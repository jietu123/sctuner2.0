from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import yaml
try:
    import psutil  # type: ignore
except Exception:
    psutil = None

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.sample_paths import resolve_sample_dir, sample_dir_candidates


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="SVTuner mainline runner (Stage1 -> Stage3 -> Stage4 -> Stage5, no Stage7)"
    )
    p.add_argument("--sample", default="real_brca", help="dataset sample id")
    p.add_argument("--project_root", default=".", help="project root")
    p.add_argument("--rscript", default=None, help="Rscript executable path (optional)")
    p.add_argument(
        "--stage1_fallback_python",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="when Stage1 Rscript crashes on Windows, fallback to Python stage1 builder",
    )
    p.add_argument(
        "--stage1_fallback_source_sample",
        default=None,
        help="optional source sample override for Stage1 Python fallback",
    )

    p.add_argument(
        "--missing_type",
        default="__NO_MISSING__",
        help="missing type name(s), comma-separated; use __NO_MISSING__ to disable CLI injection",
    )
    p.add_argument("--n_processors", type=int, default=1)
    p.add_argument("--n_subspots", type=int, default=800)
    p.add_argument("--mapping_cells_per_spot", type=int, default=None)
    p.add_argument(
        "--route2_filter_scope",
        choices=["unsupported_all", "missing_only", "missing_detected_only"],
        default="missing_only",
    )

    p.add_argument("--stage3_sc_expr_source", choices=["normalized", "data", "counts", "auto"], default="normalized")
    p.add_argument("--stage4_sc_expr_source", choices=["normalized", "data", "counts", "auto"], default="normalized")
    p.add_argument("--stage5_sim_dir", default=None, help="optional override of Stage5 truth/sim directory")
    p.add_argument(
        "--strict_stage5_truth",
        action="store_true",
        help="fail if Stage5 truth files are missing (default: skip Stage5 gracefully)",
    )

    p.add_argument("--skip_stage1", action="store_true")
    p.add_argument("--skip_stage3", action="store_true")
    p.add_argument("--skip_stage4", action="store_true")
    p.add_argument("--skip_stage5", action="store_true")
    p.add_argument("--progress_width", type=int, default=36, help="overall progress bar width")
    p.add_argument("--heartbeat_sec", type=int, default=30, help="heartbeat interval (seconds) for long-running commands")
    p.add_argument("--monitor_sec", type=float, default=1.0, help="resource polling interval in seconds")
    p.add_argument(
        "--quiet_subprocess",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="hide subprocess stdout/stderr in terminal and write logs to files",
    )
    p.add_argument("--log_dir", default=None, help="log directory (default: result/<sample>/mainline_logs)")
    return p.parse_args()


def _print_progress(done: int, total: int, stage: str, status: str, width: int = 36) -> None:
    total = max(total, 1)
    done = min(max(done, 0), total)
    width = max(width, 10)
    filled = int(round(width * done / total))
    bar = "#" * filled + "." * (width - filled)
    print(f"[PROGRESS] [{bar}] {done}/{total} | {stage} | {status}")


def _bytes_to_mb(x: int | float | None) -> float | None:
    if x is None:
        return None
    return round(float(x) / (1024.0 * 1024.0), 3)


def _sample_proc_rss_cpu(ps_proc) -> tuple[int | None, float | None]:
    """Return (rss_bytes, cpu_percent) for process + children."""
    try:
        procs = [ps_proc] + ps_proc.children(recursive=True)
    except Exception:
        procs = [ps_proc]
    rss_total = 0
    rss_seen = False
    cpu_total = 0.0
    cpu_seen = False
    for p in procs:
        try:
            mi = p.memory_info()
            rss_total += int(mi.rss)
            rss_seen = True
        except Exception:
            pass
        try:
            cpu_total += float(p.cpu_percent(interval=None))
            cpu_seen = True
        except Exception:
            pass
    return (rss_total if rss_seen else None, cpu_total if cpu_seen else None)


def _aggregate_stage_perf(cmd_metrics: list[dict[str, Any]]) -> dict[str, Any]:
    if not cmd_metrics:
        return {}
    elapsed = sum(float(m.get("elapsed_seconds", 0.0)) for m in cmd_metrics)
    peak_candidates = [m.get("peak_rss_mb") for m in cmd_metrics if m.get("peak_rss_mb") is not None]
    peak_rss_mb = max(peak_candidates) if peak_candidates else None

    cpu_weighted_sum = 0.0
    cpu_weight = 0.0
    cpu_max_candidates: list[float] = []
    for m in cmd_metrics:
        cm = m.get("cpu_percent_mean")
        cx = m.get("cpu_percent_max")
        w = float(m.get("elapsed_seconds", 0.0))
        if cm is not None and w > 0:
            cpu_weighted_sum += float(cm) * w
            cpu_weight += w
        if cx is not None:
            cpu_max_candidates.append(float(cx))
    cpu_percent_mean = (cpu_weighted_sum / cpu_weight) if cpu_weight > 0 else None
    cpu_percent_max = max(cpu_max_candidates) if cpu_max_candidates else None
    return {
        "elapsed_seconds": round(elapsed, 3),
        "peak_rss_mb": round(peak_rss_mb, 3) if peak_rss_mb is not None else None,
        "cpu_percent_mean": round(cpu_percent_mean, 3) if cpu_percent_mean is not None else None,
        "cpu_percent_max": round(cpu_percent_max, 3) if cpu_percent_max is not None else None,
        "commands": [m.get("label") for m in cmd_metrics],
    }


def _build_stage4_resource_comparison(cmd_metrics: list[dict[str, Any]]) -> dict[str, Any]:
    baseline = next((m for m in cmd_metrics if m.get("label") == "stage4_baseline"), None)
    route2 = next((m for m in cmd_metrics if m.get("label") == "stage4_route2"), None)
    if baseline is None or route2 is None:
        return {
            "available": False,
            "reason": "missing_stage4_commands",
        }

    def _metric_row(key: str) -> dict[str, Any]:
        b = baseline.get(key)
        r = route2.get(key)
        if b is None or r is None:
            return {
                "baseline": b,
                "route2": r,
                "delta_route2_minus_baseline": None,
                "route2_over_baseline": None,
            }
        b = float(b)
        r = float(r)
        return {
            "baseline": round(b, 3),
            "route2": round(r, 3),
            "delta_route2_minus_baseline": round(r - b, 3),
            "route2_over_baseline": round((r / b), 4) if b != 0 else None,
        }

    return {
        "available": True,
        "metrics": {
            "elapsed_seconds": _metric_row("elapsed_seconds"),
            "peak_rss_mb": _metric_row("peak_rss_mb"),
            "cpu_percent_mean": _metric_row("cpu_percent_mean"),
            "cpu_percent_max": _metric_row("cpu_percent_max"),
        },
    }


def _write_stage4_resource_comparison_csv(path: Path, cmp_dict: dict[str, Any]) -> None:
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "metric",
                "baseline",
                "route2",
                "delta_route2_minus_baseline",
                "route2_over_baseline",
            ],
        )
        w.writeheader()
        metrics = (cmp_dict.get("metrics") or {}) if cmp_dict.get("available") else {}
        for metric, row in metrics.items():
            w.writerow(
                {
                    "metric": metric,
                    "baseline": row.get("baseline"),
                    "route2": row.get("route2"),
                    "delta_route2_minus_baseline": row.get("delta_route2_minus_baseline"),
                    "route2_over_baseline": row.get("route2_over_baseline"),
                }
            )


def run_cmd(
    cmd: list[str],
    cwd: Path,
    *,
    label: str = "command",
    heartbeat_sec: int = 30,
    monitor_sec: float = 1.0,
    quiet_subprocess: bool = True,
    log_file: Path | None = None,
) -> dict[str, Any]:
    print(f"[RUN] {' '.join(str(x) for x in cmd)}")
    if quiet_subprocess:
        if log_file is None:
            raise ValueError("log_file is required when quiet_subprocess=True")
        log_file.parent.mkdir(parents=True, exist_ok=True)
        log_fh = log_file.open("w", encoding="utf-8", errors="ignore")
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            text=True,
        )
    else:
        log_fh = None
        proc = subprocess.Popen(cmd, cwd=str(cwd))
    start = time.time()
    heartbeat_sec = max(int(heartbeat_sec), 5)
    monitor_sec = max(float(monitor_sec), 0.5)
    last_heartbeat = start

    peak_rss_bytes: int | None = None
    cpu_sum = 0.0
    cpu_n = 0
    cpu_max: float | None = None
    monitor_enabled = False
    ps_proc = None
    if psutil is not None:
        try:
            ps_proc = psutil.Process(proc.pid)
            # prime cpu counters
            _sample_proc_rss_cpu(ps_proc)
            monitor_enabled = True
        except Exception:
            monitor_enabled = False

    while True:
        ret = proc.poll()
        if ret is not None:
            break
        now = time.time()
        if monitor_enabled and ps_proc is not None:
            rss, cpu = _sample_proc_rss_cpu(ps_proc)
            if rss is not None:
                peak_rss_bytes = rss if peak_rss_bytes is None else max(peak_rss_bytes, rss)
            if cpu is not None:
                cpu_sum += float(cpu)
                cpu_n += 1
                cpu_max = float(cpu) if cpu_max is None else max(cpu_max, float(cpu))
        if now - last_heartbeat >= heartbeat_sec:
            elapsed = int(now - start)
            print(f"[HEARTBEAT] {label} is still running... elapsed={elapsed}s")
            last_heartbeat = now
        time.sleep(monitor_sec)
    if log_fh is not None:
        log_fh.close()
    if int(ret) != 0:
        if quiet_subprocess and log_file is not None:
            print(f"[ERROR] {label} failed. log -> {log_file}")
        raise RuntimeError(f"command failed (exit={ret}): {' '.join(str(x) for x in cmd)}")
    end = time.time()
    elapsed = end - start
    # one final sample close to process end
    if monitor_enabled and ps_proc is not None:
        rss, cpu = _sample_proc_rss_cpu(ps_proc)
        if rss is not None:
            peak_rss_bytes = rss if peak_rss_bytes is None else max(peak_rss_bytes, rss)
        if cpu is not None:
            cpu_sum += float(cpu)
            cpu_n += 1
            cpu_max = float(cpu) if cpu_max is None else max(cpu_max, float(cpu))

    cpu_mean = (cpu_sum / cpu_n) if cpu_n > 0 else None
    cmd_perf = {
        "label": label,
        "elapsed_seconds": round(elapsed, 3),
        "peak_rss_mb": _bytes_to_mb(peak_rss_bytes),
        "cpu_percent_mean": round(cpu_mean, 3) if cpu_mean is not None else None,
        "cpu_percent_max": round(cpu_max, 3) if cpu_max is not None else None,
        "monitor_enabled": bool(monitor_enabled),
        "samples": int(cpu_n),
        "log_file": str(log_file) if log_file is not None else None,
    }

    print(
        f"[OK] {label} finished in {int(round(elapsed))}s"
        + (
            f" | peak_mem={cmd_perf['peak_rss_mb']}MB cpu_avg={cmd_perf['cpu_percent_mean']}% cpu_max={cmd_perf['cpu_percent_max']}%"
            if cmd_perf["monitor_enabled"]
            else ""
        )
    )
    if quiet_subprocess and log_file is not None:
        print(f"[LOG] {label} -> {log_file}")
    return cmd_perf


def load_project_cfg(project_root: Path) -> dict[str, Any]:
    cfg = project_root / "configs" / "project_config.yaml"
    if not cfg.exists():
        return {}
    return yaml.safe_load(cfg.read_text(encoding="utf-8")) or {}


def resolve_rscript_path(rscript_raw: str) -> str:
    """
    On Windows, prefer <env>\\Lib\\R\\bin\\x64\\Rscript.exe over <env>\\Scripts\\Rscript.exe
    to avoid mingw runtime pseudo-relocation crashes.
    """
    rs = Path(str(rscript_raw))
    try:
        parts_lower = [p.lower() for p in rs.parts]
    except Exception:
        return str(rscript_raw)

    if rs.name.lower() != "rscript.exe":
        return str(rscript_raw)
    if "scripts" not in parts_lower:
        return str(rscript_raw)

    try:
        idx = parts_lower.index("scripts")
        env_root = Path(*rs.parts[:idx])
        cand = env_root / "Lib" / "R" / "bin" / "x64" / "Rscript.exe"
        if cand.exists():
            print(f"[INFO] rscript redirected: {rs} -> {cand}")
            return str(cand)
    except Exception:
        pass
    return str(rscript_raw)


def infer_sim_dir_from_dataset(project_root: Path, sample: str) -> Path:
    ds_cfg = project_root / "configs" / "datasets" / f"{sample}.yaml"
    if ds_cfg.exists():
        data = yaml.safe_load(ds_cfg.read_text(encoding="utf-8")) or {}
        paths = data.get("paths") or {}
        st_expr = paths.get("st_expr")
        if st_expr:
            p = Path(st_expr)
            if p.is_absolute() and p.exists():
                return p.parent
            roots = sample_dir_candidates(project_root, sample, sim_group="real_brca")
            for r in roots:
                cand = (r / p).resolve()
                if cand.exists():
                    return cand.parent
            cand = (project_root / str(st_expr)).resolve()
            if cand.exists():
                return cand.parent

    # Fallback: infer directly by sample directory resolution.
    return resolve_sample_dir(project_root, sample, sim_group="real_brca", must_exist=True)


def resolve_stage5_truth_dir(project_root: Path, sample: str, preferred: Path | None) -> Path | None:
    """Pick a directory that contains Stage5 truth files."""
    stage1_export = project_root / "data" / "processed" / sample / "stage1_preprocess" / "exported"
    candidates = []
    if preferred is not None:
        candidates.append(preferred)
    candidates.append(stage1_export)
    for c in candidates:
        if c is None:
            continue
        if not c.exists():
            continue
        if (c / "sim_info.json").exists() and (c / "sim_truth_spot_type_fraction.csv").exists():
            return c
    return None


def _parse_missing_type_arg(raw: str | None) -> list[str]:
    if raw is None:
        return []
    txt = str(raw).strip()
    if not txt or txt == "__NO_MISSING__":
        return []
    txt = txt.replace(";", ",")
    out: list[str] = []
    for part in txt.split(","):
        t = str(part).strip()
        if not t:
            continue
        if t not in out:
            out.append(t)
    return out


def _read_sim_info(sim_dir: Path) -> dict[str, Any]:
    sim_info_path = sim_dir / "sim_info.json"
    if not sim_info_path.exists():
        return {}
    try:
        data = json.loads(sim_info_path.read_text(encoding="utf-8"))
        if isinstance(data, dict):
            return data
    except Exception:
        pass
    return {}


def _missing_types_from_sim_info(sim_info: dict[str, Any]) -> list[str]:
    out: list[str] = []
    mt_list = sim_info.get("missing_types")
    if isinstance(mt_list, list):
        for x in mt_list:
            t = str(x).strip()
            if t and t not in out:
                out.append(t)
    mt_single = sim_info.get("missing_type")
    if isinstance(mt_single, str):
        t = mt_single.strip()
        if t and t not in out:
            out.append(t)
    return out


def _infer_missing_types_from_source_chain(project_root: Path, sample: str, sim_dir: Path) -> list[str]:
    out: list[str] = []
    visited: set[str] = set()

    current_sample: str | None = str(sample).strip() or None
    first = True
    for _ in range(20):
        if not current_sample:
            break
        if current_sample in visited:
            break
        visited.add(current_sample)

        candidate_dirs: list[Path] = []
        if first:
            candidate_dirs.append(sim_dir)
        candidate_dirs.extend(sample_dir_candidates(project_root, current_sample, sim_group="real_brca"))

        sim_info: dict[str, Any] = {}
        for c in candidate_dirs:
            si = _read_sim_info(c)
            if si:
                sim_info = si
                break
        if not sim_info:
            break

        mts = _missing_types_from_sim_info(sim_info)
        for t in mts:
            if t not in out:
                out.append(t)

        src = sim_info.get("source_sample")
        if not isinstance(src, str) or not src.strip():
            break
        current_sample = src.strip()
        first = False
    return out


def _read_text_safe(path: Path, max_chars: int = 20000) -> str:
    try:
        if not path.exists():
            return ""
        txt = path.read_text(encoding="utf-8", errors="ignore")
        if len(txt) <= max_chars:
            return txt
        return txt[-max_chars:]
    except Exception:
        return ""


def _looks_like_windows_r_mingw_crash(err_msg: str, log_path: Path | None) -> bool:
    s = (err_msg or "").lower()
    if (
        "3221226505" in s
        or "3221225781" in s
        or "-1073740791" in s
        or "-1073741515" in s
    ):
        return True
    if "mingw-w64 runtime failure" in s or "pseudo relocation" in s:
        return True
    if log_path is not None:
        t = _read_text_safe(log_path).lower()
        if "mingw-w64 runtime failure" in t or "pseudo relocation" in t:
            return True
    return False


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    py = Path(sys.executable).resolve()
    monitoring_reason = "enabled"
    if psutil is None:
        monitoring_reason = "psutil_not_installed"
        print("[WARN] psutil is not installed; memory/cpu metrics are disabled (time metrics still available).")

    result_root = project_root / "result" / args.sample
    result_root.mkdir(parents=True, exist_ok=True)
    log_dir = Path(args.log_dir).resolve() if args.log_dir else (result_root / "mainline_logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    run_started = time.time()

    project_cfg = load_project_cfg(project_root)
    rscript = args.rscript or project_cfg.get("rscript_path") or "Rscript"
    if str(rscript).lower().endswith("rscript.exe"):
        rscript = resolve_rscript_path(str(rscript))
    sim_dir = Path(args.stage5_sim_dir).resolve() if args.stage5_sim_dir else infer_sim_dir_from_dataset(project_root, args.sample)
    cli_missing_types = _parse_missing_type_arg(args.missing_type)
    inferred_missing_types = _infer_missing_types_from_source_chain(project_root, args.sample, sim_dir)
    effective_missing_types: list[str] = []
    for t in cli_missing_types + inferred_missing_types:
        if t not in effective_missing_types:
            effective_missing_types.append(t)
    stage4_missing_arg = ",".join(effective_missing_types) if effective_missing_types else "__NO_MISSING__"
    effective_route2_filter_scope = args.route2_filter_scope
    if (
        effective_route2_filter_scope in {"missing_only", "missing_detected_only"}
        and not effective_missing_types
    ):
        print(
            "[WARN] no missing type found from CLI or sim_info chain; "
            f"route2_filter_scope auto-fallback: {effective_route2_filter_scope} -> unsupported_all"
        )
        effective_route2_filter_scope = "unsupported_all"

    stage_order = ["stage1", "stage3", "stage4", "stage5"]
    stage_status: dict[str, str] = {}
    cmd_perf_all: list[dict[str, Any]] = []
    stage_cmd_perf: dict[str, list[dict[str, Any]]] = {k: [] for k in stage_order}
    stage5_mode = "skipped"
    stage5_outputs: dict[str, str] = {}
    total_stages = len(stage_order)
    done_stages = 0

    try:
        _print_progress(done_stages, total_stages, "mainline", "started", args.progress_width)
        if not args.skip_stage1:
            _print_progress(done_stages, total_stages, "stage1", "running", args.progress_width)
            stage1_log = log_dir / "stage1_preprocess.log"
            try:
                m = run_cmd(
                    [
                        str(rscript),
                        str(project_root / "r_scripts" / "stage1_preprocess.R"),
                        "--sample",
                        args.sample,
                        "--project_root",
                        str(project_root),
                        "--export_csv",
                    ],
                    cwd=project_root,
                    label="stage1_preprocess",
                    heartbeat_sec=args.heartbeat_sec,
                    monitor_sec=args.monitor_sec,
                    quiet_subprocess=args.quiet_subprocess,
                    log_file=stage1_log,
                )
                cmd_perf_all.append(m)
                stage_cmd_perf["stage1"].append(m)
                stage_status["stage1"] = "ok"
            except Exception as stage1_err:
                can_fallback = bool(args.stage1_fallback_python) and _looks_like_windows_r_mingw_crash(
                    str(stage1_err),
                    stage1_log,
                )
                if not can_fallback:
                    raise
                print("[WARN] Stage1 Rscript crashed with Windows mingw runtime issue; trying Python fallback.")
                cmd_fb = [
                    str(py),
                    str(project_root / "scripts" / "prepare_stage1_from_sim_source.py"),
                    "--sample",
                    args.sample,
                    "--project_root",
                    str(project_root),
                ]
                if args.stage1_fallback_source_sample:
                    cmd_fb.extend(["--source_sample", str(args.stage1_fallback_source_sample)])
                m_fb = run_cmd(
                    cmd_fb,
                    cwd=project_root,
                    label="stage1_preprocess_fallback",
                    heartbeat_sec=args.heartbeat_sec,
                    monitor_sec=args.monitor_sec,
                    quiet_subprocess=args.quiet_subprocess,
                    log_file=log_dir / "stage1_preprocess_fallback.log",
                )
                cmd_perf_all.append(m_fb)
                stage_cmd_perf["stage1"].append(m_fb)
                stage_status["stage1"] = "ok_fallback"
        else:
            stage_status["stage1"] = "skipped"
        done_stages += 1
        _print_progress(done_stages, total_stages, "stage1", stage_status["stage1"], args.progress_width)

        if not args.skip_stage3:
            _print_progress(done_stages, total_stages, "stage3", "running", args.progress_width)
            m = run_cmd(
                [
                    str(py),
                    "-m",
                    "src.stages.stage3_type_plugin",
                    "--sample",
                    args.sample,
                    "--sc_expr_source",
                    args.stage3_sc_expr_source,
                ],
                cwd=project_root,
                label="stage3_type_plugin",
                heartbeat_sec=args.heartbeat_sec,
                monitor_sec=args.monitor_sec,
                quiet_subprocess=args.quiet_subprocess,
                log_file=log_dir / "stage3_type_plugin.log",
            )
            cmd_perf_all.append(m)
            stage_cmd_perf["stage3"].append(m)
            stage_status["stage3"] = "ok"
        else:
            stage_status["stage3"] = "skipped"
        done_stages += 1
        _print_progress(done_stages, total_stages, "stage3", stage_status["stage3"], args.progress_width)

        if not args.skip_stage4:
            _print_progress(done_stages, total_stages, "stage4", "running", args.progress_width)
            baseline_stage4_root = result_root / "stage4_cytospace_baseline"
            route2_stage4_root = result_root / "stage4_cytospace_route2"
            for d in [baseline_stage4_root, route2_stage4_root]:
                if d.exists():
                    shutil.rmtree(d)
                    print(f"[CLEAN] removed old stage4 dir: {d}")
            cmd_baseline = [
                str(py),
                "-m",
                "src.stages.stage4_cytospace",
                "--sample",
                args.sample,
                "--filter_mode",
                "none",
                "--cell_type_column",
                "sc_meta",
                "--stage4_suffix",
                "_baseline",
                "--missing_type",
                stage4_missing_arg,
                "--n_processors",
                str(args.n_processors),
                "--n_subspots",
                str(args.n_subspots),
                "--sc_expr_source",
                args.stage4_sc_expr_source,
            ]
            if args.mapping_cells_per_spot is not None:
                cmd_baseline.extend(["--mapping_cells_per_spot", str(args.mapping_cells_per_spot)])
            m = run_cmd(
                cmd_baseline,
                cwd=project_root,
                label="stage4_baseline",
                heartbeat_sec=args.heartbeat_sec,
                monitor_sec=args.monitor_sec,
                quiet_subprocess=args.quiet_subprocess,
                log_file=log_dir / "stage4_baseline.log",
            )
            cmd_perf_all.append(m)
            stage_cmd_perf["stage4"].append(m)

            cmd_route2 = [
                str(py),
                "-m",
                "src.stages.stage4_cytospace",
                "--sample",
                args.sample,
                "--filter_mode",
                "plugin_unknown",
                "--filter_scope",
                effective_route2_filter_scope,
                "--cell_type_column",
                "plugin_type",
                "--stage4_suffix",
                "_route2",
                "--missing_type",
                stage4_missing_arg,
                "--n_processors",
                str(args.n_processors),
                "--n_subspots",
                str(args.n_subspots),
                "--sc_expr_source",
                args.stage4_sc_expr_source,
            ]
            if args.mapping_cells_per_spot is not None:
                cmd_route2.extend(["--mapping_cells_per_spot", str(args.mapping_cells_per_spot)])
            m = run_cmd(
                cmd_route2,
                cwd=project_root,
                label="stage4_route2",
                heartbeat_sec=args.heartbeat_sec,
                monitor_sec=args.monitor_sec,
                quiet_subprocess=args.quiet_subprocess,
                log_file=log_dir / "stage4_route2.log",
            )
            cmd_perf_all.append(m)
            stage_cmd_perf["stage4"].append(m)
            stage_status["stage4"] = "ok"
        else:
            stage_status["stage4"] = "skipped"
        done_stages += 1
        _print_progress(done_stages, total_stages, "stage4", stage_status["stage4"], args.progress_width)

        if not args.skip_stage5:
            _print_progress(done_stages, total_stages, "stage5", "running", args.progress_width)
            stage5_truth_dir = resolve_stage5_truth_dir(project_root, args.sample, sim_dir)
            if stage5_truth_dir is None:
                if args.strict_stage5_truth:
                    msg = (
                        "[ERROR] Stage5 truth files not found (need sim_info.json + sim_truth_spot_type_fraction.csv). "
                        f"checked around: {sim_dir}"
                    )
                    raise FileNotFoundError(msg)
                print(
                    "[INFO] sim truth not found; switching to real-data Stage5 summary "
                    "(src.stages.stage5_real_summary)"
                )
                stage4_baseline = result_root / "stage4_cytospace_baseline" / "cytospace_output"
                stage4_route2 = result_root / "stage4_cytospace_route2" / "cytospace_output"
                m = run_cmd(
                    [
                        str(py),
                        "-m",
                        "src.stages.stage5_real_summary",
                        "--sample",
                        args.sample,
                        "--project_root",
                        str(project_root),
                        "--baseline_dir",
                        str(stage4_baseline),
                        "--route2_dir",
                        str(stage4_route2),
                        "--out_dir",
                        str(result_root / "stage5_real_summary"),
                    ],
                    cwd=project_root,
                    label="stage5_real_summary",
                    heartbeat_sec=args.heartbeat_sec,
                    monitor_sec=args.monitor_sec,
                    quiet_subprocess=args.quiet_subprocess,
                    log_file=log_dir / "stage5_real_summary.log",
                )
                cmd_perf_all.append(m)
                stage_cmd_perf["stage5"].append(m)
                stage_status["stage5"] = "ok_real_summary"
                stage5_mode = "real_summary"
                stage5_outputs = {
                    "stage5_real_summary_json": str(result_root / "stage5_real_summary" / "stage5_real_summary.json"),
                    "stage5_real_type_mass_csv": str(result_root / "stage5_real_summary" / "stage5_real_type_mass.csv"),
                    "stage5_real_spot_metrics_csv": str(result_root / "stage5_real_summary" / "stage5_real_spot_metrics.csv"),
                }
            else:
                stage4_baseline = result_root / "stage4_cytospace_baseline" / "cytospace_output"
                stage4_route2 = result_root / "stage4_cytospace_route2" / "cytospace_output"

                stage5_baseline = result_root / "stage5_eval_baseline"
                stage5_route2 = result_root / "stage5_eval_route2"
                stage5_compare = result_root / "stage5_compare_baseline_vs_route2.json"

                m = run_cmd(
                    [
                        str(py),
                        "-m",
                        "src.stages.stage5_route2_s0",
                        "--sample",
                        args.sample,
                        "--run_tag",
                        "baseline",
                        "--stage4_dir",
                        str(stage4_baseline),
                        "--sim_dir",
                        str(stage5_truth_dir),
                        "--out_dir",
                        str(stage5_baseline),
                    ],
                    cwd=project_root,
                    label="stage5_sim_baseline",
                    heartbeat_sec=args.heartbeat_sec,
                    monitor_sec=args.monitor_sec,
                    quiet_subprocess=args.quiet_subprocess,
                    log_file=log_dir / "stage5_sim_baseline.log",
                )
                cmd_perf_all.append(m)
                stage_cmd_perf["stage5"].append(m)

                m = run_cmd(
                    [
                        str(py),
                        "-m",
                        "src.stages.stage5_route2_s0",
                        "--sample",
                        args.sample,
                        "--run_tag",
                        "route2",
                        "--stage4_dir",
                        str(stage4_route2),
                        "--sim_dir",
                        str(stage5_truth_dir),
                        "--out_dir",
                        str(stage5_route2),
                    ],
                    cwd=project_root,
                    label="stage5_sim_route2",
                    heartbeat_sec=args.heartbeat_sec,
                    monitor_sec=args.monitor_sec,
                    quiet_subprocess=args.quiet_subprocess,
                    log_file=log_dir / "stage5_sim_route2.log",
                )
                cmd_perf_all.append(m)
                stage_cmd_perf["stage5"].append(m)

                baseline_json = stage5_baseline / "stage5_route2_s0__baseline.json"
                route2_json = stage5_route2 / "stage5_route2_s0__route2.json"
                if baseline_json.exists() and route2_json.exists():
                    m = run_cmd(
                        [
                            str(py),
                            "-m",
                            "src.stages.stage5_route2_s0",
                            "--compare_baseline",
                            str(baseline_json),
                            "--compare_route2",
                            str(route2_json),
                            "--compare_out",
                            str(stage5_compare),
                        ],
                        cwd=project_root,
                        label="stage5_sim_compare",
                        heartbeat_sec=args.heartbeat_sec,
                        monitor_sec=args.monitor_sec,
                        quiet_subprocess=args.quiet_subprocess,
                        log_file=log_dir / "stage5_sim_compare.log",
                    )
                    cmd_perf_all.append(m)
                    stage_cmd_perf["stage5"].append(m)
                stage_status["stage5"] = "ok"
                stage5_mode = "sim_truth_eval"
                stage5_outputs = {
                    "stage5_baseline_json": str(result_root / "stage5_eval_baseline" / "stage5_route2_s0__baseline.json"),
                    "stage5_route2_json": str(result_root / "stage5_eval_route2" / "stage5_route2_s0__route2.json"),
                    "stage5_compare_json": str(result_root / "stage5_compare_baseline_vs_route2.json"),
                }
        else:
            stage_status["stage5"] = "skipped"
            stage5_mode = "skipped"
        done_stages += 1
        _print_progress(done_stages, total_stages, "stage5", stage_status["stage5"], args.progress_width)

        stage_perf = {k: _aggregate_stage_perf(v) for k, v in stage_cmd_perf.items()}
        overall_perf = _aggregate_stage_perf(cmd_perf_all)
        stage4_resource_cmp = _build_stage4_resource_comparison(cmd_perf_all)
        stage4_resource_cmp_json = result_root / "stage4_resource_comparison.json"
        stage4_resource_cmp_csv = result_root / "stage4_resource_comparison.csv"
        stage4_resource_cmp_json.write_text(
            json.dumps(stage4_resource_cmp, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        _write_stage4_resource_comparison_csv(stage4_resource_cmp_csv, stage4_resource_cmp)
        summary = {
            "sample": args.sample,
            "status": "ok",
            "elapsed_seconds": round(time.time() - run_started, 2),
            "stages": stage_status,
            "result_root": str(result_root),
            "log_dir": str(log_dir),
            "stage4_baseline": str(result_root / "stage4_cytospace_baseline" / "cytospace_output"),
            "stage4_route2": str(result_root / "stage4_cytospace_route2" / "cytospace_output"),
            "stage4_missing_types_cli": cli_missing_types,
            "stage4_missing_types_inferred": inferred_missing_types,
            "stage4_missing_types_effective": effective_missing_types,
            "stage4_filter_scope_requested": args.route2_filter_scope,
            "stage4_filter_scope_effective": effective_route2_filter_scope,
            "stage4_resource_comparison_json": str(stage4_resource_cmp_json),
            "stage4_resource_comparison_csv": str(stage4_resource_cmp_csv),
            "stage5_mode": stage5_mode,
            "stage5_outputs": stage5_outputs,
            "performance": {
                "monitoring_enabled": bool(psutil is not None),
                "monitoring_reason": monitoring_reason,
                "monitor_sec": float(args.monitor_sec),
                "heartbeat_sec": int(args.heartbeat_sec),
                "commands": cmd_perf_all,
                "stages": stage_perf,
                "overall": overall_perf,
            },
            "stage4_resource_comparison": stage4_resource_cmp,
        }
        summary_path = result_root / "mainline_run_summary.json"
        summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
        if stage4_resource_cmp.get("available"):
            m = stage4_resource_cmp.get("metrics", {})
            print(
                "[COMPARE] Stage4 Baseline vs Route2 | "
                f"time(s): {m.get('elapsed_seconds', {}).get('baseline')} -> {m.get('elapsed_seconds', {}).get('route2')} | "
                f"mem(MB): {m.get('peak_rss_mb', {}).get('baseline')} -> {m.get('peak_rss_mb', {}).get('route2')} | "
                f"cpu_avg(%): {m.get('cpu_percent_mean', {}).get('baseline')} -> {m.get('cpu_percent_mean', {}).get('route2')}"
            )
            print(f"[COMPARE] details -> {stage4_resource_cmp_csv}")
        _print_progress(total_stages, total_stages, "mainline", "completed", args.progress_width)
        print(f"[DONE] mainline completed. summary -> {summary_path}")
        return 0
    except Exception as e:
        fail = {
            "sample": args.sample,
            "status": "failed",
            "error": str(e),
            "elapsed_seconds": round(time.time() - run_started, 2),
            "stages": stage_status,
        }
        fail_path = result_root / "mainline_run_summary.failed.json"
        fail_path.write_text(json.dumps(fail, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"[ERROR] {e}")
        print(f"[ERROR] failure summary -> {fail_path}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

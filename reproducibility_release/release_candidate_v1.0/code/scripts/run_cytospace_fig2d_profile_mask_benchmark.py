from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


PAIRS = {
    "melanoma_slide1": {
        "source_sample": "cytospace_fig2c_melanoma_mel1_rep2",
        "target_type": "Endothelial cells",
        "tumor_labels": "Melanoma|Melanoma cells",
    },
    "melanoma_slide2": {
        "source_sample": "cytospace_fig2c_melanoma_mel2_rep1",
        "target_type": "NK cells",
        "tumor_labels": "Melanoma|Melanoma cells",
    },
    "brca_er_her2_fresh_frozen": {
        "source_sample": "cytospace_fig2d_tme_brca_er_her2_fresh_frozen",
        "target_type": "T-cells",
        "tumor_labels": "Epithelial cells",
    },
    "brca_her2_ffpe": {
        "source_sample": "cytospace_fig2d_tme_brca_her2_ffpe",
        "target_type": "T-cells",
        "tumor_labels": "Epithelial cells",
    },
    "brca_tnbc_fresh_frozen": {
        "source_sample": "cytospace_fig2d_tme_brca_tnbc_fresh_frozen",
        "target_type": "Fibroblasts",
        "tumor_labels": "Epithelial cells",
    },
    "crc_fresh_frozen": {
        "source_sample": "cytospace_fig2d_tme_crc_fresh_frozen",
        "target_type": "Monocytes and Macrophages",
        "tumor_labels": "Epithelial cells",
    },
}


def _safe_name(label: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in str(label)).strip("_")


def _run(cmd: list[str], log_path: Path, env: dict[str, str] | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env)
    if proc.returncode != 0:
        tail = "\n".join(log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-80:])
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}\n--- log tail ---\n{tail}")


def _ensure_flat_processed_link(root: Path, sample: str, group: str = "cytospace_fig2d_profile_mask") -> None:
    target = root / "data" / "processed" / group / sample
    link = root / "data" / "processed" / sample
    if not target.exists():
        raise FileNotFoundError(target)
    if link.exists():
        if link.resolve() == target.resolve():
            return
        # Existing profile-mask links may be Windows junctions whose textual
        # metadata is locale-dependent. If a path already exists, leave it
        # untouched; this runner only needs Stage4's legacy flat lookup to find
        # the Stage3 files.
        return
    link.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["cmd", "/c", "mklink", "/J", str(link), str(target)], check=True, stdout=subprocess.DEVNULL)


def _read_stage3_detection(root: Path, sample: str, target_type: str) -> dict:
    candidates = [
        root / "result" / "cytospace_fig2d_profile_mask" / sample / "stage3_typematch" / "stage3_summary.json",
        root / "result" / sample / "stage3_typematch" / "stage3_summary.json",
    ]
    path = next((p for p in candidates if p.exists()), candidates[0])
    if not path.exists():
        return {
            "stage3_summary_path": str(path),
            "stage3_missing_types": [],
            "target_detected_by_stage3": False,
        }
    summary = json.loads(path.read_text(encoding="utf-8-sig"))
    action = summary.get("action_overview", {}) or {}
    missing = [str(x) for x in action.get("missing_types", []) or []]
    return {
        "stage3_summary_path": str(path),
        "stage3_missing_types": missing,
        "target_detected_by_stage3": str(target_type) in set(missing),
        "stage3_dropped_cells": int((action.get("by_cell_count", {}) or {}).get("Dropped", 0)),
    }


def _mirror_stage3_result_for_stage4(root: Path, sample: str, group: str = "cytospace_fig2d_profile_mask") -> None:
    src = root / "result" / group / sample / "stage3_typematch"
    dst = root / "result" / sample / "stage3_typematch"
    if not src.exists():
        return
    dst.mkdir(parents=True, exist_ok=True)
    for name in ["stage3_summary.json"]:
        p = src / name
        if p.exists():
            (dst / name).write_bytes(p.read_bytes())


def _stage4(
    root: Path,
    py: str,
    sample: str,
    suffix: str,
    mode: str,
    log_path: Path,
    cells_per_spot: int,
    n_subspots: int,
) -> Path:
    if mode == "baseline":
        missing_arg = "__NO_MISSING__"
        filter_mode = "none"
        cell_col = "sc_meta"
        scope = "unsupported_all"
    elif mode == "route2":
        missing_arg = "__AUTO__"
        filter_mode = "plugin_unknown"
        cell_col = "plugin_type"
        scope = "missing_detected_only"
    else:
        raise ValueError(mode)
    env = os.environ.copy()
    env["CYTOSPACE_SKIP_ASSIGNED_EXPRESSION"] = "1"
    cmd = [
        py,
        "-m",
        "src.stages.stage4_cytospace",
        "--sample",
        sample,
        "--project_root",
        str(root),
        "--missing_type",
        missing_arg,
        "--n_processors",
        "1",
        "--n_subspots",
        str(n_subspots),
        "--mapping_cells_per_spot",
        str(cells_per_spot),
        "--sc_expr_source",
        "normalized",
        "--filter_mode",
        filter_mode,
        "--cell_type_column",
        cell_col,
        "--filter_scope",
        scope,
        "--stage4_suffix",
        suffix,
    ]
    _run(cmd, log_path, env=env)
    return root / "result" / sample / f"stage4_cytospace{suffix}" / "cytospace_output"


def _build_assigned_expression(root: Path, py: str, sample: str, cyto_out: Path, log_path: Path) -> None:
    _run(
        [
            py,
            "scripts/build_assigned_expression_from_stage1.py",
            "--project_root",
            str(root),
            "--sample",
            sample,
            "--cytospace_output",
            str(cyto_out),
            "--sc_expr_source",
            "normalized",
        ],
        log_path,
    )


def _enrich(
    root: Path,
    py: str,
    rscript: str,
    cyto_out: Path,
    out_dir: Path,
    tumor_labels: str,
    slide_label: str,
    log_path: Path,
    backend: str,
) -> None:
    if backend == "python":
        _run(
            [
                py,
                "scripts/compute_fig2c_enrichment_python.py",
                "--project_root",
                str(root),
                "--cytospace_dir",
                str(cyto_out),
                "--out_dir",
                str(out_dir),
                "--tumor_labels",
                tumor_labels,
                "--slide_label",
                slide_label,
                "--gene_set_name",
                "Exhaustion",
                "--cell_types",
                "CD4 T cells|CD8 T cells",
                "--nperm",
                "1000",
            ],
            log_path,
        )
        return

    env = os.environ.copy()
    rscript_path = Path(rscript)
    if rscript_path.exists():
        env_prefix = rscript_path.parents[1]
        prepend = [
            str(env_prefix),
            str(env_prefix / "Library" / "mingw-w64" / "bin"),
            str(env_prefix / "Library" / "usr" / "bin"),
            str(env_prefix / "Library" / "bin"),
            str(env_prefix / "Scripts"),
            str(env_prefix / "Lib" / "R" / "bin"),
        ]
        env["PATH"] = os.pathsep.join(prepend + [env.get("PATH", "")])
    _run(
        [
            rscript,
            "scripts/compute_cytospace_fig2c_official_enrichment_seurat.R",
            str(root),
            str(cyto_out),
            str(out_dir),
            tumor_labels,
            slide_label,
            "Exhaustion",
            "CD4 T cells|CD8 T cells",
        ],
        log_path,
        env=env,
    )


def _run_other_method(
    root: Path,
    py: str,
    sample: str,
    method: str,
    mapping_dir: Path,
    log_path: Path,
) -> None:
    env = os.environ.copy()
    env.pop("PYTHONNOUSERSITE", None)
    env.setdefault("NUMBA_CACHE_DIR", str(root / ".numba_cache"))
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    Path(env["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
    if method == "tangram_marker":
        cmd = [
            py,
            "scripts/run_tangram_marker_mapping.py",
            "--project_root",
            str(root),
            "--group",
            "cytospace_fig2d_profile_mask",
            "--sample",
            sample,
            "--top_n_marker",
            "50",
            "--num_epochs",
            "200",
            "--device",
            "cpu",
            "--out_dir",
            str(mapping_dir),
        ]
    elif method == "celltrek":
        cmd = [
            py,
            "scripts/run_celltrek_mapping.py",
            "--project_root",
            str(root),
            "--group",
            "cytospace_fig2d_profile_mask",
            "--sample",
            sample,
            "--max_genes",
            "2000",
            "--n_pcs",
            "30",
            "--ntree",
            "500",
            "--out_dir",
            str(mapping_dir),
        ]
    else:
        raise ValueError(method)
    _run(cmd, log_path, env=env)


def _build_cytospace_like(root: Path, py: str, sample: str, mapping_dir: Path, out_dir: Path, log_path: Path) -> None:
    _run(
        [
            py,
            "scripts/build_cytospace_output_from_generic_mapping.py",
            "--project_root",
            str(root),
            "--sample",
            sample,
            "--mapping_dir",
            str(mapping_dir),
            "--out_dir",
            str(out_dir),
            "--sc_expr_source",
            "normalized",
        ],
        log_path,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Run non-forced profile-mask Fig.2d benchmark.")
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--pair", action="append", choices=sorted(PAIRS), default=[])
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--n_subspots", type=int, default=300)
    parser.add_argument("--cells_per_spot", type=int, default=20)
    parser.add_argument("--skip_other_methods", action="store_true")
    parser.add_argument("--enrichment_backend", choices=["python", "r"], default="python")
    parser.add_argument(
        "--target_override",
        action="append",
        default=[],
        help="Override one pair target as pair_id=Cell type for screening.",
    )
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    py = str(Path(args.python).resolve()) if Path(args.python).exists() else args.python
    result_root = root / "result" / "cytospace_fig2d_profile_mask_benchmark"
    log_root = root / "logs" / "fig2d_profile_mask_benchmark"
    pair_specs = {k: dict(v) for k, v in PAIRS.items()}
    for item in args.target_override:
        if "=" not in item:
            raise ValueError(f"--target_override must be pair_id=Cell type, got: {item}")
        pair_id, target = item.split("=", 1)
        pair_id = pair_id.strip()
        if pair_id not in pair_specs:
            raise ValueError(f"unknown pair in --target_override: {pair_id}")
        pair_specs[pair_id]["target_type"] = target.strip()
    pairs = args.pair or list(pair_specs)
    manifest_path = result_root / "profile_mask_manifest.json"
    if manifest_path.exists():
        existing = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
        manifest_by_pair = {str(row["pair_id"]): row for row in existing}
    else:
        manifest_by_pair = {}

    for pair_id in pairs:
        spec = pair_specs[pair_id]
        target_type = spec["target_type"]
        target_sample = f"{spec['source_sample']}_profile_mask_{_safe_name(target_type)}"
        pair_result = result_root / pair_id
        print(f"[STEP] {pair_id}: profile-mask target={target_type}")

        build_cmd = [
            py,
            "scripts/build_cytospace_fig2d_profile_mask_sample.py",
            "--project_root",
            str(root),
            "--source_sample",
            spec["source_sample"],
            "--target_sample",
            target_sample,
            "--target_type",
            target_type,
            "--target_profile_impact",
            "0.45",
            "--allow_unsafe_panel",
        ]
        if args.overwrite:
            build_cmd.append("--overwrite")
        _run(build_cmd, log_root / f"{pair_id}.build.log")
        _ensure_flat_processed_link(root, target_sample)

        _run(
            [
                py,
                "-m",
                "src.stages.stage3_type_plugin",
                "--sample",
                target_sample,
                "--sc_expr_source",
                "normalized",
            ],
            log_root / f"{pair_id}.stage3.log",
        )
        _mirror_stage3_result_for_stage4(root, target_sample)
        stage3_info = _read_stage3_detection(root, target_sample, target_type)

        baseline_out = _stage4(
            root,
            py,
            target_sample,
            "_baseline_profile300",
            "baseline",
            log_root / f"{pair_id}.baseline.stage4.log",
            args.cells_per_spot,
            args.n_subspots,
        )
        _build_assigned_expression(
            root,
            py,
            target_sample,
            baseline_out,
            log_root / f"{pair_id}.baseline.assigned_expression.log",
        )
        _enrich(
            root,
            py,
            args.rscript,
            baseline_out,
            pair_result / "baseline_profile300",
            spec["tumor_labels"],
            pair_id,
            log_root / f"{pair_id}.baseline.enrich.log",
            args.enrichment_backend,
        )

        route2_out = _stage4(
            root,
            py,
            target_sample,
            "_route2_profile300",
            "route2",
            log_root / f"{pair_id}.route2.stage4.log",
            args.cells_per_spot,
            args.n_subspots,
        )
        _build_assigned_expression(
            root,
            py,
            target_sample,
            route2_out,
            log_root / f"{pair_id}.route2.assigned_expression.log",
        )
        _enrich(
            root,
            py,
            args.rscript,
            route2_out,
            pair_result / "route2_profile300",
            spec["tumor_labels"],
            pair_id,
            log_root / f"{pair_id}.route2.enrich.log",
            args.enrichment_backend,
        )

        if not args.skip_other_methods:
            for method, enrich_name in [
                ("tangram_marker", "enrichment_tangram_marker"),
                ("celltrek", "enrichment_celltrek"),
            ]:
                mapping_dir = pair_result / "stage4_mapping" / method
                _run_other_method(
                    root,
                    py,
                    target_sample,
                    method,
                    mapping_dir,
                    log_root / f"{pair_id}.{method}.run.log",
                )
                cyto_like = pair_result / f"cytospace_like_{method}"
                _build_cytospace_like(
                    root,
                    py,
                    target_sample,
                    mapping_dir,
                    cyto_like,
                    log_root / f"{pair_id}.{method}.cytospace_like.log",
                )
                _enrich(
                    root,
                    py,
                    args.rscript,
                    cyto_like,
                    pair_result / enrich_name,
                    spec["tumor_labels"],
                    pair_id,
                    log_root / f"{pair_id}.{method}.enrich.log",
                    args.enrichment_backend,
                )

        manifest_by_pair[pair_id] = {
            "pair_id": pair_id,
            "source_sample": spec["source_sample"],
            "profile_mask_sample": target_sample,
            "masked_target_type": target_type,
            **stage3_info,
        }
        print(
            f"[OK] {pair_id}: target_detected_by_stage3={stage3_info['target_detected_by_stage3']} "
            f"missing={stage3_info['stage3_missing_types']}"
        )

    result_root.mkdir(parents=True, exist_ok=True)
    manifest = [manifest_by_pair[pair_id] for pair_id in pair_specs if pair_id in manifest_by_pair]
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"[OK] wrote: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

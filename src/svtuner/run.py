from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import yaml

from src.utils.sample_paths import resolve_sample_dir, sample_dir_candidates

TRUTH_FILES = (
    "sim_info.json",
    "sim_truth_query_cell_spot.csv",
    "sim_truth_spot_type_fraction.csv",
)


def _env_for_r_subprocess() -> dict:
    """Return an environment suitable for launching R subprocesses on Windows."""
    env = os.environ.copy()
    if os.name == "nt" and not env.get("CONDA_DLL_SEARCH_MODIFICATION_ENABLE"):
        env["CONDA_DLL_SEARCH_MODIFICATION_ENABLE"] = "1"
    return env


def _resolve_path(project_root: Path, p: str) -> Path:
    path = Path(p)
    if path.is_absolute():
        return path
    return (project_root / path).resolve()


def load_presets(project_root: Path) -> Dict[str, Any]:
    p = project_root / "configs" / "pipeline_presets.yaml"
    if not p.exists():
        return {}
    return yaml.safe_load(p.read_text(encoding="utf-8")) or {}


def load_dataset_yaml(project_root: Path, sample: str) -> Dict[str, Any]:
    cfg_path = project_root / "configs" / "datasets" / f"{sample}.yaml"
    if not cfg_path.exists():
        raise FileNotFoundError(f"Dataset config not found: {cfg_path}")
    return yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}


def prepare_stage1_python(project_root: Path, sim_dir: Path, sample: str) -> None:
    """Prepare Stage1 exports from simulation CSV files without running the R preprocessor."""
    export_dir = project_root / "data" / "processed" / sample / "stage1_preprocess" / "exported"
    export_dir.mkdir(parents=True, exist_ok=True)

    for f in [
        "sc_expression.csv",
        "st_expression.csv",
        "sc_metadata.csv",
        "st_coordinates.csv",
        *TRUTH_FILES,
    ]:
        src = sim_dir / f
        if src.exists():
            shutil.copy2(src, export_dir / f)

    sc = pd.read_csv(export_dir / "sc_expression.csv", sep=None, engine="python")
    sc_t = sc.set_index(sc.columns[0]).T
    sc_t.to_csv(export_dir / "sc_expression_normalized.csv")

    st = pd.read_csv(export_dir / "st_expression.csv", sep=None, engine="python")
    st_t = st.set_index(st.columns[0]).T
    st_t.to_csv(export_dir / "st_expression_normalized.csv")

    hvg_src = sim_dir / "hvg_genes.txt"
    if not hvg_src.exists():
        hvg_src = sim_dir.parent / "hvg_genes.txt"
    hvg_dst = project_root / "data" / "processed" / sample / "stage1_preprocess" / "hvg_genes.txt"
    if hvg_src.exists():
        hvg_dst.parent.mkdir(parents=True, exist_ok=True)
        if hvg_src.resolve() != hvg_dst.resolve():
            shutil.copy2(hvg_src, hvg_dst)
        else:
            print("  [info] hvg_genes.txt already exists in target sample directory; skipped copy")

    print("  [OK] Stage1 (Python) done: copied + transposed")


def copy_truth_to_export(project_root: Path, sim_dir: Path, sample: str) -> None:
    export_dir = project_root / "data" / "processed" / sample / "stage1_preprocess" / "exported"
    export_dir.mkdir(parents=True, exist_ok=True)
    for f in TRUTH_FILES:
        src = sim_dir / f
        if src.exists():
            shutil.copy2(src, export_dir / f)
            print(f"  [OK] Copied {f} to stage1 export")


def run_simgen(project_root: Path, preset: Dict[str, Any]) -> None:
    sg = preset.get("simgen") or {}
    module = sg.get("module")
    args_dict = sg.get("args") or {}
    if not module:
        raise ValueError("preset simgen.module is empty")
    cmd: List[str] = [sys.executable, "-m", module]
    for k, v in args_dict.items():
        cmd.append(f"--{k}")
        cmd.append(str(v))
    print(f"\n[SimGen] {' '.join(cmd)}")
    ret = subprocess.run(cmd, cwd=project_root)
    if ret.returncode != 0:
        raise RuntimeError("SimGen failed")


def run_pipeline(
    project_root: Path,
    sample: str,
    *,
    from_scratch: bool = False,
    skip_simgen: bool = False,
    skip_stage1: bool = False,
    use_python_stage1: bool = False,
    missing_type: Optional[str] = None,
    n_processors: Optional[int] = None,
    n_subspots: Optional[int] = None,
) -> None:
    """Run the maintained CLI pipeline through Stage4 only."""
    project_root = project_root.resolve()
    presets = load_presets(project_root)
    dcfg = load_dataset_yaml(project_root, sample)
    pipe = dcfg.get("pipeline") or {}

    mt = missing_type or pipe.get("missing_type") or "T cells CD8"
    n_proc = n_processors if n_processors is not None else int(pipe.get("n_processors", 1))
    n_sub = n_subspots if n_subspots is not None else int(pipe.get("n_subspots", 500))

    paths = dcfg.get("paths") or {}
    st_expr = paths.get("st_expr")
    if not st_expr:
        raise ValueError("configs/datasets/<sample>.yaml is missing paths.st_expr")
    st_expr_path = _resolve_path(project_root, str(st_expr))
    if st_expr_path.exists():
        sim_dir = st_expr_path.parent
    else:
        sim_dir = resolve_sample_dir(project_root, sample, sim_group="real_brca", must_exist=False)
        if not sim_dir.exists():
            for base in sample_dir_candidates(project_root, sample, sim_group="real_brca"):
                cand = (base / str(st_expr)).resolve()
                if cand.exists():
                    sim_dir = cand.parent
                    break

    if from_scratch and not skip_simgen:
        if sample not in presets or "simgen" not in (presets.get(sample) or {}):
            raise FileNotFoundError(
                f"No simgen preset found for {sample} in configs/pipeline_presets.yaml; "
                "add a preset or run simgen manually first."
            )
        run_simgen(project_root, presets[sample])
    elif from_scratch and skip_simgen:
        print("[WARN] --from-scratch and --skip-simgen were both set; skipping SimGen")

    if not sim_dir.is_dir():
        raise FileNotFoundError(f"Simulation/source directory not found: {sim_dir}")

    if not skip_stage1:
        if use_python_stage1:
            print("\n[Stage1] Python preparation ...")
            prepare_stage1_python(project_root, sim_dir, sample)
        else:
            print("\n[Stage1] R preprocess ...")
            from src.config import load_project_yaml

            proj_cfg = load_project_yaml(project_root)
            rscript_exe = proj_cfg.get("rscript_path") or "Rscript"
            rscript = project_root / "r_scripts" / "stage1_preprocess.R"
            cmd = [rscript_exe, str(rscript), "--sample", sample, "--project_root", str(project_root), "--export_csv"]
            resolved = rscript_exe if os.path.isabs(str(rscript_exe)) else shutil.which(str(rscript_exe))
            print(f"  [info] Rscript resolved to: {resolved or rscript_exe!r}")
            ret = subprocess.run(cmd, cwd=project_root, env=_env_for_r_subprocess())
            if ret.returncode != 0:
                raise RuntimeError("Stage1 (R) failed; retry with --use-python-stage1 if appropriate")
            print("  [OK] Stage1 (R) done")
    else:
        print("\n[Stage1] skipped (--skip-stage1)")

    copy_truth_to_export(project_root, sim_dir, sample)

    print("\n[Stage3] type plugin ...")
    ret = subprocess.run(
        [sys.executable, "-m", "src.stages.stage3_type_plugin", "--sample", sample],
        cwd=project_root,
    )
    if ret.returncode != 0:
        raise RuntimeError("Stage3 failed")

    print("\n[Stage4] CytoSPACE baseline ...")
    cmd_baseline = [
        sys.executable,
        "-m",
        "src.stages.stage4_cytospace",
        "--sample",
        sample,
        "--filter_mode",
        "none",
        "--stage4_suffix",
        "_baseline",
        "--missing_type",
        mt,
        "--n_processors",
        str(n_proc),
        "--n_subspots",
        str(n_sub),
    ]
    ret = subprocess.run(cmd_baseline, cwd=project_root)
    if ret.returncode != 0:
        raise RuntimeError("Stage4 baseline failed")

    print("\n[Stage4] CytoSPACE route2 ...")
    cmd_route2 = [
        sys.executable,
        "-m",
        "src.stages.stage4_cytospace",
        "--sample",
        sample,
        "--filter_mode",
        "plugin_unknown",
        "--filter_scope",
        "unsupported_all",
        "--stage4_suffix",
        "_route2",
        "--cell_type_column",
        "plugin_type",
        "--missing_type",
        "__NO_MISSING__",
        "--n_processors",
        str(n_proc),
        "--n_subspots",
        str(n_sub),
    ]
    ret = subprocess.run(cmd_route2, cwd=project_root)
    if ret.returncode != 0:
        raise RuntimeError("Stage4 route2 failed")

    base_result = project_root / "result" / sample
    baseline_out = base_result / "stage4_cytospace_baseline" / "cytospace_output"
    route2_out = base_result / "stage4_cytospace_route2" / "cytospace_output"
    print("\n[DONE] Pipeline completed. Stage4 mapping outputs:")
    print(f"  baseline: {baseline_out}")
    print(f"  route2:   {route2_out}")

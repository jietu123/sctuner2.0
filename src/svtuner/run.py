"""
核心流水线：Stage1（R 或 Python）→ 复制 sim truth → Stage3 → Stage4 baseline+route2 → Stage5。
"""
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
    """Windows 下为 R 子进程继承环境，并默认打开 conda 建议的 DLL 搜索修正（与 README 一致）。"""
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
        raise FileNotFoundError(f"未找到数据集配置: {cfg_path}")
    return yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}


def prepare_stage1_python(project_root: Path, sim_dir: Path, sample: str) -> None:
    """用 Python 准备 Stage1 导出（无需 R）：复制 sim 数据并转置表达矩阵。"""
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
            print("  [info] hvg_genes.txt 已在目标样本目录，跳过复制")

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
        raise ValueError("preset 中 simgen.module 为空")
    cmd: List[str] = [sys.executable, "-m", module]
    for k, v in args_dict.items():
        cmd.append(f"--{k}")
        cmd.append(str(v))
    print(f"\n[SimGen] {' '.join(cmd)}")
    ret = subprocess.run(cmd, cwd=project_root)
    if ret.returncode != 0:
        raise RuntimeError("SimGen 失败")


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
    project_root = project_root.resolve()
    presets = load_presets(project_root)
    dcfg = load_dataset_yaml(project_root, sample)
    pipe = dcfg.get("pipeline") or {}

    mt = missing_type or pipe.get("missing_type") or "T cells CD8"
    # 默认 1：与 stage4_cytospace 一致，降低 Windows 上多进程 + 大矩阵时的 MemoryError 风险
    n_proc = n_processors if n_processors is not None else int(pipe.get("n_processors", 1))
    n_sub = n_subspots if n_subspots is not None else int(pipe.get("n_subspots", 500))

    paths = dcfg.get("paths") or {}
    st_expr = paths.get("st_expr")
    if not st_expr:
        raise ValueError("configs/datasets/<sample>.yaml 缺少 paths.st_expr")
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
                f"configs/pipeline_presets.yaml 中无 {sample} 的 simgen 配置，无法 --from-scratch。"
                " 请添加 preset 或先手动运行 simgen。"
            )
        run_simgen(project_root, presets[sample])
    elif from_scratch and skip_simgen:
        print("[WARN] --from-scratch 与 --skip-simgen 同时指定，跳过 SimGen。")

    if not sim_dir.is_dir():
        raise FileNotFoundError(f"sim 目录不存在（由 st_expr 推导）: {sim_dir}")

    if not skip_stage1:
        if use_python_stage1:
            print("\n[Stage1] Python 准备 ...")
            prepare_stage1_python(project_root, sim_dir, sample)
        else:
            print("\n[Stage1] R preprocess ...")
            from src.config import load_project_yaml

            proj_cfg = load_project_yaml(project_root)
            rscript_exe = proj_cfg.get("rscript_path") or "Rscript"
            rscript = project_root / "r_scripts" / "stage1_preprocess.R"
            cmd = [rscript_exe, str(rscript), "--sample", sample, "--project_root", str(project_root), "--export_csv"]
            resolved = rscript_exe if os.path.isabs(str(rscript_exe)) else shutil.which(str(rscript_exe))
            print(f"  [info] Rscript 解析为: {resolved or rscript_exe!r}（与 scripts/run_s1_mt_*.py 相同调用方式）")
            ret = subprocess.run(cmd, cwd=project_root, env=_env_for_r_subprocess())
            if ret.returncode != 0:
                raise RuntimeError("Stage1 (R) 失败，可改用 --use-python-stage1")
            print("  [OK] Stage1 (R) done")
    else:
        print("\n[Stage1] 跳过 (--skip-stage1)")

    copy_truth_to_export(project_root, sim_dir, sample)

    print("\n[Stage3] type plugin ...")
    ret = subprocess.run(
        [sys.executable, "-m", "src.stages.stage3_type_plugin", "--sample", sample],
        cwd=project_root,
    )
    if ret.returncode != 0:
        raise RuntimeError("Stage3 失败")

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
        raise RuntimeError("Stage4 baseline 失败")

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
        "missing_only",
        "--stage4_suffix",
        "_route2",
        "--cell_type_column",
        "plugin_type",
        "--missing_type",
        mt,
        "--n_processors",
        str(n_proc),
        "--n_subspots",
        str(n_sub),
    ]
    ret = subprocess.run(cmd_route2, cwd=project_root)
    if ret.returncode != 0:
        raise RuntimeError("Stage4 route2 失败")

    base_result = project_root / "result" / sample
    print("\n[Stage5] evaluation ...")
    for run_tag, suffix in [("baseline", "_baseline"), ("route2", "_route2")]:
        stage4_dir = base_result / f"stage4_cytospace{suffix}" / "cytospace_output"
        out_dir = base_result / f"stage5_eval_{run_tag}"
        cmd = [
            sys.executable,
            "-m",
            "src.stages.stage5_route2_s0",
            "--sample",
            sample,
            "--run_tag",
            run_tag,
            "--stage4_dir",
            str(stage4_dir),
            "--sim_dir",
            str(sim_dir),
            "--out_dir",
            str(out_dir),
        ]
        ret = subprocess.run(cmd, cwd=project_root)
        if ret.returncode != 0:
            raise RuntimeError(f"Stage5 {run_tag} 失败")

    print(f"\n[DONE] Pipeline 完成。结果目录: result/{sample}/")

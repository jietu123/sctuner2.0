"""
Centralized defaults for SVTuner.
If you need to change default R command, paths, or stage selections, do it here.
Also hosts ProjectConfig to读取 project_config.yaml / configs/datasets/*.yaml.
"""
from pathlib import Path
from typing import Dict, Any, List, Optional

import yaml

# Default R command (Windows conda env recommended to avoid DLL issues)
DEFAULT_R_CMD = "conda run -n cytospace_v1.1.0_py310 Rscript"

# Default stage selection (can be overridden via CLI)
# Currently default only Stage1; Stage0 is optional dummy.
DEFAULT_STAGES = "1"

# Project paths (resolved at runtime)
def detect_project_root() -> Path:
    """Return project root assuming this file is under src/."""
    here = Path(__file__).resolve()
    return here.parents[1]


def project_path(*parts: str) -> Path:
    return detect_project_root().joinpath(*parts)


# Common subdirectories (use project_path(...) to resolve to absolute)
DATA_RAW = Path("data/raw")
DATA_PROCESSED = Path("data/processed")
RESULT = Path("result")
LOGS = Path("logs")


class ProjectConfig:
    """Thin helper to load project-level和dataset-level YAML 配置."""

    def __init__(self, project_root: Path):
        self.project_root = project_root
        cfg_path = project_root / "configs" / "project_config.yaml"
        self.project_cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) if cfg_path.exists() else {}
        self._dataset_cfg_map = self.project_cfg.get("dataset_config_map", {}) or {}
        self._dataset_cache: Dict[str, Dict[str, Any]] = {}

    def _dataset_cfg_name(self, sample: str) -> str:
        mapped = self._dataset_cfg_map.get(sample)
        if mapped:
            return str(mapped)
        return f"{sample}.yaml"

    def load_dataset_cfg(self, sample: str) -> Dict[str, Any]:
        if sample in self._dataset_cache:
            return self._dataset_cache[sample]
        cfg_path = self.project_root / "configs" / "datasets" / self._dataset_cfg_name(sample)
        if cfg_path.exists():
            data = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
        else:
            data = {}
        self._dataset_cache[sample] = data or {}
        return self._dataset_cache[sample]

    def dataset_cfg_path(self, sample: str) -> Path:
        return self.project_root / "configs" / "datasets" / self._dataset_cfg_name(sample)

    def get_stage_dir(self, sample: str, stage: str) -> Path:
        root = self.project_root
        if stage == "stage1_preprocess":
            return root / "data" / "processed" / sample / "stage1_preprocess"
        if stage == "stage2_svg_plugin":
            return root / "data" / "processed" / sample / "stage2_svg_plugin"
        if stage == "stage3_typematch":
            return root / "data" / "processed" / sample / "stage3_typematch"
        if stage == "stage4_mapping":
            return root / "result" / sample / "stage4_mapping"
        raise ValueError(f"未知 stage: {stage}")

    def enabled_backends(self, sample: str) -> List[str]:
        mapping_cfg = self.project_cfg.get("mapping", {})
        return mapping_cfg.get("enabled_backends", ["cytospace"])

    def _default_mapping_params(self) -> Dict[str, Any]:
        return {
            "seed": 42,
            "svg_refine_lambda": 0.5,
            "lambda_prior": 1.0,
            "prior_candidate_topk": 0,
            "prior_candidate_weight": 1.0,
            "abstain_unknown_sc_only": False,
            "eps": 1e-8,
            "umi_to_cell_norm": "median",
            "default_cells_per_spot": 1.0,
            "cells_per_spot_source": "auto",  # auto|spot_cell_counts|UMI_total|uniform
            "cells_per_spot_clip_min": 1,
            "cells_per_spot_clip_max": None,
            "cells_per_spot_rounding": "round",  # round|floor|ceil
            "distance_metric": "Pearson_correlation",
            "knn_block_size": 1024,
            "knn_max_dense_n": 5000,
            "knn_metric": "euclidean",
            "harden_topk": 5,
            "type_prior_apply_refine": True,
            "type_prior_apply_harden": True,
            "svg_refine_batch_size": None,
            "min_gene_overlap_ratio": 0.2,
            "max_cells_missing_type_prior_ratio": 0.0,
            "min_prior_row_nonzero_ratio": 0.0,
        }

    def get_mapping_config(self, sample: str, backend: str, config_id: Optional[str] = None) -> Dict[str, Any]:
        """
        加载映射配置，支持 config_id overlay。
        
        Merge 顺序：
        1. defaults
        2. project_config.yaml -> mapping -> {backend}
        3. configs/datasets/{sample}.yaml -> mapping -> {backend}
        4. configs/datasets/{config_id}.yaml -> mapping -> {backend} (如果 config_id 指定)
        5. CLI overrides (在调用处处理)
        
        返回的配置中包含 resolved_config_paths 和 resolved_config_sha1s 用于证据链追溯。
        """
        import hashlib
        
        cfg = self._default_mapping_params()
        resolved_paths = []
        resolved_sha1s = {}
        
        # 1. Project-level config
        proj_map = self.project_cfg.get("mapping", {})
        backend_cfg = proj_map.get(backend, {})
        if backend_cfg:
            proj_cfg_path = self.project_root / "configs" / "project_config.yaml"
            if proj_cfg_path.exists():
                resolved_paths.append(str(proj_cfg_path))
                resolved_sha1s[str(proj_cfg_path)] = hashlib.sha1(proj_cfg_path.read_bytes()).hexdigest()
        cfg.update(backend_cfg or {})
        
        # 2. Dataset-level base config
        ds_cfg = self.load_dataset_cfg(sample)
        ds_base_path = self.dataset_cfg_path(sample)
        if ds_base_path.exists():
            resolved_paths.append(str(ds_base_path))
            resolved_sha1s[str(ds_base_path)] = hashlib.sha1(ds_base_path.read_bytes()).hexdigest()
        
        ds_map = ds_cfg.get("mapping", {})
        if backend in ds_map:
            cfg.update(ds_map[backend] or {})
        
        # 3. Config ID overlay (如果指定)
        if config_id:
            overlay_path = self.project_root / "configs" / "datasets" / f"{config_id}.yaml"
            if not overlay_path.exists():
                raise FileNotFoundError(
                    f"Config ID overlay file not found: {overlay_path}\n"
                    f"如果指定了 --config_id={config_id}，必须存在对应的配置文件。"
                )
            
            resolved_paths.append(str(overlay_path))
            resolved_sha1s[str(overlay_path)] = hashlib.sha1(overlay_path.read_bytes()).hexdigest()
            
            overlay_cfg = yaml.safe_load(overlay_path.read_text(encoding="utf-8")) or {}
            overlay_map = overlay_cfg.get("mapping", {})
            if backend in overlay_map:
                cfg.update(overlay_map[backend] or {})
        
        # 记录配置证据链（必须在最后，确保不被后续更新覆盖）
        cfg["resolved_config_paths"] = resolved_paths
        cfg["resolved_config_sha1s"] = resolved_sha1s
        if config_id:
            cfg["resolved_config_id"] = config_id
        
        return cfg


def load_project_config(project_root: Path) -> ProjectConfig:
    return ProjectConfig(project_root)


__all__ = [
    "DEFAULT_R_CMD",
    "DEFAULT_STAGES",
    "detect_project_root",
    "project_path",
    "DATA_RAW",
    "DATA_PROCESSED",
    "RESULT",
    "LOGS",
    "ProjectConfig",
    "load_project_config",
]

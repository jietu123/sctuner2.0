from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import yaml


PAIR_TO_SOURCE_SAMPLE = {
    "melanoma_slide1": "cytospace_fig2c_melanoma_mel1_rep2",
    "melanoma_slide2": "cytospace_fig2c_melanoma_mel2_rep1",
    "brca_er_her2_fresh_frozen": "cytospace_fig2d_tme_brca_er_her2_fresh_frozen",
    "brca_her2_ffpe": "cytospace_fig2d_tme_brca_her2_ffpe",
    "brca_tnbc_fresh_frozen": "cytospace_fig2d_tme_brca_tnbc_fresh_frozen",
    "crc_fresh_frozen": "cytospace_fig2d_tme_crc_fresh_frozen",
}


def _candidate_processed_dirs(root: Path, sample: str) -> list[Path]:
    return [
        root / "data" / "processed" / sample,
        root / "data" / "processed" / "cytospace_fig2c_melanoma" / sample,
        root / "data" / "processed" / "cytospace_fig2d_tme" / sample,
    ]


def _resolve_processed_dir(root: Path, sample: str) -> Path:
    for path in _candidate_processed_dirs(root, sample):
        if (path / "stage1_preprocess" / "exported").exists():
            return path.resolve()
    raise FileNotFoundError(f"Cannot find Stage1 export for {sample}")


def _safe_rmdir(path: Path) -> None:
    if not path.exists():
        return
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink()


def _make_junction(link: Path, target: Path) -> None:
    link.parent.mkdir(parents=True, exist_ok=True)
    if link.exists():
        return
    subprocess.run(["cmd", "/c", "mklink", "/J", str(link), str(target)], check=True, stdout=subprocess.DEVNULL)


def _rewrite_config(root: Path, source_sample: str, forced_sample: str, unsupported_types: list[str]) -> None:
    src_cfg_path = root / "configs" / "datasets" / f"{source_sample}.yaml"
    cfg = yaml.safe_load(src_cfg_path.read_text(encoding="utf-8")) if src_cfg_path.exists() else {}
    cfg = cfg or {}
    cfg.setdefault("paths", {"sc_expr": "unused", "sc_meta": None, "st_expr": None, "st_meta": None, "svg_marker_whitelist": None})
    cfg.setdefault("qc", {})
    cfg.setdefault("gene_filter", {})
    stage3 = cfg.setdefault("stage3", {})
    stage3["force_unsupported_types"] = unsupported_types
    hvg_path = root / "data" / "processed" / forced_sample / "stage1_preprocess" / "hvg_genes.txt"
    if hvg_path.exists():
        stage3["plugin_genes_path"] = str(hvg_path.relative_to(root)).replace("\\", "/")
    cfg["storage"] = {}
    out = root / "configs" / "datasets" / f"{forced_sample}.yaml"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Create Fig.2d forced-unsupported Stage1 aliases without copying matrices.")
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--unsupported_type", default="B cells")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--out_manifest", default="result/cytospace_fig2d_forced_unsupported/forced_samples_manifest.json")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    unsupported_types = [x.strip() for x in str(args.unsupported_type).split(",") if x.strip()]
    if not unsupported_types:
        raise ValueError("--unsupported_type cannot be empty")

    rows = []
    for pair_id, source_sample in PAIR_TO_SOURCE_SAMPLE.items():
        forced_sample = f"{source_sample}_force_bcells_unsupported"
        source_dir = _resolve_processed_dir(root, source_sample)
        forced_dir = root / "data" / "processed" / forced_sample
        if args.overwrite:
            _safe_rmdir(forced_dir)
        forced_dir.mkdir(parents=True, exist_ok=True)
        stage1_dir = forced_dir / "stage1_preprocess"
        stage1_dir.mkdir(parents=True, exist_ok=True)

        source_stage1 = source_dir / "stage1_preprocess"
        source_export = source_stage1 / "exported"
        forced_export = stage1_dir / "exported"
        if args.overwrite and forced_export.exists():
            _safe_rmdir(forced_export)
        _make_junction(forced_export, source_export)

        for name in ["hvg_genes.txt", "hvg_gene_scores.csv"]:
            src = source_stage1 / name
            if src.exists():
                shutil.copy2(src, stage1_dir / name)

        _rewrite_config(root, source_sample, forced_sample, unsupported_types)
        rows.append(
            {
                "pair_id": pair_id,
                "source_sample": source_sample,
                "forced_sample": forced_sample,
                "processed_dir": str(forced_dir),
                "source_processed_dir": str(source_dir),
                "unsupported_types": unsupported_types,
            }
        )
        print(f"[OK] forced sample: {pair_id} -> {forced_sample}")

    out = root / args.out_manifest
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print(f"[OK] wrote: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

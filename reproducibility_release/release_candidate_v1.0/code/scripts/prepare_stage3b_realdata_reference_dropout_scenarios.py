#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path

import pandas as pd
import yaml


def slug(value: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", str(value).lower())).strip("_")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build real-data Stage3B reference-dropout scenarios from selected target types."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--selected_csv",
        default="visualizations/stage3b_realdata_candidate_scan/fig2d_source_selected_12_targets.csv",
    )
    parser.add_argument("--target_group", default="stage3b_realdata_reference_dropout")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--marker_gene_union_only",
        action="store_true",
        help="write reduced SC/ST matrices containing the union of marker genes in selected_csv; screening only",
    )
    return parser.parse_args()


def copy_if_needed(src: Path, dst: Path, overwrite: bool) -> None:
    if dst.exists() and not overwrite:
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def write_dataset_config(root: Path, sample: str, target_group: str) -> None:
    cfg = {
        "paths": {
            "sc_expr": "stage1_direct",
            "sc_meta": "stage1_direct",
            "st_expr": "stage1_direct",
            "st_meta": "stage1_direct",
            "svg_marker_whitelist": None,
        },
        "qc": {
            "sc_min_genes": 0,
            "sc_max_genes": 100000,
            "sc_max_mt": 100,
            "st_min_genes": 0,
            "st_max_genes": ".inf",
            "st_max_mt": 100,
            "hvg_nfeatures": 2000,
            "mt_pattern": "^(MT-|mt-)",
        },
        "gene_filter": {
            "min_cells_sc": 0,
            "min_cells_st": 0,
        },
        "stage3b": {
            "fdr": 0.05,
            "n_calibration": 0,
            "n_spatial_permutations": 200,
            "random_seed": 42,
            "sc_expr_source": "normalized",
            "sc_profile_source": "normalized",
            "expression_scale": "log1p",
            "sc_profile_scale": "log1p",
            "max_genes": 0,
        },
        "storage": {
            "group": target_group,
        },
    }
    cfg_path = root / "configs" / "datasets" / f"{sample}.yaml"
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    selected = pd.read_csv(root / args.selected_csv)
    marker_gene_union_by_source: dict[str, set[str]] = {}
    if args.marker_gene_union_only:
        for _, marker_row in selected.iterrows():
            source_sample = str(marker_row["source_sample"])
            genes = {
                gene
                for gene in str(marker_row.get("marker_genes_list", "")).split(";")
                if gene
            }
            marker_gene_union_by_source.setdefault(source_sample, set()).update(genes)
    rows: list[dict[str, object]] = []
    for _, row in selected.iterrows():
        pair_id = str(row["pair_id"])
        source_sample = str(row["source_sample"])
        source_group = str(row["storage_group"])
        target_type = str(row["cell_type"])
        sample = f"{source_sample}_sc_missing_{slug(target_type)}"
        print(f"[prepare] {sample}: source={source_sample} drop={target_type}", flush=True)

        source_export = (
            root
            / "data"
            / "processed"
            / source_group
            / source_sample
            / "stage1_preprocess"
            / "exported"
        )
        target_stage1 = (
            root
            / "data"
            / "processed"
            / args.target_group
            / sample
            / "stage1_preprocess"
        )
        target_export = target_stage1 / "exported"
        if target_export.exists() and args.overwrite:
            shutil.rmtree(target_stage1)
        target_export.mkdir(parents=True, exist_ok=True)

        sc_meta = pd.read_csv(source_export / "sc_metadata.csv", index_col=0)
        sc_meta.index = sc_meta.index.astype(str)
        if "cell_type" not in sc_meta.columns:
            raise ValueError(f"sc_metadata.csv lacks cell_type: {source_sample}")
        keep = ~sc_meta["cell_type"].astype(str).eq(target_type)
        removed = int((~keep).sum())
        if removed <= 0:
            raise ValueError(f"No SC cells removed for {sample}: {target_type}")
        sc_meta_out = sc_meta.loc[keep].copy()
        sc_meta_out.to_csv(target_export / "sc_metadata.csv", index_label="cell_id")

        if args.marker_gene_union_only:
            wanted = marker_gene_union_by_source.get(source_sample, set())
            if not wanted:
                raise ValueError(f"No marker-gene union available for {source_sample}")
            sc_header = pd.read_csv(source_export / "sc_expression_normalized.csv", nrows=0).columns
            sc_usecols = [sc_header[0], *[g for g in sc_header[1:] if g in wanted]]
            st_header = pd.read_csv(source_export / "st_expression_normalized.csv", nrows=0).columns
            st_usecols = [st_header[0], *[g for g in st_header[1:] if g in wanted]]
            common = sorted(set(sc_usecols[1:]).intersection(st_usecols[1:]))
            if len(common) < 5:
                raise ValueError(
                    f"Marker-gene screening subset too small for {source_sample}: {len(common)}"
                )
            sc_expr = pd.read_csv(
                source_export / "sc_expression_normalized.csv",
                index_col=0,
                usecols=[sc_header[0], *common],
            )
        else:
            sc_expr = pd.read_csv(source_export / "sc_expression_normalized.csv", index_col=0)
        sc_expr.index = sc_expr.index.astype(str)
        sc_expr.loc[sc_meta_out.index].to_csv(target_export / "sc_expression_normalized.csv")

        if args.marker_gene_union_only:
            st_expr = pd.read_csv(
                source_export / "st_expression_normalized.csv",
                index_col=0,
                usecols=[st_header[0], *common],
            )
            st_expr.index = st_expr.index.astype(str)
            st_expr.to_csv(target_export / "st_expression_normalized.csv")
        else:
            copy_if_needed(
                source_export / "st_expression_normalized.csv",
                target_export / "st_expression_normalized.csv",
                args.overwrite,
            )
        copy_if_needed(
            source_export / "st_coordinates.csv",
            target_export / "st_coordinates.csv",
            args.overwrite,
        )
        hvg = source_export.parent / "hvg_genes.txt"
        if hvg.exists():
            copy_if_needed(hvg, target_stage1 / "hvg_genes.txt", args.overwrite)

        info = {
            "sample": sample,
            "source_sample": source_sample,
            "source_group": source_group,
            "target_type": target_type,
            "sc_cells_removed": removed,
            "sc_cells_remaining": int(len(sc_meta_out)),
            "st_spots": int(row.get("st_spots", -1)),
            "marker_genes": str(row.get("marker_genes_list", "")).split(";"),
            "pair_id": pair_id,
            "marker_gene_union_only": bool(args.marker_gene_union_only),
        }
        (target_stage1 / "reference_dropout_info.json").write_text(
            json.dumps(info, indent=2),
            encoding="utf-8",
        )
        write_dataset_config(root, sample, args.target_group)
        rows.append({**info, "storage_group": args.target_group})

    manifest = pd.DataFrame(rows)
    out_path = (
        root
        / "visualizations"
        / "stage3b_realdata_candidate_scan"
        / "stage3b_realdata_reference_dropout_manifest.csv"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(out_path, index=False)
    print(f"[done] {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

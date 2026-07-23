#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rank real reference-dropout cases by blank regions induced specifically "
            "by removing the target SC type."
        )
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--candidate_csv",
        default=(
            "visualizations/stage3b_realdata_candidate_scan/"
            "stage3b_reference_dropout_13_full_matrix_validation.csv"
        ),
    )
    parser.add_argument("--target_quantile", type=float, default=0.85)
    parser.add_argument(
        "--marker_metadata_csv",
        default=(
            "visualizations/stage3b_realdata_candidate_scan/spatial_9x2/"
            "stage3b_reference_dropout_spatial_stack_recommended_6x2_metadata.csv"
        ),
    )
    parser.add_argument(
        "--out_csv",
        default=(
            "visualizations/stage3b_realdata_candidate_scan/"
            "stage3b_reference_dropout_causal_negative_control_screen.csv"
        ),
    )
    return parser.parse_args()


def read_yaml(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def processed_dir(root: Path, sample: str) -> Path:
    config = read_yaml(root / "configs" / "datasets" / f"{sample}.yaml")
    group = str((config.get("storage") or {}).get("group") or "").strip("/\\")
    base = root / "data" / "processed"
    return base / group / sample if group else base / sample


def read_blank(path: Path) -> pd.Series:
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    values = frame["is_unsupported_region"]
    if values.dtype != bool:
        values = values.astype(str).str.lower().isin(["true", "1", "yes"])
    return values.astype(bool)


def marker_percentile(path: Path, marker_genes: list[str]) -> pd.Series:
    header = pd.read_csv(path, nrows=0).columns.tolist()
    id_col = header[0]
    present = [gene for gene in marker_genes if gene in header]
    if not present:
        raise ValueError(f"No marker genes found in {path}")
    expression = pd.read_csv(path, index_col=0, usecols=[id_col, *present])
    expression.index = expression.index.astype(str)
    expression = expression.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return expression.mean(axis=1).rank(method="average", pct=True)


def infer_marker_genes(
    root: Path,
    source_sample: str,
    target_type: str,
    cache: dict[str, dict[str, object]],
    n_genes: int = 30,
) -> list[str]:
    if source_sample not in cache:
        export = (
            processed_dir(root, source_sample)
            / "stage1_preprocess"
            / "exported"
        )
        metadata = pd.read_csv(export / "sc_metadata.csv", index_col=0)
        metadata.index = metadata.index.astype(str)
        type_column = (
            "cell_type"
            if "cell_type" in metadata.columns
            else "sc_meta"
            if "sc_meta" in metadata.columns
            else None
        )
        if type_column is None:
            raise ValueError(f"No cell-type column in {export / 'sc_metadata.csv'}")
        metadata[type_column] = metadata[type_column].astype(str)
        expression = pd.read_csv(
            export / "sc_expression_normalized.csv",
            index_col=0,
        )
        expression.index = expression.index.astype(str)
        common = expression.index.intersection(metadata.index)
        expression = expression.loc[common]
        cell_types = metadata.loc[common, type_column]
        cache[source_sample] = {
            "grouped": expression.groupby(cell_types).mean(),
            "global_sum": expression.sum(axis=0),
            "counts": cell_types.value_counts(),
            "n_total": len(expression),
        }
    source = cache[source_sample]
    grouped = source["grouped"]
    counts = source["counts"]
    if target_type not in grouped.index:
        raise KeyError(
            f"Target type {target_type!r} not found in {source_sample}; "
            f"available={list(grouped.index)}"
        )
    count = int(counts[target_type])
    type_mean = grouped.loc[target_type]
    other_mean = (
        source["global_sum"] - type_mean * count
    ) / max(int(source["n_total"]) - count, 1)
    specificity = (type_mean + 1e-6) / (other_mean + 1e-6)
    table = pd.DataFrame(
        {
            "gene": type_mean.index,
            "type_mean": type_mean.to_numpy(dtype=float),
            "specificity": specificity.to_numpy(dtype=float),
        }
    )
    table = table[table["type_mean"] > 0]
    table = table.sort_values(
        ["specificity", "type_mean"],
        ascending=False,
    )
    return table["gene"].head(n_genes).tolist()


def safe_ratio(numerator: int, denominator: int) -> float:
    return float(numerator / denominator) if denominator else 0.0


def evaluate(
    root: Path,
    row: pd.Series,
    target_quantile: float,
    marker_cache: dict[str, dict[str, object]],
) -> dict[str, object]:
    source_sample = str(row["source_sample"])
    dropout_sample = str(row["sample"])
    source_processed = processed_dir(root, source_sample)
    dropout_processed = processed_dir(root, dropout_sample)
    source_scores_path = (
        source_processed
        / "stage3b_st_unsupported"
        / "spot_unsupported_scores.csv"
    )
    dropout_scores_path = (
        dropout_processed
        / "stage3b_st_unsupported"
        / "spot_unsupported_scores.csv"
    )
    if not source_scores_path.exists():
        raise FileNotFoundError(f"Missing full-reference control: {source_scores_path}")
    if not dropout_scores_path.exists():
        raise FileNotFoundError(f"Missing dropout scores: {dropout_scores_path}")

    if "marker_genes_list" in row.index and pd.notna(row["marker_genes_list"]):
        marker_genes = [
            gene for gene in str(row["marker_genes_list"]).split(";") if gene
        ]
    else:
        marker_genes = infer_marker_genes(
            root,
            source_sample,
            str(row["target_type"]),
            marker_cache,
            n_genes=int(row.get("marker_genes", 30)),
        )
    expression_path = (
        source_processed
        / "stage1_preprocess"
        / "exported"
        / "st_expression_normalized.csv"
    )
    percentile = marker_percentile(expression_path, marker_genes)
    control = read_blank(source_scores_path)
    dropout = read_blank(dropout_scores_path)
    common = percentile.index.intersection(control.index).intersection(dropout.index)
    percentile = percentile.loc[common]
    control = control.loc[common]
    dropout = dropout.loc[common]

    target = percentile >= percentile.quantile(target_quantile)
    induced = dropout & ~control
    lost = control & ~dropout
    persistent = control & dropout
    control_overlap = control & target
    dropout_overlap = dropout & target
    induced_overlap = induced & target
    lost_overlap = lost & target

    target_count = int(target.sum())
    control_count = int(control.sum())
    dropout_count = int(dropout.sum())
    induced_count = int(induced.sum())
    lost_count = int(lost.sum())
    persistent_count = int(persistent.sum())
    control_target = int(control_overlap.sum())
    dropout_target = int(dropout_overlap.sum())
    induced_target = int(induced_overlap.sum())
    lost_target = int(lost_overlap.sum())
    net_target_gain = dropout_target - control_target
    net_blank_gain = dropout_count - control_count
    target_remaining_after_control = max(target_count - control_target, 0)
    background_count = max(len(common) - target_count, 0)
    induced_background = induced_count - induced_target

    induced_precision = safe_ratio(induced_target, induced_count)
    induced_recall_total = safe_ratio(induced_target, target_count)
    induced_recall_residual = safe_ratio(
        induced_target,
        target_remaining_after_control,
    )
    control_recall = safe_ratio(control_target, target_count)
    dropout_recall = safe_ratio(dropout_target, target_count)
    dropout_precision = safe_ratio(dropout_target, dropout_count)
    background_false_positive_rate = safe_ratio(
        induced_background,
        background_count,
    )
    causal_score = (
        0.35 * induced_precision
        + 0.35 * induced_recall_residual
        + 0.20 * max(dropout_recall - control_recall, 0.0)
        + 0.10 * min(induced_count / max(target_count, 1), 1.0)
        - 0.20 * background_false_positive_rate
    )

    return {
        "pair_id": row["pair_id"],
        "source_sample": source_sample,
        "sample": dropout_sample,
        "target_type": row["target_type"],
        "spots": int(len(common)),
        "target_spots": target_count,
        "control_blank_spots": control_count,
        "dropout_blank_spots": dropout_count,
        "persistent_blank_spots": persistent_count,
        "dropout_induced_blank_spots": induced_count,
        "dropout_lost_blank_spots": lost_count,
        "net_blank_gain": net_blank_gain,
        "control_target_overlap": control_target,
        "dropout_target_overlap": dropout_target,
        "dropout_induced_target_overlap": induced_target,
        "dropout_lost_target_overlap": lost_target,
        "net_target_overlap_gain": net_target_gain,
        "control_target_recall": control_recall,
        "dropout_target_recall": dropout_recall,
        "dropout_target_precision": dropout_precision,
        "induced_target_precision": induced_precision,
        "induced_target_recall_total": induced_recall_total,
        "induced_target_recall_after_control": induced_recall_residual,
        "induced_background_spots": induced_background,
        "induced_background_false_positive_rate": background_false_positive_rate,
        "blank_mask_jaccard": safe_ratio(
            persistent_count,
            int((control | dropout).sum()),
        ),
        "causal_selection_score": causal_score,
        "suitable_primary_case": bool(
            control_recall <= 0.20
            and induced_count >= 20
            and induced_precision >= 0.70
            and induced_recall_residual >= 0.30
            and net_target_gain >= 20
        ),
    }


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    candidates = pd.read_csv(root / args.candidate_csv)
    marker_metadata_path = root / args.marker_metadata_csv
    if marker_metadata_path.exists():
        marker_metadata = pd.read_csv(marker_metadata_path)[
            ["sample", "marker_genes_list"]
        ].drop_duplicates("sample")
        candidates = candidates.merge(
            marker_metadata,
            on="sample",
            how="left",
        )
    marker_cache: dict[str, dict[str, object]] = {}
    rows = [
        evaluate(root, row, args.target_quantile, marker_cache)
        for _, row in candidates.iterrows()
    ]
    result = pd.DataFrame(rows).sort_values(
        [
            "suitable_primary_case",
            "causal_selection_score",
            "induced_target_precision",
            "net_target_overlap_gain",
        ],
        ascending=False,
    )
    out_path = root / args.out_csv
    out_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_path, index=False)
    columns = [
        "pair_id",
        "target_type",
        "control_blank_spots",
        "dropout_induced_blank_spots",
        "net_target_overlap_gain",
        "control_target_recall",
        "induced_target_precision",
        "induced_target_recall_after_control",
        "causal_selection_score",
        "suitable_primary_case",
    ]
    print(result[columns].to_string(index=False))
    print(f"[done] {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CURVE_COLOR = "#edae1a"
HIT_COLOR = "#f28e2b"
TRACK_COLOR = "#d8d8d8"
TRACK_EDGE = "#bcbcbc"
METHODS = [("CytoSPACE", "baseline"), ("SVTuner", "route2")]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _processed_dir(root: Path, sample: str) -> Path:
    candidates = [
        root / "data" / "processed" / sample,
        root / "data" / "processed" / "highres_profile_mask" / sample,
    ]
    for path in candidates:
        if (path / "stage1_preprocess" / "exported").exists():
            return path
    raise FileNotFoundError(f"Cannot resolve processed directory for {sample}")


def _read_marker_genes(processed: Path, limit: int) -> list[str]:
    panel_path = processed / "stage1_preprocess" / "fig2d_profile_mask_gene_panel.csv"
    panel = pd.read_csv(panel_path)
    expression_path = (
        processed / "stage1_preprocess" / "exported" / "st_expression_normalized.csv"
    )
    available = set(pd.read_csv(expression_path, nrows=0).columns.astype(str))
    genes: list[str] = []
    for gene in panel["gene"].astype(str):
        if gene in available and gene not in genes:
            genes.append(gene)
        if len(genes) >= limit:
            break
    if not genes:
        raise ValueError(f"No usable profile-mask marker genes in {processed}")
    return genes


def _masked_support(processed: Path, genes: list[str]) -> pd.DataFrame:
    expression = pd.read_csv(
        processed / "stage1_preprocess" / "exported" / "st_expression_normalized.csv",
        usecols=["spot_id", *genes],
    )
    expression["spot_id"] = expression["spot_id"].astype(str)
    expression["masked_support"] = expression[genes].mean(axis=1)
    support = expression[["spot_id", "masked_support"]].copy()
    low = float(support["masked_support"].min())
    high = float(support["masked_support"].max())
    support["masked_support_norm"] = (
        support["masked_support"] - low
    ) / (high - low + 1.0e-12)
    return support


def _cell_type_marker_scores(processed: Path, genes: list[str]) -> dict[str, float]:
    export = processed / "stage1_preprocess" / "exported"
    metadata = pd.read_csv(export / "sc_metadata.csv", usecols=["cell_id", "cell_type"])
    expression = pd.read_csv(
        export / "sc_expression_normalized.csv", usecols=["cell_id", *genes]
    )
    expression["target_marker_score"] = expression[genes].mean(axis=1)
    joined = expression[["cell_id", "target_marker_score"]].merge(
        metadata, on="cell_id", how="inner"
    )
    return joined.groupby("cell_type")["target_marker_score"].mean().to_dict()


def _assignment_path(root: Path, sample: str, suffix: str) -> Path:
    return (
        root
        / "result"
        / sample
        / f"stage4_cytospace_{suffix}_highres"
        / "cytospace_output"
        / "cell_type_assignments_by_spot.csv"
    )


def _predicted_target_like_score(
    root: Path,
    sample: str,
    suffix: str,
    type_scores: dict[str, float],
) -> pd.DataFrame:
    assignments = pd.read_csv(_assignment_path(root, sample, suffix))
    assignments = assignments.rename(columns={assignments.columns[0]: "spot_id"})
    assignments["spot_id"] = assignments["spot_id"].astype(str)

    weighted_score = np.zeros(len(assignments), dtype=float)
    assigned_count = np.zeros(len(assignments), dtype=float)
    for cell_type, marker_score in type_scores.items():
        if cell_type not in assignments.columns:
            continue
        counts = pd.to_numeric(assignments[cell_type], errors="coerce").fillna(0.0)
        values = counts.to_numpy(dtype=float)
        weighted_score += values * float(marker_score)
        assigned_count += values

    if "Total cells" in assignments.columns:
        reported = pd.to_numeric(assignments["Total cells"], errors="coerce").to_numpy(float)
        assigned_count = np.where(np.isfinite(reported), reported, assigned_count)
    score = weighted_score / np.maximum(assigned_count, 1.0)
    return assignments[["spot_id"]].assign(predicted_target_like_score=score)


def _running_enrichment(
    scores: pd.Series,
    fraction: float,
) -> tuple[np.ndarray, np.ndarray, float, float, int, int]:
    values = scores.to_numpy(dtype=float)
    n_total = len(values)
    target_mass = float(n_total) * fraction
    cutoff_index = min(max(int(np.ceil(target_mass)) - 1, 0), n_total - 1)
    threshold = float(np.sort(values, kind="mergesort")[cutoff_index])
    strict_hits = values < threshold
    tied = values == threshold
    remaining_mass = max(target_mass - float(strict_hits.sum()), 0.0)
    tie_weight = remaining_mass / max(int(tied.sum()), 1)
    hit_weights = strict_hits.astype(float) + tied.astype(float) * tie_weight
    hit_mass = float(hit_weights.sum())
    miss_mass = float(n_total) - hit_mass
    if hit_mass <= 0.0 or miss_mass <= 0.0:
        return (
            np.zeros(n_total, dtype=float),
            hit_weights,
            0.0,
            threshold,
            int(strict_hits.sum()),
            int(tied.sum()),
        )
    increments = hit_weights / hit_mass - (1.0 - hit_weights) / miss_mass
    running = np.cumsum(increments)
    return (
        running,
        hit_weights,
        float(running.max()),
        threshold,
        int(strict_hits.sum()),
        int(tied.sum()),
    )


def _build_spot_table(
    root: Path,
    row: pd.Series,
    marker_gene_limit: int,
    low_score_fraction: float,
) -> tuple[pd.DataFrame, dict[str, object]]:
    sample = str(row["profile_mask_sample"])
    processed = _processed_dir(root, sample)
    genes = _read_marker_genes(processed, marker_gene_limit)
    support = _masked_support(processed, genes)
    type_scores = _cell_type_marker_scores(processed, genes)

    spot_table = support
    for _, suffix in METHODS:
        method_score = _predicted_target_like_score(root, sample, suffix, type_scores)
        method_score = method_score.rename(
            columns={"predicted_target_like_score": f"{suffix}_target_like_score"}
        )
        spot_table = spot_table.merge(method_score, on="spot_id", how="inner")

    spot_table = spot_table.sort_values("masked_support_norm", kind="mergesort").reset_index(
        drop=True
    )
    spot_table["rank"] = np.arange(1, len(spot_table) + 1)
    metrics: dict[str, object] = {
        "raw_sample": str(row["raw_sample"]),
        "profile_mask_sample": sample,
        "masked_target_type": str(row["masked_target_type"]),
        "n_spatial_units": int(len(spot_table)),
        "n_marker_genes": int(len(genes)),
        "marker_genes": ";".join(genes),
        "support_sd": float(spot_table["masked_support_norm"].std(ddof=0)),
    }
    for method, suffix in METHODS:
        (
            running,
            hit_weights,
            peak_es,
            threshold,
            n_strict_hits,
            n_threshold_ties,
        ) = _running_enrichment(
            spot_table[f"{suffix}_target_like_score"],
            low_score_fraction,
        )
        spot_table[f"{suffix}_low_score_hit_weight"] = hit_weights
        spot_table[f"{suffix}_running_es"] = running
        metrics[f"{suffix}_peak_es"] = peak_es
        metrics[f"{suffix}_low_score_threshold"] = threshold
        metrics[f"{suffix}_effective_hit_mass"] = float(hit_weights.sum())
        metrics[f"{suffix}_n_strict_hits"] = n_strict_hits
        metrics[f"{suffix}_n_threshold_ties"] = n_threshold_ties
        metrics[f"{suffix}_method_label"] = method
    metrics["delta_peak_es"] = float(metrics["route2_peak_es"]) - float(
        metrics["baseline_peak_es"]
    )
    return spot_table, metrics


def _smooth(values: np.ndarray) -> np.ndarray:
    if len(values) < 50:
        return values
    window = max(5, int(round(len(values) * 0.006)))
    return pd.Series(values).rolling(window, center=True, min_periods=1).mean().to_numpy()


def _draw_method_row(
    fig: plt.Figure,
    cell,
    spot_table: pd.DataFrame,
    method_label: str,
    suffix: str,
    peak_es: float,
) -> None:
    grid = cell.subgridspec(3, 1, height_ratios=[2.35, 0.42, 0.34], hspace=0.025)
    curve_ax = fig.add_subplot(grid[0, 0])
    hit_ax = fig.add_subplot(grid[1, 0], sharex=curve_ax)
    track_ax = fig.add_subplot(grid[2, 0], sharex=curve_ax)

    x = spot_table["rank"].to_numpy(dtype=float)
    running = spot_table[f"{suffix}_running_es"].to_numpy(dtype=float)
    hit_weights = spot_table[f"{suffix}_low_score_hit_weight"].to_numpy(dtype=float)
    curve_ax.plot(x, _smooth(running), color=CURVE_COLOR, linewidth=1.75)
    curve_ax.axhline(0, color="#999999", linestyle="--", linewidth=0.7)
    curve_ax.set_xlim(1, len(spot_table))
    curve_ax.set_yticks([0.0])
    curve_ax.set_yticklabels(["0"], fontsize=7.4)
    curve_ax.set_ylabel("Enrichment", fontsize=7.2)
    curve_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    curve_ax.spines[["top", "right"]].set_visible(False)
    curve_ax.text(
        0.0,
        1.02,
        method_label,
        transform=curve_ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.8,
        color="#555555",
    )
    curve_ax.text(
        0.965,
        0.90,
        f"Peak ES = {peak_es:.3f}",
        transform=curve_ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.0,
    )

    strict_hits = hit_weights >= 1.0 - 1.0e-12
    fractional_hits = (hit_weights > 0.0) & ~strict_hits
    hit_ax.vlines(x[strict_hits], 0.05, 0.95, color=HIT_COLOR, linewidth=0.62, alpha=0.95)
    if fractional_hits.any():
        fractional_alpha = min(0.65, max(0.10, float(hit_weights[fractional_hits][0])))
        hit_ax.vlines(
            x[fractional_hits],
            0.05,
            0.95,
            color=HIT_COLOR,
            linewidth=0.48,
            alpha=fractional_alpha,
        )
    hit_ax.set_ylim(0, 1)
    hit_ax.axis("off")

    support = spot_table["masked_support_norm"].to_numpy(dtype=float)
    track_ax.fill_between(x, 0, support, color=TRACK_COLOR, linewidth=0)
    track_ax.plot(x, support, color=TRACK_EDGE, linewidth=0.7)
    track_ax.set_ylim(0, 1)
    track_ax.set_yticks([])
    track_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    track_ax.spines[["left", "right", "top"]].set_visible(False)


def _plot(
    spot_table: pd.DataFrame,
    selected_metrics: dict[str, object],
    out_dir: Path,
    prefix: str,
) -> dict[str, Path]:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig = plt.figure(figsize=(6.15, 3.65), dpi=320)
    outer = fig.add_gridspec(2, 1, hspace=0.22)
    _draw_method_row(
        fig,
        outer[0, 0],
        spot_table,
        "CytoSPACE",
        "baseline",
        float(selected_metrics["baseline_peak_es"]),
    )
    _draw_method_row(
        fig,
        outer[1, 0],
        spot_table,
        "SVTuner",
        "route2",
        float(selected_metrics["route2_peak_es"]),
    )

    fig.text(0.012, 0.975, "b", ha="left", va="top", fontsize=11.5, weight="bold")
    fig.text(
        0.055,
        0.975,
        "Spatial enrichment of mapped target-like signal",
        ha="left",
        va="top",
        fontsize=8.5,
        weight="bold",
    )
    sample_label = (
        str(selected_metrics["raw_sample"])
        .replace("Human", "Human ")
        .replace("Melanoma", "melanoma ")
        .replace("Patient", "patient ")
    )
    fig.text(
        0.055,
        0.949,
        f"{sample_label} | masked {selected_metrics['masked_target_type']}",
        ha="left",
        va="top",
        fontsize=6.9,
        color="#444444",
    )
    fig.text(0.075, 0.035, "Lower masked support", ha="left", va="bottom", fontsize=7.2)
    fig.text(0.50, 0.035, "Relative support rank", ha="center", va="bottom", fontsize=7.4)
    fig.text(0.93, 0.035, "Higher masked support", ha="right", va="bottom", fontsize=7.2)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.89, bottom=0.09)

    outputs = {
        extension: out_dir / f"{prefix}.{extension}"
        for extension in ("png", "svg", "pdf")
    }
    fig.savefig(outputs["png"], dpi=320)
    fig.savefig(outputs["svg"])
    fig.savefig(outputs["pdf"])
    plt.close(fig)
    return outputs


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Recompute the high-resolution profile-mask enrichment experiment for Fig. 3B."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--config", default="configs/fig3b_highres_profile_mask_enrichment.json"
    )
    parser.add_argument(
        "--out_dir", default="visualizations/highres_profile_mask_fig3b_enrichment"
    )
    parser.add_argument("--out_prefix", default="fig3b_highres_profile_mask_enrichment")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    config_path = root / args.config
    config = json.loads(config_path.read_text(encoding="utf-8"))
    manifest_path = (
        root
        / "visualizations"
        / "highres_profile_mask_mapping"
        / "highres_profile_mask_mapping_manifest.csv"
    )
    mapping_manifest = pd.read_csv(manifest_path)
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    source_tables: list[pd.DataFrame] = []
    metric_rows: list[dict[str, object]] = []
    experiment_input_paths: list[Path] = [
        Path(__file__).resolve(),
        config_path,
        manifest_path,
    ]
    for _, row in mapping_manifest.iterrows():
        sample = str(row["profile_mask_sample"])
        processed = _processed_dir(root, sample)
        export = processed / "stage1_preprocess" / "exported"
        experiment_input_paths.extend(
            [
                processed / "stage1_preprocess" / "fig2d_profile_mask_gene_panel.csv",
                export / "st_expression_normalized.csv",
                export / "sc_expression_normalized.csv",
                export / "sc_metadata.csv",
                _assignment_path(root, sample, "baseline"),
                _assignment_path(root, sample, "route2"),
            ]
        )
        spot_table, metrics = _build_spot_table(
            root,
            row,
            int(config["marker_gene_limit"]),
            float(config["low_score_fraction"]),
        )
        annotated = spot_table.copy()
        annotated.insert(0, "masked_target_type", metrics["masked_target_type"])
        annotated.insert(0, "raw_sample", metrics["raw_sample"])
        annotated.insert(0, "profile_mask_sample", metrics["profile_mask_sample"])
        source_tables.append(annotated)
        metric_rows.append(metrics)

    metrics_df = pd.DataFrame(metric_rows)
    all_source = pd.concat(source_tables, ignore_index=True)
    representative = str(config["representative_profile_mask_sample"])
    matches = metrics_df[metrics_df["profile_mask_sample"].astype(str).eq(representative)]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one representative row for {representative}, found {len(matches)}"
        )
    selected_metrics = matches.iloc[0].to_dict()
    selected_source = all_source[
        all_source["profile_mask_sample"].astype(str).eq(representative)
    ].copy()

    metrics_path = out_dir / f"{args.out_prefix}_metrics.csv"
    source_path = out_dir / f"{args.out_prefix}_source_values.csv"
    selected_path = out_dir / f"{args.out_prefix}_selected_source_values.csv"
    notes_path = out_dir / "Fig3B_实验与数据来源详细说明.md"
    metrics_df.to_csv(metrics_path, index=False)
    all_source.to_csv(source_path, index=False)
    selected_source.to_csv(selected_path, index=False)
    outputs = _plot(selected_source, selected_metrics, out_dir, args.out_prefix)

    output_manifest = {
        "experiment_id": config["experiment_id"],
        "experiment_definition": (
            "For each cell-resolution profile-mask dataset, spatial units are ranked from "
            "low to high residual expression of the first 30 available masked-target marker "
            "genes. A mapped target-like score is the assigned-cell-count-weighted mean "
            "reference marker score at each spatial unit. Hits carry exactly 10% effective "
            "mass within each method. Boundary ties receive equal fractional hit weight, "
            "avoiding arbitrary score-tie breaking. Peak ES is the maximum unweighted "
            "running enrichment score along the residual-support rank."
        ),
        "interpretation": (
            "Higher positive Peak ES indicates that low mapped target-like scores are "
            "concentrated toward spatial units with lower residual masked-target support. "
            "Peak ES is not NES and is not an accuracy measure against biological ground truth."
        ),
        "representative_profile_mask_sample": representative,
        "representative_selection_rule": config["representative_selection_rule"],
        "n_datasets": int(len(metrics_df)),
        "marker_gene_limit": int(config["marker_gene_limit"]),
        "low_score_fraction": float(config["low_score_fraction"]),
        "score_tie_handling": config["score_tie_handling"],
        "support_tie_handling": config["support_tie_handling"],
        "representative_values": {
            "CytoSPACE_peak_es": float(selected_metrics["baseline_peak_es"]),
            "SVTuner_peak_es": float(selected_metrics["route2_peak_es"]),
            "delta_peak_es": float(selected_metrics["delta_peak_es"]),
        },
        "inputs": {
            "config": str(config_path.relative_to(root)).replace("\\", "/"),
            "mapping_manifest": str(manifest_path.relative_to(root)).replace("\\", "/"),
        },
        "outputs": {
            "png": str(outputs["png"].relative_to(root)).replace("\\", "/"),
            "svg": str(outputs["svg"].relative_to(root)).replace("\\", "/"),
            "pdf": str(outputs["pdf"].relative_to(root)).replace("\\", "/"),
            "metrics_csv": str(metrics_path.relative_to(root)).replace("\\", "/"),
            "source_values_csv": str(source_path.relative_to(root)).replace("\\", "/"),
            "selected_source_values_csv": str(selected_path.relative_to(root)).replace(
                "\\", "/"
            ),
            "detailed_notes": str(notes_path.relative_to(root)).replace("\\", "/"),
        },
        "input_sha256": {
            str(path.relative_to(root)).replace("\\", "/"): _sha256(path)
            for path in experiment_input_paths
        },
        "output_sha256": {
            str(path.relative_to(root)).replace("\\", "/"): _sha256(path)
            for path in [
                outputs["png"],
                outputs["svg"],
                outputs["pdf"],
                metrics_path,
                source_path,
                selected_path,
                notes_path,
            ]
        },
    }
    manifest_out = out_dir / f"{args.out_prefix}_manifest.json"
    manifest_out.write_text(
        json.dumps(output_manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(output_manifest, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

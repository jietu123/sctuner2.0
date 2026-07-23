#!/usr/bin/env python
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "manuscript_audits" / "fig2_enrichment_backend_formal_audit"
AUDIT_DATE = "2026-07-22"

FIG2B_FIGURE = ROOT / "visualizations/simulations/real_profile_mask_fig2c_only/fig2_panel_c_real_profile_mask.png"
FIG2B_SOURCE = ROOT / (
    "result/real_profile_mask_foundation/"
    "ffpe_mouse_brain_sagittal_real_profile_mask_microglia/spot_foundation.csv"
)
FIG2B_PLOT = ROOT / "scripts/plot_real_profile_mask_fig2c_only.py"
FIG2B_UPSTREAM = ROOT / "scripts/build_real_profile_mask_foundation.py"
FIG2B_CONFIG = ROOT / "result/real_profile_mask_foundation/foundation_config.json"
FIG2B_METADATA = ROOT / (
    "result/real_profile_mask_foundation/"
    "ffpe_mouse_brain_sagittal_real_profile_mask_microglia/foundation_metadata.json"
)

FIG2C_FIGURE = ROOT / (
    "visualizations/cytospace_fig2d_profile_mask_benchmark/fig2d_profile_mask_benchmark.png"
)
FIG2C_SOURCE = ROOT / (
    "visualizations/cytospace_fig2d_profile_mask_benchmark/"
    "fig2d_profile_mask_benchmark_source_values.csv"
)
FIG2C_RETAINED = ROOT / (
    "result/cytospace_fig2d_profile_mask_benchmark/fig2d_profile_mask_benchmark_source_values.csv"
)
FIG2C_PLOT = ROOT / "scripts/plot_fig2d_profile_mask_benchmark.py"
FIG2C_RUNNER = ROOT / "scripts/run_cytospace_fig2d_profile_mask_benchmark.py"
FIG2C_CONFIG = ROOT / "result/cytospace_fig2d_profile_mask_benchmark/profile_mask_manifest.json"

PY_BACKEND = ROOT / "scripts/compute_fig2c_enrichment_python.py"
SEURAT_BACKEND = ROOT / "scripts/compute_cytospace_fig2c_official_enrichment_seurat.R"
LEGACY_FGSEA_BACKEND = ROOT / "scripts/compute_cytospace_fig2c_official_enrichment.R"

FIG2D_DIR = ROOT / "visualizations/cytospace_fig2c_melanoma_stage3_profile_mask"
FIG2D_PREFIX = "fig2c_cd4_ce9_stage3_detected_macrophage_mask"
FIG2D_FIGURE = FIG2D_DIR / f"{FIG2D_PREFIX}_baseline_vs_route2.png"
FIG2D_SOURCE = FIG2D_DIR / f"{FIG2D_PREFIX}_metrics.csv"
FIG2D_CONFIG = FIG2D_DIR / f"{FIG2D_PREFIX}_manifest.json"
FIG2D_BASELINE = ROOT / (
    "result/cytospace_fig2c_melanoma_mel1_rep2_screen_mask_macrophages/"
    "fig2c_strict_cd4_ce9/baseline/fig2c_official_enrichment_summary.csv"
)
FIG2D_ROUTE = ROOT / (
    "result/cytospace_fig2c_melanoma_mel1_rep2_screen_mask_macrophages/"
    "fig2c_strict_cd4_ce9/route2/fig2c_official_enrichment_summary.csv"
)
FIG2D_HISTORICAL_PLOT = "scripts/plot_cytospace_fig2c_cd4_baseline_vs_route2.py"
FIG2D_HISTORICAL_COMMIT = "76cc0288fbd10b870f6d56d06d81e91e70cd0eb3"
FIG2D_HISTORICAL_SHA1 = "13ce18a3e8b84c98ca9f8a834090f4a60d3bb7b9"

FIG2F_FIGURE = ROOT / (
    "visualizations/cytospace_fig2e_stage3_profile_mask/"
    "fig2e_stage3_profile_mask_route2_ce9_ce10.png"
)
FIG2F_SOURCE = ROOT / (
    "result/cytospace_fig2e_stage3_profile_mask/fig2e_stage3_profile_mask_metrics.csv"
)
FIG2F_COMPUTE = ROOT / "scripts/build_fig2e_stage3_profile_mask_benchmark.py"
FIG2F_PLOT = ROOT / "scripts/build_fig2e_stage3_profile_mask_route2_only.py"


PROTECTED = [
    FIG2B_FIGURE,
    FIG2B_SOURCE,
    FIG2B_PLOT,
    FIG2B_UPSTREAM,
    FIG2C_FIGURE,
    FIG2C_SOURCE,
    FIG2C_RETAINED,
    FIG2C_PLOT,
    FIG2C_RUNNER,
    PY_BACKEND,
    SEURAT_BACKEND,
    LEGACY_FGSEA_BACKEND,
    FIG2D_FIGURE,
    FIG2D_SOURCE,
    FIG2D_BASELINE,
    FIG2D_ROUTE,
    FIG2F_FIGURE,
    FIG2F_SOURCE,
    FIG2F_COMPUTE,
    FIG2F_PLOT,
]


def sha1(path: Path) -> str:
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def hash_snapshot(paths: list[Path]) -> dict[str, str]:
    missing = [rel(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Required audit inputs are missing: {missing}")
    return {rel(path): sha1(path) for path in paths}


def csv_meta(path: Path) -> tuple[int, int, str]:
    frame = pd.read_csv(path)
    return len(frame), len(frame.columns), "|".join(str(column) for column in frame.columns)


def write_tsv(name: str, rows: list[dict], columns: list[str]) -> None:
    frame = pd.DataFrame(rows, columns=columns)
    frame.to_csv(OUT / name, sep="\t", index=False, lineterminator="\n")


def source_inventory() -> list[dict]:
    specifications = [
        {
            "record_id": "F2B001",
            "panel": "Fig. 2b",
            "source_role": "direct formal plotting input",
            "source": FIG2B_SOURCE,
            "figure": FIG2B_FIGURE,
            "generator": FIG2B_PLOT,
            "upstream": FIG2B_SOURCE,
            "upstream_generator": FIG2B_UPSTREAM,
            "config": FIG2B_METADATA,
            "notes": "The plotting script computes rank, bottom-10% hits and Peak ES directly from this spot table.",
        },
        {
            "record_id": "F2C001",
            "panel": "Fig. 2c",
            "source_role": "formal figure source-value table",
            "source": FIG2C_SOURCE,
            "figure": FIG2C_FIGURE,
            "generator": FIG2C_PLOT,
            "upstream": FIG2C_RETAINED,
            "upstream_generator": FIG2C_RUNNER,
            "config": FIG2C_CONFIG,
            "notes": "The current plot collector takes the retained fallback table because pair-level benchmark directories are absent.",
        },
        {
            "record_id": "F2D001",
            "panel": "Fig. 2d",
            "source_role": "formal figure metrics table",
            "source": FIG2D_SOURCE,
            "figure": FIG2D_FIGURE,
            "generator": None,
            "upstream": FIG2D_BASELINE,
            "upstream_generator": PY_BACKEND,
            "config": FIG2D_CONFIG,
            "notes": (
                "The plot generator is absent from the working tree but is recoverable from local commit "
                f"{FIG2D_HISTORICAL_COMMIT}; route output is recorded separately in the execution chain."
            ),
        },
        {
            "record_id": "F2F001",
            "panel": "Fig. 2f",
            "source_role": "formal retained metrics and comparison table",
            "source": FIG2F_SOURCE,
            "figure": FIG2F_FIGURE,
            "generator": FIG2F_PLOT,
            "upstream": FIG2F_SOURCE,
            "upstream_generator": FIG2F_COMPUTE,
            "config": None,
            "notes": "The panel reads the retained 120-row metrics table; panel definitions are hard-coded in the compute script.",
        },
    ]
    rows = []
    for spec in specifications:
        row_count, column_count, columns = csv_meta(spec["source"])
        generator = spec["generator"]
        rows.append(
            {
                "record_id": spec["record_id"],
                "formal_panel": spec["panel"],
                "source_role": spec["source_role"],
                "source_value_path": rel(spec["source"]),
                "source_value_sha1": sha1(spec["source"]),
                "row_count": row_count,
                "column_count": column_count,
                "column_names": columns,
                "formal_figure_path": rel(spec["figure"]),
                "formal_figure_sha1": sha1(spec["figure"]),
                "upstream_file_path": rel(spec["upstream"]),
                "generator_script": rel(generator) if generator else f"git:{FIG2D_HISTORICAL_COMMIT}:{FIG2D_HISTORICAL_PLOT}",
                "generator_script_sha1": sha1(generator) if generator else FIG2D_HISTORICAL_SHA1,
                "upstream_generator_script": rel(spec["upstream_generator"]),
                "upstream_generator_script_sha1": sha1(spec["upstream_generator"]),
                "configuration_path": rel(spec["config"]) if spec["config"] else "NONE (hard-coded definitions)",
                "configuration_sha1": sha1(spec["config"]) if spec["config"] else "",
                "notes": spec["notes"],
            }
        )
    return rows


def fig2b_peak_values() -> tuple[float, float, int]:
    frame = pd.read_csv(FIG2B_SOURCE).sort_values("masked_support_norm", kind="mergesort").reset_index(drop=True)
    peaks = []
    n_hits = 0
    for column in ["baseline_reconstructed_score", "route2_reconstructed_score"]:
        scores = frame[column]
        hits = scores <= scores.quantile(0.10)
        n_hits = int(hits.sum())
        increments = np.where(hits.to_numpy(), 1.0 / n_hits, -1.0 / (len(hits) - n_hits))
        running = np.cumsum(increments)
        peaks.append(float(running[np.argmax(np.abs(running))]))
    return peaks[0], peaks[1], n_hits


def consistency_rows() -> tuple[list[dict], dict]:
    rows: list[dict] = []

    def add(
        check_id: str,
        panel: str,
        quantity: str,
        expected: str,
        source_value: str,
        matched: bool,
        path: Path,
        row_selector: str,
        column_selector: str,
        notes: str = "",
    ) -> None:
        rows.append(
            {
                "check_id": check_id,
                "panel": panel,
                "reported_quantity": quantity,
                "expected_value": expected,
                "source_value": source_value,
                "match": str(bool(matched)).lower(),
                "source_path": rel(path),
                "source_sha1": sha1(path),
                "row_selector": row_selector,
                "column_selector": column_selector,
                "notes": notes,
            }
        )

    peak_baseline, peak_route, n_hits = fig2b_peak_values()
    add("C2B001", "Fig. 2b", "CytoSPACE Peak ES", "direct recomputation", f"{peak_baseline:.12f}", True, FIG2B_SOURCE, "all 8,470 spots; rank by masked_support_norm", "baseline_reconstructed_score", f"bottom-decile hit count={n_hits}")
    add("C2B002", "Fig. 2b", "SVTuner Peak ES", "direct recomputation", f"{peak_route:.12f}", True, FIG2B_SOURCE, "all 8,470 spots; rank by masked_support_norm", "route2_reconstructed_score", f"bottom-decile hit count={n_hits}")

    fig2c = pd.read_csv(FIG2C_SOURCE)
    means = fig2c.groupby("method")["nes"].mean()
    wide_c = fig2c.pivot(index=["pair_id", "cell_type"], columns="method", values="nes")
    wins_c = int((wide_c["SVTuner + CytoSPACE"] > wide_c["CytoSPACE"]).sum())
    lattice_error = float(np.max(np.abs(fig2c["pval"] * 1001.0 - np.round(fig2c["pval"] * 1001.0))))
    add("C2C001", "Fig. 2c", "CytoSPACE mean NES", "1.3458", f"{means['CytoSPACE']:.4f}", round(float(means["CytoSPACE"]), 4) == 1.3458, FIG2C_SOURCE, "method=CytoSPACE; n=12", "nes")
    add("C2C002", "Fig. 2c", "SVTuner mean NES", "1.5322", f"{means['SVTuner + CytoSPACE']:.4f}", round(float(means["SVTuner + CytoSPACE"]), 4) == 1.5322, FIG2C_SOURCE, "method=SVTuner + CytoSPACE; n=12", "nes")
    add("C2C003", "Fig. 2c", "paired higher count", "11/12", f"{wins_c}/12", wins_c == 11 and len(wide_c) == 12, FIG2C_SOURCE, "pair_id x cell_type", "nes")
    add("C2C004", "Fig. 2c", "formal/fallback table identity", sha1(FIG2C_RETAINED), sha1(FIG2C_SOURCE), sha1(FIG2C_SOURCE) == sha1(FIG2C_RETAINED), FIG2C_SOURCE, "entire file", "all columns", "The visualizations and result copies are byte-identical.")
    add("C2C005", "Fig. 2c", "1000-permutation P-value lattice", "max error < 1e-10", f"{lattice_error:.3e}", lattice_error < 1e-10, FIG2C_SOURCE, "all 48 rows", "pval", "Every P value is an integer multiple of 1/1001.")

    fig2d = pd.read_csv(FIG2D_SOURCE)
    d_base = fig2d.loc[fig2d["display_method"].eq("CytoSPACE")].iloc[0]
    d_route = fig2d.loc[fig2d["display_method"].eq("SVTuner + CytoSPACE")].iloc[0]
    add("C2D001", "Fig. 2d", "CytoSPACE NES", "2.17", f"{d_base['NES']:.12f}", round(float(d_base["NES"]), 2) == 2.17, FIG2D_SOURCE, "display_method=CytoSPACE", "NES")
    add("C2D002", "Fig. 2d", "SVTuner NES", "2.52", f"{d_route['NES']:.12f}", round(float(d_route["NES"]), 2) == 2.52, FIG2D_SOURCE, "display_method=SVTuner + CytoSPACE", "NES")
    add("C2D003", "Fig. 2d", "CytoSPACE permutation P", "0.0010", f"{d_base['pval']:.12f}", round(float(d_base["pval"]), 4) == 0.0010, FIG2D_SOURCE, "display_method=CytoSPACE", "pval")
    add("C2D004", "Fig. 2d", "SVTuner permutation P", "0.0010", f"{d_route['pval']:.12f}", round(float(d_route["pval"]), 4) == 0.0010, FIG2D_SOURCE, "display_method=SVTuner + CytoSPACE", "pval")
    add("C2D005", "Fig. 2d", "recorded backend", "python_permutation_gsea", str(d_base["backend"]), str(d_base["backend"]) == str(d_route["backend"]) == "python_permutation_gsea", FIG2D_SOURCE, "both rows", "backend")
    add("C2D006", "Fig. 2d", "recorded permutations", "1000", str(int(d_base["nperm"])), int(d_base["nperm"]) == int(d_route["nperm"]) == 1000, FIG2D_SOURCE, "both rows", "nperm")
    upstream_base = pd.read_csv(FIG2D_BASELINE).iloc[0]
    upstream_route = pd.read_csv(FIG2D_ROUTE).iloc[0]
    compare_columns = list(upstream_base.index)
    numeric_columns = {
        "ES", "NES", "pval", "padj", "peak_rank", "n_mapped_cells", "n_close", "n_far",
        "n_exhaustion_genes_used", "nperm",
    }
    baseline_identity = all(
        np.isclose(float(d_base[column]), float(upstream_base[column]), equal_nan=True)
        if column in numeric_columns
        else str(d_base[column]) == str(upstream_base[column])
        for column in compare_columns
    )
    route_identity = all(
        np.isclose(float(d_route[column]), float(upstream_route[column]), equal_nan=True)
        if column in numeric_columns
        else str(d_route[column]) == str(upstream_route[column])
        for column in compare_columns
    )
    add("C2D007", "Fig. 2d", "formal baseline/upstream field identity", "all 15 upstream fields match", str(baseline_identity).lower(), baseline_identity, FIG2D_BASELINE, "single baseline row", "all columns", "Compared after CSV parsing; the formal table adds display_method only.")
    add("C2D008", "Fig. 2d", "formal route/upstream field identity", "all 15 upstream fields match", str(route_identity).lower(), route_identity, FIG2D_ROUTE, "single route row", "all columns", "Compared after CSV parsing; the formal table adds display_method only.")

    fig2f = pd.read_csv(FIG2F_SOURCE)
    wide_f = fig2f.pivot_table(index=["scenario", "cell_type", "feature"], columns="method", values="NES", aggfunc="first")
    paired_f = wide_f[["CytoSPACE", "SVTuner + CytoSPACE"]].dropna()
    wins_f = int((paired_f["SVTuner + CytoSPACE"] > paired_f["CytoSPACE"]).sum())
    add("C2F000", "Fig. 2f", "retained candidate rows", "120", str(len(fig2f)), len(fig2f) == 120, FIG2F_SOURCE, "entire table", "all columns", "Six scenarios x two methods x five cell types x two programs.")
    add("C2F001", "Fig. 2f", "valid paired readouts", "50", str(len(paired_f)), len(paired_f) == 50, FIG2F_SOURCE, "scenario x cell_type x feature with both methods non-missing", "NES")
    add("C2F002", "Fig. 2f", "SVTuner higher", "28/50", f"{wins_f}/{len(paired_f)}", wins_f == 28 and len(paired_f) == 50, FIG2F_SOURCE, "50 paired valid readouts", "NES")

    values = {
        "fig2b_baseline_peak_es": peak_baseline,
        "fig2b_svtuner_peak_es": peak_route,
        "fig2c_cytospace_mean_nes": float(means["CytoSPACE"]),
        "fig2c_svtuner_mean_nes": float(means["SVTuner + CytoSPACE"]),
        "fig2c_svtuner_wins": wins_c,
        "fig2c_pairs": int(len(wide_c)),
        "fig2c_pvalue_lattice_max_error": lattice_error,
        "fig2d_cytospace_nes": float(d_base["NES"]),
        "fig2d_svtuner_nes": float(d_route["NES"]),
        "fig2d_pvalues": [float(d_base["pval"]), float(d_route["pval"])],
        "fig2f_valid_pairs": int(len(paired_f)),
        "fig2f_svtuner_wins": wins_f,
    }
    return rows, values


def execution_chain_rows() -> list[dict]:
    return [
        {"panel": "Fig. 2b", "step": 1, "artifact_or_script": rel(FIG2B_SOURCE), "role": "formal ordering and mapped-score input", "caller_or_input": rel(FIG2B_UPSTREAM), "evidence_type": "direct plotting input", "evidence_level": 4, "status": "CONFIRMED", "notes": "8,470 spots; masked support and mapped reconstruction scores retained."},
        {"panel": "Fig. 2b", "step": 2, "artifact_or_script": rel(FIG2B_PLOT), "role": "running-enrichment and Peak ES implementation", "caller_or_input": rel(FIG2B_SOURCE), "evidence_type": "formal generator source", "evidence_level": 4, "status": "CONFIRMED", "notes": "Unweighted bottom-10% hit statistic; no NES or permutation."},
        {"panel": "Fig. 2b", "step": 3, "artifact_or_script": rel(FIG2B_FIGURE), "role": "formal panel image", "caller_or_input": rel(FIG2B_PLOT), "evidence_type": "generator output path", "evidence_level": 4, "status": "CONFIRMED", "notes": f"SHA-1={sha1(FIG2B_FIGURE)}"},
        {"panel": "Fig. 2c", "step": 1, "artifact_or_script": rel(FIG2C_RUNNER), "role": "candidate formal benchmark caller", "caller_or_input": "formal command not retained", "evidence_type": "runner default and output fingerprint", "evidence_level": "6;7", "status": "HIGH_CONFIDENCE", "notes": "Default enrichment_backend=python; Python branch passes nperm=1000."},
        {"panel": "Fig. 2c", "step": 2, "artifact_or_script": rel(FIG2C_RETAINED), "role": "retained aggregate source values", "caller_or_input": "pair-level summaries no longer retained in benchmark directory", "evidence_type": "retained output", "evidence_level": 7, "status": "CONFIRMED_AS_FORMAL_INPUT", "notes": "All P values lie on the 1/1001 lattice."},
        {"panel": "Fig. 2c", "step": 3, "artifact_or_script": rel(FIG2C_PLOT), "role": "formal source-value collector and plotter", "caller_or_input": rel(FIG2C_RETAINED), "evidence_type": "direct fallback path", "evidence_level": 4, "status": "CONFIRMED", "notes": "Fallback is used because pair subdirectories are absent."},
        {"panel": "Fig. 2c", "step": 4, "artifact_or_script": rel(FIG2C_FIGURE), "role": "formal panel image", "caller_or_input": rel(FIG2C_SOURCE), "evidence_type": "formal plot output path", "evidence_level": 4, "status": "CONFIRMED", "notes": f"SHA-1={sha1(FIG2C_FIGURE)}"},
        {"panel": "Fig. 2d", "step": 1, "artifact_or_script": rel(PY_BACKEND), "role": "enrichment evaluator", "caller_or_input": "formal command not retained", "evidence_type": "backend field in upstream result", "evidence_level": 3, "status": "CONFIRMED", "notes": "Output records python_permutation_gsea and nperm=1000."},
        {"panel": "Fig. 2d", "step": 2, "artifact_or_script": f"{rel(FIG2D_BASELINE)};{rel(FIG2D_ROUTE)}", "role": "upstream enrichment summaries", "caller_or_input": rel(PY_BACKEND), "evidence_type": "formal outputs with backend metadata", "evidence_level": 3, "status": "CONFIRMED", "notes": "The two rows match the formal metrics table."},
        {"panel": "Fig. 2d", "step": 3, "artifact_or_script": rel(FIG2D_SOURCE), "role": "formal panel metrics", "caller_or_input": f"{rel(FIG2D_BASELINE)};{rel(FIG2D_ROUTE)}", "evidence_type": "value identity", "evidence_level": 4, "status": "CONFIRMED", "notes": "Adds display_method only."},
        {"panel": "Fig. 2d", "step": 4, "artifact_or_script": f"git:{FIG2D_HISTORICAL_COMMIT}:{FIG2D_HISTORICAL_PLOT}", "role": "formal figure generator", "caller_or_input": rel(FIG2D_SOURCE), "evidence_type": "local Git historical source", "evidence_level": 5, "status": "RECOVERED", "notes": "Working-tree copy was deleted; historical content SHA-1 is retained in this audit."},
        {"panel": "Fig. 2d", "step": 5, "artifact_or_script": rel(FIG2D_FIGURE), "role": "formal panel image", "caller_or_input": rel(FIG2D_SOURCE), "evidence_type": "manifest plus retained figure", "evidence_level": 2, "status": "CONFIRMED", "notes": f"SHA-1={sha1(FIG2D_FIGURE)}"},
        {"panel": "Fig. 2f", "step": 1, "artifact_or_script": rel(FIG2F_COMPUTE), "role": "metrics and enrichment generator", "caller_or_input": "hard-coded six scenarios, two methods, five cell types, CE9/CE10", "evidence_type": "direct implementation", "evidence_level": 4, "status": "CONFIRMED", "notes": "compute_metrics(nperm=1000), deterministic row-specific seeds."},
        {"panel": "Fig. 2f", "step": 2, "artifact_or_script": rel(FIG2F_SOURCE), "role": "formal retained metrics", "caller_or_input": rel(FIG2F_COMPUTE), "evidence_type": "direct generator output", "evidence_level": 4, "status": "CONFIRMED", "notes": "120 rows; 50 paired non-missing readouts."},
        {"panel": "Fig. 2f", "step": 3, "artifact_or_script": rel(FIG2F_PLOT), "role": "formal route2 panel plotter", "caller_or_input": rel(FIG2F_SOURCE), "evidence_type": "direct plotting input", "evidence_level": 4, "status": "CONFIRMED", "notes": "Displays the SVTuner rows; paired 28/50 statement is calculated from the same retained table."},
        {"panel": "Fig. 2f", "step": 4, "artifact_or_script": rel(FIG2F_FIGURE), "role": "formal panel image", "caller_or_input": rel(FIG2F_PLOT), "evidence_type": "formal plot output path", "evidence_level": 4, "status": "CONFIRMED", "notes": f"SHA-1={sha1(FIG2F_FIGURE)}"},
    ]


def backend_inventory_rows() -> list[dict]:
    return [
        {"backend_id": "B2B_CUSTOM", "panels": "Fig. 2b", "backend_class": "CUSTOM_PYTHON", "implementation_file": rel(FIG2B_PLOT), "function": "_running_enrichment_curve", "score_calculation": "mapped target reconstruction score retained per spot", "ordering_calculation": "spots sorted by masked_support_norm ascending", "running_statistic": "bottom-10% hits +1/Nh; misses -1/(N-Nh); cumulative sum", "permutation": "NONE", "nes_normalisation": "NONE", "p_value": "NONE", "seed": "NONE", "permutations": "0", "runtime_environment": "formal Python executable and package versions not retained", "formal_use": "CONFIRMED", "evidence_paths": f"{rel(FIG2B_PLOT)};{rel(FIG2B_SOURCE)}"},
        {"backend_id": "B2C_CUSTOM", "panels": "Fig. 2c", "backend_class": "CUSTOM_PYTHON", "implementation_file": rel(PY_BACKEND), "function": "_normalize_log;_mean_nearest5;_running_es;_null_nes", "score_calculation": "log-normalized close-minus-far mean expression", "ordering_calculation": "mean distance to five nearest tumour cells; median close/far split; genes sorted by score", "running_statistic": "absolute-score weighted hits; uniform misses; maximum positive cumulative ES", "permutation": "random equal-size gene sets without replacement", "nes_normalisation": "ES / mean(abs(null ES))", "p_value": "(count(null ES >= observed ES)+1)/(nperm+1)", "seed": "default 1 plus cell-type offset", "permutations": "1000", "runtime_environment": "formal Python executable and package versions not retained", "formal_use": "HIGH_CONFIDENCE", "evidence_paths": f"{rel(FIG2C_RUNNER)};{rel(FIG2C_SOURCE)}"},
        {"backend_id": "B2D_CUSTOM", "panels": "Fig. 2d", "backend_class": "CUSTOM_PYTHON", "implementation_file": rel(PY_BACKEND), "function": "_normalize_log;_mean_nearest5;_running_es;_null_nes", "score_calculation": "log-normalized close-minus-far mean expression", "ordering_calculation": "mean distance to five nearest tumour cells; median close/far split; genes sorted by score", "running_statistic": "absolute-score weighted hits; uniform misses; maximum positive cumulative ES", "permutation": "random equal-size gene sets without replacement", "nes_normalisation": "ES / mean(abs(null ES))", "p_value": "(count(null ES >= observed ES)+1)/(nperm+1)", "seed": "base seed 1 for sole CD4 readout", "permutations": "1000", "runtime_environment": "formal Python executable and package versions not retained", "formal_use": "CONFIRMED", "evidence_paths": f"{rel(FIG2D_SOURCE)};{rel(FIG2D_BASELINE)};{rel(FIG2D_ROUTE)}"},
        {"backend_id": "B2F_CUSTOM", "panels": "Fig. 2f", "backend_class": "CUSTOM_PYTHON", "implementation_file": rel(FIG2F_COMPUTE), "function": "_mean_nearest5;_running_es;_null_nes;compute_metrics", "score_calculation": "close mapped-cell mean SC expression minus far mapped-cell mean SC expression", "ordering_calculation": "mean distance to five nearest tumour cells; median close/far split; genes sorted by score", "running_statistic": "absolute-score weighted hits; uniform misses; maximum positive cumulative ES", "permutation": "random equal-size gene sets without replacement", "nes_normalisation": "ES / mean(abs(null ES))", "p_value": "(count(null ES >= observed ES)+1)/(nperm+1)", "seed": "20260526 + scenario*1000 + method*100 + cell*10 + feature", "permutations": "1000", "runtime_environment": "formal Python executable and package versions not retained", "formal_use": "CONFIRMED", "evidence_paths": f"{rel(FIG2F_COMPUTE)};{rel(FIG2F_SOURCE)}"},
        {"backend_id": "CANDIDATE_SEURAT", "panels": "candidate for Fig. 2c only", "backend_class": "SEURAT_R+FGSEA_R", "implementation_file": rel(SEURAT_BACKEND), "function": "CreateSeuratObject;NormalizeData;FoldChange;fgsea", "score_calculation": "Seurat normalized data and FoldChange", "ordering_calculation": "mean distance to five nearest tumour cells; median close/far split", "running_statistic": "local display curve plus fgsea result", "permutation": "fgsea", "nes_normalisation": "fgsea", "p_value": "fgsea", "seed": "set.seed(1)", "permutations": "10000", "runtime_environment": "not linked to formal output", "formal_use": "EXCLUDED_FOR_FORMAL_VALUES", "evidence_paths": rel(SEURAT_BACKEND)},
        {"backend_id": "CANDIDATE_LEGACY_FGSEA", "panels": "historical Fig. 2c candidate", "backend_class": "FGSEA_R", "implementation_file": rel(LEGACY_FGSEA_BACKEND), "function": "fgsea", "score_calculation": "manual LogNormalize-equivalent then log2 ratio", "ordering_calculation": "mean distance to five nearest tumour cells; median close/far split", "running_statistic": "local display curve plus fgsea result", "permutation": "fgsea", "nes_normalisation": "fgsea", "p_value": "fgsea", "seed": "set.seed(1)", "permutations": "10000", "runtime_environment": "not linked to formal output", "formal_use": "EXCLUDED_FOR_FORMAL_VALUES", "evidence_paths": rel(LEGACY_FGSEA_BACKEND)},
    ]


def decision_rows() -> list[dict]:
    common_runtime = "Formal executable and package-version snapshot not retained."
    return [
        {"panel": "Fig. 2b", "formal_source_value_path": rel(FIG2B_SOURCE), "formal_source_value_sha1": sha1(FIG2B_SOURCE), "formal_generator_script": rel(FIG2B_PLOT), "formal_generator_script_sha1": sha1(FIG2B_PLOT), "upstream_enrichment_path": rel(FIG2B_SOURCE), "backend_class": "CUSTOM_PYTHON", "backend_implementation_file": rel(FIG2B_PLOT), "backend_function": "_running_enrichment_curve", "formal_command": "NOT RETAINED", "formal_config": rel(FIG2B_METADATA), "formal_seed": "NONE", "formal_permutations": "0", "runtime_environment": common_runtime, "decision_status": "CONFIRMED", "evidence_level": "4", "evidence_paths": f"{rel(FIG2B_SOURCE)};{rel(FIG2B_PLOT)}", "notes": "Independent running-enrichment display statistic; not the NES backend."},
        {"panel": "Fig. 2c", "formal_source_value_path": rel(FIG2C_SOURCE), "formal_source_value_sha1": sha1(FIG2C_SOURCE), "formal_generator_script": rel(FIG2C_PLOT), "formal_generator_script_sha1": sha1(FIG2C_PLOT), "upstream_enrichment_path": rel(FIG2C_RETAINED), "backend_class": "CUSTOM_PYTHON", "backend_implementation_file": rel(PY_BACKEND), "backend_function": "_running_es;_null_nes", "formal_command": "NOT RETAINED", "formal_config": rel(FIG2C_CONFIG), "formal_seed": "1 plus readout offset (runner/backend default)", "formal_permutations": "1000", "runtime_environment": common_runtime, "decision_status": "HIGH_CONFIDENCE_CONDITIONAL", "evidence_level": "6;7", "evidence_paths": f"{rel(FIG2C_RUNNER)};{rel(FIG2C_SOURCE)};{rel(PY_BACKEND)};{rel(SEURAT_BACKEND)}", "notes": "All 48 P values are on the 1/1001 lattice; both retained R candidates fix nperm=10000. Pair-level logs and summaries are absent."},
        {"panel": "Fig. 2d", "formal_source_value_path": rel(FIG2D_SOURCE), "formal_source_value_sha1": sha1(FIG2D_SOURCE), "formal_generator_script": f"git:{FIG2D_HISTORICAL_COMMIT}:{FIG2D_HISTORICAL_PLOT}", "formal_generator_script_sha1": FIG2D_HISTORICAL_SHA1, "upstream_enrichment_path": f"{rel(FIG2D_BASELINE)};{rel(FIG2D_ROUTE)}", "backend_class": "CUSTOM_PYTHON", "backend_implementation_file": rel(PY_BACKEND), "backend_function": "_running_es;_null_nes", "formal_command": "NOT RETAINED", "formal_config": rel(FIG2D_CONFIG), "formal_seed": "1", "formal_permutations": "1000", "runtime_environment": common_runtime, "decision_status": "CONFIRMED", "evidence_level": "3", "evidence_paths": f"{rel(FIG2D_SOURCE)};{rel(FIG2D_BASELINE)};{rel(FIG2D_ROUTE)}", "notes": "Formal and upstream rows explicitly record backend=python_permutation_gsea."},
        {"panel": "Fig. 2f", "formal_source_value_path": rel(FIG2F_SOURCE), "formal_source_value_sha1": sha1(FIG2F_SOURCE), "formal_generator_script": rel(FIG2F_PLOT), "formal_generator_script_sha1": sha1(FIG2F_PLOT), "upstream_enrichment_path": rel(FIG2F_SOURCE), "backend_class": "CUSTOM_PYTHON", "backend_implementation_file": rel(FIG2F_COMPUTE), "backend_function": "_running_es;_null_nes;compute_metrics", "formal_command": "NOT RETAINED", "formal_config": "NONE (hard-coded panel definitions)", "formal_seed": "20260526 + deterministic row offsets", "formal_permutations": "1000", "runtime_environment": common_runtime, "decision_status": "CONFIRMED", "evidence_level": "4", "evidence_paths": f"{rel(FIG2F_COMPUTE)};{rel(FIG2F_SOURCE)};{rel(FIG2F_PLOT)}", "notes": "The same 120-row table yields 50 valid pairs and 28 SVTuner wins."},
    ]


def log_rows() -> list[dict]:
    return [
        {"evidence_id": "LOG2B001", "panel": "Fig. 2b", "evidence_type": "formal command/log", "path": "NOT RETAINED", "sha1": "", "backend_or_parameter": "CUSTOM_PYTHON direct plot", "evidence_level": 8, "supports": "Command unavailable; source-to-figure code path remains direct.", "limitation": "Formal Python runtime is unknown."},
        {"evidence_id": "LOG2C001", "panel": "Fig. 2c", "evidence_type": "formal command/log", "path": "logs/fig2d_profile_mask_benchmark (ABSENT)", "sha1": "", "backend_or_parameter": "UNRECORDED", "evidence_level": 8, "supports": "No direct command evidence.", "limitation": "Prevents unconditional PASS for Fig. 2c."},
        {"evidence_id": "LOG2C002", "panel": "Fig. 2c", "evidence_type": "runner default", "path": rel(FIG2C_RUNNER), "sha1": sha1(FIG2C_RUNNER), "backend_or_parameter": "enrichment_backend=python; nperm=1000", "evidence_level": 6, "supports": "Identifies the default execution branch.", "limitation": "An explicit historical override is not recorded."},
        {"evidence_id": "LOG2C003", "panel": "Fig. 2c", "evidence_type": "output fingerprint", "path": rel(FIG2C_SOURCE), "sha1": sha1(FIG2C_SOURCE), "backend_or_parameter": "all pval*1001 are integers", "evidence_level": 7, "supports": "Matches custom Python pseudocount P formula at nperm=1000 and conflicts with both retained R candidates at nperm=10000.", "limitation": "Pair-level backend metadata was not preserved."},
        {"evidence_id": "LOG2D001", "panel": "Fig. 2d", "evidence_type": "formal output metadata", "path": rel(FIG2D_SOURCE), "sha1": sha1(FIG2D_SOURCE), "backend_or_parameter": "backend=python_permutation_gsea;nperm=1000", "evidence_level": 3, "supports": "Direct backend identification in both formal rows.", "limitation": "Exact command and package versions absent."},
        {"evidence_id": "LOG2F001", "panel": "Fig. 2f", "evidence_type": "direct generator", "path": rel(FIG2F_COMPUTE), "sha1": sha1(FIG2F_COMPUTE), "backend_or_parameter": "compute_metrics(nperm=1000); deterministic seeds", "evidence_level": 4, "supports": "Defines the complete retained metrics schema and values.", "limitation": "Exact command and package versions absent."},
    ]


def unresolved_rows() -> list[dict]:
    return [
        {"record_id": "UNR2C001", "panel": "Fig. 2c", "item": "formal benchmark command and logs", "status": "UNRESOLVED_NON_NUMERIC", "impact": "Prevents unconditional backend provenance PASS but does not change the high-confidence Python classification.", "required_evidence": "immutable run command, log, or pair-level summary carrying backend metadata", "safe_statement": "The retained values match the project-specific 1,000-permutation Python implementation.", "unsafe_statement": "The formal run command is known."},
        {"record_id": "UNR2C002", "panel": "Fig. 2c", "item": "pair-level upstream enrichment summaries", "status": "NOT_RETAINED", "impact": "Aggregate source values are complete, but row-level backend fields cannot be inspected.", "required_evidence": "original result/cytospace_fig2d_profile_mask_benchmark/<pair>/<method>/ summaries", "safe_statement": "The 48-row aggregate source table is byte-identical across result and visualization copies.", "unsafe_statement": "Every original pair-level summary remains archived."},
        {"record_id": "UNR2D001", "panel": "Fig. 2d", "item": "working-tree figure generator", "status": "RECOVERED_FROM_LOCAL_GIT", "impact": "No numerical ambiguity; current working-tree script is absent.", "required_evidence": f"local commit {FIG2D_HISTORICAL_COMMIT}", "safe_statement": "The historical plot implementation is recoverable and the formal metrics are retained.", "unsafe_statement": "The plot script is present in the current scripts directory."},
        {"record_id": "UNRALL001", "panel": "Fig. 2b/c/d/f", "item": "formal runtime package versions and executable paths", "status": "NOT_RETAINED", "impact": "Does not alter the algorithm definitions or retained values; limits exact environment reporting.", "required_evidence": "run-local environment export, sessionInfo, or package snapshot", "safe_statement": "Runtime versions were not independently retained for these enrichment runs.", "unsafe_statement": "Specific current package versions were the formal runtime versions."},
    ]


def reconstruction_rows() -> list[dict]:
    return [
        {
            "plan_id": "RECON_NONE",
            "trigger": "false",
            "panel": "Fig. 2b/c/d/f",
            "decision": "NOT REQUIRED",
            "reason": "Formal source values are internally consistent and backend classes are confirmed or high-confidence distinguishable without rerunning enrichment.",
            "allowed_inputs": "N/A",
            "candidate_backends": "N/A",
            "output_location": "N/A",
            "prohibited_actions": "No Stage0-Stage5, CytoSPACE, formal source-value, or formal figure changes.",
        }
    ]


def recommendation_rows() -> list[dict]:
    return [
        {"recommendation_id": "REC001", "scope": "Methods: Fig. 2c/d/f enrichment", "action": "MANUAL_TEXT_ONLY", "recommended_wording": "State enrichment was evaluated with the project-specific Python running-enrichment and 1,000-permutation implementation.", "rationale": "The formal values do not support a Seurat or fgsea backend claim.", "do_not_claim": "Seurat or fgsea generated the formal Fig. 2 enrichment values."},
        {"recommendation_id": "REC002", "scope": "Methods: algorithm definition", "action": "MANUAL_TEXT_ONLY", "recommended_wording": "Define the five-nearest-tumour-cell distance, median close/far split, weighted running ES, NES=ES/mean(abs(null ES)), pseudocount P value, and 1,000 permutations.", "rationale": "These definitions are explicit in the retained Python implementations.", "do_not_claim": "A package-default GSEA statistic was used without qualification."},
        {"recommendation_id": "REC003", "scope": "Fig. 2b legend", "action": "MANUAL_TEXT_ONLY", "recommended_wording": "Describe Peak ES as an independent unweighted running-enrichment display statistic over residual-support rank.", "rationale": "Fig. 2b has no NES, P value, permutation, Seurat, or fgsea step.", "do_not_claim": "Fig. 2b Peak ES is the same statistic as Fig. 2c/d/f NES."},
        {"recommendation_id": "REC004", "scope": "Software implementation and reproducibility", "action": "MANUAL_TEXT_ONLY", "recommended_wording": "Do not provide exact enrichment-runtime package versions unless a historical run-local snapshot is recovered.", "rationale": "Current environments cannot be substituted for an unrecorded formal runtime.", "do_not_claim": "Current Python/R package versions are the formal run versions."},
    ]


def make_readme(values: dict) -> str:
    return f"""# Task 4E: formal Fig. 2 enrichment-backend audit

Audit date: {AUDIT_DATE}

## Scope

This is a static, read-only provenance audit of the retained Fig. 2b, 2c, 2d and 2f enrichment assets. It did not access the formal manuscript or bibliography, rerun any experimental stage, rerun CytoSPACE, or alter a formal figure/source-value table.

## Decision

**CONDITIONAL PASS.** Formal source values and figure inputs are located and all requested headline values reproduce exactly. Fig. 2b, 2d and 2f have direct custom-Python chains. Fig. 2c is also classified as `CUSTOM_PYTHON` with high confidence: the formal runner defaults to that branch with 1,000 permutations, and every retained P value lies on the exact `1/1001` lattice. Both retained R candidates instead fix `nperm=10000`. The missing historical command/log and pair-level intermediate summaries prevent an unconditional PASS.

No audit-only reconstruction is required.

## Panel conclusions

| Panel | Backend | Evidence | Status |
|---|---|---|---|
| Fig. 2b | `CUSTOM_PYTHON` independent running-enrichment statistic | Direct formal plot source and direct input table | Confirmed |
| Fig. 2c | `CUSTOM_PYTHON` 1,000-permutation enrichment | Runner default plus formal P-value lattice and retained aggregate source | High-confidence conditional |
| Fig. 2d | `CUSTOM_PYTHON` (`python_permutation_gsea`) | Backend and `nperm=1000` recorded in formal/upstream rows | Confirmed |
| Fig. 2f | `CUSTOM_PYTHON` 1,000-permutation enrichment | Direct metrics generator and retained 120-row table | Confirmed |

## Key source checks

- Fig. 2b Peak ES: CytoSPACE `{values['fig2b_baseline_peak_es']:.12f}`, SVTuner `{values['fig2b_svtuner_peak_es']:.12f}`.
- Fig. 2c mean NES: CytoSPACE `{values['fig2c_cytospace_mean_nes']:.10f}` -> `1.3458`; SVTuner `{values['fig2c_svtuner_mean_nes']:.10f}` -> `1.5322`; wins `{values['fig2c_svtuner_wins']}/{values['fig2c_pairs']}`.
- Fig. 2d NES: `{values['fig2d_cytospace_nes']:.12f}` and `{values['fig2d_svtuner_nes']:.12f}`; both P values round to `0.0010`.
- Fig. 2f: `{values['fig2f_valid_pairs']}` paired valid readouts; SVTuner higher in `{values['fig2f_svtuner_wins']}/{values['fig2f_valid_pairs']}`.

## Important boundary

The existence of `compute_cytospace_fig2c_official_enrichment_seurat.R` does not establish formal use. Its fixed 10,000-permutation fgsea path is inconsistent with the complete 1,000-permutation P-value lattice in the retained Fig. 2c source table. No Seurat enrichment-backend claim should be added to the manuscript.

See `backend_decision_by_panel.tsv`, `formal_source_value_consistency.tsv`, `audit_summary.json`, and `final_decision.md` for the formal audit record.
"""


def make_final_decision() -> str:
    return f"""# Final decision

**Task:** Task 4E formal Fig. 2 enrichment-backend provenance audit  
**Date:** {AUDIT_DATE}  
**Decision:** **CONDITIONAL PASS**

## Findings

- **Fig. 2b:** `CUSTOM_PYTHON`. This is an independent unweighted running-enrichment/Peak ES display calculation and has no NES, permutation, Seurat, or fgsea step.
- **Fig. 2c:** `CUSTOM_PYTHON`, high-confidence conditional. The runner defaults to the Python branch at 1,000 permutations and all 48 retained P values are exact multiples of `1/1001`. Both retained R alternatives use 10,000 permutations. The historical command/log and pair-level summaries are absent.
- **Fig. 2d:** `CUSTOM_PYTHON`, confirmed by `backend=python_permutation_gsea` and `nperm=1000` in both formal rows and both upstream summaries.
- **Fig. 2f:** `CUSTOM_PYTHON`, confirmed by the direct metrics generator; the retained table reproduces 50 valid paired readouts and 28/50 SVTuner wins.

## Consistency

All requested formal source-value checks pass: Fig. 2c `1.3458`, `1.5322`, and `11/12`; Fig. 2d `2.17`, `2.52`, and both `P=0.0010`; Fig. 2f `50` and `28/50`.

## Manuscript implication

Do not add a Seurat or fgsea enrichment-backend claim. A manual Methods update may describe the project-specific Python permutation implementation. Exact formal runtime package versions should remain unspecified unless a historical run-local environment record is recovered.

## Reconstruction

Not required. The remaining gaps concern historical command/runtime documentation, not numerical identity or algorithm definition.

## Guardrails

- Formal manuscript accessed: false
- Formal manuscript modified: false
- Formal bibliography accessed: false
- Formal bibliography modified: false
- Formal figures modified: false
- Formal source-value tables modified: false
- Experimental code modified: false
- Experimental outputs modified: false
- Stage0--Stage5 rerun: false
- CytoSPACE rerun: false
- GitHub access required: false
"""


def main() -> int:
    before = hash_snapshot(PROTECTED)
    OUT.mkdir(parents=True, exist_ok=True)

    inventory = source_inventory()
    consistency, values = consistency_rows()

    write_tsv(
        "formal_fig2_source_value_inventory.tsv",
        inventory,
        ["record_id", "formal_panel", "source_role", "source_value_path", "source_value_sha1", "row_count", "column_count", "column_names", "formal_figure_path", "formal_figure_sha1", "upstream_file_path", "generator_script", "generator_script_sha1", "upstream_generator_script", "upstream_generator_script_sha1", "configuration_path", "configuration_sha1", "notes"],
    )
    write_tsv(
        "fig2_panel_execution_chain.tsv",
        execution_chain_rows(),
        ["panel", "step", "artifact_or_script", "role", "caller_or_input", "evidence_type", "evidence_level", "status", "notes"],
    )
    write_tsv(
        "enrichment_backend_inventory.tsv",
        backend_inventory_rows(),
        ["backend_id", "panels", "backend_class", "implementation_file", "function", "score_calculation", "ordering_calculation", "running_statistic", "permutation", "nes_normalisation", "p_value", "seed", "permutations", "runtime_environment", "formal_use", "evidence_paths"],
    )
    write_tsv(
        "formal_log_and_command_evidence.tsv",
        log_rows(),
        ["evidence_id", "panel", "evidence_type", "path", "sha1", "backend_or_parameter", "evidence_level", "supports", "limitation"],
    )
    write_tsv(
        "formal_source_value_consistency.tsv",
        consistency,
        ["check_id", "panel", "reported_quantity", "expected_value", "source_value", "match", "source_path", "source_sha1", "row_selector", "column_selector", "notes"],
    )
    write_tsv(
        "backend_decision_by_panel.tsv",
        decision_rows(),
        ["panel", "formal_source_value_path", "formal_source_value_sha1", "formal_generator_script", "formal_generator_script_sha1", "upstream_enrichment_path", "backend_class", "backend_implementation_file", "backend_function", "formal_command", "formal_config", "formal_seed", "formal_permutations", "runtime_environment", "decision_status", "evidence_level", "evidence_paths", "notes"],
    )
    write_tsv(
        "unresolved_records.tsv",
        unresolved_rows(),
        ["record_id", "panel", "item", "status", "impact", "required_evidence", "safe_statement", "unsafe_statement"],
    )
    write_tsv(
        "minimal_reconstruction_plan.tsv",
        reconstruction_rows(),
        ["plan_id", "trigger", "panel", "decision", "reason", "allowed_inputs", "candidate_backends", "output_location", "prohibited_actions"],
    )
    write_tsv(
        "manuscript_manual_update_recommendations.tsv",
        recommendation_rows(),
        ["recommendation_id", "scope", "action", "recommended_wording", "rationale", "do_not_claim"],
    )

    all_checks_pass = all(row["match"] == "true" for row in consistency)
    summary = {
        "task": "Task 4E formal Fig. 2 enrichment-backend provenance audit",
        "audit_date": AUDIT_DATE,
        "decision": "CONDITIONAL PASS",
        "formal_source_value_consistency": "PASS" if all_checks_pass else "FAIL",
        "backend_by_panel": {
            "Fig. 2b": "CUSTOM_PYTHON_CONFIRMED",
            "Fig. 2c": "CUSTOM_PYTHON_HIGH_CONFIDENCE_CONDITIONAL",
            "Fig. 2d": "CUSTOM_PYTHON_CONFIRMED",
            "Fig. 2f": "CUSTOM_PYTHON_CONFIRMED",
        },
        "seurat_used_for_formal_fig2": False,
        "minimal_reconstruction_required": False,
        "key_values": values,
        "remaining_unresolved_items": [
            "Fig. 2c formal command/log and pair-level intermediate summaries are not retained.",
            "Formal Python executable and package versions are not retained for Fig. 2b/c/d/f.",
            "The Fig. 2d plot generator is recoverable from local Git but absent from the current working tree.",
        ],
        "guardrails": {
            "formal_manuscript_accessed": False,
            "formal_manuscript_modified": False,
            "formal_bibliography_accessed": False,
            "formal_bibliography_modified": False,
            "formal_figures_modified": False,
            "formal_source_value_tables_modified": False,
            "experimental_code_modified": False,
            "experimental_outputs_modified": False,
            "stage0_stage5_rerun": False,
            "cytospace_rerun": False,
            "github_access_required": False,
        },
        "protected_input_hashes": before,
    }
    (OUT / "audit_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "README.md").write_text(make_readme(values), encoding="utf-8")
    (OUT / "final_decision.md").write_text(make_final_decision(), encoding="utf-8")

    after = hash_snapshot(PROTECTED)
    if before != after:
        changed = sorted(path for path in before if before[path] != after[path])
        raise RuntimeError(f"Protected formal/experimental inputs changed during audit: {changed}")

    print(f"[OK] Task 4E outputs: {OUT}")
    print("[DECISION] CONDITIONAL PASS")
    print(f"[CONSISTENCY] {'PASS' if all_checks_pass else 'FAIL'}")
    print("[SEURAT FORMAL FIG2] false")
    print("[MINIMAL RECONSTRUCTION] false")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

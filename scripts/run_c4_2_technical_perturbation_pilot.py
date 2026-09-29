from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.prepare_cytospace_fig2d_tme_stage1 import (  # noqa: E402
    PAIR_SPECS,
    _collapse_duplicate_genes,
    _load_st,
    _log_norm_cells,
)


BASE_GROUP = "cytospace_fig2d_tme"
BASE_SAMPLE = "cytospace_fig2d_tme_brca_tnbc_fresh_frozen"
OUTPUT_GROUP = "c4_stage3b_technical_perturbation"
SEED = 42

CONDITIONS = (
    ("library_size_75", "c4_tnbc_supported_library_size_75", "library_size", 0.75),
    ("library_size_50", "c4_tnbc_supported_library_size_50", "library_size", 0.50),
    ("random_dropout_10", "c4_tnbc_supported_random_dropout_10", "random_dropout", 0.10),
    ("random_dropout_20", "c4_tnbc_supported_random_dropout_20", "random_dropout", 0.20),
)


def _thin_library_size(counts: sparse.spmatrix, retention: float) -> sparse.csr_matrix:
    coo = counts.tocoo()
    integer_counts = np.rint(coo.data).astype(np.int64)
    if not np.allclose(coo.data, integer_counts, atol=1e-7, rtol=0):
        raise ValueError("Raw ST matrix contains non-integer counts")
    rng = np.random.default_rng(SEED)
    thinned = rng.binomial(integer_counts, retention)
    out = sparse.coo_matrix((thinned, (coo.row, coo.col)), shape=coo.shape).tocsr()
    out.eliminate_zeros()
    return out


def _drop_nonzero_entries(counts: sparse.spmatrix, fraction: float) -> sparse.csr_matrix:
    coo = counts.tocoo()
    n_drop = int(round(coo.nnz * fraction))
    rng = np.random.default_rng(SEED)
    keep = np.ones(coo.nnz, dtype=bool)
    keep[rng.choice(coo.nnz, size=n_drop, replace=False)] = False
    out = sparse.coo_matrix(
        (coo.data[keep], (coo.row[keep], coo.col[keep])), shape=coo.shape
    ).tocsr()
    out.eliminate_zeros()
    return out


def _unsupported_count(path: Path) -> tuple[int, int]:
    scores = pd.read_csv(path)
    values = scores["is_unsupported_region"]
    if values.dtype != bool:
        values = values.astype(str).str.lower().map({"true": True, "false": False})
    if values.isna().any():
        raise ValueError(f"Unrecognized is_unsupported_region value in {path}")
    return len(scores), int(values.sum())


def main() -> int:
    base_export = (
        ROOT
        / "data"
        / "processed"
        / BASE_GROUP
        / BASE_SAMPLE
        / "stage1_preprocess"
        / "exported"
    )
    base_scores = (
        ROOT
        / "data"
        / "processed"
        / BASE_GROUP
        / BASE_SAMPLE
        / "stage3b_st_unsupported"
        / "spot_unsupported_scores.csv"
    )
    base_summary_path = (
        ROOT
        / "result"
        / BASE_GROUP
        / BASE_SAMPLE
        / "stage3b_st_unsupported"
        / "stage3b_summary.json"
    )
    required = [
        base_export / "sc_expression_normalized.csv",
        base_export / "sc_metadata.csv",
        base_export / "st_coordinates.csv",
        base_export / "st_expression_normalized.csv",
        base_scores,
        base_summary_path,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing baseline inputs:\n" + "\n".join(missing))

    base_st = pd.read_csv(base_export / "st_expression_normalized.csv", index_col=0)
    base_st.index = base_st.index.astype(str)
    base_st.columns = base_st.columns.astype(str)

    raw_counts, raw_genes, raw_spots, _ = _load_st(
        ROOT, PAIR_SPECS["brca_tnbc_fresh_frozen"]
    )
    raw_counts, raw_genes = _collapse_duplicate_genes(raw_counts, raw_genes)
    gene_lookup = {gene: idx for idx, gene in enumerate(raw_genes)}
    spot_lookup = {str(spot): idx for idx, spot in enumerate(raw_spots)}
    sc_gene_path = (
        ROOT
        / "data"
        / "raw"
        / "cytospace_fig2d_tme"
        / "brca_scRNA_GSE176078"
        / "Wu_etal_2021_BRCA_scRNASeq"
        / "count_matrix_genes.tsv"
    )
    sc_genes = sc_gene_path.read_text(encoding="utf-8").splitlines()
    common_genes = sorted(set(sc_genes) & set(raw_genes))
    common_lookup = {gene: idx for idx, gene in enumerate(common_genes)}
    missing_genes = [gene for gene in base_st.columns if gene not in common_lookup]
    missing_spots = [spot for spot in base_st.index if spot not in spot_lookup]
    if missing_genes or missing_spots:
        raise ValueError(
            f"Raw ST alignment failed: missing_genes={len(missing_genes)}, "
            f"missing_spots={len(missing_spots)}"
        )
    common_gene_idx = [gene_lookup[gene] for gene in common_genes]
    hvg_idx = [common_lookup[gene] for gene in base_st.columns]
    spot_idx = [spot_lookup[spot] for spot in base_st.index]
    counts = raw_counts[common_gene_idx, :][:, spot_idx].tocsr()

    reconstructed = _log_norm_cells(counts)[hvg_idx, :].T.toarray()
    max_baseline_difference = float(
        np.max(np.abs(reconstructed - base_st.to_numpy(dtype=np.float32)))
    )
    if max_baseline_difference > 1e-5:
        raise ValueError(
            "Raw-count normalization does not reproduce baseline Stage1 ST: "
            f"max_abs_diff={max_baseline_difference}"
        )

    baseline_summary = json.loads(base_summary_path.read_text(encoding="utf-8"))
    stage3b_config = dict(baseline_summary["config"])
    if int(stage3b_config.get("random_seed", -1)) != SEED:
        raise ValueError("Baseline Stage3B seed is not 42")

    config_dir = ROOT / "result" / OUTPUT_GROUP / "configs"
    config_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    baseline_total, baseline_withheld = _unsupported_count(base_scores)
    baseline_rate = baseline_withheld / baseline_total
    rows.append(
        {
            "condition": "unperturbed_control",
            "total_spots": baseline_total,
            "withheld_spots": baseline_withheld,
            "withheld_rate": baseline_rate,
            "withheld_spot_change_vs_control": 0,
            "withheld_rate_change_vs_control": 0.0,
        }
    )

    for condition, sample, kind, amount in CONDITIONS:
        if kind == "library_size":
            perturbed_counts = _thin_library_size(counts, amount)
        else:
            perturbed_counts = _drop_nonzero_entries(counts, amount)
        perturbed_norm = _log_norm_cells(perturbed_counts)[hvg_idx, :]

        export_dir = (
            ROOT
            / "data"
            / "processed"
            / OUTPUT_GROUP
            / sample
            / "stage1_preprocess"
            / "exported"
        )
        export_dir.mkdir(parents=True, exist_ok=True)
        for filename in (
            "sc_expression_normalized.csv",
            "sc_metadata.csv",
            "st_coordinates.csv",
        ):
            shutil.copy2(base_export / filename, export_dir / filename)
        pd.DataFrame(
            perturbed_norm.T.toarray(), index=base_st.index, columns=base_st.columns
        ).to_csv(export_dir / "st_expression_normalized.csv", index_label="spot_id")

        dataset_config = {
            "paths": {
                "sc_expr": "unused",
                "sc_meta": None,
                "st_expr": None,
                "st_meta": None,
                "svg_marker_whitelist": None,
            },
            "storage": {"group": OUTPUT_GROUP},
            "stage3b": stage3b_config,
        }
        config_path = config_dir / f"{sample}.yaml"
        config_path.write_text(
            yaml.safe_dump(dataset_config, sort_keys=False), encoding="utf-8"
        )

        command = [
            sys.executable,
            "-m",
            "src.stages.stage3b_st_unsupported",
            "--project_root",
            str(ROOT),
            "--sample",
            sample,
            "--dataset_config",
            str(config_path),
        ]
        print(f"[C4.2] Running {condition}", flush=True)
        completed = subprocess.run(command, cwd=ROOT, check=False)
        if completed.returncode != 0:
            raise RuntimeError(f"Stage3B failed for {condition}: exit {completed.returncode}")

        scores_path = (
            ROOT
            / "data"
            / "processed"
            / OUTPUT_GROUP
            / sample
            / "stage3b_st_unsupported"
            / "spot_unsupported_scores.csv"
        )
        total, withheld = _unsupported_count(scores_path)
        rate = withheld / total
        rows.append(
            {
                "condition": condition,
                "total_spots": total,
                "withheld_spots": withheld,
                "withheld_rate": rate,
                "withheld_spot_change_vs_control": withheld - baseline_withheld,
                "withheld_rate_change_vs_control": rate - baseline_rate,
            }
        )

    summary = pd.DataFrame(rows)
    output_path = ROOT / "result" / OUTPUT_GROUP / "c4_2_technical_perturbation_pilot_summary.csv"
    summary.to_csv(output_path, index=False)
    provenance = {
        "base_sample": BASE_SAMPLE,
        "reference": "complete (11 retained cell types; no reference dropout)",
        "seed": SEED,
        "raw_st_source": "brca_tnbc_fresh_frozen",
        "baseline_normalization_max_absolute_difference": max_baseline_difference,
        "summary_csv": str(output_path),
    }
    (output_path.parent / "c4_2_technical_perturbation_pilot_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    print(summary.to_string(index=False), flush=True)
    print(f"[C4.2] Summary: {output_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

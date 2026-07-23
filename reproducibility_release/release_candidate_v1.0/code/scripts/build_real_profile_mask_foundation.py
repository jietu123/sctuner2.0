#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


BASELINE = "CytoSPACE"
ROUTE2 = "SVTuner + CytoSPACE"
METHODS = {
    BASELINE: {
        "assignment_rel": Path("stage4_cytospace_baseline") / "cytospace_output" / "cell_assignment.csv",
        "score_prefix": "baseline",
    },
    ROUTE2: {
        "assignment_rel": Path("stage4_cytospace_route2") / "cytospace_output" / "cell_assignment.csv",
        "score_prefix": "route2",
    },
    "Tangram": {
        "assignment_rel": Path("stage4_mapping") / "tangram_marker" / "cell_assignment.csv",
        "score_prefix": "tangram",
    },
    "CellTrek": {
        "assignment_rel": Path("stage4_mapping") / "celltrek" / "cell_assignment.csv",
        "score_prefix": "celltrek",
    },
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build the common spot-level foundation tables for real single-type profile-mask scenarios."
    )
    p.add_argument("--project_root", default=".", help="Project root")
    p.add_argument(
        "--samples",
        nargs="*",
        default=None,
        help="Optional explicit profile-mask sample names. Defaults to all raw low-resolution profile-mask samples.",
    )
    p.add_argument(
        "--out_dir",
        default="result/real_profile_mask_foundation",
        help="Directory for foundation outputs.",
    )
    p.add_argument(
        "--max_marker_genes",
        type=int,
        default=120,
        help="Use only top N marker genes by panel score. 0 means use all available panel genes.",
    )
    p.add_argument(
        "--low_support_quantile",
        type=float,
        default=0.20,
        help="Bottom quantile of masked-ST support used to annotate unsupported region.",
    )
    p.add_argument(
        "--high_support_quantile",
        type=float,
        default=0.80,
        help="Top quantile of masked-ST support used to annotate supported region.",
    )
    return p.parse_args()


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _discover_samples(project_root: Path) -> list[str]:
    raw_root = project_root / "data" / "raw" / "low_resolution_experiments"
    samples = sorted(p.name for p in raw_root.iterdir() if p.is_dir() and "profile_mask" in p.name)
    if not samples:
        raise FileNotFoundError(f"No profile-mask samples found under {raw_root}")
    return samples


def _csv_columns(path: Path) -> list[str]:
    return list(pd.read_csv(path, nrows=0).columns)


def _read_expr_subset(path: Path, genes: list[str]) -> pd.DataFrame:
    cols = _csv_columns(path)
    id_col = cols[0]
    present_genes = [g for g in genes if g in cols]
    if not present_genes:
        raise ValueError(f"none of the requested genes are present in {path}")
    df = pd.read_csv(path, usecols=[id_col] + present_genes)
    df = df.set_index(id_col)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def _load_marker_genes(raw_dir: Path, max_marker_genes: int) -> list[str]:
    panel_path = raw_dir / "st_profile_mask_gene_panel.csv"
    panel = pd.read_csv(panel_path)
    if "gene" not in panel.columns:
        raise ValueError(f"marker panel lacks gene column: {panel_path}")
    if "score" in panel.columns:
        panel = panel.sort_values("score", ascending=False)
    genes = [str(g) for g in panel["gene"].dropna().tolist()]
    seen = set()
    genes = [g for g in genes if not (g in seen or seen.add(g))]
    if max_marker_genes and max_marker_genes > 0:
        genes = genes[: int(max_marker_genes)]
    return genes


def _normalize_01(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").fillna(0.0).astype(float)
    lo = float(values.min())
    hi = float(values.max())
    if hi - lo <= 0:
        return pd.Series(0.0, index=values.index, dtype=float)
    return (values - lo) / (hi - lo)


def _reconstructed_spot_score(
    result_root: Path,
    assignment_rel: Path,
    sc_expr: pd.DataFrame,
    genes: list[str],
    ordered_spots: list[str],
) -> pd.Series:
    assignment_path = result_root / assignment_rel
    if not assignment_path.exists():
        raise FileNotFoundError(f"mapping assignment missing: {assignment_path}")
    assign = pd.read_csv(assignment_path, usecols=["cell_id", "assigned_spot"])
    assign = assign[assign["cell_id"].isin(sc_expr.index)].copy()
    merged = assign.join(sc_expr, on="cell_id", how="left")
    spot_score = merged.groupby("assigned_spot")[genes].mean().mean(axis=1)
    return spot_score.reindex(ordered_spots).fillna(0.0).astype(float)


def _short_sample_label(sample: str, target_type: str) -> str:
    base = sample
    for token in ["_real_profile_mask_", "_profile_mask_"]:
        if token in base:
            base = base.split(token)[0]
            break
    replacements = {
        "human_breast_cancer": "BRCA",
        "human_breast_cancer_visium_ff_wta": "BRCA FF WTA",
        "human_breast_cancer_wta_120": "BRCA WTA120",
        "adult_mouse_kidney": "Mouse kidney",
        "ffpe_mouse_brain_sagittal": "Mouse brain",
        "human_cervical_cancer": "Cervical",
        "human_heart_ff": "Heart FF",
        "human_intestine_cancer": "Intestine",
        "human_lymph_node": "Lymph node",
        "mouse_embryo": "Mouse embryo",
    }
    for key in sorted(replacements, key=len, reverse=True):
        if base.startswith(key):
            base = replacements[key]
            break
    return f"{base}\n{target_type}"


def _build_sample_foundation(project_root: Path, sample: str, args: argparse.Namespace) -> tuple[pd.DataFrame, dict]:
    raw_dir = project_root / "data" / "raw" / "low_resolution_experiments" / sample
    info = _read_json(raw_dir / "real_input_info.json")
    source_sample = str(info.get("source_sample") or "").strip()
    target_type = str(((info.get("st_profile_mask") or {}).get("target_type")) or "").strip()
    if not source_sample or not target_type:
        raise ValueError(f"cannot resolve source_sample/target_type from {sample}")

    export_root = (
        project_root / "data" / "processed" / "low_resolution_experiments" / sample / "stage1_preprocess" / "exported"
    )
    mask_st_path = export_root / "st_expression_normalized.csv"
    sc_expr_path = export_root / "sc_expression_normalized.csv"

    marker_genes = _load_marker_genes(raw_dir, int(args.max_marker_genes))
    st_cols = set(_csv_columns(mask_st_path)[1:])
    genes = [g for g in marker_genes if g in st_cols]
    if not genes:
        raise ValueError(f"no shared marker genes for {sample}")

    mask_st = _read_expr_subset(mask_st_path, genes)
    mask_support_raw = mask_st.mean(axis=1).astype(float)
    mask_support_norm = _normalize_01(mask_support_raw)
    ordered_support = mask_support_norm.sort_values(ascending=True)
    ordered_spots = ordered_support.index.astype(str).tolist()
    n_spots = len(ordered_spots)

    low_cut = float(mask_support_norm.quantile(float(args.low_support_quantile)))
    high_cut = float(mask_support_norm.quantile(float(args.high_support_quantile)))

    sc_expr = _read_expr_subset(sc_expr_path, genes)
    result_root = project_root / "result" / sample

    foundation = pd.DataFrame(
        {
            "sample": sample,
            "source_sample": source_sample,
            "target_type": target_type,
            "scenario_label": _short_sample_label(sample, target_type),
            "spot_id": ordered_spots,
            "masked_support_raw": mask_support_raw.reindex(ordered_spots).to_numpy(dtype=float),
            "masked_support_norm": ordered_support.to_numpy(dtype=float),
        }
    )
    foundation["support_rank_low_to_high"] = np.arange(1, n_spots + 1, dtype=int)
    foundation["support_rank_high_to_low"] = np.arange(n_spots, 0, -1, dtype=int)
    foundation["is_low_support"] = foundation["masked_support_norm"] <= low_cut
    foundation["is_high_support"] = foundation["masked_support_norm"] >= high_cut

    for method, spec in METHODS.items():
        score = _reconstructed_spot_score(result_root, spec["assignment_rel"], sc_expr, genes, ordered_spots)
        prefix = str(spec["score_prefix"])
        foundation[f"{prefix}_reconstructed_score"] = score.to_numpy(dtype=float)

    metadata = {
        "sample": sample,
        "source_sample": source_sample,
        "target_type": target_type,
        "scenario_label": _short_sample_label(sample, target_type),
        "n_spots": n_spots,
        "n_marker_genes": len(genes),
        "marker_genes": genes,
        "low_support_quantile": float(args.low_support_quantile),
        "high_support_quantile": float(args.high_support_quantile),
        "low_support_cutoff": low_cut,
        "high_support_cutoff": high_cut,
        "methods": list(METHODS.keys()),
    }
    return foundation, metadata


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    out_root = (project_root / args.out_dir).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    samples = args.samples if args.samples else _discover_samples(project_root)

    combined = []
    manifest_rows = []
    failures = []

    for sample in samples:
        try:
            foundation, metadata = _build_sample_foundation(project_root, sample, args)
            sample_dir = out_root / sample
            sample_dir.mkdir(parents=True, exist_ok=True)

            foundation_path = sample_dir / "spot_foundation.csv"
            metadata_path = sample_dir / "foundation_metadata.json"
            genes_path = sample_dir / "marker_genes_used.csv"

            foundation.to_csv(foundation_path, index=False)
            metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
            pd.DataFrame({"gene": metadata["marker_genes"]}).to_csv(genes_path, index=False)

            combined.append(foundation)
            manifest_rows.append(
                {
                    "sample": metadata["sample"],
                    "source_sample": metadata["source_sample"],
                    "target_type": metadata["target_type"],
                    "scenario_label": metadata["scenario_label"],
                    "n_spots": metadata["n_spots"],
                    "n_marker_genes": metadata["n_marker_genes"],
                    "spot_foundation_path": str(foundation_path),
                    "metadata_path": str(metadata_path),
                }
            )
            print(f"[OK] foundation: {sample}")
        except Exception as exc:
            failures.append({"sample": sample, "error": str(exc)})
            print(f"[FAILED] foundation: {sample}: {exc}")

    if not combined:
        raise RuntimeError("no foundation tables were built")

    combined_df = pd.concat(combined, ignore_index=True)
    manifest_df = pd.DataFrame(manifest_rows)

    combined_path = out_root / "spot_foundation_all_samples.csv"
    manifest_path = out_root / "foundation_manifest.csv"
    config_path = out_root / "foundation_config.json"

    combined_df.to_csv(combined_path, index=False)
    manifest_df.to_csv(manifest_path, index=False)
    config_path.write_text(
        json.dumps(
            {
                "samples": samples,
                "computed_samples": manifest_df["sample"].tolist() if not manifest_df.empty else [],
                "failures": failures,
                "methods": list(METHODS.keys()),
                "max_marker_genes": int(args.max_marker_genes),
                "low_support_quantile": float(args.low_support_quantile),
                "high_support_quantile": float(args.high_support_quantile),
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    print(f"[OK] wrote: {combined_path}")
    print(f"[OK] wrote: {manifest_path}")
    if failures:
        print(f"[WARN] {len(failures)} samples failed; see {config_path}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

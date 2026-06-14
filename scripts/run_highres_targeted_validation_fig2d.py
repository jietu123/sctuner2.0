from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.highres_enrichment_core import (
    add_generic_methods,
    evaluate_cytospace_pair,
    plot_fig2d,
    read_gene_sets,
    to_long,
)


def _evaluate_fixed_panel(
    root: Path,
    design: pd.DataFrame,
    nperm: int,
    seed: int,
) -> pd.DataFrame:
    mapping_manifest = pd.read_csv(
        root / "visualizations/highres_profile_mask_mapping/highres_profile_mask_mapping_manifest.csv"
    )
    gene_sets = read_gene_sets(root)
    records: list[dict[str, object]] = []

    for index, fixed in design.reset_index(drop=True).iterrows():
        sample = str(fixed["profile_mask_sample"])
        matches = mapping_manifest[
            mapping_manifest["profile_mask_sample"].astype(str).eq(sample)
        ]
        if len(matches) != 1:
            raise ValueError(f"Expected one mapping-manifest row for {sample}, found {len(matches)}.")
        gene_set_name = str(fixed["gene_set"])
        if gene_set_name not in gene_sets:
            raise KeyError(f"Unknown gene set: {gene_set_name}")
        result = evaluate_cytospace_pair(
            root,
            matches.iloc[0],
            str(fixed["readout_cell_type"]),
            gene_set_name,
            gene_sets[gene_set_name],
            nperm,
            seed + index * 100,
        )
        if result is None:
            raise RuntimeError(f"Fixed Fig2D design failed input requirements: {sample} / {gene_set_name}")
        record = {key: value for key, value in result.items() if key != "_curves"}
        record["candidate_label"] = (
            f"{record['raw_sample']} | {record['readout_cell_type']} | "
            f"{str(record['gene_set']).replace('EcoTyper_', '')}"
        )
        records.append(record)
    return pd.DataFrame(records)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the fixed high-resolution Fig2D targeted validation panel."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--design",
        default="configs/highres_targeted_validation_fig2d.csv",
    )
    parser.add_argument(
        "--out_dir",
        default="visualizations/highres_profile_mask_fig2d/targeted_validation",
    )
    parser.add_argument("--out_prefix", default="targeted_fig2d_highres_fixed_panel")
    parser.add_argument("--nperm", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=19)
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    design_path = root / args.design
    design = pd.read_csv(design_path)
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    selected = _evaluate_fixed_panel(root, design, args.nperm, args.seed)
    long_df = to_long(selected)
    long_df = add_generic_methods(root, selected, long_df, args.nperm, args.seed)

    selected_path = out_dir / f"{args.out_prefix}_fixed_profiles.csv"
    source_path = out_dir / f"{args.out_prefix}_source_values.csv"
    selected.to_csv(selected_path, index=False)
    long_df.to_csv(source_path, index=False)
    png = plot_fig2d(long_df, out_dir, args.out_prefix)

    baseline = selected["baseline_highres_NES"].astype(float)
    route2 = selected["route2_highres_NES"].astype(float)
    summary = {
        "experiment_id": "highres_targeted_validation_fig2d",
        "design_source": str(design_path.relative_to(root)).replace("\\", "/"),
        "n_fixed_profiles": int(len(selected)),
        "n_datasets": int(selected["raw_sample"].nunique()),
        "route2_improved_profiles": int((route2 > baseline).sum()),
        "mean_delta_nes": float((route2 - baseline).mean()),
        "median_delta_nes": float((route2 - baseline).median()),
        "mean_nes": long_df.groupby("method")["nes"].mean().to_dict(),
        "rows_per_method": long_df.groupby("method").size().to_dict(),
        "nperm": args.nperm,
        "seed": args.seed,
        "fixed_profiles_csv": str(selected_path.relative_to(root)).replace("\\", "/"),
        "source_values_csv": str(source_path.relative_to(root)).replace("\\", "/"),
        "png": str(png.relative_to(root)).replace("\\", "/"),
    }
    manifest_path = out_dir / f"{args.out_prefix}_manifest.json"
    manifest_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print("[OK] wrote:", manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

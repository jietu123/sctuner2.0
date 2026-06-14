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
    evaluate_cytospace_pair,
    plot_fig2c,
    read_gene_sets,
    slug,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the fixed high-resolution Fig2C targeted validation case study."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--design",
        default="configs/highres_targeted_validation_fig2c.json",
    )
    parser.add_argument(
        "--out_dir",
        default="visualizations/highres_targeted_validation_fig2c",
    )
    parser.add_argument("--nperm", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=11)
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    design_path = root / args.design
    design = json.loads(design_path.read_text(encoding="utf-8"))
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    mapping_manifest = pd.read_csv(
        root / "visualizations/highres_profile_mask_mapping/highres_profile_mask_mapping_manifest.csv"
    )
    matches = mapping_manifest[
        mapping_manifest["profile_mask_sample"].astype(str).eq(str(design["profile_mask_sample"]))
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one mapping-manifest row for {design['profile_mask_sample']}, found {len(matches)}."
        )

    gene_sets = read_gene_sets(root)
    gene_set_name = str(design["gene_set"])
    if gene_set_name not in gene_sets:
        raise KeyError(f"Unknown gene set: {gene_set_name}")

    result = evaluate_cytospace_pair(
        root,
        matches.iloc[0],
        str(design["readout_cell_type"]),
        gene_set_name,
        gene_sets[gene_set_name],
        args.nperm,
        args.seed,
    )
    if result is None:
        raise RuntimeError("The fixed Fig2C design did not meet the evaluation input requirements.")

    out_prefix = (
        f"targeted_fig2c_{slug(result['readout_cell_type'])}_"
        f"{slug(result['gene_set'])}_{slug(result['raw_sample'])}"
    )
    for suffix, curve in result["_curves"].items():
        curve.to_csv(out_dir / f"{out_prefix}_{suffix}_curve.csv", index=False)
    png, pdf = plot_fig2c(result, out_dir, out_prefix)

    summary = {key: value for key, value in result.items() if key != "_curves"}
    summary.update(
        {
            "experiment_id": design["experiment_id"],
            "design_source": str(design_path.relative_to(root)).replace("\\", "/"),
            "nperm": args.nperm,
            "seed": args.seed,
            "png": str(png.relative_to(root)).replace("\\", "/"),
            "pdf": str(pdf.relative_to(root)).replace("\\", "/"),
        }
    )
    manifest_path = out_dir / f"{out_prefix}_manifest.json"
    manifest_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print("[OK] wrote:", manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

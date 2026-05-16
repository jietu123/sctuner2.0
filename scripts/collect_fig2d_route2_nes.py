from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


PAIR_TO_PAPER = {
    "melanoma_slide1": ("Melanoma (Tirosh et al.)", "Legacy ST (slide 1, Thrane et al.)"),
    "melanoma_slide2": ("Melanoma (Tirosh et al.)", "Legacy ST (slide 2, Thrane et al.)"),
    "brca_er_her2_fresh_frozen": ("ER+HER2+ BRCA (Wu et al.)", "10x Visium (fresh-frozen)"),
    "brca_her2_ffpe": ("HER2+ BRCA (Wu et al.)", "10x Visium (FFPE)"),
    "brca_tnbc_fresh_frozen": ("TNBC BRCA (Wu et al.)", "10x Visium (fresh-frozen)"),
    "crc_fresh_frozen": ("CRC (Lee et al.)", "10x Visium (fresh-frozen)"),
}

DEFAULT_ROUTE2_SUMMARIES = {
    "melanoma_slide1": (
        "result/cytospace_fig2d_tme/melanoma_slide1/"
        "fig2d_enrichment_route2/fig2c_official_enrichment_summary.csv"
    ),
    "melanoma_slide2": (
        "result/cytospace_fig2d_tme/melanoma_slide2/"
        "fig2d_enrichment_route2/fig2c_official_enrichment_summary.csv"
    ),
    "brca_er_her2_fresh_frozen": (
        "result/cytospace_fig2d_tme/brca_er_her2_fresh_frozen/"
        "fig2d_enrichment_route2/fig2c_official_enrichment_summary.csv"
    ),
    "brca_her2_ffpe": (
        "result/cytospace_fig2d_tme/brca_her2_ffpe/"
        "fig2d_enrichment_route2/fig2c_official_enrichment_summary.csv"
    ),
    "brca_tnbc_fresh_frozen": (
        "result/cytospace_fig2d_tme/brca_tnbc_fresh_frozen/"
        "fig2d_enrichment_route2/fig2c_official_enrichment_summary.csv"
    ),
    "crc_fresh_frozen": (
        "result/cytospace_fig2d_tme/crc_fresh_frozen/"
        "fig2d_enrichment_route2/fig2c_official_enrichment_summary.csv"
    ),
}

EXPECTED_ROUTE2_ROWS_PER_PAIR = 2


def _parse_mapping(items: list[str]) -> dict[str, str]:
    mapping = dict(DEFAULT_ROUTE2_SUMMARIES)
    for item in items:
        if "=" not in item:
            raise ValueError(f"--route2_summary must be pair_id=path, got: {item}")
        pair_id, path = item.split("=", 1)
        mapping[pair_id.strip()] = path.strip()
    return mapping


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Collect route2 T-cell exhaustion NES values into Fig.2d-compatible source-data format."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--route2_summary",
        action="append",
        default=[],
        help="Optional pair_id=summary_csv. Defaults include the existing melanoma slide 1 route2 pilot.",
    )
    parser.add_argument(
        "--out_csv",
        default="result/cytospace_fig2d_tme/route2_fig2d_nes_values.csv",
    )
    parser.add_argument(
        "--allow_partial",
        action="store_true",
        help="Write available route2 rows even if not all six Fig.2d pairs have CD4/CD8 NES values.",
    )
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    rows: list[dict[str, object]] = []
    missing_or_incomplete: dict[str, str] = {}
    for pair_id, rel_path in _parse_mapping(args.route2_summary).items():
        if pair_id not in PAIR_TO_PAPER:
            raise KeyError(f"Unknown pair_id: {pair_id}")
        path = (root / rel_path).resolve()
        if not path.exists():
            print(f"[WARN] route2 summary missing for {pair_id}: {path}")
            missing_or_incomplete[pair_id] = f"missing summary: {path}"
            continue
        df = pd.read_csv(path)
        required = {"cell_type", "pathway", "NES", "pval"}
        missing = required.difference(df.columns)
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        scrna, st = PAIR_TO_PAPER[pair_id]
        sub = df[
            df["cell_type"].astype(str).isin(["CD4 T cells", "CD8 T cells"])
            & df["pathway"].astype(str).eq("Exhaustion")
        ].copy()
        if len(sub) != EXPECTED_ROUTE2_ROWS_PER_PAIR:
            missing_or_incomplete[pair_id] = (
                f"expected {EXPECTED_ROUTE2_ROWS_PER_PAIR} CD4/CD8 exhaustion rows, found {len(sub)}"
            )
        for _, row in sub.iterrows():
            rows.append(
                {
                    "pair_id": pair_id,
                    "scrna": scrna,
                    "st": st,
                    "method": "SVTuner + CytoSPACE",
                    "cell_type": str(row["cell_type"]),
                    "feature": "Exhaustion",
                    "pathway": f"{row['cell_type']} exhaustion (Zheng et al.)",
                    "nes": float(row["NES"]),
                    "pval": float(row["pval"]),
                    "source_summary": str(path.relative_to(root)),
                }
            )

    expected_total = len(PAIR_TO_PAPER) * EXPECTED_ROUTE2_ROWS_PER_PAIR
    if (missing_or_incomplete or len(rows) != expected_total) and not args.allow_partial:
        details = "\n".join(f"- {pair_id}: {reason}" for pair_id, reason in missing_or_incomplete.items())
        raise SystemExit(
            "Incomplete Fig.2d route2 coverage. "
            f"Expected {expected_total} rows ({len(PAIR_TO_PAPER)} pairs x CD4/CD8), got {len(rows)}.\n"
            f"{details}\n"
            "Use --allow_partial only for debugging/pilot figures."
        )

    out = (root / args.out_csv).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    out_df = pd.DataFrame(rows)
    out_df.to_csv(out, index=False)
    print(f"[OK] wrote: {out}")
    print(f"[INFO] route2 rows: {len(out_df)}")
    if missing_or_incomplete:
        print("[WARN] incomplete pairs:")
        for pair_id, reason in missing_or_incomplete.items():
            print(f"  - {pair_id}: {reason}")
    if not out_df.empty:
        print(out_df[["pair_id", "cell_type", "nes", "pval"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

import json
import csv
from pathlib import Path

SEEDS = [0, 1, 2, 42]

def load_json(path: Path):
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)

def main():
    rows = []
    for seed in SEEDS:
        base_tag = f"real_brca_simS0_seed{seed}"
        root = Path("result") / base_tag / "stage5_route2_s0" / "default"
        b_path = root / "stage5_route2_s0__baseline.json"
        r_path = root / "stage5_route2_s0__route2.json"
        if not b_path.exists() or not r_path.exists():
            print(f"[skip] seed {seed} missing stage5 outputs")
            continue
        b = load_json(b_path)
        r = load_json(r_path)
        filtered_fraction = None
        try:
            n_before = r["coverage"]["n_cells_before_prefilter"]
            n_filtered = r["coverage"]["n_filtered"]
            filtered_fraction = n_filtered / n_before if n_before else None
        except Exception:
            filtered_fraction = None
        row = {
            "seed": seed,
            "leak_rate_baseline": b["leakage"]["missing_leak_rate"],
            "leak_rate_route2": r["leakage"]["missing_leak_rate"],
            "L1_baseline": b["composition"]["L1_mean"],
            "L1_route2": r["composition"]["L1_mean"],
            "delta_L1": r["composition"]["L1_mean"] - b["composition"]["L1_mean"],
            "JS_baseline": b["composition"]["JS_mean"],
            "JS_route2": r["composition"]["JS_mean"],
            "delta_JS": r["composition"]["JS_mean"] - b["composition"]["JS_mean"],
            "corr_baseline": b["composition"]["corr_mean"],
            "corr_route2": r["composition"]["corr_mean"],
            "delta_corr": r["composition"]["corr_mean"] - b["composition"]["corr_mean"],
            "filtered_fraction_route2": filtered_fraction,
            "sha1_truth_spot": r["sha1"]["truth_spot_type_fraction"],
            "sha1_truth_query": r["sha1"]["truth_query_cell_spot"],
        }
        rows.append(row)

    out = Path("result") / "real_brca_simS0_seed_summary.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[OK] wrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()


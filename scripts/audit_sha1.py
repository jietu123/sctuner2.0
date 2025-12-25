import pandas as pd
from pathlib import Path


def main():
    p = Path("result/real_brca_simS0_seed_summary.csv")
    if not p.exists():
        print(f"[ERR] seed_summary not found: {p}")
        return

    df = pd.read_csv(p)
    for col in ["sha1_truth_spot", "sha1_truth_query"]:
        if col not in df.columns:
            print(f"\n{col}: column missing")
            continue
        print(f"\n{col}")
        print("unique =", df[col].nunique())
        print(df[["seed", col]].sort_values(col).to_string(index=False))


if __name__ == "__main__":
    main()


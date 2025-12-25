import argparse
import pandas as pd
from pathlib import Path


def transpose_sample(sample: str):
    base = Path("data/processed") / sample / "stage1_preprocess" / "exported"
    sc_path = base / "sc_expression_normalized.csv"
    st_path = base / "st_expression_normalized.csv"

    sc = pd.read_csv(sc_path, sep=None, engine="python")
    sc_t = sc.set_index(sc.columns[0]).T
    sc_t.to_csv(sc_path)
    print(f"[{sample}] sc_expr transposed:", sc_t.shape)

    st = pd.read_csv(st_path, sep=None, engine="python")
    st_t = st.set_index(st.columns[0]).T
    st_t.to_csv(st_path)
    print(f"[{sample}] st_expr transposed:", st_t.shape)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True, help="sample name under data/processed/")
    args = ap.parse_args()
    transpose_sample(args.sample)


if __name__ == "__main__":
    main()


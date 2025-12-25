import argparse
import pandas as pd
from pathlib import Path


def transpose(seed: int):
    base = Path(f"data/processed/real_brca_simS0_seed{seed}/stage1_preprocess/exported")
    sc_path = base / "sc_expression_normalized.csv"
    st_path = base / "st_expression_normalized.csv"

    sc = pd.read_csv(sc_path, sep=None, engine="python")
    sc_t = sc.set_index(sc.columns[0]).T
    sc_t.to_csv(sc_path)
    print(f"[seed{seed}] sc_expr transposed:", sc_t.shape)

    st = pd.read_csv(st_path, sep=None, engine="python")
    st_t = st.set_index(st.columns[0]).T
    st_t.to_csv(st_path)
    print(f"[seed{seed}] st_expr transposed:", st_t.shape)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, required=True)
    args = ap.parse_args()
    transpose(args.seed)


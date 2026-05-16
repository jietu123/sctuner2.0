from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from scipy import io, sparse

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.stages.stage1_io import load_stage1


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Rebuild CytoSPACE-style assigned_expression from Stage1 sc expression and assigned_locations.csv."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--cytospace_output", required=True)
    parser.add_argument("--sc_expr_source", default="normalized", choices=["normalized", "data", "counts", "auto"])
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    out_dir = Path(args.cytospace_output).resolve()
    assigned_path = out_dir / "assigned_locations.csv"
    if not assigned_path.exists():
        raise FileNotFoundError(assigned_path)

    assigned = pd.read_csv(assigned_path)
    required = {"UniqueCID", "OriginalCID"}
    missing = required - set(assigned.columns)
    if missing:
        raise ValueError(f"{assigned_path} missing columns: {sorted(missing)}")

    sc_expr, _, _, _ = load_stage1(root, args.sample, sc_expr_source=args.sc_expr_source)
    sc_expr.index = sc_expr.index.astype(str)
    original_ids = assigned["OriginalCID"].astype(str).tolist()
    missing_ids = sorted(set(original_ids) - set(sc_expr.index))
    if missing_ids:
        raise ValueError(f"{len(missing_ids)} assigned OriginalCID values missing from Stage1 sc expression, e.g. {missing_ids[:5]}")

    # Matrix Market convention used by CytoSPACE: genes x assigned cells.
    selected = sc_expr.loc[original_ids]
    mat = sparse.csr_matrix(selected.to_numpy(dtype="float32").T)

    expr_dir = out_dir / "assigned_expression"
    expr_dir.mkdir(parents=True, exist_ok=True)
    io.mmwrite(expr_dir / "matrix.mtx", mat)
    pd.Series(sc_expr.columns.astype(str)).to_csv(expr_dir / "genes.tsv", index=False, header=False)
    pd.Series(assigned["UniqueCID"].astype(str)).to_csv(expr_dir / "barcodes.tsv", index=False, header=False)

    print(f"[OK] wrote assigned_expression: {expr_dir}")
    print(f"[INFO] matrix genes x assigned_cells = {mat.shape}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

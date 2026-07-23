from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.build_assigned_expression_from_stage1 import main as _unused  # noqa: F401
from src.stages.stage1_io import load_stage1
from scipy import io, sparse


def _coord_columns(coords: pd.DataFrame) -> tuple[str, str]:
    if {"X", "Y"}.issubset(coords.columns):
        return "X", "Y"
    if {"row", "col"}.issubset(coords.columns):
        return "row", "col"
    if {"spatial_row", "spatial_col"}.issubset(coords.columns):
        return "spatial_row", "spatial_col"
    raise ValueError("st_coordinates.csv must contain X/Y, row/col, or spatial_row/spatial_col columns")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert generic cell_assignment.csv into CytoSPACE-like assigned_locations/assigned_expression."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--mapping_dir", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--sc_expr_source", default="normalized", choices=["normalized", "data", "counts", "auto"])
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    mapping_dir = Path(args.mapping_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    assignment_path = mapping_dir / "cell_assignment.csv"
    if not assignment_path.exists():
        raise FileNotFoundError(assignment_path)
    assign = pd.read_csv(assignment_path)
    required = {"cell_id", "assigned_spot", "cell_type"}
    missing = required.difference(assign.columns)
    if missing:
        raise ValueError(f"{assignment_path} missing columns: {sorted(missing)}")

    sc_expr, _, st_coords, _ = load_stage1(root, args.sample, sc_expr_source=args.sc_expr_source)
    sc_expr.index = sc_expr.index.astype(str)
    assign = assign.copy()
    assign["cell_id"] = assign["cell_id"].astype(str)
    assign["assigned_spot"] = assign["assigned_spot"].astype(str)
    assign["cell_type"] = assign["cell_type"].astype(str)
    assign = assign[assign["cell_id"].isin(sc_expr.index)].copy()
    if assign.empty:
        raise ValueError("no assigned cells overlap Stage1 sc expression")

    coords = st_coords.copy()
    if "spot_id" not in coords.columns:
        coords = coords.reset_index(names="spot_id")
    coords["spot_id"] = coords["spot_id"].astype(str)
    x_col, y_col = _coord_columns(coords)
    coords = coords.set_index("spot_id")

    missing_spots = sorted(set(assign["assigned_spot"]) - set(coords.index))
    if missing_spots:
        raise ValueError(f"{len(missing_spots)} assigned spots missing coordinates, e.g. {missing_spots[:5]}")

    loc = pd.DataFrame(
        {
            "UniqueCID": assign["cell_id"].to_numpy(),
            "OriginalCID": assign["cell_id"].to_numpy(),
            "CellType": assign["cell_type"].to_numpy(),
            "SpotID": assign["assigned_spot"].to_numpy(),
        }
    )
    loc["X"] = coords.loc[loc["SpotID"], x_col].to_numpy()
    loc["Y"] = coords.loc[loc["SpotID"], y_col].to_numpy()
    loc.to_csv(out_dir / "assigned_locations.csv", index=False)

    selected = sc_expr.loc[loc["OriginalCID"].astype(str)]
    mat = sparse.csr_matrix(selected.to_numpy(dtype="float32").T)
    expr_dir = out_dir / "assigned_expression"
    expr_dir.mkdir(parents=True, exist_ok=True)
    io.mmwrite(expr_dir / "matrix.mtx", mat)
    pd.Series(sc_expr.columns.astype(str)).to_csv(expr_dir / "genes.tsv", index=False, header=False)
    pd.Series(loc["UniqueCID"].astype(str)).to_csv(expr_dir / "barcodes.tsv", index=False, header=False)

    print(f"[OK] wrote CytoSPACE-like output: {out_dir}")
    print(f"[INFO] assigned cells: {len(loc)}")
    print(f"[INFO] assigned_expression genes x cells = {mat.shape}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

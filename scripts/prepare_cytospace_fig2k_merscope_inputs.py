from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


RAW_SOURCE_DEFAULT = "data/raw/cytospace_fig2k_breast_merscope/HumanBreastCancerPatient1"
PROCESSED_DEFAULT = "data/processed/cytospace_fig2k_breast_merscope/HumanBreastCancerPatient1"

PAN_T_GENES = ["CD3D", "CD3E", "CD3G", "TRAC", "PTPRC"]
CD4_GENES = ["CD4", "IL7R", "CCR7", "ICOS", "CTLA4"]
CD8_GENES = ["CD8A", "CD8B", "GZMB", "NKG7", "PRF1"]
EXHAUSTION_GENES = ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "TOX", "IFNG", "GZMB"]
TUMOR_PROXIMITY_GENES = ["EPCAM", "ERBB2", "KRT8", "KRT18", "KRT19", "MUC1", "CDH1"]


def _header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return next(csv.reader(handle))


def _count_rows(path: Path) -> int:
    n = 0
    with path.open("rb") as handle:
        for _ in handle:
            n += 1
    return max(n - 1, 0)


def _available_genes(expr_path: Path, genes: list[str]) -> list[str]:
    columns = set(_header(expr_path))
    return [g for g in genes if g in columns]


def _score(df: pd.DataFrame, genes: list[str]) -> pd.Series:
    if not genes:
        return pd.Series(0.0, index=df.index)
    return np.log1p(df[genes].astype("float32")).mean(axis=1)


def _build_candidate_table(
    expr_path: Path,
    meta_path: Path,
    out_path: Path,
    chunk_size: int,
    min_t_score_quantile: float,
) -> dict[str, object]:
    gene_groups = {
        "pan_t": _available_genes(expr_path, PAN_T_GENES),
        "cd4": _available_genes(expr_path, CD4_GENES),
        "cd8": _available_genes(expr_path, CD8_GENES),
        "exhaustion": _available_genes(expr_path, EXHAUSTION_GENES),
        "tumor_proximity": _available_genes(expr_path, TUMOR_PROXIMITY_GENES),
    }
    selected_genes = sorted({g for genes in gene_groups.values() for g in genes})
    if not selected_genes:
        raise ValueError("No marker genes found in MERSCOPE expression header.")

    usecols = ["cell", *selected_genes]
    score_chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(expr_path, usecols=usecols, chunksize=chunk_size):
        chunk["cell_id"] = chunk["cell"].astype(str)
        chunk = chunk.set_index("cell_id", drop=False)
        pan_t = _score(chunk, gene_groups["pan_t"])
        cd4 = _score(chunk, gene_groups["cd4"])
        cd8 = _score(chunk, gene_groups["cd8"])
        exhaustion = _score(chunk, gene_groups["exhaustion"])
        tumor = _score(chunk, gene_groups["tumor_proximity"])
        score_chunks.append(
            pd.DataFrame(
                {
                    "cell_id": chunk["cell_id"].to_numpy(),
                    "pan_t_score": pan_t.to_numpy(dtype=float),
                    "cd4_score": cd4.to_numpy(dtype=float),
                    "cd8_score": cd8.to_numpy(dtype=float),
                    "exhaustion_score": exhaustion.to_numpy(dtype=float),
                    "tumor_proximity_score": tumor.to_numpy(dtype=float),
                }
            )
        )
    scores = pd.concat(score_chunks, ignore_index=True)
    t_threshold = float(scores["pan_t_score"].quantile(min_t_score_quantile))
    scores["is_t_candidate"] = scores["pan_t_score"].ge(t_threshold)
    scores["is_cd4_t_candidate"] = scores["is_t_candidate"] & scores["cd4_score"].ge(scores["cd8_score"])

    meta = pd.read_csv(meta_path)
    id_col = meta.columns[0]
    meta = meta.rename(columns={id_col: "cell_id"})
    meta["cell_id"] = meta["cell_id"].astype(str)
    keep_cols = ["cell_id", "fov", "volume", "center_x", "center_y", "min_x", "max_x", "min_y", "max_y"]
    meta = meta[[c for c in keep_cols if c in meta.columns]]

    merged = scores.merge(meta, on="cell_id", how="left")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, index=False)

    return {
        "marker_gene_groups": gene_groups,
        "n_rows": int(len(merged)),
        "pan_t_score_quantile": float(min_t_score_quantile),
        "pan_t_score_threshold": t_threshold,
        "n_t_candidates": int(merged["is_t_candidate"].sum()),
        "n_cd4_t_candidates": int(merged["is_cd4_t_candidate"].sum()),
        "candidate_table": str(out_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare downloaded Vizgen MERSCOPE breast cancer inputs for Fig.2k work.")
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--raw_dir", default=RAW_SOURCE_DEFAULT)
    parser.add_argument("--out_dir", default=PROCESSED_DEFAULT)
    parser.add_argument("--chunk_size", type=int, default=50000)
    parser.add_argument("--min_t_score_quantile", type=float, default=0.90)
    parser.add_argument("--copy_raw", action="store_true", help="Copy raw CSVs into out_dir/raw instead of only referencing them.")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    raw_dir = (root / args.raw_dir).resolve()
    out_dir = (root / args.out_dir).resolve()
    expr_path = raw_dir / "cell_by_gene.csv"
    meta_path = raw_dir / "cell_metadata.csv"
    if not expr_path.exists():
        raise FileNotFoundError(expr_path)
    if not meta_path.exists():
        raise FileNotFoundError(meta_path)

    out_dir.mkdir(parents=True, exist_ok=True)
    raw_out = out_dir / "raw"
    if args.copy_raw:
        raw_out.mkdir(parents=True, exist_ok=True)
        for src in [expr_path, meta_path]:
            dst = raw_out / src.name
            if not dst.exists():
                shutil.copy2(src, dst)

    expr_header = _header(expr_path)
    meta_header = _header(meta_path)
    candidate_path = out_dir / "prepared" / "merscope_tcell_candidate_scores.csv"
    candidate_summary = _build_candidate_table(
        expr_path=expr_path,
        meta_path=meta_path,
        out_path=candidate_path,
        chunk_size=args.chunk_size,
        min_t_score_quantile=args.min_t_score_quantile,
    )

    manifest = {
        "raw_dir": str(raw_dir),
        "processed_dir": str(out_dir),
        "source_files": {
            "cell_by_gene": str(expr_path),
            "cell_metadata": str(meta_path),
        },
        "raw_copy_dir": str(raw_out) if args.copy_raw else None,
        "dimensions": {
            "cell_by_gene_rows": _count_rows(expr_path),
            "cell_by_gene_columns": len(expr_header),
            "cell_metadata_rows": _count_rows(meta_path),
            "cell_metadata_columns": len(meta_header),
        },
        "headers": {
            "cell_by_gene_first_columns": expr_header[:25],
            "cell_metadata_columns": meta_header,
        },
        "candidate_summary": candidate_summary,
        "notes": [
            "This prepares the MERSCOPE spatial single-cell side only.",
            "Strict CytoSPACE/route2 Fig.2k mapping still requires the Wu et al. BRCA scRNA-seq reference.",
            "The candidate table is for sanity checks and forced-unsupported scenario design; it is not a replacement for mapping output.",
        ],
    }
    manifest_path = out_dir / "fig2k_merscope_input_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[OK] wrote: {candidate_path}")
    print(f"[OK] wrote: {manifest_path}")
    print(json.dumps(candidate_summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw" / "cytospace_fig2c_melanoma"
GROUP = "cytospace_fig2c_melanoma"
SCRNA_FILE = "GSE72056_melanoma_single_cell_revised_v2.txt"
ST_REPLICATES = {
    "mel1_rep1": "berglund2018spatial_ST_mel1_rep1.h5ad",
    "mel1_rep2": "berglund2018spatial_ST_mel1_rep2.h5ad",
    "mel2_rep1": "berglund2018spatial_ST_mel2_rep1.h5ad",
}

NON_MALIGNANT_LABELS = {
    "1": "T cells",
    "2": "B cells",
    "3": "Macrophages",
    "4": "Endothelial cells",
    "5": "CAFs",
    "6": "NK cells",
}


def _gene_symbol(value: object) -> str:
    return str(value).strip().split()[0]


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare CytoSPACE Fig.2c melanoma data as SVTuner Stage1 exports."
    )
    p.add_argument("--project_root", default=str(PROJECT_ROOT))
    p.add_argument("--raw_dir", default=None)
    p.add_argument(
        "--max_genes",
        type=int,
        default=0,
        help="Optional cap on common genes by variance across ST spots; 0 keeps all common genes.",
    )
    p.add_argument("--overwrite", action="store_true")
    p.add_argument(
        "--no_split_t_cells",
        action="store_true",
        help="Keep broad T cells instead of marker-based CD4/CD8 split.",
    )
    p.add_argument(
        "--force_unsupported_type",
        action="append",
        default=[],
        help="Optional cell type to force as unsupported in the generated Stage3 config. Can be repeated.",
    )
    return p.parse_args()


def _resolve_st_path(project_root: Path, raw_dir: Path, filename: str) -> Path:
    candidates = [
        raw_dir / filename,
        project_root / "data" / "raw" / "cytospace_fig2d_tme" / "melanoma_ST_Thrane_legacyST" / filename,
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(candidates[0])


def _read_scrna_metadata_and_gene_index(path: Path) -> tuple[list[str], pd.DataFrame, list[str]]:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        tumor = handle.readline().rstrip("\n").split("\t")
        malignant = handle.readline().rstrip("\n").split("\t")
        non_malignant = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]

        genes: list[str] = []
        for line in handle:
            label = line.split("\t", 1)[0].strip()
            if label:
                genes.append(_gene_symbol(label))

    malignant_vals = malignant[1:]
    non_malignant_vals = non_malignant[1:]
    tumor_vals = tumor[1:]
    cell_types: list[str] = []
    for mal, nm in zip(malignant_vals, non_malignant_vals):
        if str(mal) == "2":
            cell_types.append("Melanoma")
        elif str(mal) == "1":
            cell_types.append(NON_MALIGNANT_LABELS.get(str(nm), "Non-malignant unresolved"))
        else:
            cell_types.append("Unresolved")

    meta = pd.DataFrame(
        {
            "cell_id": cell_ids,
            "tumor": tumor_vals,
            "malignant_code": malignant_vals,
            "non_malignant_code": non_malignant_vals,
            "cell_type": cell_types,
        }
    )
    return cell_ids, meta, genes


def _st_expression_and_coords(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    adata = ad.read_h5ad(path)
    genes = [_gene_symbol(name) for name in adata.var_names]
    arr = adata.X
    if hasattr(arr, "toarray"):
        arr = arr.toarray()
    st_expr = pd.DataFrame(np.asarray(arr, dtype=np.float32), index=adata.obs_names.astype(str), columns=genes)
    # Collapse duplicate gene symbols, if any, by summing their expression.
    st_expr = st_expr.T.groupby(level=0).sum().T

    spatial = adata.obsm["spatial"]
    if isinstance(spatial, pd.DataFrame):
        coords = spatial.copy()
        coords.index = adata.obs_names.astype(str)
        coord_values = coords.iloc[:, :2].to_numpy()
    else:
        coord_values = np.asarray(spatial)[:, :2]
    st_coords = pd.DataFrame(coord_values, index=adata.obs_names.astype(str), columns=["row", "col"])
    st_coords["orig.ident"] = path.stem
    st_coords["nCount_RNA"] = st_expr.sum(axis=1).to_numpy()
    st_coords["nFeature_RNA"] = (st_expr > 0).sum(axis=1).to_numpy()
    return st_expr, st_coords


def _select_common_genes(
    sc_genes: list[str],
    st_expr: pd.DataFrame,
    max_genes: int,
) -> list[str]:
    sc_set = set(sc_genes)
    common = [g for g in st_expr.columns.astype(str).tolist() if g in sc_set]
    common = sorted(set(common))
    if max_genes and len(common) > max_genes:
        var = st_expr[common].var(axis=0).sort_values(ascending=False)
        common = sorted(var.head(max_genes).index.astype(str).tolist())
    return common


def _stream_scrna_expression(path: Path, common_genes: list[str], cell_ids: list[str]) -> pd.DataFrame:
    wanted = set(common_genes)
    rows: dict[str, np.ndarray] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for _ in range(4):
            next(handle)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            gene = _gene_symbol(parts[0])
            if gene not in wanted:
                continue
            values = np.asarray(parts[1:], dtype=np.float32)
            if gene in rows:
                rows[gene] = rows[gene] + values
            else:
                rows[gene] = values
    missing = sorted(wanted - set(rows))
    if missing:
        raise ValueError(f"Missing {len(missing)} expected genes in scRNA matrix, e.g. {missing[:5]}")
    sc_expr = pd.DataFrame.from_dict(rows, orient="index", columns=cell_ids)
    return sc_expr.T.reindex(columns=common_genes).astype(np.float32)


def _split_t_cells_by_markers(sc_expr: pd.DataFrame, sc_meta: pd.DataFrame) -> pd.DataFrame:
    out = sc_meta.copy()
    t_mask = out["cell_type"].astype(str) == "T cells"
    if not t_mask.any():
        return out
    cd4_markers = [g for g in ["CD4", "IL7R", "CCR7"] if g in sc_expr.columns]
    cd8_markers = [g for g in ["CD8A", "CD8B", "GZMB"] if g in sc_expr.columns]
    if not cd4_markers or not cd8_markers:
        return out
    t_ids = out.loc[t_mask, "cell_id"].astype(str)
    expr = sc_expr.reindex(index=t_ids)
    cd4_score = expr[cd4_markers].mean(axis=1)
    cd8_score = expr[cd8_markers].mean(axis=1)
    max_score = pd.concat([cd4_score, cd8_score], axis=1).max(axis=1)
    labels = np.where(
        max_score <= 0,
        "T cells unclassified",
        np.where(cd4_score > cd8_score, "CD4 T cells", "CD8 T cells"),
    )
    label_map = dict(zip(t_ids, labels))
    out.loc[t_mask, "cell_type"] = out.loc[t_mask, "cell_id"].astype(str).map(label_map).to_numpy()
    out["tcell_cd4_marker_score"] = out["cell_id"].astype(str).map(cd4_score.to_dict())
    out["tcell_cd8_marker_score"] = out["cell_id"].astype(str).map(cd8_score.to_dict())
    return out


def _write_dataset_config(
    project_root: Path,
    sample: str,
    hvg_path: Path,
    force_unsupported_types: list[str] | None = None,
) -> None:
    cfg = {
        "paths": {
            "sc_expr": "GSE72056_melanoma_single_cell_revised_v2.txt",
            "sc_meta": None,
            "st_expr": None,
            "st_meta": None,
            "svg_marker_whitelist": None,
        },
        "qc": {
            "sc_min_genes": 0,
            "sc_max_genes": 100000,
            "sc_max_mt": 100,
            "st_min_genes": 0,
            "st_max_genes": "Inf",
            "st_max_mt": 100,
            "hvg_nfeatures": 2000,
            "mt_pattern": "^(MT-|mt-)",
        },
        "gene_filter": {"min_cells_sc": 0, "min_cells_st": 0},
        "stage3": {
            "strong_th": 0.7,
            "weak_th": 0.4,
            "st_cluster_k": 0,
            "unknown_floor": 0.3,
            "min_cells_rare_type": 20,
            "eps": 1.0e-8,
            "plugin_genes_path": str(hvg_path.relative_to(project_root)).replace("\\", "/"),
            "gene_weights_path": None,
            "auto_missing_detection": {
                "enable": False,
                "method": "adaptive_low_support",
                "action": "mark_unknown",
            },
            "masked_missing_detection": {"enable": False},
            "marker_identity_diagnostics": {"enable": True, "marker_top_n": 80},
        },
        "storage": {"group": GROUP},
    }
    if force_unsupported_types:
        cfg["stage3"]["force_unsupported_types"] = [str(x) for x in force_unsupported_types if str(x).strip()]
    cfg_path = project_root / "configs" / "datasets" / f"{sample}.yaml"
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")


def _write_stage1(
    project_root: Path,
    sample: str,
    sc_expr: pd.DataFrame,
    sc_meta: pd.DataFrame,
    st_expr: pd.DataFrame,
    st_coords: pd.DataFrame,
    common_genes: list[str],
    overwrite: bool,
    force_unsupported_types: list[str] | None = None,
) -> None:
    stage1_dir = project_root / "data" / "processed" / GROUP / sample / "stage1_preprocess"
    export = stage1_dir / "exported"
    if export.exists() and not overwrite:
        raise FileExistsError(f"{export} exists; pass --overwrite to rebuild")
    export.mkdir(parents=True, exist_ok=True)

    sc_expr = sc_expr.reindex(columns=common_genes)
    st_expr = st_expr.reindex(columns=common_genes)
    sc_expr.to_csv(export / "sc_expression_normalized.csv", index_label="cell_id")
    st_expr.to_csv(export / "st_expression_normalized.csv", index_label="spot_id")
    sc_meta.to_csv(export / "sc_metadata.csv", index=False)
    st_coords.to_csv(export / "st_coordinates.csv", index_label="spot_id")

    (stage1_dir / "common_genes.txt").write_text("\n".join(common_genes) + "\n", encoding="utf-8")
    (stage1_dir / "hvg_genes.txt").write_text("\n".join(common_genes) + "\n", encoding="utf-8")

    summary = {
        "sample": sample,
        "n_cells": int(sc_expr.shape[0]),
        "n_spots": int(st_expr.shape[0]),
        "n_common_genes": int(len(common_genes)),
        "cell_type_counts": sc_meta["cell_type"].value_counts().to_dict(),
        "outputs": {
            "stage1_exported": str(export),
            "common_genes": str(stage1_dir / "common_genes.txt"),
            "hvg_genes": str(stage1_dir / "hvg_genes.txt"),
        },
    }
    (stage1_dir / "fig2c_stage1_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _write_dataset_config(project_root, sample, stage1_dir / "hvg_genes.txt", force_unsupported_types)


def main() -> int:
    args = _parse_args()
    project_root = Path(args.project_root).resolve()
    raw_dir = Path(args.raw_dir).resolve() if args.raw_dir else project_root / RAW_DIR.relative_to(PROJECT_ROOT)
    scrna_path = raw_dir / SCRNA_FILE
    if not scrna_path.exists():
        raise FileNotFoundError(scrna_path)

    cell_ids, sc_meta, sc_genes = _read_scrna_metadata_and_gene_index(scrna_path)
    manifest: dict[str, object] = {
        "raw_dir": str(raw_dir),
        "scrna": {
            "n_cells": len(cell_ids),
            "n_gene_symbols": len(set(sc_genes)),
            "cell_type_counts": sc_meta["cell_type"].value_counts().to_dict(),
        },
        "samples": {},
    }

    for rep, filename in ST_REPLICATES.items():
        sample = f"cytospace_fig2c_melanoma_{rep}"
        st_path = _resolve_st_path(project_root, raw_dir, filename)
        st_expr, st_coords = _st_expression_and_coords(st_path)
        common_genes = _select_common_genes(sc_genes, st_expr, int(args.max_genes))
        if not common_genes:
            raise ValueError(f"No common genes for {sample}")
        sc_expr = _stream_scrna_expression(scrna_path, common_genes, cell_ids)
        sc_meta_out = sc_meta if args.no_split_t_cells else _split_t_cells_by_markers(sc_expr, sc_meta)
        _write_stage1(
            project_root=project_root,
            sample=sample,
            sc_expr=sc_expr,
            sc_meta=sc_meta_out,
            st_expr=st_expr,
            st_coords=st_coords,
            common_genes=common_genes,
            overwrite=bool(args.overwrite),
            force_unsupported_types=list(args.force_unsupported_type or []),
        )
        manifest["samples"][sample] = {
            "st_file": str(st_path),
            "n_spots": int(st_expr.shape[0]),
            "n_st_genes": int(st_expr.shape[1]),
            "n_common_genes": int(len(common_genes)),
            "stage1_dir": str(project_root / "data" / "processed" / GROUP / sample / "stage1_preprocess"),
            "config": str(project_root / "configs" / "datasets" / f"{sample}.yaml"),
        }
        print(f"[OK] prepared {sample}: cells={len(cell_ids)} spots={st_expr.shape[0]} genes={len(common_genes)}")

    out_manifest = raw_dir / "prepared" / "fig2c_stage1_manifest.json"
    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    out_manifest.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[OK] wrote: {out_manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

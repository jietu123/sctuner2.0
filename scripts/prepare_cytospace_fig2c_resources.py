from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw" / "cytospace_fig2c_melanoma"
OUT_DIR = RAW_DIR / "prepared"

SUPPLEMENTARY_XLSX = "41587_2023_1697_MOESM3_ESM.xlsx"
SCRNA_FILE = "GSE72056_melanoma_single_cell_revised_v2.txt"
ST_FILES = [
    "berglund2018spatial_ST_mel1_rep1.h5ad",
    "berglund2018spatial_ST_mel1_rep2.h5ad",
]
EXHAUSTION_COLUMN = "Tcell_exhaustion_Zheng_etal_Cell2017"


def _extract_table_s6_gene_sets(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name="Table S6", header=None)
    header_row = None
    for idx, row in raw.iterrows():
        values = [str(v).strip() for v in row.tolist() if pd.notna(v)]
        if EXHAUSTION_COLUMN in values:
            header_row = idx
            break
    if header_row is None:
        raise ValueError(f"Could not find {EXHAUSTION_COLUMN} in Table S6")

    table = raw.iloc[header_row + 1 :].copy()
    table.columns = [str(c).strip() for c in raw.iloc[header_row].tolist()]
    table = table.loc[:, [c for c in table.columns if c and c != "nan"]]
    return table


def _clean_gene_list(series: pd.Series) -> list[str]:
    genes: list[str] = []
    seen = set()
    for value in series.dropna().astype(str):
        gene = value.strip()
        if not gene or gene.lower() == "nan":
            continue
        if gene not in seen:
            genes.append(gene)
            seen.add(gene)
    return genes


def _scrna_header_info(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        meta_rows = [handle.readline().rstrip("\n").split("\t")[0] for _ in range(3)]
    return {
        "n_columns_including_label": len(header),
        "n_cells": max(len(header) - 1, 0),
        "first_cell_ids": header[1:6],
        "first_metadata_rows": meta_rows,
    }


def _gene_symbol(value: object) -> str:
    return str(value).strip().split()[0]


def _scrna_gene_overlap(path: Path, genes: list[str]) -> dict[str, object]:
    query = set(genes)
    observed: set[str] = set()
    n_gene_rows = 0
    metadata_labels = {
        "tumor",
        '"malignant(1=no,2=yes,0=unresolved)"',
        '"non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)"',
    }
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        next(handle, None)
        for line in handle:
            label = line.rstrip("\n").split("\t", 1)[0]
            if not label or label in metadata_labels:
                continue
            gene = _gene_symbol(label)
            observed.add(gene)
            n_gene_rows += 1
    overlap = sorted(query & observed)
    return {
        "n_gene_rows": int(n_gene_rows),
        "n_unique_gene_symbols": int(len(observed)),
        "n_exhaustion_genes_in_scrna": int(len(overlap)),
        "exhaustion_genes_in_scrna": overlap,
    }


def _st_info(path: Path) -> dict[str, object]:
    import anndata as ad

    adata = ad.read_h5ad(path, backed="r")
    try:
        info = {
            "shape": [int(adata.n_obs), int(adata.n_vars)],
            "obs_columns": list(map(str, adata.obs.columns[:20])),
            "var_columns": list(map(str, adata.var.columns[:20])),
            "obsm_keys": list(map(str, adata.obsm.keys())),
            "first_var_names": list(map(str, adata.var_names[:5])),
            "first_obs_names": list(map(str, adata.obs_names[:5])),
        }
    finally:
        adata.file.close()
    return info


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    gene_set_dir = OUT_DIR / "gene_sets"
    gene_set_dir.mkdir(parents=True, exist_ok=True)

    supp_path = RAW_DIR / SUPPLEMENTARY_XLSX
    scrna_path = RAW_DIR / SCRNA_FILE
    st_paths = [RAW_DIR / name for name in ST_FILES]

    for path in [supp_path, scrna_path, *st_paths]:
        if not path.exists():
            raise FileNotFoundError(path)

    gene_sets = _extract_table_s6_gene_sets(supp_path)
    exhaustion_genes = _clean_gene_list(gene_sets[EXHAUSTION_COLUMN])

    pd.DataFrame({"gene": exhaustion_genes}).to_csv(
        gene_set_dir / "tcell_exhaustion_zheng_cell2017.csv",
        index=False,
    )
    (gene_set_dir / "tcell_exhaustion_zheng_cell2017.txt").write_text(
        "\n".join(exhaustion_genes) + "\n",
        encoding="utf-8",
    )

    all_gene_sets_long = []
    for col in gene_sets.columns:
        genes = _clean_gene_list(gene_sets[col])
        for rank, gene in enumerate(genes, start=1):
            all_gene_sets_long.append({"gene_set": col, "rank": rank, "gene": gene})
    pd.DataFrame(all_gene_sets_long).to_csv(gene_set_dir / "table_s6_gene_sets_long.csv", index=False)

    scrna_info = _scrna_header_info(scrna_path)
    scrna_overlap = _scrna_gene_overlap(scrna_path, exhaustion_genes)
    st_info = {path.name: _st_info(path) for path in st_paths}

    # Use h5ad var_names for complete ST overlap checks without loading expression into memory.
    import anndata as ad

    overlap = {}
    for path in st_paths:
        adata = ad.read_h5ad(path, backed="r")
        try:
            var_symbols = {_gene_symbol(name) for name in adata.var_names}
            genes_in_st = sorted(set(exhaustion_genes) & var_symbols)
            overlap[path.name] = {
                "n_unique_gene_symbols": int(len(var_symbols)),
                "n_exhaustion_genes_in_st": int(len(genes_in_st)),
                "exhaustion_genes_in_st": genes_in_st,
            }
        finally:
            adata.file.close()

    manifest = {
        "raw_dir": str(RAW_DIR),
        "outputs": {
            "exhaustion_csv": str(gene_set_dir / "tcell_exhaustion_zheng_cell2017.csv"),
            "exhaustion_txt": str(gene_set_dir / "tcell_exhaustion_zheng_cell2017.txt"),
            "all_table_s6_gene_sets_long": str(gene_set_dir / "table_s6_gene_sets_long.csv"),
        },
        "source_files": {
            "supplementary_tables": str(supp_path),
            "scrna": str(scrna_path),
            "st": [str(path) for path in st_paths],
        },
        "tcell_exhaustion_gene_set": {
            "source": "CytoSPACE Supplementary Table S6; Zheng et al. Cell 2017 Table S4",
            "column": EXHAUSTION_COLUMN,
            "n_genes": len(exhaustion_genes),
            "first_genes": exhaustion_genes[:20],
        },
        "scrna": scrna_info,
        "scrna_exhaustion_gene_overlap": scrna_overlap,
        "st": st_info,
        "st_exhaustion_gene_overlap": overlap,
        "notes": [
            "ST files are the mel1 slide replicates from the Thrane melanoma legacy ST dataset.",
            "The exact replicate used for CytoSPACE Fig.2c still needs to be matched against the paper/code outputs.",
        ],
    }
    (OUT_DIR / "fig2c_resource_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"[OK] exhaustion genes: {len(exhaustion_genes)}")
    print(f"[OK] wrote: {gene_set_dir / 'tcell_exhaustion_zheng_cell2017.csv'}")
    print(f"[OK] wrote: {gene_set_dir / 'tcell_exhaustion_zheng_cell2017.txt'}")
    print(f"[OK] wrote: {gene_set_dir / 'table_s6_gene_sets_long.csv'}")
    print(f"[OK] wrote: {OUT_DIR / 'fig2c_resource_manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

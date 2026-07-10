# BioApp Figure V2.6 First Coordinate Source Report

## Earliest Proven Coordinate Source

The earliest proven source of the frozen BioApp plotting coordinates is:

`data/raw/新建文件夹/BreastCancer_CTA-main/ST_data.RData`

The coordinates were exported by:

`scripts/run_bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run.py`

## Evidence

Phase 2B-R2 parsed `ST_data.RData` as a Seurat object named `bcsa` and selected image slot `BCSA2TumB1`. The exported coordinate table contains `barcode`, `tissue`, `row`, `col`, `imagerow`, `imagecol`, `pxl_col_in_fullres`, and `pxl_row_in_fullres`.

The Phase 2B-R2 script applies:

```text
pxl_col_in_fullres = imagecol / 0.3
pxl_row_in_fullres = imagerow / 0.3
```

This demonstrates that the frozen `imagecol/imagerow` coordinates are processed Seurat image-slot coordinates, not a newly computed Phase 8 or figure-stage coordinate system.

## Classification

`Transformed / registered Visium coordinates from processed Seurat ST_data.RData image slot`

## Space Ranger Bundle Status

No local or downloaded official `spatial/` bundle was proven to be the original source of `ST_data.RData` / `BCSA2TumB1`. Therefore, a raw H&E overlay cannot be recovered directly from Space Ranger files in this phase.

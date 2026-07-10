# BioApp Figure V2.6 Coordinate Recovery Decision

## Final Decision

`PROVENANCE_VERIFIED_BUT_H_AND_E_NOT_DIRECTLY_RECOVERABLE`

## Earliest Coordinate Source Found

`data/raw/新建文件夹/BreastCancer_CTA-main/ST_data.RData`

Generated/exported by:

`scripts/run_bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run.py`

## Coordinate Source Classification

`Transformed / registered Visium coordinates from processed Seurat ST_data.RData image slot`

## Are These Original Space Ranger Image Coordinates?

Not proven. The coordinates are stored in a processed Seurat object image slot (`BCSA2TumB1`) and were exported from `ST_data.RData`. The original Space Ranger `spatial/` directory that produced this processed object was not verified.

## Is A Specific Spatial Bundle Verified?

No. Candidate local bundles can share Visium barcodes and sometimes match array geometry, but barcode/affine/array similarity alone is insufficient because Visium barcodes and array coordinates recur across samples. The downloaded official 10x CytAssist FFPE Human Breast Cancer spatial bundle was already audited in V2.5 and did not match the frozen BioApp barcodes.

## Is H&E Overlay Recoverable?

Not directly in this phase. Recovering a verified H&E overlay would require the original image/scalefactors/tissue_positions source for `ST_data.RData` or an explicit inverse registration transform from the processed `BCSA2TumB1` coordinate system back to a specific H&E image.

## Reason

The coordinate provenance is verified only up to the processed Seurat/ST object. The provenance chain from that object back to an official Space Ranger / Visium / CytAssist spatial bundle is missing.

## Valid Figure Basis

The V2 spot-lattice fallback remains the valid figure basis.

## Biological Results

No biological result, endpoint, metric, threshold, label, or conclusion was changed.

## Recommended Next Step

`BioApp Coordinate Registration Recovery Audit`

# High-resolution profile-mask Fig2D

This directory contains two independently runnable Fig2D experiments.

## all_candidates

Includes all 21 valid dataset/readout/gene-set combinations without outcome-based selection.

```powershell
python scripts/build_highres_profile_mask_fig2d_benchmark.py --project_root .
```

## targeted_validation

Runs the 10 combinations defined in `configs/highres_targeted_validation_fig2d.csv`.

```powershell
python scripts/run_highres_targeted_validation_fig2d.py --project_root .
```

Both scripts use `scripts/highres_enrichment_core.py` for the shared enrichment calculations and plotting.

# Fig. 3B high-resolution profile-mask enrichment

This directory contains the recomputed replacement for manuscript Fig. 3B.

## Run

```bash
python scripts/build_fig3b_highres_profile_mask_enrichment.py --project_root .
```

The fixed representative is `HumanMelanomaPatient2` with a B-cell profile mask. It is the
representative used by the historical Fig. 3B script, recovered from Git history before this
rerun. The current figure is computed from the current expression and mapping outputs; it does
not reuse the historical PDF values.

## Definition

Spatial units are ranked by residual expression of up to 30 masked-target marker genes. The
mapped target-like score is the assigned-cell-count-weighted mean reference marker score. The
low-score set carries exactly 10% effective hit mass per method. Scores tied at the boundary
receive equal fractional hit weights, so the result does not depend on arbitrary ordering within
a score tie. Spatial units otherwise retain their stable source order within equal residual-support
values, preserving the unsimplified running curve used by the visualization.

For the representative case, Peak ES is `0.074045` for CytoSPACE and `0.091479` for SVTuner
(`delta = 0.017434`). This is a representative-case result, not evidence of consistent improvement
across all five datasets: the other four recomputed deltas are negative.

The SVG stores labels as editable text and curves/ticks as vector paths. The metrics CSV contains
all five datasets; the source-values CSV contains every spatial unit and both running curves.
`fig3b_highres_profile_mask_enrichment.svg` is the canonical editable vector deliverable.

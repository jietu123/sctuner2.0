# BioApp Downstream Phase D3 spatial microenvironment context validation

Decision: `PASS`

## Scope

This phase builds the final downstream context layer using frozen spot
coordinates, CTA composition, morphology domains, boundary proxy, and SVTuner
withheld outputs. It does not rerun SVTuner, Stage3, Stage4, or CytoSPACE.

## Niche construction

- Neighborhoods: 1-hop and 2-hop spot graph neighborhoods.
- Clustering: KMeans on local neighborhood descriptors only.
- Chosen k: `6`
- Silhouette score: `0.3557`
- SVTuner withheld variables used for clustering: `false`

## Main outputs

- `panel_i_microenvironment_context.png/pdf`
- `panel_j_spatial_neighborhood_enrichment.png/pdf`
- `spot_microenvironment_table.csv`
- `niche_assignment.csv`
- `niche_summary.csv`
- `neighbor_enrichment.csv`
- `moran_statistics.csv`

## Interpretation boundary

SVTuner withheld outputs show preferential localization within specific
computational local spatial microenvironment contexts. These niches are proxy
niches derived from local neighborhood descriptors. They are not pathology
annotations and do not establish a new biological discovery.

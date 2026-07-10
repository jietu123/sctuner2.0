# BioApp Downstream Phase D0 morphology-domain and boundary freeze

Decision: `PASS`

## Frozen definitions

- Morphology-domain labels: frozen.
- Tumor-stroma boundary proxy: frozen.
- Domain type: CTA composition plus low-resolution H&E patch feature proxy.
- These are not pathology-grade segmentation labels.

## Domain counts

|                      |    0 |
|:---------------------|-----:|
| tumor_core           | 1602 |
| immune_enriched      |  312 |
| stroma_rich          |  147 |
| mixed_boundary       |   94 |
| unmapped_or_excluded |   93 |

## Boundary summary

|                                              | 0                                                                               |
|:---------------------------------------------|:--------------------------------------------------------------------------------|
| boundary_definition                          | stroma_rich/mixed_boundary spots directly adjacent to frozen tumor_core domains |
| median_nearest_neighbor_distance_image_px    | 7.238400000000013                                                               |
| neighbor_radius_image_px                     | 11.94336000000002                                                               |
| interface_band_width_image_px                | 9.048000000000016                                                               |
| n_boundary_seed_spots                        | 217                                                                             |
| n_interface_band_spots                       | 948                                                                             |
| n_tumor_region_spots_for_boundary            | 1602                                                                            |
| n_stroma_or_immune_region_spots_for_boundary | 241                                                                             |

## Guardrails

- SVTuner rerun: false
- Stage3 rerun: false
- Stage4 run: false
- CytoSPACE rerun: false
- Endpoint redefined: false
- Domain definitions use SVTuner outputs: false
- Boundary definitions use SVTuner outputs: false

## Next

Proceed to Panel G morphology-domain concordance and Panel H tumor-stroma interface topology using these frozen D0 definitions.

# BioApp downstream feasibility check

Decision: `FEASIBLE_WITH_CAUTION`

## Available inputs

- embedded_tissue_image: `True`
- spot_coordinates: `True`
- cta_composition: `True`
- endpoint_status: `True`
- baseline_svtuner_scores: `True`

## Spot/domain counts

- n_spots: `2248`
- image_width: `594`
- image_height: `600`
- coords_within_image: `True`
- primary_endpoint_counts: `{'negative': 1754, 'ambiguous': 290, 'positive': 134, 'excluded': 70}`
- tumor_endpoint_counts: `{'positive': 1690, 'negative': 246, 'ambiguous': 242, 'excluded': 70}`
- stroma_endpoint_counts: `{'negative': 1832, 'ambiguous': 250, 'positive': 96, 'excluded': 70}`
- dominant_CTA_counts: `{'Tumor': 1824, 'Immune cells': 217, 'Stroma': 137, 'unmapped': 70}`
- withheld_binary_counts: `{'False': 2099, 'True': 149}`

## H&E patch feature screen

| feature         |   endpoint_pos_vs_neg_d |   withheld_vs_not_d |   spearman_withheld_score |   spearman_immune_fraction |   spearman_tumor_fraction |   spearman_stroma_fraction |
|:----------------|------------------------:|--------------------:|--------------------------:|---------------------------:|--------------------------:|---------------------------:|
| rgb_mean        |                1.41973  |            1.13624  |                 0.223546  |                  0.371104  |                 -0.57378  |                  0.386688  |
| rgb_std         |               -0.302369 |           -0.315466 |                -0.0975983 |                 -0.0755563 |                  0.114998 |                 -0.0675923 |
| darkness        |               -1.41973  |           -1.13624  |                -0.223542  |                 -0.371107  |                  0.573781 |                 -0.386687  |
| saturation      |               -1.34735  |           -1.05788  |                -0.173519  |                 -0.359894  |                  0.527424 |                 -0.313689  |
| red_blue        |                0.481483 |            0.605122 |                 0.276747  |                  0.0404768 |                 -0.11029  |                  0.121344  |
| edge_proxy      |               -0.595493 |           -0.377975 |                -0.227374  |                 -0.0661031 |                  0.142676 |                 -0.101582  |
| tissue_coverage |              nan        |          nan        |               nan         |                nan         |                nan        |                nan         |

## Tumor-stroma interface proxy

- boundary_proxy_available: `True`
- n_boundary_proxy_spots: `1374`
- median_abs_distance_endpoint_positive: `7.249628073218654`
- median_abs_distance_withheld: `7.279999999999973`
- median_abs_distance_not_withheld: `0.0`
- spearman_abs_distance_withheld_score: `0.41567621797285176`

## Recommendation

Proceed first with a lightweight downstream Phase D0/D1 design: freeze H&E/CTA morphology-domain labels, then implement Panel G morphology-domain concordance and Panel H tumor-stroma interface topology. Do not yet run all four downstream experiments.
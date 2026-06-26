from pathlib import Path

import pandas as pd

from src.stages.stage4_cytospace import (
    load_stage3b_blank_spots,
    restore_blank_spot_rows,
)


def test_load_stage3b_blank_spots_uses_region_decision(tmp_path: Path):
    sample = "sample"
    scores_dir = (
        tmp_path
        / "data"
        / "processed"
        / sample
        / "stage3b_st_unsupported"
    )
    scores_dir.mkdir(parents=True)
    scores = pd.DataFrame(
        {
            "is_spot_candidate": [True, True, False],
            "is_unsupported_region": [True, False, False],
        },
        index=["s1", "s2", "s3"],
    )
    scores.to_csv(scores_dir / "spot_unsupported_scores.csv")

    blank, path = load_stage3b_blank_spots(tmp_path, sample, {})

    assert path == (scores_dir / "spot_unsupported_scores.csv").resolve()
    assert blank.index.tolist() == ["s1"]


def test_restore_blank_rows_adds_zero_rows_in_full_spot_order(tmp_path: Path):
    path = tmp_path / "cell_type_assignments_by_spot.csv"
    pd.DataFrame(
        {
            "spot_id": ["s1", "s3"],
            "A": [2, 0],
            "B": [0, 2],
            "Total cells": [2, 2],
        }
    ).to_csv(path, index=False)

    audit = restore_blank_spot_rows(
        path,
        ["s1", "s2", "s3"],
        {"s2"},
    )
    restored = pd.read_csv(path, index_col=0)

    assert restored.index.tolist() == ["s1", "s2", "s3"]
    assert restored.loc["s2"].sum() == 0
    assert audit == {
        "full_spots": 3,
        "blank_rows": 1,
        "blank_nonzero_rows": 0,
    }

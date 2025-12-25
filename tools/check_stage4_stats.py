import pandas as pd
import hashlib
import json
from pathlib import Path


def sha1(p: Path) -> str:
    h = hashlib.sha1()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def run_once(sample: str):
    project_root = Path(__file__).resolve().parents[1]
    root = project_root / "result" / sample / "stage4_cytospace" / "cytospace_output"
    ca = pd.read_csv(root / "cell_assignment.csv")
    al = pd.read_csv(root / "assigned_locations.csv")
    by = pd.read_csv(root / "cell_type_assignments_by_spot.csv", index_col=0)
    relabel = pd.read_csv(project_root / "data" / "processed" / sample / "stage3_typematch" / "cell_type_relabel.csv")

    print("==", sample)
    print("cell_assignment rows", len(ca))
    print("assigned_locations rows", len(al))
    print("by_spot total cells", int(by["Total cells"].sum()))

    a_cd8 = (ca["cell_type"] == "T cells CD8").sum()
    bs_cd8 = by["T cells CD8"].sum() if "T cells CD8" in by.columns else 0
    print("CD8 in assignment", int(a_cd8))
    print("CD8 total in by_spot", int(bs_cd8))
    print("columns sample", list(by.columns[:10]))
    print("cell_type unique sample", ca["cell_type"].unique()[:10])
    print("orig_type counts (top)", relabel["orig_type"].value_counts().head())
    print("orig_type CD8 rows", len(relabel[relabel["orig_type"] == "T cells CD8"]))
    ct_df = pd.read_csv(root.parent / "cytospace_input" / "cell_types_for_cytospace.csv")
    print("input cell_type counts (top)", ct_df["cell_type"].value_counts().head())
    print("input CD8 count", ct_df[ct_df["cell_type"] == "T cells CD8"].shape[0])
    print("example orig_type CD8 rows:", relabel[relabel["orig_type"] == "T cells CD8"].head())

    sha = {
        "sim_info": sha1(project_root / "data/processed/real_brca_simS0_seed42/stage1_preprocess/exported/sim_info.json"),
        "sc_metadata": sha1(project_root / "data/processed/real_brca_simS0_seed42/stage1_preprocess/exported/sc_metadata.csv"),
        "cell_assignment": sha1(root / "cell_assignment.csv"),
    }
    print("sha1", json.dumps(sha, indent=2))
    print()


def main():
    run_once("real_brca_simS0_seed42")


if __name__ == "__main__":
    main()
import pandas as pd
import hashlib
import json
from pathlib import Path


def sha1(p: Path) -> str:
    h = hashlib.sha1()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    project_root = Path(__file__).resolve().parents[1]
    root = project_root / "result" / "real_brca_simS0_seed42" / "stage4_cytospace" / "cytospace_output"
    ca = pd.read_csv(root / "cell_assignment.csv")
    al = pd.read_csv(root / "assigned_locations.csv")
    by = pd.read_csv(root / "cell_type_assignments_by_spot.csv", index_col=0)

    print("cell_assignment rows", len(ca))
    print("assigned_locations rows", len(al))
    print("by_spot total cells", int(by["Total cells"].sum()))

    a_cd8 = (ca["cell_type"] == "T cells CD8").sum()
    bs_cd8 = by["T cells CD8"].sum() if "T cells CD8" in by.columns else 0
    print("CD8 in assignment", int(a_cd8))
    print("CD8 total in by_spot", int(bs_cd8))

    sha = {
        "sim_info": sha1(project_root / "data/processed/real_brca_simS0_seed42/stage1_preprocess/exported/sim_info.json"),
        "sc_metadata": sha1(project_root / "data/processed/real_brca_simS0_seed42/stage1_preprocess/exported/sc_metadata.csv"),
        "cell_assignment": sha1(root / "cell_assignment.csv"),
    }
    print("sha1", json.dumps(sha, indent=2))


if __name__ == "__main__":
    main()


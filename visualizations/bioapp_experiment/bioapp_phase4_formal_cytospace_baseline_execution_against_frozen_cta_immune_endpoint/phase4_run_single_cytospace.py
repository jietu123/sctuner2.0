import argparse
import json
import sys
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
args = parser.parse_args()
cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))

sys.path.insert(0, cfg["cytospace_package_path"])
# datatable.fread can mis-handle non-ASCII Windows paths in this workspace.
# Force the official CytoSPACE reader to use its pandas fallback.
import cytospace.common.common as cytospace_common
cytospace_common.dt = None
from cytospace.cytospace import main_cytospace

main_cytospace(
    scRNA_path=cfg["scRNA_path"],
    cell_type_path=cfg["cell_type_path"],
    n_cells_per_spot_path=None,
    st_cell_type_path=None,
    cell_type_fraction_estimation_path=cfg["cell_type_fraction_estimation_path"],
    spaceranger_path=None,
    st_path=cfg["st_path"],
    coordinates_path=cfg["coordinates_path"],
    output_folder=cfg["output_folder"],
    output_prefix="",
    mean_cell_numbers=cfg.get("mean_cell_numbers", 5),
    downsample_off=True,
    scRNA_max_transcripts_per_cell=1500,
    solver_method=cfg.get("solver_method", "lap_CSPR"),
    distance_metric="Pearson_correlation",
    sampling_method="duplicates",
    single_cell=False,
    number_of_selected_spots=10000,
    sampling_sub_spots=False,
    number_of_selected_sub_spots=10000,
    number_of_processors=cfg.get("number_of_processors", 1),
    seed=cfg.get("seed", 20260705),
    plot_off=True,
    geometry="honeycomb",
    max_num_cells_plot=50000,
    num_column=3,
)

"""
Centralized defaults for SVTuner.
If you need to change default R command, paths, or stage selections, do it here.
"""
from pathlib import Path

# Default R command (Windows conda env recommended to avoid DLL issues)
DEFAULT_R_CMD = "conda run -n cytospace_v1.1.0_py310 Rscript"

# Default stage selection (can be overridden via CLI)
# Currently default only Stage1; Stage0 is optional dummy.
DEFAULT_STAGES = "1"

# Project paths (resolved at runtime)
def detect_project_root() -> Path:
    """Return project root assuming this file is under src/."""
    here = Path(__file__).resolve()
    return here.parents[1]


def project_path(*parts: str) -> Path:
    return detect_project_root().joinpath(*parts)


# Common subdirectories (use project_path(...) to resolve to absolute)
DATA_RAW = Path("data/raw")
DATA_PROCESSED = Path("data/processed")
RESULT = Path("result")
LOGS = Path("logs")


__all__ = [
    "DEFAULT_R_CMD",
    "DEFAULT_STAGES",
    "detect_project_root",
    "project_path",
    "DATA_RAW",
    "DATA_PROCESSED",
    "RESULT",
    "LOGS",
]

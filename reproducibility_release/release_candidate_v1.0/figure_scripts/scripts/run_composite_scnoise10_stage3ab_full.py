#!/usr/bin/env python
"""Compatibility entry point for the 10% composite Stage3A + Stage3B benchmark."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.run_composite_noise_stage3ab_full import main


if __name__ == "__main__":
    raise SystemExit(main(default_noise_fraction=0.10))

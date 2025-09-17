#!/usr/bin/env python3
"""CR2SUB pipeline launcher."""

from __future__ import annotations

import sys
from typing import List, Tuple

from scripts.quicksetup import run_full_setup


PIPELINE_STEPS: List[Tuple[str, List[str]]] = [
    ("Consolidate metadata and time series", [
        "Rscript",
        "scripts/01_consolidate_cr2sub_v1.R",
    ]),
    ("Remove groundwater outliers", [
        sys.executable,
        "-m",
        "jupyter",
        "nbconvert",
        "--to",
        "notebook",
        "--execute",
        "scripts/02_remove_outliers_cr2sub_v1.ipynb",
    ]),
    ("Derive well attributes", [
        "Rscript",
        "scripts/03_process_attributes_cr2sub_v1.R",
    ]),
]


def main() -> int:
    return run_full_setup(PIPELINE_STEPS)


if __name__ == "__main__":
    sys.exit(main())

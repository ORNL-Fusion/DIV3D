#!/usr/bin/env python3

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parents[0]))

from compare_outputs_common import run_compare_outputs


if __name__ == "__main__":
    run_compare_outputs(expected_timing_rows=50)

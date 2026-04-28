#!/usr/bin/env python3

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parents[0]))

from plot_results_common import run_plot_results


def main():
    run_plot_results(
        SCRIPT_DIR,
        description="Plot vessel-hit intersections for the axisymmetric vessel test.",
        title_2d="TOK/CAT vessel intersections",
        title_3d="TOK/CAT vessel intersections 3D",
        plot_parts=False,
        plot_3d_help="Plot vessel and intersections in Cartesian 3D.",
    )


if __name__ == "__main__":
    main()

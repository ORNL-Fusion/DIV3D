#!/usr/bin/env python3

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parents[0]))

from plot_results_common import run_plot_results


def main():
    run_plot_results(
        SCRIPT_DIR,
        description="Plot intersections, vessel, and parts for the axisymmetric vessel/parts test.",
        title_2d="TOK/CAT vessel and part intersections",
        title_3d="TOK/CAT vessel and part intersections 3D",
        plot_parts=True,
    )


if __name__ == "__main__":
    main()

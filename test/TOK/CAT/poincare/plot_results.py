#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[3]
sys.path.insert(0, str(REPO_ROOT / "postprocessing"))

from poincare_routines import collect_surfaces, find_surface_files, plot_slices


def main():
    parser = argparse.ArgumentParser(
        description="Plot Poincare surface_data.out files for a quick visual check."
    )
    parser.add_argument(
        "--ref",
        action="store_true",
        help="Plot reference files from ref/ instead of current run outputs.",
    )
    parser.add_argument(
        "--slice",
        type=int,
        default=None,
        help="Plot one slice index. By default all slices are plotted.",
    )
    parser.add_argument(
        "--save",
        metavar="PNG",
        help="Save the plot to a PNG file.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open an interactive plotting window.",
    )
    parser.add_argument(
        "--no-combined",
        action="store_true",
        help="Do not add the combined all-slices panel.",
    )
    parser.add_argument(
        "--parts",
        action="store_true",
        help="Overlay .part contours found in ../geo.",
    )
    parser.add_argument(
        "--parts-dir",
        default="../geo",
        help="Directory to search for .part files when --parts is used.",
    )
    args = parser.parse_args()

    data_dir = SCRIPT_DIR / "ref" if args.ref else SCRIPT_DIR
    files = find_surface_files(data_dir)

    if not files:
        source = "ref/" if args.ref else "the current case directory"
        raise SystemExit(
            f"No surface_data.out.#### files found in {source}. "
            "Run ./run_test first or pass --ref."
        )

    slices = collect_surfaces(files)
    if args.slice is not None and (args.slice < 0 or args.slice >= len(slices)):
        raise SystemExit(f"Slice index must be between 0 and {len(slices) - 1}.")

    part_files = None
    if args.parts:
        parts_dir = (SCRIPT_DIR / args.parts_dir).resolve()
        part_files = sorted(parts_dir.glob("*.part"))
        if not part_files:
            raise SystemExit(f"No .part files found in {parts_dir}.")

    plot_slices(
        slices,
        args.slice,
        Path(args.save) if args.save else None,
        show_plot=not args.no_show,
        include_combined=not args.no_combined,
        part_files=part_files,
    )


if __name__ == "__main__":
    main()

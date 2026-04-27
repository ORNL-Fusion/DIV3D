#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[3]
sys.path.insert(0, str(REPO_ROOT / "postprocessing"))

from div3d_routines import find_part_files, plot_part_files, read_part_files


def main():
    parser = argparse.ArgumentParser(
        description="Plot the DIV3D .part files in this geometry directory."
    )
    parser.add_argument(
        "parts",
        nargs="*",
        help="Part files to plot. Defaults to all *.part files in this directory.",
    )
    parser.add_argument(
        "--slice",
        type=int,
        default=0,
        help="Toroidal slice index to plot. Default: 0.",
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
    args = parser.parse_args()

    if args.save or args.no_show:
        import matplotlib

        matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    part_files = [Path(part) for part in args.parts] if args.parts else find_part_files(SCRIPT_DIR)
    if not part_files:
        raise SystemExit("No .part files found.")

    parts = read_part_files(part_files)
    fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
    plot_part_files(parts, ax=ax, slice_index=args.slice)
    ax.set_title("TOK/CAT part contours")

    if args.save:
        fig.savefig(args.save, dpi=200)
        print(f"Wrote {args.save}")

    if not args.no_show:
        plt.show()

    plt.close(fig)


if __name__ == "__main__":
    main()

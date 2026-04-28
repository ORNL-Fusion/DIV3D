#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path


def run_plot_results(
    script_dir,
    description,
    title_2d,
    title_3d,
    plot_parts=False,
    plot_3d_help="Plot vessel, parts, and intersections in Cartesian 3D.",
):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--ref",
        action="store_true",
        help="Plot reference files from ref/ instead of current run outputs.",
    )
    parser.add_argument(
        "--save",
        metavar="PNG",
        help="Save the plot to a PNG file.",
    )
    parser.add_argument(
        "--3D",
        dest="plot_3d",
        action="store_true",
        help=plot_3d_help,
    )
    parser.add_argument(
        "--nphi_AS_to_3D",
        dest="nphi_AS_to_3D",
        type=int,
        default=16,
        help="Number of toroidal cuts for mapping AS parts into 3D. Default: 16.",
    )
    parser.add_argument(
        "--surface_line",
        action="store_true",
        help="Plot surface_line.out as R-Z dots on 2D plots.",
    )
    parser.add_argument(
        "--launch_pts",
        action="store_true",
        help="Plot launch_pts.out as R-Z dots on 2D plots.",
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

    repo_root = script_dir.parents[3]
    sys.path.insert(0, str(repo_root / "postprocessing"))

    from div3d_routines import (
        plot_hits_3d,
        plot_hits_axisymmetric,
        plot_launch_pts_rz,
        plot_surface_line_rz,
        read_launch_pts_file,
        read_surface_line_file,
    )

    data_dir = script_dir / "ref" if args.ref else script_dir
    int_pts_file = data_dir / "int_pts.out"
    surface_line_file = data_dir / "surface_line.out"
    launch_pts_file = data_dir / "launch_pts.out"

    if not int_pts_file.exists():
        source = "ref/" if args.ref else "the current case directory"
        raise SystemExit(
            f"No int_pts.out found in {source}. Run ./run_test first or pass --ref."
        )

    if args.plot_3d:
        fig = plt.figure(figsize=(9, 8), constrained_layout=True)
        ax = fig.add_subplot(111, projection="3d")
        plot_hits_3d(
            run_settings_file=script_dir / "run_settings.nml",
            int_pts_file=int_pts_file,
            ax=ax,
            plot_parts=plot_parts,
            nphi_AS_to_3D=args.nphi_AS_to_3D,
        )
        ax.set_title(title_3d)
    else:
        fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
        plot_hits_axisymmetric(
            run_settings_file=script_dir / "run_settings.nml",
            int_pts_file=int_pts_file,
            ax=ax,
            plot_parts=plot_parts,
        )
        if args.surface_line:
            if not surface_line_file.exists():
                source = "ref/" if args.ref else "the current case directory"
                raise SystemExit(f"No surface_line.out found in {source}.")
            plot_surface_line_rz(read_surface_line_file(surface_line_file), ax=ax)
            ax.legend()
        if args.launch_pts:
            if not launch_pts_file.exists():
                source = "ref/" if args.ref else "the current case directory"
                raise SystemExit(f"No launch_pts.out found in {source}.")
            plot_launch_pts_rz(read_launch_pts_file(launch_pts_file), ax=ax)
            ax.legend()
        ax.set_title(title_2d)

    if args.save:
        fig.savefig(args.save, dpi=200)
        print(f"Wrote {args.save}")

    if not args.no_show:
        plt.show()

    plt.close(fig)

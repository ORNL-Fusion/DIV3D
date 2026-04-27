import glob
import os
import tempfile
from pathlib import Path

import numpy as np


def read_surface_data(filename):
    """Read a Poincare surface_data.out file."""
    with open(filename, "r") as f:
        lines = f.readlines()

    num_surfs, num_slices = map(int, lines[0].split())
    data = {"slices": {}}
    index = 1

    for slice_idx in range(num_slices):
        dphi, num_points = lines[index].split()
        dphi = float(dphi)
        num_points = int(num_points)
        index += 1

        data["slices"][slice_idx] = {
            "dphi": dphi,
            "num_points": num_points,
            "R": np.zeros((num_surfs, num_points)),
            "Z": np.zeros((num_surfs, num_points)),
        }

        for surf_idx in range(num_surfs):
            data["slices"][slice_idx]["R"][surf_idx, :] = np.array(
                lines[index].split(), dtype=float
            )
            index += 1
            data["slices"][slice_idx]["Z"][surf_idx, :] = np.array(
                lines[index].split(), dtype=float
            )
            index += 1

    return data


def find_surface_files(data_dir="."):
    """Return sorted rank-local Poincare surface files in a directory."""
    return sorted(Path(data_dir).glob("surface_data.out.[0-9][0-9][0-9][0-9]"))


def collect_surfaces(files):
    """Combine rank-local surface files into a list indexed by toroidal slice."""
    slices = []

    for fname in files:
        rank_data = read_surface_data(fname)
        if not slices:
            slices = [
                {"dphi": item["dphi"], "surfaces": []}
                for item in rank_data["slices"].values()
            ]

        for slice_idx, item in rank_data["slices"].items():
            for surf_idx in range(item["R"].shape[0]):
                slices[slice_idx]["surfaces"].append(
                    {
                        "rank_file": Path(fname).name,
                        "R": item["R"][surf_idx, :],
                        "Z": item["Z"][surf_idx, :],
                    }
                )

    return slices


def _prepare_matplotlib(output_file=None, show_plot=True):
    if "MPLCONFIGDIR" not in os.environ:
        mpl_config_dir = Path(tempfile.gettempdir()) / "div3d_matplotlib"
        mpl_config_dir.mkdir(exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(mpl_config_dir)

    if output_file or not show_plot:
        import matplotlib

        matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    return plt


def _plot_slice_data(ax, slice_data):
    for surf in slice_data["surfaces"]:
        valid = (
            np.isfinite(surf["R"])
            & np.isfinite(surf["Z"])
            & (surf["R"] > 0.0)
            & (surf["R"] < 10.0)
            & (np.abs(surf["Z"]) < 10.0)
        )
        ax.plot(
            surf["R"][valid],
            surf["Z"][valid],
            marker=".",
            linestyle="None",
            markersize=2,
            alpha=0.75,
        )


def _overlay_parts(ax, parts, phi_deg, label_parts=True):
    from div3d_routines import plot_parts_at_phi

    plot_parts_at_phi(
        parts,
        np.deg2rad(phi_deg),
        ax=ax,
        label_parts=label_parts,
        color="black",
        alpha=0.65,
    )


def plot_slices(
    slices,
    slice_index=None,
    output_file=None,
    show_plot=True,
    include_combined=False,
    part_files=None,
):
    """Plot one or all toroidal slices from collected Poincare surface data."""
    plt = _prepare_matplotlib(output_file, show_plot)

    parts = []
    if part_files:
        from div3d_routines import read_part_files

        parts = read_part_files(part_files)

    if slice_index is None:
        selected = list(range(len(slices)))
    else:
        selected = [slice_index]

    num_panels = len(selected)
    if include_combined and len(selected) > 1:
        num_panels += 1

    fig, axes = plt.subplots(
        1,
        num_panels,
        figsize=(6 * num_panels, 6),
        squeeze=False,
        constrained_layout=True,
    )

    for ax, islice in zip(axes[0], selected):
        slice_data = slices[islice]
        _plot_slice_data(ax, slice_data)
        if parts:
            _overlay_parts(ax, parts, slice_data["dphi"])
        ax.set_title(f"Phi = {slice_data['dphi']:.3f} deg")
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, alpha=0.3)
        if parts:
            ax.legend(fontsize=8)

    if include_combined and len(selected) > 1:
        ax = axes[0][len(selected)]
        for islice in selected:
            _plot_slice_data(ax, slices[islice])
        if parts:
            for islice in selected:
                _overlay_parts(ax, parts, slices[islice]["dphi"], label_parts=(islice == selected[0]))
        ax.set_title("All slices")
        ax.set_xlabel("R [m]")
        ax.set_ylabel("Z [m]")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, alpha=0.3)
        if parts:
            ax.legend(fontsize=8)

    if output_file:
        fig.savefig(output_file, dpi=200)
        print(f"Wrote {output_file}")

    if show_plot:
        plt.show()

    plt.close(fig)


def plot_all_surfaces_one_plot(slice_index=0, data_dir="."):
    """Plot the chosen slice from all rank-local surface files."""
    files = find_surface_files(data_dir)
    if not files:
        print("No 'surface_data.out.####' files found in current directory.")
        return

    slices = collect_surfaces(files)
    if slice_index < 0 or slice_index >= len(slices):
        print(f"Slice index must be between 0 and {len(slices) - 1}.")
        return

    plot_slices(slices, slice_index=slice_index)


if __name__ == "__main__":
    files = sorted(glob.glob("surface_data.out.[0-9][0-9][0-9][0-9]"))
    plot_slices(collect_surfaces(files), slice_index=0)

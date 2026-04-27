import pandas as pd
import os
import tempfile
from pathlib import Path

if "MPLCONFIGDIR" not in os.environ:
    mpl_config_dir = Path(tempfile.gettempdir()) / "div3d_matplotlib"
    mpl_config_dir.mkdir(exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(mpl_config_dir)

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def main():
    filename = 'hitline.out'
    hitline_data = read_hitline(filename)
    num_hitlines_plot = 10
#    plot_hitlines_3d(hitline_data, num_hitlines_plot)

    filename = 'surface_line.out'
    surf_line = read_surface_line_file(filename)

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    plot_hitlines_3d(hitline_data, num_hitlines_plot, ax=ax)  # Pass ax
    plot_surface_line(surf_line, ax=ax,num_points=2000)  # Add surface line to the same plot
    ax.legend()
    plt.show()

    print(surf_line["Z"][0:10])
    

def read_surface_line_file(filename):
    """
    Reads a surface line file written by the corresponding Fortran routine.

    Parameters:
        filename (str): Path to the surface line file.

    Returns:
        dict: Dictionary with metadata and arrays for R, Z, and phi.
    """
    with open(filename, 'r') as file:
        # Read the first line: period, nip0, ip_step
        period, nip0, ip_step = map(float, file.readline().split())

        # Read the second line: nsteps_line + 1 (number of data points)
        nline = int(file.readline().strip())

        # Read the R, Z, and phi values per point
        r_vals, z_vals, phi_vals = [], [], []
        for _ in range(nline):
            r_vals.append(float(file.readline().strip()))
            z_vals.append(float(file.readline().strip()))
            phi_vals.append(float(file.readline().strip()))

    return {
        "metadata": {"period": period, "nip0": int(nip0), "ip_step": int(ip_step), "nline": nline},
        "R": np.array(r_vals),
        "Z": np.array(z_vals),
        "phi": np.array(phi_vals),
    }

    
def plot_hitlines_3d(hitline_data, num_hitlines_plot, ax=None):
    """
    Plots the first num_hitlines_plot hitlines in 3D using Cartesian coordinates (X, Y, Z).
    
    Parameters:
        hitline_data (list): List of dictionaries containing 'r', 'z', and 'phi' arrays.
        num_hitlines_plot (int): Number of hitlines to plot.
        ax (matplotlib.axes._subplots.Axes3DSubplot, optional): Existing 3D axis to plot on.
    """
    num_hitlines_plot = min(num_hitlines_plot, len(hitline_data))  # Ensure we don't exceed available hitlines

    # Create a new figure if no existing axis is provided
    if ax is None:
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')

    # Plot each hitline
    for i in range(num_hitlines_plot):
        R_hit = hitline_data[i]["R"]
        Z_hit = hitline_data[i]["Z"]
        Phi_hit = hitline_data[i]["phi"]

        # Convert cylindrical (R, Phi, Z) to Cartesian (X, Y, Z)
        X_hit = R_hit * np.cos(Phi_hit)
        Y_hit = R_hit * np.sin(Phi_hit)

        ax.plot(X_hit, Y_hit, Z_hit, label=f"Hit Line {i+1}")

    # Set labels and title only if a new figure was created
    #    if ax is None:
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_zlabel("Z (m)")
    ax.set_title("3D Plot of Hitlines in Cartesian Coordinates")
#    ax.legend()
#    plt.show()

def plot_surface_line(surf_line, ax=None, num_points=1000):
    """
    Adds the surface line data to an existing 3D plot or creates a new one if needed.

    Parameters:
        surf_line (dict): Dictionary containing 'R', 'Z', and 'phi' arrays.
        ax (matplotlib.axes._subplots.Axes3DSubplot, optional): Existing 3D axis to plot on.
        num_points (int, optional): Number of points to plot. Defaults to 1000.
    """
    # Subsample the data to reduce plotting load
    total_points = len(surf_line["R"])
    step = max(1, total_points // num_points)  # Ensure at least 1 step size
    
    R_sampled = surf_line["R"][::step]
    Z_sampled = surf_line["Z"][::step]
    Phi_sampled = surf_line["phi"][::step]
    
    # Convert cylindrical (R, Phi, Z) to Cartesian (X, Y, Z)
    X_surf = R_sampled * np.cos(Phi_sampled)
    Y_surf = R_sampled * np.sin(Phi_sampled)
    Z_surf = Z_sampled

    # Create new figure if no axis is provided
    if ax is None:
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')

    # Plot the surface line
    ax.plot(X_surf, Y_surf, Z_surf, label="Surface Line", color='magenta',marker='.', linestyle='none')

    # Set labels only if a new figure was created
    #    if ax is None:
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_zlabel("Z (m)")
    ax.set_title("3D Plot of Hitlines and Surface Line")

#    plt.show()

def temp():    
    # Read int_pts.out
    file_path = 'int_pts.out'
    int_pts = read_int_pts_file(file_path)
    print(int_pts.head(10))

    # Read int_pts
    filename = "run_settings.nml"
    namelist_data = read_run_settings(filename)

    fname_ves = namelist_data["run_settings"]["fname_ves"]

    filename = 'hitline.out'
    hitline_data = read_hitline(filename)
    

    print(fname_ves)

    ves_part = read_part_file(fname_ves)
    label = ves_part["metadata"]["label"]
    print(label)

    """
    Plots the first slice of R,Z data from ves_part as lines and overlays R,Z points from int_pts as x's.
    """
    num_hitlines_plot = 10
    
    plt.figure(figsize=(8, 6))
    plt.plot(ves_part["coordinates"]["R"], ves_part["coordinates"]["Z"], label="Vessel Part", linestyle='-', color='blue')
    plt.scatter(int_pts["R"], int_pts["Z"], label="Int Points", color='red', marker='x')
    num_hitlines_plot = min(num_hitlines_plot, len(hitline_data))  # Ensure we don't exceed available hitlines
    for i in range(num_hitlines_plot):
        plt.plot(hitline_data[i]["R"], hitline_data[i]["Z"], color='black', linestyle='-')

    plt.xlabel("R (m)")
    plt.ylabel("Z (m)")
    plt.title("First Slice of Vessel Part and Intersection Points")
    plt.legend()
    plt.grid()
    plt.show()

def read_hitline(filename):
    """
    Reads a hitline file and returns a list of dictionaries, each containing r, z, and phi arrays.

    Parameters:
        filename (str): Path to the hitline file.

    Returns:
        list: A list of dictionaries, each with 'r', 'z', and 'phi' keys containing numpy arrays.
    """
    import numpy as np

    hitlines = []

    with open(filename, 'r') as file:
        lines = file.readlines()

    i = 0
    while i < len(lines):
        line_info = [int(val) for val in lines[i].split()]
        line_index = line_info[0]
        npts = line_info[1]
        i += 1  # Move to the data lines

        r = np.array([float(val) for val in lines[i].split()])
        i += 1
        z = np.array([float(val) for val in lines[i].split()])
        i += 1
        phi = np.array([float(val) for val in lines[i].split()])
        i += 1

        hitlines.append({"line_index": line_index, "R": r, "Z": z, "phi": phi})

    return hitlines




def read_part_file(filename):
    """
    Reads a .part file and returns a dictionary containing the metadata and a DataFrame for coordinates.
    
    Parameters:
        filename (str): Path to the .part file.
    
    Returns:
        dict: Dictionary with metadata and a DataFrame with R, Z coordinates.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    metadata = {
        "label": lines[0].strip(),
        "ntor": int(lines[1].split()[0]),
        "npol": int(lines[1].split()[1]),
        "nfp": int(lines[1].split()[2]),
        "rshift": float(lines[1].split()[3]),
        "zshift": float(lines[1].split()[4]),
        "force_non_AS": False  # Default value
    }
    
    # Handle optional force_non_AS field
    parts = lines[1].split()
    if len(parts) > 5:
        metadata["force_non_AS"] = bool(int(parts[5]))
    
    phi_value = float(lines[2].strip())
    
    # Read R, Z coordinates
    data = []
    for line in lines[3:]:
        parts = line.strip().split(',')
        if len(parts) == 2:
            data.append([float(parts[0]), float(parts[1])])
    
    df = pd.DataFrame(data, columns=["R", "Z"])/100
    
    return {"metadata": metadata, "phi": phi_value, "coordinates": df}


def read_part_slices(filename):
    """
    Reads a DIV3D .part file and returns metadata plus one coordinate table per
    toroidal slice. Coordinates are converted from cm to m.
    """
    with open(filename, 'r') as file:
        lines = [line.strip() for line in file.readlines() if line.strip()]

    parts = lines[1].split()
    metadata = {
        "filename": str(filename),
        "label": lines[0],
        "ntor": int(parts[0]),
        "npol": int(parts[1]),
        "nfp": int(parts[2]),
        "rshift": float(parts[3]),
        "zshift": float(parts[4]),
        "force_non_AS": False,
    }

    if len(parts) > 5:
        metadata["force_non_AS"] = bool(int(parts[5]))

    slices = []
    line_index = 2
    for _ in range(metadata["ntor"]):
        phi = float(lines[line_index])
        line_index += 1

        coords = []
        for _ in range(metadata["npol"]):
            r_cm, z_cm = [float(value) for value in lines[line_index].split(",")]
            coords.append([r_cm / 100.0, z_cm / 100.0])
            line_index += 1

        slices.append(
            {
                "phi": phi,
                "phi_rad": np.deg2rad(phi),
                "coordinates": pd.DataFrame(coords, columns=["R", "Z"]),
            }
        )

    metadata["is_axisymmetric"] = part_is_axisymmetric({"metadata": metadata, "slices": slices})

    return {"metadata": metadata, "slices": slices}


def read_part_files(filenames):
    """Read a list of DIV3D .part files."""
    return [read_part_slices(filename) for filename in filenames]


def part_is_axisymmetric(part, tol=1.e-8):
    """
    Return true if all toroidal slices have the same R-Z contour.

    This mirrors the DIV3D check used for part/vessel axisymmetry: compare each
    toroidal slice against the first slice in R and Z.
    """
    slices = part["slices"]
    if len(slices) < 2:
        return True

    reference = slices[0]["coordinates"][["R", "Z"]].to_numpy()
    for slice_data in slices[1:]:
        candidate = slice_data["coordinates"][["R", "Z"]].to_numpy()
        if candidate.shape != reference.shape:
            return False
        check_as = (
            np.max(np.abs(candidate[:, 0] - reference[:, 0]))
            + np.max(np.abs(candidate[:, 1] - reference[:, 1]))
        )
        if check_as >= tol:
            return False

    return True


def wrap_phi(phi, period=2.0 * np.pi):
    """Wrap a toroidal angle in radians into [0, period)."""
    return phi % period


def get_part_slice_at_phi(part, phi_rad, period=2.0 * np.pi):
    """
    Select the R-Z contour for a part at a toroidal location.

    Axisymmetric parts use the first toroidal contour, matching DIV3D's
    intersection logic. DIV3D builds triangles for non-axisymmetric parts; this
    2D plot helper does not yet mirror that triangle path.
    """
    if part["metadata"].get("is_axisymmetric", False):
        return part["slices"][0]

    return None


def plot_parts_at_phi(
    parts,
    phi_rad,
    ax=None,
    period=2.0 * np.pi,
    label_parts=True,
    linestyle="-",
    color="black",
    alpha=0.7,
):
    """Overlay DIV3D part contours selected at a toroidal angle."""
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))

    for part in parts:
        part_slice = get_part_slice_at_phi(part, phi_rad, period=period)
        if part_slice is None:
            print(
                f"Skipping non-axisymmetric part '{part['metadata']['label']}' "
                "in R-Z overlay; DIV3D uses triangles for this case."
            )
            continue

        coords = part_slice["coordinates"]
        label = part["metadata"]["label"]
        plot_label = f"{label} (AS)"

        ax.plot(
            coords["R"],
            coords["Z"],
            linestyle=linestyle,
            color=color,
            alpha=alpha,
            linewidth=1.2,
            label=plot_label,
        )

        if label_parts:
            ax.text(
                coords["R"].mean(),
                coords["Z"].mean(),
                label,
                ha="center",
                va="center",
                fontsize=8,
                color=color,
            )

    return ax


def plot_part_files(parts, ax=None, slice_index=0, label_parts=True):
    """
    Plot one toroidal slice from each DIV3D .part file in an R-Z view.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))

    for part in parts:
        metadata = part["metadata"]
        if slice_index >= len(part["slices"]):
            continue

        coords = part["slices"][slice_index]["coordinates"]
        label = metadata["label"]
        ax.plot(coords["R"], coords["Z"], marker=".", linestyle="-", label=label)

        if label_parts:
            ax.text(
                coords["R"].mean(),
                coords["Z"].mean(),
                label,
                ha="center",
                va="center",
            )

    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.3)
    ax.legend()

    return ax


def find_part_files(data_dir="."):
    """Return sorted .part files from a directory."""
    return sorted(Path(data_dir).glob("*.part"))


    
def read_run_settings(filename,verbose=True):
    """
    Returns:
        dict: Dictionary with namelist names as keys and their variables as sub-dictionaries.
    """
    import f90nml

    namelists = ["bfield_nml", "run_settings"]
    nml = f90nml.read(filename)


    # Print results
    if verbose:
        for section, variables in nml.items():
            print(f"\n[{section}]")
            for var, value in variables.items():
                print(f"{var} = {value}")
    
    return {name: nml.get(name, {}) for name in namelists}



    

def read_int_pts_file(filename):
    """
    Reads an int_pts.out file and returns a pandas DataFrame.
    
    Parameters:
        filename (str): Path to the file.
    
    Returns:
        pd.DataFrame: DataFrame containing the parsed data.
    """
    column_names = [
        "line_index", "R", "Z", "Phi", "ihit", "ipart", "itri", "i", "Lc", "sin(theta)"
    ]
    
    # Read the file, skipping the first line (header)
    df = pd.read_csv(
        filename,
        sep=r'\s+',  # Replaces deprecated delim_whitespace=True
        skiprows=1,
        names=column_names,
        engine="python"  # Ensures regex separators work correctly
    ) 

    return df


if __name__ == '__main__':
    main()

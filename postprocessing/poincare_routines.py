import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import glob

def read_surface_data(filename):
    """Reads surface data from the file and returns a structured dictionary indexed by slice and surface."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Read the first line with num_surfs and num_slices
    num_surfs, num_slices = map(int, lines[0].split())
    
    data = {
        'slices': {}
    }
    
    index = 1
    for slice_idx in range(num_slices):
        # Read dphi_slice and num_points_this_slice
        dphi, num_points = map(float, lines[index].split())
        num_points = int(num_points)
        index += 1
        
        data['slices'][slice_idx] = {
            'dphi': dphi,
            'num_points': num_points,
            'R': np.zeros((num_surfs, num_points)),
            'Z': np.zeros((num_surfs, num_points))
        }
        
        # Read surface data for each surface
        for surf_idx in range(num_surfs):
            R_values = np.array(list(map(float, lines[index].split())))  # Read R values
            index += 1
            Z_values = np.array(list(map(float, lines[index].split())))  # Read Z values
            index += 1
            
            data['slices'][slice_idx]['R'][surf_idx, :] = R_values
            data['slices'][slice_idx]['Z'][surf_idx, :] = Z_values
    
    return data

def plot_all_surfaces_one_plot(slice_index=0):
    """
    Combines the chosen slice from all 'surface_data.out.####' files into a single plot.
    """
    # Find all rank-based files (e.g. surface_data.out.0000, surface_data.out.0001, ...)
    files = sorted(glob.glob("surface_data.out.[0-9][0-9][0-9][0-9]"))
    if not files:
        print("No 'surface_data.out.####' files found in current directory.")
        return

    # Set up a single figure
    plt.figure(figsize=(8, 6))

    # Define filtering ranges
    Rmax = 10.0
    Rmin = 0.01
    Zmin = -10
    Zmax = 10

    # Loop over each file, read the chosen slice, plot onto the same figure
    for fname in files:
        data = read_surface_data(fname)

        # Ensure the chosen slice_index exists in this file
        if slice_index not in data['slices']:
            print(f"File '{fname}' does not have slice {slice_index}. Skipping.")
            continue
        
        slice_data = data['slices'][slice_index]
        dphi = slice_data['dphi']
        R_surfaces = slice_data['R']
        Z_surfaces = slice_data['Z']

        # Plot each surface from this rank
        for surf_idx in range(R_surfaces.shape[0]):
            R_vals = R_surfaces[surf_idx, :]
            Z_vals = Z_surfaces[surf_idx, :]

            # Filtering
            valid_indices = ((R_vals >= Rmin) & (R_vals <= Rmax) &
                             (Z_vals >= Zmin) & (Z_vals <= Zmax))

            # We'll label it with the filename and surface index
            plt.plot(
                R_vals[valid_indices],
                Z_vals[valid_indices],
                marker='o',
                linestyle='None',
                label=f"{fname} surf {surf_idx+1}"
            )

    plt.xlabel("R")
    plt.ylabel("Z")
    plt.axis("equal")
    plt.title(f"All Surfaces, Combined for Slice {slice_index}")
    plt.legend()
    plt.grid()
    plt.show()

# Optionally, you can just call this function directly:
if __name__ == "__main__":
    plot_all_surfaces_one_plot(slice_index=0)

#!/usr/bin/env python3
import numpy as np
import argparse
import matplotlib.pyplot as plt

def perpendicular_distance(point, start, end):
    """
    Compute the perpendicular distance from a point to a line segment defined by start and end.
    """
    # Handle the degenerate case where start and end are the same.
    if np.allclose(start, end):
        return np.linalg.norm(point - start)
    # Calculate the perpendicular distance using the cross product.
    return np.abs(np.cross(end - start, start - point)) / np.linalg.norm(end - start)

def rdp(points, epsilon):
    """
    Recursively simplify a curve using the Ramer-Douglas-Peucker algorithm.
    
    points: numpy array of shape (N, 2) where each row is a point.
    epsilon: maximum allowed deviation (tolerance) for removed points.
    
    Returns a numpy array of the simplified points.
    """
    if points.shape[0] < 3:
        return points

    # Find the point with the maximum distance from the line between the first and last point.
    distances = np.array([
        perpendicular_distance(points[i], points[0], points[-1])
        for i in range(1, points.shape[0]-1)
    ])
    max_index = np.argmax(distances) + 1  # adjust index since we skipped the first point
    max_distance = distances[max_index - 1]

    if max_distance > epsilon:
        # Recursively simplify the segments.
        first_half = rdp(points[:max_index+1], epsilon)
        second_half = rdp(points[max_index:], epsilon)
        # Concatenate the results, avoiding duplicate of the split point.
        return np.vstack((first_half[:-1], second_half))
    else:
        # The endpoints suffice to represent this segment.
        return np.array([points[0], points[-1]])

def plot_results(original, simplified):
    """
    Plot the original and simplified points.
    
    original: numpy array of original points.
    simplified: numpy array of simplified points.
    """
    plt.figure(figsize=(10, 6))
    # Plot the original points with blue circles connected by lines.
    plt.plot(original[:, 0], original[:, 1], 'bo-', label='Original', markersize=4)
    # Plot the simplified points with red markers connected by lines.
    plt.plot(simplified[:, 0], simplified[:, 1], 'r.-', label='Simplified', markersize=10)
    plt.xlabel("R (mm)")
    plt.ylabel("Z (mm)")
    plt.title("Comparison of Original and Simplified Points")
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    parser = argparse.ArgumentParser(
        description="Simplify a set of points by removing redundant points along straight segments. "
                    "The input file should contain two columns (R and Z in mm)."
    )
    parser.add_argument("-t", "--tolerance", type=float, default=1.0,
                        help="Tolerance for maximum deviation (in mm) allowed for removed points.")
    parser.add_argument("-i", "--input", type=str, default="vvfile.ogr",
                        help="Input file name containing the data (default: vvfile.ogr)")
    parser.add_argument("-o", "--output", type=str, default="simplified.ogr",
                        help="Output file name for simplified data (default: simplified.ogr)")
    parser.add_argument("--plot", action="store_true",
                        help="If set, display a plot comparing the original and simplified points.")
    args = parser.parse_args()

    # Load the data from file; expects two columns per row.
    try:
        data = np.loadtxt(args.input)
    except Exception as e:
        print(f"Error reading {args.input}: {e}")
        return

    if data.ndim != 2 or data.shape[1] != 2:
        print("Input data must have two columns (R and Z).")
        return

    # Simplify the points using the RDP algorithm.
    simplified_points = rdp(data, args.tolerance)

    # Save the simplified data to the output file.
    np.savetxt(args.output, simplified_points, fmt="%.4f")
    print(f"Simplified data written to '{args.output}'.")
    print(f"Original number of points: {data.shape[0]}")
    print(f"Simplified number of points: {simplified_points.shape[0]}")

    # Optionally, display a plot comparing the original and simplified points.
    if args.plot:
        plot_results(data, simplified_points)

if __name__ == '__main__':
    main()

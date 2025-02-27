#!/usr/bin/env python3
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Read simplified.ogr and produce a new part file with a specific format. " +
                    "Coordinates are converted from mm to cm."
    )
    parser.add_argument("-i", "--input", type=str, default="simplified.ogr",
                        help="Input file name containing the simplified points (default: simplified.ogr)")
    parser.add_argument("-o", "--output", type=str, default="new.part",
                        help="Output file name (default: new.part)")
    parser.add_argument("-l", "--label", type=str, default="label",
                        help="Label string to write as the first line (default: 'label')")
    args = parser.parse_args()

    # Read the simplified points from file.
    try:
        data = np.loadtxt(args.input)
    except Exception as e:
        print(f"Error reading {args.input}: {e}")
        return

    # Ensure the data has two columns (R and Z in mm)
    if data.ndim != 2 or data.shape[1] != 2:
        print("Input data must have two columns (R and Z).")
        return

    # Convert points from mm to cm (1 mm = 0.1 cm)
    data_cm = data * 0.1

    # Count the number of points.
    num_points = data_cm.shape[0]

    # Prepare the new file content.
    lines = []
    lines.append(args.label)
    lines.append(f"2 {num_points} 1 0.0 0.0")
    lines.append("1e-08")
    
    # Format each point as "R,Z"
    point_lines = []
    for r, z in data_cm:
        point_lines.append(f"{r},{z}")
    
    lines.extend(point_lines)

    lines.append("359.99999999")

    lines.extend(point_lines)
    
    # Write the output file.
    try:
        with open(args.output, 'w') as f:
            for line in lines:
                f.write(line + "\n")
        print(f"New part file written to '{args.output}'.")
        print(f"Label: {args.label}")
        print(f"Number of points: {num_points}")
    except Exception as e:
        print(f"Error writing {args.output}: {e}")

if __name__ == '__main__':
    main()

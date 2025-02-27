#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import argparse

def main():
    parser = argparse.ArgumentParser(description="Plot an OGR file (two columns).")
    parser.add_argument("-i", "--input", type=str, default="vvfile.ogr",
                        help="Input OGR file name (default: vvfile.ogr)")
    args = parser.parse_args()

    try:
        data = np.loadtxt(args.input)
    except Exception as e:
        print(f"Error reading {args.input}: {e}")
        return

    if data.ndim != 2 or data.shape[1] != 2:
        print("The input file must have exactly two columns.")
        return

    plt.figure(figsize=(8, 6))
    plt.plot(data[:, 0], data[:, 1], 'bo-', label="Data")
    plt.xlabel("R")
    plt.ylabel("Z")
    plt.title(f"Plot of data from {args.input}")
    plt.grid(True)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()

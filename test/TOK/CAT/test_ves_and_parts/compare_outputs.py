#!/usr/bin/env python3

import filecmp
import os
import subprocess
import sys


OUTPUT_FILES = [
    "allparts.out",
    "hitcount.out",
    "hitline.out",
    "int_pts.out",
    "launch_pts.out",
    "part_triangle_mids.out",
    "part_triangles.out",
]

TIMING_FILE = "timing_data.out"


def check_timing_file():
    if not os.path.exists(TIMING_FILE):
        print(f"Error: Missing output file {TIMING_FILE}")
        sys.exit(1)

    with open(TIMING_FILE, "r") as f:
        timing_lines = [line.split() for line in f.readlines()[1:]]

    if len(timing_lines) != 10:
        print(f"Test FAILED: {TIMING_FILE} should contain 10 timing rows.")
        sys.exit(1)

    for expected_index, fields in enumerate(timing_lines, start=1):
        if len(fields) != 3:
            print(f"Test FAILED: {TIMING_FILE} row {expected_index} should contain 3 columns.")
            sys.exit(1)
        if int(fields[0]) != expected_index:
            print(f"Test FAILED: {TIMING_FILE} index {fields[0]} should be {expected_index}.")
            sys.exit(1)


def main():
    proc = subprocess.run(["bash", "./run_test"])
    if proc.returncode != 0:
        print("Error: run_test script returned non-zero exit code.")
        sys.exit(1)

    ref_dir = "./ref"
    for fname in OUTPUT_FILES:
        actual_path = fname
        ref_path = os.path.join(ref_dir, fname)

        if not os.path.exists(actual_path):
            print(f"Error: Missing output file {actual_path}")
            sys.exit(1)
        if not os.path.exists(ref_path):
            print(f"Error: Missing reference file {ref_path}")
            sys.exit(1)

        if not filecmp.cmp(actual_path, ref_path, shallow=False):
            print(f"Test FAILED: {fname} differs from reference.")
            sys.exit(1)

    check_timing_file()

    print("Test PASSED.")
    sys.exit(0)


if __name__ == "__main__":
    main()

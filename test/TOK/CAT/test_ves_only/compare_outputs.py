#!/usr/bin/env python3

import os
import sys
import subprocess
import filecmp

OUTPUT_FILES = [
    "allparts.out",
    "hitcount.out",
    "hitline.out",
    "int_pts.out",
    "launch_pts.out",
    "part_triangle_mids.out",
    "part_triangles.out"
]

TIMING_FILE = "timing_data.out"

def main():
    # Because CTest will run us in test_ves_only/,
    # we can just do relative paths:
    run_test_script = "./run_test"
    ref_dir = "./ref"

    # 1. Run the test script
    proc = subprocess.run([run_test_script])
    if proc.returncode != 0:
        print("Error: run_test script returned non-zero exit code.")
        sys.exit(1)

    # 2. Compare output files to reference
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

    # Timing values vary from run to run, so only check the file shape and line indices.
    timing_path = TIMING_FILE
    if not os.path.exists(timing_path):
        print(f"Error: Missing output file {timing_path}")
        sys.exit(1)

    with open(timing_path, "r") as f:
        timing_lines = [line.split() for line in f.readlines()[1:]]

    if len(timing_lines) != 10:
        print(f"Test FAILED: {timing_path} should contain 10 timing rows.")
        sys.exit(1)

    for expected_index, fields in enumerate(timing_lines, start=1):
        if len(fields) != 3:
            print(f"Test FAILED: {timing_path} row {expected_index} should contain 3 columns.")
            sys.exit(1)
        if int(fields[0]) != expected_index:
            print(f"Test FAILED: {timing_path} index {fields[0]} should be {expected_index}.")
            sys.exit(1)

    print("Test PASSED.")
    sys.exit(0)

if __name__ == "__main__":
    main()

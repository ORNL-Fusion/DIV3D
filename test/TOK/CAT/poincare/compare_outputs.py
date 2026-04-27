#!/usr/bin/env python3

import filecmp
import os
import subprocess
import sys


OUTPUT_FILES = [
    "surface_data.out.0000",
    "surface_data.out.0001",
    "surface_data.out.0002",
    "surface_data.out.0003",
]


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

    print("Test PASSED.")
    sys.exit(0)


if __name__ == "__main__":
    main()

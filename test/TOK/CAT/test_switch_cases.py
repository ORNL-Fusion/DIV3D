#!/usr/bin/env python3

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
DIV3D_EXE = REPO_ROOT / "bin" / "div3d.exe"


def replace_setting(text, key, value):
    lines = text.splitlines()
    in_run_settings = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.lower().startswith("&run_settings"):
            in_run_settings = True
            continue
        if in_run_settings and stripped == "/":
            lines.insert(i, f"{key} = {value}")
            return "\n".join(lines) + "\n"
        if in_run_settings and stripped.startswith(key) and (
            len(stripped) == len(key) or stripped[len(key)] in (" ", "=")
        ):
            lines[i] = f"{key} = {value}"
            return "\n".join(lines) + "\n"
    raise ValueError(f"run_settings namelist not found while setting {key}")


def make_case(tmpdir, name):
    case_root = tmpdir / name / "TOK" / "CAT"
    shutil.copytree(SCRIPT_DIR, case_root)
    return case_root / "test_ves_only"


def run_case(case_dir):
    proc = subprocess.run(
        ["mpirun", "-np", "2", str(DIV3D_EXE)],
        cwd=case_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if proc.returncode != 0:
        print(proc.stdout)
        raise AssertionError(f"case failed in {case_dir}")
    return proc.stdout


def read_surface_count(case_dir):
    with open(case_dir / "surface_line.out", "r") as f:
        f.readline()
        return int(f.readline().strip())


def test_no_hitline(tmpdir):
    case_dir = make_case(tmpdir, "no_hitline")
    settings = (case_dir / "run_settings.nml").read_text()
    settings = replace_setting(settings, "hit_length", "0.0d0")
    (case_dir / "run_settings.nml").write_text(settings)
    run_case(case_dir)

    hitline = case_dir / "hitline.out"
    assert hitline.exists(), "hitline.out should be created"
    assert hitline.stat().st_size == 0, "hitline.out should be empty when hitline is disabled"
    assert (case_dir / "int_pts.out").exists(), "int_pts.out should still be written"
    assert (case_dir / "timing_data.out").exists(), "timing_data.out should still be written"


def test_existing_launch_points(tmpdir):
    case_dir = make_case(tmpdir, "existing_launch")
    shutil.copyfile(case_dir / "ref" / "launch_pts.out", case_dir / "launch_pts.out")

    settings = (case_dir / "run_settings.nml").read_text()
    settings = replace_setting(settings, "trace_surface_opt", ".false.")
    (case_dir / "run_settings.nml").write_text(settings)
    run_case(case_dir)

    with open(case_dir / "launch_pts.out", "r") as f:
        f.readline()
        indices = [int(line.split()[0]) for line in f]
    assert indices == list(range(1, len(indices) + 1)), "external launch indices should remain ordered"


def test_surface_thinning(tmpdir):
    case_dir = make_case(tmpdir, "surface_thinning")
    settings = (case_dir / "run_settings.nml").read_text()
    settings = replace_setting(settings, "npts_surf_out", "100")
    (case_dir / "run_settings.nml").write_text(settings)
    output = run_case(case_dir)

    assert "traced surface points" in output, "screen output should report surface point writing"
    assert read_surface_count(case_dir) == 100, "surface_line.out should contain requested thinned count"


def main():
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        test_no_hitline(tmpdir)
        test_existing_launch_points(tmpdir)
        test_surface_thinning(tmpdir)

    print("Switch-flow tests passed.")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as exc:
        print(f"FAIL: {exc}")
        sys.exit(1)

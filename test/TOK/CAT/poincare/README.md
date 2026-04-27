# TOK CAT Poincare Test

This case traces Poincare surfaces for the axisymmetric CAT tokamak field using
the g-file in `../bfield/g000001.00001_symm`.

The case starts 10 field lines along a line from `(R,Z) = (5.32,0.0)` to
`(4.5,0.0)`, follows them for 100 toroidal transits, and writes three toroidal
slices. The Poincare driver writes one `surface_data.out.####` file per MPI
rank, so the current reference test assumes the four-rank launch in `run_test`.

## Regression Check

CTest runs this case through `compare_outputs.py`, which compares the generated
rank-local output files against `ref/surface_data.out.*`.

To run it directly:

```bash
bash run_test
python3 compare_outputs.py
```

## Visual Check

Use `plot_results.py` for a local sanity plot. This is intentionally not part of
CTest.

```bash
bash run_test
python3 plot_results.py
```

For this axisymmetric case, the default plot includes one panel for each
toroidal slice plus a fourth panel combining all slices.

Useful options:

```bash
python3 plot_results.py --ref
python3 plot_results.py --ref --parts
python3 plot_results.py --slice 0
python3 plot_results.py --ref --save poincare_check.png --no-show
python3 plot_results.py --ref --parts --save poincare_with_parts.png --no-show
python3 plot_results.py --ref --no-combined
```

With `--parts`, the script looks for `.part` files in `../geo` and overlays
them on each toroidal slice. Axisymmetric parts use the first contour for all
slices, matching DIV3D's part-intersection path. Non-axisymmetric parts are not
overlaid yet because DIV3D builds and intersects 3D triangles for that case.

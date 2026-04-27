# TOK CAT Geometry Utilities

This directory contains small geometry files and helper scripts used by the
TOK/CAT tests.

The raw vessel contour is `vvfile.ogr`, a two-column `(R,Z)` contour in mm. The
DIV3D vessel/part files use the local `.part` format with coordinates in cm, so
the usual workflow is:

```bash
python3 plot_ogr.py -i vvfile.ogr
python3 plot_and_simplify_ogr.py -i vvfile.ogr -o simplified.ogr -t 1.0 --plot
python3 convert_ogr_to_AS_part.py -i simplified.ogr -o full_pfc.part -l full_pfc
python3 plot_parts.py
python3 plot_parts.py --save parts_check.png --no-show
```

## Scripts

- `plot_ogr.py`: quick visual check of a two-column OGR contour.
- `plot_parts.py`: plots all `.part` files in this directory in an `(R,Z)`
  view and labels each contour. This is useful for checking which geometry files
  are being used by the tests.
- `plot_and_simplify_ogr.py`: simplifies an OGR contour using the
  Ramer-Douglas-Peucker algorithm. The tolerance is in mm. Use `--plot` to
  compare the original and simplified contours.
- `convert_ogr_to_AS_part.py`: converts a simplified two-column contour from mm
  to cm and writes an axisymmetric `.part` file. The output contains two
  toroidal slices near `0` and `360` degrees.

## Geometry Files

- `vvfile.ogr`: raw vessel contour in mm.
- `full_pfc.part`: axisymmetric vessel contour used by the vessel-only test.
- `box.part`: simple box vessel used by the vessel-plus-parts test.
- `lo_div.part`: small divertor-like part used by `test_ves_and_parts`.

Generated intermediate files such as `simplified.ogr` are local working files
and do not need to be committed unless they are intentionally becoming a test
fixture.

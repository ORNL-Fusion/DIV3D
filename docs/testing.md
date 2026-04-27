# Testing

DIV3D tests are run through CTest from the build directory. The current test
suite includes a reference integration case, switch-flow integration checks, and
Fortran unit tests for core geometry and intersection routines.

## Running Tests

Build the code first, then run CTest:

```bash
cd build
cmake --build .
ctest --output-on-failure
```

To run a single test:

```bash
ctest -R test_intersections --output-on-failure
```

## Test Layout

- `test_ves_only`: vessel-only reference integration case comparing output
  files against files in `test_ves_only/ref`.
- `test/TOK/CAT/test_ves_and_parts`: vessel-plus-parts reference integration
  case comparing output files against files in `test_ves_and_parts/ref`.
- `test/TOK/CAT/poincare`: tokamak Poincare reference integration case
  comparing rank-local `surface_data.out.*` files against files in
  `poincare/ref`.
- `test/TOK/CAT/test_switch_cases.py`: integration checks that run modified
  copies of the `TOK/CAT/test_ves_only` case.
- `test/unit`: Fortran unit tests built against `libdiv3d_core`.

The Fortran unit test sources live outside `src` because they are test-only
programs and helper modules. Production routines are compiled into
`libdiv3d_core`, and the test executables link against that library.

## Unit Tests

### `test_math_routines`

This test covers low-level geometry routines.

`line_seg_facet_int` is tested with a right triangle in the `z = 0` plane:

```text
pa = (0, 0, 0)
pb = (1, 0, 0)
pc = (0, 1, 0)
```

The checks include a segment through the triangle interior, a segment crossing
the plane outside the triangle, a segment parallel to the plane, a vertex hit,
and the `calc_theta = .false.` path.

`int_line_curve` is tested with a closed unit square in `(R,Z)`:

```text
(0,0) -> (1,0) -> (1,1) -> (0,1) -> (0,0)
```

The checks include a single side hit, a miss, and a line crossing the square
twice.

`wrap_phi` is tested for negative angles and angles larger than one period.

### `test_inside_vessel`

This test uses a one-slice square vessel:

```text
R = [0, 1, 1, 0]
Z = [0, 0, 1, 1]
Phi = 0
```

It checks an inside point, an outside point, and a point with `phi = 3 * period`
to exercise phi wrapping in the vessel check.

### `test_read_parts`

This test covers part metadata, part loading, unit conversion, and triangle
generation. The fixture is `test/unit/fixtures/tiny.2d.jpart`, a small
`2 x 3` part with one field period:

```text
ntor = 2
npol = 3
nfp  = 1
```

The first toroidal slice is at `phi = 0 deg`, `Z = 0 cm`; the second is at
`phi = 90 deg`, `Z = 10 cm`; and `R = 100, 110, 120 cm`. The test verifies
that `query_part` reads the dimensions, `load_2d_jpart` converts centimeters to
meters and degrees to radians, and `make_triangles` creates the expected four
triangles for the `2 x 3` surface mesh.

### `test_intersections`

This test directly exercises `check_line_for_intersections`.

The part-hit case uses a large square vessel so the vessel does not interfere,
plus one triangular part in the Cartesian plane `x = 1`:

```text
x = [1, 1, 1]
y = [-0.5, 0.5, 0]
z = [-0.5, -0.5, 0.5]
```

The traced line passes through the part:

```text
R   = [0.5, 1.0, 1.5]
Z   = [0.0, 0.0, 0.0]
Phi = [0.0, 0.0, 0.0]
```

The vessel-hit case uses no parts and a unit-square vessel. The line starts
inside and exits through `R > 1`. This also checks the `nhitline = 0`,
`calc_lc = .false.`, and `calc_theta = .false.` paths.

The no-hit case uses no parts and a larger square vessel. The line remains
inside the vessel, verifying that output variables are initialized cleanly when
there is no intersection.

## Switch-Flow Integration Tests

`test_switch_cases.py` runs temporary copies of the `TOK/CAT/test_ves_only`
case with selected namelist changes:

- `hit_length = 0.0d0`: verifies that `hitline.out` is created empty while
  `int_pts.out` and `timing_data.out` are still written.
- `trace_surface_opt = .false.` with an existing `launch_pts.out`: verifies
  that externally supplied launch-point indices remain ordered.
- `npts_surf_out = 100`: verifies that `surface_line.out` is thinned to the
  requested count and that the screen output reports surface point writing.

## Poincare Integration Test

`test_poincare_tok` runs the `TOK/CAT/poincare` case with four MPI ranks. The
case traces 10 surfaces for 100 transits using the CAT g-file field, saves three
toroidal slices, and compares the four rank-local `surface_data.out.*` files
against references.

The Poincare driver writes one output file per MPI rank, so this test is tied to
the `mpirun -np 4` setting in `test/TOK/CAT/poincare/run_test`.

For a manual visual check, run the local plotting helper:

```bash
cd test/TOK/CAT/poincare
bash run_test
python3 plot_results.py
```

To plot the committed reference files or save a non-interactive PNG:

```bash
python3 plot_results.py --ref
python3 plot_results.py --ref --parts
python3 plot_results.py --ref --save poincare_check.png --no-show
```

This plotting helper is intentionally not registered with CTest.

## Adding Fortran Unit Tests

Add a new Fortran test program under `test/unit` and register it in
`test/unit/CMakeLists.txt` with:

```cmake
add_div3d_unit_test(test_name test_name.f90)
```

Use `test_assert_mod` for assertions. Keep synthetic geometry small and
explicit so failures are easy to reason about.

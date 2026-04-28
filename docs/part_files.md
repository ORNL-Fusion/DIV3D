# Part Files Documentation

This document describes the structure and functionality of the part files used in the DIV3D code. Part files define the geometry of plasma-facing components (PFCs) and their interaction with the magnetic field.

## Overview

Part files define component geometry in radial (`R`), vertical (`Z`), and
toroidal (`Phi`) coordinates. DIV3D processes these files to determine whether
components are axisymmetric (AS) or non-axisymmetric (non-AS), and to filter out
unnecessary intersection checks.

The part header contains an `msym` value. For non-axisymmetric parts this is a
mapping symmetry, not a copy count: DIV3D wraps the part's toroidal coordinates
into the base interval `0 <= Phi < 2*pi/msym` and uses that mapped
representative geometry for intersections. If a physical configuration needs
multiple non-AS copies, include those copies explicitly as separate part files
or explicit toroidal geometry. Setting `msym = 13` does not automatically create
13 limiters.

## Key Steps in Processing Part Files

### 1. Setting R, Z, and Phi Arrays

The `Rparts`, `Zparts`, and `Pparts` arrays store the geometry of each part:
- **Rparts**: Radial coordinates.
- **Zparts**: Vertical coordinates.
- **Pparts**: Toroidal (Phi) coordinates.

For standard `.part` and `.2d.jpart` files, the Phi values are shifted with
`wrap_phi(Phi, 2*pi/msym)`. This mapping is independent of the magnetic-field
period `2*pi/nfp_bfield`, which is used by the field-line following. This allows
part geometry and magnetic-field periodicity to differ, but it also means the
part is only checked where its mapped representative geometry lies.

### 2. Determining Phi Range
For each part:
- The minimum (`Pmins`) and maximum (`Pmaxs`) Phi values are calculated.
- For triangular parts (`part_type = 2`), the Phi range is set to the field period.

### 3. Axisymmetry Check

A part is considered axisymmetric (AS) if the R and Z coordinates are identical
across all toroidal cuts. Finite toroidal extent is allowed; `Pmins` and
`Pmaxs` filter AS intersection checks to the mapped phi range.

If a part is AS, it is flagged as `is_AS_part`. However, this can be overridden using the `force_non_AS` flag.

### 4. Writing Part Files
The processed part data is written to an output file, including the R, Z, and Phi coordinates for all toroidal and poloidal points.

### 5. Vessel File Processing

The vessel geometry is loaded from a separate file and processed similarly to
parts:
- The Phi values are shifted with the vessel file's `msym`.
- The vessel is checked for axisymmetry and flagged as `is_AS_ves` if
  applicable.

## Key Parameters and Flags

- **`ntor`**: Number of toroidal points.
- **`npol`**: Number of poloidal points.
- **`msym`**: Part mapping symmetry. Standard part-file phi values are wrapped
  into `0 <= Phi < 2*pi/msym`. This does not instantiate physical copies.
- **`check_AS_tol`**: Tolerance for axisymmetry checks.
- **`phi_period_tol`**: Tolerance for Phi range checks.
- **`force_non_AS`**: Flag to force non-axisymmetric treatment for parts.
- **`verbose`**: Enables detailed output for debugging.

## Example Output
For each part, the following information is printed if `verbose` is enabled:
- Part type and symmetry.
- Phi range in degrees.
- Axisymmetry status.

Warnings are issued if a part extends beyond the magnetic field period. An
error is issued if adjacent toroidal cuts cross an `msym` mapping boundary after
wrapping, because that would create artificial triangles across the base-period
boundary. Split the part at that boundary or choose an `msym`/placement that
keeps adjacent cuts monotonic after wrapping.

## Notes
- Triangular parts (`part_type = 2`) are always treated as non-axisymmetric.
- The vessel file is processed separately and follows the same axisymmetry checks as parts.
- Non-AS part copies are not generated from `msym`. To model multiple physical
  copies, list multiple files or define the full explicit geometry.

## Related Files
- **Part Files**: Define the geometry of individual components.
- **Vessel File**: Defines the geometry of the vessel.

For more details, refer to the main README file and the source code comments.

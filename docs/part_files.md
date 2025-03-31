# Part Files Documentation

This document describes the structure and functionality of the part files used in the DIV3D code. Part files define the geometry of plasma-facing components (PFCs) and their interaction with the magnetic field.

## Overview

Part files are used to define the 3D geometry of components in terms of their radial (R), vertical (Z), and toroidal (Phi) coordinates. These files are processed to determine whether components are axisymmetric (AS) or non-axisymmetric (non-AS) and to filter out unnecessary intersection checks.

## Key Steps in Processing Part Files

### 1. Setting R, Z, and Phi Arrays
The `Rparts`, `Zparts`, and `Pparts` arrays store the geometry of each part:
- **Rparts**: Radial coordinates.
- **Zparts**: Vertical coordinates.
- **Pparts**: Toroidal (Phi) coordinates.

The Phi values are shifted to the first field period using the `wrap_phi` function to ensure consistency.

### 2. Determining Phi Range
For each part:
- The minimum (`Pmins`) and maximum (`Pmaxs`) Phi values are calculated.
- For triangular parts (`part_type = 2`), the Phi range is set to the field period.

### 3. Axisymmetry Check
A part is considered axisymmetric (AS) if:
- The R and Z coordinates are identical across all toroidal cuts.
- The Phi range matches the magnetic field period.

If a part is AS, it is flagged as `is_AS_part`. However, this can be overridden using the `force_non_AS` flag.

### 4. Writing Part Files
The processed part data is written to an output file, including the R, Z, and Phi coordinates for all toroidal and poloidal points.

### 5. Vessel File Processing
The vessel geometry is loaded from a separate file and processed similarly to parts:
- The Phi values are shifted to the first field period.
- The vessel is checked for axisymmetry and flagged as `is_AS_ves` if applicable.

## Key Parameters and Flags

- **`ntor`**: Number of toroidal points.
- **`npol`**: Number of poloidal points.
- **`msym`**: Magnetic symmetry (number of field periods).
- **`check_AS_tol`**: Tolerance for axisymmetry checks.
- **`phi_period_tol`**: Tolerance for Phi range checks.
- **`force_non_AS`**: Flag to force non-axisymmetric treatment for parts.
- **`verbose`**: Enables detailed output for debugging.

## Example Output
For each part, the following information is printed if `verbose` is enabled:
- Part type and symmetry.
- Phi range in degrees.
- Axisymmetry status.

Warnings are issued if a part extends beyond the magnetic field period.

## Notes
- Triangular parts (`part_type = 2`) are always treated as non-axisymmetric.
- The vessel file is processed separately and follows the same axisymmetry checks as parts.

## Related Files
- **Part Files**: Define the geometry of individual components.
- **Vessel File**: Defines the geometry of the vessel.

For more details, refer to the main README file and the source code comments.
# Namelist Input Parameters for DIV3D

This document describes all parameters available in the two primary namelists used by DIV3D:
- `run_settings.nml`
- `bfield_nml`

---

## `run_settings.nml`

### 🔹 File Paths
| Parameter | Description | Default Value |
|----------|-------------|---------------|
| `fname_plist` | List of part filenames to load | `'parts.list'` |
| `fname_ves` | Vessel geometry part file | `'vessel.part'` |
| `fname_surf` | Output file for surface-traced lines | `'surface_line.out'` |
| `fname_launch` | File with launch points for field lines | `'launch_pts.out'` |
| `fname_parts` | Output file for field-line endpoints/intersections | `'allparts.out'` |
| `fname_hit` | Output for evenly spaced points along hitlines | `'hitline.out'` |
| `fname_intpts` | Output for field line intersections with parts | `'int_pts.out'` |
| `fname_nhit` | Output histogram of hit counts | `'hitcount.out'` |
| `fname_ptri` | Output triangles hit by field lines | `'part_triangles.out'` |
| `fname_ptri_mid` | Output midpoints of hit triangles | `'part_triangle_mids.out'` |

### 🔹 Launch Configuration
| Parameter | Description | Default Value |
|----------|-------------|---------------|
| `Rstart` | Initial radial position of field lines [m] | None |
| `Zstart` | Initial vertical position of field lines [m] | None |
| `Phistart` | Initial toroidal angle [rad] | None |
| `npts_start` | Number of field lines to launch | None |
| `randomize_start_dir` | If true, randomly set initial tracing direction | None |
| `myseed` | Seed for random number generation | None |

### 🔹 Surface Tracing
| Parameter | Description | Default Value |
|----------|-------------|---------------|
| `dphi_line_surf_deg` | Toroidal step size for surface tracing [deg] | None |
| `ntran_surf` | Number of toroidal turns to trace surface lines | None |
| `trace_surface_opt` | If true, enable surface line tracing | None |

### 🔹 Diffusion Tracing
| Parameter | Description | Default Value |
|----------|-------------|---------------|
| `dmag` | Step size for tracing field lines [m] | None |
| `dphi_line_diff_deg` | Toroidal step size for diffusion tracing [deg] | None |
| `ntran_diff` | Number of toroidal turns for diffusion tracing | None |
| `lambda_par` | Parallel step-length multiplier [m]; 0 disables flipping | None |

### 🔹 Output and Postprocessing
| Parameter | Description | Default Value |
|----------|-------------|---------------|
| `hit_length` | Length to sample evenly spaced points along hitlines; <=0 disables | None |
| `lsfi_tol` | Tolerance for last-surface intersection detection | `1e-12` |
| `calc_lc` | If true, compute connection length | `.true.` |
| `calc_theta` | If true, compute poloidal angle along line | `.false.` |
| `quiet_bfield` | Suppress magnetic field logging output | `.true.` |

### 🔹 Geometry Handling
| Parameter | Description | Default Value |
|----------|-------------|---------------|
| `vessel_is_nearest_slice` | If true, use vessel part from nearest toroidal slice | None |
| `vessel_int_is_last_point` | If true, save vessel intersection as final point | `.true.` |

---

## `bfield_nml`

### 🔹 Magnetic Field Source Configuration
| Parameter | Description | Default Value |
|----------|-------------|---------------|
| `rmp_type` | Magnetic field configuration type (e.g., `g`, `coils`, `xdr`, `vmec_coils`) | None |
| `nfp_bfield` | Number of field periods (toroidal symmetry mode) | `0` |

### 🔹 G-file Input
| Parameter | Description | Default Value |
|----------|-------------|---------------|
| `gfile_name` | Path to EFIT G-file to define magnetic equilibrium | `'none'` |

### 🔹 VMEC Coil Input
| Parameter | Description | Default Value |
|----------|-------------|---------------|
| `vmec_coils_file` | Path to VMEC coil file | `'none'` |
| `vmec_extcur_set` | External current settings for coils (array, max 100) | `0.0` |

### 🔹 XDR Format
| Parameter | Description | Default Value |
|----------|-------------|---------------|
| `xdr_fname` | Input filename for XDR format magnetic field grid | None |
| `xdr_check` | Enable/disable checking of XDR file contents | None |
| `xdr_verbose` | Enable/disable verbose XDR reader logging | None |

### 🔹 BGRID Input
| Parameter | Description | Default Value |
|----------|-------------|---------------|
| `bgrid_fname` | Magnetic field grid file (BGRID format) | `'none'` |

---

For detailed examples and execution guidelines, refer to the [main README](../README.md).

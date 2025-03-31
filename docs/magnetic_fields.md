# Magnetic Field Configuration in DIV3D

DIV3D supports several magnetic field representations. These are selected using the `rmp_type` parameter in the `bfield_nml` namelist.

---

## 🔹 Field Configuration Types

| `rmp_type` | Description |
|------------|-------------|
| `g` | Use an EFIT G-file for the magnetic equilibrium. Requires `gfile_name`. |
| `coils` | Use a legacy coil set with current definitions. |
| `vmec_coils` | Use VMEC coil file with optional external currents (`vmec_extcur_set`). |
| `vmec_coils_to_fil` | Converts VMEC coils into filament models for tracing. |
| `xdr` | Use an XDR-format magnetic field file. Requires `xdr_fname`. |
| `bgrid` | Use a binary grid (BGRID) format magnetic field. Requires `bgrid_fname`. |

---

## 🔹 Key Fields and Files

- **`gfile_name`**: EFIT format file used with `rmp_type='g'`
- **`vmec_coils_file`**: File containing VMEC coil configuration
- **`vmec_extcur_set`**: Array of current values for VMEC coils
- **`xdr_fname`**: Path to XDR field data file
- **`bgrid_fname`**: Path to binary field grid

Each method initializes internal field structures (`bfield`, `coil`, etc.) and selects appropriate tracing logic.

For more detail, see: [`setup_bfield_module.f90`](../setup_bfield.f90)

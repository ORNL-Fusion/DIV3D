# Examples for Running DIV3D

This document provides examples of how to use DIV3D for various tasks, including setting up magnetic field grids, creating Poincaré plots, and testing intersections with plasma-facing components (PFCs). The examples are divided into two directories: **TOK** (tokamak-based configurations) and **STELL** (stellarator-based configurations).

---

## TOK Examples

### 1. Magnetic Field from G-File
**Directory**: `test_ves_only`

**Description**: This example uses a g-file to define the magnetic field and performs basic field line tracing.

#### Input Files:
- **Magnetic Field File**: `bfield/g000001.00001_symm`
- **Run Settings**: `test_ves_only/run_settings.nml`

#### Command:
```bash
mpirun -np 2 ./div3d.exe
```

---

### 2. Poincaré Plot
**Directory**: `poincare`

**Description**: This example generates a Poincaré plot using the g-file-based magnetic field.

#### Input Files:
- **Magnetic Field File**: `bfield/g000001.00001_symm`
- **Run Script**: `poincare/run_test`

#### Command:
```bash
mpirun -np 4 ../../../../bin/poincare.exe
```

---

### 3. Intersections with Vessel Only
**Directory**: `test_ves_only`

**Description**: This example tests intersections using only a vessel part.

#### Input Files:
- **Vessel File**: `geo/full_pfc.part`
- **Run Settings**: `test_ves_only/run_settings.nml`

#### Command:
```bash
mpirun -np 2 ./div3d.exe
```

---

### 4. Intersections with Multiple Parts
**Directory**: `test_ves_and_parts`

**Description**: This example tests intersections using multiple parts, including a divertor.

#### Input Files:
- **Parts List**: `test_ves_and_parts/parts.list`
- **Run Settings**: `test_ves_and_parts/run_settings.nml`

#### Command:
```bash
mpirun -np 2 ./div3d.exe
```

---

## STELL Examples

### 1. Magnetic Field from Coils and Poincaré Plot
**Directory**: `poincare_coil_file`

**Description**: This example defines the magnetic field using coil data and generates a Poincaré plot.

#### Input Files:
- **Coil Data**: `coil_file.dat`
- **Run Script**: `poincare_coil_file/run_test`

#### Command:
```bash
mpirun -np 4 ../../../../bin/poincare.exe
```

---

### 2. Create B-Grid File from Coil Data
**Directory**: `make_bgrid`

**Description**: This example creates a B-grid file from coil data for use in subsequent simulations.

#### Input Files:
- **Coil Data**: `coil_file.dat`
- **Run Script**: `make_bgrid/run_test`

#### Command:
```bash
mpirun -np 1 ../../../../bin/create_bgrid_from_coils.exe
```

---

### 3. Poincaré Plot with B-Grid
**Directory**: `poincare_grid_file`

**Description**: This example uses the B-grid file to generate a Poincaré plot.

#### Input Files:
- **B-Grid File**: `bfield/Bgrid_vac_QPS_modular_only_90x150x180`
- **Run Script**: `poincare_grid_file/run_test`

#### Command:
```bash
mpirun -np 4 ../../../../bin/poincare.exe
```

---

### 4. Intersections Using B-Grid (Coming Soon)
**Directory**: `intersections_bgrid`

**Description**: This example will test intersections using the B-grid file. This feature is under development.

---

## Notes
- Ensure all required dependencies (e.g., MPI, HDF5) are installed and accessible.
- Modify the `run_settings.nml` file to adjust parameters such as starting points, step sizes, and tracing options.
- Use the `verbose` flag in the settings file for detailed output during debugging.

For more details, refer to the main [README.md](README.md) file and the documentation for specific input files.

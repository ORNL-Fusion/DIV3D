# DIV3D

DIV3D is a field-line tracing code that follows magnetic field lines in 3D magnetic equilibria and determines their interaction with complex plasma-facing component (PFC) geometry.

This code was developed to assess divertor and baffle loads and design a new high heat flux component, the "scraper element" for the W7-X stellarator.

---

## 📖 Citation

The code should be referenced with the following citiations:
[1] J.D. Lore, et al., "Design and Analysis of Divertor Scraper Elements for the W7-X Stellarator", *IEEE Trans. Plasma Sci.* 42 (2014), 539-544. DOI: [10.1109/TPS.2014.2303649](https://doi.org/10.1109/TPS.2014.2303649)

[2] J.D. Lore, et al., "Modeling and Preparation for Experimental Testing of Heat Fluxes on W7-X Divertor Scraper Elements", *IEEE Trans. Plasma Sci.* 46 (2017), 1387-1392. DOI: [10.1109/TPS.2017.2780624](https://doi.org/10.1109/TPS.2017.2780624)

Comparison to experimental data is given in:

[3] J.D. Lore, et al., "Measurement and modeling of magnetic configurations to mimic overload scenarios in the W7-X stellarator", *Nuclear Fusion* 59 (2019), 066041. DOI: [10.1088/1741-4326/ab18d1](https://doi.org/10.1088/1741-4326/ab18d1)

Other applications:

[4] A. Bader, et al., "Power and Particle Exhaust for the Infinity Two Fusion Pilot Plant", *Journal of Plasma Physics" Published online 2025:1-30. (https://doi.org/10.1017/S0022377825000273)

---

## 🔧 Building the Code

To build the code, follow these steps:

1. Navigate to the `build` directory:
   ```bash
   cd build
   ```

2. Run the `setup_cmake.sh` script to configure the build. You can specify the build type (`Debug` or `Release`) and the compiler (`GNU` or `IntelLLVM`). If no arguments are provided, the default is `Release` with the GNU compiler.

   Examples:
   - Default build (Release with GNU compiler):
     ```bash
     ./setup_cmake.sh
     ```
   - Debug build with GNU compiler:
     ```bash
     ./setup_cmake.sh Debug GNU
     ```

3. The script will automatically run `cmake` and `make` to build the project. If you need to clean the build files, use the `clean_cmake.sh` script:
   ```bash
   ./clean_cmake.sh
   ```

### Notes:
- The `setup_cmake.sh` script simplifies the build process by setting machine-specific paths and configurations.
- The code has been successfully compiled using both `gfortran` and Intel OneAPI compilers.
- Ensure that the required dependencies are installed and accessible on your system.
- If HDF5 (or LAPACK) is not detected by CMake, check the default installation path.

### Dependencies:
- Fortran compiler with MPI support (e.g., `gfortran` + `openmpi`)
- CMake
- HDF5, LAPACK, and MPI libraries

---

## ▶️ Running DIV3D

Use an MPI launcher such as:

```bash
mpirun -np 4 ./div3d.exe
```

### Notes:
Parameters are read from `run_settings.nml`, which contains the namelists `run_settings` and `bfield_nml`.

---

## `run_settings` namelist

📘 See full parameter reference: [`docs/input_namelists.md`](docs/input_namelists.md)

---

## `bfield_nml` namelist

Defines the magnetic field source and configuration. Supports EFIT g-files, VMEC coils, XDR and BGRID formats.

📘 See details: [`docs/magnetic_fields.md`](docs/magnetic_fields.md)

## Part files and parts.list

📘 See details: [`docs/part_files.md`](docs/part_files.md)

---

## 📤 Output Files

Outputs include field line endpoints, intersection points, hitline samplings, and diagnostic files.

📘 See descriptions: [`docs/outputs.md`](docs/outputs.md)

---

## 📘 Examples

For detailed examples of how to run DIV3D, including field line tracing, Poincaré plots, and intersection tests, see the [Examples Documentation](docs/examples.md).

---

## Code Flow Overview

1. **Input Parsing**: DIV3D reads input files such as `run_settings.nml`, magnetic field files, and part files.
2. **Magnetic Field Processing**: The magnetic field is loaded and processed based on the input format (e.g., g-files, VMEC coils, or BGRID).
3. **Field Line Tracing**: Field lines are traced using the specified starting points and step sizes.
4. **Intersection Testing**: Intersections with plasma-facing components (PFCs) and vessels are calculated.
5. **Output Generation**: Results such as field line endpoints, hitlines, and diagnostic data are written to output files.

---

## Screen Output

During execution, DIV3D provides output to the terminal. 
Regular "part" hit: [i,ipart,P] -> P is the intersection point in X,Y,Z (m)  
Vessel hit:         [i,P]       -> P is the intersection point in X,Y,Z (m)

---

## Input Files

- **`run_settings.nml`**: Contains simulation parameters. See [`docs/input_namelists.md`](docs/input_namelists.md).
- **Magnetic Field Files**: Supported formats include:
  - **G-Files**: EFIT equilibrium files.
  - **VMEC Coils**: Coil-based magnetic field data.
  - **BGRID**: Precomputed magnetic field grids.
  See [`docs/magnetic_fields.md`](docs/magnetic_fields.md) for details.
- **Part Files**: Define the geometry of plasma-facing components. See [`docs/part_files.md`](docs/part_files.md).

---

## Output Files

- **Field Line Endpoints**: Contains the final positions of traced field lines.
- **Intersection Points**: Lists where field lines intersect with parts or vessels.
- **Diagnostics**: Includes hitline samplings and other diagnostic data.

See [`docs/outputs.md`](docs/outputs.md) for detailed descriptions.

---

## Testing (Stub)

Testing instructions will be added here. For now, refer to the examples in [`docs/examples.md`](docs/examples.md) and compare outputs to reference files in the `test_ves_only/ref` directory.

---

## 🧑‍💻 Contributing

Contributions welcome! Please submit issues or pull requests.

---

## 📄 License

This project is licensed under the MIT license. See the [LICENSE](LICENSE) file for details.

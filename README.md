# DIV3D
DIV3D field line followning and intersection code.  
2010 - 2024  
Jeremy D. Lore

This code was developed to assess divertor and baffle loads and design a new high heat flux component, the "scraper element" for the W7-X stellarator [3].

The code should be referenced with this citiations:

[1] J.D. Lore, et al., "Design and Analysis of Divertor Scraper Elements for the W7-X Stellarator", IEEE Trans. Plasma Sci. 42 (2014), 539-544. DOI: 10.1109/TPS.2014.2303649

[2] J.D. Lore, et al., "Modeling and Preparation for Experimental Testing of Heat Fluxes on W7-X Divertor Scraper Elements", IEEE Trans. Plasma Sci. 46 (2017), 1387-1392. DOI: 10.1109/TPS.2017.2780624

Comparison to experimental data is given in:

[3] J.D. Lore, et al., "Measurement and modeling of magnetic configurations to mimic overload scenarios in the W7-X stellarator", Nuclear Fusion 59 (2019), 066041. DOI: 10.1088/1741-4326/ab18d1

To build: 
1) First build libbjdl from this repository: https://github.com/ORNL-Fusion/util-library/tree/master/fortran/bfield_library_jdl.  See build/setup_cmake.sh
2) Add MACHINE_ID to build/setup_cmake.sh. Set correct location of LIBBDIR and fortran mpi compiler.
3) from ./build run setup_cmake.sh
4) If hdf5 (or lapack) is not detected by cmake check default path (/usr/lib/x86_64-linux-gnu/hdf5/openmpi)  

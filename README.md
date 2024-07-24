# DIV3D
DIV3D field line followning and intersection code.  
2010 - 2024  
Jeremy D. Lore

## Description
This code was developed to assess divertor and baffle loads and design a new high heat flux component, the "scraper element" for the W7-X stellarator [3].

The code should be referenced with this citiations:

[1] J.D. Lore, et al., "Design and Analysis of Divertor Scraper Elements for the W7-X Stellarator", IEEE Trans. Plasma Sci. 42 (2014), 539-544. DOI: 10.1109/TPS.2014.2303649

[2] J.D. Lore, et al., "Modeling and Preparation for Experimental Testing of Heat Fluxes on W7-X Divertor Scraper Elements", IEEE Trans. Plasma Sci. 46 (2017), 1387-1392. DOI: 10.1109/TPS.2017.2780624

Comparison to experimental data is given in:

[3] J.D. Lore, et al., "Measurement and modeling of magnetic configurations to mimic overload scenarios in the W7-X stellarator", Nuclear Fusion 59 (2019), 066041. DOI: 10.1088/1741-4326/ab18d1

## Compilation
To build: 
1) Add MACHINE_ID to build/setup_cmake.sh. Set correct location of fortran mpi compiler.
2) from ./build run setup_cmake.sh
3) If hdf5 (or lapack) is not detected by cmake check default path (/usr/lib/x86_64-linux-gnu/hdf5/openmpi)  

## I/O Description

### Input files (names can be set in run_settings)
1) run_settings.nml   (run_settings namelist)
2) parts.list         (fname_plist)
3) vessel part        (fname_ves)
4) bfield description (bfield_nml namelist)

### Output files (names can be set in run_settings)
1) surface_line.out (fname_surf)
2) allparts.out     (fname_parts)
3) hitline.out      (fname_hit)
4) int_pts.out      (fname_intpts)
Row by row information of points that hit "parts".
R (m) | Z (m) | Phi (rad) | ihit | ipart | itri | i
R,Z,Phi -> Int point in mapped periodic section
ihit    -> (1) hit a "part", (2) hit the vessel, (0) no intersection.
itri    -> Index of intersected triangle
i       -> Index along field line

5) hitcount.out     (fname_nhit)
6) part_triangles.out (fname_ptri)
7) part_triangle_mids.out (fname_ptri_mid)

### Screen output
Regular "part" hit: [i,ipart,P] -> P is the intersection point in X,Y,Z (m)
Vessel hit:         [i,R,Z,P]   -> R,Z,P are the hit coordinates in (m,rad)
# DIV3D
DIV3D field line followning and intersection code.  
Jeremy D. Lore (2010 - 2024)

## Description
This code was developed to assess divertor and baffle loads and design a new high heat flux component, the "scraper element" for the W7-X stellarator.

The code should be referenced with the following citiations:  

[1] J.D. Lore, et al., "Design and Analysis of Divertor Scraper Elements for the W7-X Stellarator", IEEE Trans. Plasma Sci. 42 (2014), 539-544. DOI: [10.1109/TPS.2014.2303649](https://doi.org/10.1109/TPS.2014.2303649)

[2] J.D. Lore, et al., "Modeling and Preparation for Experimental Testing of Heat Fluxes on W7-X Divertor Scraper Elements", IEEE Trans. Plasma Sci. 46 (2017), 1387-1392. DOI: 10.1109/TPS.2017.2780624(https://doi.org/10.1109/TPS.2017.2780624)  

Comparison to experimental data is given in:  

[3] J.D. Lore, et al., "Measurement and modeling of magnetic configurations to mimic overload scenarios in the W7-X stellarator", Nuclear Fusion 59 (2019), 066041. DOI: 10.1088/1741-4326/ab18d1(https://doi.org/10.1088/1741-4326/ab18d1)    

## Compilation
To build: 
1) Add MACHINE_ID to build/setup_cmake.sh. Set correct location of fortran mpi compiler.
2) from ./build run setup_cmake.sh
3) If hdf5 (or lapack) is not detected by cmake check default path (/usr/lib/x86_64-linux-gnu/hdf5/openmpi)  

## I/O Description

### Input files 
### 1) run_settings.nml
#### This file contains namelists for the magnetic field description (bfield_nml namelist) and DIV3D run settings (run_settings namelist)
* [bfield namelist](https://github.com/ORNL-Fusion/DIV3D/README.md#bfield-namelist-description)  
* [run_settings namelist](https://github.com/ORNL-Fusion/DIV3D/README.md#run_settings-namelist-description)  

### 2) parts.list  
#### This file contains the description of the components to be checked for intersection.
### 3) vessel part  
#### This file contains the description of the vessel

### Output files  
#### The names of these files can be set in the [run_settings](https://github.com/ORNL-Fusion/DIV3D/README.md#run_settings-namelist-description) namelist
### 1) surface_line.out  
#### This file contains 
### 2) allparts.out  
#### This file contains 
### 3) hitline.out  
#### This file contains 
### 4) int_pts.out  
#### This file contains 
```
Row by row information of points that hit "parts".

R (m) | Z (m) | Phi (rad) | ihit | ipart | itri | i  
R,Z,Phi -> Int point in mapped periodic section  
ihit    -> (1) hit a "part", (2) hit the vessel, (0) no intersection.  
itri    -> Index of intersected triangle  
i       -> Index along field line  
```
### 6) hitcount.out     (fname_nhit)
#### This file contains 
### 7) part_triangles.out (fname_ptri)
#### This file contains 
### 8) part_triangle_mids.out (fname_ptri_mid)
#### This file contains 
### 9) launch_pts.out (fname_launch)
#### This file contains 
```
Row by row information of points where fieldlines are initiated. 
number of points
R,Z,Phi (m, radians)
```

### Screen output
Regular "part" hit: [i,ipart,P] -> P is the intersection point in X,Y,Z (m)  
Vessel hit:         [i,P]       -> P is the intersection point in X,Y,Z (m)

### Bfield namelist description
```&bfield_nml
! Several options are available, see bfield_module.F90
! Coils defined in filaments
! B components on grid
rmp_type = 'bgrid'
bgrid_fname = 'path/to/bgrid'
! VMEC coils with extcur
! xdr grid format
/
```

### run_settings namelist description
```
fname_plist    = 'parts.list'
fname_ves      = 'path/to/vessel.part'
fname_surf     = 'surface_line.out'
fname_launch   = 'launch_pts.out'
fname_parts    = 'allparts.out'
fname_hit      = 'hitline.out'
fname_intpts   = 'int_pts.out'
fname_nhit     = 'hitcount.out'
fname_ptri     = 'part_triangles.out'
fname_ptri_mid = 'part_triangle_mids.out'

nfp = 3

Rstart   = 4.d0 
Zstart   = 0.d0
Phistart = 0.d0
dphi_line_surf_deg = 0.1d0
ntran_surf = 20

npts_start = 10
dmag = 3.155d-4
dphi_line_diff_deg = 0.5d0
ntran_diff = 1000

trace_surface_opt=.true.
myseed = 110115

!hit_length = 1.0d0
hit_length = -1.0d0

lsfi_tol = 1.d-12

```


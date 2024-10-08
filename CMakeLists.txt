cmake_minimum_required(VERSION 2.8)
project ( div3d )
enable_language (Fortran)

# To compile without XDR routines set variable USE_FXDR to 0 below
#
# Setup is now automated by:
# 1) cd build
# 2) ./cmake_setup.sh [debug]
#
# Debug and release flags are set below.
#
# 1) cd build
# 2) cmake -DCMAKE_BUILD_TYPE=Release .. 
# 3) make install
#
# To see the commands during compilation use "make VERBOSE=1"
#
# To compile with debug options:
#
# 2) cmake -DCMAKE_BUILD_TYPE=Debug ..
# 

# XDR options
set (USE_FXDR 1)

# Compiler flags
if (CMAKE_COMPILER_IS_GNUCC)
    # gfortran
  Message ("Using gfortran flags")
  set (CMAKE_Fortran_FLAGS_RELEASE "-cpp -O3")
#  set (CMAKE_Fortran_FLAGS_DEBUG "-g -Og -frange-check -fbounds-check -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fcheck=all -fbacktrace")
  set (CMAKE_Fortran_FLAGS_DEBUG "-cpp -g -Og")  
  #-ffpe-trap=zero,underflow,overflow,invalid)
else ()
  Message ("Did not recognize compiler -- no flags set!")  
endif ()
#Message ("Compiling with these flags: ${CMAKE_Fortran_FLAGS}")

# Handle build setting
message("Build type is: ${CMAKE_BUILD_TYPE}")
if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
  message("\n\n-----------> WARNING: No build type specified!!  No compiler flags will be set!\n\n")
endif ()

find_library(H5 NAMES hdf5 PATHS /usr/lib/x86_64-linux-gnu/hdf5/openmpi)

# Set up directories and options
set (SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Handle FXDR settings
if (USE_FXDR MATCHES 1)
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrinit.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrint.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrclose.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrreal.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrdouble.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrrmat.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrdmat.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrimat.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrstring.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrreal64.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrrewind.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrrmat64.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/cxdrshort.c")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/fxdr/fxdr_module.f90")
  set (EXTRA_FILES "${EXTRA_FILES};${SRC_DIR}/xdr_routines_mod.f90")
  message ("Using FXDR")
  add_definitions (-DHAVE_FXDR)
endif ()


add_executable( div3d
  ${SRC_DIR}/bgrid_routines_mod.f90
  ${SRC_DIR}/biotsavart_module.f90  
  ${SRC_DIR}/bfield_module.f90
  ${SRC_DIR}/bspline90_22.f90  
  ${SRC_DIR}/coil_typedef.f90
  ${SRC_DIR}/diffuse_lines3.f90  
  ${SRC_DIR}/div3d_follow_and_int.F90  
  ${SRC_DIR}/fieldline_follow_mod.F90
  ${SRC_DIR}/find_nearby_triangles.f90
  ${SRC_DIR}/g3d_module.f90    
  ${SRC_DIR}/init_points_line.f90
  ${SRC_DIR}/init_random_seed.f90  
  ${SRC_DIR}/inside_vessel_mod.f90
  ${SRC_DIR}/io_unit_spec.f90  
  ${SRC_DIR}/kind_mod.f90
  ${SRC_DIR}/math_routines_mod.f90
  ${SRC_DIR}/parallel_mod.f90  
  ${SRC_DIR}/phys_const.f90
  ${SRC_DIR}/read_parts_mod.f90
  ${SRC_DIR}/run_settings_namelist.f90
  ${SRC_DIR}/trace_surface.f90
  ${SRC_DIR}/setup_bfield.f90
  ${SRC_DIR}/VMEC_routines_mod.f90  
  ${SRC_DIR}/LIBSTELL_LIGHT/tolower.f90
  ${SRC_DIR}/LIBSTELL_LIGHT/stel_kinds.f
  ${SRC_DIR}/LIBSTELL_LIGHT/safe_open_mod.f90
  ${SRC_DIR}/LIBSTELL_LIGHT/biotsavart.f
  ${SRC_DIR}/LIBSTELL_LIGHT/bsc_T.f
  ${EXTRA_FILES}  
  )

target_link_libraries (div3d ${H5} z stdc++ lapack)

# Installation
set_target_properties( div3d PROPERTIES OUTPUT_NAME "div3d.exe")
set_target_properties( div3d PROPERTIES DEBUG_OUTPUT_NAME "div3d.exe")
install( DIRECTORY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin )
install( TARGETS div3d
  RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)

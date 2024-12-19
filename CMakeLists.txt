cmake_minimum_required(VERSION 3.5)
project(div3d LANGUAGES C Fortran)

# Specify the Fortran standard
set(CMAKE_Fortran_STANDARD 90)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Default to Release build if none specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  message(STATUS "No build type specified. Defaulting to Release.")
endif()

# Compiler flags
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  message(STATUS "Using gfortran compiler flags")
  set(CMAKE_Fortran_FLAGS_RELEASE " -cpp -O3")
  set(CMAKE_Fortran_FLAGS_DEBUG " -g -cpp -Og -Wrealloc-lhs -fbounds-check -frange-check -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fcheck=all -fbacktrace")
else()
  message(WARNING "Unrecognized Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}. Compiler-specific flags will not be set.")
endif()

# Incorporate HDF5_ROOT and NetCDF_ROOT as hints if provided
set(HDF5_SEARCH_PATHS)
if(HDF5_ROOT)
  message(STATUS "HDF5 root set to: ${HDF5_ROOT}")
  list(APPEND HDF5_SEARCH_PATHS "${HDF5_ROOT}")
endif()

set(NetCDF_SEARCH_PATHS)
if(NetCDF_ROOT)
  message(STATUS "NetCDF root set to: ${NetCDF_ROOT}")
  list(APPEND NetCDF_SEARCH_PATHS "${NetCDF_ROOT}")
endif()

# Locate the HDF5 library
find_library(H5 NAMES hdf5 PATHS ${HDF5_SEARCH_PATHS} /usr/lib/x86_64-linux-gnu/hdf5/openmpi)
if(NOT H5)
  message(FATAL_ERROR "HDF5 library not found. Please install HDF5 or specify its path.")
endif()

# Locate the NetCDF library
find_library(NetCDF NAMES netcdff PATHS ${NetCDF_SEARCH_PATHS} /usr/lib/x86_64-linux-gnu)
if(NOT NetCDF)
  message(FATAL_ERROR "NetCDF library not found. Please install NetCDF or specify its path.")
endif()

# Locate the NetCDF Fortran module directory (netcdf.mod)
find_path(NetCDF_INCLUDE_DIR netcdf.mod 
  PATHS ${NetCDF_SEARCH_PATHS} /usr/include /usr/local/include /usr/lib/x86_64-linux-gnu/include)
if(NOT NetCDF_INCLUDE_DIR)
  message(FATAL_ERROR "NetCDF module directory (netcdf.mod) not found. Please install NetCDF or specify its path.")
endif()

# Add NetCDF include directory
include_directories(${NetCDF_INCLUDE_DIR})

# Directories
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(BIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../bin)

# Ensure binaries are built in ../bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

# FXDR-related files (always included)
set(FXDR_FILES
    ${SRC_DIR}/fxdr/cxdrinit.c
    ${SRC_DIR}/fxdr/cxdrint.c
    ${SRC_DIR}/fxdr/cxdrclose.c
    ${SRC_DIR}/fxdr/cxdrreal.c
    ${SRC_DIR}/fxdr/cxdrdouble.c
    ${SRC_DIR}/fxdr/cxdrrmat.c
    ${SRC_DIR}/fxdr/cxdrdmat.c
    ${SRC_DIR}/fxdr/cxdrimat.c
    ${SRC_DIR}/fxdr/cxdrstring.c
    ${SRC_DIR}/fxdr/cxdrreal64.c
    ${SRC_DIR}/fxdr/cxdrrewind.c
    ${SRC_DIR}/fxdr/cxdrrmat64.c
    ${SRC_DIR}/fxdr/cxdrshort.c
    ${SRC_DIR}/fxdr/fxdr_module.f90
    ${SRC_DIR}/xdr_routines_mod.f90
)

# Create the div3d executable
add_executable(div3d
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
  ${SRC_DIR}/initialize_bfield_div3d.F90
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
  ${FXDR_FILES}
)

# Define HAVE_FXDR since FXDR is always used
target_compile_definitions(div3d PRIVATE HAVE_FXDR)

# Link div3d with HDF5, NetCDF, and other required libs
target_link_libraries(div3d PRIVATE ${H5} z stdc++ lapack ${NetCDF})

# Create the convert_ascii_bgrid_to_netcdf executable
add_executable(convert_ascii_bgrid_to_netcdf
  ${SRC_DIR}/kind_mod.f90
  ${SRC_DIR}/bgrid_routines_mod.f90
  ${SRC_DIR}/convert_ascii_bgrid_to_netcdf.f90
  ${SRC_DIR}/parallel_mod.f90
  ${SRC_DIR}/phys_const.f90  
  )

# Set file names
set_target_properties(div3d PROPERTIES OUTPUT_NAME "div3d.exe")
set_target_properties(convert_ascii_bgrid_to_netcdf PROPERTIES OUTPUT_NAME "convert_ascii_bgrid_to_netcdf.exe")



# Link with NetCDF
target_link_libraries(convert_ascii_bgrid_to_netcdf PRIVATE ${NetCDF})

# Set the install prefix to the directory above the build directory
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/.." CACHE PATH "Install prefix" FORCE)

# Now you can install targets into 'bin' and they will end up in ../bin relative to the build directory.
install(TARGETS div3d convert_ascii_bgrid_to_netcdf
  RUNTIME DESTINATION bin)

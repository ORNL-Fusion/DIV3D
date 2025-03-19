# Set source and binary directories
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(BIN_DIR ${CMAKE_BINARY_DIR}/../bin)

# Ensure binaries are built in ../bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

# Add fxdr library
add_subdirectory(fxdr)
add_subdirectory(LIBSTELL_LIGHT)

# Define a shared library for common modules
add_library(libcommon STATIC
  bgrid_routines_mod.f90
  biotsavart_module.f90
  bfield_module.f90
  bspline90_22.f90
  coil_typedef.f90
  fieldline_follow_mod.F90
  g3d_module.f90  
  init_random_seed.f90
  initialize_bfield_div3d.F90
  io_unit_spec.f90  
  kind_mod.f90
  math_routines_mod.f90  
  parallel_mod.f90
  phys_const.f90  
  run_settings_namelist.f90
  trace_surface.f90
  timing.f90
  setup_bfield.f90
  VMEC_routines_mod.f90
  xdr_routines_mod.f90
)

# Apply the USE_MPIF08 definition globally to the shared library
if(USE_MPIF08)
    target_compile_definitions(libcommon PRIVATE USE_MPIF08)
endif()

# Include directories for shared library
target_include_directories(libcommon PUBLIC ${SRC_DIR} ${CMAKE_CURRENT_BINARY_DIR}/fxdr ${CMAKE_CURRENT_BINARY_DIR}/LIBSTELL_LIGHT)

# libcommon dependencies
target_link_libraries(libcommon PRIVATE libfxdr)
target_link_libraries(libcommon PRIVATE libstell)

# -------------------------------
# Define executables
# -------------------------------

# -------------------------------
# div3d executable
# -------------------------------
add_executable(div3d
  diffuse_lines3.f90
  div3d_follow_and_int.F90
  find_nearby_triangles.f90
  init_points_line.f90
  inside_vessel_mod.f90
  output_routines.f90
  read_parts_mod.f90
  )
set_target_properties(div3d PROPERTIES OUTPUT_NAME "div3d.exe")

# Apply USE_MPIF08 to div3d
if(USE_MPIF08)
    target_compile_definitions(div3d PRIVATE USE_MPIF08)
endif()

# Link libraries to div3d
target_link_libraries(div3d PRIVATE libcommon libfxdr libstell ${NetCDF_Fortran_LIBRARY} ${NetCDF_C_LIBRARY} z stdc++ lapack)

# --------------------------------------------
# convert_ascii_bgrid_to_netcdf executable
# -------------------------------------------
add_executable(convert_ascii_bgrid_to_netcdf
  util/convert_ascii_bgrid_to_netcdf.f90
)
set_target_properties(convert_ascii_bgrid_to_netcdf PROPERTIES OUTPUT_NAME "convert_ascii_bgrid_to_netcdf.exe")

# Apply USE_MPIF08 to convert_ascii_bgrid_to_netcdf
if(USE_MPIF08)
    target_compile_definitions(convert_ascii_bgrid_to_netcdf PRIVATE USE_MPIF08)
endif()

# Link libraries to convert_ascii_bgrid_to_netcdf
target_link_libraries(convert_ascii_bgrid_to_netcdf PRIVATE libcommon ${NetCDF_Fortran_LIBRARY} ${NetCDF_C_LIBRARY})

# -------------------------------
# poincare executable
# -------------------------------
add_executable(poincare
  poincare/poincare_driver.F90
  poincare/poincare_namelist.f90
)
set_target_properties(poincare PROPERTIES OUTPUT_NAME "poincare.exe")

# Apply USE_MPIF08 to poincare
if(USE_MPIF08)
    target_compile_definitions(poincare PRIVATE USE_MPIF08)
endif()

# Link libraries to poincare
target_link_libraries(poincare PRIVATE libcommon libfxdr libstell ${NetCDF_Fortran_LIBRARY} ${NetCDF_C_LIBRARY} z stdc++ lapack)


# -------------------------------
# create_bgrid executable
# -------------------------------
add_executable(create_bgrid_from_coils
  util/create_bgrid_from_coils.f90
  util/create_bgrid_from_coils_namelist.f90
)
set_target_properties(create_bgrid_from_coils PROPERTIES OUTPUT_NAME "create_bgrid_from_coils.exe")

# Apply USE_MPIF08 to create_bgrid_from_coils
if(USE_MPIF08)
    target_compile_definitions(create_bgrid_from_coils PRIVATE USE_MPIF08)
endif()

# Link libraries to create_bgrid_from_coils
target_link_libraries(create_bgrid_from_coils PRIVATE libcommon libstell ${NetCDF_Fortran_LIBRARY} ${NetCDF_C_LIBRARY} z stdc++ lapack)

# Define the fxdr library
add_library(libstell STATIC
  tolower.f90
  stel_kinds.f
  safe_open_mod.f90
  biotsavart.f
  bsc_T.f
)

# Set output directory for Fortran module files
set_target_properties(libstell PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Include directories for this library
target_include_directories(libstell PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

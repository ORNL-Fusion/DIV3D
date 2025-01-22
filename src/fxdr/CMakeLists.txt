# Define the fxdr library
add_library(libfxdr STATIC
  cxdrinit.c
  cxdrint.c
  cxdrclose.c
  cxdrreal.c
  cxdrdouble.c
  cxdrrmat.c
  cxdrdmat.c
  cxdrimat.c
  cxdrstring.c
  cxdrreal64.c
  cxdrrewind.c
  cxdrrmat64.c
  cxdrshort.c
  fxdr_module.f90
)

# Set output directory for Fortran module files
set_target_properties(libfxdr PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Check for rpc/rpc.h
include(CheckIncludeFiles)
check_include_files("rpc/rpc.h" HAVE_RPC_RPC_H)

if(NOT HAVE_RPC_RPC_H)
    message(STATUS "rpc/rpc.h not found. Using TIRPC...")
    target_include_directories(libfxdr PUBLIC /usr/include/tirpc)
    target_compile_definitions(libfxdr PUBLIC USE_TIRPC)
    target_link_libraries(libfxdr PUBLIC tirpc)
else()
    message(STATUS "Using system rpc/rpc.h.")
endif()

# Include directories for this library
target_include_directories(libfxdr PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

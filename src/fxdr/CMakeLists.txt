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

# Check for a usable RPC header.
include(CheckIncludeFiles)
check_include_files("rpc/rpc.h" HAVE_RPC_RPC_H)

if(HAVE_RPC_RPC_H)
    message(STATUS "Using system rpc/rpc.h.")
else()
    find_path(TIRPC_INCLUDE_DIR
      NAMES rpc/rpc.h
      PATH_SUFFIXES tirpc
    )
    find_library(TIRPC_LIBRARY NAMES tirpc)

    if(TIRPC_INCLUDE_DIR AND TIRPC_LIBRARY)
        message(STATUS "rpc/rpc.h not found. Using TIRPC from ${TIRPC_INCLUDE_DIR}.")
        target_include_directories(libfxdr PUBLIC ${TIRPC_INCLUDE_DIR})
        target_link_libraries(libfxdr PUBLIC ${TIRPC_LIBRARY})
    else()
        message(FATAL_ERROR
          "Neither system rpc/rpc.h nor a usable TIRPC installation was found. "
          "Please install the RPC/TIRPC development headers and library."
        )
    endif()
endif()

# Include directories for this library
target_include_directories(libfxdr PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

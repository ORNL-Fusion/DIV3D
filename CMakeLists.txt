cmake_minimum_required(VERSION 2.8)
project ( div3d )
enable_language (Fortran)

#
# Setup is now automated by: (replaces instructions below)
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
set (USE_FXDR 0)

# Compiler flags
if (CMAKE_COMPILER_IS_GNUCC)
    # gfortran
  Message ("Using gfortran flags")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
#  set (CMAKE_Fortran_FLAGS_DEBUG "-g -Og -frange-check -fbounds-check -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fcheck=all -fbacktrace")
  set (CMAKE_Fortran_FLAGS_DEBUG "-g -Og")  
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

# handle library from util-library
set ( LIBBJDL_LIB_PATH ${LIBBJDL_DIR}/lib ) 
set ( LIBBJDL_INC_PATH ${LIBBJDL_DIR}/include )
message ("Using libbjdl files from: " ${LIBBJDL_DIR})
find_library( BJDL_F NAMES bjdl PATHS ${LIBBJDL_LIB_PATH} NO_DEFAULT_PATH)
#  include_directories (${LIBBJDL_INC_PATH})
if (BJDL_F MATCHES NOTFOUND)
    message( "\n\nERROR: -----------------> Could not find libbjdl! -- Is it installed?\n\n")
 endif ()
include_directories(${LIBBJDL_INC_PATH})
add_definitions(-DHAVE_BJDL)

find_library(H5 NAMES hdf5 PATHS /usr/lib/x86_64-linux-gnu/hdf5/openmpi)

# Set up directories and options
set (SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

if (USE_FXDR MATCHES 1)
  message ("Using FXDR")
  add_definitions (-DHAVE_FXDR)
endif ()

add_executable( div3d
  ${SRC_DIR}/inside_vessel_mod.f90
  ${SRC_DIR}/kind_mod.f90
  ${SRC_DIR}/math_routines_mod.f90
  ${SRC_DIR}/io_unit_spec.f90
  ${SRC_DIR}/div3d_follow_and_int.F90
  ${SRC_DIR}/init_points_line.f90
  ${SRC_DIR}/trace_surface.f90
  ${SRC_DIR}/init_random_seed.f90
  ${SRC_DIR}/read_parts_mod.f90
  ${SRC_DIR}/diffuse_lines3.f90
  ${SRC_DIR}/parallel_mod.f90
  ${SRC_DIR}/find_nearby_triangles.f90
  )

target_link_libraries (div3d ${BJDL_F} ${H5} z stdc++ lapack)


# Installation
set_target_properties( div3d PROPERTIES OUTPUT_NAME "div3d.exe")
set_target_properties( div3d PROPERTIES DEBUG_OUTPUT_NAME "div3d.exe")
install( DIRECTORY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin )
install( TARGETS div3d
  RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)

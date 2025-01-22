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
    set(CMAKE_Fortran_FLAGS_DEBUG " -g -cpp -Og -Wrealloc-lhs -fbounds-check -frange-check -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fcheck=all -fbacktrace -Wshadow=local -Wpedantic -Wextra -fallow-argument-mismatch")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    message(STATUS "Using ifort compiler flags")
    set(CMAKE_Fortran_FLAGS_RELEASE " -O3")
    set(CMAKE_Fortran_FLAGS_DEBUG " -g -traceback -check all -fpe0 -warn -debug extended -assume realloc_lhs")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    message(STATUS "Using IntelLLVM compiler flags")
    set(CMAKE_Fortran_FLAGS_RELEASE " -O3")
    set(CMAKE_Fortran_FLAGS_DEBUG " -g -O0 -traceback -check bounds")
#    set(CMAKE_Fortran_FLAGS_DEBUG " -g -traceback -check all -fpe0 -warn -debug extended -assume realloc_lhs")
#    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fpp -check bounds,uninit,pointers,stack -traceback -assume realloc_lhs -ftrapuv -debug minimal -check arg_temp_created")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp")
else()
    message(WARNING "Unrecognized Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}. Compiler-specific flags will not be set.")
endif()

# Locate NetCDF
find_path(NetCDF_INCLUDE_DIR NAMES netcdf.mod PATHS ${NetCDF_INCLUDE_DIR})
find_library(NetCDF_Fortran_LIBRARY NAMES netcdff PATHS ${NetCDF_LIBRARY_DIR})
find_library(NetCDF_C_LIBRARY NAMES netcdf PATHS ${NetCDF_LIBRARY_DIR})

# Check if both libraries were found
if(NOT NetCDF_Fortran_LIBRARY OR NOT NetCDF_C_LIBRARY)
    message(FATAL_ERROR "NetCDF Fortran or C library not found. Please specify the correct paths.")
endif()

if(NOT NetCDF_INCLUDE_DIR)
    message(FATAL_ERROR "NetCDF include directory. Please specify the correct paths.")
endif()

# Add the include directories
include_directories(${NetCDF_INCLUDE_DIR})

# Add the src subdirectory
add_subdirectory(src)


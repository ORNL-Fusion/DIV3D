#!/bin/bash

# This script sets up a CMake build with machine-specific settings.
#  Machines are identified based on uname -n
#
# Usage:
#   ./setup_cmake.sh [BUILD_TYPE] [COMPILER]
#
# Arguments (optional, in any order):
#   BUILD_TYPE  - Debug or Release (default: Release)
#   COMPILER    - GNU or IntelLLVM (default: GNU)
#
# Examples:
#   ./setup_cmake.sh Debug INTELLLVM
#   ./setup_cmake.sh Release GNU
#   ./setup_cmake.sh GNU Debug
#
# If no arguments are provided, the default is Release with the GNU compiler.
#
# Other flags:
#  NETCDF_INCLUDE_PATH: Location of netcdf.mod
#  NETCDF_LIB_PATH: Location of libnetcdf and libnetcdff
#  FORTRAN_COMPILER: Name (and location) of MPI compiler (e.g., mpif90)
#  USE_MPIF08: Use modern (Fortran 2008) MPI module, if available
#
#  Other necessary libraries (assumed to be in standard locations)
#    z, stdc++, lapack
# 
#  Corresponding cleanup script: clean_cmake.sh


rm -rf CMakeFiles CMakeCache.txt

# Default values
BUILD_TYPE=Release
VERBOSE_BUILD=1
COMPILER=GNU
USE_MPIF08=0

FORTRAN_COMPILER=""
NETCDF_INCLUDE_PATH=""
NETCDF_LIB_PATH=""
NETCDF_ROOT_DIR=""

# Parse optional arguments
for ARG in "$@"; do
    ARG=$(echo "$ARG" | tr '[:lower:]' '[:upper:]')  # Convert to uppercase
    if [[ "$ARG" =~ ^(DEBUG|RELEASE)$ ]]; then
        BUILD_TYPE=$ARG
    elif [[ "$ARG" =~ ^(GNU|INTELLLVM)$ ]]; then
        COMPILER=$ARG
    else
        echo "Unknown argument: $ARG"
        echo "Usage: ./setup_cmake.sh [Debug|Release] [GNU|INTELLLVM]"
        exit 1
    fi
done


# Identify machine
MACHINE_ID=$(uname -n)
echo "Building for machine $MACHINE_ID with $COMPILER compiler and $BUILD_TYPE build type."

# Set machine-specific paths
if [[ "$MACHINE_ID" == "ultrabucky" || "$MACHINE_ID" == "fusion3" ]]; then

    USE_MPIF08=1
    
    if [[ "$COMPILER" == "GNU" ]]; then
	NETCDF_INCLUDE_PATH="/usr/include"
	NETCDF_LIB_PATH="/usr/lib/x86_64-linux-gnu"	
        FORTRAN_COMPILER=mpif90
    elif [[ "$COMPILER" == "INTELLLVM" ]]; then
	NETCDF_ROOT_DIR="/home/jjl/intel/netcdf"
	NETCDF_INCLUDE_PATH="$NETCDF_ROOT_DIR/include"
	NETCDF_LIB_PATH="$NETCDF_ROOT_DIR/lib"   
        FORTRAN_COMPILER=mpiifx
    fi

elif [[ "$MACHINE_ID" == "mac145666" ]]; then

    USE_MPIF08=1
    NETCDF_INCLUDE_PATH="/opt/local/include"
    NETCDF_LIB_PATH="/opt/local/lib"
    FORTRAN_COMPILER=mpif90
    
elif [[ "$MACHINE_ID" == "THEALTANG23" ]]; then
    NETCDF_INCLUDE_PATH="/path/to/netcdf/include"
    NETCDF_LIB_PATH="/path/to/netcdf/lib64"
    FORTRAN_COMPILER=mpif90

elif [[ "$MACHINE_ID" == "stellar-intel.princeton.edu" ]]; then
    NETCDF_INCLUDE_PATH=${NETCDFDIR}/include
    NETCDF_LIB_PATH=${NETCDFDIR}/lib64
    FORTRAN_COMPILER=mpiifx

else
    echo "$MACHINE_ID is not supported by this script."
    echo "Please add your machine."
    exit 1
fi

# Print configuration
echo "CMake configuration:"
echo "  Build type: $BUILD_TYPE"
echo "  Compiler: $FORTRAN_COMPILER"


# Run cmake
cmake -DCMAKE_Fortran_COMPILER=$FORTRAN_COMPILER \
      -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
      -DNetCDF_INCLUDE_DIR="$NETCDF_INCLUDE_PATH" \
      -DNetCDF_LIBRARY_DIR="$NETCDF_LIB_PATH" \
      -DCMAKE_PREFIX_PATH="$NETCDF_ROOT_DIR" \
      -DUSE_MPIF08="$USE_MPIF08" \
      ..

# Build the project
if [ $VERBOSE_BUILD -eq 1 ]; then
    make VERBOSE=1
else
    make
fi

make

#!/bin/bash

# This is the script to be modified to set machine-specific settings.
#
# Instructions:
#  Add your machine below, and set the variables appropriately. Machine identified using uname -n.
#
#  BUILD_TYPE: Release or Debug. This can also by controlled by adding an argument to the
#              call to this script. "./setup_cmake.sh debug"
#
#  VERBOSE_BUILD: "on" by default, so that compiler commands are displayed.
#
#  NETCDF_INCLUDE_PATH: Location of netcdf.mod
#  NETCDF_LIB_PATH: Location of libnetcdf and libnetcdff
#
#  Other necessary libraries (assumed to be in standard locations)
#    z, stdc++, lapack
# 
#  Corresponding cleanup script: clean_cmake.sh

rm -rf CMakeFiles CMakeCache.txt

# Default values
BUILD_TYPE=Release
VERBOSE_BUILD=1
NETCDF_INCLUDE_PATH=""
NETCDF_LIB_PATH=""

# If argument passed, build debug
if [ $# -eq 1 ]; then
    BUILD_TYPE=Debug
fi

# Identify machine
MACHINE_ID=$(uname -n)

echo "Building for machine $MACHINE_ID"
echo "cmake configured to generate a $BUILD_TYPE build."

# Set machine-specific paths
if [ "$MACHINE_ID" == "ultrabucky" ] || [ "$MACHINE_ID" == "fusion3" ] || [ "$MACHINE_ID" == "mac145666" ]; then
    NETCDF_INCLUDE_PATH="/usr/include"
    NETCDF_LIB_PATH="/usr/lib/x86_64-linux-gnu"
elif [ "$MACHINE_ID" == "THEALTANG23" ]; then
    NETCDF_INCLUDE_PATH="/path/to/netcdf/include"
    NETCDF_LIB_PATH="/path/to/netcdf/lib64"
elif [ "$MACHINE_ID" == "stellar-intel.princeton.edu" ]; then
    NETCDF_INCLUDE_PATH=${NETCDFDIR}/include
    NETCDF_LIB_PATH=${NETCDFDIR}/lib64
else
    echo "$MACHINE_ID is not supported by this script."
    echo "Please add your machine."
    exit 1
fi

# Run cmake
cmake -DCMAKE_Fortran_COMPILER=mpif90 \
      -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
      -DNetCDF_INCLUDE_DIR="$NETCDF_INCLUDE_PATH" \
      -DNetCDF_LIBRARY_DIR="$NETCDF_LIB_PATH" \
      ..

# Build the project
if [ $VERBOSE_BUILD -eq 1 ]; then
    make VERBOSE=1
else
    make
fi

make install

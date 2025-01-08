#!/bin/bash

BUILD_TYPE=Release
VERBOSE_BUILD=1

MACHINE_ID=$(uname -n)

echo "Building for machine $MACHINE_ID"
echo

if [ $# -eq 1 ]; then
    BUILD_TYPE=Debug
fi

echo "cmake configured to generate a $BUILD_TYPE build."

rm -rf CMakeFiles CMakeCache.txt

NETCDF_INCLUDE_PATH=""
NETCDF_LIB_PATH=""
HDF5_INCLUDE_PATH=""
HDF5_LIB_PATH=""

# Set machine-specific paths
if [ "$MACHINE_ID" == "ultrabucky" ] || [ "$MACHINE_ID" == "fusion3" ] || [ "$MACHINE_ID" == "mac145666" ]; then
    NETCDF_INCLUDE_PATH="/usr/include"
    NETCDF_LIB_PATH="/usr/lib/x86_64-linux-gnu"
    HDF5_INCLUDE_PATH="/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include"
    HDF5_LIB_PATH="/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib"
elif [ "$MACHINE_ID" == "THEALTANG23" ]; then
    NETCDF_INCLUDE_PATH="/path/to/netcdf/include"
    NETCDF_LIB_PATH="/path/to/netcdf/lib64"
    HDF5_INCLUDE_PATH="/path/to/hdf5/include"
    HDF5_LIB_PATH="/path/to/hdf5/lib"
elif [ "$MACHINE_ID" == "stellar-intel.princeton.edu" ]; then
    NETCDF_INCLUDE_PATH=${NETCDFDIR}/include
    NETCDF_LIB_PATH=${NETCDFDIR}/lib64
    HDF5_INCLUDE_PATH="/usr/local/hdf5/oneapi-2024.2/1.14.4//include"
    HDF5_LIB_PATH="/usr/local/hdf5/oneapi-2024.2/1.14.4/lib64" 
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
      -DHDF5_INCLUDE_DIR="$HDF5_INCLUDE_PATH" \
      -DHDF5_LIBRARY_DIR="$HDF5_LIB_PATH" \
      ..

# Build the project
if [ $VERBOSE_BUILD -eq 1 ]; then
    make VERBOSE=1
else
    make
fi

make install

#!/bin/bash
#
# setup_cmake chooses the cmake build command for a machine
#

BUILD_TYPE=Release
VERBOSE_BUILD=1

MACHINE_ID=`uname -n`

echo "Building for machine $MACHINE_ID"
echo

if [ $# -eq 1 ]
then
    BUILD_TYPE=Debug
fi

echo "cmake configured to generate a $BUILD_TYPE build."
if [ "$BUILD_TYPE" == "Debug" ]
then
    echo "    cmake may be reconfigured to generate a Release build by running this script with no arguments or using:"
    echo
    echo "    cmake -DCMAKE_BUILD_TYPE=Release"
else
    echo "    cmake may be reconfigured to generate a Debug build by running this script with a Debug argument or using:"
    echo
    echo "    cmake -DCMAKE_BUILD_TYPE=Debug"
fi

rm -rf CMakeFiles CMakeCache.txt

echo

HDF5_ROOT_PATH=""
NETCDF_ROOT_PATH=""

# Set machine-specific defaults for HDF5 and NetCDF
if [ "$MACHINE_ID" == "ultrabucky" ] || [ "$MACHINE_ID" == "fusion3" ] || [ "$MACHINE_ID" == "mac145666" ]; then
    HDF5_ROOT_PATH="/usr/lib/x86_64-linux-gnu/hdf5/openmpi"
    NETCDF_ROOT_PATH="/usr/lib/x86_64-linux-gnu"
elif [ "$MACHINE_ID" == "THEALTANG23" ]; then
    # Example of a different machine with different paths:
    HDF5_ROOT_PATH="/path/to/hdf5"
    NETCDF_ROOT_PATH="/path/to/netcdf"
else
    echo "$MACHINE_ID is not supported by this script."
    echo "Please add your machine."
    echo
    exit 1
fi

cmake -DCMAKE_Fortran_COMPILER=mpif90 \
      -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE \
      -DHDF5_ROOT="$HDF5_ROOT_PATH" \
      -DNetCDF_ROOT="$NETCDF_ROOT_PATH" \
      ..

if [ $VERBOSE_BUILD -eq 1 ]
then
    make VERBOSE=$VERBOSE_BUILD
else
    make
fi
make install

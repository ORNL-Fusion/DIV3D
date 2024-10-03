#!/bin/bash
#
# setup_cmake chooses the cmake build command for a machine
#
#

BUILD_TYPE=Release
VERBOSE_BUILD=1

MACHINE_ID=`uname -n`

echo Building for machine $MACHINE_ID
echo

if [ $# -eq 1 ]
then
    BUILD_TYPE=Debug
fi
echo cmake configured to generate a $BUILD_TYPE build.
if [ "$BUILD_TYPE" == Debug ]
then
	echo "    cmake may be reconfigured to generate a Release build by running this script with no arguments or using the commmand"	
	echo
	echo "    cmake -DCMAKE_BUILD_TYPE=Release"
else
	echo "    cmake may be reconfigured to generate a Debug build by running this script with a Debug argument or using the commmand"
	echo
	echo "    cmake -DCMAKE_BUILD_TYPE=Debug"
fi

rm -rf CMakeFiles CMakeCache.txt

echo 
if [ "$MACHINE_ID" == "ultrabucky" ] || [ "$MACHINE_ID" == "fusion3" ]; then
    cmake -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE ..
elif [ "$MACHINE_ID" == "THEALTANG23" ]; then
    cmake -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE ..
else
    echo $MACHINE_ID is not supported by this script.
    echo Please add your machine.
    echo
fi

if [ $VERBOSE_BUILD -eq 1 ]
then
    make VERBOSE=$VERBOSE_BUILD
else
    make
fi
make install

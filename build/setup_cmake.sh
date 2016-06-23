#!/bin/bash
#
# setup_cmake chooses the cmake build command for a machine
#
#

BUILD_TYPE=Release
VERBOSE_BUILD=0
HAVE_BJDL=1

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
if [ "$MACHINE_ID" == "megabucky" ]
then
   # megabucky 
   # module gcc should be loaded
   # cmake by default will take old gcc from /usr/bin so define manually
    cmake -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DHAVE_BJDL=$HAVE_BJDL ..
elif [ "$MACHINE_ID" == "fusion2.ornl.gov" ]
then
    # fusion2
    #
    LIBBDIR=/home/jjl/util-library/fortran/bfield_library_jdl
    cmake -DCMAKE_Fortran_COMPILER=mpif90 -DHAVE_BJDL=$HAVE_BJDL -DCMAKE_BUILD_TYPE:String=$BUILD_TYPE -DLIBBJDL_DIR:String=$LIBBDIR ..
else
    echo $MACHINE_ID is not supported by this script.
    echo Please add your machine.
    echo
fi
echo bye

if [ $VERBOSE_BUILD -eq 1 ]
then
    make VERBOSE=$VERBOSE_BUILD
else
    make
fi
make install

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

detect_fortran_compiler() {
    for compiler in mpif90 mpifort mpiifx; do
        if command -v "$compiler" >/dev/null 2>&1; then
            echo "$compiler"
            return 0
        fi
    done
    return 1
}

get_mpi_backend_compiler() {
    local mpi_compiler="$1"
    local show_output

    if ! command -v "$mpi_compiler" >/dev/null 2>&1; then
        return 1
    fi

    show_output=$("$mpi_compiler" -show 2>/dev/null || true)
    if [[ -z "$show_output" ]]; then
        show_output=$("$mpi_compiler" --showme:command 2>/dev/null || true)
    fi

    if [[ -n "$show_output" ]]; then
        echo "$show_output" | awk '{print $1}'
        return 0
    fi

    return 1
}

validate_mpi_fortran_compiler() {
    local mpi_compiler="$1"
    local backend_compiler

    if [[ -z "$mpi_compiler" ]]; then
        echo "No MPI Fortran compiler was specified."
        return 1
    fi

    if ! command -v "$mpi_compiler" >/dev/null 2>&1; then
        echo "MPI Fortran compiler '$mpi_compiler' was not found in PATH."
        return 1
    fi

    backend_compiler=$(get_mpi_backend_compiler "$mpi_compiler" || true)
    if [[ -z "$backend_compiler" ]]; then
        echo "Unable to determine the underlying Fortran compiler used by '$mpi_compiler'."
        return 1
    fi

    if ! command -v "$backend_compiler" >/dev/null 2>&1; then
        echo "MPI Fortran compiler '$mpi_compiler' points to '$backend_compiler', but '$backend_compiler' was not found in PATH."
        return 1
    fi

    return 0
}

check_mpi_fortran_link() {
    local mpi_compiler="$1"
    local test_dir
    local test_src
    local test_exe

    test_dir=$(mktemp -d 2>/dev/null || true)
    if [[ -z "$test_dir" || ! -d "$test_dir" ]]; then
        echo "Unable to create a temporary directory for the MPI compiler test."
        return 1
    fi

    test_src="$test_dir/test_mpi_compiler.f90"
    test_exe="$test_dir/test_mpi_compiler.exe"

    cat > "$test_src" <<'EOF'
program test_mpi_compiler
  print *, 'MPI compiler test'
end program test_mpi_compiler
EOF

    if ! "$mpi_compiler" "$test_src" -o "$test_exe" >/dev/null 2>&1; then
        rm -rf "$test_dir"
        return 1
    fi

    rm -rf "$test_dir"
    return 0
}

detect_mpi_launcher() {
    for launcher in mpirun mpiexec; do
        if command -v "$launcher" >/dev/null 2>&1; then
            echo "$launcher"
            return 0
        fi
    done
    return 1
}

detect_netcdf_include_path() {
    for path in /usr/include /usr/local/include /usr/include/netcdf /opt/local/include; do
        if [[ -f "$path/netcdf.mod" ]]; then
            echo "$path"
            return 0
        fi
    done
    return 1
}

detect_netcdf_lib_path() {
    for path in /usr/lib/x86_64-linux-gnu /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib /opt/local/lib; do
        if [[ -f "$path/libnetcdff.so" || -f "$path/libnetcdff.a" || -f "$path/libnetcdf.so" || -f "$path/libnetcdf.a" ]]; then
            echo "$path"
            return 0
        fi
    done
    return 1
}

# Default values
BUILD_TYPE=Release
VERBOSE_BUILD=1
COMPILER=${COMPILER:-GNU}
USE_MPIF08=${USE_MPIF08:-0}

FORTRAN_COMPILER=${FORTRAN_COMPILER:-""}
NETCDF_INCLUDE_PATH=${NETCDF_INCLUDE_PATH:-""}
NETCDF_LIB_PATH=${NETCDF_LIB_PATH:-""}
NETCDF_ROOT_DIR=${NETCDF_ROOT_DIR:-""}

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
    DETECTED_FORTRAN_COMPILER=$(detect_fortran_compiler || true)
DETECTED_NETCDF_INCLUDE_PATH=$(detect_netcdf_include_path || true)
DETECTED_NETCDF_LIB_PATH=$(detect_netcdf_lib_path || true)
DETECTED_MPI_LAUNCHER=$(detect_mpi_launcher || true)

    if [[ -z "$FORTRAN_COMPILER" && -n "$DETECTED_FORTRAN_COMPILER" ]]; then
        FORTRAN_COMPILER="$DETECTED_FORTRAN_COMPILER"
    fi
    if [[ -z "$NETCDF_INCLUDE_PATH" && -n "$DETECTED_NETCDF_INCLUDE_PATH" ]]; then
        NETCDF_INCLUDE_PATH="$DETECTED_NETCDF_INCLUDE_PATH"
    fi
    if [[ -z "$NETCDF_LIB_PATH" && -n "$DETECTED_NETCDF_LIB_PATH" ]]; then
        NETCDF_LIB_PATH="$DETECTED_NETCDF_LIB_PATH"
    fi

    echo
    echo "Machine $MACHINE_ID is not recognized by this setup script."
    echo
    echo "Possible settings detected on this system:"
    echo "  FORTRAN_COMPILER=${FORTRAN_COMPILER:-not found}"
    echo "  MPI_LAUNCHER=${DETECTED_MPI_LAUNCHER:-not found}"
    echo "  NETCDF_INCLUDE_PATH=${NETCDF_INCLUDE_PATH:-not found}"
    echo "  NETCDF_LIB_PATH=${NETCDF_LIB_PATH:-not found}"
    echo
    echo "If these look right, you can try:"
    echo "  FORTRAN_COMPILER=${FORTRAN_COMPILER:-mpif90} NETCDF_INCLUDE_PATH=${NETCDF_INCLUDE_PATH:-/path/to/netcdf/include} NETCDF_LIB_PATH=${NETCDF_LIB_PATH:-/path/to/netcdf/lib} ./setup_cmake.sh Release $COMPILER"
    echo
    echo "Environment variables override the machine-specific defaults in this script."

    if [[ -z "$FORTRAN_COMPILER" || -z "$NETCDF_INCLUDE_PATH" || -z "$NETCDF_LIB_PATH" ]]; then
        echo
        echo "Unable to determine a complete set of build paths automatically."
        echo "Please load the needed modules and set the variables above, or add this machine to the script."
        exit 1
    fi
fi

# Print configuration
echo "CMake configuration:"
echo "  Build type: $BUILD_TYPE"
echo "  Compiler: $FORTRAN_COMPILER"

if ! validate_mpi_fortran_compiler "$FORTRAN_COMPILER"; then
    echo
    echo "MPI Fortran compiler check failed."
    echo "Please load the MPI and Fortran compiler modules, or set FORTRAN_COMPILER to a working MPI compiler wrapper."
    exit 1
fi

if ! check_mpi_fortran_link "$FORTRAN_COMPILER"; then
    echo
    echo "MPI Fortran compiler '$FORTRAN_COMPILER' could not compile and link a simple test program."
    echo "Please load a complete MPI environment before running setup_cmake.sh."
    exit 1
fi

if [[ -z "$DETECTED_MPI_LAUNCHER" ]]; then
    echo
    echo "No MPI launcher was found in PATH."
    echo "Please load the MPI environment so that mpirun or mpiexec is available."
    exit 1
fi

BACKEND_FORTRAN_COMPILER=$(get_mpi_backend_compiler "$FORTRAN_COMPILER" || true)
echo "  MPI backend compiler: ${BACKEND_FORTRAN_COMPILER:-unknown}"
echo "  MPI launcher: $DETECTED_MPI_LAUNCHER"


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

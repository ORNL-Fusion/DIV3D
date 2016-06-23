#-----------------------------------------------------------------------
# This is the makefile for DIV3D
#
#  - Intended to be run from the 'build' directory
#
# targets:
#  [default] -  Makes executable
#   debug    -  Makes exec with debug options
#   clean    -  Removes *.mod, *.o from build directory
#   cleaner  -  Same as clean AND removes  executables from bin directory
#   
#
#	7/14/2011 - current JDL
#-----------------------------------------------------------------------

SRC_DIR := ../src/
BIN_DIR := ../bin/
W7X_DIR := $(SRC_DIR)W7X_routines/
FXDR_DIR:= $(SRC_DIR)fxdr/

# Location of bfield library and .mod files
LIBBJDL_DIR = /home/jjl/util-library/fortran/bfield_library_jdl/lib
LIBBJDL_INCDIR = /home/jjl/util-library/fortran/bfield_library_jdl/include

# Location of m3dc1 libraries
#LIBM3DC1_DIR = /home/jjl/m3dc1/install/lib/_X86_64 

FC := mpif90
FF := -O3 
FF_DEBUG := -g -Wall -frange-check -fbounds-check 

CC = gcc

EXEC := $(BIN_DIR)div3d_follow_and_int
EXEC_DBG := $(BIN_DIR)dbg_div3d_follow_and_int

#-----------------------------------------------------------------------

# Libraries needed to build the executable
FXDRLIB = libfxdr.a
#LIBDIR := -L$(LIBBJDL_DIR) -L$(LIBM3DC1_DIR) -L.
#FL := -lfxdr -lbjdl -lfusionio -lm3dc1 -lz -lhdf5 -lstdc++
LIBDIR := -L$(LIBBJDL_DIR) -L.
FL := -lfxdr -lbjdl -llapack
LIBS := $(FXDRLIB) 

# Objects needed to build the executable
OBJS90 :=  \
follow_fieldlines_rzphi.o \
inside_vessel_mod.o  \
kind_mod.o \
math_routines_mod.o \
io_unit_spec.o \
div3d_follow_and_int.o \
init_points_line.o \
trace_surface.o \
init_random_seed.o \
read_parts_mod.o \
diffuse_lines3.o \
parallel_mod.o \
integrator_routines_mod.o \
fieldline_following_mod.o \
find_nearby_triangles.o 

OBJS90_W7X := bfield_xdr.o
OBJS90_FXDR := fxdr_module.o

# Modules needed to build the executable
MODS90 := \
kind_mod.mod \
io_unit_spec.mod \
read_parts_mod.mod \
parallel_mod.mod \
integrator_routines_mod.mod \
fieldline_following_mod.mod \
math_routines_mod.mod 

MODS90_W7X := bfield_xdr.mod 
MODS90_FXDR := fxdr_module.mod

# list of c object files for fxdr library
OBJS_C_FXDR := cxdrinit.o cxdrint.o cxdrclose.o cxdrreal.o cxdrdouble.o \
	       cxdrrmat.o cxdrdmat.o cxdrimat.o cxdrstring.o		\
	       cxdrreal64.o cxdrrewind.o cxdrrmat64.o cxdrshort.o

#-----------------------------------------------------------------------

$(EXEC): $(OBJS90) $(OBJS90_W7X) $(OBJS90_FXDR) $(MODS90) $(MODS90_W7X) $(MODS90_FXDR) $(LIBS)
	@echo "Creating "$@" in directory "$(PWD)"."
	${FC} -o $@ $(FF) $(OBJS90) $(OBJS90_W7X) $(OBJS90_FXDR) $(LIBDIR) $(FL)

.PHONY.: debug
debug: FF=$(FF_DEBUG)
debug: $(EXEC)

#-----------------------------------------------------------------------
# Static patterns

$(OBJS_C_FXDR) : %.o : $(FXDR_DIR)%.c
	$(CC) $(CF) -c $<

$(OBJS90) : %.o:  $(SRC_DIR)%.f90
	$(FC) $(FF) -c -I$(LIBBJDL_INCDIR) $(LIBDIR) $(FL) $<

$(OBJS90_W7X) : %.o:  $(W7X_DIR)%.f90
	$(FC) $(FF) -c $<

$(OBJS90_FXDR) : %.o:  $(FXDR_DIR)%.f90
	$(FC) $(FF) -c $<

$(MODS90) : %.mod:  $(SRC_DIR)%.f90
	$(FC) $(FF) -c -I$(LIBBJDL_INCDIR) $(LIBDIR) $(FL) $<

$(MODS90_W7X) : %.mod:  $(W7X_DIR)%.f90
	$(FC) $(FF) -c $<

$(MODS90_FXDR) : %.mod:  $(FXDR_DIR)%.f90
	$(FC) $(FF) -c $<

$(FXDRLIB): $(OBJS_C_FXDR)
	@echo "Creating "$(FXDRLIB)" in directory "$(PWD)"."
	$(AR) rc $(FXDRLIB) $(OBJS_C_FXDR)
#	-$(RANLIB) $(FXDRLIB)


# dependencies
div3d_follow_and_int.o : kind_mod.mod io_unit_spec.mod read_parts_mod.mod parallel_mod.mod
follow_fieldlines_rzphi.o: kind_mod.mod bfield_xdr.mod
parallel_mod.o : kind_mod.mod
math_routines_mod.mod: kind_mod.mod
inside_vessel_mod.mod: kind_mod.mod
read_parts_mod.mod: kind_mod.mod io_unit_spec.mod  
integrator_routines_mod.o: kind_mod.mod
integrator_routines_mod.mod: kind_mod.mod
diffuse_lines3.o: kind_mod.mod integrator_routines_mod.mod math_routines_mod.mod \
fieldline_following_mod.mod
fieldline_following_mod.o: kind_mod.mod bfield_xdr.mod
fieldline_following_mod.mod: kind_mod.mod
find_nearby_triangles.o: read_parts_mod.mod math_routines_mod.mod
bfield_xdr.o: kind_mod.mod fxdr_module.mod
bfield_xdr.mod: fxdr_module.mod

clean:
	-rm *.o *.mod *.a
	@echo "Cleaned object, mod and library files from "$(PWD)

cleaner:
	-rm *.o *.mod *.a $(EXEC) $(EXEC_DBG)
	@echo "Cleaned object, mod, library and executable files from "$(PWD)

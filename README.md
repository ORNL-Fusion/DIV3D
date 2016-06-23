# DIV3D
DIV3D field line following and intersection code.
Written 2010 - current JDL

#
# COMPILATION
#
0) Build process uses cmake.  GCC compilers have been tested. LAPACK required.  MPI required.
1) In directory build edit setup_cmake.sh.  Add your system to have the proper
   flags sent to cmake calls, or call "cmake .." with appropriate options from
   the build directory. Important options are HAVE_BJDL, set to 1 if libbjdl is
   available.  Note this is only required for tokamak cases that have bfield files
   read from EFIT equdsk files. Other options should be self-explanatory, but note
   that the fortran compiler/linker must be able to compile with MPI, otherwise
   link libraries manually.
2) Check that compiler flags set in CmakeLists.txt are appropriate (should be fine for GCC).
   Note that if linking fails you may need to add system-dependent libraries in the
   target_link_libraries call (where lapack is listed).
3) enter build and run "./setup_cmake.sh [debug]"
   After initial setup you can just run "make install" unless you are switching between
   debug and release builds.

#
# SETUP AND RUNNING
#
1) See test directory
2) Input files: parts.list, run_settings.nml

   parts.list just contains locations of geometry files, which can be either type "part" or type "jpart".
   ==========================================
   See test/geo_files for examples and info.

   run_settings.nml is the control namelist.
   ===========================================
   R,Z,phistart: Initial position of field line used to trace closed surface with dphi_line_surf_deg and ntran_surf steps. 

   npts_start lines are followed in both directions with dphi_line_diff_deg and ntran_diff steps, with dmag [m**2/m] perpendicular diffusion.

   xdr_check: Set to true to perform the check of xdr file for "bad points".  Will then process and make new file. TURN OFF FOR VACUUM CASES! Ask J. Geiger for more details.

   working_master: set to true to have main process do work instead of just control.

   trace_surface_opt: Set to false to skip generation of closed surface if restarting run

   hit_length: Length of line segment saved for post-processing [m].  (Used for scraper element design).

   lsfi_tol: Tolerance on checking for intersection of line segments with facets.


This product includes software developed by the
University of California, San Diego and its contributor
# DIV3D
DIV3D field line following and intersection code.
Written 2010 - current JDL

#
# COMPILATION
#
0) Build process uses cmake.  GCC compilers have been tested. LAPACK required if linking to libbjdl.  MPI required.
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


#
# OUTPUT
#

%%%%%%%%%%%%%%
launch_pts.out
%%%%%%%%%%%%%%
-- Contains list of points from which field lines are initialized
Write(Int*1) : NUM_PTS 
for i = 1:NUM_PTS
    Write(Real*3) : R  Z  Phi
end

%%%%%%%%%%%%%%
hitcount.out
%%%%%%%%%%%%%%
-- Contains number of lines that intersected the listed components
Write(Int*1) : NHIT

%%%%%%%%%%%%%%
allparts.out
%%%%%%%%%%%%%%
-- Contains component geometry
-- NTOR_MAX and NPOL_MAX are used to initialize part arrays to the same size
Write(Int*3) : NPARTS NTOR_MAX NPOL_MAX    ! Number of components, max number of toroidal points, max number of poloidal points
for ipart = 1:nparts
    Write(Int*2) : NTOR(ipart) NPOL(ipart)
    for i = 1:ntor_max
        for j = 1:npol_max
	    Write(Real*3) : R(ipart,i,j)  Z(ipart,i,j)  Phi(ipart,i,j)
        end
    end
end

%%%%%%%%%%%%%%%%%%%
part_triangles.out
%%%%%%%%%%%%%%%%%%%
-- Contains triangle verticies
-- Each segement of each part (toroidal and poloidal) is broken into two triangles
Write(Int*1) : NPARTS
for i = 1:nparts
    ITRI = 1
    Write(Int*2) : IPART  NTRI_PARTS(ipart)
    for itor = 1:ntor-1
      for jpol = 1:npol-1
        Write(Int*1) : ITRI
	Write(Real*3) : XTRI(IPART,ITRI,1)  YTRI(IPART,ITRI,1)  ZTRI(IPART,ITRI,1)  ! 1st triangle point 1
        Write(Real*3) : XTRI(IPART,ITRI,2)  YTRI(IPART,ITRI,2)  ZTRI(IPART,ITRI,2)  ! 1st triangle point 2
        Write(Real*3) : XTRI(IPART,ITRI,3)  YTRI(IPART,ITRI,3)  ZTRI(IPART,ITRI,3)  ! 1st triangle point 3
	ITRI = ITRI + 1
        Write(Int*1) : ITRI
	Write(Real*3) : XTRI(IPART,ITRI,1)  YTRI(IPART,ITRI,1)  ZTRI(IPART,ITRI,1)  ! 2nd triangle point 1
        Write(Real*3) : XTRI(IPART,ITRI,2)  YTRI(IPART,ITRI,2)  ZTRI(IPART,ITRI,2)  ! 2nd triangle point 2
        Write(Real*3) : XTRI(IPART,ITRI,3)  YTRI(IPART,ITRI,3)  ZTRI(IPART,ITRI,3)  ! 2nd triangle point 3
	ITRI = ITRI + 1
      end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
part_triangle_mids.out
%%%%%%%%%%%%%%%%%%%%%%% 
-- Contains triangle midpoints
Write(Int) : NPARTS
for i = 1:nparts
    ITRI = 1
    Write(Int) : IPART  NTRI_PARTS(ipart)
    for itor = 1:ntor-1
      for jpol = 1:npol-1
        Write(Int*1,Real*4) : ITRI  XMID(IPART,ITRI)  YMID(IPART,ITRI)  ZMID(IPART,ITRI)  DMID(IPART,ITRI)  ! 1st triangle: Index, x,y,z of midpoint, Max distance from midpoint to vertex
        ITRI = ITRI + 1
	Write(Int*1,Real*4) : ITRI  XMID(IPART,ITRI)  YMID(IPART,ITRI)  ZMID(IPART,ITRI)  DMID(IPART,ITRI)  ! 2nd triangle: Index, x,y,z of midpoint, Max distance from midpoint to vertex
      end
    end
end
    
%%%%%%%%%%%%%%
int_pts.out
%%%%%%%%%%%%%%
-- Contains intersection info
Column 1: R_hit
Column 2: Z_hit
Column 3: Phi_hit (radians)
Column 4: ihit               ! = 1 if it hit a part
Column 5: ipart              ! Index of the part hit (from allparts)
Column 6: itri               ! Index of the triangle intersected (from part_triangles)
Column 7: iline              ! Index of the line (from launch_points)


%%%%%%%%%%%%%%
hitline.out
%%%%%%%%%%%%%%
-- Contains diffused field line segment for each intersecting line along with intersection data
For i = 1,nhit
  Write(Real*3,Int*4) : Rhit  Zhit  Phihit  IOUT  ! Intersection point, IOUT = [IHIT,IPART,ITRI,I] (where I is index of fieldline integration)
  Write(Int*1) : nhitline  ! Number of points along line that will be output
  Write(Real*nhitline) : r_hitline       ! Line geometry
  Write(Real*nhitline) : z_hitline 
  Write(Real*nhitline) : phi_hitline
end

%%%%%%%%%%%%%%
surface_line.out
%%%%%%%%%%%%%%
-- Contains initializing surface data
Write(Real*1,Int*2) : Period, nip0, ip_step   ! Field line info
Write(Int*1) : Num_pts ! Number of points on field line
For i = 1:Num_pts
  Write(Real*1) : r(i)
  Write(Real*1) : z(i)
  Write(Real*1) : phi(i)
end
This product includes software developed by the
University of California, San Diego and its contributor
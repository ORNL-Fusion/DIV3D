!-----------------------------------------------------------------------------
!+ Program to calculate intersections of fieldlines initiated on a flux surface
!  and diffused until they intersect components
!-----------------------------------------------------------------------------
program div3d_follow_and_int
! Author(s): J.D. Lore - 07/14/2011 - xxx
!
!
! To do:
!  - Add default namelist variables or check for missing
!  - Vessel intersection is just nearest phi cut
!  - Need to add exiting routine to deallocate and finalize mpi
!
! Modules used:
Use run_settings_namelist
Use parallel_mod, Only : rank, nprocs, init_mpi, fin_mpi
Use read_parts_mod, Only : read_parts, make_triangles
Use setup_bfield_module
Use diffusion, Only : diffuse_lines3, diffuse_lines3_worker
Implicit none

Logical :: verbose = .false.

!----------------------------------------------------------
! 0. Setup
! --Initialize MPI
! --Read namelist
!----------------------------------------------------------
Call init_mpi
If (rank .eq. 0) verbose = .true. 
If (verbose) Write(*,'(/A)') '-------------------------------------------------------------------------'
setup_bfield_verbose = verbose

! Read namelists
Call read_run_settings_namelist(verbose)

!----------------------------------------------------------
! 1. Initialize magnetic field
!----------------------------------------------------------
! Setup rmp field
Select Case (rmp_type)
  Case ('g')
    Call setup_bfield_g3d
#ifdef HAVE_FXDR 
  Case ('xdr')
    Call setup_bfield_xdr
#endif 
  Case ('vmec_coils')
    Call setup_bfield_vmec_coils
  Case ('vmec_coils_to_fil')
    Call setup_bfield_vmec_coils_to_fil
  Case ('bgrid')
    Call setup_bfield_bgrid        
  Case Default
    If (rank == 0) Then
      Write(*,*) 'Unknown rmp_type in div3d!'
      Write(*,*) 'Current options are:'
      Write(*,*) '''g'''      
      Write(*,*) '''vmec_coils'''
      Write(*,*) '''vmec_coils_to_fil'''
      Write(*,*) '''bgrid'''
#ifdef HAVE_FXDR 
      Write(*,*) '''xdr'''
#endif
    Endif
    Stop      
End Select

!----------------------------------------------------------------
! 2. Load intersection components
!  -- read parts.list file for file names
!  -- read part files
!  -- make trianges from parts
!----------------------------------------------------------------
If (verbose) Write(*,'(A)') ' Reading parts list and part files:'
Call read_parts(fname_plist,fname_parts,fname_ves,verbose)
  
! Make triangles from 2d parts
If (verbose) Write(*,'(/A/)') ' Generating 2d part triangles'
Call make_triangles(fname_ptri,fname_ptri_mid)

!----------------------------------------------------------------
! Master node traces initial surface and defines starting points
! Master then calls diffuse_lines3 and distributes field lines
! to be followed
!----------------------------------------------------------------
If (rank .eq. 0) Then

  Write(*,'(A,I0,A)') ' Master node (rank ',rank,') is initializing run'
  Write(*,'(A,I0,A)') ' Total number of processes: ',nprocs

  If (npts_start .lt. nprocs) Then
    Write(*,*) 'nprocs must be less than npts_start!!!'
    Stop "Cannot handle this" 
  Endif

  !----------------------------------------------------------
  ! 3. Trace out a surface (no diffusion)
  !----------------------------------------------------------
  If (trace_surface_opt .eqv. .true.) Then
    Write(*,'(/A,3(F8.2))') ' Tracing initial surface from (R,Z,Phi) = ',Rstart,Zstart,Phistart
    If (.true.) Call trace_surface(Rstart,Zstart,Phistart,dphi_line_surf,ns_line_surf,period,fname_surf)

    !----------------------------------------------------------
    ! 4. Initialize points on surface to carry heat
    !----------------------------------------------------------
    Write(*,*) 'Initializing points along initial surface line'
    If (.true.) Call init_points_line(fname_surf,npts_start,fname_launch)
  Else
    Write(*,*) 'Skipping trace_surface, loading file.'
  Endif

  !----------------------------------------------------------------------------------------
  ! 5. Follow fieldlines from init points and check for intersections 
  !----------------------------------------------------------------------------------------

  Write(*,'(/A)') '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
  write(*,*) 'Beginning MPI line diffusion'
  Call diffuse_lines3(fname_launch,dmag,ns_line_diff,fname_hit,fname_intpts,fname_nhit,nhitline)

Else

 !----------------------------------------------------------------
 ! Additional nodes wait for lines and follow them
 !----------------------------------------------------------------
 Call diffuse_lines3_worker(dmag,dphi_line_diff,nhitline,period,calc_lc,calc_theta,lsfi_tol)
Endif 


! Finialize MPI
Call fin_mpi(.false.) ! False means this is a non-error exit

end program div3d_follow_and_int






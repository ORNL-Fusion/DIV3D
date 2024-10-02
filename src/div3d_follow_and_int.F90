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
Use kind_mod, Only : real64, int32
Use parallel_mod
Use io_unit_spec, Only: iu_nl 
Use read_parts_mod
Use setup_bfield_module
Use phys_const, Only : pi
Use diffusion, Only : diffuse_lines3, diffuse_lines3_worker
Use init_random, Only : init_random_seed
Implicit none

! Local scalars
Real(real64) :: Rstart, Zstart, Phistart
Real(real64) :: dmag, period, lsfi_tol
Real(real64) :: dphi_line_surf, dphi_line_diff
Real(real64) :: dphi_line_surf_deg, dphi_line_diff_deg, hit_length

Integer(int32) :: myseed, nhitline
Integer(int32) :: npts_start, nfp
Integer(int32) :: iocheck
Integer(int32) :: ntran_surf, ns_line_surf
Integer(int32) :: ntran_diff, ns_line_diff

Character(len=300) :: fname_hit, fname_ptri, fname_ptri_mid
Character(len=300) :: fname_launch,fname_surf, fname_parts, fname_intpts, fname_ves
Character(len=300) :: fname_plist, fname_nhit

Logical :: verbose, trace_surface_opt, calc_lc, calc_theta

! Namelists
Namelist / run_settings / fname_plist, fname_ves,  &
  fname_surf, fname_launch, fname_parts, fname_hit, fname_intpts, nfp, &
  Rstart, Zstart, Phistart, dphi_line_surf_deg, ntran_surf, &
  npts_start, dmag, dphi_line_diff_deg, ntran_diff, myseed, &
  fname_nhit, hit_length, lsfi_tol, trace_surface_opt, &
  fname_ptri, fname_ptri_mid, calc_lc, calc_theta

!- End of header -------------------------------------------------------------

!----------------------------------------------------------
! 0. Setup
! --Initialize MPI
! --Read namelist
! --Init random number generator
! --basic conversions
!----------------------------------------------------------
Call init_mpi
verbose = .false.
If (rank .eq. 0) verbose = .true. 
If (verbose) Write(6,'(/A)') '-------------------------------------------------------------------------'

! Defaults
calc_lc = .true.
calc_theta = .false.

! Read the run settings namelist file
setup_bfield_verbose = verbose
If (verbose) Write(*,*) 'Reading run settings from run_settings.nml'
Open(iu_nl,file="run_settings.nml",status="old",form="formatted",iostat=iocheck)
If ( iocheck .ne. 0 ) Then
  Write(6,*) 'Error opening namelist file'
  Stop 'Exiting: I/O Error in div3d_follow_and_int.f90'
Endif
Read(iu_nl,nml=run_settings)
Rewind(iu_nl)
Read(iu_nl,nml=bfield_nml)
Close(iu_nl)

If (verbose) Write(*,*) 'Initializing random number with base seed:',myseed
Call init_random_seed(myseed*(rank+1))

period = 2.d0*pi/Real(nfp,real64)
dphi_line_surf = dphi_line_surf_deg * pi/180.d0
dphi_line_diff = dphi_line_diff_deg * pi/180.d0
ns_line_surf = Floor(ntran_surf*2.d0*pi/Abs(dphi_line_surf))
ns_line_diff = Floor(ntran_diff*2.d0*pi/Abs(dphi_line_diff))

If (verbose) Write(*,'(/A,G0.3)') 'hit_length: ',hit_length
If (hit_length .le. 0.d0) Then
   If (verbose) Write(*,'(/A)') 'Turning off hitline because hit_length <= 0'
   nhitline = 0
Else
   nhitline = Floor(hit_length/Rstart/abs(dphi_line_diff))
   If (verbose) Write(*,'(A,I0,A)') ' Returning ',nhitline,' points along intersecting lines'
Endif

If (verbose) Then
   If (calc_lc) Then
      Write(*,*) 'Computing one-directional connection length'
   Else
      Write(*,*) 'Not computing one-directional connection length'   
   Endif
Endif

!----------------------------------------------------------
! 1. Initialize magnetic field
!----------------------------------------------------------
! Setup rmp field
if (verbose .AND. rank .EQ. 0) write(*,*) 'Bfield method is ',rmp_type
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
#ifdef HAVE_FXDR 
      Write(*,*) '''xdr'''
#endif
      Write(*,*) '''bgrid'''            
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






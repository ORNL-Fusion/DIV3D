!-----------------------------------------------------------------------------
!+ Program to calculate intersections of fieldlines initiated on a flux surface
!  and diffused until they intersect components
!-----------------------------------------------------------------------------
program div3d_follow_and_int
!
! Description: 
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/14/2011  
!  1.1     12/07/2011  -- Added 3D part reading and triangle intersection. JDL
!  1.2     12/14/2011  -- Added MPI. JDL
!  1.3     02/15/2012  -- Updated to allow reading of parts not defined in toroidal planes.
! Author(s): J.D. Lore - 07/14/2011 - xxx
!
!
! To do:
!  - Add default namelist variables or check for missing
!  - Vessel intersection file not wiped if no intersections?
!  - Vessel intersection is just nearest phi cut
!  - Check allocation/deallocation
!  - Need to add exiting routine to deallocate and finalize mpi
!
! Modules used:
Use kind_mod
Use parallel_mod
Use bfield_xdr, Only: &
! Imported subroutines
readbgrid
Use io_unit_spec, Only: &
iu_nl,    &  ! Run settings namelist file unit (run_settings.nml,input)
iu_plist, &  ! Parts filename list file (input) 
iu_parts, &  ! Parts data file (output)
iu_nhit,  & 
iu_surf,  &
iu_hit,   &
iu_int
Use read_parts_mod
#ifdef HAVE_BJDL
Use g3d_module, Only : readg_g3d
Use setup_bfield_module
#endif
Implicit none

! Local scalars
Real(rknd) :: Rstart, Zstart, Phistart
Real(rknd) :: dmag, period, lsfi_tol
Real(rknd) :: dphi_line_surf, dphi_line_diff
Real(rknd) :: dphi_line_surf_deg, dphi_line_diff_deg, hit_length

Integer(iknd) :: myseed, nhitline
Integer(iknd) :: npts_start, nfp
Integer(iknd) :: iocheck
Integer(iknd) :: ntran_surf, ns_line_surf
Integer(iknd) :: ntran_diff, ns_line_diff
Integer(iknd) :: my_numl, ifl, ierr_follow, iline, dest, num_myjobs, source, tag
Integer(iknd) :: div3d_bfield_method

Character(len=300) :: fname_hit, fname_bfile, fname_ptri, fname_ptri_mid
Character(len=300) :: fname_launch,fname_surf, fname_parts, fname_intpts, fname_ves
Character(len=300) :: fname_plist, fname_nhit

Logical :: verbose, trace_surface_opt

! Local arrays

Real(rknd), Dimension(3) :: pint
Real(rknd), Dimension(:), Allocatable :: r_hitline,z_hitline,phi_hitline

Real(rknd), Dimension(6) :: line_start_data_r
Integer(iknd), Dimension(3) :: line_start_data_i
Integer(iknd), Dimension(4) :: iout
Real(rknd), Dimension(:), Allocatable :: line_done_data_r
Integer(iknd), Dimension(5) :: line_done_data_i

! Local Parameters
Real(rknd), Parameter :: pi = 3.141592653589793238462643383279502_rknd

! Namelists
Namelist / run_settings / fname_plist, fname_ves,  &
  fname_surf, fname_launch, fname_parts, fname_hit, fname_intpts, nfp, &
  Rstart, Zstart, Phistart, dphi_line_surf_deg, ntran_surf, &
  npts_start, dmag, dphi_line_diff_deg, ntran_diff, myseed, &
  fname_nhit, hit_length, lsfi_tol, working_master, trace_surface_opt, &
  fname_ptri, fname_ptri_mid

!- End of header -------------------------------------------------------------

!----------------------------------------------------------
! 0. Setup
! --Initialize MPI
! --Read namelist
! --Init random number generator
! --basic conversions
!----------------------------------------------------------
Call init_mpi()
verbose = .false.
If (rank .eq. 0) verbose = .true. 
If (verbose) Write(6,'(/A)') '-------------------------------------------------------------------------'

! Read the run settings namelist file
If (verbose) Write(6,*) 'Reading run settings from run_settings.nml'
Open(iu_nl,file="run_settings.nml",status="old",form="formatted",iostat=iocheck)
If ( iocheck .ne. 0 ) Then
  Write(6,*) 'Error opening namelist file'
  Stop 'Exiting: I/O Error in div3d_follow_and_int.f90'
Endif
Read(iu_nl,nml=run_settings)
Rewind(iu_nl)
Read(iu_nl,nml=bfield_nml)
Close(iu_nl)

If (verbose) Write(6,*) 'Initializing random number with base seed:',myseed
Call init_random_seed(myseed*(rank+1))

period = 2.d0*pi/Real(nfp,rknd)
dphi_line_surf = dphi_line_surf_deg * pi/180.d0
dphi_line_diff = dphi_line_diff_deg * pi/180.d0
ns_line_surf = Floor(ntran_surf*2.d0*pi/Abs(dphi_line_surf))
ns_line_diff = Floor(ntran_diff*2.d0*pi/Abs(dphi_line_diff))

!----------------------------------------------------------
! 1. Initialize magnetic field
!  W7X: 
!     Read the bfield file
!     -- The loaded variables are passed via the bfield_xdr
!        module to subroutine bint (called during fieldline 
!        following)
!  NSTX: 
!     Load the gfile
!     --> RMP etc not implemented yet
!----------------------------------------------------------
! Setup rmp field
Select Case (rmp_type)
  Case ('g3d')
    Call setup_bfield_g3d
    div3d_bfield_method = 1
    if (verbose .AND. rank .EQ. 0) write(*,*) 'Bfield method is g3d'
  Case ('xdr')
    Call setup_bfield_xdr
    div3d_bfield_method = 0
    if (verbose .AND. rank .EQ. 0) write(*,*) 'Bfield method is xdr'
  Case Default
    If (rank == 0) Then
      Write(*,*) 'Unknown rmp_type in div3d!'
      Write(*,*) 'Current options are:'
      Write(*,*) '''g3d'''
    Endif
    Stop      
End Select

!----------------------------------------------------------------
! 2. Load intersection components
!  -- read parts.list file for file names
!  -- read part files
!  -- make trianges from parts
!----------------------------------------------------------------
If (verbose) Write(6,'(A)') ' Reading parts list and part files:'
Call read_parts(fname_plist,fname_parts,fname_ves,verbose)
  
! Make triangles from 2d parts
If (verbose) Write(6,'(/A/)') ' Generating 2d part triangles'
Call make_triangles(fname_ptri,fname_ptri_mid)

!----------------------------------------------------------------
! Master node traces initial surface and defines starting points
! Master then calls diffuse_lines3 and distributes field
!  lines to be followed
!----------------------------------------------------------------
If (rank .eq. 0) Then


  Write(6,'(A,I0,A)') ' Master node (rank ',rank,') is initializing run'
  Write(6,'(A,I0,A)') ' Total number of processes: ',nprocs

  If (npts_start .lt. nprocs) Then
    Write(*,*) 'For now nprocs must be less than npts_start!!!'
  Endif

  If ((nprocs .lt. 2) .AND. (working_master .eqv. .false.)) Then
    Write(6,*) ' Must use at least two processors if working_master is false!!!'
    Stop
  Endif

  !----------------------------------------------------------
  ! 3. Trace out a surface (no diffusion)
  !----------------------------------------------------------
  If (trace_surface_opt .eqv. .true.) Then
    Write(6,'(/A,3(F8.2))') ' Tracing initial surface from (R,Z,Phi) = ',Rstart,Zstart,Phistart
    If (.true.) Call trace_surface(Rstart,Zstart,Phistart,dphi_line_surf,ns_line_surf,period,fname_surf,&
         div3d_bfield_method)

    !----------------------------------------------------------
    ! 4. Initialize points on surface to carry heat
    !----------------------------------------------------------
    Write(6,*) 'Initializing points along initial surface line'
    If (.true.) Call init_points_line(fname_surf,npts_start,fname_launch)
  Else
    Write(6,*) 'Skipping trace_surface, loading file.'
  Endif

  !----------------------------------------------------------------------------------------
  ! 5. Follow fieldlines from init points and check for intersections 
  !----------------------------------------------------------------------------------------
  nhitline = Floor(hit_length/Rstart/abs(dphi_line_diff))
  Write(6,'(A,I0,A)') ' Returning ',nhitline,' points along intersecting lines'

  Write(6,'(/A)') '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
  write(6,*) 'Beginning MPI line diffusion'
  Call diffuse_lines3(fname_launch,dmag,dphi_line_diff,ns_line_diff, &
       fname_hit,period,fname_intpts,fname_nhit,nhitline, &
       lsfi_tol)

  Deallocate(ntri_parts)
  Deallocate(xtri,ytri,ztri,check_tri)
  Deallocate(xmid,ymid,zmid)
  Deallocate(R_ves,Z_ves,P_ves)
  Deallocate(dmid)

Endif ! rank == 0

!----------------------------------------------------------------
! Additional nodes wait for lines and follow them
!----------------------------------------------------------------
If (rank .gt. 0) Then

  num_myjobs = 0
  line_start_data_i = 0
!  Write(*,*) 'Process ',rank,' reporting as READY'
  Do While (line_start_data_i(1) .ne. -1) 

    ! Wait for line data
    source = 0
    dest = 0
    tag = rank
    Call MPI_RECV(line_start_data_i,3,MPI_INTEGER         ,source,tag,MPI_COMM_WORLD,status,ierr_mpi)   

    ! Check for kill signal    
    if (line_start_data_i(1) .ne. -1 ) Then
      ! QQ == diffuse_lines should really send a real array of zeros or something
      ! Get rest of start data
      Call MPI_RECV(line_start_data_r,6,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,ierr_mpi)

      ! Follow the line
      nhitline = line_start_data_i(2)
      Allocate(r_hitline(nhitline))
      Allocate(z_hitline(nhitline))
      Allocate(phi_hitline(nhitline))      
      ! QQ -- clean this up
!      Rstart_local = 
!Subroutine line_follow_and_int(Rstart,Zstart,Phistart,dphi_line,nsteps_line,dmag,period,pint,iout, &
!r_hitline,z_hitline,phi_hitline,nhitline,linenum,lsfi_tol)
      Call line_follow_and_int(line_start_data_r(1),line_start_data_r(2), &
           line_start_data_r(3),line_start_data_r(4),line_start_data_i(1),&
           line_start_data_r(5),line_start_data_r(6),pint,iout,r_hitline, &
           z_hitline,phi_hitline,nhitline,line_start_data_i(3),lsfi_tol,div3d_bfield_method)      
      ierr_follow = 0

      ! Compile results and send data back to master

      ! first send handshake signal (this_job_done)
      Call MPI_SEND(1,1,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,status,ierr_mpi)

      Allocate(line_done_data_r(3+3*nhitline))
      line_done_data_i(1) = ierr_follow
      line_done_data_i(2:5) = iout
      line_done_data_r(1:3) = pint
      line_done_data_r(3+1+0*nhitline:3+1*nhitline) = r_hitline
      line_done_data_r(3+1+1*nhitline:3+2*nhitline) = z_hitline
      line_done_data_r(3+1+2*nhitline:3+3*nhitline) = phi_hitline
      Call MPI_SEND(line_done_data_i,5,MPI_INTEGER                    ,dest,tag,MPI_COMM_WORLD,status,ierr_mpi)
      Call MPI_SEND(line_done_data_r,3+nhitline*3,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,status,ierr_mpi)
      Deallocate(r_hitline,z_hitline,phi_hitline)
      Deallocate(line_done_data_r)
      num_myjobs = num_myjobs + 1
    Endif ! kill signal check

  EndDo ! while mywork ne -1

  Write(*,*) 'Process ',rank,'received signal that all work is complete'
  Write(*,*) 'Process ',rank,' completed ',num_myjobs,' jobs'

Endif ! rank > 0


! Finialize MPI
Call fin_mpi()

end program div3d_follow_and_int






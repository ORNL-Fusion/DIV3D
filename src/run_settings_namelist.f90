! --Read namelist
! --Init random number generator
! --basic conversions
Module run_settings_namelist
  Use kind_mod, Only : real64, int32
  Implicit None

  Real(real64) :: Rstart, Zstart, Phistart
  Real(real64) :: dmag, lsfi_tol
  Real(real64) :: dphi_line_surf, dphi_line_diff, period
  Real(real64) :: dphi_line_surf_deg, dphi_line_diff_deg, hit_length

  Integer(int32) :: myseed, nhitline
  Integer(int32) :: npts_start

  Integer(int32) :: ntran_surf, ns_line_surf
  Integer(int32) :: ntran_diff, ns_line_diff
  
  Logical :: trace_surface_opt, calc_lc, calc_theta, quiet_bfield
  Logical :: vessel_int_is_last_point
  Logical :: vessel_is_nearest_slice
  
  Character(len=300) :: fname_hit, fname_ptri, fname_ptri_mid
  Character(len=300) :: fname_launch,fname_surf, fname_parts, fname_intpts, fname_ves
  Character(len=300) :: fname_plist, fname_nhit  

  
  Namelist / run_settings / fname_plist, fname_ves,  &
       fname_surf, fname_launch, fname_parts, fname_hit, fname_intpts, &
       Rstart, Zstart, Phistart, dphi_line_surf_deg, ntran_surf, &
       npts_start, dmag, dphi_line_diff_deg, ntran_diff, myseed, &
       fname_nhit, hit_length, lsfi_tol, trace_surface_opt, &
       fname_ptri, fname_ptri_mid, calc_lc, calc_theta, quiet_bfield, &
       vessel_is_nearest_slice, vessel_int_is_last_point
  
Contains
  

  Subroutine read_run_settings_namelist(verbose)    
    Use io_unit_spec, Only: iu_nl 
    Use phys_const, Only : pi
    Use parallel_mod, Only : fin_mpi, rank, nprocs
    Use setup_bfield_module, Only : bfield_nml, nfp_bfield
    Use init_random, Only : init_random_seed
    Use bfield, Only : verbose_bfield
    Implicit None
    Logical, Intent(In) :: verbose
    Character(len=256) :: iomsg
    Integer(int32) :: iocheck

    ! --------------------------------
    ! Defaults
    ! --------------------------------
    calc_lc = .true.
    calc_theta = .false.
    quiet_bfield = .true.
    vessel_is_nearest_slice = .true.
    vessel_int_is_last_point = .true.

    fname_plist    = 'parts.list'
    fname_ves      = 'vessel.part'
    fname_surf     = 'surface_line.out'
    fname_launch   = 'launch_pts.out'
    fname_parts    = 'allparts.out'
    fname_hit      = 'hitline.out'
    fname_intpts   = 'int_pts.out'
    fname_nhit     = 'hitcount.out'
    fname_ptri     = 'part_triangles.out'
    fname_ptri_mid = 'part_triangle_mids.out'

    lsfi_tol = 1.d-12
    hit_length = -1._real64
    
    ! -------------------------------
    ! Read settings namelists
    ! -------------------------------
    
    ! Open file
    If (verbose) Write(*,*) 'Reading run settings from run_settings.nml'
    Open(iu_nl,file="run_settings.nml",status="old",form="formatted",iostat=iocheck)
    If ( iocheck .ne. 0 ) Then
       if (rank .eq. 0) Write(*,*) 'Error opening namelist file'
       Call fin_mpi(.true.)
    Endif

    ! Read run_settings namelist
    Read(iu_nl,nml=run_settings,iostat=iocheck,iomsg=iomsg)
    If ( iocheck .ne. 0 ) Then
       if (rank .eq. 0) Write(*,*) 'Error reading run_settings namelist: ',trim(iomsg)
       Call fin_mpi(.true.)
    Endif

    ! Read bfield namelist
    Rewind(iu_nl)
    Read(iu_nl,nml=bfield_nml,iostat=iocheck,iomsg=iomsg)
    If ( iocheck .ne. 0 ) Then
       if (rank .eq. 0) Write(*,*) 'Error reading bfield_nml namelist: ',trim(iomsg)
       Call fin_mpi(.true.)
    Endif
    Close(iu_nl)

    ! -------------------------    
    ! Check inputs
    ! -------------------------

    ! vessel_is_nearest_slice
    If (vessel_is_nearest_slice) Then
       If (verbose) Write(*,*) 'Vessel representation is: Nearest slice'
    Else
       If (verbose) Write(*,*) 'Vessel representation is: Interpolated at constant phi index'
       ! ToDo: need to check at least 2 slices or some minimum phi resolution
    End If
    
    ! npts_start
    If (npts_start .lt. nprocs) Then
       if (rank .eq. 0) Write(*,*) 'Error: nprocs must be less than npts_start!!!'
       Call fin_mpi(.true.)
    Endif   

    ! nfp_bfield, default 0 comes from setup_bfield
    If (nfp_bfield .eq. 0) Then
       if (rank .eq. 0) Write(*,*) 'Error: nfp_bfield must be set in bfield_nml!'
       Call fin_mpi(.true.)
    Endif


    
    ! --------------------------
    ! Process inputs
    ! --------------------------
    
    period = 2.d0*pi/Real(nfp_bfield,real64)
    If (verbose) Write(*,*) 'Magnetic field set to have tor. symm. mode number ',nfp_bfield,period*180./pi,' deg.'
    dphi_line_surf = dphi_line_surf_deg * pi/180.d0
    dphi_line_diff = dphi_line_diff_deg * pi/180.d0
    ns_line_surf = Floor(ntran_surf*2.d0*pi/Abs(dphi_line_surf))
    ns_line_diff = Floor(ntran_diff*2.d0*pi/Abs(dphi_line_diff))

    If (verbose) Write(*,*) 'Initializing random number with base seed:',myseed
    Call init_random_seed(myseed*(rank+1))
    
    If (verbose) Write(*,'(A,G0.3)') ' hit_length: ',hit_length
    If (hit_length .le. 0.d0) Then
       If (verbose) Write(*,'(A)') ' Turning off hitline because hit_length <= 0'
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
    
    verbose_bfield = .not. quiet_bfield
    
  End Subroutine read_run_settings_namelist
  
End Module run_settings_namelist

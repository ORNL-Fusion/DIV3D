! --Read namelist
! --Init random number generator
! --basic conversions
Module run_settings_namelist
  Use kind_mod, Only : real64, int32
  Implicit None

  Real(real64) :: Rstart, Zstart, Phistart
  Real(real64) :: dmag, period, lsfi_tol
  Real(real64) :: dphi_line_surf, dphi_line_diff
  Real(real64) :: dphi_line_surf_deg, dphi_line_diff_deg, hit_length

  Integer(int32) :: myseed, nhitline
  Integer(int32) :: npts_start, nfp

  Integer(int32) :: ntran_surf, ns_line_surf
  Integer(int32) :: ntran_diff, ns_line_diff
  Integer(int32), Private :: iocheck
  
  Logical :: trace_surface_opt, calc_lc, calc_theta, quiet_bfield
  
  Character(len=300) :: fname_hit, fname_ptri, fname_ptri_mid
  Character(len=300) :: fname_launch,fname_surf, fname_parts, fname_intpts, fname_ves
  Character(len=300) :: fname_plist, fname_nhit  

  
  Namelist / run_settings / fname_plist, fname_ves,  &
       fname_surf, fname_launch, fname_parts, fname_hit, fname_intpts, nfp, &
       Rstart, Zstart, Phistart, dphi_line_surf_deg, ntran_surf, &
       npts_start, dmag, dphi_line_diff_deg, ntran_diff, myseed, &
       fname_nhit, hit_length, lsfi_tol, trace_surface_opt, &
       fname_ptri, fname_ptri_mid, calc_lc, calc_theta, quiet_bfield

  
Contains
  

  Subroutine read_run_settings_namelist(verbose)    
    Use io_unit_spec, Only: iu_nl 
    Use phys_const, Only : pi
    Use parallel_mod, Only : fin_mpi, rank
    Use setup_bfield_module, Only : bfield_nml
    Use init_random, Only : init_random_seed
    Use bfield, Only : verbose_bfield
    Implicit None
    Logical, Intent(In) :: verbose
    
    ! Defaults
    calc_lc = .true.
    calc_theta = .false.
    quiet_bfield = .true.

    verbose_bfield = .not. quiet_bfield
    
    If (verbose) Write(*,*) 'Reading run settings from run_settings.nml'

    Open(iu_nl,file="run_settings.nml",status="old",form="formatted",iostat=iocheck)
    If ( iocheck .ne. 0 ) Then
       Write(*,*) 'Error opening namelist file'
       Call fin_mpi(.true.)
    Endif
    Read(iu_nl,nml=run_settings,iostat=iocheck)
    If ( iocheck .ne. 0 ) Then
       Write(*,*) 'Error reading run_settings namelist'
       Call fin_mpi(.true.)
    Endif    
    Rewind(iu_nl)
    Read(iu_nl,nml=bfield_nml,iostat=iocheck)
    If ( iocheck .ne. 0 ) Then
       Write(*,*) 'Error reading bfield_nml namelist'
       Call fin_mpi(.true.)
    Endif
    Close(iu_nl)

    
    period = 2.d0*pi/Real(nfp,real64)
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
    
    
  End Subroutine read_run_settings_namelist
  
End Module run_settings_namelist

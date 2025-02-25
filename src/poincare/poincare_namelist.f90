Module poincare_namelist
  Use kind_mod, Only : real64, int32
  Implicit None

  Real(real64) :: period, dphi_line_poincare, dphi_slice
  Integer(int32) :: ind_poin, nsteps
  
  ! Namelist variables:
  Real(real64) :: phistart_deg_poincare, rstart_poincare, rend_poincare, zstart_poincare, zend_poincare, dphi_line_deg_poincare
  Integer(int32) :: num_surfs, ntransits, num_slices, num_points_per_surf
  Logical :: follow_both_ways, quiet_bfield
  
  ! Namelist files
  Namelist / poincare_settings_nml / phistart_deg_poincare, Rstart_poincare, Rend_poincare, &
       Zstart_poincare, Zend_poincare, num_surfs, ntransits, dphi_line_deg_poincare, follow_both_ways, &
       num_slices, quiet_bfield

Contains

  Subroutine read_poincare_namelist(verbose)
    Use parallel_mod, Only : fin_mpi, rank
    Use setup_bfield_module, Only : bfield_nml, nfp_bfield
    Use bfield_module, Only : verbose_bfield
    Use phys_const, Only : pi    
    Implicit None
    Logical, Intent(In) :: verbose
    Character(len=256) :: iomsg
    Integer(int32) :: iocheck, i
    Real(real64) :: Adphirat

    ! --------------------------------
    ! Defaults
    ! --------------------------------
    num_surfs = 2
    ntransits = 1
    follow_both_ways = .false.
    quiet_bfield = .true.
    num_slices = 1

    ! -------------------------------
    ! Read settings namelists
    ! -------------------------------    
    
    ! Open file
    If (verbose) Write(*,*) 'Reading poincare settings from poincare_settings.nml'
    Open(99,file="run_settings.nml",status="old",form="formatted",iostat=iocheck)
    If ( iocheck .ne. 0 ) Then
       if (rank .eq. 0) Write(*,*) 'Error opening namelist file'
       Call fin_mpi(.true.)
    Endif

    ! Read run_settings namelist
    Read(99,nml=poincare_settings_nml,iostat=iocheck,iomsg=iomsg)
    If ( iocheck .ne. 0 ) Then
       if (rank .eq. 0) Write(*,*) 'Error reading run_settings namelist: ',trim(iomsg)
       Call fin_mpi(.true.)
    Endif

    ! Read bfield namelist
    Rewind(99)
    Read(99,nml=bfield_nml,iostat=iocheck,iomsg=iomsg)
    If ( iocheck .ne. 0 ) Then
       if (rank .eq. 0) Write(*,*) 'Error reading bfield_nml namelist: ',trim(iomsg)
       Call fin_mpi(.true.)
    Endif
    Close(99)

    ! --------------------------
    ! Process inputs
    ! --------------------------
    
    period = 2.d0*pi/Real(nfp_bfield,real64)
    If (verbose) Write(*,*) 'Magnetic field set to have toroidal symmetry mode number ',nfp_bfield,period*180./pi,' deg.'    

    verbose_bfield = .not. quiet_bfield

    ! Step size in radians
    dphi_line_poincare = dphi_line_deg_poincare*pi/180.d0
    
    ! Number of steps to complete the required number of transits. ntransits is full toroidal transits.    
    nsteps = Floor(ntransits*2.d0*pi/Abs(dphi_line_poincare))
    
    ! Check that step size fits into the symmetric toroidal angle, --> also account for slices
    Adphirat = Abs(360.d0/dphi_line_deg_poincare/Real(nfp_bfield,real64)/Real(num_slices,real64))
    If (Adphirat - Real(Nint(Adphirat)) > 1.d-8) Then
       Write(*,*) 'Error!: 2*pi/nfp_bfield/num_slices must be an integer multiple of dphi_line_deg'
       Write(*,*) 'Exiting'
       Call fin_mpi(.true.)
    Endif

    ! This is the stride required to stay at constant toroidal angle
    ind_poin = Nint(Adphirat*num_slices)

    ! Step size between slices
    dphi_slice = 2.d0*pi/Real(nfp_bfield,real64)/Real(num_slices,real64)

    ! Number of points returns per surface per slice    
    num_points_per_surf = ntransits*nfp_bfield

    If (verbose) Then
       Write(*,'(/1a,i0,a)')      'Fieldlines are traced for ',ntransits,' transits.'    
       Write(*,'(1a,i0,a,i0,a)')  'Tracing ',num_surfs,' surfaces with ',nsteps,' steps'
       Write(*,'(a,g0.3,a)')      'Toroidal step size for tracing is ',dphi_line_deg_poincare,' degrees.'
       Write(*,'(a,i0,a,i0,a)')   'With nfp = ',nfp_bfield,' there will be ~',num_points_per_surf,' points per surface'
       Write(*,'(a,i0)')          'Stride for poincare plot is ',ind_poin
       Write(*,'(a,i0,a,i0,a)')   'Saving data for ',num_surfs,' surfaces with ',num_slices,' slices (toroidal angles)'
       If (num_slices .lt. 20) Then
          Write(*,'(a,g0.3,a)')      'Toroidal angles for slices with dphi_slice = ',dphi_slice*180.d0/pi,' degrees'       
          Do i = 0,num_slices-1
             Write(*,'(f12.3)') dphi_slice*i*180.d0/pi
          End Do
       End If
    End If
        
  End Subroutine read_poincare_namelist

End Module poincare_namelist

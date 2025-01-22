Module poincare_namelist
  Use kind_mod, Only : real64, int32
  Implicit None

  Real(real64) :: period
  
  ! Namelist variables:
  Real(real64) :: phistart_deg_poincare, rstart_poincare, rend_poincare, zstart_poincare, zend_poincare, dphi_line_deg_poincare
  Integer(int32) :: num_pts, ntransits
  Logical :: follow_both_ways, quiet_bfield
  
  ! Namelist files
  Namelist / poincare_settings_nml / phistart_deg_poincare, Rstart_poincare, Rend_poincare, &
       Zstart_poincare, Zend_poincare, num_pts, ntransits, dphi_line_deg_poincare, follow_both_ways, quiet_bfield

Contains

  Subroutine read_poincare_namelist(verbose)
    Use parallel_mod, Only : fin_mpi, rank, nprocs
    Use setup_bfield_module, Only : bfield_nml, nfp_bfield
    Use bfield_module, Only : verbose_bfield
    Use phys_const, Only : pi    
    Implicit None
    Logical, Intent(In) :: verbose
    Character(len=256) :: iomsg
    Integer(int32) :: iocheck    

    ! --------------------------------
    ! Defaults
    ! --------------------------------
    num_pts = 2
    ntransits = 1
    follow_both_ways = .false.
    quiet_bfield = .true.

    ! -------------------------------
    ! Read settings namelists
    ! -------------------------------    
    
    ! Open file
    If (verbose) Write(*,*) 'Reading poincare settings from poincare_settings.nml'
    Open(99,file="poincare_settings.nml",status="old",form="formatted",iostat=iocheck)
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


  End Subroutine read_poincare_namelist

End Module poincare_namelist

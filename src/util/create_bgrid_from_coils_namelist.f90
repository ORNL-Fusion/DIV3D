Module create_bgrid_from_coils_namelist
  Use kind_mod, Only : real64, int32
  Implicit None

  Real(real64) :: period

  ! Namelist variables
  Real(real64) :: Rmin, Rmax, Zmin, Zmax
  Integer(int32) :: nr, nz, nphi
  Logical :: quiet_bfield
  Character(len=256) :: outfile_prefix

  ! Namelist file
  Namelist / create_bgrid_settings_nml / Rmin, Rmax, Zmin, Zmax, &
       nr, nz, nphi, quiet_bfield, outfile_prefix

Contains

  Subroutine read_create_bgrid_from_coils_namelist(verbose)
    Use parallel_mod, Only : fin_mpi, rank
    Use setup_bfield_module, Only : bfield_nml, nfp_bfield
    Use bfield_module, Only : verbose_bfield
    Use phys_const, Only : pi
    Implicit None
    Logical, Intent(In) :: verbose
    Character(len=256) :: iomsg
    Integer(int32) :: iocheck, i

    ! --------------------------------
    ! Defaults
    ! --------------------------------
    Rmin = 1.d0
    Rmax = 3.d0
    Zmin = -1.d0
    Zmax = 1.d0
    
    nr = 3
    nz = 3
    ! Number of planes including 0 and 2*pi/nfp_bfield.
    ! The last (redundant) slice is not written to the grid file
    nphi = 3

    outfile_prefix = 'none'

    ! -------------------------------
    ! Read settings namelists
    ! -------------------------------    
    
    ! Open file
    If (verbose) Write(*,*) 'Reading create_bgrid settings from run_settings.nml'
    Open(99,file="run_settings.nml",status="old",form="formatted",iostat=iocheck)
    If ( iocheck .ne. 0 ) Then
       if (rank .eq. 0) Write(*,*) 'Error opening namelist file'
       Call fin_mpi(.true.)
    Endif

    ! Read run_settings namelist
    Read(99,nml=create_bgrid_settings_nml,iostat=iocheck,iomsg=iomsg)
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
    If (verbose) Write(*,*) 'Magnetic field set to have toroidal symmetry mode number ',nfp_bfield,period*180./pi,'deg.'

  End Subroutine read_create_bgrid_from_coils_namelist

End Module create_bgrid_from_coils_namelist

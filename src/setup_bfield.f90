Module setup_bfield_module
  Use kind_mod, Only: int32, real64
  Use bfield_module, Only : bfield_type, coil_type, g_type
  
  Implicit None
  
  Type(bfield_type) :: bfield
  Type(g_type) :: g  
  Type(coil_type) :: coil
  
  Integer(int32), Parameter :: max_extcur = 100
  Logical :: setup_bfield_verbose = .true. ! To be used by calling routines to supress output (particularly in MPI runs)
  ! ------------ BFIELD NAMELIST VARIABLES ----------------
  
  Real(real64) :: &
       vmec_extcur_set(max_extcur)            = 0.d0

  Integer(int32) :: &
       nfp_bfield = 0
  
  Character(Len=300) :: &
       rmp_type                         = 'none', &
       gfile_name                       = 'none', &       
       vmec_coils_file                  = 'none', &
       xdr_fname                        = 'none', &
       bgrid_fname                      = 'none'

  Logical :: &
       xdr_check   = .true.,   &
       xdr_verbose = .true.
  
  Namelist / bfield_nml / &
       rmp_type, &
       gfile_name, &       
       vmec_coils_file, vmec_extcur_set, &
       xdr_fname, xdr_check, xdr_verbose, &
       bgrid_fname, nfp_bfield
  
  ! ------------ BFIELD NAMELIST VARIABLES ----------------
  
!  Private

Contains
  ! *********************************************
  ! *************** GFILE ONLY ******************
  ! *********************************************
  Subroutine setup_bfield_g3d
    Use g3d_module, Only : readg_g3d
    Implicit None
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS G3D'
    bfield%method      = 0
    bfield%method_2d   = 0
    bfield%method_pert = -1
    Call readg_g3d(gfile_name,g)
    bfield%g = g
  End Subroutine setup_bfield_g3d

  ! *********************************************
  ! *************** VMEC COILS ******************
  ! *********************************************
  Subroutine setup_bfield_vmec_coils
    ! Requires vmec_coils_file and vmec_extcur_set are set
    Use vmec_routines_mod, Only : read_vmec_coils_file, vmec_extcur, vmec_nextcur    
    Implicit None
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS VMEC COILS'
    bfield%method      = 14
    bfield%method_2d   = -1
    bfield%method_pert = -1    
    Call read_vmec_coils_file(vmec_coils_file)
    vmec_extcur(1:vmec_nextcur) = vmec_extcur_set(1:vmec_nextcur)

  End Subroutine setup_bfield_vmec_coils

  ! *********************************************
  ! *************** VMEC COILS TO FIL************
  ! *********************************************
  Subroutine setup_bfield_vmec_coils_to_fil
    ! Requires vmec_coils_file and vmec_extcur_set are set
    Use vmec_routines_mod, Only : read_vmec_coils_file, vmec_extcur, &
         vmec_nextcur, vmec_coil_new, convert_vmec_coils_to_filaments
    Implicit None
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS VMEC COILS TO FIL'
    bfield%method      =  6
    bfield%method_2d   = -1
    bfield%method_pert = -1    
    Call read_vmec_coils_file(vmec_coils_file)
    vmec_extcur(1:vmec_nextcur) = vmec_extcur_set(1:vmec_nextcur)
    Call convert_vmec_coils_to_filaments
    bfield%coil = vmec_coil_new
  End Subroutine setup_bfield_vmec_coils_to_fil
    
#ifdef HAVE_FXDR    
  ! *********************************************
  ! ***************  XDR  ***********************
  ! *********************************************    
  Subroutine setup_bfield_xdr
    Use xdr_routines_mod, Only : readbgrid_xdr
    Implicit None
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS XDR'
    Call readbgrid_xdr(xdr_fname,xdr_check,xdr_verbose)
    bfield%method      = 15
    bfield%method_2d   = -1
    bfield%method_pert = -1
  End Subroutine setup_bfield_xdr
#endif
  
  ! *********************************************
  ! *************** BGRID ***********************
  ! *********************************************    
  Subroutine setup_bfield_bgrid
    Use bgrid_module, Only: open_bgrid_fields
    Implicit None
    
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS BGRID'
      bfield%method      = 16
      bfield%method_2d   = -1
      bfield%method_pert = -1
    Call open_bgrid_fields(bgrid_fname,setup_bfield_verbose)
  End Subroutine setup_bfield_bgrid

  
End Module setup_bfield_module
  

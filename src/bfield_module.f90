Module bfield_typedef
  Use kind_mod, Only : int32, real64
  Use coil_typedef, Only : coil_type
  Implicit None
  Private
  Type, Public :: bfield_type
    Integer(int32) :: method = -1
    Type(coil_type) :: coil
    Integer(int32) :: method_2d = -1    ! Method corresponding to AS fields
    Integer(int32) :: method_pert = -1  ! Method corresponding to pert only
    Integer(int32) :: method_save = -1  ! Used to save standard method
    Logical :: method_switched = .false.
  End Type bfield_type
End Module bfield_typedef

!   Set bfield_method to control the fieldline deriviative calls.
!   Appropriate loading must be done before calls to any fieldline following routine 
!
!    bfield%method == 
!                     6 -- Just coils
!                    14 -- VMEC coils file with extcur
!                    15 -- xdr file
!                    16 -- bgrid file (FLARE format)
! 


Module bfield
  Use bfield_typedef, Only : bfield_type
  Use coil_typedef, Only : coil_type  
  Implicit None
  Private
  Public :: bfield_type, coil_type
  Public :: calc_B_rzphi_general

Contains

  Subroutine calc_B_rzphi_general(bfield,r,z,phi,n,br,bz,bphi,ierr_out)
    Use kind_mod, Only: real64, int32
    Use biotsavart_module, Only : bfield_bs_cyl
    Use VMEC_routines_mod, Only : bfield_vmec_coils
#ifdef HAVE_FXDR    
    Use xdr_routines_mod, Only : bint_xdr_n
#endif    
    Use bgrid_module, Only : bfield_bgrid
    Implicit None
    Type(bfield_type), Intent(In) :: bfield
    Integer(int32), Intent(In) :: n
    Real(real64), Intent(In) :: r(n),z(n),phi(n)
    Real(real64), Intent(Out) :: br(n),bz(n),bphi(n)
    Integer(int32), Intent(Out), Optional :: ierr_out
    Real(real64) :: btmp(n,3)
    Integer(int32) :: ierr

    ierr = 0
    If (Present(ierr_out)) ierr_out = 0

    br   = 0._real64
    bz   = 0._real64
    bphi = 0._real64
    btmp = 0._real64
    
    Select Case (bfield%method)
    Case (6) ! just coils
      Call bfield_bs_cyl(r,phi,z,n,bfield%coil,br,bphi,bz)
    Case (14) ! VMEC coils
      Call bfield_vmec_coils(r,phi,z,n,btmp,ierr)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)
    Case (15) ! Xdr
#ifdef HAVE_FXDR         
      Call bint_xdr_n(r,phi,z,n,btmp,ierr)
      br   = btmp(:,1)
      bz   = btmp(:,3)  ! Note order!
      bphi = btmp(:,2)
#else
      Stop "Compiled without fxdr support"
#endif      
    Case (16) ! bgrid
      Call bfield_bgrid(r,phi,z,n,btmp,ierr)
      br   = btmp(:,1)
      bz   = btmp(:,2)
      bphi = btmp(:,3)
    Case Default
      Write(*,*) 'Unknown bfield%method:',bfield%method
      Stop "Exiting from bfield general"
    End Select
    If (Present(ierr_out)) ierr_out = ierr
  End Subroutine calc_B_rzphi_general
End Module bfield


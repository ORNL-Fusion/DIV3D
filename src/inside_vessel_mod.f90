Module inside_vessel_mod
Use kind_mod, Only : real64, int32
Implicit None

Real(real64), Allocatable, Private :: Rslice(:), Zslice(:)
  
Contains

!-----------------------------------------------------------------------------  
Subroutine init_find_vessel_intersection
Use read_parts_mod, Only : npol_ves
Implicit None
!- End of header -------------------------------------------------------------
Allocate(Rslice(npol_ves),Zslice(npol_ves))
End Subroutine init_find_vessel_intersection


!-----------------------------------------------------------------------------  
Subroutine fin_find_vessel_intersection
Implicit None
!- End of header -------------------------------------------------------------
Deallocate(Rslice,Zslice)
End Subroutine fin_find_vessel_intersection


!-----------------------------------------------------------------------------
Function inside_vessel(R,Z,Phiin,Rves,Zves,Pves,ntor,npol) &
Result(inside)

Use kind_mod, Only : real64, int32
Use math_routines_mod, Only : wrap_phi, inside_poly
Use run_settings_namelist, Only : period
Implicit None

Real(real64), Intent(in) :: R, Z, Phiin
Logical :: inside
Integer(int32), Intent(in) :: ntor,npol
Real(real64), Dimension(ntor), Intent(in) :: Pves
Real(real64),Dimension(ntor,npol), Intent(in) :: Rves, Zves

Real(real64) :: Rs(npol), Zs(npol), Phi
!- End of header -------------------------------------------------------------

Phi = Phiin
Call wrap_phi(Phi,period)

! Get vessel at phi
Call get_vessel_at_phi(Phi, Pves, Rves, Zves, Rs, Zs, ntor, npol)

! Check polygon -- could drop one point but probably small impact
inside = inside_poly(R,Z,Rs,Zs,npol)

End Function inside_vessel

!-----------------------------------------------------------------------------
Subroutine get_vessel_at_phi(Phi, Pves, Rves, Zves, Rs, Zs, ntor, npol)
Use kind_mod, Only : real64, int32
Use read_parts_mod, Only : is_AS_ves
Use run_settings_namelist, Only : vessel_is_nearest_slice
Use parallel_mod, Only : fin_mpi
Implicit None

Real(real64), Intent(in) :: Phi
Integer(int32), Intent(in) :: ntor, npol
Real(real64), Dimension(ntor), Intent(in) :: Pves
Real(real64), Dimension(ntor, npol), Intent(in) :: Rves, Zves
Real(real64), Dimension(npol), Intent(out) :: Rs, Zs
Integer(int32) :: ind_near
!- End of header -------------------------------------------------------------

If (vessel_is_nearest_slice) Then
   If (is_AS_ves) Then
      Rs = Rves(1,:)
      Zs = Zves(1,:)
   Else
      ind_near = Minloc(Abs(Phi-Pves),1)
      Rs = Rves(ind_near,:)
      Zs = Zves(ind_near,:)
   End If
Else
   Write(*,*) "Error: vessel interpolation not yet implemented"
   Call fin_mpi(.true.)
End If

End Subroutine get_vessel_at_phi

!-----------------------------------------------------------------------------
!
! Must call init_find_vessel_intersection first!
! This all assumes that the vessel has a constant npol
! JDL 12/2024
Subroutine find_vessel_intersection(R1,Z1,R2,Z2,Phi,Rint,Zint)
Use kind_mod, Only : real64, int32
Use read_parts_mod, Only : R_ves, Z_ves, P_ves, npol_ves, ntor_ves
Use math_routines_mod, Only : int_line_curve
Implicit None
Real(real64), Intent(in) :: R1,Z1,R2,Z2,Phi
Real(real64), Intent(out) :: Rint, Zint
Real(real64) :: pint2D(2), Uint
Integer(int32) :: ierr_pint
!- End of header -------------------------------------------------------------

! This is not perfect because it just uses the vessel at Phi1
Call get_vessel_at_phi(Phi, P_ves, R_ves, Z_ves, Rslice, Zslice, ntor_ves, npol_ves)
Call int_line_curve((/R1,Z1/),(/R2,Z2/),Rslice,Zslice,.true.,pint2D,ierr_pint,Uint)

Rint = pint2D(1)
Zint = pint2D(2)

End Subroutine find_vessel_intersection


End Module inside_vessel_mod

Module inside_vessel_mod

Contains

Function inside_vessel(Rin,Zin,Phiin,Rves,Zves,Pves,ntor,npol,msym) &
Result(inside)

Use kind_mod
Use phys_const, Only : pi
Implicit None

Real(real64), Intent(in) :: Rin,Zin,Phiin
Integer(int32) :: inside
Integer(int32), Intent(in) :: ntor,npol
Integer(int32),Intent(in) :: msym
Real(real64), Dimension(ntor), Intent(in) :: &
  Pves
Real(real64),Dimension(ntor,npol), Intent(in) :: &
  Rves, Zves

Real(real64) :: R,Z,Phi
Real(real64) :: Rs(npol), Zs(npol)
Integer(int32) :: ind_near
!- End of header -------------------------------------------------------------



R=Rin
Z=Zin
Phi=Phiin

Do While (Phi .lt. 0.d0) 
  Phi = Phi + 2._real64*pi/Real(msym,real64)
Enddo
Phi = Mod(Phi,2._real64*pi/Real(msym,real64))

ind_near = Minloc(Dabs(Phi-Pves),1)
!write(*,*) ind_near,phi*180.d0/pi,pves(ind_near)*180.d0/pi


Rs = Rves(ind_near,:)
Zs = Zves(ind_near,:)

inside = inside_poly(R,Z,Rs(1:npol-1),Zs(1:npol-1),npol-1)

Endfunction inside_vessel


!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

Function inside_poly(x,y,px,py,npoly) &
Result(inside)
!
! Inputs:
!   npoly - size of vectors px and py (integer)
!                                ;  x - The x coordinate of the point.
!                                ;  y - The y coordinate of the point.
!                                ; px - The x coordinates of the polygon.
!                                ; py - The y coordinates of the polygon.
!                                ;
!                                ; The return value of the function is 1 if the point is inside the
!                                ; polygon and 0 if it is outside the polygon.

!
! Polygon is closed by connecting first and last points
!

Use kind_mod
Use phys_const, Only : pi
Implicit none

Real(real64), Intent(in) :: x,y
Integer(int32), Intent(in) :: npoly
Real(real64), Dimension(npoly), Intent(in) :: px,py

Integer(int32) :: inside

! Local variables
Real(real64), Dimension(npoly+1) :: tmp_px, tmp_py
Real(real64), Dimension(npoly) :: theta,dp,cp

<<<<<<< Updated upstream
Real(real64),parameter :: pi = 3.14159265358979323846_real64  

=======
>>>>>>> Stashed changes
!- End of header -------------------------------------------------------------

! Close polygon
tmp_px(1:npoly) = px
tmp_px(npoly+1) = px(1)
tmp_py(1:npoly) = py
tmp_py(npoly+1) = py(1)

dp = (tmp_px(1:npoly)-x)*(tmp_px(2:npoly+1)-x) + (tmp_py(1:npoly)-y)*(tmp_py(2:npoly+1)-y) ! dot product
cp = (tmp_px(1:npoly)-x)*(tmp_py(2:npoly+1)-y) - (tmp_py(1:npoly)-y)*(tmp_px(2:npoly+1)-x) ! cross product

theta = atan2(cp,dp)

inside = 0
If ( Abs(Sum(theta)) .gt. pi ) Then
  inside = 1
Endif

Endfunction inside_poly

Endmodule inside_vessel_mod

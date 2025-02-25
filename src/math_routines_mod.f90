!-----------------------------------------------------------------------------
!+ Contains special math functions
!  dist_2pts_cyl
!  wrap_phi
!  int_two_lines
!  rlinspace
!  line_seg_facet_int
!  inside_poly
!  int_line_curve
!-----------------------------------------------------------------------------
Module math_routines_mod
! Author(s): J. Lore 07/2009 - xxx 
Implicit None

Contains

!------------------------------------------------------------------------------
!+ Linear interpolation of curve C to find intersection with line segment L
!------------------------------------------------------------------------------
Subroutine int_line_curve(p1,p2,rline,zline,first,pint1,ierr,u1,found_ind,int_count_out)
! Curve is defined by array of points, L by p1,p2
!------------------------------------------------------------------------------
!+ Linear interpolation of curve C to find intersection with line L
!
! Inputs:
!   p1, p2       - Line endpoints defining line segment L (2D, R,Z).
!   rline, zline - Arrays defining curve C in cylindrical coordinates.
!   first        - Logical flag to return the first intersection found.
!
! Outputs:
!   pint1        - Intersection point (2D) between L and C.
!   ierr         - Error code:
!                    0 - Success
!                    1 - No intersection found
!   u1           - Parameter defining intersection point along L.
!   found_ind    - (Optional) Index of the segment of C where intersection occurs.
!   int_count_out - (Optional) Number of intersections found.
!
! Notes:
!   rline and zline must have the same size.
!------------------------------------------------------------------------------
!
! JL 2/2011-2024
Use kind_mod, Only: real64, int32
Implicit None
Logical, Intent(in) :: first
Real(real64), Intent(in) :: p1(2), p2(2), rline(:), zline(:)
Real(real64), Intent(out) :: pint1(2), u1
Integer(int32), Intent(out):: ierr
Integer(int32), Intent(out), Optional :: found_ind, int_count_out
Logical :: test
Integer(int32) :: ii, ierr2, int_count, nline
Real(real64) :: p3(2), p4(2), u2

nline = Size(rline)

int_count = 0
If (Present(found_ind)) found_ind = 0
Do ii = 1,nline-1
  p3 = [rline(ii),zline(ii)]
  p4 = [rline(ii+1),zline(ii+1)]
  Call int_two_lines(p1,p2,p3,p4,u1,u2,ierr2)

  If (ierr2 .ne. 0) Cycle

  If (ii == nline -1) Then
    test = ((u1 .ge. 0._real64) .AND. (u1 .le. 1._real64) .AND. &
         (u2 .ge. 0._real64) .AND. (u2 .le. 1._real64))
  Else
    test = ((u1 .ge. 0._real64) .AND. (u1 .le. 1._real64) .AND. &
         (u2 .ge. 0._real64) .AND. (u2 .lt. 1._real64))
  Endif

  If (test) Then
    int_count = int_count + 1
    pint1 = p1 + u1*(p2-p1)
    ierr = 0
    If (Present(found_ind)) found_ind = ii
    If (first) Then
      If (Present(int_count_out)) int_count_out = int_count
      Return
    Endif
  Endif  
Enddo

If (int_count .eq. 0) Then
    pint1 = [0.d0,0.d0]
    ierr = 1
Endif
If (int_count > 1) Then
    Write(*,*) 'Warning: More than one intersection found in int_line_curve. Returning last point'
Endif
If (Present(int_count_out)) int_count_out = int_count

End Subroutine int_line_curve


!-----------------------------------------------------------------------------
!+ Distance formula in cylindrical coordinates
!-----------------------------------------------------------------------------
Function dist_2pts_cyl(R1, R2, Z1, Z2, P1, P2) Result(distance)
  Use kind_mod, Only : real64
  Implicit None

  Real(real64), Intent(In) :: R1, R2  ! Radial components
  Real(real64), Intent(In) :: Z1, Z2  ! Height components
  Real(real64), Intent(In) :: P1, P2  ! Azimuthal angle components (in radians)
  Real(real64) :: distance  ! Result: the distance between two points

  ! Distance formula in cylindrical coordinates
  distance = sqrt((R1**2 + R2**2 - 2.0_real64*R1*R2*cos(P1-P2)) + (Z1-Z2)**2)

End Function dist_2pts_cyl


!-----------------------------------------------------------------------------
!+ Wrap phi given a symmetry periodicity
!-----------------------------------------------------------------------------
Subroutine wrap_phi(Phi, period)
   Use kind_mod, Only : real64
   Implicit None

   Real(real64), Intent(InOut) :: Phi
   Real(real64), Intent(In) :: period

   ! Wrap Phi into the range [0, period)
   Phi = Mod(Phi, period)
   If (Phi .lt. 0._real64) Then
      Phi = Phi + period
   End If
End Subroutine wrap_phi

!-----------------------------------------------------------------------------
!+ Checks for the intersection of two line segments
!-----------------------------------------------------------------------------
Subroutine int_two_lines(p1,p2,p3,p4,u1,u2,ierr,pint)
! Calculates the intersection of lines p1 to p2, and p3 to p4.  Each
! point is an x-y pair.  Returns two values, which are the distance
! between p1 and p2 of the intersection, and between p3 and p4.  This
! distance is normalized to the length of the lines.

Use kind_mod, Only: real64, int32
implicit none

Real(real64), Intent(in),Dimension(2) :: p1,p2,p3,p4
Real(real64), intent(out) :: u1,u2
Real(real64) :: denom
Real(real64), Intent(out), Optional :: pint(2)
Integer(int32), Intent(out) :: ierr

denom = (p4(2)-p3(2))*(p2(1)-p1(1)) - (p4(1)-p3(1))*(p2(2)-p1(2))

if ( ABS(denom) < epsilon(1.0_real64) ) then
  u1 = 1.d30
  u2 = 1.d30
  If (Present(pint)) pint = 0._real64
  ierr = 1
else
    u1 = ((p4(1)-p3(1))*(p1(2)-p3(2)) - (p4(2)-p3(2))*(p1(1)-p3(1)))/denom
    u2 = ((p2(1)-p1(1))*(p1(2)-p3(2)) - (p2(2)-p1(2))*(p1(1)-p3(1)))/denom
    If (Present(pint)) pint = p1 + u1*(p2-p1)
    ierr = 0
endif

Endsubroutine int_two_lines
!-----------------------------------------------------------------------------

!------------------------------------------------------------------------------
!+ Returns a linearly spaced real vector given endpoints and number of elements
!------------------------------------------------------------------------------
Function rlinspace(xstart,xend,numel)  & 
     Result(rlinvec)
  !
  ! Description: 
  !   This function returns a real vector of length(numel) with linearly spaced
  !   values from xstart to xend.  Similar to the Matlab function.
  !
  !   If numel == 1, then xstart is returned (unlike Matlab which returns xend)
  !
  ! Inputs:
  !  xstart,xend: Values of the first and last points of the array [real]
  !  numel: Number of elements in the array [integer]
  ! Outputs:
  !  rlinvec: The linearly spaced array
  !
  ! Author(s): J. Lore 07/2009 - 7/18/2011
  !
  ! Modules used:
  Use kind_mod, Only : real64, int32
  
  Implicit None
  
  ! Input/output
  Real(real64),    Intent(in) :: xstart
  Real(real64),    Intent(in) :: xend
  Integer(int32), Intent(in) :: numel
  Real(real64)                :: rlinvec(numel)
  
  ! Local scalars
  Integer(int32)   ::  ii
  !- End of header -------------------------------------------------------------
  
  If (numel == 1) Then
     rlinvec(1) = xstart
     Return
  Else   
     Do ii = 1,numel
        rlinvec(ii) = ( Real(ii,real64) - 1._real64 ) * ( xend - xstart ) &
             / ( numel - 1._real64 ) + xstart
     Enddo
  End If

End Function rlinspace


!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------
Subroutine line_seg_facet_int(pa,pb,pc,p1,p2,ithit,p,tol,calc_theta,sin_theta)
!
! Description: 
!  Based on code from Paul Bourke
! Inputs: 
!
! Outputs:
!   none
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     12/07/2011   JDL
! Author(s): J.D. Lore - 12/07/2011 - xxx

! Modules used:
Use kind_mod, Only : real64, int32
Use phys_const, Only : pi
Implicit none

! Input/output
Real(real64), Dimension(3), Intent(In) :: pa,pb,pc,p1,p2
Real(real64), Dimension(3), Intent(Out) :: p
Integer(int32), Intent(Out) :: ithit
real(real64), Intent(In) :: tol
Real(real64), Intent(Out) :: sin_theta
Logical, Intent(In) :: calc_theta


!local variables and arrays
Real(real64), Dimension(3) :: n, pa1, pa2, pa3, v
Real(real64) :: d, denom, mu, total, a1, a2, a3, n2

Real(real64) :: mag_v, dot_nv


!- End of header -------------------------------------------------------------


!tol = 1.d-12

! Calc unit vector normal to plane of pa-pc (gives plane components A-C)
! Cross product of AB and AC
n(1) = (pb(2) - pa(2))*(pc(3) - pa(3)) - (pb(3) - pa(3))*(pc(2) - pa(2))
n(2) = (pb(3) - pa(3))*(pc(1) - pa(1)) - (pb(1) - pa(1))*(pc(3) - pa(3))
n(3) = (pb(1) - pa(1))*(pc(2) - pa(2)) - (pb(2) - pa(2))*(pc(1) - pa(1))
n2=Sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3))
n2 = 1._real64/n2
n=n*n2
!n=n/Sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3))

! Calculate plane component D
d = - n(1)*pa(1) - n(2)*pa(2) - n(3)*pa(3)

! Calculate the position on the line that intersects the plane
denom = n(1)*(p2(1) - p1(1)) + n(2)*(p2(2) - p1(2)) + n(3)*(p2(3) - p1(3))

p(:) = 0._real64

If (abs(denom) .lt. tol) Then
   ! line is parallel to plane
    ithit = 0
else
    ! how far along line intersection occurs [0,1]
    mu = - (d + n(1) * p1(1) + n(2) * p1(2) + n(3) * p1(3)) / denom
    p(1) = p1(1) + mu * (p2(1) - p1(1))
    p(2) = p1(2) + mu * (p2(2) - p1(2))
    p(3) = p1(3) + mu * (p2(3) - p1(3))
        
    if ((mu .lt. 0._real64) .or. (mu .gt. 1._real64)) Then   !Intersection not along line segment
        ithit = 0
    else        
        !  Determine whether or not the intersection point is bounded by pa,pb,pc
        pa1(1) = pa(1) - p(1)
        pa1(2) = pa(2) - p(2)
        pa1(3) = pa(3) - p(3)
        pa1 = pa1/Sqrt(pa1(1)*pa1(1) + pa1(2)*pa1(2) + pa1(3)*pa1(3))
        pa2(1) = pb(1) - p(1)
        pa2(2) = pb(2) - p(2)
        pa2(3) = pb(3) - p(3)
        pa2 = pa2/Sqrt(pa2(1)*pa2(1) + pa2(2)*pa2(2) + pa2(3)*pa2(3))
        pa3(1) = pc(1) - p(1)
        pa3(2) = pc(2) - p(2)
        pa3(3) = pc(3) - p(3)
        pa3 = pa3/Sqrt(pa3(1)*pa3(1) + pa3(2)*pa3(2) + pa3(3)*pa3(3))
        a1 = pa1(1)*pa2(1) + pa1(2)*pa2(2) + pa1(3)*pa2(3)
        a2 = pa2(1)*pa3(1) + pa2(2)*pa3(2) + pa2(3)*pa3(3)
        a3 = pa3(1)*pa1(1) + pa3(2)*pa1(2) + pa3(3)*pa1(3)
        ! all the min/max is to avoid acos(x) where abs(x) > 1 due to truncation error
        total = acos(max(-1._real64,min(1._real64,a1))) + acos(max(-1._real64,min(1._real64,a2))) &
             + acos(max(-1._real64,min(1._real64,a3)))

        if (abs(total - 2._real64*pi) .gt. tol) Then
            ithit = 0_int32
        else
            ithit = 1_int32
        endif
    endif
endif

sin_theta = 0._real64
If ((ithit .eq. 1_int32) .and. (calc_theta)) Then
!   Write(*,*) 'computing angle between the line and the plane'

   ! Below assumes |n| = 1, as done above
   
   ! Compute direction vector of the line
   v(1) = p2(1) - p1(1)
   v(2) = p2(2) - p1(2)
   v(3) = p2(3) - p1(3)

   ! Compute dot product between the normalized normal vector n and the direction vector v
   dot_nv = n(1) * v(1) + n(2) * v(2) + n(3) * v(3)

   ! Compute magnitudes of n and v
   mag_v = Sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))

   ! Calculate sin of the angle between the line and the plane
   sin_theta = dot_nv / (mag_v)

   ! Ensure sin_theta is within valid range due to potential floating-point inaccuracies
   sin_theta = max(-1.0_real64, min(1.0_real64, sin_theta))

Endif


End Subroutine line_seg_facet_int

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

Function inside_poly(x,y,px,py,npoly) &
Result(inside)
!
! Inputs:
!   npoly - size of vectors px and py (integer)
!   x - The x coordinate of the point.
!   y - The y coordinate of the point.
!   px - The x coordinates of the polygon.
!   py - The y coordinates of the polygon.
!

!
! Polygon is closed by connecting first and last points
!

Use kind_mod, Only : real64, int32
Use phys_const, Only : pi
Implicit none

Real(real64), Intent(in) :: x,y
Integer(int32), Intent(in) :: npoly
Real(real64), Dimension(npoly), Intent(in) :: px,py

Logical :: inside

! Local variables
Real(real64), Dimension(npoly+1) :: tmp_px, tmp_py
Real(real64), Dimension(npoly) :: theta,dp,cp
!- End of header -------------------------------------------------------------

! Close polygon
tmp_px(1:npoly) = px
tmp_px(npoly+1) = px(1)
tmp_py(1:npoly) = py
tmp_py(npoly+1) = py(1)

dp = (tmp_px(1:npoly)-x)*(tmp_px(2:npoly+1)-x) + (tmp_py(1:npoly)-y)*(tmp_py(2:npoly+1)-y) ! dot product
cp = (tmp_px(1:npoly)-x)*(tmp_py(2:npoly+1)-y) - (tmp_py(1:npoly)-y)*(tmp_px(2:npoly+1)-x) ! cross product

theta = atan2(cp,dp)

inside = .false.
If ( Abs(Sum(theta)) .gt. pi ) Then
  inside = .true.
Endif

Endfunction inside_poly


EndModule math_routines_mod

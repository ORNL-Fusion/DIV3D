!-----------------------------------------------------------------------------
!+ Contains special math functions for w7 routines
!-----------------------------------------------------------------------------
Module math_routines_mod
!
! Description:
!
! History:
! Version   Date      Comment
! -------   ----      -------
!  1.0     07/18/2011  Ported from PENTA  JDL
!  1.1     12/01/2011  Added int_two_lines JDL
! Author(s): J. Lore 07/2009 - xxx 
!
! Modules used:
Use kind_mod                ! Import real64, int32 specifications       

Implicit None

Contains


!-----------------------------------------------------------------------------
!+ Checks for the intersection of two line segments
!-----------------------------------------------------------------------------
Subroutine int_two_lines(p1,p2,p3,p4,u1,u2)
!; Calculates the intersection of lines p1 to p2, and p3 to p4.  Each
!; point is an x-y pair.  Returns two values, which are the distance
!; between p1 and p2 of the intersection, and between p3 and p4.  This
!; distance is normalized to the length of the lines.

Use kind_mod
implicit none

Real(real64), Intent(in),Dimension(2) :: p1,p2,p3,p4
Real(real64), intent(out) :: u1,u2
Real(real64) :: denom

denom = (p4(2)-p3(2))*(p2(1)-p1(1)) - (p4(1)-p3(1))*(p2(2)-p1(2))

if ( denom .eq. 0._real64 ) then 
  u1 = 1.d30
  u2 = 1.d30
else
    u1 = ((p4(1)-p3(1))*(p1(2)-p3(2)) - (p4(2)-p3(2))*(p1(1)-p3(1)))/denom
    u2 = ((p2(1)-p1(1))*(p1(2)-p3(2)) - (p2(2)-p1(2))*(p1(1)-p3(1)))/denom
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
! Inputs:
!  xstart,xend: Values of the first and last points of the array [real]
!  numel: Number of elements in the array [integer]
! Outputs:
!  rlinvec: The linearly spaced array
!         
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/22/2009  Original Code.  JL
!  1.1     09/01/2010  Updated for PENTA3. JL
!  1.2     07/18/2011  Ported to w7 routines. JL
!
! Author(s): J. Lore 07/2009 - 7/18/2011
!
! Modules used:
Use kind_mod                ! Import real64, int32 specifications

Implicit None

! Input/output                      !See above for descriptions
Real(real64),    Intent(in) :: xstart  
Real(real64),    Intent(in) :: xend
Integer(int32), Intent(in) :: numel
Real(real64)                :: rlinvec(numel)

! Local scalars
Integer(int32)   ::  ii
!- End of header -------------------------------------------------------------

Do ii = 1,numel
  rlinvec(ii) = ( Real(ii,real64) - 1._real64 ) * ( xend - xstart ) &
       / ( numel - 1._real64 ) + xstart
Enddo

EndFunction rlinspace


!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------
Subroutine line_seg_facet_int(pa,pb,pc,p1,p2,ithit,p,tol)
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
Use kind_mod

Implicit none

! Input/output
Real(real64), Dimension(3), Intent(In) :: pa,pb,pc,p1,p2
Real(real64), Dimension(3), Intent(Out) :: p
Integer(int32), Intent(Out) :: ithit
real(real64), Intent(In) :: tol


!local variables and arrays
Real(real64), Dimension(3) :: n, pa1, pa2, pa3
Real(real64) :: d, denom, mu, total, a1, a2, a3

! Local Parameters
Real(real64), Parameter :: pi = 3.141592653589793238462643383279502_real64
!- End of header -------------------------------------------------------------


!tol = 1.d-12

! Calc unit vector normal to plane of Pa-c (gives plane components A-C)
n(1) = (pb(2) - pa(2))*(pc(3) - pa(3)) - (pb(3) - pa(3))*(pc(2) - pa(2))
n(2) = (pb(3) - pa(3))*(pc(1) - pa(1)) - (pb(1) - pa(1))*(pc(3) - pa(3))
n(3) = (pb(1) - pa(1))*(pc(2) - pa(2)) - (pb(2) - pa(2))*(pc(1) - pa(1))
n=n/Sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3))

! Calculate plane component D
d = - n(1)*pa(1) - n(2)*pa(2) - n(3)*pa(3)

! Calculate the position on the line that intersects the plane
denom = n(1)*(p2(1) - p1(1)) + n(2)*(p2(2) - p1(2)) + n(3)*(p2(3) - p1(3))

p(:) = 0._real64

If (abs(denom) .lt. tol) Then
    ithit = 0
else    
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



End Subroutine line_seg_facet_int



EndModule math_routines_mod

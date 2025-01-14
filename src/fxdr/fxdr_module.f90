module fxdr_module
implicit none

integer, parameter :: FXDRERR_WRNEGNELS=-10, FXDRERR_WRITEERR=-11
integer, parameter :: FXDRERR_READERR=-12, FXDRERR_READWRONGNELS=-13
integer, parameter :: FXDRERR_REWIND=-14

contains

!!       -------------------------------------------------
!!       Returns < 0 on error, otherwise the xdrid used to
!!       access the xdr file 
!!
!!       Args:
!!           -   filename: character*(*): name of file to 
!!                       access
!!           -   mode: character*1: either 'r', 'w', or 'a' 
!!                       (for read, write, and append)
!!           -   returnonerror: logical: if .TRUE., then
!!                       routines return even if there is
!!                       an I/O error.  If .FALSE., then
!!                       routines halt on I/O error, 
!!                       printing out a useful message 
!!                       (this is like what fortran does,
!!                       so set to .FALSE. for most fortran-
!!                       like behavior)
!!       -------------------------------------------------
function initxdr( filename, mode, returnonerror ) result(ixdr)

implicit none

character(len=*), intent(in) :: filename
character(1), intent(in) :: mode
logical, intent(in) :: returnonerror
integer :: ixdr

integer :: iretval, cxdrinit, iroe, imode, lenf, istrlen

select case(mode)
 case("r","R")
  imode = 1
 case("w","W")
  imode = 2
 case("a","A")
  imode = 3
 case default
  write(0,*) 'Error: fxdr library, initxdr called '
  write(0,*) 'with unknown mode (should be r, w, or a):', mode
  stop 'fxdr library: initxdr routine, bad mode'
end select

lenf    = len(filename)
istrlen = ifxdrstrlen( filename )
if( istrlen .lt. lenf ) lenf = istrlen

if( returnonerror ) then
 iroe = 1
else
 iroe = 0
endif

iretval = cxdrinit( lenf, filename, imode, iroe )
if( iretval .lt. 0 ) then
 ixdr = iretval
 return
endif

!-----------------------------------------------------------
! Add one to the returned value because the C library returns 
! starting with zero, which is awkward considering that 
! unitialized variables often start with zero in Fortran.
!-----------------------------------------------------------
ixdr = iretval + 1

return

contains

!========================================================================

function ifxdrstrlen( s ) result(ifxdrsl)
implicit none

character(len=*), intent(in) :: s
integer :: ifxdrsl
integer :: i

i = 1
do
 if( s(i:i) .eq. ' ' ) exit
 i = i + 1
 if( i .gt. len(s) ) exit 
enddo
ifxdrsl = i-1

return
end function ifxdrstrlen

end function initxdr
!! -----------------------------
!! Returns < 0 if error, 0 if OK
!! -----------------------------
function ixdrclose( ixdr ) result(ierr)

implicit none

integer, intent(in) :: ixdr
integer :: ixdrm1, ierr

ixdrm1 = ixdr - 1
call cxdrclose( ixdrm1, ierr )

return
end function ixdrclose

!========================================================================
!! -----------------------------
!! Returns < 0 if error, 0 if OK
!! -----------------------------
function ixdrdmat( ixdrs, nels, dval ) result(ierr)

implicit none

integer, intent(in) :: ixdrs, nels
real(kind(1.0d0)), dimension(nels), intent(inout) :: dval(nels)
integer :: ixdrsm1, ierr

ixdrsm1 = ixdrs - 1
call cxdrdmat( ixdrsm1, nels, dval, ierr )

return
end function ixdrdmat

!! -----------------------------
!! Returns < 0 if error, 0 if OK
!! -----------------------------
function ixdrdouble( ixdrs, dval ) result(ierr)

implicit none

integer, intent(in) :: ixdrs
real(kind(1.0d0)), intent(inout) :: dval
integer :: ixdrsm1, ierr

ixdrsm1 = ixdrs - 1

!print *, "Got :", dval, " for processing"
call cxdrdouble( ixdrsm1, dval, ierr )
!print *, "Have:", dval, " for processing"

return
end function ixdrdouble


!! -----------------------------
!! Returns < 0 if error, 0 if OK
!! -----------------------------
function ixdrimat( ixdrs, nels, ival ) result(ierr)

implicit none

integer, intent(in) :: ixdrs, nels
integer, dimension(nels), intent(inout) :: ival(nels)
integer :: ixdrsm1, ierr

ixdrsm1 = ixdrs - 1
call cxdrimat( ixdrsm1, nels, ival, ierr )

return
end function ixdrimat

!! -----------------------------
!! Returns < 0 if error, 0 if OK
!! -----------------------------
function ixdrint( ixdrs, ival ) result(ierr)

implicit none

integer, intent(in) :: ixdrs
integer, intent(inout) :: ival
integer :: ixdrsm1, ierr

ixdrsm1 = ixdrs - 1
call cxdrint( ixdrsm1, ival, ierr )

return
end function ixdrint

!! -----------------------------
!! Returns < 0 if error, 0 if OK
!! -----------------------------
function ixdrreal64( ixdrs, rval ) result(ierr)
use, intrinsic :: iso_fortran_env, only: real64
implicit none

integer, intent(in) :: ixdrs
real(real64), intent(inout) :: rval
integer :: ixdrsm1, ierr

ixdrsm1 = ixdrs - 1

call cxdrreal64( ixdrsm1, rval, ierr )

return
end function ixdrreal64

!! -----------------------------
!! Returns < 0 if error, 0 if OK
!! -----------------------------
function ixdrreal( ixdrs, rval ) result(ierr)

implicit none

integer, intent(in)    :: ixdrs
real(kind(1.0e0)),    intent(inout) :: rval
integer                :: ixdrsm1, ierr

ixdrsm1 = ixdrs - 1

call cxdrreal( ixdrsm1, rval, ierr )

return
end function ixdrreal

!! -----------------------------
!! Returns < 0 if error, 0 if OK
!! -----------------------------
function ixdrrewind( ixdrs ) result(ierr)

implicit none

integer, intent(in) :: ixdrs
integer :: ixdrsm1, ierr

ixdrsm1 = ixdrs - 1

call cxdrrewind( ixdrsm1, ierr )

return
end function ixdrrewind

!! -----------------------------
!! Returns < 0 if error, 0 if OK
!! -----------------------------
function ixdrrmat64( ixdrs, nels, rval ) result(ierr)
use, intrinsic :: iso_fortran_env, only: real64
implicit none

integer, intent(in) :: ixdrs, nels
real(real64), dimension(nels), intent(inout) :: rval
integer :: ixdrsm1, ierr

ixdrsm1 = ixdrs - 1

call cxdrrmat64( ixdrsm1, nels, rval, ierr )

return
end function ixdrrmat64

!! -----------------------------
!! Returns < 0 if error, 0 if OK
!! -----------------------------
function ixdrrmat( ixdrs, nels, rval ) result(ierr)

implicit none

integer, intent(in) :: ixdrs, nels
real(kind(1.0e0)), dimension(nels), intent(inout) :: rval
integer :: ixdrsm1, ierr

ixdrsm1 = ixdrs - 1

call cxdrrmat( ixdrsm1, nels, rval, ierr )

return
end function ixdrrmat

!! -----------------------------
!! Returns < 0 if error, 0 if OK
!! -----------------------------
function ixdrshort (ixdrs, num) result(ierr)

implicit none

integer, intent(in) :: ixdrs,num
integer :: ixdrsm1,ierr

ixdrsm1=ixdrs - 1

call cxdrshort( ixdrsm1, num, ierr )

return
end function ixdrshort

!! -----------------------------
!! Returns < 0 if error, 0 if OK
!! -----------------------------
function ixdrstring( ixdrs, string ) result(ierr)

implicit none

integer, intent(in) :: ixdrs
character(len=*), intent(inout) :: string
integer :: ixdrsm1, ierr, ilen

ixdrsm1 = ixdrs - 1
        
ilen = len(string)
call cxdrstring( ixdrsm1, ilen, string, ierr )

return
end function ixdrstring

end module fxdr_module

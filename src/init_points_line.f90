!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------
Subroutine init_points_line(fname_surf,numl,fname_launch)
!
! Description: 
!
! Inputs: 
!  numl -- number of points to define along line
!
! Outputs:
!   none
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/26/2011   JDL
! Author(s): J.D. Lore - 07/26/2011 - xxx

! Modules used:
Use kind_mod
Use io_unit_spec, Only: &
iu_launch, iu_surf
Implicit none

! Input/output
Character(len=100),Intent(in) :: fname_surf, fname_launch
Integer(int32), Intent(in) :: numl

! Local scalars
Integer(int32) :: ii, rand_ind
Integer(int32) :: npts_line,nip0, ip_step
Real(real64) :: period, P1, rnum

! Local arrays
Real(real64),Allocatable,Dimension(:) :: rsurf,zsurf,phisurf

! Local parameters
Real(real64), Parameter :: pi = 3.141592653589793238462643383279502_real64

!- End of header -------------------------------------------------------------

! Read surface file
Write(6,*) 'Reading surface data from ',Trim(Adjustl(fname_surf))
Open(iu_surf,file=fname_surf)
Read(iu_surf,*) period, nip0, ip_step
Read(iu_surf,*) npts_line
Allocate(rsurf(npts_line),zsurf(npts_line),phisurf(npts_line))
Do ii = 1,npts_line
  Read(iu_surf,*) rsurf(ii)
  Read(iu_surf,*) zsurf(ii)
  Read(iu_surf,*) phisurf(ii)
Enddo
Close(iu_surf)

Write(6,*) 'Writing launch point data to ',Trim(Adjustl(fname_launch))
Open(iu_launch,file=fname_launch)
Write(iu_launch,*) numl

Do ii = 1,numl
  ! Choose an integer between 1 and npts_line
!  rand_ind = Nint(npts_line*rand())
  Call Random_number(rnum)
  rand_ind = Nint(npts_line*rnum)
  P1 = phisurf(rand_ind)
  Do While (P1 .lt. 0.d0)
    P1 = P1 + period
  Enddo
  P1 = Mod(P1,period)
  Write(iu_launch,*) rsurf(rand_ind),zsurf(rand_ind),P1
Enddo
Close(iu_launch)

Endsubroutine init_points_line
!-----------------------------------------------------------------------------

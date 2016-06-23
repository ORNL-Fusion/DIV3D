!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------
Subroutine trace_surface(Rstart,Zstart,Phistart,dphi_line,nsteps_line,period,fname_surf,div3d_bfield_method)
!
! Description: 
!
! Inputs: 
!   Rstart,Zstart: (meters)
!   Phistart: (radians)
!   dphi_line: Integration step size in phi (radians)
!   nsteps_line: Number of integration steps
!   period: Toroidal extent of field period (degrees)
!   fname_surf: Output file name
!   
! Outputs:
!   none
!
! Surface data:
!   rsurf,zsurf: Array(nsteps_line+1) (meters)
!   phisurf: Array(nsteps_line+1) (radians)
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/15/2011   JDL
! Author(s): J.D. Lore - 07/15/2011 - xxx

! Modules used:
Use kind_mod
Use io_unit_spec, Only: &
iu_surf
Implicit none

! Input/output
Real(rknd), Intent(in) :: Rstart, Zstart, Phistart
Real(rknd), Intent(in) :: dphi_line, Period
Integer(iknd), Intent(in) :: nsteps_line
Integer(iknd), Intent(in) :: div3d_bfield_method
Character(len=100), Intent(in) :: fname_surf

! Local scalars
Integer(iknd) :: ifail, ii, ip_step, nip0
Real(rknd) :: adp

! Local arrays
Real(rknd), Dimension(nsteps_line+1) :: &
  rsurf,zsurf,phisurf

! Local parameters
Real(rknd), Parameter :: pi = 3.141592653589793238462643383279502_rknd

!- End of header -------------------------------------------------------------


! Check periodicity

adp = Abs(dphi_line)
nip0 = floor( nsteps_line*adp/period) + 1
ip_step = Nint(period/adp)
If (Real(period/adp,rknd) - Real(ip_step,rknd) .gt. 1.d-12 ) Then
  Write(6,*) 'Choose dphi_line such that ',period*180.d0/pi,' degrees is divisible'
  Write(*,*) 'Dphi_line (deg):',dphi_line*180.d0/pi
  Write(6,*) 'Should be equal:',Real(period/adp,rknd),ip_step,period/adp
  Write(6,*) 'Periodicity, (ip_step,nip0) = ', ip_step, nip0
  Stop
Else
  Write(6,*) 'Periodicity ok, (ip_step,nip0) = ', ip_step, nip0
Endif

! Trace out surface
!Write(6,*) 'Following fieldline from [R,Z,phi] = ' &
!  ,Rstart,Zstart,phistart*180./pi
Call follow_fieldline_rzphi(Rstart,Zstart,Phistart,dphi_line,nsteps_line,rsurf,zsurf,phisurf,.false.,0.d0,ifail&
     ,div3d_bfield_method)

If (ifail .ne. 0 ) Then
  Write(6,*) 'fieldline following error',ifail
  Stop
Endif

! Write line data
Write(6,*) 'Writing surface data to ',Trim(Adjustl(fname_surf))
Open(iu_surf,file=fname_surf)
Write(iu_surf,*) period, nip0, ip_step
Write(iu_surf,*) nsteps_line + 1
Do ii = 1,nsteps_line + 1 
  Write(iu_surf,*) rsurf(ii)
  Write(iu_surf,*) zsurf(ii)
  Write(iu_surf,*) phisurf(ii)
Enddo
Close(iu_surf)


Endsubroutine trace_surface

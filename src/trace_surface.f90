!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------
Subroutine trace_surface(Rstart,Zstart,Phistart,dphi_line,nsteps_line,period,fname_surf)
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
Use io_unit_spec, Only : iu_surf
Use fieldline_follow_mod, Only : follow_fieldlines_rzphi
Use setup_bfield_module, Only : bfield
Use parallel_mod
Use phys_const, Only : pi
Implicit none

! Input/output
Real(real64), Intent(in) :: Rstart, Zstart, Phistart
Real(real64), Intent(in) :: dphi_line, Period
Integer(int32), Intent(in) :: nsteps_line
Character(len=100), Intent(in) :: fname_surf

! Local scalars
Integer(int32) :: ifail, ii, ip_step, nip0, ierr
Real(real64) :: adp

! Local arrays
Real(real64), Dimension(nsteps_line+1) :: &
  rsurf,zsurf,phisurf

!- End of header -------------------------------------------------------------


! Check periodicity

adp = Abs(dphi_line)
nip0 = floor( nsteps_line*adp/period) + 1
ip_step = Nint(period/adp)
If (abs(Real(period/adp,real64) - Real(ip_step,real64)) .gt. 1.d-12 ) Then
  Write(*,*) 'Choose dphi_line such that ',period*180.d0/pi,' degrees is divisible'
  Write(*,*) 'Dphi_line (deg):',dphi_line*180.d0/pi
  Write(*,*) 'Should be equal:',Real(period/adp,real64),ip_step,period/adp
  Write(*,*) 'Periodicity, (ip_step,nip0) = ', ip_step, nip0
  Stop
Else
  Write(*,*) 'Periodicity ok. (ip_step,nip0) = ', ip_step, nip0
Endif

! Trace out surface
!Write(6,*) 'Following fieldline from [R,Z,phi] = ',Rstart,Zstart,phistart*180./pi
Call follow_fieldlines_rzphi(bfield,Rstart,Zstart,Phistart,dphi_line,nsteps_line,rsurf,zsurf,phisurf,ierr,ifail)
If (ierr .ne. 0 ) Then
  Write(*,*) 'fieldline following error',ierr
  Call flush()
  Stop
Endif

! Write line data
Write(*,*) 'Writing surface data to ',Trim(Adjustl(fname_surf))
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

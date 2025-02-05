Module surface_mod
  Implicit None
  Public :: trace_surface

Contains
!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------
Subroutine trace_surface(nsteps_line)
!
!
! Author(s): J.D. Lore - 07/15/2011 - xxx
Use kind_mod, Only : real64, int32
Use io_unit_spec, Only : iu_surf
Use fieldline_follow_mod, Only : follow_fieldlines_rzphi
Use setup_bfield_module, Only : bfield, rmp_type
Use g3d_module, Only : get_psiN_bicub
Use parallel_mod, Only : fin_mpi
Use phys_const, Only : pi
Use run_settings_namelist, Only : Rstart, Zstart, Phistart, dphi_line_surf, &
     period, fname_surf
Implicit none

Integer(int32), Intent(in) :: nsteps_line
Integer(int32) :: ifail, ii, ip_step, nip0, ierr
Real(real64) :: adp
Real(real64), Dimension(nsteps_line+1) :: rsurf,zsurf,phisurf
Real(real64) :: psiN(1)

!- End of header -------------------------------------------------------------

Write(*,'(/A,3(F8.2))') ' Tracing initial surface from (R,Z,Phi) = ',Rstart,Zstart,Phistart

Select Case (rmp_type)
Case ('g')
   Call get_psiN_bicub(bfield%g,(/Rstart/),(/Zstart/),1,psiN,ierr)
   Write(*,*) "Axis position [R,Z]",bfield%g%rmaxis,bfield%g%zmaxis
   Write(*,*) "Psi_N of start point",psiN(1)
End Select

! Check periodicity
adp = Abs(dphi_line_surf)
nip0 = floor( nsteps_line*adp/period) + 1
ip_step = Nint(period/adp)
If (abs(Real(period/adp,real64) - Real(ip_step,real64)) .gt. 1.d-12 ) Then
  Write(*,*) 'Choose dphi_line such that ',period*180.d0/pi,' degrees is divisible'
  Write(*,*) 'Dphi_line (deg):',dphi_line_surf*180.d0/pi
  Write(*,*) 'Should be equal:',Real(period/adp,real64),ip_step,period/adp
  Write(*,*) 'Periodicity, (ip_step,nip0) = ', ip_step, nip0
  Stop
Else
  Write(*,*) 'Periodicity ok. (ip_step,nip0) = ', ip_step, nip0
Endif

! Trace out surface
!Write(6,*) 'Following fieldline from [R,Z,phi] = ',Rstart,Zstart,phistart*180./pi
Call follow_fieldlines_rzphi(bfield,Rstart,Zstart,Phistart,dphi_line_surf,nsteps_line,rsurf,zsurf,phisurf,ierr,ifail)
If (ierr .ne. 0 ) Then
  Write(*,*) 'Error: fieldline following error',ierr
  Call fin_mpi(.true.)
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

End Subroutine trace_surface

End Module surface_mod

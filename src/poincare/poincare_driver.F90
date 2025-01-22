!-----------------------------------------------------------------------------
!+ Driver for poincare in DIV3D
!-----------------------------------------------------------------------------
Program poincare_driver
  !
  ! Author(s): J.D. Lore 
  !
  Use kind_mod, Only: int32, real64
  Use poincare_namelist, Only : read_poincare_namelist
  Use initialize_bfield_div3d, Only : init_bfield
  Use parallel_mod, Only : rank, nprocs, init_mpi, fin_mpi
  Use timing_mod, Only : get_elapsed_time, init_timing  
  !Use g3d_module, Only : get_psi_bicub
  !Use math_geo_module, Only : rlinspace
  !Use fieldline_follow_mod, Only: follow_fieldlines_rzphi
  !Use util_routines, Only: get_psin_2d
  !Use phys_const, Only: pi
  !Use setup_bfield_module

  Implicit none

  Integer :: tstart ! Timing variable
  Integer(int32) :: iocheck
  Logical :: verbose = .false.  

  !  Integer(int32) :: i, nstart_fl, nsteps
  !  Integer(int32) :: iocheck, ind_poin, ierr
  !  Real(real64), Allocatable :: r1d(:), z1d(:),phistart_arr(:), fl_r(:,:), fl_z(:,:), fl_p(:,:), &
  !       psiout(:), psiNout(:), fl_r2(:,:), fl_z2(:,:), fl_p2(:,:)
  !  Integer(int32), Allocatable :: ilg(:), fl_ierr(:), ilg2(:), fl_ierr2(:)
  !  Real(real64) :: dphi_line, Adphirat
  !  Real(Kind=4)  :: tarray(2),tres,tres0

  
  !---------------------------------------------------------------------------
  ! 0. Initialization
  !---------------------------------------------------------------------------  
  Call init_mpi
  Call init_timing
  
  If (rank .eq. 0) Then
     call system_clock(tstart)
     verbose = .true.
     Write(*,'(/A)') '-------------------------------------------------------------------------'
     Write(*,'(a)') " Starting Poincare driver"     
  End If
  
  ! Read namelists
  Call read_poincare_namelist(verbose)
  
  !----------------------------------------------------------
  ! 1. Initialize magnetic field
  !----------------------------------------------------------
  Call init_bfield(verbose)


  ! !
  ! ! 2) Magnetic equilibrium data
  ! !
  ! ! Setup rmp field
  ! Select Case (rmp_type)
  ! Case ('g3d')
  !    Call setup_bfield_g3d
  ! Case ('g3d+rmpcoil')
  !    Call setup_bfield_g_and_rmp
  ! Case ('g3d+m3dc1')
  !    Call setup_bfield_g_and_m3dc1
  ! Case ('m3dc1_full_field')
  !    Call setup_bfield_m3dc1_full
  ! Case ('m3dc1_as')
  !    Call setup_bfield_m3dc1_as
  ! Case ('ipec')
  !    Call setup_bfield_ipec
  ! Case ('xpand')
  !    Call setup_bfield_xpand
  ! Case ('vmec_coils')
  !    Call setup_bfield_vmec_coils
  ! Case ('vmec_coils_to_fil')
  !    Call setup_bfield_vmec_coils_to_fil
  ! Case ('bgrid')
  !    Call setup_bfield_bgrid
  ! Case ('xdr')
  !    Call setup_bfield_xdr    
  ! Case Default
  !    Write(*,*) 'Unknown rmp_type in poincare_driver!'
  !    Write(*,*) 'Current options are:'
  !    Write(*,*) '''g3d'''
  !    Write(*,*) '''g3d+rmpcoil'''
  !    Write(*,*) '''g3d+m3dc1'''
  !    Write(*,*) '''m3dc1_full_field'''
  !    Write(*,*) '''m3dc1_as'''
  !    Write(*,*) '''ipec'''
  !    Write(*,*) '''xpand'''
  !    Write(*,*) '''vmec_coils'''
  !    Write(*,*) '''vmec_coils_to_fil'''
  !    Write(*,*) '''bgrid'''
  !    Write(*,*) '''xdr'''    
  !    Stop      
  ! End Select


  ! ! 
  ! ! 3) Follow fls
  ! !

  ! Allocate(r1d(num_pts),z1d(num_pts))
  ! r1d = rlinspace(rstart,rend,num_pts)
  ! z1d = rlinspace(zstart,zend,num_pts)

  ! dphi_line = dphi_line_deg*pi/180.d0
  ! nsteps = Floor(ntransits*2.d0*pi/Abs(dphi_line))

  ! Adphirat = Abs(360.d0/dphi_line_deg/Nsym)
  ! If (Adphirat - Real(Nint(Adphirat)) > 1.d-8) Then
  !    Write(*,*) 'Error!: 2*pi/Nsym must be an integer multiple of dphi_line_deg'
  !    Write(*,*) 'Exiting'
  !    Stop
  ! Endif
  ! ind_poin = Nint(Adphirat)


  ! Write(*,'(/1a,i0,a,i0,a)') 'Following ',num_pts,' fls for ',ntransits,' transits.'
  ! Write(*,'(a,f12.3,a)') 'Toroidal step size is ',dphi_line_deg,' degrees'
  ! Write(*,'(a,i0)') 'Poincare plot step size index is ',ind_poin
  ! Write(*,'(a,f12.3)') 'Starting fl at phi = ',phistart_deg
  ! Write(*,'(a,i0)') 'Number of steps = ',nsteps

  ! nstart_fl = num_pts
  ! Allocate(ilg(nstart_fl),fl_ierr(nstart_fl),phistart_arr(nstart_fl))
  ! Allocate(fl_r(nstart_fl,nsteps+1),fl_z(nstart_fl,nsteps+1),fl_p(nstart_fl,nsteps+1))
  ! fl_r = 0.d0; fl_z = 0.d0; fl_p = 0.d0
  ! phistart_arr = phistart_deg*pi/180.d0

  ! Call follow_fieldlines_rzphi(bfield,r1d,z1d,phistart_arr,nstart_fl, dphi_line,nsteps,fl_r,fl_z,fl_p,fl_ierr,ilg)

  ! If (follow_both_ways) Then
  !    Allocate(ilg2(nstart_fl),fl_ierr2(nstart_fl))
  !    Allocate(fl_r2(nstart_fl,nsteps+1),fl_z2(nstart_fl,nsteps+1),fl_p2(nstart_fl,nsteps+1))
  !    fl_r2 = 0.d0; fl_z2 = 0.d0; fl_p2 = 0.d0 
  !    Call follow_fieldlines_rzphi(bfield,r1d,z1d,phistart_arr,nstart_fl,-dphi_line,nsteps,fl_r2,fl_z2,fl_p2,fl_ierr2,ilg2)
  ! Endif

  ! Open(99,file="poincare_output.out",status="unknown",form="formatted",iostat=iocheck)
  ! If ( iocheck /= 0 ) Then
  !    Write(*,*) 'Error opening output file'
  !    Stop 'Exiting: I/O Error in poincare_driver.f90'
  ! Endif
  ! Write(99,*) phistart_deg,nstart_fl,nsteps/ind_poin+1
  ! !Write(*,*) phistart_deg,nstart_fl,nsteps/ind_poin+1

  ! Do i = 1,nstart_fl
  !    Write(99,*) i
  !    Write(99,'(6e20.12)') fl_r(i,1:nsteps+1:ind_poin)
  !    Write(99,'(6e20.12)') fl_z(i,1:nsteps+1:ind_poin)
  ! Enddo
  ! Close(99)
  ! Deallocate(r1d,z1d)

  ! If (follow_both_ways) Then
  !    Open(99,file="poincare_output2.out",status="unknown",form="formatted",iostat=iocheck)
  !    If ( iocheck /= 0 ) Then
  !       Write(*,*) 'Error opening output file'
  !       Stop 'Exiting: I/O Error in poincare_driver.f90'
  !    Endif
  !    Write(99,*) phistart_deg,nstart_fl,nsteps/ind_poin+1

  !    Do i = 1,nstart_fl
  !       Write(99,*) i
  !       Write(99,'(6e20.12)') fl_r2(i,1:nsteps+1:ind_poin)
  !       Write(99,'(6e20.12)') fl_z2(i,1:nsteps+1:ind_poin)
  !    Enddo
  !    Close(99)
  ! Endif

  ! !
  ! ! Calculate min psi_N

  ! !
  ! If (calc_psiN_min) Then
  !    Write(*,*) 'Calculating minimum psi_N'
  !    Allocate(psiout(nsteps+1),psiNout(nsteps+1))

  !    Open(99,file="psiN_min_output.out",status="unknown",form="formatted",iostat=iocheck)
  !    If ( iocheck /= 0 ) Then
  !       Write(*,*) 'Error opening output file'
  !       Stop 'Exiting: I/O Error in poincare_driver.f90 (2)'
  !    Endif
  !    Write(99,*) nstart_fl

  !    Do i = 1,nstart_fl

  !       !Call get_psi_bicub(fl_r(i,:),fl_z(i,:),nsteps+1,psiout,psiNout,ierr)
  !       psiNout = get_psiN_2d(bfield,fl_r(i,:),fl_z(i,:),nsteps+1,ierr)

  !       Where (psiNout < 1.e-3) psiNout = 1000000.d0
  !       Write(99,*) Minval(psiNout)
  !    Enddo

  !    Close(99)
  !    Deallocate(psiout,psiNout)
  ! Endif

  ! If (follow_both_ways) Then
  !    ! Calculate min psi_N
  !    If (calc_psiN_min) Then
  !       Write(*,*) 'Calculating minimum psi_N'
  !       Allocate(psiout(nsteps+1),psiNout(nsteps+1))

  !       Open(99,file="psiN_min_output2.out",status="unknown",form="formatted",iostat=iocheck)
  !       If ( iocheck /= 0 ) Then
  !          Write(*,*) 'Error opening output file'
  !          Stop 'Exiting: I/O Error in poincare_driver.f90 (2)'
  !       Endif
  !       Write(99,*) nstart_fl

  !       Do i = 1,nstart_fl
  !          psiNout = get_psiN_2d(bfield,fl_r2(i,:),fl_z2(i,:),nsteps+1,ierr)
  !          Where (psiNout < 1.e-3) psiNout = 1000000.d0
  !          Write(99,*) Minval(psiNout)
  !       Enddo

  !       Close(99)
  !       Deallocate(psiout,psiNout)
  !    Endif
  ! Endif


  ! Deallocate(ilg,fl_ierr,fl_r,fl_z,fl_p,phistart_arr)

  ! If (follow_both_ways) Then
  !    Deallocate(ilg2,fl_ierr2,fl_r2,fl_z2,fl_p2)
  ! Endif

  ! Call Etime(tarray,tres)
  ! Write(*,*) ' Poincare_driver took ',tres-tres0,' seconds'


End Program poincare_driver






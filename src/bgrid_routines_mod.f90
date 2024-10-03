!-----------------------------------------------------------------------------
!
!   Routines/modules related to bgrid evaluation
!
!
!   Contains:
!-----------------------------------------------------------------------------
Module bgrid_module
Use kind_mod, Only: real64, int32
Implicit None

Real(real64), Allocatable, Save :: bgrid_r(:), bgrid_z(:), bgrid_phi(:)
Real(real64), Allocatable, Dimension(:,:,:), Save :: &
     bgrid_br,bgrid_bz, bgrid_bphi
Integer(int32), Save :: bgrid_nr, bgrid_nz, bgrid_nphi, nsym

Contains

!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------
Subroutine read_bgrid_field_file(fname_base,verbose)
  ! Author(s): J.D. Lore
  Use kind_mod, Only: int32
  Use phys_const, Only: pi
  Use parallel_mod, Only : fin_mpi
  Implicit None
  Character(len=*), Intent(In) :: fname_base
  Logical, Intent(In) :: verbose
  Integer(int32) :: iocheck, nphi,nr,nz, i, j, k
  Character(len=1000) :: fname
  Real(real64), Allocatable :: data(:,:,:)
  Real(real64) :: R_min, R_max, Z_min, Z_max
  !- End of header -------------------------------------------------------------

  fname = Trim(fname_base)//'_layout.dat'
  If (verbose) Write(*,'(a,a)') ' Reading ', Trim(fname)
  open(99,file=fname,IOSTAT=iocheck,status='old')
  If (iocheck /= 0) Then
     Write(*,*) 'Error opening file: ',Trim(fname)
     Call fin_mpi(.true.)
  Endif

  read(99, *) nr, nz, nphi, nsym, R_min, R_max, Z_min, Z_max
  close(99)
  If (verbose) Then
     write(*,*) 'nr   = ',nr
     write(*,*) 'nphi = ',nphi
     write(*,*) 'nz   = ',nz
     write(*,*) 'nsym = ',nsym
     write(*,*) 'R_(min,max) =',R_min, R_max
     write(*,*) 'Z_(min,max) =',Z_min, Z_max
  Endif

  nphi = nphi + 1 ! Will add first slice to 2*pi

  bgrid_nr = nr
  bgrid_nz = nz
  bgrid_nphi = nphi

  Allocate(bgrid_r(nr))
  Allocate(bgrid_z(nz))
  Allocate(bgrid_phi(nphi))

  do i=1,nr
     bgrid_r(i) = R_min + 1.d0*(i-1)/(nr-1) * (R_max - R_min)
  enddo
  do i=1,nz
     bgrid_z(i) = Z_min + 1.d0*(i-1)/(nz-1) * (Z_max - Z_min)
  enddo
  do i=1,nphi
     bgrid_phi(i) = 2.d0*pi/nsym*(i-1)/(nphi-1)
  enddo
  Allocate(bgrid_br(nr,nz,nphi))
  Allocate(bgrid_bz(nr,nz,nphi))
  Allocate(bgrid_bphi(nr,nz,nphi))

  bgrid_br = 0.d0
  bgrid_bz = 0.d0
  bgrid_bphi = 0.d0

  Allocate(data(nr,nz,nphi))

  fname = Trim(fname_base)//'_r.dat'
  If (verbose) Write(*,'(a,a)') ' Reading ', Trim(fname)
  open(99,file=fname,IOSTAT=iocheck,status='old')
  If (iocheck /= 0) Then
     Write(*,*) 'Error opening file: ',Trim(fname)
     Call fin_mpi(.true.)
  Endif
  read(99, *) (((data(i,j,k), k=1,nphi-1), j=1,nz), i=1,nr)
  close(99)
  bgrid_br = data

  fname = Trim(fname_base)//'_z.dat'
  If (verbose) Write(*,'(a,a)') ' Reading ', Trim(fname)
  open(99,file=fname,IOSTAT=iocheck,status='old')
  If (iocheck /= 0) Then
     Write(*,*) 'Error opening file: ',Trim(fname)
     Call fin_mpi(.true.)
  Endif
  read(99, *) (((data(i,j,k), k=1,nphi-1), j=1,nz), i=1,nr)
  close(99)
  bgrid_bz = data

  fname = Trim(fname_base)//'_phi.dat'
  If (verbose)   Write(*,'(a,a)') ' Reading ', Trim(fname)
  open(99,file=fname,IOSTAT=iocheck,status='old')
  If (iocheck /= 0) Then
     Write(*,*) 'Error opening file: ',Trim(fname)
     Call fin_mpi(.true.)
  Endif
  read(99, *) (((data(i,j,k), k=1,nphi-1), j=1,nz), i=1,nr)
  close(99)
  bgrid_bphi = data

  Deallocate(data)

  bgrid_br(:,:,nphi)   = bgrid_br(:,:,1)
  bgrid_bz(:,:,nphi)   = bgrid_bz(:,:,1)
  bgrid_bphi(:,:,nphi) = bgrid_bphi(:,:,1)

End Subroutine read_bgrid_field_file

Subroutine open_bgrid_fields(fname,verbose)
  Implicit None
  Character(Len=*), Intent(In) :: fname
  Logical, Intent(In) :: verbose
  Call read_bgrid_field_file(fname,verbose)
End Subroutine open_bgrid_fields

Subroutine close_bgrid_fields
  Implicit None
  Deallocate(bgrid_r,bgrid_z,bgrid_phi)
  Deallocate(bgrid_br,bgrid_bz,bgrid_bphi)
End Subroutine close_bgrid_fields

!-----------------------------------------------------------------------------
!+ Evaluate B(r,phi,z) using Equilibrium only BGRID fields
!-----------------------------------------------------------------------------
Subroutine bfield_bgrid(r,phi,z,Npts,Bout,ierr)
  !   Bout = (:,[Br,Bz,Bt])
  Use kind_mod, Only: int32, real64
  Use phys_const, Only: pi
  Use parallel_mod, Only : fin_mpi
  Implicit None
  Real(Real64), Intent(In), Dimension(Npts) :: r, z, phi
  Integer(int32), Intent(In) :: Npts
  Real(real64), Intent(Out), Dimension(Npts,3) :: Bout
  Integer(int32), Intent(Out) :: ierr
  ! Local variables
  Real(real64) :: phi_tmp, phi_fac, dphi_grid
  Real(real64) :: dz_grid, dr_grid, dr1, dr2, dz1, dz2, QQ1(2,2), QQ2(2,2)
  Integer(int32) :: i, ir, iz, iphi
  ierr = 0

  Do i=1,Npts

     Bout(i,1:3) = 0._real64
     If (r(i) .lt. bgrid_r(1) .OR. r(i) .gt. bgrid_r(bgrid_nr-1) &
          .OR. z(i) .lt. bgrid_z(1) .OR. z(i) .gt. bgrid_z(bgrid_nz-1)) Then
        ierr = 1
        Cycle
     Endif
  
     ir = Floor((r(i) - bgrid_r(1))/(bgrid_r(2)-bgrid_r(1))) + 1
     iz = Floor((z(i) - bgrid_z(1))/(bgrid_z(2)-bgrid_z(1))) + 1
     phi_tmp = phi(i)
     Do While (phi_tmp .lt. 0._real64)
        phi_tmp = phi_tmp + 2._real64*pi
     Enddo
     phi_tmp = Mod(phi_tmp,2._real64*pi/nsym)
     iphi = Floor((phi_tmp - bgrid_phi(1))/(bgrid_phi(2)-bgrid_phi(1))) + 1
     If (iphi .lt. 1 .OR. iphi .gt. bgrid_nphi - 1) Then
        Write(*,*) 'iphi out of range... should not happen',iphi
        Write(*,*) (iphi .lt. 1)
        Write(*,*) bgrid_nphi -1
        Write(*,*) "Error: giving up"
        Call fin_mpi(.true.)
     Endif
  
     dr_grid = bgrid_r(ir+1) - bgrid_r(ir)
     dz_grid = bgrid_z(iz+1) - bgrid_z(iz)
     dphi_grid = bgrid_phi(iphi+1) - bgrid_phi(iphi)

     phi_fac = (phi_tmp - bgrid_phi(iphi))/dphi_grid
  
     dr2 = bgrid_r(ir+1) - r(i)
     dr1 = dr_grid - dr2
     dz2 = bgrid_z(iz+1) - z(i)
     dz1 = dz_grid - dz2

     QQ1 = bgrid_br(ir:ir+1,iz:iz+1,iphi)
     QQ2 = bgrid_br(ir:ir+1,iz:iz+1,iphi+1)
     Bout(i,1) = ((1._real64-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1) + &
          phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1))/(dr_grid*dz_grid)
     QQ1 = bgrid_bz(ir:ir+1,iz:iz+1,iphi)
     QQ2 = bgrid_bz(ir:ir+1,iz:iz+1,iphi+1)
     Bout(i,2) = ((1._real64-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1) + &
          phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1))/(dr_grid*dz_grid)
     QQ1 = bgrid_bphi(ir:ir+1,iz:iz+1,iphi)
     QQ2 = bgrid_bphi(ir:ir+1,iz:iz+1,iphi+1)
     Bout(i,3) = ((1._real64-phi_fac)*(QQ1(1,1)*dr2*dz2 + QQ1(2,1)*dr1*dz2 + QQ1(1,2)*dr2*dz1 + QQ1(2,2)*dr1*dz1) + &
          phi_fac*(QQ2(1,1)*dr2*dz2 + QQ2(2,1)*dr1*dz2 + QQ2(1,2)*dr2*dz1 + QQ2(2,2)*dr1*dz1))/(dr_grid*dz_grid)
  Enddo

End Subroutine bfield_bgrid

End Module bgrid_module

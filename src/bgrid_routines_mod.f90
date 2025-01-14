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

Private

Real(real64), Public, Allocatable, Save :: bgrid_r(:), bgrid_z(:), bgrid_phi(:)
Real(real64), Public, Allocatable, Dimension(:,:,:), Save :: &
     bgrid_br,bgrid_bz, bgrid_bphi
Integer(int32), Public, Save :: bgrid_nr, bgrid_nz, bgrid_nphi, nsym

Public :: bfield_bgrid, open_bgrid_fields, close_bgrid_fields, read_ascii_version

Contains


Subroutine read_ascii_version(fname_base,verbose,R_min,R_max,Z_min,Z_max,nphi_minus_one)
  ! Author(s): J.D. Lore
  Use kind_mod, Only: real64, int32
  Use parallel_mod, Only : fin_mpi
  Implicit None
  Character(len=*), Intent(In) :: fname_base
  Logical, Intent(In) :: verbose
  Real(real64), Intent(Out) :: R_min, R_max, Z_min, Z_max
  Integer(int32), Intent(Out) :: nphi_minus_one
  Real(real64), Allocatable :: data(:,:,:)
  Character(len=1000) :: fname
  Integer(int32) :: iocheck, i, j, k
  !- End of header -------------------------------------------------------------

  fname = Trim(fname_base)//'_layout.dat'
  If (verbose) Write(*,'(a,a)') ' Reading ', Trim(fname)
  open(99,file=fname,IOSTAT=iocheck,status='old')
  If (iocheck /= 0) Then
     Write(*,*) 'Error opening file: ',Trim(fname)
     Call fin_mpi(.true.)
  Endif
  read(99, *) bgrid_nr, bgrid_nz, nphi_minus_one, nsym, R_min, R_max, Z_min, Z_max
  close(99)
  If (verbose) Then
     write(*,*) 'nr   = ',bgrid_nr
     write(*,*) 'nphi = ',nphi_minus_one
     write(*,*) 'nz   = ',bgrid_nz
     write(*,*) 'nsym = ',nsym
     write(*,*) 'R_(min,max) =',R_min, R_max
     write(*,*) 'Z_(min,max) =',Z_min, Z_max
  Endif

  ! Add one phi slice here for final arrays
  bgrid_nphi = nphi_minus_one + 1
  Allocate(bgrid_br(bgrid_nr,bgrid_nz,bgrid_nphi),   source=0._real64)
  Allocate(bgrid_bz(bgrid_nr,bgrid_nz,bgrid_nphi),   source=0._real64)
  Allocate(bgrid_bphi(bgrid_nr,bgrid_nz,bgrid_nphi), source=0._real64)
  
  ! Data array read from file is one slice smaller
  Allocate(data(bgrid_nr,bgrid_nz,nphi_minus_one),   source=0._real64)

  fname = Trim(fname_base)//'_r.dat'
  If (verbose) Write(*,'(a,a)') ' Reading ', Trim(fname)
  open(99,file=fname,IOSTAT=iocheck,status='old')
  If (iocheck /= 0) Then
     Write(*,*) 'Error opening file: ',Trim(fname)
     Call fin_mpi(.true.)
  Endif
  read(99, *) (((data(i,j,k), k=1,nphi_minus_one), j=1,bgrid_nz), i=1,bgrid_nr)
  close(99)
  bgrid_br(:,:,1:nphi_minus_one) = data
  
  fname = Trim(fname_base)//'_z.dat'
  If (verbose) Write(*,'(a,a)') ' Reading ', Trim(fname)
  open(99,file=fname,IOSTAT=iocheck,status='old')
  If (iocheck /= 0) Then
     Write(*,*) 'Error opening file: ',Trim(fname)
     Call fin_mpi(.true.)
  Endif
  read(99, *) (((data(i,j,k), k=1,nphi_minus_one), j=1,bgrid_nz), i=1,bgrid_nr)
  close(99)
  bgrid_bz(:,:,1:nphi_minus_one) = data

  fname = Trim(fname_base)//'_phi.dat'
  If (verbose)   Write(*,'(a,a)') ' Reading ', Trim(fname)
  open(99,file=fname,IOSTAT=iocheck,status='old')
  If (iocheck /= 0) Then
     Write(*,*) 'Error opening file: ',Trim(fname)
     Call fin_mpi(.true.)
  Endif
  read(99, *) (((data(i,j,k), k=1,nphi_minus_one), j=1,bgrid_nz), i=1,bgrid_nr)
  close(99)
  bgrid_bphi(:,:,1:nphi_minus_one) = data

  Deallocate(data)

  ! Write in last slice
  bgrid_br(:,:,bgrid_nphi)   = bgrid_br(:,:,1)
  bgrid_bz(:,:,bgrid_nphi)   = bgrid_bz(:,:,1)
  bgrid_bphi(:,:,bgrid_nphi) = bgrid_bphi(:,:,1)

  
End Subroutine read_ascii_version
!-----------------------------------------------------------------------------

Subroutine read_netcdf_version(fname,verbose,R_min, R_max, Z_min, Z_max, nphi_minus_one)
  ! Author(s): J.D. Lore
  Use kind_mod, Only: int32
  Use parallel_mod, Only : fin_mpi
  Use netcdf
  Implicit None
  Character(len=*), Intent(In) :: fname
  Logical, Intent(In) :: verbose
  Real(real64), Intent(Out) :: R_min, R_max, Z_min, Z_max
  Integer(int32), Intent(Out) :: nphi_minus_one

  ! Local variables
  Integer(int32) :: ncid, varid_br, varid_bz, varid_bphi
  Integer(int32) :: varid_nr, varid_nz, varid_nphi, varid_nsym
  Integer(int32) :: varid_rmin, varid_rmax, varid_zmin, varid_zmax
!  Integer(int32) :: dimid_nr, dimid_nz, dimid_nphi
  Integer(int32) :: ierr

  !- End of header -------------------------------------------------------------
  

  ! Open the NetCDF file for reading
  ierr = nf90_open(fname, nf90_nowrite, ncid)
  if (ierr /= nf90_noerr) then
     If (verbose) write(*,*) "Error opening NetCDF file:", trim(fname)
     If (verbose) write(*,*) nf90_strerror(ierr)
     Call fin_mpi(.true.)
  end if

  ! Inquire variable IDs
  ierr = nf90_inq_varid(ncid, "bgrid_br", varid_br)
  ierr = nf90_inq_varid(ncid, "bgrid_bz", varid_bz)
  ierr = nf90_inq_varid(ncid, "bgrid_bphi", varid_bphi)
  
  ! Inquire Scalars
  ierr = nf90_inq_varid(ncid, "nr", varid_nr)
  ierr = nf90_inq_varid(ncid, "nz", varid_nz)
  ierr = nf90_inq_varid(ncid, "nphi", varid_nphi)
  ierr = nf90_inq_varid(ncid, "nsym", varid_nsym)
  ierr = nf90_inq_varid(ncid, "R_min", varid_rmin)
  ierr = nf90_inq_varid(ncid, "R_max", varid_rmax)
  ierr = nf90_inq_varid(ncid, "Z_min", varid_zmin)
  ierr = nf90_inq_varid(ncid, "Z_max", varid_zmax)

  ! Read scalars
  ierr = nf90_get_var(ncid, varid_nr, bgrid_nr)
  ierr = nf90_get_var(ncid, varid_nz, bgrid_nz)
  ierr = nf90_get_var(ncid, varid_nphi, nphi_minus_one)
  ierr = nf90_get_var(ncid, varid_nsym, nsym)
  ierr = nf90_get_var(ncid, varid_rmin, R_min)
  ierr = nf90_get_var(ncid, varid_rmax, R_max)
  ierr = nf90_get_var(ncid, varid_zmin, Z_min)
  ierr = nf90_get_var(ncid, varid_zmax, Z_max)

  If (verbose) Then
     write(*,*) 'nr   = ',bgrid_nr
     write(*,*) 'nphi = ',nphi_minus_one
     write(*,*) 'nz   = ',bgrid_nz
     write(*,*) 'nsym = ',nsym
     write(*,*) 'R_(min,max) =',R_min, R_max
     write(*,*) 'Z_(min,max) =',Z_min, Z_max
  Endif
  

  ! Add one phi slice here for final arrays
  bgrid_nphi = nphi_minus_one + 1
  Allocate(bgrid_br(bgrid_nr,bgrid_nz,bgrid_nphi),   source=0._real64)
  Allocate(bgrid_bz(bgrid_nr,bgrid_nz,bgrid_nphi),   source=0._real64)
  Allocate(bgrid_bphi(bgrid_nr,bgrid_nz,bgrid_nphi), source=0._real64)

  ! Read arrays
  ierr = nf90_get_var(ncid, varid_br, bgrid_br(:,:,1:nphi_minus_one))
  ierr = nf90_get_var(ncid, varid_bz, bgrid_bz(:,:,1:nphi_minus_one))
  ierr = nf90_get_var(ncid, varid_bphi, bgrid_bphi(:,:,1:nphi_minus_one))

  ! Close the NetCDF file
  ierr = nf90_close(ncid)

  ! Write in last slice
  bgrid_br(:,:,bgrid_nphi)   = bgrid_br(:,:,1)
  bgrid_bz(:,:,bgrid_nphi)   = bgrid_bz(:,:,1)
  bgrid_bphi(:,:,bgrid_nphi) = bgrid_bphi(:,:,1)
  


  
End Subroutine read_netcdf_version



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
  Integer(int32) :: i, nphi_minus_one
  Real(real64) :: R_min,R_max,Z_min,Z_max
  Character(len=1000) :: fname_nc
  Logical :: nc_exists
  !- End of header -------------------------------------------------------------

  ! Read file, nc or ascii.
  ! This will allocate arrays and slices
  fname_nc = Trim(fname_base)//'.nc'
  Inquire(file=fname_nc, exist=nc_exists)
  If (nc_exists) Then
     If (verbose) Write(*,*) "Reading .nc version of bgrid: ", Trim(fname_nc)
     Call read_netcdf_version(fname_nc,verbose,R_min,R_max,Z_min,Z_max,nphi_minus_one)
  Else
     If (verbose) Write(*,*) "Reading ascii version of bgrid"
     Call read_ascii_version(fname_base,verbose,R_min,R_max,Z_min,Z_max,nphi_minus_one)
  End If


  ! Set up 1D arrays for interpolation
  Allocate(bgrid_r(bgrid_nr), source=0._real64)
  Allocate(bgrid_z(bgrid_nz), source=0._real64)
  Allocate(bgrid_phi(bgrid_nphi), source=0._real64)

  Do i=1,bgrid_nr
     bgrid_r(i) = R_min + 1.d0*(i-1)/(bgrid_nr-1) * (R_max - R_min)
  End Do
  Do i=1,bgrid_nz
     bgrid_z(i) = Z_min + 1.d0*(i-1)/(bgrid_nz-1) * (Z_max - Z_min)
  End Do
  Do i=1,bgrid_nphi
     bgrid_phi(i) = 2.d0*pi/nsym*(i-1)/(bgrid_nphi-1)
  End Do

End Subroutine read_bgrid_field_file

!-----------------------------------------------------------------------------
Subroutine open_bgrid_fields(fname,verbose)
  Implicit None
  Character(Len=*), Intent(In) :: fname
  Logical, Intent(In) :: verbose
  Call read_bgrid_field_file(fname,verbose)
End Subroutine open_bgrid_fields
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
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
  Integer(int32), Intent(In) :: Npts  
  Real(Real64), Intent(In), Dimension(Npts) :: r, z, phi
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

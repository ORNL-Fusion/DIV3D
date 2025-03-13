program write_to_netcdf
  Use kind_mod, Only : int32, real64
  Use bgrid_module
  use netcdf
  implicit none

  Character(len=1000) :: fname_base = 'Bgrid_thea_C640'
  Character(len=1000) :: fname

  
  ! Variables for NetCDF file and error handling
  integer(int32) :: ncid, varid_br, varid_bz, varid_bphi
  integer(int32) :: dimid_nr, dimid_nz, dimid_nphi, dimids(3)
  integer(int32) :: varid_nr, varid_nz, varid_nphi, varid_nsym, varid_rmin, varid_rmax, varid_zmin, varid_zmax
  integer(int32) :: ierr

  Real(real64) :: R_min,R_max,Z_min,Z_max
  Integer(int32) :: nphi_minus_one

  
  Call read_ascii_version(fname_base,.true.,R_min,R_max,Z_min,Z_max,nphi_minus_one)
  
  ! Create NetCDF file
  fname = Trim(fname_base)//'.nc'
  ierr = nf90_create(fname, nf90_clobber, ncid)
  if (ierr /= nf90_noerr) stop "Error creating NetCDF file"
  ! Define dimensions
  ierr = nf90_def_dim(ncid, "nr_dim", bgrid_nr, dimid_nr)
  ierr = nf90_def_dim(ncid, "nz_dim", bgrid_nz, dimid_nz)
  ierr = nf90_def_dim(ncid, "nphi_dim", nphi_minus_one, dimid_nphi)

  ! Define scalar variables
  ierr = nf90_def_var(ncid, "nr",    nf90_int,    varid=varid_nr)
  if (ierr /= nf90_noerr) stop "Error def_var nr"
  ierr = nf90_def_var(ncid, "nz",    nf90_int,    varid=varid_nz)
  ierr = nf90_def_var(ncid, "nphi",  nf90_int,    varid=varid_nphi)
  ierr = nf90_def_var(ncid, "nsym",  nf90_int,    varid=varid_nsym)
  ierr = nf90_def_var(ncid, "R_min", nf90_double, varid=varid_rmin)
  if (ierr /= nf90_noerr) stop "Error def_var R_min"  
  ierr = nf90_def_var(ncid, "R_max", nf90_double, varid=varid_rmax)
  ierr = nf90_def_var(ncid, "Z_min", nf90_double, varid=varid_zmin)
  ierr = nf90_def_var(ncid, "Z_max", nf90_double, varid=varid_zmax)
  
  
  ! Define grid variables
  dimids = (/dimid_nr, dimid_nz, dimid_nphi/)
  ierr = nf90_def_var(ncid, "bgrid_br",   nf90_double, dimids, varid_br)
  ierr = nf90_def_var(ncid, "bgrid_bz",   nf90_double, dimids, varid_bz)
  ierr = nf90_def_var(ncid, "bgrid_bphi", nf90_double, dimids, varid_bphi)
  
  ! End definition mode
  ierr = nf90_enddef(ncid)

  ! Write scalar values
  ierr = nf90_put_var(ncid, varid_nr, bgrid_nr)
  ierr = nf90_put_var(ncid, varid_nz, bgrid_nz)
  ierr = nf90_put_var(ncid, varid_nphi, nphi_minus_one)
  ierr = nf90_put_var(ncid, varid_nsym, nsym)
    write(*,*) 'jdl 1',r_min
  ierr = nf90_put_var(ncid, varid_rmin, R_min)
  if (ierr /= nf90_noerr) stop "Error put_var  R_min"  
  ierr = nf90_put_var(ncid, varid_rmax, R_max)
  ierr = nf90_put_var(ncid, varid_zmin, Z_min)
  ierr = nf90_put_var(ncid, varid_zmax, Z_max)

  write(*,*) 'jdl 2',r_min
  
  ! Write grid data
  ierr = nf90_put_var(ncid, varid_br, bgrid_br(:,:,1:nphi_minus_one))
  ierr = nf90_put_var(ncid, varid_bz, bgrid_bz(:,:,1:nphi_minus_one))
  ierr = nf90_put_var(ncid, varid_bphi, bgrid_bphi(:,:,1:nphi_minus_one))

  ! Close NetCDF file
  ierr = nf90_close(ncid)
  if (ierr /= nf90_noerr) stop "Error closing NetCDF file"

  write(*,*) "NetCDF file written successfully"
end program write_to_netcdf
  

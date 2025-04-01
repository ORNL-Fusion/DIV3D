Module output_routines
  Implicit None

Contains

  Subroutine init_hitline_netcdf(fname,nhitline)
    use netcdf
    implicit none

    character(len=*), intent(in) :: fname
    integer,          intent(in) :: nhitline

    integer :: ncid, ierr
    integer :: dimid_snapshot, dimid_hit_length
    integer :: varid_r, varid_z, varid_phi

    ! Master process creates the file (serial, no MPI communicator needed)
    ierr = nf90_create(trim(fname), IOR(NF90_CLOBBER, NF90_NETCDF4), ncid)
    if (ierr /= NF90_NOERR) stop "Error creating NetCDF file"

    ! Define dimensions
    ierr = nf90_def_dim(ncid, "snapshot", NF90_UNLIMITED, dimid_snapshot)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    ierr = nf90_def_dim(ncid, "hit_length", nhitline, dimid_hit_length)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    ! Define variables with dimensions [snapshot, hit_length]
    ierr = nf90_def_var(ncid, "r_hitline", NF90_DOUBLE, [dimid_snapshot, dimid_hit_length], varid_r)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    ierr = nf90_def_var(ncid, "z_hitline", NF90_DOUBLE, [dimid_snapshot, dimid_hit_length], varid_z)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    ierr = nf90_def_var(ncid, "phi_hitline", NF90_DOUBLE, [dimid_snapshot, dimid_hit_length], varid_phi)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    ierr = nf90_enddef(ncid)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    ierr = nf90_close(ncid)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

  end subroutine init_hitline_netcdf


  subroutine write_hitline_data_netcdf(fname, r_hitline, z_hitline, phi_hitline)
    use netcdf
    Use kind_mod, Only : real64
    implicit none

    character(len=*), intent(in) :: fname
    character(len=NF90_MAX_NAME) :: dim_name
    real(real64), intent(in) :: r_hitline(:)
    real(real64), intent(in) :: z_hitline(:)
    real(real64), intent(in) :: phi_hitline(:)

    integer :: ncid, ierr
    integer :: varid_r, varid_z, varid_phi
    integer :: dimid_snapshot
    integer(kind=4) :: snapshot_len
    integer :: start(2), count(2)

    ! Master process reopens file in write mode
    ierr = nf90_open(trim(fname), NF90_WRITE, ncid)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    ! Get dimension IDs and var IDs
    ierr = nf90_inq_dimid(ncid, "snapshot", dimid_snapshot)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    ! Current length of snapshot dimension
    ierr = nf90_inquire_dimension(ncid, dimid_snapshot, dim_name, snapshot_len)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    ierr = nf90_inq_varid(ncid, "r_hitline", varid_r)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    ierr = nf90_inq_varid(ncid, "z_hitline", varid_z)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    ierr = nf90_inq_varid(ncid, "phi_hitline", varid_phi)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    start = [snapshot_len+1, 1]
    count = [1, size(r_hitline)]

    ! write data
    ierr = nf90_put_var(ncid, varid_r, r_hitline, start=start, count=count)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    ierr = nf90_put_var(ncid, varid_z, z_hitline, start=start, count=count)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    ierr = nf90_put_var(ncid, varid_phi, phi_hitline, start=start, count=count)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

    ierr = nf90_close(ncid)
    If (ierr /= NF90_NOERR) Then
       Write(*,*) nf90_strerror(ierr)
       Stop
    End If

  End Subroutine write_hitline_data_netcdf

End Module output_routines

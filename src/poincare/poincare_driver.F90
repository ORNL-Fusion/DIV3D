!-----------------------------------------------------------------------------
!+ Driver for poincare in DIV3D
!-----------------------------------------------------------------------------
Program poincare_driver
  !
  ! Author(s): J.D. Lore 
  !
  Use kind_mod, Only: int32, real64
  Use poincare_namelist
  Use initialize_bfield_div3d, Only : init_bfield
  Use parallel_mod
  Use timing_mod, Only : get_elapsed_time, init_timing  
  Use math_routines_mod, Only : rlinspace
  Use fieldline_follow_mod, Only: follow_fieldlines_rzphi
  Use phys_const, Only: pi
  Use setup_bfield_module, Only : bfield

  Implicit none

  Integer(int32), Parameter :: print_progress_every = 1
  
  Integer :: tstart ! Timing variable
  Integer(int32) :: iocheck, i, j, offset, num_points_this_slice, ierr, &
       i_start, i_end, chunk_size, remainder, num_surfs_local, local_idx
  Logical :: verbose
  Real(real64) :: Rs, Zs, Ps
  Real(real64), Allocatable :: r1d(:), z1d(:), phistart_arr(:)
  Real(real64), Allocatable :: fl_r(:,:), fl_z(:,:), fl_p(:,:)
  Integer(int32), Allocatable :: ilg(:), fl_ierr(:)

  ! Temporary arrays for each fieldline call
  Real(real64), Allocatable :: rtemp(:), ztemp(:), ptemp(:)
  Integer(int32) :: ierrtemp, ilgtemp

  ! For output filename
  Character(len=256) :: outfilename

  !---------------------------------------------------------------------------
  ! 0. Initialization
  !---------------------------------------------------------------------------
  Call init_mpi
  Call init_timing

  verbose = .false.
  If (rank == 0) Then
     Call system_clock(tstart)
     verbose = .true.
     Write(*,'(/A)') '-------------------------------------------------------------------------'
     Write(*,'(a)') " Starting Poincare driver"
  End If

  ! Read namelists
  Call read_poincare_namelist(verbose)

  ! Determine workload distribution
  chunk_size = num_surfs / nprocs
  remainder  = mod(num_surfs, nprocs)

  ! Assign loop bounds for each rank
  i_start = rank * chunk_size + min(rank, remainder) + 1
  i_end   = i_start + chunk_size - 1
  If (rank < remainder) i_end = i_end + 1

  ! Number of surfaces this rank is responsible for
  num_surfs_local = i_end - i_start + 1

  Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !----------------------------------------------------------
  ! 1. Initialize magnetic field
  !----------------------------------------------------------
  Call init_bfield(verbose)

  !----------------------------------------------------------
  ! 2. Allocate variables and set start points
  !----------------------------------------------------------
  Allocate(r1d(num_surfs))
  Allocate(z1d(num_surfs))
  Allocate(phistart_arr(num_surfs))

  r1d(:)          = rlinspace(Rstart_poincare, Rend_poincare, num_surfs)
  z1d(:)          = rlinspace(Zstart_poincare, Zend_poincare, num_surfs)
  phistart_arr(:) = phistart_deg_poincare*pi/180.d0

  ! Local arrays to hold our subset of surfaces
  Allocate(fl_r(num_surfs_local, nsteps+1))
  Allocate(fl_z(num_surfs_local, nsteps+1))
  Allocate(fl_p(num_surfs_local, nsteps+1))
  Allocate(fl_ierr(num_surfs_local))
  Allocate(ilg(num_surfs_local))

  ! Allocate temp arrays for the fieldline routine
  Allocate(rtemp(nsteps+1), ztemp(nsteps+1), ptemp(nsteps+1))

  !----------------------------------------------------------
  ! 3. Follow fieldlines
  !----------------------------------------------------------
  Do i = i_start, i_end

     local_idx = i - i_start + 1

     Rs = r1d(i)
     Zs = z1d(i)
     Ps = phistart_arr(i)

     Call follow_fieldlines_rzphi(bfield, Rs, Zs, Ps, dphi_line_poincare, &
          nsteps, rtemp, ztemp, ptemp, ierrtemp, ilgtemp)

     ! Copy results to our local arrays
     fl_r(local_idx,:)   = rtemp
     fl_z(local_idx,:)   = ztemp
     fl_p(local_idx,:)   = ptemp
     fl_ierr(local_idx)  = ierrtemp
     ilg(local_idx)       = ilgtemp


     ! Optional progress message every so many surfaces
     If ( mod(i - i_start, print_progress_every) == 0 ) Then
        Write(*,*) 'Rank ', rank, ': processed surface ', local_idx, ' of ', num_surfs_local
     End If
  End Do

  Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  Write(*,*) 'Rank ', rank, ' has completed processing.'

  !----------------------------------------------------------
  ! 4. Write local data to unique file
  !----------------------------------------------------------
  ! surface_data.out.0000 for rank=0, etc.
  Write(outfilename, '(A,I4.4)') "surface_data.out.", rank
  Open(unit=99, file=outfilename, status="unknown", form="formatted", iostat=iocheck)
  If (iocheck /= 0) Then
     Write(*,*) 'Error opening output file on rank=', rank
     Stop 'Exiting: I/O Error in poincare_driver.F90'
  End If

  ! We can write exactly the same "sliced" data structure for this local subset
  Write(99,*) num_surfs_local, num_slices

  Do i = 0, num_slices-1
     offset = i*ind_poin/num_slices + 1
     num_points_this_slice = Int(nsteps + 1 - offset) / ind_poin + 1
     Write(99,*) dphi_slice*i*180.d0/pi, num_points_this_slice

     ! Now write fl_r, fl_z for each local surface
     Do j = 1, num_surfs_local
        Write(99,*) fl_r(j, offset:nsteps+1:ind_poin)
        Write(99,*) fl_z(j, offset:nsteps+1:ind_poin)
     End Do
  End Do

  Close(99)

  !----------------------------------------------------------
  ! 5. Print timing
  !----------------------------------------------------------
  If (rank == 0) Then
     Write(*,*) "Time spent in total: ", get_elapsed_time(tstart), " seconds"
     Write(*,*) "Time spent per line avg: ", get_elapsed_time(tstart)/num_surfs, " seconds"
  End If

  Call fin_mpi(.false.)

End Program poincare_driver

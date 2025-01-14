!-----------------------------------------------------------------------------
!+ Contains routines for MPI computing
!
! Use precompiler flag USE_MPIF08 to control MPI interface
!-----------------------------------------------------------------------------
Module parallel_mod
  
#ifdef USE_MPIF08
  Use mpi_f08      ! Modern MPI interface
#else
  Use mpi          ! Legacy MPI interface
#endif
  
  Implicit None

  Integer :: nprocs, ierr_mpi, rank, ret_code

#ifdef USE_MPIF08
  Type(MPI_Status) :: status
  Type(MPI_Request) :: request
#else
  Integer :: request
  Integer :: status(MPI_STATUS_SIZE)
#endif

Contains

  !-----------------------------------------------------------------------------
  !+ Initializes mpi computing
  !-----------------------------------------------------------------------------
  Subroutine init_mpi

    Implicit None

    ! initialize mpi
    call MPI_INIT(ierr_mpi)
    If (ierr_mpi .ne. MPI_SUCCESS) then
       Write(*,*) 'Error starting MPI program. Terminating.'
       Call MPI_ABORT(MPI_COMM_WORLD, ret_code, ierr_mpi)
    End If

    ! Determine rank of calling process. Rank = 0:nprocs-1
    Call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr_mpi)
    ! Determine # of processes
    Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr_mpi)

  End Subroutine init_mpi

  !-----------------------------------------------------------------------------
  !+ Finalizes mpi computing
  !-----------------------------------------------------------------------------
  Subroutine fin_mpi(iserror)
    Implicit None
    Logical, Intent(In) :: iserror

    If (iserror) Then
       !   Write(*,*) "Exiting due to error on rank: ", rank
       Call MPI_ABORT(MPI_COMM_WORLD,1,ierr_mpi)
    Else
       Call MPI_FINALIZE(ierr_mpi)
    End If

  End Subroutine fin_mpi

End Module parallel_mod

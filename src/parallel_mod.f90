!-----------------------------------------------------------------------------
!+ Contains routines for MPI computing
!-----------------------------------------------------------------------------
Module parallel_mod
  Use mpi
  Implicit None

  Integer :: nprocs, ierr_mpi, rank, ret_code, request
  Integer :: status(MPI_STATUS_SIZE)

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

!-----------------------------------------------------------------------------
!+ Contains routines for MPI computing
!-----------------------------------------------------------------------------
Module parallel_mod
!
! Description:
!
! History:
! Version   Date      Comment
! -------   ----      -------
!  1.0     12/14/2011  Begun JDL
! Author(s): J.D. Lore 12/14/2011 - xxx 
!
! Modules used:
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
if (ierr_mpi .ne. MPI_SUCCESS) then
  Write(6,*) 'Error starting MPI program. Terminating.'
  call MPI_ABORT(MPI_COMM_WORLD, ret_code, ierr_mpi)
end if

! Determine rank of calling process. Rank = 0:nprocs-1
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr_mpi)
! Determine # of processes
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr_mpi)

End Subroutine init_mpi

!-----------------------------------------------------------------------------
!+ Finalizses mpi computing
!-----------------------------------------------------------------------------
Subroutine fin_mpi(iserror)
Implicit None
Logical, Intent(In) :: iserror  

If (iserror) Then
   Write(*,*) "Exiting due to error on rank: ", rank
   Call MPI_ABORT(MPI_COMM_WORLD,1,ierr_mpi)
Else
   Call MPI_FINALIZE(ierr_mpi)
Endif

End Subroutine fin_mpi

End Module parallel_mod

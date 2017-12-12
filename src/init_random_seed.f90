Subroutine init_random_seed(myseed)
  Use kind_mod
  Implicit none
  Integer(int32) :: i, n, clock
  Integer(int32),Intent(in) :: myseed
  Integer(int32), Dimension(:), Allocatable :: seed
          
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  i=1          
  If (myseed .ne. 0) Then
!    Write(6,*) 'Using seed: ',myseed
    seed(:) = myseed
  Else
    Write(6,*) 'Using clock seed'
    Call SYSTEM_CLOCK(COUNT=clock)          
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  Endif

  Call RANDOM_SEED(PUT = seed)
          
  Deallocate(seed)
End Subroutine init_random_seed

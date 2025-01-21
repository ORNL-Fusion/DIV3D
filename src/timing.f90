Module timing_mod
  Use kind_mod, Only : real64, int32
  Implicit None

  ! Variables for timing
  Integer(int32), Private :: timing_count_rate = 1

Contains

  Subroutine init_timing
    Implicit None
    Call system_clock(count_rate = timing_count_rate)
  End Subroutine init_timing

  Function get_elapsed_time(tstart) Result(elapsed_time)
    Implicit None
    Integer, Intent(In) :: tstart
    Real(real64) :: elapsed_time
    Integer(int32) :: tend
    Call system_clock(tend)
    elapsed_time = Real(tend-tstart)/Real(timing_count_rate)
  End Function get_elapsed_time

End Module timing_mod

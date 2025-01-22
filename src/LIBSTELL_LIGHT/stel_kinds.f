      MODULE stel_kinds
!----------------------------------------------------------------------
!     Kind specifications
!----------------------------------------------------------------------
      Use Iso_fortran_env , Only: int16, int32, int64, real32, real64
      Implicit none
      Private
      Integer, Public, parameter :: rprec = real64
      Integer, Public, parameter :: iprec = int32

      END MODULE stel_kinds

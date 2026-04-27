Module test_assert_mod
  Use kind_mod, Only : int32, real64
  Implicit None

  Integer(int32) :: tests_failed = 0

Contains

  Subroutine assert_true(condition,message)
    Logical, Intent(in) :: condition
    Character(len=*), Intent(in) :: message

    If (.not. condition) Then
       tests_failed = tests_failed + 1
       Write(*,*) 'FAIL: ',Trim(message)
    End If
  End Subroutine assert_true


  Subroutine assert_equal_int(actual,expected,message)
    Integer(int32), Intent(in) :: actual, expected
    Character(len=*), Intent(in) :: message

    If (actual .ne. expected) Then
       tests_failed = tests_failed + 1
       Write(*,*) 'FAIL: ',Trim(message),' actual=',actual,' expected=',expected
    End If
  End Subroutine assert_equal_int


  Subroutine assert_near(actual,expected,tol,message)
    Real(real64), Intent(in) :: actual, expected, tol
    Character(len=*), Intent(in) :: message

    If (Abs(actual - expected) .gt. tol) Then
       tests_failed = tests_failed + 1
       Write(*,*) 'FAIL: ',Trim(message),' actual=',actual,' expected=',expected
    End If
  End Subroutine assert_near


  Subroutine finish_tests()
    If (tests_failed .ne. 0) Then
       Write(*,*) tests_failed,' unit test checks failed.'
       Stop 1
    End If

    Write(*,*) 'All unit test checks passed.'
  End Subroutine finish_tests

End Module test_assert_mod

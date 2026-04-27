Program test_inside_vessel
  Use kind_mod, Only : int32, real64
  Use inside_vessel_mod, Only : inside_vessel
  Use phys_const, Only : pi
  Use read_parts_mod, Only : is_AS_ves
  Use run_settings_namelist, Only : period, vessel_is_nearest_slice
  Use test_assert_mod
  Implicit None

  Integer(int32), Parameter :: ntor = 1_int32, npol = 4_int32
  Real(real64) :: Rves(ntor,npol), Zves(ntor,npol), Pves(ntor)

  period = 2._real64*pi
  vessel_is_nearest_slice = .true.
  is_AS_ves = .true.

  Pves = [0._real64]
  Rves(1,:) = [0._real64,1._real64,1._real64,0._real64]
  Zves(1,:) = [0._real64,0._real64,1._real64,1._real64]

  Call assert_true(inside_vessel(0.5_real64,0.5_real64,0._real64,Rves,Zves,Pves,ntor,npol), &
       'point inside square vessel')
  Call assert_true(.not. inside_vessel(1.5_real64,0.5_real64,0._real64,Rves,Zves,Pves,ntor,npol), &
       'point outside square vessel')
  Call assert_true(inside_vessel(0.5_real64,0.5_real64,3._real64*period,Rves,Zves,Pves,ntor,npol), &
       'phi wrapping preserves inside test')

  Call finish_tests

End Program test_inside_vessel

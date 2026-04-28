Program test_read_parts
  Use kind_mod, Only : int32, real64
  Use phys_const, Only : pi
  Use read_parts_mod
  Use run_settings_namelist, Only : fname_ptri, fname_ptri_mid, period
  Use test_assert_mod
  Implicit None

  Call test_load_2d_jpart
  Call test_make_triangles
  Call finish_tests

Contains

  Subroutine test_load_2d_jpart
    Character(len=300) :: label, fname
    Integer(int32) :: ntor, npol, msym
    Real(real64) :: Rpart(2,3), Zpart(2,3), Phipart(2,3)
    Logical :: force_non_AS_local

    fname = 'fixtures/tiny.2d.jpart'

    Call query_part(fname,ntor,npol,msym)
    Call assert_equal_int(ntor,2_int32,'query_part ntor')
    Call assert_equal_int(npol,3_int32,'query_part npol')
    Call assert_equal_int(msym,1_int32,'query_part nfp')

    Call load_2d_jpart(fname,label,ntor,npol,msym, &
         Rpart,Zpart,Phipart,force_non_AS_local)

    Call assert_near(Rpart(1,1),1._real64,1.e-12_real64,'load_2d_jpart converts R cm to m')
    Call assert_near(Zpart(2,1),0.1_real64,1.e-12_real64,'load_2d_jpart converts Z cm to m')
    Call assert_near(Phipart(2,1),0.5_real64*pi,1.e-12_real64,'load_2d_jpart converts Phi deg to rad')
    Call assert_true(.not. force_non_AS_local,'load_2d_jpart default force_non_AS')
  End Subroutine test_load_2d_jpart


  Subroutine test_make_triangles
    period = 2._real64*pi
    fname_ptri = 'unit_part_triangles.out'
    fname_ptri_mid = 'unit_part_triangle_mids.out'

    nparts = 1
    nt_max = 2
    np_max = 3

    Allocate(nt_parts(nparts),np_parts(nparts),part_type(nparts))
    Allocate(is_AS_part(nparts),force_non_AS(nparts))
    Allocate(Rparts(nparts,nt_max,np_max),Zparts(nparts,nt_max,np_max),Pparts(nparts,nt_max,np_max))
    Allocate(Pmins(nparts),Pmaxs(nparts))

    nt_parts = [2_int32]
    np_parts = [3_int32]
    part_type = [0_int32]
    is_AS_part = .false.
    force_non_AS = .true.

    Rparts(1,1,:) = [1._real64,1.1_real64,1.2_real64]
    Zparts(1,1,:) = [0._real64,0._real64,0._real64]
    Pparts(1,1,:) = [0._real64,0._real64,0._real64]

    Rparts(1,2,:) = [1._real64,1.1_real64,1.2_real64]
    Zparts(1,2,:) = [0.1_real64,0.1_real64,0.1_real64]
    Pparts(1,2,:) = [0.5_real64*pi,0.5_real64*pi,0.5_real64*pi]
    Pmins = [0._real64]
    Pmaxs = [0.5_real64*pi]

    Call make_triangles(.false.)
    Call assert_equal_int(ntri_parts(1),4_int32,'make_triangles creates expected triangle count')
    Call assert_near(xtri(1,1,1),1._real64,1.e-12_real64,'first triangle x coordinate')
    Call assert_near(ztri(1,1,1),0._real64,1.e-12_real64,'first triangle z coordinate')
  End Subroutine test_make_triangles

End Program test_read_parts

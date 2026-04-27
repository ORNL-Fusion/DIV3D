Program test_intersections
  Use kind_mod, Only : int32, real64
  Use diffusion, Only : check_line_for_intersections
  Use phys_const, Only : pi
  Use read_parts_mod
  Use run_settings_namelist, Only : period, vessel_int_is_last_point, vessel_is_nearest_slice
  Use test_assert_mod
  Implicit None

  Call test_part_hit
  Call test_vessel_hit
  Call test_no_hit_initializes_outputs
  Call finish_tests

Contains

  Subroutine reset_geometry
    If (Allocated(is_AS_part)) Deallocate(is_AS_part)
    If (Allocated(Pmins)) Deallocate(Pmins)
    If (Allocated(Pmaxs)) Deallocate(Pmaxs)
    If (Allocated(ntri_parts)) Deallocate(ntri_parts)
    If (Allocated(xtri)) Deallocate(xtri)
    If (Allocated(ytri)) Deallocate(ytri)
    If (Allocated(ztri)) Deallocate(ztri)
    If (Allocated(pmintri)) Deallocate(pmintri)
    If (Allocated(pmaxtri)) Deallocate(pmaxtri)
    If (Allocated(R_ves)) Deallocate(R_ves)
    If (Allocated(Z_ves)) Deallocate(Z_ves)
    If (Allocated(P_ves)) Deallocate(P_ves)
  End Subroutine reset_geometry


  Subroutine setup_vessel(rmin,rmax,zmin,zmax)
    Real(real64), Intent(in) :: rmin, rmax, zmin, zmax

    ntor_ves = 1
    npol_ves = 4
    Allocate(R_ves(ntor_ves,npol_ves),Z_ves(ntor_ves,npol_ves),P_ves(ntor_ves,npol_ves))
    R_ves(1,:) = [rmin,rmax,rmax,rmin]
    Z_ves(1,:) = [zmin,zmin,zmax,zmax]
    P_ves(1,:) = 0._real64
    is_AS_ves = .true.
    vessel_is_nearest_slice = .true.
    vessel_int_is_last_point = .true.
  End Subroutine setup_vessel


  Subroutine setup_no_parts
    nparts = 0
    Allocate(is_AS_part(0),Pmins(0),Pmaxs(0),ntri_parts(0))
    Allocate(xtri(0,0,3),ytri(0,0,3),ztri(0,0,3),pmintri(0,0),pmaxtri(0,0))
  End Subroutine setup_no_parts


  Subroutine test_part_hit
    Real(real64) :: rout(3), zout(3), phiout(3), pint(3), totL, theta
    Real(real64) :: r_hitline(2), z_hitline(2), phi_hitline(2)
    Integer(int32) :: iout(4)

    Call reset_geometry
    period = 2._real64*pi
    Call setup_vessel(0._real64,2._real64,-1._real64,1._real64)

    nparts = 1
    Allocate(is_AS_part(1),Pmins(1),Pmaxs(1),ntri_parts(1))
    Allocate(xtri(1,1,3),ytri(1,1,3),ztri(1,1,3),pmintri(1,1),pmaxtri(1,1))
    is_AS_part = .false.
    Pmins = [-1._real64]
    Pmaxs = [ 1._real64]
    ntri_parts = [1_int32]
    pmintri = -1._real64
    pmaxtri =  1._real64
    xtri(1,1,:) = [1._real64,1._real64,1._real64]
    ytri(1,1,:) = [-0.5_real64,0.5_real64,0._real64]
    ztri(1,1,:) = [-0.5_real64,-0.5_real64,0.5_real64]

    rout = [0.5_real64,1._real64,1.5_real64]
    zout = [0._real64,0._real64,0._real64]
    phiout = [0._real64,0._real64,0._real64]
    Call check_line_for_intersections(pint,iout,r_hitline,z_hitline,phi_hitline,2_int32, &
         1_int32,2_int32,rout,zout,phiout,2_int32,totL,.true.,.true.,theta)

    Call assert_equal_int(iout(1),1_int32,'check_line detects part hit')
    Call assert_equal_int(iout(2),1_int32,'part hit reports part index')
    Call assert_equal_int(iout(4),1_int32,'part hit reports segment index')
    Call assert_near(pint(1),1._real64,1.e-12_real64,'part hit x coordinate')
    Call assert_near(theta,1._real64,1.e-12_real64,'part hit theta')
  End Subroutine test_part_hit


  Subroutine test_vessel_hit
    Real(real64) :: rout(3), zout(3), phiout(3), pint(3), totL, theta
    Real(real64) :: r_hitline(0), z_hitline(0), phi_hitline(0)
    Integer(int32) :: iout(4)

    Call reset_geometry
    period = 2._real64*pi
    Call setup_no_parts
    Call setup_vessel(0._real64,1._real64,0._real64,1._real64)

    rout = [0.5_real64,0.75_real64,1.5_real64]
    zout = [0.5_real64,0.5_real64,0.5_real64]
    phiout = [0._real64,0._real64,0._real64]
    Call check_line_for_intersections(pint,iout,r_hitline,z_hitline,phi_hitline,0_int32, &
         1_int32,2_int32,rout,zout,phiout,2_int32,totL,.false.,.false.,theta)

    Call assert_equal_int(iout(1),2_int32,'check_line detects vessel hit')
    Call assert_near(totL,0._real64,1.e-12_real64,'calc_lc false keeps total length zero')
    Call assert_near(theta,0._real64,1.e-12_real64,'vessel hit theta initialized to zero')
  End Subroutine test_vessel_hit


  Subroutine test_no_hit_initializes_outputs
    Real(real64) :: rout(3), zout(3), phiout(3), pint(3), totL, theta
    Real(real64) :: r_hitline(0), z_hitline(0), phi_hitline(0)
    Integer(int32) :: iout(4)

    Call reset_geometry
    period = 2._real64*pi
    Call setup_no_parts
    Call setup_vessel(0._real64,2._real64,0._real64,2._real64)

    rout = [0.5_real64,0.75_real64,1._real64]
    zout = [0.5_real64,0.5_real64,0.5_real64]
    phiout = [0._real64,0._real64,0._real64]
    Call check_line_for_intersections(pint,iout,r_hitline,z_hitline,phi_hitline,0_int32, &
         1_int32,2_int32,rout,zout,phiout,2_int32,totL,.false.,.false.,theta)

    Call assert_equal_int(iout(1),0_int32,'no hit leaves ihit zero')
    Call assert_near(pint(1),0._real64,1.e-12_real64,'no hit initializes point')
    Call assert_near(theta,0._real64,1.e-12_real64,'no hit initializes theta')
  End Subroutine test_no_hit_initializes_outputs

End Program test_intersections

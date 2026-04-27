Program test_math_routines
  Use kind_mod, Only : int32, real64
  Use math_routines_mod, Only : int_line_curve, line_seg_facet_int, wrap_phi
  Use phys_const, Only : pi
  Use test_assert_mod
  Implicit None

  Call test_line_seg_facet_int
  Call test_int_line_curve
  Call test_wrap_phi
  Call finish_tests

Contains

  Subroutine test_line_seg_facet_int
    Real(real64) :: pa(3), pb(3), pc(3), p1(3), p2(3), p(3), sin_theta
    Integer(int32) :: ithit

    pa = [0._real64,0._real64,0._real64]
    pb = [1._real64,0._real64,0._real64]
    pc = [0._real64,1._real64,0._real64]

    p1 = [0.25_real64,0.25_real64,-1._real64]
    p2 = [0.25_real64,0.25_real64, 1._real64]
    Call line_seg_facet_int(pa,pb,pc,p1,p2,ithit,p,1.e-12_real64,.true.,sin_theta)
    Call assert_equal_int(ithit,1_int32,'segment hits triangle interior')
    Call assert_near(p(1),0.25_real64,1.e-12_real64,'hit x coordinate')
    Call assert_near(p(2),0.25_real64,1.e-12_real64,'hit y coordinate')
    Call assert_near(p(3),0._real64,1.e-12_real64,'hit z coordinate')
    Call assert_near(sin_theta,1._real64,1.e-12_real64,'theta for normal incidence')

    p1 = [1.25_real64,1.25_real64,-1._real64]
    p2 = [1.25_real64,1.25_real64, 1._real64]
    Call line_seg_facet_int(pa,pb,pc,p1,p2,ithit,p,1.e-12_real64,.true.,sin_theta)
    Call assert_equal_int(ithit,0_int32,'segment crosses plane outside triangle')

    p1 = [0.25_real64,0.25_real64,1._real64]
    p2 = [0.75_real64,0.25_real64,1._real64]
    Call line_seg_facet_int(pa,pb,pc,p1,p2,ithit,p,1.e-12_real64,.true.,sin_theta)
    Call assert_equal_int(ithit,0_int32,'parallel segment misses triangle plane')

    p1 = [0._real64,0._real64,-1._real64]
    p2 = [0._real64,0._real64, 1._real64]
    Call line_seg_facet_int(pa,pb,pc,p1,p2,ithit,p,1.e-12_real64,.false.,sin_theta)
    Call assert_equal_int(ithit,1_int32,'segment hits triangle vertex')
    Call assert_near(sin_theta,0._real64,1.e-12_real64,'theta disabled returns zero')
  End Subroutine test_line_seg_facet_int


  Subroutine test_int_line_curve
    Real(real64) :: rline(5), zline(5), p1(2), p2(2), pint(2), u
    Integer(int32) :: ierr, found_ind, int_count

    rline = [0._real64,1._real64,1._real64,0._real64,0._real64]
    zline = [0._real64,0._real64,1._real64,1._real64,0._real64]

    p1 = [-0.5_real64,0.5_real64]
    p2 = [ 0.5_real64,0.5_real64]
    Call int_line_curve(p1,p2,rline,zline,.true.,pint,ierr,u,found_ind,int_count)
    Call assert_equal_int(ierr,0_int32,'line intersects square')
    Call assert_equal_int(found_ind,4_int32,'first intersection segment index')
    Call assert_equal_int(int_count,1_int32,'first mode reports one intersection')
    Call assert_near(pint(1),0._real64,1.e-12_real64,'curve intersection R')
    Call assert_near(pint(2),0.5_real64,1.e-12_real64,'curve intersection Z')

    p1 = [-0.5_real64,1.5_real64]
    p2 = [ 0.5_real64,1.5_real64]
    Call int_line_curve(p1,p2,rline,zline,.true.,pint,ierr,u,found_ind,int_count)
    Call assert_equal_int(ierr,1_int32,'line misses square')
    Call assert_equal_int(int_count,0_int32,'miss reports zero intersections')

    p1 = [-0.5_real64,0.5_real64]
    p2 = [ 1.5_real64,0.5_real64]
    Call int_line_curve(p1,p2,rline,zline,.false.,pint,ierr,u,found_ind,int_count)
    Call assert_equal_int(ierr,0_int32,'line intersects square twice')
    Call assert_equal_int(int_count,2_int32,'all mode counts two intersections')
    Call assert_near(pint(1),0._real64,1.e-12_real64,'all mode returns last curve-order intersection R')
  End Subroutine test_int_line_curve


  Subroutine test_wrap_phi
    Real(real64) :: phi, period

    period = 2._real64*pi
    phi = -0.5_real64*pi
    Call wrap_phi(phi,period)
    Call assert_near(phi,1.5_real64*pi,1.e-12_real64,'negative phi wraps into period')

    phi = 2.5_real64*pi
    Call wrap_phi(phi,period)
    Call assert_near(phi,0.5_real64*pi,1.e-12_real64,'large phi wraps into period')
  End Subroutine test_wrap_phi

End Program test_math_routines

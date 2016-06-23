!-----------------------------------------------------------------------------
!+ Fieldline following routines
!-----------------------------------------------------------------------------
Module fieldline_following_mod
!
! Description:

! History:
! Version   Date      Comment
! -------   ----      -------
!  1.0    07/14/2011  
! 
! Author(s): J. Lore 7/2011 - xxx
!

Implicit None

Contains

!------------------------------------------------------------------------------------------
Subroutine follow_fieldline_rzphi_ci(rstart,zstart,phistart,dphi,nsteps,method, &
     rtol,atol,dphimin,nmax_step,r,z,phi,ifail,diffuse,dmag)
! General fieldline following routine.
Use kind_mod
Use bfield_xdr, Only: bint
Use integrator_routines_mod
Implicit None
! Input/output                      !See above for descriptions
Real(rknd),Intent(in) :: rstart,zstart,phistart
Real(rknd),Intent(in) :: dphi, dphimin, rtol, atol
Integer(iknd), Intent(in) :: method, nsteps, nmax_step
Real(rknd),Intent(out),Dimension(nsteps+1) :: &
  r,z,phi
Integer(iknd),Intent(out) :: ifail

! Local variables
Integer(iknd), Parameter :: n = 2
Real(rknd), Dimension(n) :: y
Real(rknd), Dimension(nsteps+1) :: xout
Real(rknd), Dimension(n,nsteps+1) :: yout
Real(rknd) :: x, dx

Integer(iknd) :: nok, nbad, i
Real(rknd) :: xtmp, x1, x2
Real(rknd), Dimension(n) :: ytmp

Real(rknd) :: xvec(3), bval_bint(3), perpdir1(3), perpdir2(3), alpha, &
bphi, br, bz, dca, dL, delta_x, dmag, dsa, rnum
Integer(iknd) :: idiv

Logical, intent(in) :: diffuse
Real(rknd), Parameter :: pi = 3.141592653589793238462643383279502_rknd
!External fl_derivs_dphi_bint

ifail = 1
r = 0._rknd
z = 0._rknd
phi =  0._rknd

If (method .eq. 1) Then  ! Fixed step size RK4 integration
!  y(1) = Rstart
!  y(2) = Zstart
!  x = phistart
!  dx = dphi
!
!  Call rk45_fixed_step_integrate(y,n,x,dx,nsteps,fl_derivs_dphi_bint,yout,xout)
!
!  r = yout(1,:)
!  z = yout(2,:)
!  phi = xout

  write(*,*) 'not implemented'
  stop
Elseif (method .eq. 2) Then ! Adaptive step size RK45 integration, fixed output step sizes

  ! Save initial point
  r(1) = Rstart ;   z(1) = Zstart ;   phi(1) = Phistart

  ! first step
  
  ytmp(1) = Rstart
  ytmp(2) = Zstart
  xtmp = Phistart

  Do i = 1,nsteps

!    write(*,*) i,nsteps
    
    y = ytmp
    x1 = xtmp
    x2 = xtmp + dphi

    Call rk45_adp_step_integrate(y,n,x1,x2,rtol,atol,dphi,dphimin,nok,nbad,fl_derivs_dphi_bint,ytmp,xtmp,nmax_step)
    
    r(i+1) = ytmp(1)
    z(i+1) = ytmp(2)
    phi(i+1) = xtmp

    if (diffuse) then 

    dL = sqrt(r(i+1)*r(i+1) + r(i)*r(i) - 2._rknd*r(i+1)*r(i)*cos(phi(i+1)-phi(i)) & 
         + z(i+1)*z(i+1) + z(i)*z(i) - 2._rknd*z(i+1)*z(i))

    xvec(1) = r(i+1)
    xvec(2) = phi(i+1)
    xvec(3) = z(i+1)
    call bint(xvec,bval_bint,idiv)  
    if (idiv .ne. 0) then
      write(*,*) 'adfadsfdasf'
      stop
      exit
    endif
 
    br   = bval_bint(1)
    bphi = bval_bint(2)
    bz   = bval_bint(3)


    ! B cross z^hat
    perpdir1(1) = bphi   !r 
    perpdir1(2) = -br    !phi
    perpdir1(3) = 0.d0   !z
    perpdir1 = perpdir1/dsqrt(bphi*bphi + br*br)

    ! B cross r^hat
    perpdir2(1) = 0.d0
    perpdir2(2) = bz
    perpdir2(3) = -bphi
    perpdir2 = perpdir2/dsqrt(bz*bz + bphi*bphi)

!    rnum = rand()
    Call Random_number(rnum)

    alpha = 2.d0*pi*(-1.d0 + 2.d0*rnum) ! -2pi to 2pi kick    
    dca = dcos(alpha)
    dsa = dsin(alpha)

    delta_x = dsqrt(dmag*dL)
    r(i+1)   = r(i+1)   + delta_x*(dca*perpdir1(1) + dsa*perpdir2(1))
    phi(i+1) = phi(i+1) + delta_x*(dca*perpdir1(2) + dsa*perpdir2(2))
    z(i+1)   = z(i+1)   + delta_x*(dca*perpdir1(3) + dsa*perpdir2(3))
  endif

  Enddo

Else
  Write(6,*) 'Unknown method in follow_fieldline_rzphi_ci'
  Stop
Endif

ifail = 0
Return

End Subroutine follow_fieldline_rzphi_ci
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


!------------------------------------------------------------------------------------------
Subroutine fl_derivs_dphi_bint(n,phi,RZ,df)
Use kind_mod
Use bfield_xdr, Only: bint
!Use vmec_coils_mod
Implicit None
Real(rknd), Intent(In) :: phi
Integer(iknd), Intent(In) :: n
Real(rknd), Intent(In), Dimension(n) :: RZ
Real(rknd), Intent(Out), Dimension(n) :: df

Real(rknd), Dimension(3) :: xvec, bval
Integer(iknd) :: idiv
!Real(rknd) :: Bx, By, Bz

idiv = 0
bval = 0._rknd

!if (.true.) Then
  xvec(1) = RZ(1)
  xvec(2) = phi
  xvec(3) = RZ(2)
  Call bint(xvec,bval,idiv)
!Else
!  Call bfield_bs_jdl(RZ(1)*cos(phi),RZ(1)*sin(phi),RZ(2),&
!       coil,current,nfil,Bx,By,Bz)
!  bval(1) = Bx*cos(phi) + By*sin(phi)
!  bval(2) = -Bx*sin(phi) + By*cos(phi)
!  bval(3) = Bz
!Endif

If (idiv .ne. 0) Then
  Write(6,*) 'Error in fl_derivs_dphi_bint'
  Stop
Endif
df(1) = RZ(1)*bval(1)/bval(2)
df(2) = RZ(1)*bval(3)/bval(2)
Return
End Subroutine fl_derivs_dphi_bint
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

End Module fieldline_following_mod

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

!-----------------------------------------------------------------------------
!+ Follows fieldlines in cylindrical coords and checks for intersection at each step
!-----------------------------------------------------------------------------
Subroutine follow_fieldline_rzphi_and_check(rstart,zstart,phistart,dphi,nsteps,r,z,phi,diffuse,drat,ifail,&
period,imin,lsfi_tol,linenum,pint,iout)

Use kind_mod                ! Import rknd, iknd specifications
use bfield_xdr, Only: &
! Imported subroutines
bint
Use read_parts_mod
Implicit None

! Input/output                      !See above for descriptions
Integer(iknd),Intent(in) :: nsteps,linenum
Real(rknd),Intent(in) :: rstart,zstart,phistart
Real(rknd),Intent(in) :: dphi
Real(rknd),Intent(in) :: drat, period, lsfi_tol
logical,Intent(in) :: diffuse
Real(rknd),Intent(out),Dimension(nsteps+1) :: &
  r,z,phi
Integer(iknd),Intent(out) :: ifail

! Local scalars
Real(rknd) :: br,bz,bphi,rnum, dL_tmp
Real(Rknd) :: k1r,k2r,k3r,k4r,k1z,k2z,k3z,k4z
Integer(iknd) :: ii,idiv, imin, iclose,ipart_close, imin_tmp, iout(4),jj
real(rknd) :: alpha, delta_x, dL, dca, dsa, adp

! Local arrays (1D)
Real(rknd) :: bval_bint(3), dmin, pint(3)
Real(rknd) :: xvec(3), perpdir1(3),perpdir2(3)

! Local Parameters
Real(rknd), parameter :: pi = 3.1415926535897932384626433832795_rknd

!- End of header -------------------------------------------------------------

! Initialize arrays
r(:) = 0._rknd
z(:) = 0._rknd
phi(:) = 0._rknd
r(1)   = rstart
z(1)   = zstart
phi(1) = phistart

adp = abs(dphi)
iclose = 0
ifail = 0
imin = 1
Do ii=2,nsteps+1
    phi(ii) = phi(ii-1) + dphi
   
    xvec(1) = r(ii-1)
    xvec(2) = phi(ii-1)
    xvec(3) = z(ii-1)
    call bint(xvec,bval_bint,idiv)  
    if (idiv .ne. 0) then
!    Do jj = ii+1,nsteps+1
!      Call Random_number(rnum)
!    Enddo
!      write(6,*),'ff ifail1',xvec,ii
      exit
    endif

    k1r = dphi*r(ii-1)*bval_bint(1)/bval_bint(2)
    k1z = dphi*r(ii-1)*bval_bint(3)/bval_bint(2)

    xvec(1) = r(ii-1) + 0.5_rknd*k1r
    xvec(2) = phi(ii-1) + 0.5_rknd*dphi
    xvec(3) = z(ii-1) + 0.5_rknd*k1z
    call bint(xvec,bval_bint,idiv)  
    if (idiv .ne. 0) then
!    Do jj = ii+1,nsteps+1
!      Call Random_number(rnum)
!    Enddo
!      write(6,*),'ff ifail2',xvec,ii
      exit
    endif

    k2r = dphi*(r(ii-1)+0.5_rknd*k1r)*bval_bint(1)/bval_bint(2)
    k2z = dphi*(r(ii-1)+0.5_rknd*k1r)*bval_bint(3)/bval_bint(2)

    xvec(1) = r(ii-1) + 0.5_rknd*k2r
    xvec(2) = phi(ii-1) + 0.5_rknd*dphi
    xvec(3) = z(ii-1) + 0.5_rknd*k2z
    call bint(xvec,bval_bint,idiv)  
    if (idiv .ne. 0) then
!    Do jj = ii+1,nsteps+1
!      Call Random_number(rnum)
!    Enddo
!      write(6,*),'ff ifail3',xvec,ii
      exit
    endif

    k3r = dphi*(r(ii-1)+0.5_rknd*k2r)*bval_bint(1)/bval_bint(2)
    k3z = dphi*(r(ii-1)+0.5_rknd*k2r)*bval_bint(3)/bval_bint(2) 

    xvec(1) = r(ii-1) + k3r
    xvec(2) = phi(ii-1) + dphi
    xvec(3) = z(ii-1) + k3z
    call bint(xvec,bval_bint,idiv)  
 
    if (idiv .ne. 0) then
!      write(6,*),'ff ifail4',xvec,ii
!    Do jj = ii+1,nsteps+1
!      Call Random_number(rnum)
!    Enddo
      exit
    endif

    k4r = dphi*(r(ii-1)+k3r)*bval_bint(1)/bval_bint(2)
    k4z = dphi*(r(ii-1)+k3r)*bval_bint(3)/bval_bint(2)

    r(ii) = r(ii-1) + k1r/6._rknd + k2r/3._rknd + k3r/3._rknd + k4r/6._rknd
    z(ii) = z(ii-1) + k1z/6._rknd + k2z/3._rknd + k3z/3._rknd + k4z/6._rknd




  if (diffuse) then 

!    dL = r(ii)*Dabs(phi(ii)-phi(ii-1))
    dL = sqrt(r(ii)*r(ii) + r(ii-1)*r(ii-1) - 2._rknd*r(ii)*r(ii-1)*cos(phi(ii)-phi(ii-1)) & 
         + z(ii)*z(ii) + z(ii-1)*z(ii-1) - 2._rknd*z(ii)*z(ii-1))

    xvec(1) = r(ii)
    xvec(2) = phi(ii)
    xvec(3) = z(ii)
    call bint(xvec,bval_bint,idiv)  

    if (idiv .ne. 0) then
!      write(6,*),'ff ifail5',xvec,ii
!    Do jj = ii+1,nsteps+1
!      Call Random_number(rnum)
!    Enddo
      exit
    endif
 
    br   = bval_bint(1)
    bz   = bval_bint(3)
    bphi = bval_bint(2)

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

    Call Random_number(rnum)

    alpha = 2.d0*pi*(-1.d0 + 2.d0*rnum) ! -2pi to 2pi kick    
    dca = dcos(alpha)
    dsa = dsin(alpha)

    delta_x = dsqrt(drat*dL)
    r(ii)   = r(ii)   + delta_x*(dca*perpdir1(1) + dsa*perpdir2(1))
    phi(ii) = phi(ii) + delta_x*(dca*perpdir1(2) + dsa*perpdir2(2))
    z(ii)   = z(ii)   + delta_x*(dca*perpdir1(3) + dsa*perpdir2(3))
  endif



  

  ! Check for nearby triangles
  Call find_nearby_triangles_v2(r(ii),phi(ii),z(ii),dL,period,dmin,imin_tmp,iclose,r(ii-1),phi(ii-1),z(ii-1),&
       lsfi_tol,ipart_close,linenum)

  ! Check each triangle identified for an intersection
  Call check_line_segment_for_intersections(period,pint,iout, &
       linenum,lsfi_tol,1_iknd,r(ii-1:ii),z(ii-1:ii),phi(ii-1:ii),0_iknd)
!       r_hitline,z_hitline,phi_hitline,nhitline, &


  deallocate(near_part,near_tri)
  ! If intersection found, quit
 
  If (iout(1) .ne. 0) Then
    iout(4) = ii - 1
    ! SET UP HITLINE HERE!! QQ
    Exit
  Endif
Enddo ! nsteps

If (idiv .ne. 0 ) ifail = ii - 1

EndSubroutine follow_fieldline_rzphi_and_check



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
Real(rknd) :: Bx, By, Bz

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

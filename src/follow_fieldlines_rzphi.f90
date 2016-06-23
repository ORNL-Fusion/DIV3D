!-----------------------------------------------------------------------------
!+ Follows fieldlines in cylindrical coords
!-----------------------------------------------------------------------------
Subroutine follow_fieldline_rzphi(rstart,zstart,phistart,dphi,nsteps,r,z,phi,diffuse,drat,ifail,method)
!b_routine)
!
! Description: 
!
! Inputs:
!   rstart
!   zstart
!   phistart
!   dphi
!   nsteps
!   diffuse
!   drat
!
! Outputs:
!   r
!   z
!   phi
!   ifail
!
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/20/2011  Ported from g3d.  JDL
!  1.1     07/15/2011  Modified for W7X routines. JDL
! Author(s): J.D. Lore - 04/20/2011 - xxx
!
! Modules used:
Use kind_mod                ! Import rknd, iknd specifications
use bfield_xdr, Only: &
! Imported subroutines
bint
Use g3d_module, Only : bfield_geq_bicub
Implicit None

! Input/output                      !See above for descriptions
Integer(iknd),Intent(in) :: nsteps
Real(rknd),Intent(in) :: rstart,zstart,phistart
Real(rknd),Intent(in) :: dphi
Real(rknd),Intent(in) :: drat
logical,Intent(in) :: diffuse
Real(rknd),Intent(out),Dimension(nsteps+1) :: &
  r,z,phi
Integer(iknd),Intent(out) :: ifail
Integer(iknd),Intent(In) :: method

! Local scalars
Real(rknd) :: br,bz,bphi,rnum
Real(Rknd) :: k1r,k2r,k3r,k4r,k1z,k2z,k3z,k4z
Integer(iknd) :: ii,idiv
real(rknd) :: alpha, delta_x, dL, dca, dsa

! Local arrays (1D)
Real(rknd) :: bval_bint(3), Bval(1,3)
Real(rknd) :: xvec(3), perpdir1(3),perpdir2(3)

! Local Parameters
Real(rknd), parameter :: pi = 3.1415926535897932384626433832795_rknd

!interface
!subroutine bint(xvec,bval,idiv)
!real(kind(1.d+0)), intent(in) :: xvec(3)
!real(kind(1.d+0)), intent(out) :: bval(3)
!integer, intent(out) :: idiv
!end subroutine bint
!end interface

!- End of header -------------------------------------------------------------

! Initialize arrays
r(:) = 0._rknd
z(:) = 0._rknd
phi(:) = 0._rknd
r(1)   = rstart
z(1)   = zstart
phi(1) = phistart

ifail = 0
Do ii=2,nsteps+1
    phi(ii) = phi(ii-1) + dphi
   
    xvec(1) = r(ii-1)
    xvec(2) = phi(ii-1)
    xvec(3) = z(ii-1)

    If (method .eq. 0) Then
      call bint(xvec,bval_bint,idiv)  
    Else
      Call bfield_geq_bicub((/xvec(1)/),(/xvec(3)/),1,Bval,idiv)     
      bval_bint(1) = Bval(1,1)
      bval_bint(3) = Bval(1,2)
      bval_bint(2) = Bval(1,3)
    Endif

    if (idiv .ne. 0) then
!      write(6,*),'ff ifail1',xvec,ii
      exit
    endif

    k1r = dphi*r(ii-1)*bval_bint(1)/bval_bint(2)
    k1z = dphi*r(ii-1)*bval_bint(3)/bval_bint(2)

    xvec(1) = r(ii-1) + 0.5_rknd*k1r
    xvec(2) = phi(ii-1) + 0.5_rknd*dphi
    xvec(3) = z(ii-1) + 0.5_rknd*k1z

    If (method .eq. 0) Then
      call bint(xvec,bval_bint,idiv)  
    Else
      Call bfield_geq_bicub((/xvec(1)/),(/xvec(3)/),1,Bval,idiv)     
      bval_bint(1) = Bval(1,1)
      bval_bint(3) = Bval(1,2)
      bval_bint(2) = Bval(1,3)
    Endif

    if (idiv .ne. 0) then
!      write(6,*),'ff ifail2',xvec,ii
      exit
    endif

    k2r = dphi*(r(ii-1)+0.5_rknd*k1r)*bval_bint(1)/bval_bint(2)
    k2z = dphi*(r(ii-1)+0.5_rknd*k1r)*bval_bint(3)/bval_bint(2)

    xvec(1) = r(ii-1) + 0.5_rknd*k2r
    xvec(2) = phi(ii-1) + 0.5_rknd*dphi
    xvec(3) = z(ii-1) + 0.5_rknd*k2z

    If (method .eq. 0) Then
      call bint(xvec,bval_bint,idiv)  
    Else
      Call bfield_geq_bicub((/xvec(1)/),(/xvec(3)/),1,Bval,idiv)     
      bval_bint(1) = Bval(1,1)
      bval_bint(3) = Bval(1,2)
      bval_bint(2) = Bval(1,3)
    Endif

    if (idiv .ne. 0) then
!      write(6,*),'ff ifail3',xvec,ii
      exit
    endif

    k3r = dphi*(r(ii-1)+0.5_rknd*k2r)*bval_bint(1)/bval_bint(2)
    k3z = dphi*(r(ii-1)+0.5_rknd*k2r)*bval_bint(3)/bval_bint(2) 

    xvec(1) = r(ii-1) + k3r
    xvec(2) = phi(ii-1) + dphi
    xvec(3) = z(ii-1) + k3z

    If (method .eq. 0) Then
      call bint(xvec,bval_bint,idiv)  
    Else
      Call bfield_geq_bicub((/xvec(1)/),(/xvec(3)/),1,Bval,idiv)     
      bval_bint(1) = Bval(1,1)
      bval_bint(3) = Bval(1,2)
      bval_bint(2) = Bval(1,3)
    Endif
 
    if (idiv .ne. 0) then
!      write(6,*),'ff ifail4',xvec,ii
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

    If (method .eq. 0) Then
      call bint(xvec,bval_bint,idiv)  
    Else
      Call bfield_geq_bicub((/xvec(1)/),(/xvec(3)/),1,Bval,idiv)     
      bval_bint(1) = Bval(1,1)
      bval_bint(3) = Bval(1,2)
      bval_bint(2) = Bval(1,3)
    Endif

    if (idiv .ne. 0) then
!      write(6,*),'ff ifail5',xvec,ii
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

!    rnum = rand()
    Call Random_number(rnum)

    alpha = 2.d0*pi*(-1.d0 + 2.d0*rnum) ! -2pi to 2pi kick    
    dca = dcos(alpha)
    dsa = dsin(alpha)

    delta_x = dsqrt(drat*dL)
    
!    delta_x = dsqrt(drat*dL*(1.+0.4d0*sin(3.*phi(ii)))*1.55d0/sqrt(br*br+bz*bz+bphi*bphi)**2)
!    delta_x = dsqrt(drat*dL*(1.+0.99*sin(3.*phi(ii)))*0.3d0**2/sqrt(br*br+bz*bz+bphi*bphi)**2)
    
!    delta_x = 0.d0
!    If ((r(ii) .gt. 1.d0) .AND. (abs(z(ii)) .lt. 0.1d0)) Then
!      delta_x = sqrt(drat*dL*(1.+0.5*sin(3.*phi(ii))))
!    Endif
    


    r(ii)   = r(ii)   + delta_x*(dca*perpdir1(1) + dsa*perpdir2(1))
    phi(ii) = phi(ii) + delta_x*(dca*perpdir1(2) + dsa*perpdir2(2))
    z(ii)   = z(ii)   + delta_x*(dca*perpdir1(3) + dsa*perpdir2(3))
  endif
Enddo


If (idiv .ne. 0 ) ifail = ii - 1

EndSubroutine follow_fieldline_rzphi



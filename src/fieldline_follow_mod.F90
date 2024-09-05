!-----------------------------------------------------------------------------
!
!   Routines/modules related to fieldline following
!
!   Contains:
!     Subroutine follow_fieldlines_rzphi
!     Subroutine fl_derivs_fun
!     Subroutine rk45_fixed_step_integrate
!     Subroutine rk4_core
!   
!-----------------------------------------------------------------------------

Module fieldline_follow_mod
  Use kind_mod, Only: int32
  Use bfield, Only : bfield_type
  Implicit None
  Private
  Save

  Public :: follow_fieldlines_rzphi
  Public :: follow_fieldlines_rzphi_diffuse
  Public :: follow_fieldlines_rzphi_dz
  Private :: fl_derivs_fun, fl_derivs_fun_dz

  Interface follow_fieldlines_rzphi_dz
    Module Procedure follow_fieldlines_rzphi_dz_Npts
    Module Procedure follow_fieldlines_rzphi_dz_1pt
  End Interface follow_fieldlines_rzphi_dz

  Interface follow_fieldlines_rzphi
    Module Procedure follow_fieldlines_rzphi_Npts
    Module Procedure follow_fieldlines_rzphi_1pt
  End Interface follow_fieldlines_rzphi
 
Contains


!-----------------------------------------------------------------------------
!+ Follows fieldlines in cylindrical coords with diffusion
!-----------------------------------------------------------------------------
Subroutine follow_fieldlines_rzphi_diffuse(bfield,rstart,zstart,phistart,Npts,dphi,nsteps,r,z,phi,ierr,i_last_good,dmag)
!
! Description: 
!  Follows fieldlines by integrating along toroidal angle in cylindrical coordinates. 
!
!  Integration method is RK45 with fixed step size
!
! Input:
!  rstart, zstart, phistart : real64(Npts) : Launching points for fieldlines (m,m,radians)
!  Npts : int32 : Number of starting points
!  dphi : real64 : Integration step size (radians)
!  nsteps : int32 : Number of integration steps to take
!  dmag : magnetic diffusivity (m^2/m)
! Output:
!  r,z,phi : real64(Npts,nsteps+1) : Field line trajectories (m,m,radians)
!  ierr : int32(Npts) : Error flag for each fl. 0 indicates no error, 1 indicates an error from bfield_geq_bicub
!  i_last_good : int32(Npts) : If an error occured for the fl, this array gives the last 'good' index, i.e., before the 
!                             error occured
!
! Calls:
!  Subroutine rk45_fixed_step_integrate
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/20/2011  Original code
!
! Author(s): J.D Lore - 04/20/2011
!
! Modules used:
  Use kind_mod, Only: real64, int32
  Implicit None
  ! Input/output                      !See above for descriptions
  Type(bfield_type), Intent(In) :: bfield
  Integer(int32),Intent(in) :: Npts, nsteps
  Real(real64),Intent(in),Dimension(Npts) :: rstart(Npts),zstart(Npts),phistart(Npts)
  Real(real64),Intent(in) :: dphi, dmag
  
Real(real64),Intent(out),Dimension(Npts,nsteps+1) :: r,z,phi
Integer(int32),Intent(out) :: ierr(Npts), i_last_good(Npts)
! Local variables
Integer(int32), Parameter :: n = 2
Real(real64) :: y(n),x,dx,xout(nsteps+1),yout(n,nsteps+1)
Integer(int32) :: ierr_rk45, i_last_good_rk45, ipt

!- End of header -------------------------------------------------------------

r = 0._real64
z = 0._real64
phi = 0._real64
dx = dphi
Do ipt = 1,Npts 
  y(1) = rstart(ipt)
  y(2) = zstart(ipt)
  x = phistart(ipt)
  Call rk45_fixed_step_integrate_diffuse(bfield,y,n,x,dx,nsteps,fl_derivs_fun,yout,xout,ierr_rk45,i_last_good_rk45,dmag)
  r(ipt,1:nsteps+1)   = yout(1,1:nsteps+1)
  z(ipt,1:nsteps+1)   = yout(2,1:nsteps+1)
  phi(ipt,1:nsteps+1) = xout
  ierr(ipt) = ierr_rk45
  i_last_good(ipt) = i_last_good_rk45
Enddo

EndSubroutine follow_fieldlines_rzphi_diffuse

!-----------------------------------------------------------------------------
!+ Follows fieldlines in cylindrical coords (from one point)
!-----------------------------------------------------------------------------
Subroutine follow_fieldlines_rzphi_1pt(bfield,rstart,zstart,phistart,dphi,nsteps,r,z,phi,ierr,i_last_good)
  Use kind_mod, Only: real64, int32
Implicit None
! Input/output
  Type(bfield_type), Intent(In) :: bfield
Integer(int32),Intent(in) :: nsteps
Real(real64),Intent(in) :: rstart,zstart,phistart
Real(real64),Intent(in) :: dphi
Real(real64),Intent(out),Dimension(nsteps+1) :: r,z,phi
Integer(int32),Intent(out) :: ierr, i_last_good

! Local
Integer(int32), Parameter :: n = 2
Real(real64) :: y(n),x,dx,xout(nsteps+1),yout(n,nsteps+1)
Integer(int32) :: ierr_rk45, i_last_good_rk45

r = 0._real64
z = 0._real64
phi = 0._real64
dx = dphi
y(1) = rstart
y(2) = zstart
x = phistart
Call rk45_fixed_step_integrate(bfield,y,n,x,dx,nsteps,fl_derivs_fun,yout,xout,ierr_rk45,i_last_good_rk45)
r(1:nsteps+1)   = yout(1,1:nsteps+1)
z(1:nsteps+1)   = yout(2,1:nsteps+1)
phi(1:nsteps+1) = xout
ierr = ierr_rk45
i_last_good = i_last_good_rk45

End Subroutine follow_fieldlines_rzphi_1pt

!-----------------------------------------------------------------------------
!+ Follows fieldlines in cylindrical coords (from one point)
!-----------------------------------------------------------------------------
Subroutine follow_fieldlines_rzphi_dz_1pt(bfield,rstart,zstart,phistart,dz,nsteps,r,z,phi,ierr,i_last_good)
  Use kind_mod, Only: real64, int32
  Implicit None
  ! Input/output 
  Type(bfield_type), Intent(In) :: bfield
Integer(int32),Intent(in) :: nsteps
Real(real64),Intent(in) :: rstart,zstart,phistart
Real(real64),Intent(in) :: dz

Real(real64),Intent(out),Dimension(nsteps+1) :: r,z,phi
Integer(int32),Intent(out) :: ierr, i_last_good

! Local
Integer(int32), Parameter :: n = 2
Real(real64) :: y(n),x,dx,xout(nsteps+1),yout(n,nsteps+1)
Integer(int32) :: ierr_rk45, i_last_good_rk45

r = 0._real64
z = 0._real64
phi = 0._real64
dx = dz
y(1) = rstart
y(2) = phistart
x = zstart
Call rk45_fixed_step_integrate(bfield,y,n,x,dx,nsteps,fl_derivs_fun_dz,yout,xout,ierr_rk45,i_last_good_rk45)
r(1:nsteps+1)   = yout(1,1:nsteps+1)
phi(1:nsteps+1) = yout(2,1:nsteps+1)
z(1:nsteps+1)   = xout
ierr = ierr_rk45
i_last_good = i_last_good_rk45

End Subroutine follow_fieldlines_rzphi_dz_1pt

!-----------------------------------------------------------------------------
!+ Follows fieldlines dphi in cylindrical coords (from multiple points)
!-----------------------------------------------------------------------------
Subroutine follow_fieldlines_rzphi_Npts(bfield,rstart,zstart,phistart,Npts,dphi,nsteps,r,z,phi,ierr,i_last_good,verbose)
!
! Description: 
!  Follows fieldlines by integrating along toroidal angle in cylindrical coordinates. 
!
!  Integration method is RK45 with fixed step size
!
! Input:
!  rstart, zstart, phistart : real64(Npts) : Launching points for fieldlines (m,m,radians)
!  Npts : int32 : Number of starting points
!  dphi : real64 : Integration step size (radians)
!  nsteps : int32 : Number of integration steps to take
! verbose : logcl : Display fl counter  
! Output:
!  r,z,phi : real64(Npts,nsteps+1) : Field line trajectories (m,m,radians)
!  ierr : int32(Npts) : Error flag for each fl. 0 indicates no error, 1 indicates an error from bfield_geq_bicub
!  i_last_good : int32(Npts) : If an error occured for the fl, this array gives the last 'good' index, i.e., before the 
!                             error occured
!
! Calls:
!  Subroutine rk45_fixed_step_integrate
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/20/2011  Original code
!
! Author(s): J.D Lore - 04/20/2011
!
! Modules used:
Use kind_mod, Only: real64, int32
Implicit None
! Input/output                      !See above for descriptions
  Type(bfield_type), Intent(In) :: bfield
Integer(int32),Intent(in) :: Npts, nsteps
Real(real64),Intent(in),Dimension(Npts) :: rstart(Npts),zstart(Npts),phistart(Npts)
Real(real64),Intent(in) :: dphi
Logical, Intent(in), Optional :: verbose
Real(real64),Intent(out),Dimension(Npts,nsteps+1) :: r,z,phi
Integer(int32),Intent(out) :: ierr(Npts), i_last_good(Npts)
! Local variables
Integer(int32), Parameter :: n = 2
Real(real64) :: y(n),x,dx,xout(nsteps+1),yout(n,nsteps+1)
Integer(int32) :: ierr_rk45, i_last_good_rk45, ipt
Logical :: isverbose
!- End of header -------------------------------------------------------------

isverbose = .false.
If (Present(verbose)) isverbose = verbose
r = 0._real64
z = 0._real64
phi = 0._real64
dx = dphi
Do ipt = 1,Npts
  If (isverbose) Write(*,'(2(a,i0))') ' Following fl ',ipt,' of ',Npts
!  Call follow_fieldlines_rzphi_1pt(rstart(ipt),zstart(ipt),phistart(ipt),dphi,nsteps,r,z,phi,ierr,i_last_good)
  y(1) = rstart(ipt)
  y(2) = zstart(ipt)
  x = phistart(ipt)
  Call rk45_fixed_step_integrate(bfield,y,n,x,dx,nsteps,fl_derivs_fun,yout,xout,ierr_rk45,i_last_good_rk45)
  r(ipt,1:nsteps+1)   = yout(1,1:nsteps+1)
  z(ipt,1:nsteps+1)   = yout(2,1:nsteps+1)
  phi(ipt,1:nsteps+1) = xout
  ierr(ipt) = ierr_rk45
  i_last_good(ipt) = i_last_good_rk45
Enddo

EndSubroutine follow_fieldlines_rzphi_Npts

!-----------------------------------------------------------------------------
!+ Follows fieldlines dz in cylindrical coords (from multiple points)
!-----------------------------------------------------------------------------
Subroutine follow_fieldlines_rzphi_dz_Npts(bfield,rstart,zstart,phistart,Npts,dz,nsteps,r,z,phi,ierr,i_last_good,verbose)
  !
! Description: 
!  Follows fieldlines by integrating along toroidal angle in cylindrical coordinates. 
!
!  Integration method is RK45 with fixed step size
!
! Input:
!  rstart, zstart, phistart : real64(Npts) : Launching points for fieldlines (m,m,radians)
!  Npts : int32 : Number of starting points
!  dphi : real64 : Integration step size (radians)
!  nsteps : int32 : Number of integration steps to take
! verbose : logcl : Display fl counter  
! Output:
!  r,z,phi : real64(Npts,nsteps+1) : Field line trajectories (m,m,radians)
!  ierr : int32(Npts) : Error flag for each fl. 0 indicates no error, 1 indicates an error from bfield_geq_bicub
!  i_last_good : int32(Npts) : If an error occured for the fl, this array gives the last 'good' index, i.e., before the 
!                             error occured
!
! Calls:
!  Subroutine rk45_fixed_step_integrate
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/20/2011  Original code
!
! Author(s): J.D Lore - 04/20/2011
!
! Modules used:
Use kind_mod, Only: real64, int32
Implicit None
! Input/output                      !See above for descriptions
Type(bfield_type), Intent(In) :: bfield
Integer(int32),Intent(in) :: Npts, nsteps
Real(real64),Intent(in),Dimension(Npts) :: rstart(Npts),zstart(Npts),phistart(Npts)
Real(real64),Intent(in) :: dz
Logical, Intent(in), Optional :: verbose
Real(real64),Intent(out),Dimension(Npts,nsteps+1) :: r,z,phi
Integer(int32),Intent(out) :: ierr(Npts), i_last_good(Npts)
! Local variables
Integer(int32), Parameter :: n = 2
Real(real64) :: y(n),x,dx,xout(nsteps+1),yout(n,nsteps+1)
Integer(int32) :: ierr_rk45, i_last_good_rk45, ipt
Logical :: isverbose
!- End of header -------------------------------------------------------------

isverbose = .false.
If (Present(verbose)) isverbose = verbose
r = 0._real64
z = 0._real64
phi = 0._real64
dx = dz
Do ipt = 1,Npts
  If (isverbose) Write(*,'(2(a,i0))') ' Following fl ',ipt,' of ',Npts
  y(1) = rstart(ipt)
  y(2) = phistart(ipt)
  x = zstart(ipt)
  Call rk45_fixed_step_integrate(bfield,y,n,x,dx,nsteps,fl_derivs_fun_dz,yout,xout,ierr_rk45,i_last_good_rk45)
  r(ipt,1:nsteps+1)   = yout(1,1:nsteps+1)
  phi(ipt,1:nsteps+1) = yout(2,1:nsteps+1)
  z(ipt,1:nsteps+1)   = xout
  ierr(ipt) = ierr_rk45
  i_last_good(ipt) = i_last_good_rk45
Enddo

EndSubroutine follow_fieldlines_rzphi_dz_Npts

!-----------------------------------------------------------------------------
!+ Routines field line equation derivatives (dphi)
!-----------------------------------------------------------------------------
Subroutine fl_derivs_fun(bfield,n,phi,RZ,df,ierr)
!
! Description: 
!  Evaluates field line deriviatives (based on bfield%method). Should be easy to generalize to different
!  integration methods, but right now it is hard-coded to assume two simultaneous equations evaluated at 1 pt.
!
! Input:
!  n : int32 : Number of equations
!  phi : real64 : Integration variable (radians)
!  RZ : real64(n) : Evaluation points (solution vector)  RZ(1) = R, RZ(2) = Z in meters
! 
! Output:
!  df : real64(n) : derivative evaluation
!  ierr : int32 : Error flag (0 = no error)
!
! Calls:
!  Subroutine bfield_geq_bicub
!  Subroutine bfield_bs_cyl
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/20/2011  Original code
!
! Author(s): J.D Lore - 04/20/2011
!
! Modules used:
Use kind_mod, Only: real64, int32
Use bfield, Only : calc_B_rzphi_general
Implicit None
Type(bfield_type), Intent(In) :: bfield
Real(real64), Intent(In) :: phi
Integer(int32), Intent(In) :: n
Integer(int32), Intent(Out) :: ierr
Real(real64), Intent(In), Dimension(n) :: RZ
Real(real64), Intent(Out), Dimension(n) :: df

Integer(int32),Parameter :: Npts = 1
Real(real64) :: bval(Npts,3)
Integer(int32) :: ierr_b
Real(real64) :: Bz, Br, Bphi, phi_tmp(Npts)
!- End of header -------------------------------------------------------------

bval = 0._real64
phi_tmp = phi
Call calc_B_rzphi_general(bfield,RZ(1),RZ(2),phi_tmp,Npts,bval(1,1),bval(1,2),bval(1,3),ierr_b) 
Br   = bval(1,1)
Bz   = bval(1,2)
Bphi = bval(1,3)

If (ierr_b .ne. 0) Then
  ierr = 1
  df = 0._real64
  Return
Else
  ierr = 0
Endif
df(1) = RZ(1)*Br/Bphi
df(2) = RZ(1)*Bz/Bphi

End Subroutine fl_derivs_fun


!-----------------------------------------------------------------------------
!+ Routines field line equation derivatives (dz)
!-----------------------------------------------------------------------------
Subroutine fl_derivs_fun_dz(bfield,n,Z,RP,df,ierr)
!
! Description: 
!  Evaluates field line deriviatives (based on bfield%method). Should be easy to generalize to different
!  integration methods, but right now it is hard-coded to assume two simultaneous equations evaluated at 1 pt.
!
! Input:
!  n : int32 : Number of equations
!  z : real64 : Integration variable (m)
!  RP : real64(n) : Evaluation points (solution vector)  RP(1) = R, RP(2) = Phi in radians
! 
! Output:
!  df : real64(n) : derivative evaluation
!  ierr : int32 : Error flag (0 = no error)
!
! Calls:
!  Subroutine bfield_geq_bicub
!  Subroutine bfield_bs_cyl
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/20/2011  Original code
!
! Author(s): J.D Lore - 04/20/2011
!
! Modules used:
Use kind_mod, Only: real64, int32
Use biotsavart_module, Only : bfield_bs_cyl
Implicit None
Type(bfield_type), Intent(In) :: bfield
Real(real64), Intent(In) :: Z
Integer(int32), Intent(In) :: n
Integer(int32), Intent(Out) :: ierr
Real(real64), Intent(In), Dimension(n) :: RP
Real(real64), Intent(Out), Dimension(n) :: df

Real(real64) :: Bz, Br, Bphi
!- End of header -------------------------------------------------------------

ierr = 0
If (bfield%method == 6) Then     ! just coils
  Call bfield_bs_cyl(RP(1),RP(2),Z,bfield%coil,Br,Bphi,Bz)
Else
  Write(*,*) 'Unknown bfield%method in fl_derivs_fun_dz'
  Write(*,*) 'The following are supported'
  Write(*,*) '(6) just coils'
  stop
Endif
df(1) = Br/Bz
df(2) = Bphi/(Bz*RP(1))

End Subroutine fl_derivs_fun_dz


!-----------------------------------------------------------------------------
!+ Main routine for RK45 fixed step integration
!-----------------------------------------------------------------------------
Subroutine rk45_fixed_step_integrate(bfield,y0,n,x0,dx,nsteps,odefun,yout,xout,ierr,i_last_good)
!
! Description: 
!  Should be a general implementation of RK45 fixed step integration
!
! Input:
!  y0     : real64(n) : initial values
!  n      : int32    : number of initial values
!  x0     : real64    : Location of initial values
!  dx     : real64    : Step size
!  nsteps : int32    : Number of integration steps
!  odefun : External function that evaluates derivatives
! 
! Output:
!  yout : real64(n,nsteps+1) : Solution
!  xout : real64(nsteps+1)   : Solution evaluation locations
!  ierr : int32             : Error flag (0 = no error)
!  i_last_good : int32      : Index of last good evaluation before error
!
! Calls:
!  Subroutine odefun
!  Subroutine rk4_core
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/20/2011  Original code
!
! Author(s): J.D Lore - 04/20/2011
!
! Modules used:
Use kind_mod, Only: real64, int32
Implicit None
Type(bfield_type), Intent(In) :: bfield
Real(real64), Intent(In), Dimension(n) :: y0
Integer(int32), Intent(In) :: n, nsteps
Real(real64), Intent(In) :: x0, dx
Real(real64), Intent(Out), Dimension(n,nsteps+1) :: yout
Real(real64), Intent(Out), Dimension(nsteps+1) :: xout
Integer(int32), Intent(Out) :: ierr, i_last_good

Integer(int32) :: i, ierr_odefun, ierr_rk4core
Real(real64), Dimension(n) :: y, dydx, ytmp
Real(real64) :: x

Interface
  Subroutine odefun(bfield,n,x,y,dydx,ierr)
    Use kind_mod, Only: int32, real64
    Use bfield_typedef, Only : bfield_type
    Type(bfield_type), Intent(In) :: bfield
    Real(real64), Intent(In) :: x
    Integer(int32), Intent(In) :: n
    Integer(int32), Intent(Out) :: ierr
    Real(real64), Intent(In), Dimension(n) :: y
    Real(real64), Intent(Out), Dimension(n) :: dydx
  End Subroutine odefun
End Interface
!- End of header -------------------------------------------------------------

! Store initial point
yout(:,:) = 0._real64
xout(:) = 0._real64
yout(:,1) = y0
xout(1) = x0

y = y0
x = x0
ierr = 0
i_last_good = nsteps+1
Do i=1,nsteps
  Call odefun(bfield,n,x,y,dydx,ierr_odefun)
  If (ierr_odefun == 1) Then
    ierr = 1
    i_last_good = i
    Return
  Endif
  Call rk4_core(bfield,y,dydx,n,x,dx,odefun,ytmp,ierr_rk4core)
  If (ierr_rk4core == 1) Then
    ierr = 1
    i_last_good = i
    Return
  Endif

  x = x + dx

  yout(:,i+1) = ytmp
  xout(i+1) = x
  y = ytmp
Enddo

End Subroutine rk45_fixed_step_integrate

!-----------------------------------------------------------------------------
!+ Main routine for RK45 fixed step integration with diffusion
!-----------------------------------------------------------------------------
Subroutine rk45_fixed_step_integrate_diffuse(bfield,y0,n,x0,dx,nsteps,odefun,yout,xout,ierr,i_last_good,dmag)
!
! Description: 
!  Should be a general implementation of RK45 fixed step integration
!
! Input:
!  y0     : real64(n) : initial values
!  n      : int32    : number of initial values
!  x0     : real64    : Location of initial value
!  dx     : real64    : Step size
!  dmag   : real64    : Magnetic diffusivity (m^2/m)
!  nsteps : int32    : Number of integration steps
!  odefun : External function that evaluates derivatives
! 
! Output:
!  yout : real64(n,nsteps+1) : Solution
!  xout : real64(nsteps+1)   : Solution evaluation locations
!  ierr : int32             : Error flag (0 = no error)
!  i_last_good : int32      : Index of last good evaluation before error
!
! Calls:
!  Subroutine odefun
!  Subroutine rk4_core
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     04/20/2011  Original code
!
! Author(s): J.D Lore - 04/20/2011
!
! Modules used:
Use kind_mod, Only: real64, int32
Use bfield, Only : calc_B_rzphi_general
Use phys_const, Only : pi
Implicit None

Type(bfield_type), Intent(In) :: bfield
Real(real64), Intent(In), Dimension(n) :: y0
Integer(int32), Intent(In) :: n, nsteps
Real(real64), Intent(In) :: x0, dx, dmag
Real(real64), Intent(Out), Dimension(n,nsteps+1) :: yout
Real(real64), Intent(Out), Dimension(nsteps+1) :: xout
Integer(int32), Intent(Out) :: ierr, i_last_good
Integer(int32) :: i, ierr_odefun, ierr_rk4core
Real(real64), Dimension(n) :: y, dydx, ytmp
Real(real64) :: x, RZ(2), perpdir1(3), perpdir2(3), alpha, dca, dsa, delta_x, dL
Real(real64) :: bval(1,3), phi_tmp(1), phi, rnum
Integer(int32) :: ierr_b
Real(real64) :: Bz, Br, Bphi

Interface
  Subroutine odefun(bfield,n,x,y,dydx,ierr)
    Use kind_mod, Only: int32, real64
    Use bfield_typedef, Only : bfield_type
    Type(bfield_type), Intent(In) :: bfield
    Real(real64), Intent(In) :: x
    Integer(int32), Intent(In) :: n
    Integer(int32), Intent(Out) :: ierr
    Real(real64), Intent(In), Dimension(n) :: y
    Real(real64), Intent(Out), Dimension(n) :: dydx
  End Subroutine odefun
End Interface
!- End of header -------------------------------------------------------------

! Store initial point
yout(:,:) = 0._real64
xout(:) = 0._real64
yout(:,1) = y0
xout(1) = x0

y = y0
x = x0
ierr = 0
i_last_good = nsteps+1
Do i=1,nsteps
  Call odefun(bfield,n,x,y,dydx,ierr_odefun)
  If (ierr_odefun == 1) Then
    ierr = 1
    i_last_good = i
    Return
  Endif
  Call rk4_core(bfield,y,dydx,n,x,dx,odefun,ytmp,ierr_rk4core)
  If (ierr_rk4core == 1) Then
    ierr = 1
    i_last_good = i
    Return
  Endif

  x = x + dx

  !
  ! Diffuse
  ! 

  RZ(1) = ytmp(1)
  RZ(2) = ytmp(2)
  phi = x
  bval = 0._real64
  ierr_b = 0
  phi_tmp(1) = phi
  Call calc_B_rzphi_general(bfield,RZ(1),RZ(2),phi_tmp,1,bval(1,1),bval(1,2),bval(1,3),ierr_b) 
  Br   = bval(1,1)
  Bz   = bval(1,2)
  Bphi = bval(1,3)
  If (ierr_b .ne. 0)Then
    ierr = 1
    i_last_good = i
    Return
  Endif
  
  dL = sqrt(ytmp(1)*ytmp(1) + y(1)*y(1) - 2._real64*ytmp(1)*y(1)*cos(dx) & 
          + ytmp(2)*ytmp(2) + y(2)*y(2) - 2._real64*ytmp(2)*y(2))
  
  ! B cross z^hat
  perpdir1(1) = Bphi   !r 
  perpdir1(2) = -Br    !phi
  perpdir1(3) = 0.d0   !z
  perpdir1 = perpdir1/Sqrt(bphi*bphi + br*br)
  
  ! B cross r^hat
  perpdir2(1) = 0.d0
  perpdir2(2) = bz
  perpdir2(3) = -bphi
  perpdir2 = perpdir2/Sqrt(bz*bz + bphi*bphi)

  Call Random_number(rnum)
  
  alpha = 2.d0*pi*(-1.d0 + 2.d0*rnum) ! -2pi to 2pi kick    
  dca = Cos(alpha)
  dsa = Sin(alpha)
  
  delta_x = Sqrt(dmag*dL)
  
  
  ytmp(1)   = ytmp(1)   + delta_x*(dca*perpdir1(1) + dsa*perpdir2(1))
  x         = x         + delta_x*(dca*perpdir1(2) + dsa*perpdir2(2))
  ytmp(2)   = ytmp(2)   + delta_x*(dca*perpdir1(3) + dsa*perpdir2(3))

  yout(:,i+1) = ytmp
  xout(i+1) = x
  y = ytmp


Enddo





End Subroutine rk45_fixed_step_integrate_diffuse


!-----------------------------------------------------------------------------
!+ RK4 stepper
!-----------------------------------------------------------------------------
Subroutine rk4_core(bfield,y,dydx_in,n,x,dx,odefun,yout,ierr)
! Advance y(x) to yout=y(x+dx) given dydx(x) using 
! the RK4 method. y, dydx_in are vectors of length n, 
! function odefun(n,x,y,dxdy) returns dydx. 
! JDL 5/2012 
Use kind_mod, Only: real64, int32
Implicit None
Type(bfield_type), Intent(In) :: bfield
Real(real64), Intent(In), Dimension(n) :: y, dydx_in
Real(real64), Intent(Out), Dimension(n) :: yout
Real(real64), Intent(In) :: x, dx
Integer(int32), Intent(In) :: n
Integer(int32), Intent(Out) :: ierr

Real(real64), Parameter :: TWO = 2._real64
Real(real64), Parameter :: SIX = 6._real64
Integer(int32) :: ierr_odefun
Real(real64), Dimension(n) :: d1,d2,d3,d4,dydx

Interface
  Subroutine odefun(bfield,n,x,y,dydx,ierr)
    Use kind_mod, Only: int32, real64
    Use bfield_typedef, Only : bfield_type
    Type(bfield_type), Intent(In) :: bfield
    Real(real64), Intent(In) :: x
    Integer(int32), Intent(In) :: n
    Integer(int32), Intent(Out) :: ierr
    Real(real64), Intent(In), Dimension(n) :: y
    Real(real64), Intent(Out), Dimension(n) :: dydx
  End Subroutine odefun
End Interface
!- End of header -------------------------------------------------------------

dydx = dydx_in ! Do not want to overwrite input
! First step (uses supplied derivatives)
d1 = dx*dydx
ierr = 0
! Second step
Call odefun(bfield,n,x+dx/TWO,y+d1/TWO,dydx,ierr_odefun)
If (ierr_odefun == 1) Then
  ierr = 1
  yout = 0._real64
  Return
Endif
d2 = dx*dydx

! Third step
Call odefun(bfield,n,x+dx/TWO,y+d2/TWO,dydx,ierr_odefun)
If (ierr_odefun == 1) Then
  ierr = 1
  yout = 0._real64
  Return
Endif
d3 = dx*dydx

! Fourth step
Call odefun(bfield,n,x+dx,y+d3,dydx,ierr_odefun)
If (ierr_odefun == 1) Then
  ierr = 1
  yout = 0._real64
  Return
Endif
d4 = dx*dydx

yout = y + (d1 + TWO*d2 + TWO*d3 + d4)/SIX

Return
End Subroutine rk4_core


End Module fieldline_follow_mod

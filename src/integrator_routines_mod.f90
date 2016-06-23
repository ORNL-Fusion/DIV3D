!-----------------------------------------------------------------------------
!+ Integrator core and stepping routines
!-----------------------------------------------------------------------------
Module integrator_routines_mod
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
!+
!------------------------------------------------------------------------------------------
Subroutine rk45_adp_step_integrate(y0,n,x1,x2,rtol,atol,dx1,dxmin,nok,nbad,odefun,yout,xout,nmax_step)
Use kind_mod
Implicit None
Real(rknd), Intent(In), Dimension(n) :: y0
Integer(iknd), Intent(In) :: n, nmax_step
Real(rknd), Intent(In) :: x1, x2, rtol, atol, dx1, dxmin
Real(rknd), Intent(Out), Dimension(n) :: yout
Real(rknd), Intent(Out) :: xout
Integer(iknd), Intent(Out) :: nok, nbad

Integer(iknd) :: i
Real(rknd) :: x, dx, dxused, dxnext, xtmp
Real(rknd), Dimension(n) :: y, yscal, dydx, ytmp
Real(rknd ), Parameter :: tinyx = tiny(1._rknd)
External odefun

x = x1
dx = sign(dx1,x2-x1)
nok = 0
nbad = 0
y = y0

Do i = 1,nmax_step
  
!write(*,*) 'here'
  Call odefun(n,x,y,dydx)
  
  ! General form of error scaling
  yscal = abs(y) + abs(dx*dydx) + tinyx

  ! Don't allow stepsize to overshoot
  If ( (x+dx-x2)*(x+dx-x1) .gt. 0._rknd ) dx = x2-x
 
!  Call rk45_adp_step(y,dydx,n,x,dx,tol,yscal,dxused,dxnext,odefun,ytmp,xtmp)
!  Call rk45_adp_step2(y,dydx,n,x,dx,rtol,atol,dxused,dxnext,odefun,ytmp,xtmp)
  Call Dopri853_adp_step(y,dydx,n,x,dx,rtol,atol,dxused,dxnext,odefun,ytmp,xtmp)

  x = xtmp
  y = ytmp
  If ( dxused .eq. dx ) Then
    nok = nok + 1
  Else
    nbad = nbad + 1
  Endif

  ! Check if done
  If ( (x-x2)*(x2-x1) .ge. 0._rknd ) Then
    yout = ytmp
    xout = xtmp
    Return
  Endif

  If (abs(dxnext) .lt. dxmin) Then
    Write(6,*) 'stepsize smaller than minimum allowed in rk45_adp_step_integrate',dxmin
    Stop
  Endif

  dx = dxnext

Enddo

Write(6,*) 'Max number of steps exceeded in rk45_adp_step_integrate',nmax_step
Stop


End Subroutine rk45_adp_step_integrate
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


!------------------------------------------------------------------------------------------
!+
!------------------------------------------------------------------------------------------
Subroutine Dopri853_adp_step(y,dydx,n,x,dxtry,rtol,atol,dxused,dxnext,odefun,yout,xout)
! DP853 stepper routine. Attempts a stepsize of dxtry and, based on the error of that
! step, uses a stepsize dxused and recommends a new stepsize dxnext. 
! JDL 5/2012 
Use kind_mod
Implicit None
Real(rknd), Intent(In), Dimension(n) :: y, dydx
Real(rknd), Intent(In) :: x, dxtry, rtol, atol
Integer(iknd), Intent(In) :: n
Real(rknd), Intent(Out), Dimension(n) :: yout
Real(rknd), Intent(Out) :: dxused, dxnext, xout


Real(rknd) :: dx, err, err2, xnew, sk, scale, denom
Real(rknd), Dimension(n) :: ytemp, yerr, yerr2 !, dydxnext
Integer(iknd) :: iter, i, reject

Integer(iknd), Parameter :: iter_max = 100_iknd
Real(rknd), Parameter :: &
Safety = 0.9_rknd, &
errcon = 1.0e-4_rknd, & 
maxscale = 6._rknd, &
minscale = 0.333_rknd, &
beta = 0._rknd, &  ! Change to 0.04 or 0.08 for PI control. (also uncomment lines below)
alpha = -0.125_rknd + beta*0.2_rknd

External odefun

! Initial stepsize
dx = dxtry
reject = 0
! errold = errcon

Do iter = 1,iter_max 

  ! Take a step

!  ytemp=0._rknd
  Call Dopri853_core(y,dydx,n,x,dx,odefun,ytemp,yerr,yerr2)
!  write(*,*) '1',ytemp
!  Call Dopri5_core(y,dydx,n,x,dx,odefun,ytemp,yerr,dydxnext)
!  write(*,*) '2',ytemp
!  stop



  ! Evaluate scaled error
  err = 0._rknd
  err2 = 0._rknd
  Do i = 1,n
    sk = atol + rtol*Max(Abs(y(i)),Abs(ytemp(i)))
    err = err + (yerr2(i)/sk)*(yerr2(i)/sk)
    err2 = err2 + (yerr(i)/sk)*(yerr(i)/sk)
  Enddo

  denom = err+0.01_rknd*err2
  If (denom .le. 0._rknd) Then
    denom = 1._rknd
  Endif
  err = abs(dx)*err*Sqrt(1._rknd/(denom*Real(n,rknd)))


  If (err .gt. 1) Then  
    ! Step failed
    ! Error too large, reduce stepsize and re-evaluate
    scale = Max(safety*(err**alpha),minscale)
    dx = dx*scale

    ! Check for dx underflow
    xnew = x + dx
    If (xnew .eq. x) Then 
      Write(6,*) 'Stepsize (dx) underflow in rk4f_adp_step2'
    Endif
    reject = 1
  Else 
    ! Step succeeded
    ! Determine next step size
    If (err .eq. 0._rknd) Then
      scale = maxscale
    Else
!     scale = safety*(err**alpha)*(errold^beta)  USE FOR PI CONTROL (would need to implement errold too)
      scale = safety*(err**alpha)
      If (scale .lt. minscale) scale = minscale
      If (scale .gt. maxscale) scale = maxscale
    Endif 

    ! If last step was too large don't increase this time
    If (reject .eq. 1) Then
      dxnext = dx*Min(scale,1._rknd)
    Else
      dxnext = dx*scale
    Endif

!   errold = Max(err,errcon)
    reject = 0
    dxused = dx
    xout = x + dx
    yout = ytemp
    Return
  Endif

Enddo
End Subroutine Dopri853_adp_step
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!------------------------------------------------------------------------------------------
!+
!------------------------------------------------------------------------------------------
Subroutine Dopri853_core(y,dydx,n,x,dx,odefun,yout,yerr,yerr2)
! Routine to take a Dormand-Prince Runge-Kutta step. Coefficients from
! NR v3 3.02
! JDL 5/2012
Use kind_mod
Implicit None
Real(rknd), Intent(In), Dimension(n) :: y, dydx
Real(rknd), Intent(Out), Dimension(n) :: yout, yerr, yerr2
Real(rknd), Intent(In) :: x, dx
Integer(iknd), Intent(In) :: n
Real(rknd), Dimension(n) ::  ak2, ak3, ak4, ak5, ak6, &
ak7, ak8, ak9, ak10, ak11, ak12, ak13, ytemp

Real(rknd), Parameter :: &
C2  = 0.526001519587677318785587544488e-01_rknd, &
C3  = 0.789002279381515978178381316732e-01_rknd, &
C4  = 0.118350341907227396726757197510e+00_rknd, &
C5  = 0.281649658092772603273242802490e+00_rknd, &
C6  = 0.333333333333333333333333333333e+00_rknd, &
C7  = 0.25e+00_rknd, &
C8  = 0.307692307692307692307692307692e+00_rknd, &
C9  = 0.651282051282051282051282051282e+00_rknd, &
C10 = 0.6e+00_rknd, &
C11 = 0.857142857142857142857142857142e+00_rknd, &
C14 = 0.1e+00_rknd, &
C15 = 0.2e+00_rknd, &
C16 = 0.777777777777777777777777777778e+00_rknd, &
B1 = 5.42937341165687622380535766363e-2_rknd, &
B6 = 4.45031289275240888144113950566e0_rknd, &
B7 = 1.89151789931450038304281599044e0_rknd, &
B8 = -5.8012039600105847814672114227e0_rknd, &
B9 = 3.1116436695781989440891606237e-1_rknd, &
B10 = -1.52160949662516078556178806805e-1_rknd, &
B11 = 2.01365400804030348374776537501e-1_rknd, &
B12 = 4.47106157277725905176885569043e-2_rknd, &
BHH1 = 0.244094488188976377952755905512e+00_rknd, &
BHH2 = 0.733846688281611857341361741547e+00_rknd, &
BHH3 = 0.220588235294117647058823529412e-01_rknd, &
ER1 = 0.1312004499419488073250102996e-01_rknd, &
ER6 = -0.1225156446376204440720569753e+01_rknd, &
ER7 = -0.4957589496572501915214079952e+00_rknd, &
ER8 = 0.1664377182454986536961530415e+01_rknd, &
ER9 = -0.3503288487499736816886487290e+00_rknd, &
ER10 = 0.3341791187130174790297318841e+00_rknd, &
ER11 = 0.8192320648511571246570742613e-01_rknd, &
ER12 = -0.2235530786388629525884427845e-01_rknd, &
A21 = 5.26001519587677318785587544488e-2_rknd, &
A31 = 1.97250569845378994544595329183e-2_rknd, &
A32 = 5.91751709536136983633785987549e-2_rknd, &
A41 = 2.95875854768068491816892993775e-2_rknd, &
A43 = 8.87627564304205475450678981324e-2_rknd, &
A51 = 2.41365134159266685502369798665e-1_rknd, &
A53 = -8.84549479328286085344864962717e-1_rknd, &
A54 = 9.24834003261792003115737966543e-1_rknd, &
A61 = 3.7037037037037037037037037037e-2_rknd, &
A64 = 1.70828608729473871279604482173e-1_rknd, &
A65 = 1.25467687566822425016691814123e-1_rknd, &
A71 = 3.7109375e-2_rknd, &
A74 = 1.70252211019544039314978060272e-1_rknd, &
A75 = 6.02165389804559606850219397283e-2_rknd, &
A76 = -1.7578125e-2_rknd, &
A81 = 3.70920001185047927108779319836e-2_rknd, &
A84 = 1.70383925712239993810214054705e-1_rknd, &
A85 = 1.07262030446373284651809199168e-1_rknd, &
A86 = -1.53194377486244017527936158236e-2_rknd, &
A87 = 8.27378916381402288758473766002e-3_rknd, &
A91 = 6.24110958716075717114429577812e-1_rknd, &
A94 = -3.36089262944694129406857109825e0_rknd, &
A95 = -8.68219346841726006818189891453e-1_rknd, &
A96 = 2.75920996994467083049415600797e1_rknd, &
A97 = 2.01540675504778934086186788979e1_rknd, &
A98 = -4.34898841810699588477366255144e1_rknd, &
A101 = 4.77662536438264365890433908527e-1_rknd, &
A104 = -2.48811461997166764192642586468e0_rknd, &
A105 = -5.90290826836842996371446475743e-1_rknd, &
A106 = 2.12300514481811942347288949897e1_rknd, &
A107 = 1.52792336328824235832596922938e1_rknd, &
A108 = -3.32882109689848629194453265587e1_rknd, &
A109 = -2.03312017085086261358222928593e-2_rknd, &
A111 = -9.3714243008598732571704021658e-1_rknd, &
A114 = 5.18637242884406370830023853209e0_rknd, &
A115 = 1.09143734899672957818500254654e0_rknd, &
A116 = -8.14978701074692612513997267357e0_rknd, &
A117 = -1.85200656599969598641566180701e1_rknd, &
A118 = 2.27394870993505042818970056734e1_rknd, &
A119 = 2.49360555267965238987089396762e0_rknd, &
A1110 = -3.0467644718982195003823669022e0_rknd, &
A121 = 2.27331014751653820792359768449e0_rknd, &
A124 = -1.05344954667372501984066689879e1_rknd, &
A125 = -2.00087205822486249909675718444e0_rknd, &
A126 = -1.79589318631187989172765950534e1_rknd, &
A127 = 2.79488845294199600508499808837e1_rknd, &
A128 = -2.85899827713502369474065508674e0_rknd, &
A129 = -8.87285693353062954433549289258e0_rknd, &
A1210 = 1.23605671757943030647266201528e1_rknd, &
A1211 = 6.43392746015763530355970484046e-1_rknd, &
A141 = 5.61675022830479523392909219681e-2_rknd, &
A147 = 2.53500210216624811088794765333e-1_rknd, &
A148 = -2.46239037470802489917441475441e-1_rknd, &
A149 = -1.24191423263816360469010140626e-1_rknd, &
A1410 = 1.5329179827876569731206322685e-1_rknd, &
A1411 = 8.20105229563468988491666602057e-3_rknd, &
A1412 = 7.56789766054569976138603589584e-3_rknd, &
A1413 = -8.298e-3_rknd, &
A151 = 3.18346481635021405060768473261e-2_rknd, &
A156 = 2.83009096723667755288322961402e-2_rknd, &
A157 = 5.35419883074385676223797384372e-2_rknd, &
A158 = -5.49237485713909884646569340306e-2_rknd, &
A1511 = -1.08347328697249322858509316994e-4_rknd, &
A1512 = 3.82571090835658412954920192323e-4_rknd, &
A1513 = -3.40465008687404560802977114492e-4_rknd, &
A1514 = 1.41312443674632500278074618366e-1_rknd, &
A161 = -4.28896301583791923408573538692e-1_rknd, &
A166 = -4.69762141536116384314449447206e0_rknd, &
A167 = 7.68342119606259904184240953878e0_rknd, &
A168 = 4.06898981839711007970213554331e0_rknd, &
A169 = 3.56727187455281109270669543021e-1_rknd, &
A1613 = -1.39902416515901462129418009734e-3_rknd, &
A1614 = 2.9475147891527723389556272149e0_rknd, &
A1615 = -9.15095847217987001081870187138e0_rknd, &
D41 = -0.84289382761090128651353491142e+01_rknd, &
D46 = 0.56671495351937776962531783590e+00_rknd, &
D47 = -0.30689499459498916912797304727e+01_rknd, &
D48 = 0.23846676565120698287728149680e+01_rknd, &
D49 = 0.21170345824450282767155149946e+01_rknd, &
D410 = -0.87139158377797299206789907490e+00_rknd, &
D411 = 0.22404374302607882758541771650e+01_rknd, &
D412 = 0.63157877876946881815570249290e+00_rknd, &
D413 = -0.88990336451333310820698117400e-01_rknd, &
D414 = 0.18148505520854727256656404962e+02_rknd, &
D415 = -0.91946323924783554000451984436e+01_rknd, &
D416 = -0.44360363875948939664310572000e+01_rknd, &
D51 = 0.10427508642579134603413151009e+02_rknd, &
D56 = 0.24228349177525818288430175319e+03_rknd, &
D57 = 0.16520045171727028198505394887e+03_rknd, &
D58 = -0.37454675472269020279518312152e+03_rknd, &
D59 = -0.22113666853125306036270938578e+02_rknd, &
D510 = 0.77334326684722638389603898808e+01_rknd, &
D511 = -0.30674084731089398182061213626e+02_rknd, &
D512 = -0.93321305264302278729567221706e+01_rknd, &
D513 = 0.15697238121770843886131091075e+02_rknd, &
D514 = -0.31139403219565177677282850411e+02_rknd, &
D515 = -0.93529243588444783865713862664e+01_rknd, &
D516 = 0.35816841486394083752465898540e+02_rknd, &
D61 = 0.19985053242002433820987653617e+02_rknd, &
D66 = -0.38703730874935176555105901742e+03_rknd, &
D67 = -0.18917813819516756882830838328e+03_rknd, &
D68 = 0.52780815920542364900561016686e+03_rknd, &
D69 = -0.11573902539959630126141871134e+02_rknd, &
D610 = 0.68812326946963000169666922661e+01_rknd, &
D611 = -0.10006050966910838403183860980e+01_rknd, &
D612 = 0.77771377980534432092869265740e+00_rknd, &
D613 = -0.27782057523535084065932004339e+01_rknd, &
D614 = -0.60196695231264120758267380846e+02_rknd, &
D615 = 0.84320405506677161018159903784e+02_rknd, &
D616 = 0.11992291136182789328035130030e+02_rknd, &
D71 = -0.25693933462703749003312586129e+02_rknd, &
D76 = -0.15418974869023643374053993627e+03_rknd, &
D77 = -0.23152937917604549567536039109e+03_rknd, &
D78 = 0.35763911791061412378285349910e+03_rknd, &
D79 = 0.93405324183624310003907691704e+02_rknd, &
D710 = -0.37458323136451633156875139351e+02_rknd, &
D711 = 0.10409964950896230045147246184e+03_rknd, &
D712 = 0.29840293426660503123344363579e+02_rknd, &
D713 = -0.43533456590011143754432175058e+02_rknd, &
D714 = 0.96324553959188282948394950600e+02_rknd, &
D715 = -0.39177261675615439165231486172e+02_rknd, &
D716 = -0.14972683625798562581422125276e+03_rknd

External odefun

! 12 steps
ytemp = y + A21*dx*dydx
Call odefun(n,x+C2*dx,ytemp,ak2)
ytemp = y + dx*(A31*dydx + A32*ak2)
Call odefun(n,x+C3*dx,ytemp,ak3)
ytemp = y + dx*(A41*dydx + A43*ak3)
Call odefun(n,x+C4*dx,ytemp,ak4)
ytemp = y + dx*(A51*dydx + A53*ak3 + A54*ak4)
Call odefun(n,x+C5*dx,ytemp,ak5)
ytemp = y + dx*(A61*dydx + A64*ak4 + A65*ak5)
Call odefun(n,x+C6*dx,ytemp,ak6)
ytemp = y + dx*(A71*dydx + A74*ak4 + A75*ak5 + A76*ak6)
Call odefun(n,x+C7*dx,ytemp,ak7)
ytemp = y + dx*(A81*dydx + A84*ak4 + A85*ak5 + A86*ak6 + A87*ak7)
Call odefun(n,x+C8*dx,ytemp,ak8)
ytemp = y + dx*(A91*dydx + A94*ak4 + A95*ak5 + A96*ak6 + A97*ak7 + A98*ak8)
Call odefun(n,x+C9*dx,ytemp,ak9)
ytemp = y + dx*(A101*dydx + A104*ak4 + A105*ak5 + A106*ak6 + A107*ak7 + A108*ak8 + A109*ak9)
Call odefun(n,x+C10*dx,ytemp,ak10)
ytemp = y + dx*(A111*dydx + A114*ak4 + A115*ak5 + A116*ak6 + A117*ak7 + A118*ak8 + A119*ak9 + A1110*ak10)
Call odefun(n,x+C11*dx,ytemp,ak11)
ytemp = y + dx*(A121*dydx + A124*ak4 + A125*ak5 + A126*ak6 + A127*ak7 + A128*ak8 + A129*ak9 + A1210*ak10 + A1211*ak11)

Call odefun(n,x+dx,ytemp,ak12)
ak13 = B1*dydx + B6*ak6 + B7*ak7 + B8*ak8 + B9*ak9 + B10*ak10 + B11*ak11 + B12*ak12
yout = y + dx*ak13

! Estimate error
yerr = (ak13 - BHH1*dydx - BHH2*ak9 - BHH3*ak12)
yerr2 = (ER1*dydx + ER6*ak6 + ER7*ak7 + ER8*ak8 + ER9*ak9 + ER10*ak10 + ER11*ak11 + ER12*ak12)


Return
End Subroutine Dopri853_core
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

End Module integrator_routines_mod

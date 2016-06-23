!-----------------------------------------------------------------------------
!+ Module to read bfield (xdr) file and perform interpolation
!-----------------------------------------------------------------------------
module bfield_xdr
use fxdr_module
use kind_mod
Implicit none

integer :: nperio,ialfa,knull
integer :: k2, iald21
real(rknd) :: rnull,ronull,eta
real(rknd), dimension(:,:,:), allocatable :: brg
real(rknd), dimension(:,:,:), allocatable :: bfg
real(rknd), dimension(:,:,:), allocatable :: bzg

contains


!--------------------------------------------------
! Reads the xdr file containing bfield info
!--------------------------------------------------
Subroutine Readbgrid(filename,xdr_check,verbose)
Use kind_mod
Implicit None

character(len=300),intent(in) :: filename
logical,Intent(in) :: xdr_check, verbose

character(len=300) :: newfilename
logical :: l_val=.false., l_check=.false.
integer(iknd) :: ierr, ixdrs,iganz
integer(iknd) :: i1, i2
integer(iknd) :: iz,ir,i 
real(rknd) :: bmean, fnull,bfak,ampby0,bz0
logical :: debug=.false.

!logical :: verbose = .true.

!
! Read the xdr file (using list from Geiger)
!
If (verbose) Write(6,*) 'Reading xdr file ',Trim(Adjustl(filename))

ixdrs = initxdr(Trim(Adjustl(filename)),'r',l_val)
ierr = ixdrint(ixdrs,i1)
ierr = ixdrint(ixdrs,i2)
ierr = ixdrdouble(ixdrs,rnull)    !rmin + ronull (midpt in R)
ierr = ixdrdouble(ixdrs,ronull)   !0.5*(rmax-rmin)
ierr = ixdrdouble(ixdrs,eta)      !0.5*(zmax-zmin)/ronull
ierr = ixdrdouble(ixdrs,fnull)    !dummy (=1)
ierr = ixdrint(ixdrs,nperio)      !# field periods
ierr = ixdrint(ixdrs,ialfa)       !# intervals/period
ierr = ixdrint(ixdrs,knull)       !#of intervals within ronull
ierr = ixdrint(ixdrs,iganz)       !not used (=0)
ierr = ixdrdouble(ixdrs,ampby0)   !not used (=0)
ierr = ixdrdouble(ixdrs,bz0)      !not used (=0)
ierr = ixdrdouble(ixdrs,bfak)     !not used (=1?) (normalizer for B)

! Note that B indices are (i,j,k)= (1-ialfa/2+1, 1-2*knull+1,1-2*knull+1)

k2=2*knull+1
iald21=ialfa/2+1

! Sanity checks
If (fnull .ne. 1 .OR. iganz .ne. 0 .OR. ampby0 .ne. 0 .OR. bz0 .ne. 0 .OR. bfak .ne. 1) Then
  Write(6,*) 'Bad assumptions on xdr values!'
  Write(6,*) 'fnull: ',fnull
  Write(6,*) 'iganz: ',iganz
  Write(6,*) 'ampby0: ',ampby0
  Write(6,*) 'bz0: ',bz0
  Write(6,*) 'bfak: ',bfak
  Stop
Endif



If (debug) Then
  write(6,*) "Parameters read:"
  write(6,*) "i1 : ", i1
  write(6,*) "i2 : ", i2
  write(6,*) "rnull : ", rnull
  write(6,*) "ronull: ", ronull
  write(6,*) "eta   : ", eta
  write(6,*) "fnull : ", fnull
  write(6,*) "nperio: ", nperio
  write(6,*) "ialfa : ", ialfa
  write(6,*) "knull : ", knull
  write(6,*) "Number of R-grid points: ",k2
  write(6,*) "Number of z-grid points: ",k2
Endif


allocate(brg(iald21,k2,k2),bfg(iald21,k2,k2),bzg(iald21,k2,k2))
brg = 0.d+0 ; bfg = 0.d+0 ; bzg = 0.d+0
do iz=1,k2
  do ir=1,k2
    do i=1,iald21
      ierr = ixdrdouble(ixdrs,bfg(i,ir,iz))
    enddo
  enddo
enddo
do iz=1,k2
  do ir=1,k2
    do i=1,iald21
      ierr = ixdrdouble(ixdrs,brg(i,ir,iz))
    enddo
  enddo
enddo
do iz=1,k2
  do ir=1,k2
    do i=1,iald21
      ierr = ixdrdouble(ixdrs,bzg(i,ir,iz))
    enddo
  enddo
enddo
ierr = ixdrclose(ixdrs)

!
! Check for bad points in grid
!

If (xdr_check .eqv. .true.) Then
  If (verbose) Write(6,*) 'Checking for bad points'
  do i=1,iald21
    do iz=2,k2-1
      do ir=2,k2-1
        if( bfg(i,ir-1,iz) /= 0.d+0 .and. &
             bfg(i,ir+1,iz) /= 0.d+0 .and. &
             bfg(i,ir,iz-1) /= 0.d+0 .and. &
             bfg(i,ir,iz+1) /= 0.d+0 ) then
          bmean=0.25d+0*(bfg(i,ir-1,iz)+bfg(i,ir+1,iz)+ &
               bfg(i,ir,iz-1)+bfg(i,ir,iz+1))
          if( abs(bfg(i,ir,iz)) < 0.5*abs(bmean) .or. &
               abs(bfg(i,ir,iz)) > 1.5*abs(bmean) ) then
            write(6,*)"Point detected:"
            write(6,*)"bfg (",i,",",ir,",",iz,") = ",bfg(i,ir,iz)
            write(6,*)"bmean(",i,",",ir,",",iz,") = ",bmean
            brg(i,ir,iz)=0.25d+0*(brg(i,ir-1,iz)+brg(i,ir+1,iz)+ &
                 brg(i,ir,iz-1)+brg(i,ir,iz+1))
            bfg(i,ir,iz)=0.25d+0*(bfg(i,ir-1,iz)+bfg(i,ir+1,iz)+ &
                 bfg(i,ir,iz-1)+bfg(i,ir,iz+1))
            bzg(i,ir,iz)=0.25d+0*(bzg(i,ir-1,iz)+bzg(i,ir+1,iz)+ &
                 bzg(i,ir,iz-1)+bzg(i,ir,iz+1))
            l_check=.true.
          endif
        endif
      enddo
    enddo
  enddo
Else 
  If (verbose) Write(6,'(/A/)') ' WARNING: >>>>>>>>>>>>>>>> Skipping bfield check of xdr file!.  <<<<<<<<<<<<<<<<<<<<'
Endif


!
! Use average of surrounding points, if necessary
!

if(l_check) then
  Write(*,*) 'Bad points found'
  write(6,*)"Write corrected field to fieldn.data"
  newfilename = 'fieldn.data'
  ixdrs = initxdr(newfilename,'w',l_val)
  ierr = ixdrint(ixdrs,i1)
  ierr = ixdrint(ixdrs,i2)
  ierr = ixdrdouble(ixdrs,rnull)
  ierr = ixdrdouble(ixdrs,ronull)
  ierr = ixdrdouble(ixdrs,eta)
  ierr = ixdrdouble(ixdrs,fnull)
  ierr = ixdrint(ixdrs,nperio)
  ierr = ixdrint(ixdrs,ialfa)
  ierr = ixdrint(ixdrs,knull)
  ierr = ixdrint(ixdrs,iganz)
  ierr = ixdrdouble(ixdrs,ampby0)
  ierr = ixdrdouble(ixdrs,bz0)
  ierr = ixdrdouble(ixdrs,bfak)
  do iz=1,k2
    do ir=1,k2
      do i=1,iald21
        ierr = ixdrdouble(ixdrs,bfg(i,ir,iz))
      enddo
    enddo
  enddo
  do iz=1,k2
    do ir=1,k2
      do i=1,iald21
        ierr = ixdrdouble(ixdrs,brg(i,ir,iz))
      enddo
    enddo
  enddo
  do iz=1,k2
    do ir=1,k2
      do i=1,iald21
        ierr = ixdrdouble(ixdrs,bzg(i,ir,iz))
      enddo
    enddo
  enddo
  ierr = ixdrclose(ixdrs)
  stop
else
  If (verbose .AND. xdr_check) Write(*,*) 'No bad points found in xdr file'
endif

end subroutine readbgrid


!***********************************************************************
subroutine bint(xvec,bval,idiv)
!***********************************************************************
!  interpolation routine
!  for   linear and  toroidal magnetic fields without derivatives,
!               ***                           *******
!  "welters"  magnetic fields for program "gourdon" .
!  original program from Lotz ==>                zol:weg.grad2121
!                                    modified by jag, 27.Nov.1995
!***********************************************************************
! input : xvec(3), where
!
!      xvec(1)  = R          xvec(2) = phi         xvec(3) = z
!
! output: bval(3), where
!
!      bval(1)  = br         bval(2)  = bphi       bval(3)  = bz
!      idiv =   0     if bmod.ne.0 and xvec inside of int.region
!               1     if bmod.eq.0
!                  or if xvec is out of the interpolation region
!***********************************************************************

Use kind_mod
Implicit none

!********* list of variables *********************

Real(rknd), intent(in) :: xvec(3)
Real(rknd), intent(out) :: bval(3)
Integer(iknd), intent(out) :: idiv

!********* list of local variables ***************

Real(rknd) :: bx,by,bz,bri,bsi,bfi
Real(rknd) :: phi,cosphi,sinphi,cosphi0,sinphi0
Real(rknd) :: ray,dx,rq,v1,v2,sq, dphi
Real(rknd) :: dr,dz
Real(rknd) :: fq,fq1,fq1v1,fqv2, dphi_period
Real(rknd) :: ar,as,af,a1,a2,a3,a4
Real(rknd), dimension(6) :: bxt,byt,bzt,xpt,ypt,zpt
Real(rknd) :: x,y,z,xeval,yeval,zeval

integer(iknd) :: i,ng_phi, ng_rz, ig_ax,  nf, nc_phi
integer(iknd) :: lf,lf1,lr,lr1,ls,ls1,ls11,ls2,ls21

Real(rknd), parameter :: pi = 3.1415926535897932384626433832795_rknd

!*************************************************

integer(iknd) :: icall = 0
save :: icall,ng_rz,ng_phi,ig_ax,dr,dz,dx,dphi,dphi_period

!*************************************************

If ( icall .eq. 0 ) Then
  ng_rz = k2       
  ng_phi = iald21
  ig_ax = knull + 1
  nc_phi = ialfa     ! number of cells for full period
  dphi_period = 2._rknd*pi/Real(nperio,rknd)

  ! Set grid step sizes
  dphi = dphi_period/Real(nc_phi,rknd)
  dr = ronull/Real(knull,rknd)
  dz = dr*eta
  dx = 0.2d0*dr ! Step size for points around evalution point +/i

  icall = 1
Endif

! Convert to cartesian coordinates
cosphi0 = cos(xvec(2))
sinphi0 = sin(xvec(2))
xeval = xvec(1)*cosphi0
yeval = xvec(1)*sinphi0
zeval = xvec(3)

!*************************************************
idiv = 0
!*************************************************

!
! Set interpolation grid points
!
xpt(1) = xeval + dx
ypt(1) = yeval
zpt(1) = zeval

xpt(2) = xeval - dx
ypt(2) = yeval
zpt(2) = zeval

xpt(3) = xeval
ypt(3) = yeval + dx
zpt(3) = zeval

xpt(4) = xeval
ypt(4) = yeval - dx
zpt(4) = zeval

xpt(5) = xeval
ypt(5) = yeval
zpt(5) = zeval + dx

xpt(6) = xeval
ypt(6) = yeval
zpt(6) = zeval - dx

Do i=1,6
  x = xpt(i)
  y = ypt(i)
  z = zpt(i)

  ! Map this point to the fp
  phi    = atan2(y,x)
  phi = phi + (pi - sign(pi,phi))
  nf     = floor(phi/dphi_period)
  phi      = phi - nf*dphi_period
  ray    = sqrt(x*x + y*y)
  cosphi  = x/ray
  sinphi  = y/ray

  ! Find the grid indices around this point for interpolation
  af     = phi/dphi + 1._rknd
  lf     = af
  fq     = af - lf
  lf1    = lf + 1
  if ( lf1 .gt. ialfa )  lf1 = 1
  ar     = (ray-rnull)/dr + ig_ax
  lr     = Floor(ar)
  rq     = ar - lr
  lr1    = lr + 1
  as     = z/dz + ig_ax
  ls     = Floor(as)
  sq     = as - ls
!*************************************************
  if (idiv .eq. 1 .or. lr*(lr-ng_rz).ge.0 .or. ls*(ls-ng_rz).ge.0)then
    idiv = 1  ! Outside grid
  else
    v1     = 1._rknd
    v2     = 1._rknd
    ls1    = ls
    ls11   = ls1 + 1
    ls2    = ls1
    ls21   = ls11 
    If ( lf .gt. ng_phi ) Then ! apply stellarator symmetry to the vectors
      lf     = ialfa - lf + 2
      v1     = -1._rknd
      ls1    = -ls1  + ng_rz + 1
      ls11   = -ls11 + ng_rz + 1
    Endif
    if ( lf1 .gt. ng_phi ) then
      lf1    = ialfa - lf1 + 2
      v2     = -1._rknd
      ls2    = -ls2  + ng_rz + 1
      ls21   = -ls21 + ng_rz + 1
    endif
    !*************************************************
    !  interpolation : bilinear interpolation in r and z
    !                  and linear in phi for each point
    !*************************************************
    fq1   = 1._rknd - fq  
    fqv2  = fq*v2
    fq1v1 = fq1*v1
    a1    = (1._rknd - rq)*(1._rknd - sq)
    a2    = (1._rknd - rq)*sq
    a3    = rq*(1._rknd - sq)
    a4    = rq*sq
    !*************************************************

    bri =  a1*( fq1v1*brg(lf ,lr ,ls1 ) +fqv2*brg(lf1,lr ,ls2 ) ) &
          +a2*( fq1v1*brg(lf ,lr ,ls11) +fqv2*brg(lf1,lr ,ls21) ) &
          +a3*( fq1v1*brg(lf ,lr1,ls1 ) +fqv2*brg(lf1,lr1,ls2 ) ) &
          +a4*( fq1v1*brg(lf ,lr1,ls11) +fqv2*brg(lf1,lr1,ls21) )
    bsi =  a1*( fq1  *bzg(lf ,lr ,ls1 ) +fq  *bzg(lf1,lr ,ls2 ) ) &
          +a2*( fq1  *bzg(lf ,lr ,ls11) +fq  *bzg(lf1,lr ,ls21) ) &
          +a3*( fq1  *bzg(lf ,lr1,ls1 ) +fq  *bzg(lf1,lr1,ls2 ) ) &
          +a4*( fq1  *bzg(lf ,lr1,ls11) +fq  *bzg(lf1,lr1,ls21) )
    bfi =  a1*( fq1  *bfg(lf ,lr ,ls1 ) +fq  *bfg(lf1,lr ,ls2 ) ) &
          +a2*( fq1  *bfg(lf ,lr ,ls11) +fq  *bfg(lf1,lr ,ls21) ) &
          +a3*( fq1  *bfg(lf ,lr1,ls1 ) +fq  *bfg(lf1,lr1,ls2 ) ) &
          +a4*( fq1  *bfg(lf ,lr1,ls11) +fq  *bfg(lf1,lr1,ls21) )
    !*************************************************
    bxt(i) = -bfi*sinphi +bri*cosphi
    byt(i) = +bfi*cosphi +bri*sinphi
    bzt(i) =  bsi                 
  endif
enddo  ! i=1,6
!*************************************************

!
! Calculate interpolated bvector
!
bx=0._rknd
by=0._rknd
bz=0._rknd
do i=1,6
  bx = bx+bxt(i)
  by = by+byt(i)
  bz = bz+bzt(i)
enddo
bx=bx/6._rknd
by=by/6._rknd
bz=bz/6._rknd

! Calculate bval

!*************************************************
if ( bx .eq. 0 .AND. by .eq. 0 .AND. bz .eq. 0 )then
  idiv = 1
!  write(*,*) 'idiv set 2 --> bmod == 0'
else
!*************************************************
  bval(1)  = bx*cosphi0+by*sinphi0    ! br
  bval(2)  =-bx*sinphi0+by*cosphi0    ! bphi
  bval(3)  = bz
endif

return
end subroutine bint


end module bfield_xdr

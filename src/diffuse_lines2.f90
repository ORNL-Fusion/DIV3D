!-----------------------------------------------------------------------------
!+ Main subroutine for following fieldlines and calculating intersections
!-----------------------------------------------------------------------------
Subroutine diffuse_lines2(fname_launch,dmag,dphi_line,nsteps_line,fname_hit, &
period,fname_intpts,fname_nhit,nhitline,lsfi_tol)
!
! Description: 
!
! Inputs: 
!
! Outputs:
!   none
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/15/2011   JDL
!  1.1     09/28/2012   Switched to buffered calls for load balancing.
! Author(s): J.D. Lore - 07/15/2011 - xxx

! Modules used:
Use kind_mod
Use parallel_mod
Use inside_vessel_mod, Only: &
inside_vessel
Use io_unit_spec, Only: &
iu_hit, iu_launch, iu_nhit, iu_vhit, iu_int
Use read_parts_mod
Use integrator_routines_mod
Implicit none

! Input/output
Character(len=100),Intent(in) :: fname_launch, fname_hit, fname_intpts, fname_nhit
Real(rknd), Intent(in) :: dmag, dphi_line, period,lsfi_tol
Integer(iknd), Intent(in) :: nsteps_line, nhitline

! Local scalars
Integer(iknd) :: numl, iline, ii, hitcount, &
  ihit, iocheck, dest, source
Real(rknd) :: Rstart, Zstart, Phistart
Integer(iknd) :: numl_per_proc, numl_extra
Integer(iknd) :: ierr_follow
Integer(iknd) :: WM_OFFSET

! Local arrays
Real(rknd), Dimension(3) :: pint 
Integer(iknd), Dimension(4) :: iout
Real(rknd), Allocatable :: R0(:),Z0(:),Phi0(:)
Real(rknd), Dimension(nhitline) :: r_hitline,z_hitline,phi_hitline
Real(rknd), Dimension(6) :: line_start_data_r
Integer(iknd), Dimension(3) :: line_start_data_i


!- End of header -------------------------------------------------------------

! Read init point data
Write(6,'(A,A)') ' Reading launch point data from ',Trim(Adjustl(fname_launch))
Open(iu_launch,file=fname_launch)
Read(iu_launch,*) numl
Allocate(R0(numl),Z0(numl),Phi0(numl))
Do ii = 1,numl
  Read(iu_launch,*) R0(ii),Z0(ii),Phi0(ii)
Enddo
Close(iu_launch)

If (working_master) Then
  nprocs_working = nprocs
  WM_OFFSET = 0
Else
  nprocs_working = nprocs - 1
  WM_OFFSET = -1
Endif

numl_extra = Mod(numl,nprocs_working)
numl_per_proc = numl/nprocs_working

Write(6,*) 'Total number of fieldlines to follow:',numl
Write(6,*) 'Each working proc gets',numl_per_proc,'plus master gets extra',numl_extra
Write(6,*) 'Diffusing fieldlines with D_mag (m**2/m) = ',dmag,nparts
Write(6,*) 

! Tell each proc how many fls to wait for
Do dest = 1,nprocs - 1
!  Write(6,*) 'master sending numl',numl_per_proc, ' to rank',dest
  call MPI_SEND(numl_per_proc,1,MPI_INTEGER,dest,tag1,MPI_COMM_WORLD,ierr_mpi)
Enddo

hitcount = 0
Open(iu_int,file=fname_intpts,iostat=iocheck)
Open(iu_hit,file=fname_hit,iostat=iocheck)   

Do iline = 1, numl - numl_extra,nprocs_working

  ! Send lines to be followed
  Do dest = 1,nprocs - 1
    Rstart   = R0(iline+dest + WM_OFFSET)
    Zstart   = Z0(iline+dest + WM_OFFSET)
    Phistart = Phi0(iline+dest + WM_OFFSET)
    ! send line number
    call MPI_SEND(iline+dest + WM_OFFSET,1,MPI_INTEGER,dest,tag1,MPI_COMM_WORLD,ierr_mpi)
    ! send [R,Z,P] start
    line_start_data_r(1) = Rstart
    line_start_data_r(2) = Zstart
    line_start_data_r(3) = Phistart
    line_start_data_r(4) = dphi_line
    line_start_data_r(5) = dmag
    line_start_data_r(6) = period

    line_start_data_i(1) = nsteps_line
    line_start_data_i(2) = nhitline
    line_start_data_i(3) = iline+dest+WM_OFFSET

    call MPI_SEND(line_start_data_r,6,MPI_DOUBLE_PRECISION,dest,tag1,MPI_COMM_WORLD,ierr_mpi)
    call MPI_SEND(line_start_data_i,3,MPI_INTEGER,dest,tag1,MPI_COMM_WORLD,ierr_mpi)

!!!!!!!!1    Write(6,*) 'Sent fl ',iline+dest,'to proc ',dest
  Enddo

  If (working_master) Then
    ! Master follows its line
    !!!!!!!!!  Write(6,*) 'Line ',iline,' handled by master'
    Rstart   = R0(iline)
    Zstart   = Z0(iline)
    Phistart = Phi0(iline)
    Call line_follow_and_int(Rstart,Zstart,Phistart,dphi_line,nsteps_line,dmag,period,pint,iout, &
         r_hitline,z_hitline,phi_hitline,nhitline,iline,lsfi_tol)
    ihit = iout(1)
    if ( ihit .ge. 1 ) Then
      hitcount = hitcount + 1
      ! Write to file
      Write(iu_int,*) sqrt(pint(1)*pint(1)+pint(2)*pint(2)),pint(3),atan2(pint(2),pint(1)), iout
      Write(iu_hit,*) nhitline
      Write(iu_hit,*) r_hitline
      Write(iu_hit,*) z_hitline
      Write(iu_hit,*) phi_hitline
    Endif
  Endif !working master

  ! Wait to receive results
  Do source = 1, nprocs -1
    call MPI_RECV(ierr_follow,1,MPI_INTEGER,source,tag1,MPI_COMM_WORLD,status,ierr_mpi)
    call MPI_RECV(iout,4,MPI_INTEGER,source,tag1,MPI_COMM_WORLD,status,ierr_mpi)
    call MPI_RECV(pint,3,MPI_DOUBLE_PRECISION,source,tag1,MPI_COMM_WORLD,status,ierr_mpi)
    call MPI_RECV(r_hitline,nhitline,MPI_DOUBLE_PRECISION,source,tag1,MPI_COMM_WORLD,status,ierr_mpi)
    call MPI_RECV(z_hitline,nhitline,MPI_DOUBLE_PRECISION,source,tag1,MPI_COMM_WORLD,status,ierr_mpi)
    call MPI_RECV(phi_hitline,nhitline,MPI_DOUBLE_PRECISION,source,tag1,MPI_COMM_WORLD,status,ierr_mpi)

!    Write(6,*) 'Master received fl ',iline+source,' from ',source
    ihit = iout(1)
    if ( ihit .ge. 1 ) Then 
      hitcount = hitcount + 1
      Write(iu_int,*)   sqrt(pint(1)*pint(1)+pint(2)*pint(2)),pint(3),atan2(pint(2),pint(1)), iout
      Write(iu_hit,*) nhitline
      Write(iu_hit,*) r_hitline
      Write(iu_hit,*) z_hitline
      Write(iu_hit,*) phi_hitline
    Endif
  Enddo

Enddo ! Line index


! follow extra lines
Do iline = numl-numl_extra+1,numl
  ! Master follows its line
!!!!!!  Write(6,*) 'Master following (extra) fl ',iline
  Rstart   = R0(iline)
  Zstart   = Z0(iline)
  Phistart = Phi0(iline)
  Call line_follow_and_int(Rstart,Zstart,Phistart,dphi_line,nsteps_line,dmag,period,pint,iout, &
       r_hitline,z_hitline,phi_hitline,nhitline,iline,lsfi_tol)

  ihit = iout(1)
  if ( ihit .ge. 1 ) Then 
    hitcount = hitcount + 1
    Write(iu_int,*)   sqrt(pint(1)*pint(1)+pint(2)*pint(2)),pint(3),atan2(pint(2),pint(1)), iout
    Write(iu_hit,*) nhitline
    Write(iu_hit,*) r_hitline
    Write(iu_hit,*) z_hitline
    Write(iu_hit,*) phi_hitline
  Endif

!QQ Write to file

Enddo

Close(iu_hit)
close(iu_int)

Write(6,*) ' We had ',hitcount,' lines -hit-'

Open(iu_nhit,file=fname_nhit)
Write(iu_nhit,*) hitcount
Close(iu_nhit)
Deallocate(R0,Z0,Phi0)

Endsubroutine diffuse_lines2
!-----------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------
Subroutine line_follow_and_int(Rstart,Zstart,Phistart,dphi_line,nsteps_line,dmag,period,pint,iout, &
r_hitline,z_hitline,phi_hitline,nhitline,linenum,lsfi_tol)

Use kind_mod
Use read_parts_mod
Use math_routines_mod, Only: line_seg_facet_int
Implicit None

Real(rknd), Intent(in) :: Rstart, Zstart, Phistart, dmag, dphi_line, period, lsfi_tol
Integer(iknd), Intent(in) :: nsteps_line, linenum !, ntri_max
Integer(iknd), Intent(out), Dimension(4) :: iout
Real(rknd), Dimension(3), Intent(out) :: pint
Integer(iknd), Intent(in) :: nhitline
Real(rknd), Intent(out),Dimension(nhitline) :: r_hitline,z_hitline,phi_hitline


Real(rknd), Dimension(nsteps_line+1) :: rout,zout,phiout
Integer(iknd) :: ifail, imin

!Real(rknd) :: rtol, atol, dphimin
!Integer(iknd) :: nmax_step, method, imin

!- End of header -------------------------------------------------------------




if ( .true.) Then

!
!
!  FIRST FOLLOW LINE THEN CHECK FOR INTERSECTIONS
!
!


!Call follow_fieldline_rzphi_interrupt(Rstart,Zstart,Phistart,dphi_line,nsteps_line,rout,zout,phiout,.true.,dmag,&
!ifail,period,imin,lsfi_tol,linenum)


!write(*,*) ' im here!!'
!stop

Call follow_fieldline_rzphi(Rstart,Zstart,Phistart,dphi_line,nsteps_line,rout,zout,phiout,.true.,dmag,ifail)

!write(*,*) 'rout',rout(1:5)
!  rtol = 1.e-3_rknd ; atol = 1.e-6_rknd ; nmax_step = 100000 ; dphimin = 1.e-9_rknd
!  method = 2
!  Call follow_fieldline_rzphi_ci(Rstart,Zstart,Phistart,dphi_line,nsteps_line,method,rtol,atol, &
!       dphimin,nmax_step,rout,zout,phiout,ifail,.true.,dmag)
!write(*,*) '2out',rout(1:5)
!stop


! Now check for intersections
!Write(6,*) 'Checking line for intersections'
!Write(6,*) 'This line left the grid at point ',ifail

Call check_line_for_intersections(period,pint,iout, &
r_hitline,z_hitline,phi_hitline,nhitline,linenum,lsfi_tol,nsteps_line,rout,zout,phiout,ifail)


else
!
!
!  CHECK FOR INTERSECTIONS AT EACH FIELDLINE FOLLOWING STEP
!
!

Call follow_fieldline_rzphi_and_check(Rstart,Zstart,Phistart,dphi_line,nsteps_line,rout,zout,phiout,.true.,dmag,&
ifail,period,imin,lsfi_tol,linenum,pint,iout)

r_hitline = 0._rknd
phi_hitline = 0._rknd
z_hitline = 0._rknd


endif


End Subroutine line_follow_and_int
!-----------------------------------------------------------------------------






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



  

!  ! Check for nearby triangles
!  Call find_nearby_triangles_v2(r(ii),phi(ii),z(ii),dL,period,dmin,imin_tmp,iclose,r(ii-1),phi(ii-1),z(ii-1),&
!       lsfi_tol,ipart_close,linenum)

!  if (ic_near .ne. 0) Then
!    write(*,*) 'this parts',near_part,near_tri
!    stop
!  endif

  ! Check each triangle identified for an intersection
  Call check_line_segment_for_intersections(period,pint,iout, &
       linenum,lsfi_tol,1_iknd,r(ii-1:ii),z(ii-1:ii),phi(ii-1:ii),0_iknd)
!       r_hitline,z_hitline,phi_hitline,nhitline, &


  deallocate(near_part,near_tri)
  ! If intersection found, quit
!  write(*,*) 'iout!',iout(1)
!  write(*,*) 'r',r(ii-1:ii)
!  write(*,*) 'p',phi(ii-1:ii)
!  write(*,*) 'z',z(ii-1:ii)
 
  If (iout(1) .ne. 0) Then
!    Write(*,*) 'an intersection was found at ii = ',ii
    iout(4) = ii - 1
    ! SET UP HITLINE HERE!!

!    Do jj = ii+1,nsteps+1
!      Call Random_number(rnum)
!    Enddo


    Exit
  endif

Enddo


If (idiv .ne. 0 ) ifail = ii - 1

EndSubroutine follow_fieldline_rzphi_and_check




Subroutine find_nearby_triangles_v2(rr,pp,zz,dL,period)

Use kind_mod
Use read_parts_mod
Use math_routines_mod, Only: line_seg_facet_int
Implicit None
Real(rknd),Intent(in) :: rr,pp,zz,dL, period
Real(rknd) :: xm,ym,zm,R2,P2,Z2,X2,Y2,dd
Integer(iknd) :: ipart, itri, icount, itmp,linenum
Real(rknd), allocatable :: near_part_tmp(:), near_tri_tmp(:)
!Real(rknd), Intent(out) :: dmin
!Integer(iknd), Intent(out) :: imin, iclose, ipart_close
Integer(iknd) :: ii


itmp = sum(ntri_parts(1:nparts))
allocate(near_part_tmp(itmp))
allocate(near_tri_tmp(itmp))

R2 = rr
P2 = pp
Z2 = zz
Do While (pp .lt. 0.d0)
  P2 = P2 + period
Enddo
P2 = mod(P2,period)
X2 = R2*cos(P2)
Y2 = R2*sin(P2)

icount = 1
ic_near = 0
Do ipart=1,nparts  
  Do itri = 1,ntri_parts(ipart)
    xm = xmid(ipart,itri)
    ym = ymid(ipart,itri)
    zm = zmid(ipart,itri)
    
    dd = sqrt( (X2-xm)*(X2-xm) + (Y2-ym)*(Y2-ym) + (Z2-zm)*(Z2-zm) )
    
!    write(*,*) 'this:',dd,dL
!    stop

    if (dd .le. dL) Then    
      ic_near = ic_near + 1
      near_part_tmp(ic_near) = ipart    
      near_tri_tmp(ic_near) = itri
    endif
    
    icount = icount + 1    
  Enddo ! itri
Enddo ! ipart


allocate(near_part(ic_near),near_tri(ic_near))
near_part = near_part_tmp(1:ic_near)
near_tri = near_tri_tmp(1:ic_near)

deallocate(near_part_tmp,near_tri_tmp)

!write(*,*) 'hjeeree>>>>',ic_near
!stop



End Subroutine find_nearby_triangles_v2


!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------

Subroutine check_line_segment_for_intersections(period,pint,iout, &
linenum,lsfi_tol,nsteps_line,rout,zout,phiout,ifail)

!r_hitline,z_hitline,phi_hitline,nhitline,& 

Use kind_mod
Use read_parts_mod
Use inside_vessel_mod, Only: &
inside_vessel
Use math_routines_mod, Only: line_seg_facet_int
Implicit None

Real(rknd), Intent(in) :: period, lsfi_tol
Integer(iknd), Intent(in) :: linenum, nsteps_line,ifail
Integer(iknd), Intent(out), Dimension(4) :: iout
Real(rknd), Dimension(3), Intent(out) :: pint
Real(rknd), Dimension(nsteps_line+1),intent(in) :: rout,zout,phiout

!Integer(iknd), Intent(in) :: nhitline
!Real(rknd), Intent(out),Dimension(nhitline) :: r_hitline,z_hitline,phi_hitline

Integer(iknd) :: iseg, nparts_tmp, ipart_ind, itri_ind
Integer(iknd) :: npts_line, ihit, i, twofer, inphi, ntri, ihit_tmp, ipart, itri, inside_it
Integer(iknd), Dimension(1) :: ind_min, ind_max
Real(rknd) :: R1, Z1, P1, P1a, P2a, X3, Y3, Z3, R3, mu, Aplane, Bplane, denom
Real(rknd) :: p_start, x_start, y_start, z_start, p_end, x_end, y_end, z_end
Real(Rknd), Dimension(2) :: Rtmp, Ztmp, Ytmp, Xtmp
Real(rknd) :: R2,Z2,P2, X1, Y1, X2, Y2, dphi1, dphi2, Pmin, Pmax, X1a, Y1a, Z1a, X2a, Y2a, Z2a


Real(rknd), Dimension(3) :: Pt1, pt2, pa, pb, pc

Real(rknd) :: rtol, atol, dphimin
Integer(iknd) :: nmax_step, method, imin


logical :: close_part_check = .false.

!- End of header -------------------------------------------------------------

!r_hitline = 0.d0
!z_hitline = 0.d0
!phi_hitline = 0.d0


if (close_part_check) Then
  nparts_tmp = ic_near
Else
  nparts_tmp = nparts
Endif

npts_line = nsteps_line + 1
If (ifail .ne. 0) npts_line = ifail - 1

ihit = 0
pint(:) = 0.
iout(:) = -1
iout(1) = ihit
Do i=1,npts_line - 1

  ! Current point along line
  R1 = rout(i)
  Z1 = zout(i)
  P1 = phiout(i)

  ! next point
  R2 = rout(i+1)
  Z2 = zout(i+1)
  P2 = phiout(i+1)
  
  ! Convert to Cartesian coordinates
  dphi1 = P2-P1

  Do While (P1 .lt. 0.d0)
    P1 = P1 + period
  Enddo
  P1 = Mod(P1,period)
  X1 = R1*cos(P1)
  Y1 = R1*sin(P1)

  Do While (P2 .lt. 0.d0)
    P2 = P2 + period
  Enddo
  P2 = Mod(P2,period)
  X2 = R2*cos(P2)
  Y2 = R2*sin(P2)

  dphi2 = P2-P1

  ! Check for lines that cross symmetry (field period) plane
  twofer = 0
  If ( abs(dphi1-dphi2) .gt. 1.d-12) Then 
    twofer = 1
    P1a = Minval([P1,P2])
    P2a = Maxval([P1,P2])
    ind_min = Minloc([P1,P2])
    ind_max = Maxloc([P1,P2])
    Rtmp = [R1,R2]
    Ztmp = [Z1,Z2]
    Xtmp = [X1,X2]
    Ytmp = [Y1,Y2]

    ! Define adjusted two original points
    X1a = Xtmp(ind_min(1))
    Y1a = Ytmp(ind_min(1))
    Z1a = Ztmp(ind_min(1))
    X2a = Xtmp(ind_max(1))
    Y2a = Ytmp(ind_max(1))
    Z2a = Ztmp(ind_max(1))


    !---------------------------------------------
    ! Find intersection of line segment with plane
    !---------------------------------------------
    pt1(1) = Rtmp(ind_min(1))*cos(P1a+period)
    pt1(2) = Rtmp(ind_min(1))*sin(P1a+period)
    pt1(3) = Ztmp(ind_min(1))
    pt2(1) = Xtmp(ind_max(1))
    pt2(2) = Ytmp(ind_max(1))
    pt2(3) = Ztmp(ind_max(1))

    ! Define explicit form of plane at 2*pi/Nfp 
    pa = [0.d0,0.d0,0.d0]
    pb = [0.d0,0.d0,1.d0]
    pc = [cos(period),sin(period),0.d0]
    
    ! Calc unit vector normal to plane of Pa-c (gives plane components A-C)
    Aplane=-sin(period)
    Bplane= cos(period)

    ! Calculate the position on the line that intersects the plane
    denom = Aplane*(pt2(1) - pt1(1)) + Bplane*(pt2(2) - pt1(2))
    If (abs(denom) .lt. 1.d-15) Then
      Write(6,*) 'Error: Line did not hit plane at 2*pi/Nfp'
      Stop
    Endif
    mu = - (Aplane*pt1(1) + Bplane*pt1(2)) / denom
    X3 = pt1(1) + mu * (pt2(1) - pt1(1))
    Y3 = pt1(2) + mu * (pt2(2) - pt1(2))
    Z3 = pt1(3) + mu * (pt2(3) - pt1(3))
    R3 = sqrt(X3*X3+Y3*Y3)
    
    ! Define two new points 
    P1 = 0._rknd
    X1 = R3*cos(P1)
    Y1 = R3*sin(P1)
    Z1 = Z3
    P2 = period
    X2 = R3*cos(P2)
    Y2 = R3*sin(P2)
    Z2 = Z3
  endif

  Do iseg = 1,1+twofer

  if (twofer .eq. 0) Then
    p_start = P1
    x_start = X1
    y_start = Y1
    z_start = Z1
    p_end   = P2
    x_end   = X2
    y_end   = Y2
    z_end   = Z2
  endif
  if (twofer .eq. 1) Then
    if ( iseg .eq. 1) Then
        p_start = P1
        x_start = X1
        y_start = Y1
        z_start = Z1
        p_end   = P1a
        x_end   = X1a
        y_end   = Y1a
        z_end   = Z1a
    else
        p_start = P2a
        x_start = X2a
        y_start = Y2a
        z_start = Z2a
        p_end   = P2
        x_end   = X2
        y_end   = Y2
        z_end   = Z2
    endif
  endif

    
  ihit = 0
  Do ipart_ind=1,nparts_tmp

    if (close_part_check) Then
      ipart = near_part(ipart_ind)
    Else
      ipart = ipart_ind
    Endif

    !Check if this line segment is within the phi bounds of the part
    Pmin = Pmins(ipart)
    Pmax = Pmaxs(ipart)

    inphi = 0

    If ( p_start .ge. Pmin .AND. p_start .le. Pmax ) inphi = 1
    If ( p_end .ge. Pmin .AND. p_end .le. Pmax ) inphi = 1

    If (inphi .eq. 1 ) Then

      if (close_part_check) Then
        ntri = ic_near
      Else          
        ntri = ntri_parts(ipart_ind)
      Endif
      Do itri_ind = 1,ntri
          
      if (close_part_check) Then
        itri = near_tri(itri_ind)
      Else          
        itri = itri_ind
      Endif


!        If (check_tri(ipart,itri) .eq. 1 ) Then

          pa = [xtri(ipart,itri,1),ytri(ipart,itri,1),ztri(ipart,itri,1)]
          pb = [xtri(ipart,itri,2),ytri(ipart,itri,2),ztri(ipart,itri,2)]
          pc = [xtri(ipart,itri,3),ytri(ipart,itri,3),ztri(ipart,itri,3)]
          pt1 = [x_start,y_start,z_start]
          pt2 = [x_end,y_end,z_end]

          Call line_seg_facet_int(pa,pb,pc,pt1,pt2,ihit_tmp,pint,lsfi_tol)          

          If (ihit_tmp .eq. 1) Then
            ihit = 1
            Write(6,'(A,I0,A,I0,A,I0,3(F8.3))') ' Line ',linenum,' hit! [i,ipart,P] ',i,' ',ipart,pint(1),pint(2),pint(3)
            iout(1) = ihit
            iout(2) = ipart
            iout(3) = itri
            iout(4) = i
!            if ( (i - nhitline+1) .lt. 1 ) Then
!              Write(*,*) 'Truncating hitline'
!              r_hitline = 0.d0
!              z_hitline = 0.d0
!              phi_hitline = 0.d0
!              r_hitline(1:i) = rout(1:i)
!              z_hitline(1:i) = zout(1:i)
!              phi_hitline(1:i) = phiout(1:i)
!            Else
!              r_hitline = rout(i-nhitline+1:i)
!              z_hitline = zout(i-nhitline+1:i)
!              phi_hitline = phiout(i-nhitline+1:i)
!            Endif
            Exit ! stop looking for triangle intersections
          Endif

!        Else ! check tri
!          Write(6,*) 'Should not be here unless 3d parts have been implemented'
!          stop
!        Endif
          
      Enddo !triangle loop

      If (ihit .eq. 1 ) Exit  ! Stop looking for part intersections
    Endif !inphi check
  Enddo ! part index

      If (ihit .eq. 1 ) Exit  ! Stop looking for segment intersections
  Enddo ! seg index (twofer)

  If (ihit .eq. 1 ) Exit  ! Quit fieldline if it hit a part
Enddo ! points along line (i)


if (.false.) Then
If ( ihit .eq. 0 ) Then
  Write(6,'(A,I0,A)') ' No part intersections found for line ',linenum,', checking for vessel intersection'
  Do i=1,npts_line
    ! Current point along line
    R1 = rout(i)
    Z1 = zout(i)
    P1 = phiout(i)

    inside_it = inside_vessel(R1,Z1,P1,R_ves,Z_ves,P_ves,ntor_ves,npol_ves,msym_ves)
    If (inside_it .eq. 0 ) Then
      Write(6,'(A,I0,A,I0,3(F8.3))') ' Line ',linenum,' did hit the vessel at [i,R,Z,P]=',i,R1,Z1,P1
      ihit = 2
      iout(1) = ihit
      iout(2) = -2
      iout(3) = -2
      iout(4) = i
      Exit
    Endif
  Enddo
Endif
Endif

End subroutine check_line_segment_for_intersections







!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------

Subroutine check_line_for_intersections(period,pint,iout, &
r_hitline,z_hitline,phi_hitline,nhitline,linenum,lsfi_tol,nsteps_line,rout,zout,phiout,ifail)

Use kind_mod
Use read_parts_mod
Use inside_vessel_mod, Only: &
inside_vessel
Use math_routines_mod, Only: line_seg_facet_int
Implicit None

Real(rknd), Intent(in) :: period, lsfi_tol
Integer(iknd), Intent(in) :: linenum, nsteps_line,ifail
Integer(iknd), Intent(out), Dimension(4) :: iout
Real(rknd), Dimension(3), Intent(out) :: pint
Integer(iknd), Intent(in) :: nhitline
Real(rknd), Intent(out),Dimension(nhitline) :: r_hitline,z_hitline,phi_hitline

Integer(iknd) :: iseg
Integer(iknd) :: npts_line, ihit, i, twofer, inphi, ntri, ihit_tmp, ipart, itri, inside_it
Integer(iknd), Dimension(1) :: ind_min, ind_max
Real(rknd) :: R1, Z1, P1, P1a, P2a, X3, Y3, Z3, R3, mu, Aplane, Bplane, denom
Real(rknd) :: p_start, x_start, y_start, z_start, p_end, x_end, y_end, z_end
Real(Rknd), Dimension(2) :: Rtmp, Ztmp, Ytmp, Xtmp
Real(rknd) :: R2,Z2,P2, X1, Y1, X2, Y2, dphi1, dphi2, Pmin, Pmax, X1a, Y1a, Z1a, X2a, Y2a, Z2a

Real(rknd), Dimension(nsteps_line+1),intent(in) :: rout,zout,phiout
Real(rknd), Dimension(3) :: Pt1, pt2, pa, pb, pc

Real(rknd) :: rtol, atol, dphimin
Integer(iknd) :: nmax_step, method, imin

!- End of header -------------------------------------------------------------

r_hitline = 0.d0
z_hitline = 0.d0
phi_hitline = 0.d0


npts_line = nsteps_line + 1
If (ifail .ne. 0) npts_line = ifail - 1

ihit = 0
pint(:) = 0.
iout(:) = -1
iout(1) = ihit
Do i=1,npts_line - 1
!write(*,*) 'here:>>>>>>>>>>>',imin
!Do i=imin,npts_line - 1
!  Write(*,*) i
  ! Current point along line
  R1 = rout(i)
  Z1 = zout(i)
  P1 = phiout(i)

  ! next point
  R2 = rout(i+1)
  Z2 = zout(i+1)
  P2 = phiout(i+1)
  
  ! Convert to Cartesian coordinates
  dphi1 = P2-P1

  Do While (P1 .lt. 0.d0)
    P1 = P1 + period
  Enddo
  P1 = Mod(P1,period)
  X1 = R1*cos(P1)
  Y1 = R1*sin(P1)

  Do While (P2 .lt. 0.d0)
    P2 = P2 + period
  Enddo
  P2 = Mod(P2,period)
  X2 = R2*cos(P2)
  Y2 = R2*sin(P2)

  dphi2 = P2-P1

  ! Check for lines that cross symmetry (field period) plane
  twofer = 0
  If ( abs(dphi1-dphi2) .gt. 1.d-12) Then 
    twofer = 1
    P1a = Minval([P1,P2])
    P2a = Maxval([P1,P2])
    ind_min = Minloc([P1,P2])
    ind_max = Maxloc([P1,P2])
    Rtmp = [R1,R2]
    Ztmp = [Z1,Z2]
    Xtmp = [X1,X2]
    Ytmp = [Y1,Y2]

    ! Define adjusted two original points
    X1a = Xtmp(ind_min(1))
    Y1a = Ytmp(ind_min(1))
    Z1a = Ztmp(ind_min(1))
    X2a = Xtmp(ind_max(1))
    Y2a = Ytmp(ind_max(1))
    Z2a = Ztmp(ind_max(1))


    !---------------------------------------------
    ! Find intersection of line segment with plane
    !---------------------------------------------
    pt1(1) = Rtmp(ind_min(1))*cos(P1a+period)
    pt1(2) = Rtmp(ind_min(1))*sin(P1a+period)
    pt1(3) = Ztmp(ind_min(1))
    pt2(1) = Xtmp(ind_max(1))
    pt2(2) = Ytmp(ind_max(1))
    pt2(3) = Ztmp(ind_max(1))

    ! Define explicit form of plane at 2*pi/Nfp 
    pa = [0.d0,0.d0,0.d0]
    pb = [0.d0,0.d0,1.d0]
    pc = [cos(period),sin(period),0.d0]
    
    ! Calc unit vector normal to plane of Pa-c (gives plane components A-C)
    Aplane=-sin(period)
    Bplane= cos(period)

    ! Calculate the position on the line that intersects the plane
    denom = Aplane*(pt2(1) - pt1(1)) + Bplane*(pt2(2) - pt1(2))
    If (abs(denom) .lt. 1.d-15) Then
      Write(6,*) 'Error: Line did not hit plane at 2*pi/Nfp'
      Stop
    Endif
    mu = - (Aplane*pt1(1) + Bplane*pt1(2)) / denom
    X3 = pt1(1) + mu * (pt2(1) - pt1(1))
    Y3 = pt1(2) + mu * (pt2(2) - pt1(2))
    Z3 = pt1(3) + mu * (pt2(3) - pt1(3))
    R3 = sqrt(X3*X3+Y3*Y3)
    
    ! Define two new points 
    P1 = 0._rknd
    X1 = R3*cos(P1)
    Y1 = R3*sin(P1)
    Z1 = Z3
    P2 = period
    X2 = R3*cos(P2)
    Y2 = R3*sin(P2)
    Z2 = Z3
  endif

  Do iseg = 1,1+twofer

  if (twofer .eq. 0) Then
    p_start = P1
    x_start = X1
    y_start = Y1
    z_start = Z1
    p_end   = P2
    x_end   = X2
    y_end   = Y2
    z_end   = Z2
  endif
  if (twofer .eq. 1) Then
    if ( iseg .eq. 1) Then
        p_start = P1
        x_start = X1
        y_start = Y1
        z_start = Z1
        p_end   = P1a
        x_end   = X1a
        y_end   = Y1a
        z_end   = Z1a
    else
        p_start = P2a
        x_start = X2a
        y_start = Y2a
        z_start = Z2a
        p_end   = P2
        x_end   = X2
        y_end   = Y2
        z_end   = Z2
    endif
  endif

    
  ihit = 0
  Do ipart=1,nparts

    !Check if this line segment is within the phi bounds of the part
    Pmin = Pmins(ipart)
    Pmax = Pmaxs(ipart)

    inphi = 0
!    If ( twofer .eq. 1 ) Then  ! This line crossed the f.p. plane (2pi/Nfp)
!      If ( P1a .ge. Pmin .AND. P1a .le. Pmax ) inphi = 1
!      If ( P2a .ge. Pmin .AND. P2a .le. Pmax ) inphi = 1
!      If (inphi .eq. 1 ) Then
!        Write(6,*) 'NEED TO WRITE THIS PART!!!!'
!        WRite(6,*) 'ipart',ipart
!        Write(6,*) P1a,P2a,Pmin,Pmax,dphi1,dphi2,abs(dphi1-dphi2)
!        stop
!      Endif
!    Endif

    If ( p_start .ge. Pmin .AND. p_start .le. Pmax ) inphi = 1
    If ( p_end .ge. Pmin .AND. p_end .le. Pmax ) inphi = 1

    If (inphi .eq. 1 ) Then
      ntri = ntri_parts(ipart)
      Do itri = 1,ntri
          
        If (check_tri(ipart,itri) .eq. 1 ) Then

          pa = [xtri(ipart,itri,1),ytri(ipart,itri,1),ztri(ipart,itri,1)]
          pb = [xtri(ipart,itri,2),ytri(ipart,itri,2),ztri(ipart,itri,2)]
          pc = [xtri(ipart,itri,3),ytri(ipart,itri,3),ztri(ipart,itri,3)]
          pt1 = [x_start,y_start,z_start]
          pt2 = [x_end,y_end,z_end]

          Call line_seg_facet_int(pa,pb,pc,pt1,pt2,ihit_tmp,pint,lsfi_tol)          

!          If ( ihit .eq. 1 ) Then
!            Write(6,*) 'Something is wrong'
!          Endif

          If (ihit_tmp .eq. 1) Then
            ihit = 1
            Write(6,'(A,I0,A,I0,A,I0,3(F8.3))') ' Line ',linenum,' hit! [i,ipart,P] ',i,' ',ipart,pint(1),pint(2),pint(3)
            iout(1) = ihit
            iout(2) = ipart
            iout(3) = itri
            iout(4) = i
            if ( (i - nhitline+1) .lt. 1 ) Then
              !Write(6,*) 'Write something to handle this', i, nhitline
              !Stop
              Write(*,*) 'Truncating hitline'
              r_hitline = 0.d0
              z_hitline = 0.d0
              phi_hitline = 0.d0
              r_hitline(1:i) = rout(1:i)
              z_hitline(1:i) = zout(1:i)
              phi_hitline(1:i) = phiout(1:i)
            Else
              r_hitline = rout(i-nhitline+1:i)
              z_hitline = zout(i-nhitline+1:i)
              phi_hitline = phiout(i-nhitline+1:i)
            Endif
            Exit ! stop looking for triangle intersections
          Endif

        Else ! check tri
          Write(6,*) 'Should not be here unless 3d parts have been implemented'
          stop
        Endif
          
      Enddo !triangle loop

      If (ihit .eq. 1 ) Exit  ! Stop looking for part intersections
    Endif !inphi check
  Enddo ! part index

      If (ihit .eq. 1 ) Exit  ! Stop looking for segment intersections
  Enddo ! seg index (twofer)

  If (ihit .eq. 1 ) Exit  ! Quit fieldline if it hit a part
Enddo ! points along line (i)



If ( ihit .eq. 0 ) Then
  Write(6,'(A,I0,A)') ' No part intersections found for line ',linenum,', checking for vessel intersection'
  Do i=1,npts_line
    ! Current point along line
    R1 = rout(i)
    Z1 = zout(i)
    P1 = phiout(i)

    inside_it = inside_vessel(R1,Z1,P1,R_ves,Z_ves,P_ves,ntor_ves,npol_ves,msym_ves)
    If (inside_it .eq. 0 ) Then
      Write(6,'(A,I0,A,I0,3(F8.3))') ' Line ',linenum,' did hit the vessel at [i,R,Z,P]=',i,R1,Z1,P1
      ihit = 2
      iout(1) = ihit
      iout(2) = -2
      iout(3) = -2
      iout(4) = i
      Exit
    Endif
  Enddo
Endif


End subroutine check_line_for_intersections

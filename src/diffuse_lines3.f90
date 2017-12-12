!-----------------------------------------------------------------------------
!+ Main subroutine for following fieldlines and calculating intersections
!-----------------------------------------------------------------------------
Subroutine diffuse_lines3(fname_launch,dmag,nsteps_line,fname_hit, &
fname_intpts,fname_nhit,nhitline)
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
iu_hit, iu_launch, iu_nhit, iu_int
Use read_parts_mod
Implicit none

! Input/output
Character(len=100),Intent(in) :: fname_launch, fname_hit, fname_intpts, fname_nhit
Real(real64), Intent(in) :: dmag
Integer(int32), Intent(in) :: nsteps_line, nhitline

! Local scalars
Integer(int32) :: numl, iline, ii, hitcount, &
  ihit, iocheck, dest
Real(real64) :: Rstart, Zstart, Phistart

Integer(int32) :: tag, flag, work_done, work_done_count

! Local arrays
Real(real64), Dimension(3) :: pint 
Integer(int32), Dimension(4) :: iout
Real(real64), Allocatable :: R0(:),Z0(:),Phi0(:)
Real(real64), Dimension(nhitline) :: r_hitline,z_hitline,phi_hitline
Real(real64), Dimension(3) :: line_start_data_r
Integer(int32), Dimension(3) :: line_start_data_i
Real(real64), Dimension(3) :: line_done_data_r
Real(real64), Dimension(3*nhitline) :: line_done_data_r2
Integer(int32), Dimension(5) :: line_done_data_i

Integer(int32), Dimension(:), Allocatable :: this_job_done, is_working_arr, req_arr

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

Write(6,*) 'Total number of fieldlines to follow:',numl
Write(6,*) 'Diffusing fieldlines with D_mag (m**2/m) = ',dmag
Write(6,*) 

Allocate(this_job_done(nprocs-1),is_working_arr(nprocs-1),req_arr(nprocs-1))

!-----------------------------------------------------------------------------
! Send initial fls and request results
!-----------------------------------------------------------------------------
iline = 0
work_done_count = 0
req_arr = 0
is_working_arr = 0
Do dest = 1,nprocs - 1
   
  tag = dest

  ! Send line start data
  iline = iline + 1
  Rstart   = R0(iline)
  Zstart   = Z0(iline)
  Phistart = Phi0(iline)
  
  line_start_data_i(1) = nsteps_line
  line_start_data_i(2) = nhitline
  line_start_data_i(3) = iline
  Call MPI_SEND(line_start_data_i,3,MPI_INTEGER         ,dest,tag,MPI_COMM_WORLD,ierr_mpi)

  line_start_data_r(1) = Rstart
  line_start_data_r(2) = Zstart
  line_start_data_r(3) = Phistart  
  Call MPI_SEND(line_start_data_r,3,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,ierr_mpi)

  ! Request signal that job is finished
  Call MPI_IRECV(this_job_done(dest),1,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,request,ierr_mpi)
  req_arr(dest) = request
  is_working_arr(dest) = 1
  
Enddo

!-----------------------------------------------------------
! Loop over processors until work is complete
!-----------------------------------------------------------
work_done = 0
hitcount = 0

Open(iu_int,file=fname_intpts,iostat=iocheck)
Open(iu_hit,file=fname_hit   ,iostat=iocheck)   
Do While ( work_done .ne. 1 )
  Do dest = 1,nprocs - 1

    ! If we are still waiting on a job from this process, check if it is finished
    If (is_working_arr(dest) == 1) Then
      request = req_arr(dest)
      Call MPI_TEST(request,flag,status,ierr_mpi) 
      If (flag == 1) Then
        is_working_arr(dest) = 0
        work_done_count = work_done_count + 1

        ! Request results
        tag = dest
        Call MPI_RECV(line_done_data_i ,5         ,MPI_INTEGER         ,dest,tag,MPI_COMM_WORLD,status,ierr_mpi)
        Call MPI_RECV(line_done_data_r ,3         ,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,status,ierr_mpi)

        If (nhitline .gt. 0) Then
           Call MPI_RECV(line_done_data_r2,nhitline*3,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,status,ierr_mpi)
        Endif
        ihit        = line_done_data_i(2)
        if ( ihit .ge. 1 ) Then 
          hitcount = hitcount + 1
          iout = line_done_data_i(2:5)
          pint = line_done_data_r(1:3)
          Write(iu_int,*)   sqrt(pint(1)*pint(1)+pint(2)*pint(2)),pint(3),atan2(pint(2),pint(1)),iout
          If (nhitline .gt. 0) Then
             r_hitline   = line_done_data_r2(1+0*nhitline:1*nhitline)
             z_hitline   = line_done_data_r2(1+1*nhitline:2*nhitline)
             phi_hitline = line_done_data_r2(1+2*nhitline:3*nhitline)
             Write(iu_hit,*) nhitline
             Write(iu_hit,*) r_hitline
             Write(iu_hit,*) z_hitline
             Write(iu_hit,*) phi_hitline
          Endif
        Endif
      Endif ! flag == 1

    ! If the process is not working, send a new job (if there are still jobs to do)           
    Elseif (is_working_arr(dest) == 0 ) Then        
!      write(*,*) 'process ',dest,' is reporting idle, sending new job'
      If (iline .lt. numl) Then       
        tag = dest
        iline = iline + 1
        Rstart   = R0(iline)
        Zstart   = Z0(iline)
        Phistart = Phi0(iline)
        line_start_data_r(1) = Rstart
        line_start_data_r(2) = Zstart
        line_start_data_r(3) = Phistart
        line_start_data_i(1) = nsteps_line
        line_start_data_i(2) = nhitline
        line_start_data_i(3) = iline
        Call MPI_SEND(line_start_data_i,3,MPI_INTEGER         ,dest,tag,MPI_COMM_WORLD,ierr_mpi)
        Call MPI_SEND(line_start_data_r,3,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,ierr_mpi)
        is_working_arr(dest) = 1
        
        ! Initiate request for completed job
        call MPI_IRECV(this_job_done(dest),1,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,request,ierr_mpi)       
        req_arr(dest) = request
      Else
        ! If there are no more jobs, turn this process 'off'
        is_working_arr(dest) = -1
      Endif
      
    Endif
    
  Enddo ! nprocs

  if ( work_done_count .ge. numl ) Then
    work_done = 1
  Endif
  
Enddo ! while NOT work done

! 
! Send signal to each process that we are done
!

Do dest = 1,nprocs-1
  tag = dest
  Write(*,*) 'Master sending kill signal to process',dest
  line_start_data_i = -1 ! Kill signal
  Call MPI_SEND(line_start_data_i,3,MPI_INTEGER         ,dest,tag,MPI_COMM_WORLD,ierr_mpi) 
Enddo


Close(iu_hit)
close(iu_int)

Write(6,*) ' We had ',hitcount,' lines -hit-'

Open(iu_nhit,file=fname_nhit)
Write(iu_nhit,*) hitcount
Close(iu_nhit)
Deallocate(R0,Z0,Phi0)

End Subroutine diffuse_lines3
!-----------------------------------------------------------------------------



!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------
Subroutine line_follow_and_int(Rstart,Zstart,Phistart,dphi_line,nsteps_line,dmag,period,pint,iout, &
r_hitline,z_hitline,phi_hitline,nhitline,linenum,lsfi_tol)

Use kind_mod
Use read_parts_mod
Use math_routines_mod, Only: line_seg_facet_int
Use fieldline_follow_mod, Only : follow_fieldlines_rzphi_diffuse
Use setup_bfield_module, Only : bfield
Implicit None

Real(real64), Intent(in) :: Rstart, Zstart, Phistart, dmag, dphi_line, period, lsfi_tol
Integer(int32), Intent(in) :: nsteps_line, linenum
Integer(int32), Intent(out), Dimension(4) :: iout
Real(real64), Dimension(3), Intent(out) :: pint
Integer(int32), Intent(in) :: nhitline
Real(real64), Intent(out),Dimension(nhitline) :: r_hitline,z_hitline,phi_hitline


Real(real64), Dimension(nsteps_line+1) :: rout,zout,phiout
Integer(int32) :: ifail, ierr(1),ilg(1)

!- End of header -------------------------------------------------------------


Call follow_fieldlines_rzphi_diffuse(bfield,(/Rstart/),(/Zstart/),(/Phistart/),1,&
     dphi_line,nsteps_line,rout,zout,phiout,ierr,ilg,dmag)
ifail = ilg(1)

Call check_line_for_intersections(period,pint,iout, &
     r_hitline,z_hitline,phi_hitline,nhitline,linenum,lsfi_tol,nsteps_line,rout,zout,phiout,ifail)


End Subroutine line_follow_and_int
!-----------------------------------------------------------------------------




!-----------------------------------------------------------------------------
!+ 
!-----------------------------------------------------------------------------

Subroutine check_line_segment_for_intersections(period,pint,iout, &
linenum,lsfi_tol,nsteps_line,rout,zout,phiout,ifail)

Use kind_mod
Use read_parts_mod
Use inside_vessel_mod, Only: &
inside_vessel
Use math_routines_mod, Only: line_seg_facet_int
Implicit None

Real(real64), Intent(in) :: period, lsfi_tol
Integer(int32), Intent(in) :: linenum, nsteps_line,ifail
Integer(int32), Intent(out), Dimension(4) :: iout
Real(real64), Dimension(3), Intent(out) :: pint
Real(real64), Dimension(nsteps_line+1),intent(in) :: rout,zout,phiout

Integer(int32) :: iseg, nparts_tmp, ipart_ind, itri_ind
Integer(int32) :: npts_line, ihit, i, twofer, inphi, ntri, ihit_tmp, ipart, itri, inside_it
Integer(int32), Dimension(1) :: ind_min, ind_max
Real(real64) :: R1, Z1, P1, P1a, P2a, X3, Y3, Z3, R3, mu, Aplane, Bplane, denom
Real(real64) :: p_start, x_start, y_start, z_start, p_end, x_end, y_end, z_end
Real(Real64), Dimension(2) :: Rtmp, Ztmp, Ytmp, Xtmp
Real(real64) :: R2,Z2,P2, X1, Y1, X2, Y2, dphi1, dphi2, Pmin, Pmax, X1a, Y1a, Z1a, X2a, Y2a, Z2a


Real(real64), Dimension(3) :: Pt1, pt2, pa, pb, pc

logical :: close_part_check = .false.

!- End of header -------------------------------------------------------------


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
    P1 = 0._real64
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
            Exit ! stop looking for triangle intersections
          Endif
          
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

Real(real64), Intent(in) :: period, lsfi_tol
Integer(int32), Intent(in) :: linenum, nsteps_line,ifail
Integer(int32), Intent(out), Dimension(4) :: iout
Real(real64), Dimension(3), Intent(out) :: pint
Integer(int32), Intent(in) :: nhitline
Real(real64), Intent(out),Dimension(nhitline) :: r_hitline,z_hitline,phi_hitline

Integer(int32) :: iseg
Integer(int32) :: npts_line, ihit, i, twofer, inphi, ntri, ihit_tmp, ipart, itri, inside_it
Integer(int32), Dimension(1) :: ind_min, ind_max
Real(real64) :: R1, Z1, P1, P1a, P2a, X3, Y3, Z3, R3, mu, Aplane, Bplane, denom
Real(real64) :: p_start, x_start, y_start, z_start, p_end, x_end, y_end, z_end
Real(Real64), Dimension(2) :: Rtmp, Ztmp, Ytmp, Xtmp
Real(real64) :: R2,Z2,P2, X1, Y1, X2, Y2, dphi1, dphi2, Pmin, Pmax, X1a, Y1a, Z1a, X2a, Y2a, Z2a

Real(real64), Dimension(nsteps_line+1),intent(in) :: rout,zout,phiout
Real(real64), Dimension(3) :: Pt1, pt2, pa, pb, pc

!- End of header -------------------------------------------------------------

r_hitline = 0.d0
z_hitline = 0.d0
phi_hitline = 0.d0


npts_line = nsteps_line + 1
If (ifail .ne. nsteps_line) npts_line = ifail - 1

ihit = 0
pint(:) = 0.
iout(:) = -1
iout(1) = ihit
Do i=1,npts_line - 1
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
    P1 = 0._real64
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

Module diffusion
  Implicit None
  Public :: diffuse_lines3, line_follow_and_int, diffuse_lines3_worker

  ! Message array sizes
  Integer, Private, Parameter :: size_rdata = 7
  Integer, Private, Parameter :: size_idata = 5

Contains

  !-----------------------------------------------------------------------------
  !+ Follow fieldlines and then compute intersections, used by worker nodes
  !-----------------------------------------------------------------------------
  Subroutine line_follow_and_int(Rstart,Zstart,Phistart,pint,iout, &
       r_hitline,z_hitline,phi_hitline,nhitline,linenum,totL,theta, &
       etime_follow,etime_int)

    Use kind_mod, Only : real64, int32
    Use fieldline_follow_mod, Only : follow_fieldlines_rzphi_diffuse
    Use setup_bfield_module, Only : bfield
    Use run_settings_namelist, Only : dmag, dphi_line_diff, calc_lc, calc_theta, ns_line_diff, &
         lambda_par, randomize_start_dir
    Use timing_mod, Only : get_elapsed_time
    Implicit None

    Real(real64), Intent(in) :: Rstart, Zstart, Phistart
    Integer(int32), Intent(in) :: linenum, nhitline
    Integer(int32), Intent(out), Dimension(4) :: iout
    Real(real64), Dimension(3), Intent(out) :: pint
    Real(real64), Intent(out) :: totL, theta
    Real(real64), Intent(out), Dimension(nhitline) :: r_hitline,z_hitline,phi_hitline
    Real(real64), Intent(out) :: etime_follow, etime_int
    Real(real64), Allocatable :: rout(:) ,zout(:) ,phiout(:)
    Integer(int32) :: ifail, ierr(1),ilg(1)
    Integer :: tstart
    !- End of header -------------------------------------------------------------

    Allocate(rout(ns_line_diff+1),zout(ns_line_diff+1),phiout(ns_line_diff+1))
    
    ! Follow
    Call system_clock(tstart)
    Call follow_fieldlines_rzphi_diffuse(bfield,[Rstart],[Zstart],[Phistart],1,&
         dphi_line_diff,ns_line_diff,rout,zout,phiout,ierr,ilg,dmag,lambda_par,randomize_start_dir)
    ifail = ilg(1)
    etime_follow = get_elapsed_time(tstart)

    ! Check intersection
    Call system_clock(tstart)
    Call check_line_for_intersections(pint,iout, &
         r_hitline,z_hitline,phi_hitline,nhitline,linenum, &
         ns_line_diff,rout,zout,phiout,ifail,totL,calc_lc,calc_theta,theta)
    etime_int = get_elapsed_time(tstart)

    Deallocate(rout,zout,phiout)

  End Subroutine line_follow_and_int
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !+ Worker node subroutine for following fieldlines and calculating intersections
  !  Counterpart to diffuse_lines3
  !-----------------------------------------------------------------------------
  Subroutine diffuse_lines3_worker
    Use kind_mod, Only : real64, int32, int16
    Use run_settings_namelist, Only : dmag, nhitline,calc_lc, calc_theta
    Use parallel_mod!, Only : ierr_mpi, rank, status, &
!         MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_PRECISION
    Implicit None

    Real(real64), Dimension(3) :: line_start_data_r
    Integer(int32), Dimension(3) :: line_start_data_i
    Real(real64) :: Rstart_local, Zstart_local, Pstart_local
    Integer(int32) :: iline_local !, nsteps_line_local
    Integer :: dest, source, tag
    Real(real64), Dimension(:), Allocatable :: r_hitline,z_hitline,phi_hitline
    Integer(int32) :: ierr_follow
    Integer(int32) :: num_myjobs
    Real(real64) :: totL, theta, etime_follow, etime_int

    Real(real64), Dimension(size_rdata) :: line_done_data_r
    Integer(int32), Dimension(size_idata) :: line_done_data_i

    Real(real64), Dimension(:), Allocatable :: line_done_data_r2
    Integer(int32), Dimension(4) :: iout
    Real(real64), Dimension(3) :: pint
    Integer(int16) :: buffer

    ! Initialize status
    num_myjobs = 0
    line_start_data_i = 0

    !  Write(*,*) 'Process ',rank,' reporting as READY'
    Do While (line_start_data_i(1) .ne. -1)

       ! Wait for line data
       source = 0
       dest = 0
       tag = rank
       Call MPI_RECV(line_start_data_i,3,MPI_INTEGER,source,tag,MPI_COMM_WORLD,status,ierr_mpi)

       ! Check for kill signal
       if (line_start_data_i(1) .ne. -1 ) Then

          ! Get rest of start data
          Call MPI_RECV(line_start_data_r,3,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,status,ierr_mpi)
          Rstart_local    = line_start_data_r(1)
          Zstart_local    = line_start_data_r(2)
          Pstart_local    = line_start_data_r(3)
!          nsteps_line_local = line_start_data_i(1)
          iline_local       = line_start_data_i(3)

          ! Follow the line
          Allocate(r_hitline(nhitline))
          Allocate(z_hitline(nhitline))
          Allocate(phi_hitline(nhitline))

          Call line_follow_and_int(Rstart_local,Zstart_local,Pstart_local,pint,iout, &
               r_hitline,z_hitline,phi_hitline,nhitline,iline_local, &
               totL,theta,etime_follow,etime_int)
          ierr_follow = 0

          ! Compile results and send data back to master

          ! first send handshake signal (this_job_done)
          buffer = 1
          Call MPI_SEND(buffer,1,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,ierr_mpi)

          line_done_data_i(1) = ierr_follow
          line_done_data_i(2:5) = iout
          Call MPI_SEND(line_done_data_i,size_idata,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,ierr_mpi)

          line_done_data_r(1:3) = pint
          line_done_data_r(4) = totL
          line_done_data_r(5) = theta
          line_done_data_r(6) = etime_follow
          line_done_data_r(7) = etime_int
          Call MPI_SEND(line_done_data_r,size_rdata,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,ierr_mpi)

          If (nhitline .gt. 0) Then
             Allocate(line_done_data_r2(3*nhitline))
             line_done_data_r2(1+0*nhitline:1*nhitline) = r_hitline
             line_done_data_r2(1+1*nhitline:2*nhitline) = z_hitline
             line_done_data_r2(1+2*nhitline:3*nhitline) = phi_hitline
             Call MPI_SEND(line_done_data_r2,nhitline*3,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,ierr_mpi)
             Deallocate(line_done_data_r2)
          Endif
          Deallocate(r_hitline,z_hitline,phi_hitline)

          num_myjobs = num_myjobs + 1
       Endif ! kill signal check

    EndDo ! while mywork ne -1

    Write(*,*) 'Process ',rank,'received signal that all work is complete'
    Write(*,*) 'Process ',rank,' completed ',num_myjobs,' jobs'
  End Subroutine diffuse_lines3_worker

  !-----------------------------------------------------------------------------
  !+ Main subroutine for following fieldlines and calculating intersections
  !-----------------------------------------------------------------------------
  Subroutine diffuse_lines3
    ! Author(s): J.D. Lore - 07/15/2011 - xxx

    ! Modules used:
    Use kind_mod, Only : int32, real64
    Use parallel_mod
    Use io_unit_spec, Only: iu_hit, iu_launch, iu_nhit, iu_int
    Use output_routines, Only : init_hitline_netcdf, write_hitline_data_netcdf
    Use run_settings_namelist, Only : dmag, fname_launch, ns_line_diff, &
         fname_hit, fname_intpts, fname_nhit, nhitline, lambda_par
    Implicit none

    Logical, Parameter :: write_hitline_to_netcdf = .false.

    Integer(int32) :: numl, iline, ii, hitcount, ihit, iocheck
    Real(real64) :: Rstart, Zstart, Phistart

    Integer(int32) :: work_done, work_done_count
    Integer :: tag, dest
    Logical :: flag

    Real(real64), Dimension(3) :: pint
    Real(real64) :: totL, theta
    Integer(int32), Dimension(4) :: iout
    Real(real64), Allocatable :: R0(:),Z0(:),Phi0(:)
    Real(real64), Dimension(nhitline) :: r_hitline,z_hitline,phi_hitline
    Real(real64), Dimension(3) :: line_start_data_r
    Integer(int32), Dimension(3) :: line_start_data_i
    Real(real64), Dimension(size_rdata) :: line_done_data_r
    Real(real64), Dimension(3*nhitline) :: line_done_data_r2
    Integer(int32), Dimension(size_idata) :: line_done_data_i
    Real(real64) :: etime_follow, etime_int
    Integer(int32), Dimension(:), Allocatable :: this_job_done, is_working_arr
#ifdef USE_MPIF08
    Type(MPI_Request), Dimension(:), Allocatable :: req_arr
#else
    Integer, Dimension(:), Allocatable :: req_arr
#endif
    !- End of header -------------------------------------------------------------

    ! Read init point data
    Write(*,'(A,A)') ' Reading launch point data from ',Trim(Adjustl(fname_launch))
    Open(iu_launch,file=fname_launch)
    Read(iu_launch,*) numl
    Allocate(R0(numl),Z0(numl),Phi0(numl))
    Do ii = 1,numl
       Read(iu_launch,*) R0(ii),Z0(ii),Phi0(ii)
    End Do
    Close(iu_launch)

    Write(*,*) 'Total number of fieldlines to follow:',numl
    Write(*,*) 'Diffusing fieldlines with D_mag (m**2/m) = ',dmag
    If (lambda_par .ge. 0._real64) Then       
       Write(*,*) 'Parallel scattering mean free path (m) = ',lambda_par
    Else
       Write(*,*) 'Parallel scattering deactivated'
    End If
    Write(*,*)

    Allocate(this_job_done(nprocs-1),is_working_arr(nprocs-1),req_arr(nprocs-1))

    !-----------------------------------------------------------------------------
    ! Send initial fls and request results
    !-----------------------------------------------------------------------------
    iline = 0
    work_done_count = 0
    is_working_arr = 0
#ifdef USE_MPIF08
    req_arr(:) = MPI_REQUEST_NULL
#else
    req_arr(:) = 0
#endif
    
    Do dest = 1,nprocs - 1

       tag = dest

       ! Send line start data
       iline = iline + 1
       Rstart   = R0(iline)
       Zstart   = Z0(iline)
       Phistart = Phi0(iline)

       line_start_data_i(1) = ns_line_diff
       line_start_data_i(2) = nhitline
       line_start_data_i(3) = iline
       Call MPI_SEND(line_start_data_i,3,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,ierr_mpi)

       line_start_data_r(1) = Rstart
       line_start_data_r(2) = Zstart
       line_start_data_r(3) = Phistart
       Call MPI_SEND(line_start_data_r,3,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,ierr_mpi)

       ! Request signal that job is finished
       Call MPI_IRECV(this_job_done(dest),1,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,req_arr(dest),ierr_mpi)
       is_working_arr(dest) = 1

    End Do


    ! Set up output files
    Open(iu_int,file=fname_intpts,iostat=iocheck)
    Write(iu_int,*) '# R (m) | Z (m) | Phi (rad) | ihit | ipart | itri | i | Lc | sin(theta) | t_follow (s) | t_int (s)'

    ! Open hitline file
    If (write_hitline_to_netcdf) Then
       Call init_hitline_netcdf(fname_hit,nhitline)
    Else
       Open(iu_hit,file=fname_hit,iostat=iocheck)
    End If

    !-----------------------------------------------------------
    ! Loop over processors until work is complete
    !-----------------------------------------------------------
    work_done = 0
    hitcount = 0

    Do While ( work_done .ne. 1 )
       Do dest = 1,nprocs - 1

          ! If we are still waiting on a job from this process, check if it is finished
          If (is_working_arr(dest) == 1) Then
             Call MPI_TEST(req_arr(dest),flag,status,ierr_mpi)
             If (flag) Then
                is_working_arr(dest) = 0
                work_done_count = work_done_count + 1

                ! Request results
                tag = dest
                Call MPI_RECV(line_done_data_i,size_idata,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,status,ierr_mpi)
                Call MPI_RECV(line_done_data_r,size_rdata,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,status,ierr_mpi)

                If (nhitline .gt. 0) Then
                   Call MPI_RECV(line_done_data_r2,nhitline*3,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,status,ierr_mpi)
                Endif
                ihit        = line_done_data_i(2)
                etime_follow = line_done_data_r(6)
                etime_int = line_done_data_r(7)

                if ( ihit .ge. 1 ) Then
                   hitcount = hitcount + 1
                   iout = line_done_data_i(2:5)
                   pint = line_done_data_r(1:3) ! X,Y,Z
                   totL = line_done_data_r(4)
                   theta = line_done_data_r(5)

                   !R,Z,phi,ihit,ipart,itri,i,totL,theta
                   Write(iu_int,*)   sqrt(pint(1)*pint(1)+pint(2)*pint(2)),pint(3),atan2(pint(2),pint(1)),iout,totL,theta, &
                        etime_follow,etime_int


                   If (nhitline .gt. 0) Then
                      r_hitline   = line_done_data_r2(1+0*nhitline:1*nhitline)
                      z_hitline   = line_done_data_r2(1+1*nhitline:2*nhitline)
                      phi_hitline = line_done_data_r2(1+2*nhitline:3*nhitline)

                      If (write_hitline_to_netcdf) Then
                         Call write_hitline_data_netcdf(fname_hit, r_hitline, z_hitline, phi_hitline)
                      Else
                         Write(iu_hit,*) nhitline
                         Write(iu_hit,*) r_hitline
                         Write(iu_hit,*) z_hitline
                         Write(iu_hit,*) phi_hitline
                      Endif
                   End If
                Endif
             Endif ! flag true

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
                line_start_data_i(1) = ns_line_diff
                line_start_data_i(2) = nhitline
                line_start_data_i(3) = iline
                Call MPI_SEND(line_start_data_i,3,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,ierr_mpi)
                Call MPI_SEND(line_start_data_r,3,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,ierr_mpi)
                is_working_arr(dest) = 1

                ! Initiate request for completed job
                call MPI_IRECV(this_job_done(dest),1,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,req_arr(dest),ierr_mpi)
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
       Call MPI_SEND(line_start_data_i,3,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,ierr_mpi)
    Enddo


    Close(iu_hit)
    close(iu_int)

    Write(*,*) ' We had ',hitcount,' lines -hit-'

    Open(iu_nhit,file=fname_nhit)
    Write(iu_nhit,*) hitcount
    Close(iu_nhit)
    Deallocate(R0,Z0,Phi0)

  End Subroutine diffuse_lines3
  !-----------------------------------------------------------------------------




  !-----------------------------------------------------------------------------
  !+
  !-----------------------------------------------------------------------------
  ! Todo: cleanup
  Subroutine check_line_for_intersections(pint,iout, &
       r_hitline,z_hitline,phi_hitline,nhitline,linenum,nsteps_line,rout,zout,phiout,ifail,totL, &
       calc_lc,calc_theta,theta)

    Use kind_mod, Only : real64, int32
    Use read_parts_mod
    Use inside_vessel_mod, Only : inside_vessel, find_vessel_intersection, &
         init_find_vessel_intersection, fin_find_vessel_intersection
    Use math_routines_mod, Only : line_seg_facet_int, dist_2pts_cyl, wrap_phi, int_line_curve
    Use parallel_mod, Only : fin_mpi
    Use run_settings_namelist, Only : period, lsfi_tol, vessel_int_is_last_point
    Implicit None

    Integer(int32), Intent(in) :: linenum, nsteps_line,ifail
    Integer(int32), Intent(out), Dimension(4) :: iout
    Real(real64), Dimension(3), Intent(out) :: pint
    Real(real64), Intent(out) :: totL, theta
    Integer(int32), Intent(in) :: nhitline
    Real(real64), Intent(out),Dimension(nhitline) :: r_hitline,z_hitline,phi_hitline
    Real(real64), Dimension(nsteps_line+1),intent(in) :: rout,zout,phiout
    Logical, Intent(In) :: calc_lc, calc_theta

    Integer(int32) :: iseg, ierr_pint
    Integer(int32) :: npts_line, ihit, i, twofer, inphi, ntri, ihit_tmp, ipart, itri
    Logical :: inside_it
    Integer(int32), Dimension(1) :: ind_min, ind_max
    Real(real64) :: R1, Z1, P1, P1a, P2a, X3, Y3, Z3, R3, mu, Aplane, Bplane, denom, pint2D(2), rint,zint,uint
    Real(real64) :: p_start, x_start, y_start, z_start, p_end, x_end, y_end, z_end,r_start,r_end,R2a,R1a
    Real(Real64), Dimension(2) :: Rtmp, Ztmp, Ytmp, Xtmp
    Real(real64) :: R2,Z2,P2, X1, Y1, X2, Y2, dphi1, dphi2, X1a, Y1a, Z1a, X2a, Y2a, Z2a
    Real(real64) :: RLast, ZLast, PLast
    Real(real64), Dimension(3) :: Pt1, pt2, pa, pb, pc
    Real(real64) :: Pmin_line, Pmax_line, Pmin_test, Pmax_test
    !- End of header -------------------------------------------------------------

    ! Set length of line to check. I
    npts_line = nsteps_line + 1
    If (ifail .ne. nsteps_line) npts_line = ifail - 1

    ! Initialize variables 
    ihit = 0
    pint(:) = 0._real64
    iout(:) = -1
    iout(1) = ihit
    totL = 0._real64
    r_hitline(:) = 0._real64
    z_hitline(:) = 0._real64
    phi_hitline(:) = 0._real64

    ! Loop over fieldline points
    Do i=1,npts_line - 1

       ! Current point along line
       R1 = rout(i)
       Z1 = zout(i)
       P1 = phiout(i)

       ! next point
       R2 = rout(i+1)
       Z2 = zout(i+1)
       P2 = phiout(i+1)

       ! Initial dphi used to check if line segment crosses symmetry plane
       dphi1 = P2-P1

       ! Distance
       If (calc_lc) Then
          totL = totL + dist_2pts_cyl(R1, R2, Z1, Z2, P1, P2)
       Endif

       ! Convert to Cartesian coordinates
       Call wrap_phi(P1,period)
       X1 = R1*cos(P1)
       Y1 = R1*sin(P1)

       Call wrap_phi(P2,period)
       X2 = R2*cos(P2)
       Y2 = R2*sin(P2)

       ! Post wrapping dphi
       dphi2 = P2-P1

       ! Check for lines that cross symmetry (field period) plane
       twofer = 0
       If ( abs(dphi1-dphi2) .gt. 1.d-4) Then
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
          R1a = Rtmp(ind_min(1))
          
          X2a = Xtmp(ind_max(1))
          Y2a = Ytmp(ind_max(1))
          Z2a = Ztmp(ind_max(1))
          R2a = Rtmp(ind_max(1))

          !---------------------------------------------
          ! Find intersection of line segment with plane
          !---------------------------------------------
          pt1(1) = Rtmp(ind_min(1))*cos(P1a+period)
          pt1(2) = Rtmp(ind_min(1))*sin(P1a+period)
          pt1(3) = Ztmp(ind_min(1))
          pt2(1) = Xtmp(ind_max(1))
          pt2(2) = Ytmp(ind_max(1))
          pt2(3) = Ztmp(ind_max(1))

          ! Define explicit form of plane at 2*pi/nfp (period)
          pa = [0.d0,0.d0,0.d0]
          pb = [0.d0,0.d0,1.d0]
          pc = [cos(period),sin(period),0.d0]

          ! Calc unit vector normal to plane of Pa-c (gives plane components A-C)
          Aplane=-sin(period)
          Bplane= cos(period)

          ! Calculate the position on the line that intersects the plane
          denom = Aplane*(pt2(1) - pt1(1)) + Bplane*(pt2(2) - pt1(2))
          If (abs(denom) .lt. 1.d-15) Then
             Write(*,*) 'Error: Line did not hit plane at 2*pi/nfp (is parallel based on dot product)'
             Call fin_mpi(.true.)
          Endif
          mu = - (Aplane*pt1(1) + Bplane*pt1(2)) / denom
          X3 = pt1(1) + mu * (pt2(1) - pt1(1))
          Y3 = pt1(2) + mu * (pt2(2) - pt1(2))
          Z3 = pt1(3) + mu * (pt2(3) - pt1(3))
          R3 = sqrt(X3*X3+Y3*Y3)

          ! Define two new points
          P1 = 0._real64
          R1 = R3
          X1 = R3*cos(P1)
          Y1 = R3*sin(P1)
          Z1 = Z3

          P2 = period
          R2 = R3
          X2 = R3*cos(P2)
          Y2 = R3*sin(P2)
          Z2 = Z3
       End If

       Do iseg = 1,1+twofer

          If (twofer .eq. 0) Then
             p_start = P1
             x_start = X1
             y_start = Y1
             z_start = Z1
             r_start = R1
             p_end   = P2
             x_end   = X2
             y_end   = Y2
             z_end   = Z2
             r_end   = R2
          Else
             if ( iseg .eq. 1) Then
                p_start = P1
                x_start = X1
                y_start = Y1
                z_start = Z1
                r_start = R1
                p_end   = P1a
                x_end   = X1a
                y_end   = Y1a
                z_end   = Z1a
                r_end   = R1a
             else
                p_start = P2a
                x_start = X2a
                y_start = Y2a
                z_start = Z2a
                r_start = R2a
                p_end   = P2
                x_end   = X2
                y_end   = Y2
                z_end   = Z2
                r_end   = R2
             End If
          End If

          Pmin_line = Min(p_start,p_end)
          Pmax_line = Max(p_start,p_end)

          ihit = 0
          Do ipart = 1,nparts

             ! Check if this line segment is within the phi bounds of the part
             ! The phi segments overlap if Max(Min(phi_line),Min(phi_test)) <= Min(Max(phi_line),Max(phi_test))
             Pmin_test = Min(Pmins(ipart),Pmaxs(ipart))
             Pmax_test = Max(Pmins(ipart),Pmaxs(ipart))

             inphi = 0
             If (Max(Pmin_line, Pmin_test) <= Min(Pmax_line, Pmax_test)) Then
                inphi = 1
             End If
             If (inphi .eq. 1 ) Then

                If (is_AS_part(ipart)) Then
                   ! Part is AS (allowing for it to have a finite toroidal extent)
                   Call int_line_curve((/r_start,z_start/),(/r_end,z_end/), &
                        Rparts(ipart,1,1:np_parts(ipart)),Zparts(ipart,1,1:np_parts(ipart)), &
                        .true.,pint2D,ierr_pint,uint)
                   If (ierr_pint .eq. 0) Then
                      ihit = 1
                      pint(1) = pint2D(1)*cos(uint*(p_end-p_start)+p_start)
                      pint(2) = pint2D(1)*sin(uint*(p_end-p_start)+p_start)
                      pint(3) = pint2D(2)
                      itri = -1
                   End If
                Else
                   ! If not AS then check triangles
                   ntri = ntri_parts(ipart)
                   Do itri = 1,ntri

                      ! Check if this line segment is within the phi bounds of the triangle
                      ! The phi segments overlap if Max(Min(phi_line),Min(phi_test)) <= Min(Max(phi_line),Max(phi_test))
                      Pmin_test = Min(Pmintri(ipart,itri),Pmaxtri(ipart,itri))
                      Pmax_test = Max(Pmintri(ipart,itri),Pmaxtri(ipart,itri))

                      If (Max(Pmin_line, Pmin_test) <= Min(Pmax_line, Pmax_test)) Then

                         pa = [xtri(ipart,itri,1),ytri(ipart,itri,1),ztri(ipart,itri,1)]
                         pb = [xtri(ipart,itri,2),ytri(ipart,itri,2),ztri(ipart,itri,2)]
                         pc = [xtri(ipart,itri,3),ytri(ipart,itri,3),ztri(ipart,itri,3)]
                         pt1 = [x_start,y_start,z_start]
                         pt2 = [x_end,y_end,z_end]

                         Call line_seg_facet_int(pa,pb,pc,pt1,pt2,ihit_tmp,pint,lsfi_tol,calc_theta,theta)

                         If (ihit_tmp .eq. 1) Then
                            ihit = 1
                            Exit ! stop looking for triangle intersections
                         End If

                      End If ! Is in triangle phi range
                   End Do !triangle loop
                End If ! is AS part

                ! If it hit then write pint and hitline
                If (ihit .eq. 1 ) Then
                   Write(*,'(A,I0,A,I0,A,I0,1X,4(G0.3,1X))') ' Line ',linenum,' hit! [i,ipart,Px,Py,Pz,Lc] ',&
                        i,' ',ipart,pint(1),pint(2),pint(3),totL
                   iout(1) = ihit
                   iout(2) = ipart
                   iout(3) = itri
                   iout(4) = i

                   If ( (i - nhitline+1) .lt. 1 ) Then
                      Write(*,*) 'Truncating hitline'
                      r_hitline(1:i)   =   rout(1:i)
                      z_hitline(1:i)   =   zout(1:i)
                      phi_hitline(1:i) = phiout(1:i)
                   Else
                      r_hitline(1:nhitline)   =   rout(i-nhitline+1:i)
                      z_hitline(1:nhitline)   =   zout(i-nhitline+1:i)
                      phi_hitline(1:nhitline) = phiout(i-nhitline+1:i)
                   End If
                   Exit  ! Stop looking for part intersections
                End If

             End If !inphi check
          Enddo ! part index

          If (ihit .eq. 1 ) Exit  ! Stop looking for segment intersections
       Enddo ! seg index (twofer)

       If (ihit .eq. 1 ) Exit  ! Quit fieldline if it hit a part
    Enddo ! points along line (i)



    If ( ihit .eq. 0 ) Then
       Write(*,'(A,I0,A)') ' No part intersections found for line ',linenum,', checking for vessel intersection'

       ! Check if first point is outside vessel, if so quit
       i = 1
       RLast = rout(i)
       ZLast = zout(i)
       PLast = phiout(i)
       inside_it = inside_vessel(RLast,ZLast,PLast,R_ves,Z_ves,P_ves,ntor_ves,npol_ves)
       If (.not. inside_it) Then
          Write(*,*) 'Error: initial point on fl',linenum,' is already outside vessel!'
          Write(*,*) 'Point (R,Z,Phi) = ',RLast,ZLast,PLast
          Call fin_mpi(.true.)
       Endif

       ! Allocate slice arrays for find_vessel_intersection
       If (.not. vessel_int_is_last_point) Call init_find_vessel_intersection

       ! Loop over curve points
       totL = 0._real64
       Do i=2,npts_line

          ! Current point along line
          R1 = rout(i)
          Z1 = zout(i)
          P1 = phiout(i)

          ! Calculate connection length
          If (calc_lc) Then
             totL = totL + dist_2pts_cyl(R1, RLast, Z1, ZLast, P1, PLast)
          Endif

          ! Check if current point is inside vessel
          inside_it = inside_vessel(R1,Z1,P1,R_ves,Z_ves,P_ves,ntor_ves,npol_ves)

          ! Three cases
          ! 1) Both points outside. Could intersect poly but we assume the fl started inside.
          ! 2) Both points inside. Cannot intersect
          ! 3) One in and one out. Intersects. Since we started inside we only need to check inside_it.

          ! If there is an intersection find int point and define hitline
          If (.not. inside_it) Then

             ! Define vessel intersection point
             If (vessel_int_is_last_point) Then
                ! Just use first point on fl outside of vessel
                Call wrap_phi(P1,period)
                pint(1) = R1*cos(P1)
                pint(2) = R1*sin(P1)
                pint(3) = Z1
             Else
                ! Check for intersection with vessel
                ! This is not perfect because it just uses the vessel at Phi1
                Call find_vessel_intersection(R1,Z1,RLast,ZLast,P1,Rint,Zint)
                pint(1) = Rint*cos(P1)
                pint(2) = Rint*sin(P1)
                pint(3) = Zint
             End If
             Write(*,*) ' Line ',linenum,' did hit the vessel at [i,P,Lc]=',i,pint,totL

             ! Set the MPI integer output
             ihit = 2
             iout(1) = ihit
             iout(2) = -2
             iout(3) = -2
             iout(4) = i

             ! Set hitline
             if ( (i - nhitline+1) .lt. 1 ) Then
                r_hitline(1:i) = rout(1:i)
                z_hitline(1:i) = zout(1:i)
                phi_hitline(1:i) = phiout(1:i)
             Else
                r_hitline = rout(i-nhitline+1:i)
                z_hitline = zout(i-nhitline+1:i)
                phi_hitline = phiout(i-nhitline+1:i)
             Endif
             Exit
          Endif ! not inside_it

          ! Save last point
          RLast = R1
          ZLast = Z1
          PLast = P1

       Enddo ! npts line

       ! Allocate slice arrays for find_vessel_intersection
       If (.not. vessel_int_is_last_point) Call fin_find_vessel_intersection

    Endif ! did not hit part

  End subroutine check_line_for_intersections


End Module diffusion

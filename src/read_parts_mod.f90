!-----------------------------------------------------------------------------
!+ Contains routines used to read and interact with part files.
!-----------------------------------------------------------------------------
Module read_parts_mod
  ! Description:
  !
  ! Author(s): J.D. Lore - 07/14/2011 - xxx
  !

  Use kind_mod, Only : real64, int32

  Implicit None

  ! Parts
  Real(real64), Allocatable :: Rparts(:,:,:),Zparts(:,:,:),Pparts(:,:,:)
  Real(real64), Allocatable, Dimension(:) :: Pmaxs, Pmins
  Logical, Allocatable, Dimension(:) :: is_AS_part, force_non_AS
  Integer(int32) :: nparts, nt_max, np_max
  Real(real64), Allocatable :: near_part(:)
  Integer(int32),Allocatable,Dimension(:) :: nt_parts, np_parts, part_type
  Character(len=300),Allocatable,Dimension(:) :: part_names
  Integer(int32) :: ic_near

  ! Vessel
  Real(real64), Allocatable :: R_ves(:,:),Z_ves(:,:),P_ves(:,:)
  Integer(int32) :: ntor_ves, npol_ves
  Logical :: is_AS_ves, force_non_AS_ves

  ! Triangles
  Integer(int32) :: ntri_max
  Real(real64), Allocatable, Dimension(:,:) :: xmid, ymid, zmid, dmid, pmintri, pmaxtri
  Real(real64), Allocatable,Dimension(:,:,:) :: xtri,ytri,ztri
  Integer(int32), Allocatable,Dimension(:) :: ntri_parts
  Real(real64), Allocatable :: near_tri(:)

  !- End of header -------------------------------------------------------------

Contains

  !-----------------------------------------------------------------------------
  !+ Makes triangles from 2d parts
  !-----------------------------------------------------------------------------
  Subroutine make_triangles(verbose)
    Use kind_mod, Only : real64, int32
    Use io_unit_spec, Only : iu_ptri, iu_ptmid
    Use run_settings_namelist, Only : period, fname_ptri, fname_ptri_mid
    Use math_routines_mod, Only : wrap_phi
    Use phys_const, Only : pi
    Use parallel_mod
    Implicit None
    Logical, Intent(in) :: verbose
    Real(real64) :: R1,Z1,P1,Pt1(3)
    Real(real64) :: R2,Z2,P2,Pt2(3)
    Real(real64) :: R3,Z3,P3,Pt3(3)
    Real(real64) :: R4,Z4,P4,Pt4(3), dmids(3)
    Integer(int32) :: ipart, jpol, itor, itri, npol, ntor, ifacet, itri_tot
    Real(real64), Allocatable :: rtri_part(:,:), ptri_part(:,:)
    Real(real64), Allocatable :: xtri_tmp(:,:),ytri_tmp(:,:),ztri_tmp(:,:)
    !- End of header -------------------------------------------------------------

    ! First have to get max number of triangles for allocation
    Allocate(ntri_parts(nparts))
    Do ipart=1,nparts
       If (part_type(ipart) .eq. 2) Then
          Call query_tri_part(part_names(ipart),ntri_parts(ipart))
       Else
          ntor = nt_parts(ipart)
          npol = np_parts(ipart)
          ntri_parts(ipart) = (ntor-1)*(npol-1)*2
       End If
    End Do
    ntri_max = Maxval(ntri_parts)

    Allocate(xtri(nparts,ntri_max,3) ,source=0._real64)
    Allocate(ytri(nparts,ntri_max,3) ,source=0._real64)
    Allocate(ztri(nparts,ntri_max,3) ,source=0._real64)
    Allocate(pmintri(nparts,ntri_max),source=0._real64)
    Allocate(pmaxtri(nparts,ntri_max),source=0._real64)

    ! Create triangles
    itri_tot = 0
    Do ipart = 1,nparts

       If (part_type(ipart) .eq. 2) Then
          ! For triangle parts
          Allocate(xtri_tmp(ntri_parts(ipart),3))
          Allocate(ytri_tmp(ntri_parts(ipart),3))
          Allocate(ztri_tmp(ntri_parts(ipart),3))
          Call load_tri_part(part_names(ipart),ntri_parts(ipart),xtri_tmp,ytri_tmp,ztri_tmp)

          xtri(ipart,1:ntri_parts(ipart),1:3) = xtri_tmp
          ytri(ipart,1:ntri_parts(ipart),1:3) = ytri_tmp
          ztri(ipart,1:ntri_parts(ipart),1:3) = ztri_tmp
          Deallocate(xtri_tmp,ytri_tmp,ztri_tmp)

          itri_tot = itri_tot + ntri_parts(ipart)

          ! Convert to cylindrical coordinates
          Allocate(rtri_part(ntri_max,3) ,source=0._real64)
          Allocate(ptri_part(ntri_max,3) ,source=0._real64)
          rtri_part(:,:) = Sqrt(xtri(ipart,:,:)**2 + ytri(ipart,:,:)**2)
          ptri_part(:,:) = Atan2(ytri(ipart,:,:),xtri(ipart,:,:))

          ! This section (and back conversion) commented out because how do we set periodicity
          ! for triangle parts? For now assume this has been done.
          !          ! Map to first period
          !          Do itri = 1,ntri_parts(ipart)
          !             Do j = 1,3
          !                Call wrap_phi(ptri_part(itri,j),2._real64*pi/Real(msym,real64))
          !             End Do
          !          End Do

          ! Get min/max of Phi for filtering out intersection checks
          Pmins(ipart) = Minval(ptri_part(1:ntri_parts(ipart),:))
          Pmaxs(ipart) = Maxval(ptri_part(1:ntri_parts(ipart),:))
          Do itri = 1,ntri_parts(ipart)
             pmintri(ipart,itri) = Minval(ptri_part(itri,:))
             pmaxtri(ipart,itri) = Maxval(ptri_part(itri,:))
          End Do

          !          ! Convert back to cartesian
          !          xtri(ipart,:,:) = rtri_part*Cos(ptri_part)
          !          ytri(ipart,:,:) = rtri_part*Sin(ptri_part)

          If (verbose) Then
             Write(*,*) 'Triangle part ',ipart,'extends from Phi = ',Pmins(ipart)*180./pi,' to ',Pmaxs(ipart)*180./pi,' deg.'
             If (Pmaxs(ipart) .gt. period) Then
                Write(*,*) 'Warning: Part extends beyond Bfield period, extra range is not used!'
             End If
          End If
          Deallocate(rtri_part)
          Deallocate(ptri_part)

       Else
          ! For "jparts" and "parts"
          ntor = nt_parts(ipart)
          npol = np_parts(ipart)
          itri = 1
          Do itor=1,ntor-1
             Do jpol=1,npol-1

                ! define 4 points for this face
                R1 = Rparts(ipart,itor,jpol)
                Z1 = Zparts(ipart,itor,jpol)
                P1 = Pparts(ipart,itor,jpol)
                Pt1(1:3) = [R1*cos(P1),R1*sin(P1),Z1]
                R2 = Rparts(ipart,itor+1,jpol)
                Z2 = Zparts(ipart,itor+1,jpol)
                P2 = Pparts(ipart,itor+1,jpol)
                Pt2(1:3) = [R2*cos(P2),R2*sin(P2),Z2]
                R3 = Rparts(ipart,itor,jpol+1)
                Z3 = Zparts(ipart,itor,jpol+1)
                P3 = Pparts(ipart,itor,jpol)
                Pt3(1:3) = [R3*cos(P3),R3*sin(P3),Z3]
                R4 = Rparts(ipart,itor+1,jpol+1)
                Z4 = Zparts(ipart,itor+1,jpol+1)
                P4 = Pparts(ipart,itor+1,jpol)
                Pt4(1:3) = [R4*cos(P4),R4*sin(P4),Z4]

                ! define triangles for this part

                ! Loop over the two facets from each quad
                Do ifacet = 1,2
                   If (ifacet .eq. 1) Then
                      ! 1st triangle
                      xtri(ipart,itri,1:3) = [Pt1(1),Pt2(1),Pt3(1)]
                      ytri(ipart,itri,1:3) = [Pt1(2),Pt2(2),Pt3(2)]
                      ztri(ipart,itri,1:3) = [Pt1(3),Pt2(3),Pt3(3)]
                      ! Set phi range of triangle for checking
                      pmintri(ipart,itri) = Minval([P1,P2,P3])
                      pmaxtri(ipart,itri) = Maxval([P1,P2,P3])
                   Else
                      ! 2nd triangle
                      xtri(ipart,itri,1:3) = [Pt4(1),Pt3(1),Pt2(1)]
                      ytri(ipart,itri,1:3) = [Pt4(2),Pt3(2),Pt2(2)]
                      ztri(ipart,itri,1:3) = [Pt4(3),Pt3(3),Pt2(3)]
                      ! Set phi range of triangle for checking
                      pmintri(ipart,itri) = Minval([P4,P3,P2])
                      pmaxtri(ipart,itri) = Maxval([P4,P3,P2])
                   End If

                   itri_tot = itri_tot + 1
                   itri = itri + 1
                End Do ! ifacet
             End Do !jpol
          End Do ! ntor
       End If ! part type
    End Do !part

    If (verbose) Write(*,*) "Total number of triangles:",itri_tot

    ! Write full triangle file
    Open(iu_ptri,file=fname_ptri)
    Write(iu_ptri,*) nparts
    Do ipart = 1,nparts
       write(iu_ptri,*) ipart,ntri_parts(ipart)
       Do itri = 1,ntri_parts(ipart)
          Write(iu_ptri,*) itri
          Write(iu_ptri,*) xtri(ipart,itri,1),ytri(ipart,itri,1),ztri(ipart,itri,1)
          Write(iu_ptri,*) xtri(ipart,itri,2),ytri(ipart,itri,2),ztri(ipart,itri,2)
          Write(iu_ptri,*) xtri(ipart,itri,3),ytri(ipart,itri,3),ztri(ipart,itri,3)
       End Do
    End Do
    Close(iu_ptri)

    ! Compute and write triangle midpoints
    Allocate(xmid(nparts,ntri_max)   ,source=0._real64)
    Allocate(ymid(nparts,ntri_max)   ,source=0._real64)
    Allocate(zmid(nparts,ntri_max)   ,source=0._real64)
    Allocate(dmid(nparts,ntri_max)   ,source=0._real64)

    Open(iu_ptmid,file=fname_ptri_mid)
    Write(iu_ptmid,*) nparts
    Do ipart = 1,nparts
       Write(iu_ptmid,*) ipart,ntri_parts(ipart)
       Do itri = 1,ntri_parts(ipart)
          ! Define triangle midpoint
          xmid(ipart,itri) = Sum(xtri(ipart,itri,1:3))/3._real64
          ymid(ipart,itri) = Sum(ytri(ipart,itri,1:3))/3._real64
          zmid(ipart,itri) = Sum(ztri(ipart,itri,1:3))/3._real64
          dmids(1) = Sqrt( (xtri(ipart,itri,1) - xmid(ipart,itri))**2 + &
               (ytri(ipart,itri,1) - ymid(ipart,itri))**2 + (ztri(ipart,itri,1) - zmid(ipart,itri))**2 )
          dmids(2) = Sqrt( (xtri(ipart,itri,2) - xmid(ipart,itri))**2 + &
               (ytri(ipart,itri,2) - ymid(ipart,itri))**2 + (ztri(ipart,itri,2) - zmid(ipart,itri))**2 )
          dmids(3) = Sqrt( (xtri(ipart,itri,3) - xmid(ipart,itri))**2 + &
               (ytri(ipart,itri,3) - ymid(ipart,itri))**2 + (ztri(ipart,itri,3) - zmid(ipart,itri))**2 )
          dmid(ipart,itri) = Maxval(dmids,1)
          write(iu_ptmid,*) itri,xmid(ipart,itri),ymid(ipart,itri),zmid(ipart,itri),dmid(ipart,itri)
       End Do
    End Do
    Close(iu_ptmid)



  End Subroutine make_triangles
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !+ Reads 2d part file (w7 format)
  !-----------------------------------------------------------------------------
  Subroutine load_w7_part(fname,label,ntor,npol,msym,Rpart,Zpart,Phipart,force_non_AS)
    Use kind_mod
    Use io_unit_spec, Only : iu_thispart
    Use phys_const, Only : pi
    Implicit none
    ! Dummy variables
    Character(len=300), Intent(in) :: fname
    Integer(int32), Intent(in) :: ntor, npol
    Character(len=300), Intent(out) :: label
    Integer(int32),Intent(out) :: msym
    Real(real64),Dimension(ntor,npol), Intent(out) :: Rpart, Zpart, Phipart
    Logical, Intent(out) :: force_non_AS
    Integer(int32) :: itor, ipol, ntor_dum, npol_dum, iostat
    Real(real64) :: rshift, zshift, Phitmp
    Character(len=300) :: line
    !- End of header -------------------------------------------------------------

    open(iu_thispart,file=fname)
    Read(iu_thispart,*) label

    ! Parse header, allowing for optional logical force_non_AS
    force_non_AS = .false.
    Read(iu_thispart, '(A)') line
    Read(line, *, IOSTAT=iostat) ntor_dum, npol_dum, msym, rshift, zshift, force_non_AS

    ! When read [R,Z] = cm and Phi = degrees
    Do itor = 1,ntor
       Read(iu_thispart,*) Phitmp
       Phipart(itor,:) = Phitmp
       Do ipol = 1,npol
          Read(iu_thispart,*) Rpart(itor,ipol), Zpart(itor,ipol)
       Enddo
    Enddo
    close(iu_thispart)

    ! Convert to meters and radians
    Rpart = (Rpart+rshift)*0.01_real64
    Zpart = (Zpart+zshift)*0.01_real64
    Phipart = Phipart*pi/180._real64

  End Subroutine load_w7_part
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !+ Queries triangle part file for number of triangles, ignoring blank lines
  !-----------------------------------------------------------------------------
  Subroutine query_tri_part(fname,num_tri)
    Use kind_mod, Only : int32

    Use parallel_mod, Only : fin_mpi
    Use io_unit_spec, Only : iu_thispart
    Implicit None
    Character(len=300), Intent(in) :: fname
    Integer(int32), Intent(out) :: num_tri
    Integer(int32) :: iostat
    Character(len=500) :: dummy_line
    !- End of header -------------------------------------------------------------

    Open(iu_thispart,file=fname,status='old',iostat=iostat)
    If (iostat /= 0) Then
       Write(*,*) "Error: Unable to open file: ",trim(fname)
       Call fin_mpi(.true.) ! True means this is an exit-on-error
    Endif

    num_tri = 0

    Do
       Read(iu_thispart, '(A)', iostat=iostat) dummy_line
       If (iostat /= 0) Exit
       If (len_trim(dummy_line) == 0) Cycle
       num_tri = num_tri + 1
    End Do
    Close(iu_thispart)

  End Subroutine query_tri_part
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !+ Reads triangle part file
  !-----------------------------------------------------------------------------
  Subroutine load_tri_part(fname,num_tri,xtri,ytri,ztri)
    Use kind_mod, Only : int32, real64
    Use parallel_mod, Only : fin_mpi
    Use io_unit_spec, Only : iu_thispart
    Implicit None
    Character(len=300), Intent(in) :: fname
    Integer(int32), Intent(in) :: num_tri
    Real(real64), Intent(out) :: xtri(num_tri,3),ytri(num_tri,3),ztri(num_tri,3)
    Integer(int32) :: iostat, i, id, pos, pos_open, pos_close, comma1, comma2, itri, sp
    character(len=500)  :: line, token
    Real(real64) :: v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z
    real(8), dimension(3,3) :: v  ! v(i,1:3) for each vertex
    !- End of header -------------------------------------------------------------

    Open(iu_thispart,file=fname,status='old',iostat=iostat)
    If (iostat /= 0) Then
       Write(*,*) "Error: Unable to open file: ",trim(fname)
       Call fin_mpi(.true.) ! True means this is an exit-on-error
    Endif

    Do itri = 1, num_tri

       Read(iu_thispart, '(A)', iostat=iostat) line

       ! Extract the integer (identified as preceeding first "(")
       sp = index(line, '(')
       if (sp == 0) then
          print *, "Error: No ( found to separate the integer."
          stop
       endif
       read(line(1:sp-1), *) id

       ! Loop over each vertex group
       pos = sp
       do i = 1, 3
          ! Find the next "("
          pos_open = index(line(pos:), '(')
          if (pos_open == 0) then
             print *, "Error: Missing '(' for vertex ", i
             stop
          endif
          pos_open = pos + pos_open - 1

          ! Find the corresponding ")"
          pos_close = index(line(pos_open:), ')')
          if (pos_close == 0) then
             print *, "Error: Missing ')' for vertex ", i
             stop
          endif
          pos_close = pos_open + pos_close - 1

          ! Extract the substring inside the parentheses.
          token = trim(line(pos_open+1:pos_close-1))
          read(token, *) v(i,:)

          ! Move pos beyond the closing parenthesis for next iteration.
          pos = pos_close + 1
       end do

       ! Store each vertex
       xtri(itri,1) = v(1,1);  xtri(itri,2) = v(2,1);  xtri(itri,3) = v(3,1)
       ytri(itri,1) = v(1,2);  ytri(itri,2) = v(2,2);  ytri(itri,3) = v(3,2)
       ztri(itri,1) = v(1,3);  ztri(itri,2) = v(2,3);  ztri(itri,3) = v(3,3)
    End Do
    Close(iu_thispart)


  End Subroutine load_tri_part
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !+ Reads 2d jpart file
  !-----------------------------------------------------------------------------
  Subroutine load_2d_jpart(fname,label,ntor,npol,msym,Rpart,Zpart,Phipart,force_non_AS)
    Use kind_mod, Only : int32, real64
    Use parallel_mod, Only : fin_mpi
    Use io_unit_spec, Only : iu_thispart
    Use phys_const, Only : pi
    Implicit none

    ! Dummy variables
    Character(len=300), Intent(in) :: fname
    Integer(int32), Intent(in) :: ntor,npol
    Character(len=300), Intent(out) :: label
    Integer(int32),Intent(out) :: msym
    Real(real64),Dimension(ntor,npol), Intent(out) :: Rpart, Zpart, Phipart
    Logical, Intent(out) :: force_non_AS

    ! Local vars
    Integer(int32) :: itor, ipol, ntor_dum, npol_dum, iostat
    Real(real64) :: rshift, zshift
    Character(len=256) :: line
    !- End of header -------------------------------------------------------------

    Open(iu_thispart,file=fname,status='old',iostat=iostat)
    If (iostat /= 0) Then
       Write(*,*) "Error: Unable to open file: ",trim(fname)
       Call fin_mpi(.true.) ! True means this is an exit-on-error
    Endif
    Read(iu_thispart,*) label

    !Read(iu_thispart,*) ntor_dum, npol_dum, msym, rshift, zshift
    force_non_AS = .false.
    Read(iu_thispart, '(A)') line
    Read(line, *, IOSTAT=iostat) ntor_dum, npol_dum, msym, rshift, zshift, force_non_AS

    ! When read [R,Z] = cm and Phi = degrees
    Do itor = 1,ntor
       Do ipol = 1,npol
          Read(iu_thispart,*) Rpart(itor,ipol), Zpart(itor,ipol), Phipart(itor,ipol)
       Enddo
    Enddo
    close(iu_thispart)

    ! Convert to meters and radians
    Rpart = (Rpart+rshift)*0.01_real64
    Zpart = (Zpart+zshift)*0.01_real64
    Phipart = Phipart*pi/180._real64

  Endsubroutine load_2d_jpart
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------




  !-----------------------------------------------------------------------------
  !+ Get dimensions from part file
  !-----------------------------------------------------------------------------
  Subroutine query_part(fname,ntor,npol,nfp_part)
    Use kind_mod, Only : int32, real64
    Use io_unit_spec, Only : iu_thispart
    Use parallel_mod, Only : fin_mpi
    Implicit none
    Character(len=300), Intent(in) :: fname
    Integer(int32), Intent(out) :: ntor,npol,nfp_part
    Integer(int32) :: iostat
    Real(real64) :: rshift, zshift
    Character(len=300) :: label
    !- End of header -------------------------------------------------------------
    Open(iu_thispart,file=fname,status='old',iostat=iostat)
    If (iostat /= 0) Then
       Write(*,*) "Error: Unable to open file: ",trim(fname)
       Call fin_mpi(.true.) ! True means this is an exit-on-error
    Endif
    Read(iu_thispart,*) label
    Read(iu_thispart,*) ntor, npol, nfp_part, rshift, zshift
    Close(iu_thispart)
  Endsubroutine query_part
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------




  !-----------------------------------------------------------------------------
  !+ Closes a part with only two points
  !-----------------------------------------------------------------------------
  Subroutine close_2pt_part(Rin,Zin,ntor,npol,dir,step,Rout,Zout)
    Use kind_mod
    Implicit none
    Integer(int32), Intent(in) :: ntor, npol, dir
    Real(real64), Intent(in) :: step
    Real(real64), Dimension(ntor,npol), Intent(in) :: Rin,Zin
    Real(real64), Dimension(ntor,npol+1), Intent(out) :: Rout, Zout
    Integer(int32) :: i
    Real(real64) :: R1,Z1,R2,Z2,dL,Rmid,Zmid
    Real(real64), Dimension(2) ::  V1,unit_out
    !- End of header -------------------------------------------------------------

    Rout(1:ntor,1:npol) = Rin
    Zout(1:ntor,1:npol) = Zin

    Do i=1,ntor
       R1 = Rin(i,1)
       Z1 = Zin(i,1)

       R2 = Rin(i,npol)
       Z2 = Zin(i,npol)
       dL = Sqrt( (R1-R2)*(R1-R2) + (Z1-Z2)*(Z1-Z2) )

       V1(1) = R2-R1
       V1(2) = Z2-Z1
       V1 = V1/dL

       Rmid = R1 + V1(1)*dL/2._real64
       Zmid = Z1 + V1(2)*dL/2._real64

       unit_out(1) = -V1(2)
       unit_out(2) =  V1(1)
       unit_out = unit_out/Sqrt( V1(1)*V1(1) + V1(2)*V1(2) )

       Rout(i,npol+1) = Rmid + real(dir,real64)*unit_out(1)*step
       Zout(i,npol+1) = Zmid + real(dir,real64)*unit_out(2)*step
       !  Rout(i,npol+1) = Rmid
       !  Zout(i,npol+1) = Zmid

       !  Write(*,*) '------'
       !  Write(*,*) Rout(i,:)
       !  Write(*,*) Zout(i,:)
       !  Write(*,*) '------'

    Enddo
  End Subroutine close_2pt_part
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !+ Reads parts list and load part coordinates
  !-----------------------------------------------------------------------------
  Subroutine read_parts(verbose)
    !
    ! Description:
    !
    !  Notes:
    !    part_type == 0 indicates 'w7' type. (R,Z) points at toroidal angles
    !    part_type == 1 indicates '2d.jpart' type. (R,Z,phi) points
    !    part_type == 2 indicates 'triangle type', inferred from .txt extension
    !
    ! Author(s): J.D. Lore - 07/14/2011 - xxx
    !
    Use kind_mod, Only : int32, real64
    Use io_unit_spec, Only: iu_plist, iu_parts
    Use phys_const, Only : pi
    Use run_settings_namelist, Only : period, fname_plist, fname_parts, fname_ves
    Use math_routines_mod, Only : wrap_phi
    Use parallel_mod, Only : fin_mpi
    Implicit none
    Logical, Intent(in) :: verbose
    Real(real64),Allocatable :: Rpart(:,:),Zpart(:,:),Ppart(:,:)
    Real(real64) :: check_AS
    Integer(int32) :: i, j, ipart, ntor, npol, nfp_part, msym, msym_ves, iostat, name_len
    Character(len=300) :: part_name, label

    ! Parameters
    ! Tolerance on checking phi limits against bfield periodicity.
    Real(real64), Parameter :: phi_period_tol = 1.e-3_real64
    ! Tolerance on part R,Z contour comparison
    Real(real64), Parameter :: check_AS_tol   = 1.e-8_real64

    !- End of header -------------------------------------------------------------

    ! Read parts list file and query each part for dimensions and file type
    Open(iu_plist,file=fname_plist,status='old',iostat=iostat)
    If (iostat /= 0) Then
       Write(*,*) "Error: Unable to open file: ",Trim(fname_plist)
       Call fin_mpi(.true.) ! True means this is an exit-on-error
    End If
    Read(iu_plist,*) nparts
    Allocate(nt_parts(nparts),np_parts(nparts))
    Allocate(part_names(nparts),part_type(nparts))

    ! Determine file format from extension
    Do ipart = 1,nparts
       Read(iu_plist,*) part_name
       part_names(ipart) = part_name
       name_len = Len_trim(part_name)
       If ( part_name(name_len-8:name_len) == '.2d.jpart' ) Then
          part_type(ipart) = 1
       Else If ( part_name(name_len-3:name_len) == '.txt' ) Then
          part_type(ipart) = 2
       Else
          part_type(ipart) = 0
       End If
       ! Set up quads for triangulation
       If (part_type(ipart) .le. 1) Then
          ! For "part" or "jpart" files
          Call query_part(part_name,ntor,npol,nfp_part)
       Else
          ! For pre-triangulated files (.txt)
          ntor = 0
          npol = 0
          nfp_part = 1
       End If
       nt_parts(ipart) = ntor
       np_parts(ipart) = npol

       If (verbose) Then
          Write(*,*) 'Part',ipart,Trim(Adjustl(part_names(ipart)))
       End If
    End Do
    Close(iu_plist)

    ! Allocate arrays for part coordinates
    nt_max = Maxval(nt_parts)
    np_max = Maxval(np_parts)
    Allocate(Rparts(nparts,nt_max,np_max))
    Allocate(Zparts(nparts,nt_max,np_max))
    Allocate(Pparts(nparts,nt_max,np_max))
    Allocate(Pmins(nparts),Pmaxs(nparts))
    Allocate(is_AS_part(nparts),force_non_AS(nparts))

    Rparts = 0._real64
    Zparts = 0._real64
    Pparts = 0._real64
    Pmins  = 0._real64
    Pmaxs  = 0._real64
    is_AS_part = .false.

    ! Part coordinates are written to file
    Open(iu_parts,file=fname_parts)
    Write(iu_parts,*) nparts,nt_max,np_max

    ! Load part coordinates
    Do ipart = 1,nparts
       ntor = nt_parts(ipart)
       npol = np_parts(ipart)

       Allocate( Ppart(ntor,npol),Rpart(ntor,npol),Zpart(ntor,npol) )
       ! Read part
       If (part_type(ipart) .EQ. 0) Then
          Call load_w7_part(part_names(ipart),label,ntor,npol,msym,Rpart,Zpart,Ppart,force_non_AS(ipart))
       Else If (part_type(ipart) .EQ. 1) Then
          Call load_2d_jpart(part_names(ipart),label,ntor,npol,msym,Rpart,Zpart,Ppart,force_non_AS(ipart))
       Else If (part_type(ipart) .EQ. 2) Then
          ! For triangle parts here we just want to know how many triangles (lines) each has
          !          Call query_tri_part(part_names(ipart),ntri_parts(ipart))
       Else
          If (verbose) Write(*,*) 'Did not recognize part_type',part_type(ipart),'for part',ipart
       Endif

       ! Set R,Z arrays
       Rparts(ipart,1:ntor,1:npol) = Rpart
       Zparts(ipart,1:ntor,1:npol) = Zpart

       ! Shift phi to first field period and set Phi array
       Do i=1,ntor
          Do j=1,npol
             Call wrap_phi(Ppart(i,j),2._real64*pi/Real(msym,real64))
          Enddo
       Enddo
       Pparts(ipart,1:ntor,1:npol) = Ppart

       ! Get min/max of Phi for filtering out intersection checks and to determine if part is AS
       If (part_type(ipart) .eq. 2) Then
          ! For tri parts this will get updated later
          Pmins(ipart) = 0.d0
          Pmaxs(ipart) = period
       Else
          Pmins(ipart) = Minval(Ppart)
          Pmaxs(ipart) = Maxval(Ppart)
       End If

       ! Extra screen output
       If (verbose) Then
          If (part_type(ipart) .eq. 2) Then
             Write(*,*)'Part ',ipart,'is a triangle part defined with nsym',nfp_part
             force_non_AS(ipart) = .true.
          Else
             Write(*,*)'Part ',ipart,'[nt,np] =',ntor,npol,'defined with nsym',nfp_part
             Write(*,*)'  Extends from Phi = ',Pmins(ipart)*180./pi,' to ',Pmaxs(ipart)*180./pi,' deg.'
             Write(*,*)'  Has part_type',part_type(ipart)
          End If
       End If

       If (Pmaxs(ipart) .gt. period) Then
          Write(*,*) 'Warning: Part extends beyond Bfield period, extra range is not used!'
       End If

       ! Part is AS if the points are the same across each phi cut
       ! AND the phi range is equal to the Bfield period
       check_AS = 0._real64
       Do i=2,ntor
          check_AS = Max(Maxval(Abs(Rpart(i,1:npol) - Rpart(1,1:npol))) &
               + Maxval(Abs(Zpart(i,1:npol) - Zpart(1,1:npol))),check_AS)
       Enddo
       If (       (check_AS .lt. check_AS_tol) &
                                !       .and. (Pmins(ipart)          .le. phi_period_tol) &
                                !       .and. (period - Pmaxs(ipart) .le. phi_period_tol) &
            ) Then
          If (force_non_AS(ipart)) Then
             If (part_type(ipart) .eq. 2) Then
                If (verbose) Write(*,*) '  Triangle parts are always treated as non-axistymmetric'
             Else
                If (verbose) Write(*,*) '  Part appears AS but AS treatment overridden in part file'
             End If
          Else
             If (verbose) Write(*,*) '  Part will be treated as axisymmetric over this phi range!'
             is_AS_part(ipart) = .true.
          End If
       Endif

       ! Clean up
       Deallocate(Rpart,Zpart,Ppart)

       ! Write all parts file
       Write(iu_parts,*) ntor,npol
       Do i = 1,nt_max
          Do j = 1,np_max
             Write(iu_parts,*) Rparts(ipart,i,j),Zparts(ipart,i,j),Pparts(ipart,i,j)
          Enddo
       Enddo
    Enddo

    Close(iu_parts) !all parts file

    ! Load vessel
    If (verbose) Write(6,*) 'Loading vessel file: ',Trim(Adjustl(fname_ves))
    Call query_part(fname_ves,ntor_ves,npol_ves,msym_ves)
    Allocate(R_ves(ntor_ves,npol_ves),Z_ves(ntor_ves,npol_ves),P_ves(ntor_ves,npol_ves))
    Call load_w7_part(fname_ves,label,ntor_ves,npol_ves,msym_ves,R_ves,Z_ves,P_ves,force_non_AS_ves)

    ! Shift phi to first field period and set Phi array
    Do i=1,ntor_ves
       Do j=1,npol_ves
          Call wrap_phi(P_ves(i,j),2._real64*pi/Real(msym_ves,real64))
       Enddo
    Enddo

    ! Check if vessel is AS
    is_AS_ves = .false.
    check_AS = 0._real64
    Do i=2,ntor_ves
       check_AS = Max(Maxval(Abs(R_ves(i,1:npol_ves) - R_ves(1,1:npol_ves))) &
            + Maxval(Abs(Z_ves(i,1:npol_ves) - Z_ves(1,1:npol_ves))),check_AS)
    Enddo
    If (       (check_AS .lt. check_AS_tol) &
         .and. (Minval(P_ves)          .le. phi_period_tol) &
         .and. (period - Maxval(P_ves) .le. phi_period_tol) &
         ) Then
       If (force_non_AS_ves) Then
          If (verbose) Write(*,*) '  Vessel appears AS but AS treatment overridden in part file'
       Else
          If (verbose) Write(*,*) '  Vessel will be treated as axisymmetric!'
          is_AS_ves = .true.
       End If
    Endif

  End Subroutine read_parts
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

End Module read_parts_mod

!-----------------------------------------------------------------------------
!
!   Routines/modules related to g3d/efit
!
!     Module g3d_module
!       Subroutine bfield_geq_bicub
!       Subroutine readg_g3d
!       Subroutine get_psi_bicub_coeffs
!       Subroutine get_becub_mat
!       Subroutine get_psi_derivs_bicub
!       Subroutine close_gfile
!       Function psi_bi
!       Function dsdr_bi
!       Function dsdz_bi
!-----------------------------------------------------------------------------

Module g_typedef
  Use kind_mod, Only : int32, real64
  Implicit None
  Private
  Public :: get_g_bspl_ord
  Type, Public :: g_type
    Integer(int32) :: mw, mh, nbdry, limitr
    Character(Len=6) :: ecase
    Real(real64) :: xdim, zdim, rzero, rgrid1, zmid, rmaxis, zmaxis, &
         ssimag, ssibry, bcentr, cpasma, dr, dz, ip_sign
    Real(real64), Allocatable, Dimension(:) :: fpol, pres, ffprim, &
         pprime, qpsi, r, z, pn
    Real(real64), Allocatable, Dimension(:,:) :: bdry, lim, psirz
    Real(real64), Allocatable, Dimension(:) :: pnknot, fpol_bscoef    
    Real(real64), Allocatable, Dimension(:,:,:) :: bicub_coeffs
    ! Private
    Integer, Private :: bspl_ord = 3
  End Type g_type

Contains

  Function get_g_bspl_ord(g) Result(bspl_ord)
    Use kind_mod, Only : int32
    Implicit None
    Type(g_type), Intent(In) :: g
    Integer(int32) :: bspl_ord
    bspl_ord = g%bspl_ord
  End Function get_g_bspl_ord
End Module g_typedef

!-----------------------------------------------------------------------------
!+ Module containing routines for gfile psi, B interpolation
!-----------------------------------------------------------------------------
Module g3d_module
  Use kind_mod, Only : int32, real64
  Use g_typedef, Only : g_type, get_g_bspl_ord
  Implicit None
  Private
  Public :: g_type
  Public :: readg_g3d, close_gfile
  Public :: bfield_geq_bicub, get_psi_derivs_bicub
  Public :: get_psiN_bicub
Contains
  
  !-----------------------------------------------------------------------------
  !+ returns Bcyl from gfile bfield at R,Z
  !-----------------------------------------------------------------------------
  Subroutine bfield_geq_bicub(g,R1,Z1,Npts,Bout,ierr,verbose)
    ! Description: 
    !  
    ! Output:
    !   Bout = (:,[Br,Bz,Bt])
    !
    !         
    ! History:
    !  Version   Date      Comment
    !  -------   ----      -------
    !  1.0     04/14/2011  Ported from Canik's idl routines.  JDL
    ! 
    ! Author(s): J.D. Lore 04/14/2011 
    !
    ! Modules used:
    Use kind_mod, Only : int32, real64
    Use bspline, Only : dbsval
    Implicit None

    ! Input/output                      !See above for descriptions
    Type(g_type), Intent(In) :: g
    Integer(int32),Intent(in) :: Npts
    Real(real64),Dimension(Npts),Intent(in)  :: R1,Z1
    Real(real64),Dimension(Npts,3),Intent(out) :: Bout(Npts,3)
    Integer(int32),Intent(out) :: ierr
    Logical, Optional :: verbose
    ! Local Scalars
    Integer(int32) :: ir,iz,index,ii
    Real(real64) :: dir,diz
    Real(real64) :: psi1,dsdr1,dsdz1,br1,bz1,psiN
    Real(real64) :: fpol,bt1
    Logical :: myverbose 
    !- End of header -------------------------------------------------------------
    myverbose = .true.
    If (present(verbose)) Then
      myverbose = verbose
    Endif

    !If (.not. Allocated(g%r) ) Then
    !  Write(*,*) 'G VARIABLES NOT ALLOCATED, EXITING FROM bfield_geq_bicub!'
    !  Stop
    !Endif

    ierr = 0
    Do ii = 1,Npts 

      ! update for vec.
      ir = Floor( (R1(ii)-g%r(1))/g%dr ) + 1
      iz = Floor( (Z1(ii)-g%z(1))/g%dz ) + 1

      ! Check for points off grid
      If ( (ir .le. 2) .or. (ir .ge. g%mw - 1) ) Then
        If (myverbose) &
             Write(*,'(3(a,f12.3),a)') 'bfield_geq: Point off grid in R: R = ',R1(ii), &
             '. [Rmin,Rmax] = [',g%r(3),',',g%r(g%mw)-1,']'
        ierr = 1
        Bout = 0.d0
        return
      Endif
      If ( (iz .le. 1) .or. (iz .ge. g%mh - 1) ) Then
        If (myverbose) &
             Write(*,'(3(a,f12.3),a)') 'bfield_geq: Point off grid in Z: Z = ',Z1(ii), &
             '. [Zmin,Zmax] = [',g%z(3),',',g%z(g%mh)-1,']'
        ierr = 1
        Bout = 0.d0
        return
      Endif

      dir = (R1(ii) - g%r(ir))/g%dr
      diz = (Z1(ii) - g%z(iz))/g%dz

      index = iz + g%mh*(ir-1)
      psi1 = psi_bi(g,index,dir,diz)
      dsdr1 = dsdr_bi(g,index,dir,diz)
      dsdz1 = dsdz_bi(g,index,dir,diz)

      br1 = -dsdz1/R1(ii)
      bz1 =  dsdr1/R1(ii)

      psiN = (psi1*g%ip_sign - g%ssimag)/(g%ssibry-g%ssimag)

      ! Toroidal field
      If (psiN .ge. 0._real64 .AND. psiN .le. 1._real64) Then 
        fpol = dbsval(psiN,get_g_bspl_ord(g),g%pnknot,g%mw,g%fpol_bscoef)
        bt1 = fpol/R1(ii)
      Else
        bt1 = g%bcentr*g%rzero/R1(ii)
      Endif

      Bout(ii,1) = br1
      Bout(ii,2) = bz1
      Bout(ii,3) = bt1
    Enddo

  End Subroutine bfield_geq_bicub

  !-----------------------------------------------------------------------------
  !+ reads a g file
  !-----------------------------------------------------------------------------
  Subroutine readg_g3d(filename,g)
    !
    ! Description: 
    !  Abbreviated function to read gfiles.  Assumes formatted file
    !
    ! History:
    !  Version   Date      Comment
    !  -------   ----      -------
    !  1.0     04/12/2011  Ported from Canik's idl routines.  JDL
    !
    ! Author(s): J.D. Lore - 04/12/2011
    !
    ! Modules used:
    Use kind_mod, Only: real64, int32
    Use bspline
    Implicit None

    ! Input/output                      !See above for descriptions
    Character(len=*), Intent(In) :: filename
    Type(g_type), Intent(Out) :: g
    
    ! Local scalars
    Integer(int32) :: iocheck,idum,i,j
    Real(real64) :: xdum

    ! Local arrays (1D)
    Character(len=100) :: sjunk

    Logical :: debug = .false.
    ! Allocatable arrays

    !- End of header -------------------------------------------------------------

    Write(*,*) 'Reading gfile: ',trim(filename)

    Open(UNIT=99,FILE=filename,STATUS="old",IOSTAT=iocheck)
    If ( iocheck /= 0 ) Then
      Write(*,*) 'Error opening gfile: ', filename
      Stop 'Exiting: I/O Error in subroutine readg_g3d'
    Endif
    If (debug) Write(*,*) 'Debugging readg_g3d'
    Read(99,'(a8,a42,i3,2i4)') g%ecase,sjunk,idum,g%mw,g%mh
    If (debug) Then
      Write(*,*) 'g%ecase: ',g%ecase
      Write(*,*) 'g%mw,g%mh: ',g%mw,g%mh
    Endif

    Allocate(g%fpol(g%mw),g%pres(g%mw),g%ffprim(g%mw))
    Allocate(g%pprime(g%mw),g%psirz(g%mw,g%mh),g%qpsi(g%mw))

    read(99,'(5e16.9)') g%xdim,g%zdim,g%rzero,g%rgrid1,g%zmid
    If (debug) Then
      Write(*,*) 'g%xdim',g%xdim
      Write(*,*) 'g%zdim',g%zdim
      Write(*,*) 'g%rzero',g%rzero
      Write(*,*) 'g%rgrid1',g%rgrid1
      Write(*,*) 'g%zmid',g%zmid
    Endif
    read(99,'(5e16.9)') g%rmaxis,g%zmaxis,g%ssimag,g%ssibry,g%bcentr
    If (debug) Then
      Write(*,*) 'g%rmaxis',g%rmaxis
      Write(*,*) 'g%zmaxis',g%zmaxis
      Write(*,*) 'g%ssimag',g%ssimag
      Write(*,*) 'g%ssibry',g%ssibry
      Write(*,*) 'g%bcentr',g%bcentr
    Endif
    read(99,'(1e16.9)') g%cpasma
    If (debug) Write(*,*) 'g%cpasma',g%cpasma
    read(99,'(5e16.9)') xdum

    read(99,'(5e16.9)') (g%fpol(i),i=1,g%mw)
    If (debug) Write(*,*) 'g%fpol([1,g%mw])',g%fpol(1),g%fpol(g%mw)
    read(99,'(5e16.9)') (g%pres(i),i=1,g%mw)
    read(99,'(5e16.9)') (g%ffprim(i),i=1,g%mw)
    read(99,'(5e16.9)') (g%pprime(i),i=1,g%mw)

    read(99,'(5e16.9)') ((g%psirz(i,j),i=1,g%mw),j=1,g%mh)
    read(99,'(5e16.9)') (g%qpsi(i),i=1,g%mw)
    read(99,'(2i5)')    g%nbdry,g%limitr

    Allocate(g%bdry(2,g%nbdry))
    Allocate(g%lim(2,g%limitr))
    read(99,'(5e16.9)') ((g%bdry(i,j),i=1,2),j=1,g%nbdry)
    read(99,'(5e16.9)') ((g%lim(i,j),i=1,2),j=1,g%limitr)

    Close(99)



    !
    ! Postprocessing
    !

    g%dR = g%xdim/(g%mw-1)
    g%dZ = g%zdim/(g%mh-1)

    Allocate(g%r(g%mw),g%z(g%mh),g%pn(g%mw))

    Do i=0,g%mw-1
      g%r(i+1) = g%rgrid1 + g%dR*i
      g%pn(i+1) = Real(i,real64)/(g%mw-1)
    Enddo

    Do i=0,g%mh-1
      g%z(i+1) = g%zmid - 0.5_real64*g%zdim + g%dZ*i
    Enddo

    g%ip_sign = -g%cpasma/abs(g%cpasma)

    Allocate(g%bicub_coeffs(g%mw*g%mh,4,4))
    g%bicub_coeffs = get_psi_bicub_coeffs(g)

    ! B-Spline fit poloidal current function
    Allocate(g%pnknot(get_g_bspl_ord(g)+g%mw))
    Allocate(g%fpol_bscoef(g%mw))
    Call dbsnak(g%mw,g%pn,get_g_bspl_ord(g),g%pnknot)
    Call dbsint(g%mw,g%pn,g%fpol,get_g_bspl_ord(g),g%pnknot,g%fpol_bscoef)

  End Subroutine readg_g3d


  !-----------------------------------------------------------------------------
  !+ Returns array of coefficients for bicubic interpolation
  !-----------------------------------------------------------------------------
  Function get_psi_bicub_coeffs(g)  &
       Result(psi_bicub_coeffs)
    !
    ! Description: 
    !  
    ! Function arguments:
    ! Output:
    !
    !         
    ! History:
    !  Version   Date      Comment
    !  -------   ----      -------
    !  1.0     04/15/2011  Ported from Canik's idl routines.  JDL
    ! 
    ! Author(s): J.D. Lore 04/15/2011 
    !
    ! Modules used:
    Use kind_mod, Only: real64, int32

    Implicit None

    ! Input/output                      !See above for descriptions
    Type(g_type), Intent(In) :: g
    Real(real64),dimension(g%mw*g%mh,4,4)   :: psi_bicub_coeffs
    ! Local Scalars
    Integer(int32) :: nr,nz,ir,iz,index,inv_err
    ! Local arrays 
    Real(real64),dimension(g%mw,g%mh) :: psi2d,dsdr,dsdz,d2sdrdz
    Real(real64),dimension(16,16) :: bicub_mat,bicub_mat_fac
    Real(real64),dimension(16) :: b
    Integer(int32),dimension(16) :: ipiv

    Interface
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        CHARACTER          TRANS
        INTEGER            INFO, LDA, LDB, N, NRHS
        INTEGER            IPIV( * )
        DOUBLE PRECISION   A( lda, * ), B( ldb, * )
      End SUBROUTINE DGETRS
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
        INTEGER INFO, LDA, M, N
        INTEGER            IPIV( * )
        DOUBLE PRECISION   A( lda, * )
      End SUBROUTINE DGETRF
    End Interface

    ! Local parameters                 


    !- End of header -------------------------------------------------------------

    nr = g%mw
    nz = g%mh
    psi2d = g%ip_sign * g%psirz

    dsdr = (cshift(psi2d,shift=1,dim=1) - cshift(psi2d,shift=-1,dim=1))/(2._real64*g%dr)
    dsdz = (cshift(psi2d,shift=1,dim=2) - cshift(psi2d,shift=-1,dim=2))/(2._real64*g%dz)
    d2sdrdz = (cshift(cshift(psi2d,shift=1,dim=1),shift=1,dim=2) & 
         - cshift(cshift(psi2d,shift=-1,dim=1),shift=1,dim=2) &
         - cshift(cshift(psi2d,shift=1,dim=1),shift=-1,dim=2) &
         + cshift(cshift(psi2d,shift=-1,dim=1),shift=-1,dim=2))/(4._real64*g%dr*g%dz)

    bicub_mat = get_bicub_mat()

    ! Perform LU factorization
    bicub_mat_fac = bicub_mat
    Call DGETRF(16,16,bicub_mat_fac,16,ipiv,inv_err)
    If (inv_err .ne. 0) Then
      Stop "Error from DGETRF factorizing bicub_mat in get_psi_bicub_coeffs"
    Endif

    psi_bicub_coeffs = 0._real64
    Do ir = 1,nr-1
      Do iz = 1,nz-1
        b = (/psi2d(ir,iz),            psi2d(ir+1,iz),            psi2d(ir,iz+1),            psi2d(ir+1,iz+1),     &
             dsdr(ir,iz)*g%dr,        dsdr(ir+1,iz)*g%dr,        dsdr(ir,iz+1)*g%dr,        dsdr(ir+1,iz+1)*g%dr, &
             dsdz(ir,iz)*g%dz,        dsdz(ir+1,iz)*g%dz,        dsdz(ir,iz+1)*g%dz,        dsdz(ir+1,iz+1)*g%dz, &
             d2sdrdz(ir,iz)*g%dr*g%dz,d2sdrdz(ir+1,iz)*g%dr*g%dz,d2sdrdz(ir,iz+1)*g%dr*g%dz,d2sdrdz(ir+1,iz+1)*g%dr*g%dz/)
        ! Solve Ax=B
        Call DGETRS('N',16,1,bicub_mat_fac,16,ipiv,b,16,inv_err)
        If (inv_err .ne. 0) Then
          Stop "Error from DGETRS solving AX=B for X in get_psi_bicub_coeffs"
        Endif
        index = iz + nz*(ir-1)
        psi_bicub_coeffs(index,1:4,1) = b(1:4)
        psi_bicub_coeffs(index,1:4,2) = b(5:8)
        psi_bicub_coeffs(index,1:4,3) = b(9:12)
        psi_bicub_coeffs(index,1:4,4) = b(13:16)
      Enddo
    Enddo

  EndFunction get_psi_bicub_coeffs

  !-----------------------------------------------------------------------------
  !+ Returns bicubic interpolation matrix
  !-----------------------------------------------------------------------------
  Function get_bicub_mat()  &
       Result(bicub_mat)
    !
    ! Description: 
    !  
    ! Function arguments:
    ! Output:
    !
    !         
    ! History:
    !  Version   Date      Comment
    !  -------   ----      -------
    !  1.0     04/19/2011  Ported from Canik's idl routines.  JDL
    ! 
    ! Author(s): J.D. Lore 04/19/2011 
    !
    ! Modules used:
    Use kind_mod, Only: real64

    Implicit None

    ! Input/output                      !See above for descriptions
    Real(real64),dimension(16,16)        :: bicub_mat
    
    ! Local Parameters
    Real(real64), Parameter :: One   = 1._real64, &
         Two   = 2._real64, &
         Three = 3._real64, &
         Four  = 4._real64, &
         Six   = 6._real64, &
         Nine  = 9._real64

    !- End of header -------------------------------------------------------------

    bicub_mat(:,:) = 0._real64

    ! Function values at corners
    bicub_mat(1,1) = One
    bicub_mat(2,(/1,2,3,4/)) = (/One,One,One,One/)
    bicub_mat(3,(/1,5,9,13/)) = (/One,One,One,One/)
    bicub_mat(4,:) = One


    ! 1st derivatives at corners: x direction
    bicub_mat(5,2) = One
    bicub_mat(6,(/2,3,4/)) = (/One,Two,Three/)
    bicub_mat(7,(/2,6,10,14/)) = (/One,One,One,One/)
    bicub_mat(8,(/2,3,4/)) = (/One,Two,Three/)
    bicub_mat(8,(/6,7,8/)) = (/One,Two,Three/)
    bicub_mat(8,(/10,11,12/)) = (/One,Two,Three/)
    bicub_mat(8,(/14,15,16/)) = (/One,Two,Three/)

    ! 1st derivatives at corners: y direction
    bicub_mat(9,5) = One
    bicub_mat(10,(/5,6,7,8/)) = (/One,One,One,One/)
    bicub_mat(11,(/5,9,13/)) = (/One,Two,Three/)
    bicub_mat(12,(/5,9,13/)) = (/One,Two,Three/)
    bicub_mat(12,(/6,10,14/)) = (/One,Two,Three/)
    bicub_mat(12,(/7,11,15/)) = (/One,Two,Three/)
    bicub_mat(12,(/8,12,16/)) = (/One,Two,Three/)

    ! cross derivatives at corners
    bicub_mat(13,6) = One
    bicub_mat(14,(/6,7,8/)) = (/One,Two,Three/)
    bicub_mat(15,(/6,10,14/)) = (/One,Two,Three/)
    bicub_mat(16,(/6,10,14/)) = (/One,Two,Three/)
    bicub_mat(16,(/7,11,15/)) = (/Two,Four,Six/)
    bicub_mat(16,(/8,12,16/)) =(/Three,Six,Nine/)

  EndFunction get_bicub_mat


  !-----------------------------------------------------------------------------
  !+ returns psi and derivatives from gfile bfield at R,Z
  !-----------------------------------------------------------------------------
  Subroutine get_psi_derivs_bicub(g,R1,Z1,Npts,psiout,dpsidr,dpsidz,ierr)
    Use kind_mod, Only: real64, int32
    Implicit None

    ! Input/output                      !See above for descriptions
    Type(g_type), Intent(In) :: g
    Integer(int32),Intent(in) :: Npts
    Real(real64),Dimension(Npts),Intent(in)  :: R1,Z1
    Real(real64),Dimension(Npts),Intent(out) :: psiout, dpsidr, dpsidz
    Integer(int32),Intent(out) :: ierr

    ! Local Scalars
    Integer(int32) :: ir,iz,index,ii
    Real(real64) :: dir,diz

    !- End of header -------------------------------------------------------------

    If (.not. Allocated(g%r) ) Then
      Write(*,*) 'G VARIABLES NOT ALLOCATED, EXITING FROM get_psi_derivs_bicub!'
      Stop
    Endif

    ierr = 0
    Do ii = 1,Npts 

      ! update for vec.
      ir = Floor( (R1(ii)-g%r(1))/g%dr ) + 1
      iz = Floor( (Z1(ii)-g%z(1))/g%dz ) + 1

      ! Check for points off grid
      If ( (ir .le. 2) .or. (ir .ge. g%mw - 1) ) Then
        Write(*,'(3(a,f12.3),a)') 'psi deriv: Point off grid in R: R = ',R1(ii),&
             '. [Rmin,Rmax] = [',g%r(1),',',g%r(g%mw),']'
        ierr = 1
        psiout = 0.d0
        return
      Endif
      If ( (iz .le. 1) .or. (iz .ge. g%mh - 1) ) Then
        Write(*,'(3(a,f12.3),a)') 'psi deriv: Point off grid in Z: Z = ',Z1(ii),&
             '. [Zmin,Zmax] = [',g%z(1),',',g%z(g%mh),']'
        ierr = 1
        psiout = 0.d0
        return
      Endif

      dir = (R1(ii) - g%r(ir))/g%dr
      diz = (Z1(ii) - g%z(iz))/g%dz
      index = iz + g%mh*(ir-1)
      psiout(ii) = psi_bi(g,index,dir,diz)
      dpsidr(ii) = dsdr_bi(g,index,dir,diz)
      dpsidz(ii) = dsdz_bi(g,index,dir,diz)
    Enddo

  End Subroutine get_psi_derivs_bicub



  !-----------------------------------------------------------------------------
  !+ 
  !-----------------------------------------------------------------------------
  Function psi_bi(g,index,dir,diz)
    Use kind_mod, Only: real64, int32
    Implicit None
    Type(g_type), Intent(In) :: g
    Integer(int32), Intent(In) :: index
    Real(real64), Intent(In) :: dir, diz
    Real(real64) :: psi_bi
    psi_bi = g%bicub_coeffs(index,1,1)   + g%bicub_coeffs(index,2,1)*dir       + &
         g%bicub_coeffs(index,3,1)*dir*dir       + g%bicub_coeffs(index,4,1)*dir*dir*dir      + &
         g%bicub_coeffs(index,1,2)*diz   + g%bicub_coeffs(index,2,2)*dir*diz   + &
         g%bicub_coeffs(index,3,2)*dir*dir*diz   + g%bicub_coeffs(index,4,2)*dir*dir*dir*diz  + &
         g%bicub_coeffs(index,1,3)*diz*diz  + g%bicub_coeffs(index,2,3)*dir*diz*diz  + &
         g%bicub_coeffs(index,3,3)*dir*dir*diz*diz  + g%bicub_coeffs(index,4,3)*dir*dir*dir*diz*diz + &
         g%bicub_coeffs(index,1,4)*diz*diz*diz  + g%bicub_coeffs(index,2,4)*dir*diz*diz*diz  + &
         g%bicub_coeffs(index,3,4)*dir*dir*diz*diz*diz  + g%bicub_coeffs(index,4,4)*dir*dir*dir*diz*diz*diz
  End Function psi_bi


  !-----------------------------------------------------------------------------
  !+ 
  !-----------------------------------------------------------------------------
  Function dsdr_bi(g,index,dir,diz)
    Use kind_mod, Only: real64, int32
    Implicit None
    Type(g_type), Intent(In) :: g
    Integer(int32), Intent(In) :: index
    Real(real64), Intent(In) :: dir, diz
    Real(real64) :: dsdr_bi
    Real(real64), Parameter :: TWO = 2._real64, THREE = 3._real64  
    dsdr_bi = (g%bicub_coeffs(index,2,1)       + TWO*g%bicub_coeffs(index,3,1)*dir +    &
         THREE*g%bicub_coeffs(index,4,1)*dir*dir      + &
         g%bicub_coeffs(index,2,2)*diz   + TWO*g%bicub_coeffs(index,3,2)*dir*diz   + &
         THREE*g%bicub_coeffs(index,4,2)*dir*dir*diz  + &
         g%bicub_coeffs(index,2,3)*diz*diz  + TWO*g%bicub_coeffs(index,3,3)*dir*diz*diz  + &
         THREE*g%bicub_coeffs(index,4,3)*dir*dir*diz*diz + &
         g%bicub_coeffs(index,2,4)*diz*diz*diz  + TWO*g%bicub_coeffs(index,3,4)*dir*diz*diz*diz + &
         THREE*g%bicub_coeffs(index,4,4)*dir*dir*diz*diz*diz)/g%dr
  End Function dsdr_bi


  !-----------------------------------------------------------------------------
  !+ 
  !-----------------------------------------------------------------------------
  Function dsdz_bi(g,index,dir,diz)
    Use kind_mod, Only : int32, real64
    Implicit None
    Type(g_type), Intent(In) :: g
    Integer(int32), Intent(In) :: index
    Real(real64), Intent(In) :: dir, diz
    Real(real64) :: dsdz_bi
    Real(real64), Parameter :: TWO = 2._real64, THREE = 3._real64  
    dsdz_bi = (g%bicub_coeffs(index,1,2)                 + g%bicub_coeffs(index,2,2)*dir                 + &
               g%bicub_coeffs(index,3,2)*dir*dir         + g%bicub_coeffs(index,4,2)*dir*dir*dir         + &
           TWO*g%bicub_coeffs(index,1,3)*diz             + TWO*g%bicub_coeffs(index,2,3)*dir*diz         + &
           TWO*g%bicub_coeffs(index,3,3)*dir*dir*diz     + TWO*g%bicub_coeffs(index,4,3)*dir*dir*dir*diz + &
         THREE*g%bicub_coeffs(index,1,4)*diz*diz         + THREE*g%bicub_coeffs(index,2,4)*dir*diz*diz   + &
         THREE*g%bicub_coeffs(index,3,4)*dir*dir*diz*diz + THREE*g%bicub_coeffs(index,4,4)*dir*dir*dir*diz*diz)/g%dz
  End Function dsdz_bi

  Subroutine close_gfile(g)
    Implicit None
    Type(g_type), Intent(InOut) :: g
    Write(*,*) 'Deallocating gfile variables'
    Deallocate(g%fpol,g%pres,g%ffprim)
    Deallocate(g%pprime,g%psirz,g%qpsi)
    Deallocate(g%bdry)
    Deallocate(g%lim)
    Deallocate(g%r,g%z,g%pn)
    Deallocate(g%bicub_coeffs)
    Deallocate(g%pnknot)
    Deallocate(g%fpol_bscoef)
  End Subroutine close_gfile

  !-----------------------------------------------------------------------------
  !+ returns psiN from gfile bfield at R,Z
  !-----------------------------------------------------------------------------
  Subroutine get_psiN_bicub(g,R1,Z1,Npts,psiNout,ierr)
    Use kind_mod, Only: real64, int32
    Use bspline
    Implicit None

    ! Input/output                      !See above for descriptions
    Type(g_type), Intent(In) :: g
    Integer(int32),Intent(in) :: Npts
    Real(real64),Dimension(Npts),Intent(in)  :: R1,Z1
    Real(real64),Dimension(Npts),Intent(out) :: psiNout
    Integer(int32),Intent(out) :: ierr
    ! Local variables
    Real(real64),Dimension(Npts) :: psi
    Integer(int32) :: ierr_tmp

    ierr = 0
    Call get_psi_bicub(g,R1,Z1,Npts,psi,ierr_tmp)
    psiNout = (psi - g%ssimag)/(g%ssibry-g%ssimag)
  End Subroutine get_psiN_bicub

  !-----------------------------------------------------------------------------
  !+ returns psi from gfile bfield at R,Z
  !-----------------------------------------------------------------------------
  Subroutine get_psi_bicub(g,R1,Z1,Npts,psiout,ierr)
    Use kind_mod, Only: real64, int32
    Use bspline
    Implicit None

    ! Input/output                      !See above for descriptions
    Type(g_type), Intent(In) :: g
    Integer(int32),Intent(in) :: Npts
    Real(real64),Dimension(Npts),Intent(in)  :: R1,Z1
    Real(real64),Dimension(Npts),Intent(out) :: psiout
    Integer(int32),Intent(out) :: ierr

    ! Local Scalars
    Integer(int32) :: ir,iz,ii
    Real(real64) :: dir,diz

    !- End of header -------------------------------------------------------------

    If (.not. Allocated(g%r) ) Then
      Write(*,*) 'G VARIABLES NOT ALLOCATED, EXITING FROM get_psi_bicub!'
      Stop
    Endif

    ierr = 0
    Do ii = 1,Npts 

      ! update for vec.
      ir = Floor( (R1(ii)-g%r(1))/g%dr ) + 1
      iz = Floor( (Z1(ii)-g%z(1))/g%dz ) + 1

      ! Check for points off grid
      If ( (ir .le. 2) .or. (ir .ge. g%mw - 1) ) Then
        Write(*,'(3(a,f12.3),a)') 'psi_geq: Point off grid in R: R = ',R1(ii),&
             '. [Rmin,Rmax] = [',g%r(1),',',g%r(g%mw),']'
        ierr = 1
        psiout(ii:Npts) = 0.d0
        return
      Endif
      If ( (iz .le. 1) .or. (iz .ge. g%mh - 1) ) Then
        Write(*,'(3(a,f12.3),a)') 'psi_geq: Point off grid in Z: Z = ',Z1(ii),&
             '. [Zmin,Zmax] = [',g%z(1),',',g%z(g%mh),']'
        ierr = 1
        psiout(ii:Npts) = 0.d0
        return
      Endif

      dir = (R1(ii) - g%r(ir))/g%dr
      diz = (Z1(ii) - g%z(iz))/g%dz
      psiout(ii) = g%ip_sign*psi_bi(g,iz + g%mh*(ir-1),dir,diz)


    Enddo

  End Subroutine get_psi_bicub
  
  
End Module g3d_module




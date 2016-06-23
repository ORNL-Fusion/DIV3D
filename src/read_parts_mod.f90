!-----------------------------------------------------------------------------
!+ Contains routines used to read and interact with part files.
!-----------------------------------------------------------------------------
Module read_parts_mod
! Description: 
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/14/2011  
!  1.1     12/14/2011  Combined with load parts JDL
! Author(s): J.D. Lore - 07/14/2011 - xxx
!

Use kind_mod

Implicit None
Integer(iknd),Allocatable,Dimension(:) :: nt_parts,np_parts
Integer(iknd), Allocatable, Dimension(:) :: pol_dirs, tor_dirs, part_type
Character(len=300),Allocatable,Dimension(:) :: part_names
Real(rknd),Allocatable :: Rparts(:,:,:),Zparts(:,:,:),Pparts(:,:,:)
Real(rknd), Allocatable, Dimension(:) :: Pmaxs, Pmins
Real(rknd), Allocatable :: R_ves(:,:),Z_ves(:,:),P_ves(:,:)
Real(rknd), Allocatable, Dimension(:,:) :: xmid, ymid, zmid, dmid

Integer(iknd) :: ntor_ves, npol_ves, msym_ves, msym
Integer(iknd) :: nparts
Integer(iknd) :: nt_max, np_max

Integer(iknd) :: ntri_max, ic_near
Real(rknd),Allocatable,Dimension(:,:,:) :: xtri,ytri,ztri
Integer(iknd),Allocatable,Dimension(:,:) :: check_tri
Integer(iknd),Allocatable,Dimension(:) :: ntri_parts
Real(rknd), allocatable :: near_part(:), near_tri(:)
!- End of header -------------------------------------------------------------

Contains

!-----------------------------------------------------------------------------
!+ Makes triangles from 2d parts
!-----------------------------------------------------------------------------
Subroutine make_triangles(fname_ptri,fname_ptri_mid)
Use kind_mod
Use io_unit_spec, Only : iu_ptri, iu_ptmid
Implicit None

Character(len=300), Intent(in) :: fname_ptri,fname_ptri_mid
Real(rknd) :: R1,Z1,P1,Pt1(3)
Real(rknd) :: R2,Z2,P2,Pt2(3)
Real(rknd) :: R3,Z3,P3,Pt3(3)
Real(rknd) :: R4,Z4,P4,Pt4(3), dmids(3)
Integer(iknd) :: ipart, jpol, itor, itri, npol, ntor
!- End of header -------------------------------------------------------------

!first have to get max number of triangles for allocation
Allocate(ntri_parts(nparts))
Do ipart=1,nparts
  ntor = nt_parts(ipart)
  npol = np_parts(ipart)
  ntri_parts(ipart) = (ntor-1)*(npol-1)*2_iknd
Enddo
ntri_max = Maxval(ntri_parts)

Allocate(xtri(nparts,ntri_max,3))
Allocate(ytri(nparts,ntri_max,3))
Allocate(ztri(nparts,ntri_max,3))
Allocate(xmid(nparts,ntri_max))
Allocate(ymid(nparts,ntri_max))
Allocate(zmid(nparts,ntri_max))
Allocate(dmid(nparts,ntri_max))
Allocate(check_tri(nparts,ntri_max))

xmid(:,:) = 0._rknd
ymid(:,:) = 0._rknd
zmid(:,:) = 0._rknd
xtri(:,:,:) = 0._rknd
ytri(:,:,:) = 0._rknd
ztri(:,:,:) = 0._rknd
check_tri(:,:) = 0_iknd

Open(iu_ptri,file=fname_ptri)
Write(iu_ptri,*) nparts
Open(iu_ptmid,file=fname_ptri_mid)
Write(iu_ptmid,*) nparts

Do ipart=1,nparts
  write(iu_ptri,*) ipart,ntri_parts(ipart)
  write(iu_ptmid,*) ipart,ntri_parts(ipart)
  ntor = nt_parts(ipart)
  npol = np_parts(ipart)
  itri = 1_iknd
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
      xtri(ipart,itri,1:3) = [Pt1(1),Pt2(1),Pt3(1)]
      ytri(ipart,itri,1:3) = [Pt1(2),Pt2(2),Pt3(2)]
      ztri(ipart,itri,1:3) = [Pt1(3),Pt2(3),Pt3(3)]
      check_tri(ipart,itri) = 1  !all 2d triangles are checked
      write(iu_ptri,*) itri
      write(iu_ptri,*) xtri(ipart,itri,1),ytri(ipart,itri,1),ztri(ipart,itri,1)
      write(iu_ptri,*) xtri(ipart,itri,2),ytri(ipart,itri,2),ztri(ipart,itri,2)
      write(iu_ptri,*) xtri(ipart,itri,3),ytri(ipart,itri,3),ztri(ipart,itri,3)
      ! Define triangle midpoint
      xmid(ipart,itri) = sum(xtri(ipart,itri,1:3))/3._rknd
      ymid(ipart,itri) = sum(ytri(ipart,itri,1:3))/3._rknd
      zmid(ipart,itri) = sum(ztri(ipart,itri,1:3))/3._rknd
      dmids(1) = sqrt( (xtri(ipart,itri,1) - xmid(ipart,itri))**2 + &
           (ytri(ipart,itri,1) - ymid(ipart,itri))**2 + (ztri(ipart,itri,1) - zmid(ipart,itri))**2 )
      dmids(2) = sqrt( (xtri(ipart,itri,2) - xmid(ipart,itri))**2 + &
           (ytri(ipart,itri,2) - ymid(ipart,itri))**2 + (ztri(ipart,itri,2) - zmid(ipart,itri))**2 )
      dmids(3) = sqrt( (xtri(ipart,itri,3) - xmid(ipart,itri))**2 + &
           (ytri(ipart,itri,3) - ymid(ipart,itri))**2 + (ztri(ipart,itri,3) - zmid(ipart,itri))**2 )
      dmid(ipart,itri) = maxval(dmids,1)
      write(iu_ptmid,*) itri,xmid(ipart,itri),ymid(ipart,itri),zmid(ipart,itri),dmid(ipart,itri)

      ! 2nd triangle
      itri = itri+1_iknd
      xtri(ipart,itri,1:3) = [Pt4(1),Pt2(1),Pt3(1)]
      ytri(ipart,itri,1:3) = [Pt4(2),Pt2(2),Pt3(2)]
      ztri(ipart,itri,1:3) = [Pt4(3),Pt2(3),Pt3(3)]
      check_tri(ipart,itri) = 1  !all 2d triangles are checked
      write(iu_ptri,*) itri
      write(iu_ptri,*) xtri(ipart,itri,1),ytri(ipart,itri,1),ztri(ipart,itri,1)
      write(iu_ptri,*) xtri(ipart,itri,2),ytri(ipart,itri,2),ztri(ipart,itri,2)
      write(iu_ptri,*) xtri(ipart,itri,3),ytri(ipart,itri,3),ztri(ipart,itri,3)
      ! Define triangle midpoint
      xmid(ipart,itri) = sum(xtri(ipart,itri,1:3))/3._rknd
      ymid(ipart,itri) = sum(ytri(ipart,itri,1:3))/3._rknd
      zmid(ipart,itri) = sum(ztri(ipart,itri,1:3))/3._rknd
      dmids(1) = sqrt( (xtri(ipart,itri,1) - xmid(ipart,itri))**2 + &
           (ytri(ipart,itri,1) - ymid(ipart,itri))**2 + (ztri(ipart,itri,1) - zmid(ipart,itri))**2 )
      dmids(2) = sqrt( (xtri(ipart,itri,2) - xmid(ipart,itri))**2 + &
           (ytri(ipart,itri,2) - ymid(ipart,itri))**2 + (ztri(ipart,itri,2) - zmid(ipart,itri))**2 )
      dmids(3) = sqrt( (xtri(ipart,itri,3) - xmid(ipart,itri))**2 + &
           (ytri(ipart,itri,3) - ymid(ipart,itri))**2 + (ztri(ipart,itri,3) - zmid(ipart,itri))**2 )
      dmid(ipart,itri) = maxval(dmids,1)
      write(iu_ptmid,*) itri,xmid(ipart,itri),ymid(ipart,itri),zmid(ipart,itri),dmid(ipart,itri)

      itri = itri+1_iknd

    Enddo
  Enddo
Enddo
Close(iu_ptri)
Close(iu_ptmid)

End Subroutine make_triangles
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!+ Reads 2d part file (w7 format)
!-----------------------------------------------------------------------------
Subroutine load_w7_part(fname,label,ntor,npol,msym,Rpart,Zpart,Phipart)
Use kind_mod
Use io_unit_spec, Only : iu_thispart
Implicit none
! Dummy variables
Character(len=100), Intent(in) :: fname
Integer(iknd), Intent(in) :: ntor, npol
Character(len=100), Intent(out) :: label
Integer(iknd),Intent(out) :: msym  
Real(rknd),Dimension(ntor,npol), Intent(out) :: &
  Rpart, Zpart, Phipart
Integer(iknd) :: itor, ipol, ntor_dum, npol_dum
! Local variables
Real(rknd) :: rshift, zshift, Phitmp
! Local parameters
Real(rknd), Parameter :: pi = 3.141592653589793238462643383279502_rknd
!- End of header -------------------------------------------------------------

open(iu_thispart,file=fname)
Read(iu_thispart,*) label
Read(iu_thispart,*) ntor_dum, npol_dum, msym, rshift, zshift

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
Rpart = (Rpart+rshift)*0.01_rknd
Zpart = (Zpart+zshift)*0.01_rknd
Phipart = Phipart*pi/180._rknd

Endsubroutine load_w7_part
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------





!-----------------------------------------------------------------------------
!+ Reads 2d jpart file 
!-----------------------------------------------------------------------------
Subroutine load_2d_jpart(fname,label,ntor,npol,msym,Rpart,Zpart,Phipart)
Use kind_mod
Use io_unit_spec, Only : iu_thispart
Implicit none
! Dummy variables
Character(len=100), Intent(in) :: fname
Integer(iknd), Intent(in) :: ntor,npol
Character(len=100), Intent(out) :: label
Integer(iknd),Intent(out) :: msym
Real(rknd),Dimension(ntor,npol), Intent(out) :: &
  Rpart, Zpart, Phipart
Integer(iknd) :: itor, ipol, ntor_dum, npol_dum
! Local variables
Real(rknd) :: rshift, zshift
! Local parameters
Real(rknd), Parameter :: pi = 3.141592653589793238462643383279502_rknd
!- End of header -------------------------------------------------------------

open(iu_thispart,file=fname)
Read(iu_thispart,*) label
Read(iu_thispart,*) ntor_dum, npol_dum, msym, rshift, zshift

! When read [R,Z] = cm and Phi = degrees
Do itor = 1,ntor 
  Do ipol = 1,npol
    Read(iu_thispart,*) Rpart(itor,ipol), Zpart(itor,ipol), Phipart(itor,ipol)
  Enddo
Enddo
close(iu_thispart)

! Convert to meters and radians
Rpart = (Rpart+rshift)*0.01_rknd
Zpart = (Zpart+zshift)*0.01_rknd
Phipart = Phipart*pi/180._rknd

Endsubroutine load_2d_jpart
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------




!-----------------------------------------------------------------------------
!+ Get dimensions from part file
!-----------------------------------------------------------------------------
Subroutine query_part(fname,ntor,npol)
Use kind_mod
Use io_unit_spec, Only : iu_thispart
Implicit none
Character(len=100), Intent(in) :: fname
Integer(iknd), Intent(out) :: ntor,npol
Integer(iknd) :: nfp
Real(rknd) :: rshift, zshift
Character(len=100) :: label
!- End of header -------------------------------------------------------------
Open(iu_thispart,file=fname)
Read(iu_thispart,*) label
Read(iu_thispart,*) ntor, npol, nfp, rshift, zshift
Close(iu_thispart)
Endsubroutine
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------




!-----------------------------------------------------------------------------
!+ Makes triangles in RZ plane from 2d part
!-----------------------------------------------------------------------------
Subroutine close_2pt_part(Rin,Zin,ntor,npol,dir,step,Rout,Zout)
Use kind_mod
Implicit none
Integer(iknd), Intent(in) :: ntor, npol, dir
Real(rknd), Intent(in) :: step
Real(rknd), Dimension(ntor,npol), Intent(in) :: Rin,Zin
Real(rknd), Dimension(ntor,npol+1), Intent(out) :: Rout, Zout
Integer(iknd) :: i
Real(rknd) :: R1,Z1,R2,Z2,dL,Rmid,Zmid
Real(rknd), Dimension(2) ::  V1,unit_out
!- End of header -------------------------------------------------------------

Rout(1:ntor,1:npol) = Rin
Zout(1:ntor,1:npol) = Zin

Do i=1,ntor
  R1 = Rin(i,1)
  Z1 = Zin(i,1)

  R2 = Rin(i,npol)
  Z2 = Zin(i,npol)
  dL = Dsqrt( (R1-R2)*(R1-R2) + (Z1-Z2)*(Z1-Z2) )
  
  V1(1) = R2-R1
  V1(2) = Z2-Z1
  V1 = V1/dL

  Rmid = R1 + V1(1)*dL/2._rknd
  Zmid = Z1 + V1(2)*dL/2._rknd

  unit_out(1) = -V1(2)
  unit_out(2) =  V1(1)
  unit_out = unit_out/Dsqrt( V1(1)*V1(1) + V1(2)*V1(2) )

  Rout(i,npol+1) = Rmid + real(dir,rknd)*unit_out(1)*step
  Zout(i,npol+1) = Zmid + real(dir,rknd)*unit_out(2)*step
!  Rout(i,npol+1) = Rmid 
!  Zout(i,npol+1) = Zmid 

!  Write(*,*) '------'
!  Write(*,*) Rout(i,:)
!  Write(*,*) Zout(i,:)
!  Write(*,*) '------'

Enddo
Endsubroutine close_2pt_part
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------




!-----------------------------------------------------------------------------
!+ Reads '3d' .jpart files
!-----------------------------------------------------------------------------
Subroutine read_3d_part()
!
! Description: 
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     12/07/2011  
! Author(s): J.D. Lore - 12/07/2011 - xxx
!
! Modules used:
Use kind_mod

Implicit none
!- End of header -------------------------------------------------------------
End Subroutine read_3d_part
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!+ Reads parts list and load part coordinates
!-----------------------------------------------------------------------------
Subroutine read_parts(fname_plist,fname_parts,fname_ves,verbose)
!
! Description: 
!
! History:
!  Version   Date      Comment
!  -------   ----      -------
!  1.0     07/14/2011  
!  1.1     02/14/2012  Updated to allow for general 2d parts and removed part directions
!
!  Notes:
!    part_type == 0 indicates 'w7' type. (R,Z) points at toroidal angles
!    part_type == 1 indicates '2d.jpart' type. (R,Z,phi) points
!  
!
!
! Author(s): J.D. Lore - 07/14/2011 - xxx
!
! Modules used:
Use kind_mod
Use io_unit_spec, Only: &
iu_plist, &    ! Parts filename list file (input) 
iu_parts
Implicit none
Character(len=300), Intent(in) ::  fname_plist, fname_parts, fname_ves
Logical, Intent(in) :: verbose
Real(rknd),Allocatable :: Rpart(:,:),Zpart(:,:),Ppart(:,:)
Integer(iknd) :: i, j
Integer(iknd) :: ipart
Integer(iknd) :: ntor, npol
Character(len=300) :: part_name
Character(len=300) :: label
! Local Parameters
Real(rknd), Parameter :: pi = 3.141592653589793238462643383279502_rknd
!- End of header -------------------------------------------------------------

! Read parts list file and query each part for dimensions
Open(iu_plist,file=fname_plist)
Read(iu_plist,*) nparts
Allocate(nt_parts(nparts),np_parts(nparts))
Allocate(part_names(nparts),part_type(nparts))
Do ipart = 1,nparts
  Read(iu_plist,*) part_name
!  WRite(6,*) 'fix this!'
  If ( Verify('.2d.jpart',part_name(10:len(part_name))) .EQ. 0 ) Then
    part_type(ipart) = 1
  Else
    part_type(ipart)= 0
  Endif
  Call query_part(part_name,ntor,npol)
  part_names(ipart) = part_name
  nt_parts(ipart) = ntor
  np_parts(ipart) = npol  
  If (verbose) Write(6,'(A,A,A,I4,I4)') ' ',Trim(Adjustl(part_names(ipart))),' [nt,np] =',ntor,npol
Enddo
Close(iu_plist)

! Allocate arrays for part coordinates
nt_max = Maxval(nt_parts)
np_max = Maxval(np_parts)
Allocate(Rparts(nparts,nt_max,np_max))
Allocate(Zparts(nparts,nt_max,np_max))
Allocate(Pparts(nparts,nt_max,np_max))
Allocate(Pmins(nparts),Pmaxs(nparts))

! Part coordinates are written to file
Open(iu_parts,file=fname_parts)
Write(iu_parts,*) nparts,nt_max,np_max

! Load part coordinates
Do ipart = 1,nparts
  ntor = nt_parts(ipart)
  npol = np_parts(ipart)
  Allocate( Ppart(ntor,npol),Rpart(ntor,npol),Zpart(ntor,npol) )

  If (part_type(ipart) .EQ. 0) Then
    Call load_w7_part(part_names(ipart),label,ntor,npol,msym,Rpart,Zpart,Ppart)
  Else
    Call load_2d_jpart(part_names(ipart),label,ntor,npol,msym,Rpart,Zpart,Ppart)
  Endif

  Rparts(ipart,1:ntor,1:npol) = Rpart
  Zparts(ipart,1:ntor,1:npol) = Zpart

  ! Shift phi to first field period
  Do i=1,ntor
    Do j=1,npol
      Do While (Ppart(i,j) .lt. 0.d0) 
        Ppart(i,j) = Ppart(i,j) + 2._rknd*pi/Real(msym,rknd)
      Enddo
      Ppart(i,j) = Mod(Ppart(i,j),2._rknd*pi/Real(msym,rknd))
    Enddo
  Enddo
  Pparts(ipart,1:ntor,1:npol) = Ppart

  Pmins(ipart) = Minval(Ppart)
  Pmaxs(ipart) = Maxval(Ppart)

  Deallocate(Rpart,Zpart,Ppart)

  Write(iu_parts,*) ntor,npol
  Do i = 1,nt_max
    Do j = 1,np_max
      Write(iu_parts,*) Rparts(ipart,i,j),Zparts(ipart,i,j),Pparts(ipart,i,j)
    Enddo
  Enddo

Enddo

Deallocate(part_names,part_type)
Close(iu_parts) !all parts file

! Load vessel
If (verbose) Write(6,*) 'Loading vessel file: ',Trim(Adjustl(fname_ves))
Call query_part(fname_ves,ntor_ves,npol_ves)
Allocate(R_ves(ntor_ves,npol_ves),Z_ves(ntor_ves,npol_ves),P_ves(ntor_ves,npol_ves))
Call load_w7_part(fname_ves,label,ntor_ves,npol_ves,msym_ves,R_ves,Z_ves,P_ves)

End Subroutine read_parts
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

End Module read_parts_mod


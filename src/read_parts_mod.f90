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
Integer(int32),Allocatable,Dimension(:) :: nt_parts,np_parts
Integer(int32), Allocatable, Dimension(:) :: pol_dirs, tor_dirs, part_type
Character(len=300),Allocatable,Dimension(:) :: part_names
Real(real64),Allocatable :: Rparts(:,:,:),Zparts(:,:,:),Pparts(:,:,:)
Real(real64), Allocatable, Dimension(:) :: Pmaxs, Pmins
Real(real64), Allocatable :: R_ves(:,:),Z_ves(:,:),P_ves(:,:)
Real(real64), Allocatable, Dimension(:,:) :: xmid, ymid, zmid, dmid

Integer(int32) :: ntor_ves, npol_ves, msym_ves, msym
Integer(int32) :: nparts
Integer(int32) :: nt_max, np_max

Integer(int32) :: ntri_max, ic_near
Real(real64),Allocatable,Dimension(:,:,:) :: xtri,ytri,ztri
Integer(int32),Allocatable,Dimension(:,:) :: check_tri
Integer(int32),Allocatable,Dimension(:) :: ntri_parts
Real(real64), allocatable :: near_part(:), near_tri(:)
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
Real(real64) :: R1,Z1,P1,Pt1(3)
Real(real64) :: R2,Z2,P2,Pt2(3)
Real(real64) :: R3,Z3,P3,Pt3(3)
Real(real64) :: R4,Z4,P4,Pt4(3), dmids(3)
Integer(int32) :: ipart, jpol, itor, itri, npol, ntor
!- End of header -------------------------------------------------------------

!first have to get max number of triangles for allocation
Allocate(ntri_parts(nparts))
Do ipart=1,nparts
  ntor = nt_parts(ipart)
  npol = np_parts(ipart)
  ntri_parts(ipart) = (ntor-1)*(npol-1)*2_int32
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

xmid(:,:) = 0._real64
ymid(:,:) = 0._real64
zmid(:,:) = 0._real64
xtri(:,:,:) = 0._real64
ytri(:,:,:) = 0._real64
ztri(:,:,:) = 0._real64
check_tri(:,:) = 0_int32

Open(iu_ptri,file=fname_ptri)
Write(iu_ptri,*) nparts
Open(iu_ptmid,file=fname_ptri_mid)
Write(iu_ptmid,*) nparts

Do ipart=1,nparts
  write(iu_ptri,*) ipart,ntri_parts(ipart)
  write(iu_ptmid,*) ipart,ntri_parts(ipart)
  ntor = nt_parts(ipart)
  npol = np_parts(ipart)
  itri = 1_int32
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
      xmid(ipart,itri) = sum(xtri(ipart,itri,1:3))/3._real64
      ymid(ipart,itri) = sum(ytri(ipart,itri,1:3))/3._real64
      zmid(ipart,itri) = sum(ztri(ipart,itri,1:3))/3._real64
      dmids(1) = sqrt( (xtri(ipart,itri,1) - xmid(ipart,itri))**2 + &
           (ytri(ipart,itri,1) - ymid(ipart,itri))**2 + (ztri(ipart,itri,1) - zmid(ipart,itri))**2 )
      dmids(2) = sqrt( (xtri(ipart,itri,2) - xmid(ipart,itri))**2 + &
           (ytri(ipart,itri,2) - ymid(ipart,itri))**2 + (ztri(ipart,itri,2) - zmid(ipart,itri))**2 )
      dmids(3) = sqrt( (xtri(ipart,itri,3) - xmid(ipart,itri))**2 + &
           (ytri(ipart,itri,3) - ymid(ipart,itri))**2 + (ztri(ipart,itri,3) - zmid(ipart,itri))**2 )
      dmid(ipart,itri) = maxval(dmids,1)
      write(iu_ptmid,*) itri,xmid(ipart,itri),ymid(ipart,itri),zmid(ipart,itri),dmid(ipart,itri)

      ! 2nd triangle
      itri = itri+1_int32
      xtri(ipart,itri,1:3) = [Pt4(1),Pt2(1),Pt3(1)]
      ytri(ipart,itri,1:3) = [Pt4(2),Pt2(2),Pt3(2)]
      ztri(ipart,itri,1:3) = [Pt4(3),Pt2(3),Pt3(3)]
      check_tri(ipart,itri) = 1  !all 2d triangles are checked
      write(iu_ptri,*) itri
      write(iu_ptri,*) xtri(ipart,itri,1),ytri(ipart,itri,1),ztri(ipart,itri,1)
      write(iu_ptri,*) xtri(ipart,itri,2),ytri(ipart,itri,2),ztri(ipart,itri,2)
      write(iu_ptri,*) xtri(ipart,itri,3),ytri(ipart,itri,3),ztri(ipart,itri,3)
      ! Define triangle midpoint
      xmid(ipart,itri) = sum(xtri(ipart,itri,1:3))/3._real64
      ymid(ipart,itri) = sum(ytri(ipart,itri,1:3))/3._real64
      zmid(ipart,itri) = sum(ztri(ipart,itri,1:3))/3._real64
      dmids(1) = sqrt( (xtri(ipart,itri,1) - xmid(ipart,itri))**2 + &
           (ytri(ipart,itri,1) - ymid(ipart,itri))**2 + (ztri(ipart,itri,1) - zmid(ipart,itri))**2 )
      dmids(2) = sqrt( (xtri(ipart,itri,2) - xmid(ipart,itri))**2 + &
           (ytri(ipart,itri,2) - ymid(ipart,itri))**2 + (ztri(ipart,itri,2) - zmid(ipart,itri))**2 )
      dmids(3) = sqrt( (xtri(ipart,itri,3) - xmid(ipart,itri))**2 + &
           (ytri(ipart,itri,3) - ymid(ipart,itri))**2 + (ztri(ipart,itri,3) - zmid(ipart,itri))**2 )
      dmid(ipart,itri) = maxval(dmids,1)
      write(iu_ptmid,*) itri,xmid(ipart,itri),ymid(ipart,itri),zmid(ipart,itri),dmid(ipart,itri)

      itri = itri+1_int32

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
Integer(int32), Intent(in) :: ntor, npol
Character(len=100), Intent(out) :: label
Integer(int32),Intent(out) :: msym  
Real(real64),Dimension(ntor,npol), Intent(out) :: &
  Rpart, Zpart, Phipart
Integer(int32) :: itor, ipol, ntor_dum, npol_dum
! Local variables
Real(real64) :: rshift, zshift, Phitmp
! Local parameters
Real(real64), Parameter :: pi = 3.141592653589793238462643383279502_real64
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
Rpart = (Rpart+rshift)*0.01_real64
Zpart = (Zpart+zshift)*0.01_real64
Phipart = Phipart*pi/180._real64

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
Integer(int32), Intent(in) :: ntor,npol
Character(len=100), Intent(out) :: label
Integer(int32),Intent(out) :: msym
Real(real64),Dimension(ntor,npol), Intent(out) :: &
  Rpart, Zpart, Phipart
Integer(int32) :: itor, ipol, ntor_dum, npol_dum
! Local variables
Real(real64) :: rshift, zshift
! Local parameters
Real(real64), Parameter :: pi = 3.141592653589793238462643383279502_real64
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
Subroutine query_part(fname,ntor,npol)
Use kind_mod
Use io_unit_spec, Only : iu_thispart
Implicit none
Character(len=100), Intent(in) :: fname
Integer(int32), Intent(out) :: ntor,npol
Integer(int32) :: nfp
Real(real64) :: rshift, zshift
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
  dL = Dsqrt( (R1-R2)*(R1-R2) + (Z1-Z2)*(Z1-Z2) )
  
  V1(1) = R2-R1
  V1(2) = Z2-Z1
  V1 = V1/dL

  Rmid = R1 + V1(1)*dL/2._real64
  Zmid = Z1 + V1(2)*dL/2._real64

  unit_out(1) = -V1(2)
  unit_out(2) =  V1(1)
  unit_out = unit_out/Dsqrt( V1(1)*V1(1) + V1(2)*V1(2) )

  Rout(i,npol+1) = Rmid + real(dir,real64)*unit_out(1)*step
  Zout(i,npol+1) = Zmid + real(dir,real64)*unit_out(2)*step
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
Real(real64),Allocatable :: Rpart(:,:),Zpart(:,:),Ppart(:,:)
Integer(int32) :: i, j
Integer(int32) :: ipart
Integer(int32) :: ntor, npol
Character(len=300) :: part_name
Character(len=300) :: label
! Local Parameters
Real(real64), Parameter :: pi = 3.141592653589793238462643383279502_real64
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
        Ppart(i,j) = Ppart(i,j) + 2._real64*pi/Real(msym,real64)
      Enddo
      Ppart(i,j) = Mod(Ppart(i,j),2._real64*pi/Real(msym,real64))
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


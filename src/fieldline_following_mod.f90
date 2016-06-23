!-----------------------------------------------------------------------------
!+ Fieldline following routines
!-----------------------------------------------------------------------------
Module fieldline_following_mod
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
Subroutine fl_derivs_dphi_bint(n,phi,RZ,df)
Use kind_mod
Use bfield_xdr, Only: bint
!Use vmec_coils_mod
Implicit None
Real(rknd), Intent(In) :: phi
Integer(iknd), Intent(In) :: n
Real(rknd), Intent(In), Dimension(n) :: RZ
Real(rknd), Intent(Out), Dimension(n) :: df

Real(rknd), Dimension(3) :: xvec, bval
Integer(iknd) :: idiv
!Real(rknd) :: Bx, By, Bz

idiv = 0
bval = 0._rknd

  xvec(1) = RZ(1)
  xvec(2) = phi
  xvec(3) = RZ(2)
  Call bint(xvec,bval,idiv)

If (idiv .ne. 0) Then
  Write(6,*) 'Error in fl_derivs_dphi_bint'
  Stop
Endif
df(1) = RZ(1)*bval(1)/bval(2)
df(2) = RZ(1)*bval(3)/bval(2)
Return
End Subroutine fl_derivs_dphi_bint
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

End Module fieldline_following_mod

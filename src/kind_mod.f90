!-----------------------------------------------------------------------------
!+ Kind specifications
!-----------------------------------------------------------------------------
Module kind_mod
!
! Description:
!   This module contains the kind specifications

! History:
! Version   Date      Comment
! -------   ----      -------
!  1.0    07/14/2011  From PENTA
! 
! Author(s): J. Lore 7/2009 - 5/24/2010    
!
  Implicit none
  Integer, parameter :: rknd = selected_real_kind(15,307) 
  Integer, parameter :: iknd = selected_int_kind(8)       
End module kind_mod
!- End of header -------------------------------------------------------------

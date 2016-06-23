!-----------------------------------------------------------------------------
!+ Sets file I/O unit numbers
!-----------------------------------------------------------------------------
Module io_unit_spec

!
! Description:
!   This module sets the file i/o numbers to avoid conflicts.
!
! History:
! Version   Date        Comment
! -------   ----        -------
!  1.0     01/07/2009   Original Code.  JL
! 
! Author(s): J. Lore 7/2009 - ***
!

  Implicit none
  
  Integer, parameter ::   &
    iu_nl=10,             &   ! Run settings namelist file (input)
    iu_plist=11,          &   ! Parts list file (input)
    iu_thispart=12,       &   ! Used to open each part file (input)
    iu_parts=30,          &   ! Parts data file (output)
    iu_nhit=31,           &   ! Number of hit lines (output)
    iu_vhit=31,           &   ! Vessel intersection data (output)
    iu_surf=33,           &   ! Surface data (output)
    iu_hit=34,            &   ! Hit line data (output)
    iu_int=35,            &   ! Intersection data (output)
    iu_launch=36,         &   ! Launch points (output)
    iu_axis=37,           &   ! Axis data file (output)
    iu_ptri=38,           &   ! Part_triangles (output)
    iu_ptmid=39               ! Part triangle midpoints (output)

End module io_unit_spec
!- End of header -------------------------------------------------------------

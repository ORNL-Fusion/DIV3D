!-----------------------------------------------------------------------------
!+ Module for physical constants used in several routines
!-----------------------------------------------------------------------------
Module phys_const
!
! Description:
!   This module contains physical constant parameters.
! 
! Author(s): J. Lore 7/2009 - 12/7/2015
!     
  ! Modules used:
  Use kind_mod, Only : real64
  
  Implicit none
  Public 
  ! Scalar parameters
  Real(real64),parameter :: p_mass = 1.672621637e-27_real64      !proton mass
  Real(real64),parameter :: e_mass = 9.10938215e-31_real64       !electron mass
  Real(real64),parameter :: elem_charge = 1.602176487e-19_real64 !elem. charge 
  Real(real64),parameter :: eps0 = 8.854187817e-12_real64        !electric const
  Real(real64),parameter :: pi = 3.14159265358979323846_real64   !pi
End module phys_const

!- End of header -------------------------------------------------------------

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!    Definition of some global variables for the FIST code

MODULE fist_global
  
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: set_fist_global, fistpar
  
  TYPE ( global_environment_type ) :: fistpar
  
CONTAINS

!******************************************************************************

SUBROUTINE set_fist_global ( glopar )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( IN ) :: glopar
  
!------------------------------------------------------------------------------
  
  fistpar = glopar
  
END SUBROUTINE set_fist_global

!******************************************************************************

END MODULE fist_global

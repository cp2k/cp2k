!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!***** cp2k/ao_types [1.0] *
!!
!!   NAME
!!     ao_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!   NOTES
!!
!!
!!   SOURCE
!******************************************************************************

MODULE ao_types
  
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ao_type, ao_structure_factor_type

  TYPE ao_structure_factor_type
     COMPLEX (KIND=dp), DIMENSION ( : ), POINTER :: ex, ey, ez
     INTEGER, DIMENSION ( : ), POINTER :: centre
  END TYPE ao_structure_factor_type

  TYPE ao_type
     REAL (KIND=dp), DIMENSION ( : ), POINTER :: cr
     REAL (KIND=dp), DIMENSION ( : ), POINTER :: norm
     REAL (KIND=dp), DIMENSION ( : ), POINTER :: mass
     INTEGER, DIMENSION ( : ), POINTER :: coef_to_basis 
     INTEGER, DIMENSION ( : ), POINTER :: coef_to_set 
     INTEGER, DIMENSION ( : ), POINTER :: coef_to_part 
     INTEGER, DIMENSION ( : ), POINTER :: grid_type_id
     INTEGER, DIMENSION ( : ), POINTER :: grid_id
     TYPE ( ao_structure_factor_type ), DIMENSION ( : ),  POINTER :: eigr
  END TYPE ao_type


END MODULE ao_types

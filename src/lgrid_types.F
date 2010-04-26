!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Scratch space for multiple threads writing to rs grids (see
!>        qs_collocate_density.F for an example
!> \par History
!>      IAB 26-Apr-2010 : initial version - moved out of qs_collocate_density.F
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2010 on behalf of the HECToR project
!> \author IAB
! *****************************************************************************

MODULE lgrid_types

  USE kinds,                           ONLY: dp

  IMPLICIT NONE

  PUBLIC :: lgrid_type, lgrid_p_type

  TYPE lgrid_type
     INTEGER :: ldim
     REAL(dp), DIMENSION(:), POINTER :: r
  END TYPE lgrid_type

  TYPE lgrid_p_type
     TYPE(lgrid_type), POINTER :: l
  END TYPE lgrid_p_type

END MODULE


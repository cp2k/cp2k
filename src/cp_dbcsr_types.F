!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR types in CP2K
!> \author  VW
!> \date    2010
!> \version 0.1
!>
!> <b>Modification history:</b>
!> - Created 2010
! *****************************************************************************
MODULE cp_dbcsr_types
  USE dbcsr_dist_operations
  USE dbcsr_io
  USE dbcsr_types,                     ONLY: dbcsr_obj

  !$ USE OMP_LIB
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_types'


  PUBLIC :: cp_dbcsr_obj_p_type

!> \brief Emulation of real_matrix_p_type
  TYPE cp_dbcsr_obj_p_type
     TYPE(dbcsr_obj), POINTER :: matrix
  END TYPE cp_dbcsr_obj_p_type

END MODULE cp_dbcsr_types

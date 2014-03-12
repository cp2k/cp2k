!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR types in CP2K
!> \author  VW
!> \date    2010
!> \version 0.1
!>
!> <b>Modification history:</b>
!> - Created 2010
!> - Moved basic routines from _operations to
!> - 2010-08-12 [UB] Added buffer definitions and mapping type.
!> - 2011-11-14 [UB] Removed buffers and maps
! *****************************************************************************
MODULE cp_dbcsr_types
  USE dbcsr_types,                     ONLY: cp_dbcsr_iterator=>dbcsr_iterator,&
                                             dbcsr_obj

  !$ USE OMP_LIB
#include "../cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_types'

  PUBLIC :: cp_dbcsr_p_type
  PUBLIC :: cp_dbcsr_iterator
  PUBLIC :: cp_dbcsr_type

  PRIVATE

  TYPE cp_dbcsr_p_type
     TYPE(cp_dbcsr_type), POINTER :: matrix
  END TYPE cp_dbcsr_p_type

  TYPE cp_dbcsr_type
     INTEGER                              :: ref_count
     TYPE(dbcsr_obj)                      :: matrix
  END TYPE cp_dbcsr_type

END MODULE cp_dbcsr_types

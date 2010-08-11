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
  USE dbcsr_types,                     ONLY: cp_dbcsr_iterator=>dbcsr_iterator,&
                                             dbcsr_obj

  !$ USE OMP_LIB
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_types'

  PUBLIC :: cp_dbcsr_p_type
  PUBLIC :: cp_dbcsr_iterator
  PUBLIC :: cp_dbcsr_type

  TYPE cp_dbcsr_p_type
     TYPE(cp_dbcsr_type), POINTER :: matrix
  END TYPE cp_dbcsr_p_type

  TYPE cp_dbcsr_type
     INTEGER                      :: ref_count
     TYPE(dbcsr_obj)              :: matrix
     !TYPE(mapping_type)           :: mapping
     !TYPE(dbcsr_block_buffer_obj) :: buffer
     !TYPE(dbcsr_distribution_obj) :: dbcsr_distribution
     !TYPE(dbcsr_distribution_obj) :: cp_distribution
     !type(array_i1d_obj)          :: cp_col_offsets, cp_row_offsets
     !type(array_i1d_obj)          :: cp_col_sizes, cp_row_sizes
  END TYPE cp_dbcsr_type

END MODULE cp_dbcsr_types

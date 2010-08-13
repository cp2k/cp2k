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
!> - Moved basic routines from _operations to 
!> - 2010-08-12 [UB] Added buffer definitions and mapping type.
! *****************************************************************************
MODULE cp_dbcsr_types
  USE cluster_types,                   ONLY: cluster_type
  USE dbcsr_types,                     ONLY: cp_dbcsr_iterator=>dbcsr_iterator,&
                                             dbcsr_block_buffer_obj,&
                                             dbcsr_obj

  !$ USE OMP_LIB
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_types'

  PUBLIC :: cp_dbcsr_p_type
  PUBLIC :: cp_dbcsr_iterator
  PUBLIC :: cp_dbcsr_type


  TYPE cp_dbcsr_mapping_type
     INTEGER :: refcount
     LOGICAL :: defined, active
     TYPE(cluster_type) :: clusters
  END TYPE cp_dbcsr_mapping_type


  TYPE cp_dbcsr_p_type
     TYPE(cp_dbcsr_type), POINTER :: matrix
  END TYPE cp_dbcsr_p_type

  TYPE cp_dbcsr_type
     INTEGER                              :: ref_count
     TYPE(dbcsr_obj)                      :: matrix
     TYPE(cp_dbcsr_mapping_type), POINTER :: row_mapping, col_mapping
     TYPE(dbcsr_block_buffer_obj)         :: buffers
     !TYPE(dbcsr_distribution_obj) :: dbcsr_distribution
     !TYPE(dbcsr_distribution_obj) :: cp_distribution
     !type(array_i1d_obj)          :: cp_col_offsets, cp_row_offsets
     !type(array_i1d_obj)          :: cp_col_sizes, cp_row_sizes
  END TYPE cp_dbcsr_type

END MODULE cp_dbcsr_types

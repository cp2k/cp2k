!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Methods for the CP2K DBCSR type
!> \author  Urban Borstnik
!> \date    2010-08-11
!> \version 0.1
!>
!> <b>Modification history:</b>
!> - Created 2010-08-11 [UB]
! *****************************************************************************
MODULE cp_dbcsr_methods
  USE array_types,                     ONLY: array_i1d_obj
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type

  !$ USE OMP_LIB
#include "cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_methods'

  PRIVATE

!  PUBLIC :: cp_dbcsr_get_mapping, cp_dbcsr_get_dbcsr

  PUBLIC :: cp_dbcsr_col_block_sizes, cp_dbcsr_row_block_sizes
  PUBLIC :: cp_dbcsr_col_block_offsets, cp_dbcsr_row_block_offsets

CONTAINS

  FUNCTION cp_dbcsr_row_block_sizes (matrix) RESULT (row_blk_sizes)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: row_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_row_block_sizes', &
      routineP = moduleN//':'//routineN

    row_blk_sizes = matrix%matrix%m%row_blk_size
  END FUNCTION cp_dbcsr_row_block_sizes

  FUNCTION cp_dbcsr_col_block_sizes (matrix) RESULT (col_blk_sizes)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: col_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_col_block_sizes', &
      routineP = moduleN//':'//routineN

    col_blk_sizes = matrix%matrix%m%col_blk_size
  END FUNCTION cp_dbcsr_col_block_sizes

  FUNCTION cp_dbcsr_col_block_offsets (matrix) RESULT (col_blk_offsets)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: col_blk_offsets

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_col_block_offsets', &
      routineP = moduleN//':'//routineN

    col_blk_offsets = matrix%matrix%m%col_blk_offset
  END FUNCTION cp_dbcsr_col_block_offsets

  FUNCTION cp_dbcsr_row_block_offsets (matrix) RESULT (row_blk_offsets)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: row_blk_offsets

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_row_block_offsets', &
      routineP = moduleN//':'//routineN

    row_blk_offsets = matrix%matrix%m%row_blk_offset
  END FUNCTION cp_dbcsr_row_block_offsets

END MODULE cp_dbcsr_methods

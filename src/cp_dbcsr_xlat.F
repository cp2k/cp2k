!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Utility for clustering translations and mapping.
!> \author  Urban Borstnik
!> \date    2010-08-11
!> \version 0.1
!>
!> <b>Modification history:</b>
!> - Created 2010-08-11 [UB]
! *****************************************************************************
MODULE cp_dbcsr_xlat
  USE array_types,                     ONLY: array_data
  USE cluster_buffers,                 ONLY: cluster_buffer_flush,&
                                             cluster_buffer_set_ptr
  USE cluster_methods,                 ONLY: cluster_part2cluster
  USE cp_dbcsr_methods,                ONLY: cp_dbcsr_col_block_offsets,&
                                             cp_dbcsr_col_block_sizes,&
                                             cp_dbcsr_has_mapping,&
                                             cp_dbcsr_row_block_offsets,&
                                             cp_dbcsr_row_block_sizes
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_mapping_type,&
                                             cp_dbcsr_type
  USE dbcsr_block_buffers,             ONLY: dbcsr_buffers_mark_dirty,&
                                             dbcsr_buffers_resize
  USE dbcsr_error_handling,            ONLY: dbcsr_error_type
  USE dbcsr_methods,                   ONLY: dbcsr_col_block_offsets,&
                                             dbcsr_row_block_offsets,&
                                             dbcsr_row_block_sizes
  USE dbcsr_ptr_util,                  ONLY: dbcsr_ptr_remapping,&
                                             pointer_d_rank_remap2
  USE dbcsr_util,                      ONLY: swap
  USE kinds,                           ONLY: real_8

  !$ USE OMP_LIB
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_xlat'

  LOGICAL, PRIVATE, PARAMETER :: careful_mod = .TRUE.

  PUBLIC :: set_block_pointer_to_cluster

CONTAINS


! *****************************************************************************
!> \brief Sets a block pointer to a region of a clustered block.
!> \par Notes
!>      This routine does all the necessary lookups of blocks sizes
!>      and offsets and then calls another routine for the actual
!>      pointer assignments.
!> \par History
!>      - created 2010-08-10
!> \author UB
! *****************************************************************************
  SUBROUTINE set_block_pointer_to_cluster(matrix,&
       high_row, high_col, high_block, high_tr, contiguous_pointers, read_only,&
       low_block, low_row, low_col, low_tr,&
       high_row_size, high_col_size, low_col_size, low_row_size,&
       error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: high_row, high_col
    REAL(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: high_block
    LOGICAL, INTENT(IN)                      :: high_tr, contiguous_pointers, &
                                                read_only
    REAL(kind=real_8), DIMENSION(:), POINTER :: low_block
    INTEGER, INTENT(IN), OPTIONAL            :: low_row, low_col
    LOGICAL, INTENT(IN)                      :: low_tr
    INTEGER, INTENT(IN), OPTIONAL            :: high_row_size, high_col_size, &
                                                low_col_size, low_row_size
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_block_pointer_to_cluster', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: actual_cs, actual_rs, &
                                                col_offset, high_cs, high_rs, &
                                                low_c, low_cs, low_r, low_rs, &
                                                row_offset
    LOGICAL                                  :: failure
    REAL(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: low_block_2d
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL get4blksizes (matrix, high_row, high_col, low_row, low_col,&
       high_row_size, high_col_size, low_row_size, low_col_size,&
       high_rs, high_cs, low_r, low_c, low_rs, low_cs, row_offset, col_offset,&
       error=error)
    !
    CALL cluster_buffer_flush (matrix%buffers, error=error)
    actual_rs = high_rs
    actual_cs = high_cs
    IF (high_tr) CALL swap (actual_rs, actual_cs)
    IF (.NOT. contiguous_pointers .AND. dbcsr_ptr_remapping) THEN
       ! point high pointer into cluster block
       !high_block => low_block(low_row_offset, low_col_offset)
       CALL pointer_d_rank_remap2 (low_block_2d,&
            actual_rs, actual_cs, low_block)
       high_block => low_block_2d(row_offset:row_offset+high_rs-1,&
                     col_offset:col_offset+high_cs-1)
    ELSE
       ! Resize the buffers to correct row size
       CALL dbcsr_buffers_resize (matrix%buffers,&
            row_size=actual_rs, col_size=actual_cs,&
            minimum=.NOT.contiguous_pointers, error=dbcsr_error)
       CALL cluster_buffer_set_ptr (high_block, high_row, high_col,&
            high_rs, high_cs, high_tr, &
            low_rs, low_cs, low_tr, row_offset, col_offset,&
            matrix%buffers, error=error)
       IF (read_only) THEN
          CALL dbcsr_buffers_mark_dirty (matrix%buffers, dirty=.FALSE.,&
               error=dbcsr_error)
       ENDIF
    ENDIF
  END SUBROUTINE set_block_pointer_to_cluster


  SUBROUTINE get4blksizes (matrix, high_row, high_col, low_row, low_col,&
       high_row_size, high_col_size,&
       low_row_size, low_col_size,&
       high_rs, high_cs, low_r, low_c, low_rs, low_cs, row_off, col_off, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, INTENT(IN)                      :: high_row, high_col
    INTEGER, INTENT(IN), OPTIONAL            :: low_row, low_col, &
                                                high_row_size, high_col_size, &
                                                low_row_size, low_col_size
    INTEGER, INTENT(OUT)                     :: high_rs, high_cs, low_r, &
                                                low_c, low_rs, low_cs, &
                                                row_off, col_off
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get4blksizes', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(:), POINTER           :: cbs, high_offset, low_offset, &
                                                rbs

    IF (.NOT. PRESENT (high_row_size)) THEN
       rbs => array_data (cp_dbcsr_row_block_sizes (matrix))
       high_rs = rbs(high_row)
    ELSE
       high_rs = high_row_size
    ENDIF
    IF (.NOT. PRESENT (high_col_size)) THEN
       cbs => array_data (cp_dbcsr_col_block_sizes (matrix))
       high_cs = cbs(high_col)
    ELSE
       high_cs = high_col_size
    ENDIF
    !
    IF (.NOT. PRESENT (low_row)) THEN
       IF (careful_mod) THEN
          CALL cp_assert (ASSOCIATED (matrix%row_mapping), cp_fatal_level,&
               cp_internal_error, routineN, "Row mapping nonexistent.",&
               error=error)
       ENDIF
       low_r = cp_part2cluster (matrix%row_mapping, high_row, error)
    ELSE
       low_r = low_row
    ENDIF
    IF (.NOT. PRESENT (low_col)) THEN
       IF (careful_mod) THEN
          CALL cp_assert (ASSOCIATED (matrix%col_mapping), cp_fatal_level,&
               cp_internal_error, routineN, "Column mapping nonexistent.",&
               error=error)
       ENDIF
       low_c = cp_part2cluster (matrix%col_mapping, high_col, error)
    ELSE
       low_c = low_col
    ENDIF
    !
    IF (.NOT. PRESENT (low_row_size)) THEN
       rbs => array_data (dbcsr_row_block_sizes (matrix%matrix))
       low_rs = rbs(low_r)
    ELSE
       low_rs = low_row_size
    ENDIF
    IF (.NOT. PRESENT (low_col_size)) THEN
       cbs => array_data (dbcsr_row_block_sizes (matrix%matrix))
       low_cs = cbs(low_c)
    ELSE
       low_cs = low_col_size
    ENDIF
    !
    low_offset => array_data (dbcsr_row_block_offsets (matrix%matrix))
    high_offset => array_data (cp_dbcsr_row_block_offsets (matrix))
    row_off = high_offset(high_row) - low_offset(low_r) + 1
    !
    low_offset => array_data (dbcsr_col_block_offsets (matrix%matrix))
    high_offset => array_data (cp_dbcsr_col_block_offsets (matrix))
    col_off = high_offset(high_col) - low_offset(low_c) + 1
  END SUBROUTINE get4blksizes


  FUNCTION cp_part2cluster (mapping, particle, error) RESULT (cluster)
    TYPE(cp_dbcsr_mapping_type), INTENT(IN)  :: mapping
    INTEGER, INTENT(IN)                      :: particle
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER                                  :: cluster

    IF (cp_dbcsr_has_mapping (mapping)) THEN
       cluster = cluster_part2cluster (mapping%clusters, particle, error)
    ELSE
       cluster = particle
    ENDIF
  END FUNCTION cp_part2cluster

END MODULE cp_dbcsr_xlat

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Cluster pointer and buffer operations.
!> \par History
!> \author UB (2010-08-10)
! *****************************************************************************
MODULE cluster_buffers

  USE dbcsr_block_buffers,             ONLY: dbcsr_buffers_get_area,&
                                             dbcsr_buffers_get_num,&
                                             dbcsr_buffers_is_dirty,&
                                             dbcsr_buffers_new,&
                                             dbcsr_buffers_release,&
                                             dbcsr_buffers_valid,&
                                             dbcsr_buffers_which_is_my
  USE dbcsr_block_operations,          ONLY: dbcsr_block_partial_copy
  USE dbcsr_data_methods,              ONLY: dbcsr_data_get_sizes,&
                                             dbcsr_data_get_type,&
                                             dbcsr_data_query_type,&
                                             dbcsr_data_reset_type,&
                                             dbcsr_data_set_pointer
  USE dbcsr_error_handling,            ONLY: dbcsr_error_type
  USE dbcsr_types,                     ONLY: dbcsr_block_buffer_obj,&
                                             dbcsr_data_obj,&
                                             dbcsr_type_real_8_2d
  USE dbcsr_util,                      ONLY: swap
  USE kinds,                           ONLY: dp,&
                                             real_8

  !$ USE OMP_LIB

#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cluster_buffers'

  LOGICAL, PARAMETER, PRIVATE          :: careful_mod = .TRUE.,&
                                          debug_mod = .TRUE.

! *** Public subroutines ***

  PUBLIC :: cluster_buffer_new, cluster_buffer_release
  PUBLIC :: cluster_buffer_flush
  PUBLIC :: cluster_buffer_set_ptr

CONTAINS

  SUBROUTINE cluster_buffer_new (buffers, data_area, error)
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(OUT)                            :: buffers
    TYPE(dbcsr_data_obj), INTENT(IN)         :: data_area
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_buffers_new (buffers, data_area, error=dbcsr_error)
  END SUBROUTINE cluster_buffer_new

  SUBROUTINE cluster_buffer_release (buffers, error)
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(INOUT)                          :: buffers
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    TYPE(dbcsr_error_type)                   :: dbcsr_error

! Buffers don't exist if the matrix was not created.

    IF (dbcsr_buffers_valid (buffers)) THEN
       CALL dbcsr_buffers_release (buffers, error=dbcsr_error)
    ENDIF
  END SUBROUTINE cluster_buffer_release

! *****************************************************************************
!> \brief Sets a block pointer to a region of a clustered block.
!> \par History
!>      - created 2010-08-10
!> \author UB
! *****************************************************************************
  SUBROUTINE cluster_buffer_set_ptr (high_block, high_row, high_col,&
       high_row_size, high_col_size, high_tr,&
       low_row_size, low_col_size, low_tr, row_offset, col_offset,&
       buffers, error)
    REAL(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: high_block
    INTEGER, INTENT(IN)                      :: high_row, high_col, &
                                                high_row_size, high_col_size
    LOGICAL, INTENT(IN)                      :: high_tr
    INTEGER, INTENT(IN)                      :: low_row_size, low_col_size
    LOGICAL, INTENT(IN)                      :: low_tr
    INTEGER, INTENT(IN)                      :: row_offset, col_offset
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(INOUT)                          :: buffers
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_buffer_set_ptr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: buff_cs, buff_rs, cs, ibuff, &
                                                rs, transposes
    LOGICAL                                  :: failure, valid
    REAL(kind=real_8), DIMENSION(:), POINTER :: low_block
    TYPE(dbcsr_data_obj)                     :: backing, buff_area
    TYPE(dbcsr_error_type)                   :: dbcsr_error

!   ---------------------------------------------------------------------------

    ibuff = dbcsr_buffers_which_is_my (buffers)
    buff_area = dbcsr_buffers_get_area (buffers, ibuff, error=dbcsr_error)
    !
    !CALL dbcsr_assert (low_tr, "EQV", buff_tr, dbcsr_fatal_level, &
    !     dbcsr_wrong_args_error, routineN,&
    !     "Source and buffer should have same transposed status.",&
    !     __LINE__, error=error)
    IF (careful_mod) THEN
       failure = .FALSE.
       CALL cp_assert (ASSOCIATED (buff_area%d%r2_dp),&
            cp_fatal_level, cp_wrong_args_error, routineN,&
            "Buffer not associated.", failure=failure, error=error)
       failure = .TRUE.
       CALL cp_assert (ASSOCIATED (buffers%b%backing),&
            cp_fatal_level, cp_internal_error, routineN,&
            "Individual backing not associated.", failure=failure, error=error)
    ENDIF
    CALL dbcsr_data_get_sizes (buff_area, rs, cs, valid, error=dbcsr_error)
    CALL cp_assert (valid, cp_fatal_level, cp_internal_error,&
         routineN, "Buffer not setup", error=error)
    !
    buff_rs = high_row_size
    buff_cs = high_col_size
    IF (high_tr) CALL swap (buff_rs, buff_cs)
    !
    backing = buffers%b%backing(ibuff)
    IF (dbcsr_data_get_type (backing) .NE. dbcsr_data_query_type (low_block)) &
         CALL dbcsr_data_reset_type (backing, dbcsr_data_query_type (low_block))
    CALL dbcsr_data_set_pointer (backing, low_block)
    CALL dbcsr_block_partial_copy (dst=buff_area,&
         dst_rs=high_row_size, dst_cs=high_col_size, dst_tr=high_tr,&
         src=backing, src_offset=0,&
         src_rs=low_row_size, src_cs=low_col_size, src_tr=low_tr,&
         dst_r_lb=1, dst_c_lb=1,&
         src_r_lb=row_offset, src_c_lb=col_offset,&
         nrow = high_row_size, ncol = high_col_size)
    !high_block => buff_area%d%[base1]2_[prec1](1:buff_rs, 1:buff_cs)
    CALL cp_assert (dbcsr_data_get_type (buff_area) .EQ. dbcsr_type_real_8_2d,&
         cp_fatal_level, cp_internal_error, routineN, "type NYI", error=error)
    high_block => buff_area%d%r2_dp(1:buff_rs, 1:buff_cs)
    transposes = 0
    IF (high_tr) transposes = IBSET (transposes, 0)
    IF (low_tr) transposes = IBSET (transposes, 1)
    buffers%b%rcb(1, ibuff) = high_row
    buffers%b%rcb(2, ibuff) = high_col
    buffers%b%rcb(3, ibuff) = 1 ! offset
    buffers%b%rcb(4, ibuff) = high_row_size
    buffers%b%rcb(5, ibuff) = high_col_size
    buffers%b%rcb(6, ibuff) = transposes
    buffers%b%rcb(7, ibuff) = low_row_size
    buffers%b%rcb(8, ibuff) = low_col_size
    buffers%b%rcb(9, ibuff) = row_offset
    buffers%b%rcb(10, ibuff) = col_offset
    buffers%b%dirty(ibuff) = .TRUE.
  END SUBROUTINE cluster_buffer_set_ptr


! *****************************************************************************
!> \brief Flushes buffers
!> \param[in,out] buffers     Buffers to flush
!> \param[in,out] error       error
!> \author UB
! *****************************************************************************
  SUBROUTINE cluster_buffer_flush (buffers, error)
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(INOUT)                          :: buffers
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cluster_buffer_flush', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ibuff, nthreads
    LOGICAL                                  :: reap
    TYPE(dbcsr_error_type)                   :: dbcsr_error

!   ---------------------------------------------------------------------------

    IF (.NOT. dbcsr_buffers_valid (buffers)) THEN
       RETURN
    ENDIF
    reap = .FALSE.
    nthreads = 1
!$  nthreads = omp_get_num_threads()
    reap = nthreads .LT. dbcsr_buffers_get_num(buffers)
!$  reap = reap .OR. .NOT. omp_in_parallel()
    IF (reap) THEN
!$omp master
       DO ibuff = 1, dbcsr_buffers_get_num(buffers)
          IF (buffers%b%dirty(ibuff)) THEN
             CALL buffer_flush_single (buffers, ibuff, error=error)
          ENDIF
       ENDDO
!$omp end master
    ELSE
       ibuff = dbcsr_buffers_which_is_my (buffers)
       IF (dbcsr_buffers_is_dirty (buffers, ibuff, error=dbcsr_error)) THEN
          CALL buffer_flush_single (buffers, ibuff, error=error)
       ENDIF
    ENDIF
  END SUBROUTINE cluster_buffer_flush



! *****************************************************************************
!> \brief Copies data from the buffer to the backing store.
!> \author UB
! *****************************************************************************
  SUBROUTINE buffer_flush_single (buffers, which_buffer, error)
    TYPE(dbcsr_block_buffer_obj), &
      INTENT(INOUT)                          :: buffers
    INTEGER, INTENT(IN)                      :: which_buffer
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'buffer_flush_single', &
      routineP = moduleN//':'//routineN

    INTEGER :: col, col_offset, dst_cs, dst_rs, high_col_size, high_row_size, &
      low_col_size, low_row_size, offset, row, row_offset
    LOGICAL                                  :: high_tr, low_tr, valid
    TYPE(dbcsr_data_obj)                     :: area, backing

!   ---------------------------------------------------------------------------

    IF (debug_mod) THEN
       CALL cp_assert (ASSOCIATED (buffers%b), cp_fatal_level,&
            cp_caller_error, routineN, "Buffers not valid", error=error)
    ENDIF
    IF (careful_mod) THEN
       CALL cp_assert (buffers%b%dirty(which_buffer), cp_warning_level,&
            cp_caller_error, routineN, "Should not flush clean buffer.",&
            error=error)
       CALL cp_assert (ASSOCIATED (buffers%b%backing), cp_fatal_level,&
            cp_internal_error, routineN, "Individual backing nonexistent.",&
            error=error)
    ENDIF
    IF (buffers%b%dirty(which_buffer)) THEN
       backing       = buffers%b%backing(which_buffer)
       area          = buffers%b%buffers(which_buffer)
       row           = buffers%b%rcb(1, which_buffer)
       col           = buffers%b%rcb(2, which_buffer)
       offset        = buffers%b%rcb(3, which_buffer)
       high_row_size = buffers%b%rcb(4, which_buffer)
       high_col_size = buffers%b%rcb(5, which_buffer)
       low_row_size  = buffers%b%rcb(7, which_buffer)
       low_col_size  = buffers%b%rcb(8, which_buffer)
       row_offset    = buffers%b%rcb(9, which_buffer)
       col_offset    = buffers%b%rcb(10, which_buffer)
       low_tr = BTEST (buffers%b%rcb(6, which_buffer), 1)
       high_tr = BTEST (buffers%b%rcb(6, which_buffer), 0)
       CALL dbcsr_block_partial_copy (dst=backing,&
            dst_offset=ABS(offset)-1,&
            dst_rs = low_row_size, dst_cs = low_col_size, dst_tr = low_tr,&
            src = area, src_rs = high_row_size, src_cs = high_col_size,&
            src_tr = high_tr,&
            dst_r_lb=row_offset, dst_c_lb=col_offset, src_r_lb=1, src_c_lb=1,&
            nrow = high_row_size, ncol = high_col_size)
       buffers%b%dirty(which_buffer) = .FALSE.
    ENDIF
  END SUBROUTINE buffer_flush_single


END MODULE cluster_buffers

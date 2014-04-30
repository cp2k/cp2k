!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Routines for basic block transformations.
!> \author  Urban Borstnik
!> \date    2010-02-18
!> \version 0.9
!>
!> <b>Modification history:</b>
!> - 2009-05-12 moved from dbcsr_util
! *****************************************************************************
MODULE dbcsr_block_operations

  USE acc_devmem,                      ONLY: acc_devmem_allocated,&
                                             acc_devmem_setzero_bytes
  USE dbcsr_data_methods_low,          ONLY: &
       dbcsr_data_exists, dbcsr_data_get_size, &
       dbcsr_data_get_size_referenced, dbcsr_data_get_type, &
       dbcsr_data_verify_bounds, dbcsr_get_data_p_2d_c, &
       dbcsr_get_data_p_2d_d, dbcsr_get_data_p_2d_s, dbcsr_get_data_p_2d_z, &
       dbcsr_get_data_p_c, dbcsr_get_data_p_d, dbcsr_get_data_p_s, &
       dbcsr_get_data_p_z, dbcsr_type_2d_to_1d, dbcsr_type_is_2d
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: dp,&
                                             real_4,&
                                             real_8,&
                                             sp
  USE dbcsr_types,                     ONLY: &
       dbcsr_data_obj, dbcsr_datatype_sizeof, dbcsr_scalar_type, &
       dbcsr_type_complex_4, dbcsr_type_complex_4_2d, dbcsr_type_complex_8, &
       dbcsr_type_complex_8_2d, dbcsr_type_real_4, dbcsr_type_real_4_2d, &
       dbcsr_type_real_8, dbcsr_type_real_8_2d

  !$ USE OMP_LIB
  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_block_operations'

  PUBLIC :: dbcsr_block_transpose, dbcsr_data_set
  PUBLIC :: dbcsr_block_copy, dbcsr_block_partial_copy
  PUBLIC :: set_block_diagonal, dbcsr_data_clear
  PUBLIC :: set_block2d_diagonal, get_block2d_diagonal
  PUBLIC :: dbcsr_block_scale, dbcsr_block_conjg, dbcsr_block_real_neg

  PUBLIC :: block_set, dbcsr_data_copy
  PUBLIC :: block_add_on_diag
  PUBLIC :: block_chol_inv
  PUBLIC :: dbcsr_block_norm_frob

  PUBLIC :: block_add

  ! For quick access
  PUBLIC :: block_copy_s, block_copy_d,&
            block_copy_c, block_copy_z

  INTERFACE dbcsr_block_transpose
     MODULE PROCEDURE block_transpose_inplace_s, block_transpose_inplace_d,&
                      block_transpose_inplace_c, block_transpose_inplace_z
     MODULE PROCEDURE block_transpose_copy_d, block_transpose_copy_s,&
                      block_transpose_copy_z, block_transpose_copy_c
     MODULE PROCEDURE block_transpose_copy_2d1d_d,&
                      block_transpose_copy_2d1d_s,&
                      block_transpose_copy_2d1d_z,&
                      block_transpose_copy_2d1d_c
     MODULE PROCEDURE block_transpose_copy_1d2d_d,&
                      block_transpose_copy_1d2d_s,&
                      block_transpose_copy_1d2d_z,&
                      block_transpose_copy_1d2d_c
     MODULE PROCEDURE dbcsr_block_transpose_aa, dbcsr_block_transpose_a
  END INTERFACE

  INTERFACE dbcsr_block_copy
     MODULE PROCEDURE block_copy_2d1d_s, block_copy_2d1d_d,&
                      block_copy_2d1d_c, block_copy_2d1d_z
     MODULE PROCEDURE block_copy_1d2d_s, block_copy_1d2d_d,&
                      block_copy_1d2d_c, block_copy_1d2d_z
     MODULE PROCEDURE block_copy_1d1d_s, block_copy_1d1d_d,&
                      block_copy_1d1d_c, block_copy_1d1d_z
     MODULE PROCEDURE block_copy_2d2d_s, block_copy_2d2d_d,&
                      block_copy_2d2d_c, block_copy_2d2d_z
  END INTERFACE


  INTERFACE dbcsr_block_partial_copy
     MODULE PROCEDURE block_partial_copy_a
     MODULE PROCEDURE block_partial_copy_d,&
                      block_partial_copy_s,&
                      block_partial_copy_z,&
                      block_partial_copy_c
     MODULE PROCEDURE block_partial_copy_1d2d_d,&
                      block_partial_copy_1d2d_s,&
                      block_partial_copy_1d2d_z,&
                      block_partial_copy_1d2d_c
     MODULE PROCEDURE block_partial_copy_2d1d_d,&
                      block_partial_copy_2d1d_s,&
                      block_partial_copy_2d1d_z,&
                      block_partial_copy_2d1d_c
     MODULE PROCEDURE block_partial_copy_2d2d_d,&
                      block_partial_copy_2d2d_s,&
                      block_partial_copy_2d2d_z,&
                      block_partial_copy_2d2d_c
  END INTERFACE

  INTERFACE dbcsr_data_clear
     MODULE PROCEDURE dbcsr_data_clear_nt
     MODULE PROCEDURE dbcsr_data_clear0
  END INTERFACE

  ! Supports copy between two data areas, or to a data area from a
  ! given explicit array.
  INTERFACE dbcsr_data_set
     MODULE PROCEDURE dbcsr_data_copy_aa, dbcsr_data_set_as,&
          dbcsr_data_set_ad, dbcsr_data_set_ac, dbcsr_data_set_az
  END INTERFACE

  INTERFACE dbcsr_data_copy
     MODULE PROCEDURE dbcsr_data_copy_aa2, dbcsr_data_set_as,&
          dbcsr_data_set_ad, dbcsr_data_set_ac, dbcsr_data_set_az
  END INTERFACE


  INTERFACE set_block_diagonal
     MODULE PROCEDURE set_block_diagonal_s, set_block_diagonal_d,&
          set_block_diagonal_c, set_block_diagonal_z
  END INTERFACE

  INTERFACE get_block2d_diagonal
     MODULE PROCEDURE get_block2d_diagonal_a
     MODULE PROCEDURE get_block2d_diagonal_s, get_block2d_diagonal_d,&
                      get_block2d_diagonal_c, get_block2d_diagonal_z
  END INTERFACE

  INTERFACE set_block2d_diagonal
     MODULE PROCEDURE set_block2d_diagonal_a
     MODULE PROCEDURE set_block2d_diagonal_s, set_block2d_diagonal_d,&
                      set_block2d_diagonal_c, set_block2d_diagonal_z
  END INTERFACE

  INTERFACE block_set
     MODULE PROCEDURE block_2d_set_d, block_2d_set_s,&
                      block_2d_set_z, block_2d_set_c
     MODULE PROCEDURE block_set_d, block_set_s,&
                      block_set_z, block_set_c
  END INTERFACE

  INTERFACE block_add_on_diag
     MODULE PROCEDURE block_add_on_diag_anytype
     MODULE PROCEDURE block_2d_add_on_diag_d, block_2d_add_on_diag_s,&
                      block_2d_add_on_diag_z, block_2d_add_on_diag_c
     MODULE PROCEDURE block_add_on_diag_d, block_add_on_diag_s,&
                      block_add_on_diag_z, block_add_on_diag_c
  END INTERFACE

  INTERFACE block_chol_inv
     MODULE PROCEDURE block_2d_chol_inv_d, block_2d_chol_inv_s,&
                      block_2d_chol_inv_z, block_2d_chol_inv_c
     MODULE PROCEDURE block_chol_inv_d, block_chol_inv_s,&
                      block_chol_inv_z, block_chol_inv_c
  END INTERFACE


  INTERFACE dbcsr_block_norm_frob
     MODULE PROCEDURE block_norm_frob_anytype_s
     MODULE PROCEDURE block_norm_frob_s_s, block_norm_frob_d_s,&
                      block_norm_frob_c_s, block_norm_frob_z_s
  END INTERFACE

  INTERFACE block_add
     MODULE PROCEDURE block_add_anytype
     MODULE PROCEDURE block_add_s, block_add_d, block_add_c, block_add_z
  END INTERFACE block_add


  LOGICAL, PARAMETER :: debug_mod    = .FALSE.
  LOGICAL, PARAMETER :: careful_mod  = .FALSE.

CONTAINS

! *****************************************************************************
!> \brief Copy data from one data area to another.
!>
!> There are no checks done for correctness!
!> \param[in] dst        destination data area
!> \param[in] src        source data area
!> \param[in] row_size   row size of existing block
!> \param[in] col_size   column size of existing block
!> \param[in] lb         (optional) lower bound for destination (and source if
!>                       not given explicity)
!> \param[in] source_lb  (optional) lower bound of source
!> \param[in] scale      (optional) scale data
!> \param[in] lb2        (optional) lower bound of 2nd dimension for target
!> \param[in] source_lb2 (optional) lower bound of 2nd dimension for source
! *****************************************************************************
  SUBROUTINE dbcsr_block_transpose_aa (dst, src, &
       row_size, col_size, lb, source_lb, scale, lb2, source_lb2, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dst
    TYPE(dbcsr_data_obj), INTENT(IN)         :: src
    INTEGER, INTENT(IN)                      :: row_size, col_size
    INTEGER, INTENT(IN), OPTIONAL            :: lb, source_lb
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: scale
    INTEGER, INTENT(IN), OPTIONAL            :: lb2, source_lb2
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_block_transpose_aa', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_size, lb2_s, lb2_t, &
                                                lb_s, lb_t, ub_s, ub_t

!   ---------------------------------------------------------------------------

    IF (debug_mod) THEN
       CALL dbcsr_assert (ASSOCIATED(dst%d), 'AND', ASSOCIATED(src%d),&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Data areas must be setup.",__LINE__,error)
       CALL dbcsr_assert (dst%d%data_type, 'EQ', src%d%data_type, dbcsr_fatal_level,&
            dbcsr_caller_error, routineN, "Data type must be the same.",__LINE__,error)
       CALL dbcsr_assert (dst%d%data_type.EQ.dbcsr_type_real_8&
            .OR.dst%d%data_type.EQ.dbcsr_type_real_4&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_8&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_4&
            .OR.dst%d%data_type.EQ.dbcsr_type_real_8_2d&
            .OR.dst%d%data_type.EQ.dbcsr_type_real_4_2d&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_8_2d&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_4_2d,&
            dbcsr_warning_level, dbcsr_caller_error, routineN, "Incorrect data type.",__LINE__,error)
    ENDIF
    IF (PRESENT (scale)) THEN
       IF (dbcsr_type_is_2d (src%d%data_type)) THEN
          CALL dbcsr_assert (scale%data_type, &
               'EQ', dbcsr_type_2d_to_1d (src%d%data_type),&
               dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
               "Incompatible data types",__LINE__,error)
       ELSE
          CALL dbcsr_assert (scale%data_type, 'EQ', src%d%data_type, dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN, "Incompatible data types",__LINE__,error)
       ENDIF
    ENDIF
    data_size = row_size*col_size
    lb_t = 1
    IF (PRESENT (lb)) lb_t = lb
    ub_t = lb_t + data_size - 1
    IF (PRESENT (source_lb)) THEN
       lb_s = source_lb
       ub_s = source_lb + data_size-1
    ELSE
       lb_s = lb_t
       ub_s = ub_t
    ENDIF
    lb2_t = 1
    IF (PRESENT (lb2)) lb2_t = lb2
    IF (PRESENT (source_lb2)) THEN
       lb2_s = source_lb2
    ELSE
       lb2_s = lb2_t
    ENDIF
    SELECT CASE (src%d%data_type)
    CASE (dbcsr_type_real_8)
       IF (PRESENT (scale)) THEN
          CALL dbcsr_block_transpose (dst%d%r_dp(lb_t:ub_t),&
               src%d%r_dp(lb_s:ub_s) * scale%r_dp,&
               row_size, col_size)
       ELSE
          CALL dbcsr_block_transpose (dst%d%r_dp(lb_t:ub_t),&
               src%d%r_dp(lb_s:ub_s),&
               row_size, col_size)
       ENDIF
    CASE (dbcsr_type_real_4)
       IF (PRESENT (scale)) THEN
          CALL dbcsr_block_transpose (dst%d%r_sp(lb_t:ub_t),&
               src%d%r_sp(lb_s:ub_s) * scale%r_sp,&
               row_size, col_size)
       ELSE
          CALL dbcsr_block_transpose (dst%d%r_sp(lb_t:ub_t),&
               src%d%r_sp(lb_s:ub_s),&
               row_size, col_size)
       ENDIF
    CASE (dbcsr_type_complex_8)
       IF (PRESENT (scale)) THEN
          CALL dbcsr_block_transpose (dst%d%c_dp(lb_t:ub_t),&
               src%d%c_dp(lb_s:ub_s) * scale%c_dp,&
               row_size, col_size)
       ELSE
          CALL dbcsr_block_transpose (dst%d%c_dp(lb_t:ub_t),&
               src%d%c_dp(lb_s:ub_s),&
               row_size, col_size)
       ENDIF
    CASE (dbcsr_type_complex_4)
       IF (PRESENT (scale)) THEN
          CALL dbcsr_block_transpose (dst%d%c_sp(lb_t:ub_t),&
               src%d%c_sp(lb_s:ub_s) * scale%c_sp,&
               row_size, col_size)
       ELSE
          CALL dbcsr_block_transpose (dst%d%c_sp(lb_t:ub_t),&
               src%d%c_sp(lb_s:ub_s),&
               row_size, col_size)
       ENDIF
    CASE (dbcsr_type_real_8_2d)
       IF (PRESENT (scale)) THEN
          dst%d%r2_dp(lb_t:lb_t+col_size-1,lb2_t:lb2_t+row_size-1) =&
               TRANSPOSE(&
               src%d%r2_dp(lb_s:lb_s+row_size-1,lb2_s:lb2_s+col_size-1) &
               * scale%r_dp)
       ELSE
          dst%d%r2_dp(lb_t:lb_t+col_size-1,lb2_t:lb2_t+row_size-1) =&
               TRANSPOSE(&
               src%d%r2_dp(lb_s:lb_s+row_size-1,lb2_s:lb2_s+col_size-1))
       ENDIF
    CASE (dbcsr_type_real_4_2d)
       IF (PRESENT (scale)) THEN
          dst%d%r2_sp(lb_t:lb_t+col_size-1,lb2_t:lb2_t+row_size-1) =&
               TRANSPOSE(&
               src%d%r2_sp(lb_s:lb_s+row_size-1,lb2_s:lb2_s+col_size-1) &
               * scale%r_sp)
       ELSE
          dst%d%r2_sp(lb_t:lb_t+col_size-1,lb2_t:lb2_t+row_size-1) =&
               TRANSPOSE(&
               src%d%r2_sp(lb_s:lb_s+row_size-1,lb2_s:lb2_s+col_size-1))
       ENDIF
    CASE (dbcsr_type_complex_8_2d)
       IF (PRESENT (scale)) THEN
          dst%d%c2_dp(lb_t:lb_t+col_size-1,lb2_t:lb2_t+row_size-1) =&
               TRANSPOSE(&
               src%d%c2_dp(lb_s:lb_s+row_size-1,lb2_s:lb2_s+col_size-1) &
               * scale%c_dp)
       ELSE
          dst%d%c2_dp(lb_t:lb_t+col_size-1,lb2_t:lb2_t+row_size-1) =&
               TRANSPOSE(&
               src%d%c2_dp(lb_s:lb_s+row_size-1,lb2_s:lb2_s+col_size-1))
       ENDIF
    CASE (dbcsr_type_complex_4_2d)
       IF (PRESENT (scale)) THEN
          dst%d%c2_sp(lb_t:lb_t+col_size-1,lb2_t:lb2_t+row_size-1) =&
               TRANSPOSE(&
               src%d%c2_sp(lb_s:lb_s+row_size-1,lb2_s:lb2_s+col_size-1) &
               * scale%c_sp)
       ELSE
          dst%d%c2_sp(lb_t:lb_t+col_size-1,lb2_t:lb2_t+row_size-1) =&
               TRANSPOSE(&
               src%d%c2_sp(lb_s:lb_s+row_size-1,lb2_s:lb2_s+col_size-1))
       ENDIF
    CASE default
       CALL dbcsr_assert (.FALSE.,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN, "Incorrect data type.",__LINE__,error)
    END SELECT
  END SUBROUTINE dbcsr_block_transpose_aa


! *****************************************************************************
!> \brief Scale a data area.
!>
!> There are no checks done for correctness!
!> \param[in] src        source data area
!> \param[in] row_size   row size of existing block
!> \param[in] col_size   column size of existing block
!> \param[in] lb         (optional) lower bound for destination (and source if
!>                       not given explicity)
!> \param[in] source_lb  (optional) lower bound of source
!> \param[in] scale      (optional) scale data
!> \param[in] lb2        (optional) lower bound of 2nd dimension for target
!> \param[in] source_lb2 (optional) lower bound of 2nd dimension for source
!> \par History
!> - 2010-09 [??] Copied from block_transpose?
!> - 2010-09-20 [UB] Swaps/corrects row/column definitions for 2-D bounds
! *****************************************************************************
  SUBROUTINE dbcsr_block_scale (dst, scale,&
       row_size, col_size, lb, lb2, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dst
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: scale
    INTEGER, INTENT(IN)                      :: row_size, col_size
    INTEGER, INTENT(IN), OPTIONAL            :: lb, lb2
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_block_scale', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_size, error_handler, &
                                                lb2_t, lb_t, ub_t

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handler, error)
    IF (debug_mod) THEN
       CALL dbcsr_assert (ASSOCIATED(dst%d),&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Data area must be setup.",__LINE__,error)
       CALL dbcsr_assert (dst%d%data_type.EQ.dbcsr_type_real_8&
            .OR.dst%d%data_type.EQ.dbcsr_type_real_4&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_8&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_4&
            .OR.dst%d%data_type.EQ.dbcsr_type_real_8_2d&
            .OR.dst%d%data_type.EQ.dbcsr_type_real_4_2d&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_8_2d&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_4_2d,&
            dbcsr_warning_level, dbcsr_caller_error, routineN, "Incorrect data type.",__LINE__,error)
    ENDIF
    CALL dbcsr_assert (scale%data_type, &
         'EQ', dbcsr_type_2d_to_1d (dst%d%data_type),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Incompatible data types",__LINE__,error)

    data_size = row_size*col_size
    lb_t = 1
    IF (PRESENT (lb)) lb_t = lb
    ub_t = lb_t + data_size - 1
    lb2_t = 1
    IF (PRESENT (lb2)) lb2_t = lb2
    SELECT CASE (dst%d%data_type)
    CASE (dbcsr_type_real_8)
       dst%d%r_dp(lb_t:ub_t) = dst%d%r_dp(lb_t:ub_t) * scale%r_dp
    CASE (dbcsr_type_real_4)
       dst%d%r_sp(lb_t:ub_t) = dst%d%r_sp(lb_t:ub_t) * scale%r_sp
    CASE (dbcsr_type_complex_8)
       dst%d%c_dp(lb_t:ub_t) = dst%d%c_dp(lb_t:ub_t) * scale%c_dp
    CASE (dbcsr_type_complex_4)
       dst%d%c_sp(lb_t:ub_t) = dst%d%c_sp(lb_t:ub_t) * scale%c_sp
    CASE (dbcsr_type_real_8_2d)
       dst%d%r2_dp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) =&
       dst%d%r2_dp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) * scale%r_dp
    CASE (dbcsr_type_real_4_2d)
       dst%d%r2_sp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) =&
       dst%d%r2_sp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) * scale%r_sp
    CASE (dbcsr_type_complex_8_2d)
       dst%d%c2_dp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) =&
       dst%d%c2_dp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) * scale%c_dp
    CASE (dbcsr_type_complex_4_2d)
       dst%d%c2_sp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) =&
       dst%d%c2_sp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) * scale%c_sp
    CASE default
       CALL dbcsr_assert (.FALSE.,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN, "Incorrect data type.",__LINE__,error)
    END SELECT
    IF (careful_mod) CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_block_scale

! *****************************************************************************
!> \brief Negates the real part of a block
!>
!> There are no checks done for correctness!
!> \param[in] src        source data area
!> \param[in] row_size   row size of existing block
!> \param[in] col_size   column size of existing block
!> \param[in] lb         (optional) lower bound for destination (and source if
!>                       not given explicity)
!> \param[in] source_lb  (optional) lower bound of source
!> \param[in] scale      (optional) scale data
!> \param[in] lb2        (optional) lower bound of 2nd dimension for target
!> \param[in] source_lb2 (optional) lower bound of 2nd dimension for source
!> \par History
!> - 2010-09 [??] Copied from block_transpose?
!> - 2010-09-20 [UB] Swaps/corrects row/column definitions for 2-D bounds
! *****************************************************************************
  SUBROUTINE dbcsr_block_real_neg (dst,&
       row_size, col_size, lb, lb2, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dst
    INTEGER, INTENT(IN)                      :: row_size, col_size
    INTEGER, INTENT(IN), OPTIONAL            :: lb, lb2
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_block_real_neg', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_size, error_handler, &
                                                lb2_t, lb_t, ub_t

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handler, error)
    IF (debug_mod) THEN
       CALL dbcsr_assert (ASSOCIATED(dst%d),&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Data area must be setup.",__LINE__,error)
       CALL dbcsr_assert (dst%d%data_type.EQ.dbcsr_type_real_8&
            .OR.dst%d%data_type.EQ.dbcsr_type_real_4&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_8&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_4&
            .OR.dst%d%data_type.EQ.dbcsr_type_real_8_2d&
            .OR.dst%d%data_type.EQ.dbcsr_type_real_4_2d&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_8_2d&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_4_2d,&
            dbcsr_warning_level, dbcsr_caller_error, routineN, "Incorrect data type.",__LINE__,error)
    ENDIF

    data_size = row_size*col_size
    lb_t = 1
    IF (PRESENT (lb)) lb_t = lb
    ub_t = lb_t + data_size - 1
    lb2_t = 1
    IF (PRESENT (lb2)) lb2_t = lb2
    SELECT CASE (dst%d%data_type)
    CASE (dbcsr_type_real_8)
       dst%d%r_dp(lb_t:ub_t) = -dst%d%r_dp(lb_t:ub_t)
    CASE (dbcsr_type_real_4)
       dst%d%r_sp(lb_t:ub_t) = -dst%d%r_sp(lb_t:ub_t)
    CASE (dbcsr_type_complex_8)
       dst%d%c_dp(lb_t:ub_t) = CMPLX (&
            REAL(dst%d%c_dp(lb_t:ub_t), KIND=real_8),&
            AIMAG(dst%d%c_dp(lb_t:ub_t)),&
            KIND=real_8)
    CASE (dbcsr_type_complex_4)
       dst%d%c_sp(lb_t:ub_t) = CMPLX(&
            REAL(dst%d%c_sp(lb_t:ub_t), KIND=real_4),&
            AIMAG(dst%d%c_sp(lb_t:ub_t)),&
            KIND=real_4)
    CASE (dbcsr_type_real_8_2d)
       dst%d%r2_dp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) =&
      -dst%d%r2_dp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1)
    CASE (dbcsr_type_real_4_2d)
       dst%d%r2_sp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) =&
      -dst%d%r2_sp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1)
    CASE (dbcsr_type_complex_8_2d)
       dst%d%c2_dp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) =&
            CMPLX(&
            REAL(dst%d%c2_dp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1), KIND=real_8),&
            AIMAG(dst%d%c2_dp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1)),&
            KIND=real_8)
    CASE (dbcsr_type_complex_4_2d)
       dst%d%c2_sp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) =&
            CMPLX(&
            REAL(dst%d%c2_sp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1), KIND=real_4),&
            AIMAG(dst%d%c2_sp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1)),&
            KIND=real_4)
    CASE default
       CALL dbcsr_assert (.FALSE.,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN, &
            "Incorrect data type.",__LINE__,error)
    END SELECT
    IF (careful_mod) CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_block_real_neg

! *****************************************************************************
!> \brief Conjugate a data area.
!>
!> There are no checks done for correctness!
!> \param[in] src        source data area
!> \param[in] row_size   row size of existing block
!> \param[in] col_size   column size of existing block
!> \param[in] lb         (optional) lower bound for destination (and source if
!>                       not given explicity)
!> \param[in] source_lb  (optional) lower bound of source
!> \param[in] scale      (optional) scale data
!> \param[in] lb2        (optional) lower bound of 2nd dimension for target
!> \param[in] source_lb2 (optional) lower bound of 2nd dimension for source
!> \par History
!> - 2010-09 [??] Copied from block_transpose?
!> - 2010-09-20 [UB] Swaps/corrects row/column definitions for 2-D bounds
! *****************************************************************************
  SUBROUTINE dbcsr_block_conjg (dst,&
       row_size, col_size, lb, lb2, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dst
    INTEGER, INTENT(IN)                      :: row_size, col_size
    INTEGER, INTENT(IN), OPTIONAL            :: lb, lb2
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_block_conjg', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_size, error_handler, &
                                                lb2_t, lb_t, ub_t

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handler, error)
    IF (debug_mod) THEN
       CALL dbcsr_assert (ASSOCIATED(dst%d),&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Data area must be setup.",__LINE__,error)
       CALL dbcsr_assert (dst%d%data_type.EQ.dbcsr_type_real_8&
            .OR.dst%d%data_type.EQ.dbcsr_type_real_4&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_8&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_4&
            .OR.dst%d%data_type.EQ.dbcsr_type_real_8_2d&
            .OR.dst%d%data_type.EQ.dbcsr_type_real_4_2d&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_8_2d&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_4_2d,&
            dbcsr_warning_level, dbcsr_caller_error, routineN, "Incorrect data type.",__LINE__,error)
    ENDIF

    data_size = row_size*col_size
    lb_t = 1
    IF (PRESENT (lb)) lb_t = lb
    ub_t = lb_t + data_size - 1
    lb2_t = 1
    IF (PRESENT (lb2)) lb2_t = lb2
    SELECT CASE (dst%d%data_type)
    CASE (dbcsr_type_real_8)
       dst%d%r_dp(lb_t:ub_t) = dst%d%r_dp(lb_t:ub_t)
    CASE (dbcsr_type_real_4)
       dst%d%r_sp(lb_t:ub_t) = dst%d%r_sp(lb_t:ub_t)
    CASE (dbcsr_type_complex_8)
       dst%d%c_dp(lb_t:ub_t) = CONJG(dst%d%c_dp(lb_t:ub_t))
    CASE (dbcsr_type_complex_4)
       dst%d%c_sp(lb_t:ub_t) = CONJG(dst%d%c_sp(lb_t:ub_t))
    CASE (dbcsr_type_real_8_2d)
       dst%d%r2_dp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) =&
       dst%d%r2_dp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1)
    CASE (dbcsr_type_real_4_2d)
       dst%d%r2_sp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) =&
       dst%d%r2_sp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1)
    CASE (dbcsr_type_complex_8_2d)
       dst%d%c2_dp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) =&
       CONJG(dst%d%c2_dp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1))
    CASE (dbcsr_type_complex_4_2d)
       dst%d%c2_sp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1) =&
       CONJG(dst%d%c2_sp(lb_t:lb_t+row_size-1,lb2_t:lb2_t+col_size-1))
    CASE default
       CALL dbcsr_assert (.FALSE.,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN, &
            "Incorrect data type.",__LINE__,error)
    END SELECT
    IF (careful_mod) CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_block_conjg

! *****************************************************************************
!> \brief In-place transpose of encapsulated data
!>
!> There are no checks done for correctness!
!> \param[in] area       encapsulated data area
!> \param[in] row_size   number of rows in existing block
!> \param[in] col_size   number of columns in existing block
! *****************************************************************************
  SUBROUTINE dbcsr_block_transpose_a (area, row_size, col_size, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, INTENT(IN)                      :: row_size, col_size
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_block_transpose_a', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

!   ---------------------------------------------------------------------------

    IF (careful_mod) &
         CALL dbcsr_error_set(routineN, error_handler, error)
    SELECT CASE (area%d%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_block_transpose (area%d%r_dp,&
            row_size, col_size)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_block_transpose (area%d%r_sp,&
            row_size, col_size)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_block_transpose (area%d%c_dp,&
            row_size, col_size)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_block_transpose (area%d%c_sp,&
            row_size, col_size)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
    IF (careful_mod) &
         CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_block_transpose_a


! *****************************************************************************
!> \brief Copy data from one data area to another.
!>
!> There are no checks done for correctness!
!> \param[in] dst        destination data area
!> \param[in] lb         lower bound for destination (and source if
!>                       not given explicity)
!> \param[in] data_size  number of elements to copy
!> \param[in] src        source data area
!> \param[in] source_lb  (optional) lower bound of source
!> \param[in] scale      (optional) scale by this factor
!> \param[in] lb2        (optional) 2nd dimension lower bound
!> \param[in] data_size2 (optional) 2nd dimension data size
!> \param[in] source_lb2 (optional) 2nd dimension lower bound for source
! *****************************************************************************
  SUBROUTINE dbcsr_data_copy_aa (dst, lb, data_size, src, source_lb, scale,&
       lb2, data_size2, source_lb2)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dst
    INTEGER, INTENT(IN)                      :: lb, data_size
    TYPE(dbcsr_data_obj), INTENT(IN)         :: src
    INTEGER, INTENT(IN), OPTIONAL            :: source_lb
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: scale
    INTEGER, INTENT(IN), OPTIONAL            :: lb2, data_size2, source_lb2

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_copy_aa', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lb2_s, lb_s, ub, ub2, ub2_s, &
                                                ub_s
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (debug_mod) THEN
       CALL dbcsr_assert (ASSOCIATED(dst%d) .AND. ASSOCIATED(src%d),&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Data areas must be setup.",__LINE__,error)
       CALL dbcsr_assert (dst%d%data_type .EQ. src%d%data_type, dbcsr_fatal_level,&
            dbcsr_caller_error, routineN,&
            "Data type must be the same.",__LINE__,error)
       !CALL dbcsr_assert (dst%d%data_type.EQ.dbcsr_type_real_8&
       !     .OR.dst%d%data_type.EQ.dbcsr_type_real_4&
       !     .OR.dst%d%data_type.EQ.dbcsr_type_complex_8&
       !     .OR.dst%d%data_type.EQ.dbcsr_type_complex_4,&
       !     dbcsr_warning_level, dbcsr_caller_error, routineN,&
       !     "Incorrect data types: "//dst%d%data_type//", "//src%d%data_type)
    ENDIF
    IF (PRESENT (scale) .AND. careful_mod) THEN
       IF (dbcsr_type_is_2d(src%d%data_type)) THEN
          CALL dbcsr_assert (scale%data_type&
               .EQ. dbcsr_type_2d_to_1d (src%d%data_type),&
               dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
               "Incomptable data types",__LINE__,error)
       ELSE
          CALL dbcsr_assert (scale%data_type .EQ. src%d%data_type, dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN, "Incomptable data types",__LINE__,error)
       ENDIF
    ENDIF
    ub = lb + data_size - 1
    IF (PRESENT (source_lb)) THEN
       lb_s = source_lb
       ub_s = source_lb + data_size-1
    ELSE
       lb_s = lb
       ub_s = ub
    ENDIF
    IF (careful_mod) THEN
      CALL dbcsr_assert (.NOT.dbcsr_type_is_2d (src%d%data_type)&
           .OR. PRESENT (lb2),&
           dbcsr_fatal_level, dbcsr_caller_error, routineN,&
           "Must specify lb2 for 2-D data area",__LINE__,error)
      CALL dbcsr_assert (.NOT.dbcsr_type_is_2d (src%d%data_type)&
           .OR. PRESENT (data_size2),&
           dbcsr_fatal_level, dbcsr_caller_error, routineN,&
           "Must specify data_size2 for 2-D data area",__LINE__,error)
    ENDIF
    IF (PRESENT (lb2)) THEN
       IF (careful_mod) THEN
          CALL dbcsr_assert (dbcsr_type_is_2d (src%d%data_type),&
               dbcsr_warning_level, dbcsr_caller_error, routineN,&
               "Should not specify lb2 for 1-D data",__LINE__,error)
       ENDIF
       ub2 = lb2 + data_size2 - 1
       IF (PRESENT (source_lb2)) THEN
          lb2_s = source_lb2
          ub2_s = source_lb2 + data_size2 - 1
       ELSE
          lb2_s = lb2
          ub2_s = ub2
       ENDIF
       !write(*,*)routineN//" lb,ub2D <- S",lb2, ub2, lb2_s, ub2_s
    ENDIF
    SELECT CASE (src%d%data_type)
    CASE (dbcsr_type_real_4)
       IF (debug_mod) THEN
          CALL dbcsr_assert(ASSOCIATED(dst%d%r_sp),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "associated(dst%d%r_sp)",__LINE__,error)
          CALL dbcsr_assert(ASSOCIATED(src%d%r_sp),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "associated(src%d%r_sp)",__LINE__,error)
          CALL dbcsr_assert(lb, "GE", LBOUND(dst%d%r_sp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lb dst%d%r_sp",__LINE__,error)
          CALL dbcsr_assert(ub, "LE", UBOUND(dst%d%r_sp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "ub dst%d%r_sp",__LINE__,error)
          CALL dbcsr_assert(lb_s, "GE", LBOUND(src%d%r_sp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lb src%d%r_sp",__LINE__,error)
          CALL dbcsr_assert(ub_s, "LE", UBOUND(src%d%r_sp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "ub src%d%r_sp",__LINE__,error)
       ENDIF
       IF (PRESENT (scale)) THEN
          dst%d%r_sp(lb:ub) = scale%r_sp * src%d%r_sp(lb_s:ub_s)
       ELSE
          dst%d%r_sp(lb:ub) = src%d%r_sp(lb_s:ub_s)
       ENDIF
    CASE (dbcsr_type_real_8)
       IF (debug_mod) THEN
          CALL dbcsr_assert(ASSOCIATED(dst%d%r_dp),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "associated(dst%d%r_dp)",__LINE__,error)
          CALL dbcsr_assert(ASSOCIATED(src%d%r_dp),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "associated(src%d%r_dp)",__LINE__,error)
          CALL dbcsr_assert(lb, "GE", LBOUND(dst%d%r_dp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lb dst%d%r_dp",__LINE__,error)
          CALL dbcsr_assert(ub, "LE", UBOUND(dst%d%r_dp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "ub dst%d%r_dp",__LINE__,error)
          CALL dbcsr_assert(lb_s, "GE", LBOUND(src%d%r_dp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lb src%d%r_dp",__LINE__,error)
          CALL dbcsr_assert(ub_s, "LE", UBOUND(src%d%r_dp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "ub src%d%r_dp",__LINE__,error)
       ENDIF
       IF (PRESENT (scale)) THEN
          dst%d%r_dp(lb:ub) = scale%r_dp * src%d%r_dp(lb_s:ub_s)
       ELSE
          dst%d%r_dp(lb:ub) = src%d%r_dp(lb_s:ub_s)
       ENDIF
    CASE (dbcsr_type_complex_4)
       IF (debug_mod) THEN
          CALL dbcsr_assert(ASSOCIATED(dst%d%c_sp),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "associated(dst%d%c_sp)",__LINE__,error)
          CALL dbcsr_assert(ASSOCIATED(src%d%c_sp),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "associated(src%d%c_sp)",__LINE__,error)
          CALL dbcsr_assert(lb, "GE", LBOUND(dst%d%c_sp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lb dst%d%c_sp",__LINE__,error)
          CALL dbcsr_assert(ub, "LE", UBOUND(dst%d%c_sp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "ub dst%d%c_sp",__LINE__,error)
          CALL dbcsr_assert(lb_s, "GE", LBOUND(src%d%c_sp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lb src%d%c_sp",__LINE__,error)
          CALL dbcsr_assert(ub_s, "LE", UBOUND(src%d%c_sp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "ub src%d%c_sp",__LINE__,error)
       ENDIF
       IF (PRESENT (scale)) THEN
          dst%d%c_sp(lb:ub) = scale%c_sp * src%d%c_sp(lb_s:ub_s)
       ELSE
          dst%d%c_sp(lb:ub) = src%d%c_sp(lb_s:ub_s)
       ENDIF
    CASE (dbcsr_type_complex_8)
       IF (debug_mod) THEN
          CALL dbcsr_assert(ASSOCIATED(dst%d%c_dp),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "associated(dst%d%c_dp)",__LINE__,error)
          CALL dbcsr_assert(ASSOCIATED(src%d%c_dp),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "associated(src%d%c_dp)",__LINE__,error)
          CALL dbcsr_assert(lb, "GE", LBOUND(dst%d%c_dp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lb dst%d%c_dp",__LINE__,error)
          CALL dbcsr_assert(ub, "LE", UBOUND(dst%d%c_dp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "ub dst%d%c_dp",__LINE__,error)
          CALL dbcsr_assert(lb_s, "GE", LBOUND(src%d%c_dp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lb src%d%c_dp",__LINE__,error)
          CALL dbcsr_assert(ub_s, "LE", UBOUND(src%d%c_dp,1),dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "ub src%d%c_dp",__LINE__,error)
       ENDIF
       IF (PRESENT (scale)) THEN
          dst%d%c_dp(lb:ub) = scale%c_dp * src%d%c_dp(lb_s:ub_s)
       ELSE
          dst%d%c_dp(lb:ub) = src%d%c_dp(lb_s:ub_s)
       ENDIF
    CASE (dbcsr_type_real_4_2d)
       IF (PRESENT (scale)) THEN
          dst%d%r2_sp(lb:ub, lb2:ub2) =&
               scale%r_sp * src%d%r2_sp(lb_s:ub_s, lb2_s:ub2_s)
       ELSE
          dst%d%r2_sp(lb:ub, lb2:ub2) = src%d%r2_sp(lb_s:ub_s, lb2_s:ub2_s)
       ENDIF
    CASE (dbcsr_type_real_8_2d)
       IF (PRESENT (scale)) THEN
          dst%d%r2_dp(lb:ub, lb2:ub2) =&
               scale%r_dp * src%d%r2_dp(lb_s:ub_s, lb2_s:ub2_s)
       ELSE
          dst%d%r2_dp(lb:ub, lb2:ub2) = src%d%r2_dp(lb_s:ub_s, lb2_s:ub2_s)
       ENDIF
    CASE (dbcsr_type_complex_4_2d)
       IF (PRESENT (scale)) THEN
          dst%d%c2_sp(lb:ub, lb2:ub2) =&
               scale%c_sp * src%d%c2_sp(lb_s:ub_s, lb2_s:ub2_s)
       ELSE
          dst%d%c2_sp(lb:ub, lb2:ub2) = src%d%c2_sp(lb_s:ub_s, lb2_s:ub2_s)
       ENDIF
    CASE (dbcsr_type_complex_8_2d)
       IF (PRESENT (scale)) THEN
          dst%d%c2_dp(lb:ub, lb2:ub2) =&
               scale%c_dp * src%d%c2_dp(lb_s:ub_s, lb2_s:ub2_s)
       ELSE
          dst%d%c2_dp(lb:ub, lb2:ub2) = src%d%c2_dp(lb_s:ub_s, lb2_s:ub2_s)
       ENDIF
    CASE default
       CALL dbcsr_assert(.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
  END SUBROUTINE dbcsr_data_copy_aa


! *****************************************************************************
!> \brief Copy data from one data area to another, the most basic form.
!>
!> There are no checks done for correctness!
!> \param[in,out] dst    destination data area
!> \param[in] dst_lb     lower bounds for destination
!> \param[in] dst_sizes  sizes for destination
!> \param[in] src        source data area
!> \param[in] src_lb     lower bounds for source
!> \param[in] src_sizes  sizes for source
!> \param[in,out] error  error
! *****************************************************************************
  SUBROUTINE dbcsr_data_copy_aa2 (dst, dst_lb, dst_sizes,&
       src, src_lb, src_sizes, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dst
    INTEGER, DIMENSION(:), INTENT(IN)        :: dst_lb, dst_sizes
    TYPE(dbcsr_data_obj), INTENT(IN)         :: src
    INTEGER, DIMENSION(:), INTENT(IN)        :: src_lb, src_sizes
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_copy_aa2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dst_d, dst_dt, error_handler, &
                                                src_d, src_dt
    INTEGER, DIMENSION(2)                    :: dst_ub, src_ub

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    !
    src_dt = dbcsr_data_get_type (src)
    dst_dt = dbcsr_data_get_type (dst)
    IF (dbcsr_type_is_2d (src_dt)) THEN
       src_d = 2
    ELSE
       src_d = 1
    ENDIF
    IF (dbcsr_type_is_2d (dst_dt)) THEN
       dst_d = 2
    ELSE
       dst_d = 1
    ENDIF
    src_ub(1:src_d) = src_lb(1:src_d) + src_sizes(1:src_d) - 1
    dst_ub(1:dst_d) = dst_lb(1:dst_d) + dst_sizes(1:dst_d) - 1
    IF (careful_mod) THEN
       CALL dbcsr_assert (dbcsr_data_exists (dst, error), dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN, "Invalid target data area",&
            __LINE__, error=error)
       CALL dbcsr_assert (dbcsr_data_exists (src, error), dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN, "Invalid source data area",&
            __LINE__, error=error)
       CALL dbcsr_assert (dbcsr_type_2d_to_1d (src_dt), "EQ",&
            dbcsr_type_2d_to_1d (dst_dt), dbcsr_fatal_level, &
            dbcsr_caller_error, routineN,&
            "Data types must be comptable: ", __LINE__, error=error)
       IF (dbcsr_type_is_2d (dst_dt)) THEN
          CALL dbcsr_assert (SIZE (dst_lb), "EQ", 2, dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "size must be 2 for 2-d dst_lb", __LINE__, error=error)
          CALL dbcsr_assert (SIZE (dst_sizes), "EQ", 2, dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "size must be 2 for 2-d dst_sizes", __LINE__, error=error)
       ELSE
          CALL dbcsr_assert (SIZE (dst_lb), "EQ", 1, dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "size must be 1 for 1-d dst_lb", __LINE__, error=error)
          CALL dbcsr_assert (SIZE (dst_sizes), "EQ", 1, dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "size must be 1 for 1-d dst_sizes", __LINE__, error=error)
       ENDIF
       IF (dbcsr_type_is_2d (src_dt)) THEN
          CALL dbcsr_assert (SIZE (src_lb), "EQ", 2, dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "size must be 2 for 2-d src_lb", __LINE__, error=error)
          CALL dbcsr_assert (SIZE (src_sizes), "EQ", 2, dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "size must be 2 for 2-d src_sizes", __LINE__, error=error)
       ELSE
          CALL dbcsr_assert (SIZE (src_lb), "EQ", 1, dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "size must be 1 for 1-d src_lb", __LINE__, error=error)
          CALL dbcsr_assert (SIZE (src_sizes), "EQ", 1, dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "size must be 1 for 1-d src_sizes", __LINE__, error=error)
       ENDIF
       IF (debug_mod) THEN
          CALL dbcsr_data_verify_bounds (dst, dst_lb(1:dst_d), dst_ub(1:dst_d),&
               error)
          CALL dbcsr_data_verify_bounds (src, src_lb(1:src_d), src_ub(1:src_d),&
               error)
       ENDIF
    ENDIF
    !
    SELECT CASE (src_dt)
    CASE (dbcsr_type_real_4)
       IF (dbcsr_type_is_2d (dst_dt)) THEN
          CALL dbcsr_block_copy (dst%d%r2_sp(dst_lb(1):dst_ub(1),&
                                             dst_lb(2):dst_ub(2)),&
                                 src%d%r_sp(src_lb(1):src_ub(1)),&
                                 src_sizes(1), 1)
       ELSE
          CALL dbcsr_block_copy (dst%d%r_sp(dst_lb(1):dst_ub(1)),&
                                 src%d%r_sp(src_lb(1):src_ub(1)),&
                                 src_sizes(1), 1)
       ENDIF
    CASE (dbcsr_type_real_8)
       IF (dbcsr_type_is_2d (dst_dt)) THEN
          CALL dbcsr_block_copy (dst%d%r2_dp(dst_lb(1):dst_ub(1),&
                                             dst_lb(2):dst_ub(2)),&
                                 src%d%r_dp(src_lb(1):src_ub(1)),&
                                 src_sizes(1), 1)
       ELSE
          CALL dbcsr_block_copy (dst%d%r_dp(dst_lb(1):dst_ub(1)),&
                                 src%d%r_dp(src_lb(1):src_ub(1)),&
                                 src_sizes(1), 1)
       ENDIF
    CASE (dbcsr_type_complex_4)
       IF (dbcsr_type_is_2d (dst_dt)) THEN
          CALL dbcsr_block_copy (dst%d%c2_sp(dst_lb(1):dst_ub(1),&
                                             dst_lb(2):dst_ub(2)),&
                                 src%d%c_sp(src_lb(1):src_ub(1)),&
                                 src_sizes(1), 1)
       ELSE
          CALL dbcsr_block_copy (dst%d%c_sp(dst_lb(1):dst_ub(1)),&
                                 src%d%c_sp(src_lb(1):src_ub(1)),&
                                 src_sizes(1), 1)
       ENDIF
    CASE (dbcsr_type_complex_8)
       IF (dbcsr_type_is_2d (dst_dt)) THEN
          CALL dbcsr_block_copy (dst%d%c2_dp(dst_lb(1):dst_ub(1),&
                                             dst_lb(2):dst_ub(2)),&
                                 src%d%c_dp(src_lb(1):src_ub(1)),&
                                 src_sizes(1), 1)
       ELSE
          CALL dbcsr_block_copy (dst%d%c_dp(dst_lb(1):dst_ub(1)),&
                                 src%d%c_dp(src_lb(1):src_ub(1)),&
                                 src_sizes(1), 1)
       ENDIF
    CASE (dbcsr_type_real_4_2d)
       IF (dbcsr_type_is_2d (dst_dt)) THEN
          CALL dbcsr_block_copy (dst%d%r2_sp(dst_lb(1):dst_ub(1),&
                                             dst_lb(2):dst_ub(2)),&
                                 src%d%r2_sp(src_lb(1):src_ub(1),&
                                             src_lb(2):src_ub(2)),&
                                 dst_sizes(1), dst_sizes(2))
       ELSE
          CALL dbcsr_block_copy (dst%d%r_sp(dst_lb(1):dst_ub(1)),&
                                 src%d%r2_sp(src_lb(1):src_ub(1),&
                                             src_lb(2):src_ub(2)),&
                                 dst_sizes(1), dst_sizes(2))
       ENDIF
    CASE (dbcsr_type_real_8_2d)
       IF (dbcsr_type_is_2d (dst_dt)) THEN
          CALL dbcsr_block_copy (dst%d%r2_dp(dst_lb(1):dst_ub(1),&
                                             dst_lb(2):dst_ub(2)),&
                                 src%d%r2_dp(src_lb(1):src_ub(1),&
                                             src_lb(2):src_ub(2)),&
                                 dst_sizes(1), dst_sizes(2))
       ELSE
          CALL dbcsr_block_copy (dst%d%r_dp(dst_lb(1):dst_ub(1)),&
                                 src%d%r2_dp(src_lb(1):src_ub(1),&
                                             src_lb(2):src_ub(2)),&
                                 dst_sizes(1), dst_sizes(2))
       ENDIF
    CASE (dbcsr_type_complex_4_2d)
       IF (dbcsr_type_is_2d (dst_dt)) THEN
          CALL dbcsr_block_copy (dst%d%c2_sp(dst_lb(1):dst_ub(1),&
                                             dst_lb(2):dst_ub(2)),&
                                 src%d%c2_sp(src_lb(1):src_ub(1),&
                                             src_lb(2):src_ub(2)),&
                                 dst_sizes(1), dst_sizes(2))
       ELSE
          CALL dbcsr_block_copy (dst%d%c_sp(dst_lb(1):dst_ub(1)),&
                                 src%d%c2_sp(src_lb(1):src_ub(1),&
                                             src_lb(2):src_ub(2)),&
                                 dst_sizes(1), dst_sizes(2))
       ENDIF
    CASE (dbcsr_type_complex_8_2d)
       IF (dbcsr_type_is_2d (dst_dt)) THEN
          CALL dbcsr_block_copy (dst%d%c2_dp(dst_lb(1):dst_ub(1),&
                                             dst_lb(2):dst_ub(2)),&
                                 src%d%c2_dp(src_lb(1):src_ub(1),&
                                             src_lb(2):src_ub(2)),&
                                 dst_sizes(1), dst_sizes(2))
       ELSE
          CALL dbcsr_block_copy (dst%d%c_dp(dst_lb(1):dst_ub(1)),&
                                 src%d%c2_dp(src_lb(1):src_ub(1),&
                                             src_lb(2):src_ub(2)),&
                                 dst_sizes(1), dst_sizes(2))
       ENDIF
    CASE default
       CALL dbcsr_assert(.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_data_copy_aa2


! *****************************************************************************
!> \brief Clears a data area, possibly transposed.
! *****************************************************************************
  SUBROUTINE dbcsr_data_clear_nt (area, lb, ub, value, lb2, ub2, tr)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, INTENT(IN), OPTIONAL            :: lb, ub
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: value
    INTEGER, INTENT(IN), OPTIONAL            :: lb2, ub2
    LOGICAL, INTENT(in)                      :: tr

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_clear_nt', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF (tr) THEN
       CALL dbcsr_data_clear0 (area, lb=lb2, ub=ub2, value=value, lb2=lb, ub2=ub)
    ELSE
       CALL dbcsr_data_clear0 (area, lb=lb, ub=ub, value=value, lb2=lb2, ub2=ub2)
    ENDIF
  END SUBROUTINE dbcsr_data_clear_nt

! *****************************************************************************
!> \brief Clears a data area
!> \param[in,out] area   area with encapsulated data
!> \param[in] lb, lb2    (optional) lower bound for clearing
!> \param[in] ub, ub2    (optional) upper bound for clearing
!> \param[in] value      (optional) value to use for clearing
! *****************************************************************************
  SUBROUTINE dbcsr_data_clear0 (area, lb, ub, value, lb2, ub2)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, INTENT(IN), OPTIONAL            :: lb, ub
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: value
    INTEGER, INTENT(IN), OPTIONAL            :: lb2, ub2

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_clear0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: l, l2, s, u, u2
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------
! CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (ASSOCIATED(area%d),&
         dbcsr_fatal_level, dbcsr_caller_error, routineN,&
         "Data area must be setup.",__LINE__,error)
    IF (PRESENT (value)) THEN
       CALL dbcsr_assert (area%d%data_type .EQ. value%data_type,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Incompatible data types",__LINE__,error)
    ENDIF

    IF(acc_devmem_allocated(area%d%acc_devmem)) THEN
      IF (PRESENT(value)) STOP "dbcsr_data_clear0 with value not implemented for acc_devmem"
      s = dbcsr_datatype_sizeof(area%d%data_type)
      CALL acc_devmem_setzero_bytes(area%d%acc_devmem, s*lb, s*ub, area%d%memory_type%acc_stream)
    END IF

    SELECT CASE (area%d%data_type)
    CASE (dbcsr_type_real_4)
       l = LBOUND(area%d%r_sp, 1)
       u = UBOUND(area%d%r_sp, 1)
       IF (PRESENT (lb)) THEN
          CALL dbcsr_assert(lb, "GE", l, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lower bound too low",__LINE__,error)
          l = lb
       ENDIF
       IF (PRESENT (ub)) THEN
          CALL dbcsr_assert(ub, "LE", u, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "upper bound too high",__LINE__,error)
          u = ub
       ENDIF
       IF (PRESENT (value)) THEN
          area%d%r_sp(l:u) = value%r_sp
       ELSE
          area%d%r_sp(l:u) = 0.0_real_4
       ENDIF
    CASE (dbcsr_type_real_8)
       l = LBOUND(area%d%r_dp, 1)
       u = UBOUND(area%d%r_dp, 1)
       IF (PRESENT (lb)) THEN
          CALL dbcsr_assert(lb, "GE", l, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lower bound too low",__LINE__,error)
          l = lb
       ENDIF
       IF (PRESENT (ub)) THEN
          CALL dbcsr_assert(ub, "LE", u, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "upper bound too high",__LINE__,error)
          u = ub
       ENDIF
       IF (PRESENT (value)) THEN
          area%d%r_dp(l:u) = value%r_dp
       ELSE
          area%d%r_dp(l:u) = 0.0_real_8
       ENDIF
    CASE (dbcsr_type_complex_4)
       l = LBOUND(area%d%c_sp, 1)
       u = UBOUND(area%d%c_sp, 1)
       IF (PRESENT (lb)) THEN
          CALL dbcsr_assert(lb, "GE", l, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lower bound too low",__LINE__,error)
          l = lb
       ENDIF
       IF (PRESENT (ub)) THEN
          CALL dbcsr_assert(ub, "LE", u, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "upper bound too high",__LINE__,error)
          u = ub
       ENDIF
       IF (PRESENT (value)) THEN
          area%d%c_sp(l:u) = value%c_sp
       ELSE
          area%d%c_sp(l:u) = CMPLX(0.0, 0.0, real_4)
       ENDIF
    CASE (dbcsr_type_complex_8)
       l = LBOUND(area%d%c_dp, 1)
       u = UBOUND(area%d%c_dp, 1)
       IF (PRESENT (lb)) THEN
          CALL dbcsr_assert(lb, "GE", l, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lower bound too low",__LINE__,error)
          l = lb
       ENDIF
       IF (PRESENT (ub)) THEN
          CALL dbcsr_assert(ub, "LE", u, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "upper bound too high",__LINE__,error)
          u = ub
       ENDIF
       IF (PRESENT (value)) THEN
          area%d%c_dp(l:u) = value%c_dp
       ELSE
          area%d%c_dp(l:u) = CMPLX(0.0, 0.0, real_8)
       ENDIF
    CASE (dbcsr_type_real_4_2d)
       l = LBOUND(area%d%r2_sp, 1)
       u = UBOUND(area%d%r2_sp, 1)
       l2 = LBOUND(area%d%r2_sp, 2)
       u2 = UBOUND(area%d%r2_sp, 2)
       IF (PRESENT (lb)) THEN
          CALL dbcsr_assert(lb, "GE", l, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lower bound too low",__LINE__,error)
          l = lb
       ENDIF
       IF (PRESENT (ub)) THEN
          CALL dbcsr_assert(ub, "LE", u, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "upper bound too high",__LINE__,error)
          u = ub
       ENDIF
       IF (PRESENT (lb2)) THEN
          CALL dbcsr_assert(lb2, "GE", l2, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lower2 bound too low",__LINE__,error)
          l2 = lb2
       ENDIF
       IF (PRESENT (ub2)) THEN
          CALL dbcsr_assert(ub2, "LE", u2, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "upper2 bound too high",__LINE__,error)
          u2 = ub2
       ENDIF
       IF (PRESENT (value)) THEN
          area%d%r2_sp(l:u, l2:u2) = value%r_sp
       ELSE
          area%d%r2_sp(l:u, l2:u2) = 0.0_real_4
       ENDIF
    CASE (dbcsr_type_real_8_2d)
       l = LBOUND(area%d%r2_dp, 1)
       u = UBOUND(area%d%r2_dp, 1)
       l2 = LBOUND(area%d%r2_dp, 2)
       u2 = UBOUND(area%d%r2_dp, 2)
       IF (PRESENT (lb)) THEN
          CALL dbcsr_assert(lb, "GE", l, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lower bound too low",__LINE__,error)
          l = lb
       ENDIF
       IF (PRESENT (ub)) THEN
          CALL dbcsr_assert(ub, "LE", u, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "upper bound too high",__LINE__,error)
          u = ub
       ENDIF
       IF (PRESENT (lb2)) THEN
          CALL dbcsr_assert(lb2, "GE", l2, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lower2 bound too low",__LINE__,error)
          l2 = lb2
       ENDIF
       IF (PRESENT (ub2)) THEN
          CALL dbcsr_assert(ub2, "LE", u2, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "upper2 bound too high",__LINE__,error)
          u2 = ub2
       ENDIF
       IF (PRESENT (value)) THEN
          area%d%r2_dp(l:u, l2:u2) = value%r_dp
       ELSE
          area%d%r2_dp(l:u, l2:u2) = 0.0_real_8
       ENDIF
    CASE (dbcsr_type_complex_4_2d)
       l = LBOUND(area%d%c2_sp, 1)
       u = UBOUND(area%d%c2_sp, 1)
       l2 = LBOUND(area%d%c2_sp, 2)
       u2 = UBOUND(area%d%c2_sp, 2)
       IF (PRESENT (lb)) THEN
          CALL dbcsr_assert(lb, "GE", l, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lower bound too low",__LINE__,error)
          l = lb
       ENDIF
       IF (PRESENT (ub)) THEN
          CALL dbcsr_assert(ub, "LE", u, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "upper bound too high",__LINE__,error)
          u = ub
       ENDIF
       IF (PRESENT (lb2)) THEN
          CALL dbcsr_assert(lb2, "GE", l2, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lower2 bound too low",__LINE__,error)
          l2 = lb2
       ENDIF
       IF (PRESENT (ub2)) THEN
          CALL dbcsr_assert(ub2, "LE", u2, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "upper2 bound too high",__LINE__,error)
          u2 = ub2
       ENDIF
       IF (PRESENT (value)) THEN
          area%d%c2_sp(l:u, l2:u2) = value%c_sp
       ELSE
          area%d%c2_sp(l:u, l2:u2) = CMPLX(0.0, 0.0, real_4)
       ENDIF
    CASE (dbcsr_type_complex_8_2d)
       l = LBOUND(area%d%c2_dp, 1)
       u = UBOUND(area%d%c2_dp, 1)
       l2 = LBOUND(area%d%c2_dp, 2)
       u2 = UBOUND(area%d%c2_dp, 2)
       IF (PRESENT (lb)) THEN
          CALL dbcsr_assert(lb, "GE", l, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lower bound too low",__LINE__,error)
          l = lb
       ENDIF
       IF (PRESENT (ub)) THEN
          CALL dbcsr_assert(ub, "LE", u, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "upper bound too high",__LINE__,error)
          u = ub
       ENDIF
       IF (PRESENT (lb2)) THEN
          CALL dbcsr_assert(lb2, "GE", l2, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "lower2 bound too low",__LINE__,error)
          l2 = lb2
       ENDIF
       IF (PRESENT (ub2)) THEN
          CALL dbcsr_assert(ub2, "LE", u2, dbcsr_fatal_level,&
               dbcsr_caller_error, routineN, "upper2 bound too high",__LINE__,error)
          u2 = ub2
       ENDIF
       IF (PRESENT (value)) THEN
          area%d%c2_dp(l:u, l2:u2) = value%c_dp
       ELSE
          area%d%c2_dp(l:u, l2:u2) = CMPLX(0.0, 0.0, real_8)
       ENDIF
    CASE default
       CALL dbcsr_assert(.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid or unsupported data type.",__LINE__,error)
    END SELECT

    ! CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_data_clear0


! *****************************************************************************
!> \brief Copies a block subset
!> \param[in,out] dst  target data area
!> \param[in] dst_rs   target block row size (logical)
!> \param[in] dst_cs   target block column size (logical)
!> \param[in] dst_tr   whether target block is transposed
!> \param[in] src      source data area
!> \param[in] src_rs   source block row size (logical)
!> \param[in] src_cs   source block column size (logical)
!> \param[in] src_tr   whether source block is transposed
!> \param[in] dst_r_lb   first row in target
!> \param[in] dst_c_lb   first column in target
!> \param[in] src_r_lb   first_row in source
!> \param[in] src_r_lb   first column in target
!> \param[in] nrow       number of rows to copy
!> \param[in] ncol       number of columns to copy
!> \param[in] dst_offset offset in target
!> \param[in] src_offset offset in source
! *****************************************************************************
  SUBROUTINE block_partial_copy_a(dst, dst_rs, dst_cs, dst_tr,&
       src, src_rs, src_cs, src_tr,&
       dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol,&
       dst_offset, src_offset)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dst
    INTEGER, INTENT(IN)                      :: dst_rs, dst_cs
    LOGICAL                                  :: dst_tr
    TYPE(dbcsr_data_obj), INTENT(IN)         :: src
    INTEGER, INTENT(IN)                      :: src_rs, src_cs
    LOGICAL                                  :: src_tr
    INTEGER, INTENT(IN)                      :: dst_r_lb, dst_c_lb, src_r_lb, &
                                                src_c_lb, nrow, ncol
    INTEGER, INTENT(IN), OPTIONAL            :: dst_offset, src_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'block_partial_copy_a', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: verification = careful_mod

    INTEGER                                  :: dst_o, src_o
    LOGICAL                                  :: src_is_2d
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (careful_mod) &
         CALL dbcsr_assert (dbcsr_type_2d_to_1d (dbcsr_data_get_type (dst)),&
         "EQ", dbcsr_type_2d_to_1d (dbcsr_data_get_type(src)),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Incompatible data types.",__LINE__,error)
    dst_o = 0 ; src_o = 0
    IF (PRESENT (dst_offset)) dst_o = dst_offset
    IF (PRESENT (src_offset)) src_o = src_offset
    IF (verification) THEN
       CALL dbcsr_assert (dst_r_lb+nrow-1, "LE", dst_rs, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN, "Incompatible dst row sizes",__LINE__,error)
       CALL dbcsr_assert (dst_c_lb+ncol-1, "LE", dst_cs, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN, "Incompatible dst col sizes",__LINE__,error)
       CALL dbcsr_assert (src_r_lb+nrow-1, "LE", src_rs, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN, "Incompatible src row sizes",__LINE__,error)
       CALL dbcsr_assert (src_c_lb+ncol-1, "LE", src_cs, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN, "Incompatible src col sizes",__LINE__,error)
    ENDIF
    !
    src_is_2d = dbcsr_type_is_2d (dbcsr_data_get_type (src))
    SELECT CASE (dbcsr_data_get_type (dst))
    CASE (dbcsr_type_real_4)
       IF (src_is_2d) THEN
          CALL dbcsr_block_partial_copy(dst%d%r_sp, dst_rs, dst_cs, dst_tr,&
               src%d%r2_sp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol,&
               dst_offset=dst_o)
       ELSE
          CALL dbcsr_block_partial_copy(dst%d%r_sp, dst_rs, dst_cs, dst_tr,&
               src%d%r_sp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol,&
               dst_offset=dst_o, src_offset=src_o)
       ENDIF
    CASE (dbcsr_type_real_8)
       IF (src_is_2d) THEN
          CALL dbcsr_block_partial_copy(dst%d%r_dp, dst_rs, dst_cs, dst_tr,&
               src%d%r2_dp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol,&
               dst_offset=dst_o)
       ELSE
          CALL dbcsr_block_partial_copy(dst%d%r_dp, dst_rs, dst_cs, dst_tr,&
               src%d%r_dp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol,&
               dst_offset=dst_o, src_offset=src_o)
       ENDIF
    CASE (dbcsr_type_complex_4)
       IF (src_is_2d) THEN
          CALL dbcsr_block_partial_copy(dst%d%c_sp, dst_rs, dst_cs, dst_tr,&
               src%d%c2_sp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol,&
               dst_offset=dst_o)
       ELSE
          CALL dbcsr_block_partial_copy(dst%d%c_sp, dst_rs, dst_cs, dst_tr,&
               src%d%c_sp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol,&
               dst_offset=dst_o, src_offset=src_o)
       ENDIF
    CASE (dbcsr_type_complex_8)
       IF (src_is_2d) THEN
          CALL dbcsr_block_partial_copy(dst%d%c_dp, dst_rs, dst_cs, dst_tr,&
               src%d%c2_dp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol,&
               dst_offset=dst_o)
       ELSE
          CALL dbcsr_block_partial_copy(dst%d%c_dp, dst_rs, dst_cs, dst_tr,&
               src%d%c_dp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol,&
               dst_offset=dst_o, src_offset=src_o)
       ENDIF
    CASE (dbcsr_type_real_4_2d)
       IF (src_is_2d) THEN
          CALL dbcsr_block_partial_copy(dst%d%r2_sp, dst_rs, dst_cs, dst_tr,&
               src%d%r2_sp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol)
       ELSE
          CALL dbcsr_block_partial_copy(dst%d%r2_sp, dst_rs, dst_cs, dst_tr,&
               src%d%r_sp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol,&
               src_offset=src_o)
       ENDIF
    CASE (dbcsr_type_real_8_2d)
       IF (src_is_2d) THEN
          CALL dbcsr_block_partial_copy(dst%d%r2_dp, dst_rs, dst_cs, dst_tr,&
               src%d%r2_dp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol)
       ELSE
          CALL dbcsr_block_partial_copy(dst%d%r2_dp, dst_rs, dst_cs, dst_tr,&
               src%d%r_dp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol,&
               src_offset=src_o)
       ENDIF
    CASE (dbcsr_type_complex_4_2d)
       IF (src_is_2d) THEN
          CALL dbcsr_block_partial_copy(dst%d%c2_sp, dst_rs, dst_cs, dst_tr,&
               src%d%c2_sp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol)
       ELSE
          CALL dbcsr_block_partial_copy(dst%d%c2_sp, dst_rs, dst_cs, dst_tr,&
               src%d%c_sp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol,&
               src_offset=src_o)
       ENDIF
    CASE (dbcsr_type_complex_8_2d)
       IF (src_is_2d) THEN
          CALL dbcsr_block_partial_copy(dst%d%c2_dp, dst_rs, dst_cs, dst_tr,&
               src%d%c2_dp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol)
       ELSE
          CALL dbcsr_block_partial_copy(dst%d%c2_dp, dst_rs, dst_cs, dst_tr,&
               src%d%c_dp, src_rs, src_cs, src_tr,&
               dst_r_lb, dst_c_lb, src_r_lb, src_c_lb, nrow, ncol,&
               src_offset=src_o)
       ENDIF
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid data type.",__LINE__,error)
    END SELECT
  END SUBROUTINE block_partial_copy_a


! *****************************************************************************
!> \brief Sets the diagonal of a square data block
!> \par Off-diagonal values
!>      Other values are untouched.
!> \param[in,out] block_data  sets the diagonal of this data block
!> \param[in] diagonal        set diagonal of block_data to these values
!> \param[in] d               block size
!> \param[in] error           error
! *****************************************************************************
  SUBROUTINE set_block2d_diagonal_a (block_data, diagonal, d, &
       error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: block_data
    TYPE(dbcsr_data_obj), INTENT(IN)         :: diagonal
    INTEGER, INTENT(IN)                      :: d
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_block2d_diagonal_a', &
      routineP = moduleN//':'//routineN

    COMPLEX(kind=dp), DIMENSION(:, :), &
      POINTER                                :: c2_dp
    COMPLEX(kind=sp), DIMENSION(:, :), &
      POINTER                                :: c2_sp
    INTEGER                                  :: dt_block, dt_diag
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: r2_dp
    REAL(kind=sp), DIMENSION(:, :), POINTER  :: r2_sp

!   ---------------------------------------------------------------------------

    dt_block = dbcsr_data_get_type(block_data)
    dt_diag = dbcsr_data_get_type(diagonal)
    CALL dbcsr_assert (dbcsr_type_is_2d (dt_block), dbcsr_fatal_level,&
         dbcsr_internal_error, routineN, "Block must be 2-d", __LINE__,&
         error=error)
    CALL dbcsr_assert ("NOT", dbcsr_type_is_2d (dt_diag), dbcsr_fatal_level,&
         dbcsr_internal_error, routineN, "Diagonal can not be 2-d", __LINE__,&
         error=error)
    SELECT CASE (dt_diag)
    CASE (dbcsr_type_real_4)
       r2_sp => dbcsr_get_data_p_2d_s(block_data)
       CALL set_block2d_diagonal (r2_sp,&
            dbcsr_get_data_p_s(diagonal), d)
    CASE (dbcsr_type_real_8)
       r2_dp => dbcsr_get_data_p_2d_d(block_data)
       CALL set_block2d_diagonal (r2_dp,&
            dbcsr_get_data_p_d(diagonal), d)
    CASE (dbcsr_type_complex_4)
       c2_sp => dbcsr_get_data_p_2d_c(block_data)
       CALL set_block2d_diagonal (c2_sp,&
            dbcsr_get_data_p_c(diagonal), d)
    CASE (dbcsr_type_complex_8)
       c2_dp => dbcsr_get_data_p_2d_z(block_data)
       CALL set_block2d_diagonal (c2_dp,&
            dbcsr_get_data_p_z(diagonal), d)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid data type!",__LINE__,error)
    END SELECT
  END SUBROUTINE set_block2d_diagonal_a


! *****************************************************************************
!> \brief Gets the diagonal of a square data block
!> \param[in] block_data      gets the diagonal of this data block
!> \param[out] diagonal       values of the diagonal elements
!> \param[in] d               dimension of block
!> \param[in] error           error
! *****************************************************************************
  SUBROUTINE get_block2d_diagonal_a (block_data, diagonal, d, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: block_data, diagonal
    INTEGER, INTENT(IN)                      :: d
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_block2d_diagonal_a', &
      routineP = moduleN//':'//routineN

    COMPLEX(kind=dp), DIMENSION(:), POINTER  :: c_dp
    COMPLEX(kind=dp), DIMENSION(:, :), &
      POINTER                                :: c2_dp
    COMPLEX(kind=sp), DIMENSION(:), POINTER  :: c_sp
    COMPLEX(kind=sp), DIMENSION(:, :), &
      POINTER                                :: c2_sp
    INTEGER                                  :: dt_block, dt_diag
    REAL(kind=dp), DIMENSION(:), POINTER     :: r_dp
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: r2_dp
    REAL(kind=sp), DIMENSION(:), POINTER     :: r_sp
    REAL(kind=sp), DIMENSION(:, :), POINTER  :: r2_sp

!   ---------------------------------------------------------------------------

    dt_block = dbcsr_data_get_type(block_data)
    dt_diag = dbcsr_data_get_type(diagonal)
    CALL dbcsr_assert (dbcsr_type_is_2d (dt_block), dbcsr_fatal_level,&
         dbcsr_internal_error, routineN, "Block must be 2-d", __LINE__,&
         error=error)
    CALL dbcsr_assert ("NOT", dbcsr_type_is_2d (dt_diag), dbcsr_fatal_level,&
         dbcsr_internal_error, routineN, "Diagonal can not be 2-d", __LINE__,&
         error=error)
    SELECT CASE (dt_diag)
    CASE (dbcsr_type_real_4)
       r2_sp => dbcsr_get_data_p_2d_s(block_data)
       r_sp => dbcsr_get_data_p_s(diagonal)
       CALL get_block2d_diagonal_s (r2_sp, r_sp, d)
    CASE (dbcsr_type_real_8)
       r2_dp => dbcsr_get_data_p_2d_d(block_data)
       r_dp => dbcsr_get_data_p_d(diagonal)
       CALL get_block2d_diagonal_d (r2_dp, r_dp, d)
    CASE (dbcsr_type_complex_4)
       c2_sp => dbcsr_get_data_p_2d_c(block_data)
       c_sp => dbcsr_get_data_p_c(diagonal)
       CALL get_block2d_diagonal_c (c2_sp, c_sp, d)
    CASE (dbcsr_type_complex_8)
       c2_dp => dbcsr_get_data_p_2d_z(block_data)
       c_dp => dbcsr_get_data_p_z(diagonal)
       CALL get_block2d_diagonal_z (c2_dp, c_dp, d)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid data type!",__LINE__,error)
    END SELECT
  END SUBROUTINE get_block2d_diagonal_a

! *****************************************************************************
!> \brief Adds a value to the diagonal of the block
!> \param[in,out] block_a      Sets the diagonal of this block
!> \param[in] alpha_scalar     Adds this value to the diagonal elements of the
!>                             block
!> \param[in] block_dimension  Dimension of the block
! *****************************************************************************
  SUBROUTINE block_add_on_diag_anytype (block_a, alpha_scalar, block_dimension,&
       imin, imax, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: block_a
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: alpha_scalar
    INTEGER, INTENT(IN)                      :: block_dimension
    INTEGER, INTENT(in), OPTIONAL            :: imin, imax
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'block_add_on_diag_anytype', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (dbcsr_type_2d_to_1d(dbcsr_data_get_type(block_a)), "EQ",&
         alpha_scalar%data_type, dbcsr_fatal_level, dbcsr_wrong_args_error,&
         routineN, "Mismatched data types.", __LINE__, error=error)
    SELECT CASE (dbcsr_data_get_type(block_a))
    CASE (dbcsr_type_real_4)
       CALL block_add_on_diag (&
            block_dimension, block_a%d%r_sp, alpha_scalar%r_sp, imin, imax)
    CASE (dbcsr_type_real_8)
       CALL block_add_on_diag (&
            block_dimension, block_a%d%r_dp, alpha_scalar%r_dp, imin, imax)
    CASE (dbcsr_type_complex_4)
       CALL block_add_on_diag (&
            block_dimension, block_a%d%c_sp, alpha_scalar%c_sp, imin, imax)
    CASE (dbcsr_type_complex_8)
       CALL block_add_on_diag (&
            block_dimension, block_a%d%c_dp, alpha_scalar%c_dp, imin, imax)
    CASE (dbcsr_type_real_4_2d)
       CALL block_add_on_diag (&
            block_dimension, block_a%d%r2_sp, alpha_scalar%r_sp, imin, imax)
    CASE (dbcsr_type_real_8_2d)
       CALL block_add_on_diag (&
            block_dimension, block_a%d%r2_dp, alpha_scalar%r_dp, imin, imax)
    CASE (dbcsr_type_complex_4_2d)
       CALL block_add_on_diag (&
            block_dimension, block_a%d%c2_sp, alpha_scalar%c_sp, imin, imax)
    CASE (dbcsr_type_complex_8_2d)
       CALL block_add_on_diag (&
            block_dimension, block_a%d%c2_dp, alpha_scalar%c_dp, imin, imax)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid data type!",__LINE__,error)
    END SELECT
  END SUBROUTINE block_add_on_diag_anytype

! *****************************************************************************
!> \brief Adds two blocks
!> \param[in,out] block_a      Block to add to
!> \param[in] block_b          Block to add to block_a
!> \param[inout] error         error
! *****************************************************************************
  SUBROUTINE block_add_anytype (block_a, block_b, len, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: block_a
    TYPE(dbcsr_data_obj), INTENT(IN)         :: block_b
    INTEGER, INTENT(IN), OPTIONAL            :: len
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'block_add_anytype', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n

!   ---------------------------------------------------------------------------

    IF (careful_mod) &
         CALL dbcsr_assert (dbcsr_data_get_type(block_a), "EQ",&
         dbcsr_data_get_type(block_a),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Mismatched data types.", __LINE__, error=error)
    IF (PRESENT (len)) THEN
       n = len
       CALL dbcsr_assert (dbcsr_data_get_size (block_b), "GE", n,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Block B too small.", __LINE__, error=error)
    ELSE
       n = dbcsr_data_get_size_referenced (block_b)
    ENDIF
    CALL dbcsr_assert (dbcsr_data_get_size (block_a), "GE", n,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Block A too small.", __LINE__, error=error)
    SELECT CASE (dbcsr_data_get_type(block_a))
    CASE (dbcsr_type_real_4)
       CALL block_add_s (block_a%d%r_sp, block_b%d%r_sp, n)
    CASE (dbcsr_type_real_8)
       CALL block_add_d (block_a%d%r_dp, block_b%d%r_dp, n)
    CASE (dbcsr_type_complex_4)
       CALL block_add_c (block_a%d%c_sp, block_b%d%c_sp, n)
    CASE (dbcsr_type_complex_8)
       CALL block_add_z (block_a%d%c_dp, block_b%d%c_dp, n)
    CASE (dbcsr_type_real_4_2d)
       CALL block_add_s (block_a%d%r2_sp, block_b%d%r2_sp, n)
    CASE (dbcsr_type_real_8_2d)
       CALL block_add_d (block_a%d%r2_dp, block_b%d%r2_dp, n)
    CASE (dbcsr_type_complex_4_2d)
       CALL block_add_c (block_a%d%c2_sp, block_b%d%c2_sp, n)
    CASE (dbcsr_type_complex_8_2d)
       CALL block_add_z (block_a%d%c2_dp, block_b%d%c2_dp, n)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid data type!",__LINE__,error)
    END SELECT
  END SUBROUTINE block_add_anytype



! *****************************************************************************
!> \brief Calculates frobenious norm of a block
!>
!> \param[in] area       encapsulated data area
!> \param[in,out] error  error
!> \result norm     single-precision real frobenius norm
! *****************************************************************************
  SUBROUTINE block_norm_frob_anytype_s (area, norm, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    REAL(kind=sp), INTENT(OUT)               :: norm
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'block_norm_frob_anytype_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handler, error)

    SELECT CASE (area%d%data_type)
    CASE (dbcsr_type_real_8)
       CALL block_norm_frob_d_s (SIZE(area%d%r_dp), area%d%r_dp, norm)
    CASE (dbcsr_type_real_4)
       CALL block_norm_frob_s_s (SIZE(area%d%r_sp), area%d%r_sp, norm)
    CASE (dbcsr_type_complex_8)
       CALL block_norm_frob_z_s (SIZE(area%d%c_dp), area%d%c_dp, norm)
    CASE (dbcsr_type_complex_4)
       CALL block_norm_frob_c_s (SIZE(area%d%c_sp), area%d%c_sp, norm)
    CASE (dbcsr_type_real_8_2d)
       CALL block_norm_frob_d_s (SIZE(area%d%r2_dp), area%d%r2_dp, norm)
    CASE (dbcsr_type_real_4_2d)
       CALL block_norm_frob_s_s (SIZE(area%d%r2_sp), area%d%r2_sp, norm)
    CASE (dbcsr_type_complex_8_2d)
       CALL block_norm_frob_z_s (SIZE(area%d%c2_dp), area%d%c2_dp, norm)
    CASE (dbcsr_type_complex_4_2d)
       CALL block_norm_frob_c_s (SIZE(area%d%c2_sp), area%d%c2_sp, norm)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
    IF (careful_mod) CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE block_norm_frob_anytype_s

  SUBROUTINE block_norm_frob_s_s (m, blk, norm)
    INTEGER, INTENT(IN)                      :: m
    REAL(KIND=real_4), DIMENSION(m), &
      INTENT(IN)                             :: blk
    REAL(KIND=real_4), INTENT(OUT)           :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'block_norm_frob_s_s', &
      routineP = moduleN//':'//routineN

    norm = SQRT (REAL (SUM (blk**2), KIND=real_4))
  END SUBROUTINE block_norm_frob_s_s
  SUBROUTINE block_norm_frob_d_s (m, blk, norm)
    INTEGER, INTENT(IN)                      :: m
    REAL(KIND=real_8), DIMENSION(m), &
      INTENT(IN)                             :: blk
    REAL(KIND=real_4), INTENT(OUT)           :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'block_norm_frob_d_s', &
      routineP = moduleN//':'//routineN

    norm = SQRT (REAL (SUM (blk**2), KIND=real_4))
  END SUBROUTINE block_norm_frob_d_s
  SUBROUTINE block_norm_frob_c_s (m, blk, norm)
    INTEGER, INTENT(IN)                      :: m
    COMPLEX(KIND=real_4), DIMENSION(m), &
      INTENT(IN)                             :: blk
    REAL(KIND=real_4), INTENT(OUT)           :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'block_norm_frob_c_s', &
      routineP = moduleN//':'//routineN

    norm = SQRT (REAL (SUM (CONJG(blk) * blk), KIND=real_4))
  END SUBROUTINE block_norm_frob_c_s
  SUBROUTINE block_norm_frob_z_s (m, blk, norm)
    INTEGER, INTENT(IN)                      :: m
    COMPLEX(KIND=real_8), DIMENSION(m), &
      INTENT(IN)                             :: blk
    REAL(KIND=real_4), INTENT(OUT)           :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'block_norm_frob_z_s', &
      routineP = moduleN//':'//routineN

    norm = SQRT (REAL (SUM (CONJG(blk)* blk), KIND=real_4))
  END SUBROUTINE block_norm_frob_z_s



#include "dbcsr_block_operations_d.f90"
#include "dbcsr_block_operations_z.f90"
#include "dbcsr_block_operations_s.f90"
#include "dbcsr_block_operations_c.f90"

END MODULE dbcsr_block_operations

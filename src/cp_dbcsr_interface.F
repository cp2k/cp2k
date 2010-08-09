!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR interface in CP2K
!> \author  VE
!> \date    2010
!> \version 0.1
!>
!> <b>Modification history:</b>
!> - Created 2010
! *****************************************************************************
MODULE cp_dbcsr_interface
  USE array_types,                     ONLY: array_i1d_obj
  USE dbcsr_dist_operations
  USE dbcsr_error_handling,            ONLY: dbcsr_error_type
  USE dbcsr_io
  USE dbcsr_methods,                   ONLY: dbcsr_get_matrix_type,&
                                             dbcsr_init,&
                                             dbcsr_init_p,&
                                             dbcsr_nfullcols_total,&
                                             dbcsr_nfullrows_total,&
                                             dbcsr_release,&
                                             dbcsr_release_p
  USE dbcsr_operations,                ONLY: &
       dbcsr_add, dbcsr_add_on_diag, dbcsr_btriu, dbcsr_copy, dbcsr_filter, &
       dbcsr_get_block_diag, dbcsr_get_diag, dbcsr_hadamard_product, &
       dbcsr_multiply, dbcsr_norm, dbcsr_scale, dbcsr_scale_by_vector, &
       dbcsr_set, dbcsr_set_diag, dbcsr_sum_replicated, dbcsr_trace
  USE dbcsr_transformations,           ONLY: dbcsr_distribute,&
                                             dbcsr_new_transposed,&
                                             dbcsr_replicate_all
  USE dbcsr_types,                     ONLY: dbcsr_2d_array_type,&
                                             dbcsr_distribution_obj,&
                                             dbcsr_obj
  USE dbcsr_util,                      ONLY: dbcsr_checksum,&
                                             dbcsr_verify_matrix
  USE dbcsr_work_operations,           ONLY: dbcsr_create,&
                                             dbcsr_finalize,&
                                             dbcsr_work_create
  USE kinds,                           ONLY: dp,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop

  !$ USE OMP_LIB
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_interface'

  !
  ! Interface to libdbcsr (contains cp2k timing and error)
  PUBLIC :: cp_dbcsr_set
  PUBLIC :: cp_dbcsr_add
  PUBLIC :: cp_dbcsr_scale
  PUBLIC :: cp_dbcsr_scale_by_vector
  PUBLIC :: cp_dbcsr_hadamard_product
  PUBLIC :: cp_dbcsr_transposed
  PUBLIC :: cp_dbcsr_multiply
  PUBLIC :: cp_dbcsr_copy
  PUBLIC :: cp_dbcsr_add_on_diag
  PUBLIC :: cp_dbcsr_get_block_diag
  PUBLIC :: cp_dbcsr_set_diag
  PUBLIC :: cp_dbcsr_get_diag
  PUBLIC :: cp_dbcsr_filter
  PUBLIC :: cp_dbcsr_finalize
  PUBLIC :: cp_dbcsr_create
  PUBLIC :: cp_dbcsr_work_create
  PUBLIC :: cp_dbcsr_verify_matrix
  PUBLIC :: cp_dbcsr_btriu
  PUBLIC :: cp_dbcsr_sum_replicated
  PUBLIC :: cp_dbcsr_checksum
  PUBLIC :: cp_dbcsr_trace
  PUBLIC :: cp_dbcsr_print
  PUBLIC :: cp_dbcsr_init
  PUBLIC :: cp_dbcsr_init_p
  PUBLIC :: cp_dbcsr_release
  PUBLIC :: cp_dbcsr_release_p
  PUBLIC :: cp_dbcsr_distribute
  PUBLIC :: cp_dbcsr_replicate_all
  PUBLIC :: cp_dbcsr_norm

  INTERFACE cp_dbcsr_trace
     MODULE PROCEDURE cp_dbcsr_trace_ab_d
     MODULE PROCEDURE cp_dbcsr_trace_ab_s
     MODULE PROCEDURE cp_dbcsr_trace_a_d
     MODULE PROCEDURE cp_dbcsr_trace_a_s
  END INTERFACE

  INTERFACE cp_dbcsr_set
     MODULE PROCEDURE cp_dbcsr_set_d
     MODULE PROCEDURE cp_dbcsr_set_s
     MODULE PROCEDURE cp_dbcsr_set_c
     MODULE PROCEDURE cp_dbcsr_set_z
  END INTERFACE

  INTERFACE cp_dbcsr_add
     MODULE PROCEDURE cp_dbcsr_add_d
     MODULE PROCEDURE cp_dbcsr_add_s
     MODULE PROCEDURE cp_dbcsr_add_c
     MODULE PROCEDURE cp_dbcsr_add_z
  END INTERFACE

  INTERFACE cp_dbcsr_scale
     MODULE PROCEDURE cp_dbcsr_scale_d
     MODULE PROCEDURE cp_dbcsr_scale_s
     MODULE PROCEDURE cp_dbcsr_scale_c
     MODULE PROCEDURE cp_dbcsr_scale_z
  END INTERFACE

  INTERFACE cp_dbcsr_scale_by_vector
     MODULE PROCEDURE cp_dbcsr_scale_by_vector_d
     MODULE PROCEDURE cp_dbcsr_scale_by_vector_s
     MODULE PROCEDURE cp_dbcsr_scale_by_vector_c
     MODULE PROCEDURE cp_dbcsr_scale_by_vector_z
  END INTERFACE

  INTERFACE cp_dbcsr_multiply
     MODULE PROCEDURE cp_dbcsr_multiply_d
     MODULE PROCEDURE cp_dbcsr_multiply_s
     MODULE PROCEDURE cp_dbcsr_multiply_c
     MODULE PROCEDURE cp_dbcsr_multiply_z
  END INTERFACE

  PRIVATE

CONTAINS

! *****************************************************************************
!> \brief cp2k interface to libdbcsr.
!>     The cp2k error is passed as an argument.
!>     Those procedure can de timed with the cp2k timer.    
!> \author vw
! *****************************************************************************
  SUBROUTINE cp_dbcsr_norm(matrix, which_norm, norm_scalar, norm_vector, error)

    TYPE(dbcsr_obj), INTENT(INOUT), TARGET   :: matrix
    INTEGER, INTENT(IN)                      :: which_norm
    REAL(dp), INTENT(OUT), OPTIONAL          :: norm_scalar
    REAL(dp), DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: norm_vector
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_norm', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    IF (PRESENT (norm_scalar)) THEN
       CALL dbcsr_norm(matrix, which_norm, norm_scalar=norm_scalar,&
            error=dbcsr_error)
    ELSEIF (PRESENT (norm_vector)) THEN
       CALL dbcsr_norm (matrix, which_norm, norm_vector=norm_vector,&
            error=dbcsr_error)
    ELSE
       CALL cp_assert (.FALSE., cp_wrong_args_error, cp_internal_error,&
            routineN, "Must pass either scalar or vector norm.", error=error)
    ENDIF

  END SUBROUTINE cp_dbcsr_norm

  SUBROUTINE cp_dbcsr_replicate_all(matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_replicate_all', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_replicate_all(matrix, dbcsr_error)
  END SUBROUTINE cp_dbcsr_replicate_all

  SUBROUTINE cp_dbcsr_distribute(matrix, fast, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    LOGICAL, INTENT(in), OPTIONAL            :: fast
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_distribute', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_distribute(matrix, fast, dbcsr_error)
  END SUBROUTINE cp_dbcsr_distribute

  SUBROUTINE cp_dbcsr_release_p (matrix, error)
    TYPE(dbcsr_obj), POINTER                 :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_release_p', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_release_p (matrix,dbcsr_error)

  END SUBROUTINE cp_dbcsr_release_p

  SUBROUTINE cp_dbcsr_release (matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_release', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_release (matrix)

  END SUBROUTINE cp_dbcsr_release

  SUBROUTINE cp_dbcsr_init (matrix, error)
    TYPE(dbcsr_obj), INTENT(OUT)             :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_init', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_init(matrix)

  END SUBROUTINE cp_dbcsr_init

  SUBROUTINE cp_dbcsr_init_p (matrix, error)
    TYPE(dbcsr_obj), POINTER                 :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_init_p', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_init_p(matrix,dbcsr_error)

  END SUBROUTINE cp_dbcsr_init_p

  SUBROUTINE cp_dbcsr_print(matrix, nodata, matlab_format, variable_name, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: nodata, matlab_format
    CHARACTER(*), INTENT(in), OPTIONAL       :: variable_name
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_print', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_print(matrix, nodata, matlab_format, variable_name, dbcsr_error)
  END SUBROUTINE cp_dbcsr_print

  SUBROUTINE cp_dbcsr_trace_a_d(matrix_a, trace, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(kind=real_8), INTENT(INOUT)         :: trace
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_trace_a_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_trace(matrix_a, trace, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_trace_a_d

  SUBROUTINE cp_dbcsr_trace_a_s(matrix_a, trace, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(kind=real_4), INTENT(INOUT)         :: trace
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_trace_a_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_trace(matrix_a, trace, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_trace_a_s

  SUBROUTINE cp_dbcsr_trace_ab_d(matrix_a, matrix_b, trace, trans_a, trans_b, local_sum, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a, matrix_b
    REAL(kind=real_8), INTENT(INOUT)         :: trace
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: trans_a, trans_b
    LOGICAL, INTENT(IN), OPTIONAL            :: local_sum
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_trace_ab_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_trace(matrix_a, matrix_b, trace, trans_a, trans_b, local_sum, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_trace_ab_d

  SUBROUTINE cp_dbcsr_trace_ab_s(matrix_a, matrix_b, trace, trans_a, trans_b, local_sum, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a, matrix_b
    REAL(kind=real_4), INTENT(INOUT)         :: trace
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: trans_a, trans_b
    LOGICAL, INTENT(IN), OPTIONAL            :: local_sum
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_trace_ab_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_trace(matrix_a, matrix_b, trace, trans_a, trans_b, local_sum, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_trace_ab_s

  FUNCTION cp_dbcsr_checksum(matrix, local, error) RESULT(checksum)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: local
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    REAL(KIND=dp)                            :: checksum

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_checksum', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    checksum = dbcsr_checksum(matrix, local, error=dbcsr_error)
  END FUNCTION cp_dbcsr_checksum

  SUBROUTINE cp_dbcsr_sum_replicated (matrix, error)
    TYPE(dbcsr_obj), INTENT(inout)           :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_sum_replicated', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_sum_replicated (matrix, dbcsr_error)
  END SUBROUTINE cp_dbcsr_sum_replicated

  SUBROUTINE cp_dbcsr_btriu(matrix_b, matrix_a, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_btriu', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_btriu(matrix_b, matrix_a, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_btriu

  SUBROUTINE cp_dbcsr_verify_matrix(m, error, verbosity, local)
    TYPE(dbcsr_obj), INTENT(IN)              :: m
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER, INTENT(IN), OPTIONAL            :: verbosity
    LOGICAL, INTENT(IN), OPTIONAL            :: local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_verify_matrix', &
      r = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_verify_matrix(m, verbosity, local, dbcsr_error)
  END SUBROUTINE cp_dbcsr_verify_matrix

  SUBROUTINE cp_dbcsr_finalize(matrix, resort, reshuffle, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: resort, reshuffle
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_finalize(matrix, resort, reshuffle, error=dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_finalize

  SUBROUTINE cp_dbcsr_work_create(matrix, nblks_guess, sizedata_guess, n,&
       error, work_mutable)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN), OPTIONAL            :: nblks_guess, sizedata_guess, n
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(in), OPTIONAL            :: work_mutable

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_work_create', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_work_create(matrix, nblks_guess, sizedata_guess, n,&
         work_mutable, dbcsr_error)
  END SUBROUTINE cp_dbcsr_work_create

  SUBROUTINE cp_dbcsr_create(matrix, name, dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, data_type, special, reuse,&
       mutable_work, replication_type, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    CHARACTER(len=*), INTENT(IN)             :: name
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    CHARACTER, INTENT(IN)                    :: matrix_type
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_size, col_blk_size
    INTEGER, INTENT(IN), OPTIONAL            :: nblks, nze
    CHARACTER, INTENT(IN), OPTIONAL          :: data_type
    LOGICAL, INTENT(IN), OPTIONAL            :: special, reuse, mutable_work
    CHARACTER, INTENT(IN), OPTIONAL          :: replication_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_create', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_create(matrix, name, dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, data_type, special, reuse,&
       mutable_work, replication_type, dbcsr_error)
  END SUBROUTINE cp_dbcsr_create

  SUBROUTINE cp_dbcsr_filter(matrix, eps, method, use_absolute, filter_diag, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(dp), INTENT(IN)                     :: eps
    INTEGER, INTENT(IN), OPTIONAL            :: method
    LOGICAL, INTENT(in), OPTIONAL            :: use_absolute, filter_diag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_filter', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_filter(matrix, eps, method, use_absolute, filter_diag, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_filter

  SUBROUTINE cp_dbcsr_set_diag(matrix, diag, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(dp), DIMENSION(:), INTENT(IN)       :: diag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_set_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_set_diag(matrix, diag, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_set_diag

  SUBROUTINE cp_dbcsr_get_diag(matrix, diag, error)

    TYPE(dbcsr_obj), INTENT(IN), TARGET      :: matrix
    REAL(dp), DIMENSION(:), INTENT(INOUT), &
      TARGET                                 :: diag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_get_diag(matrix, diag, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_get_diag

  SUBROUTINE cp_dbcsr_get_block_diag(matrix, diag, error)

    TYPE(dbcsr_obj), INTENT(IN), TARGET      :: matrix
    TYPE(dbcsr_obj), INTENT(INOUT), TARGET   :: diag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_block_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_get_block_diag(matrix, diag, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_get_block_diag

  SUBROUTINE cp_dbcsr_add_on_diag(matrix, alpha_scalar, first_row, last_row, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(KIND=dp), INTENT(IN)                :: alpha_scalar
    INTEGER, INTENT(in), OPTIONAL            :: first_row, last_row
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_add_on_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_add_on_diag(matrix, alpha_scalar, first_row, last_row, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_add_on_diag

  SUBROUTINE cp_dbcsr_copy(matrix_b, matrix_a, name, error, keep_sparsity,&
       shallow_data, keep_imaginary)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity, shallow_data, &
                                                keep_imaginary

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_copy(matrix_b, matrix_a, name, keep_sparsity,&
         shallow_data, keep_imaginary, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_copy

  SUBROUTINE cp_dbcsr_multiply_s(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, left_set, right_set, error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    REAL(real_4), INTENT(IN)                 :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    REAL(real_4), INTENT(IN)                 :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_2d_array_type), OPTIONAL, &
      POINTER                                :: left_set, right_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(int_8), INTENT(OUT), OPTIONAL    :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_s', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: shape_a, shape_b, trans_a, &
                                                trans_b
    INTEGER                                  :: timing_handle, &
                                                timing_handle_detail
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    trans_a = transa
    trans_b = transb
    CALL uppercase(trans_a)
    CALL uppercase(trans_b)
    shape_a='R'
    IF(dbcsr_nfullcols_total(matrix_a).EQ.dbcsr_nfullrows_total(matrix_a)) shape_a='S'
    shape_b='R'
    IF(dbcsr_nfullcols_total(matrix_b).EQ.dbcsr_nfullrows_total(matrix_b)) shape_b='S'
    CALL timeset('dbcsr_mult_'//trans_a//shape_a//dbcsr_get_matrix_type(matrix_a)//'_'&
         //trans_b//shape_b//dbcsr_get_matrix_type(matrix_b), timing_handle_detail)
    CALL dbcsr_multiply(transa, transb,&
         alpha, matrix_a, matrix_b, beta, matrix_c,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, left_set=left_set, right_set=right_set, error=dbcsr_error, flop=flop)
    CALL timestop(timing_handle_detail)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_multiply_s

  SUBROUTINE cp_dbcsr_multiply_d(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, left_set, right_set, error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    REAL(real_8), INTENT(IN)                 :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    REAL(real_8), INTENT(IN)                 :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_2d_array_type), OPTIONAL, &
      POINTER                                :: left_set, right_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(int_8), INTENT(OUT), OPTIONAL    :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_d', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: shape_a, shape_b, trans_a, &
                                                trans_b
    INTEGER                                  :: timing_handle, &
                                                timing_handle_detail
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)

    trans_a = transa
    trans_b = transb
    CALL uppercase(trans_a)
    CALL uppercase(trans_b)
    shape_a='R'
    IF(dbcsr_nfullcols_total(matrix_a).EQ.dbcsr_nfullrows_total(matrix_a)) shape_a='S'
    shape_b='R'
    IF(dbcsr_nfullcols_total(matrix_b).EQ.dbcsr_nfullrows_total(matrix_b)) shape_b='S'
    CALL timeset('dbcsr_mult_'//trans_a//shape_a//dbcsr_get_matrix_type(matrix_a)//'_'&
         //trans_b//shape_b//dbcsr_get_matrix_type(matrix_b), timing_handle_detail)
    CALL dbcsr_multiply(transa, transb,&
         alpha, matrix_a, matrix_b, beta, matrix_c,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, left_set, right_set, dbcsr_error, flop)
    CALL timestop(timing_handle_detail)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_multiply_d

  SUBROUTINE cp_dbcsr_multiply_c(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, left_set, right_set, error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    COMPLEX(real_4), INTENT(IN)              :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    COMPLEX(real_4), INTENT(IN)              :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_2d_array_type), OPTIONAL, &
      POINTER                                :: left_set, right_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(int_8), INTENT(OUT), OPTIONAL    :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_c', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: shape_a, shape_b, trans_a, &
                                                trans_b
    INTEGER                                  :: timing_handle, &
                                                timing_handle_detail
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)

    trans_a = transa
    trans_b = transb
    CALL uppercase(trans_a)
    CALL uppercase(trans_b)
    shape_a='R'
    IF(dbcsr_nfullcols_total(matrix_a).EQ.dbcsr_nfullrows_total(matrix_a)) shape_a='S'
    shape_b='R'
    IF(dbcsr_nfullcols_total(matrix_b).EQ.dbcsr_nfullrows_total(matrix_b)) shape_b='S'
    CALL timeset('dbcsr_mult_'//trans_a//shape_a//dbcsr_get_matrix_type(matrix_a)//'_'&
         //trans_b//shape_b//dbcsr_get_matrix_type(matrix_b), timing_handle_detail)
    CALL dbcsr_multiply(transa, transb,&
         alpha, matrix_a, matrix_b, beta, matrix_c,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, left_set=left_set, right_set=right_set, error=dbcsr_error, flop=flop)
    CALL timestop(timing_handle_detail)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_multiply_c

  SUBROUTINE cp_dbcsr_multiply_z(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, left_set, right_set, error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    COMPLEX(real_8), INTENT(IN)              :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    COMPLEX(real_8), INTENT(IN)              :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_2d_array_type), OPTIONAL, &
      POINTER                                :: left_set, right_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(int_8), INTENT(OUT), OPTIONAL    :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_z', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: shape_a, shape_b, trans_a, &
                                                trans_b
    INTEGER                                  :: timing_handle, &
                                                timing_handle_detail
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)

    trans_a = transa
    trans_b = transb
    CALL uppercase(trans_a)
    CALL uppercase(trans_b)
    shape_a='R'
    IF(dbcsr_nfullcols_total(matrix_a).EQ.dbcsr_nfullrows_total(matrix_a)) shape_a='S'
    shape_b='R'
    IF(dbcsr_nfullcols_total(matrix_b).EQ.dbcsr_nfullrows_total(matrix_b)) shape_b='S'
    CALL timeset('dbcsr_mult_'//trans_a//shape_a//dbcsr_get_matrix_type(matrix_a)//'_'&
         //trans_b//shape_b//dbcsr_get_matrix_type(matrix_b), timing_handle_detail)
    CALL dbcsr_multiply(transa, transb,&
         alpha, matrix_a, matrix_b, beta, matrix_c,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, left_set=left_set, right_set=right_set, error=dbcsr_error, flop=flop)
    CALL timestop(timing_handle_detail)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_multiply_z

  SUBROUTINE cp_dbcsr_transposed (transposed, normal, shallow_data_copy,&
       transpose_data, transpose_distribution, use_distribution, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: transposed
    TYPE(dbcsr_obj), INTENT(IN)              :: normal
    LOGICAL, INTENT(IN), OPTIONAL            :: shallow_data_copy, &
                                                transpose_data, &
                                                transpose_distribution
    TYPE(dbcsr_distribution_obj), &
      INTENT(IN), OPTIONAL                   :: use_distribution
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_transposed', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_new_transposed(transposed, normal, shallow_data_copy,&
         transpose_data, transpose_distribution,&
         use_distribution=use_distribution, error=dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_transposed

  SUBROUTINE cp_dbcsr_hadamard_product(matrix_a, matrix_b, matrix_c, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_hadamard_product', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_hadamard_product(matrix_a, matrix_b, matrix_c, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_hadamard_product


  SUBROUTINE cp_dbcsr_scale_by_vector_d(matrix_a, alpha, side, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(real_8), DIMENSION(:), INTENT(IN), &
      TARGET                                 :: alpha
    CHARACTER(LEN=*), INTENT(IN)             :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_by_vector_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale_by_vector(matrix_a, alpha, side, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_by_vector_d

  SUBROUTINE cp_dbcsr_scale_by_vector_s(matrix_a, alpha, side, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(real_4), DIMENSION(:), INTENT(IN), &
      TARGET                                 :: alpha
    CHARACTER(LEN=*), INTENT(IN)             :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_by_vector_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale_by_vector(matrix_a, alpha, side, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_by_vector_s

  SUBROUTINE cp_dbcsr_scale_by_vector_z(matrix_a, alpha, side, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    COMPLEX(real_8), DIMENSION(:), &
      INTENT(IN), TARGET                     :: alpha
    CHARACTER(LEN=*), INTENT(IN)             :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_by_vector_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale_by_vector(matrix_a, alpha, side, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_by_vector_z

  SUBROUTINE cp_dbcsr_scale_by_vector_c(matrix_a, alpha, side, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    COMPLEX(real_4), DIMENSION(:), &
      INTENT(IN), TARGET                     :: alpha
    CHARACTER(LEN=*), INTENT(IN)             :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_by_vector_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale_by_vector(matrix_a, alpha, side, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_by_vector_c

  SUBROUTINE cp_dbcsr_scale_d(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(real_8), INTENT(IN)                 :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale(matrix_a, alpha_scalar, last_column, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_d

  SUBROUTINE cp_dbcsr_scale_s(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(real_4), INTENT(IN)                 :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale(matrix_a, alpha_scalar, last_column, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_s

  SUBROUTINE cp_dbcsr_scale_z(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    COMPLEX(real_8), INTENT(IN)              :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale(matrix_a, alpha_scalar, last_column, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_z

  SUBROUTINE cp_dbcsr_scale_c(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    COMPLEX(real_4), INTENT(IN)              :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale(matrix_a, alpha_scalar, last_column, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_c

  SUBROUTINE cp_dbcsr_set_d(matrix, alpha, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(real_8), INTENT(IN)                 :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_set_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_set(matrix, alpha, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_set_d

  SUBROUTINE cp_dbcsr_set_s(matrix, alpha, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(real_4), INTENT(IN)                 :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_set_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_set(matrix, alpha, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_set_s

  SUBROUTINE cp_dbcsr_set_z(matrix, alpha, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    COMPLEX(real_8), INTENT(IN)              :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_set_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_set(matrix, alpha, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_set_z

  SUBROUTINE cp_dbcsr_set_c(matrix, alpha, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    COMPLEX(real_4), INTENT(IN)              :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_set_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_set(matrix, alpha, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_set_c

  SUBROUTINE cp_dbcsr_add_d(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    REAL(real_8), INTENT(IN)                 :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_add_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_add(matrix_a, matrix_b, alpha_scalar, beta_scalar, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_add_d

  SUBROUTINE cp_dbcsr_add_s(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    REAL(real_4), INTENT(IN)                 :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_add_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_add(matrix_a, matrix_b, alpha_scalar, beta_scalar, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_add_s

  SUBROUTINE cp_dbcsr_add_z(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    COMPLEX(real_8), INTENT(IN)              :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_add_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_add(matrix_a, matrix_b, alpha_scalar, beta_scalar, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_add_z

  SUBROUTINE cp_dbcsr_add_c(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    COMPLEX(real_4), INTENT(IN)              :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_add_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_add(matrix_a, matrix_b, alpha_scalar, beta_scalar, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_add_c

END MODULE cp_dbcsr_interface

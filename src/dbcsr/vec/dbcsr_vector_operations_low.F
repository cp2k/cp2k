!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief operations for skinny matrices/vectors expressed in dbcsr form
!> \par History
!>       2014.10 created [Florian Schiffmann]
!> \author Florian Schiffmann
! *****************************************************************************

MODULE dbcsr_vector_operations_low
  USE array_types,                     ONLY: array_data
  USE dbcsr_data_methods,              ONLY: dbcsr_get_data_p
  USE dbcsr_dist_methods,              ONLY: dbcsr_distribution_col_dist, &
                                             dbcsr_distribution_row_dist
  USE dbcsr_error_handling,            ONLY: dbcsr_assert, &
                                             dbcsr_caller_error, &
                                             dbcsr_error_type, &
                                             dbcsr_fatal_level, &
                                             dbcsr_error_set, &
                                             dbcsr_error_stop
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left, &
                                             dbcsr_iterator_next_block, &
                                             dbcsr_iterator_start, &
                                             dbcsr_iterator_stop
  USE dbcsr_methods,                   ONLY: dbcsr_distribution, &
                                             dbcsr_get_matrix_type
  USE dbcsr_operations,                ONLY: dbcsr_get_info, &
                                             dbcsr_set
  USE dbcsr_types,                     ONLY: &
       dbcsr_distribution_obj, dbcsr_iterator, dbcsr_obj, dbcsr_scalar_type, &
       dbcsr_type_antisymmetric, dbcsr_type_complex_4, dbcsr_type_complex_8, &
       dbcsr_type_no_symmetry, dbcsr_type_real_4, dbcsr_type_real_8, &
       dbcsr_type_symmetric
  USE kinds,                           ONLY: real_4, &
                                             real_8
  USE message_passing,                 ONLY: mp_bcast, &
                                             mp_sum
! types for mapping matrix blocks

  !$ USE OMP_LIB

  IMPLICIT NONE

  TYPE block_ptr_d
    REAL(real_8), DIMENSION(:, :), POINTER          :: ptr =>NULL()
    INTEGER                                         :: assigned_thread
  END TYPE
  TYPE block_ptr_s
    REAL(real_4), DIMENSION(:, :), POINTER          :: ptr =>NULL()
    INTEGER                                         :: assigned_thread
  END TYPE
  TYPE block_ptr_c
    COMPLEX(real_4), DIMENSION(:, :), POINTER       :: ptr =>NULL()
    INTEGER                                         :: assigned_thread
  END TYPE
  TYPE block_ptr_z
    COMPLEX(real_8), DIMENSION(:, :), POINTER       :: ptr =>NULL()
    INTEGER                                         :: assigned_thread
  END TYPE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_vector_operations_low'

  PUBLIC :: dbcsr_matrix_colvec_multiply_low

INTERFACE  dbcsr_matrix_colvec_multiply_low
   MODULE PROCEDURE dbcsr_matrix_colvec_multiply_low_any
   MODULE PROCEDURE dbcsr_matrix_colvec_multiply_low_d, dbcsr_matrix_colvec_multiply_low_s, &
                    dbcsr_matrix_colvec_multiply_low_c, dbcsr_matrix_colvec_multiply_low_z
END INTERFACE


CONTAINS

! finally we are here so we can make a final decision which argument types to pass, ... long way to here

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param vec_in ...
!> \param vec_out ...
!> \param alpha ...
!> \param beta ...
!> \param work_row ...
!> \param work_col ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_matrix_colvec_multiply_low_any(matrix, vec_in, vec_out, alpha, beta, work_row, work_col, error)
    TYPE(dbcsr_obj)                          :: matrix, vec_in, vec_out
    TYPE(dbcsr_scalar_type)                  :: alpha, beta
    TYPE(dbcsr_obj)                          :: work_row, work_col
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'dbcsr_matrix_colvec_multiply_low_any', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: matrix_type, handle

    CALL dbcsr_error_set(routineN, handle, error)

    matrix_type=matrix%m%data_type
    SELECT CASE(matrix_type)
    CASE(dbcsr_type_real_8)
       CALL dbcsr_matrix_colvec_multiply_low(matrix, vec_in, vec_out, alpha%r_dp, beta%r_dp, work_row, work_col, error)
    CASE(dbcsr_type_real_4)
       CALL dbcsr_matrix_colvec_multiply_low(matrix, vec_in, vec_out, alpha%r_sp, beta%r_sp, work_row, work_col, error)
    CASE(dbcsr_type_complex_8)
       CALL dbcsr_matrix_colvec_multiply_low(matrix, vec_in, vec_out, alpha%c_dp, beta%c_dp, work_row, work_col, error)
    CASE(dbcsr_type_complex_4)
       CALL dbcsr_matrix_colvec_multiply_low(matrix, vec_in, vec_out, alpha%c_sp, beta%c_sp, work_row, work_col, error)
    END SELECT

    CALL dbcsr_error_stop(handle, error)

  END SUBROUTINE dbcsr_matrix_colvec_multiply_low_any

#include "dbcsr_vector_operations_low_d.f90"
#include "dbcsr_vector_operations_low_z.f90"
#include "dbcsr_vector_operations_low_s.f90"
#include "dbcsr_vector_operations_low_c.f90"

END MODULE dbcsr_vector_operations_low

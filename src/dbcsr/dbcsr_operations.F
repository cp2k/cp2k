!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Higher-level operations on DBCSR matrices.
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.9
!>
!> <b>Modification history:</b>
!  - Created 2009-05-12
! *****************************************************************************

MODULE dbcsr_operations
  USE array_types,                     ONLY: array_i1d_obj,&
                                             array_release
  USE dbcsr_config,                    ONLY: is_initialized
  USE dbcsr_data_methods,              ONLY: dbcsr_scalar
  USE dbcsr_dist_operations,           ONLY: create_bl_distribution
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type,&
                                             dbcsr_failure_level,&
                                             dbcsr_fatal_level,&
                                             dbcsr_internal_error
  USE dbcsr_methods,                   ONLY: &
       dbcsr_distribution, dbcsr_distribution_mp, dbcsr_distribution_new, &
       dbcsr_distribution_release, dbcsr_distribution_row_dist, &
       dbcsr_get_data_type, dbcsr_get_info, dbcsr_init, dbcsr_mp_npcols, &
       dbcsr_release
  USE dbcsr_mm_cannon,                 ONLY: dbcsr_mm_cannon_clear_mempools,&
                                             dbcsr_mm_cannon_lib_finalize,&
                                             dbcsr_mm_cannon_lib_init,&
                                             dbcsr_mm_cannon_multiply
  USE dbcsr_operations_low,            ONLY: &
       dbcsr_add, dbcsr_add_on_diag, dbcsr_btriu, dbcsr_copy, &
       dbcsr_copy_into_existing, dbcsr_filter, dbcsr_frobenius_norm, &
       dbcsr_gershgorin_norm, dbcsr_get_block_diag, dbcsr_get_diag, &
       dbcsr_hadamard_product, dbcsr_init_random, dbcsr_maxabs, dbcsr_norm, &
       dbcsr_scale, dbcsr_scale_by_vector, dbcsr_scale_mat, dbcsr_set, &
       dbcsr_set_diag, dbcsr_sum_replicated, dbcsr_trace, dbcsr_triu
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_obj,&
                                             dbcsr_type_no_symmetry,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE dbcsr_work_operations,           ONLY: dbcsr_create,&
                                             dbcsr_finalize
  USE kinds,                           ONLY: dp,&
                                             int_1_size,&
                                             int_2_size,&
                                             int_4_size,&
                                             int_8,&
                                             int_8_size,&
                                             real_4,&
                                             real_8

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_operations'

  PUBLIC :: dbcsr_init_lib, dbcsr_finalize_lib,&
            dbcsr_multiply,&
            dbcsr_trace, dbcsr_add_on_diag,&
            dbcsr_set, dbcsr_scale, dbcsr_scale_mat, dbcsr_add, dbcsr_copy,&
            dbcsr_copy_into_existing,&
            dbcsr_get_diag, dbcsr_set_diag,&
            dbcsr_get_block_diag, dbcsr_hadamard_product,&
            dbcsr_filter, dbcsr_scale_by_vector,&
            dbcsr_btriu,&
            dbcsr_triu,&
            dbcsr_init_random,&
            dbcsr_sum_replicated, dbcsr_norm,&
            dbcsr_gershgorin_norm, dbcsr_maxabs, dbcsr_frobenius_norm,&
            dbcsr_clear_mempools

  INTERFACE dbcsr_multiply
     MODULE PROCEDURE dbcsr_mm_cannon_multiply
     MODULE PROCEDURE dbcsr_multiply_s, dbcsr_multiply_d,&
                      dbcsr_multiply_c, dbcsr_multiply_z
  END INTERFACE





 CONTAINS


! *****************************************************************************
!> \brief Initialize the DBCSR library
!>
!> Prepares the DBCSR library for use.
!> \param group ...
!> \param[in,out] error     error
! *****************************************************************************
  SUBROUTINE dbcsr_init_lib (group, error)
    INTEGER, INTENT(IN)                      :: group
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_lib', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!n_stack_buffers,mem_type, n_stack_mem_regions, stack_size, nstacks,&
!INTEGER, DIMENSION(3)                    :: nxstacks
!   ---------------------------------------------------------------------------
!TODO: problem: init/finalize are called by cp2k_runs AND f77_interface

    IF (is_initialized) RETURN
    CALL dbcsr_error_set(routineN, error_handle, error)
    !
    CALL dbcsr_assert (int_1_size, "EQ", 1,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Incorrect assumption of an 8-bit integer size!",&
         __LINE__, error=error)
    CALL dbcsr_assert (int_2_size, "EQ", 2,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Incorrect assumption of a 16-bit integer size!",&
         __LINE__, error=error)
    CALL dbcsr_assert (int_4_size, "EQ", 4,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Incorrect assumption of a 32-bit integer size!",&
         __LINE__, error=error)
    CALL dbcsr_assert (int_8_size, "EQ", 8,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Incorrect assumption of a 64-bit integer size!",&
         __LINE__, error=error)

    !$omp parallel default(none)  shared(error)
    CALL dbcsr_mm_cannon_lib_init(error)
    !$omp end parallel

    is_initialized = .TRUE.
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_init_lib


! *****************************************************************************
!> \brief Finalize the DBCSR library
!>
!> Cleans up after the DBCSR library.  Used to deallocate persistent objects.
!> \param group ...
!> \param output_unit ...
!> \param[in,out] error     error
! *****************************************************************************
  SUBROUTINE dbcsr_finalize_lib (group, output_unit, error)
    INTEGER, INTENT(IN)                      :: group, output_unit
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_finalize_lib', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!TODO: problem: init/finalize are called by cp2k_runs AND f77_interface

    IF (.NOT. is_initialized) RETURN
    CALL dbcsr_error_set(routineN, error_handle, error)

    IF(output_unit>0)THEN
       WRITE (UNIT=output_unit,FMT="(/,T2,A)") REPEAT("-",79)
       WRITE (UNIT=output_unit,FMT="(T2,A,T80,A)") "-","-"
       WRITE (UNIT=output_unit,FMT="(T2,A,T35,A,T80,A)") "-","DBCSR STATISTICS","-"
       WRITE (UNIT=output_unit,FMT="(T2,A,T80,A)") "-","-"
       WRITE (UNIT=output_unit,FMT="(T2,A)") REPEAT("-",79)
    END IF

    !$omp parallel default(none) shared(output_unit, group, error)
    CALL dbcsr_mm_cannon_lib_finalize(group, output_unit, error)
    !$omp end parallel
    IF(output_unit>0) WRITE (UNIT=output_unit,FMT="(T2,A)") REPEAT("-",79)

    is_initialized = .FALSE.
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_finalize_lib


! *****************************************************************************
!> \brief  Deallocate memory contained in mempools
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_clear_mempools(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    !$omp parallel default(none) shared( error)
    CALL dbcsr_mm_cannon_clear_mempools(error)
    !$omp end parallel
  END SUBROUTINE dbcsr_clear_mempools

! *****************************************************************************
!> \brief ...
!> \param transa ...
!> \param transb ...
!> \param alpha ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param beta ...
!> \param matrix_c ...
!> \param first_row ...
!> \param last_row ...
!> \param first_column ...
!> \param last_column ...
!> \param first_k ...
!> \param last_k ...
!> \param retain_sparsity ...
!> \param filter_eps ...
!> \param error ...
!> \param flop ...
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_s(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, filter_eps,&
       error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    REAL(KIND=real_4), INTENT(IN)            :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    REAL(KIND=real_4), INTENT(IN)            :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(KIND=real_8), INTENT(IN), OPTIONAL  :: filter_eps
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CALL dbcsr_mm_cannon_multiply(transa, transb,&
         dbcsr_scalar(alpha), matrix_a, matrix_b, dbcsr_scalar(beta), matrix_c,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, &
         filter_eps=filter_eps,&
         error=error, flop=flop)
  END SUBROUTINE dbcsr_multiply_s


! *****************************************************************************
!> \brief ...
!> \param transa ...
!> \param transb ...
!> \param alpha ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param beta ...
!> \param matrix_c ...
!> \param first_row ...
!> \param last_row ...
!> \param first_column ...
!> \param last_column ...
!> \param first_k ...
!> \param last_k ...
!> \param retain_sparsity ...
!> \param filter_eps ...
!> \param error ...
!> \param flop ...
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_d(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, filter_eps,&
       error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    REAL(KIND=real_8), INTENT(IN)            :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    REAL(KIND=real_8), INTENT(IN)            :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(KIND=real_8), INTENT(IN), OPTIONAL  :: filter_eps
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_d', &
      routineP = moduleN//':'//routineN

    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4 .AND.&
       dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_4 .AND.&
       dbcsr_get_data_type(matrix_c) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_mm_cannon_multiply(transa, transb,&
            dbcsr_scalar(REAL(alpha,real_4)), matrix_a, matrix_b, &
            dbcsr_scalar(REAL(beta,real_4)), matrix_c,&
            first_row, last_row, first_column, last_column, first_k, last_k,&
            retain_sparsity, &
            filter_eps=filter_eps,&
            error=error, flop=flop)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_8 .AND.&
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_8 .AND.&
           dbcsr_get_data_type(matrix_c) .EQ. dbcsr_type_real_8) THEN
       CALL dbcsr_mm_cannon_multiply(transa, transb,&
            dbcsr_scalar(alpha), matrix_a, matrix_b, dbcsr_scalar(beta), matrix_c,&
            first_row, last_row, first_column, last_column, first_k, last_k,&
            retain_sparsity, &
            filter_eps=filter_eps,&
            error=error, flop=flop)
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_internal_error,&
            routineP, "This combination of data types NYI",__LINE__, error)
    ENDIF
  END SUBROUTINE dbcsr_multiply_d

! *****************************************************************************
!> \brief ...
!> \param transa ...
!> \param transb ...
!> \param alpha ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param beta ...
!> \param matrix_c ...
!> \param first_row ...
!> \param last_row ...
!> \param first_column ...
!> \param last_column ...
!> \param first_k ...
!> \param last_k ...
!> \param retain_sparsity ...
!> \param filter_eps ...
!> \param error ...
!> \param flop ...
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_c(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, filter_eps,&
       error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    COMPLEX(KIND=real_4), INTENT(IN)         :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    COMPLEX(KIND=real_4), INTENT(IN)         :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(KIND=real_8), INTENT(IN), OPTIONAL  :: filter_eps
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CALL dbcsr_mm_cannon_multiply(transa, transb,&
         dbcsr_scalar(alpha), matrix_a, matrix_b, dbcsr_scalar(beta), matrix_c,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, &
         filter_eps=filter_eps,&
         error=error, flop=flop)
  END SUBROUTINE dbcsr_multiply_c


! *****************************************************************************
!> \brief ...
!> \param transa ...
!> \param transb ...
!> \param alpha ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param beta ...
!> \param matrix_c ...
!> \param first_row ...
!> \param last_row ...
!> \param first_column ...
!> \param last_column ...
!> \param first_k ...
!> \param last_k ...
!> \param retain_sparsity ...
!> \param filter_eps ...
!> \param error ...
!> \param flop ...
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_z(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, filter_eps,&
       error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    COMPLEX(KIND=real_8), INTENT(IN)         :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    COMPLEX(KIND=real_8), INTENT(IN)         :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(KIND=real_8), INTENT(IN), OPTIONAL  :: filter_eps
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CALL dbcsr_mm_cannon_multiply(transa, transb,&
         dbcsr_scalar(alpha), matrix_a, matrix_b, dbcsr_scalar(beta), matrix_c,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, &
         filter_eps=filter_eps,&
         error=error, flop=flop)
  END SUBROUTINE dbcsr_multiply_z

END MODULE dbcsr_operations

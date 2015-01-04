!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Facade for matrix multiplication
!> \author  Ole Schuett
!> \date    10-2014
! *****************************************************************************
MODULE dbcsr_multiply_api
  USE dbcsr_data_methods,              ONLY: dbcsr_scalar
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_error_type,&
                                             dbcsr_failure_level,&
                                             dbcsr_internal_error
  USE dbcsr_methods,                   ONLY: dbcsr_get_data_type
  USE dbcsr_mm_cannon,                 ONLY: dbcsr_mm_cannon_multiply,&
                                             dbcsr_multiply_clear_mempools => dbcsr_mm_cannon_clear_mempools,&
                                             dbcsr_multiply_lib_finalize  => dbcsr_mm_cannon_lib_finalize,&
                                             dbcsr_multiply_lib_init => dbcsr_mm_cannon_lib_init
  USE dbcsr_types,                     ONLY: dbcsr_obj,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE kinds,                           ONLY: int_8,&
                                             real_4,&
                                             real_8

  !$ USE OMP_LIB

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_multiply_api'

  PUBLIC :: dbcsr_multiply
  PUBLIC :: dbcsr_multiply_clear_mempools
  PUBLIC :: dbcsr_multiply_lib_finalize, dbcsr_multiply_lib_init

  INTERFACE dbcsr_multiply
     MODULE PROCEDURE dbcsr_mm_cannon_multiply
     MODULE PROCEDURE dbcsr_multiply_s, dbcsr_multiply_d,&
                      dbcsr_multiply_c, dbcsr_multiply_z
  END INTERFACE


 CONTAINS

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

END MODULE dbcsr_multiply_api

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief operations for skinny matrices/vectors expressed in dbcsr form
!> \par History
!>       2014.10 created [Florian Schiffmann]
!> \author Florian Schiffmann
! *****************************************************************************

MODULE dbcsr_vector_operations
  USE dbcsr_data_methods,              ONLY: dbcsr_scalar,&
                                             dbcsr_scalar_fill_all
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_caller_error,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level
  USE dbcsr_types,                     ONLY: dbcsr_obj,&
                                             dbcsr_scalar_type,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE dbcsr_vector_operations_low,     ONLY: create_col_vec_from_matrix,&
                                             create_replicated_col_vec_from_matrix,&
                                             create_replicated_row_vec_from_matrix,&
                                             create_row_vec_from_matrix,&
                                             dbcsr_matrix_colvec_multiply_low
  USE kinds,                           ONLY: real_4,&
                                             real_8

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_vector_operations'

  PUBLIC :: dbcsr_matrix_colvec_multiply,&
            create_col_vec_from_matrix,&
            create_row_vec_from_matrix,&
            create_replicated_col_vec_from_matrix,&
            create_replicated_row_vec_from_matrix

  INTERFACE dbcsr_matrix_colvec_multiply
     MODULE PROCEDURE dbcsr_matrix_colvec_multiply_d, dbcsr_matrix_colvec_multiply_s,&
                      dbcsr_matrix_colvec_multiply_z, dbcsr_matrix_colvec_multiply_c
  END INTERFACE
CONTAINS

! *****************************************************************************
!> \brief multiplies a dbcsr matrix with a column vector like dbcsr matrix.
!>        v_out=beta*v_out+alpha*M*V
!>        IMPORTANT: vector have to be created via the vec create routines:
!>                   create_col_vec_from_matrix,
!>                   create_row_vec_from_matrix,
!>                   create_rep_col_vec_from_matrix,
!>                   create_rep_row_vec_from_matrix
!>        WARNING:   Do not filter the vectors as they are assumed to be non
!>                   sparse in the underlying routines. If your vector is
!>                   sparse, fill it!!!
!> \param matrix a dbcsr matrix
!> \param vec_in the vector to be multiplied (only available on proc_col 0)
!> \param vec_out the result vector (only available on proc_col 0)
!> \param alpha  as described in formula
!> \param beta  as described in formula
!> \param work_row a work row vector replicated on all proc_cols.
!> \param work_col a work col vector replicated on all proc_rows.
!> \param error ...
! *****************************************************************************

  SUBROUTINE dbcsr_matrix_colvec_multiply_d(matrix,vec_in,vec_out,alpha,beta,work_row,work_col,error)
    TYPE(dbcsr_obj)                          :: matrix, vec_in, vec_out
    REAL(real_8)                             :: alpha, beta
    TYPE(dbcsr_obj)                          :: work_row, work_col
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'dbcsr_matrix_colvec_multiply_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: matrix_type
    TYPE(dbcsr_scalar_type)                  :: a, b

    matrix_type=matrix%m%data_type
    SELECT CASE (matrix_type)
    CASE (dbcsr_type_real_8,dbcsr_type_real_4,dbcsr_type_complex_4,dbcsr_type_complex_8)
        ! We can convert the prefactors if necessary
    CASE DEFAULT
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Unknown data type in matrix",__LINE__,error)
    END SELECT

    a=dbcsr_scalar(alpha); b=dbcsr_scalar(beta)
    CALL dbcsr_scalar_fill_all (a); CALL dbcsr_scalar_fill_all (b)

    CALL dbcsr_matrix_colvec_multiply_low(matrix,vec_in,vec_out,a,b,work_row,work_col,error)

  END SUBROUTINE dbcsr_matrix_colvec_multiply_d

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
  SUBROUTINE dbcsr_matrix_colvec_multiply_s(matrix,vec_in,vec_out,alpha,beta,work_row,work_col,error)
    TYPE(dbcsr_obj)                          :: matrix, vec_in, vec_out
    REAL(real_4)                             :: alpha, beta
    TYPE(dbcsr_obj)                          :: work_row, work_col
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'dbcsr_matrix_colvec_multiply_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: matrix_type
    TYPE(dbcsr_scalar_type)                  :: a, b

    matrix_type=matrix%m%data_type
    SELECT CASE (matrix_type)
    CASE (dbcsr_type_real_8,dbcsr_type_real_4,dbcsr_type_complex_4,dbcsr_type_complex_8)
        ! We can convert the prefactors if necessary
    CASE DEFAULT
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Unknown data type in matrix",__LINE__,error)
    END SELECT

    a=dbcsr_scalar(alpha); b=dbcsr_scalar(beta)
    CALL dbcsr_scalar_fill_all (a); CALL dbcsr_scalar_fill_all (b)

    CALL dbcsr_matrix_colvec_multiply_low(matrix,vec_in,vec_out,a,b,work_row,work_col,error)

  END SUBROUTINE dbcsr_matrix_colvec_multiply_s

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
  SUBROUTINE dbcsr_matrix_colvec_multiply_c(matrix,vec_in,vec_out,alpha,beta,work_row,work_col,error)
    TYPE(dbcsr_obj)                          :: matrix, vec_in, vec_out
    COMPLEX(real_4)                          :: alpha, beta
    TYPE(dbcsr_obj)                          :: work_row, work_col
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'dbcsr_matrix_colvec_multiply_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: matrix_type
    TYPE(dbcsr_scalar_type)                  :: a, b

    matrix_type=matrix%m%data_type
    SELECT CASE (matrix_type)
    CASE (dbcsr_type_real_8,dbcsr_type_real_4)
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Complex prefactors dissallowed for real matrices",__LINE__,error)
    CASE (dbcsr_type_complex_4,dbcsr_type_complex_8)
        ! We can convert the prefactors if necessary
    CASE DEFAULT
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Unknown data type in matrix",__LINE__,error)
    END SELECT  
    
    a=dbcsr_scalar(alpha); b=dbcsr_scalar(beta)
    CALL dbcsr_scalar_fill_all (a); CALL dbcsr_scalar_fill_all (b)

    CALL dbcsr_matrix_colvec_multiply_low(matrix,vec_in,vec_out,a,b,work_row,work_col,error)

  END SUBROUTINE dbcsr_matrix_colvec_multiply_c

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
  SUBROUTINE dbcsr_matrix_colvec_multiply_z(matrix,vec_in,vec_out,alpha,beta,work_row,work_col,error)
    TYPE(dbcsr_obj)                          :: matrix, vec_in, vec_out
    COMPLEX(real_8)                          :: alpha, beta
    TYPE(dbcsr_obj)                          :: work_row, work_col
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'dbcsr_matrix_colvec_multiply_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: matrix_type
    TYPE(dbcsr_scalar_type)                  :: a, b

    matrix_type=matrix%m%data_type
    SELECT CASE (matrix_type)
    CASE (dbcsr_type_real_8,dbcsr_type_real_4)
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Complex prefactors dissallowed for real matrices",__LINE__,error)
    CASE (dbcsr_type_complex_4,dbcsr_type_complex_8)
        ! We can convert the prefactors if necessary
    CASE DEFAULT
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Unknown data type in matrix",__LINE__,error)
    END SELECT

    a=dbcsr_scalar(alpha); b=dbcsr_scalar(beta)
    CALL dbcsr_scalar_fill_all (a); CALL dbcsr_scalar_fill_all (b)

    CALL dbcsr_matrix_colvec_multiply_low(matrix,vec_in,vec_out,a,b,work_row,work_col,error)

  END SUBROUTINE dbcsr_matrix_colvec_multiply_z

END MODULE dbcsr_vector_operations

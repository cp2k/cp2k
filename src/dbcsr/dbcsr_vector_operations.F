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

MODULE dbcsr_vector_operations
  USE array_types,                     ONLY: array_i1d_obj,&
                                             array_new,&
                                             array_release
  USE dbcsr_block_access,              ONLY: dbcsr_reserve_all_blocks
  USE dbcsr_data_methods,              ONLY: dbcsr_scalar,&
                                             dbcsr_scalar_fill_all
  USE dbcsr_dist_methods,              ONLY: dbcsr_distribution_col_dist,&
                                             dbcsr_distribution_mp,&
                                             dbcsr_distribution_new,&
                                             dbcsr_distribution_release,&
                                             dbcsr_distribution_row_dist
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_caller_error,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level
  USE dbcsr_methods,                   ONLY: dbcsr_col_block_sizes,&
                                             dbcsr_distribution,&
                                             dbcsr_get_info,&
                                             dbcsr_init,&
                                             dbcsr_row_block_sizes
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_obj,&
                                             dbcsr_scalar_type,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_no_symmetry,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE dbcsr_vector_operations_low,     ONLY: dbcsr_matrix_colvec_multiply_low
  USE dbcsr_work_operations,           ONLY: dbcsr_create
  USE kinds,                           ONLY: real_4,&
                                             real_8

  IMPLICIT NONE

  PRIVATE

  TYPE block_ptr_d
    REAL(real_8),DIMENSION(:,:), POINTER          :: ptr =>NULL()
  END TYPE
  TYPE block_ptr_s
    REAL(real_4),DIMENSION(:,:), POINTER          :: ptr =>NULL()
  END TYPE
  TYPE block_ptr_c
    COMPLEX(real_4),DIMENSION(:,:), POINTER       :: ptr =>NULL()
  END TYPE
  TYPE block_ptr_z
    COMPLEX(real_8),DIMENSION(:,:), POINTER       :: ptr =>NULL()
  END TYPE

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

! *****************************************************************************
!> \brief creates a dbcsr col vector like object which lives on proc_col 0
!>        and has the same row dist as the template matrix
!>        the returned matrix is fully allocated and all blocks are set to 0
!>        this is not a sparse object (and must never be)
!> \param dbcsr_vec  the vector object to create must be allocated but not initialized
!> \param matrix a dbcsr matrix used as template
!> \param ncol number of vectors in the dbcsr_object (1 for vector, n for skinny matrix)
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_col_vec_from_matrix(dbcsr_vec,matrix,ncol,error) 
    TYPE(dbcsr_obj)                          :: dbcsr_vec, matrix
    INTEGER                                  :: ncol
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_col_vec_from_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: npcol
    INTEGER, DIMENSION(:), POINTER           :: col_sizes, dist_col
    TYPE(array_i1d_obj)                      :: bin_dist, col_size
    TYPE(dbcsr_distribution_obj)             :: dist_col_vec, distri

    CALL dbcsr_get_info(matrix=matrix,distribution=distri)
    npcol=SIZE(distri%d%mp_env%mp%pgrid,2)
    ALLOCATE(dist_col(1));  ALLOCATE(col_sizes(1))
    dist_col=0; col_sizes=ncol

    CALL array_new (bin_dist, dist_col, gift=.FALSE.)
    CALL array_new (col_size, col_sizes, gift=.FALSE.)

    CALL dbcsr_distribution_new (dist_col_vec,&
         dbcsr_distribution_mp (dbcsr_distribution(matrix)),&
         dbcsr_distribution_row_dist (dbcsr_distribution(matrix)),&
         bin_dist)
    CALL dbcsr_init(dbcsr_vec)
    CALL dbcsr_create(dbcsr_vec, "D", dist_col_vec, dbcsr_type_no_symmetry, &
            dbcsr_row_block_sizes(matrix), col_size,&
            0, 0, data_type=matrix%m%data_type, error=error)    

    CALL dbcsr_reserve_all_blocks(dbcsr_vec,error)    

    CALL dbcsr_distribution_release(dist_col_vec)
    DEALLOCATE(dist_col,col_sizes)
    CALL array_release(bin_dist)
    CALL array_release(col_size)

  END SUBROUTINE create_col_vec_from_matrix

! *****************************************************************************
!> \brief creates a dbcsr row vector like object which lives on proc_row 0
!>        and has the same row dist as the template matrix
!>        the returned matrix is fully allocated and all blocks are set to 0
!>        this is not a sparse object (and must never be)
!> \param dbcsr_vec ...
!> \param matrix a dbcsr matrix used as template
!> \param nrow number of vectors in the dbcsr_object (1 for vector, n for skinny matrix)
!> \param error ... 
! *****************************************************************************
  SUBROUTINE create_row_vec_from_matrix(dbcsr_vec,matrix,nrow,error)
    TYPE(dbcsr_obj)                          :: dbcsr_vec, matrix
    INTEGER                                  :: nrow
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_row_vec_from_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nprow
    INTEGER, DIMENSION(:), POINTER           :: dist_row, row_sizes
    TYPE(array_i1d_obj)                      :: bin_dist, row_size
    TYPE(dbcsr_distribution_obj)             :: dist_row_vec, distri

    CALL dbcsr_get_info(matrix=matrix,distribution=distri)
    nprow=SIZE(distri%d%mp_env%mp%pgrid,1)
    ALLOCATE(dist_row(1)); ALLOCATE(row_sizes(1))
    dist_row=0; row_sizes=nrow
    CALL array_new (bin_dist, dist_row, gift=.FALSE.)
    CALL array_new (row_size, row_sizes, gift=.FALSE.)

    CALL dbcsr_distribution_new (dist_row_vec,&
         dbcsr_distribution_mp (dbcsr_distribution(matrix)),&
         bin_dist,&
         dbcsr_distribution_col_dist (dbcsr_distribution(matrix)))
    
    CALL dbcsr_init(dbcsr_vec)
    CALL dbcsr_create(dbcsr_vec, "D", dist_row_vec, dbcsr_type_no_symmetry, &
            row_size, dbcsr_col_block_sizes(matrix), &
            0, 0, data_type=matrix%m%data_type, error=error)    

    CALL dbcsr_reserve_all_blocks(dbcsr_vec,error) 

    CALL dbcsr_distribution_release(dist_row_vec)
    DEALLOCATE(dist_row,row_sizes)
    CALL array_release(bin_dist)
    CALL array_release(row_size)

  END SUBROUTINE create_row_vec_from_matrix

! *****************************************************************************
!> \brief creates a col vector like object whose blocks can be replicated 
!>        along the processor row and has the same row dist as the template matrix
!>        the returned matrix is fully allocated and all blocks are set to 0
!>        this is not a sparse object (and must never be)
!> \param dbcsr_vec the vector object to create must be allocated but not initialized
!> \param matrix a dbcsr matrix used as template
!> \param ncol number of vectors in the dbcsr_object (1 for vector, n for skinny matrix)
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_replicated_col_vec_from_matrix(dbcsr_vec,matrix,ncol,error)
    TYPE(dbcsr_obj)                          :: dbcsr_vec, matrix
    INTEGER                                  :: ncol
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'create_replicated_col_vec_from_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, npcol
    INTEGER, DIMENSION(:), POINTER           :: col_sizes, dist_col
    TYPE(array_i1d_obj)                      :: bin_dist, col_size
    TYPE(dbcsr_distribution_obj)             :: dist_col_vec, distri

    CALL dbcsr_get_info(matrix=matrix,distribution=distri)
    npcol=SIZE(distri%d%mp_env%mp%pgrid,2)
    ALLOCATE(dist_col(npcol));  ALLOCATE(col_sizes(npcol))
    col_sizes=ncol
    DO i=0,npcol-1
      dist_col(i+1)=i
    END DO

    CALL array_new (bin_dist, dist_col, gift=.FALSE.)
    CALL array_new (col_size, col_sizes, gift=.FALSE.)

    CALL dbcsr_distribution_new (dist_col_vec,&
         dbcsr_distribution_mp (dbcsr_distribution(matrix)),&
         dbcsr_distribution_row_dist (dbcsr_distribution(matrix)),&
         bin_dist)

    CALL dbcsr_init(dbcsr_vec)
    CALL dbcsr_create(dbcsr_vec, "D", dist_col_vec, dbcsr_type_no_symmetry, &
            dbcsr_row_block_sizes(matrix), col_size,&
            0, 0, data_type=matrix%m%data_type, error=error)

    CALL dbcsr_reserve_all_blocks(dbcsr_vec,error)

    CALL dbcsr_distribution_release(dist_col_vec)
    DEALLOCATE(dist_col,col_sizes)
    CALL array_release(bin_dist)
    CALL array_release(col_size)

  END SUBROUTINE create_replicated_col_vec_from_matrix

! *****************************************************************************
!> \brief creates a row vector like object whose blocks can be replicated 
!>        along the processor col and has the same col dist as the template matrix
!>        the returned matrix is fully allocated and all blocks are set to 0
!>        this is not a sparse object (and must never be)
!> \param dbcsr_vec the vector object to create must be allocated but not initialized
!> \param matrix a dbcsr matrix used as template
!> \param nrow number of vectors in the dbcsr_object (1 for vector, n for skinny matrix)
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_replicated_row_vec_from_matrix(dbcsr_vec,matrix,nrow,error)
    TYPE(dbcsr_obj)                          :: dbcsr_vec, matrix
    INTEGER                                  :: nrow
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'create_replicated_row_vec_from_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nprow
    INTEGER, DIMENSION(:), POINTER           :: dist_row, row_sizes
    TYPE(array_i1d_obj)                      :: bin_dist, row_size
    TYPE(dbcsr_distribution_obj)             :: dist_row_vec, distri

    CALL dbcsr_get_info(matrix=matrix,distribution=distri)
    nprow=SIZE(distri%d%mp_env%mp%pgrid,1)
    ALLOCATE(dist_row(nprow)); ALLOCATE(row_sizes(nprow))
    row_sizes=nrow
    DO i=0,nprow-1
      dist_row(i+1)=i
    END DO
 
    CALL array_new (bin_dist, dist_row, gift=.FALSE.)
    CALL array_new (row_size, row_sizes, gift=.FALSE.)

    CALL dbcsr_distribution_new (dist_row_vec,&
         dbcsr_distribution_mp (dbcsr_distribution(matrix)),&
         bin_dist,&
         dbcsr_distribution_col_dist (dbcsr_distribution(matrix)))
   
    CALL dbcsr_init(dbcsr_vec)
    CALL dbcsr_create(dbcsr_vec, "D", dist_row_vec, dbcsr_type_no_symmetry, &
            row_size, dbcsr_col_block_sizes(matrix), &
            0, 0, data_type=matrix%m%data_type, error=error)   

    CALL dbcsr_reserve_all_blocks(dbcsr_vec,error)

    CALL dbcsr_distribution_release(dist_row_vec)
    DEALLOCATE(dist_row,row_sizes)
    CALL array_release(bin_dist)
    CALL array_release(row_size)

  END SUBROUTINE create_replicated_row_vec_from_matrix

END MODULE dbcsr_vector_operations

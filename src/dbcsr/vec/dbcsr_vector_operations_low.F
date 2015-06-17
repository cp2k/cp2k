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

MODULE dbcsr_vector_operations_low
  USE dbcsr_block_access,              ONLY: dbcsr_reserve_all_blocks
  USE dbcsr_data_methods,              ONLY: dbcsr_get_data_p
  USE dbcsr_dist_methods,              ONLY: dbcsr_distribution_col_clusters,&
                                             dbcsr_distribution_col_dist,&
                                             dbcsr_distribution_mp,&
                                             dbcsr_distribution_new,&
                                             dbcsr_distribution_release,&
                                             dbcsr_distribution_row_clusters,&
                                             dbcsr_distribution_row_dist
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_caller_error,&
                                             dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_methods,                   ONLY: dbcsr_col_block_sizes,&
                                             dbcsr_distribution,&
                                             dbcsr_get_matrix_type,&
                                             dbcsr_init,&
                                             dbcsr_release,&
                                             dbcsr_row_block_sizes
  USE dbcsr_operations,                ONLY: dbcsr_copy,&
                                             dbcsr_get_info,&
                                             dbcsr_set
  USE dbcsr_types,                     ONLY: &
       dbcsr_distribution_obj, dbcsr_iterator, dbcsr_obj, dbcsr_scalar_type, &
       dbcsr_type_antisymmetric, dbcsr_type_complex_4, dbcsr_type_complex_8, &
       dbcsr_type_no_symmetry, dbcsr_type_real_4, dbcsr_type_real_8, &
       dbcsr_type_symmetric
  USE dbcsr_work_operations,           ONLY: dbcsr_create
  USE kinds,                           ONLY: real_4,&
                                             real_8
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum

  !$ USE OMP_LIB

  IMPLICIT NONE

! the following types provide fast access to distributed dbcsr vectors
INCLUDE 'hash_table_types.f90'

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

  TYPE fast_vec_access_type
    TYPE(hash_table_type) :: hash_table
    TYPE(block_ptr_d), DIMENSION(:), ALLOCATABLE :: blk_map_d
    TYPE(block_ptr_s), DIMENSION(:), ALLOCATABLE :: blk_map_s
    TYPE(block_ptr_c), DIMENSION(:), ALLOCATABLE :: blk_map_c
    TYPE(block_ptr_z), DIMENSION(:), ALLOCATABLE :: blk_map_z
  END TYPE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_vector_operations_low'

  PUBLIC :: dbcsr_matrix_colvec_multiply_low, &
            create_col_vec_from_matrix,&
            create_row_vec_from_matrix,&
            create_replicated_col_vec_from_matrix,&
            create_replicated_row_vec_from_matrix

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

    INTEGER                                  :: handle, matrix_type

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

    INTEGER                                  :: handle, npcol
    INTEGER, DIMENSION(:), POINTER           :: col_sizes, dist_col
    TYPE(dbcsr_distribution_obj)             :: dist_col_vec, distri

    CALL dbcsr_error_set(routineN, handle, error)

    CALL dbcsr_get_info(matrix=matrix,distribution=distri)
    npcol=SIZE(distri%d%mp_env%mp%pgrid,2)
    ALLOCATE(dist_col(1));  ALLOCATE(col_sizes(1))
    dist_col=0; col_sizes=ncol

    CALL dbcsr_distribution_new (dist_col_vec,&
         dbcsr_distribution_mp (dbcsr_distribution(matrix)),&
         dbcsr_distribution_row_dist (dbcsr_distribution(matrix)),&
         dist_col,&
         row_dist_cluster=dbcsr_distribution_row_clusters (dbcsr_distribution(matrix)))
    CALL dbcsr_init(dbcsr_vec)
    CALL dbcsr_create(dbcsr_vec, "D", dist_col_vec, dbcsr_type_no_symmetry, &
            dbcsr_row_block_sizes(matrix), col_sizes,&
            0, 0, data_type=matrix%m%data_type, error=error)    

    CALL dbcsr_reserve_all_blocks(dbcsr_vec,error)    

    CALL dbcsr_distribution_release(dist_col_vec)
    DEALLOCATE(dist_col,col_sizes)

    CALL dbcsr_error_stop(handle, error)

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

    INTEGER                                  :: handle, nprow
    INTEGER, DIMENSION(:), POINTER           :: dist_row, row_sizes
    TYPE(dbcsr_distribution_obj)             :: dist_row_vec, distri

    CALL dbcsr_error_set(routineN, handle, error)

    CALL dbcsr_get_info(matrix=matrix,distribution=distri)
    nprow=SIZE(distri%d%mp_env%mp%pgrid,1)
    ALLOCATE(dist_row(1)); ALLOCATE(row_sizes(1))
    dist_row=0; row_sizes=nrow

    CALL dbcsr_distribution_new (dist_row_vec,&
         dbcsr_distribution_mp (dbcsr_distribution(matrix)),&
         dist_row,&
         dbcsr_distribution_col_dist (dbcsr_distribution(matrix)),&
         col_dist_cluster=dbcsr_distribution_col_clusters (dbcsr_distribution(matrix)))
    
    CALL dbcsr_init(dbcsr_vec)
    CALL dbcsr_create(dbcsr_vec, "D", dist_row_vec, dbcsr_type_no_symmetry, &
            row_sizes, dbcsr_col_block_sizes(matrix), &
            0, 0, data_type=matrix%m%data_type, error=error)    

    CALL dbcsr_reserve_all_blocks(dbcsr_vec,error) 

    CALL dbcsr_distribution_release(dist_row_vec)
    DEALLOCATE(dist_row,row_sizes)

    CALL dbcsr_error_stop(handle, error)

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

    INTEGER                                  :: handle, i, npcol
    INTEGER, DIMENSION(:), POINTER           :: col_sizes, dist_col
    TYPE(dbcsr_distribution_obj)             :: dist_col_vec, distri

    CALL dbcsr_error_set(routineN, handle, error)

    CALL dbcsr_get_info(matrix=matrix,distribution=distri)
    npcol=SIZE(distri%d%mp_env%mp%pgrid,2)
    ALLOCATE(dist_col(npcol));  ALLOCATE(col_sizes(npcol))
    col_sizes=ncol
    DO i=0,npcol-1
      dist_col(i+1)=i
    END DO

    CALL dbcsr_distribution_new (dist_col_vec,&
         dbcsr_distribution_mp (dbcsr_distribution(matrix)),&
         dbcsr_distribution_row_dist (dbcsr_distribution(matrix)),&
         dist_col,&
         row_dist_cluster=dbcsr_distribution_row_clusters (dbcsr_distribution(matrix)))

    CALL dbcsr_init(dbcsr_vec)
    CALL dbcsr_create(dbcsr_vec, "D", dist_col_vec, dbcsr_type_no_symmetry, &
            dbcsr_row_block_sizes(matrix), col_sizes,&
            0, 0, data_type=matrix%m%data_type, error=error)

    CALL dbcsr_reserve_all_blocks(dbcsr_vec,error)

    CALL dbcsr_distribution_release(dist_col_vec)
    DEALLOCATE(dist_col,col_sizes)

    CALL dbcsr_error_stop(handle, error)

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

    INTEGER                                  :: handle, i, nprow
    INTEGER, DIMENSION(:), POINTER           :: dist_row, row_sizes
    TYPE(dbcsr_distribution_obj)             :: dist_row_vec, distri

    CALL dbcsr_error_set(routineN, handle, error)

    CALL dbcsr_get_info(matrix=matrix,distribution=distri)
    nprow=SIZE(distri%d%mp_env%mp%pgrid,1)
    ALLOCATE(dist_row(nprow)); ALLOCATE(row_sizes(nprow))
    row_sizes=nrow
    DO i=0,nprow-1
      dist_row(i+1)=i
    END DO
 
    CALL dbcsr_distribution_new (dist_row_vec,&
         dbcsr_distribution_mp (dbcsr_distribution(matrix)),&
         dist_row,&
         dbcsr_distribution_col_dist (dbcsr_distribution(matrix)),&
         col_dist_cluster=dbcsr_distribution_col_clusters (dbcsr_distribution(matrix)))
   
    CALL dbcsr_init(dbcsr_vec)
    CALL dbcsr_create(dbcsr_vec, "D", dist_row_vec, dbcsr_type_no_symmetry, &
            row_sizes, dbcsr_col_block_sizes(matrix), &
            0, 0, data_type=matrix%m%data_type, error=error)   

    CALL dbcsr_reserve_all_blocks(dbcsr_vec,error)

    CALL dbcsr_distribution_release(dist_row_vec)
    DEALLOCATE(dist_row,row_sizes)

    CALL dbcsr_error_stop(handle, error)

  END SUBROUTINE create_replicated_row_vec_from_matrix

! *****************************************************************************
!> \brief given a column vector, prepare the fast_vec_access container
!> \param vec ...
!> \param fast_vec_access ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_fast_col_vec_access(vec, fast_vec_access, error)
    TYPE(dbcsr_obj)                          :: vec
    TYPE(fast_vec_access_type)               :: fast_vec_access
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_fast_col_vec_access', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL dbcsr_error_set(routineN, handle, error)

    SELECT CASE(vec%m%data_type)
    CASE(dbcsr_type_real_8)
      CALL create_fast_col_vec_access_d(vec, fast_vec_access, error)
    CASE(dbcsr_type_real_4)
      CALL create_fast_col_vec_access_s(vec, fast_vec_access, error)
    CASE(dbcsr_type_complex_8)
      CALL create_fast_col_vec_access_z(vec, fast_vec_access, error)
    CASE(dbcsr_type_complex_4)
      CALL create_fast_col_vec_access_c(vec, fast_vec_access, error)
    END SELECT

    CALL dbcsr_error_stop(handle, error)

  END SUBROUTINE create_fast_col_vec_access

! *****************************************************************************
!> \brief given a row vector, prepare the fast_vec_access_container
!> \param vec ...
!> \param fast_vec_access ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_fast_row_vec_access(vec, fast_vec_access, error)
    TYPE(dbcsr_obj)                          :: vec
    TYPE(fast_vec_access_type)               :: fast_vec_access
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_fast_row_vec_access', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL dbcsr_error_set(routineN, handle, error)

    SELECT CASE(vec%m%data_type)
    CASE(dbcsr_type_real_8)
       CALL create_fast_row_vec_access_d(vec, fast_vec_access, error)
    CASE(dbcsr_type_real_4)
       CALL create_fast_row_vec_access_s(vec, fast_vec_access, error)
    CASE(dbcsr_type_complex_8)
       CALL create_fast_row_vec_access_c(vec, fast_vec_access, error)
    CASE(dbcsr_type_complex_4)
       CALL create_fast_row_vec_access_z(vec, fast_vec_access, error)
    END SELECT

    CALL dbcsr_error_stop(handle, error)

  END SUBROUTINE create_fast_row_vec_access

! *****************************************************************************
!> \brief release all memory associated with the fast_vec_access type
!> \param fast_vec_access ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE release_fast_vec_access(fast_vec_access, error)
    TYPE(fast_vec_access_type)               :: fast_vec_access
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'release_fast_vec_access', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL dbcsr_error_set(routineN, handle, error)

    CALL hash_table_release(fast_vec_access%hash_table)
    IF (ALLOCATED(fast_vec_access%blk_map_d)) DEALLOCATE(fast_vec_access%blk_map_d)
    IF (ALLOCATED(fast_vec_access%blk_map_s)) DEALLOCATE(fast_vec_access%blk_map_s)
    IF (ALLOCATED(fast_vec_access%blk_map_c)) DEALLOCATE(fast_vec_access%blk_map_c)
    IF (ALLOCATED(fast_vec_access%blk_map_z)) DEALLOCATE(fast_vec_access%blk_map_z)

    CALL dbcsr_error_stop(handle, error)

  END SUBROUTINE release_fast_vec_access

INCLUDE 'hash_table.f90'

#include "dbcsr_vector_operations_low_d.f90"
#include "dbcsr_vector_operations_low_z.f90"
#include "dbcsr_vector_operations_low_s.f90"
#include "dbcsr_vector_operations_low_c.f90"

END MODULE dbcsr_vector_operations_low

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief operations for skinny matrices/vectors expressed in dbcsr form
!> \par History
!>       2014.10 created [Florian Schiffmann]
!> \author Florian Schiffmann
! **************************************************************************************************

MODULE dbcsr_vector
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_distribution_get,&
                                              dbcsr_distribution_new,&
                                              dbcsr_distribution_release,&
                                              dbcsr_distribution_type,&
                                              dbcsr_get_info,&
                                              dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop,&
                                              dbcsr_iterator_type,&
                                              dbcsr_release,&
                                              dbcsr_reserve_all_blocks,&
                                              dbcsr_set,dbcsr_get_data_p,&
                                              dbcsr_type,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_complex_4,&
                                              dbcsr_type_complex_4,&
                                              dbcsr_type_complex_8,&
                                              dbcsr_type_complex_8,&
                                              dbcsr_type_no_symmetry,&
                                              dbcsr_type_real_4,&
                                              dbcsr_type_real_4,&
                                              dbcsr_type_real_8,&
                                              dbcsr_type_real_8,&
                                              dbcsr_type_symmetric
   USE kinds,                           ONLY: real_4,&
                                              real_8
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_sum


#include "../base/base_uses.f90"

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_vector_operations'

!    ! the following types provide fast access to distributed dbcsr vectors
#include "hash_table_types.f90"

   TYPE block_ptr_d
      REAL(real_8), DIMENSION(:, :), POINTER          :: ptr => NULL()
      INTEGER                                         :: assigned_thread
   END TYPE
   TYPE block_ptr_s
      REAL(real_4), DIMENSION(:, :), POINTER          :: ptr => NULL()
      INTEGER                                         :: assigned_thread
   END TYPE
   TYPE block_ptr_c
      COMPLEX(real_4), DIMENSION(:, :), POINTER       :: ptr => NULL()
      INTEGER                                         :: assigned_thread
   END TYPE
   TYPE block_ptr_z
      COMPLEX(real_8), DIMENSION(:, :), POINTER       :: ptr => NULL()
      INTEGER                                         :: assigned_thread
   END TYPE

   TYPE fast_vec_access_type
      TYPE(hash_table_type) :: hash_table
      TYPE(block_ptr_d), DIMENSION(:), ALLOCATABLE :: blk_map_d
      TYPE(block_ptr_s), DIMENSION(:), ALLOCATABLE :: blk_map_s
      TYPE(block_ptr_c), DIMENSION(:), ALLOCATABLE :: blk_map_c
      TYPE(block_ptr_z), DIMENSION(:), ALLOCATABLE :: blk_map_z
   END TYPE

   PUBLIC :: dbcsr_matrix_colvec_multiply, &
             create_col_vec_from_matrix, &
             create_row_vec_from_matrix, &
             create_replicated_col_vec_from_matrix, &
             create_replicated_row_vec_from_matrix

   INTERFACE dbcsr_matrix_colvec_multiply
      MODULE PROCEDURE dbcsr_matrix_colvec_multiply_d
      MODULE PROCEDURE dbcsr_matrix_colvec_multiply_s
      MODULE PROCEDURE dbcsr_matrix_colvec_multiply_z
      MODULE PROCEDURE dbcsr_matrix_colvec_multiply_c
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief creates a dbcsr col vector like object which lives on proc_col 0
!>        and has the same row dist as the template matrix
!>        the returned matrix is fully allocated and all blocks are set to 0
!>        this is not a sparse object (and must never be)
!> \param dbcsr_vec  the vector object to create must be allocated but not initialized
!> \param matrix a dbcsr matrix used as template
!> \param ncol number of vectors in the dbcsr_typeect (1 for vector, n for skinny matrix)
! **************************************************************************************************
   SUBROUTINE create_col_vec_from_matrix(dbcsr_vec, matrix, ncol)
      TYPE(dbcsr_type)                                   :: dbcsr_vec, matrix
      INTEGER                                            :: ncol

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_col_vec_from_matrix'

      INTEGER                                            :: handle, npcols, data_type
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size, col_blk_size, row_dist, col_dist
      TYPE(dbcsr_distribution_type)                      :: dist_col_vec, dist

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix, data_type=data_type, row_blk_size=row_blk_size, distribution=dist)
      CALL dbcsr_distribution_get(dist, npcols=npcols, row_dist=row_dist)

      ALLOCATE (col_dist(1), col_blk_size(1))
      col_dist(1) = 0
      col_blk_size(1) = ncol
      CALL dbcsr_distribution_new(dist_col_vec, template=dist, row_dist=row_dist, col_dist=col_dist)

      CALL dbcsr_create(dbcsr_vec, "D", dist_col_vec,&
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=row_blk_size,&
                        col_blk_size=col_blk_size, &
                        nze=0, data_type=data_type)
      CALL dbcsr_reserve_all_blocks(dbcsr_vec)

      CALL dbcsr_distribution_release(dist_col_vec)
      DEALLOCATE (col_dist, col_blk_size)
      CALL timestop(handle)

   END SUBROUTINE create_col_vec_from_matrix

! **************************************************************************************************
!> \brief creates a dbcsr row vector like object which lives on proc_row 0
!>        and has the same row dist as the template matrix
!>        the returned matrix is fully allocated and all blocks are set to 0
!>        this is not a sparse object (and must never be)
!> \param dbcsr_vec ...
!> \param matrix a dbcsr matrix used as template
!> \param nrow number of vectors in the dbcsr_typeect (1 for vector, n for skinny matrix)
! **************************************************************************************************
SUBROUTINE create_row_vec_from_matrix(dbcsr_vec, matrix, nrow)
      TYPE(dbcsr_type)                                   :: dbcsr_vec, matrix
      INTEGER                                            :: nrow

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_row_vec_from_matrix'

      INTEGER                                            :: handle, nprows, data_type
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size, col_blk_size, row_dist, col_dist
      TYPE(dbcsr_distribution_type)                      :: dist_row_vec, dist

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix, data_type=data_type, col_blk_size=col_blk_size, distribution=dist)
      CALL dbcsr_distribution_get(dist, nprows=nprows, col_dist=col_dist)

      ALLOCATE (row_dist(1), row_blk_size(1))
      row_dist(1) = 0
      row_blk_size(1) = nrow
      CALL dbcsr_distribution_new(dist_row_vec, template=dist, row_dist=row_dist, col_dist=col_dist)

      CALL dbcsr_create(dbcsr_vec, "D", dist_row_vec,&
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=row_blk_size,&
                        col_blk_size=col_blk_size, &
                        nze=0, data_type=data_type)
      CALL dbcsr_reserve_all_blocks(dbcsr_vec)

      CALL dbcsr_distribution_release(dist_row_vec)
      DEALLOCATE (row_dist, row_blk_size)
      CALL timestop(handle)

   END SUBROUTINE create_row_vec_from_matrix


! **************************************************************************************************
!> \brief creates a col vector like object whose blocks can be replicated
!>        along the processor row and has the same row dist as the template matrix
!>        the returned matrix is fully allocated and all blocks are set to 0
!>        this is not a sparse object (and must never be)
!> \param dbcsr_vec the vector object to create must be allocated but not initialized
!> \param matrix a dbcsr matrix used as template
!> \param ncol number of vectors in the dbcsr_typeect (1 for vector, n for skinny matrix)
! **************************************************************************************************
   SUBROUTINE create_replicated_col_vec_from_matrix(dbcsr_vec, matrix, ncol)
      TYPE(dbcsr_type)                                   :: dbcsr_vec, matrix
      INTEGER                                            :: ncol

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_replicated_col_vec_from_matrix'

      INTEGER                                            :: handle, npcols, data_type, i
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size, col_blk_size, row_dist, col_dist
      TYPE(dbcsr_distribution_type)                      :: dist_col_vec, dist
      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix, data_type=data_type, row_blk_size=row_blk_size, distribution=dist)
      CALL dbcsr_distribution_get(dist, npcols=npcols, row_dist=row_dist)

      ALLOCATE (col_dist(npcols), col_blk_size(npcols))
      col_blk_size(:) = ncol
      DO i = 0, npcols-1
         col_dist(i+1) = i
      END DO
      CALL dbcsr_distribution_new(dist_col_vec, template=dist, row_dist=row_dist, col_dist=col_dist)

      CALL dbcsr_create(dbcsr_vec, "D", dist_col_vec, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=row_blk_size,&
                        col_blk_size=col_blk_size, &
                        nze=0, data_type=data_type)
      CALL dbcsr_reserve_all_blocks(dbcsr_vec)

      CALL dbcsr_distribution_release(dist_col_vec)
      DEALLOCATE (col_dist, col_blk_size)
      CALL timestop(handle)

   END SUBROUTINE create_replicated_col_vec_from_matrix

! **************************************************************************************************
!> \brief creates a row vector like object whose blocks can be replicated
!>        along the processor col and has the same col dist as the template matrix
!>        the returned matrix is fully allocated and all blocks are set to 0
!>        this is not a sparse object (and must never be)
!> \param dbcsr_vec the vector object to create must be allocated but not initialized
!> \param matrix a dbcsr matrix used as template
!> \param nrow number of vectors in the dbcsr_typeect (1 for vector, n for skinny matrix)
! **************************************************************************************************
   SUBROUTINE create_replicated_row_vec_from_matrix(dbcsr_vec, matrix, nrow)
      TYPE(dbcsr_type)                                   :: dbcsr_vec
      TYPE(dbcsr_type)                                   :: matrix
      INTEGER                                            :: nrow

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_replicated_row_vec_from_matrix'

      INTEGER                                            :: handle, i, nprows, data_type
      INTEGER, DIMENSION(:), POINTER                     :: row_dist, col_dist, row_blk_size, col_blk_size
      TYPE(dbcsr_distribution_type)                      :: dist_row_vec, dist

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix, distribution=dist, col_blk_size=col_blk_size, data_type=data_type)
      CALL dbcsr_distribution_get(dist, nprows=nprows, col_dist=col_dist)

      ALLOCATE (row_dist(nprows), row_blk_size(nprows))
      row_blk_size(:) = nrow
      DO i = 0, nprows-1
         row_dist(i+1) = i
      END DO
      CALL dbcsr_distribution_new(dist_row_vec, template=dist, row_dist=row_dist, col_dist=col_dist)

      CALL dbcsr_create(dbcsr_vec, "D", dist_row_vec, dbcsr_type_no_symmetry, &
                        row_blk_size=row_blk_size, col_blk_size=col_blk_size, &
                        nze=0, data_type=data_type)
      CALL dbcsr_reserve_all_blocks(dbcsr_vec)

      CALL dbcsr_distribution_release(dist_row_vec)
      DEALLOCATE (row_dist, row_blk_Size)
      CALL timestop(handle)

   END SUBROUTINE create_replicated_row_vec_from_matrix

! **************************************************************************************************
!> \brief given a column vector, prepare the fast_vec_access container
!> \param vec ...
!> \param fast_vec_access ...
! **************************************************************************************************
   SUBROUTINE create_fast_col_vec_access(vec, fast_vec_access)
      TYPE(dbcsr_type)                                   :: vec
      TYPE(fast_vec_access_type)                         :: fast_vec_access

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_fast_col_vec_access'

      INTEGER                                            :: handle, data_type

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(vec, data_type=data_type)

      SELECT CASE (data_type)
      CASE (dbcsr_type_real_8)
         CALL create_fast_col_vec_access_d(vec, fast_vec_access)
      CASE (dbcsr_type_real_4)
         CALL create_fast_col_vec_access_s(vec, fast_vec_access)
      CASE (dbcsr_type_complex_8)
         CALL create_fast_col_vec_access_z(vec, fast_vec_access)
      CASE (dbcsr_type_complex_4)
         CALL create_fast_col_vec_access_c(vec, fast_vec_access)
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE create_fast_col_vec_access

! **************************************************************************************************
!> \brief given a row vector, prepare the fast_vec_access_container
!> \param vec ...
!> \param fast_vec_access ...
! **************************************************************************************************
   SUBROUTINE create_fast_row_vec_access(vec, fast_vec_access)
      TYPE(dbcsr_type)                                   :: vec
      TYPE(fast_vec_access_type)                         :: fast_vec_access

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_fast_row_vec_access'

      INTEGER                                            :: handle, data_type

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(vec, data_type=data_type)

      SELECT CASE (data_type)
      CASE (dbcsr_type_real_8)
         CALL create_fast_row_vec_access_d(vec, fast_vec_access)
      CASE (dbcsr_type_real_4)
         CALL create_fast_row_vec_access_s(vec, fast_vec_access)
      CASE (dbcsr_type_complex_8)
         CALL create_fast_row_vec_access_c(vec, fast_vec_access)
      CASE (dbcsr_type_complex_4)
         CALL create_fast_row_vec_access_z(vec, fast_vec_access)
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE create_fast_row_vec_access

! **************************************************************************************************
!> \brief release all memory associated with the fast_vec_access type
!> \param fast_vec_access ...
! **************************************************************************************************
   SUBROUTINE release_fast_vec_access(fast_vec_access)
      TYPE(fast_vec_access_type)                         :: fast_vec_access

      CHARACTER(LEN=*), PARAMETER :: routineN = 'release_fast_vec_access'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL hash_table_release(fast_vec_access%hash_table)
      IF (ALLOCATED(fast_vec_access%blk_map_d)) DEALLOCATE (fast_vec_access%blk_map_d)
      IF (ALLOCATED(fast_vec_access%blk_map_s)) DEALLOCATE (fast_vec_access%blk_map_s)
      IF (ALLOCATED(fast_vec_access%blk_map_c)) DEALLOCATE (fast_vec_access%blk_map_c)
      IF (ALLOCATED(fast_vec_access%blk_map_z)) DEALLOCATE (fast_vec_access%blk_map_z)

      CALL timestop(handle)

   END SUBROUTINE release_fast_vec_access

#include "hash_table.f90"

#:set instances = [ ('s', 'REAL(kind=real_4)',    '0.0_real_4'), &
                    ('d', 'REAL(kind=real_8)',    '0.0_real_8'), &
                    ('c', 'COMPLEX(kind=real_4)', 'CMPLX(0.0, 0.0, real_4)'), &
                    ('z', 'COMPLEX(kind=real_8)', 'CMPLX(0.0, 0.0, real_8)') ]

#:for nametype, type, zero in instances

! **************************************************************************************************
!> \brief the real driver routine for the multiply, not all symmetries implemented yet
!> \param matrix ...
!> \param vec_in ...
!> \param vec_out ...
!> \param alpha ...
!> \param beta ...
!> \param work_row ...
!> \param work_col ...
! **************************************************************************************************
  SUBROUTINE dbcsr_matrix_colvec_multiply_${nametype}$(matrix, vec_in, vec_out, alpha, beta, work_row, work_col)
    TYPE(dbcsr_type)                          :: matrix, vec_in, vec_out
    ${type}$                                  :: alpha, beta
    TYPE(dbcsr_type)                          :: work_row, work_col

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_matrix_colvec_multiply_low', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: matrix_type

    CALL dbcsr_get_info(matrix, matrix_type=matrix_type)

    SELECT CASE(matrix_type)
    CASE(dbcsr_type_no_symmetry)
       CALL dbcsr_matrix_vector_mult_${nametype}$(matrix, vec_in, vec_out, alpha, beta, work_row, work_col)
    CASE(dbcsr_type_symmetric)
       CALL dbcsr_sym_matrix_vector_mult_${nametype}$(matrix, vec_in, vec_out, alpha, beta, work_row, work_col)
    CASE(dbcsr_type_antisymmetric)
        ! Not yet implemented, should mainly be some prefactor magic, but who knows how antisymmetric matrices are stored???
       CPABORT("NYI, antisymmetric matrix not permitted")
    CASE DEFAULT
       CPABORT("Unknown matrix type, ...")
    END SELECT

  END SUBROUTINE dbcsr_matrix_colvec_multiply_${nametype}$

! **************************************************************************************************
!> \brief low level routines for matrix vector multiplies
!> \param matrix ...
!> \param vec_in ...
!> \param vec_out ...
!> \param alpha ...
!> \param beta ...
!> \param work_row ...
!> \param work_col ...
! **************************************************************************************************
  SUBROUTINE dbcsr_matrix_vector_mult_${nametype}$(matrix, vec_in, vec_out, alpha, beta, work_row, work_col)
    TYPE(dbcsr_type)                          :: matrix, vec_in, vec_out
    ${type}$                                  :: alpha, beta
    TYPE(dbcsr_type)                          :: work_row, work_col

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_matrix_vector_mult'

    INTEGER                                  :: col, mypcol, &
                                                myprow, &
                                                ncols, pcol_group, nrows, &
                                                prow_group, row, &
                                                handle, handle1, ithread 
    ${type}$, DIMENSION(:), POINTER          :: data_vec
    ${type}$, DIMENSION(:, :), POINTER       :: data_d, vec_res
    TYPE(dbcsr_distribution_type)            :: dist
    TYPE(dbcsr_iterator_type)                :: iter
    TYPE(fast_vec_access_type)               :: fast_vec_row, fast_vec_col
    INTEGER                                  :: prow, pcol

    CALL timeset(routineN, handle)
    ithread=0

! Collect some data about the parallel environment. We will use them later to move the vector around
    CALL dbcsr_get_info(matrix, distribution=dist)
    CALL dbcsr_distribution_get(dist, prow_group=prow_group, pcol_group=pcol_group, myprow=myprow, mypcol=mypcol)

    CALL create_fast_row_vec_access(work_row, fast_vec_row)
    CALL create_fast_col_vec_access(work_col, fast_vec_col)

! Transfer the correct parts of the input vector to the correct locations so we can do a local multiply
    CALL dbcsr_col_vec_to_rep_row_${nametype}$(vec_in, work_col, work_row, fast_vec_col)

! Set the work vector for the results to 0
    CALL dbcsr_set(work_col, ${zero}$)

! Perform the local multiply. Here we exploit, that we have the blocks replicated on the mpi processes
! It is important to note, that the input and result vector are distributed differently (row wise, col wise respectively)
    CALL timeset(routineN//"_local_mm", handle1)

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(row,col,iter,data_d,ithread,pcol,prow) &
!$OMP          SHARED(matrix,fast_vec_col,fast_vec_row)
    !$ ithread = omp_get_thread_num ()
    CALL dbcsr_iterator_start(iter, matrix, shared=.FALSE.)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_d)
       prow=hash_table_get(fast_vec_col%hash_table,row)
       IF(fast_vec_col%blk_map_${nametype}$(prow)%assigned_thread .NE. ithread ) CYCLE
       pcol=hash_table_get(fast_vec_row%hash_table,col)
       fast_vec_col%blk_map_${nametype}$(prow)%ptr=fast_vec_col%blk_map_${nametype}$(prow)%ptr+&
            MATMUL(data_d,TRANSPOSE(fast_vec_row%blk_map_${nametype}$(pcol)%ptr))
    END DO
    CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

    CALL timestop(handle1)

! sum all the data onto the first processor col where the original vector is stored
    data_vec => dbcsr_get_data_p(work_col, select_data_type=${zero}$)
    CALL dbcsr_get_info(matrix=work_col, nfullrows_local=nrows, nfullcols_local=ncols)
    CALL mp_sum(data_vec(1:nrows*ncols), prow_group)

! Local copy on the first mpi col (as this is the localtion of the vec_res blocks) of the result vector
! from the replicated to the original vector. Let's play it safe and use the iterator
    CALL dbcsr_iterator_start(iter, vec_out)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_res)
       prow=hash_table_get(fast_vec_col%hash_table,row)
       IF(ASSOCIATED(fast_vec_col%blk_map_${nametype}$(prow)%ptr)) THEN
          vec_res(:, :)= beta*vec_res(:, :)+alpha*fast_vec_col%blk_map_${nametype}$(prow)%ptr(:,:)
       ELSE
          vec_res(:, :)= beta*vec_res(:, :)
       END IF
    END DO
    CALL dbcsr_iterator_stop(iter)

    CALL release_fast_vec_access(fast_vec_row)
    CALL release_fast_vec_access(fast_vec_col)

    CALL timestop(handle)

  END SUBROUTINE dbcsr_matrix_vector_mult_${nametype}$

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param vec_in ...
!> \param vec_out ...
!> \param alpha ...
!> \param beta ...
!> \param work_row ...
!> \param work_col ...
!> \param skip_diag ...
! **************************************************************************************************
  SUBROUTINE dbcsr_matrixT_vector_mult_${nametype}$(matrix, vec_in, vec_out, alpha, beta, work_row, work_col, skip_diag)
    TYPE(dbcsr_type)                          :: matrix, vec_in, vec_out
    ${type}$                                  :: alpha, beta
    TYPE(dbcsr_type)                          :: work_row, work_col
    LOGICAL                                   :: skip_diag

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_matrixT_vector_mult'

    INTEGER                                  :: col, col_size, mypcol, &
                                                myprow, &
                                                ncols, pcol_group, nrows, &
                                                prow_group, row, row_size, &
                                                handle, handle1, ithread
    ${type}$, DIMENSION(:), POINTER          :: data_vec
    ${type}$, DIMENSION(:, :), POINTER       :: data_d, vec_bl, vec_res
    TYPE(dbcsr_distribution_type)            :: dist
    TYPE(dbcsr_iterator_type)                :: iter

    TYPE(fast_vec_access_type)               :: fast_vec_row, fast_vec_col
    INTEGER                                  :: prow, pcol

    CALL timeset(routineN, handle)
    ithread=0

! Collect some data about the parallel environment. We will use them later to move the vector around
    CALL dbcsr_get_info(matrix, distribution=dist)
    CALL dbcsr_distribution_get(dist, prow_group=prow_group, pcol_group=pcol_group, myprow=myprow, mypcol=mypcol)

    CALL create_fast_row_vec_access(work_row, fast_vec_row)
    CALL create_fast_col_vec_access(work_col, fast_vec_col)

! Set the work vector for the results to 0
    CALL dbcsr_set(work_row, ${zero}$)

! Transfer the correct parts of the input vector to the replicated vector on proc_col 0
    CALL dbcsr_iterator_start(iter, vec_in)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_bl, row_size=row_size, col_size=col_size)
       prow=hash_table_get(fast_vec_col%hash_table,row)
       fast_vec_col%blk_map_${nametype}$(prow)%ptr(1:row_size, 1:col_size)= vec_bl(1:row_size, 1:col_size)
    END DO
    CALL dbcsr_iterator_stop(iter)
! Replicate the data on all processore in the row
    data_vec => dbcsr_get_data_p(work_col, select_data_type=${zero}$)
    CALL mp_bcast(data_vec, 0, prow_group)

! Perform the local multiply. Here it is obvious why the vectors are replicated on the mpi rows and cols
    CALL timeset(routineN//"local_mm", handle1)
    CALL dbcsr_get_info(matrix=work_col, nfullcols_local=ncols)
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(row,col,iter,data_d,row_size,col_size,ithread,prow,pcol) &
!$OMP          SHARED(matrix,fast_vec_row,fast_vec_col,skip_diag,ncols)
    !$ ithread = omp_get_thread_num ()
    CALL dbcsr_iterator_start(iter, matrix, shared=.FALSE.)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_d, row_size=row_size, col_size=col_size)
       IF(skip_diag.AND.col==row)CYCLE
       prow=hash_table_get(fast_vec_col%hash_table,row)
       pcol=hash_table_get(fast_vec_row%hash_table,col)
       IF ( ASSOCIATED(fast_vec_row%blk_map_${nametype}$(pcol)%ptr) .AND. &
            ASSOCIATED(fast_vec_col%blk_map_${nametype}$(prow)%ptr) )THEN
          IF(fast_vec_row%blk_map_${nametype}$(pcol)%assigned_thread .NE. ithread ) CYCLE
          fast_vec_row%blk_map_${nametype}$(pcol)%ptr=fast_vec_row%blk_map_${nametype}$(pcol)%ptr+&
               MATMUL(TRANSPOSE(fast_vec_col%blk_map_${nametype}$(prow)%ptr),data_d)
       ELSE
          prow=hash_table_get(fast_vec_row%hash_table,row)
          pcol=hash_table_get(fast_vec_col%hash_table,col)
          IF(fast_vec_row%blk_map_${nametype}$(prow)%assigned_thread .NE. ithread ) CYCLE
          fast_vec_row%blk_map_${nametype}$(prow)%ptr=fast_vec_row%blk_map_${nametype}$(prow)%ptr+&
             MATMUL(TRANSPOSE(fast_vec_col%blk_map_${nametype}$(pcol)%ptr),TRANSPOSE(data_d))
       END IF
    END DO
    CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

    CALL timestop(handle1)

! sum all the data within a processor column to obtain the replicated result
    data_vec => dbcsr_get_data_p(work_row, select_data_type=${zero}$)
! we use the replicated vector but the final answer is only summed to proc_col 0 for efficiency
    CALL dbcsr_get_info(matrix=work_row, nfullrows_local=nrows, nfullcols_local=ncols)
    CALL mp_sum(data_vec(1:nrows*ncols), pcol_group)

! Convert the result to a column wise distribution
    CALL dbcsr_rep_row_to_rep_col_vec_${nametype}$(work_col, work_row, fast_vec_row)
    
! Create_the final vector by summing it to the result vector which lives on proc_col 0
    CALL dbcsr_iterator_start(iter, vec_out)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_res, row_size=row_size)
       prow=hash_table_get(fast_vec_col%hash_table,row)
       IF(ASSOCIATED(fast_vec_col%blk_map_${nametype}$(prow)%ptr)) THEN
          vec_res(:, :)= beta*vec_res(:, :)+alpha*fast_vec_col%blk_map_${nametype}$(prow)%ptr(:,:)
       ELSE
          vec_res(:, :)= beta*vec_res(:, :)
       END IF
    END DO
    CALL dbcsr_iterator_stop(iter)

    CALL timestop(handle)

  END SUBROUTINE dbcsr_matrixT_vector_mult_${nametype}$

! **************************************************************************************************
!> \brief ...
!> \param vec_in ...
!> \param rep_col_vec ...
!> \param rep_row_vec ...
!> \param fast_vec_col ...
! **************************************************************************************************
  SUBROUTINE dbcsr_col_vec_to_rep_row_${nametype}$(vec_in, rep_col_vec, rep_row_vec, fast_vec_col)
    TYPE(dbcsr_type)                          :: vec_in, rep_col_vec, &
                                                rep_row_vec
    TYPE(fast_vec_access_type), INTENT(IN)   :: fast_vec_col

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_col_vec_to_rep_row'

    INTEGER                                  :: col, mypcol, myprow, ncols, &
                                                nrows, pcol_group, &
                                                prow_group, row, handle
    INTEGER, DIMENSION(:), POINTER           :: local_cols, row_dist
    ${type}$, DIMENSION(:), POINTER          :: data_vec, data_vec_rep
    ${type}$, DIMENSION(:, :), POINTER       :: vec_row
    TYPE(dbcsr_distribution_type)            :: dist_in, dist_rep_col
    TYPE(dbcsr_iterator_type)                :: iter

    CALL timeset(routineN, handle)

! get information about the parallel environment
    CALL dbcsr_get_info(vec_in, distribution=dist_in)
    CALL dbcsr_distribution_get(dist_in, &
                                prow_group=prow_group, &
                                pcol_group=pcol_group, &
                                myprow=myprow, &
                                mypcol=mypcol)

! Get the vector which tells us which blocks are local to which processor row in the col vec
    CALL dbcsr_get_info(rep_col_vec, distribution=dist_rep_col)
    CALL dbcsr_distribution_get(dist_rep_col, row_dist=row_dist)

! Copy the local vector to the replicated on the first processor column (this is where vec_in lives)
    CALL dbcsr_get_info(matrix=rep_col_vec, nfullrows_local=nrows, nfullcols_local=ncols)
    data_vec_rep => dbcsr_get_data_p(rep_col_vec, select_data_type=${zero}$)
    data_vec => dbcsr_get_data_p(vec_in, select_data_type=${zero}$)
    IF(mypcol==0)data_vec_rep(1:nrows*ncols)=data_vec(1:nrows*ncols)
! Replicate the data along the row
    CALL mp_bcast(data_vec_rep(1:nrows*ncols), 0, prow_group)

! Here it gets a bit tricky as we are dealing with two different parallel layouts:
! The rep_col_vec contains all blocks local to the row distribution of the vector. 
! The rep_row_vec only needs the fraction which is local to the col distribution.
! However in most cases this won't the complete set of block which can be obtained from col_vector p_row i
! Anyway, as the blocks don't repeat in the col_vec, a different fraction of the row vec will be available
! on every replica in the processor column, by summing along the column we end up with the complete vector everywhere
! Hope this clarifies the idea
    CALL dbcsr_set(rep_row_vec, ${zero}$)
    CALL dbcsr_get_info(matrix=rep_row_vec, nfullrows_local=nrows, local_cols=local_cols, nfullcols_local=ncols)
    CALL dbcsr_iterator_start(iter, rep_row_vec)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_row)
       IF(row_dist(col)==myprow)THEN
          vec_row=TRANSPOSE(fast_vec_col%blk_map_${nametype}$(hash_table_get(fast_vec_col%hash_table,col))%ptr)
       END IF
    END DO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_get_info(matrix=rep_row_vec, nfullrows_local=nrows, nfullcols_local=ncols)
    data_vec_rep => dbcsr_get_data_p(rep_row_vec, select_data_type=${zero}$)
    CALL mp_sum(data_vec_rep(1:ncols*nrows), pcol_group)

    CALL timestop(handle)

  END SUBROUTINE dbcsr_col_vec_to_rep_row_${nametype}$

! **************************************************************************************************
!> \brief ...
!> \param rep_col_vec ...
!> \param rep_row_vec ...
!> \param fast_vec_row ...
!> \param fast_vec_col_add ...
! **************************************************************************************************
  SUBROUTINE dbcsr_rep_row_to_rep_col_vec_${nametype}$(rep_col_vec, rep_row_vec, fast_vec_row, fast_vec_col_add)
    TYPE(dbcsr_type)                          :: rep_col_vec, rep_row_vec
    TYPE(fast_vec_access_type), OPTIONAL     :: fast_vec_col_add
    TYPE(fast_vec_access_type)               :: fast_vec_row

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_rep_row_to_rep_col_vec'

    INTEGER                                  :: col, mypcol, myprow, ncols, &
                                                nrows, pcol_group, &
                                                prow_group, row, handle
    INTEGER, DIMENSION(:), POINTER           :: col_dist
    ${type}$, DIMENSION(:), POINTER          :: data_vec_rep
    ${type}$, DIMENSION(:, :), POINTER       :: vec_col
    TYPE(dbcsr_distribution_type)            :: dist_rep_row, dist_rep_col
    TYPE(dbcsr_iterator_type)                :: iter

    CALL timeset(routineN, handle)

! get information about the parallel environment
    CALL dbcsr_get_info(matrix=rep_col_vec, distribution=dist_rep_col)
    CALL dbcsr_distribution_get(dist_rep_col, &
                                prow_group=prow_group, &
                                pcol_group=pcol_group, &
                                myprow=myprow, &
                                mypcol=mypcol)

! Get the vector which tells us which blocks are local to which processor col in the row vec
    CALL dbcsr_get_info(matrix=rep_row_vec, distribution=dist_rep_row)
    CALL dbcsr_distribution_get(dist_rep_row, col_dist=col_dist)


! The same trick as described above with opposite direction
    CALL dbcsr_set(rep_col_vec, ${zero}$)
    CALL dbcsr_iterator_start(iter, rep_col_vec)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_col)
       IF(col_dist(row)==mypcol)THEN
          vec_col=TRANSPOSE(fast_vec_row%blk_map_${nametype}$(hash_table_get(fast_vec_row%hash_table,row))%ptr)
       END IF
       ! this one is special and allows to add the elements of a not yet summed replicated 
       ! column vector as it appears in M*V(row_rep) as result. Save an mp_sum in the symmetric case
       IF(PRESENT(fast_vec_col_add))vec_col=vec_col+&
            fast_vec_col_add%blk_map_${nametype}$(hash_table_get(fast_vec_col_add%hash_table,row))%ptr(:,:)
    END DO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_get_info(matrix=rep_col_vec, nfullrows_local=nrows, nfullcols_local=ncols)
    data_vec_rep => dbcsr_get_data_p(rep_col_vec, select_data_type=${zero}$)
    CALL mp_sum(data_vec_rep(1:nrows*ncols), prow_group)

    CALL timestop(handle)

  END SUBROUTINE dbcsr_rep_row_to_rep_col_vec_${nametype}$


! **************************************************************************************************
!> \brief given a column vector, prepare the fast_vec_access container
!> \param vec ...
!> \param fast_vec_access ...
! **************************************************************************************************
  SUBROUTINE create_fast_col_vec_access_${nametype}$(vec, fast_vec_access)
    TYPE(dbcsr_type)                          :: vec
    TYPE(fast_vec_access_type)               :: fast_vec_access

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_fast_col_vec_access_${nametype}$'

    INTEGER                                  :: handle, nblk_local
    INTEGER                                  :: col, row, iblock, nthreads
    ${type}$, DIMENSION(:, :), POINTER       :: vec_bl
    TYPE(dbcsr_iterator_type)                :: iter

    CALL timeset(routineN, handle)

    ! figure out the number of threads
    nthreads = 1
!$OMP PARALLEL DEFAULT(NONE) SHARED(nthreads)
!$OMP MASTER
    !$ nthreads = OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

    CALL dbcsr_get_info(matrix=vec, nblkrows_local=nblk_local)
    ! 4 times makes sure the table is big enough to limit collisions.
    CALL hash_table_create(fast_vec_access%hash_table,4*nblk_local)
    ! include zero for effective dealing with values not in the hash table (will return 0)
    ALLOCATE(fast_vec_access%blk_map_${nametype}$(0:nblk_local))

    CALL dbcsr_get_info(matrix=vec, nblkcols_local=col)
    IF (col.GT.1) CPABORT("BUG")

    ! go through the blocks of the vector
    iblock=0
    CALL dbcsr_iterator_start(iter, vec)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_bl)
       iblock=iblock+1
       CALL hash_table_add(fast_vec_access%hash_table,row,iblock)
       fast_vec_access%blk_map_${nametype}$(iblock)%ptr=>vec_bl
       fast_vec_access%blk_map_${nametype}$(iblock)%assigned_thread=MOD(iblock,nthreads)
    END DO
    CALL dbcsr_iterator_stop(iter)

    CALL timestop(handle)

  END SUBROUTINE create_fast_col_vec_access_${nametype}$

! **************************************************************************************************
!> \brief given a row vector, prepare the fast_vec_access_container
!> \param vec ...
!> \param fast_vec_access ...
! **************************************************************************************************
  SUBROUTINE create_fast_row_vec_access_${nametype}$(vec, fast_vec_access)
    TYPE(dbcsr_type)                          :: vec
    TYPE(fast_vec_access_type)                :: fast_vec_access

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_fast_row_vec_access_${nametype}$'

    INTEGER                                  :: handle, nblk_local
    INTEGER                                  :: col, row, iblock, nthreads
    ${type}$, DIMENSION(:, :), POINTER       :: vec_bl
    TYPE(dbcsr_iterator_type)                :: iter

    CALL timeset(routineN, handle)

    ! figure out the number of threads
    nthreads = 1
!$OMP PARALLEL DEFAULT(NONE) SHARED(nthreads)
!$OMP MASTER
    !$ nthreads = OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

    CALL dbcsr_get_info(matrix=vec, nblkcols_local=nblk_local)
    ! 4 times makes sure the table is big enough to limit collisions.
    CALL hash_table_create(fast_vec_access%hash_table,4*nblk_local)
    ! include zero for effective dealing with values not in the hash table (will return 0)
    ALLOCATE(fast_vec_access%blk_map_${nametype}$(0:nblk_local))

    ! sanity check
    CALL dbcsr_get_info(matrix=vec, nblkrows_local=row)
    IF (row.GT.1) CPABORT("BUG")

    ! go through the blocks of the vector
    iblock=0
    CALL dbcsr_iterator_start(iter, vec)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_bl)
       iblock=iblock+1
       CALL hash_table_add(fast_vec_access%hash_table,col,iblock)
       fast_vec_access%blk_map_${nametype}$(iblock)%ptr=>vec_bl
       fast_vec_access%blk_map_${nametype}$(iblock)%assigned_thread=MOD(iblock,nthreads)
    END DO
    CALL dbcsr_iterator_stop(iter)

    CALL timestop(handle)

  END SUBROUTINE create_fast_row_vec_access_${nametype}$

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param vec_in ...
!> \param vec_out ...
!> \param alpha ...
!> \param beta ...
!> \param work_row ...
!> \param work_col ...
! **************************************************************************************************
  SUBROUTINE dbcsr_sym_matrix_vector_mult_${nametype}$(matrix, vec_in, vec_out, alpha, beta, work_row, work_col)
    TYPE(dbcsr_type)                          :: matrix, vec_in, vec_out
    ${type}$                                  :: alpha, beta
    TYPE(dbcsr_type)                          :: work_row, work_col

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_sym_m_v_mult'

    INTEGER                                  :: col, mypcol, &
                                                myprow, &
                                                pcol_group, nrows, ncols,&
                                                prow_group, row, &
                                                handle, handle1, ithread, vec_dim
    ${type}$, DIMENSION(:), POINTER          :: data_vec
    ${type}$, DIMENSION(:, :), POINTER       :: data_d, vec_res
    TYPE(dbcsr_distribution_type)            :: dist
    TYPE(dbcsr_iterator_type)                :: iter
    TYPE(dbcsr_type)                         :: result_row, result_col

    TYPE(fast_vec_access_type)               :: fast_vec_row, fast_vec_col, res_fast_vec_row, res_fast_vec_col
    INTEGER                                  :: prow, pcol, rprow, rpcol

    CALL timeset(routineN, handle)
    ithread=0
! We need some work matrices as we try to exploit operations on the replicated vectors which are duplicated otherwise
    CALL dbcsr_get_info(matrix=vec_in,nfullcols_total=vec_dim)
! This is a performance hack as the new creation of a replicated vector is a fair bit more expensive
    CALL dbcsr_set(work_col, ${zero}$)
    CALL dbcsr_copy(result_col, work_col)
    CALL dbcsr_set(work_row, ${zero}$)
    CALL dbcsr_copy(result_row, work_row)

! Collect some data about the parallel environment. We will use them later to move the vector around
    CALL dbcsr_get_info(matrix=matrix, distribution=dist)
    CALL dbcsr_distribution_get(dist, prow_group=prow_group, pcol_group=pcol_group, myprow=myprow, mypcol=mypcol)

    CALL create_fast_row_vec_access(work_row, fast_vec_row)
    CALL create_fast_col_vec_access(work_col, fast_vec_col)
    CALL create_fast_row_vec_access(result_row, res_fast_vec_row)
    CALL create_fast_col_vec_access(result_col, res_fast_vec_col)

! Transfer the correct parts of the input vector to the correct locations so we can do a local multiply
    CALL dbcsr_col_vec_to_rep_row_${nametype}$(vec_in, work_col, work_row, fast_vec_col)

! Probably I should rename the routine above as it delivers both the replicated row and column vector

! Perform the local multiply. Here we exploit, that we have the blocks replicated on the mpi processes
! It is important to note, that the input and result vector are distributed differently (row wise, col wise respectively)
    CALL timeset(routineN//"_local_mm", handle1)

!------ perform the multiplication, we have to take car to take the correct blocks ----------

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(row,col,iter,data_d,ithread,pcol,prow,rpcol,rprow) &
!$OMP          SHARED(matrix,fast_vec_row,res_fast_vec_col,res_fast_vec_row,fast_vec_col)
    !$ ithread = omp_get_thread_num ()
    CALL dbcsr_iterator_start(iter, matrix, shared=.FALSE.)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_d)
       pcol=hash_table_get(fast_vec_row%hash_table,col)
       rprow=hash_table_get(res_fast_vec_col%hash_table,row)
       IF(ASSOCIATED(fast_vec_row%blk_map_${nametype}$(pcol)%ptr) .AND.&
          ASSOCIATED(res_fast_vec_col%blk_map_${nametype}$(rprow)%ptr))THEN
          IF(res_fast_vec_col%blk_map_${nametype}$(rprow)%assigned_thread .EQ. ithread ) THEN
             res_fast_vec_col%blk_map_${nametype}$(rprow)%ptr=res_fast_vec_col%blk_map_${nametype}$(rprow)%ptr+&
               MATMUL(data_d,TRANSPOSE(fast_vec_row%blk_map_${nametype}$(pcol)%ptr))
          END IF
          prow=hash_table_get(fast_vec_col%hash_table,row)
          rpcol=hash_table_get(res_fast_vec_row%hash_table,col)
          IF(res_fast_vec_row%blk_map_${nametype}$(rpcol)%assigned_thread .EQ. ithread .AND. row .NE. col) THEN
             res_fast_vec_row%blk_map_${nametype}$(rpcol)%ptr=res_fast_vec_row%blk_map_${nametype}$(rpcol)%ptr+&
                MATMUL(TRANSPOSE(fast_vec_col%blk_map_${nametype}$(prow)%ptr),data_d)
          END IF
       ELSE
          rpcol=hash_table_get(res_fast_vec_col%hash_table,col)
          prow=hash_table_get(fast_vec_row%hash_table,row)
          IF(res_fast_vec_col%blk_map_${nametype}$(rpcol)%assigned_thread .EQ. ithread ) THEN
             res_fast_vec_col%blk_map_${nametype}$(rpcol)%ptr=res_fast_vec_col%blk_map_${nametype}$(rpcol)%ptr+&
                TRANSPOSE(MATMUL(fast_vec_row%blk_map_${nametype}$(prow)%ptr,data_d))
          END IF
          rprow=hash_table_get(res_fast_vec_row%hash_table,row)
          pcol=hash_table_get(fast_vec_col%hash_table,col)
          IF(res_fast_vec_row%blk_map_${nametype}$(rprow)%assigned_thread .EQ. ithread  .AND. row .NE. col ) THEN
             res_fast_vec_row%blk_map_${nametype}$(rprow)%ptr=res_fast_vec_row%blk_map_${nametype}$(rprow)%ptr+&
                TRANSPOSE(MATMUL(data_d,fast_vec_col%blk_map_${nametype}$(pcol)%ptr))
          END IF
       END IF 
    END DO
    CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

    CALL timestop(handle1)

    ! sum all the data within a processor column to obtain the replicated result from lower
    data_vec => dbcsr_get_data_p(result_row, select_data_type=${zero}$)
    CALL dbcsr_get_info(matrix=result_row, nfullrows_local=nrows, nfullcols_local=ncols)

    CALL mp_sum(data_vec(1:nrows*ncols), pcol_group)
!
!! Convert the results to a column wise distribution, this is a bit involved as the result_row is fully replicated
!! While the result_col still has the partial results in parallel. The routine below takes care of that and saves an
!! mp_sum. Of the res_row vectors are created only taking the approriate element (0 otherwise) while the res_col
!! parallel bits are locally added. The mp_sum magically creates the correct vector
    CALL dbcsr_rep_row_to_rep_col_vec_${nametype}$(work_col, result_row, res_fast_vec_row, res_fast_vec_col)

!    ! Create_the final vector by summing it to the result vector which lives on proc_col 0 lower
    CALL dbcsr_iterator_start(iter, vec_out)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, vec_res)
       prow=hash_table_get(fast_vec_col%hash_table,row)
       IF(ASSOCIATED(fast_vec_col%blk_map_${nametype}$(prow)%ptr))THEN
          vec_res(:, :)= beta*vec_res(:, :)+alpha*(fast_vec_col%blk_map_${nametype}$(prow)%ptr(:, :)) 
       ELSE
          vec_res(:, :)= beta*vec_res(:, :)
       END IF
    END DO
    CALL dbcsr_iterator_stop(iter)

    CALL release_fast_vec_access(fast_vec_row)
    CALL release_fast_vec_access(fast_vec_col)
    CALL release_fast_vec_access(res_fast_vec_row)
    CALL release_fast_vec_access(res_fast_vec_col)

    CALL dbcsr_release(result_row); CALL dbcsr_release(result_col)

    CALL timestop(handle)

  END SUBROUTINE dbcsr_sym_matrix_vector_mult_${nametype}$

#:endfor

END MODULE dbcsr_vector

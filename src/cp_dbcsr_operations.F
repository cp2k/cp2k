!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   DBCSR operations in CP2K
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.8
!>
!> <b>Modification history:</b>
!> - Created 2009-05-12
! **************************************************************************************************
MODULE cp_dbcsr_operations
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type,&
                                              get_blacs_info
   USE cp_cfm_types,                    ONLY: cp_cfm_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        convert_sizes_to_offsets, dbcsr_add, dbcsr_col_block_sizes, dbcsr_complete_redistribute, &
        dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, dbcsr_distribution, &
        dbcsr_distribution_col_clusters, dbcsr_distribution_col_dist, dbcsr_distribution_mp, &
        dbcsr_distribution_ncols, dbcsr_distribution_new, dbcsr_distribution_obj, &
        dbcsr_distribution_release, dbcsr_distribution_row_clusters, dbcsr_distribution_row_dist, &
        dbcsr_get_data_type, dbcsr_get_info, dbcsr_get_matrix_type, dbcsr_init, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_mp_hold, dbcsr_mp_new, dbcsr_mp_npcols, &
        dbcsr_mp_nprows, dbcsr_mp_obj, dbcsr_mp_release, dbcsr_multiply, dbcsr_name, dbcsr_norm, &
        dbcsr_norm_frobenius, dbcsr_release, dbcsr_reserve_all_blocks, dbcsr_row_block_sizes, &
        dbcsr_scale, dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_complex_8, &
        dbcsr_type_no_symmetry, dbcsr_type_real_8, dbcsr_type_symmetric, dbcsr_valid_index, &
        dbcsr_verify_matrix
   USE distribution_2d_types,           ONLY: distribution_2d_get,&
                                              distribution_2d_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: gcd,&
                                              lcm

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
#include "base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_operations'
   LOGICAL, PARAMETER :: debug_mod = .FALSE.

   INTEGER, SAVE, PUBLIC :: max_elements_per_block = 32

   PUBLIC :: dbcsr_multiply_local

   ! CP2K API emulation
   PUBLIC :: copy_fm_to_dbcsr, copy_dbcsr_to_fm, &
             copy_dbcsr_to_cfm, copy_cfm_to_dbcsr, &
             cp_dbcsr_sm_fm_multiply, cp_dbcsr_plus_fm_fm_t, &
             copy_dbcsr_to_fm_bc, copy_fm_to_dbcsr_bc, cp_fm_to_dbcsr_row_template, &
             cp_dbcsr_m_by_n_from_template, cp_dbcsr_m_by_n_from_row_template, &
             dbcsr_create_dist_r_unrot

   ! distribution_2d_type compatibility
   PUBLIC :: cp_dbcsr_dist2d_to_dist

   PUBLIC :: dbcsr_copy_columns_hack

   

   PRIVATE

CONTAINS




! **************************************************************************************************
!> \brief   Copy a BLACS matrix to a dbcsr matrix.
!>
!>          real_matrix=beta*real_matrix+alpha*fm
!>          beta defaults to 0, alpha to 1
!> \param[in] fm              full matrix
!> \param[out] matrix         DBCSR matrix
!> \param[in] keep_sparsity   (optional) retains the sparsity of the input
!>                            matrix
!> \date    2009-10-13
!> \par History
!>          2009-10-13 rewritten based on copy_dbcsr_to_fm
!> \author  Urban Borstnik
!> \version 2.0
! **************************************************************************************************
   SUBROUTINE copy_fm_to_dbcsr(fm, matrix, keep_sparsity)
      TYPE(cp_fm_type), POINTER                          :: fm
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_sparsity

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'copy_fm_to_dbcsr'

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: bc_mat

      CALL timeset(routineN, handle)

      CALL copy_fm_to_dbcsr_bc(fm, bc_mat)
      CALL dbcsr_complete_redistribute(bc_mat, matrix, keep_sparsity=keep_sparsity)
      CALL dbcsr_release(bc_mat)

      CALL timestop(handle)
   END SUBROUTINE copy_fm_to_dbcsr


! **************************************************************************************************
!> \brief   Copy a BLACS matrix to a dbcsr matrix with a special block-cyclic distribution,
!>           which requires no complete redistribution.
!> \param fm ...
!> \param bc_mat ...
! **************************************************************************************************
   SUBROUTINE copy_fm_to_dbcsr_bc(fm, bc_mat)
      TYPE(cp_fm_type), POINTER                          :: fm
      TYPE(dbcsr_type)                                   :: bc_mat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_fm_to_dbcsr_bc'

      INTEGER                                            :: col, group, handle, mypcol, myproc, &
                                                            myprow, ncol_block, ncol_global, &
                                                            nrow_block, nrow_global, numproc, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_col, first_row, last_col, last_row
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dbcsr_block, fm_block
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(dbcsr_distribution_obj)                       :: bc_dist
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_mp_obj)                                 :: mp_env

      CALL timeset(routineN, handle)

      IF (fm%use_sp) CPABORT("copy_fm_to_dbcsr_bc: single precision not supported")

      ! Create processor grid
      context => fm%matrix_struct%context
      pgrid => context%blacs2mpi
      group = fm%matrix_struct%para_env%group
      myproc = fm%matrix_struct%para_env%mepos
      numproc = fm%matrix_struct%para_env%num_pe
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      CALL dbcsr_mp_new(mp_env, pgrid, group, myproc, numproc, myprow, mypcol)

      ! Create a block-cyclic distribution compatible with the FM matrix.
      nrow_block = fm%matrix_struct%nrow_block
      ncol_block = fm%matrix_struct%ncol_block
      nrow_global = fm%matrix_struct%nrow_global
      ncol_global = fm%matrix_struct%ncol_global
      NULLIFY (col_blk_size, row_blk_size)
      CALL dbcsr_create_dist_block_cyclic(bc_dist, &
                                          nrows=nrow_global, ncolumns=ncol_global, & ! Actual full matrix size
                                          nrow_block=nrow_block, ncol_block=ncol_block, & ! BLACS parameters
                                          mp_env=mp_env, &
                                          row_blk_sizes=row_blk_size, col_blk_sizes=col_blk_size) ! block-cyclic row/col sizes
      CALL dbcsr_mp_release(mp_env)

      ! Create the block-cyclic DBCSR matrix
      CALL dbcsr_init(bc_mat)
      CALL dbcsr_create(bc_mat, "Block-cyclic ", bc_dist, &
                           dbcsr_type_no_symmetry, row_blk_size, col_blk_size, nze=0, &
                           reuse_arrays=.TRUE., data_type=dbcsr_type_real_8)
      CALL dbcsr_distribution_release(bc_dist)

      ! allocate all blocks
      CALL dbcsr_reserve_all_blocks(bc_mat)

      CALL calculate_fm_block_ranges(bc_mat, first_row, last_row, first_col, last_col)

      ! Copy the FM data to the block-cyclic DBCSR matrix.  This step
      ! could be skipped with appropriate DBCSR index manipulation.
      fm_block => fm%local_data
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(iter, row, col, dbcsr_block) &
!$OMP SHARED(bc_mat, last_row, first_row, last_col, first_col, fm_block)
       CALL dbcsr_iterator_start(iter, bc_mat)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          CALL dbcsr_iterator_next_block(iter, row, col, dbcsr_block)
          dbcsr_block(:, :) = fm_block(first_row(row):last_row(row), first_col(col):last_col(col))
       ENDDO
       CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)
   END SUBROUTINE copy_fm_to_dbcsr_bc


! **************************************************************************************************
!> \brief Copy a DBCSR matrix to a BLACS matrix
!> \param[in] matrix          DBCSR matrix
!> \param[out] fm             full matrix
! **************************************************************************************************
   SUBROUTINE copy_dbcsr_to_fm(matrix, fm)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(cp_fm_type), POINTER                          :: fm

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'copy_dbcsr_to_fm'

      INTEGER                                            :: handle, ncol_block, nfullcols_total, &
                                                            nfullrows_total, nrow_block
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      TYPE(dbcsr_distribution_obj)                       :: bc_dist
      TYPE(dbcsr_type)                                   :: bc_mat, matrix_nosym

      CALL timeset(routineN, handle)

      ! check compatibility
      CALL dbcsr_get_info(matrix, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      CPASSERT(fm%matrix_struct%nrow_global==nfullrows_total)
      CPASSERT(fm%matrix_struct%ncol_global==nfullcols_total)

      ! info about the full matrix
      nrow_block = fm%matrix_struct%nrow_block
      ncol_block = fm%matrix_struct%ncol_block

      ! Convert DBCSR to a block-cyclic
      NULLIFY (col_blk_size, row_blk_size)
      CALL dbcsr_create_dist_block_cyclic(bc_dist, &
                                          nfullrows_total, nfullcols_total, &
                                          nrow_block, ncol_block, &
                                          dbcsr_distribution_mp(dbcsr_distribution(matrix)), &
                                          row_blk_size, col_blk_size)

      CALL dbcsr_init(bc_mat)
      CALL dbcsr_create(bc_mat, "Block-cyclic"//dbcsr_name(matrix), bc_dist, &
                           dbcsr_type_no_symmetry, row_blk_size, col_blk_size, &
                           nze=0, data_type=dbcsr_get_data_type(matrix), &
                           reuse_arrays=.TRUE.)
      CALL dbcsr_distribution_release(bc_dist)

      CALL dbcsr_init(matrix_nosym)
      CALL dbcsr_create(matrix_nosym, template=matrix, matrix_type="N")
      CALL dbcsr_desymmetrize(matrix, matrix_nosym)
      CALL dbcsr_complete_redistribute(matrix_nosym, bc_mat)
      CALL dbcsr_release(matrix_nosym)

      CALL copy_dbcsr_to_fm_bc(bc_mat, fm)

      CALL dbcsr_release(bc_mat)

      CALL timestop(handle)
   END SUBROUTINE copy_dbcsr_to_fm


! **************************************************************************************************
!> \brief Copy a DBCSR_BLACS matrix to a BLACS matrix
!> \param bc_mat DBCSR matrix
!> \param[out] fm             full matrix
! **************************************************************************************************
   SUBROUTINE copy_dbcsr_to_fm_bc(bc_mat, fm)
      TYPE(dbcsr_type), INTENT(IN)                       :: bc_mat
      TYPE(cp_fm_type), POINTER                          :: fm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_dbcsr_to_fm_bc'

      INTEGER                                            :: col, handle, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_col, first_row, last_col, last_row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dbcsr_block, fm_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      IF (fm%use_sp) CPABORT("copy_dbcsr_to_fm_bc: single precision not supported")

      CALL calculate_fm_block_ranges(bc_mat, first_row, last_row, first_col, last_col)

      ! Now copy data to the FM matrix
      fm_block => fm%local_data
      fm_block = REAL(0.0, KIND=dp)
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(iter, row, col, dbcsr_block) &
!$OMP SHARED(bc_mat, last_row, first_row, last_col, first_col, fm_block)
      CALL dbcsr_iterator_start(iter, bc_mat)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, dbcsr_block)
         fm_block(first_row(row):last_row(row), first_col(col):last_col(col)) = dbcsr_block(:, :)
      ENDDO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)
   END SUBROUTINE copy_dbcsr_to_fm_bc


! **************************************************************************************************
!> \brief   Copy a BLACS matrix to a dbcsr matrix.
!>
!>          real_matrix=beta*real_matrix+alpha*fm
!>          beta defaults to 0, alpha to 1
!> \param[in] fm              full matrix
!> \param[out] matrix         DBCSR matrix
!> \param[in] keep_sparsity   (optional) retains the sparsity of the input
!>                            matrix
!> \date    2009-10-13
!> \par History
!>          2009-10-13 rewritten based on copy_dbcsr_to_fm
!> \author  Urban Borstnik
!> \version 2.0
! **************************************************************************************************
   SUBROUTINE copy_cfm_to_dbcsr(fm, matrix, keep_sparsity)
      TYPE(cp_cfm_type), POINTER                         :: fm
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_sparsity

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'copy_cfm_to_dbcsr'

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: bc_mat

      CALL timeset(routineN, handle)

      CALL copy_cfm_to_dbcsr_bc(fm, bc_mat)
      CALL dbcsr_complete_redistribute(bc_mat, matrix, keep_sparsity=keep_sparsity)
      CALL dbcsr_release(bc_mat)

      CALL timestop(handle)
   END SUBROUTINE copy_cfm_to_dbcsr


! **************************************************************************************************
!> \brief   Copy a BLACS matrix to a dbcsr matrix with a special block-cyclic distribution,
!>           which requires no complete redistribution.
!> \param fm ...
!> \param bc_mat ...
! **************************************************************************************************
   SUBROUTINE copy_cfm_to_dbcsr_bc(fm, bc_mat)
      TYPE(cp_cfm_type), POINTER                         :: fm
      TYPE(dbcsr_type)                                   :: bc_mat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_cfm_to_dbcsr_bc'

      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: dbcsr_block, fm_block
      INTEGER                                            :: col, group, handle, mypcol, myproc, &
                                                            myprow, ncol_block, ncol_global, &
                                                            nrow_block, nrow_global, numproc, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_col, first_row, last_col, last_row
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(dbcsr_distribution_obj)                       :: bc_dist
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_mp_obj)                                 :: mp_env

      CALL timeset(routineN, handle)


      ! Create processor grid
      context => fm%matrix_struct%context
      pgrid => context%blacs2mpi
      group = fm%matrix_struct%para_env%group
      myproc = fm%matrix_struct%para_env%mepos
      numproc = fm%matrix_struct%para_env%num_pe
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      CALL dbcsr_mp_new(mp_env, pgrid, group, myproc, numproc, myprow, mypcol)

      ! Create a block-cyclic distribution compatible with the FM matrix.
      nrow_block = fm%matrix_struct%nrow_block
      ncol_block = fm%matrix_struct%ncol_block
      nrow_global = fm%matrix_struct%nrow_global
      ncol_global = fm%matrix_struct%ncol_global
      NULLIFY (col_blk_size, row_blk_size)
      CALL dbcsr_create_dist_block_cyclic(bc_dist, &
                                          nrows=nrow_global, ncolumns=ncol_global, & ! Actual full matrix size
                                          nrow_block=nrow_block, ncol_block=ncol_block, & ! BLACS parameters
                                          mp_env=mp_env, &
                                          row_blk_sizes=row_blk_size, col_blk_sizes=col_blk_size) ! block-cyclic row/col sizes
      CALL dbcsr_mp_release(mp_env)

      ! Create the block-cyclic DBCSR matrix
      CALL dbcsr_init(bc_mat)
      CALL dbcsr_create(bc_mat, "Block-cyclic ", bc_dist, &
                           dbcsr_type_no_symmetry, row_blk_size, col_blk_size, nze=0, &
                           reuse_arrays=.TRUE., data_type=dbcsr_type_complex_8)
      CALL dbcsr_distribution_release(bc_dist)

      ! allocate all blocks
      CALL dbcsr_reserve_all_blocks(bc_mat)

      CALL calculate_fm_block_ranges(bc_mat, first_row, last_row, first_col, last_col)

      ! Copy the FM data to the block-cyclic DBCSR matrix.  This step
      ! could be skipped with appropriate DBCSR index manipulation.
      fm_block => fm%local_data
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(iter, row, col, dbcsr_block) &
!$OMP SHARED(bc_mat, last_row, first_row, last_col, first_col, fm_block)
       CALL dbcsr_iterator_start(iter, bc_mat)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          CALL dbcsr_iterator_next_block(iter, row, col, dbcsr_block)
          dbcsr_block(:, :) = fm_block(first_row(row):last_row(row), first_col(col):last_col(col))
       ENDDO
       CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)
   END SUBROUTINE copy_cfm_to_dbcsr_bc


! **************************************************************************************************
!> \brief Copy a DBCSR matrix to a BLACS matrix
!> \param[in] matrix          DBCSR matrix
!> \param[out] fm             full matrix
! **************************************************************************************************
   SUBROUTINE copy_dbcsr_to_cfm(matrix, fm)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(cp_cfm_type), POINTER                         :: fm

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'copy_dbcsr_to_cfm'

      INTEGER                                            :: handle, ncol_block, nfullcols_total, &
                                                            nfullrows_total, nrow_block
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      TYPE(dbcsr_distribution_obj)                       :: bc_dist
      TYPE(dbcsr_type)                                   :: bc_mat, matrix_nosym

      CALL timeset(routineN, handle)

      ! check compatibility
      CALL dbcsr_get_info(matrix, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      CPASSERT(fm%matrix_struct%nrow_global==nfullrows_total)
      CPASSERT(fm%matrix_struct%ncol_global==nfullcols_total)

      ! info about the full matrix
      nrow_block = fm%matrix_struct%nrow_block
      ncol_block = fm%matrix_struct%ncol_block

      ! Convert DBCSR to a block-cyclic
      NULLIFY (col_blk_size, row_blk_size)
      CALL dbcsr_create_dist_block_cyclic(bc_dist, &
                                          nfullrows_total, nfullcols_total, &
                                          nrow_block, ncol_block, &
                                          dbcsr_distribution_mp(dbcsr_distribution(matrix)), &
                                          row_blk_size, col_blk_size)

      CALL dbcsr_init(bc_mat)
      CALL dbcsr_create(bc_mat, "Block-cyclic"//dbcsr_name(matrix), bc_dist, &
                           dbcsr_type_no_symmetry, row_blk_size, col_blk_size, &
                           nze=0, data_type=dbcsr_get_data_type(matrix), &
                           reuse_arrays=.TRUE.)
      CALL dbcsr_distribution_release(bc_dist)

      CALL dbcsr_init(matrix_nosym)
      CALL dbcsr_create(matrix_nosym, template=matrix, matrix_type="N")
      CALL dbcsr_desymmetrize(matrix, matrix_nosym)
      CALL dbcsr_complete_redistribute(matrix_nosym, bc_mat)
      CALL dbcsr_release(matrix_nosym)

      CALL copy_dbcsr_to_cfm_bc(bc_mat, fm)

      CALL dbcsr_release(bc_mat)

      CALL timestop(handle)
   END SUBROUTINE copy_dbcsr_to_cfm


! **************************************************************************************************
!> \brief Copy a DBCSR_BLACS matrix to a BLACS matrix
!> \param bc_mat DBCSR matrix
!> \param[out] fm             full matrix
! **************************************************************************************************
   SUBROUTINE copy_dbcsr_to_cfm_bc(bc_mat, fm)
      TYPE(dbcsr_type), INTENT(IN)                       :: bc_mat
      TYPE(cp_cfm_type), POINTER                         :: fm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_dbcsr_to_cfm_bc'

      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: dbcsr_block, fm_block
      INTEGER                                            :: col, handle, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_col, first_row, last_col, last_row
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)


      CALL calculate_fm_block_ranges(bc_mat, first_row, last_row, first_col, last_col)

      ! Now copy data to the FM matrix
      fm_block => fm%local_data
      fm_block = CMPLX(0.0, KIND=dp)
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(iter, row, col, dbcsr_block) &
!$OMP SHARED(bc_mat, last_row, first_row, last_col, first_col, fm_block)
      CALL dbcsr_iterator_start(iter, bc_mat)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, dbcsr_block)
         fm_block(first_row(row):last_row(row), first_col(col):last_col(col)) = dbcsr_block(:, :)
      ENDDO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)
   END SUBROUTINE copy_dbcsr_to_cfm_bc


! **************************************************************************************************
!> \brief Helper routine used to copy blocks from DBCSR into FM matrices and vice versa
!> \param bc_mat ...
!> \param first_row ...
!> \param last_row ...
!> \param first_col ...
!> \param last_col ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE calculate_fm_block_ranges(bc_mat, first_row, last_row, first_col, last_col)
      TYPE(dbcsr_type), INTENT(IN)                       :: bc_mat
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_row, last_row, first_col, last_col

      INTEGER                                            :: col, nblkcols_local, nblkcols_total, &
                                                            nblkrows_local, nblkrows_total, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: local_col_sizes, local_row_sizes
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, local_cols, local_rows, &
                                                            row_blk_size

      CALL dbcsr_get_info(bc_mat, &
                             nblkrows_total=nblkrows_total, &
                             nblkcols_total=nblkcols_total, &
                             nblkrows_local=nblkrows_local, &
                             nblkcols_local=nblkcols_local, &
                             local_rows=local_rows, &
                             local_cols=local_cols, &
                             row_blk_size=row_blk_size, &
                             col_blk_size=col_blk_size )

      ! calculate first_row and last_row
      ALLOCATE (local_row_sizes(nblkrows_total))
      local_row_sizes(:) = 0
      IF (nblkrows_local .GE. 1) THEN
         DO row = 1, nblkrows_local
            local_row_sizes(local_rows(row)) = row_blk_size(local_rows(row))
         END DO
      ENDIF
      ALLOCATE (first_row(nblkrows_total), last_row(nblkrows_total))
      CALL convert_sizes_to_offsets(local_row_sizes, first_row, last_row)

      ! calculate first_col and last_col
      ALLOCATE (local_col_sizes(nblkcols_total))
      local_col_sizes(:) = 0
      IF (nblkcols_local .GE. 1) THEN
         DO col = 1, nblkcols_local
            local_col_sizes(local_cols(col)) = col_blk_size(local_cols(col))
         END DO
      ENDIF
      ALLOCATE (first_col(nblkcols_total), last_col(nblkcols_total))
      CALL convert_sizes_to_offsets(local_col_sizes, first_col, last_col)

  END SUBROUTINE calculate_fm_block_ranges


! **************************************************************************************************
!> \brief hack for dbcsr_copy_columns
!> \param matrix_b ...
!> \param matrix_a ...
!> \param ncol ...
!> \param source_start ...
!> \param target_start ...
!> \param para_env ...
!> \param blacs_env ...
!> \author vw
! **************************************************************************************************
   SUBROUTINE dbcsr_copy_columns_hack(matrix_b, matrix_a, &
                                         ncol, source_start, target_start, para_env, blacs_env)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_b
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a
      INTEGER, INTENT(IN)                                :: ncol, source_start, target_start
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_copy_columns_hack', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: nfullcols_total, nfullrows_total
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_a, fm_matrix_b

      NULLIFY (fm_matrix_a, fm_matrix_b, fm_struct)
      CALL dbcsr_get_info(matrix_a, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nfullrows_total, &
                               ncol_global=nfullcols_total, para_env=para_env)
      CALL cp_fm_create(fm_matrix_a, fm_struct, name="fm_matrix_a")
      CALL cp_fm_struct_release(fm_struct)

      CALL dbcsr_get_info(matrix_b, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nfullrows_total, &
                               ncol_global=nfullcols_total, para_env=para_env)
      CALL cp_fm_create(fm_matrix_b, fm_struct, name="fm_matrix_b")
      CALL cp_fm_struct_release(fm_struct)

      CALL copy_dbcsr_to_fm(matrix_a, fm_matrix_a)
      CALL copy_dbcsr_to_fm(matrix_b, fm_matrix_b)

      CALL cp_fm_to_fm(fm_matrix_a, fm_matrix_b, ncol, source_start, target_start)

      CALL copy_fm_to_dbcsr(fm_matrix_b, matrix_b)

      CALL cp_fm_release(fm_matrix_a)
      CALL cp_fm_release(fm_matrix_b)

   END SUBROUTINE dbcsr_copy_columns_hack

! **************************************************************************************************
!> \brief Creates a DBCSR distribution from a distribution_2d
!> \param[in] dist2d          distribution_2d
!> \param[out] dist           DBCSR distribution
!> \param mp_obj ...
!> \par History
!>    move form dbcsr_operation 01.2010
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_dist2d_to_dist(dist2d, dist, mp_obj)
      TYPE(distribution_2d_type), INTENT(IN), TARGET     :: dist2d
      TYPE(dbcsr_distribution_obj), INTENT(OUT)          :: dist
      TYPE(dbcsr_mp_obj), INTENT(IN), OPTIONAL           :: mp_obj

      INTEGER                                            :: mypcol, myproc, myprow, numproc
      INTEGER, DIMENSION(:, :), POINTER                  :: col_dist_data, pgrid, row_dist_data
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_mp_obj)                                 :: mp_env
      TYPE(distribution_2d_type), POINTER                :: dist2d_p

!
!   ---------------------------------------------------------------------------

      dist2d_p => dist2d
      CALL distribution_2d_get(dist2d_p, &
                               row_distribution=row_dist_data, col_distribution=col_dist_data, &
                               blacs_env=blacs_env)
      CALL get_blacs_info(blacs_env, para_env=para_env, &
                          my_process_row=myprow, my_process_column=mypcol, &
                          blacs2mpi=pgrid)
      myproc = para_env%mepos
      numproc = para_env%num_pe
      IF (PRESENT(mp_obj)) THEN
         mp_env = mp_obj
         CALL dbcsr_mp_hold(mp_env)
      ELSE
         CALL dbcsr_mp_new(mp_env, pgrid, para_env%group, myproc, numproc, &
                              myprow, mypcol)
      ENDIF
      CALL dbcsr_distribution_new(dist, mp_env, row_dist_data, col_dist_data)
      CALL dbcsr_mp_release(mp_env)
   END SUBROUTINE cp_dbcsr_dist2d_to_dist

! **************************************************************************************************
!> \brief multiply a dbcsr with a replicated array
!>        c = alpha_scalar * A (dbscr) * b + c
!> \param[in] matrix_a DBSCR matrxx
!> \param[in]  vec_b        vectors b
!> \param[inout] vec_c      vectors c
!> \param[in]  ncol         nbr of columns
!> \param[in]  alpha        alpha
!>
! **************************************************************************************************
   SUBROUTINE dbcsr_multiply_local(matrix_a, vec_b, vec_c, ncol, alpha)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: vec_b
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: vec_c
      INTEGER, INTENT(in), OPTIONAL                      :: ncol
      REAL(dp), INTENT(IN), OPTIONAL                     :: alpha

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_multiply_local', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, col, coloff, my_ncol, row, rowoff, &
                                                            timing_handle
      LOGICAL                                            :: has_symm
      REAL(dp)                                           :: my_alpha, my_alpha2
      REAL(dp), DIMENSION(:, :), POINTER                 :: data_d
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, timing_handle)

      my_alpha = 1.0_dp
      IF (PRESENT(alpha)) my_alpha = alpha

      my_ncol = SIZE(vec_b, 2)
      IF (PRESENT(ncol)) my_ncol = ncol

      my_alpha2 = 0.0_dp
      IF (dbcsr_get_matrix_type(matrix_a) .EQ. dbcsr_type_symmetric) my_alpha2 = my_alpha
      IF (dbcsr_get_matrix_type(matrix_a) .EQ. dbcsr_type_antisymmetric) my_alpha2 = -my_alpha

      has_symm = (dbcsr_get_matrix_type(matrix_a) .EQ. dbcsr_type_symmetric .OR. &
                  dbcsr_get_matrix_type(matrix_a) .EQ. dbcsr_type_antisymmetric)

!$OMP     PARALLEL DEFAULT(NONE) SHARED(matrix_a,vec_b,vec_c,ncol,my_alpha2,my_alpha,my_ncol,has_symm) &
!$OMP              PRIVATE(iter,row,col,data_d,blk,rowoff,coloff)
      CALL dbcsr_iterator_start(iter, matrix_a, read_only=.TRUE., dynamic=.TRUE., dynamic_byrows=.TRUE.)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, data_d, blk, row_offset=rowoff, col_offset=coloff)
         IF (my_ncol .NE. 1) THEN
            CALL dgemm('N', 'N', &
                       SIZE(data_d, 1), my_ncol, SIZE(data_d, 2), &
                       my_alpha, data_d(1, 1), SIZE(data_d, 1), &
                       vec_b(coloff, 1), SIZE(vec_b, 1), &
                       1.0_dp, vec_c(rowoff, 1), SIZE(vec_c, 1))
         ELSE
            CALL dgemv('N', SIZE(data_d, 1), SIZE(data_d, 2), &
                       my_alpha, data_d(1, 1), SIZE(data_d, 1), &
                       vec_b(coloff, 1), 1, &
                       1.0_dp, vec_c(rowoff, 1), 1)
         ENDIF
      ENDDO
      CALL dbcsr_iterator_stop(iter)
!$OMP     END PARALLEL

      ! FIXME ... in the symmetric case, the writes to vec_c depend on the column, not the row. This makes OMP-ing more difficult
      ! needs e.g. a buffer for vec_c and a reduction of that buffer.
      IF (has_symm) THEN
         CALL dbcsr_iterator_start(iter, matrix_a)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row, col, data_d, blk, row_offset=rowoff, col_offset=coloff)
            IF (row .NE. col) THEN
               IF (my_ncol .NE. 1) THEN
                  CALL dgemm('T', 'N', &
                             SIZE(data_d, 2), my_ncol, SIZE(data_d, 1), &
                             my_alpha2, data_d(1, 1), SIZE(data_d, 1), &
                             vec_b(rowoff, 1), SIZE(vec_b, 1), &
                             1.0_dp, vec_c(coloff, 1), SIZE(vec_c, 1))
               ELSE
                  CALL dgemv('T', SIZE(data_d, 1), SIZE(data_d, 2), &
                             my_alpha2, data_d(1, 1), SIZE(data_d, 1), &
                             vec_b(rowoff, 1), 1, &
                             1.0_dp, vec_c(coloff, 1), 1)
               ENDIF
            ENDIF
         ENDDO
         CALL dbcsr_iterator_stop(iter)
      ENDIF

      CALL timestop(timing_handle)
   END SUBROUTINE dbcsr_multiply_local


! **************************************************************************************************
!> \brief multiply a dbcsr with a fm matrix
!>
!> For backwards compatibility with BLAS XGEMM, this routine supports
!> the multiplication of matrices with incompatible dimensions.
!>
!> \param[in]  matrix         DBCSR matrix
!> \param fm_in full matrix
!> \param fm_out full matrix
!> \param[in]  ncol           nbr of columns
!> \param[in]  alpha          alpha
!> \param[in]  beta           beta
!>
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_sm_fm_multiply(matrix, fm_in, fm_out, ncol, alpha, beta)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(cp_fm_type), POINTER                          :: fm_in, fm_out
      INTEGER, INTENT(IN)                                :: ncol
      REAL(dp), INTENT(IN), OPTIONAL                     :: alpha, beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_sm_fm_multiply', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: k_in, k_out, timing_handle, &
                                                            timing_handle_mult
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size_right_in, &
                                                            col_blk_size_right_out
      TYPE(dbcsr_type)                                   :: in, out
      TYPE(dbcsr_distribution_obj)                       :: dist_right_in, product_dist
      REAL(dp)                                           :: my_alpha, my_beta

      CALL timeset(routineN, timing_handle)

      my_alpha = 1.0_dp
      my_beta = 0.0_dp
      IF (PRESENT(alpha)) my_alpha = alpha
      IF (PRESENT(beta)) my_beta = beta

      CALL cp_fm_get_info(fm_out, ncol_global=k_out)

      CALL cp_fm_get_info(fm_in, ncol_global=k_in)
      !write(*,*)routineN//" -----------------------------------"
      !IF (k_in .NE. k_out) &
      !   WRITE(*,'(3(A,I5,1X),2(A,F5.2,1X))')&
      !   routineN//" ncol", ncol,'k_in',k_in,'k_out',k_out,&
      !   'alpha',my_alpha,'beta',my_beta

      IF (ncol .GT. 0 .AND. k_out .GT. 0 .AND. k_in .GT. 0) THEN
         CALL dbcsr_create_dist_r_unrot(dist_right_in, dbcsr_distribution(matrix), k_in, &
                                           col_blk_size_right_in)

         CALL dbcsr_init(in)
         CALL dbcsr_create(in, "D", dist_right_in, dbcsr_type_no_symmetry, &
                              dbcsr_row_block_sizes(matrix), col_blk_size_right_in, &
                              nze=0)

         CALL dbcsr_init(out)
         CALL dbcsr_distribution_new(product_dist, &
                                        dbcsr_distribution_mp(dbcsr_distribution(matrix)), &
                                        dbcsr_distribution_row_dist(dbcsr_distribution(matrix)), &
                                        dbcsr_distribution_col_dist(dist_right_in), &
                                        dbcsr_distribution_row_clusters(dbcsr_distribution(matrix)), &
                                        dbcsr_distribution_col_clusters(dist_right_in))
         ALLOCATE (col_blk_size_right_out(SIZE(col_blk_size_right_in)))
         col_blk_size_right_out = col_blk_size_right_in
         CALL match_col_sizes(col_blk_size_right_out, col_blk_size_right_in, k_out)

         !if (k_in .ne. k_out) then
         !   write(*,*)routineN//" in cs", col_blk_size_right_in
         !   write(*,*)routineN//" out cs", col_blk_size_right_out
         !endif

         CALL dbcsr_create(out, "D", product_dist, dbcsr_type_no_symmetry, &
                              dbcsr_row_block_sizes(matrix), col_blk_size_right_out, &
                              nze=0)

         CALL copy_fm_to_dbcsr(fm_in, in)
         IF (ncol .NE. k_out .OR. my_beta .NE. 0.0_dp) &
            CALL copy_fm_to_dbcsr(fm_out, out)

         CALL timeset(routineN//'_core', timing_handle_mult)
         CALL dbcsr_multiply("N", "N", my_alpha, matrix, in, my_beta, out, &
                                last_column=ncol, match_matrix_sizes=.TRUE.)
         CALL timestop(timing_handle_mult)

         CALL copy_dbcsr_to_fm(out, fm_out)

         CALL dbcsr_release(in)
         CALL dbcsr_release(out)
         DEALLOCATE (col_blk_size_right_in, col_blk_size_right_out)
         CALL dbcsr_distribution_release(dist_right_in)
         CALL dbcsr_distribution_release(product_dist)

      ENDIF

      CALL timestop(timing_handle)

   END SUBROUTINE cp_dbcsr_sm_fm_multiply

! **************************************************************************************************
!> \brief ...
!> \param sizes1 ...
!> \param sizes2 ...
!> \param full_num ...
! **************************************************************************************************
   SUBROUTINE match_col_sizes(sizes1, sizes2, full_num)
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: sizes1
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes2
      INTEGER, INTENT(IN)                                :: full_num

      INTEGER                                            :: left, n1, n2, p, rm, used

      n1 = SIZE(sizes1)
      n2 = SIZE(sizes2)
      IF (n1 .NE. n2) &
         CPABORT("distributions must be equal!")
      sizes1(1:n1) = sizes2(1:n1)
      used = SUM(sizes1(1:n1))
      ! If sizes1 does not cover everything, then we increase the
      ! size of the last block; otherwise we reduce the blocks
      ! (from the end) until it is small enough.
      IF (used .LT. full_num) THEN
         sizes1(n1) = sizes1(n1)+full_num-used
      ELSE
         left = used-full_num
         p = n1
         DO WHILE (left .GT. 0 .AND. p .GT. 0)
            rm = MIN(left, sizes1(p))
            sizes1(p) = sizes1(p)-rm
            left = left-rm
            p = p-1
         ENDDO
      ENDIF
   END SUBROUTINE match_col_sizes

! **************************************************************************************************
!> \brief performs the multiplication sparse_matrix+dense_mat*dens_mat^T
!>        if matrix_g is not explicitly given, matrix_v^T will be used
!>        this can be important to save the necessary redistribute for a
!>        different matrix_g and increase performance.
!> \param sparse_matrix ...
!> \param matrix_v ...
!> \param matrix_g ...
!> \param ncol ...
!> \param alpha ...
!> \param keep_sparsity Determines if the sparsity of sparse_matrix is retained
!>        by default it is TRUE
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_plus_fm_fm_t(sparse_matrix, matrix_v, matrix_g, ncol, &
                                           alpha, keep_sparsity)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: sparse_matrix
      TYPE(cp_fm_type), POINTER                          :: matrix_v
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: matrix_g
      INTEGER, INTENT(IN)                                :: ncol
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: alpha
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_sparsity

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_plus_fm_fm_t', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: k, nao, timing_handle
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size_left, col_dist_cluster, &
                                                            col_dist_left
      LOGICAL                                            :: check_product, my_keep_sparsity
      REAL(KIND=dp)                                      :: my_alpha, norm
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: fm_matrix
      TYPE(dbcsr_distribution_obj)                       :: dist_left
      TYPE(dbcsr_mp_obj)                                 :: mp
      TYPE(dbcsr_type)                                   :: mat_g, mat_v, sparse_matrix2, &
                                                            sparse_matrix3

      check_product = .FALSE.

      CALL timeset(routineN, timing_handle)

      my_keep_sparsity = .TRUE.
      IF (PRESENT(keep_sparsity)) my_keep_sparsity = keep_sparsity
      NULLIFY (col_dist_left)

      IF (ncol .GT. 0) THEN
         IF (.NOT. dbcsr_valid_index(sparse_matrix)) &
            CPABORT("sparse_matrix must pre-exist")
         !
         ! Setup matrix_v
         CALL cp_fm_get_info(matrix_v, ncol_global=k)
         !WRITE(*,*)routineN//'truncated mult k, ncol',k,ncol,' PRESENT (matrix_g)',PRESENT (matrix_g)
         mp = dbcsr_distribution_mp(dbcsr_distribution(sparse_matrix))
         CALL create_bl_distribution(col_dist_left, col_blk_size_left, &
                                     k, dbcsr_mp_npcols(mp))
         NULLIFY (col_dist_cluster)
         CALL dbcsr_distribution_new(dist_left, mp, &
                                        dbcsr_distribution_row_dist(dbcsr_distribution(sparse_matrix)), &
                                        col_dist_left, &
                                        dbcsr_distribution_row_clusters(dbcsr_distribution(sparse_matrix)), &
                                        col_dist_cluster)
         DEALLOCATE (col_dist_left)
         CALL dbcsr_init(mat_v)
         CALL dbcsr_create(mat_v, "DBCSR matrix_v", dist_left, dbcsr_type_no_symmetry, &
                              dbcsr_row_block_sizes(sparse_matrix), col_blk_size_left, nze=0, &
                              data_type=dbcsr_get_data_type(sparse_matrix))
         CALL copy_fm_to_dbcsr(matrix_v, mat_v)
         CALL dbcsr_verify_matrix(mat_v)
         !
         ! Setup matrix_g
         IF (PRESENT(matrix_g)) THEN
            CALL dbcsr_init(mat_g)
            CALL dbcsr_create(mat_g, "DBCSR matrix_g", dist_left, &
                                 dbcsr_type_no_symmetry, &
                                 dbcsr_row_block_sizes(sparse_matrix), &
                                 dbcsr_col_block_sizes(mat_v), &
                                 data_type=dbcsr_get_data_type(sparse_matrix))
            CALL copy_fm_to_dbcsr(matrix_g, mat_g)
         ENDIF
         !
         DEALLOCATE (col_blk_size_left)
         CALL dbcsr_distribution_release(dist_left)
         !
         !
         IF (check_product) THEN
            NULLIFY (fm_matrix)
            CALL cp_fm_get_info(matrix_v, nrow_global=nao)
            CALL cp_fm_struct_create(fm_struct_tmp, context=matrix_v%matrix_struct%context, nrow_global=nao, &
                                     ncol_global=nao, para_env=matrix_v%matrix_struct%para_env)
            CALL cp_fm_create(fm_matrix, fm_struct_tmp, name="fm matrix")
            CALL cp_fm_struct_release(fm_struct_tmp)
            CALL copy_dbcsr_to_fm(sparse_matrix, fm_matrix)
            CALL dbcsr_init(sparse_matrix3)
            CALL dbcsr_copy(sparse_matrix3, sparse_matrix)
         ENDIF
         !
         my_alpha = 1.0_dp
         IF (PRESENT(alpha)) my_alpha = alpha
         IF (PRESENT(matrix_g)) THEN
            CALL dbcsr_multiply("N", "T", my_alpha, mat_v, mat_g, &
                                   1.0_dp, sparse_matrix, &
                                   retain_sparsity=my_keep_sparsity, &
                                   last_k=ncol)
         ELSE
            CALL dbcsr_multiply("N", "T", my_alpha, mat_v, mat_v, &
                                   1.0_dp, sparse_matrix, &
                                   retain_sparsity=my_keep_sparsity, &
                                   last_k=ncol)
         ENDIF

         IF (check_product) THEN
            IF (PRESENT(matrix_g)) THEN
               CALL cp_fm_gemm("N", "T", nao, nao, ncol, my_alpha, matrix_v, matrix_g, &
                               1.0_dp, fm_matrix)
            ELSE
               CALL cp_fm_gemm("N", "T", nao, nao, ncol, my_alpha, matrix_v, matrix_v, &
                               1.0_dp, fm_matrix)
            ENDIF

            CALL dbcsr_init(sparse_matrix2)
            CALL dbcsr_copy(sparse_matrix2, sparse_matrix)
            CALL dbcsr_scale(sparse_matrix2, alpha_scalar=0.0_dp)
            CALL copy_fm_to_dbcsr(fm_matrix, sparse_matrix2, keep_sparsity=my_keep_sparsity)
            CALL dbcsr_add(sparse_matrix2, sparse_matrix, alpha_scalar=1.0_dp, &
                              beta_scalar=-1.0_dp)
            CALL dbcsr_norm(sparse_matrix2, which_norm=dbcsr_norm_frobenius, &
                               norm_scalar=norm)
            WRITE (*, *) 'nao=', nao, ' k=', k, ' ncol=', ncol, ' my_alpha=', my_alpha
            WRITE (*, *) 'PRESENT (matrix_g)', PRESENT(matrix_g)
            WRITE (*, *) 'matrix_type=', dbcsr_get_matrix_type(sparse_matrix)
            WRITE (*, *) 'norm(sm+alpha*v*g^t - fm+alpha*v*g^t)/n=', norm/REAL(nao, dp)
            IF (norm/REAL(nao, dp) .GT. 1e-12_dp) THEN
               !WRITE(*,*) 'fm_matrix'
               !DO j=1,SIZE(fm_matrix%local_data,2)
               !   DO i=1,SIZE(fm_matrix%local_data,1)
               !      WRITE(*,'(A,I3,A,I3,A,E26.16,A)') 'a(',i,',',j,')=',fm_matrix%local_data(i,j),';'
               !   ENDDO
               !ENDDO
               !WRITE(*,*) 'mat_v'
               !CALL dbcsr_print(mat_v,matlab_format=.TRUE.)
               !WRITE(*,*) 'mat_g'
               !CALL dbcsr_print(mat_g,matlab_format=.TRUE.)
               !WRITE(*,*) 'sparse_matrix'
               !CALL dbcsr_print(sparse_matrix,matlab_format=.TRUE.)
               !WRITE(*,*) 'sparse_matrix2 (-sm + sparse(fm))'
               !CALL dbcsr_print(sparse_matrix2,matlab_format=.TRUE.)
               !WRITE(*,*) 'sparse_matrix3 (copy of sm input)'
               !CALL dbcsr_print(sparse_matrix3,matlab_format=.TRUE.)
               !stop
            ENDIF
            CALL dbcsr_release(sparse_matrix2)
            CALL dbcsr_release(sparse_matrix3)
            CALL cp_fm_release(fm_matrix)
         ENDIF
         CALL dbcsr_release(mat_v)
         IF (PRESENT(matrix_g)) CALL dbcsr_release(mat_g)
      ENDIF
      CALL timestop(timing_handle)

   END SUBROUTINE cp_dbcsr_plus_fm_fm_t

! **************************************************************************************************
!> \brief Utility function to copy a specially shaped fm to dbcsr_matrix
!>        The result matrix will be the matrix in dbcsr format
!>        with the row blocks sizes according to the block_sizes of the template
!>        and the col blocks sizes evenly blocked with the internal dbcsr conversion
!>        size (32 is the current default)
!> \param matrix ...
!> \param fm_in ...
!> \param template ...
! **************************************************************************************************
   SUBROUTINE cp_fm_to_dbcsr_row_template(matrix, fm_in, template)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      TYPE(cp_fm_type), POINTER                          :: fm_in
      TYPE(dbcsr_type), INTENT(IN)                       :: template

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_to_dbcsr_row_template', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: k_in
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size_right_in
      TYPE(dbcsr_distribution_obj)                       :: dist_right_in

      CALL cp_fm_get_info(fm_in, ncol_global=k_in)

      CALL dbcsr_create_dist_r_unrot(dist_right_in, dbcsr_distribution(template), k_in, &
                                        col_blk_size_right_in)
      CALL dbcsr_init(matrix)
      CALL dbcsr_create(matrix, "D", dist_right_in, dbcsr_type_no_symmetry, &
                           dbcsr_row_block_sizes(template), col_blk_size_right_in, &
                           nze=0, data_type=dbcsr_get_data_type(template))

      CALL copy_fm_to_dbcsr(fm_in, matrix)
      DEALLOCATE (col_blk_size_right_in)
      CALL dbcsr_distribution_release(dist_right_in)

   END SUBROUTINE cp_fm_to_dbcsr_row_template

! **************************************************************************************************
!> \brief Utility function to create an arbitrary shaped dbcsr matrix
!>        with the same processor grid as the template matrix
!>        both row sizes and col sizes are evenly blocked with the internal
!>        dbcsr_conversion size (32 is the current default)
!> \param matrix dbcsr matrix to be created
!> \param template template dbcsr matrix giving its mp_env
!> \param m global row size of output matrix
!> \param n global col size of output matrix
!> \param sym ...
!> \param data_type ...
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_m_by_n_from_template(matrix, template, m, n, sym, data_type)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix, template
      INTEGER                                            :: m, n
      CHARACTER, OPTIONAL                                :: sym
      INTEGER, OPTIONAL                                  :: data_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_m_by_n_from_template', &
         routineP = moduleN//':'//routineN

      CHARACTER                                          :: mysym
      INTEGER                                            :: my_data_type
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, col_cluster_dist, &
                                                            col_dist, row_blk_size, &
                                                            row_cluster_dist, row_dist
      TYPE(dbcsr_distribution_obj)                       :: dist_m_n

      mysym = dbcsr_get_matrix_type(template)
      IF (PRESENT(sym)) mysym = sym
      my_data_type = dbcsr_get_data_type(template)
      IF (PRESENT(data_type)) my_data_type = data_type

      NULLIFY (row_dist, col_dist)
      NULLIFY (row_blk_size, col_blk_size)
      NULLIFY (row_cluster_dist, col_cluster_dist)

      CALL create_bl_distribution(row_dist, row_blk_size, m, &
                                  dbcsr_mp_nprows(dbcsr_distribution_mp(dbcsr_distribution(template))))
      CALL create_bl_distribution(col_dist, col_blk_size, n, &
                                  dbcsr_mp_npcols(dbcsr_distribution_mp(dbcsr_distribution(template))))
      CALL dbcsr_distribution_new(dist_m_n, dbcsr_distribution_mp(dbcsr_distribution(template)), &
                                     row_dist, col_dist, &
                                     row_cluster_dist, col_cluster_dist, &
                                     reuse_arrays=.TRUE.)

      CALL dbcsr_create(matrix, "m_n_template", dist_m_n, mysym, &
                           row_blk_size, col_blk_size, nze=0, data_type=my_data_type, &
                           reuse_arrays=.TRUE.)
      CALL dbcsr_distribution_release(dist_m_n)

   END SUBROUTINE cp_dbcsr_m_by_n_from_template

! **************************************************************************************************
!> \brief Utility function to create dbcsr matrix, m x n matrix (n arbitrary)
!>        with the same processor grid and row distribution  as the template matrix
!>        col sizes are evenly blocked with the internal
!>        dbcsr_conversion size (32 is the current default)
!> \param matrix dbcsr matrix to be created
!> \param template template dbcsr matrix giving its mp_env
!> \param n global col size of output matrix
!> \param sym ...
!> \param data_type ...
! **************************************************************************************************
   SUBROUTINE cp_dbcsr_m_by_n_from_row_template(matrix, template, n, sym, data_type)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix, template
      INTEGER                                            :: n
      CHARACTER, OPTIONAL                                :: sym
      INTEGER, OPTIONAL                                  :: data_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_m_by_n_from_row_template', &
         routineP = moduleN//':'//routineN

      CHARACTER                                          :: mysym
      INTEGER                                            :: my_data_type
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, col_cluster_dist, col_dist
      TYPE(dbcsr_distribution_obj)                       :: dist_m_n

      mysym = dbcsr_get_matrix_type(template)
      IF (PRESENT(sym)) mysym = sym
      my_data_type = dbcsr_get_data_type(template)
      IF (PRESENT(data_type)) my_data_type = data_type

      NULLIFY (col_dist, col_blk_size)

      CALL create_bl_distribution(col_dist, col_blk_size, n, &
                                  dbcsr_mp_npcols(dbcsr_distribution_mp(dbcsr_distribution(template))))
      NULLIFY (col_cluster_dist)
      CALL dbcsr_distribution_new(dist_m_n, dbcsr_distribution_mp(dbcsr_distribution(template)), &
                                     dbcsr_distribution_row_dist(dbcsr_distribution(template)), col_dist, &
                                     dbcsr_distribution_row_clusters(dbcsr_distribution(template)), &
                                     col_cluster_dist)

      CALL dbcsr_create(matrix, "m_n_template", dist_m_n, mysym, &
                           dbcsr_row_block_sizes(template), col_blk_size, nze=0, data_type=my_data_type)

      DEALLOCATE (col_dist, col_blk_size)
      CALL dbcsr_distribution_release(dist_m_n)

   END SUBROUTINE cp_dbcsr_m_by_n_from_row_template

! **************************************************************************************************
!> \brief Distributes elements into blocks and into bins
!>
!> \param[out] block_distribution       block distribution to bins
!> \param[out] block_size       sizes of blocks
!> \param[in] nelements number of elements to bin
!> \param[in] nbins             number of bins
!> \par Term clarification
!>      An example: blocks are atom blocks and bins are process rows/colums.
! **************************************************************************************************
  SUBROUTINE create_bl_distribution (block_distribution,&
       block_size, nelements, nbins)
      INTEGER, DIMENSION(:), INTENT(OUT), POINTER        :: block_distribution, block_size
      INTEGER, INTENT(IN)                                :: nelements, nbins

      CHARACTER(len=*), PARAMETER :: routineN = 'create_bl_distribution', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: bin, blk_layer, element_stack, els, &
                                                            estimated_blocks, max_blocks_per_bin, &
                                                            nblks, nblocks, stat
      INTEGER, DIMENSION(:), POINTER                     :: blk_dist, blk_sizes

!   ---------------------------------------------------------------------------

    NULLIFY (block_distribution)
    NULLIFY (block_size)
    ! Define the sizes on which we build the distribution.
    IF (nelements .GT. 0) THEN

       nblocks = CEILING(REAL(nelements, KIND=dp)/REAL(max_elements_per_block, KIND=dp))
       max_blocks_per_bin = CEILING(REAL(nblocks, KIND=dp)/REAL(nbins, KIND=dp))

       IF (debug_mod) THEN
          WRITE(*,'(1X,A,1X,A,I7,A,I7,A)')routineP,"For",nelements,&
               " elements and",nbins," bins"
          WRITE(*,'(1X,A,1X,A,I7,A)')routineP,"There are",&
               max_elements_per_block," max elements per block"
          WRITE(*,'(1X,A,1X,A,I7,A)')routineP,"There are",&
               nblocks," blocks"
          WRITE(*,'(1X,A,1X,A,I7,A)')routineP,"There are",&
               max_blocks_per_bin," max blocks/bin"
       ENDIF

       estimated_blocks = max_blocks_per_bin * nbins
       ALLOCATE (blk_dist (estimated_blocks), stat=stat)
       IF(stat /= 0)&
          CPABORT("blk_dist")
       ALLOCATE (blk_sizes (estimated_blocks), stat=stat)
       IF(stat /= 0)&
          CPABORT("blk_sizes")
       element_stack = 0
       nblks = 0
       DO blk_layer = 1, max_blocks_per_bin
          DO bin = 0, nbins-1
             els = MIN (max_elements_per_block, nelements-element_stack)
             IF (els .GT. 0) THEN
                element_stack = element_stack + els
                nblks = nblks + 1
                blk_dist(nblks) = bin
                blk_sizes(nblks) = els
                IF (debug_mod) WRITE(*,'(1X,A,I5,A,I5,A,I5)')routineP//" Assigning",&
                     els," elements as block",nblks," to bin",bin
             ENDIF
          ENDDO
       ENDDO
       ! Create the output arrays.
       IF (nblks .EQ. estimated_blocks) THEN
          block_distribution => blk_dist
          block_size => blk_sizes
       ELSE
          ALLOCATE (block_distribution (nblks), stat=stat)
          IF(stat /= 0)&
             CPABORT("blk_dist")
          block_distribution(:)=blk_dist(1:nblks)
          DEALLOCATE (blk_dist)
          ALLOCATE (block_size (nblks), stat=stat)
          IF(stat /= 0)&
             CPABORT("blk_sizes")
          block_size(:)=blk_sizes(1:nblks)
          DEALLOCATE(blk_sizes)
       ENDIF
    ELSE
       ALLOCATE (block_distribution (0), stat=stat)
       IF(stat /= 0)&
          CPABORT("blk_dist")
       ALLOCATE (block_size (0), stat=stat)
       IF(stat /= 0)&
          CPABORT("blk_sizes")
    ENDIF
1579 FORMAT (I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5)
    IF (debug_mod) THEN
       WRITE(*,'(1X,A,A)')routineP//" Distribution"
       WRITE(*,1579) block_distribution(:)
       WRITE(*,'(1X,A,A)')routineP//" Sizes"
       WRITE(*,1579) block_size(:)
    ENDIF
  END SUBROUTINE create_bl_distribution

! **************************************************************************************************
!> \brief Creates a new distribution for the right matrix in a matrix
!>        multiplication with unrotated grid.
!> \param[out] dist_right     new distribution for the right matrix
!> \param[in] dist_left       the distribution of the left matrix
!> \param[in] ncolumns        number of columns in right matrix
!> \param[out] right_col_blk_sizes      sizes of blocks in the created column
!> \par The new row distribution for the right matrix is the same as the row
!>      distribution of the left matrix, while the column distribution is
!>      created so that it is appropriate to the parallel environment.
! **************************************************************************************************
  SUBROUTINE dbcsr_create_dist_r_unrot(dist_right, dist_left, ncolumns,&
       right_col_blk_sizes)
      TYPE(dbcsr_distribution_obj), INTENT(OUT)          :: dist_right
      TYPE(dbcsr_distribution_obj), INTENT(IN)           :: dist_left
      INTEGER, INTENT(IN)                                :: ncolumns
      INTEGER, DIMENSION(:), INTENT(OUT), POINTER        :: right_col_blk_sizes

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_r_unrot', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: multiplicity, new_npcols, new_nprows, &
                                                            nimages
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: tmp_images
      INTEGER, DIMENSION(:), POINTER                     :: dummy, right_col_dist, right_row_dist
      TYPE(dbcsr_mp_obj)                                 :: right_mp_env

!   ---------------------------------------------------------------------------

    right_mp_env = dbcsr_distribution_mp (dist_left)
    ! Create the column distribution
    CALL create_bl_distribution (right_col_dist,&
         right_col_blk_sizes, ncolumns, dbcsr_mp_npcols (right_mp_env))
    ! Create an even row distribution.
    ALLOCATE (right_row_dist(dbcsr_distribution_ncols (dist_left)))
    ALLOCATE (tmp_images(dbcsr_distribution_ncols (dist_left)))
    new_nprows = dbcsr_mp_nprows (right_mp_env)
    new_npcols = dbcsr_mp_npcols (right_mp_env)
    nimages = lcm (new_nprows, new_npcols) / new_nprows
    multiplicity = new_nprows / gcd (new_nprows, new_npcols)
    CALL rebin_distribution (right_row_dist, tmp_images,& 
         dbcsr_distribution_col_dist (dist_left),&
         new_nprows, multiplicity, nimages)

    NULLIFY(dummy)
    CALL dbcsr_distribution_new (dist_right,&
                                    mp_env=right_mp_env,&
                                    row_dist_block=right_row_dist,&
                                    col_dist_block=right_col_dist,&
                                    row_dist_cluster=dummy,&
                                    col_dist_cluster=dummy,&
                                    reuse_arrays=.TRUE.)
    DEALLOCATE(tmp_images)
  END SUBROUTINE dbcsr_create_dist_r_unrot

! **************************************************************************************************
!> \brief Makes new distribution with decimation and multiplicity
!> \param[out] new_bins      new real distribution
!> \param[out] images        new image distribution
!> \param[in] source_bins    Basis for the new distribution and images
!> \param[in] nbins          number of bins in the new real distribution
!> \param[in] multiplicity   multiplicity
!> \param[in] nimages        number of images in the new distribution
!> \par Definition of multiplicity and nimages
!>      Multiplicity and decimation (number of images) are used to
!>      match process grid coordinates on non-square process
!>      grids. Given source_nbins and target_nbins, their relation is
!>                source_nbins * target_multiplicity
!>              = target_nbins * target_nimages.
!>      It is best when both multiplicity and nimages are small. To
!>      get these two factors, then, one can use the following formulas:
!>          nimages      = lcm(source_nbins, target_nbins) / target_nbins
!>          multiplicity = target_nbins / gcd(source_nbins, target_nbins)
!>      from the target's point of view (nimages = target_nimages).
!> \par Mapping
!>      The new distribution comprises of real bins and images within
!>      bins. These can be view as target_nbins*nimages virtual
!>      columns. These same virtual columns are also
!>      source_nbins*multiplicity in number. Therefore these virtual
!>      columns are mapped from source_nbins*multiplicity onto
!>      target_bins*nimages (each target bin has nimages images):
!>      Source 4: |1 2 3|4 5 6|7 8 9|A B C| (4*3)
!>      Target 6: |1 2|3 4|5 6|7 8|9 A|B C| (6*2)
!>      multiplicity=3, nimages=2, 12 virtual columns (1-C).
!>      Source bin elements are evenly mapped into one of multiplicity
!>      virtual columns. Other (non-even, block-size aware) mappings
!>      could be better.
! **************************************************************************************************
  SUBROUTINE rebin_distribution (new_bins, images, source_bins,&
       nbins, multiplicity, nimages)
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: new_bins, images
      INTEGER, DIMENSION(:), INTENT(IN)                  :: source_bins
      INTEGER, INTENT(IN)                                :: nbins, multiplicity, nimages

      CHARACTER(len=*), PARAMETER :: routineN = 'rebin_distribution', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: bin, i, old_nbins, virtual_bin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bin_multiplier

!   ---------------------------------------------------------------------------

    IF(MOD(nbins*nimages,multiplicity).NE.0)&
       CPWARN("mulitplicity is not divisor of new process grid coordinate")
    old_nbins = (nbins * nimages) / multiplicity
    ALLOCATE (bin_multiplier(0:old_nbins-1))
    bin_multiplier(:) = 0
    DO i = 1, SIZE (new_bins)
       IF (i .LE. SIZE (source_bins)) THEN
          bin = source_bins (i)
       ELSE
          ! Fill remainder with a cyclic distribution
          bin = MOD (i, old_nbins)
       ENDIF
       virtual_bin = bin*multiplicity+bin_multiplier(bin)
       new_bins(i) = virtual_bin / nimages
       images(i) = 1 + MOD (virtual_bin, nimages)
       bin_multiplier(bin) = bin_multiplier(bin)+1
       IF (bin_multiplier(bin) .GE. multiplicity) THEN
          bin_multiplier(bin) = 0
       ENDIF
    ENDDO
  END SUBROUTINE rebin_distribution

! **************************************************************************************************
!> \brief Creates a block-cyclic compatible distribution
!>
!>        All blocks in a dimension, except for possibly the last
!>        block, have the same size.
!> \param[out] dist           the elemental distribution
!> \param[in] nrows           number of full rows
!> \param[in] ncolumns        number of full columns
!> \param[in] nrow_block      size of row blocks
!> \param[in] ncol_block      size of column blocks
!> \param[in] mp_env          multiprocess environment
!> \param[out] row_blk_sizes  row block sizes
!> \param[out] col_blk_sizes  column block sizes
! **************************************************************************************************
  SUBROUTINE dbcsr_create_dist_block_cyclic (dist, nrows, ncolumns,&
       nrow_block, ncol_block, mp_env, row_blk_sizes, col_blk_sizes)
      TYPE(dbcsr_distribution_obj), INTENT(OUT)          :: dist
      INTEGER, INTENT(IN)                                :: nrows, ncolumns, nrow_block, ncol_block
      TYPE(dbcsr_mp_obj), INTENT(IN)                     :: mp_env
      INTEGER, DIMENSION(:), INTENT(OUT), POINTER        :: row_blk_sizes, col_blk_sizes

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_block_cyclic', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: nblkcols, nblkrows, npcols, nprows, &
                                                            pdim, sz
      INTEGER, DIMENSION(:), POINTER                     :: cd_data, dummy, rd_data

!   ---------------------------------------------------------------------------
! Row sizes

    IF(nrow_block.EQ.0) THEN
       nblkrows = 0
       sz = 0
    ELSE
       nblkrows = nrows / nrow_block
       sz = MOD(nrows, nrow_block)
    ENDIF
    IF (sz .GT. 0) nblkrows = nblkrows + 1
    ALLOCATE (row_blk_sizes(nblkrows), rd_data(nblkrows))
    row_blk_sizes = nrow_block
    IF (sz .NE. 0) row_blk_sizes(nblkrows) = sz

    ! Column sizes
    IF(ncol_block.EQ.0) THEN
       nblkcols = 0
       sz = 0
    ELSE
       nblkcols = ncolumns / ncol_block
       sz = MOD(ncolumns, ncol_block)
    ENDIF
    IF (sz .GT. 0) nblkcols = nblkcols + 1
    ALLOCATE (col_blk_sizes(nblkcols), cd_data(nblkcols))
    col_blk_sizes = ncol_block
    IF (sz .NE. 0) col_blk_sizes(nblkcols) = sz
    !
    IF(debug_mod) THEN
       WRITE(*,*)routineN//" nrows,nrow_block,nblkrows=",&
            nrows,nrow_block,nblkrows
       WRITE(*,*)routineN//" ncols,ncol_block,nblkcols=",&
            ncolumns,ncol_block,nblkcols
    ENDIF
    ! Calculate process row distribution
    nprows = dbcsr_mp_nprows (mp_env)
    DO pdim = 0 , MIN(nprows-1,nblkrows-1)
       rd_data(1+pdim:nblkrows:nprows) = pdim
    END DO
    ! Calculate process column distribution
    npcols = dbcsr_mp_npcols (mp_env)
    DO pdim = 0 , MIN(npcols-1, nblkcols-1)
       cd_data(1+pdim:nblkcols:npcols) = pdim
    END DO
    !
    IF(debug_mod) THEN
       WRITE(*,*)routineN//" row_dist",&
            rd_data
       WRITE(*,*)routineN//" col_dist",&
            cd_data
    ENDIF
    !
    NULLIFY(dummy)
    CALL dbcsr_distribution_new (dist,&
                                    mp_env=mp_env,&
                                    row_dist_block=rd_data,&
                                    col_dist_block=cd_data,&
                                    row_dist_cluster=dummy,&
                                    col_dist_cluster=dummy,&
                                    reuse_arrays=.TRUE.)


  END SUBROUTINE dbcsr_create_dist_block_cyclic

END MODULE cp_dbcsr_operations

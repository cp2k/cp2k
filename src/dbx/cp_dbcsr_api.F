!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE cp_dbcsr_api
   USE dbcsr_api,                       ONLY: &
        convert_csr_to_dbcsr_prv => dbcsr_convert_csr_to_dbcsr, &
        convert_dbcsr_to_csr_prv => dbcsr_convert_dbcsr_to_csr, dbcsr_add_prv => dbcsr_add, &
        dbcsr_binary_read_prv => dbcsr_binary_read, dbcsr_binary_write_prv => dbcsr_binary_write, &
        dbcsr_clear_mempools, dbcsr_clear_prv => dbcsr_clear, &
        dbcsr_complete_redistribute_prv => dbcsr_complete_redistribute, &
        dbcsr_convert_offsets_to_sizes, dbcsr_convert_sizes_to_offsets, &
        dbcsr_copy_prv => dbcsr_copy, dbcsr_create_prv => dbcsr_create, dbcsr_csr_create, &
        dbcsr_csr_create_from_dbcsr_prv => dbcsr_csr_create_from_dbcsr, &
        dbcsr_csr_dbcsr_blkrow_dist, dbcsr_csr_destroy, dbcsr_csr_eqrow_floor_dist, &
        dbcsr_csr_p_type, dbcsr_csr_print_sparsity, dbcsr_csr_type, &
        dbcsr_csr_type_real_8 => dbcsr_type_real_8, dbcsr_csr_write, &
        dbcsr_desymmetrize_prv => dbcsr_desymmetrize, dbcsr_distribute_prv => dbcsr_distribute, &
        dbcsr_distribution_get_num_images, dbcsr_distribution_get_prv => dbcsr_distribution_get, &
        dbcsr_distribution_hold_prv => dbcsr_distribution_hold, &
        dbcsr_distribution_new_prv => dbcsr_distribution_new, &
        dbcsr_distribution_release_prv => dbcsr_distribution_release, &
        dbcsr_distribution_type_prv => dbcsr_distribution_type, dbcsr_dot_prv => dbcsr_dot, &
        dbcsr_filter_prv => dbcsr_filter, dbcsr_finalize_lib, &
        dbcsr_finalize_prv => dbcsr_finalize, dbcsr_get_block_p_prv => dbcsr_get_block_p, &
        dbcsr_get_data_p_prv => dbcsr_get_data_p, dbcsr_get_data_size_prv => dbcsr_get_data_size, &
        dbcsr_get_default_config, dbcsr_get_info_prv => dbcsr_get_info, &
        dbcsr_get_matrix_type_prv => dbcsr_get_matrix_type, &
        dbcsr_get_num_blocks_prv => dbcsr_get_num_blocks, &
        dbcsr_get_occupation_prv => dbcsr_get_occupation, &
        dbcsr_get_stored_coordinates_prv => dbcsr_get_stored_coordinates, &
        dbcsr_has_symmetry_prv => dbcsr_has_symmetry, dbcsr_init_lib, &
        dbcsr_iterator_blocks_left_prv => dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block_prv => dbcsr_iterator_next_block, &
        dbcsr_iterator_start_prv => dbcsr_iterator_start, &
        dbcsr_iterator_stop_prv => dbcsr_iterator_stop, &
        dbcsr_iterator_type_prv => dbcsr_iterator_type, &
        dbcsr_mp_grid_setup_prv => dbcsr_mp_grid_setup, dbcsr_multiply_prv => dbcsr_multiply, &
        dbcsr_no_transpose, dbcsr_print_config, dbcsr_print_statistics, &
        dbcsr_put_block_prv => dbcsr_put_block, dbcsr_release_prv => dbcsr_release, &
        dbcsr_replicate_all_prv => dbcsr_replicate_all, &
        dbcsr_reserve_blocks_prv => dbcsr_reserve_blocks, dbcsr_reset_randmat_seed, &
        dbcsr_run_tests, dbcsr_scale_prv => dbcsr_scale, dbcsr_set_config, &
        dbcsr_set_prv => dbcsr_set, dbcsr_sum_replicated_prv => dbcsr_sum_replicated, &
        dbcsr_test_mm, dbcsr_transpose, dbcsr_transposed_prv => dbcsr_transposed, &
        dbcsr_type_antisymmetric, dbcsr_type_complex_8, dbcsr_type_no_symmetry, &
        dbcsr_type_prv => dbcsr_type, dbcsr_type_real_8, dbcsr_type_symmetric, &
        dbcsr_valid_index_prv => dbcsr_valid_index, &
        dbcsr_verify_matrix_prv => dbcsr_verify_matrix, dbcsr_work_create_prv => dbcsr_work_create
   USE dbm_api,                         ONLY: &
        dbm_add, dbm_clear, dbm_copy, dbm_distribution_obj, dbm_iterator, dbm_redistribute, &
        dbm_scale, dbm_type, dbm_zero
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_comm_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! constants
   PUBLIC :: dbcsr_type_no_symmetry
   PUBLIC :: dbcsr_type_symmetric
   PUBLIC :: dbcsr_type_antisymmetric
   PUBLIC :: dbcsr_transpose
   PUBLIC :: dbcsr_no_transpose

   ! types
   PUBLIC :: dbcsr_type
   PUBLIC :: dbcsr_p_type
   PUBLIC :: dbcsr_distribution_type
   PUBLIC :: dbcsr_iterator_type

   ! lib init/finalize
   PUBLIC :: dbcsr_clear_mempools
   PUBLIC :: dbcsr_init_lib
   PUBLIC :: dbcsr_finalize_lib
   PUBLIC :: dbcsr_set_config
   PUBLIC :: dbcsr_get_default_config
   PUBLIC :: dbcsr_print_config
   PUBLIC :: dbcsr_reset_randmat_seed
   PUBLIC :: dbcsr_mp_grid_setup
   PUBLIC :: dbcsr_print_statistics

   ! create / release
   PUBLIC :: dbcsr_distribution_hold
   PUBLIC :: dbcsr_distribution_release
   PUBLIC :: dbcsr_distribution_new
   PUBLIC :: dbcsr_create
   PUBLIC :: dbcsr_init_p
   PUBLIC :: dbcsr_release
   PUBLIC :: dbcsr_release_p
   PUBLIC :: dbcsr_deallocate_matrix

   ! primitive matrix operations
   PUBLIC :: dbcsr_set
   PUBLIC :: dbcsr_add
   PUBLIC :: dbcsr_scale
   PUBLIC :: dbcsr_transposed
   PUBLIC :: dbcsr_multiply
   PUBLIC :: dbcsr_copy
   PUBLIC :: dbcsr_desymmetrize
   PUBLIC :: dbcsr_filter
   PUBLIC :: dbcsr_complete_redistribute
   PUBLIC :: dbcsr_reserve_blocks
   PUBLIC :: dbcsr_put_block
   PUBLIC :: dbcsr_get_block_p
   PUBLIC :: dbcsr_get_readonly_block_p
   PUBLIC :: dbcsr_clear

   ! iterator
   PUBLIC :: dbcsr_iterator_start
   PUBLIC :: dbcsr_iterator_readonly_start
   PUBLIC :: dbcsr_iterator_stop
   PUBLIC :: dbcsr_iterator_blocks_left
   PUBLIC :: dbcsr_iterator_next_block

   ! getters
   PUBLIC :: dbcsr_get_info
   PUBLIC :: dbcsr_distribution_get
   PUBLIC :: dbcsr_get_matrix_type
   PUBLIC :: dbcsr_get_occupation
   PUBLIC :: dbcsr_get_num_blocks
   PUBLIC :: dbcsr_get_data_size
   PUBLIC :: dbcsr_has_symmetry
   PUBLIC :: dbcsr_get_stored_coordinates
   PUBLIC :: dbcsr_valid_index

   ! work operations
   PUBLIC :: dbcsr_work_create
   PUBLIC :: dbcsr_verify_matrix
   PUBLIC :: dbcsr_get_data_p
   PUBLIC :: dbcsr_finalize

   ! replication
   PUBLIC :: dbcsr_replicate_all
   PUBLIC :: dbcsr_sum_replicated
   PUBLIC :: dbcsr_distribute

   ! misc
   PUBLIC :: dbcsr_distribution_get_num_images
   PUBLIC :: dbcsr_convert_offsets_to_sizes
   PUBLIC :: dbcsr_convert_sizes_to_offsets
   PUBLIC :: dbcsr_run_tests
   PUBLIC :: dbcsr_test_mm
   PUBLIC :: dbcsr_dot_threadsafe

   ! csr conversion
   PUBLIC :: dbcsr_csr_type
   PUBLIC :: dbcsr_csr_p_type
   PUBLIC :: dbcsr_convert_csr_to_dbcsr
   PUBLIC :: dbcsr_convert_dbcsr_to_csr
   PUBLIC :: dbcsr_csr_create_from_dbcsr
   PUBLIC :: dbcsr_csr_destroy
   PUBLIC :: dbcsr_csr_create
   PUBLIC :: dbcsr_csr_eqrow_floor_dist
   PUBLIC :: dbcsr_csr_dbcsr_blkrow_dist
   PUBLIC :: dbcsr_csr_print_sparsity
   PUBLIC :: dbcsr_csr_write
   PUBLIC :: dbcsr_csr_create_and_convert_complex
   PUBLIC :: dbcsr_csr_type_real_8

   ! binary io
   PUBLIC :: dbcsr_binary_write
   PUBLIC :: dbcsr_binary_read

   TYPE dbcsr_p_type
      TYPE(dbcsr_type), POINTER :: matrix => Null()
   END TYPE

   TYPE dbcsr_type
      PRIVATE
      TYPE(dbcsr_type_prv)  :: dbcsr = dbcsr_type_prv()
      TYPE(dbm_type) :: dbm = dbm_type()
   END TYPE dbcsr_type

   TYPE dbcsr_distribution_type
      PRIVATE
      TYPE(dbcsr_distribution_type_prv) :: dbcsr = dbcsr_distribution_type_prv()
      TYPE(dbm_distribution_obj)        :: dbm = dbm_distribution_obj()
   END TYPE dbcsr_distribution_type

   TYPE dbcsr_iterator_type
      PRIVATE
      TYPE(dbcsr_iterator_type_prv)  :: dbcsr = dbcsr_iterator_type_prv()
      TYPE(dbm_iterator) :: dbm = dbm_iterator()
   END TYPE dbcsr_iterator_type

   INTERFACE dbcsr_create
      MODULE PROCEDURE dbcsr_create_new, dbcsr_create_template
   END INTERFACE

   LOGICAL, PARAMETER, PRIVATE :: USE_DBCSR_BACKEND = .TRUE.

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_init_p(matrix)
      TYPE(dbcsr_type), POINTER                          :: matrix

      IF (ASSOCIATED(matrix)) THEN
         CALL dbcsr_release(matrix)
         DEALLOCATE (matrix)
      END IF

      ALLOCATE (matrix)
   END SUBROUTINE dbcsr_init_p

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_release_p(matrix)
      TYPE(dbcsr_type), POINTER                          :: matrix

      IF (ASSOCIATED(matrix)) THEN
         CALL dbcsr_release(matrix)
         DEALLOCATE (matrix)
      END IF
   END SUBROUTINE dbcsr_release_p

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_deallocate_matrix(matrix)
      TYPE(dbcsr_type), POINTER                          :: matrix

      CALL dbcsr_release(matrix)
      IF (dbcsr_valid_index(matrix)) &
         CALL cp_abort(__LOCATION__, &
                       'You should not "deallocate" a referenced matrix. '// &
                       'Avoid pointers to DBCSR matrices.')
      DEALLOCATE (matrix)
   END SUBROUTINE dbcsr_deallocate_matrix

! **************************************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param alpha_scalar ...
!> \param beta_scalar ...
! **************************************************************************************************
   SUBROUTINE dbcsr_add(matrix_a, matrix_b, alpha_scalar, beta_scalar)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_a
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_b
      REAL(kind=dp), INTENT(IN)                          :: alpha_scalar, beta_scalar

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_add_prv(matrix_a%dbcsr, matrix_b%dbcsr, alpha_scalar, beta_scalar)
      ELSE
         IF (alpha_scalar /= 1.0_dp .OR. beta_scalar /= 1.0_dp) CPABORT("Not yet implemented for DBM.")
         CALL dbm_add(matrix_a%dbm, matrix_b%dbm)
      END IF
   END SUBROUTINE dbcsr_add

! **************************************************************************************************
!> \brief ...
!> \param filepath ...
!> \param distribution ...
!> \param matrix_new ...
! **************************************************************************************************
   SUBROUTINE dbcsr_binary_read(filepath, distribution, matrix_new)
      CHARACTER(len=*), INTENT(IN)                       :: filepath
      TYPE(dbcsr_distribution_type), INTENT(IN)          :: distribution
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_new

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_binary_read_prv(filepath, distribution%dbcsr, matrix_new%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_binary_read

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param filepath ...
! **************************************************************************************************
   SUBROUTINE dbcsr_binary_write(matrix, filepath)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      CHARACTER(LEN=*), INTENT(IN)                       :: filepath

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_binary_write_prv(matrix%dbcsr, filepath)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_binary_write

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_clear(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_clear_prv(matrix%dbcsr)
      ELSE
         CALL dbm_clear(matrix%dbm)
      END IF
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param redist ...
! **************************************************************************************************
   SUBROUTINE dbcsr_complete_redistribute(matrix, redist)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(dbcsr_type), INTENT(INOUT)                    :: redist

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_complete_redistribute_prv(matrix%dbcsr, redist%dbcsr)
      ELSE
         CALL dbm_redistribute(matrix%dbm, redist%dbm)
      END IF
   END SUBROUTINE dbcsr_complete_redistribute

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_mat ...
!> \param csr_mat ...
! **************************************************************************************************
   SUBROUTINE dbcsr_convert_csr_to_dbcsr(dbcsr_mat, csr_mat)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: dbcsr_mat
      TYPE(dbcsr_csr_type), INTENT(INOUT)                :: csr_mat

      IF (USE_DBCSR_BACKEND) THEN
         CALL convert_csr_to_dbcsr_prv(dbcsr_mat%dbcsr, csr_mat)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_convert_csr_to_dbcsr

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_mat ...
!> \param csr_mat ...
! **************************************************************************************************
   SUBROUTINE dbcsr_convert_dbcsr_to_csr(dbcsr_mat, csr_mat)
      TYPE(dbcsr_type), INTENT(IN)                       :: dbcsr_mat
      TYPE(dbcsr_csr_type), INTENT(INOUT)                :: csr_mat

      IF (USE_DBCSR_BACKEND) THEN
         CALL convert_dbcsr_to_csr_prv(dbcsr_mat%dbcsr, csr_mat)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_convert_dbcsr_to_csr

! **************************************************************************************************
!> \brief ...
!> \param matrix_b ...
!> \param matrix_a ...
!> \param name ...
!> \param keep_sparsity ...
!> \param keep_imaginary ...
! **************************************************************************************************
   SUBROUTINE dbcsr_copy(matrix_b, matrix_a, name, keep_sparsity, keep_imaginary)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_b
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: name
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_sparsity, keep_imaginary

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_copy_prv(matrix_b%dbcsr, matrix_a%dbcsr, name=name, &
                             keep_sparsity=keep_sparsity, keep_imaginary=keep_imaginary)
      ELSE
         IF (PRESENT(name) .OR. PRESENT(keep_sparsity) .OR. PRESENT(keep_imaginary)) THEN
            CPABORT("Not yet implemented for DBM.")
         END IF
         CALL dbm_copy(matrix_b%dbm, matrix_a%dbm)
      END IF
   END SUBROUTINE dbcsr_copy

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param name ...
!> \param dist ...
!> \param matrix_type ...
!> \param row_blk_size ...
!> \param col_blk_size ...
!> \param reuse_arrays ...
!> \param mutable_work ...
! **************************************************************************************************
   SUBROUTINE dbcsr_create_new(matrix, name, dist, matrix_type, row_blk_size, col_blk_size, &
                               reuse_arrays, mutable_work)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      CHARACTER(len=*), INTENT(IN)                       :: name
      TYPE(dbcsr_distribution_type), INTENT(IN)          :: dist
      CHARACTER, INTENT(IN)                              :: matrix_type
      INTEGER, DIMENSION(:), INTENT(INOUT), POINTER      :: row_blk_size, col_blk_size
      LOGICAL, INTENT(IN), OPTIONAL                      :: reuse_arrays, mutable_work

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_create_prv(matrix=matrix%dbcsr, name=name, dist=dist%dbcsr, &
                               matrix_type=matrix_type, row_blk_size=row_blk_size, &
                               col_blk_size=col_blk_size, nze=0, data_type=dbcsr_type_real_8, &
                               reuse_arrays=reuse_arrays, mutable_work=mutable_work)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_create_new

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param name ...
!> \param template ...
!> \param dist ...
!> \param matrix_type ...
!> \param row_blk_size ...
!> \param col_blk_size ...
!> \param reuse_arrays ...
!> \param mutable_work ...
! **************************************************************************************************
   SUBROUTINE dbcsr_create_template(matrix, name, template, dist, matrix_type, &
                                    row_blk_size, col_blk_size, reuse_arrays, mutable_work)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: name
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(dbcsr_distribution_type), INTENT(IN), &
         OPTIONAL                                        :: dist
      CHARACTER, INTENT(IN), OPTIONAL                    :: matrix_type
      INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL, &
         POINTER                                         :: row_blk_size, col_blk_size
      LOGICAL, INTENT(IN), OPTIONAL                      :: reuse_arrays, mutable_work

      IF (USE_DBCSR_BACKEND) THEN
         IF (PRESENT(dist)) THEN
            CALL dbcsr_create_prv(matrix=matrix%dbcsr, name=name, template=template%dbcsr, &
                                  dist=dist%dbcsr, matrix_type=matrix_type, &
                                  row_blk_size=row_blk_size, col_blk_size=col_blk_size, &
                                  nze=0, data_type=dbcsr_type_real_8, reuse_arrays=reuse_arrays, &
                                  mutable_work=mutable_work)
         ELSE
            CALL dbcsr_create_prv(matrix=matrix%dbcsr, name=name, template=template%dbcsr, &
                                  matrix_type=matrix_type, &
                                  row_blk_size=row_blk_size, col_blk_size=col_blk_size, &
                                  nze=0, data_type=dbcsr_type_real_8, reuse_arrays=reuse_arrays, &
                                  mutable_work=mutable_work)
         END IF
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_create_template

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_mat ...
!> \param csr_mat ...
!> \param dist_format ...
!> \param csr_sparsity ...
!> \param numnodes ...
! **************************************************************************************************
   SUBROUTINE dbcsr_csr_create_from_dbcsr(dbcsr_mat, csr_mat, dist_format, csr_sparsity, numnodes)

      TYPE(dbcsr_type), INTENT(IN)                       :: dbcsr_mat
      TYPE(dbcsr_csr_type), INTENT(OUT)                  :: csr_mat
      INTEGER                                            :: dist_format
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: csr_sparsity
      INTEGER, INTENT(IN), OPTIONAL                      :: numnodes

      IF (USE_DBCSR_BACKEND) THEN
         IF (PRESENT(csr_sparsity)) THEN
            CALL dbcsr_csr_create_from_dbcsr_prv(dbcsr_mat%dbcsr, csr_mat, dist_format, &
                                                 csr_sparsity%dbcsr, numnodes)
         ELSE
            CALL dbcsr_csr_create_from_dbcsr_prv(dbcsr_mat%dbcsr, csr_mat, &
                                                 dist_format, numnodes=numnodes)
         END IF
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_csr_create_from_dbcsr

! **************************************************************************************************
!> \brief Combines csr_create_from_dbcsr and convert_dbcsr_to_csr to produce a complex CSR matrix.
!> \param rmatrix Real part of the matrix.
!> \param imatrix Imaginary part of the matrix.
!> \param csr_mat The resulting CSR matrix.
!> \param dist_format ...
! **************************************************************************************************
   SUBROUTINE dbcsr_csr_create_and_convert_complex(rmatrix, imatrix, csr_mat, dist_format)
      TYPE(dbcsr_type), INTENT(IN)                       :: rmatrix, imatrix
      TYPE(dbcsr_csr_type), INTENT(INOUT)                :: csr_mat
      INTEGER                                            :: dist_format

      COMPLEX(KIND=dp), PARAMETER :: ione = CMPLX(0.0_dp, 1.0_dp, KIND=dp), &
         rone = CMPLX(1.0_dp, 0.0_dp, KIND=dp)

      TYPE(dbcsr_type)                                   :: cmatrix, tmp_matrix

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_create_prv(tmp_matrix%dbcsr, template=rmatrix%dbcsr, data_type=dbcsr_type_complex_8)
         CALL dbcsr_create_prv(cmatrix%dbcsr, template=rmatrix%dbcsr, data_type=dbcsr_type_complex_8)
         CALL dbcsr_copy_prv(cmatrix%dbcsr, rmatrix%dbcsr)
         CALL dbcsr_copy_prv(tmp_matrix%dbcsr, imatrix%dbcsr)
         CALL dbcsr_add_prv(cmatrix%dbcsr, tmp_matrix%dbcsr, rone, ione)
         CALL dbcsr_release_prv(tmp_matrix%dbcsr)
         ! Convert to csr
         CALL dbcsr_csr_create_from_dbcsr_prv(cmatrix%dbcsr, csr_mat, dist_format)
         CALL convert_dbcsr_to_csr_prv(cmatrix%dbcsr, csr_mat)
         CALL dbcsr_release_prv(cmatrix%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_csr_create_and_convert_complex

! **************************************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
! **************************************************************************************************
   SUBROUTINE dbcsr_desymmetrize(matrix_a, matrix_b)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_b

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_desymmetrize_prv(matrix_a%dbcsr, matrix_b%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_desymmetrize

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_distribute(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_distribute_prv(matrix%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_distribute

! **************************************************************************************************
!> \brief ...
!> \param dist ...
!> \param row_dist ...
!> \param col_dist ...
!> \param nrows ...
!> \param ncols ...
!> \param has_threads ...
!> \param group ...
!> \param mynode ...
!> \param numnodes ...
!> \param nprows ...
!> \param npcols ...
!> \param myprow ...
!> \param mypcol ...
!> \param pgrid ...
!> \param subgroups_defined ...
!> \param prow_group ...
!> \param pcol_group ...
! **************************************************************************************************
   SUBROUTINE dbcsr_distribution_get(dist, row_dist, col_dist, nrows, ncols, has_threads, &
                                     group, mynode, numnodes, nprows, npcols, myprow, mypcol, &
                                     pgrid, subgroups_defined, prow_group, pcol_group)
      TYPE(dbcsr_distribution_type), INTENT(IN)          :: dist
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: row_dist, col_dist
      INTEGER, INTENT(OUT), OPTIONAL                     :: nrows, ncols
      LOGICAL, INTENT(OUT), OPTIONAL                     :: has_threads
      INTEGER, INTENT(OUT), OPTIONAL                     :: group, mynode, numnodes, nprows, npcols, &
                                                            myprow, mypcol
      INTEGER, DIMENSION(:, :), OPTIONAL, POINTER        :: pgrid
      LOGICAL, INTENT(OUT), OPTIONAL                     :: subgroups_defined
      INTEGER, INTENT(OUT), OPTIONAL                     :: prow_group, pcol_group

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_distribution_get_prv(dist%dbcsr, row_dist, col_dist, nrows, ncols, has_threads, &
                                         group, mynode, numnodes, nprows, npcols, myprow, mypcol, &
                                         pgrid, subgroups_defined, prow_group, pcol_group)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_distribution_get

! **************************************************************************************************
!> \brief ...
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE dbcsr_distribution_hold(dist)
      TYPE(dbcsr_distribution_type)                      :: dist

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_distribution_hold_prv(dist%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_distribution_hold

! **************************************************************************************************
!> \brief ...
!> \param dist ...
!> \param template ...
!> \param group ...
!> \param pgrid ...
!> \param row_dist ...
!> \param col_dist ...
!> \param reuse_arrays ...
! **************************************************************************************************
   SUBROUTINE dbcsr_distribution_new(dist, template, group, pgrid, row_dist, col_dist, reuse_arrays)
      TYPE(dbcsr_distribution_type), INTENT(OUT)         :: dist
      TYPE(dbcsr_distribution_type), INTENT(IN), &
         OPTIONAL                                        :: template
      INTEGER, INTENT(IN), OPTIONAL                      :: group
      INTEGER, DIMENSION(:, :), OPTIONAL, POINTER        :: pgrid
      INTEGER, DIMENSION(:), INTENT(INOUT), POINTER      :: row_dist, col_dist
      LOGICAL, INTENT(IN), OPTIONAL                      :: reuse_arrays

      IF (USE_DBCSR_BACKEND) THEN
         IF (PRESENT(template)) THEN
            CALL dbcsr_distribution_new_prv(dist%dbcsr, template%dbcsr, group, pgrid, &
                                            row_dist, col_dist, reuse_arrays)
         ELSE
            CALL dbcsr_distribution_new_prv(dist%dbcsr, group=group, pgrid=pgrid, &
                                            row_dist=row_dist, col_dist=col_dist, &
                                            reuse_arrays=reuse_arrays)
         END IF
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_distribution_new

! **************************************************************************************************
!> \brief ...
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE dbcsr_distribution_release(dist)
      TYPE(dbcsr_distribution_type)                      :: dist

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_distribution_release_prv(dist%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_distribution_release

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eps ...
! **************************************************************************************************
   SUBROUTINE dbcsr_filter(matrix, eps)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(dp), INTENT(IN)                               :: eps

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_filter_prv(matrix%dbcsr, eps)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_filter

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_finalize(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_finalize_prv(matrix%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_finalize

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param found ...
!> \param row_size ...
!> \param col_size ...
! **************************************************************************************************
   SUBROUTINE dbcsr_get_block_p(matrix, row, col, block, found, row_size, col_size)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN)                                :: row, col
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: block
      LOGICAL, INTENT(OUT)                               :: found
      INTEGER, INTENT(OUT), OPTIONAL                     :: row_size, col_size

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_get_block_p_prv(matrix%dbcsr, row, col, block, found, row_size, col_size)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_get_block_p

! **************************************************************************************************
!> \brief Like dbcsr_get_block_p() but with matrix being INTENT(IN).
!>        When invoking this routine, the caller promises not to modify the returned block.
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param found ...
!> \param row_size ...
!> \param col_size ...
! **************************************************************************************************
   SUBROUTINE dbcsr_get_readonly_block_p(matrix, row, col, block, found, row_size, col_size)
      TYPE(dbcsr_type), INTENT(IN), TARGET               :: matrix
      INTEGER, INTENT(IN)                                :: row, col
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: block
      LOGICAL, INTENT(OUT)                               :: found
      INTEGER, INTENT(OUT), OPTIONAL                     :: row_size, col_size

      TYPE(dbcsr_type), POINTER                          :: matrix_p

      MARK_USED(matrix)
      MARK_USED(row)
      MARK_USED(col)
      MARK_USED(block)
      MARK_USED(found)
      MARK_USED(row_size)
      MARK_USED(col_size)
      IF (USE_DBCSR_BACKEND) THEN
         matrix_p => matrix ! Hacky workaround to shake the INTENT(IN).
         CALL dbcsr_get_block_p_prv(matrix_p%dbcsr, row, col, block, found, row_size, col_size)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_get_readonly_block_p

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param lb ...
!> \param ub ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_get_data_p(matrix, lb, ub) RESULT(res)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER, INTENT(IN), OPTIONAL                      :: lb, ub
      REAL(kind=dp), DIMENSION(:), POINTER               :: res

      IF (USE_DBCSR_BACKEND) THEN
         res => dbcsr_get_data_p_prv(matrix%dbcsr, select_data_type=0.0_dp, lb=lb, ub=ub)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END FUNCTION dbcsr_get_data_p

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_get_data_size(matrix) RESULT(data_size)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: data_size

      IF (USE_DBCSR_BACKEND) THEN
         data_size = dbcsr_get_data_size_prv(matrix%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END FUNCTION dbcsr_get_data_size

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param nblkrows_total ...
!> \param nblkcols_total ...
!> \param nfullrows_total ...
!> \param nfullcols_total ...
!> \param nblkrows_local ...
!> \param nblkcols_local ...
!> \param nfullrows_local ...
!> \param nfullcols_local ...
!> \param my_prow ...
!> \param my_pcol ...
!> \param local_rows ...
!> \param local_cols ...
!> \param proc_row_dist ...
!> \param proc_col_dist ...
!> \param row_blk_size ...
!> \param col_blk_size ...
!> \param row_blk_offset ...
!> \param col_blk_offset ...
!> \param distribution ...
!> \param name ...
!> \param matrix_type ...
!> \param group ...
! **************************************************************************************************
   SUBROUTINE dbcsr_get_info(matrix, nblkrows_total, nblkcols_total, &
                             nfullrows_total, nfullcols_total, nblkrows_local, nblkcols_local, &
                             nfullrows_local, nfullcols_local, my_prow, my_pcol, &
                             local_rows, local_cols, proc_row_dist, proc_col_dist, &
                             row_blk_size, col_blk_size, row_blk_offset, col_blk_offset, &
                             distribution, name, matrix_type, group)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER, INTENT(OUT), OPTIONAL :: nblkrows_total, nblkcols_total, nfullrows_total, &
         nfullcols_total, nblkrows_local, nblkcols_local, nfullrows_local, nfullcols_local, &
         my_prow, my_pcol
      INTEGER, DIMENSION(:), OPTIONAL, POINTER :: local_rows, local_cols, proc_row_dist, &
         proc_col_dist, row_blk_size, col_blk_size, row_blk_offset, col_blk_offset
      TYPE(dbcsr_distribution_type), INTENT(OUT), &
         OPTIONAL                                        :: distribution
      CHARACTER(len=*), INTENT(OUT), OPTIONAL            :: name
      CHARACTER, INTENT(OUT), OPTIONAL                   :: matrix_type
      TYPE(mp_comm_type), INTENT(OUT), OPTIONAL          :: group

      INTEGER                                            :: group_handle
      TYPE(dbcsr_distribution_type_prv)                  :: my_distribution

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_get_info_prv(matrix=matrix%dbcsr, &
                                 nblkrows_total=nblkrows_total, &
                                 nblkcols_total=nblkcols_total, &
                                 nfullrows_total=nfullrows_total, &
                                 nfullcols_total=nfullcols_total, &
                                 nblkrows_local=nblkrows_local, &
                                 nblkcols_local=nblkcols_local, &
                                 nfullrows_local=nfullrows_local, &
                                 nfullcols_local=nfullcols_local, &
                                 my_prow=my_prow, &
                                 my_pcol=my_pcol, &
                                 local_rows=local_rows, &
                                 local_cols=local_cols, &
                                 proc_row_dist=proc_row_dist, &
                                 proc_col_dist=proc_col_dist, &
                                 row_blk_size=row_blk_size, &
                                 col_blk_size=col_blk_size, &
                                 row_blk_offset=row_blk_offset, &
                                 col_blk_offset=col_blk_offset, &
                                 distribution=my_distribution, &
                                 name=name, &
                                 matrix_type=matrix_type, &
                                 group=group_handle)

         IF (PRESENT(distribution)) distribution%dbcsr = my_distribution
         IF (PRESENT(group)) CALL group%set_handle(group_handle)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_get_info

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_get_matrix_type(matrix) RESULT(matrix_type)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      CHARACTER                                          :: matrix_type

      IF (USE_DBCSR_BACKEND) THEN
         matrix_type = dbcsr_get_matrix_type_prv(matrix%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END FUNCTION dbcsr_get_matrix_type

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_get_num_blocks(matrix) RESULT(num_blocks)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: num_blocks

      IF (USE_DBCSR_BACKEND) THEN
         num_blocks = dbcsr_get_num_blocks_prv(matrix%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END FUNCTION dbcsr_get_num_blocks

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_get_occupation(matrix) RESULT(occupation)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      REAL(KIND=dp)                                      :: occupation

      IF (USE_DBCSR_BACKEND) THEN
         occupation = dbcsr_get_occupation_prv(matrix%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END FUNCTION dbcsr_get_occupation

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param row ...
!> \param column ...
!> \param processor ...
! **************************************************************************************************
   SUBROUTINE dbcsr_get_stored_coordinates(matrix, row, column, processor)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER, INTENT(IN)                                :: row, column
      INTEGER, INTENT(OUT)                               :: processor

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_get_stored_coordinates_prv(matrix%dbcsr, row, column, processor)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_get_stored_coordinates

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_has_symmetry(matrix) RESULT(has_symmetry)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      LOGICAL                                            :: has_symmetry

      IF (USE_DBCSR_BACKEND) THEN
         has_symmetry = dbcsr_has_symmetry_prv(matrix%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END FUNCTION dbcsr_has_symmetry

! **************************************************************************************************
!> \brief ...
!> \param iterator ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_iterator_blocks_left(iterator) RESULT(blocks_left)
      TYPE(dbcsr_iterator_type), INTENT(IN)              :: iterator
      LOGICAL                                            :: blocks_left

      IF (USE_DBCSR_BACKEND) THEN
         blocks_left = dbcsr_iterator_blocks_left_prv(iterator%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END FUNCTION dbcsr_iterator_blocks_left

! **************************************************************************************************
!> \brief ...
!> \param iterator ...
!> \param row ...
!> \param column ...
!> \param block ...
!> \param block_number_argument_has_been_removed ...
!> \param row_size ...
!> \param col_size ...
!> \param row_offset ...
!> \param col_offset ...
! **************************************************************************************************
   SUBROUTINE dbcsr_iterator_next_block(iterator, row, column, block, &
                                        block_number_argument_has_been_removed, &
                                        row_size, col_size, &
                                        row_offset, col_offset)
      TYPE(dbcsr_iterator_type), INTENT(INOUT)           :: iterator
      INTEGER, INTENT(OUT), OPTIONAL                     :: row, column
      REAL(kind=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: block
      LOGICAL, OPTIONAL :: block_number_argument_has_been_removed
      INTEGER, INTENT(OUT), OPTIONAL                     :: row_size, col_size, row_offset, &
                                                            col_offset

      INTEGER                                            :: my_column, my_row
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: my_block

      CPASSERT(.NOT. PRESENT(block_number_argument_has_been_removed))

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_iterator_next_block_prv(iterator%dbcsr, row=my_row, column=my_column, &
                                            block=my_block, row_size=row_size, col_size=col_size, &
                                            row_offset=row_offset, col_offset=col_offset)
         IF (PRESENT(block)) block => my_block
         IF (PRESENT(row)) row = my_row
         IF (PRESENT(column)) column = my_column
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_iterator_next_block

! **************************************************************************************************
!> \brief ...
!> \param iterator ...
!> \param matrix ...
!> \param shared ...
!> \param dynamic ...
!> \param dynamic_byrows ...
! **************************************************************************************************
   SUBROUTINE dbcsr_iterator_start(iterator, matrix, shared, dynamic, dynamic_byrows)
      TYPE(dbcsr_iterator_type), INTENT(OUT)             :: iterator
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: shared, dynamic, dynamic_byrows

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_iterator_start_prv(iterator%dbcsr, matrix%dbcsr, shared, dynamic, dynamic_byrows)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_iterator_start

! **************************************************************************************************
!> \brief Like dbcsr_iterator_start() but with matrix being INTENT(IN).
!>        When invoking this routine, the caller promises not to modify the returned blocks.
!> \param iterator ...
!> \param matrix ...
!> \param shared ...
!> \param dynamic ...
!> \param dynamic_byrows ...
! **************************************************************************************************
   SUBROUTINE dbcsr_iterator_readonly_start(iterator, matrix, shared, dynamic, dynamic_byrows)
      TYPE(dbcsr_iterator_type), INTENT(OUT)             :: iterator
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: shared, dynamic, dynamic_byrows

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_iterator_start_prv(iterator%dbcsr, matrix%dbcsr, shared, dynamic, &
                                       dynamic_byrows, read_only=.TRUE.)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_iterator_readonly_start

! **************************************************************************************************
!> \brief ...
!> \param iterator ...
! **************************************************************************************************
   SUBROUTINE dbcsr_iterator_stop(iterator)
      TYPE(dbcsr_iterator_type), INTENT(INOUT)           :: iterator

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_iterator_stop_prv(iterator%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_iterator_stop

! **************************************************************************************************
!> \brief ...
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE dbcsr_mp_grid_setup(dist)
      TYPE(dbcsr_distribution_type), INTENT(INOUT)       :: dist

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_mp_grid_setup_prv(dist%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_mp_grid_setup

! **************************************************************************************************
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
!> \param flop ...
! **************************************************************************************************
   SUBROUTINE dbcsr_multiply(transa, transb, alpha, matrix_a, matrix_b, beta, &
                             matrix_c, first_row, last_row, &
                             first_column, last_column, first_k, last_k, &
                             retain_sparsity, filter_eps, flop)
      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb
      REAL(kind=dp), INTENT(IN)                          :: alpha
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a, matrix_b
      REAL(kind=dp), INTENT(IN)                          :: beta
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_c
      INTEGER, INTENT(IN), OPTIONAL                      :: first_row, last_row, first_column, &
                                                            last_column, first_k, last_k
      LOGICAL, INTENT(IN), OPTIONAL                      :: retain_sparsity
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: filter_eps
      INTEGER(int_8), INTENT(OUT), OPTIONAL              :: flop

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_multiply_prv(transa, transb, alpha, matrix_a%dbcsr, matrix_b%dbcsr, beta, &
                                 matrix_c%dbcsr, first_row, last_row, first_column, last_column, &
                                 first_k, last_k, retain_sparsity, filter_eps=filter_eps, flop=flop)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_multiply

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param summation ...
! **************************************************************************************************
   SUBROUTINE dbcsr_put_block(matrix, row, col, block, summation)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN)                                :: row, col
      REAL(kind=dp), DIMENSION(:, :), INTENT(IN)         :: block
      LOGICAL, INTENT(IN), OPTIONAL                      :: summation

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_put_block_prv(matrix%dbcsr, row, col, block, summation=summation)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_put_block

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_release(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_release_prv(matrix%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_release

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_replicate_all(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_replicate_all_prv(matrix%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_replicate_all

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param rows ...
!> \param cols ...
! **************************************************************************************************
   SUBROUTINE dbcsr_reserve_blocks(matrix, rows, cols)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, DIMENSION(:), INTENT(IN)                  :: rows, cols

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_reserve_blocks_prv(matrix%dbcsr, rows, cols)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_reserve_blocks

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param alpha_scalar ...
! **************************************************************************************************
   SUBROUTINE dbcsr_scale(matrix, alpha_scalar)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(kind=dp), INTENT(IN)                          :: alpha_scalar

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_scale_prv(matrix%dbcsr, alpha_scalar)
      ELSE
         CALL dbm_scale(matrix%dbm, alpha_scalar)
      END IF
   END SUBROUTINE dbcsr_scale

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param alpha ...
! **************************************************************************************************
   SUBROUTINE dbcsr_set(matrix, alpha)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(kind=dp), INTENT(IN)                          :: alpha

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_set_prv(matrix%dbcsr, alpha)
      ELSE
         IF (alpha == 0.0_dp) THEN
            CALL dbm_zero(matrix%dbm)
         ELSE
            CPABORT("Not yet implemented for DBM.")
         END IF
      END IF
   END SUBROUTINE dbcsr_set

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_sum_replicated(matrix)
      TYPE(dbcsr_type), INTENT(inout)                    :: matrix

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_sum_replicated_prv(matrix%dbcsr)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_sum_replicated

! **************************************************************************************************
!> \brief ...
!> \param transposed ...
!> \param normal ...
!> \param shallow_data_copy ...
!> \param transpose_distribution ...
!> \param use_distribution ...
! **************************************************************************************************
   SUBROUTINE dbcsr_transposed(transposed, normal, shallow_data_copy, transpose_distribution, &
                               use_distribution)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: transposed
      TYPE(dbcsr_type), INTENT(IN)                       :: normal
      LOGICAL, INTENT(IN), OPTIONAL                      :: shallow_data_copy, transpose_distribution
      TYPE(dbcsr_distribution_type), INTENT(IN), &
         OPTIONAL                                        :: use_distribution

      IF (USE_DBCSR_BACKEND) THEN
         IF (PRESENT(use_distribution)) THEN
            CALL dbcsr_transposed_prv(transposed%dbcsr, normal%dbcsr, &
                                      shallow_data_copy=shallow_data_copy, &
                                      transpose_distribution=transpose_distribution, &
                                      use_distribution=use_distribution%dbcsr)
         ELSE
            CALL dbcsr_transposed_prv(transposed%dbcsr, normal%dbcsr, &
                                      shallow_data_copy=shallow_data_copy, &
                                      transpose_distribution=transpose_distribution)
         END IF
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_transposed

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \return ...
! **************************************************************************************************
   FUNCTION dbcsr_valid_index(matrix) RESULT(valid_index)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      LOGICAL                                            :: valid_index

      IF (USE_DBCSR_BACKEND) THEN
         valid_index = dbcsr_valid_index_prv(matrix%dbcsr)
      ELSE
         valid_index = .TRUE. ! Does not apply to DBM.
      END IF
   END FUNCTION dbcsr_valid_index

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param verbosity ...
!> \param local ...
! **************************************************************************************************
   SUBROUTINE dbcsr_verify_matrix(matrix, verbosity, local)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER, INTENT(IN), OPTIONAL                      :: verbosity
      LOGICAL, INTENT(IN), OPTIONAL                      :: local

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_verify_matrix_prv(matrix%dbcsr, verbosity, local)
      ELSE
         ! Does not apply to DBM.
      END IF
   END SUBROUTINE dbcsr_verify_matrix

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param nblks_guess ...
!> \param sizedata_guess ...
!> \param n ...
!> \param work_mutable ...
! **************************************************************************************************
   SUBROUTINE dbcsr_work_create(matrix, nblks_guess, sizedata_guess, n, work_mutable)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN), OPTIONAL                      :: nblks_guess, sizedata_guess, n
      LOGICAL, INTENT(in), OPTIONAL                      :: work_mutable

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_work_create_prv(matrix%dbcsr, nblks_guess, sizedata_guess, n, work_mutable)
      ELSE
         ! Does not apply to DBM.
      END IF
   END SUBROUTINE dbcsr_work_create

! **************************************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param RESULT ...
! **************************************************************************************************
   SUBROUTINE dbcsr_dot_threadsafe(matrix_a, matrix_b, RESULT)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a, matrix_b
      REAL(kind=dp), INTENT(INOUT)                       :: result

      IF (USE_DBCSR_BACKEND) THEN
         CALL dbcsr_dot_prv(matrix_a%dbcsr, matrix_b%dbcsr, RESULT)
      ELSE
         CPABORT("Not yet implemented for DBM.")
      END IF
   END SUBROUTINE dbcsr_dot_threadsafe

END MODULE cp_dbcsr_api

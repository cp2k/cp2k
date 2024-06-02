!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE cp_dbcsr_api
   USE kinds, ONLY: dp, int_8
   USE dbcsr_api, ONLY: &
      dbcsr_add_prv => dbcsr_add, &
      dbcsr_add_block_node_prv => dbcsr_add_block_node, &
      dbcsr_add_on_diag_prv => dbcsr_add_on_diag, &
      dbcsr_binary_read_prv => dbcsr_binary_read, &
      dbcsr_binary_write_prv => dbcsr_binary_write, &
      dbcsr_checksum_prv => dbcsr_checksum, &
      dbcsr_clear_prv => dbcsr_clear, &
      dbcsr_clear_mempools, &
      dbcsr_complete_redistribute_prv => dbcsr_complete_redistribute, &
      convert_csr_to_dbcsr_prv => dbcsr_convert_csr_to_dbcsr, &
      convert_dbcsr_to_csr_prv => dbcsr_convert_dbcsr_to_csr, &
      dbcsr_convert_offsets_to_sizes, &
      dbcsr_convert_sizes_to_offsets, &
      dbcsr_copy_prv => dbcsr_copy, &
      dbcsr_create_prv => dbcsr_create, &
      dbcsr_csr_create, &
      dbcsr_csr_create_from_dbcsr_prv => dbcsr_csr_create_from_dbcsr, &
      dbcsr_csr_dbcsr_blkrow_dist, dbcsr_csr_destroy, dbcsr_csr_eqrow_floor_dist, &
      dbcsr_csr_p_type, dbcsr_csr_print_sparsity, dbcsr_csr_type, dbcsr_csr_write, &
      dbcsr_desymmetrize_prv => dbcsr_desymmetrize, &
      dbcsr_distribute_prv => dbcsr_distribute, &
      dbcsr_distribution_get_prv => dbcsr_distribution_get, &
      dbcsr_distribution_get_num_images, &
      dbcsr_distribution_hold_prv => dbcsr_distribution_hold, &
      dbcsr_distribution_new_prv => dbcsr_distribution_new, &
      dbcsr_distribution_release_prv => dbcsr_distribution_release, &
      dbcsr_distribution_type_prv => dbcsr_distribution_type, &
      dbcsr_dot_prv => dbcsr_dot, &
      dbcsr_filter_prv => dbcsr_filter, &
      dbcsr_finalize_prv => dbcsr_finalize, &
      dbcsr_finalize_lib, &
      dbcsr_frobenius_norm_prv => dbcsr_frobenius_norm, &
      dbcsr_func_dtanh, dbcsr_func_inverse, dbcsr_func_tanh, &
      dbcsr_function_of_elements_prv => dbcsr_function_of_elements, &
      dbcsr_gershgorin_norm_prv => dbcsr_gershgorin_norm, &
      dbcsr_get_block_diag_prv => dbcsr_get_block_diag, &
      dbcsr_get_block_p_prv => dbcsr_get_block_p, &
      dbcsr_get_data_p_prv => dbcsr_get_data_p, &
      dbcsr_get_data_size_prv => dbcsr_get_data_size, &
      dbcsr_get_data_type_prv => dbcsr_get_data_type, &
      dbcsr_get_default_config, &
      dbcsr_get_diag_prv => dbcsr_get_diag, &
      dbcsr_get_info_prv => dbcsr_get_info, &
      dbcsr_get_matrix_type_prv => dbcsr_get_matrix_type, &
      dbcsr_get_num_blocks_prv => dbcsr_get_num_blocks, &
      dbcsr_get_occupation_prv => dbcsr_get_occupation, &
      dbcsr_get_stored_coordinates_prv => dbcsr_get_stored_coordinates, &
      dbcsr_hadamard_product_prv => dbcsr_hadamard_product, &
      dbcsr_has_symmetry_prv => dbcsr_has_symmetry, &
      dbcsr_init_lib, &
      dbcsr_init_random_prv => dbcsr_init_random, &
      dbcsr_iterator_blocks_left_prv => dbcsr_iterator_blocks_left, &
      dbcsr_iterator_next_block_prv => dbcsr_iterator_next_block, &
      dbcsr_iterator_start_prv => dbcsr_iterator_start, &
      dbcsr_iterator_stop_prv => dbcsr_iterator_stop, &
      dbcsr_iterator_type_prv => dbcsr_iterator_type, &
      dbcsr_maxabs_prv => dbcsr_maxabs, &
      dbcsr_mp_grid_setup_prv => dbcsr_mp_grid_setup, &
      dbcsr_multiply_prv => dbcsr_multiply, &
      dbcsr_nblkcols_total_prv => dbcsr_nblkcols_total, &
      dbcsr_nblkrows_total_prv => dbcsr_nblkrows_total, &
      dbcsr_nfullcols_total_prv => dbcsr_nfullcols_total, &
      dbcsr_nfullrows_total_prv => dbcsr_nfullrows_total, &
      dbcsr_no_transpose, &
      dbcsr_norm_prv => dbcsr_norm, &
      dbcsr_norm_column, dbcsr_norm_frobenius, dbcsr_norm_maxabsnorm, &
      dbcsr_print_prv => dbcsr_print, &
      dbcsr_print_block_sum_prv => dbcsr_print_block_sum, &
      dbcsr_print_config, dbcsr_print_statistics, &
      dbcsr_put_block_prv => dbcsr_put_block, &
      dbcsr_release_prv => dbcsr_release, &
      dbcsr_replicate_all_prv => dbcsr_replicate_all, &
      dbcsr_reserve_all_blocks_prv => dbcsr_reserve_all_blocks, &
      dbcsr_reserve_block2d_prv => dbcsr_reserve_block2d, &
      dbcsr_reserve_blocks_prv => dbcsr_reserve_blocks, &
      dbcsr_reserve_diag_blocks_prv => dbcsr_reserve_diag_blocks, &
      dbcsr_reset_randmat_seed, dbcsr_run_tests, &
      dbcsr_scale_prv => dbcsr_scale, &
      dbcsr_scale_by_vector_prv => dbcsr_scale_by_vector, &
      dbcsr_set_prv => dbcsr_set, &
      dbcsr_set_config, &
      dbcsr_set_diag_prv => dbcsr_set_diag, &
      dbcsr_setname_prv => dbcsr_setname, &
      dbcsr_sum_replicated_prv => dbcsr_sum_replicated, &
      dbcsr_test_binary_io, dbcsr_test_mm, &
      dbcsr_trace_prv => dbcsr_trace, &
      dbcsr_transpose, &
      dbcsr_transposed_prv => dbcsr_transposed, &
      dbcsr_triu_prv => dbcsr_triu, &
      dbcsr_type_prv => dbcsr_type, &
      dbcsr_type_antisymmetric, dbcsr_type_complex_8, &
      dbcsr_type_complex_default, dbcsr_type_no_symmetry, dbcsr_type_real_8, &
      dbcsr_type_real_default, dbcsr_type_symmetric, &
      dbcsr_valid_index_prv => dbcsr_valid_index, &
      dbcsr_verify_matrix_prv => dbcsr_verify_matrix, &
      dbcsr_work_create_prv => dbcsr_work_create

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! constants
   PUBLIC :: dbcsr_type_no_symmetry
   PUBLIC :: dbcsr_type_symmetric
   PUBLIC :: dbcsr_type_antisymmetric
   PUBLIC :: dbcsr_transpose
   PUBLIC :: dbcsr_no_transpose
   PUBLIC :: dbcsr_type_complex_8
   PUBLIC :: dbcsr_type_real_8
   PUBLIC :: dbcsr_type_complex_default
   PUBLIC :: dbcsr_type_real_default

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
   PUBLIC :: dbcsr_scale_by_vector
   PUBLIC :: dbcsr_transposed
   PUBLIC :: dbcsr_multiply
   PUBLIC :: dbcsr_copy
   PUBLIC :: dbcsr_desymmetrize
   PUBLIC :: dbcsr_add_on_diag
   PUBLIC :: dbcsr_get_block_diag
   PUBLIC :: dbcsr_set_diag
   PUBLIC :: dbcsr_get_diag
   PUBLIC :: dbcsr_filter
   PUBLIC :: dbcsr_trace
   PUBLIC :: dbcsr_dot
   PUBLIC :: dbcsr_complete_redistribute
   PUBLIC :: dbcsr_get_block_p
   PUBLIC :: dbcsr_clear

   ! block reservation
   PUBLIC :: dbcsr_reserve_diag_blocks
   PUBLIC :: dbcsr_reserve_block2d
   PUBLIC :: dbcsr_reserve_blocks
   PUBLIC :: dbcsr_reserve_all_blocks

   ! iterator
   PUBLIC :: dbcsr_iterator_start
   PUBLIC :: dbcsr_iterator_stop
   PUBLIC :: dbcsr_iterator_blocks_left
   PUBLIC :: dbcsr_iterator_next_block

   ! getters / setters
   PUBLIC :: dbcsr_get_info
   PUBLIC :: dbcsr_distribution_get
   PUBLIC :: dbcsr_setname
   PUBLIC :: dbcsr_get_matrix_type
   PUBLIC :: dbcsr_get_occupation
   PUBLIC :: dbcsr_nblkrows_total
   PUBLIC :: dbcsr_nblkcols_total
   PUBLIC :: dbcsr_get_num_blocks
   PUBLIC :: dbcsr_get_data_size
   PUBLIC :: dbcsr_has_symmetry
   PUBLIC :: dbcsr_nfullrows_total
   PUBLIC :: dbcsr_nfullcols_total
   PUBLIC :: dbcsr_get_stored_coordinates
   PUBLIC :: dbcsr_valid_index
   PUBLIC :: dbcsr_get_data_type

   ! work operations
   PUBLIC :: dbcsr_add_block_node
   PUBLIC :: dbcsr_put_block
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

   ! high level matrix functions
   PUBLIC :: dbcsr_norm_frobenius
   PUBLIC :: dbcsr_norm_maxabsnorm
   PUBLIC :: dbcsr_norm_column
   PUBLIC :: dbcsr_hadamard_product
   PUBLIC :: dbcsr_func_dtanh
   PUBLIC :: dbcsr_func_inverse
   PUBLIC :: dbcsr_func_tanh
   PUBLIC :: dbcsr_print
   PUBLIC :: dbcsr_print_block_sum
   PUBLIC :: dbcsr_checksum
   PUBLIC :: dbcsr_maxabs
   PUBLIC :: dbcsr_norm
   PUBLIC :: dbcsr_gershgorin_norm
   PUBLIC :: dbcsr_frobenius_norm
   PUBLIC :: dbcsr_init_random
   PUBLIC :: dbcsr_function_of_elements
   PUBLIC :: dbcsr_triu

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

   ! binary io
   PUBLIC :: dbcsr_binary_write
   PUBLIC :: dbcsr_binary_read
   PUBLIC :: dbcsr_test_binary_io

   TYPE dbcsr_p_type
      TYPE(dbcsr_type), POINTER :: matrix => Null()
   END TYPE

   TYPE dbcsr_type
      TYPE(dbcsr_type_prv), PRIVATE :: prv
   END TYPE dbcsr_type

   TYPE dbcsr_distribution_type
      TYPE(dbcsr_distribution_type_prv), PRIVATE        :: prv
   END TYPE dbcsr_distribution_type

   TYPE dbcsr_iterator_type
      TYPE(dbcsr_iterator_type_prv), PRIVATE :: prv
   END TYPE dbcsr_iterator_type

   INTERFACE dbcsr_add
      MODULE PROCEDURE dbcsr_add_d, dbcsr_add_z
   END INTERFACE

   INTERFACE dbcsr_add_on_diag
      MODULE PROCEDURE dbcsr_add_on_diag_d, dbcsr_add_on_diag_z
   END INTERFACE

   INTERFACE dbcsr_create
      MODULE PROCEDURE dbcsr_create_new, dbcsr_create_template
   END INTERFACE

   INTERFACE dbcsr_dot
      MODULE PROCEDURE dbcsr_dot_d, dbcsr_dot_z
   END INTERFACE

   INTERFACE dbcsr_get_block_p
      MODULE PROCEDURE dbcsr_get_block_p_d, dbcsr_get_block_p_z
      MODULE PROCEDURE dbcsr_get_2d_block_p_d, dbcsr_get_2d_block_p_z
   END INTERFACE

   INTERFACE dbcsr_get_data_p
      MODULE PROCEDURE dbcsr_get_data_d, dbcsr_get_data_z
   END INTERFACE

   INTERFACE dbcsr_get_diag
      MODULE PROCEDURE dbcsr_get_diag_d, dbcsr_get_diag_z
   END INTERFACE

   INTERFACE dbcsr_iterator_next_block
      MODULE PROCEDURE dbcsr_iterator_next_block_index
      MODULE PROCEDURE dbcsr_iterator_next_1d_block_d, dbcsr_iterator_next_1d_block_z
      MODULE PROCEDURE dbcsr_iterator_next_2d_block_d, dbcsr_iterator_next_2d_block_z
   END INTERFACE

   INTERFACE dbcsr_multiply
      MODULE PROCEDURE dbcsr_multiply_d, dbcsr_multiply_z
   END INTERFACE

   INTERFACE dbcsr_norm
      MODULE PROCEDURE dbcsr_norm_scalar, dbcsr_norm_vector
   END INTERFACE

   INTERFACE dbcsr_put_block
      MODULE PROCEDURE dbcsr_put_block_d, dbcsr_put_block_z
      MODULE PROCEDURE dbcsr_put_block2d_d, dbcsr_put_block2d_z
   END INTERFACE

   INTERFACE dbcsr_reserve_block2d
      MODULE PROCEDURE dbcsr_reserve_block2d_d, dbcsr_reserve_block2d_z
   END INTERFACE

   INTERFACE dbcsr_scale
      MODULE PROCEDURE dbcsr_scale_d, dbcsr_scale_z
   END INTERFACE

   INTERFACE dbcsr_scale_by_vector
      MODULE PROCEDURE dbcsr_scale_by_vector_d, dbcsr_scale_by_vector_z
   END INTERFACE

   INTERFACE dbcsr_set
      MODULE PROCEDURE dbcsr_set_d, dbcsr_set_z
   END INTERFACE

   INTERFACE dbcsr_set_diag
      MODULE PROCEDURE dbcsr_set_diag_d, dbcsr_set_diag_z
   END INTERFACE

   INTERFACE dbcsr_trace
      MODULE PROCEDURE dbcsr_trace_d, dbcsr_trace_z
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief ...
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
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_add_${nametype1}$ (matrix_a, matrix_b, alpha_scalar, beta_scalar)
         TYPE(dbcsr_type), INTENT(INOUT)            :: matrix_a
         TYPE(dbcsr_type), INTENT(IN)               :: matrix_b
         ${type1}$, INTENT(IN)                      :: alpha_scalar, beta_scalar

         CALL dbcsr_add_prv(matrix_a%prv, matrix_b%prv, alpha_scalar, beta_scalar)
      END SUBROUTINE dbcsr_add_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_add_block_node(matrix, block_row, block_col, block)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN)                                :: block_row, block_col
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block

      CALL dbcsr_add_block_node_prv(matrix%prv, block_row, block_col, block)
   END SUBROUTINE dbcsr_add_block_node

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_add_on_diag_${nametype1}$ (matrix, alpha_scalar)
         TYPE(dbcsr_type), INTENT(INOUT)                      :: matrix
         ${type1}$, INTENT(IN)                                :: alpha_scalar

         CALL dbcsr_add_on_diag_prv(matrix%prv, alpha_scalar)
      END SUBROUTINE dbcsr_add_on_diag_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_binary_read(filepath, distribution, matrix_new)
      CHARACTER(len=*), INTENT(IN)                       :: filepath
      TYPE(dbcsr_distribution_type), INTENT(IN)          :: distribution
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_new

      CALL dbcsr_binary_read_prv(filepath, distribution%prv, matrix_new%prv)
   END SUBROUTINE dbcsr_binary_read

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_binary_write(matrix, filepath)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      CHARACTER(LEN=*), INTENT(IN)                       :: filepath

      CALL dbcsr_binary_write_prv(matrix%prv, filepath)
   END SUBROUTINE dbcsr_binary_write

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   FUNCTION dbcsr_checksum(matrix, pos) RESULT(checksum)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: pos
      REAL(KIND=dp)                                      :: checksum

      checksum = dbcsr_checksum_prv(matrix%prv, pos=pos)
   END FUNCTION dbcsr_checksum

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_clear(matrix)
      TYPE(dbcsr_type), INTENT(INOUT) :: matrix

      CALL dbcsr_clear_prv(matrix%prv)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_complete_redistribute(matrix, redist)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(dbcsr_type), INTENT(INOUT)                    :: redist

      CALL dbcsr_complete_redistribute_prv(matrix%prv, redist%prv)
   END SUBROUTINE dbcsr_complete_redistribute

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_convert_csr_to_dbcsr(dbcsr_mat, csr_mat)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: dbcsr_mat
      TYPE(dbcsr_csr_type), INTENT(INOUT)                :: csr_mat

      CALL convert_csr_to_dbcsr_prv(dbcsr_mat%prv, csr_mat)
   END SUBROUTINE dbcsr_convert_csr_to_dbcsr

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_convert_dbcsr_to_csr(dbcsr_mat, csr_mat)
      TYPE(dbcsr_type), INTENT(IN)                       :: dbcsr_mat
      TYPE(dbcsr_csr_type), INTENT(INOUT)                :: csr_mat

      CALL convert_dbcsr_to_csr_prv(dbcsr_mat%prv, csr_mat)
   END SUBROUTINE dbcsr_convert_dbcsr_to_csr

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_copy(matrix_b, matrix_a, name, keep_sparsity, keep_imaginary)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_b
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: name
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_sparsity, keep_imaginary

      CALL dbcsr_copy_prv(matrix_b%prv, matrix_a%prv, name=name, keep_sparsity=keep_sparsity, &
                          keep_imaginary=keep_imaginary)
   END SUBROUTINE dbcsr_copy

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_create_new(matrix, name, dist, matrix_type, row_blk_size, col_blk_size, nze, &
                               data_type, reuse_arrays, mutable_work)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      CHARACTER(len=*), INTENT(IN)                       :: name
      TYPE(dbcsr_distribution_type), INTENT(IN)          :: dist
      CHARACTER, INTENT(IN)                              :: matrix_type
      INTEGER, DIMENSION(:), INTENT(INOUT), POINTER      :: row_blk_size, col_blk_size
      INTEGER, INTENT(IN), OPTIONAL                      :: nze, data_type
      LOGICAL, INTENT(IN), OPTIONAL                      :: reuse_arrays, mutable_work

      CALL dbcsr_create_prv(matrix=matrix%prv, name=name, dist=dist%prv, matrix_type=matrix_type, &
                            row_blk_size=row_blk_size, col_blk_size=col_blk_size, nze=nze, &
                            data_type=data_type, reuse_arrays=reuse_arrays, &
                            mutable_work=mutable_work)
   END SUBROUTINE dbcsr_create_new

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_create_template(matrix, name, template, dist, matrix_type, &
                                    row_blk_size, col_blk_size, nze, data_type, &
                                    reuse_arrays, mutable_work)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: name
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(dbcsr_distribution_type), INTENT(IN), &
         OPTIONAL                                        :: dist
      CHARACTER, INTENT(IN), OPTIONAL                    :: matrix_type
      INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL, &
         POINTER                                         :: row_blk_size, col_blk_size
      INTEGER, INTENT(IN), OPTIONAL                      :: nze, data_type
      LOGICAL, INTENT(IN), OPTIONAL                      :: reuse_arrays, mutable_work

      CALL dbcsr_create_prv(matrix=matrix%prv, name=name, template=template%prv, dist=dist%prv, &
                            matrix_type=matrix_type, &
                            row_blk_size=row_blk_size, col_blk_size=col_blk_size, &
                            nze=nze, data_type=data_type, reuse_arrays=reuse_arrays, &
                            mutable_work=mutable_work)
   END SUBROUTINE dbcsr_create_template

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_csr_create_from_dbcsr(dbcsr_mat, csr_mat, dist_format, csr_sparsity, numnodes)

      TYPE(dbcsr_type), INTENT(IN)                       :: dbcsr_mat
      TYPE(dbcsr_csr_type), INTENT(OUT)                  :: csr_mat
      INTEGER                                            :: dist_format
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: csr_sparsity
      INTEGER, INTENT(IN), OPTIONAL                      :: numnodes

      IF (PRESENT(csr_sparsity)) THEN
         CALL dbcsr_csr_create_from_dbcsr_prv(dbcsr_mat%prv, csr_mat, dist_format, &
                                              csr_sparsity%prv, numnodes)
      ELSE
         CALL dbcsr_csr_create_from_dbcsr_prv(dbcsr_mat%prv, csr_mat, dist_format, numnodes=numnodes)
      END IF
   END SUBROUTINE dbcsr_csr_create_from_dbcsr

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_desymmetrize(matrix_a, matrix_b)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_b

      CALL dbcsr_desymmetrize_prv(matrix_a%prv, matrix_b%prv)
   END SUBROUTINE dbcsr_desymmetrize

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_distribute(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CALL dbcsr_distribute_prv(matrix%prv)
   END SUBROUTINE dbcsr_distribute

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_distribution_get(dist, row_dist, col_dist, nrows, ncols, has_threads, &
                                     group, mynode, numnodes, nprows, npcols, myprow, mypcol, &
                                     pgrid, subgroups_defined, prow_group, pcol_group)
      TYPE(dbcsr_distribution_type), INTENT(IN)          :: dist
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: row_dist, col_dist
      INTEGER, INTENT(OUT), OPTIONAL                     :: nrows, ncols
      LOGICAL, INTENT(OUT), OPTIONAL                     :: has_threads
      INTEGER, INTENT(OUT), OPTIONAL                     :: group, mynode, numnodes, &
                                                            nprows, npcols, myprow, mypcol
      INTEGER, DIMENSION(:, :), OPTIONAL, POINTER        :: pgrid
      LOGICAL, INTENT(OUT), OPTIONAL                     :: subgroups_defined
      INTEGER, INTENT(OUT), OPTIONAL                     :: prow_group, pcol_group

      call dbcsr_distribution_get_prv(dist%prv, row_dist, col_dist, nrows, ncols, has_threads, &
                                      group, mynode, numnodes, nprows, npcols, myprow, mypcol, &
                                      pgrid, subgroups_defined, prow_group, pcol_group)
   END SUBROUTINE dbcsr_distribution_get

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_distribution_hold(dist)
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL dbcsr_distribution_hold_prv(dist%prv)
   END SUBROUTINE dbcsr_distribution_hold

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_distribution_new(dist, template, group, pgrid, row_dist, col_dist, reuse_arrays)
      TYPE(dbcsr_distribution_type), INTENT(OUT)         :: dist
      TYPE(dbcsr_distribution_type), INTENT(IN), &
         OPTIONAL                                        :: template
      INTEGER, INTENT(IN), OPTIONAL                      :: group
      INTEGER, DIMENSION(:, :), OPTIONAL, POINTER        :: pgrid
      INTEGER, DIMENSION(:), INTENT(INOUT), POINTER      :: row_dist, col_dist
      LOGICAL, INTENT(IN), OPTIONAL                      :: reuse_arrays

      CALL dbcsr_distribution_new_prv(dist%prv, template%prv, group, pgrid, &
                                      row_dist, col_dist, reuse_arrays)
   END SUBROUTINE dbcsr_distribution_new

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_distribution_release(dist)
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL dbcsr_distribution_release_prv(dist%prv)
   END SUBROUTINE dbcsr_distribution_release

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_dot_${nametype1}$ (matrix_a, matrix_b, result)
         TYPE(dbcsr_type), INTENT(IN)              :: matrix_a, matrix_b
         ${type1}$, INTENT(INOUT)                  :: result

         CALL dbcsr_dot_prv(matrix_a%prv, matrix_b%prv, result)
      END SUBROUTINE dbcsr_dot_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_filter(matrix, eps)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(dp), INTENT(IN)                               :: eps

      CALL dbcsr_filter_prv(matrix%prv, eps)
   END SUBROUTINE dbcsr_filter

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_finalize(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CALL dbcsr_finalize_prv(matrix%prv)
   END SUBROUTINE dbcsr_finalize

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   FUNCTION dbcsr_frobenius_norm(matrix) RESULT(norm)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp)                                      :: norm

      norm = dbcsr_frobenius_norm_prv(matrix%prv)
   END FUNCTION dbcsr_frobenius_norm

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_function_of_elements(matrix, func, a0, a1, a2)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN)                                :: func
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: a0, a1, a2

      CALL dbcsr_function_of_elements_prv(matrix%prv, func, a0, a1, a2)
   END SUBROUTINE dbcsr_function_of_elements

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   FUNCTION dbcsr_gershgorin_norm(matrix) RESULT(norm)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp)                                      :: norm

      norm = dbcsr_gershgorin_norm_prv(matrix%prv)
   END FUNCTION dbcsr_gershgorin_norm

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_get_block_diag(matrix, diag)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(dbcsr_type), INTENT(INOUT)                    :: diag

      CALL dbcsr_get_block_diag_prv(matrix%prv, diag%prv)
   END SUBROUTINE dbcsr_get_block_diag

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_get_block_p_${nametype1}$ (matrix, row, col, block, found, row_size, col_size)
         TYPE(dbcsr_type), INTENT(IN)              :: matrix
         INTEGER, INTENT(IN)                       :: row, col
         ${type1}$, DIMENSION(:), POINTER          :: block
         LOGICAL, INTENT(OUT)                      :: found
         INTEGER, INTENT(OUT), OPTIONAL            :: row_size, col_size

         CALL dbcsr_get_block_p_prv(matrix%prv, row, col, block, found, row_size, col_size)
      END SUBROUTINE dbcsr_get_block_p_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_get_2d_block_p_${nametype1}$ (matrix, row, col, block, found, row_size, col_size)
         TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
         INTEGER, INTENT(IN)                      :: row, col
         ${type1}$, DIMENSION(:, :), POINTER      :: block
         LOGICAL, INTENT(OUT)                     :: found
         INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

         CALL dbcsr_get_block_p_prv(matrix%prv, row, col, block, found, row_size, col_size)
      END SUBROUTINE dbcsr_get_2d_block_p_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      FUNCTION dbcsr_get_data_${nametype1}$ (matrix, select_data_type, lb, ub) RESULT(res)
         TYPE(dbcsr_type), INTENT(IN)             :: matrix
         ${type1}$, INTENT(IN)                    :: select_data_type
         ${type1}$, DIMENSION(:), POINTER         :: res
         INTEGER, INTENT(IN), OPTIONAL            :: lb, ub

         res => dbcsr_get_data_p_prv(matrix%prv, select_data_type, lb, ub)
      END FUNCTION dbcsr_get_data_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   FUNCTION dbcsr_get_data_size(matrix) RESULT(data_size)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: data_size

      data_size = dbcsr_get_data_size_prv(matrix%prv)
   END FUNCTION dbcsr_get_data_size

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   PURE FUNCTION dbcsr_get_data_type(matrix) RESULT(data_type)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: data_type

      data_type = dbcsr_get_data_type_prv(matrix%prv)
   END FUNCTION dbcsr_get_data_type

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_get_diag_${nametype1}$ (matrix, diag)
         TYPE(dbcsr_type), INTENT(IN)                       :: matrix
         ${type1}$, DIMENSION(:), INTENT(OUT)                 :: diag

         CALL dbcsr_get_diag_prv(matrix%prv, diag)
      END SUBROUTINE dbcsr_get_diag_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_get_info(matrix, nblkrows_total, nblkcols_total, &
                             nfullrows_total, nfullcols_total, nblkrows_local, nblkcols_local, &
                             nfullrows_local, nfullcols_local, my_prow, my_pcol, &
                             local_rows, local_cols, proc_row_dist, proc_col_dist, &
                             row_blk_size, col_blk_size, row_blk_offset, col_blk_offset, &
                             distribution, name, matrix_type, data_type, group)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER, INTENT(OUT), OPTIONAL :: nblkrows_total, nblkcols_total, nfullrows_total, &
                                        nfullcols_total, nblkrows_local, nblkcols_local, &
                                        nfullrows_local, nfullcols_local, my_prow, my_pcol
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: local_rows, local_cols, &
                                                            proc_row_dist, proc_col_dist, &
                                                            row_blk_size, col_blk_size, &
                                                            row_blk_offset, col_blk_offset
      TYPE(dbcsr_distribution_type), INTENT(OUT), &
         OPTIONAL                                        :: distribution
      CHARACTER(len=*), INTENT(OUT), OPTIONAL            :: name
      CHARACTER, INTENT(OUT), OPTIONAL                   :: matrix_type
      INTEGER, INTENT(OUT), OPTIONAL                     :: data_type, group

      TYPE(dbcsr_distribution_type_prv)                  :: my_distribution

      CALL dbcsr_get_info_prv(matrix=matrix%prv, &
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
                              data_type=data_type, &
                              group=group)

      IF (PRESENT(distribution)) distribution%prv = my_distribution

   END SUBROUTINE dbcsr_get_info

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   PURE FUNCTION dbcsr_get_matrix_type(matrix) RESULT(matrix_type)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      CHARACTER                                          :: matrix_type

      matrix_type = dbcsr_get_matrix_type_prv(matrix%prv)
   END FUNCTION dbcsr_get_matrix_type

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   PURE FUNCTION dbcsr_get_num_blocks(matrix) RESULT(num_blocks)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: num_blocks

      num_blocks = dbcsr_get_num_blocks_prv(matrix%prv)
   END FUNCTION dbcsr_get_num_blocks

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   FUNCTION dbcsr_get_occupation(matrix) RESULT(occupation)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      REAL(KIND=dp)                                      :: occupation

      occupation = dbcsr_get_occupation_prv(matrix%prv)
   END FUNCTION dbcsr_get_occupation

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_get_stored_coordinates(matrix, row, column, processor)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER, INTENT(IN)                                :: row, column
      INTEGER, INTENT(OUT)                               :: processor

      CALL dbcsr_get_stored_coordinates_prv(matrix%prv, row, column, processor)
   END SUBROUTINE dbcsr_get_stored_coordinates

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_hadamard_product(matrix_a, matrix_b, matrix_c, b_assume_value)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a, matrix_b
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_c
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: b_assume_value

      CALL dbcsr_hadamard_product_prv(matrix_a%prv, matrix_b%prv, matrix_c%prv, b_assume_value)
   END SUBROUTINE dbcsr_hadamard_product

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   PURE FUNCTION dbcsr_has_symmetry(matrix) RESULT(has_symmetry)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      LOGICAL                                            :: has_symmetry

      has_symmetry = dbcsr_has_symmetry_prv(matrix%prv)
   END FUNCTION dbcsr_has_symmetry

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_init_random(matrix, keep_sparsity)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      LOGICAL, OPTIONAL                                  :: keep_sparsity

      CALL dbcsr_init_random_prv(matrix%prv, keep_sparsity=keep_sparsity)
   END SUBROUTINE dbcsr_init_random

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   PURE FUNCTION dbcsr_iterator_blocks_left(iterator) RESULT(blocks_left)
      TYPE(dbcsr_iterator_type), INTENT(IN)              :: iterator
      LOGICAL                                            :: blocks_left

      blocks_left = dbcsr_iterator_blocks_left_prv(iterator%prv)
   END FUNCTION dbcsr_iterator_blocks_left

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_iterator_next_block_index(iterator, row, column, blk)
      TYPE(dbcsr_iterator_type), INTENT(INOUT)           :: iterator
      INTEGER, INTENT(OUT)                               :: row, column, blk

      CALL dbcsr_iterator_next_block_prv(iterator%prv, row=row, column=column, blk=blk)
   END SUBROUTINE dbcsr_iterator_next_block_index

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_iterator_next_1d_block_${nametype1}$ (iterator, row, column, block, &
                                                             block_number, row_size, col_size, &
                                                             row_offset, col_offset)
         TYPE(dbcsr_iterator_type), INTENT(INOUT)        :: iterator
         INTEGER, INTENT(OUT)                            :: row, column
         ${type1}$, DIMENSION(:), POINTER                :: block
         INTEGER, INTENT(OUT), OPTIONAL                  :: block_number, row_size, col_size, &
                                                            row_offset, col_offset

         CALL dbcsr_iterator_next_block_prv(iterator%prv, row, column, block, block_number, &
                                            row_size, col_size, row_offset, col_offset)
      END SUBROUTINE dbcsr_iterator_next_1d_block_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_iterator_next_2d_block_${nametype1}$ (iterator, row, column, block, &
                                                             block_number, row_size, col_size, &
                                                             row_offset, col_offset)
         TYPE(dbcsr_iterator_type), INTENT(INOUT)        :: iterator
         INTEGER, INTENT(OUT)                            :: row, column
         ${type1}$, DIMENSION(:, :), POINTER             :: block
         INTEGER, INTENT(OUT), OPTIONAL                  :: block_number, row_size, col_size, &
                                                            row_offset, col_offset

         CALL dbcsr_iterator_next_block_prv(iterator%prv, row, column, block, block_number, &
                                            row_size, col_size, row_offset, col_offset)
      END SUBROUTINE dbcsr_iterator_next_2d_block_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_iterator_start(iterator, matrix, shared, dynamic, dynamic_byrows, read_only)
      TYPE(dbcsr_iterator_type), INTENT(OUT)             :: iterator
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: shared, dynamic, dynamic_byrows, &
                                                            read_only

      CALL dbcsr_iterator_start_prv(iterator%prv, matrix%prv, shared, dynamic, dynamic_byrows, &
                                    read_only)
   END SUBROUTINE dbcsr_iterator_start

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_iterator_stop(iterator)
      TYPE(dbcsr_iterator_type), INTENT(INOUT)           :: iterator

      CALL dbcsr_iterator_stop_prv(iterator%prv)
   END SUBROUTINE dbcsr_iterator_stop

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   FUNCTION dbcsr_maxabs(matrix) RESULT(norm)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp)                                      :: norm

      norm = dbcsr_maxabs_prv(matrix%prv)
   END FUNCTION dbcsr_maxabs

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_mp_grid_setup(dist)
      TYPE(dbcsr_distribution_type), INTENT(INOUT)       :: dist

      CALL dbcsr_mp_grid_setup_prv(dist%prv)
   END SUBROUTINE dbcsr_mp_grid_setup

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_multiply_${nametype1}$ (transa, transb, alpha, matrix_a, matrix_b, beta, &
                                               matrix_c, first_row, last_row, &
                                               first_column, last_column, first_k, last_k, &
                                               retain_sparsity, filter_eps, flop)
         CHARACTER(LEN=1), INTENT(IN)                    :: transa, transb
         ${type1}$, INTENT(IN)                           :: alpha
         TYPE(dbcsr_type), INTENT(IN)                    :: matrix_a, matrix_b
         ${type1}$, INTENT(IN)                           :: beta
         TYPE(dbcsr_type), INTENT(INOUT)                 :: matrix_c
         INTEGER, INTENT(IN), OPTIONAL                   :: first_row, last_row, &
                                                            first_column, last_column, &
                                                            first_k, last_k
         LOGICAL, INTENT(IN), OPTIONAL                   :: retain_sparsity
         REAL(kind=dp), INTENT(IN), OPTIONAL             :: filter_eps
         INTEGER(int_8), INTENT(OUT), OPTIONAL           :: flop

         CALL dbcsr_multiply_prv(transa, transb, alpha, matrix_a%prv, matrix_b%prv, beta, &
                                 matrix_c%prv, first_row, last_row, first_column, last_column, &
                                 first_k, last_k, retain_sparsity, filter_eps=filter_eps, flop=flop)
      END SUBROUTINE dbcsr_multiply_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   FUNCTION dbcsr_nblkcols_total(matrix) RESULT(nblkcols_total)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: nblkcols_total

      nblkcols_total = dbcsr_nblkcols_total_prv(matrix%prv)
   END FUNCTION dbcsr_nblkcols_total

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   FUNCTION dbcsr_nblkrows_total(matrix) RESULT(nblkrows_total)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: nblkrows_total

      nblkrows_total = dbcsr_nblkrows_total_prv(matrix%prv)
   END FUNCTION dbcsr_nblkrows_total

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   FUNCTION dbcsr_nfullcols_total(matrix) RESULT(nfullcols_total)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: nfullcols_total

      nfullcols_total = dbcsr_nfullcols_total_prv(matrix%prv)
   END FUNCTION dbcsr_nfullcols_total

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   FUNCTION dbcsr_nfullrows_total(matrix) RESULT(nfullrows_total)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: nfullrows_total

      nfullrows_total = dbcsr_nfullrows_total_prv(matrix%prv)
   END FUNCTION dbcsr_nfullrows_total

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_norm_scalar(matrix, which_norm, norm_scalar)
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: matrix
      INTEGER, INTENT(IN)                                :: which_norm
      REAL(dp), INTENT(OUT)                              :: norm_scalar

      CALL dbcsr_norm_prv(matrix%prv, which_norm, norm_scalar)
   END SUBROUTINE dbcsr_norm_scalar

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_norm_vector(matrix, which_norm, norm_vector)
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: matrix
      INTEGER, INTENT(IN)                                :: which_norm
      REAL(dp), DIMENSION(:), INTENT(OUT)                :: norm_vector

      CALL dbcsr_norm_prv(matrix%prv, which_norm, norm_vector)
   END SUBROUTINE dbcsr_norm_vector

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_print(matrix, nodata, matlab_format, variable_name, unit_nr)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: nodata, matlab_format
      CHARACTER(*), INTENT(in), OPTIONAL                 :: variable_name
      INTEGER, OPTIONAL                                  :: unit_nr

      CALL dbcsr_print_prv(matrix%prv, nodata, matlab_format, variable_name, unit_nr)
   END SUBROUTINE dbcsr_print

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_print_block_sum(matrix, unit_nr)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER, OPTIONAL                                  :: unit_nr

      CALL dbcsr_print_block_sum_prv(matrix%prv, unit_nr)
   END SUBROUTINE dbcsr_print_block_sum

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_put_block2d_${nametype1}$ (matrix, row, col, block, summation)
         TYPE(dbcsr_type), INTENT(INOUT)                 :: matrix
         INTEGER, INTENT(IN)                             :: row, col
         ${type1}$, DIMENSION(:, :), INTENT(IN)          :: block
         LOGICAL, INTENT(IN), OPTIONAL                   :: summation

         CALL dbcsr_put_block_prv(matrix%prv, row, col, block, summation=summation)
      END SUBROUTINE dbcsr_put_block2d_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_put_block_${nametype1}$ (matrix, row, col, block, summation)
         TYPE(dbcsr_type), INTENT(INOUT)                 :: matrix
         INTEGER, INTENT(IN)                             :: row, col
         ${type1}$, DIMENSION(:), INTENT(IN)             :: block
         LOGICAL, INTENT(IN), OPTIONAL                   :: summation

         CALL dbcsr_put_block_prv(matrix%prv, row, col, block, summation=summation)
      END SUBROUTINE dbcsr_put_block_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_release(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CALL dbcsr_release_prv(matrix%prv)
   END SUBROUTINE dbcsr_release

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_replicate_all(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CALL dbcsr_replicate_all_prv(matrix%prv)
   END SUBROUTINE dbcsr_replicate_all

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_reserve_all_blocks(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CALL dbcsr_reserve_all_blocks_prv(matrix%prv)
   END SUBROUTINE dbcsr_reserve_all_blocks

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_reserve_block2d_${nametype1}$ (matrix, row, col, block)
         TYPE(dbcsr_type), INTENT(INOUT)                 :: matrix
         INTEGER, INTENT(IN)                             :: row, col
         ${type1}$, DIMENSION(:, :), POINTER             :: block

         CALL dbcsr_reserve_block2d_prv(matrix%prv, row, col, block)
      END SUBROUTINE dbcsr_reserve_block2d_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_reserve_blocks(matrix, rows, cols)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, DIMENSION(:), INTENT(IN)                  :: rows, cols

      CALL dbcsr_reserve_blocks_prv(matrix%prv, rows, cols)
   END SUBROUTINE dbcsr_reserve_blocks

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_reserve_diag_blocks(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CALL dbcsr_reserve_diag_blocks_prv(matrix%prv)
   END SUBROUTINE dbcsr_reserve_diag_blocks

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_scale_${nametype1}$ (matrix, alpha_scalar)
         TYPE(dbcsr_type), INTENT(INOUT)                 :: matrix
         ${type1}$, INTENT(IN)                           :: alpha_scalar

         CALL dbcsr_scale_prv(matrix%prv, alpha_scalar)
      END SUBROUTINE dbcsr_scale_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_scale_by_vector_${nametype1}$ (matrix, alpha, side)
         TYPE(dbcsr_type), INTENT(INOUT)                 :: matrix
         ${type1}$, DIMENSION(:), INTENT(IN), TARGET     :: alpha
         CHARACTER(LEN=*), INTENT(IN)                    :: side

         CALL dbcsr_scale_by_vector_prv(matrix%prv, alpha, side)
      END SUBROUTINE dbcsr_scale_by_vector_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_set_${nametype1}$ (matrix, alpha)
         TYPE(dbcsr_type), INTENT(INOUT)                 :: matrix
         ${type1}$, INTENT(IN)                           :: alpha

         CALL dbcsr_set_prv(matrix%prv, alpha)
      END SUBROUTINE dbcsr_set_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_set_diag_${nametype1}$ (matrix, diag)
         TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
         ${type1}$, DIMENSION(:), INTENT(IN)                  :: diag

         CALL dbcsr_set_diag_prv(matrix%prv, diag)
      END SUBROUTINE dbcsr_set_diag_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_setname(matrix, newname)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      CHARACTER(len=*), INTENT(IN)                       :: newname

      CALL dbcsr_setname_prv(matrix%prv, newname)
   END SUBROUTINE dbcsr_setname

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_sum_replicated(matrix)
      TYPE(dbcsr_type), INTENT(inout)                    :: matrix

      CALL dbcsr_sum_replicated_prv(matrix%prv)
   END SUBROUTINE dbcsr_sum_replicated

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   #:for nametype1, type1 in [('d', 'REAL(kind=dp)'), ('z', 'COMPLEX(kind=dp)')]
      SUBROUTINE dbcsr_trace_${nametype1}$ (matrix, trace)
         TYPE(dbcsr_type), INTENT(IN)                    :: matrix
         ${type1}$, INTENT(OUT)                          :: trace

         CALL dbcsr_trace_prv(matrix%prv, trace)
      END SUBROUTINE dbcsr_trace_${nametype1}$
   #:endfor

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_transposed(transposed, normal, shallow_data_copy, transpose_distribution, &
                               use_distribution)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: transposed
      TYPE(dbcsr_type), INTENT(IN)                       :: normal
      LOGICAL, INTENT(IN), OPTIONAL                      :: shallow_data_copy, transpose_distribution
      TYPE(dbcsr_distribution_type), INTENT(IN), &
         OPTIONAL                                        :: use_distribution

      IF (PRESENT(use_distribution)) THEN
         CALL dbcsr_transposed_prv(transposed%prv, normal%prv, shallow_data_copy=shallow_data_copy, &
                                   transpose_distribution=transpose_distribution, &
                                   use_distribution=use_distribution%prv)
      ELSE
         CALL dbcsr_transposed_prv(transposed%prv, normal%prv, shallow_data_copy=shallow_data_copy, &
                                   transpose_distribution=transpose_distribution)
      END IF
   END SUBROUTINE dbcsr_transposed

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_triu(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CALL dbcsr_triu_prv(matrix%prv)
   END SUBROUTINE dbcsr_triu

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   PURE FUNCTION dbcsr_valid_index(matrix) RESULT(valid_index)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      LOGICAL                                            :: valid_index

      valid_index = dbcsr_valid_index_prv(matrix%prv)
   END FUNCTION dbcsr_valid_index

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_verify_matrix(matrix, verbosity, local)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER, INTENT(IN), OPTIONAL                      :: verbosity
      LOGICAL, INTENT(IN), OPTIONAL                      :: local

      CALL dbcsr_verify_matrix_prv(matrix%prv, verbosity, local)
   END SUBROUTINE dbcsr_verify_matrix

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE dbcsr_work_create(matrix, nblks_guess, sizedata_guess, n, work_mutable)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN), OPTIONAL                      :: nblks_guess, sizedata_guess, n
      LOGICAL, INTENT(in), OPTIONAL                      :: work_mutable

      CALL dbcsr_work_create_prv(matrix%prv, nblks_guess, sizedata_guess, n, work_mutable)
   END SUBROUTINE dbcsr_work_create

END MODULE cp_dbcsr_api

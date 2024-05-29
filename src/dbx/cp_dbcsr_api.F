!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE cp_dbcsr_api
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_block_node, dbcsr_add_on_diag, dbcsr_add_work_coordinate, &
        dbcsr_binary_read, dbcsr_binary_write, dbcsr_checksum, dbcsr_clear, dbcsr_clear_mempools, &
        dbcsr_complete_redistribute, dbcsr_convert_csr_to_dbcsr, dbcsr_convert_dbcsr_to_csr, &
        dbcsr_convert_offsets_to_sizes, dbcsr_convert_sizes_to_offsets, dbcsr_copy, &
        dbcsr_copy_into_existing, dbcsr_create, dbcsr_csr_create, dbcsr_csr_create_from_dbcsr, &
        dbcsr_csr_dbcsr_blkrow_dist, dbcsr_csr_destroy, dbcsr_csr_eqrow_ceil_dist, &
        dbcsr_csr_eqrow_floor_dist, dbcsr_csr_p_type, dbcsr_csr_print_sparsity, dbcsr_csr_type, &
        dbcsr_csr_write, dbcsr_deallocate_matrix, dbcsr_desymmetrize, dbcsr_distribute, &
        dbcsr_distribution_get, dbcsr_distribution_get_num_images, dbcsr_distribution_hold, &
        dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_dot, &
        dbcsr_filter, dbcsr_finalize, dbcsr_finalize_lib, dbcsr_frobenius_norm, dbcsr_func_artanh, &
        dbcsr_func_dtanh, dbcsr_func_inverse, dbcsr_func_tanh, dbcsr_function_of_elements, &
        dbcsr_gershgorin_norm, dbcsr_get_block_diag, dbcsr_get_block_p, dbcsr_get_data_p, &
        dbcsr_get_data_size, dbcsr_get_data_type, dbcsr_get_default_config, dbcsr_get_diag, &
        dbcsr_get_info, dbcsr_get_matrix_type, dbcsr_get_num_blocks, dbcsr_get_occupation, &
        dbcsr_get_stored_coordinates, dbcsr_get_wms_data_p, dbcsr_hadamard_product, &
        dbcsr_has_symmetry, dbcsr_init_lib, dbcsr_init_p, dbcsr_init_random, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_maxabs, dbcsr_mp_grid_setup, &
        dbcsr_multiply, dbcsr_nblkcols_local, dbcsr_nblkcols_total, dbcsr_nblkrows_local, &
        dbcsr_nblkrows_total, dbcsr_nfullcols_total, dbcsr_nfullrows_total, dbcsr_no_transpose, &
        dbcsr_norm, dbcsr_norm_column, dbcsr_norm_frobenius, dbcsr_norm_maxabsnorm, dbcsr_p_type, &
        dbcsr_print, dbcsr_print_block_sum, dbcsr_print_config, dbcsr_print_statistics, &
        dbcsr_put_block, dbcsr_release, dbcsr_release_p, dbcsr_replicate_all, &
        dbcsr_reserve_all_blocks, dbcsr_reserve_block2d, dbcsr_reserve_blocks, &
        dbcsr_reserve_diag_blocks, dbcsr_reset_randmat_seed, dbcsr_run_tests, dbcsr_scalar, &
        dbcsr_scalar_type, dbcsr_scale, dbcsr_scale_by_vector, dbcsr_set, dbcsr_set_config, &
        dbcsr_set_diag, dbcsr_set_work_size, dbcsr_setname, dbcsr_sum_replicated, &
        dbcsr_test_binary_io, dbcsr_test_mm, dbcsr_to_csr_filter, dbcsr_trace, dbcsr_transpose, &
        dbcsr_transposed, dbcsr_triu, dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_complex_4, &
        dbcsr_type_complex_8, dbcsr_type_complex_default, dbcsr_type_no_symmetry, &
        dbcsr_type_real_4, dbcsr_type_real_8, dbcsr_type_real_default, dbcsr_type_symmetric, &
        dbcsr_valid_index, dbcsr_verify_matrix, dbcsr_work_create

   IMPLICIT NONE

   ! constants
   PUBLIC :: dbcsr_type_no_symmetry
   PUBLIC :: dbcsr_type_symmetric
   PUBLIC :: dbcsr_type_antisymmetric
   PUBLIC :: dbcsr_transpose
   PUBLIC :: dbcsr_no_transpose
   PUBLIC :: dbcsr_type_complex_8
   PUBLIC :: dbcsr_type_real_4
   PUBLIC :: dbcsr_type_real_8
   PUBLIC :: dbcsr_type_complex_4
   PUBLIC :: dbcsr_type_complex_default
   PUBLIC :: dbcsr_type_real_default

   ! types
   PUBLIC :: dbcsr_type
   PUBLIC :: dbcsr_p_type
   PUBLIC :: dbcsr_distribution_type
   PUBLIC :: dbcsr_iterator_type
   PUBLIC :: dbcsr_scalar_type

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
   PUBLIC :: dbcsr_copy_into_existing
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
   PUBLIC :: dbcsr_nblkrows_local
   PUBLIC :: dbcsr_nblkcols_local
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
   PUBLIC :: dbcsr_add_work_coordinate
   PUBLIC :: dbcsr_get_wms_data_p
   PUBLIC :: dbcsr_get_data_p
   PUBLIC :: dbcsr_set_work_size
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
   PUBLIC :: dbcsr_scalar

   ! high level matrix functions
   PUBLIC :: dbcsr_norm_frobenius
   PUBLIC :: dbcsr_norm_maxabsnorm
   PUBLIC :: dbcsr_norm_column
   PUBLIC :: dbcsr_hadamard_product
   PUBLIC :: dbcsr_func_artanh
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
   PUBLIC :: dbcsr_csr_eqrow_ceil_dist
   PUBLIC :: dbcsr_csr_dbcsr_blkrow_dist
   PUBLIC :: dbcsr_csr_print_sparsity
   PUBLIC :: dbcsr_to_csr_filter
   PUBLIC :: dbcsr_csr_write

   ! binary io
   PUBLIC :: dbcsr_binary_write
   PUBLIC :: dbcsr_binary_read
   PUBLIC :: dbcsr_test_binary_io

END MODULE cp_dbcsr_api

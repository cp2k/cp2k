!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   This is the start of a dbcsr_api, all publically needed functions
!>          are exported here. The others remain private to the library.
!>          Currently, this is the CP2K used set.
!>          Ultimately, a reduced subset and well defined api will remain,
!>          possibly grouped in to standard and expert api.
!>          Currently, this is work in progress.
!> \author  Joost VandeVondele
! *****************************************************************************
MODULE dbcsr_api
  USE array_types,                     ONLY: array_data,&
                                             array_hold,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release,&
                                             array_size
  USE dbcsr_arnoldi,                   ONLY: dbcsr_arnoldi_ev
  USE dbcsr_arnoldi_data_methods,      ONLY: dbcsr_arnoldi_data,&
                                             deallocate_arnoldi_data,&
                                             get_selected_ritz_val,&
                                             get_selected_ritz_vector,&
                                             set_arnoldi_initial_vector,&
                                             setup_arnoldi_data
  USE dbcsr_block_access,              ONLY: dbcsr_get_block_p,&
                                             dbcsr_put_block,&
                                             dbcsr_reserve_all_blocks,&
                                             dbcsr_reserve_block2d,&
                                             dbcsr_reserve_blocks,&
                                             dbcsr_reserve_diag_blocks
  USE dbcsr_config,                    ONLY: &
       accdrv_avoid_after_busy, accdrv_binning_binsize, accdrv_binning_nbins, &
       accdrv_do_inhomogenous, accdrv_min_flop_process, accdrv_min_flop_sort, &
       accdrv_posterior_buffers, accdrv_posterior_streams, &
       accdrv_priority_buffers, accdrv_priority_streams, &
       dbcsr_get_conf_combtypes, dbcsr_get_conf_comm_thread_load, &
       dbcsr_get_conf_max_ele_block, dbcsr_get_conf_max_elements_images, &
       dbcsr_get_conf_mm_driver, dbcsr_get_conf_mm_stacksize, &
       dbcsr_get_conf_mpi_mem, dbcsr_get_conf_nstacks, &
       dbcsr_get_conf_subcomm, dbcsr_get_conf_use_comm_thread, &
       dbcsr_set_conf_combtypes, dbcsr_set_conf_comm_thread_load, &
       dbcsr_set_conf_external_timing, dbcsr_set_conf_max_ele_block, &
       dbcsr_set_conf_max_elements_images, dbcsr_set_conf_mm_driver, &
       dbcsr_set_conf_mm_stacksize, dbcsr_set_conf_mpi_mem, &
       dbcsr_set_conf_nstacks, dbcsr_set_conf_subcomm, &
       dbcsr_set_conf_use_comm_thread, dbcsr_set_default_config, has_acc, &
       has_mpi, mm_driver_acc, mm_driver_blas, mm_driver_matmul, &
       mm_driver_smm, mm_name_acc, mm_name_blas, mm_name_matmul, mm_name_smm, &
       multrec_limit
  USE dbcsr_csr_conversions,           ONLY: &
       convert_csr_to_dbcsr, convert_dbcsr_to_csr, csr_create, &
       csr_create_from_dbcsr, csr_create_new, csr_create_template, &
       csr_dbcsr_blkrow_dist, csr_destroy, csr_eqrow_ceil_dist, &
       csr_eqrow_floor_dist, csr_type, csr_write
  USE dbcsr_data_methods,              ONLY: &
       dbcsr_data_clear_pointer, dbcsr_data_init, dbcsr_data_new, &
       dbcsr_data_release, dbcsr_get_data_p, dbcsr_scalar, &
       dbcsr_scalar_fill_all, dbcsr_scalar_get_type, dbcsr_scalar_get_value, &
       dbcsr_scalar_set_type, dbcsr_scalar_zero
  USE dbcsr_dist_methods,              ONLY: &
       dbcsr_distribution_col_dist, dbcsr_distribution_has_threads, &
       dbcsr_distribution_init, dbcsr_distribution_local_cols, &
       dbcsr_distribution_local_rows, dbcsr_distribution_mp, &
       dbcsr_distribution_ncols, dbcsr_distribution_new, &
       dbcsr_distribution_nlocal_cols, dbcsr_distribution_nlocal_rows, &
       dbcsr_distribution_nrows, dbcsr_distribution_release, &
       dbcsr_distribution_row_dist
  USE dbcsr_dist_operations,           ONLY: create_bl_distribution,&
                                             dbcsr_create_dist_block_cyclic,&
                                             dbcsr_create_dist_r_unrot,&
                                             dbcsr_get_stored_coordinates
  USE dbcsr_error_handling,            ONLY: dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type
  USE dbcsr_io,                        ONLY: dbcsr_binary_read,&
                                             dbcsr_binary_write,&
                                             dbcsr_print,&
                                             dbcsr_print_block_sum
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_lib,                       ONLY: dbcsr_clear_mempools,&
                                             dbcsr_finalize_lib,&
                                             dbcsr_init_lib
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_sizes, dbcsr_distribution, dbcsr_get_data_size, &
       dbcsr_get_data_type, dbcsr_get_matrix_type, dbcsr_get_num_blocks, &
       dbcsr_get_nze, dbcsr_has_symmetry, dbcsr_init, dbcsr_max_col_size, &
       dbcsr_max_row_size, dbcsr_name, dbcsr_nblkcols_local, &
       dbcsr_nblkcols_total, dbcsr_nblkrows_local, dbcsr_nblkrows_total, &
       dbcsr_nfullcols_total, dbcsr_nfullrows_total, dbcsr_release, &
       dbcsr_row_block_sizes, dbcsr_setname, dbcsr_valid_index, &
       dbcsr_wm_use_mutable
  USE dbcsr_mp_methods,                ONLY: dbcsr_mp_group,&
                                             dbcsr_mp_hold,&
                                             dbcsr_mp_mynode,&
                                             dbcsr_mp_new,&
                                             dbcsr_mp_npcols,&
                                             dbcsr_mp_nprows,&
                                             dbcsr_mp_numnodes,&
                                             dbcsr_mp_release
  USE dbcsr_multiply_api,              ONLY: dbcsr_multiply
  USE dbcsr_operations,                ONLY: &
       dbcsr_add, dbcsr_add_on_diag, dbcsr_copy, dbcsr_copy_into_existing, &
       dbcsr_filter, dbcsr_frobenius_norm, dbcsr_function_of_elements, &
       dbcsr_gershgorin_norm, dbcsr_get_block_diag, dbcsr_get_diag, &
       dbcsr_get_info, dbcsr_get_occupation, dbcsr_hadamard_product, &
       dbcsr_init_random, dbcsr_maxabs, dbcsr_norm, dbcsr_scale, &
       dbcsr_scale_by_vector, dbcsr_set, dbcsr_set_diag, &
       dbcsr_sum_replicated, dbcsr_trace, dbcsr_triu
  USE dbcsr_tests,                     ONLY: dbcsr_run_tests,&
                                             dbcsr_test_arnoldi,&
                                             dbcsr_test_binary_io,&
                                             dbcsr_test_mm
  USE dbcsr_toollib,                   ONLY: swap
  USE dbcsr_transformations,           ONLY: dbcsr_complete_redistribute,&
                                             dbcsr_desymmetrize_deep,&
                                             dbcsr_distribute,&
                                             dbcsr_new_transposed,&
                                             dbcsr_replicate,&
                                             dbcsr_replicate_all
  USE dbcsr_types,                     ONLY: &
       dbcsr_data_obj, dbcsr_distribution_obj, dbcsr_func_artanh, &
       dbcsr_func_dtanh, dbcsr_func_inverse, dbcsr_func_tanh, dbcsr_iterator, &
       dbcsr_mp_obj, dbcsr_no_transpose, dbcsr_norm_column, &
       dbcsr_norm_frobenius, dbcsr_norm_maxabsnorm, dbcsr_obj, &
       dbcsr_obj_type_p, dbcsr_repl_full, dbcsr_scalar_type, &
       dbcsr_type_antisymmetric, dbcsr_type_complex_4, dbcsr_type_complex_8, &
       dbcsr_type_complex_default, dbcsr_type_no_symmetry, dbcsr_type_real_4, &
       dbcsr_type_real_8, dbcsr_type_real_default, dbcsr_type_symmetric, &
       dbcsr_work_type
  USE dbcsr_util,                      ONLY: convert_offsets_to_sizes,&
                                             convert_sizes_to_offsets,&
                                             dbcsr_checksum,&
                                             dbcsr_verify_matrix
  USE dbcsr_vector_operations,         ONLY: create_col_vec_from_matrix,&
                                             create_replicated_col_vec_from_matrix,&
                                             create_replicated_row_vec_from_matrix,&
                                             create_row_vec_from_matrix
  USE dbcsr_work_operations,           ONLY: add_work_coordinate,&
                                             dbcsr_create,&
                                             dbcsr_finalize,&
                                             dbcsr_work_create
  USE min_heap,                        ONLY: heap_fill,&
                                             heap_get_first,&
                                             heap_new,&
                                             heap_release,&
                                             heap_reset_first,&
                                             heap_t

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: heap_fill,&
            heap_get_first,&
            heap_new,&
            heap_release,&
            heap_reset_first,&
            heap_t,&
            array_data,&
            array_hold,&
            array_i1d_obj,&
            array_new,&
            array_release,&
            array_size,&
            array_nullify,&
            add_work_coordinate,&
            convert_offsets_to_sizes,&
            convert_sizes_to_offsets,&
            dbcsr_run_tests,&
            dbcsr_test_mm,&
            dbcsr_test_arnoldi,&
            dbcsr_test_binary_io,&
            create_bl_distribution,&
            accdrv_avoid_after_busy,&
            accdrv_binning_binsize,&
            accdrv_binning_nbins,&
            accdrv_do_inhomogenous,&
            accdrv_min_flop_process,&
            accdrv_min_flop_sort,&
            accdrv_posterior_buffers,&
            accdrv_posterior_streams,&
            accdrv_priority_buffers,&
            accdrv_priority_streams,&
            dbcsr_add,&
            dbcsr_add_on_diag,&
            dbcsr_arnoldi_data,&
            setup_arnoldi_data,&
            set_arnoldi_initial_vector,&
            deallocate_arnoldi_data,&
            dbcsr_arnoldi_ev,&
            get_selected_ritz_val,&
            get_selected_ritz_vector,&
            dbcsr_binary_read,&
            dbcsr_binary_write,&
            dbcsr_triu,&
            dbcsr_checksum,&
            dbcsr_clear_mempools,&
            dbcsr_complete_redistribute,&
            dbcsr_copy,&
            dbcsr_copy_into_existing,&
            dbcsr_create,&
            dbcsr_create_dist_block_cyclic,&
            dbcsr_create_dist_r_unrot,&
            dbcsr_data_clear_pointer,&
            dbcsr_data_init,&
            dbcsr_data_new,&
            dbcsr_data_obj,&
            dbcsr_data_release,&
            dbcsr_desymmetrize_deep,&
            dbcsr_distribute,&
            dbcsr_distribution,&
            dbcsr_distribution_col_dist,&
            dbcsr_distribution_has_threads,&
            dbcsr_distribution_init,&
            dbcsr_distribution_local_cols,&
            dbcsr_distribution_local_rows,&
            dbcsr_distribution_mp,&
            dbcsr_distribution_ncols,&
            dbcsr_distribution_new,&
            dbcsr_distribution_nlocal_cols,&
            dbcsr_distribution_nlocal_rows,&
            dbcsr_distribution_nrows,&
            dbcsr_distribution_obj,&
            dbcsr_distribution_release,&
            dbcsr_distribution_row_dist,&
            dbcsr_error_set,&
            dbcsr_error_stop,&
            dbcsr_error_type,&
            dbcsr_filter,&
            dbcsr_finalize,&
            dbcsr_finalize_lib,&
            dbcsr_frobenius_norm,&
            dbcsr_func_artanh,&
            dbcsr_func_dtanh,&
            dbcsr_func_inverse,&
            dbcsr_func_tanh,&
            dbcsr_function_of_elements,&
            dbcsr_gershgorin_norm,&
            dbcsr_get_block_diag,&
            dbcsr_get_block_p,&
            dbcsr_get_conf_combtypes,&
            dbcsr_get_conf_comm_thread_load,&
            dbcsr_get_conf_max_ele_block,&
            dbcsr_get_conf_mm_driver,&
            dbcsr_get_conf_max_elements_images,&
            dbcsr_get_conf_mm_stacksize,&
            dbcsr_get_conf_mpi_mem,&
            dbcsr_get_conf_nstacks,&
            dbcsr_get_conf_subcomm,&
            dbcsr_get_conf_use_comm_thread,&
            dbcsr_get_data_p,&
            dbcsr_get_data_size,&
            dbcsr_get_data_type,&
            dbcsr_get_diag,&
            dbcsr_get_info,&
            dbcsr_get_matrix_type,&
            dbcsr_get_num_blocks,&
            dbcsr_get_occupation,&
            dbcsr_get_stored_coordinates,&
            dbcsr_hadamard_product,&
            dbcsr_init,&
            dbcsr_init_lib,&
            dbcsr_init_random,&
            dbcsr_iterator,&
            dbcsr_iterator_blocks_left,&
            dbcsr_iterator_next_block,&
            dbcsr_iterator_start,&
            dbcsr_iterator_stop,&
            dbcsr_maxabs,&
            dbcsr_max_col_size,&
            dbcsr_max_row_size,&
            dbcsr_mp_group,&
            dbcsr_mp_hold,&
            dbcsr_mp_mynode,&
            dbcsr_mp_new,&
            dbcsr_mp_npcols,&
            dbcsr_mp_nprows,&
            dbcsr_mp_numnodes,&
            dbcsr_mp_obj,&
            dbcsr_mp_release,&
            dbcsr_multiply,&
            dbcsr_name,&
            dbcsr_setname,&
            dbcsr_nblkcols_local,&
            dbcsr_nblkcols_total,&
            dbcsr_nblkrows_local,&
            dbcsr_nblkrows_total,&
            dbcsr_new_transposed,&
            dbcsr_nfullcols_total,&
            dbcsr_nfullrows_total,&
            dbcsr_norm,&
            dbcsr_norm_column,&
            dbcsr_norm_frobenius,&
            dbcsr_norm_maxabsnorm,&
            dbcsr_no_transpose,&
            dbcsr_obj,&
            dbcsr_obj_type_p,&
            dbcsr_print,&
            dbcsr_print_block_sum,&
            dbcsr_put_block,&
            dbcsr_release,&
            dbcsr_repl_full,&
            dbcsr_replicate,&
            dbcsr_replicate_all,&
            dbcsr_reserve_all_blocks,&
            dbcsr_reserve_block2d,&
            dbcsr_reserve_blocks,&
            dbcsr_reserve_diag_blocks,&
            dbcsr_scalar,&
            dbcsr_scalar_fill_all,&
            dbcsr_scalar_get_type,&
            dbcsr_scalar_get_value,&
            dbcsr_scalar_set_type,&
            dbcsr_scalar_type,&
            dbcsr_work_type,&
            dbcsr_scalar_zero,&
            dbcsr_scale,&
            dbcsr_scale_by_vector,&
            dbcsr_set,&
            dbcsr_set_conf_combtypes,&
            dbcsr_set_conf_comm_thread_load,&
            dbcsr_set_conf_external_timing,&
            dbcsr_set_conf_max_ele_block,&
            dbcsr_set_conf_mm_driver,&
            dbcsr_set_conf_mm_stacksize,&
            dbcsr_set_conf_mpi_mem,&
            dbcsr_set_conf_max_elements_images,&
            dbcsr_set_conf_nstacks,&
            dbcsr_set_conf_subcomm,&
            dbcsr_set_conf_use_comm_thread,&
            dbcsr_set_diag,&
            dbcsr_sum_replicated,&
            dbcsr_trace,&
            dbcsr_type_antisymmetric,&
            dbcsr_type_complex_4,&
            dbcsr_type_complex_8,&
            dbcsr_type_complex_default,&
            dbcsr_type_no_symmetry,&
            dbcsr_type_real_4,&
            dbcsr_type_real_8,&
            dbcsr_type_real_default,&
            dbcsr_type_symmetric,&
            dbcsr_valid_index,&
            dbcsr_verify_matrix,&
            dbcsr_wm_use_mutable,&
            dbcsr_work_create,&
            create_col_vec_from_matrix,&
            create_row_vec_from_matrix,&
            create_replicated_col_vec_from_matrix,&
            create_replicated_row_vec_from_matrix,&
            has_acc,&
            has_mpi,&
            mm_driver_blas,&
            mm_driver_acc,&
            mm_driver_matmul,&
            mm_driver_smm,&
            mm_name_blas,&
            mm_name_acc,&
            mm_name_matmul,&
            mm_name_smm,&
            multrec_limit,&
            swap,&
            dbcsr_set_default_config,&
            dbcsr_row_block_sizes,&
            dbcsr_col_block_sizes,&
            dbcsr_has_symmetry,&
            convert_csr_to_dbcsr,&
            convert_dbcsr_to_csr,&
            csr_create_from_dbcsr,&
            csr_destroy,&
            csr_type,&
            csr_create,&
            csr_create_new,&
            csr_create_template,&
            csr_write,&
            csr_dbcsr_blkrow_dist,& 
            csr_eqrow_ceil_dist,&
            csr_eqrow_floor_dist,&
            dbcsr_get_nze
END MODULE dbcsr_api

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2014  CP2K developers group   !
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
  USE dbcsr_block_access,              ONLY: dbcsr_access_flush,&
                                             dbcsr_get_block,&
                                             dbcsr_get_block_p,&
                                             dbcsr_put_block,&
                                             dbcsr_reserve_all_blocks,&
                                             dbcsr_reserve_block2d,&
                                             dbcsr_reserve_blocks,&
                                             dbcsr_reserve_diag_blocks
  USE dbcsr_block_operations,          ONLY: block_add_on_diag,&
                                             block_chol_inv,&
                                             block_set
  USE dbcsr_config,                    ONLY: &
       cudadrv_avoid_after_busy, cudadrv_binning_binsize, &
       cudadrv_binning_nbins, cudadrv_do_inhomogenous, &
       cudadrv_min_flop_process, cudadrv_min_flop_sort, &
       cudadrv_posterior_buffers, cudadrv_posterior_streams, &
       cudadrv_priority_buffers, cudadrv_priority_streams, &
       dbcsr_get_conf_combtypes, dbcsr_get_conf_comm_thread_load, &
       dbcsr_get_conf_max_ele_block, dbcsr_get_conf_mm_driver, &
       dbcsr_get_conf_mm_stacksize, dbcsr_get_conf_mpi_mem, &
       dbcsr_get_conf_nstacks, dbcsr_get_conf_subcomm, &
       dbcsr_get_conf_use_comm_thread, dbcsr_set_conf_combtypes, &
       dbcsr_set_conf_comm_thread_load, dbcsr_set_conf_external_timing, &
       dbcsr_set_conf_max_ele_block, dbcsr_set_conf_mm_driver, &
       dbcsr_set_conf_mm_stacksize, dbcsr_set_conf_mpi_mem, &
       dbcsr_set_conf_nstacks, dbcsr_set_conf_subcomm, &
       dbcsr_set_conf_use_comm_thread, has_cuda, has_mpi, mm_driver_blas, &
       mm_driver_cuda, mm_driver_matmul, mm_driver_plasma, mm_driver_smm, &
       mm_name_blas, mm_name_cuda, mm_name_matmul, mm_name_plasma, &
       mm_name_smm, multrec_limit
  USE dbcsr_data_methods,              ONLY: &
       dbcsr_data_clear_pointer, dbcsr_data_init, dbcsr_data_new, &
       dbcsr_data_release, dbcsr_get_data, dbcsr_get_data_p, dbcsr_scalar, &
       dbcsr_scalar_fill_all, dbcsr_scalar_get_type, dbcsr_scalar_get_value, &
       dbcsr_scalar_set_type, dbcsr_scalar_zero
  USE dbcsr_dist_operations,           ONLY: create_bl_distribution,&
                                             dbcsr_create_dist_block_cyclic,&
                                             dbcsr_create_dist_r_unrot,&
                                             dbcsr_get_stored_coordinates
  USE dbcsr_error_handling,            ONLY: dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type
  USE dbcsr_index_operations,          ONLY: dbcsr_index_compact
  USE dbcsr_io,                        ONLY: dbcsr_binary_read,&
                                             dbcsr_binary_write,&
                                             dbcsr_print,&
                                             dbcsr_print_block_sum
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_methods,                   ONLY: &
       dbcsr_blk_col_offset, dbcsr_blk_column_size, dbcsr_blk_row_offset, &
       dbcsr_blk_row_size, dbcsr_distribution, dbcsr_distribution_col_dist, &
       dbcsr_distribution_has_threads, dbcsr_distribution_init, &
       dbcsr_distribution_local_cols, dbcsr_distribution_local_rows, &
       dbcsr_distribution_mp, dbcsr_distribution_ncols, &
       dbcsr_distribution_new, dbcsr_distribution_nlocal_cols, &
       dbcsr_distribution_nlocal_rows, dbcsr_distribution_nrows, &
       dbcsr_distribution_release, dbcsr_distribution_row_dist, &
       dbcsr_get_data_size, dbcsr_get_data_type, dbcsr_get_info, &
       dbcsr_get_matrix_type, dbcsr_get_num_blocks, dbcsr_get_occupation, &
       dbcsr_init, dbcsr_max_col_size, dbcsr_max_row_size, &
       dbcsr_mp_get_coordinates, dbcsr_mp_group, dbcsr_mp_hold, &
       dbcsr_mp_mynode, dbcsr_mp_mypcol, dbcsr_mp_myprow, dbcsr_mp_new, &
       dbcsr_mp_npcols, dbcsr_mp_nprows, dbcsr_mp_numnodes, dbcsr_mp_pgrid, &
       dbcsr_mp_release, dbcsr_name, dbcsr_nblkcols_local, &
       dbcsr_nblkcols_total, dbcsr_nblkrows_local, dbcsr_nblkrows_total, &
       dbcsr_nfullcols_local, dbcsr_nfullcols_total, dbcsr_nfullrows_local, &
       dbcsr_nfullrows_total, dbcsr_release, dbcsr_valid_index, &
       dbcsr_wm_use_mutable
  USE dbcsr_operations,                ONLY: &
       dbcsr_add, dbcsr_add_on_diag, dbcsr_btriu, dbcsr_clear_mempools, &
       dbcsr_copy, dbcsr_copy_into_existing, dbcsr_filter, &
       dbcsr_finalize_lib, dbcsr_frobenius_norm, dbcsr_gershgorin_norm, &
       dbcsr_get_block_diag, dbcsr_get_diag, dbcsr_hadamard_product, &
       dbcsr_init_lib, dbcsr_init_random, dbcsr_maxabs, dbcsr_multiply, &
       dbcsr_norm, dbcsr_scale, dbcsr_scale_by_vector, dbcsr_scale_mat, &
       dbcsr_set, dbcsr_set_diag, dbcsr_sum_replicated, dbcsr_trace, &
       dbcsr_triu
  USE dbcsr_operations_low,            ONLY: dbcsr_function_of_elements
  USE dbcsr_tests,                     ONLY: cp_test_multiplies
  USE dbcsr_toollib,                   ONLY: swap
  USE dbcsr_transformations,           ONLY: dbcsr_complete_redistribute,&
                                             dbcsr_desymmetrize_deep,&
                                             dbcsr_distribute,&
                                             dbcsr_make_untransposed_blocks,&
                                             dbcsr_new_transposed,&
                                             dbcsr_replicate,&
                                             dbcsr_replicate_all
  USE dbcsr_types,                     ONLY: &
       dbcsr_data_obj, dbcsr_distribution_obj, dbcsr_func_artanh, &
       dbcsr_func_dtanh, dbcsr_func_inverse, dbcsr_func_tanh, dbcsr_iterator, &
       dbcsr_mp_obj, dbcsr_no_transpose, dbcsr_norm_column, &
       dbcsr_norm_frobenius, dbcsr_norm_maxabsnorm, dbcsr_obj, &
       dbcsr_repl_full, dbcsr_scalar_type, dbcsr_type_antisymmetric, &
       dbcsr_type_complex_4, dbcsr_type_complex_8, &
       dbcsr_type_complex_default, dbcsr_type_no_symmetry, dbcsr_type_real_4, &
       dbcsr_type_real_8, dbcsr_type_real_default, dbcsr_type_symmetric
  USE dbcsr_util,                      ONLY: convert_offsets_to_sizes,&
                                             convert_sizes_to_offsets,&
                                             dbcsr_checksum,&
                                             dbcsr_verify_matrix
  USE dbcsr_work_operations,           ONLY: add_work_coordinate,&
                                             dbcsr_create,&
                                             dbcsr_finalize,&
                                             dbcsr_work_create

  PUBLIC :: array_data,&
           array_hold,&
           array_i1d_obj,&
           array_new,&
           array_release,&
           array_size,&
           array_nullify,&
           add_work_coordinate ,&
           block_add_on_diag ,&
           block_chol_inv ,&
           block_set ,&
           convert_offsets_to_sizes ,&
           convert_sizes_to_offsets ,&
           cp_test_multiplies ,&
           create_bl_distribution ,&
           cudadrv_avoid_after_busy ,&
           cudadrv_binning_binsize ,&
           cudadrv_binning_nbins ,&
           cudadrv_do_inhomogenous ,&
           cudadrv_min_flop_process ,&
           cudadrv_min_flop_sort ,&
           cudadrv_posterior_buffers ,&
           cudadrv_posterior_streams ,&
           cudadrv_priority_buffers ,&
           cudadrv_priority_streams ,&
           dbcsr_access_flush ,&
           dbcsr_add ,&
           dbcsr_add_on_diag ,&
           dbcsr_binary_read ,&
           dbcsr_binary_write ,&
           dbcsr_blk_col_offset ,&
           dbcsr_blk_column_size ,&
           dbcsr_blk_row_offset ,&
           dbcsr_blk_row_size ,&
           dbcsr_btriu ,&
           dbcsr_checksum ,&
           dbcsr_clear_mempools ,&
           dbcsr_complete_redistribute ,&
           dbcsr_copy ,&
           dbcsr_copy_into_existing ,&
           dbcsr_create ,&
           dbcsr_create_dist_block_cyclic,&
           dbcsr_create_dist_r_unrot,&
           dbcsr_data_clear_pointer ,&
           dbcsr_data_init ,&
           dbcsr_data_new ,&
           dbcsr_data_obj ,&
           dbcsr_data_release ,&
           dbcsr_desymmetrize_deep ,&
           dbcsr_distribute ,&
           dbcsr_distribution ,&
           dbcsr_distribution_col_dist ,&
           dbcsr_distribution_has_threads ,&
           dbcsr_distribution_init ,&
           dbcsr_distribution_local_cols ,&
           dbcsr_distribution_local_rows ,&
           dbcsr_distribution_mp ,&
           dbcsr_distribution_ncols ,&
           dbcsr_distribution_new ,&
           dbcsr_distribution_nlocal_cols ,&
           dbcsr_distribution_nlocal_rows ,&
           dbcsr_distribution_nrows ,&
           dbcsr_distribution_obj ,&
           dbcsr_distribution_release ,&
           dbcsr_distribution_row_dist ,&
           dbcsr_error_set ,&
           dbcsr_error_stop ,&
           dbcsr_error_type ,&
           dbcsr_filter ,&
           dbcsr_finalize ,&
           dbcsr_finalize_lib ,&
           dbcsr_frobenius_norm ,&
           dbcsr_func_artanh ,&
           dbcsr_func_dtanh ,&
           dbcsr_func_inverse ,&
           dbcsr_func_tanh ,&
           dbcsr_function_of_elements ,&
           dbcsr_gershgorin_norm ,&
           dbcsr_get_block ,&
           dbcsr_get_block_diag ,&
           dbcsr_get_block_p ,&
           dbcsr_get_conf_combtypes ,&
           dbcsr_get_conf_comm_thread_load ,&
           dbcsr_get_conf_max_ele_block ,&
           dbcsr_get_conf_mm_driver ,&
           dbcsr_get_conf_mm_stacksize ,&
           dbcsr_get_conf_mpi_mem ,&
           dbcsr_get_conf_nstacks ,&
           dbcsr_get_conf_subcomm ,&
           dbcsr_get_conf_use_comm_thread ,&
           dbcsr_get_data ,&
           dbcsr_get_data_p ,&
           dbcsr_get_data_size ,&
           dbcsr_get_data_type ,&
           dbcsr_get_diag ,&
           dbcsr_get_info ,&
           dbcsr_get_matrix_type ,&
           dbcsr_get_num_blocks ,&
           dbcsr_get_occupation ,&
           dbcsr_get_stored_coordinates ,&
           dbcsr_hadamard_product ,&
           dbcsr_index_compact ,&
           dbcsr_init ,&
           dbcsr_init_lib ,&
           dbcsr_init_random ,&
           dbcsr_iterator ,&
           dbcsr_iterator_blocks_left ,&
           dbcsr_iterator_next_block ,&
           dbcsr_iterator_start ,&
           dbcsr_iterator_stop ,&
           dbcsr_make_untransposed_blocks ,&
           dbcsr_maxabs ,&
           dbcsr_max_col_size ,&
           dbcsr_max_row_size ,&
           dbcsr_mp_get_coordinates ,&
           dbcsr_mp_group ,&
           dbcsr_mp_hold ,&
           dbcsr_mp_mynode ,&
           dbcsr_mp_mypcol ,&
           dbcsr_mp_myprow ,&
           dbcsr_mp_new ,&
           dbcsr_mp_npcols ,&
           dbcsr_mp_nprows ,&
           dbcsr_mp_numnodes ,&
           dbcsr_mp_obj ,&
           dbcsr_mp_pgrid ,&
           dbcsr_mp_release ,&
           dbcsr_multiply ,&
           dbcsr_name ,&
           dbcsr_nblkcols_local ,&
           dbcsr_nblkcols_total ,&
           dbcsr_nblkrows_local ,&
           dbcsr_nblkrows_total ,&
           dbcsr_new_transposed ,&
           dbcsr_nfullcols_local ,&
           dbcsr_nfullcols_total ,&
           dbcsr_nfullrows_local ,&
           dbcsr_nfullrows_total ,&
           dbcsr_norm ,&
           dbcsr_norm_column ,&
           dbcsr_norm_frobenius ,&
           dbcsr_norm_maxabsnorm ,&
           dbcsr_no_transpose ,&
           dbcsr_obj ,&
           dbcsr_print ,&
           dbcsr_print_block_sum ,&
           dbcsr_put_block ,&
           dbcsr_release ,&
           dbcsr_repl_full ,&
           dbcsr_replicate ,&
           dbcsr_replicate_all ,&
           dbcsr_reserve_all_blocks ,&
           dbcsr_reserve_block2d ,&
           dbcsr_reserve_blocks ,&
           dbcsr_reserve_diag_blocks ,&
           dbcsr_scalar ,&
           dbcsr_scalar_fill_all ,&
           dbcsr_scalar_get_type ,&
           dbcsr_scalar_get_value ,&
           dbcsr_scalar_set_type ,&
           dbcsr_scalar_type ,&
           dbcsr_scalar_zero ,&
           dbcsr_scale ,&
           dbcsr_scale_by_vector ,&
           dbcsr_scale_mat ,&
           dbcsr_set ,&
           dbcsr_set_conf_combtypes ,&
           dbcsr_set_conf_comm_thread_load ,&
           dbcsr_set_conf_external_timing ,&
           dbcsr_set_conf_max_ele_block ,&
           dbcsr_set_conf_mm_driver ,&
           dbcsr_set_conf_mm_stacksize ,&
           dbcsr_set_conf_mpi_mem ,&
           dbcsr_set_conf_nstacks ,&
           dbcsr_set_conf_subcomm ,&
           dbcsr_set_conf_use_comm_thread ,&
           dbcsr_set_diag ,&
           dbcsr_sum_replicated ,&
           dbcsr_trace ,&
           dbcsr_triu ,&
           dbcsr_type_antisymmetric ,&
           dbcsr_type_complex_4 ,&
           dbcsr_type_complex_8 ,&
           dbcsr_type_complex_default ,&
           dbcsr_type_no_symmetry ,&
           dbcsr_type_real_4 ,&
           dbcsr_type_real_8 ,&
           dbcsr_type_real_default ,&
           dbcsr_type_symmetric ,&
           dbcsr_valid_index ,&
           dbcsr_verify_matrix ,&
           dbcsr_wm_use_mutable ,&
           dbcsr_work_create ,&
           has_cuda ,&
           has_mpi ,&
           mm_driver_blas ,&
           mm_driver_cuda ,&
           mm_driver_matmul ,&
           mm_driver_plasma ,&
           mm_driver_smm ,&
           mm_name_blas ,&
           mm_name_cuda ,&
           mm_name_matmul ,&
           mm_name_plasma ,&
           mm_name_smm ,&
           multrec_limit ,&
           swap

END MODULE dbcsr_api

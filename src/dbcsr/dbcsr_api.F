!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   This is the start of a dbcsr_api, all publically needed functions
!>          are exported here. The others remain private to the library.
!>          Currently, this is the CP2K used set.
!>          Ultimately, a reduced subset and well defined api will remain,
!>          possibly grouped in to standard and expert api.
!>          Currently, this is work in progress.
!> \author  Joost VandeVondele
! **************************************************************************************************
MODULE dbcsr_api
   USE dbcsr_block_access,              ONLY: dbcsr_get_block_p_prv         => dbcsr_get_block_p,&
                                              dbcsr_put_block_prv           => dbcsr_put_block,&
                                              dbcsr_reserve_all_blocks_prv  => dbcsr_reserve_all_blocks,&
                                              dbcsr_reserve_block2d_prv     => dbcsr_reserve_block2d,&
                                              dbcsr_reserve_blocks_prv      => dbcsr_reserve_blocks,&
                                              dbcsr_reserve_diag_blocks_prv => dbcsr_reserve_diag_blocks
   USE dbcsr_config,                    ONLY: dbcsr_get_default_config,&
                                              dbcsr_print_config,&
                                              dbcsr_set_config
   USE dbcsr_csr_conversions,           ONLY: &
        convert_csr_to_dbcsr_prv            => convert_csr_to_dbcsr, &
        convert_dbcsr_to_csr_prv            => convert_dbcsr_to_csr, &
        csr_create_from_dbcsr_prv           => csr_create_from_dbcsr, csr_create_new, &
        csr_create_template, csr_dbcsr_blkrow_dist, csr_destroy, csr_eqrow_ceil_dist, &
        csr_eqrow_floor_dist, csr_p_type, csr_print_sparsity, csr_type, &
        dbcsr_to_csr_filter_prv             => dbcsr_to_csr_filter
   USE dbcsr_data_methods,              ONLY: dbcsr_get_data_p_prv => dbcsr_get_data_p,&
                                              dbcsr_scalar,&
                                              dbcsr_scalar_fill_all,&
                                              dbcsr_scalar_get_type,&
                                              dbcsr_scalar_get_value,&
                                              dbcsr_scalar_set_type,&
                                              dbcsr_scalar_zero
   USE dbcsr_dist_methods,              ONLY: &
        dbcsr_distribution_col_clusters, dbcsr_distribution_col_dist, &
        dbcsr_distribution_get_num_images => dbcsr_distribution_get_num_images_1d, &
        dbcsr_distribution_has_col_clusters, dbcsr_distribution_has_row_clusters, &
        dbcsr_distribution_has_threads, dbcsr_distribution_hold, dbcsr_distribution_local_cols, &
        dbcsr_distribution_local_rows, dbcsr_distribution_mp, dbcsr_distribution_ncols, &
        dbcsr_distribution_new_prv => dbcsr_distribution_new, dbcsr_distribution_nlocal_cols, &
        dbcsr_distribution_nlocal_rows, dbcsr_distribution_nrows, dbcsr_distribution_release, &
        dbcsr_distribution_row_clusters, dbcsr_distribution_row_dist
   USE dbcsr_dist_operations,           ONLY: dbcsr_get_stored_coordinates_prv => dbcsr_get_stored_coordinates
   USE dbcsr_io,                        ONLY: dbcsr_binary_read_prv            => dbcsr_binary_read,&
                                              dbcsr_binary_write_prv           => dbcsr_binary_write,&
                                              dbcsr_print_block_sum_prv        => dbcsr_print_block_sum,&
                                              dbcsr_print_prv                  => dbcsr_print
   USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left_prv   => dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block_prv    => dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start_prv         => dbcsr_iterator_start,&
                                              dbcsr_iterator_stop_prv          => dbcsr_iterator_stop
   USE dbcsr_lib,                       ONLY: dbcsr_clear_mempools,&
                                              dbcsr_finalize_lib,&
                                              dbcsr_init_lib
   USE dbcsr_methods,                   ONLY: &
        dbcsr_distribution_prv              => dbcsr_distribution, &
        dbcsr_get_data_size_prv             => dbcsr_get_data_size, &
        dbcsr_get_data_type_prv             => dbcsr_get_data_type, &
        dbcsr_get_matrix_type_prv           => dbcsr_get_matrix_type, &
        dbcsr_get_num_blocks_prv            => dbcsr_get_num_blocks, &
        dbcsr_has_symmetry_prv              => dbcsr_has_symmetry, &
        dbcsr_max_col_size_prv              => dbcsr_max_col_size, &
        dbcsr_max_row_size_prv              => dbcsr_max_row_size, &
        dbcsr_name_prv                      => dbcsr_name, &
        dbcsr_nblkcols_total_prv            => dbcsr_nblkcols_total, &
        dbcsr_nblkrows_total_prv            => dbcsr_nblkrows_total, &
        dbcsr_nfullcols_total_prv           => dbcsr_nfullcols_total, &
        dbcsr_nfullrows_total_prv           => dbcsr_nfullrows_total, &
        dbcsr_release_prv                   => dbcsr_release, &
        dbcsr_setname_prv                   => dbcsr_setname, &
        dbcsr_valid_index_prv               => dbcsr_valid_index, dbcsr_wm_use_mutable
   USE dbcsr_mp_methods,                ONLY: &
        dbcsr_mp_get_process, dbcsr_mp_grid_setup, dbcsr_mp_group, dbcsr_mp_hold, dbcsr_mp_mynode, &
        dbcsr_mp_new, dbcsr_mp_npcols, dbcsr_mp_nprows, dbcsr_mp_numnodes, dbcsr_mp_release
   USE dbcsr_multiply_api,              ONLY: dbcsr_multiply_prv => dbcsr_multiply
   USE dbcsr_operations,                ONLY: &
        dbcsr_add_on_diag_prv               => dbcsr_add_on_diag, &
        dbcsr_add_prv                       => dbcsr_add, &
        dbcsr_copy_into_existing_prv        => dbcsr_copy_into_existing, &
        dbcsr_copy_prv                      => dbcsr_copy, &
        dbcsr_filter_prv                    => dbcsr_filter, &
        dbcsr_frobenius_norm_prv            => dbcsr_frobenius_norm, &
        dbcsr_function_of_elements_prv      => dbcsr_function_of_elements, &
        dbcsr_gershgorin_norm_prv           => dbcsr_gershgorin_norm, &
        dbcsr_get_block_diag_prv            => dbcsr_get_block_diag, &
        dbcsr_get_diag_prv                  => dbcsr_get_diag, &
        dbcsr_get_info_prv                  => dbcsr_get_info, &
        dbcsr_get_occupation_prv            => dbcsr_get_occupation, &
        dbcsr_hadamard_product_prv          => dbcsr_hadamard_product, &
        dbcsr_init_random_prv               => dbcsr_init_random, &
        dbcsr_maxabs_prv                    => dbcsr_maxabs, &
        dbcsr_norm_prv                      => dbcsr_norm, &
        dbcsr_scale_by_vector_prv           => dbcsr_scale_by_vector, &
        dbcsr_scale_prv                     => dbcsr_scale, &
        dbcsr_set_diag_prv                  => dbcsr_set_diag, &
        dbcsr_set_prv                       => dbcsr_set, &
        dbcsr_sum_replicated_prv            => dbcsr_sum_replicated, &
        dbcsr_trace_prv                     => dbcsr_trace, &
        dbcsr_triu_prv                      => dbcsr_triu
   USE dbcsr_test_methods,              ONLY: dbcsr_reset_randmat_seed
   USE dbcsr_tests,                     ONLY: dbcsr_run_tests,&
                                              dbcsr_test_binary_io,&
                                              dbcsr_test_mm
   USE dbcsr_toollib,                   ONLY: uppercase
   USE dbcsr_transformations,           ONLY: dbcsr_complete_redistribute_prv => dbcsr_complete_redistribute,&
                                              dbcsr_desymmetrize_deep_prv     => dbcsr_desymmetrize_deep,&
                                              dbcsr_distribute_prv            => dbcsr_distribute,&
                                              dbcsr_new_transposed_prv        => dbcsr_new_transposed,&
                                              dbcsr_replicate_all_prv         => dbcsr_replicate_all
   USE dbcsr_types,                     ONLY: &
        dbcsr_distribution_obj, dbcsr_func_artanh, dbcsr_func_dtanh, dbcsr_func_inverse, &
        dbcsr_func_tanh, dbcsr_iterator_prv => dbcsr_iterator, dbcsr_mp_obj, dbcsr_no_transpose, &
        dbcsr_norm_column, dbcsr_norm_frobenius, dbcsr_norm_maxabsnorm, &
        dbcsr_prv_type => dbcsr_type, dbcsr_scalar_type, dbcsr_type_antisymmetric, &
        dbcsr_type_complex_4, dbcsr_type_complex_8, dbcsr_type_complex_default, &
        dbcsr_type_no_symmetry, dbcsr_type_real_4, dbcsr_type_real_8, dbcsr_type_real_default, &
        dbcsr_type_symmetric
   USE dbcsr_util,                      ONLY: convert_offsets_to_sizes,&
                                              convert_sizes_to_offsets,&
                                              dbcsr_checksum_prv      => dbcsr_checksum,&
                                              dbcsr_verify_matrix_prv => dbcsr_verify_matrix
   USE dbcsr_vector_operations,         ONLY: create_col_vec_from_matrix_prv            => create_col_vec_from_matrix,&
                                              create_replicated_col_vec_from_matrix_prv => create_replicated_col_vec_from_matrix,&
                                              create_replicated_row_vec_from_matrix_prv => create_replicated_row_vec_from_matrix,&
                                              create_row_vec_from_matrix_prv            => create_row_vec_from_matrix,&
                                              dbcsr_matrix_colvec_multiply_prv          => dbcsr_matrix_colvec_multiply
   USE dbcsr_work_operations,           ONLY: add_work_coordinate_prv => add_work_coordinate,&
                                              dbcsr_create_prv        => dbcsr_create,&
                                              dbcsr_finalize_prv      => dbcsr_finalize,&
                                              dbcsr_work_create_prv   => dbcsr_work_create
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8,&
                                              real_4,&
                                              real_8

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
#include "../base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_api'

   ! constants
   PUBLIC :: dbcsr_type_no_symmetry
   PUBLIC :: dbcsr_type_symmetric
   PUBLIC :: dbcsr_type_antisymmetric
   PUBLIC :: dbcsr_type_complex_8
   PUBLIC :: dbcsr_type_real_4
   PUBLIC :: dbcsr_type_real_8
   PUBLIC :: dbcsr_type_complex_4
   PUBLIC :: dbcsr_type_complex_default
   PUBLIC :: dbcsr_type_real_default

   ! types
   PUBLIC :: dbcsr_type
   PUBLIC :: dbcsr_p_type
   PUBLIC :: dbcsr_mp_obj
   PUBLIC :: dbcsr_distribution_obj
   PUBLIC :: dbcsr_iterator_type

   ! lib init/finalize
   PUBLIC :: dbcsr_clear_mempools
   PUBLIC :: dbcsr_init_lib
   PUBLIC :: dbcsr_finalize_lib
   PUBLIC :: dbcsr_set_config
   PUBLIC :: dbcsr_get_default_config
   PUBLIC :: dbcsr_print_config
   PUBLIC :: dbcsr_reset_randmat_seed

   ! create / release
   PUBLIC :: dbcsr_mp_new
   PUBLIC :: dbcsr_mp_release
   PUBLIC :: dbcsr_distribution_release
   PUBLIC :: dbcsr_distribution_new
   PUBLIC :: dbcsr_create
   PUBLIC :: dbcsr_init
   PUBLIC :: dbcsr_init_p
   PUBLIC :: dbcsr_release
   PUBLIC :: dbcsr_release_p
   PUBLIC :: dbcsr_deallocate_matrix
   PUBLIC :: dbcsr_mp_hold

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
   PUBLIC :: dbcsr_complete_redistribute
   PUBLIC :: dbcsr_get_block_p

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
   PUBLIC :: dbcsr_distribution
   PUBLIC :: dbcsr_distribution_local_cols
   PUBLIC :: dbcsr_distribution_hold
   PUBLIC :: dbcsr_distribution_nlocal_rows
   PUBLIC :: dbcsr_distribution_nrows
   PUBLIC :: dbcsr_distribution_local_rows
   PUBLIC :: dbcsr_distribution_ncols
   PUBLIC :: dbcsr_distribution_nlocal_cols
   PUBLIC :: dbcsr_distribution_get_num_images
   PUBLIC :: dbcsr_distribution_row_dist
   PUBLIC :: dbcsr_distribution_col_dist
   PUBLIC :: dbcsr_distribution_row_clusters
   PUBLIC :: dbcsr_distribution_col_clusters
   PUBLIC :: dbcsr_distribution_has_row_clusters
   PUBLIC :: dbcsr_distribution_has_col_clusters
   PUBLIC :: dbcsr_distribution_has_threads
   PUBLIC :: dbcsr_distribution_mp
   PUBLIC :: dbcsr_mp_group
   PUBLIC :: dbcsr_mp_mynode
   PUBLIC :: dbcsr_mp_numnodes
   PUBLIC :: dbcsr_mp_nprows
   PUBLIC :: dbcsr_mp_npcols
   PUBLIC :: dbcsr_mp_get_process
   PUBLIC :: dbcsr_mp_grid_setup
   PUBLIC :: dbcsr_setname
   PUBLIC :: dbcsr_name
   PUBLIC :: dbcsr_max_row_size
   PUBLIC :: dbcsr_max_col_size
   PUBLIC :: dbcsr_col_block_sizes
   PUBLIC :: dbcsr_row_block_sizes
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
   PUBLIC :: dbcsr_add_work_coordinate
   PUBLIC :: dbcsr_get_wms_data_p
   PUBLIC :: dbcsr_get_data_p
   PUBLIC :: dbcsr_set_work_size
   PUBLIC :: dbcsr_finalize

   ! replication
   PUBLIC :: dbcsr_replicate_all
   PUBLIC :: dbcsr_sum_replicated
   PUBLIC :: dbcsr_distribute

   ! matrix set
   PUBLIC :: dbcsr_allocate_matrix_set
   PUBLIC :: dbcsr_deallocate_matrix_set

   ! misc
   PUBLIC :: convert_offsets_to_sizes
   PUBLIC :: convert_sizes_to_offsets
   PUBLIC :: dbcsr_run_tests
   PUBLIC :: dbcsr_test_mm

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
   PUBLIC :: csr_type
   PUBLIC :: csr_p_type
   PUBLIC :: convert_csr_to_dbcsr
   PUBLIC :: convert_dbcsr_to_csr
   PUBLIC :: csr_create_from_dbcsr
   PUBLIC :: csr_destroy
   PUBLIC :: csr_create
   PUBLIC :: csr_eqrow_floor_dist
   PUBLIC :: csr_eqrow_ceil_dist
   PUBLIC :: csr_dbcsr_blkrow_dist
   PUBLIC :: csr_print_sparsity
   PUBLIC :: dbcsr_to_csr_filter

   ! binary io
   PUBLIC :: dbcsr_binary_write
   PUBLIC :: dbcsr_binary_read
   PUBLIC :: dbcsr_test_binary_io

   ! matrix vector multiplication
   PUBLIC :: dbcsr_matrix_colvec_multiply
   PUBLIC :: create_col_vec_from_matrix
   PUBLIC :: create_row_vec_from_matrix
   PUBLIC :: create_replicated_col_vec_from_matrix
   PUBLIC :: create_replicated_row_vec_from_matrix

   ! -----------------------------------------------------------------------------------------------
   TYPE dbcsr_p_type
      TYPE(dbcsr_type), POINTER :: matrix => Null()
   END TYPE dbcsr_p_type

   ! the components of this type must remain private to encapsulate better the internals
   ! of the dbcsr library.
   TYPE dbcsr_type
      INTEGER, PRIVATE                     :: ref_count = 0
      TYPE(dbcsr_prv_type), PRIVATE        :: prv
   END TYPE dbcsr_type

   TYPE dbcsr_iterator_type
      TYPE(dbcsr_iterator_prv), PRIVATE        :: prv
   END TYPE dbcsr_iterator_type

   INTERFACE dbcsr_distribution_new
      MODULE PROCEDURE dbcsr_distribution_new_noclusters
      MODULE PROCEDURE dbcsr_distribution_new_clusters1
      MODULE PROCEDURE dbcsr_distribution_new_clusters2
   END INTERFACE

   INTERFACE dbcsr_create
      MODULE PROCEDURE dbcsr_create_new
      MODULE PROCEDURE dbcsr_create_template
   END INTERFACE

   INTERFACE dbcsr_conform_scalar
      MODULE PROCEDURE make_conformant_scalar_d, &
         make_conformant_scalar_s, &
         make_conformant_scalar_c, &
         make_conformant_scalar_z
   END INTERFACE

   INTERFACE dbcsr_trace
      MODULE PROCEDURE dbcsr_trace_ab_d, &
         dbcsr_trace_ab_s
      MODULE PROCEDURE dbcsr_trace_a_d, &
         dbcsr_trace_a_s
      MODULE PROCEDURE dbcsr_trace_ab_z, &
         dbcsr_trace_ab_c
      MODULE PROCEDURE dbcsr_trace_a_z, &
         dbcsr_trace_a_c
   END INTERFACE

   INTERFACE dbcsr_set
      MODULE PROCEDURE dbcsr_set_d, &
         dbcsr_set_s, &
         dbcsr_set_c, &
         dbcsr_set_z
   END INTERFACE

   INTERFACE dbcsr_add
      MODULE PROCEDURE dbcsr_add_d, &
         dbcsr_add_s, &
         dbcsr_add_c, &
         dbcsr_add_z
   END INTERFACE

   INTERFACE dbcsr_scale
      MODULE PROCEDURE dbcsr_scale_d, &
         dbcsr_scale_s, &
         dbcsr_scale_c, &
         dbcsr_scale_z
   END INTERFACE

   INTERFACE dbcsr_scale_by_vector
      MODULE PROCEDURE dbcsr_scale_by_vector_d, &
         dbcsr_scale_by_vector_s, &
         dbcsr_scale_by_vector_c, &
         dbcsr_scale_by_vector_z
   END INTERFACE

   INTERFACE dbcsr_multiply
      MODULE PROCEDURE dbcsr_multiply_d, &
         dbcsr_multiply_s, &
         dbcsr_multiply_c, &
         dbcsr_multiply_z
   END INTERFACE

   INTERFACE dbcsr_get_block_p
      MODULE PROCEDURE dbcsr_get_block_p_d, &
         dbcsr_get_block_p_s, &
         dbcsr_get_block_p_z, &
         dbcsr_get_block_p_c
      MODULE PROCEDURE dbcsr_get_2d_block_p_d, &
         dbcsr_get_2d_block_p_s, &
         dbcsr_get_2d_block_p_z, &
         dbcsr_get_2d_block_p_c
   END INTERFACE

   INTERFACE dbcsr_put_block
      MODULE PROCEDURE dbcsr_put_block_d, &
         dbcsr_put_block_s, &
         dbcsr_put_block_z, &
         dbcsr_put_block_c
      MODULE PROCEDURE dbcsr_put_block2d_d, &
         dbcsr_put_block2d_s, &
         dbcsr_put_block2d_z, &
         dbcsr_put_block2d_c
   END INTERFACE

   INTERFACE dbcsr_iterator_next_block
      MODULE PROCEDURE dbcsr_iterator_next_block_index
      MODULE PROCEDURE dbcsr_iterator_next_2d_block_d, &
         dbcsr_iterator_next_2d_block_s, &
         dbcsr_iterator_next_2d_block_c, &
         dbcsr_iterator_next_2d_block_z
      MODULE PROCEDURE dbcsr_iterator_next_1d_block_d, &
         dbcsr_iterator_next_1d_block_s, &
         dbcsr_iterator_next_1d_block_c, &
         dbcsr_iterator_next_1d_block_z
   END INTERFACE

   INTERFACE dbcsr_reserve_block2d
      MODULE PROCEDURE dbcsr_reserve_block2d_d, &
         dbcsr_reserve_block2d_s, &
         dbcsr_reserve_block2d_c, &
         dbcsr_reserve_block2d_z
   END INTERFACE

   INTERFACE dbcsr_allocate_matrix_set
      MODULE PROCEDURE allocate_dbcsr_matrix_set, allocate_dbcsr_matrix_set_2d, allocate_dbcsr_matrix_set_3d
   END INTERFACE

   INTERFACE dbcsr_deallocate_matrix_set
      MODULE PROCEDURE deallocate_dbcsr_matrix_set, &
         deallocate_dbcsr_matrix_set_2d, deallocate_dbcsr_matrix_set_3d
   END INTERFACE

   INTERFACE csr_create
      MODULE PROCEDURE csr_create_new, &
         csr_create_template
   END INTERFACE

   INTERFACE dbcsr_get_wms_data_p
      MODULE PROCEDURE dbcsr_get_wms_data_s, dbcsr_get_wms_data_c, &
         dbcsr_get_wms_data_d, dbcsr_get_wms_data_z
   END INTERFACE

   INTERFACE dbcsr_get_data_p
      MODULE PROCEDURE dbcsr_get_data_s, dbcsr_get_data_c, &
         dbcsr_get_data_d, dbcsr_get_data_z
   END INTERFACE

   INTERFACE dbcsr_matrix_colvec_multiply
      MODULE PROCEDURE dbcsr_matrix_colvec_multiply_s, dbcsr_matrix_colvec_multiply_c, &
         dbcsr_matrix_colvec_multiply_d, dbcsr_matrix_colvec_multiply_z
   END INTERFACE

   PRIVATE

CONTAINS

! **************************************************************************************************
!> \brief cp2k interface to libdbcsr.
!>     The cp2k error is passed as an argument.
!>     Those procedure can de timed with the cp2k timer.
!> \param matrix ...
!> \retval name ...
!> \author vw
! **************************************************************************************************
   FUNCTION dbcsr_name(matrix) RESULT(name)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      CHARACTER(len=default_string_length)               :: name

      name = dbcsr_name_prv(matrix%prv)
   END FUNCTION dbcsr_name

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param newname ...
! **************************************************************************************************
   SUBROUTINE dbcsr_setname(matrix, newname)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      CHARACTER(len=*), INTENT(IN)                       :: newname

      CALL dbcsr_setname_prv(matrix%prv, newname)
   END SUBROUTINE dbcsr_setname

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval norm ...
! **************************************************************************************************
   FUNCTION dbcsr_gershgorin_norm(matrix) RESULT(norm)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=real_8)                                  :: norm

      norm = dbcsr_gershgorin_norm_prv(matrix%prv)
   END FUNCTION dbcsr_gershgorin_norm

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param local ...
!> \retval norm ...
! **************************************************************************************************
   FUNCTION dbcsr_frobenius_norm(matrix, local) RESULT(norm)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      LOGICAL, INTENT(in), OPTIONAL                      :: local
      REAL(KIND=real_8)                                  :: norm

      norm = dbcsr_frobenius_norm_prv(matrix%prv, local)
   END FUNCTION dbcsr_frobenius_norm

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval norm ...
! **************************************************************************************************
   FUNCTION dbcsr_maxabs(matrix) RESULT(norm)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=real_8)                                  :: norm

      norm = dbcsr_maxabs_prv(matrix%prv)
   END FUNCTION dbcsr_maxabs

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param redist ...
!> \param keep_sparsity ...
! **************************************************************************************************
   SUBROUTINE dbcsr_complete_redistribute(matrix, redist, keep_sparsity)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(dbcsr_type), INTENT(INOUT)                    :: redist
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_sparsity

      CALL dbcsr_complete_redistribute_prv(matrix%prv, redist%prv, keep_sparsity)
   END SUBROUTINE dbcsr_complete_redistribute

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param rows ...
!> \param cols ...
!> \param blk_pointers ...
! **************************************************************************************************
   SUBROUTINE dbcsr_reserve_blocks(matrix, rows, cols, blk_pointers)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, DIMENSION(:), INTENT(IN)                  :: rows, cols
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: blk_pointers

      CALL dbcsr_reserve_blocks_prv(matrix%prv, rows, cols, blk_pointers)
   END SUBROUTINE dbcsr_reserve_blocks

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_reserve_all_blocks(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CALL dbcsr_reserve_all_blocks_prv(matrix%prv)
   END SUBROUTINE dbcsr_reserve_all_blocks

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_reserve_diag_blocks(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CALL dbcsr_reserve_diag_blocks_prv(matrix%prv)
   END SUBROUTINE dbcsr_reserve_diag_blocks

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param index_matrix ...
!> \param row ...
!> \param col ...
!> \param blk ...
!> \param index ...
! **************************************************************************************************
   SUBROUTINE dbcsr_add_work_coordinate(matrix, index_matrix, row, col, blk, index)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN)                                :: index_matrix, row, col
      INTEGER, INTENT(IN), OPTIONAL                      :: blk
      INTEGER, INTENT(OUT), OPTIONAL                     :: index

      CALL add_work_coordinate_prv(matrix%prv%wms(index_matrix), row, col, blk, index)
   END SUBROUTINE dbcsr_add_work_coordinate

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param index_matrix ...
!> \param newvalue ...
! **************************************************************************************************
   SUBROUTINE dbcsr_set_work_size(matrix, index_matrix, newvalue)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN)                                :: index_matrix, newvalue

      matrix%prv%wms(index_matrix)%datasize = newvalue
   END SUBROUTINE dbcsr_set_work_size

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param keep_sparsity ...
! **************************************************************************************************
   SUBROUTINE dbcsr_init_random(matrix, keep_sparsity)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      LOGICAL, OPTIONAL                                  :: keep_sparsity

      CALL dbcsr_init_random_prv(matrix%prv, keep_sparsity=keep_sparsity)
   END SUBROUTINE dbcsr_init_random

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval data_type ...
! **************************************************************************************************
   PURE FUNCTION dbcsr_get_data_type(matrix) RESULT(data_type)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: data_type

      data_type = dbcsr_get_data_type_prv(matrix%prv)
   END FUNCTION dbcsr_get_data_type

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval valid_index ...
! **************************************************************************************************
   PURE FUNCTION dbcsr_valid_index(matrix) RESULT(valid_index)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      LOGICAL                                            :: valid_index

      valid_index = dbcsr_valid_index_prv(matrix%prv)
   END FUNCTION dbcsr_valid_index

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
      INTEGER, INTENT(OUT), OPTIONAL                     :: processor

      CALL dbcsr_get_stored_coordinates_prv(matrix%prv, row, column, processor)
   END SUBROUTINE dbcsr_get_stored_coordinates

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval num_blocks ...
! **************************************************************************************************
   PURE FUNCTION dbcsr_get_num_blocks(matrix) RESULT(num_blocks)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: num_blocks

      num_blocks = dbcsr_get_num_blocks_prv(matrix%prv)
   END FUNCTION dbcsr_get_num_blocks

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval data_size ...
! **************************************************************************************************
   FUNCTION dbcsr_get_data_size(matrix) RESULT(data_size)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: data_size

      data_size = dbcsr_get_data_size_prv(matrix%prv)
   END FUNCTION dbcsr_get_data_size

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval matrix_type ...
! **************************************************************************************************
   PURE FUNCTION dbcsr_get_matrix_type(matrix) RESULT(matrix_type)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      CHARACTER                                          :: matrix_type

      matrix_type = dbcsr_get_matrix_type_prv(matrix%prv)
   END FUNCTION dbcsr_get_matrix_type

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval occupation ...
! **************************************************************************************************
   FUNCTION dbcsr_get_occupation(matrix) RESULT(occupation)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      REAL(KIND=real_8)                                  :: occupation

      occupation = dbcsr_get_occupation_prv(matrix%prv)
   END FUNCTION dbcsr_get_occupation

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nblkrows_total ...
! **************************************************************************************************
   FUNCTION dbcsr_nblkrows_total(matrix) RESULT(nblkrows_total)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: nblkrows_total

      nblkrows_total = dbcsr_nblkrows_total_prv(matrix%prv)
   END FUNCTION dbcsr_nblkrows_total

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nblkcols_total ...
! **************************************************************************************************
   FUNCTION dbcsr_nblkcols_total(matrix) RESULT(nblkcols_total)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: nblkcols_total

      nblkcols_total = dbcsr_nblkcols_total_prv(matrix%prv)
   END FUNCTION dbcsr_nblkcols_total

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nfullrows_total ...
! **************************************************************************************************
   FUNCTION dbcsr_nfullrows_total(matrix) RESULT(nfullrows_total)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: nfullrows_total

      nfullrows_total = dbcsr_nfullrows_total_prv(matrix%prv)
   END FUNCTION dbcsr_nfullrows_total

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nfullcols_total ...
! **************************************************************************************************
   FUNCTION dbcsr_nfullcols_total(matrix) RESULT(nfullcols_total)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: nfullcols_total

      nfullcols_total = dbcsr_nfullcols_total_prv(matrix%prv)
   END FUNCTION dbcsr_nfullcols_total

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval max_row_size ...
! **************************************************************************************************
   FUNCTION dbcsr_max_row_size(matrix) RESULT(max_row_size)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: max_row_size

      max_row_size = dbcsr_max_row_size_prv(matrix%prv)
   END FUNCTION dbcsr_max_row_size

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval max_col_size ...
! **************************************************************************************************
   FUNCTION dbcsr_max_col_size(matrix) RESULT(max_col_size)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER                                            :: max_col_size

      max_col_size = dbcsr_max_col_size_prv(matrix%prv)
   END FUNCTION dbcsr_max_col_size

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval distribution ...
! **************************************************************************************************
   FUNCTION dbcsr_distribution(matrix) RESULT(distribution)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(dbcsr_distribution_obj)                       :: distribution

      distribution = dbcsr_distribution_prv(matrix%prv)
   END FUNCTION dbcsr_distribution

! **************************************************************************************************
!> \brief ...
!> \param iterator ...
!> \retval blocks_left ...
! **************************************************************************************************
   PURE FUNCTION dbcsr_iterator_blocks_left(iterator) RESULT(blocks_left)
      TYPE(dbcsr_iterator_type), INTENT(IN)              :: iterator
      LOGICAL                                            :: blocks_left

      blocks_left = dbcsr_iterator_blocks_left_prv(iterator%prv)
   END FUNCTION dbcsr_iterator_blocks_left

! **************************************************************************************************
!> \brief ...
!> \param iterator ...
! **************************************************************************************************
   SUBROUTINE dbcsr_iterator_stop(iterator)
      TYPE(dbcsr_iterator_type), INTENT(INOUT)           :: iterator

      CALL dbcsr_iterator_stop_prv(iterator%prv)
   END SUBROUTINE dbcsr_iterator_stop

! **************************************************************************************************
!> \brief ...
!> \param iterator ...
!> \param matrix ...
!> \param shared ...
!> \param dynamic ...
!> \param dynamic_byrows ...
!> \param contiguous_pointers ...
!> \param read_only ...
! **************************************************************************************************
   SUBROUTINE dbcsr_iterator_start(iterator, matrix, shared, dynamic, &
                                   dynamic_byrows, contiguous_pointers, read_only)
      TYPE(dbcsr_iterator_type), INTENT(OUT)             :: iterator
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: shared, dynamic, dynamic_byrows, &
                                                            contiguous_pointers, read_only

      CALL dbcsr_iterator_start_prv(iterator%prv, matrix%prv, shared, dynamic, &
                                    dynamic_byrows, contiguous_pointers, read_only)
   END SUBROUTINE dbcsr_iterator_start

! **************************************************************************************************
!> \brief Gets the index information of the next block, no data.
!> \param[in,out] iterator   the iterator
!> \param[out] row           row of the data block
!> \param[out] column        column of the data block
!> \param[out] blk           block number
!> \param[out] blk_p         (optional) index into block data array
! **************************************************************************************************
   SUBROUTINE dbcsr_iterator_next_block_index(iterator, row, column, blk, blk_p)
      TYPE(dbcsr_iterator_type), INTENT(INOUT)           :: iterator
      INTEGER, INTENT(OUT)                               :: row, column, blk
      INTEGER, INTENT(OUT), OPTIONAL                     :: blk_p

      CALL dbcsr_iterator_next_block_prv(iterator%prv, row=row, column=column, blk=blk, blk_p=blk_p)
   END SUBROUTINE dbcsr_iterator_next_block_index

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
!> \param data_type ...
! **************************************************************************************************
   SUBROUTINE dbcsr_get_info(matrix, nblkrows_total, nblkcols_total, &
                             nfullrows_total, nfullcols_total, &
                             nblkrows_local, nblkcols_local, &
                             nfullrows_local, nfullcols_local, &
                             my_prow, my_pcol, &
                             local_rows, local_cols, proc_row_dist, proc_col_dist, &
                             row_blk_size, col_blk_size, row_blk_offset, col_blk_offset, &
                             distribution, name, matrix_type, data_type)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER, INTENT(OUT), OPTIONAL :: nblkrows_total, nblkcols_total, nfullrows_total, &
         nfullcols_total, nblkrows_local, nblkcols_local, nfullrows_local, nfullcols_local, &
         my_prow, my_pcol
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: local_rows, local_cols, proc_row_dist, &
                                                            proc_col_dist
      INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL, &
         POINTER                                         :: row_blk_size, col_blk_size, &
                                                            row_blk_offset, col_blk_offset
      TYPE(dbcsr_distribution_obj), INTENT(OUT), &
         OPTIONAL                                        :: distribution
      CHARACTER(len=*), INTENT(OUT), OPTIONAL            :: name
      CHARACTER, OPTIONAL                                :: matrix_type
      INTEGER, OPTIONAL                                  :: data_type

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
                              distribution=distribution, &
                              name=name, &
                              matrix_type=matrix_type, &
                              data_type=data_type)
   END SUBROUTINE dbcsr_get_info

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param which_norm ...
!> \param norm_scalar ...
!> \param norm_vector ...
! **************************************************************************************************
   SUBROUTINE dbcsr_norm(matrix, which_norm, norm_scalar, norm_vector)

      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: matrix
      INTEGER, INTENT(IN)                                :: which_norm
      REAL(dp), INTENT(OUT), OPTIONAL                    :: norm_scalar
      REAL(dp), DIMENSION(:), INTENT(OUT), OPTIONAL      :: norm_vector

      IF (PRESENT(norm_scalar)) THEN
         CALL dbcsr_norm_prv(matrix%prv, which_norm, norm_scalar=norm_scalar)
      ELSEIF (PRESENT(norm_vector)) THEN
         CALL dbcsr_norm_prv(matrix%prv, which_norm, norm_vector=norm_vector)
      ELSE
         CPABORT("Must pass either scalar or vector norm.")
      ENDIF
   END SUBROUTINE dbcsr_norm

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_replicate_all(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CALL dbcsr_replicate_all_prv(matrix%prv)
   END SUBROUTINE dbcsr_replicate_all

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param fast ...
! **************************************************************************************************
   SUBROUTINE dbcsr_distribute(matrix, fast)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      LOGICAL, INTENT(in), OPTIONAL                      :: fast

      CALL dbcsr_distribute_prv(matrix%prv, fast)
   END SUBROUTINE dbcsr_distribute

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_release_p(matrix)
      TYPE(dbcsr_type), POINTER                          :: matrix

      IF (ASSOCIATED(matrix)) THEN
         CALL dbcsr_release(matrix)
         DEALLOCATE (matrix)
      ENDIF

   END SUBROUTINE dbcsr_release_p

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_release(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CALL dbcsr_release_prv(matrix%prv)

   END SUBROUTINE dbcsr_release

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_init(matrix)
      TYPE(dbcsr_type), INTENT(OUT)                      :: matrix

      matrix%ref_count = 0
   END SUBROUTINE dbcsr_init

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_init_p(matrix)
      TYPE(dbcsr_type), POINTER                          :: matrix

      IF (ASSOCIATED(matrix)) THEN
         CALL dbcsr_release(matrix)
         DEALLOCATE (matrix)
      ENDIF

      ALLOCATE (matrix)
      CALL dbcsr_init(matrix)
   END SUBROUTINE dbcsr_init_p

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param nodata ...
!> \param matlab_format ...
!> \param variable_name ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE dbcsr_print(matrix, nodata, matlab_format, variable_name, unit_nr)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: nodata, matlab_format
      CHARACTER(*), INTENT(in), OPTIONAL                 :: variable_name
      INTEGER, OPTIONAL                                  :: unit_nr

      CALL dbcsr_print_prv(matrix%prv, nodata, matlab_format, variable_name, unit_nr)
   END SUBROUTINE dbcsr_print

! **************************************************************************************************
!> \brief Prints the sum of the elements in each block
!> \param matrix ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE dbcsr_print_block_sum(matrix, unit_nr)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER, OPTIONAL                                  :: unit_nr

      CALL dbcsr_print_block_sum_prv(matrix%prv, unit_nr)
   END SUBROUTINE dbcsr_print_block_sum

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param local ...
!> \param pos ...
!> \retval checksum ...
! **************************************************************************************************
   FUNCTION dbcsr_checksum(matrix, local, pos) RESULT(checksum)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: local, pos
      REAL(KIND=dp)                                      :: checksum

      checksum = dbcsr_checksum_prv(matrix%prv, local=local, pos=pos)
   END FUNCTION dbcsr_checksum

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_sum_replicated(matrix)
      TYPE(dbcsr_type), INTENT(inout)                    :: matrix

      CALL dbcsr_sum_replicated_prv(matrix%prv)
   END SUBROUTINE dbcsr_sum_replicated

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE dbcsr_triu(matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      CALL dbcsr_triu_prv(matrix%prv)
   END SUBROUTINE dbcsr_triu

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

      CALL dbcsr_verify_matrix_prv(matrix%prv, verbosity, local)
   END SUBROUTINE dbcsr_verify_matrix

! **************************************************************************************************
!> \brief Creates new distribution from block and cluster distributions
!> \param[out] dist           distribution
!> \param[in] mp_env          multiprocessing environment
!> \param row_dist_block ...
!> \param col_dist_block ...
!> \param row_dist_cluster ...
!> \param col_dist_cluster ...
!> \param reuse_arrays ...
! **************************************************************************************************
   SUBROUTINE dbcsr_distribution_new_clusters1(dist, mp_env, &
                                               row_dist_block, col_dist_block, &
                                               row_dist_cluster, col_dist_cluster, &
                                               reuse_arrays)
      TYPE(dbcsr_distribution_obj), INTENT(OUT)          :: dist
      TYPE(dbcsr_mp_obj), INTENT(IN)                     :: mp_env
      INTEGER, DIMENSION(:), INTENT(IN), POINTER         :: row_dist_block, col_dist_block, &
                                                            row_dist_cluster, col_dist_cluster
      LOGICAL, INTENT(IN), OPTIONAL                      :: reuse_arrays

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_distribution_new_clusters1', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: timing_handle
      INTEGER, DIMENSION(:), POINTER                     :: col_dist_cluster2, row_dist_cluster2

      CALL timeset(routineN, timing_handle)

      row_dist_cluster2 => row_dist_cluster
      col_dist_cluster2 => col_dist_cluster
      IF (ASSOCIATED(row_dist_cluster2)) THEN
         IF (ANY(row_dist_cluster2 .EQ. -1)) NULLIFY (row_dist_cluster2)
      END IF
      IF (ASSOCIATED(col_dist_cluster2)) THEN
         IF (ANY(col_dist_cluster2 .EQ. -1)) NULLIFY (col_dist_cluster2)
      END IF

      CALL dbcsr_distribution_new_prv(dist, mp_env, row_dist_block, col_dist_block, &
                                      row_dist_cluster=row_dist_cluster2, &
                                      col_dist_cluster=col_dist_cluster2, &
                                      reuse_arrays=reuse_arrays)

      CALL timestop(timing_handle)

   END SUBROUTINE dbcsr_distribution_new_clusters1

! **************************************************************************************************
!> \brief Creates new distribution from block and cluster distributions
!> \param[out] dist           distribution
!> \param[in] mp_env          multiprocessing environment
!> \param row_dist row distribution
!> \param col_dist column distribution
!> \param reuse_arrays ...
! **************************************************************************************************
   SUBROUTINE dbcsr_distribution_new_clusters2(dist, mp_env, row_dist, col_dist, reuse_arrays)
      TYPE(dbcsr_distribution_obj), INTENT(OUT)          :: dist
      TYPE(dbcsr_mp_obj), INTENT(IN)                     :: mp_env
      INTEGER, DIMENSION(:, :), INTENT(IN), POINTER      :: row_dist, col_dist
      LOGICAL, INTENT(IN), OPTIONAL                      :: reuse_arrays

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_distribution_new_clusters2', &
         routineP = moduleN//':'//routineN

      INTEGER, DIMENSION(:), POINTER                     :: col_dist_block, col_dist_cluster, &
                                                            row_dist_block, row_dist_cluster

      row_dist_block => row_dist(:, 1)
      col_dist_block => col_dist(:, 1)
      row_dist_cluster => row_dist(:, 2)
      col_dist_cluster => col_dist(:, 2)

      CALL dbcsr_distribution_new_clusters1(dist, mp_env, &
                                            row_dist_block, col_dist_block, &
                                            row_dist_cluster, col_dist_cluster, &
                                            reuse_arrays)

   END SUBROUTINE dbcsr_distribution_new_clusters2

! **************************************************************************************************
!> \brief Creates new distribution from blockr distributions
!> \param[out] dist           distribution
!> \param[in] mp_env          multiprocessing environment
!> \param row_dist_block ...
!> \param col_dist_block ...
!> \param reuse_arrays ...
! **************************************************************************************************
   SUBROUTINE dbcsr_distribution_new_noclusters(dist, mp_env, row_dist_block, col_dist_block, reuse_arrays)
      TYPE(dbcsr_distribution_obj), INTENT(OUT)          :: dist
      TYPE(dbcsr_mp_obj), INTENT(IN)                     :: mp_env
      INTEGER, DIMENSION(:), INTENT(IN), POINTER         :: row_dist_block, col_dist_block
      LOGICAL, INTENT(IN), OPTIONAL                      :: reuse_arrays

      CALL dbcsr_distribution_new_clusters1(dist, mp_env, &
                                            row_dist_block=row_dist_block, &
                                            col_dist_block=col_dist_block, &
                                            row_dist_cluster=Null(), &
                                            col_dist_cluster=Null(), &
                                            reuse_arrays=reuse_arrays)

   END SUBROUTINE dbcsr_distribution_new_noclusters

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param reshuffle ...
! **************************************************************************************************
   SUBROUTINE dbcsr_finalize(matrix, reshuffle)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: reshuffle

      CALL dbcsr_finalize_prv(matrix%prv, reshuffle)
   END SUBROUTINE dbcsr_finalize

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

      CALL dbcsr_work_create_prv(matrix%prv, nblks_guess, sizedata_guess, n, work_mutable)
   END SUBROUTINE dbcsr_work_create

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param name ...
!> \param dist ...
!> \param matrix_type ...
!> \param row_blk_size ...
!> \param col_blk_size ...
!> \param nze ...
!> \param data_type ...
!> \param reuse ...
!> \param reuse_arrays ...
!> \param mutable_work ...
!> \param replication_type ...
! **************************************************************************************************
   SUBROUTINE dbcsr_create_new(matrix, name, dist, matrix_type, &
                               row_blk_size, col_blk_size, nze, data_type, reuse, &
                               reuse_arrays, mutable_work, replication_type)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      CHARACTER(len=*), INTENT(IN)                       :: name
      TYPE(dbcsr_distribution_obj), INTENT(IN)           :: dist
      CHARACTER, INTENT(IN)                              :: matrix_type
      INTEGER, DIMENSION(:), INTENT(IN), POINTER         :: row_blk_size, col_blk_size
      INTEGER, INTENT(IN), OPTIONAL                      :: nze, data_type
      LOGICAL, INTENT(IN), OPTIONAL                      :: reuse, reuse_arrays, mutable_work
      CHARACTER, INTENT(IN), OPTIONAL                    :: replication_type

      CALL dbcsr_create_prv(matrix%prv, name, dist, &
                            matrix_type, &
                            row_blk_size, col_blk_size, nze=nze, &
                            data_type=data_type, reuse=reuse, &
                            reuse_arrays=reuse_arrays, &
                            mutable_work=mutable_work, replication_type=replication_type)

      matrix%ref_count = 1

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
!> \param nze ...
!> \param data_type ...
!> \param reuse_arrays ...
!> \param mutable_work ...
!> \param replication_type ...
! **************************************************************************************************
   SUBROUTINE dbcsr_create_template(matrix, name, template, &
                                    dist, matrix_type, &
                                    row_blk_size, col_blk_size, nze, data_type, &
                                    reuse_arrays, mutable_work, replication_type)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: name
      TYPE(dbcsr_type), INTENT(IN)                       :: template
      TYPE(dbcsr_distribution_obj), INTENT(IN), OPTIONAL :: dist
      CHARACTER, INTENT(IN), OPTIONAL                    :: matrix_type
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL, &
         POINTER                                         :: row_blk_size, col_blk_size
      INTEGER, INTENT(IN), OPTIONAL                      :: nze, data_type
      LOGICAL, INTENT(IN), OPTIONAL                      :: reuse_arrays, mutable_work
      CHARACTER, INTENT(IN), OPTIONAL                    :: replication_type

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_template', &
         routineP = moduleN//':'//routineN

      CALL dbcsr_create_prv(matrix%prv, template%prv, &
                            name=name, dist=dist, matrix_type=matrix_type, &
                            row_blk_size=row_blk_size, col_blk_size=col_blk_size, &
                            nze=nze, data_type=data_type, &
                            reuse_arrays=reuse_arrays, mutable_work=mutable_work, &
                            replication_type=replication_type)

      matrix%ref_count = 1
   END SUBROUTINE dbcsr_create_template

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eps ...
!> \param method ...
!> \param use_absolute ...
!> \param filter_diag ...
!> \param thorough ...
! **************************************************************************************************
   SUBROUTINE dbcsr_filter(matrix, eps, method, use_absolute, filter_diag, &
                           thorough)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(dp), INTENT(IN)                               :: eps
      INTEGER, INTENT(IN), OPTIONAL                      :: method
      LOGICAL, INTENT(in), OPTIONAL                      :: use_absolute, filter_diag, thorough

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_filter', routineP = moduleN//':'//routineN

      LOGICAL                                            :: quick

      IF (PRESENT(thorough)) THEN
         quick = .NOT. thorough
      ELSE
         quick = .FALSE.
      ENDIF
      CALL dbcsr_filter_prv(matrix%prv, dbcsr_conform_scalar(eps, matrix), &
                            method, use_absolute, filter_diag, &
                            quick=quick)
   END SUBROUTINE dbcsr_filter

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
   SUBROUTINE dbcsr_set_diag(matrix, diag)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: diag

      CALL dbcsr_set_diag_prv(matrix%prv, diag)
   END SUBROUTINE dbcsr_set_diag

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
   SUBROUTINE dbcsr_get_diag(matrix, diag)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      REAL(dp), DIMENSION(:), INTENT(INOUT)              :: diag

      CALL dbcsr_get_diag_prv(matrix%prv, diag)
   END SUBROUTINE dbcsr_get_diag

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param diag ...
! **************************************************************************************************
   SUBROUTINE dbcsr_get_block_diag(matrix, diag)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(dbcsr_type), INTENT(INOUT)                    :: diag

      CALL dbcsr_create(diag, "Diagonal of "//TRIM(dbcsr_name(matrix)), template=matrix)
      CALL dbcsr_get_block_diag_prv(matrix%prv, diag%prv)
   END SUBROUTINE dbcsr_get_block_diag

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param alpha_scalar ...
!> \param first_row ...
!> \param last_row ...
! **************************************************************************************************
   SUBROUTINE dbcsr_add_on_diag(matrix, alpha_scalar, first_row, last_row)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), INTENT(IN)                          :: alpha_scalar
      INTEGER, INTENT(in), OPTIONAL                      :: first_row, last_row

      CALL dbcsr_add_on_diag_prv(matrix%prv, &
                                 dbcsr_conform_scalar(alpha_scalar, matrix), &
                                 first_row, last_row)
   END SUBROUTINE dbcsr_add_on_diag

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param filepath ...
! **************************************************************************************************
   SUBROUTINE dbcsr_binary_write(matrix, filepath)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      CHARACTER(LEN=*), INTENT(IN)                       :: filepath

      CALL dbcsr_binary_write_prv(matrix%prv, filepath)
   END SUBROUTINE dbcsr_binary_write

! **************************************************************************************************
!> \brief ...
!> \param filepath ...
!> \param distribution ...
!> \param groupid ...
!> \param matrix_new ...
! **************************************************************************************************
   SUBROUTINE dbcsr_binary_read(filepath, distribution, groupid, matrix_new)
      CHARACTER(len=*), INTENT(IN)                       :: filepath
      TYPE(dbcsr_distribution_obj), INTENT(IN)           :: distribution
      INTEGER, INTENT(IN), OPTIONAL                      :: groupid
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_new

      CALL dbcsr_binary_read_prv(filepath, distribution, groupid, matrix_new%prv)
   END SUBROUTINE dbcsr_binary_read

! **************************************************************************************************
!> \brief ...
!> \param matrix_b ...
!> \param matrix_a ...
!> \param name ...
!> \param keep_sparsity ...
!> \param shalprv_data ...
!> \param keep_imaginary ...
!> \param matrix_type ...
! **************************************************************************************************
   SUBROUTINE dbcsr_copy(matrix_b, matrix_a, name, keep_sparsity, &
                         shalprv_data, keep_imaginary, matrix_type)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_b
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: name
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_sparsity, shalprv_data, &
                                                            keep_imaginary
      CHARACTER, INTENT(IN), OPTIONAL                    :: matrix_type

      IF (matrix_b%ref_count .EQ. 0) THEN
         CALL dbcsr_create(matrix_b, template=matrix_a)
      ENDIF
      CALL dbcsr_copy_prv(matrix_b%prv, matrix_a%prv, name, keep_sparsity, &
                          shalprv_data, keep_imaginary, matrix_type)
   END SUBROUTINE dbcsr_copy

! **************************************************************************************************
!> \brief ...
!> \param matrix_b ...
!> \param matrix_a ...
! **************************************************************************************************
   SUBROUTINE dbcsr_copy_into_existing(matrix_b, matrix_a)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_b
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a

      IF (matrix_b%ref_count .EQ. 0) THEN
         CALL dbcsr_create(matrix_b, template=matrix_a)
      ENDIF
      CALL dbcsr_copy_into_existing_prv(matrix_b%prv, matrix_a%prv)
   END SUBROUTINE dbcsr_copy_into_existing

! **************************************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
! **************************************************************************************************
   SUBROUTINE dbcsr_desymmetrize(matrix_a, matrix_b)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_b

      IF (matrix_b%ref_count .EQ. 0) THEN
         CALL dbcsr_create(matrix_b, template=matrix_a)
      ENDIF
      CALL dbcsr_desymmetrize_deep_prv(matrix_a%prv, matrix_b%prv, untransposed_data=.TRUE.)
   END SUBROUTINE dbcsr_desymmetrize

! **************************************************************************************************
!> \brief ...
!> \param transposed ...
!> \param normal ...
!> \param shallow_data_copy ...
!> \param transpose_data ...
!> \param transpose_distribution ...
!> \param use_distribution ...
! **************************************************************************************************
   SUBROUTINE dbcsr_transposed(transposed, normal, shallow_data_copy, &
                               transpose_data, transpose_distribution, use_distribution)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: transposed
      TYPE(dbcsr_type), INTENT(IN)                       :: normal
      LOGICAL, INTENT(IN), OPTIONAL                      :: shallow_data_copy, transpose_data, &
                                                            transpose_distribution
      TYPE(dbcsr_distribution_obj), INTENT(IN), OPTIONAL :: use_distribution

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_transposed', &
         routineP = moduleN//':'//routineN

      LOGICAL                                            :: myshallow_data_copy, &
                                                            mytranspose_distribution
      TYPE(dbcsr_distribution_obj)                       :: myuse_distribution

!   set some defaults to make usage a bit less painful (fschiff)

      myshallow_data_copy = .FALSE.
      myuse_distribution = dbcsr_distribution(normal)
      mytranspose_distribution = .FALSE.
      IF (PRESENT(shallow_data_copy)) myshallow_data_copy = shallow_data_copy
      IF (PRESENT(use_distribution)) myuse_distribution = use_distribution
      IF (PRESENT(transpose_distribution)) mytranspose_distribution = transpose_distribution

      CALL dbcsr_new_transposed_prv(transposed%prv, normal%prv, myshallow_data_copy, &
                                    transpose_data, mytranspose_distribution, &
                                    use_distribution=myuse_distribution)
   END SUBROUTINE dbcsr_transposed

! **************************************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param func ...
!> \param a0 ...
!> \param a1 ...
!> \param a2 ...
! **************************************************************************************************
   SUBROUTINE dbcsr_function_of_elements(matrix_a, func, a0, a1, a2)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_a
      INTEGER, INTENT(IN)                                :: func
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: a0, a1, a2

      CALL dbcsr_function_of_elements_prv(matrix_a%prv, func, a0, a1, a2)
   END SUBROUTINE dbcsr_function_of_elements

! **************************************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param matrix_c ...
!> \param b_assume_value ...
! **************************************************************************************************
   SUBROUTINE dbcsr_hadamard_product(matrix_a, matrix_b, matrix_c, b_assume_value)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_a, matrix_b
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_c
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: b_assume_value

      CALL dbcsr_hadamard_product_prv(matrix_a%prv, matrix_b%prv, matrix_c%prv, b_assume_value)
   END SUBROUTINE dbcsr_hadamard_product

! **************************************************************************************************
!> \brief     Adjusts matrices to be compatibly for multiplication.
!>
!> In BLAS xGEMM or PxGEMM style, CP2K alprvs with differing
!> dimensions to be multiplied.  This compatibilyt has been removed
!> from DBCSR.  Therefore, matrices must be matched for
!> multiplication:
!> <ul>
!> <li> The rows of matrix_b are adjusted to fit the columns of matrix_a,
!> <li> The columns of matrix_b are adjusted to fit the columns of matrix_c, and
!> <li> The rows of matrix_a are adjusted to fit the rows of matrix_c.
!> </ul>
!> Transformations are skipped if the sizes are already compatible.
!> However, in any case a copy of matrix_a and matrix_b is made.
!> It is the responsibility of the caller to release these new matrices.
!>
!> Some copies could be avoided with more code.
!> \param matrix_c ...
!> \param matrix_a ...
!> \param tr_a ...
!> \param matrix_b ...
!> \param tr_b ...
!> \param new_a ...
!> \param new_b ...
! **************************************************************************************************
   SUBROUTINE matrix_match_sizes(matrix_c, matrix_a, tr_a, matrix_b, tr_b, &
                                 new_a, new_b)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_c, matrix_a
      CHARACTER, INTENT(IN)                              :: tr_a
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_b
      CHARACTER, INTENT(IN)                              :: tr_b
      TYPE(dbcsr_type), INTENT(OUT)                      :: new_a, new_b

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_match_sizes', &
         routineP = moduleN//':'//routineN

      CHARACTER                                          :: tr_a_l, tr_b_l
      INTEGER                                            :: a_cols_total, a_rows_total, &
                                                            b_cols_total, b_rows_total, &
                                                            c_cols_total, c_rows_total, tmp
      INTEGER, DIMENSION(:), POINTER                     :: new_col_clusters, new_col_dist, &
                                                            new_col_size, new_row_clusters, &
                                                            new_row_dist, new_row_size
      LOGICAL                                            :: atr, btr, create_col_dist, &
                                                            create_row_dist, new_a_rows, &
                                                            new_b_cols, new_b_rows
      TYPE(dbcsr_distribution_obj)                       :: a_dist, b_dist, c_dist, new_dist

!   ---------------------------------------------------------------------------

      tr_a_l = tr_a
      CALL uppercase(tr_a_l)
      tr_b_l = tr_b
      CALL uppercase(tr_b_l)

      create_col_dist = .FALSE.
      create_row_dist = .FALSE.

      btr = tr_b_l .NE. dbcsr_no_transpose
      atr = tr_a_l .NE. dbcsr_no_transpose
      !
      c_rows_total = dbcsr_nfullrows_total(matrix_c)
      c_cols_total = dbcsr_nfullcols_total(matrix_c)
      b_rows_total = dbcsr_nfullrows_total(matrix_b)
      b_cols_total = dbcsr_nfullcols_total(matrix_b)
      a_rows_total = dbcsr_nfullrows_total(matrix_a)
      a_cols_total = dbcsr_nfullcols_total(matrix_a)
      IF (atr) THEN
         tmp = a_cols_total
         a_cols_total = a_rows_total
         a_rows_total = tmp
      ENDIF
      IF (btr) THEN
         tmp = b_cols_total
         b_cols_total = b_rows_total
         b_rows_total = tmp
      ENDIF

      NULLIFY (new_row_clusters, new_col_clusters)

      !
      ! Process matrix B.
      !
      ! This check is faster than explicity comparing blocked row and
      ! column sizes.
      new_b_cols = c_cols_total .NE. b_cols_total
      new_b_rows = a_cols_total .NE. b_rows_total
      a_dist = dbcsr_distribution(matrix_a)
      b_dist = dbcsr_distribution(matrix_b)
      c_dist = dbcsr_distribution(matrix_c)
      IF (new_b_rows .OR. new_b_cols) THEN
         CALL dbcsr_init(new_b)
         IF (.NOT. btr) THEN
            IF (new_b_cols) THEN
               CALL match_1_dist(new_col_dist, &
                                 dbcsr_distribution_col_dist(b_dist), &
                                 dbcsr_distribution_col_dist(c_dist), &
                                 new_col_clusters, &
                                 dbcsr_distribution_col_clusters(b_dist), &
                                 dbcsr_distribution_col_clusters(c_dist))
               create_col_dist = .TRUE.
            ELSE
               new_col_dist => dbcsr_distribution_col_dist(b_dist)
               new_col_clusters => dbcsr_distribution_col_clusters(b_dist)
            ENDIF
            new_col_size => dbcsr_col_block_sizes(matrix_c)
            IF (new_b_rows) THEN
               IF (.NOT. atr) THEN
                  CALL match_1_dist(new_row_dist, &
                                    dbcsr_distribution_row_dist(b_dist), &
                                    dbcsr_distribution_col_dist(a_dist), &
                                    new_row_clusters, &
                                    dbcsr_distribution_row_clusters(b_dist), &
                                    dbcsr_distribution_col_clusters(a_dist))
                  create_row_dist = .TRUE.
                  new_row_size => dbcsr_col_block_sizes(matrix_a)
               ELSE
                  CALL match_1_dist(new_row_dist, &
                                    dbcsr_distribution_row_dist(b_dist), &
                                    dbcsr_distribution_row_dist(a_dist), &
                                    new_row_clusters, &
                                    dbcsr_distribution_row_clusters(b_dist), &
                                    dbcsr_distribution_row_clusters(a_dist))
                  create_row_dist = .TRUE.
                  new_row_size => dbcsr_row_block_sizes(matrix_a)
               ENDIF
            ELSE
               new_row_dist => dbcsr_distribution_row_dist(b_dist)
               new_row_clusters => dbcsr_distribution_row_clusters(b_dist)
               new_row_size => dbcsr_row_block_sizes(matrix_b)
            ENDIF
         ELSE
            IF (new_b_cols) THEN
               CALL match_1_dist(new_row_dist, &
                                 dbcsr_distribution_row_dist(b_dist), &
                                 dbcsr_distribution_col_dist(c_dist), &
                                 new_row_clusters, &
                                 dbcsr_distribution_row_clusters(b_dist), &
                                 dbcsr_distribution_col_clusters(c_dist))
               create_row_dist = .TRUE.
            ELSE
               new_row_dist => dbcsr_distribution_row_dist(b_dist)
               new_row_clusters => dbcsr_distribution_row_clusters(b_dist)
            ENDIF
            new_row_size => dbcsr_col_block_sizes(matrix_c)
            IF (new_b_rows) THEN
               IF (.NOT. atr) THEN
                  CALL match_1_dist(new_col_dist, &
                                    dbcsr_distribution_col_dist(b_dist), &
                                    dbcsr_distribution_col_dist(a_dist), &
                                    new_col_clusters, &
                                    dbcsr_distribution_col_clusters(b_dist), &
                                    dbcsr_distribution_col_clusters(a_dist))
                  create_col_dist = .TRUE.
                  new_col_size => dbcsr_col_block_sizes(matrix_a)
               ELSE
                  CALL match_1_dist(new_col_dist, &
                                    dbcsr_distribution_col_dist(b_dist), &
                                    dbcsr_distribution_row_dist(a_dist), &
                                    new_col_clusters, &
                                    dbcsr_distribution_col_clusters(b_dist), &
                                    dbcsr_distribution_row_clusters(a_dist))
                  create_col_dist = .TRUE.
                  new_col_size => dbcsr_row_block_sizes(matrix_a)
               ENDIF
            ELSE
               new_col_dist => dbcsr_distribution_col_dist(b_dist)
               new_col_clusters => dbcsr_distribution_col_clusters(b_dist)
               new_col_size => dbcsr_col_block_sizes(matrix_b)
            ENDIF
         ENDIF
         CALL dbcsr_distribution_new(new_dist, dbcsr_distribution_mp(c_dist), &
                                     new_row_dist, new_col_dist, &
                                     row_dist_cluster=new_row_clusters, &
                                     col_dist_cluster=new_col_clusters)
         IF (create_col_dist) THEN
            DEALLOCATE (new_col_dist)
            IF (ASSOCIATED(new_col_clusters)) DEALLOCATE (new_col_clusters)
         ENDIF
         IF (create_row_dist) THEN
            DEALLOCATE (new_row_dist)
            IF (ASSOCIATED(new_row_clusters)) DEALLOCATE (new_row_clusters)
         ENDIF
         CALL dbcsr_create(new_b, template=matrix_b, &
                           dist=new_dist, &
                           row_blk_size=new_row_size, &
                           col_blk_size=new_col_size)
         CALL dbcsr_distribution_release(new_dist)
         CALL dbcsr_complete_redistribute(matrix_b, new_b)
      ELSE
         CALL dbcsr_init(new_b)
         CALL dbcsr_copy(new_b, matrix_b, shalprv_data=.TRUE.)
      ENDIF

      create_col_dist = .FALSE.
      create_row_dist = .FALSE.
      NULLIFY (new_row_clusters, new_col_clusters)
      !
      ! Process matrix A
      new_a_rows = a_rows_total .NE. c_rows_total
      IF (new_a_rows) THEN
         CALL dbcsr_init(new_a)
         IF (atr) THEN
            new_row_dist => dbcsr_distribution_row_dist(a_dist)
            new_row_clusters => dbcsr_distribution_row_clusters(a_dist)
            new_row_size => dbcsr_row_block_sizes(matrix_a)
            CALL match_1_dist(new_col_dist, &
                              dbcsr_distribution_col_dist(a_dist), &
                              dbcsr_distribution_row_dist(c_dist), &
                              new_col_clusters, &
                              dbcsr_distribution_col_clusters(a_dist), &
                              dbcsr_distribution_row_clusters(c_dist))
            create_col_dist = .TRUE.
            new_col_size => dbcsr_row_block_sizes(matrix_c)
         ELSE
            CALL match_1_dist(new_row_dist, &
                              dbcsr_distribution_row_dist(a_dist), &
                              dbcsr_distribution_row_dist(c_dist), &
                              new_row_clusters, &
                              dbcsr_distribution_row_clusters(a_dist), &
                              dbcsr_distribution_row_clusters(c_dist))
            create_row_dist = .TRUE.
            new_row_size => dbcsr_row_block_sizes(matrix_c)
            new_col_dist => dbcsr_distribution_col_dist(a_dist)
            new_col_clusters => dbcsr_distribution_col_clusters(a_dist)
            new_col_size => dbcsr_col_block_sizes(matrix_a)
         ENDIF
         CALL dbcsr_distribution_new(new_dist, dbcsr_distribution_mp(a_dist), &
                                     new_row_dist, new_col_dist, &
                                     new_row_clusters, new_col_clusters)
         IF (create_col_dist) THEN
            DEALLOCATE (new_col_dist)
            DEALLOCATE (new_col_clusters)
         ENDIF
         IF (create_row_dist) THEN
            DEALLOCATE (new_row_dist)
            DEALLOCATE (new_row_clusters)
         ENDIF
         CALL dbcsr_create(new_a, template=matrix_a, &
                           dist=new_dist, &
                           row_blk_size=new_row_size, &
                           col_blk_size=new_col_size)
         CALL dbcsr_distribution_release(new_dist)
         CALL dbcsr_complete_redistribute(matrix_a, new_a)
      ELSE
         CALL dbcsr_init(new_a)
         CALL dbcsr_copy(new_a, matrix_a, shalprv_data=.TRUE.)
      ENDIF
   END SUBROUTINE matrix_match_sizes

! **************************************************************************************************
!> \brief Matches a 1-D distribution
!>
!> Forms a new distribution new_dist of size template_dist.  As much data as
!> possible is taken from old_dist.
!> \param new_dist ...
!> \param old_dist ...
!> \param template_dist ...
!> \param new_clusters ...
!> \param old_clusters ...
!> \param template_clusters ...
! **************************************************************************************************
   SUBROUTINE match_1_dist(new_dist, old_dist, template_dist, &
                           new_clusters, old_clusters, template_clusters)
      INTEGER, DIMENSION(:), INTENT(out), POINTER        :: new_dist
      INTEGER, DIMENSION(:), INTENT(in), POINTER         :: old_dist, template_dist
      INTEGER, DIMENSION(:), INTENT(out), POINTER        :: new_clusters
      INTEGER, DIMENSION(:), INTENT(in), POINTER         :: old_clusters, template_clusters

      INTEGER                                            :: i, max_bin, new_size, old_size, sz

      old_size = SIZE(old_dist)
      new_size = SIZE(template_dist)
      sz = MIN(old_size, new_size)
      ALLOCATE (new_dist(new_size))
      new_dist = template_dist
      new_dist(1:sz) = old_dist(1:sz)
      max_bin = MAXVAL(old_dist)
      IF (max_bin .GT. 0) THEN
         DO i = sz+1, new_size
            new_dist(i) = MOD(i, max_bin)
         END DO
      ELSE
         new_dist(sz+1:new_size) = 0
      ENDIF
      IF (ASSOCIATED(template_clusters) .AND. ASSOCIATED(old_clusters)) THEN
         ALLOCATE (new_clusters(new_size))
         new_clusters = template_clusters
         new_clusters(1:sz) = old_clusters(1:sz)
         max_bin = MAXVAL(old_clusters)
         DO i = sz+1, new_size
            new_clusters(i) = MOD(i, max_bin)
         END DO
      ENDIF
   END SUBROUTINE match_1_dist

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval row_blk_sizes ...
! **************************************************************************************************
   FUNCTION dbcsr_row_block_sizes(matrix) RESULT(row_blk_sizes)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes

      CALL dbcsr_get_info(matrix, row_blk_size=row_blk_sizes)
   END FUNCTION dbcsr_row_block_sizes

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval col_blk_sizes ...
! **************************************************************************************************
   FUNCTION dbcsr_col_block_sizes(matrix) RESULT(col_blk_sizes)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes

      CALL dbcsr_get_info(matrix, col_blk_size=col_blk_sizes)
   END FUNCTION dbcsr_col_block_sizes

! **************************************************************************************************
!> \brief Deallocates a DBCSR matrix for compatibility with CP2K
!> \param[in,out] matrix      DBCSR matrix
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
!> \brief   Allocate and initialize a real matrix 1-dimensional set.
!> \param[in,out] matrix_set  Set containing the DBCSR matrices
!> \param[in] nmatrix         Size of set
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! **************************************************************************************************
   SUBROUTINE allocate_dbcsr_matrix_set(matrix_set, nmatrix)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_set
      INTEGER, INTENT(IN)                                :: nmatrix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_dbcsr_matrix_set', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: imatrix

      IF (ASSOCIATED(matrix_set)) CALL dbcsr_deallocate_matrix_set(matrix_set)
      ALLOCATE (matrix_set(nmatrix))
      DO imatrix = 1, nmatrix
         NULLIFY (matrix_set(imatrix)%matrix)
      END DO
   END SUBROUTINE allocate_dbcsr_matrix_set

! **************************************************************************************************
!> \brief   Allocate and initialize a real matrix 2-dimensional set.
!> \param[in,out] matrix_set  Set containing the DBCSR matrix pointer type
!> \param[in] nmatrix         Size of set
!> \param mmatrix ...
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! **************************************************************************************************
   SUBROUTINE allocate_dbcsr_matrix_set_2d(matrix_set, nmatrix, mmatrix)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_set
      INTEGER, INTENT(IN)                                :: nmatrix, mmatrix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_dbcsr_matrix_set_2d', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: imatrix, jmatrix

      IF (ASSOCIATED(matrix_set)) CALL dbcsr_deallocate_matrix_set(matrix_set)
      ALLOCATE (matrix_set(nmatrix, mmatrix))
      DO jmatrix = 1, mmatrix
         DO imatrix = 1, nmatrix
            NULLIFY (matrix_set(imatrix, jmatrix)%matrix)
         END DO
      END DO
   END SUBROUTINE allocate_dbcsr_matrix_set_2d

! **************************************************************************************************
!> \brief   Allocate and initialize a real matrix 2-dimensional set.
!> \param[in,out] matrix_set  Set containing the DBCSR matrix pointer type
!> \param[in] nmatrix         Size of set
!> \param mmatrix ...
!> \param pmatrix ...
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! **************************************************************************************************
   SUBROUTINE allocate_dbcsr_matrix_set_3d(matrix_set, nmatrix, mmatrix, pmatrix)
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: matrix_set
      INTEGER, INTENT(IN)                                :: nmatrix, mmatrix, pmatrix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_dbcsr_matrix_set_3d', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: imatrix, jmatrix, kmatrix

      IF (ASSOCIATED(matrix_set)) CALL dbcsr_deallocate_matrix_set(matrix_set)
      ALLOCATE (matrix_set(nmatrix, mmatrix, pmatrix))
      DO kmatrix = 1, pmatrix
         DO jmatrix = 1, mmatrix
            DO imatrix = 1, nmatrix
               NULLIFY (matrix_set(imatrix, jmatrix, kmatrix)%matrix)
            END DO
         END DO
      END DO
   END SUBROUTINE allocate_dbcsr_matrix_set_3d

! **************************************************************************************************
!> \brief Deallocate a real matrix set and release all of the member matrices.
!> \param[in,out] matrix_set  Set containing the DBCSR matrix pointer type
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! **************************************************************************************************
   SUBROUTINE deallocate_dbcsr_matrix_set(matrix_set)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_dbcsr_matrix_set', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: imatrix

      IF (ASSOCIATED(matrix_set)) THEN
         DO imatrix = 1, SIZE(matrix_set)
            CALL dbcsr_deallocate_matrix(matrix_set(imatrix)%matrix)
         END DO
         DEALLOCATE (matrix_set)
      END IF

   END SUBROUTINE deallocate_dbcsr_matrix_set

! **************************************************************************************************
!> \brief Deallocate a real matrix set and release all of the member matrices.
!> \param[in,out] matrix_set  Set containing the DBCSR matrix pointer type
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! **************************************************************************************************
   SUBROUTINE deallocate_dbcsr_matrix_set_2d(matrix_set)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_dbcsr_matrix_set_2d', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: imatrix, jmatrix

      IF (ASSOCIATED(matrix_set)) THEN
         DO jmatrix = 1, SIZE(matrix_set, 2)
            DO imatrix = 1, SIZE(matrix_set, 1)
               CALL dbcsr_deallocate_matrix(matrix_set(imatrix, jmatrix)%matrix)
            END DO
         END DO
         DEALLOCATE (matrix_set)
      END IF
   END SUBROUTINE deallocate_dbcsr_matrix_set_2d

! **************************************************************************************************
!> \brief Deallocate a real matrix set and release all of the member matrices.
!> \param[in,out] matrix_set  Set containing the DBCSR matrix pointer type
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! **************************************************************************************************
   SUBROUTINE deallocate_dbcsr_matrix_set_3d(matrix_set)

      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: matrix_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_dbcsr_matrix_set_3d', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: imatrix, jmatrix, kmatrix

      IF (ASSOCIATED(matrix_set)) THEN
         DO kmatrix = 1, SIZE(matrix_set, 3)
            DO jmatrix = 1, SIZE(matrix_set, 2)
               DO imatrix = 1, SIZE(matrix_set, 1)
                  CALL dbcsr_deallocate_matrix(matrix_set(imatrix, jmatrix, kmatrix)%matrix)
               END DO
            END DO
         END DO
         DEALLOCATE (matrix_set)
      END IF
   END SUBROUTINE deallocate_dbcsr_matrix_set_3d

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval has_symmetry ...
! **************************************************************************************************
   PURE FUNCTION dbcsr_has_symmetry(matrix) RESULT(has_symmetry)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      LOGICAL                                            :: has_symmetry

      has_symmetry = dbcsr_has_symmetry_prv(matrix%prv)
   END FUNCTION dbcsr_has_symmetry

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_mat ...
!> \param csr_mat ...
!> \param dist_format ...
!> \param csr_sparsity ...
!> \param numnodes ...
! **************************************************************************************************
   SUBROUTINE csr_create_from_dbcsr(dbcsr_mat, csr_mat, dist_format, csr_sparsity, numnodes)

      TYPE(dbcsr_type), INTENT(IN)                       :: dbcsr_mat
      TYPE(csr_type), INTENT(OUT)                        :: csr_mat
      INTEGER                                            :: dist_format
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: csr_sparsity
      INTEGER, INTENT(IN), OPTIONAL                      :: numnodes

      IF (PRESENT(csr_sparsity)) THEN
         CALL csr_create_from_dbcsr_prv(dbcsr_mat%prv, csr_mat, dist_format, csr_sparsity%prv, numnodes)
      ELSE
         CALL csr_create_from_dbcsr_prv(dbcsr_mat%prv, csr_mat, dist_format, numnodes=numnodes)
      ENDIF
   END SUBROUTINE csr_create_from_dbcsr

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_mat ...
!> \param csr_mat ...
! **************************************************************************************************
   SUBROUTINE convert_csr_to_dbcsr(dbcsr_mat, csr_mat)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: dbcsr_mat
      TYPE(csr_type), INTENT(INOUT)                      :: csr_mat

      CALL convert_csr_to_dbcsr_prv(dbcsr_mat%prv, csr_mat)

   END SUBROUTINE convert_csr_to_dbcsr

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_mat ...
!> \param csr_mat ...
! **************************************************************************************************
   SUBROUTINE convert_dbcsr_to_csr(dbcsr_mat, csr_mat)
      TYPE(dbcsr_type), INTENT(IN)                       :: dbcsr_mat
      TYPE(csr_type), INTENT(INOUT)                      :: csr_mat

      CALL convert_dbcsr_to_csr_prv(dbcsr_mat%prv, csr_mat)
   END SUBROUTINE convert_dbcsr_to_csr

! **************************************************************************************************
!> \brief Apply filtering threshold eps to DBCSR blocks in order to improve
!>        CSR sparsity (currently only used for testing purposes)
!> \param dbcsr_mat ...
!> \param csr_sparsity ...
!> \param eps ...
! **************************************************************************************************
   SUBROUTINE dbcsr_to_csr_filter(dbcsr_mat, csr_sparsity, eps)
      TYPE(dbcsr_type), INTENT(IN)                       :: dbcsr_mat
      TYPE(dbcsr_type), INTENT(OUT)                      :: csr_sparsity
      REAL(kind=real_8), INTENT(IN)                      :: eps

      CALL dbcsr_to_csr_filter_prv(dbcsr_mat%prv, csr_sparsity%prv, eps)
   END SUBROUTINE dbcsr_to_csr_filter

! **************************************************************************************************
!> \brief Emulation of sparse_matrix_types/add_block_node mapped
!>        to add_real_matrix_block.... should not be used any longer
!>
!> It adds a block to the dbcsr matrix and returns a rank-2 pointer to the
!> block. Currently it only and always uses the mutable data.
!> \param[in,out] matrix      DBCSR matrix
!> \param[in]  block_row      the row
!> \param[in]  block_col      the column
!> \param[in]  block          the block to put
! **************************************************************************************************
   SUBROUTINE dbcsr_add_block_node(matrix, block_row, block_col, block)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, INTENT(IN)                                :: block_row, block_col
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_add_block_node', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: c, ithread, p, r
      LOGICAL                                            :: dbg, existed, is_there, tr

!   ---------------------------------------------------------------------------

      dbg = .FALSE.

      ithread = 0
!$    ithread = omp_get_thread_num()
      IF (.NOT. ASSOCIATED(matrix%prv%wms)) THEN
         CALL dbcsr_work_create(matrix, work_mutable=.TRUE.)
         matrix%prv%valid = .FALSE.
      ENDIF
!$    IF (SIZE(matrix%prv%wms) .LT. omp_get_num_threads()) &
!$       CPABORT("Too few threads.")
      IF (.NOT. dbcsr_wm_use_mutable(matrix%prv%wms(ithread+1))) &
         CPABORT("Data loss due to no conversion of appendable to mutable data")
      is_there = ASSOCIATED(block)
      !r = row ; c = col ; tr = .FALSE.
      !CALL dbcsr_get_stored_coordinates (matrix, r, c, tr)
      !CALL dbcsr_reserve_block2d (matrix, row, col, block)
      !write(*,*) 'add_block_node: block_row',block_row,' block_col',block_col
      CALL dbcsr_reserve_block2d(matrix, block_row, block_col, block, &
                                 existed=existed)
!
      IF (dbg) THEN
         r = block_row; c = block_col; tr = .FALSE.
         CALL dbcsr_get_stored_coordinates(matrix, r, c, p)
         IF (p .NE. dbcsr_mp_mynode(dbcsr_distribution_mp(dbcsr_distribution(matrix)))) &
            CPWARN("Adding non-local element")
      ENDIF
      IF (existed) CPWARN("You should not add existing blocks according to old API.")
      IF (.NOT. is_there) block(:, :) = 0.0_dp
   END SUBROUTINE dbcsr_add_block_node

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_vec ...
!> \param matrix ...
!> \param ncol ...
! **************************************************************************************************
   SUBROUTINE create_col_vec_from_matrix(dbcsr_vec, matrix, ncol)
      TYPE(dbcsr_type)                                   :: dbcsr_vec, matrix
      INTEGER                                            :: ncol

      CALL create_col_vec_from_matrix_prv(dbcsr_vec%prv, matrix%prv, ncol)
   END SUBROUTINE create_col_vec_from_matrix

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_vec ...
!> \param matrix ...
!> \param nrow ...
! **************************************************************************************************
   SUBROUTINE create_row_vec_from_matrix(dbcsr_vec, matrix, nrow)
      TYPE(dbcsr_type)                                   :: dbcsr_vec, matrix
      INTEGER                                            :: nrow

      CALL create_row_vec_from_matrix_prv(dbcsr_vec%prv, matrix%prv, nrow)
   END SUBROUTINE create_row_vec_from_matrix

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_vec ...
!> \param matrix ...
!> \param ncol ...
! **************************************************************************************************
   SUBROUTINE create_replicated_col_vec_from_matrix(dbcsr_vec, matrix, ncol)
      TYPE(dbcsr_type)                                   :: dbcsr_vec, matrix
      INTEGER                                            :: ncol

      CALL create_replicated_col_vec_from_matrix_prv(dbcsr_vec%prv, matrix%prv, ncol)
   END SUBROUTINE create_replicated_col_vec_from_matrix

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_vec ...
!> \param matrix ...
!> \param nrow ...
! **************************************************************************************************
   SUBROUTINE create_replicated_row_vec_from_matrix(dbcsr_vec, matrix, nrow)
      TYPE(dbcsr_type)                                   :: dbcsr_vec, matrix
      INTEGER                                            :: nrow

      CALL create_replicated_row_vec_from_matrix_prv(dbcsr_vec%prv, matrix%prv, nrow)
   END SUBROUTINE create_replicated_row_vec_from_matrix

#include "dbcsr_api_s.f90"
#include "dbcsr_api_d.f90"
#include "dbcsr_api_c.f90"
#include "dbcsr_api_z.f90"

END MODULE dbcsr_api

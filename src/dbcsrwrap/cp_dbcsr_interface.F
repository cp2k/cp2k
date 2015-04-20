!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR interface in CP2K
!> \author  VE
!> \date    2010
!> \version 0.1
!>
!> <b>Modification history:</b>
!> - Created 2010
! *****************************************************************************
MODULE cp_dbcsr_interface
  USE bibliography,                    ONLY: Borstnik2014,&
                                             cite_reference
  USE dbcsr_api,                       ONLY: &
       accdrv_avoid_after_busy, accdrv_binning_binsize, accdrv_binning_nbins, &
       accdrv_do_inhomogenous, accdrv_min_flop_process, accdrv_min_flop_sort, &
       accdrv_posterior_buffers, accdrv_posterior_streams, &
       accdrv_priority_buffers, accdrv_priority_streams, add_work_coordinate, &
       array_data, convert_csr_to_dbcsr, convert_dbcsr_to_csr, &
       convert_offsets_to_sizes, convert_sizes_to_offsets, &
       cp_dbcsr_iterator=>dbcsr_iterator, &
       cp_distribution_get_num_images=>dbcsr_distribution_get_num_images_1d, &
       create_bl_distribution, create_col_vec_from_matrix, &
       create_replicated_col_vec_from_matrix, &
       create_replicated_row_vec_from_matrix, create_row_vec_from_matrix, &
       csr_create_from_dbcsr, csr_create_new, csr_create_template, &
       csr_dbcsr_blkrow_dist, csr_destroy, csr_eqrow_floor_dist, &
       csr_print_sparsity, csr_type, dbcsr_add, dbcsr_add_on_diag, &
       dbcsr_arnoldi_data, dbcsr_arnoldi_ev, dbcsr_binary_read, &
       dbcsr_binary_write, dbcsr_checksum, dbcsr_clear_mempools, &
       dbcsr_complete_redistribute, dbcsr_copy, dbcsr_copy_into_existing, &
       dbcsr_create, dbcsr_create_dist_block_cyclic, &
       dbcsr_create_dist_r_unrot, dbcsr_data_clear_pointer, dbcsr_data_init, &
       dbcsr_data_new, dbcsr_data_obj, dbcsr_data_release, &
       dbcsr_desymmetrize_deep, dbcsr_distribute, dbcsr_distribution, &
       dbcsr_distribution_col_dist, dbcsr_distribution_has_threads, &
       dbcsr_distribution_init, dbcsr_distribution_local_cols, &
       dbcsr_distribution_local_rows, dbcsr_distribution_mp, &
       dbcsr_distribution_ncols, dbcsr_distribution_new, &
       dbcsr_distribution_nlocal_cols, dbcsr_distribution_nlocal_rows, &
       dbcsr_distribution_nrows, dbcsr_distribution_obj, &
       dbcsr_distribution_release, dbcsr_distribution_row_dist, &
       dbcsr_error_set, dbcsr_error_stop, dbcsr_error_type, dbcsr_filter, &
       dbcsr_finalize, dbcsr_finalize_lib, dbcsr_frobenius_norm, &
       dbcsr_func_artanh, dbcsr_func_dtanh, dbcsr_func_inverse, &
       dbcsr_func_tanh, dbcsr_function_of_elements, dbcsr_gershgorin_norm, &
       dbcsr_get_block_diag, dbcsr_get_block_p, &
       dbcsr_get_conf_avg_elements_images, dbcsr_get_conf_comm_thread_load, &
       dbcsr_get_conf_max_ele_block, dbcsr_get_conf_mm_driver, &
       dbcsr_get_conf_mm_stacksize, dbcsr_get_conf_mpi_mem, &
       dbcsr_get_conf_nstacks, dbcsr_get_conf_subcomm, &
       dbcsr_get_conf_use_comm_thread, dbcsr_get_conf_use_mpi_rma, &
       dbcsr_get_data_p, dbcsr_get_data_size, dbcsr_get_data_type, &
       dbcsr_get_diag, dbcsr_get_info, dbcsr_get_matrix_type, &
       dbcsr_get_num_blocks, dbcsr_get_occupation, &
       dbcsr_get_stored_coordinates, dbcsr_hadamard_product, &
       dbcsr_has_symmetry, dbcsr_init, dbcsr_init_lib, dbcsr_init_random, &
       dbcsr_iterator, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
       dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_maxabs, &
       dbcsr_mp_group, dbcsr_mp_hold, dbcsr_mp_mynode, dbcsr_mp_new, &
       dbcsr_mp_npcols, dbcsr_mp_nprows, dbcsr_mp_numnodes, dbcsr_mp_obj, &
       dbcsr_mp_release, dbcsr_multiply, dbcsr_name, dbcsr_nblkcols_local, &
       dbcsr_nblkcols_total, dbcsr_nblkrows_local, dbcsr_nblkrows_total, &
       dbcsr_new_transposed, dbcsr_nfullcols_total, dbcsr_nfullrows_total, &
       dbcsr_no_transpose, dbcsr_norm, dbcsr_norm_column, &
       dbcsr_norm_frobenius, dbcsr_norm_maxabsnorm, dbcsr_obj, &
       dbcsr_obj_type_p, dbcsr_print, dbcsr_print_block_sum, dbcsr_put_block, &
       dbcsr_release, dbcsr_replicate, dbcsr_replicate_all, &
       dbcsr_reserve_all_blocks, dbcsr_reserve_block2d, dbcsr_reserve_blocks, &
       dbcsr_reserve_diag_blocks, dbcsr_run_tests, dbcsr_scalar, &
       dbcsr_scalar_fill_all, dbcsr_scalar_get_type, dbcsr_scalar_get_value, &
       dbcsr_scalar_set_type, dbcsr_scalar_type, dbcsr_scalar_zero, &
       dbcsr_scale, dbcsr_scale_by_vector, dbcsr_set, &
       dbcsr_set_conf_avg_elements_images, dbcsr_set_conf_comm_thread_load, &
       dbcsr_set_conf_external_timing, dbcsr_set_conf_max_ele_block, &
       dbcsr_set_conf_mm_driver, dbcsr_set_conf_mm_stacksize, &
       dbcsr_set_conf_nstacks, dbcsr_set_conf_subcomm, &
       dbcsr_set_conf_use_comm_thread, dbcsr_set_conf_use_mpi_rma, &
       dbcsr_set_diag, dbcsr_setname, dbcsr_sum_replicated, &
       dbcsr_test_arnoldi, dbcsr_test_binary_io, dbcsr_test_mm, dbcsr_trace, &
       dbcsr_triu, dbcsr_type_antisymmetric, dbcsr_type_complex_4, &
       dbcsr_type_complex_8, dbcsr_type_complex_default, &
       dbcsr_type_no_symmetry, dbcsr_type_real_4, dbcsr_type_real_8, &
       dbcsr_type_real_default, dbcsr_type_symmetric, dbcsr_valid_index, &
       dbcsr_verify_matrix, dbcsr_wm_use_mutable, dbcsr_work_create, &
       deallocate_arnoldi_data, get_selected_ritz_val, &
       get_selected_ritz_vector, has_acc, has_mpi, heap_fill, heap_get_first, &
       heap_new, heap_release, heap_reset_first, heap_t, mm_driver_acc, &
       mm_driver_blas, mm_driver_matmul, mm_driver_smm, mm_name_acc, &
       mm_name_blas, mm_name_matmul, mm_name_smm, multrec_limit, &
       set_arnoldi_initial_vector, setup_arnoldi_data, swap
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE string_utilities,                ONLY: uppercase

  !$ USE OMP_LIB
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_interface'

  !
  ! Interface to libdbcsr (contains cp2k timing and error) 
  PUBLIC :: dbcsr_error_set, dbcsr_error_stop, dbcsr_get_conf_comm_thread_load,&
            dbcsr_get_conf_avg_elements_images,&
            dbcsr_get_conf_mm_stacksize, dbcsr_get_conf_nstacks,&
            dbcsr_set_conf_external_timing, dbcsr_set_conf_mm_stacksize,&
            dbcsr_get_conf_use_mpi_rma
  PUBLIC :: dbcsr_set_conf_comm_thread_load, dbcsr_set_conf_max_ele_block, dbcsr_set_conf_mm_driver,&
            dbcsr_set_conf_avg_elements_images, dbcsr_set_conf_nstacks, dbcsr_set_conf_subcomm, dbcsr_set_conf_use_comm_thread,&
            has_mpi, dbcsr_set_conf_use_mpi_rma
  PUBLIC :: has_acc, dbcsr_create_dist_block_cyclic, dbcsr_create_dist_r_unrot, dbcsr_distribution_local_cols,&
            dbcsr_distribution_nlocal_rows, dbcsr_distribution_nrows, dbcsr_distribution_local_rows,&
            dbcsr_distribution_ncols, dbcsr_distribution_nlocal_cols,&
            dbcsr_mp_hold, cp_distribution_get_num_images
  PUBLIC :: dbcsr_type_no_symmetry, dbcsr_type_symmetric, dbcsr_type_antisymmetric
  PUBLIC :: dbcsr_type_complex_8, dbcsr_type_real_4, dbcsr_type_real_8, dbcsr_type_complex_4
  PUBLIC :: dbcsr_type_complex_default, dbcsr_type_real_default
  PUBLIC :: convert_offsets_to_sizes, dbcsr_distribution_obj, dbcsr_distribution_row_dist
  PUBLIC :: dbcsr_get_conf_use_comm_thread, dbcsr_data_clear_pointer, dbcsr_data_new, dbcsr_data_init
  PUBLIC :: dbcsr_mp_group, dbcsr_mp_mynode, dbcsr_mp_numnodes, dbcsr_mp_nprows, dbcsr_mp_npcols
  PUBLIC :: dbcsr_iterator_stop, dbcsr_data_obj, dbcsr_iterator, dbcsr_data_release
  PUBLIC :: dbcsr_iterator_blocks_left, dbcsr_iterator_next_block
  PUBLIC :: dbcsr_error_type, dbcsr_mp_obj, dbcsr_norm_frobenius
  PUBLIC :: dbcsr_run_tests, dbcsr_test_mm, dbcsr_test_arnoldi, dbcsr_norm_maxabsnorm,&
            dbcsr_norm_column, dbcsr_distribution_has_threads, dbcsr_test_binary_io
  PUBLIC :: dbcsr_distribution_col_dist, dbcsr_distribution_init, dbcsr_distribution_release
  PUBLIC :: dbcsr_func_artanh, dbcsr_func_dtanh, dbcsr_func_inverse, dbcsr_func_tanh
  PUBLIC :: dbcsr_distribution_mp
  PUBLIC :: convert_sizes_to_offsets
  PUBLIC :: cp_dbcsr_iterator
  PUBLIC :: cp_dbcsr_type
  PUBLIC :: cp_dbcsr_p_type
  PUBLIC :: cp_dbcsr_init_lib, cp_dbcsr_finalize_lib
  PUBLIC :: cp_dbcsr_set
  PUBLIC :: cp_dbcsr_add
  PUBLIC :: cp_dbcsr_scale
  PUBLIC :: cp_dbcsr_scale_by_vector
  PUBLIC :: cp_dbcsr_hadamard_product
  PUBLIC :: cp_dbcsr_transposed
  PUBLIC :: cp_dbcsr_multiply
  PUBLIC :: cp_dbcsr_copy
  PUBLIC :: cp_dbcsr_copy_into_existing
  PUBLIC :: cp_dbcsr_desymmetrize
  PUBLIC :: cp_dbcsr_add_on_diag
  PUBLIC :: cp_dbcsr_get_block_diag
  PUBLIC :: cp_dbcsr_set_diag
  PUBLIC :: cp_dbcsr_get_diag
  PUBLIC :: cp_dbcsr_filter
  PUBLIC :: cp_dbcsr_finalize
  PUBLIC :: cp_dbcsr_distribution_new
  PUBLIC :: cp_dbcsr_create
  PUBLIC :: cp_dbcsr_work_create
  PUBLIC :: cp_dbcsr_verify_matrix
  PUBLIC :: cp_dbcsr_sum_replicated
  PUBLIC :: cp_dbcsr_checksum
  PUBLIC :: cp_dbcsr_trace
  PUBLIC :: cp_dbcsr_print
  PUBLIC :: cp_dbcsr_print_block_sum
  PUBLIC :: cp_dbcsr_init
  PUBLIC :: cp_dbcsr_init_p
  PUBLIC :: cp_dbcsr_release
  PUBLIC :: cp_dbcsr_release_p
  PUBLIC :: cp_dbcsr_distribute
  PUBLIC :: cp_dbcsr_replicate_all
  PUBLIC :: cp_dbcsr_norm
  PUBLIC :: cp_dbcsr_get_info
  PUBLIC :: cp_dbcsr_get_block_p
  PUBLIC :: cp_dbcsr_put_block
  PUBLIC :: cp_dbcsr_iterator_start
  PUBLIC :: cp_dbcsr_iterator_stop
  PUBLIC :: cp_dbcsr_iterator_next_block
  PUBLIC :: cp_dbcsr_mp_new
  PUBLIC :: cp_dbcsr_mp_release
  PUBLIC :: cp_dbcsr_iterator_blocks_left
  PUBLIC :: cp_dbcsr_distribution_release
  PUBLIC :: cp_dbcsr_col_block_sizes
  PUBLIC :: cp_dbcsr_row_block_sizes
  PUBLIC :: cp_create_bl_distribution
  PUBLIC :: cp_dbcsr_get_matrix_type
  PUBLIC :: cp_dbcsr_get_occupation
  PUBLIC :: cp_dbcsr_distribution
  PUBLIC :: cp_dbcsr_nblkrows_total
  PUBLIC :: cp_dbcsr_nblkcols_total
  PUBLIC :: cp_dbcsr_get_num_blocks
  PUBLIC :: cp_dbcsr_get_data_size
  PUBLIC :: cp_dbcsr_has_symmetry
  PUBLIC :: cp_dbcsr_nfullrows_total
  PUBLIC :: cp_dbcsr_nfullcols_total
  !         1234567890123456789012345678901
  PUBLIC :: cp_dbcsr_get_stored_coordinates
  PUBLIC :: cp_dbcsr_valid_index
  PUBLIC :: cp_dbcsr_get_data_type
  PUBLIC :: cp_dbcsr_reserve_block2d
  PUBLIC :: cp_dbcsr_reserve_blocks
  PUBLIC :: cp_dbcsr_reserve_all_blocks
  PUBLIC :: cp_dbcsr_complete_redistribute
  PUBLIC :: cp_dbcsr_gershgorin_norm
  PUBLIC :: cp_dbcsr_frobenius_norm
  PUBLIC :: cp_dbcsr_name
  PUBLIC :: cp_dbcsr_binary_write
  PUBLIC :: cp_dbcsr_binary_read
  PUBLIC :: cp_dbcsr_function_of_elements
  PUBLIC :: cp_dbcsr_triu
  PUBLIC :: cp_dbcsr_clear_mempools
  PUBLIC :: accdrv_avoid_after_busy, accdrv_binning_binsize, accdrv_binning_nbins,&
            accdrv_do_inhomogenous, accdrv_min_flop_process, accdrv_min_flop_sort,&
            accdrv_posterior_buffers, accdrv_posterior_streams,&
            accdrv_priority_buffers, accdrv_priority_streams,&
            dbcsr_get_conf_max_ele_block,&
            dbcsr_get_conf_mm_driver, dbcsr_get_conf_mpi_mem,&
            dbcsr_get_conf_subcomm, mm_driver_acc,&
            mm_driver_blas, mm_driver_matmul, mm_driver_smm,&
            mm_name_acc, mm_name_blas, mm_name_matmul,&
            mm_name_smm, multrec_limit
  PUBLIC :: heap_fill,&
            heap_get_first,&
            heap_new,&
            heap_release,&
            heap_reset_first,&
            heap_t
  PUBLIC :: cp_dbcsr_allocate_matrix_set, cp_dbcsr_deallocate_matrix_set, cp_dbcsr_deallocate_matrix
  PUBLIC :: cp_dbcsr_add_block_node
  PUBLIC :: cp_dbcsr_maxabs
  PUBLIC :: cp_dbcsr_arnoldi_ev, dbcsr_arnoldi_data, cp_dbcsr_setup_arnoldi_data, deallocate_arnoldi_data,&
            get_selected_ritz_val, get_selected_ritz_vec, cp_set_arnoldi_initial_vector
  PUBLIC :: csr_type, cp_convert_csr_to_dbcsr, cp_convert_dbcsr_to_csr,&
            cp_csr_create_from_dbcsr, cp_csr_destroy, cp_csr_create,&
            csr_dbcsr_blkrow_dist, csr_eqrow_floor_dist,&
            cp_csr_print_sparsity
  PUBLIC :: cp_add_work_coordinate, cp_dbcsr_get_data_p,&
            cp_dbcsr_init_random, cp_dbcsr_setname, cp_set_work_size
  TYPE cp_dbcsr_p_type
     TYPE(cp_dbcsr_type), POINTER :: matrix
  END TYPE cp_dbcsr_p_type

  ! the components of this type must remain private to encapsulate better the internals 
  ! of the dbcsr library.
  TYPE cp_dbcsr_type
     INTEGER, PRIVATE                     :: ref_count
     TYPE(dbcsr_obj), PRIVATE             :: matrix
  END TYPE cp_dbcsr_type

  INTERFACE cp_dbcsr_distribution_new
     MODULE PROCEDURE cp_dbcsr_distribution_new_blocks
     MODULE PROCEDURE cp_dbcsr_distribution_new_clusters1
  END INTERFACE

  INTERFACE cp_dbcsr_create
     MODULE PROCEDURE cp_dbcsr_create_new
     MODULE PROCEDURE cp_dbcsr_create_template
  END INTERFACE

  INTERFACE cp_dbcsr_conform_scalar
     MODULE PROCEDURE make_conformant_scalar_d,&
                      make_conformant_scalar_s,&
                      make_conformant_scalar_c,&
                      make_conformant_scalar_z
  END INTERFACE

  INTERFACE cp_dbcsr_trace
     MODULE PROCEDURE cp_dbcsr_trace_ab_d,&
                      cp_dbcsr_trace_ab_s
     MODULE PROCEDURE cp_dbcsr_trace_a_d,&
                      cp_dbcsr_trace_a_s
     MODULE PROCEDURE cp_dbcsr_trace_ab_z,&
                      cp_dbcsr_trace_ab_c
     MODULE PROCEDURE cp_dbcsr_trace_a_z,&
                      cp_dbcsr_trace_a_c
  END INTERFACE

  INTERFACE cp_dbcsr_set
     MODULE PROCEDURE cp_dbcsr_set_d,&
                      cp_dbcsr_set_s,&
                      cp_dbcsr_set_c,&
                      cp_dbcsr_set_z
  END INTERFACE

  INTERFACE cp_dbcsr_add
     MODULE PROCEDURE cp_dbcsr_add_d,&
                      cp_dbcsr_add_s,&
                      cp_dbcsr_add_c,&
                      cp_dbcsr_add_z
  END INTERFACE

  INTERFACE cp_dbcsr_scale
     MODULE PROCEDURE cp_dbcsr_scale_d,&
                      cp_dbcsr_scale_s,&
                      cp_dbcsr_scale_c,&
                      cp_dbcsr_scale_z
  END INTERFACE

  INTERFACE cp_dbcsr_scale_by_vector
     MODULE PROCEDURE cp_dbcsr_scale_by_vector_d,&
                      cp_dbcsr_scale_by_vector_s,&
                      cp_dbcsr_scale_by_vector_c,&
                      cp_dbcsr_scale_by_vector_z
  END INTERFACE

  INTERFACE cp_dbcsr_multiply
     MODULE PROCEDURE cp_dbcsr_multiply_d,&
                      cp_dbcsr_multiply_s,&
                      cp_dbcsr_multiply_c,&
                      cp_dbcsr_multiply_z
  END INTERFACE

  INTERFACE cp_dbcsr_get_block_p
     MODULE PROCEDURE cp_dbcsr_get_block_p_d,&
                      cp_dbcsr_get_block_p_s,&
                      cp_dbcsr_get_block_p_z,&
                      cp_dbcsr_get_block_p_c
     MODULE PROCEDURE cp_dbcsr_get_2d_block_p_d,&
                      cp_dbcsr_get_2d_block_p_s,&
                      cp_dbcsr_get_2d_block_p_z,&
                      cp_dbcsr_get_2d_block_p_c
  END INTERFACE

  INTERFACE cp_dbcsr_put_block
     MODULE PROCEDURE cp_dbcsr_put_block_area
     MODULE PROCEDURE cp_dbcsr_put_block_d,&
                      cp_dbcsr_put_block_s,&
                      cp_dbcsr_put_block_z,&
                      cp_dbcsr_put_block_c
     MODULE PROCEDURE cp_dbcsr_put_block2d_d,&
                      cp_dbcsr_put_block2d_s,&
                      cp_dbcsr_put_block2d_z,&
                      cp_dbcsr_put_block2d_c
  END INTERFACE

  INTERFACE cp_dbcsr_iterator_next_block
     MODULE PROCEDURE cp_iterator_next_2d_block_d,&
                      cp_iterator_next_2d_block_s,&
                      cp_iterator_next_2d_block_c,&
                      cp_iterator_next_2d_block_z
     MODULE PROCEDURE cp_iterator_next_1d_block_d,&
                      cp_iterator_next_1d_block_s,&
                      cp_iterator_next_1d_block_c,&
                      cp_iterator_next_1d_block_z
  END INTERFACE

  INTERFACE cp_dbcsr_reserve_block2d
     MODULE PROCEDURE cp_dbcsr_reserve_block2d_d,&
                      cp_dbcsr_reserve_block2d_s,&
                      cp_dbcsr_reserve_block2d_c,&
                      cp_dbcsr_reserve_block2d_z
  END INTERFACE

  INTERFACE cp_dbcsr_allocate_matrix_set
     MODULE PROCEDURE allocate_dbcsr_matrix_set, allocate_dbcsr_matrix_set_2d
  END INTERFACE

  INTERFACE cp_dbcsr_deallocate_matrix_set
     MODULE PROCEDURE deallocate_dbcsr_matrix_set,&
                      deallocate_dbcsr_matrix_set_2d
  END INTERFACE

  INTERFACE cp_csr_create
     MODULE PROCEDURE cp_csr_create_new,&
                      cp_csr_create_template
  END INTERFACE

  INTERFACE cp_dbcsr_get_data_p
     MODULE PROCEDURE cp_dbcsr_get_data_c_s, cp_dbcsr_get_data_c_c,&
                      cp_dbcsr_get_data_c_d, cp_dbcsr_get_data_c_z
  END INTERFACE


  PRIVATE

CONTAINS

! *****************************************************************************
!> \brief Initializes DBCSR
!> \param group ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_init_lib (group, error)
    INTEGER, INTENT(IN)                      :: group
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_init_lib', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL cite_reference(Borstnik2014)

    CALL dbcsr_init_lib (group, dbcsr_error)
  END SUBROUTINE cp_dbcsr_init_lib

! *****************************************************************************
!> \brief Finalizes DBCSR
!> \param group ...
!> \param output_unit ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_finalize_lib (group, output_unit, error)
    INTEGER, INTENT(IN)                      :: group, output_unit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_finalize_lib', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_finalize_lib (group, output_unit, dbcsr_error)
  END SUBROUTINE cp_dbcsr_finalize_lib

! *****************************************************************************
!> \brief Deallocate memory contained in mempools
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_clear_mempools(error)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_clear_mempools(dbcsr_error)
  END SUBROUTINE cp_dbcsr_clear_mempools

! *****************************************************************************
!> \brief cp2k interface to libdbcsr.
!>     The cp2k error is passed as an argument.
!>     Those procedure can de timed with the cp2k timer.
!> \param matrix ...
!> \retval name ...
!> \author vw
! *****************************************************************************

  FUNCTION cp_dbcsr_name (matrix) RESULT (name)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    CHARACTER(len=default_string_length)     :: name

    name = dbcsr_name (matrix%matrix)
  END FUNCTION cp_dbcsr_name

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param newname ...
! *****************************************************************************

  SUBROUTINE cp_dbcsr_setname (matrix,newname)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    CHARACTER(len=*), INTENT(IN)             :: newname

    CALL dbcsr_setname(matrix%matrix,newname)
  END SUBROUTINE cp_dbcsr_setname

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval norm ...
! *****************************************************************************
  FUNCTION cp_dbcsr_gershgorin_norm(matrix) RESULT (norm)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    REAL(KIND=real_8)                        :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_gershgorin_norm', &
      routineP = moduleN//':'//routineN

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    norm = dbcsr_gershgorin_norm(matrix%matrix)

  END FUNCTION cp_dbcsr_gershgorin_norm

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param local ...
!> \retval norm ...
! *****************************************************************************
  FUNCTION cp_dbcsr_frobenius_norm(matrix, local) RESULT (norm)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    LOGICAL, INTENT(in), OPTIONAL            :: local
    REAL(KIND=real_8)                        :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_frobenius_norm', &
      routineP = moduleN//':'//routineN

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    norm = dbcsr_frobenius_norm(matrix%matrix, local)

  END FUNCTION cp_dbcsr_frobenius_norm

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval norm ...
! *****************************************************************************
  FUNCTION cp_dbcsr_maxabs(matrix) RESULT (norm)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    REAL(KIND=real_8)                        :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_maxabs', &
      routineP = moduleN//':'//routineN

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    norm = dbcsr_maxabs(matrix%matrix)

  END FUNCTION cp_dbcsr_maxabs

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param redist ...
!> \param keep_sparsity ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_complete_redistribute(matrix, redist, keep_sparsity, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: redist
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_complete_redistribute', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_complete_redistribute(matrix%matrix, redist%matrix, keep_sparsity, dbcsr_error)

  END SUBROUTINE cp_dbcsr_complete_redistribute

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param rows ...
!> \param cols ...
!> \param blk_pointers ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_reserve_blocks(matrix, rows, cols, blk_pointers, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, DIMENSION(:), INTENT(IN)        :: rows, cols
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: blk_pointers
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_reserve_blocks', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_reserve_blocks(matrix%matrix, rows, cols, blk_pointers,&
         dbcsr_error)
  END SUBROUTINE cp_dbcsr_reserve_blocks

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_reserve_all_blocks(matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_reserve_all_blocks', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_reserve_all_blocks(matrix%matrix, dbcsr_error)
  END SUBROUTINE cp_dbcsr_reserve_all_blocks

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_reserve_diag_blocks(matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_reserve_diag_blocks', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_reserve_diag_blocks(matrix%matrix, dbcsr_error)
  END SUBROUTINE cp_dbcsr_reserve_diag_blocks

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param index_matrix ...
!> \param row ...
!> \param col ...
!> \param blk ...
!> \param index ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_add_work_coordinate(matrix, index_matrix, row, col, blk, index, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: index_matrix, row, col
    INTEGER, INTENT(IN), OPTIONAL            :: blk
    INTEGER, INTENT(OUT), OPTIONAL           :: index
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_add_work_coordinate', &
      routineP = moduleN//':'//routineN

    CALL add_work_coordinate(matrix%matrix%m%wms(index_matrix), row, col, blk, index, error)
  END SUBROUTINE cp_add_work_coordinate

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param index_matrix ...
!> \param newvalue ...
! *****************************************************************************
  SUBROUTINE cp_set_work_size(matrix, index_matrix, newvalue)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: index_matrix, newvalue

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_set_work_size', &
      routineP = moduleN//':'//routineN

    matrix%matrix%m%wms(index_matrix)%datasize = newvalue
  END SUBROUTINE cp_set_work_size

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param keep_sparsity ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_init_random (matrix, keep_sparsity, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    LOGICAL, OPTIONAL                        :: keep_sparsity
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_init_random', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_init_random(matrix%matrix,keep_sparsity=keep_sparsity,error=error)
  END SUBROUTINE cp_dbcsr_init_random

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval data_type ...
! *****************************************************************************
  PURE FUNCTION cp_dbcsr_get_data_type (matrix) RESULT (data_type)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: data_type

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_data_type', &
      routineP = moduleN//':'//routineN

    data_type = dbcsr_get_data_type (matrix%matrix)
  END FUNCTION cp_dbcsr_get_data_type

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval valid_index ...
! *****************************************************************************
  PURE FUNCTION cp_dbcsr_valid_index (matrix) RESULT (valid_index)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    LOGICAL                                  :: valid_index

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_valid_index', &
      routineP = moduleN//':'//routineN

    valid_index = dbcsr_valid_index(matrix%matrix)

  END FUNCTION cp_dbcsr_valid_index

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param row ...
!> \param column ...
!> \param transpose ...
!> \param processor ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_get_stored_coordinates(matrix, row, column, transpose, processor)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, INTENT(INOUT)                   :: row, column
    LOGICAL, INTENT(INOUT)                   :: transpose
    INTEGER, INTENT(OUT), OPTIONAL           :: processor

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_get_stored_coordinates', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_get_stored_coordinates(matrix%matrix, row, column, transpose, processor)

  END SUBROUTINE cp_dbcsr_get_stored_coordinates

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval num_blocks ...
! *****************************************************************************
  PURE FUNCTION cp_dbcsr_get_num_blocks (matrix) RESULT (num_blocks)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: num_blocks

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_num_blocks', &
      routineP = moduleN//':'//routineN

    num_blocks = dbcsr_get_num_blocks (matrix%matrix)

  END FUNCTION cp_dbcsr_get_num_blocks

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval data_size ...
! *****************************************************************************
  FUNCTION cp_dbcsr_get_data_size (matrix) RESULT (data_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: data_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_data_size', &
      routineP = moduleN//':'//routineN

    data_size = dbcsr_get_data_size(matrix%matrix)

  END FUNCTION cp_dbcsr_get_data_size

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval matrix_type ...
! *****************************************************************************
  PURE FUNCTION cp_dbcsr_get_matrix_type (matrix) RESULT (matrix_type)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    CHARACTER                                :: matrix_type

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_matrix_type', &
      routineP = moduleN//':'//routineN

    matrix_type = dbcsr_get_matrix_type (matrix%matrix)

  END FUNCTION cp_dbcsr_get_matrix_type

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval occupation ...
! *****************************************************************************
  FUNCTION cp_dbcsr_get_occupation (matrix) RESULT (occupation)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    REAL(KIND=real_8)                        :: occupation

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_occupation', &
      routineP = moduleN//':'//routineN

    occupation = dbcsr_get_occupation (matrix%matrix)

  END FUNCTION cp_dbcsr_get_occupation


! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nblkrows_total ...
! *****************************************************************************
  FUNCTION cp_dbcsr_nblkrows_total(matrix) RESULT (nblkrows_total)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nblkrows_total

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nblkrows_total', &
      routineP = moduleN//':'//routineN

    nblkrows_total = dbcsr_nblkrows_total(matrix%matrix)

  END FUNCTION cp_dbcsr_nblkrows_total

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nblkcols_total ...
! *****************************************************************************
  FUNCTION cp_dbcsr_nblkcols_total(matrix) RESULT (nblkcols_total)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nblkcols_total

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nblkcols_total', &
      routineP = moduleN//':'//routineN

    nblkcols_total = dbcsr_nblkcols_total(matrix%matrix)
  END FUNCTION cp_dbcsr_nblkcols_total

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nfullrows_total ...
! *****************************************************************************
  FUNCTION cp_dbcsr_nfullrows_total(matrix) RESULT (nfullrows_total)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nfullrows_total

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nfullrows_total', &
      routineP = moduleN//':'//routineN

    nfullrows_total = dbcsr_nfullrows_total(matrix%matrix)

  END FUNCTION cp_dbcsr_nfullrows_total

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nfullcols_total ...
! *****************************************************************************
  FUNCTION cp_dbcsr_nfullcols_total(matrix) RESULT (nfullcols_total)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nfullcols_total

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nfullcols_total', &
      routineP = moduleN//':'//routineN

    nfullcols_total = dbcsr_nfullcols_total(matrix%matrix)
  END FUNCTION cp_dbcsr_nfullcols_total

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nblkrows_local ...
! *****************************************************************************
  FUNCTION cp_dbcsr_nblkrows_local(matrix) RESULT (nblkrows_local)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nblkrows_local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nblkrows_local', &
      routineP = moduleN//':'//routineN

    nblkrows_local = dbcsr_nblkrows_local(matrix%matrix)
  END FUNCTION cp_dbcsr_nblkrows_local

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nblkcols_local ...
! *****************************************************************************
  FUNCTION cp_dbcsr_nblkcols_local(matrix) RESULT (nblkcols_local)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nblkcols_local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nblkcols_local', &
      routineP = moduleN//':'//routineN

    nblkcols_local = dbcsr_nblkcols_local(matrix%matrix)
  END FUNCTION cp_dbcsr_nblkcols_local

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval distribution ...
! *****************************************************************************
  FUNCTION cp_dbcsr_distribution (matrix) RESULT (distribution)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(dbcsr_distribution_obj)             :: distribution

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_distribution', &
      routineP = moduleN//':'//routineN

    distribution = dbcsr_distribution (matrix%matrix)

  END FUNCTION cp_dbcsr_distribution

! *****************************************************************************
!> \brief ...
!> \param block_distribution ...
!> \param block_size ...
!> \param nelements ...
!> \param nbins ...
! *****************************************************************************
  SUBROUTINE cp_create_bl_distribution (block_distribution,&
       block_size, nelements, nbins)
    INTEGER, DIMENSION(:), INTENT(OUT), &
      POINTER                                :: block_distribution, block_size
    INTEGER, INTENT(IN)                      :: nelements, nbins

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_create_bl_distribution', &
      routineP = moduleN//':'//routineN

    CALL create_bl_distribution (block_distribution,&
       block_size, nelements, nbins)

  END SUBROUTINE cp_create_bl_distribution

! *****************************************************************************
!> \brief ...
!> \param dist ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_distribution_release(dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_distribution_release', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_distribution_release(dist)

  END SUBROUTINE cp_dbcsr_distribution_release

! *****************************************************************************
!> \brief ...
!> \param iterator ...
!> \retval blocks_left ...
! *****************************************************************************
  PURE FUNCTION cp_dbcsr_iterator_blocks_left (iterator) RESULT (blocks_left)
    TYPE(cp_dbcsr_iterator), INTENT(IN)      :: iterator
    LOGICAL                                  :: blocks_left

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_iterator_blocks_left', &
      routineP = moduleN//':'//routineN

    blocks_left = dbcsr_iterator_blocks_left (iterator)

  END FUNCTION cp_dbcsr_iterator_blocks_left

! *****************************************************************************
!> \brief ...
!> \param mp_env ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_mp_release(mp_env)
    TYPE(dbcsr_mp_obj), INTENT(INOUT)        :: mp_env

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_mp_release', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_mp_release(mp_env)

  END SUBROUTINE cp_dbcsr_mp_release

! *****************************************************************************
!> \brief ...
!> \param mp_env ...
!> \param pgrid ...
!> \param mp_group ...
!> \param mynode ...
!> \param numnodes ...
!> \param myprow ...
!> \param mypcol ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_mp_new(mp_env, pgrid, mp_group, mynode, numnodes, myprow,&
       mypcol)
    TYPE(dbcsr_mp_obj), INTENT(OUT)          :: mp_env
    INTEGER, DIMENSION(0:, 0:), INTENT(IN)   :: pgrid
    INTEGER, INTENT(IN)                      :: mp_group, mynode
    INTEGER, INTENT(IN), OPTIONAL            :: numnodes, myprow, mypcol

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_mp_new', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_mp_new(mp_env, pgrid, mp_group, mynode, numnodes, myprow,&
       mypcol)

  END SUBROUTINE cp_dbcsr_mp_new

! *****************************************************************************
!> \brief ...
!> \param iterator ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_iterator_stop (iterator)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)   :: iterator

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_iterator_stop', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_iterator_stop (iterator)

  END SUBROUTINE cp_dbcsr_iterator_stop

! *****************************************************************************
!> \brief ...
!> \param iterator ...
!> \param matrix ...
!> \param shared ...
!> \param dynamic ...
!> \param dynamic_byrows ...
!> \param contiguous_pointers ...
!> \param read_only ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_iterator_start (iterator, matrix, shared, dynamic,&
       dynamic_byrows, contiguous_pointers, read_only)
    TYPE(cp_dbcsr_iterator), INTENT(OUT)     :: iterator
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: shared, dynamic, &
                                                dynamic_byrows, &
                                                contiguous_pointers, read_only

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_iterator_start', &
      routineP = moduleN//':'//routineN

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL dbcsr_iterator_start (iterator, matrix%matrix, shared, dynamic,&
       dynamic_byrows, contiguous_pointers, read_only)

  END SUBROUTINE cp_dbcsr_iterator_start


! *****************************************************************************
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
!> \param data_area ...
!> \param matrix_type ...
!> \param data_type ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_get_info(matrix, nblkrows_total, nblkcols_total,&
       nfullrows_total, nfullcols_total,&
       nblkrows_local, nblkcols_local,&
       nfullrows_local, nfullcols_local,&
       my_prow, my_pcol,&
       local_rows, local_cols, proc_row_dist, proc_col_dist,&
       row_blk_size, col_blk_size, row_blk_offset, col_blk_offset, distribution, name, data_area,&
       matrix_type, data_type)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, INTENT(OUT), OPTIONAL :: nblkrows_total, nblkcols_total, &
      nfullrows_total, nfullcols_total, nblkrows_local, nblkcols_local, &
      nfullrows_local, nfullcols_local, my_prow, my_pcol
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: local_rows, local_cols, &
                                                proc_row_dist, proc_col_dist
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL, POINTER                      :: row_blk_size, col_blk_size, &
                                                row_blk_offset, col_blk_offset
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT), OPTIONAL                  :: distribution
    CHARACTER(len=*), INTENT(OUT), OPTIONAL  :: name
    TYPE(dbcsr_data_obj), INTENT(OUT), &
      OPTIONAL                               :: data_area
    CHARACTER, OPTIONAL                      :: matrix_type
    INTEGER, OPTIONAL                        :: data_type

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_info', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_get_info(matrix%matrix, nblkrows_total, nblkcols_total,&
         nfullrows_total, nfullcols_total,&
         nblkrows_local, nblkcols_local,&
         nfullrows_local, nfullcols_local,&
         my_prow, my_pcol,&
         local_rows, local_cols, proc_row_dist, proc_col_dist,&
         row_blk_size, col_blk_size, row_blk_offset, col_blk_offset, distribution, name, data_area,&
         matrix_type, data_type)

  END SUBROUTINE cp_dbcsr_get_info

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param which_norm ...
!> \param norm_scalar ...
!> \param norm_vector ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_norm(matrix, which_norm, norm_scalar, norm_vector, error)

    TYPE(cp_dbcsr_type), INTENT(INOUT), &
      TARGET                                 :: matrix
    INTEGER, INTENT(IN)                      :: which_norm
    REAL(dp), INTENT(OUT), OPTIONAL          :: norm_scalar
    REAL(dp), DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: norm_vector
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_norm', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    IF (PRESENT (norm_scalar)) THEN
       CALL dbcsr_norm(matrix%matrix, which_norm, norm_scalar=norm_scalar,&
            error=dbcsr_error)
    ELSEIF (PRESENT (norm_vector)) THEN
       CALL dbcsr_norm (matrix%matrix, which_norm, norm_vector=norm_vector,&
            error=dbcsr_error)
    ELSE
       CALL cp_assert (.FALSE., cp_wrong_args_error, cp_internal_error,&
            routineN, "Must pass either scalar or vector norm.", error=error)
    ENDIF

  END SUBROUTINE cp_dbcsr_norm

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_replicate_all(matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_replicate_all', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_replicate_all(matrix%matrix, dbcsr_error)
  END SUBROUTINE cp_dbcsr_replicate_all

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param replicate_rows ...
!> \param replicate_columns ...
!> \param restrict_source ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_replicate(matrix, replicate_rows, replicate_columns,&
       restrict_source, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    LOGICAL, INTENT(IN)                      :: replicate_rows, &
                                                replicate_columns
    INTEGER, INTENT(IN), OPTIONAL            :: restrict_source
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_replicate', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_replicate(matrix%matrix, replicate_rows, replicate_columns,&
         restrict_source, dbcsr_error)
  END SUBROUTINE cp_dbcsr_replicate

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param fast ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_distribute(matrix, fast, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    LOGICAL, INTENT(in), OPTIONAL            :: fast
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_distribute', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_distribute(matrix%matrix, fast, dbcsr_error)

  END SUBROUTINE cp_dbcsr_distribute

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_release_p (matrix, error)
    TYPE(cp_dbcsr_type), POINTER             :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_release_p', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(matrix)) THEN
       CALL cp_dbcsr_release (matrix, error)
       DEALLOCATE(matrix)
    ENDIF


  END SUBROUTINE cp_dbcsr_release_p

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_release (matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_release', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_release (matrix%matrix)

  END SUBROUTINE cp_dbcsr_release

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_init (matrix, error)
    TYPE(cp_dbcsr_type), INTENT(OUT)         :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_init', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_init(matrix%matrix)
    matrix%ref_count = 0

  END SUBROUTINE cp_dbcsr_init

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_init_p (matrix, error)
    TYPE(cp_dbcsr_type), POINTER             :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_init_p', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(matrix)) THEN
       CALL cp_dbcsr_release(matrix, error)
       DEALLOCATE(matrix)
    ENDIF

    ALLOCATE(matrix)
    CALL cp_dbcsr_init (matrix, error)

  END SUBROUTINE cp_dbcsr_init_p

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param nodata ...
!> \param matlab_format ...
!> \param variable_name ...
!> \param unit_nr ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_print(matrix, nodata, matlab_format, variable_name, unit_nr, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: nodata, matlab_format
    CHARACTER(*), INTENT(in), OPTIONAL       :: variable_name
    INTEGER, OPTIONAL                        :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_print', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_print(matrix%matrix, nodata, matlab_format, variable_name, unit_nr, dbcsr_error)

  END SUBROUTINE cp_dbcsr_print

! *****************************************************************************
!> \brief Prints the sum of the elements in each block
!> \param matrix ...
!> \param unit_nr ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_print_block_sum(matrix, unit_nr, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, OPTIONAL                        :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_print_block_sum', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_print_block_sum(matrix%matrix, unit_nr, dbcsr_error)

  END SUBROUTINE cp_dbcsr_print_block_sum

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param local ...
!> \param pos ...
!> \param error ...
!> \retval checksum ...
! *****************************************************************************
  FUNCTION cp_dbcsr_checksum(matrix, local, pos, error) RESULT(checksum)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: local, pos
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    REAL(KIND=dp)                            :: checksum

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_checksum', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    checksum = dbcsr_checksum(matrix%matrix,&
         local=local, pos=pos, error=dbcsr_error)
  END FUNCTION cp_dbcsr_checksum

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_sum_replicated (matrix, error)
    TYPE(cp_dbcsr_type), INTENT(inout)       :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_sum_replicated', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_sum_replicated (matrix%matrix, dbcsr_error)
  END SUBROUTINE cp_dbcsr_sum_replicated

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_triu(matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_triu', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

!CALL cluster_buffer_flush (matrix%buffers, error=error)

    CALL dbcsr_triu(matrix%matrix, dbcsr_error)
    !CALL dbcsr_index_compact(matrix%matrix, dbcsr_error)

  END SUBROUTINE cp_dbcsr_triu

! *****************************************************************************
!> \brief ...
!> \param m ...
!> \param error ...
!> \param verbosity ...
!> \param local ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_verify_matrix(m, error, verbosity, local)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: m
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER, INTENT(IN), OPTIONAL            :: verbosity
    LOGICAL, INTENT(IN), OPTIONAL            :: local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_verify_matrix', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_verify_matrix(m%matrix, verbosity, local, dbcsr_error)

  END SUBROUTINE cp_dbcsr_verify_matrix

! *****************************************************************************
!> \brief Creates new distribution from only block distributions
!> \param[out] dist           distribution
!> \param[in] mp_env          multiprocessing environment
!> \param row_dist row distribution
!> \param col_dist column distribution
!> \param reuse_arrays ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_distribution_new_blocks(dist, mp_env, row_dist, col_dist, reuse_arrays)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, DIMENSION(:), INTENT(IN), &
      POINTER                                :: row_dist, col_dist
    LOGICAL, INTENT(IN), OPTIONAL            :: reuse_arrays

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_distribution_new_blocks', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_distribution_new(dist, mp_env, row_dist, col_dist, reuse_arrays=reuse_arrays)

  END SUBROUTINE cp_dbcsr_distribution_new_blocks

! *****************************************************************************
!> \brief Creates new distribution from block and cluster distributions
!> \param[out] dist           distribution
!> \param[in] mp_env          multiprocessing environment
!> \param row_dist row distribution
!> \param col_dist column distribution
!> \param reuse_arrays ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_distribution_new_clusters1(dist, mp_env, row_dist, col_dist, reuse_arrays)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, DIMENSION(:, :), INTENT(IN), &
      POINTER                                :: row_dist, col_dist
    LOGICAL, INTENT(IN), OPTIONAL            :: reuse_arrays

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_distribution_new_clusters1', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(:), POINTER           :: col_dist_block, &
                                                col_dist_cluster, &
                                                row_dist_block, &
                                                row_dist_cluster

    row_dist_block => row_dist(:,1)
    col_dist_block => col_dist(:,1)
    row_dist_cluster => row_dist(:,2)
    col_dist_cluster => col_dist(:,2)
    IF (ANY(row_dist_cluster.eq.-1)) NULLIFY(row_dist_cluster)
    IF (ANY(col_dist_cluster.eq.-1)) NULLIFY(col_dist_cluster)
    CALL dbcsr_distribution_new(dist, mp_env, row_dist_block, col_dist_block, &
                                row_dist_cluster=row_dist_cluster, &
                                col_dist_cluster=col_dist_cluster, &
                                reuse_arrays=reuse_arrays)

  END SUBROUTINE cp_dbcsr_distribution_new_clusters1

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param reshuffle ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_finalize(matrix, reshuffle, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: reshuffle
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_finalize', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_finalize(matrix%matrix, reshuffle, error=dbcsr_error)
  END SUBROUTINE cp_dbcsr_finalize

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param nblks_guess ...
!> \param sizedata_guess ...
!> \param n ...
!> \param error ...
!> \param work_mutable ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_work_create(matrix, nblks_guess, sizedata_guess, n,&
       error, work_mutable)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN), OPTIONAL            :: nblks_guess, sizedata_guess, n
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(in), OPTIONAL            :: work_mutable

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_work_create', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_work_create(matrix%matrix, nblks_guess, sizedata_guess, n,&
         work_mutable, error=dbcsr_error)
  END SUBROUTINE cp_dbcsr_work_create

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param name ...
!> \param dist ...
!> \param matrix_type ...
!> \param row_blk_size ...
!> \param col_blk_size ...
!> \param nblks ...
!> \param nze ...
!> \param data_type ...
!> \param reuse ...
!> \param reuse_arrays ...
!> \param mutable_work ...
!> \param replication_type ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_create_new(matrix, name, dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, data_type, reuse,&
       reuse_arrays, mutable_work, replication_type, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    CHARACTER(len=*), INTENT(IN)             :: name
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    CHARACTER, INTENT(IN)                    :: matrix_type
    INTEGER, DIMENSION(:), INTENT(IN), &
      POINTER                                :: row_blk_size, col_blk_size
    INTEGER, INTENT(IN), OPTIONAL            :: nblks, nze, data_type
    LOGICAL, INTENT(IN), OPTIONAL            :: reuse, reuse_arrays, &
                                                mutable_work
    CHARACTER, INTENT(IN), OPTIONAL          :: replication_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_create_new', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_create(matrix%matrix, name, dist,&
         matrix_type,&
         row_blk_size, col_blk_size, nblks, nze, data_type, reuse=reuse,&
         reuse_arrays=reuse_arrays,&
         mutable_work=mutable_work, replication_type=replication_type,&
         error=dbcsr_error)

    matrix%ref_count = 1

  END SUBROUTINE cp_dbcsr_create_new

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param name ...
!> \param template ...
!> \param dist ...
!> \param matrix_type ...
!> \param row_blk_size ...
!> \param col_blk_size ...
!> \param nblks ...
!> \param nze ...
!> \param data_type ...
!> \param reuse ...
!> \param reuse_arrays ...
!> \param mutable_work ...
!> \param replication_type ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_create_template(matrix, name, template,&
       dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, data_type, reuse,&
       reuse_arrays, mutable_work, replication_type, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: name
    TYPE(cp_dbcsr_type), INTENT(IN)          :: template
    TYPE(dbcsr_distribution_obj), &
      INTENT(IN), OPTIONAL                   :: dist
    CHARACTER, INTENT(IN), OPTIONAL          :: matrix_type
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL, POINTER                      :: row_blk_size, col_blk_size
    INTEGER, INTENT(IN), OPTIONAL            :: nblks, nze, data_type
    LOGICAL, INTENT(IN), OPTIONAL            :: reuse, reuse_arrays, &
                                                mutable_work
    CHARACTER, INTENT(IN), OPTIONAL          :: replication_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_create_template', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_create(matrix%matrix, template%matrix,&
         name=name, dist=dist, matrix_type=matrix_type,&
         row_blk_size=row_blk_size, col_blk_size=col_blk_size,&
         nblks=nblks, nze=nze, data_type=data_type,&
         reuse=reuse, reuse_arrays=reuse_arrays, mutable_work=mutable_work,&
         replication_type=replication_type,&
         error=dbcsr_error)
    !
    matrix%ref_count = 1
  END SUBROUTINE cp_dbcsr_create_template

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eps ...
!> \param method ...
!> \param use_absolute ...
!> \param filter_diag ...
!> \param thorough ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_filter(matrix, eps, method, use_absolute, filter_diag, &
       thorough, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    REAL(dp), INTENT(IN)                     :: eps
    INTEGER, INTENT(IN), OPTIONAL            :: method
    LOGICAL, INTENT(in), OPTIONAL            :: use_absolute, filter_diag, &
                                                thorough
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_filter', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: quick
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    IF (PRESENT (thorough)) THEN
       quick = .NOT. thorough
    ELSE
       quick = .FALSE.
    ENDIF
    CALL dbcsr_filter(matrix%matrix, cp_dbcsr_conform_scalar (eps, matrix, error),&
                      method, use_absolute, filter_diag,&
         quick=quick, error=dbcsr_error)
  END SUBROUTINE cp_dbcsr_filter

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param diag ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_set_diag(matrix, diag, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    REAL(dp), DIMENSION(:), INTENT(IN)       :: diag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_set_diag', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_set_diag(matrix%matrix, diag, dbcsr_error)
  END SUBROUTINE cp_dbcsr_set_diag

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param diag ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_get_diag(matrix, diag, error)

    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: diag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_diag', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_get_diag(matrix%matrix, diag, dbcsr_error)
  END SUBROUTINE cp_dbcsr_get_diag

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param diag ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_get_block_diag(matrix, diag, error)

    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: diag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_block_diag', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL cp_dbcsr_create (diag, "Diagonal of "//TRIM(cp_dbcsr_name (matrix)),&
         template=matrix,&
         error=error)
    CALL dbcsr_get_block_diag(matrix%matrix, diag%matrix, dbcsr_error)
  END SUBROUTINE cp_dbcsr_get_block_diag

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param alpha_scalar ...
!> \param first_row ...
!> \param last_row ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_add_on_diag(matrix, alpha_scalar, first_row, last_row, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    REAL(KIND=dp), INTENT(IN)                :: alpha_scalar
    INTEGER, INTENT(in), OPTIONAL            :: first_row, last_row
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_add_on_diag', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_add_on_diag(matrix%matrix,&
         cp_dbcsr_conform_scalar (alpha_scalar, matrix, error),&
         first_row, last_row, dbcsr_error)
  END SUBROUTINE cp_dbcsr_add_on_diag

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param filepath ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_binary_write(matrix, filepath, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    CHARACTER(LEN=*), INTENT(IN)             :: filepath
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_binary_write', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_binary_write(matrix%matrix, filepath, dbcsr_error)

  END SUBROUTINE cp_dbcsr_binary_write

! *****************************************************************************
!> \brief ...
!> \param filepath ...
!> \param distribution ...
!> \param groupid ...
!> \param matrix_new ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_binary_read(filepath, distribution, groupid, matrix_new, error)
    CHARACTER(len=*), INTENT(IN)             :: filepath
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: distribution
    INTEGER, INTENT(IN), OPTIONAL            :: groupid
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_new
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_binary_read', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_binary_read(filepath, distribution, groupid, matrix_new%matrix, dbcsr_error)

  END SUBROUTINE cp_dbcsr_binary_read

! *****************************************************************************
!> \brief ...
!> \param matrix_b ...
!> \param matrix_a ...
!> \param name ...
!> \param error ...
!> \param keep_sparsity ...
!> \param shallow_data ...
!> \param keep_imaginary ...
!> \param matrix_type ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_copy(matrix_b, matrix_a, name, error, keep_sparsity,&
       shallow_data, keep_imaginary, matrix_type)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_b
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity, shallow_data, &
                                                keep_imaginary
    CHARACTER, INTENT(IN), OPTIONAL          :: matrix_type

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_copy', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

!call cp_assert (matrix_b%ref_count .gt. 0,&
!     cp_warning_level, cp_caller_error, routineN,&
!     "Matrix not created.", error=error)

    IF (matrix_b%ref_count .EQ. 0) THEN
       CALL cp_dbcsr_create (matrix_b, template=matrix_a,&
            error=error)
    ENDIF
    CALL dbcsr_copy(matrix_b%matrix, matrix_a%matrix, name, keep_sparsity,&
         shallow_data, keep_imaginary, matrix_type, dbcsr_error)
  END SUBROUTINE cp_dbcsr_copy

! *****************************************************************************
!> \brief ...
!> \param matrix_b ...
!> \param matrix_a ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_copy_into_existing(matrix_b, matrix_a, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_b
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_copy_into_existing', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    IF (matrix_b%ref_count .EQ. 0) THEN
       CALL cp_dbcsr_create (matrix_b, template=matrix_a,&
            error=error)
    ENDIF
    CALL dbcsr_copy_into_existing(matrix_b%matrix, matrix_a%matrix, dbcsr_error)
  END SUBROUTINE cp_dbcsr_copy_into_existing


! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_desymmetrize(matrix_a, matrix_b, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_b
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_desymmetrize', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

!call cp_assert (matrix_b%ref_count .gt. 0,&
!     cp_warning_level, cp_caller_error, routineN,&
!     "Matrix not created.", error=error)

    IF (matrix_b%ref_count .EQ. 0) THEN
       CALL cp_dbcsr_create (matrix_b, template=matrix_a,&
            error=error)
    ENDIF
    CALL dbcsr_desymmetrize_deep(matrix_a%matrix, matrix_b%matrix,&
         untransposed_data = .TRUE., error=dbcsr_error)
  END SUBROUTINE cp_dbcsr_desymmetrize


! *****************************************************************************
!> \brief ...
!> \param transposed ...
!> \param normal ...
!> \param shallow_data_copy ...
!> \param transpose_data ...
!> \param transpose_distribution ...
!> \param use_distribution ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_transposed (transposed, normal, shallow_data_copy,&
       transpose_data, transpose_distribution, use_distribution, error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: transposed
    TYPE(cp_dbcsr_type), INTENT(IN)          :: normal
    LOGICAL, INTENT(IN), OPTIONAL            :: shallow_data_copy, &
                                                transpose_data, &
                                                transpose_distribution
    TYPE(dbcsr_distribution_obj), &
      INTENT(IN), OPTIONAL                   :: use_distribution
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_transposed', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: myshallow_data_copy, &
                                                mytranspose_distribution
    TYPE(dbcsr_distribution_obj)             :: myuse_distribution
    TYPE(dbcsr_error_type)                   :: dbcsr_error

!   set some defaults to make usage a bit less painful (fschiff)

    myshallow_data_copy=.FALSE.
    myuse_distribution=cp_dbcsr_distribution(normal)
    mytranspose_distribution=.FALSE.
    IF(PRESENT(shallow_data_copy)) myshallow_data_copy=shallow_data_copy
    IF(PRESENT(use_distribution))myuse_distribution=use_distribution
    IF(PRESENT(transpose_distribution))mytranspose_distribution=transpose_distribution

    CALL dbcsr_new_transposed(transposed%matrix, normal%matrix, myshallow_data_copy,&
         transpose_data, mytranspose_distribution,&
         use_distribution=myuse_distribution, error=dbcsr_error)
  END SUBROUTINE cp_dbcsr_transposed

! *****************************************************************************
!> \brief creates a dbcsr col vector like object which lives on proc_col 0
!>        and has the same row dist as the template matrix
!>        the returned matrix is fully allocated and all blocks are set to 0
!>        this is not a sparse object (and must never be)
!> \param vec  the vector object to create must be allocated but not initialized
!> \param matrix a dbcsr matrix used as template
!> \param ncol number of vectors in the dbcsr_object (1 for vector, n for skinny matrix)
!> \param error ...
! *****************************************************************************
   SUBROUTINE cp_dbcsr_create_col_vec_from_matrix(vec,matrix,ncol,error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: vec
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: ncol
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_create_col_vec_from_matrix', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL create_col_vec_from_matrix(vec%matrix,matrix%matrix,ncol,dbcsr_error)

   END SUBROUTINE cp_dbcsr_create_col_vec_from_matrix

! *****************************************************************************
!> \brief creates a dbcsr row vector like object which lives on proc_row 0
!>        and has the same row dist as the template matrix
!>        the returned matrix is fully allocated and all blocks are set to 0
!>        this is not a sparse object (and must never be)
!> \param vec  the vector object to create must be allocated but not initialized
!> \param matrix a dbcsr matrix used as template
!> \param nrow number of vectors in the dbcsr_object (1 for vector, n for skinny matrix)
!> \param error ...
! *****************************************************************************
   SUBROUTINE cp_dbcsr_create_row_vec_from_matrix(vec,matrix,nrow,error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: vec
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nrow
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_create_row_vec_from_matrix', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL create_row_vec_from_matrix(vec%matrix,matrix%matrix,nrow,dbcsr_error)

   END SUBROUTINE cp_dbcsr_create_row_vec_from_matrix

! *****************************************************************************
!> \brief creates a col vector like object whose blocks can be replicated
!>        along the processor row and has the same row dist as the template matrix
!>        the returned matrix is fully allocated and all blocks are set to 0
!>        this is not a sparse object (and must never be)
!> \param vec  the vector object to create must be allocated but not initialized
!> \param matrix a dbcsr matrix used as template
!> \param ncol number of vectors in the dbcsr_object (1 for vector, n for skinny matrix)
!> \param error ...
! *****************************************************************************
   SUBROUTINE cp_dbcsr_create_rep_col_vec_from_matrix(vec,matrix,ncol,error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: vec
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: ncol
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_create_rep_col_vec_from_matrix', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL create_replicated_col_vec_from_matrix(vec%matrix,matrix%matrix,ncol,dbcsr_error)

  END SUBROUTINE  cp_dbcsr_create_rep_col_vec_from_matrix

! *****************************************************************************
!> \brief creates a row vector like object whose blocks can be replicated
!>        along the processor col and has the same col dist as the template matrix
!>        the returned matrix is fully allocated and all blocks are set to 0
!>        this is not a sparse object (and must never be)
!> \param vec  the vector object to create must be allocated but not initialized
!> \param matrix a dbcsr matrix used as template
!> \param nrow number of vectors in the dbcsr_object (1 for vector, n for skinny matrix)
!> \param error ...
! *****************************************************************************
   SUBROUTINE cp_dbcsr_create_rep_row_vec_from_matrix(vec,matrix,nrow,error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: vec
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nrow
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_create_rep_row_vec_from_matrix', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL create_replicated_row_vec_from_matrix(vec%matrix,matrix%matrix,nrow,dbcsr_error)

  END SUBROUTINE  cp_dbcsr_create_rep_row_vec_from_matrix

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param arnoldi_data ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_arnoldi_ev(matrix,arnoldi_data,error)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      INTENT(IN)                             :: matrix
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_arnoldi_ev', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_obj_type_p), ALLOCATABLE, &
      DIMENSION(:)                           :: m_ptrs

! convert cp_dbcsr_p_type to dbcsr_obj_type_p

    ALLOCATE(m_ptrs(SIZE(matrix)))
    DO i=1,SIZE(matrix)
       m_ptrs(i)%matrix=>matrix(i)%matrix%matrix
    END DO

    CALL dbcsr_arnoldi_ev(m_ptrs,arnoldi_data,dbcsr_error)    
    
    DEALLOCATE(m_ptrs)
    
  END SUBROUTINE cp_dbcsr_arnoldi_ev  

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param matrix ...
!> \param max_iter ...
!> \param threshold ...
!> \param selection_crit ...
!> \param nval_request ...
!> \param nrestarts ...
!> \param generalized_ev ...
!> \param iram ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_setup_arnoldi_data(arnoldi_data,matrix,max_iter, threshold, selection_crit, &
                                         nval_request, nrestarts, generalized_ev, iram)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(cp_dbcsr_p_type), DIMENSION(:)      :: matrix
    INTEGER                                  :: max_iter
    REAL(real_8)                             :: threshold
    INTEGER                                  :: selection_crit, nval_request, &
                                                nrestarts
    LOGICAL                                  :: generalized_ev, iram

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_setup_arnoldi_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    TYPE(dbcsr_obj_type_p), ALLOCATABLE, &
      DIMENSION(:)                           :: m_ptrs

! convert cp_dbcsr_p_type to dbcsr_obj_type_p

    ALLOCATE(m_ptrs(SIZE(matrix)))
    DO i=1,SIZE(matrix)
       m_ptrs(i)%matrix=>matrix(i)%matrix%matrix
    END DO

    CALL setup_arnoldi_data(arnoldi_data,m_ptrs ,max_iter, threshold,selection_crit, & 
                            nval_request, nrestarts, generalized_ev, iram)
      

    DEALLOCATE(m_ptrs)
  END SUBROUTINE cp_dbcsr_setup_arnoldi_data

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param vector ...
! *****************************************************************************
  SUBROUTINE cp_set_arnoldi_initial_vector(arnoldi_data,vector)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    TYPE(cp_dbcsr_type)                      :: vector

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_set_arnoldi_initial_vector', &
      routineP = moduleN//':'//routineN

    CALL set_arnoldi_initial_vector(arnoldi_data,vector%matrix)

  END SUBROUTINE cp_set_arnoldi_initial_vector

! *****************************************************************************
!> \brief ...
!> \param arnoldi_data ...
!> \param ind ...
!> \param matrix ...
!> \param vector ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_selected_ritz_vec(arnoldi_data,ind,matrix,vector,error)
    TYPE(dbcsr_arnoldi_data)                 :: arnoldi_data
    INTEGER                                  :: ind
    TYPE(cp_dbcsr_type)                      :: matrix, vector
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_selected_ritz_vec', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL get_selected_ritz_vector(arnoldi_data,ind,matrix%matrix,vector%matrix,dbcsr_error)

  END SUBROUTINE get_selected_ritz_vec

! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param func ...
!> \param a0 ...
!> \param a1 ...
!> \param a2 ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_function_of_elements(matrix_a, func, a0, a1, a2, error) 
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    INTEGER, INTENT(IN)                      :: func
    REAL(kind=dp), INTENT(IN), OPTIONAL      :: a0, a1, a2
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_function_of_elements', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_function_of_elements(matrix_a%matrix,func,a0,a1,a2,dbcsr_error) 
  END SUBROUTINE cp_dbcsr_function_of_elements 
 
! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param matrix_c ...
!> \param b_assume_value ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_hadamard_product(matrix_a, matrix_b, matrix_c, b_assume_value,& 
    error) 
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a, matrix_b
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_c
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: b_assume_value
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_hadamard_product', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_hadamard_product(matrix_a%matrix, matrix_b%matrix, matrix_c%matrix,&
            b_assume_value, dbcsr_error)
  END SUBROUTINE cp_dbcsr_hadamard_product

! *****************************************************************************
!> \brief     Adjusts matrices to be compatibly for multiplication.
!>
!> In BLAS xGEMM or PxGEMM style, CP2K allows with differing
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
!> \param error ...
! *****************************************************************************
  SUBROUTINE matrix_match_sizes (matrix_c, matrix_a, tr_a, matrix_b, tr_b,&
       new_a, new_b, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_c, matrix_a
    CHARACTER, INTENT(IN)                    :: tr_a
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_b
    CHARACTER, INTENT(IN)                    :: tr_b
    TYPE(cp_dbcsr_type), INTENT(OUT)         :: new_a, new_b
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_match_sizes', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: tr_a_l, tr_b_l
    INTEGER                                  :: a_cols_total, a_rows_total, &
                                                b_cols_total, b_rows_total, &
                                                c_cols_total, c_rows_total
    INTEGER, DIMENSION(:), POINTER           :: new_col_dist, new_col_size, &
                                                new_row_dist, new_row_size
    LOGICAL                                  :: atr, btr, create_col_dist, &
                                                create_row_dist, new_a_rows, &
                                                new_b_cols, new_b_rows
    TYPE(dbcsr_distribution_obj)             :: a_dist, b_dist, c_dist, &
                                                new_dist

!   ---------------------------------------------------------------------------

    tr_a_l = tr_a
    CALL uppercase (tr_a_l)
    tr_b_l = tr_b
    CALL uppercase (tr_b_l)

    create_col_dist = .FALSE.
    create_row_dist = .FALSE.

    btr = tr_b_l .NE. dbcsr_no_transpose
    atr = tr_a_l .NE. dbcsr_no_transpose
    !
    c_rows_total = cp_dbcsr_nfullrows_total (matrix_c)
    c_cols_total = cp_dbcsr_nfullcols_total (matrix_c)
    b_rows_total = cp_dbcsr_nfullrows_total (matrix_b)
    b_cols_total = cp_dbcsr_nfullcols_total (matrix_b)
    a_rows_total = cp_dbcsr_nfullrows_total (matrix_a)
    a_cols_total = cp_dbcsr_nfullcols_total (matrix_a)
    IF (atr) CALL swap (a_cols_total, a_rows_total)
    IF (btr) CALL swap (b_cols_total, b_rows_total)
    !
    ! Process matrix B.
    !
    ! This check is faster than explicity comparing blocked row and
    ! column sizes.
    new_b_cols = c_cols_total .NE. b_cols_total
    new_b_rows = a_cols_total .NE. b_rows_total
    a_dist = cp_dbcsr_distribution (matrix_a)
    b_dist = cp_dbcsr_distribution (matrix_b)
    c_dist = cp_dbcsr_distribution (matrix_c)
    IF (new_b_rows .OR. new_b_cols) THEN
       CALL cp_dbcsr_init (new_b, error=error)
       IF (.NOT. btr) THEN
          IF (new_b_cols) THEN
             CALL match_1_dist (new_col_dist,&
                  dbcsr_distribution_col_dist (b_dist),&
                  dbcsr_distribution_col_dist (c_dist))
             create_col_dist = .TRUE.
          ELSE
             new_col_dist => dbcsr_distribution_col_dist (b_dist)
          ENDIF
          new_col_size => cp_dbcsr_col_block_sizes (matrix_c)
          IF (new_b_rows) THEN
             IF (.NOT. atr) THEN
                CALL match_1_dist (new_row_dist,&
                     dbcsr_distribution_row_dist(b_dist),&
                     dbcsr_distribution_col_dist(a_dist))
                create_row_dist = .TRUE.
                new_row_size => cp_dbcsr_col_block_sizes(matrix_a)
             ELSE
                CALL match_1_dist (new_row_dist,&
                     dbcsr_distribution_row_dist(b_dist),&
                     dbcsr_distribution_row_dist(a_dist))
                create_row_dist = .TRUE.
                new_row_size => cp_dbcsr_row_block_sizes(matrix_a)
             ENDIF
          ELSE
             new_row_dist => dbcsr_distribution_row_dist(b_dist)
             new_row_size => cp_dbcsr_row_block_sizes(matrix_b)
          ENDIF
       ELSE
          IF (new_b_cols) THEN
             CALL match_1_dist (new_row_dist,&
                  dbcsr_distribution_row_dist (b_dist),&
                  dbcsr_distribution_col_dist (c_dist))
             create_row_dist = .TRUE.
          ELSE
             new_row_dist => dbcsr_distribution_row_dist (b_dist)
          ENDIF
          new_row_size => cp_dbcsr_col_block_sizes(matrix_c)
          IF (new_b_rows) THEN
             IF (.not.atr) THEN
                CALL match_1_dist (new_col_dist,&
                     dbcsr_distribution_col_dist (b_dist),&
                     dbcsr_distribution_col_dist (a_dist))
                create_col_dist = .TRUE.
                new_col_size => cp_dbcsr_col_block_sizes(matrix_a)
             ELSE
                CALL match_1_dist (new_col_dist,&
                     dbcsr_distribution_col_dist (b_dist),&
                     dbcsr_distribution_row_dist (a_dist))
                create_col_dist = .TRUE.
                new_col_size => cp_dbcsr_row_block_sizes(matrix_a)
             ENDIF
          ELSE
             new_col_dist => dbcsr_distribution_col_dist(b_dist)
             new_col_size => cp_dbcsr_col_block_sizes(matrix_b)
          ENDIF
       ENDIF
       CALL dbcsr_distribution_new (new_dist, dbcsr_distribution_mp (c_dist),&
                                    new_row_dist, new_col_dist)
       IF (create_col_dist) DEALLOCATE(new_col_dist)
       IF (create_row_dist) DEALLOCATE(new_row_dist)
       CALL cp_dbcsr_create (new_b, template=matrix_b,&
            dist = new_dist,&
            row_blk_size = new_row_size,&
            col_blk_size = new_col_size,&
            error=error)
       CALL dbcsr_distribution_release (new_dist)
       CALL cp_dbcsr_complete_redistribute (matrix_b, new_b, error=error)
    ELSE
       CALL cp_dbcsr_init (new_b, error=error)
       CALL cp_dbcsr_copy (new_b, matrix_b, shallow_data=.TRUE., error=error)
    ENDIF

    create_col_dist = .FALSE.
    create_row_dist = .FALSE.
    !
    ! Process matrix A
    new_a_rows = a_rows_total .NE. c_rows_total
    IF (new_a_rows) THEN
       CALL cp_dbcsr_init (new_a, error=error)
       IF (atr) THEN
          new_row_dist => dbcsr_distribution_row_dist (a_dist)
          new_row_size => cp_dbcsr_row_block_sizes(matrix_a)
          CALL match_1_dist (new_col_dist,&
               dbcsr_distribution_col_dist (a_dist),&
               dbcsr_distribution_row_dist (c_dist))
          create_col_dist = .TRUE.
          new_col_size => cp_dbcsr_row_block_sizes(matrix_c)
       ELSE
          CALL match_1_dist (new_row_dist,&
               dbcsr_distribution_row_dist (a_dist),&
               dbcsr_distribution_row_dist (c_dist))
          create_row_dist = .TRUE.
          new_row_size => cp_dbcsr_row_block_sizes(matrix_c)
          new_col_dist => dbcsr_distribution_col_dist (a_dist)
          new_col_size => cp_dbcsr_col_block_sizes(matrix_a)
       ENDIF
       CALL dbcsr_distribution_new (new_dist, dbcsr_distribution_mp (a_dist),&
                                    new_row_dist, new_col_dist)
       IF (create_col_dist) DEALLOCATE(new_col_dist)
       IF (create_row_dist) DEALLOCATE(new_row_dist)       
       CALL cp_dbcsr_create (new_a, template=matrix_a,&
            dist=new_dist,&
            row_blk_size = new_row_size,&
            col_blk_size = new_col_size,&
            error=error)
       CALL dbcsr_distribution_release (new_dist)
       CALL cp_dbcsr_complete_redistribute (matrix_a, new_a, error=error)
    ELSE
       CALL cp_dbcsr_init (new_a, error=error)
       CALL cp_dbcsr_copy (new_a, matrix_a, shallow_data=.TRUE., error=error)
    ENDIF
  END SUBROUTINE matrix_match_sizes

! *****************************************************************************
!> \brief Matches a 1-D distribution
!>
!> Forms a new distribution new_dist of size template_dist.  As much data as
!> possible is taken from old_dist.
!> \param new_dist ...
!> \param old_dist ...
!> \param template_dist ...
! *****************************************************************************
  SUBROUTINE match_1_dist (new_dist, old_dist, template_dist)
    INTEGER, DIMENSION(:), INTENT(out), &
      POINTER                                :: new_dist
    INTEGER, DIMENSION(:), INTENT(in), &
      POINTER                                :: old_dist, template_dist

    INTEGER                                  :: i, max_bin, new_size, &
                                                old_size, sz

    old_size = SIZE(old_dist)
    new_size = SIZE(template_dist)
    sz = MIN(old_size, new_size)
    ALLOCATE(new_dist(new_size))
    new_dist = template_dist
    new_dist(1:sz) = old_dist(1:sz)
    max_bin = MAXVAL (old_dist)
    IF (max_bin .GT. 0) THEN
       DO i = sz+1 , new_size
          new_dist(i) = MOD (i, max_bin)
       END DO
    ELSE
       new_dist(sz+1:new_size) = 0
    ENDIF
  END SUBROUTINE match_1_dist

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval row_blk_sizes ...
! *****************************************************************************
  FUNCTION cp_dbcsr_row_block_sizes (matrix) RESULT (row_blk_sizes)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, DIMENSION(:), POINTER           :: row_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_row_block_sizes', &
      routineP = moduleN//':'//routineN

    row_blk_sizes => array_data(matrix%matrix%m%row_blk_size)
  END FUNCTION cp_dbcsr_row_block_sizes

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval col_blk_sizes ...
! *****************************************************************************
  FUNCTION cp_dbcsr_col_block_sizes (matrix) RESULT (col_blk_sizes)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, DIMENSION(:), POINTER           :: col_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_col_block_sizes', &
      routineP = moduleN//':'//routineN

    col_blk_sizes => array_data(matrix%matrix%m%col_blk_size)
  END FUNCTION cp_dbcsr_col_block_sizes

! *****************************************************************************
!> \brief Deallocates a DBCSR matrix for compatibility with CP2K
!> \param[in,out] matrix      DBCSR matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE cp_dbcsr_deallocate_matrix(matrix, error)
    TYPE(cp_dbcsr_type), POINTER             :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_deallocate_matrix', &
      routineP = moduleN//':'//routineN

    CALL cp_dbcsr_release (matrix, error=error)
    CALL cp_assert (.NOT. cp_dbcsr_valid_index(matrix), cp_warning_level,&
         cp_caller_error, routineN,&
         'You should not "deallocate" a referenced matrix. '//&
         'Avoid pointers to DBCSR matrices.')
    DEALLOCATE (matrix)
    NULLIFY (matrix)

  END SUBROUTINE cp_dbcsr_deallocate_matrix

! *****************************************************************************
! CP2k-compatible matrix sets
! *****************************************************************************


! *****************************************************************************
!> \brief   Allocate and initialize a real matrix 1-dimensional set.
!> \param[in,out] matrix_set  Set containing the DBCSR matrices
!> \param[in] nmatrix         Size of set
!> \param[in,out] error       cp2k error
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! *****************************************************************************
  SUBROUTINE allocate_dbcsr_matrix_set(matrix_set, nmatrix, error)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_set
    INTEGER, INTENT(IN)                      :: nmatrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_dbcsr_matrix_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix

    IF (ASSOCIATED(matrix_set)) CALL cp_dbcsr_deallocate_matrix_set(matrix_set,error=error)
    ALLOCATE (matrix_set(nmatrix))
    DO imatrix=1,nmatrix
       NULLIFY (matrix_set(imatrix)%matrix)
       !ALLOCATE (matrix_set(imatrix)%matrix)
       !CALL cp_dbcsr_init (matrix_set(imatrix)%matrix)
    END DO
  END SUBROUTINE allocate_dbcsr_matrix_set

! *****************************************************************************
!> \brief   Allocate and initialize a real matrix 2-dimensional set.
!> \param[in,out] matrix_set  Set containing the DBCSR matrix pointer type
!> \param[in] nmatrix         Size of set
!> \param mmatrix ...
!> \param[in,out] error       cp2k error
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! *****************************************************************************
  SUBROUTINE allocate_dbcsr_matrix_set_2d(matrix_set,nmatrix,mmatrix,error)
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrix_set
    INTEGER, INTENT(IN)                      :: nmatrix, mmatrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_dbcsr_matrix_set_2d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix, jmatrix

    IF (ASSOCIATED(matrix_set)) CALL cp_dbcsr_deallocate_matrix_set(matrix_set,error=error)
    ALLOCATE (matrix_set(nmatrix,mmatrix))
    DO jmatrix=1,mmatrix
      DO imatrix=1,nmatrix
         NULLIFY (matrix_set(imatrix,jmatrix)%matrix)
         !ALLOCATE (matrix_set(imatrix,jmatrix)%matrix)
         !CALL cp_dbcsr_init (matrix_set(imatrix,jmatrix)%matrix)
      END DO
    END DO
  END SUBROUTINE allocate_dbcsr_matrix_set_2d


! *****************************************************************************
!> \brief Deallocate a real matrix set and release all of the member matrices.
!> \param[in,out] matrix_set  Set containing the DBCSR matrix pointer type
!> \param[in,out] error       cp2k error
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! *****************************************************************************
  SUBROUTINE deallocate_dbcsr_matrix_set(matrix_set,error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_dbcsr_matrix_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix

    IF (ASSOCIATED(matrix_set)) THEN
      DO imatrix=1,SIZE(matrix_set)
        CALL cp_dbcsr_deallocate_matrix(matrix_set(imatrix)%matrix,error=error)
      END DO
      DEALLOCATE (matrix_set)
    END IF

  END SUBROUTINE deallocate_dbcsr_matrix_set

! *****************************************************************************
!> \brief Deallocate a real matrix set and release all of the member matrices.
!> \param[in,out] matrix_set  Set containing the DBCSR matrix pointer type
!> \param[in,out] error       cp2k error
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! *****************************************************************************
  SUBROUTINE deallocate_dbcsr_matrix_set_2d(matrix_set,error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrix_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'deallocate_dbcsr_matrix_set_2d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix, jmatrix

    IF (ASSOCIATED(matrix_set)) THEN
      DO jmatrix=1,SIZE(matrix_set,2)
        DO imatrix=1,SIZE(matrix_set,1)
          CALL cp_dbcsr_deallocate_matrix(matrix_set(imatrix,jmatrix)%matrix,error=error)
        END DO
      END DO
      DEALLOCATE (matrix_set)
    END IF
  END SUBROUTINE deallocate_dbcsr_matrix_set_2d

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval has_symmetry ...
! *****************************************************************************
  PURE FUNCTION cp_dbcsr_has_symmetry (matrix) RESULT (has_symmetry)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    LOGICAL                                  :: has_symmetry

    has_symmetry = dbcsr_has_symmetry(matrix%matrix)
  END FUNCTION cp_dbcsr_has_symmetry

! *****************************************************************************
!> \brief ...
!> \param dbcsr_mat ...
!> \param csr_mat ...
!> \param dist_format ...
!> \param csr_sparsity ...
!> \param numnodes ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_csr_create_from_dbcsr(dbcsr_mat, csr_mat, dist_format, csr_sparsity, numnodes, error)

    TYPE(cp_dbcsr_type), INTENT(IN)          :: dbcsr_mat
    TYPE(csr_type), INTENT(OUT)              :: csr_mat
    INTEGER                                  :: dist_format
    TYPE(cp_dbcsr_type), INTENT(IN), &
      OPTIONAL                               :: csr_sparsity
    INTEGER, INTENT(IN), OPTIONAL            :: numnodes
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_csr_create_from_dbcsr', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    IF (PRESENT(csr_sparsity)) THEN
      CALL csr_create_from_dbcsr(dbcsr_mat%matrix, csr_mat, dist_format, &
                                 csr_sparsity%matrix, numnodes, dbcsr_error)
    ELSE
      CALL csr_create_from_dbcsr(dbcsr_mat%matrix, csr_mat, dist_format, &
                                 numnodes=numnodes, error=dbcsr_error)
    ENDIF
  
  END SUBROUTINE cp_csr_create_from_dbcsr

! *****************************************************************************
!> \brief ...
!> \param dbcsr_mat ...
!> \param csr_mat ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_convert_csr_to_dbcsr(dbcsr_mat, csr_mat, error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: dbcsr_mat
    TYPE(csr_type), INTENT(INOUT)            :: csr_mat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_convert_csr_to_dbcsr', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL convert_csr_to_dbcsr(dbcsr_mat%matrix, csr_mat, dbcsr_error)

  END SUBROUTINE cp_convert_csr_to_dbcsr

! *****************************************************************************
!> \brief ...
!> \param dbcsr_mat ...
!> \param csr_mat ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_convert_dbcsr_to_csr(dbcsr_mat, csr_mat, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: dbcsr_mat
    TYPE(csr_type), INTENT(INOUT)            :: csr_mat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_convert_dbcsr_to_csr', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL convert_dbcsr_to_csr(dbcsr_mat%matrix, csr_mat, dbcsr_error)

  END SUBROUTINE cp_convert_dbcsr_to_csr

! *****************************************************************************
!> \brief ...
!> \param csr_mat ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_csr_destroy(csr_mat, error)
    TYPE(csr_type), INTENT(INOUT)            :: csr_mat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_csr_destroy', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL csr_destroy(csr_mat, dbcsr_error)

  END SUBROUTINE cp_csr_destroy


! *****************************************************************************
!> \brief ...
!> \param csr_mat ...
!> \param nrows_total ...
!> \param ncols_total ...
!> \param nze_total ...
!> \param nze_local ...
!> \param nrows_local ...
!> \param mp_group ...
!> \param data_type ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_csr_create_new(csr_mat, nrows_total, ncols_total, nze_total,&
                               nze_local, nrows_local, mp_group, data_type, error)

    TYPE(csr_type), INTENT(OUT)              :: csr_mat
    INTEGER, INTENT(IN)                      :: nrows_total, ncols_total
    INTEGER(KIND=int_8)                      :: nze_total
    INTEGER, INTENT(IN)                      :: nze_local, nrows_local, &
                                                mp_group
    INTEGER, INTENT(IN), OPTIONAL            :: data_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_csr_create_new', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL csr_create_new(csr_mat, nrows_total, ncols_total, nze_total,&
                        nze_local, nrows_local, mp_group, data_type, dbcsr_error)
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param matrix_b ...
!> \param matrix_a ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_csr_create_template(matrix_b, matrix_a, error)
    TYPE(csr_type), INTENT(OUT)              :: matrix_b
    TYPE(csr_type), INTENT(IN)               :: matrix_a
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_csr_create_template', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL csr_create_template(matrix_b, matrix_a, dbcsr_error)

  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param csr_mat ...
!> \param unit_nr ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE cp_csr_print_sparsity(csr_mat, unit_nr, error)
    TYPE(csr_type), INTENT(IN)               :: csr_mat
    INTEGER, INTENT(IN)                      :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL csr_print_sparsity(csr_mat, unit_nr, dbcsr_error)

  END SUBROUTINE 

! *****************************************************************************
!> \brief Emulation of sparse_matrix_types/add_block_node mapped
!>        to add_real_matrix_block.... should not be used any longer
!>
!> It adds a block to the dbcsr matrix and returns a rank-2 pointer to the
!> block. Currently it only and always uses the mutable data.
!> \param[in,out] matrix      DBCSR matrix
!> \param[in]  block_row      the row
!> \param[in]  block_col      the column
!> \param[in]  block          the block to put
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE cp_dbcsr_add_block_node (matrix, block_row, block_col, block, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_dbcsr_add_block_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: c, ithread, p, r
    LOGICAL                                  :: dbg, existed, is_there, tr

!   ---------------------------------------------------------------------------

    dbg = .FALSE.

    ithread = 0
!$  ithread = omp_get_thread_num()
    IF (.NOT. ASSOCIATED (matrix%matrix%m%wms)) THEN
       CALL cp_dbcsr_work_create (matrix, work_mutable=.TRUE., error=error)
       matrix%matrix%m%valid = .FALSE.
    ENDIF
!$  CALL cp_assert (SIZE (matrix%matrix%m%wms) .GE. omp_get_num_threads(),&
!$       cp_fatal_level, cp_wrong_args_error, routineN,&
!$       "Too few threads.", error=error)
    CALL cp_assert (dbcsr_wm_use_mutable (matrix%matrix%m%wms(ithread+1)),&
         cp_warning_level,&
         cp_unimplemented_error_nr, routineN,&
         "Data loss due to no conversion of appendable to mutable data")
    is_there = ASSOCIATED(block)
    !r = row ; c = col ; tr = .FALSE.
    !CALL dbcsr_get_stored_coordinates (matrix, r, c, tr)
    !CALL dbcsr_reserve_block2d (matrix, row, col, block)
    !write(*,*) 'add_block_node: block_row',block_row,' block_col',block_col
    CALL cp_dbcsr_reserve_block2d (matrix, block_row, block_col, block,&
         existed=existed)
!
    IF (dbg) THEN
       r = block_row ; c = block_col ; tr = .FALSE.
       CALL cp_dbcsr_get_stored_coordinates (matrix, r, c, tr, p)
       CALL cp_assert (p .EQ. dbcsr_mp_mynode (dbcsr_distribution_mp (&
            cp_dbcsr_distribution(matrix))),&
            cp_warning_level, cp_internal_error, routineN,&
            "Adding non-local element", error=error)
    ENDIF
    CALL cp_assert (.NOT.existed, cp_warning_level, cp_wrong_args_error,&
         routineN, "You should not add existing blocks according to old API.")
    IF(.NOT.is_there) block(:,:) = 0.0_dp
  END SUBROUTINE cp_dbcsr_add_block_node

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param row ...
!> \param col ...
!> \param block ...
!> \param summation ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE cp_dbcsr_put_block_area (matrix, row, col, block,&
       summation, scale)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    TYPE(dbcsr_data_obj)                     :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: summation
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_put_block_area', &
      routineP = moduleN//':'//routineN

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL dbcsr_put_block(matrix%matrix, row, col, block,&
                         summation=summation, scale=scale)

  END SUBROUTINE cp_dbcsr_put_block_area

#include "cp_dbcsr_interface_s.f90"
#include "cp_dbcsr_interface_d.f90"
#include "cp_dbcsr_interface_c.f90"
#include "cp_dbcsr_interface_z.f90"

END MODULE cp_dbcsr_interface

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE array_types,                     ONLY: array_data,&
                                             array_hold,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_release,&
                                             array_size
  USE cp_dbcsr_methods,                ONLY: cp_dbcsr_col_block_offsets,&
                                             cp_dbcsr_col_block_sizes,&
                                             cp_dbcsr_row_block_offsets,&
                                             cp_dbcsr_row_block_sizes
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_iterator,&
                                             cp_dbcsr_type
  USE dbcsr_api,                       ONLY: &
       create_bl_distribution, dbcsr_add, dbcsr_add_on_diag, &
       dbcsr_binary_read, dbcsr_binary_write, dbcsr_blk_col_offset, &
       dbcsr_blk_column_size, dbcsr_blk_row_offset, dbcsr_blk_row_size, &
       dbcsr_btriu, dbcsr_checksum, dbcsr_clear_mempools, &
       dbcsr_complete_redistribute, dbcsr_copy, dbcsr_copy_into_existing, &
       dbcsr_create, dbcsr_data_obj, dbcsr_desymmetrize_deep, &
       dbcsr_distribute, dbcsr_distribution, dbcsr_distribution_col_dist, &
       dbcsr_distribution_mp, dbcsr_distribution_new, dbcsr_distribution_obj, &
       dbcsr_distribution_release, dbcsr_distribution_row_dist, &
       dbcsr_error_type, dbcsr_filter, dbcsr_finalize, dbcsr_finalize_lib, &
       dbcsr_frobenius_norm, dbcsr_function_of_elements, &
       dbcsr_gershgorin_norm, dbcsr_get_block, dbcsr_get_block_diag, &
       dbcsr_get_block_p, dbcsr_get_data_size, dbcsr_get_data_type, &
       dbcsr_get_diag, dbcsr_get_info, dbcsr_get_matrix_type, &
       dbcsr_get_num_blocks, dbcsr_get_occupation, &
       dbcsr_get_stored_coordinates, dbcsr_hadamard_product, &
       dbcsr_index_compact, dbcsr_init, dbcsr_init_lib, &
       dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
       dbcsr_iterator_start, dbcsr_iterator_stop, &
       dbcsr_make_untransposed_blocks, dbcsr_max_col_size, &
       dbcsr_max_row_size, dbcsr_maxabs, dbcsr_mp_new, dbcsr_mp_obj, &
       dbcsr_mp_release, dbcsr_multiply, dbcsr_name, dbcsr_nblkcols_local, &
       dbcsr_nblkcols_total, dbcsr_nblkrows_local, dbcsr_nblkrows_total, &
       dbcsr_new_transposed, dbcsr_nfullcols_local, dbcsr_nfullcols_total, &
       dbcsr_nfullrows_local, dbcsr_nfullrows_total, dbcsr_no_transpose, &
       dbcsr_norm, dbcsr_print, dbcsr_print_block_sum, dbcsr_put_block, &
       dbcsr_release, dbcsr_replicate, dbcsr_replicate_all, &
       dbcsr_reserve_all_blocks, dbcsr_reserve_block2d, dbcsr_reserve_blocks, &
       dbcsr_reserve_diag_blocks, dbcsr_scalar, dbcsr_scalar_fill_all, &
       dbcsr_scalar_get_type, dbcsr_scalar_get_value, dbcsr_scalar_set_type, &
       dbcsr_scalar_type, dbcsr_scalar_zero, dbcsr_scale, &
       dbcsr_scale_by_vector, dbcsr_scale_mat, dbcsr_set, dbcsr_set_diag, &
       dbcsr_sum_replicated, dbcsr_trace, dbcsr_triu, dbcsr_type_complex_4, &
       dbcsr_type_complex_8, dbcsr_valid_index, dbcsr_verify_matrix, &
       dbcsr_work_create, swap
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop

  !$ USE OMP_LIB
#include "common/cp_common_uses.f90"
#define CHECKUSAGE CALL cp_assert(.FALSE.,cp_fatal_level,cp_unimplemented_error_nr,routineN,"OBS")
#define BUFFERME ! buffering code

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_interface'

  !
  ! Interface to libdbcsr (contains cp2k timing and error)
  PUBLIC :: cp_dbcsr_conform_scalar
  PUBLIC :: cp_dbcsr_init_lib, cp_dbcsr_finalize_lib
  PUBLIC :: cp_dbcsr_untranspose_blocks
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
  PUBLIC :: cp_dbcsr_create
  PUBLIC :: cp_dbcsr_work_create
  PUBLIC :: cp_dbcsr_verify_matrix
  PUBLIC :: cp_dbcsr_btriu
  PUBLIC :: cp_dbcsr_sum_replicated
  PUBLIC :: cp_dbcsr_checksum
  PUBLIC :: cp_dbcsr_trace
  PUBLIC :: cp_dbcsr_print
  PUBLIC :: cp_dbcsr_init
  PUBLIC :: cp_dbcsr_init_p
  PUBLIC :: cp_dbcsr_release
  PUBLIC :: cp_dbcsr_release_p
  PUBLIC :: cp_dbcsr_distribute
  PUBLIC :: cp_dbcsr_replicate_all
  PUBLIC :: cp_dbcsr_replicate
  PUBLIC :: cp_dbcsr_norm
  PUBLIC :: cp_dbcsr_get_info
  PUBLIC :: cp_dbcsr_get_block
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
  PUBLIC :: cp_dbcsr_nblkrows_local
  PUBLIC :: cp_dbcsr_nblkcols_local
  PUBLIC :: cp_dbcsr_nblkrows_total
  PUBLIC :: cp_dbcsr_nblkcols_total
  PUBLIC :: cp_dbcsr_get_num_blocks
  PUBLIC :: cp_dbcsr_get_data_size
  PUBLIC :: cp_dbcsr_col_block_offsets
  PUBLIC :: cp_dbcsr_row_block_offsets
  PUBLIC :: cp_dbcsr_nfullrows_total
  PUBLIC :: cp_dbcsr_nfullcols_total
  PUBLIC :: cp_dbcsr_nfullrows_local
  PUBLIC :: cp_dbcsr_nfullcols_local
  !         1234567890123456789012345678901
  PUBLIC :: cp_dbcsr_get_stored_coordinates
  PUBLIC :: cp_dbcsr_valid_index
  PUBLIC :: cp_dbcsr_get_data_type
  PUBLIC :: cp_dbcsr_reserve_block2d
  PUBLIC :: cp_dbcsr_reserve_blocks
  PUBLIC :: cp_dbcsr_reserve_all_blocks
  PUBLIC :: cp_dbcsr_reserve_diag_blocks
  PUBLIC :: cp_dbcsr_complete_redistribute
  PUBLIC :: cp_dbcsr_gershgorin_norm
  PUBLIC :: cp_dbcsr_frobenius_norm
  PUBLIC :: cp_dbcsr_maxabs
  PUBLIC :: cp_dbcsr_name
  PUBLIC :: cp_dbcsr_binary_write
  PUBLIC :: cp_dbcsr_binary_read
  PUBLIC :: cp_dbcsr_function_of_elements
  PUBLIC :: cp_dbcsr_print_block_sum
  PUBLIC :: cp_dbcsr_triu
  PUBLIC :: cp_dbcsr_clear_mempools

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
     MODULE PROCEDURE cp_dbcsr_scale_d_m,&
                      cp_dbcsr_scale_s_m,&
                      cp_dbcsr_scale_c_m,&
                      cp_dbcsr_scale_z_m
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

  INTERFACE cp_dbcsr_get_block
     MODULE PROCEDURE cp_dbcsr_get_block_d,&
                      cp_dbcsr_get_block_s,&
                      cp_dbcsr_get_block_z,&
                      cp_dbcsr_get_block_c
     MODULE PROCEDURE cp_dbcsr_get_2d_block_d,&
                      cp_dbcsr_get_2d_block_s,&
                      cp_dbcsr_get_2d_block_z,&
                      cp_dbcsr_get_2d_block_c
  END INTERFACE

  INTERFACE cp_dbcsr_put_block
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
     MODULE PROCEDURE cp_iterator_next_block_index
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

  PRIVATE

CONTAINS

! *****************************************************************************
!> \brief Initializes DBCSR
! *****************************************************************************
  SUBROUTINE cp_dbcsr_init_lib (group, error)
    INTEGER, INTENT(IN)                      :: group
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_init_lib', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_init_lib (group, dbcsr_error)
  END SUBROUTINE cp_dbcsr_init_lib

! *****************************************************************************
!> \brief Finalizes DBCSR
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
!> \author vw
! *****************************************************************************

  FUNCTION cp_dbcsr_name (matrix) RESULT (name)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    CHARACTER(len=default_string_length)     :: name

    name = dbcsr_name (matrix%matrix)
  END FUNCTION cp_dbcsr_name


  FUNCTION cp_dbcsr_gershgorin_norm(matrix) RESULT (norm)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    REAL(KIND=real_8)                        :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_gershgorin_norm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_error_type)                      :: error

    CALL timeset(routineN,handle)

    CALL cp_error_init (error)
    norm = dbcsr_gershgorin_norm(matrix%matrix)

    CALL timestop(handle)

  END FUNCTION cp_dbcsr_gershgorin_norm

  FUNCTION cp_dbcsr_frobenius_norm(matrix, local) RESULT (norm)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    LOGICAL, INTENT(in), OPTIONAL            :: local
    REAL(KIND=real_8)                        :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_frobenius_norm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_error_type)                      :: error

    CALL timeset(routineN,handle)

    CALL cp_error_init (error)
    norm = dbcsr_frobenius_norm(matrix%matrix, local)

    CALL timestop(handle)
  END FUNCTION cp_dbcsr_frobenius_norm

  FUNCTION cp_dbcsr_maxabs(matrix) RESULT (norm)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    REAL(KIND=real_8)                        :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_maxabs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_error_type)                      :: error

    CALL timeset(routineN,handle)

    CALL cp_error_init (error)
    norm = dbcsr_maxabs(matrix%matrix)

    CALL timestop(handle)
  END FUNCTION cp_dbcsr_maxabs



  SUBROUTINE cp_dbcsr_complete_redistribute(matrix, redist, keep_sparsity, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: redist
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_complete_redistribute', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN,handle)

    CALL dbcsr_complete_redistribute(matrix%matrix, redist%matrix, keep_sparsity, dbcsr_error)

    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_complete_redistribute


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

  SUBROUTINE cp_dbcsr_reserve_all_blocks(matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_reserve_all_blocks', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_reserve_all_blocks(matrix%matrix, dbcsr_error)
  END SUBROUTINE cp_dbcsr_reserve_all_blocks

  SUBROUTINE cp_dbcsr_reserve_diag_blocks(matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_reserve_diag_blocks', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_reserve_diag_blocks(matrix%matrix, dbcsr_error)
  END SUBROUTINE cp_dbcsr_reserve_diag_blocks

  PURE FUNCTION cp_dbcsr_get_data_type (matrix) RESULT (data_type)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: data_type

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_data_type', &
      routineP = moduleN//':'//routineN

    data_type = dbcsr_get_data_type (matrix%matrix)

  END FUNCTION cp_dbcsr_get_data_type

  PURE FUNCTION cp_dbcsr_valid_index (matrix) RESULT (valid_index)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    LOGICAL                                  :: valid_index

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_valid_index', &
      routineP = moduleN//':'//routineN

    valid_index = dbcsr_valid_index(matrix%matrix)

  END FUNCTION cp_dbcsr_valid_index

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

  PURE FUNCTION cp_dbcsr_get_num_blocks (matrix) RESULT (num_blocks)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: num_blocks

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_num_blocks', &
      routineP = moduleN//':'//routineN

    num_blocks = dbcsr_get_num_blocks (matrix%matrix)

  END FUNCTION cp_dbcsr_get_num_blocks

  FUNCTION cp_dbcsr_get_data_size (matrix) RESULT (data_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: data_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_data_size', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    data_size = dbcsr_get_data_size(matrix%matrix)

    CALL timestop(handle)

  END FUNCTION cp_dbcsr_get_data_size

  PURE FUNCTION cp_dbcsr_get_matrix_type (matrix) RESULT (matrix_type)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    CHARACTER                                :: matrix_type

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_matrix_type', &
      routineP = moduleN//':'//routineN

    matrix_type = dbcsr_get_matrix_type (matrix%matrix)

  END FUNCTION cp_dbcsr_get_matrix_type

  FUNCTION cp_dbcsr_get_occupation (matrix) RESULT (occupation)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    REAL(KIND=real_8)                        :: occupation

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_occupation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    occupation = dbcsr_get_occupation (matrix%matrix)

    CALL timestop(handle)

  END FUNCTION cp_dbcsr_get_occupation


  FUNCTION cp_dbcsr_nblkrows_total(matrix) RESULT (nblkrows_total)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nblkrows_total

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nblkrows_total', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    nblkrows_total = dbcsr_nblkrows_total(matrix%matrix)

    CALL timestop(handle)
  END FUNCTION cp_dbcsr_nblkrows_total

  FUNCTION cp_dbcsr_nblkcols_total(matrix) RESULT (nblkcols_total)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nblkcols_total

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nblkcols_total', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    nblkcols_total = dbcsr_nblkcols_total(matrix%matrix)
    CALL timestop(handle)
  END FUNCTION cp_dbcsr_nblkcols_total

  FUNCTION cp_dbcsr_nfullrows_total(matrix) RESULT (nfullrows_total)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nfullrows_total

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nfullrows_total', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    nfullrows_total = dbcsr_nfullrows_total(matrix%matrix)

    CALL timestop(handle)
  END FUNCTION cp_dbcsr_nfullrows_total

  FUNCTION cp_dbcsr_nfullcols_total(matrix) RESULT (nfullcols_total)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nfullcols_total

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nfullcols_total', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    nfullcols_total = dbcsr_nfullcols_total(matrix%matrix)
    CALL timestop(handle)
  END FUNCTION cp_dbcsr_nfullcols_total

  FUNCTION cp_dbcsr_nblkrows_local(matrix) RESULT (nblkrows_local)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nblkrows_local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nblkrows_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    nblkrows_local = dbcsr_nblkrows_local(matrix%matrix)
    CALL timestop(handle)
  END FUNCTION cp_dbcsr_nblkrows_local

  FUNCTION cp_dbcsr_nblkcols_local(matrix) RESULT (nblkcols_local)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nblkcols_local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nblkcols_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    nblkcols_local = dbcsr_nblkcols_local(matrix%matrix)
    CALL timestop(handle)
  END FUNCTION cp_dbcsr_nblkcols_local

  FUNCTION cp_dbcsr_nfullrows_local(matrix) RESULT (nfullrows_local)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nfullrows_local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nfullrows_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    nfullrows_local = dbcsr_nfullrows_local(matrix%matrix)

    CALL timestop(handle)

  END FUNCTION cp_dbcsr_nfullrows_local

  FUNCTION cp_dbcsr_nfullcols_local(matrix) RESULT (nfullcols_local)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nfullcols_local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nfullcols_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    nfullcols_local = dbcsr_nfullcols_local(matrix%matrix)

    CALL timestop(handle)

  END FUNCTION cp_dbcsr_nfullcols_local

  FUNCTION cp_dbcsr_max_row_size(matrix) RESULT (max_row_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: max_row_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_max_row_size', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    max_row_size = dbcsr_max_row_size(matrix%matrix)

    CALL timestop(handle)
  END FUNCTION cp_dbcsr_max_row_size

  FUNCTION cp_dbcsr_max_col_size(matrix) RESULT (max_col_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: max_col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_max_col_size', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    max_col_size = dbcsr_max_col_size(matrix%matrix)
    CALL timestop(handle)

  END FUNCTION cp_dbcsr_max_col_size

  FUNCTION cp_dbcsr_distribution (matrix) RESULT (distribution)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(dbcsr_distribution_obj)             :: distribution

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_distribution', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    distribution = dbcsr_distribution (matrix%matrix)
    CALL timestop(handle)

  END FUNCTION cp_dbcsr_distribution

  SUBROUTINE cp_create_bl_distribution (block_distribution,&
       block_size, nelements, nbins)
    TYPE(array_i1d_obj), INTENT(OUT)         :: block_distribution, block_size
    INTEGER, INTENT(IN)                      :: nelements, nbins

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_create_bl_distribution', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    CALL create_bl_distribution (block_distribution,&
       block_size, nelements, nbins)

    CALL timestop(handle)

  END SUBROUTINE cp_create_bl_distribution

  SUBROUTINE cp_dbcsr_distribution_release(dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_distribution_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    CALL dbcsr_distribution_release(dist)
    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_distribution_release

  PURE FUNCTION cp_dbcsr_iterator_blocks_left (iterator) RESULT (blocks_left)
    TYPE(cp_dbcsr_iterator), INTENT(IN)      :: iterator
    LOGICAL                                  :: blocks_left

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_iterator_blocks_left', &
      routineP = moduleN//':'//routineN

    blocks_left = dbcsr_iterator_blocks_left (iterator)

  END FUNCTION cp_dbcsr_iterator_blocks_left

  SUBROUTINE cp_dbcsr_mp_release(mp_env)
    TYPE(dbcsr_mp_obj), INTENT(INOUT)        :: mp_env

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_mp_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    CALL dbcsr_mp_release(mp_env)
    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_mp_release

  SUBROUTINE cp_dbcsr_mp_new(mp_env, pgrid, mp_group, mynode, numnodes, myprow,&
       mypcol)
    TYPE(dbcsr_mp_obj), INTENT(OUT)          :: mp_env
    INTEGER, DIMENSION(0:, 0:), INTENT(IN)   :: pgrid
    INTEGER, INTENT(IN)                      :: mp_group, mynode
    INTEGER, INTENT(IN), OPTIONAL            :: numnodes, myprow, mypcol

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_mp_new', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    CALL dbcsr_mp_new(mp_env, pgrid, mp_group, mynode, numnodes, myprow,&
       mypcol)

    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_mp_new

  SUBROUTINE cp_iterator_next_block_index (iterator, row, column, blk,&
       blk_p, row_size, col_size, row_offset, col_offset)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)   :: iterator
    INTEGER, INTENT(OUT)                     :: row, column, blk
    INTEGER, INTENT(OUT), OPTIONAL           :: blk_p, row_size, col_size, &
                                                row_offset, col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iterator_next_block_index', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: transposed

    CALL dbcsr_iterator_next_block (iterator, row, column, blk,&
       transposed, blk_p, row_size, col_size, row_offset, col_offset)

  END SUBROUTINE cp_iterator_next_block_index


  SUBROUTINE cp_dbcsr_iterator_stop (iterator)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)   :: iterator

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_iterator_stop', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    CALL dbcsr_iterator_stop (iterator)
    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_iterator_stop

  SUBROUTINE cp_dbcsr_iterator_start (iterator, matrix, shared, dynamic,&
       dynamic_byrows, contiguous_pointers, read_only)
    TYPE(cp_dbcsr_iterator), INTENT(OUT)     :: iterator
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: shared, dynamic, &
                                                dynamic_byrows, &
                                                contiguous_pointers, read_only

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_iterator_start', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_error_type)                      :: error

    CALL timeset(routineN,handle)
    CALL cp_error_init (error)
    CALL dbcsr_iterator_start (iterator, matrix%matrix, shared, dynamic,&
       dynamic_byrows, contiguous_pointers, read_only)
    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_iterator_start


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
    TYPE(array_i1d_obj), INTENT(OUT), &
      OPTIONAL                               :: row_blk_size, col_blk_size, &
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

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    CALL dbcsr_get_info(matrix%matrix, nblkrows_total, nblkcols_total,&
         nfullrows_total, nfullcols_total,&
         nblkrows_local, nblkcols_local,&
         nfullrows_local, nfullcols_local,&
         my_prow, my_pcol,&
         local_rows, local_cols, proc_row_dist, proc_col_dist,&
         row_blk_size, col_blk_size, row_blk_offset, col_blk_offset, distribution, name, data_area,&
         matrix_type, data_type)

    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_get_info

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

    INTEGER                                  :: handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN,handle)

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
    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_norm

  SUBROUTINE cp_dbcsr_replicate_all(matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_replicate_all', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN,handle)
    CALL dbcsr_replicate_all(matrix%matrix, dbcsr_error)
    CALL timestop(handle)
  END SUBROUTINE cp_dbcsr_replicate_all

  SUBROUTINE cp_dbcsr_replicate(matrix, replicate_rows, replicate_columns,&
       restrict_source, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    LOGICAL, INTENT(IN)                      :: replicate_rows, &
                                                replicate_columns
    INTEGER, INTENT(IN), OPTIONAL            :: restrict_source
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_replicate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN,handle)
    CALL dbcsr_replicate(matrix%matrix, replicate_rows, replicate_columns,&
         restrict_source, dbcsr_error)
    CALL timestop(handle)
  END SUBROUTINE cp_dbcsr_replicate

  SUBROUTINE cp_dbcsr_distribute(matrix, fast, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    LOGICAL, INTENT(in), OPTIONAL            :: fast
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_distribute', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN,handle)
    CALL dbcsr_distribute(matrix%matrix, fast, dbcsr_error)
    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_distribute

  SUBROUTINE cp_dbcsr_release_p (matrix, error)
    TYPE(cp_dbcsr_type), POINTER             :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_release_p', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    IF(ASSOCIATED(matrix)) THEN
       CALL cp_dbcsr_release (matrix, error)
       DEALLOCATE(matrix)
    ENDIF
    CALL timestop(handle)


  END SUBROUTINE cp_dbcsr_release_p

  SUBROUTINE cp_dbcsr_release (matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    CALL dbcsr_release (matrix%matrix)

    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_release

  SUBROUTINE cp_dbcsr_init (matrix, error)
    TYPE(cp_dbcsr_type), INTENT(OUT)         :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    CALL dbcsr_init(matrix%matrix)
    matrix%ref_count = 0

    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_init

  SUBROUTINE cp_dbcsr_init_p (matrix, error)
    TYPE(cp_dbcsr_type), POINTER             :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_init_p', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    IF(ASSOCIATED(matrix)) THEN
       CALL cp_dbcsr_release(matrix, error)
       DEALLOCATE(matrix)
    ENDIF

    ALLOCATE(matrix)
    CALL cp_dbcsr_init (matrix, error)

    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_init_p

  SUBROUTINE cp_dbcsr_print(matrix, nodata, matlab_format, variable_name, unit_nr, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: nodata, matlab_format
    CHARACTER(*), INTENT(in), OPTIONAL       :: variable_name
    INTEGER, OPTIONAL                        :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_print', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN,handle)

    CALL dbcsr_print(matrix%matrix, nodata, matlab_format, variable_name, unit_nr, dbcsr_error)
    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_print


  FUNCTION cp_dbcsr_checksum(matrix, local, pos, error) RESULT(checksum)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: local, pos
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    REAL(KIND=dp)                            :: checksum

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_checksum', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN,handle)
    checksum = dbcsr_checksum(matrix%matrix,&
         local=local, pos=pos, error=dbcsr_error)
    CALL timestop(handle)
  END FUNCTION cp_dbcsr_checksum

  SUBROUTINE cp_dbcsr_sum_replicated (matrix, error)
    TYPE(cp_dbcsr_type), INTENT(inout)       :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_sum_replicated', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN,handle)
    CALL dbcsr_sum_replicated (matrix%matrix, dbcsr_error)
    CALL timestop(handle)
  END SUBROUTINE cp_dbcsr_sum_replicated

  SUBROUTINE cp_dbcsr_btriu(matrix_b, matrix_a, error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_b
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_btriu', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL cp_dbcsr_create (matrix_b, name="BTRIU of"//cp_dbcsr_name (matrix_a),&
         template=matrix_a, error=error)
    CALL dbcsr_btriu(matrix_b%matrix, matrix_a%matrix, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_btriu

  SUBROUTINE cp_dbcsr_verify_matrix(m, error, verbosity, local)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: m
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER, INTENT(IN), OPTIONAL            :: verbosity
    LOGICAL, INTENT(IN), OPTIONAL            :: local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_verify_matrix', &
      r = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN,handle)
    CALL dbcsr_verify_matrix(m%matrix, verbosity, local, dbcsr_error)

    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_verify_matrix

  SUBROUTINE cp_dbcsr_finalize(matrix, resort, reshuffle, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: resort, reshuffle
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_finalize(matrix%matrix, resort, reshuffle, error=dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_finalize

  SUBROUTINE cp_dbcsr_work_create(matrix, nblks_guess, sizedata_guess, n,&
       error, work_mutable)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN), OPTIONAL            :: nblks_guess, sizedata_guess, n
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(in), OPTIONAL            :: work_mutable

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_work_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN,handle)

    CALL dbcsr_work_create(matrix%matrix, nblks_guess, sizedata_guess, n,&
         work_mutable, error=dbcsr_error)
    CALL timestop(handle)
  END SUBROUTINE cp_dbcsr_work_create

  SUBROUTINE cp_dbcsr_create_new(matrix, name, dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, data_type, reuse,&
       mutable_work, replication_type, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    CHARACTER(len=*), INTENT(IN)             :: name
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    CHARACTER, INTENT(IN)                    :: matrix_type
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_size, col_blk_size
    INTEGER, INTENT(IN), OPTIONAL            :: nblks, nze, data_type
    LOGICAL, INTENT(IN), OPTIONAL            :: reuse, mutable_work
    CHARACTER, INTENT(IN), OPTIONAL          :: replication_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_create_new', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset (routineN, handle)

    CALL dbcsr_create(matrix%matrix, name, dist,&
         matrix_type,&
         row_blk_size, col_blk_size, nblks, nze, data_type, reuse=reuse,&
         mutable_work=mutable_work, replication_type=replication_type,&
         error=dbcsr_error)

    matrix%ref_count = 1

    CALL timestop (handle)
  END SUBROUTINE cp_dbcsr_create_new

  SUBROUTINE cp_dbcsr_create_template(matrix, name, template,&
       dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, data_type, reuse,&
       mutable_work, replication_type, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: name
    TYPE(cp_dbcsr_type), INTENT(IN)          :: template
    TYPE(dbcsr_distribution_obj), &
      INTENT(IN), OPTIONAL                   :: dist
    CHARACTER, INTENT(IN), OPTIONAL          :: matrix_type
    TYPE(array_i1d_obj), INTENT(IN), &
      OPTIONAL                               :: row_blk_size, col_blk_size
    INTEGER, INTENT(IN), OPTIONAL            :: nblks, nze, data_type
    LOGICAL, INTENT(IN), OPTIONAL            :: reuse, mutable_work
    CHARACTER, INTENT(IN), OPTIONAL          :: replication_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_create_template', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset (routineN, handle)
    CALL dbcsr_create(matrix%matrix, template%matrix,&
         name=name, dist=dist, matrix_type=matrix_type,&
         row_blk_size=row_blk_size, col_blk_size=col_blk_size,&
         nblks=nblks, nze=nze, data_type=data_type,&
         reuse=reuse, mutable_work=mutable_work,&
         replication_type=replication_type,&
         error=dbcsr_error)
    !
    matrix%ref_count = 1
    CALL timestop (handle)
  END SUBROUTINE cp_dbcsr_create_template


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

    INTEGER                                  :: timing_handle
    LOGICAL                                  :: quick
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    IF (PRESENT (thorough)) THEN
       quick = .NOT. thorough
    ELSE
       quick = .FALSE.
    ENDIF
    CALL dbcsr_filter(matrix%matrix, cp_dbcsr_conform_scalar (eps, matrix, error),&
                      method, use_absolute, filter_diag,&
         quick=quick, error=dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_filter

  SUBROUTINE cp_dbcsr_set_diag(matrix, diag, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    REAL(dp), DIMENSION(:), INTENT(IN)       :: diag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_set_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_set_diag(matrix%matrix, diag, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_set_diag

  SUBROUTINE cp_dbcsr_get_diag(matrix, diag, error)

    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: diag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_get_diag(matrix%matrix, diag, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_get_diag

  SUBROUTINE cp_dbcsr_get_block_diag(matrix, diag, error)

    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: diag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_block_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL cp_dbcsr_create (diag, "Diagonal of "//TRIM(cp_dbcsr_name (matrix)),&
         template=matrix,&
         error=error)
    CALL dbcsr_get_block_diag(matrix%matrix, diag%matrix, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_get_block_diag

  SUBROUTINE cp_dbcsr_add_on_diag(matrix, alpha_scalar, first_row, last_row, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    REAL(KIND=dp), INTENT(IN)                :: alpha_scalar
    INTEGER, INTENT(in), OPTIONAL            :: first_row, last_row
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_add_on_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_add_on_diag(matrix%matrix,&
         cp_dbcsr_conform_scalar (alpha_scalar, matrix, error),&
         first_row, last_row, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_add_on_diag

  SUBROUTINE cp_dbcsr_binary_write(matrix, filepath, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    CHARACTER(LEN=*), INTENT(IN)             :: filepath
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_binary_write', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_binary_write(matrix%matrix, filepath, dbcsr_error)
    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_binary_write

  SUBROUTINE cp_dbcsr_binary_read(filepath, distribution, groupid, matrix_new, error)
    CHARACTER(len=*), INTENT(IN)             :: filepath
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: distribution
    INTEGER, INTENT(IN), OPTIONAL            :: groupid
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_new
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_binary_read', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_binary_read(filepath, distribution, groupid, matrix_new%matrix, dbcsr_error)
    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_binary_read

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

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    !call cp_assert (matrix_b%ref_count .gt. 0,&
    !     cp_warning_level, cp_caller_error, routineN,&
    !     "Matrix not created.", error=error)
    IF (matrix_b%ref_count .EQ. 0) THEN
       CALL cp_dbcsr_create (matrix_b, template=matrix_a,&
            error=error)
    ENDIF
    CALL dbcsr_copy(matrix_b%matrix, matrix_a%matrix, name, keep_sparsity,&
         shallow_data, keep_imaginary, matrix_type, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_copy

  SUBROUTINE cp_dbcsr_copy_into_existing(matrix_b, matrix_a, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_b
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_copy_into_existing', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    IF (matrix_b%ref_count .EQ. 0) THEN
       CALL cp_dbcsr_create (matrix_b, template=matrix_a,&
            error=error)
    ENDIF
    CALL dbcsr_copy_into_existing(matrix_b%matrix, matrix_a%matrix, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_copy_into_existing


  SUBROUTINE cp_dbcsr_desymmetrize(matrix_a, matrix_b, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_b
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_desymmetrize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    !call cp_assert (matrix_b%ref_count .gt. 0,&
    !     cp_warning_level, cp_caller_error, routineN,&
    !     "Matrix not created.", error=error)
    IF (matrix_b%ref_count .EQ. 0) THEN
       CALL cp_dbcsr_create (matrix_b, template=matrix_a,&
            error=error)
    ENDIF
    CALL dbcsr_desymmetrize_deep(matrix_a%matrix, matrix_b%matrix,&
         untransposed_data = .TRUE., error=dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_desymmetrize


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

    INTEGER                                  :: timing_handle
    LOGICAL                                  :: myshallow_data_copy, &
                                                mytranspose_distribution
    TYPE(dbcsr_distribution_obj)             :: myuse_distribution
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
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
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_transposed

  SUBROUTINE cp_dbcsr_function_of_elements(matrix_a, func, a0, a1, a2, a3, a4,& 
    a5, a6, a7, a8, a9, error) 
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    INTEGER, INTENT(IN)                      :: func
    REAL(kind=dp), INTENT(IN), OPTIONAL      :: a0, a1, a2, a3, a4, a5, a6, &
                                                a7, a8, a9
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_function_of_elements', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle) 
    CALL dbcsr_function_of_elements(matrix_a%matrix,func,a0,a1,a2,a3,a4,a5,a6,& 
            a7,a8,a9,dbcsr_error) 
    CALL timestop(timing_handle) 
  END SUBROUTINE cp_dbcsr_function_of_elements 
 
  SUBROUTINE cp_dbcsr_hadamard_product(matrix_a, matrix_b, matrix_c, b_assume_value,& 
    error) 
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a, matrix_b
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_c
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: b_assume_value
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_hadamard_product', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_hadamard_product(matrix_a%matrix, matrix_b%matrix, matrix_c%matrix,&
            b_assume_value, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_hadamard_product

  ! ultimately this routine should not be needed
  SUBROUTINE cp_dbcsr_untranspose_blocks(matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_untranspose_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_make_untransposed_blocks(matrix%matrix,dbcsr_error)
    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_untranspose_blocks


  SUBROUTINE fill_sizes (matrix, row, col,&
       row_size, col_size, row_offset, col_offset)
    TYPE(cp_dbcsr_type), INTENT(in)          :: matrix
    INTEGER, INTENT(in)                      :: row, col
    INTEGER, INTENT(out), OPTIONAL           :: row_size, col_size, &
                                                row_offset, col_offset

    IF (PRESENT (row_size)) THEN
       row_size = dbcsr_blk_row_size (matrix%matrix, row)
    ENDIF
    IF (PRESENT (col_size)) THEN
       col_size = dbcsr_blk_column_size (matrix%matrix, col)
    ENDIF
    IF (PRESENT (row_offset)) THEN
       row_offset = dbcsr_blk_row_offset (matrix%matrix, row)
    ENDIF
    IF (PRESENT (col_offset)) THEN
       col_offset = dbcsr_blk_col_offset (matrix%matrix, col)
    ENDIF
  END SUBROUTINE fill_sizes


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
!>
!> Some copies could be avoided with more code.
! *****************************************************************************
  SUBROUTINE matrix_match_sizes (matrix_c, matrix_a, tr_a, matrix_b, tr_b,&
       new_a, new_b, new_a_is_new, new_b_is_new,&
       error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_c, matrix_a
    CHARACTER, INTENT(IN)                    :: tr_a
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_b
    CHARACTER, INTENT(IN)                    :: tr_b
    TYPE(cp_dbcsr_type), INTENT(OUT)         :: new_a, new_b
    LOGICAL, INTENT(OUT)                     :: new_a_is_new, new_b_is_new
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_match_sizes', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: tr_a_l, tr_b_l
    INTEGER :: a_cols_total, a_rows_total, b_cols_total, b_rows_total, &
      c_cols_total, c_rows_total, handle
    LOGICAL                                  :: atr, btr, new_a_rows, &
                                                new_b_cols, new_b_rows
    TYPE(array_i1d_obj)                      :: new_col_dist, new_col_size, &
                                                new_row_dist, new_row_size
    TYPE(dbcsr_distribution_obj)             :: a_dist, b_dist, c_dist, &
                                                new_dist

!type(dbcsr_imagedistribution_obj) :: b_imgdist
!   ---------------------------------------------------------------------------

    CALL timeset (routineN, handle)
    tr_a_l = tr_a
    CALL uppercase (tr_a_l)
    tr_b_l = tr_b
    CALL uppercase (tr_b_l)
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
       new_b_is_new = .TRUE.
       CALL cp_dbcsr_init (new_b, error=error)
       IF (.NOT. btr) THEN
          IF (new_b_cols) THEN
             CALL match_1_dist (new_col_dist,&
                  dbcsr_distribution_col_dist (b_dist),&
                  dbcsr_distribution_col_dist (c_dist))
          ELSE
             new_col_dist = dbcsr_distribution_col_dist (b_dist)
             CALL array_hold (new_col_dist)
          ENDIF
          new_col_size = cp_dbcsr_col_block_sizes (matrix_c)
          IF (new_b_rows) THEN
             IF (.NOT. atr) THEN
                CALL match_1_dist (new_row_dist,&
                     dbcsr_distribution_row_dist(b_dist),&
                     dbcsr_distribution_col_dist(a_dist))
                new_row_size = cp_dbcsr_col_block_sizes(matrix_a)
             ELSE
                CALL match_1_dist (new_row_dist,&
                     dbcsr_distribution_row_dist(b_dist),&
                     dbcsr_distribution_row_dist(a_dist))
                new_row_size = cp_dbcsr_row_block_sizes(matrix_a)
             ENDIF
          ELSE
             new_row_dist = dbcsr_distribution_row_dist(b_dist)
             CALL array_hold (new_row_dist)
             new_row_size = cp_dbcsr_row_block_sizes (matrix_b)
          ENDIF
       ELSE
          IF (new_b_cols) THEN
             CALL match_1_dist (new_row_dist,&
                  dbcsr_distribution_row_dist (b_dist),&
                  dbcsr_distribution_col_dist (c_dist))
          ELSE
             new_row_dist = dbcsr_distribution_row_dist (b_dist)
             CALL array_hold (new_row_dist)
          ENDIF
          new_row_size = cp_dbcsr_col_block_sizes (matrix_c)
          IF (new_b_rows) THEN
             IF (.not.atr) THEN
                CALL match_1_dist (new_col_dist,&
                     dbcsr_distribution_col_dist (b_dist),&
                     dbcsr_distribution_col_dist (a_dist))
                new_col_size = cp_dbcsr_col_block_sizes (matrix_a)
             ELSE
                CALL match_1_dist (new_col_dist,&
                     dbcsr_distribution_col_dist (b_dist),&
                     dbcsr_distribution_row_dist (a_dist))
                new_col_size = cp_dbcsr_row_block_sizes (matrix_a)
             ENDIF
          ELSE
             new_col_dist = dbcsr_distribution_col_dist(b_dist)
             CALL array_hold (new_col_dist)
             new_col_size = cp_dbcsr_col_block_sizes (matrix_b)
          ENDIF
       ENDIF
       CALL dbcsr_distribution_new (new_dist, dbcsr_distribution_mp (c_dist),&
            row_dist=new_row_dist,&
            col_dist=new_col_dist)
       CALL array_release (new_row_dist)
       CALL array_release (new_col_dist)
       CALL cp_dbcsr_create (new_b, template=matrix_b,&
            dist = new_dist,&
            row_blk_size = new_row_size,&
            col_blk_size = new_col_size,&
            error=error)
       CALL dbcsr_distribution_release (new_dist)
       CALL cp_dbcsr_complete_redistribute (matrix_b, new_b, error=error)
    ELSE
       !new_b_is_new = .FALSE.
       !new_b = matrix_b
       new_b_is_new = .TRUE.
       CALL cp_dbcsr_init (new_b, error=error)
       CALL cp_dbcsr_copy (new_b, matrix_b, shallow_data=.TRUE., error=error)
    ENDIF
    !
    ! Process matrix A
    new_a_rows = a_rows_total .NE. c_rows_total
    IF (new_a_rows) THEN
       CALL cp_dbcsr_init (new_a, error=error)
       IF (atr) THEN
          new_row_dist = dbcsr_distribution_row_dist (a_dist)
          CALL array_hold (new_row_dist)
          new_row_size = cp_dbcsr_row_block_sizes (matrix_a)
          CALL match_1_dist (new_col_dist,&
               dbcsr_distribution_col_dist (a_dist),&
               dbcsr_distribution_row_dist (c_dist))
          new_col_size = cp_dbcsr_row_block_sizes (matrix_c)
       ELSE
          CALL match_1_dist (new_row_dist,&
               dbcsr_distribution_row_dist (a_dist),&
               dbcsr_distribution_row_dist (c_dist))
          new_row_size = cp_dbcsr_row_block_sizes (matrix_c)
          new_col_dist = dbcsr_distribution_col_dist (a_dist)
          CALL array_hold (new_col_dist)
          new_col_size = cp_dbcsr_col_block_sizes (matrix_a)
       ENDIF
       CALL dbcsr_distribution_new (new_dist, dbcsr_distribution_mp (a_dist),&
            row_dist = new_row_dist,&
            col_dist = new_col_dist)
       CALL array_release (new_row_dist)
       CALL array_release (new_col_dist)
       CALL cp_dbcsr_create (new_a, template=matrix_a,&
            dist=new_dist,&
            row_blk_size = new_row_size,&
            col_blk_size = new_col_size,&
            error=error)
       CALL dbcsr_distribution_release (new_dist)
       CALL cp_dbcsr_complete_redistribute (matrix_a, new_a, error=error)
    ELSE
       !new_a_is_new = .FALSE.
       !new_a = matrix_a
       new_a_is_new = .TRUE.
       CALL cp_dbcsr_init (new_a, error=error)
       CALL cp_dbcsr_copy (new_a, matrix_a, shallow_data=.TRUE., error=error)
    ENDIF
    CALL timestop (handle)
  END SUBROUTINE matrix_match_sizes

! *****************************************************************************
!> \brief Matches a 1-D distribution
!>
!> Forms a new distribution new_dist of size template_dist.  As much data as
!> possible is taken from old_dist.
! *****************************************************************************
  SUBROUTINE match_1_dist (new_dist, old_dist, template_dist)
    TYPE(array_i1d_obj), INTENT(out)         :: new_dist
    TYPE(array_i1d_obj), INTENT(in)          :: old_dist, template_dist

    INTEGER                                  :: i, max_bin, new_size, &
                                                old_size, sz
    INTEGER, DIMENSION(:), POINTER           :: new_data, old_data

    old_size = array_size (old_dist)
    new_size = array_size (template_dist)
    sz = MIN(old_size, new_size)
    old_data => array_data (old_dist)
    CALL array_new (new_dist, array_data (template_dist), lb=1)
    new_data => array_data (new_dist)
    new_data(1:sz) = old_data(1:sz)
    max_bin = MAXVAL (old_data)
    IF (max_bin .GT. 0) THEN
       FORALL (i = sz+1 : new_size)
          new_data(i) = MOD (i, max_bin)
       END FORALL
    ELSE
       new_data(sz+1:new_size) = 0
    ENDIF
  END SUBROUTINE match_1_dist

  SUBROUTINE cp_dbcsr_print_block_sum(matrix, unit_nr, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, OPTIONAL                        :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_print_block_sum', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN,handle)

    CALL dbcsr_print_block_sum(matrix%matrix, unit_nr, dbcsr_error)
    CALL timestop(handle)

  END SUBROUTINE cp_dbcsr_print_block_sum

  SUBROUTINE cp_dbcsr_triu(matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_triu', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    !CALL cluster_buffer_flush (matrix%buffers, error=error)
    CALL dbcsr_triu(matrix%matrix, dbcsr_error)
    CALL dbcsr_index_compact(matrix%matrix, dbcsr_error)
    CALL timestop(timing_handle)

  END SUBROUTINE cp_dbcsr_triu

#include "cp_dbcsr_interface_s.f90"
#include "cp_dbcsr_interface_d.f90"
#include "cp_dbcsr_interface_c.f90"
#include "cp_dbcsr_interface_z.f90"

END MODULE cp_dbcsr_interface

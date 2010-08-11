!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
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
  USE array_types,                     ONLY: array_i1d_obj
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_iterator,&
                                             cp_dbcsr_type
  USE dbcsr_block_access,              ONLY: dbcsr_get_block,&
                                             dbcsr_get_block_p,&
                                             dbcsr_put_block,&
                                             dbcsr_reserve_block2d
  USE dbcsr_dist_operations,           ONLY: create_bl_distribution,&
                                             dbcsr_get_stored_coordinates
  USE dbcsr_error_handling,            ONLY: dbcsr_error_type
  USE dbcsr_io,                        ONLY: dbcsr_print
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_offsets, dbcsr_col_block_sizes, dbcsr_distribution, &
       dbcsr_distribution_release, dbcsr_get_data_size, dbcsr_get_data_type, &
       dbcsr_get_info, dbcsr_get_matrix_type, dbcsr_get_num_blocks, &
       dbcsr_get_occupation, dbcsr_init, dbcsr_max_col_size, &
       dbcsr_max_row_size, dbcsr_mp_new, dbcsr_mp_release, dbcsr_name, &
       dbcsr_nblkcols_local, dbcsr_nblkcols_total, dbcsr_nblkrows_local, &
       dbcsr_nblkrows_total, dbcsr_nfullcols_local, dbcsr_nfullcols_total, &
       dbcsr_nfullrows_local, dbcsr_nfullrows_total, dbcsr_release, &
       dbcsr_row_block_offsets, dbcsr_row_block_sizes, &
       dbcsr_uses_special_memory, dbcsr_valid_index
  USE dbcsr_operations,                ONLY: &
       dbcsr_add, dbcsr_add_on_diag, dbcsr_btriu, dbcsr_copy, dbcsr_filter, &
       dbcsr_frobenius_norm, dbcsr_gershgorin_norm, dbcsr_get_block_diag, &
       dbcsr_get_diag, dbcsr_hadamard_product, dbcsr_multiply, dbcsr_norm, &
       dbcsr_scale, dbcsr_scale_by_vector, dbcsr_set, dbcsr_set_diag, &
       dbcsr_sum_replicated, dbcsr_trace
  USE dbcsr_transformations,           ONLY: dbcsr_complete_redistribute,&
                                             dbcsr_distribute,&
                                             dbcsr_new_transposed,&
                                             dbcsr_replicate_all
  USE dbcsr_types,                     ONLY: dbcsr_2d_array_type,&
                                             dbcsr_data_obj,&
                                             dbcsr_distribution_obj,&
                                             dbcsr_mp_obj
  USE dbcsr_util,                      ONLY: dbcsr_checksum,&
                                             dbcsr_verify_matrix
  USE dbcsr_work_operations,           ONLY: dbcsr_create,&
                                             dbcsr_finalize,&
                                             dbcsr_work_create
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop

  !$ USE OMP_LIB
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_interface'

  !
  ! Interface to libdbcsr (contains cp2k timing and error)
  PUBLIC :: cp_dbcsr_set
  PUBLIC :: cp_dbcsr_add
  PUBLIC :: cp_dbcsr_scale
  PUBLIC :: cp_dbcsr_scale_by_vector
  PUBLIC :: cp_dbcsr_hadamard_product
  PUBLIC :: cp_dbcsr_transposed
  PUBLIC :: cp_dbcsr_multiply
  PUBLIC :: cp_dbcsr_copy
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
  PUBLIC :: cp_dbcsr_uses_special_memory
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
  PUBLIC :: cp_dbcsr_complete_redistribute
  PUBLIC :: cp_dbcsr_gershgorin_norm
  PUBLIC :: cp_dbcsr_frobenius_norm
  PUBLIC :: cp_dbcsr_name


  INTERFACE cp_dbcsr_trace
     MODULE PROCEDURE cp_dbcsr_trace_ab_d
     MODULE PROCEDURE cp_dbcsr_trace_ab_s
     MODULE PROCEDURE cp_dbcsr_trace_a_d
     MODULE PROCEDURE cp_dbcsr_trace_a_s
  END INTERFACE

  INTERFACE cp_dbcsr_set
     MODULE PROCEDURE cp_dbcsr_set_d
     MODULE PROCEDURE cp_dbcsr_set_s
     MODULE PROCEDURE cp_dbcsr_set_c
     MODULE PROCEDURE cp_dbcsr_set_z
  END INTERFACE

  INTERFACE cp_dbcsr_add
     MODULE PROCEDURE cp_dbcsr_add_d
     MODULE PROCEDURE cp_dbcsr_add_s
     MODULE PROCEDURE cp_dbcsr_add_c
     MODULE PROCEDURE cp_dbcsr_add_z
  END INTERFACE

  INTERFACE cp_dbcsr_scale
     MODULE PROCEDURE cp_dbcsr_scale_d
     MODULE PROCEDURE cp_dbcsr_scale_s
     MODULE PROCEDURE cp_dbcsr_scale_c
     MODULE PROCEDURE cp_dbcsr_scale_z
     MODULE PROCEDURE cp_dbcsr_scale_d_m
     MODULE PROCEDURE cp_dbcsr_scale_s_m
     MODULE PROCEDURE cp_dbcsr_scale_c_m
     MODULE PROCEDURE cp_dbcsr_scale_z_m
  END INTERFACE

  INTERFACE cp_dbcsr_scale_by_vector
     MODULE PROCEDURE cp_dbcsr_scale_by_vector_d
     MODULE PROCEDURE cp_dbcsr_scale_by_vector_s
     MODULE PROCEDURE cp_dbcsr_scale_by_vector_c
     MODULE PROCEDURE cp_dbcsr_scale_by_vector_z
  END INTERFACE

  INTERFACE cp_dbcsr_multiply
     MODULE PROCEDURE cp_dbcsr_multiply_d
     MODULE PROCEDURE cp_dbcsr_multiply_s
     MODULE PROCEDURE cp_dbcsr_multiply_c
     MODULE PROCEDURE cp_dbcsr_multiply_z
  END INTERFACE

  INTERFACE cp_dbcsr_get_block_p
     MODULE PROCEDURE cp_dbcsr_get_block_p_d, cp_dbcsr_get_block_p_s,&
                      cp_dbcsr_get_block_p_z, cp_dbcsr_get_block_p_c
     MODULE PROCEDURE cp_dbcsr_get_2d_block_p_d, cp_dbcsr_get_2d_block_p_s,&
                      cp_dbcsr_get_2d_block_p_z, cp_dbcsr_get_2d_block_p_c
  END INTERFACE

  INTERFACE cp_dbcsr_get_block
     MODULE PROCEDURE cp_dbcsr_get_block_d, cp_dbcsr_get_block_s,&
                      cp_dbcsr_get_block_z, cp_dbcsr_get_block_c
     MODULE PROCEDURE cp_dbcsr_get_2d_block_d, cp_dbcsr_get_2d_block_s,&
                      cp_dbcsr_get_2d_block_z, cp_dbcsr_get_2d_block_c
  END INTERFACE

  INTERFACE cp_dbcsr_put_block
     MODULE PROCEDURE cp_dbcsr_put_block_d, cp_dbcsr_put_block_s,&
                      cp_dbcsr_put_block_z, cp_dbcsr_put_block_c
     MODULE PROCEDURE cp_dbcsr_put_block2d_d, cp_dbcsr_put_block2d_s,&
                      cp_dbcsr_put_block2d_z, cp_dbcsr_put_block2d_c
  END INTERFACE

  INTERFACE cp_dbcsr_iterator_next_block
     MODULE PROCEDURE cp_iterator_next_block_index
     MODULE PROCEDURE cp_iterator_next_2d_block_d,&
                      cp_iterator_next_2d_block_s,&
                      cp_iterator_next_2d_block_c,&
                      cp_iterator_next_2d_block_z,&
                      cp_iterator_next_1d_block_d,&
                      cp_iterator_next_1d_block_s,&
                      cp_iterator_next_1d_block_c,&
                      cp_iterator_next_1d_block_z
  END INTERFACE

  INTERFACE cp_dbcsr_reserve_block2d
     MODULE PROCEDURE cp_dbcsr_reserve_block2d_d
     MODULE PROCEDURE cp_dbcsr_reserve_block2d_s
     MODULE PROCEDURE cp_dbcsr_reserve_block2d_c
     MODULE PROCEDURE cp_dbcsr_reserve_block2d_z
  END INTERFACE

  PRIVATE

CONTAINS

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

    norm = dbcsr_gershgorin_norm(matrix%matrix)

  END FUNCTION cp_dbcsr_gershgorin_norm

  FUNCTION cp_dbcsr_frobenius_norm(matrix, local) RESULT (norm)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    LOGICAL, INTENT(in), OPTIONAL            :: local
    REAL(KIND=real_8)                        :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_frobenius_norm', &
      routineP = moduleN//':'//routineN

    norm = dbcsr_frobenius_norm(matrix%matrix, local)
  END FUNCTION cp_dbcsr_frobenius_norm


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

  SUBROUTINE cp_dbcsr_reserve_block2d_d(matrix, row, col, block,&
       transposed, existed)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed
    LOGICAL, INTENT(OUT), OPTIONAL           :: existed

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_reserve_block2d_d', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_reserve_block2d(matrix%matrix, row, col, block,&
         transposed, existed)

  END SUBROUTINE cp_dbcsr_reserve_block2d_d

  SUBROUTINE cp_dbcsr_reserve_block2d_s(matrix, row, col, block,&
       transposed, existed)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_4), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed
    LOGICAL, INTENT(OUT), OPTIONAL           :: existed

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_reserve_block2d_s', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_reserve_block2d(matrix%matrix, row, col, block,&
         transposed, existed)

  END SUBROUTINE cp_dbcsr_reserve_block2d_s

  SUBROUTINE cp_dbcsr_reserve_block2d_z(matrix, row, col, block,&
       transposed, existed)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed
    LOGICAL, INTENT(OUT), OPTIONAL           :: existed

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_reserve_block2d_z', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_reserve_block2d(matrix%matrix, row, col, block,&
         transposed, existed)

  END SUBROUTINE cp_dbcsr_reserve_block2d_z

  SUBROUTINE cp_dbcsr_reserve_block2d_c(matrix, row, col, block,&
       transposed, existed)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed
    LOGICAL, INTENT(OUT), OPTIONAL           :: existed

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_reserve_block2d_c', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_reserve_block2d(matrix%matrix, row, col, block,&
         transposed, existed)

  END SUBROUTINE cp_dbcsr_reserve_block2d_c

  PURE FUNCTION cp_dbcsr_get_data_type (matrix) RESULT (data_type)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    CHARACTER                                :: data_type

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

    data_size = dbcsr_get_data_size(matrix%matrix)

  END FUNCTION cp_dbcsr_get_data_size

  FUNCTION cp_dbcsr_row_block_sizes (matrix) RESULT (row_blk_sizes)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: row_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_row_block_sizes', &
      routineP = moduleN//':'//routineN

    row_blk_sizes = dbcsr_row_block_sizes (matrix%matrix)
  END FUNCTION cp_dbcsr_row_block_sizes

  FUNCTION cp_dbcsr_col_block_sizes (matrix) RESULT (col_blk_sizes)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: col_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_col_block_sizes', &
      routineP = moduleN//':'//routineN

    col_blk_sizes = dbcsr_col_block_sizes (matrix%matrix)
  END FUNCTION cp_dbcsr_col_block_sizes

  FUNCTION cp_dbcsr_col_block_offsets (matrix) RESULT (col_blk_offsets)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: col_blk_offsets

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_col_block_offsets', &
      routineP = moduleN//':'//routineN

    col_blk_offsets = dbcsr_col_block_offsets (matrix%matrix)
  END FUNCTION cp_dbcsr_col_block_offsets

  FUNCTION cp_dbcsr_row_block_offsets (matrix) RESULT (row_blk_offsets)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(array_i1d_obj)                      :: row_blk_offsets

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_row_block_offsets', &
      routineP = moduleN//':'//routineN

    row_blk_offsets = dbcsr_row_block_offsets (matrix%matrix)
  END FUNCTION cp_dbcsr_row_block_offsets

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

    occupation = dbcsr_get_occupation (matrix%matrix)

  END FUNCTION cp_dbcsr_get_occupation


  FUNCTION cp_dbcsr_nblkrows_total(matrix) RESULT (nblkrows_total)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nblkrows_total

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nblkrows_total', &
      routineP = moduleN//':'//routineN

    nblkrows_total = dbcsr_nblkrows_total(matrix%matrix)
  END FUNCTION cp_dbcsr_nblkrows_total

  FUNCTION cp_dbcsr_nblkcols_total(matrix) RESULT (nblkcols_total)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nblkcols_total

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nblkcols_total', &
      routineP = moduleN//':'//routineN

    nblkcols_total = dbcsr_nblkcols_total(matrix%matrix)
  END FUNCTION cp_dbcsr_nblkcols_total

  FUNCTION cp_dbcsr_nfullrows_total(matrix) RESULT (nfullrows_total)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nfullrows_total

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nfullrows_total', &
      routineP = moduleN//':'//routineN

    nfullrows_total = dbcsr_nfullrows_total(matrix%matrix)
  END FUNCTION cp_dbcsr_nfullrows_total

  FUNCTION cp_dbcsr_nfullcols_total(matrix) RESULT (nfullcols_total)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nfullcols_total

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nfullcols_total', &
      routineP = moduleN//':'//routineN

    nfullcols_total = dbcsr_nfullcols_total(matrix%matrix)
  END FUNCTION cp_dbcsr_nfullcols_total

  FUNCTION cp_dbcsr_nblkrows_local(matrix) RESULT (nblkrows_local)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nblkrows_local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nblkrows_local', &
      routineP = moduleN//':'//routineN

    nblkrows_local = dbcsr_nblkrows_local(matrix%matrix)
  END FUNCTION cp_dbcsr_nblkrows_local

  FUNCTION cp_dbcsr_nblkcols_local(matrix) RESULT (nblkcols_local)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nblkcols_local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nblkcols_local', &
      routineP = moduleN//':'//routineN

    nblkcols_local = dbcsr_nblkcols_local(matrix%matrix)
  END FUNCTION cp_dbcsr_nblkcols_local

  FUNCTION cp_dbcsr_nfullrows_local(matrix) RESULT (nfullrows_local)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nfullrows_local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nfullrows_local', &
      routineP = moduleN//':'//routineN

    nfullrows_local = dbcsr_nfullrows_local(matrix%matrix)
  END FUNCTION cp_dbcsr_nfullrows_local

  FUNCTION cp_dbcsr_nfullcols_local(matrix) RESULT (nfullcols_local)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: nfullcols_local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_nfullcols_local', &
      routineP = moduleN//':'//routineN

    nfullcols_local = dbcsr_nfullcols_local(matrix%matrix)
  END FUNCTION cp_dbcsr_nfullcols_local

  FUNCTION cp_dbcsr_max_row_size(matrix) RESULT (max_row_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: max_row_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_max_row_size', &
      routineP = moduleN//':'//routineN

    max_row_size = dbcsr_max_row_size(matrix%matrix)
  END FUNCTION cp_dbcsr_max_row_size

  FUNCTION cp_dbcsr_max_col_size(matrix) RESULT (max_col_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER                                  :: max_col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_max_col_size', &
      routineP = moduleN//':'//routineN

    max_col_size = dbcsr_max_col_size(matrix%matrix)
  END FUNCTION cp_dbcsr_max_col_size

  FUNCTION cp_dbcsr_distribution (matrix) RESULT (distribution)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(dbcsr_distribution_obj)             :: distribution

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_distribution', &
      routineP = moduleN//':'//routineN

    distribution = dbcsr_distribution (matrix%matrix)
  END FUNCTION cp_dbcsr_distribution

  PURE FUNCTION cp_dbcsr_uses_special_memory (matrix) RESULT (uses_special)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    LOGICAL                                  :: uses_special

    uses_special = dbcsr_uses_special_memory(matrix%matrix)
  END FUNCTION cp_dbcsr_uses_special_memory

  SUBROUTINE cp_create_bl_distribution (block_distribution,&
       block_size, nelements, nbins)
    TYPE(array_i1d_obj), INTENT(OUT)         :: block_distribution, block_size
    INTEGER, INTENT(IN)                      :: nelements, nbins

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_create_bl_distribution', &
      routineP = moduleN//':'//routineN

    CALL create_bl_distribution (block_distribution,&
       block_size, nelements, nbins)

  END SUBROUTINE cp_create_bl_distribution

  SUBROUTINE cp_dbcsr_distribution_release(dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_dbcsr_distribution_release', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_distribution_release(dist)

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

    CALL dbcsr_mp_release(mp_env)

  END SUBROUTINE cp_dbcsr_mp_release

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

  SUBROUTINE cp_iterator_next_block_index (iterator, row, column, blk,&
       transposed, blk_p, row_size, col_size, row_offset, col_offset)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)   :: iterator
    INTEGER, INTENT(OUT)                     :: row, column, blk
    LOGICAL, INTENT(OUT), OPTIONAL           :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: blk_p, row_size, col_size, &
                                                row_offset, col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iterator_next_block_index', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_iterator_next_block (iterator, row, column, blk,&
       transposed, blk_p, row_size, col_size, row_offset, col_offset)

  END SUBROUTINE cp_iterator_next_block_index

  SUBROUTINE cp_iterator_next_2d_block_d (iterator, row, column,&
       block, transposed,&
       block_number, row_size, col_size, row_offset, col_offset)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)   :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    REAL(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number, row_size, &
                                                col_size, row_offset, &
                                                col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iterator_next_2d_block_d', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_iterator_next_block (iterator, row, column,&
       block, transposed,&
       block_number, row_size, col_size, row_offset, col_offset)

  END SUBROUTINE cp_iterator_next_2d_block_d

  SUBROUTINE cp_iterator_next_1d_block_d (iterator, row, column, block,&
       transposed, block_number, row_size, col_size, row_offset, col_offset)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)   :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    REAL(kind=real_8), DIMENSION(:), POINTER :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number, row_size, &
                                                col_size, row_offset, &
                                                col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iterator_next_1d_block_d', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_iterator_next_block (iterator, row, column, block,&
       transposed, block_number, row_size, col_size, row_offset, col_offset)

  END SUBROUTINE cp_iterator_next_1d_block_d

  SUBROUTINE cp_iterator_next_2d_block_s (iterator, row, column,&
       block, transposed,&
       block_number, row_size, col_size, row_offset, col_offset)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)   :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    REAL(kind=real_4), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number, row_size, &
                                                col_size, row_offset, &
                                                col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iterator_next_2d_block_s', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_iterator_next_block (iterator, row, column,&
       block, transposed,&
       block_number, row_size, col_size, row_offset, col_offset)

  END SUBROUTINE cp_iterator_next_2d_block_s

  SUBROUTINE cp_iterator_next_1d_block_s (iterator, row, column, block,&
       transposed, block_number, row_size, col_size, row_offset, col_offset)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)   :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    REAL(kind=real_4), DIMENSION(:), POINTER :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number, row_size, &
                                                col_size, row_offset, &
                                                col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iterator_next_1d_block_s', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_iterator_next_block (iterator, row, column, block,&
       transposed, block_number, row_size, col_size, row_offset, col_offset)

  END SUBROUTINE cp_iterator_next_1d_block_s

  SUBROUTINE cp_iterator_next_2d_block_z (iterator, row, column,&
       block, transposed,&
       block_number, row_size, col_size, row_offset, col_offset)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)   :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    COMPLEX(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number, row_size, &
                                                col_size, row_offset, &
                                                col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iterator_next_2d_block_z', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_iterator_next_block (iterator, row, column,&
       block, transposed,&
       block_number, row_size, col_size, row_offset, col_offset)

  END SUBROUTINE cp_iterator_next_2d_block_z

  SUBROUTINE cp_iterator_next_1d_block_z (iterator, row, column, block,&
       transposed, block_number, row_size, col_size, row_offset, col_offset)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)   :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    COMPLEX(kind=real_8), DIMENSION(:), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number, row_size, &
                                                col_size, row_offset, &
                                                col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iterator_next_1d_block_z', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_iterator_next_block (iterator, row, column, block,&
       transposed, block_number, row_size, col_size, row_offset, col_offset)

  END SUBROUTINE cp_iterator_next_1d_block_z

  SUBROUTINE cp_iterator_next_2d_block_c (iterator, row, column,&
       block, transposed,&
       block_number, row_size, col_size, row_offset, col_offset)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)   :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    COMPLEX(kind=real_4), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number, row_size, &
                                                col_size, row_offset, &
                                                col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iterator_next_2d_block_c', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_iterator_next_block (iterator, row, column,&
       block, transposed,&
       block_number, row_size, col_size, row_offset, col_offset)

  END SUBROUTINE cp_iterator_next_2d_block_c

  SUBROUTINE cp_iterator_next_1d_block_c (iterator, row, column, block,&
       transposed, block_number, row_size, col_size, row_offset, col_offset)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)   :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    COMPLEX(kind=real_4), DIMENSION(:), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number, row_size, &
                                                col_size, row_offset, &
                                                col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iterator_next_1d_block_c', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_iterator_next_block (iterator, row, column, block,&
       transposed, block_number, row_size, col_size, row_offset, col_offset)

  END SUBROUTINE cp_iterator_next_1d_block_c

  SUBROUTINE cp_dbcsr_iterator_stop (iterator)
    TYPE(cp_dbcsr_iterator), INTENT(INOUT)   :: iterator

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_iterator_stop', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_iterator_stop (iterator)

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

    CALL dbcsr_iterator_start (iterator, matrix%matrix, shared, dynamic,&
       dynamic_byrows, contiguous_pointers, read_only)

  END SUBROUTINE cp_dbcsr_iterator_start

  SUBROUTINE cp_dbcsr_put_block2d_d(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:, :), &
      INTENT(IN)                             :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    REAL(kind=real_8), INTENT(IN), OPTIONAL  :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_put_block2d_d', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_put_block(matrix%matrix, row, col, block, transposed,&
       summation, scale)

  END SUBROUTINE cp_dbcsr_put_block2d_d

  SUBROUTINE cp_dbcsr_put_block2d_s(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_4), DIMENSION(:, :), &
      INTENT(IN)                             :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    REAL(kind=real_4), INTENT(IN), OPTIONAL  :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_put_block2d_s', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_put_block(matrix%matrix, row, col, block, transposed,&
       summation, scale)

  END SUBROUTINE cp_dbcsr_put_block2d_s

  SUBROUTINE cp_dbcsr_put_block2d_z(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:, :), &
      INTENT(IN)                             :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    COMPLEX(kind=real_8), INTENT(IN), &
      OPTIONAL                               :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_put_block2d_z', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_put_block(matrix%matrix, row, col, block, transposed,&
       summation, scale)

  END SUBROUTINE cp_dbcsr_put_block2d_z

  SUBROUTINE cp_dbcsr_put_block2d_c(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:, :), &
      INTENT(IN)                             :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    COMPLEX(kind=real_4), INTENT(IN), &
      OPTIONAL                               :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_put_block2d_c', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_put_block(matrix%matrix, row, col, block, transposed,&
       summation, scale)

  END SUBROUTINE cp_dbcsr_put_block2d_c

  SUBROUTINE cp_dbcsr_put_block_d(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:), &
      INTENT(IN)                             :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    REAL(kind=real_8), INTENT(IN), OPTIONAL  :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_put_block_d', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_put_block(matrix%matrix, row, col, block, transposed,&
       summation, scale)

  END SUBROUTINE cp_dbcsr_put_block_d

  SUBROUTINE cp_dbcsr_put_block_s(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_4), DIMENSION(:), &
      INTENT(IN)                             :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    REAL(kind=real_4), INTENT(IN), OPTIONAL  :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_put_block_s', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_put_block(matrix%matrix, row, col, block, transposed,&
       summation, scale)

  END SUBROUTINE cp_dbcsr_put_block_s

  SUBROUTINE cp_dbcsr_put_block_z(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:), &
      INTENT(IN)                             :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    COMPLEX(kind=real_8), INTENT(IN), &
      OPTIONAL                               :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_put_block_z', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_put_block(matrix%matrix, row, col, block, transposed,&
       summation, scale)

  END SUBROUTINE cp_dbcsr_put_block_z

  SUBROUTINE cp_dbcsr_put_block_c(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:), &
      INTENT(IN)                             :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    COMPLEX(kind=real_4), INTENT(IN), &
      OPTIONAL                               :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_put_block_c', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_put_block(matrix%matrix, row, col, block, transposed,&
       summation, scale)

  END SUBROUTINE cp_dbcsr_put_block_c

  SUBROUTINE cp_dbcsr_get_block_d(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:), &
      INTENT(OUT)                            :: block
    LOGICAL, INTENT(IN)                      :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_block_d', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_get_block(matrix%matrix,row,col,block,tr,found,&
       row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_block_d

  SUBROUTINE cp_dbcsr_get_block_s(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_4), DIMENSION(:), &
      INTENT(OUT)                            :: block
    LOGICAL, INTENT(IN)                      :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_block_s', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_get_block(matrix%matrix,row,col,block,tr,found,&
       row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_block_s

  SUBROUTINE cp_dbcsr_get_block_c(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:), &
      INTENT(OUT)                            :: block
    LOGICAL, INTENT(IN)                      :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_block_c', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_get_block(matrix%matrix,row,col,block,tr,found,&
       row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_block_c

  SUBROUTINE cp_dbcsr_get_block_z(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:), &
      INTENT(OUT)                            :: block
    LOGICAL, INTENT(IN)                      :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_block_z', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_get_block(matrix%matrix,row,col,block,tr,found,&
       row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_block_z

  SUBROUTINE cp_dbcsr_get_2d_block_d(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:, :), &
      INTENT(OUT)                            :: block
    LOGICAL, INTENT(IN)                      :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CALL dbcsr_get_block(matrix%matrix,row,col,block,tr,found,&
         row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_2d_block_d

  SUBROUTINE cp_dbcsr_get_2d_block_s(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_4), DIMENSION(:, :), &
      INTENT(OUT)                            :: block
    LOGICAL, INTENT(IN)                      :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CALL dbcsr_get_block(matrix%matrix,row,col,block,tr,found,&
         row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_2d_block_s

  SUBROUTINE cp_dbcsr_get_2d_block_c(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:, :), &
      INTENT(OUT)                            :: block
    LOGICAL, INTENT(IN)                      :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CALL dbcsr_get_block(matrix%matrix,row,col,block,tr,found,&
         row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_2d_block_c

  SUBROUTINE cp_dbcsr_get_2d_block_z(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:, :), &
      INTENT(OUT)                            :: block
    LOGICAL, INTENT(IN)                      :: tr
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CALL dbcsr_get_block(matrix%matrix,row,col,block,tr,found,&
         row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_2d_block_z

  SUBROUTINE cp_dbcsr_get_2d_block_p_d(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: tr, found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_2d_block_p_d', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_get_block_p(matrix%matrix,row,col,block,tr,found,&
       row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_2d_block_p_d

  SUBROUTINE cp_dbcsr_get_2d_block_p_s(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_4), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: tr, found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_2d_block_p_s', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_get_block_p(matrix%matrix,row,col,block,tr,found,&
       row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_2d_block_p_s

  SUBROUTINE cp_dbcsr_get_2d_block_p_c(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: tr, found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_2d_block_p_c', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_get_block_p(matrix%matrix,row,col,block,tr,found,&
       row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_2d_block_p_c

  SUBROUTINE cp_dbcsr_get_2d_block_p_z(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: tr, found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_2d_block_p_z', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_get_block_p(matrix%matrix,row,col,block,tr,found,&
       row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_2d_block_p_z


  SUBROUTINE cp_dbcsr_get_block_p_d(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_8), DIMENSION(:), POINTER :: block
    LOGICAL, INTENT(OUT)                     :: tr, found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_block_p_d', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_get_block_p(matrix%matrix,row,col,block,tr,found,&
       row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_block_p_d

  SUBROUTINE cp_dbcsr_get_block_p_s(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    REAL(kind=real_4), DIMENSION(:), POINTER :: block
    LOGICAL, INTENT(OUT)                     :: tr, found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_block_p_s', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_get_block_p(matrix%matrix,row,col,block,tr,found,&
       row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_block_p_s

  SUBROUTINE cp_dbcsr_get_block_p_c(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_4), DIMENSION(:), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: tr, found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_block_p_c', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_get_block_p(matrix%matrix,row,col,block,tr,found,&
       row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_block_p_c

  SUBROUTINE cp_dbcsr_get_block_p_z(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    COMPLEX(kind=real_8), DIMENSION(:), &
      POINTER                                :: block
    LOGICAL, INTENT(OUT)                     :: tr, found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_get_block_p_z', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_get_block_p(matrix%matrix,row,col,block,tr,found,&
       row_size, col_size)

  END SUBROUTINE cp_dbcsr_get_block_p_z


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
    CHARACTER, OPTIONAL                      :: matrix_type, data_type

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

  SUBROUTINE cp_dbcsr_replicate_all(matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_replicate_all', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_replicate_all(matrix%matrix, dbcsr_error)
  END SUBROUTINE cp_dbcsr_replicate_all

  SUBROUTINE cp_dbcsr_distribute(matrix, fast, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    LOGICAL, INTENT(in), OPTIONAL            :: fast
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_distribute', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_distribute(matrix%matrix, fast, dbcsr_error)
  END SUBROUTINE cp_dbcsr_distribute

  SUBROUTINE cp_dbcsr_release_p (matrix, error)
    TYPE(cp_dbcsr_type), POINTER             :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_release_p', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(matrix)) THEN
       CALL dbcsr_release (matrix%matrix)
       DEALLOCATE(matrix)
    ENDIF

  END SUBROUTINE cp_dbcsr_release_p

  SUBROUTINE cp_dbcsr_release (matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_release', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_release (matrix%matrix)

  END SUBROUTINE cp_dbcsr_release

  SUBROUTINE cp_dbcsr_init (matrix, error)
    TYPE(cp_dbcsr_type), INTENT(OUT)         :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_init', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_init(matrix%matrix)

  END SUBROUTINE cp_dbcsr_init

  SUBROUTINE cp_dbcsr_init_p (matrix, error)
    TYPE(cp_dbcsr_type), POINTER             :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_init_p', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(matrix)) THEN
       CALL dbcsr_release(matrix%matrix)
       DEALLOCATE(matrix)
    ENDIF

    ALLOCATE(matrix)
    CALL dbcsr_init(matrix%matrix)

  END SUBROUTINE cp_dbcsr_init_p

  SUBROUTINE cp_dbcsr_print(matrix, nodata, matlab_format, variable_name, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: nodata, matlab_format
    CHARACTER(*), INTENT(in), OPTIONAL       :: variable_name
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_print', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_print(matrix%matrix, nodata, matlab_format, variable_name, dbcsr_error)
  END SUBROUTINE cp_dbcsr_print

  SUBROUTINE cp_dbcsr_trace_a_d(matrix_a, trace, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(kind=real_8), INTENT(INOUT)         :: trace
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_trace_a_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_trace(matrix_a%matrix, trace, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_trace_a_d

  SUBROUTINE cp_dbcsr_trace_a_s(matrix_a, trace, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(kind=real_4), INTENT(INOUT)         :: trace
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_trace_a_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_trace(matrix_a%matrix, trace, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_trace_a_s

  SUBROUTINE cp_dbcsr_trace_ab_d(matrix_a, matrix_b, trace, trans_a, trans_b, local_sum, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a, matrix_b
    REAL(kind=real_8), INTENT(INOUT)         :: trace
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: trans_a, trans_b
    LOGICAL, INTENT(IN), OPTIONAL            :: local_sum
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_trace_ab_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_trace(matrix_a%matrix, matrix_b%matrix, trace, trans_a, trans_b, local_sum, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_trace_ab_d

  SUBROUTINE cp_dbcsr_trace_ab_s(matrix_a, matrix_b, trace, trans_a, trans_b, local_sum, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a, matrix_b
    REAL(kind=real_4), INTENT(INOUT)         :: trace
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: trans_a, trans_b
    LOGICAL, INTENT(IN), OPTIONAL            :: local_sum
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_trace_ab_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_trace(matrix_a%matrix, matrix_b%matrix, trace, trans_a, trans_b, local_sum, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_trace_ab_s

  FUNCTION cp_dbcsr_checksum(matrix, local, error) RESULT(checksum)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: local
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    REAL(KIND=dp)                            :: checksum

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_checksum', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    checksum = dbcsr_checksum(matrix%matrix, local, error=dbcsr_error)
  END FUNCTION cp_dbcsr_checksum

  SUBROUTINE cp_dbcsr_sum_replicated (matrix, error)
    TYPE(cp_dbcsr_type), INTENT(inout)       :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_sum_replicated', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_sum_replicated (matrix%matrix, dbcsr_error)
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

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_verify_matrix(m%matrix, verbosity, local, dbcsr_error)
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

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_work_create(matrix%matrix, nblks_guess, sizedata_guess, n,&
         work_mutable, dbcsr_error)
  END SUBROUTINE cp_dbcsr_work_create

  SUBROUTINE cp_dbcsr_create(matrix, name, dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, data_type, special, reuse,&
       mutable_work, replication_type, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    CHARACTER(len=*), INTENT(IN)             :: name
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    CHARACTER, INTENT(IN)                    :: matrix_type
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_size, col_blk_size
    INTEGER, INTENT(IN), OPTIONAL            :: nblks, nze
    CHARACTER, INTENT(IN), OPTIONAL          :: data_type
    LOGICAL, INTENT(IN), OPTIONAL            :: special, reuse, mutable_work
    CHARACTER, INTENT(IN), OPTIONAL          :: replication_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_create', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL dbcsr_create(matrix%matrix, name, dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, data_type, special, reuse,&
       mutable_work, replication_type, dbcsr_error)
  END SUBROUTINE cp_dbcsr_create

  SUBROUTINE cp_dbcsr_filter(matrix, eps, method, use_absolute, filter_diag, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    REAL(dp), INTENT(IN)                     :: eps
    INTEGER, INTENT(IN), OPTIONAL            :: method
    LOGICAL, INTENT(in), OPTIONAL            :: use_absolute, filter_diag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_filter', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_filter(matrix%matrix, eps, method, use_absolute, filter_diag, dbcsr_error)
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
    CALL dbcsr_add_on_diag(matrix%matrix, alpha_scalar, first_row, last_row, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_add_on_diag

  SUBROUTINE cp_dbcsr_copy(matrix_b, matrix_a, name, error, keep_sparsity,&
       shallow_data, keep_imaginary)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_b
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity, shallow_data, &
                                                keep_imaginary

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_copy(matrix_b%matrix, matrix_a%matrix, name, keep_sparsity,&
         shallow_data, keep_imaginary, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_copy

  SUBROUTINE cp_dbcsr_multiply_s(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, left_set, right_set, error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    REAL(real_4), INTENT(IN)                 :: alpha
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a, matrix_b
    REAL(real_4), INTENT(IN)                 :: beta
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_2d_array_type), OPTIONAL, &
      POINTER                                :: left_set, right_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(int_8), INTENT(OUT), OPTIONAL    :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_s', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: shape_a, shape_b, trans_a, &
                                                trans_b
    INTEGER                                  :: timing_handle, &
                                                timing_handle_detail
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    trans_a = transa
    trans_b = transb
    CALL uppercase(trans_a)
    CALL uppercase(trans_b)
    shape_a='R'
    IF(cp_dbcsr_nfullcols_total(matrix_a).EQ.cp_dbcsr_nfullrows_total(matrix_a)) shape_a='S'
    shape_b='R'
    IF(cp_dbcsr_nfullcols_total(matrix_b).EQ.cp_dbcsr_nfullrows_total(matrix_b)) shape_b='S'
    CALL timeset('dbcsr_mult_'//trans_a//shape_a//cp_dbcsr_get_matrix_type(matrix_a)//'_'&
         //trans_b//shape_b//cp_dbcsr_get_matrix_type(matrix_b), timing_handle_detail)
    CALL dbcsr_multiply(transa, transb,&
         alpha, matrix_a%matrix, matrix_b%matrix, beta, matrix_c%matrix,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, left_set=left_set, right_set=right_set, error=dbcsr_error, flop=flop)
    CALL timestop(timing_handle_detail)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_multiply_s

  SUBROUTINE cp_dbcsr_multiply_d(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, left_set, right_set, error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    REAL(real_8), INTENT(IN)                 :: alpha
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a, matrix_b
    REAL(real_8), INTENT(IN)                 :: beta
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_2d_array_type), OPTIONAL, &
      POINTER                                :: left_set, right_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(int_8), INTENT(OUT), OPTIONAL    :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_d', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: shape_a, shape_b, trans_a, &
                                                trans_b
    INTEGER                                  :: timing_handle, &
                                                timing_handle_detail
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)

    trans_a = transa
    trans_b = transb
    CALL uppercase(trans_a)
    CALL uppercase(trans_b)
    shape_a='R'
    IF(cp_dbcsr_nfullcols_total(matrix_a).EQ.cp_dbcsr_nfullrows_total(matrix_a)) shape_a='S'
    shape_b='R'
    IF(cp_dbcsr_nfullcols_total(matrix_b).EQ.cp_dbcsr_nfullrows_total(matrix_b)) shape_b='S'
    CALL timeset('dbcsr_mult_'//trans_a//shape_a//cp_dbcsr_get_matrix_type(matrix_a)//'_'&
         //trans_b//shape_b//cp_dbcsr_get_matrix_type(matrix_b), timing_handle_detail)
    CALL dbcsr_multiply(transa, transb,&
         alpha, matrix_a%matrix, matrix_b%matrix, beta, matrix_c%matrix,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, left_set, right_set, dbcsr_error, flop)
    CALL timestop(timing_handle_detail)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_multiply_d

  SUBROUTINE cp_dbcsr_multiply_c(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, left_set, right_set, error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    COMPLEX(real_4), INTENT(IN)              :: alpha
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a, matrix_b
    COMPLEX(real_4), INTENT(IN)              :: beta
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_2d_array_type), OPTIONAL, &
      POINTER                                :: left_set, right_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(int_8), INTENT(OUT), OPTIONAL    :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_c', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: shape_a, shape_b, trans_a, &
                                                trans_b
    INTEGER                                  :: timing_handle, &
                                                timing_handle_detail
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)

    trans_a = transa
    trans_b = transb
    CALL uppercase(trans_a)
    CALL uppercase(trans_b)
    shape_a='R'
    IF(cp_dbcsr_nfullcols_total(matrix_a).EQ.cp_dbcsr_nfullrows_total(matrix_a)) shape_a='S'
    shape_b='R'
    IF(cp_dbcsr_nfullcols_total(matrix_b).EQ.cp_dbcsr_nfullrows_total(matrix_b)) shape_b='S'
    CALL timeset('dbcsr_mult_'//trans_a//shape_a//cp_dbcsr_get_matrix_type(matrix_a)//'_'&
         //trans_b//shape_b//cp_dbcsr_get_matrix_type(matrix_b), timing_handle_detail)
    CALL dbcsr_multiply(transa, transb,&
         alpha, matrix_a%matrix, matrix_b%matrix, beta, matrix_c%matrix,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, left_set=left_set, right_set=right_set, error=dbcsr_error, flop=flop)
    CALL timestop(timing_handle_detail)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_multiply_c

  SUBROUTINE cp_dbcsr_multiply_z(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, left_set, right_set, error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    COMPLEX(real_8), INTENT(IN)              :: alpha
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a, matrix_b
    COMPLEX(real_8), INTENT(IN)              :: beta
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_2d_array_type), OPTIONAL, &
      POINTER                                :: left_set, right_set
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(int_8), INTENT(OUT), OPTIONAL    :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_multiply_z', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: shape_a, shape_b, trans_a, &
                                                trans_b
    INTEGER                                  :: timing_handle, &
                                                timing_handle_detail
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)

    trans_a = transa
    trans_b = transb
    CALL uppercase(trans_a)
    CALL uppercase(trans_b)
    shape_a='R'
    IF(cp_dbcsr_nfullcols_total(matrix_a).EQ.cp_dbcsr_nfullrows_total(matrix_a)) shape_a='S'
    shape_b='R'
    IF(cp_dbcsr_nfullcols_total(matrix_b).EQ.cp_dbcsr_nfullrows_total(matrix_b)) shape_b='S'
    CALL timeset('dbcsr_mult_'//trans_a//shape_a//cp_dbcsr_get_matrix_type(matrix_a)//'_'&
         //trans_b//shape_b//cp_dbcsr_get_matrix_type(matrix_b), timing_handle_detail)
    CALL dbcsr_multiply(transa, transb,&
         alpha, matrix_a%matrix, matrix_b%matrix, beta, matrix_c%matrix,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, left_set=left_set, right_set=right_set, error=dbcsr_error, flop=flop)
    CALL timestop(timing_handle_detail)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_multiply_z

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
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_new_transposed(transposed%matrix, normal%matrix, shallow_data_copy,&
         transpose_data, transpose_distribution,&
         use_distribution=use_distribution, error=dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_transposed

  SUBROUTINE cp_dbcsr_hadamard_product(matrix_a, matrix_b, matrix_c, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a, matrix_b
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_c
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_hadamard_product', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_hadamard_product(matrix_a%matrix, matrix_b%matrix, matrix_c%matrix, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_hadamard_product


  SUBROUTINE cp_dbcsr_scale_by_vector_d(matrix_a, alpha, side, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(real_8), DIMENSION(:), INTENT(IN), &
      TARGET                                 :: alpha
    CHARACTER(LEN=*), INTENT(IN)             :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_by_vector_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale_by_vector(matrix_a%matrix, alpha, side, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_by_vector_d

  SUBROUTINE cp_dbcsr_scale_by_vector_s(matrix_a, alpha, side, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(real_4), DIMENSION(:), INTENT(IN), &
      TARGET                                 :: alpha
    CHARACTER(LEN=*), INTENT(IN)             :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_by_vector_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale_by_vector(matrix_a%matrix, alpha, side, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_by_vector_s

  SUBROUTINE cp_dbcsr_scale_by_vector_z(matrix_a, alpha, side, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    COMPLEX(real_8), DIMENSION(:), &
      INTENT(IN), TARGET                     :: alpha
    CHARACTER(LEN=*), INTENT(IN)             :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_by_vector_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale_by_vector(matrix_a%matrix, alpha, side, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_by_vector_z

  SUBROUTINE cp_dbcsr_scale_by_vector_c(matrix_a, alpha, side, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    COMPLEX(real_4), DIMENSION(:), &
      INTENT(IN), TARGET                     :: alpha
    CHARACTER(LEN=*), INTENT(IN)             :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_by_vector_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale_by_vector(matrix_a%matrix, alpha, side, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_by_vector_c

  SUBROUTINE cp_dbcsr_scale_d(matrix_a, alpha_scalar, last_column, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(real_8), INTENT(IN)                 :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale(matrix_a%matrix, alpha_scalar, last_column, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_d

  SUBROUTINE cp_dbcsr_scale_s(matrix_a, alpha_scalar, last_column, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(real_4), INTENT(IN)                 :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale(matrix_a%matrix, alpha_scalar, last_column, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_s

  SUBROUTINE cp_dbcsr_scale_z(matrix_a, alpha_scalar, last_column, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    COMPLEX(real_8), INTENT(IN)              :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale(matrix_a%matrix, alpha_scalar, last_column, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_z

  SUBROUTINE cp_dbcsr_scale_c(matrix_a, alpha_scalar, last_column, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    COMPLEX(real_4), INTENT(IN)              :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale(matrix_a%matrix, alpha_scalar, last_column, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_c


  SUBROUTINE cp_dbcsr_scale_d_m(matrix_a, alpha_matrix,&
       last_column, side, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(kind=real_8), DIMENSION(:), &
      INTENT(IN), TARGET                     :: alpha_matrix
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: side
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_d_m', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale(matrix_a%matrix, alpha_matrix, last_column, side, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_d_m

  SUBROUTINE cp_dbcsr_scale_s_m(matrix_a, alpha_matrix,&
       last_column, side, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    REAL(kind=real_4), DIMENSION(:), &
      INTENT(IN), TARGET                     :: alpha_matrix
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: side
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_s_m', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale(matrix_a%matrix, alpha_matrix, last_column, side, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_s_m

  SUBROUTINE cp_dbcsr_scale_z_m(matrix_a, alpha_matrix,&
       last_column, side, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    COMPLEX(kind=real_8), DIMENSION(:), &
      INTENT(IN), TARGET                     :: alpha_matrix
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: side
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_z_m', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale(matrix_a%matrix, alpha_matrix, last_column, side, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_z_m

  SUBROUTINE cp_dbcsr_scale_c_m(matrix_a, alpha_matrix,&
       last_column, side, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    COMPLEX(kind=real_4), DIMENSION(:), &
      INTENT(IN), TARGET                     :: alpha_matrix
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: side
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_scale_c_m', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_scale(matrix_a%matrix, alpha_matrix, last_column, side, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_scale_c_m


  SUBROUTINE cp_dbcsr_set_d(matrix, alpha, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    REAL(real_8), INTENT(IN)                 :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_set_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_set(matrix%matrix, alpha, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_set_d

  SUBROUTINE cp_dbcsr_set_s(matrix, alpha, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    REAL(real_4), INTENT(IN)                 :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_set_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_set(matrix%matrix, alpha, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_set_s

  SUBROUTINE cp_dbcsr_set_z(matrix, alpha, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    COMPLEX(real_8), INTENT(IN)              :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_set_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_set(matrix%matrix, alpha, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_set_z

  SUBROUTINE cp_dbcsr_set_c(matrix, alpha, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    COMPLEX(real_4), INTENT(IN)              :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_set_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_set(matrix%matrix, alpha, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_set_c

  SUBROUTINE cp_dbcsr_add_d(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_b
    REAL(real_8), INTENT(IN)                 :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_add_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_add(matrix_a%matrix, matrix_b%matrix, alpha_scalar, beta_scalar, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_add_d

  SUBROUTINE cp_dbcsr_add_s(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_b
    REAL(real_4), INTENT(IN)                 :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_add_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_add(matrix_a%matrix, matrix_b%matrix, alpha_scalar, beta_scalar, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_add_s

  SUBROUTINE cp_dbcsr_add_z(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_b
    COMPLEX(real_8), INTENT(IN)              :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_add_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_add(matrix_a%matrix, matrix_b%matrix, alpha_scalar, beta_scalar, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_add_z

  SUBROUTINE cp_dbcsr_add_c(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_a
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_b
    COMPLEX(real_4), INTENT(IN)              :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_add_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL dbcsr_add(matrix_a%matrix, matrix_b%matrix, alpha_scalar, beta_scalar, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_add_c

END MODULE cp_dbcsr_interface

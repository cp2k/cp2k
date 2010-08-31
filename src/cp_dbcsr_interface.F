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
  USE array_types,                     ONLY: array_data,&
                                             array_hold,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release,&
                                             array_size
  USE cluster_buffers,                 ONLY: cluster_buffer_flush,&
                                             cluster_buffer_release
  USE cluster_types,                   ONLY: cluster_type
  USE cp_dbcsr_methods,                ONLY: &
       cp_dbcsr_col_block_offsets, cp_dbcsr_col_block_sizes, &
       cp_dbcsr_has_mapping, cp_dbcsr_mapping_activate, &
       cp_dbcsr_mapping_create, cp_dbcsr_mapping_hold, &
       cp_dbcsr_mapping_release, cp_dbcsr_row_block_offsets, &
       cp_dbcsr_row_block_sizes
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_iterator,&
                                             cp_dbcsr_mapping_type,&
                                             cp_dbcsr_type
  USE cp_dbcsr_xlat,                   ONLY: blksizes_blk_to_cluster,&
                                             blksizes_blk_to_cluster4,&
                                             dist_blk_to_cluster,&
                                             get4blksizes,&
                                             set_block_pointer_to_cluster
  USE dbcsr_block_access,              ONLY: dbcsr_get_block,&
                                             dbcsr_get_block_p,&
                                             dbcsr_put_block,&
                                             dbcsr_reserve_block2d
  USE dbcsr_block_buffers,             ONLY: dbcsr_buffers_init,&
                                             dbcsr_buffers_new
  USE dbcsr_block_operations,          ONLY: dbcsr_block_partial_copy,&
                                             dbcsr_data_clear
  USE dbcsr_data_methods,              ONLY: dbcsr_data_clear_pointer,&
                                             dbcsr_data_init,&
                                             dbcsr_data_new,&
                                             dbcsr_data_release,&
                                             dbcsr_data_set_pointer,&
                                             dbcsr_type_1d_to_2d
  USE dbcsr_dist_operations,           ONLY: create_bl_distribution,&
                                             dbcsr_get_stored_coordinates
  USE dbcsr_error_handling,            ONLY: dbcsr_error_type
  USE dbcsr_io,                        ONLY: dbcsr_print
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_methods,                   ONLY: &
       dbcsr_blk_col_offset, dbcsr_blk_column_size, dbcsr_blk_row_offset, &
       dbcsr_blk_row_size, dbcsr_distribution, &
       dbcsr_distribution_add_col_map, dbcsr_distribution_add_row_map, &
       dbcsr_distribution_del_col_map, dbcsr_distribution_del_row_map, &
       dbcsr_distribution_hold, dbcsr_distribution_init, &
       dbcsr_distribution_release, dbcsr_get_data_size, dbcsr_get_data_type, &
       dbcsr_get_info, dbcsr_get_matrix_type, dbcsr_get_num_blocks, &
       dbcsr_get_occupation, dbcsr_has_symmetry, dbcsr_init, &
       dbcsr_max_col_size, dbcsr_max_row_size, dbcsr_mp_new, &
       dbcsr_mp_release, dbcsr_name, dbcsr_nblkcols_local, &
       dbcsr_nblkcols_total, dbcsr_nblkrows_local, dbcsr_nblkrows_total, &
       dbcsr_nfullcols_local, dbcsr_nfullcols_total, dbcsr_nfullrows_local, &
       dbcsr_nfullrows_total, dbcsr_release, dbcsr_uses_special_memory, &
       dbcsr_valid_index
  USE dbcsr_operations,                ONLY: &
       dbcsr_add, dbcsr_add_on_diag, dbcsr_btriu, dbcsr_copy, dbcsr_filter, &
       dbcsr_frobenius_norm, dbcsr_gershgorin_norm, dbcsr_get_block_diag, &
       dbcsr_get_diag, dbcsr_hadamard_product, dbcsr_multiply, dbcsr_norm, &
       dbcsr_scale, dbcsr_scale_by_vector, dbcsr_set, dbcsr_set_diag, &
       dbcsr_sum_replicated, dbcsr_trace
  USE dbcsr_transformations,           ONLY: dbcsr_complete_redistribute,&
                                             dbcsr_distribute,&
                                             dbcsr_new_transposed,&
                                             dbcsr_redistribute,&
                                             dbcsr_replicate_all
  USE dbcsr_types,                     ONLY: dbcsr_2d_array_type,&
                                             dbcsr_block_buffer_obj,&
                                             dbcsr_data_obj,&
                                             dbcsr_distribution_obj,&
                                             dbcsr_iterator,&
                                             dbcsr_mp_obj
  USE dbcsr_util,                      ONLY: convert_sizes_to_offsets,&
                                             dbcsr_checksum,&
                                             dbcsr_verify_matrix,&
                                             swap
  USE dbcsr_work_operations,           ONLY: dbcsr_create,&
                                             dbcsr_finalize,&
                                             dbcsr_work_create
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8,&
                                             int_size,&
                                             real_4,&
                                             real_8
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  !$ USE OMP_LIB
#include "cp_common_uses.h"
#define CHECKUSAGE CALL cp_assert(.FALSE.,cp_fatal_level,cp_unimplemented_error_nr,routineN,"OBS")
#define BUFFERME ! buffering code

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


  INTERFACE cp_dbcsr_create
     MODULE PROCEDURE cp_dbcsr_create_new
     MODULE PROCEDURE cp_dbcsr_create_template
  END INTERFACE

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

  TYPE(cp_dbcsr_mapping_type), POINTER, SAVE :: saved_ao_map => NULL()
  INTEGER, SAVE                     :: map_size = -1
  PUBLIC :: save_map, trash_map

CONTAINS

  SUBROUTINE save_map (clusters, aosize, error)
    TYPE(cluster_type), INTENT(IN)           :: clusters
    INTEGER, INTENT(in)                      :: aosize
    TYPE(cp_error_type), INTENT(INOUT)       :: error

!

    NULLIFY (saved_ao_map)
    CALL cp_dbcsr_mapping_create (saved_ao_map, clusters, error=error)
    map_size = aosize
  END SUBROUTINE save_map

  SUBROUTINE trash_map (error)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

!

    CALL cp_dbcsr_mapping_release (saved_ao_map, error=error)
    NULLIFY (saved_ao_map)
  END SUBROUTINE trash_map

  FUNCTION okay_size (matrix)
    TYPE(cp_dbcsr_type), INTENT(in)          :: matrix
    LOGICAL                                  :: okay_size

!okay_size = (cp_dbcsr_nblkrows_total (matrix) .eq. map_size .and.&
!     cp_dbcsr_nblkcols_total (matrix) .eq. map_size)

    IF (okay_size) &
         WRITE(*,*)"Testing", map_size,&
         cp_dbcsr_nblkrows_total (matrix),&
         cp_dbcsr_nblkcols_total (matrix)
    okay_size = .TRUE.
  END FUNCTION okay_size

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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
    norm = dbcsr_gershgorin_norm(matrix%matrix)

  END FUNCTION cp_dbcsr_gershgorin_norm

  FUNCTION cp_dbcsr_frobenius_norm(matrix, local) RESULT (norm)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    LOGICAL, INTENT(in), OPTIONAL            :: local
    REAL(KIND=real_8)                        :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_frobenius_norm', &
      routineP = moduleN//':'//routineN

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    buffers = matrix%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    CALL cluster_buffer_flush (redist%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error
    TYPE(dbcsr_block_buffer_obj)             :: buffers

    CALL cp_error_init (error)
    buffers = matrix%buffers
    CALL cluster_buffer_flush (buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error
    TYPE(dbcsr_block_buffer_obj)             :: buffers

    CALL cp_error_init (error)
    buffers = matrix%buffers
    CALL cluster_buffer_flush (buffers, error=error)
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

    TYPE(cp_error_type)                      :: error
    TYPE(dbcsr_block_buffer_obj)             :: buffers

    CALL cp_error_init (error)
    buffers = matrix%buffers
    CALL cluster_buffer_flush (buffers, error=error)
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

    TYPE(cp_error_type)                      :: error
    TYPE(dbcsr_block_buffer_obj)             :: buffers

    CALL cp_error_init (error)
    buffers = matrix%buffers
    CALL cluster_buffer_flush (buffers, error=error)
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

    TYPE(cp_error_type)                      :: error
    TYPE(dbcsr_block_buffer_obj)             :: buffers

    CALL cp_error_init (error)
    buffers = matrix%buffers
    CALL cluster_buffer_flush (buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    INTEGER                                  :: low_col, low_row
    LOGICAL                                  :: low_tr
    REAL(kind=real_8), DIMENSION(:), POINTER :: low_block
    TYPE(cp_error_type)                      :: error

    IF (.TRUE.) THEN
       ! Direct
       CALL dbcsr_get_block_p(matrix%matrix,row,col,block,tr,found,&
            row_size, col_size)
    ELSE
       CALL cp_error_init (error)
       CALL cluster_buffer_flush (matrix%buffers, error=error)
       low_row = row
       low_col = col
       CALL dbcsr_get_block_p(matrix%matrix,low_row,low_col,low_block,low_tr,found)
       tr = .FALSE.
       CALL set_block_pointer_to_cluster (matrix, row, col, block, tr,&
            contiguous_pointers = .FALSE., read_only = .FALSE.,&
            low_block=low_block, low_tr=low_tr, error=error)
       CALL fill_sizes (matrix, row, col, row_size, col_size)
    ENDIF
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
    BUFFERME
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
    BUFFERME
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

    TYPE(cp_error_type)                      :: error

    CALL cp_error_init (error)
    CALL cluster_buffer_flush (matrix%buffers, error=error)
    BUFFERME
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

    TYPE(cp_error_type)                      :: error
    TYPE(dbcsr_block_buffer_obj)             :: buffers

    CALL cp_error_init (error)
    buffers = matrix%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    BUFFERME
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

    TYPE(cp_error_type)                      :: error
    TYPE(dbcsr_block_buffer_obj)             :: buffers

    CALL cp_error_init (error)
    buffers = matrix%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    BUFFERME
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

    TYPE(cp_error_type)                      :: error
    TYPE(dbcsr_block_buffer_obj)             :: buffers

    CALL cp_error_init (error)
    buffers = matrix%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    BUFFERME
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

    TYPE(cp_error_type)                      :: error
    TYPE(dbcsr_block_buffer_obj)             :: buffers

    CALL cp_error_init (error)
    buffers = matrix%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    BUFFERME
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

    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

    CALL cluster_buffer_flush (matrix%buffers, error=error)
    CALL dbcsr_replicate_all(matrix%matrix, dbcsr_error)
  END SUBROUTINE cp_dbcsr_replicate_all

  SUBROUTINE cp_dbcsr_distribute(matrix, fast, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    LOGICAL, INTENT(in), OPTIONAL            :: fast
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_distribute', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL cluster_buffer_flush (matrix%buffers, error=error)
    CALL dbcsr_distribute(matrix%matrix, fast, dbcsr_error)
  END SUBROUTINE cp_dbcsr_distribute

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

  SUBROUTINE cp_dbcsr_release (matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_release', &
      routineP = moduleN//':'//routineN

    CALL cluster_buffer_release (matrix%buffers, error=error)
    CALL cp_dbcsr_mapping_release (matrix%row_mapping, error=error)
    CALL cp_dbcsr_mapping_release (matrix%col_mapping, error=error)
    CALL dbcsr_distribution_release (matrix%dbcsr_distribution)
    CALL dbcsr_distribution_release (matrix%cp_distribution)
    !
    CALL array_release (matrix%cp_row_sizes)
    CALL array_release (matrix%cp_col_sizes)
    CALL array_release (matrix%cp_row_offsets)
    CALL array_release (matrix%cp_col_offsets)
    !
    CALL dbcsr_release (matrix%matrix)

  END SUBROUTINE cp_dbcsr_release

  SUBROUTINE cp_dbcsr_init (matrix, error)
    TYPE(cp_dbcsr_type), INTENT(OUT)         :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_init', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_init(matrix%matrix)
    CALL dbcsr_buffers_init (matrix%buffers)
    NULLIFY (matrix%row_mapping, matrix%col_mapping)
    CALL dbcsr_distribution_init (matrix%dbcsr_distribution)
    CALL dbcsr_distribution_init (matrix%cp_distribution)
    CALL array_nullify (matrix%cp_row_sizes)
    CALL array_nullify (matrix%cp_col_sizes)
    CALL array_nullify (matrix%cp_row_offsets)
    CALL array_nullify (matrix%cp_col_offsets)
    matrix%ref_count = 0
  END SUBROUTINE cp_dbcsr_init

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

  SUBROUTINE cp_dbcsr_print(matrix, nodata, matlab_format, variable_name, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: nodata, matlab_format
    CHARACTER(*), INTENT(in), OPTIONAL       :: variable_name
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_print', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    buffers = matrix%buffers
    CALL cluster_buffer_flush (buffers, error=error)
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
    CALL cluster_buffer_flush (matrix_a%buffers, error=error)
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
    CALL cluster_buffer_flush (matrix_a%buffers, error=error)
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
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    buffers = matrix_a%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    buffers = matrix_b%buffers
    CALL cluster_buffer_flush (buffers, error=error)
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
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    buffers = matrix_a%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    buffers = matrix_b%buffers
    CALL cluster_buffer_flush (buffers, error=error)
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

    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    buffers = matrix%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    checksum = dbcsr_checksum(matrix%matrix, local, error=dbcsr_error)
  END FUNCTION cp_dbcsr_checksum

  SUBROUTINE cp_dbcsr_sum_replicated (matrix, error)
    TYPE(cp_dbcsr_type), INTENT(inout)       :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_sum_replicated', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL cluster_buffer_flush (matrix%buffers, error=error)
    CALL dbcsr_sum_replicated (matrix%matrix, dbcsr_error)
  END SUBROUTINE cp_dbcsr_sum_replicated

  SUBROUTINE cp_dbcsr_btriu(matrix_b, matrix_a, error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_b
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_a
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_btriu', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL cluster_buffer_flush (matrix_b%buffers, error=error)
    buffers = matrix_a%buffers
    CALL cluster_buffer_flush (buffers, error=error)
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
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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

  SUBROUTINE cp_dbcsr_create_new(matrix, name, dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, data_type, special, reuse,&
       mutable_work, replication_type, error, row_mapping, col_mapping,round_blk_size_up)
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
    TYPE(cp_dbcsr_mapping_type), OPTIONAL, &
      POINTER                                :: row_mapping, col_mapping
    LOGICAL, INTENT(IN), OPTIONAL            :: round_blk_size_up

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_create_new', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: as, handle, stat
    INTEGER, DIMENSION(:), POINTER           :: offsets
    LOGICAL                                  :: col_mapping_active, round_up, &
                                                row_mapping_active
    TYPE(array_i1d_obj)                      :: low_cbs, low_rbs
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset (routineN, handle)
    IF (PRESENT (round_blk_size_up)) THEN
       round_up = round_blk_size_up
    ELSE
       round_up = .FALSE.
    ENDIF
    IF (ASSOCIATED (matrix%row_mapping)) THEN
       CALL cp_dbcsr_mapping_release (matrix%row_mapping, error)
    ENDIF
    IF (ASSOCIATED (matrix%col_mapping)) THEN
       CALL cp_dbcsr_mapping_release (matrix%col_mapping, error)
    ENDIF
    NULLIFY (matrix%row_mapping, matrix%col_mapping)
    IF (PRESENT (row_mapping)) THEN
       matrix%row_mapping => row_mapping
       CALL cp_dbcsr_mapping_hold (matrix%row_mapping, error=error)
       row_mapping_active = cp_dbcsr_has_mapping (row_mapping)
    ELSE
       CALL cp_dbcsr_mapping_create (matrix%row_mapping, error=error)
       row_mapping_active = .FALSE.
    ENDIF
    IF (PRESENT (col_mapping)) THEN
       matrix%col_mapping => col_mapping
       CALL cp_dbcsr_mapping_hold (matrix%col_mapping, error=error)
       col_mapping_active = cp_dbcsr_has_mapping (col_mapping)
    ELSE       
       CALL cp_dbcsr_mapping_create (matrix%col_mapping, error=error)
       col_mapping_active = .FALSE.
    ENDIF
    !
    IF (.NOT. col_mapping_active .AND. .NOT. row_mapping_active) THEN
       matrix%dbcsr_distribution = dist
       CALL dbcsr_distribution_hold (matrix%dbcsr_distribution)
    ELSE
       ! Create a clustered distribution.
       CALL dist_blk_to_cluster (dist, matrix%row_mapping, matrix%col_mapping,&
            matrix%dbcsr_distribution, error=error)
    ENDIF
    matrix%cp_distribution = dist
    CALL dbcsr_distribution_hold (matrix%cp_distribution)
    !
    IF (row_mapping_active) THEN
       IF (round_up) THEN
          CALL blksizes_blk_to_cluster4 (row_blk_size, matrix%row_mapping,&
               low_rbs, error)
       ELSE
          CALL blksizes_blk_to_cluster (row_blk_size, matrix%row_mapping,&
               low_rbs, error)
       ENDIF
    ELSE
       low_rbs = row_blk_size
       CALL array_hold (low_rbs)
    ENDIF
    IF (col_mapping_active) THEN
       IF (round_up) THEN
          CALL blksizes_blk_to_cluster4 (col_blk_size, matrix%col_mapping,&
               low_cbs, error)
       ELSE
          CALL blksizes_blk_to_cluster (col_blk_size, matrix%col_mapping,&
               low_cbs, error)
       ENDIF
    ELSE
       low_cbs = col_blk_size
       CALL array_hold (low_cbs)
    ENDIF
    CALL dbcsr_create(matrix%matrix, name, matrix%dbcsr_distribution,&
         matrix_type,&
       low_rbs, low_cbs, nblks, nze, data_type, special, reuse,&
       mutable_work, replication_type, dbcsr_error)
    CALL array_release (low_rbs)
    CALL array_release (low_cbs)
    CALL dbcsr_buffers_init (matrix%buffers)
    CALL dbcsr_buffers_new (matrix%buffers, matrix%matrix%m%data_area,&
         dbcsr_error)
    !
    matrix%cp_row_sizes = row_blk_size
    CALL array_hold (matrix%cp_row_sizes)
    matrix%cp_col_sizes = col_blk_size
    CALL array_hold (matrix%cp_col_sizes)
    IF (row_mapping_active) THEN
       as = array_size(row_blk_size)+1
       ALLOCATE (offsets(as), stat=stat)
       IF (stat /= 0) CALL stop_memory (routineN, moduleN, __LINE__,&
            "offsets", int_size*as)
       CALL convert_sizes_to_offsets (array_data (row_blk_size), offsets)
       CALL array_new (matrix%cp_row_offsets, offsets, gift=.TRUE.)
    ELSE
       matrix%cp_row_offsets = matrix%matrix%m%row_blk_offset
       CALL array_hold (matrix%cp_row_offsets)
    ENDIF
    IF (col_mapping_active) THEN
       as = array_size(col_blk_size)+1
       ALLOCATE (offsets(as), stat=stat)
       IF (stat /= 0) CALL stop_memory (routineN, moduleN, __LINE__,&
            "offsets", int_size*as)
       CALL convert_sizes_to_offsets (array_data (col_blk_size), offsets)
       CALL array_new (matrix%cp_col_offsets, offsets, gift=.TRUE.)
    ELSE
       matrix%cp_col_offsets = matrix%matrix%m%col_blk_offset
       CALL array_hold (matrix%cp_col_offsets)
    ENDIF
    !
    matrix%ref_count = 1
    CALL timestop (handle)
  END SUBROUTINE cp_dbcsr_create_new

  SUBROUTINE cp_dbcsr_create_template(matrix, name, template,&
       dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, data_type, special, reuse,&
       mutable_work, replication_type, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: name
    TYPE(cp_dbcsr_type), INTENT(IN)          :: template
    TYPE(dbcsr_distribution_obj), &
      INTENT(IN), OPTIONAL                   :: dist
    CHARACTER, INTENT(IN), OPTIONAL          :: matrix_type
    TYPE(array_i1d_obj), INTENT(IN), &
      OPTIONAL                               :: row_blk_size, col_blk_size
    INTEGER, INTENT(IN), OPTIONAL            :: nblks, nze
    CHARACTER, INTENT(IN), OPTIONAL          :: data_type
    LOGICAL, INTENT(IN), OPTIONAL            :: special, reuse, mutable_work
    CHARACTER, INTENT(IN), OPTIONAL          :: replication_type
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_create_template', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    INTEGER, DIMENSION(:), POINTER           :: offsets
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset (routineN, handle)
    CALL dbcsr_create(matrix%matrix, template%matrix,&
         name=name, dist=dist, matrix_type=matrix_type,&
         row_blk_size=row_blk_size, col_blk_size=col_blk_size,&
         nblks=nblks, nze=nze, data_type=data_type, special=special,&
         reuse=reuse, mutable_work=mutable_work,&
         replication_type=replication_type,&
         error=dbcsr_error)
    CALL dbcsr_buffers_init (matrix%buffers)
    CALL dbcsr_buffers_new (matrix%buffers, matrix%matrix%m%data_area,&
         dbcsr_error)
    !
    IF (ASSOCIATED (matrix%row_mapping)) THEN
       CALL cp_dbcsr_mapping_release (matrix%row_mapping, error)
    ENDIF
    IF (ASSOCIATED (matrix%col_mapping)) THEN
       CALL cp_dbcsr_mapping_release (matrix%col_mapping, error)
    ENDIF
    IF (ASSOCIATED (template%row_mapping)) THEN
       matrix%row_mapping => template%row_mapping
       CALL cp_dbcsr_mapping_hold (matrix%row_mapping, error=error)
    ENDIF
    IF (ASSOCIATED (template%col_mapping)) THEN
       matrix%col_mapping => template%col_mapping
       CALL cp_dbcsr_mapping_hold (matrix%col_mapping, error=error)
    ENDIF
    !
    IF (PRESENT (row_blk_size)) THEN
       matrix%cp_row_sizes = row_blk_size
       ALLOCATE (offsets (cp_dbcsr_nblkrows_total(template)+1), stat=stat)
       IF (stat /= 0) CALL stop_memory (routineN, moduleN, __LINE__, "offsets",&
            int_size*(cp_dbcsr_nblkrows_total(template)+1))
       CALL convert_sizes_to_offsets (array_data (matrix%cp_row_sizes), offsets)
       CALL array_new (matrix%cp_row_offsets, offsets, gift=.TRUE.)
    ELSE
       matrix%cp_row_sizes = cp_dbcsr_row_block_sizes (template)
       matrix%cp_row_offsets = cp_dbcsr_row_block_offsets (template)
       CALL array_hold (matrix%cp_row_offsets)
    ENDIF
    CALL array_hold (matrix%cp_row_sizes)
    !
    IF (PRESENT (col_blk_size)) THEN
       matrix%cp_col_sizes = col_blk_size
       ALLOCATE (offsets (cp_dbcsr_nblkcols_total(template)+1), stat=stat)
       IF (stat /= 0) CALL stop_memory (routineN, moduleN, __LINE__, "offsets",&
            int_size*(cp_dbcsr_nblkcols_total(template)+1))
       CALL convert_sizes_to_offsets (array_data (matrix%cp_col_sizes), offsets)
       CALL array_new (matrix%cp_col_offsets, offsets, gift=.TRUE.)
    ELSE
       matrix%cp_col_sizes = cp_dbcsr_col_block_sizes (template)
       matrix%cp_col_offsets = cp_dbcsr_col_block_offsets (template)
       CALL array_hold (matrix%cp_col_offsets)
    ENDIF
    CALL array_hold (matrix%cp_col_sizes)
    !
    IF (PRESENT (dist)) THEN
       matrix%cp_distribution = dist
    ELSE
       matrix%cp_distribution = cp_dbcsr_distribution (template)
    ENDIF
    CALL dbcsr_distribution_hold (matrix%cp_distribution)
    matrix%dbcsr_distribution = dbcsr_distribution (matrix%matrix)
    CALL dbcsr_distribution_hold (matrix%dbcsr_distribution)
    !
    matrix%ref_count = 1
    CALL timestop (handle)
  END SUBROUTINE cp_dbcsr_create_template


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
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    buffers = matrix%buffers
    CALL cluster_buffer_flush (buffers, error=error)
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
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    buffers = matrix%buffers
    CALL cluster_buffer_flush (buffers, error=error)
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
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    !call cp_assert (matrix_b%ref_count .gt. 0,&
    !     cp_warning_level, cp_caller_error, routineN,&
    !     "Matrix not created.", error=error)
    CALL cluster_buffer_flush (matrix_b%buffers, error)
    IF (matrix_b%ref_count .EQ. 0) THEN
       CALL cp_dbcsr_create (matrix_b, template=matrix_a,&
            error=error)
    ENDIF
    buffers = matrix_a%buffers
    CALL cluster_buffer_flush (buffers, error)
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
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    buffers = matrix_a%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    buffers = matrix_b%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    CALL cluster_buffer_flush (matrix_c%buffers, error=error)
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
    LOGICAL, PARAMETER                       :: prnt = .FALSE., &
                                                verify = .FALSE.

    CHARACTER(LEN=1)                         :: shape_a, shape_b, trans_a, &
                                                trans_b
    INTEGER                                  :: timing_handle, &
                                                timing_handle_clustering, &
                                                timing_handle_detail
    REAL(kind=dp)                            :: cs_b, cs_c
    TYPE(cp_dbcsr_type)                      :: clustered_a, clustered_b, &
                                                clustered_c
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    buffers = matrix_a%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    buffers = matrix_b%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    CALL cluster_buffer_flush (matrix_c%buffers, error=error)
    !
    IF (.FALSE.) THEN
       IF (prnt) &
            WRITE(*,*)routineN//" Multiply======================", alpha, beta
       CALL cp_dbcsr_init (clustered_a, error)
       CALL cp_dbcsr_docluster (matrix_a, clustered_a, error=error)
       IF (prnt) THEN
          CALL cp_dbcsr_print (matrix_a, matlab_format=.TRUE.,&
               variable_name="mao", error=error)
          CALL cp_dbcsr_print (clustered_a, matlab_format=.TRUE.,&
               variable_name="man", error=error)
       ENDIF
       CALL cp_dbcsr_init (clustered_b, error)
       CALL cp_dbcsr_docluster (matrix_b, clustered_b, error=error)
       IF (prnt) THEN
          CALL cp_dbcsr_print (matrix_b, matlab_format=.TRUE.,&
               variable_name="mbo", error=error)
          CALL cp_dbcsr_print (clustered_b, matlab_format=.TRUE.,&
               variable_name="mbn", error=error)
       ENDIF
       CALL cp_dbcsr_init (clustered_c, error)
       CALL cp_dbcsr_docluster (matrix_c, clustered_c, error=error)
       IF (prnt) THEN
          CALL cp_dbcsr_print (matrix_c, matlab_format=.TRUE.,&
               variable_name="mco", error=error)
          CALL cp_dbcsr_print (clustered_c, matlab_format=.TRUE.,&
               variable_name="mcn", error=error)
       ENDIF
       CALL timeset ("multiply_clustering", timing_handle_clustering)
       CALL dbcsr_multiply (transa, transb,&
            alpha, clustered_a%matrix, clustered_b%matrix, beta, clustered_c%matrix,&
            first_row, last_row, first_column, last_column, first_k, last_k,&
            retain_sparsity, left_set, right_set, dbcsr_error, flop)
       IF (verify) cs_c = cp_dbcsr_checksum (clustered_c, error=error)
       IF (prnt) THEN
          CALL cp_dbcsr_print (clustered_c, matlab_format=.TRUE.,&
               variable_name="mpn", error=error)
       ENDIF
       CALL cp_dbcsr_release (clustered_a, error)
       CALL cp_dbcsr_release (clustered_b, error)
       CALL cp_dbcsr_release (clustered_c, error)
       CALL timestop (timing_handle_clustering)
       !!
       IF (prnt) &
            WRITE(*,*)routineN//" Multiply2======================", alpha, beta
       CALL cp_dbcsr_init (clustered_a, error)
       CALL cp_dbcsr_docluster (matrix_a, clustered_a, round_blk_size_up=.TRUE., error=error)
       IF (prnt) THEN
          CALL cp_dbcsr_print (matrix_a, matlab_format=.TRUE.,&
               variable_name="mao", error=error)
          CALL cp_dbcsr_print (clustered_a, matlab_format=.TRUE.,&
               variable_name="man", error=error)
       ENDIF
       CALL cp_dbcsr_init (clustered_b, error)
       CALL cp_dbcsr_docluster (matrix_b, clustered_b, round_blk_size_up=.TRUE., error=error)
       IF (prnt) THEN
          CALL cp_dbcsr_print (matrix_b, matlab_format=.TRUE.,&
               variable_name="mbo", error=error)
          CALL cp_dbcsr_print (clustered_b, matlab_format=.TRUE.,&
               variable_name="mbn", error=error)
       ENDIF
       CALL cp_dbcsr_init (clustered_c, error)
       CALL cp_dbcsr_docluster (matrix_c, clustered_c, round_blk_size_up=.TRUE., error=error)
       IF (prnt) THEN
          CALL cp_dbcsr_print (matrix_c, matlab_format=.TRUE.,&
               variable_name="mco", error=error)
          CALL cp_dbcsr_print (clustered_c, matlab_format=.TRUE.,&
               variable_name="mcn", error=error)
       ENDIF
       CALL timeset ("multiply_clustering4", timing_handle_clustering)
       CALL dbcsr_multiply (transa, transb,&
            alpha, clustered_a%matrix, clustered_b%matrix, beta, clustered_c%matrix,&
            first_row, last_row, first_column, last_column, first_k, last_k,&
            retain_sparsity, left_set, right_set, dbcsr_error, flop)
       IF (verify) cs_c = cp_dbcsr_checksum (clustered_c, error=error)
       IF (prnt) THEN
          CALL cp_dbcsr_print (clustered_c, matlab_format=.TRUE.,&
               variable_name="mpn", error=error)
       ENDIF
       CALL cp_dbcsr_release (clustered_a, error)
       CALL cp_dbcsr_release (clustered_b, error)
       CALL cp_dbcsr_release (clustered_c, error)
       CALL timestop (timing_handle_clustering)
    ENDIF
    !!
    CALL timeset(routineN, timing_handle)
    !
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
    IF (prnt) THEN
       CALL cp_dbcsr_print (matrix_c, matlab_format=.TRUE.,&
            variable_name="mpo", error=error)
    ENDIF
    IF (verify) cs_b = cp_dbcsr_checksum (matrix_c, error=error)
    CALL timestop(timing_handle_detail)
    CALL timestop(timing_handle)
    IF (verify) THEN
       WRITE(*,'(A,4(1X,E9.3))')routineN//" checksums", cs_c, cs_b,&
            cs_c-cs_b, ABS(cs_c-cs_b)/cs_b
       WRITE(*,*)routineN//" multiply type",&
            trans_a//shape_a//cp_dbcsr_get_matrix_type(matrix_a)//'_'&
            //trans_b//shape_b//cp_dbcsr_get_matrix_type(matrix_b)
       
       IF (ABS(cs_c-cs_b) .GT. 0.00001) STOP "Bad multiply"
    ENDIF
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
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    buffers = matrix_a%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    buffers = matrix_b%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    CALL cluster_buffer_flush (matrix_c%buffers, error=error)
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
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    buffers = matrix_a%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    buffers = matrix_b%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    CALL cluster_buffer_flush (matrix_c%buffers, error=error)
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
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    buffers = normal%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    CALL cluster_buffer_flush (transposed%buffers, error=error)
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
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    buffers = matrix_a%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    buffers = matrix_b%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    CALL cluster_buffer_flush (matrix_c%buffers, error=error)
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
    CALL cluster_buffer_flush (matrix_a%buffers, error=error)
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
    CALL cluster_buffer_flush (matrix_a%buffers, error=error)
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
    CALL cluster_buffer_flush (matrix_a%buffers, error=error)
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
    CALL cluster_buffer_flush (matrix_a%buffers, error=error)
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
    CALL cluster_buffer_flush (matrix_a%buffers, error=error)
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
    CALL cluster_buffer_flush (matrix_a%buffers, error=error)
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
    CALL cluster_buffer_flush (matrix_a%buffers, error=error)
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
    CALL cluster_buffer_flush (matrix_a%buffers, error=error)
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
    TYPE(cp_error_type)                      :: cp_error
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL cp_error_init (cp_error)
    dbcsr_error = error
    CALL cluster_buffer_flush (matrix_a%buffers, error=cp_error)
    CALL dbcsr_scale(matrix_a%matrix, alpha_matrix, last_column, side, dbcsr_error)
    error = dbcsr_error
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
    TYPE(cp_error_type)                      :: cp_error
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL cp_error_init (cp_error)
    dbcsr_error = error
    CALL cluster_buffer_flush (matrix_a%buffers, error=cp_error)
    CALL dbcsr_scale(matrix_a%matrix, alpha_matrix, last_column, side, dbcsr_error)
    error = dbcsr_error
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
    TYPE(cp_error_type)                      :: cp_error
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    dbcsr_error = error
    CALL cp_error_init (cp_error)
    CALL cluster_buffer_flush (matrix_a%buffers, error=cp_error)
    CALL dbcsr_scale(matrix_a%matrix, alpha_matrix, last_column, side, dbcsr_error)
    error = dbcsr_error
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
    TYPE(cp_error_type)                      :: cp_error
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL cp_error_init (cp_error)
    dbcsr_error = error
    CALL cluster_buffer_flush (matrix_a%buffers, error=cp_error)
    CALL dbcsr_scale(matrix_a%matrix, alpha_matrix, last_column, side, dbcsr_error)
    error = dbcsr_error
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
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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
    CALL cluster_buffer_flush (matrix%buffers, error=error)
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
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL cluster_buffer_flush (matrix_a%buffers, error=error)
    buffers = matrix_b%buffers
    CALL cluster_buffer_flush (buffers, error=error)
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
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL cluster_buffer_flush (matrix_a%buffers, error=error)
    buffers = matrix_b%buffers
    CALL cluster_buffer_flush (buffers, error=error)
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
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL cluster_buffer_flush (matrix_a%buffers, error=error)
    buffers = matrix_b%buffers
    CALL cluster_buffer_flush (buffers, error=error)
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
    TYPE(dbcsr_block_buffer_obj)             :: buffers
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN, timing_handle)
    CALL cluster_buffer_flush (matrix_a%buffers, error=error)
    buffers = matrix_b%buffers
    CALL cluster_buffer_flush (buffers, error=error)
    CALL dbcsr_add(matrix_a%matrix, matrix_b%matrix, alpha_scalar, beta_scalar, dbcsr_error)
    CALL timestop(timing_handle)
  END SUBROUTINE cp_dbcsr_add_c


  ! Buffer routines

  SUBROUTINE cp_dbcsr_flush_buffers (matrix, error)
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CALL cluster_buffer_flush (matrix%buffers, error)
  END SUBROUTINE cp_dbcsr_flush_buffers

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
!> \brief Converts a non-clustered matrix to a clustered.
! *****************************************************************************
  SUBROUTINE cp_dbcsr_docluster (matrix, clustered_matrix, round_blk_size_up, error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix
    TYPE(cp_dbcsr_type), INTENT(OUT)         :: clustered_matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: round_blk_size_up
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dbcsr_docluster', &
      routineP = moduleN//':'//routineN

    INTEGER :: block_col, block_col_size, block_row, block_row_size, &
      cluster_col, cluster_col_size, cluster_row, cluster_row_size, &
      col_offset, icopy, max_blk_size, ncopies, row_offset, timing_handle
    LOGICAL                                  :: block_tr, cluster_found, &
                                                cluster_tr, mirrored
    TYPE(array_i1d_obj)                      :: col_map, row_map
    TYPE(cp_dbcsr_mapping_type), POINTER     :: cl_col_mapping, cl_row_mapping
    TYPE(cp_dbcsr_type)                      :: predistributed
    TYPE(dbcsr_data_obj)                     :: block_data, buffer_data, &
                                                cluster_data
    TYPE(dbcsr_distribution_obj)             :: block_distribution
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_iterator)                     :: block_iter

!   ---------------------------------------------------------------------------

    CALL timeset (routineN, timing_handle)
    !
    NULLIFY (cl_row_mapping, cl_col_mapping)
    !CALL cp_dbcsr_mapping_activate (matrix%row_mapping, cl_row_mapping,&
    !     activation=.TRUE., error=error)
    !CALL cp_dbcsr_mapping_activate (matrix%col_mapping, cl_col_mapping,&
    !     activation=.TRUE., error=error)
    !write(*,*)"row, col", cp_dbcsr_nblkrows_total (matrix) .eq. map_size,&
    !     cp_dbcsr_nblkcols_total (matrix) .eq. map_size
    IF (cp_dbcsr_nblkrows_total (matrix) .EQ. map_size) THEN
       CALL cp_dbcsr_mapping_activate (saved_ao_map, cl_row_mapping,&
            activation=.TRUE., error=error)
    ELSE
       CALL cp_dbcsr_mapping_create (cl_row_mapping, error=error)
    ENDIF
    IF (cp_dbcsr_nblkcols_total (matrix) .EQ. map_size) THEN
       CALL cp_dbcsr_mapping_activate (saved_ao_map, cl_col_mapping,&
            activation=.TRUE., error=error)
    ELSE
       CALL cp_dbcsr_mapping_create (cl_col_mapping, error=error)
    ENDIF
    !
    block_distribution = dbcsr_distribution (matrix%matrix)
    IF (cp_dbcsr_has_mapping (cl_row_mapping)) THEN
       CALL array_new (row_map, cl_row_mapping%clusters%particle_to_cluster,&
            lb=1)
       CALL dbcsr_distribution_add_row_map (block_distribution, row_map)
       CALL array_release (row_map)
    ENDIF
    IF (cp_dbcsr_has_mapping (cl_col_mapping)) THEN
       CALL array_new (col_map, cl_col_mapping%clusters%particle_to_cluster,&
            lb=1)
       CALL dbcsr_distribution_add_col_map (block_distribution, col_map)
       CALL array_release (col_map)
    ENDIF
    !
    CALL cp_dbcsr_init (predistributed, error=error)
    CALL cp_dbcsr_create (predistributed, dbcsr_name (matrix%matrix),&
         block_distribution,&
         dbcsr_get_matrix_type (matrix%matrix),&
         cp_dbcsr_row_block_sizes (matrix),&
         cp_dbcsr_col_block_sizes (matrix),&
         data_type = cp_dbcsr_get_data_type (matrix),&
         error=error)
    !call cp_dbcsr_redistribute (matrix, predistributed, error=error)
    CALL dbcsr_redistribute (matrix%matrix, predistributed%matrix, error=dbcsr_error)
    IF (cp_dbcsr_has_mapping (cl_row_mapping)) THEN
       CALL dbcsr_distribution_del_row_map (block_distribution)
    ENDIF
    IF (cp_dbcsr_has_mapping (cl_col_mapping)) THEN
       CALL dbcsr_distribution_del_col_map (block_distribution)
    ENDIF
    !
    !
    !CALL dbcsr_distribution_init (clustered_distribution)
    !CALL blksizes_blk_to_cluster (matrix%cp_row_sizes,&
    !     matrix%row_mapping, clustered_row_blk_sizes, error=error)
    !CALL blksizes_blk_to_cluster (matrix%cp_col_sizes,&
    !     matrix%col_mapping, clustered_col_blk_sizes, error=error)
    !CALL dist_blk_to_cluster (matrix%cp_distribution, matrix%row_mapping,&
    !     matrix%col_mapping, clustered_distribution, error=error)
    !!
    !CALL cp_dbcsr_init (clustered_matrix, error=error)
    !CALL cp_dbcsr_create (clustered_matrix, dbcsr_name (matrix%matrix),&
    !     clustered_distribution, dbcsr_get_matrix_type (matrix%matrix),&
    !     clustered_row_blk_sizes, clustered_col_blk_sizes,&
    !     nze=dbcsr_get_data_size (matrix%matrix),&
    !     data_type=dbcsr_get_data_type(matrix%matrix),&
    !     row_mapping = cl_row_mapping, col_mapping = cl_col_mapping,&
    !     error=error)
    !CALL cp_dbcsr_mapping_release (cl_row_mapping, error)
    !CALL cp_dbcsr_mapping_release (cl_col_mapping, error)
    !CALL array_release (clustered_row_blk_sizes)
    !CALL array_release (clustered_col_blk_sizes)
    !CALL dbcsr_distribution_release (clustered_distribution)
    CALL cp_dbcsr_create (clustered_matrix, dbcsr_name (matrix%matrix),&
         matrix%cp_distribution, dbcsr_get_matrix_type (matrix%matrix),&
         matrix%cp_row_sizes, matrix%cp_col_sizes,&
         nze=dbcsr_get_data_size (matrix%matrix),&
         data_type=dbcsr_get_data_type(matrix%matrix),&
         row_mapping = cl_row_mapping, col_mapping = cl_col_mapping,&
         round_blk_size_up=round_blk_size_up,&
         error=error)
    CALL cp_dbcsr_finalize (clustered_matrix, error=error)
    CALL dbcsr_work_create (clustered_matrix%matrix, work_mutable=.TRUE., error=dbcsr_error)
    !
    CALL dbcsr_data_init (block_data)
    CALL dbcsr_data_init (cluster_data)
    CALL dbcsr_data_init (buffer_data)
    CALL dbcsr_data_new (block_data,&
         dbcsr_type_1d_to_2d (dbcsr_get_data_type (matrix%matrix)))
    CALL dbcsr_data_new (cluster_data,&
         dbcsr_type_1d_to_2d (dbcsr_get_data_type (matrix%matrix)))
    max_blk_size = MAX (dbcsr_max_row_size (clustered_matrix%matrix),&
            dbcsr_max_col_size (clustered_matrix%matrix))
    CALL dbcsr_data_new (buffer_data,&
         dbcsr_type_1d_to_2d (dbcsr_get_data_type (matrix%matrix)),&
         max_blk_size, max_blk_size)
    !!@@@
    !write(*,*)routineN//" high rbs", array_data(dbcsr_row_block_sizes (matrix%matrix))
    !write(*,*)routineN//" high cbs", array_data(dbcsr_col_block_sizes (matrix%matrix))
    !write(*,*)routineN//" low rbs", array_data(dbcsr_row_block_sizes (clustered_matrix%matrix))
    !write(*,*)routineN//" low cbs", array_data(dbcsr_col_block_sizes (clustered_matrix%matrix))
    CALL dbcsr_iterator_start (block_iter, predistributed%matrix, read_only=.TRUE.,&
         contiguous_pointers=.FALSE.)
    block_iteration: DO WHILE (dbcsr_iterator_blocks_left (block_iter))
       CALL dbcsr_iterator_next_block (block_iter, &
            row=block_row, column=block_col, transposed=block_tr,&
            block=block_data,&
            row_size=block_row_size, col_size=block_col_size)
       ! Warnings about aliasing here are OK.
       CALL get4blksizes (clustered_matrix,&
            high_row=block_row, high_col=block_col,&
            low_r=cluster_row, low_c=cluster_col,&
            low_rs=cluster_row_size, low_cs=cluster_col_size,&
            row_off=row_offset,  col_off=col_offset, error=error)
       !
       mirrored = .FALSE.
       ncopies = 1
       IF (dbcsr_has_symmetry (clustered_matrix%matrix)) THEN
          IF (cluster_row .GT. cluster_col) THEN
             mirrored = .TRUE.
             CALL swap (cluster_row, cluster_col)
          ELSEIF (cluster_row .EQ. cluster_col .AND. block_row .NE. block_col) THEN
             ncopies = 2
          ENDIF
       ENDIF
       !
       cluster_tr = .FALSE.
       CALL dbcsr_get_block_p (clustered_matrix%matrix,&
            cluster_row, cluster_col, cluster_data, cluster_tr, cluster_found,&
            row_size=cluster_row_size, col_size=cluster_col_size)
       IF (.NOT. cluster_found) THEN
          CALL dbcsr_data_set_pointer (cluster_data,&
               rsize=cluster_row_size, csize=cluster_col_size,&
               pointee=buffer_data)
          CALL dbcsr_data_clear (cluster_data)
          !cluster_tr = mirrored
       ENDIF
       IF (mirrored) THEN
          CALL swap (cluster_row_size, cluster_col_size)
       ENDIF
       cluster_tr = cluster_tr .NEQV. mirrored
       DO icopy = 1, ncopies
          IF (icopy .EQ. 2) THEN
             CALL swap (block_row_size, block_col_size)
             CALL swap (row_offset, col_offset)
             block_tr = .NOT. block_tr
          ENDIF
          !!@@@
          !!write(*,*)routineN//" blk row,col", block_row, block_col
          !!write(*,*)routineN//" cls row,col",&
          !!     cluster_row, cluster_col, cluster_found, row_offset, col_offset,&
          !!     cluster_row_size, cluster_col_size, cluster_tr, mirrored
          !!@@@
          !write(*,'(A,4(1X,I5),2(1X,L1))')"Block row, col -> cluster row, col, m, f",&
          !     block_row, block_col, cluster_row, cluster_col, mirrored, cluster_found
          !write(*,'(A,6(1X,I5),1X,L1)')"Block, cluster row, col size, offsets, tr",&
          !     block_row_size, block_col_size,&
          !     cluster_row_size, cluster_col_size,&
          !     row_offset, col_offset,&
          !     cluster_tr
          !call dbcsr_print2dmat (block_data%d%r2_dp)
          CALL dbcsr_block_partial_copy(&
               dst=cluster_data, dst_tr=cluster_tr,&
               dst_rs=cluster_row_size, dst_cs=cluster_col_size,&
               dst_r_lb=row_offset, dst_c_lb=col_offset,&
               src=block_data, src_rs=block_row_size, src_cs=block_col_size,&
               src_tr=block_tr, src_r_lb=1, src_c_lb=1,&
               nrow=block_row_size, ncol=block_col_size)
          !call dbcsr_print2dmat (cluster_data%d%r2_dp)
       ENDDO
       IF (.NOT. cluster_found) THEN
          cluster_tr = .FALSE.
          CALL dbcsr_put_block (clustered_matrix%matrix,&
               cluster_row, cluster_col, cluster_data, cluster_tr)
       ENDIF
    ENDDO block_iteration
    CALL dbcsr_iterator_stop (block_iter)
    CALL cp_dbcsr_release (predistributed, error=error)
    CALL dbcsr_data_clear_pointer (block_data)
    CALL dbcsr_data_clear_pointer (cluster_data)
    CALL dbcsr_data_release (block_data)
    CALL dbcsr_data_release (cluster_data)
    CALL dbcsr_data_release (buffer_data)
    CALL dbcsr_finalize (clustered_matrix%matrix, error=dbcsr_error)
    ! Can't use this because blocks are not in the same positions.
    !!!call cp_dbcsr_complete_redistribute (matrix, clustered_matrix, error=error)
    !@@@
    !!write(*,*)routineN//" sizesO",&
    !!     dbcsr_get_data_size_referenced (matrix%matrix),&
    !!     dbcsr_get_data_size(matrix%matrix), cp_dbcsr_get_occupation (matrix)
    !!write(*,*)routineN//" sizesN",&
    !!     dbcsr_get_data_size_referenced (clustered_matrix%matrix),&
    !!     dbcsr_get_data_size(clustered_matrix%matrix),&
    !!     cp_dbcsr_get_occupation (clustered_matrix)
    !!!a = REAL(dbcsr_get_data_size_referenced (matrix%matrix), KIND=dp)
    !!!b = REAL(dbcsr_get_data_size_referenced (clustered_matrix%matrix), KIND=dp)
    !!a = cp_dbcsr_get_occupation (matrix)
    !!b = cp_dbcsr_get_occupation (clustered_matrix)
    !!if (a .gt. 0.0001) THEN
    !!   size_inflation = (b-a)/a
    !!else
    !!   size_inflation = 0.0_dp
    !!endif
    !!write(*,'(1X,A,3(1X,F9.3))')routineN//" Inflation", size_inflation, a, b
    !call cp_dbcsr_print (matrix, nodata=.true., error=error)
    !call cp_dbcsr_print (clustered_matrix, nodata=.true., error=error)
    CALL timestop (timing_handle)
  END SUBROUTINE cp_dbcsr_docluster


END MODULE cp_dbcsr_interface

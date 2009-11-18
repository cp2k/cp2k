!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR sparse matrix utility routines
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.9
!>
!> <b>Modification history:</b>
!> - Created 2009-05-12
! *****************************************************************************
MODULE dbcsr_util

  USE array_types,                     ONLY: array_data,&
                                             array_hold,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release,&
                                             array_size
  USE btree_I8_k_dp2d_v,               ONLY: btree_2d_data_d => dp2d,&
                                             btree_destroy_d => btree_delete,&
                                             btree_size_d => btree_get_entries
  USE cp_array_i_utils,                ONLY: cp_1d_i_p_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release,&
                                             cp_blacs_env_type,&
                                             cp_blacs_env_write,&
                                             get_blacs_info
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE dbcsr_methods,                   ONLY: &
       dbcsr_data_copyall, dbcsr_data_hold, dbcsr_data_init, dbcsr_data_new, &
       dbcsr_data_release, dbcsr_data_set_pointer, dbcsr_destroy, &
       dbcsr_distribution_col_dist, dbcsr_distribution_hold, &
       dbcsr_distribution_init, dbcsr_distribution_local_cols, &
       dbcsr_distribution_local_rows, dbcsr_distribution_mp, &
       dbcsr_distribution_ncols, dbcsr_distribution_new, &
       dbcsr_distribution_nlocal_cols, dbcsr_distribution_nlocal_rows, &
       dbcsr_distribution_nrows, dbcsr_distribution_processor, &
       dbcsr_distribution_release, dbcsr_distribution_row_dist, &
       dbcsr_get_data, dbcsr_get_data_p, dbcsr_get_data_size, &
       dbcsr_get_data_type, dbcsr_get_info, dbcsr_get_num_blocks, dbcsr_init, &
       dbcsr_is_initialized, dbcsr_iterator_blocks_left, &
       dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, &
       dbcsr_modify_lock, dbcsr_modify_unlock, dbcsr_mp_group, dbcsr_mp_new, &
       dbcsr_mp_new_transposed, dbcsr_mp_npcols, dbcsr_mp_nprows, &
       dbcsr_mp_numnodes, dbcsr_mp_pgrid, dbcsr_mp_release, &
       dbcsr_mutable_destroy, dbcsr_mutable_init, dbcsr_mutable_instantiated, &
       dbcsr_mutable_new, dbcsr_mutable_release, dbcsr_valid_index
  USE dbcsr_types,                     ONLY: &
       dbcsr_2d_array_type, dbcsr_array_type, dbcsr_data_obj, &
       dbcsr_distribution_obj, dbcsr_imagedistribution_type, dbcsr_iterator, &
       dbcsr_magic_number, dbcsr_meta_size, dbcsr_mp_obj, dbcsr_num_slots, &
       dbcsr_obj, dbcsr_repl_col, dbcsr_repl_full, dbcsr_repl_none, &
       dbcsr_repl_row, dbcsr_slot_blk_p, dbcsr_slot_col_i, dbcsr_slot_nblks, &
       dbcsr_slot_nze, dbcsr_slot_row_p, dbcsr_slot_size, dbcsr_type, &
       dbcsr_type_antisymmetric, dbcsr_type_complex_4, dbcsr_type_complex_8, &
       dbcsr_type_normal, dbcsr_type_real_4, dbcsr_type_real_8, &
       dbcsr_type_real_default, dbcsr_type_symmetric, dbcsr_type_transposed, &
       dbcsr_work_type
  USE distribution_2d_types,           ONLY: distribution_2d_create,&
                                             distribution_2d_get,&
                                             distribution_2d_release,&
                                             distribution_2d_type,&
                                             distribution_2d_write
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8,&
                                             sp
  USE machine,                         ONLY: m_walltime
  USE mathlib,                         ONLY: gcd,&
                                             lcm
  USE message_passing,                 ONLY: mp_allocate,&
                                             mp_alltoall,&
                                             mp_deallocate,&
                                             mp_sum

  !$ USE OMP_LIB
  USE termination,                     ONLY: stop_memory,stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "cp_common_uses.h"
  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_util'

  INTERFACE dbcsr_get_stored_block_info
     MODULE PROCEDURE get_stored_block_info_type, get_stored_block_info_obj
  END INTERFACE

  INTERFACE dbcsr_get_stored_coordinates
     MODULE PROCEDURE get_stored_coordinates_type, get_stored_coordinates_obj
  END INTERFACE

  INTERFACE dbcsr_get_block_index
     MODULE PROCEDURE get_block_index_type, get_block_index_obj
  END INTERFACE

  INTERFACE dbcsr_blk_row_size
     MODULE PROCEDURE dbcsr_blk_row_size_type, dbcsr_blk_row_size_obj
  END INTERFACE

  INTERFACE dbcsr_blk_column_size
     MODULE PROCEDURE dbcsr_blk_column_size_type, dbcsr_blk_column_size_obj
  END INTERFACE

  INTERFACE swap
     MODULE PROCEDURE iswap
  END INTERFACE

  INTERFACE dbcsr_printmat
     MODULE PROCEDURE printmat_s, printmat_d, printmat_c, printmat_z
  END INTERFACE

  INTERFACE addto_array
     MODULE PROCEDURE addto_array_i, addto_array_d, addto_array_ad,&
          addto_array_as, addto_array_ac, addto_array_az
  END INTERFACE

  INTERFACE ensure_array_size
     MODULE PROCEDURE ensure_array_size_i, ensure_array_size_i_2d,&
          ensure_array_size_s, ensure_array_size_d, ensure_array_size_c,&
          ensure_array_size_z
  END INTERFACE

  INTERFACE dbcsr_block_transpose
     MODULE PROCEDURE block_transpose_inplace,&
          block_transpose_copy_d, block_transpose_copy_s,&
          block_transpose_copy_z, block_transpose_copy_c,&
          block_transpose_copy_2d1d_d, block_transpose_copy_2d1d_s,&
          block_transpose_copy_2d1d_z, block_transpose_copy_2d1d_c,&
          block_transpose_copy_1d2d_d, block_transpose_copy_1d2d_s,&
          block_transpose_copy_1d2d_z, block_transpose_copy_1d2d_c,&
          dbcsr_block_transpose_aa
  END INTERFACE

  INTERFACE dbcsr_block_copy
     MODULE PROCEDURE block_copy_2d1d_s, block_copy_2d1d_d,&
          block_copy_2d1d_c, block_copy_2d1d_z,&
          block_copy_1d2d_s, block_copy_1d2d_d,&
          block_copy_1d2d_c, block_copy_1d2d_z
  END INTERFACE

  ! Supports copy between two data areas, or to a data area from a
  ! given explicit array.
  INTERFACE dbcsr_data_set
     MODULE PROCEDURE dbcsr_data_copy_aa, dbcsr_data_set_as,&
          dbcsr_data_set_ad, dbcsr_data_set_ac, dbcsr_data_set_az
  END INTERFACE

  ! Ugly fortran hack
  INTERFACE pointer_view
     MODULE PROCEDURE pointer_view_s, pointer_view_d,&
          pointer_view_c, pointer_view_z
  END INTERFACE

  INTERFACE pointer_set
     MODULE PROCEDURE pointer_set_r2_from_r1_s, pointer_set_r2_from_r1_d,&
                      pointer_set_r2_from_r1_c, pointer_set_r2_from_r1_z
  END INTERFACE

  INTERFACE pointer_replace
     MODULE PROCEDURE pointer_replace_i
  END INTERFACE

  INTERFACE set_block_diagonal
     MODULE PROCEDURE set_block_diagonal_s, set_block_diagonal_d,&
          set_block_diagonal_c, set_block_diagonal_z
  END INTERFACE

  INTERFACE dbcsr_set_diagonal
     MODULE PROCEDURE dbcsr_set_diagonal_s, dbcsr_set_diagonal_d,&
          dbcsr_set_diagonal_z, dbcsr_set_diagonal_c
  END INTERFACE

  INTERFACE dbcsr_destroy_array
     MODULE PROCEDURE dbcsr_destroy_2d_array, dbcsr_destroy_1d_array
  END INTERFACE

  INTERFACE dbcsr_assert
     MODULE PROCEDURE bcsr_assert, dbcsr_int_assert, dbcsr_logical_assert
  END INTERFACE

  ! Main
  PUBLIC :: dbcsr_create, dbcsr_work_create, dbcsr_finalize,&
       dbcsr_checksum, dbcsr_ssum, dbcsr_verify_matrix,&
       dbcsr_repoint_index,&
       dbcsr_work_destroy, dbcsr_print,&
       dbcsr_get_block_index, checker_square_proc,&
       transpose_index_local, dbcsr_make_index_canonical
  ! Distributions
  PUBLIC :: dbcsr_create_dist_r_rot, dbcsr_create_dist_r_unrot,&
       dbcsr_create_dist_l_rot, dbcsr_create_dist_l_unrot,&
       dbcsr_transpose_distribution, dbcsr_create_dist_elemental,&
       dbcsr_create_dist_fullest, dbcsr_create_dist_block_cyclic,&
       create_bl_distribution, dbcsr_transpose_dims, dbcsr_print_dist,&
       dbcsr_reblocking_targets
  ! Images and their distributions
  PUBLIC ::&
       dbcsr_create_col_image_dist,&
       dbcsr_create_image_dist, dbcsr_destroy_image_dist,&
       dbcsr_destroy_array,&
       dbcsr_make_images
  ! distribution_2d_type compatibility
  PUBLIC :: dbcsr_create_dist2d_r, dbcsr_create_dist2d_r_unrot,&
            dbcsr_create_dist2d_transposed, dbcsr_dist2d_to_dist
  ! Low-level printing
  PUBLIC :: dbcsr_printmat, dbcsr_print2dmat
  ! Utility routines
  PUBLIC :: dbcsr_addto_index_array,&
       addto_array, swap, ensure_array_size,&
       dbcsr_set_debug, dbcsr_sleep,&
       checker_tr, bcsr_assert, dbcsr_assert,&
       dbcsr_get_stored_coordinates, dbcsr_get_stored_block_info,&
       dbcsr_find_column,&
       add_work_coordinate,&
       dbcsr_copy_block_sizes, dbcsr_block_transpose, dbcsr_block_copy,&
       dbcsr_copy_sort_data,&
       dbcsr_pack_meta,&
       dbcsr_unpack_i8_2i4,&
       pointer_view, pointer_set, pointer_replace, dbcsr_data_ensure_size,&
       dbcsr_data_set,&
       dbcsr_use_mutable, dbcsr_wm_use_mutable,&
       convert_sizes_to_offsets, convert_offsets_to_sizes,&
       dbcsr_blk_row_size, dbcsr_blk_column_size,&
       dbcsr_add_wm_from_matrix, absmod
  ! Some higher-level operations
  PUBLIC :: dbcsr_set_diagonal



#define DBG IF (dbg) WRITE(*,*)routineN//" ",
#define DBGV IF (bcsr_verbose) WRITE(*,*)routineN//" ",
#define DBGI IF (info) WRITE(*,*)routineN//" ",
#define DEBUG_HEADER        LOGICAL :: dbg, info
#define DEBUG_BODY        dbg = .FALSE. ; CALL dbcsr_set_debug(dbg, info=info)
#define DEBUG_BODYY       dbg = .TRUE. ; CALL dbcsr_set_debug(dbg, info=info)

  LOGICAL, PARAMETER :: bcsr_debug =   .TRUE.
  LOGICAL, PARAMETER :: bcsr_info =    .FALSE.
  LOGICAL, PARAMETER :: bcsr_verbose = .FALSE.

CONTAINS

! *****************************************************************************
!> \brief Determines whether a transpose must be applied
!> \param[in] row   The absolute matrix row.
!> \param[in] column          The absolute matrix column.
! *****************************************************************************
  ELEMENTAL FUNCTION checker_tr(row, column) RESULT(transpose)
    INTEGER, INTENT(IN)                      :: row, column
    LOGICAL                                  :: transpose

    transpose = BTEST(column+row, 0) .EQV. column.GT.row

  END FUNCTION checker_tr

! *****************************************************************************
!> \brief Determines the home process for a given logical matrix element.
!> \note This routine is a more low-level version of
!>       dbcsr_get_stored_coordinate without the setting the row and column
!>       to the stored position.
!> \note It assumes a square matrix.
!> \param[in] row             logical row
!> \param[in] col             logical column
!> \param[in] pgrid           process grid
!> \param[in] row_dist        row distribution
!> \param[in] col_dist        column distribution
!> \retval process            home process of the given element
! *****************************************************************************
  PURE FUNCTION checker_square_proc(row, col, pgrid,&
       row_dist, col_dist) RESULT(process)
    INTEGER, INTENT(IN)                      :: row, col
    INTEGER, DIMENSION(0:, 0:), INTENT(IN)   :: pgrid
    INTEGER, DIMENSION(:), INTENT(IN)        :: row_dist, col_dist
    INTEGER                                  :: process

    IF (.NOT. checker_tr (row, col)) THEN
       process = pgrid (row_dist (row), col_dist(col))
    ELSE
       process = pgrid (row_dist (col), col_dist(row))
    ENDIF
  END FUNCTION checker_square_proc


! *****************************************************************************
!> \brief Sets the correct source matrix, row, column and possible data
!>        transposition for a given matrix and matrix logical row and
!>        column.
!> \param[in] matrix          DBCSR matrix
!> \param[in,out] row         input is logical row, output is lookup row
!> \param[in,out] column      input is logical column, output is lookup
!>                            column
!> \param[in,out] transpose   input: current block state; output: whether the
!>                            data should be considered transposed
!> \param[out] processor      (optional) returns the processor on which this
!>                            block resides
! *****************************************************************************
  PURE SUBROUTINE get_stored_coordinates_type(matrix, row, column, transpose,&
       processor)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(INOUT)                   :: row, column
    LOGICAL, INTENT(INOUT)                   :: transpose
    INTEGER, INTENT(OUT), OPTIONAL           :: processor

!   ---------------------------------------------------------------------------
! SM-compatible way

    IF (PRESENT (processor)) THEN
       IF ((matrix%matrix_type .EQ. dbcsr_type_symmetric&
            .OR. matrix%matrix_type .EQ. dbcsr_type_antisymmetric)&
            .AND. checker_tr(row, column)) THEN
          ! The transposed way.
          processor = dbcsr_distribution_processor (matrix%dist, column, row)
       ELSE
          ! The default way.
          processor = dbcsr_distribution_processor (matrix%dist, row, column)
       ENDIF
    ENDIF
  END SUBROUTINE get_stored_coordinates_type

! \brief Canonical logic
  PURE SUBROUTINE get_stored_canonical(matrix, row, column,&
       transpose, processor)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(INOUT)                   :: row, column
    LOGICAL, INTENT(INOUT)                   :: transpose
    INTEGER, INTENT(OUT), OPTIONAL           :: processor

    INTEGER                                  :: tmp
    LOGICAL                                  :: straight

! The old way

    straight = matrix%matrix_type .NE. dbcsr_type_normal
    straight = (.NOT.&
         (matrix%matrix_type.EQ.dbcsr_type_symmetric&
         .OR.matrix%matrix_type.EQ.dbcsr_type_antisymmetric))&
         .OR.&
         (straight .EQV. .NOT.checker_tr(row, column))
    !transpose = .NOT. straight! .NEQV. transpose
    transpose =  straight .EQV. transpose
    IF (.NOT. straight) THEN
       tmp = row
       row = column
       column = tmp
    ENDIF
    IF (PRESENT (processor)) THEN
       processor = dbcsr_distribution_processor (matrix%dist, row, column)
    END IF
  END SUBROUTINE get_stored_canonical

! *****************************************************************************
!> \brief Sets the correct source matrix, row, column and possible data
!>        transposition for a given matrix and matrix logical row and
!>        column.
!> \param[in] matrix          DBCSR matrix
!> \param[in,out] row         input is logical row, output is lookup row
!> \param[in,out] column      input is logical column, output is lookup
!>                            column
!> \param[out] transpose      whether the data must be transposed
!> \param[out] processor      (optional) returns the processor on which this
!>                            block resides
! *****************************************************************************
  SUBROUTINE get_stored_coordinates_obj(matrix, row, column, transpose, processor)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(INOUT)                   :: row, column
    LOGICAL, INTENT(INOUT)                   :: transpose
    INTEGER, INTENT(OUT), OPTIONAL           :: processor

    CHARACTER(len=*), PARAMETER :: routineN = 'get_stored_coordinates_obj', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF(.NOT.matrix%m%initialized.EQ.dbcsr_magic_number)&
         CALL stop_program (routineP, 'matrix%m not associated') 

    IF (PRESENT (processor)) THEN
       CALL dbcsr_get_stored_coordinates (matrix%m, row, column, transpose, processor)
    ELSE
       CALL dbcsr_get_stored_coordinates (matrix%m, row, column, transpose)
    ENDIF
  END SUBROUTINE get_stored_coordinates_obj


! *****************************************************************************
!> \brief Looks up a block's index given logical coordinates.
!> \param[in] matrix          DBCSR matrix
!> \param[in] row             logical row
!> \param[in] column          logical column
!> \param[out] stored_row     row where block is actually stored
!> \param[out] stored_column  column where block is actually stored
!> \param[out] transpose      whether the data must be transposed
!> \param[out] found          whether the block was found
!> \param[out] block_number   returns the block number of the row
!>                            and column
!> \param[out] data_offset    (optional) data offset for the block; 0 if
!>                            nonexistent
!> \param[out] processor      (optional) returns the processor on which this
!>                            block resides
! *****************************************************************************
  PURE SUBROUTINE get_block_index_type(matrix, row, column, stored_row,&
       stored_column, transpose, found, block_number, data_offset)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(IN)                      :: row, column
    INTEGER, INTENT(OUT)                     :: stored_row, stored_column
    LOGICAL, INTENT(OUT)                     :: transpose, found
    INTEGER, INTENT(OUT)                     :: block_number
    INTEGER, INTENT(OUT), OPTIONAL           :: data_offset

!   ---------------------------------------------------------------------------

    stored_row = row
    stored_column = column
    transpose = .FALSE.
    CALL dbcsr_get_stored_coordinates (matrix, stored_row, stored_column,&
         transpose)
    IF (PRESENT (data_offset)) THEN
       CALL get_stored_block_info_type(matrix, stored_row, stored_column,&
            found, block_number, data_offset)
    ELSE
       CALL get_stored_block_info_type(matrix, stored_row, stored_column,&
            found, block_number)
    ENDIF
  END SUBROUTINE get_block_index_type

! *****************************************************************************
!> \brief Sets the correct source matrix, row, column and possible data
!>        transposition for a given matrix and matrix logical row and
!>        column.
!> \param[in] matrix          DBCSR matrix
!> \param[in] row             logical row
!> \param[in] column          logical column
!> \param[out] stored_row     row where block is actually stored
!> \param[out] stored_column  column where block is actually stored
!> \param[out] transpose      whether the data must be transposed
!> \param[out] found          whether the block was found
!> \param[out] block_number   returns the block number of the row
!>                            and column
!> \param[out] data_offset    (optional) data offset for the block; 0 if
!>                            nonexistent
!> \param[out] processor      (optional) returns the processor on which this
!>                            block resides
! *****************************************************************************
  PURE SUBROUTINE get_block_index_obj(matrix, row, column, stored_row,&
       stored_column, transpose, found, block_number, data_offset)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, column
    INTEGER, INTENT(OUT)                     :: stored_row, stored_column
    LOGICAL, INTENT(OUT)                     :: transpose, found
    INTEGER, INTENT(OUT)                     :: block_number
    INTEGER, INTENT(OUT), OPTIONAL           :: data_offset

!   ---------------------------------------------------------------------------

    IF (PRESENT (data_offset)) THEN
       CALL dbcsr_get_block_index (matrix%m, row, column, stored_row,&
            stored_column, transpose, found, block_number,&
            data_offset=data_offset)
    ELSE
       CALL dbcsr_get_block_index (matrix%m, row, column, stored_row,&
            stored_column, transpose, found, block_number)
    ENDIF
  END SUBROUTINE get_block_index_obj

! *****************************************************************************
!> \brief Returns the index to a queried block, given a real (stored) row and
!>        column
!> \param[in] matrix          bcsr matrix
!> \param[in] row             input is logical row, output is lookup row
!> \param[in] column          input is logical column, output is lookup
!>                            column
!> \param[out] found          whether the block was found
!> \param[out] block_number   returns the block number of the row
!>                            and column
!> \param[out] data_offset    (optional) data offset for the block; 0 if
!>                            nonexistent
!> \param[out] transposed     (optional) whether the block is stored transposed
!>                            according to its position
! *****************************************************************************
  PURE SUBROUTINE get_stored_block_info_type(matrix, row, column,&
       found, block_number, data_offset, transposed)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(IN)                      :: row, column
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT)                     :: block_number
    INTEGER, INTENT(OUT), OPTIONAL           :: data_offset
    LOGICAL, INTENT(OUT), OPTIONAL           :: transposed

    INTEGER                                  :: blk_last, offset

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (matrix%row_p)) THEN
       blk_last = matrix%row_p(row+1)
       IF (blk_last .GT. 0) THEN
          CALL dbcsr_find_column(column, matrix%row_p(row)+1, blk_last,&
               matrix%col_i(1:blk_last), matrix%blk_p(1:blk_last),&
               block_number, found)
       ELSE
          found = .FALSE.
       ENDIF
    ELSE
       found = .FALSE.
    ENDIF
    IF (found) THEN
       IF (PRESENT (data_offset) .OR. PRESENT (transposed)) THEN
          offset = matrix%blk_p(block_number)
       ENDIF
       IF (PRESENT (data_offset)) THEN
          data_offset = ABS (offset)
       ENDIF
       IF (PRESENT (transposed)) THEN
          transposed = offset .LT. 0
       ENDIF
    ELSE
       IF (PRESENT (data_offset)) THEN
          data_offset = 0
       ENDIF
       IF (PRESENT (transposed)) THEN
          transposed = .FALSE.
       ENDIF
    ENDIF
  END SUBROUTINE get_stored_block_info_type

! *****************************************************************************
!> \brief Returns the index to a queried block, given a real (stored) row and
!>        column
!> \param[in] matrix          bcsr matrix
!> \param[in] row             input is logical row, output is lookup row
!> \param[in] column          input is logical column, output is lookup
!>                            column
!> \param[out] found          whether the block was found
!> \param[out] block_number   returns the block number of the row
!>                            and column; 0 if nonexistent
!> \param[out] data_offset    (optional) data offset for the block; 0 if
!>                            nonexistent
!> \param[out] transposed     (optional) whether the block is stored transposed
!>                            according to its position
! *****************************************************************************
  PURE SUBROUTINE get_stored_block_info_obj(matrix, row, column,&
       found, block_number, data_offset, transposed)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, column
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT)                     :: block_number
    INTEGER, INTENT(OUT), OPTIONAL           :: data_offset
    LOGICAL, INTENT(OUT), OPTIONAL           :: transposed

!   ---------------------------------------------------------------------------

    IF (PRESENT (data_offset)) THEN
       IF (PRESENT (transposed)) THEN
          CALL dbcsr_get_stored_block_info (matrix%m, row, column, found,&
               block_number, data_offset, transposed)
       ELSE
          CALL dbcsr_get_stored_block_info (matrix%m, row, column, found,&
               block_number, data_offset)
       ENDIF
    ELSE
       IF (PRESENT (transposed)) THEN
          CALL dbcsr_get_stored_block_info (matrix%m, row, column, found,&
               block_number, transposed=transposed)
       ELSE
          CALL dbcsr_get_stored_block_info (matrix%m, row, column, found,&
               block_number)
       ENDIF
    ENDIF
  END SUBROUTINE get_stored_block_info_obj


! *****************************************************************************
!> \brief Finds the block that has the given column.
!>
!>        If the block having the queried column is found, the blk parameter
!>        is set to this block number and the found parameter is true.
!>        Otherwise found is false and the block number is invalid.
!> \par Index validity
!>      The blk_p array of block pointers is a required parameter to enable
!>      the detection of deleted blocks.
!> \param[in] find_col        column to find
!> \param[in] frst_blk        first block number in row
!> \param[in] last_blk        last block number in row
!> \param[in] col_i           col indices
!> \param[in] blk_p           block pointers
!> \param[out] blk            block number with searched-for column
!> \param[out] found          flag specified whether a block that has the
!>                            correct column was found
! *****************************************************************************
  PURE SUBROUTINE dbcsr_find_column(find_col,frst_blk,last_blk,col_i,blk_p,&
       blk,found)
    INTEGER, INTENT(IN)                      :: find_col, frst_blk, last_blk
    INTEGER, DIMENSION(last_blk), INTENT(IN) :: col_i, blk_p
    INTEGER, INTENT(OUT)                     :: blk
    LOGICAL, INTENT(OUT)                     :: found

    INTEGER                                  :: col, high, low

    found = .FALSE.
    low = frst_blk
    high = last_blk
    blk = (low + high) / 2
    DO WHILE (blk .GE. low .AND. blk .LE. high)
       col = col_i(blk)
       IF(col.EQ.find_col) THEN
          found =  blk_p(blk) .NE. 0
          EXIT
       ELSEIF (col .LT. find_col) THEN
          low = blk+1
       ELSE
          high = blk-1
       ENDIF
       blk = (low + high) / 2
    ENDDO
  END SUBROUTINE dbcsr_find_column


! *****************************************************************************
!> \brief Returns the blocked row size of a row
!>
!> This routine is optimized for speed and no checks are performed.
!> \param[in] matrix          DBCSR matrix
!> \param[in] row             row number  
!> \retval row_size           blocked row size
! *****************************************************************************
  PURE FUNCTION dbcsr_blk_row_size_type (matrix, row) RESULT (row_size)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(IN)                      :: row
    INTEGER                                  :: row_size

    row_size = matrix%row_blk_size%low%data(row)
  END FUNCTION dbcsr_blk_row_size_type
  PURE FUNCTION dbcsr_blk_row_size_obj (matrix, row) RESULT (row_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row
    INTEGER                                  :: row_size

    row_size = matrix%m%row_blk_size%low%data(row)
  END FUNCTION dbcsr_blk_row_size_obj

! *****************************************************************************
!> \brief Returns the blocked column size of a column
!>
!> This routine is optimized for speed and no checks are performed.
!> \param[in] matrix          DBCSR matrix
!> \param[in] column          column number  
!> \retval column_size        blocked row size
! *****************************************************************************
  PURE FUNCTION dbcsr_blk_column_size_type (matrix, column) RESULT (column_size)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(IN)                      :: column
    INTEGER                                  :: column_size

    column_size = matrix%col_blk_size%low%data(column)
  END FUNCTION dbcsr_blk_column_size_type
  PURE FUNCTION dbcsr_blk_column_size_obj (matrix, column) RESULT (column_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: column
    INTEGER                                  :: column_size

    column_size = matrix%m%col_blk_size%low%data(column)
  END FUNCTION dbcsr_blk_column_size_obj


! *****************************************************************************
!> \brief Creates a matrix, allocating the essentials.
!> \par The matrix itself is allocated, as well as the essential parts of
!>      the index. When passed the nze argument, the data is also allocated
!>      to that size.
!> \param[in,out] matrix      new matrix
!> \param[in] dist            distribution_2d distribution
!> \param[in] matrix_type     'N' for normal, 'T' for transposed, 'S' for
!>                            symmetric, and 'A' for antisymmetric
!> \param[in] nblks           (optional) number of blocks
!> \param[in] nze             (optional) number of elements
!> \param[in] data_type       type of data from [rRcC] for single/double
!>                            precision real/complex, default is 'R'
!> \param[in] special         (optional) allocate indices and data using
!>                            special memory
!> \param[in] reuse           (optional) reuses an existing matrix, default
!>                            is to create a fresh one
!> \param[in] mutable_work    uses the mutable data for working and not the
!>                            append-only data; default is append-only
!> \param[in] replication_type     replication to be used for this matrix;
!>                                 default is dbcsr_repl_none
!>                                 \see dbcsr_types.F
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_create(matrix, name, dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, data_type, special, reuse,&
       mutable_work, replication_type, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    CHARACTER(len=*), INTENT(IN)             :: name
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    CHARACTER, INTENT(IN)                    :: matrix_type
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_size, col_blk_size
    INTEGER, INTENT(IN), OPTIONAL            :: nblks, nze
    CHARACTER, INTENT(IN), OPTIONAL          :: data_type
    LOGICAL, INTENT(IN), OPTIONAL            :: special, reuse, mutable_work
    CHARACTER, INTENT(IN), OPTIONAL          :: replication_type
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: mt
    INTEGER                                  :: my_nblks, my_nze, stat
    INTEGER, DIMENSION(dbcsr_meta_size)      :: new_meta
    LOGICAL                                  :: hijack

!   ---------------------------------------------------------------------------

    hijack = .FALSE.
    IF (PRESENT (reuse)) THEN
       hijack = reuse
    ELSE
       IF (matrix%m%initialized.EQ.dbcsr_magic_number) THEN
          ! Reuse matrix only if has actually been allocated.
          IF (ASSOCIATED (matrix%m%index)) THEN
             hijack = .TRUE.
          ELSE
             hijack = .FALSE.
          ENDIF
       ELSE
          CALL cp_assert (matrix%m%initialized .EQ. 0,&
               cp_fatal_level, cp_caller_error, routineN,&
               "Matrix may not have been initialized with dbcsr_init")
          hijack = .FALSE.
       ENDIF
    ENDIF
    IF (.NOT.hijack) THEN
       CALL dbcsr_init (matrix%m)
       matrix%m%refcount = 1
    ENDIF
    ! Mark matrix index as having an invalid index.
    matrix%m%valid = .FALSE.
    matrix%m%name = name
    ! Sets the type of matrix building/modifying work structures.
    IF (PRESENT (mutable_work)) THEN
       matrix%m%work_mutable = mutable_work
    ELSE
       matrix%m%work_mutable = .FALSE.
    ENDIF
    ! Sets the correct data type.
    IF (PRESENT (data_type)) THEN
       SELECT CASE (data_type)
       CASE (dbcsr_type_real_4)
          matrix%m%data_type = dbcsr_type_real_4
          CALL cp_unimplemented_error (routineN,&
               "Single precision real matrices not yet tested or fully supported",&
               error_level=cp_warning_level)
       CASE (dbcsr_type_real_8)
          matrix%m%data_type = dbcsr_type_real_8
       CASE (dbcsr_type_complex_4)
          matrix%m%data_type = dbcsr_type_complex_4
          CALL cp_unimplemented_error (routineN,&
               "Complex matrices not yet tested or fully supported",&
               error_level=cp_warning_level)
       CASE (dbcsr_type_complex_8)
          matrix%m%data_type = dbcsr_type_complex_8
          CALL cp_unimplemented_error (routineN,&
               "Complex matrices not yet tested or fully supported",&
               error_level=cp_warning_level)
       CASE DEFAULT
         CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error,&
              routineN, "Invalid matrix type")
      END SELECT
    ELSE
       matrix%m%data_type = dbcsr_type_real_default
    ENDIF
    IF (hijack) THEN
       ! Release/deallocate elements that are replaced or not needed
       ! by the new matrix. This is similar to what dbcsr_destroy
       ! does, except that it keeps the index and data.
       CALL array_release (matrix%m%row_blk_size)
       CALL array_release (matrix%m%col_blk_size)
       CALL dbcsr_distribution_release (matrix%m%dist)
       IF (ASSOCIATED (matrix%m%dist2d)) &
            CALL distribution_2d_release (matrix%m%dist2d, error)
       IF (ASSOCIATED (matrix%m%wms)) THEN
          CALL dbcsr_work_destroy_all(matrix%m)
       ENDIF
    ELSE
       ! Invalidate index
       NULLIFY(matrix%m%index)
       ! Invalidate data
       CALL dbcsr_data_init (matrix%m%data_area)
    ENDIF
    ! These are always invalidated.
    NULLIFY(matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p)
    matrix%m%row_blk_size = row_blk_size
    CALL array_hold (matrix%m%row_blk_size)
    IF (array_size (matrix%m%row_blk_size) .GT. 0) THEN
       matrix%m%max_rbs = MAXVAL (array_data (matrix%m%row_blk_size))
    ELSE
       matrix%m%max_rbs = 0
    ENDIF
    matrix%m%col_blk_size = col_blk_size
    CALL array_hold (matrix%m%col_blk_size)
    IF (array_size (matrix%m%col_blk_size) .GT. 0) THEN
       matrix%m%max_cbs = MAXVAL (array_data (matrix%m%col_blk_size))
    ELSE
       matrix%m%max_cbs = 0
    ENDIF
    matrix%m%dist = dist
    CALL dbcsr_distribution_hold (matrix%m%dist)
    ! Set up some data.
    CALL meta_from_dist (new_meta, dist, row_blk_size, col_blk_size)
    matrix%m%nblkrows_total  = new_meta(3)
    matrix%m%nblkcols_total  = new_meta(4)
    matrix%m%nfullrows_total = new_meta(5)
    matrix%m%nfullcols_total = new_meta(6)
    matrix%m%nblkrows_local  = new_meta(7)
    matrix%m%nblkcols_local  = new_meta(8)
    matrix%m%nfullrows_local = new_meta(9)
    matrix%m%nfullcols_local = new_meta(10)
    my_nze = 0; IF (PRESENT (nze)) my_nze = nze
    my_nblks = 0; IF (PRESENT (nblks)) my_nblks = nblks
    matrix%m%nblks = 0
    matrix%m%nze = 0
    IF (PRESENT (special)) THEN
       matrix%m%special_memory = special
    ELSE
       matrix%m%special_memory = .FALSE.
    ENDIF
    IF (PRESENT (replication_type)) THEN
       CALL cp_assert (replication_type .EQ. dbcsr_repl_none&
            .OR. replication_type .EQ. dbcsr_repl_full&
            .OR. replication_type .EQ. dbcsr_repl_row&
            .OR. replication_type .EQ. dbcsr_repl_col,&
            cp_failure_level, cp_wrong_args_error, routineN,&
            "Invalid replication type '"//replication_type//"'")
       CALL cp_assert (replication_type .NE. dbcsr_repl_row&
            .AND. replication_type .NE. dbcsr_repl_col,&
            cp_warning_level, cp_unimplemented_error_nr, routineN,&
            "Row and column replication not fully supported")
       matrix%m%replication_type = replication_type
    ELSE
       matrix%m%replication_type = dbcsr_repl_none
    ENDIF
    ! CP2K SM compatibility
    NULLIFY(matrix%m%dist2d)
    !
    ! Setup the data
    IF (.NOT. hijack) THEN
       CALL dbcsr_data_new (matrix%m%data_area, matrix%m%data_type, my_nze,&
            special_memory=matrix%m%special_memory)
       IF (matrix%m%special_memory) THEN
          CALL mp_allocate(matrix%m%index, dbcsr_num_slots, stat=stat)
          IF (stat /= 0) CALL stop_memory (routineP, 'matrix%m%index', my_nze)
       ELSE
          ALLOCATE(matrix%m%index(dbcsr_num_slots), stat=stat)
          IF (stat /= 0) CALL stop_memory (routineP, 'matrix%m%index', my_nze)
       ENDIF
    ENDIF
    CALL cp_assert (LBOUND (matrix%m%index, 1) .LE. 1&
         .AND. UBOUND (matrix%m%index, 1) .GE. dbcsr_num_slots,&
         cp_fatal_level, cp_internal_error, routineN,&
         "Index is not large enough")
    matrix%m%index(1:dbcsr_num_slots) = 0
    matrix%m%index(dbcsr_slot_size) = dbcsr_num_slots
    !
    matrix%m%sym_negation = .FALSE.
    SELECT CASE (matrix_type)
    CASE ('N', 'n')
       mt = dbcsr_type_normal
    CASE ('T', 't')
       mt = dbcsr_type_transposed
    CASE ('S', 's')
       mt = dbcsr_type_symmetric
    CASE ('A', 'a')
       mt = dbcsr_type_antisymmetric
       matrix%m%sym_negation = .TRUE.
    CASE DEFAULT
       CALL cp_assert(.FALSE., cp_failure_level,&
            cp_wrong_args_error, routineP, "Invalid matrix type.")
    END SELECT
    matrix%m%matrix_type = mt
!$  CALL OMP_INIT_LOCK (matrix%m%modification_lock)
  END SUBROUTINE dbcsr_create


! *****************************************************************************
!> \brief Initializes one work matrix
!> \param[out] wm             initialized work matrix
!> \param[in] nblks_guess     (optional) estimated number of blocks
!> \param[in] sizedata_guess  (optional) estimated size of data
! *****************************************************************************
  SUBROUTINE dbcsr_init_wm (wm, data_type, nblks_guess, sizedata_guess)
    TYPE(dbcsr_work_type), INTENT(OUT)       :: wm
    CHARACTER, INTENT(IN)                    :: data_type
    INTEGER, INTENT(IN), OPTIONAL            :: nblks_guess, sizedata_guess

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_wm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nblks, stat

!   ---------------------------------------------------------------------------

    wm%lastblk = 0
    wm%datasize = 0
    ! Index
    IF(PRESENT(nblks_guess)) THEN
       nblks = nblks_guess
       ALLOCATE(wm%row_i(nblks), stat=stat)
       IF (stat /= 0) CALL stop_memory(routineP, 'wm%row_i', nblks)
       ALLOCATE(wm%col_i(nblks), stat=stat)
       IF (stat /= 0) CALL stop_memory(routineP, 'wm%col_i', nblks)
       ALLOCATE(wm%blk_p(nblks), stat=stat)
       IF (stat /= 0) CALL stop_memory(routineP, 'wm%blk_p', nblks)
    ELSE
       NULLIFY (wm%row_i, wm%col_i, wm%blk_p)
       !nblks = CEILING (REAL (matrix%nblkrows_local * matrix%nblkcols_local)&
       !     / REAL (dbcsr_mp_numnodes (dbcsr_distribution_mp (matrix%dist))))
    ENDIF
    ! Data
    CALL dbcsr_data_init (wm%data_area)
    IF(PRESENT(sizedata_guess)) THEN
       CALL dbcsr_data_new (wm%data_area, data_type,&
            data_size=sizedata_guess)
    ELSE
       CALL dbcsr_data_new (wm%data_area, data_type)
    ENDIF
    CALL dbcsr_mutable_init (wm%mutable)
  END SUBROUTINE dbcsr_init_wm


! *****************************************************************************
!> \brief Creates a the working matrix(es) for a DBCSR matrix.
!> \param[out] matrix         new matrix
!> \param[in] nblks_guess     (optional) estimated number of blocks
!> \param[in] sizedata_guess  (optional) estimated size of data
!> \param[in] n               (optional) number work matrices to create,
!>                            default is 1
!> \param[in] work_mutable    (optional) use mutable work type, default is
!>                            what was specified in create
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_work_create(matrix, nblks_guess, sizedata_guess, n,&
       error, work_mutable)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN), OPTIONAL            :: nblks_guess, sizedata_guess, n
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(in), OPTIONAL            :: work_mutable

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_work_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iw, nw, ow
    LOGICAL                                  :: wms_new, wms_realloc
    TYPE(dbcsr_work_type), DIMENSION(:), &
      POINTER                                :: wms

!   ---------------------------------------------------------------------------

    CALL cp_assert (matrix%m%initialized.EQ.dbcsr_magic_number,&
         cp_fatal_level, cp_caller_error,&
         routineN, "Can not create work matrix for nonexisting matrix object.")
    IF (PRESENT (n)) THEN
       nw = n
    ELSE
       nw = 1
!$     nw = omp_get_max_threads()
    ENDIF
!$  CALL dbcsr_modify_lock (matrix)
    wms_new = .NOT. ASSOCIATED (matrix%m%wms)
    wms_realloc = .FALSE.
    nw = 1
!$    nw = omp_get_max_threads()
    IF (ASSOCIATED (matrix%m%wms)) THEN
       ow = SIZE(matrix%m%wms)
       CALL cp_assert (ow .GE. nw, cp_warning_level,&
            cp_internal_error, routineN,&
            "Number of work matrices less than threads.")
       IF (ow .LT. nw) wms_realloc = .TRUE.
    ENDIF
    IF (PRESENT (work_mutable)) THEN
       matrix%m%work_mutable = work_mutable
    ENDIF
    IF (wms_realloc) THEN
       ALLOCATE (wms(nw))
       wms(1:ow) = matrix%m%wms(1:ow)
       DEALLOCATE (matrix%m%wms)
       matrix%m%wms => wms
       DO iw = ow+1, nw
          CALL dbcsr_init_wm (matrix%m%wms(iw), matrix%m%data_type,&
               nblks_guess=nblks_guess, sizedata_guess=sizedata_guess)
          IF (matrix%m%work_mutable) &
               CALL dbcsr_mutable_new (matrix%m%wms(iw)%mutable,&
               dbcsr_get_data_type (matrix))
       END DO
    ENDIF
    IF (wms_new) THEN
       ALLOCATE (matrix%m%wms(nw))
       DO iw = 1, nw
          CALL dbcsr_init_wm (matrix%m%wms(iw), matrix%m%data_type,&
               nblks_guess=nblks_guess, sizedata_guess=sizedata_guess)
          IF (matrix%m%work_mutable) &
               CALL dbcsr_mutable_new (matrix%m%wms(iw)%mutable,&
               dbcsr_get_data_type (matrix))
       END DO
    ENDIF
!$OMP FLUSH (matrix)
!$  CALL dbcsr_modify_unlock (matrix)
  END SUBROUTINE dbcsr_work_create


! *****************************************************************************
!> \brief Creates the final dbcsr_type matrix from the working matrix.
!>
!> If a matrix has a valid index, then nothing is done.
!> \note The columns and rows must be ordered! We don't support appending
!>       data to an existing matrix.
!> \param[in,out] matrix      final matrix
!> \param[in] resort          whether the indices should be sorted, default
!>                            is true
!> \param[in] reshuffle       whether the data should be reshuffled,
!>                            default is false
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_finalize(matrix, resort, reshuffle, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: resort, reshuffle
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: all_nblks, all_nze, i, j, &
                                                timing_handle
    LOGICAL                                  :: multidata, sort_data, &
                                                sort_index

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL cp_assert (dbcsr_is_initialized (matrix),&
         cp_fatal_level, cp_caller_error,&
         routineN, "Can not finalize uninitialized matrix.")
    IF (dbcsr_valid_index(matrix)) THEN
       IF (dbg) WRITE(*,'(2(1X,A))')&
            "No need to finalize a valid matrix, skipping."
!vw       CALL cp_assert (.NOT. ASSOCIATED (matrix%m%wms), cp_warning_level,&
!vw            cp_internal_error, routineN,&
!vw            "Valid matrix should not contain work data")
       IF (ASSOCIATED (matrix%m%wms)) &
            CALL dbcsr_work_destroy_all(matrix%m)
       RETURN
    ENDIF
    CALL timeset(routineN, timing_handle)
    IF (dbg) WRITE(*,'(A,1X,A,2(1X,I9))')routineN//' finalizing matrix ',&
         matrix%m%name, matrix%m%nze, matrix%m%nblks
    !CALL cp_assert (matrix%m%nze.EQ.0 .EQV. matrix%m%nblks.EQ.0,&
    !     cp_warning_level, cp_unimplemented_error_nr, routineN,&
    !     "Matrix data size not consistent with number of blocks.")
    ! If there is data in the matrix proper, create a work matrix from it.
    IF (matrix%m%nze .GT. 0) THEN
       CALL dbcsr_add_wm_from_matrix (matrix, 1)
    ENDIF
    IF (PRESENT (reshuffle)) THEN
       sort_data = reshuffle
    ELSE
       sort_data = .FALSE.
    ENDIF
    sort_index = .TRUE.
    IF (PRESENT (resort)) sort_index = resort
    IF (ASSOCIATED (matrix%m%wms)) THEN
       multidata = SIZE (matrix%m%wms) .GT. 1
    ELSE
       multidata = .FALSE.
    ENDIF
    IF (.NOT. (matrix%m%nze .GT. 0 .OR. ASSOCIATED (matrix%m%wms))) THEN
       IF (dbg) THEN
          CALL cp_assert (matrix%m%nze .GT. 0 .OR. ASSOCIATED (matrix%m%wms),&
               cp_warning_level, cp_caller_error, routineN,&
               "Nothing to be done.")
       ENDIF
       ! Should short-circut and just create an empty index.
       CALL dbcsr_work_create (matrix)
    ENDIF
    CALL cp_assert (ASSOCIATED (matrix%m%wms), cp_fatal_level, cp_caller_error,&
         routineN, "Can not finalize matrix containing no work matrices.")
    CALL cp_assert (SIZE (matrix%m%wms) .GT. 0, cp_failure_level,&
         cp_caller_error, routineN,&
         "Can not finalize matrix containing no work matrices.")
    ! If there are multiple work matrices, the indices have to be sorted to
    ! enable merging the index.
    IF (sort_index .OR. multidata) THEN
       IF (multidata) THEN
          CALL dbcsr_sort_many_indices (matrix%m)
          sort_data = .TRUE.
       ELSE
          IF (.NOT.dbcsr_wm_use_mutable (matrix%m%wms(1))&
               .AND. matrix%m%wms(1)%lastblk .GT. 0) THEN
             CALL dbcsr_sort_indices(matrix%m%wms(1)%lastblk,&
                  matrix%m%wms(1)%row_i,&
                  matrix%m%wms(1)%col_i, matrix%m%wms(1)%blk_p)
          ENDIF
       ENDIF
    ENDIF
    ! Count how large the new matrix will have to be.
    all_nblks = 0
    all_nze = 0
    DO i = 1, SIZE (matrix%m%wms)
       all_nblks = all_nblks + matrix%m%wms(i)%lastblk
       all_nze = all_nze + matrix%m%wms(i)%datasize
    ENDDO
    IF (dbg) THEN
       DO i = 1, SIZE (matrix%m%wms)
          WRITE(*,'(A,4(1X,I4))')"Sorted work matrix",i, matrix%m%wms(i)%lastblk,&
               matrix%m%wms(i)%datasize,&
               dbcsr_get_data_size (matrix%m%wms(i)%data_area)
          IF (.NOT.dbcsr_wm_use_mutable (matrix%m%wms(i))) THEN
             DO j = 1, matrix%m%wms(i)%lastblk
                WRITE(*,'(A,I4,":",3(1X,I5))')"Sorted",j,&
                     matrix%m%wms(i)%row_i(j),&
                     matrix%m%wms(i)%col_i(j), matrix%m%wms(i)%blk_p(j)
             ENDDO
          ENDIF
       ENDDO
    ENDIF
    matrix%m%nblks = all_nblks
    matrix%m%nze = all_nze
    matrix%m%index(dbcsr_slot_nblks) = matrix%m%nblks
    matrix%m%index(dbcsr_slot_nze) = matrix%m%nze
    DBG routineP//' index',matrix%m%index(:dbcsr_num_slots)
    DBG routineP//' # blocks=',matrix%m%nblks, all_nblks
    DBG routinep//' # nze=',matrix%m%nze, all_nze
    ! We first figure out how large the indices should be and allocate them.
    CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_row_p,&
         reservation=matrix%m%nblkrows_total+1, extra=all_nblks*2)
    CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_col_i,&
         reservation=all_nblks)
    CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_blk_p,&
         reservation=all_nblks)
    CALL dbcsr_repoint_index(matrix%m)
    matrix%m%row_p(1) = 0
    !CALL dbcsr_data_ensure_size (matrix%m%data_area, all_nze)
    ! Now merge the index and data.
    matrix%m%nze = all_nze ! Specifies size of data array to be allocated
    CALL dbcsr_merge_data (matrix%m, sort_data)
    !
    ! Clean up.
    CALL dbcsr_work_destroy_all(matrix%m)
    matrix%m%valid = .TRUE.
    IF (dbg) CALL dbcsr_verify_matrix (matrix)
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_finalize


! *****************************************************************************
!> \brief Creates a work matrix from the data present in a finalized matrix.
!> \param[in,out] matrix      DBCSR matrix
!> \param[in] nwms            Number of work matrices (i.e., threads)
! *****************************************************************************
  SUBROUTINE dbcsr_add_wm_from_matrix(matrix, nwms)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: nwms

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_wm_from_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n, timing_handle
    TYPE(dbcsr_work_type), DIMENSION(:), &
      POINTER                                :: new_wms

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    DBG "Adding to an existing matrix, wee!"
    IF (ASSOCIATED (matrix%m%wms)) THEN
       n = SIZE (matrix%m%wms)
       IF (n .EQ. 0) DEALLOCATE (matrix%m%wms)
    ELSE
       n = 0
    ENDIF
    IF (n .EQ. 0) THEN
       ! Just create it
       CALL dbcsr_work_create (matrix, matrix%m%nblks, matrix%m%nze, n=nwms,&
            work_mutable=.FALSE.)
       CALL dbcsr_fill_wm_from_matrix (matrix%m%wms, matrix%m, nwms=nwms,&
            offset=0)
    ELSE
       ALLOCATE (new_wms (1:n+nwms))
       new_wms(1:n) = matrix%m%wms(1:n)
       DEALLOCATE (matrix%m%wms)
       matrix%m%wms => new_wms
       DO i = 1, nwms
          CALL dbcsr_init_wm (matrix%m%wms(n+i), dbcsr_get_data_type (matrix),&
               nblks_guess=dbcsr_get_num_blocks (matrix))
       ENDDO
       CALL dbcsr_fill_wm_from_matrix (matrix%m%wms, matrix%m,&
            nwms=nwms, offset=n)
    ENDIF
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_add_wm_from_matrix


! *****************************************************************************
!> \brief Fills index and data of the work matrix from the
!>        previously-finalized one.
!> \param[out] wm             the work matrix to fill
!> \param[in,out] matrix      DBCSR matrix
!> \param[in] nwms            number of work matrices
!> \param[in] offset          first work matrix to use
! *****************************************************************************
  SUBROUTINE dbcsr_fill_wm_from_matrix(wm, matrix, nwms, offset)
    TYPE(dbcsr_work_type), DIMENSION(:), &
      INTENT(INOUT)                          :: wm
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, INTENT(IN)                      :: nwms, offset

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_fill_wm_from_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_p, col, i, nblks, &
                                                row, timing_handle, which_wm
    INTEGER, DIMENSION(:), POINTER           :: tdist
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_obj)                          :: m

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    nblks = matrix%nblks
    m%m = matrix
    IF (nwms .GT. 1) THEN
       tdist => array_data (matrix%dist%d%thread_dist)
       CALL cp_assert (ASSOCIATED (tdist), cp_fatal_level, cp_internal_error,&
            routineN, "thread distribution must be defined in threaded env.")
       WRITE(*,*)routineN//" Threading fill", tdist
    ENDIF
    CALL dbcsr_iterator_start (iter, m)
    IF (dbg) WRITE(*,*)routineN//" There are",nblks,"blocks."
    i = 1
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       IF (dbg) WRITE(*,*)routineN//" Getting block", i,"/",nblks
       CALL dbcsr_iterator_next_block (iter, row, col, blk=blk, blk_p=blk_p)
       IF (dbg) WRITE(*,*)routineN//" With coor", row, col, "at", blk_p
       !which_wm = MOD (row, nwms) + 1 + offset
       IF (nwms .GT. 1) THEN
          which_wm = tdist (row) + offset
       ELSE
          which_wm = offset+1
       ENDIF
       CALL add_work_coordinate (wm(which_wm), row, col, blk_p)
       i = i+1
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    !wm%data_area = matrix%data_area
    !CALL dbcsr_data_hold (wm%data_area)
    DO which_wm = 1+offset, nwms+offset
       CALL dbcsr_data_ensure_size (wm(which_wm)%data_area,&
            dbcsr_get_data_size (matrix%data_area))
       CALL dbcsr_data_copyall (wm(which_wm)%data_area, matrix%data_area)
       wm(which_wm)%datasize = dbcsr_get_data_size (wm(which_wm)%data_area)
    ENDDO
    IF (dbg) WRITE(*,*)routineN//" Size is",wm%datasize, wm%lastblk
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_fill_wm_from_matrix


! *****************************************************************************
!> \brief Deallocates and destroys a work matrix.
!> \param[in,out] wm          work matrix
! *****************************************************************************
  SUBROUTINE dbcsr_work_destroy(wm)
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: wm

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_work_destroy', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (wm%row_i)) THEN
       DEALLOCATE(wm%row_i)
       NULLIFY (wm%row_i)
    ENDIF
    IF (ASSOCIATED (wm%col_i)) THEN
       DEALLOCATE(wm%col_i)
       NULLIFY (wm%row_i)
    ENDIF
    IF (ASSOCIATED (wm%blk_p)) THEN
       DEALLOCATE(wm%blk_p)
       NULLIFY (wm%blk_p)
    ENDIF
    CALL dbcsr_data_release (wm%data_area)
    CALL dbcsr_mutable_destroy (wm%mutable)
  END SUBROUTINE dbcsr_work_destroy


! *****************************************************************************
!> \brief Deallocates and destroys a work matrix.
!> \param[in,out] wm          work matrix
!> \param[in,out] error       cp2k error
!> \param keepdata    do not deallocate data
!> \param keepfinal   do not destroy the final, non-work matrix
!> \param keepfinaldata       do not destroy the data in the final,
!>                            non-work matrix
! *****************************************************************************
  SUBROUTINE dbcsr_work_destroy_all(m)
    TYPE(dbcsr_type), INTENT(INOUT)          :: m

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_work_destroy_all', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (m%wms), cp_warning_level, cp_caller_error,&
         routineN, "Want to destroy nonexisting work matrices.")
    IF (ASSOCIATED (m%wms)) THEN
       DO i = 1, SIZE (m%wms)
          CALL dbcsr_work_destroy (m%wms(i))
       ENDDO
       DEALLOCATE (m%wms)
       NULLIFY (m%wms)
    ENDIF
  END SUBROUTINE dbcsr_work_destroy_all


! *****************************************************************************
!> \brief Converts sizes to offsets
!>
!> \param[in] sizes           array with sizes
!> \param[out] offsets_start  offsets of starts
!> \param[out] offsets_stop   (optional) offsets of ends
! *****************************************************************************

  PURE SUBROUTINE convert_sizes_to_offsets (sizes,&
       offsets_start, offsets_stop)
    INTEGER, DIMENSION(:), INTENT(IN)        :: sizes
    INTEGER, DIMENSION(:), INTENT(OUT)       :: offsets_start
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: offsets_stop

    CHARACTER(len=*), PARAMETER :: routineN = 'convert_sizes_to_offsets', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n

!   ---------------------------------------------------------------------------

    n = SIZE (sizes)
    !IF (bcsr_debug) THEN
    !   CALL cp_assert (SIZE (offsets_start) .GE. n, cp_failure_level,&
    !        cp_caller_error, routineN, "Start offset array too small")
    !   CALL cp_assert (SIZE (offsets_start) .EQ. n, cp_warning_level,&
    !        cp_caller_error, routineN, "Start offset array too big")
    !   IF (PRESENT (offsests_stop)) THEN
    !      CALL cp_assert (SIZE (offsets_stop) .GE. n, cp_failure_level,&
    !           cp_caller_error, routineN, "Stop offset array too small")
    !      CALL cp_assert (SIZE (offsets_stop) .EQ. n, cp_warning_level,&
    !           cp_caller_error, routineN, "Stop offset array too big")
    !   ENDIF
    !ENDIF
    IF (n .GT. 0) THEN
       offsets_start(1) = 1
       IF (PRESENT (offsets_stop)) offsets_stop(1) = sizes(1)
       IF (.NOT. PRESENT (offsets_stop)) THEN
          DO i = 2, n
             offsets_start(i) = offsets_start(i-1) + sizes(i-1) + 1
          ENDDO
       ELSE
          DO i = 2, n
             offsets_start(i) = offsets_start(i-1) + sizes(i-1)
             offsets_stop(i) = offsets_stop(i-1) + sizes(i)
          ENDDO
       ENDIF
    ENDIF
  END SUBROUTINE convert_sizes_to_offsets

! *****************************************************************************
!> \brief Converts offsets to sizes
!>
!> If the offsets of ends are not given, then the array of sizes is assumed
!> to be one greater than the desired sizes.
!>
!> \param[in] offsets_start  offsets of starts
!> \param[out] sizes         array with sizes
!> \param[in] offsets_stop   (optional) offsets of ends
! *****************************************************************************

  PURE SUBROUTINE convert_offsets_to_sizes (offsets_start, sizes, offsets_stop)
    INTEGER, DIMENSION(:), INTENT(IN)        :: offsets_start
    INTEGER, DIMENSION(:), INTENT(OUT)       :: sizes
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: offsets_stop

    CHARACTER(len=*), PARAMETER :: routineN = 'convert_offsets_to_sizes', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n

!   ---------------------------------------------------------------------------

    n = SIZE (offsets_start)
    IF (PRESENT (offsets_stop)) THEN
       sizes(:) = offsets_stop(:) - offsets_start(:) + 1
    ELSE
       IF (n .GT. 1) THEN
          DO i = 1, n-1
             sizes(i) = sizes(i+1) - sizes(i)
          ENDDO
       ENDIF
    ENDIF
  END SUBROUTINE convert_offsets_to_sizes


! *****************************************************************************
!> \brief Determines the correct transposed type of a DBCSR matrix.
!> \param[out] new_type       new matrix type
!> \param[in] old_type        current matrix_type
! *****************************************************************************
  ELEMENTAL SUBROUTINE dbcsr_transposed_type (new_type, old_type)
    CHARACTER, INTENT(OUT)                   :: new_type
    CHARACTER, INTENT(IN)                    :: old_type

    SELECT CASE (old_type)
    CASE ('N')
       new_type = 'T'
    CASE ('T')
       new_type = 'N'
    CASE DEFAULT
       new_type = old_type
    END SELECT
  END SUBROUTINE dbcsr_transposed_type


! *****************************************************************************
!> \brief Re-indexes row_p and blk_i according to columns.
!> \par The re-indexing is equivalent to a local-only transpose.
!> \param[out] new_col_p      new column pointer
!> \param[out] new_row_i      new row index
!> \param[in] old_row_p       old row pointer
!> \param[in] old_col_i       old column index
!> \param[out] new_blk_p      (optional) new block pointer
!> \param[in] old_blk_p       (optional) old block pointer
! *****************************************************************************
  SUBROUTINE transpose_index_local (new_col_p, new_row_i, old_row_p,&
       old_col_i, new_blk_p, old_blk_p)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: new_col_p, new_row_i
    INTEGER, DIMENSION(:), INTENT(IN)        :: old_row_p, old_col_i
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: new_blk_p
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: old_blk_p

    CHARACTER(len=*), PARAMETER :: routineN = 'transpose_index_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: curcol, mincol, nblks, &
                                                new_ncols, old_max_col, &
                                                old_nrows, old_row, prev_col, &
                                                prev_nblks
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: row_p
    LOGICAL                                  :: column_le

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY

1323 FORMAT (I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5)
    new_ncols = SIZE (new_col_p)-1
    old_nrows = SIZE (old_row_p)-1
    ALLOCATE (row_p (SIZE (old_row_p)))
    row_p(:) = old_row_p(:)+1
    IF (SIZE (old_col_i) .GT. 0) THEN
       old_max_col = MAXVAL (old_col_i)
    ELSE
       old_max_col = 0
    ENDIF
    mincol = old_max_col
    IF (dbg) THEN
       WRITE(*,*)'old index'
       WRITE(*,1323)old_row_p
       WRITE(*,1323)old_col_i
       IF (PRESENT (old_blk_p)) THEN
          WRITE(*,*)'old blkp'
          WRITE(*,1323)old_blk_p
       ENDIF
    ENDIF
    ! Scan through the columns of the given matrix.
    ! row_p is an index, for each row, of the current column scanning position
    ! appropriate to the column being scanned.
    ! mincol is use to jump ahead for very sparse matrices
    nblks = 0
    prev_nblks = 0
    new_col_p(1) = 0
    prev_col = 0
    curcol = 1
    DO WHILE (curcol .LE. new_ncols)
       DBG 'curcol', curcol
       IF (dbg) THEN
          WRITE(*,*)'rowp'
          WRITE(*,1323)row_p
       ENDIF
       DO old_row = 1, old_nrows
          ! Condition 1:
          ! row_p(old_row) .LE. old_row_p(old_row+1)
          !           => there are still blocks in the row
          ! Condition 2:
          ! old_col_i(row_p(old_row)) .LE. curcol
          !           => the examined column in this row is <= to the
          !              curcol current column
          column_le = row_p(old_row) .LE. old_row_p(old_row+1)
          IF (column_le) column_le = old_col_i(row_p(old_row)) .LE. curcol
          ! In every row, try to find the column appropriate to the current
          ! column.
          DO WHILE (column_le)
             IF (old_col_i(row_p(old_row)) .EQ. curcol) THEN
                nblks = nblks + 1
                DBG 'Adding block',nblks,'at old row/new col',old_row
                !new_col_p(curcol+1) = i = new_col_p(curcol+1) + 1
                new_row_i(nblks) = old_row
                IF (PRESENT (new_blk_p) .AND. PRESENT (old_blk_p)) THEN
                   new_blk_p (nblks) = old_blk_p(row_p(old_row))
                ENDIF
             ENDIF
             row_p(old_row) = row_p(old_row) + 1
             ! We've come to the end of this row
             column_le = row_p(old_row) .LE. old_row_p(old_row+1)
             IF (column_le) column_le = old_col_i(row_p(old_row)) .LE. curcol
          ENDDO
          IF (row_p(old_row) .LE. old_row_p(old_row+1)) THEN
             mincol = MIN (mincol, old_col_i(row_p(old_row)))
             !DBG 'for row',old_row,'leftoff column is',old_col_i(row_p(old_row))
          ELSE
             !DBG 'for row',old_row,'forced column is ',old_col_i(MIN(row_p(old_row), old_max_col))
          ENDIF
       ENDDO
       !! Fill in missing indices.
       DBG 'prev_col, curcol, mincol, nblks', prev_col, curcol, mincol, nblks
       IF (prev_col+1 .LE. mincol) THEN
          DBG 'Skipped!'
       ENDIF
       !new_col_p(prev_col+2:MIN(mincol+1,new_ncols+1)) = nblks
       new_col_p(prev_col+1:MIN(mincol+1,new_ncols+1)-1) = prev_nblks
       new_col_p(MIN(mincol+1,new_ncols+1)) = nblks
       prev_nblks = nblks
       prev_col = curcol
       curcol = MAX (mincol, curcol+1)
       !curcol = mincol
       mincol = mincol+1
       IF (.NOT. (prev_col .LT. curcol)) THEN
          DBG 'BAD new state: prev_col, curcol, mincol',&
               prev_col, curcol, mincol
          STOP
       ENDIF
       !IF(curcol+2 .LE. mincol .AND. curcol+2 .LE. old_max_col+1) THEN
       !   new_col_p(curcol+2:mincol) = new_col_p(curcol+1)
       !END IF
       !curcol = MAX(mincol,curcol+1)
       !mincol = old_max_col
    ENDDO
    new_col_p(new_ncols+1) = nblks
    IF (dbg) THEN
       WRITE(*,*)'new index'
       WRITE(*,1323)new_col_p
       WRITE(*,1323)new_row_i
       IF (PRESENT (new_blk_p)) THEN
          WRITE(*,*)'new blkp'
          WRITE(*,1323)new_blk_p
       ENDIF
    ENDIF
  END SUBROUTINE transpose_index_local



! *****************************************************************************
!> \brief Makes a canonical index to the distribution.
!> \par Canonical means that it respects the distribution.
!> \param[inout] matrix       matrix for which to make canonical index
! *****************************************************************************
  SUBROUTINE dbcsr_make_index_canonical (matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    INTEGER                                  :: nb, nc, nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: new_blk_p, new_col_i, &
                                                new_row_p

!   ---------------------------------------------------------------------------

    nr = SIZE(matrix%m%row_p)
    ALLOCATE (new_row_p (nr))
    nc = SIZE(matrix%m%col_i)
    ALLOCATE (new_col_i (nc))
    nb = SIZE(matrix%m%blk_p)
    ALLOCATE (new_blk_p (nb))
    CALL make_index_canonical (new_row_p, new_col_i, new_blk_p,&
         matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p, matrix%m)
    matrix%m%row_p(:) = new_row_p
    matrix%m%col_i(:) = new_col_i
    matrix%m%blk_p(:) = new_blk_p
  END SUBROUTINE dbcsr_make_index_canonical


! *****************************************************************************
!> \brief Makes a canonical index given the index arrays
!>
!> \note
!> This routine uses hard-coded logic as to what constitutes a
!> canonical ordering
!> \par Description of canonical ordering
!> A non-(anti)symmetric matrix is left as is. Otherwise, the row and column
!> are stored in the position prescribed by the distribution.
! *****************************************************************************
  SUBROUTINE make_index_canonical (new_row_p, new_col_i, new_blk_p,&
       old_row_p, old_col_i, old_blk_p, matrix)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: new_row_p, new_col_i, &
                                                new_blk_p
    INTEGER, DIMENSION(:), INTENT(IN)        :: old_row_p, old_col_i, &
                                                old_blk_p
    TYPE(dbcsr_type), INTENT(IN)             :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'make_index_canonical', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_p, col, nblks, row, &
                                                stored_col, stored_row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: row_i
    LOGICAL                                  :: tr

!   ---------------------------------------------------------------------------

    nblks = SIZE(old_blk_p)
    ALLOCATE(row_i (nblks))
    !write(*,*)"old row_p", old_row_p
    !write(*,*)"old col_i", old_col_i
    !write(*,*)"old blk_p", old_blk_p
    DO row = 1, SIZE(old_row_p)-1
       DO blk = old_row_p(row)+1, old_row_p(row+1)
          col = old_col_i(blk)
          stored_row = row
          stored_col = col
          tr = .FALSE.
          CALL get_stored_canonical (matrix, stored_row, stored_col, tr)
          !write(*,'(A,2(1X,I5),A,2(1X,I5),";",I7)')&
          !     routineN//" X->",row,col,"->",&
          !     stored_row, stored_col,blk
          row_i(blk) = stored_row
          new_col_i(blk) = stored_col
          IF (.NOT. tr) THEN
             new_blk_p(blk) = old_blk_p(blk)
          ELSE
             new_blk_p(blk) = -old_blk_p(blk)
          ENDIF
       ENDDO
    ENDDO
    CALL dbcsr_sort_indices(nblks, row_i, new_col_i, blk_p=new_blk_p)
    ! Re-create the index
    CALL dbcsr_make_dbcsr_index (new_row_p, row_i, SIZE(new_row_p)-1, nblks)
    !write(*,*)"new row_p", new_row_p
    !write(*,*)"new row_i", row_i
    !write(*,*)"new col_i", new_col_i
    !write(*,*)"new blk_p", new_blk_p
  END SUBROUTINE make_index_canonical


  SUBROUTINE dbcsr_make_dbcsr_index (row_p, row_i, nrows, nblks)
    INTEGER, INTENT(in)                      :: nrows, nblks
    INTEGER, DIMENSION(1:nrows+1), &
      INTENT(out)                            :: row_p
    INTEGER, DIMENSION(1:nblks), INTENT(in)  :: row_i

    INTEGER                                  :: blk, row

!
!

    row_p(1) = 0
    row_p(nrows+1) = nblks
    row = 1
    blk = 1
    DO WHILE (row .LE. nrows)
       IF (blk .LE. nblks) THEN
          DO WHILE (row_i(blk) .EQ. row)
             blk = blk+1
             IF (blk .GT. nblks) THEN
                row_p(row+1) = nblks-1
                EXIT
             ENDIF
          ENDDO
       ENDIF
       row_p(row+1) = blk-1
       row = row+1
    ENDDO
  END SUBROUTINE dbcsr_make_dbcsr_index



! *****************************************************************************
!> \brief Fills meta information from a given distribution_2d
!> \param[out] meta           meta information array to fill
!> \param dist                processor distribution
!> \param[in] row_blk_size    row block sizes
!> \param[in] col_blk_size    column block sizes
! *****************************************************************************
  SUBROUTINE meta_from_dist (meta, dist, row_blk_size, col_blk_size)
    INTEGER, DIMENSION(dbcsr_meta_size), &
      INTENT(OUT)                            :: meta
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_size, col_blk_size

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_from_dist', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nfullcols_local, &
                                                nfullcols_total, &
                                                nfullrows_local, &
                                                nfullrows_total
    INTEGER, DIMENSION(:), POINTER           :: blkcols_local, blkrows_local, &
                                                cbs, rbs

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    blkrows_local => array_data (dbcsr_distribution_local_rows (dist))
    blkcols_local => array_data (dbcsr_distribution_local_cols (dist))
    rbs => array_data (row_blk_size)
    cbs => array_data (col_blk_size)
    nfullrows_total = SUM (rbs)
    nfullcols_total = SUM (cbs)
    nfullrows_local = 0
    nfullcols_local = 0
    DO i = 1, dbcsr_distribution_nlocal_rows (dist)
       nfullrows_local = nfullrows_local + rbs(blkrows_local(i))
    ENDDO
    DO i = 1, dbcsr_distribution_nlocal_cols (dist)
       nfullcols_local = nfullcols_local + cbs(blkcols_local(i))
    ENDDO
    meta(:) = 0
    meta(3)  = dbcsr_distribution_nrows (dist)
    meta(4)  = dbcsr_distribution_ncols (dist)
    meta(5)  = nfullrows_total
    meta(6)  = nfullcols_total
    meta(7)  = dbcsr_distribution_nlocal_rows (dist)
    meta(8)  = dbcsr_distribution_nlocal_cols (dist)
    meta(9)  = nfullrows_local
    meta(10) = nfullcols_local
  END SUBROUTINE meta_from_dist


! *****************************************************************************
!> \brief Test routine for create_bl_distribution.
! *****************************************************************************
  SUBROUTINE test_create_bl_distribution()

    INTEGER, PARAMETER                       :: rpt = 10

    INTEGER                                  :: bins, elements, i
    REAL                                     :: rn
    TYPE(array_i1d_obj)                      :: block_distribution, block_size

    DO i = 1, rpt
       CALL RANDOM_NUMBER (rn)
       elements = INT(rn * 500)
       IF (i .EQ. 1) elements = 0
       CALL RANDOM_NUMBER (rn)
       bins = MAX(1, INT(rn*8))
       CALL create_bl_distribution (block_distribution,&
            block_size, elements, bins)
       CALL array_release (block_distribution)
       CALL array_release (block_size)
    ENDDO
  END SUBROUTINE test_create_bl_distribution

! *****************************************************************************
!> \brief Distributes elements into blocks and into bins
!>
!> \par Term clarification
!>      An example: blocks are atom blocks and bins are process rows/colums.
!> \param[out] block_distribution       block distribution to bins
!> \param[out] block_size       sizes of blocks
!> \param[in] elements          number of elements to bin
!> \param[in] nbins             number of bins
! *****************************************************************************
  SUBROUTINE create_bl_distribution (block_distribution,&
       block_size, nelements, nbins)
    TYPE(array_i1d_obj), INTENT(OUT)         :: block_distribution, block_size
    INTEGER, INTENT(IN)                      :: nelements, nbins

    CHARACTER(len=*), PARAMETER :: routineN = 'create_bl_distribution', &
      routineP = moduleN//':'//routineN

    INTEGER :: bin, blk_layer, element_stack, els, estimated_blocks, &
      max_block_size = 32, max_blocks_per_bin, max_elements_per_bin, &
      max_elements_per_block, nblks, stat
    INTEGER, DIMENSION(:), POINTER           :: blk_dist, blk_sizes

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL array_nullify (block_distribution)
    CALL array_nullify (block_size)
    ! Define the sizes on which we build the distribution.
    IF (nelements .GT. 0) THEN
       max_elements_per_bin = CEILING (REAL(nelements)/REAL(nbins))
       max_elements_per_block = MIN (max_elements_per_bin, max_block_size)
       max_blocks_per_bin = CEILING (REAL(max_elements_per_bin)/REAL(max_elements_per_block))


       IF (dbg) bin = max_elements_per_block
       max_elements_per_block = CEILING (REAL(max_elements_per_bin)/REAL(max_blocks_per_bin))

       IF (max_elements_per_block*max_blocks_per_bin*(nbins-1) .GE. nelements) THEN
          max_elements_per_block = MAX(1, max_elements_per_block-1)
          max_blocks_per_bin = max_blocks_per_bin+1
          DBG " small number FIX"
       ENDIF

       IF (dbg) THEN
          WRITE(*,'(1X,A,1X,A,I7,A,I7,A)')routineP,"For",nelements,&
               " elements and",nbins," bins"
          WRITE(*,'(1X,A,1X,A,I7,A)')routineP,"There are",&
               max_elements_per_bin," max elements/bin"
          WRITE(*,'(1X,A,1X,A,I7,A)')routineP,"There are",&
               bin," true max elements/block"
          WRITE(*,'(1X,A,1X,A,I7,A)')routineP,"There are",&
               max_blocks_per_bin," max blocks/bin"
          WRITE(*,'(1X,A,1X,A,I7,A)')routineP,"There are",&
               max_elements_per_block," max elements/block"
       ENDIF
       ! Prepare block distribution and sizes array, which may be
       ! oversized at this point.
       estimated_blocks = max_blocks_per_bin * nbins
       ALLOCATE (blk_dist (estimated_blocks), stat=stat)
       IF (stat /= 0) CALL stop_memory (routineP, "blk_dist", estimated_blocks)
       ALLOCATE (blk_sizes (estimated_blocks), stat=stat)
       IF (stat /= 0) CALL stop_memory (routineP, "blk_sizes", estimated_blocks)
       !
       element_stack = 0
       nblks = 0
       DO blk_layer = 1, max_blocks_per_bin
          DO bin = 0, nbins-1
             els = MIN (max_elements_per_block, nelements-element_stack)
             IF (els .GT. 0) THEN
                element_stack = element_stack + els
                nblks = nblks + 1
                blk_dist(nblks) = bin
                blk_sizes(nblks) = els
                IF (dbg) WRITE(*,'(1X,A,I5,A,I5,A,I5)')routineP//" Assigning",&
                     els," elements as block",nblks," to bin",bin
             ENDIF
          ENDDO
       ENDDO
       ! Create the output arrays.
       DBG "Wound up with",nblks,'blocks instead of',estimated_blocks
       IF (nblks .EQ. estimated_blocks) THEN
          CALL array_new (block_distribution, blk_dist, gift=.TRUE.)
          CALL array_new (block_size, blk_sizes, gift=.TRUE.)
       ELSE
          CALL array_new (block_distribution, blk_dist(1:nblks), lb=1)
          CALL array_new (block_size, blk_sizes(1:nblks), lb=1)
          DEALLOCATE (blk_dist)
          DEALLOCATE (blk_sizes)
       ENDIF
    ELSE
       ALLOCATE (blk_dist (0), stat=stat)
       IF (stat /= 0) CALL stop_memory (routineP, "blk_dist", estimated_blocks)
       ALLOCATE (blk_sizes (0), stat=stat)
       IF (stat /= 0) CALL stop_memory (routineP, "blk_sizes", estimated_blocks)
       CALL array_new (block_distribution, blk_dist, gift=.TRUE.)
       CALL array_new (block_size, blk_sizes, gift=.TRUE.)
    ENDIF
1579 FORMAT (I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5)
    IF (dbg) THEN
       WRITE(*,'(1X,A,A)')routineP//" Distribution"
       WRITE(*,1579)array_data (block_distribution)
       WRITE(*,'(1X,A,A)')routineP//" Sizes"
       WRITE(*,1579)array_data (block_size)
    ENDIF
  END SUBROUTINE create_bl_distribution


! *****************************************************************************
!> \brief Create the distribution_2d for the right matrix in a matrix-matrix
!>        multiplication, given the distribution_2d of the left matrix.
!> \par The new row distribution for the right matrix mirrors the column
!>      distribution of the left matrix, while the column distribution is
!>      created so that it is appropriate to the parallel environment. The
!>      columns are blocked together and each block is assigned to its own
!>      atom type.
!> \param[out] dist2d         new distribution_2d for the right matrix
!> \param[in] dist2d_left     the distribution_2d of the left matrix
!> \param[in] ncolumns        number of columns in right matrix
!> \param row_blk_sizes, col_blk_sizes  (optional) sizes of blocks for the
!>                                      row or column. Calculated if
!>                                      specified but not associated.
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist2d_r(dist2d, dist2d_left, ncolumns,&
       col_blk_size, error)
    TYPE(distribution_2d_type), POINTER      :: dist2d, dist2d_left
    INTEGER, INTENT(IN)                      :: ncolumns
    TYPE(array_i1d_obj), INTENT(OUT), &
      OPTIONAL                               :: col_blk_size
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist2d_r', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: blksize_default = 3

    INTEGER :: at, at_high, at_low, column_stack, i, k, mynod, mypcol, &
      mypcol_left, myprow, myprow_left, nblkcols, nblkcols_left, nblkrows, &
      nblks, ncols, npcols, npcols_left, nprows, nprows_left, pcol, pcsize, &
      stat
    INTEGER, DIMENSION(:), POINTER :: cbs, col_dist, col_dist_left, &
      n_at_local_cols, n_at_local_cols_left, n_at_local_rows, row_dist, &
      row_dist_left
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi, mpi2blacs
    LOGICAL                                  :: dbg
    TYPE(cp_1d_i_p_type), DIMENSION(:), &
      POINTER                                :: at_local_cols, &
                                                at_local_cols_left, &
                                                at_local_rows
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env, blacs_env_right
    TYPE(cp_para_env_type), POINTER          :: para_env

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    IF (dbg) CALL distribution_2d_write(dist2d_left, 6,&
         long_description=.TRUE., error=error, local=.TRUE.)
    IF (dbg) &
         CALL distribution_2d_write(dist2d_left, 6, long_description=.TRUE.,&
         error=error, local=.TRUE.)
    CALL distribution_2d_get(dist2d_left, error=error,&
         row_distribution=row_dist_left, col_distribution=col_dist_left,&
         local_cols=at_local_cols_left, n_local_cols=n_at_local_cols_left,&
         n_col_distribution=nblkcols_left, blacs_env=blacs_env)
    CALL get_blacs_info(blacs_env, my_process_row=myprow_left,&
         my_process_column=mypcol_left, number_of_process_rows=nprows_left,&
         number_of_process_columns=npcols_left, para_env=para_env,&
         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
    !  myprow = mypcol_left
    !  mypcol = myprow_left
    !  nprows = npcols_left
    !  npcols = nprows_left

    !  myprow = myprow_left
    !  mypcol = mypcol_left
    !  nprows = nprows_left
    !  npcols = npcols_left

    !
    !> Create the transposed blacs_env
    IF (dbg) WRITE(*,*)routineP//' LEFT prows, pcols:', nprows_left, npcols_left
    CALL cp_transpose_blacs_env(blacs_env, blacs_env_right, error)
    CALL get_blacs_info(blacs_env_right, my_process_row=myprow,&
         my_process_column=mypcol, number_of_process_rows=nprows,&
         number_of_process_columns=npcols,&
         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
    IF (dbg) WRITE(*,*)routineP//' RIGHT prows, pcols:', nprows, npcols
    IF (dbg) WRITE(*,*)routineP//' myprow, mypcol',myprow,'/',nprows,';',&
         mypcol,'/',npcols
    IF (dbg) CALL cp_blacs_env_write(blacs_env_right, 6, error)
    !
    mynod = para_env%mepos
    !> Handle row distribution. This is based on the column distribution of the
    !> left matrix.
    ALLOCATE(at_local_rows(LBOUND(at_local_cols_left,1):&
         UBOUND(at_local_cols_left,1)), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'at_local_rows',&
         SIZE(at_local_cols_left,1))
    at_low = LBOUND(n_at_local_cols_left,1) 
    at_high = UBOUND(n_at_local_cols_left,1)
    ALLOCATE(n_at_local_rows(at_low:at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'n_at_local_rows',&
         SIZE(n_at_local_cols_left,1))
    nblkrows = nblkcols_left
    ALLOCATE(row_dist(nblkrows), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'row_dist', nblkrows)
    n_at_local_rows(:) = n_at_local_cols_left(:)
    !> The rows and columns of the right matrix are the transposes of the left,
    !> so we can safely use the right-matrix's column data for our row data.
    !  IF(nprows.EQ.npcols) THEN
    !at_local_rows%array(:) = at_local_cols_left%array(:)
    DO i = LBOUND(at_local_cols_left,1), UBOUND(at_local_cols_left,1)
       ALLOCATE(at_local_rows(i)%array(&
            LBOUND(at_local_cols_left(i)%array,1):&
            UBOUND(at_local_cols_left(i)%array,1)), stat=stat)
       IF(stat /= 0) CALL stop_memory(routineP, 'at_local_rows(?)%array',&
            SIZE(at_local_cols_left(i)%array))
       at_local_rows(i)%array(:) = at_local_cols_left(i)%array(:)
    ENDDO
    n_at_local_rows(:) = n_at_local_cols_left(:)
    row_dist(:) = col_dist_left(:)
    ! Handle column distribution.
    ! Estimate the number of full columns per processor column, but it
    ! must be at least 1.
    pcsize = CEILING((1.0*ncolumns)/(1.0*npcols)) ! Full columns per processor column.
    column_stack = ncolumns
    nblkcols = 0
    DO pcol = 0, npcols-1
       at_low = nblkcols+1
       at_high = nblkcols
       IF (column_stack .GT. 0) THEN
          ncols = MIN(column_stack, pcsize)
          column_stack = column_stack - ncols
          nblks = CEILING ((1.0*ncols)/(1.0*blksize_default))
          nblkcols = nblkcols + nblks
          at_high = nblkcols
       ENDIF
       IF (dbg) WRITE(*,*)routineP//' pcol',pcol,'has blkcols',at_low,':',&
            at_high,'col stack is',column_stack
    ENDDO

    IF (dbg) WRITE(*,*)routineP//' There are',nblkcols,'blocked columns in',&
         npcols,'pcols for',ncolumns,'full columns.'

    at_low = 1
    at_high = nblkcols
    ALLOCATE(col_dist(nblkcols), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'col_dist', nblkcols)
    col_dist(:) = -999
    ALLOCATE(n_at_local_cols(at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'n_at_local_cols', at_high)
    ALLOCATE(at_local_cols(at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols', at_high)
    n_at_local_cols(:) = 0

    IF (PRESENT(col_blk_size)) THEN
       ALLOCATE(cbs(nblkcols), stat=stat)
       IF (stat /= 0) CALL stop_memory(routineP,&
            'cbs', nblkcols)
    ENDIF
    column_stack = ncolumns
    nblkcols = 0
    i = 1
    DO pcol = 0, npcols-1
       at_low = nblkcols+1
       at_high = nblkcols
       nblks = 0
       IF (column_stack .GT. 0) THEN
          ncols = MIN(column_stack, pcsize)
          column_stack = column_stack - ncols
          nblks = CEILING ((1.0*ncols)/(1.0*blksize_default))
          nblkcols = nblkcols + nblks
          at_high = nblkcols
       ENDIF
       IF (dbg) WRITE(*,*)routineP//' pcol',pcol,'has blkcols',at_low,':',&
            at_high,'col stack is',column_stack
       col_dist(at_low:at_high) = pcol
       DO at = at_low, at_high
          IF (pcol.EQ.mypcol) THEN
             ALLOCATE(at_local_cols(at)%array(1), stat=stat)
             IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols(?)', 1)
             at_local_cols(at)%array(1) = at
             n_at_local_cols(at) = 1
          ELSE
             ALLOCATE(at_local_cols(at)%array(0), stat=stat)
             IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols(?)', 1)
          ENDIF
          IF (PRESENT(col_blk_size)) THEN
             k = MIN(ncols,blksize_default)
             cbs(at) = k
             i = i + k
             ncols = ncols - k
             IF (dbg) WRITE(*,*)routineP//' col',at,'has',k-1,&
                  'full colums; ncols, i',ncols, i,&
                  'nblks',nblks
          ENDIF
       ENDDO
    ENDDO

    IF(dbg) WRITE(*,*)routineP//' Col. distribution',col_dist
    IF (PRESENT (col_blk_size)) THEN
       IF(dbg) WRITE(*,*)routineP//' Col. block sizes',cbs
       CALL array_new (col_blk_size, cbs, gift=.TRUE.)
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Col. at_local_cols=',n_at_local_cols
    IF(dbg) WRITE(*,*)routineP//' Col. n_at_local_cols=',n_at_local_cols

    NULLIFY(dist2d)
    CALL distribution_2d_create(dist2d, blacs_env_right, error=error,&
         local_rows_ptr=at_local_rows, n_local_rows=n_at_local_rows,&
         local_cols_ptr=at_local_cols, n_local_cols=n_at_local_cols,&
         row_distribution_ptr=row_dist, col_distribution_ptr=col_dist)
    DEALLOCATE(n_at_local_rows)
    DEALLOCATE(n_at_local_cols)
    CALL cp_blacs_env_release(blacs_env_right, error)
    IF (dbg) CALL distribution_2d_write(dist2d, 6, long_description=.TRUE.,&
        error=error, local=.TRUE.)
    IF (dbg) &
         CALL distribution_2d_write(dist2d, 6, long_description=.TRUE.,&
         error=error, local=.TRUE.)
  END SUBROUTINE dbcsr_create_dist2d_r

! *****************************************************************************
!> \brief Create a transposed distribution_2d
!> \param[out] dist2d         transposed distribution_2d
!> \param[in] dist2d_normal   the base distribution_2d
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist2d_transposed(dist2d, dist2d_normal,  error)
    TYPE(distribution_2d_type), POINTER      :: dist2d, dist2d_normal
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_create_dist2d_transposed', &
      routineP = moduleN//':'//routineN

    INTEGER :: at_high, at_low, i, mynod, mypcol_normal, mypcol_tr, &
      myprow_normal, myprow_tr, nblkcols_normal, nblkcols_tr, &
      nblkrows_normal, nblkrows_tr, npcols_normal, npcols_tr, nprows_normal, &
      nprows_tr, stat
    INTEGER, DIMENSION(:), POINTER :: col_dist_normal, col_dist_tr, &
      n_at_local_cols_normal, n_at_local_cols_tr, n_at_local_rows_normal, &
      n_at_local_rows_tr, row_dist_normal, row_dist_tr
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi, mpi2blacs
    TYPE(cp_1d_i_p_type), DIMENSION(:), &
      POINTER                                :: at_local_cols_normal, &
                                                at_local_cols_tr, &
                                                at_local_rows_normal, &
                                                at_local_rows_tr
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env_normal, blacs_env_tr
    TYPE(cp_para_env_type), POINTER          :: para_env

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    IF (dbg) CALL distribution_2d_write(dist2d_normal, 6,&
         long_description=.TRUE.,&
         error=error, local=.TRUE.)
    CALL dbcsr_set_debug(dbg)
    CALL distribution_2d_get(dist2d_normal, error=error,&
         row_distribution=row_dist_normal, col_distribution=col_dist_normal,&
         local_rows=at_local_rows_normal, n_local_rows=n_at_local_rows_normal,&
         local_cols=at_local_cols_normal, n_local_cols=n_at_local_cols_normal,&
         n_row_distribution=nblkrows_normal,&
         n_col_distribution=nblkcols_normal,&
         blacs_env=blacs_env_normal)
    CALL get_blacs_info(blacs_env_normal, my_process_row=myprow_normal,&
         my_process_column=mypcol_normal, number_of_process_rows=nprows_normal,&
         number_of_process_columns=npcols_normal, para_env=para_env,&
         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
    ! Create the transposed blacs_env
    CALL cp_transpose_blacs_env(blacs_env_normal, blacs_env_tr, error)
    CALL get_blacs_info(blacs_env_tr, my_process_row=myprow_tr,&
         my_process_column=mypcol_tr, number_of_process_rows=nprows_tr,&
         number_of_process_columns=npcols_tr,&
         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
    mynod = para_env%mepos
    !> Handle row distribution. This is based on the column distribution of the
    !> given matrix.
    at_low = LBOUND(n_at_local_cols_normal,1) 
    at_high = UBOUND(n_at_local_cols_normal,1)
    ALLOCATE(at_local_rows_tr(at_low:at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'at_local_rows_tr',&
         at_high-at_low+1)
    ALLOCATE(n_at_local_rows_tr(at_low:at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'n_at_local_rows_tr',&
         at_high-at_low+1)
    nblkrows_tr = nblkcols_normal
    ALLOCATE(row_dist_tr(nblkrows_tr), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'row_dist_tr', nblkrows_tr)
    n_at_local_rows_tr(:) = n_at_local_cols_normal(:)
    DO i = at_low, at_high
       ALLOCATE(at_local_rows_tr(i)%array(&
            LBOUND(at_local_cols_normal(i)%array,1):&
            UBOUND(at_local_cols_normal(i)%array,1)), stat=stat)
       IF(stat /= 0) CALL stop_memory(routineP, 'at_local_rows_tr(?)%array',&
            SIZE(at_local_cols_normal(i)%array))
       at_local_rows_tr(i)%array(:) = at_local_cols_normal(i)%array(:)
    END DO
    n_at_local_rows_tr(:) = n_at_local_cols_normal(:)
    row_dist_tr(:) = col_dist_normal(:)
    !
    ! Handle column distribution.
    at_low = LBOUND(n_at_local_rows_normal,1) 
    at_high = UBOUND(n_at_local_rows_normal,1)
    ALLOCATE(at_local_cols_tr(at_low:at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols_tr',&
         at_high-at_low+1)
    ALLOCATE(n_at_local_cols_tr(at_low:at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'n_at_local_cols_tr',&
         at_high-at_low+1)
    nblkcols_tr = nblkrows_normal
    ALLOCATE(col_dist_tr(nblkcols_tr), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'col_dist_tr', nblkcols_tr)
    n_at_local_cols_tr(:) = n_at_local_rows_normal(:)
    DO i = at_low, at_high
       ALLOCATE(at_local_cols_tr(i)%array(&
            LBOUND(at_local_rows_normal(i)%array,1):&
            UBOUND(at_local_rows_normal(i)%array,1)), stat=stat)
       IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols_tr(?)%array',&
            SIZE(at_local_rows_normal(i)%array))
       at_local_cols_tr(i)%array(:) = at_local_rows_normal(i)%array(:)
    END DO
    n_at_local_cols_tr(:) = n_at_local_rows_normal(:)
    col_dist_tr(:) = row_dist_normal(:)
    !
    NULLIFY(dist2d)
    CALL distribution_2d_create(dist2d, blacs_env_tr, error=error,&
         local_rows_ptr=at_local_rows_tr, n_local_rows=n_at_local_rows_tr,&
         local_cols_ptr=at_local_cols_tr, n_local_cols=n_at_local_cols_tr,&
         row_distribution_ptr=row_dist_tr, col_distribution_ptr=col_dist_tr)
    DEALLOCATE(n_at_local_rows_tr)
    DEALLOCATE(n_at_local_cols_tr)
    IF (dbg) CALL distribution_2d_write(dist2d, 6, long_description=.TRUE.,&
        error=error, local=.TRUE.)
  END SUBROUTINE dbcsr_create_dist2d_transposed



! *****************************************************************************
!> \brief Create the distribution_2d for the right matrix in a matrix-matrix
!>        multiplication, given the distribution_2d of the left matrix.
!> \par The new row distribution for the right matrix is the same as the row
!>      distribution of the left matrix, while the column distribution is
!>      created so that it is appropriate to the parallel environment. The
!>      columns are blocked together and each block is assigned to its own
!>      atom type.
!> \param[out] dist2d         new distribution_2d for the right matrix
!> \param[in] dist2d_left     the distribution_2d of the left matrix
!> \param[in] ncolumns        number of columns in right matrix
!> \param row_blk_sizes, col_blk_sizes  (optional) sizes of blocks for the
!>                                      row or column. Calculated if
!>                                      specified but not associated.
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist2d_r_unrot(dist2d, dist2d_left, ncolumns,&
       col_blk_size, error)
    TYPE(distribution_2d_type), POINTER      :: dist2d, dist2d_left
    INTEGER, INTENT(IN)                      :: ncolumns
    TYPE(array_i1d_obj), INTENT(OUT), &
      OPTIONAL                               :: col_blk_size
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist2d_r_unrot', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: blksize_default = 3

    INTEGER :: at, at_high, at_low, column_stack, i, k, mynod, mypcol, &
      mypcol_left, myprow, myprow_left, nblkcols, nblkrows, nblkrows_left, &
      nblks, ncols, npcols, npcols_left, nprows, nprows_left, pcol, pcsize, &
      stat
    INTEGER, DIMENSION(:), POINTER :: cbs, col_dist, col_dist_left, &
      n_at_local_cols, n_at_local_rows, n_at_local_rows_left, row_dist, &
      row_dist_left
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi, mpi2blacs
    TYPE(cp_1d_i_p_type), DIMENSION(:), &
      POINTER                                :: at_local_cols, at_local_rows, &
                                                at_local_rows_left
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env, blacs_env_right
    TYPE(cp_para_env_type), POINTER          :: para_env

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    IF (dbg) CALL distribution_2d_write(dist2d_left, 6,&
         long_description=.TRUE., error=error, local=.TRUE.)
    IF (dbg) &
         CALL distribution_2d_write(dist2d_left, 6, long_description=.TRUE.,&
         error=error, local=.TRUE.)
    CALL distribution_2d_get(dist2d_left, error=error,&
         row_distribution=row_dist_left, col_distribution=col_dist_left,&
         local_rows=at_local_rows_left, n_local_rows=n_at_local_rows_left,&
         n_row_distribution=nblkrows_left, blacs_env=blacs_env)
    CALL get_blacs_info(blacs_env, my_process_row=myprow_left,&
         my_process_column=mypcol_left, number_of_process_rows=nprows_left,&
         number_of_process_columns=npcols_left, para_env=para_env,&
         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
    !
    IF (dbg) WRITE(*,*)routineP//' LEFT prows, pcols:', nprows_left, npcols_left
    blacs_env_right => blacs_env
    !CALL cp_blacs_env_retain(blacs_env, error)
    CALL get_blacs_info(blacs_env_right, my_process_row=myprow,&
         my_process_column=mypcol, number_of_process_rows=nprows,&
         number_of_process_columns=npcols,&
         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
    IF (dbg) WRITE(*,*)routineP//' RIGHT prows, pcols:', nprows, npcols
    IF (dbg) WRITE(*,*)routineP//' myprow, mypcol',myprow,'/',nprows,';',&
         mypcol,'/',npcols
    IF (dbg) CALL cp_blacs_env_write(blacs_env_right, 6, error)
    !
    mynod = para_env%mepos
    !> Handle row distribution. This is based on the column distribution of the
    !> left matrix.
    ALLOCATE(at_local_rows(LBOUND(at_local_rows_left,1):&
         UBOUND(at_local_rows_left,1)), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'at_local_rows',&
         SIZE(at_local_rows_left,1))
    at_low = LBOUND(n_at_local_rows_left,1) 
    at_high = UBOUND(n_at_local_rows_left,1)
    ALLOCATE(n_at_local_rows(at_low:at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'n_at_local_rows',&
         SIZE(n_at_local_rows_left,1))
    nblkrows = nblkrows_left
    ALLOCATE(row_dist(nblkrows), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'row_dist', nblkrows)
    n_at_local_rows(:) = n_at_local_rows_left(:)
    DO i = LBOUND(at_local_rows_left,1), UBOUND(at_local_rows_left,1)
       ALLOCATE(at_local_rows(i)%array(&
            LBOUND(at_local_rows_left(i)%array,1):&
            UBOUND(at_local_rows_left(i)%array,1)), stat=stat)
       IF(stat /= 0) CALL stop_memory(routineP, 'at_local_rows(?)%array',&
            SIZE(at_local_rows_left(i)%array))
       at_local_rows(i)%array(:) = at_local_rows_left(i)%array(:)
    ENDDO
    n_at_local_rows(:) = n_at_local_rows_left(:)
    row_dist(:) = row_dist_left(:)
    ! Handle column distribution.
    ! Estimate the number of full columns per processor column, but it
    ! must be at least 1.
    pcsize = CEILING((1.0*ncolumns)/(1.0*npcols)) ! Full columns per processor column.
    column_stack = ncolumns
    nblkcols = 0
    DO pcol = 0, npcols-1
       at_low = nblkcols+1
       at_high = nblkcols
       IF (column_stack .GT. 0) THEN
          ncols = MIN(column_stack, pcsize)
          column_stack = column_stack - ncols
          nblks = CEILING ((1.0*ncols)/(1.0*blksize_default))
          nblkcols = nblkcols + nblks
          at_high = nblkcols
       ENDIF
       IF (dbg) WRITE(*,*)routineP//' pcol',pcol,'has blkcols',at_low,':',&
            at_high,'col stack is',column_stack
    ENDDO

    IF (dbg) WRITE(*,*)routineP//' There are',nblkcols,'blocked columns in',&
         npcols,'pcols for',ncolumns,'full columns.'

    at_low = 1
    at_high = nblkcols
    ALLOCATE(col_dist(nblkcols), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'col_dist', nblkcols)
    col_dist(:) = -999
    ALLOCATE(n_at_local_cols(at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'n_at_local_cols', at_high)
    ALLOCATE(at_local_cols(at_high), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols', at_high)
    n_at_local_cols(:) = 0

    IF (PRESENT(col_blk_size)) THEN
       ALLOCATE(cbs(nblkcols), stat=stat)
       IF (stat /= 0) CALL stop_memory(routineP,&
            'cbs', nblkcols)
    ENDIF
    column_stack = ncolumns
    nblkcols = 0
    i = 1
    DO pcol = 0, npcols-1
       at_low = nblkcols+1
       at_high = nblkcols
       nblks = 0
       IF (column_stack .GT. 0) THEN
          ncols = MIN(column_stack, pcsize)
          column_stack = column_stack - ncols
          nblks = CEILING ((1.0*ncols)/(1.0*blksize_default))
          nblkcols = nblkcols + nblks
          at_high = nblkcols
       ENDIF
       IF (dbg) WRITE(*,*)routineP//' pcol',pcol,'has blkcols',at_low,':',&
            at_high,'col stack is',column_stack
       col_dist(at_low:at_high) = pcol
       DO at = at_low, at_high
          IF (pcol.EQ.mypcol) THEN
             ALLOCATE(at_local_cols(at)%array(1), stat=stat)
             IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols(?)', 1)
             at_local_cols(at)%array(1) = at
             n_at_local_cols(at) = 1
          ELSE
             ALLOCATE(at_local_cols(at)%array(0), stat=stat)
             IF(stat /= 0) CALL stop_memory(routineP, 'at_local_cols(?)', 1)
          ENDIF
          IF (PRESENT(col_blk_size)) THEN
             k = MIN(ncols,blksize_default)
             cbs(at) = k
             i = i + k
             ncols = ncols - k
             IF (dbg) WRITE(*,*)routineP//' col',at,'has',k-1,&
                  'full colums; ncols, i',ncols, i,&
                  'nblks',nblks
          ENDIF
       ENDDO
    ENDDO

    IF(dbg) WRITE(*,*)routineP//' Col. distribution',col_dist
    IF (PRESENT (col_blk_size)) THEN
       IF(dbg) WRITE(*,*)routineP//' Col. block sizes',cbs
       CALL array_new (col_blk_size, cbs, gift=.TRUE.)
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Col. at_local_cols=',n_at_local_cols
    IF(dbg) WRITE(*,*)routineP//' Col. n_at_local_cols=',n_at_local_cols

    NULLIFY(dist2d)
    CALL distribution_2d_create(dist2d, blacs_env_right, error=error,&
         local_rows_ptr=at_local_rows, n_local_rows=n_at_local_rows,&
         local_cols_ptr=at_local_cols, n_local_cols=n_at_local_cols,&
         row_distribution_ptr=row_dist, col_distribution_ptr=col_dist)
    DEALLOCATE(n_at_local_rows)
    DEALLOCATE(n_at_local_cols)
    !CALL cp_blacs_env_release(blacs_env_right, error)
    IF (dbg) CALL distribution_2d_write(dist2d, 6, long_description=.TRUE.,&
        error=error, local=.TRUE.)
    IF (dbg) &
         CALL distribution_2d_write(dist2d, 6, long_description=.TRUE.,&
         error=error, local=.TRUE.)
  END SUBROUTINE dbcsr_create_dist2d_r_unrot



! *****************************************************************************
!> \brief Creates a new distribution for the right matrix in a matrix
!>        multiplication with unrotated grid.
!> \par The new row distribution for the right matrix is the same as the row
!>      distribution of the left matrix, while the column distribution is
!>      created so that it is appropriate to the parallel environment.
!> \param[out] dist_right     new distribution for the right matrix
!> \param[in] dist_left       the distribution of the left matrix
!> \param[in] ncolumns        number of columns in right matrix
!> \param[out] right_col_blk_sizes      sizes of blocks in the created column
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_r_unrot(dist_right, dist_left, ncolumns,&
       right_col_blk_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_right
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_left
    INTEGER, INTENT(IN)                      :: ncolumns
    TYPE(array_i1d_obj), INTENT(OUT)         :: right_col_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_r_unrot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: multiplicity, new_npcols, &
                                                new_nprows, nimages
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_images
    INTEGER, DIMENSION(:), POINTER           :: right_rd
    TYPE(array_i1d_obj)                      :: left_col_dist, &
                                                right_col_dist, right_row_dist
    TYPE(dbcsr_mp_obj)                       :: right_mp_env

!   ---------------------------------------------------------------------------

    right_mp_env = dbcsr_distribution_mp (dist_left)
    ! Create the column distribution
    CALL create_bl_distribution (right_col_dist,&
         right_col_blk_sizes, ncolumns, dbcsr_mp_npcols (right_mp_env))
    ! Create an even row distribution.
    left_col_dist = dbcsr_distribution_col_dist (dist_left)
    CALL array_nullify (right_row_dist)
    ALLOCATE (right_rd(array_size (left_col_dist)))
    ALLOCATE (tmp_images(array_size (left_col_dist)))
    new_nprows = dbcsr_mp_nprows (right_mp_env)
    new_npcols = dbcsr_mp_npcols (right_mp_env)
    nimages = lcm (new_nprows, new_npcols) / new_nprows
    multiplicity = new_nprows / gcd (new_nprows, new_npcols)
    CALL rebin_distribution (right_rd, tmp_images, array_data (left_col_dist),&
         new_nprows, multiplicity, nimages)
    CALL array_new (right_row_dist, right_rd, gift=.TRUE.)
    CALL dbcsr_distribution_new (dist_right, right_mp_env,&
         right_row_dist, right_col_dist)
    CALL array_release (right_col_dist)
    CALL array_release (right_row_dist)
  END SUBROUTINE dbcsr_create_dist_r_unrot

! *****************************************************************************
!> \brief Creates a new distribution for the right matrix in a matrix
!>        multiplication with rotated grid.
!> \par The new row distribution for the right matrix is the same as the
!>      column distribution of the left matrix, while the column distribution
!>      is created so that it is appropriate to the parallel environment.
!> \param[out] dist_right     new distribution for the right matrix
!> \param[in] dist_left       the distribution of the left matrix
!> \param[in] ncolumns        number of columns in right matrix
!> \param[out] right_col_blk_sizes      sizes of blocks in the created column
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_r_rot(dist_right, dist_left, ncolumns,&
       right_col_blk_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_right
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_left
    INTEGER, INTENT(IN)                      :: ncolumns
    TYPE(array_i1d_obj), INTENT(OUT)         :: right_col_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_r_rot', &
      routineP = moduleN//':'//routineN

    TYPE(array_i1d_obj)                      :: right_col_dist
    TYPE(dbcsr_mp_obj)                       :: right_mp_env

!   ---------------------------------------------------------------------------
!CALL test_create_bl_distribution (error)

    CALL dbcsr_mp_new_transposed (right_mp_env,&
         dbcsr_distribution_mp (dist_left))
    CALL create_bl_distribution (right_col_dist,&
         right_col_blk_sizes, ncolumns, dbcsr_mp_npcols (right_mp_env))
    CALL dbcsr_distribution_new (dist_right, right_mp_env,&
         dbcsr_distribution_col_dist (dist_left), right_col_dist,&
         local_rows=dbcsr_distribution_local_cols (dist_left))
    CALL dbcsr_mp_release (right_mp_env)
    CALL array_release (right_col_dist)
  END SUBROUTINE dbcsr_create_dist_r_rot

! *****************************************************************************
!> \brief Creates a new distribution for the left matrix in a matrix
!>        multiplication with unrotated grid.
!> \par The new column distribution for the left matrix is the same as the
!>      column distribution of the right matrix, while the row distribution is
!>      created so that it is appropriate to the parallel environment.
!> \param[out] dist_left      new distribution for the left matrix
!> \param[in] dist_right      old distribution of the right matrix
!> \param[in] nrows           number of rows in the left matrix
!> \param[out] left_row_blk_sizes       sizes of blocks in the created rows
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_l_unrot(dist_left, dist_right, nrows,&
       left_row_blk_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_left
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_right
    INTEGER, INTENT(IN)                      :: nrows
    TYPE(array_i1d_obj), INTENT(OUT)         :: left_row_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_l_unrot', &
      routineP = moduleN//':'//routineN

    TYPE(array_i1d_obj)                      :: left_row_dist
    TYPE(dbcsr_mp_obj)                       :: left_mp_env

!   ---------------------------------------------------------------------------

    left_mp_env = dbcsr_distribution_mp (dist_left)
    CALL create_bl_distribution (left_row_dist,&
         left_row_blk_sizes, nrows, dbcsr_mp_nprows (left_mp_env))
    CALL dbcsr_distribution_new (dist_left, left_mp_env,&
         left_row_dist, dbcsr_distribution_col_dist (dist_right),&
         local_cols=dbcsr_distribution_local_cols (dist_right))
    CALL array_release (left_row_dist)
  END SUBROUTINE dbcsr_create_dist_l_unrot

! *****************************************************************************
!> \brief Creates a new distribution for the left matrix in a matrix
!>        multiplication with rotated grid.
!> \par The new column distribution for the left matrix is the same as the
!>      row distribution of the right matrix, while the row distribution is
!>      created so that it is appropriate to the parallel environment.
!> \param[out] dist_left      new distribution for the left matrix
!> \param[in] dist_right      the distribution of the right matrix
!> \param[in] nrows           (optional) number of rows in left matrix
!> \param[out] left_row_blk_sizes      sizes of blocks in the created rows
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_l_rot(dist_right, dist_left, nrows,&
       left_row_blk_sizes)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_right
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_left
    INTEGER, INTENT(IN)                      :: nrows
    TYPE(array_i1d_obj), INTENT(OUT)         :: left_row_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_l_rot', &
      routineP = moduleN//':'//routineN

    TYPE(array_i1d_obj)                      :: left_row_dist
    TYPE(dbcsr_mp_obj)                       :: left_mp_env

!   ---------------------------------------------------------------------------

    CALL dbcsr_mp_new_transposed (left_mp_env,&
         dbcsr_distribution_mp (dist_right))
    CALL create_bl_distribution (left_row_dist,&
         left_row_blk_sizes, nrows, dbcsr_mp_nprows (left_mp_env))
    CALL dbcsr_distribution_new (dist_left, left_mp_env,&
         left_row_dist, dbcsr_distribution_row_dist (dist_right),&
         local_cols=dbcsr_distribution_local_rows (dist_right))
    CALL dbcsr_mp_release (left_mp_env)
    CALL array_release (left_row_dist)
  END SUBROUTINE dbcsr_create_dist_l_rot

! *****************************************************************************
!> \brief Creates a distribution where all block sizes are 1.
!> \param[out] dist           the elemental distribution
!> \param[in] nrows           number of full rows
!> \param[in] ncolumns        number of full columns
!> \param[in] mp_env          mp environment
!> \param[out] row_blk_sizes  (optional) row block sizes (all 1)
!> \param[out] col_blk_sizes  (optional) column block sizes (all 1)
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_elemental(dist, nrows, ncolumns,&
       mp_env, row_blk_sizes, col_blk_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    INTEGER, INTENT(IN)                      :: nrows, ncolumns
    TYPE(dbcsr_mp_obj)                       :: mp_env
    TYPE(array_i1d_obj), INTENT(OUT), &
      OPTIONAL                               :: row_blk_sizes, col_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_elemental', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: npdims, pdim, sz
    INTEGER, DIMENSION(:), POINTER           :: cbs_data, cd_data, rbs_data, &
                                                rd_data
    TYPE(array_i1d_obj)                      :: cbs, col_dist, rbs, row_dist

!   ---------------------------------------------------------------------------

    ALLOCATE (rbs_data(nrows), cbs_data(ncolumns))
    rbs_data(:) = 1
    cbs_data(:) = 1
    CALL array_new (rbs, rbs_data, gift=.TRUE.)
    CALL array_new (cbs, cbs_data, gift=.TRUE.)
    ALLOCATE (rd_data(nrows), cd_data(ncolumns))
    ! Calclate process row distribution
    npdims = dbcsr_mp_nprows (mp_env)
    sz = CEILING (REAL(nrows) / REAL(npdims))
    FORALL (pdim = 0 : npdims-1)
       rd_data(pdim*sz+1 : MIN((pdim+1)*sz, nrows)) = pdim
    END FORALL
    CALL array_new (row_dist, rd_data, gift=.TRUE.)
    ! Calclate process column distribution
    npdims = dbcsr_mp_npcols (mp_env)
    sz = CEILING (REAL(ncolumns) / REAL(npdims))
    FORALL (pdim = 0 : npdims-1)
       cd_data(pdim*sz+1 : MIN((pdim+1)*sz, ncolumns)) = pdim
    END FORALL
    CALL array_new (col_dist, cd_data, gift=.TRUE.)
    !
    CALL dbcsr_distribution_new (dist, mp_env,&
         row_dist, col_dist)
    CALL array_release (row_dist)
    CALL array_release (col_dist)
    IF (PRESENT (row_blk_sizes)) THEN
       row_blk_sizes = rbs
    ELSE
       CALL array_release (rbs)
    ENDIF
    IF (PRESENT (col_blk_sizes)) THEN
       col_blk_sizes = cbs
    ELSE
       CALL array_release (cbs)
    ENDIF
  END SUBROUTINE dbcsr_create_dist_elemental


! *****************************************************************************
!> \brief Creates a distribution where all block sizes are as large as
!>        possible for the given distribution.
!> \param[out] dist           the elemental distribution
!> \param[in] nrows           number of full rows
!> \param[in] ncolumns        number of full columns
!> \param[in] mp_env          mp environment
!> \param[out] row_blk_sizes  row block sizes
!> \param[out] col_blk_sizes  column block sizes
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_fullest(dist, nrows, ncolumns,&
       mp_env, row_blk_sizes, col_blk_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    INTEGER, INTENT(IN)                      :: nrows, ncolumns
    TYPE(dbcsr_mp_obj)                       :: mp_env
    TYPE(array_i1d_obj), INTENT(OUT)         :: row_blk_sizes, col_blk_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_dist_fullest', &
      routineP = moduleN//':'//routineN

    TYPE(array_i1d_obj)                      :: col_dist, row_dist

!   ---------------------------------------------------------------------------

    CALL create_bl_distribution(row_dist, row_blk_sizes,&
         nrows, dbcsr_mp_nprows (mp_env))
    CALL create_bl_distribution(col_dist, col_blk_sizes,&
         ncolumns, dbcsr_mp_npcols (mp_env))
    !
    CALL dbcsr_distribution_new (dist, mp_env,&
         row_dist, col_dist)
    !
    CALL array_release (row_dist)
    CALL array_release (col_dist)
  END SUBROUTINE dbcsr_create_dist_fullest

! *****************************************************************************
!> \brief Creates a block-cyclic compatible distribution
!>
!>        All blocks in a dimension, except for possibly the last
!>        block, have the same size.
!> \param[out] dist           the elemental distribution
!> \param[in] nrows           number of full rows
!> \param[in] ncolumns        number of full columns
!> \param[in] nrow_block      size of row blocks
!> \param[in] ncol_block      size of column blocks
!> \param[in] mp_env          multiprocess environment
!> \param[out] row_blk_sizes  row block sizes
!> \param[out] col_blk_sizes  column block sizes
! *****************************************************************************
  SUBROUTINE dbcsr_create_dist_block_cyclic (dist, nrows, ncolumns,&
       nrow_block, ncol_block, mp_env, row_blk_sizes, col_blk_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    INTEGER, INTENT(IN)                      :: nrows, ncolumns, nrow_block, &
                                                ncol_block
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    TYPE(array_i1d_obj), INTENT(OUT)         :: row_blk_sizes, col_blk_sizes

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_create_dist_block_cyclic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nblkcols, nblkrows, npcols, &
                                                nprows, pdim, sz
    INTEGER, DIMENSION(:), POINTER           :: cbs_data, cd_data, rbs_data, &
                                                rd_data
    TYPE(array_i1d_obj)                      :: col_dist, row_dist

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY

    ! Row sizes
    nblkrows = nrows / nrow_block
    sz = MOD(nrows, nrow_block)
    IF (sz .GT. 0) nblkrows = nblkrows + 1
    ALLOCATE (rbs_data(nblkrows), rd_data(nblkrows))
    rbs_data(1:nblkrows) = nrow_block
    IF (sz .NE. 0) rbs_data(nblkrows) = sz
    ! Column sizes
    nblkcols = ncolumns / ncol_block
    sz = MOD(ncolumns, ncol_block)
    IF (sz .GT. 0) nblkcols = nblkcols + 1
    ALLOCATE (cbs_data(nblkcols), cd_data(nblkcols))
    cbs_data(1:nblkcols) = ncol_block
    IF (sz .NE. 0) cbs_data(nblkcols) = sz
    CALL array_new (row_blk_sizes, rbs_data, gift=.TRUE.)
    CALL array_new (col_blk_sizes, cbs_data, gift=.TRUE.)
    !
    IF(dbg) THEN
       WRITE(*,*)routineN//" nrows,nrow_block,nblkrows=",&
            nrows,nrow_block,nblkrows
       WRITE(*,*)routineN//" ncols,ncol_block,nblkcols=",&
            ncolumns,ncol_block,nblkcols
    ENDIF
    ! Calculate process row distribution
    nprows = dbcsr_mp_nprows (mp_env)
    FORALL (pdim = 0 : nprows-1)
       rd_data(1+pdim:nblkrows:nprows) = pdim
    END FORALL
    CALL array_new (row_dist, rd_data, gift=.TRUE.)
    ! Calculate process column distribution
    npcols = dbcsr_mp_npcols (mp_env)
    FORALL (pdim = 0 : npcols-1)
       cd_data(1+pdim:nblkcols:npcols) = pdim
    END FORALL
    CALL array_new (col_dist, cd_data, gift=.TRUE.)
    !
    IF(dbg) THEN
       WRITE(*,*)routineN//" row_dist",&
            rd_data
       WRITE(*,*)routineN//" col_dist",&
            cd_data
    ENDIF
    !
    CALL dbcsr_distribution_new (dist, mp_env,&
         row_dist, col_dist)
    !
    CALL array_release (row_dist)
    CALL array_release (col_dist)
  END SUBROUTINE dbcsr_create_dist_block_cyclic



! *****************************************************************************
!> \brief Creates a transposed BLACS environment.
!> \par This routine assumes that all BLACS environments are created equal.
!>      It just fiddles with the generated data.
!> \param[in] blacs_env       original BLACS environment.
!> \param[out] blacs_env_t    transposed BLACS environment.
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE cp_transpose_blacs_env(blacs_env, blacs_env_t, error)
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env, blacs_env_t
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_transpose_blacs_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi, mpi2blacs
    LOGICAL                                  :: dbg
    TYPE(cp_para_env_type), POINTER          :: para_env

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    CALL get_blacs_info(blacs_env, para_env=para_env,&
         blacs2mpi=blacs2mpi, mpi2blacs=mpi2blacs)
    !
    !> Create a default
    NULLIFY(blacs_env_t)
    CALL cp_blacs_env_create(blacs_env_t, para_env, blacs_repeatable=blacs_env%repeatable, error=error)
    DEALLOCATE(blacs_env_t%blacs2mpi, blacs_env_t%mpi2blacs)
    blacs_env_t%num_pe(1:2) = blacs_env_t%num_pe(2:1:-1)
    blacs_env_t%mepos(1:2) = blacs_env_t%mepos(2:1:-1)
    !
    ALLOCATE(blacs_env_t%blacs2mpi(LBOUND(blacs2mpi,2):UBOUND(blacs2mpi,2),&
         &                         LBOUND(blacs2mpi,1):UBOUND(blacs2mpi,1)),&
         stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'blacs2mpi_t', SIZE(blacs2mpi))
    blacs_env_t%blacs2mpi = TRANSPOSE(blacs2mpi)
    !
    IF (dbg) WRITE(*,*)'blacs2mpi new',SIZE(blacs_env_t%blacs2mpi,1),&
         SIZE(blacs_env_t%blacs2mpi,2),&
         ';',blacs_env_t%blacs2mpi
    !
    ALLOCATE(blacs_env_t%mpi2blacs(LBOUND(mpi2blacs,2):UBOUND(mpi2blacs,2),&
         &                         LBOUND(mpi2blacs,1):UBOUND(mpi2blacs,1)),&
         stat=stat)
    IF (stat /= 0) CALL stop_memory(routineP, 'mpi2blacs_t', SIZE(mpi2blacs))
    blacs_env_t%mpi2blacs = TRANSPOSE(mpi2blacs)
  END SUBROUTINE cp_transpose_blacs_env


! *****************************************************************************
!> \brief Adds data to the index. Incleases the index size when neccessary.
!> \param[in,out] matrix      bcsr matrix
!> \param[in] slot  which index array to add (e.g., dbcsr_slot_row_blk_sizes)
!> \param[in] data  (optional) array holding the index data to add to the index
!>                  array
!> \param[in] reservation     (optional) only reserve space for subsequent
!>                            array
!> \param[in] extra (optional) reserve extra space for later additions
! *****************************************************************************
  SUBROUTINE dbcsr_addto_index_array(matrix, slot, DATA, reservation, extra)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, INTENT(IN)                      :: slot
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: DATA
    INTEGER, INTENT(IN), OPTIONAL            :: reservation, extra

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_addto_index_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: deplus, space, timing_handle, &
                                                ub, ub_new
    LOGICAL                                  :: dbg

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    CALL timeset(routineN, timing_handle)
    CALL cp_assert (ASSOCIATED (matrix%index), cp_fatal_level,&
         cp_internal_error, routineP,&
         "Index must be preallocated.")
    CALL cp_assert (UBOUND(matrix%index,1).GE.dbcsr_num_slots,&
         cp_failure_level, cp_internal_error, routineP,&
         "Actual index size less than declared size")
    CALL cp_assert(PRESENT(DATA) .OR. PRESENT(reservation), cp_fatal_level,&
         cp_caller_error, routineP,&
         'Either an array or its size must be specified.')
    IF(dbg) WRITE(*,*)routineP//' index',&
         matrix%index(:dbcsr_num_slots)
    IF (PRESENT (reservation)) THEN
       space = reservation
    ELSE
       space = SIZE(DATA)
    ENDIF
    IF (PRESENT (extra)) THEN
       deplus = extra
    ELSE
       deplus = 0
    ENDIF
    ub = UBOUND(matrix%index,1)
    !> The data area was not defined or the new area is greater than the old.
    IF (matrix%index(slot).EQ.0 .OR.&
         space.GT.matrix%index(slot+1)-matrix%index(slot)+1) THEN
       IF(dbg) WRITE(*,*)routineP//' Slot',slot,'not filled, adding at',&
            matrix%index(dbcsr_slot_size)+1,'sized',space
       matrix%index(slot) = matrix%index(dbcsr_slot_size)+1
       matrix%index(slot+1) = matrix%index(slot) + space - 1
       matrix%index(dbcsr_slot_size) = matrix%index(slot+1)
    ENDIF
    ! Shorten an index entry.
    IF (space .LT. matrix%index(slot+1) - matrix%index(slot)+1) THEN
       IF(dbg) WRITE(*,*)routineP//' Shortening index to ubound', ub_new
       matrix%index(slot+1) = matrix%index(slot) + space -1
       CALL dbcsr_repoint_index(matrix, slot)
    ENDIF
    ub_new = matrix%index(slot+1) + deplus
    IF(dbg) WRITE(*,*)routineP//' need',space,'at',matrix%index(slot),&
         'to',matrix%index(slot+1),'(',ub_new,')','have',ub
    IF (ub_new .GT. ub) THEN
       IF(dbg) WRITE(*,*)routineP//' Reallocating index to ubound', ub_new
       !CALL reallocate(matrix%index, 1, ub_new)
       CALL ensure_array_size(matrix%index, 1, ub_new, factor=1.0,&
            nocopy=.FALSE., special=matrix%special_memory)
       CALL dbcsr_repoint_index(matrix)
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Adding slot',slot,'at',&
         matrix%index(slot),'sized',space
    CALL dbcsr_repoint_index(matrix, slot)
    IF (PRESENT(DATA)) &
         matrix%index(matrix%index(slot):matrix%index(slot+1)) = DATA(:)
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_addto_index_array


! *****************************************************************************
!> \brief Updates the index pointers of a bcsr matrix
!> \param[in,out] m           matrix for which index pointers are updated
!> \param[in] slot            (optional) only repoint this index
! *****************************************************************************
  SUBROUTINE dbcsr_repoint_index(m, slot)
    TYPE(dbcsr_type), INTENT(INOUT)          :: m
    INTEGER, INTENT(IN), OPTIONAL            :: slot

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_repoint_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: s
    LOGICAL                                  :: all

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    DBG 'Repointing index for matrix ',m%name,'with #blocks',&
         m%nblks,m%index(dbcsr_slot_nblks)
    IF (m%nblks .NE. m%index(dbcsr_slot_nblks)) THEN
       m%nblks = m%index(dbcsr_slot_nblks)
       m%nze = m%index(dbcsr_slot_nze)
    ENDIF
    all = .TRUE.
    IF (PRESENT (slot)) THEN
       all = .FALSE.
       s = slot
       DBG 'only slot',slot, m%index(slot)
    ELSE
       s = 0
    ENDIF
    DBGV 'index',m%index(1:dbcsr_num_slots)
    DBG 'all?',all,'s=',s

    IF (m%index(dbcsr_slot_row_p).GT.0.AND.all .OR.&
         s.EQ.dbcsr_slot_row_p) THEN
       DBG ' setting pointer row_p'
       m%row_p => m%index(m%index(dbcsr_slot_row_p):&
            &                  m%index(dbcsr_slot_row_p+1))
       DBG ' set pointer row_p, size',&
            m%index(dbcsr_slot_row_p+1)-m%index(dbcsr_slot_row_p)+1
    ENDIF
    IF (m%index(dbcsr_slot_col_i).GT.0.AND.all .OR.&
         s.EQ.dbcsr_slot_col_i) THEN
       DBG ' setting pointer col_i'
       m%col_i => m%index(m%index(dbcsr_slot_col_i):&
            &                  m%index(dbcsr_slot_col_i+1))
       DBG ' set pointer col_i, size',&
            m%index(dbcsr_slot_col_i+1)-m%index(dbcsr_slot_col_i)+1
    ENDIF
    IF (m%index(dbcsr_slot_blk_p).GT.0.AND.all .OR.&
         s.EQ.dbcsr_slot_blk_p) THEN
       DBG ' setting pointer blk_p'
       m%blk_p => m%index(m%index(dbcsr_slot_blk_p):&
            &                  m%index(dbcsr_slot_blk_p+1))
       DBG ' set pointer blk_p, size',&
            m%index(dbcsr_slot_blk_p+1)-m%index(dbcsr_slot_blk_p)+1
    ENDIF
    IF (all) THEN
       m%index(dbcsr_slot_nblks) = m%nblks
       m%index(dbcsr_slot_nze) = m%nze
    ENDIF
  END SUBROUTINE dbcsr_repoint_index


! *****************************************************************************
!> \brief Calculates the checksum of a DBCSR matrix.
!> \param[in] matrix          matrix
!> \param[out] chksum         calculated checksum
!> \param[in,out] error       cp2k error
!> \param[in] local           (optional) no global communication
! *****************************************************************************
  FUNCTION dbcsr_checksum(matrix, error, local) RESULT(checksum)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(cp_error_type), INTENT(IN)          :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: local
    REAL(KIND=dp)                            :: checksum

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_checksum', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: c_dp
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: c_sp
    INTEGER                                  :: bc, blk, blk_p, br, m, mn, n, &
                                                timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: nocomm
    REAL(KIND=dp)                            :: blk_cs, local_cs
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_dp
    REAL(KIND=sp), DIMENSION(:), POINTER     :: r_sp
    TYPE(cp_error_type)                      :: e

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    CALL cp_assert (dbcsr_valid_index(matrix),&
         cp_warning_level, cp_caller_error,&
         routineN, "Invalid matrix.")
    nocomm = .FALSE.
    IF (PRESENT (local)) nocomm = local
    e = error
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    local_cs = 0.0_dp
    SELECT CASE (matrix%m%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (matrix%m%data_area, r_dp)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (matrix%m%data_area, r_sp)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (matrix%m%data_area, c_dp)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (matrix%m%data_area, c_sp)
    END SELECT
    DO br = 1, matrix%m%nblkrows_total
       m = row_blk_size(br)
       DO blk = matrix%m%row_p(br)+1, matrix%m%row_p(br+1)
          bc = matrix%m%col_i(blk)
          n = col_blk_size(bc)
          mn = m*n
          blk_p = ABS(matrix%m%blk_p(blk))
          IF (blk_p .NE. 0) THEN
             ! Calculate DDOT
             SELECT CASE (matrix%m%data_type)
             CASE (dbcsr_type_real_8)
                blk_cs = DOT_PRODUCT(r_dp(blk_p:blk_p+mn-1),&
                     &               r_dp(blk_p:blk_p+mn-1))
             CASE (dbcsr_type_real_4)
                blk_cs = DOT_PRODUCT(r_sp(blk_p:blk_p+mn-1),&
                     &               r_sp(blk_p:blk_p+mn-1))
             CASE (dbcsr_type_complex_8)
                blk_cs = DOT_PRODUCT(c_dp(blk_p:blk_p+mn-1),&
                     &               c_dp(blk_p:blk_p+mn-1))
             CASE (dbcsr_type_complex_4)
                blk_cs = DOT_PRODUCT(c_sp(blk_p:blk_p+mn-1),&
                     &               c_sp(blk_p:blk_p+mn-1))
             END SELECT
             DBGV ' chksum for (',br,',',bc,') at', blk_p,'l',mn,'= ', blk_cs
             local_cs = local_cs + blk_cs
          ENDIF
       ENDDO
    ENDDO
    DBG 'local chksum=',local_cs
    checksum = local_cs
    IF (.NOT. nocomm) THEN
       CALL mp_sum(local_cs, dbcsr_mp_group (dbcsr_distribution_mp (&
            matrix%m%dist)))
       checksum = local_cs
       DBGI 'local chksum=',local_cs,'global=',checksum
    ENDIF
    CALL timestop(timing_handle)
  END FUNCTION dbcsr_checksum

! *****************************************************************************
!> \brief Calculates the sum of all DBCSR matrix elements.
!> \param[in] matrix          matrix
!> \param[out] chksum         calculated sum
!> \param[in,out] error       cp2k error
!> \param[in] local           (optional) no global communication
! *****************************************************************************
  FUNCTION dbcsr_ssum(matrix, local) RESULT(ssum)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: local
    REAL(KIND=dp)                            :: ssum

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ssum', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: c_dp
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: c_sp
    INTEGER                                  :: bc, blk, blk_p, br, m, mn, n, &
                                                timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: nocomm
    REAL(KIND=dp)                            :: blk_cs, local_cs
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_dp
    REAL(KIND=sp), DIMENSION(:), POINTER     :: r_sp

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    CALL cp_assert (dbcsr_valid_index(matrix),&
         cp_warning_level, cp_caller_error,&
         routineN, "Invalid matrix.")
    nocomm = .FALSE.
    IF (PRESENT (local)) nocomm = local
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    local_cs = 0.0_dp
    SELECT CASE (matrix%m%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (matrix%m%data_area, r_dp)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (matrix%m%data_area, r_sp)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (matrix%m%data_area, c_dp)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (matrix%m%data_area, c_sp)
    END SELECT
    DO br = 1, matrix%m%nblkrows_total
       m = row_blk_size(br)
       DO blk = matrix%m%row_p(br)+1, matrix%m%row_p(br+1)
          bc = matrix%m%col_i(blk)
          n = col_blk_size(bc)
          mn = m*n
          blk_p = ABS(matrix%m%blk_p(blk))
          IF (blk_p .NE. 0) THEN
             ! Calculate DDOT
             SELECT CASE (matrix%m%data_type)
             CASE (dbcsr_type_real_8)
                blk_cs = SUM(r_dp(blk_p:blk_p+mn-1))
             CASE (dbcsr_type_real_4)
                blk_cs = SUM(r_sp(blk_p:blk_p+mn-1))
             CASE (dbcsr_type_complex_8)
                blk_cs = SUM(c_dp(blk_p:blk_p+mn-1))
             CASE (dbcsr_type_complex_4)
                blk_cs = SUM(c_sp(blk_p:blk_p+mn-1))
             END SELECT
             DBGV ' chksum for (',br,',',bc,') at', blk_p,'l',mn,'= ', blk_cs
             local_cs = local_cs + blk_cs
          ENDIF
       ENDDO
    ENDDO
    DBG 'local chksum=',local_cs
    ssum = local_cs
    IF (.NOT. nocomm) THEN
       CALL mp_sum(local_cs, dbcsr_mp_group (dbcsr_distribution_mp (&
            matrix%m%dist)))
       ssum = local_cs
       DBGI 'local chksum=',local_cs,'global=',ssum
    ENDIF
    CALL timestop(timing_handle)
  END FUNCTION dbcsr_ssum


! *****************************************************************************
!> \brief Prints a BCSR matrix (block-style, not full)
!> \param[in] matrix          matrix
!> \param[in,out] error       cp2k error
!> \param[in] nodata          (optional) don't print actual data
! *****************************************************************************
  SUBROUTINE dbcsr_print(matrix, error, nodata)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: nodata

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_print', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: c_dp
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: c_sp
    INTEGER                                  :: bc, blk, blk_p, br, m, mn, n, &
                                                timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: tr, yesprint
    REAL(KIND=dp)                            :: blk_cs
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_dp
    REAL(KIND=sp), DIMENSION(:), POINTER     :: r_sp

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    CALL cp_assert (dbcsr_valid_index (matrix),&
         cp_warning_level, cp_caller_error,&
         routineN, "Can not print invalid matrix.")
    yesprint = .TRUE.
    IF (PRESENT (nodata)) yesprint = .NOT. nodata
    WRITE(*,*)routineP//' Contents of matrix named ', matrix%m%name
    WRITE(*,*)routineP//' Flags ', matrix%m%matrix_type
    WRITE(*,'(1X,A,I5,A,I5)')routineP//" Full size",&
         matrix%m%nfullrows_total, "x",matrix%m%nfullcols_total
    WRITE(*,'(1X,A,I5,A,I5)')routineP//" Blocked size",&
         matrix%m%nblkrows_total, "x",matrix%m%nblkcols_total
    SELECT CASE (matrix%m%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (matrix%m%data_area, r_dp)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (matrix%m%data_area, r_sp)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (matrix%m%data_area, r_dp)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (matrix%m%data_area, r_sp)
    END SELECT
    DBGV 'index', matrix%m%index
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    IF (matrix%m%nblks .GT. 0) THEN
       DO br = 1, matrix%m%nblkrows_total
          m = row_blk_size(br)
          DO blk = matrix%m%row_p(br)+1, matrix%m%row_p(br+1)
             bc = matrix%m%col_i(blk)
             n = col_blk_size(bc)
             mn = m*n
             blk_p = ABS(matrix%m%blk_p(blk))
             tr = matrix%m%blk_p(blk) .LT. 0
             block_exists: IF (blk_p .NE. 0) THEN
                SELECT CASE (matrix%m%data_type)
                CASE (dbcsr_type_real_8)
                   blk_cs = DOT_PRODUCT(r_dp(blk_p:blk_p+mn-1),&
                        &               r_dp(blk_p:blk_p+mn-1))
                   !CALL &
                   !     dbcsr_printmat(r_dp(blk_p:blk_p+mn-1),m,n, tr=tr)
                CASE (dbcsr_type_real_4)
                   blk_cs = DOT_PRODUCT(r_sp(blk_p:blk_p+mn-1),&
                        &               r_sp(blk_p:blk_p+mn-1))
                   !CALL &
                   !     dbcsr_printmat(r_sp(blk_p:blk_p+mn-1),m,n, tr=tr)
                CASE (dbcsr_type_complex_8)
                   blk_cs = DOT_PRODUCT(c_dp(blk_p:blk_p+mn-1),&
                        &               c_dp(blk_p:blk_p+mn-1))
                   !CALL &
                   !     dbcsr_printmat(c_dp(blk_p:blk_p+mn-1),m,n, tr=tr)
                CASE (dbcsr_type_complex_4)
                   blk_cs = DOT_PRODUCT(c_sp(blk_p:blk_p+mn-1),&
                        &               c_sp(blk_p:blk_p+mn-1))
                   !CALL &
                   !     dbcsr_printmat(c_sp(blk_p:blk_p+mn-1),m,n, tr=tr)
                END SELECT
                !WRITE(*,*)routineP//' chksum for (',br,',',bc,') at',&
                !     blk_p,'l',mn,'= ', blk_cs,'size',m,n
                WRITE(*,'(A,I6,",",I6,A,I7,A,I6,I6,"=",I7,A,E12.3)')&
                     !" Checksum for (",br,bc,") at ",blk_p," size ",m,n,mn,&
                     " Checksum for (",br,bc,") at ",matrix%m%blk_p(blk)," size ",m,n,mn,&
                     " checksum=",blk_cs
                IF (yesprint .AND. blk_p .NE. 0) THEN
                   SELECT CASE (matrix%m%data_type)
                   CASE (dbcsr_type_real_8)
                      !WRITE(*,'(10(1X,F7.2))')r_dp(blk_p:blk_p+mn-1)
                      CALL dbcsr_printmat(r_dp(blk_p:blk_p+mn-1),m,n, tr=tr)
                   CASE (dbcsr_type_real_4)
                      CALL dbcsr_printmat(r_sp(blk_p:blk_p+mn-1),m,n, tr=tr)
                   CASE (dbcsr_type_complex_8)
                      CALL dbcsr_printmat(c_dp(blk_p:blk_p+mn-1),m,n, tr=tr)
                   CASE (dbcsr_type_complex_4)
                      CALL dbcsr_printmat(c_sp(blk_p:blk_p+mn-1),m,n, tr=tr)
                   END SELECT
                ENDIF
             ENDIF block_exists
          ENDDO
       ENDDO
    ENDIF
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_print



! *****************************************************************************
!> \brief Verify the correctness of a BCSR matrix.
!> \param[in] m     bcsr matrix
!> \param[in,out] error     cp2k error
!> \param[in] verbosity     how detailed errors are; 0=nothing; 1=summary at
!>                          end if matrix not consistent; 2=also individual
!>                          errors; 3=always print info about matrix; >3=even
!>                          more info
!> \param[in] local         (optional) no global communication
! *****************************************************************************
  SUBROUTINE dbcsr_verify_matrix(m, error, verbosity, local)
    TYPE(dbcsr_obj), INTENT(IN)              :: m
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    INTEGER, INTENT(IN), OPTIONAL            :: verbosity
    LOGICAL, INTENT(IN), OPTIONAL            :: local

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_verify_matrix', &
      r = moduleN//':'//routineN

    INTEGER :: bc, blk, blk_p, br, dbg, i, mb, mn, n, n_full_blocks_total, &
      n_have_blocks_local, n_have_blocks_total, prev_br, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: nocomm
    REAL(KIND=dp)                            :: sparsity_total

!   ---------------------------------------------------------------------------

    CALL cp_assert (m%m%initialized .EQ. dbcsr_magic_number,&
         cp_fatal_level, cp_caller_error,&
         routineN, "Can not verify non-existing matrix object.")
    dbg = 2
    nocomm = .FALSE.
    IF (PRESENT (local)) nocomm = local
    IF (PRESENT (verbosity)) dbg = verbosity
    CALL timeset(routineN, timing_handle)
    IF (dbg.GE.3) WRITE(*,'(1X,A,A,A,A)')r//'Matrix name: ', m%m%name,&
         " of type ",m%m%matrix_type
    IF (dbg.GE.3) THEN
       WRITE(*,'(1X,A,I5,"x",I5,A,I5,"x",I5)')r//' Size blocked',&
            m%m%nblkrows_total, m%m%nblkcols_total, ", full ",&
            m%m%nfullrows_total, m%m%nfullcols_total
    ENDIF
    row_blk_size => array_data (m%m%row_blk_size)
    col_blk_size => array_data (m%m%col_blk_size)
    !
    IF(m%m%matrix_type.EQ.'N' .OR. m%m%matrix_type.EQ.'T') THEN
       n_full_blocks_total = m%m%nblkrows_total * m%m%nblkcols_total
    ELSE
       CALL bcsr_assert(dbg, m%m%nblkrows_total.EQ.m%m%nblkcols_total, r,&
            'Symmetric matrix is not square')
       n_full_blocks_total = m%m%nblkrows_total * (m%m%nblkrows_total+1)/2
    ENDIF
    n_have_blocks_local = m%m%nblks
2045 FORMAT (I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5)
2047 FORMAT (I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7,1X,I7)
    IF (dbg.GE.4) THEN
       WRITE(*,'(1X,A)')r//' index='
       WRITE(*,2045)m%m%index(:dbcsr_num_slots)
    ENDIF
    CALL bcsr_assert(dbg, m%m%index(1).GT.0, r, 'Index size 0')
    DO i = dbcsr_slot_row_p, dbcsr_num_slots
       CALL bcsr_assert(dbg, m%m%index(i) .GT. 0, r, 'Index member is 0')
       IF (.NOT. (i.EQ.dbcsr_slot_col_i.OR.i.EQ.dbcsr_slot_blk_p)) THEN
            CALL bcsr_assert(dbg,m%m%index(i) .LE. m%m%index(1), r,&
            'Index member is greater than size')
         ENDIF
    ENDDO
    !
    IF(dbg.GE.4) WRITE(*,*)r//' row_p extents',m%m%index(dbcsr_slot_row_p+1),&
         m%m%index(dbcsr_slot_row_p), SIZE(m%m%row_p)
    CALL bcsr_assert(dbg,m%m%index(dbcsr_slot_row_p+1)-m%m%index(dbcsr_slot_row_p)+1&
         .EQ.m%m%nblkrows_total+1, r,&
         'Size of row_p index inconsistent with number of rows')
    CALL bcsr_assert(dbg,SIZE(m%m%row_p).EQ.m%m%nblkrows_total+1, r,&
         'Size of row_p inconsistent with number of rows')
    !
    IF(dbg.GE.4) WRITE(*,*)r//' col_i extents',m%m%index(dbcsr_slot_col_i+1),&
         m%m%index(dbcsr_slot_col_i), SIZE(m%m%col_i)
    CALL bcsr_assert(dbg,m%m%index(dbcsr_slot_col_i+1)-m%m%index(dbcsr_slot_col_i)+1&
         .EQ.m%m%nblks, r,&
         'Size of col_i index inconsistent with number of blocks')
    CALL bcsr_assert(dbg,SIZE(m%m%col_i).EQ.m%m%nblks, r,&
         'Size of col inconsistent with number of blocks')
    !
    IF(dbg.GE.4) WRITE(*,*)r//' blk_p extents',m%m%index(dbcsr_slot_blk_p+1),&
         m%m%index(dbcsr_slot_blk_p), SIZE(m%m%blk_p)
    CALL bcsr_assert(dbg,m%m%index(dbcsr_slot_blk_p+1)-m%m%index(dbcsr_slot_blk_p)+1&
         .EQ.m%m%nblks, r,&
         'Size of blk_p index inconsistent with number of blocks')
    CALL bcsr_assert(dbg,SIZE(m%m%col_i).EQ.m%m%nblks, r,&
         'Size of blk_p inconsistent with number of blocks')
    !
    CALL bcsr_assert(dbg,SIZE(row_blk_size).EQ.m%m%nblkrows_total, r,&
         'Row block size array inconsistent with number of blocked rows')
    CALL bcsr_assert(dbg,SIZE(col_blk_size).EQ.m%m%nblkcols_total, r,&
         'Column block size array inconsistent with number of blocked columns')
    !
    IF(dbg.GE.4) THEN
       WRITE(*,'(1X,A,I7,A,I7)')r//' nze=',m%m%nze,'data size',&
            dbcsr_get_data_size(m%m%data_area)
    ENDIF
    CALL bcsr_assert(dbg,dbcsr_get_data_size(m%m%data_area).GE.m%m%nze, r,&
         'Data storage may be too small.')
    IF (dbg.GE.5) THEN
       WRITE (*,'(1X,A,I7,A)') r//' size=',SIZE(m%m%row_p),' row_p='
       WRITE (*,2047) m%m%row_p(1:m%m%nblkrows_total+1)
       WRITE (*,'(1X,A)') r//' col_i='
       WRITE (*,2047) m%m%col_i(1:m%m%nblks)
       WRITE (*,'(1X,A)') r//' blk_p='
       WRITE (*,2047) m%m%blk_p(1:m%m%nblks)
    ENDIF
    prev_br = 0
    DO br = 1, m%m%nblkrows_total
       CALL bcsr_assert(dbg,m%m%row_p(br).GE.0, r, 'row_p less than zero')
       IF(br.GT.1) CALL bcsr_assert(dbg,m%m%row_p(br).GE.m%m%row_p(prev_br),&
            r, 'row_p decreases')
       mb = row_blk_size(br)
       CALL bcsr_assert(dbg,mb.GT.0, r, 'Row blocked size is zero')
       DO blk = m%m%row_p(br)+1, m%m%row_p(br+1)
          CALL bcsr_assert(dbg,blk.GT.0, r, 'Block number is zero')
          CALL bcsr_assert(dbg,blk.LE.m%m%nblks, r, 'Block number too high')
          bc = m%m%col_i(blk)
          IF (dbg.GE.5) THEN
             WRITE(*,'(1X,A,I7,"(",I5,",",I5,")")')r//' block',blk,br,bc
          ENDIF
          CALL bcsr_assert(dbg,bc.GT.0, r, 'col_i is zero')
          CALL bcsr_assert(dbg,bc.LE.m%m%nblkcols_total, r, 'col_i too high')
          n = col_blk_size(bc)
          CALL bcsr_assert(dbg,n.GT.0, r, 'Column blocked size is zero')
          blk_p = m%m%blk_p(blk)
          mn = mb*n
          !CALL bcsr_assert(dbg,blk_p.GT.0, r, 'Block pointer is negative')
          CALL bcsr_assert(dbg,blk_p.LE.m%m%nze, r, 'Block pointer too large')
          CALL bcsr_assert(dbg,blk_p+mn-1.LE.m%m%nze, r, 'Block extends too far')
       ENDDO
       prev_br = br
    ENDDO
    IF (dbg.GE.3.AND..NOT.nocomm) THEN
       CALL mp_sum(n_have_blocks_local, dbcsr_mp_group (dbcsr_distribution_mp (&
            m%m%dist)))
       n_have_blocks_total = n_have_blocks_local
       sparsity_total = REAL(n_have_blocks_total)&
            / REAL(n_full_blocks_total)*100.0
       !WRITE(*,FMT='(30A,F5.1,A)')r//' Sparsity: ', sparsity_total,'%'
       WRITE(*,FMT='(1X,A,F5.1,A)')r//' Non-sparsity: ',&
            sparsity_total,'%'
    ENDIF
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_verify_matrix


! *****************************************************************************
!> \brief Glue for the real assert code
!> \param[in] l     Assertion
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
! *****************************************************************************
  SUBROUTINE bcsr_assert(level, l, routine, msg)
    INTEGER, INTENT(IN)                      :: level
    LOGICAL, INTENT(IN)                      :: l
    CHARACTER(len=*), INTENT(IN)             :: routine, msg

!   ---------------------------------------------------------------------------

    IF (level .GE. 2) THEN
       CALL cp_assert(l, cp_warning_level, cp_internal_error, routine, msg)
    ENDIF
  END SUBROUTINE bcsr_assert

! *****************************************************************************
!> \brief Assertion
!> \param[in] left            left value
!> \param[in] rel             relation
!> \param[in] right           right value
!> \param[in] level           error level
!> \param[in] etype           error type
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
! *****************************************************************************
  SUBROUTINE dbcsr_int_assert(left, rel, right, level, etype, routine, msg)
    INTEGER, INTENT(IN)                      :: left
    CHARACTER(len=2), INTENT(IN)             :: rel
    INTEGER, INTENT(IN)                      :: right, level, etype
    CHARACTER(len=*), INTENT(IN)             :: routine, msg

    LOGICAL                                  :: l

!   ---------------------------------------------------------------------------

    SELECT CASE (rel)
    CASE ("EQ")
       l = left .EQ. right
    CASE ("LT")
       l = left .LT. right
    CASE ("LE")
       l = left .LE. right
    CASE ("GT")
       l = left .GT. right
    CASE ("GE")
       l = left .GE. right
    CASE ("NE")
       l = left .NE. right
    CASE default
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error,&
            "dbcsr_int_assert", "Invalid relation specified: "//rel)
    END SELECT
    IF (.NOT. l) THEN
       WRITE(*,'(1X,A,1X,I9,A,I9)')"ASSERTION FAILED:",&
            left, "."//rel//".", right
    ENDIF
    CALL cp_assert(l, cp_fatal_level, cp_internal_error, routine, msg)
  END SUBROUTINE dbcsr_int_assert

! *****************************************************************************
!> \brief Assertion
!> \param[in] left            left value
!> \param[in] rel             relation
!> \param[in] right           right value
!> \param[in] level           error level
!> \param[in] etype           error type
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
! *****************************************************************************
  SUBROUTINE dbcsr_logical_assert(left, rel, right, level, etype, routine, msg)
    LOGICAL, INTENT(IN)                      :: left
    CHARACTER(len=*), INTENT(IN)             :: rel
    LOGICAL, INTENT(IN)                      :: right
    INTEGER, INTENT(IN)                      :: level, etype
    CHARACTER(len=*), INTENT(IN)             :: routine, msg

    LOGICAL                                  :: l

!   ---------------------------------------------------------------------------

    SELECT CASE (rel)
    CASE ("EQV")
       l = left .EQV. right
    CASE ("NEQV")
       l = left .NEQV. right
    CASE ("OR")
       l = left .OR. right
    CASE ("AND")
       l = left .AND. right
    CASE default
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error,&
            "dbcsr_int_assert", "Invalid relation specified: "//rel)
    END SELECT
    IF (.NOT. l) THEN
       WRITE(*,'(1X,A,1X,L1,A,L1)')"ASSERTION FAILED:",&
            left, "."//rel//".", right
    ENDIF
    CALL cp_assert(l, cp_fatal_level, cp_internal_error, routine, msg)
    END SUBROUTINE dbcsr_logical_assert

! *****************************************************************************
!> \brief Assertion
!> \param[in] rel             relation
!> \param[in] right           right value
!> \param[in] level           error level
!> \param[in] etype           error type
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
! *****************************************************************************
  SUBROUTINE dbcsr_not_assert(rel, right, level, etype, routine, msg)
    CHARACTER(len=3), INTENT(IN)             :: rel
    LOGICAL, INTENT(IN)                      :: right
    INTEGER, INTENT(IN)                      :: level, etype
    CHARACTER(len=*), INTENT(IN)             :: routine, msg

    LOGICAL                                  :: l

!   ---------------------------------------------------------------------------

    SELECT CASE (rel)
    CASE ("NOT")
       l = .NOT. right
    CASE default
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error,&
            "dbcsr_int_assert", "Invalid relation specified: "//rel)
    END SELECT
    IF (.NOT. l) THEN
       WRITE(*,'(1X,A,1X,A,L1)')"ASSERTION FAILED:",&
            "."//rel//".", right
    ENDIF
    CALL cp_assert(l, cp_fatal_level, cp_internal_error, routine, msg)
  END SUBROUTINE dbcsr_not_assert


!! *****************************************************************************
!!> \brief Makes a matrix canonical: sorts the indices and makes the data
!!>        sequential.
!!> \param[in,out] matrix      matrix in which to sort the data blocks
!!> \param[in,out] error       cp2k error
!!> \param[in] onlydata        (optional) only sort data
!! *****************************************************************************
!  SUBROUTINE dbcsr_order(matrix, error, onlydata)
!    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
!    TYPE(cp_error_type), INTENT(INOUT)       :: error
!    LOGICAL, INTENT(IN), OPTIONAL            :: onlydata
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_order', &
!      routineP = moduleN//':'//routineN
!
!    INTEGER                                  :: blk, col_size, row, row_size, &
!                                                stat, timing_handle
!    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sizes
!    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
!    REAL(KIND=dp), DIMENSION(:), POINTER     :: new_data
!
!    DEBUG_HEADER
!!   ---------------------------------------------------------------------------
!    DEBUG_BODY
!    CALL timeset(routineN, timing_handle)
!    ! We should also check that the indices actually exist.
!    CALL cp_assert(ASSOCIATED(matrix%data), cp_failure_level, cp_caller_error,&
!         routineN, "Matrix has no data.", error)
!    row_blk_size => array_data (matrix%row_blk_size)
!    col_blk_size => array_data (matrix%col_blk_size)
!    !
!    IF (PRESENT (onlydata)) THEN
!       IF (onlydata) THEN
!          CALL dbcsr_sort_indices(matrix%nblks, matrix%row_p, matrix%col_i,&
!               matrix%blk_p, error=error)
!       ENDIF
!    ENDIF
!    ALLOCATE(sizes(matrix%nblks), stat=stat)
!    IF (stat /= 0) CALL stop_memory(routineN, "sizes", matrix%nblks)
!    !$OMP PARALLEL DO PRIVATE(row, blk, row_size, col_size)
!    DO row = 1, matrix%nblkrows_total
!       row_size = row_blk_size(row)
!       DO blk = matrix%row_p(row)+1, matrix%row_p(row+1)
!          col_size = col_blk_size(matrix%col_i(blk))
!          sizes(blk) = row_size*col_size
!       ENDDO
!    ENDDO
!    !$OMP END PARALLEL DO
!    IF (matrix%special_memory) THEN
!       CALL mp_allocate(new_data, SIZE(matrix%data), stat=stat)
!    ELSE
!       ALLOCATE(new_data(SIZE(matrix%data)), stat=stat)
!    ENDIF
!    IF (stat /= 0) CALL stop_memory(routineN, "new_data", SIZE(matrix%data))
!    CALL dbcsr_sort_data(matrix%blk_p, matrix%blk_p, sizes,&
!         new_data, matrix%data, error=error)
!    IF (matrix%special_memory) THEN
!       CALL mp_deallocate(matrix%data)
!    ELSE
!       DEALLOCATE(matrix%data)
!    ENDIF
!    matrix%data => new_data
!    DEALLOCATE(sizes)
!    CALL timestop(timing_handle)
!  END SUBROUTINE dbcsr_order


! *****************************************************************************
!> \brief Sorts the rows & columns of a work matrix
!> \par Description
!>      Sorts the row and column indices so that the rows monotonically
!>      increase and the columns monotonically increase within each row.
!>      Passing the blk_p array rearranges the block pointers accordingly.
!>      This must be done if they are pointing to valid data, otherwise
!>      they become invalid.
!> \param[in] n               number of blocks (elements) to sort
!> \param[in,out] row_i       row indices
!> \param[in,out] col_i       column indices
!> \param[in,out] blk_p       (optional) block pointers
!> \param[in,out] blk_d       (optional) data storage
! *****************************************************************************
  SUBROUTINE dbcsr_sort_indices(n, row_i, col_i, blk_p, blk_d)
    INTEGER, INTENT(IN)                      :: n
    INTEGER, DIMENSION(1:n), INTENT(INOUT)   :: row_i, col_i
    INTEGER, DIMENSION(1:n), INTENT(INOUT), &
      OPTIONAL                               :: blk_p, blk_d

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sort_indices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: buf, buf_d

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    IF (SIZE (row_i) .EQ. 0) RETURN
    CALL timeset(routineN, timing_handle)
    CALL cp_assert(SIZE(row_i).GE.n, cp_failure_level,&
         cp_caller_error, routineP, 'row_i too small')
    CALL cp_assert(SIZE(col_i).GE.n, cp_failure_level,&
         cp_caller_error, routineP, 'col_i too small')
    IF (PRESENT (blk_p)) CALL cp_assert(SIZE(blk_p).GE.n, cp_failure_level,&
         cp_caller_error, routineP, 'blk_p too small')
    IF (PRESENT (blk_p)) THEN
       ALLOCATE(buf(n), stat=stat)
       IF (stat /= 0) CALL stop_memory(routineP, 'buf', n)
       buf(1:n) = blk_p(1:n)
    ENDIF
    IF (PRESENT (blk_d)) THEN
       ALLOCATE(buf_d(n), stat=stat)
       IF (stat /= 0) CALL stop_memory(routineP, 'buf_d', n)
       buf_d(1:n) = blk_d(1:n)
    ENDIF
    !> Create an ordering for both rows and columns while. If the blk_p must be
    !> rearranged, then the col_i array will be used as a permutation vector.
    row_i(:) = ISHFT(row_i(:), 16) + col_i
    IF (PRESENT (blk_p)) col_i(1:n) = (/ (i, i=1,n) /)
    !> Now do a nice quicksort.
    DBG ' sort, row_i=',row_i
    !col_i = (/ (i, i=1, n) /)
    CALL sort(row_i, n, col_i)
    DO i = 1, n
       IF (PRESENT (blk_p)) THEN
          blk_p(i) = buf(col_i(i))
       ENDIF
       IF (PRESENT (blk_d)) THEN
          blk_d(i) = buf_d(col_i(i))
       ENDIF
       col_i(i) = IAND(row_i(i), (ISHFT(1,16)-1))
       row_i(i) = ISHFT(row_i(i), -16)
    ENDDO
    DBG ' sort, row_i=',row_i
    DBG ' sort, col_i=',col_i
    IF(dbg.AND.bcsr_verbose.AND.PRESENT(blk_p))&
         WRITE(*,*)routineP//' sort, blk_p =',blk_p
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_sort_indices

! *****************************************************************************
!> \brief Sorts all work matrices.
!> \param[in,out] matrix      sort work matrices in this matrix
! *****************************************************************************
  SUBROUTINE dbcsr_sort_many_indices(matrix)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sort_many_indices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, timing_handle

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    DO i = 1, SIZE (matrix%wms)
       IF (.NOT. dbcsr_wm_use_mutable (matrix%wms(i))) THEN
          IF (matrix%wms(i)%lastblk .GT. 0) THEN
             CALL dbcsr_sort_indices (matrix%wms(i)%lastblk,&
                  matrix%wms(i)%row_i,&
                  matrix%wms(i)%col_i, matrix%wms(i)%blk_p)
          ENDIF
       ENDIF
    ENDDO
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_sort_many_indices


! *****************************************************************************
!> \brief Calculates explicit sizes for all data blocks.
!> \param[out] sizes          sizes of all data blocks
!> \param[in] row_p, col_i    index structure
!> \param[in] rbs, cbs        row block sizes and column block sizes
! *****************************************************************************
  SUBROUTINE calc_blk_sizes(sizes, row_p, col_i, rbs, cbs)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: sizes
    INTEGER, DIMENSION(:), INTENT(IN)        :: row_p, col_i, rbs, cbs

    INTEGER                                  :: blk, nrows, row, row_size

    nrows = SIZE (row_p)-1
    !$OMP PARALLEL DO PRIVATE (row, blk, row_size)
    DO row = 1, nrows
       row_size = rbs(row)
       FORALL (blk = row_p(row)+1 : row_p(row+1))
          sizes(blk) = row_size * cbs(col_i(blk))
       END FORALL
    ENDDO
    !$OMP END PARALLEL DO
  END SUBROUTINE calc_blk_sizes


! *****************************************************************************
!> \brief Sorts the data in a matrix so that the data blocks follow
!>        sequentially.
!> \param[out] blk_p          re-arranged block pointers reflecting the new
!>                            data order
!> \param[in] old_blk_p       current block pointers
!> \param[in] sizes           sizes of the data blocks
!> \param[out] dst            sorted data
!> \param[in] src             existing unordered data
!> \param[in] srcs            (optional) multiple source areas
!> \param[in] blk_d           (optional) mapping from blocks to source areas
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_sort_data(blk_p, old_blk_p, sizes, dsts, src,&
       srcs, old_blk_d, error)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: blk_p
    INTEGER, DIMENSION(:), INTENT(IN)        :: old_blk_p, sizes
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dsts
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: src
    TYPE(dbcsr_data_obj), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: srcs
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: old_blk_d
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sort_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nblks, timing_handle
    LOGICAL                                  :: multidata

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    multidata = PRESENT (srcs) .AND. PRESENT (old_blk_d)
    nblks = SIZE(old_blk_p)
    DBG 'blk sizes', SIZE(blk_p), nblks
    DBG 'dst data sizes', dbcsr_get_data_size(dsts)
    CALL cp_assert(SIZE(blk_p).GE.nblks, cp_failure_level,&
         cp_wrong_args_error, routineN, 'Destination blk_p too small.', error)
    ! If the source is deliberately bigger, we don't want any warnings.
    !CALL cp_assert(SIZE(dst).GE.SIZE(src), cp_warning_level,&
    !     cp_wrong_args_error, routineN,&
    !     'Destination data space may be too small.',&
    !     error)
    DBGV 'old_blk_p',old_blk_p
    DBGV 'sizes', sizes
    IF (SIZE (old_blk_p) .GT. 0) THEN
       blk_p(1) = SIGN(1, old_blk_p(1))
       DO i = 2, nblks
          blk_p(i) = SIGN(ABS(blk_p(i-1)) + sizes(i-1), old_blk_p(i))
       ENDDO
       !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i)
       DO i = 1, nblks
          IF (.NOT. multidata) THEN
             CALL dbcsr_data_set (dsts,&
                  ABS(blk_p(i)), sizes(i),&
                  src, source_lb=ABS(old_blk_p(i)))
             !dst(ABS(blk_p(i)):ABS(blk_p(i))+sizes(i)-1) =&
             !     src(ABS(old_blk_p(i)):ABS(old_blk_p(i))+sizes(i)-1)
          ELSE
             CALL dbcsr_data_set (dsts,&
                  ABS(blk_p(i)), sizes(i),&
                  srcs(old_blk_d(i)), source_lb=ABS(old_blk_p(i)))
             !dst(ABS(blk_p(i)):ABS(blk_p(i))+sizes(i)-1) =&
             !     srcs(old_blk_d(i))%d&
             !     %r_dp(ABS(old_blk_p(i)):ABS(old_blk_p(i))+sizes(i)-1)
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF
    DBGV 'new blk_p =',blk_p
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_sort_data


! *****************************************************************************
!> \brief Determines the sign for a block pointer.
!> \param[in] n          block offset
!> \param[in] oldsign    sign for the new block offset
!> \param[in] x          reverse the sign of the new block offset
! *****************************************************************************
  ELEMENTAL FUNCTION sgn (n, oldsign, x) RESULT (val)
    INTEGER, INTENT(IN)                      :: n, oldsign
    LOGICAL, INTENT(IN)                      :: x
    INTEGER                                  :: val

    IF (.NOT.x) THEN
       val = SIGN (n, oldsign)
    ELSE
       val = -SIGN (n, oldsign)
    ENDIF
  END FUNCTION sgn


! *****************************************************************************
!> \brief Sorts the data in a matrix so that the data blocks follow
!>        sequentially and does various transposing options.
!>
!> As opposed to dbcsr_sort_data, this routine calculates block sizes
!> \param[out] blk_p          re-arranged block pointers reflecting the new
!>                            data order
!> \param[in] old_blk_p       current block pointers
!> \param[in] rbs, cbs        sizes of the blocked rows and columns
!> \param[in] row_p, col_i    index
!> \param[out] dst            sorted data
!> \param[in] src             existing unordered data
!> \param[in] mark_transposed           mark data as transposed by negating the
!>                                      blk_p index entries
!> \param[in] transpose_blocks          transpose data blocks
! *****************************************************************************
  SUBROUTINE dbcsr_copy_sort_data(blk_p, old_blk_p, row_p, col_i, rbs, cbs,&
       dst, src, mark_transposed, transpose_blocks)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: blk_p
    INTEGER, DIMENSION(:), INTENT(IN)        :: old_blk_p, row_p, col_i, rbs, &
                                                cbs
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dst
    TYPE(dbcsr_data_obj), INTENT(IN)         :: src
    LOGICAL, INTENT(IN), OPTIONAL            :: mark_transposed, &
                                                transpose_blocks

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_sort_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col_size, nblks, nrows, &
                                                nze, nze_prev, row, row_size
    LOGICAL                                  :: mark, trb

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    ! Analyze parameters
    mark = .FALSE.
    IF (PRESENT (mark_transposed)) mark = mark_transposed
    trb = .FALSE.
    IF (PRESENT (transpose_blocks)) trb = transpose_blocks
    !
    nblks = SIZE (old_blk_p)
    nrows = SIZE (row_p)-1
    CALL cp_assert(SIZE(blk_p).GE.nblks, cp_failure_level,&
         cp_wrong_args_error, routineN, 'Destination blk_p too small.')
    CALL cp_assert(dbcsr_get_data_size(dst).GE.dbcsr_get_data_size(src),&
         cp_warning_level, cp_wrong_args_error, routineN,&
         'Destination data space may be too small.')
    IF (nblks .GE. 1) &
         blk_p(1) = SGN(1, old_blk_p(1), mark)
    nze_prev = 0
    DO row = 1, nrows
       row_size = rbs (row)
       DO blk = row_p(row)+1, row_p(row+1)
          IF (old_blk_p(blk) .NE. 0) THEN
             col_size = cbs(col_i(blk))
             nze = row_size * col_size
             IF (blk .GT. 1) THEN
                blk_p(blk) = SGN(ABS(blk_p(blk-1)) + nze_prev, old_blk_p(blk),&
                     mark)
             ENDIF
             IF (.NOT. trb) THEN
                CALL dbcsr_data_set (dst, ABS(blk_p(blk)), nze,&
                     src, source_lb=ABS(old_blk_p(blk)))
                !dst(ABS(blk_p(blk)):ABS(blk_p(blk))+nze-1) =&
                !     src(ABS(old_blk_p(blk)):ABS(old_blk_p(blk))+nze-1)
             ELSE
                CALL dbcsr_block_transpose (dst, src, ABS(blk_p(blk)),&
                     !row_size, col_size, ABS(old_blk_p(blk)))
                     col_size, row_size, ABS(old_blk_p(blk)))
                !CALL dbcsr_block_transpose (&
                !     dst(ABS(blk_p(blk)):ABS(blk_p(blk))+nze-1),&
                !     src(ABS(old_blk_p(blk)):ABS(old_blk_p(blk))+nze-1),&
                !     row_size, col_size)
             ENDIF
             nze_prev = nze
          ENDIF ! blk exists
       ENDDO ! blk
    ENDDO ! row
  END SUBROUTINE dbcsr_copy_sort_data


! *****************************************************************************
!> \brief Merge data from sorted work matrices indices into the final matrix.
!>
!> \param[in,out] matrix      matrix to work on
!> \param[in] sort_data       (optional) whether to sort data or copy/move
! *****************************************************************************
  SUBROUTINE dbcsr_merge_data(matrix, sort_data)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: sort_data

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_merge_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_p, i, nwms, nze, &
                                                o_blk_p, prev_row, row, &
                                                smallest, smallest_col, &
                                                smallest_row, treesize
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: keys
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: wm_blk_p
    INTEGER, DIMENSION(:), POINTER           :: cbs, rbs
    LOGICAL                                  :: srt
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: finished, tree_based
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_2d
    TYPE(btree_2d_data_d), ALLOCATABLE, &
      DIMENSION(:)                           :: values_d

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    srt = .FALSE.
    IF (PRESENT (sort_data)) srt = sort_data
    rbs => array_data (matrix%row_blk_size)
    cbs => array_data (matrix%col_blk_size)
    nwms = SIZE (matrix%wms)
    ALLOCATE (wm_blk_p(nwms), finished(nwms), tree_based(nwms))
    wm_blk_p(:) = 1
    ! Convert mutable data type into standard work arrays
    DO i = 1, nwms
       tree_based(i) = dbcsr_wm_use_mutable (matrix%wms(i))
       IF (tree_based(i)) THEN
          IF (dbg) WRITE(*,'(1X,A,I5,A)')routineN//" Work matrix",i,&
               " has mutable data"
          ! srt = .TRUE. ! Not needed because of the copy
          treesize = btree_size_d(matrix%wms(i)%mutable%m%btree_d)
          IF (dbg) WRITE(*,*)routineN//" Declared size",matrix%wms(i)%lastblk,&
               "present blocks",treesize
          CALL cp_assert(matrix%wms(i)%lastblk .EQ. treesize,&
               cp_fatal_level, cp_internal_error, routineN,&
               "Mismatch in number of blocks")
          ALLOCATE (keys(treesize), values_d(treesize))
          CALL btree_destroy_d (matrix%wms(i)%mutable%m%btree_d, keys, values_d)
          CALL ensure_array_size (matrix%wms(i)%row_i, ub=treesize)
          CALL ensure_array_size (matrix%wms(i)%col_i, ub=treesize)
          CALL dbcsr_unpack_i8_2i4 (keys, matrix%wms(i)%row_i,&
               matrix%wms(i)%col_i)
          ! For now we also fill the data, sloooowly, but this should
          ! be avoided and the data should be copied directly from the
          ! source in the subroutine's main loop.
          CALL ensure_array_size (matrix%wms(i)%blk_p, ub=treesize)
          CALL dbcsr_data_ensure_size (matrix%wms(i)%data_area,&
               matrix%wms(i)%datasize)
          blk_p = 1
          DO blk = 1, treesize
             block_2d => values_d(blk)%p
             !write(*,*)"block", blk
             !write(*,'(10(1X,F5.1))')block_2d

             IF (.NOT. values_d(blk)%tr) THEN
                matrix%wms(i)%blk_p(blk) = blk_p
             ELSE
                matrix%wms(i)%blk_p(blk) = -blk_p
             ENDIF

             CALL dbcsr_data_set (matrix%wms(i)%data_area, blk_p,&
                  SIZE(block_2d),&
                  RESHAPE (block_2d, (/SIZE(block_2d)/)), source_lb=1)
             blk_p = blk_p + SIZE(block_2d)
             DEALLOCATE (block_2d)
          ENDDO
          DEALLOCATE (keys, values_d)
          CALL dbcsr_mutable_release (matrix%wms(i)%mutable)
          !write(*,*)"final data",dbcsr_get_data_size (matrix%wms(i)%data_area),&
          !     treesize
          !write(*,'(10(1X,F5.1))')&
          !dbcsr_get_data_p (matrix%wms(i)%data_area, 0.0_dp)
       ENDIF
    ENDDO
    FORALL (i = 1:nwms)
       finished(i) = wm_blk_p(i) .GT. matrix%wms(i)%lastblk
    END FORALL
    IF (nwms .GT. 1) srt = .TRUE.
    IF (matrix%special_memory) srt = .TRUE.
    IF (.NOT. srt) THEN
       CALL dbcsr_data_release (matrix%data_area)
       matrix%data_area = matrix%wms(1)%data_area
       CALL dbcsr_data_hold (matrix%data_area)
    ELSE
       CALL dbcsr_data_ensure_size (matrix%data_area, matrix%nze)
    ENDIF
    IF (dbg) THEN
       DO i = 1, nwms
          WRITE(*,'(A,1X,I4,1X,A,3(1X,I7))')routineN//" size of wm",i,"is",&
               matrix%wms(i)%lastblk,&
               matrix%wms(i)%datasize,&
               dbcsr_get_data_size (matrix%wms(i)%data_area)
          WRITE(*,'(A)')routineP//' row_i='
          IF (ASSOCIATED (matrix%wms(i)%row_i))&
               WRITE(*,'(10(1X,I7))')matrix%wms(i)%row_i(:)
          WRITE(*,'(A)')routineP//' col_i='
          IF (ASSOCIATED (matrix%wms(i)%col_i))&
               WRITE(*,'(10(1X,I7))')matrix%wms(i)%col_i(:)
          WRITE(*,'(A)')routineP//' blk_p='
          IF (ASSOCIATED (matrix%wms(i)%blk_p))&
               WRITE(*,'(10(1X,I7))')matrix%wms(i)%blk_p(:)
       ENDDO
    ENDIF
    matrix%row_p(:) = -1
    matrix%row_p(1) = 0
    blk_p = 1
    blk = 1
    row = 1
    prev_row = 1
    smallest_row = 1
    DO WHILE (.NOT. ALL (finished))
       ! Find the (still-valid) source block that has the smallest column in
       ! the smallest row. Start with any valid source block.
       smallest = MINLOC (wm_blk_p, 1, .NOT.finished)
       smallest_row = matrix%wms(smallest)%row_i(wm_blk_p(smallest))
       smallest_col = matrix%wms(smallest)%col_i(wm_blk_p(smallest))
       DO i = 1, nwms
          IF (.NOT. finished(i)) THEN
             IF (matrix%wms(i)%row_i(wm_blk_p(i)) .LE. smallest_row) THEN
                IF (matrix%wms(i)%row_i(wm_blk_p(i)) .LT. smallest_row) THEN
                   smallest_row = matrix%wms(i)%row_i(wm_blk_p(i))
                   smallest_col = matrix%wms(i)%col_i(wm_blk_p(i))
                   smallest = i
                ELSE
                   IF (matrix%wms(i)%col_i(wm_blk_p(i)) .LT. smallest_col) THEN
                      smallest_col = matrix%wms(i)%col_i(wm_blk_p(i))
                      smallest = i
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
       ENDDO
       matrix%col_i(blk) = smallest_col
       o_blk_p = matrix%wms(smallest)%blk_p(wm_blk_p(smallest))
       nze = rbs(smallest_row) * cbs(smallest_col)
       IF (dbg) THEN
          WRITE(*,'(A,"(",I5,",",I5,") at ",I7,1X,A,I5,1X,I7,1X,A,1X,I7)')&
               routineN//" Adding",&
               smallest_row, smallest_col, blk_p, "from",smallest, o_blk_p,&
               "size", nze
       END IF
       IF (srt) THEN
          matrix%blk_p(blk) = SIGN (blk_p, o_blk_p)
          CALL dbcsr_data_set (matrix%data_area, ABS(blk_p), nze,&
               matrix%wms(smallest)%data_area, ABS(o_blk_p))
       ENDIF
       IF (smallest_row .GT. prev_row) THEN
          matrix%row_p(prev_row+1:smallest_row) = blk-1
          IF (dbg) WRITE(*,'(15(1X,I4))')matrix%row_p(:)
       ENDIF
       prev_row = smallest_row
       blk = blk + 1
       blk_p = blk_p + nze
       wm_blk_p(smallest) = wm_blk_p(smallest)+1
       finished(smallest) = wm_blk_p(smallest) .GT. matrix%wms(smallest)%lastblk
    ENDDO
    matrix%row_p(smallest_row+1:matrix%nblkrows_total+1) = blk-1
    matrix%nze = blk_p - 1
    IF (.NOT. srt .AND. ASSOCIATED (matrix%wms(1)%blk_p)) THEN
       matrix%blk_p(1:matrix%wms(1)%lastblk)&
            = matrix%wms(1)%blk_p(1:matrix%wms(1)%lastblk)
    ENDIF
    IF (dbg) WRITE(*,*)routineP//'row_p=',matrix%row_p(:)
    IF (dbg) WRITE(*,*)routineP//'col_i=',matrix%col_i(:)
    IF (dbg) WRITE(*,*)routineP//'blk_p=',matrix%blk_p(:)
  END SUBROUTINE dbcsr_merge_data


! *****************************************************************************
!> \brief Creates row and column images of a matrix.
!> \param[in] source          input matrix
!> \param[in,out] normalized  image array of the normalized matrix
!> \param[in] target_image_dist          normalize to this image distribution
!> \param[in] predistribute   (optional) predistribute data for multiplication
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_make_images(source, normalized, target_image_dist,&
       predistribute, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: source
    TYPE(dbcsr_2d_array_type), INTENT(OUT)   :: normalized
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(IN)                             :: target_image_dist
    CHARACTER, INTENT(IN), OPTIONAL          :: predistribute
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_images', &
      routineP = moduleN//':'//routineN

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL cp_assert (dbcsr_valid_index (source),&
         cp_fatal_level, cp_caller_error, routineN,&
         "Matrix not initialized.")
    DBG 'Make images for matrix', source%m%matrix_type, source%m%name
    normalized%image_dist = target_image_dist
    SELECT CASE (source%m%matrix_type)
    CASE ('N', 'S', 'A')
       CALL make_images(source, normalized,&
            target_image_dist, desymmetrize=.TRUE.,&
            predistribute=predistribute, error=error)
    CASE ('T')
       CALL cp_assert(.FALSE., cp_failure_level,&
            cp_unimplemented_error_nr, routineP,&
            "Making images of transposed matrices not "//&
            "supported yet.", error)
    CASE DEFAULT
       CALL cp_assert(.FALSE., cp_failure_level,&
            cp_wrong_args_error, routineP, "Invalid matrix type.", error)
    END SELECT
    normalized%image_dist = target_image_dist
    DBG 'parent ums(1) name', normalized%mats(1,1)%m%name
  END SUBROUTINE dbcsr_make_images


! *****************************************************************************
!> \brief Makes column-based and row-based images of a matrix.
!> \param[in] sm    input symmetric matrix
!> \param[in,out] ums         normalied matrices
!> \param[in] target_imgdist  image distribution to normalize to
!> \param[in] desymmetrize    (optional) desymmetrize a symmetric matrix
!> \param[in] predistribute   (optional) predistribute data for multiplication
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE make_images(ism, ums, target_imgdist, desymmetrize, predistribute,&
       error)
    TYPE(dbcsr_obj), INTENT(IN)              :: ism
    TYPE(dbcsr_2d_array_type), INTENT(OUT)   :: ums
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(IN)                             :: target_imgdist
    LOGICAL, INTENT(IN), OPTIONAL            :: desymmetrize
    CHARACTER, INTENT(IN), OPTIONAL          :: predistribute
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_images', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 5

    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: c_sp
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: c_dp
    INTEGER :: blk, blk_l, blk_p, bp, col, col_img, data_p, dst_p, i, &
      ithread, j, mp_group, ncol_images, nrow_images, nsymmetries, nthreads, &
      numproc, nze, pcol, prev_blk_p, prev_dst_p, prev_row, prow, row, &
      row_img, sd_pos, sm_pos, src_p, stored_blk_p, stored_col, stored_row, &
      symmetry_i, timing_handle, vcol, vrow
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: lsdp, lsmp, rd_disp, &
                                                recv_meta, rm_disp, sd_disp, &
                                                sdp, send_meta, sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: blk_ps, blks, &
                                                total_recv_count, &
                                                total_send_count
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: row_dst_counts, &
                                                row_dst_offsets
    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: myt_send_count, recv_count, &
                                                send_count
    INTEGER, DIMENSION(:), POINTER :: col_blk_size, col_dist, col_img_dist, &
      row_blk_size, row_dist, row_img_dist, stored_local_rows
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    LOGICAL                                  :: predist, same_dst_p, tr
    REAL(KIND=dp)                            :: tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: send_data
    REAL(KIND=dp), DIMENSION(:), POINTER     :: recv_data
    REAL(KIND=real_4), DIMENSION(:), POINTER :: r_sp
    REAL(KIND=real_8), DIMENSION(:), POINTER :: r_dp
    TYPE(dbcsr_data_obj)                     :: received_data_area
    TYPE(dbcsr_distribution_obj)             :: target_dist
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_type)                         :: sm

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    tstart = m_walltime ()
    sm = ism%m
    nsymmetries = 1
    IF (PRESENT (desymmetrize)) THEN
       IF (desymmetrize .AND. (&
            sm%matrix_type .EQ. dbcsr_type_symmetric&
            .OR. sm%matrix_type .EQ. dbcsr_type_antisymmetric)) THEN
          nsymmetries = 2
       ENDIF
    ENDIF
    predist = .FALSE.
    IF (PRESENT (predistribute)) predist = .TRUE.
    SELECT CASE (sm%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (sm%data_area, r_dp)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (sm%data_area, r_sp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (sm%data_area, c_dp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (sm%data_area, c_sp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    END SELECT
    row_blk_size => array_data (sm%row_blk_size)
    col_blk_size => array_data (sm%col_blk_size)
    target_dist = target_imgdist%main
    row_dist => array_data (dbcsr_distribution_row_dist (target_dist))
    col_dist => array_data (dbcsr_distribution_col_dist (target_dist))
    nrow_images = target_imgdist%row_decimation
    IF (nrow_images .GT. 1) THEN
       row_img_dist => array_data (target_imgdist%row_image)
    ELSE
       NULLIFY (row_img_dist)
    ENDIF
    ncol_images = target_imgdist%col_decimation
    IF (ncol_images .GT. 1) THEN
       col_img_dist => array_data (target_imgdist%col_image)
    ELSE
       NULLIFY (col_img_dist)
    ENDIF
    stored_local_rows => array_data (dbcsr_distribution_local_rows (sm%dist))
    mp_obj = dbcsr_distribution_mp (target_dist)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    numproc = dbcsr_mp_numnodes (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    DBG 'row_dist', row_dist
    DBG 'col_dist', col_dist
    DBG 'row, col maxval:', MAXVAL(row_dist), MAXVAL(col_dist)
    DBG 'blacs2mpi',blacs2mpi
    DBG 'blacs2mpi LB',LBOUND(blacs2mpi)
    DBG 'blacs2mpi UB',UBOUND(blacs2mpi)
    DBG 'row,col images',nrow_images,ncol_images
    ALLOCATE (ums%mats(nrow_images,ncol_images))
    IF (sm%matrix_type .EQ. 'S' .OR. sm%matrix_type .EQ. 'A') THEN
       CALL cp_assert(SIZE(row_dist) .EQ. SIZE(col_dist), cp_warning_level,&
            cp_wrong_args_error, routineN,&
            'Unequal row and column distributions for symmetric matrix.', error)
    ENDIF
    CALL cp_assert(MAXVAL(row_dist).LE.UBOUND(blacs2mpi,1), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows', error)
    IF (dbg) &
         CALL cp_assert(MAXVAL(row_dist).EQ.UBOUND(blacs2mpi,1), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Range of row distribution not equal to processor rows', error)
    CALL cp_assert(MAXVAL(col_dist).LE.UBOUND(blacs2mpi,2), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols', error)
    IF (dbg) &
         CALL cp_assert(MAXVAL(col_dist).EQ.UBOUND(blacs2mpi,2), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Range of col distribution not equal to processor cols', error)
    ALLOCATE (send_count(2, nrow_images, ncol_images, 0:numproc-1))
    ALLOCATE (recv_count(2, nrow_images, ncol_images, 0:numproc-1))
    ALLOCATE (total_send_count(2, 0:numproc-1))
    ALLOCATE (total_recv_count(2, 0:numproc-1))
    ALLOCATE (sdp(0:numproc-1))
    ALLOCATE (sd_disp(0:numproc-1))
    ALLOCATE (smp(0:numproc-1))
    ALLOCATE (sm_disp(0:numproc-1))
    ALLOCATE (rd_disp(0:numproc-1))
    ALLOCATE (rm_disp(0:numproc-1))
    ALLOCATE (row_dst_counts(2, 0:numproc-1, sm%nblkrows_total))
    ALLOCATE (row_dst_offsets(2, 0:numproc-1, sm%nblkrows_total))
    ithread = 0
    nthreads = 1
    ! Create and allocate the imaged matrices.
    DO row_img = 1, nrow_images
       DO col_img = 1, ncol_images
          CALL dbcsr_init (ums%mats(row_img, col_img), error)
          CALL dbcsr_create(ums%mats(row_img, col_img), "imaged "//sm%name,&
               target_dist,&
               dbcsr_type_normal, sm%row_blk_size, sm%col_blk_size,&
               0,0, sm%data_type, error=error, special=.TRUE.)
          IF (sm%matrix_type .EQ. dbcsr_type_antisymmetric) THEN
             ums%mats(row_img, col_img)%m%sym_negation = .TRUE.
          ENDIF
       ENDDO
    ENDDO
    send_count(:,:,:,:) = 0
    row_dst_counts(:,:,:) = 0
    ! C
    ALLOCATE (myt_send_count(2, nrow_images, ncol_images, 0:numproc-1))
    myt_send_count(:,:,:,:) = 0
    prev_dst_p = -1
    ! Count sizes for sending.
    CALL dbcsr_iterator_start(iter, ism)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk)
       DO symmetry_i = 1, nsymmetries
          IF (symmetry_i .EQ. 1) THEN
             stored_row = row ; stored_col = col
          ELSE
             IF (row .EQ. col) CYCLE
             stored_row = col ; stored_col = row
          ENDIF
          ! Where do we send this block?
          prow = row_dist(stored_row)
          pcol = col_dist(stored_col)
          row_img = 1
          IF (nrow_images .GT. 1) row_img = row_img_dist (stored_row)
          col_img = 1
          IF (ncol_images .GT. 1) col_img = col_img_dist (stored_col)
          IF (predist) THEN
             IF (dbg) THEN
                WRITE(*,'(1X,A,2(1X,I5))')routineN//" Block", stored_row, stored_col
             ENDIF
             vcol = pcol*target_imgdist%col_decimation + col_img - 1
             vrow = prow*target_imgdist%row_decimation + row_img - 1
             i = vrow ; j = vcol
             IF (predistribute.eq."R") THEN
                ! shift = vpcol
                ! vrow = vrow - shift
                vrow = vrow - pcol*target_imgdist%row_multiplicity
                vrow = absmod(vrow, dbcsr_mp_nprows(mp_obj)*target_imgdist%row_decimation)
                prow = vrow / target_imgdist%row_decimation
                row_img = 1+absmod (vrow, target_imgdist%row_decimation)
                IF (dbg) &
                     WRITE(*,'(2(1X,A,2(1X,I5)))')routineN//" R shift from",&
                     i,j, "to", vrow, vcol
             ELSEIF (predistribute.eq."L") THEN
                ! shift = vprow
                ! vcol = vcol - shift
                vcol = vcol - prow*target_imgdist%col_multiplicity
                vcol = absmod(vcol, dbcsr_mp_npcols(mp_obj)*target_imgdist%col_decimation)
                pcol = vcol / target_imgdist%col_decimation
                col_img = 1+absmod (vcol, target_imgdist%col_decimation)
                IF (dbg) &
                     WRITE(*,'(2(1X,A,2(1X,I5)))')routineN//" L shift from",&
                     i,j, "to", vrow, vcol
             ENDIF
          ENDIF
          dst_p = blacs2mpi(prow, pcol)
          same_dst_p = prev_dst_p .EQ. dst_p
          prev_dst_p = dst_p
          ! These counts are meant for the thread that processes this row.
          row_dst_counts(1, dst_p, row) =&
               row_dst_counts(1, dst_p, row) + metalen
          myt_send_count(1, row_img, col_img, dst_p) =&
               myt_send_count(1, row_img, col_img, dst_p) + 1
          ! Data can be duplicated if the transpose is destined to the same
          ! process.
          IF (.NOT. same_dst_p .OR. symmetry_i .EQ. 1) THEN
             nze = row_blk_size(stored_row) * col_blk_size(stored_col)
             row_dst_counts(2, dst_p, row) =&
                  row_dst_counts(2, dst_p, row) + nze
             myt_send_count(2, row_img, col_img, dst_p) =&
                  myt_send_count(2, row_img, col_img, dst_p) + nze
          ENDIF
       ENDDO ! symmetry_i
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    ! Merge the send counts.
    send_count(:,:,:,:) = send_count(:,:,:,:) + myt_send_count(:,:,:,:)
    CALL mp_alltoall(send_count, recv_count, 2*nrow_images*ncol_images,&
         mp_group)
    ! Fill in the meta data structures and copy the data.
    DO dst_p = 0, numproc-1
       total_send_count(1, dst_p) = SUM (send_count (1, :, :, dst_p))
       total_send_count(2, dst_p) = SUM (send_count (2, :, :, dst_p))
       total_recv_count(1, dst_p) = SUM (recv_count (1, :, :, dst_p))
       total_recv_count(2, dst_p) = SUM (recv_count (2, :, :, dst_p))
    ENDDO
    ! Allocate data structures needed for data exchange.
    ALLOCATE (recv_data(SUM(recv_count(2, :, :, :))))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :, :, :))))
    ALLOCATE (send_data(SUM(send_count(2, :, :, :))))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :, :, :))))
    ! Calculate displacements for processors needed for the exchanges.
    sd_disp = -1 ; sm_disp = -1
    rd_disp = -1 ; rm_disp = -1
    sd_disp(0) = 1 ; sm_disp(0) = 1
    rd_disp(0) = 1 ; rm_disp(0) = 1
    DO dst_p = 1, numproc-1
       sm_disp(dst_p) = sm_disp(dst_p-1)&
                        + metalen*total_send_count(1, dst_p-1)
       sd_disp(dst_p) = sd_disp(dst_p-1)&
                        + total_send_count(2, dst_p-1)
       rm_disp(dst_p) = rm_disp(dst_p-1)&
                        + metalen*total_recv_count(1, dst_p-1)
       rd_disp(dst_p) = rd_disp(dst_p-1)&
                        + total_recv_count(2, dst_p-1)
    ENDDO
    ! Calculate the displacements for rows on processors needed for
    ! the exchanges. (Used to avoid critical sections in the following
    ! DO loop that copies data into the send buffer.)
    row_dst_offsets(:,0,1) = 1
    DO dst_p = 1, numproc-1
       row_dst_offsets(1, dst_p, 1) = row_dst_offsets(1, dst_p-1, 1)&
            + metalen*total_send_count(1, dst_p-1)
       row_dst_offsets(2, dst_p, 1) = row_dst_offsets(2, dst_p-1, 1)&
            + total_send_count(2, dst_p-1)
    ENDDO
    DO row = 2, sm%nblkrows_total
       row_dst_offsets(:, :, row) =&
            + row_dst_offsets(:, :, row-1)&
            + row_dst_counts(:, :, row-1)
    ENDDO
    ! Prepares some indices needed for the last DO loop that copies
    ! from buffer to local space. Placed early to take advantage of
    ! the SECTIONS.
    ALLOCATE (blk_ps(nrow_images, ncol_images))
    ALLOCATE (blks (nrow_images, ncol_images))
    blk_ps(:,:) = 1
    blks(:,:) = 1
    ! Prepares the work matrices used in the last DO loop. Placed
    ! early.
    CALL cp_assert (sm%data_type .EQ. dbcsr_type_real_8, cp_fatal_level,&
         cp_unimplemented_error_nr, routineN, "Only real double support.")
    CALL dbcsr_data_init (received_data_area)
    CALL dbcsr_data_new (received_data_area, sm%data_type,&
         special_memory=.FALSE.)
    ! This is an ugly cuckoo.
    CALL dbcsr_data_set_pointer (received_data_area, recv_data)
    !received_data_area%d%r_dp => recv_data
    DO row_img = 1, nrow_images
       DO col_img = 1, ncol_images
          CALL dbcsr_work_create(ums%mats(row_img, col_img),&
               SUM(recv_count(1,row_img,col_img,:)), error=error)
          CALL dbcsr_data_hold (received_data_area)
          CALL dbcsr_data_release (ums%mats(row_img,col_img)%m%wms(1)%data_area)
          ums%mats(row_img,col_img)%m%wms(1)%data_area = received_data_area
       ENDDO
    ENDDO
    ALLOCATE (lsmp(0:numproc-1), lsdp(0:numproc-1))
    prev_dst_p = -1
    ! Copies metadata and actual data to be sent into the send buffers.
    CALL dbcsr_iterator_start(iter, ism)
    prev_row = 0
    prev_blk_p = 0
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk, blk_p=blk_p)
       ! New row
       IF (row .NE. prev_row) THEN
          lsmp = row_dst_offsets(1, :, row)
          lsdp = row_dst_offsets(2, :, row)
          prev_row = row
       ENDIF
       bp = ABS(blk_p)
       DO symmetry_i = 1, nsymmetries
          IF (symmetry_i .EQ. 1) THEN
             stored_row = row ; stored_col = col; tr = blk_p .LT. 0
          ELSE                                                
             IF (row .EQ. col) CYCLE
             stored_row = col ; stored_col = row; tr = blk_p .GT. 0
          ENDIF
          ! Where do we send this block?
          prow = row_dist(stored_row)
          pcol = col_dist(stored_col)
          row_img = 1
          IF (nrow_images .GT. 1) row_img = row_img_dist (stored_row)
          col_img = 1
          IF (ncol_images .GT. 1) col_img = col_img_dist (stored_col)
          IF (predist) THEN
             vcol = pcol*target_imgdist%col_decimation + col_img - 1
             vrow = prow*target_imgdist%row_decimation + row_img - 1
             IF (predistribute.eq."R") THEN
                ! shift = vpcol
                ! vrow = vrow + shift
                vrow = vrow - pcol*target_imgdist%row_multiplicity
                vrow = absmod(vrow, dbcsr_mp_nprows(mp_obj)*target_imgdist%row_decimation)
                prow = vrow / target_imgdist%row_decimation
                row_img = 1+absmod (vrow, target_imgdist%row_decimation)
             ELSEIF (predistribute.eq."L") THEN
                ! shift = vprow
                ! vcol = vcol + shift
                vcol = vcol - prow*target_imgdist%col_multiplicity
                vcol = absmod(vcol, dbcsr_mp_npcols(mp_obj)*target_imgdist%col_decimation)
                pcol = vcol / target_imgdist%col_decimation
                col_img = 1+absmod (vcol, target_imgdist%col_decimation)
             ENDIF
          ENDIF
          dst_p = blacs2mpi(prow, pcol)
          same_dst_p = dst_p .EQ. prev_dst_p
          prev_dst_p = dst_p
          sm_pos = lsmp(dst_p)
          lsmp(dst_p) = lsmp(dst_p) + metalen
          send_meta(sm_pos) = stored_row
          IF (.NOT. same_dst_p .OR. symmetry_i .EQ. 1) THEN
             nze = row_blk_size(stored_row) * col_blk_size(stored_col)
             sd_pos = lsdp(dst_p)
             lsdp(dst_p) = lsdp(dst_p) + nze
             send_data(sd_pos:sd_pos+nze-1) = r_dp(bp:bp+nze-1)
             send_meta(sm_pos+1) = stored_col
             send_meta(sm_pos+2) = SGN (&
                  sd_pos-sd_disp(dst_p)+1, blk_p, tr)
             send_meta(sm_pos+3) = row_img
             send_meta(sm_pos+4) = col_img
             prev_blk_p = send_meta(sm_pos+2)
          ELSE
             send_meta(sm_pos+1) = -stored_col
             send_meta(sm_pos+2) = -prev_blk_p
             send_meta(sm_pos+3) = row_img
             send_meta(sm_pos+4) = col_img
          ENDIF
       ENDDO ! symmetry_i
    ENDDO ! iterator
    CALL dbcsr_iterator_stop(iter)
    ! Exchange the data and metadata structures.
    CALL mp_alltoall(send_data(:), total_send_count(2,:), sd_disp(:)-1,&
         recv_data(:), total_recv_count(2,:), rd_disp(:)-1, mp_group)
    CALL mp_alltoall(send_meta(:), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), metalen*total_recv_count(1,:), rm_disp(:)-1, mp_group)
    ! Now create the work index and/or copy the relevant data from the
    ! receive buffer into the local indices.
    prev_blk_p = 0
    DO src_p = 0, numproc-1
       data_p = 0
       DO blk_l = 1, total_recv_count(1, src_p)
          stored_row = recv_meta(rm_disp(src_p)+metalen*(blk_l-1))
          stored_col = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+1)
          stored_blk_p = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+2)
          row_img = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+3)
          col_img = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+4)
          nze = row_blk_size(ABS(stored_row))&
               * col_blk_size(ABS(stored_col))
          blk = blks(row_img,col_img)
          blks(row_img,col_img) = blks(row_img,col_img) + 1
          blk_ps(row_img,col_img) = blk_ps(row_img,col_img) + nze
          IF (stored_col .GT. 0) THEN
             blk_p = data_p
             data_p = data_p + nze
          ELSE
             blk_p = prev_blk_p
          ENDIF
          ums%mats(row_img,col_img)%m%wms(1)%row_i(blk) = ABS(stored_row)
          ums%mats(row_img,col_img)%m%wms(1)%col_i(blk) = ABS(stored_col)
          ums%mats(row_img,col_img)%m%wms(1)%blk_p(blk) =&
               SIGN(rd_disp(src_p) + ABS(stored_blk_p)-1, stored_blk_p)
          prev_blk_p = blk_p
       ENDDO
    ENDDO
    ! Finalize the actual imaged matrices from the work matrices.
    DO row_img = 1, nrow_images
       DO col_img = 1, ncol_images
          ums%mats(row_img,col_img)%m%wms(1)%lastblk = blks(row_img,col_img) - 1
          ums%mats(row_img,col_img)%m%wms(1)%datasize = blk_ps(row_img,col_img) - 1
          IF (nrow_images.EQ.1 .AND. ncol_images.eq.1) THEN
             CALL dbcsr_finalize(ums%mats(row_img,col_img),reshuffle=.FALSE.)
          ELSE
             CALL dbcsr_finalize(ums%mats(row_img,col_img),reshuffle=.TRUE.)
          ENDIF
       ENDDO
    ENDDO
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_data)
    DEALLOCATE(send_meta)
    ! Get rid of the cuckoo.
    CALL dbcsr_data_release (received_data_area)
    tstop = m_walltime ()
    !WRITE(*,'(1X,A,1X,EN12.4)')'make_column_images time:',tstop-tstart
    DBG 'ums(1) name',ums%mats(1,1)%m%name
    DBG 'Done making images'
    CALL timestop(timing_handle)
  END SUBROUTINE make_images


! *****************************************************************************
!> \brief A modulo operation that always returns non-negative numbers.
!> \param[in] n               take modulo of this number
!> \param[in] divisor         modulo divisor
!> \result absmod             resulting remainder
! *****************************************************************************
  pure FUNCTION absmod(n, divisor)
    INTEGER, INTENT(IN)                      :: n, divisor
    INTEGER                                  :: absmod

    absmod = MOD(n+divisor, divisor)
  END FUNCTION absmod


! *****************************************************************************
!> \brief Copies an array into another array, ensuring that the target array
!>        is large enough.
!> \par Makes use of FORTRAN's array bounds for array sections.
!> \param to        array to copy to (target array)
!> \param[in] from  array to copy from (source array)
!> \param[in] offset          starting position in the target array
!> \param[in] factor          factor for size increase
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE addto_array_i(to, from, offset, factor, special, error)
    INTEGER, DIMENSION(:), POINTER           :: to
    INTEGER, DIMENSION(:), INTENT(IN)        :: from
    INTEGER, INTENT(IN), OPTIONAL            :: offset
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: special
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'addto_array_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: last, lb_t, s_f, s_t, start, &
                                                timing_handle, ub_t
    LOGICAL                                  :: spec
    REAL                                     :: f

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    spec = .FALSE.
    IF (PRESENT (special)) spec = special
    s_f = SIZE(from, 1)
    ub_t = UBOUND(to, 1)
    lb_t = LBOUND(to, 1)
    s_t = SIZE(to)
    start = lb_t
    IF (PRESENT (offset)) start = offset
    ! assert offset > lb_t
    last = start+s_f-1
    IF(dbg) WRITE(*,*)routineP//' target array bounds',lb_t,':',ub_t,&
         'starting at',start,'last',last,'source size', s_f
    IF (ub_t .LT. last) THEN
       IF (PRESENT (factor)) THEN
          f = factor
       ELSE
          f = 1.4142136
       ENDIF
       ub_t = MAX(start + INT(s_t * f), last)
       IF(dbg) WRITE(*,*)routineP//' expanding array from :',UBOUND(to),&
            'to',ub_t
       !CALL reallocate(to, lb_t, ub_t)
       CALL ensure_array_size(to, lb_t, ub_t, factor=1.0, nocopy=.FALSE.,&
            special=spec, error=error)
    ENDIF
    to(start:last) = from(:)
    CALL timestop(timing_handle)
  END SUBROUTINE addto_array_i


! *****************************************************************************
!> \brief Copies an array into another array, ensuring that the target array
!>        large enough.
!> \par Makes use of FORTRAN's array bounds for array sections.
!> \param to        array to copy to
!> \param[in] from  array to copy from
!> \param[in] offset          starting position in the target array
!> \param[in] factor          factor for size increase
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE addto_array_d(to, from, offset, factor, special, error)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: to
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: from
    INTEGER, INTENT(IN), OPTIONAL            :: offset
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: special
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'addto_array_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: last, lb_t, s_f, s_t, start, &
                                                timing_handle, ub_t
    LOGICAL                                  :: dbg, spec
    REAL                                     :: f

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    CALL timeset(routineN, timing_handle)
    spec = .FALSE.
    IF (PRESENT (special)) spec = special
    s_f = SIZE(from, 1)
    IF (dbg) THEN
       WRITE(*,*)routineP//' from size',SIZE(from),';', SIZE(from,1)
       WRITE(*,*)routineP//' from lbound',LBOUND(from),';', LBOUND(from,1)
       WRITE(*,*)routineP//' from ubound',UBOUND(from),';', UBOUND(from,1)
    ENDIF
    ub_t = UBOUND(to, 1)
    lb_t = LBOUND(to, 1)
    s_t = SIZE(to)
    start = lb_t
    IF (PRESENT (offset)) start = offset
    ! assert offset > lb_t
    last = start+s_f-1
    IF(dbg) WRITE(*,*)routineP//' target array bounds',lb_t,':',ub_t,&
         'starting at',start,'last',last,'source size', s_f
    IF (ub_t .LT. last) THEN
       IF (PRESENT (factor)) THEN
          f = factor
       ELSE
          f = 1.4142136
       ENDIF
       !ub_t = MAX(lb_t + INT(s_t * f), lb_t+s_f-1)
       ub_t = MAX(start + INT(s_t * f), last)
       IF(dbg) WRITE(*,*)routineP//' expanding array from :',UBOUND(to),&
            'to',ub_t
       !CALL reallocate(to, lb_t, ub_t)
       CALL ensure_array_size(to, lb_t, ub_t, factor=1.0, nocopy=.FALSE.,&
            special=spec)
    ENDIF
    to(start:last) = from(:)
    CALL timestop(timing_handle)
  END SUBROUTINE addto_array_d


! *****************************************************************************
!> \brief Ensures that an array is appropriately large.
!> \param[in,out] array       array to verify and possibly resize
!> \param[in] lb    (optional) desired array lower bound
!> \param[in] ub    desired array upper bound
!> \param[in] factor          (optional) factor by which to exagerrate
!>                            enlargements
!> \param[in] nocopy          (optional) copy array on enlargement; default
!>                            is to copy
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE ensure_array_size_i(array, lb, ub, factor, nocopy, special, error)
    INTEGER, DIMENSION(:), POINTER           :: array
    INTEGER, INTENT(IN), OPTIONAL            :: lb
    INTEGER, INTENT(IN)                      :: ub
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, special
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ensure_array_size_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lb_new, lb_orig, stat, &
                                                ub_new, ub_orig
    INTEGER, DIMENSION(:), POINTER           :: newarray
    LOGICAL                                  :: dbg, docopy = .TRUE., &
                                                failure = .TRUE., info, spec

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    info = .FALSE.
    CALL dbcsr_set_debug(dbg, info=info)

    docopy = .TRUE.
    IF (PRESENT (nocopy)) docopy = .NOT. nocopy
    spec = .FALSE.
    IF (PRESENT (special)) spec = special
    lb_new = 1
    IF (PRESENT (lb)) lb_new = lb
    !> Creates a new array if it doesn't yet exist.
    IF (.NOT.ASSOCIATED(array)) THEN
       !CALL timeset(routineN, timing_handle)
       failure = .FALSE.
       !CPPrecondition(.NOT.docopy, cp_warning_level, routineP, error, failure)
       failure = .TRUE.
       IF (spec) THEN
          CALL cp_assert (lb_new.EQ.1, cp_warning_level,&
               cp_internal_error,&
               routineN, "Special memory can only start at 1.")
          CALL mp_allocate(array, ub-lb_new+1, stat=stat)
       ELSE
          ALLOCATE(array(lb_new:ub), stat=stat)
       ENDIF
       IF (stat /= 0) CALL stop_memory(routineP, 'array', ub-lb_new+1)
       !CALL timestop(timing_handle)
       RETURN
    ENDIF
    lb_orig = LBOUND(array,1)
    ub_orig = UBOUND(array,1)
    IF(dbg) WRITE(*,*)routineP//' Current bounds are',lb_orig,':',ub_orig
    ! The existing array is big enough.
    IF (lb_orig.LE.lb_new .AND. ub_orig.GE.ub) THEN
       RETURN
    ENDIF
    ! A reallocation must be performed.
    !CALL timeset(routineN, timing_handle)
    IF (lb_orig.GT.lb_new) THEN
       IF (PRESENT(factor)) THEN
          lb_new = MIN(lb_orig, lb_new - NINT(REAL(lb_orig - lb_new)*factor))
       ENDIF
    ENDIF
    IF (ub_orig.LT.ub) THEN
       IF (PRESENT(factor)) THEN
          ub_new = MAX(ub_orig, ub + NINT(REAL(ub-ub_orig)*factor))
       ELSE
          ub_new = ub
       ENDIF
    ELSE
       ub_new = ub
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Resizing to bounds',lb_new,':',ub_new
    IF(.NOT.docopy) THEN
       IF (spec) THEN
          CALL mp_deallocate(array)
       ELSE
          DEALLOCATE(array)
       ENDIF
    ENDIF
    IF (spec) THEN
       CALL cp_assert (lb_new.EQ.1,&
            cp_warning_level, cp_internal_error,&
            routineN, "Special memory can only start at 1.")
       CALL mp_allocate(newarray, ub_new-lb_new+1, stat=stat)
    ELSE
       ALLOCATE(newarray(lb_new:ub_new), stat=stat)
    ENDIF
    IF(stat /= 0) CALL stop_memory(routineP, 'new', ub_new-lb_new+1)
    IF(docopy) THEN
       IF(dbg) CALL cp_assert(lb_new.LE.lb_orig .AND. ub_new.GE.ub_orig,&
            cp_failure_level, cp_internal_error, routineP,&
            "Old extent exceeds the new one.")
       newarray(lb_orig:ub_orig) = array(lb_orig:ub_orig)
       IF (spec) THEN
          CALL mp_deallocate(array)
       ELSE
          DEALLOCATE(array)
       ENDIF
    END IF
    array => newarray
    !CALL timestop(timing_handle)
  END SUBROUTINE ensure_array_size_i


! *****************************************************************************
!> \brief Ensures that an array is appropriately large.
!> \param[in,out] array       array to verify and possibly resize
!> \param[in] lb    (optional) desired array lower bound
!> \param[in] ub    desired array upper bound
!> \param[in] factor          (optional) factor by which to exagerrate
!>                            enlargements
!> \param[in] nocopy          (optional) copy array on enlargement; default
!>                            is to copy
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE ensure_array_size_i_2d(array, dim, lb, ub, factor, nocopy,&
       special, error)
    INTEGER, DIMENSION(:, :), POINTER        :: array
    INTEGER, INTENT(IN)                      :: dim
    INTEGER, INTENT(IN), OPTIONAL            :: lb
    INTEGER, INTENT(IN)                      :: ub
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, special
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ensure_array_size_i_2d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lb_new, lb_orig, osize, stat, &
                                                timing_handle, ub_new, ub_orig
    INTEGER, DIMENSION(:, :), POINTER        :: newarray
    LOGICAL                                  :: dbg, docopy = .TRUE., &
                                                failure = .TRUE., info

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    info = .FALSE.
    CALL dbcsr_set_debug(dbg, info=info)
    docopy = .TRUE.
    IF (PRESENT (nocopy)) docopy = .NOT. nocopy
    lb_new = 1
    IF (PRESENT (lb)) lb_new = lb
    !> Creates a new array if it doesn't yet exist.
    CALL cp_assert(ASSOCIATED(array), cp_failure_level, cp_internal_error, routineP,&
         'cannot resize empty multi-dimensional array', error=error, failure=failure)
    lb_orig = LBOUND(array,dim)
    ub_orig = UBOUND(array,dim)
    IF(dbg) WRITE(*,*)routineP//' Current bounds are',lb_orig,':',ub_orig
    ! The existing array is big enough.
    IF (lb_orig.LE.lb_new .AND. ub_orig.GE.ub) THEN
       RETURN
    ENDIF
    ! A reallocation must be performed.
    CALL timeset(routineN,timing_handle)
    IF (lb_orig.GT.lb_new) THEN
       IF (PRESENT(factor)) THEN
          lb_new = MIN(lb_orig, lb_new - NINT(REAL(lb_orig - lb_new)*factor))
       ENDIF
    ENDIF
    IF (ub_orig.LT.ub) THEN
       IF (PRESENT(factor)) THEN
          ub_new = MAX(ub_orig, ub + NINT(REAL(ub-ub_orig)*factor))
       ELSE
          ub_new = ub
       ENDIF
    ELSE
       ub_new = ub
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Resizing to bounds',lb_new,':',ub_new
    IF(.NOT.docopy) DEALLOCATE(array)
    CALL cp_assert(dim.EQ.1, cp_failure_level, cp_unimplemented_error_nr,&
         routineP, 'Only dimension 1 can be specified.', error=error,&
         failure=failure)
    osize = SIZE(array,2)
    ALLOCATE(newarray(lb_new:ub_new,osize), stat=stat)
    IF(stat /= 0) CALL stop_memory(routineP, 'new', ub_new-lb_new+1)
    IF(docopy) THEN
       IF(dbg) CALL cp_assert(lb_new.LE.lb_orig .AND. ub_new.GE.ub_orig,&
            cp_failure_level, cp_internal_error, routineP,&
            "Old extent exceeds the new one.")
       newarray(lb_orig:ub_orig,:) = array(lb_orig:ub_orig,:)
       DEALLOCATE(array)
    END IF
    array => newarray
    CALL timestop(timing_handle)
  END SUBROUTINE ensure_array_size_i_2d


! *****************************************************************************
!> \brief Adds a coordinate (or other data) into a work matrix's row_i and
!>        col_i arrays and returns its position.
!> \note  Uses the matrix%lastblk to keep track of the current position.
!> \param[in,out] matrix      work matrix
!> \param[in] row,col         row, col data to add
!> \param[in] blk   (optional) block pointer to add
!> \param[out] index          (optional) saved position
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE add_work_coordinate(matrix, row, col, blk, index, error)
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    INTEGER, INTENT(IN), OPTIONAL            :: blk
    INTEGER, INTENT(OUT), OPTIONAL           :: index
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_work_coordinate', &
      routineP = moduleN//':'//routineN

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    dbg = .FALSE.
    !CALL timeset(routineN, timing_handle)
    matrix%lastblk = matrix%lastblk+1
    CALL ensure_array_size(matrix%row_i, ub=matrix%lastblk)
    CALL ensure_array_size(matrix%col_i, ub=matrix%lastblk)
    matrix%row_i(matrix%lastblk) = row
    matrix%col_i(matrix%lastblk) = col
    IF (PRESENT(blk)) THEN
       CALL ensure_array_size(matrix%blk_p, ub=matrix%lastblk)
       matrix%blk_p(matrix%lastblk) = blk
    ENDIF
    IF(dbg.AND.PRESENT(blk))&
         WRITE(*,*)routineP//' Adding',row,col,blk,'at',matrix%lastblk
    IF (dbg.AND.bcsr_verbose) THEN
       WRITE(*,*)routineP//' row_i=',matrix%row_i(1:matrix%lastblk)
       WRITE(*,*)routineP//' col_i=',matrix%col_i(1:matrix%lastblk)
       WRITE(*,*)routineP//' blk_p=',matrix%blk_p(1:matrix%lastblk)
    ENDIF
    IF (PRESENT (index)) index = matrix%lastblk
    !CALL timestop(timing_handle)
  END SUBROUTINE add_work_coordinate


!! *****************************************************************************
!!> \brief Ensures that the data array is appropriately large.
!!> \param[in,out] matrix      ensure data size of this matrix
!!> \param[in] size            desired data size
!!> \param[in] factor          (optional) factor by which to exagerrate
!!>                            enlargements, default is no exagerration
!!> \param[in] nocopy          (optional) copy array on enlargement; default
!!>                            is to copy
!!> \param error     cp2k error
!! *****************************************************************************
!  SUBROUTINE dbcsr_ensure_data_size (matrix, size, factor, nocopy, error)
!    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
!    INTEGER, INTENT(IN)                      :: size
!    REAL, INTENT(IN), OPTIONAL               :: factor
!    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy
!    TYPE(cp_error_type), INTENT(INOUT)       :: error
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_ensure_data_size', &
!      routineP = moduleN//':'//routineN
!
!!   ---------------------------------------------------------------------------
!
!    CALL ensure_array_size (matrix%data, ub=size, factor=factor, nocopy=nocopy,&
!         special=matrix%special_memory, error=error)
!    SELECT CASE (matrix%data_type)
!    CASE (dbcsr_type_real_8)
!       CALL ensure_array_size (matrix%data_area%d%r_dp,&
!            ub=size, factor=factor, nocopy=nocopy,&
!            special=matrix%data_area%d%special, error=error)
!    END SELECT
!  END SUBROUTINE dbcsr_ensure_data_size


! *****************************************************************************
!> \brief Ensures a minimum size of a previously-setup data area.
!>
!> The data area must have been previously setup with dbcsr_data_new.
!> \param[inout] area         data area
!> \param[in] data_size       (optional) allocate this much data
!> \param[in] nocopy          (optional) do not keep potentially existing data,
!>                            default is to keep it
!> \param[in] zero_pad        (optional) pad new data with zeros
! *****************************************************************************
  SUBROUTINE dbcsr_data_ensure_size (area, data_size, nocopy, zero_pad)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, INTENT(IN), OPTIONAL            :: data_size
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, zero_pad

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_ensure_size', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: current_size
    LOGICAL                                  :: nocp, pad

!   ---------------------------------------------------------------------------

    CALL cp_assert(ASSOCIATED (area%d), cp_fatal_level, cp_caller_error,&
         routineN, "Data area must be setup.")
    current_size = dbcsr_get_data_size (area)
    IF (current_size .GT. data_size) RETURN
    !
    nocp = .FALSE.
    IF (PRESENT (nocopy)) nocp = nocopy
    pad = .FALSE.
    IF (PRESENT (zero_pad)) pad = zero_pad
    !
    SELECT CASE (area%d%data_type)
    CASE (dbcsr_type_real_8)
       IF (.NOT. ASSOCIATED (area%d%r_dp)) THEN
          IF (area%d%special) THEN
             CALL mp_allocate (area%d%r_dp, data_size)
          ELSE
             ALLOCATE (area%d%r_dp(data_size))
          ENDIF
          IF (pad) area%d%r_dp(:) = 0.0_dp
       ELSE
          CALL ensure_array_size (area%d%r_dp, ub=data_size,&
               special=area%d%special, nocopy=nocp, zero_pad=zero_pad)
       ENDIF
    CASE (dbcsr_type_real_4)
       IF (.NOT. ASSOCIATED (area%d%r_sp)) THEN
          IF (area%d%special) THEN
             CALL mp_allocate (area%d%r_sp, data_size)
          ELSE
             ALLOCATE (area%d%r_sp(data_size))
          ENDIF
          IF (pad) area%d%r_sp(:) = 0.0_sp
       ELSE
          CALL ensure_array_size (area%d%r_sp, ub=data_size,&
               special=area%d%special, nocopy=nocp, zero_pad=zero_pad)
       ENDIF
    CASE (dbcsr_type_complex_8)
       IF (.NOT. ASSOCIATED (area%d%c_dp)) THEN
          IF (area%d%special) THEN
             CALL mp_allocate (area%d%c_dp, data_size)
          ELSE
             ALLOCATE (area%d%c_dp(data_size))
          ENDIF
          IF (pad) area%d%c_dp(:) = 0.0_dp
       ELSE
          CALL ensure_array_size (area%d%c_dp, ub=data_size,&
               special=area%d%special, nocopy=nocp, zero_pad=zero_pad)
       ENDIF
    CASE (dbcsr_type_complex_4)
       IF (.NOT. ASSOCIATED (area%d%c_sp)) THEN
          IF (area%d%special) THEN
             CALL mp_allocate (area%d%c_sp, data_size)
          ELSE
             ALLOCATE (area%d%c_sp(data_size))
          ENDIF
          IF (pad) area%d%c_sp(:) = 0.0_sp
       ELSE
          CALL ensure_array_size (area%d%c_sp, ub=data_size,&
               special=area%d%special, nocopy=nocp, zero_pad=zero_pad)
       ENDIF
    CASE default
       CALL cp_unimplemented_error (routineN,&
            "Currently only double precision reals are supported")
    END SELECT
  END SUBROUTINE dbcsr_data_ensure_size

! *****************************************************************************
!> \brief Prints the elements of a matrix.
!> \param[in] rows,columns    is the logical (possibly detransposed) matrix
!>                            size, not the stored size.
!> \param[in] tr    (optional) specifies whether the elements are stored
!>                  transposed
! *****************************************************************************
  SUBROUTINE printmat_s(matrix, rows, cols, title, tr)
    REAL(KIND=real_4), DIMENSION(:), &
      INTENT(IN)                             :: matrix
    INTEGER, INTENT(IN)                      :: rows, cols
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

    CHARACTER(30)                            :: f
    INTEGER                                  :: m, n, r
    LOGICAL                                  :: t
    REAL(KIND=dp)                            :: bit_bucket

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    !IF (.NOT.dbg.OR..NOT.bcsr_verbose) RETURN
    m = rows
    n = cols
    t = .FALSE.
    IF(PRESENT(title)) WRITE(*,*)title
    IF(PRESENT(tr)) THEN
       IF(tr) THEN
          t = .TRUE.
          m = cols
          n = rows
       ENDIF
    ENDIF
    IF (dbg) WRITE(*,*)'Print Matrix 1D:',rows,cols,t,'range',&
         LBOUND(matrix),':',UBOUND(matrix)
    DO r = LBOUND(matrix,1), UBOUND(matrix,1)
       bit_bucket = matrix(r)
    ENDDO
    IF (dbg) WRITE(*,*)'Range OK'
    bit_bucket = 0.0_dp
    DO r = LBOUND(matrix,1), UBOUND(matrix,1)
       !IF (dbg.AND.bcsr_verbose) WRITE(*,*)':@',r,'=',matrix(r)
       bit_bucket = bit_bucket + matrix(r)
    ENDDO
    IF(m.GT.10000) m=0
    IF(n.GT.10000) n=0
    IF(m*n.LT.1 .OR. m*n.GT.SIZE(matrix)) RETURN
    !WRITE(f,FMT="((I4),(F9.4)))")n
    WRITE(f,FMT="('(',I4,'(F9.4))')")cols
    !matrix(:) = (/ (real(r), r=1,m*n) /)
    DO r = 1, rows
       !WRITE(*,FMT=f)matrix(r,1:(n-1)*m+r:n)
       IF(.NOT.t) THEN
          WRITE(*,FMT=f)matrix(r:r+(cols-1)*rows:rows)
       ELSE
          WRITE(*,FMT=f)matrix((r-1)*cols+1:r*cols)
       ENDIF
    ENDDO
  END SUBROUTINE printmat_s

  SUBROUTINE printmat_d(matrix, rows, cols, title, tr)
    REAL(KIND=real_8), DIMENSION(:), &
      INTENT(IN)                             :: matrix
    INTEGER, INTENT(IN)                      :: rows, cols
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

    IF (PRESENT (title)) THEN
       IF (PRESENT (tr)) THEN
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, title, tr)
       ELSE
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, title)
       ENDIF
    ELSE
       IF (PRESENT (tr)) THEN
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, tr=tr)
       ELSE
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols)
       ENDIF
    ENDIF
  END SUBROUTINE printmat_d

  SUBROUTINE printmat_c(matrix, rows, cols, title, tr)
    COMPLEX(KIND=real_4), DIMENSION(:), &
      INTENT(IN)                             :: matrix
    INTEGER, INTENT(IN)                      :: rows, cols
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

    IF (PRESENT (title)) THEN
       IF (PRESENT (tr)) THEN
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, title, tr)
       ELSE
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, title)
       ENDIF
    ELSE
       IF (PRESENT (tr)) THEN
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, tr=tr)
       ELSE
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols)
       ENDIF
    ENDIF
  END SUBROUTINE printmat_c

  SUBROUTINE printmat_z(matrix, rows, cols, title, tr)
    COMPLEX(KIND=real_8), DIMENSION(:), &
      INTENT(IN)                             :: matrix
    INTEGER, INTENT(IN)                      :: rows, cols
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

    IF (PRESENT (title)) THEN
       IF (PRESENT (tr)) THEN
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, title, tr)
       ELSE
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, title)
       ENDIF
    ELSE
       IF (PRESENT (tr)) THEN
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols, tr=tr)
       ELSE
          CALL printmat_s(REAL(matrix, KIND=sp), rows, cols)
       ENDIF
    ENDIF
  END SUBROUTINE printmat_z


! *****************************************************************************
!> \brief Prints the elements of a matrix.
!> rows and columns is the matrix size. tr specifies whether
!> the elements are stored transposed
! *****************************************************************************
  SUBROUTINE dbcsr_print2dmat(matrix, title, tr)
    REAL(KIND=real_8), DIMENSION(:, :), &
      INTENT(IN)                             :: matrix
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    LOGICAL, INTENT(IN), OPTIONAL            :: tr

    CHARACTER(30)                            :: f
    INTEGER                                  :: cols, m, n, r, rows
    LOGICAL                                  :: dbg, t

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg, on=.FALSE.)
    !IF (.NOT.dbg.OR..NOT.bcsr_verbose) RETURN
    IF (.NOT.dbg) RETURN
    rows = SIZE(matrix,1)
    cols = SIZE(matrix,2)
    m = rows
    n = cols
    t = .FALSE.
    IF(PRESENT(title)) WRITE(*,*)title
    IF(PRESENT(tr)) THEN
       IF(tr) THEN
          t = .TRUE.
          m = cols
          n = rows
       ENDIF
    ENDIF
    WRITE(*,*)'Print Matrix 2D:',rows,cols,t,'range',&
         LBOUND(matrix),':',UBOUND(matrix)
    IF(m.GT.10000) m=0
    IF(n.GT.10000) n=0
    IF(m*n.LT.1 .OR. m*n.GT.SIZE(matrix)) RETURN
    WRITE(f,FMT="('(',I4,'(F9.4))')")cols
    DO r = 1, rows
       !WRITE(*,FMT=f)matrix(r,1:(n-1)*m+r:n)
       IF(.NOT.t) THEN
          WRITE(*,FMT=f)matrix(r,:)
       ELSE
          WRITE(*,FMT=f)matrix(:,r)
       ENDIF
    ENDDO
  END SUBROUTINE dbcsr_print2dmat


! *****************************************************************************
!> \brief Swaps two integers
!> \param[in,out] a,b         Integers to swap
! *****************************************************************************
  SUBROUTINE iswap(a, b)
    INTEGER, INTENT(INOUT)                   :: a, b

    INTEGER                                  :: tmp

!   ---------------------------------------------------------------------------

    tmp = a
    a = b
    b = tmp
  END SUBROUTINE iswap


! *****************************************************************************
!> \brief In-place block transpose.
!> \param[in,out] extent      Matrix in the form of a 1-d array
!> \param[in] rows, columns   Input matrix size
!> \param error     cp2k error
! *****************************************************************************
  PURE SUBROUTINE block_transpose_inplace(extent, rows, columns)
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: extent
    INTEGER, INTENT(IN)                      :: rows, columns

    CHARACTER(len=*), PARAMETER :: routineN = 'block_transpose_inplace', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    extent = RESHAPE(TRANSPOSE(&
         RESHAPE(extent, (/rows, columns/))), (/rows*columns/))
  END SUBROUTINE block_transpose_inplace

! *****************************************************************************
!> \brief Copy data from one data area to another.
!>
!> There are no checks done for correctness!
!> \param[in] dst        destination data area
!> \param[in] src        source data area
!> \param[in] lb         lower bound for destination (and source if
!>                       not given explicity)
!> \param[in] data_size  number of elements to copy
!> \param[in] source_lb  (optional) lower bound of source
! *****************************************************************************
  SUBROUTINE dbcsr_block_transpose_aa (dst, src, lb,& 
       row_size, col_size, source_lb)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dst
    TYPE(dbcsr_data_obj), INTENT(IN)         :: src
    INTEGER, INTENT(IN)                      :: lb, row_size, col_size
    INTEGER, INTENT(IN), OPTIONAL            :: source_lb

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_block_transpose_aa', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_size, lb_s, ub, ub_s

!   ---------------------------------------------------------------------------

    IF (bcsr_debug) THEN
       CALL cp_assert (ASSOCIATED(dst%d) .AND. ASSOCIATED(src%d),&
            cp_fatal_level, cp_caller_error, routineN,&
            "Data areas must be setup.")
       CALL cp_assert (dst%d%data_type .EQ. src%d%data_type, cp_fatal_level,&
            cp_caller_error, routineN, "Data type must be the same.")
       CALL cp_assert (dst%d%data_type.EQ.dbcsr_type_real_8.OR.dst%d%data_type.EQ.dbcsr_type_real_4&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_8.OR.dst%d%data_type.EQ.dbcsr_type_complex_4,&
            cp_warning_level, cp_caller_error, routineN, "Incorrect data type.")
    ENDIF
    data_size = row_size*col_size
    ub = lb + data_size - 1
    IF (PRESENT (source_lb)) THEN
       lb_s = source_lb
       ub_s = source_lb + data_size-1
    ELSE
       lb_s = lb
       ub_s = ub
    ENDIF
    SELECT CASE (src%d%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_block_transpose (dst%d%r_dp(lb:ub), src%d%r_dp(lb_s:ub_s),&
            row_size, col_size)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_block_transpose (dst%d%r_sp(lb:ub), src%d%r_sp(lb_s:ub_s),&
            row_size, col_size)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_block_transpose (dst%d%c_dp(lb:ub), src%d%c_dp(lb_s:ub_s),&
            row_size, col_size)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_block_transpose (dst%d%c_sp(lb:ub), src%d%c_sp(lb_s:ub_s),&
            row_size, col_size)
    END SELECT
  END SUBROUTINE dbcsr_block_transpose_aa


! *****************************************************************************
!> \brief Copy data from one data area to another.
!>
!> There are no checks done for correctness!
!> \param[in] dst        destination data area
!> \param[in] lb         lower bound for destination (and source if
!>                       not given explicity)
!> \param[in] data_size  number of elements to copy
!> \param[in] src        source data area
!> \param[in] source_lb  (optional) lower bound of source
! *****************************************************************************
  SUBROUTINE dbcsr_data_copy_aa (dst, lb, data_size, src, source_lb)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dst
    INTEGER, INTENT(IN)                      :: lb, data_size
    TYPE(dbcsr_data_obj), INTENT(IN)         :: src
    INTEGER, INTENT(IN), OPTIONAL            :: source_lb

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_copy_aa', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lb_s, ub, ub_s

!   ---------------------------------------------------------------------------

    IF (bcsr_debug) THEN
       CALL cp_assert (ASSOCIATED(dst%d) .AND. ASSOCIATED(src%d),&
            cp_fatal_level, cp_caller_error, routineN,&
            "Data areas must be setup.")
       CALL cp_assert (dst%d%data_type .EQ. src%d%data_type, cp_fatal_level,&
            cp_caller_error, routineN, "Data type must be the same.")
       CALL cp_assert (dst%d%data_type.EQ.dbcsr_type_real_8.OR.dst%d%data_type.EQ.dbcsr_type_real_4&
            .OR.dst%d%data_type.EQ.dbcsr_type_complex_8.OR.dst%d%data_type.EQ.dbcsr_type_complex_4,&
            cp_warning_level, cp_caller_error, routineN, "Incorrect data type.")
    ENDIF
    ub = lb + data_size - 1
    IF (PRESENT (source_lb)) THEN
       lb_s = source_lb
       ub_s = source_lb + data_size-1
    ELSE
       lb_s = lb
       ub_s = ub
    ENDIF
    SELECT CASE (src%d%data_type)
    CASE (dbcsr_type_real_4)
       dst%d%r_sp(lb:ub) = src%d%r_sp(lb_s:ub_s)
    CASE (dbcsr_type_real_8)
       CALL cp_assert(ASSOCIATED(dst%d%r_dp),cp_fatal_level,&
            cp_caller_error, routineN, "associated(dst%d%r_dp)")
       CALL cp_assert(ASSOCIATED(src%d%r_dp),cp_fatal_level,&
            cp_caller_error, routineN, "associated(src%d%r_dp)")
       CALL cp_assert(lb.GE.LBOUND(dst%d%r_dp,1),cp_fatal_level,&
            cp_caller_error, routineN, "lb dst%d%r_dp")
       CALL cp_assert(ub.LE.UBOUND(dst%d%r_dp,1),cp_fatal_level,&
            cp_caller_error, routineN, "ub dst%d%r_dp")
       CALL cp_assert(lb_s.GE.LBOUND(src%d%r_dp,1),cp_fatal_level,&
            cp_caller_error, routineN, "lb src%d%r_dp")
       CALL cp_assert(ub_s.LE.UBOUND(src%d%r_dp,1),cp_fatal_level,&
            cp_caller_error, routineN, "ub src%d%r_dp")
       dst%d%r_dp(lb:ub) = src%d%r_dp(lb_s:ub_s)
    CASE (dbcsr_type_complex_4)
       dst%d%c_sp(lb:ub) = src%d%c_sp(lb_s:ub_s)
    CASE (dbcsr_type_complex_8)
       dst%d%c_dp(lb:ub) = src%d%c_dp(lb_s:ub_s)
    END SELECT
  END SUBROUTINE dbcsr_data_copy_aa


! *****************************************************************************
!> \brief Copies metadata into an array.
!> \param[in] matrix      Matrix
!> \param[out] meta       Metadata elements
! *****************************************************************************
  SUBROUTINE dbcsr_pack_meta(matrix, meta)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, DIMENSION(dbcsr_meta_size), &
      INTENT(OUT)                            :: meta

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_pack_meta', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    meta(1)  = matrix%nblks
    meta(2)  = matrix%nze
    meta(3)  = matrix%nblkrows_total
    meta(4)  = matrix%nblkcols_total
    meta(5)  = matrix%nfullrows_total
    meta(6)  = matrix%nfullcols_total
    meta(7)  = matrix%nblkrows_local
    meta(8)  = matrix%nblkcols_local
    meta(9)  = matrix%nfullrows_local
    meta(10) = matrix%nfullcols_local
    SELECT CASE (matrix%matrix_type)
    CASE ('N')
       meta(11) = 0
    CASE ('T')
       meta(11) = 1
    CASE ('S')
       meta(11) = 2
    CASE ('A')
       meta(11) = 3
    CASE DEFAULT
       meta(11) = -1
       CALL bcsr_assert(5, .FALSE., routineP, "Invalid matrix type.")
    END SELECT
  END SUBROUTINE dbcsr_pack_meta


! *****************************************************************************
!> \brief Copies row and column block sizes from another matrix.
!> \param[in,out] matrix      target matrix
!> \param[in] meta            source matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_copy_block_sizes(dst, src)
    TYPE(dbcsr_type), INTENT(INOUT)          :: dst
    TYPE(dbcsr_type), INTENT(IN)             :: src

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_block_sizes', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    dst%row_blk_size = src%row_blk_size
    CALL array_hold (dst%row_blk_size)
    dst%col_blk_size = src%col_blk_size
    CALL array_hold (dst%col_blk_size)
  END SUBROUTINE dbcsr_copy_block_sizes


! *****************************************************************************
!> \brief Sets debugging for an individual subroutine depending on the
!>        value of the bcsr_debug module variable, default initialization in
!>        the declared subroutine, or 
!> \param[in,out] dbg         debug flag
!> \param[in] on    (optional) turn debugging on or off in this subroutine
!> \param[in,out] info        print informational messages (less than debug)
! *****************************************************************************
  PURE SUBROUTINE dbcsr_set_debug(dbg, on, info)
    LOGICAL, INTENT(INOUT)                   :: dbg
    LOGICAL, INTENT(IN), OPTIONAL            :: on
    LOGICAL, INTENT(INOUT), OPTIONAL         :: info

!   ---------------------------------------------------------------------------

    IF (PRESENT (info)) info = info .AND. bcsr_info
    IF (PRESENT (on)) THEN
       dbg = on
       RETURN
    ENDIF
    IF (.NOT. bcsr_debug) THEN
       dbg = .FALSE.
       RETURN
    ENDIF
  END SUBROUTINE dbcsr_set_debug


! *****************************************************************************
!> \brief Busy-loop sleep.
!> \param[in] duration        time to sleep/s.
! *****************************************************************************
  SUBROUTINE dbcsr_sleep(duration)
    REAL(KIND=dp), INTENT(IN)                :: duration

    REAL(KIND=dp)                            :: t_now, t_start

!   ---------------------------------------------------------------------------

    t_start = m_walltime ()
    t_now = t_start
    DO WHILE (t_now .LT. t_start + duration)
       t_now = m_walltime ()
    ENDDO
  END SUBROUTINE dbcsr_sleep


! *****************************************************************************
!> \brief Replaces an existing pointer with a new one, freeing memory as
!>        required.
!> \param[in,out] original_p  original pointer, to be replaced
!> \param[in] new_p           replacement pointer
! *****************************************************************************
  SUBROUTINE pointer_replace_i (original_p, new_p)
    INTEGER, DIMENSION(:), POINTER           :: original_p, new_p

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (original_p)) DEALLOCATE (original_p)
    original_p => new_p
  END SUBROUTINE pointer_replace_i

! *****************************************************************************
!> \brief Creates a DBCSR distribution from a distribution_2d
!> \param[in] dist2d          distribution_2d
!> \param[out] dist           DBCSR distribution
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_dist2d_to_dist(dist2d, dist, error)
    TYPE(distribution_2d_type), INTENT(IN), &
      TARGET                                 :: dist2d
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: mypcol, myproc, myprow, &
                                                numproc
    INTEGER, DIMENSION(:), POINTER           :: col_dist_data, row_dist_data
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    TYPE(array_i1d_obj)                      :: cd, rd
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_mp_obj)                       :: mp_env
    TYPE(distribution_2d_type), POINTER      :: dist2d_p

!
!   ---------------------------------------------------------------------------

    dist2d_p => dist2d
    CALL distribution_2d_get(dist2d_p, error=error,&
         row_distribution=row_dist_data, col_distribution=col_dist_data,&
         blacs_env=blacs_env)
    CALL get_blacs_info(blacs_env, para_env=para_env,&
         my_process_row=myprow, my_process_column=mypcol,&
         blacs2mpi=pgrid)
    myproc = para_env%mepos
    numproc = para_env%num_pe
    CALL dbcsr_mp_new(mp_env, pgrid, para_env%group, myproc, numproc,&
         myprow, mypcol)
    CALL array_nullify (rd)
    CALL array_nullify (cd)
    CALL array_new(rd, row_dist_data)
    CALL array_new(cd, col_dist_data)
    CALL dbcsr_distribution_new(dist, mp_env, rd, cd)
    CALL dbcsr_mp_release (mp_env)
    CALL array_release (rd)
    CALL array_release (cd)
  END SUBROUTINE dbcsr_dist2d_to_dist


! *****************************************************************************
!> \brief Prints a distribution
!> \param[in] dist  distribution to print
! *****************************************************************************
  SUBROUTINE dbcsr_print_dist(dist)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist

    TYPE(dbcsr_mp_obj)                       :: mp

    mp = dbcsr_distribution_mp (dist)
    WRITE(*,'(1X,A,I5,"x",I5)')"Distribution process grid size",&
         dbcsr_mp_nprows (mp), dbcsr_mp_npcols (mp)
    IF (array_size(dist%d%row_dist).gt.0) &
         WRITE(*,*)"Distribution: Row dist", array_data(dist%d%row_dist)
    IF (array_size(dist%d%col_dist).gt.0) &
         WRITE(*,*)"Distribution: Col dist", array_data(dist%d%col_dist)
    IF (array_size(dist%d%local_rows).gt.0) &
         WRITE(*,*)"Distribution: Local rows", array_data(dist%d%local_rows)
    IF (array_size(dist%d%local_cols).gt.0) &
    WRITE(*,*)"Distribution: Local cols", array_data(dist%d%local_cols)
  END SUBROUTINE dbcsr_print_dist


! *****************************************************************************
!> \brief Transposes a distribution
!> \param[out] dist_tr        transposed distribution
!> \param[in] dist2d          current distribution
! *****************************************************************************
  SUBROUTINE dbcsr_transpose_distribution(dist_tr, dist_normal)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_tr
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_normal

    TYPE(array_i1d_obj)                      :: cd, lc, lr, rd
    TYPE(dbcsr_mp_obj)                       :: mp_env_tr

!   ---------------------------------------------------------------------------

    rd = dbcsr_distribution_col_dist (dist_normal)
    cd = dbcsr_distribution_row_dist (dist_normal)
    lr = dbcsr_distribution_local_cols (dist_normal)
    lc = dbcsr_distribution_local_rows (dist_normal)
    CALL dbcsr_mp_new_transposed (mp_env_tr, dbcsr_distribution_mp (&
         dist_normal))
    CALL dbcsr_distribution_new(dist_tr, mp_env_tr, rd, cd, lr, lc)
    CALL dbcsr_mp_release (mp_env_tr)
  END SUBROUTINE dbcsr_transpose_distribution

! *****************************************************************************
!> \brief Transposes a distribution but keeps the same mp_env
!> \param[out] dist_tr        transposed distribution
!> \param[in] dist2d          current distribution
! *****************************************************************************
  SUBROUTINE dbcsr_transpose_dims(dist_tr, dist_normal)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist_tr
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist_normal

    INTEGER, DIMENSION(:), POINTER           :: cda, rda
    TYPE(array_i1d_obj)                      :: cd, new_cd, new_rd, rd
    TYPE(dbcsr_mp_obj)                       :: mp

!   ---------------------------------------------------------------------------

    mp = dbcsr_distribution_mp (dist_normal)
    rd = dbcsr_distribution_row_dist (dist_normal)
    cd = dbcsr_distribution_col_dist (dist_normal)
    ALLOCATE (rda (array_size (cd)))
    CALL evenly_bin (rda, dbcsr_mp_nprows (mp))
    CALL array_new (new_rd, rda, gift=.TRUE.)

    ALLOCATE (cda (array_size (rd)))
    CALL evenly_bin (cda, dbcsr_mp_npcols (mp))
    CALL array_new (new_cd, cda, gift=.TRUE.)

    CALL dbcsr_distribution_new(dist_tr,&
         mp, new_rd, new_cd)
    CALL array_release (new_rd)
    CALL array_release (new_cd)
  END SUBROUTINE dbcsr_transpose_dims


! *****************************************************************************
!> \brief Evenly distributes elements into bins
!> \param[out] array          elements to bin
!> \param[in] nbins           number of bins into which elemnts are distributed
! *****************************************************************************
  PURE SUBROUTINE evenly_bin (array, nbins)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: array
    INTEGER, INTENT(IN)                      :: nbins

    INTEGER                                  :: bin, first, last, maxbinsize

!   ---------------------------------------------------------------------------

    maxbinsize = CEILING (REAL(SIZE (array)) / REAL(nbins))
    DO bin = 1, nbins
       first = 1 + (bin-1)*maxbinsize
       last = MIN (1 + (bin)*maxbinsize-1, SIZE (array))
       array (first : last) = bin - 1
    ENDDO
  END SUBROUTINE evenly_bin


! *****************************************************************************
!> \brief Makes new distribution with decimation and multiplicity
!> \par Definition of multiplicity and nimages
!>      Multiplicity and number of images are used to match process
!>      grid coordinates on non-square process grids. Given
!>      source_nbins and target_nbins, their relation is
!>          source_nbins * multiplicity = target_nbins * nimages.
!>      It is best when both multiplicity and nimages are small. To
!>      get these two factors, then, one can use the following formulas:
!>          nimages      = lcm(source_nbins, target_nbins) / target_nbins
!>          multiplicity = target_nbins / gcd(source_nbins, target_nbins)
!> \par Mapping
!>      The new distribution comprises of real bins and images within
!>      bins. These can be view as target_nbins*nimages virtual
!>      columns. These same virtual columns are also
!>      source_nbins*multiplicity in number. Therefore these virtual
!>      columns are mapped from source_nbins*multiplicity onto
!>      target_bins*nimages (each target bin has nimages images):
!>      Source 4: |1 2 3|4 5 6|7 8 9|A B C| (4*3)
!>      Target 6: |1 2|3 4|5 6|7 8|9 A|B C| (6*2)
!>      multiplicity=3, nimages=2, 12 virtual columns (1-C).
!>      Source bin elements are evenly mapped into one of multiplicity
!>      virtual columns. Other (non-even, block-size aware) mappings
!>      could be better.
!> \param[out] new_bins      new real distribution
!> \param[out] images        new image distribution
!> \param[in] source_bins    Basis for the new distribution and images
!> \param[in] nbins          number of bins in the new real distribution
!> \param[in] multiplicity   multiplicity
!> \param[in] nimages        number of images in the new distribution
! *****************************************************************************
  SUBROUTINE rebin_distribution (new_bins, images, source_bins,&
       nbins, multiplicity, nimages)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: new_bins, images
    INTEGER, DIMENSION(:), INTENT(IN)        :: source_bins
    INTEGER, INTENT(IN)                      :: nbins, multiplicity, nimages

    CHARACTER(len=*), PARAMETER :: routineN = 'rebin_distribution', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 2

    INTEGER                                  :: bin, i, old_nbins, virtual_bin
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bin_multiplier

!   ---------------------------------------------------------------------------

    CALL bcsr_assert (debug_level, SIZE(source_bins).EQ.SIZE(new_bins),&
         routineN, "array sizes not equal")
    CALL bcsr_assert (debug_level, SIZE(source_bins).EQ.SIZE(images),&
         routineN, "array sizes not equal")
    CALL bcsr_assert (debug_level, MOD (nbins, multiplicity) .EQ. 0,&
         routineN, "mulitplicity is not divisor of new process grid coordinate")
    old_nbins = nbins * nimages / multiplicity
    ALLOCATE (bin_multiplier(0:old_nbins-1))
    bin_multiplier(:) = 0
    DO i = 1, SIZE (source_bins)
       bin = source_bins (i)
       bin_multiplier(bin) = bin_multiplier(bin)+1
       IF (bin_multiplier(bin) .GE. multiplicity) THEN
          bin_multiplier(bin) = 0
       ENDIF
       virtual_bin = bin*multiplicity+bin_multiplier(bin)
       new_bins(i) = virtual_bin / nimages
       images(i) = 1 + MOD (virtual_bin, nimages)
    ENDDO
  END SUBROUTINE rebin_distribution

! *****************************************************************************
!> \brief Makes new distribution with decimation and multiplicity
!> \par Definition of multiplicity and nimages
!>      Multiplicity and number of images are used to match process
!>      grid coordinates on non-square process grids. Given
!>      source_nbins and target_nbins, their relation is
!>          source_nbins * multiplicity = target_nbins * nimages.
!>      It is best when both multiplicity and nimages are small. To
!>      get these two factors, then, one can use the following formulas:
!>          nimages      = lcm(source_nbins, target_nbins) / target_nbins
!>          multiplicity = target_nbins / gcd(source_nbins, target_nbins)
!> \par Mapping
!>      The new distribution comprises of real bins and images within
!>      bins. These can be view as target_nbins*nimages virtual
!>      columns. These same virtual columns are also
!>      source_nbins*multiplicity in number. Therefore these virtual
!>      columns are mapped from source_nbins*multiplicity onto
!>      target_bins*nimages (each target bin has nimages images):
!>      Source 4: |1 2 3|4 5 6|7 8 9|A B C| (4*3)
!>      Target 6: |1 2|3 4|5 6|7 8|9 A|B C| (6*2)
!>      multiplicity=3, nimages=2, 12 virtual columns (1-C).
!>      Source bin elements are evenly mapped into one of multiplicity
!>      virtual columns. Other (non-even, block-size aware) mappings
!>      could be better.
!> \param[out] images        new image distribution
!> \param[in] my_bins        Basis for the new images
!> \param[in] nbins          number of bins in the new real distribution
!> \param[in] multiplicity   multiplicity
!> \param[in] nimages        number of images in the new distribution
! *****************************************************************************
  SUBROUTINE reimage_distribution (images, my_bins,&
       nbins, multiplicity, nimages)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: images
    INTEGER, DIMENSION(:), INTENT(IN)        :: my_bins
    INTEGER, INTENT(IN)                      :: nbins, multiplicity, nimages

    CHARACTER(len=*), PARAMETER :: routineN = 'reimage_distribution', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 2

    INTEGER                                  :: bin, i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bin_multiplier

!   ---------------------------------------------------------------------------

    CALL bcsr_assert (debug_level, SIZE(my_bins).EQ.SIZE(images),&
         routineN, "array sizes not equal")
    ALLOCATE (bin_multiplier(0:nbins-1))
    bin_multiplier(:) = 0
    DO i = 1, SIZE (my_bins)
       bin = my_bins (i)
       bin_multiplier(bin) = bin_multiplier(bin)+1
       IF (bin_multiplier(bin) .GE. nimages) THEN
          bin_multiplier(bin) = 0
       ENDIF
       images(i) = 1 + bin_multiplier(bin)
    ENDDO
  END SUBROUTINE reimage_distribution


! *****************************************************************************
!> \brief Makes new distribution with multiplicity
!> \par Definition of multiplicity and nimages
!>      Multiplicity and number of images are used to match process
!>      grid coordinates on non-square process grids. Given
!>      source_nbins and target_nbins, their relation is
!>          source_nbins * multiplicity = target_nbins * nimages.
!>      It is best when both multiplicity and nimages are small. To
!>      get these two factors, then, one can use the following formulas:
!>          nimages      = lcm(source_nbins, target_nbins) / target_nbins
!>          multiplicity = target_nbins / gcd(source_nbins, target_nbins)
!> \par Mapping
!>      The new distribution comprises of real bins and images within
!>      bins. These can be view as target_nbins*nimages virtual
!>      columns. These same virtual columns are also
!>      source_nbins*multiplicity in number. Therefore these virtual
!>      columns are mapped from source_nbins*multiplicity onto
!>      target_bins*nimages (each target bin has nimages images):
!>      Source 4: |1 2 3|4 5 6|7 8 9|A B C| (4*3)
!>      Target 6: |1 2|3 4|5 6|7 8|9 A|B C| (6*2)
!>      multiplicity=3, nimages=2, 12 virtual columns (1-C).
!>      Source bin elements are evenly mapped into one of multiplicity
!>      virtual columns. Other (non-even, block-size aware) mappings
!>      could be better.
!> \param[out] new_bins      new real distribution
!> \param[out] images        new image distribution
!> \param[in] source_bins    Basis for the new distribution and images
!> \param[in] source_images  Basis for the new distribution and images
!> \param[in] nbins          number of bins in the new real distribution
!> \param[in] multiplicity   multiplicity
!> \param[in] nimages        number of images in the new distribution
! *****************************************************************************
  SUBROUTINE rebin_imaged_distribution (new_bins, images,&
       source_bins, source_images, nbins, multiplicity, nimages)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: new_bins, images
    INTEGER, DIMENSION(:), INTENT(IN)        :: source_bins, source_images
    INTEGER, INTENT(IN)                      :: nbins, multiplicity, nimages

    CHARACTER(len=*), PARAMETER :: routineN = 'rebin_imaged_distribution', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 2

    INTEGER                                  :: bin, i, virtual_bin

!   ---------------------------------------------------------------------------

    CALL bcsr_assert (debug_level, SIZE(source_bins).EQ.SIZE(new_bins),&
         routineN, "array sizes not equal")
    CALL bcsr_assert (debug_level, SIZE(source_bins).EQ.SIZE(images),&
         routineN, "array sizes not equal")
    CALL bcsr_assert (debug_level, MOD (nbins, multiplicity) .EQ. 0,&
         routineN, "mulitplicity is not divisor of new process grid coordinate")
    DO i = 1, SIZE (source_bins)
       virtual_bin = source_bins(i)*multiplicity + source_images(i)-1
       new_bins(i) = virtual_bin / nimages
       images(i) = 1 + MOD (virtual_bin, nimages)
       CALL cp_assert (new_bins(i).lt.nbins, cp_fatal_level, cp_internal_error,&
            routineN, "Wrong bin calculation")
       CALL cp_assert (images(i).le.nimages, cp_fatal_level, cp_internal_error,&
            routineN, "Wrong image calculation")
    ENDDO
  END SUBROUTINE rebin_imaged_distribution



! *****************************************************************************
!> \brief Creates an image distribution with column images.
!>
!> The column images make the column distribution compatible with the
!> row distribution of the supplied other_dist distribution.
!>
!> This routine is used to make a left matrix (with new_row_dist)
!> column-comptable with a given right matrix (other_dist) for
!> multiplication with non-transposed processor grids.
!>
!> \param[out] imgdist        image distribution with columns
!> \param[in] new_row_dist    use this row distribution
!> \param[in] other_dist      makes the distribution compatible with and based
!>                            on this one
! *****************************************************************************
  SUBROUTINE dbcsr_create_col_image_dist(imgdist, new_row_dist, other_dist)
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(OUT)                            :: imgdist
    TYPE(array_i1d_obj), INTENT(IN)          :: new_row_dist
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: other_dist

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_col_image_dist', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: debug_level, new_ncols, &
                                                npcols, nprows
    INTEGER, DIMENSION(:), POINTER           :: new_col_dist_data, &
                                                new_col_img_data, rd
    TYPE(array_i1d_obj)                      :: new_col_dist, new_col_image
    TYPE(dbcsr_distribution_obj)             :: new_dist
    TYPE(dbcsr_mp_obj)                       :: other_mp_env

    DEBUG_HEADER
!
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    IF (dbg) debug_level = 5
    CALL dbcsr_distribution_init (new_dist)
    CALL dbcsr_distribution_init (imgdist%main)
    CALL array_nullify (imgdist%col_image)
    CALL array_nullify (imgdist%row_image)
    other_mp_env = dbcsr_distribution_mp (other_dist)
    ! Determine the factors.
    new_ncols = dbcsr_distribution_nrows (other_dist) ! They must be compatible
    nprows = dbcsr_mp_nprows (other_mp_env)
    npcols = dbcsr_mp_npcols (other_mp_env)
    imgdist%col_decimation = lcm (nprows, npcols) / npcols
    imgdist%col_multiplicity = npcols / gcd (nprows, npcols)
    DBG 'image col decimation, multiplicity',&
         imgdist%col_decimation, imgdist%col_multiplicity
    imgdist%row_decimation = 1
    imgdist%row_multiplicity = 1
    ! Create the new column distribution and column image distribution
    ALLOCATE (new_col_dist_data (new_ncols))
    ALLOCATE (new_col_img_data (new_ncols))
    rd => array_data (dbcsr_distribution_row_dist (other_dist))
    CALL rebin_distribution (new_col_dist_data, new_col_img_data, rd,&
         npcols, imgdist%col_multiplicity, imgdist%col_decimation)
    CALL array_new (new_col_dist, new_col_dist_data, gift=.TRUE.)
    CALL array_new (new_col_image, new_col_img_data, gift=.TRUE.)
    ! Copy the row distribution from old distribution
    CALL dbcsr_distribution_new(new_dist,&
         other_mp_env,&
         new_row_dist, new_col_dist)
    CALL array_release (new_col_dist)
    ! Now finish the distribution image.
    imgdist%main = new_dist
    imgdist%col_image = new_col_image
  END SUBROUTINE dbcsr_create_col_image_dist


! *****************************************************************************
!> \brief Creates an image distribution given the other compatibility images
!> \param[out] imgdist                  distribution repetition
!> \param[in] dist                      distribution for which to form the
!>                                      image distribution
!> \param[in] match_row_pdist           (optional) match the new row
!>                                      distribution to this row distribution
!> \param[in] match_row_idist           (optional) match the row
!>                                      distribution to these row images
!> \param[in] match_row_nbins           number of bins in the distribution to
!>                                      match the local rows
!> \param[in] match_col_pdist           (optional) match the new column
!>                                      distribution to this column
!>                                      distribution
!> \param[in] match_col_idist           (optional) match the new column
!>                                      distribution to these column images
!> \param[in] match_col_nbins           number of bins in the distribution to
!>                                      match the local columns
! *****************************************************************************
  SUBROUTINE dbcsr_create_image_dist(imgdist, dist,&
       match_row_pdist, match_row_idist, match_row_nbins,&
       match_col_pdist, match_col_idist, match_col_nbins)
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(OUT)                            :: imgdist
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL, TARGET                       :: match_row_pdist, &
                                                match_row_idist
    INTEGER, INTENT(IN)                      :: match_row_nbins
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: match_col_pdist, &
                                                match_col_idist
    INTEGER, INTENT(IN)                      :: match_col_nbins

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_image_dist', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ncols, npcols, nprows, nrows
    INTEGER, DIMENSION(:), POINTER           :: col_dist_data, col_img_data, &
                                                row_dist_data, row_img_data
    LOGICAL                                  :: dbg
    TYPE(array_i1d_obj)                      :: col_dist, col_image, &
                                                row_dist, row_image
    TYPE(dbcsr_distribution_obj)             :: new_dist
    TYPE(dbcsr_mp_obj)                       :: mp_env

!   ---------------------------------------------------------------------------

    dbg = .FALSE.

    CALL dbcsr_distribution_init (new_dist)
    CALL dbcsr_distribution_init (imgdist%main)
    CALL array_nullify (imgdist%col_image)
    CALL array_nullify (imgdist%row_image)
    mp_env = dbcsr_distribution_mp (dist)
    ! Determine the factors.
    nrows = dbcsr_distribution_nrows (dist)
    ncols = dbcsr_distribution_ncols (dist)
    nprows = dbcsr_mp_nprows (mp_env)
    npcols = dbcsr_mp_npcols (mp_env)
    ! Create the new row distribution and row image distribution
    imgdist%row_decimation = lcm (nprows, match_row_nbins) / nprows
    imgdist%row_multiplicity = nprows / gcd (nprows, match_row_nbins)
    IF (dbg) WRITE (*,*)routineN//'row decimation, multiplicity',&
         imgdist%row_decimation, imgdist%row_multiplicity
    IF (dbg) WRITE(*,*)routineN//" nprows, match prows", nprows, match_row_nbins
    IF (match_row_nbins .EQ. nprows .AND. .NOT. PRESENT (match_row_pdist)) THEN
       row_dist = dbcsr_distribution_row_dist (dist)
       CALL array_hold (row_dist)
       ALLOCATE (row_img_data (nrows))
       row_img_data(:) = 1
       CALL array_new (row_image, row_img_data, gift=.TRUE.)
    ELSE
       ALLOCATE (row_img_data (nrows))
       IF (PRESENT (match_row_pdist)) THEN
          ALLOCATE (row_dist_data (nrows))
          IF (PRESENT (match_row_idist)) THEN
             CALL rebin_imaged_distribution (row_dist_data, row_img_data,&
                  match_row_pdist, match_row_idist,&
                  nprows,&
                  imgdist%row_multiplicity, imgdist%row_decimation)
          ELSE
             CALL rebin_distribution (row_dist_data, row_img_data,&
                  match_row_pdist, &
                  nprows,&
                  imgdist%row_multiplicity, imgdist%row_decimation)
          ENDIF
          CALL array_new (row_dist, row_dist_data, gift=.TRUE.)
       ELSE
          row_dist = dbcsr_distribution_row_dist (dist)
          CALL array_hold (row_dist)
          CALL reimage_distribution (row_img_data,&
               array_data (row_dist),&
               nprows,&
               imgdist%row_multiplicity, imgdist%row_decimation)
       ENDIF
       CALL array_new (row_image, row_img_data, gift=.TRUE.)
    ENDIF
    ! Create the new column distribution and column image distribution
    imgdist%col_decimation = lcm (npcols, match_col_nbins) / npcols
    imgdist%col_multiplicity = npcols / gcd (npcols, match_col_nbins)
    IF (dbg) WRITE (*,*)routineN//'col decimation, multiplicity',&
         imgdist%col_decimation, imgdist%col_multiplicity
    IF (dbg) WRITE(*,*)routineN//" npcols, match pcols", npcols, match_col_nbins
    IF (match_col_nbins .EQ. npcols .AND. .NOT. PRESENT (match_col_pdist)) THEN
       col_dist = dbcsr_distribution_col_dist (dist)
       CALL array_hold (col_dist)
       ALLOCATE (col_img_data (ncols))
       col_img_data(:) = 1
       CALL array_new (col_image, col_img_data, gift=.TRUE.)
    ELSE
       ALLOCATE (col_img_data (ncols))
       IF (PRESENT (match_col_pdist)) THEN
          ALLOCATE (col_dist_data (ncols))
          IF (PRESENT (match_col_idist)) THEN
             CALL rebin_imaged_distribution (col_dist_data, col_img_data,&
                  match_col_pdist, match_col_idist,&
                  npcols,&
                  imgdist%col_multiplicity, imgdist%col_decimation)
          ELSE
             CALL rebin_distribution (col_dist_data, col_img_data,&
                  match_col_pdist,&
                  npcols,&
                  imgdist%col_multiplicity, imgdist%col_decimation)
          ENDIF
          CALL array_new (col_dist, col_dist_data, gift=.TRUE.)
       ELSE
          col_dist = dbcsr_distribution_col_dist (dist)
          CALL array_hold (col_dist)
          CALL reimage_distribution (col_img_data,&
               array_data (col_dist),&
               npcols,&
               imgdist%col_multiplicity, imgdist%col_decimation)
       ENDIF
       CALL array_new (col_image, col_img_data, gift=.TRUE.)
    ENDIF
    ! Copy the row distribution from old distribution
    CALL dbcsr_distribution_new(new_dist,&
         mp_env,&
         row_dist, col_dist)
    CALL array_release (row_dist)
    CALL array_release (col_dist)
    ! Now finish the distribution image.
    imgdist%main = new_dist
    imgdist%col_image = col_image
    imgdist%row_image = row_image
  END SUBROUTINE dbcsr_create_image_dist

! *****************************************************************************
!> \brief Destroys a DBCSR distribution for a matrix multiplication based on
!>        the right matrix
!> \param[inout] imgdist_product        product distribution repetition
! *****************************************************************************
  SUBROUTINE dbcsr_destroy_image_dist(imgdist)
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(INOUT)                          :: imgdist

!
!   ---------------------------------------------------------------------------

    CALL array_release (imgdist%row_image)
    CALL array_release (imgdist%col_image)
    CALL dbcsr_distribution_release (imgdist%main)
  END SUBROUTINE dbcsr_destroy_image_dist

! *****************************************************************************
!> \brief Splits an array of int8 values into two int4 arrays.
!> \param[in] merged          array of merged values
!> \param[out] array_upper    array to fill with the upper bytes of the merged
!>                            values
!> \param[out] array_lower    array to fill with the lower bytes of the merged
!>                            values
! *****************************************************************************
  SUBROUTINE dbcsr_unpack_i8_2i4 (merged, array_upper, array_lower)
    INTEGER(KIND=int_8), DIMENSION(:), &
      INTENT(IN)                             :: merged
    INTEGER(KIND=int_4), DIMENSION(:), &
      INTENT(OUT)                            :: array_upper, array_lower

    INTEGER(KIND=int_8), PARAMETER           :: lmask8 = 4294967295_int_8

    INTEGER                                  :: i

!
!   ---------------------------------------------------------------------------
! Lmask is used to filter in the lower 4 bytes and so its lower 32 bits are
! set to 1: lmask8 = 2^32-1.
! Umask is used to filter in the higher 4 bytes and so its higher 32 bits
! are set to 1: umask8 = 2^32-1 << 32
!lmask8 = 4294967295 ! 2^32-1
!umask8 = 18446744069414584320 ! (2^32-1) * 2^32 = (2^64-1)-(2^32-1)

!$OMP DO SCHEDULE (STATIC)
    DO i = 1, SIZE (merged)
       array_upper(i) = ISHFT (merged(i), -32)
       array_lower(i) = IAND (merged(i), lmask8)
    ENDDO
!$OMP END DO
  END SUBROUTINE dbcsr_unpack_i8_2i4


! *****************************************************************************
!> \brief Releases all matrices in a 1-d arrray.
!> \param[in] source          input matrix
!> \param[in,out] marray      matrix array
! *****************************************************************************
  SUBROUTINE dbcsr_destroy_1d_array(marray)
    TYPE(dbcsr_array_type), INTENT(INOUT)    :: marray

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_destroy_1d_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    DO i = LBOUND(marray%mats,1), UBOUND (marray%mats,1)
       CALL dbcsr_destroy (marray%mats(i), force=.TRUE.)
    ENDDO
    DEALLOCATE (marray%mats)
    !CALL dbcsr_destroy_image_dist(marray%image_dist)
  END SUBROUTINE dbcsr_destroy_1d_array


! *****************************************************************************
!> \brief Releases all matrices in 2-d arrray.
!> \param[in] source          input matrix
!> \param[in,out] marray      matrix array
! *****************************************************************************
  SUBROUTINE dbcsr_destroy_2d_array(marray)
    TYPE(dbcsr_2d_array_type), INTENT(INOUT) :: marray

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_destroy_2d_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, row

!   ---------------------------------------------------------------------------

    DO row = LBOUND(marray%mats,1), UBOUND (marray%mats,1)
       DO col = LBOUND(marray%mats,2), UBOUND (marray%mats,2)
          CALL dbcsr_destroy (marray%mats(row, col), force=.TRUE.)
       ENDDO
    ENDDO
    DEALLOCATE (marray%mats)
  END SUBROUTINE dbcsr_destroy_2d_array


!! *****************************************************************************
!!> \brief Extracts a column from a matrix into a new matrix, local data only.
!!> \param[inout] vector_m     matrix into which the column is written
!!> \param[in] matrix          matrix in which the column resides
!!> \param[in] column          column to extract
!!> \param[in,out] error       cp2k error
!! *****************************************************************************
!  SUBROUTINE dbcsr_extract_column (vector_m, matrix, column, error)
!    TYPE(dbcsr_type), INTENT(INOUT)          :: vector_m
!    TYPE(dbcsr_type), INTENT(IN)             :: matrix
!    INTEGER, INTENT(IN)                      :: column
!    TYPE(cp_error_type), INTENT(INOUT)       :: error
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_extract_column', &
!      routineP = moduleN//':'//routineN
!
!    INTEGER                                  :: block_nze, col, nblks, &
!                                                next_data_p, row, src_blk, &
!                                                src_blk_p, src_nblk_rows_total
!    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
!    LOGICAL                                  :: symmetric_src, transposed_src
!
!!   ---------------------------------------------------------------------------
!
!    symmetric_src = (matrix%matrix_type.EQ.'S'.OR.matrix%matrix_type.EQ.'A')
!    transposed_src = (matrix%matrix_type.EQ.'T')
!    row_blk_size => array_data (matrix%row_blk_size)
!    col_blk_size => array_data (matrix%col_blk_size)
!    !> We first ensure there is enough space in the vector matrix,
!    !> then we build it up.
!    CALL dbcsr_get_info (matrix,&
!         nblkrows_total = src_nblk_rows_total)
!    !> Now go through the rows and find the appropriate column. If it
!    !> exists, place the data into the vector matrix and update the
!    !> index.
!    vector_m%nze = 0
!    nblks = 0
!    next_data_p = 1
!    CALL cp_assert(dbcsr_valid_index (vector_m), cp_warning_level,&
!         cp_wrong_args_error, routineP, "Vector matrix must be valid.",&
!         error=error)
!    vector_m%row_p(1) = 0
!    DO row = 1, src_nblk_rows_total
!       vector_m%row_p(row+1) = nblks
!       DO src_blk = matrix%row_p(row)+1, matrix%row_p(row+1)
!          col = matrix%col_i(src_blk)
!          !WRITE(*,*)col,'vs.',column,transposed_src,symmetric_src
!          IF ((.NOT. transposed_src .AND. col .EQ. column) .OR.&
!               (symmetric_src .AND. row .EQ. column) .OR.&
!               (transposed_src .AND. row .EQ. column)) THEN
!             nblks = nblks + 1
!             vector_m%row_p(row+1) = nblks
!             vector_m%blk_p(nblks) = next_data_p
!             vector_m%col_i(nblks) = 1
!             block_nze = row_blk_size(row) * col_blk_size(col)
!             src_blk_p = matrix%blk_p(src_blk)
!             CALL dbcsr_data_set(vector_m%data_area,&
!                  next_data_p, block_nze, matrix%data_area, src_blk_p)
!             !vector_m%data(next_data_p : next_data_p+block_nze-1) =&
!             !     matrix%data(src_blk_p : src_blk_p+block_nze-1)
!             next_data_p = next_data_p + block_nze
!             vector_m%nze = vector_m%nze + block_nze
!             WRITE(*,*)'added block',nblks,'in row',row
!          ENDIF ! coordinates match
!       ENDDO ! src_blk
!    ENDDO ! blk_row
!    !WRITE(*,*)nblks,'blocks',lbound(vector_m%row_p),ubound(vector_m%row_p)
!    !WRITE(*,*)'rowp=',vector_m%row_p
!    !WRITE(*,*)'col_i=',vector_m%col_i(1:nblks)
!    !WRITE(*,*)'blk_p=',vector_m%blk_p(1:nblks)
!    !vector_m%row_p(src_nblk_rows_total+1) = nblks
!    vector_m%nblks = nblks
!    vector_m%index(dbcsr_slot_nblks) = nblks
!    vector_m%index(dbcsr_slot_nze) = vector_m%nze
!    ! Are there any other elements we have to update, especially
!    ! regarding the number of blocks and the data size?
!  END SUBROUTINE dbcsr_extract_column


! *****************************************************************************
!> \brief Returns whether this work matrix uses the mutable type
!> \param[in] wm              work matrix
!> \param[out] use_mutable    use the mutable and not append-only working
!>                            structures
! *****************************************************************************
  PURE FUNCTION dbcsr_wm_use_mutable (wm) RESULT (use_mutable)
    TYPE(dbcsr_work_type), INTENT(IN)        :: wm
    LOGICAL                                  :: use_mutable

!   ---------------------------------------------------------------------------

    use_mutable = dbcsr_mutable_instantiated (wm%mutable)
  END FUNCTION dbcsr_wm_use_mutable

! *****************************************************************************
!> \brief Returns whether work matrices should use the mutable data type
!> \param[in] matrix          matrix
!> \param[out] use_mutable    use the mutable and not append-only working
!>                            structures
! *****************************************************************************
  PURE FUNCTION dbcsr_use_mutable (matrix) RESULT (use_mutable)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    LOGICAL                                  :: use_mutable

!   ---------------------------------------------------------------------------

    use_mutable = matrix%work_mutable
  END FUNCTION dbcsr_use_mutable


! *****************************************************************************
!> \brief Calculates the intersections of blocks
!> \par nints output format
!>      The ints array should be up to twice as large as the number of
!>      intersecting blocks. Each entry is comprised of the target
!>      block and the common length.
!> \par n_src_dsts format
!>      This arrays stored the number of intersecting blocks in common
!>      (position 2) and the offset of the first common intersecting
!>      block (position 1).
!> \note This routine is used in the counting and sending loops in
!>       dbcsr_complete_redistribute
!> \param[out] ints           intersections of blocks
!> \param[in] num_ints        maximum number of expected intersections
!> \param[out] n_src_dsts     offset and number intersections belonging
!>                            to source blocks
!> \param[in] src_sizes       sizes of source blocks
!> \param[in] dst_sizes       sizes of target blocks
! *****************************************************************************
  SUBROUTINE dbcsr_reblocking_targets (ints, numints, n_src_dsts,&
       src_sizes, dst_sizes)
    INTEGER, INTENT(INOUT)                   :: numints
    INTEGER, DIMENSION(2, numints), &
      INTENT(OUT)                            :: ints
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: n_src_dsts
    INTEGER, DIMENSION(:), INTENT(IN)        :: src_sizes, dst_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reblocking_targets', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: common_extent, current_dst, &
                                                current_int, current_src, &
                                                n_dst, n_src, s_dst, s_src

!   ---------------------------------------------------------------------------

    n_src = SIZE(src_sizes)
    n_dst = SIZE(dst_sizes)
    current_int = 1
    current_src = 1
    current_dst = 1
    s_src = src_sizes (current_src)
    s_dst = dst_sizes (current_dst)
    n_src_dsts(:,1) = (/ 1, 0 /)
    DO WHILE (current_src .LE. n_src .AND. current_dst .LE. n_dst)
       ! Calculate many elements do the current blocks have in common
       ! and record these as going to the current target block.
       common_extent = MIN (s_src, s_dst)
       ints(1, current_int) = current_dst ! target block
       ints(2, current_int) = common_extent
       ! We've used up the common extents.
       s_src = s_src - common_extent
       s_dst = s_dst - common_extent
       n_src_dsts(2,current_src) = n_src_dsts(2,current_src) + 1
       IF (s_src .EQ. 0) THEN
          current_src = current_src + 1
          IF (current_src .LE. n_src) THEN
             s_src = src_sizes (current_src)
             n_src_dsts(:, current_src) = (/ current_int+1, 0 /)
          ENDIF
       ENDIF
       IF (s_dst .EQ. 0) THEN
          current_dst = current_dst + 1
          IF (current_dst .LE. n_dst) s_dst = dst_sizes (current_dst)
       ENDIF
       current_int = current_int + 1
    ENDDO
    numints = current_int - 1
  END SUBROUTINE dbcsr_reblocking_targets


! *****************************************************************************
!> \brief Calls the functions that are called in the included files but not
!>        in thin one. Needed to make prettify think the functions ary really
!>        used.
!> \param[in] matrix          data_area
!> \param[out] temp           pointer to the data
! *****************************************************************************
  SUBROUTINE prettify_fooler (data_area, temp)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: data_area
    REAL(KIND=dp), DIMENSION(:), POINTER     :: temp

    TYPE(dbcsr_obj)                          :: matrix

    IF (.FALSE.) THEN
       temp => dbcsr_get_data_p (data_area, 0.0_dp)
       CALL dbcsr_get_info (matrix)
    ENDIF
  END SUBROUTINE prettify_fooler

!  INCLUDE "dbcsr_util_d_.F"
!  INCLUDE "dbcsr_util_z_.F"
#include "dbcsr_util_d_.F"
#include "dbcsr_util_z_.F"
!#if !defined (__SGL)
!  INCLUDE "dbcsr_util_s_.F"
!  INCLUDE "dbcsr_util_c_.F"
#include "dbcsr_util_s_.F"
#include "dbcsr_util_c_.F"
!#endif


END MODULE dbcsr_util

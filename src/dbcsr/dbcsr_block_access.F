!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR block access
!> \author  Urban Borstnik
!> \date    2010-02-18
!> \version 0.9
!>
!> <b>Modification history:</b>
!  - 2010-02-18 Moved from dbcsr_operations and dbcsr_methods
!  - 2010-04-22 Added block buffer operations
! *****************************************************************************
MODULE dbcsr_block_access
  USE array_types,                     ONLY: array_data
  USE btree_I8_k_cp2d_v,               ONLY: btree_2d_data_c => cp2d,&
                                             btree_add_c => btree_add,&
                                             btree_get_c => btree_find
  USE btree_I8_k_dp2d_v,               ONLY: btree_2d_data_d => dp2d,&
                                             btree_add_d => btree_add,&
                                             btree_get_d => btree_find
  USE btree_I8_k_sp2d_v,               ONLY: btree_2d_data_s => sp2d,&
                                             btree_add_s => btree_add,&
                                             btree_get_s => btree_find
  USE btree_I8_k_zp2d_v,               ONLY: btree_2d_data_z => zp2d,&
                                             btree_add_z => btree_add,&
                                             btree_get_z => btree_find
  USE dbcsr_block_operations,          ONLY: dbcsr_block_copy,&
                                             dbcsr_block_transpose,&
                                             dbcsr_data_clear,&
                                             dbcsr_data_set
  USE dbcsr_config,                    ONLY: default_resize_factor
  USE dbcsr_data_methods,              ONLY: dbcsr_data_clear_pointer,&
                                             dbcsr_data_ensure_size,&
                                             dbcsr_data_set_pointer,&
                                             dbcsr_get_data,&
                                             dbcsr_get_data_p
  USE dbcsr_dist_operations,           ONLY: dbcsr_get_block_index,&
                                             dbcsr_get_stored_block_info,&
                                             dbcsr_get_stored_coordinates
  USE dbcsr_error_handling,            ONLY: &
       dbcsr_assert, dbcsr_caller_error, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_failure_level, dbcsr_fatal_level, &
       dbcsr_internal_error, dbcsr_warning_level, dbcsr_wrong_args_error
  USE dbcsr_index_operations,          ONLY: dbcsr_addto_index_array,&
                                             dbcsr_clearfrom_index_array,&
                                             dbcsr_expand_row_index,&
                                             dbcsr_make_dbcsr_index,&
                                             dbcsr_sort_indices,&
                                             merge_index_arrays
  USE dbcsr_methods,                   ONLY: &
       dbcsr_blk_column_size, dbcsr_blk_row_size, dbcsr_distribution, &
       dbcsr_distribution_mp, dbcsr_get_data_size_referenced, &
       dbcsr_get_data_type, dbcsr_get_num_blocks, dbcsr_mp_mynode, &
       dbcsr_mutable_instantiated, dbcsr_mutable_new, dbcsr_nblkcols_total, &
       dbcsr_nblkrows_total, dbcsr_use_mutable, dbcsr_wm_use_mutable
  USE dbcsr_ptr_util,                  ONLY: pointer_c_rank_remap2,&
                                             pointer_d_rank_remap2,&
                                             pointer_s_rank_remap2,&
                                             pointer_view,&
                                             pointer_z_rank_remap2
  USE dbcsr_toollib,                   ONLY: make_coordinate_tuple,&
                                             swap
  USE dbcsr_types,                     ONLY: &
       dbcsr_data_obj, dbcsr_obj, dbcsr_scalar_type, dbcsr_slot_blk_p, &
       dbcsr_slot_col_i, dbcsr_slot_nblks, dbcsr_slot_nze, &
       dbcsr_type_complex_4, dbcsr_type_complex_4_2d, dbcsr_type_complex_8, &
       dbcsr_type_complex_8_2d, dbcsr_type_real_4, dbcsr_type_real_4_2d, &
       dbcsr_type_real_8, dbcsr_type_real_8_2d
  USE dbcsr_work_operations,           ONLY: add_work_coordinate,&
                                             dbcsr_work_create
  USE kinds,                           ONLY: dp,&
                                             real_4,&
                                             real_8

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_block_access'

  PUBLIC :: dbcsr_get_block, dbcsr_get_block_p,&
            dbcsr_put_block, dbcsr_remove_block

  PUBLIC :: dbcsr_reserve_block2d,&
            dbcsr_reserve_blocks, dbcsr_reserve_all_blocks, dbcsr_reserve_diag_blocks


  INTERFACE dbcsr_get_block
     MODULE PROCEDURE dbcsr_get_block_d, dbcsr_get_block_s,&
                      dbcsr_get_block_z, dbcsr_get_block_c
     MODULE PROCEDURE dbcsr_get_2d_block_d, dbcsr_get_2d_block_s,&
                      dbcsr_get_2d_block_z, dbcsr_get_2d_block_c
  END INTERFACE

  INTERFACE dbcsr_get_block_p
     MODULE PROCEDURE dbcsr_get_block_p_d, dbcsr_get_block_p_s,&
                      dbcsr_get_block_p_z, dbcsr_get_block_p_c
     MODULE PROCEDURE dbcsr_get_2d_block_p_d, dbcsr_get_2d_block_p_s,&
                      dbcsr_get_2d_block_p_z, dbcsr_get_2d_block_p_c
     MODULE PROCEDURE dbcsr_get_block_p_area
  END INTERFACE

  INTERFACE dbcsr_put_block
     MODULE PROCEDURE dbcsr_put_block_area
     MODULE PROCEDURE dbcsr_put_block_d, dbcsr_put_block_s,&
                      dbcsr_put_block_z, dbcsr_put_block_c
     MODULE PROCEDURE dbcsr_put_block2d_d, dbcsr_put_block2d_s,&
                      dbcsr_put_block2d_z, dbcsr_put_block2d_c
  END INTERFACE

  INTERFACE dbcsr_reserve_block2d
     MODULE PROCEDURE dbcsr_reserve_block2d_s, dbcsr_reserve_block2d_d,&
                      dbcsr_reserve_block2d_c, dbcsr_reserve_block2d_z
  END INTERFACE

  INTERFACE dbcsr_get_tmp_block2d
     MODULE PROCEDURE dbcsr_get_tmp_block2d_s, dbcsr_get_tmp_block2d_d,&
                      dbcsr_get_tmp_block2d_c, dbcsr_get_tmp_block2d_z
  END INTERFACE

  INTERFACE dbcsr_set_block_pointer
     MODULE PROCEDURE dbcsr_set_block_pointer_any
     MODULE PROCEDURE dbcsr_set_block_pointer_2d_s,&
                      dbcsr_set_block_pointer_2d_d,&
                      dbcsr_set_block_pointer_2d_c,&
                      dbcsr_set_block_pointer_2d_z
  END INTERFACE



  LOGICAL, PARAMETER :: careful_mod = .FALSE.
  LOGICAL, PARAMETER :: debug_mod = .FALSE.


  INTEGER, PARAMETER, PRIVATE :: rpslot_owner = 1
  INTEGER, PARAMETER, PRIVATE :: rpslot_addblks = 2
  INTEGER, PARAMETER, PRIVATE :: rpslot_addoffset = 3
  INTEGER, PARAMETER, PRIVATE :: rpslot_oldblks = 4
  INTEGER, PARAMETER, PRIVATE :: rpslot_oldoffset = 5
  INTEGER, PARAMETER, PRIVATE :: rpslot_totaloffset = 6
  INTEGER, PARAMETER, PRIVATE :: rpnslots = 6


  LOGICAL, PARAMETER, PRIVATE :: detailed_timing = .FALSE.

  TYPE block_parameters
     LOGICAL :: tr
     INTEGER :: logical_rows, logical_cols
     INTEGER :: offset, nze
  END TYPE block_parameters

  TYPE dgemm_join
     INTEGER :: p_a, p_b, p_c
     INTEGER :: last_k, last_n
     TYPE(dbcsr_scalar_type) :: alpha, beta
  END TYPE dgemm_join

CONTAINS

! *****************************************************************************
!> \brief Marks a block for removal from a DBCSR matrix. Handles
!>        symmetric matrices.
!> \param[in]  matrix         DBCSR matrix
!> \param[in]  row            row of block to remove
!> \param[in]  col            column of block to remove
!> \param[in]  block_number   (optional) the block number, if it is known
! *****************************************************************************
  SUBROUTINE dbcsr_remove_block(matrix, row, col, block_nze, block_number)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col, block_nze
    INTEGER, INTENT(IN), OPTIONAL            :: block_number

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_remove_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: b, c, error_handle, r
    LOGICAL                                  :: found, tr
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set (routineN, error_handle, error)
    IF (PRESENT (block_number)) THEN
       b = block_number
       CALL dbcsr_assert (block_number .LE. matrix%m%nblks, dbcsr_failure_level,&
            dbcsr_caller_error, routineN, "Block number too big.",__LINE__,error)
       found = .TRUE.
    ELSE
       CALL dbcsr_get_block_index (matrix, row, col, r, c, tr, found, b)
    ENDIF
    b = ABS (b)
    IF (found .AND. b .GT. 0) THEN
       ! Mark the block for deletion.
       matrix%m%blk_p(b) = 0
       matrix%m%valid = .FALSE.
       ! update nze accordingly
       matrix%m%nze = matrix%m%nze - block_nze
       IF (debug_mod) THEN
          CALL dbcsr_assert (matrix%m%nze, "GE", 0, dbcsr_failure_level,&
             dbcsr_caller_error, routineN, "nze < 0!",__LINE__,error)
       ENDIF
    ELSE
       IF (debug_mod) THEN
          CALL dbcsr_assert (b, "NE", 0,&
               dbcsr_warning_level, dbcsr_caller_error, routineN,&
               "Block does not exist or is already deleted.",__LINE__,error)
       ENDIF
    ENDIF
    IF (careful_mod) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_remove_block

! *****************************************************************************
!> \brief Gets a block from a dbcsr matrix as a data area
!> \par Data area
!>      The pointer encapsulated in the data area points to data stored in the
!>      matrix. It must be 2-dimensional.
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get
!> \param[in] tr      whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE dbcsr_get_block_p_area(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: block
    LOGICAL, INTENT(OUT)                     :: tr, found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_p_area', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, csize, error_handle, iw, &
                                                offset, rsize, stored_col, &
                                                stored_row
    LOGICAL                                  :: stored_tr
    TYPE(btree_2d_data_c)                    :: data_block_c
    TYPE(btree_2d_data_d)                    :: data_block_d
    TYPE(btree_2d_data_s)                    :: data_block_s
    TYPE(btree_2d_data_z)                    :: data_block_z
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set (routineN, error_handle, error)
    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)

    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    CALL dbcsr_data_clear_pointer (block)
    IF(found) THEN
       CALL dbcsr_set_block_pointer (matrix, block, stored_row, stored_col,&
            rsize, csize, stored_tr, offset, buffer_tr=stored_tr,&
            contiguous_pointers=.TRUE., error=error)
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       iw = 1
!$     iw = omp_get_thread_num()+1
       CALL dbcsr_assert (dbcsr_use_mutable (matrix%m), dbcsr_failure_level,&
            dbcsr_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.",__LINE__,error)
       IF (dbcsr_mutable_instantiated(matrix%m%wms(iw)%mutable)) THEN
          SELECT CASE (block%d%data_type)
          CASE (dbcsr_type_real_4_2d)
             CALL btree_get_s (&
                  matrix%m%wms(iw)%mutable%m%btree_s,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block_s, found)
             IF (found) THEN
                CALL dbcsr_data_set_pointer (block, data_block_s%p)
             ENDIF
          CASE (dbcsr_type_real_8_2d)
             CALL btree_get_d (&
                  matrix%m%wms(iw)%mutable%m%btree_d,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block_d, found)
             IF (found) THEN
                CALL dbcsr_data_set_pointer (block, data_block_d%p)
             ENDIF
          CASE (dbcsr_type_complex_4_2d)
             CALL btree_get_c (&
                  matrix%m%wms(iw)%mutable%m%btree_c,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block_c, found)
             IF (found) THEN
                CALL dbcsr_data_set_pointer (block, data_block_c%p)
             ENDIF
          CASE (dbcsr_type_complex_8_2d)
             CALL btree_get_z (&
                  matrix%m%wms(iw)%mutable%m%btree_z,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block_z, found)
             IF (found) THEN
                CALL dbcsr_data_set_pointer (block, data_block_z%p)
             ENDIF
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
                  routineN, "Only 2-D data for block pointers!",__LINE__,error)
          END SELECT
       ENDIF
    ENDIF
    IF (careful_mod) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_get_block_p_area


! *****************************************************************************
!> \brief
!>          We allow :
!>                  matrix(dp) [+]= [scale(dp)] * block(dp)
!>                  matrix(dp) [+]= [scale(dp)] * block(sp)
!>                  matrix(sp) [+]= [scale(dp)] * block(sp)
!> \param[in]
!> \param[out]
!>
! *****************************************************************************
  SUBROUTINE dbcsr_put_block_area(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    TYPE(dbcsr_data_obj)                     :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block_area', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_type_m, error_handle
    LOGICAL                                  :: do_scale
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set (routineN, error_handle, error=error)
    data_type_m = dbcsr_get_data_type (matrix)
    do_scale = PRESENT (scale)
    IF (do_scale) THEN
       !CALL dbcsr_assert (data_type_m .EQ. scale%data_type, dbcsr_fatal_level,&
       !     dbcsr_wrong_args_error, routineN, "Incompatible data types matrix="//&
       !     data_type_m//" scale="//scale%data_type)
    ENDIF
    CALL dbcsr_assert (ASSOCIATED (block%d), dbcsr_fatal_level, dbcsr_wrong_args_error,&
         routineN, "Can only add valid data block!",__LINE__,error)
    SELECT CASE(block%d%data_type)
    CASE (dbcsr_type_real_4)
       IF (do_scale) THEN
          IF(data_type_m.EQ.dbcsr_type_real_4) THEN
             CALL dbcsr_put_block(matrix, row, col, block%d%r_sp, transposed,&
                  summation, scale=scale%r_sp)
          ELSEIF(data_type_m.EQ.dbcsr_type_real_8) THEN
             CALL dbcsr_put_block(matrix, row, col, REAL(block%d%r_sp,real_8), transposed,&
                  summation, scale=REAL(scale%r_sp,real_8))
          ENDIF
       ELSE
          IF(data_type_m.EQ.dbcsr_type_real_4) THEN
             CALL dbcsr_put_block(matrix, row, col, block%d%r_sp, transposed,&
                  summation)
          ELSEIF(data_type_m.EQ.dbcsr_type_real_8) THEN
             CALL dbcsr_put_block(matrix, row, col, REAL(block%d%r_sp,real_8), transposed,&
                  summation)
          ENDIF
       ENDIF
    CASE (dbcsr_type_real_8)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%r_dp, transposed,&
               summation, scale=scale%r_dp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%r_dp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_complex_4)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%c_sp, transposed,&
               summation, scale=scale%c_sp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%c_sp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_complex_8)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%c_dp, transposed,&
               summation, scale=scale%c_dp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%c_dp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_real_4_2d)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%r2_sp, transposed,&
               summation, scale=scale%r_sp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%r2_sp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_real_8_2d)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%r2_dp, transposed,&
               summation, scale=scale%r_dp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%r2_dp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_complex_4_2d)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%c2_sp, transposed,&
               summation, scale=scale%c_sp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%c2_sp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_complex_8_2d)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%c2_dp, transposed,&
               summation, scale=scale%c_dp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%c2_dp, transposed,&
               summation)
       ENDIF
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid data type",__LINE__,error)
    END SELECT
    IF (careful_mod) CALL dbcsr_error_stop (error_handle, error=error)
  END SUBROUTINE dbcsr_put_block_area

! *****************************************************************************
!> \brief Inserts all blocks of a dbcsr matrix to make it a full matrix.
!>        Thus obviously not linear scaling.
!> \param[in,out] matrix      Matrix into which blocks should be added.
!> \param[in,out] error       Error
! *****************************************************************************
  SUBROUTINE dbcsr_reserve_all_blocks(matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_all_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_count, col, col_s, &
                                                myrank, rank, row, row_s
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: columns, rows
    LOGICAL                                  :: tr

    myrank = dbcsr_mp_mynode (dbcsr_distribution_mp (dbcsr_distribution (matrix)))

    blk_count=0
    ! should be possible to loop only over the local blockrows/blockcols
    DO row = 1, dbcsr_nblkrows_total(matrix)
     DO col = 1, dbcsr_nblkcols_total(matrix)
        tr = .FALSE.
        row_s=row ; col_s=col
        CALL dbcsr_get_stored_coordinates (matrix, row_s, col_s, tr, rank)
        ! is that the correct condition for symmetric matrices ?
        IF (rank.EQ.myrank .AND. row_s.EQ.row .AND. col_s.EQ.col) blk_count=blk_count+1
     ENDDO
    ENDDO

    ALLOCATE(rows(blk_count),columns(blk_count))

    blk_count=0
    DO row = 1, dbcsr_nblkrows_total(matrix)
     DO col = 1, dbcsr_nblkcols_total(matrix)
        tr = .FALSE.
        row_s=row ; col_s=col
        CALL dbcsr_get_stored_coordinates (matrix, row_s, col_s, tr, rank)
        IF (rank.EQ.myrank .AND. row_s.EQ.row .AND. col_s.EQ.col) THEN
           blk_count=blk_count+1 
           rows(blk_count)=row
           columns(blk_count)=col
        ENDIF
     ENDDO
    ENDDO

    CALL dbcsr_reserve_blocks(matrix,rows,columns,error=error)

  END SUBROUTINE dbcsr_reserve_all_blocks

! *****************************************************************************
!> \brief Inserts diagonal blocks of a dbcsr matrix to make it a matrix with at least all diagonal blocks present
!> \param[in,out] matrix      Matrix into which blocks should be added.
!> \param[in,out] error       Error
! *****************************************************************************
  SUBROUTINE dbcsr_reserve_diag_blocks(matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_diag_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_count, col, col_s, &
                                                myrank, rank, row, row_s
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: columns, rows
    LOGICAL                                  :: tr

    myrank = dbcsr_mp_mynode (dbcsr_distribution_mp (dbcsr_distribution (matrix)))

    blk_count=0
    ! should be possible to loop only over the local blockrows/blockcols
    DO row = 1, dbcsr_nblkrows_total(matrix)
       col = row
       tr = .FALSE.
       row_s=row ; col_s=col
       CALL dbcsr_get_stored_coordinates (matrix, row_s, col_s, tr, rank)
       IF (rank.EQ.myrank .AND. row_s.EQ.row .AND. col_s.EQ.col) blk_count=blk_count+1
    ENDDO

    ALLOCATE(rows(blk_count),columns(blk_count))

    blk_count=0
    DO row = 1, dbcsr_nblkrows_total(matrix)
       col = row
       tr = .FALSE.
       row_s=row ; col_s=col
       CALL dbcsr_get_stored_coordinates (matrix, row_s, col_s, tr, rank)
       IF (rank.EQ.myrank .AND. row_s.EQ.row .AND. col_s.EQ.col) THEN
          blk_count=blk_count+1
          rows(blk_count)=row
          columns(blk_count)=col
       ENDIF
    ENDDO

    CALL dbcsr_reserve_blocks(matrix,rows,columns,error=error)

  END SUBROUTINE dbcsr_reserve_diag_blocks

! *****************************************************************************
!> \brief Inserts block reservations into a matrix, avoiding the work matrix.
!> \par Data
!>      No data can be specified; instead, space is reserved and zeroed. To
!>      add data, call dbcsr_put_block afterwards.
!> \par Reserving existing blocks
!>      Duplicates are not added, but allocations may be greater than
!>      the minimum necessary.
!> \par blk_pointers
!>      When blk_pointers is passed, the newly added blocks use these pointers.
!>      No data is cleared in this case
!> \param[in,out] matrix      Matrix into which blocks should be added.
!> \param[in] rows            Rows of the blocks to add
!> \param[in] columns         Columns of the blocks to add
!> \param[in] blk_pointers    (optional) block pointers to use for new blocks
!> \param[in,out] error       Error
! *****************************************************************************
  SUBROUTINE dbcsr_reserve_blocks(matrix, rows, columns, blk_pointers, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, DIMENSION(:), INTENT(IN)        :: rows, columns
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: blk_pointers
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reserve_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, blk_p, data_size_new, data_size_old, error_handler, &
      nblkrows, nblks_actual_added, nblks_added, nblks_new, nblks_old, &
      new_data_sizes, nze
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: add_blkp, add_cols, add_rows, &
                                                added_sizes, new_blk_p, &
                                                new_col_i, new_row_i, &
                                                old_row_i
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: added_blk_info

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    CALL dbcsr_assert (SIZE(rows), "EQ", SIZE(columns), dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN,&
         "Size of rows and columns array must match.", __LINE__, error=error)
    IF (PRESENT (blk_pointers)) THEN
       CALL dbcsr_assert (SIZE(rows), "EQ", SIZE(blk_pointers),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Size of rows and block pointecs arrays must match.",&
            __LINE__, error=error)
       data_size_old = 0
    ELSE
       ! Get current data size
       data_size_old = dbcsr_get_data_size_referenced (matrix)
    ENDIF
    !> Ensures that the rows and columns are sorted.
    nblks_added = SIZE(rows)
    ALLOCATE (add_rows (nblks_added))
    add_rows(:) = rows(:)
    ALLOCATE (add_cols (nblks_added))
    add_cols(:) = columns(:)
    IF (PRESENT(blk_pointers)) THEN
       ALLOCATE (add_blkp (nblks_added))
       add_blkp(:) = blk_pointers(:)
       CALL dbcsr_sort_indices (nblks_added, add_rows, add_cols,&
            blk_p = add_blkp)
    ELSE
       CALL dbcsr_sort_indices (nblks_added, add_rows, add_cols)
    ENDIF
    nblks_old = dbcsr_get_num_blocks (matrix)
    nblkrows = dbcsr_nblkrows_total(matrix)
    CALL dbcsr_assert (SIZE(rows) .GT. 0, "IMP",&
         nblkrows .GT. 0, dbcsr_fatal_level,&
         dbcsr_internal_error, routineN,&
         "Can not add blocks to matrix with no rows.", __LINE__, error=error)
    ! Adjust the index.
    ! Get the old row indices
    ALLOCATE (old_row_i (nblks_old))
    CALL dbcsr_expand_row_index(matrix%m%row_p, old_row_i,&
         nblkrows, nblks_old)
    ! Calculate new block pointers. Possibly high estimates.
    new_data_sizes = 0
    blk_p = data_size_old + 1   ! New blocks start at the end of the old
    ALLOCATE (added_blk_info (3, nblks_added))
    ALLOCATE (added_sizes (nblks_added))
    DO blk = 1, nblks_added
       IF (PRESENT (blk_pointers)) THEN
          blk_p = add_blkp(blk)
       ENDIF
       added_blk_info(1:3,blk) = (/ add_rows(blk), add_cols(blk), blk_p /)
       nze = dbcsr_blk_row_size (matrix, add_rows(blk)) &
            * dbcsr_blk_column_size (matrix, add_cols(blk))
       added_sizes(blk) = nze
       blk_p = blk_p + nze
    ENDDO
    DEALLOCATE (add_rows)
    DEALLOCATE (add_cols)
    IF (PRESENT (blk_pointers)) DEALLOCATE (add_blkp)
    !
    nblks_new = nblks_old + nblks_added ! Possibly high estimate
    ALLOCATE (new_row_i (nblks_new))
    ALLOCATE (new_col_i (nblks_new))
    ALLOCATE (new_blk_p (nblks_new))
    ! Merge the two indices
    IF (PRESENT (blk_pointers)) THEN
       CALL merge_index_arrays (new_row_i, new_col_i, new_blk_p, nblks_new,&
            old_row_i, matrix%m%col_i, matrix%m%blk_p, nblks_old,&
            added_blk_info, nblks_added, added_nblks=nblks_actual_added,&
            error=error)
       data_size_new = 0
    ELSE
       CALL merge_index_arrays (new_row_i, new_col_i, new_blk_p, nblks_new,&
            old_row_i, matrix%m%col_i, matrix%m%blk_p, nblks_old,&
            added_blk_info, nblks_added, added_nblks=nblks_actual_added,&
            added_sizes=added_sizes, added_size_offset=data_size_old+1,&
            added_size=data_size_new,&
            error=error)
    ENDIF
    nblks_new = nblks_actual_added + nblks_old
    ! Free some memory
    DEALLOCATE (added_blk_info)
    DEALLOCATE (added_sizes)
    DEALLOCATE (old_row_i)
    ! We can skip this if no block was actually added.
    IF (nblks_actual_added .GT. 0) THEN
       ! Write the new index
       matrix%m%nblks = nblks_new
       matrix%m%nze = matrix%m%nze + data_size_new
       matrix%m%index(dbcsr_slot_nblks) = matrix%m%nblks
       matrix%m%index(dbcsr_slot_nze) = matrix%m%index(dbcsr_slot_nze)
       CALL dbcsr_clearfrom_index_array (matrix%m, dbcsr_slot_col_i)
       CALL dbcsr_clearfrom_index_array (matrix%m, dbcsr_slot_blk_p)
       CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_col_i,&
            new_col_i(1:nblks_new),&
            extra=nblks_new, error=error)
       CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_blk_p,&
            new_blk_p(1:nblks_new),&
            error=error)
       CALL dbcsr_make_dbcsr_index (matrix%m%row_p, new_row_i(1:nblks_new),&
            nblkrows, nblks_new)
       IF (.NOT. PRESENT (blk_pointers)) THEN
          ! Resize data area to fit the new blocks.
          CALL dbcsr_data_ensure_size (matrix%m%data_area,&
               data_size = matrix%m%nze, error=error)
          ! Zero the new data blocks.
          CALL dbcsr_data_clear (matrix%m%data_area,&
               lb=data_size_old+1, ub=matrix%m%nze)
       ENDIF
    ENDIF
    CALL dbcsr_error_stop (error_handler, error)
  END SUBROUTINE dbcsr_reserve_blocks


! *****************************************************************************
!> \brief Sets a pointer, possibly using the buffers.
!> \param[in] matrix           Matrix to use
!> \param[in,out] pointer      The pointer to set
!> \param[in] row, col         Row and column of block to point to
!> \param[in] rsize, csize     Row and column sizes of block to point to
!> \param[in] main_tr          Whether block is transposed in the matrix
!> \param[in] base_offset      The block pointer
!> \param[in] buffer_tr        Whether buffer should be transposed
!> \param[in] contiguous_pointers  (optional) Whether pointers should be made
!>                                 contiguous
!> \param[in] read_only        (optional) User promise not to change data
!> \param[in,out] error        error
! *****************************************************************************
  SUBROUTINE dbcsr_set_block_pointer_any (matrix, pointer_any, row, col,&
       rsize, csize, main_tr, base_offset, buffer_tr, contiguous_pointers,&
       read_only, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: pointer_any
    INTEGER, INTENT(IN)                      :: row, col, rsize, csize
    LOGICAL, INTENT(IN)                      :: main_tr
    INTEGER, INTENT(IN)                      :: base_offset
    LOGICAL, INTENT(IN)                      :: buffer_tr
    LOGICAL, INTENT(IN), OPTIONAL            :: contiguous_pointers, read_only
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_block_pointer_any', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

       IF (main_tr) THEN
          CALL dbcsr_data_set_pointer (pointer_any, csize, rsize,&
               matrix%m%data_area, source_lb = base_offset)
       ELSE
          CALL dbcsr_data_set_pointer (pointer_any, rsize, csize,&
               matrix%m%data_area, source_lb = base_offset)
       ENDIF
  END SUBROUTINE dbcsr_set_block_pointer_any



#include "dbcsr_block_access_d.f90"
#include "dbcsr_block_access_z.f90"
#include "dbcsr_block_access_s.f90"
#include "dbcsr_block_access_c.f90"


END MODULE dbcsr_block_access

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR block access
!> \author  Urban Borstnik
!> \date    2010-02-18
!> \version 0.9
!>
!> <b>Modification history:</b>
!  - 2010-02-18 Moved from dbcsr_operations and dbcsr_methods
! *****************************************************************************
MODULE dbcsr_block_access
  USE array_types,                     ONLY: array_data
  USE btree_I8_k_cp2d_v,               ONLY: btree_2d_data_c => cp2d,&
                                             btree_add_c => btree_add,&
                                             btree_c => btree,&
                                             btree_get_c => btree_find
  USE btree_I8_k_dp2d_v,               ONLY: btree_2d_data_d => dp2d,&
                                             btree_add_d => btree_add,&
                                             btree_d => btree,&
                                             btree_get_d => btree_find
  USE btree_I8_k_sp2d_v,               ONLY: btree_2d_data_s => sp2d,&
                                             btree_add_s => btree_add,&
                                             btree_get_s => btree_find,&
                                             btree_s => btree
  USE btree_I8_k_zp2d_v,               ONLY: btree_2d_data_z => zp2d,&
                                             btree_add_z => btree_add,&
                                             btree_get_z => btree_find,&
                                             btree_z => btree
  USE dbcsr_block_operations,          ONLY: dbcsr_block_copy,&
                                             dbcsr_block_transpose
  USE dbcsr_dist_operations,           ONLY: dbcsr_get_block_index,&
                                             dbcsr_get_stored_block_info,&
                                             dbcsr_get_stored_coordinates
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: real_4,&
                                             real_8
  USE dbcsr_methods,                   ONLY: &
       dbcsr_blk_column_size, dbcsr_blk_row_size, dbcsr_data_clear_pointer, &
       dbcsr_data_ensure_size, dbcsr_data_set_pointer, dbcsr_get_data_p, &
       dbcsr_get_data_type, dbcsr_mutable_instantiated, dbcsr_mutable_new, &
       dbcsr_use_mutable, dbcsr_wm_use_mutable
  USE dbcsr_ptr_util,                  ONLY: addto_array,&
                                             pointer_view
  USE dbcsr_types,                     ONLY: &
       dbcsr_data_obj, dbcsr_obj, dbcsr_scalar_type, dbcsr_type_complex_4, &
       dbcsr_type_complex_4_2d, dbcsr_type_complex_8, &
       dbcsr_type_complex_8_2d, dbcsr_type_real_4, dbcsr_type_real_4_2d, &
       dbcsr_type_real_8, dbcsr_type_real_8_2d
  USE dbcsr_util,                      ONLY: dbcsr_set_debug,&
                                             make_coordinate_tuple,&
                                             swap
  USE dbcsr_work_operations,           ONLY: add_work_coordinate,&
                                             dbcsr_work_create

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_block_access'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034


  PUBLIC :: dbcsr_get_block, dbcsr_get_block_p,&
            dbcsr_put_block, dbcsr_remove_block

  PUBLIC :: dbcsr_reserve_block2d, dbcsr_get_tmp_block2d


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





#define DBG IF (dbg) WRITE(*,*)routineN//" ",
#define DBGV IF (bcsr_verbose) WRITE(*,*)routineN//" ",
#define DBGI IF (info) WRITE(*,*)routineN//" ",
#define DEBUG_HEADER        LOGICAL :: dbg, info
#define DEBUG_BODY        dbg = .FALSE. ; CALL dbcsr_set_debug(dbg, info=info)
#define DEBUG_BODYY       dbg = .TRUE. ; CALL dbcsr_set_debug(dbg, info=info)
  LOGICAL, PARAMETER :: bcsr_debug =   .TRUE.
  LOGICAL, PARAMETER :: bcsr_info =    .FALSE.
  LOGICAL, PARAMETER :: bcsr_verbose = .FALSE.


#define temp_transpose(v, r, c) RESHAPE(TRANSPOSE(RESHAPE(v,(/r,c/))),(/r*c/))

  INTEGER, PARAMETER, PRIVATE :: rpslot_owner = 1
  INTEGER, PARAMETER, PRIVATE :: rpslot_addblks = 2
  INTEGER, PARAMETER, PRIVATE :: rpslot_addoffset = 3
  INTEGER, PARAMETER, PRIVATE :: rpslot_oldblks = 4
  INTEGER, PARAMETER, PRIVATE :: rpslot_oldoffset = 5
  INTEGER, PARAMETER, PRIVATE :: rpslot_totaloffset = 6
  INTEGER, PARAMETER, PRIVATE :: rpnslots = 6


  LOGICAL, PARAMETER, PRIVATE :: detailed_timing = .FALSE.

  INTEGER, PARAMETER :: dgemm_stack_size = 32

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
! Block access functions
! *****************************************************************************

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

    INTEGER                                  :: b, c, r
    LOGICAL                                  :: found, tr
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

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
       CALL dbcsr_assert (matrix%m%nze .GE. 0, dbcsr_failure_level,&
            dbcsr_caller_error, routineN, "nze < 0!",__LINE__,error)
    ELSE
       CALL dbcsr_assert (b .NE. 0, dbcsr_warning_level, dbcsr_caller_error, routineN,&
            "Block does not exist or is already deleted.",__LINE__,error)
    ENDIF
  END SUBROUTINE dbcsr_remove_block


! *****************************************************************************
!> \brief Gets a block from a dbcsr matrix as a data area
!> \par Data area
!>      The pointer encapsulated in the data area points to data stored in the
!>      matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get
!> \param[in] tr      whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
! *****************************************************************************
!> \brief Gets a 2-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-2 array)
!> \param[out] tr     whether the data is transposed
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

    INTEGER                                  :: blk, csize, iw, nze, offset, &
                                                rsize, stored_col, stored_row
    LOGICAL                                  :: stored_tr
    TYPE(btree_2d_data_c)                    :: data_block_c
    TYPE(btree_2d_data_d)                    :: data_block_d
    TYPE(btree_2d_data_s)                    :: data_block_s
    TYPE(btree_2d_data_z)                    :: data_block_z
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)

    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    CALL dbcsr_data_clear_pointer (block)
    IF(found) THEN
       nze = rsize*csize
       !IF(nze.eq.0) THEN
       !   found = .FALSE.
       !ELSE
       CALL dbcsr_data_set_pointer (block, rsize, csize, matrix%m%data_area,&
            source_lb = offset)
       !ENDIF
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
                  routineN, "Invalid data type",__LINE__,error)
          END SELECT
       ENDIF
    ENDIF
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

    CHARACTER                                :: data_type_m
    LOGICAL                                  :: do_scale
    TYPE(dbcsr_error_type)                   :: error

    data_type_m = dbcsr_get_data_type (matrix)
    do_scale = PRESENT (scale)
    IF (do_scale) THEN
       !CALL dbcsr_assert (data_type_m .EQ. scale%data_type, dbcsr_fatal_level,&
       !     dbcsr_wrong_args_error, routineN, "Incompatible data types matrix="//&
       !     data_type_m//" scale="//scale%data_type)
    ENDIF
    CALL dbcsr_assert (ASSOCIATED (block%d), dbcsr_fatal_level, dbcsr_wrong_args_error,&
         routineN, "Can only add valid data block!",__LINE__,error)
    !CALL dbcsr_assert (data_type_m .EQ. dbcsr_type_2d_to_1d(block%d%data_type),&
    !     dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
    !     "Incompatible data types matrix="//&
    !        data_type_m//" block="//block%d%data_type)
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
  END SUBROUTINE dbcsr_put_block_area



#include "dbcsr_block_access_d.F"
#include "dbcsr_block_access_z.F"
#include "dbcsr_block_access_s.F"
#include "dbcsr_block_access_c.F"


END MODULE dbcsr_block_access

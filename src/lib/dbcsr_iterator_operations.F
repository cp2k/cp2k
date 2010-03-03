!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR iterator operations
!> \author  Urban Borstnik
!> \date    2010-02-18
!> \version 0.9
!>
!> <b>Modification history:</b>
!  - 2010-02-18 Moved from dbcsr_operations and dbcsr_methods
! *****************************************************************************
MODULE dbcsr_iterator_operations
  USE array_types,                     ONLY: array_data,&
                                             array_size
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
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: int_8,&
                                             real_4,&
                                             real_8
  USE dbcsr_message_passing,           ONLY: mp_bcast,&
                                             mp_sum
  USE dbcsr_methods,                   ONLY: dbcsr_data_hold,&
                                             dbcsr_data_init,&
                                             dbcsr_data_new,&
                                             dbcsr_data_release,&
                                             dbcsr_data_set_pointer,&
                                             dbcsr_distribution_thread_dist,&
                                             dbcsr_get_data_p,&
                                             dbcsr_get_data_type
  USE dbcsr_ptr_util,                  ONLY: addto_array
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_iterator,&
                                             dbcsr_obj,&
                                             dbcsr_scalar_type,&
                                             dbcsr_type
  USE dbcsr_util,                      ONLY: dbcsr_set_debug,&
                                             swap

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_iterator_operations'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034


  PUBLIC :: dbcsr_iterator_start, dbcsr_iterator_next_block,&
            dbcsr_iterator_stop, dbcsr_iterator_blocks_left,&
            iterator_next_2d_block_d


  INTERFACE dbcsr_iterator_next_block
     MODULE PROCEDURE iterator_next_block_index,&
                      iterator_next_area_block
     MODULE PROCEDURE iterator_next_2d_block_d,&
                      iterator_next_2d_block_s,&
                      iterator_next_2d_block_c,&
                      iterator_next_2d_block_z,&
                      iterator_next_1d_block_d,&
                      iterator_next_1d_block_s,&
                      iterator_next_1d_block_c,&
                      iterator_next_1d_block_z
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
! Iterator functions
! *****************************************************************************


! *****************************************************************************
!> \brief Sets up an iterator
!> \param[out] iterator   the iterator
!> \param[in] matrix      DBCSR matrix
!> \param[in] shared      (optional) Threads do not iterator through
!>                        the same matrix in a multi-threaded
!>                        environment; default is TRUE.
!> \param[in] contiguous_pointers  (optional) Whether returned pointers need
!>                                 to be contiguous; default is FALSE.
!>
!> Contiguous pointers may incur reallocation penalties but enable quick
!> passing of arrays to routines with unspecified interfaces (i.e., direct
!> calls to BLACS or MPI).
! *****************************************************************************
  SUBROUTINE dbcsr_iterator_start (iterator, matrix, shared,&
       contiguous_pointers)
    TYPE(dbcsr_iterator), INTENT(OUT)        :: iterator
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: shared, contiguous_pointers

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_iterator_start', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ithread, max_c_size, &
                                                max_r_size, max_size
    LOGICAL                                  :: share

!   ---------------------------------------------------------------------------

    share = .TRUE.
    IF (PRESENT (shared)) share = shared
    share = share .AND.&
         array_size (dbcsr_distribution_thread_dist (matrix%m%dist)) .GT. 0
    iterator%shared = share
    IF (PRESENT (contiguous_pointers)) THEN
       iterator%contiguous_pointers = contiguous_pointers
    ELSE
       iterator%contiguous_pointers = .FALSE.
    ENDIF
    iterator%row = 0
    iterator%pos = 0
    iterator%rbs => array_data (matrix%m%row_blk_size)
    iterator%cbs => array_data (matrix%m%col_blk_size)
    iterator%roff => array_data (matrix%m%row_blk_offset)
    iterator%coff => array_data (matrix%m%col_blk_offset)

    iterator%nblks = matrix%m%nblks
    iterator%nblkrows_total = matrix%m%nblkrows_total
    iterator%max_r_size = matrix%m%max_rbs
    iterator%max_c_size = matrix%m%max_cbs

    ! Can't do that because fortran forces us lose the matrix reference.
    !CALL dbcsr_hold (matrix)
    iterator%row_p => matrix%m%row_p
    iterator%col_i => matrix%m%col_i
    iterator%blk_p => matrix%m%blk_p
    iterator%data_area = matrix%m%data_area
    CALL dbcsr_data_hold (iterator%data_area)
    iterator%row_size = 0
    IF (iterator%shared) THEN
       iterator%tdist => array_data (dbcsr_distribution_thread_dist (matrix%m%dist))
    ELSE
       NULLIFY (iterator%tdist)
    ENDIF
    IF (matrix%m%nblks .GE. 1) THEN
       iterator%row = 1
       iterator%pos = 1
       IF (.NOT. iterator%shared) THEN
          CALL find_proper_position(iterator%pos,&
               iterator%row, iterator%nblks, iterator%nblkrows_total,&
               iterator%blk_p, iterator%row_p)
       ELSE
          ithread = 0
!$        ithread = OMP_GET_THREAD_NUM()
          CALL find_proper_position(iterator%pos,&
               iterator%row, iterator%nblks, iterator%nblkrows_total,&
               iterator%blk_p, iterator%row_p,&
               tdist=iterator%tdist, tid=ithread)
       ENDIF
       IF (iterator%row .GT. 0) THEN
          iterator%row_size = iterator%rbs(iterator%row)
          iterator%row_offset = iterator%roff(iterator%row)
          max_r_size = iterator%max_r_size
          max_c_size = iterator%max_c_size
       ENDIF
    ENDIF
    max_size = MAX (matrix%m%max_rbs, matrix%m%max_cbs)
    CALL dbcsr_data_init(iterator%buffer_2d)
    CALL dbcsr_data_new(iterator%buffer_2d,&
         dbcsr_get_data_type(matrix),&
         max_size, max_size)
  END SUBROUTINE dbcsr_iterator_start

! *****************************************************************************
!> \brief Stops up an iterator
!> \param[out] iterator   the iterator
!> \param[in] matrix      DBCSR matrix
! *****************************************************************************
  SUBROUTINE dbcsr_iterator_stop (iterator)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_iterator_stop', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------
!CALL dbcsr_release (iterator%matrix)

    iterator%row = 0
    iterator%pos = 0

    CALL dbcsr_data_release (iterator%buffer_2d)
    NULLIFY (iterator%tdist)
    !CALL dbcsr_release (iterator%matrix)
    CALL dbcsr_data_release (iterator%data_area)
  END SUBROUTINE dbcsr_iterator_stop


! *****************************************************************************
!> \brief Finds the first valid block, inclusive from the current position.
!>
!> If there is no valid block, pos is set to 0
!> \param[in,out] pos        input: current position; output: next valid
!>                           position or 0
!> \param[in] maxpos         maximal allowed position
!> \param[in] blk_p          block pointers, used to check validity
! *****************************************************************************
  PURE SUBROUTINE find_first_valid_block (pos, maxpos, blk_p)
    INTEGER, INTENT(INOUT)                   :: pos
    INTEGER, INTENT(IN)                      :: maxpos
    INTEGER, DIMENSION(1:maxpos), INTENT(IN) :: blk_p

    CHARACTER(len=*), PARAMETER :: routineN = 'find_first_valid_block', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF (pos .LT. 1) pos = 1
    DO WHILE (pos .LE. maxpos)
       IF (blk_p(pos) .EQ. 0) THEN
          !WRITE(*,*)routineN//" Skipping deleted block."
          pos = pos+1
       ELSE
          EXIT
       ENDIF
    ENDDO
    IF (pos .GT. maxpos) pos = 0
  END SUBROUTINE find_first_valid_block

! *****************************************************************************
!> \brief Finds the row to which the current block belongs
!>
!> If there is no valid block, pos is set to 0
!> \param[in pos             current position
!> \param[in,out] row        input: current row; output: the row corresponding
!>                           to the position
!> \param[in] maxrows        maxmimum row
!> \param[in] row_p          row pointers
! *****************************************************************************
  PURE SUBROUTINE find_proper_row (pos, row, maxrows, row_p)
    INTEGER, INTENT(IN)                      :: pos
    INTEGER, INTENT(INOUT)                   :: row
    INTEGER, INTENT(IN)                      :: maxrows
    INTEGER, DIMENSION(1:maxrows+1), &
      INTENT(IN)                             :: row_p

    CHARACTER(len=*), PARAMETER :: routineN = 'find_proper_row', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF (pos .GT. 0) THEN
       IF (row .LT. 1) THEN
          row = 1
       ELSEIF (row .GT. maxrows) THEN
          row = maxrows
       ENDIF
       DO WHILE (row_p(row+1) .LT. pos)
          row = row+1
          IF (row .GT. maxrows) THEN
             row = 0
             EXIT
          ENDIF
       ENDDO
    ELSE
       row = 0
    ENDIF
  END SUBROUTINE find_proper_row

! *****************************************************************************
!> \brief Finds the next proper position accounting for threads
!>
!> First time: pos and row are set to 0.
!> If there is no valid block, pos is set to 0
!> \param[inout] pos         current position and updated position
!> \param[in,out] row        input: current row; output: the row corresponding
!>                           to the next proper position
!> \param[in] maxpos         maximum allowable position
!> \param[in] maxrows        maxmimum row
!> \param[in] blk_p          block pointercs
!> \param[in] row_p          row pointers
!> \param[in] tdist          (optional) thread distribution
!> \param[in] tid            my thread number
! *****************************************************************************
  PURE SUBROUTINE find_proper_position (pos, row, maxpos, maxrows,&
       blk_p, row_p, tdist, tid)
    INTEGER, INTENT(INOUT)                   :: pos, row
    INTEGER, INTENT(IN)                      :: maxpos, maxrows
    INTEGER, DIMENSION(1:maxpos), INTENT(IN) :: blk_p
    INTEGER, DIMENSION(1:maxrows+1), &
      INTENT(IN)                             :: row_p
    INTEGER, DIMENSION(1:maxrows), &
      INTENT(IN), OPTIONAL                   :: tdist
    INTEGER, INTENT(IN), OPTIONAL            :: tid

    CHARACTER(len=*), PARAMETER :: routineN = 'find_proper_position', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: row_inrange, row_ok

!   ---------------------------------------------------------------------------

    IF (maxpos .GE. 1) THEN
       IF (pos.EQ.0) pos = 1
       CALL find_first_valid_block (pos, maxpos, blk_p)
       CALL find_proper_row(pos, row, maxrows, row_p)
       row_inrange = row .NE. 0 .AND. row .LE. maxrows
       row_ok = row_inrange
!$     IF (present (tdist) .AND. PRESENT (tid) .AND. row_inrange) THEN
!$        row_ok = tdist(row) .EQ. tid
!$     ENDIF
       DO WHILE (row_inrange .AND. .NOT. row_ok)
          row = row + 1
          pos = row_p(row)+1
          IF (row .GT. maxrows) THEN
             row = 0
             EXIT
          ENDIF
          CALL find_first_valid_block (pos, maxpos, blk_p)
          CALL find_proper_row(pos, row, maxrows, row_p)
          row_inrange = row .NE. 0
          row_ok = row_inrange
!$     IF (present (tdist) .AND. PRESENT (tid) .AND. row_inrange) THEN
!$        row_ok = tdist(row) .EQ. tid
!$     ENDIF
       ENDDO
       IF (row .EQ. 0) pos = 0
    ELSE
       pos = 0
       row = 0
    ENDIF
  END SUBROUTINE find_proper_position

! *****************************************************************************
!> \brief Gets the index information of the next block, no data.
!> \param[in,out] iterator   the iterator
!> \param[out] row           row of the data block
!> \param[out] column        column of the data block
!> \param[out] blk           block number
!> \param[out] transposed    (optional) whether block is transposed
!> \param[out] blk_p         (optional) index into block data array
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
!> \param[out] row_offset    (optional) logical row offset of block
!> \param[out] col_offset    (optional) logical column offset of block
! *****************************************************************************
  SUBROUTINE iterator_next_block_index (iterator, row, column, blk,&
       transposed, blk_p, row_size, col_size, row_offset, col_offset)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator
    INTEGER, INTENT(OUT)                     :: row, column, blk
    LOGICAL, INTENT(OUT), OPTIONAL           :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: blk_p, row_size, col_size, &
                                                row_offset, col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'iterator_next_block_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bp, ithread

!   ---------------------------------------------------------------------------

    IF (iterator%pos .LE. iterator%nblks&
         .AND. iterator%pos .NE. 0) THEN
       row = iterator%row
       column = iterator%col_i(iterator%pos)
       blk = iterator%pos
       IF (PRESENT (row_size)) row_size = iterator%row_size
       IF (PRESENT (col_size)) col_size = iterator%cbs(column)
       IF (PRESENT (row_offset)) row_offset = iterator%row_offset
       IF (PRESENT (col_offset)) col_offset = iterator%coff(column)
       IF (PRESENT (blk_p) .OR. PRESENT (transposed)) THEN
          bp = iterator%blk_p(iterator%pos)
          IF (PRESENT (blk_p)) blk_p = bp
          IF (PRESENT (transposed)) transposed = bp .LT. 0
       ENDIF
       iterator%pos = iterator%pos+1
       IF (.NOT.iterator%shared) THEN
          CALL find_proper_position(iterator%pos,&
               iterator%row, iterator%nblks, iterator%nblkrows_total,&
               iterator%blk_p, iterator%row_p)
       ELSE
          ithread = 0
!$        ithread = OMP_GET_THREAD_NUM()
          CALL find_proper_position(iterator%pos,&
               iterator%row, iterator%nblks, iterator%nblkrows_total,&
               iterator%blk_p, iterator%row_p,&
               tdist=iterator%tdist, tid=ithread)
       ENDIF
       !CALL find_first_valid_block (iterator%pos, iterator%nblks,&
       !     iterator%blk_p)
       !CALL find_proper_row (iterator%pos, iterator%row,&
       !     iterator%nblkrows_total, iterator%row_p)
       IF (iterator%row .GT. 0) iterator%row_size = iterator%rbs(iterator%row)
       IF (iterator%row .GT. 0) iterator%row_offset = iterator%roff(iterator%row)
    ELSE
       row = 0
       column = 0
    ENDIF
    
  END SUBROUTINE iterator_next_block_index

! *****************************************************************************
!> \brief Gets the next data block encapsulated in an object.
!> \param[in,out] iterator   the iterator
!> \param[out] row           row of the data block
!> \param[out] column        column of the data block
!> \param[out] block         encapsulated data
!> \param[out] transposed    whether the block data is transposed
!> \param[out] block_number  (optional) block number
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
!> \param[out] row_offset    (optional) logical row offset of block
!> \param[out] col_offset    (optional) logical column offset of block
! *****************************************************************************
  SUBROUTINE iterator_next_area_block (iterator, row, column, block,&
       transposed, block_number, row_size, col_size, row_offset, col_offset)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: block
    LOGICAL, INTENT(OUT)                     :: transposed
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number, row_size, &
                                                col_size, row_offset, &
                                                col_offset

    CHARACTER(len=*), PARAMETER :: routineN = 'iterator_next_area_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_p, bp, csize, ithread, &
                                                nze, rsize

!   ---------------------------------------------------------------------------

    IF (iterator%pos .LE. iterator%nblks&
         .AND. iterator%pos .NE. 0) THEN
       row = iterator%row
       column = iterator%col_i(iterator%pos)
       blk_p = iterator%blk_p(iterator%pos)
       transposed = blk_p .LT. 0
       bp = ABS (blk_p)
       rsize = iterator%row_size
       csize = iterator%cbs(column)
       nze = rsize * csize
       IF (PRESENT (row_size)) row_size = rsize
       IF (PRESENT (col_size)) col_size = csize
       IF (PRESENT (row_offset)) row_offset = iterator%row_offset
       IF (PRESENT (col_offset)) col_offset = iterator%coff(column)
       ! Redirect the encapsulated pointer to the correct pointer here.
       IF (transposed) CALL swap (rsize, csize)
       CALL dbcsr_data_set_pointer (block, rsize, csize, iterator%data_area,&
            source_lb=bp)
       IF (PRESENT (block_number)) block_number = iterator%pos
       ! Move to the next non-deleted position.
       iterator%pos = iterator%pos+1
       IF (.NOT. iterator%shared) THEN
          CALL find_proper_position(iterator%pos,&
               iterator%row, iterator%nblks, iterator%nblkrows_total,&
               iterator%blk_p, iterator%row_p)
       ELSE
          ithread = 0
!$        ithread = OMP_GET_THREAD_NUM()
          CALL find_proper_position(iterator%pos,&
               iterator%row, iterator%nblks, iterator%nblkrows_total,&
               iterator%blk_p, iterator%row_p,&
               tdist=iterator%tdist, tid=ithread)
       ENDIF
       IF (iterator%row .GT. 0) iterator%row_size = iterator%rbs(iterator%row)
    ELSE
       row = 0
       column = 0
       IF (PRESENT (block_number)) block_number = 0
    ENDIF
  END SUBROUTINE iterator_next_area_block



! *****************************************************************************
!> \brief Returns whether there any blocks left in the iterator.
!> \param[in] iterator       the iterator
!> \param[in]
! *****************************************************************************
  PURE FUNCTION dbcsr_iterator_blocks_left (iterator) RESULT (blocks_left)
    TYPE(dbcsr_iterator), INTENT(IN)         :: iterator
    LOGICAL                                  :: blocks_left

    blocks_left = iterator%pos .NE. 0
  END FUNCTION dbcsr_iterator_blocks_left


! *****************************************************************************
!> \brief Calls the functions that are called in the included files but not
!>        in thin one. Needed to make prettify think the functions ary really
!>        used.
!> \param[in] matrix          matrix
!> \param[out] temp1, temp2   pointers to the data
! *****************************************************************************
  SUBROUTINE prettify_fooler (matrix, temp1, temp2, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    REAL(KIND=real_4), DIMENSION(:), POINTER :: temp1
    REAL(KIND=real_8), DIMENSION(:), POINTER :: temp2
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    INTEGER                                  :: i
    INTEGER(KIND=int_8)                      :: key
    TYPE(btree_2d_data_c)                    :: bt_dc
    TYPE(btree_2d_data_d)                    :: bt_dd
    TYPE(btree_2d_data_s)                    :: bt_ds
    TYPE(btree_2d_data_z)                    :: bt_dz
    TYPE(btree_c)                            :: bt_c
    TYPE(btree_d)                            :: bt_d
    TYPE(btree_s)                            :: bt_s
    TYPE(btree_z)                            :: bt_z

    IF (.FALSE.) THEN
       temp1 => dbcsr_get_data_p (matrix%data_area, 0.0_real_4)
       temp2 => dbcsr_get_data_p (matrix%data_area, 0.0_real_8)
       CALL mp_sum (i, 0)
       CALL mp_bcast (temp1, 0, 0)
       CALL addto_array (matrix%data_area, temp1, error=error)
       CALL btree_add_s (bt_s, key, bt_ds)
       CALL btree_add_d (bt_d, key, bt_dd)
       CALL btree_add_c (bt_c, key, bt_dc)
       CALL btree_add_z (bt_z, key, bt_dz)
       CALL btree_get_s (bt_s, key, bt_ds)
       CALL btree_get_d (bt_d, key, bt_dd)
       CALL btree_get_c (bt_c, key, bt_dc)
       CALL btree_get_z (bt_z, key, bt_dz)
    ENDIF
  END SUBROUTINE prettify_fooler


#include "dbcsr_iterator_operations_d.F"
#include "dbcsr_iterator_operations_z.F"
#include "dbcsr_iterator_operations_s.F"
#include "dbcsr_iterator_operations_c.F"


END MODULE dbcsr_iterator_operations

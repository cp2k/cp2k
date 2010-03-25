!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Operations on the DBCSR index
!> \author  Urban Borstnik
!> \date    2010-02-18
!> \version 0.9
!>
!> <b>Modification history:</b>
!> - Created 2010-02-18
! *****************************************************************************
MODULE dbcsr_index_operations

  USE array_types,                     ONLY: array_data
  USE dbcsr_dist_operations,           ONLY: get_stored_canonical
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: int_4
  USE dbcsr_methods,                   ONLY: dbcsr_distribution_local_cols,&
                                             dbcsr_distribution_local_rows,&
                                             dbcsr_distribution_ncols,&
                                             dbcsr_distribution_nlocal_cols,&
                                             dbcsr_distribution_nlocal_rows,&
                                             dbcsr_distribution_nrows,&
                                             dbcsr_wm_use_mutable
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size
  USE dbcsr_types,                     ONLY: &
       dbcsr_distribution_obj, dbcsr_meta_size, dbcsr_num_slots, dbcsr_obj, &
       dbcsr_slot_blk_p, dbcsr_slot_col_i, dbcsr_slot_nblkcols_local, &
       dbcsr_slot_nblkcols_total, dbcsr_slot_nblkrows_local, &
       dbcsr_slot_nblkrows_total, dbcsr_slot_nblks, &
       dbcsr_slot_nfullcols_local, dbcsr_slot_nze, dbcsr_slot_row_p, &
       dbcsr_slot_size, dbcsr_type
  USE dbcsr_util,                      ONLY: dbcsr_set_debug,&
                                             swap

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_index_operations'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034


  ! Index transformations
  PUBLIC :: dbcsr_make_index_canonical, make_index_triangular,&
            transpose_index_local
  ! Dense/Sparse
  PUBLIC :: make_dense_index, make_undense_index
  ! Working with DBCSR and linear indices
  PUBLIC :: dbcsr_make_dbcsr_index, dbcsr_sort_indices,&
            dbcsr_sort_many_indices, merge_index_arrays
  ! Index array manipulation
  PUBLIC :: dbcsr_addto_index_array, dbcsr_clearfrom_index_array,&
            dbcsr_repoint_index


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

    INTEGER                                  :: blk, col, nblks, row, &
                                                stored_col, stored_row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: row_i
    LOGICAL                                  :: dbg = .FALSE., tr

!   ---------------------------------------------------------------------------

    nblks = SIZE(old_blk_p)
    ALLOCATE(row_i (nblks))
    IF (dbg) THEN
       WRITE(*,*)"old row_p", old_row_p
       WRITE(*,*)"old col_i", old_col_i
       WRITE(*,*)"old blk_p", old_blk_p
    ENDIF
    DO row = 1, SIZE(old_row_p)-1
       DO blk = old_row_p(row)+1, old_row_p(row+1)
          col = old_col_i(blk)
          stored_row = row
          stored_col = col
          tr = .FALSE.
          CALL get_stored_canonical (matrix, stored_row, stored_col, tr)
          IF (dbg) &
               WRITE(*,'(A,2(1X,I5),A,2(1X,I5),";",I7,1X,L1)')&
               routineN//" X->",row,col,"->",&
               stored_row, stored_col,blk,tr
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
    IF (dbg) THEN
       WRITE(*,*)"new row_p", new_row_p
       WRITE(*,*)"new row_i", row_i
       WRITE(*,*)"new col_i", new_col_i
       WRITE(*,*)"new blk_p", new_blk_p
    ENDIF
  END SUBROUTINE make_index_canonical


! *****************************************************************************
!> \brief Makes a CP2K triangular index given the index arrays
!>
!> \note
!> This routine uses hard-coded logic as to what constitutes a
!> canonical ordering
!> \par Description of canonical ordering
!> A non-(anti)symmetric matrix is left as is. Otherwise, the row and column
!> are stored in the position prescribed by the distribution.
! *****************************************************************************
  SUBROUTINE make_index_triangular (new_row_p, new_col_i, new_blk_p,&
       old_row_p, old_col_i, old_blk_p, matrix)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: new_row_p, new_col_i, &
                                                new_blk_p
    INTEGER, DIMENSION(:), INTENT(IN)        :: old_row_p, old_col_i, &
                                                old_blk_p
    TYPE(dbcsr_type), INTENT(IN)             :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'make_index_triangular', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, nblks, row, &
                                                stored_col, stored_row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: row_i
    LOGICAL                                  :: dbg = .FALSE., tr

!   ---------------------------------------------------------------------------

    nblks = SIZE(old_blk_p)
    ALLOCATE(row_i (nblks))
    IF (dbg) THEN
       WRITE(*,*)"old row_p", old_row_p
       WRITE(*,*)"old col_i", old_col_i
       WRITE(*,*)"old blk_p", old_blk_p
    ENDIF
    DO row = 1, SIZE(old_row_p)-1
       DO blk = old_row_p(row)+1, old_row_p(row+1)
          col = old_col_i(blk)
          stored_row = row
          stored_col = col
          tr = .FALSE.
          CALL get_stored_canonical (matrix, stored_row, stored_col, tr)
          IF (stored_row .GT. stored_col) THEN
             CALL swap (stored_row, stored_col)
             tr = .NOT. tr
          ENDIF
          IF (dbg) &
               WRITE(*,'(A,2(1X,I5),A,2(1X,I5),";",I7,1X,L1)')&
               routineN//" X->",row,col,"->",&
               stored_row, stored_col,blk,tr
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
    IF (dbg) THEN
       WRITE(*,*)"new row_p", new_row_p
       WRITE(*,*)"new row_i", row_i
       WRITE(*,*)"new col_i", new_col_i
       WRITE(*,*)"new blk_p", new_blk_p
    ENDIF
  END SUBROUTINE make_index_triangular


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
!> \brief Adds data to the index. Incleases the index size when neccessary.
!> \param[in,out] matrix      bcsr matrix
!> \param[in] slot  which index array to add (e.g., dbcsr_slot_row_blk_sizes)
!> \param[in] data  (optional) array holding the index data to add to the index
!>                  array
!> \param[in] reservation     (optional) only reserve space for subsequent
!>                            array
!> \param[in] extra (optional) reserve extra space for later additions
! *****************************************************************************
  SUBROUTINE dbcsr_addto_index_array(matrix, slot, DATA, reservation, extra, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, INTENT(IN)                      :: slot
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: DATA
    INTEGER, INTENT(IN), OPTIONAL            :: reservation, extra
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_addto_index_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: deplus, error_handler, space, &
                                                ub, ub_new
    LOGICAL                                  :: dbg

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    dbg = .FALSE.
    CALL dbcsr_assert (ASSOCIATED (matrix%index), dbcsr_fatal_level,&
         dbcsr_internal_error, routineP,&
         "Index must be preallocated.",__LINE__,error)
    CALL dbcsr_assert (UBOUND(matrix%index,1),'GE',dbcsr_num_slots,&
         dbcsr_failure_level, dbcsr_internal_error, routineP,&
         "Actual index size less than declared size",__LINE__,error)
    CALL dbcsr_assert(PRESENT(DATA), 'OR', PRESENT(reservation), dbcsr_fatal_level,&
         dbcsr_caller_error, routineP,&
         'Either an array or its size must be specified.',__LINE__,error)
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
       IF(dbg) WRITE(*,*)routineP//' Shortening index'
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
            nocopy=.FALSE., special=matrix%special_memory,error=error)
       CALL dbcsr_repoint_index(matrix)
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Adding slot',slot,'at',&
         matrix%index(slot),'sized',space
    CALL dbcsr_repoint_index(matrix, slot)
    IF (PRESENT(DATA)) &
         matrix%index(matrix%index(slot):matrix%index(slot+1)) = DATA(:)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_addto_index_array

! *****************************************************************************
!> \brief Removes data from the index.
!> \param[in,out] matrix      bcsr matrix
!> \param[in] slot  which index array to remove (e.g., dbcsr_slot_row_blk_sizes)
! *****************************************************************************
  SUBROUTINE dbcsr_clearfrom_index_array(matrix, slot)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, INTENT(IN)                      :: slot

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_clearfrom_index_array', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: space
    INTEGER, DIMENSION(3)                    :: max_extents
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (matrix%index), dbcsr_fatal_level,&
         dbcsr_internal_error, routineP,&
         "Index must be preallocated.",__LINE__,error)
    CALL dbcsr_assert (UBOUND(matrix%index,1),'GE',dbcsr_num_slots,&
         dbcsr_failure_level, dbcsr_internal_error, routineP,&
         "Actual index size less than declared size",__LINE__,error)
    IF(dbg) WRITE(*,*)routineP//' index',&
         matrix%index(:dbcsr_num_slots)
    ! Clear index entry pointer
    matrix%index(slot) = 0
    matrix%index(slot+1) = 0
    ! Update the declared index size
    max_extents = (/ &
         matrix%index(dbcsr_slot_row_p+1),&
         matrix%index(dbcsr_slot_col_i+1),&
         matrix%index(dbcsr_slot_blk_p+1) /)
    space = MAX (MAXVAL (max_extents), dbcsr_num_slots)
    matrix%index (dbcsr_slot_size) = space
  END SUBROUTINE dbcsr_clearfrom_index_array



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

    INTEGER                                  :: i, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: buf, buf_d
    TYPE(dbcsr_error_type)                   :: error

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    IF (SIZE (row_i) .EQ. 0) RETURN
    CALL dbcsr_assert(SIZE(row_i),'GE',n, dbcsr_failure_level,&
         dbcsr_caller_error, routineP, 'row_i too small',__LINE__,error)
    CALL dbcsr_assert(SIZE(col_i),'GE',n, dbcsr_failure_level,&
         dbcsr_caller_error, routineP, 'col_i too small',__LINE__,error)
    IF (PRESENT (blk_p)) CALL dbcsr_assert(SIZE(blk_p),'GE',n, dbcsr_failure_level,&
         dbcsr_caller_error, routineP, 'blk_p too small',__LINE__,error)
    IF (PRESENT (blk_p)) THEN
       ALLOCATE(buf(n), stat=stat)
       CALL dbcsr_assert(stat == 0, dbcsr_failure_level,&
            dbcsr_caller_error, routineP, 'buf',__LINE__,error)
       buf(1:n) = blk_p(1:n)
    ENDIF
    IF (PRESENT (blk_d)) THEN
       ALLOCATE(buf_d(n), stat=stat)
       CALL dbcsr_assert(stat == 0, dbcsr_failure_level,&
            dbcsr_caller_error, routineP, 'buf_d',__LINE__,error)
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
  END SUBROUTINE dbcsr_sort_indices

! *****************************************************************************
!> \brief Sorts all work matrices.
!> \param[in,out] matrix      sort work matrices in this matrix
! *****************************************************************************
  SUBROUTINE dbcsr_sort_many_indices(matrix)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sort_many_indices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    DO i = 1, SIZE (matrix%wms)
       IF (.NOT. dbcsr_wm_use_mutable (matrix%wms(i))) THEN
          IF (matrix%wms(i)%lastblk .GT. 0) THEN
             CALL dbcsr_sort_indices (matrix%wms(i)%lastblk,&
                  matrix%wms(i)%row_i,&
                  matrix%wms(i)%col_i, matrix%wms(i)%blk_p)
          ENDIF
       ENDIF
    ENDDO
  END SUBROUTINE dbcsr_sort_many_indices



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
!> \param[in] cp2k            (optional) make CP2K triangular index
!>                            from canonical; default is false
! *****************************************************************************
  SUBROUTINE dbcsr_make_index_canonical (matrix, cp2k)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: cp2k

    INTEGER                                  :: nb, nc, nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: new_blk_p, new_col_i, &
                                                new_row_p
    LOGICAL                                  :: rev

!   ---------------------------------------------------------------------------

    rev = .FALSE.
    IF (PRESENT (cp2k)) rev = cp2k
    nr = SIZE(matrix%m%row_p)
    ALLOCATE (new_row_p (nr))
    nc = SIZE(matrix%m%col_i)
    ALLOCATE (new_col_i (nc))
    nb = SIZE(matrix%m%blk_p)
    ALLOCATE (new_blk_p (nb))
    IF (rev) THEN
       CALL make_index_triangular (new_row_p, new_col_i, new_blk_p,&
            matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p, matrix%m)
    ELSE
       CALL make_index_canonical (new_row_p, new_col_i, new_blk_p,&
            matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p, matrix%m)
    ENDIF
    matrix%m%row_p(:) = new_row_p
    matrix%m%col_i(:) = new_col_i
    matrix%m%blk_p(:) = new_blk_p
  END SUBROUTINE dbcsr_make_index_canonical


! *****************************************************************************
!> \brief Makes the index for a dense matrix
!> \param[out] row_p, col_i, blk_p            Storage for new index
!> \param[in,out] distribution                Distribution for matrix
!> \param[in,out] matrix                      Matrix
!> \param[in,out] meta                        Metadata updates for new index
!> \param[in] nblkrows_total                  Total blocked rows
!> \param[in] nblkcols_total                  Total blocked columns
!> \param[in] myblkrows                       List of blocked rows in my
!>                                            process row
!> \param[in] myblkcols                       List of blocked columns in my
!>                                            process column
!> \param[in] make_tr                         Dense blocks are transposed
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE make_dense_index (row_p, col_i, blk_p,&
       nblkrows_total, nblkcols_total, myblkrows, myblkcols,&
       row_blk_offsets, col_blk_offsets, meta, make_tr, error)

    !INTEGER, DIMENSION(:), INTENT(OUT)       :: row_p, col_i, blk_p
    INTEGER, INTENT(IN)                      :: nblkrows_total
    INTEGER, DIMENSION(:), INTENT(OUT)       :: blk_p, col_i
    INTEGER, DIMENSION(1:nblkrows_total+1), &
      INTENT(OUT)                            :: row_p
    INTEGER, INTENT(IN)                      :: nblkcols_total
    INTEGER, DIMENSION(:), INTENT(IN)        :: myblkrows, myblkcols, &
                                                row_blk_offsets, &
                                                col_blk_offsets
    INTEGER, DIMENSION(dbcsr_meta_size), &
      INTENT(INOUT)                          :: meta
    LOGICAL, INTENT(IN), OPTIONAL            :: make_tr
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_dense_index', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: blk, c, col_l, mynblkcols, &
                                                mynblkrows, nblks, nze, &
                                                prev_row, row, row_l, &
                                                sign_carrier, sz

!   ---------------------------------------------------------------------------

    sign_carrier = 1
    IF (PRESENT (make_tr)) THEN
       IF (make_tr) sign_carrier = -1
    ENDIF
    mynblkrows = SIZE(myblkrows)
    mynblkcols = SIZE(myblkcols)
    meta(dbcsr_slot_nblkrows_local) = mynblkrows
    meta(dbcsr_slot_nblkcols_local) = mynblkcols
    nblks = mynblkrows * mynblkcols
    nze = 1
    IF (nblks .EQ. 0) THEN
       row_p(1:) = 0
    ELSE
       row_p(1) = 0
       !row_p(nrows+1) = nblks
       prev_row = 1
       blk = 0
       DO row_l = 1, mynblkrows
          row = myblkrows(row_l)
          row_p(prev_row+1:row) = blk
          DO col_l = 1, mynblkcols
             c = myblkcols(col_l)
             col_i(blk+col_l) = c
             sz = (row_blk_offsets(row+1)-row_blk_offsets(row))*&
                  (col_blk_offsets(c+1)-col_blk_offsets(c))
             IF (sz .GT. 0) THEN
                blk_p(blk+col_l) = SIGN (nze, sign_carrier)
                nze = nze + sz
             ELSE
                blk_p(blk+col_l) = 0
             ENDIF
          ENDDO
          prev_row = row
          blk = blk + mynblkcols
       END DO
       CALL dbcsr_assert (blk, "EQ", nblks, dbcsr_fatal_level,&
            dbcsr_internal_error, routineN, "Block mismatch", __LINE__,&
            error=error)
       row_p(prev_row+1:nblkrows_total+1) = nblks
    ENDIF
    IF (dbg) THEN
       WRITE(*,*)routineN//" new index"
       WRITE(*,*)"row_p=",row_p
       WRITE(*,*)"col_i=",col_i
       WRITE(*,*)"blk_p=",blk_p
    ENDIF
    meta(dbcsr_slot_nblkrows_total) = nblkrows_total
    meta(dbcsr_slot_nblkcols_total) = nblkcols_total
  END SUBROUTINE make_dense_index

! *****************************************************************************
!> \brief Makes a blocked index from a dense matrix
!> \param[out] row_p, col_i, blk_p            Storage for new index
!> \param[in] distribution     Blocked distribution
!> \param[in] row_blk_offsets  Row block offsets
!> \param[in] col_blk_offsets  Column block offsets
!> \param[in,out] meta         Metadata updates for new index
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE make_undense_index (&
       row_p, col_i, blk_p,&
       distribution,  local_row_offsets, local_col_offsets,&
       meta)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: row_p, col_i, blk_p
    TYPE(dbcsr_distribution_obj)             :: distribution
    INTEGER, DIMENSION(:), INTENT(IN)        :: local_row_offsets, &
                                                local_col_offsets
    INTEGER, DIMENSION(dbcsr_meta_size), &
      INTENT(INOUT)                          :: meta

    CHARACTER(len=*), PARAMETER :: routineN = 'make_undense_index', &
      routineP = moduleN//':'//routineN

    INTEGER :: col, lr, lrow, nblkcols_local, nblkrows_local, nblkrows_total, &
      nfullcols_local, prev_row, row
    INTEGER, DIMENSION(:), POINTER           :: local_cols, local_rows

!   ---------------------------------------------------------------------------

    local_cols => array_data (dbcsr_distribution_local_cols (distribution))
    local_rows => array_data (dbcsr_distribution_local_rows (distribution))
    meta(dbcsr_slot_nblkrows_total) = dbcsr_distribution_nrows (distribution)
    meta(dbcsr_slot_nblkcols_total) = dbcsr_distribution_ncols (distribution)
    meta(dbcsr_slot_nblkrows_local) = dbcsr_distribution_nlocal_rows (distribution)
    meta(dbcsr_slot_nblkcols_local) = dbcsr_distribution_nlocal_cols (distribution)
    nblkrows_total = meta(dbcsr_slot_nblkrows_total)
    nblkcols_local = meta(dbcsr_slot_nblkcols_local)
    nblkrows_local = meta(dbcsr_slot_nblkrows_local)
    nfullcols_local = meta(dbcsr_slot_nfullcols_local)
    ! Fill the row_p array.
    lr = 0
    row_p(1) = 0
    prev_row = 1
    DO lrow = 1, nblkrows_local
       row = local_rows(lrow)
       row_p(prev_row+1:row) = lr
       lr = lr + nblkcols_local
       row_p(row+1) = lr
       prev_row = row
    ENDDO
    row_p(prev_row+1:nblkrows_total+1) = lr
    !
    FORALL (row = 1 : nblkrows_local)
       FORALL (col = 1 : nblkcols_local)
          col_i(nblkcols_local*(row-1)+col) = local_cols(col)
          blk_p(nblkcols_local*(row-1)+col) = 1 + &
               (local_row_offsets(row)-1)*nfullcols_local&
               + (local_col_offsets(col)-1)*&
               (local_row_offsets(row+1)-local_row_offsets(row))
       END FORALL
    END FORALL
  END SUBROUTINE make_undense_index



! *****************************************************************************
!> \brief Merges two indices
!> \note Used in local multiply
!>
!> Assumes they are both pre-sorted
!> \param[in,out] new_row_i, new_col_i, new_blk_p      merged result
!> \param[in] new_size                  size of merged index
!> \param[in,out] old_row_i, old_col_i, old_blk_p      current index
!> \param[in] old_size                  size of current index
!> \param[in] add_ip                    index to add into the current index
!> \param[in] add_size                  size of index to add into the current
!>                                      index
! *****************************************************************************
  SUBROUTINE merge_index_arrays (new_row_i, new_col_i, new_blk_p, new_size,&
       old_row_i, old_col_i, old_blk_p, old_size,&
       add_ip, add_size, new_blk_d, old_blk_d, error)
    INTEGER, INTENT(IN)                      :: new_size
    INTEGER, DIMENSION(new_size), &
      INTENT(OUT)                            :: new_blk_p, new_col_i, &
                                                new_row_i
    INTEGER, INTENT(IN)                      :: old_size
    INTEGER, DIMENSION(old_size), INTENT(IN) :: old_blk_p, old_col_i, &
                                                old_row_i
    INTEGER, INTENT(IN)                      :: add_size
    INTEGER, DIMENSION(3, add_size), &
      INTENT(IN)                             :: add_ip
    INTEGER, DIMENSION(new_size), &
      INTENT(OUT), OPTIONAL                  :: new_blk_d
    INTEGER, DIMENSION(old_size), &
      INTENT(IN), OPTIONAL                   :: old_blk_d
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    INTEGER                                  :: add_blk, i, merge_from_whom, &
                                                new_blk, old_blk
    LOGICAL                                  :: multidata

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    multidata = PRESENT (old_blk_d) .AND. PRESENT (new_blk_d)
    CALL dbcsr_assert(old_size+add_size .EQ. new_size, &
         dbcsr_warning_level, dbcsr_internal_error ,"merge_arrays",&
         "Mismatch of new and old size",__LINE__,error)
    IF (dbg) THEN
       WRITE (*,*) " Old array", old_size
       DO i = 1, old_size
          WRITE(*,'(I7,2X,I7,2X,I7)')old_row_i(i),old_col_i(i),old_blk_p(i)
       ENDDO
       WRITE (*,*) " Add array", add_size
       DO i = 1, add_size
          WRITE(*,'(I7,2X,I7,2X,I7)')add_ip (1:3, i)
       ENDDO
    ENDIF
    IF (add_size .GT. 0) THEN
       old_blk = 1
       add_blk = 1
       new_blk = 1
       IF (old_size .EQ. 0) THEN
          new_row_i(1:add_size) = add_ip(1, 1:add_size)
          new_col_i(1:add_size) = add_ip(2, 1:add_size)
          new_blk_p(1:add_size) = add_ip(3, 1:add_size)
          !IF (multidata) new_blk_d(1:add_size) = add_ip(4, 1:add_size)
       ELSE
          DO WHILE (new_blk .LE. new_size)
             merge_from_whom = 0
             IF (old_blk .LE. old_size .AND. add_blk .LE. add_size) THEN
                IF (add_ip(1, add_blk) .EQ. old_row_i(old_blk)&
                     .AND.add_ip(2, add_blk) .EQ. old_col_i(old_blk)) THEN
                   !WRITE(*,*)"Duplicate block! addblk", add_blk, "oldblk", old_blk
                ENDIF
                ! Rows come first
                IF (add_ip(1, add_blk) .LT. old_row_i(old_blk)) THEN
                   merge_from_whom = 2
                ELSEIF (add_ip(1, add_blk) .GT. old_row_i(old_blk)) THEN
                   merge_from_whom = 1
                ELSE ! Same rows, so now come the columns
                   IF (add_ip(2, add_blk) .LE. old_col_i(old_blk)) THEN
                      ! Merges from the add array
                      merge_from_whom = 2
                   ELSE
                      ! Merges from the old array
                      merge_from_whom = 1
                   ENDIF
                ENDIF
             ELSE
                IF (add_blk .LE. add_size) THEN
                   ! Merges from the add array
                   merge_from_whom = 2
                ELSEIF (old_blk .LE. old_size) THEN
                   ! Merges from the old array
                   merge_from_whom = 1
                ELSE
                   ! Hmmm, nothing to merge...
                   !WRITE(*,*)"Error: Ran out of data to merge"
                ENDIF
             ENDIF
             SELECT CASE (merge_from_whom)
             CASE (2)
                ! Merges from the add array
                new_row_i(new_blk) = add_ip(1, add_blk)
                new_col_i(new_blk) = add_ip(2, add_blk)
                new_blk_p(new_blk) = add_ip(3, add_blk)
                !IF (multidata) new_blk_d(new_blk) = add_ip(4, add_blk)
                add_blk = add_blk + 1
             CASE (1)
                ! Merges from the old array
                new_row_i(new_blk) = old_row_i(old_blk)
                new_col_i(new_blk) = old_col_i(old_blk)
                new_blk_p(new_blk) = old_blk_p(old_blk)
                IF (multidata) new_blk_p(new_blk) = old_blk_d(old_blk)
                old_blk = old_blk + 1
             CASE DEFAULT
                !WRITE(*,*)"Error: Nothing to merge"
             END SELECT
             new_blk = new_blk + 1
          ENDDO
       ENDIF
    ELSE
       new_row_i(1:old_size) = old_row_i(1:old_size)
       new_col_i(1:old_size) = old_col_i(1:old_size)
       new_blk_p(1:old_size) = old_blk_p(1:old_size)
       IF (multidata) new_blk_d(1:old_size) = old_blk_d(1:old_size)
    ENDIF
    IF (dbg) THEN
       WRITE (*,*) " New array"
       DO i = 1, new_size
          WRITE(*,'(4(2X,I7))')new_row_i(i),new_col_i(i),new_blk_p(i)
       ENDDO
    ENDIF
  END SUBROUTINE merge_index_arrays

  SUBROUTINE sort ( arr, n, index )
    INTEGER, INTENT(IN)                      :: n
    INTEGER(KIND=int_4), INTENT(INOUT)       :: arr(1:n)
    INTEGER, INTENT(OUT)                     :: INDEX(1:n)

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sort', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: m = 7, nstack = 50

    INTEGER                                  :: i, ib, ir, istack(1:nstack), &
                                                itemp, j, jstack, k, l
    INTEGER(kind=int_4)                      :: a, temp

    IF (n==0) RETURN
    DO i = 1, n
       INDEX(i) = i
    END DO
    !
    IF (ALL(arr==arr(1))) RETURN ! Nothing to order..
    jstack = 0
    l = 1
    ir = n
1   IF (ir-l<m) THEN
       DO j = l + 1, ir
          a = arr(j)
          ib = INDEX(j)
          DO i = j - 1, 1, -1
             IF (arr(i)<=a) GO TO 2
             arr(i+1) = arr(i)
             INDEX(i+1) = INDEX(i)
          END DO
          i = 0
2         arr(i+1) = a
          INDEX(i+1) = ib
       END DO
       IF (jstack==0) RETURN
       ir = istack(jstack)
       l = istack(jstack-1)
       jstack = jstack - 2
    ELSE
       k = (l+ir)/2
       temp = arr(k)
       arr(k) = arr(l+1)
       arr(l+1) = temp
       itemp = INDEX(k)
       INDEX(k) = INDEX(l+1)
       INDEX(l+1) = itemp
       IF (arr(l+1)>arr(ir)) THEN
          temp = arr(l+1)
          arr(l+1) = arr(ir)
          arr(ir) = temp
          itemp = INDEX(l+1)
          INDEX(l+1) = INDEX(ir)
          INDEX(ir) = itemp
       END IF
       IF (arr(l)>arr(ir)) THEN
          temp = arr(l)
          arr(l) = arr(ir)
          arr(ir) = temp
          itemp = INDEX(l)
          INDEX(l) = INDEX(ir)
          INDEX(ir) = itemp
       END IF
       IF (arr(l+1)>arr(l)) THEN
          temp = arr(l+1)
          arr(l+1) = arr(l)
          arr(l) = temp
          itemp = INDEX(l+1)
          INDEX(l+1) = INDEX(l)
          INDEX(l) = itemp
       END IF
       i = l + 1
       j = ir
       a = arr(l)
       ib = INDEX(l)
3      CONTINUE
       i = i + 1
       IF (arr(i)<a) GO TO 3
4      CONTINUE
       j = j - 1
       IF (arr(j)>a) GO TO 4
       IF (j<i) GO TO 5
       temp = arr(i)
       arr(i) = arr(j)
       arr(j) = temp
       itemp = INDEX(i)
       INDEX(i) = INDEX(j)
       INDEX(j) = itemp
       GO TO 3
5      arr(l) = arr(j)
       arr(j) = a
       INDEX(l) = INDEX(j)
       INDEX(j) = ib
       jstack = jstack + 2
       IF (jstack>nstack) STOP ' Nstack too small in sortr'
       IF (ir-i+1>=j-l) THEN
          istack(jstack) = ir
          istack(jstack-1) = i
          ir = j - 1
       ELSE
          istack(jstack) = j - 1
          istack(jstack-1) = l
          l = i
       END IF
    END IF
    GO TO 1
  END SUBROUTINE sort



END MODULE dbcsr_index_operations

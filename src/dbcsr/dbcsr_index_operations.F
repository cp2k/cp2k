!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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

  USE array_types,                     ONLY: array_data,&
                                             array_size
  USE dbcsr_config,                    ONLY: default_resize_factor
  USE dbcsr_dist_operations,           ONLY: get_stored_canonical
  USE dbcsr_error_handling,            ONLY: &
       dbcsr_assert, dbcsr_caller_error, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_failure_level, dbcsr_fatal_level, &
       dbcsr_internal_error, dbcsr_unimplemented_error_nr, &
       dbcsr_warning_level, dbcsr_wrong_args_error
  USE dbcsr_methods,                   ONLY: &
       dbcsr_distribution, dbcsr_distribution_local_cols, &
       dbcsr_distribution_local_rows, dbcsr_distribution_ncols, &
       dbcsr_distribution_nlocal_cols, dbcsr_distribution_nlocal_rows, &
       dbcsr_distribution_nrows, dbcsr_distribution_thread_dist, &
       dbcsr_get_index_memory_type, dbcsr_nblkcols_total, &
       dbcsr_nblkrows_total, dbcsr_wm_use_mutable
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size,&
                                             memory_allocate,&
                                             memory_deallocate
  USE dbcsr_toollib,                   ONLY: joaat_hash,&
                                             sort,&
                                             swap
  USE dbcsr_types,                     ONLY: &
       dbcsr_distribution_obj, dbcsr_meta_size, dbcsr_num_slots, dbcsr_obj, &
       dbcsr_slot_blk_p, dbcsr_slot_col_i, dbcsr_slot_coo_l, &
       dbcsr_slot_home_prow, dbcsr_slot_home_vprow, &
       dbcsr_slot_nblkcols_local, dbcsr_slot_nblkcols_total, &
       dbcsr_slot_nblkrows_local, dbcsr_slot_nblkrows_total, &
       dbcsr_slot_nblks, dbcsr_slot_nfullcols_local, dbcsr_slot_nze, &
       dbcsr_slot_row_p, dbcsr_slot_size, dbcsr_slot_thr_c, dbcsr_type
  USE kinds,                           ONLY: int_4,&
                                             int_8

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_index_operations'


  LOGICAL, PARAMETER :: careful_mod = .FALSE.
  LOGICAL, PARAMETER :: debug_mod = .FALSE.




  ! Index transformations
  PUBLIC :: dbcsr_make_index_canonical, make_index_triangular,&
            transpose_index_local, dbcsr_order_matrix_index
  PUBLIC :: dbcsr_make_index_local_row, dbcsr_has_local_row_index
  PUBLIC :: dbcsr_make_index_list
  ! Dense/Sparse
  PUBLIC :: make_dense_index, make_undense_index
  ! Working with DBCSR and linear indices
  PUBLIC :: dbcsr_make_dbcsr_index, dbcsr_sort_indices,&
            dbcsr_sort_many_indices, merge_index_arrays,&
            dbcsr_expand_row_index,&
            dbcsr_count_row_index, dbcsr_build_row_index,&
            dbcsr_index_prune_deleted,&
            dbcsr_index_compact
  PUBLIC :: dbcsr_index_checksum
  ! Index array manipulation
  PUBLIC :: dbcsr_addto_index_array, dbcsr_clearfrom_index_array,&
            dbcsr_repoint_index, dbcsr_make_index_exist

  INTERFACE dbcsr_count_row_index
     MODULE PROCEDURE dbcsr_count_row_index_copy,&
                      dbcsr_count_row_index_inplace
  END INTERFACE

  INTERFACE dbcsr_build_row_index
     MODULE PROCEDURE dbcsr_build_row_index_copy,&
                      dbcsr_build_row_index_inplace
  END INTERFACE

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
    LOGICAL                                  :: tr

!   ---------------------------------------------------------------------------

    nblks = SIZE(old_blk_p)
    ALLOCATE(row_i (nblks))
    IF (debug_mod) THEN
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
          IF (debug_mod) &
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
    IF (debug_mod) THEN
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
    LOGICAL                                  :: tr

!   ---------------------------------------------------------------------------

    nblks = SIZE(old_blk_p)
    ALLOCATE(row_i (nblks))
    IF (debug_mod) THEN
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
          IF (debug_mod) &
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
    IF (debug_mod) THEN
       WRITE(*,*)"new row_p", new_row_p
       WRITE(*,*)"new row_i", row_i
       WRITE(*,*)"new col_i", new_col_i
       WRITE(*,*)"new blk_p", new_blk_p
    ENDIF
  END SUBROUTINE make_index_triangular

! *****************************************************************************
!> \brief Collapses a row_p index
! *****************************************************************************
  SUBROUTINE dbcsr_make_dbcsr_index (row_p, row_i, nrows, nblks)
    INTEGER, INTENT(in)                      :: nrows, nblks
    INTEGER, DIMENSION(1:nrows+1), &
      INTENT(out)                            :: row_p
    INTEGER, DIMENSION(1:nblks), INTENT(in)  :: row_i

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_dbcsr_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, error_handle, row
    TYPE(dbcsr_error_type)                   :: error

    CALL dbcsr_error_set(routineN, error_handle, error)

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

    CALL dbcsr_error_stop(error_handle, error)

  END SUBROUTINE dbcsr_make_dbcsr_index

! *****************************************************************************
!> \brief Expands a row_p index
! *****************************************************************************
  PURE SUBROUTINE dbcsr_expand_row_index (row_p, row_i, nrows, nblks)
    INTEGER, INTENT(in)                      :: nrows, nblks
    INTEGER, DIMENSION(1:nrows+1), &
      INTENT(IN)                             :: row_p
    INTEGER, DIMENSION(1:nblks), INTENT(out) :: row_i

    INTEGER                                  :: row

    DO row = 1, nrows
       row_i(row_p(row)+1 : row_p(row+1)) = row
    ENDDO
  END SUBROUTINE dbcsr_expand_row_index


! *****************************************************************************
!> \brief Counts columns-per-row count from row index array, in-place.
!> \param[in] nrows           number of rows
!> \param[in,out] rows        the row_p index (input); the count of the number
!>                            of columns per row (output)
! *****************************************************************************
  PURE SUBROUTINE dbcsr_count_row_index_inplace (rows, nrows)
    INTEGER, INTENT(IN)                      :: nrows
    INTEGER, DIMENSION(1:nrows+1), &
      INTENT(INOUT)                          :: rows

    INTEGER                                  :: row

    DO row = 1, nrows
       rows(row) = rows(row+1) - rows(row)
    ENDDO
    rows(nrows+1) = 0
  END SUBROUTINE dbcsr_count_row_index_inplace

! *****************************************************************************
!> \brief Counts columns-per-row count from row index array.
!> \param[in] nrows           number of rows
!> \param[in] rows            the row_p index (input)
!> \param[out] counts         the count of the number of columns per row
! *****************************************************************************
  PURE SUBROUTINE dbcsr_count_row_index_copy (rows, counts, nrows)
    INTEGER, INTENT(IN)                      :: nrows
    INTEGER, DIMENSION(1:nrows), INTENT(OUT) :: counts
    INTEGER, DIMENSION(1:nrows+1), &
      INTENT(IN)                             :: rows

    INTEGER                                  :: row

    FORALL (row = 1:nrows)
       counts(row) = rows(row+1) - rows(row)
    END FORALL
  END SUBROUTINE dbcsr_count_row_index_copy

! *****************************************************************************
!> \brief Builds row index array from a columns-per-row count, in-place.
!> \param[in] nrows           number of rows
!> \param[in,out] rows        count of the number of colums per row (input);
!>                            the row_p index (output)
! *****************************************************************************
  PURE SUBROUTINE dbcsr_build_row_index_inplace (rows, nrows)
    INTEGER, INTENT(IN)                      :: nrows
    INTEGER, DIMENSION(1:nrows+1), &
      INTENT(INOUT)                          :: rows

    INTEGER                                  :: o, old_count, row

    old_count = rows(1)
    rows(1) = 0
    IF (nrows .GE. 1) THEN
       DO row = 2, nrows+1
          o = rows(row)
          rows(row) = rows(row-1) + old_count
          old_count = o
       ENDDO
    ENDIF
  END SUBROUTINE dbcsr_build_row_index_inplace

! *****************************************************************************
!> \brief Builds row index array from a columns-per-row count.
!> \param[in] nrows           number of rows
!> \param[in] counts          count of the number of colums per row
!> \param[out] rows           count of the number of colums per row (input);
!>                            the row_p index (output)
! *****************************************************************************
  PURE SUBROUTINE dbcsr_build_row_index_copy (counts, rows, nrows)
    INTEGER, INTENT(IN)                      :: nrows
    INTEGER, DIMENSION(1:nrows+1), &
      INTENT(OUT)                            :: rows
    INTEGER, DIMENSION(1:nrows), INTENT(IN)  :: counts

!WTF?!rows(1) = 0
!WTF?!IF (nrows .GE. 1) THEN
!WTF?!   DO row = 2, nrows+1
!WTF?!      rows(row) = rows(row-1) + counts(rows-1)
!WTF?!   ENDDO
!WTF?!ENDIF

    rows(1:nrows) = counts(1:nrows)
    CALL dbcsr_build_row_index_inplace (rows, nrows)
  END SUBROUTINE dbcsr_build_row_index_copy


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

    INTEGER                                  :: deplus, space, ub, ub_new

!   ---------------------------------------------------------------------------

    IF (debug_mod) THEN
        CALL dbcsr_assert (ASSOCIATED (matrix%index), dbcsr_fatal_level,&
             dbcsr_internal_error, routineP,&
             "Index must be preallocated.",__LINE__,error)
        CALL dbcsr_assert (UBOUND(matrix%index,1),'GE',dbcsr_num_slots,&
             dbcsr_failure_level, dbcsr_internal_error, routineP,&
             "Actual index size less than declared size",__LINE__,error)
        CALL dbcsr_assert(PRESENT(DATA), 'OR', PRESENT(reservation), dbcsr_fatal_level,&
             dbcsr_caller_error, routineP,&
             'Either an array or its size must be specified.',__LINE__,error)
        WRITE(*,*)routineP//' index', matrix%index(:dbcsr_num_slots)
    ENDIF
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
       IF(debug_mod) WRITE(*,*)routineP//' Slot',slot,'not filled, adding at',&
            matrix%index(dbcsr_slot_size)+1,'sized',space
       matrix%index(slot) = matrix%index(dbcsr_slot_size)+1
       matrix%index(slot+1) = matrix%index(slot) + space - 1
       matrix%index(dbcsr_slot_size) = matrix%index(slot+1)
    ENDIF
    ! Shorten an index entry.
    IF (space .LT. matrix%index(slot+1) - matrix%index(slot)+1) THEN
       IF(debug_mod) WRITE(*,*)routineP//' Shortening index'
       matrix%index(slot+1) = matrix%index(slot) + space -1
       CALL dbcsr_repoint_index(matrix, slot)
    ENDIF
    ub_new = matrix%index(slot+1) + deplus
    IF(debug_mod) WRITE(*,*)routineP//' need',space,'at',matrix%index(slot),&
         'to',matrix%index(slot+1),'(',ub_new,')','have',ub
    IF (ub_new .GT. ub) THEN
       IF(debug_mod) WRITE(*,*)routineP//' Reallocating index to ubound', ub_new
       !CALL reallocate(matrix%index, 1, ub_new)
       CALL ensure_array_size(matrix%index, lb=1, ub=ub_new,&
            factor=default_resize_factor,&
            nocopy=.FALSE.,&
            memory_type=matrix%index_memory_type, error=error)
       CALL dbcsr_repoint_index(matrix)
    ENDIF
    IF(debug_mod) WRITE(*,*)routineP//' Adding slot',slot,'at',&
         matrix%index(slot),'sized',space
    CALL dbcsr_repoint_index(matrix, slot)
    IF (PRESENT(DATA)) &
         matrix%index(matrix%index(slot):matrix%index(slot+1)) = DATA(:)
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

    INTEGER                                  :: space
    INTEGER, DIMENSION(5)                    :: max_extents
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (matrix%index), dbcsr_fatal_level,&
         dbcsr_internal_error, routineP,&
         "Index must be preallocated.",__LINE__,error)
    CALL dbcsr_assert (UBOUND(matrix%index,1),'GE',dbcsr_num_slots,&
         dbcsr_failure_level, dbcsr_internal_error, routineP,&
         "Actual index size less than declared size",__LINE__,error)
    IF(debug_mod) WRITE(*,*)routineP//' index',&
         matrix%index(:dbcsr_num_slots)
    ! Clear index entry pointer
    matrix%index(slot) = 1
    matrix%index(slot+1) = 0
    CALL dbcsr_repoint_index(matrix, slot)
    ! Update the declared index size
    max_extents = (/ &
         matrix%index(dbcsr_slot_row_p+1),&
         matrix%index(dbcsr_slot_col_i+1),&
         matrix%index(dbcsr_slot_blk_p+1),&
         matrix%index(dbcsr_slot_thr_c+1),&
         matrix%index(dbcsr_slot_coo_l+1) /)
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

!   ---------------------------------------------------------------------------

    IF (m%nblks .NE. m%index(dbcsr_slot_nblks)) THEN
       m%nblks = m%index(dbcsr_slot_nblks)
       m%nze = m%index(dbcsr_slot_nze)
    ENDIF
    all = .TRUE.
    IF (PRESENT (slot)) THEN
       all = .FALSE.
       s = slot
    ELSE
       s = 0
    ENDIF

    IF (m%index(dbcsr_slot_row_p).GT.0.AND.all .OR.&
         s.EQ.dbcsr_slot_row_p) THEN
       IF (m%index(dbcsr_slot_row_p).GT.0) THEN
          m%row_p => m%index(m%index(dbcsr_slot_row_p):&
               &             m%index(dbcsr_slot_row_p+1))
       ELSE
          NULLIFY (m%row_p)
       ENDIF
    ENDIF
    IF (m%index(dbcsr_slot_col_i).GT.0.AND.all .OR.&
         s.EQ.dbcsr_slot_col_i) THEN
       IF (m%index(dbcsr_slot_col_i).GT.0) THEN
          m%col_i => m%index(m%index(dbcsr_slot_col_i):&
               &             m%index(dbcsr_slot_col_i+1))
       ELSE
          NULLIFY (m%col_i)
       ENDIF
    ENDIF
    IF (m%index(dbcsr_slot_blk_p).GT.0.AND.all .OR.&
         s.EQ.dbcsr_slot_blk_p) THEN
       IF (m%index(dbcsr_slot_blk_p).GT.0) THEN
          m%blk_p => m%index(m%index(dbcsr_slot_blk_p):&
               &             m%index(dbcsr_slot_blk_p+1))
       ELSE
          NULLIFY (m%blk_p)
       ENDIF
    ENDIF
    IF (m%index(dbcsr_slot_thr_c).GT.0.AND.all .OR.&
         s.EQ.dbcsr_slot_thr_c) THEN
       IF (m%index(dbcsr_slot_thr_c).GT.0) THEN
          m%thr_c => m%index(m%index(dbcsr_slot_thr_c):&
               &             m%index(dbcsr_slot_thr_c+1))
       ELSE
          NULLIFY (m%thr_c)
       ENDIF
    ENDIF
    IF (m%index(dbcsr_slot_coo_l).GT.0.AND.all .OR.&
         s.EQ.dbcsr_slot_coo_l) THEN
       IF (m%index(dbcsr_slot_coo_l).GT.0) THEN
          m%coo_l => m%index(m%index(dbcsr_slot_coo_l):&
               &             m%index(dbcsr_slot_coo_l+1))
       ELSE
          NULLIFY (m%coo_l)
       ENDIF
    ENDIF
    IF (all) THEN
       m%index(dbcsr_slot_nblks) = m%nblks
       m%index(dbcsr_slot_nze) = m%nze
    ENDIF
  END SUBROUTINE dbcsr_repoint_index


! *****************************************************************************
!> \brief Creates an index for the matrix.

!> \param[in,out] m           Create index for this matrix
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE dbcsr_make_index_exist(m, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: m
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_index_exist', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    CALL dbcsr_assert (ASSOCIATED (m%index),&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Index array does not yet exist.", __LINE__, error=error)
    IF (.NOT. ASSOCIATED (m%row_p)) THEN
       CALL dbcsr_addto_index_array (m, dbcsr_slot_row_p,&
            reservation=m%nblkrows_total+1, error=error)
       m%row_p(:) = 0
    ENDIF
    IF (.NOT. ASSOCIATED (m%col_i)) THEN
       CALL dbcsr_addto_index_array (m, dbcsr_slot_col_i,&
            reservation=0, error=error)
    ENDIF
    IF (.NOT. ASSOCIATED (m%blk_p)) THEN
       CALL dbcsr_addto_index_array (m, dbcsr_slot_blk_p,&
            reservation=0, error=error)
    ENDIF
    CALL dbcsr_repoint_index (m)
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_make_index_exist



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
    INTEGER, DIMENSION(1:), INTENT(INOUT)    :: row_i, col_i
    INTEGER, DIMENSION(1:), INTENT(INOUT), &
      OPTIONAL                               :: blk_p, blk_d

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sort_indices', &
      routineP = moduleN//':'//routineN
    INTEGER(KIND=int_8), PARAMETER           :: lmask8 = 4294967295_int_8

    INTEGER                                  :: error_handle, i, stat
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: sort_keys
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: buf, buf_d
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (SIZE (row_i) .EQ. 0) RETURN

    CALL dbcsr_error_set(routineN, error_handle, error)

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
    !> Create an ordering for both rows and columns. If the blk_p must
    !> be rearranged, then the col_i array will be used as a
    !> permutation vector.
    ALLOCATE (sort_keys(n), stat=stat)
    CALL dbcsr_assert (stat, "EQ", 0, dbcsr_fatal_level, dbcsr_internal_error,&
         routineN, "Could not allocate memory for sorting buffer.", __LINE__,&
         error=error)
    sort_keys(:) = IOR (ISHFT(INT(row_i(1:n),int_8), 32), INT(col_i(1:n),int_8))
    IF (PRESENT (blk_p)) col_i(1:n) = (/ (i, i=1,n) /)
    !> Now do a nice quicksort.
    CALL sort(sort_keys, n, col_i)
    ! Since blk_d is usually not present we can have two loops that
    ! are essentially the same.
    IF (PRESENT (blk_p)) THEN
       FORALL (i = 1:n)
          blk_p(i) = buf(col_i(i))
       END FORALL
       DEALLOCATE (buf)
    END IF
    IF (PRESENT (blk_d)) THEN
       FORALL (i = 1:n)
          blk_d(i) = buf_d(col_i(i))
       END FORALL
       DEALLOCATE (buf_d)
    ENDIF
    FORALL (i = 1:n)
       col_i(i) = INT (IAND(sort_keys(i), lmask8), int_4)
       row_i(i) = INT (ISHFT(sort_keys(i), -32), int_4)
    END FORALL
    DEALLOCATE (sort_keys)
    IF(debug_mod.AND.PRESENT(blk_p))&
         WRITE(*,*)routineP//' sort, blk_p =',blk_p
    CALL dbcsr_error_stop (error_handle, error)

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
!> \brief Sorts the index of a DBCSR
!> \par Description
!>      Sorts the row and column indices so that the rows monotonically
!>      increase and the columns monotonically increase within each row.
!> \param[in,out] matrix      matrix for which to sort
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE dbcsr_order_matrix_index(matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_order_matrix_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, first_blk, i, &
                                                last_blk, ncols, row, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: buff, permutation
    INTEGER, DIMENSION(:), POINTER           :: blk_p, col_i, full_blk_p, &
                                                full_col_i, full_row_p

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handle, error)

    full_row_p => matrix%m%row_p
    full_col_i => matrix%m%col_i
    full_blk_p => matrix%m%blk_p
!test!$OMP PARALLEL DEFAULT (none) &
!test!$OMP PRIVATE (permutation, buff, stat, &
!test!$             row, first_blk, last_blk, ncols, &
!test!$             col_i, blk_p, i) &
!test!$OMP SHARED (matrix, full_row_p, full_col_i, full_blk_p)
    ALLOCATE (permutation(dbcsr_nblkcols_total (matrix)), stat=stat)
    ALLOCATE (buff(dbcsr_nblkcols_total (matrix)), stat=stat)
    ! Goes through all rows and sorts the column in each row. It also
    ! reorders the blk_p correspondingly.
!test!$OMP DO
    DO row = 1, dbcsr_nblkrows_total (matrix)
       first_blk = full_row_p(row)+1
       last_blk = full_row_p(row+1)
       ncols = last_blk - first_blk + 1
       col_i => full_col_i(first_blk : last_blk)
       blk_p => full_blk_p(first_blk : last_blk)
       permutation(1:ncols) = (/ (i, i=1,ncols) /)
       buff(1:ncols) = blk_p(1:ncols)
       CALL sort(col_i, ncols, permutation)
       DO i = 1, ncols
          blk_p(i) = buff(permutation(i))
       ENDDO
    ENDDO
!test!OMP END DO NOWAIT
    DEALLOCATE (permutation, buff)
!test!$OMP END PARALLEL
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_order_matrix_index


! *****************************************************************************
!> \brief Removes the deleted blocks from the index.
!> \par Description
!> \param[in,out] matrix      Prune the index of this matrix.
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE dbcsr_index_prune_deleted(matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_index_prune_deleted', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, nblks_max, &
                                                new_blk, nrows, old_blk, row, &
                                                stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: new_blk_p, new_col_i, &
                                                new_row_count
    INTEGER, DIMENSION(:), POINTER           :: old_blk_p, old_col_i, &
                                                old_row_p

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handle, error)
    !
    old_row_p => matrix%m%row_p
    old_col_i => matrix%m%col_i
    old_blk_p => matrix%m%blk_p
    !
    nrows = matrix%m%nblkrows_total
    nblks_max = old_row_p(nrows+1)
    ALLOCATE (new_row_count(nrows), stat=stat)
    CALL dbcsr_assert (stat, "EQ", 0, dbcsr_fatal_level, dbcsr_internal_error,&
         routineN, "Can not allocate memory.", __LINE__, error=error)
    ALLOCATE (new_col_i(nblks_max), stat=stat)
    CALL dbcsr_assert (stat, "EQ", 0, dbcsr_fatal_level, dbcsr_internal_error,&
         routineN, "Can not allocate memory.", __LINE__, error=error)
    ALLOCATE (new_blk_p(nblks_max), stat=stat)
    CALL dbcsr_assert (stat, "EQ", 0, dbcsr_fatal_level, dbcsr_internal_error,&
         routineN, "Can not allocate memory.", __LINE__, error=error)
    !
    ! Build up the new index from all non-deleted blocks in the
    ! existing index.
    new_blk = 0
    DO row = 1, nrows
       new_row_count(row) = 0
       DO old_blk = old_row_p(row)+1, old_row_p(row+1)
          IF (old_blk_p(old_blk) .GT. 0) THEN
             new_blk = new_blk + 1
             new_row_count(row) = new_row_count(row) + 1
             new_col_i(new_blk) = old_col_i(old_blk)
             new_blk_p(new_blk) = old_blk_p(old_blk)
          ENDIF
       ENDDO
    ENDDO
    !
    CALL dbcsr_clearfrom_index_array(matrix%m, dbcsr_slot_row_p)
    CALL dbcsr_clearfrom_index_array(matrix%m, dbcsr_slot_col_i)
    CALL dbcsr_clearfrom_index_array(matrix%m, dbcsr_slot_blk_p)
    CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_row_p,&
         reservation=nrows+1, extra=2*new_blk, error=error)
    old_row_p => matrix%m%row_p
    CALL dbcsr_build_row_index(counts=new_row_count, rows=old_row_p,&
         nrows=nrows)
    CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_col_i, DATA=new_col_i(1:new_blk),&
         error=error)
    CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_blk_p, DATA=new_blk_p(1:new_blk),&
         error=error)
    matrix%m%nblks = new_blk
    matrix%m%index(dbcsr_slot_nblks) = new_blk
    !
    DEALLOCATE (new_col_i, new_blk_p, new_row_count)
    !
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_index_prune_deleted


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

    INTEGER                                  :: error_handle, nblks, &
                                                ncols_new, nrows_old
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: new_col_i
    LOGICAL                                  :: blks
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    blks = PRESENT (new_blk_p) .AND. PRESENT (old_blk_p)
    nblks = SIZE (old_col_i)
    nrows_old = SIZE (old_row_p)-1
    ncols_new = SIZE (new_col_p)-1
    IF (blks) new_blk_p(:) = old_blk_p(:)
    ALLOCATE (new_col_i(nblks))
    CALL dbcsr_expand_row_index(old_row_p, new_row_i, nrows_old, nblks)
    new_col_i(:) = old_col_i(:)
    CALL dbcsr_sort_indices (nblks, new_col_i, new_row_i, new_blk_p)
    CALL dbcsr_make_dbcsr_index (new_col_p, new_col_i, ncols_new, nblks)
    DEALLOCATE(new_col_i)
    CALL dbcsr_error_stop (error_handle, error)
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
    IF (debug_mod) THEN
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
!> \par Added sizes
!>      added_size_offset and added_sizes can be optionally
!>      specified. This is meant for cases where the added blocks may
!>      be duplicates of existing blocks. In this way it is possible
!>      to recalculate new block pointers to avoid wasted space.
!> \param[in,out] new_row_i, new_col_i, new_blk_p      merged result
!> \param[in] new_size                  size of merged index
!> \param[in,out] old_row_i, old_col_i, old_blk_p      current index
!> \param[in] old_size                  size of current index
!> \param[in] add_ip                    index to add into the current index
!> \param[in] add_size                  size of index to add into the current
!>                                      index
!> \param[in] added_size_offset         (optional) specify base of added sizes
!> \param[in] added_sizes               (optional) specify sizes of added
!>                                      blocks
!> \param[out] added_size               (optional) counts number of sizes of
!>                                      added blocks
!> \param[out] added_nblks              (optional) actual number of new
!>                                      elements
!> \param[in,out] error                 error
! *****************************************************************************
  SUBROUTINE merge_index_arrays (new_row_i, new_col_i, new_blk_p, new_size,&
       old_row_i, old_col_i, old_blk_p, old_size,&
       add_ip, add_size, new_blk_d, old_blk_d,&
       added_size_offset, added_sizes, added_size, added_nblks, error)
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
    INTEGER, INTENT(IN), OPTIONAL            :: added_size_offset
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: added_sizes
    INTEGER, INTENT(OUT), OPTIONAL           :: added_size, added_nblks
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'merge_index_arrays', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: add_blk, bp, i, &
                                                merge_from_whom, new_blk, &
                                                old_blk
    LOGICAL                                  :: multidata

!   ---------------------------------------------------------------------------

    multidata = PRESENT (old_blk_d) .AND. PRESENT (new_blk_d)
    CALL dbcsr_assert(old_size+add_size .EQ. new_size, &
         dbcsr_warning_level, dbcsr_internal_error ,"merge_arrays",&
         "Mismatch of new and old size",__LINE__,error)
    CALL dbcsr_assert (PRESENT (added_size_offset), "EQV",&
         PRESENT (added_sizes), dbcsr_fatal_level, dbcsr_wrong_args_error,&
         routineN, "Must specify a set of arguments", __LINE__, error=error)
    CALL dbcsr_assert (PRESENT (added_sizes), "EQV", PRESENT (added_size),&
         dbcsr_fatal_level, dbcsr_wrong_args_error,&
         routineN, "Must specify a set of arguments", __LINE__, error=error)
    IF (debug_mod) THEN
       WRITE (*,*) " Old array", old_size
       DO i = 1, old_size
          WRITE(*,'(I7,2X,I7,2X,I7)')old_row_i(i),old_col_i(i),old_blk_p(i)
       ENDDO
       WRITE (*,*) " Add array", add_size
       DO i = 1, add_size
          WRITE(*,'(I7,2X,I7,2X,I7)')add_ip (1:3, i)
       ENDDO
    ENDIF
    IF (PRESENT (added_nblks)) added_nblks = 0
    IF (PRESENT (added_size)) THEN
       added_size = 0
       bp = added_size_offset
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
          IF (PRESENT (added_nblks)) added_nblks = add_size
          IF (PRESENT (added_size)) added_size = SUM (added_sizes)
       ELSE
          DO WHILE (new_blk .LE. new_size)
             merge_from_whom = 0
             IF (old_blk .LE. old_size .AND. add_blk .LE. add_size) THEN
                IF (add_ip(1, add_blk) .EQ. old_row_i(old_blk)&
                     .AND.add_ip(2, add_blk) .EQ. old_col_i(old_blk)) THEN
                   IF (debug_mod) THEN
                      WRITE(*,*)"Duplicate block! addblk",&
                           add_blk, "oldblk", old_blk
                   ENDIF
                ENDIF
                ! Rows come first
                IF (add_ip(1, add_blk) .LT. old_row_i(old_blk)) THEN
                   merge_from_whom = 2
                ELSEIF (add_ip(1, add_blk) .GT. old_row_i(old_blk)) THEN
                   merge_from_whom = 1
                ELSE ! Same rows, so now come the columns
                   IF (add_ip(2, add_blk) .LT. old_col_i(old_blk)) THEN
                      ! Merges from the add array
                      merge_from_whom = 2
                   ELSEIF (add_ip(2, add_blk) .GT. old_col_i(old_blk)) THEN
                      ! Merges from the old array
                      merge_from_whom = 1
                   ELSE
                      ! Merge from old array and skip one in the new array
                      IF (debug_mod) THEN
                         WRITE(*,*)"Duplicate, keeping old",&
                         add_ip(1, add_blk), add_ip(2, add_blk)
                      ENDIF
                      merge_from_whom = 1
                      add_blk = add_blk + 1
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
                   merge_from_whom = 0
                   !WRITE(*,*)"Ran out of data to merge"
                ENDIF
             ENDIF
             SELECT CASE (merge_from_whom)
             CASE (2)
                ! Merges from the add array
                new_row_i(new_blk) = add_ip(1, add_blk)
                new_col_i(new_blk) = add_ip(2, add_blk)
                new_blk_p(new_blk) = add_ip(3, add_blk)
                !IF (multidata) new_blk_d(new_blk) = add_ip(4, add_blk)
                IF (PRESENT (added_nblks)) added_nblks = added_nblks + 1
                IF (PRESENT (added_sizes)) THEN
                   new_blk_p(new_blk) = bp
                   bp = bp + added_sizes(add_blk)
                   added_size = added_size + added_sizes(add_blk)
                ENDIF
                add_blk = add_blk + 1
             CASE (1)
                ! Merges from the old array
                new_row_i(new_blk) = old_row_i(old_blk)
                new_col_i(new_blk) = old_col_i(old_blk)
                new_blk_p(new_blk) = old_blk_p(old_blk)
                IF (multidata) new_blk_p(new_blk) = old_blk_d(old_blk)
                old_blk = old_blk + 1
             CASE DEFAULT
                !WRITE(*,*)"Nothing to merge"
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
    IF (debug_mod) THEN
       WRITE (*,*) " New array"
       DO i = 1, new_size
          WRITE(*,'(4(2X,I7))')new_row_i(i),new_col_i(i),new_blk_p(i)
       ENDDO
    ENDIF
  END SUBROUTINE merge_index_arrays


! *****************************************************************************
!> \brief Converts BCSR global row index to local row index.
!> \param[in,out] matrix       matrix for which to make canonical index
!> \param[in] prow             Use this process row
!> \param[in,out] error        error
! *****************************************************************************
  SUBROUTINE dbcsr_make_index_local_row (matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_index_local_row', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, lrow, &
                                                nlocal_rows, ntotal_rows, &
                                                prow, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: local_row_p
    INTEGER, DIMENSION(:), POINTER           :: local_rows

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error=error)
    CALL dbcsr_assert (ASSOCIATED (matrix%m%row_p),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "The row index must be initialized.", __LINE__, error=error)
    CALL dbcsr_assert ("NOT", matrix%m%bcsc,&
         dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
         "Not support for BCSC yet.", __LINE__, error=error)
    !
    prow = matrix%m%index(dbcsr_slot_home_vprow)
    IF (prow .LT. 0) THEN
       prow = matrix%m%index(dbcsr_slot_home_prow)
    ENDIF
    nlocal_rows = matrix%m%nblkrows_local
    ALLOCATE (local_row_p(nlocal_rows+1), stat=stat)
    CALL dbcsr_assert (stat, "EQ", 0,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Could not allocate memory.", __LINE__, error=error)
    !> The existing row_p is converted from an indexing array into a
    !> counting array.  Because it is later discarded, the counting is
    !> done in-place.
    ntotal_rows = matrix%m%nblkrows_total
    CALL dbcsr_count_row_index (matrix%m%row_p, ntotal_rows)
    ! We first have to find the local rows for the given prow.
    local_rows => array_data (matrix%m%local_rows)
    CALL dbcsr_assert (SIZE(local_rows), "EQ", nlocal_rows,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Mismatch in the number of local rows.", __LINE__, error=error)
    !> The counts are mapped to local rows,
    FORALL (lrow = 1:nlocal_rows)
       local_row_p(lrow) = matrix%m%row_p(local_rows(lrow))
    END FORALL
    CALL dbcsr_assert (SUM (matrix%m%row_p(1:ntotal_rows)), "EQ",&
         SUM(local_row_p(1:nlocal_rows)),&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Inconsistent row counts. Perhaps non-local rows contain data?.",&
         __LINE__, error=error)
    !> then converted into an index.
    CALL dbcsr_build_row_index (local_row_p, nlocal_rows)
    !> The local row index replaces the global one.
    CALL dbcsr_clearfrom_index_array(matrix%m, dbcsr_slot_row_p)
    CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_row_p, DATA=local_row_p,&
         error=error)
    !> Finally the matrix is designated as having a local-based index.
    matrix%m%local_indexing = .TRUE.
    IF (careful_mod) THEN
       CALL dbcsr_assert (array_size(matrix%m%local_rows), "EQ", nlocal_rows,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Inconsistent local row counts.",&
            __LINE__, error=error)
       CALL dbcsr_assert (array_size(matrix%m%global_rows), "EQ", ntotal_rows,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Inconsistent global row counts.",&
            __LINE__, error=error)
       IF (array_size (matrix%m%global_rows) .EQ. 0) THEN
          CALL dbcsr_assert (nlocal_rows, "EQ", 0,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Invalid number of local or global rows.",&
               __LINE__, error=error)
       ENDIF
    ENDIF
    DEALLOCATE (local_row_p)
    !
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_make_index_local_row


! *****************************************************************************
!> \brief Converts BCSR index into list indices (similar to work matrices)
!> \param[in,out] matrix       matrix for which to make canonical index
!> \param[in] thread_redist    make a thread subdistribution
!> \param[in,out] error        error
! *****************************************************************************
  SUBROUTINE dbcsr_make_index_list (matrix, thread_redist, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    LOGICAL, INTENT(IN)                      :: thread_redist
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_index_list', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: blk, error_handle, ithread, &
                                                my_cnt, nblks, nrows, nthreads
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: blki
    INTEGER, DIMENSION(0)                    :: zero_len_array
    INTEGER, DIMENSION(:), POINTER           :: global_cols, local_rows, td, &
                                                thr_c

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error=error)
    CALL dbcsr_assert (ASSOCIATED (matrix%m%row_p),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "The row index must be initialized.", __LINE__, error=error)
    CALL dbcsr_assert ("NOT", matrix%m%list_indexing,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "List index already exists?", __LINE__, error=error)
    CALL dbcsr_assert ("NOT", matrix%m%bcsc,&
         dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
         "Not support for BCSC yet.", __LINE__, error=error)
    CALL dbcsr_assert (matrix%m%nblks, "EQ", SIZE(matrix%m%col_i),&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Block count mismatch.", __LINE__, error=error)
    CALL dbcsr_assert (matrix%m%nblks, "EQ", SIZE(matrix%m%blk_p),&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Block count mismatch", __LINE__, error=error)
    !
    IF (matrix%m%local_indexing) THEN
       CALL dbcsr_assert (SIZE(matrix%m%row_p)-1, "EQ", matrix%m%nblkrows_local,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Local row index incorrectly sized.", __LINE__, error=error)
    ELSE
       CALL dbcsr_assert (SIZE(matrix%m%row_p)-1, "EQ", matrix%m%nblkrows_total,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Global row index incorrectly sized", __LINE__, error=error)
    ENDIF
    !
    matrix%m%list_indexing = .TRUE.
    !
    IF (matrix%m%local_indexing) THEN
       global_cols => array_data (matrix%m%global_cols)
       nrows = matrix%m%nblkrows_local
    ELSE
       nrows = matrix%m%nblkrows_total
    ENDIF
    !
    nblks = matrix%m%nblks
    ALLOCATE (blki(3, nblks))
    CALL dbcsr_expand_row_index (matrix%m%row_p, blki(1,:), nrows, nblks)
    IF (matrix%m%local_indexing) THEN
       ! If local indexing is enabled, then the rows but not the
       ! columns are already localized
       IF (dbg) THEN
          WRITE(*,*)routineN//" Making local columns"
          WRITE(*,'(10(1X,i7))')global_cols
          WRITE(*,*)'local'
          WRITE(*,'(10(1X,i7))')array_data(matrix%m%local_cols)
       ENDIF
       FORALL (blk = 1 : nblks)
          blki(2,blk) = global_cols(matrix%m%col_i(blk))
          blki(3,blk) = matrix%m%blk_p(blk)
       END FORALL
    ELSE
       IF (dbg) WRITE(*,*)routineN//" Not making local columns"
       FORALL (blk = 1 : nblks)
          blki(2,blk) = matrix%m%col_i(blk)
          blki(3,blk) = matrix%m%blk_p(blk)
       END FORALL
    ENDIF
    CALL dbcsr_clearfrom_index_array (matrix%m, dbcsr_slot_row_p)
    CALL dbcsr_clearfrom_index_array (matrix%m, dbcsr_slot_col_i)
    CALL dbcsr_clearfrom_index_array (matrix%m, dbcsr_slot_blk_p)
    nthreads = 0
    !$ nthreads = OMP_GET_MAX_THREADS ()
    !
    ! Reshuffle according to threads
    IF (nthreads .GT. 0 .AND. thread_redist) THEN
       td => array_data (dbcsr_distribution_thread_dist (dbcsr_distribution (matrix)))
       IF (matrix%m%local_indexing) THEN
          local_rows => array_data (matrix%m%local_rows)
       ENDIF
       !$OMP PARALLEL DEFAULT (none) &
       !$OMP PRIVATE (my_cnt, ithread, blk) &
       !$OMP SHARED (td, blki, nthreads, thr_c, nblks, matrix, error, local_rows)
       !
       ithread = 0
       !$ ithread = OMP_GET_THREAD_NUM ()
       !$OMP MASTER
       !$ nthreads = OMP_GET_NUM_THREADS ()
       CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_thr_c,&
            reservation=nthreads+1, extra=3*nblks, error=error)
       thr_c => matrix%m%thr_c
       CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_coo_l,&
            reservation=3*nblks, error=error)
       !$OMP END MASTER
       !$OMP BARRIER
       my_cnt = 0
       IF (matrix%m%local_indexing) THEN
          my_cnt = COUNT (td(local_rows(blki(1,:))) .EQ. ithread)
       ELSE
          my_cnt = COUNT (td(blki(1,:)) .EQ. ithread)
       ENDIF
       !DO blk = 1, nblks
       !   IF (td(blki(1, blk)) .EQ. ithread) my_cnt = my_cnt+1
       !ENDDO
       thr_c(ithread+1) = my_cnt
       !$OMP BARRIER
       !$OMP MASTER
       CALL dbcsr_build_row_index_inplace (thr_c, nthreads)
       !$OMP END MASTER
       !$OMP BARRIER
       my_cnt = (thr_c(ithread+1)+1) * 3 -2
       IF (matrix%m%local_indexing) THEN
          DO blk = 1, nblks
             IF (td(local_rows(blki(1, blk))) .EQ. ithread) THEN
                matrix%m%coo_l(my_cnt:my_cnt+2) = blki(1:3, blk)
                my_cnt = my_cnt + 3
             ENDIF
          ENDDO
       ELSE
          DO blk = 1, nblks
             IF (td(blki(1, blk)) .EQ. ithread) THEN
                matrix%m%coo_l(my_cnt:my_cnt+2) = blki(1:3, blk)
                my_cnt = my_cnt + 3
             ENDIF
          ENDDO
       ENDIF
       !$OMP END PARALLEL
    ELSE
       ! Small price to pay for avoiding infinite recursions.
       DO blk = 2, nblks
          IF (       blki(1, blk) .EQ. blki(1, blk-1)&
               .AND. blki(2, blk) .EQ. blki(2, blk-1)) &
               ! Weird assertion to give some idea of the two blocks.
               CALL dbcsr_assert (-blki(1,blk), "EQ", blki(2,blk),&
               dbcsr_fatal_level, dbcsr_caller_error, routineN,&
      "Should not have duplicate matrix blocks. (-row, col) is duplicated.",&
               __LINE__, error=error)
       ENDDO
       !
       IF (nblks>0) THEN
          CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_coo_l,&
               DATA=RESHAPE(blki, (/3*nblks/)), error=error)
       ELSE
          CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_coo_l,&
               DATA=zero_len_array, error=error)
       ENDIF
    ENDIF

    DEALLOCATE (blki)
    !
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_make_index_list



! *****************************************************************************
!> \brief Compacts an index.
!> \param[in,out] matrix       matrix for which to make canonical index
!> \param[in,out] error        error
! *****************************************************************************
  SUBROUTINE dbcsr_index_compact(matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_index_compact', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, new_size, &
                                                size_blk_p, size_col_i, &
                                                size_coo_l, size_row_p, &
                                                size_thr_c
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: blk_p, col_i, coo_l, meta, &
                                                row_p, thr_c
    LOGICAL                                  :: compact, has_blk_p, &
                                                has_col_i, has_coo_l, &
                                                has_row_p, has_thr_c

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error=error)
    !> Ensures the index pointers are set.
    CALL dbcsr_repoint_index (matrix%m)
    ! Check that compaction is even needed.
    has_row_p = ASSOCIATED (matrix%m%row_p)
    IF (has_row_p) THEN
       size_row_p = SIZE(matrix%m%row_p)
    ELSE
       size_row_p = 0
    ENDIF
    has_col_i = ASSOCIATED (matrix%m%col_i)
    IF (has_col_i) THEN
       size_col_i = SIZE(matrix%m%col_i)
    ELSE
       size_col_i = 0
    ENDIF
    has_blk_p = ASSOCIATED (matrix%m%blk_p)
    IF (has_blk_p) THEN
       size_blk_p = SIZE(matrix%m%blk_p)
    ELSE
       size_blk_p = 0
    ENDIF
    has_thr_c = ASSOCIATED (matrix%m%thr_c)
    IF (has_thr_c) THEN
       size_thr_c = SIZE(matrix%m%thr_c)
    ELSE
       size_thr_c = 0
    ENDIF
    has_coo_l = ASSOCIATED (matrix%m%coo_l)
    IF (has_coo_l) THEN
       size_coo_l = SIZE(matrix%m%coo_l)
    ELSE
       size_coo_l = 0
    ENDIF
    !
    new_size = dbcsr_num_slots +&
         size_row_p + size_col_i + size_blk_p + size_thr_c + size_coo_l
    compact = new_size .LT. SIZE(matrix%m%index)
    IF (compact) THEN
       ! Store old index arrays.
       IF (has_row_p) THEN
          ALLOCATE (row_p(size_row_p))
          row_p(:) = matrix%m%row_p(:)
       ENDIF
       IF (has_col_i) THEN
          ALLOCATE (col_i(size_col_i))
          col_i(:) = matrix%m%col_i(:)
       ENDIF
       IF (has_blk_p) THEN
          ALLOCATE (blk_p(size_blk_p))
          blk_p(:) = matrix%m%blk_p(:)
       ENDIF
       IF (has_thr_c) THEN
          ALLOCATE (thr_c(size_thr_c))
          thr_c(:) = matrix%m%thr_c(:)
       ENDIF
       IF (has_coo_l) THEN
          ALLOCATE (coo_l(size_coo_l))
          coo_l(:) = matrix%m%coo_l(:)
       ENDIF
       ALLOCATE (meta(dbcsr_num_slots))
       meta(:) = matrix%m%index(1:dbcsr_num_slots)
       ! Clear the index.
       CALL memory_deallocate (matrix%m%index,&
            dbcsr_get_index_memory_type (matrix), error=error)
       NULLIFY (matrix%m%index)
       CALL memory_allocate (matrix%m%index, new_size,&
            dbcsr_get_index_memory_type (matrix), error=error)
       !
       ! Now copy the old index arrays into the index. We must not
       ! copy the positions of the old pointers.
       matrix%m%index(1:dbcsr_meta_size) = meta(1:dbcsr_meta_size)
       matrix%m%index(dbcsr_meta_size+1:) = 0
       matrix%m%index(dbcsr_slot_size) = dbcsr_num_slots
       IF (has_thr_c) THEN
          CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_thr_c, thr_c,&
               error=error)
          DEALLOCATE (thr_c)
       ENDIF
       IF (has_row_p) THEN
          CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_row_p, row_p,&
               error=error)
          DEALLOCATE (row_p)
       ENDIF
       IF (has_col_i) THEN
          CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_col_i, col_i,&
               error=error)
          DEALLOCATE (col_i)
       ENDIF
       IF (has_blk_p) THEN
          CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_blk_p, blk_p,&
               error=error)
          DEALLOCATE (blk_p)
       ENDIF
       IF (has_coo_l) THEN
          CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_coo_l, coo_l,&
               error=error)
          DEALLOCATE (coo_l)
       ENDIF
       DEALLOCATE (meta)
       IF (careful_mod) THEN
          ! This is pretty strong but it should be true.
          CALL dbcsr_assert (matrix%m%index(dbcsr_slot_size), "EQ", new_size,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Unexpected index size.", __LINE__, error=error)
          CALL dbcsr_assert (SIZE(matrix%m%index), "EQ", new_size,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Unexpected index size.", __LINE__, error=error)
       ENDIF
       CALL dbcsr_repoint_index (matrix%m)
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_index_compact

! *****************************************************************************
!> \brief Calculates the checksum of an index.
!> \param[in] matrix           matrix for which to make canonical index
!> \param[in,out] error        error
! *****************************************************************************
  SUBROUTINE dbcsr_index_checksum(matrix, checksum, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(OUT)                     :: checksum
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_index_checksum', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error=error)
    !
    checksum = joaat_hash ( (/ joaat_hash (matrix%m%row_p),&
                               joaat_hash (matrix%m%col_i),&
                               joaat_hash (matrix%m%blk_p) /) )
    !
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_index_checksum

  FUNCTION dbcsr_has_local_row_index (matrix) RESULT (local_indexing)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    LOGICAL                                  :: local_indexing

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_has_local_row_index', &
      routineP = moduleN//':'//routineN

    local_indexing = matrix%m%local_indexing
  END FUNCTION dbcsr_has_local_row_index

END MODULE dbcsr_index_operations

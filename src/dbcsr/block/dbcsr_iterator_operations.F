!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
                                             array_exists
  USE dbcsr_data_methods,              ONLY: dbcsr_data_hold,&
                                             dbcsr_data_release,&
                                             dbcsr_data_set_pointer,&
                                             dbcsr_get_data
  USE dbcsr_dist_methods,              ONLY: dbcsr_distribution_has_threads,&
                                             dbcsr_distribution_num_threads,&
                                             dbcsr_distribution_thread_dist
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level,&
                                             dbcsr_internal_error,&
                                             dbcsr_warning_level,&
                                             dbcsr_wrong_args_error
  USE dbcsr_methods,                   ONLY: dbcsr_distribution
  USE dbcsr_ptr_util,                  ONLY: pointer_rank_remap2
  USE dbcsr_toollib,                   ONLY: swap
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_distribution_obj,&
                                             dbcsr_iterator,&
                                             dbcsr_obj,&
                                             dbcsr_scalar_type
  USE kinds,                           ONLY: real_4,&
                                             real_8

  !$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads, omp_in_parallel

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_iterator_operations'

  INTEGER, PRIVATE, POINTER, SAVE      :: common_int_pointer


  PUBLIC :: dbcsr_iterator_start, dbcsr_iterator_stop

  PUBLIC :: dbcsr_iterator_blocks_left, dbcsr_iterator_next_block


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





  LOGICAL, PARAMETER :: careful_mod  = .FALSE.

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
! Iterator functions
! *****************************************************************************


! *****************************************************************************
!> \brief Sets up an iterator
!> \param[out] iterator   the iterator
!> \param[in] matrix      DBCSR matrix
!> \param[in] shared      (optional) The matrix is shared between several
!>                        iterators. Default is .TRUE.
!> \param[in] dynamic     (optional) Threads are given blocks regardless of
!>                        the thread distribution; default is .FALSE.
!> \param[in] dynamic_byrows  (optional) Threads are given blocks regardless
!>                            of the thread distribution, but still grouped
!>                            by rows; default is .FALSE.
!> \param[in] contiguous_pointers  (optional) Whether returned pointers need
!>                                 to be contiguous; default is FALSE.
!> \param[in] read_only            (optional) User promises not to change
!>                                 returned data; default is FALSE
!>
!> \par Contiguous pointers
!> Contiguous pointers may incur reallocation penalties but enable quick
!> passing of arrays to routines with unspecified interfaces (i.e., direct
!> calls to BLACS or MPI).
!> \par Threading
!>      The TYPE(dbcsr_iterator) variable should be thread-private.
!> \par
!>      The iterator has several modes of operation when used with
!>      OpenMP. Two options can be set to influence the behavior.
!> \par Threading: shared vs. non-shared
!>      The "shared" flag specifies that several threads will be
!>      iterating through the same matrix.
!>      - Sharing is the default when called from an active parallel
!>        region. In the shared mode no two threads will receive the
!>        same block; i.e., the work is split among the threads.
!>      - If each (or one) thread needs to iterator through all blocks
!>        then shared should be set to .FALSE.. (E.g., when called
!>        from an enclosing MASTER region or when each thread has its
!>        own matrix.)
!>      - It is safe to use an iterator in non-shared mode with only
!>        one thread.  No thread synchronization constructs are used
!>        in this case)
!> \par Threading in shared mode
!>      When in shared mode there are three possibilities to select
!>      how the blocks are distributed to the threads.
!>      <DL>
!>      <DT>Thread distribution</DT>
!>      <DD>The default is to use the thread distribution. The thread
!>          distribution statically maps rows to threads and should be
!>          used whenever retaining a consistent mapping among
!>          subsequent iterations is important.</DD>
!>      <DT>Dynamic scheduling</DT>
!>      <DD>If the dynamic flag is .TRUE., then blocks are given to
!>          threads dynamically. By default the assignment is grouped
!>          by rows (to minimize synchronization); however, if the
!>          dynamic_byrows flag is .FALSE. then every block is
!>          assigned dynamically.</DD></DL>
! *****************************************************************************
  SUBROUTINE dbcsr_iterator_start (iterator, matrix, shared, dynamic,&
       dynamic_byrows, contiguous_pointers, read_only)
    TYPE(dbcsr_iterator), INTENT(OUT)        :: iterator
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: shared, dynamic, &
                                                dynamic_byrows, &
                                                contiguous_pointers, read_only

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_iterator_start', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    iterator%shared = .TRUE.
!$  iterator%shared = omp_in_parallel()
!ub!$  CALL dbcsr_assert (omp_in_parallel(), "IMP", PRESENT (shared),&
!ub!$       dbcsr_warning_level, dbcsr_caller_error, routineN,&
!ub!$       "SHARED must be specified in OMP invocation.", __LINE__, error)
    IF (PRESENT (shared)) iterator%shared = shared
    iterator%dynamic = .TRUE.
!$  iterator%dynamic = .FALSE.
!$  IF (PRESENT (dynamic)) iterator%dynamic = dynamic
    IF (PRESENT (dynamic_byrows)) THEN
       iterator%dynamic_byrows = dynamic_byrows
       IF (iterator%dynamic_byrows) iterator%dynamic = .TRUE.
    ELSE
       iterator%dynamic_byrows = iterator%dynamic
       !$iterator%dynamic_byrows = iterator%dynamic
    ENDIF
!$  IF (.NOT.iterator%shared) THEN
!$     iterator%dynamic = .FALSE.
!$  ENDIF
    dist = dbcsr_distribution (matrix)
!$  CALL dbcsr_assert (dbcsr_distribution_has_threads(dist),&
!$       dbcsr_warning_level, dbcsr_internal_error, routineN,&
!$       "Thread distribution should be defined for OpenMP.",&
!$       __LINE__, error=error)
    CALL dbcsr_assert (.NOT.iterator%dynamic, "IMP",&
         dbcsr_distribution_has_threads(dist),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Thread distribution must be defined for non-dynamic iterator.",&
         __LINE__, error=error)
!$  IF (omp_in_parallel ()) &
!$     CALL dbcsr_assert (omp_get_num_threads(), "EQ",&
!$          dbcsr_distribution_num_threads(dist),&
!$          dbcsr_fatal_level, dbcsr_internal_error, routineN,&
!$          "Number of threads has changed!", __LINE__, error=error)
    !Synchronize the positions
    NULLIFY (iterator%common_pos)
    IF (iterator%dynamic) THEN
       ! All threads point into the master thread's data space
       ! (temporarily using the common_int_pointer variable). This is
       ! not the nicest OpenMP way of doing this but it is also not
       ! explicity forbidden.
       !
       !$OMP BARRIER
       !$OMP MASTER
       ALLOCATE (iterator%common_pos)
       common_int_pointer => iterator%common_pos
       common_int_pointer = 0
       !$OMP FLUSH (common_int_pointer)
       !$OMP END MASTER
       !$OMP BARRIER
       IF (.NOT. ASSOCIATED (iterator%common_pos)) THEN
          iterator%common_pos => common_int_pointer
       ENDIF
       !$OMP BARRIER
    ENDIF
    !
    IF (PRESENT (contiguous_pointers)) THEN
       iterator%contiguous_pointers = contiguous_pointers
    ELSE
       iterator%contiguous_pointers = .TRUE.
    ENDIF
    IF (PRESENT (read_only)) THEN
       iterator%read_only = read_only
    ELSE
       iterator%read_only = .FALSE.
    ENDIF
    iterator%row = 0
    iterator%pos = 0
    iterator%rbs => array_data (matrix%m%row_blk_size)
    iterator%cbs => array_data (matrix%m%col_blk_size)
    iterator%roff => array_data (matrix%m%row_blk_offset)
    iterator%coff => array_data (matrix%m%col_blk_offset)

    iterator%local_indexing = matrix%m%local_indexing
    !CALL dbcsr_assert (iterator%local_indexing, "IMP",&
    !     iterator%dynamic,&
    !     dbcsr_fatal_level, dbcsr_caller_error, routineN,&
    !     "Locally-indexed matrices can only have a dynamic iterator.",&
    !     __LINE__, error=error)
    CALL dbcsr_assert (iterator%local_indexing, "IMP",&
         array_exists (matrix%m%local_rows),&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Local rows mapping array should exist when local indexing is used.",&
         __LINE__, error=error)
    CALL dbcsr_assert (iterator%local_indexing, "IMP",&
         array_exists (matrix%m%global_rows),&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Global rows mapping array should exist when local indexing is used.",&
         __LINE__, error=error)
    iterator%global_rows => array_data(matrix%m%global_rows)
    iterator%local_rows => array_data(matrix%m%local_rows)

    iterator%transpose = .FALSE.!matrix%m%transpose
    iterator%nblks = matrix%m%nblks
    IF (iterator%transpose) THEN
       iterator%nblkrows_total = matrix%m%nblkcols_total
    ELSE
       iterator%nblkrows_total = matrix%m%nblkrows_total
    ENDIF

    iterator%row_p => matrix%m%row_p
    iterator%col_i => matrix%m%col_i
    iterator%blk_p => matrix%m%blk_p
    !$OMP CRITICAL (crit_data)
    iterator%data_area = matrix%m%data_area
    CALL dbcsr_data_hold (iterator%data_area)
    !$OMP END CRITICAL (crit_data)
    iterator%row_size = 0
    IF (.NOT. iterator%dynamic) THEN
       iterator%tdist => array_data (dbcsr_distribution_thread_dist (dist))
    ELSE
       NULLIFY (iterator%tdist)
    ENDIF
    !$ IF (iterator%dynamic) THEN
          !$OMP SINGLE
    !$    IF (iterator%dynamic_byrows) THEN
    !$       iterator%common_pos = omp_get_num_threads()
    !$    ENDIF
          !$OMP END SINGLE
    !$    CALL dbcsr_iterator_seek (iterator, omp_get_thread_num()+1)
    !$ ELSE
    CALL dbcsr_iterator_seek (iterator, 1)
    !$ ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_iterator_start

! *****************************************************************************
!> \brief Stops up an iterator
!> \param[out] iterator   the iterator
! *****************************************************************************
  SUBROUTINE dbcsr_iterator_stop (iterator)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_iterator_stop', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    iterator%row = 0
    iterator%pos = 0

    NULLIFY (iterator%tdist)
    !$OMP CRITICAL (crit_data)
    CALL dbcsr_data_release (iterator%data_area)
    !$OMP END CRITICAL (crit_data)
    IF (iterator%dynamic) THEN
       !$OMP BARRIER
       !$OMP MASTER
       common_int_pointer => iterator%common_pos
       DEALLOCATE (common_int_pointer)
       !$OMP FLUSH (common_int_pointer)
       !$OMP END MASTER
       NULLIFY (iterator%common_pos)
       !$OMP BARRIER
    ENDIF
    IF (iterator%local_indexing) THEN
       NULLIFY (iterator%local_rows)
       NULLIFY (iterator%global_rows)
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
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
    INTEGER, DIMENSION(:), INTENT(IN)        :: blk_p

    CHARACTER(len=*), PARAMETER :: routineN = 'find_first_valid_block', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------
!IF (pos .LT. 1) pos = 1

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
!> \param[in] pos            current position
!> \param[in,out] row        input: current row; output: the row corresponding
!>                           to the position
!> \param[in] maxrows        maxmimum row
!> \param[in] row_p          row pointers
! *****************************************************************************
  PURE SUBROUTINE find_proper_row (pos, row, maxrows, row_p)
    INTEGER, INTENT(IN)                      :: pos
    INTEGER, INTENT(INOUT)                   :: row
    INTEGER, INTENT(IN)                      :: maxrows
    INTEGER, DIMENSION(:), INTENT(IN)        :: row_p

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
!> \param local2global ...
! *****************************************************************************
  PURE SUBROUTINE find_proper_position (pos, row, maxpos, maxrows,&
       blk_p, row_p, tdist, tid, local2global)
    INTEGER, INTENT(INOUT)                   :: pos, row
    INTEGER, INTENT(IN)                      :: maxpos, maxrows
    INTEGER, DIMENSION(:), INTENT(IN)        :: blk_p, row_p
    INTEGER, DIMENSION(1:maxrows), &
      INTENT(IN), OPTIONAL                   :: tdist
    INTEGER, INTENT(IN), OPTIONAL            :: tid
    INTEGER, DIMENSION(1:*), INTENT(IN), &
      OPTIONAL                               :: local2global

    CHARACTER(len=*), PARAMETER :: routineN = 'find_proper_position', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: local, row_inrange, row_ok

!   ---------------------------------------------------------------------------

    local = PRESENT (local2global)
    IF (maxpos .GE. 1) THEN
       !IF (pos.EQ.0) pos = 1
       CALL find_first_valid_block (pos, maxpos, blk_p)
       CALL find_proper_row(pos, row, maxrows, row_p)
       row_inrange = row .NE. 0 .AND. row .LE. maxrows
       row_ok = row_inrange
!$     IF (present (tdist) .AND. PRESENT (tid) .AND. row_inrange) THEN
!$        IF (.NOT. local) THEN
!$           row_ok = tdist(row) .EQ. tid
!$        ELSE
!$           row_ok = tdist(local2global(row)) .EQ. tid
!$        ENDIF
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
!$        IF (.NOT. local) THEN
!$           row_ok = tdist(row) .EQ. tid
!$        ELSE
!$           row_ok = tdist(local2global(row)) .EQ. tid
!$        ENDIF
!$     ENDIF
       ENDDO
       IF (row .EQ. 0) pos = 0
    ELSE
       pos = 0
       row = 0
    ENDIF
  END SUBROUTINE find_proper_position

! *****************************************************************************
!> \brief Advances to the next block
!> \param[in,out] iterator     the iterator
! *****************************************************************************
  SUBROUTINE iterator_advance (iterator)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator

    INTEGER                                  :: ithread, my_old_row, &
                                                next_row, p
    LOGICAL                                  :: advance, jumped_row

!   ---------------------------------------------------------------------------

    ithread = 0
!$  ithread = omp_get_thread_num()
    IF (iterator%dynamic .AND. iterator%shared) THEN
       IF (iterator%dynamic_byrows) THEN
          ! In this case common_row holds the last assigned row.
          !
          ! The idea is to advance in this thread's row. If it gets bumped
          ! into the next row, then we have to find the correct one.
          advance = .TRUE.
          DO WHILE (advance)
             iterator%pos = iterator%pos + 1
             my_old_row = iterator%row
             CALL find_proper_position_caller(iterator)
             jumped_row = iterator%row .GT. my_old_row
             advance = jumped_row
             IF (jumped_row) THEN
                !$OMP CRITICAL (crit_common_pos)
                ! Set the common_pos to the next available row.
                next_row = MAX (iterator%common_pos + 1,  iterator%row)
                next_row = next_row - iterator%common_pos
                iterator%common_pos = iterator%common_pos + next_row
                iterator%row = iterator%common_pos
                !$OMP END CRITICAL (crit_common_pos)
                IF (iterator%row .GT. iterator%nblkrows_total) THEN
                   iterator%pos = 0
                   iterator%row = 0
                   advance = .FALSE.
                ELSE
                   ! To be incremented in the next loop.
                   IF (.NOT. iterator%local_indexing) THEN
                      iterator%pos = iterator%row_p(iterator%row)
                   ELSE
                      iterator%pos = iterator%row_p(iterator%global_rows(iterator%row))
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ELSE
          ! In this case common_pos holds the last-assigned block.
          !
          !$OMP CRITICAL (crit_common_pos)
          iterator%common_pos = iterator%common_pos + 1
          iterator%pos = iterator%common_pos
          CALL find_proper_position_caller(iterator)
          p = iterator%pos
          iterator%common_pos = MAX(iterator%common_pos, p)
          !$OMP END CRITICAL (crit_common_pos)
       ENDIF
    ELSEIF (iterator%shared) THEN
       iterator%pos = iterator%pos+1
       ithread = 0
!$     ithread = OMP_GET_THREAD_NUM()
       CALL find_proper_position_caller(iterator, use_ithread=ithread)
    ELSE
       iterator%pos = iterator%pos+1
       CALL find_proper_position_caller(iterator)
    ENDIF
  END SUBROUTINE iterator_advance


! *****************************************************************************
!> \brief ...
!> \param iterator ...
!> \param use_ithread ...
! *****************************************************************************
  SUBROUTINE find_proper_position_caller (iterator, use_ithread)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator
    INTEGER, INTENT(in), OPTIONAL            :: use_ithread

    INTEGER                                  :: post_pos, post_row, pre_pos, &
                                                pre_row, r

    pre_row = iterator%row
    pre_pos = iterator%pos
    IF (.NOT. iterator%local_indexing) THEN
       IF (PRESENT (use_ithread)) THEN
          CALL find_proper_position(iterator%pos,&
               iterator%row, iterator%nblks, iterator%nblkrows_total,&
               iterator%blk_p, iterator%row_p,&
               tdist=iterator%tdist, tid=use_ithread)
       ELSE
          CALL find_proper_position(iterator%pos,&
               iterator%row, iterator%nblks, iterator%nblkrows_total,&
               iterator%blk_p, iterator%row_p)
       ENDIF
    ELSE
       IF (iterator%row .GT. 0) THEN
          r = iterator%global_rows(iterator%row)
       ELSE
          r = 0
       ENDIF
       IF (PRESENT(use_ithread)) THEN
          CALL find_proper_position(iterator%pos,&
               r, iterator%nblks, SIZE(iterator%local_rows),&
               iterator%blk_p, iterator%row_p,&
               local2global=iterator%local_rows,&
               tdist=iterator%tdist, tid=use_ithread)
       ELSE
          CALL find_proper_position(iterator%pos,&
               r, iterator%nblks, SIZE(iterator%local_rows),&
               iterator%blk_p, iterator%row_p,&
               local2global=iterator%local_rows)
       ENDIF
       IF (r .GT. 0) THEN
          iterator%row = iterator%local_rows(r)
       ELSE
          iterator%row = 0
       ENDIF
    ENDIF
    post_row = iterator%row
    post_pos = iterator%pos
  END SUBROUTINE find_proper_position_caller

! *****************************************************************************
!> \brief Updates the row info stored in the iterator
!> \param[in,out] iterator     the iterator
!> \note Added to handle the complexity introduced with the transpose
! *****************************************************************************
  PURE SUBROUTINE update_row_info (iterator)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator

    IF (iterator%row .GT. 0) THEN
       IF (iterator%transpose) THEN
          iterator%row_size = iterator%cbs(iterator%row)
          iterator%row_offset = iterator%coff(iterator%row)
       ELSE
          iterator%row_size = iterator%rbs(iterator%row)
          iterator%row_offset = iterator%roff(iterator%row)
       ENDIF
    ENDIF
  END SUBROUTINE update_row_info

! *****************************************************************************
!> \brief Places the iterator to the desired row.
!> \param[in,out] iterator     the iterator
!> \param[in] row              seek to this row
! *****************************************************************************
  SUBROUTINE dbcsr_iterator_seek (iterator, row)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator
    INTEGER, INTENT(in)                      :: row

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_iterator_seek', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    IF (iterator%nblks .GT. 0 .AND. row .LE. iterator%nblkrows_total) THEN
       iterator%row = row
       ! This line is replaced because iterator_advance increments the block
       ! number
       !iterator%pos = iterator%row_p(row)+1
       iterator%pos = iterator%row_p(row)  ! +1-1
       CALL iterator_advance (iterator)
       CALL update_row_info (iterator)
    ELSE
       iterator%row = 0
       iterator%pos = 0
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_iterator_seek

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

    INTEGER                                  :: bp, error_handle
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set (routineN, error_handle, error)
    IF (iterator%pos .LE. iterator%nblks&
         .AND. iterator%pos .NE. 0) THEN
       row = iterator%row
       column = iterator%col_i(iterator%pos)
       IF (iterator%transpose) CALL swap (row, column)
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
       CALL iterator_advance (iterator)
       CALL update_row_info (iterator)
    ELSE
       row = 0
       column = 0
    ENDIF
    IF (careful_mod) CALL dbcsr_error_stop (error_handle, error)
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

    INTEGER                                  :: blk_p, block_col_size, &
                                                block_row_size, bp, csize, &
                                                error_handle, nze, rsize
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set (routineN, error_handle, error)

    IF (iterator%pos .LE. iterator%nblks&
         .AND. iterator%pos .NE. 0) THEN
       row = iterator%row
       column = iterator%col_i(iterator%pos)
       IF (iterator%transpose) CALL swap (row, column)
       blk_p = iterator%blk_p(iterator%pos)
       transposed = blk_p .LT. 0
       bp = ABS (blk_p)
       rsize = iterator%row_size
       csize = iterator%cbs(column)
       block_row_size = rsize
       block_col_size = csize
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
       CALL iterator_advance (iterator)
       CALL update_row_info (iterator)
    ELSE
       row = 0
       column = 0
       IF (PRESENT (block_number)) block_number = 0
    ENDIF

    IF (careful_mod) CALL dbcsr_error_stop (error_handle, error)

  END SUBROUTINE iterator_next_area_block



! *****************************************************************************
!> \brief Returns whether there any blocks left in the iterator.
!> \param[in] iterator       the iterator
!> \retval blocks_left ...
!> \param[in]
! *****************************************************************************
  PURE FUNCTION dbcsr_iterator_blocks_left (iterator) RESULT (blocks_left)
    TYPE(dbcsr_iterator), INTENT(IN)         :: iterator
    LOGICAL                                  :: blocks_left

    blocks_left = iterator%pos .NE. 0
  END FUNCTION dbcsr_iterator_blocks_left




#include "dbcsr_iterator_operations_d.f90"
#include "dbcsr_iterator_operations_z.f90"
#include "dbcsr_iterator_operations_s.f90"
#include "dbcsr_iterator_operations_c.f90"


END MODULE dbcsr_iterator_operations

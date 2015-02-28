!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Routines related to DBCSR distributions
!> \author  Ole Schuett
!> \date    10-2014
! *****************************************************************************
MODULE dbcsr_dist_methods
  USE array_types,                     ONLY: array_data,&
                                             array_hold,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release,&
                                             array_size
  USE dbcsr_config,                    ONLY: comm_thread_load
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level,&
                                             dbcsr_internal_error,&
                                             dbcsr_wrong_args_error
  USE dbcsr_methods,                   ONLY: dbcsr_col_block_sizes,&
                                             dbcsr_distribution,&
                                             dbcsr_distribution_init,&
                                             dbcsr_distribution_release,&
                                             dbcsr_nblkcols_total,&
                                             dbcsr_nblkrows_total,&
                                             dbcsr_row_block_sizes
  USE dbcsr_mp_methods,                ONLY: dbcsr_mp_hold,&
                                             dbcsr_mp_mypcol,&
                                             dbcsr_mp_myprow,&
                                             dbcsr_mp_npcols,&
                                             dbcsr_mp_nprows
  USE dbcsr_toollib,                   ONLY: sort
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_mp_obj,&
                                             dbcsr_obj
  USE kinds,                           ONLY: int_8,&
                                             sp
  USE min_heap,                        ONLY: heap_fill,&
                                             heap_get_first,&
                                             heap_new,&
                                             heap_release,&
                                             heap_reset_first,&
                                             heap_t

  !$ USE OMP_LIB
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_dist_methods'

  PUBLIC :: dbcsr_distribution_new, dbcsr_distribution_hold,&
            dbcsr_distribution_release, dbcsr_distribution_init,&
            dbcsr_distribution_mp, dbcsr_distribution_processor,&
            dbcsr_distribution_nrows, dbcsr_distribution_ncols,&
            dbcsr_distribution_row_dist, dbcsr_distribution_col_dist,&
            dbcsr_distribution_nlocal_rows, dbcsr_distribution_nlocal_cols,&
            dbcsr_distribution_local_rows, dbcsr_distribution_local_cols,&
            dbcsr_distribution_thread_dist, dbcsr_distribution_has_threads,&
            dbcsr_distribution_make_threads, dbcsr_distribution_no_threads,&
            dbcsr_distribution_num_threads, dbcsr_distribution_data_max_size_allproc

CONTAINS

! *****************************************************************************
!> \brief Creates new distribution
!> \param[out] dist           distribution
!> \param[in] mp_env          multiprocessing environment
!> \param row_dist row distribution
!> \param col_dist column distribution
!> \param local_rows ...
!> \param local_cols ...
! *****************************************************************************
  SUBROUTINE dbcsr_distribution_new(dist, mp_env, row_dist, col_dist,&
       local_rows, local_cols)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    TYPE(array_i1d_obj), INTENT(IN)          :: row_dist, col_dist
    TYPE(array_i1d_obj), INTENT(IN), &
      OPTIONAL                               :: local_rows, local_cols

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_distribution_new', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, mypcoor, seq
    INTEGER, DIMENSION(:), POINTER           :: dd, ld
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    ALLOCATE (dist%d)
    dist%d%refcount = 1
    dist%d%row_dist = row_dist
    CALL array_hold (dist%d%row_dist)
    dist%d%col_dist = col_dist
    CALL array_hold (dist%d%col_dist)
    dist%d%mp_env = mp_env
    CALL dbcsr_mp_hold (dist%d%mp_env)
    ! Verify given process row distribution.
    dd => array_data (row_dist)
    CALL dbcsr_assert (MAXVAL (dd), "LT", dbcsr_mp_nprows (mp_env),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "A process row is too big for process grid",&
         __LINE__, error=error)
    ! Verify given process column distribution.
    dd => array_data (col_dist)
    CALL dbcsr_assert (MAXVAL (dd), "LT", dbcsr_mp_npcols (mp_env),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "A process column is too big for process grid",&
         __LINE__, error=error)
    IF (PRESENT (local_rows)) THEN
       dist%d%local_rows = local_rows
       CALL array_hold (dist%d%local_rows)
    ELSE
       dd => array_data (row_dist)
       mypcoor = dbcsr_mp_myprow (mp_env)
       i = COUNT (dd.EQ.mypcoor)
       NULLIFY (ld)
       ALLOCATE (ld(i))
       seq = 1
       DO i = 1, array_size (row_dist)
          IF (dd(i) .EQ. mypcoor) THEN
             ld(seq) = i
             seq = seq+1
          ENDIF
       ENDDO
       CALL array_new (dist%d%local_rows, ld, gift=.TRUE.)
    ENDIF
    IF (PRESENT (local_cols)) THEN
       dist%d%local_cols = local_cols
       CALL array_hold (dist%d%local_cols)
    ELSE
       dd => array_data (col_dist)
       mypcoor = dbcsr_mp_mypcol (mp_env)
       i = COUNT (dd.EQ.mypcoor)
       NULLIFY (ld)
       ALLOCATE (ld(i))
       seq = 1
       DO i = 1, array_size (col_dist)
          IF (dd(i) .EQ. mypcoor) THEN
             ld(seq) = i
             seq = seq+1
          ENDIF
       ENDDO
       CALL array_new (dist%d%local_cols, ld, gift=.TRUE.)
    ENDIF
    dist%d%num_threads = 1
!$  dist%d%num_threads = OMP_GET_MAX_THREADS()
    dist%d%has_thread_dist = .FALSE.
    CALL array_nullify (dist%d%thread_dist)
    CALL array_nullify (dist%d%row_map)
    CALL array_nullify (dist%d%col_map)
    NULLIFY (dist%d%other_l_rows)
    NULLIFY (dist%d%other_l_cols)
    dist%d%has_other_l_rows = .FALSE.
    dist%d%has_other_l_cols = .FALSE.
    CALL array_nullify (dist%d%global_row_map)
    CALL array_nullify (dist%d%global_col_map)
    dist%d%has_global_row_map = .FALSE.
    dist%d%has_global_col_map = .FALSE.
  END SUBROUTINE dbcsr_distribution_new

! *****************************************************************************
!> \brief Marks another use of the distribution
!> \param dist ...
! *****************************************************************************
  SUBROUTINE dbcsr_distribution_hold(dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist

!   ---------------------------------------------------------------------------

    dist%d%refcount = dist%d%refcount + 1
  END SUBROUTINE dbcsr_distribution_hold

! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \retval mp_env ...
! *****************************************************************************
  FUNCTION dbcsr_distribution_mp(dist) RESULT (mp_env)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(dbcsr_mp_obj)                       :: mp_env

!   ---------------------------------------------------------------------------

    mp_env = dist%d%mp_env
  END FUNCTION dbcsr_distribution_mp
! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \retval nrows ...
! *****************************************************************************
  PURE FUNCTION dbcsr_distribution_nrows(dist) RESULT (nrows)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER                                  :: nrows

    nrows = array_size (dist%d%row_dist)
  END FUNCTION dbcsr_distribution_nrows
! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \retval ncols ...
! *****************************************************************************
  PURE FUNCTION dbcsr_distribution_ncols(dist) RESULT (ncols)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER                                  :: ncols

    ncols = array_size (dist%d%col_dist)
  END FUNCTION dbcsr_distribution_ncols
! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \retval row_dist ...
! *****************************************************************************
  FUNCTION dbcsr_distribution_row_dist(dist) RESULT (row_dist)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: row_dist

!   ---------------------------------------------------------------------------

    row_dist = dist%d%row_dist
  END FUNCTION dbcsr_distribution_row_dist
! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \retval col_dist ...
! *****************************************************************************
  FUNCTION dbcsr_distribution_col_dist(dist) RESULT (col_dist)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: col_dist

!   ---------------------------------------------------------------------------

    col_dist = dist%d%col_dist
  END FUNCTION dbcsr_distribution_col_dist

! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \retval nlocalrows ...
! *****************************************************************************
  PURE FUNCTION dbcsr_distribution_nlocal_rows(dist) RESULT (nlocalrows)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER                                  :: nlocalrows

    nlocalrows = array_size (dist%d%local_rows)
  END FUNCTION dbcsr_distribution_nlocal_rows
! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \retval nlocalcols ...
! *****************************************************************************
  PURE FUNCTION dbcsr_distribution_nlocal_cols(dist) RESULT (nlocalcols)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER                                  :: nlocalcols

    nlocalcols = array_size (dist%d%local_cols)
  END FUNCTION dbcsr_distribution_nlocal_cols
! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \retval local_rows ...
! *****************************************************************************
  FUNCTION dbcsr_distribution_local_rows(dist) RESULT (local_rows)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: local_rows

    local_rows = dist%d%local_rows
  END FUNCTION dbcsr_distribution_local_rows
! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \retval local_cols ...
! *****************************************************************************
  FUNCTION dbcsr_distribution_local_cols(dist) RESULT (local_cols)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: local_cols

    local_cols = dist%d%local_cols
  END FUNCTION dbcsr_distribution_local_cols
  !
! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \param row ...
!> \param col ...
!> \retval processor ...
! *****************************************************************************
  PURE FUNCTION dbcsr_distribution_processor(dist, row, col)&
       RESULT (processor)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER, INTENT(IN)                      :: row, col
    INTEGER                                  :: processor

    INTEGER                                  :: c, r

    IF (ASSOCIATED (dist%d%row_map%low)) THEN ! instead of array_exists
       r = dist%d%row_map%low%data(row)
    ELSE
       r = row
    ENDIF
    IF (ASSOCIATED (dist%d%col_map%low)) THEN ! instead of array_exists
       c = dist%d%col_map%low%data(col)
    ELSE
       c = col
    ENDIF
    processor = dist%d%mp_env%mp%pgrid(dist%d%row_dist%low%data(r),&
         dist%d%col_dist%low%data(c))
  END FUNCTION dbcsr_distribution_processor

! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \retval thread_dist ...
! *****************************************************************************
  FUNCTION dbcsr_distribution_thread_dist(dist) RESULT (thread_dist)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: thread_dist

!   ---------------------------------------------------------------------------

    thread_dist = dist%d%thread_dist
  END FUNCTION dbcsr_distribution_thread_dist


! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \retval has_thread_dist ...
! *****************************************************************************
  PURE FUNCTION dbcsr_distribution_has_threads(dist) RESULT (has_thread_dist)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    LOGICAL                                  :: has_thread_dist

!   ---------------------------------------------------------------------------

    has_thread_dist = dist%d%has_thread_dist
  END FUNCTION dbcsr_distribution_has_threads

! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \retval num_threads ...
! *****************************************************************************
  PURE FUNCTION dbcsr_distribution_num_threads(dist) RESULT (num_threads)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER                                  :: num_threads

!   ---------------------------------------------------------------------------

    num_threads = dist%d%num_threads
  END FUNCTION dbcsr_distribution_num_threads

! *****************************************************************************
!> \brief Count max possible data storage for all processors used for matrix data
!> \param[in] matrix          Count data of this matrix
!> \retval data_size          Max data size used by matrix  for all processors
! *****************************************************************************
  FUNCTION dbcsr_distribution_data_max_size_allproc (matrix) RESULT (data_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER(KIND=int_8)                      :: data_size

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_distribution_data_max_size_allproc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, error_handler, row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: counts_cols, counts_rows
    INTEGER, DIMENSION(:), POINTER           :: col_blk_sizes, pcol_blk, &
                                                prow_blk, row_blk_sizes
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(dbcsr_error_type)                   :: error
    TYPE(dbcsr_mp_obj)                       :: mp_env

    CALL dbcsr_error_set(routineN, error_handler, error)

    dist = dbcsr_distribution(matrix)
    mp_env = dbcsr_distribution_mp(dist)
    ALLOCATE(counts_rows(0:dbcsr_mp_nprows(mp_env)-1),&
             counts_cols(0:dbcsr_mp_npcols(mp_env)-1))
    counts_rows = 0
    counts_cols = 0
    prow_blk => array_data (dbcsr_distribution_row_dist (dist))
    pcol_blk => array_data (dbcsr_distribution_col_dist (dist))
    row_blk_sizes => array_data (dbcsr_row_block_sizes (matrix))
    col_blk_sizes => array_data (dbcsr_col_block_sizes (matrix))

    DO row = 1, dbcsr_nblkrows_total(matrix)
       counts_rows(prow_blk(row)) = & 
            counts_rows(prow_blk(row)) + row_blk_sizes(row)
    ENDDO

    DO col = 1, dbcsr_nblkcols_total(matrix)
       counts_cols(pcol_blk(col)) = & 
            counts_cols(pcol_blk(col)) + col_blk_sizes(col)
    ENDDO

    data_size = MAXVAL(counts_rows)*MAXVAL(counts_cols)

    DEALLOCATE(counts_rows)
    DEALLOCATE(counts_cols)
    CALL dbcsr_error_stop(error_handler, error)
  END FUNCTION dbcsr_distribution_data_max_size_allproc


! *****************************************************************************
!> \brief Creates a distribution for threads
!> \param[in,out] dist   Add thread distribution to this distribution
!> \param[in] row_sizes  (optional) row block sizes
! *****************************************************************************
  SUBROUTINE dbcsr_distribution_make_threads(dist, row_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT), TARGET                  :: dist
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: row_sizes

    TYPE(dbcsr_distribution_obj), POINTER    :: dist_p

!   ---------------------------------------------------------------------------

    dist_p => dist
    !$ IF (.NOT. OMP_IN_PARALLEL ()) THEN
       !$OMP PARALLEL DEFAULT(NONE) &
       !$OMP          SHARED(dist_p,row_sizes)
       !$    CALL make_threads (dist_p, row_sizes=row_sizes)
       !$OMP END PARALLEL
    !$ ELSE
       CALL make_threads (dist_p, row_sizes=row_sizes)
       !$OMP BARRIER
    !$ ENDIF
  END SUBROUTINE dbcsr_distribution_make_threads

! *****************************************************************************
!> \brief Creates a distribution for threads
!> \param[in,out] dist   Add thread distribution to this distribution
!> \param[in] row_sizes  (optional) row block sizes
!> \par Presence of row_sizes
!>      When row_sizes is present then the thread distribution
!>      attempts to distribute rows to threads such that the sum of
!>      delegated row sizes is approximately matched for all rows.
!>
!>      When row_sizes is not present then a random distribution is chosen.
! *****************************************************************************
  SUBROUTINE make_threads(dist, row_sizes)
    TYPE(dbcsr_distribution_obj), POINTER    :: dist
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: row_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'make_threads', &
      routineP = moduleN//':'//routineN

    INTEGER :: block_size, block_size0, cur_block, group_size, i, last_row, &
      nlrows, nrows, nthreads, row, t, t_cnt
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: itemp1, itemp2, reorder, &
                                                sorted_row_sizes
    INTEGER, DIMENSION(:), POINTER           :: lrows, td
    LOGICAL                                  :: assigned, found, heap_error
    REAL(kind=sp)                            :: load_fraction, rn, soft_thr

!$  TYPE(dbcsr_error_type)                   :: error
    TYPE(heap_t)                             :: t_heap

!   ---------------------------------------------------------------------------

    nthreads = 1
!$  nthreads = OMP_GET_NUM_THREADS () ;

!$  CALL dbcsr_assert (dist%d%num_threads, "EQ", nthreads,&
!$                     dbcsr_fatal_level, dbcsr_internal_error, routineN,&
!$                     "Thread number has changed", __LINE__, error=error)
    nrows = dbcsr_distribution_nrows (dist)
    nlrows = dbcsr_distribution_nlocal_rows (dist)
    lrows => array_data (dbcsr_distribution_local_rows (dist))

    !$OMP BARRIER
    !$OMP MASTER


    load_fraction = REAL(comm_thread_load)/100.0
    IF(nthreads==1) load_fraction = 1.0

    IF (.NOT. dist%d%has_thread_dist) THEN
       dist%d%num_threads = nthreads
       group_size = 0 ; cur_block = 0

       ALLOCATE (td(nrows))
       dist%d%has_thread_dist = .TRUE.
       CALL array_new (dist%d%thread_dist, td, gift=.TRUE.)
       td => array_data (dist%d%thread_dist)

       IF (PRESENT (row_sizes)) THEN
          ! The goal is to distribute rows to threads as equally as
          ! possible. The row sizes are first sorted. Each group of
          ! equally sized rows (group_size rows of size cur_block) is
          ! distributed to threads (keeping consecutive rows
          ! together). The group is divided into equally-sized blocks
          ! (block_size0, block_size).  Leftover rows (those that can
          ! not be equally distributed to threads) are then assigned
          ! to threads so that each thread's commulative load attempts
          ! to be equal. This distribution is achieved using a heap.
          !
          ! The heap is used to distribute "leftover"rows to threads.
          ! Lefotver rows are those of the same size that can not be
          ! evenly distributed among all threads.
          CALL heap_new (t_heap, nthreads-1, heap_error)
          ! We do not want thread 0 to be in the heap.
          ALLOCATE(itemp1(1:nthreads-1))
          ALLOCATE(itemp2(1:nthreads-1))
          DO i=1,nthreads-1
             itemp1(i)=i
             itemp2(i)=0
          ENDDO
          CALL heap_fill (t_heap,itemp1,itemp2,heap_error)
          DEALLOCATE(itemp1,itemp2)
          ALLOCATE (sorted_row_sizes (nrows))
          ALLOCATE (reorder (nrows))
          sorted_row_sizes(:) = row_sizes(:)
          CALL sort (sorted_row_sizes, nrows, reorder)

          row = 1
          DO WHILE ( row .LE. nrows)
             cur_block = sorted_row_sizes(nrows-row+1)
             assigned = .FALSE.
             group_size = 0

             last_row = nrows-row+1
             DO i = last_row, 1, -1
                IF ( cur_block == sorted_row_sizes(i) ) THEN
                   group_size = group_size + 1
                   row = row + 1
                ELSE
                   EXIT
                END IF
             END DO

             soft_thr = load_fraction + nthreads - 1
             block_size0 = INT(load_fraction*(group_size/soft_thr))
             block_size = INT(group_size/soft_thr)

             !blocks for master thread
             IF (block_size0 > 0) &
                td(reorder(last_row:last_row-block_size0+1:-1)) = 0

             !Other threads
             IF (block_size > 0) THEN
                DO t=1, nthreads-1
                   td(reorder(last_row-block_size0-(t-1)*block_size:&
                        last_row-block_size0-(t)*block_size+1:-1)) = t
                END DO
             END IF

             !Leftover bocks
             DO i=last_row-block_size0-(nthreads-1)*block_size, last_row+1-group_size, -1
                CALL heap_get_first (t_heap, t, t_cnt, found,heap_error)
                t_cnt = t_cnt + cur_block
                CALL heap_reset_first (t_heap, t_cnt, heap_error)
                td(reorder(i)) = t
             END DO

          END DO
          CALL heap_release (t_heap, heap_error)
          DEALLOCATE (sorted_row_sizes)
          DEALLOCATE (reorder)
       ELSE
          DO t = 1, nrows
             IF (.FALSE.) THEN
                td(t) = MOD(t-1, nthreads)
             ELSE
                CALL RANDOM_NUMBER (rn)
                ! Makes sure the numbers are in the proper integer range.
                td(t) = MOD (INT (rn*REAL(nthreads)), nthreads)
             ENDIF
          END DO
       ENDIF
    ENDIF
    !$OMP END MASTER
  END SUBROUTINE make_threads

!> \brief Removes the thread distribution from a distribution
!> \param dist ...
! *****************************************************************************
  SUBROUTINE dbcsr_distribution_no_threads(dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist

!$OMP MASTER
    CALL array_release (dist%d%thread_dist)
    dist%d%has_thread_dist = .FALSE.
!$OMP END MASTER
  END SUBROUTINE dbcsr_distribution_no_threads

END MODULE dbcsr_dist_methods

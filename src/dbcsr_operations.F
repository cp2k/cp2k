!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   [DESCRIPTION OF THE MODULE]
!> \author  [AUTHOR]
!> \date    [DATE]
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - [ADD MODIFICATION HISTORY or CREATION]
! *****************************************************************************
MODULE dbcsr_operations
  USE btree_I8_k_I4_v,                 ONLY: btree_I8_I4 => btree,&
                                             btree_add,&
                                             btree_delete,&
                                             btree_find,&
                                             btree_get_entries,&
                                             btree_new
  USE dbcsr_methods,                   ONLY: dbcsr_distribution_mp,&
                                             dbcsr_mp_group,&
                                             dbcsr_mp_mypcol,&
                                             dbcsr_mp_myproc,&
                                             dbcsr_mp_myprow,&
                                             dbcsr_mp_nprows,&
                                             dbcsr_mp_numproc,&
                                             dbcsr_mp_pgrid
  USE dbcsr_types,                     ONLY: dbcsr_array_type,&
                                             dbcsr_meta_size,&
                                             dbcsr_mp_obj,&
                                             dbcsr_type,&
                                             dbcsr_work_type
  USE dbcsr_util,                      ONLY: &
       add_work_coordinate, dbcsr_create, dbcsr_destroy_array, &
       dbcsr_finalize, dbcsr_pack_meta, dbcsr_print, dbcsr_repoint_index, &
       dbcsr_set_debug, dbcsr_unpack_i8_2i4, ensure_array_size
  USE kinds,                           ONLY: dp,&
                                             int_4,&
                                             int_8
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_isendrecv,&
                                             mp_rma_get,&
                                             mp_wait,&
                                             mp_window_create,&
                                             mp_window_destroy,&
                                             mp_window_fence
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_operations'

  PUBLIC :: dbcsr_multiply_sets, dbcsr_multiply_sets_sqrt, dbcsr_nn_mult

  PRIVATE

#define DBG IF (dbg) WRITE(*,*)routineN//" ",
#define DBGV IF (bcsr_verbose) WRITE(*,*)routineN//" ",
#define DBGI IF (info) WRITE(*,*)routineN//" ",
#define DEBUG_HEADER        LOGICAL :: dbg, info
#define DEBUG_BODY        dbg = .FALSE. ; CALL dbcsr_set_debug(dbg, info=info)
#define DEBUG_BODYY       dbg = .TRUE. ; CALL dbcsr_set_debug(dbg, info=info)

#define temp_transpose(v, r, c) RESHAPE(TRANSPOSE(RESHAPE(v,(/r,c/))),(/r*c/))


CONTAINS

! *****************************************************************************
!> \brief Multiplies two DBCSR matrices, using the current local data.
!> \param[in] left_set, right_set  left and right DBCSR matrices
!> \param[out] product             DBCSR product workmatrix
!> \param[in,out] error            cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_sets (left_set, right_set, product_matrix,&
       product_wm, error)
    TYPE(dbcsr_array_type), INTENT(INOUT)    :: left_set, right_set
    TYPE(dbcsr_type), INTENT(INOUT)          :: product_matrix
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: product_wm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_sets', &
      routineP = moduleN//':'//routineN

    INTEGER :: entries, flop, flop_total, i, left_nprows, left_prow, &
      left_prow_i, left_rep, left_reps, left_src_p_calc, left_src_p_comm, &
      mp_group, mypcol, myproc, myprow, numproc, rep, right_nprows, &
      right_prow, right_prow_i, right_rep, right_reps, right_src_p_calc, &
      right_src_p_comm, timing_handle
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: keys
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: my_sizes
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: all_sizes
    INTEGER, DIMENSION(:, :), POINTER        :: left_pgrid, product_pgrid, &
                                                right_pgrid
    INTEGER, DIMENSION(dbcsr_meta_size)      :: new_meta
    REAL(KIND=dp)                            :: trun, tstart, tstop
    TYPE(btree_I8_I4)                        :: windex
    TYPE(dbcsr_array_type), POINTER          :: lbuffer_set_calc, &
                                                lbuffer_set_comm, &
                                                rbuffer_set_calc, &
                                                rbuffer_set_comm
    TYPE(dbcsr_array_type), TARGET           :: lbuffer_set1, lbuffer_set2, &
                                                rbuffer_set1, rbuffer_set2
    TYPE(dbcsr_mp_obj)                       :: mp_obj

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset (routineN, timing_handle)
    CALL btree_new (windex)
    tstart = 0.0d0 ; tstop = 0.0d0 ; trun = 0.0d0
    left_reps = SIZE (left_set%mats)
    right_reps = SIZE (right_set%mats)
    mp_obj = dbcsr_distribution_mp (product_matrix%dist)
    numproc = dbcsr_mp_numproc (mp_obj)
    myproc = dbcsr_mp_myproc (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    product_pgrid => dbcsr_mp_pgrid (mp_obj)
    left_pgrid => dbcsr_mp_pgrid (dbcsr_distribution_mp (&
         left_set%rdist%main))
    left_nprows = dbcsr_mp_nprows (dbcsr_distribution_mp (&
         left_set%rdist%main))
    right_pgrid => dbcsr_mp_pgrid (dbcsr_distribution_mp (&
         right_set%mats(1)%dist))
    DBG 'product pgrid shape', SHAPE (dbcsr_mp_pgrid (mp_obj))
    DBG 'left pgrid shape', SHAPE (left_pgrid),';',left_nprows,left_reps
    DBG 'right pgrid shape', SHAPE (right_pgrid)
    ALLOCATE (lbuffer_set1%mats(0:left_nprows-1))
    ALLOCATE (lbuffer_set2%mats(0:left_nprows-1))
    ALLOCATE (rbuffer_set1%mats(right_reps))
    ALLOCATE (rbuffer_set2%mats(right_reps))
    ALLOCATE (left_set%data_windows(left_reps))
    ALLOCATE (left_set%index_windows(left_reps))
    ALLOCATE (right_set%data_windows(right_reps))
    ALLOCATE (right_set%index_windows(right_reps))
    ALLOCATE (my_sizes(4, MAX(left_reps, right_reps)))
    my_sizes(:,:) = 0
    DO rep = 1, left_reps
       my_sizes(1, rep) = SIZE (left_set%mats(rep)%data)
       my_sizes(2, rep) = SIZE (left_set%mats(rep)%index)
    ENDDO
    DO rep = 1, right_reps
       my_sizes(3, rep) = SIZE (right_set%mats(rep)%data)
       my_sizes(4, rep) = SIZE (right_set%mats(rep)%index)
    ENDDO
    ALLOCATE (all_sizes(4, LBOUND(my_sizes,2):UBOUND(my_sizes,2), 0:numproc-1))
    DO i = 1, 4
       DBG 'my sizes', my_sizes(i,:)
    ENDDO
    CALL mp_allgather(my_sizes, all_sizes, mp_group)
    DBG 'data sizes', all_sizes
    ! Setup the buffer matrices
    DO rep = 1, left_reps
       CALL mp_window_create (left_set%data_windows(rep),&
            left_set%mats(rep)%data, gid=mp_group)
       CALL mp_window_create (left_set%index_windows(rep),&
            left_set%mats(rep)%index, gid=mp_group)
       !CALL mp_window_fence (left_set%data_windows(rep))
       !CALL mp_window_fence (left_set%index_windows(rep))
    ENDDO
    DO left_prow = 0, left_nprows-1
       DBG 'setting up left buffer for repetition', left_prow
       CALL dbcsr_pack_meta(left_set%mats(1), new_meta, error)
       DBG 'meta for right matrix',new_meta
       CALL dbcsr_create (lbuffer_set1%mats(left_prow),&
            'lbufferl1'//left_set%mats(1)%name,&
            left_set%mats(1)%dist2d,&
            'N', new_meta, 0, 0, special=.FALSE., error=error)
       lbuffer_set1%mats(left_prow)%row_blk_size =>&
            left_set%mats(1)%row_blk_size
       lbuffer_set1%mats(left_prow)%col_blk_size =>&
            left_set%mats(1)%col_blk_size
       CALL ensure_array_size (lbuffer_set1%mats(left_prow)%data,&
            ub=MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
       CALL ensure_array_size (lbuffer_set1%mats(left_prow)%index,&
            ub=MAXVAL (all_sizes(2, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
       !
       CALL dbcsr_create (lbuffer_set2%mats(left_prow),&
            'lbuffer2'//left_set%mats(1)%name,&
            left_set%mats(1)%dist2d,&
            'N', new_meta, 0, 0, special=.FALSE., error=error)
       lbuffer_set2%mats(left_prow)%row_blk_size =>&
            left_set%mats(1)%row_blk_size
       lbuffer_set2%mats(left_prow)%col_blk_size =>&
            left_set%mats(1)%col_blk_size
       CALL ensure_array_size (lbuffer_set2%mats(left_prow)%data,&
            ub=MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
       CALL ensure_array_size (lbuffer_set2%mats(left_prow)%index,&
            ub=MAXVAL (all_sizes(2, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
    ENDDO
    lbuffer_set_calc => lbuffer_set1
    lbuffer_set_comm => lbuffer_set2
    DO rep = 1, right_reps
       DBG 'setting up right buffer for repetition', rep
       CALL dbcsr_pack_meta(right_set%mats(rep), new_meta, error)
       CALL dbcsr_create (rbuffer_set1%mats(rep),&
            'rbuffer1'//right_set%mats(rep)%name, right_set%mats(rep)%dist2d,&
            'N', new_meta, 0, 0, special=.FALSE., error=error)
       rbuffer_set1%mats(rep)%row_blk_size => right_set%mats(rep)%row_blk_size
       rbuffer_set1%mats(rep)%col_blk_size => right_set%mats(rep)%col_blk_size
       CALL ensure_array_size (rbuffer_set1%mats(rep)%data,&
            ub=MAXVAL (all_sizes(3, rep, :)), nocopy=.TRUE., special=.FALSE.,&
            error=error)
       CALL ensure_array_size (rbuffer_set1%mats(rep)%index,&
            ub=MAXVAL (all_sizes(4, rep, :)), nocopy=.TRUE., special=.FALSE.,&
            error=error)
       !
       CALL dbcsr_create (rbuffer_set2%mats(rep),&
            'rbuffer2'//right_set%mats(rep)%name, right_set%mats(rep)%dist2d,&
            'N', new_meta, 0, 0, special=.FALSE., error=error)
       rbuffer_set2%mats(rep)%row_blk_size => right_set%mats(rep)%row_blk_size
       rbuffer_set2%mats(rep)%col_blk_size => right_set%mats(rep)%col_blk_size
       CALL ensure_array_size (rbuffer_set2%mats(rep)%data,&
            ub=MAXVAL (all_sizes(3, rep, :)), nocopy=.TRUE., special=.FALSE.,&
            error=error)
       CALL ensure_array_size (rbuffer_set2%mats(rep)%index,&
            ub=MAXVAL (all_sizes(4, rep, :)), nocopy=.TRUE., special=.FALSE.,&
            error=error)
       !
       DBG 'RIGHT data win size', SIZE (right_set%mats(rep)%data)
       DBG 'RIGHT indx win size', SIZE (right_set%mats(rep)%index)
       CALL mp_window_create (right_set%data_windows(rep),&
            right_set%mats(rep)%data, gid=mp_group)
       CALL mp_window_create (right_set%index_windows(rep),&
            right_set%mats(rep)%index, gid=mp_group)
       !CALL mp_window_fence (right_set%data_windows(rep))
       !CALL mp_window_fence (right_set%index_windows(rep))
    ENDDO
    rbuffer_set_calc => rbuffer_set1
    rbuffer_set_comm => rbuffer_set2
    !
    right_nprows = dbcsr_mp_nprows (dbcsr_distribution_mp (&
         right_set%mats(1)%dist))
    right_src_p_calc = -7
    left_src_p_calc = -7
    k_index: DO right_prow_i = 0, right_nprows
       right_prow = MOD (right_prow_i + myprow, right_nprows)
       DBG 'RIGHT_PROW from/to', right_prow_i, right_prow, myprow, right_nprows
       ! No fetches in the last loop
       right_rep = 1
       CALL mp_window_fence (right_set%data_windows(right_rep),&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
       CALL mp_window_fence (right_set%index_windows(right_rep),&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
       xfer_case_right: IF (right_prow_i .LT. right_nprows) THEN
          right_src_p_comm = right_pgrid (right_prow, mypcol)
          IF (dbg) WRITE(*,'(1X,A,"(",I3,",",I3,")")')&
               routineN//' Getting RIGHT matrix from pgrid coor ',&
               right_prow, mypcol
          IF (dbg) WRITE(*,'(1X,A,I4,A,"(",I3,",",I3,") sizes",I7,I7)')&
               routineN//' Getting RIGHT matrix from',right_src_p_comm,&
               ' pgrid coor ',right_prow, mypcol,&
               all_sizes(3:4, right_rep, right_src_p_comm)
          CALL mp_rma_get (right_set%data_windows(right_rep),&
               right_src_p_comm, 0, all_sizes(3, right_rep, right_src_p_comm),&
               rbuffer_set_comm%mats(right_rep)%data)
          CALL mp_rma_get (right_set%index_windows(right_rep),&
               right_src_p_comm, 0, all_sizes(4, right_rep, right_src_p_comm),&
               rbuffer_set_comm%mats(right_rep)%index)
       ENDIF xfer_case_right
       calc_case_right: IF (right_prow_i .GT. 0) THEN
          !right_src_p_calc = right_pgrid (right_prow-1, mypcol)
          right_rep = 1
          IF (all_sizes(4, right_rep, right_src_p_calc) .EQ. 0) THEN
             DBG 'ZERO SIZED INDEX'
          ENDIF
          CALL dbcsr_repoint_index (rbuffer_set_calc%mats(right_rep),&
               error=error)
       ENDIF calc_case_right
       right_src_p_calc = right_src_p_comm
       DO left_rep = 1, left_reps
          CALL mp_window_fence (left_set%data_windows(left_rep),&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
          CALL mp_window_fence (left_set%index_windows(left_rep),&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
       ENDDO
       left_rep = myprow+1 ! Because reps are offset 1 from myprow
       left_row_aggregation: DO left_prow_i = 0, left_nprows-1
          left_prow = MOD (left_prow_i + myprow, left_nprows)
          xfer_case_left: IF (right_prow_i .LT. right_nprows) THEN
             IF (dbg) WRITE(*,'(1X,A,I3,A,"(",I3,",",I3,")")')&
                  routineN//' Getting LEFT matrix repetition',left_rep,&
                  ' pgrid coor ',left_prow, right_prow
             left_src_p_comm = left_pgrid (left_prow, right_prow)
             IF (dbg) WRITE(*,'(1X,A,I4,A,I3,A,"(",I3,",",I3,")")')&
                  routineN//' Getting LEFT matrix from',left_src_p_comm,&
                  ' repetition',left_rep,&
                  ' pgrid coor ',left_prow, right_prow
             CALL mp_rma_get (left_set%data_windows(left_rep),&
                  left_src_p_comm, 0,&
                  all_sizes(1, left_rep, left_src_p_comm),&
                  lbuffer_set_comm%mats(left_prow)%data)
             CALL mp_rma_get (left_set%index_windows(left_rep),&
                  left_src_p_comm, 0,&
                  all_sizes(2, left_rep, left_src_p_comm),&
                  lbuffer_set_comm%mats(left_prow)%index)
          ENDIF xfer_case_left
          calc_case_left: IF (right_prow_i .GT. 0) THEN
             !left_src_p_calc = left_pgrid (left_prow, right_prow-1)
             CALL dbcsr_repoint_index (lbuffer_set_calc%mats(left_prow),&
                  error=error)
             CALL dbcsr_print (lbuffer_set_calc%mats(left_prow), error)
             CALL dbcsr_print (rbuffer_set_calc%mats(1), error)
             tstart = m_walltime ()
             CALL dbcsr_nn_mult(lbuffer_set_calc%mats(left_prow),&
                  rbuffer_set_calc%mats(1),&
                  product_matrix, product_wm, error, windex, flop=flop)
             flop_total = flop_total + flop
             tstop = m_walltime ()
             DBG 'Segment Local Multiplication time',tstop-tstart,'and',&
                  (REAL(flop, dp)/1000000.0d0)/(tstop-tstart),'MFLOPS'
             trun = trun + (tstop - tstart)
          ENDIF calc_case_left
       ENDDO left_row_aggregation
       left_src_p_calc = left_src_p_comm
       CALL dbcsr_switch_sets(lbuffer_set_calc, lbuffer_set_comm)
       CALL dbcsr_switch_sets(rbuffer_set_calc, rbuffer_set_comm)
    ENDDO k_index
    WRITE(*,'(1X,A,F9.3,A,F9.3,A,I9)')&
         routineN//' Total Local Multiplication time',&
         trun,' and',&
         (REAL(flop_total, dp)/1000000.0d0)/trun,' MFLOPS',flop_total
    DO rep = 1, left_reps
       CALL mp_window_destroy (left_set%data_windows(rep))
       CALL mp_window_destroy (left_set%index_windows(rep))
    ENDDO
    DO rep = 1, right_reps
       CALL mp_window_destroy (right_set%data_windows(rep))
       CALL mp_window_destroy (right_set%index_windows(rep))
    ENDDO
    !
    entries = btree_get_entries (windex)
    CALL ensure_array_size (product_wm%blk_p, ub=entries, error=error)
    ALLOCATE (keys (entries))
    CALL btree_delete (windex, keys, product_wm%blk_p)
    ! Recreate the row, column arrays from the B-tree.
    !lmask = 0 ; umask = 0
    !DO i = 0, 31
    !   lmask = IOR (lmask, IBSET (lmask, i))
    !   umask = IOR (umask, IBSET (umask, i+32))
    !ENDDO
    !lmask8 = lmask
    CALL ensure_array_size(product_wm%row_i, ub=entries, error=error)
    CALL ensure_array_size(product_wm%col_i, ub=entries, error=error)
    !DO i = 1, entries
    !   product_wm%row_i(i) = ISHFT (keys(i), -32)
    !   product_wm%col_i(i) = IAND (keys(i), lmask8)
    !ENDDO
    CALL dbcsr_unpack_i8_2i4 (keys, product_wm%row_i, product_wm%col_i)
    DEALLOCATE (keys)
    CALL dbcsr_finalize(product_wm, product_matrix, .TRUE., .FALSE., error)
    !
    !DO left_prow = 0, left_nprows - 1
    !   CALL dbcsr_destroy (lbuffer_set1%mats(left_prow), error)
    !   CALL dbcsr_destroy (lbuffer_set2%mats(left_prow), error)
    !ENDDO
    !DO right_prow = 0, right_nprows - 1
    !   CALL dbcsr_destroy (rbuffer_set1%mats(right_prow), error)
    !   CALL dbcsr_destroy (rbuffer_set2%mats(right_prow), error)
    !ENDDO
    CALL dbcsr_destroy_array (rbuffer_set1, error)
    CALL dbcsr_destroy_array (rbuffer_set2, error)
    CALL dbcsr_destroy_array (lbuffer_set1, error)
    CALL dbcsr_destroy_array (lbuffer_set2, error)
    DEALLOCATE (left_set%data_windows)
    DEALLOCATE (left_set%index_windows)
    DEALLOCATE (right_set%data_windows)
    DEALLOCATE (right_set%index_windows)
    DEALLOCATE (my_sizes)
    CALL timestop (timing_handle)
  END SUBROUTINE dbcsr_multiply_sets


! *****************************************************************************
!> \brief Switches pointers between two matrix sets
!> \param[in,out] set1p, set2p
! *****************************************************************************
  SUBROUTINE dbcsr_switch_sets (set1p, set2p)
    TYPE(dbcsr_array_type), POINTER          :: set1p, set2p

    TYPE(dbcsr_array_type), POINTER          :: tmp_set

!   ---------------------------------------------------------------------------

    tmp_set => set1p
    set1p => set2p
    set2p => tmp_set
  END SUBROUTINE dbcsr_switch_sets


! *****************************************************************************
!> \brief Makes an MPI tag
!> \param[in,out] set1p, set2p
! *****************************************************************************
  SUBROUTINE make_tag (tag, to, from, seq)
    INTEGER, INTENT(OUT)                     :: tag
    INTEGER, INTENT(IN)                      :: to, from, seq

    INTEGER, PARAMETER                       :: s = 8

!   ---------------------------------------------------------------------------

    tag = ISHFT (to, s) + from + ISHFT (seq, 2*s)
  END SUBROUTINE make_tag


! *****************************************************************************
!> \brief Multiplies two DBCSR matrices, passing data only among processor
!>        rows
!> \param[in] left_set, right_set  left and right DBCSR matrices
!> \param[out] product             DBCSR product workmatrix
!> \param[in,out] error            cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_sets_sqrt (left_set, right_set, product_matrix,&
       product_wm, error)
    TYPE(dbcsr_array_type), INTENT(INOUT)    :: left_set, right_set
    TYPE(dbcsr_type), INTENT(INOUT)          :: product_matrix
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: product_wm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_sets_sqrt', &
      routineP = moduleN//':'//routineN

    INTEGER :: entries, flop, flop_total, i, left_nprows, left_prow, &
      left_prow_i, left_rep, left_reps, left_src_p_calc, left_src_p_comm, &
      mp_group, mypcol, myproc, myprow, numproc, rep, right_dst_p_comm, &
      right_nprows, right_prow, right_prow_i, right_rep, right_reps, &
      right_src_p_calc, right_src_p_comm, timing_handle
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: keys
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: right_data_rr, right_data_sr, &
                                                right_index_rr, right_index_sr
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: my_sizes
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: all_sizes
    INTEGER, DIMENSION(:, :), POINTER        :: left_pgrid, product_pgrid, &
                                                right_pgrid
    INTEGER, DIMENSION(dbcsr_meta_size)      :: new_meta
    REAL(KIND=dp)                            :: trun, tstart, tstop
    TYPE(btree_I8_I4)                        :: windex
    TYPE(dbcsr_array_type), POINTER          :: lbuffer_set_calc, &
                                                lbuffer_set_comm, &
                                                rbuffer_set_calc, &
                                                rbuffer_set_comm
    TYPE(dbcsr_array_type), TARGET           :: lbuffer_set1, lbuffer_set2, &
                                                rbuffer_set1, rbuffer_set2
    TYPE(dbcsr_mp_obj)                       :: mp_obj

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset (routineN, timing_handle)
    CALL btree_new (windex)
    !
    ! An overview of the algorithm: a processor multiplies the local
    ! left matrices (i.e., the left_set) and the local right matrix
    ! the (right_set). The data must then be sent to the appropriate
    ! processor, which is in the same processor column. The
    ! destination processor then sums up the matrices.
    !
    ! Instead of prefetching all local left matrices, they could be
    ! fetched in turn (through the process rows).
    !
    tstart = 0.0d0 ; tstop = 0.0d0 ; trun = 0.0d0
    tstart = m_walltime ()
    left_reps = SIZE (left_set%mats)
    right_reps = SIZE (right_set%mats)
    mp_obj = dbcsr_distribution_mp (product_matrix%dist)
    numproc = dbcsr_mp_numproc (mp_obj)
    myproc = dbcsr_mp_myproc (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    product_pgrid => dbcsr_mp_pgrid (mp_obj)
    left_pgrid => dbcsr_mp_pgrid (dbcsr_distribution_mp (&
         left_set%rdist%main))
    left_nprows = dbcsr_mp_nprows (dbcsr_distribution_mp (&
         left_set%rdist%main))
    right_pgrid => dbcsr_mp_pgrid (dbcsr_distribution_mp (&
         right_set%mats(1)%dist))
    DBG 'product pgrid shape', SHAPE (dbcsr_mp_pgrid (mp_obj))
    DBG 'left pgrid shape', SHAPE (left_pgrid),';',left_nprows,left_reps
    DBG 'right pgrid shape', SHAPE (right_pgrid)
    ALLOCATE (lbuffer_set1%mats(0:left_nprows-1))
    ALLOCATE (lbuffer_set2%mats(0:left_nprows-1))
    ALLOCATE (rbuffer_set1%mats(right_reps))
    ALLOCATE (rbuffer_set2%mats(right_reps))
    ALLOCATE (right_data_sr(right_reps))
    ALLOCATE (right_data_rr(right_reps))
    ALLOCATE (right_index_sr(right_reps))
    ALLOCATE (right_index_rr(right_reps))
    ALLOCATE (left_set%data_windows(left_reps))
    ALLOCATE (left_set%index_windows(left_reps))
    ALLOCATE (right_set%data_windows(right_reps))
    ALLOCATE (right_set%index_windows(right_reps))
    ALLOCATE (my_sizes(4, MAX(left_reps, right_reps)))
    my_sizes(:,:) = 0
    DO rep = 1, left_reps
       my_sizes(1, rep) = SIZE (left_set%mats(rep)%data)
       my_sizes(2, rep) = SIZE (left_set%mats(rep)%index)
    ENDDO
    DO rep = 1, right_reps
       my_sizes(3, rep) = SIZE (right_set%mats(rep)%data)
       my_sizes(4, rep) = SIZE (right_set%mats(rep)%index)
    ENDDO
    ALLOCATE (all_sizes(4, LBOUND(my_sizes,2):UBOUND(my_sizes,2), 0:numproc-1))
    DO i = 1, 4
       DBG 'my sizes', my_sizes(i,:)
    ENDDO
    CALL mp_allgather(my_sizes, all_sizes, mp_group)
    DBG 'data sizes', all_sizes
    ! Setup the buffer matrices
    DO rep = 1, left_reps
       CALL mp_window_create (left_set%data_windows(rep),&
            left_set%mats(rep)%data, gid=mp_group)
       CALL mp_window_create (left_set%index_windows(rep),&
            left_set%mats(rep)%index, gid=mp_group)
       !CALL mp_window_fence (left_set%data_windows(rep))
       !CALL mp_window_fence (left_set%index_windows(rep))
    ENDDO
    DO left_prow = 0, left_nprows-1
       DBG 'setting up left buffer for repetition', left_prow
       CALL dbcsr_pack_meta(left_set%mats(1), new_meta, error)
       DBG 'meta for right matrix',new_meta
       CALL dbcsr_create (lbuffer_set1%mats(left_prow),&
            'lbufferl1'//left_set%mats(1)%name,&
            left_set%mats(1)%dist2d,&
            'N', new_meta, 0, 0, special=.FALSE., error=error)
       lbuffer_set1%mats(left_prow)%row_blk_size =>&
            left_set%mats(1)%row_blk_size
       lbuffer_set1%mats(left_prow)%col_blk_size =>&
            left_set%mats(1)%col_blk_size
       CALL ensure_array_size (lbuffer_set1%mats(left_prow)%data,&
            ub=MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
       CALL ensure_array_size (lbuffer_set1%mats(left_prow)%index,&
            ub=MAXVAL (all_sizes(2, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
       !
       CALL dbcsr_create (lbuffer_set2%mats(left_prow),&
            'lbuffer2'//left_set%mats(1)%name,&
            left_set%mats(1)%dist2d,&
            'N', new_meta, 0, 0, special=.FALSE., error=error)
       lbuffer_set2%mats(left_prow)%row_blk_size =>&
            left_set%mats(1)%row_blk_size
       lbuffer_set2%mats(left_prow)%col_blk_size =>&
            left_set%mats(1)%col_blk_size
       CALL ensure_array_size (lbuffer_set2%mats(left_prow)%data,&
            ub=MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
       CALL ensure_array_size (lbuffer_set2%mats(left_prow)%index,&
            ub=MAXVAL (all_sizes(2, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
    ENDDO
    lbuffer_set_calc => lbuffer_set1
    lbuffer_set_comm => lbuffer_set2
    DO rep = 1, right_reps
       DBG 'setting up right buffer for repetition', rep
       CALL dbcsr_pack_meta(right_set%mats(rep), new_meta, error)
       CALL dbcsr_create (rbuffer_set1%mats(rep),&
            'rbuffer1'//right_set%mats(rep)%name, right_set%mats(rep)%dist2d,&
            'N', new_meta, 0, 0, special=.FALSE., error=error)
       rbuffer_set1%mats(rep)%row_blk_size => right_set%mats(rep)%row_blk_size
       rbuffer_set1%mats(rep)%col_blk_size => right_set%mats(rep)%col_blk_size
       CALL ensure_array_size (rbuffer_set1%mats(rep)%data,&
            ub=MAXVAL (all_sizes(3, rep, :)), nocopy=.TRUE., special=.FALSE.,&
            error=error)
       CALL ensure_array_size (rbuffer_set1%mats(rep)%index,&
            ub=MAXVAL (all_sizes(4, rep, :)), nocopy=.TRUE., special=.FALSE.,&
            error=error)
       !
       CALL dbcsr_create (rbuffer_set2%mats(rep),&
            'rbuffer2'//right_set%mats(rep)%name, right_set%mats(rep)%dist2d,&
            'N', new_meta, 0, 0, special=.FALSE., error=error)
       rbuffer_set2%mats(rep)%row_blk_size => right_set%mats(rep)%row_blk_size
       rbuffer_set2%mats(rep)%col_blk_size => right_set%mats(rep)%col_blk_size
       CALL ensure_array_size (rbuffer_set2%mats(rep)%data,&
            ub=MAXVAL (all_sizes(3, rep, :)), nocopy=.TRUE., special=.FALSE.,&
            error=error)
       CALL ensure_array_size (rbuffer_set2%mats(rep)%index,&
            ub=MAXVAL (all_sizes(4, rep, :)), nocopy=.TRUE., special=.FALSE.,&
            error=error)
       !
       DBG 'RIGHT data win size', SIZE (right_set%mats(rep)%data)
       DBG 'RIGHT indx win size', SIZE (right_set%mats(rep)%index)
       CALL mp_window_create (right_set%data_windows(rep),&
            right_set%mats(rep)%data, gid=mp_group)
       CALL mp_window_create (right_set%index_windows(rep),&
            right_set%mats(rep)%index, gid=mp_group)
       !CALL mp_window_fence (right_set%data_windows(rep))
       !CALL mp_window_fence (right_set%index_windows(rep))
    ENDDO
    rbuffer_set_calc => rbuffer_set1
    rbuffer_set_comm => rbuffer_set2
    !
    right_nprows = dbcsr_mp_nprows (dbcsr_distribution_mp (&
         right_set%mats(1)%dist))
    right_src_p_calc = -7
    left_src_p_calc = -7
    ! RMA fetch the right matrix
    right_prow = MOD (myprow, right_nprows)
    DBG 'RIGHT_PROW from/to', right_prow_i, right_prow, myprow, right_nprows
    right_rep = 1
    CALL mp_window_fence (right_set%data_windows(right_rep),&
         nothing_stored=.TRUE., no_put_follows=.TRUE.)
    CALL mp_window_fence (right_set%index_windows(right_rep),&
         nothing_stored=.TRUE., no_put_follows=.TRUE.)
    right_src_p_comm = right_pgrid (right_prow, mypcol)
    IF (dbg) WRITE(*,'(1X,A,"(",I3,",",I3,")")')&
         routineN//' Getting RIGHT matrix from pgrid coor ',&
         right_prow, mypcol
    IF (dbg) WRITE(*,'(1X,A,I4,A,"(",I3,",",I3,") sizes",I7,I7)')&
         routineN//' Getting RIGHT matrix from',right_src_p_comm,&
         ' pgrid coor ',right_prow, mypcol,&
         all_sizes(3:4, right_rep, right_src_p_comm)
    CALL mp_rma_get (right_set%data_windows(right_rep),&
         right_src_p_comm, 0, all_sizes(3, right_rep, right_src_p_comm),&
         rbuffer_set_comm%mats(right_rep)%data)
    CALL mp_rma_get (right_set%index_windows(right_rep),&
         right_src_p_comm, 0, all_sizes(4, right_rep, right_src_p_comm),&
         rbuffer_set_comm%mats(right_rep)%index)
    CALL mp_window_fence (right_set%data_windows(right_rep),&
         nothing_stored=.TRUE., no_put_follows=.TRUE.)
    CALL mp_window_fence (right_set%index_windows(right_rep),&
         nothing_stored=.TRUE., no_put_follows=.TRUE.)
    CALL dbcsr_repoint_index (rbuffer_set_comm%mats(right_rep),&
         error=error)
    right_src_p_calc = right_src_p_comm
    ! RMA fetch the left matrix
    DO left_rep = 1, left_reps
       CALL mp_window_fence (left_set%data_windows(left_rep),&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
       CALL mp_window_fence (left_set%index_windows(left_rep),&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
    ENDDO
    left_rep = myprow+1 ! Because reps are offset 1 from myprow
    left_row_aggregation: DO left_prow_i = 0, left_nprows-1
       left_prow = MOD (left_prow_i + myprow, left_nprows)
       IF (dbg) WRITE(*,'(1X,A,I3,A,"(",I3,",",I3,")")')&
            routineN//' Getting LEFT matrix repetition',left_rep,&
            ' pgrid coor ',left_prow, right_prow
       left_src_p_comm = left_pgrid (left_prow, right_prow)
       IF (dbg) WRITE(*,'(1X,A,I4,A,I3,A,"(",I3,",",I3,")")')&
            routineN//' Getting LEFT matrix from',left_src_p_comm,&
            ' repetition',left_rep,&
            ' pgrid coor ',left_prow, right_prow
       CALL mp_rma_get (left_set%data_windows(left_rep),&
            left_src_p_comm, 0,&
            all_sizes(1, left_rep, left_src_p_comm),&
            lbuffer_set_comm%mats(left_prow)%data)
       CALL mp_rma_get (left_set%index_windows(left_rep),&
            left_src_p_comm, 0,&
            all_sizes(2, left_rep, left_src_p_comm),&
            lbuffer_set_comm%mats(left_prow)%index)
    ENDDO left_row_aggregation
    DO left_rep = 1, left_reps
       CALL mp_window_fence (left_set%data_windows(left_rep),&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
       CALL mp_window_fence (left_set%index_windows(left_rep),&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
    ENDDO
    left_row_aggregation_repoint: DO left_prow_i = 0, left_nprows-1
       CALL dbcsr_repoint_index (lbuffer_set_comm%mats(left_prow),&
            error=error)
    ENDDO left_row_aggregation_repoint
    CALL dbcsr_switch_sets(lbuffer_set_calc, lbuffer_set_comm)
    CALL dbcsr_switch_sets(rbuffer_set_calc, rbuffer_set_comm)
    tstop = m_walltime ()
    WRITE (*,'(1X,A,F9.3)')routineN//' Setup time',tstop-tstart
    !
    ! Performs calculations and exchanges.
    k_index: DO right_prow_i = 1, right_nprows
       right_prow = MOD (right_prow_i + myprow, right_nprows)
       ! Right
       DBG 'right_prow from/to', right_prow_i, right_prow, myprow, right_nprows
       right_rep = 1
       wait_case_right: IF (right_prow_i .GT. 1) THEN
          DBG 'waiting right_data_sr'
          CALL mp_wait (right_data_sr(right_rep))
          DBG 'waiting right_data_rr'
          CALL mp_wait (right_data_rr(right_rep))
          DBG 'waiting right_index_sr'
          CALL mp_wait (right_index_sr(right_rep))
          DBG 'waiting right_index_rr'
          CALL mp_wait (right_index_rr(right_rep))
       ENDIF wait_case_right
       xfer_case_right: IF (right_prow_i .LT. right_nprows) THEN
          ! Sends what I will also be calculating now and receives what I will
          ! calculate in the next round.
          !right_src_p_comm = right_pgrid (right_prow, mypcol)
          right_src_p_comm = right_pgrid (&
               MOD (myprow+1+right_nprows, right_nprows), mypcol)
          right_dst_p_comm = right_pgrid (&
               MOD (myprow-1+right_nprows, right_nprows), mypcol)
          IF (dbg) WRITE(*,'(1X,A,"(",I3,",",I3,")")')&
               routineN//' Getting RIGHT matrix from pgrid coor ',&
               right_prow, mypcol
          IF (dbg) WRITE(*,'(1X,A,I4,A,"(",I3,",",I3,") sizes",I7,I7)')&
               routineN//' Getting RIGHT matrix from',right_src_p_comm,&
               ' pgrid coor ',right_prow, mypcol,&
               all_sizes(3:4, right_rep, right_src_p_comm)
          IF (dbg) WRITE(*,'(1X,A,I4,A,I4)')routineN//' EXCH Send to',&
               right_dst_p_comm,' recv_from',right_src_p_comm
          DBG 'sizes',SIZE(rbuffer_set_calc%mats(right_rep)%data),&
               SIZE(rbuffer_set_comm%mats(right_rep)%data)
          DBG 'sizes',SIZE(rbuffer_set_calc%mats(right_rep)%index),&
               SIZE(rbuffer_set_comm%mats(right_rep)%index)
          CALL mp_isendrecv (rbuffer_set_calc%mats(right_rep)%data,&
               right_dst_p_comm, rbuffer_set_comm%mats(right_rep)%data,&
               right_src_p_comm, mp_group, right_data_sr(right_rep),&
               right_data_rr(right_rep))
          CALL mp_isendrecv (rbuffer_set_calc%mats(right_rep)%index,&
               right_dst_p_comm, rbuffer_set_comm%mats(right_rep)%index,&
               right_src_p_comm, mp_group, right_index_sr(right_rep),&
               right_index_rr(right_rep))
          !CALL mp_sendrecv (rbuffer_set_calc%mats(right_rep)%data,&
          !     right_dst_p_comm, rbuffer_set_comm%mats(right_rep)%data,&
          !     right_src_p_comm, mp_group)
          !CALL mp_sendrecv (rbuffer_set_calc%mats(right_rep)%index,&
          !     right_dst_p_comm, rbuffer_set_comm%mats(right_rep)%index,&
          !     right_src_p_comm, mp_group)
       ENDIF xfer_case_right
       calc_case_right: IF (right_prow_i .GT. 0) THEN
          right_rep = 1
          CALL dbcsr_repoint_index (rbuffer_set_calc%mats(right_rep),&
               error=error)
       ENDIF calc_case_right
       right_src_p_calc = right_src_p_comm
       ! Left
       DO left_rep = 1, left_reps
          CALL mp_window_fence (left_set%data_windows(left_rep),&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
          CALL mp_window_fence (left_set%index_windows(left_rep),&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
       ENDDO
       left_rep = myprow+1 ! Because reps are offset 1 from myprow
       left_row_walk: DO left_prow_i = 0, left_nprows-1
          left_prow = MOD (left_prow_i + myprow, left_nprows)
          xfer_case_left: IF (right_prow_i .LT. right_nprows) THEN
             IF (dbg) WRITE(*,'(1X,A,I3,A,"(",I3,",",I3,")")')&
                  routineN//' Getting LEFT matrix repetition',left_rep,&
                  ' pgrid coor ',left_prow, right_prow
             left_src_p_comm = left_pgrid (left_prow, right_prow)
             IF (dbg) WRITE(*,'(1X,A,I4,A,I3,A,"(",I3,",",I3,")")')&
                  routineN//' Getting LEFT matrix from',left_src_p_comm,&
                  ' repetition',left_rep,&
                  ' pgrid coor ',left_prow, right_prow
             CALL mp_rma_get (left_set%data_windows(left_rep),&
                  left_src_p_comm, 0,&
                  all_sizes(1, left_rep, left_src_p_comm),&
                  lbuffer_set_comm%mats(left_prow)%data)
             CALL mp_rma_get (left_set%index_windows(left_rep),&
                  left_src_p_comm, 0,&
                  all_sizes(2, left_rep, left_src_p_comm),&
                  lbuffer_set_comm%mats(left_prow)%index)
          ENDIF xfer_case_left
          calc_case_left: IF (right_prow_i .GT. 0) THEN
             !left_src_p_calc = left_pgrid (left_prow, right_prow-1)
             CALL dbcsr_repoint_index (lbuffer_set_calc%mats(left_prow),&
                  error=error)
             CALL dbcsr_print (lbuffer_set_calc%mats(left_prow), error)
             CALL dbcsr_print (rbuffer_set_calc%mats(1), error)
             tstart = m_walltime ()
             CALL dbcsr_nn_mult(lbuffer_set_calc%mats(left_prow),&
                  rbuffer_set_calc%mats(1),&
                  product_matrix, product_wm, error, windex, flop=flop)
             flop_total = flop_total + flop
             tstop = m_walltime ()
             DBG 'Segment Local Multiplication time',tstop-tstart,'and',&
                  (REAL(flop, dp)/1000000.0d0)/(tstop-tstart),'MFLOPS'
             trun = trun + (tstop - tstart)
          ENDIF calc_case_left
       ENDDO left_row_walk
       left_src_p_calc = left_src_p_comm
       CALL dbcsr_switch_sets(lbuffer_set_calc, lbuffer_set_comm)
       CALL dbcsr_switch_sets(rbuffer_set_calc, rbuffer_set_comm)
    ENDDO k_index
    !
    !
    WRITE(*,'(1X,A,F9.3,A,F9.3,A,I9)')&
         routineN//' Total Local Multiplication time',&
         trun,' and',&
         (REAL(flop_total, dp)/1000000.0d0)/trun,' MFLOPS',flop_total
    DO rep = 1, left_reps
       CALL mp_window_destroy (left_set%data_windows(rep))
       CALL mp_window_destroy (left_set%index_windows(rep))
    ENDDO
    DO rep = 1, right_reps
       CALL mp_window_destroy (right_set%data_windows(rep))
       CALL mp_window_destroy (right_set%index_windows(rep))
    ENDDO
    !
    entries = btree_get_entries (windex)
    CALL ensure_array_size(product_wm%blk_p, ub=entries, error=error)
    ALLOCATE (keys(entries))
    CALL btree_delete (windex, keys, product_wm%blk_p)
    !lmask = 0 ; umask = 0
    !DO i = 0, 31
    !   lmask = IOR (lmask, IBSET (lmask, i))
    !   umask = IOR (umask, IBSET (umask, i+32))
    !ENDDO
    CALL ensure_array_size(product_wm%row_i, ub=entries, error=error)
    CALL ensure_array_size(product_wm%col_i, ub=entries, error=error)
    !DO i = 1, entries
    !   product_wm%row_i(i) = ISHFT (keys(i), -32)
    !   product_wm%col_i(i) = IAND (keys(i), lmask)
    !ENDDO
    CALL dbcsr_unpack_i8_2i4 (keys, product_wm%row_i, product_wm%col_i)
    DEALLOCATE (keys)
    CALL dbcsr_finalize(product_wm, product_matrix, .TRUE., .FALSE., error)
    CALL timestop (timing_handle)
  END SUBROUTINE dbcsr_multiply_sets_sqrt


! *****************************************************************************
!> \brief Multiplies two DBCSR matrices, using the current local data.
!> \param[in] left, right     left and right DBCSR matrices
!> \param[in,out] product     resulting DBCSR product matrix
!> \param[in,out] wm          DBCSR product working matrix
!> \param[in,out] error       cp2k error
!> \param[out] flop           (optional) number of effective double-precision
!>                            floating point operations performed
! *****************************************************************************
  SUBROUTINE dbcsr_nn_mult(left, right, product, wm, error, index, flop)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: wm
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    TYPE(btree_I8_I4), INTENT(INOUT)         :: index
    INTEGER, INTENT(OUT), OPTIONAL           :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_nn_mult', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: alpha = 1.0d0, &
                                                beta_add = 1.0d0, &
                                                beta_new = 0.0d0

    CHARACTER                                :: a_type, b_type
    INTEGER :: a_blk, a_col, a_col_size, a_nze, a_row, a_row_size, b_blk, &
      b_col, b_col_size, b_nze, b_row, b_row_size, c_blk_p, c_col, c_nze, &
      c_row, lda, ldb, lflop, timing_handle
    INTEGER(KIND=int_4)                      :: offset
    INTEGER(KIND=int_8)                      :: coordinate
    LOGICAL                                  :: block_exists
    REAL(KIND=dp)                            :: beta, tstart, tstop
    REAL(KIND=dp), DIMENSION(:), POINTER     :: a_data, b_data, c_data

!.!INTEGER, DIMENSION(:), ALLOCATABLE   :: c_blk_ps

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    !IF (dbg) tstart = m_walltime ()
    tstart = m_walltime ()
    DBG "Doing NN matrix multiplication",left%name,' x ',right%name
    DBG "Pre-mult blk sizes, data", wm%lastblk, wm%datasize
    ! c_blk_ps contains the block pointers in the product's current row
    !.!ALLOCATE (c_blk_ps(product%nblkcols_total), stat=stat)
    !IF (stat /= 0) CALL stop_memory (routineN, "c_blk_ps",&
    !     product%nblkcols_total)
    lflop = 0
    c_blk_p = wm%datasize + 1
    DBG 'A rows:',left%nblkrows_total
    IF (left%nblks .GT. 0 .AND. right%nblks .GT. 0) THEN
    !CALL ensure_array_size(wm%data,&
    !     ub=wm%datasize+left%nfullrows_total*right%nfullcols_total, error=error)
    DO a_row = 1, left%nblkrows_total
       DBG 'A row', a_row,' with #cols=',left%row_p(a_row+1)-left%row_p(a_row)
       a_row_size = left%row_blk_size(a_row)
       CALL ensure_array_size(wm%data,&
            ub=wm%datasize+a_row_size*right%nfullcols_total &
            * MAX(1, (left%nblkrows_total-a_row)/4), error=error)
       DO a_blk = left%row_p(a_row)+1, left%row_p(a_row+1)
          a_col = left%col_i(a_blk)
          DBG 'A col', a_col
          b_row = a_col ; c_row = a_row
          a_col_size = left%col_blk_size(a_col)
          a_nze = a_row_size * a_col_size
          b_row_size = a_col_size
          a_data => left%data(ABS(left%blk_p(a_blk))&
                             :ABS(left%blk_p(a_blk))+a_nze-1)
          IF (left%blk_p(a_blk) .LT. 0) THEN
             a_type = 'T'; lda = a_col_size
          ELSE
             a_type = 'N'; lda = a_row_size
          ENDIF
          DO b_blk = right%row_p(b_row)+1, right%row_p(b_row+1)
             b_col = right%col_i(b_blk) ; c_col = b_col
             DBG 'B col', b_col
             b_col_size = right%col_blk_size(b_col)
             b_nze = b_row_size * b_col_size
             b_data => right%data(&
                  ABS(right%blk_p(b_blk)):ABS(right%blk_p(b_blk))+b_nze-1)
             IF (right%blk_p(b_blk) .LT. 0) THEN
                b_type = 'T'; ldb = b_col_size
             ELSE
                b_type = 'N'; ldb = b_row_size
             ENDIF
             c_nze = a_row_size * b_col_size
             ! If this c column does not yet exist, "allocate" it.
             coordinate = make_coordinate_tuple (c_row, c_col)
             CALL btree_find (index, coordinate, offset, block_exists)
             !.!IF (c_blk_ps(c_col) .EQ. 0) THEN
             IF (.NOT. block_exists) THEN
                !.!c_blk_ps(c_col) = c_blk_p
                !.!CALL add_work_coordinate(wm, c_row, c_col, c_blk_p, error=error)
                offset = c_blk_p
                CALL btree_add (index, coordinate, offset)
                wm%lastblk = wm%lastblk + 1
                c_blk_p = c_blk_p + c_nze
                wm%datasize = wm%datasize + c_nze
                beta = beta_new
             ELSE
                beta = beta_add
             ENDIF
             !.!c_data => wm%data(c_blk_ps(c_col):c_blk_ps(c_col)+c_nze-1)
             DBG 'offset:',offset,offset+c_nze-1,'of',wm%datasize,UBOUND(wm%data),'ex',block_exists
             c_data => wm%data(offset:offset+c_nze-1)
             IF (dbg) WRITE(*,&
                  "(1X,A,A1,1X,A1,A,I6,A,F4.1,' (',I3,',',I3,')=(',I3,',',I3,')x(',I3,',',I3,')')")&
                  !routineN//' DGEMM ',a_type,b_type,' at',c_blk_ps(c_col),' beta=', beta,&
                  routineN//' DGEMM ',a_type,b_type,' at',offset,' beta=', beta,&
                  c_row,c_col,a_row,a_col,b_row,b_col
             CALL DGEMM(a_type, b_type, a_row_size, b_col_size, a_col_size,&
                  alpha,&
                  a_data, lda,&
                  b_data, ldb, beta,&
                  c_data, a_row_size)
             lflop = lflop + a_row_size * b_col_size * a_row_size * 2.0d0
          ENDDO
       ENDDO ! a_col
    ENDDO ! a_row
    ELSE
       DBG 'One of the matrices is empty!'
    ENDIF
    !.!DEALLOCATE (c_blk_ps)
    DBG "Post-mult blk sizes, data", wm%lastblk, wm%datasize
    dbg = .TRUE.
    IF (dbg) THEN
       tstop = m_walltime ()
       DBG tstop-tstart,'s',&
            (REAL(lflop, dp)/1000000.0d0)/(tstop-tstart),'MFLOPS'
    ENDIF
    IF (PRESENT (flop)) flop = lflop
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_nn_mult

  FUNCTION make_coordinate_tuple (most, least) RESULT (tuple)
    INTEGER, INTENT(IN)                      :: most, least
    INTEGER(KIND=int_8)                      :: tuple

!tuple = IOR (ISHFT (most, 32), least)

    tuple = most
    tuple = IOR (ISHFT(tuple, 32), INT(least, int_8))
  END FUNCTION make_coordinate_tuple
    


!! *****************************************************************************
!!> \brief Temporary transpose of a matrix represented with 1-d data
!!> \param[in] vector          matrix data
!!> \param[in] rows            number of matrix rows
!!> \param[in] columns         number of matrix columns
!! *****************************************************************************
!
!  FUNCTION temp_transpose (vector, rows, columns) RESULT vector_tr
!    REAL(KIND=dp)                           :: vector_tr
!    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vector
!    INTEGER, INTENT(IN)                     :: rows, columns
!!   ---------------------------------------------------------------------------
!    vector_tr =
!    RESHAPE (&
!            TRANSPOSE (&
!                    RESHAPE (vector, (/ rows, columns /))),&
!            (/ rows*columns /))
!  END FUNCTION temp_transpose


! *****************************************************************************
!> \brief Adds two DBCSR work matrices
!> \param[in,out] wm          resulting DBCSR sum work matrix
!> \param[in] addee           DBCSR work matrix to add to the wm matrix
!> \param[out] flop           (optional) number of effective double-precision
!>                            floating point operations performed
! *****************************************************************************
  SUBROUTINE dbcsr_wm_add(wm, addee, error, flop)
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: wm
    TYPE(dbcsr_type), INTENT(IN)             :: addee
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER, INTENT(OUT), OPTIONAL           :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_wm_add', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: alpha = 1.0d0, &
                                                beta_add = 1.0d0, &
                                                beta_new = 0.0d0

    CHARACTER                                :: addee_type, sum_type
    INTEGER :: addee_blk, addee_col, addee_col_size, addee_nze, addee_row, &
      addee_row_size, ldaddee, ldsum, lflop, sum_blk, sum_blk_p, &
      sum_col_size, sum_row_size, timing_handle
    REAL(KIND=dp)                            :: tstart, tstop
    REAL(KIND=dp), DIMENSION(:), POINTER     :: addee_data, sum_data

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    IF (dbg) tstart = m_walltime ()
    DBG "Doing WM matrix sum",' += ',addee%name
    DBG "Pre-mult blk sizes, data", wm%lastblk, wm%datasize
    DBG "Addee blk sizes, data", addee%nblks, addee%nze
    !
    ! The algorithm is a follows: go through all of addee's
    ! blocks. For every block, see if a corresponding block already
    ! exists in wm. If such a block exists, then add the two
    ! blocks. In this case, the index does not change. If a
    ! corresponding block does not exist, then simply add tho new
    ! block. The data goes to the end of the data array and the index
    ! information is appended to wm's index structure.
    !
    ! c_blk_ps contains the block pointers in the product's current row
    lflop = 0
    IF (wm%lastblk .GT. 0 .OR. addee%nblks .GT. 0) THEN
    DO addee_row = 1, addee%nblkrows_total
       DBG 'Addee row', addee_row,' with #cols=',&
            addee%row_p(addee_row+1)-addee%row_p(addee_row)
       ! Find the blocks present in addee_row.
       addee_row_size = addee%row_blk_size(addee_row)
       CALL ensure_array_size(wm%data,&
            ub=wm%datasize+addee_row_size*addee%nfullcols_total &
            * MAX(1, (addee%nblkrows_total-addee_row)/4), error=error)
       DO addee_blk = addee%row_p(addee_row)+1, addee%row_p(addee_row+1)
          addee_col = addee%col_i(addee_blk)
          DBG 'Addee col', addee_col
          addee_col_size = addee%col_blk_size(addee_col)
          addee_nze = addee_row_size * addee_col_size
          addee_data => addee%data(ABS(addee%blk_p(addee_blk))&
                             :ABS(addee%blk_p(addee_blk))+addee_nze-1)
          IF (addee%blk_p(addee_blk) .LT. 0) THEN
             addee_type = 'T'; ldaddee = addee_col_size
          ELSE
             addee_type = 'N'; ldaddee = addee_row_size
          ENDIF
          ! Attempt to find co-lying block in the sum wm matrix.
          sum_blk_p = 0 ! Signifies colier not present.
          find_sum_colier: DO sum_blk = 1, wm%lastblk
             IF (wm%row_i(sum_blk) .EQ. addee_row .AND.&
                  wm%col_i(sum_blk) .EQ. addee_col) THEN
                sum_blk_p = wm%blk_p(sum_blk)
                EXIT find_sum_colier
             ENDIF
          ENDDO find_sum_colier
          ! If the colier is present, then sum the data, otherwise
          ! augment the wm sum matrix with the data and coordinates.
          IF (sum_blk_p .GT. 0) THEN
             sum_data => wm%data(&
                  ABS(sum_blk_p):ABS(sum_blk_p)+addee_nze-1)
             IF (wm%blk_p(sum_blk) .LT. 0) THEN
                sum_type = 'T'; ldsum = sum_col_size
             ELSE
                sum_type = 'N'; ldsum = sum_row_size
             ENDIF
             IF (sum_type .EQ. addee_type) THEN
                CALL DAXPY (addee_nze, 1.0d0, addee_data, 1, sum_data, 1)
             ELSE
                CALL DAXPY (addee_nze, 1.0d0,&
temp_transpose(addee_data,addee_row_size,addee_col_size),&
                     1, sum_data, 1)
             ENDIF
             lflop = lflop + addee_nze
          ELSE
             sum_blk_p = wm%datasize + 1
             sum_data => wm%data(&
                  ABS(sum_blk_p):ABS(sum_blk_p)+addee_nze-1)
             sum_data = addee_data
             CALL add_work_coordinate(wm, addee_row, addee_col, sum_blk_p,&
                  error=error)
             wm%datasize = wm%datasize + addee_nze
          END IF
       ENDDO ! addee_blk
    ENDDO ! addee_row
    ELSE
       DBG 'Both matrices were empty!'
    ENDIF
    DBG "Post-add blk sizes, data", wm%lastblk, wm%datasize
    IF (dbg) THEN
       tstop = m_walltime ()
       DBG tstop-tstart,'s',&
            (REAL(lflop, dp)/1000000.0d0)/(tstop-tstart),'MFLOPS'
    ENDIF
    IF (PRESENT (flop)) flop = lflop
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_wm_add


END MODULE dbcsr_operations

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Higher-level operations on DBCSR matrices.
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.9
!>
!> <b>Modification history:</b>
!  - Created 2009-05-12
! *****************************************************************************
MODULE dbcsr_operations
  USE array_types,                     ONLY: array_data
  USE btree_I8_k_I4_v,                 ONLY: btree_add,&
                                             btree_delete,&
                                             btree_find,&
                                             btree_get_entries,&
                                             btree_new,&
                                             btree_verify
  USE dbcsr_methods,                   ONLY: &
       dbcsr_distribution_mp, dbcsr_mp_group, dbcsr_mp_mynode, &
       dbcsr_mp_mypcol, dbcsr_mp_myprow, dbcsr_mp_npcols, dbcsr_mp_nprows, &
       dbcsr_mp_numnodes, dbcsr_mp_pgrid
  USE dbcsr_types,                     ONLY: dbcsr_array_type,&
                                             dbcsr_mp_obj,&
                                             dbcsr_type,&
                                             dbcsr_work_type
  USE dbcsr_util,                      ONLY: &
       add_work_coordinate, bcsr_assert, dbcsr_create, dbcsr_destroy, &
       dbcsr_destroy_array, dbcsr_finalize, dbcsr_init, dbcsr_repoint_index, &
       dbcsr_set_debug, dbcsr_unpack_i8_2i4, ensure_array_size, swap
  USE kinds,                           ONLY: dp,&
                                             int_4,&
                                             int_8
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_rma_get,&
                                             mp_sum,&
                                             mp_window_create,&
                                             mp_window_destroy,&
                                             mp_window_fence,&
                                             mp_window_type

  !$ USE OMP_LIB
  USE termination,                     ONLY: stop_memory,stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_operations'

  PUBLIC :: dbcsr_multiply_from_left, dbcsr_nn_mult, dbcsr_trace, dbcsr_scale

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
!> \brief Multiplies two DBCSR matrices.
!> \param[in] left_set             set of left images
!> \param[in] right_matrix         right matrix
!> \param[out] product             DBCSR product workmatrix
!> \param[in,out] error            cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_from_left (left_set, right_matrix, product_matrix,&
       error)
    TYPE(dbcsr_array_type), INTENT(INOUT)    :: left_set
    TYPE(dbcsr_type), INTENT(IN)             :: right_matrix
    TYPE(dbcsr_type), INTENT(INOUT)          :: product_matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_from_left', &
      routineP = moduleN//':'//routineN

    INTEGER :: debug_level, entries, flop, flop_total, i, img, left_image, &
      left_images, left_multiplicity, left_nprows, left_prow, &
      left_src_p_calc, left_src_p_comm, mp_group, mult, mynode, mypcol, &
      myprow, nblks, npcols, nprows, numnodes, nze, output_unit, &
      right_nprows, right_prow, right_prow_i, right_src_p_calc, &
      right_src_p_comm, timing_handle
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: keys
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: my_sizes
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: all_sizes
    INTEGER, DIMENSION(:, :), POINTER        :: left_pgrid, product_pgrid, &
                                                right_pgrid
    REAL(KIND=dp)                            :: trun, trun_t, tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: my_waits
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: all_waits
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_array_type), POINTER          :: lbuffer_calc, lbuffer_comm
    TYPE(dbcsr_array_type), TARGET           :: lbuffer_1, lbuffer_2
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_type), POINTER                :: rbuffer_calc, rbuffer_comm
    TYPE(dbcsr_type), TARGET                 :: rbuffer_1, rbuffer_2
    TYPE(mp_window_type)                     :: right_data_window, &
                                                right_index_window

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    debug_level=0
    IF (dbg) debug_level = 2
    CALL timeset (routineN, timing_handle)
    CALL btree_new (product_matrix%w%tree)
    tstart = 0.0_dp ; tstop = 0.0_dp ; trun = 0.0_dp
    flop_total = 0
    trun_t = m_walltime ()
    left_images = SIZE (left_set%mats)
    CALL bcsr_assert (debug_level,&
         left_images .EQ. left_set%image_dist%row_decimation, routineP,&
         "Non-compatible images")
    left_multiplicity = left_set%image_dist%row_multiplicity
    mp_obj = dbcsr_distribution_mp (product_matrix%dist)
    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    npcols = dbcsr_mp_npcols (mp_obj)
    nprows = dbcsr_mp_nprows (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    product_pgrid => dbcsr_mp_pgrid (mp_obj)
    left_pgrid => dbcsr_mp_pgrid (dbcsr_distribution_mp (&
         left_set%image_dist%main))
    left_nprows = dbcsr_mp_nprows (dbcsr_distribution_mp (&
         left_set%image_dist%main))
    right_pgrid => dbcsr_mp_pgrid (dbcsr_distribution_mp (&
         right_matrix%dist))
    DBG 'product pgrid shape', SHAPE (dbcsr_mp_pgrid (mp_obj))
    DBG 'left pgrid shape', SHAPE (left_pgrid),';',left_nprows,left_images
    DBG 'right pgrid shape', SHAPE (right_pgrid)
    DBG 'left images, multiplicity', left_images, left_multiplicity
    CALL bcsr_assert (debug_level,&
         left_multiplicity * npcols .EQ. left_nprows, routineP,&
         "Mulitplicity mismatch and processer grid.")
    ALLOCATE (left_set%data_windows(left_images))
    ALLOCATE (left_set%index_windows(left_images))
    ALLOCATE (my_sizes(4, left_images))
    my_sizes(:,:) = 0
    DO img = 1, left_images
       my_sizes(1, img) = SIZE (left_set%mats(img)%data)
       my_sizes(2, img) = SIZE (left_set%mats(img)%index)
    ENDDO
    my_sizes(3, 1) = SIZE (right_matrix%data)
    my_sizes(4, 1) = SIZE (right_matrix%index)
    ALLOCATE (all_sizes(4, LBOUND(my_sizes,2):UBOUND(my_sizes,2), 0:numnodes-1))
    DO i = 1, 4
       DBG 'my sizes', my_sizes(i,:)
    ENDDO
    CALL mp_allgather(my_sizes, all_sizes, mp_group)
    nblks = MAXVAL (all_sizes (2, :, :))
    nze = MAXVAL (all_sizes(1, :, :))
    DBG 'data sizes', all_sizes
    ! Setup the buffer matrices
    DO img = 1, left_images
       CALL mp_window_create (left_set%data_windows(img),&
            left_set%mats(img)%data, gid=mp_group)
       CALL mp_window_create (left_set%index_windows(img),&
            left_set%mats(img)%index, gid=mp_group)
    ENDDO
    ALLOCATE (lbuffer_1%mats(left_multiplicity))
    ALLOCATE (lbuffer_2%mats(left_multiplicity))
    DO mult = 1, left_multiplicity
       CALL dbcsr_init (lbuffer_1%mats(mult), error)
       CALL dbcsr_create (lbuffer_1%mats(mult),&
            'lbufferl1'//left_set%mats(1)%name, left_set%mats(1)%dist, 'N',&
            left_set%mats(1)%row_blk_size, left_set%mats(1)%col_blk_size,&
            nblks, nze, special=.FALSE., error=error)
       CALL ensure_array_size (lbuffer_1%mats(mult)%data,&
            ub=MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
       CALL ensure_array_size (lbuffer_1%mats(mult)%index,&
            ub=MAXVAL (all_sizes(2, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
       !
       CALL dbcsr_init (lbuffer_2%mats(mult), error)
       CALL dbcsr_create (lbuffer_2%mats(mult),&
            'lbuffer2'//left_set%mats(1)%name, left_set%mats(1)%dist, 'N',&
            left_set%mats(1)%row_blk_size, left_set%mats(1)%col_blk_size,&
            nblks, nze, special=.FALSE., error=error)
       CALL ensure_array_size (lbuffer_2%mats(mult)%data,&
            ub=MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
       CALL ensure_array_size (lbuffer_2%mats(mult)%index,&
            ub=MAXVAL (all_sizes(2, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
    ENDDO
    lbuffer_calc => lbuffer_1
    lbuffer_comm => lbuffer_2
    ! Setup buffers for right matrix
    nblks = MAXVAL (all_sizes (4, 1, :))
    nze = MAXVAL (all_sizes(3, 1, :))
    CALL dbcsr_init (rbuffer_1, error)
    CALL dbcsr_create (rbuffer_1,&
         'rbuffer1'//right_matrix%name, right_matrix%dist,&
         'N',&
         right_matrix%row_blk_size,&
         right_matrix%col_blk_size,&
         0, 0, special=.FALSE., error=error)
    CALL ensure_array_size (rbuffer_1%data,&
         ub=MAXVAL (all_sizes(3, 1, :)), nocopy=.TRUE., special=.FALSE.,&
         error=error)
    CALL ensure_array_size (rbuffer_1%index,&
         ub=MAXVAL (all_sizes(4, 1, :)), nocopy=.TRUE., special=.FALSE.,&
         error=error)
    !
    CALL dbcsr_init (rbuffer_2, error)
    CALL dbcsr_create (rbuffer_2,&
         'rbuffer2'//right_matrix%name, right_matrix%dist,&
         'N',&
         right_matrix%row_blk_size,&
         right_matrix%col_blk_size,&
         0, 0, special=.FALSE., error=error)
    CALL ensure_array_size (rbuffer_2%data,&
         ub=MAXVAL (all_sizes(3, 1, :)), nocopy=.TRUE., special=.FALSE.,&
         error=error)
    CALL ensure_array_size (rbuffer_2%index,&
         ub=MAXVAL (all_sizes(4, 1, :)), nocopy=.TRUE., special=.FALSE.,&
         error=error)
    ! Create windows for the right matrices
    DBG 'RIGHT data win size', SIZE (right_matrix%data)
    DBG 'RIGHT indx win size', SIZE (right_matrix%index)
    CALL mp_window_create (right_data_window, right_matrix%data, gid=mp_group)
    CALL mp_window_create (right_index_window, right_matrix%index, gid=mp_group)
    rbuffer_calc => rbuffer_1
    rbuffer_comm => rbuffer_2
    !
    right_nprows = dbcsr_mp_nprows (dbcsr_distribution_mp (&
         right_matrix%dist))
    right_src_p_calc = -7
    left_src_p_calc = -7
    ! Do the actual transfers and call local multiplication
    ALLOCATE (my_waits(4, 0:right_nprows))
    my_waits(:,:) = 0
    k_index: DO right_prow_i = 0, right_nprows
       ! Right matrix
       right_prow = MOD (right_prow_i+myprow+mypcol, right_nprows)
       DBG 'RIGHT_PROW from/to', right_prow_i, right_prow, myprow, right_nprows
       my_waits(1, right_prow_i) = m_walltime()
       !CALL dbcsr_sleep (10.0_dp)
       CALL mp_window_fence (right_data_window,&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
       my_waits(1, right_prow_i) = m_walltime() - my_waits(1, right_prow_i)
       my_waits(2, right_prow_i) = m_walltime()
       !CALL dbcsr_sleep (10.0_dp)
       CALL mp_window_fence (right_index_window,&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
       my_waits(2, right_prow_i) = m_walltime() - my_waits(2, right_prow_i)
       ! No fetches in the last loop
       xfer_case_right: IF (right_prow_i .LT. right_nprows) THEN
          right_src_p_comm = right_pgrid (right_prow, mypcol)
          IF (dbg) WRITE(*,'(1X,A,"(",I3,",",I3,")")')&
               routineN//' Getting RIGHT matrix from pgrid coor ',&
               right_prow, mypcol
          IF (dbg) WRITE(*,'(1X,A,I4,A,"(",I3,",",I3,") sizes",I7,I7)')&
               routineN//' Getting RIGHT matrix from',right_src_p_comm,&
               ' pgrid coor ',right_prow, mypcol,&
               all_sizes(3:4, 1, right_src_p_comm)
          CALL mp_rma_get (right_data_window,&
               right_src_p_comm, 0, all_sizes(3, 1, right_src_p_comm),&
               rbuffer_comm%data)
          CALL mp_rma_get (right_index_window,&
               right_src_p_comm, 0, all_sizes(4, 1, right_src_p_comm),&
               rbuffer_comm%index)
       ENDIF xfer_case_right
       ! No calculations in the first loop
       calc_case_right: IF (right_prow_i .GT. 0) THEN
          IF (all_sizes(4, 1, right_src_p_calc) .EQ. 0) THEN
             DBG 'ZERO SIZED INDEX'
          ENDIF
          CALL dbcsr_repoint_index (rbuffer_calc, error=error)
       ENDIF calc_case_right
       ! Now we transfer/multiply with the left matrix
       right_src_p_calc = right_src_p_comm
       DO left_image = 1, left_images
          my_waits(3, right_prow_i) = m_walltime()
          !CALL dbcsr_sleep (10.0_dp)
          CALL mp_window_fence (left_set%data_windows(left_image),&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
          my_waits(3, right_prow_i) = m_walltime() - my_waits(3, right_prow_i)
          my_waits(4, right_prow_i) = m_walltime()
          !CALL dbcsr_sleep (10.0_dp)
          CALL mp_window_fence (left_set%index_windows(left_image),&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
       my_waits(4, right_prow_i) = m_walltime() - my_waits(4, right_prow_i)
       ENDDO
       left_prow = myprow / left_images
       left_image = 1 + MOD (myprow, left_images)
       xfer_case_left: IF (right_prow_i .LT. right_nprows) THEN
          DO mult = 1, left_multiplicity
             IF (dbg) WRITE(*,'(1X,A,I3,A,"(",I3,",",I3,")")')&
                  routineN//' Getting LEFT matrix image',left_image,&
                  ' pgrid coor ',left_prow, right_prow
             IF (dbg) THEN
                WRITE(*,*)routineP//'left pgrid shape', SHAPE (left_pgrid),';'
                WRITE(*,*)routineP, left_prow,&
                     (mult-1)*(left_nprows/left_multiplicity),&
                     '=>',left_prow + (mult-1)*(left_nprows/left_multiplicity)
             ENDIF
             left_src_p_comm = left_pgrid (&
                  left_prow + (mult-1)*(left_nprows/left_multiplicity),&
                  right_prow)
             IF (dbg) WRITE(*,'(1X,A,I4,A,I3,A,"(",I3,",",I3,")")')&
                  routineN//' Getting LEFT matrix from',left_src_p_comm,&
                  ' image',left_image,&
                  ' pgrid coor ',left_prow, right_prow
             CALL mp_rma_get (left_set%data_windows(left_image),&
                  left_src_p_comm, 0,&
                  all_sizes(1, left_image, left_src_p_comm),&
                  lbuffer_comm%mats(mult)%data)
             CALL mp_rma_get (left_set%index_windows(left_image),&
                  left_src_p_comm, 0,&
                  all_sizes(2, left_image, left_src_p_comm),&
                  lbuffer_comm%mats(mult)%index)
          ENDDO
       ENDIF xfer_case_left
       calc_case_left: IF (right_prow_i .GT. 0) THEN
             !left_src_p_calc = left_pgrid (left_prow, right_prow-1)
          DO mult = 1, left_multiplicity
             CALL dbcsr_repoint_index (lbuffer_calc%mats(mult), error=error)
             tstart = m_walltime ()
             CALL dbcsr_nn_mult(lbuffer_calc%mats(mult), rbuffer_calc,&
                  product_matrix, error, flop=flop)
             flop_total = flop_total + flop
             tstop = m_walltime ()
             IF (tstop-tstart .NE. 0) THEN
                DBG 'Segment Local Multiplication time',tstop-tstart,'and',&
                     (REAL(flop, dp)/1000000.0_dp)/(tstop-tstart),'MFLOPS'
             ENDIF
             trun = trun + (tstop - tstart)
          ENDDO
       ENDIF calc_case_left
       left_src_p_calc = left_src_p_comm
       CALL dbcsr_switch_sets (lbuffer_calc, lbuffer_comm)
       CALL dbcsr_switch_m_ptrs (rbuffer_calc, rbuffer_comm)
    ENDDO k_index
    trun_t = m_walltime () - trun_t
    IF (trun_t .EQ. 0) trun_t = 0.001_dp
    IF (dbg) THEN
       WRITE(*,'(1X,A,F9.3,A,F9.3,A,I9)')&
            routineN//' Total Local Multiplication time',&
            trun_t,' and',&
            (REAL(flop_total, dp)/1000000.0_dp)/trun_t,' MFLOPS',flop_total
    ENDIF
    DO img = 1, left_images
       CALL mp_window_destroy (left_set%data_windows(img))
       CALL mp_window_destroy (left_set%index_windows(img))
    ENDDO
    DEALLOCATE (left_set%data_windows)
    DEALLOCATE (left_set%index_windows)
    CALL mp_window_destroy (right_data_window)
    CALL mp_window_destroy (right_index_window)
    !
    CALL dbcsr_destroy (rbuffer_1, error)
    CALL dbcsr_destroy (rbuffer_2, error)
    CALL dbcsr_destroy_array (lbuffer_1, error)
    CALL dbcsr_destroy_array (lbuffer_2, error)
    DEALLOCATE (my_sizes)
    !
    ALLOCATE (all_waits(4, 0:right_nprows, 0:numnodes-1))
    CALL mp_allgather (my_waits, all_waits, mp_group)
    logger => cp_error_get_logger (error)
    IF (logger%para_env%mepos == logger%para_env%source) THEN
       output_unit = cp_logger_get_default_unit_nr (logger)
       WRITE(output_unit,'(A)')"Waiting times per place per step"
       DO right_prow_i = 0, right_nprows
          WRITE(output_unit,'(4(F9.4),5X,4(F9.4))')&
               (MAXVAL(all_waits(i,right_prow_i,:)), i=1, 4),&
               (SUM(all_waits(i,right_prow_i,:))/numnodes, i=1, 4)
          WRITE(*,'(A,I4)')"Step ", right_prow_i
          DO i = 0, numnodes-1
             WRITE(*,'(I4,1X,4(F9.4))')i, all_waits(:, right_prow_i, i)
          ENDDO
       ENDDO
    ENDIF
    !
    entries = btree_get_entries (product_matrix%w%tree)
    CALL ensure_array_size (product_matrix%w%blk_p, ub=entries, error=error)
    ALLOCATE (keys (entries))
    CALL btree_delete (product_matrix%w%tree, keys, product_matrix%w%blk_p)
    CALL ensure_array_size(product_matrix%w%row_i, ub=entries, error=error)
    CALL ensure_array_size(product_matrix%w%col_i, ub=entries, error=error)
    CALL dbcsr_unpack_i8_2i4 (keys, product_matrix%w%row_i,&
         product_matrix%w%col_i)
    DEALLOCATE (keys)
    CALL dbcsr_finalize(product_matrix, .TRUE., .FALSE., error)
    !
    CALL timestop (timing_handle)
  END SUBROUTINE dbcsr_multiply_from_left


! *****************************************************************************
!> \brief Switches pointers between two matrices
!> \param[in,out] set1p, set2p
! *****************************************************************************
  SUBROUTINE dbcsr_switch_m_ptrs (m1p, m2p)
    TYPE(dbcsr_type), POINTER                :: m1p, m2p

    TYPE(dbcsr_type), POINTER                :: tmp_p

!   ---------------------------------------------------------------------------

    tmp_p => m1p
    m1p => m2p
    m2p => tmp_p
  END SUBROUTINE dbcsr_switch_m_ptrs


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
  ELEMENTAL SUBROUTINE make_tag (tag, to, from, seq)
    INTEGER, INTENT(OUT)                     :: tag
    INTEGER, INTENT(IN)                      :: to, from, seq

    INTEGER, PARAMETER                       :: s = 8

!   ---------------------------------------------------------------------------

    tag = ISHFT (to, s) + from + ISHFT (seq, 2*s)
  END SUBROUTINE make_tag


! *****************************************************************************
!> \brief Multiplies two DBCSR matrices, using the current local data.
!> \param[in] left, right     left and right DBCSR matrices
!> \param[in,out] product     resulting DBCSR product matrix
!> \param[in,out] wm          DBCSR product working matrix
!> \param[in,out] error       cp2k error
!> \param[out] flop           (optional) number of effective double-precision
!>                            floating point operations performed
! *****************************************************************************
  SUBROUTINE dbcsr_nn_mult(left, right, product, error, flop)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER, INTENT(OUT), OPTIONAL           :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_nn_mult', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: alpha = 1.0_dp, &
                                                beta_add = 1.0_dp, &
                                                beta_new = 0.0_dp

    CHARACTER                                :: a_type, b_type
    INTEGER :: a_blk, a_col, a_col_log, a_col_size, a_col_size_logical, &
      a_nze, a_row, a_row_log, a_row_size, a_row_size_logical, b_blk, b_col, &
      b_col_size, b_nze, b_row, b_row_size, c_blk, c_blk_p, c_col, c_nze, &
      c_row, lda, ldb, ldc, lflop, num_symmetries, stat, symmetry, &
      timing_handle
    INTEGER(KIND=int_4)                      :: offset
    INTEGER(KIND=int_8)                      :: coordinate
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: c_blk_ps
    INTEGER, DIMENSION(:), POINTER           :: left_col_blk_size, &
                                                left_row_blk_size, &
                                                right_col_blk_size, &
                                                right_row_blk_size
    LOGICAL                                  :: block_exists, tr_a, tr_b
    REAL(KIND=dp)                            :: beta, tstart, tstop
    REAL(KIND=dp), DIMENSION(:), POINTER     :: a_data, b_data, c_data
    TYPE(dbcsr_work_type)                    :: wm

    !$ INTEGER                           :: ithread, nthreads
    !$ INTEGER(KIND=omp_lock_kind)       :: data_lock, index_lock
    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    wm = product%w
    DBG "Doing NN matrix multiplication",left%name,' x ',right%name
    DBG "Pre-mult blk sizes, data", wm%lastblk, wm%datasize
    ! c_blk_ps contains the block pointers in the product's current row
    !.!ALLOCATE (c_blk_ps(product%nblkcols_total), stat=stat)
    !IF (stat /= 0) CALL stop_memory (routineN, "c_blk_ps",&
    !     product%nblkcols_total)
    lflop = 0
    num_symmetries = 1
    left_row_blk_size => array_data (left%row_blk_size)
    left_col_blk_size => array_data (left%col_blk_size)
    right_row_blk_size => array_data (right%row_blk_size)
    right_col_blk_size => array_data (right%col_blk_size)
    c_blk_p = wm%datasize + 1
    DBG 'A rows:',left%nblkrows_total
    IF (left%nblks .GT. 0 .AND. right%nblks .GT. 0) THEN
    CALL ensure_array_size(wm%data,&
         ub=MAX(wm%datasize,left%nfullrows_total*right%nfullcols_total), error=error)
       !CALL ensure_array_size(wm%data,&
       !     ub=wm%datasize+a_row_size*right%nfullcols_total &
       !     * MAX(1, (left%nblkrows_total-a_row)/4), error=error)
    !$OMP PARALLEL DEFAULT (SHARED) PRIVATE(a_row, a_col, a_blk,&
    !$OMP                                   b_row, b_col, b_blk,&
    !$OMP                                   c_row, c_col, c_blk,&
    !$OMP                                   a_row_log, a_col_log,&
    !$OMP                                   a_row_size_logical,&
    !$OMP                                   a_col_size_logical,&
    !$OMP                                   a_row_size, a_col_size, a_nze,&
    !$OMP                                   b_row_size, b_col_size, b_nze, c_nze,&
    !$OMP                                   lda, ldb, ldc,&
    !$OMP                                   tr_a, tr_b,&
    !$OMP                                   c_blk_ps,&
    !$OMP                                   a_data, b_data, c_data,&
    !$OMP                                   beta,&
    !$OMP                                   a_type, b_type,&
    !$OMP                                   symmetry,&
    !$OMP                                   coordinate,&
    !$OMP                                   offset,&
    !$OMP                                   block_exists,&
    !$OMP                                   tstart, tstop,&
    !$OMP                                   ithread)&
    !$OMP          REDUCTION (+: lflop)
    !
    !$ CALL OMP_INIT_LOCK (index_lock)
    !$ CALL OMP_INIT_LOCK (data_lock)
    !$ ithread = OMP_GET_THREAD_NUM ()
    !!$ WRITE (*,*)"Hello! I'm thread number", ithread
    !$OMP SINGLE
    !$ nthreads = OMP_GET_NUM_THREADS ()
    !!$ WRITE (*,*)'And there are', nthreads, 'threads, from', ithread
    !$OMP FLUSH (nthreads)
    !$OMP END SINGLE
    ! c_blk_ps contains the block pointers in the product's current row
    ALLOCATE (c_blk_ps(product%nblkcols_total), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineN, "c_blk_ps",&
         product%nblkcols_total)
    !
    !$OMP DO
    DO a_row = 1, left%nblkrows_total
       a_row_log = a_row
       DBG 'A row', a_row,' with #cols=',left%row_p(a_row+1)-left%row_p(a_row)
       a_row_size = left_row_blk_size(a_row) ! stored
       c_blk_ps(:) = 0
       DO a_blk = left%row_p(a_row)+1, left%row_p(a_row+1)
          DO symmetry = 1, num_symmetries
             a_col = left%col_i(a_blk)
             IF (a_col .EQ. a_row .AND. symmetry .GT. 1) CYCLE
             a_col_size = left_col_blk_size(a_col) !stored
             a_row_log = a_row
             a_col_log = a_col
             a_row_size_logical = a_row_size
             a_col_size_logical = a_col_size
             lda = a_row_size
             ! For regular symmetric matrices
             IF (symmetry .GT. 1) THEN
                CALL swap (a_row_log, a_col_log)
                CALL swap (a_row_size_logical, a_col_size_logical)
                tr_a = .TRUE.
                a_type = 'T'
             ELSE
                tr_a = .FALSE.
                a_type = 'N'
             ENDIF
             DBG 'A col', a_col, '=>',a_col_log
             b_row = a_col_log
             c_row = a_row_log
             a_nze = a_row_size * a_col_size
             a_data => left%data(ABS(left%blk_p(a_blk))&
                  :ABS(left%blk_p(a_blk))+a_nze-1)
             b_row_size = a_col_size_logical
             ! For index-augmented matrices.
             IF (left%blk_p(a_blk) .LT. 0) THEN
                lda = a_col_size
                a_type = 'T'
             ENDIF
             DO b_blk = right%row_p(b_row)+1, right%row_p(b_row+1)
                b_col = right%col_i(b_blk)
                c_col = b_col
                DBG 'B col', b_col
                b_col_size = right_col_blk_size(b_col)
                b_nze = b_row_size * b_col_size
                b_data => right%data(&
                     ABS(right%blk_p(b_blk)):ABS(right%blk_p(b_blk))+b_nze-1)
                IF (right%blk_p(b_blk) .LT. 0) THEN
                   b_type = 'T'
                   ldb = b_col_size
                ELSE
                   b_type = 'N'
                   ldb = b_row_size
                ENDIF
                c_nze = a_row_size * b_col_size
                ! If this c column does not yet exist, "allocate" it.
                coordinate = make_coordinate_tuple (c_row, c_col)
                IF (c_blk_ps(c_col) .EQ. 0 .AND. .NOT. tr_a) THEN
                   CALL btree_find (wm%tree, coordinate, offset, block_exists)
                ELSE
                   block_exists = .TRUE.
                   offset = c_blk_ps(c_col)
                ENDIF
                IF (.NOT. block_exists) THEN
                   !.!CALL add_work_coordinate(wm, c_row, c_col, c_blk_p, error=error)
                   !$ CALL OMP_SET_LOCK (data_lock)
                   offset = c_blk_p
                   c_blk_p = c_blk_p + c_nze
                   wm%datasize = wm%datasize + c_nze
                   wm%lastblk = wm%lastblk + 1
                   !$ CALL OMP_UNSET_LOCK (data_lock)
                   CALL btree_add (wm%tree, coordinate, offset, block_exists,&
                        existing_value = offset)
                   !CALL btree_verify (index)
                   IF (block_exists) THEN
                      IF (dbg) WRITE(*,*)routineP//' EXISTS'
                      ! We can not rewind c_blk_p!
                      !$OMP ATOMIC
                      wm%datasize = wm%datasize - c_nze
                      !$OMP ATOMIC
                      wm%lastblk = wm%lastblk - 1
                   ENDIF
                   c_blk_ps(c_col) = offset
                   beta = beta_new
                ELSE
                   beta = beta_add
                ENDIF
                !offset = 1 + ithread * 4096
                !.!c_data => wm%data(c_blk_ps(c_col):c_blk_ps(c_col)+c_nze-1)
                DBG 'offset:',offset,offset+c_nze-1,'of',wm%datasize,UBOUND(wm%data),'ex',block_exists
                !% CALL OMP_SET_LOCK (data_lock)
                c_data => wm%data(offset:offset+c_nze-1)
                !% CALL OMP_UNSET_LOCK (data_lock)
                IF (dbg) WRITE(*,&
                     "(1X,A,A1,1X,A1,A,I6,A,F4.1,' (',I3,',',I3,')=(',I3,',',I3,')x(',I3,',',I3,')')")&
                                !routineN//' DGEMM ',a_type,b_type,' at',c_blk_ps(c_col),' beta=', beta,&
                     routineN//' DGEMM ',a_type,b_type,' at',offset,' beta=', beta,&
                     c_row,c_col,a_row,a_col,b_row,b_col
                CALL DGEMM(a_type, b_type,&
                     a_row_size_logical, b_col_size, a_col_size_logical,&
                     alpha,&
                     a_data, lda,&
                     b_data, ldb, beta,&
                     c_data, a_row_size_logical)
                !write(6,*)'NEWTIME',a_row_size,b_col_size,a_col_size,t2-t1,&
                !     REAL(t2-t1)/(REAL(2*a_row_size*b_col_size*a_col_size))
                lflop = lflop + a_row_size * b_col_size * a_col_size * 2.0_dp
             ENDDO ! b
          ENDDO ! symmetry
       ENDDO ! a_col
    ENDDO ! a_row

!    !$OMP END DO
    !$ CALL OMP_DESTROY_LOCK (index_lock)
    !$ CALL OMP_DESTROY_LOCK (data_lock)
    DEALLOCATE (c_blk_ps)
    !$OMP END PARALLEL
    ELSE
       DBG 'One of the matrices is empty!'
    ENDIF
    IF (dbg) CALL btree_verify (wm%tree)
    !DEALLOCATE (c_blk_ps)
    !WRITE(*,*)'NEW TOTAL TIME', nanotot, REAL(nanotot)/(2.3*REAL(10**6)),'ms'
    DBG "Post-mult blk sizes, data", wm%lastblk, wm%datasize
    !dbg = .TRUE.
    IF (dbg) THEN
       tstop = m_walltime ()
       IF (tstop-tstart .NE. 0) THEN
          DBG tstop-tstart,'s',&
               (REAL(lflop, dp)/1000000.0_dp)/(tstop-tstart),'MFLOPS'
       ENDIF
    ENDIF
    IF (PRESENT (flop)) flop = lflop
    product%w = wm
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_nn_mult


  ELEMENTAL FUNCTION make_coordinate_tuple (most, least) RESULT (tuple)
    INTEGER, INTENT(IN)                      :: most, least
    INTEGER(KIND=int_8)                      :: tuple

!tuple = IOR (ISHFT (most, 32), least)

    tuple = most
    tuple = IOR (ISHFT(tuple, 32), INT(least, int_8))
  END FUNCTION make_coordinate_tuple

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
    REAL(KIND=dp), PARAMETER                 :: alpha = 1.0_dp, &
                                                beta_add = 1.0_dp, &
                                                beta_new = 0.0_dp

    CHARACTER                                :: addee_type, sum_type
    INTEGER :: addee_blk, addee_col, addee_col_size, addee_nze, addee_row, &
      addee_row_size, ldaddee, ldsum, lflop, sum_blk, sum_blk_p, &
      sum_col_size, sum_row_size, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: addee_col_blk_size, &
                                                addee_row_blk_size
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
    addee_row_blk_size => array_data (addee%row_blk_size)
    addee_col_blk_size => array_data (addee%col_blk_size)
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
       addee_row_size = addee_row_blk_size(addee_row)
       CALL ensure_array_size(wm%data,&
            ub=wm%datasize+addee_row_size*addee%nfullcols_total &
            * MAX(1, (addee%nblkrows_total-addee_row)/4), error=error)
       DO addee_blk = addee%row_p(addee_row)+1, addee%row_p(addee_row+1)
          addee_col = addee%col_i(addee_blk)
          DBG 'Addee col', addee_col
          addee_col_size = addee_col_blk_size(addee_col)
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
                CALL DAXPY (addee_nze, 1.0_dp, addee_data, 1, sum_data, 1)
             ELSE
                CALL DAXPY (addee_nze, 1.0_dp,&
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
       IF (tstop-tstart .NE. 0) THEN
          DBG tstop-tstart,'s',&
               (REAL(lflop, dp)/1000000.0_dp)/(tstop-tstart),'MFLOPS'
       ENDIF
    ENDIF
    IF (PRESENT (flop)) flop = lflop
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_wm_add


! *****************************************************************************
!> \brief Scales a DBCSR matrix by alpha
!> \param[inout] matrix_a       DBCSR matrix
!> \param[in] alpha_scalar      (optional) a scalar
!> \param[in] alpha_vector      (optional) a vector
!> \param[in] scale_column      (optional) do we scale the columns?
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_scale(matrix_a, alpha_scalar, alpha_vector, scale_column, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_a
    REAL(dp), INTENT(IN), OPTIONAL           :: alpha_scalar
    REAL(dp), DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: alpha_vector
    LOGICAL, INTENT(IN), OPTIONAL            :: scale_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a_blk, a_col, a_col_size, &
                                                a_nze, a_row, a_row_size, &
                                                timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: my_scale_column
    REAL(dp), DIMENSION(:), POINTER          :: a_data

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    !
    !
    row_blk_size => array_data (matrix_a%row_blk_size)
    col_blk_size => array_data (matrix_a%col_blk_size)

    !
    ! set vars
    my_scale_column = .TRUE.
    IF(PRESENT(scale_column)) my_scale_column = scale_column
    !
    ! check if we have enough infos
    IF(.NOT.PRESENT(alpha_scalar).AND.&
       .NOT.PRESENT(alpha_vector).AND.&
       .NOT.PRESENT(scale_column)) THEN
       CALL stop_program("dbcsr_scale","missing an argument...")
    ENDIF
    !
    ! check if we have too many infos
    IF(PRESENT(alpha_scalar).AND.PRESENT(alpha_vector)) THEN
       CALL stop_program("dbcsr_scale","incompatible arguments...")
    ENDIF
    !
    !
    ! let's go
    DO a_row = 1, matrix_a%nblkrows_total
       a_row_size = row_blk_size(a_row)
       DO a_blk = matrix_a%row_p(a_row)+1,matrix_a%row_p(a_row+1)
          a_col = matrix_a%col_i(a_blk)
          a_col_size = col_blk_size(a_col)
          a_nze = a_row_size * a_col_size
          a_data => matrix_a%data(ABS(matrix_a%blk_p(a_blk)):ABS(matrix_a%blk_p(a_blk))+a_nze-1)
          !
          ! let's scale
          IF(PRESENT(alpha_scalar)) THEN
             !a_data = a_data * alpha_scalar
             CALL dscal(a_nze,alpha_scalar,a_data(1),1)
          ELSEIF(PRESENT(alpha_vector)) THEN
             CALL stop_program("dbcsr_scale","NYI")
             IF(my_scale_column) THEN
                !col_global = ...
                !DO col = 1,a_col_size
                !   !a_data(:,col) = a_data(:,col) * alpha_vector(col_global+col)
                !   CALL dscal(a_row_size,alpha_vector(col_global+col),a_data((col-1)*a_row_size+1),1)
                !ENDDO
             ELSE
                !row_global = ...
                !DO row = 1,a_row_size
                !   !a_data(row,:) = a_data(row,:) * alpha_vector(row_global+row)
                !   CALL dscal(a_col_size,alpha_vector(row_global+row),a_data(row),a_col_size)
                !ENDDO
             ENDIF
          ELSE
             CALL stop_program("dbcsr_scale","wrong logic")
          ENDIF
       ENDDO ! a_col
    ENDDO ! a_row
    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_scale

! *****************************************************************************
!> \brief traces a DBCSR matrix
!> \param[in] matrix_a       DBCSR matrix
!> \param[out] trace         the trace of the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_trace(matrix_a, trace, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_a
    REAL(dp), INTENT(INOUT)                  :: trace
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a_blk, a_col, a_col_size, &
                                                a_nze, a_row, a_row_size, &
                                                row, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    REAL(dp), DIMENSION(:), POINTER          :: a_data

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    row_blk_size => array_data (matrix_a%row_blk_size)
    col_blk_size => array_data (matrix_a%col_blk_size)
    !
    ! let's go
    trace = 0.0_dp
    DO a_row = 1, matrix_a%nblkrows_total
       a_row_size = row_blk_size(a_row)
       DO a_blk = matrix_a%row_p(a_row)+1,matrix_a%row_p(a_row+1)
          a_col = matrix_a%col_i(a_blk)
          IF(a_col.ne.a_row) CYCLE
          a_col_size = col_blk_size(a_col)
          IF(a_row_size.NE.a_col_size) CALL stop_program("dbcsr_trace","is that a square matrix?")
          a_nze = a_row_size**2
          a_data => matrix_a%data(ABS(matrix_a%blk_p(a_blk)):ABS(matrix_a%blk_p(a_blk))+a_nze-1)
          !
          ! let's trace the blockw
          DO row = 1,a_row_size
             trace = trace + a_data((row-1)*a_row_size+row)
          ENDDO
       ENDDO ! a_col
    ENDDO ! a_row
    !
    ! summe
    CALL mp_sum(trace,dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_trace

END MODULE dbcsr_operations

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
       dbcsr_mp_numnodes, dbcsr_mp_pgrid, dbcsr_mp_pgrid_equal
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
  USE message_passing,                 ONLY: &
       mp_allgather, mp_isendrecv, mp_rma_get, mp_sum, mp_wait, &
       mp_window_create, mp_window_destroy, mp_window_fence, mp_window_type

  !$ USE OMP_LIB
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_operations'

  PUBLIC :: dbcsr_mult_e_ec_m, dbcsr_multiply_from_left, dbcsr_trace,&
       dbcsr_scale

  PRIVATE

  INTERFACE dbcsr_trace
     MODULE PROCEDURE dbcsr_trace_a, dbcsr_trace_ab
  END INTERFACE

#define DBG IF (dbg) WRITE(*,*)routineN//" ",
#define DBGV IF (bcsr_verbose) WRITE(*,*)routineN//" ",
#define DBGI IF (info) WRITE(*,*)routineN//" ",
#define DEBUG_HEADER        LOGICAL :: dbg, info
#define DEBUG_BODY        dbg = .FALSE. ; CALL dbcsr_set_debug(dbg, info=info)
#define DEBUG_BODYY       dbg = .TRUE. ; CALL dbcsr_set_debug(dbg, info=info)

#define temp_transpose(v, r, c) RESHAPE(TRANSPOSE(RESHAPE(v,(/r,c/))),(/r*c/))


CONTAINS


! *****************************************************************************
!> \brief Multiplies two DBCSR matrices using non-RMA transfers.
!>
!> All matrices have the same processor grid distribution.
!>
!> \param[in] left_set             set of left column images
!> \param[in] right_matrix         right matrix
!> \param[out] product             DBCSR product matrix
!> \param[in,out] error            cp2k error
!> \param[out] flop                (optional) effective flop
! *****************************************************************************
  SUBROUTINE dbcsr_mult_e_ec_m (left_set, right_matrix, product_matrix,&
       error, flop)
    TYPE(dbcsr_array_type), INTENT(INOUT)    :: left_set
    TYPE(dbcsr_type), INTENT(IN)             :: right_matrix
    TYPE(dbcsr_type), INTENT(INOUT)          :: product_matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mult_e_ec_m', &
      routineP = moduleN//':'//routineN

    INTEGER :: debug_level, entries, i, img, left_dst_image, left_dst_p, &
      left_dst_pcol, left_dst_prow, left_multiplicity, left_mypcol, &
      left_myprow, left_nimages, left_npcols, left_nprows, left_p_recv, &
      left_p_send, left_src_image, left_src_p, left_src_pcol, left_src_prow, &
      mp_group, mult, mynode, mypcol, myprow, nblks, npcols, nprows, &
      numnodes, nze, right_data_rr, right_data_sr, right_dst_p, &
      right_dst_prow, right_index_rr, right_index_sr, right_mypcol, &
      right_myprow, right_npcols, right_nprows, right_p_recv, right_p_send, &
      right_prow_i, right_src_p, right_src_p_calc, right_src_prow, &
      timing_handle
    INTEGER(KIND=int_8)                      :: flop_single, flop_total
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: keys
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: left_data_rr, left_data_sr, &
                                                left_index_rr, left_index_sr
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: my_sizes, tmp_index
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: all_sizes
    INTEGER, DIMENSION(:), POINTER           :: left_index_rp, left_index_sp, &
                                                right_index_rp, right_index_sp
    INTEGER, DIMENSION(:, :), POINTER        :: left_pgrid, pgrid, right_pgrid
    REAL(KIND=dp)                            :: t_all, t_dgemm, trun, trun_t, &
                                                tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: my_waits
    REAL(KIND=dp), DIMENSION(:), POINTER     :: left_data_rp, left_data_sp, &
                                                right_data_rp, right_data_sp
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
    t_all = 0.0_dp
    t_dgemm = 0.0_dp
    IF (dbg) debug_level = 2
    IF (dbcsr_mp_nprows (dbcsr_distribution_mp (right_matrix%dist)) .NE.&
         dbcsr_mp_npcols (dbcsr_distribution_mp (right_matrix%dist))) THEN
       CALL dbcsr_mult_e_ec_m_rma (left_set, right_matrix,&
            product_matrix, error, flop=flop_total)
       IF (PRESENT (flop)) flop = flop_total
       RETURN
    ENDIF
    CALL timeset (routineN, timing_handle)
    CALL btree_new (product_matrix%w%tree)
    tstart = 0.0_dp ; tstop = 0.0_dp ; trun = 0.0_dp
    flop_total = 0
    trun_t = m_walltime ()
    left_nimages = left_set%image_dist%col_decimation
    left_multiplicity = left_set%image_dist%col_multiplicity
    CALL bcsr_assert (debug_level,&
         SIZE (left_set%mats) .EQ. left_set%image_dist%row_decimation,&
         routineP, "Declared images not equal to present images")
    mp_obj = dbcsr_distribution_mp (product_matrix%dist)
    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    npcols = dbcsr_mp_npcols (mp_obj)
    nprows = dbcsr_mp_nprows (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    pgrid => dbcsr_mp_pgrid (mp_obj)
    left_myprow = myprow
    left_mypcol = mypcol
    left_nprows = nprows
    left_npcols = npcols
    left_pgrid => pgrid
    right_myprow = myprow
    right_mypcol = mypcol
    right_nprows = nprows
    right_npcols = npcols
    right_pgrid => pgrid
    CALL bcsr_assert (debug_level,&
         left_multiplicity * npcols .EQ. left_nimages*nprows, routineP,&
         "Mulitplicity and image  mismatch with processer grid.")
    ALLOCATE (tmp_index (3, product_matrix%nblkcols_total*product_matrix%nblkrows_total))
    CALL ensure_array_size(product_matrix%w%row_i,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%col_i,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%blk_p,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    DBG 'pgrid shape', SHAPE (pgrid)
    IF (dbg .AND. mynode .EQ. 0) THEN
       WRITE(*,*) 'left images, multiplicity', left_nimages, left_multiplicity
    ENDIF
    ALLOCATE (left_set%data_windows(left_nimages))
    ALLOCATE (left_set%index_windows(left_nimages))
    ALLOCATE (my_sizes(4, left_nimages))
    my_sizes(:,:) = 0
    DO img = 1, left_nimages
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
    DO img = 1, left_nimages
       CALL mp_window_create (left_set%data_windows(img),&
            left_set%mats(img)%data, gid=mp_group)
       CALL mp_window_create (left_set%index_windows(img),&
            left_set%mats(img)%index, gid=mp_group)
    ENDDO
    ALLOCATE (lbuffer_1%mats(left_multiplicity))
    ALLOCATE (lbuffer_2%mats(left_multiplicity))
    ALLOCATE (left_data_sr(left_multiplicity), left_data_rr(left_multiplicity))
    ALLOCATE (left_index_sr(left_multiplicity), left_index_rr(left_multiplicity))
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
    ! To which processor we normally send
    right_p_send = right_pgrid (MOD (right_nprows&
         + right_myprow - 1, right_nprows),&
         right_mypcol)
    ! From which processor we normally receive.
    right_p_recv = right_pgrid (MOD (&
         right_myprow + 1, right_nprows),&
         right_mypcol)
    ! Do the actual transfers and call local multiplication
    ALLOCATE (my_waits(4, 0:right_nprows))
    my_waits(:,:) = 0
    k_index: DO right_prow_i = 0, right_nprows
       ! Right matrix
       ! right_src_prow is the matrix "source" processor (as if a
       ! direct mpi_get) would be called.
       right_src_prow = MOD (&
            right_prow_i+right_myprow + right_mypcol/left_multiplicity,&
            right_nprows)
       ! right_dst_prow is the matrix "destination" processor (i.e.,
       ! the matrix data that we currently have)
       right_dst_prow = MOD ( right_prow_i - 1 + left_multiplicity*right_nprows&
            + right_myprow + right_mypcol/left_multiplicity,&
            right_nprows)
       DBG 'RIGHT_PROW from/to', right_prow_i, right_src_prow
       !
       ! No wait in first loop iteration.
       IF (right_prow_i .GT. 0) THEN
          DBG 'It', right_prow_i,'waiting for right'
          my_waits(1, right_prow_i) = m_walltime()
          CALL mp_wait (right_data_sr)
          CALL mp_wait (right_index_sr)
          my_waits(1, right_prow_i) = m_walltime() - my_waits(1, right_prow_i)
          my_waits(2, right_prow_i) = m_walltime()
          CALL mp_wait (right_data_rr)
          CALL mp_wait (right_index_rr)
          my_waits(2, right_prow_i) = m_walltime() - my_waits(2, right_prow_i)
       ENDIF
       !
       ! No fetches in the last loop
       xfer_case_right: IF (right_prow_i .LT. right_nprows) THEN
          ! First time the meaning of "send" is different
          IF (right_prow_i .EQ. 0) THEN
             right_dst_prow = MOD (left_multiplicity*right_nprows&
                  + right_myprow - right_mypcol/left_multiplicity,&
             right_nprows)
          ENDIF
          right_src_p = right_pgrid (right_src_prow, right_mypcol)
          right_dst_p = right_pgrid (right_dst_prow, right_mypcol)
          IF (dbg) WRITE(*,'(1X,A,I4,A,"(",I3,",",I3,") sizes",I7,I7)')&
               routineN//' Getting RIGHT matrix from',right_src_p,&
               ' pgrid coor ',right_src_prow, right_mypcol,&
               all_sizes(3:4, 1, right_src_p)
          IF (dbg) WRITE(*,'(1X,A,I4,A,"(",I3,",",I3,") sizes",I7,I7)')&
               routineN//' Sending RIGHT matrix to  ',right_dst_p,&
               ' pgrid coor ',right_dst_prow, right_mypcol,&
               all_sizes(3:4, 1, right_dst_p)
          right_data_rp => rbuffer_comm%data(1:all_sizes(3,1,right_src_p))
          right_index_rp => rbuffer_comm%index(1:all_sizes(4,1,right_src_p))
          right_data_sp => rbuffer_calc%data(1:all_sizes(3,1,right_dst_p))
          right_index_sp => rbuffer_calc%index(1:all_sizes(4,1,right_dst_p))
          IF (right_prow_i .EQ. 0) THEN
             CALL mp_isendrecv (right_matrix%data, right_dst_p,&
                  right_data_rp, right_src_p, &
                  mp_group, right_data_sr, right_data_rr, tag=3)
             CALL mp_isendrecv (right_matrix%index, right_dst_p, &
                  right_index_rp,  right_src_p,&
                  mp_group, right_index_sr, right_index_rr, tag=4)
          ELSE
             CALL mp_isendrecv (right_data_sp, right_p_send,&
                  right_data_rp, right_p_recv, &
                  mp_group, right_data_sr, right_data_rr, tag=3)
             CALL mp_isendrecv (right_index_sp, right_p_send, &
                  right_index_rp,  right_p_recv,&
                  mp_group, right_index_sr, right_index_rr, tag=4)
          ENDIF
       ENDIF xfer_case_right
       ! No calculations in the first loop
       calc_case_right: IF (right_prow_i .GT. 0) THEN
          IF (all_sizes(4, 1, right_src_p_calc) .EQ. 0) THEN
             DBG 'ZERO SIZED INDEX'
          ENDIF
          CALL dbcsr_repoint_index (rbuffer_calc, error=error)
       ENDIF calc_case_right
       !
       ! No waits in first loop iteration.
       IF (right_prow_i .GT. 0) THEN
          ! Now we transfer/multiply with the left matrix
          DO mult = 1, left_multiplicity
             DBG 'It', right_prow_i,'waiting for left, multiplicity', mult
             my_waits(3, right_prow_i) = m_walltime()
             CALL mp_wait (left_data_sr(mult))
             CALL mp_wait (left_index_sr(mult))
             my_waits(3, right_prow_i) = m_walltime() - my_waits(3, right_prow_i)
             my_waits(4, right_prow_i) = m_walltime()
             CALL mp_wait (left_data_rr(mult))
             CALL mp_wait (left_index_rr(mult))
             my_waits(4, right_prow_i) = m_walltime() - my_waits(4, right_prow_i)
          ENDDO
       ENDIF
       left_src_prow = left_myprow
       left_dst_prow = left_myprow
       left_src_image = 1 + MOD (right_src_prow, left_nimages)
       left_dst_image = 1 + MOD (right_dst_prow, left_nimages)
       xfer_case_left: IF (right_prow_i .LT. nprows) THEN
          DO mult = 1, left_multiplicity
             !left_src_pcol = INT (right_src_prow / left_nimages)&
             !     + (mult-1)*(right_nprows/left_nimages)
             !left_dst_pcol = INT (right_dst_prow / left_nimages)&
             !     + (mult-1)*(right_nprows/left_nimages)
             
             !left_dst_pcol = MOD ((right_src_prow-1) / left_nimages&
             !     + (mult-1)*(right_nprows/left_nimages)&
             !     +left_npcols, left_npcols)
             left_src_pcol = MOD (&
                  right_prow_i+right_myprow + right_mypcol/left_multiplicity&
                  + right_nprows,&
                  left_npcols)
             left_dst_pcol = MOD (&
                  right_prow_i-1+right_myprow + right_mypcol/left_multiplicity&
                  + left_npcols,&
                  left_npcols)
             ! First time "dst" has a different meaning.
             IF (right_prow_i .EQ. 0) THEN
                left_dst_pcol = MOD (left_multiplicity*right_nprows&
                     - right_myprow + right_mypcol/left_multiplicity,&
                     right_nprows)
                !left_dst_pcol = INT (right_dst_prow / left_nimages)&
                !     + (mult-1)*(right_nprows/left_nimages)
             ENDIF
             !left_src_pcol = right_src_prow
             !left_dst_pcol = right_dst_prow
             left_src_p = left_pgrid (left_src_prow, left_src_pcol)
             left_dst_p = left_pgrid (left_dst_prow, left_dst_pcol)
             left_p_send = left_pgrid (left_myprow,&
                  MOD (left_mypcol - 1+left_npcols, left_npcols))
             left_p_recv = left_pgrid (left_myprow,&
                  MOD (left_mypcol + 1, left_npcols))
             IF (dbg) WRITE(*,'(1X,A,I4,A,I3,A,"(",I3,",",I3,")")')&
                  routineN//' Getting LEFT matrix from',left_src_p,&
                  ' image',left_src_image,&
                  ' pgrid coor ',left_src_prow, left_src_pcol
             IF (dbg) WRITE(*,'(1X,A,I4,A,I3,A,"(",I3,",",I3,")")')&
                  routineN//' Sending LEFT matrix to  ',left_dst_p,&
                  ' image',left_dst_image,&
                  ' pgrid coor ',left_dst_prow, left_dst_pcol
             left_data_rp => lbuffer_comm%mats(mult)%data(&
                  1:all_sizes(1, left_src_image, left_src_p))
             left_index_rp => lbuffer_comm%mats(mult)%index(&
                  1:all_sizes(2, left_src_image, left_src_p))
             left_data_sp => lbuffer_calc%mats(mult)%data(&
                  1:all_sizes(1, left_dst_image, left_dst_p))
             left_index_sp => lbuffer_calc%mats(mult)%index(&
                  1:all_sizes(2, left_dst_image, left_dst_p))
             IF (right_prow_i .EQ. 0) THEN
                CALL mp_isendrecv (left_set%mats(mult)%data, left_dst_p,&
                     left_data_rp, left_src_p, mp_group,&
                     left_data_sr(mult), left_data_rr (mult), tag=1)
                CALL mp_isendrecv (left_set%mats(mult)%index, left_dst_p,&
                     left_index_rp, left_src_p, mp_group,&
                     left_index_sr(mult), left_index_rr (mult), tag=2)
             ELSE
                CALL mp_isendrecv (left_data_sp, left_p_send,&
                     left_data_rp, left_p_recv, mp_group,&
                     left_data_sr(mult), left_data_rr (mult), tag=1)
                CALL mp_isendrecv (left_index_sp, left_p_send,&
                     left_index_rp, left_p_recv, mp_group,&
                     left_index_sr(mult), left_index_rr (mult), tag=2)
             ENDIF
          ENDDO
       ENDIF xfer_case_left
       calc_case_left: IF (right_prow_i .GT. 0) THEN
             !left_src_p_calc = left_pgrid (left_prow, right_prow-1)
          DO mult = 1, left_multiplicity
             CALL dbcsr_repoint_index (lbuffer_calc%mats(mult), error=error)
             tstart = m_walltime ()
             CALL dbcsr_nn_mult_lin(lbuffer_calc%mats(mult), rbuffer_calc,&
                  product_matrix, error, flop=flop_single,&
                  t_all=t_all, t_dgemm=t_dgemm, new_index=tmp_index)
             flop_total = flop_total + flop_single
             tstop = m_walltime ()
             IF (dbg .AND. tstop-tstart .NE. 0) THEN
                DBG 'Segment Local Multiplication time',tstop-tstart,'and',&
                     (REAL(flop_total, dp)/1000000.0_dp)/(tstop-tstart),'MFLOPS'
             ENDIF
             trun = trun + (tstop - tstart)
          ENDDO
       ENDIF calc_case_left
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
    DO img = 1, left_nimages
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
    !IF (dbg) THEN
    !   ALLOCATE (all_waits(4, 0:right_nprows, 0:numnodes-1))
    !   CALL mp_allgather (my_waits, all_waits, mp_group)
    !   logger => cp_error_get_logger (error)
    !   IF (logger%para_env%mepos == logger%para_env%source) THEN
    !      output_unit = cp_logger_get_default_unit_nr (logger)
    !      WRITE(output_unit,'(A)')"Waiting times per place per step"
    !      DO right_prow_i = 0, right_nprows
    !         WRITE(output_unit,'(4(F9.4),5X,4(F9.4))')&
    !              (MAXVAL(all_waits(i,right_prow_i,:)), i=1, 4),&
    !              (SUM(all_waits(i,right_prow_i,:))/numnodes, i=1, 4)
    !         WRITE(*,'(A,I4)')"Step ", right_prow_i
    !         DO i = 0, numnodes-1
    !            WRITE(*,'(I4,1X,4(F9.4))')i, all_waits(:, right_prow_i, i)
    !         ENDDO
    !      ENDDO
    !   ENDIF
    !ENDIF
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
    IF (PRESENT (flop)) flop = flop_total
    CALL timestop (timing_handle)
  END SUBROUTINE dbcsr_mult_e_ec_m



! *****************************************************************************
!> \brief Multiplies two DBCSR matrices with general RMA.
!>
!> All matrices have the same processor grid distribution and column images
!> in the left matrix are used to match the processor grid to the
!> right matrix rows.
!>
!> \param[in] left_set             set of left column images
!> \param[in] right_matrix         right matrix
!> \param[out] product             DBCSR product matrix
!> \param[in,out] error            cp2k error
!> \param[out] flop                (optional) effective flop
! *****************************************************************************
  SUBROUTINE dbcsr_mult_e_ec_m_rma (left_set, right_matrix, product_matrix,&
       error, flop)
    TYPE(dbcsr_array_type), INTENT(INOUT)    :: left_set
    TYPE(dbcsr_type), INTENT(IN)             :: right_matrix
    TYPE(dbcsr_type), INTENT(INOUT)          :: product_matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mult_e_ec_m_rma', &
      routineP = moduleN//':'//routineN

    INTEGER :: debug_level, entries, i, img, left_multiplicity, left_mypcol, &
      left_myprow, left_nimages, left_npcols, left_nprows, left_src_image, &
      left_src_p_calc, left_src_p_comm, left_src_pcol, left_src_prow, &
      mp_group, mult, mynode, mypcol, myprow, nblks, npcols, nprows, &
      numnodes, nze, right_mypcol, right_myprow, right_npcols, right_nprows, &
      right_prow_i, right_src_p_calc, right_src_p_comm, right_src_prow, &
      timing_handle
    INTEGER(KIND=int_8)                      :: flop_single, flop_total
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: keys
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: my_sizes, tmp_index
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: all_sizes
    INTEGER, DIMENSION(:, :), POINTER        :: left_pgrid, pgrid, right_pgrid
    REAL(KIND=dp)                            :: t_all, t_dgemm, trun, trun_t, &
                                                tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: my_waits
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
    t_all = 0.0_dp
    t_dgemm = 0.0_dp
    IF (dbg) debug_level = 2
    CALL timeset (routineN, timing_handle)
    CALL btree_new (product_matrix%w%tree)
    tstart = 0.0_dp ; tstop = 0.0_dp ; trun = 0.0_dp
    flop_total = 0
    trun_t = m_walltime ()
    left_nimages = left_set%image_dist%col_decimation
    left_multiplicity = left_set%image_dist%col_multiplicity
    CALL bcsr_assert (debug_level,&
         SIZE (left_set%mats) .EQ. left_set%image_dist%row_decimation,&
         routineP, "Declared images not equal to present images")
    mp_obj = dbcsr_distribution_mp (product_matrix%dist)
    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    npcols = dbcsr_mp_npcols (mp_obj)
    nprows = dbcsr_mp_nprows (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    pgrid => dbcsr_mp_pgrid (mp_obj)
    left_myprow = myprow
    left_mypcol = mypcol
    left_nprows = nprows
    left_npcols = npcols
    left_pgrid => pgrid
    right_myprow = myprow
    right_mypcol = mypcol
    right_nprows = nprows
    right_npcols = npcols
    right_pgrid => pgrid
    CALL bcsr_assert (debug_level,&
         left_multiplicity * npcols .EQ. left_nimages*nprows, routineP,&
         "Mulitplicity and image  mismatch with processer grid.")
    ALLOCATE (tmp_index (3, product_matrix%nblkcols_total*product_matrix%nblkrows_total))
    CALL ensure_array_size(product_matrix%w%row_i,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%col_i,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%blk_p,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    DBG 'pgrid shape', SHAPE (pgrid)
    DBG 'left images, multiplicity', left_nimages, left_multiplicity
    ALLOCATE (left_set%data_windows(left_nimages))
    ALLOCATE (left_set%index_windows(left_nimages))
    ALLOCATE (my_sizes(4, left_nimages))
    my_sizes(:,:) = 0
    DO img = 1, left_nimages
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
    DO img = 1, left_nimages
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
    ! Do the actual transfers and call local multiplication
    ALLOCATE (my_waits(4, 0:right_nprows))
    my_waits(:,:) = 0
    k_index: DO right_prow_i = 0, right_nprows
       ! Right matrix
       right_src_prow = MOD (right_prow_i+right_myprow+right_mypcol,&
            right_nprows)
       DBG 'RIGHT_PROW from/to', right_prow_i, right_src_prow
       my_waits(1, right_prow_i) = m_walltime()
       CALL mp_window_fence (right_data_window,&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
       my_waits(1, right_prow_i) = m_walltime() - my_waits(1, right_prow_i)
       my_waits(2, right_prow_i) = m_walltime()
       CALL mp_window_fence (right_index_window,&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
       my_waits(2, right_prow_i) = m_walltime() - my_waits(2, right_prow_i)
       ! No fetches in the last loop
       xfer_case_right: IF (right_prow_i .LT. right_nprows) THEN
          right_src_p_comm = right_pgrid (right_src_prow, right_mypcol)
          IF (dbg) WRITE(*,'(1X,A,"(",I3,",",I3,")")')&
               routineN//' Getting RIGHT matrix from pgrid coor ',&
               right_src_prow, right_mypcol
          IF (dbg) WRITE(*,'(1X,A,I4,A,"(",I3,",",I3,") sizes",I7,I7)')&
               routineN//' Getting RIGHT matrix from',right_src_p_comm,&
               ' pgrid coor ',right_src_prow, right_mypcol,&
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
       DO left_src_image = 1, left_nimages
          my_waits(3, right_prow_i) = m_walltime()
          !CALL dbcsr_sleep (10.0_dp)
          CALL mp_window_fence (left_set%data_windows(left_src_image),&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
          my_waits(3, right_prow_i) = m_walltime() - my_waits(3, right_prow_i)
          my_waits(4, right_prow_i) = m_walltime()
          !CALL dbcsr_sleep (10.0_dp)
          CALL mp_window_fence (left_set%index_windows(left_src_image),&
            nothing_stored=.TRUE., no_put_follows=.TRUE.)
       my_waits(4, right_prow_i) = m_walltime() - my_waits(4, right_prow_i)
       ENDDO
       left_src_prow = left_myprow
       left_src_image = 1 + MOD (right_src_prow, left_nimages)
       xfer_case_left: IF (right_prow_i .LT. nprows) THEN
          DO mult = 1, left_multiplicity
             left_src_pcol = INT (right_src_prow / left_nimages)
             IF (dbg) WRITE(*,'(1X,A,I3,A,"(",I3,",",I3,")")')&
                  routineN//' Getting LEFT matrix image',left_src_image,&
                  ' pgrid coor ',left_src_prow, left_src_pcol
             left_src_p_comm = left_pgrid (left_src_prow, left_src_pcol)
             IF (dbg) WRITE(*,'(1X,A,I4,A,I3,A,"(",I3,",",I3,")")')&
                  routineN//' Getting LEFT matrix from',left_src_p_comm,&
                  ' image',left_src_image,&
                  ' pgrid coor ',left_src_prow, left_src_pcol
             CALL mp_rma_get (left_set%data_windows(left_src_image),&
                  left_src_p_comm, 0,&
                  all_sizes(1, left_src_image, left_src_p_comm),&
                  lbuffer_comm%mats(mult)%data)
             CALL mp_rma_get (left_set%index_windows(left_src_image),&
                  left_src_p_comm, 0,&
                  all_sizes(2, left_src_image, left_src_p_comm),&
                  lbuffer_comm%mats(mult)%index)
          ENDDO
       ENDIF xfer_case_left
       calc_case_left: IF (right_prow_i .GT. 0) THEN
             !left_src_p_calc == left_pgrid (left_prow, right_prow-1)
          DO mult = 1, left_multiplicity
             CALL dbcsr_repoint_index (lbuffer_calc%mats(mult), error=error)
             tstart = m_walltime ()
             CALL dbcsr_nn_mult_lin(lbuffer_calc%mats(mult), rbuffer_calc,&
                  product_matrix, error, flop=flop_single,&
                  t_all=t_all, t_dgemm=t_dgemm, new_index=tmp_index)
             flop_total = flop_total + flop_single
             tstop = m_walltime ()
             IF (dbg .AND. tstop-tstart .NE. 0) THEN
                DBG 'Segment Local Multiplication time',tstop-tstart,'and',&
                     (REAL(flop_total, dp)/1000000.0_dp)/(tstop-tstart),'MFLOPS'
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
    DO img = 1, left_nimages
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
    !IF (dbg) THEN
    !   ALLOCATE (all_waits(4, 0:right_nprows, 0:numnodes-1))
    !   CALL mp_allgather (my_waits, all_waits, mp_group)
    !   logger => cp_error_get_logger (error)
    !   IF (logger%para_env%mepos == logger%para_env%source) THEN
    !      output_unit = cp_logger_get_default_unit_nr (logger)
    !      WRITE(output_unit,'(A)')"Waiting times per place per step"
    !      DO right_prow_i = 0, right_nprows
    !         WRITE(output_unit,'(4(F9.4),5X,4(F9.4))')&
    !              (MAXVAL(all_waits(i,right_prow_i,:)), i=1, 4),&
    !              (SUM(all_waits(i,right_prow_i,:))/numnodes, i=1, 4)
    !         WRITE(*,'(A,I4)')"Step ", right_prow_i
    !         DO i = 0, numnodes-1
    !            WRITE(*,'(I4,1X,4(F9.4))')i, all_waits(:, right_prow_i, i)
    !         ENDDO
    !      ENDDO
    !   ENDIF
    !ENDIF
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
    IF (PRESENT (flop)) flop = flop_total
    CALL timestop (timing_handle)
  END SUBROUTINE dbcsr_mult_e_ec_m_rma


! *****************************************************************************
!> \brief Multiplies two DBCSR matrices.
!> \param[in] left_set             set of left images
!> \param[in] right_matrix         right matrix
!> \param[out] product             DBCSR product workmatrix
!> \param[in,out] error            cp2k error
!> \param[out] flop                (optional) effective flop
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_from_left (left_set, right_matrix, product_matrix,&
       error, flop)
    TYPE(dbcsr_array_type), INTENT(INOUT)    :: left_set
    TYPE(dbcsr_type), INTENT(IN)             :: right_matrix
    TYPE(dbcsr_type), INTENT(INOUT)          :: product_matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_from_left', &
      routineP = moduleN//':'//routineN

    INTEGER :: debug_level, entries, i, img, left_image, left_images, &
      left_multiplicity, left_nprows, left_prow, left_src_p_calc, &
      left_src_p_comm, mp_group, mult, mynode, mypcol, myprow, nblks, npcols, &
      nprows, numnodes, nze, output_unit, right_nprows, right_prow, &
      right_prow_i, right_src_p_calc, right_src_p_comm, timing_handle
    INTEGER(KIND=int_8)                      :: flop_single, flop_total
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
    IF (dbcsr_mp_pgrid_equal (dbcsr_distribution_mp (left_set%image_dist%main),&
         dbcsr_distribution_mp (right_matrix%dist))) THEN
       CALL dbcsr_multiply_from_left_sqpg (left_set, right_matrix,&
            product_matrix, error, flop=flop_total)
       IF (PRESENT (flop)) flop = flop_total
       RETURN
    ENDIF
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
                  product_matrix, error, flop=flop_single)
             flop_total = flop_total + flop_single
             tstop = m_walltime ()
             IF (dbg .AND. tstop-tstart .NE. 0) THEN
                DBG 'Segment Local Multiplication time',tstop-tstart,'and',&
                     (REAL(flop_total, dp)/1000000.0_dp)/(tstop-tstart),'MFLOPS'
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
    IF (dbg) THEN
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
    IF (PRESENT (flop)) flop = flop_total
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
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_nn_mult', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: alpha = 1.0_dp, &
                                                beta_add = 1.0_dp, &
                                                beta_new = 0.0_dp

    CHARACTER                                :: a_type, b_type
    INTEGER :: a_blk, a_col, a_col_log, a_col_size, a_col_size_logical, &
      a_nze, a_row, a_row_log, a_row_size, a_row_size_logical, b_blk, b_col, &
      b_col_size, b_nze, b_row, b_row_size, c_blk, c_blk_p, c_col, c_nze, &
      c_row, lda, ldb, ldc, num_symmetries, stat, symmetry, timing_handle
    INTEGER(KIND=int_4)                      :: offset
    INTEGER(KIND=int_8)                      :: coordinate, lflop
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
                lflop = lflop + a_row_size * b_col_size * a_col_size * 2
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
!> \brief Multiplies two DBCSR matrices with square processor grids.
!> \param[in] left_set             set of left images
!> \param[in] right_matrix         right matrix
!> \param[out] product             DBCSR product workmatrix
!> \param[in,out] error            cp2k error
!> \param[out] flop                (optional) effective flop
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_from_left_sqpg (left_set, right_matrix,&
       product_matrix, error, flop)
    TYPE(dbcsr_array_type), INTENT(INOUT)    :: left_set
    TYPE(dbcsr_type), INTENT(IN)             :: right_matrix
    TYPE(dbcsr_type), INTENT(INOUT)          :: product_matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_multiply_from_left_sqpg', &
      routineP = moduleN//':'//routineN

    INTEGER :: debug_level, i, img, left_dst_p, left_dst_pcol, left_dst_prow, &
      left_image, left_images, left_multiplicity, left_mypcol, left_myprow, &
      left_npcols, left_nprows, left_p_recv, left_p_send, left_src_p, &
      left_src_pcol, left_src_prow, mp_group, mult, mynode, mypcol, myprow, &
      nblks, npcols, nprows, numnodes, nze, right_data_rr, right_data_sr, &
      right_dst_p, right_dst_prow, right_index_rr, right_index_sr, &
      right_mypcol, right_myprow, right_npcols, right_nprows, right_p_recv, &
      right_p_send, right_prow_i, right_src_p, right_src_prow, tag, &
      timing_handle
    INTEGER(KIND=int_8)                      :: flop_single, flop_total
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: left_data_rr, left_data_sr, &
                                                left_index_rr, left_index_sr
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: my_sizes, tmp_index
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: all_sizes
    INTEGER, DIMENSION(:), POINTER           :: left_index_rp, left_index_sp, &
                                                right_index_rp, right_index_sp
    INTEGER, DIMENSION(:, :), POINTER        :: left_pgrid, product_pgrid, &
                                                right_pgrid
    REAL(KIND=dp)                            :: t_all, t_dgemm, trun, trun_t, &
                                                tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: my_waits
    REAL(KIND=dp), DIMENSION(:), POINTER     :: left_data_rp, left_data_sp, &
                                                right_data_rp, right_data_sp
    TYPE(dbcsr_array_type), POINTER          :: lbuffer_calc, lbuffer_comm
    TYPE(dbcsr_array_type), TARGET           :: lbuffer_1, lbuffer_2
    TYPE(dbcsr_mp_obj)                       :: left_mp, mp_obj, right_mp
    TYPE(dbcsr_type), POINTER                :: rbuffer_calc, rbuffer_comm
    TYPE(dbcsr_type), TARGET                 :: rbuffer_1, rbuffer_2

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    debug_level=0
    IF (dbg) debug_level = 2
    CALL timeset (routineN, timing_handle)
    t_all = 0.0_dp
    t_dgemm = 0.0_dp
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
    mp_group = dbcsr_mp_group (mp_obj)
    left_mp = dbcsr_distribution_mp (left_set%image_dist%main)
    right_mp = dbcsr_distribution_mp (right_matrix%dist)
    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    npcols = dbcsr_mp_npcols (mp_obj)
    nprows = dbcsr_mp_nprows (mp_obj)
    product_pgrid => dbcsr_mp_pgrid (mp_obj)
    right_myprow = dbcsr_mp_myprow (right_mp)
    right_mypcol = dbcsr_mp_mypcol (right_mp)
    right_nprows = dbcsr_mp_nprows (right_mp)
    right_npcols = dbcsr_mp_npcols (right_mp)
    right_pgrid => dbcsr_mp_pgrid (right_mp)
    left_pgrid => dbcsr_mp_pgrid (left_mp)
    left_myprow = dbcsr_mp_myprow (left_mp)
    left_mypcol = dbcsr_mp_mypcol (left_mp)
    left_nprows = dbcsr_mp_nprows (left_mp)
    left_npcols = dbcsr_mp_npcols (left_mp)
    ALLOCATE (tmp_index (3, product_matrix%nblkcols_total*product_matrix%nblkrows_total))
    CALL ensure_array_size(product_matrix%w%row_i,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%col_i,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%blk_p,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    DBG 'product pgrid shape', SHAPE (dbcsr_mp_pgrid (mp_obj))
    DBG 'left pgrid shape', SHAPE (left_pgrid),';',left_nprows,left_images
    DBG 'right pgrid shape', SHAPE (right_pgrid)
    DBG 'left images, multiplicity', left_images, left_multiplicity
    CALL bcsr_assert (debug_level, left_npcols .EQ. right_nprows, routineP,&
         "Left and right matrix processor grid dimension mismatch")
    CALL bcsr_assert (debug_level,&
         left_multiplicity * npcols .EQ. left_nprows, routineP,&
         "Mulitplicity mismatch and processer grid.")
    DBG 'left pgrid =', left_pgrid
    DBG 'right pgrid =', right_pgrid
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
    !! Setup the buffer matrices
    !DO img = 1, left_images
    !   CALL mp_window_create (left_set%data_windows(img),&
    !        left_set%mats(img)%data, gid=mp_group)
    !   CALL mp_window_create (left_set%index_windows(img),&
    !        left_set%mats(img)%index, gid=mp_group)
    !ENDDO
    ALLOCATE (lbuffer_1%mats(left_multiplicity))
    ALLOCATE (lbuffer_2%mats(left_multiplicity))
    ALLOCATE (left_data_sr(left_multiplicity), left_data_rr(left_multiplicity))
    ALLOCATE (left_index_sr(left_multiplicity), left_index_rr(left_multiplicity))
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
    !CALL mp_window_create (right_data_window, right_matrix%data, gid=mp_group)
    !CALL mp_window_create (right_index_window, right_matrix%index, gid=mp_group)
    rbuffer_calc => rbuffer_1
    rbuffer_comm => rbuffer_2
    ! To which processor we normally send
    right_p_send = right_pgrid (MOD (right_nprows&
         + right_myprow - 1, right_nprows),&
         right_mypcol)
    ! From which processor we normally receive.
    right_p_recv = right_pgrid (MOD (&
         right_myprow + 1, right_nprows),&
         right_mypcol)
    !
    ! Do the actual transfers and call local multiplication
    ALLOCATE (my_waits(4, 0:right_nprows))
    my_waits(:,:) = 0    
    k_index: DO right_prow_i = 0, right_nprows
       ! Right matrix
       !
       ! right_src_prow is the matrix "source" processor (as if a
       ! direct mpi_get) would be called.
       right_src_prow = MOD ( right_prow_i &
            + right_myprow + right_mypcol/left_multiplicity, right_nprows)
       ! right_dst_prow is the matrix "destination" processor (i.e.,
       ! the matrix data that we currently have)
       right_dst_prow = MOD ( right_prow_i - 1 + 4*right_nprows&
            + right_myprow + right_mypcol/left_multiplicity, right_nprows)

       DBG 'RIGHT_PROW from/to', right_prow_i, right_src_prow, right_dst_prow,&
            myprow, right_nprows
       ! No wait in first loop iteration.
       IF (right_prow_i .GT. 0) THEN
          DBG 'It', right_prow_i,'waiting for right'
          my_waits(1, right_prow_i) = m_walltime()
          CALL mp_wait (right_data_sr)
          CALL mp_wait (right_index_sr)
          my_waits(1, right_prow_i) = m_walltime() - my_waits(1, right_prow_i)
          my_waits(2, right_prow_i) = m_walltime()
          CALL mp_wait (right_data_rr)
          CALL mp_wait (right_index_rr)
          my_waits(2, right_prow_i) = m_walltime() - my_waits(2, right_prow_i)
       ENDIF
       ! No fetches in the last loop
       xfer_case_right: IF (right_prow_i .LT. right_nprows) THEN
          ! First time the meaning of "send" is different
          IF (right_prow_i .EQ. 0) THEN
             right_dst_prow = MOD (4*right_nprows&
                  + right_myprow - right_mypcol/left_multiplicity, right_nprows)
          ENDIF
          right_src_p = right_pgrid (right_src_prow, right_mypcol)
          right_dst_p = right_pgrid (right_dst_prow, right_mypcol)
          IF (dbg) WRITE(*,'(1X,A,"(",I3,",",I3,")")')&
               routineN//' Getting RIGHT matrix pgrid coor ',&
               right_src_prow, mypcol
          IF (dbg) WRITE(*,'(1X,A,"(",I3,",",I3,")")')&
               routineN//' Sending RIGHT matrix pgrid coor ',&
               right_dst_prow, mypcol
          IF (dbg) WRITE(*,'(1X,A,I4,A,"(",I3,",",I3,") sizes",I7,I7)')&
               routineN//' Getting RIGHT matrix from',right_src_p,&
               ' pgrid coor ',right_src_prow, mypcol
          IF (dbg) WRITE(*,'(1X,A,I4,A,"(",I3,",",I3,") sizes",I7,I7)')&
               routineN//' Sending RIGHT matrix to  ',right_dst_p,&
               ' pgrid coor ',right_dst_prow, mypcol
          IF (dbg) WRITE(*,'(1X,A,2(I4,1X),A,2(I4,1X))')routineN//" Right src ",&
               right_src_p, all_sizes(3, 1, right_src_p),&
               " dst", right_dst_p,SIZE(right_matrix%data)
          right_data_rp => rbuffer_comm%data(1:all_sizes(3,1,right_src_p))
          right_index_rp => rbuffer_comm%index(1:all_sizes(4,1,right_src_p))
          right_data_sp => rbuffer_calc%data(1:all_sizes(3,1,right_dst_p))
          right_index_sp => rbuffer_calc%index(1:all_sizes(4,1,right_dst_p))
          IF (right_prow_i .EQ. 0) THEN
             CALL mp_isendrecv (right_matrix%data, right_dst_p,&
                  right_data_rp, right_src_p, &
                  mp_group, right_data_sr, right_data_rr, tag=3)
             CALL mp_isendrecv (right_matrix%index, right_dst_p, &
                  right_index_rp,  right_src_p,&
                  mp_group, right_index_sr, right_index_rr, tag=4)
          ELSE
             CALL mp_isendrecv (right_data_sp, right_p_send,&
                  right_data_rp, right_p_recv, &
                  mp_group, right_data_sr, right_data_rr, tag=3)
             CALL mp_isendrecv (right_index_sp, right_p_send, &
                  right_index_rp,  right_p_recv,&
                  mp_group, right_index_sr, right_index_rr, tag=4)
          ENDIF
       ENDIF xfer_case_right
       ! No calculations in the first loop
       calc_case_right: IF (right_prow_i .GT. 0) THEN
          CALL dbcsr_repoint_index (rbuffer_calc, error=error)
       ENDIF calc_case_right
       ! No waits in first loop iteration.
       IF (right_prow_i .GT. 0) THEN
          ! Now we transfer/multiply with the left matrix
          DO mult = 1, left_multiplicity
             DBG 'It', right_prow_i,'waiting for left, multiplicity', mult
             my_waits(3, right_prow_i) = m_walltime()
             CALL mp_wait (left_data_sr(mult))
             CALL mp_wait (left_index_sr(mult))
             my_waits(3, right_prow_i) = m_walltime() - my_waits(3, right_prow_i)
             my_waits(4, right_prow_i) = m_walltime()
             CALL mp_wait (left_data_rr(mult))
             CALL mp_wait (left_index_rr(mult))
             my_waits(4, right_prow_i) = m_walltime() - my_waits(4, right_prow_i)
          ENDDO
       ENDIF
       left_image = 1 + MOD (right_myprow, left_images)
       xfer_case_left: IF (right_prow_i .LT. right_nprows) THEN
          DO mult = 1, left_multiplicity
             left_src_prow = (myprow / left_images)
             left_src_prow = left_src_prow&
                  + (mult-1)*(left_nprows/left_multiplicity)
             left_dst_prow = (myprow / left_images)
             left_dst_prow = left_dst_prow&
                  + (mult-1)*(left_nprows/left_multiplicity)
             !left_dst_prow = left_src_prow
             !!+...left_myprow/left_multiplicity
             left_src_pcol = MOD ( right_prow_i&
                  + left_mypcol + (myprow/left_images), left_npcols)
             left_dst_pcol = MOD ( right_prow_i-1&
                  + left_mypcol + (myprow/left_images), left_npcols)
             !left_dst_pcol = MOD (-right_prow_i + 4*left_npcols&
             !     + left_mypcol - (myprow/left_images), left_npcols)
             IF (right_prow_i .EQ. 0) THEN
                left_dst_pcol = MOD (4*left_npcols&
                     + left_mypcol - (myprow/left_images), left_npcols)
             ENDIF
             !
             IF (dbg) WRITE(*,'(1X,A,I3,A,"(",I3,",",I3,")")')&
                  routineN//' Getting LEFT matrix image',left_image,&
                  ' pgrid coor ',left_src_prow, left_src_pcol
             IF (dbg) WRITE(*,'(1X,A,I3,A,"(",I3,",",I3,")")')&
                  routineN//' Sending LEFT matrix image',left_image,&
                  ' pgrid coor ',left_dst_prow, left_dst_pcol
             left_src_p = left_pgrid (left_src_prow, left_src_pcol)
             left_dst_p = left_pgrid (left_dst_prow, left_dst_pcol)
             !
             left_p_send = left_pgrid (left_src_prow,&
                  MOD (left_mypcol - 1+left_npcols, left_npcols))
             left_p_recv = left_pgrid (left_src_prow,&
                  MOD (left_mypcol + 1, left_npcols))
             !
             !
             IF (dbg) WRITE(*,'(1X,A,I4,A,I3,A,"(",I3,",",I3,")")')&
                  routineN//' Getting LEFT matrix from',left_src_p,&
                  ' image',left_image,&
                  ' pgrid coor ',left_src_prow, left_src_pcol
             IF (dbg) WRITE(*,'(1X,A,I4,A,I3,A,"(",I3,",",I3,")")')&
                  routineN//' Sending LEFT matrix to  ',left_dst_p,&
                  ' image',left_image,&
                  ' pgrid coor ',left_dst_prow, left_dst_pcol
             IF (dbg) WRITE(*,'(1X,A,I4,A,I4,A,I4)')routineN//"  Left mult. ",&
                  mult," src ", left_src_p, " dst ",right_dst_p
             left_data_rp => lbuffer_comm%mats(mult)%data(&
                  1:all_sizes(1, left_image, left_src_p))
             left_index_rp => lbuffer_comm%mats(mult)%index(&
                  1:all_sizes(2, left_image, left_src_p))
             left_data_sp => lbuffer_calc%mats(mult)%data(&
                  1:all_sizes(1, left_image, left_dst_p))
             left_index_sp => lbuffer_calc%mats(mult)%index(&
                  1:all_sizes(2, left_image, left_dst_p))
             IF (right_prow_i .EQ. 0) THEN
                CALL mp_isendrecv (left_set%mats(mult)%data, left_dst_p,&
                     left_data_rp, left_src_p, mp_group,&
                     left_data_sr(mult), left_data_rr (mult), tag=1)
                CALL mp_isendrecv (left_set%mats(mult)%index, left_dst_p,&
                     left_index_rp, left_src_p, mp_group,&
                     left_index_sr(mult), left_index_rr (mult), tag=2)
             ELSE
                CALL mp_isendrecv (left_data_sp, left_p_send,&
                     left_data_rp, left_p_recv, mp_group,&
                     left_data_sr(mult), left_data_rr (mult), tag=1)
                CALL mp_isendrecv (left_index_sp, left_p_send,&
                     left_index_rp, left_p_recv, mp_group,&
                     left_index_sr(mult), left_index_rr (mult), tag=2)
             ENDIF
          ENDDO
       ENDIF xfer_case_left
       calc_case_left: IF (right_prow_i .GT. 0) THEN
          DO mult = 1, left_multiplicity
             CALL dbcsr_repoint_index (lbuffer_calc%mats(mult), error=error)
             tstart = m_walltime ()
             CALL dbcsr_nn_mult_lin(lbuffer_calc%mats(mult), rbuffer_calc,&
                  product_matrix, error, flop=flop_single, new_index=tmp_index,&
                  t_all=t_all, t_dgemm=t_dgemm)
             flop_total = flop_total + flop_single
             tstop = m_walltime ()
             IF (tstop-tstart .NE. 0) THEN
                DBG 'Segment Local Multiplication time',tstop-tstart,'and',&
                     (REAL(flop_total, dp)/1000000.0_dp)/(tstop-tstart),'MFLOPS'
             ENDIF
             trun = trun + (tstop - tstart)
          ENDDO
       ENDIF calc_case_left
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
    !DO img = 1, left_images
    !   CALL mp_window_destroy (left_set%data_windows(img))
    !   CALL mp_window_destroy (left_set%index_windows(img))
    !ENDDO
    DEALLOCATE (left_set%data_windows)
    DEALLOCATE (left_set%index_windows)
    !CALL mp_window_destroy (right_data_window)
    !CALL mp_window_destroy (right_index_window)
    !
    CALL dbcsr_destroy (rbuffer_1, error)
    CALL dbcsr_destroy (rbuffer_2, error)
    CALL dbcsr_destroy_array (lbuffer_1, error)
    CALL dbcsr_destroy_array (lbuffer_2, error)
    DEALLOCATE (my_sizes)
    !
    !ALLOCATE (all_waits(4, 0:right_nprows, 0:numnodes-1))
    !CALL mp_allgather (my_waits, all_waits, mp_group)
    !logger => cp_error_get_logger (error)
    !IF (logger%para_env%mepos == logger%para_env%source) THEN
    !   output_unit = cp_logger_get_default_unit_nr (logger)
    !   WRITE(output_unit,'(A)')"Waiting times per place per step"
    !   DO right_prow_i = 0, right_nprows
    !      WRITE(output_unit,'(4(F9.4),5X,4(F9.4))')&
    !           (MAXVAL(all_waits(i,right_prow_i,:)), i=1, 4),&
    !           (SUM(all_waits(i,right_prow_i,:))/numnodes, i=1, 4)
    !      WRITE(*,'(A,I4)')"Step ", right_prow_i
    !      DO i = 0, numnodes-1
    !         WRITE(*,'(I4,1X,4(F9.4))')i, all_waits(:, right_prow_i, i)
    !      ENDDO
    !   ENDDO
    !ENDIF
    !
    !entries = btree_get_entries (product_matrix%w%tree)
    !CALL ensure_array_size (product_matrix%w%blk_p, ub=entries, error=error)
    !ALLOCATE (keys (entries))
    !CALL btree_delete (product_matrix%w%tree, keys, product_matrix%w%blk_p)
    !CALL ensure_array_size(product_matrix%w%row_i, ub=entries, error=error)
    !CALL ensure_array_size(product_matrix%w%col_i, ub=entries, error=error)
    !CALL dbcsr_unpack_i8_2i4 (keys, product_matrix%w%row_i,&
    !     product_matrix%w%col_i)
    !DEALLOCATE (keys)
    !CALL mp_sync (mp_group)
    CALL dbcsr_finalize(product_matrix, .TRUE., .FALSE., error)
    !
    IF (PRESENT (flop)) flop = flop_total
    CALL timestop (timing_handle)
  END SUBROUTINE dbcsr_multiply_from_left_sqpg



! *****************************************************************************
!> \brief Multiplies two DBCSR matrices, using the current local data.
!>
!> Attempts to block consecutive columns in rows of the right matrix.
!> \param[in] left, right     left and right DBCSR matrices
!> \param[in,out] product     resulting DBCSR product matrix
!> \param[in,out] wm          DBCSR product working matrix
!> \param[in,out] error       cp2k error
!> \param[out] flop           (optional) number of effective double-precision
!>                            floating point operations performed
! *****************************************************************************
  SUBROUTINE dbcsr_nn_mult_lin(left, right, product, error, flop, new_index,&
       t_all, t_dgemm)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop
    INTEGER, DIMENSION(:, :), INTENT(INOUT)  :: new_index
    REAL(KIND=dp), INTENT(INOUT)             :: t_all, t_dgemm

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_nn_mult_lin', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: alpha = 1.0_dp, &
                                                beta_add = 1.0_dp, &
                                                beta_new = 0.0_dp

    CHARACTER                                :: a_type, b_type, use_a_type, &
                                                use_b_type
    INTEGER :: a_blk, a_col, a_col_log, a_col_size, a_col_size_logical, &
      a_nze, a_row, a_row_log, a_row_size, a_row_size_logical, b_blk, b_col, &
      b_col_size, b_nze, b_row, b_row_size, c_blk_p, c_col, c_nze, c_row, i, &
      lda, ldb, new_blk, old_blk, old_row_begin_blk, old_size, prev_b_col, &
      row_begin_blk, row_blk, size_c_blk_ps, stat, timing_handle, &
      use_a_col_size_logical, use_a_offset, use_a_row_size_logical, &
      use_b_col_size, use_b_offset, use_c_offset, use_lda, use_ldb
    INTEGER(KIND=int_4)                      :: offset
    INTEGER(KIND=int_8)                      :: coordinate, lflop
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: c_blk_ps, c_blk_ps_p
    INTEGER, DIMENSION(:), POINTER           :: left_col_blk_size, &
                                                left_row_blk_size, &
                                                right_col_blk_size, &
                                                right_row_blk_size
    LOGICAL                                  :: block_exists, &
                                                start_fresh_extent, tr_a, &
                                                use_accumulated_data
    REAL(KINd=dp)                            :: beta, epoch, use_alpha, &
                                                use_beta

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    epoch = m_walltime()
    DBG "Doing NN matrix multiplication",left%name,' x ',right%name
    DBG "Pre-mult blk sizes, data", product%w%lastblk, product%w%datasize
    lflop = 0
    left_row_blk_size => array_data (left%row_blk_size)
    left_col_blk_size => array_data (left%col_blk_size)
    right_row_blk_size => array_data (right%row_blk_size)
    right_col_blk_size => array_data (right%col_blk_size)
    c_blk_p = product%w%datasize + 1
    DBG 'A rows:',left%nblkrows_total
    old_blk = 1
    old_size = product%w%lastblk
    new_blk = 0
    IF (left%nblks .GT. 0 .AND. right%nblks .GT. 0) THEN
    CALL ensure_array_size(product%w%data,&
         ub=MAX(product%w%datasize,left%nfullrows_total*right%nfullcols_total), error=error)
    ! c_blk_ps contains the block pointers in the product's current row
    ALLOCATE (c_blk_ps(product%nblkcols_total), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineN, "c_blk_ps",&
         product%nblkcols_total)
    c_blk_ps(:) = 0
    ALLOCATE (c_blk_ps_p(product%nblkcols_total), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineN, "c_blk_ps_p",&
         product%nblkcols_total)
    size_c_blk_ps = 0
    !
    use_accumulated_data = .FALSE.
    !
    DO a_row = 1, left%nblkrows_total
       a_row_log = a_row
       DBG '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
       DBG 'A row', a_row,' with #cols=',left%row_p(a_row+1)-left%row_p(a_row)
       a_row_size = left_row_blk_size(a_row) ! stored
       !c_blk_ps(:) = 0
       row_begin_blk = new_blk+1
       IF (old_blk .LE. old_size) THEN
          DO WHILE (product%w%row_i(old_blk) .LT. a_row)
             old_blk = old_blk+1
             IF (old_blk .GT. old_size) THEN
                old_blk = old_blk-1
                EXIT
             ENDIF
          ENDDO
       ENDIF
       old_row_begin_blk = old_blk
       ! Clear the c_blk_ps array
       FORALL (i = 1 : size_c_blk_ps)
          c_blk_ps(c_blk_ps_p(i)) = 0
       END FORALL
       size_c_blk_ps = 0
       ! Fill in the c_blk_ps array
       IF (old_blk .LE. old_size) THEN
          ! Old index
          DO WHILE (product%w%row_i(old_blk) .LE. a_row) ! EQ
             DBG "Adding c_blk_ps column", product%w%col_i(old_blk)
             c_blk_ps(product%w%col_i(old_blk)) = product%w%blk_p(old_blk)
             size_c_blk_ps = size_c_blk_ps + 1
             c_blk_ps_p(size_c_blk_ps) = product%w%col_i(old_blk)
             old_blk = old_blk+1
             IF (old_blk .GT. old_size) THEN
                old_blk = old_blk-1
                EXIT
             ENDIF
          ENDDO
       ENDIF
       DO a_blk = left%row_p(a_row)+1, left%row_p(a_row+1)
          row_blk = row_begin_blk
          old_blk = old_row_begin_blk
          a_col = left%col_i(a_blk)
          a_col_size = left_col_blk_size(a_col) !stored
          a_row_log = a_row
          a_col_log = a_col
          a_row_size_logical = a_row_size
          a_col_size_logical = a_col_size
          lda = a_row_size
          tr_a = .FALSE.
          a_type = 'N'
          DBG '---------------'
          DBG 'A col', a_col, '=>',a_col_log
          b_row = a_col_log
          c_row = a_row_log
          a_nze = a_row_size * a_col_size
          !a_data => left%data(ABS(left%blk_p(a_blk))&
          !     :ABS(left%blk_p(a_blk))+a_nze-1)
          b_row_size = a_col_size_logical
          ! Handle transposed blocks.
          IF (left%blk_p(a_blk) .LT. 0) THEN
             lda = a_col_size
             a_type = 'T'
          ENDIF
          ldb = b_row_size
          start_fresh_extent = .TRUE.
          DO b_blk = right%row_p(b_row)+1, right%row_p(b_row+1)
             b_col = right%col_i(b_blk)
             c_col = b_col
             DBG 'B col', b_col
             b_col_size = right_col_blk_size(b_col)
             b_nze = b_row_size * b_col_size
             !b_data => right%data(&
             !     ABS(right%blk_p(b_blk)):ABS(right%blk_p(b_blk))+b_nze-1)
             IF (right%blk_p(b_blk) .LT. 0) THEN
                b_type = 'T'
                ldb = b_col_size
             ELSE
                b_type = 'N'
                ldb = b_row_size
             ENDIF
             c_nze = a_row_size * b_col_size
             offset = c_blk_ps(c_col)
             block_exists = c_blk_ps(c_col) .GT. 0
             DBG "Found old?", block_exists!, product%w%col_i(old_blk), c_col, old_blk
             ! Try to find the current coordinate in the current row index
             IF (.NOT. block_exists) THEN
                new_blk = new_blk+1
                new_index(1, new_blk) = c_row
                new_index(2, new_blk) = c_col
                new_index(3, new_blk) = c_blk_p
                DBG 'Adding coordinate', c_row, c_col, c_blk_p
                offset = c_blk_p
                c_blk_ps (c_col) = offset
                size_c_blk_ps = size_c_blk_ps + 1
                c_blk_ps_p(size_c_blk_ps) = c_col
                c_blk_p = c_blk_p + c_nze
                product%w%datasize = product%w%datasize + c_nze
                product%w%lastblk = product%w%lastblk + 1
                beta = beta_new
             ELSE
                beta = beta_add
             ENDIF
             !DBG 'offset:',offset,offset+c_nze-1,'of',product%w%datasize,UBOUND(product%w%data),'ex',block_exists
             !c_data => product%w%data(offset:offset+c_nze-1)
             !IF (dbg) WRITE(*,&
             !     "(1X,A,A1,1X,A1,A,I6,A,F4.1,' (',I3,',',I3,')=(',I3,',',I3,')x(',I3,',',I3,')')")&
             !                   !routineN//' DGEMM ',a_type,b_type,' at',c_blk_ps(c_col),' beta=', beta,&
             !     routineN//' DGEMM ',a_type,b_type,' at',offset,' beta=', beta,&
             !     c_row,c_col,a_row,a_col,b_row,b_col
             !
             ! If we start afresh, then set the variables.
             new_extent: IF (start_fresh_extent) THEN
                use_a_type = a_type
                use_b_type = b_type
                use_a_row_size_logical = a_row_size_logical
                use_b_col_size = b_col_size
                use_a_col_size_logical = a_col_size_logical
                use_alpha = alpha
                use_a_offset = ABS (left%blk_p(a_blk))
                use_lda = lda
                use_b_offset = ABS (right%blk_p(b_blk))
                use_ldb = ldb
                use_beta = beta
                use_c_offset = offset
                prev_b_col = b_col
                start_fresh_extent = .FALSE.
             ELSE
                use_accumulated_data = .TRUE.
                ! Check to see whether the right matrix row extent can
                ! be expanded: is the previous column of the extent
                ! directly before the current one?
                IF (b_blk .LE. right%row_p(b_row+1) &
                     .AND. prev_b_col .EQ. b_col-1 .AND. b_type .EQ. 'N') THEN
                   use_accumulated_data = .FALSE.
                   prev_b_col = b_col
                ENDIF
                use_b_col_size = use_b_col_size + b_col_size
             ENDIF new_extent
             end_of_row: IF (b_blk .EQ. right%row_p(b_row+1)) THEN
                use_accumulated_data = .TRUE.
             ENDIF end_of_row
             flush_dgemm: IF (use_accumulated_data) THEN
                t_dgemm = t_dgemm - (m_walltime() - epoch)
                CALL DGEMM(use_a_type, use_b_type,&
                     use_a_row_size_logical, use_b_col_size,&
                     use_a_col_size_logical,&
                     use_alpha,&
                     left%data(use_a_offset), use_lda,&
                     right%data(use_b_offset), use_ldb, use_beta,&
                     product%w%data(use_c_offset), use_a_row_size_logical)
                !CALL DGEMM(a_type, b_type,&
                !     a_row_size_logical, b_col_size, a_col_size_logical,&
                !     alpha,&
                !     left%data(ABS(left%blk_p(a_blk))), lda,&
                !     right%data(ABS(right%blk_p(b_blk))), ldb, beta,&
                !     product%w%data(offset), a_row_size_logical)
                t_dgemm = t_dgemm + (m_walltime() - epoch)
                lflop = lflop + INT(use_a_row_size_logical, int_8)&
                     * INT(use_b_col_size, int_8)&
                     * INT(use_a_col_size_logical, int_8) * INT(2, int_8)
                use_accumulated_data = .FALSE.
                use_b_col_size = 0
                start_fresh_extent = .TRUE.
             ENDIF flush_dgemm
          ENDDO ! b
       ENDDO ! a_col
    ENDDO ! a_row
    IF (.NOT. start_fresh_extent) THEN
       t_dgemm = t_dgemm - (m_walltime() - epoch)
       CALL DGEMM(use_a_type, use_b_type,&
            use_a_row_size_logical, use_b_col_size,&
            use_a_col_size_logical,&
            use_alpha,&
            left%data(use_a_offset), use_lda,&
            right%data(use_b_offset), use_ldb, use_beta,&
            product%w%data(use_c_offset), use_a_row_size_logical)
       !CALL DGEMM(a_type, b_type,&
       !     a_row_size_logical, b_col_size, a_col_size_logical,&
       !     alpha,&
       !     left%data(ABS(left%blk_p(a_blk))), lda,&
       !     right%data(ABS(right%blk_p(b_blk))), ldb, beta,&
       !     product%w%data(offset), a_row_size_logical)
       t_dgemm = t_dgemm + (m_walltime() - epoch)
       lflop = lflop + a_row_size * b_col_size * a_col_size * 2
       use_accumulated_data = .FALSE.
    ENDIF

    DEALLOCATE (c_blk_ps, c_blk_ps_p)
    ELSE
       DBG 'One of the matrices is empty!'
    ENDIF
    ! Now merge the two indices
    CALL merge_arrays (product%w%row_i, product%w%col_i, product%w%blk_p, old_size, new_index, new_blk)
    DBG "Post-mult blk sizes, data", product%w%lastblk, product%w%datasize
    IF (PRESENT (flop)) flop = lflop
    t_all = t_all + (m_walltime() - epoch)
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_nn_mult_lin


! *****************************************************************************
!> \brief Merges two sorted indices
!> \param[in,out] old_row_i, old_col_i, old_blk_p      current index and
!>                                                     merged result
!> \param[in] add_ip                    index to add into the current index
!> \param[in] old_size                  size of current index
!> \param[in] add_size                  size of index to add into the current
!>                                      index
! *****************************************************************************
  SUBROUTINE merge_arrays (old_row_i, old_col_i, old_blk_p, old_size,&
       add_ip, add_size)
    INTEGER, DIMENSION(:), POINTER           :: old_row_i, old_col_i, &
                                                old_blk_p
    INTEGER, INTENT(IN)                      :: old_size
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: add_ip
    INTEGER, INTENT(IN)                      :: add_size

    INTEGER                                  :: add_blk, merge_from_whom, &
                                                new_blk, new_size, old_blk
    INTEGER, DIMENSION(:), POINTER           :: new_blk_p, new_col_i, &
                                                new_row_i, tmp_p

!

    new_size = old_size+add_size
    !WRITE (*,*) " Old array", old_size
    !DO i = 1, old_size
    !   WRITE(*,'(I7,2X,I7,2X,I7)')old_row_i(i),old_col_i(i),old_blk_p(i)
    !ENDDO
    !WRITE (*,*) " Add array", add_size
    !DO i = 1, add_size
    !   WRITE(*,'(I7,2X,I7,2X,I7)')add_ip (1:3, i)
    !ENDDO
    IF (add_size .GT. 0) THEN
       ALLOCATE (new_row_i (new_size))
       ALLOCATE (new_col_i (new_size))
       ALLOCATE (new_blk_p (new_size))
       old_blk = 1
       add_blk = 1
       new_blk = 1
       IF (old_size .EQ. 0) THEN
          new_row_i(1:add_size) = add_ip(1, 1:add_size)
          new_col_i(1:add_size) = add_ip(2, 1:add_size)
          new_blk_p(1:add_size) = add_ip(3, 1:add_size)
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
                add_blk = add_blk + 1
             CASE (1)
                ! Merges from the old array
                new_row_i(new_blk) = old_row_i(old_blk)
                new_col_i(new_blk) = old_col_i(old_blk)
                new_blk_p(new_blk) = old_blk_p(old_blk)
                old_blk = old_blk + 1
             CASE DEFAULT
                !WRITE(*,*)"Error: Nothing to merge"
             END SELECT
             new_blk = new_blk + 1
          ENDDO
       ENDIF
       tmp_p => old_row_i
       old_row_i => new_row_i
       DEALLOCATE (tmp_p)
       tmp_p => old_col_i
       old_col_i => new_col_i
       DEALLOCATE (tmp_p)
       tmp_p => old_blk_p
       old_blk_p => new_blk_p
       DEALLOCATE (tmp_p)
    ENDIF
    !WRITE (*,*) " New array"
    !DO i = 1, new_size
    !   WRITE(*,'(I7,2X,I7,2X,I7)')old_row_i(i),old_col_i(i),old_blk_p(i)
    !ENDDO
  END SUBROUTINE merge_arrays


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
       CALL stop_program(routineN,"missing an argument...")
    ENDIF
    !
    ! check if we have too many infos
    IF(PRESENT(alpha_scalar).AND.PRESENT(alpha_vector)) THEN
       CALL stop_program(routineN,"incompatible arguments...")
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
             CALL stop_program(routineN,"NYI")
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
             CALL stop_program(routineN,"wrong logic")
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
  SUBROUTINE dbcsr_trace_a(matrix_a, trace, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_a
    REAL(dp), INTENT(INOUT)                  :: trace
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_a', &
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
          IF(a_row_size.NE.a_col_size) CALL stop_program(routineN,"is that a square matrix?")
          a_nze = a_row_size**2
          a_data => matrix_a%data(ABS(matrix_a%blk_p(a_blk)):ABS(matrix_a%blk_p(a_blk))+a_nze-1)
          !
          ! let's trace the block
          DO row = 1,a_row_size
             trace = trace + a_data((row-1)*a_row_size+row)
          ENDDO
       ENDDO ! a_col
    ENDDO ! a_row
    !
    ! summe
    CALL mp_sum(trace,dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_trace_a

! *****************************************************************************
!> \brief traces a product of DBCSR matrices
!> \param[in] matrix_a,matrix_b  DBCSR matrices
!> \param[out] trace             the trace of the product of the matrices
!> \param[in] trans_a            (optional) is matrix_a transposed or not?
!> \param[in] trans_b            (optional) is matrix_b transposed or not?
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_trace_ab(matrix_a, matrix_b, trace, trans_a, trans_b, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_a, matrix_b
    REAL(dp), INTENT(INOUT)                  :: trace
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: trans_a, trans_b
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_ab', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1)                         :: my_trans_a, my_trans_b
    INTEGER :: a_blk, a_col, a_col_size, a_row_size, b_blk, b_col, &
      b_col_size, b_last_blk, b_row_size, nze, row, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: a_col_blk_size, &
                                                a_row_blk_size, &
                                                b_col_blk_size, b_row_blk_size
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:), POINTER          :: a_data, b_data
    REAL(dp), EXTERNAL                       :: DDOT

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    a_row_blk_size => array_data (matrix_a%row_blk_size)
    a_col_blk_size => array_data (matrix_a%col_blk_size)
    b_row_blk_size => array_data (matrix_b%row_blk_size)
    b_col_blk_size => array_data (matrix_b%col_blk_size)
    
    my_trans_a = 'T'
    IF(PRESENT(trans_a)) my_trans_a = trans_a
    my_trans_b = 'N'
    IF(PRESENT(trans_b)) my_trans_b = trans_b
    IF(my_trans_a.NE.'T'.OR.my_trans_b.NE.'N') CALL stop_program(routineN,"NYI")
    !
    ! let's go
    trace = 0.0_dp
    IF(matrix_a%nblkrows_total.NE.matrix_b%nblkrows_total) THEN
       CALL stop_program(routineN,"matrices not consistent")
    ENDIF
    DO row = 1, matrix_a%nblkrows_total
       a_row_size = a_row_blk_size(row)
       b_row_size = b_row_blk_size(row)
       IF(a_row_size.NE.b_row_size) CALL stop_program(routineN,"matrices not consistent")
       b_blk = matrix_b%row_p(row)+1
       b_last_blk = matrix_b%row_p(row+1)
       DO a_blk = matrix_a%row_p(row)+1,matrix_a%row_p(row+1)
          a_col = matrix_a%col_i(a_blk)
          a_col_size = a_col_blk_size(a_col)
          !
          ! find the b_blk we assume here that the colums are ordered !
          found = .FALSE.
          DO
             IF(b_blk.GT.b_last_blk) EXIT
             b_col = matrix_b%col_i(b_blk)
             IF(b_col.EQ.a_col) THEN
                found = .TRUE.
                EXIT
             ENDIF
             b_blk = b_blk+1
          ENDDO
          IF(.NOT.found) CYCLE
          !
          b_col_size = b_col_blk_size(b_col)
          IF(a_col_size.NE.b_col_size) CALL stop_program(routineN,"matrices not consistent")
          !
          nze = a_row_size*a_col_size
          a_data => matrix_a%data(ABS(matrix_a%blk_p(a_blk)):ABS(matrix_a%blk_p(a_blk))+nze-1)
          b_data => matrix_b%data(ABS(matrix_b%blk_p(b_blk)):ABS(matrix_b%blk_p(b_blk))+nze-1)
          !
          ! let's trace the blocks
          !trace = trace + DOT_PRODUCT(a_data,b_data)
          trace = trace + ddot(nze,a_data(1),1,b_data(1),1)
       ENDDO ! a_col
    ENDDO ! a_row
    !
    ! summe
    CALL mp_sum(trace,dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_trace_ab

END MODULE dbcsr_operations

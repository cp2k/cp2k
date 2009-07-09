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
       mp_allgather, mp_bcast, mp_isendrecv, mp_rma_get, mp_sum, mp_wait, &
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
       dbcsr_scale, dbcsr_add, dbcsr_copy, dbcsr_get_block, dbcsr_copy_vec,&
       dbcsr_put_block, dbcsr_multiply_vec, dbcsr_pack_vec, dbcsr_unpack_vec,&
       dbcsr_btriu, dbcsr_add_on_diag, packed_vec_scale, packed_vec_bif_tech,&
       block_add_on_diag, block_set, block_chol_inv

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

  INTEGER, PARAMETER, PRIVATE :: rpslot_owner = 1
  INTEGER, PARAMETER, PRIVATE :: rpslot_addblks = 2
  INTEGER, PARAMETER, PRIVATE :: rpslot_addoffset = 3
  INTEGER, PARAMETER, PRIVATE :: rpslot_oldblks = 4
  INTEGER, PARAMETER, PRIVATE :: rpslot_oldoffset = 5
  INTEGER, PARAMETER, PRIVATE :: rpslot_totaloffset = 6
  INTEGER, PARAMETER, PRIVATE :: rpnslots = 6


  LOGICAL, PARAMETER, PRIVATE :: detailed_timing = .FALSE.


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

    INTEGER :: debug_level, i, img, left_dst_image, left_dst_p, &
      left_dst_pcol, left_dst_prow, left_multiplicity, left_mypcol, &
      left_myprow, left_nimages, left_npcols, left_nprows, left_p_recv, &
      left_p_send, left_src_image, left_src_p, left_src_pcol, left_src_prow, &
      mp_group, mult, mynode, mypcol, myprow, nblks, npcols, nprows, &
      numnodes, nze, output_unit, right_data_rr, right_data_sr, right_dst_p, &
      right_dst_prow, right_index_rr, right_index_sr, right_mypcol, &
      right_myprow, right_npcols, right_nprows, right_p_recv, right_p_send, &
      right_prow_i, right_src_p, right_src_p_calc, right_src_prow, &
      timing_handle
    INTEGER(KIND=int_8)                      :: flop_single, flop_total
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: left_data_rr, left_data_sr, &
                                                left_index_rr, left_index_sr
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: my_sizes
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: all_sizes
    INTEGER, DIMENSION(:), POINTER           :: left_index_rp, left_index_sp, &
                                                right_index_rp, right_index_sp
    INTEGER, DIMENSION(:, :), POINTER        :: left_pgrid, pgrid, right_pgrid
    LOGICAL                                  :: ionode
    REAL(KIND=dp)                            :: t_all, t_dgemm, trun, trun_t, &
                                                tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: my_waits
    REAL(KIND=dp), DIMENSION(:), POINTER     :: left_data_rp, left_data_sp, &
                                                right_data_rp, right_data_sp
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
    logger => cp_error_get_logger (error)
    ionode = logger%para_env%mepos == logger%para_env%source
    output_unit = 0
    IF (ionode) output_unit = cp_logger_get_default_unit_nr (logger)
    tstart = 0.0_dp ; tstop = 0.0_dp ; trun = 0.0_dp
    t_dgemm = 0.0_dp ; t_all = 0.0_dp
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
    CALL ensure_array_size(product_matrix%w%row_i,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%col_i,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%blk_p,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%blk_d,&
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
       right_src_p_calc = right_src_p
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
                  t_all=t_all, t_dgemm=t_dgemm)
             flop_total = flop_total + flop_single
             tstop = m_walltime ()
             IF (tstop-tstart .EQ. 0) tstop = tstart+0.000001_dp
             !IF (dbg .AND. tstop-tstart .NE. 0) THEN
             IF (ionode .AND. detailed_timing) THEN
                WRITE(output_unit,'(1X,A,F9.4,A,EN12.4,A,1X,EN12.4,1X,A)')&
                     "Segment Local Multiplication time ",tstop-tstart," and ",&
                     (REAL(flop_single, dp)/1000000.0_dp)/(tstop-tstart),&
                     " MFLOPS", REAL(flop_single), "FLOP"
             ENDIF
             !ENDIF
             trun = trun + (tstop - tstart)
          ENDDO
       ENDIF calc_case_left
       CALL dbcsr_switch_sets (lbuffer_calc, lbuffer_comm)
       CALL dbcsr_switch_m_ptrs (rbuffer_calc, rbuffer_comm)
    ENDDO k_index
    trun_t = m_walltime () - trun_t
    IF (trun_t .EQ. 0) trun_t = 0.000001_dp
    IF (trun .EQ. 0) trun = 0.000001_dp
    !IF (dbg) THEN
    IF (ionode .AND. detailed_timing) THEN
       WRITE(output_unit,'(1X,A,F9.4,1X,A,EN12.4,1X,A)')&
            "  Total Local Multiplication time ",trun,"and ",&
            (REAL(flop_total, dp)/1000000.0_dp)/trun,&
            "MFLOPS"
       WRITE(output_unit,'(1X,A,F9.4,1X,A,F9.4,1X,F9.4,"%")')&
            "  Total Local DGEMM time          ",t_dgemm,&
            "index time ", t_all - t_dgemm, 100.0_dp*(t_all-t_dgemm)/t_all
       WRITE(output_unit,'(1X,A,F9.4,1X,A,EN12.4,1X,A)')&
            "        Total Multiplication time ",trun_t,"and ",&
            (REAL(flop_total, dp)/1000000.0_dp)/trun_t,&
            "MFLOPS"
    ENDIF
    !ENDIF
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
    CALL dbcsr_finalize(product_matrix, .TRUE., .FALSE., error=error)
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

    INTEGER :: debug_level, i, img, left_multiplicity, left_mypcol, &
      left_myprow, left_nimages, left_npcols, left_nprows, left_src_image, &
      left_src_p_calc, left_src_p_comm, left_src_pcol, left_src_prow, &
      mp_group, mult, mynode, mypcol, myprow, nblks, npcols, nprows, &
      numnodes, nze, output_unit, right_mypcol, right_myprow, right_npcols, &
      right_nprows, right_prow_i, right_src_p_calc, right_src_p_comm, &
      right_src_prow, timing_handle
    INTEGER(KIND=int_8)                      :: flop_single, flop_total
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: my_sizes
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: all_sizes
    INTEGER, DIMENSION(:, :), POINTER        :: left_pgrid, pgrid, right_pgrid
    LOGICAL                                  :: ionode
    REAL(KIND=dp)                            :: t_all, t_dgemm, trun, trun_t, &
                                                tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: my_waits
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
    t_all = 0.0_dp
    t_dgemm = 0.0_dp
    IF (dbg) debug_level = 2
    CALL timeset (routineN, timing_handle)
    logger => cp_error_get_logger (error)
    ionode = logger%para_env%mepos == logger%para_env%source
    output_unit = 0
    IF (ionode) output_unit = cp_logger_get_default_unit_nr (logger)
    tstart = 0.0_dp ; tstop = 0.0_dp ; trun = 0.0_dp
    t_dgemm = 0.0_dp ; t_all = 0.0_dp
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
    CALL ensure_array_size(product_matrix%w%row_i,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%col_i,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%blk_p,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%blk_d,&
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
                  t_all=t_all, t_dgemm=t_dgemm)
             flop_total = flop_total + flop_single
             tstop = m_walltime ()
             IF (tstop-tstart .EQ. 0) tstop = tstart+0.000001_dp
             !IF (dbg .AND. tstop-tstart .NE. 0) THEN
             IF (ionode .AND. detailed_timing) THEN
                WRITE(output_unit,'(1X,A,F9.4,A,EN12.4,A,1X,EN12.4,1X,A)')&
                     "Segment Local Multiplication time ",tstop-tstart," and ",&
                     (REAL(flop_single, dp)/1000000.0_dp)/(tstop-tstart),&
                     " MFLOPS", REAL(flop_single), "FLOP"
             ENDIF
             !ENDIF
             trun = trun + (tstop - tstart)
          ENDDO
       ENDIF calc_case_left
       left_src_p_calc = left_src_p_comm
       CALL dbcsr_switch_sets (lbuffer_calc, lbuffer_comm)
       CALL dbcsr_switch_m_ptrs (rbuffer_calc, rbuffer_comm)
    ENDDO k_index
    trun_t = m_walltime () - trun_t
    IF (trun_t .EQ. 0) trun_t = 0.000001_dp
    IF (trun .EQ. 0) trun = 0.000001_dp
    !IF (dbg) THEN
    IF (ionode .AND. detailed_timing) THEN
       WRITE(output_unit,'(1X,A,F9.4,1X,A,EN12.4,1X,A)')&
            "  Total Local Multiplication time ",trun,"and ",&
            (REAL(flop_total, dp)/1000000.0_dp)/trun,&
            "MFLOPS"
       WRITE(output_unit,'(1X,A,F9.4,1X,A,F9.4,1X,F9.4,"%")')&
            "  Total Local DGEMM time          ",t_dgemm,&
            "index time ", t_all - t_dgemm, 100.0_dp*(t_all-t_dgemm)/t_all
       WRITE(output_unit,'(1X,A,F9.4,1X,A,EN12.4,1X,A)')&
            "        Total Multiplication time ",trun_t,"and ",&
            (REAL(flop_total, dp)/1000000.0_dp)/trun_t,&
            "MFLOPS"
    ENDIF
    !ENDIF
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
    CALL dbcsr_finalize(product_matrix, .TRUE., .FALSE., error=error)
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
    IF (dbg .AND. detailed_timing) THEN
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
    CALL dbcsr_finalize(product_matrix, .TRUE., .FALSE., error=error)
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
    IF (dbg .AND. detailed_timing) THEN
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
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: my_sizes
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
    CALL ensure_array_size(product_matrix%w%row_i,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%col_i,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%blk_p,&
         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
    CALL ensure_array_size(product_matrix%w%blk_d,&
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
                  product_matrix, error, flop=flop_single,&
                  t_all=t_all, t_dgemm=t_dgemm)
             flop_total = flop_total + flop_single
             tstop = m_walltime ()
             IF (tstop-tstart .NE. 0 .AND. detailed_timing) THEN
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
    IF (dbg .AND. detailed_timing) THEN
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
    CALL dbcsr_finalize(product_matrix, .TRUE., .FALSE., error=error)
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
  SUBROUTINE dbcsr_nn_mult_lin(left, right, product, error, flop,&
       t_all, t_dgemm)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop
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
      b_col_size, b_nze, b_row, b_row_size, c_blk_p, c_blk_pt, c_col, c_nze, &
      c_row, datasize, i, ithread, lastblk, lda, ldb, new_blk, new_row_blks, &
      nthreads, old_blk, old_row_begin_blk, old_size, prev_b_col, &
      row_begin_blk, row_blk, size_c_blk_ps, stat, timing_handle, &
      use_a_col_size_logical, use_a_offset, use_a_row_size_logical, &
      use_b_col_size, use_b_offset, use_c_offset, use_lda, use_ldb
    INTEGER(KIND=int_4)                      :: offset
    INTEGER(KIND=int_8)                      :: coordinate, lflop
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: all_nnew_blk, c_blk_ps, &
                                                c_blk_ps_p
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: mt_row_processing
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: new_index
    INTEGER, DIMENSION(:), POINTER :: left_col_blk_size, left_row_blk_size, &
      new_blk_d, new_blk_p, new_col_i, new_row_i, right_col_blk_size, &
      right_row_blk_size
    LOGICAL                                  :: block_exists, &
                                                start_fresh_extent, tr_a, &
                                                use_accumulated_data
    REAL(KIND=dp)                            :: beta, epoch, mt_t_dgemm, &
                                                use_alpha, use_beta

    !$ INTEGER(KIND=omp_lock_kind)           :: data_lock, index_lock
    INTEGER, DIMENSION(:,:), ALLOCATABLE     :: exchange

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    epoch = m_walltime()
    DBG "Doing NN matrix multiplication",left%name,' x ',right%name
    DBG "Pre-mult blk sizes, data", product%w%lastblk, product%w%datasize
    left_row_blk_size => array_data (left%row_blk_size)
    left_col_blk_size => array_data (left%col_blk_size)
    right_row_blk_size => array_data (right%row_blk_size)
    right_col_blk_size => array_data (right%col_blk_size)
    c_blk_p = product%w%datasize + 1
    DBG 'A rows:',left%nblkrows_total
    old_size = product%w%lastblk
    mt_t_dgemm = 0.0_dp
    IF (left%nblks .GT. 0 .AND. right%nblks .GT. 0) THEN
    IF (left%nfullrows_local*right%nfullcols_local .GT.&
         SIZE (product%w%data)) THEN
       CALL ensure_array_size(product%w%data,&
            ub=MAX(product%w%datasize,&
            left%nfullrows_local*right%nfullcols_local),&
            error=error)
    ENDIF
    ALLOCATE (mt_row_processing(rpnslots, left%nblkrows_total))
    !
    !$OMP PARALLEL DEFAULT (NONE) &
    !$OMP          PRIVATE(a_row, a_col, a_blk,&
    !$OMP                  b_row, b_col, b_blk,&
    !$OMP                  c_row, c_col,&
    !$OMP                  a_row_log, a_col_log,&
    !$OMP                  a_row_size_logical,&
    !$OMP                  a_col_size_logical,&
    !$OMP                  a_row_size, a_col_size, a_nze,&
    !$OMP                  b_row_size, b_col_size, b_nze, c_nze,&
    !$OMP                  lda, ldb,&
    !$OMP                  tr_a,&
    !$OMP                  beta,&
    !$OMP                  a_type, b_type,&
    !$OMP                  coordinate,&
    !$OMP                  offset,&
    !$OMP                  block_exists,&
    !$OMP                  c_blk_ps, c_blk_ps_p, size_c_blk_ps,&
    !$OMP                  new_blk, row_blk, row_begin_blk, new_row_blks,&
    !$OMP                  old_blk, old_row_begin_blk,&
    !$OMP                  start_fresh_extent, use_accumulated_data,&
    !$OMP                  use_a_type, use_b_type,&
    !$OMP                  use_a_row_size_logical, use_a_col_size_logical,&
    !$OMP                  use_b_col_size,&
    !$OMP                  use_alpha, use_beta,&
    !$OMP                  use_a_offset, use_b_offset, use_c_offset,&
    !$OMP                  use_lda, use_ldb,&
    !$OMP                  prev_b_col,&
    !$OMP                  lastblk, datasize, c_blk_pt,&
    !$OMP                  stat, i,&
    !$OMP                  ithread)&
    !$OMP          REDUCTION (+: lflop, mt_t_dgemm) &
    !$OMP          SHARED (nthreads, mt_row_processing, epoch, dbg, error,&
    !$OMP                  new_index, all_nnew_blk, exchange,&
    !$OMP                  left, right, product,&
    !$OMP                  left_row_blk_size, left_col_blk_size,&
    !$OMP                  right_row_blk_size, right_col_blk_size,&
    !$OMP                  old_size,&
    !$OMP                  c_blk_p,&
    !$OMP                  new_row_i, new_col_i, new_blk_p, new_blk_d)
    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM ()
    !$OMP SINGLE
    nthreads = 1
    !$ nthreads = OMP_GET_NUM_THREADS ()
    ALLOCATE (new_index (4, left%nfullrows_local*right%nfullcols_local,&
         0:nthreads-1))
    ALLOCATE (all_nnew_blk(0:nthreads-1))
    ALLOCATE (exchange(2, 0:nthreads-1))
    ! Bad assumption for now: Assume alway the same number of threads.
    IF (.NOT. ASSOCIATED (product%w%datas)) THEN
       ALLOCATE (product%w%datas(nthreads))
       ALLOCATE (product%w%lastblks(nthreads))
       product%w%lastblks(:) = 0
       ALLOCATE (product%w%datasizes(nthreads))
       product%w%datasizes(:) = 0
       DO i = 1, nthreads
          NULLIFY (product%w%datas(i)%r, product%w%datas(i)%c)
       ENDDO
    ENDIF
    !$OMP END SINGLE
    CALL ensure_array_size(product%w%datas(ithread+1)%r,&
         ub=MAX(product%w%datasize,&
         left%nfullrows_local*right%nfullcols_local),&
         error=error)
    ! c_blk_ps contains the block pointers in the product's current row
    ALLOCATE (c_blk_ps(product%nblkcols_total), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineN, "c_blk_ps",&
         product%nblkcols_total)
    c_blk_ps(:) = 0
    ALLOCATE (c_blk_ps_p(product%nblkcols_total), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineN, "c_blk_ps_p",&
         product%nblkcols_total)
    size_c_blk_ps = 0
    c_blk_pt = 1
    lastblk = 0
    datasize = 0
    !
    mt_t_dgemm = 0.0_dp
    lflop = 0
    old_blk = 1
    use_accumulated_data = .FALSE.
    new_blk = 0
    start_fresh_extent = .TRUE.
    !$OMP DO
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
       mt_row_processing(rpslot_oldblks, a_row) = size_c_blk_ps
       mt_row_processing(rpslot_addoffset, a_row) = new_blk+1
       new_row_blks = 0
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
             DBG "Found old?", block_exists, c_blk_ps(c_col), offset
             ! Try to find the current coordinate in the current row index
             IF (.NOT. block_exists) THEN
                new_blk = new_blk+1
                new_index(1, new_blk, ithread) = c_row
                new_index(2, new_blk, ithread) = c_col
                new_index(3, new_blk, ithread) = c_blk_pt
                new_index(4, new_blk, ithread) = ithread+1
                new_row_blks = new_row_blks + 1
                DBG 'Adding coordinate', c_row, c_col, c_blk_p
                offset = c_blk_pt
                c_blk_pt = c_blk_pt + c_nze
                !WRITE(*,'("(",I4,"): ",A,1X,I7)')ithread,"GOT offset", offset
                datasize = datasize + c_nze
                lastblk = lastblk + 1
                c_blk_ps (c_col) = offset
                size_c_blk_ps = size_c_blk_ps + 1
                c_blk_ps_p(size_c_blk_ps) = c_col
                beta = beta_new
             ELSE
                DBG 'Using existing coordinate at', offset
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
                DBG "New extent", a_row_log, b_row, b_col, ';', a_blk, b_blk
             ELSE
                use_accumulated_data = .TRUE.
                ! Check to see whether the right matrix row extent can
                ! be expanded: is the previous column of the extent
                ! directly before the current one?
                IF (b_blk .LE. right%row_p(b_row+1) &
                     .AND. prev_b_col .EQ. b_col-1 .AND. b_type .EQ. 'N') THEN
                   use_accumulated_data = .FALSE.
                   prev_b_col = b_col
                   DBG "Continuing at b_col", b_col
                ENDIF
                use_b_col_size = use_b_col_size + b_col_size
             ENDIF new_extent
             end_of_row: IF (b_blk .EQ. right%row_p(b_row+1)) THEN
                use_accumulated_data = .TRUE.
             ENDIF end_of_row
             flush_dgemm: IF (use_accumulated_data) THEN
                mt_t_dgemm = mt_t_dgemm - (m_walltime() - epoch)
                IF (dbg) WRITE(*,&
                     "(1X,A,A1,1X,A1,A,I6,A,F4.1)")&
                     routineN//' DGEMM ',use_a_type,use_b_type,&
                     ' at',use_c_offset,' beta=', use_beta
                CALL DGEMM(use_a_type, use_b_type,&
                     use_a_row_size_logical, use_b_col_size,&
                     use_a_col_size_logical,&
                     use_alpha,&
                     left%data(use_a_offset), use_lda,&
                     right%data(use_b_offset), use_ldb, use_beta,&
                     !product%w%data(use_c_offset), use_a_row_size_logical)
                     product%w%datas(ithread+1)%r(use_c_offset), use_a_row_size_logical)
                !CALL DGEMM(a_type, b_type,&
                !     a_row_size_logical, b_col_size, a_col_size_logical,&
                !     alpha,&
                !     left%data(ABS(left%blk_p(a_blk))), lda,&
                !     right%data(ABS(right%blk_p(b_blk))), ldb, beta,&
                !     product%w%data(offset), a_row_size_logical)
                mt_t_dgemm = mt_t_dgemm + (m_walltime() - epoch)
                lflop = lflop + INT(use_a_row_size_logical, int_8)&
                     * INT(use_b_col_size, int_8)&
                     * INT(use_a_col_size_logical, int_8) * INT(2, int_8)
                use_accumulated_data = .FALSE.
                use_b_col_size = 0
                start_fresh_extent = .TRUE.
             ENDIF flush_dgemm
          ENDDO ! b
       ENDDO ! a_col
       mt_row_processing(rpslot_owner, a_row) = ithread
       mt_row_processing(rpslot_addblks, a_row) = new_row_blks
    ENDDO ! a_row
    !$OMP END DO
    DEALLOCATE (c_blk_ps, c_blk_ps_p)
    ! Now merge the new indices with the old index.
    all_nnew_blk(ithread) = new_blk
    product%w%lastblks(ithread+1) = lastblk
    product%w%datasizes(ithread+1) = datasize
    !$OMP BARRIER
    
    CALL accumulate_row_offsets (mt_row_processing, exchange,&
         SIZE(mt_row_processing,2))
    ! The upper replaces this commented block.
    !!$OMP SINGLE
    !mt_row_processing(rpslot_totaloffset, 1) = 1
    !mt_row_processing(rpslot_oldoffset, 1) = 1
    !DO a_row = 2, left%nblkrows_total
    !   mt_row_processing(rpslot_totaloffset, a_row)&
    !        = mt_row_processing(rpslot_totaloffset, a_row-1)&
    !        + mt_row_processing(rpslot_addblks, a_row-1)&
    !        + mt_row_processing(rpslot_oldblks, a_row-1)
    !   mt_row_processing(rpslot_oldoffset, a_row)&
    !        = mt_row_processing(rpslot_oldoffset, a_row-1)&
    !        + mt_row_processing(rpslot_oldblks, a_row-1)
    !ENDDO
    !!$OMP END SINGLE

    !$OMP SINGLE
    product%w%lastblk = product%w%lastblk + SUM(product%w%lastblks(:))
    product%w%datasize = product%w%datasize + SUM(product%w%datasizes(:))
    !WRITE(*,*)"New total # blocks:", product%w%lastblk
    !WRITE(*,*)"New all # blocks:  ", product%w%lastblks
    !WRITE(*,*)"New total size:", product%w%datasize
    !WRITE(*,*)"New all sizes: ", product%w%datasizes
    i = old_size + SUM (all_nnew_blk)
    ALLOCATE (new_row_i(i), new_col_i(i), new_blk_p(i), new_blk_d(i))
    !$OMP END SINGLE
    !$OMP BARRIER
    CALL merge_many_arrays (new_row_i, new_col_i, new_blk_p,&
         product%w%row_i, product%w%col_i, product%w%blk_p,&
         old_size, new_index, all_nnew_blk, mt_row_processing,&
         new_blk_d, product%w%blk_d)
    !$OMP END PARALLEL
    DEALLOCATE (new_index)
    ELSE
       DBG 'One of the matrices is empty!'
    ENDIF
    DBG "Post-mult blk sizes, data", product%w%lastblk, product%w%datasize
    IF (PRESENT (flop)) flop = lflop
    t_all = t_all + (m_walltime() - epoch)
    t_dgemm = t_dgemm + mt_t_dgemm
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_nn_mult_lin


! *****************************************************************************
!> \brief Merges two arrays
!> \param[in,out] new_row_i, new_col_i, new_blk_p      merged result
!> \param[in] new_size                  size of merged index
!> \param[in,out] old_row_i, old_col_i, old_blk_p      current index
!> \param[in] old_size                  size of current index
!> \param[in] add_ip                    index to add into the current index
!> \param[in] add_size                  size of index to add into the current
!>                                      index
! *****************************************************************************
  SUBROUTINE merge_arrays (new_row_i, new_col_i, new_blk_p, new_size,&
       old_row_i, old_col_i, old_blk_p, old_size,&
       add_ip, add_size, new_blk_d, old_blk_d)
    INTEGER, INTENT(IN)                      :: new_size
    INTEGER, DIMENSION(new_size), &
      INTENT(OUT)                            :: new_blk_p, new_col_i, &
                                                new_row_i
    INTEGER, INTENT(IN)                      :: old_size
    INTEGER, DIMENSION(old_size), INTENT(IN) :: old_blk_p, old_col_i, &
                                                old_row_i
    INTEGER, INTENT(IN)                      :: add_size
    INTEGER, DIMENSION(4, add_size), &
      INTENT(IN)                             :: add_ip
    INTEGER, DIMENSION(new_size), &
      INTENT(OUT), OPTIONAL                  :: new_blk_d
    INTEGER, DIMENSION(old_size), &
      INTENT(IN), OPTIONAL                   :: old_blk_d

    INTEGER                                  :: add_blk, merge_from_whom, &
                                                new_blk, old_blk
    LOGICAL                                  :: multidata

!
!

    multidata = PRESENT (old_blk_d) .AND. PRESENT (new_blk_d)
    CALL bcsr_assert(2, old_size+add_size .EQ. new_size, "merge_arrays",&
         "Mismatch of new and old size")
    !WRITE (*,*) " Old array", old_size
    !DO i = 1, old_size
    !   WRITE(*,'(I7,2X,I7,2X,I7)')old_row_i(i),old_col_i(i),old_blk_p(i)
    !ENDDO
    !WRITE (*,*) " Add array", add_size
    !DO i = 1, add_size
    !   WRITE(*,'(I7,2X,I7,2X,I7)')add_ip (1:3, i)
    !ENDDO
    IF (add_size .GT. 0) THEN
       old_blk = 1
       add_blk = 1
       new_blk = 1
       IF (old_size .EQ. 0) THEN
          new_row_i(1:add_size) = add_ip(1, 1:add_size)
          new_col_i(1:add_size) = add_ip(2, 1:add_size)
          new_blk_p(1:add_size) = add_ip(3, 1:add_size)
          IF (multidata) new_blk_d(1:add_size) = add_ip(4, 1:add_size)
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
                IF (multidata) new_blk_d(new_blk) = add_ip(4, add_blk)
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
    !WRITE (*,*) " New array"
    !DO i = 1, new_size
    !   WRITE(*,'(4(2X,I7))')new_row_i(i),new_col_i(i),new_blk_p(i)
    !ENDDO
  END SUBROUTINE merge_arrays


! *****************************************************************************
!> \brief Merges many sorted indices
!> \param[in,out] new_row_i, new_col_i, new_blk_p      merge index
!> \param[in,out] old_row_i, old_col_i, old_blk_p      current index
!> \param[in] add_ip                    index to add into the current index
!> \param[in] old_size                  size of current index
!> \param[in] add_sizes                 size of index to add into the current
!>                                      index
!> \param[in] row_metadata              who processed a row, how many elements
!>                                      were added, and how many were in the
!>                                      old row
!> \param[in] new_blk_d, old_blk_d      (optional) new and old location of data
!>                                      blocks
! *****************************************************************************
  SUBROUTINE merge_many_arrays (&
       new_row_i, new_col_i, new_blk_p,&
       old_row_i, old_col_i, old_blk_p, old_size,&
       add_ip, add_sizes, row_metadata,&
       new_blk_d, old_blk_d)
    INTEGER, DIMENSION(:), POINTER           :: new_row_i, new_col_i, &
                                                new_blk_p, old_row_i, &
                                                old_col_i, old_blk_p
    INTEGER, INTENT(IN)                      :: old_size
    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: add_ip
    INTEGER, DIMENSION(0:), INTENT(IN)       :: add_sizes
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: row_metadata
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: new_blk_d, old_blk_d

    INTEGER                                  :: a1, a2, ablks, add_size, &
                                                ithread, n1, n2, new_size, &
                                                nrows, nthreads, o1, o2, &
                                                oblks, row
    INTEGER, DIMENSION(:), POINTER           :: tmp_p
    LOGICAL                                  :: multidata

!

    multidata = PRESENT (old_blk_d) .AND. PRESENT (new_blk_d)
    ithread = 0
    nthreads = 1
    !$ ithread = OMP_GET_THREAD_NUM ()
    !$ nthreads = OMP_GET_NUM_THREADS ()
    !IF (ithread .EQ. 0) THEN
    !   WRITE (*,*) " MANY Old array", old_size
    !   DO n1 = 1, old_size
    !      WRITE(*,'(I7,2X,I7,2X,I7)')old_row_i(n1),old_col_i(n1),old_blk_p(n1)
    !   ENDDO
    !   DO n1 = 0, nthreads-1
    !      WRITE (*,*) " MANY Add array from", n1, ":", add_sizes(n1)
    !      DO n2 = 1, add_sizes(n1)
    !         WRITE(*,'(4(2X,I7))')add_ip (1:4, n2, n1)
    !      ENDDO
    !   ENDDO
    !   WRITE(*,*)"Metadata:"
    !   DO n1 = 1, SIZE (row_metadata, 2)
    !      WRITE(*,'(6(1X,I5))') row_metadata(:, n1)
    !   ENDDO
    !ENDIF
    !!$OMP BARRIER
    nrows = SIZE(row_metadata, 2)
    add_size = SUM (add_sizes)
    new_size = old_size+add_size
    ! There's nothing to do if nothing was added
    IF (add_size .GT. 0) THEN
       ! Process my rows.
       DO row = 1, nrows
          IF (row_metadata(rpslot_owner, row) .EQ. ithread) THEN
             oblks = row_metadata(rpslot_oldblks, row)
             o1 = row_metadata(rpslot_oldoffset, row)
             o2 = o1 + oblks - 1
             ablks = row_metadata(rpslot_addblks, row)
             a1 = row_metadata(rpslot_addoffset, row)
             a2 = a1 + ablks - 1
             n1 = row_metadata(rpslot_totaloffset, row)
             n2 = n1 + oblks + ablks - 1
             IF (n2 .GE. n1) THEN
                CALL merge_arrays(new_row_i(n1:n2), new_col_i(n1:n2),&
                     new_blk_p(n1:n2), oblks + ablks,&
                     old_row_i(o1:o2), old_col_i(o1:o2),&
                     old_blk_p(o1:o2), oblks,&
                     add_ip(:,a1:a2,ithread), ablks,&
                     new_blk_d(n1:n2), old_blk_d(o1:o2))
             ENDIF
          END IF
       ENDDO
       ! Now set the pointers to the arrays with the new values.
       IF (ithread .EQ. 0) THEN
          tmp_p => old_row_i
          old_row_i => new_row_i
          DEALLOCATE (tmp_p)
          tmp_p => old_col_i
          old_col_i => new_col_i
          DEALLOCATE (tmp_p)
          tmp_p => old_blk_p
          old_blk_p => new_blk_p
          DEALLOCATE (tmp_p)
          IF (multidata) THEN
             tmp_p => old_blk_d
             old_blk_d => new_blk_d
             DEALLOCATE (tmp_p)
          ENDIF
       ENDIF
    ELSE
       IF (ithread .EQ. 0) THEN
          DEALLOCATE (new_row_i, new_col_i, new_blk_p)
          IF (multidata) DEALLOCATE (new_blk_d)
       ENDIF
    ENDIF
    !!!$OMP BARRIER
    !IF (ithread .EQ. 0) THEN
    !   WRITE (*,*) " MANY New array"
    !   DO i = 1, new_size
    !      WRITE(*,'(4(2X,I7))')old_row_i(i),old_col_i(i),&
    !           old_blk_p(i),old_blk_d(i)
    !   ENDDO
    !ENDIF
  END SUBROUTINE merge_many_arrays


! *****************************************************************************
!> \brief Accumulates sums.
!>
!> accumulates[x,i],i>1 = \sum_{j=1,i-1}terms[x,j]
!> \param[in] terms                     terms to accumulate
!> \param[out] accumulates              accumulations of terms
!> \par Must be a pointer to force pass by reference.
!> \param exchange                      for exchanging data
! *****************************************************************************
  SUBROUTINE accumulate_row_offsets (row_meta, exchange, n)
    !INTEGER, DIMENSION(:, :), INTENT(INOUT), VOLATILE :: row_meta
    !INTEGER, DIMENSION(:, :), POINTER                 :: row_meta
    INTEGER, DIMENSION(rpnslots, *)          :: row_meta
    INTEGER, DIMENSION(:, 0:), &
      INTENT(INOUT), VOLATILE                :: exchange
    INTEGER, INTENT(IN)                      :: n

    INTEGER                                  :: i, ithread, my_first, &
                                                my_last, nthreads, s_add, &
                                                s_old, segment_size

    ithread = 0 ; nthreads = 1
    !$ ithread = OMP_GET_THREAD_NUM ()
    !$ nthreads = OMP_GET_NUM_THREADS ()
    !segment_size = CEILING (n / nthreads)
    segment_size = n / nthreads
    IF (nthreads * segment_size .LT. n) segment_size = segment_size + 1
    my_first = 1 + ithread * segment_size
    my_last = MIN (n, (ithread+1) * segment_size)
    exchange(2, ithread) = SUM (row_meta(rpslot_addblks, my_first:my_last))
    exchange(1, ithread) = exchange(2, ithread)&
         + SUM (row_meta(rpslot_oldblks, my_first:my_last))
    !$OMP BARRIER
    !$OMP SINGLE
    DO i = 1, nthreads-2
       exchange(:, i) = exchange(:, i-1) + exchange(:, i)
    ENDDO
    exchange(:, 1:nthreads-1) = exchange(:, 0:nthreads-2) + 1
    exchange(:, 0) = 1
    !$OMP END SINGLE
    IF (my_first .LE. my_last) THEN
       row_meta(rpslot_totaloffset, my_first) = exchange(1, ithread)
       row_meta(rpslot_oldoffset, my_first) = exchange(2, ithread)
    ENDIF
    !IF (ithread .EQ. 0) THEN
    !   WRITE(*,*)"row_meta preassign"
    !   WRITE(*,'(6(1X,I5))')row_meta(:,1:n)
    !ENDIF
    DO i = my_first+1, my_last
       s_add = row_meta(rpslot_addblks, i-1)
       s_old = row_meta(rpslot_oldblks, i-1)
       !IF (ithread .EQ. 0) THEN
       !   WRITE(*,'(A,I3,1X,I3,1X,A,2(1X,I5),A,6(1X,I5))')'me',ithread,i,&
       !        'add, old', s_add, s_old, ';',row_meta(:,i-1)
       !ENDIF
       row_meta(rpslot_totaloffset, i) = row_meta(rpslot_totaloffset, i-1)&
            + s_add+s_old
       row_meta(rpslot_oldoffset, i) = row_meta(rpslot_oldoffset, i-1)&
            + s_old
    ENDDO
    !$OMP FLUSH (row_meta)
  END SUBROUTINE accumulate_row_offsets



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
    IF (dbg .AND. detailed_timing) THEN
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
!> \param[in] alpha_matrix      (optional) a list of matrices
!> \param[in] scale_column      (optional) do we scale the columns?
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_scale(matrix_a, alpha_scalar, alpha_vector, alpha_matrix, side, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_a
    REAL(dp), INTENT(IN), OPTIONAL           :: alpha_scalar
    REAL(dp), DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: alpha_vector, alpha_matrix
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a_blk, a_col, a_col_size, &
                                                a_nze, a_offset, a_row, &
                                                a_row_size, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: m_offset
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: right
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    !
    !
    row_blk_size => array_data (matrix_a%row_blk_size)
    col_blk_size => array_data (matrix_a%col_blk_size)
    !
    ! set vars
    right = .TRUE.
    IF(PRESENT(side)) THEN
       SELECT CASE(side)
       CASE('right');right = .TRUE.
       CASE('left' );right = .FALSE.
       CASE DEFAULT
          WRITE(*,*) 'side='//side
          CALL stop_program(routineN,"wrong side...")
       END SELECT
    ENDIF
    !
    ! check if we have enough infos
    IF(.NOT.PRESENT(alpha_scalar).AND.&
       .NOT.PRESENT(alpha_vector).AND.&
       .NOT.PRESENT(alpha_matrix).AND.&
       .NOT.PRESENT(side)) THEN
       CALL stop_program(routineN,"missing an argument...")
    ENDIF
    !
    ! check if we have too many infos
    IF(PRESENT(alpha_scalar).AND.PRESENT(alpha_vector).OR.&
       PRESENT(alpha_scalar).AND.PRESENT(alpha_matrix).OR.&
       PRESENT(alpha_vector).AND.PRESENT(alpha_matrix)) THEN
       CALL stop_program(routineN,"incompatible arguments...")
    ENDIF
    !
    !
    IF(PRESENT(alpha_matrix)) THEN
       ALLOCATE(buff(MAXVAL(row_blk_size)*MAXVAL(col_blk_size)),&
            &   m_offset(matrix_a%nblkcols_total))
       buff = 0.0_dp
       m_offset(1) = 1
       IF(right) THEN
          DO a_col = 2,matrix_a%nblkcols_total
             a_col_size = col_blk_size(a_col-1)
             m_offset(a_col) = m_offset(a_col-1) + a_col_size**2
          ENDDO
       ELSE
          DO a_row = 2,matrix_a%nblkrows_total
             a_row_size = row_blk_size(a_row-1)
             m_offset(a_row) = m_offset(a_row-1) + a_row_size**2
          ENDDO          
       ENDIF
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
          a_offset = ABS(matrix_a%blk_p(a_blk))
          !
          ! let's scale
          IF(PRESENT(alpha_scalar)) THEN
             !a_data = a_data * alpha_scalar
             CALL dscal(a_nze,alpha_scalar,matrix_a%data(a_offset),1)
          ELSEIF(PRESENT(alpha_vector)) THEN
             CALL stop_program(routineN,"NYI")
             IF(right) THEN
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
          ELSEIF(PRESENT(alpha_matrix)) THEN
             IF(right) THEN
                !A = A * alpha
                CALL dgemm('N','N',a_row_size,a_col_size,a_col_size,&
                     &     1.0_dp,matrix_a%data(a_offset),a_row_size,&
                     &            alpha_matrix(m_offset(a_col)),a_col_size,&
                     &     0.0_dp,buff,a_row_size)
                CALL dcopy(a_nze,buff,1,matrix_a%data(a_offset),1)
             ELSE
                !A = alpha * A
                CALL dgemm('N','N',a_row_size,a_col_size,a_row_size,&
                     &     1.0_dp,alpha_matrix(m_offset(a_row)),a_row_size,&
                     &            matrix_a%data(a_offset),a_row_size,&
                     &     0.0_dp,buff,a_row_size)
                CALL dcopy(a_nze,buff,1,matrix_a%data(a_offset),1)
             ENDIF
          ELSE
             CALL stop_program(routineN,"wrong logic")
          ENDIF
       ENDDO ! a_col
    ENDDO ! a_row

    IF(PRESENT(alpha_matrix)) THEN
       DEALLOCATE(buff,m_offset)
    ENDIF

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
                                                a_nze, a_row, a_row_size, i, &
                                                timing_handle
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
          DO i = 1,a_row_size
             trace = trace + a_data((i-1)*a_row_size+i)
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
    INTEGER :: a_blk, a_col, a_col_size, a_row_size, b_blk, b_col_size, &
      b_frst_blk, b_last_blk, b_row_size, nze, row, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: a_col_blk_size, &
                                                a_row_blk_size, &
                                                b_col_blk_size, b_row_blk_size
    LOGICAL                                  :: found
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
       b_frst_blk = matrix_b%row_p(row)+1
       b_last_blk = matrix_b%row_p(row+1)
       DO a_blk = matrix_a%row_p(row)+1,matrix_a%row_p(row+1)
          a_col = matrix_a%col_i(a_blk)
          a_col_size = a_col_blk_size(a_col)
          !
          ! find the b_blk we assume here that the colums are ordered !
          CALL find_column(a_col,b_frst_blk,b_last_blk,matrix_b%col_i,b_blk,found)
          IF(found) THEN
             b_col_size = b_col_blk_size(a_col)
             IF(a_col_size.NE.b_col_size) CALL stop_program(routineN,"matrices not consistent")
             !
             nze = a_row_size*a_col_size
             !
             ! let's trace the blocks
             !trace = trace + DOT_PRODUCT(a_data,b_data)
             trace = trace + ddot(nze,matrix_a%data(ABS(matrix_a%blk_p(a_blk))),1,&
                  matrix_b%data(ABS(matrix_b%blk_p(b_blk))),1)
          ENDIF
       ENDDO ! a_col
    ENDDO ! a_row
    !
    ! summe
    CALL mp_sum(trace,dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_trace_ab

  SUBROUTINE find_column(find_col,frst_blk,last_blk,col_i,blk,found)
    INTEGER, INTENT(IN)                      :: find_col, frst_blk, last_blk
    INTEGER, DIMENSION(last_blk), INTENT(IN) :: col_i
    INTEGER, INTENT(OUT)                     :: blk
    LOGICAL, INTENT(OUT)                     :: found

    INTEGER                                  :: col

    found = .FALSE.
    blk = frst_blk
    DO
       IF(blk.GT.last_blk) EXIT
       col = col_i(blk)
       IF(col.EQ.find_col) THEN
          found = .TRUE.
          EXIT
       ENDIF
       blk = blk+1
    ENDDO
  END SUBROUTINE find_column

! *****************************************************************************
!> \brief add and scale matrices 
!>    A = alpha*A + beta*B or 
!> \param[inout] matrix_a    DBCSR matrix
!> \param[in] matrix_b       DBCSR matrix
!> \param[in] alpha_scalar   (optional) 
!> \param[in] beta_scalar    (optional)
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_add(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_a
    TYPE(dbcsr_type), INTENT(IN)             :: matrix_b
    REAL(dp), INTENT(IN), OPTIONAL           :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add', &
      routineP = moduleN//':'//routineN

    INTEGER :: a_blk, a_frst_blk, a_last_blk, a_row_size, b_blk, b_col, &
      b_col_size, b_row_size, nze, row, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: a_col_blk_size, &
                                                a_row_blk_size, &
                                                b_col_blk_size, b_row_blk_size
    LOGICAL                                  :: found
    REAL(dp)                                 :: my_alpha_scalar, &
                                                my_beta_scalar

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    a_row_blk_size => array_data (matrix_a%row_blk_size)
    a_col_blk_size => array_data (matrix_a%col_blk_size)
    b_row_blk_size => array_data (matrix_b%row_blk_size)
    b_col_blk_size => array_data (matrix_b%col_blk_size)
    
    my_alpha_scalar = 1.0_dp
    IF(PRESENT(alpha_scalar)) my_alpha_scalar = alpha_scalar
    my_beta_scalar = 1.0_dp
    IF(PRESENT(beta_scalar)) my_beta_scalar = beta_scalar
    !
    ! let's go
    IF(matrix_a%nblkrows_total.NE.matrix_b%nblkrows_total) THEN
       CALL stop_program(routineN,"matrices not consistent")
    ENDIF

    CALL dbcsr_scale(matrix_a,alpha_scalar=my_alpha_scalar,error=error)

    IF(my_beta_scalar.NE.0.0_dp) THEN
       DO row = 1, matrix_b%nblkrows_total
          a_row_size = a_row_blk_size(row)
          b_row_size = b_row_blk_size(row)
          IF(a_row_size.NE.b_row_size) CALL stop_program(routineN,"matrices not consistent")
          a_blk = matrix_a%row_p(row)+1
          a_last_blk = matrix_a%row_p(row+1)
          DO b_blk = matrix_b%row_p(row)+1,matrix_b%row_p(row+1)
             b_col = matrix_b%col_i(b_blk)
             ! find the a_blk we assume here that the colums are ordered !
             a_frst_blk = a_blk
             CALL find_column(b_col,a_frst_blk,a_last_blk,matrix_a%col_i,a_blk,found)
             IF(found) THEN
                !
                b_col_size = b_col_blk_size(b_col)
                !
                nze = a_row_size*b_col_size
                !
                ! let's scale the block
                CALL daxpy(nze,my_beta_scalar,matrix_b%data(ABS(matrix_a%blk_p(a_blk))),1,&
                     &                        matrix_a%data(ABS(matrix_b%blk_p(b_blk))),1)
             ELSE
                CALL stop_program(routineN,"need to add the block")
             ENDIF
          ENDDO ! a_col
       ENDDO ! a_row
    ENDIF

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_add

! *****************************************************************************
!> \brief add a constant to the diagonal of a matrix
!> \param[inout] matrix       DBCSR matrix
!> \param[in]    alpha_scalar scalar
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_add_on_diag(matrix, alpha_scalar, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    REAL(dp), INTENT(IN)                     :: alpha_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_on_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, frst_blk, i, last_blk, &
                                                nze, offset, row, row_size, &
                                                timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:), POINTER          :: DATA

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    row_blk_size => array_data (matrix%row_blk_size)
    col_blk_size => array_data (matrix%col_blk_size)
    
    DO row = 1, matrix%nblkrows_total
       row_size = row_blk_size(row)
       blk = matrix%row_p(row)+1
       frst_blk = matrix%row_p(row)+1
       last_blk = matrix%row_p(row+1)
       CALL find_column(row,frst_blk,last_blk,matrix%col_i,blk,found)
       IF(found) THEN
          offset = ABS(matrix%blk_p(blk))
          nze = row_size**2
          DATA => matrix%data(offset:offset+nze-1)
          !
          ! let's add on diag
          DO i = 1,row_size
             DATA((i-1)*row_size+i) = DATA((i-1)*row_size+i) + alpha_scalar
          ENDDO
       ELSE
          ! need to create the block
          CALL stop_program(routineN,"need ot create the block, NYI")
       ENDIF
    ENDDO ! row

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_add_on_diag

! *****************************************************************************
!> \brief copy a matrix 
!> \param[in]    matrix_a  DBCSR matrix
!> \param[inout] matrix_b  DBCSR matrix
!> \param[in]    name      name of the new matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_copy(matrix_b, matrix_a, name, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_b
    TYPE(dbcsr_type), INTENT(IN)             :: matrix_a
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nze, timing_handle

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL dbcsr_init(matrix_b, error)
    IF(PRESENT(name)) THEN
       CALL dbcsr_create(matrix_b, name, &
            matrix_a%dist, matrix_a%matrix_type, matrix_a%row_blk_size,&
            matrix_a%col_blk_size, matrix_a%nblks, matrix_a%nze, &
            matrix_a%special_memory, error=error)
    ELSE
       CALL dbcsr_create(matrix_b,'copy of '//matrix_a%name, &
            matrix_a%dist, matrix_a%matrix_type, matrix_a%row_blk_size,&
            matrix_a%col_blk_size, matrix_a%nblks, matrix_a%nze, &
            matrix_a%special_memory, error=error)
    ENDIF
    CALL ensure_array_size(matrix_b%index, ub=SIZE(matrix_a%index), error=error)
    CALL ensure_array_size(matrix_b%data, ub=SIZE(matrix_a%data), error=error)
    !
    ! copy index and data
    matrix_b%index(:) = matrix_a%index(:)
    CALL dcopy(SIZE(matrix_a%data),matrix_a%data,1,matrix_b%data,1)
    !
    ! the row_p, col_i and blk_p ...
    CALL dbcsr_repoint_index(matrix_b, error=error)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_copy

! *****************************************************************************
!> \brief copy a vector (row or column) from a matrix to another matrix
!> \param[inout] matrix_b  matrix
!> \param[in]    matrix_a  matrix
!> \param[in]    what      row or column
!> \param[in]    ivec      the column/row index
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_copy_vec(matrix_a, matrix_b, what, ivec, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_a
    TYPE(dbcsr_type), INTENT(IN)             :: matrix_b
    CHARACTER(len=*), INTENT(in)             :: what
    INTEGER, INTENT(IN)                      :: ivec
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_vec', &
      routineP = moduleN//':'//routineN

    INTEGER :: a_blk, a_col_size, a_frst_blk, a_last_blk, a_offset, &
      a_row_size, b_blk, b_col, b_col_size, b_offset, b_row_size, col, nze, &
      row, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: a_col_blk_size, &
                                                a_row_blk_size, &
                                                b_col_blk_size, b_row_blk_size
    LOGICAL                                  :: found

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    a_row_blk_size => array_data (matrix_a%row_blk_size)
    a_col_blk_size => array_data (matrix_a%col_blk_size)
    b_row_blk_size => array_data (matrix_b%row_blk_size)
    b_col_blk_size => array_data (matrix_b%col_blk_size)
    !
    ! let's go
    SELECT CASE(what)
    CASE('column')
       col = ivec
       CALL stop_program(routineN, "NYI")
    CASE('row')
       row = ivec
       a_row_size = a_row_blk_size(row)
       b_row_size = b_row_blk_size(row)
       IF(a_row_size.NE.b_row_size) CALL stop_program(routineN, "matrices not consistent")
       a_frst_blk = matrix_a%row_p(row)+1
       a_last_blk = matrix_a%row_p(row+1)
       DO b_blk = matrix_b%row_p(row)+1,matrix_b%row_p(row+1)
          b_col = matrix_b%col_i(b_blk)
          CALL find_column(b_col, a_frst_blk, a_last_blk, matrix_a%col_i, a_blk, found)
          IF(found) THEN
             a_col_size = a_col_blk_size(b_col)
             b_col_size = b_col_blk_size(b_col)
             IF(a_col_size.NE.b_col_size) CALL stop_program(routineN, "matrices not consistent")
             nze = a_row_size * a_col_size
             !
             ! let's copy the block
             a_offset = ABS(matrix_a%blk_p(a_blk))
             b_offset = ABS(matrix_b%blk_p(b_blk))
             CALL dcopy(nze, matrix_b%data(b_offset), 1, matrix_a%data(a_offset), 1)
          ELSE
             CALL stop_program(routineN, "need ot create the block, NYI")
          ENDIF
       ENDDO ! a_col
    CASE DEFAULT
       CALL stop_program(routineN,"copy what?")
    END SELECT
    !
    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_copy_vec

! *****************************************************************************
!> \brief block triu of a dbcsr matrix
!> \param[in] matrix_b  triu of the matrix
!> \param[in] matrix_a  the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_btriu(matrix_b, matrix_a, error)

    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix_b
    TYPE(dbcsr_type), INTENT(IN)             :: matrix_a
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_btriu', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_size, nze, &
                                                offset, row, row_size, &
                                                timing_handle
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL dbcsr_copy(matrix_b, matrix_a, name="triu of "//matrix_a%name, error=error)

    row_blk_size => array_data (matrix_a%row_blk_size)
    col_blk_size => array_data (matrix_a%col_blk_size)

    !
    ! This is not the way to do it !
    DO row = 1, matrix_b%nblkrows_total
       row_size = row_blk_size(row)
       DO blk = matrix_b%row_p(row)+1,matrix_b%row_p(row+1)
          col = matrix_b%col_i(blk)
          IF(row.GT.col) THEN
             col_size = col_blk_size(col)
             offset = ABS(matrix_b%blk_p(blk))
             nze = row_size * col_size
             matrix_b%data(offset:offset+nze-1) = 0.0_dp
          ENDIF
       ENDDO
    ENDDO

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_btriu

! *****************************************************************************
!> \brief get a block in a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  block  the block to get
!> \param[in]  row    the row
!> \param[in]  col    the column
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_get_block(matrix,block,row,col,found,error)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: block
    INTEGER, INTENT(IN)                      :: row, col
    LOGICAL, INTENT(OUT)                     :: found
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col_size, frst_blk, &
                                                last_blk, nze, offset, &
                                                row_size
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size

    row_blk_size => array_data (matrix%row_blk_size)
    col_blk_size => array_data (matrix%col_blk_size)

    row_size = row_blk_size(row)
    col_size = col_blk_size(col)
    blk = matrix%row_p(row)+1
    frst_blk = matrix%row_p(row)+1
    last_blk = matrix%row_p(row+1)
    CALL find_column(col,frst_blk,last_blk,matrix%col_i,blk,found)
    IF(found) THEN
       nze = row_size*col_size
       !
       ! let's copy the block
       offset = ABS(matrix%blk_p(blk))
       CALL dcopy(nze,matrix%data(offset),1,block(1),1)
    ENDIF

  END SUBROUTINE dbcsr_get_block

! *****************************************************************************
!> \brief put a block in a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  block  the block to put
!> \param[in]  row    the row
!> \param[in]  col    the column
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_put_block(matrix, block, row, col, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    REAL(dp), DIMENSION(:), INTENT(IN)       :: block
    INTEGER, INTENT(IN)                      :: row, col
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col_size, frst_blk, &
                                                last_blk, nze, offset, &
                                                row_size
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found

! IN PARALLEL WE SHOULD ENSURE THAT THE BLOCK GOES TO THE RIGHT PLACE !

    row_blk_size => array_data (matrix%row_blk_size)
    col_blk_size => array_data (matrix%col_blk_size)

    row_size = row_blk_size(row)
    col_size = col_blk_size(col)
    blk = matrix%row_p(row)+1
    frst_blk = matrix%row_p(row)+1
    last_blk = matrix%row_p(row+1)
    CALL find_column(col,frst_blk,last_blk,matrix%col_i,blk,found)
    IF(found) THEN
       nze = row_size*col_size
       !
       ! let's copy the block
       offset = ABS(matrix%blk_p(blk))
       CALL dcopy(nze,block(1),1,matrix%data(offset),1)
    ELSE
       CALL stop_program(routineN,"need ot create the block, NYI")
    ENDIF

  END SUBROUTINE dbcsr_put_block

! *****************************************************************************
!> \brief multiply a dbcsr with a column vector
!>        c (packed column) = A (dbscr) * B_b (dbcsr)
!> \param[in]  matrix_a  DBCSR matrix
!> \param[in]  matrix_b  DBCSR matrix
!> \param[in]  b_col     the column of the matrix B
!> \param[out] pkd       the result in a packed format
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_vec(matrix_a, matrix_b, b_col, pkd, error)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix_a, matrix_b
    INTEGER, INTENT(in)                      :: b_col
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pkd
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_vec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: pkd_b

    CALL timeset(routineN, timing_handle)

    ALLOCATE(pkd_b(SIZE(pkd)))
    !
    ! packed the b_col
    CALL dbcsr_pack_vec(matrix_b, b_col, pkd_b, 'column', error)
    !
    ! send the packed col to the right guy
    !CALL packed_vec_bcast(pkd_b, source, , mp_obj, error)
    !
    ! local multiply
    CALL dbcsr_mult_pack_vec(matrix_a, pkd_b, b_col, pkd, error)
    !
    ! redistribute the result
    !CALL ...

    DEALLOCATE(pkd_b)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_multiply_vec

! *****************************************************************************
!> \brief multiply a dbcsr with a packed vector as
!>        c (packed column) = A (dbscr) * b (packed column)
!> \param[in]  matrix      DBCSR matrix
!> \param[in]  pkd_b       packed vector b
!> \param[in]  b_col       the column of the packed b (this should be changed in b_col_size)
!> \param[out] pkd_c       packed vector c
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_mult_pack_vec(matrix_a, pkd_b, b_col, pkd_c, error)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix_a
    REAL(dp), DIMENSION(:), INTENT(in)       :: pkd_b
    INTEGER, INTENT(in)                      :: b_col
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pkd_c
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mult_pack_vec', &
      routineP = moduleN//':'//routineN

    INTEGER :: a_blk, a_col, a_col_size, a_offset, a_row, a_row_size, &
      b_col_size, b_offset, c_offset, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: a_col_blk_size, a_row_blk_size
    LOGICAL                                  :: new_c

    CALL timeset(routineN, timing_handle)

    a_row_blk_size => array_data (matrix_a%row_blk_size)
    a_col_blk_size => array_data (matrix_a%col_blk_size)


    IF(pkd_b(matrix_a%nblkrows_total+1).LE.0) CALL stop_program(routineN,"1pkd_b(n+1).LE.0")

    pkd_c(:) = 0.0_dp ! if not done
    c_offset = matrix_a%nblkrows_total + 2
    b_col_size = a_col_blk_size(b_col)


    DO a_row = 1, matrix_a%nblkrows_total
       a_row_size = a_row_blk_size(a_row)
       new_c = .FALSE.
       DO a_blk = matrix_a%row_p(a_row)+1,matrix_a%row_p(a_row+1)
          a_col = matrix_a%col_i(a_blk)
          a_offset = ABS(matrix_a%blk_p(a_blk))
          a_col_size = a_col_blk_size(a_col)
          b_offset = INT(pkd_b(a_col))
          IF(b_offset.LE.0) CYCLE
          !
          pkd_c(a_row) = REAL(c_offset,dp)
          CALL dgemm('N','N',a_row_size,b_col_size,a_col_size,&
               &     1.0_dp,matrix_a%data(a_offset),a_row_size,&
               &            pkd_b(b_offset),a_col_size,&
               &     1.0_dp,pkd_c(c_offset),a_row_size)
          new_c = .TRUE.
       ENDDO ! a_col
       IF(new_c) c_offset = c_offset + a_row_size*b_col_size
    ENDDO ! a_row
    pkd_c(matrix_a%nblkrows_total+1) = REAL(c_offset,dp)

    IF(pkd_b(matrix_a%nblkrows_total+1).LE.0) CALL stop_program(routineN,"2pkd_b(n+1).LE.0")
    IF(pkd_c(matrix_a%nblkrows_total+1).LE.0) CALL stop_program(routineN,"2pkd_c(n+1).LE.0")

    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_mult_pack_vec

! *****************************************************************************
!> \brief pack a row or column of a dbcsr
!> \param[in]  matrix      DBCSR matrix
!> \param[in]  ivec        the vector to pack
!> \param[out] pkd_vec     packed vector
!> \param[in]  what        what to pack (row or column)
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_pack_vec(matrix, ivec, pkd_vec, what, error)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(IN)                      :: ivec
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_vec
    CHARACTER(LEN=*), INTENT(IN)             :: what
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_pack_vec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_size, frst_blk, &
                                                last_blk, nze, row, row_size, &
                                                timing_handle, v_offset
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    row_blk_size => array_data (matrix%row_blk_size)
    col_blk_size => array_data (matrix%col_blk_size)

    pkd_vec(:) = 0.0_dp ! should be big enough to hold all the data
    !
    ! let's go
    SELECT CASE(what)
    CASE('column')
       v_offset = matrix%nblkrows_total+2 ! we need nblkrows_total+1 to store rows 
       !                                    and 1 extra for adding new data
       col = ivec
       col_size = col_blk_size(col)
       DO row = 1, matrix%nblkrows_total
          row_size = row_blk_size(row)
          blk = matrix%row_p(row)+1
          frst_blk = matrix%row_p(row)+1
          last_blk = matrix%row_p(row+1)
          CALL find_column(col,frst_blk,last_blk,matrix%col_i,blk,found)
          IF(found) THEN
             nze = row_size*col_size
             !
             ! let's copy the block
             pkd_vec(row) = REAL(v_offset,dp)
             CALL dcopy(nze,matrix%data(ABS(matrix%blk_p(blk))),1,pkd_vec(v_offset),1)
             v_offset = v_offset + nze
          ENDIF
       ENDDO ! row
       pkd_vec(matrix%nblkrows_total+1) = REAL(v_offset,dp)
    CASE('row')
       v_offset = matrix%nblkcols_total+2 ! we need nblkcols_total+1 to store rows 
       !                                    and 1 extra for the size of the data
       row = ivec
       row_size = row_blk_size(row)
       DO blk = matrix%row_p(row)+1,matrix%row_p(row+1)
          col = matrix%col_i(blk)
          col_size = col_blk_size(col)
          nze = row_size*col_size
          !
          ! let's copy the block
          pkd_vec(col) = REAL(v_offset,dp)
          CALL dcopy(nze,matrix%data(ABS(matrix%blk_p(blk))),1,pkd_vec(v_offset),1)
          v_offset = v_offset + nze
       ENDDO ! col
       pkd_vec(matrix%nblkcols_total+1) = REAL(v_offset,dp)
    CASE DEFAULT
       CALL stop_program(routineN,"pack what?")       
    END SELECT

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_pack_vec

! *****************************************************************************
!> \brief unpack a row or column of a dbcsr
!> \param[inout]  matrix      DBCSR matrix
!> \param[in]     ivec        the vector to pack
!> \param[in]     pkd_vec     packed vector
!> \param[in]     what        what to unpack (row or column)
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_unpack_vec(matrix, ivec, pkd_vec, what, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, INTENT(IN)                      :: ivec
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_vec
    CHARACTER(LEN=*), INTENT(IN)             :: what
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_unpack_vec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_size, frst_blk, &
                                                last_blk, nze, offset, row, &
                                                row_size, timing_handle, &
                                                v_offset
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    row_blk_size => array_data (matrix%row_blk_size)
    col_blk_size => array_data (matrix%col_blk_size)

    !
    ! let's go
    SELECT CASE(what)
    CASE('column')
       IF(pkd_vec(matrix%nblkrows_total+1).LE.0) CALL stop_program(routineN,"1pkd_vec(n+1).LE.0")
       col = ivec
       col_size = col_blk_size(col)
       DO row = 1, matrix%nblkrows_total
          v_offset = INT(pkd_vec(row))
          IF(v_offset.LE.0) CYCLE
          row_size = row_blk_size(row)
          blk = matrix%row_p(row)+1
          frst_blk = matrix%row_p(row)+1
          last_blk = matrix%row_p(row+1)
          CALL find_column(col,frst_blk,last_blk,matrix%col_i,blk,found)
          IF(found) THEN
             nze = row_size*col_size
             !
             ! let's scale the block
             offset = ABS(matrix%blk_p(blk))
             ! we just dcopy here, could aslo do a daxpy
             !CALL dcopy(nze,pkd_vec(v_offset),1,matrix%data(offset),1)
            matrix%data(offset:offset+nze-1)=pkd_vec(v_offset:v_offset+nze-1)
          ELSE
             !
             ! need to add it
             WRITE(*,*) routineN//' row',row,' col',col
             CALL stop_program(routineN,"some more work here")
          ENDIF
       ENDDO ! row
    CASE('row')
       IF(pkd_vec(matrix%nblkcols_total+1).LE.0) CALL stop_program(routineN,"2pkd_vec(n+1).LE.0")
       row = ivec
       row_size = row_blk_size(row)
       blk = matrix%row_p(row)+1
       frst_blk = matrix%row_p(row)+1
       last_blk = matrix%row_p(row+1)
       DO col = 1,matrix%nblkcols_total
          v_offset = INT(pkd_vec(col))
          IF(v_offset.LE.0) CYCLE
          col_size = col_blk_size(col)
          ! we assume that the columns are ordered !
          CALL find_column(col,frst_blk,last_blk,matrix%col_i,blk,found)
          IF(found) THEN
             nze = row_size*col_size
             offset = ABS(matrix%blk_p(blk))
             ! we just dcopy here, could aslo do a daxpy
             CALL dcopy(nze,pkd_vec(v_offset),1,matrix%data(offset),1)
          ELSE
             !
             ! need to add it
             CALL stop_program(routineN,"some more work here")
          ENDIF
       ENDDO ! row
    CASE DEFAULT
       CALL stop_program(routineN,"unpack what?")
    END SELECT

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_unpack_vec


! *****************************************************************************
!> \brief Does one of the operations:
!>      left)  pkd_ji = alpha * block_jj (block) * pkd_ji for j = 1,..., ivec - 1
!>      right) pkd_ij = alpha * pkd_ij * block_jj (block) for j = 1,..., ivec - 1
!>
! *****************************************************************************
  SUBROUTINE packed_vec_scale(alpha, block, pkd, ivec, vec_blk_size, side, error)
    REAL(dp), INTENT(IN)                     :: alpha
    REAL(dp), DIMENSION(:), INTENT(in)       :: block
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pkd
    INTEGER, INTENT(in)                      :: ivec
    INTEGER, DIMENSION(:), INTENT(in)        :: vec_blk_size
    CHARACTER(LEN=*), INTENT(IN)             :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_scale', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iblk, k, m, n, offset
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff

!   ---------------------------------------------------------------------------

    ALLOCATE(buff(SIZE(pkd)))!too long buff allocated here...
    buff=0.0_dp
    SELECT CASE(side)
    CASE('left')
       iblk = 1
       DO i = 1,ivec-1
          offset = INT(pkd(i))
          IF(offset.GT.0) THEN
             m = vec_blk_size(i)
             n = vec_blk_size(ivec)
             k = vec_blk_size(i)
             CALL dgemm('N','N',m,n,k,alpha,block(iblk),m,pkd(offset),k,0.0_dp,buff(1),m)
             CALL dcopy(m*n,buff(1),1,pkd(offset),1)
          ENDIF
          iblk = iblk + vec_blk_size(i)**2
       ENDDO
    CASE('right')
       iblk = 1
       DO i = 1,ivec-1
          offset = INT(pkd(i))
          IF(offset.GT.0) THEN
             m = vec_blk_size(ivec)
             n = vec_blk_size(i)
             k = vec_blk_size(i)
             CALL dgemm('N','N',m,n,k,alpha,pkd(offset),m,block(iblk),k,0.0_dp,buff(1),m)
             CALL dcopy(m*n,buff(1),1,pkd(offset),1)
          ENDIF
          iblk = iblk + vec_blk_size(i)**2
       ENDDO
    CASE DEFAULT
       CALL stop_program(routineN,"side?")
    END SELECT
    DEALLOCATE(buff)

  END SUBROUTINE packed_vec_scale

! *****************************************************************************
!> \brief Does the operations:
!>      v(k,:) = v(k,:) - uAk(i,k)' * inv(d(i,i)) * v(i,:) / s
!>      u(k,:) = u(k,:) -   v(i,k)' * inv(d(i,i)) * u(i,:) / s
!>
! *****************************************************************************
  SUBROUTINE packed_vec_bif_tech(mat_v, mat_u, pkd_v_fac, pkd_u_fac, ivec, pkd_v, pkd_u, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: mat_v, mat_u
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_v_fac, pkd_u_fac
    INTEGER, INTENT(IN)                      :: ivec
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_v, pkd_u
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_bif_tech', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, col, col_size, k_row, k_row_size, offset, row, row_size, &
      timing_handle, u_offset, u_offset_last, ufac_offset, v_offset, &
      v_offset_last, vfac_offset
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    row_blk_size => array_data (mat_u%row_blk_size)
    col_blk_size => array_data (mat_u%col_blk_size)

    IF(pkd_u(mat_u%nblkrows_total+1).LE.0) CALL stop_program(routineN,"1pkd_u(n+1).LE.0")
    IF(pkd_v(mat_v%nblkrows_total+1).LE.0) CALL stop_program(routineN,"1pkd_v(n+1).LE.0")
    !
    !v(k,:) (pkd) = v(k,:) (pkd) - (inv(d(i,i)) * uAk(i,k))' / s (pkd) * v(i,:) (bcsr) i = 1,ivec-1
    k_row = ivec
    k_row_size = row_blk_size(k_row)
    v_offset_last = INT(pkd_v(mat_v%nblkrows_total+1))
    DO row = 1, mat_v%nblkrows_total
       IF(row.GE.k_row) EXIT !i = 1,ivec-1
       row_size = row_blk_size(row)
       ufac_offset = INT(pkd_u_fac(row))
       IF(ufac_offset.LE.0) CYCLE
       DO blk = mat_v%row_p(row)+1,mat_v%row_p(row+1)
          col = mat_v%col_i(blk)
          col_size = col_blk_size(col)
          offset = ABS(mat_v%blk_p(blk))
          v_offset = INT(pkd_v(col))
          IF(v_offset.LE.0) THEN
             v_offset = v_offset_last
             v_offset_last = v_offset_last + k_row_size * col_size
             pkd_v(v_offset:v_offset+k_row_size * col_size-1) = 0.0_dp
             pkd_v(mat_v%nblkrows_total+1) = REAL(v_offset_last,dp)
          ENDIF
          !
          ! let's multiply and add
          CALL dgemm('T', 'N', k_row_size, col_size, row_size,&
               &    -1.0_dp, pkd_u_fac(ufac_offset), row_size,&
               &             mat_v%data(offset), row_size,&
                  &     1.0_dp, pkd_v(v_offset), k_row_size)
       ENDDO
    ENDDO ! row
    !
    !u(k,:) (pkd) = u(k,:) (pkd) - (inv(d(i,i)) * v(i,k))' / s (pkd) * u(i,:) (bcsr) i = 1,ivec-1
    u_offset_last = INT(pkd_u(mat_u%nblkrows_total+1))
    DO row = 1, mat_u%nblkrows_total
       IF(row.GE.k_row) EXIT !i = 1,ivec-1
       row_size = row_blk_size(row)
       vfac_offset = INT(pkd_v_fac(row))
       IF(vfac_offset.LE.0) CYCLE
       DO blk = mat_u%row_p(row)+1,mat_u%row_p(row+1)
          col = mat_u%col_i(blk)
          col_size = col_blk_size(col)
          offset = ABS(mat_u%blk_p(blk))
          u_offset = INT(pkd_u(col))
          IF(u_offset.LE.0) THEN
             u_offset = u_offset_last
             u_offset_last = u_offset_last + k_row_size * col_size
             pkd_u(u_offset:u_offset+k_row_size * col_size-1) = 0.0_dp
             pkd_u(mat_u%nblkrows_total+1) = REAL(u_offset_last,dp)
          ENDIF
          !
          ! let's multiply and add
          CALL dgemm('T', 'N', k_row_size, col_size, row_size,&
               &    -1.0_dp, pkd_v_fac(vfac_offset), row_size,&
               &             mat_u%data(offset), row_size,&
               &     1.0_dp, pkd_u(u_offset), k_row_size)
       ENDDO
    ENDDO ! row

    IF(pkd_u(mat_u%nblkrows_total+1).LE.0) CALL stop_program(routineN,"2pkd_u(n+1).LE.0")
    IF(pkd_v(mat_v%nblkrows_total+1).LE.0) CALL stop_program(routineN,"2pkd_v(n+1).LE.0")

    CALL timestop(timing_handle)

  END SUBROUTINE packed_vec_bif_tech

  SUBROUTINE packed_vec_bcast(pkd_vec, source, scope, mp_obj, error)
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_vec
    INTEGER, INTENT(IN)                      :: source
    CHARACTER(LEN=*), INTENT(IN)             :: scope
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_bcast', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icol, irow, mp_group, mynode, &
                                                mypcol, myprow, npcols, &
                                                nprows, numnodes, src, &
                                                timing_handle
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    npcols = dbcsr_mp_npcols (mp_obj)
    nprows = dbcsr_mp_nprows (mp_obj)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)

    ALLOCATE(buff(SIZE(pkd_vec)))
    SELECT CASE(scope)
    CASE('rowise')
       !
       ! simple hack
       DO icol = 0,npcols-1
          buff = pkd_vec
          src = blacs2mpi(source,icol)
          CALL mp_bcast(buff,src,mp_group)
          IF(mypcol.EQ.icol) THEN
             pkd_vec = buff
          ENDIF
       ENDDO
    CASE('columnwise')
       !
       ! simple hack
       DO irow = 0,nprows-1
          buff = pkd_vec
          src = blacs2mpi(irow,source)
          CALL mp_bcast(buff,src,mp_group)
          IF(myprow.EQ.irow) THEN
             pkd_vec = buff
          ENDIF
       ENDDO
    CASE('all')
       CALL mp_bcast(pkd_vec,source,mp_group)
    CASE DEFAULT
       CALL stop_program(routineN,"how do you wanna bcast")
    END SELECT
    DEALLOCATE(buff)

    CALL timestop(timing_handle)

  END SUBROUTINE packed_vec_bcast

  SUBROUTINE packed_vec_reduce(pkd_vec, to, scope, vec_blk_size, n, mp_obj, error)
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_vec
    INTEGER, INTENT(IN)                      :: to
    CHARACTER(LEN=*), INTENT(IN)             :: scope
    INTEGER, DIMENSION(:), INTENT(IN)        :: vec_blk_size
    INTEGER, INTENT(IN)                      :: n
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_reduce', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icol, irow, mp_group, mynode, &
                                                mypcol, myprow, npcols, &
                                                nprows, numnodes, src, &
                                                timing_handle
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff, sums

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    npcols = dbcsr_mp_npcols (mp_obj)
    nprows = dbcsr_mp_nprows (mp_obj)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)

    IF(pkd_vec(n+1).LE.0) CALL stop_program(routineN,"a_offset_last.LE.0!!!!!!!")

    ALLOCATE(buff(SIZE(pkd_vec)),sums(SIZE(pkd_vec)))
    sums = 0.0_dp
    SELECT CASE(scope)
    CASE('rowise')
       CALL stop_program(routineN,"more work here")
       !
       ! simple hack
       DO icol = 0,npcols-1
          buff = pkd_vec
          src = blacs2mpi(to,icol)
          CALL mp_bcast(buff,src,mp_group)
          IF(mypcol.EQ.icol) THEN
             pkd_vec = buff
             CALL add_packed_vecs(sums,buff,vec_blk_size,n,error)
          ENDIF
       ENDDO
    CASE('columnwise')
       !CALL stop_program(routineN,"more work here")
       !
       ! simple hack
       DO irow = 0,nprows-1
          buff = pkd_vec
          src = blacs2mpi(irow,to)
          CALL mp_bcast(buff,src,mp_group)
          IF(myprow.EQ.irow) THEN
             CALL add_packed_vecs(sums,buff,vec_blk_size,n,error)
          ENDIF
       ENDDO
    CASE('all')
       CALL stop_program(routineN,"more work here")
    CASE DEFAULT
       CALL stop_program(routineN,"how do you wanna reduce")
    END SELECT
    pkd_vec = sums
    DEALLOCATE(buff,sums)

    IF(pkd_vec(n+1).LE.0) CALL stop_program(routineN,"a_offset_last.LE.0!!!!!!!")

    CALL timestop(timing_handle)

  END SUBROUTINE packed_vec_reduce

  SUBROUTINE add_packed_vecs(pkd_vec_a,pkd_vec_b,vec_blk_size,n,error)
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_vec_a
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_vec_b
    INTEGER, DIMENSION(:), INTENT(IN)        :: vec_blk_size
    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_packed_vecs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a_offset, a_offset_last, &
                                                b_offset, i

!   ---------------------------------------------------------------------------
! this points to the last empty entry

    a_offset_last = INT(pkd_vec_a(n+1)) 
    IF(a_offset_last.LE.0) CALL stop_program(routineN,"a_offset_last.LE.0!!!!!!!")
    DO i = 1,n
       b_offset = INT(pkd_vec_b(i))
       a_offset = INT(pkd_vec_a(i))
       IF(b_offset.LE.0) CYCLE
       IF(a_offset.GT.0) THEN
          ! the block exsits in a, just add
          CALL daxpy(vec_blk_size(i),1.0_dp,pkd_vec_b(b_offset),1,pkd_vec_a(a_offset),1)
       ELSE
          ! the block doesnt exsit in a, copy b at the end
          CALL dcopy(vec_blk_size(i),pkd_vec_b(b_offset),1,pkd_vec_a(a_offset_last),1)
          a_offset_last = a_offset_last + vec_blk_size(i)
       ENDIF
    ENDDO
    ! reset the last empty entry if needed
    pkd_vec_a(n+1) = REAL(a_offset_last,dp)
  END SUBROUTINE add_packed_vecs

  SUBROUTINE bcast_block(blk, source, scope, mp_obj, error)
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: blk
    INTEGER, INTENT(IN)                      :: source
    CHARACTER(LEN=*), INTENT(IN)             :: scope
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'bcast_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icol, irow, mp_group, mynode, &
                                                mypcol, myprow, npcols, &
                                                nprows, numnodes, src
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff

!   ---------------------------------------------------------------------------

    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    npcols = dbcsr_mp_npcols (mp_obj)
    nprows = dbcsr_mp_nprows (mp_obj)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)

    ALLOCATE(buff(SIZE(blk)))
    SELECT CASE(scope)
    CASE('rowise')
       !
       ! simple hack
       DO icol = 0,npcols-1
          buff = blk
          src = blacs2mpi(source,icol)
          CALL mp_bcast(buff,src,mp_group)
          IF(mypcol.EQ.icol) THEN
             blk = buff
          ENDIF
       ENDDO
    CASE('columnwise')
       !
       ! simple hack
       DO irow = 0,nprows-1
          buff = blk
          src = blacs2mpi(irow,source)
          CALL mp_bcast(buff,src,mp_group)
          IF(myprow.EQ.irow) THEN
             blk = buff
          ENDIF
       ENDDO
    CASE('all')
       CALL stop_program(routineN,"more work here")
    CASE DEFAULT
       CALL stop_program(routineN,"how do you wanna bcast")
    END SELECT
    DEALLOCATE(buff)

  END SUBROUTINE bcast_block

  SUBROUTINE block_add_on_diag( m, blk, alpha)
    INTEGER, INTENT(in)                      :: m
    REAL(dp), INTENT(inout)                  :: blk(m,m)
    REAL(dp), INTENT(in)                     :: alpha

    CHARACTER(len=*), PARAMETER :: routineN = 'block_add_on_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    DO i=1,m
       blk(i,i) = blk(i,i) + alpha
    ENDDO
  END SUBROUTINE block_add_on_diag

  SUBROUTINE block_set( m, n, blk, alpha, beta)
    INTEGER                                  :: m, n
    REAL(dp)                                 :: blk(m,n)
    REAL(dp), OPTIONAL                       :: alpha, beta

    CHARACTER(len=*), PARAMETER :: routineN = 'block_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(dp)                                 :: my_alpha, my_beta

!   ---------------------------------------------------------------------------

    my_beta = 0.0_dp
    my_alpha = 0.0_dp
    IF(PRESENT(alpha)) my_alpha = alpha
    IF(PRESENT(beta)) my_beta = beta
    blk(:,:) = my_beta
    IF(m.EQ.n) THEN
       DO i=1,m
          blk(i,i)=my_alpha
       ENDDO
    ENDIF
  END SUBROUTINE block_set

  SUBROUTINE block_chol_inv(m, blk)
    INTEGER                                  :: m
    REAL(dp)                                 :: blk(m,m)

    CHARACTER(len=*), PARAMETER :: routineN = 'block_chol_inv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, info, j

!   ---------------------------------------------------------------------------

    CALL dpotrf( 'U', m, blk, m, info )
    IF(info.NE.0)CALL stop_program(routineN,"error in dpotrf")
    CALL dpotri( 'U', m, blk, m, info )
    IF(info.NE.0)CALL stop_program(routineN,"error in dpotri")
    !
    ! symmetrize
    DO i=1,m
       DO j=i,m
          blk(j,i)=blk(i,j)
       ENDDO
    ENDDO
  END SUBROUTINE block_chol_inv

END MODULE dbcsr_operations

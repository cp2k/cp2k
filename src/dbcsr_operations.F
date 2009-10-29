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
  USE dbcsr_methods
  USE dbcsr_types
  USE dbcsr_util
  USE distribution_2d_types,           ONLY: distribution_2d_retain
  USE kinds,                           ONLY: dp,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8,&
                                             sp
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: &
       mp_allgather, mp_alltoall, mp_bcast, mp_isendrecv, mp_max, mp_rma_get, &
       mp_sum, mp_wait, mp_window_create, mp_window_destroy, mp_window_fence, &
       mp_window_type

  !$ USE OMP_LIB
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_operations'

  PUBLIC :: dbcsr_mult_e_ec_m, dbcsr_trace, dbcsr_add_on_diag,&
            dbcsr_set, dbcsr_scale, dbcsr_add, dbcsr_copy, dbcsr_get_diag, &
            dbcsr_set_diag, dbcsr_get_block, dbcsr_get_block_p,&
            dbcsr_get_block_diag,&
            dbcsr_put_block, dbcsr_filter,&
            dbcsr_replace_blocks,&
            dbcsr_btriu, block_add_on_diag,&
            block_chol_inv, block_set, dbcsr_remove_block, dbcsr_norm
  PUBLIC :: dbcsr_desymmetrize_deep,&
            dbcsr_new_transposed,&
            dbcsr_complete_redistribute,&
            dbcsr_redistribute, dbcsr_make_untransposed_blocks
  PUBLIC :: dbcsr_replicate_all, dbcsr_distribute, dbcsr_sum_replicated
  ! Direct access to the variants because array reshaping and
  ! typematching is incompatible.
  PUBLIC :: block_add_on_diag_d, block_add_on_diag_s,&
            block_add_on_diag_z, block_add_on_diag_c
  PUBLIC :: block_chol_inv_d, block_chol_inv_s,&
            block_chol_inv_z, block_chol_inv_c
  PUBLIC :: block_set_d, block_set_s,&
            block_set_z, block_set_c

  PUBLIC :: dbcsr_reserve_block2d, dbcsr_get_tmp_block2d

  PRIVATE


! The interfaces for the generic routines found in the generated
! generic files.

  INTERFACE dbcsr_trace
     MODULE PROCEDURE dbcsr_trace_a_s, dbcsr_trace_a_d,&
          dbcsr_trace_a_c, dbcsr_trace_a_z,&
          dbcsr_trace_ab_s, dbcsr_trace_ab_d,&
          dbcsr_trace_ab_c, dbcsr_trace_ab_z
  END INTERFACE

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
  END INTERFACE

  INTERFACE dbcsr_put_block
     MODULE PROCEDURE dbcsr_put_block_d, dbcsr_put_block_s,&
                      dbcsr_put_block_z, dbcsr_put_block_c,&
                      dbcsr_put_block2d_d, dbcsr_put_block2d_s,&
                      dbcsr_put_block2d_z, dbcsr_put_block2d_c
  END INTERFACE

  INTERFACE block_add_on_diag
     MODULE PROCEDURE block_add_on_diag_d, block_add_on_diag_s,&
          block_add_on_diag_z, block_add_on_diag_c
  END INTERFACE

  INTERFACE block_chol_inv
     MODULE PROCEDURE block_chol_inv_d, block_chol_inv_s,&
          block_chol_inv_z, block_chol_inv_c
  END INTERFACE

  INTERFACE block_set
     MODULE PROCEDURE block_set_d, block_set_s,&
          block_set_z, block_set_c
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


CONTAINS


! *****************************************************************************
!> \brief Multiplies two DBCSR matrices using non-RMA transfers.
!>
!> All matrices have the same processor grid distribution.
!>
!> \param[in] left_set             set of left column images
!> \param[in] right_matrix         right matrix
!> \param[out] product             DBCSR product matrix
!> \param[in] submatrix_selection  select submatrix
!> \param[in,out] error            cp2k error
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix; default is no
!> \param[in] alpha_d              (optional) product multiplication factor
!> \param[in] beta_d               (optional) existing target multiplication
!>                                 factor
!> \param[out] flop                (optional) effective flop
! *****************************************************************************
  SUBROUTINE dbcsr_mult_e_ec_m (left_set, right_matrix, product_matrix,&
       submatrix_selection, error, retain_sparsity, alpha_d, beta_d, flop)
    TYPE(dbcsr_array_type), INTENT(INOUT)    :: left_set
    TYPE(dbcsr_obj), INTENT(IN)              :: right_matrix
    TYPE(dbcsr_obj), INTENT(INOUT), TARGET   :: product_matrix
    INTEGER, DIMENSION(6), INTENT(in)        :: submatrix_selection
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: alpha_d, beta_d
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mult_e_ec_m', &
      routineP = moduleN//':'//routineN

    INTEGER :: debug_level, i, img, ithread, left_dst_image, left_dst_p, &
      left_dst_pcol, left_dst_prow, left_multiplicity, left_mypcol, &
      left_myprow, left_nimages, left_npcols, left_nprows, left_p_recv, &
      left_p_send, left_src_image, left_src_p, left_src_pcol, left_src_prow, &
      mp_group, mult, mynode, mypcol, myprow, myt_mult, nblks, npcols, &
      nprows, nthreads, numnodes, nze, output_unit, right_data_rr, &
      right_data_sr, right_dst_p, right_dst_prow, right_index_rr, &
      right_index_sr, right_mypcol, right_myprow, right_npcols, right_nprows, &
      right_p_recv, right_p_send, right_prow_i, right_src_p, &
      right_src_p_calc, right_src_prow, timing_handle
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
    TYPE(dbcsr_data_obj)                     :: right_matrix_data
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_obj), POINTER                 :: product_matrix_p, &
                                                rbuffer_calc, rbuffer_comm
    TYPE(dbcsr_obj), TARGET                  :: rbuffer_1, rbuffer_2

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    debug_level=0
    t_all = 0.0_dp
    t_dgemm = 0.0_dp
    IF (dbg) debug_level = 2
    IF (dbcsr_mp_nprows (dbcsr_distribution_mp (right_matrix%m%dist)) .NE.&
         dbcsr_mp_npcols (dbcsr_distribution_mp (right_matrix%m%dist))) THEN
       CALL dbcsr_mult_e_ec_m_rma (left_set, right_matrix,&
            product_matrix, submatrix_selection, error=error,&
            alpha_d=alpha_d, beta_d=beta_d, retain_sparsity=retain_sparsity,&
            flop=flop_total)
       IF (PRESENT (flop)) flop = flop_total
       RETURN
    ENDIF
    CALL timeset (routineN, timing_handle)
    product_matrix_p => product_matrix
    !!$OMP PARALLEL DEFAULT (SHARED) &
    !!$OMP          PRIVATE (myt_mult, right_prow_i, ithread)
    nthreads = 1
    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM ()
    !$OMP MASTER
    !$ nthreads = OMP_GET_NUM_THREADS ()
    IF (.NOT. ASSOCIATED (product_matrix%m%wms)) THEN
       CALL dbcsr_work_create(product_matrix, n=nthreads, error=error)
    ENDIF
    !$OMP FLUSH (product_matrix)
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
    mp_obj = dbcsr_distribution_mp (dbcsr_distribution (product_matrix))
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
    DO i = 1, nthreads
       CALL dbcsr_data_ensure_size(product_matrix%m%wms(i)%data_area,&
            MAX(product_matrix%m%wms(i)%datasize,&
            dbcsr_nfullrows_local (product_matrix)&
            *dbcsr_nfullcols_local (product_matrix)))
            !dbcsr_get_nfullrowsleft_set%m(1)%nfullrows_local*right_matrix%nfullcols_local))
       CALL ensure_array_size(product_matrix%m%wms(i)%row_i,&
            ub=product_matrix%m%nblkcols_total*product_matrix%m%nblkrows_total, error=error)
       CALL ensure_array_size(product_matrix%m%wms(i)%col_i,&
            ub=product_matrix%m%nblkcols_total*product_matrix%m%nblkrows_total, error=error)
       CALL ensure_array_size(product_matrix%m%wms(i)%blk_p,&
            ub=product_matrix%m%nblkcols_total*product_matrix%m%nblkrows_total, error=error)

       !left_data_rp => dbcsr_get_data_p (product_matrix%m%wms(i)%data_area, 0.0_dp)
       !left_data_rp(:) = 0.0_dp

    ENDDO
    DBG 'pgrid shape', SHAPE (pgrid)
    IF (dbg .AND. mynode .EQ. 0) THEN
       WRITE(*,*) 'left images, multiplicity', left_nimages, left_multiplicity
    ENDIF
    ALLOCATE (my_sizes(4, left_nimages))
    my_sizes(:,:) = 0
    DO img = 1, left_nimages
       my_sizes(1, img) = dbcsr_get_data_size (left_set%mats(img)%m%data_area)
       my_sizes(2, img) = SIZE (left_set%mats(img)%m%index)
    ENDDO
    my_sizes(3, 1) = dbcsr_get_data_size (right_matrix%m%data_area)
    my_sizes(4, 1) = SIZE (right_matrix%m%index)
    ALLOCATE (all_sizes(4, LBOUND(my_sizes,2):UBOUND(my_sizes,2), 0:numnodes-1))
    DO i = 1, 4
       DBG 'my sizes', my_sizes(i,:)
    ENDDO
    CALL mp_allgather(my_sizes, all_sizes, mp_group)
    nblks = MAXVAL (all_sizes (2, :, :))
    nze = MAXVAL (all_sizes(1, :, :))
    DBG 'data sizes', all_sizes
    ! Setup the buffer matrices
    ALLOCATE (lbuffer_1%mats(left_multiplicity))
    ALLOCATE (lbuffer_2%mats(left_multiplicity))
    ALLOCATE (left_data_sr(left_multiplicity), left_data_rr(left_multiplicity))
    ALLOCATE (left_index_sr(left_multiplicity), left_index_rr(left_multiplicity))
    DO mult = 1, left_multiplicity
       CALL dbcsr_init (lbuffer_1%mats(mult), error)
       CALL dbcsr_create (lbuffer_1%mats(mult),&
            'lbufferl1'//left_set%mats(1)%m%name, left_set%mats(1)%m%dist, 'N',&
            left_set%mats(1)%m%row_blk_size, left_set%mats(1)%m%col_blk_size,&
            nblks, nze, special=.FALSE., error=error)
            lbuffer_1%mats(mult)%m%sym_negation = left_set%mats(1)%m%sym_negation
       CALL dbcsr_data_ensure_size (lbuffer_1%mats(mult)%m%data_area,&
            MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.)
       !CALL ensure_array_size (lbuffer_1%mats(mult)%m%data,&
       !     ub=MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.,&
       !     special=.FALSE., error=error)
       CALL ensure_array_size (lbuffer_1%mats(mult)%m%index,&
            ub=MAXVAL (all_sizes(2, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
       !
       CALL dbcsr_init (lbuffer_2%mats(mult), error)
       CALL dbcsr_create (lbuffer_2%mats(mult),&
            'lbuffer2'//left_set%mats(1)%m%name, left_set%mats(1)%m%dist, 'N',&
            left_set%mats(1)%m%row_blk_size, left_set%mats(1)%m%col_blk_size,&
            nblks, nze, special=.FALSE., error=error)
            lbuffer_2%mats(mult)%m%sym_negation = left_set%mats(1)%m%sym_negation
       CALL dbcsr_data_ensure_size (lbuffer_2%mats(mult)%m%data_area,&
            MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.)
       !CALL ensure_array_size (lbuffer_2%mats(mult)%m%data,&
       !     ub=MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.,&
       !     special=.FALSE., error=error)
       CALL ensure_array_size (lbuffer_2%mats(mult)%m%index,&
            ub=MAXVAL (all_sizes(2, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
    ENDDO
    lbuffer_calc => lbuffer_1
    lbuffer_comm => lbuffer_2
    ! Setup buffers for right matrix
    nblks = MAXVAL (all_sizes (4, 1, :))
    nze = MAXVAL (all_sizes(3, 1, :))
    right_matrix_data = right_matrix%m%data_area
    CALL dbcsr_data_ensure_size (right_matrix_data,&
         MAX(dbcsr_get_data_size(right_matrix_data), 0))
    CALL dbcsr_init (rbuffer_1, error)
    CALL dbcsr_create (rbuffer_1,&
         'rbuffer1'//right_matrix%m%name, right_matrix%m%dist,&
         'N',&
         right_matrix%m%row_blk_size,&
         right_matrix%m%col_blk_size,&
         0, 0, special=.FALSE., error=error)
         rbuffer_1%m%sym_negation = right_matrix%m%sym_negation
    CALL dbcsr_data_ensure_size (rbuffer_1%m%data_area,&
         MAXVAL (all_sizes(3, 1, :)), nocopy=.TRUE.)
    !CALL ensure_array_size (rbuffer_1%data,&
    !     ub=MAXVAL (all_sizes(3, 1, :)), nocopy=.TRUE., special=.FALSE.,&
    !     error=error)
    CALL ensure_array_size (rbuffer_1%m%index,&
         ub=MAXVAL (all_sizes(4, 1, :)), nocopy=.TRUE., special=.FALSE.,&
         error=error)
    !
    CALL dbcsr_init (rbuffer_2, error)
    CALL dbcsr_create (rbuffer_2,&
         'rbuffer2'//right_matrix%m%name, right_matrix%m%dist,&
         'N',&
         right_matrix%m%row_blk_size,&
         right_matrix%m%col_blk_size,&
         0, 0, special=.FALSE., error=error)
         rbuffer_2%m%sym_negation = right_matrix%m%sym_negation
    CALL dbcsr_data_ensure_size (rbuffer_2%m%data_area,&
         MAXVAL (all_sizes(3, 1, :)), nocopy=.TRUE.)
    !CALL ensure_array_size (rbuffer_2%data,&
    !     ub=MAXVAL (all_sizes(3, 1, :)), nocopy=.TRUE., special=.FALSE.,&
    !     error=error)
    CALL ensure_array_size (rbuffer_2%m%index,&
         ub=MAXVAL (all_sizes(4, 1, :)), nocopy=.TRUE., special=.FALSE.,&
         error=error)
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
    !$OMP FLUSH (right_nprows)
    !$OMP END MASTER
    !$OMP BARRIER
    !WRITE(*,*)ithread,'going zwischen',0,right_nprows
    k_index: DO right_prow_i = 0, right_nprows
       !$OMP MASTER
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
          right_data_rp => pointer_view (&
               dbcsr_get_data_p (rbuffer_comm%m%data_area, 0.0_dp),&
               1, all_sizes(3,1,right_src_p))
          right_index_rp => rbuffer_comm%m%index(1:all_sizes(4,1,right_src_p))
          right_data_sp => pointer_view (&
               dbcsr_get_data_p (rbuffer_calc%m%data_area, 0.0_dp),&
               1, all_sizes(3,1,right_dst_p))
          right_index_sp => rbuffer_calc%m%index(1:all_sizes(4,1,right_dst_p))
          IF (right_prow_i .EQ. 0) THEN
             CALL mp_isendrecv (dbcsr_get_data_p (right_matrix%m%data_area, 0.0_dp),&
                  right_dst_p,&
                  right_data_rp, right_src_p, &
                  mp_group, right_data_sr, right_data_rr, tag=3)
             CALL mp_isendrecv (right_matrix%m%index, right_dst_p, &
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
          CALL dbcsr_repoint_index (rbuffer_calc%m)
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
             left_data_rp => pointer_view (&
                     dbcsr_get_data_p (lbuffer_comm%mats(mult)%m%data_area, 0.0_dp),&
                     1, all_sizes(1, left_src_image, left_src_p))
             left_index_rp => lbuffer_comm%mats(mult)%m%index(&
                  1:all_sizes(2, left_src_image, left_src_p))
             left_data_sp => pointer_view (&
                     dbcsr_get_data_p (lbuffer_calc%mats(mult)%m%data_area, 0.0_dp),&
                     1, all_sizes(1, left_dst_image, left_dst_p))
             left_index_sp => lbuffer_calc%mats(mult)%m%index(&
                  1:all_sizes(2, left_dst_image, left_dst_p))
             IF (right_prow_i .EQ. 0) THEN
                CALL mp_isendrecv (&
                     dbcsr_get_data_p (left_set%mats(mult)%m%data_area, 0.0_dp),&
                     left_dst_p,&
                     left_data_rp, left_src_p, mp_group,&
                     left_data_sr(mult), left_data_rr (mult), tag=1)
                CALL mp_isendrecv (left_set%mats(mult)%m%index, left_dst_p,&
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
       !$OMP END MASTER
       !$OMP BARRIER
       calc_case_left: IF (right_prow_i .GT. 0) THEN
             !left_src_p_calc = left_pgrid (left_prow, right_prow-1)
          DO myt_mult = 1, left_multiplicity
             !$OMP MASTER
             mult = myt_mult
             CALL dbcsr_repoint_index (lbuffer_calc%mats(myt_mult)%m)
             tstart = m_walltime ()
             !$OMP END MASTER
             !$OMP BARRIER
             !WRITE(*,*)ithread,'running local multiply', right_prow_i, myt_mult, nthreads
             CALL dbcsr_nn_mult_lin(lbuffer_calc%mats(myt_mult)%m, rbuffer_calc%m,&
                  product_matrix_p%m, submatrix_box=submatrix_selection,&
                  retain_sparsity=retain_sparsity,&
                  alpha_d=alpha_d,&
                  flop=flop_single,&
                  t_all=t_all, t_dgemm=t_dgemm)
             !$OMP BARRIER
             !$OMP MASTER
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
             !$OMP END MASTER
          ENDDO
       ENDIF calc_case_left
       !$OMP BARRIER
       !$OMP MASTER
       CALL dbcsr_switch_sets (lbuffer_calc, lbuffer_comm)
       CALL dbcsr_switch_m_ptrs (rbuffer_calc%m, rbuffer_comm%m)
       !$OMP END MASTER
       !$BARRIER
    ENDDO k_index
    !$BARRIER
    !$OMP MASTER
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
    !
    CALL dbcsr_destroy (rbuffer_1, force=.TRUE.)
    CALL dbcsr_destroy (rbuffer_2, force=.TRUE.)
    CALL dbcsr_destroy_array (lbuffer_1)
    CALL dbcsr_destroy_array (lbuffer_2)
    DEALLOCATE (my_sizes)
    !
    !$OMP END MASTER
    !!$OMP END PARALLEL
    CALL dbcsr_finalize(product_matrix)
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
!> \param[in] submatrix_selection  select submatrix
!> \param[in,out] error            cp2k error
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix, default is no
!> \param[in] alpha_d              (optional) product multiplication factor
!> \param[in] beta_d               (optional) existing target multiplication
!>                                 factor
!> \param[out] flop                (optional) effective flop
! *****************************************************************************
  SUBROUTINE dbcsr_mult_e_ec_m_rma (left_set, right_matrix, product_matrix,&
       submatrix_selection, error, retain_sparsity, alpha_d, beta_d, flop)
    TYPE(dbcsr_array_type), INTENT(INOUT)    :: left_set
    TYPE(dbcsr_obj), INTENT(IN)              :: right_matrix
    TYPE(dbcsr_obj), INTENT(INOUT), TARGET   :: product_matrix
    INTEGER, DIMENSION(6), INTENT(in)        :: submatrix_selection
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: alpha_d, beta_d
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mult_e_ec_m_rma', &
      routineP = moduleN//':'//routineN

    INTEGER :: debug_level, i, img, ithread, left_multiplicity, left_mypcol, &
      left_myprow, left_nimages, left_npcols, left_nprows, left_src_image, &
      left_src_p_calc, left_src_p_comm, left_src_pcol, left_src_prow, &
      left_virt_src_pcol, mp_group, mult, mynode, mypcol, myprow, nblks, &
      npcols, nprows, nthreads, numnodes, nze, output_unit, right_mypcol, &
      right_myprow, right_npcols, right_nprows, right_prow_i, &
      right_src_p_calc, right_src_p_comm, right_src_prow, timing_handle
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
    TYPE(dbcsr_data_obj)                     :: data_modifier
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_obj), POINTER                 :: product_matrix_p, &
                                                rbuffer_calc, rbuffer_comm
    TYPE(dbcsr_obj), TARGET                  :: rbuffer_1, rbuffer_2
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
    product_matrix_p => product_matrix
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
    mp_obj = dbcsr_distribution_mp (product_matrix%m%dist)
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
    nthreads = 1
    ithread = 0
    !!$ ithread = OMP_GET_THREAD_NUM ()
    !!$ nthreads = OMP_GET_NUM_THREADS ()
    IF (.NOT. ASSOCIATED (product_matrix%m%wms)) THEN
       CALL dbcsr_work_create(product_matrix, n=nthreads, error=error)
    ENDIF
    DO i = 1, nthreads
       CALL dbcsr_data_ensure_size(product_matrix%m%wms(i)%data_area,&
            MAX(product_matrix%m%wms(i)%datasize,&
            dbcsr_nfullrows_local (product_matrix)&
            *dbcsr_nfullcols_local (product_matrix)))
            !dbcsr_get_nfullrowsleft_set%m(1)%nfullrows_local*right_matrix%m%nfullcols_local))
       CALL ensure_array_size(product_matrix%m%wms(i)%row_i,&
            ub=product_matrix%m%nblkcols_total*product_matrix%m%nblkrows_total, error=error)
       CALL ensure_array_size(product_matrix%m%wms(i)%col_i,&
            ub=product_matrix%m%nblkcols_total*product_matrix%m%nblkrows_total, error=error)
       CALL ensure_array_size(product_matrix%m%wms(i)%blk_p,&
            ub=product_matrix%m%nblkcols_total*product_matrix%m%nblkrows_total, error=error)

       !left_data_rp => dbcsr_get_data_p (product_matrix%m%wms(i)%data_area, 0.0_dp)
       !left_data_rp(:) = 0.0_dp

    ENDDO
    DBG 'pgrid shape', SHAPE (pgrid)
    DBG 'left images, multiplicity', left_nimages, left_multiplicity
    ALLOCATE (left_set%data_windows(left_nimages))
    ALLOCATE (left_set%index_windows(left_nimages))
    ALLOCATE (my_sizes(4, left_nimages))
    my_sizes(:,:) = 0
    DO img = 1, left_nimages
       my_sizes(1, img) = dbcsr_get_data_size (left_set%mats(img)%m%data_area)
       my_sizes(2, img) = SIZE (left_set%mats(img)%m%index)
       data_modifier = left_set%mats(img)%m%data_area
       CALL dbcsr_data_ensure_size (data_modifier, 0)
    ENDDO
    my_sizes(3, 1) = dbcsr_get_data_size (right_matrix%m%data_area)
    my_sizes(4, 1) = SIZE (right_matrix%m%index)
    data_modifier = right_matrix%m%data_area
    CALL dbcsr_data_ensure_size (data_modifier, 0)
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
            dbcsr_get_data_p (left_set%mats(img)%m%data_area, 0.0_dp), gid=mp_group)
       CALL mp_window_create (left_set%index_windows(img),&
            left_set%mats(img)%m%index, gid=mp_group)
    ENDDO
    ALLOCATE (lbuffer_1%mats(left_multiplicity))
    ALLOCATE (lbuffer_2%mats(left_multiplicity))
    DO mult = 1, left_multiplicity
       CALL dbcsr_init (lbuffer_1%mats(mult), error)
       CALL dbcsr_create (lbuffer_1%mats(mult),&
            'lbufferl1'//left_set%mats(1)%m%name, left_set%mats(1)%m%dist, 'N',&
            left_set%mats(1)%m%row_blk_size, left_set%mats(1)%m%col_blk_size,&
            nblks, nze, special=.FALSE., error=error)
            lbuffer_1%mats(mult)%m%sym_negation = left_set%mats(1)%m%sym_negation
       CALL dbcsr_data_ensure_size (lbuffer_1%mats(mult)%m%data_area,&
            MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.)
       !CALL ensure_array_size (lbuffer_1%mats(mult)%m%data,&
       !     ub=MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.,&
       !     special=.FALSE., error=error)
       CALL ensure_array_size (lbuffer_1%mats(mult)%m%index,&
            ub=MAXVAL (all_sizes(2, :, :)), nocopy=.TRUE.,&
            special=.FALSE., error=error)
       !
       CALL dbcsr_init (lbuffer_2%mats(mult), error)
       CALL dbcsr_create (lbuffer_2%mats(mult),&
            'lbuffer2'//left_set%mats(1)%m%name, left_set%mats(1)%m%dist, 'N',&
            left_set%mats(1)%m%row_blk_size, left_set%mats(1)%m%col_blk_size,&
            nblks, nze, special=.FALSE., error=error)
            lbuffer_2%mats(mult)%m%sym_negation = left_set%mats(1)%m%sym_negation
       CALL dbcsr_data_ensure_size (lbuffer_2%mats(mult)%m%data_area,&
            MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.)
       !CALL ensure_array_size (lbuffer_2%mats(mult)%m%data,&
       !     ub=MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.,&
       !     special=.FALSE., error=error)
       CALL ensure_array_size (lbuffer_2%mats(mult)%m%index,&
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
         'rbuffer1'//right_matrix%m%name, right_matrix%m%dist,&
         'N',&
         right_matrix%m%row_blk_size,&
         right_matrix%m%col_blk_size,&
         0, 0, special=.FALSE., error=error)
         rbuffer_1%m%sym_negation = right_matrix%m%sym_negation
    CALL dbcsr_data_ensure_size (rbuffer_1%m%data_area,&
         MAXVAL (all_sizes(3, 1, :)), nocopy=.TRUE.)
    !CALL ensure_array_size (rbuffer_1%data,&
    !     ub=MAXVAL (all_sizes(3, 1, :)), nocopy=.TRUE., special=.FALSE.,&
    !     error=error)
    CALL ensure_array_size (rbuffer_1%m%index,&
         ub=MAXVAL (all_sizes(4, 1, :)), nocopy=.TRUE., special=.FALSE.,&
         error=error)
    !
    CALL dbcsr_init (rbuffer_2, error)
    CALL dbcsr_create (rbuffer_2,&
         'rbuffer2'//right_matrix%m%name, right_matrix%m%dist,&
         'N',&
         right_matrix%m%row_blk_size,&
         right_matrix%m%col_blk_size,&
         0, 0, special=.FALSE., error=error)
         rbuffer_2%m%sym_negation = right_matrix%m%sym_negation
    CALL dbcsr_data_ensure_size (rbuffer_2%m%data_area,&
         MAXVAL (all_sizes(3, 1, :)), nocopy=.TRUE.)
    !CALL ensure_array_size (rbuffer_2%data,&
    !     ub=MAXVAL (all_sizes(3, 1, :)), nocopy=.TRUE., special=.FALSE.,&
    !     error=error)
    CALL ensure_array_size (rbuffer_2%m%index,&
         ub=MAXVAL (all_sizes(4, 1, :)), nocopy=.TRUE., special=.FALSE.,&
         error=error)
    ! Create windows for the right matrices
    CALL mp_window_create (right_data_window,&
         dbcsr_get_data_p (right_matrix%m%data_area, 0.0_dp), gid=mp_group)
    CALL mp_window_create (right_index_window, right_matrix%m%index, gid=mp_group)
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
               dbcsr_get_data_p (rbuffer_comm%m%data_area, 0.0_dp))
          CALL mp_rma_get (right_index_window,&
               right_src_p_comm, 0, all_sizes(4, 1, right_src_p_comm),&
               rbuffer_comm%m%index)
       ENDIF xfer_case_right
       ! No calculations in the first loop
       calc_case_right: IF (right_prow_i .GT. 0) THEN
          IF (all_sizes(4, 1, right_src_p_calc) .EQ. 0) THEN
             DBG 'ZERO SIZED INDEX'
          ENDIF
          CALL dbcsr_repoint_index (rbuffer_calc%m)
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
       !left_src_image = 1 + MOD (right_src_prow, left_nimages)
       xfer_case_left: IF (right_prow_i .LT. nprows) THEN
          DO mult = 1, left_multiplicity
             left_virt_src_pcol = right_src_prow*left_multiplicity+mult-1
             left_src_pcol = INT (left_virt_src_pcol / left_nimages)
             left_src_image = 1 + MOD (left_virt_src_pcol, left_nimages)
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
                  dbcsr_get_data_p (lbuffer_comm%mats(mult)%m%data_area, 0.0_dp))
             CALL mp_rma_get (left_set%index_windows(left_src_image),&
                  left_src_p_comm, 0,&
                  all_sizes(2, left_src_image, left_src_p_comm),&
                  lbuffer_comm%mats(mult)%m%index)
          ENDDO
       ENDIF xfer_case_left
       calc_case_left: IF (right_prow_i .GT. 0) THEN
             !left_src_p_calc == left_pgrid (left_prow, right_prow-1)
          DO mult = 1, left_multiplicity
             CALL dbcsr_repoint_index (lbuffer_calc%mats(mult)%m)
             tstart = m_walltime ()
             CALL dbcsr_nn_mult_lin(lbuffer_calc%mats(mult)%m, rbuffer_calc%m,&
                  product_matrix_p%m, submatrix_box=submatrix_selection,&
                  retain_sparsity=retain_sparsity,&
                  alpha_d=alpha_d,&
                  flop=flop_single,&
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
       CALL dbcsr_switch_m_ptrs (rbuffer_calc%m, rbuffer_comm%m)
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
    CALL dbcsr_destroy (rbuffer_1, force=.TRUE.)
    CALL dbcsr_destroy (rbuffer_2, force=.TRUE.)
    CALL dbcsr_destroy_array (lbuffer_1)
    CALL dbcsr_destroy_array (lbuffer_2)
    DEALLOCATE (my_sizes)
    !
    CALL dbcsr_finalize(product_matrix, .TRUE., .FALSE.)
    !
    IF (PRESENT (flop)) flop = flop_total
    CALL timestop (timing_handle)
  END SUBROUTINE dbcsr_mult_e_ec_m_rma


!! *****************************************************************************
!!> \brief Multiplies two DBCSR matrices.
!!> \param[in] left_set             set of left images
!!> \param[in] right_matrix         right matrix
!!> \param[out] product             DBCSR product workmatrix
!!> \param[in,out] error            cp2k error
!!> \param[out] flop                (optional) effective flop
!! *****************************************************************************
!  SUBROUTINE dbcsr_multiply_from_left (left_set, right_matrix, product_matrix,&
!       error, flop)
!    TYPE(dbcsr_array_type), INTENT(INOUT)    :: left_set
!    TYPE(dbcsr_type), INTENT(IN)             :: right_matrix
!    TYPE(dbcsr_type), INTENT(INOUT)          :: product_matrix
!    TYPE(cp_error_type), INTENT(INOUT)       :: error
!    INTEGER(KIND=int_8), INTENT(OUT), &
!      OPTIONAL                               :: flop
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_from_left', &
!      routineP = moduleN//':'//routineN
!
!    INTEGER :: debug_level, entries, i, img, left_image, left_images, &
!      left_multiplicity, left_nprows, left_prow, left_src_p_calc, &
!      left_src_p_comm, mp_group, mult, mynode, mypcol, myprow, nblks, npcols, &
!      nprows, numnodes, nze, output_unit, right_nprows, right_prow, &
!      right_prow_i, right_src_p_calc, right_src_p_comm, timing_handle
!    INTEGER(KIND=int_8)                      :: flop_single, flop_total
!    INTEGER(KIND=int_8), ALLOCATABLE, &
!      DIMENSION(:)                           :: keys
!    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: my_sizes
!    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: all_sizes
!    INTEGER, DIMENSION(:, :), POINTER        :: left_pgrid, product_pgrid, &
!                                                right_pgrid
!    REAL(KIND=dp)                            :: trun, trun_t, tstart, tstop
!    REAL(KIND=dp), ALLOCATABLE, &
!      DIMENSION(:, :)                        :: my_waits
!    REAL(KIND=dp), ALLOCATABLE, &
!      DIMENSION(:, :, :)                     :: all_waits
!    TYPE(cp_logger_type), POINTER            :: logger
!    TYPE(dbcsr_array_type), POINTER          :: lbuffer_calc, lbuffer_comm
!    TYPE(dbcsr_array_type), TARGET           :: lbuffer_1, lbuffer_2
!    TYPE(dbcsr_mp_obj)                       :: mp_obj
!    TYPE(dbcsr_type), POINTER                :: rbuffer_calc, rbuffer_comm
!    TYPE(dbcsr_type), TARGET                 :: rbuffer_1, rbuffer_2
!    TYPE(mp_window_type)                     :: right_data_window, &
!                                                right_index_window
!
!    DEBUG_HEADER
!!   ---------------------------------------------------------------------------
!    DEBUG_BODY
!    debug_level=0
!    IF (dbg) debug_level = 2
!    IF (dbcsr_mp_pgrid_equal (dbcsr_distribution_mp (left_set%image_dist%main),&
!         dbcsr_distribution_mp (right_matrix%dist))) THEN
!       CALL dbcsr_multiply_from_left_sqpg (left_set, right_matrix,&
!            product_matrix, error, flop=flop_total)
!       IF (PRESENT (flop)) flop = flop_total
!       RETURN
!    ENDIF
!    CALL timeset (routineN, timing_handle)
!    CALL btree_new (product_matrix%w%tree)
!    tstart = 0.0_dp ; tstop = 0.0_dp ; trun = 0.0_dp
!    flop_total = 0
!    trun_t = m_walltime ()
!    left_images = SIZE (left_set%mats)
!    CALL bcsr_assert (debug_level,&
!         left_images .EQ. left_set%image_dist%row_decimation, routineP,&
!         "Non-compatible images")
!    left_multiplicity = left_set%image_dist%row_multiplicity
!    mp_obj = dbcsr_distribution_mp (product_matrix%dist)
!    numnodes = dbcsr_mp_numnodes (mp_obj)
!    mynode = dbcsr_mp_mynode (mp_obj)
!    myprow = dbcsr_mp_myprow (mp_obj)
!    mypcol = dbcsr_mp_mypcol (mp_obj)
!    npcols = dbcsr_mp_npcols (mp_obj)
!    nprows = dbcsr_mp_nprows (mp_obj)
!    mp_group = dbcsr_mp_group (mp_obj)
!    product_pgrid => dbcsr_mp_pgrid (mp_obj)
!    left_pgrid => dbcsr_mp_pgrid (dbcsr_distribution_mp (&
!         left_set%image_dist%main))
!    left_nprows = dbcsr_mp_nprows (dbcsr_distribution_mp (&
!         left_set%image_dist%main))
!    right_pgrid => dbcsr_mp_pgrid (dbcsr_distribution_mp (&
!         right_matrix%dist))
!    DBG 'product pgrid shape', SHAPE (dbcsr_mp_pgrid (mp_obj))
!    DBG 'left pgrid shape', SHAPE (left_pgrid),';',left_nprows,left_images
!    DBG 'right pgrid shape', SHAPE (right_pgrid)
!    DBG 'left images, multiplicity', left_images, left_multiplicity
!    CALL bcsr_assert (debug_level,&
!         left_multiplicity * npcols .EQ. left_nprows, routineP,&
!         "Mulitplicity mismatch and processer grid.")
!    ALLOCATE (left_set%data_windows(left_images))
!    ALLOCATE (left_set%index_windows(left_images))
!    ALLOCATE (my_sizes(4, left_images))
!    my_sizes(:,:) = 0
!    DO img = 1, left_images
!       my_sizes(1, img) = SIZE (left_set%mats(img)%data)
!       my_sizes(2, img) = SIZE (left_set%mats(img)%index)
!    ENDDO
!    my_sizes(3, 1) = SIZE (right_matrix%data)
!    my_sizes(4, 1) = SIZE (right_matrix%index)
!    ALLOCATE (all_sizes(4, LBOUND(my_sizes,2):UBOUND(my_sizes,2), 0:numnodes-1))
!    DO i = 1, 4
!       DBG 'my sizes', my_sizes(i,:)
!    ENDDO
!    CALL mp_allgather(my_sizes, all_sizes, mp_group)
!    nblks = MAXVAL (all_sizes (2, :, :))
!    nze = MAXVAL (all_sizes(1, :, :))
!    DBG 'data sizes', all_sizes
!    ! Setup the buffer matrices
!    DO img = 1, left_images
!       CALL mp_window_create (left_set%data_windows(img),&
!            left_set%mats(img)%data, gid=mp_group)
!       CALL mp_window_create (left_set%index_windows(img),&
!            left_set%mats(img)%index, gid=mp_group)
!    ENDDO
!    ALLOCATE (lbuffer_1%mats(left_multiplicity))
!    ALLOCATE (lbuffer_2%mats(left_multiplicity))
!    DO mult = 1, left_multiplicity
!       CALL dbcsr_init (lbuffer_1%mats(mult), error)
!       CALL dbcsr_create (lbuffer_1%mats(mult),&
!            'lbufferl1'//left_set%mats(1)%m%name, left_set%mats(1)%m%dist, 'N',&
!            left_set%mats(1)%m%row_blk_size, left_set%mats(1)%m%col_blk_size,&
!            nblks, nze, special=.FALSE., error=error)
!       CALL ensure_array_size (lbuffer_1%mats(mult)%data,&
!            ub=MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.,&
!            special=.FALSE., error=error)
!       CALL ensure_array_size (lbuffer_1%mats(mult)%index,&
!            ub=MAXVAL (all_sizes(2, :, :)), nocopy=.TRUE.,&
!            special=.FALSE., error=error)
!       !
!       CALL dbcsr_init (lbuffer_2%mats(mult), error)
!       CALL dbcsr_create (lbuffer_2%mats(mult),&
!            'lbuffer2'//left_set%mats(1)%m%name, left_set%mats(1)%m%dist, 'N',&
!            left_set%mats(1)%m%row_blk_size, left_set%mats(1)%m%col_blk_size,&
!            nblks, nze, special=.FALSE., error=error)
!       CALL ensure_array_size (lbuffer_2%mats(mult)%data,&
!            ub=MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.,&
!            special=.FALSE., error=error)
!       CALL ensure_array_size (lbuffer_2%mats(mult)%index,&
!            ub=MAXVAL (all_sizes(2, :, :)), nocopy=.TRUE.,&
!            special=.FALSE., error=error)
!    ENDDO
!    lbuffer_calc => lbuffer_1
!    lbuffer_comm => lbuffer_2
!    ! Setup buffers for right matrix
!    nblks = MAXVAL (all_sizes (4, 1, :))
!    nze = MAXVAL (all_sizes(3, 1, :))
!    CALL dbcsr_init (rbuffer_1, error)
!    CALL dbcsr_create (rbuffer_1,&
!         'rbuffer1'//right_matrix%name, right_matrix%dist,&
!         'N',&
!         right_matrix%row_blk_size,&
!         right_matrix%col_blk_size,&
!         0, 0, special=.FALSE., error=error)
!    CALL ensure_array_size (rbuffer_1%data,&
!         ub=MAXVAL (all_sizes(3, 1, :)), nocopy=.TRUE., special=.FALSE.,&
!         error=error)
!    CALL ensure_array_size (rbuffer_1%index,&
!         ub=MAXVAL (all_sizes(4, 1, :)), nocopy=.TRUE., special=.FALSE.,&
!         error=error)
!    !
!    CALL dbcsr_init (rbuffer_2, error)
!    CALL dbcsr_create (rbuffer_2,&
!         'rbuffer2'//right_matrix%name, right_matrix%dist,&
!         'N',&
!         right_matrix%row_blk_size,&
!         right_matrix%col_blk_size,&
!         0, 0, special=.FALSE., error=error)
!    CALL ensure_array_size (rbuffer_2%data,&
!         ub=MAXVAL (all_sizes(3, 1, :)), nocopy=.TRUE., special=.FALSE.,&
!         error=error)
!    CALL ensure_array_size (rbuffer_2%index,&
!         ub=MAXVAL (all_sizes(4, 1, :)), nocopy=.TRUE., special=.FALSE.,&
!         error=error)
!    ! Create windows for the right matrices
!    DBG 'RIGHT data win size', SIZE (right_matrix%data)
!    DBG 'RIGHT indx win size', SIZE (right_matrix%index)
!    CALL mp_window_create (right_data_window, right_matrix%data, gid=mp_group)
!    CALL mp_window_create (right_index_window, right_matrix%index, gid=mp_group)
!    rbuffer_calc => rbuffer_1
!    rbuffer_comm => rbuffer_2
!    !
!    right_nprows = dbcsr_mp_nprows (dbcsr_distribution_mp (&
!         right_matrix%dist))
!    right_src_p_calc = -7
!    left_src_p_calc = -7
!    ! Do the actual transfers and call local multiplication
!    ALLOCATE (my_waits(4, 0:right_nprows))
!    my_waits(:,:) = 0
!    k_index: DO right_prow_i = 0, right_nprows
!       ! Right matrix
!       right_prow = MOD (right_prow_i+myprow+mypcol, right_nprows)
!       DBG 'RIGHT_PROW from/to', right_prow_i, right_prow, myprow, right_nprows
!       my_waits(1, right_prow_i) = m_walltime()
!       !CALL dbcsr_sleep (10.0_dp)
!       CALL mp_window_fence (right_data_window,&
!            nothing_stored=.TRUE., no_put_follows=.TRUE.)
!       my_waits(1, right_prow_i) = m_walltime() - my_waits(1, right_prow_i)
!       my_waits(2, right_prow_i) = m_walltime()
!       !CALL dbcsr_sleep (10.0_dp)
!       CALL mp_window_fence (right_index_window,&
!            nothing_stored=.TRUE., no_put_follows=.TRUE.)
!       my_waits(2, right_prow_i) = m_walltime() - my_waits(2, right_prow_i)
!       ! No fetches in the last loop
!       xfer_case_right: IF (right_prow_i .LT. right_nprows) THEN
!          right_src_p_comm = right_pgrid (right_prow, mypcol)
!          IF (dbg) WRITE(*,'(1X,A,"(",I3,",",I3,")")')&
!               routineN//' Getting RIGHT matrix from pgrid coor ',&
!               right_prow, mypcol
!          IF (dbg) WRITE(*,'(1X,A,I4,A,"(",I3,",",I3,") sizes",I7,I7)')&
!               routineN//' Getting RIGHT matrix from',right_src_p_comm,&
!               ' pgrid coor ',right_prow, mypcol,&
!               all_sizes(3:4, 1, right_src_p_comm)
!          CALL mp_rma_get (right_data_window,&
!               right_src_p_comm, 0, all_sizes(3, 1, right_src_p_comm),&
!               rbuffer_comm%data)
!          CALL mp_rma_get (right_index_window,&
!               right_src_p_comm, 0, all_sizes(4, 1, right_src_p_comm),&
!               rbuffer_comm%index)
!       ENDIF xfer_case_right
!       ! No calculations in the first loop
!       calc_case_right: IF (right_prow_i .GT. 0) THEN
!          IF (all_sizes(4, 1, right_src_p_calc) .EQ. 0) THEN
!             DBG 'ZERO SIZED INDEX'
!          ENDIF
!          CALL dbcsr_repoint_index (rbuffer_calc, error=error)
!       ENDIF calc_case_right
!       ! Now we transfer/multiply with the left matrix
!       right_src_p_calc = right_src_p_comm
!       DO left_image = 1, left_images
!          my_waits(3, right_prow_i) = m_walltime()
!          !CALL dbcsr_sleep (10.0_dp)
!          CALL mp_window_fence (left_set%data_windows(left_image),&
!            nothing_stored=.TRUE., no_put_follows=.TRUE.)
!          my_waits(3, right_prow_i) = m_walltime() - my_waits(3, right_prow_i)
!          my_waits(4, right_prow_i) = m_walltime()
!          !CALL dbcsr_sleep (10.0_dp)
!          CALL mp_window_fence (left_set%index_windows(left_image),&
!            nothing_stored=.TRUE., no_put_follows=.TRUE.)
!       my_waits(4, right_prow_i) = m_walltime() - my_waits(4, right_prow_i)
!       ENDDO
!       left_prow = myprow / left_images
!       left_image = 1 + MOD (myprow, left_images)
!       xfer_case_left: IF (right_prow_i .LT. right_nprows) THEN
!          DO mult = 1, left_multiplicity
!             IF (dbg) WRITE(*,'(1X,A,I3,A,"(",I3,",",I3,")")')&
!                  routineN//' Getting LEFT matrix image',left_image,&
!                  ' pgrid coor ',left_prow, right_prow
!             IF (dbg) THEN
!                WRITE(*,*)routineP//'left pgrid shape', SHAPE (left_pgrid),';'
!                WRITE(*,*)routineP, left_prow,&
!                     (mult-1)*(left_nprows/left_multiplicity),&
!                     '=>',left_prow + (mult-1)*(left_nprows/left_multiplicity)
!             ENDIF
!             left_src_p_comm = left_pgrid (&
!                  left_prow + (mult-1)*(left_nprows/left_multiplicity),&
!                  right_prow)
!             IF (dbg) WRITE(*,'(1X,A,I4,A,I3,A,"(",I3,",",I3,")")')&
!                  routineN//' Getting LEFT matrix from',left_src_p_comm,&
!                  ' image',left_image,&
!                  ' pgrid coor ',left_prow, right_prow
!             CALL mp_rma_get (left_set%data_windows(left_image),&
!                  left_src_p_comm, 0,&
!                  all_sizes(1, left_image, left_src_p_comm),&
!                  lbuffer_comm%mats(mult)%data)
!             CALL mp_rma_get (left_set%index_windows(left_image),&
!                  left_src_p_comm, 0,&
!                  all_sizes(2, left_image, left_src_p_comm),&
!                  lbuffer_comm%mats(mult)%index)
!          ENDDO
!       ENDIF xfer_case_left
!       calc_case_left: IF (right_prow_i .GT. 0) THEN
!             !left_src_p_calc = left_pgrid (left_prow, right_prow-1)
!          DO mult = 1, left_multiplicity
!             CALL dbcsr_repoint_index (lbuffer_calc%mats(mult), error=error)
!             tstart = m_walltime ()
!             CALL dbcsr_nn_mult(lbuffer_calc%mats(mult), rbuffer_calc,&
!                  product_matrix, error, flop=flop_single)
!             flop_total = flop_total + flop_single
!             tstop = m_walltime ()
!             IF (dbg .AND. tstop-tstart .NE. 0) THEN
!                DBG 'Segment Local Multiplication time',tstop-tstart,'and',&
!                     (REAL(flop_total, dp)/1000000.0_dp)/(tstop-tstart),'MFLOPS'
!             ENDIF
!             trun = trun + (tstop - tstart)
!          ENDDO
!       ENDIF calc_case_left
!       left_src_p_calc = left_src_p_comm
!       CALL dbcsr_switch_sets (lbuffer_calc, lbuffer_comm)
!       CALL dbcsr_switch_m_ptrs (rbuffer_calc, rbuffer_comm)
!    ENDDO k_index
!    trun_t = m_walltime () - trun_t
!    IF (trun_t .EQ. 0) trun_t = 0.001_dp
!    IF (dbg .AND. detailed_timing) THEN
!       WRITE(*,'(1X,A,F9.3,A,F9.3,A,I9)')&
!            routineN//' Total Local Multiplication time',&
!            trun_t,' and',&
!            (REAL(flop_total, dp)/1000000.0_dp)/trun_t,' MFLOPS',flop_total
!    ENDIF
!    DO img = 1, left_images
!       CALL mp_window_destroy (left_set%data_windows(img))
!       CALL mp_window_destroy (left_set%index_windows(img))
!    ENDDO
!    DEALLOCATE (left_set%data_windows)
!    DEALLOCATE (left_set%index_windows)
!    CALL mp_window_destroy (right_data_window)
!    CALL mp_window_destroy (right_index_window)
!    !
!    CALL dbcsr_destroy (rbuffer_1, error)
!    CALL dbcsr_destroy (rbuffer_2, error)
!    CALL dbcsr_destroy_array (lbuffer_1, error)
!    CALL dbcsr_destroy_array (lbuffer_2, error)
!    DEALLOCATE (my_sizes)
!    !
!    IF (dbg) THEN
!       ALLOCATE (all_waits(4, 0:right_nprows, 0:numnodes-1))
!       CALL mp_allgather (my_waits, all_waits, mp_group)
!       logger => cp_error_get_logger (error)
!       IF (logger%para_env%mepos == logger%para_env%source) THEN
!          output_unit = cp_logger_get_default_unit_nr (logger)
!          WRITE(output_unit,'(A)')"Waiting times per place per step"
!          DO right_prow_i = 0, right_nprows
!             WRITE(output_unit,'(4(F9.4),5X,4(F9.4))')&
!                  (MAXVAL(all_waits(i,right_prow_i,:)), i=1, 4),&
!                  (SUM(all_waits(i,right_prow_i,:))/numnodes, i=1, 4)
!             WRITE(*,'(A,I4)')"Step ", right_prow_i
!             DO i = 0, numnodes-1
!                WRITE(*,'(I4,1X,4(F9.4))')i, all_waits(:, right_prow_i, i)
!             ENDDO
!          ENDDO
!       ENDIF
!    ENDIF
!    !
!    entries = btree_get_entries (product_matrix%w%tree)
!    CALL ensure_array_size (product_matrix%w%blk_p, ub=entries, error=error)
!    ALLOCATE (keys (entries))
!    CALL btree_delete (product_matrix%w%tree, keys, product_matrix%w%blk_p)
!    CALL ensure_array_size(product_matrix%w%row_i, ub=entries, error=error)
!    CALL ensure_array_size(product_matrix%w%col_i, ub=entries, error=error)
!    CALL dbcsr_unpack_i8_2i4 (keys, product_matrix%w%row_i,&
!         product_matrix%w%col_i)
!    DEALLOCATE (keys)
!    CALL dbcsr_finalize(product_matrix, .TRUE., .FALSE., error=error)
!    !
!    IF (PRESENT (flop)) flop = flop_total
!    CALL timestop (timing_handle)
!  END SUBROUTINE dbcsr_multiply_from_left


! *****************************************************************************
!> \brief Switches pointers between two matrices
!> \param[in,out] set1p, set2p
! *****************************************************************************
  SUBROUTINE dbcsr_switch_m_ptrs (m1p, m2p)
    TYPE(dbcsr_type)                         :: m1p, m2p

    TYPE(dbcsr_type)                         :: tmp_p

!   ---------------------------------------------------------------------------

    tmp_p = m1p
    m1p = m2p
    m2p = tmp_p
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


!! *****************************************************************************
!!> \brief Multiplies two DBCSR matrices, using the current local data.
!!> \param[in] left, right     left and right DBCSR matrices
!!> \param[in,out] product     resulting DBCSR product matrix
!!> \param[in,out] wm          DBCSR product working matrix
!!> \param[in,out] error       cp2k error
!!> \param[out] flop           (optional) number of effective double-precision
!!>                            floating point operations performed
!! *****************************************************************************
!  SUBROUTINE dbcsr_nn_mult(left, right, product, error, flop)
!    TYPE(dbcsr_type), INTENT(IN)             :: left, right
!    TYPE(dbcsr_type), INTENT(INOUT)          :: product
!    TYPE(cp_error_type), INTENT(INOUT)       :: error
!    INTEGER(KIND=int_8), INTENT(OUT), &
!      OPTIONAL                               :: flop
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_nn_mult', &
!      routineP = moduleN//':'//routineN
!    REAL(KIND=dp), PARAMETER                 :: alpha = 1.0_dp, &
!                                                beta_add = 1.0_dp, &
!                                                beta_new = 0.0_dp
!
!    CHARACTER                                :: a_type, b_type
!    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: left_data_z, right_data_z
!    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: left_data_c, right_data_c
!    INTEGER :: a_blk, a_col, a_col_log, a_col_size, a_col_size_logical, &
!      a_nze, a_row, a_row_log, a_row_size, a_row_size_logical, b_blk, b_col, &
!      b_col_size, b_nze, b_row, b_row_size, c_blk, c_blk_p, c_col, c_nze, &
!      c_row, lda, ldb, ldc, num_symmetries, stat, symmetry, timing_handle
!    INTEGER(KIND=int_4)                      :: offset
!    INTEGER(KIND=int_8)                      :: coordinate, lflop
!    INTEGER, ALLOCATABLE, DIMENSION(:)       :: c_blk_ps
!    INTEGER, DIMENSION(:), POINTER           :: left_col_blk_size, &
!                                                left_row_blk_size, &
!                                                right_col_blk_size, &
!                                                right_row_blk_size
!    LOGICAL                                  :: block_exists, tr_a, tr_b
!    REAL(KIND=dp)                            :: beta, tstart, tstop
!    REAL(KIND=dp), DIMENSION(:), POINTER     :: a_data, b_data, c_data, &
!                                                left_data_d, right_data_d
!    REAL(KIND=sp), DIMENSION(:), POINTER     :: left_data_r, right_data_r
!    TYPE(dbcsr_work_type)                    :: wm
!
!    !$ INTEGER                           :: ithread, nthreads
!    !$ INTEGER(KIND=omp_lock_kind)       :: data_lock, index_lock
!    DEBUG_HEADER
!!   ---------------------------------------------------------------------------
!    DEBUG_BODY
!    CALL timeset(routineN, timing_handle)
!    wm = product%wms(1)
!    DBG "Doing NN matrix multiplication",left%name,' x ',right%name
!    DBG "Pre-mult blk sizes, data", wm%lastblk, wm%datasize
!    ! c_blk_ps contains the block pointers in the product's current row
!    !.!ALLOCATE (c_blk_ps(product%nblkcols_total), stat=stat)
!    !IF (stat /= 0) CALL stop_memory (routineN, "c_blk_ps",&
!    !     product%nblkcols_total)
!    lflop = 0
!    num_symmetries = 1
!    left_row_blk_size => array_data (left%row_blk_size)
!    left_col_blk_size => array_data (left%col_blk_size)
!    right_row_blk_size => array_data (right%row_blk_size)
!    right_col_blk_size => array_data (right%col_blk_size)
!    SELECT CASE (left%data_type)
!    CASE ('r')
!       CALL dbcsr_get_data (left%data_area, left_data_r)
!       CALL dbcsr_get_data (right%data_area, right_data_r)
!    CASE ('R')
!       CALL dbcsr_get_data (left%data_area, left_data_d)
!       CALL dbcsr_get_data (right%data_area, right_data_d)
!    CASE ('c')
!       CALL dbcsr_get_data (left%data_area, left_data_c)
!       CALL dbcsr_get_data (right%data_area, right_data_c)
!    CASE ('C')
!       CALL dbcsr_get_data (left%data_area, left_data_z)
!       CALL dbcsr_get_data (right%data_area, right_data_z)
!    END SELECT
!    c_blk_p = wm%datasize + 1
!    DBG 'A rows:',left%nblkrows_total
!    IF (left%nblks .GT. 0 .AND. right%nblks .GT. 0) THEN
!    CALL dbcsr_data_ensure_size(wm%data_area,&
!         MAX(wm%datasize,left%nfullrows_total*right%nfullcols_total))
!    !CALL ensure_array_size(wm%data,&
!    !     ub=MAX(wm%datasize,left%nfullrows_total*right%nfullcols_total), error=error)
!       !CALL ensure_array_size(wm%data,&
!       !     ub=wm%datasize+a_row_size*right%nfullcols_total &
!       !     * MAX(1, (left%nblkrows_total-a_row)/4), error=error)
!    !$OMP PARALLEL DEFAULT (SHARED) PRIVATE(a_row, a_col, a_blk,&
!    !$OMP                                   b_row, b_col, b_blk,&
!    !$OMP                                   c_row, c_col, c_blk,&
!    !$OMP                                   a_row_log, a_col_log,&
!    !$OMP                                   a_row_size_logical,&
!    !$OMP                                   a_col_size_logical,&
!    !$OMP                                   a_row_size, a_col_size, a_nze,&
!    !$OMP                                   b_row_size, b_col_size, b_nze, c_nze,&
!    !$OMP                                   lda, ldb, ldc,&
!    !$OMP                                   tr_a, tr_b,&
!    !$OMP                                   c_blk_ps,&
!    !$OMP                                   a_data, b_data, c_data,&
!    !$OMP                                   beta,&
!    !$OMP                                   a_type, b_type,&
!    !$OMP                                   symmetry,&
!    !$OMP                                   coordinate,&
!    !$OMP                                   offset,&
!    !$OMP                                   block_exists,&
!    !$OMP                                   tstart, tstop,&
!    !$OMP                                   ithread)&
!    !$OMP          REDUCTION (+: lflop)
!    !
!    !$ CALL OMP_INIT_LOCK (index_lock)
!    !$ CALL OMP_INIT_LOCK (data_lock)
!    !$ ithread = OMP_GET_THREAD_NUM ()
!    !!$ WRITE (*,*)"Hello! I'm thread number", ithread
!    !$OMP SINGLE
!    !$ nthreads = OMP_GET_NUM_THREADS ()
!    !!$ WRITE (*,*)'And there are', nthreads, 'threads, from', ithread
!    !$OMP FLUSH (nthreads)
!    !$OMP END SINGLE
!    ! c_blk_ps contains the block pointers in the product's current row
!    ALLOCATE (c_blk_ps(product%nblkcols_total), stat=stat)
!    IF (stat /= 0) CALL stop_memory (routineN, "c_blk_ps",&
!         product%nblkcols_total)
!    !
!    !$OMP DO
!    DO a_row = 1, left%nblkrows_total
!       a_row_log = a_row
!       DBG 'A row', a_row,' with #cols=',left%row_p(a_row+1)-left%row_p(a_row)
!       a_row_size = left_row_blk_size(a_row) ! stored
!       c_blk_ps(:) = 0
!       DO a_blk = left%row_p(a_row)+1, left%row_p(a_row+1)
!          DO symmetry = 1, num_symmetries
!             a_col = left%col_i(a_blk)
!             IF (a_col .EQ. a_row .AND. symmetry .GT. 1) CYCLE
!             a_col_size = left_col_blk_size(a_col) !stored
!             a_row_log = a_row
!             a_col_log = a_col
!             a_row_size_logical = a_row_size
!             a_col_size_logical = a_col_size
!             lda = a_row_size
!             ! For regular symmetric matrices
!             IF (symmetry .GT. 1) THEN
!                CALL swap (a_row_log, a_col_log)
!                CALL swap (a_row_size_logical, a_col_size_logical)
!                tr_a = .TRUE.
!                a_type = 'T'
!             ELSE
!                tr_a = .FALSE.
!                a_type = 'N'
!             ENDIF
!             DBG 'A col', a_col, '=>',a_col_log
!             b_row = a_col_log
!             c_row = a_row_log
!             a_nze = a_row_size * a_col_size
!             a_data => left_data_d(ABS(left%blk_p(a_blk))&
!                  :ABS(left%blk_p(a_blk))+a_nze-1)
!             !a_data => left%data(ABS(left%blk_p(a_blk))&
!             !     :ABS(left%blk_p(a_blk))+a_nze-1)
!             b_row_size = a_col_size_logical
!             ! For index-augmented matrices.
!             IF (left%blk_p(a_blk) .LT. 0) THEN
!                lda = a_col_size
!                a_type = 'T'
!             ENDIF
!             DO b_blk = right%row_p(b_row)+1, right%row_p(b_row+1)
!                b_col = right%col_i(b_blk)
!                c_col = b_col
!                DBG 'B col', b_col
!                b_col_size = right_col_blk_size(b_col)
!                b_nze = b_row_size * b_col_size
!                b_data => right_data_d(&
!                     ABS(right%blk_p(b_blk)):ABS(right%blk_p(b_blk))+b_nze-1)
!                !b_data => right%data(&
!                !     ABS(right%blk_p(b_blk)):ABS(right%blk_p(b_blk))+b_nze-1)
!                IF (right%blk_p(b_blk) .LT. 0) THEN
!                   b_type = 'T'
!                   ldb = b_col_size
!                ELSE
!                   b_type = 'N'
!                   ldb = b_row_size
!                ENDIF
!                c_nze = a_row_size * b_col_size
!                ! If this c column does not yet exist, "allocate" it.
!                coordinate = make_coordinate_tuple (c_row, c_col)
!                IF (c_blk_ps(c_col) .EQ. 0 .AND. .NOT. tr_a) THEN
!                   CALL btree_find (wm%tree, coordinate, offset, block_exists)
!                ELSE
!                   block_exists = .TRUE.
!                   offset = c_blk_ps(c_col)
!                ENDIF
!                IF (.NOT. block_exists) THEN
!                   !.!CALL add_work_coordinate(wm, c_row, c_col, c_blk_p, error=error)
!                   !$ CALL OMP_SET_LOCK (data_lock)
!                   offset = c_blk_p
!                   c_blk_p = c_blk_p + c_nze
!                   wm%datasize = wm%datasize + c_nze
!                   wm%lastblk = wm%lastblk + 1
!                   !$ CALL OMP_UNSET_LOCK (data_lock)
!                   CALL btree_add (wm%tree, coordinate, offset, block_exists,&
!                        existing_value = offset)
!                   !CALL btree_verify (index)
!                   IF (block_exists) THEN
!                      IF (dbg) WRITE(*,*)routineP//' EXISTS'
!                      ! We can not rewind c_blk_p!
!                      !$OMP ATOMIC
!                      wm%datasize = wm%datasize - c_nze
!                      !$OMP ATOMIC
!                      wm%lastblk = wm%lastblk - 1
!                   ENDIF
!                   c_blk_ps(c_col) = offset
!                   beta = beta_new
!                ELSE
!                   beta = beta_add
!                ENDIF
!                !offset = 1 + ithread * 4096
!                !.!c_data => wm%data(c_blk_ps(c_col):c_blk_ps(c_col)+c_nze-1)
!                DBG 'offset:',offset,offset+c_nze-1,'of',wm%datasize,'ex',block_exists
!                !% CALL OMP_SET_LOCK (data_lock)
!                c_data => pointer_view (&
!                     dbcsr_get_data_p (wm%data_area, 0.0_dp),&
!                     offset, offset+c_nze-1)
!                !c_data => wm%data(offset:offset+c_nze-1)
!                !% CALL OMP_UNSET_LOCK (data_lock)
!                IF (dbg) WRITE(*,&
!                     "(1X,A,A1,1X,A1,A,I6,A,F4.1,' (',I3,',',I3,')=(',I3,',',I3,')x(',I3,',',I3,')')")&
!                                !routineN//' DGEMM ',a_type,b_type,' at',c_blk_ps(c_col),' beta=', beta,&
!                     routineN//' DGEMM ',a_type,b_type,' at',offset,' beta=', beta,&
!                     c_row,c_col,a_row,a_col,b_row,b_col
!                CALL DGEMM(a_type, b_type,&
!                     a_row_size_logical, b_col_size, a_col_size_logical,&
!                     alpha,&
!                     a_data, lda,&
!                     b_data, ldb, beta,&
!                     c_data, a_row_size_logical)
!                !write(6,*)'NEWTIME',a_row_size,b_col_size,a_col_size,t2-t1,&
!                !     REAL(t2-t1)/(REAL(2*a_row_size*b_col_size*a_col_size))
!                lflop = lflop + a_row_size * b_col_size * a_col_size * 2
!             ENDDO ! b
!          ENDDO ! symmetry
!       ENDDO ! a_col
!    ENDDO ! a_row
!
!!    !$OMP END DO
!    !$ CALL OMP_DESTROY_LOCK (index_lock)
!    !$ CALL OMP_DESTROY_LOCK (data_lock)
!    DEALLOCATE (c_blk_ps)
!    !$OMP END PARALLEL
!    ELSE
!       DBG 'One of the matrices is empty!'
!    ENDIF
!    IF (dbg) CALL btree_verify (wm%tree)
!    !DEALLOCATE (c_blk_ps)
!    !WRITE(*,*)'NEW TOTAL TIME', nanotot, REAL(nanotot)/(2.3*REAL(10**6)),'ms'
!    DBG "Post-mult blk sizes, data", wm%lastblk, wm%datasize
!    !dbg = .TRUE.
!    IF (dbg .AND. detailed_timing) THEN
!       tstop = m_walltime ()
!       IF (tstop-tstart .NE. 0) THEN
!          DBG tstop-tstart,'s',&
!               (REAL(lflop, dp)/1000000.0_dp)/(tstop-tstart),'MFLOPS'
!       ENDIF
!    ENDIF
!    IF (PRESENT (flop)) flop = lflop
!    product%wms(1) = wm
!    CALL timestop(timing_handle)
!  END SUBROUTINE dbcsr_nn_mult


  ELEMENTAL FUNCTION make_coordinate_tuple (most, least) RESULT (tuple)
    INTEGER, INTENT(IN)                      :: most, least
    INTEGER(KIND=int_8)                      :: tuple

!tuple = IOR (ISHFT (most, 32), least)

    tuple = most
    tuple = IOR (ISHFT(tuple, 32), INT(least, int_8))
  END FUNCTION make_coordinate_tuple


!! *****************************************************************************
!!> \brief Multiplies two DBCSR matrices with square processor grids.
!!> \param[in] left_set             set of left images
!!> \param[in] right_matrix         right matrix
!!> \param[out] product             DBCSR product workmatrix
!!> \param[in,out] error            cp2k error
!!> \param[out] flop                (optional) effective flop
!! *****************************************************************************
!  SUBROUTINE dbcsr_multiply_from_left_sqpg (left_set, right_matrix,&
!       product_matrix, error, flop)
!    TYPE(dbcsr_array_type), INTENT(INOUT)    :: left_set
!    TYPE(dbcsr_type), INTENT(IN)             :: right_matrix
!    TYPE(dbcsr_type), INTENT(INOUT)          :: product_matrix
!    TYPE(cp_error_type), INTENT(INOUT)       :: error
!    INTEGER(KIND=int_8), INTENT(OUT), &
!      OPTIONAL                               :: flop
!
!    CHARACTER(len=*), PARAMETER :: &
!      routineN = 'dbcsr_multiply_from_left_sqpg', &
!      routineP = moduleN//':'//routineN
!
!    INTEGER :: debug_level, i, img, left_dst_p, left_dst_pcol, left_dst_prow, &
!      left_image, left_images, left_multiplicity, left_mypcol, left_myprow, &
!      left_npcols, left_nprows, left_p_recv, left_p_send, left_src_p, &
!      left_src_pcol, left_src_prow, mp_group, mult, mynode, mypcol, myprow, &
!      nblks, npcols, nprows, numnodes, nze, right_data_rr, right_data_sr, &
!      right_dst_p, right_dst_prow, right_index_rr, right_index_sr, &
!      right_mypcol, right_myprow, right_npcols, right_nprows, right_p_recv, &
!      right_p_send, right_prow_i, right_src_p, right_src_prow, tag, &
!      timing_handle
!    INTEGER(KIND=int_8)                      :: flop_single, flop_total
!    INTEGER, ALLOCATABLE, DIMENSION(:)       :: left_data_rr, left_data_sr, &
!                                                left_index_rr, left_index_sr
!    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: my_sizes
!    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: all_sizes
!    INTEGER, DIMENSION(:), POINTER           :: left_index_rp, left_index_sp, &
!                                                right_index_rp, right_index_sp
!    INTEGER, DIMENSION(:, :), POINTER        :: left_pgrid, product_pgrid, &
!                                                right_pgrid
!    REAL(KIND=dp)                            :: t_all, t_dgemm, trun, trun_t, &
!                                                tstart, tstop
!    REAL(KIND=dp), ALLOCATABLE, &
!      DIMENSION(:, :)                        :: my_waits
!    REAL(KIND=dp), DIMENSION(:), POINTER     :: left_data_rp, left_data_sp, &
!                                                right_data_rp, right_data_sp
!    TYPE(dbcsr_array_type), POINTER          :: lbuffer_calc, lbuffer_comm
!    TYPE(dbcsr_array_type), TARGET           :: lbuffer_1, lbuffer_2
!    TYPE(dbcsr_mp_obj)                       :: left_mp, mp_obj, right_mp
!    TYPE(dbcsr_type), POINTER                :: rbuffer_calc, rbuffer_comm
!    TYPE(dbcsr_type), TARGET                 :: rbuffer_1, rbuffer_2
!
!    DEBUG_HEADER
!!   ---------------------------------------------------------------------------
!    DEBUG_BODY
!    debug_level=0
!    IF (dbg) debug_level = 2
!    CALL timeset (routineN, timing_handle)
!    t_all = 0.0_dp
!    t_dgemm = 0.0_dp
!    CALL btree_new (product_matrix%w%tree)
!    tstart = 0.0_dp ; tstop = 0.0_dp ; trun = 0.0_dp
!    flop_total = 0
!    trun_t = m_walltime ()
!    left_images = SIZE (left_set%mats)
!    CALL bcsr_assert (debug_level,&
!         left_images .EQ. left_set%image_dist%row_decimation, routineP,&
!         "Non-compatible images")
!    left_multiplicity = left_set%image_dist%row_multiplicity
!    mp_obj = dbcsr_distribution_mp (product_matrix%dist)
!    mp_group = dbcsr_mp_group (mp_obj)
!    left_mp = dbcsr_distribution_mp (left_set%image_dist%main)
!    right_mp = dbcsr_distribution_mp (right_matrix%dist)
!    numnodes = dbcsr_mp_numnodes (mp_obj)
!    mynode = dbcsr_mp_mynode (mp_obj)
!    myprow = dbcsr_mp_myprow (mp_obj)
!    mypcol = dbcsr_mp_mypcol (mp_obj)
!    npcols = dbcsr_mp_npcols (mp_obj)
!    nprows = dbcsr_mp_nprows (mp_obj)
!    product_pgrid => dbcsr_mp_pgrid (mp_obj)
!    right_myprow = dbcsr_mp_myprow (right_mp)
!    right_mypcol = dbcsr_mp_mypcol (right_mp)
!    right_nprows = dbcsr_mp_nprows (right_mp)
!    right_npcols = dbcsr_mp_npcols (right_mp)
!    right_pgrid => dbcsr_mp_pgrid (right_mp)
!    left_pgrid => dbcsr_mp_pgrid (left_mp)
!    left_myprow = dbcsr_mp_myprow (left_mp)
!    left_mypcol = dbcsr_mp_mypcol (left_mp)
!    left_nprows = dbcsr_mp_nprows (left_mp)
!    left_npcols = dbcsr_mp_npcols (left_mp)
!    CALL ensure_array_size(product_matrix%w%row_i,&
!         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
!    CALL ensure_array_size(product_matrix%w%col_i,&
!         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
!    CALL ensure_array_size(product_matrix%w%blk_p,&
!         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
!    CALL ensure_array_size(product_matrix%w%blk_d,&
!         ub=product_matrix%nblkcols_total*product_matrix%nblkrows_total, error=error)
!    DBG 'product pgrid shape', SHAPE (dbcsr_mp_pgrid (mp_obj))
!    DBG 'left pgrid shape', SHAPE (left_pgrid),';',left_nprows,left_images
!    DBG 'right pgrid shape', SHAPE (right_pgrid)
!    DBG 'left images, multiplicity', left_images, left_multiplicity
!    CALL bcsr_assert (debug_level, left_npcols .EQ. right_nprows, routineP,&
!         "Left and right matrix processor grid dimension mismatch")
!    CALL bcsr_assert (debug_level,&
!         left_multiplicity * npcols .EQ. left_nprows, routineP,&
!         "Mulitplicity mismatch and processer grid.")
!    DBG 'left pgrid =', left_pgrid
!    DBG 'right pgrid =', right_pgrid
!    ALLOCATE (left_set%data_windows(left_images))
!    ALLOCATE (left_set%index_windows(left_images))
!    ALLOCATE (my_sizes(4, left_images))
!    my_sizes(:,:) = 0
!    DO img = 1, left_images
!       my_sizes(1, img) = SIZE (left_set%mats(img)%data)
!       my_sizes(2, img) = SIZE (left_set%mats(img)%index)
!    ENDDO
!    my_sizes(3, 1) = SIZE (right_matrix%data)
!    my_sizes(4, 1) = SIZE (right_matrix%index)
!    ALLOCATE (all_sizes(4, LBOUND(my_sizes,2):UBOUND(my_sizes,2), 0:numnodes-1))
!    DO i = 1, 4
!       DBG 'my sizes', my_sizes(i,:)
!    ENDDO
!    CALL mp_allgather(my_sizes, all_sizes, mp_group)
!    nblks = MAXVAL (all_sizes (2, :, :))
!    nze = MAXVAL (all_sizes(1, :, :))
!    DBG 'data sizes', all_sizes
!    !! Setup the buffer matrices
!    !DO img = 1, left_images
!    !   CALL mp_window_create (left_set%data_windows(img),&
!    !        left_set%mats(img)%data, gid=mp_group)
!    !   CALL mp_window_create (left_set%index_windows(img),&
!    !        left_set%mats(img)%index, gid=mp_group)
!    !ENDDO
!    ALLOCATE (lbuffer_1%mats(left_multiplicity))
!    ALLOCATE (lbuffer_2%mats(left_multiplicity))
!    ALLOCATE (left_data_sr(left_multiplicity), left_data_rr(left_multiplicity))
!    ALLOCATE (left_index_sr(left_multiplicity), left_index_rr(left_multiplicity))
!    DO mult = 1, left_multiplicity
!       CALL dbcsr_init (lbuffer_1%mats(mult), error)
!       CALL dbcsr_create (lbuffer_1%mats(mult),&
!            'lbufferl1'//left_set%mats(1)%m%name, left_set%mats(1)%m%dist, 'N',&
!            left_set%mats(1)%m%row_blk_size, left_set%mats(1)%m%col_blk_size,&
!            nblks, nze, special=.FALSE., error=error)
!       CALL ensure_array_size (lbuffer_1%mats(mult)%data,&
!            ub=MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.,&
!            special=.FALSE., error=error)
!       CALL ensure_array_size (lbuffer_1%mats(mult)%index,&
!            ub=MAXVAL (all_sizes(2, :, :)), nocopy=.TRUE.,&
!            special=.FALSE., error=error)
!       !
!       CALL dbcsr_init (lbuffer_2%mats(mult), error)
!       CALL dbcsr_create (lbuffer_2%mats(mult),&
!            'lbuffer2'//left_set%mats(1)%m%name, left_set%mats(1)%m%dist, 'N',&
!            left_set%mats(1)%m%row_blk_size, left_set%mats(1)%m%col_blk_size,&
!            nblks, nze, special=.FALSE., error=error)
!       CALL ensure_array_size (lbuffer_2%mats(mult)%data,&
!            ub=MAXVAL (all_sizes(1, :, :)), nocopy=.TRUE.,&
!            special=.FALSE., error=error)
!       CALL ensure_array_size (lbuffer_2%mats(mult)%index,&
!            ub=MAXVAL (all_sizes(2, :, :)), nocopy=.TRUE.,&
!            special=.FALSE., error=error)
!    ENDDO
!    lbuffer_calc => lbuffer_1
!    lbuffer_comm => lbuffer_2
!    ! Setup buffers for right matrix
!    nblks = MAXVAL (all_sizes (4, 1, :))
!    nze = MAXVAL (all_sizes(3, 1, :))
!    CALL dbcsr_init (rbuffer_1, error)
!    CALL dbcsr_create (rbuffer_1,&
!         'rbuffer1'//right_matrix%name, right_matrix%dist,&
!         'N',&
!         right_matrix%row_blk_size,&
!         right_matrix%col_blk_size,&
!         0, 0, special=.FALSE., error=error)
!    CALL ensure_array_size (rbuffer_1%data,&
!         ub=MAXVAL (all_sizes(3, 1, :)), nocopy=.TRUE., special=.FALSE.,&
!         error=error)
!    CALL ensure_array_size (rbuffer_1%index,&
!         ub=MAXVAL (all_sizes(4, 1, :)), nocopy=.TRUE., special=.FALSE.,&
!         error=error)
!    !
!    CALL dbcsr_init (rbuffer_2, error)
!    CALL dbcsr_create (rbuffer_2,&
!         'rbuffer2'//right_matrix%name, right_matrix%dist,&
!         'N',&
!         right_matrix%row_blk_size,&
!         right_matrix%col_blk_size,&
!         0, 0, special=.FALSE., error=error)
!    CALL ensure_array_size (rbuffer_2%data,&
!         ub=MAXVAL (all_sizes(3, 1, :)), nocopy=.TRUE., special=.FALSE.,&
!         error=error)
!    CALL ensure_array_size (rbuffer_2%index,&
!         ub=MAXVAL (all_sizes(4, 1, :)), nocopy=.TRUE., special=.FALSE.,&
!         error=error)
!    ! Create windows for the right matrices
!    DBG 'RIGHT data win size', SIZE (right_matrix%data)
!    DBG 'RIGHT indx win size', SIZE (right_matrix%index)
!    !CALL mp_window_create (right_data_window, right_matrix%data, gid=mp_group)
!    !CALL mp_window_create (right_index_window, right_matrix%index, gid=mp_group)
!    rbuffer_calc => rbuffer_1
!    rbuffer_comm => rbuffer_2
!    ! To which processor we normally send
!    right_p_send = right_pgrid (MOD (right_nprows&
!         + right_myprow - 1, right_nprows),&
!         right_mypcol)
!    ! From which processor we normally receive.
!    right_p_recv = right_pgrid (MOD (&
!         right_myprow + 1, right_nprows),&
!         right_mypcol)
!    !
!    ! Do the actual transfers and call local multiplication
!    ALLOCATE (my_waits(4, 0:right_nprows))
!    my_waits(:,:) = 0    
!    k_index: DO right_prow_i = 0, right_nprows
!       ! Right matrix
!       !
!       ! right_src_prow is the matrix "source" processor (as if a
!       ! direct mpi_get) would be called.
!       right_src_prow = MOD ( right_prow_i &
!            + right_myprow + right_mypcol/left_multiplicity, right_nprows)
!       ! right_dst_prow is the matrix "destination" processor (i.e.,
!       ! the matrix data that we currently have)
!       right_dst_prow = MOD ( right_prow_i - 1 + 4*right_nprows&
!            + right_myprow + right_mypcol/left_multiplicity, right_nprows)
!
!       DBG 'RIGHT_PROW from/to', right_prow_i, right_src_prow, right_dst_prow,&
!            myprow, right_nprows
!       ! No wait in first loop iteration.
!       IF (right_prow_i .GT. 0) THEN
!          DBG 'It', right_prow_i,'waiting for right'
!          my_waits(1, right_prow_i) = m_walltime()
!          CALL mp_wait (right_data_sr)
!          CALL mp_wait (right_index_sr)
!          my_waits(1, right_prow_i) = m_walltime() - my_waits(1, right_prow_i)
!          my_waits(2, right_prow_i) = m_walltime()
!          CALL mp_wait (right_data_rr)
!          CALL mp_wait (right_index_rr)
!          my_waits(2, right_prow_i) = m_walltime() - my_waits(2, right_prow_i)
!       ENDIF
!       ! No fetches in the last loop
!       xfer_case_right: IF (right_prow_i .LT. right_nprows) THEN
!          ! First time the meaning of "send" is different
!          IF (right_prow_i .EQ. 0) THEN
!             right_dst_prow = MOD (4*right_nprows&
!                  + right_myprow - right_mypcol/left_multiplicity, right_nprows)
!          ENDIF
!          right_src_p = right_pgrid (right_src_prow, right_mypcol)
!          right_dst_p = right_pgrid (right_dst_prow, right_mypcol)
!          IF (dbg) WRITE(*,'(1X,A,"(",I3,",",I3,")")')&
!               routineN//' Getting RIGHT matrix pgrid coor ',&
!               right_src_prow, mypcol
!          IF (dbg) WRITE(*,'(1X,A,"(",I3,",",I3,")")')&
!               routineN//' Sending RIGHT matrix pgrid coor ',&
!               right_dst_prow, mypcol
!          IF (dbg) WRITE(*,'(1X,A,I4,A,"(",I3,",",I3,") sizes",I7,I7)')&
!               routineN//' Getting RIGHT matrix from',right_src_p,&
!               ' pgrid coor ',right_src_prow, mypcol
!          IF (dbg) WRITE(*,'(1X,A,I4,A,"(",I3,",",I3,") sizes",I7,I7)')&
!               routineN//' Sending RIGHT matrix to  ',right_dst_p,&
!               ' pgrid coor ',right_dst_prow, mypcol
!          IF (dbg) WRITE(*,'(1X,A,2(I4,1X),A,2(I4,1X))')routineN//" Right src ",&
!               right_src_p, all_sizes(3, 1, right_src_p),&
!               " dst", right_dst_p,SIZE(right_matrix%data)
!          right_data_rp => rbuffer_comm%data(1:all_sizes(3,1,right_src_p))
!          right_index_rp => rbuffer_comm%index(1:all_sizes(4,1,right_src_p))
!          right_data_sp => rbuffer_calc%data(1:all_sizes(3,1,right_dst_p))
!          right_index_sp => rbuffer_calc%index(1:all_sizes(4,1,right_dst_p))
!          IF (right_prow_i .EQ. 0) THEN
!             CALL mp_isendrecv (right_matrix%data, right_dst_p,&
!                  right_data_rp, right_src_p, &
!                  mp_group, right_data_sr, right_data_rr, tag=3)
!             CALL mp_isendrecv (right_matrix%index, right_dst_p, &
!                  right_index_rp,  right_src_p,&
!                  mp_group, right_index_sr, right_index_rr, tag=4)
!          ELSE
!             CALL mp_isendrecv (right_data_sp, right_p_send,&
!                  right_data_rp, right_p_recv, &
!                  mp_group, right_data_sr, right_data_rr, tag=3)
!             CALL mp_isendrecv (right_index_sp, right_p_send, &
!                  right_index_rp,  right_p_recv,&
!                  mp_group, right_index_sr, right_index_rr, tag=4)
!          ENDIF
!       ENDIF xfer_case_right
!       ! No calculations in the first loop
!       calc_case_right: IF (right_prow_i .GT. 0) THEN
!          CALL dbcsr_repoint_index (rbuffer_calc, error=error)
!       ENDIF calc_case_right
!       ! No waits in first loop iteration.
!       IF (right_prow_i .GT. 0) THEN
!          ! Now we transfer/multiply with the left matrix
!          DO mult = 1, left_multiplicity
!             DBG 'It', right_prow_i,'waiting for left, multiplicity', mult
!             my_waits(3, right_prow_i) = m_walltime()
!             CALL mp_wait (left_data_sr(mult))
!             CALL mp_wait (left_index_sr(mult))
!             my_waits(3, right_prow_i) = m_walltime() - my_waits(3, right_prow_i)
!             my_waits(4, right_prow_i) = m_walltime()
!             CALL mp_wait (left_data_rr(mult))
!             CALL mp_wait (left_index_rr(mult))
!             my_waits(4, right_prow_i) = m_walltime() - my_waits(4, right_prow_i)
!          ENDDO
!       ENDIF
!       left_image = 1 + MOD (right_myprow, left_images)
!       xfer_case_left: IF (right_prow_i .LT. right_nprows) THEN
!          DO mult = 1, left_multiplicity
!             left_src_prow = (myprow / left_images)
!             left_src_prow = left_src_prow&
!                  + (mult-1)*(left_nprows/left_multiplicity)
!             left_dst_prow = (myprow / left_images)
!             left_dst_prow = left_dst_prow&
!                  + (mult-1)*(left_nprows/left_multiplicity)
!             !left_dst_prow = left_src_prow
!             !!+...left_myprow/left_multiplicity
!             left_src_pcol = MOD ( right_prow_i&
!                  + left_mypcol + (myprow/left_images), left_npcols)
!             left_dst_pcol = MOD ( right_prow_i-1&
!                  + left_mypcol + (myprow/left_images), left_npcols)
!             !left_dst_pcol = MOD (-right_prow_i + 4*left_npcols&
!             !     + left_mypcol - (myprow/left_images), left_npcols)
!             IF (right_prow_i .EQ. 0) THEN
!                left_dst_pcol = MOD (4*left_npcols&
!                     + left_mypcol - (myprow/left_images), left_npcols)
!             ENDIF
!             !
!             IF (dbg) WRITE(*,'(1X,A,I3,A,"(",I3,",",I3,")")')&
!                  routineN//' Getting LEFT matrix image',left_image,&
!                  ' pgrid coor ',left_src_prow, left_src_pcol
!             IF (dbg) WRITE(*,'(1X,A,I3,A,"(",I3,",",I3,")")')&
!                  routineN//' Sending LEFT matrix image',left_image,&
!                  ' pgrid coor ',left_dst_prow, left_dst_pcol
!             left_src_p = left_pgrid (left_src_prow, left_src_pcol)
!             left_dst_p = left_pgrid (left_dst_prow, left_dst_pcol)
!             !
!             left_p_send = left_pgrid (left_src_prow,&
!                  MOD (left_mypcol - 1+left_npcols, left_npcols))
!             left_p_recv = left_pgrid (left_src_prow,&
!                  MOD (left_mypcol + 1, left_npcols))
!             !
!             !
!             IF (dbg) WRITE(*,'(1X,A,I4,A,I3,A,"(",I3,",",I3,")")')&
!                  routineN//' Getting LEFT matrix from',left_src_p,&
!                  ' image',left_image,&
!                  ' pgrid coor ',left_src_prow, left_src_pcol
!             IF (dbg) WRITE(*,'(1X,A,I4,A,I3,A,"(",I3,",",I3,")")')&
!                  routineN//' Sending LEFT matrix to  ',left_dst_p,&
!                  ' image',left_image,&
!                  ' pgrid coor ',left_dst_prow, left_dst_pcol
!             IF (dbg) WRITE(*,'(1X,A,I4,A,I4,A,I4)')routineN//"  Left mult. ",&
!                  mult," src ", left_src_p, " dst ",right_dst_p
!             left_data_rp => lbuffer_comm%mats(mult)%data(&
!                  1:all_sizes(1, left_image, left_src_p))
!             left_index_rp => lbuffer_comm%mats(mult)%index(&
!                  1:all_sizes(2, left_image, left_src_p))
!             left_data_sp => lbuffer_calc%mats(mult)%data(&
!                  1:all_sizes(1, left_image, left_dst_p))
!             left_index_sp => lbuffer_calc%mats(mult)%index(&
!                  1:all_sizes(2, left_image, left_dst_p))
!             IF (right_prow_i .EQ. 0) THEN
!                CALL mp_isendrecv (left_set%mats(mult)%data, left_dst_p,&
!                     left_data_rp, left_src_p, mp_group,&
!                     left_data_sr(mult), left_data_rr (mult), tag=1)
!                CALL mp_isendrecv (left_set%mats(mult)%index, left_dst_p,&
!                     left_index_rp, left_src_p, mp_group,&
!                     left_index_sr(mult), left_index_rr (mult), tag=2)
!             ELSE
!                CALL mp_isendrecv (left_data_sp, left_p_send,&
!                     left_data_rp, left_p_recv, mp_group,&
!                     left_data_sr(mult), left_data_rr (mult), tag=1)
!                CALL mp_isendrecv (left_index_sp, left_p_send,&
!                     left_index_rp, left_p_recv, mp_group,&
!                     left_index_sr(mult), left_index_rr (mult), tag=2)
!             ENDIF
!          ENDDO
!       ENDIF xfer_case_left
!       calc_case_left: IF (right_prow_i .GT. 0) THEN
!          DO mult = 1, left_multiplicity
!             CALL dbcsr_repoint_index (lbuffer_calc%mats(mult), error=error)
!             tstart = m_walltime ()
!             CALL dbcsr_nn_mult_lin(lbuffer_calc%mats(mult), rbuffer_calc,&
!                  product_matrix, error, flop=flop_single,&
!                  t_all=t_all, t_dgemm=t_dgemm)
!             flop_total = flop_total + flop_single
!             tstop = m_walltime ()
!             IF (tstop-tstart .NE. 0 .AND. detailed_timing) THEN
!                DBG 'Segment Local Multiplication time',tstop-tstart,'and',&
!                     (REAL(flop_total, dp)/1000000.0_dp)/(tstop-tstart),'MFLOPS'
!             ENDIF
!             trun = trun + (tstop - tstart)
!          ENDDO
!       ENDIF calc_case_left
!       CALL dbcsr_switch_sets (lbuffer_calc, lbuffer_comm)
!       CALL dbcsr_switch_m_ptrs (rbuffer_calc, rbuffer_comm)
!    ENDDO k_index
!    trun_t = m_walltime () - trun_t
!    IF (trun_t .EQ. 0) trun_t = 0.001_dp
!    IF (dbg .AND. detailed_timing) THEN
!       WRITE(*,'(1X,A,F9.3,A,F9.3,A,I9)')&
!            routineN//' Total Local Multiplication time',&
!            trun_t,' and',&
!            (REAL(flop_total, dp)/1000000.0_dp)/trun_t,' MFLOPS',flop_total
!    ENDIF
!    !DO img = 1, left_images
!    !   CALL mp_window_destroy (left_set%data_windows(img))
!    !   CALL mp_window_destroy (left_set%index_windows(img))
!    !ENDDO
!    DEALLOCATE (left_set%data_windows)
!    DEALLOCATE (left_set%index_windows)
!    !CALL mp_window_destroy (right_data_window)
!    !CALL mp_window_destroy (right_index_window)
!    !
!    CALL dbcsr_destroy (rbuffer_1, error)
!    CALL dbcsr_destroy (rbuffer_2, error)
!    CALL dbcsr_destroy_array (lbuffer_1, error)
!    CALL dbcsr_destroy_array (lbuffer_2, error)
!    DEALLOCATE (my_sizes)
!    !
!    !ALLOCATE (all_waits(4, 0:right_nprows, 0:numnodes-1))
!    !CALL mp_allgather (my_waits, all_waits, mp_group)
!    !logger => cp_error_get_logger (error)
!    !IF (logger%para_env%mepos == logger%para_env%source) THEN
!    !   output_unit = cp_logger_get_default_unit_nr (logger)
!    !   WRITE(output_unit,'(A)')"Waiting times per place per step"
!    !   DO right_prow_i = 0, right_nprows
!    !      WRITE(output_unit,'(4(F9.4),5X,4(F9.4))')&
!    !           (MAXVAL(all_waits(i,right_prow_i,:)), i=1, 4),&
!    !           (SUM(all_waits(i,right_prow_i,:))/numnodes, i=1, 4)
!    !      WRITE(*,'(A,I4)')"Step ", right_prow_i
!    !      DO i = 0, numnodes-1
!    !         WRITE(*,'(I4,1X,4(F9.4))')i, all_waits(:, right_prow_i, i)
!    !      ENDDO
!    !   ENDDO
!    !ENDIF
!    !
!    !entries = btree_get_entries (product_matrix%w%tree)
!    !CALL ensure_array_size (product_matrix%w%blk_p, ub=entries, error=error)
!    !ALLOCATE (keys (entries))
!    !CALL btree_delete (product_matrix%w%tree, keys, product_matrix%w%blk_p)
!    !CALL ensure_array_size(product_matrix%w%row_i, ub=entries, error=error)
!    !CALL ensure_array_size(product_matrix%w%col_i, ub=entries, error=error)
!    !CALL dbcsr_unpack_i8_2i4 (keys, product_matrix%w%row_i,&
!    !     product_matrix%w%col_i)
!    !DEALLOCATE (keys)
!    !CALL mp_sync (mp_group)
!    CALL dbcsr_finalize(product_matrix, .TRUE., .FALSE., error=error)
!    !
!    IF (PRESENT (flop)) flop = flop_total
!    CALL timestop (timing_handle)
!  END SUBROUTINE dbcsr_multiply_from_left_sqpg


! *****************************************************************************
!> \brief Multiplies two DBCSR matrices, using the current local data.
!>
!> Attempts to block consecutive columns in rows of the right matrix.
!> \param[in] left, right     left and right DBCSR matrices
!> \param[in,out] product     resulting DBCSR product matrix
!> \param[in] submatrix_box   select submatrix
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix, default is no
!> \param[in] alpha_d              (optional) product multiplication factor
!> \param[in] beta_d               (optional) existing target multiplication
!>                                 factor
!> \param[out] flop           (optional) number of effective double-precision
!>                            floating point operations performed
! *****************************************************************************
  SUBROUTINE dbcsr_nn_mult_lin(left, right, product, submatrix_box, flop,&
       retain_sparsity, alpha_d, beta_d,&
       t_all, t_dgemm)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    INTEGER, DIMENSION(6), INTENT(in)        :: submatrix_box
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: alpha_d, beta_d
    REAL(KIND=dp), INTENT(INOUT)             :: t_all, t_dgemm

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_nn_mult_lin', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: beta_add = 1.0_dp, &
                                                beta_new = 0.0_dp

    CHARACTER                                :: a_type, b_type, use_a_type, &
                                                use_b_type
    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: left_data_z, right_data_z
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: left_data_c, right_data_c
    INTEGER :: a_blk, a_col, a_col_log, a_col_size, a_col_size_logical, &
      a_nze, a_row, a_row_l, a_row_log, a_row_size, a_row_size_logical, &
      b_blk, b_col, b_col_size, b_col_size_logical, b_nze, b_row, b_row_size, &
      b_row_size_logical, c_blk_p, c_blk_pt, c_col, c_col_size_logical, &
      c_nze, c_row, c_row_size_logical, clamp_k_len, clamp_k_offset, &
      datasize, first_k, first_row, i, ithread, j, last_k, last_row, lastblk, &
      lda, ldb, new_blk, new_row_blks, nlda, nldb, nthreads, old_blk, &
      old_row_begin_blk, old_size, prev_b_col, row_begin_blk, row_blk, rpt, &
      size_c_blk_ps, stat, submat_f_col, submat_f_row, submat_l_col, &
      submat_l_row
    INTEGER :: use_a_col_size, use_a_col_size_logical, use_a_offset, &
      use_a_row_size, use_a_row_size_logical, use_b_col_size, &
      use_b_col_size_logical, use_b_offset, use_b_row_size, &
      use_b_row_size_logical, use_c_col_size_logical, use_c_nze, &
      use_c_offset, use_c_row_size_logical, use_lda, use_ldb, use_nlda, &
      use_nldb
    INTEGER(KIND=int_4)                      :: offset
    INTEGER(KIND=int_8)                      :: coordinate, lflop
    INTEGER, ALLOCATABLE, DIMENSION(:) :: c_blk_ps, c_blk_ps_p, &
      left_col_blk_end, left_col_blk_start, left_row_blk_end, &
      left_row_blk_start, right_col_blk_end, right_col_blk_start
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: mt_row_processing, new_index
    INTEGER, DIMENSION(:), POINTER :: left_col_blk_size, left_local_rows, &
      left_row_blk_size, new_blk_p, new_col_i, new_row_i, &
      product_col_blk_size, right_col_blk_size, right_row_blk_size
    LOGICAL :: accumulated_data_valid, block_exists, found, keep_sparsity, &
      left_limit_rows, limit_k, right_limit_cols, start_fresh_extent, tr_a, &
      tr_b, tr_p, use_accumulated_data, use_tr_p
    REAL(KIND=dp)                            :: beta, epoch, mt_t_dgemm, &
                                                my_alpha, my_beta, use_alpha, &
                                                use_beta, which_beta
    REAL(KIND=dp), DIMENSION(:), POINTER     :: left_data_d, product_data_d, &
                                                right_data_d
    REAL(KIND=sp), DIMENSION(:), POINTER     :: left_data_r, right_data_r

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    my_alpha = 1.0_dp ; IF (PRESENT (alpha_d)) my_alpha = alpha_d
    my_beta = 0.0_dp ; IF (PRESENT (beta_d)) my_beta = beta_d
    keep_sparsity = .FALSE.
    IF (PRESENT (retain_sparsity)) THEN
       !WRITE(*,*)routineN//" Retaining sparsity!"
       keep_sparsity = .TRUE.
    END IF
    !CALL timeset(routineN, timing_handle)
    !!$OMP END MASTER
    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM ()
    nthreads = 1
    !$ nthreads = OMP_GET_NUM_THREADS ()
    !WRITE(*,*)'Hello! I am thread',ithread
!  write(*,*)routineN//" got data", my_beta
!  write(*,'(10(1X,F7.2))')dbcsr_get_data_p_d (&
!                     product%wms(ithread+1)%data_area)
    epoch = m_walltime()
    DBG "Doing NN matrix multiplication",left%name,' x ',right%name
    DBG "Pre-mult blk sizes, data", product%wms(ithread+1)%lastblk,&
         product%wms(ithread+1)%datasize
    left_row_blk_size => array_data (left%row_blk_size)
    left_col_blk_size => array_data (left%col_blk_size)
    right_row_blk_size => array_data (right%row_blk_size)
    right_col_blk_size => array_data (right%col_blk_size)
    product_col_blk_size => array_data (product%col_blk_size)
    !left_thread_dist => array_data (dbcsr_distribution_thread_dist (product%dist))
    left_local_rows => array_data (dbcsr_distribution_local_rows (left%dist))
    ! Submatrix handling
    CALL cp_assert (submatrix_box(3) .EQ. 0,&
         cp_warning_level, cp_unimplemented_error_nr, routineN,&
         "Limiting first column untested")
    CALL cp_assert (submatrix_box(5) .EQ. 0,&
         cp_warning_level, cp_unimplemented_error_nr, routineN,&
         "Limiting first k untested")
    left_limit_rows = MAXVAL (submatrix_box(1:2)) .GT. 0
    right_limit_cols = MAXVAL (submatrix_box(3:4)) .GT. 0
    limit_k = MAXVAL (submatrix_box(5:6)) .GT. 0
    IF (left_limit_rows) THEN
       ALLOCATE (left_row_blk_start (left%nblkrows_total),&
                 left_row_blk_end (left%nblkrows_total))
       CALL convert_sizes_to_offsets (left_row_blk_size,&
            left_row_blk_start, left_row_blk_end)
       submat_f_row = submatrix_box(1)
       submat_l_row = submatrix_box(2)
       !write(*,*)"limiting rows",submatrix_box(1:2)
    ENDIF
    IF (right_limit_cols) THEN
       ALLOCATE (right_col_blk_start (right%nblkcols_total),&
                 right_col_blk_end (right%nblkcols_total))
       CALL convert_sizes_to_offsets (right_col_blk_size,&
            right_col_blk_start, right_col_blk_end)
       submat_f_col = submatrix_box(3)
       submat_l_col = submatrix_box(4)
       !write(*,*)"limiting cols",submatrix_box(3:4)
    ENDIF
    IF (limit_k) THEN
       ALLOCATE (left_col_blk_start (left%nblkcols_total),&
                 left_col_blk_end (left%nblkcols_total))
       CALL convert_sizes_to_offsets (left_col_blk_size,&
            left_col_blk_start, left_col_blk_end)
       first_k = 1
       last_k = left%nfullcols_total
       IF (submatrix_box(5) .GT. 0) THEN
          first_k = submatrix_box(5)
       ENDIF
       IF (submatrix_box(6) .GT. 0) THEN
          last_k = submatrix_box(6)
       ENDIF
       !write(*,*)routineN//" Limiting k", submatrix_box(5:6)
    ENDIF
    !
    SELECT CASE (left%data_type)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (left%data_area, left_data_r)
       CALL dbcsr_get_data (right%data_area, right_data_r)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (left%data_area, left_data_d)
       CALL dbcsr_get_data (right%data_area, right_data_d)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (left%data_area, left_data_c)
       CALL dbcsr_get_data (right%data_area, right_data_c)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (left%data_area, left_data_z)
       CALL dbcsr_get_data (right%data_area, right_data_z)
    END SELECT
    c_blk_p = product%wms(ithread+1)%datasize + 1
    DBG 'A rows:',left%nblkrows_total
    old_size = product%wms(ithread+1)%lastblk
    mt_t_dgemm = 0.0_dp
    IF (left%nblks .GT. 0 .AND. right%nblks .GT. 0) THEN
    !IF (left%nfullrows_local*right%nfullcols_local .GT.&
    !     dbcsr_get_data_size (product%wms(ithread+1)%data_area)) THEN
    !   !!$OMP CRITICAL
    !   !CALL dbcsr_data_ensure_size(product%wms(ithread+1)%data_area,&
    !   !     MAX(product%wms(ithread+1)%datasize,&
    !   !     left%nfullrows_local*right%nfullcols_local))
    !   !!$OMP END CRITICAL
    !   !CALL ensure_array_size(product%wms(1)%data,&
    !   !     ub=MAX(product%wms(1)%datasize,&
    !   !     left%nfullrows_local*right%nfullcols_local),&
    !   !     error=error)
    !ENDIF
    ALLOCATE (mt_row_processing(rpnslots, left%nblkrows_total))
    ALLOCATE (new_index (4, left%nfullrows_local*right%nfullcols_local))
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
    accumulated_data_valid = .FALSE.
    new_blk = 0
    start_fresh_extent = .TRUE.
    !write(*,*)ithread,':',left_thread_dist
    rpt = CEILING (REAL(left%nblkrows_local) / REAL(nthreads))
    first_row = MIN (ithread*rpt+1, left%nblkrows_local)
    last_row = MIN ((ithread+1)*rpt, left%nblkrows_local)
    !write(*,*)ithread,'Doing rows', first_row, last_row
    !DO a_row_l = first_row, last_row
    DO a_row_l = 1, left%nblkrows_local
       a_row = left_local_rows (a_row_l)
       !write(*,*)ithread,'testing lrow',a_row_l,a_row
       !DO a_row = 1, left%nblkrows_total
       IF (left_limit_rows) THEN
          IF (left_row_blk_start (a_row) .GT. submat_l_row .OR.&
              left_row_blk_end (a_row) .LT. submat_f_row) THEN
             !write(*,*)"Skipping row", a_row
             CYCLE
          ELSE
             CALL cp_assert(&
                    MIN(submat_l_row, left_row_blk_end (a_row))&
                  - MAX(submat_f_row, left_row_blk_start (a_row))+1&
                  .EQ. left_row_blk_size (a_row),&
                  cp_fatal_level, cp_unimplemented_error_nr,&
                  routineN, "Submatrix must encompass complete block rows")
          ENDIF
       ENDIF
       a_row_log = a_row
       DBG '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
       DBG 'A row', a_row,' with #cols=',left%row_p(a_row+1)-left%row_p(a_row)
       a_row_size_logical = left_row_blk_size(a_row) ! stored
       c_row_size_logical = a_row_size_logical
       !c_blk_ps(:) = 0
       row_begin_blk = new_blk+1
       IF (old_blk .LE. old_size) THEN
          DO WHILE (product%wms(ithread+1)%row_i(old_blk) .LT. a_row)
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
          DO WHILE (product%wms(ithread+1)%row_i(old_blk) .LE. a_row) ! EQ
             DBG "Adding c_blk_ps column",&
                  product%wms(ithread+1)%col_i(old_blk),&
                  product%wms(ithread+1)%blk_p(old_blk)
             c_col = product%wms(ithread+1)%col_i(old_blk)
             found = .FALSE.
             !IF (keep_sparsity) THEN
             !   ! Check to see whether the block exists in the
             !   ! sparsity pattern.
             !   tr = .FALSE.
             !   CALL dbcsr_get_block_index (product, a_row, a_col,&
             !        b_row, b_col, tr, found, i)
             !ENDIF
             !use_block: IF (.NOT. keep_sparsity) THEN
                c_blk_ps(c_col)&
                     = product%wms(ithread+1)%blk_p(old_blk)
                size_c_blk_ps = size_c_blk_ps + 1
                c_blk_ps_p(size_c_blk_ps) = product%wms(ithread+1)%col_i(old_blk)
             !ELSEIF (dbg) THEN
             !   WRITE(*,*)routineN//" Skipping sparse block", a_row, c_col
             !ENDIF use_block
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
          IF (limit_k) THEN
             IF (left_col_blk_start (a_col) .GT. last_k .OR.&
                 left_col_blk_end (a_col) .LT. first_k) THEN
                !write(*,'(A,I5,"(",I5,"-",I5,") not in ",I5,"-",I5)')&
                !     routineN//"Skipping a col", a_col,&
                !     left_col_blk_start (a_col), left_col_blk_end (a_col),&
                !     first_k, last_k
                CYCLE
             ELSE
                i = MAX (first_k, left_col_blk_start (a_col))
                j = MIN (last_k, left_col_blk_end (a_col))
                clamp_k_offset = i - left_col_blk_start (a_col)
                clamp_k_len = j - i + 1
                !write(*,'(A,I5,"(",I5,"-",I5,") not in ",I5,"-",I5)')&
                !     routineN//"Not skipping a col", a_col,&
                !     left_col_blk_start (a_col), left_col_blk_end (a_col),&
                !     first_k, last_k
                !write(*,*)routineN//" Clamping", clamp_k_offset, clamp_k_len
             ENDIF
          ENDIF
          a_row_log = a_row
          a_col_log = a_col
          a_col_size_logical = left_col_blk_size(a_col_log) !stored
          a_row_size = a_row_size_logical
          a_col_size = a_col_size_logical
          lda = a_row_size
          nlda = a_col_size
          tr_a = .FALSE.
          a_type = 'N'
          DBG '---------------'
          DBG 'A col', a_col, '=>',a_col_log
          b_row = a_col_log
          c_row = a_row_log
          a_nze = a_row_size_logical * a_col_size_logical
          !a_data => left%data(ABS(left%blk_p(a_blk))&
          !     :ABS(left%blk_p(a_blk))+a_nze-1)
          b_row_size_logical = a_col_size_logical
          c_row_size_logical = a_row_size_logical
          ! Handle transposed blocks.
          IF (left%blk_p(a_blk) .LT. 0) THEN
             CALL swap (lda, nlda)
             a_type = 'T'
             CALL swap (a_row_size, a_col_size)
          ENDIF
          ldb = b_row_size_logical
          start_fresh_extent = .TRUE.
          DO b_blk = right%row_p(b_row)+1, right%row_p(b_row+1)
             b_col = right%col_i(b_blk)
             IF (right_limit_cols) THEN
                IF (right_col_blk_start (b_col) .GT. submat_l_col .OR.&
                    right_col_blk_end (b_col) .LT. submat_f_col) THEN
                   !WRITE(*,*)"Skipping column", b_col
                   CYCLE
                ENDIF
             ENDIF
             c_col = b_col
             DBG 'B col', b_col
             b_col_size_logical = right_col_blk_size(b_col)
             IF (right_limit_cols) THEN
                c_col_size_logical = product_col_blk_size(c_col)
                !write(*,*)routineN//" c_col_size_log", c_col_size_logical, b_col_size_logical
             ELSE
                c_col_size_logical = b_col_size_logical
             ENDIF
             b_col_size = b_col_size_logical
             b_nze = b_row_size_logical * b_col_size_logical
             !b_data => right%data(&
             !     ABS(right%blk_p(b_blk)):ABS(right%blk_p(b_blk))+b_nze-1)
             tr_b = right%blk_p(b_blk) .LT. 0
             b_row_size = b_row_size_logical
             IF (right%blk_p(b_blk) .LT. 0) THEN
                b_type = 'T'
                b_col_size = b_row_size_logical
                ldb = b_col_size_logical
                nldb = b_row_size_logical
                b_row_size = b_col_size_logical
                !b_col_size = b_row_size_logical
             ELSEIF (right%blk_p(b_blk) .GT. 0) THEN
                b_type = 'N'
                ldb = b_row_size_logical
                nldb = b_col_size_logical
             ELSE
                b_type = '0'
                ldb = 0
                nldb = 0
                !CYCLE
             ENDIF
             ! This must be adjusted afterwards so that if tr_b, then 
             ! b_row_size retains the full value!
             IF (right_limit_cols) THEN
                b_col_size_logical =&
                     MIN (right_col_blk_end (b_col), submat_l_col) &
                     - MAX (right_col_blk_start (b_col), submat_f_col)+1
                !write(*,*)routineN//" b_col_size_log", b_col_size_logical
             ENDIF
             c_nze = c_row_size_logical * c_col_size_logical
             offset = ABS(c_blk_ps(c_col))
             tr_p = c_blk_ps(c_col) .LT. 0
             IF (tr_p) use_accumulated_data = .TRUE.
             block_exists = c_blk_ps(c_col) .NE. 0
             DBG "Found old?", block_exists, c_blk_ps(c_col), offset
             ! Try to find the current coordinate in the current row index
             sparsity_enforcement: IF &
                  (keep_sparsity .AND. .NOT. block_exists .OR. b_type.eq.'0') THEN
                use_accumulated_data = .TRUE.
             ELSE
                IF (.NOT. block_exists) THEN
                   new_blk = new_blk+1
                   new_index(1, new_blk) = c_row
                   new_index(2, new_blk) = c_col
                   new_index(3, new_blk) = c_blk_pt
                   new_index(4, new_blk) = ithread+1
                   new_row_blks = new_row_blks + 1
                   DBG 'Adding coordinate', c_row, c_col, c_blk_p
                   !write(*,*)routineN//' Adding coordinate', c_row, c_col, c_blk_p
                   offset = c_blk_pt
                   c_blk_pt = c_blk_pt + c_nze
                   !WRITE(*,'("(",I4,"): ",A,1X,I7)')ithread,"GOT offset", offset
                   datasize = datasize + c_nze
                   lastblk = lastblk + 1
                   c_blk_ps (c_col) = offset
                   size_c_blk_ps = size_c_blk_ps + 1
                   c_blk_ps_p(size_c_blk_ps) = c_col
                   which_beta = beta_new
                   tr_p = .FALSE.
                ELSE
                   DBG 'Using existing coordinate at', offset, tr_p
                   !write(*,*)routineN//' Using existing coordinate at', offset,&
                   !     pointer_view (dbcsr_get_data_p_d (&
                   !     product%wms(ithread+1)%data_area),&
                   !     offset, offset+c_nze-1)
                   which_beta = beta_add
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
                   use_a_row_size = a_row_size
                   use_a_col_size = a_col_size
                   use_a_row_size_logical = a_row_size_logical
                   use_a_col_size_logical = a_col_size_logical
                   use_b_row_size = b_row_size
                   use_b_col_size = b_col_size
                   use_b_row_size_logical = b_row_size_logical
                   use_b_col_size_logical = b_col_size_logical
                   IF (right_limit_cols) THEN
                      use_b_col_size_logical =&
                           MIN (right_col_blk_end (b_col), submat_l_col) &
                           - MAX (right_col_blk_start (b_col), submat_f_col)+1
                      use_accumulated_data = .TRUE.
                   ENDIF
                   IF (.NOT. limit_k) THEN
                      use_a_col_size_logical = a_col_size_logical
                   ELSE
                      use_a_col_size_logical = clamp_k_len
                   ENDIF
                   IF (use_a_type .EQ. 'T' .AND. left%sym_negation) THEN
                      use_alpha = -my_alpha
                   ELSE
                      use_alpha = my_alpha
                   ENDIF
                   IF (use_b_type .EQ. 'T' .AND. right%sym_negation) THEN
                      use_alpha = -use_alpha
                   ENDIF
                   IF (tr_p .AND. product%sym_negation) THEN
                      use_alpha = -use_alpha
                   ENDIF
                   use_a_offset = ABS (left%blk_p(a_blk))
                   use_lda = lda
                   use_nlda = nlda
                   use_b_offset = ABS (right%blk_p(b_blk))
                   use_ldb = ldb
                   use_nldb = nldb
                   use_beta = which_beta
                   use_c_offset = offset
                   use_c_nze = use_c_row_size_logical * use_c_col_size_logical
                   use_c_col_size_logical = c_col_size_logical
                   use_c_row_size_logical = c_row_size_logical
                   prev_b_col = b_col
                   use_tr_p = tr_p
                   start_fresh_extent = .FALSE.
                   accumulated_data_valid = .TRUE.
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
                   use_b_col_size_logical = use_b_col_size_logical + b_col_size_logical
                   use_c_col_size_logical = use_c_col_size_logical + c_col_size_logical
                ENDIF new_extent
                end_of_row: IF (b_blk .EQ. right%row_p(b_row+1)) THEN
                   use_accumulated_data = .TRUE.
                ENDIF end_of_row
             ENDIF sparsity_enforcement
!vw this doesnt work for right sparse matrix !
!vw   flush_dgemm: IF (use_accumulated_data) THEN
!vw reset to true for the moment !
             flush_dgemm: IF (accumulated_data_valid) THEN
!vw
!
                mt_t_dgemm = mt_t_dgemm - (m_walltime() - epoch)
                IF (dbg) WRITE(*,&
                     "(1X,A,A1,1X,A1,A,I6,A,F4.1)")&
                     routineN//' DGEMM ',use_a_type,use_b_type,&
                     ' at',use_c_offset,' beta=', use_beta
!WRITE(*,*) 'sum A',SUM(left_data_d(use_a_offset : use_a_offset + use_a_row_size_logical*use_a_col_size_logical-1))
!WRITE(*,*) 'sum B',SUM(right_data_d(use_b_offset : use_b_offset + use_a_col_size_logical * use_b_col_size-1))
!WRITE(*,*) 'sum C',SUM(pointer_view (dbcsr_get_data_p_d (&
!                            product%wms(ithread+1)%data_area),&
!                            use_c_offset, use_c_offset+c_nze-1) )
                product_data_d => pointer_view (dbcsr_get_data_p_d (&
                     product%wms(ithread+1)%data_area),&
                     use_c_offset, use_c_offset+c_nze-1)
                ! Avoid nagging compilers
                IF (use_beta .EQ. 0.0_dp) THEN
                   product_data_d(1&
                        :use_c_row_size_logical*use_c_col_size_logical) = 0.0_dp
                ENDIF
                IF (.NOT. use_tr_p) THEN
                   IF (dbg) &
                        WRITE(*,'("NRC ",A,A,2X,3(I5,1X),F5.1,1X,2(I9,1X,I5,1X),F5.1,1X,I9,1X,I5)')&
                        use_a_type, use_b_type,&
                        use_c_row_size_logical,&
                        use_c_col_size_logical,&
                        use_a_col_size_logical,&
                        use_alpha,&
                        use_a_offset, use_a_row_size,&
                        use_b_offset, use_b_row_size, use_beta,&
                        use_c_offset,&
                        use_c_row_size_logical
                   use_c_col_size_logical = MIN(&
                        use_b_col_size_logical,&
                        use_c_col_size_logical)
                   IF(use_c_row_size_logical.GT.0.AND.use_c_col_size_logical.GT.0&
                        .AND.use_a_col_size_logical.GT.0) THEN
                   CALL DGEMM(use_a_type, use_b_type,&
                        use_c_row_size_logical, use_c_col_size_logical,&
                        use_a_col_size_logical,&
                        use_alpha,&
                        left_data_d(use_a_offset), use_a_row_size,&
                        right_data_d(use_b_offset), use_b_row_size, use_beta,&
                        product_data_d,&
                        use_c_row_size_logical)
                   ENDIF
                ELSE
                   !WRITE(*,*)'Flipping'
                   CALL flip_type (use_a_type)
                   CALL flip_type (use_b_type)
                   use_c_col_size_logical = MIN(&
                        use_b_col_size_logical,&
                        use_c_col_size_logical)
                   IF (dbg) &
                        WRITE(*,'("TRC ",A,A,2X,3(I5,1X),F5.1,1X,2(I9,1X,I5,1X),F5.1,1X,I9,1X,I5)')&
                        use_b_type, use_a_type,&
                        use_c_col_size_logical,&
                        use_c_row_size_logical,&
                        use_a_col_size_logical,&
                        use_alpha,&
                        use_b_offset, use_b_row_size,&
                        use_a_offset, use_a_row_size, use_beta,&
                        use_c_offset,&
                        use_c_col_size_logical
                   IF(use_c_col_size_logical.GT.0.AND.use_c_row_size_logical.GT.0&
                        .AND.use_a_col_size_logical.GT.0) THEN
                   CALL DGEMM(use_b_type, use_a_type,&
                        use_c_col_size_logical, use_c_row_size_logical,&
                        use_a_col_size_logical,&
                        !use_a_row_size_logical, use_b_col_size,&
                        !use_a_col_size_logical,&
                        use_alpha,&
                        !use_alpha,&
                        right_data_d(use_b_offset), use_b_row_size,&
                        !left_data_d(use_a_offset), use_lda,&
                        left_data_d(use_a_offset), use_a_row_size, use_beta,&
                        !right_data_d(use_b_offset), use_ldb, use_beta,&
                        product_data_d,&
                        use_c_col_size_logical)
                        !use_a_row_size_logical)
                   ENDIF
                ENDIF
                mt_t_dgemm = mt_t_dgemm + (m_walltime() - epoch)
                lflop = lflop + INT(use_a_row_size_logical, int_8)&
                     * INT(use_b_col_size, int_8)&
                     * INT(use_a_col_size_logical, int_8) * INT(2, int_8)
             ENDIF flush_dgemm
             reset_accumulation: IF (.TRUE.) THEN
                accumulated_data_valid = .FALSE.
                use_accumulated_data = .FALSE.
                use_b_col_size = 0
                start_fresh_extent = .TRUE.
             ENDIF reset_accumulation
          ENDDO ! b
          row_reset_accumulation: IF (.TRUE.) THEN
             accumulated_data_valid = .FALSE.
             use_accumulated_data = .FALSE.
             use_b_col_size = 0
             start_fresh_extent = .TRUE.
          ENDIF row_reset_accumulation
       ENDDO ! a_col
       mt_row_processing(rpslot_owner, a_row) = ithread
       mt_row_processing(rpslot_addblks, a_row) = new_row_blks
    ENDDO ! a_row Final flush This is a copy of the flush appearing in
    ! the main loop. They should be factored out.
    IF (accumulated_data_valid) THEN
       mt_t_dgemm = mt_t_dgemm - (m_walltime() - epoch)
       IF (dbg) WRITE(*,&
            "(1X,A,A1,1X,A1,A,I6,A,F4.1)")&
            routineN//' DGEMM ',use_a_type,use_b_type,&
            ' at',use_c_offset,' beta=', use_beta
       product_data_d => pointer_view (dbcsr_get_data_p_d (&
            product%wms(ithread+1)%data_area),&
            use_c_offset, use_c_offset+c_nze-1)
       ! Avoid nagging compilers
       IF (use_beta .EQ. 0.0_dp) THEN
          product_data_d(1&
               :use_c_row_size_logical*use_c_col_size_logical) = 0.0_dp
       ENDIF
       IF (.NOT. use_tr_p) THEN
          IF (dbg) &
               WRITE(*,'("NRC ",A,A,2X,3(I5,1X),F5.1,1X,2(I9,1X,I5,1X),F5.1,1X,I9,1X,I5)')&
               use_a_type, use_b_type,&
               use_c_row_size_logical,&
               use_c_col_size_logical,&
               use_a_col_size_logical,&
               use_alpha,&
               use_a_offset, use_a_row_size,&
               use_b_offset, use_b_row_size, use_beta,&
               use_c_offset,&
               use_c_row_size_logical
          use_c_col_size_logical = MIN(&
               use_b_col_size_logical,&
               use_c_col_size_logical)
          CALL DGEMM(use_a_type, use_b_type,&
               use_c_row_size_logical, use_c_col_size_logical,&
               use_a_col_size_logical,&
               use_alpha,&
               left_data_d(use_a_offset), use_a_row_size,&
               right_data_d(use_b_offset), use_b_row_size, use_beta,&
               product_data_d,&
               use_c_row_size_logical)
       ELSE
          CALL flip_type (use_a_type)
          CALL flip_type (use_b_type)
          IF (dbg) &
               WRITE(*,'("TRC ",A,A,2X,3(I5,1X),F5.1,1X,2(I9,1X,I5,1X),F5.1,1X,I9,1X,I5)')&
               use_b_type, use_a_type,&
               use_c_col_size_logical,&
               use_c_row_size_logical,&
               use_a_col_size_logical,&
               use_alpha,&
               use_b_offset, use_b_row_size,&
               use_a_offset, use_a_row_size, use_beta,&
               use_c_offset,&
               use_c_col_size_logical
          use_c_col_size_logical = MIN(&
               use_b_col_size_logical,&
               use_c_col_size_logical)
          CALL DGEMM(use_b_type, use_a_type,&
               use_c_col_size_logical, use_c_row_size_logical,&
               use_a_col_size_logical,&
               !use_a_row_size_logical, use_b_col_size,&
               !use_a_col_size_logical,&
               use_alpha,&
               !use_alpha,&
               right_data_d(use_b_offset), use_b_row_size,&
               !left_data_d(use_a_offset), use_lda,&
               left_data_d(use_a_offset), use_a_row_size, use_beta,&
               !right_data_d(use_b_offset), use_ldb, use_beta,&
               product_data_d,&
               use_c_col_size_logical)
          !use_a_row_size_logical)
       ENDIF
       mt_t_dgemm = mt_t_dgemm + (m_walltime() - epoch)
       lflop = lflop + INT(use_a_row_size_logical, int_8)&
            * INT(use_b_col_size, int_8)&
            * INT(use_a_col_size_logical, int_8) * INT(2, int_8)
    ENDIF
    !
    DEALLOCATE (c_blk_ps, c_blk_ps_p)
    ! Now merge the new indices with the old index.
    !!CALL accumulate_row_offsets (mt_row_processing, exchange,&
    !!     SIZE(mt_row_processing,2))
    ! The upper replaces this commented block.
    mt_row_processing(rpslot_totaloffset, 1) = 1
    mt_row_processing(rpslot_oldoffset, 1) = 1
    DO a_row = 2, left%nblkrows_total
       mt_row_processing(rpslot_totaloffset, a_row)&
            = mt_row_processing(rpslot_totaloffset, a_row-1)&
            + mt_row_processing(rpslot_addblks, a_row-1)&
            + mt_row_processing(rpslot_oldblks, a_row-1)
       mt_row_processing(rpslot_oldoffset, a_row)&
            = mt_row_processing(rpslot_oldoffset, a_row-1)&
            + mt_row_processing(rpslot_oldblks, a_row-1)
    ENDDO

    product%wms(ithread+1)%lastblk&
         = product%wms(ithread+1)%lastblk + lastblk
    product%wms(ithread+1)%datasize&
         = product%wms(ithread+1)%datasize + datasize
    !WRITE(*,*)"New total # blocks:", product%w%lastblk
    !WRITE(*,*)"New all # blocks:  ", product%w%lastblks
    !WRITE(*,*)"New total size:", product%w%datasize
    !WRITE(*,*)"New all sizes: ", product%w%datasizes
    i = old_size + new_blk
    ALLOCATE (new_row_i(i), new_col_i(i), new_blk_p(i))
    CALL merge_arrays (new_row_i, new_col_i, new_blk_p, i,&
         product%wms(ithread+1)%row_i, product%wms(ithread+1)%col_i,&
         product%wms(ithread+1)%blk_p,&
         old_size,&
         new_index, new_blk)
    ! Replace the index.
    CALL pointer_replace (product%wms(ithread+1)%row_i, new_row_i)
    CALL pointer_replace (product%wms(ithread+1)%col_i, new_col_i)
    CALL pointer_replace (product%wms(ithread+1)%blk_p, new_blk_p)
    DEALLOCATE (new_index)
    ELSE
       DBG 'One of the matrices is empty!'
    ENDIF
    DBG "Post-mult blk sizes, data", product%wms(ithread+1)%lastblk,&
         product%wms(ithread+1)%datasize
    IF (PRESENT (flop)) flop = lflop
    t_all = t_all + (m_walltime() - epoch)
    t_dgemm = t_dgemm + mt_t_dgemm
    !WRITE(*,*)'Goodbye! I am thread',ithread
    !!$OMP MASTER
    !CALL timestop(timing_handle)
    !!$OMP END MASTER
    !!$BARRIER
  END SUBROUTINE dbcsr_nn_mult_lin


  elemental SUBROUTINE flip_type (t)
    CHARACTER, INTENT(INOUT)                 :: t

    SELECT CASE (t)
    CASE ('N')
       t = 'T'
    CASE ('T')
       t = 'N'
    CASE DEFAULT
       t = '@'
    END SELECT
  END SUBROUTINE flip_type


! *****************************************************************************
!> \brief Merges two indices
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

    INTEGER                                  :: add_blk, i, merge_from_whom, &
                                                new_blk, old_blk
    LOGICAL                                  :: multidata

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    multidata = PRESENT (old_blk_d) .AND. PRESENT (new_blk_d)
    CALL bcsr_assert(2, old_size+add_size .EQ. new_size, "merge_arrays",&
         "Mismatch of new and old size")
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
    IF (dbg) THEN
       WRITE (*,*) " New array"
       DO i = 1, new_size
          WRITE(*,'(4(2X,I7))')new_row_i(i),new_col_i(i),new_blk_p(i)
       ENDDO
    ENDIF
  END SUBROUTINE merge_arrays


!! *****************************************************************************
!!> \brief Merges many sorted indices
!!> \param[in,out] new_row_i, new_col_i, new_blk_p      merge index
!!> \param[in,out] old_row_i, old_col_i, old_blk_p      current index
!!> \param[in] add_ip                    index to add into the current index
!!> \param[in] old_size                  size of current index
!!> \param[in] add_sizes                 size of index to add into the current
!!>                                      index
!!> \param[in] row_metadata              who processed a row, how many elements
!!>                                      were added, and how many were in the
!!>                                      old row
!!> \param[in] new_blk_d, old_blk_d      (optional) new and old location of data
!!>                                      blocks
!! *****************************************************************************
!  SUBROUTINE merge_many_arrays (&
!       new_row_i, new_col_i, new_blk_p,&
!       old_row_i, old_col_i, old_blk_p, old_size,&
!       add_ip, add_sizes, row_metadata,&
!       new_blk_d, old_blk_d)
!    INTEGER, DIMENSION(:), POINTER           :: new_row_i, new_col_i, &
!                                                new_blk_p, old_row_i, &
!                                                old_col_i, old_blk_p
!    INTEGER, INTENT(IN)                      :: old_size
!    INTEGER, DIMENSION(:, :, 0:), INTENT(IN) :: add_ip
!    INTEGER, DIMENSION(0:), INTENT(IN)       :: add_sizes
!    INTEGER, DIMENSION(:, :), INTENT(IN)     :: row_metadata
!    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: new_blk_d, old_blk_d
!
!    INTEGER                                  :: a1, a2, ablks, add_size, &
!                                                ithread, n1, n2, new_size, &
!                                                nrows, nthreads, o1, o2, &
!                                                oblks, row
!    INTEGER, DIMENSION(:), POINTER           :: tmp_p
!    LOGICAL                                  :: multidata
!
!!
!
!    multidata = PRESENT (old_blk_d) .AND. PRESENT (new_blk_d)
!    ithread = 0
!    nthreads = 1
!    !$ ithread = OMP_GET_THREAD_NUM ()
!    !$ nthreads = OMP_GET_NUM_THREADS ()
!    !IF (ithread .EQ. 0) THEN
!    !   WRITE (*,*) " MANY Old array", old_size
!    !   DO n1 = 1, old_size
!    !      WRITE(*,'(I7,2X,I7,2X,I7)')old_row_i(n1),old_col_i(n1),old_blk_p(n1)
!    !   ENDDO
!    !   DO n1 = 0, nthreads-1
!    !      WRITE (*,*) " MANY Add array from", n1, ":", add_sizes(n1)
!    !      DO n2 = 1, add_sizes(n1)
!    !         WRITE(*,'(4(2X,I7))')add_ip (1:4, n2, n1)
!    !      ENDDO
!    !   ENDDO
!    !   WRITE(*,*)"Metadata:"
!    !   DO n1 = 1, SIZE (row_metadata, 2)
!    !      WRITE(*,'(6(1X,I5))') row_metadata(:, n1)
!    !   ENDDO
!    !ENDIF
!    !!$OMP BARRIER
!    nrows = SIZE(row_metadata, 2)
!    add_size = SUM (add_sizes)
!    new_size = old_size+add_size
!    ! There's nothing to do if nothing was added
!    IF (add_size .GT. 0) THEN
!       ! Process my rows.
!       DO row = 1, nrows
!          IF (row_metadata(rpslot_owner, row) .EQ. ithread) THEN
!             oblks = row_metadata(rpslot_oldblks, row)
!             o1 = row_metadata(rpslot_oldoffset, row)
!             o2 = o1 + oblks - 1
!             ablks = row_metadata(rpslot_addblks, row)
!             a1 = row_metadata(rpslot_addoffset, row)
!             a2 = a1 + ablks - 1
!             n1 = row_metadata(rpslot_totaloffset, row)
!             n2 = n1 + oblks + ablks - 1
!             IF (n2 .GE. n1) THEN
!                CALL merge_arrays(new_row_i(n1:n2), new_col_i(n1:n2),&
!                     new_blk_p(n1:n2), oblks + ablks,&
!                     old_row_i(o1:o2), old_col_i(o1:o2),&
!                     old_blk_p(o1:o2), oblks,&
!                     add_ip(:,a1:a2,ithread), ablks,&
!                     new_blk_d(n1:n2), old_blk_d(o1:o2))
!             ENDIF
!          END IF
!       ENDDO
!       ! Now set the pointers to the arrays with the new values.
!       IF (ithread .EQ. 0) THEN
!          tmp_p => old_row_i
!          old_row_i => new_row_i
!          DEALLOCATE (tmp_p)
!          tmp_p => old_col_i
!          old_col_i => new_col_i
!          DEALLOCATE (tmp_p)
!          tmp_p => old_blk_p
!          old_blk_p => new_blk_p
!          DEALLOCATE (tmp_p)
!          IF (multidata) THEN
!             tmp_p => old_blk_d
!             old_blk_d => new_blk_d
!             DEALLOCATE (tmp_p)
!          ENDIF
!       ENDIF
!    ELSE
!       IF (ithread .EQ. 0) THEN
!          DEALLOCATE (new_row_i, new_col_i, new_blk_p)
!          IF (multidata) DEALLOCATE (new_blk_d)
!       ENDIF
!    ENDIF
!    !!!$OMP BARRIER
!    !IF (ithread .EQ. 0) THEN
!    !   WRITE (*,*) " MANY New array"
!    !   DO i = 1, new_size
!    !      WRITE(*,'(4(2X,I7))')old_row_i(i),old_col_i(i),&
!    !           old_blk_p(i),old_blk_d(i)
!    !   ENDDO
!    !ENDIF
!  END SUBROUTINE merge_many_arrays


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
    ithread = 0
    nthreads = 1
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
    TYPE(dbcsr_obj), INTENT(IN)              :: addee
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER, INTENT(OUT), OPTIONAL           :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_wm_add', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: alpha = 1.0_dp, &
                                                beta_add = 1.0_dp, &
                                                beta_new = 0.0_dp

    CHARACTER                                :: addee_type, sum_type
    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: addee_data_z
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: addee_data_c
    INTEGER :: addee_blk, addee_col, addee_col_size, addee_nze, addee_row, &
      addee_row_size, ldaddee, ldsum, lflop, sum_blk, sum_blk_p, &
      sum_col_size, sum_row_size, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: addee_col_blk_size, &
                                                addee_row_blk_size
    REAL(KIND=dp)                            :: tstart, tstop
    REAL(KIND=dp), DIMENSION(:), POINTER     :: addee_data, addee_data_d, &
                                                sum_data
    REAL(KIND=sp), DIMENSION(:), POINTER     :: addee_data_r

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    IF (dbg) tstart = m_walltime ()
    DBG "Doing WM matrix sum",' += ',addee%m%name
    DBG "Pre-mult blk sizes, data", wm%lastblk, wm%datasize
    DBG "Addee blk sizes, data", addee%m%nblks, addee%m%nze
    addee_row_blk_size => array_data (addee%m%row_blk_size)
    addee_col_blk_size => array_data (addee%m%col_blk_size)
    SELECT CASE (addee%m%data_type)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (addee%m%data_area, addee_data_r)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (addee%m%data_area, addee_data_d)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (addee%m%data_area, addee_data_c)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (addee%m%data_area, addee_data_z)
    END SELECT
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
    IF (wm%lastblk .GT. 0 .OR. addee%m%nblks .GT. 0) THEN
    DO addee_row = 1, addee%m%nblkrows_total
       DBG 'Addee row', addee_row,' with #cols=',&
            addee%m%row_p(addee_row+1)-addee%m%row_p(addee_row)
       ! Find the blocks present in addee_row.
       addee_row_size = addee_row_blk_size(addee_row)
       CALL dbcsr_data_ensure_size(wm%data_area,&
            wm%datasize+addee_row_size*addee%m%nfullcols_total &
            * MAX(1, (addee%m%nblkrows_total-addee_row)/4))
       !CALL ensure_array_size(wm%data,&
       !     ub=wm%datasize+addee_row_size*addee%m%nfullcols_total &
       !     * MAX(1, (addee%m%nblkrows_total-addee_row)/4), error=error)
       DO addee_blk = addee%m%row_p(addee_row)+1, addee%m%row_p(addee_row+1)
          addee_col = addee%m%col_i(addee_blk)
          DBG 'Addee col', addee_col
          addee_col_size = addee_col_blk_size(addee_col)
          addee_nze = addee_row_size * addee_col_size
          addee_data_d => addee_data_d(ABS(addee%m%blk_p(addee_blk))&
               :ABS(addee%m%blk_p(addee_blk))+addee_nze-1)
          !addee_data => addee%m%data(ABS(addee%m%blk_p(addee_blk))&
          !                   :ABS(addee%m%blk_p(addee_blk))+addee_nze-1)
          IF (addee%m%blk_p(addee_blk) .LT. 0) THEN
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
             sum_data => pointer_view (dbcsr_get_data_p (wm%data_area, 0.0_dp),&
                  ABS(sum_blk_p), ABS(sum_blk_p)+addee_nze-1)
             !sum_data => wm%data(&
             !     ABS(sum_blk_p):ABS(sum_blk_p)+addee_nze-1)
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
             sum_data => pointer_view (dbcsr_get_data_p (wm%data_area, 0.0_dp),&
                  ABS(sum_blk_p), ABS(sum_blk_p)+addee_nze-1)
             !sum_data => wm%data(&
             !     ABS(sum_blk_p):ABS(sum_blk_p)+addee_nze-1)
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
!> \param[in] side              (optional) ?
!> \param[in] last_column       (optional) only scale up to this full column
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_scale(matrix_a, alpha_scalar, alpha_vector, alpha_matrix,&
       last_column, side, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(dp), INTENT(IN), OPTIONAL           :: alpha_scalar
    REAL(dp), DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: alpha_vector, alpha_matrix
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: data_z
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: data_c
    INTEGER                                  :: a_blk, a_col, a_col_size, &
                                                a_nze, a_offset, a_row, &
                                                a_row_size, i, last_col, &
                                                timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: col_end, col_start, m_offset
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: do_scale, limit_col, right, tr
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_d
    REAL(KIND=sp), DIMENSION(:), POINTER     :: data_r
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    !
    !
    row_blk_size => array_data (matrix_a%m%row_blk_size)
    col_blk_size => array_data (matrix_a%m%col_blk_size)
    ! last_column => alpha_scalar
    CALL cp_assert (.NOT. PRESENT (last_column) .OR. PRESENT (alpha_scalar),&
         cp_fatal_level, cp_wrong_args_error, routineN,&
         "Column limiting supported only with scalar scaling")
    last_col = 0
    IF (PRESENT (last_column)) THEN
       CALL cp_assert (last_column .GE. 0 .AND. last_column .LE.&
            dbcsr_nfullcols_total (matrix_a), cp_fatal_level,&
            cp_wrong_args_error, routineN,&
            "Specified last column is out of bounds.")
       last_col = last_column
       ALLOCATE (col_start(SIZE (col_blk_size)))
       ALLOCATE (col_end(SIZE (col_blk_size)))
       CALL convert_sizes_to_offsets(col_blk_size, col_start, col_end)
    ENDIF
    SELECT CASE (matrix_a%m%data_type)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (matrix_a%m%data_area, data_r)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (matrix_a%m%data_area, data_d)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (matrix_a%m%data_area, data_c)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (matrix_a%m%data_area, data_z)
    END SELECT
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
            &   m_offset(matrix_a%m%nblkcols_total))
       buff = 0.0_dp
       m_offset(1) = 1
       IF(right) THEN
          DO a_col = 2,matrix_a%m%nblkcols_total
             a_col_size = col_blk_size(a_col-1)
             m_offset(a_col) = m_offset(a_col-1) + a_col_size**2
          ENDDO
       ELSE
          DO a_row = 2,matrix_a%m%nblkrows_total
             a_row_size = row_blk_size(a_row-1)
             m_offset(a_row) = m_offset(a_row-1) + a_row_size**2
          ENDDO          
       ENDIF
    ENDIF
    !
    ! quick return if possible
    do_scale = .TRUE.
    IF(PRESENT(alpha_scalar)) THEN
       IF(alpha_scalar.EQ.1.0_dp) do_scale = .FALSE.
    ENDIF
    !
    ! let's go
    IF(do_scale) THEN
       CALL dbcsr_iterator_start(iter, matrix_a)
       iterations: DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block(iter, a_row, a_col, a_blk, tr, a_offset)
          a_row_size = row_blk_size(a_row)
          a_col_size = col_blk_size(a_col)
          a_nze = a_row_size * a_col_size
          !
          ! let's scale
          IF(PRESENT(alpha_scalar)) THEN
             !a_data = a_data * alpha_scalar
             SELECT CASE (dbcsr_get_data_type(matrix_a))
                !CASE (dbcsr_type_real_4)
             CASE (dbcsr_type_real_8)
                limit_col = last_col .GT. 0
                IF (limit_col) limit_col = col_end(a_col).gt.last_col
                IF (limit_col) THEN
                   IF (col_start(a_col) .LE. last_col) THEN
                      i = MIN(last_col, col_end(a_col)) - col_start(a_col) + 1
                      CALL dscal(a_row_size*i, alpha_scalar,&
                           data_d(a_offset), 1)
                   ELSE
                      ! Nothing to do, the block is beyond the boundry.
                   ENDIF
                ELSE
                   CALL dscal(a_nze,alpha_scalar,data_d(a_offset),1)
                ENDIF
                !CASE (dbcsr_type_complex_4)
                !CASE (dbcsr_type_complex_8)
             CASE default
                CALL cp_unimplemented_error (routineN,&
                     "unknown precision")
             END SELECT
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
                SELECT CASE (dbcsr_get_data_type(matrix_a))
                   !CASE (dbcsr_type_real_4)
                CASE (dbcsr_type_real_8)
                   CALL dgemm('N','N',a_row_size,a_col_size,a_col_size,&
                        &     1.0_dp,data_d(a_offset),a_row_size,&
                        &            alpha_matrix(m_offset(a_col)),a_col_size,&
                        &     0.0_dp,buff,a_row_size)
                   CALL dcopy(a_nze,buff,1,data_d(a_offset),1)
                   !CASE (dbcsr_type_complex_4)
                   !CASE (dbcsr_type_complex_8)
                CASE default
                   CALL cp_unimplemented_error (routineN,&
                        "Only real double precision")
                END SELECT
             ELSE
                !A = alpha * A
                SELECT CASE (dbcsr_get_data_type(matrix_a))
                   !CASE (dbcsr_type_real_4)
                CASE (dbcsr_type_real_8)
                   CALL dgemm('N','N',a_row_size,a_col_size,a_row_size,&
                        &     1.0_dp,alpha_matrix(m_offset(a_row)),a_row_size,&
                        &            data_d(a_offset),a_row_size,&
                        &     0.0_dp,buff,a_row_size)
                   CALL dcopy(a_nze,buff,1,data_d(a_offset),1)
                   !CASE (dbcsr_type_complex_4)
                   !CASE (dbcsr_type_complex_8)
                CASE default
                   CALL cp_unimplemented_error (routineN,&
                        "Only real double precision")
                END SELECT
             ENDIF
          ELSE
             CALL stop_program(routineN,"wrong logic")
          ENDIF
       ENDDO iterations
       CALL dbcsr_iterator_stop(iter)

    ENDIF

    IF(PRESENT(alpha_matrix)) THEN
       DEALLOCATE(buff,m_offset)
    ENDIF

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_scale


  SUBROUTINE dbcsr_make_untransposed_blocks(matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_make_untransposed_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_size, row, &
                                                row_size
    INTEGER, DIMENSION(:), POINTER           :: cbs, rbs
    LOGICAL                                  :: sym_negation, tr
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_dp
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    rbs => array_data (dbcsr_row_block_sizes (matrix))
    cbs => array_data (dbcsr_col_block_sizes (matrix))
    CALL dbcsr_iterator_start (iter, matrix)
    sym_negation = matrix%m%sym_negation
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, r_dp, transposed=tr,&
            block_number=blk)
       IF (tr) THEN
          row_size = rbs(row)
          col_size = cbs(col)
          CALL dbcsr_block_transpose(r_dp, col_size, row_size)
          IF (sym_negation) THEN
             r_dp(:) = -r_dp(:)
          ENDIF
          matrix%m%blk_p(blk) = -matrix%m%blk_p(blk)
       ENDIF
    ENDDO
    CALL dbcsr_iterator_stop (iter)
  END SUBROUTINE dbcsr_make_untransposed_blocks


! *****************************************************************************
!> \brief Set a DBCSR matrix by alpha
!> \param[inout] matrix       DBCSR matrix
!> \param[in] alpha           a scalar
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_set(matrix, alpha, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(dp), INTENT(IN)                     :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, row, timing_handle
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    !
    !
    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       SELECT CASE (dbcsr_get_data_type(matrix))
       !CASE (dbcsr_type_real_4)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
          data_d(:,:) = alpha
       !CASE (dbcsr_type_complex_4)
       !CASE (dbcsr_type_complex_8)
       CASE DEFAULT
          CALL cp_unimplemented_error (routineN,"Only real double precision")
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_set

! *****************************************************************************
!> \brief Marks a block for removal from a DBCSR matrix. Handles
!>        symmetric matrices.
!> \param[in]  matrix         DBCSR matrix
!> \param[in]  row            row of block to remove
!> \param[in]  col            column of block to remove
!> \param[in]  block_number   (optional) the block number, if it is known
! *****************************************************************************
  SUBROUTINE dbcsr_remove_block(matrix, row, col, block_number)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    INTEGER, INTENT(IN), OPTIONAL            :: block_number

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_remove_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: b, c, r
    LOGICAL                                  :: found, tr

!   ---------------------------------------------------------------------------
!r = row
!c = col
! Figure out where the desired row and col is actually stored.
!CALL dbcsr_get_stored_coordinates (matrix, r, c, tr)
!IF (PRESENT (block_number)) THEN
!   b = block_number
!   CALL cp_assert (block_number .LE. matrix%m%nblks, cp_failure_level,&
!        cp_caller_error, routineN, "Block number too big.")
!   found = .TRUE.
!ELSE
!   CALL dbcsr_get_stored_block_info (matrix, r, c, found, b)
!ENDIF

    IF (PRESENT (block_number)) THEN
       b = block_number
       CALL cp_assert (block_number .LE. matrix%m%nblks, cp_failure_level,&
            cp_caller_error, routineN, "Block number too big.")
       found = .TRUE.
    ELSE
       CALL dbcsr_get_block_index (matrix, row, col, r, c, tr, found, b)
    ENDIF
    b = ABS (b)
    IF (found .AND. b .GT. 0) THEN
       ! Mark the block for deletion.
       matrix%m%blk_p(b) = 0
       matrix%m%valid = .FALSE.
    ELSE
       CALL cp_assert (b .NE. 0, cp_warning_level, cp_caller_error, routineN,&
            "Block does not exist or is already deleted.")
    ENDIF
  END SUBROUTINE dbcsr_remove_block

! *****************************************************************************
!> \brief add and scale matrices 
!>    A = alpha*A + beta*B or 
!> \param[in,out] matrix_a   DBCSR matrix
!> \param[in] matrix_b       DBCSR matrix
!> \param[in] alpha_scalar   (optional) 
!> \param[in] beta_scalar    (optional)
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_add(matrix_a, matrix_b, alpha_scalar, beta_scalar, &
                       first_column, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    REAL(dp), INTENT(IN), OPTIONAL           :: alpha_scalar, beta_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: first_column, last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, ibeg, iend, &
                                                my_first_column, &
                                                my_last_column, row, &
                                                timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: offset
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size
    LOGICAL                                  :: is_allocated, tr
    REAL(dp)                                 :: my_alpha_scalar, &
                                                my_beta_scalar
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d, data_d_tmp
    TYPE(dbcsr_iterator)                     :: iter

!real(KIND=sp), DIMENSION(:, :), POINTER  :: data_r
!complex(KIND=sp), DIMENSION(:, :), POINTER  :: data_z
!complex(KIND=dp), DIMENSION(:, :), POINTER  :: data_c
!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    IF(dbcsr_get_data_type(matrix_a).NE.dbcsr_get_data_type(matrix_b)) &
         CALL stop_program(routineN,"data types not consistent")

    !> this should be precomputed somewhere else
    col_blk_size => array_data (matrix_b%m%col_blk_size)
    ALLOCATE(offset(matrix_a%m%nblkcols_total+1))
    offset(1) = 1
    DO col = 2,matrix_a%m%nblkcols_total+1
       offset(col) = offset(col-1) + col_blk_size(col-1)
    ENDDO
    !<

    !
    my_alpha_scalar = 1.0_dp
    IF(PRESENT(alpha_scalar)) my_alpha_scalar = alpha_scalar
    my_beta_scalar = 1.0_dp
    IF(PRESENT(beta_scalar)) my_beta_scalar = beta_scalar
    !
    my_first_column = 1
    IF(PRESENT(first_column)) my_first_column = first_column
    my_last_column = dbcsr_nfullcols_total(matrix_a)
    IF(PRESENT(last_column)) my_last_column = last_column
    IF(dbcsr_nfullcols_total(matrix_a).NE.offset(matrix_a%m%nblkcols_total+1)-1) &
         CALL stop_program(routineN,"something wrong here 1")
    IF(my_first_column.LT.1) CALL stop_program(routineN,"first_column less than 1")
    IF(my_last_column.GT.dbcsr_nfullcols_total(matrix_a)) &
         CALL stop_program(routineN,"last_column greater than nfullcols_total")
    !
    ! let's go
    IF(matrix_a%m%nblkrows_total.NE.matrix_b%m%nblkrows_total) THEN
       CALL stop_program(routineN,"matrices not consistent")
    ENDIF

    CALL dbcsr_scale(matrix_a,alpha_scalar=my_alpha_scalar,error=error)

    IF(my_beta_scalar.NE.0.0_dp) THEN

       CALL dbcsr_iterator_start(iter, matrix_b)

       DO WHILE (dbcsr_iterator_blocks_left(iter))

          SELECT CASE (dbcsr_get_data_type(matrix_a))
          !CASE (dbcsr_type_real_4)
             !CALL dbcsr_iterator_next_block(iter, row, col, data_r, blk)
             !CALL dbcsr_put_block(matrix_a, row, col, data_r)
          CASE (dbcsr_type_real_8)
             CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)

             IF(offset(col+1).LT.my_first_column.OR.offset(col).GT.my_last_column) CYCLE

             data_d_tmp => data_d
             is_allocated = .FALSE.
             !
             IF(offset(col).LT.my_first_column.AND.offset(col+1).GT.my_first_column.AND.&
                offset(col+1).LT.my_last_column) THEN
                NULLIFY(data_d_tmp)
                ALLOCATE(data_d_tmp(SIZE(data_d,1),SIZE(data_d,2)))
                is_allocated = .TRUE.
                data_d_tmp = 0.0_dp
                ibeg = my_first_column - offset(col) + 1
                iend = SIZE(data_d_tmp,2)
                data_d_tmp(:,ibeg:iend) = data_d(:,ibeg:iend)
             ELSEIF(offset(col).LT.my_last_column.AND.offset(col+1).GT.my_last_column.AND.&
                offset(col).GT.my_first_column) THEN
                NULLIFY(data_d_tmp)
                ALLOCATE(data_d_tmp(SIZE(data_d,1),SIZE(data_d,2)))
                is_allocated = .TRUE.
                data_d_tmp = 0.0_dp
                ibeg = 1
                iend = my_last_column - offset(col) + 1
                data_d_tmp(:,ibeg:iend) = data_d(:,ibeg:iend)
             ELSEIF(offset(col).LE.my_first_column.AND.offset(col+1).GT.my_last_column.OR.&
                offset(col).LT.my_first_column.AND.offset(col+1).GE.my_last_column) THEN
                NULLIFY(data_d_tmp)
                ALLOCATE(data_d_tmp(SIZE(data_d,1),SIZE(data_d,2)))
                is_allocated = .TRUE.
                data_d_tmp = 0.0_dp
                ibeg = my_first_column - offset(col) + 1
                iend = my_last_column - offset(col) + 1
                data_d_tmp(:,ibeg:iend) = data_d(:,ibeg:iend)
             ENDIF
             !
             IF (my_beta_scalar.EQ.1.0_dp) THEN
                CALL dbcsr_put_block(matrix_a, row, col, data_d_tmp, tr,&
                     summation=.TRUE.)
             ELSE
                CALL dbcsr_put_block(matrix_a, row, col, data_d_tmp, tr,&
                     summation=.TRUE., scale=my_beta_scalar)
             ENDIF

             IF(is_allocated) DEALLOCATE(data_d_tmp)
             !CASE (dbcsr_type_complex_4)
             !   CALL dbcsr_iterator_next_block(iter, row, col, data_c, blk)
             !   CALL dbcsr_put_block(matrix_a, row, col, data_c)
             !CASE (dbcsr_type_complex_8)
             !   CALL dbcsr_iterator_next_block(iter, row, col, data_z, blk)
             !   CALL dbcsr_put_block(matrix_a, row, col, data_z)
          CASE DEFAULT
             CALL cp_unimplemented_error (routineN,"Only real double precision")
          END SELECT

       ENDDO

       CALL dbcsr_iterator_stop(iter)
       CALL dbcsr_finalize (matrix_a, error=error)

    ENDIF

    DEALLOCATE(offset)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_add

! *****************************************************************************
!> \brief Replaces blocks with blocks from another matrix.
!>
!> Every block in matrix B is copied over to matrix A, replacing any possibly
!> pre-existing blocks.
!> \param[in,out] matrix_a   replace blocks in this DBCSR matrix
!> \param[in] matrix_b       replace blocks taken from the DBCSR matrix
! *****************************************************************************
  SUBROUTINE dbcsr_replace_blocks(matrix_a, matrix_b)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_replace_blocks', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: data_c
    COMPLEX(KIND=sp), DIMENSION(:, :), &
      POINTER                                :: data_z
    INTEGER                                  :: blk, col, row, timing_handle
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: data_r
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    IF(dbcsr_get_data_type(matrix_a).NE.dbcsr_get_data_type(matrix_b)) &
         CALL stop_program(routineN,"data types not consistent")

    !
    ! let's go
    IF(matrix_a%m%nblkrows_total.NE.matrix_b%m%nblkrows_total) THEN
       CALL stop_program(routineN,"matrices not consistent")
    ENDIF

    CALL dbcsr_iterator_start(iter, matrix_b)

    DO WHILE (dbcsr_iterator_blocks_left(iter))

       SELECT CASE (dbcsr_get_data_type(matrix_a))
       CASE (dbcsr_type_real_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_r, tr, blk)
          CALL dbcsr_put_block(matrix_a, row, col, data_r)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
          CALL dbcsr_put_block(matrix_a, row, col, data_d, tr,&
               summation=.FALSE.)
       CASE (dbcsr_type_complex_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_c, tr, blk)
          CALL dbcsr_put_block(matrix_a, row, col, data_c)
       CASE (dbcsr_type_complex_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_z, tr, blk)
          CALL dbcsr_put_block(matrix_a, row, col, data_z)
       CASE DEFAULT
          CALL cp_unimplemented_error (routineN,"Unkown precision")
       END SELECT

    ENDDO

    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_finalize (matrix_a)
    !
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_replace_blocks

! *****************************************************************************
!> \brief add a constant to the diagonal of a matrix
!> \param[inout] matrix       DBCSR matrix
!> \param[in]    alpha_scalar scalar
!> \note  We assume the sparse matrix has at least its diagonal blocks.
! *****************************************************************************
  SUBROUTINE dbcsr_add_on_diag(matrix, alpha_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(dp), INTENT(IN)                     :: alpha_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_on_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, row, timing_handle
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL dbcsr_iterator_start(iter, matrix)

    DO WHILE (dbcsr_iterator_blocks_left(iter))

       SELECT CASE (dbcsr_get_data_type(matrix))
          !CASE (dbcsr_type_real_4)
          !   CALL dbcsr_iterator_next_block(iter, row, col, data_r, blk)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
          IF(row.EQ.col) CALL block_add_on_diag(SIZE(data_d,1),data_d,alpha_scalar)
          !CASE (dbcsr_type_complex_4)
          !   CALL dbcsr_iterator_next_block(iter, row, col, data_c, blk)
          !CASE (dbcsr_type_complex_8)
          !   CALL dbcsr_iterator_next_block(iter, row, col, data_z, blk)
       CASE DEFAULT
          CALL cp_unimplemented_error (routineN,"Only real double precision")
       END SELECT

    ENDDO

    CALL dbcsr_iterator_stop(iter)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_add_on_diag

! *****************************************************************************
!> \brief get the diagonal of a dbcsr matrix
!> \param[in] matrix    the matrix
!> \param[inout] diag   the diagonal
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_get_block_diag(matrix, diag, error)

    TYPE(dbcsr_obj), INTENT(IN), TARGET      :: matrix
    TYPE(dbcsr_obj), INTENT(INOUT), TARGET   :: diag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, row, timing_handle
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL dbcsr_create(diag, 'diag of '//matrix%m%name, &
         matrix%m%dist, matrix%m%matrix_type, matrix%m%row_blk_size,&
         matrix%m%col_blk_size, matrix%m%nblks, matrix%m%nze,&
         matrix%m%data_type, matrix%m%special_memory, error=error)

    CALL dbcsr_iterator_start(iter, matrix, contiguous_pointers=.TRUE.)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       SELECT CASE (dbcsr_get_data_type(matrix))
          !CASE (dbcsr_type_real_4)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
          IF(row.EQ.col) CALL dbcsr_put_block(diag, row, col, data_d, tr)
          !CASE (dbcsr_type_complex_4)
          !CASE (dbcsr_type_complex_8)
       CASE DEFAULT
          CALL cp_unimplemented_error (routineN,&
               "Only real double precision")
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_finalize(diag, error=error)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_get_block_diag

! *****************************************************************************
!> \brief get the diagonal of a dbcsr matrix
!> \param[in] matrix    the matrix
!> \param[inout] diag   diagonal
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_get_diag(matrix, diag, error)

    TYPE(dbcsr_obj), INTENT(IN), TARGET      :: matrix
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: diag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, i, row, &
                                                timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: offset
    INTEGER, DIMENSION(:), POINTER           :: row_blk_size
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    !> this should be precomputed somewhere else
    row_blk_size => array_data (matrix%m%row_blk_size)
    ALLOCATE(offset(matrix%m%nblkrows_total))
    offset(1) = 1
    DO row = 2,matrix%m%nblkrows_total
       offset(row) = offset(row-1) + row_blk_size(row-1)
    ENDDO
    !<
    IF(SUM(row_blk_size).GT.SIZE(diag)) CALL stop_program(routineN,"SIZE(diag) too small") 

    CALL dbcsr_iterator_start(iter, matrix)

    DO WHILE (dbcsr_iterator_blocks_left(iter))

       SELECT CASE (dbcsr_get_data_type(matrix))
          !CASE (dbcsr_type_real_4)
          !CALL dbcsr_iterator_next_block(iter, row, col, data_r, tr, blk)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
          IF(row.EQ.col) THEN
             DO i = 1,SIZE(data_d,1)
                diag(offset(row)+i-1) = data_d(i,i)
             ENDDO
          ENDIF
          !CASE (dbcsr_type_complex_4)
          !CALL dbcsr_iterator_next_block(iter, row, col, data_c, tr, blk)
          !CASE (dbcsr_type_complex_8)
          !CALL dbcsr_iterator_next_block(iter, row, col, data_z, tr, blk)
       CASE DEFAULT
          CALL cp_unimplemented_error (routineN,"Only real double precision")
       END SELECT

    ENDDO

    CALL dbcsr_iterator_stop(iter)

    DEALLOCATE(offset)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_get_diag

! *****************************************************************************
!> \brief set the diagonal of a dbcsr matrix
!> \param[in] matrix    the matrix
!> \param[inout] diag   diagonal
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_set_diag(matrix, diag, error)

    TYPE(dbcsr_obj), INTENT(INOUT), TARGET   :: matrix
    REAL(dp), DIMENSION(:), INTENT(IN)       :: diag
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, i, row, &
                                                timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: offset
    INTEGER, DIMENSION(:), POINTER           :: row_blk_size
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    !> this should be precomputed somewhere else
    row_blk_size => array_data (matrix%m%row_blk_size)
    ALLOCATE(offset(matrix%m%nblkrows_total))
    offset(1) = 1
    DO row = 2,matrix%m%nblkrows_total
       offset(row) = offset(row-1) + row_blk_size(row-1)
    ENDDO
    !<
    IF(SUM(row_blk_size).GT.SIZE(diag)) CALL stop_program(routineN,"SIZE(diag) too small")

    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       SELECT CASE (dbcsr_get_data_type(matrix))
          !CASE (dbcsr_type_real_4)
          !CALL dbcsr_iterator_next_block(iter, row, col, data_s, tr, blk)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
          IF(row.EQ.col) THEN
             DO i = 1,SIZE(data_d,1)
                data_d(i,i) = diag(offset(row)+i-1)
             ENDDO
          ENDIF
          !CASE (dbcsr_type_complex_4)
          !CALL dbcsr_iterator_next_block(iter, row, col, data_c, tr, blk)
          !CASE (dbcsr_type_complex_8)
          !CALL dbcsr_iterator_next_block(iter, row, col, data_z, tr, blk)
       CASE DEFAULT
          CALL cp_unimplemented_error (routineN,"Only real double precision")
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_set_diag

! *****************************************************************************
!> \brief copy a matrix 
!> \param[in]    matrix_a  DBCSR matrix
!> \param[inout] matrix_b  DBCSR matrix
!> \param[in]    name      (optional) name of the new matrix
!> \param[in,out] error    cp2k error
!> \param[in]    shallow_data  (optional) shallow data copy
! *****************************************************************************
  SUBROUTINE dbcsr_copy(matrix_b, matrix_a, name, error, shallow_data)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: shallow_data

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    LOGICAL                                  :: shallow

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    shallow = .FALSE. ; IF (PRESENT (shallow_data)) shallow = shallow_data

    IF(PRESENT(name)) THEN
       CALL dbcsr_create(matrix_b, name, &
            matrix_a%m%dist, matrix_a%m%matrix_type, matrix_a%m%row_blk_size,&
            matrix_a%m%col_blk_size, nblks=matrix_a%m%nblks,&
            data_type=matrix_a%m%data_type, special=matrix_a%m%special_memory,&
            error=error)
    ELSE
       CALL dbcsr_create(matrix_b,'copy of '//matrix_a%m%name, &
            matrix_a%m%dist, matrix_a%m%matrix_type, matrix_a%m%row_blk_size,&
            matrix_a%m%col_blk_size, nblks=matrix_a%m%nblks, &
            data_type=matrix_a%m%data_type, special=matrix_a%m%special_memory,&
            error=error)
    ENDIF
    CALL ensure_array_size(matrix_b%m%index, ub=SIZE(matrix_a%m%index), error=error)
    !
    ! copy index and data
    matrix_b%m%index(1:SIZE(matrix_a%m%index)) = matrix_a%m%index(:)
    IF (.NOT. shallow) THEN
       CALL dbcsr_data_ensure_size (matrix_b%m%data_area, matrix_a%m%nze)
       CALL dbcsr_data_copyall (matrix_b%m%data_area, matrix_a%m%data_area)
    ELSE
       CALL dbcsr_data_release (matrix_b%m%data_area)
       matrix_b%m%data_area = matrix_a%m%data_area
       CALL dbcsr_data_hold (matrix_b%m%data_area)
    ENDIF
    !
    ! the row_p, col_i and blk_p ...
    CALL dbcsr_repoint_index(matrix_b%m)
    matrix_b%m%nze = matrix_a%m%nze
    matrix_b%m%nblks = matrix_b%m%nblks
    matrix_b%m%valid = .TRUE.

    ! Now do cp2k special stuff
    IF (ASSOCIATED (matrix_a%m%dist2d)) THEN
       matrix_b%m%dist2d => matrix_a%m%dist2d
       CALL distribution_2d_retain (matrix_b%m%dist2d, error)
    ENDIF
    matrix_b%m%sparsity_id = matrix_a%m%sparsity_id

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_copy


! *****************************************************************************
!> \brief block triu of a dbcsr matrix
!> \param[in] matrix_b  triu of the matrix
!> \param[in] matrix_a  the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_btriu(matrix_b, matrix_a, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_btriu', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, row, timing_handle
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL dbcsr_copy(matrix_b, matrix_a, name="triu of "//matrix_a%m%name, error=error)

    CALL dbcsr_iterator_start(iter, matrix_b)

    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, blk)
       IF(row.GT.col) CALL dbcsr_remove_block(matrix_b, row, col, blk)
    ENDDO

    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_finalize (matrix_b, error=error)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_btriu

! *****************************************************************************
!> \brief filter a dbcsr matrix
!> \param[inout] matrix_a  the matrix
!> \param[in] eps          the threshold
!> \param[in] method       how the matrix is filtered 
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_filter(matrix, eps, method, error)

    TYPE(dbcsr_obj), INTENT(INOUT), TARGET   :: matrix
    REAL(dp), INTENT(IN)                     :: eps
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: method
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_filter', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=3)                         :: my_method
    INTEGER                                  :: blk, col, row, timing_handle
    LOGICAL                                  :: tr
    REAL(dp)                                 :: ddot, nrm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    my_method = 'fro'
    IF(PRESENT(method)) my_method = method

    CALL dbcsr_iterator_start(iter, matrix, contiguous_pointers=.TRUE.)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       SELECT CASE (dbcsr_get_data_type(matrix))
          !CASE (dbcsr_type_real_4)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
          SELECT CASE(my_method)
          CASE('fro')
             !
             ! Frobenius based
             nrm = SQRT(DDOT(SIZE(data_d), data_d(1), 1, data_d(1), 1))
             IF(nrm.LT.eps) CALL dbcsr_remove_block(matrix, row, col, blk)
          CASE DEFAULT
             CALL cp_unimplemented_error (routineN,&
                  "Only Frobenius based filtering")
          END SELECT
          !CASE (dbcsr_type_complex_4)
          !CASE (dbcsr_type_complex_8)
       CASE DEFAULT
          CALL cp_unimplemented_error (routineN,&
               "Only real double precision")
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_finalize (matrix, error=error)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_filter

! *****************************************************************************
!> \brief compute a norm of a dbcsr matrix
!> \param[in] matrix  the matrix
!> \param[out] norm
!> \param[in] which_norm
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_norm(matrix, which_norm, norm_scalar, norm_vector, error)

    TYPE(dbcsr_obj), INTENT(INOUT), TARGET   :: matrix
    CHARACTER(LEN=*), INTENT(IN)             :: which_norm
    REAL(dp), INTENT(OUT), OPTIONAL          :: norm_scalar
    REAL(dp), DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: norm_vector
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_norm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, i, j, nc, nr, row, &
                                                timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: c_offset, r_offset
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: tr
    REAL(dp)                                 :: fac
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff_d
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    col_blk_size => array_data (matrix%m%col_blk_size)
    row_blk_size => array_data (matrix%m%row_blk_size)

    !> this should be precomputed somewhere else
    ALLOCATE(r_offset(matrix%m%nblkrows_total),c_offset(matrix%m%nblkcols_total))
    r_offset(1) = 1
    DO row = 2,matrix%m%nblkrows_total
       r_offset(row) = r_offset(row-1) + row_blk_size(row-1)
    ENDDO
    nr = SUM(row_blk_size)
    c_offset(1) = 1
    DO col = 2,matrix%m%nblkcols_total
       c_offset(col) = c_offset(col-1) + col_blk_size(col-1)
    ENDDO
    nc = SUM(col_blk_size)
    !<

    SELECT CASE(which_norm)

    CASE('Frobenius')

       IF(PRESENT(norm_scalar)) THEN
          norm_scalar = 0.0_dp
          CALL dbcsr_iterator_start(iter, matrix)
          DO WHILE (dbcsr_iterator_blocks_left(iter))
             SELECT CASE (dbcsr_get_data_type(matrix))
                !CASE (dbcsr_type_real_4)
             CASE (dbcsr_type_real_8)
                CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
                fac = 1.0_dp
                IF((dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_symmetric.OR.&
                    dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_antisymmetric).AND.&
                    row.NE.col) fac = 2.0_dp
                norm_scalar = norm_scalar + fac * SUM(data_d**2)
                !CASE (dbcsr_type_complex_4)
                !CASE (dbcsr_type_complex_8)
             CASE DEFAULT
                CALL cp_unimplemented_error (routineN,"Only real double precision")
             END SELECT
          ENDDO
          CALL dbcsr_iterator_stop(iter)
          CALL mp_max(norm_scalar,dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))
          norm_scalar = SQRT(norm_scalar)
       ENDIF

    CASE('MaxAbsNorm')

       IF(PRESENT(norm_scalar)) THEN
          norm_scalar = 0.0_dp
          CALL dbcsr_iterator_start(iter, matrix)
          DO WHILE (dbcsr_iterator_blocks_left(iter))
             SELECT CASE (dbcsr_get_data_type(matrix))
                !CASE (dbcsr_type_real_4)
             CASE (dbcsr_type_real_8)
                CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
                norm_scalar = MAX(norm_scalar,MAXVAL(ABS(data_d)))
                !CASE (dbcsr_type_complex_4)
                !CASE (dbcsr_type_complex_8)
             CASE DEFAULT
                CALL cp_unimplemented_error (routineN,"Only real double precision")
             END SELECT
          ENDDO
          CALL dbcsr_iterator_stop(iter)
          CALL mp_max(norm_scalar,dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))
       ENDIF

    CASE('Gershgorin')

       IF(PRESENT(norm_scalar)) THEN
          IF(nr.NE.nc) CALL stop_program(routineN,"not a square matrix")
          norm_scalar = 0.0_dp
          ALLOCATE(buff_d(nr))
          buff_d = 0.0_dp
          CALL dbcsr_iterator_start(iter, matrix)
          DO WHILE (dbcsr_iterator_blocks_left(iter))
             SELECT CASE (dbcsr_get_data_type(matrix))
             !CASE (dbcsr_type_real_4)
             CASE (dbcsr_type_real_8)
                CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
                IF((dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_symmetric.OR.&
                    dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_antisymmetric).AND.&
                    row.NE.col) THEN
                   DO j=1,SIZE(data_d,2)
                   DO i=1,SIZE(data_d,1)
                      buff_d(r_offset(row)+i-1) = buff_d(r_offset(row)+i-1) + ABS(data_d(i,j))
                      buff_d(c_offset(col)+j-1) = buff_d(c_offset(col)+j-1) + ABS(data_d(i,j))
                   ENDDO
                   ENDDO
                ELSE
                   DO j=1,SIZE(data_d,2)
                   DO i=1,SIZE(data_d,1)
                      buff_d(r_offset(row)+i-1) = buff_d(r_offset(row)+i-1) + ABS(data_d(i,j))
                   ENDDO
                   ENDDO
                ENDIF
             !CASE (dbcsr_type_complex_4)
             !CASE (dbcsr_type_complex_8)
             CASE DEFAULT
                CALL cp_unimplemented_error (routineN,"Only real double precision")
             END SELECT
          ENDDO
          CALL dbcsr_iterator_stop(iter)
          CALL mp_sum(buff_d,dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))       
          norm_scalar = MAXVAL(buff_d)
          DEALLOCATE(buff_d)
       ENDIF

    CASE('Column')

       IF(PRESENT(norm_vector)) THEN
          norm_vector(:) = 0.0_dp
          CALL dbcsr_iterator_start(iter, matrix)
          DO WHILE (dbcsr_iterator_blocks_left(iter))
             SELECT CASE (dbcsr_get_data_type(matrix))
             !CASE (dbcsr_type_real_4)
             CASE (dbcsr_type_real_8)
                CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
                IF((dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_symmetric.OR.&
                    dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_antisymmetric).AND.&
                    row.NE.col) THEN
                   DO j=1,SIZE(data_d,2)
                   DO i=1,SIZE(data_d,1)
                      norm_vector(c_offset(col)+j-1) = norm_vector(c_offset(col)+j-1) &
                                                     + data_d(i,j)**2
                      norm_vector(r_offset(row)+i-1) = norm_vector(r_offset(row)+i-1) &
                                                     + data_d(i,j)**2
                   ENDDO
                   ENDDO
                ELSE
                   DO j=1,SIZE(data_d,2)
                   DO i=1,SIZE(data_d,1)
                      norm_vector(c_offset(col)+j-1) = norm_vector(c_offset(col)+j-1) &
                                                     + data_d(i,j) * data_d(i,j)
                   ENDDO
                   ENDDO
                ENDIF
             !CASE (dbcsr_type_complex_4)
             !CASE (dbcsr_type_complex_8)
             CASE DEFAULT
                CALL cp_unimplemented_error (routineN,"Only real double precision")
             END SELECT
          ENDDO
          CALL dbcsr_iterator_stop(iter)
          CALL mp_sum(norm_vector,dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))
          norm_vector = SQRT(norm_vector)
       ENDIF

    CASE DEFAULT

       CALL stop_program(routineN,"this norm is NYI")

    END SELECT

    DEALLOCATE(r_offset,c_offset)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_norm


! *****************************************************************************
!> \brief Replicates a DBCSR on all processors.
!> \param[in,out] matrix      matrix to replicate
! *****************************************************************************
  SUBROUTINE dbcsr_replicate_all(matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_replicate_all', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 2

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: c_dp
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: c_sp
    INTEGER :: blk, blk_l, blk_ps, blks, col, dst_p, mp_group, numnodes, nze, &
      row, src_p, stored_col, stored_row, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rd_disp, recv_meta, rm_disp, &
                                                sd_disp, sdp, send_meta, &
                                                sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_count, send_count, &
                                                total_recv_count, &
                                                total_send_count
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, col_dist, &
                                                row_blk_size, row_dist, &
                                                tmp_index
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    LOGICAL                                  :: tr
    REAL(KIND=dp)                            :: tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: send_data
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: recv_data
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_dp
    REAL(KIND=sp), DIMENSION(:), POINTER     :: r_sp
    TYPE(dbcsr_data_obj)                     :: tmp_data
    TYPE(dbcsr_distribution_obj)             :: target_dist
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_obj)                          :: replicated

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    tstart = m_walltime ()
    CALL cp_assert (dbcsr_valid_index(matrix%m),&
         cp_fatal_level, cp_caller_error,&
         routineN, "Matrix not initialized.")
    CALL cp_assert (matrix%m%replication_type .EQ. dbcsr_repl_none,&
         cp_warning_level, cp_caller_error, routineN,&
         "Replicating a non-distributed matrix makes no sense.")
    SELECT CASE (matrix%m%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (matrix%m%data_area, r_dp)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (matrix%m%data_area, r_sp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (matrix%m%data_area, c_dp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (matrix%m%data_area, c_sp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    END SELECT
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    target_dist = matrix%m%dist
    row_dist => array_data (dbcsr_distribution_row_dist (target_dist))
    col_dist => array_data (dbcsr_distribution_col_dist (target_dist))
    mp_obj = dbcsr_distribution_mp (target_dist)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    numnodes = dbcsr_mp_numnodes (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    DBGV 'row_dist', row_dist
    DBGV 'col_dist', col_dist
    DBG 'row, col maxval:', MAXVAL(row_dist), MAXVAL(col_dist)
    DBG 'blacs2mpi',blacs2mpi
    DBG 'blacs2mpi LB',LBOUND(blacs2mpi)
    DBG 'blacs2mpi UB',UBOUND(blacs2mpi)
    CALL cp_assert(MAXVAL(row_dist).LE.UBOUND(blacs2mpi,1), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows')
    IF (dbg) THEN
       CALL cp_assert(MAXVAL(row_dist).EQ.UBOUND(blacs2mpi,1),&
            cp_warning_level, cp_wrong_args_error, routineN,&
            'Range of row distribution not equal to processor rows')
    ENDIF
    CALL cp_assert(MAXVAL(col_dist).LE.UBOUND(blacs2mpi,2), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols')
    IF (dbg) THEN
       CALL cp_assert(MAXVAL(col_dist).EQ.UBOUND(blacs2mpi,2),&
            cp_warning_level, cp_wrong_args_error, routineN,&
            'Range of col distribution not equal to processor cols')
    ENDIF
    ALLOCATE (send_count(2, 0:numnodes-1))
    ALLOCATE (recv_count(2, 0:numnodes-1))
    ALLOCATE (total_send_count(2, 0:numnodes-1))
    ALLOCATE (total_recv_count(2, 0:numnodes-1))
    ALLOCATE (sdp(0:numnodes-1))
    ALLOCATE (sd_disp(0:numnodes-1))
    ALLOCATE (smp(0:numnodes-1))
    ALLOCATE (sm_disp(0:numnodes-1))
    ALLOCATE (rd_disp(0:numnodes-1))
    ALLOCATE (rm_disp(0:numnodes-1))
    CALL dbcsr_init (replicated)
    CALL dbcsr_create(replicated, 'Replicated '//matrix%m%name, matrix%m%dist,&
         'N', matrix%m%row_blk_size, matrix%m%col_blk_size,&
         0, 0, matrix%m%data_type, special=matrix%m%special_memory,&
         replication_type=dbcsr_repl_full)
    ! Count initial sizes for sending.
    DBG 'A'
    send_count(:,:) = 0
    CALL dbcsr_iterator_start (iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, r_dp, tr, blk)
       nze = SIZE(r_dp)
       DO dst_p = 0, numnodes-1
          send_count(1, dst_p) = send_count(1, dst_p) + 1
          send_count(2, dst_p) = send_count(2, dst_p) + nze
       END DO
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    DBG 'B'
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    DBG 'C'
    DBG 'send counts',send_count
    DBG 'recv counts',recv_count
    ! Allocate data structures needed for data exchange.
    ALLOCATE (recv_data(SUM(recv_count(2, :))))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :))))
    ALLOCATE (send_data(SUM(send_count(2, :))))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :))))
    DBG 'send data len',SIZE (send_data)
    DBG 'send meta len',SIZE (send_meta)
    DBG 'recv data len',SIZE (recv_data)
    DBG 'recv meta len',SIZE (recv_meta)
    ! Fill in the meta data structures and copy the data.
    DO dst_p = 0, numnodes-1
       total_send_count(1, dst_p) = send_count (1, dst_p)
       total_send_count(2, dst_p) = send_count (2, dst_p)
       total_recv_count(1, dst_p) = recv_count (1, dst_p)
       total_recv_count(2, dst_p) = recv_count (2, dst_p)
    ENDDO
    sd_disp = -1 ; sm_disp = -1
    rd_disp = -1 ; rm_disp = -1
    sd_disp(0) = 1 ; sm_disp(0) = 1
    rd_disp(0) = 1 ; rm_disp(0) = 1
    DO dst_p = 1, numnodes-1
       sm_disp(dst_p) = sm_disp(dst_p-1)&
                        + metalen*total_send_count(1, dst_p-1)
       sd_disp(dst_p) = sd_disp(dst_p-1)&
                        + total_send_count(2, dst_p-1)
       rm_disp(dst_p) = rm_disp(dst_p-1)&
                        + metalen*total_recv_count(1, dst_p-1)
       rd_disp(dst_p) = rd_disp(dst_p-1)&
                        + total_recv_count(2, dst_p-1)
    ENDDO
    !DO dst_p = 0, numnodes-1
    !      sm_disp (dst_p) = sm_disp (dst_p) + metalen*send_count(1, dst_p)
    !      sd_disp (dst_p) = sd_disp (dst_p) + send_count(2, dst_p)
    !      rm_disp (dst_p) = rm_disp (dst_p) + metalen*recv_count(1, dst_p)
    !      rd_disp (dst_p) = rd_disp (dst_p) + recv_count(2, dst_p)
    !ENDDO
    sdp(:) = sd_disp
    smp(:) = sm_disp
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    DBG 'D'
    CALL dbcsr_iterator_start (iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, r_dp, tr, blk)
       nze = SIZE(r_dp)
       DO dst_p = 0, numnodes-1
          send_meta(smp(dst_p)+1) = col
          IF (.NOT. tr) THEN
             send_meta(smp(dst_p)) = row
             send_data(sdp(dst_p):sdp(dst_p)+nze-1) =&
                  r_dp(:)
          ELSE
             send_meta(smp(dst_p)) = -row ! Mark block transposed
             send_data(sdp(dst_p):sdp(dst_p)+nze-1) =&
                  r_dp(:)
          ENDIF
          smp(dst_p) = smp(dst_p) + metalen
          sdp(dst_p) = sdp(dst_p) + nze
       END DO
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    ! Exchange the data and metadata structures.
    DBG 'E'
    CALL mp_alltoall(send_data(:), total_send_count(2,:), sd_disp(:)-1,&
         recv_data(:), total_recv_count(2,:), rd_disp(:)-1, mp_group)
    DBG 'F'
    DBG 'send_meta',send_meta
    CALL mp_alltoall(send_meta(:), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), metalen*total_recv_count(1,:), rm_disp(:)-1, mp_group)
    DBG 'recv_meta',recv_meta
    DBG 'G'
    ! Now fill in the data.
    CALL dbcsr_work_create(replicated,&
            SUM(recv_count(1,:)),&
            SUM(recv_count(2,:)))
    !DBG 'send disps data',sd_disp
    !DBG 'send disps meta',sm_disp
    !DBG 'recv disps data',rd_disp
    !DBG 'recv disps meta',rm_disp
    !blk_p = 1
    !blk = 1
    blk_ps = 1
    blks = 1
    DO src_p = 0, numnodes-1
       nze = recv_count(2, src_p)
       CALL dbcsr_data_set (replicated%m%wms(1)%data_area, blk_ps, nze,&
            recv_data, rd_disp(src_p))
       DO blk_l = 1, recv_count(1, src_p)
          stored_row = recv_meta(rm_disp(src_p)+metalen*(blk_l-1))
          stored_col = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+1)
          DBG 'blk: p',blks,src_p,'; row, col',&
               stored_row,stored_col,&
               '; blk_ps',blk_ps
          replicated%m%wms(1)%row_i(blks) = ABS(stored_row)
          replicated%m%wms(1)%col_i(blks) = stored_col
          replicated%m%wms(1)%blk_p(blks) = SIGN(blk_ps, stored_row)
          nze = row_blk_size(ABS(stored_row))&
               * col_blk_size(stored_col)
          blk_ps = blk_ps + nze
          blks = blks + 1
       ENDDO
    ENDDO
    !
    replicated%m%wms(1)%lastblk = blks - 1
    replicated%m%wms(1)%datasize = blk_ps - 1
    DBG 'Finalizing normalization'
    CALL dbcsr_finalize(replicated, resort=.TRUE., reshuffle=.TRUE.)
    tstop = m_walltime ()
    DBG 'time:',tstop-tstart
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_data)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_data)
    DEALLOCATE(send_meta)
    matrix%m%replication_type = dbcsr_repl_full
    ! Now replace the data and index
    tmp_data = matrix%m%data_area
    matrix%m%data_area = replicated%m%data_area
    replicated%m%data_area = tmp_data
    tmp_index => matrix%m%index
    matrix%m%index => replicated%m%index
    replicated%m%index => tmp_index
    CALL dbcsr_repoint_index (matrix%m)
    matrix%m%nze = replicated%m%nze
    matrix%m%nblks = replicated%m%nblks
    CALL dbcsr_release (replicated)
    DBG 'Done replicating'
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_replicate_all


! *****************************************************************************
!> \brief Finds block to which a full element belongs.
!> \par Assumptions
!>        It is assumed that block_start and block_end are sorted and
!>        that hint is in the range [0, nblocks].
!> \param[in] full            full element
!> \param[out] block          block to which full belongs
!> \param[in] block_start     starting full elements of blocks
!> \param[out] block_end      ending full elements of blocks
!> \param[in] hint            where to start looking; ignored if 0
! *****************************************************************************
  SUBROUTINE find_block_of_element(full, block, nblocks,&
       block_start, block_end, hint)
    INTEGER, INTENT(in)                      :: full
    INTEGER, INTENT(out)                     :: block
    INTEGER, INTENT(in)                      :: nblocks
    INTEGER, DIMENSION(1:nblocks), &
      INTENT(in)                             :: block_start, block_end
    INTEGER, INTENT(in)                      :: hint

    INTEGER                                  :: count

    IF (hint .NE. 0) THEN
       block = hint
    ELSE
       block = MAX(1,(nblocks+1)/2)
    ENDIF
    count = 0
    DO WHILE (block_start(block) .GT. full .OR. block_end(block) .LT. full)
       IF (block_start(block) .GT. full) THEN
          block = block-1
       ELSEIF (block_end(block) .LT. full) THEN
          block = block+1
       ENDIF
       count = count+1
       IF (count .GT. nblocks .OR. block .LT. 1 .OR. block .GT. nblocks) THEN
          WRITE(*,'(1X,A,I9,A,I9,A)')"Want to find block",&
               block," of",nblocks," blocks"
          !WRITE(*,'(10(1X,I7))')block_start
          !WRITE(*,'(10(1X,I7))')block_end
          !WRITE(*,*)"full, nblocks, count, hint", full, nblocks, count, hint
          CALL cp_assert (count .LE. nblocks, cp_failure_level,&
               cp_internal_error, "find_block_of_element",&
               "Too many searches")
       ENDIF
    ENDDO
  END SUBROUTINE find_block_of_element


! *****************************************************************************
!> \brief Fully redistributes a DBCSR matrix.
!>
!>        The new distribution may be arbitrary as long as the total
!>        number full rows and columns matches that of the existing
!>        matrix.
!> \param[in] matrix          matrix to redistribute
!> \param[in,out] redist      redistributed matrix
!> \param[in] alpha           (optional) scaling of input matrix
!> \param[in] keep_sparsity   (optional) retains the sparsity of the redist
!>                            matrix
! *****************************************************************************
  SUBROUTINE dbcsr_complete_redistribute(matrix, redist, alpha, keep_sparsity)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_obj), INTENT(INOUT)           :: redist
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: alpha
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_complete_redistribute', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 7

    CHARACTER(len=80)                        :: ifmt
    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: c_dp
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: c_sp
    INTEGER :: blk, blk_col_new, blk_col_prev_new, blk_ps, blk_row_new, &
      blk_row_prev_new, blks, col, col_offset_old, col_rle, data_offset_l, &
      dst_p, full_col, full_row, i, meta_l, mp_group, nblkcols_total_new, &
      nblkrows_total_new, numnodes, rel_full_col_new, rel_full_row_new, row, &
      row_offset_old, row_prev, row_rle, src_p, stored_col_new, &
      stored_row_new, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:) :: col_end_new, col_end_old, &
      col_start_new, col_start_old, rd_disp, recv_meta, rm_disp, row_end_new, &
      row_end_old, row_start_new, row_start_old, sd_disp, sdp, send_meta, &
      sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_count, send_count, &
                                                total_recv_count, &
                                                total_send_count
    INTEGER, DIMENSION(:), POINTER :: col_blk_size_new, col_blk_size_old, &
      col_dist_new, col_dist_old, row_blk_size_new, row_blk_size_old, &
      row_dist_new, row_dist_old
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    LOGICAL                                  :: found, my_keep_sparsity, &
                                                new_meta, tr, valid_block
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: send_data
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: recv_data
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_dp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: r2_dp
    REAL(KIND=sp), DIMENSION(:), POINTER     :: r_sp
    TYPE(dbcsr_distribution_obj)             :: dist_new, dist_old
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj_new
    INTEGER, DIMENSION(:), ALLOCATABLE       :: map_row_el_blk, map_col_el_blk

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    CALL cp_assert (dbcsr_valid_index(matrix%m),&
         cp_fatal_level, cp_caller_error,&
         routineN, "Input not valid.")
    CALL cp_assert (dbcsr_is_initialized(redist%m),&
         cp_fatal_level, cp_caller_error,&
         routineN, "Output matrix not initialized.")
    CALL cp_assert (matrix%m%replication_type .EQ. dbcsr_repl_none,&
         cp_warning_level, cp_caller_error, routineN,&
         "Can not redistribute replicated matrix.")
    SELECT CASE (matrix%m%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (matrix%m%data_area, r_dp)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (matrix%m%data_area, r_sp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (matrix%m%data_area, c_dp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (matrix%m%data_area, c_sp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    END SELECT
    !
    my_keep_sparsity = .FALSE.
    IF (PRESENT (keep_sparsity)) my_keep_sparsity = my_keep_sparsity
    ! Get row and column start and end positions
    ! Old matrix
    dist_old = dbcsr_distribution (matrix)
    row_blk_size_old => array_data (matrix%m%row_blk_size)
    col_blk_size_old => array_data (matrix%m%col_blk_size)
    ALLOCATE (row_start_old (dbcsr_nblkrows_total (matrix)),&
         row_end_old (dbcsr_nblkrows_total (matrix)),&
         col_start_old (dbcsr_nblkcols_total (matrix)),&
         col_end_old (dbcsr_nblkcols_total (matrix)))
    CALL convert_sizes_to_offsets (row_blk_size_old,&
         row_start_old, row_end_old)
    CALL convert_sizes_to_offsets (col_blk_size_old,&
         col_start_old, col_end_old)
    row_dist_old => array_data (dbcsr_distribution_row_dist (dist_old))
    col_dist_old => array_data (dbcsr_distribution_col_dist (dist_old))
    ! New matrix
    dist_new = dbcsr_distribution (redist)
    row_blk_size_new => array_data (redist%m%row_blk_size)
    col_blk_size_new => array_data (redist%m%col_blk_size)
    ALLOCATE (row_start_new (dbcsr_nblkrows_total (redist)),&
         row_end_new (dbcsr_nblkrows_total (redist)),&
         col_start_new (dbcsr_nblkcols_total (redist)),&
         col_end_new (dbcsr_nblkcols_total (redist)))
    CALL convert_sizes_to_offsets (row_blk_size_new,&
         row_start_new, row_end_new)
    CALL convert_sizes_to_offsets (col_blk_size_new,&
         col_start_new, col_end_new)
    row_dist_new => array_data (dbcsr_distribution_row_dist (dist_new))
    col_dist_new => array_data (dbcsr_distribution_col_dist (dist_new))
    nblkrows_total_new = dbcsr_nblkrows_total (redist)
    nblkcols_total_new = dbcsr_nblkcols_total (redist)
    ! Create mappings of elements to new blocks.
    ALLOCATE (map_row_el_blk (dbcsr_nfullrows_total (redist)))
    CALL fill_map (map_row_el_blk, dbcsr_nfullrows_total (redist),&
         row_start_new, nblkrows_total_new,&
         row_blk_size_new(nblkrows_total_new))
    ALLOCATE (map_col_el_blk (dbcsr_nfullcols_total (redist)))
    CALL fill_map (map_col_el_blk, dbcsr_nfullcols_total (redist),&
         col_start_new, nblkcols_total_new,&
         col_blk_size_new(nblkcols_total_new))
    !
    mp_obj_new = dbcsr_distribution_mp (dist_new)
    pgrid => dbcsr_mp_pgrid (mp_obj_new)
    numnodes = dbcsr_mp_numnodes (mp_obj_new)
    mp_group = dbcsr_mp_group (mp_obj_new)
    !
    CALL cp_assert(MAXVAL(row_dist_new).LE.UBOUND(pgrid,1), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows')
    IF (dbg) THEN
       CALL cp_assert(MAXVAL(row_dist_new).EQ.UBOUND(pgrid,1),&
            cp_warning_level, cp_wrong_args_error, routineN,&
            'Range of row distribution not equal to processor rows')
    ENDIF
    CALL cp_assert(MAXVAL(col_dist_new).LE.UBOUND(pgrid,2), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols')
    IF (dbg) THEN
       CALL cp_assert(MAXVAL(col_dist_new).EQ.UBOUND(pgrid,2),&
            cp_warning_level, cp_wrong_args_error, routineN,&
            'Range of col distribution not equal to processor cols')
    ENDIF
    ALLOCATE (send_count(2, 0:numnodes-1))
    ALLOCATE (recv_count(2, 0:numnodes-1))
    ALLOCATE (total_send_count(2, 0:numnodes-1))
    ALLOCATE (total_recv_count(2, 0:numnodes-1))
    ALLOCATE (sdp(0:numnodes-1))
    ALLOCATE (sd_disp(0:numnodes-1))
    ALLOCATE (smp(0:numnodes-1))
    ALLOCATE (sm_disp(0:numnodes-1))
    ALLOCATE (rd_disp(0:numnodes-1))
    ALLOCATE (rm_disp(0:numnodes-1))
    ! Count initial sizes for sending.
    !
    ! We go through every element of every local block and determine
    ! to which processor it must be sent. It could be more efficient,
    ! but at least the index data are run-length encoded.
    DBG 'A'
    send_count(:,:) = 0
    CALL dbcsr_iterator_start (iter, matrix)
    blk_col_prev_new = 0
    blk_row_prev_new = 0
    row_prev = 0
    dst_p = -1
    new_meta = .TRUE.
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk)
       IF (row .NE. row_prev) blk_col_prev_new = 0 ! Reset the "hint"
       ! For all new blocked columns that overlap the old block colum,
       ! count the data that should be sent to that column.
       DO full_col = col_start_old(col), col_end_old(col)
          blk_col_new = map_col_el_blk(full_col)
          DO full_row = row_start_old(row), row_end_old(row)
             blk_row_new = map_row_el_blk(full_row)
             IF (blk_col_new .NE. blk_col_prev_new&
                  .OR. blk_row_new .NE. blk_row_prev_new .OR. new_meta) THEN
                CALL dbcsr_get_stored_coordinates (redist%m,&
                     blk_row_new, blk_col_new, tr, dst_p)
                send_count(1, dst_p) = send_count(1, dst_p) + 1
             ENDIF
             send_count(2, dst_p) = send_count(2, dst_p) + 1
          ENDDO
          new_meta = .TRUE.
       ENDDO
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    DBG 'B'
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    DBG 'C'
    DBG 'send counts',send_count
    DBG 'recv counts',recv_count
    ! Allocate data structures needed for data exchange.
    ALLOCATE (recv_data(SUM(recv_count(2, :))))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :))))
    ALLOCATE (send_data(SUM(send_count(2, :))))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :))))
    DBG 'send data len',SIZE (send_data)
    DBG 'send meta len',SIZE (send_meta)
    DBG 'recv data len',SIZE (recv_data)
    DBG 'recv meta len',SIZE (recv_meta)
    ! Fill in the meta data structures and copy the data.
    DO dst_p = 0, numnodes-1
       total_send_count(1, dst_p) = send_count (1, dst_p)
       total_send_count(2, dst_p) = send_count (2, dst_p)
       total_recv_count(1, dst_p) = recv_count (1, dst_p)
       total_recv_count(2, dst_p) = recv_count (2, dst_p)
    ENDDO
    sd_disp = -1 ; sm_disp = -1
    rd_disp = -1 ; rm_disp = -1
    sd_disp(0) = 1 ; sm_disp(0) = 1
    rd_disp(0) = 1 ; rm_disp(0) = 1
    DO dst_p = 1, numnodes-1
       sm_disp(dst_p) = sm_disp(dst_p-1)&
                        + metalen*total_send_count(1, dst_p-1)
       sd_disp(dst_p) = sd_disp(dst_p-1)&
                        + total_send_count(2, dst_p-1)
       rm_disp(dst_p) = rm_disp(dst_p-1)&
                        + metalen*total_recv_count(1, dst_p-1)
       rd_disp(dst_p) = rd_disp(dst_p-1)&
                        + total_recv_count(2, dst_p-1)
    ENDDO
    sdp(:) = sd_disp     ! sdp points to the the next place to store
                         ! data. It is postincremented.
    smp(:) = sm_disp - metalen  ! But smp points to the "working" data, not
                                ! the next. It is pre-incremented, so we must
                                ! first rewind it.
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    DBG 'D'
    IF (dbg) THEN
       WRITE(*,*)'row_start_old',row_start_old
       WRITE(*,*)'row_start_new',row_start_new
       WRITE(*,*)'row_blk_size_new',row_blk_size_new
    ENDIF
    CALL dbcsr_iterator_start (iter, matrix)
    blk_col_prev_new = 0
    blk_row_prev_new = 0
    row_prev = 0
    dst_p = -1
    new_meta = .TRUE.
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, r2_dp, tr, blk)
       IF (row .NE. row_prev) blk_col_prev_new = 0 ! Reset the "hint"
       !IF (tr) WRITE(*,*)"block at",row,col," is transposed"
       col_offset_old = col_start_old(col) - 1
       DO full_col = col_start_old(col), col_end_old(col)
          blk_col_new = map_col_el_blk(full_col)
          row_offset_old = row_start_old(row) - 1
          IF (dbg) THEN
             WRITE(*,*)'col',full_col,'of',&
                  col_start_old(col),col_end_old(col),'/',SIZE(r2_dp,2)
             WRITE(*,*)'rows',0,'of',&
                  row_start_old(row),row_end_old(row),'/',SIZE(r2_dp,1)
             DO i = 1, SIZE(r2_dp,1)
                WRITE(ifmt,'("(",I2,"(1X,F5.2))")')SIZE(r2_dp,2)
                WRITE(*,ifmt)r2_dp(i,:)
             ENDDO
          ENDIF
          DO full_row = row_start_old(row), row_end_old(row)
             blk_row_new = map_row_el_blk(full_row)
             IF (dbg) &
                  WRITE(*,*)'using block',blk_row_new,'x',blk_col_new
             IF (blk_col_new .NE. blk_col_prev_new&
                  .OR. blk_row_new .NE. blk_row_prev_new .OR. new_meta) THEN
                ! Start a new RLE run
                CALL dbcsr_get_stored_coordinates (redist%m,&
                     blk_row_new, blk_col_new, tr, dst_p)
                send_count(1, dst_p) = send_count(1, dst_p) + 1
                rel_full_row_new = full_row - row_start_new (blk_row_new)+1
                rel_full_col_new = full_col - col_start_new (blk_col_new)+1
                smp(dst_p) = smp(dst_p) + metalen
                send_meta(smp(dst_p)) = blk_row_new   ! new blocked row
                send_meta(smp(dst_p)+1) = blk_col_new ! new blocked column
                send_meta(smp(dst_p)+2) = rel_full_row_new  ! row in new block
                send_meta(smp(dst_p)+3) = rel_full_col_new  ! col in new block
                send_meta(smp(dst_p)+4) = 1           ! RLE rows
                send_meta(smp(dst_p)+5) = 1           ! RLE columns
                send_meta(smp(dst_p)+6) = -7          ! Debugging marker
             ELSE
                ! We are continuing the RLE run
                IF (blk_row_new .EQ. blk_row_prev_new) THEN
                   send_meta(smp(dst_p)+2) = send_meta(smp(dst_p)+2) + 1
                ELSEIF (blk_col_new .EQ. blk_col_prev_new) THEN
                   send_meta(smp(dst_p)+3) = send_meta(smp(dst_p)+3) + 1
                ELSE
                   CALL cp_assert(.FALSE., cp_failure_level, cp_internal_error,&
                        routineN, "logic error")
                ENDIF
             ENDIF
             IF (.NOT. tr) THEN
                IF (.NOT. PRESENT (alpha)) THEN
                   send_data(sdp(dst_p))&
                        = r2_dp(full_row-row_offset_old, full_col-col_offset_old)
                   ELSE
                      send_data(sdp(dst_p)) = alpha * &
                           r2_dp(full_row-row_offset_old, full_col-col_offset_old)

                   ENDIF
             ELSE
                IF (.NOT. PRESENT (alpha)) THEN
                   send_data(sdp(dst_p))&
                        = r2_dp(full_col-col_offset_old, full_row-row_offset_old)
                ELSE
                   send_data(sdp(dst_p)) = alpha * &
                        r2_dp(full_col-col_offset_old, full_row-row_offset_old)
                ENDIF
             ENDIF
             sdp(dst_p) = sdp(dst_p)+1
             send_count(2, dst_p) = send_count(2, dst_p) + 1
          ENDDO
          new_meta = .TRUE.
       ENDDO
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    ! Exchange the data and metadata structures.
    DBG 'E'
    CALL mp_alltoall(send_data(:), total_send_count(2,:), sd_disp(:)-1,&
         recv_data(:), total_recv_count(2,:), rd_disp(:)-1, mp_group)
    DBG 'F'
    DBG 'send_meta'
    WRITE(ifmt,'("(",I2,"(1X,I5))")') metalen
    IF (dbg) WRITE(*,ifmt) send_meta
    DBG 'send_data'
    IF (dbg) WRITE(*,'(10(1X,F5.2))') send_data
    CALL mp_alltoall(send_meta(:), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), metalen*total_recv_count(1,:), rm_disp(:)-1, mp_group)
    DBG 'recv_meta'
    IF (dbg) WRITE(*,ifmt) recv_meta
    DBG 'recv_data'
    IF (dbg) WRITE(*,'(10(1X,F5.2))') recv_data
    DBG 'G'
    ! Now fill in the data.
    CALL dbcsr_work_create(redist,&
            SUM(recv_count(1,:)),&
            SUM(recv_count(2,:)), work_mutable=.TRUE.)
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    !blk_p = 1
    !blk = 1
    blk_ps = 0
    blks = 0
    DO src_p = 0, numnodes-1
       data_offset_l = rd_disp(src_p)
       DO meta_l = 1, recv_count(1, src_p)
          stored_row_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1))
          stored_col_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+1)
          rel_full_row_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+2)
          rel_full_col_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+3)
          row_rle = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+4)
          col_rle = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+5)
          DBG 'blk: p',blks,src_p,'; row, col',&
               stored_row_new,stored_col_new,&
               '; blk_ps',blk_ps
          NULLIFY (r2_dp)
          CALL dbcsr_get_block_p(redist, stored_row_new, stored_col_new,&
               r2_dp, tr, found)
          valid_block = found
          IF (.NOT. found .AND. .NOT. my_keep_sparsity) THEN
             ALLOCATE (r2_dp(row_blk_size_new (stored_row_new),&
                  col_blk_size_new (stored_col_new)))
             r2_dp(:,:) = 0.0_dp ! Needed for nagging compilers.
             tr = .FALSE.
             !redist%m%wms(1)%lastblk = redist%m%wms(1)%lastblk+1
             !redist%m%wms(1)%datasize = redist%m%wms(1)%datasize&
             !     + row_rle + col+col_rle
             blks = blks + 1
             blk_ps = blk_ps + row_blk_size_new (stored_row_new) * &
                               col_blk_size_new (stored_col_new)
             valid_block = .TRUE.
          ENDIF
          IF (valid_block) THEN
             DO col = rel_full_col_new, rel_full_col_new + col_rle - 1
                DO row = rel_full_row_new, rel_full_row_new + row_rle - 1
                   IF (dbg) &
                        WRITE(*,*)'adding point',row,col,&
                        data_offset_l,recv_data(data_offset_l)
                   IF (.NOT. tr) THEN
                      r2_dp(row, col) = recv_data(data_offset_l)
                   ELSE
                      r2_dp(col, row) = recv_data(data_offset_l)
                   ENDIF
                   data_offset_l = data_offset_l + 1
                ENDDO
             ENDDO
          ENDIF
          IF (.NOT. found .AND. valid_block) THEN
             IF (dbg) WRITE(*,*)routineN//" Adding new block at",&
                  stored_row_new, stored_col_new
             CALL dbcsr_put_block(redist, stored_row_new, stored_col_new,&
                  r2_dp, tr)
             DEALLOCATE (r2_dp)
          ELSE
             IF (.NOT. my_keep_sparsity .AND. dbg) &
                  WRITE(*,*)routineN//" Reusing block at",&
                  stored_row_new, stored_col_new
          ENDIF
       ENDDO
    ENDDO
    !
    IF (dbg) THEN
       WRITE(*,*)routineN//" Declared blocks=",redist%m%wms(1)%lastblk,&
            "actual=",blks
       WRITE(*,*)routineN//" Declared data size=",redist%m%wms(1)%datasize,&
            "actual=",blk_ps
    ENDIF
    DBG 'Finalizing redistribution'
    CALL dbcsr_finalize(redist)
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_data)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_data)
    DEALLOCATE(send_meta)
    DBG 'Done redistributing'
    CALL timestop(timing_handle)
  CONTAINS
    SUBROUTINE fill_map (map, nfull, offsets, nblks, last_size)
    INTEGER, INTENT(in)                      :: nfull
    INTEGER, DIMENSION(nfull), INTENT(out)   :: map
    INTEGER, INTENT(in)                      :: nblks
    INTEGER, DIMENSION(nblks), INTENT(in)    :: offsets
    INTEGER, INTENT(in)                      :: last_size

    INTEGER                                  :: blk

!

      FORALL (blk = 1:nblks-1)
         map(offsets(blk):offsets(blk+1)-1) = blk
      END FORALL
      map(offsets(nblks):offsets(nblks)+last_size-1) = nblks
    END SUBROUTINE fill_map
  END SUBROUTINE dbcsr_complete_redistribute

! *****************************************************************************
!> \brief Redistributes a DBCSR matrix.
!>
!>        The new distribution should have compatible row and column blocks.
!> \param[in] matrix          matrix to redistribute
!> \param[in,out] redist      redistributed matrix, which should already be
!>                            created
! *****************************************************************************
  SUBROUTINE dbcsr_redistribute(matrix, redist)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_obj), INTENT(INOUT)           :: redist

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_redistribute', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 2

    CHARACTER(len=80)                        :: ifmt
    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: c_dp
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: c_sp
    INTEGER :: blk, blk_ps, blks, col, dst_p, meta_l, mp_group, numnodes, &
      nze, row, row_prev, src_p, stored_col_new, stored_row_new, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rd_disp, recv_meta, rm_disp, &
                                                sd_disp, sdp, send_meta, &
                                                sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_count, send_count, &
                                                total_recv_count, &
                                                total_send_count
    INTEGER, DIMENSION(:), POINTER :: col_blk_size_new, col_blk_size_old, &
      col_dist_new, col_dist_old, row_blk_size_new, row_blk_size_old, &
      row_dist_new, row_dist_old
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    LOGICAL                                  :: tr
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: send_data
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: recv_data
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_dp
    REAL(KIND=sp), DIMENSION(:), POINTER     :: r_sp
    TYPE(dbcsr_distribution_obj)             :: dist_new, dist_old
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj_new

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    IF (dbg) CALL dbcsr_print (matrix)
    !call dbcsr_print_dist (matrix%m%dist)
    !call dbcsr_print_dist (redist%m%dist)
    CALL cp_assert (dbcsr_valid_index(matrix%m),&
         cp_fatal_level, cp_caller_error,&
         routineN, "Input not valid.")
    CALL cp_assert (dbcsr_is_initialized(redist%m),&
         cp_fatal_level, cp_caller_error,&
         routineN, "Output matrix not initialized.")
    CALL cp_assert (matrix%m%replication_type .EQ. dbcsr_repl_none,&
         cp_warning_level, cp_caller_error, routineN,&
         "Can not redistribute replicated matrix.")
    SELECT CASE (matrix%m%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (matrix%m%data_area, r_dp)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (matrix%m%data_area, r_sp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (matrix%m%data_area, c_dp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (matrix%m%data_area, c_sp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    END SELECT
    ! Get row and column start and end positions
    ! Old matrix
    dist_old = dbcsr_distribution (matrix)
    row_blk_size_old => array_data (matrix%m%row_blk_size)
    col_blk_size_old => array_data (matrix%m%col_blk_size)
    row_dist_old => array_data (dbcsr_distribution_row_dist (dist_old))
    col_dist_old => array_data (dbcsr_distribution_col_dist (dist_old))
    ! New matrix
    dist_new = dbcsr_distribution (redist)
    row_blk_size_new => array_data (redist%m%row_blk_size)
    col_blk_size_new => array_data (redist%m%col_blk_size)
    row_dist_new => array_data (dbcsr_distribution_row_dist (dist_new))
    col_dist_new => array_data (dbcsr_distribution_col_dist (dist_new))
    !
    mp_obj_new = dbcsr_distribution_mp (dist_new)
    pgrid => dbcsr_mp_pgrid (mp_obj_new)
    numnodes = dbcsr_mp_numnodes (mp_obj_new)
    mp_group = dbcsr_mp_group (mp_obj_new)
    !
    CALL cp_assert(MAXVAL(row_dist_new).LE.UBOUND(pgrid,1), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows')
    IF (dbg) THEN
       CALL cp_assert(MAXVAL(row_dist_new).EQ.UBOUND(pgrid,1),&
            cp_warning_level, cp_wrong_args_error, routineN,&
            'Range of row distribution not equal to processor rows')
    ENDIF
    CALL cp_assert(MAXVAL(col_dist_new).LE.UBOUND(pgrid,2), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols')
    IF (dbg) THEN
       CALL cp_assert(MAXVAL(col_dist_new).EQ.UBOUND(pgrid,2),&
            cp_warning_level, cp_wrong_args_error, routineN,&
            'Range of col distribution not equal to processor cols')
    ENDIF
    ALLOCATE (send_count(2, 0:numnodes-1))
    ALLOCATE (recv_count(2, 0:numnodes-1))
    ALLOCATE (total_send_count(2, 0:numnodes-1))
    ALLOCATE (total_recv_count(2, 0:numnodes-1))
    ALLOCATE (sdp(0:numnodes-1))
    ALLOCATE (sd_disp(0:numnodes-1))
    ALLOCATE (smp(0:numnodes-1))
    ALLOCATE (sm_disp(0:numnodes-1))
    ALLOCATE (rd_disp(0:numnodes-1))
    ALLOCATE (rm_disp(0:numnodes-1))
    ! Count initial sizes for sending.
    !
    DBG 'A'
    send_count(:,:) = 0
    CALL dbcsr_iterator_start (iter, matrix)
    row_prev = 0
    dst_p = -1
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, r_dp, tr)
       CALL dbcsr_get_stored_coordinates (redist%m,&
            row, col, tr, dst_p)
       IF (dbg) WRITE(*,*)routinen//" Sending",row,col,"to",dst_p
       nze = SIZE(r_dp)
       send_count(1, dst_p) = send_count(1, dst_p) + 1
       send_count(2, dst_p) = send_count(2, dst_p) + SIZE(r_dp)
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    DBG 'B'
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    DBG 'C'
    DBG 'send counts',send_count
    DBG 'recv counts',recv_count
    ! Allocate data structures needed for data exchange.
    ALLOCATE (recv_data(SUM(recv_count(2, :))))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :))))
    ALLOCATE (send_data(SUM(send_count(2, :))))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :))))
    DBG 'send data len',SIZE (send_data)
    DBG 'send meta len',SIZE (send_meta)
    DBG 'recv data len',SIZE (recv_data)
    DBG 'recv meta len',SIZE (recv_meta)
    ! Fill in the meta data structures and copy the data.
    DO dst_p = 0, numnodes-1
       total_send_count(1, dst_p) = send_count (1, dst_p)
       total_send_count(2, dst_p) = send_count (2, dst_p)
       total_recv_count(1, dst_p) = recv_count (1, dst_p)
       total_recv_count(2, dst_p) = recv_count (2, dst_p)
    ENDDO
    sd_disp = -1 ; sm_disp = -1
    rd_disp = -1 ; rm_disp = -1
    sd_disp(0) = 1 ; sm_disp(0) = 1
    rd_disp(0) = 1 ; rm_disp(0) = 1
    DO dst_p = 1, numnodes-1
       sm_disp(dst_p) = sm_disp(dst_p-1)&
                        + metalen*total_send_count(1, dst_p-1)
       sd_disp(dst_p) = sd_disp(dst_p-1)&
                        + total_send_count(2, dst_p-1)
       rm_disp(dst_p) = rm_disp(dst_p-1)&
                        + metalen*total_recv_count(1, dst_p-1)
       rd_disp(dst_p) = rd_disp(dst_p-1)&
                        + total_recv_count(2, dst_p-1)
    ENDDO
    sdp(:) = sd_disp     ! sdp points to the the next place to store
                         ! data. It is postincremented.
    smp(:) = sm_disp - metalen  ! But smp points to the "working" data, not
                                ! the next. It is pre-incremented, so we must
                                ! first rewind it.
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    DBG 'D'
    IF (dbg) THEN
       WRITE(*,*)'row_blk_size_new',row_blk_size_new
    ENDIF
    CALL dbcsr_iterator_start (iter, matrix)
    dst_p = -1
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, r_dp, tr, blk)
       !IF (tr) WRITE(*,*)"block at",row,col," is transposed"
       CALL dbcsr_get_stored_coordinates (redist%m,&
            row, col, tr, dst_p)
       smp(dst_p) = smp(dst_p) + metalen
       IF (tr) THEN
          send_meta(smp(dst_p)) = -row
       ELSE
          send_meta(smp(dst_p)) = row
       ENDIF
       send_meta(smp(dst_p)+1) = col ! new blocked column
       send_data(sdp(dst_p):sdp(dst_p)+SIZE(r_dp)-1) &
            = r_dp(:)
       sdp(dst_p) = sdp(dst_p) + SIZE(r_dp)
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    ! Exchange the data and metadata structures.
    DBG 'E'
    CALL mp_alltoall(send_data(:), total_send_count(2,:), sd_disp(:)-1,&
         recv_data(:), total_recv_count(2,:), rd_disp(:)-1, mp_group)
    DBG 'F'
    DBG 'send_meta'
    IF (dbg) WRITE(ifmt,'("(",I2,"(1X,I5))")') metalen
    IF (dbg) WRITE(*,ifmt) send_meta
    DBG 'send_data'
    IF (dbg) WRITE(*,'(10(1X,F5.2))') send_data
    CALL mp_alltoall(send_meta(:), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), metalen*total_recv_count(1,:), rm_disp(:)-1, mp_group)
    DBG 'recv_meta'
    IF (dbg) WRITE(*,ifmt) recv_meta
    DBG 'recv_data'
    IF (dbg) WRITE(*,'(10(1X,F5.2))') recv_data
    DBG 'G'
    ! Now fill in the data.
    CALL dbcsr_work_create(redist,&
            SUM(recv_count(1,:)),&
            SUM(recv_count(2,:)), work_mutable=.FALSE.)
    !
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    blk_ps = 1
    blks = 0
    DO src_p = 0, numnodes-1
       !data_offset_l = rd_disp(src_p)
       DO meta_l = 1, recv_count(1, src_p)
          row = recv_meta(rm_disp(src_p)+metalen*(meta_l-1))
          tr = row .LT. 0
          stored_row_new = ABS(row)
          stored_col_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+1)
          nze = row_blk_size_new(stored_row_new) * col_blk_size_new(stored_col_new)
          r_dp => recv_data(blk_ps:blk_ps+nze-1)
          CALL dbcsr_put_block(redist, stored_row_new, stored_col_new, r_dp, tr)
          blk_ps = blk_ps + nze
          blks = blks + 1
       ENDDO
    ENDDO
    IF (dbg) WRITE(*,*)routineN//" blk_p", redist%m%wms(1)%blk_p
    !
    IF (dbg) THEN
       WRITE(*,*)routineN//" Declared blocks=",redist%m%wms(1)%lastblk,&
            "actual=",blks
       WRITE(*,*)routineN//" Declared data size=",redist%m%wms(1)%datasize,&
            "actual=",blk_ps
    ENDIF
    DBG 'Finalizing redistribution'
    CALL dbcsr_finalize(redist)
    IF (dbg) CALL dbcsr_print (redist)
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_data)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_data)
    DEALLOCATE(send_meta)
    DBG 'Done redistributing'
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_redistribute


! *****************************************************************************
!> \brief Transposes a DBCSR matrix.
!> \param[out] transposed     transposed DBCSR matrix
!> \param[in] normal          input DBCSR matrix
!> \param[in] shallow_data_copy         (optional) only shallow data_copy;
!>                                      default is no; if set, the
!>                                      transpose_data option is ignored
!> \param[in] transpose_data  (optional) transpose data blocks, default is True
!> \param[in] transpose_distribution    (optional) transpose the distribution
!>                                      from the input matrix, default is True
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_new_transposed (transposed, normal, shallow_data_copy,&
       transpose_data, transpose_distribution, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: transposed
    TYPE(dbcsr_obj), INTENT(IN)              :: normal
    LOGICAL, INTENT(IN), OPTIONAL            :: shallow_data_copy, &
                                                transpose_data, &
                                                transpose_distribution
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_new_transposed', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: new_type
    INTEGER                                  :: stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: blk_p
    LOGICAL                                  :: shallow, tr_blocks, tr_dist
    TYPE(dbcsr_distribution_obj)             :: new_dist
    TYPE(dbcsr_obj)                          :: t2

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL cp_assert (dbcsr_valid_index(normal), cp_fatal_level, cp_caller_error,&
         routineN, "Matrix does not exist.")
    ! Internalize options
    shallow = .FALSE.
    IF (PRESENT (shallow_data_copy)) shallow = shallow_data_copy
    tr_blocks = .TRUE.
    IF (PRESENT (transpose_data)) tr_blocks = transpose_data
    tr_dist = .TRUE.
    IF (PRESENT (transpose_distribution)) tr_dist = transpose_distribution
    ! Prepare the distribution for the transposed matrix
    IF (tr_dist) THEN
       !write(*,*)' transposing mp dist...'
       CALL dbcsr_transpose_distribution (new_dist, normal%m%dist)
    ELSE
       !write(*,*)' transposing dist not mp...'
       CALL dbcsr_transpose_dims (new_dist, normal%m%dist)
    ENDIF
    ! Create the transposed matrix
    !CALL dbcsr_transposed_type (new_type, normal%m%matrix_type)
    new_type = normal%m%matrix_type
    CALL dbcsr_create (transposed, "transposed "//normal%m%name, new_dist,&
         new_type, normal%m%col_blk_size, normal%m%row_blk_size,&
         data_type=normal%m%data_type, special=normal%m%special_memory)
    CALL dbcsr_distribution_release (new_dist)
    ! Reserve the space for the new indices.
    CALL dbcsr_addto_index_array (transposed%m, dbcsr_slot_row_p,&
         reservation=transposed%m%nblkrows_total+1, extra=transposed%m%nblks*2)
    CALL dbcsr_addto_index_array (transposed%m, dbcsr_slot_col_i,&
         reservation=normal%m%nblks)
    CALL dbcsr_addto_index_array (transposed%m, dbcsr_slot_blk_p,&
         reservation=normal%m%nblks)
    CALL dbcsr_repoint_index (transposed%m)
    IF (.NOT. shallow) THEN
       CALL dbcsr_data_ensure_size (transposed%m%data_area, normal%m%nze,&
            nocopy=.TRUE.)
    ENDIF
    !
    transposed%m%nblks = normal%m%nblks
    transposed%m%nze = normal%m%nze
    transposed%m%index(dbcsr_slot_nblks) = normal%m%nblks
    transposed%m%index(dbcsr_slot_nze) = normal%m%nze
    ! Transpose the local index.
    DBG 'Original size', normal%m%nblks, normal%m%nze
    DBG 'meta', transposed%m%nblkrows_total+1, normal%m%nblks
    DBG 'size', SIZE (transposed%m%row_p), SIZE (transposed%m%col_i)
    ALLOCATE (blk_p (normal%m%nblks), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineP, "blk_p", normal%m%nblks)
    CALL transpose_index_local (transposed%m%row_p, transposed%m%col_i,&
         normal%m%row_p, normal%m%col_i, blk_p, normal%m%blk_p)
    IF (dbg) THEN
       WRITE(*,*)'orig. row_p',normal%m%row_p
       WRITE(*,*)'orig. col_i',normal%m%col_i
       WRITE(*,*)'orig. blk_p',normal%m%blk_p
       WRITE(*,*)'new . row_p',transposed%m%row_p
       WRITE(*,*)'new . col_i',transposed%m%col_i
       WRITE(*,*)'new . blk_p',blk_p!transposed%m%blk_p
    ENDIF
    ! Copy the data
    IF (shallow) THEN
       CALL dbcsr_data_release (transposed%m%data_area)
       transposed%m%data_area = normal%m%data_area
       CALL dbcsr_data_hold (transposed%m%data_area)
       transposed%m%blk_p(1:transposed%m%nblks) =&
            -blk_p(1:transposed%m%nblks)
    ELSE
       CALL dbcsr_copy_sort_data (transposed%m%blk_p, blk_p, transposed%m%row_p,&
            transposed%m%col_i, array_data (transposed%m%row_blk_size),&
            array_data (transposed%m%col_blk_size),&
            transposed%m%data_area, normal%m%data_area,&
            mark_transposed=.not.tr_blocks,&
            transpose_blocks=tr_blocks)
    ENDIF
    transposed%m%valid = .TRUE.
    !CALL dbcsr_copy_sort_data (transposed%m%blk_p, blk_p, transposed%m%row_p,&
    !     transposed%m%col_i, array_data (transposed%m%row_blk_size),&
    !     array_data (transposed%m%col_blk_size),&
    !     transposed%m%data_area, normal%m%data_area,&
    !     transpose_blocks=.TRUE.)
    !
1315 FORMAT (I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5)
    IF (dbg) THEN
       WRITE(*,*)'new FINAL index'
       WRITE(*,1315)transposed%m%row_p
       WRITE(*,1315)transposed%m%col_i
       WRITE(*,1315)transposed%m%blk_p
    ENDIF
    !
    DEALLOCATE (blk_p)
    !
    ! Redistribute if mpgrid was not transposed
    IF (.NOT. tr_dist) THEN
       !write (*,*)routineN//" redistributing"
       CALL dbcsr_init (t2)
       CALL dbcsr_create (t2, transposed%m%name, transposed%m%dist,&
            'N', transposed%m%row_blk_size, transposed%m%col_blk_size,&
            error=error)
       CALL dbcsr_redistribute (transposed, t2)
       CALL dbcsr_release (transposed)
       transposed = t2
    ENDIF
  END SUBROUTINE dbcsr_new_transposed

! *****************************************************************************
!> \brief Duplicates data in symmetric matrix to make it normal (w.r.t. data
!>        structure
!> \param[in] sm              input symmetric matrix
!> \param[out] desm           desymmetrized matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_desymmetrize_deep(sm, desm, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: sm
    TYPE(dbcsr_obj), INTENT(OUT)             :: desm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_desymmetrize_deep', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 2

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: c_dp
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: c_sp
    INTEGER :: blk, blk_l, blk_p, blk_ps, blks, col, dst_p, mp_group, &
      nsymmetries, numproc, nze, pcol, prow, row, src_p, stored_col, &
      stored_row, symmetry_i, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rd_disp, recv_meta, rm_disp, &
                                                sd_disp, sdp, send_meta, &
                                                sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_count, send_count, &
                                                total_recv_count, &
                                                total_send_count
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, col_dist, &
                                                row_blk_size, row_dist
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    LOGICAL                                  :: tr
    REAL(KIND=dp)                            :: tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: send_data
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: recv_data
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_dp
    REAL(KIND=sp), DIMENSION(:), POINTER     :: r_sp
    TYPE(dbcsr_distribution_obj)             :: target_dist
    TYPE(dbcsr_mp_obj)                       :: mp_obj

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    tstart = m_walltime ()
    CALL cp_assert (dbcsr_valid_index (sm),&
         cp_fatal_level, cp_caller_error, routineN,&
         "Matrix not initialized.")
    nsymmetries = 1
    IF (sm%m%matrix_type .EQ. 'S' .OR. sm%m%matrix_type .EQ. 'A') THEN
       nsymmetries = 2
    ENDIF
    SELECT CASE (sm%m%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (sm%m%data_area, r_dp)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (sm%m%data_area, r_sp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (sm%m%data_area, c_dp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (sm%m%data_area, c_sp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    END SELECT
    row_blk_size => array_data (sm%m%row_blk_size)
    col_blk_size => array_data (sm%m%col_blk_size)
    target_dist = sm%m%dist
    row_dist => array_data (dbcsr_distribution_row_dist (target_dist))
    col_dist => array_data (dbcsr_distribution_col_dist (target_dist))
    mp_obj = dbcsr_distribution_mp (target_dist)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    numproc = dbcsr_mp_numnodes (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    DBGV 'row_dist', row_dist
    DBGV 'col_dist', col_dist
    DBG 'row, col maxval:', MAXVAL(row_dist), MAXVAL(col_dist)
    DBG 'blacs2mpi',blacs2mpi
    DBG 'blacs2mpi LB',LBOUND(blacs2mpi)
    DBG 'blacs2mpi UB',UBOUND(blacs2mpi)
    IF (sm%m%matrix_type .EQ. 'S' .OR. sm%m%matrix_type .EQ. 'A') THEN
       CALL cp_assert(SIZE(row_dist) .EQ. SIZE(col_dist), cp_warning_level,&
            cp_wrong_args_error, routineN,&
            'Unequal row and column distributions for symmetric matrix.', error)
    ENDIF
    CALL cp_assert(MAXVAL(row_dist).LE.UBOUND(blacs2mpi,1), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows', error)
    IF (dbg) &
         CALL cp_assert(MAXVAL(row_dist).EQ.UBOUND(blacs2mpi,1), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Range of row distribution not equal to processor rows', error)
    CALL cp_assert(MAXVAL(col_dist).LE.UBOUND(blacs2mpi,2), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols', error)
    IF (dbg) &
         CALL cp_assert(MAXVAL(col_dist).EQ.UBOUND(blacs2mpi,2), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Range of col distribution not equal to processor cols', error)
    ALLOCATE (send_count(2, 0:numproc-1))
    ALLOCATE (recv_count(2, 0:numproc-1))
    ALLOCATE (total_send_count(2, 0:numproc-1))
    ALLOCATE (total_recv_count(2, 0:numproc-1))
    ALLOCATE (sdp(0:numproc-1))
    ALLOCATE (sd_disp(0:numproc-1))
    ALLOCATE (smp(0:numproc-1))
    ALLOCATE (sm_disp(0:numproc-1))
    ALLOCATE (rd_disp(0:numproc-1))
    ALLOCATE (rm_disp(0:numproc-1))
    CALL dbcsr_create(desm, 'desym '//sm%m%name, sm%m%dist,&
         dbcsr_type_normal, sm%m%row_blk_size, sm%m%col_blk_size,&
         0,0, dbcsr_get_data_type(sm), error=error, special=sm%m%special_memory)
    IF (sm%m%matrix_type .EQ. dbcsr_type_antisymmetric) THEN
       desm%m%sym_negation = .TRUE.
    ENDIF
    ! Count initial sizes for sending.
    DBG 'A'
    send_count(:,:) = 0
    DO row = 1, sm%m%nblkrows_total
       DO blk = sm%m%row_p(row)+1, sm%m%row_p(row+1)
          col = sm%m%col_i(blk)
          DO symmetry_i = 1, nsymmetries
             IF (symmetry_i .EQ. 1) THEN
                stored_row = row ; stored_col = col
             ELSE
                IF (row .EQ. col) CYCLE
                stored_row = col ; stored_col = row
             ENDIF
             ! Where do we send this block?
             prow = row_dist(stored_row)
             pcol = col_dist(stored_col)
             dst_p = blacs2mpi(prow, pcol)
             nze = row_blk_size(stored_row) * col_blk_size(stored_col)
             send_count(1, dst_p) = send_count(1, dst_p) + 1
             send_count(2, dst_p) = send_count(2, dst_p) + nze
          ENDDO ! symmetry_i
       ENDDO ! col_i
    ENDDO ! row
    DBG 'B'
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    DBG 'C'
    DBG 'send counts',send_count
    DBG 'recv counts',recv_count
    ! Allocate data structures needed for data exchange.
    ALLOCATE (recv_data(SUM(recv_count(2, :))))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :))))
    ALLOCATE (send_data(SUM(send_count(2, :))))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :))))
    DBG 'send data len',SIZE (send_data)
    DBG 'send meta len',SIZE (send_meta)
    DBG 'recv data len',SIZE (recv_data)
    DBG 'recv meta len',SIZE (recv_meta)
    ! Fill in the meta data structures and copy the data.
    DO dst_p = 0, numproc-1
       total_send_count(1, dst_p) = send_count (1, dst_p)
       total_send_count(2, dst_p) = send_count (2, dst_p)
       total_recv_count(1, dst_p) = recv_count (1, dst_p)
       total_recv_count(2, dst_p) = recv_count (2, dst_p)
    ENDDO
    sd_disp = -1 ; sm_disp = -1
    rd_disp = -1 ; rm_disp = -1
    sd_disp(0) = 1 ; sm_disp(0) = 1
    rd_disp(0) = 1 ; rm_disp(0) = 1
    DO dst_p = 1, numproc-1
       sm_disp(dst_p) = sm_disp(dst_p-1)&
                        + metalen*total_send_count(1, dst_p-1)
       sd_disp(dst_p) = sd_disp(dst_p-1)&
                        + total_send_count(2, dst_p-1)
       rm_disp(dst_p) = rm_disp(dst_p-1)&
                        + metalen*total_recv_count(1, dst_p-1)
       rd_disp(dst_p) = rd_disp(dst_p-1)&
                        + total_recv_count(2, dst_p-1)
    ENDDO
    !DO dst_p = 0, numproc-1
    !      sm_disp (dst_p) = sm_disp (dst_p) + metalen*send_count(1, dst_p)
    !      sd_disp (dst_p) = sd_disp (dst_p) + send_count(2, dst_p)
    !      rm_disp (dst_p) = rm_disp (dst_p) + metalen*recv_count(1, dst_p)
    !      rd_disp (dst_p) = rd_disp (dst_p) + recv_count(2, dst_p)
    !ENDDO
    sdp(:) = sd_disp
    smp(:) = sm_disp
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    DBG 'D'
    DO row = 1, sm%m%nblkrows_total
       DO blk = sm%m%row_p(row)+1, sm%m%row_p(row+1)
          col = sm%m%col_i(blk)
          blk_p = sm%m%blk_p(blk)
          DO symmetry_i = 1, nsymmetries
             IF (symmetry_i .EQ. 1) THEN
                stored_row = row ; stored_col = col; tr = .FALSE.
             ELSE                                                
                IF (row .EQ. col) CYCLE
                stored_row = col ; stored_col = row; tr = .TRUE.
             ENDIF
             ! Where do we send this block?
             prow = row_dist(stored_row)
             pcol = col_dist(stored_col)
             dst_p = blacs2mpi(prow, pcol)
             nze = row_blk_size(stored_row) * col_blk_size(stored_col)
             send_meta(smp(dst_p)) = stored_row
             send_meta(smp(dst_p)+1) = stored_col
             DBG stored_row, stored_col, blk_p, tr, nze
             IF (.NOT. tr) THEN
                send_data(sdp(dst_p):sdp(dst_p)+nze-1) =&
                     r_dp(blk_p:blk_p+nze-1)
                !send_data(sdp(dst_p):sdp(dst_p)+nze-1) =&
                !     sm%m%data(blk_p:blk_p+nze-1)
             ELSE
                send_meta(smp(dst_p)) = -stored_row
                send_data(sdp(dst_p):sdp(dst_p)+nze-1) =&
                     r_dp(blk_p:blk_p+nze-1)
                !send_data(sdp(dst_p):sdp(dst_p)+nze-1) =&
                !     sm%m%data(blk_p:blk_p+nze-1)
             ENDIF
             smp(dst_p) = smp(dst_p) + metalen
             sdp(dst_p) = sdp(dst_p) + nze
          ENDDO ! symmetry_i
       ENDDO ! col_i
    ENDDO ! row
    ! Exchange the data and metadata structures.
    DBG 'E'
    CALL mp_alltoall(send_data(:), total_send_count(2,:), sd_disp(:)-1,&
         recv_data(:), total_recv_count(2,:), rd_disp(:)-1, mp_group)
    DBG 'F'
    DBG 'send_meta',send_meta
    CALL mp_alltoall(send_meta(:), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), metalen*total_recv_count(1,:), rm_disp(:)-1, mp_group)
    DBG 'recv_meta',recv_meta
    DBG 'G'
    ! Now fill in the data.
    CALL dbcsr_work_create(desm,&
            SUM(recv_count(1,:)),&
            SUM(recv_count(2,:)), error=error)
    !DBG 'send disps data',sd_disp
    !DBG 'send disps meta',sm_disp
    !DBG 'recv disps data',rd_disp
    !DBG 'recv disps meta',rm_disp
    !blk_p = 1
    !blk = 1
    blk_ps = 1
    blks = 1
    DO src_p = 0, numproc-1
       !wm%data(blk_p : blk_p+recv_count(2,src_p)-1) =&
       !     recv_data(rd_disp(src_p) :&
       !               rd_disp(src_p)+recv_count(2,src_p)-1)
       nze = recv_count(2, src_p)
       CALL dbcsr_data_set (desm%m%wms(1)%data_area, blk_ps, nze,&
            recv_data, rd_disp(src_p))
       !desm%m%w%data(blk_ps:blk_ps+nze-1) =&
       !     recv_data(rd_disp(src_p):rd_disp(src_p)+nze-1)
       DBG 'processing from',src_p,'counts',send_count(1, src_p)
       DO blk_l = 1, recv_count(1, src_p)
          stored_row = recv_meta(rm_disp(src_p)+metalen*(blk_l-1))
          stored_col = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+1)
          !DBG 'blk: p',blk,src_p,'; row, col',stored_row,stored_col,&
          !     '; blk_p',blk_p,'; row_rep, imgdist',row_rep,&
          !     row_rep_dist (ABS (stored_row))
          DBG 'blk: p',blks,src_p,'; row, col',&
               stored_row,stored_col,&
               '; blk_ps',blk_ps
          !row_rep = row_rep_dist (ABS(stored_row))
          desm%m%wms(1)%row_i(blks) = ABS(stored_row)
          desm%m%wms(1)%col_i(blks) = stored_col
          desm%m%wms(1)%blk_p(blks) = SIGN(blk_ps, stored_row)
          nze = row_blk_size(ABS(stored_row))&
               * col_blk_size(stored_col)
          !blk_p = blk_p + nze
          !blk = blk + 1
          blk_ps = blk_ps + nze
          blks = blks + 1
       ENDDO
    ENDDO
    !
    desm%m%wms(1)%lastblk = blks - 1
    desm%m%wms(1)%datasize = blk_ps - 1
    DBG 'Finalizing normalization'
    CALL dbcsr_finalize(desm)
    tstop = m_walltime ()
    DBG 'time:',tstop-tstart
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_data)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_data)
    DEALLOCATE(send_meta)
    DBG 'desym name',desm%m%name
    DBG 'Done desymmetrizing!'
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_desymmetrize_deep


! *****************************************************************************
!> \brief Distributes a matrix that is currently replicated.
!> \param[in,out] matrix      matrix to replicate
!> \param[in] fast            change just the index, don't touch the data
! *****************************************************************************
  SUBROUTINE dbcsr_distribute(matrix, fast)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    LOGICAL, INTENT(in), OPTIONAL            :: fast

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_distribute', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 2

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: c_dp
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: c_sp
    INTEGER                                  :: blk, col, mp_group, mynode, &
                                                numnodes, p, row
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, col_dist, &
                                                row_blk_size, row_dist, &
                                                tmp_index
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    LOGICAL                                  :: mini, tr
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_dp
    REAL(KIND=sp), DIMENSION(:), POINTER     :: r_sp
    TYPE(dbcsr_data_obj)                     :: tmp_data
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_obj)                          :: distributed

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL cp_assert (dbcsr_valid_index(matrix%m),&
         cp_fatal_level, cp_caller_error,&
         routineN, "Matrix not initialized.")
    CALL cp_assert (matrix%m%replication_type .NE. dbcsr_repl_none,&
         cp_warning_level, cp_caller_error, routineN,&
         "Distributing a non-replicated matrix makes no sense.")
    IF (PRESENT (fast)) THEN
       mini=fast
    ELSE
       mini=.FALSE.
    ENDIF
    SELECT CASE (matrix%m%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (matrix%m%data_area, r_dp)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (matrix%m%data_area, r_sp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (matrix%m%data_area, c_dp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (matrix%m%data_area, c_sp)
       CALL cp_unimplemented_error(routineN, "Only double precision real.")
    END SELECT
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    dist = dbcsr_distribution (matrix)
    row_dist => array_data (dbcsr_distribution_row_dist (dist))
    col_dist => array_data (dbcsr_distribution_col_dist (dist))
    mp_obj = dbcsr_distribution_mp (dist)
    numnodes = dbcsr_mp_numnodes (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    mynode = dbcsr_mp_mynode (dbcsr_distribution_mp (dist))
    pgrid => dbcsr_mp_pgrid (dbcsr_distribution_mp (dist))

    !
    IF (mini) THEN
       ! We just mark the blocks as deleted.
       CALL dbcsr_iterator_start (iter, matrix)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, row, col, r_dp, tr, blk)
          tr = .FALSE.
          CALL dbcsr_get_stored_coordinates (matrix%m, row, col, tr, p)
          IF (mynode .EQ. p) THEN
             matrix%m%blk_p(blk) = 0
          ENDIF
       ENDDO
       CALL dbcsr_iterator_stop (iter)
       matrix%m%replication_type = dbcsr_repl_none
    ELSE
       CALL dbcsr_init(distributed)
       CALL dbcsr_create(distributed, 'Distributed '//matrix%m%name,&
            matrix%m%dist,&
            'N', matrix%m%row_blk_size, matrix%m%col_blk_size,&
            0, 0, matrix%m%data_type, special=matrix%m%special_memory)
       CALL dbcsr_work_create(distributed)
       distributed%m%replication_type = dbcsr_repl_none
       CALL dbcsr_iterator_start (iter, matrix)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, row, col, r_dp, tr, blk)
          tr = .FALSE.
          CALL dbcsr_get_stored_coordinates (matrix%m, row, col, tr, p)
          IF (mynode .EQ. p) THEN
             CALL dbcsr_put_block (distributed, row, col, r_dp, tr)
          ENDIF
       ENDDO
       CALL dbcsr_iterator_stop (iter)
       CALL dbcsr_finalize (distributed)
       ! Now replace the data and index
       tmp_data = matrix%m%data_area
       matrix%m%data_area = distributed%m%data_area
       distributed%m%data_area = tmp_data
       tmp_index => matrix%m%index
       matrix%m%index => distributed%m%index
       distributed%m%index => tmp_index
       CALL dbcsr_repoint_index (matrix%m)
       matrix%m%nze = distributed%m%nze
       matrix%m%nblks = distributed%m%nblks
       CALL dbcsr_release (distributed)
    ENDIF
    DBG 'Done replicating'
  END SUBROUTINE dbcsr_distribute


! *****************************************************************************
!> \brief Sums blocks in a replicated dbcsr matrix.
!> \param[in,out] matrix      dbcsr matrix to operate on
! *****************************************************************************
  SUBROUTINE dbcsr_sum_replicated (matrix)
    TYPE(dbcsr_obj), INTENT(inout)           :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sum_replicated', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: comm, max_data_size, &
                                                max_index_size, mynode, &
                                                numnodes, p
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: send_index
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: all_sizes
    INTEGER, DIMENSION(2)                    :: my_sizes
    INTEGER, DIMENSION(:), POINTER           :: recv_index
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: send_dp
    REAL(kind=dp), DIMENSION(:), POINTER     :: recv_dp
    TYPE(cp_error_type)                      :: error
    TYPE(dbcsr_mp_obj)                       :: mp
    TYPE(dbcsr_obj)                          :: recv_matrix

!   ---------------------------------------------------------------------------

    mp = dbcsr_distribution_mp (dbcsr_distribution (matrix))
    comm = dbcsr_mp_group (mp)
    numnodes = dbcsr_mp_numnodes (mp)
    mynode = dbcsr_mp_mynode (mp)
    ! Count data sizes
    my_sizes(1) = SIZE (matrix%m%index)
    recv_dp => dbcsr_get_data_p (matrix%m%data_area, 0.0_dp)
    my_sizes(2) = SIZE (recv_dp)
    ! Allocate receive buffers
    ALLOCATE (all_sizes (2, 0:numnodes-1))
    CALL mp_allgather (my_sizes, all_sizes, comm)
    max_index_size = MAXVAL (all_sizes(1,:))
    max_data_size = MAXVAL (all_sizes(2,:))
    CALL dbcsr_init (recv_matrix)
    CALL cp_error_init (error)
    CALL dbcsr_create(recv_matrix, "Temporary receiving buffer",&
         dbcsr_distribution(matrix), dbcsr_get_matrix_type(matrix),&
         dbcsr_row_block_sizes(matrix), dbcsr_col_block_sizes(matrix),&
         0, 0, dbcsr_get_data_type(matrix), error=error)
    CALL ensure_array_size (recv_matrix%m%index, ub=max_index_size)
    CALL dbcsr_data_ensure_size (recv_matrix%m%data_area, max_data_size)
    ! Setup data to be sent.
    ALLOCATE (send_dp(my_sizes(2)))
    send_dp(:) = recv_dp(:)
    ALLOCATE (send_index(my_sizes(1)))
    send_index(:) = matrix%m%index(:)
    ! Broadcasts a matrix and every processor sums it. Inefficient and
    ! slow, but it gets the job done with little memory use.
    DO p = 0, numnodes-1
       IF (p .EQ. mynode) THEN
          CALL mp_bcast (send_index, p, comm)
          CALL mp_bcast (send_dp, p, comm)
       ELSE
          recv_index => recv_matrix%m%index(1:all_sizes(1,p))
          recv_dp => pointer_view ( dbcsr_get_data_p (&
               recv_matrix%m%data_area, 0.0_dp), 1, all_sizes(2,p))
          CALL mp_bcast (recv_index, p, comm)
          CALL mp_bcast (recv_dp, p, comm)
          CALL dbcsr_repoint_index(recv_matrix%m)
          CALL dbcsr_add (matrix, recv_matrix, 1.0_dp, 1.0_dp, error=error)
        ENDIF
    ENDDO
    CALL dbcsr_destroy (recv_matrix, force=.TRUE.)
  END SUBROUTINE dbcsr_sum_replicated

! *****************************************************************************
!> \brief Calls the functions that are called in the included files but not
!>        in thin one. Needed to make prettify think the functions ary really
!>        used.
!> \param[in] matrix          matrix
!> \param[out] temp1, temp2   pointers to the data
! *****************************************************************************
  SUBROUTINE prettify_fooler (matrix, temp1, temp2)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    REAL(KIND=real_4), DIMENSION(:), POINTER :: temp1
    REAL(KIND=real_8), DIMENSION(:), POINTER :: temp2

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
       CALL addto_array (matrix%data_area, temp1)
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

!  INCLUDE "dbcsr_operations_d_.F"
!  INCLUDE "dbcsr_operations_z_.F"
#include "dbcsr_operations_d_.F"
#include "dbcsr_operations_z_.F"
!#if !defined (__SGL)
!  INCLUDE "dbcsr_operations_s_.F"
!  INCLUDE "dbcsr_operations_c_.F"
#include "dbcsr_operations_s_.F"
#include "dbcsr_operations_c_.F"
!#endif

END MODULE dbcsr_operations

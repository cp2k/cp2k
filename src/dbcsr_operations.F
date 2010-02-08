!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
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
  USE dbcsr_kinds,                     ONLY: dp,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8,&
                                             sp
  USE dbcsr_machine,                   ONLY: m_walltime
  USE dbcsr_message_passing,           ONLY: dmp_max,&
                                             mp_allgather,&
                                             mp_alltoall,&
                                             mp_bcast,&
                                             mp_irecv,&
                                             mp_isend,&
                                             mp_sum,&
                                             mp_waitall
!dmp_max, this fixes a bug in g95
! -> ambigous interface for mp_max
!    in module dbcsr_message_passing and wrongly
!    exported from module timings
  USE dbcsr_methods
  USE dbcsr_types
  USE dbcsr_util

  !$ USE OMP_LIB
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_operations'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034


  PUBLIC :: dbcsr_mult_m_e_e, dbcsr_trace, dbcsr_add_on_diag,&
            dbcsr_set, dbcsr_scale, dbcsr_add, dbcsr_copy, dbcsr_get_diag, &
            dbcsr_set_diag, dbcsr_get_block, dbcsr_get_block_p,&
            dbcsr_get_block_diag, dbcsr_hadamard_product, &
            dbcsr_put_block, dbcsr_filter, dbcsr_scale_by_vector,&
            dbcsr_replace_blocks, dbcsr_conjg,&
            dbcsr_btriu, dbcsr_triu, dbcsr_tril, block_add_on_diag,&
            block_chol_inv, block_set, dbcsr_remove_block, dbcsr_norm, &
            dbcsr_gershgorin_norm, dbcsr_maxabs_norm, dbcsr_frobenius_norm, &
            dbcsr_symmetrize_block_diag, dbcsr_copy_columns
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
          dbcsr_trace_ab_s, dbcsr_trace_a_b_d,&!dbcsr_trace_ab_d,&
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
     MODULE PROCEDURE dbcsr_get_block_p_area
  END INTERFACE

  INTERFACE dbcsr_put_block
     MODULE PROCEDURE dbcsr_put_block_area
     MODULE PROCEDURE dbcsr_put_block_d, dbcsr_put_block_s,&
                      dbcsr_put_block_z, dbcsr_put_block_c
     MODULE PROCEDURE dbcsr_put_block2d_d, dbcsr_put_block2d_s,&
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

  INTERFACE dbcsr_scale
     MODULE PROCEDURE dbcsr_scale_anytype
     MODULE PROCEDURE dbcsr_scale_s, dbcsr_scale_d,&
                      dbcsr_scale_c, dbcsr_scale_z
     MODULE PROCEDURE dbcsr_scale_s_m, dbcsr_scale_d_m,&
                      dbcsr_scale_c_m, dbcsr_scale_z_m
  END INTERFACE

  INTERFACE dbcsr_scale_by_vector
     MODULE PROCEDURE dbcsr_scale_by_vector_anytype
     MODULE PROCEDURE dbcsr_scale_by_vector_s, dbcsr_scale_by_vector_d,&
                      dbcsr_scale_by_vector_c, dbcsr_scale_by_vector_z
  END INTERFACE

  INTERFACE dbcsr_set
     MODULE PROCEDURE dbcsr_set_anytype
     MODULE PROCEDURE dbcsr_set_s, dbcsr_set_d, dbcsr_set_c, dbcsr_set_z
  END INTERFACE

  INTERFACE dbcsr_add
     MODULE PROCEDURE dbcsr_add_anytype
     MODULE PROCEDURE dbcsr_add_s, dbcsr_add_d,&
                      dbcsr_add_c, dbcsr_add_z
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

  TYPE dgemm_parameters
     CHARACTER     :: type_a, type_b, type_c
     INTEGER       :: lrows_a, lcols_a, lrows_b, lcols_b, lrows_c, lcols_c
     INTEGER       :: last_k, last_n
     TYPE(dbcsr_scalar_type) :: alpha, beta
     INTEGER       :: offset_a, nze_a, offset_b, nze_b, offset_c, nze_c
  END TYPE dgemm_parameters

CONTAINS

! *****************************************************************************
!> \brief Multiplies two DBCSR matrices
!>
!> \param[in] left_set             set of imaged left matrices
!> \param[in] right_set            set of imaged right matrices
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
  SUBROUTINE dbcsr_mult_m_e_e (left_set, right_set, product_matrix,&
       submatrix_selection, error, retain_sparsity, alpha, beta, flop)
    TYPE(dbcsr_2d_array_type), &
      INTENT(INOUT), TARGET                  :: left_set, right_set
    TYPE(dbcsr_obj), INTENT(INOUT), TARGET   :: product_matrix
    INTEGER, DIMENSION(6), INTENT(IN)        :: submatrix_selection
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: alpha, beta
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mult_m_e_e', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=80), PARAMETER :: &
      fdata = '(A,1X,I4,"(",2(I3),"x",2(I3),")","(",I3,"x",I3,")")', fxfer = &
      '(A,1X,I4,"->",I4,2(1X,"(",I3,"x",I3,")"),1X,"IM (",I3,"x",I3,")"))'
    INTEGER, PARAMETER                       :: idata = 1, imeta = 2, &
                                                M_L = 2, M_P = 1, M_R = 3, &
                                                RC_C = 2, RC_R = 1

    CHARACTER                                :: data_type
    INTEGER :: debug_level, grp, i, ithread, left_col_image, left_col_mult, &
      left_col_nimages, left_col_v_offset, left_dst_icol, left_dst_irow, &
      left_dst_p, left_dst_pcol, left_dst_prow, left_dst_v_col, &
      left_dst_v_row, left_recv_icol, left_recv_irow, left_recv_p, &
      left_recv_pcol, left_recv_prow, left_recv_v_col, left_recv_v_row, &
      left_row_image, left_row_mult, left_row_nimages, left_send_icol, &
      left_send_irow, left_send_p, left_send_pcol, left_send_prow, &
      left_send_v_col, left_send_v_row, left_src_icol, left_src_irow, &
      left_src_p, left_src_pcol, left_src_prow, left_src_v_col, &
      left_src_v_row, max_nblks, max_nze, metronome
    INTEGER :: min_nimages, mp_group, mynode, nsteps_k, nthreads, numnodes, &
      nvirt_k, output_unit, right_col_image, right_col_mult, &
      right_col_nimages, right_dst_icol, right_dst_irow, right_dst_p, &
      right_dst_pcol, right_dst_prow, right_dst_v_col, right_dst_v_row, &
      right_recv_icol, right_recv_irow, right_recv_p, right_recv_pcol, &
      right_recv_prow, right_recv_v_col, right_recv_v_row, right_row_image, &
      right_row_mult, right_row_nimages, right_row_v_offset, right_send_icol, &
      right_send_irow, right_send_p, right_send_pcol, right_send_prow, &
      right_send_v_col, right_send_v_row, right_src_icol, right_src_irow, &
      right_src_p, right_src_pcol
    INTEGER :: right_src_prow, right_src_v_col, right_src_v_row, &
      timing_handle, v_k, v_ki
    INTEGER(KIND=int_8)                      :: flop_single, flop_total
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: infos
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: left_data_rr, left_data_sr, &
      left_index_rr, left_index_sr, right_data_rr, right_data_sr, &
      right_index_rr, right_index_sr
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: left_sizes, my_sizes, &
                                                right_sizes
    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: all_sizes
    INTEGER, DIMENSION(3, 2)                 :: mp_rc_groups
    INTEGER, DIMENSION(:), POINTER           :: left_index_rp, left_index_sp, &
                                                right_index_rp, right_index_sp
    INTEGER, DIMENSION(:, :), POINTER        :: left_pgrid, product_pgrid, &
                                                right_pgrid
    REAL(KIND=dp)                            :: checksum, t_all, t_dgemm, &
                                                trun, trun_t, tstart, tstop
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_2d_array_type), POINTER       :: left_buffer_calc, &
                                                left_buffer_comm, &
                                                right_buffer_calc, &
                                                right_buffer_comm
    TYPE(dbcsr_2d_array_type), TARGET        :: left_buffer_1, left_buffer_2, &
                                                right_buffer_1, right_buffer_2
    TYPE(dbcsr_data_obj)                     :: left_data_rp, left_data_sp, &
                                                right_data_rp, right_data_sp
    TYPE(dbcsr_mp_obj)                       :: left_mp_obj, product_mp_obj, &
                                                right_mp_obj

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    debug_level=0
    t_all = 0.0_dp
    t_dgemm = 0.0_dp
    IF (dbg) debug_level = 2
    CALL timeset (routineN, timing_handle)
    ithread= 0 ; nthreads = 1
    !$ ithread = OMP_GET_THREAD_NUM () ; nthreads = OMP_GET_NUM_THREADS ()
    !IF (.NOT. ASSOCIATED (product_matrix%m%wms)) THEN
       CALL dbcsr_work_create(product_matrix, error=error,&
            work_mutable=.FALSE.)
    !ENDIF
    logger => cp_error_get_logger (error)
    output_unit = cp_logger_get_default_io_unit (logger)
    tstart = 0.0_dp ; tstop = 0.0_dp ; trun = 0.0_dp
    t_dgemm = 0.0_dp ; t_all = 0.0_dp
    flop_total = 0
    trun_t = m_walltime ()
    ! Set up variables
    data_type = dbcsr_get_data_type (product_matrix)
    left_row_nimages =  left_set%image_dist%row_decimation
    left_row_mult =     left_set%image_dist%row_multiplicity
    left_col_nimages =  left_set%image_dist%col_decimation
    left_col_mult =     left_set%image_dist%col_multiplicity
    right_row_nimages = right_set%image_dist%row_decimation
    right_row_mult =    right_set%image_dist%row_multiplicity
    right_col_nimages = right_set%image_dist%col_decimation
    right_col_mult =    right_set%image_dist%col_multiplicity
    left_mp_obj = dbcsr_distribution_mp (left_set%image_dist%main)
    right_mp_obj = dbcsr_distribution_mp (right_set%image_dist%main)
    product_mp_obj = dbcsr_distribution_mp (product_matrix%m%dist)
    numnodes = dbcsr_mp_numnodes (product_mp_obj)
    mynode = dbcsr_mp_mynode (product_mp_obj)
    mp_group = dbcsr_mp_group (product_mp_obj)
    left_pgrid => dbcsr_mp_pgrid (left_mp_obj)
    right_pgrid => dbcsr_mp_pgrid (right_mp_obj)
    product_pgrid => dbcsr_mp_pgrid (product_mp_obj)
    CALL dbcsr_mp_grid_setup (product_mp_obj)
    CALL dbcsr_assert (dbcsr_mp_has_subgroups (product_mp_obj),&
         cp_fatal_level, cp_internal_error, routineN,&
         "Row/column subgroups must be defined")
    mp_rc_groups(M_P, 1:2) = (/ dbcsr_mp_my_row_group (product_mp_obj),&
                                dbcsr_mp_my_col_group (product_mp_obj) /)
    CALL dbcsr_mp_grid_setup (left_mp_obj)
    CALL dbcsr_assert (dbcsr_mp_has_subgroups (left_mp_obj),&
         cp_fatal_level, cp_internal_error, routineN,&
         "Row/column subgroups must be defined")
    mp_rc_groups(M_L, 1:2) = (/ dbcsr_mp_my_row_group (left_mp_obj),&
                                dbcsr_mp_my_col_group (left_mp_obj) /)
    CALL dbcsr_mp_grid_setup (right_mp_obj)
    CALL dbcsr_assert (dbcsr_mp_has_subgroups (right_mp_obj),&
         cp_fatal_level, cp_internal_error, routineN,&
         "Row/column subgroups must be defined")
    mp_rc_groups(M_R, 1:2) = (/ dbcsr_mp_my_row_group (right_mp_obj),&
                                dbcsr_mp_my_col_group (right_mp_obj) /)
    ! Dummy checks
    ! left/right matching
    CALL cp_assert (left_col_nimages .EQ. right_row_mult&
         .AND. left_col_mult .EQ. right_row_nimages, cp_fatal_level,&
         cp_internal_error, routineN, "Left/Right image mismatch", error=error)
    CALL cp_assert (left_col_nimages * dbcsr_mp_npcols (left_mp_obj) &
         .EQ. right_row_nimages * dbcsr_mp_nprows (right_mp_obj), &
         cp_fatal_level, cp_internal_error, routineN,&
         "Left/Right total mismatch", error=error)
    ! product/left matching
    CALL cp_assert (left_row_mult * dbcsr_mp_nprows (product_mp_obj) &
         .EQ. left_row_nimages * dbcsr_mp_nprows (left_mp_obj), &
         cp_fatal_level, cp_internal_error, routineN,&
         "Product/Left total mismatch", error=error)
    ! product/left matching
    CALL cp_assert (right_col_mult * dbcsr_mp_npcols (product_mp_obj) &
         .EQ. right_col_nimages * dbcsr_mp_npcols (right_mp_obj), &
         cp_fatal_level, cp_internal_error, routineN,&
         "Product/Right total mismatch", error=error)
    ! Limitations
    CALL cp_assert (left_row_nimages.eq.1 .AND. left_row_mult.eq.1,&
         cp_fatal_level, cp_unimplemented_error_nr, routineN,&
         "Product/Left matrix process grid mismatch", error=error)
    CALL cp_assert (right_col_nimages.eq.1 .AND. right_col_mult.eq.1,&
         cp_fatal_level, cp_unimplemented_error_nr, routineN,&
         "Product/Right matrix process grid mismatch", error=error)
    !
    ! Setup product work areas
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
    ENDDO
    !
    ! Exchange size data
    ALLOCATE (my_sizes(4, MAX (left_row_nimages, right_row_nimages),&
         MAX (left_col_nimages, right_col_nimages)))
    my_sizes(:,:,:) = 0
    DO left_row_image = 1, left_row_nimages
       DO left_col_image = 1, left_col_nimages
          my_sizes(idata, left_row_image, left_col_image) &
               = dbcsr_get_data_size (&
               left_set%mats(left_row_image, left_col_image)%m%data_area)
          my_sizes(imeta, left_row_image, left_col_image) = SIZE (&
               left_set%mats(left_row_image, left_col_image)%m%index)
       ENDDO
    ENDDO
    DO right_row_image = 1, right_row_nimages
       DO right_col_image = 1, right_col_nimages
          my_sizes(idata+2, right_row_image, right_col_image) &
               = dbcsr_get_data_size (&
               right_set%mats(right_row_image, right_col_image)%m%data_area)
          my_sizes(imeta+2, right_row_image, right_col_image) = SIZE (&
               right_set%mats(right_row_image, right_col_image)%m%index)
       ENDDO
    ENDDO
    ALLOCATE (all_sizes(4, LBOUND(my_sizes,2):UBOUND(my_sizes,2),&
         LBOUND(my_sizes,3):UBOUND(my_sizes,3), 0:numnodes-1))
    CALL mp_allgather(my_sizes, all_sizes, mp_group)
    !
    !
    ! The main transfer loop goes through the virtual
    ! rows/columns. They are grouped together according to the process
    ! coordinate with the smallest number of images.
    min_nimages = MIN (left_col_nimages, right_row_nimages)
    nvirt_k = dbcsr_mp_npcols(left_mp_obj) * left_col_nimages
    nsteps_k = nvirt_k / min_nimages
    left_col_v_offset = MOD (dbcsr_mp_mypcol(left_mp_obj)*left_col_nimages,&
         left_col_mult)
    right_row_v_offset = MOD (dbcsr_mp_myprow(right_mp_obj)*right_row_nimages,&
         right_row_mult)
    !
    ! Translate the all_sizes to account for pre-distribution
    ALLOCATE (left_sizes(2, 0:dbcsr_mp_nprows(left_mp_obj)*left_row_nimages-1, 0:nvirt_k-1))
    left_sizes = -1
    DO left_src_pcol = 0, dbcsr_mp_npcols(left_mp_obj)-1
       DO left_col_image = 1, left_col_nimages
          DO left_src_prow = 0, dbcsr_mp_nprows(left_mp_obj)-1
             left_src_p = left_pgrid (left_src_prow, left_src_pcol)
             DO left_row_image = 1, left_row_nimages
                left_src_v_row = left_src_prow*left_row_nimages + left_row_image-1
                left_src_v_col = left_src_pcol*left_col_nimages + left_col_image-1
                ! Rewind the column
                left_src_v_col = left_src_v_col + left_src_prow*left_col_mult
                left_src_v_col = absmod(left_src_v_col + nvirt_k, nvirt_k)
                left_sizes(idata, left_src_v_row, left_src_v_col) = all_sizes(&
                     idata, left_row_image, left_col_image, left_src_p)
                left_sizes(imeta, left_src_v_row, left_src_v_col) = all_sizes(&
                     imeta, left_row_image, left_col_image, left_src_p)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !write(*,*)"left sizes"
    !write(*,'(2(1X,I7))')left_sizes
    ALLOCATE (right_sizes(2, 0:nvirt_k-1, 0:dbcsr_mp_npcols(right_mp_obj)*right_col_nimages-1))
    right_sizes = -1
    DO right_src_pcol = 0, dbcsr_mp_npcols(right_mp_obj)-1
       DO right_col_image = 1, right_col_nimages
          DO right_src_prow = 0, dbcsr_mp_nprows(right_mp_obj)-1
             right_src_p = right_pgrid (right_src_prow, right_src_pcol)
             DO right_row_image = 1, right_row_nimages
                right_src_v_row = right_src_prow*right_row_nimages + right_row_image-1
                right_src_v_col = right_src_pcol*right_col_nimages + right_col_image-1
                ! Rewind the row
                right_src_v_row = right_src_v_row + right_src_pcol*right_row_mult
                right_src_v_row = absmod(right_src_v_row + nvirt_k, nvirt_k)
                right_sizes(idata, right_src_v_row, right_src_v_col) =&
                     all_sizes(&
                     idata+2, right_row_image, right_col_image, right_src_p)
                right_sizes(imeta, right_src_v_row, right_src_v_col) =&
                     all_sizes(&
                     imeta+2, right_row_image, right_col_image, right_src_p)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !write(*,*)"right sizes"
    !write(*,'(2(1X,I7))')right_sizes
    !
    ! Setup the left buffer matrices
    max_nblks = MAXVAL (all_sizes (imeta, :, :, :))
    max_nze = MAXVAL (all_sizes(idata, :, :, :))
    CALL setup_buffer_matrices (left_buffer_1, left_row_mult, left_col_nimages,&
         left_set%mats(1,1), index_size=max_nblks,&
         data_size=max_nze, error=error)
    CALL setup_buffer_matrices (left_buffer_2, left_row_mult, left_col_nimages,&
         left_set%mats(1,1), index_size=max_nblks,&
         data_size=max_nze, error=error)
    !left_buffer_calc => left_buffer_1
    left_buffer_calc => left_set
    left_buffer_comm => left_buffer_2
    ALLOCATE (left_data_sr  (left_row_nimages, left_col_nimages))
    ALLOCATE (left_index_sr (left_row_nimages, left_col_nimages))
    ALLOCATE (left_data_rr  (left_row_mult, left_col_nimages))
    ALLOCATE (left_index_rr (left_row_mult, left_col_nimages))
    ! Setup buffers for right matrix
    max_nblks = MAXVAL (all_sizes (imeta+2, :,:,:))
    max_nze = MAXVAL (all_sizes(idata+2, :,:,:))
    CALL setup_buffer_matrices (right_buffer_1, right_row_nimages, right_col_mult,&
         right_set%mats(1,1), index_size=max_nblks, data_size=max_nze,&
         error=error)
    CALL setup_buffer_matrices (right_buffer_2, right_row_nimages, right_col_mult,&
         right_set%mats(1,1), index_size=max_nblks, data_size=max_nze,&
         error=error)
    !right_buffer_calc => right_buffer_1
    right_buffer_calc => right_set
    right_buffer_comm => right_buffer_2
    ALLOCATE (right_data_sr  (right_row_nimages, right_col_nimages))
    ALLOCATE (right_index_sr (right_row_nimages, right_col_nimages))
    ALLOCATE (right_data_rr  (right_row_nimages, right_col_mult))
    ALLOCATE (right_index_rr (right_row_nimages, right_col_mult))
    !
    ! Setup the send/receive data pointers
    CALL dbcsr_data_init(left_data_sp)
    CALL dbcsr_data_init(left_data_rp)
    CALL dbcsr_data_init(right_data_sp)
    CALL dbcsr_data_init(right_data_rp)
    CALL dbcsr_data_new(left_data_sp, data_type)
    CALL dbcsr_data_new(left_data_rp, data_type)
    CALL dbcsr_data_new(right_data_sp, data_type)
    CALL dbcsr_data_new(right_data_rp, data_type)
    !
    ALLOCATE (infos (numnodes-1))
    !
    ! In the first loop iteration, the data is fetched from the
    ! sources. In the remaining iterations, the data are exchanged
    ! among neighbors.
    grouped_k_index: DO metronome = 1, nsteps_k
       IF (dbg) WRITE(*,'(1X,A,3(1X,A,1X,I5))')routineN,&
            "step",metronome,&
            "first k",metronome*min_nimages,&
            "last k",(metronome+1)*min_nimages-1
       ! Wait for right matrix transfer completion. Wait in all but
       ! the first loop iteration.
       wait_right: IF (metronome .GT. 1) THEN
          IF (dbg) WRITE (*,'(1X,A)')routineN//" waiting for right"
          CALL mp_waitall (right_data_sr)
          CALL mp_waitall (right_data_rr)
          CALL mp_waitall (right_index_sr)
          CALL mp_waitall (right_index_rr)
       ENDIF wait_right
       ! Right matrix transfer. Transfer in all but the last loop
       ! iteration.
       xfer_right: IF (metronome .LT. nsteps_k) THEN
          DO v_ki = 0, right_row_nimages-1
             v_k = metronome*min_nimages + v_ki
             ! Calculate the process to send to.
             right_send_v_row = v_ki &
                  + dbcsr_mp_myprow(right_mp_obj)*right_row_nimages &
                  - min_nimages
             right_send_v_row = MOD (right_send_v_row + nvirt_k, nvirt_k)
             right_send_prow = right_send_v_row / right_row_nimages
             right_send_irow = MOD (right_send_v_row, right_row_nimages)
             CALL dbcsr_assert (&
                  right_send_prow, "LT", dbcsr_mp_nprows(right_mp_obj),&
                  cp_fatal_level, cp_internal_error, routineN,&
                  "Wrong right send prow")

             right_send_pcol = dbcsr_mp_mypcol(right_mp_obj)
             right_send_icol = 1
             right_send_v_col = dbcsr_mp_mypcol(right_mp_obj)
             ! Calculate which data I send.
             right_dst_v_row = v_k &
                  + dbcsr_mp_myprow(right_mp_obj)*right_row_nimages &
                  + dbcsr_mp_mypcol(right_mp_obj)*right_row_mult &
                  - min_nimages
             right_dst_v_row = MOD (right_dst_v_row + nvirt_k, nvirt_k)
             right_dst_prow = right_dst_v_row / right_row_nimages
             right_dst_irow = 1 + MOD (right_dst_v_row, right_row_nimages)
             CALL dbcsr_assert(&
                  right_dst_prow, "LT", dbcsr_mp_nprows(right_mp_obj),&
                  cp_fatal_level, cp_internal_error, routineN,&
                  "Wrong right dst prow")
             CALL dbcsr_assert(&
                  right_dst_irow, "LE", right_row_nimages,&
                  cp_fatal_level, cp_internal_error, routineN,&
                  "Wrong right dst row image")
             right_dst_pcol = dbcsr_mp_mypcol(right_mp_obj)
             right_dst_icol = 1
             right_dst_v_col =&
                  dbcsr_mp_mypcol(right_mp_obj)*right_col_nimages
             !
             right_dst_p = right_pgrid(right_dst_prow, right_dst_pcol)
             CALL dbcsr_data_set_pointer(&
                  area=right_data_sp,&
                  rsize=right_sizes(idata, right_dst_v_row, right_dst_v_col),&
                  csize=1,&
                  pointee=right_buffer_calc%mats(v_ki+1, 1)%m%data_area)
             !right_data_sp => pointer_view (&
             !     dbcsr_get_data_p (right_buffer_calc%mats(&
             !     v_ki+1, 1&
             !     )%m%data_area, 0.0_dp),&
             !     1, right_sizes(idata, right_dst_v_row, right_dst_v_col))
             right_index_sp => right_buffer_calc%mats(&
                  v_ki+1, 1&
                  )%m%index(1:&
                  right_sizes(imeta, right_dst_v_row, right_dst_v_col))
             !
             ! Calculate the process to receive from
             right_recv_v_row = v_ki &
                  + dbcsr_mp_myprow(right_mp_obj)*right_row_nimages &
                  + min_nimages
             right_recv_v_row = MOD(right_recv_v_row + nvirt_k, nvirt_k)
             right_recv_prow = right_recv_v_row / right_row_nimages
             right_recv_irow = MOD(right_recv_v_row, right_row_nimages)
             CALL dbcsr_assert (&
                  right_recv_prow, "LT", dbcsr_mp_nprows(right_mp_obj),&
                  cp_fatal_level, cp_internal_error, routineN,&
                  "Wrong right receive prow")
             right_recv_pcol = dbcsr_mp_mypcol (right_mp_obj)
             right_recv_icol = 1
             right_recv_v_col = dbcsr_mp_mypcol(right_mp_obj)*right_col_nimages
             ! Calculate which data I receive
             right_src_v_row = v_k &
                  + dbcsr_mp_myprow(right_mp_obj)*right_row_nimages &
                  + dbcsr_mp_mypcol(right_mp_obj)*right_row_mult
             right_src_v_row = MOD (right_src_v_row, nvirt_k)
             right_src_prow = right_src_v_row / right_row_nimages
             right_src_irow = 1 + MOD (right_src_v_row, right_row_nimages)
             CALL dbcsr_assert(&
                  right_src_prow, "LT", dbcsr_mp_nprows(right_mp_obj),&
                  cp_fatal_level, cp_internal_error, routineN,&
                  "Wrong right src prow")
             CALL dbcsr_assert(&
                  right_src_irow, "LE", right_row_nimages,&
                  cp_fatal_level, cp_internal_error, routineN,&
                  "Wrong right src row image")
             right_src_pcol = dbcsr_mp_mypcol(right_mp_obj)
             right_src_icol = 1
             right_src_v_col = dbcsr_mp_mypcol(right_mp_obj)*right_col_nimages
             !
             right_src_p = right_pgrid(right_src_prow, right_src_pcol)
             CALL dbcsr_data_set_pointer(&
                  area=right_data_rp,&
                  rsize=right_sizes(idata, right_src_v_row, right_src_v_col),&
                  csize=1,&
                  pointee=right_buffer_comm%mats(v_ki+1, 1)%m%data_area)
             !right_data_rp => pointer_view (&
             !     dbcsr_get_data_p (right_buffer_comm%mats(&
             !        v_ki+1, 1&
             !     )%m%data_area, 0.0_dp),&
             !     1, right_sizes(idata, right_src_v_row, right_src_v_col))
             right_index_rp => right_buffer_comm%mats(&
                     v_ki+1, 1&
                  )%m%index(1:&
                     right_sizes(imeta, right_src_v_row, right_src_v_col))
             !
             right_send_p = right_pgrid (right_send_prow, right_send_pcol)
             right_recv_p = right_pgrid (right_recv_prow, right_recv_pcol)
             IF (dbg) THEN
                CALL all_print (right_send_p, "R SEND P", product_mp_obj)
                CALL all_print (right_recv_p, "R RECV P", product_mp_obj)
                CALL all_print (right_src_p, "R SRC P ", product_mp_obj)
                CALL all_print (right_dst_p, "R DST P ", product_mp_obj)
                WRITE(*,fxfer)"SEND RIGHT",mynode,right_send_p,&
                     dbcsr_mp_myprow(right_mp_obj),dbcsr_mp_mypcol(right_mp_obj),&
                     right_send_prow, right_send_pcol,&
                     right_send_v_row, right_send_v_col
                WRITE(*,fdata)"WHAT RIGHT",right_dst_p,&
                     right_dst_prow, right_dst_irow,&
                     right_dst_pcol, right_dst_icol
                WRITE(*,fxfer)"RECV RIGHT",right_recv_p,mynode,&
                     right_recv_prow, right_recv_pcol,&
                     dbcsr_mp_myprow(right_mp_obj),dbcsr_mp_mypcol(right_mp_obj),&
                     right_recv_v_row, right_recv_v_col
                WRITE(*,fdata)"WHAT RIGHT",right_src_p,&
                     right_src_prow, right_src_irow,&
                     right_src_pcol, right_src_icol
             ENDIF
             IF (dbg) &
                  WRITE(*,fdata)"CALC RIGHT",right_dst_p,&
                  right_dst_prow, right_dst_irow,&
                  right_dst_pcol, right_dst_icol,&
                  right_dst_v_row, right_dst_v_col
             ! These are column-communicator relative
             right_send_p = right_send_prow
             right_recv_p = right_recv_prow
             grp = mp_rc_groups(M_R, RC_C)
             SELECT CASE (data_type)
             CASE (dbcsr_type_real_4)
                CALL mp_isend (right_data_sp%d%r_sp, right_send_p,&
                     grp, right_data_sr(v_ki+1, 1), tag=right_dst_v_row)
                CALL mp_irecv (right_data_rp%d%r_sp, right_recv_p,&
                     grp, right_data_rr(v_ki+1, 1), tag=right_src_v_row)
             CASE (dbcsr_type_real_8)
                CALL mp_isend (right_data_sp%d%r_dp, right_send_p,&
                     grp, right_data_sr(v_ki+1, 1), tag=right_dst_v_row)
                CALL mp_irecv (right_data_rp%d%r_dp, right_recv_p,&
                     grp, right_data_rr(v_ki+1, 1), tag=right_src_v_row)
             CASE (dbcsr_type_complex_4)
                CALL mp_isend (right_data_sp%d%c_sp, right_send_p,&
                     grp, right_data_sr(v_ki+1, 1), tag=right_dst_v_row)
                CALL mp_irecv (right_data_rp%d%c_sp, right_recv_p,&
                     grp, right_data_rr(v_ki+1, 1), tag=right_src_v_row)
             CASE (dbcsr_type_complex_8)
                CALL mp_isend (right_data_sp%d%c_dp, right_send_p,&
                     grp, right_data_sr(v_ki+1, 1), tag=right_dst_v_row)
                CALL mp_irecv (right_data_rp%d%c_dp, right_recv_p,&
                     grp, right_data_rr(v_ki+1, 1), tag=right_src_v_row)
             END SELECT
             CALL mp_isend (right_index_sp, right_send_p,&
                  grp, right_index_sr(v_ki+1, 1), tag=right_dst_v_row)
             CALL mp_irecv (right_index_rp, right_recv_p,&
                  grp, right_index_rr(v_ki+1, 1), tag=right_src_v_row)
          ENDDO
       ENDIF xfer_right
       !
       ! Repoint indices of right matrices
       calc_case_right: IF (metronome .GT. 1) THEN
          DO v_ki = 0, right_row_nimages-1
             CALL dbcsr_repoint_index (right_buffer_calc%mats(v_ki+1,1)%m)
             right_buffer_calc%mats(v_ki+1,1)%m%valid = .TRUE.
          ENDDO
       ENDIF calc_case_right
       !
       ! Wait for left matrix transfer completion. Wait in all but
       ! the first loop iteration.
       wait_left: IF (metronome .GT. 1) THEN
          IF (dbg) WRITE (*,'(1X,A)')routineN//" waiting for left"
          CALL mp_waitall (left_data_sr)
          CALL mp_waitall (left_data_rr)
          CALL mp_waitall (left_index_sr)
          CALL mp_waitall (left_index_rr)
       ENDIF wait_left
       ! Left matrix transfer. Transfer in all but the last loop
       ! iteration.
       xfer_left: IF (metronome .LT. nsteps_k) THEN
          DO v_ki = 0, left_col_nimages-1
             v_k = metronome*min_nimages + v_ki
             ! Calculate the process to send to.
             left_send_prow = dbcsr_mp_myprow(left_mp_obj)
             left_send_irow = 1
             left_send_v_row = dbcsr_mp_myprow(left_mp_obj)
             !
             left_send_v_col = dbcsr_mp_mypcol(left_mp_obj)*left_col_nimages &
                  - min_nimages + v_ki
             left_send_v_col = MOD (left_send_v_col + nvirt_k, nvirt_k)
             left_send_pcol = left_send_v_col / left_col_nimages
             left_send_icol = 1
             CALL dbcsr_assert (&
                  left_send_pcol, "LT", dbcsr_mp_npcols(left_mp_obj),&
                  cp_fatal_level, cp_internal_error, routineN,&
                  "Wrong left send pcol")
             ! Calculate which data I send.
             left_dst_prow = dbcsr_mp_myprow(left_mp_obj)
             left_dst_irow = 1
             left_dst_v_row = dbcsr_mp_myprow(left_mp_obj)*left_row_nimages
             !
             left_dst_v_col = v_k &
                  + dbcsr_mp_mypcol(left_mp_obj)*left_col_nimages &
                  + dbcsr_mp_myprow(left_mp_obj)*left_col_mult &
                  - min_nimages
             left_dst_v_col = MOD (left_dst_v_col + nvirt_k, nvirt_k)
             left_dst_pcol = left_dst_v_col / left_col_nimages
             left_dst_icol = 1 + MOD (left_dst_v_col, left_col_nimages)
             CALL dbcsr_assert(&
                  left_dst_pcol, "LT", dbcsr_mp_npcols(left_mp_obj),&
                  cp_fatal_level, cp_internal_error, routineN,&
                  "Wrong left dst pcol")
             CALL dbcsr_assert(&
                  left_dst_icol, "LE", left_col_nimages,&
                  cp_fatal_level, cp_internal_error, routineN,&
                  "Wrong left dst col image")
             !
             left_dst_p = left_pgrid (left_dst_prow, left_dst_pcol)
             CALL dbcsr_data_set_pointer(&
                  area=left_data_sp,&
                  rsize=left_sizes(idata, left_dst_v_row, left_dst_v_col),&
                  csize=1,&
                  pointee=left_buffer_calc%mats(1, v_ki+1)%m%data_area)
             !left_data_sp => pointer_view (&
             !     dbcsr_get_data_p (left_buffer_calc%mats(&
             !        1, v_ki+1&
             !     )%m%data_area, 0.0_dp),&
             !     1, left_sizes(idata, left_dst_v_row, left_dst_v_col))
             left_index_sp => left_buffer_calc%mats(&
                     1, v_ki+1&
                  )%m%index(1:&
                     left_sizes(imeta, left_dst_v_row, left_dst_v_col))
             !
             ! Calculate the process to receive from
             left_recv_prow = dbcsr_mp_myprow (left_mp_obj)
             left_recv_irow = 1
             left_recv_v_row = dbcsr_mp_myprow (left_mp_obj)
             !
             left_recv_v_col = dbcsr_mp_mypcol(left_mp_obj)*left_col_nimages &
                  + min_nimages + v_ki
             left_recv_v_col = MOD(left_recv_v_col + nvirt_k, nvirt_k)
             left_recv_pcol = left_recv_v_col / left_col_nimages
             left_recv_icol = 1
             CALL dbcsr_assert (&
                  left_recv_pcol, "LT", dbcsr_mp_npcols(left_mp_obj),&
                  cp_fatal_level, cp_internal_error, routineN,&
                  "Wrong receive pcol")
             ! Calculate which data I receive
             left_src_prow = dbcsr_mp_myprow (left_mp_obj)
             left_src_irow = 1
             left_src_v_row = dbcsr_mp_myprow(left_mp_obj)*left_row_nimages
             !
             left_src_v_col = v_k &
                  !+ dbcsr_mp_myprow(right_mp_obj)*right_row_nimages &
                  !+ dbcsr_mp_mypcol(right_mp_obj)*right_row_mult
                  + dbcsr_mp_mypcol(left_mp_obj)*left_col_nimages &
                  + dbcsr_mp_myprow(left_mp_obj)*left_col_mult
             left_src_v_col = MOD(left_src_v_col, nvirt_k)
             left_src_pcol = left_src_v_col / left_col_nimages
             left_src_icol = 1 + MOD (left_src_v_col, left_col_nimages)
             CALL dbcsr_assert(&
                  left_src_pcol, "LT", dbcsr_mp_npcols(left_mp_obj),&
                  cp_fatal_level, cp_internal_error, routineN,&
                  "Wrong left src pcol")
             CALL dbcsr_assert(&
                  left_src_icol, "LE", left_col_nimages,&
                  cp_fatal_level, cp_internal_error, routineN,&
                  "Wrong left src column image")
             !
             left_src_p = left_pgrid (left_src_prow, left_src_pcol)
             CALL dbcsr_data_set_pointer(&
                  area=left_data_rp,&
                  rsize=left_sizes(idata, left_src_v_row, left_src_v_col),&
                  csize=1,&
                  pointee=left_buffer_comm%mats(1, v_ki+1)%m%data_area)
             !left_data_rp => pointer_view (&
             !     dbcsr_get_data_p (left_buffer_comm%mats(&
             !        1, v_ki+1&
             !     )%m%data_area, 0.0_dp),&
             !     1, left_sizes(idata, left_src_v_row, left_src_v_col))
             left_index_rp => left_buffer_comm%mats(&
                     1, v_ki+1&
                  )%m%index(1:&
                     left_sizes(imeta, left_src_v_row, left_src_v_col))
             !
             left_send_p = left_pgrid (left_send_prow, left_send_pcol)
             left_recv_p = left_pgrid (left_recv_prow, left_recv_pcol)
             IF (dbg) THEN
                CALL all_print (left_send_p, "L SEND P", product_mp_obj)
                CALL all_print (left_recv_p, "L RECV P", product_mp_obj)
                CALL all_print (left_src_p, "L SRC P ", product_mp_obj)
                CALL all_print (left_dst_p, "L DST P ", product_mp_obj)
                WRITE(*,fxfer)"SEND LEFT ",mynode,left_send_p,&
                     left_send_prow, left_send_pcol,&
                     dbcsr_mp_myprow(left_mp_obj),dbcsr_mp_mypcol(left_mp_obj),&
                     left_send_v_row, left_send_v_col
                WRITE(*,fdata)"WHAT LEFT ",left_dst_p,&
                     left_dst_prow, left_dst_irow,&
                     left_dst_pcol, left_dst_icol
                WRITE(*,fxfer)"RECV LEFT ",left_recv_p,mynode,&
                     left_recv_prow,left_recv_pcol,&
                     dbcsr_mp_myprow(left_mp_obj),dbcsr_mp_mypcol(left_mp_obj),&
                     left_recv_v_row, left_recv_v_col
                WRITE(*,fdata)"WHAT LEFT ",left_src_p,&
                     left_src_prow, left_src_irow,&
                     left_src_pcol, left_src_icol
             ENDIF
             IF (dbg) &
                  WRITE(*,fdata)"CALC LEFT ",left_dst_p,&
                  left_dst_prow, left_dst_irow,&
                  left_dst_pcol, left_dst_icol,&
                  left_dst_v_row, left_dst_v_col
             ! These are column-communicator relative
             left_send_p = left_send_pcol
             left_recv_p = left_recv_pcol
             grp = mp_rc_groups(M_L, RC_R)
             SELECT CASE (data_type)
             CASE (dbcsr_type_real_4)
                CALL mp_isend (left_data_sp%d%r_sp, left_send_p,&
                     grp, left_data_sr(1, v_ki+1), tag=left_dst_v_col)
                CALL mp_irecv (left_data_rp%d%r_sp, left_recv_p,&
                     grp, left_data_rr(1, v_ki+1), tag=left_src_v_col)
             CASE (dbcsr_type_real_8)
                CALL mp_isend (left_data_sp%d%r_dp, left_send_p,&
                     grp, left_data_sr(1, v_ki+1), tag=left_dst_v_col)
                CALL mp_irecv (left_data_rp%d%r_dp, left_recv_p,&
                     grp, left_data_rr(1, v_ki+1), tag=left_src_v_col)
             CASE (dbcsr_type_complex_4)
                CALL mp_isend (left_data_sp%d%c_sp, left_send_p,&
                     grp, left_data_sr(1, v_ki+1), tag=left_dst_v_col)
                CALL mp_irecv (left_data_rp%d%c_sp, left_recv_p,&
                     grp, left_data_rr(1, v_ki+1), tag=left_src_v_col)
             CASE (dbcsr_type_complex_8)
                CALL mp_isend (left_data_sp%d%c_dp, left_send_p,&
                     grp, left_data_sr(1, v_ki+1), tag=left_dst_v_col)
                CALL mp_irecv (left_data_rp%d%c_dp, left_recv_p,&
                     grp, left_data_rr(1, v_ki+1), tag=left_src_v_col)
             END SELECT
             CALL mp_isend (left_index_sp, left_send_p,&
                  grp, left_index_sr(1, v_ki+1), tag=left_dst_v_col)
             CALL mp_irecv (left_index_rp, left_recv_p,&
                  grp, left_index_rr(1, v_ki+1), tag=left_src_v_col)
          ENDDO
       ENDIF xfer_left
       !
       ! Repoint indices of left matrices and do the multiplications.
       calc_case_left: IF (metronome .GT. 0) THEN
          IF (metronome .GT. 1) THEN
             DO v_ki = 0, left_col_nimages-1
                CALL dbcsr_repoint_index (left_buffer_calc%mats(1,v_ki+1)%m)
                left_buffer_calc%mats(1, v_ki+1)%m%valid=.TRUE.
             ENDDO
          ENDIF
          DO v_ki = 0, min_nimages-1
             IF (dbg) THEN
                CALL dbcsr_print(left_buffer_calc%mats(1, v_ki+1), nodata=.TRUE.)
                CALL dbcsr_print(right_buffer_calc%mats(v_ki+1, 1), nodata=.TRUE.)
             ENDIF
             tstart = m_walltime ()
             CALL dbcsr_nn_mult_lin(&
                  left_buffer_calc%mats(1, v_ki+1)%m,&
                  right_buffer_calc%mats(v_ki+1, 1)%m,&
                  product_matrix%m, submatrix_box=submatrix_selection,&
                  retain_sparsity=retain_sparsity,&
                  alpha=alpha,&
                  flop=flop_single,&
                  t_all=t_all, t_dgemm=t_dgemm)
             flop_total = flop_total + flop_single
             tstop = m_walltime ()
             IF (tstop-tstart .EQ. 0) tstop = tstart+0.000001_dp
             IF ((output_unit>0) .AND. detailed_timing) THEN
                WRITE(output_unit,'(1X,A,F9.4,A,EN12.4,A,1X,EN12.4,1X,A)')&
                     "Segment Local Multiplication time ",tstop-tstart," and ",&
                     (REAL(flop_single, dp)/1000000.0_dp)/(tstop-tstart),&
                     " MFLOPS", REAL(flop_single), "FLOP"
             ENDIF
             trun = trun + (tstop - tstart)
          ENDDO
       ENDIF calc_case_left
       IF (metronome .EQ. 1) THEN
          left_buffer_calc => left_buffer_1
          right_buffer_calc => right_buffer_1
       ENDIF
       CALL dbcsr_switch_sets (left_buffer_calc, left_buffer_comm)
       CALL dbcsr_switch_sets (right_buffer_calc, right_buffer_comm)
    ENDDO grouped_k_index
    trun_t = m_walltime () - trun_t
    IF (trun_t .EQ. 0) trun_t = 0.000001_dp
    IF (trun .EQ. 0) trun = 0.000001_dp
    IF ((output_unit>0) .AND. detailed_timing) THEN
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
    !
    CALL dbcsr_destroy_array (right_buffer_1)
    CALL dbcsr_destroy_array (right_buffer_2)
    CALL dbcsr_destroy_array (left_buffer_1)
    CALL dbcsr_destroy_array (left_buffer_2)
    DEALLOCATE (my_sizes)
    !
    CALL dbcsr_data_clear_pointer(left_data_sp)
    CALL dbcsr_data_clear_pointer(left_data_rp)
    CALL dbcsr_data_clear_pointer(right_data_sp)
    CALL dbcsr_data_clear_pointer(right_data_rp)
    CALL dbcsr_data_release(left_data_sp)
    CALL dbcsr_data_release(left_data_rp)
    CALL dbcsr_data_release(right_data_sp)
    CALL dbcsr_data_release(right_data_rp)
    !
    CALL dbcsr_finalize(product_matrix)
    !
    IF (dbg) THEN
       checksum = dbcsr_checksum (product_matrix, error=error)
       IF ((output_unit>0)) THEN
          WRITE(output_unit,'(1X,A,1X,F9.4)')"Product Checksum=",checksum
       ENDIF
    ENDIF
    !
    IF (PRESENT (flop)) flop = flop_total
    CALL timestop (timing_handle)
  END SUBROUTINE dbcsr_mult_m_e_e


  SUBROUTINE setup_buffer_matrices (buffer_set, buff_rows, buff_cols,&
       source_matrix, index_size, data_size, error)
    TYPE(dbcsr_2d_array_type), INTENT(OUT)   :: buffer_set
    INTEGER, INTENT(IN)                      :: buff_rows, buff_cols
    TYPE(dbcsr_obj), INTENT(IN)              :: source_matrix
    INTEGER, INTENT(IN)                      :: index_size, data_size
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: col_image, row_image

    ALLOCATE (buffer_set%mats(buff_rows, buff_cols))
    DO row_image = 1, buff_rows
       DO col_image = 1, buff_cols
          CALL dbcsr_init(buffer_set%mats(row_image, col_image), error)
          CALL dbcsr_create(buffer_set%mats(row_image, col_image),&
               "Buffer of "//source_matrix%m%name,&
               dbcsr_distribution (source_matrix),&
               dbcsr_get_matrix_type (source_matrix),&
               dbcsr_row_block_sizes (source_matrix),&
               dbcsr_col_block_sizes (source_matrix),&
               nblks=index_size, nze=data_size, &
               data_type=dbcsr_get_data_type(source_matrix),&
               special=.TRUE., error=error)
          CALL dbcsr_data_ensure_size (&
               buffer_set%mats(row_image, col_image)%m%data_area,&
               data_size, nocopy=.TRUE.)
          CALL ensure_array_size (&
               buffer_set%mats(row_image, col_image)%m%index,&
               ub=index_size, nocopy=.TRUE., special=.TRUE.,&
               error=error)
          buffer_set%mats(row_image, col_image)%m%sym_negation&
               = source_matrix%m%sym_negation
       ENDDO
    ENDDO
  END SUBROUTINE setup_buffer_matrices

  SUBROUTINE all_print (val, desc, mp)
    INTEGER, INTENT(IN)                      :: val
    CHARACTER(LEN=*), INTENT(IN)             :: desc
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp

    CHARACTER(LEN=80)                        :: fmt
    INTEGER                                  :: n
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: infos

!

    n = dbcsr_mp_numnodes (mp)
    WRITE(fmt,'("(1X,A,",I5,"(1X,I4))")')n
    ALLOCATE (infos(n))
    CALL mp_allgather(val, infos, dbcsr_mp_group(mp))
    WRITE(*,fmt)desc,infos
  END SUBROUTINE all_print


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
    TYPE(dbcsr_2d_array_type), POINTER       :: set1p, set2p

    TYPE(dbcsr_2d_array_type), POINTER       :: tmp_set

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



  ELEMENTAL FUNCTION make_coordinate_tuple (most, least) RESULT (tuple)
    INTEGER, INTENT(IN)                      :: most, least
    INTEGER(KIND=int_8)                      :: tuple

!tuple = IOR (ISHFT (most, 32), least)

    tuple = most
    tuple = IOR (ISHFT(tuple, 32), INT(least, int_8))
  END FUNCTION make_coordinate_tuple


! *****************************************************************************
!> \brief Multiplies two DBCSR matrices, using the current local data.
!>
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
       retain_sparsity, alpha, beta,&
       t_all, t_dgemm)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    INTEGER, DIMENSION(6), INTENT(in)        :: submatrix_box
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: alpha, beta
    REAL(KIND=dp), INTENT(INOUT)             :: t_all, t_dgemm

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_nn_mult_lin', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: a_type, b_type, data_type
    INTEGER :: a_blk, a_col, a_col_logical, a_col_size, a_col_size_logical, &
      a_nze, a_offset, a_row, a_row_l, a_row_logical, a_row_size, &
      a_row_size_logical, b_blk, b_col, b_col_logical, b_col_size, &
      b_col_size_logical, b_nze, b_offset, b_row, b_row_logical, b_row_size, &
      b_row_size_logical, blk_end, blk_size, blk_start, c_blk_p, c_blk_pt, &
      c_col, c_col_logical, c_col_size, c_col_size_logical, c_nze, c_offset, &
      c_row_logical, c_row_size, c_row_size_logical, clamp_k_len, datasize, &
      first_k, first_row, first_working_sp, i, ithread, j, last_k, last_row, &
      lastblk, new_blk, new_row_blks, nthreads, old_blk, old_row_begin_blk
    INTEGER :: old_size, row_begin_blk, row_blk, rpt, size_c_blk_ps, &
      stack_pointer, stat, submat_f_col, submat_f_row, submat_l_col, &
      submat_l_row, this_right_limit
    INTEGER(KIND=int_4)                      :: offset
    INTEGER(KIND=int_8)                      :: coordinate, lflop, one_flop
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: c_blk_ps, c_blk_ps_p
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: new_index
    INTEGER, DIMENSION(:), POINTER :: left_col_blk_offset, left_local_rows, &
      left_row_blk_offset, new_blk_p, new_col_i, new_row_i, &
      product_col_blk_offset, right_col_blk_offset
    LOGICAL :: block_exists, found, keep_sparsity, left_limit_rows, limit_k, &
      negate_alpha, product_is_symmetric, right_limit_cols, tr_a, tr_b, tr_p
    REAL(KIND=dp)                            :: epoch, mt_t_dgemm
    TYPE(dbcsr_scalar_type)                  :: beta_add, beta_new, my_alpha, &
                                                my_beta, which_beta
    TYPE(dgemm_parameters), DIMENSION(16)    :: dgemm_parameter_stack

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    data_type = product%data_type
    my_alpha = dbcsr_scalar_one (data_type)
    IF (PRESENT (alpha)) my_alpha = alpha
    my_beta = dbcsr_scalar_zero (data_type)
    IF (PRESENT (beta)) my_beta = beta
    beta_new = dbcsr_scalar_zero (data_type)
    beta_add = dbcsr_scalar_one (data_type)
    keep_sparsity = .FALSE.
    IF (PRESENT (retain_sparsity)) THEN
       !WRITE(*,*)routineN//" Retaining sparsity!"
       keep_sparsity = .TRUE.
    END IF
    product_is_symmetric = (product%matrix_type .EQ. dbcsr_type_symmetric &
         .OR. product%matrix_type .EQ. dbcsr_type_antisymmetric)
    ithread = 0 ; nthreads = 1
    !$ ithread = OMP_GET_THREAD_NUM () ; nthreads = OMP_GET_NUM_THREADS ()
    !WRITE(*,*)'Hello! I am thread',ithread
    epoch = m_walltime()
    DBG "Doing NN matrix multiplication",left%name,' x ',right%name
    DBG "Pre-mult blk sizes, data", product%wms(ithread+1)%lastblk,&
         product%wms(ithread+1)%datasize
    left_row_blk_offset => array_data (left%row_blk_offset)
    left_col_blk_offset => array_data (left%col_blk_offset)
    right_col_blk_offset => array_data (right%col_blk_offset)
    product_col_blk_offset => array_data (product%col_blk_offset)
    !left_thread_dist => array_data (dbcsr_distribution_thread_dist (product%dist))
    left_local_rows => array_data (dbcsr_distribution_local_rows (left%dist))
    ! Submatrix handling
    CALL cp_assert (submatrix_box(1) .EQ. 0,&
         cp_fatal_level, cp_unimplemented_error_nr, routineN,&
         "Limiting first row not supported")
    CALL cp_assert (submatrix_box(3) .EQ. 0,&
         cp_fatal_level, cp_unimplemented_error_nr, routineN,&
         "Limiting first column not supported")
    CALL cp_assert (submatrix_box(5) .EQ. 0,&
         cp_fatal_level, cp_unimplemented_error_nr, routineN,&
         "Limiting first k not supported")
    left_limit_rows = MAXVAL (submatrix_box(1:2)) .GT. 0
    right_limit_cols = MAXVAL (submatrix_box(3:4)) .GT. 0
    limit_k = MAXVAL (submatrix_box(5:6)) .GT. 0
    IF (left_limit_rows) THEN
       submat_f_row = submatrix_box(1)
       submat_l_row = submatrix_box(2)
       !write(*,*)"limiting rows",submatrix_box(1:2)
    ENDIF
    IF (right_limit_cols) THEN
       submat_f_col = submatrix_box(3)
       submat_l_col = submatrix_box(4)
       !write(*,*)"limiting cols",submatrix_box(3:4)
    ENDIF
    IF (limit_k) THEN
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
    c_blk_p = product%wms(ithread+1)%datasize + 1
    DBG 'A rows:',left%nblkrows_total
    old_size = product%wms(ithread+1)%lastblk
    mt_t_dgemm = 0.0_dp
    IF (left%nblks .GT. 0 .AND. right%nblks .GT. 0) THEN
       !
    ALLOCATE (new_index (3, left%nfullrows_local*right%nfullcols_local))
    ! c_blk_ps contains the block pointers in the product's current row
    ALLOCATE (c_blk_ps(product%nblkcols_total), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineN, "c_blk_ps",&
         product%nblkcols_total)
    c_blk_ps(:) = 0
    ALLOCATE (c_blk_ps_p(product%nblkcols_total), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineN, "c_blk_ps_p",&
         product%nblkcols_total)
    size_c_blk_ps = 0
    c_blk_pt = product%wms(ithread+1)%datasize + 1
    lastblk = 0
    datasize = 0
    !
    mt_t_dgemm = 0.0_dp
    lflop = 0
    old_blk = 1
    first_working_sp = 1
    stack_pointer = 0
    new_blk = 0
    !write(*,*)ithread,':',left_thread_dist
    rpt = CEILING (REAL(left%nblkrows_local) / REAL(nthreads))
    first_row = MIN (ithread*rpt+1, left%nblkrows_local)
    last_row = MIN ((ithread+1)*rpt, left%nblkrows_local)
    !write(*,*)ithread,'Doing rows', first_row, last_row
    DO a_row_l = 1, left%nblkrows_local
       a_row = left_local_rows (a_row_l)
       a_row_logical = a_row
       IF (dbg) WRITE(*,*)routineN//" A row", a_row
       blk_start = left_row_blk_offset(a_row_logical)
       blk_end = left_row_blk_offset(a_row_logical+1)-1
       blk_size = left_row_blk_offset(a_row_logical+1) - blk_start
       IF (left_limit_rows) THEN
          IF (left_row_blk_offset (a_row) .GT. submat_l_row .OR.&
              left_row_blk_offset (a_row+1)-1 .LT. submat_f_row) THEN
             !write(*,*)"Skipping row", a_row
             CYCLE
          ELSE
             CALL cp_assert(&
                    MIN(submat_l_row, blk_start)&
                  - MAX(submat_f_row, blk_end)+1&
                  .EQ. blk_size,&
                  cp_fatal_level, cp_unimplemented_error_nr,&
                  routineN, "Submatrix must encompass complete block rows")
          ENDIF
       ENDIF
       a_row_size_logical = blk_size ! stored
       c_row_size_logical = blk_size
       row_begin_blk = new_blk+1
       IF (old_blk .LE. old_size) THEN
          DO WHILE (product%wms(ithread+1)%row_i(old_blk) .LT. a_row)
             old_blk = old_blk+1
             IF (old_blk .GT. old_size) THEN
                !old_blk = old_blk-1
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
       IF (dbg) WRITE(*,*)routineN//" old_blk, old_size, old_row_begin_blk",&
            old_blk,old_size,old_row_begin_blk
       IF (old_blk .LE. old_size) THEN
          ! Old index
          DO WHILE (product%wms(ithread+1)%row_i(old_blk) .LE. a_row) ! EQ
             DBG "Adding c_blk_ps column",&
                  product%wms(ithread+1)%col_i(old_blk),&
                  product%wms(ithread+1)%blk_p(old_blk)
             c_col = product%wms(ithread+1)%col_i(old_blk)
             found = .FALSE.
             c_blk_ps(c_col)&
                  = product%wms(ithread+1)%blk_p(old_blk)
             size_c_blk_ps = size_c_blk_ps + 1
             c_blk_ps_p(size_c_blk_ps) = product%wms(ithread+1)%col_i(old_blk)
             old_blk = old_blk+1
             IF (old_blk .GT. old_size) THEN
                !old_blk = old_blk-1
                EXIT
             ENDIF
          ENDDO
       ENDIF
       IF (dbg) THEN
          WRITE(*,*)routineN//" c_blk_ps="
          WRITE(*,'(10(1X,I7))')c_blk_ps
       ENDIF
       new_row_blks = 0
       DO a_blk = left%row_p(a_row)+1, left%row_p(a_row+1)
          IF (left%blk_p(a_blk) .EQ. 0) THEN
             !write(*,*)"Skipping deleted block"
             CYCLE
          ENDIF
          row_blk = row_begin_blk
          old_blk = old_row_begin_blk
          a_offset = ABS (left%blk_p(a_blk))
          a_col = left%col_i(a_blk)
          IF (dbg) WRITE(*,*)routineN//" A col", a_col,";",a_row
          a_row_logical = a_row
          a_col_logical = a_col
          blk_start = left_col_blk_offset (a_col_logical)
          blk_end = left_col_blk_offset (a_col_logical+1)-1
          blk_size = left_col_blk_offset (a_col_logical+1) - blk_start
          a_col_size_logical = blk_size
          a_row_size = a_row_size_logical
          a_col_size = a_col_size_logical
          tr_a = left%blk_p(a_blk) .LT. 0
          a_type = blas_mat_type (tr_a)
          IF (tr_a) CALL swap (a_row_size, a_col_size)
          IF (limit_k) THEN
             IF (blk_start .GT. last_k .OR.&
                 blk_end .LT. first_k) THEN
                !write(*,'(A,I5,"(",I5,"-",I5,") not in ",I5,"-",I5)')&
                !     routineN//"Skipping a col", a_col,&
                !     left_col_blk_start (a_col), left_col_blk_end (a_col),&
                !     first_k, last_k
                CYCLE
             ELSE
                !i = MAX (first_k, left_col_blk_start (a_col_logical))
                i = blk_start
                j = MIN (last_k, blk_end)
                !clamp_k_offset = i - left_col_blk_start (a_col_logical)
                clamp_k_len = j - i + 1
                !write(*,'(A,I5,"(",I5,"-",I5,") not in ",I5,"-",I5)')&
                !     routineN//"Not skipping a col", a_col,&
                !     left_col_blk_start (a_col), left_col_blk_end (a_col),&
                !     first_k, last_k
                !write(*,*)routineN//" Clamping", clamp_k_offset, clamp_k_len
             ENDIF
          ENDIF
          b_row_logical = a_col_logical
          b_row = b_row_logical
          c_row_logical = a_row_logical
          a_nze = a_row_size_logical * a_col_size_logical
          b_row_size_logical = a_col_size_logical
          c_row_size_logical = a_row_size_logical
          DO b_blk = right%row_p(b_row_logical)+1, right%row_p(b_row_logical+1)
             IF (right%blk_p(b_blk) .EQ. 0) THEN
                !write(*,*)"Skipping deleted block"
                CYCLE
             ENDIF
             b_col = right%col_i(b_blk)
             c_col_logical = b_col
             ! Don't calculate symmetric blocks.
             IF (product_is_symmetric) THEN
                IF (c_row_logical .NE. c_col_logical&
                     .AND. checker_tr (c_row_logical, c_col_logical)) THEN
                   !write(*,*)"Skipping symmetric block!", c_row_logical,&
                   !     c_col_logical
                   CYCLE
                ENDIF
             ENDIF
             !IF (dbg) WRITE(*,*)routineN//" B col", b_col, ";",a_row, a_col
             b_col_logical = b_col
             b_offset = ABS (right%blk_p(b_blk))
             blk_start = right_col_blk_offset (b_col_logical)
             blk_end = right_col_blk_offset (b_col_logical+1)-1
             blk_size = blk_end - blk_start + 1
             IF (right_limit_cols) THEN
                !IF (right_col_blk_start (b_col_logical) .GT. submat_l_col .OR.&
                !    right_col_blk_end (b_col_logical) .LT. submat_f_col) THEN
                IF (blk_start .GT. submat_l_col) THEN
                   !WRITE(*,*)"Skipping column", b_col
                   CYCLE
                ELSEIF (blk_end .GT. submat_l_col) THEN
                   this_right_limit =&
                        submat_l_col - blk_start + 1
                ELSE
                   this_right_limit = 0
                ENDIF
             ELSE
                this_right_limit = 0
             ENDIF
             c_col_logical = b_col_logical
             b_col_size_logical = blk_size
             c_col_size_logical = &
                  product_col_blk_offset (c_col_logical+1) -&
                  product_col_blk_offset (c_col_logical)
             c_row_size = c_row_size_logical
             c_col_size = c_col_size_logical
             !
             b_col_size = b_col_size_logical
             b_nze = b_row_size_logical * b_col_size_logical
             tr_b = right%blk_p(b_blk) .LT. 0
             b_row_size = b_row_size_logical
             b_type = blas_mat_type (tr_b)
             IF (tr_b) CALL swap (b_col_size, b_row_size)
             !
             c_nze = c_row_size_logical * c_col_size_logical
             offset = ABS(c_blk_ps(c_col_logical))
             tr_p = c_blk_ps(c_col_logical) .LT. 0
             block_exists = c_blk_ps(c_col_logical) .NE. 0
             IF (tr_p) CALL swap (c_row_size, c_col_size)
             !DBG "Found old?", block_exists, c_blk_ps(c_col_logical), offset
             ! Try to find the current coordinate in the current row index
             sparsity_enforcement: IF &
                  (keep_sparsity .AND. .NOT. block_exists) THEN
                !use_accumulated_data = .TRUE.
                CYCLE
             ENDIF sparsity_enforcement
             IF (.NOT. block_exists) THEN
                offset = c_blk_pt
                ! CP2K can't handle transposed blocks. Blocks that are
                ! canonically under the diagonal should therefore be
                ! set at transposed.
                IF (product_is_symmetric&
                     .AND. c_row_logical .GT. c_col_logical) THEN
                   offset = -offset
                ENDIF
                new_blk = new_blk+1
                new_index(1, new_blk) = c_row_logical
                new_index(2, new_blk) = c_col_logical
                new_index(3, new_blk) = offset
                new_row_blks = new_row_blks + 1
                !DBG 'Adding coordinate', c_row_logical, c_col_logical, c_blk_pt
                c_blk_pt = c_blk_pt + c_nze
                datasize = datasize + c_nze
                lastblk = lastblk + 1
                c_blk_ps (c_col_logical) = offset
                size_c_blk_ps = size_c_blk_ps + 1
                c_blk_ps_p(size_c_blk_ps) = c_col_logical
                which_beta = beta_new
                tr_p = offset .LT. 0
             ELSE
                DBG 'Using existing coordinate for', c_row_logical, c_col_logical, offset, tr_p
                which_beta = beta_add
             ENDIF
             c_offset = ABS(offset)
             stack_pointer = stack_pointer + 1
             dgemm_parameter_stack(stack_pointer)%type_a = a_type
             dgemm_parameter_stack(stack_pointer)%type_b = b_type
             dgemm_parameter_stack(stack_pointer)%type_c = blas_mat_type(tr_p)
             dgemm_parameter_stack(stack_pointer)%lrows_a = a_row_size_logical
             dgemm_parameter_stack(stack_pointer)%lcols_a = a_col_size_logical
             dgemm_parameter_stack(stack_pointer)%lrows_b = b_row_size_logical
             dgemm_parameter_stack(stack_pointer)%lcols_b = b_col_size_logical
             dgemm_parameter_stack(stack_pointer)%lrows_c = c_row_size_logical
             dgemm_parameter_stack(stack_pointer)%lcols_c = c_col_size_logical
             dgemm_parameter_stack(stack_pointer)%offset_a = a_offset
             dgemm_parameter_stack(stack_pointer)%nze_a = a_nze
             dgemm_parameter_stack(stack_pointer)%offset_b = b_offset
             dgemm_parameter_stack(stack_pointer)%nze_b = b_nze
             dgemm_parameter_stack(stack_pointer)%offset_c = ABS(c_offset)
             dgemm_parameter_stack(stack_pointer)%nze_c = c_nze
             dgemm_parameter_stack(stack_pointer)%beta = which_beta
             IF (.NOT. limit_k) THEN
                dgemm_parameter_stack(stack_pointer)%last_k = a_col_size_logical
             ELSE
                dgemm_parameter_stack(stack_pointer)%last_k = clamp_k_len
             ENDIF
             IF (right_limit_cols .AND. this_right_limit.gt.0) THEN
                dgemm_parameter_stack(stack_pointer)%last_n = this_right_limit
             ELSE
                dgemm_parameter_stack(stack_pointer)%last_n = b_col_size_logical
                ! assert b_col_size_logical .eq. c_col_size_logical
             ENDIF
             negate_alpha = a_type .EQ. 'T' .AND. left%sym_negation
             negate_alpha = negate_alpha .NEQV. (b_type .EQ. 'T' .AND. right%sym_negation)
             negate_alpha = negate_alpha .NEQV. (tr_p .AND. product%sym_negation)
             IF (negate_alpha) THEN
                dgemm_parameter_stack(stack_pointer)%alpha &
                     = dbcsr_scalar_negative (my_alpha)
             ELSE
                dgemm_parameter_stack(stack_pointer)%alpha = my_alpha
             ENDIF
             i = UBOUND (dgemm_parameter_stack,1)
             i = MIN (i, 3*i/4)
             IF (dbg) WRITE(*,*)routineN//" Before merge, stack pointer is at",stack_pointer
             IF (stack_pointer .GE. i) THEN
                CALL merge_dgemm_stack(dgemm_parameter_stack,&
                     first_working_sp, stack_pointer)
             ENDIF
             IF (dbg) WRITE(*,*)routineN//"  After merge, stack pointer is at",stack_pointer
             IF (stack_pointer .GE. i) THEN
                IF (dbg) WRITE(*,*)routineN//" Processing stack at", stack_pointer
                mt_t_dgemm = mt_t_dgemm - (m_walltime() - epoch)
                CALL process_dgemm_stack(dgemm_parameter_stack(1:stack_pointer),&
                     stack_pointer,&
                     left%data_area, right%data_area,&
                     product%wms(ithread+1)%data_area,&
                     one_flop)
                mt_t_dgemm = mt_t_dgemm + (m_walltime() - epoch)
                lflop = lflop + one_flop
                stack_pointer = 0
                first_working_sp = 1
             ELSE
                !WRITE(*,*)routineN//" Skipping DGEMM"
             ENDIF
          ENDDO ! b
          ! Conclude the row, merge everything, and prepare for the new row.
          CALL merge_dgemm_stack(dgemm_parameter_stack,&
               first_working_sp, stack_pointer)
          first_working_sp = stack_pointer+1
       ENDDO ! a_col
    ENDDO ! a_row
    IF (stack_pointer .GT. 0) THEN
       IF (dbg) WRITE(*,*)routineN//" Final processing at", stack_pointer
       mt_t_dgemm = mt_t_dgemm - (m_walltime() - epoch)
       CALL process_dgemm_stack(dgemm_parameter_stack(1:stack_pointer),&
            stack_pointer,&
            left%data_area, right%data_area,&
            product%wms(ithread+1)%data_area,&
            one_flop)
       mt_t_dgemm = mt_t_dgemm + (m_walltime() - epoch)
       lflop = lflop + one_flop
       stack_pointer = 0
       first_working_sp = 1
    ENDIF
    !
    DEALLOCATE (c_blk_ps, c_blk_ps_p)
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
  END SUBROUTINE dbcsr_nn_mult_lin


! *****************************************************************************
!> \brief Joins compatible DGEMM parameters into a single DGEMM parameter
!>
!> \param[in,out] params      Stack of parameters
!> \param[in,out] first       First element considered joinable
!> \par first
!>      This is the first element in the list (and highest in the
!>      stack) should be considered joinable with following
!>      elements. It is used to avoid rescanning the elements at the
!>      bottom of the stack that have already been determined to be
!>      non-joinable.
!> \par[in,out] n             Number of elements (stack size)
! *****************************************************************************
  PURE SUBROUTINE merge_dgemm_stack(params, first, n)
    INTEGER, INTENT(INOUT)                   :: n, first
    TYPE(dgemm_parameters), DIMENSION(1:n), &
      INTENT(INOUT)                          :: params

    CHARACTER(len=*), PARAMETER :: routineN = 'merge_dgemm_stack', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: f, last, sp
    LOGICAL                                  :: merge_this

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY

    last = n
    f = first
    !IF (dbg) WRITE(*,*)routineN//" merge attempt",first,"to",last
    DO sp = last, f+1, -1
       !IF (dbg) THEN
       !   WRITE(*,*)routineN//" testing merge of", sp, sp-1
       !   CALL print_dgemm_parameters(params(sp))
       !   CALL print_dgemm_parameters(params(sp-1))
       !ENDIF
       merge_this = .FALSE.
       merge_this = merge_this &
            .AND. params(sp-1)%type_b .EQ. 'N'.AND. params(sp)%type_b .EQ. 'N'&
            .AND. params(sp-1)%type_c .EQ. 'N' .AND. params(sp)%type_c .EQ. 'N'

       merge_this = merge_this&
            .AND. params(sp-1)%lcols_b .EQ. params(sp-1)%last_n &
            .AND. params(sp-1)%lcols_c .EQ. params(sp-1)%last_n &
            .AND. params(sp-1)%offset_a .EQ. params(sp)%offset_a
       IF (merge_this) THEN
            merge_this = merge_this &
            .AND. dbcsr_scalar_are_equal(params(sp-1)%beta, params(sp)%beta) &
            .AND. dbcsr_scalar_are_equal(params(sp-1)%alpha, params(sp)%alpha)
       ENDIF
!       IF (params(sp-1)%type_b .EQ. 'N'.AND. params(sp)%type_b .EQ. 'N' &
!            .AND. params(sp-1)%type_c .EQ. 'N' .AND. params(sp)%type_c .EQ. 'N'&
!            .AND. dbcsr_scalar_are_equal(params(sp-1)%beta, params(sp)%beta) &
!            .AND. dbcsr_scalar_are_equal(params(sp-1)%alpha, params(sp)%alpha) &
!            .AND. params(sp-1)%lcols_b .EQ. params(sp-1)%last_n &
!            .AND. params(sp-1)%lcols_c .EQ. params(sp-1)%last_n &
!            .AND. params(sp-1)%offset_a .EQ. params(sp)%offset_a &
!            .AND. params(sp-1)%offset_b + params(sp-1)%nze_b .EQ. params(sp)%offset_b &
!            .AND. params(sp-1)%offset_c + params(sp-1)%nze_c .EQ. params(sp)%offset_c) THEN
       IF (merge_this) THEN
          !IF (dbg) WRITE(*,*)routineN//" Merging into:"
          params(sp-1)%last_n = params(sp-1)%last_n + params(sp)%last_n
          params(sp-1)%lcols_b = params(sp-1)%lcols_b + params(sp)%lcols_b
          params(sp-1)%lcols_c = params(sp-1)%lcols_c + params(sp)%lcols_c
          params(sp-1)%nze_b = params(sp-1)%nze_b + params(sp)%nze_b
          params(sp-1)%nze_c = params(sp-1)%nze_c + params(sp)%nze_c
          params(sp)%type_c = '0'
          !IF (dbg) CALL print_dgemm_parameters(params(sp-1))
          IF (n .GT. first) THEN
             n = n - 1
          ENDIF
       ELSE
          first = MAX(first, sp)
          !IF (dbg) WRITE(*,*)routineN//" New first is ",first
       ENDIF
    ENDDO
    !IF (dbg) WRITE(*,*)routineN//" The new first and n are ",first, n
  END SUBROUTINE merge_dgemm_stack


! *****************************************************************************
!> \brief Issues actual DGEMM calls.
!>
!> \param[in] params           Stack of DGEMM parameters
!> \param[in] n                Number of parameters
!> \param[in] left_data_area   Left-matrix data
!> \param[in] right_data_area  Right-matrix data
!> \param[in,out] product_data_area  Data for results
!> \param[out] lflop           (optional) Number of FLOPs used by DGEMM
! *****************************************************************************
  SUBROUTINE process_dgemm_stack(params, n,&
       left_data_area, right_data_area, product_data_area, lflop)
    INTEGER, INTENT(IN)                      :: n
    TYPE(dgemm_parameters), DIMENSION(1:n), &
      INTENT(IN)                             :: params
    TYPE(dbcsr_data_obj), INTENT(IN)         :: left_data_area, &
                                                right_data_area
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: product_data_area
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: lflop

    CHARACTER(len=*), PARAMETER :: routineN = 'process_dgemm_stack', &
      routineP = moduleN//':'//routineN
    REAL, PARAMETER                          :: resize_factor = 1.618034

    INTEGER                                  :: c, lda, ldb, ldc, maxs, r, sp
    LOGICAL                                  :: do_resize
    TYPE(dbcsr_data_obj)                     :: left_data, product_data, &
                                                right_data
    TYPE(dbcsr_scalar_type)                  :: any_zero

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY

    SELECT CASE (product_data_area%d%data_type)
    CASE (dbcsr_type_real_4)
       CALL process_dgemm_stack_s (params, n,&
       left_data_area%d%r_sp, right_data_area%d%r_sp, product_data_area,&
       lflop=lflop)
    CASE (dbcsr_type_real_8)
       CALL process_dgemm_stack_d (params, n,&
       left_data_area%d%r_dp, right_data_area%d%r_dp, product_data_area,&
       lflop=lflop)
    CASE (dbcsr_type_complex_4)
       CALL process_dgemm_stack_c (params, n,&
       left_data_area%d%c_sp, right_data_area%d%c_sp, product_data_area,&
       lflop=lflop)
    CASE (dbcsr_type_complex_8)
       CALL process_dgemm_stack_z (params, n,&
       left_data_area%d%c_dp, right_data_area%d%c_dp, product_data_area,&
       lflop=lflop)
    CASE default
       CALL cp_assert (.FALSE., cp_fatal_level, cp_caller_error,&
            routineN, "Invalid data type")
    END SELECT
    
    RETURN

    IF (PRESENT (lflop)) lflop = INT(0, int_8)
    ! Increase product data area size if necessary.
    do_resize = .FALSE.
    maxs = dbcsr_get_data_size(product_data_area)
    DO sp = 1, n
       IF (params(sp)%offset_c+params(sp)%nze_c-1 .GT. maxs) THEN
          maxs = params(sp)%offset_c+params(sp)%nze_c-1
          do_resize = .TRUE.
       ENDIF
    ENDDO
    IF (do_resize) THEN
       IF (dbg) &
            WRITE(*,*)routineN//" Data area too small, increasing!", &
            maxs, dbcsr_get_data_size(product_data_area)
       CALL dbcsr_data_ensure_size (product_data_area,&
            maxs, factor=resize_factor)
    ENDIF
    ! Setup encapsulated data area
    CALL dbcsr_data_init (left_data)
    CALL dbcsr_data_init (right_data)
    CALL dbcsr_data_init (product_data)
    CALL dbcsr_data_new (left_data, left_data_area%d%data_type)
    CALL dbcsr_data_new (right_data, right_data_area%d%data_type)
    CALL dbcsr_data_new (product_data, product_data_area%d%data_type)
    any_zero = dbcsr_scalar_zero(product_data%d%data_type)
    DO sp = 1, n
       IF (dbg) WRITE(*,*)routineN//" Stack pointer", sp
       IF (dbg) CALL print_dgemm_parameters(params(sp))
       IF (params(sp)%type_c .EQ. '0') THEN
          IF (dbg) WRITE(*,*)routineN//" Skipping", sp
          CYCLE
       ENDIF
       IF (params(sp)%nze_c .EQ. 0&
            .or.params(sp)%nze_a .EQ. 0&
            .or.params(sp)%nze_b .EQ. 0) THEN
          IF (dbg) WRITE(*,*)routineN//" Skipping empty blocks"
          CYCLE
       ENDIF

       !IF(params(sp)%offset_c+params(sp)%nze_c-1.GT.&
       !     SIZE(dbcsr_get_data_p_d(product_data_area))) THEN
          !WRITE(*,*) 'SIZE(dbcsr_get_data_p_d(product_data_area))',&
          !     SIZE(dbcsr_get_data_p_d(product_data_area))
          !WRITE(*,*) 'params(sp)%offset_c+params(sp)%nze_c-1',&
          !     params(sp)%offset_c+params(sp)%nze_c-1
          !stop
       !ENDIF
       !IF(params(sp)%offset_c.GT.&
       !     SIZE(dbcsr_get_data_p_d(product_data_area))) THEN
          !WRITE(*,*) 'SIZE(dbcsr_get_data_p_d(product_data_area))',&
          !     SIZE(dbcsr_get_data_p_d(product_data_area))
          !WRITE(*,*) 'params(sp)%offset_c',&
          !     params(sp)%offset_c
          !stop
       !ENDIF
       CALL dbcsr_data_set_pointer (product_data,&
            rsize=params(sp)%nze_c, csize=1,&
            source_lb=params(sp)%offset_c,&
            pointee=product_data_area)
       CALL dbcsr_data_set_pointer (left_data,&
            rsize=params(sp)%nze_a, csize=1,&
            source_lb=params(sp)%offset_a,&
            pointee=left_data_area)
       CALL dbcsr_data_set_pointer (right_data,&
            rsize=params(sp)%nze_b, csize=1,&
            source_lb=params(sp)%offset_b,&
            pointee=right_data_area)
       ! Keep the commented template below if needed for speed
       !product_data_d => pointer_view (dbcsr_get_data_p_d (&
       !     product_data_area),&
       !     params(sp)%offset_c,&
       !     params(sp)%offset_c+params(sp)%nze_c-1)
       !left_data_d => pointer_view (dbcsr_get_data_p_d (&
       !     left_data_area),&
       !     params(sp)%offset_a,&
       !     params(sp)%offset_a+params(sp)%nze_a-1)
       !right_data_d => pointer_view (dbcsr_get_data_p_d (&
       !     right_data_area),&
       !     params(sp)%offset_b,&
       !     params(sp)%offset_b+params(sp)%nze_b-1)
       !IF (dbg) THEN
       !   WRITE(*,*)'A'
       !   CALL dbcsr_printmat (left_data_d,&
       !        params(sp)%lrows_a, params(sp)%lcols_a,&
       !        tr=params(sp)%type_a.eq.'T')
       !   WRITE(*,*)'B'
       !   CALL dbcsr_printmat (right_data_d,&
       !        params(sp)%lrows_b, params(sp)%lcols_b,&
       !        tr=params(sp)%type_b.eq.'T')
       !   WRITE(*,*)'Cin'
       !   CALL dbcsr_printmat (product_data_d,&
       !        params(sp)%lrows_c, params(sp)%lcols_c,&
       !        tr=params(sp)%type_c.eq.'T')
       !ENDIF
       ! Avoid nagging compilers
       IF (dbcsr_scalar_are_equal(params(sp)%beta, any_zero)) THEN
          IF (params(sp)%type_c .EQ. 'N'&
               .OR. params(sp)%last_n .EQ. params(sp)%lcols_c) THEN
             CALL dbcsr_data_clear (product_data,&
                  ub=params(sp)%lrows_c*params(sp)%last_n)
             !product_data_d(1&
             !     :params(sp)%lrows_c*params(sp)%last_n) = 0.0_dp
          ELSE
             ! It must be transposed.
             !### optimize out the inner loop (turn into a range)
             DO r = 1, params(sp)%lrows_c
                DO c = 1, params(sp)%last_n
                   CALL dbcsr_data_clear(product_data,&
                        lb=(r-1)*params(sp)%lcols_c+c,&
                        ub=(r-1)*params(sp)%lcols_c+c)
                ENDDO
             ENDDO
             !FORALL (r = 1:params(sp)%lrows_c, c = 1:params(sp)%last_n)
             !   product_data_d((r-1)*params(sp)%lcols_c+c) = 0.0_dp
             !END FORALL
          ENDIF
       ENDIF
       IF (params(sp)%type_c .EQ. 'N') THEN
          IF (dbg) &
               WRITE(*,'("DGN",A,1X,A,A,2X,3(I5,1X),F5.1,1X,2(I9,1X,I5,1X),F5.1,1X,I9,1X,I5)')&
               params(sp)%type_c,&
               params(sp)%type_a, params(sp)%type_b,&
               params(sp)%lrows_c, params(sp)%last_n,& !m, n
               params(sp)%last_k,& ! k
               params(sp)%alpha,&
               params(sp)%offset_a, params(sp)%nze_a,&
               params(sp)%offset_b, params(sp)%nze_b,&
               params(sp)%beta,&
               params(sp)%offset_c, params(sp)%nze_c
          lda = select_n_or_t (params(sp)%type_a,&
               params(sp)%lrows_a, params(sp)%lcols_a)
          ldb = select_n_or_t (params(sp)%type_b,&
               params(sp)%lrows_b, params(sp)%lcols_b)
          ldc = select_n_or_t (params(sp)%type_c,&
               params(sp)%lrows_c, params(sp)%lcols_c)
          CALL dbcsr_blas_gemm(&
               params(sp)%type_a, params(sp)%type_b,&
               params(sp)%lrows_c, params(sp)%last_n,& !m, n
               params(sp)%last_k,& ! k
               params(sp)%alpha,&
               left_data, lda,&
               right_data, ldb,&
               params(sp)%beta,&
               product_data, ldc)
          !CALL DGEMM(params(sp)%type_a, params(sp)%type_b,&
          !     params(sp)%lrows_c, params(sp)%last_n,& !m, n
          !     params(sp)%last_k,& ! k
          !     params(sp)%alpha,&
          !     left_data_d, lda,&
          !     right_data_d, ldb,&
          !     params(sp)%beta,&
          !     product_data_d, ldc)
          IF (PRESENT (lflop)) &
               lflop = lflop + INT(2, int_8) * &
               INT(params(sp)%lrows_c &
               * params(sp)%last_n * params(sp)%last_k, int_8)
       ELSE
          IF (dbg) &
               WRITE(*,'("DGT",A,1X,A,A,2X,3(I5,1X),F5.1,1X,2(I9,1X,I5,1X),F5.1,1X,I9,1X,I5)')&
               params(sp)%type_c,&
               flip_type(params(sp)%type_a),&
               flip_type(params(sp)%type_b),&
               params(sp)%lcols_c, params(sp)%lrows_c,& !m, n (switched)
               params(sp)%last_k,& ! k
               params(sp)%alpha,&
               params(sp)%offset_a, params(sp)%nze_a,&
               params(sp)%offset_b, params(sp)%nze_b,&
               params(sp)%beta,&
               params(sp)%offset_c, params(sp)%nze_c
          lda = select_n_or_t (params(sp)%type_a,&
               params(sp)%lrows_a, params(sp)%lcols_a)
          ldb = select_n_or_t (params(sp)%type_b,&
               params(sp)%lrows_b, params(sp)%lcols_b)
          ldc = select_n_or_t (params(sp)%type_c,&
               params(sp)%lrows_c, params(sp)%lcols_c)
          CALL dbcsr_blas_GEMM(flip_type(params(sp)%type_b),&
               flip_type(params(sp)%type_a),&
               params(sp)%lcols_c, params(sp)%lrows_c,& !m, n (switched)
               params(sp)%last_k,& ! k
               params(sp)%alpha,&
               right_data, ldb,&
               left_data, lda,&
               params(sp)%beta,&
               product_data, ldc)
          !CALL DGEMM(flip_type(params(sp)%type_b),&
          !     flip_type(params(sp)%type_a),&
          !     params(sp)%lcols_c, params(sp)%lrows_c,& !m, n (switched)
          !     params(sp)%last_k,& ! k
          !     params(sp)%alpha,&
          !     right_data_d, ldb,&
          !     left_data_d, lda,&
          !     params(sp)%beta,&
          !     product_data_d, ldc)
          IF (PRESENT (lflop)) &
               lflop = lflop + INT(2, int_8) * &
               INT(params(sp)%last_n &
               * params(sp)%lcols_c * params(sp)%last_k, int_8)
       ENDIF
       !IF (dbg) THEN
       !   WRITE(*,*)'Cout'
       !   CALL dbcsr_printmat (product_data_d,&
       !        params(sp)%lrows_c, params(sp)%lcols_c,&
       !        tr=params(sp)%type_c.eq.'T')
       !ENDIF
    ENDDO
    CALL dbcsr_data_clear_pointer (product_data)
    CALL dbcsr_data_clear_pointer (left_data)
    CALL dbcsr_data_clear_pointer (right_data)
    CALL dbcsr_data_release (product_data)
    CALL dbcsr_data_release (left_data)
    CALL dbcsr_data_release (right_data)
  END SUBROUTINE process_dgemm_stack

  SUBROUTINE print_dgemm_parameters(params)
    TYPE(dgemm_parameters), INTENT(in)       :: params

    SELECT CASE (params%alpha%data_type)
    CASE (dbcsr_type_real_4)
       WRITE(*,'(1X,A,1X,A,1X,A,A,3(1X,2(1X,I4),";"),".",2(1X,I4),2(1X,F5.1),/,10X,"@",3(1X,2(1X,I7),";"))')&
            "DGEMM PARAMETERS",&
            params%type_c,&
            params%type_a, params%type_b,&
            params%lrows_a, params%lcols_a,&
            params%lrows_b, params%lcols_b,&
            params%lrows_c, params%lcols_c,&
            params%last_k, params%last_n,&
            params%alpha%r_sp, params%beta%r_sp,&
            params%offset_a, params%nze_a,&
            params%offset_b, params%nze_b,&
            params%offset_c, params%nze_c
    CASE (dbcsr_type_real_8)
       WRITE(*,'(1X,A,1X,A,1X,A,A,3(1X,2(1X,I4),";"),".",2(1X,I4),2(1X,F5.1),/,10X,"@",3(1X,2(1X,I7),";"))')&
            "DGEMM PARAMETERS",&
            params%type_c,&
            params%type_a, params%type_b,&
            params%lrows_a, params%lcols_a,&
            params%lrows_b, params%lcols_b,&
            params%lrows_c, params%lcols_c,&
            params%last_k, params%last_n,&
            params%alpha%r_dp, params%beta%r_dp,&
            params%offset_a, params%nze_a,&
            params%offset_b, params%nze_b,&
            params%offset_c, params%nze_c
    CASE (dbcsr_type_complex_4)
       WRITE(*,'(1X,A,1X,A,1X,A,A,3(1X,2(1X,I4),";"),".",2(1X,I4),2(1X,F5.1),/,10X,"@",3(1X,2(1X,I7),";"))')&
            "DGEMM PARAMETERS",&
            params%type_c,&
            params%type_a, params%type_b,&
            params%lrows_a, params%lcols_a,&
            params%lrows_b, params%lcols_b,&
            params%lrows_c, params%lcols_c,&
            params%last_k, params%last_n,&
            REAL(params%alpha%c_sp), REAL(params%beta%c_sp),&
            params%offset_a, params%nze_a,&
            params%offset_b, params%nze_b,&
            params%offset_c, params%nze_c
    CASE (dbcsr_type_complex_8)
       WRITE(*,'(1X,A,1X,A,1X,A,A,3(1X,2(1X,I4),";"),".",2(1X,I4),2(1X,F5.1),/,10X,"@",3(1X,2(1X,I7),";"))')&
            "DGEMM PARAMETERS",&
            params%type_c,&
            params%type_a, params%type_b,&
            params%lrows_a, params%lcols_a,&
            params%lrows_b, params%lcols_b,&
            params%lrows_c, params%lcols_c,&
            params%last_k, params%last_n,&
            REAL(params%alpha%c_dp), REAL(params%beta%c_dp),&
            params%offset_a, params%nze_a,&
            params%offset_b, params%nze_b,&
            params%offset_c, params%nze_c
    END SELECT
   END SUBROUTINE print_dgemm_parameters



  ELEMENTAL FUNCTION blas_mat_type (t)
    LOGICAL, INTENT(IN)                      :: t
    CHARACTER                                :: blas_mat_type

    IF (t) THEN
       blas_mat_type = 'T'
    ELSE
       blas_mat_type = 'N'
    ENDIF
  END FUNCTION blas_mat_type


  ELEMENTAL FUNCTION flip_type (t)
    CHARACTER, INTENT(IN)                    :: t
    CHARACTER                                :: flip_type

    SELECT CASE (t)
    CASE ('N')
       flip_type = 'T'
    CASE ('T')
       flip_type = 'N'
    CASE DEFAULT
       flip_type = '@'
    END SELECT
  END FUNCTION flip_type

  ELEMENTAL FUNCTION select_n_or_t (t, n1, n2) RESULT (val)
    CHARACTER, INTENT(in)                    :: t
    INTEGER, INTENT(in)                      :: n1, n2
    INTEGER                                  :: val

    IF (t .EQ. 'N') THEN
       val = n1
    ELSE
       val = n2
    ENDIF
  END FUNCTION select_n_or_t


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
    INTEGER, DIMENSION(3, add_size), &
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
  END SUBROUTINE merge_arrays


!! *****************************************************************************
!!> \brief Adds two DBCSR work matrices
!!> \param[in,out] wm          resulting DBCSR sum work matrix
!!> \param[in] addee           DBCSR work matrix to add to the wm matrix
!!> \param[out] flop           (optional) number of effective double-precision
!!>                            floating point operations performed
!! *****************************************************************************
!  SUBROUTINE dbcsr_wm_add(wm, addee, error, flop)
!    TYPE(dbcsr_work_type), INTENT(INOUT)     :: wm
!    TYPE(dbcsr_obj), INTENT(IN)              :: addee
!    TYPE(cp_error_type), INTENT(INOUT)       :: error
!    INTEGER, INTENT(OUT), OPTIONAL           :: flop
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_wm_add', &
!      routineP = moduleN//':'//routineN
!    REAL(KIND=dp), PARAMETER                 :: alpha = 1.0_dp, &
!                                                beta_add = 1.0_dp, &
!                                                beta_new = 0.0_dp
!
!    CHARACTER                                :: addee_type, sum_type
!    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: addee_data_z
!    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: addee_data_c
!    INTEGER :: addee_blk, addee_col, addee_col_size, addee_nze, addee_row, &
!      addee_row_size, ldaddee, ldsum, lflop, sum_blk, sum_blk_p, &
!      sum_col_size, sum_row_size, timing_handle
!    INTEGER, DIMENSION(:), POINTER           :: addee_col_blk_size, &
!                                                addee_row_blk_size
!    REAL(KIND=dp)                            :: tstart, tstop
!    REAL(KIND=dp), DIMENSION(:), POINTER     :: addee_data, addee_data_d, &
!                                                sum_data
!    REAL(KIND=sp), DIMENSION(:), POINTER     :: addee_data_r
!
!    DEBUG_HEADER
!!   ---------------------------------------------------------------------------
!    DEBUG_BODY
!    CALL timeset(routineN, timing_handle)
!    IF (dbg) tstart = m_walltime ()
!    DBG "Doing WM matrix sum",' += ',addee%m%name
!    DBG "Pre-mult blk sizes, data", wm%lastblk, wm%datasize
!    DBG "Addee blk sizes, data", addee%m%nblks, addee%m%nze
!    addee_row_blk_size => array_data (addee%m%row_blk_size)
!    addee_col_blk_size => array_data (addee%m%col_blk_size)
!    SELECT CASE (addee%m%data_type)
!    CASE (dbcsr_type_real_4)
!       CALL dbcsr_get_data (addee%m%data_area, addee_data_r)
!    CASE (dbcsr_type_real_8)
!       CALL dbcsr_get_data (addee%m%data_area, addee_data_d)
!    CASE (dbcsr_type_complex_4)
!       CALL dbcsr_get_data (addee%m%data_area, addee_data_c)
!    CASE (dbcsr_type_complex_8)
!       CALL dbcsr_get_data (addee%m%data_area, addee_data_z)
!    END SELECT
!    !
!    ! The algorithm is a follows: go through all of addee's
!    ! blocks. For every block, see if a corresponding block already
!    ! exists in wm. If such a block exists, then add the two
!    ! blocks. In this case, the index does not change. If a
!    ! corresponding block does not exist, then simply add tho new
!    ! block. The data goes to the end of the data array and the index
!    ! information is appended to wm's index structure.
!    !
!    ! c_blk_ps contains the block pointers in the product's current row
!    lflop = 0
!    IF (wm%lastblk .GT. 0 .OR. addee%m%nblks .GT. 0) THEN
!    DO addee_row = 1, addee%m%nblkrows_total
!       DBG 'Addee row', addee_row,' with #cols=',&
!            addee%m%row_p(addee_row+1)-addee%m%row_p(addee_row)
!       ! Find the blocks present in addee_row.
!       addee_row_size = addee_row_blk_size(addee_row)
!       CALL dbcsr_data_ensure_size(wm%data_area,&
!            wm%datasize+addee_row_size*addee%m%nfullcols_total &
!            * MAX(1, (addee%m%nblkrows_total-addee_row)/4))
!       !CALL ensure_array_size(wm%data,&
!       !     ub=wm%datasize+addee_row_size*addee%m%nfullcols_total &
!       !     * MAX(1, (addee%m%nblkrows_total-addee_row)/4), error=error)
!       DO addee_blk = addee%m%row_p(addee_row)+1, addee%m%row_p(addee_row+1)
!          addee_col = addee%m%col_i(addee_blk)
!          DBG 'Addee col', addee_col
!          addee_col_size = addee_col_blk_size(addee_col)
!          addee_nze = addee_row_size * addee_col_size
!          addee_data_d => addee_data_d(ABS(addee%m%blk_p(addee_blk))&
!               :ABS(addee%m%blk_p(addee_blk))+addee_nze-1)
!          !addee_data => addee%m%data(ABS(addee%m%blk_p(addee_blk))&
!          !                   :ABS(addee%m%blk_p(addee_blk))+addee_nze-1)
!          IF (addee%m%blk_p(addee_blk) .LT. 0) THEN
!             addee_type = 'T'; ldaddee = addee_col_size
!          ELSE
!             addee_type = 'N'; ldaddee = addee_row_size
!          ENDIF
!          ! Attempt to find co-lying block in the sum wm matrix.
!          sum_blk_p = 0 ! Signifies colier not present.
!          find_sum_colier: DO sum_blk = 1, wm%lastblk
!             IF (wm%row_i(sum_blk) .EQ. addee_row .AND.&
!                  wm%col_i(sum_blk) .EQ. addee_col) THEN
!                sum_blk_p = wm%blk_p(sum_blk)
!                EXIT find_sum_colier
!             ENDIF
!          ENDDO find_sum_colier
!          ! If the colier is present, then sum the data, otherwise
!          ! augment the wm sum matrix with the data and coordinates.
!          IF (sum_blk_p .GT. 0) THEN
!             sum_data => pointer_view (dbcsr_get_data_p (wm%data_area, 0.0_dp),&
!                  ABS(sum_blk_p), ABS(sum_blk_p)+addee_nze-1)
!             !sum_data => wm%data(&
!             !     ABS(sum_blk_p):ABS(sum_blk_p)+addee_nze-1)
!             IF (wm%blk_p(sum_blk) .LT. 0) THEN
!                sum_type = 'T'; ldsum = sum_col_size
!             ELSE
!                sum_type = 'N'; ldsum = sum_row_size
!             ENDIF
!             IF (sum_type .EQ. addee_type) THEN
!                CALL DAXPY (addee_nze, 1.0_dp, addee_data, 1, sum_data, 1)
!             ELSE
!                CALL DAXPY (addee_nze, 1.0_dp,&
!temp_transpose(addee_data,addee_row_size,addee_col_size),&
!                     1, sum_data, 1)
!             ENDIF
!             lflop = lflop + addee_nze
!          ELSE
!             sum_blk_p = wm%datasize + 1
!             sum_data => pointer_view (dbcsr_get_data_p (wm%data_area, 0.0_dp),&
!                  ABS(sum_blk_p), ABS(sum_blk_p)+addee_nze-1)
!             !sum_data => wm%data(&
!             !     ABS(sum_blk_p):ABS(sum_blk_p)+addee_nze-1)
!             sum_data = addee_data
!             CALL add_work_coordinate(wm, addee_row, addee_col, sum_blk_p,&
!                  error=error)
!             wm%datasize = wm%datasize + addee_nze
!          END IF
!       ENDDO ! addee_blk
!    ENDDO ! addee_row
!    ELSE
!       DBG 'Both matrices were empty!'
!    ENDIF
!    DBG "Post-add blk sizes, data", wm%lastblk, wm%datasize
!    IF (dbg .AND. detailed_timing) THEN
!       tstop = m_walltime ()
!       IF (tstop-tstart .NE. 0) THEN
!          DBG tstop-tstart,'s',&
!               (REAL(lflop, dp)/1000000.0_dp)/(tstop-tstart),'MFLOPS'
!       ENDIF
!    ENDIF
!    IF (PRESENT (flop)) flop = lflop
!    CALL timestop(timing_handle)
!  END SUBROUTINE dbcsr_wm_add


! *****************************************************************************
!> \brief Scales a DBCSR matrix by alpha
!> \param[inout] matrix_a       DBCSR matrix
!> \param[in] alpha_scalar      (optional) a scalar
!> \param[in] alpha_matrix      (optional) a list of matrices
!> \param[in] side              (optional) ?
!> \param[in] last_column       (optional) only scale up to this full column
!> \param[in] scale type        corresponds to scalar, vector, matrix
! *****************************************************************************
  SUBROUTINE dbcsr_scale_anytype(matrix_a, last_column, error, scale_type,&
         side,&
         alpha_scalar, alpha_matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    CHARACTER, INTENT(IN)                    :: scale_type
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: side
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: alpha_scalar
    TYPE(dbcsr_data_obj), INTENT(IN), &
      OPTIONAL                               :: alpha_matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_anytype', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER :: a_blk, a_col, a_col_size, a_nze, a_row, a_row_size, &
      col_offset, i, last_col, row_offset, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: col_end, col_start, m_offset
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: do_scale, limit_col, right, tr
    TYPE(dbcsr_data_obj)                     :: buffer, data_any, mat_tmp
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    !
    ! last_column => alpha_scalar
    CALL cp_assert (.NOT. PRESENT (last_column) .OR. scale_type .EQ. 'S', &
         cp_fatal_level, cp_wrong_args_error, routineN,&
         "Column limiting supported only with scalar scaling")
    last_col = 0
    IF (PRESENT (last_column)) THEN
       CALL cp_assert (last_column .GE. 0 .AND. last_column .LE.&
            dbcsr_nfullcols_total (matrix_a), cp_fatal_level,&
            cp_wrong_args_error, routineN,&
            "Specified last column is out of bounds.")
       last_col = last_column
       ALLOCATE (col_start(dbcsr_nblkcols_total(matrix_a)))
       ALLOCATE (col_end(dbcsr_nblkcols_total(matrix_a)))
       CALL convert_sizes_to_offsets(array_data (matrix_a%m%col_blk_size), &
            col_start, col_end)
    ENDIF
    data_type = dbcsr_get_data_type (matrix_a)
    !
    ! set vars
    right = .TRUE.
    IF(PRESENT(side)) THEN
       SELECT CASE(side)
       CASE('right');right = .TRUE.
       CASE('left' );right = .FALSE.
       CASE DEFAULT
          CALL cp_assert (.FALSE., cp_fatal_level,&
               cp_wrong_args_error, routineN,&
               "wrong side="//side)
       END SELECT
    ENDIF
    IF(scale_type .EQ. 'M') THEN     ! Matrix scaling
       CALL dbcsr_data_init (buffer)
       ! Clean up these maxvalues
       CALL dbcsr_data_new (buffer, dbcsr_get_data_type(matrix_a),&
            data_size=MAXVAL(array_data (matrix_a%m%row_blk_size)) &
            *MAXVAL(array_data (matrix_a%m%col_blk_size)))
       CALL dbcsr_data_clear (buffer) ! Zero the buffer
       !
       ! TODO: Clean up these offsets using the matrix-supplied ones
       row_blk_size => array_data (matrix_a%m%row_blk_size)
       col_blk_size => array_data (matrix_a%m%col_blk_size)
       ALLOCATE(m_offset(dbcsr_nblkrows_total(matrix_a)))
       m_offset(1) = 1
       IF(right) THEN
          DO a_col = 2,dbcsr_nblkcols_total(matrix_a)
             a_col_size = col_blk_size(a_col-1)
             m_offset(a_col) = m_offset(a_col-1) + a_col_size**2
          ENDDO
       ELSE
          DO a_row = 2,dbcsr_nblkrows_total(matrix_a)
             a_row_size = row_blk_size(a_row-1)
             m_offset(a_row) = m_offset(a_row-1) + a_row_size**2
          ENDDO
       ENDIF
       !
       ! This data area holds the current column
       CALL dbcsr_data_init (mat_tmp)
       CALL dbcsr_data_new (mat_tmp, dbcsr_get_data_type (alpha_matrix))
    ENDIF
    !
    ! quick return if possible
    do_scale = .TRUE.
    IF (PRESENT (alpha_scalar)) THEN
       SELECT CASE (alpha_scalar%data_type)
       CASE (dbcsr_type_real_4)
          IF (alpha_scalar%r_sp .EQ. 1.0_real_4) do_scale = .FALSE.
       CASE (dbcsr_type_real_8)
          IF (alpha_scalar%r_dp .EQ. 1.0_real_8) do_scale = .FALSE.
       CASE (dbcsr_type_complex_4)
          IF (alpha_scalar%c_sp .EQ. CMPLX(1.0, 0.0, real_4)) do_scale = .FALSE.
       CASE (dbcsr_type_complex_8)
          IF (alpha_scalar%c_dp .EQ. CMPLX(1.0, 0.0, real_8)) do_scale = .FALSE.
       END SELECT
    END IF
    !
    ! let's go
    IF(do_scale) THEN
       CALL dbcsr_data_init (data_any)
       CALL dbcsr_data_new (data_any, dbcsr_get_data_type (matrix_a))
       CALL dbcsr_iterator_start(iter, matrix_a)
       iterations: DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block(iter, a_row, a_col, data_any, tr,&
               block_number=a_blk,&
               row_size=a_row_size, col_size=a_col_size, &
               row_offset=row_offset, col_offset=col_offset)
          ! ### UB: Find where the transposed blocks come from.
          !CALL cp_assert (.NOT. tr, cp_fatal_level, cp_unimplemented_error_nr,&
          !     routineN, "No support for transposed blocks.")
          a_nze = a_row_size * a_col_size
          IF (a_nze .EQ. 0) CYCLE ! Skip empty blocks
          !
          ! let's scale
          IF(scale_type .EQ. 'S') THEN
             limit_col = last_col .GT. 0
             IF (limit_col) limit_col = col_end(a_col).gt.last_col
             IF (limit_col) THEN
                IF (col_start(a_col) .LE. last_col) THEN
                   i = MIN(last_col, col_end(a_col)) - col_start(a_col) + 1
                   a_nze = a_row_size*i
                ELSE
                   ! Nothing to do, the block is beyond the boundry.
                   a_nze = 0
                ENDIF
             ENDIF
             IF (a_nze .GT. 0) THEN
                CALL dbcsr_blas_scal (a_nze, alpha_scalar, data_any, 1)
             ENDIF
          ELSEIF(scale_type .EQ. 'M') THEN
             IF(right) THEN
                !A = A * alpha
                mat_tmp = pointer_view (mat_tmp, alpha_matrix, m_offset(a_col))
                CALL dbcsr_blas_gemm ('N','N',&
                     a_row_size, a_col_size, a_col_size,&
                     dbcsr_scalar_one (alpha_matrix%d%data_type),&
                     data_any, a_row_size,&
                     mat_tmp, a_col_size,&
                     dbcsr_scalar_zero (alpha_matrix%d%data_type),&
                     buffer, a_row_size)
                !CALL dgemm('N','N',a_row_size,a_col_size,a_col_size,&
                !     1.0_dp,data_any%d%r_dp(1),a_row_size,&
                !     alpha_matrix%d%r_dp(m_offset(a_col)),a_col_size,&
                !     0.0_dp,buffer%d%r_dp(1),a_row_size)
                CALL dbcsr_blas_copy (a_nze, buffer, 1, data_any, 1)
             ELSE
                !A = alpha * A
                mat_tmp = pointer_view (mat_tmp, alpha_matrix, m_offset(a_row))
                CALL dbcsr_blas_gemm ('N','N',&
                     a_row_size, a_col_size, a_row_size,&
                     dbcsr_scalar_one (alpha_matrix%d%data_type),&
                     mat_tmp, a_row_size,&
                     data_any, a_row_size,&
                     dbcsr_scalar_zero (alpha_matrix%d%data_type),&
                     buffer, a_row_size)
                !CALL dgemm('N','N',a_row_size,a_col_size,a_row_size,&
                !        &     1.0_dp,alpha_matrix%d%r_dp(m_offset(a_col)),a_row_size,&
                !        &            data_any%d%r_dp(1),a_row_size,&
                !        &     0.0_dp,buffer%d%r_dp(1),a_row_size)
                CALL dbcsr_blas_copy(a_nze, buffer, 1, data_any, 1)
             ENDIF
          ELSE
             CALL cp_assert (.FALSE., cp_fatal_level,&
                  cp_wrong_args_error, routineN,&
                  "wrong logic")
          ENDIF
       ENDDO iterations
       CALL dbcsr_iterator_stop(iter)
       CALL dbcsr_data_clear_pointer (data_any)
       CALL dbcsr_data_release (data_any)
    ENDIF

    IF(scale_type .EQ. 'M') THEN
       DEALLOCATE(m_offset)
       !@@@
       WRITE(*,*)"refs:", buffer%d%refcount, buffer%d%data_type
       CALL dbcsr_data_release (buffer)
       CALL dbcsr_data_clear_pointer (mat_tmp)
       CALL dbcsr_data_release (mat_tmp)
    ENDIF

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_scale_anytype

  !> \brief Interface for matrix scaling by a scalar
  SUBROUTINE dbcsr_scale_d(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(real_8), INTENT(IN)                 :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_d', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: scale_type

    scale_type = 'S'
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_8) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (alpha_scalar))
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (REAL(alpha_scalar,real_4)))
    ELSE
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error, routineN,&
            "Invalid data type")
    ENDIF
  END SUBROUTINE dbcsr_scale_d

  SUBROUTINE dbcsr_scale_s(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(real_4), INTENT(IN)                 :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_s', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: scale_type

    scale_type = 'S'
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (alpha_scalar))
    ELSE
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error, routineN,&
            "Invalid data type")
    ENDIF
  END SUBROUTINE dbcsr_scale_s
  
SUBROUTINE dbcsr_scale_z(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    COMPLEX(real_8), INTENT(IN)              :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_z', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: scale_type

    scale_type = 'S'
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_8) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (alpha_scalar))
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_4) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (CMPLX(alpha_scalar,KIND=real_4)))
    ELSE
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error, routineN,&
            "Invalid data type")
    ENDIF
  END SUBROUTINE dbcsr_scale_z

  SUBROUTINE dbcsr_scale_c(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    COMPLEX(real_4), INTENT(IN)              :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_c', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: scale_type

    scale_type = 'S'
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_4) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (alpha_scalar))
    ELSE
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error, routineN,&
            "Invalid data type")
    ENDIF
  END SUBROUTINE dbcsr_scale_c

! *****************************************************************************
!> \brief Scales a DBCSR matrix by alpha
!> \param[inout] matrix_a       DBCSR matrix
!> \param[in] alpha             the scaling vector
!> \param[in] side              apply the scaling from the side
! *****************************************************************************
  SUBROUTINE dbcsr_scale_by_vector_anytype(matrix_a, alpha, side, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_data_obj), INTENT(IN), &
      OPTIONAL                               :: alpha
    CHARACTER(LEN=*), INTENT(IN)             :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_scale_by_vector_anytype', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER                                  :: a_blk, a_col, a_col_size, &
                                                a_nze, a_row, a_row_size, &
                                                col_offset, i, row_offset, &
                                                timing_handle
    LOGICAL                                  :: right, tr
    TYPE(dbcsr_data_obj)                     :: data_any
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    !
    ! check that alpha and matrix have the same data type
    CALL cp_assert (dbcsr_get_data_type (matrix_a).EQ.alpha%d%data_type, cp_fatal_level,&
         cp_wrong_args_error, routineN, "wrong data type matrix_a="//&
         dbcsr_get_data_type (matrix_a)//' alpha='//alpha%d%data_type)
    !
    ! set vars
    right = .TRUE.
    SELECT CASE(side)
    CASE('right');right = .TRUE.
    CASE('left' );right = .FALSE.
    CASE DEFAULT
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error, &
            routineN, "wrong side="//side)
    END SELECT
    !
    ! let's go
    data_type = dbcsr_get_data_type (matrix_a)
    CALL dbcsr_data_init(data_any)
    CALL dbcsr_data_new(data_any, dbcsr_get_data_type (matrix_a))
    CALL dbcsr_iterator_start(iter, matrix_a)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block(iter, a_row, a_col, data_any, tr,&
            block_number=a_blk,&
            row_size=a_row_size, col_size=a_col_size, &
            row_offset=row_offset, col_offset=col_offset)
       a_nze = a_row_size * a_col_size
       IF (a_nze .EQ. 0) CYCLE ! Skip empty blocks
       !
       ! let's scale
       IF(right) THEN
          DO i = 1,a_col_size
             SELECT CASE (data_type)
             CASE (dbcsr_type_real_4)
                CALL sscal(a_row_size,alpha%d%r_sp(col_offset+i-1),&
                     data_any%d%r_sp((i-1)*a_row_size+1),1)
             CASE (dbcsr_type_real_8)
                CALL dscal(a_row_size,alpha%d%r_dp(col_offset+i-1),&
                     data_any%d%r_dp((i-1)*a_row_size+1),1)
             CASE (dbcsr_type_complex_4)
                CALL cscal(a_row_size,alpha%d%c_sp(col_offset+i-1),&
                     data_any%d%c_sp((i-1)*a_row_size+1),1)
             CASE (dbcsr_type_complex_8)
                CALL zscal(a_row_size,alpha%d%c_dp(col_offset+i-1),&
                     data_any%d%c_dp((i-1)*a_row_size+1),1)
             END SELECT
          ENDDO
       ELSE
          DO i = 1,a_row_size
             SELECT CASE (data_type)
             CASE (dbcsr_type_real_4)
                CALL sscal(a_col_size,alpha%d%r_sp(row_offset+i-1),&
                     data_any%d%r_sp(i),a_col_size)
             CASE (dbcsr_type_real_8)
                CALL dscal(a_col_size,alpha%d%r_dp(row_offset+i-1),&
                     data_any%d%r_dp(i),a_col_size)
             CASE (dbcsr_type_complex_4)
                CALL cscal(a_col_size,alpha%d%c_sp(row_offset+i-1),&
                     data_any%d%c_sp(i),a_col_size)
             CASE (dbcsr_type_complex_8)
                CALL zscal(a_col_size,alpha%d%c_dp(row_offset+i-1),&
                     data_any%d%c_dp(i),a_col_size)
             END SELECT
          ENDDO
       ENDIF
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_data_clear_pointer (data_any)
    CALL dbcsr_data_release (data_any)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_scale_by_vector_anytype


! *****************************************************************************
!> \brief Detransposes all blocks in a matrix
!> \param[in,out] matrix      DBCSR matrix
! *****************************************************************************
  SUBROUTINE dbcsr_make_untransposed_blocks(matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_make_untransposed_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_size, row, &
                                                row_size
    INTEGER, DIMENSION(:), POINTER           :: cbs, rbs
    LOGICAL                                  :: sym_negation, tr
    TYPE(dbcsr_data_obj)                     :: block_data
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    rbs => array_data (dbcsr_row_block_sizes (matrix))
    cbs => array_data (dbcsr_col_block_sizes (matrix))
    CALL dbcsr_data_init (block_data)
    CALL dbcsr_data_new (block_data, dbcsr_get_data_type(matrix))
    CALL dbcsr_iterator_start (iter, matrix)
    sym_negation = matrix%m%sym_negation
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, block_data,&
            transposed=tr,&
            block_number=blk)
       IF (tr) THEN
          row_size = rbs(row)
          col_size = cbs(col)
          CALL dbcsr_block_transpose(block_data, col_size, row_size)
          IF (sym_negation) THEN
             SELECT CASE (block_data%d%data_type)
             CASE (dbcsr_type_real_4)
                block_data%d%r_sp(:) = -block_data%d%r_sp(:)
             CASE (dbcsr_type_real_8)
                block_data%d%r_dp(:) = -block_data%d%r_dp(:)
             CASE (dbcsr_type_complex_4)
                block_data%d%c_sp(:) = -block_data%d%c_sp(:)
             CASE (dbcsr_type_complex_8)
                block_data%d%c_dp(:) = -block_data%d%c_dp(:)
             END SELECT
          ENDIF
          matrix%m%blk_p(blk) = -matrix%m%blk_p(blk)
       ENDIF
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    CALL dbcsr_data_clear_pointer (block_data)
    CALL dbcsr_data_release (block_data)
  END SUBROUTINE dbcsr_make_untransposed_blocks


! *****************************************************************************
!> \brief Set a DBCSR matrix by alpha
!> \param[inout] matrix       DBCSR matrix
!> \param[in] alpha           a scalar
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_set_anytype(matrix, alpha, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_anytype', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, row, timing_handle
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_block
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    !
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_new (data_block, dbcsr_get_data_type (matrix))
    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_block, tr, blk)
       CALL dbcsr_data_clear (data_block, value=alpha)
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_data_clear_pointer (data_block)
    CALL dbcsr_data_release (data_block)
    !
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_set_anytype

! *****************************************************************************
!> \brief Conjugate a DBCSR matrix
!> \param[inout] matrix       DBCSR matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_conjg(matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_conjg', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER                                  :: blk, col, row, timing_handle
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_any
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    !
    data_type = dbcsr_get_data_type(matrix)
    CALL dbcsr_data_init (data_any)
    CALL dbcsr_data_new (data_any, data_type)
    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_any, tr, blk)
       SELECT CASE (data_type)
       CASE (dbcsr_type_complex_4)
          data_any%d%c_sp = CONJG(data_any%d%c_sp)
       CASE (dbcsr_type_complex_8)
          data_any%d%c_dp = CONJG(data_any%d%c_dp)
       CASE DEFAULT
          ! needed for g95
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_data_clear_pointer(data_any)
    CALL dbcsr_data_release(data_any)
    !
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_conjg

! *****************************************************************************
!> \brief Marks a block for removal from a DBCSR matrix. Handles
!>        symmetric matrices.
!> \param[in]  matrix         DBCSR matrix
!> \param[in]  row            row of block to remove
!> \param[in]  col            column of block to remove
!> \param[in]  block_number   (optional) the block number, if it is known
! *****************************************************************************
  SUBROUTINE dbcsr_remove_block(matrix, row, col, block_nze, block_number)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col, block_nze
    INTEGER, INTENT(IN), OPTIONAL            :: block_number

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_remove_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: b, c, r
    LOGICAL                                  :: found, tr

!   ---------------------------------------------------------------------------

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
       ! update nze accordingly
       matrix%m%nze = matrix%m%nze - block_nze
       CALL cp_assert (matrix%m%nze .GE. 0, cp_failure_level,&
            cp_caller_error, routineN, "nze < 0!")
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
  SUBROUTINE dbcsr_add_anytype(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_anytype', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type_a, data_type_b
    INTEGER                                  :: blk, col, row, size_a, &
                                                size_b, timing_handle
    LOGICAL                                  :: do_scale, tr
    TYPE(dbcsr_data_obj)                     :: data_block
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_scalar_type)                  :: my_alpha_scalar, &
                                                my_beta_scalar

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    data_type_a = dbcsr_get_data_type(matrix_a)
    data_type_b = dbcsr_get_data_type(matrix_b)
    !
    my_alpha_scalar = dbcsr_scalar_one (data_type_a)
    IF(PRESENT(alpha_scalar)) my_alpha_scalar = alpha_scalar
    my_beta_scalar = dbcsr_scalar_one (data_type_b)
    IF(PRESENT(beta_scalar)) my_beta_scalar = beta_scalar
    !
    ! let's go
    CALL cp_assert (dbcsr_nblkrows_total(matrix_a).EQ.dbcsr_nblkrows_total(matrix_b), &
         cp_fatal_level, cp_wrong_args_error, routineN, "matrices not consistent")

    do_scale = dbcsr_scalar_are_equal (&
         my_beta_scalar, dbcsr_scalar_one (data_type_b))

    CALL dbcsr_scale(matrix_a,scale_type='S',&
         alpha_scalar=my_alpha_scalar,error=error)

    ! Pre-size work arrays of matrix_a to avoid continuous reallocation.
    size_a = dbcsr_get_data_size (matrix_a)
    size_b = dbcsr_get_data_size (matrix_b)
    IF (size_b .GT. size_a) THEN
       CALL dbcsr_work_create (matrix_a,&
            nblks_guess=matrix_b%m%nblks-matrix_a%m%nblks,&
            sizedata_guess=size_b-size_a)
    ENDIF

    IF(.NOT.dbcsr_scalar_are_equal(my_beta_scalar,&
         dbcsr_scalar_zero(data_type_b))) THEN
       CALL dbcsr_data_init (data_block)
       CALL dbcsr_data_new (data_block, data_type_b)
       CALL dbcsr_iterator_start(iter, matrix_b)

       DO WHILE (dbcsr_iterator_blocks_left(iter))

          CALL dbcsr_iterator_next_block(iter, row, col, data_block, tr, blk)

             IF (do_scale) THEN
                CALL dbcsr_put_block(matrix_a, row, col, data_block, tr,&
                     summation=.TRUE.)
             ELSE
                CALL dbcsr_put_block(matrix_a, row, col, data_block, tr,&
                     summation=.TRUE., scale=my_beta_scalar)
             ENDIF

       ENDDO

       CALL dbcsr_iterator_stop(iter)
       CALL dbcsr_finalize (matrix_a, error=error)
       CALL dbcsr_data_clear_pointer (data_block)
       CALL dbcsr_data_release (data_block)

    ENDIF

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_add_anytype

  !> \brief Interface for dbcsr_add
  SUBROUTINE dbcsr_add_d(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    REAL(real_8), INTENT(IN)                 :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_d', &
      routineP = moduleN//':'//routineN

    IF(    dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_8 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_8) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar),&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(REAL(alpha_scalar,real_4)),&
            beta_scalar=dbcsr_scalar(REAL(beta_scalar,real_4)), error=error)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_8 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar),&
            beta_scalar=dbcsr_scalar(REAL(beta_scalar,real_4)), error=error)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_8) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(REAL(alpha_scalar,real_4)),&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ELSE
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI")
    ENDIF
  END SUBROUTINE dbcsr_add_d

  SUBROUTINE dbcsr_add_s(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    REAL(real_4), INTENT(IN)                 :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_s', &
      routineP = moduleN//':'//routineN

    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4 .AND. &
       dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar),&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ELSE
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI")
    ENDIF
  END SUBROUTINE dbcsr_add_s

  SUBROUTINE dbcsr_add_z(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    COMPLEX(real_8), INTENT(IN)              :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_z', &
      routineP = moduleN//':'//routineN

    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_8 .AND. &
       dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_complex_8) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar),&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_4 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_complex_4) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(CMPLX(alpha_scalar,KIND=real_4)),&
            beta_scalar=dbcsr_scalar(CMPLX(beta_scalar,KIND=real_4)), error=error)
    ELSE
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI")
    ENDIF
  END SUBROUTINE dbcsr_add_z

  SUBROUTINE dbcsr_add_c(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    COMPLEX(real_4), INTENT(IN)              :: alpha_scalar, beta_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_c', &
      routineP = moduleN//':'//routineN

    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_4 .AND. &
       dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_complex_4) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar),&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ELSE
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI")
    ENDIF
  END SUBROUTINE dbcsr_add_c

! *****************************************************************************
!> \brief Hadamard product
!>    C = A . B (C needs to be different from A and B)
!> \param[in,out] matrix_c        DBCSR matrix
!> \param[in] matrix_a, matrix_b  DBCSR matrix
! *****************************************************************************
  SUBROUTINE dbcsr_hadamard_product(matrix_a, matrix_b, matrix_c, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_hadamard_product', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER                                  :: blk, col, col_size, nze, row, &
                                                row_size, timing_handle
    LOGICAL                                  :: found, tr_a, tr_b
    TYPE(dbcsr_data_obj)                     :: a_data, b_data, c_data
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL cp_assert (dbcsr_get_data_type(matrix_a).EQ.dbcsr_get_data_type(matrix_b).AND.&
         dbcsr_get_data_type(matrix_a).EQ.dbcsr_get_data_type(matrix_c), &
         cp_fatal_level, cp_wrong_args_error, routineN, &
         "data types not consistent, need to fix that")

    CALL cp_assert (dbcsr_nblkrows_total(matrix_a).EQ.dbcsr_nblkrows_total(matrix_b).AND.&
         dbcsr_nblkrows_total(matrix_c).EQ.dbcsr_nblkrows_total(matrix_a), &
         cp_fatal_level, cp_wrong_args_error, routineN, &
         "matrices not consistent")

    data_type = dbcsr_get_data_type(matrix_a)
    CALL dbcsr_data_init (c_data)
    CALL dbcsr_data_new (c_data, data_type,&
         data_size=dbcsr_max_row_size(matrix_a)*dbcsr_max_col_size(matrix_a))
    CALL dbcsr_set(matrix_c, dbcsr_scalar_zero(data_type), error=error)
    CALL dbcsr_data_init (a_data)
    CALL dbcsr_data_new (a_data, data_type)
    CALL dbcsr_data_init (b_data)
    CALL dbcsr_data_new (b_data, data_type)
   CALL dbcsr_iterator_start(iter, matrix_a)
   DO WHILE (dbcsr_iterator_blocks_left(iter))
       SELECT CASE (dbcsr_get_data_type(matrix_a))
          !CASE (dbcsr_type_real_4)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, a_data, tr_a, blk, &
               row_size=row_size, col_size=col_size)
          nze = row_size * col_size
          CALL dbcsr_get_block_p(matrix_b, row, col, b_data, tr_b, found)
          CALL cp_assert (tr_a.EQV.tr_b, cp_fatal_level, cp_wrong_args_error, routineN, &
               "tr not consistent, need to fix that")
          IF(found) THEN
             SELECT CASE (data_type)
             CASE (dbcsr_type_real_4)
                c_data%d%r_sp(1:nze) = a_data%d%r_sp(1:nze) * b_data%d%r_sp(1:nze)
             CASE (dbcsr_type_real_8)
                c_data%d%r_dp(1:nze) = a_data%d%r_dp(1:nze) * b_data%d%r_dp(1:nze)
             CASE (dbcsr_type_complex_4)
                c_data%d%c_sp(1:nze) = a_data%d%c_sp(1:nze) * b_data%d%c_sp(1:nze)
             CASE (dbcsr_type_complex_8)
                c_data%d%c_dp(1:nze) = a_data%d%c_dp(1:nze) * b_data%d%c_dp(1:nze)
             END SELECT
             CALL dbcsr_put_block(matrix_c, row, col, c_data, tr_a, &
                  summation=.FALSE.)
          ENDIF
          !CASE (dbcsr_type_complex_4)
          !CASE (dbcsr_type_complex_8)
       CASE DEFAULT
          CALL cp_unimplemented_error (routineN,"Only real double precision")
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_finalize (matrix_c, error=error)
    CALL dbcsr_data_clear_pointer (a_data)
    CALL dbcsr_data_clear_pointer (b_data)
    CALL dbcsr_data_release (c_data)
    CALL dbcsr_data_release (a_data)
    CALL dbcsr_data_release (b_data)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_hadamard_product

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

    CALL cp_assert (dbcsr_get_data_type(matrix_a).EQ.dbcsr_get_data_type(matrix_b), &
         cp_fatal_level, cp_wrong_args_error, routineN, "data types not consistent")
    !
    ! let's go
    CALL cp_assert (matrix_a%m%nblkrows_total.EQ.matrix_b%m%nblkrows_total, &
         cp_fatal_level, cp_wrong_args_error, routineN, "matrices not consistent")

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
! *****************************************************************************
  SUBROUTINE dbcsr_add_on_diag(matrix, alpha_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(KIND=real_8), INTENT(IN)            :: alpha_scalar
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_on_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: hold, mynode, n, row, &
                                                timing_handle
    INTEGER, DIMENSION(:), POINTER           :: row_blk_size
    LOGICAL                                  :: found, tr
    REAL(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: buff, data_d

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL cp_assert (dbcsr_nblkrows_total(matrix).EQ.dbcsr_nblkcols_total(matrix).AND.&
         dbcsr_nfullrows_total(matrix).EQ.dbcsr_nfullrows_total(matrix), &
         cp_fatal_level, cp_wrong_args_error, routineN, "matrices not consistent")

    row_blk_size => array_data(dbcsr_row_block_sizes(matrix))
    mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(dbcsr_distribution(matrix)))

    CALL dbcsr_work_create(matrix, work_mutable=.TRUE., error=error)

    DO row = 1,dbcsr_nblkrows_total(matrix)
       CALL dbcsr_get_stored_coordinates (matrix, row, row, tr, hold)
       IF(hold.EQ.mynode) THEN
          CALL dbcsr_get_block_p(matrix, row, row, data_d, tr, found)
          IF(found) THEN
             CALL block_add_on_diag(SIZE(data_d,1),data_d,alpha_scalar)
          ELSE
             n = row_blk_size(row)
             ALLOCATE(buff(n,n))
             buff(:,:) = 0.0_dp
             CALL block_add_on_diag(SIZE(buff,1),buff,alpha_scalar)
             CALL dbcsr_put_block(matrix, row, row, buff)
             DEALLOCATE(buff)
          ENDIF
       ENDIF
    ENDDO

    CALL dbcsr_finalize(matrix, error=error)

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

    INTEGER                                  :: blk, col, i, row, row_offset, &
                                                timing_handle
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL cp_assert (dbcsr_nfullrows_total(matrix).LE.SIZE(diag), &
         cp_fatal_level, cp_wrong_args_error, routineN, "SIZE(diag) too small")

    diag = 0.0_dp

    CALL dbcsr_iterator_start(iter, matrix)

    DO WHILE (dbcsr_iterator_blocks_left(iter))

       SELECT CASE (dbcsr_get_data_type(matrix))
          !CASE (dbcsr_type_real_4)
          !CALL dbcsr_iterator_next_block(iter, row, col, data_r, tr, blk)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk, &
               row_offset=row_offset)
          IF(row.EQ.col) THEN
             DO i = 1,SIZE(data_d,1)
                diag(row_offset+i-1) = data_d(i,i)
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

    INTEGER                                  :: blk, col, i, row, row_offset, &
                                                timing_handle
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL cp_assert (dbcsr_nfullrows_total(matrix).LE.SIZE(diag), &
         cp_fatal_level, cp_wrong_args_error, routineN, "SIZE(diag) too small")

    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       SELECT CASE (dbcsr_get_data_type(matrix))
          !CASE (dbcsr_type_real_4)
          !CALL dbcsr_iterator_next_block(iter, row, col, data_s, tr, blk)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk, &
               row_offset=row_offset)
          IF(row.EQ.col) THEN
             DO i = 1,SIZE(data_d,1)
                data_d(i,i) = diag(row_offset+i-1)
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
!> \param[inout] matrix_b       target DBCSR matrix
!> \param[in]    matrix_a       source DBCSR matrix
!> \param[in]    name           (optional) name of the new matrix
!> \param[in,out] error         cp2k error
!> \param[in]    keep_sparsity  (optional) keep the target matrix sparsity;
!>                              default is False.
!> \param[in]    shallow_data   (optional) shallow data copy
!> \param[in]    keep_imaginary  (optional) when copy from complex to real,&
!>                               the default is to keep only the real part; if
!>                               this flag is set, the imaginary part is used
! *****************************************************************************
  SUBROUTINE dbcsr_copy(matrix_b, matrix_a, name, error, keep_sparsity,&
       shallow_data, keep_imaginary)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity, shallow_data, &
                                                keep_imaginary

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: new_type
    INTEGER                                  :: timing_handle
    LOGICAL                                  :: keep_sparse, shallow

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    shallow = .FALSE. ; IF (PRESENT (shallow_data)) shallow = shallow_data
    keep_sparse = .FALSE.
    IF (PRESENT (keep_sparsity)) keep_sparse = keep_sparsity

    !CALL cp_assert (dbcsr_is_initialized (matrix_b),&
    !     cp_fatal_level, cp_wrong_args_error, routineN,&
    !     "Target matrix must be initialized", error=error)
    CALL cp_assert (.not.keep_sparse.or.dbcsr_valid_index(matrix_b),&
         cp_fatal_level, cp_wrong_args_error, routineN,&
         "Target matrix must be valid to keep its sparsity", error=error)
    CALL cp_assert (.not.keep_sparse.or..not.shallow, cp_warning_level,&
         cp_wrong_args_error, routineN,&
         "Shallow copy not compatibly with sparsity retainment", error=error)
    IF (keep_sparse) THEN
       IF (PRESENT (name)) matrix_b%m%name = name
       CALL dbcsr_copy_into_existing(matrix_b, matrix_a, error)
    ELSE
       IF (dbcsr_is_initialized (matrix_b)) THEN
          new_type = dbcsr_get_data_type (matrix_b)
       ELSE
          new_type = dbcsr_get_data_type (matrix_a)
       ENDIF
       IF(PRESENT(name)) THEN
          CALL dbcsr_create(matrix_b, name, &
               matrix_a%m%dist, matrix_a%m%matrix_type,&
               matrix_a%m%row_blk_size, matrix_a%m%col_blk_size,&
               nblks=matrix_a%m%nblks, data_type=new_type,&
               special=matrix_a%m%special_memory, error=error)
       ELSE
          CALL dbcsr_create(matrix_b,'copy of '//matrix_a%m%name, &
               matrix_a%m%dist, matrix_a%m%matrix_type,&
               matrix_a%m%row_blk_size, matrix_a%m%col_blk_size,&
               nblks=matrix_a%m%nblks, data_type=new_type,&
               special=matrix_a%m%special_memory, error=error)
       ENDIF
       CALL ensure_array_size(matrix_b%m%index, ub=SIZE(matrix_a%m%index),&
            error=error)
       !
       ! copy index and data
       matrix_b%m%index(1:SIZE(matrix_a%m%index)) = matrix_a%m%index(:)
       IF (.NOT. shallow) THEN
          CALL cp_assert (matrix_a%m%nze .LE. dbcsr_get_data_size(matrix_a),&
               cp_fatal_level, cp_internal_error, routineN,&
               "Source matrix sizes not consistent!", error=error)
          CALL dbcsr_data_ensure_size (matrix_b%m%data_area,&
               dbcsr_get_data_size(matrix_a))
          IF (dbcsr_get_data_type(matrix_a) .EQ. dbcsr_get_data_type(matrix_b))&
               THEN
             CALL dbcsr_data_copyall (matrix_b%m%data_area,&
                  matrix_a%m%data_area)
          ELSE
             CALL dbcsr_data_convert (matrix_b%m%data_area,&
                  matrix_a%m%data_area, drop_real=keep_imaginary)
          ENDIF
       ELSE
          CALL cp_assert (dbcsr_get_data_type(matrix_a) &
               .EQ. dbcsr_get_data_type(matrix_b), cp_fatal_level,&
               cp_wrong_args_error, routineN,&
               "Shallow copy only possible when retaining data type.",&
               error=error)
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

       matrix_b%m%sparsity_id = matrix_a%m%sparsity_id
    ENDIF
    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_copy

! *****************************************************************************
!> \brief copy a matrix, retaining current sparsity
!> \param[inout] matrix_b       target DBCSR matrix
!> \param[in]    matrix_a       source DBCSR matrix
!> \param[in,out] error         cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_copy_into_existing(matrix_b, matrix_a, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_into_existing', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER                                  :: col_size, cs, dst_col, &
                                                dst_row, rel, row_size, rs, &
                                                src_col, src_row, &
                                                timing_handle
    LOGICAL                                  :: dst_tr, src_tr
    TYPE(dbcsr_data_obj)                     :: dst_data, src_data
    TYPE(dbcsr_iterator)                     :: dst_iter, src_iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL cp_assert (dbcsr_get_data_type(matrix_b)&
         .EQ. dbcsr_get_data_type(matrix_a), cp_fatal_level,&
         cp_wrong_args_error, routineN, "Matrices have different data types.",&
         error=error)
    data_type = dbcsr_get_data_type (matrix_b)

    CALL dbcsr_data_init (src_data)
    CALL dbcsr_data_init (dst_data)
    CALL dbcsr_data_new (src_data, data_type)
    CALL dbcsr_data_new (dst_data, data_type)
    CALL dbcsr_iterator_start (src_iter, matrix_a)
    CALL dbcsr_iterator_start (dst_iter, matrix_b)
    ! Iterate through the blocks of the source and destination
    ! matrix. There are three possibilites: 1. copy the data for
    ! blocks present in both; 2 skip source blocks not present in the
    ! target; 3 zero blocks not present in the source.
    IF (dbcsr_iterator_blocks_left (src_iter)) THEN
       CALL dbcsr_iterator_next_block (src_iter, src_row, src_col, src_data,&
            src_tr)
    ELSE
       src_row = 0 ; src_col = 0
    ENDIF
    DO WHILE (dbcsr_iterator_blocks_left (dst_iter))
       CALL dbcsr_iterator_next_block (dst_iter, dst_row, dst_col, dst_data,&
            dst_tr, row_size=row_size, col_size=col_size)
       ! Now find the source position that is greater or equal to the
       ! target one. I.e, skip blocks that the target doesn't have.
       rel = pos_relation (dst_row, dst_col, src_row, src_col)
       DO WHILE (rel .EQ. 1 .AND. dbcsr_iterator_blocks_left (src_iter))
          CALL dbcsr_iterator_next_block (src_iter, src_row, src_col,&
               src_data, src_tr)
          rel = pos_relation (dst_row, dst_col, src_row, src_col)
       ENDDO
       SELECT CASE (rel)
       CASE (-1, 1)
           ! Target lags source or ran out of source
          CALL dbcsr_data_clear (dst_data)
          !dst_r_dp(:) = 0.0_dp
       CASE (0)
          ! Copy the data
          CALL cp_assert (dbcsr_get_data_size (src_data)&
               .EQ. dbcsr_get_data_size (dst_data), cp_fatal_level,&
               cp_internal_error, routineN, "Block sizes not equal!",&
               error=error)
          IF (src_tr .EQV. dst_tr) THEN
             CALL dbcsr_data_copyall (dst_data, src_data)
             !dst_r_dp(:) = src_r_dp(:)
          ELSE
             rs = row_size
             cs = col_size
             IF (src_tr) CALL swap (rs, cs)
             CALL dbcsr_block_transpose (dst_data, src_data, rs, cs)
             !CALL dbcsr_block_transpose (dst_r_dp, src_r_dp, rs, cs)
          ENDIF
       CASE default
          CALL cp_assert (.FALSE., cp_fatal_level, cp_internal_error,&
               routineN, "Trouble syncing iterators", error=error)
       END SELECT
    END DO
    CALL dbcsr_iterator_stop (src_iter)
    CALL dbcsr_iterator_stop (dst_iter)
    CALL dbcsr_data_clear_pointer (src_data)
    CALL dbcsr_data_clear_pointer (dst_data)
    CALL dbcsr_data_release (src_data)
    CALL dbcsr_data_release (dst_data)

    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_copy_into_existing


! *****************************************************************************
!> \brief Copy only a subset of matrix columns (element-wise, not block-size)
!> \par Column block compatibility
!>      The process column distributions of the copied columns must be
!>      equal. Therefore, the source and target offsets must also be
!>      the same.
!> \param[in,out] matrix_b      target DBCSR matrix
!> \param[in]    matrix_a       source DBCSR matrix
!> \param[in]    ncol           number of elemental column to copy
!> \param[in]    source_start   starting elemental column in source
!> \param[in]    target_start   starting elemental column in target
!> \param[in,out] error         cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_copy_columns(matrix_b, matrix_a,&
       ncol, source_start, target_start, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    INTEGER, INTENT(IN)                      :: ncol, source_start, &
                                                target_start
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_columns', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER :: col_size, dst_col, dst_row, first_col, last_col, num_col, rel, &
      row_size, src_col, src_row, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: dst_col_offsets, &
                                                dst_col_sizes, &
                                                src_col_offsets, src_col_sizes
    LOGICAL                                  :: dst_tr, src_tr
    TYPE(dbcsr_data_obj)                     :: dst_data, src_data, tmp_buffer
    TYPE(dbcsr_iterator)                     :: dst_iter, src_iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL cp_assert (dbcsr_get_data_type(matrix_b)&
         .EQ. dbcsr_get_data_type(matrix_a), cp_fatal_level,&
         cp_wrong_args_error, routineN, "Matrices have different data types.",&
         error=error)
    data_type = dbcsr_get_data_type (matrix_b)

    CALL cp_assert (source_start .EQ. target_start, cp_fatal_level,&
         cp_unimplemented_error_nr, routineN, "Column shifting not supported",&
         error=error)

    src_col_offsets => array_data (matrix_a%m%col_blk_offset)
    dst_col_offsets => array_data (matrix_b%m%col_blk_offset)
    src_col_sizes => array_data (matrix_a%m%col_blk_size)
    dst_col_sizes => array_data (matrix_b%m%col_blk_size)

    CALL dbcsr_data_init (src_data)
    CALL dbcsr_data_init (dst_data)
    CALL dbcsr_data_init (tmp_buffer)
    CALL dbcsr_data_new (src_data, data_type)
    CALL dbcsr_data_new (dst_data, data_type)
    CALL dbcsr_data_new (tmp_buffer, data_type,&
         matrix_a%m%max_rbs*matrix_a%m%max_cbs)
    CALL dbcsr_iterator_start (src_iter, matrix_a)
    CALL dbcsr_iterator_start (dst_iter, matrix_b)
    ! Iterate through the blocks of the source and destination
    ! matrix. There are three possibilites: 1 copy data for blocks
    ! present in both that are fully within the specified range; 2
    ! copy partial data for blacks partially within the specified
    ! range; 3 add (partial or full) data for blocks present in source
    ! but not in the target
    IF (dbcsr_iterator_blocks_left (src_iter)) THEN
       CALL dbcsr_iterator_next_block (src_iter, src_row, src_col, src_data,&
            src_tr)
    ELSE
       src_row = 0 ; src_col = 0
    ENDIF
    DO WHILE (dbcsr_iterator_blocks_left (dst_iter))
       CALL dbcsr_iterator_next_block (dst_iter, dst_row, dst_col, dst_data,&
            dst_tr, row_size=row_size, col_size=col_size)
       ! Now find the source position that is greater or equal to the
       ! target one. I.e, skip blocks that the target doesn't have.
       rel = pos_relation (dst_row, dst_col, src_row, src_col)
       DO WHILE (rel .EQ. 1 .AND. dbcsr_iterator_blocks_left (src_iter))
          CALL dbcsr_iterator_next_block (src_iter, src_row, src_col,&
               src_data, src_tr)
          rel = pos_relation (dst_row, dst_col, src_row, src_col)
       ENDDO
       SELECT CASE (rel)
       CASE (0) ! Equal block coordinates; block exists in both src and dst
          ! Check if the full block is being copied
          CALL cp_assert (src_col_sizes(src_col) .EQ. dst_col_sizes(dst_col),&
               cp_fatal_level, cp_internal_error, routineN,&
               "Unequal column sizes", error=error)
          !CALL cp_unimplemented_error (routineN, "Partial block copy not tested",&
          !     error=error, error_level=cp_warning_level)
          CALL cp_assert (src_col_offsets(src_col).eq.dst_col_offsets(dst_col),&
               cp_fatal_level, cp_internal_error, routineN,&
               "Block offsets must be equal", error=error)
          first_col = 1 + MAX (src_col_offsets(src_col), source_start)&
               - src_col_offsets(src_col)
          last_col = MIN (src_col_offsets(src_col)+ncol-1, source_start+ncol-1) - &
               src_col_offsets(src_col+1)-1 + 1
          num_col = last_col - first_col + 1
          IF (num_col .GT. 0) THEN
             CALL dbcsr_block_partial_copy(&
                  dst_data, row_size, col_size, dst_tr,&
                  src_data, row_size, col_size, src_tr,&
                  1, first_col, 1, first_col,& ! offsets
                  row_size, num_col) ! sizes
          ENDIF
       CASE (-1) ! Block exists in dst but must be added to dst
          IF (source_start .GE. src_col_offsets(src_col) .AND. &
              src_col_offsets(src_col)+ncol-1 .LE. source_start+ncol-1) THEN
             CALL dbcsr_put_block (matrix_b, src_row, src_col, src_data,&
                  transposed=src_tr)
          ELSE
             CALL cp_unimplemented_error (routineN,&
                  "Partial block addition not yet tested",&
                  error=error, error_level=cp_warning_level)
             CALL dbcsr_data_clear (tmp_buffer)
             CALL dbcsr_block_partial_copy (&
                  tmp_buffer, row_size, col_size, src_tr,&
                  src_data, row_size, col_size, src_tr,&
                  1, first_col, 1, first_col,&
                  row_size, num_col)
             CALL dbcsr_put_block (matrix_b, src_row, src_col, tmp_buffer,&
                  src_tr)
          ENDIF
       CASE (1)
          ! Ran out of source. Do nothing.
       CASE default
          CALL cp_assert (.FALSE., cp_fatal_level, cp_internal_error,&
               routineN, "Trouble syncing iterators", error=error)
       END SELECT
    END DO
    CALL dbcsr_iterator_stop (src_iter)
    CALL dbcsr_iterator_stop (dst_iter)
    CALL dbcsr_data_clear_pointer (src_data)
    CALL dbcsr_data_clear_pointer (dst_data)
    CALL dbcsr_data_release (src_data)
    CALL dbcsr_data_release (dst_data)

    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_copy_columns


! *****************************************************************************
!> \brief Determines the relation between two matrix positions.
!> \retval relation  Relation between positions 1 and 2.
!>                    0: same
!>                   -1: pos1 < pos2
!>                    1: pos1 > pos2
! *****************************************************************************
  ELEMENTAL FUNCTION pos_relation (row1, col1, row2, col2) RESULT (relation)
    INTEGER, INTENT(IN)                      :: row1, col1, row2, col2
    INTEGER                                  :: relation

    IF (row1 .LT. row2) THEN
       relation = -1
    ELSEIF (row1 .GT. row2) THEN
       relation = 1
    ELSE ! rows are equal, check column
       IF (col1 .LT. col2) THEN
          relation = -1
       ELSEIF (col1 .GT. col2) THEN
          relation = 1
       ELSE
          relation = 0
       ENDIF
    ENDIF
  END FUNCTION pos_relation


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

    INTEGER                                  :: blk, blk_nze, col, col_size, &
                                                row, row_size, timing_handle
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL dbcsr_copy(matrix_b, matrix_a, name="triu of "//matrix_a%m%name, error=error)

    CALL dbcsr_iterator_start(iter, matrix_b)

    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, blk=blk,&
            row_size=row_size, col_size=col_size)
       blk_nze = row_size * col_size
       IF(row.GT.col) CALL dbcsr_remove_block(matrix_b, row, col, blk_nze, blk)
    ENDDO

    CALL dbcsr_iterator_stop(iter)

    CALL dbcsr_finalize (matrix_b, error=error)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_btriu

! *****************************************************************************
!> \brief  triu of a dbcsr matrix
!> \param[inout] matrix_a  the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_triu(matrix_a, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_triu', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_nze, col, col_size, &
                                                i, j, row, row_size, &
                                                timing_handle
    LOGICAL                                  :: tr
    REAL(dp), DIMENSION(:, :), POINTER       :: DATA
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL dbcsr_iterator_start(iter, matrix_a)

    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, DATA, tr, &
            block_number=blk, row_size=row_size, col_size=col_size)
       blk_nze = row_size * col_size
       IF(row.GT.col) CALL dbcsr_remove_block(matrix_a, row, col, blk_nze, blk)
       IF(row.EQ.col) THEN
          DO j = 1,col_size
          DO i = j+1,row_size
             DATA(i,j) = 0.0_dp
          ENDDO
          ENDDO
       ENDIF
    ENDDO

    CALL dbcsr_iterator_stop(iter)

    CALL dbcsr_finalize (matrix_a, error=error)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_triu

! *****************************************************************************
!> \brief  symmetrize the block diagonal of a dbcsr matrix
!> \param[inout] matrix_a  the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_symmetrize_block_diag(matrix_a, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_symmetrize_block_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_size, i, j, &
                                                row, row_size, timing_handle
    LOGICAL                                  :: tr
    REAL(dp)                                 :: dum
    REAL(dp), DIMENSION(:, :), POINTER       :: DATA
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL dbcsr_iterator_start(iter, matrix_a)

    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, DATA, tr, &
            block_number=blk, row_size=row_size, col_size=col_size)
       IF(col.NE.row) CYCLE
       DO j = 1,col_size
       DO i = j+1,row_size
          dum = (DATA(i,j) + DATA(j,i)) / 2.0_dp
          DATA(i,j) = dum
          DATA(j,i) = dum
       ENDDO
       ENDDO
    ENDDO

    CALL dbcsr_iterator_stop(iter)

    CALL dbcsr_finalize(matrix_a, error=error)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_symmetrize_block_diag

! *****************************************************************************
!> \brief  tril of a dbcsr matrix
!> \param[inout] matrix_a  the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_tril(matrix_a, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_tril', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_nze, col, col_size, &
                                                i, j, row, row_size, &
                                                timing_handle
    LOGICAL                                  :: tr
    REAL(dp), DIMENSION(:, :), POINTER       :: DATA
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    CALL dbcsr_iterator_start(iter, matrix_a)

    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, block=DATA, transposed=tr, &
            block_number=blk, row_size=row_size, col_size=col_size)
       blk_nze = row_size * col_size
       IF(row.GT.col) CALL dbcsr_remove_block(matrix_a, row, col, blk_nze, blk)
       IF(row.EQ.col) THEN
          DO j = 1,col_size
          DO i = 1,j-1
             DATA(i,j) = 0.0_dp
          ENDDO
          ENDDO
       ENDIF
    ENDDO

    CALL dbcsr_iterator_stop(iter)

    CALL dbcsr_finalize (matrix_a, error=error)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_tril

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
    INTEGER, INTENT(IN), OPTIONAL            :: method
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_filter', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_nze, col, col_size, &
                                                my_method, row, row_size, &
                                                timing_handle
    LOGICAL                                  :: tr
    REAL(dp)                                 :: ddot, nrm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    my_method = dbcsr_filter_frobenius
    IF(PRESENT(method)) my_method = method

    IF(eps.GT.0.0_dp) THEN

       CALL dbcsr_iterator_start(iter, matrix, contiguous_pointers=.TRUE.)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          SELECT CASE (dbcsr_get_data_type(matrix))
             !CASE (dbcsr_type_real_4)
          CASE (dbcsr_type_real_8)
             CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk, &
                  row_size, col_size)
             blk_nze = row_size * col_size
             IF (blk_nze .EQ. 0) CYCLE ! Skip empty blocks
             SELECT CASE(my_method)
             CASE(dbcsr_filter_frobenius)
                !
                ! Frobenius based
                nrm = SQRT(DDOT(SIZE(data_d), data_d(1), 1, data_d(1), 1))
                IF(nrm.LT.eps) CALL dbcsr_remove_block(matrix, row, col, blk_nze, blk)
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

    ENDIF

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
    INTEGER, INTENT(IN)                      :: which_norm
    REAL(dp), INTENT(OUT), OPTIONAL          :: norm_scalar
    REAL(dp), DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: norm_vector
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_norm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_offset, i, j, &
                                                row, row_offset
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    SELECT CASE(which_norm)

    CASE(dbcsr_norm_frobenius)

       IF(PRESENT(norm_scalar)) norm_scalar = dbcsr_frobenius_norm(matrix)

    CASE(dbcsr_norm_maxabsnorm)

       IF(PRESENT(norm_scalar)) norm_scalar = dbcsr_maxabs_norm(matrix)

    CASE(dbcsr_norm_gershgorin)

       IF(PRESENT(norm_scalar)) norm_scalar = dbcsr_gershgorin_norm(matrix)

    CASE(dbcsr_norm_column)

       IF(PRESENT(norm_vector)) THEN
          norm_vector(:) = 0.0_dp
          CALL dbcsr_iterator_start(iter, matrix)
          DO WHILE (dbcsr_iterator_blocks_left(iter))
             SELECT CASE (dbcsr_get_data_type(matrix))
             !CASE (dbcsr_type_real_4)
             CASE (dbcsr_type_real_8)
                CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk, &
                     row_offset=row_offset, col_offset=col_offset)
                IF((dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_symmetric.OR.&
                    dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_antisymmetric).AND.&
                    row.NE.col) THEN
                   DO j=1,SIZE(data_d,2)
                   DO i=1,SIZE(data_d,1)
                      norm_vector(col_offset+j-1) = norm_vector(col_offset+j-1) &
                                                     + data_d(i,j)**2
                      norm_vector(row_offset+i-1) = norm_vector(row_offset+i-1) &
                                                     + data_d(i,j)**2
                   ENDDO
                   ENDDO
                ELSE
                   DO j=1,SIZE(data_d,2)
                   DO i=1,SIZE(data_d,1)
                      norm_vector(col_offset+j-1) = norm_vector(col_offset+j-1) &
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

       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error, routineN, "this norm is NYI")

    END SELECT

  END SUBROUTINE dbcsr_norm

! *****************************************************************************
!> \brief compute a norm of a dbcsr matrix
!> \param[in] matrix  the matrix
!> \param[out] norm
!> 
! *****************************************************************************
  FUNCTION dbcsr_gershgorin_norm(matrix)

    TYPE(dbcsr_obj), INTENT(INOUT), TARGET   :: matrix
    REAL(dp)                                 :: dbcsr_gershgorin_norm

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_gershgorin_norm', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: data_z
    COMPLEX(KIND=sp), DIMENSION(:, :), &
      POINTER                                :: data_c
    INTEGER                                  :: blk, col, col_offset, i, j, &
                                                nc, nr, row, row_offset, &
                                                timing_handle
    LOGICAL                                  :: any_sym, tr
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff_d
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: data_r
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    nr = dbcsr_nfullrows_total(matrix)
    nc = dbcsr_nfullcols_total(matrix)

    any_sym = dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_symmetric.OR.&
              dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_antisymmetric

    CALL cp_assert (nr.EQ.nc, cp_fatal_level, cp_wrong_args_error, &
         routineN, "not a square matrix")

    dbcsr_gershgorin_norm = 0.0_dp
    ALLOCATE(buff_d(nr))
    buff_d = 0.0_dp
    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       SELECT CASE (dbcsr_get_data_type(matrix))
       CASE (dbcsr_type_real_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_r, tr, blk, &
               row_offset=row_offset, col_offset=col_offset)
          DO j=1,SIZE(data_r,2)
          DO i=1,SIZE(data_r,1)
             buff_d(row_offset+i-1) = buff_d(row_offset+i-1) + ABS(data_r(i,j))
             IF(any_sym.AND.row.NE.col) &
                  buff_d(col_offset+j-1) = buff_d(col_offset+j-1) + ABS(data_r(i,j))
          ENDDO
          ENDDO
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk, &
               row_offset=row_offset, col_offset=col_offset)
          DO j=1,SIZE(data_d,2)
          DO i=1,SIZE(data_d,1)
             buff_d(row_offset+i-1) = buff_d(row_offset+i-1) + ABS(data_d(i,j))
             IF(any_sym.AND.row.NE.col) &
                  buff_d(col_offset+j-1) = buff_d(col_offset+j-1) + ABS(data_d(i,j))
          ENDDO
          ENDDO
       CASE (dbcsr_type_complex_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_c, tr, blk, &
               row_offset=row_offset, col_offset=col_offset)
          DO j=1,SIZE(data_c,2)
          DO i=1,SIZE(data_c,1)
             buff_d(row_offset+i-1) = buff_d(row_offset+i-1) + ABS(data_c(i,j))
             IF(any_sym.AND.row.NE.col) &
                  CALL cp_unimplemented_error(routineN,"Only nonsymmetric matrix so far")
             !     buff_d(col_offset+j-1) = buff_d(col_offset+j-1) + ABS(data_c(i,j))
          ENDDO
          ENDDO
       CASE (dbcsr_type_complex_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_z, tr, blk, &
               row_offset=row_offset, col_offset=col_offset)
          DO j=1,SIZE(data_z,2)
          DO i=1,SIZE(data_z,1)
             buff_d(row_offset+i-1) = buff_d(row_offset+i-1) + ABS(data_z(i,j))
             IF(any_sym.AND.row.NE.col) &
                  CALL cp_unimplemented_error(routineN,"Only nonsymmetric matrix so far")
             !     buff_d(col_offset+j-1) = buff_d(col_offset+j-1) + ABS(data_z(i,j))
          ENDDO
          ENDDO
       CASE DEFAULT
          CALL cp_unimplemented_error (routineN,"Wrong data type")
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL mp_sum(buff_d,dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))       
    dbcsr_gershgorin_norm = MAXVAL(buff_d)
    DEALLOCATE(buff_d)

    CALL timestop(timing_handle)

  END FUNCTION dbcsr_gershgorin_norm

! *****************************************************************************
!> \brief compute a norm of a dbcsr matrix
!> \param[in] matrix  the matrix
!> \param[out] norm
!> 
! *****************************************************************************
  FUNCTION dbcsr_maxabs_norm(matrix)
    TYPE(dbcsr_obj), INTENT(INOUT), TARGET   :: matrix
    REAL(dp)                                 :: dbcsr_maxabs_norm

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_maxabs_norm', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: data_z
    COMPLEX(KIND=sp), DIMENSION(:, :), &
      POINTER                                :: data_c
    INTEGER                                  :: blk, col, row, timing_handle
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: data_r
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    dbcsr_maxabs_norm = 0.0_dp
    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       SELECT CASE (dbcsr_get_data_type(matrix))
       CASE (dbcsr_type_real_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_r, tr, blk)
          dbcsr_maxabs_norm = MAX(dbcsr_maxabs_norm,REAL(MAXVAL(ABS(data_r)),dp))
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
          dbcsr_maxabs_norm = MAX(dbcsr_maxabs_norm,MAXVAL(ABS(data_d)))
       CASE (dbcsr_type_complex_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_c, tr, blk)
          dbcsr_maxabs_norm = MAX(dbcsr_maxabs_norm,REAL(MAXVAL(ABS(data_c)),dp))
       CASE (dbcsr_type_complex_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_z, tr, blk)
          dbcsr_maxabs_norm = MAX(dbcsr_maxabs_norm,MAXVAL(ABS(data_z)))
       CASE DEFAULT
          CALL cp_unimplemented_error (routineN,"Wrong data type")
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)

    !dmp_max, this fixes a bug in g95
    ! -> ambigous interface for mp_max
    !    in module dbcsr_message_passing and wrongly
    !    exported from the module timings
    CALL dmp_max( dbcsr_maxabs_norm, dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))

    CALL timestop(timing_handle)

  END FUNCTION dbcsr_maxabs_norm

! *****************************************************************************
!> \brief compute a norm of a dbcsr matrix
!> \param[in] matrix  the matrix
!> \param[out] norm
!> 
! *****************************************************************************
  FUNCTION dbcsr_frobenius_norm(matrix)

    TYPE(dbcsr_obj), INTENT(INOUT), TARGET   :: matrix
    REAL(dp)                                 :: dbcsr_frobenius_norm

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_frobenius_norm', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: data_z
    COMPLEX(KIND=sp), DIMENSION(:, :), &
      POINTER                                :: data_c
    INTEGER                                  :: blk, col, row, timing_handle
    LOGICAL                                  :: any_sym, tr
    REAL(dp)                                 :: fac
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: data_r
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    any_sym = dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_symmetric.OR.&
              dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_antisymmetric

    dbcsr_frobenius_norm = 0.0_dp
    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       SELECT CASE (dbcsr_get_data_type(matrix))
       CASE (dbcsr_type_real_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_r, tr, blk)
          fac = 1.0_dp
          IF(any_sym.AND.row.NE.col) fac = 2.0_dp
          dbcsr_frobenius_norm = dbcsr_frobenius_norm + fac * SUM(data_r**2)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
          fac = 1.0_dp
          IF(any_sym.AND.row.NE.col) fac = 2.0_dp
          dbcsr_frobenius_norm = dbcsr_frobenius_norm + fac * SUM(data_d**2)
       CASE (dbcsr_type_complex_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_c, tr, blk)
          fac = 1.0_dp
          IF(any_sym.AND.row.NE.col) &
               CALL cp_unimplemented_error(routineN,"Only nonsymmetric matrix so far")
          dbcsr_frobenius_norm = dbcsr_frobenius_norm + fac * SUM(CONJG(data_c)*data_c)
       CASE (dbcsr_type_complex_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_z, tr, blk)
          fac = 1.0_dp
          IF(any_sym.AND.row.NE.col) &
               CALL cp_unimplemented_error(routineN,"Only nonsymmetric matrix so far")
          dbcsr_frobenius_norm = dbcsr_frobenius_norm + fac * SUM(CONJG(data_z)*data_z)
       CASE DEFAULT
          CALL cp_unimplemented_error(routineN,"Wrong data type")
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL mp_sum(dbcsr_frobenius_norm,dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))
    dbcsr_frobenius_norm = SQRT(dbcsr_frobenius_norm)

    CALL timestop(timing_handle)

  END FUNCTION dbcsr_frobenius_norm

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
         dbcsr_type_normal, matrix%m%row_blk_size, matrix%m%col_blk_size,&
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
            SUM(recv_count(2,:)), n=1)
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
!> \param[in] keep_sparsity   (optional) retains the sparsity of the redist
!>                            matrix
! *****************************************************************************
  SUBROUTINE dbcsr_complete_redistribute(matrix, redist, keep_sparsity)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_obj), INTENT(INOUT)           :: redist
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_complete_redistribute', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 7

    CHARACTER                                :: data_type
    INTEGER :: blk, blk_col_new, blk_col_prev_new, blk_ps, blk_row_new, blks, &
      cnt_fnd, cnt_new, cnt_skip, col, col_int, col_offset_new, &
      col_offset_old, col_rle, col_size, data_offset_l, dst_p, i, meta_l, &
      mp_group, nblkcols_total_new, nblkrows_total_new, numnodes, nze_rle, &
      row, row_int, row_offset_new, row_offset_old, row_prev, row_rle, &
      row_size, src_p, stored_col_new, stored_row_new, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:) :: col_end_new, col_end_old, &
      col_start_new, col_start_old, rd_disp, recv_meta, rm_disp, row_end_new, &
      row_end_old, row_start_new, row_start_old, sd_disp, sdp, send_meta, &
      sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: col_reblocks, n_col_reblocks, &
      n_row_reblocks, recv_count, row_reblocks, send_count, total_recv_count, &
      total_send_count
    INTEGER, DIMENSION(:), POINTER :: col_blk_size_new, col_blk_size_old, &
      col_dist_new, row_blk_size_new, row_blk_size_old, row_dist_new
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    LOGICAL                                  :: found, my_keep_sparsity, &
                                                new_meta, sym, tr, valid_block
    REAL(kind=dp)                            :: t_all, t_blk, t_blk2, &
                                                t_count, t_fill, t_prolog, &
                                                t_unpack, t_xfer1, t_xfer2
    TYPE(dbcsr_data_obj)                     :: buff_data, data_block, &
                                                recv_data, send_data, src_2d
    TYPE(dbcsr_distribution_obj)             :: dist_new, dist_old
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj_new

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    t_all = m_walltime()
    t_prolog = t_all
    CALL cp_assert (dbcsr_valid_index(matrix%m),&
         cp_fatal_level, cp_caller_error,&
         routineN, "Input not valid.")
    CALL cp_assert (dbcsr_is_initialized(redist%m),&
         cp_fatal_level, cp_caller_error,&
         routineN, "Output matrix not initialized.")
    CALL cp_assert (matrix%m%replication_type .EQ. dbcsr_repl_none,&
         cp_warning_level, cp_caller_error, routineN,&
         "Can not redistribute replicated matrix.")
    !
    my_keep_sparsity = .FALSE.
    IF (PRESENT (keep_sparsity)) my_keep_sparsity = my_keep_sparsity
    !
    sym = dbcsr_get_matrix_type(redist) .EQ. dbcsr_type_symmetric .OR.&
         dbcsr_get_matrix_type(redist) .EQ. dbcsr_type_antisymmetric
    data_type = matrix%m%data_type
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
    ! Create mappings
    i = dbcsr_nfullrows_total (redist)
    ALLOCATE (row_reblocks (2,i))
    ALLOCATE (n_row_reblocks (2, dbcsr_nblkrows_total (matrix)))
    CALL dbcsr_reblocking_targets (row_reblocks, i, n_row_reblocks,&
         row_blk_size_old, row_blk_size_new)
    i = dbcsr_nfullcols_total (redist)
    ALLOCATE (col_reblocks (2,i))
    ALLOCATE (n_col_reblocks (2, dbcsr_nblkcols_total (matrix)))
    CALL dbcsr_reblocking_targets (col_reblocks, i, n_col_reblocks,&
         col_blk_size_old, col_blk_size_new)
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
    IF (dbg) CALL dbcsr_print(matrix)
    !cs1 = dbcsr_checksum (matrix)
    !cs1 = dbcsr_checksum (matrix)
    !call dbcsr_print(matrix)
    !
    t_count = m_walltime()
    t_prolog = t_count - t_prolog
    !
    ! Count initial sizes for sending.
    !
    ! We go through every element of every local block and determine
    ! to which processor it must be sent. It could be more efficient,
    ! but at least the index data are run-length encoded.
    DBG 'A'
    send_count(:,:) = 0
    CALL dbcsr_iterator_start (iter, matrix)
    dst_p = -1
    new_meta = .TRUE.
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk)
       DO col_int = n_col_reblocks (1,col),&
            n_col_reblocks(1,col) + n_col_reblocks(2,col)-1
          blk_col_new = col_reblocks(1, col_int)
          DO row_int = n_row_reblocks (1, row),&
               n_row_reblocks(1,row) + n_row_reblocks(2,row)-1
             blk_row_new = row_reblocks(1, row_int)
             IF (.NOT. sym .OR. blk_col_new .GE. blk_row_new) THEN
                tr = .FALSE.
                CALL dbcsr_get_stored_coordinates (redist%m,&
                     blk_row_new, blk_col_new, tr, dst_p)
                send_count(1, dst_p) = send_count(1, dst_p) + 1
                send_count(2, dst_p) = send_count(2, dst_p) +&
                     col_reblocks(2,col_int)*row_reblocks(2,row_int)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    !
    t_xfer1 = m_walltime()
    t_count = t_xfer1 - t_count
    !
    DBG 'B'
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    DBG 'C'
    DBG 'send counts',send_count
    DBG 'recv counts',recv_count
    ! Allocate data structures needed for data exchange.
    CALL dbcsr_data_init (recv_data)
    CALL dbcsr_data_new (recv_data, data_type, SUM(recv_count(2, :)))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :))))
    CALL dbcsr_data_init (send_data)
    CALL dbcsr_data_new (send_data, data_type, SUM(send_count(2, :)))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :))))
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
    !
    t_fill = m_walltime()
    t_xfer1 = t_fill - t_xfer1
    !
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_new (data_block, data_type)
    CALL dbcsr_iterator_start (iter, matrix)
    row_prev = 0
    dst_p = -1
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, data_block, tr, blk,&
            row_size=row_size, col_size=col_size)
       IF (row .NE. row_prev) blk_col_prev_new = 0 ! Reset the "hint"
       !IF (tr) WRITE(*,*)"block at",row,col," is transposed"
       col_offset_old = 1
       DO col_int = n_col_reblocks (1,col),&
            n_col_reblocks(1,col) + n_col_reblocks(2,col)-1
          blk_col_new = col_reblocks(1, col_int)
          row_offset_old = 1
          DO row_int = n_row_reblocks (1, row),&
               n_row_reblocks(1,row) + n_row_reblocks(2,row)-1
             blk_row_new = row_reblocks(1, row_int)
             loc_ok: IF (.NOT. sym .OR. blk_col_new .GE. blk_row_new) THEN
                IF (dbg) &
                     WRITE(*,*)'using block',blk_row_new,'x',blk_col_new
                ! Start a new RLE run
                tr = .FALSE.
                CALL dbcsr_get_stored_coordinates (redist%m,&
                     blk_row_new, blk_col_new, tr, dst_p)
                row_rle = row_reblocks(2,row_int)
                col_rle = col_reblocks(2,col_int)
                row_offset_new = row_offset_old + row_start_old(row) - &
                     row_start_new(blk_row_new)
                col_offset_new = col_offset_old + col_start_old(col) - &
                     col_start_new(blk_col_new)
                smp(dst_p) = smp(dst_p) + metalen
                send_meta(smp(dst_p)) = blk_row_new   ! new blocked row
                send_meta(smp(dst_p)+1) = blk_col_new ! new blocked column
                send_meta(smp(dst_p)+2) = row_offset_new  ! row in new block
                send_meta(smp(dst_p)+3) = col_offset_new  ! col in new block
                send_meta(smp(dst_p)+4) = row_rle ! RLE rows
                send_meta(smp(dst_p)+5) = col_rle ! RLE columns
                send_meta(smp(dst_p)+6) = sdp(dst_p)-sd_disp(dst_p) ! Offset in data
                nze_rle = row_rle * col_rle
                ! Copy current block into the send buffer
                CALL dbcsr_block_partial_copy(&
                     send_data, dst_offset=sdp(dst_p)-1,&
                     dst_rs=row_rle, dst_cs=col_rle, dst_tr=.FALSE.,&
                     dst_r_lb=1, dst_c_lb=1,&
                     src=data_block,&
                     src_rs=row_size, src_cs=col_size, src_tr=tr,&
                     src_r_lb=row_offset_old, src_c_lb=col_offset_old,&
                     nrow=row_rle, ncol=col_rle)
                !IF (.NOT. tr) THEN
                !      CALL dbcsr_block_copy(&
                !           send_data(sdp(dst_p):sdp(dst_p)+nze_rle-1),&
                !           r2_dp(&
                !           row_offset_old:row_offset_old+row_rle-1,&
                !           col_offset_old:col_offset_old+col_rle-1),&
                !           row_rle, col_rle)
                !ELSE
                !   CALL dbcsr_block_transpose(&
                !        send_data(sdp(dst_p):sdp(dst_p)+nze_rle-1),&
                !        r2_dp(&
                !        col_offset_old:col_offset_old+col_rle-1,&
                !        row_offset_old:row_offset_old+row_rle-1),&
                !        row_rle, col_rle)
                !ENDIF
                sdp(dst_p) = sdp(dst_p)+nze_rle
             ENDIF loc_ok
             row_offset_old = row_offset_old + row_reblocks(2, row_int)
          ENDDO ! row_int
          col_offset_old = col_offset_old + col_reblocks(2, col_int)
       ENDDO ! col_int
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    CALL dbcsr_data_clear_pointer (data_block)
    CALL dbcsr_data_release (data_block)

    ! Exchange the data and metadata structures.
    !
    t_xfer2 = m_walltime()
    t_fill = t_xfer2 - t_fill
    !
    DBG 'E'
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       CALL hybrid_alltoall_s1(&
            send_data%d%r_sp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%r_sp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_real_8)
       !CALL mp_alltoall(&
       !     send_data%d%r_dp(:), total_send_count(2,:), sd_disp(:)-1,&
       !     recv_data%d%r_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
       !     mp_group)
       CALL hybrid_alltoall_d1 (&
            send_data%d%r_dp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%r_dp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_complex_4)
       CALL hybrid_alltoall_c1(&
            send_data%d%c_sp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%c_sp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_complex_8)
       CALL hybrid_alltoall_z1(&
            send_data%d%c_dp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%c_dp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE default
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error,&
            routineN, "Invalid matrix type")
    END SELECT
    DBG 'F'
    DBG 'send_data'
    CALL hybrid_alltoall_i1(send_meta(:), SUM(metalen*total_send_count(1,:)), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), SUM(metalen*total_recv_count(1,:)), metalen*total_recv_count(1,:), rm_disp(:)-1, mp_obj_new)
    DBG 'G'
    !
    t_unpack = m_walltime()
    t_xfer2 = t_unpack - t_xfer2
    !
    ! Now fill in the data.
    CALL dbcsr_work_create(redist,&
         nblks_guess=SUM(recv_count(1,:)),&
         sizedata_guess=SUM(recv_count(2,:)), work_mutable=.TRUE.)
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    CALL dbcsr_data_init (buff_data)
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_init (src_2d)
    CALL dbcsr_data_new (buff_data, dbcsr_type_1d_to_2d (data_type), &
         redist%m%max_rbs, redist%m%max_cbs)
    CALL dbcsr_data_new (data_block, dbcsr_type_1d_to_2d (data_type))
    CALL dbcsr_data_new (src_2d, dbcsr_type_1d_to_2d (data_type))

    !blk_p = 1
    !blk = 1
    blk_ps = 0
    blks = 0
    t_blk2 = 0.0_dp
    cnt_fnd = 0 ; cnt_new = 0 ; cnt_skip = 0
    DO src_p = 0, numnodes-1
       data_offset_l = rd_disp(src_p)
       DO meta_l = 1, recv_count(1, src_p)
          stored_row_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1))
          stored_col_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+1)
          row_offset_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+2)
          col_offset_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+3)
          row_rle = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+4)
          col_rle = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+5)
          data_offset_l = rd_disp(src_p)&
               + recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+6)
          DBG 'blk: p',blks,src_p,'; row, col',&
               stored_row_new,stored_col_new,&
               '; blk_ps',blk_ps

          CALL dbcsr_data_clear_pointer (data_block)
          t_blk = m_walltime()
          CALL dbcsr_get_block_p(redist, stored_row_new, stored_col_new,&
               data_block, tr, found)
          t_blk2 = t_blk2 + (m_walltime() - t_blk)
          valid_block = found

          IF (found) cnt_fnd = cnt_fnd + 1
          IF (.NOT. found .AND. .NOT. my_keep_sparsity) THEN
             ! We have to set up a buffer block
             CALL dbcsr_data_set_pointer (data_block,&
                  rsize=row_blk_size_new (stored_row_new),&
                  csize=col_blk_size_new (stored_col_new),&
                  pointee=buff_data)
             CALL dbcsr_data_clear (data_block)
             !r2_dp => r2_dp_buff(1:row_blk_size_new (stored_row_new),&
             !     1:col_blk_size_new (stored_col_new))
             !r2_dp(:,:) = 0.0_dp
             tr = .FALSE.
             blks = blks + 1
             blk_ps = blk_ps + row_blk_size_new (stored_row_new) * &
                  col_blk_size_new (stored_col_new)
             valid_block = .TRUE.
             cnt_new = cnt_new + 1
          ENDIF
          nze_rle = row_rle*col_rle

          IF (valid_block) THEN
             CALL dbcsr_data_set_pointer (src_2d,&
                  rsize=row_rle, csize=col_rle,&
                  pointee=recv_data, source_lb=data_offset_l)

             IF (.NOT. tr) THEN
                CALL dbcsr_data_set (dst=data_block,&
                     lb=row_offset_new, data_size=row_rle,&
                     lb2=col_offset_new, data_size2=col_rle,&
                     src=src_2d,&
                     source_lb=1, source_lb2=1)
                !CALL dbcsr_block_copy (&
                !     r2_dp(&
                !     row_offset_new:row_offset_new+row_rle-1,&
                !     col_offset_new:col_offset_new+col_rle-1),&
                !     recv_data(data_offset_l:data_offset_l+nze_rle-1),&
                !     row_rle, col_rle)
             ELSE
                CALL dbcsr_block_transpose(&
                     dst=data_block,&
                     lb=col_offset_new, lb2=row_offset_new,&
                     src=src_2d,&
                     row_size=row_rle, col_size=col_rle)
                !CALL dbcsr_block_transpose (&
                !     r2_dp(&
                !     col_offset_new:col_offset_new+col_rle-1,&
                !     row_offset_new:row_offset_new+row_rle-1),&
                !     recv_data(data_offset_l:data_offset_l+nze_rle-1),&
                !     row_rle, col_rle)
             ENDIF
          ELSE
             cnt_skip = cnt_skip+1
          ENDIF

          data_offset_l = data_offset_l + nze_rle
          IF (.NOT. found .AND. valid_block) THEN
             IF (dbg) WRITE(*,*)routineN//" Adding new block at",&
                  stored_row_new, stored_col_new
             CALL dbcsr_put_block(redist, stored_row_new, stored_col_new,&
                  data_block, tr)
             !DEALLOCATE (r2_dp)
          ELSE
             IF (.NOT. my_keep_sparsity .AND. dbg) &
                  WRITE(*,*)routineN//" Reusing block at",&
                  stored_row_new, stored_col_new
          ENDIF
       ENDDO
    ENDDO

    CALL dbcsr_data_clear_pointer(data_block)
    CALL dbcsr_data_clear_pointer(src_2d)
    CALL dbcsr_data_release (buff_data)
    CALL dbcsr_data_release (data_block)
    CALL dbcsr_data_release (src_2d)
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

    CALL dbcsr_data_release (recv_data)
    CALL dbcsr_data_release (send_data)

    DEALLOCATE(recv_meta)
    DEALLOCATE(send_meta)
    DBG 'Done redistributing'
    t_unpack = m_walltime() - t_unpack
    t_all = MAX(m_walltime() - t_all, 0.001_dp)
    !write(*,'(A,2(1X,F9.3),F5.1,"%",4(1X,I7))')&
    !     " Times total, lookup", t_all, t_blk2,&
    !     100.0_dp*t_blk2/t_all,&
    !     cnt_fnd, cnt_new, cnt_skip, redist%m%nblks
    !write(*,'(A,6(1X,F9.3))')" Times              ",&
    !     t_prolog, t_count, t_xfer1, t_fill, t_xfer2, t_unpack
    !write(*,*)" "
    !if (dbg) call dbcsr_print(redist)
    !call dbcsr_print(redist)
    !cs2 = dbcsr_checksum (redist)
    !WRITE(*,*)routineN//" Checksums=",cs1, cs2, cs1-cs2
    !CALL cp_assert (cs1-cs2 .LT. 0.00001, cp_fatal_level, cp_internal_error,&
    !     routineN, "Mangled data!")
    CALL timestop(timing_handle)
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

    CHARACTER                                :: data_type
    CHARACTER(len=80)                        :: ifmt
    INTEGER :: blk, blk_ps, blks, col, col_size, dst_p, meta_l, mp_group, &
      numnodes, nze, row, row_prev, row_size, src_p, stored_col_new, &
      stored_row_new, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rd_disp, recv_meta, rm_disp, &
                                                sd_disp, sdp, send_meta, &
                                                sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_count, send_count, &
                                                total_recv_count, &
                                                total_send_count
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size_new, &
                                                col_dist_new, &
                                                row_blk_size_new, row_dist_new
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_block, recv_data, &
                                                send_data
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
    data_type = matrix%m%data_type
    ! Get row and column start and end positions
    ! Old matrix
    dist_old = dbcsr_distribution (matrix)
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
       CALL dbcsr_iterator_next_block (iter, row, col, blk, tr,&
            row_size=row_size, col_size=col_size)
       CALL dbcsr_get_stored_coordinates (redist%m,&
            row, col, tr, dst_p)
       IF (dbg) WRITE(*,*)routinen//" Sending",row,col,"to",dst_p
       nze = row_size*col_size
       send_count(1, dst_p) = send_count(1, dst_p) + 1
       send_count(2, dst_p) = send_count(2, dst_p) + nze
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    DBG 'B'
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    DBG 'C'
    DBG 'send counts',send_count
    DBG 'recv counts',recv_count
    ! Allocate data structures needed for data exchange.
    CALL dbcsr_data_init (recv_data)
    CALL dbcsr_data_new (recv_data, data_type, SUM(recv_count(2, :)))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :))))
    CALL dbcsr_data_init (send_data)
    CALL dbcsr_data_new (send_data, data_type, SUM(send_count(2, :)))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :))))
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
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_new (data_block, data_type)
    CALL dbcsr_iterator_start (iter, matrix)
    dst_p = -1
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, data_block, tr, blk)
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
       nze = dbcsr_get_data_size (data_block)
       CALL dbcsr_data_set (send_data, lb=sdp(dst_p), data_size=nze,&
            src=data_block, source_lb=1)
       !send_data(sdp(dst_p):sdp(dst_p)+SIZE(r_dp)-1) &
       !     = r_dp(:)
       sdp(dst_p) = sdp(dst_p) + nze
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    CALL dbcsr_data_clear_pointer (data_block)
    ! Exchange the data and metadata structures.
    DBG 'E'
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       CALL hybrid_alltoall_s1(&
            send_data%d%r_sp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%r_sp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_real_8)
       !CALL mp_alltoall(&
       !     send_data%d%r_dp(:), total_send_count(2,:), sd_disp(:)-1,&
       !     recv_data%d%r_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
       !     mp_group)
       CALL hybrid_alltoall_d1 (&
            send_data%d%r_dp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%r_dp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_complex_4)
       CALL hybrid_alltoall_c1(&
            send_data%d%c_sp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%c_sp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_complex_8)
       CALL hybrid_alltoall_z1(&
            send_data%d%c_dp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%c_dp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    END SELECT
    !CALL mp_alltoall(send_data(:), total_send_count(2,:), sd_disp(:)-1,&
    !     recv_data(:), total_recv_count(2,:), rd_disp(:)-1, mp_group)
    DBG 'F'
    DBG 'send_meta'
    IF (dbg) WRITE(ifmt,'("(",I2,"(1X,I5))")') metalen
    IF (dbg) WRITE(*,ifmt) send_meta
    DBG 'send_data'
    !IF (dbg) WRITE(*,'(10(1X,F5.2))') send_data
    CALL hybrid_alltoall_i1(send_meta(:), SIZE(send_meta), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), SIZE(recv_meta), metalen*total_recv_count(1,:), rm_disp(:)-1, mp_obj_new)
    DBG 'recv_meta'
    IF (dbg) WRITE(*,ifmt) recv_meta
    DBG 'recv_data'
    !IF (dbg) WRITE(*,'(10(1X,F5.2))') recv_data
    DBG 'G'
    ! Now fill in the data.
    CALL dbcsr_work_create(redist,&
            SUM(recv_count(1,:)),&
            SUM(recv_count(2,:)), work_mutable=.FALSE., n=1)
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
          !r_dp => recv_data(blk_ps:blk_ps+nze-1)
          !CALL dbcsr_put_block(redist, stored_row_new, stored_col_new, r_dp, tr)
          !### this should be changed to be like the make images (i.e., copy data in finalize, not here & now)
          data_block = pointer_view (data_block, recv_data, blk_ps, nze)
          CALL dbcsr_put_block(redist, stored_row_new, stored_col_new, data_block, tr)
          blk_ps = blk_ps + nze
          blks = blks + 1
       ENDDO
    ENDDO
    IF (dbg) WRITE(*,*)routineN//" blk_p", redist%m%wms(1)%blk_p
    CALL dbcsr_data_clear_pointer (data_block)
    CALL dbcsr_data_release (data_block)
    !
    IF (dbg) THEN
       WRITE(*,*)routineN//" Declared blocks=",redist%m%wms(1)%lastblk,&
            "actual=",blks
       WRITE(*,*)routineN//" Declared data size=",redist%m%wms(1)%datasize,&
            "actual=",blk_ps
    ENDIF
    DBG 'Finalizing redistribution'
    CALL dbcsr_finalize(redist)
    CALL dbcsr_data_release (recv_data)
    CALL dbcsr_data_release (send_data)
    IF (dbg) CALL dbcsr_print (redist)
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_meta)
    DBG 'Done redistributing'
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_redistribute


! *****************************************************************************
!> \brief Transposes a DBCSR matrix.
!> \par Distribution options
!>      By default the distribution is transposed. If transpose_distribution
!>      is false, then an undetermined distribution is created that is
!>      compatible with the same process grid.
!> \param[out] transposed     transposed DBCSR matrix
!> \param[in] normal          input DBCSR matrix
!> \param[in] shallow_data_copy         (optional) only shallow data_copy;
!>                                      default is no; if set, the
!>                                      transpose_data option is ignored
!> \param[in] transpose_data  (optional) transpose data blocks, default is True
!> \param[in] transpose_distribution    (optional) transpose the distribution
!>                                      from the input matrix, default is True
!> \param[in] use_distribution          (optional) use this distribution
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_new_transposed (transposed, normal, shallow_data_copy,&
       transpose_data, transpose_distribution, use_distribution, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: transposed
    TYPE(dbcsr_obj), INTENT(IN)              :: normal
    LOGICAL, INTENT(IN), OPTIONAL            :: shallow_data_copy, &
                                                transpose_data, &
                                                transpose_distribution
    TYPE(dbcsr_distribution_obj), &
      INTENT(IN), OPTIONAL                   :: use_distribution
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_new_transposed', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: new_type
    INTEGER                                  :: stat, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: blk_p
    LOGICAL                                  :: shallow, tr_blocks, tr_dist
    TYPE(dbcsr_distribution_obj)             :: new_dist
    TYPE(dbcsr_obj)                          :: t2

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
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
    IF (PRESENT (use_distribution)) THEN
       CALL dbcsr_assert (dbcsr_distribution_nrows (use_distribution),&
            "EQ", dbcsr_distribution_ncols (normal%m%dist), cp_fatal_level,&
            cp_wrong_args_error, routineN,&
            "Given distribution must be compatible"&
            //" with the current distribution")
       CALL dbcsr_assert (dbcsr_distribution_ncols (use_distribution),&
            "EQ", dbcsr_distribution_nrows (normal%m%dist), cp_fatal_level,&
            cp_wrong_args_error, routineN,&
            "Given distribution must be compatible"&
            //" with the current distribution")
       new_dist = use_distribution
       CALL dbcsr_distribution_hold (new_dist)
    ELSE
       IF (tr_dist) THEN
          CALL dbcsr_transpose_distribution (new_dist, normal%m%dist)
       ELSE
          CALL dbcsr_transpose_dims (new_dist, normal%m%dist)
       ENDIF
    ENDIF
    ! Create the transposed matrix
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
            dbcsr_type_normal, transposed%m%row_blk_size, transposed%m%col_blk_size,&
            0, 0, transposed%m%data_type, error=error)
       CALL dbcsr_redistribute (transposed, t2)
       CALL dbcsr_release (transposed)
       transposed = t2
    ENDIF
    CALL timestop(timing_handle)
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
            SUM(recv_count(2,:)), error=error, n=1)
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
                                                nblks, numnodes, nze, p, row
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size, &
                                                tmp_index
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
    mp_obj = dbcsr_distribution_mp (dist)
    numnodes = dbcsr_mp_numnodes (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    mynode = dbcsr_mp_mynode (dbcsr_distribution_mp (dist))
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
            dbcsr_type_normal, matrix%m%row_blk_size, matrix%m%col_blk_size,&
            0, 0, matrix%m%data_type, special=matrix%m%special_memory)
       distributed%m%replication_type = dbcsr_repl_none
       ! First count how many blocks are local.
       nze = 0
       nblks = 0
       CALL dbcsr_iterator_start (iter, matrix)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, row, col, r_dp, tr, blk)
          tr = .FALSE.
          CALL dbcsr_get_stored_coordinates (matrix%m, row, col, tr, p)
          IF (mynode .EQ. p) THEN
             nze = nze + row_blk_size(row) * col_blk_size(col)
             nblks = nblks + 1
          ENDIF
       ENDDO
       CALL dbcsr_iterator_stop (iter)
       ! Preallocate the erray
       CALL dbcsr_work_create(distributed, nblks_guess=nblks,&
            sizedata_guess=nze, work_mutable=.FALSE.)
       ! Now actually do the work
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
!> \brief Gets a block from a dbcsr matrix as a data area
!> \par Data area
!>      The pointer encapsulated in the data area points to data stored in the
!>      matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get
!> \param[in] tr      whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
! *****************************************************************************
!> \brief Gets a 2-d block from a dbcsr matrix
!> \param[in]  matrix DBCSR matrix
!> \param[in]  row    the row
!> \param[in]  col    the column
!> \param[out] block  the block to get (rank-2 array)
!> \param[out] tr     whether the data is transposed
!> \param[out] found  whether the block exists in the matrix
!> \param[out] row_size      (optional) logical row size of block
!> \param[out] col_size      (optional) logical column size of block
! *****************************************************************************
  SUBROUTINE dbcsr_get_block_p_area(matrix,row,col,block,tr,found,&
       row_size, col_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: block
    LOGICAL, INTENT(OUT)                     :: tr, found
    INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_p_area', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, csize, iw, nze, offset, &
                                                rsize, stored_col, stored_row
    LOGICAL                                  :: stored_tr
    TYPE(btree_2d_data_c)                    :: data_block_c
    TYPE(btree_2d_data_d)                    :: data_block_d
    TYPE(btree_2d_data_s)                    :: data_block_s
    TYPE(btree_2d_data_z)                    :: data_block_z

!   ---------------------------------------------------------------------------

    CALL dbcsr_get_block_index (matrix, row, col, stored_row, stored_col,&
         stored_tr, found, blk, offset)

    tr = stored_tr

    rsize = dbcsr_blk_row_size (matrix%m, stored_row)
    csize = dbcsr_blk_column_size (matrix%m, stored_col)
    IF (PRESENT (row_size)) row_size = rsize
    IF (PRESENT (col_size)) col_size = csize

    CALL dbcsr_data_clear_pointer (block)
    IF(found) THEN
       nze = rsize*csize
       !IF(nze.eq.0) THEN
       !   found = .FALSE.
       !ELSE
       CALL dbcsr_data_set_pointer (block, rsize, csize, matrix%m%data_area,&
            source_lb = offset)
       !ENDIF
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       iw = 1
!$     iw = omp_get_thread_num()+1
       CALL cp_assert (dbcsr_use_mutable (matrix%m), cp_failure_level,&
            cp_caller_error, routineN,&
            "Can not retrieve blocks from non-mutable work matrices.")
       IF (dbcsr_mutable_instantiated(matrix%m%wms(iw)%mutable)) THEN
          SELECT CASE (block%d%data_type)
          CASE (dbcsr_type_real_4_2d)
             CALL btree_get_s (&
                  matrix%m%wms(iw)%mutable%m%btree_s,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block_s, found)
             IF (found) THEN
                CALL dbcsr_data_set_pointer (block, data_block_s%p)
             ENDIF
          CASE (dbcsr_type_real_8_2d)
             CALL btree_get_d (&
                  matrix%m%wms(iw)%mutable%m%btree_d,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block_d, found)
             IF (found) THEN
                CALL dbcsr_data_set_pointer (block, data_block_d%p)
             ENDIF
          CASE (dbcsr_type_complex_4_2d)
             CALL btree_get_c (&
                  matrix%m%wms(iw)%mutable%m%btree_c,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block_c, found)
             IF (found) THEN
                CALL dbcsr_data_set_pointer (block, data_block_c%p)
             ENDIF
          CASE (dbcsr_type_complex_8_2d)
             CALL btree_get_z (&
                  matrix%m%wms(iw)%mutable%m%btree_z,&
                  make_coordinate_tuple(stored_row, stored_col),&
                  data_block_z, found)
             IF (found) THEN
                CALL dbcsr_data_set_pointer (block, data_block_z%p)
             ENDIF
          CASE default
             CALL cp_assert (.FALSE., cp_fatal_level, cp_internal_error,&
                  routineN, "Invalid data type")
          END SELECT
       ENDIF
    ENDIF
  END SUBROUTINE dbcsr_get_block_p_area


! *****************************************************************************
!> \brief 
!>          We allow : 
!>                  matrix(dp) [+]= [scale(dp)] * block(dp)
!>                  matrix(dp) [+]= [scale(dp)] * block(sp)
!>                  matrix(sp) [+]= [scale(dp)] * block(sp)
!> \param[in]  
!> \param[out] 
!>
! *****************************************************************************
  SUBROUTINE dbcsr_put_block_area(matrix, row, col, block, transposed,&
       summation, scale)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    TYPE(dbcsr_data_obj)                     :: block
    LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: scale

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_put_block_area', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type_m
    LOGICAL                                  :: do_scale

    data_type_m = dbcsr_get_data_type (matrix)
    do_scale = PRESENT (scale)
    IF (do_scale) THEN
       !CALL cp_assert (data_type_m .EQ. scale%data_type, cp_fatal_level,&
       !     cp_wrong_args_error, routineN, "Incompatible data types matrix="//&
       !     data_type_m//" scale="//scale%data_type)
    ENDIF
    CALL cp_assert (ASSOCIATED (block%d), cp_fatal_level, cp_wrong_args_error,&
         routineN, "Can only add valid data block!")
    !CALL cp_assert (data_type_m .EQ. dbcsr_type_2d_to_1d(block%d%data_type),&
    !     cp_fatal_level, cp_wrong_args_error, routineN,&
    !     "Incompatible data types matrix="//&
    !        data_type_m//" block="//block%d%data_type)
    SELECT CASE(block%d%data_type)
    CASE (dbcsr_type_real_4)
       IF (do_scale) THEN
          IF(data_type_m.EQ.dbcsr_type_real_4) THEN
             CALL dbcsr_put_block(matrix, row, col, block%d%r_sp, transposed,&
                  summation, scale=scale%r_sp)
          ELSEIF(data_type_m.EQ.dbcsr_type_real_8) THEN
             CALL dbcsr_put_block(matrix, row, col, REAL(block%d%r_sp,real_8), transposed,&
                  summation, scale=REAL(scale%r_sp,real_8))
          ENDIF
       ELSE
          IF(data_type_m.EQ.dbcsr_type_real_4) THEN
             CALL dbcsr_put_block(matrix, row, col, block%d%r_sp, transposed,&
                  summation)
          ELSEIF(data_type_m.EQ.dbcsr_type_real_8) THEN
             CALL dbcsr_put_block(matrix, row, col, REAL(block%d%r_sp,real_8), transposed,&
                  summation)
          ENDIF
       ENDIF
    CASE (dbcsr_type_real_8)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%r_dp, transposed,&
               summation, scale=scale%r_dp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%r_dp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_complex_4)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%c_sp, transposed,&
               summation, scale=scale%c_sp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%c_sp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_complex_8)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%c_dp, transposed,&
               summation, scale=scale%c_dp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%c_dp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_real_4_2d)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%r2_sp, transposed,&
               summation, scale=scale%r_sp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%r2_sp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_real_8_2d)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%r2_dp, transposed,&
               summation, scale=scale%r_dp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%r2_dp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_complex_4_2d)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%c2_sp, transposed,&
               summation, scale=scale%c_sp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%c2_sp, transposed,&
               summation)
       ENDIF
    CASE (dbcsr_type_complex_8_2d)
       IF (do_scale) THEN
          CALL dbcsr_put_block(matrix, row, col, block%d%c2_dp, transposed,&
               summation, scale=scale%c_dp)
       ELSE
          CALL dbcsr_put_block(matrix, row, col, block%d%c2_dp, transposed,&
               summation)
       ENDIF
    CASE default
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error, routineN,&
            "Invalid data type")
    END SELECT
  END SUBROUTINE dbcsr_put_block_area


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

  SUBROUTINE dbcsr_trace_a_b_d(matrix_a, matrix_b, trace, trans_a, trans_b, local_sum, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a, matrix_b
    REAL(kind=real_8), INTENT(INOUT)         :: trace
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: trans_a, trans_b
    LOGICAL, INTENT(IN), OPTIONAL            :: local_sum
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_a_b_d', &
      routineP = moduleN//':'//routineN

    REAL(kind=real_4)                        :: trace_4

    IF(    dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_8 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_8 .OR. &
           dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_8 .OR. &
           dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_8 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_trace_ab_d(matrix_a, matrix_b, trace, trans_a, trans_b, local_sum, error)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_4) THEN
       trace_4 = 0.0_real_4
       CALL dbcsr_trace_ab_s(matrix_a, matrix_b, trace_4, trans_a, trans_b, local_sum, error)
       trace = REAL(trace_4,real_8)
    ELSE
       CALL cp_assert (.FALSE., cp_fatal_level, cp_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI")
    ENDIF
  END SUBROUTINE dbcsr_trace_a_b_d


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

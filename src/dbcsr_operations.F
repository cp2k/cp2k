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
  USE dbcsr_blas_operations,           ONLY: dbcsr_blas_copy,&
                                             dbcsr_blas_gemm,&
                                             dbcsr_blas_scal
  USE dbcsr_block_access,              ONLY: dbcsr_get_block_p,&
                                             dbcsr_put_block,&
                                             dbcsr_remove_block
  USE dbcsr_block_operations,          ONLY: dbcsr_block_partial_copy,&
                                             dbcsr_block_transpose,&
                                             dbcsr_data_clear
  USE dbcsr_data_operations,           ONLY: dbcsr_data_convert,&
                                             dbcsr_data_copyall
  USE dbcsr_dist_operations,           ONLY: checker_square_proc,&
                                             checker_tr,&
                                             dbcsr_find_column,&
                                             dbcsr_get_stored_coordinates
  USE dbcsr_error_handling
  USE dbcsr_index_operations,          ONLY: dbcsr_repoint_index
  USE dbcsr_io,                        ONLY: dbcsr_print
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_kinds,                     ONLY: dp,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8,&
                                             sp
  USE dbcsr_machine,                   ONLY: default_output_unit,&
                                             m_walltime
  USE dbcsr_message_passing,           ONLY: &
       dmp_max, mp_allgather, mp_bcast, mp_irecv, mp_isend, mp_sum, &
       mp_type_descriptor_type, mp_type_free, mp_type_make, mp_waitall
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_sizes, dbcsr_data_clear_pointer, &
       dbcsr_data_ensure_size, dbcsr_data_hold, dbcsr_data_init, &
       dbcsr_data_new, dbcsr_data_release, dbcsr_data_set_pointer, &
       dbcsr_destroy, dbcsr_destroy_array, dbcsr_distribution, &
       dbcsr_distribution_col_dist, dbcsr_distribution_local_rows, &
       dbcsr_distribution_mp, dbcsr_distribution_row_dist, dbcsr_get_data, &
       dbcsr_get_data_p, dbcsr_get_data_size, dbcsr_get_data_type, &
       dbcsr_get_matrix_type, dbcsr_init, dbcsr_is_initialized, &
       dbcsr_max_col_size, dbcsr_max_row_size, dbcsr_mp_grid_setup, &
       dbcsr_mp_group, dbcsr_mp_has_subgroups, dbcsr_mp_my_col_group, &
       dbcsr_mp_my_row_group, dbcsr_mp_mynode, dbcsr_mp_mypcol, &
       dbcsr_mp_myprow, dbcsr_mp_npcols, dbcsr_mp_nprows, dbcsr_mp_numnodes, &
       dbcsr_mp_pgrid, dbcsr_nblkcols_total, dbcsr_nblkrows_total, &
       dbcsr_nfullcols_local, dbcsr_nfullcols_total, dbcsr_nfullrows_local, &
       dbcsr_nfullrows_total, dbcsr_row_block_sizes, dbcsr_scalar, &
       dbcsr_scalar_are_equal, dbcsr_scalar_negative, dbcsr_scalar_one, &
       dbcsr_scalar_zero, dbcsr_valid_index
  USE dbcsr_ptr_util,                  ONLY: addto_array,&
                                             ensure_array_size,&
                                             pointer_replace,&
                                             pointer_view
  USE dbcsr_types,                     ONLY: &
       dbcsr_2d_array_type, dbcsr_data_obj, dbcsr_distribution_obj, &
       dbcsr_filter_frobenius, dbcsr_iterator, dbcsr_mp_obj, &
       dbcsr_norm_column, dbcsr_norm_frobenius, dbcsr_norm_gershgorin, &
       dbcsr_norm_maxabsnorm, dbcsr_obj, dbcsr_repl_full, dbcsr_repl_none, &
       dbcsr_scalar_type, dbcsr_type, dbcsr_type_antisymmetric, &
       dbcsr_type_complex_4, dbcsr_type_complex_8, dbcsr_type_real_4, &
       dbcsr_type_real_8, dbcsr_type_symmetric
  USE dbcsr_util,                      ONLY: absmod,&
                                             convert_sizes_to_offsets,&
                                             dbcsr_checksum,&
                                             dbcsr_set_debug,&
                                             swap
  USE dbcsr_work_operations,           ONLY: dbcsr_create,&
                                             dbcsr_finalize,&
                                             dbcsr_work_create

  !$ USE OMP_LIB
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_operations'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034


  PUBLIC :: dbcsr_mult_m_e_e, dbcsr_trace, dbcsr_add_on_diag,&
            dbcsr_set, dbcsr_scale, dbcsr_add, dbcsr_copy, dbcsr_get_diag, &
            dbcsr_set_diag, &
            dbcsr_get_block_diag, dbcsr_hadamard_product, &
            dbcsr_filter, dbcsr_scale_by_vector,&
            dbcsr_replace_blocks, dbcsr_conjg,&
            dbcsr_btriu, dbcsr_triu, dbcsr_tril, block_add_on_diag,&
            block_chol_inv, block_set, dbcsr_norm, &
            dbcsr_gershgorin_norm, dbcsr_maxabs_norm, dbcsr_frobenius_norm, &
            dbcsr_symmetrize_block_diag, dbcsr_copy_columns
  PUBLIC :: dbcsr_sum_replicated
  ! Direct access to the variants because array reshaping and
  ! typematching is incompatible.
  PUBLIC :: block_add_on_diag_d, block_add_on_diag_s,&
            block_add_on_diag_z, block_add_on_diag_c
  PUBLIC :: block_chol_inv_d, block_chol_inv_s,&
            block_chol_inv_z, block_chol_inv_c
  PUBLIC :: block_set_d, block_set_s,&
            block_set_z, block_set_c


! The interfaces for the generic routines found in the generated
! generic files.

  INTERFACE dbcsr_trace
     MODULE PROCEDURE dbcsr_trace_a_s, dbcsr_trace_a_d,&
          dbcsr_trace_a_c, dbcsr_trace_a_z,&
          dbcsr_trace_ab_s, dbcsr_trace_a_b_d,&!dbcsr_trace_ab_d,&
          dbcsr_trace_ab_c, dbcsr_trace_ab_z
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

  INTEGER, PARAMETER :: dgemm_stack_size = 32

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
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
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
    LOGICAL, PARAMETER                       :: use_combined_types = .FALSE.

    CHARACTER                                :: data_type
    INTEGER :: debug_level, error_handler, grp, i, ithread, left_col_image, &
      left_col_mult, left_col_nimages, left_col_v_offset, left_dst_icol, &
      left_dst_irow, left_dst_p, left_dst_pcol, left_dst_prow, &
      left_dst_v_col, left_dst_v_row, left_recv_icol, left_recv_irow, &
      left_recv_p, left_recv_pcol, left_recv_prow, left_recv_v_col, &
      left_recv_v_row, left_row_image, left_row_mult, left_row_nimages, &
      left_send_icol, left_send_irow, left_send_p, left_send_pcol, &
      left_send_prow, left_send_v_col, left_send_v_row, left_src_icol, &
      left_src_irow, left_src_p, left_src_pcol, left_src_prow, &
      left_src_v_col, left_src_v_row, max_nblks
    INTEGER :: max_nze, metronome, min_nimages, mp_group, mynode, nsteps_k, &
      nthreads, numnodes, nvirt_k, output_unit, right_col_image, &
      right_col_mult, right_col_nimages, right_dst_icol, right_dst_irow, &
      right_dst_p, right_dst_pcol, right_dst_prow, right_dst_v_col, &
      right_dst_v_row, right_recv_icol, right_recv_irow, right_recv_p, &
      right_recv_pcol, right_recv_prow, right_recv_v_col, right_recv_v_row, &
      right_row_image, right_row_mult, right_row_nimages, right_row_v_offset, &
      right_send_icol, right_send_irow, right_send_p, right_send_pcol, &
      right_send_prow, right_send_v_col, right_send_v_row, right_src_icol, &
      right_src_irow
    INTEGER :: right_src_p, right_src_pcol, right_src_prow, right_src_v_col, &
      right_src_v_row, timing_handle, v_k, v_ki
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
    TYPE(mp_type_descriptor_type), &
      ALLOCATABLE, DIMENSION(:, :)           :: left_recv_type, &
                                                left_send_type, &
                                                right_recv_type, &
                                                right_send_type
    TYPE(mp_type_descriptor_type), &
      DIMENSION(2)                           :: left_recv_subtypes, &
                                                left_send_subtypes, &
                                                right_recv_subtypes, &
                                                right_send_subtypes

!TYPE(cp_logger_type), POINTER            :: logger

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    debug_level=0
    t_all = 0.0_dp
    t_dgemm = 0.0_dp
    IF (dbg) debug_level = 2
    CALL timeset (routineN, timing_handle)
    ithread= 0 ; nthreads = 1
    !$ ithread = OMP_GET_THREAD_NUM () ; nthreads = OMP_GET_NUM_THREADS ()
    !IF (.NOT. ASSOCIATED (product_matrix%m%wms)) THEN
       CALL dbcsr_work_create(product_matrix,&
            work_mutable=.FALSE., error=error)
    !ENDIF
    !logger => cp_error_get_logger (error)
    output_unit = default_output_unit
    !cp_logger_get_default_io_unit (logger)
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
    CALL dbcsr_mp_grid_setup (left_mp_obj)
    CALL dbcsr_mp_grid_setup (right_mp_obj)
    IF (dbcsr_mp_has_subgroups (product_mp_obj)) THEN
       mp_rc_groups(M_P, 1:2) = (/ dbcsr_mp_my_row_group (product_mp_obj),&
            dbcsr_mp_my_col_group (product_mp_obj) /)
    ENDIF
    IF (dbcsr_mp_has_subgroups (left_mp_obj)) THEN
       mp_rc_groups(M_L, 1:2) = (/ dbcsr_mp_my_row_group (left_mp_obj),&
            dbcsr_mp_my_col_group (left_mp_obj) /)
    ENDIF
    IF (dbcsr_mp_has_subgroups (right_mp_obj)) THEN
       mp_rc_groups(M_R, 1:2) = (/ dbcsr_mp_my_row_group (right_mp_obj),&
            dbcsr_mp_my_col_group (right_mp_obj) /)
    ENDIF
    !CALL dbcsr_mp_grid_setup (product_mp_obj)
    !CALL dbcsr_assert (dbcsr_mp_has_subgroups (product_mp_obj),&
    !     dbcsr_fatal_level, dbcsr_internal_error, routineN,&
    !     "Row/column subgroups must be defined",__LINE__,error)
    !mp_rc_groups(M_P, 1:2) = (/ dbcsr_mp_my_row_group (product_mp_obj),&
    !                            dbcsr_mp_my_col_group (product_mp_obj) /)
    !CALL dbcsr_mp_grid_setup (left_mp_obj)
    !CALL dbcsr_assert (dbcsr_mp_has_subgroups (left_mp_obj),&
    !     dbcsr_fatal_level, dbcsr_internal_error, routineN,&
    !     "Row/column subgroups must be defined",__LINE__,error)
    !mp_rc_groups(M_L, 1:2) = (/ dbcsr_mp_my_row_group (left_mp_obj),&
    !                            dbcsr_mp_my_col_group (left_mp_obj) /)
    !CALL dbcsr_mp_grid_setup (right_mp_obj)
    !CALL dbcsr_assert (dbcsr_mp_has_subgroups (right_mp_obj),&
    !     dbcsr_fatal_level, dbcsr_internal_error, routineN,&
    !     "Row/column subgroups must be defined",__LINE__,error)
    !mp_rc_groups(M_R, 1:2) = (/ dbcsr_mp_my_row_group (right_mp_obj),&
    !                            dbcsr_mp_my_col_group (right_mp_obj) /)
    ! Dummy checks
    ! left/right matching
    CALL dbcsr_assert (left_col_nimages .EQ. right_row_mult&
         .AND. left_col_mult .EQ. right_row_nimages, dbcsr_fatal_level,&
         dbcsr_internal_error, routineN, "Left/Right image mismatch",__LINE__,error)
    CALL dbcsr_assert (left_col_nimages * dbcsr_mp_npcols (left_mp_obj) &
         .EQ. right_row_nimages * dbcsr_mp_nprows (right_mp_obj), &
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Left/Right total mismatch",__LINE__,error)
    ! product/left matching
    CALL dbcsr_assert (left_row_mult * dbcsr_mp_nprows (product_mp_obj) &
         .EQ. left_row_nimages * dbcsr_mp_nprows (left_mp_obj), &
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Product/Left total mismatch",__LINE__,error)
    ! product/left matching
    CALL dbcsr_assert (right_col_mult * dbcsr_mp_npcols (product_mp_obj) &
         .EQ. right_col_nimages * dbcsr_mp_npcols (right_mp_obj), &
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Product/Right total mismatch",__LINE__,error)
    ! Limitations
    CALL dbcsr_assert (left_row_nimages.eq.1 .AND. left_row_mult.eq.1,&
         dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
         "Product/Left matrix process grid mismatch",__LINE__,error)
    CALL dbcsr_assert (right_col_nimages.eq.1 .AND. right_col_mult.eq.1,&
         dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
         "Product/Right matrix process grid mismatch",__LINE__,error)
    !
    ! Setup product work areas
    DO i = 1, nthreads
       CALL dbcsr_data_ensure_size(product_matrix%m%wms(i)%data_area,&
            MAX(product_matrix%m%wms(i)%datasize,&
            dbcsr_nfullrows_local (product_matrix)&
            *dbcsr_nfullcols_local (product_matrix)),error=error)
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
    ALLOCATE (left_send_type (left_row_nimages, left_col_nimages))
    ALLOCATE (left_recv_type (left_row_nimages, left_col_nimages))
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
    ALLOCATE (right_send_type (right_row_nimages, right_col_nimages))
    ALLOCATE (right_recv_type (right_row_nimages, right_col_nimages))
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
          IF (use_combined_types) THEN
             DO v_ki = 1, right_row_nimages
                CALL mp_type_free (right_recv_type(v_ki, 1))
                CALL mp_type_free (right_send_type(v_ki, 1))
             ENDDO
          ELSE
             CALL mp_waitall (right_index_sr)
             CALL mp_waitall (right_index_rr)
          ENDIF
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
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Wrong right send prow",__LINE__,error)

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
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Wrong right dst prow",__LINE__,error)
             CALL dbcsr_assert(&
                  right_dst_irow, "LE", right_row_nimages,&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Wrong right dst row image",__LINE__,error)
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
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Wrong right receive prow",__LINE__,error)
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
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Wrong right src prow",__LINE__,error)
             CALL dbcsr_assert(&
                  right_src_irow, "LE", right_row_nimages,&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Wrong right src row image",__LINE__,error)
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
             IF (dbcsr_mp_has_subgroups (right_mp_obj)) THEN
                right_send_p = right_send_prow
                right_recv_p = right_recv_prow
                grp = dbcsr_mp_my_col_group (right_mp_obj)
             ELSE
                grp = dbcsr_mp_group (right_mp_obj)
             ENDIF
             !
             IF (use_combined_types) THEN
                right_send_subtypes(1) = dbcsr_mp_type_from_anytype (right_data_sp)
                right_send_subtypes(2) = mp_type_make (right_index_sp)
                right_recv_subtypes(1) = dbcsr_mp_type_from_anytype (right_data_rp)
                right_recv_subtypes(2) = mp_type_make (right_index_rp)
                right_send_type(v_ki+1, 1) = mp_type_make (right_send_subtypes)
                right_recv_type(v_ki+1, 1) = mp_type_make (right_recv_subtypes)
                CALL mp_isend (right_send_type(v_ki+1, 1), right_send_p,&
                     grp, right_data_sr(v_ki+1, 1), tag=right_dst_v_row)
                CALL mp_irecv (right_recv_type(v_ki+1, 1), right_recv_p,&
                     grp, right_data_rr(v_ki+1, 1), tag=right_src_v_row)
             ELSE
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
             ENDIF
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
          IF (use_combined_types) THEN
             DO v_ki = 1, left_col_nimages
                CALL mp_type_free (left_send_type(1, v_ki))
                CALL mp_type_free (left_recv_type(1, v_ki))
             ENDDO
          ELSE
             CALL mp_waitall (left_index_sr)
             CALL mp_waitall (left_index_rr)
          ENDIF
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
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Wrong left send pcol",__LINE__,error)
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
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Wrong left dst pcol",__LINE__,error)
             CALL dbcsr_assert(&
                  left_dst_icol, "LE", left_col_nimages,&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Wrong left dst col image",__LINE__,error)
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
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Wrong receive pcol",__LINE__,error)
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
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Wrong left src pcol",__LINE__,error)
             CALL dbcsr_assert(&
                  left_src_icol, "LE", left_col_nimages,&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Wrong left src column image",__LINE__,error)
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
             IF (dbcsr_mp_has_subgroups (left_mp_obj)) THEN
                left_send_p = left_send_pcol
                left_recv_p = left_recv_pcol
                grp = dbcsr_mp_my_row_group (left_mp_obj)
             ELSE
                grp = dbcsr_mp_group (left_mp_obj)
             ENDIF
             !
             IF (use_combined_types) THEN
                left_send_subtypes(1) = dbcsr_mp_type_from_anytype (left_data_sp)
                left_send_subtypes(2) = mp_type_make (left_index_sp)
                left_recv_subtypes(1) = dbcsr_mp_type_from_anytype (left_data_rp)
                left_recv_subtypes(2) = mp_type_make (left_index_rp)
                left_send_type(1, v_ki+1) = mp_type_make (left_send_subtypes)
                left_recv_type(1, v_ki+1) = mp_type_make (left_recv_subtypes)
                CALL mp_isend (left_send_type(1, v_ki+1), left_send_p,&
                     grp, left_data_sr(1, v_ki+1), tag=left_dst_v_col)
                CALL mp_irecv (left_recv_type(1, v_ki+1), left_recv_p,&
                     grp, left_data_rr(1, v_ki+1), tag=left_src_v_col)
             ELSE
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
             ENDIF
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
                CALL dbcsr_print(left_buffer_calc%mats(1, v_ki+1), nodata=.TRUE., error=error)
                CALL dbcsr_print(right_buffer_calc%mats(v_ki+1, 1), nodata=.TRUE., error=error)
             ENDIF
             tstart = m_walltime ()
             CALL dbcsr_nn_mult_lin(&
                  left_buffer_calc%mats(1, v_ki+1)%m,&
                  right_buffer_calc%mats(v_ki+1, 1)%m,&
                  product_matrix%m, submatrix_box=submatrix_selection,&
                  retain_sparsity=retain_sparsity,&
                  alpha=alpha,&
                  flop=flop_single,&
                  t_all=t_all, t_dgemm=t_dgemm, error=error)
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
    CALL dbcsr_destroy_array (right_buffer_1, error=error)
    CALL dbcsr_destroy_array (right_buffer_2, error=error)
    CALL dbcsr_destroy_array (left_buffer_1, error=error)
    CALL dbcsr_destroy_array (left_buffer_2, error=error)
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
    CALL dbcsr_finalize(product_matrix, error=error)
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
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_mult_m_e_e


  SUBROUTINE setup_buffer_matrices (buffer_set, buff_rows, buff_cols,&
       source_matrix, index_size, data_size, error)
    TYPE(dbcsr_2d_array_type), INTENT(OUT)   :: buffer_set
    INTEGER, INTENT(IN)                      :: buff_rows, buff_cols
    TYPE(dbcsr_obj), INTENT(IN)              :: source_matrix
    INTEGER, INTENT(IN)                      :: index_size, data_size
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_buffer_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col_image, error_handler, &
                                                row_image

    CALL dbcsr_error_set(routineN, error_handler, error)
    ALLOCATE (buffer_set%mats(buff_rows, buff_cols))
    DO row_image = 1, buff_rows
       DO col_image = 1, buff_cols
          CALL dbcsr_init(buffer_set%mats(row_image, col_image))
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
               data_size, nocopy=.TRUE.,error=error)
          CALL ensure_array_size (&
               buffer_set%mats(row_image, col_image)%m%index,&
               ub=index_size, nocopy=.TRUE., special=.TRUE.,&
               error=error)
          buffer_set%mats(row_image, col_image)%m%sym_negation&
               = source_matrix%m%sym_negation
       ENDDO
    ENDDO
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE setup_buffer_matrices


  FUNCTION dbcsr_mp_type_from_anytype(data_area) RESULT (mp_type)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: data_area
    TYPE(mp_type_descriptor_type)            :: mp_type

    SELECT CASE (data_area%d%data_type)
    CASE (dbcsr_type_real_4)
       mp_type = mp_type_make (data_area%d%r_sp)
    CASE (dbcsr_type_real_8)
       mp_type = mp_type_make (data_area%d%r_dp)
    CASE (dbcsr_type_complex_4)
       mp_type = mp_type_make (data_area%d%c_sp)
    CASE (dbcsr_type_complex_8)
       mp_type = mp_type_make (data_area%d%c_dp)
    END SELECT
  END FUNCTION dbcsr_mp_type_from_anytype


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
       t_all, t_dgemm, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    INTEGER, DIMENSION(6), INTENT(in)        :: submatrix_box
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: alpha, beta
    REAL(KIND=dp), INTENT(INOUT)             :: t_all, t_dgemm
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_nn_mult_lin', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER :: a_blk, a_col_logical, a_col_size_logical, a_row_l, &
      a_row_logical, a_row_size_logical, b_blk, b_col_logical, &
      b_col_size_logical, b_row_logical, b_row_size_logical, blk_end, &
      blk_size, blk_start, c_blk_p, c_blk_pt, c_col_logical, &
      c_col_size_logical, c_nze, c_row_logical, c_row_size_logical, &
      clamp_k_len, datasize, error_handler, first_k, first_mergeable, &
      first_row, i, ithread, j, last_k, last_row, lastblk, new_blk, &
      new_row_blks, nthreads, old_blk, old_row_begin_blk, old_size, &
      row_begin_blk, row_blk, rpt, size_c_blk_ps, stack_p, stack_p_a, &
      stack_p_b, stack_p_c, stat, submat_f_col, submat_f_row, submat_l_col
    INTEGER :: submat_l_row, this_right_limit
    INTEGER(KIND=int_4)                      :: offset
    INTEGER(KIND=int_8)                      :: lflop, one_flop
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: c_blk_ps, c_blk_ps_p
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: new_index
    INTEGER, DIMENSION(:), POINTER :: left_col_blk_offset, left_local_rows, &
      left_row_blk_offset, new_blk_p, new_col_i, new_row_i, &
      product_col_blk_offset, right_col_blk_offset
    LOGICAL :: block_exists, flush_stack, found, keep_sparsity, &
      left_limit_rows, limit_k, negate_alpha, product_is_symmetric, &
      right_limit_cols, tr_p
    REAL(KIND=dp)                            :: epoch, mt_t_dgemm
    REAL(KIND=real_8)                        :: one_time
    TYPE(block_parameters), &
      DIMENSION(dgemm_stack_size)            :: parameters_a, parameters_b, &
                                                parameters_c
    TYPE(dbcsr_scalar_type)                  :: beta_add, beta_new, my_alpha, &
                                                my_beta, which_beta
    TYPE(dgemm_join), &
      DIMENSION(dgemm_stack_size)            :: dgemm_joins

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
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
    CALL dbcsr_assert (submatrix_box(1) .EQ. 0,&
         dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
         "Limiting first row not supported",__LINE__,error)
    CALL dbcsr_assert (submatrix_box(3) .EQ. 0,&
         dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
         "Limiting first column not supported",__LINE__,error)
    CALL dbcsr_assert (submatrix_box(5) .EQ. 0,&
         dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
         "Limiting first k not supported",__LINE__,error)
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
    one_time = 0.0_dp
    old_blk = 1
    first_mergeable = 1
    stack_p = 0 ; stack_p_a = 0 ; stack_p_b = 0 ; stack_p_c = 0
    new_blk = 0
    !write(*,*)ithread,':',left_thread_dist
    rpt = CEILING (REAL(left%nblkrows_local) / REAL(nthreads))
    first_row = MIN (ithread*rpt+1, left%nblkrows_local)
    last_row = MIN ((ithread+1)*rpt, left%nblkrows_local)
    !write(*,*)ithread,'Doing rows', first_row, last_row
    DO a_row_l = 1, left%nblkrows_local
       a_row_logical = left_local_rows (a_row_l)
       IF (dbg) WRITE(*,*)routineN//" A row", a_row_logical
       blk_start = left_row_blk_offset(a_row_logical)
       blk_end = left_row_blk_offset(a_row_logical+1)-1
       blk_size = left_row_blk_offset(a_row_logical+1) - blk_start
       IF (left_limit_rows) THEN
          IF (left_row_blk_offset (a_row_logical) .GT. submat_l_row .OR.&
              left_row_blk_offset (a_row_logical+1)-1 .LT. submat_f_row) THEN
             CYCLE
          ELSE
             CALL dbcsr_assert(&
                    MIN(submat_l_row, blk_start)&
                  - MAX(submat_f_row, blk_end)+1&
                  .EQ. blk_size,&
                  dbcsr_fatal_level, dbcsr_unimplemented_error_nr,&
                  routineN, "Submatrix must encompass complete block rows",__LINE__,error)
          ENDIF
       ENDIF
       a_row_size_logical = blk_size ! stored
       c_row_size_logical = blk_size
       row_begin_blk = new_blk+1
       IF (old_blk .LE. old_size) THEN
          DO WHILE (product%wms(ithread+1)%row_i(old_blk) .LT. a_row_logical)
             old_blk = old_blk+1
             IF (old_blk .GT. old_size) THEN
                !old_blk = old_blk-1
                EXIT
             ENDIF
          ENDDO
       ENDIF
       old_row_begin_blk = old_blk
       ! Clear the c_blk_ps array
       DO i = 1, size_c_blk_ps
          c_blk_ps(c_blk_ps_p(i)) = 0
       ENDDO
       size_c_blk_ps = 0
       ! Fill in the c_blk_ps array
       IF (dbg) WRITE(*,*)routineN//" old_blk, old_size, old_row_begin_blk",&
            old_blk,old_size,old_row_begin_blk
       IF (old_blk .LE. old_size) THEN
          ! Old index
          DO WHILE (product%wms(ithread+1)%row_i(old_blk) .LE. a_row_logical) ! EQ
             DBG "Adding c_blk_ps column",&
                  product%wms(ithread+1)%col_i(old_blk),&
                  product%wms(ithread+1)%blk_p(old_blk)
             c_col_logical = product%wms(ithread+1)%col_i(old_blk)
             found = .FALSE.
             c_blk_ps(c_col_logical)&
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
       c_row_logical = a_row_logical
       new_row_blks = 0
       DO a_blk = left%row_p(a_row_logical)+1, left%row_p(a_row_logical+1)
          IF (left%blk_p(a_blk) .EQ. 0) THEN
             CYCLE
          ENDIF
          row_blk = row_begin_blk
          old_blk = old_row_begin_blk
          a_col_logical = left%col_i(a_blk)
          IF (dbg) WRITE(*,*)routineN//" A col", a_col_logical,";",a_row_logical
          blk_start = left_col_blk_offset (a_col_logical)
          blk_end = left_col_blk_offset (a_col_logical+1)-1
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
          a_col_size_logical = left_col_blk_offset (a_col_logical+1) - blk_start
          stack_p_a = stack_p_a + 1
          parameters_a(stack_p_a)%offset = ABS(left%blk_p(a_blk))
          parameters_a(stack_p_a)%tr = left%blk_p(a_blk) .LT. 0
          parameters_a(stack_p_a)%logical_rows = a_row_size_logical
          parameters_a(stack_p_a)%logical_cols = a_col_size_logical
          parameters_a(stack_p_a)%nze = a_row_size_logical*a_col_size_logical
          !
          b_row_logical = a_col_logical
          b_row_size_logical = a_col_size_logical
          DO b_blk = right%row_p(b_row_logical)+1, right%row_p(b_row_logical+1)
             IF (right%blk_p(b_blk) .EQ. 0) THEN
                CYCLE
             ENDIF
             b_col_logical = right%col_i(b_blk)
             c_col_logical = b_col_logical
             block_exists = c_blk_ps(c_col_logical) .NE. 0
             sparsity_enforcement: IF &
                  (keep_sparsity .AND. .NOT. block_exists) THEN
                CYCLE
             ENDIF sparsity_enforcement
             ! Don't calculate symmetric blocks.
             IF (product_is_symmetric) THEN
                IF (c_row_logical .NE. c_col_logical&
                     .AND. checker_tr (c_row_logical, c_col_logical)) THEN
                   !write(*,*)"Skipping symmetric block!", c_row_logical,&
                   !     c_col_logical
                   CYCLE
                ENDIF
             ENDIF
             !IF (dbg) WRITE(*,*)routineN//" B col", b_col_logical, ";",a_row_logical, a_col_logical
             blk_start = right_col_blk_offset (b_col_logical)
             blk_end = right_col_blk_offset (b_col_logical+1)-1
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
             b_col_size_logical = right_col_blk_offset(b_col_logical+1) - blk_start
             offset = ABS(c_blk_ps(c_col_logical))
             tr_p = c_blk_ps(c_col_logical) .LT. 0
             stack_p_b = stack_p_b + 1
             parameters_b(stack_p_b)%offset = ABS(right%blk_p(b_blk))
             parameters_b(stack_p_b)%tr = right%blk_p(b_blk) .LT. 0

             parameters_b(stack_p_b)%logical_rows = b_row_size_logical
             parameters_b(stack_p_b)%logical_cols = b_col_size_logical
             parameters_b(stack_p_b)%nze = b_row_size_logical * b_col_size_logical
             !
             c_col_size_logical = &
                  product_col_blk_offset (c_col_logical+1) -&
                  product_col_blk_offset (c_col_logical)
             c_nze = c_row_size_logical * c_col_size_logical
             !
             !
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
                !write(*,*)'Adding coordinate', c_row_logical, c_col_logical, c_blk_pt
                c_blk_pt = c_blk_pt + c_nze
                datasize = datasize + c_nze
                lastblk = lastblk + 1
                c_blk_ps (c_col_logical) = offset
                size_c_blk_ps = size_c_blk_ps + 1
                c_blk_ps_p(size_c_blk_ps) = c_col_logical
                which_beta = beta_new
                tr_p = offset .LT. 0
             ELSE
                !DBG 'Using existing coordinate for', c_row_logical, c_col_logical, offset, tr_p
                !write(*,*) 'Using existing coordinate for', c_row_logical, c_col_logical, offset, tr_p
                which_beta = beta_add
             ENDIF
             stack_p_c = stack_p_c + 1
             parameters_c(stack_p_c)%logical_rows = c_row_size_logical
             parameters_c(stack_p_c)%logical_cols = c_col_size_logical
             parameters_c(stack_p_c)%nze = c_nze
             parameters_c(stack_p_c)%tr = tr_p
             parameters_c(stack_p_c)%offset = ABS(offset)
             !write(*,*)'@@@',stack_p_c, parameters_c(stack_p_c)%offset
             !write(*,*)'@@@',stack_p_c, c_row_size_logical, c_col_size_logical, c_nze
             !
             stack_p = stack_p + 1
             dgemm_joins(stack_p)%p_a = stack_p_a
             dgemm_joins(stack_p)%p_b = stack_p_b
             dgemm_joins(stack_p)%p_c = stack_p_c
             dgemm_joins(stack_p)%beta = which_beta
             IF (.NOT. limit_k) THEN
                dgemm_joins(stack_p)%last_k = b_row_size_logical
             ELSE
                dgemm_joins(stack_p)%last_k = clamp_k_len
             ENDIF
             IF (right_limit_cols .AND. this_right_limit.gt.0) THEN
                dgemm_joins(stack_p)%last_n = this_right_limit
             ELSE
                dgemm_joins(stack_p)%last_n = b_col_size_logical
                ! assert b_col_size_logical .eq. c_col_size_logical
             ENDIF
             negate_alpha = parameters_a(stack_p_a)%tr .AND. left%sym_negation
             negate_alpha = negate_alpha .NEQV. (parameters_b(stack_p_b)%tr .AND. right%sym_negation)
             !write(*,*)routineN//" tr_p:", tr_p
             !write(*,*)routineN//" p.s.n", product%sym_negation
             negate_alpha = negate_alpha .NEQV. (tr_p .AND. product%sym_negation)
             IF (negate_alpha) THEN
                dgemm_joins(stack_p)%alpha &
                     = dbcsr_scalar_negative (my_alpha)
             ELSE
                dgemm_joins(stack_p)%alpha = my_alpha
             ENDIF
             !i = SIZE (dgemm_joins,1) / 4
             i = dgemm_stack_size / 8
             !IF (dbg) WRITE(*,*)routineN//" Before merge, stack pointer is at",stack_p
             !IF (stack_p .GE. i .and. stack_p .gt. first_mergeable) THEN
             IF (stack_p .GE. i + first_mergeable .OR. stack_p .GE. dgemm_stack_size) THEN
                CALL merge_dgemm_stack(dgemm_joins,&
                     parameters_a, parameters_b, parameters_c,&
                     first_mergeable,&
                     stack_p,&
                     stack_p_a, stack_p_b, stack_p_c)
             ENDIF
             IF (dbg) WRITE(*,*)routineN//"  After merge, stack pointer is at",stack_p
             flush_stack = (stack_p .GE. dgemm_stack_size-1)
             IF (flush_stack) THEN
                !IF (dbg) WRITE(*,*)routineN//" Processing stack at", stack_p
                !mt_t_dgemm = mt_t_dgemm - (m_walltime() - epoch)
                CALL process_dgemm_stack(dgemm_joins,&
                     parameters_a, parameters_b, parameters_c,&
                     stack_p,&
                     stack_p_a, stack_p_b, stack_p_c,&
                     left%data_area, right%data_area,&
                     product%wms(ithread+1)%data_area,&
                     !lflop=one_flop,ltime=one_time,error=error)
                     lflop=one_flop,error=error)
                !mt_t_dgemm = mt_t_dgemm + (m_walltime() - epoch)
                mt_t_dgemm = mt_t_dgemm + one_time
                lflop = lflop + one_flop
                first_mergeable = 1
                stack_p = 0
                stack_p_b = 0
                stack_p_c = 0
                !parameters_a(1) = parameters_a(stack_p_a)
                !stack_p_a = 1
             ELSE
                !WRITE(*,*)routineN//" Skipping DGEMM"
                !stack_p = stack_p+1
                !stack_p_b = stack_p_b+1
                !stack_p_c = stack_p_c+1
             ENDIF
          ENDDO ! b
          !@@@
          !write(*,*)routineN//" End of ROW"
          ! Conclude the row, merge everything, and prepare for the new row.
          IF (stack_p_a .GE. SIZE (parameters_a)) THEN
             !stack_p = stack_p - 1
             !stack_p_b = stack_p_b - 1
             !stack_p_c = stack_p_c - 1
             CALL merge_dgemm_stack(dgemm_joins,&
                  parameters_a, parameters_b, parameters_c,&
                  first_mergeable,&
                  stack_p,&
                  stack_p_a, stack_p_b, stack_p_c)
             CALL process_dgemm_stack(dgemm_joins,&
                  parameters_a, parameters_b, parameters_c,&
                  stack_p,&
                  stack_p_a, stack_p_b, stack_p_c,&
                  left%data_area, right%data_area,&
                  product%wms(ithread+1)%data_area,&
                  !lflop=one_flop,ltime=one_time,error=error)
                  lflop=one_flop,error=error)
             mt_t_dgemm = mt_t_dgemm + one_time
             stack_p = 0
             stack_p_a = 0
             stack_p_b = 0
             stack_p_c = 0
             first_mergeable = 1
          ENDIF
          first_mergeable = MAX (stack_p, 1)
       ENDDO ! a_col
    ENDDO ! a_row
    IF (stack_p .GE. 1) THEN
       IF (dbg) WRITE(*,*)routineN//" Final processing at", stack_p
       CALL merge_dgemm_stack(dgemm_joins,&
            parameters_a, parameters_b, parameters_c,&
            first_mergeable,&
            stack_p,&
            stack_p_a, stack_p_b, stack_p_c)
       !mt_t_dgemm = mt_t_dgemm - (m_walltime() - epoch)
       CALL process_dgemm_stack(dgemm_joins,&
            parameters_a, parameters_b, parameters_c,&
            stack_p,&
            stack_p_a, stack_p_b, stack_p_c,&
            left%data_area, right%data_area,&
            product%wms(ithread+1)%data_area,&
            !lflop=one_flop,ltime=one_time,error=error)
            lflop=one_flop,error=error)
       !mt_t_dgemm = mt_t_dgemm + (m_walltime() - epoch)
       mt_t_dgemm = mt_t_dgemm + one_time
       lflop = lflop + one_flop
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
         new_index, new_blk, error=error)
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
    CALL dbcsr_error_stop(error_handler, error)
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
  PURE SUBROUTINE merge_dgemm_stack(params, params_a, params_b, params_c,&
       first_mergeable,&
       stack_size,&
       stack_size_a, stack_size_b, stack_size_c)
    INTEGER, INTENT(INOUT)                   :: stack_size, first_mergeable, &
                                                stack_size_c, stack_size_b, &
                                                stack_size_a
    TYPE(dgemm_join), &
      DIMENSION(1:stack_size), INTENT(INOUT) :: params
    TYPE(block_parameters), &
      DIMENSION(1:stack_size_a), &
      INTENT(INOUT)                          :: params_a
    TYPE(block_parameters), &
      DIMENSION(1:stack_size_b), &
      INTENT(INOUT)                          :: params_b
    TYPE(block_parameters), &
      DIMENSION(1:stack_size_c), &
      INTENT(INOUT)                          :: params_c

    CHARACTER(len=*), PARAMETER :: routineN = 'merge_dgemm_stack', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: f, last, sp, sp_a, sp_b, sp_c
    LOGICAL                                  :: merge_this, remainder

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY

    last = stack_size
    f = first_mergeable
    !IF (dbg) THEN
    !   WRITE(*,*)routineN//" merge attempt",f,"to",last,";",&
    !        stack_size_a, stack_size_b, stack_size_c
    !ENDIF
    !call dbcsr_assert (last, "LE", size(params), cp_fatal_level,&
    !     cp_internal_error, routineN, "sp too large")
    ! This is used to flag merging of blocks prior to stack.
    remainder = .FALSE.
    DO sp = last, f+1, -1
       !IF (dbg) THEN
       !   WRITE(*,*)routineN//" testing merge of", sp-1, sp
       !   CALL print_dgemm_parameters(params(sp-1:sp),&
       !        params_a, params_b, params_c)
       !ENDIF
       sp_a = params(sp)%p_a
       sp_b = params(sp)%p_b
       sp_c = params(sp)%p_c
       merge_this = sp_b .GT. 1 .AND. sp_c .GT. 1
       !call dbcsr_assert (sp_b, "LE", stack_size_b, cp_fatal_level,&
       !     cp_internal_error, routineN, "spb too large")
       !call dbcsr_assert (sp_c, "LE", stack_size_c, cp_fatal_level,&
       !     cp_internal_error, routineN, "spc too large")
       IF (merge_this) THEN
          merge_this = &
                  (.NOT. params_b(sp_b-1)%tr) .AND. (.NOT. params_b(sp_b)%tr) &
            .AND. (.NOT. params_c(sp_c-1)%tr) .AND. (.NOT. params_c(sp_c)%tr)
       ENDIF
       IF (merge_this) THEN
          merge_this = &
                        params_b(sp_b-1)%logical_cols .EQ. params(sp-1)%last_n &
                  .AND. params_c(sp_c-1)%logical_cols .EQ. params(sp-1)%last_n &
                  .AND. params(sp-1)%p_a .EQ. sp_a &
                  .AND. params_b(sp_b-1)%offset + params_b(sp_b-1)%nze .EQ.&
                        params_b(sp_b)%offset &
                  .AND. params_c(sp_c-1)%offset + params_c(sp_c-1)%nze .EQ.&
                        params_c(sp_c)%offset
       ENDIF
       IF (merge_this) THEN
          merge_this = &
               dbcsr_scalar_are_equal(params(sp-1)%beta, params(sp)%beta) &
               .AND. dbcsr_scalar_are_equal(params(sp-1)%alpha, params(sp)%alpha)
       ENDIF
       IF (merge_this) THEN
          !IF (dbg) WRITE(*,*)routineN//" Merging into:"
          params(sp-1)%last_n = params(sp-1)%last_n + params(sp)%last_n
          params(sp)%last_k = 0
          params_b(sp_b-1)%logical_cols =&
               params_b(sp_b-1)%logical_cols + params_b(sp_b)%logical_cols
          params_b(sp_b-1)%nze = params_b(sp_b-1)%nze + params_b(sp_b)%nze
          params_c(sp_c-1)%logical_cols =&
               params_c(sp_c-1)%logical_cols + params_c(sp_c)%logical_cols
          params_c(sp_c-1)%nze = params_c(sp_c-1)%nze + params_c(sp_c)%nze
          !IF (dbg) CALL print_dgemm_parameters(params(sp-1:sp-1),&
          !     params_a, params_b, params_c)
          IF (.not.remainder) THEN
             stack_size = stack_size - 1
             stack_size_b = stack_size_b - 1
             stack_size_c = stack_size_c - 1
          ENDIF
       ELSE
          remainder = .TRUE.
          IF (sp .GE. first_mergeable) THEN
             first_mergeable = sp
          ENDIF
          !IF (dbg) WRITE(*,*)routineN//" New first is ",first_mergeable
       ENDIF
    ENDDO
    !IF (dbg) THEN
    !   WRITE(*,*)routineN//" The new first and n are ",&
    !        first_mergeable, stack_size,";",&
    !        stack_size_a, stack_size_b, stack_size_c
    !ENDIF
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
!> \param[out] ltime           (optional) Time used by DGEMM
! *****************************************************************************
  SUBROUTINE process_dgemm_stack(params, params_a, params_b, params_c,&
       stack_size, stack_size_a, stack_size_b, stack_size_c,&
       left_data_area, right_data_area, product_data_area, lflop, ltime, error)
    INTEGER, INTENT(INOUT)                   :: stack_size, stack_size_c, &
                                                stack_size_b, stack_size_a
    TYPE(dgemm_join), &
      DIMENSION(1:stack_size), INTENT(IN)    :: params
    TYPE(block_parameters), &
      DIMENSION(1:stack_size_a), INTENT(IN)  :: params_a
    TYPE(block_parameters), &
      DIMENSION(1:stack_size_b), INTENT(IN)  :: params_b
    TYPE(block_parameters), &
      DIMENSION(1:stack_size_c), INTENT(IN)  :: params_c
    TYPE(dbcsr_data_obj), INTENT(IN)         :: left_data_area, &
                                                right_data_area
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: product_data_area
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: lflop
    REAL(KIND=real_8), INTENT(OUT), OPTIONAL :: ltime
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'process_dgemm_stack', &
      routineP = moduleN//':'//routineN
    REAL, PARAMETER                          :: resize_factor = 1.618034

    INTEGER                                  :: error_handler

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    !WRITE(*,*)routineN//"========== Mulitplying"
    !CALL print_dgemm_parameters (params(1:stack_size),&
    !     params_a, params_b, params_c)
    !WRITE(*,*)routineN//"========== END of multiplies"
    SELECT CASE (product_data_area%d%data_type)
    CASE (dbcsr_type_real_4)
       CALL process_dgemm_stack_s (params, params_a, params_b, params_c,&
       stack_size, stack_size_a, stack_size_b, stack_size_c,&
       left_data_area%d%r_sp, right_data_area%d%r_sp, product_data_area,&
       lflop=lflop, ltime=ltime, error=error)
    CASE (dbcsr_type_real_8)
       CALL process_dgemm_stack_d (params, params_a, params_b, params_c,&
       stack_size, stack_size_a, stack_size_b, stack_size_c,&
       left_data_area%d%r_dp, right_data_area%d%r_dp, product_data_area,&
       lflop=lflop, ltime=ltime, error=error)
    CASE (dbcsr_type_complex_4)
       CALL process_dgemm_stack_c (params, params_a, params_b, params_c,&
       stack_size, stack_size_a, stack_size_b, stack_size_c,&
       left_data_area%d%c_sp, right_data_area%d%c_sp, product_data_area,&
       lflop=lflop, ltime=ltime, error=error)
    CASE (dbcsr_type_complex_8)
       CALL process_dgemm_stack_z (params, params_a, params_b, params_c,&
       stack_size, stack_size_a, stack_size_b, stack_size_c,&
       left_data_area%d%c_dp, right_data_area%d%c_dp, product_data_area,&
       lflop=lflop, ltime=ltime, error=error)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
    stack_size = 0
    CALL dbcsr_error_stop(error_handler, error)
    RETURN
  END SUBROUTINE process_dgemm_stack

  SUBROUTINE print_dgemm_parameters(params, params_a, params_b, params_c)
    TYPE(dgemm_join), DIMENSION(:), &
      INTENT(in)                             :: params
    TYPE(block_parameters), DIMENSION(:), &
      INTENT(IN)                             :: params_a, params_b, params_c

    INTEGER                                  :: sp, sp_a, sp_b, sp_c

    DO sp = 1, SIZE(params)
       sp_a = params(sp)%p_a
       sp_b = params(sp)%p_b
       sp_c = params(sp)%p_c
       WRITE(*,*)"SP",sp, sp_a, sp_b, sp_c
       SELECT CASE (params(sp)%alpha%data_type)
       CASE (dbcsr_type_real_4)
          WRITE(*,'(1X,A,1X,L1,1X,2(L1),3(1X,2(1X,I4),";"),".",2(1X,I4),2(1X,F5.1),/,10X,"@",3(1X,2(1X,I7),";"))')&
               "DGEMM PARAMETERS",&
               params_c(sp_c)%tr,&
               params_a(sp_a)%tr, params_b(sp_b)%tr,&
               params_a(sp_a)%logical_rows, params_a(sp_a)%logical_cols,&
               params_b(sp_b)%logical_rows, params_b(sp_b)%logical_cols,&
               params_c(sp_c)%logical_rows, params_c(sp_c)%logical_cols,&
               params(sp)%last_k, params(sp)%last_n,&
               params(sp)%alpha%r_sp, params(sp)%beta%r_sp,&
               params_a(sp_a)%offset, params_a(sp_a)%nze,&
               params_b(sp_b)%offset, params_b(sp_b)%nze,&
               params_c(sp_c)%offset, params_c(sp_c)%nze
       CASE (dbcsr_type_real_8)
          WRITE(*,'(1X,A,1X,L1,1X,2(L1),3(1X,2(1X,I4),";"),".",2(1X,I4),2(1X,F5.1),/,10X,"@",3(1X,2(1X,I7),";"))')&
               "DGEMM PARAMETERS",&
               params_c(sp_c)%tr,&
               params_a(sp_a)%tr, params_b(sp_b)%tr,&
               params_a(sp_a)%logical_rows, params_a(sp_a)%logical_cols,&
               params_b(sp_b)%logical_rows, params_b(sp_b)%logical_cols,&
               params_c(sp_c)%logical_rows, params_c(sp_c)%logical_cols,&
               params(sp)%last_k, params(sp)%last_n,&
               params(sp)%alpha%r_dp, params(sp)%beta%r_dp,&
               params_a(sp_a)%offset, params_a(sp_a)%nze,&
               params_b(sp_b)%offset, params_b(sp_b)%nze,&
               params_c(sp_c)%offset, params_c(sp_c)%nze
       CASE (dbcsr_type_complex_4)
          WRITE(*,'(1X,A,1X,L1,1X,2(L1),3(1X,2(1X,I4),";"),".",2(1X,I4),2(1X,F5.1),/,10X,"@",3(1X,2(1X,I7),";"))')&
               "DGEMM PARAMETERS",&
               params_c(sp_c)%tr,&
               params_a(sp_a)%tr, params_b(sp_b)%tr,&
               params_a(sp_a)%logical_rows, params_a(sp_a)%logical_cols,&
               params_b(sp_b)%logical_rows, params_b(sp_b)%logical_cols,&
               params_c(sp_c)%logical_rows, params_c(sp_c)%logical_cols,&
               params(sp)%last_k, params(sp)%last_n,&
               REAL(params(sp)%alpha%c_sp), REAL(params(sp)%beta%c_sp),&
               params_a(sp_a)%offset, params_a(sp_a)%nze,&
               params_b(sp_b)%offset, params_b(sp_b)%nze,&
               params_c(sp_c)%offset, params_c(sp_c)%nze
       CASE (dbcsr_type_complex_8)
          WRITE(*,'(1X,A,1X,L1,1X,2(L1),3(1X,2(1X,I4),";"),".",2(1X,I4),2(1X,F5.1),/,10X,"@",3(1X,2(1X,I7),";"))')&
               "DGEMM PARAMETERS",&
               params_c(sp_c)%tr,&
               params_a(sp_a)%tr, params_b(sp_b)%tr,&
               params_a(sp_a)%logical_rows, params_a(sp_a)%logical_cols,&
               params_b(sp_b)%logical_rows, params_b(sp_b)%logical_cols,&
               params_c(sp_c)%logical_rows, params_c(sp_c)%logical_cols,&
               params(sp)%last_k, params(sp)%last_n,&
               REAL(params(sp)%alpha%c_dp), REAL(params(sp)%beta%c_dp),&
               params_a(sp_a)%offset, params_a(sp_a)%nze,&
               params_b(sp_b)%offset, params_b(sp_b)%nze,&
               params_c(sp_c)%offset, params_c(sp_c)%nze
       END SELECT
    ENDDO
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
    LOGICAL, INTENT(in)                      :: t
    INTEGER, INTENT(in)                      :: n1, n2
    INTEGER                                  :: val

    IF (.NOT. t) THEN
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
  END SUBROUTINE merge_arrays



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
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
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
      col_offset, error_handler, i, last_col, row_offset
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: col_end, col_start, m_offset
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: do_scale, limit_col, right, tr
    TYPE(dbcsr_data_obj)                     :: buffer, data_any, mat_tmp
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (.NOT. PRESENT (last_column) .OR. scale_type .EQ. 'S', &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Column limiting supported only with scalar scaling",__LINE__,error)
    last_col = 0
    IF (PRESENT (last_column)) THEN
       CALL dbcsr_assert (last_column .GE. 0 .AND. last_column .LE.&
            dbcsr_nfullcols_total (matrix_a), dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Specified last column is out of bounds.",__LINE__,error)
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
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "wrong side="//side,__LINE__,error)
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
             CALL dbcsr_assert (.FALSE., dbcsr_fatal_level,&
                  dbcsr_wrong_args_error, routineN,&
                  "wrong logic",__LINE__,error)
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
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_scale_anytype

  !> \brief Interface for matrix scaling by a scalar
  SUBROUTINE dbcsr_scale_d(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(real_8), INTENT(IN)                 :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_d', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: scale_type
    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
    scale_type = 'S'
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_8) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (alpha_scalar))
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (REAL(alpha_scalar,real_4)))
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid data type",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_scale_d

  SUBROUTINE dbcsr_scale_s(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(real_4), INTENT(IN)                 :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_s', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: scale_type
    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
    scale_type = 'S'
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (alpha_scalar))
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid data type",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_scale_s
  
  SUBROUTINE dbcsr_scale_z(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    COMPLEX(real_8), INTENT(IN)              :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_z', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: scale_type
    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
    scale_type = 'S'
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_8) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (alpha_scalar))
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_4) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (CMPLX(alpha_scalar,KIND=real_4)))
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid data type",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_scale_z

  SUBROUTINE dbcsr_scale_c(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    COMPLEX(real_4), INTENT(IN)              :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_c', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: scale_type
    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
    scale_type = 'S'
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_4) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (alpha_scalar))
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid data type",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
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
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_scale_by_vector_anytype', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER                                  :: a_blk, a_col, a_col_size, &
                                                a_nze, a_row, a_row_size, &
                                                col_offset, error_handler, i, &
                                                row_offset
    LOGICAL                                  :: right, tr
    TYPE(dbcsr_data_obj)                     :: data_any
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

! check that alpha and matrix have the same data type
    CALL dbcsr_assert (dbcsr_get_data_type (matrix_a).EQ.alpha%d%data_type, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "wrong data type matrix_a="//&
         dbcsr_get_data_type (matrix_a)//' alpha='//alpha%d%data_type,__LINE__,error)
    !
    ! set vars
    right = .TRUE.
    SELECT CASE(side)
    CASE('right');right = .TRUE.
    CASE('left' );right = .FALSE.
    CASE DEFAULT
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, &
            routineN, "wrong side="//side,__LINE__,error)
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
    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_scale_by_vector_anytype


! *****************************************************************************
!> \brief Set a DBCSR matrix by alpha
!> \param[inout] matrix       DBCSR matrix
!> \param[in] alpha           a scalar
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_set_anytype(matrix, alpha, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: alpha
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_anytype', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, error_handler, row
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_block
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------
!

    CALL dbcsr_error_set(routineN, error_handler, error)
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
    CALL dbcsr_error_stop(error_handler, error)
    !
  END SUBROUTINE dbcsr_set_anytype

! *****************************************************************************
!> \brief Conjugate a DBCSR matrix
!> \param[inout] matrix       DBCSR matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_conjg(matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_conjg', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER                                  :: blk, col, error_handler, row
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_any
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------
!

    CALL dbcsr_error_set(routineN, error_handler, error)
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
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_conjg

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
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_anytype', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type_a, data_type_b
    INTEGER                                  :: blk, col, error_handler, row, &
                                                size_a, size_b
    LOGICAL                                  :: do_scale, tr
    TYPE(dbcsr_data_obj)                     :: data_block
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_scalar_type)                  :: my_alpha_scalar, &
                                                my_beta_scalar

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    data_type_a = dbcsr_get_data_type(matrix_a)
    data_type_b = dbcsr_get_data_type(matrix_b)
    !
    my_alpha_scalar = dbcsr_scalar_one (data_type_a)
    IF(PRESENT(alpha_scalar)) my_alpha_scalar = alpha_scalar
    my_beta_scalar = dbcsr_scalar_one (data_type_b)
    IF(PRESENT(beta_scalar)) my_beta_scalar = beta_scalar
    !
    ! let's go
    CALL dbcsr_assert (dbcsr_nblkrows_total(matrix_a).EQ.dbcsr_nblkrows_total(matrix_b), &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, "matrices not consistent",__LINE__,error)

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
            sizedata_guess=size_b-size_a,error=error)
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
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_add_anytype

  !> \brief Interface for dbcsr_add
  SUBROUTINE dbcsr_add_d(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    REAL(real_8), INTENT(IN)                 :: alpha_scalar, beta_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
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
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_add_d

  SUBROUTINE dbcsr_add_s(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    REAL(real_4), INTENT(IN)                 :: alpha_scalar, beta_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4 .AND. &
       dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar),&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_add_s

  SUBROUTINE dbcsr_add_z(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    COMPLEX(real_8), INTENT(IN)              :: alpha_scalar, beta_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
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
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_add_z

  SUBROUTINE dbcsr_add_c(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    COMPLEX(real_4), INTENT(IN)              :: alpha_scalar, beta_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_4 .AND. &
       dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_complex_4) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar),&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
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
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_hadamard_product', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER                                  :: blk, col, col_size, &
                                                error_handler, nze, row, &
                                                row_size
    LOGICAL                                  :: found, tr_a, tr_b
    TYPE(dbcsr_data_obj)                     :: a_data, b_data, c_data
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_get_data_type(matrix_a).EQ.dbcsr_get_data_type(matrix_b).AND.&
         dbcsr_get_data_type(matrix_a).EQ.dbcsr_get_data_type(matrix_c), &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, &
         "data types not consistent, need to fix that",__LINE__,error)

    CALL dbcsr_assert (dbcsr_nblkrows_total(matrix_a).EQ.dbcsr_nblkrows_total(matrix_b).AND.&
         dbcsr_nblkrows_total(matrix_c).EQ.dbcsr_nblkrows_total(matrix_a), &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, &
         "matrices not consistent",__LINE__,error)

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
          CALL dbcsr_assert (tr_a.EQV.tr_b, dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, &
               "tr not consistent, need to fix that",__LINE__,error)
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
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Only real double precision",__LINE__,error)
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_finalize (matrix_c, error=error)
    CALL dbcsr_data_clear_pointer (a_data)
    CALL dbcsr_data_clear_pointer (b_data)
    CALL dbcsr_data_release (c_data)
    CALL dbcsr_data_release (a_data)
    CALL dbcsr_data_release (b_data)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_hadamard_product

! *****************************************************************************
!> \brief Replaces blocks with blocks from another matrix.
!>
!> Every block in matrix B is copied over to matrix A, replacing any possibly
!> pre-existing blocks.
!> \param[in,out] matrix_a   replace blocks in this DBCSR matrix
!> \param[in] matrix_b       replace blocks taken from the DBCSR matrix
! *****************************************************************************
  SUBROUTINE dbcsr_replace_blocks(matrix_a, matrix_b, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_replace_blocks', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: data_c
    COMPLEX(KIND=sp), DIMENSION(:, :), &
      POINTER                                :: data_z
    INTEGER                                  :: blk, col, error_handler, row
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: data_r
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (dbcsr_get_data_type(matrix_a).EQ.dbcsr_get_data_type(matrix_b), &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, "data types not consistent",__LINE__,error)
    !
    ! let's go
    CALL dbcsr_assert (matrix_a%m%nblkrows_total.EQ.matrix_b%m%nblkrows_total, &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, "matrices not consistent",__LINE__,error)

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
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Unkown precision",__LINE__,error)
       END SELECT

    ENDDO

    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_finalize (matrix_a, error=error)
    !
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_replace_blocks

! *****************************************************************************
!> \brief add a constant to the diagonal of a matrix
!> \param[inout] matrix       DBCSR matrix
!> \param[in]    alpha_scalar scalar
! *****************************************************************************
  SUBROUTINE dbcsr_add_on_diag(matrix, alpha_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(KIND=real_8), INTENT(IN)            :: alpha_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_on_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, hold, mynode, &
                                                n, row
    INTEGER, DIMENSION(:), POINTER           :: row_blk_size
    LOGICAL                                  :: found, tr
    REAL(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: buff, data_d

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_nblkrows_total(matrix).EQ.dbcsr_nblkcols_total(matrix).AND.&
         dbcsr_nfullrows_total(matrix).EQ.dbcsr_nfullrows_total(matrix), &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, "matrices not consistent",__LINE__,error)

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
    CALL dbcsr_error_stop(error_handler, error)
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
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, error_handler, row
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
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
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Only real double precision",__LINE__,error)
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_finalize(diag, error=error)
    CALL dbcsr_error_stop(error_handler, error)
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
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, error_handler, i, &
                                                row, row_offset
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (dbcsr_nfullrows_total(matrix).LE.SIZE(diag), &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, "SIZE(diag) too small",__LINE__,error)

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
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Only real double precision",__LINE__,error)
       END SELECT

    ENDDO

    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_error_stop(error_handler, error)

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
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, error_handler, i, &
                                                row, row_offset
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (dbcsr_nfullrows_total(matrix).LE.SIZE(diag), &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, "SIZE(diag) too small",__LINE__,error)

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
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Only real double precision",__LINE__,error)
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_error_stop(error_handler, error)
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
  SUBROUTINE dbcsr_copy(matrix_b, matrix_a, name, keep_sparsity,&
       shallow_data, keep_imaginary, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity, shallow_data, &
                                                keep_imaginary
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: new_type
    INTEGER                                  :: error_handler
    LOGICAL                                  :: keep_sparse, shallow

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    shallow = .FALSE. ; IF (PRESENT (shallow_data)) shallow = shallow_data
    keep_sparse = .FALSE.
    IF (PRESENT (keep_sparsity)) keep_sparse = keep_sparsity

    !CALL dbcsr_assert (dbcsr_is_initialized (matrix_b),&
    !     dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
    !     "Target matrix must be initialized", error=error)
    CALL dbcsr_assert (.not.keep_sparse.or.dbcsr_valid_index(matrix_b),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Target matrix must be valid to keep its sparsity",__LINE__,error)
    CALL dbcsr_assert (.not.keep_sparse.or..not.shallow, dbcsr_warning_level,&
         dbcsr_wrong_args_error, routineN,&
         "Shallow copy not compatibly with sparsity retainment",__LINE__,error)
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
       CALL ensure_array_size(matrix_b%m%index, ub=SIZE(matrix_a%m%index), error=error)
       !
       ! copy index and data
       matrix_b%m%index(1:SIZE(matrix_a%m%index)) = matrix_a%m%index(:)
       IF (.NOT. shallow) THEN
          CALL dbcsr_assert (matrix_a%m%nze .LE. dbcsr_get_data_size(matrix_a),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Source matrix sizes not consistent!",__LINE__,error)
          CALL dbcsr_data_ensure_size (matrix_b%m%data_area,&
               dbcsr_get_data_size(matrix_a), error=error)
          IF (dbcsr_get_data_type(matrix_a) .EQ. dbcsr_get_data_type(matrix_b))&
               THEN
             CALL dbcsr_data_copyall (matrix_b%m%data_area,&
                  matrix_a%m%data_area)
          ELSE
             CALL dbcsr_data_convert (matrix_b%m%data_area,&
                  matrix_a%m%data_area, drop_real=keep_imaginary)
          ENDIF
       ELSE
          CALL dbcsr_assert (dbcsr_get_data_type(matrix_a) &
               .EQ. dbcsr_get_data_type(matrix_b), dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "Shallow copy only possible when retaining data type.", __LINE__, error)
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
    CALL dbcsr_error_stop(error_handler, error)
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
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_into_existing', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER                                  :: col_size, cs, dst_col, &
                                                dst_row, error_handler, rel, &
                                                row_size, rs, src_col, src_row
    LOGICAL                                  :: dst_tr, src_tr
    TYPE(dbcsr_data_obj)                     :: dst_data, src_data
    TYPE(dbcsr_iterator)                     :: dst_iter, src_iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_get_data_type(matrix_b)&
         .EQ. dbcsr_get_data_type(matrix_a), dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Matrices have different data types.",__LINE__,error)
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
          CALL dbcsr_assert (dbcsr_get_data_size (src_data)&
               .EQ. dbcsr_get_data_size (dst_data), dbcsr_fatal_level,&
               dbcsr_internal_error, routineN, "Block sizes not equal!",__LINE__,error)
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
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
               routineN, "Trouble syncing iterators",__LINE__,error)
       END SELECT
    END DO
    CALL dbcsr_iterator_stop (src_iter)
    CALL dbcsr_iterator_stop (dst_iter)
    CALL dbcsr_data_clear_pointer (src_data)
    CALL dbcsr_data_clear_pointer (dst_data)
    CALL dbcsr_data_release (src_data)
    CALL dbcsr_data_release (dst_data)
    CALL dbcsr_error_stop(error_handler, error)

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
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_columns', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER :: col_size, dst_col, dst_row, error_handler, first_col, &
      last_col, num_col, rel, row_size, src_col, src_row
    INTEGER, DIMENSION(:), POINTER           :: dst_col_offsets, &
                                                dst_col_sizes, &
                                                src_col_offsets, src_col_sizes
    LOGICAL                                  :: dst_tr, src_tr
    TYPE(dbcsr_data_obj)                     :: dst_data, src_data, tmp_buffer
    TYPE(dbcsr_iterator)                     :: dst_iter, src_iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_get_data_type(matrix_b)&
         .EQ. dbcsr_get_data_type(matrix_a), dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Matrices have different data types.",__LINE__,error)
    data_type = dbcsr_get_data_type (matrix_b)

    CALL dbcsr_assert (source_start .EQ. target_start, dbcsr_fatal_level,&
         dbcsr_unimplemented_error_nr, routineN, "Column shifting not supported",__LINE__,error)

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
          CALL dbcsr_assert (src_col_sizes(src_col) .EQ. dst_col_sizes(dst_col),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Unequal column sizes",__LINE__,error)
          !CALL dbcsr_unimplemented_error (routineN, "Partial block copy not tested",&
          !     error=error, error_level=dbcsr_warning_level)
          CALL dbcsr_assert (src_col_offsets(src_col).eq.dst_col_offsets(dst_col),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Block offsets must be equal",__LINE__,error)
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
             CALL dbcsr_assert (.FALSE., dbcsr_fatal_level,  dbcsr_internal_error, routineN,&
                  "Partial block addition not yet tested",__LINE__,error)
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
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
               routineN, "Trouble syncing iterators",__LINE__,error)
       END SELECT
    END DO
    CALL dbcsr_iterator_stop (src_iter)
    CALL dbcsr_iterator_stop (dst_iter)
    CALL dbcsr_data_clear_pointer (src_data)
    CALL dbcsr_data_clear_pointer (dst_data)
    CALL dbcsr_data_release (src_data)
    CALL dbcsr_data_release (dst_data)
    CALL dbcsr_error_stop(error_handler, error)
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
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_btriu', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_nze, col, col_size, &
                                                error_handler, row, row_size
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
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
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_btriu

! *****************************************************************************
!> \brief  triu of a dbcsr matrix
!> \param[inout] matrix_a  the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_triu(matrix_a, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_triu', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_nze, col, col_size, &
                                                error_handler, i, j, row, &
                                                row_size
    LOGICAL                                  :: tr
    REAL(dp), DIMENSION(:, :), POINTER       :: DATA
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
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
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_triu

! *****************************************************************************
!> \brief  symmetrize the block diagonal of a dbcsr matrix
!> \param[inout] matrix_a  the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_symmetrize_block_diag(matrix_a, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_symmetrize_block_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_size, &
                                                error_handler, i, j, row, &
                                                row_size
    LOGICAL                                  :: tr
    REAL(dp)                                 :: dum
    REAL(dp), DIMENSION(:, :), POINTER       :: DATA
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
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
    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_symmetrize_block_diag

! *****************************************************************************
!> \brief  tril of a dbcsr matrix
!> \param[inout] matrix_a  the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_tril(matrix_a, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_tril', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_nze, col, col_size, &
                                                error_handler, i, j, row, &
                                                row_size
    LOGICAL                                  :: tr
    REAL(dp), DIMENSION(:, :), POINTER       :: DATA
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
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
    CALL dbcsr_error_stop(error_handler, error)
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
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_filter', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_nze, col, col_size, &
                                                error_handler, my_method, &
                                                row, row_size
    LOGICAL                                  :: tr
    REAL(dp)                                 :: ddot, nrm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
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
                CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
                     routineN,"Only Frobenius based filtering",__LINE__,error)
             END SELECT
             !CASE (dbcsr_type_complex_4)
             !CASE (dbcsr_type_complex_8)
          CASE DEFAULT
             CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
                  routineN,"Only real double precision",__LINE__,error)
          END SELECT
       ENDDO
       CALL dbcsr_iterator_stop(iter)
       CALL dbcsr_finalize (matrix, error=error)

    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
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
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_norm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_offset, &
                                                error_handler, i, j, row, &
                                                row_offset
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

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
                CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
                     routineN,"Only real double precision",__LINE__,error)
             END SELECT
          ENDDO
          CALL dbcsr_iterator_stop(iter)
          CALL mp_sum(norm_vector,dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))
          norm_vector = SQRT(norm_vector)
       ENDIF

    CASE DEFAULT

       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, &
            "this norm is NYI",__LINE__,error)

    END SELECT
    CALL dbcsr_error_stop(error_handler, error)
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
                                                nc, nr, row, row_offset
    LOGICAL                                  :: any_sym, tr
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff_d
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: data_r
    TYPE(dbcsr_error_type)                   :: error
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    nr = dbcsr_nfullrows_total(matrix)
    nc = dbcsr_nfullcols_total(matrix)

    any_sym = dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_symmetric.OR.&
              dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_antisymmetric

    CALL dbcsr_assert (nr.EQ.nc, dbcsr_fatal_level, dbcsr_wrong_args_error, &
         routineN, "not a square matrix",__LINE__,error)

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
                  CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
                  routineN,"Only nonsymmetric matrix so far",__LINE__,error)
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
                  CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
                  routineN,"Only nonsymmetric matrix so far",__LINE__,error)
             !     buff_d(col_offset+j-1) = buff_d(col_offset+j-1) + ABS(data_z(i,j))
          ENDDO
          ENDDO
       CASE DEFAULT
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Wrong data type",__LINE__,error)
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL mp_sum(buff_d,dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))       
    dbcsr_gershgorin_norm = MAXVAL(buff_d)
    DEALLOCATE(buff_d)

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
    INTEGER                                  :: blk, col, row
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: data_r
    TYPE(dbcsr_error_type)                   :: error
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

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
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Wrong data type",__LINE__,error)
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)

    !dmp_max, this fixes a bug in g95
    ! -> ambigous interface for mp_max
    !    in module dbcsr_message_passing and wrongly
    !    exported from the module timings
    CALL dmp_max( dbcsr_maxabs_norm, dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))

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
    INTEGER                                  :: blk, col, row
    LOGICAL                                  :: any_sym, tr
    REAL(dp)                                 :: fac
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: data_r
    TYPE(dbcsr_error_type)                   :: error
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

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
               CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Only nonsymmetric matrix so far",__LINE__,error)
          dbcsr_frobenius_norm = dbcsr_frobenius_norm + fac * SUM(CONJG(data_c)*data_c)
       CASE (dbcsr_type_complex_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_z, tr, blk)
          fac = 1.0_dp
          IF(any_sym.AND.row.NE.col) &
               CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Only nonsymmetric matrix so far",__LINE__,error)
          dbcsr_frobenius_norm = dbcsr_frobenius_norm + fac * SUM(CONJG(data_z)*data_z)
       CASE DEFAULT
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Wrong data type",__LINE__,error)
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL mp_sum(dbcsr_frobenius_norm,dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))
    dbcsr_frobenius_norm = SQRT(dbcsr_frobenius_norm)

  END FUNCTION dbcsr_frobenius_norm


! *****************************************************************************
!> \brief Sums blocks in a replicated dbcsr matrix.
!> \param[in,out] matrix      dbcsr matrix to operate on
! *****************************************************************************
  SUBROUTINE dbcsr_sum_replicated (matrix, error)
    TYPE(dbcsr_obj), INTENT(inout)           :: matrix
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sum_replicated', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: comm, error_handler, &
                                                max_data_size, &
                                                max_index_size, mynode, &
                                                numnodes, p
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: send_index
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: all_sizes
    INTEGER, DIMENSION(2)                    :: my_sizes
    INTEGER, DIMENSION(:), POINTER           :: recv_index
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: send_dp
    REAL(kind=dp), DIMENSION(:), POINTER     :: recv_dp
    TYPE(dbcsr_mp_obj)                       :: mp
    TYPE(dbcsr_obj)                          :: recv_matrix

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
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
    CALL dbcsr_create(recv_matrix, "Temporary receiving buffer",&
         dbcsr_distribution(matrix), dbcsr_get_matrix_type(matrix),&
         dbcsr_row_block_sizes(matrix), dbcsr_col_block_sizes(matrix),&
         0, 0, dbcsr_get_data_type(matrix), error=error)
    CALL ensure_array_size (recv_matrix%m%index, ub=max_index_size, error=error)
    CALL dbcsr_data_ensure_size (recv_matrix%m%data_area, max_data_size, error=error)
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
    CALL dbcsr_destroy (recv_matrix, error=error,force=.TRUE.)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_sum_replicated


  SUBROUTINE dbcsr_trace_a_b_d(matrix_a, matrix_b, trace, trans_a, trans_b, local_sum, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a, matrix_b
    REAL(kind=real_8), INTENT(INOUT)         :: trace
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: trans_a, trans_b
    LOGICAL, INTENT(IN), OPTIONAL            :: local_sum
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_a_b_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    REAL(kind=real_4)                        :: trace_4

    CALL dbcsr_error_set(routineN, error_handler, error)
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
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_trace_a_b_d




! *****************************************************************************
!> \brief Calls the functions that are called in the included files but not
!>        in thin one. Needed to make prettify think the functions ary really
!>        used.
!> \param[in] matrix          matrix
!> \param[out] temp1, temp2   pointers to the data
! *****************************************************************************
  SUBROUTINE prettify_fooler (matrix, temp1, temp2, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    REAL(KIND=real_4), DIMENSION(:), POINTER :: temp1
    REAL(KIND=real_8), DIMENSION(:), POINTER :: temp2
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    INTEGER                                  :: i
    INTEGER(KIND=int_8)                      :: key
    INTEGER, DIMENSION(:), POINTER           :: ip
    LOGICAL                                  :: tr
    TYPE(btree_2d_data_c)                    :: bt_dc
    TYPE(btree_2d_data_d)                    :: bt_dd
    TYPE(btree_2d_data_s)                    :: bt_ds
    TYPE(btree_2d_data_z)                    :: bt_dz
    TYPE(btree_c)                            :: bt_c
    TYPE(btree_d)                            :: bt_d
    TYPE(btree_s)                            :: bt_s
    TYPE(btree_z)                            :: bt_z
    TYPE(dbcsr_data_obj)                     :: md
    TYPE(dbcsr_distribution_obj)             :: dist

    IF (.FALSE.) THEN
       temp1 => dbcsr_get_data_p (matrix%data_area, 0.0_real_4)
       temp2 => dbcsr_get_data_p (matrix%data_area, 0.0_real_8)
       ip => array_data (dbcsr_distribution_row_dist (matrix%dist))
       ip => array_data (dbcsr_distribution_col_dist (matrix%dist))
       CALL mp_sum (i, 0)
       CALL mp_bcast (temp1, 0, 0)
       CALL addto_array (matrix%data_area, temp1, error=error)
       CALL btree_add_s (bt_s, key, bt_ds)
       CALL btree_add_d (bt_d, key, bt_dd)
       CALL btree_add_c (bt_c, key, bt_dc)
       CALL btree_add_z (bt_z, key, bt_dz)
       CALL btree_get_s (bt_s, key, bt_ds)
       CALL btree_get_d (bt_d, key, bt_dd)
       CALL btree_get_c (bt_c, key, bt_dc)
       CALL btree_get_z (bt_z, key, bt_dz)
       WRITE(*,*)dbcsr_repl_none, dbcsr_repl_full
       WRITE(*,*) checker_tr (1, 1)
       WRITE(*,*)checker_square_proc(1,1,RESHAPE((/ 1 /),(/1,1/)),(/0/),(/0/))
       dist = matrix%dist
       CALL dbcsr_get_data (md, temp2)
       CALL dbcsr_find_column (i, ip(1), ip(2), (/ 1 /), (/ 1 /), ip(3), tr)
    ENDIF
  END SUBROUTINE prettify_fooler

#include "dbcsr_operations_d_.F"
#include "dbcsr_operations_z_.F"
#include "dbcsr_operations_s_.F"
#include "dbcsr_operations_c_.F"

END MODULE dbcsr_operations

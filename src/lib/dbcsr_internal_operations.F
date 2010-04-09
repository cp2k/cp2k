!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Library-internal subroutines for DBCSR matrix operations.
!> \author  Urban Borstnik
!> \date    2010-02-23
!> \version 0.9
!>
!> <b>Modification history:</b>
!  - 2010-02-23 Moved from dbcsr_operations
! *****************************************************************************
MODULE dbcsr_internal_operations
  USE array_types,                     ONLY: array_data
  USE dbcsr_block_operations,          ONLY: dbcsr_data_clear
  USE dbcsr_dist_operations,           ONLY: checker_tr
  USE dbcsr_error_handling
  USE dbcsr_index_operations,          ONLY: dbcsr_repoint_index,&
                                             merge_index_arrays
  USE dbcsr_io,                        ONLY: all_print,&
                                             dbcsr_print
  USE dbcsr_kinds,                     ONLY: dp,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8,&
                                             sp
  USE dbcsr_machine,                   ONLY: default_output_unit,&
                                             m_walltime
  USE dbcsr_message_passing,           ONLY: mp_allgather,&
                                             mp_irecv,&
                                             mp_isend,&
                                             mp_type_descriptor_type,&
                                             mp_type_free,&
                                             mp_type_make,&
                                             mp_waitall
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_sizes, dbcsr_data_clear_pointer, &
       dbcsr_data_ensure_size, dbcsr_data_init, dbcsr_data_new, &
       dbcsr_data_release, dbcsr_data_set_pointer, dbcsr_destroy_array, &
       dbcsr_distribution, dbcsr_distribution_local_rows, &
       dbcsr_distribution_mp, dbcsr_get_data_size, dbcsr_get_data_type, &
       dbcsr_get_matrix_type, dbcsr_init, dbcsr_mp_grid_setup, &
       dbcsr_mp_group, dbcsr_mp_has_subgroups, dbcsr_mp_my_col_group, &
       dbcsr_mp_my_row_group, dbcsr_mp_mynode, dbcsr_mp_mypcol, &
       dbcsr_mp_myprow, dbcsr_mp_npcols, dbcsr_mp_nprows, dbcsr_mp_numnodes, &
       dbcsr_mp_pgrid, dbcsr_nfullcols_local, dbcsr_nfullrows_local, &
       dbcsr_row_block_sizes, dbcsr_scalar_are_equal, dbcsr_scalar_negative, &
       dbcsr_scalar_one, dbcsr_scalar_zero

!$ USE dbcsr_methods, ONLY: dbcsr_distribution_thread_dist
  USE dbcsr_mp_operations,             ONLY: dbcsr_irecv_any,&
                                             dbcsr_isend_any,&
                                             dbcsr_mp_type_from_anytype
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size,&
                                             pointer_replace
  USE dbcsr_types,                     ONLY: &
       dbcsr_2d_array_type, dbcsr_data_obj, dbcsr_mp_obj, dbcsr_obj, &
       dbcsr_scalar_type, dbcsr_type, dbcsr_type_complex_4, &
       dbcsr_type_complex_8, dbcsr_type_real_4, dbcsr_type_real_8
  USE dbcsr_util,                      ONLY: dbcsr_checksum,&
                                             dbcsr_set_debug
  USE dbcsr_work_operations,           ONLY: dbcsr_create,&
                                             dbcsr_finalize,&
                                             dbcsr_work_create

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_internal_operations'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034


  PUBLIC :: dbcsr_mult_m_e_e


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
    TYPE(dbcsr_obj), INTENT(INOUT)           :: product_matrix
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
      '(A,1X,I4,"->",I4,2(1X,"(",I3,"x",I3,")"),1X,"IM (",I3,"x",I3,")")'
    INTEGER, PARAMETER                       :: idata = 1, imeta = 2, &
                                                M_L = 2, M_P = 1, M_R = 3, &
                                                RC_C = 2, RC_R = 1
    LOGICAL, PARAMETER                       :: use_combined_types = .FALSE.

    CHARACTER                                :: data_type
    INTEGER :: error_handler, grp, i, left_col_image, left_col_mult, &
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
    INTEGER :: right_src_prow, right_src_v_col, right_src_v_row, v_k, v_ki
    INTEGER(KIND=int_8)                      :: flop_single, flop_total
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
    t_all = 0.0_dp
    t_dgemm = 0.0_dp
!$omp parallel default (none) &
!$omp shared (nthreads, product_matrix, error)
!$omp single
    nthreads = 1
    !$ nthreads = OMP_GET_NUM_THREADS ()
!$omp end single
    CALL dbcsr_work_create(product_matrix,&
         work_mutable=.FALSE., n=nthreads, error=error)
!$omp end parallel
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
                left_src_v_col = MODULO(left_src_v_col, nvirt_k)
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
                right_src_v_row = MODULO(right_src_v_row, nvirt_k)
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
             right_send_v_row = MODULO (right_send_v_row, nvirt_k)
             right_send_prow = right_send_v_row / right_row_nimages
             right_send_irow = MODULO (right_send_v_row, right_row_nimages)
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
             right_dst_v_row = MODULO (right_dst_v_row, nvirt_k)
             right_dst_prow = right_dst_v_row / right_row_nimages
             right_dst_irow = 1 + MODULO (right_dst_v_row, right_row_nimages)
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
             right_src_v_row = MODULO (right_src_v_row, nvirt_k)
             right_src_prow = right_src_v_row / right_row_nimages
             right_src_irow = 1 + MODULO (right_src_v_row, right_row_nimages)
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
                CALL dbcsr_isend_any (right_data_sp, right_send_p,&
                     grp, right_data_sr(v_ki+1, 1), tag=right_dst_v_row,&
                     error=error)
                CALL dbcsr_irecv_any (right_data_rp, right_recv_p,&
                     grp, right_data_rr(v_ki+1, 1), tag=right_src_v_row,&
                     error=error)
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
             left_send_v_col = MODULO (left_send_v_col, nvirt_k)
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
             left_dst_v_col = MODULO (left_dst_v_col, nvirt_k)
             left_dst_pcol = left_dst_v_col / left_col_nimages
             left_dst_icol = 1 + MODULO (left_dst_v_col, left_col_nimages)
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
             left_src_v_col = MODULO (left_src_v_col, nvirt_k)
             left_src_pcol = left_src_v_col / left_col_nimages
             left_src_icol = 1 + MODULO (left_src_v_col, left_col_nimages)
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
                CALL dbcsr_isend_any (left_data_sp, left_send_p,&
                     grp, left_data_sr(1, v_ki+1), tag=left_dst_v_col,&
                     error=error)
                CALL dbcsr_irecv_any (left_data_rp, left_recv_p,&
                     grp, left_data_rr(1, v_ki+1), tag=left_src_v_col,&
                     error=error)
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
             flop_single = 0
!$omp parallel default (none) &
!$omp shared (left_buffer_calc, right_buffer_calc, &
!$omp         v_ki, &
!$omp         product_matrix, submatrix_selection, &
!$omp         retain_sparsity, alpha, error) &
!$omp reduction (+: flop_single, t_all, t_dgemm)
             CALL dbcsr_nn_mult_lin(&
                  left_buffer_calc%mats(1, v_ki+1)%m,&
                  right_buffer_calc%mats(v_ki+1, 1)%m,&
                  product_matrix%m, submatrix_box=submatrix_selection,&
                  retain_sparsity=retain_sparsity,&
                  alpha=alpha,&
                  flop=flop_single,&
                  t_all=t_all, t_dgemm=t_dgemm, error=error)
!$omp end parallel
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
!$omp parallel default(none) &
!$omp shared (product_matrix, error)
    CALL dbcsr_finalize(product_matrix, error=error)
!$omp end parallel
    !
    IF (dbg) THEN
       checksum = dbcsr_checksum (product_matrix, error=error)
       IF ((output_unit>0)) THEN
          WRITE(output_unit,'(1X,A,1X,F9.4)')"Product Checksum=",checksum
       ENDIF
    ENDIF
    !
    IF (PRESENT (flop)) flop = flop_total
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
          buffer_set%mats(row_image, col_image)%m%transpose&
               = source_matrix%m%transpose
          buffer_set%mats(row_image, col_image)%m%negate_real&
               = source_matrix%m%negate_real
          buffer_set%mats(row_image, col_image)%m%negate_imaginary&
               = source_matrix%m%negate_imaginary
       ENDDO
    ENDDO
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE setup_buffer_matrices



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
       retain_sparsity, alpha,&
       t_all, t_dgemm, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    INTEGER, DIMENSION(6), INTENT(in)        :: submatrix_box
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: alpha
    REAL(KIND=dp), INTENT(INOUT)             :: t_all, t_dgemm
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_nn_mult_lin', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER :: a_blk, a_col_logical, a_col_size_logical, a_row_l, &
      a_row_logical, a_row_size_logical, b_blk, b_col_logical, &
      b_col_size_logical, b_row_logical, b_row_size_logical, blk_end, &
      blk_size, blk_start, c_blk_pt, c_col_logical, c_col_size_logical, &
      c_nze, c_row_logical, c_row_size_logical, clamp_k_len, datasize, &
      error_handler, first_k, first_mergeable, i, ithread, j, last_k, &
      lastblk, new_blk, new_row_blks, nthreads, old_blk, old_row_begin_blk, &
      old_size, row_begin_blk, size_c_blk_ps, stack_p, stack_p_a, stack_p_b, &
      stack_p_c, stat, submat_f_col, submat_f_row, submat_l_col, &
      submat_l_row, this_right_limit
    INTEGER(KIND=int_4)                      :: offset
    INTEGER(KIND=int_8)                      :: lflop, one_flop
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: c_blk_ps, c_blk_ps_p
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: new_index
    INTEGER, DIMENSION(:), POINTER :: left_col_blk_offset, left_local_rows, &
      left_row_blk_offset, new_blk_p, new_col_i, new_row_i, &
      product_col_blk_offset, right_col_blk_offset

!$  INTEGER, DIMENSION(:), POINTER           :: product_thread_dist
    LOGICAL :: block_exists, flush_stack, keep_sparsity, &
      left_limit_rows, limit_k, negate_alpha, product_is_symmetric, &
      right_limit_cols, tr_p
    REAL(KIND=dp)                            :: epoch, mt_t_dgemm
    REAL(KIND=real_8)                        :: one_time
    TYPE(block_parameters), &
      DIMENSION(dgemm_stack_size)            :: parameters_a, parameters_b, &
                                                parameters_c
    TYPE(dbcsr_scalar_type)                  :: beta_add, beta_new, my_alpha, &
                                                which_beta
    TYPE(dgemm_join), &
      DIMENSION(dgemm_stack_size)            :: dgemm_joins

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    data_type = product%data_type
    my_alpha = dbcsr_scalar_one (data_type)
    IF (PRESENT (alpha)) my_alpha = alpha
    beta_new = dbcsr_scalar_zero (data_type)
    beta_add = dbcsr_scalar_one (data_type)
    keep_sparsity = .FALSE.
    IF (PRESENT (retain_sparsity)) THEN
       !WRITE(*,*)routineN//" Retaining sparsity!"
       keep_sparsity = .TRUE.
    END IF
    product_is_symmetric = product%symmetry
    ithread = 0 ; nthreads = 1
    !$ ithread = OMP_GET_THREAD_NUM () ; nthreads = OMP_GET_NUM_THREADS ()
    epoch = m_walltime()
    DBG "Doing NN matrix multiplication",left%name,' x ',right%name
    DBG "Pre-mult blk sizes, data", product%wms(ithread+1)%lastblk,&
         product%wms(ithread+1)%datasize
    left_row_blk_offset => array_data (left%row_blk_offset)
    left_col_blk_offset => array_data (left%col_blk_offset)
    right_col_blk_offset => array_data (right%col_blk_offset)
    product_col_blk_offset => array_data (product%col_blk_offset)
!$  NULLIFY (product_thread_dist)    
!$  product_thread_dist => array_data (&
!$       dbcsr_distribution_thread_dist (product%dist))
    left_local_rows => array_data (dbcsr_distribution_local_rows (left%dist))
    ! Submatrix handling
    CALL dbcsr_assert (submatrix_box(1) .EQ. 0,&
         dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
         "Limiting first row not supported",__LINE__,error)
!vw    CALL dbcsr_assert (submatrix_box(3) .EQ. 0,&
!vw         dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
!vw         "Limiting first column not supported",__LINE__,error)
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
    DBG 'A rows:',left%nblkrows_total
    old_size = product%wms(ithread+1)%lastblk
    mt_t_dgemm = 0.0_dp
    lflop = 0
    one_time = 0.0_dp
    IF (left%nblks .GT. 0 .AND. right%nblks .GT. 0) THEN
       !
    ALLOCATE (new_index (3, left%nblkrows_local*right%nblkcols_local))
    ! c_blk_ps contains the block pointers in the product's current row
    ALLOCATE (c_blk_ps(product%nblkcols_total), stat=stat)
    CALL dbcsr_assert (stat == 0, dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "c_blk_ps",__LINE__,error)
    c_blk_ps(:) = 0
    ALLOCATE (c_blk_ps_p(product%nblkcols_total), stat=stat)
    CALL dbcsr_assert (stat == 0, dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "c_blk_ps_p",__LINE__,error)
    size_c_blk_ps = 0
    c_blk_pt = product%wms(ithread+1)%datasize + 1
    lastblk = 0
    datasize = 0
    !
    old_blk = 1
    first_mergeable = 1
    stack_p = 0 ; stack_p_a = 0 ; stack_p_b = 0 ; stack_p_c = 0
    new_blk = 0
    DO a_row_l = 1, left%nblkrows_local
       a_row_logical = left_local_rows (a_row_l)
!$     IF (product_thread_dist(a_row_logical) .ne. ithread) CYCLE
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

                IF (blk_end .LT. submat_f_col) THEN
                   CYCLE
                ELSEIF (blk_end .GE. submat_f_col.and. blk_start.LT. submat_f_col) THEN
                   WRITE(*,*) routinen//' blk_end',blk_end,' submat_f_col', submat_f_col,&
                        ' blk_start',blk_start,' submat_f_col',submat_f_col,' b_col_logical',b_col_logical
                   ! some more work here and further with none maching col/block
                   CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr,&
                        routineN, "none maching first column/block NYI", __LINE__, error)
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
             negate_alpha = parameters_a(stack_p_a)%tr .AND. left%negate_real
             negate_alpha = negate_alpha .NEQV. (parameters_b(stack_p_b)%tr .AND. right%negate_real)
             !write(*,*)routineN//" tr_p:", tr_p
             !write(*,*)routineN//" p.s.n", product%sym_negation
             negate_alpha = negate_alpha .NEQV. (tr_p .AND. product%negate_real)
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
    CALL merge_index_arrays (new_row_i, new_col_i, new_blk_p, i,&
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
    !call dbcsr_assert (last, "LE", size(params), dbcsr_fatal_level,&
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
       !call dbcsr_assert (sp_b, "LE", stack_size_b, dbcsr_fatal_level,&
       !     cp_internal_error, routineN, "spb too large")
       !call dbcsr_assert (sp_c, "LE", stack_size_c, dbcsr_fatal_level,&
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
! The following routines are helped here to help the compiler optimize them
! out.
! *****************************************************************************

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



#include "dbcsr_internal_operations_d.F"
#include "dbcsr_internal_operations_z.F"
#include "dbcsr_internal_operations_s.F"
#include "dbcsr_internal_operations_c.F"


END MODULE dbcsr_internal_operations

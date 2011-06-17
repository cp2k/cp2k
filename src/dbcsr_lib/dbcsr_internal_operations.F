!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
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
  USE array_types,                     ONLY: array_data,&
                                             array_equality,&
                                             array_exists,&
                                             array_hold
  USE dbcsr_block_access,              ONLY: dbcsr_put_block,&
                                             dbcsr_reserve_blocks
  USE dbcsr_block_operations,          ONLY: dbcsr_data_clear
  USE dbcsr_config,                    ONLY: &
       detailed_timing, kernel_timing, mm_driver, mm_driver_blas, &
       mm_driver_cuda, mm_driver_matmul, mm_driver_plasma, mm_driver_smm, &
       mm_stack_size, use_CUDA_host_pinned_memory, use_MPI_memory, &
       use_combined_types
  USE dbcsr_cuda_memory
  USE dbcsr_cuda_methods
  USE dbcsr_cuda_operations
  USE dbcsr_cuda_types
  USE dbcsr_data_methods,              ONLY: &
       dbcsr_data_clear_pointer, dbcsr_data_ensure_size, dbcsr_data_get_size, &
       dbcsr_data_get_size_referenced, dbcsr_data_hold, dbcsr_data_init, &
       dbcsr_data_new, dbcsr_data_release, dbcsr_data_set_pointer, &
       dbcsr_data_set_size_referenced, dbcsr_get_data_p_c, &
       dbcsr_get_data_p_d, dbcsr_get_data_p_s, dbcsr_get_data_p_z
  USE dbcsr_dist_operations,           ONLY: dbcsr_reset_vlocals,&
                                             image_calculator
  USE dbcsr_error_handling
  USE dbcsr_index_operations,          ONLY: dbcsr_count_row_index,&
                                             dbcsr_expand_row_index,&
                                             dbcsr_has_local_row_index,&
                                             dbcsr_repoint_index
  USE dbcsr_io,                        ONLY: dbcsr_print,&
                                             print_xfer_timings
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_kinds,                     ONLY: &
       dp, int_4, int_4_size, int_8, real_4, real_4_size, real_8, &
       real_8_size, sp
  USE dbcsr_machine,                   ONLY: default_output_unit,&
                                             m_flush,&
                                             m_walltime
  USE dbcsr_message_passing,           ONLY: mp_allgather,&
                                             mp_irecv,&
                                             mp_isend,&
                                             mp_sum,&
                                             mp_type_descriptor_type,&
                                             mp_type_free,&
                                             mp_type_make,&
                                             mp_waitall
  USE dbcsr_methods,                   ONLY: &
       dbcsr_destroy_array, dbcsr_distribution_col_dist, &
       dbcsr_distribution_has_threads, dbcsr_distribution_local_cols, &
       dbcsr_distribution_local_rows, dbcsr_distribution_mp, &
       dbcsr_distribution_row_dist, dbcsr_distribution_thread_dist, &
       dbcsr_get_data_type, dbcsr_get_index_memory_type, &
       dbcsr_get_num_blocks, dbcsr_image_dist_init, dbcsr_init, &
       dbcsr_mp_grid_setup, dbcsr_mp_group, dbcsr_mp_has_subgroups, &
       dbcsr_mp_my_col_group, dbcsr_mp_my_row_group, dbcsr_mp_mynode, &
       dbcsr_mp_mypcol, dbcsr_mp_myprow, dbcsr_mp_npcols, dbcsr_mp_nprows, &
       dbcsr_mp_numnodes, dbcsr_mp_pgrid, dbcsr_nblkcols_local, &
       dbcsr_nblkrows_local, dbcsr_nblkrows_total, dbcsr_nfullcols_local, &
       dbcsr_nfullrows_local, dbcsr_valid_index
  USE dbcsr_mp_operations,             ONLY: dbcsr_irecv_any,&
                                             dbcsr_isend_any,&
                                             dbcsr_mp_type_from_anytype
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size
  USE dbcsr_types,                     ONLY: &
       dbcsr_2d_array_type, dbcsr_data_obj, dbcsr_iterator, &
       dbcsr_memory_CUDA_host_pinned, dbcsr_memory_MPI, dbcsr_memory_default, &
       dbcsr_mp_obj, dbcsr_obj, dbcsr_slot_size, dbcsr_type, &
       dbcsr_type_complex_4, dbcsr_type_complex_8, dbcsr_type_real_4, &
       dbcsr_type_real_8, dbcsr_work_type
  USE dbcsr_util,                      ONLY: count_bins,&
                                             dbcsr_checksum
  USE dbcsr_work_operations,           ONLY: dbcsr_create,&
                                             dbcsr_finalize

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_internal_operations'

  CHARACTER(len=*), PARAMETER, PRIVATE :: int_print = "(10(1X,I7))"

  REAL, PARAMETER                      :: default_resize_factor = 1.618034

  INTEGER :: multrec_calls = 0


#if defined(__HAS_NO_OMP_3)
#define __COLLAPSE2
#else
#define __COLLAPSE2 collapse(2)
#endif



  PUBLIC :: dbcsr_mult_m_e_e
  PUBLIC :: dbcsr_insert_blocks


  LOGICAL, PARAMETER :: debug_mod  = .FALSE.
  LOGICAL, PARAMETER :: careful_mod = .FALSE.


#define temp_transpose(v, r, c) RESHAPE(TRANSPOSE(RESHAPE(v,(/r,c/))),(/r*c/))

  INTEGER, PARAMETER, PRIVATE :: rpslot_owner = 1
  INTEGER, PARAMETER, PRIVATE :: rpslot_addblks = 2
  INTEGER, PARAMETER, PRIVATE :: rpslot_addoffset = 3
  INTEGER, PARAMETER, PRIVATE :: rpslot_oldblks = 4
  INTEGER, PARAMETER, PRIVATE :: rpslot_oldoffset = 5
  INTEGER, PARAMETER, PRIVATE :: rpslot_totaloffset = 6
  INTEGER, PARAMETER, PRIVATE :: rpnslots = 6

  INTEGER, PARAMETER, PRIVATE :: n_mult_params = 7
  INTEGER, PARAMETER, PRIVATE :: p_m = 1
  INTEGER, PARAMETER, PRIVATE :: p_n = 2
  INTEGER, PARAMETER, PRIVATE :: p_k = 3
  INTEGER, PARAMETER, PRIVATE :: p_a_first = 4
  INTEGER, PARAMETER, PRIVATE :: p_b_first = 5
  INTEGER, PARAMETER, PRIVATE :: p_c_first = 6
  INTEGER, PARAMETER, PRIVATE :: p_c_blk = 7


  LOGICAL, PARAMETER, PRIVATE :: verbose_acc = .FALSE.

  ! Types needed for the hashtable.
  TYPE ele_type
     INTEGER :: c=0
     INTEGER :: p=0
  END TYPE ele_type
  TYPE hash_table_type
     TYPE(ele_type), DIMENSION(:), POINTER :: table
     INTEGER :: nele=0
     INTEGER :: nmax=0
     INTEGER :: prime=0
  END TYPE hash_table_type

!> \var mult_params    array of ((m, n, k, a_first, b_first, c_first, c_blk),*)
  TYPE carrier_type
     INTEGER :: id
     LOGICAL :: c_has_symmetry, keep_sparsity, use_eps
     LOGICAL :: local_indexing
     TYPE(hash_table_type), DIMENSION(:), POINTER :: c_hashes
     INTEGER, DIMENSION(:,:), POINTER              :: mult_params
     INTEGER :: mm_stack_p
     INTEGER :: zero_first, zero_last
     INTEGER, DIMENSION(:), POINTER :: m_sizes, n_sizes, k_sizes
     INTEGER, DIMENSION(:), POINTER :: c_local_rows, c_local_cols, k_locals,&
                                       c_global_rows, c_global_cols
     REAL(KIND=sp), DIMENSION(:), POINTER :: row_max_epss, a_norms, b_norms
     REAL(KIND=real_8)     :: eps
     TYPE(dbcsr_work_type) :: my_wm
     INTEGER               :: lastblk, datasize
     INTEGER               :: original_lastblk
     INTEGER(kind=int_8)   :: flop
     REAL(KIND=dp)         :: t_index, t_gemm
     TYPE(dbcsr_data_obj)  :: a_data, b_data, c_data
  END TYPE carrier_type



  TYPE(dbcsr_cuda_mem_type) :: c_dev
  TYPE(dbcsr_cuda_mem_type), DIMENSION(2) :: a_dev, b_dev
  TYPE(dbcsr_cuda_mem_type_i4), DIMENSION(2) :: params_dev
  TYPE(dbcsr_cuda_mem_type_i4)               :: c_locks_dev
  INTEGER :: data_dev_calc, param_dev_calc, data_dev_xfer, param_dev_xfer,&
       param_host_fill, param_host_xfer

  REAL(kind=dp), DIMENSION(:), POINTER :: a_pin, b_pin, c_pin
  INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_p
  REAL(kind=dp) :: t_xfer_c_in, t_xfer_c_out, t_xfer_a, t_xfer_b,&
       t_calc, t_xfer_params, t_devdata_resize, t_calc_step, t_process_stack
  INTEGER :: len_xfer_a, len_xfer_b, len_xfer_c_in, len_xfer_c_out
  REAL :: t_r
  REAL(kind=dp), PARAMETER :: factor_mega = 1E6_dp, factor_giga = 1E9_dp


CONTAINS


! *****************************************************************************
!> \brief Multiplies two DBCSR matrices
!>
!> \param[in] left_set             set of imaged left matrices
!> \param[in] right_set            set of imaged right matrices
!> \param[out] product             DBCSR product matrix
!> \param[in,out] error            cp2k error
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix; default is no
!> \param[out] flop                (optional) effective flop
! *****************************************************************************
  SUBROUTINE dbcsr_mult_m_e_e (left_set, right_set, product_matrix,&
       error, retain_sparsity, &
       filter_eps, flop)
    TYPE(dbcsr_2d_array_type), POINTER       :: left_set, right_set
    TYPE(dbcsr_obj), INTENT(INOUT)           :: product_matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(kind=real_8), INTENT(in), OPTIONAL  :: filter_eps
    INTEGER(KIND=int_8), INTENT(OUT)         :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mult_m_e_e', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=80), PARAMETER :: &
      fdata = '(A,1X,I4,"(",2(I3),"x",2(I3),")","(",I3,"x",I3,")")', fxfer = &
      '(A,1X,I4,"->",I4,2(1X,"(",I3,"x",I3,")"),1X,"IM (",I3,"x",I3,")")'
    INTEGER, PARAMETER :: id_bytes = 3, id_recv = 2, id_send = 1, &
      id_time = 1, id_waittime = 2, idata = 1, ileft = 0, imeta = 2, &
      iright = 2, M_L = 2, M_P = 1, M_R = 3, RC_C = 2, RC_R = 1
    LOGICAL, PARAMETER                       :: excessive_output = .FALSE., &
                                                time_xfers = .FALSE.

    INTEGER :: data_type, error_handler, error_handler2, flop_metronome, grp, &
      i, ithread, left_col_image, left_col_mult, left_col_nimages, &
      left_data_recv_size, left_data_send_size, left_dst_icol, left_dst_irow, &
      left_dst_p, left_dst_pcol, left_dst_prow, left_dst_vcol, left_dst_vrow, &
      left_index_recv_size, left_index_send_size, left_max_nblks, &
      left_max_nze, left_myfirstvcol, left_myfirstvrow, left_mypcol, &
      left_myprow, left_npcols, left_nprows, left_recv_icol, left_recv_irow, &
      left_recv_p, left_recv_pcol, left_recv_prow, left_recv_vcol, &
      left_recv_vrow, left_row_image, left_row_mult, left_row_nimages, &
      left_send_icol
    INTEGER :: left_send_irow, left_send_p, left_send_pcol, left_send_prow, &
      left_send_vcol, left_send_vrow, left_src_icol, left_src_irow, &
      left_src_p, left_src_pcol, left_src_prow, left_src_vcol, left_src_vrow, &
      metronome, min_nimages, mp_group, mynode, nblkrows_total, &
      nblkrows_used, nsteps_k, nthreads, numnodes, nvirt_k, output_unit, &
      right_col_image, right_col_mult, right_col_nimages, &
      right_data_recv_size, right_data_send_size, right_dst_icol, &
      right_dst_irow, right_dst_p, right_dst_pcol, right_dst_prow, &
      right_dst_vcol, right_dst_vrow, right_index_recv_size, &
      right_index_send_size, right_max_nblks, right_max_nze
    INTEGER :: right_myfirstvcol, right_myfirstvrow, right_mypcol, &
      right_myprow, right_npcols, right_nprows, right_recv_icol, &
      right_recv_irow, right_recv_p, right_recv_pcol, right_recv_prow, &
      right_recv_vcol, right_recv_vrow, right_row_image, right_row_mult, &
      right_row_nimages, right_send_icol, right_send_irow, right_send_p, &
      right_send_pcol, right_send_prow, right_send_vcol, right_send_vrow, &
      right_src_icol, right_src_irow, right_src_p, right_src_pcol, &
      right_src_prow, right_src_vcol, right_src_vrow, row, size_guess, stat, &
      v_k, v_ki
    INTEGER(KIND=int_8)                      :: flop_single, flop_total, size8
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: row_counts, total_row_counts
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: left_data_rr, left_data_sr, &
      left_index_rr, left_index_sr, right_data_rr, right_data_sr, &
      right_index_rr, right_index_sr
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: left_sizes, my_sizes, &
                                                right_sizes
    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: all_sizes
    INTEGER, DIMENSION(3, 2)                 :: mp_rc_groups
    INTEGER, DIMENSION(:), POINTER           :: left_index_rp, left_index_sp, &
                                                local_rows, right_index_rp, &
                                                right_index_sp
    INTEGER, DIMENSION(:, :), POINTER        :: left_pgrid, product_pgrid, &
                                                right_pgrid
    INTEGER, SAVE                            :: mult_id = 0
    LOGICAL                                  :: keep_sparsity, list_indexing, &
                                                my_use_plasma, otf_filtering
    REAL(KIND=dp)                            :: checksum, t_all, t_gemm, &
                                                trun, trun_t, tstart, tstop
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: xfer_timings
    REAL(KIND=real_8)                        :: fill_guess, left_fill, &
                                                right_fill
    REAL(kind=sp), ALLOCATABLE, DIMENSION(:) :: left_norms, right_norms, &
                                                row_max_epss
    TYPE(carrier_type), SAVE                 :: carrier
    TYPE(dbcsr_2d_array_type), POINTER :: left_buffer_1, left_buffer_2, &
      left_buffer_calc, left_buffer_comm, right_buffer_1, right_buffer_2, &
      right_buffer_calc, right_buffer_comm
    TYPE(dbcsr_data_obj)                     :: left_data_rp, left_data_sp, &
                                                right_data_rp, right_data_sp
    TYPE(dbcsr_error_type)                   :: t_error
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

!$OMP THREADPRIVATE (carrier)


!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    !
    ALLOCATE (left_buffer_1, left_buffer_2, right_buffer_1, right_buffer_2)
    t_all = 0.0_dp
    t_gemm = 0.0_dp
    mult_id=mult_id+1

    my_use_plasma = mm_driver_blas .EQ. mm_driver_plasma
    IF (PRESENT (retain_sparsity)) THEN
       keep_sparsity = retain_sparsity
    ELSE
       keep_sparsity = .FALSE.
    ENDIF
    otf_filtering = PRESENT (filter_eps)

!$omp parallel if( .NOT. my_use_plasma ) &
!$omp default (none) &
!$omp shared (nthreads, product_matrix, error)
!$omp master
    nthreads = 1
    !$ nthreads = OMP_GET_NUM_THREADS ()
    CALL dbcsr_assert (ASSOCIATED (product_matrix%m%wms),&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Work matrices do not exist",__LINE__,error)
    CALL dbcsr_assert (SIZE (product_matrix%m%wms), "EQ", nthreads,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Work matrices not correctly sized.",__LINE__,error)
!$omp end master
!$omp end parallel

    output_unit = default_output_unit
    tstart = 0.0_dp ; tstop = 0.0_dp ; trun = 0.0_dp
    t_gemm = 0.0_dp ; t_all = 0.0_dp
    flop_total = 0
    flop_metronome=0
    left_index_send_size=0
    right_index_send_size=0
    left_data_send_size=0
    right_data_send_size=0
    left_index_recv_size=0
    right_index_recv_size=0
    left_data_recv_size=0
    right_data_recv_size=0
    trun_t = m_walltime ()
    ! Set up variables
    data_type = dbcsr_get_data_type (product_matrix)
    left_row_nimages =  left_set%image_dist%i%row_decimation
    left_row_mult =     left_set%image_dist%i%row_multiplicity
    left_col_nimages =  left_set%image_dist%i%col_decimation
    left_col_mult =     left_set%image_dist%i%col_multiplicity
    right_row_nimages = right_set%image_dist%i%row_decimation
    right_row_mult =    right_set%image_dist%i%row_multiplicity
    right_col_nimages = right_set%image_dist%i%col_decimation
    right_col_mult =    right_set%image_dist%i%col_multiplicity
    left_mp_obj    = dbcsr_distribution_mp (left_set%image_dist%i%main)
    right_mp_obj   = dbcsr_distribution_mp (right_set%image_dist%i%main)
    product_mp_obj = dbcsr_distribution_mp (product_matrix%m%dist)
    numnodes          = dbcsr_mp_numnodes (product_mp_obj)
    mynode            = dbcsr_mp_mynode (product_mp_obj)
    left_nprows       = dbcsr_mp_nprows(left_mp_obj)
    left_npcols       = dbcsr_mp_npcols(left_mp_obj)
    left_myprow       = dbcsr_mp_myprow(left_mp_obj)
    left_mypcol       = dbcsr_mp_mypcol(left_mp_obj)
    left_myfirstvrow  = dbcsr_mp_myprow(left_mp_obj)*left_row_nimages
    left_myfirstvcol  = dbcsr_mp_mypcol(left_mp_obj)*left_col_nimages
    right_nprows      = dbcsr_mp_nprows(right_mp_obj)
    right_npcols      = dbcsr_mp_npcols(right_mp_obj)
    right_myprow      = dbcsr_mp_myprow(right_mp_obj)
    right_mypcol      = dbcsr_mp_mypcol(right_mp_obj)
    right_myfirstvrow = dbcsr_mp_myprow(right_mp_obj)*right_row_nimages
    right_myfirstvcol = dbcsr_mp_mypcol(right_mp_obj)*right_col_nimages
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
    !
    ! Dummy checks
    ! left/right matching
    CALL dbcsr_assert (left_col_nimages, "EQ", right_row_mult,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Left/Right image mismatch",__LINE__,error)
    CALL dbcsr_assert (left_col_mult, "EQ", right_row_nimages,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Left/Right image mismatch",__LINE__,error)
    CALL dbcsr_assert (left_col_nimages * left_npcols,&
         "EQ", right_row_nimages * right_nprows, &
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Left/Right total mismatch",__LINE__,error)
    ! product/left matching
    CALL dbcsr_assert (left_row_mult * dbcsr_mp_nprows (product_mp_obj), &
         "EQ", left_row_nimages * left_nprows, &
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Product/Left total mismatch",__LINE__,error)
    ! product/left matching
    CALL dbcsr_assert (right_col_mult * dbcsr_mp_npcols (product_mp_obj), &
         "EQ", right_col_nimages * right_npcols, &
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Product/Right total mismatch",__LINE__,error)
    ! Limitations
    CALL dbcsr_assert (left_row_nimages, "EQ", 1,&
         dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
         "Product/Left matrix process grid mismatch",__LINE__,error)
    CALL dbcsr_assert (left_row_mult, "EQ", 1,&
         dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
         "Product/Left matrix process grid mismatch",__LINE__,error)
    CALL dbcsr_assert (right_col_nimages, "EQ", 1,&
         dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
         "Product/Right matrix process grid mismatch",__LINE__,error)
    CALL dbcsr_assert (right_col_mult, "EQ", 1,&
         dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
         "Product/Right matrix process grid mismatch",__LINE__,error)
    !
    ! Exchange size data
    ALLOCATE (my_sizes(4, MAX (left_row_nimages, right_row_nimages),&
         MAX (left_col_nimages, right_col_nimages)))
    my_sizes(:,:,:) = 0
    DO left_row_image = 1, left_row_nimages
       DO left_col_image = 1, left_col_nimages
          my_sizes(idata+ileft, left_row_image, left_col_image) &
               = dbcsr_data_get_size_referenced (&
               left_set%mats(left_row_image, left_col_image)%m%data_area)
          my_sizes(imeta+ileft, left_row_image, left_col_image) = &
               left_set%mats(left_row_image, left_col_image)%m%index&
               (dbcsr_slot_size)
       ENDDO
    ENDDO
    DO right_row_image = 1, right_row_nimages
       DO right_col_image = 1, right_col_nimages
          my_sizes(idata+iright, right_row_image, right_col_image) &
               = dbcsr_data_get_size_referenced (&
               right_set%mats(right_row_image, right_col_image)%m%data_area)
          my_sizes(imeta+iright, right_row_image, right_col_image) = &
               right_set%mats(right_row_image, right_col_image)%m%index&
               (dbcsr_slot_size)
       ENDDO
    ENDDO
    ALLOCATE (all_sizes(4, LBOUND(my_sizes,2):UBOUND(my_sizes,2),&
         LBOUND(my_sizes,3):UBOUND(my_sizes,3), 0:numnodes-1))
    CALL mp_allgather(my_sizes, all_sizes, mp_group)
    !
    ! Count the maximum possible multiplies per row for on-the-fly
    ! filtering.
    per_row_eps: IF (.NOT.otf_filtering) THEN
       ! These arrays must be valid when passed to called subroutines.
       ALLOCATE(left_norms(0),right_norms(0),row_max_epss(0), stat=stat)
       CALL dbcsr_assert (stat, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Could not allocate memory",&
            __LINE__, error=error)
    ELSE
       IF (careful_mod) THEN
          CALL dbcsr_assert ("NOT", left_set%mats(1, 1)%m%bcsc,&
               dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
               "Can not do on-the-fly filtering with CSC-indexed matrices.",&
               __LINE__, error=error)
       ENDIF
       IF (dbcsr_has_local_row_index (left_set%mats(1, 1))) THEN
          nblkrows_used = dbcsr_nblkrows_local (left_set%mats(1, 1))
       ELSE
          nblkrows_used = dbcsr_nblkrows_total (left_set%mats(1, 1))
       ENDIF
       nblkrows_total = dbcsr_nblkrows_total (left_set%mats(1, 1))
       ALLOCATE (row_max_epss (nblkrows_total), stat=stat)
       CALL dbcsr_assert (stat, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Could not allocate memory for left epsilons",&
            __LINE__, error=error)
       ALLOCATE (row_counts (nblkrows_used), stat=stat)
       CALL dbcsr_assert (stat, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Could not allocate memory for left row counts",&
            __LINE__, error=error)
       ! The summation could be done prow-locally but it would
       ! complicate the pre-row eps calculation.
       ALLOCATE (total_row_counts (nblkrows_total), stat=stat)
       CALL dbcsr_assert (stat, "EQ", 0, dbcsr_fatal_level,&
            dbcsr_internal_error,&
            routineN, "Could not allocate memory for left row counts",&
            __LINE__, error=error)
       ! Each prow member matrix (npcols * row_images) counts the
       ! blocks present in each of its rows.
       total_row_counts(:) = 0
       DO left_row_image = 1, left_row_nimages
          DO left_col_image = 1, left_col_nimages
             list_indexing =&
                  left_set%mats(left_row_image, left_col_image)%m%list_indexing
             IF (careful_mod) THEN
                IF (list_indexing) THEN
                   CALL dbcsr_assert ((left_set%mats(left_row_image, left_col_image)%m%nblks)*3, "EQ",&
                        SIZE(left_set%mats(left_row_image, left_col_image)%m%row_p),&
                        dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                        "Row count mismatch", __LINE__, error=error)
                ELSE
                   CALL dbcsr_assert (nblkrows_used+1, "EQ",&
                        SIZE(left_set%mats(left_row_image, left_col_image)%m%row_p),&
                        dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                        "Row count mismatch", __LINE__, error=error)
                ENDIF
             ENDIF
             IF (list_indexing) THEN
                CALL count_bins (&
                     left_set%mats(left_row_image, left_col_image)%m%nblks,&
                     left_set%mats(left_row_image, left_col_image)%m%row_p(1::3),&
                     nblkrows_used, row_counts)
             ELSE
                CALL dbcsr_count_row_index (&
                     left_set%mats(left_row_image, left_col_image)%m%row_p,&
                     row_counts, nblkrows_used)
             ENDIF
             IF (dbcsr_has_local_row_index (left_set%mats(left_row_image, left_col_image))) THEN
                local_rows => array_data (left_set%mats(left_row_image, left_col_image)%m%local_rows)
                CALL dbcsr_assert (SIZE(local_rows), "EQ", SIZE(row_counts),&
                     dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                     "Mismatch in number of local rows.", __LINE__, error=error)
                total_row_counts(local_rows) = total_row_counts(local_rows)&
                     + row_counts(1:nblkrows_used)
             ELSE
                total_row_counts(:) = total_row_counts(:)&
                     + row_counts(:)
             ENDIF
          ENDDO
       ENDDO
       ! The counted blocks are then summed up
       CALL mp_sum(total_row_counts, mp_group)
       ! and used to determine the maximum per-block epsilon.
       FORALL (row = 1 : nblkrows_total)
          row_max_epss (row) &
               = REAL(filter_eps&
               / REAL(MAX(1, total_row_counts(row)), KIND=KIND(row_max_epss)),&
               KIND=KIND(row_max_epss))
       END FORALL
       !
       DEALLOCATE (row_counts, stat=stat)
       CALL dbcsr_assert (stat, "EQ", 0, dbcsr_fatal_level,&
            dbcsr_internal_error, routineN,&
            "Could not deallocate memory for right matrix row counts",&
            __LINE__, error=error)
       DEALLOCATE (total_row_counts, stat=stat)
       CALL dbcsr_assert (stat, "EQ", 0, dbcsr_fatal_level,&
            dbcsr_internal_error, routineN,&
            "Could not deallocate memory for right matrix row counts",&
            __LINE__, error=error)
    ENDIF per_row_eps
    !
    ! The main transfer loop goes through the virtual rows/columns.
    ! The number of steps may be smaller if the grid dimension is very
    ! non-optimal (both left column images and right row images are >
    ! 1).
    min_nimages = MIN (left_col_nimages, right_row_nimages)
    nvirt_k = left_npcols * left_col_nimages
    nsteps_k = nvirt_k / min_nimages
    !
    ! Timings
    ALLOCATE (xfer_timings(4, 2, 3, nsteps_k))
    xfer_timings = 0.0_dp
    !
    ! Translate the all_sizes to account for pre-distribution.  This
    ! is just done to simplify lookups.
    ALLOCATE (left_sizes(2, 0:left_nprows*left_row_nimages-1, 0:nvirt_k-1))
    left_sizes = -1
    DO left_src_vcol = 0, left_col_nimages*left_npcols-1
       DO left_src_vrow = 0, left_row_nimages*left_nprows-1
          ! Calculate what was shifted.  The left_src_v{row,col} are
          ! the "source" rows/columns; the left_dst are the shifted
          ! targets where the data was placed in make_images.
          CALL image_calculator(left_set%image_dist,&
               prow = left_dst_prow, pcol = left_dst_pcol,&
               rowi = left_dst_irow, coli = left_dst_icol,&
               myvprow = left_src_vrow, myvpcol = left_src_vcol,&
               shifting = 'l', error=error)
          left_dst_p = left_pgrid (left_dst_prow, left_dst_pcol)
          left_sizes(idata, left_src_vrow, left_src_vcol) =&
               all_sizes(&
               idata+ileft, left_dst_irow, left_dst_icol, left_dst_p)
          left_sizes(imeta, left_src_vrow, left_src_vcol) =&
               all_sizes(&
               imeta+ileft, left_dst_irow, left_dst_icol, left_dst_p)
       ENDDO
    ENDDO
    !
    ALLOCATE (right_sizes(2, 0:nvirt_k-1, 0:right_npcols*right_col_nimages-1))
    right_sizes = -1
    DO right_src_vcol = 0, right_col_nimages*right_npcols-1
       DO right_src_vrow = 0, right_row_nimages*right_nprows-1
          ! Calculate what was shifted.  The right_src_v{row,col} are
          ! the "source" rows/columns; the right_dst are the shifted
          ! targets where the data was placed in make_images.
          CALL image_calculator(right_set%image_dist,&
               prow = right_dst_prow, pcol = right_dst_pcol,&
               rowi = right_dst_irow, coli = right_dst_icol,&
               myvprow = right_src_vrow, myvpcol = right_src_vcol,&
               shifting = 'r', error=error)
          right_dst_p = right_pgrid (right_dst_prow, right_dst_pcol)
          right_sizes(idata, right_src_vrow, right_src_vcol) =&
               all_sizes(&
               idata+iright, right_dst_irow, right_dst_icol, right_dst_p)
          right_sizes(imeta, right_src_vrow, right_src_vcol) =&
               all_sizes(&
               imeta+iright, right_dst_irow, right_dst_icol, right_dst_p)
       ENDDO
    ENDDO
    !
    ! Setup product work areas
    left_max_nze    = MAXVAL (all_sizes(idata+ileft, :, :, :))
    left_max_nblks  = MAXVAL (all_sizes(imeta+ileft, :, :, :))
    right_max_nze   = MAXVAL (all_sizes(idata+iright, :, :, :))
    right_max_nblks = MAXVAL (all_sizes(imeta+iright, :, :, :))
    !!
    ithread = 0
!$omp parallel default(none) &
!$omp          private (i, size_guess, size8, fill_guess, &
!$omp                   left_fill, right_fill, ithread) &
!$omp          shared (product_matrix, left_max_nze, right_max_nze) &
!$omp          shared (left_set, right_set, &
!$omp                 left_col_nimages, right_row_nimages) &
!$omp          shared (error, nthreads, keep_sparsity, mynode)
    !
    !$ ithread = OMP_GET_THREAD_NUM()
    ! The work arrays have to be setup.
    i = ithread + 1
    IF (keep_sparsity) THEN
       size_guess = product_matrix%m%wms(i)%datasize ! Should be minimal
    ELSE
       ! First we calculate the sparsities
       size8 = INT(dbcsr_nfullrows_local (left_set%mats(1,1)), KIND=int_8)&
            * INT(dbcsr_nfullcols_local (left_set%mats(1,1)), KIND=int_8)
       size8 = MAX(1_int_8,size8)
       left_fill = REAL(left_max_nze, KIND=real_8)*REAL(left_col_nimages, KIND=real_8) &
            / REAL(size8, KIND=real_8)
       size8 = INT(dbcsr_nfullrows_local (right_set%mats(1,1)), KIND=int_8)&
            * INT(dbcsr_nfullcols_local (right_set%mats(1,1)), KIND=int_8)
       size8 = MAX(1_int_8,size8)
       right_fill = REAL(right_max_nze,KIND=real_8)*REAL(right_row_nimages, KIND=real_8) &
            / REAL(size8, KIND=real_8)
       size8 = INT(dbcsr_nfullrows_local (product_matrix), KIND=int_8)&
            * INT(dbcsr_nfullcols_local (product_matrix), KIND=int_8)
       size8 = MAX(1_int_8,size8)
       IF (debug_mod .AND. mynode .EQ. 0) THEN
          WRITE(*,'(1X,A,2(1X,F12.3))')routineN//" fill orig =",&
               left_fill, right_fill
       ENDIF
       ! Old guess: fill_guess = 7 * MAX(left_fill,right_fill)
       fill_guess = 7 * MAX(left_fill,right_fill)
       ! New guess: fill_guess = 2.4_real_8 * MAX(left_fill,right_fill)
       fill_guess = MIN (1.0_real_8, MAX (0.0_real_8, fill_guess))
       IF (nthreads .GT. 1) THEN
          fill_guess = fill_guess * 3.0_real_8 / REAL(2*nthreads, KIND=real_8)
       ENDIF
       IF (debug_mod .AND. mynode .EQ. 0) THEN
          WRITE(*,*)routineN//" fill guess=", fill_guess
       ENDIF
       size_guess = MAX(product_matrix%m%wms(i)%datasize,&
            INT(REAL(size8, KIND=real_8) * fill_guess, KIND=int_4))
    ENDIF
    IF (debug_mod) &
         WRITE(*,'(1X,A,2(1X,F12.3))')routineN//" Using size guess",&
         LOG(REAL(size_guess))/LOG(10.0), LOG(REAL(size8))/LOG(10.0)
    CALL dbcsr_data_ensure_size(product_matrix%m%wms(i)%data_area,&
         size_guess,error=error)
    CALL dbcsr_data_set_size_referenced (product_matrix%m%wms(i)%data_area,&
         product_matrix%m%wms(i)%datasize)
    ! XXXXXXX a quick fix right now, allocation with size 1 might actually not be needed at all,
    !         but something expects this to be associated
    CALL ensure_array_size(product_matrix%m%wms(i)%row_i, ub=1, error=error)
    CALL ensure_array_size(product_matrix%m%wms(i)%col_i, ub=1, error=error)
    CALL ensure_array_size(product_matrix%m%wms(i)%blk_p, ub=1, error=error)
!$omp end parallel

    ! Allocate unresizable data on GPU: Multiplication parameters and
    ! A and B data.
    cuda_firstalloc: IF (mm_driver .EQ. mm_driver_cuda) THEN
       CALL dbcsr_cuda_dev_mem_alloc (params_dev(1), mm_stack_size*n_mult_params,&
            error=error)
       CALL dbcsr_cuda_dev_mem_alloc (params_dev(2), mm_stack_size*n_mult_params,&
            error=error)
       param_host_fill = 1
       param_host_xfer = 1
       param_dev_xfer = 1
       param_dev_calc = 1
       !
       CALL dbcsr_cuda_dev_mem_setup (a_dev(1), data_type)
       CALL dbcsr_cuda_dev_mem_setup (a_dev(2), data_type)
       CALL dbcsr_cuda_dev_mem_setup (b_dev(1), data_type)
       CALL dbcsr_cuda_dev_mem_setup (b_dev(2), data_type)
       CALL dbcsr_cuda_dev_mem_alloc (a_dev(1), left_max_nze, error=error)
       CALL dbcsr_cuda_dev_mem_alloc (a_dev(2), left_max_nze, error=error)
       CALL dbcsr_cuda_dev_mem_alloc (b_dev(1), right_max_nze, error=error)
       CALL dbcsr_cuda_dev_mem_alloc (b_dev(2), right_max_nze, error=error)
       data_dev_calc = 1!2
       data_dev_xfer = 1
       !
       ! Allocate C locks on GPU.  Locks should be cleared as soon as
       ! allocated.
       CALL dbcsr_cuda_dev_mem_alloc (c_locks_dev,&
            product_matrix%m%wms(1)%lastblk, error=error)
       CALL dbcsr_cuda_dev_mem_zero (c_locks_dev,&
            first=1, last=dbcsr_cuda_dev_mem_get_alloc(c_locks_dev),&
            error=error)
       !
       ! Send C data to GPU
       CALL dbcsr_assert (nthreads, "EQ", 1,&
            dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
            "Multiple threads not supported with GPU.", __LINE__, error=error)
       i = dbcsr_data_get_size (product_matrix%m%wms(1)%data_area)
       IF (verbose_acc) WRITE(*,*)routineN//" c data size", i
       !
       ! Time CUDA device memory allocation
       IF (verbose_acc) WRITE(*,*)'allocating c_dev'
       CALL dbcsr_cuda_dev_mem_setup (c_dev, data_type)
       t_xfer_c_in = -m_walltime()
       CALL dbcsr_cuda_dev_mem_alloc (c_dev, MAX(1,i), error=error)
       t_xfer_c_in = t_xfer_c_in + m_walltime()
       IF (verbose_acc) THEN
          WRITE(*,'(1X,A20,1(1X,EN12.3,1X,A))')"Device Allocation", t_xfer_c_in, "s"
       ENDIF
       !
       ! Time transfer from regular memory to GPU device memory.
       t_xfer_c_in = -m_walltime()
       IF (verbose_acc) WRITE(*,*)routineN//" copy in c", i
       CALL dbcsr_cuda_cp_host_to_dev (&
            product_matrix%m%wms(1)%data_area,&
            c_dev, error=error)
       IF (i .EQ. 0) THEN
          IF (verbose_acc) WRITE(*,*)'zeroing c_dev'
          CALL dbcsr_cuda_dev_mem_zero(c_dev, first=1, last=1, error=error)
       ENDIF
       len_xfer_c_in = i
       t_xfer_c_in = m_walltime() + t_xfer_c_in
       IF (verbose_acc) THEN
          WRITE(*,'(1X,A20,3(1X,EN12.3,1X,A))')"Cin regular", t_xfer_c_in, "s",&
               REAL(len_xfer_c_in*8,kind=dp)/t_xfer_c_in/factor_mega, "MB/s",&
               REAL(len_xfer_c_in*8), "B"
       ENDIF
       !
    ENDIF cuda_firstalloc
    !-------------
    !
    IF (debug_mod .AND. mynode .EQ. 0) THEN
       WRITE(*,*)routineN//" All sizes"
       WRITE(*,'(1X,F12.3)') LOG(REAL(all_sizes(idata, :,:,:)))/LOG(10.0)
    ENDIF
    !
    ! Setup the left buffer matrices
    !
    IF (debug_mod .AND. mynode .EQ. 0) THEN
       WRITE(*,*)routineN//" All sizes"
       WRITE(*,'(1X,F12.3)') LOG(REAL(all_sizes(idata, :,:,:)))/LOG(10.0)
    ENDIF
    CALL setup_buffer_matrices (left_buffer_1, left_row_mult, left_col_nimages,&
         left_set%mats(1,1), index_size=left_max_nblks,&
         data_size=left_max_nze, error=error)
    CALL setup_buffer_matrices (left_buffer_2, left_row_mult, left_col_nimages,&
         left_set%mats(1,1), index_size=left_max_nblks,&
         data_size=left_max_nze, error=error)
    IF (otf_filtering) THEN
       ALLOCATE (left_norms (left_max_nblks), stat=stat)
       CALL dbcsr_assert (stat, "EQ", 0, dbcsr_fatal_level,&
            dbcsr_internal_error,&
            routineN, "Could not allocate memory for left norms", __LINE__,&
            error=error)
       IF (stat .NE. 0) otf_filtering = .FALSE.
    ENDIF
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
    CALL setup_buffer_matrices (right_buffer_1, right_row_nimages, right_col_mult,&
         right_set%mats(1,1), index_size=right_max_nblks, data_size=right_max_nze,&
         error=error)
    CALL setup_buffer_matrices (right_buffer_2, right_row_nimages, right_col_mult,&
         right_set%mats(1,1), index_size=right_max_nblks, data_size=right_max_nze,&
         error=error)
    IF (otf_filtering) THEN
       ALLOCATE (right_norms (right_max_nblks), stat=stat)
       CALL dbcsr_assert (stat, "EQ", 0, dbcsr_warning_level,&
            dbcsr_internal_error,&
            routineN, "Could not allocate memory for right norms", __LINE__,&
            error=error)
       IF (stat .NE. 0) otf_filtering = .FALSE.
    ENDIF
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
!$omp parallel &
!$omp default (none) &
!$omp shared (left_buffer_calc, right_buffer_calc, product_matrix,&
!$omp         keep_sparsity, filter_eps, row_max_epss, error)
    CALL dbcsr_multrec_init(&
         left_buffer_calc%mats(1, 1)%m,&
         right_buffer_calc%mats(1, 1)%m,&
         product_matrix%m,&
         carrier,&
         keep_sparsity=keep_sparsity,&
         eps=filter_eps,&
         row_max_epss = row_max_epss,&
         error=error)
!$omp end parallel
    !
    ! Setup indexing
    CALL setup_rec_index (left_set, error)
    CALL setup_rec_index (right_set, error)
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
    ! Here is the main loop.
    !
    ! In the first loop iteration, the data is fetched from the
    ! sources. In the remaining iterations, the data are exchanged
    ! among neighbors.  In the last loop only calculations take place.
    grouped_k_index: DO metronome = 1, nsteps_k
       IF (debug_mod) WRITE(*,'(1X,A,3(1X,A,1X,I5))')routineN,&
            "step",metronome,&
            "first k",metronome*min_nimages,&
            "last k",(metronome+1)*min_nimages-1
       ! Wait for right matrix transfer completion. Wait in all but
       ! the first loop iteration.
       CALL dbcsr_error_set(routineN//"_metrocomm", error_handler2, error)
       wait_right: IF (metronome .GT. 1) THEN
          IF (debug_mod) WRITE (*,'(1X,A)')routineN//" waiting for right"
          !
          CALL xtime_set (xfer_timings(&
               idata+iright, id_send, id_waittime, metronome-1))
          CALL mp_waitall (right_data_sr)
          CALL xtime_stop (xfer_timings(&
               idata+iright, id_send, id_waittime, metronome-1))
          CALL xtime_stop (xfer_timings(&
               idata+iright, id_send, id_time, metronome-1))
          !
          CALL xtime_set (xfer_timings(&
               idata+iright, id_recv, id_waittime, metronome-1))
          CALL mp_waitall (right_data_rr)
          CALL xtime_stop (xfer_timings(&
               idata+iright, id_recv, id_waittime, metronome-1))
          CALL xtime_stop (xfer_timings(&
               idata+iright, id_recv, id_time, metronome-1))
          IF (use_combined_types) THEN
             DO v_ki = 1, right_row_nimages
                CALL mp_type_free (right_recv_type(v_ki, 1))
                CALL mp_type_free (right_send_type(v_ki, 1))
             ENDDO
          ELSE
             CALL xtime_set (xfer_timings(&
                  imeta+iright, id_send, id_waittime, metronome-1))
             CALL mp_waitall (right_index_sr)
             CALL xtime_stop (xfer_timings(&
                  imeta+iright, id_send, id_waittime, metronome-1))
             CALL xtime_stop (xfer_timings(&
                  imeta+iright, id_send, id_time, metronome-1))
             CALL xtime_set (xfer_timings(&
                  imeta+iright, id_recv, id_waittime, metronome-1))
             CALL mp_waitall (right_index_rr)
             CALL xtime_stop (xfer_timings(&
                  imeta+iright, id_recv, id_waittime, metronome-1))
             CALL xtime_stop (xfer_timings(&
                  imeta+iright, id_recv, id_time, metronome-1))
          ENDIF
       ENDIF wait_right
       CALL dbcsr_error_stop(error_handler2, error)
       ! Right matrix transfer. Transfer in all but the last loop
       ! iteration.
       xfer_right: IF (metronome .LT. nsteps_k) THEN
          DO v_ki = 0, right_row_nimages-1
             v_k = metronome*min_nimages + v_ki
             ! Calculate the process to send to.  It's the virtual
             ! process row -min_nimages up (i.e., smaller row number)
             ! from me.
             CALL image_calculator (right_set%image_dist,&
                  prow=right_send_prow, rowi=right_send_irow,&   ! output
                  pcol=right_send_pcol, coli=right_send_icol,&   ! output
                  vprow=right_send_vrow, vpcol=right_send_vcol,& ! output
                  ! myvprow goes through all of my (process row) images
                  myvprow=v_ki+right_myfirstvrow,&
                  myvpcol=right_myfirstvcol,& ! nothing happens in the columns
                  ! send to process min_nimages up in the grid
                  vprow_shift=-min_nimages,&
                  shifting='0', error=error)
             ! Calculate which data I send.
             CALL image_calculator (right_set%image_dist,&
                  prow=right_dst_prow, rowi=right_dst_irow,&
                  pcol=right_dst_pcol, coli=right_dst_icol,&
                  vprow=right_dst_vrow, vpcol=right_dst_vcol,&
                  ! myvprows goes through all of my (process row) images
                  myvprow=v_ki+right_myfirstvrow,&
                  myvpcol=right_myfirstvcol,& ! nothing happens in the columns
                  ! send what I got from min_nimages down, appropriate
                  ! to the metronome tick
                  vprow_shift=-min_nimages + metronome*min_nimages,&
                  ! This is with relative shifting.
                  shifting='R', error=error)
             right_dst_p = right_pgrid(right_dst_prow, right_dst_pcol)
             CALL dbcsr_data_set_pointer(&
                  area=right_data_sp,&
                  rsize=right_sizes(idata, right_dst_vrow, right_dst_vcol),&
                  csize=1,&
                  pointee=right_buffer_calc%mats(v_ki+1, 1)%m%data_area)
             right_index_sp => right_buffer_calc%mats(&
                  v_ki+1, 1&
                  )%m%index(1:&
                  right_sizes(imeta, right_dst_vrow, right_dst_vcol))
             !
             ! Calculate the process to receive from
             CALL image_calculator (right_set%image_dist,&
                  prow=right_recv_prow, rowi=right_recv_irow,&
                  pcol=right_recv_pcol, coli=right_recv_icol,&
                  vprow=right_recv_vrow, vpcol=right_recv_vcol,&
                  myvprow=v_ki+right_myfirstvrow,&
                  myvpcol=right_myfirstvcol,&
                  vprow_shift=+min_nimages,& ! just the opposite as "send to"
                  shifting='0', error=error)
             ! Calculate which data I receive
             CALL image_calculator (right_set%image_dist,&
                  prow=right_src_prow, rowi=right_src_irow,&
                  pcol=right_src_pcol, coli=right_src_icol,&
                  vprow=right_src_vrow, vpcol=right_src_vcol,&
                  myvprow=v_ki+right_myfirstvrow,&
                  myvpcol=right_myfirstvcol,&
                  ! receive window moves with the metronome
                  vprow_shift=metronome*min_nimages,&
                  shifting='R', error=error)
             !
             right_src_p = right_pgrid(right_src_prow, right_src_pcol)
             CALL dbcsr_data_set_pointer(&
                  area=right_data_rp,&
                  rsize=right_sizes(idata, right_src_vrow, right_src_vcol),&
                  csize=1,&
                  pointee=right_buffer_comm%mats(v_ki+1, 1)%m%data_area)
             right_index_rp => right_buffer_comm%mats(&
                     v_ki+1, 1&
                  )%m%index(1:&
                     right_sizes(imeta, right_src_vrow, right_src_vcol))
             !
             right_send_p = right_pgrid (right_send_prow, right_send_pcol)
             right_recv_p = right_pgrid (right_recv_prow, right_recv_pcol)
             ! These are column-communicator relative
             IF (dbcsr_mp_has_subgroups (right_mp_obj)) THEN
                right_send_p = right_send_prow
                right_recv_p = right_recv_prow
                grp = dbcsr_mp_my_col_group (right_mp_obj)
             ELSE
                grp = dbcsr_mp_group (right_mp_obj)
             ENDIF
             !
             CALL dbcsr_error_set(routineN//"_metrocomm", error_handler2, error)
             IF (use_combined_types) THEN
                right_send_subtypes(1) = dbcsr_mp_type_from_anytype (right_data_sp)
                right_send_subtypes(2) = mp_type_make (right_index_sp)
                right_recv_subtypes(1) = dbcsr_mp_type_from_anytype (right_data_rp)
                right_recv_subtypes(2) = mp_type_make (right_index_rp)
                right_send_type(v_ki+1, 1) = mp_type_make (right_send_subtypes)
                right_recv_type(v_ki+1, 1) = mp_type_make (right_recv_subtypes)
                CALL xtime_set(xfer_timings(idata+iright, id_send, id_time, metronome))
                CALL mp_isend (right_send_type(v_ki+1, 1), right_send_p,&
                     grp, right_data_sr(v_ki+1, 1), tag=right_dst_vrow)
                CALL xtime_set(xfer_timings(idata+iright, id_recv, id_time, metronome))
                CALL mp_irecv (right_recv_type(v_ki+1, 1), right_recv_p,&
                     grp, right_data_rr(v_ki+1, 1), tag=right_src_vrow)
             ELSE
                CALL xtime_set(xfer_timings(idata+iright, id_send, id_time, metronome))
                CALL dbcsr_isend_any (right_data_sp, right_send_p,&
                     grp, right_data_sr(v_ki+1, 1), tag=right_dst_vrow,&
                     error=error)
                CALL xtime_set(xfer_timings(idata+iright, id_recv, id_time, metronome))
                CALL dbcsr_irecv_any (right_data_rp, right_recv_p,&
                     grp, right_data_rr(v_ki+1, 1), tag=right_src_vrow,&
                     error=error)
                CALL xtime_set(xfer_timings(imeta+iright, id_send, id_time, metronome))
                CALL mp_isend (right_index_sp, right_send_p,&
                     grp, right_index_sr(v_ki+1, 1), tag=right_dst_vrow)
                CALL xtime_set(xfer_timings(imeta+iright, id_recv, id_time, metronome))
                CALL mp_irecv (right_index_rp, right_recv_p,&
                     grp, right_index_rr(v_ki+1, 1), tag=right_src_vrow)
             ENDIF
             xfer_timings(idata+iright, id_send, id_bytes, metronome) = &
                xfer_timings(idata+iright, id_send, id_bytes, metronome) + &
                dbcsr_data_get_size(right_data_sp)
             xfer_timings(idata+iright, id_recv, id_bytes, metronome) = &
                xfer_timings(idata+iright, id_recv, id_bytes, metronome) + &
                dbcsr_data_get_size(right_data_rp)
             xfer_timings(imeta+iright, id_send, id_bytes, metronome) = &
                xfer_timings(imeta+iright, id_send, id_bytes, metronome) + &
                SIZE(right_index_sp)
             xfer_timings(imeta+iright, id_recv, id_bytes, metronome) = &
                xfer_timings(imeta+iright, id_recv, id_bytes, metronome) + &
                SIZE(right_index_rp)
             IF (excessive_output) THEN
                right_data_send_size = right_data_send_size +&
                     dbcsr_data_get_size(right_data_sp)
                right_data_recv_size = right_data_send_size +&
                     dbcsr_data_get_size(right_data_rp)
                right_index_send_size = right_index_send_size +&
                     SIZE(right_index_sp)
                right_index_recv_size = right_index_send_size +&
                     SIZE(right_index_rp)
             ENDIF
             CALL dbcsr_error_stop(error_handler2, error)
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
       CALL dbcsr_error_set(routineN//"_metrocomm", error_handler2, error)
       wait_left: IF (metronome .GT. 1) THEN
          IF (debug_mod) WRITE (*,'(1X,A)')routineN//" waiting for left"
          CALL xtime_set (xfer_timings(&
               idata+ileft, id_send, id_waittime, metronome-1))
          CALL mp_waitall (left_data_sr)
          CALL xtime_stop (xfer_timings(&
               idata+ileft, id_send, id_waittime, metronome-1))
          CALL xtime_stop (xfer_timings(&
               idata+ileft, id_send, id_time, metronome-1))
          !
          CALL xtime_set (xfer_timings(&
               idata+ileft, id_recv, id_waittime, metronome-1))
          CALL mp_waitall (left_data_rr)
          CALL xtime_stop (xfer_timings(&
               idata+ileft, id_recv, id_waittime, metronome-1))
          CALL xtime_stop (xfer_timings(&
               idata+ileft, id_recv, id_time, metronome-1))
          IF (use_combined_types) THEN
             DO v_ki = 1, left_col_nimages
                CALL mp_type_free (left_send_type(1, v_ki))
                CALL mp_type_free (left_recv_type(1, v_ki))
             ENDDO
          ELSE
             CALL xtime_set (xfer_timings(&
                  imeta+ileft, id_send, id_waittime, metronome-1))
             CALL mp_waitall (left_index_sr)
             CALL xtime_stop (xfer_timings(&
                  imeta+ileft, id_send, id_waittime, metronome-1))
             CALL xtime_stop (xfer_timings(&
                  imeta+ileft, id_send, id_time, metronome-1))
             CALL xtime_set (xfer_timings(&
                  imeta+ileft, id_recv, id_waittime, metronome-1))
             CALL mp_waitall (left_index_rr)
             CALL xtime_stop (xfer_timings(&
                  imeta+ileft, id_recv, id_waittime, metronome-1))
             CALL xtime_stop (xfer_timings(&
                  imeta+ileft, id_recv, id_time, metronome-1))
          ENDIF
       ENDIF wait_left
       CALL dbcsr_error_stop(error_handler2, error)
       ! Left matrix transfer. Transfer in all but the last loop
       ! iteration.
       xfer_left: IF (metronome .LT. nsteps_k) THEN
          DO v_ki = 0, left_col_nimages-1
             v_k = metronome*min_nimages + v_ki
             ! Calculate the process to send to.
             CALL image_calculator (left_set%image_dist,&
                  prow=left_send_prow, rowi=left_send_irow,&   ! output
                  pcol=left_send_pcol, coli=left_send_icol,&   ! output
                  vprow=left_send_vrow, vpcol=left_send_vcol,& ! output
                  myvprow=left_myfirstvrow,& ! nothing happens in the rows
                  ! go through all my column images
                  myvpcol=v_ki+left_myfirstvcol,& 
                  ! send to process min_nimages left in the grid
                  vpcol_shift=-min_nimages,&
                  shifting='0', error=error)
             ! Calculate which data I send.
             CALL image_calculator (left_set%image_dist,&
                  prow=left_dst_prow, rowi=left_dst_irow,&
                  pcol=left_dst_pcol, coli=left_dst_icol,&
                  vprow=left_dst_vrow, vpcol=left_dst_vcol,&
                  myvprow=left_myfirstvrow,&
                  ! go through all my column images
                  myvpcol=v_ki+left_myfirstvcol,&
                  ! send what I got from min_nimages left, appropriate
                  ! to the metronome tick
                  vpcol_shift=-min_nimages + metronome*min_nimages,&
                  ! This is with relative shifting.
                  shifting='L', error=error)
             !
             left_dst_p = left_pgrid (left_dst_prow, left_dst_pcol)
             CALL dbcsr_data_set_pointer(&
                  area=left_data_sp,&
                  rsize=left_sizes(idata, left_dst_vrow, left_dst_vcol),&
                  csize=1,&
                  pointee=left_buffer_calc%mats(1, v_ki+1)%m%data_area)
             left_index_sp => left_buffer_calc%mats(&
                     1, v_ki+1&
                  )%m%index(1:&
                     left_sizes(imeta, left_dst_vrow, left_dst_vcol))
             !
             ! Calculate the process to receive from
             CALL image_calculator (left_set%image_dist,&
                  prow=left_recv_prow, rowi=left_recv_irow,&
                  pcol=left_recv_pcol, coli=left_recv_icol,&
                  vprow=left_recv_vrow, vpcol=left_recv_vcol,&
                  myvprow=left_myfirstvrow,&
                  myvpcol=v_ki+left_myfirstvcol,&
                  vpcol_shift=+min_nimages,& ! just the opposite as "send to"
                  shifting='0', error=error)
             ! Calculate which data I receive
             CALL image_calculator (left_set%image_dist,&
                  prow=left_src_prow, rowi=left_src_irow,&
                  pcol=left_src_pcol, coli=left_src_icol,&
                  vprow=left_src_vrow, vpcol=left_src_vcol,&
                  myvprow=left_myfirstvrow,&
                  myvpcol=v_ki+left_myfirstvcol,&
                  ! receive window moves with the metronome
                  vpcol_shift=metronome*min_nimages,&
                  shifting='L', error=error)
             !
             left_src_p = left_pgrid (left_src_prow, left_src_pcol)
             CALL dbcsr_data_set_pointer(&
                  area=left_data_rp,&
                  rsize=left_sizes(idata, left_src_vrow, left_src_vcol),&
                  csize=1,&
                  pointee=left_buffer_comm%mats(1, v_ki+1)%m%data_area)
             left_index_rp => left_buffer_comm%mats(&
                     1, v_ki+1&
                  )%m%index(1:&
                     left_sizes(imeta, left_src_vrow, left_src_vcol))
             !
             left_send_p = left_pgrid (left_send_prow, left_send_pcol)
             left_recv_p = left_pgrid (left_recv_prow, left_recv_pcol)
             ! These are column-communicator relative
             IF (dbcsr_mp_has_subgroups (left_mp_obj)) THEN
                left_send_p = left_send_pcol
                left_recv_p = left_recv_pcol
                grp = dbcsr_mp_my_row_group (left_mp_obj)
             ELSE
                grp = dbcsr_mp_group (left_mp_obj)
             ENDIF
             !
             CALL dbcsr_error_set(routineN//"_metrocomm", error_handler2, error)
             IF (use_combined_types) THEN
                left_send_subtypes(1) = dbcsr_mp_type_from_anytype (left_data_sp)
                left_send_subtypes(2) = mp_type_make (left_index_sp)
                left_recv_subtypes(1) = dbcsr_mp_type_from_anytype (left_data_rp)
                left_recv_subtypes(2) = mp_type_make (left_index_rp)
                left_send_type(1, v_ki+1) = mp_type_make (left_send_subtypes)
                left_recv_type(1, v_ki+1) = mp_type_make (left_recv_subtypes)
                CALL xtime_set(xfer_timings(idata+ileft, id_send, id_time, metronome))
                CALL mp_isend (left_send_type(1, v_ki+1), left_send_p,&
                     grp, left_data_sr(1, v_ki+1), tag=left_dst_vcol)
                CALL xtime_set(xfer_timings(idata+ileft, id_recv, id_time, metronome))
                CALL mp_irecv (left_recv_type(1, v_ki+1), left_recv_p,&
                     grp, left_data_rr(1, v_ki+1), tag=left_src_vcol)
             ELSE
                CALL xtime_set(xfer_timings(idata+ileft, id_send, id_time, metronome))
                CALL dbcsr_isend_any (left_data_sp, left_send_p,&
                     grp, left_data_sr(1, v_ki+1), tag=left_dst_vcol,&
                     error=error)
                CALL xtime_set(xfer_timings(idata+ileft, id_recv, id_time, metronome))
                CALL dbcsr_irecv_any (left_data_rp, left_recv_p,&
                     grp, left_data_rr(1, v_ki+1), tag=left_src_vcol,&
                     error=error)
                CALL xtime_set(xfer_timings(imeta+ileft, id_send, id_time, metronome))
                CALL mp_isend (left_index_sp, left_send_p,&
                     grp, left_index_sr(1, v_ki+1), tag=left_dst_vcol)
                CALL xtime_set(xfer_timings(imeta+ileft, id_recv, id_time, metronome))
                CALL mp_irecv (left_index_rp, left_recv_p,&
                     grp, left_index_rr(1, v_ki+1), tag=left_src_vcol)
             ENDIF
             xfer_timings(idata+ileft, id_send, id_bytes, metronome) = &
                xfer_timings(idata+ileft, id_send, id_bytes, metronome) + &
                dbcsr_data_get_size(left_data_sp)
             xfer_timings(idata+ileft, id_recv, id_bytes, metronome) = &
                xfer_timings(idata+ileft, id_recv, id_bytes, metronome) + &
                dbcsr_data_get_size(left_data_rp)
             xfer_timings(imeta+ileft, id_send, id_bytes, metronome) = &
                xfer_timings(imeta+ileft, id_send, id_bytes, metronome) + &
                SIZE(left_index_sp)
             xfer_timings(imeta+ileft, id_recv, id_bytes, metronome) = &
                xfer_timings(imeta+ileft, id_recv, id_bytes, metronome) + &
                SIZE(left_index_rp)
             IF (excessive_output) THEN
                left_data_send_size = left_data_send_size +&
                     dbcsr_data_get_size(left_data_sp)
                left_data_recv_size = left_data_send_size +&
                     dbcsr_data_get_size(left_data_rp)
                left_index_send_size = left_index_send_size +&
                     SIZE(left_index_sp)
                left_index_recv_size = left_index_send_size +&
                     SIZE(left_index_rp)
             ENDIF
             CALL dbcsr_error_stop(error_handler2, error)
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
             IF (debug_mod) THEN
                CALL dbcsr_print(left_buffer_calc%mats(1, v_ki+1), nodata=.TRUE., error=error)
                CALL dbcsr_print(right_buffer_calc%mats(v_ki+1, 1), nodata=.TRUE., error=error)
             ENDIF
             !
             tmp_aix_bix: IF (mm_driver .EQ. mm_driver_cuda) THEN
                t_xfer_a = t_xfer_a - m_walltime()
                IF (verbose_acc) WRITE(*,*)'copy in a', v_ki, left_max_nze
                CALL dbcsr_cuda_cp_host_to_dev(&
                     left_buffer_calc%mats(1,v_ki+1)%m%data_area,&
                     a_dev(data_dev_xfer), error=error)
                len_xfer_a = len_xfer_a + left_max_nze
                t_xfer_a = t_xfer_a + m_walltime()
                !
                t_xfer_b = t_xfer_b - m_walltime()
                IF (verbose_acc) WRITE(*,*)'copy in b', v_ki, right_max_nze
                CALL dbcsr_cuda_cp_host_to_dev(&
                     right_buffer_calc%mats(v_ki+1,1)%m%data_area,&
                     b_dev(data_dev_xfer), error=error)
                len_xfer_b = len_xfer_b + right_max_nze
                t_xfer_b = t_xfer_b + m_walltime()
             ENDIF tmp_aix_bix
             !
             IF (.FALSE.) WRITE(*,*)routineN//" TICK", v_ki
             IF (right_buffer_calc%mats(v_ki+1, 1)%m%local_indexing) THEN
                ! Since the right matrix is shifted vertically, the
                ! received data always has different notions of "local
                ! rows".  Thus the local_rows and global_rows must be
                ! recalculated.
                CALL dbcsr_reset_vlocals (right_buffer_calc%mats(v_ki+1, 1),&
                     right_set%image_dist, error=error)
             ENDIF
             IF (left_buffer_calc%mats(1, v_ki+1)%m%local_indexing) THEN
                ! Since the right matrix is shifted vertically, the
                ! received data always has different notions of "local
                ! rows".  Thus the local_rows and global_rows must be
                ! recalculated.
                CALL dbcsr_reset_vlocals (left_buffer_calc%mats(1, v_ki+1),&
                     left_set%image_dist, error=error)
             ENDIF
             ! Sets the local right-matrix columns
             IF (otf_filtering) THEN
                left_norms(:) = HUGE(left_norms(1))
                right_norms(:) = HUGE(right_norms(1))
                CALL calculate_norms(right_buffer_calc%mats(v_ki+1, 1),&
                     right_norms, error=error)
                CALL calculate_norms(left_buffer_calc%mats(1, v_ki+1),&
                     left_norms, error=error)
             ENDIF
             !
             tstart = m_walltime ()
             flop_single = 0
!$omp parallel if( .NOT. my_use_plasma ) &
!$omp default (none) &
!$omp shared (left_buffer_calc, right_buffer_calc, &
!$omp         v_ki, &
!$omp         product_matrix, &
!$omp         filter_eps, right_norms, left_norms, row_max_epss, &
!$omp         keep_sparsity, error) &
!$omp reduction (+: flop_single, t_all, t_gemm)
             CALL dbcsr_multrec_sparse(&
                  left_buffer_calc%mats(1, v_ki+1)%m,&
                  right_buffer_calc%mats(v_ki+1, 1)%m,&
                  product_matrix%m,&
                  carrier,&
                  keep_sparsity=keep_sparsity,&
                  flop=flop_single,&
                  a_norms=left_norms, b_norms=right_norms,&
                  eps=filter_eps,&
                  row_max_epss = row_max_epss, &
                  t_all=t_all, t_gemm=t_gemm, error=error)
!$omp end parallel
             flop_total = flop_total + flop_single
             flop_metronome=flop_metronome+flop_single
             tstop = m_walltime ()
             IF (tstop-tstart .EQ. 0) tstop = tstart+EPSILON(tstart)
             IF ((output_unit>0) .AND. detailed_timing) THEN
                WRITE(output_unit,'(1X,A,F9.4,A,EN12.4,A,1X,EN12.4,1X,A)')&
                     "Segment Local Multiplication time ",tstop-tstart," and ",&
                     (REAL(flop_single, dp)/1000000.0_dp)/(tstop-tstart),&
                     " MFLOP/s", REAL(flop_single), "FLOP"
                t_process_stack = MAX (EPSILON (t_process_stack), t_process_stack)
                WRITE(output_unit,'(1X,A,F9.4,A,EN12.4,A,1X,EN12.4,1X,A)')&
                     " Stack Process Multiplication time ",t_process_stack," and ",&
                     (REAL(flop_single, dp)/1000000.0_dp)/(t_process_stack),&
                     " MFLOP/s", REAL(flop_single), "FLOP"
                t_calc_step = MAX (EPSILON (t_calc_step), t_calc_step)
                WRITE(output_unit,'(1X,A,F9.4,A,EN12.4,A,1X,EN12.4,1X,A)')&
                     "Kernel Process Multiplication time ",t_calc_step," and ",&
                     (REAL(flop_single, dp)/1000000.0_dp)/(t_calc_step),&
                     " MFLOP/s", REAL(flop_single), "FLOP"
             ENDIF
             trun = trun + (tstop - tstart)
          ENDDO

          IF (excessive_output) THEN
             WRITE(1000000+mynode,*) mult_id,&
                  metronome,flop_metronome,&
                  left_index_send_size,right_index_send_size, &
                  left_data_send_size,right_data_send_size
          ENDIF
          flop_metronome=0
          left_index_send_size=0
          right_index_send_size=0
          left_data_send_size=0
          right_data_send_size=0

       ENDIF calc_case_left
       IF (metronome .EQ. 1) THEN
          left_buffer_calc => left_buffer_1
          right_buffer_calc => right_buffer_1
       ENDIF
       CALL dbcsr_switch_sets (left_buffer_calc, left_buffer_comm)
       CALL dbcsr_switch_sets (right_buffer_calc, right_buffer_comm)
    ENDDO grouped_k_index
    trun_t = m_walltime () - trun_t
    IF (excessive_output) CALL m_flush(1000000+mynode)
    IF (trun_t .EQ. 0) trun_t = 0.000001_dp
    IF (trun .EQ. 0) trun = 0.000001_dp
    IF ((output_unit>0) .AND. detailed_timing) THEN
       WRITE(output_unit,'(1X,A,F9.4,1X,A,EN12.4,1X,A)')&
            "  Total Local Multiplication time ",trun,"and ",&
            (REAL(flop_total, dp)/1000000.0_dp)/trun,&
            "MFLOP/s"
       !WRITE(output_unit,'(1X,A,F9.4,1X,A,F9.4,1X,F9.4,"%")')&
       !     "  Total Local xGEMM time          ",t_gemm,&
       !     "index time ", t_all - t_gemm, 100.0_dp*(t_all-t_gemm)/t_all
       WRITE(output_unit,'(1X,A,F9.4,1X,A,EN12.4,1X,A)')&
            "        Total Multiplication time ",trun_t,"and ",&
            (REAL(flop_total, dp)/1000000.0_dp)/trun_t,&
            "MFLOP/s"
    ENDIF
    !
    ! Finalize multiplication
    !
!$omp parallel &
!$omp default (none) &
!$omp shared (product_matrix, error)
    CALL dbcsr_multrec_finalize (carrier, product_matrix%m, error)
!$omp end parallel
    !
    ! Transfer data from GPU to host.
    cuda_copyout_c: IF (mm_driver .EQ. mm_driver_cuda) THEN
       CALL dbcsr_data_ensure_size(product_matrix%m%wms(1)%data_area,&
            carrier%datasize, error=error)
       t_xfer_c_out = -m_walltime()
       CALL dbcsr_cuda_cp_dev_to_host(c_dev,&
            product_matrix%m%wms(1)%data_area,&
            error=error)
       len_xfer_c_out = carrier%datasize
       t_xfer_c_out = t_xfer_c_out + m_walltime()
       IF (verbose_acc) THEN
          WRITE(*,'(1X,A20,3(1X,EN12.3,1X,A))')"Cout regular",&
               t_xfer_c_out, "s",&
               REAL(len_xfer_c_out*8,kind=dp)/t_xfer_c_out/factor_mega, "MB/s",&
               REAL(len_xfer_c_out*8), "B"
       ENDIF
    ENDIF cuda_copyout_c
    !
    IF (ALLOCATED (right_norms)) THEN
       DEALLOCATE (right_norms, stat=stat)
       CALL dbcsr_assert (stat, "EQ", 0, dbcsr_fatal_level,&
            dbcsr_internal_error,&
            routineN, "Could not deallocate memory for right norms", __LINE__,&
            error=error)
    ENDIF
    IF (ALLOCATED (left_norms)) THEN
       DEALLOCATE (left_norms, stat=stat)
       CALL dbcsr_assert (stat, "EQ", 0, dbcsr_fatal_level,&
            dbcsr_internal_error,&
            routineN, "Could not deallocate memory for left norms", __LINE__,&
            error=error)
    ENDIF
    IF (ALLOCATED (row_max_epss)) THEN
       DEALLOCATE (row_max_epss, stat=stat)
       CALL dbcsr_assert (stat, "EQ", 0, dbcsr_fatal_level,&
            dbcsr_internal_error,&
            routineN, "Could not deallocate memory for row block epsilons",&
            __LINE__,&
            error=error)
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
    cuda_final_release: IF (mm_driver .EQ. mm_driver_cuda) THEN
       CALL dbcsr_cuda_dev_mem_dealloc (c_locks_dev, error=error)
       CALL dbcsr_cuda_dev_mem_dealloc (c_dev, error=error)
       !
       CALL dbcsr_cuda_dev_mem_dealloc (params_dev(1), error=error)
       CALL dbcsr_cuda_dev_mem_dealloc (params_dev(2), error=error)
       !
       CALL dbcsr_cuda_dev_mem_dealloc (a_dev(1), error=error)
       CALL dbcsr_cuda_dev_mem_dealloc (a_dev(2), error=error)
       CALL dbcsr_cuda_dev_mem_dealloc (b_dev(1), error=error)
       CALL dbcsr_cuda_dev_mem_dealloc (b_dev(2), error=error)
    ENDIF cuda_final_release
    !

    !
    t_error = error
    !
    IF (debug_mod) THEN
       v_ki = 0
       DO i = 1, product_matrix%m%nblks
          v_ki = MAX(v_ki, ABS(product_matrix%m%blk_p(i)))
       ENDDO
       WRITE(*,*)routineN//" Actual final size",&
            LOG(REAL(dbcsr_data_get_size(product_matrix%m%data_area)))/LOG(10.0),&
            LOG(REAL(v_ki))/LOG(10.0)
    ENDIF
    IF (debug_mod) THEN
       checksum = dbcsr_checksum (product_matrix, error=error)
       IF ((output_unit>0)) THEN
          WRITE(output_unit,'(1X,A,1X,F9.4)')"Product Checksum=",checksum
       ENDIF
    ENDIF
    !
    flop = flop_total
    DEALLOCATE (left_buffer_1, left_buffer_2, right_buffer_1, right_buffer_2)
    !
    ! Handle timings
    xfer_timings(imeta+ileft,:,id_bytes,:) =&
         xfer_timings(imeta+ileft,:,id_bytes,:) * int_4_size
    xfer_timings(imeta+iright,:,id_bytes,:) =&
         xfer_timings(imeta+iright,:,id_bytes,:) * int_4_size
    SELECT CASE (dbcsr_get_data_type (product_matrix))
    CASE (dbcsr_type_real_4)
       i = real_4_size
    CASE (dbcsr_type_real_8)
       i = real_8_size
    CASE (dbcsr_type_complex_4)
       i = real_4_size * 2
    CASE (dbcsr_type_complex_8)
       i = real_8_size * 2
    END SELECT
    xfer_timings(idata+ileft,:,id_bytes,:) =&
         xfer_timings(idata+ileft,:,id_bytes,:) * i
    xfer_timings(idata+iright,:,id_bytes,:) =&
         xfer_timings(idata+iright,:,id_bytes,:) * i
    IF (time_xfers) &
         CALL print_xfer_timings(xfer_timings)
    DEALLOCATE (xfer_timings)
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

    INTEGER                                  :: col_image, data_mem_type, &
                                                error_handler, &
                                                index_mem_type, row_image
    INTEGER, DIMENSION(:), POINTER           :: i1
    LOGICAL                                  :: use_acc_mem

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    !
    ! Determine which type of memory to use for the buffers.  The
    ! first priority is given to CUDA, then MPI-allocated.
    IF (mm_driver .EQ. mm_driver_cuda) THEN
       use_acc_mem = use_CUDA_host_pinned_memory
    ELSE
       use_acc_mem = .FALSE.
    ENDIF
    CALL dbcsr_assert (use_acc_mem, "IMP",&
         mm_driver .EQ. mm_driver_cuda,&
         dbcsr_fatal_level, dbcsr_caller_error, routineN,&
         "Should not use CUDA host pinned memory without CUDA support.",&
         __LINE__, error=error)
    IF (use_MPI_memory) THEN
       data_mem_type = dbcsr_memory_MPI
       index_mem_type = dbcsr_memory_MPI
    ELSE
       data_mem_type = dbcsr_memory_default
       index_mem_type = dbcsr_memory_default
    ENDIF
    IF (use_acc_mem) THEN
       data_mem_type = dbcsr_memory_CUDA_host_pinned
    ENDIF
    !
    CALL dbcsr_image_dist_init (buffer_set%image_dist, error=error)
    ALLOCATE (buffer_set%mats(buff_rows, buff_cols))
    DO row_image = 1, buff_rows
       DO col_image = 1, buff_cols
          CALL dbcsr_init(buffer_set%mats(row_image, col_image))
          ! Dummy allocation only needed for NAG (at least for 5.1(327))
          ALLOCATE(i1(10000))
          CALL dbcsr_create(buffer_set%mats(row_image, col_image),&
               template = source_matrix,&
               name = TRIM("Buffer of "//TRIM(source_matrix%m%name)),&
               nblks = index_size, nze = data_size, &
               data_memory_type = data_mem_type,&
               index_memory_type = index_mem_type,&
               error = error)
          ! Dummy allocation only needed for NAG (at least for 5.1(327))
          DEALLOCATE(i1)
          CALL dbcsr_data_ensure_size (&
               buffer_set%mats(row_image, col_image)%m%data_area,&
               data_size, nocopy=.TRUE.,error=error)
          CALL ensure_array_size (&
               buffer_set%mats(row_image, col_image)%m%index,&
               ub=index_size, nocopy=.TRUE.,&
               memory_type=dbcsr_get_index_memory_type(buffer_set%mats(row_image, col_image)),&
               error=error)
          buffer_set%mats(row_image, col_image)%m%negate_real&
               = source_matrix%m%negate_real
          buffer_set%mats(row_image, col_image)%m%negate_imaginary&
               = source_matrix%m%negate_imaginary
          buffer_set%mats(row_image, col_image)%m%local_indexing &
               = source_matrix%m%local_indexing
          buffer_set%mats(row_image, col_image)%m%list_indexing &
               = source_matrix%m%list_indexing
          !
          IF (source_matrix%m%has_local_rows) THEN
             buffer_set%mats(row_image, col_image)%m%local_rows &
                                   = source_matrix%m%local_rows
             CALL array_hold (buffer_set%mats(row_image, col_image)%m%local_rows)
             buffer_set%mats(row_image, col_image)%m%has_local_rows = .TRUE.
          ENDIF
          IF (source_matrix%m%has_global_rows) THEN
             buffer_set%mats(row_image, col_image)%m%global_rows &
                                   = source_matrix%m%global_rows
             CALL array_hold (buffer_set%mats(row_image, col_image)%m%global_rows)
             buffer_set%mats(row_image, col_image)%m%has_global_rows = .TRUE.
          ENDIF
          IF (source_matrix%m%has_local_cols) THEN
             buffer_set%mats(row_image, col_image)%m%local_cols &
                                   = source_matrix%m%local_cols
             CALL array_hold (buffer_set%mats(row_image, col_image)%m%local_cols)
             buffer_set%mats(row_image, col_image)%m%has_local_cols = .TRUE.
          ENDIF
          IF (source_matrix%m%has_global_cols) THEN
             buffer_set%mats(row_image, col_image)%m%global_cols &
                                   = source_matrix%m%global_cols
             CALL array_hold (buffer_set%mats(row_image, col_image)%m%global_cols)
             buffer_set%mats(row_image, col_image)%m%has_global_cols = .TRUE.
          ENDIF
          IF (source_matrix%m%local_indexing .AND. careful_mod) THEN
             CALL dbcsr_assert (array_exists (source_matrix%m%local_rows),&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Local rows should exist.", __LINE__, error=error)
             CALL dbcsr_assert (array_exists (source_matrix%m%global_rows),&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Global rows should exist.", __LINE__, error=error)
             !
             CALL dbcsr_assert (array_exists (source_matrix%m%local_cols),&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Local cols should exist.", __LINE__, error=error)
             CALL dbcsr_assert (array_exists (source_matrix%m%global_cols),&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Global cols should exist.", __LINE__, error=error)
          ENDIF
       ENDDO
    ENDDO
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE setup_buffer_matrices


! *****************************************************************************
!> \brief Multiplies two DBCSR matrices using recursive algorithm
!>
!> This routine sets up the multiplication.  Specifically, it <ul>
!> <li> verifies input sanity
!> <li> converts everything into "local indexing"
!> </ul>
!>
!> \param[in] left, right     left and right DBCSR matrices
!> \param[in,out] product     resulting DBCSR product matrix
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix, default is no
!> \param[in,out] t_all       accumulated time spent in entire routine
!> \param[in,out] t_gemm      accumulated time spent just for multiplication
!> \param[in] a_norms         (optional) norms of left-matrix blocks
!> \param[in] b_norms         (optional) norms of right-matrix blocks
!> \param[in] eps             (optional) on-the-fly filtering epsilon
!> \param[in] row_max_epss  (optional) on-the-fly filtering epsilon per block
!> \param[out] flop           (optional) number of effective double-precision
!>                            floating point operations performed
! *****************************************************************************
  SUBROUTINE dbcsr_multrec_sparse(left, right, product, carrier, flop,&
       keep_sparsity, &
       t_all, t_gemm, a_norms, b_norms, eps, row_max_epss, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(carrier_type), INTENT(inout)        :: carrier
    INTEGER(KIND=int_8), INTENT(OUT)         :: flop
    LOGICAL, INTENT(IN)                      :: keep_sparsity
    REAL(KIND=dp), INTENT(INOUT)             :: t_all, t_gemm
    REAL(kind=sp), DIMENSION(:), &
      INTENT(in), TARGET                     :: a_norms, b_norms
    REAL(kind=real_8), INTENT(in), OPTIONAL  :: eps
    REAL(kind=sp), DIMENSION(:), &
      INTENT(in), TARGET                     :: row_max_epss
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multrec_sparse', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handler, ithread, &
                                                nthreads, t_a_f, t_a_l, &
                                                t_b_f, t_b_l
    INTEGER, DIMENSION(:), POINTER           :: k_locals
    LOGICAL                                  :: local_cols, local_indexing, &
                                                local_rows
    REAL(KIND=dp)                            :: epoch

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    !
    ithread = 0 ; nthreads = 1
    !$ ithread = OMP_GET_THREAD_NUM () ; nthreads = OMP_GET_NUM_THREADS ()
    epoch = m_walltime()
    carrier%flop = 0
    t_calc_step = 0.0_dp
    t_process_stack = 0.0_dp
    !
    local_cols = right%local_indexing
    local_rows = left%local_indexing
    local_indexing = local_rows
    !
    ! Setup the carrier data that changes in each multiplication step.
    !$OMP CRITICAL (crit_data)
    carrier%a_data = left%data_area
    CALL dbcsr_data_hold (carrier%a_data)
    carrier%b_data = right%data_area
    CALL dbcsr_data_hold (carrier%b_data)
    !$OMP END CRITICAL (crit_data)
    !
    ! Find out the local A columns / B rows and sizes
    ! The right%local_rows is setup by the communication engine.
    k_locals => array_data (right%local_rows)
    carrier%k_locals => k_locals
    CALL ensure_array_size (carrier%k_sizes, ub=SIZE(k_locals), error=error)
    CALL local_filter(array_data(right%row_blk_size), SIZE(k_locals),&
         k_locals, carrier%k_sizes)    
    ! Setup the block norms
    carrier%a_norms => a_norms
    carrier%b_norms => b_norms
    !
    ! Start local multiplication
    IF (.TRUE.) THEN
       IF (dbg) THEN
          WRITE(*,*)"Calling sparse_multrec", left%nblkrows_local,&
               right%nblkcols_local, SIZE(k_locals), left%nblks, right%nblks
       ENDIF
       t_a_f = 1
       t_a_l = left%nblks
       t_b_f = 1
       t_b_l = right%nblks
       !$ ithread = OMP_GET_THREAD_NUM()
       !$ t_a_f = left%thr_c(ithread+1)+1
       !$ t_a_l = left%thr_c(ithread+2)
!write(*,*)routineN, ithread, t_a_f, t_a_l, "/", left%nblks
       CALL sparse_multrec(&
            1, left%nblkrows_local,&
            1, right%nblkcols_local,&
            1, SIZE(k_locals),&
            t_a_f, t_a_l, left%row_p,&
            t_b_f, t_b_l, right%row_p,&
            carrier, error, 0)
    ELSE
       CALL csr_multiply_unwrap(&
            1, left%nblkrows_local,&
            1, right%nblkrows_local,&
            1, SIZE(k_locals),&
            1, left%nblks, left%row_p,&
            1, right%nblks, right%row_p,&
            carrier, error)
    ENDIF
    !
    ! Flush the MM stack
    IF (carrier%mm_stack_p .GE. 1) THEN
       CALL process_mm_stack(carrier%mult_params,&
            carrier%mm_stack_p,&
            carrier%a_data, carrier%b_data,&
            carrier%c_data,&
            carrier%zero_first, carrier%zero_last,&
            carrier%lastblk,&
            error=error)
    ENDIF
    !
    flop = carrier%flop
    !
    !$OMP CRITICAL (crit_data)
    CALL dbcsr_data_release (carrier%a_data)
    CALL dbcsr_data_release (carrier%b_data)
    !$OMP END CRITICAL (crit_data)
    !
    IF (ASSOCIATED (carrier%k_sizes)) DEALLOCATE (carrier%k_sizes)
    !
    product%wms(ithread+1)%lastblk = carrier%lastblk
    product%wms(ithread+1)%datasize = carrier%datasize
    !
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_multrec_sparse


  SUBROUTINE setup_rec_index (matrix_set, error)
    TYPE(dbcsr_2d_array_type), INTENT(INOUT) :: matrix_set
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_rec_index', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handler, i_col, i_row, &
                                                ithread, t_f, t_l, t_size
    LOGICAL                                  :: thread_redist

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    DO i_row = 1, SIZE (matrix_set%mats, 1)
       DO i_col = 1, SIZE (matrix_set%mats, 2)
          IF (.FALSE.) &
               CALL dbcsr_reset_vlocals (matrix_set%mats(i_row, i_col),&
               matrix_set%image_dist, error=error)
          IF (dbg) THEN
             WRITE(*,*)routineN//" m, n, size",&
                  SIZE(matrix_set%mats(i_row, i_col)%m%row_p),&
                  dbcsr_nblkrows_local(matrix_set%mats(i_row, i_col)),&
                  dbcsr_nblkcols_local(matrix_set%mats(i_row, i_col))
             WRITE(*,'(3(1X,I7))')matrix_set%mats(i_row, i_col)%m%row_p
          ENDIF
          IF (careful_mod) THEN
             CALL dbcsr_assert (SIZE(matrix_set%mats(i_row, i_col)%m%row_p),&
                  "EQ", matrix_set%mats(i_row, i_col)%m%nblks*3,&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Block count mismatch.", __LINE__, error=error)
          ENDIF
          thread_redist = ASSOCIATED (matrix_set%mats(i_row, i_col)%m%thr_c)
          t_size = SIZE(matrix_set%mats(i_row, i_col)%m%row_p)/3
          t_f = 1
          t_l = t_size
          !$OMP PARALLEL IF (thread_redist) DEFAULT (none) &
          !$OMP PRIVATE (ithread) &
          !$OMP FIRSTPRIVATE (t_f, t_l, t_size) &
          !$OMP SHARED (matrix_set, i_row, i_col, thread_redist, error)
          !$ ithread = OMP_GET_THREAD_NUM()
          !$ IF (thread_redist) THEN
          !$    t_f = matrix_set%mats(i_row, i_col)%m%thr_c(ithread+1)+1
          !$    t_l = matrix_set%mats(i_row, i_col)%m%thr_c(ithread+2)
          !$ ENDIF
          t_size =  t_l - t_f + 1
          !$OMP BARRIER
          IF (t_size .GT. 0) THEN
             CALL call_rec_sort_index (&
                  dbcsr_nblkrows_local(matrix_set%mats(i_row, i_col)),&
                  dbcsr_nblkcols_local(matrix_set%mats(i_row, i_col)),&
                  t_size,&
                  matrix_set%mats(i_row, i_col)%m%row_p((t_f*3-2):(t_l*3)), error=error)
          ENDIF
          !$OMP END PARALLEL
       ENDDO
    ENDDO
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE setup_rec_index

! *****************************************************************************
!> \brief Sets up recursive multiplication
!>
!>
!> \param[in] left, right     left and right DBCSR matrices
!> \param[in,out] product     resulting DBCSR product matrix
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix, default is no
!> \param[in] eps             (optional) on-the-fly filtering epsilon
! *****************************************************************************
  SUBROUTINE dbcsr_multrec_init(left, right, product, carrier,&
       keep_sparsity, eps, row_max_epss, error)
    TYPE(dbcsr_type), INTENT(IN)             :: left, right
    TYPE(dbcsr_type), INTENT(INOUT)          :: product
    TYPE(carrier_type), INTENT(out)          :: carrier
    LOGICAL, INTENT(IN)                      :: keep_sparsity
    REAL(kind=real_8), INTENT(in), OPTIONAL  :: eps
    REAL(kind=sp), DIMENSION(:), INTENT(IN)  :: row_max_epss
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multrec_init', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: block_estimate, &
                                                c_nblkcols_local, &
                                                c_nblkrows_local, &
                                                error_handler, ithread, &
                                                nthreads
    INTEGER, DIMENSION(:), POINTER           :: c_local_cols, c_local_rows, &
                                                product_thread_dist
    LOGICAL                                  :: local_cols, local_indexing, &
                                                local_rows

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    IF (dbg) THEN
       multrec_calls = multrec_calls + 1
       WRITE(*,*)routineN//" multrec calls=", multrec_calls
    ENDIF
    !
    ithread = 0 ; nthreads = 1
    !$ ithread = OMP_GET_THREAD_NUM () ; nthreads = OMP_GET_NUM_THREADS ()
    carrier%id = ithread
    !
    ! Ensures that the index is correctly defined.
    CALL dbcsr_assert (left%list_indexing,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Must use list indexing for this routine.", __LINE__, error=error)
    !
    CALL dbcsr_assert ("NOT", left%bcsc,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Wrong routine for BCSC matrices.", __LINE__, error=error)
    CALL dbcsr_assert ("NOT", right%bcsc,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Wrong routine for BCSC matrices.", __LINE__, error=error)
    local_cols = right%local_indexing
    local_rows = left%local_indexing
    CALL dbcsr_assert (local_cols, "EQV", local_rows,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Local index useage must be consistent.", __LINE__, error=error)
    local_indexing = local_rows
    IF (local_cols) THEN
       CALL dbcsr_assert (left%local_indexing,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Wrong left format for local_cols.", __LINE__, error=error)
       CALL dbcsr_assert (right%local_indexing,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Wrong right format for local_cols.", __LINE__, error=error)
    ELSE
       CALL dbcsr_assert ("NOT",left%local_indexing,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Wrong left format for not local_cols.", __LINE__, error=error)
       CALL dbcsr_assert ("NOT",right%local_indexing,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Wrong right format for not local_cols.", __LINE__, error=error)
    ENDIF
    CALL dbcsr_assert (local_indexing,&
            dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
            "Global indexing not implemented.", __LINE__, error=error)
    !
    ! Fill carrier data structure.
    carrier%local_indexing = local_indexing
    carrier%keep_sparsity = keep_sparsity
    carrier%c_has_symmetry = product%symmetry
    carrier%use_eps = PRESENT (eps)
    carrier%my_wm = product%wms(ithread+1)
    carrier%lastblk = product%wms(ithread+1)%lastblk
    carrier%original_lastblk = carrier%lastblk
    carrier%datasize = product%wms(ithread+1)%datasize
    carrier%flop = INT(0, int_8)
    carrier%t_index = 0.0_dp
    carrier%t_gemm = 0.0_dp
    IF (PRESENT (eps)) THEN
       carrier%eps = eps
    ELSE
       carrier%eps = 0.0_real_8
    ENDIF
    !
    carrier%c_data = product%wms(ithread+1)%data_area
    CALL dbcsr_data_hold (carrier%c_data)
    !
    !
    !$ NULLIFY (product_thread_dist)    
    !$ CALL dbcsr_assert (dbcsr_distribution_has_threads (product%dist),&
    !$      dbcsr_fatal_level, dbcsr_internal_error, routineN,&
    !$      "Missing thread distribution.", __LINE__, error=error)
    !$ product_thread_dist => array_data (&
    !$      dbcsr_distribution_thread_dist (product%dist))
    !
    ! Find out the C/A rows and C/B columns and sizes.
    c_nblkrows_local = product%nblkrows_local
    c_local_rows => array_data (product%local_rows)
    c_nblkcols_local = product%nblkcols_local
    c_local_cols => array_data (product%local_cols)
    carrier%c_local_rows => c_local_rows
    carrier%c_local_cols => c_local_cols
    IF (dbg) WRITE(*,*)"setting up for product", product%name
    IF (careful_mod) THEN
       IF (.NOT. array_equality (dbcsr_distribution_local_rows (product%dist),&
                                 product%local_rows)) THEN
          WRITE(*,*)"row dist", array_data(dbcsr_distribution_row_dist(product%dist))
          WRITE(*,*)"dist local rows", array_data(dbcsr_distribution_local_rows (product%dist))
          WRITE(*,*)" mat local rows", array_data(product%local_rows)
          CALL dbcsr_assert (.FALSE., &
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Array mismatch.", __LINE__, error=error)
       ENDIF
       IF (.NOT. array_equality (dbcsr_distribution_local_cols (product%dist),&
                                 product%local_cols)) THEN
          WRITE(*,*)"col dist", array_data(dbcsr_distribution_col_dist(product%dist))
          WRITE(*,*)"dist local cols", array_data(dbcsr_distribution_local_cols (product%dist))
          WRITE(*,*)" mat local cols", array_data(product%local_cols)
          CALL dbcsr_assert (.FALSE., &
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Array mismatch.", __LINE__, error=error)
       ENDIF
       CALL dbcsr_assert (SIZE(c_local_rows), "EQ", c_nblkrows_local,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Row count mismatch.", __LINE__, error=error)
       CALL dbcsr_assert (SIZE(c_local_cols), "EQ", c_nblkcols_local,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Column count mismatch.", __LINE__, error=error)
    ENDIF
    !
    ! And the k epsilons
    ALLOCATE (carrier%row_max_epss(c_nblkrows_local))
    IF (carrier%use_eps) THEN
       CALL local_filter_sp(row_max_epss, c_nblkrows_local, c_local_rows,&
            carrier%row_max_epss)
    ELSE
       carrier%row_max_epss(:) = -HUGE(0.0_sp)
    ENDIF
    !
    ALLOCATE(carrier%m_sizes(c_nblkrows_local))
    CALL local_filter(array_data (product%row_blk_size), SIZE(c_local_rows),&
         c_local_rows, carrier%m_sizes)
    ALLOCATE(carrier%n_sizes(c_nblkcols_local))
    CALL local_filter(array_data (product%col_blk_size), SIZE(c_local_cols),&
         c_local_cols, carrier%n_sizes)
    !
    NULLIFY (carrier%k_locals)
    NULLIFY (carrier%k_sizes)
    !
    ! Setup the hash tables if needed
    ALLOCATE (carrier%c_hashes (product%nblkrows_local))
    block_estimate=MAX(product%nblks,left%nblks,right%nblks)/nthreads
    IF (local_indexing) THEN
       CALL fill_hash_tables (carrier%c_hashes, product,block_estimate,&
            row_map=array_data(product%global_rows),&
            col_map=array_data(product%global_cols),&
            error=error)
    ELSE
       CALL fill_hash_tables (carrier%c_hashes, product,block_estimate,&
            error=error)
    ENDIF
    !
    ! Setup the MM stack
    ALLOCATE (carrier%mult_params(n_mult_params, mm_stack_size))
    carrier%mm_stack_p = 0
    carrier%zero_first = carrier%datasize+1
    carrier%zero_last = carrier%zero_first-1
    !
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_multrec_init


! *****************************************************************************
!> \brief Sets up recursive multiplication
!>
!>
!> \param[in] left, right     left and right DBCSR matrices
!> \param[in,out] product     resulting DBCSR product matrix
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix, default is no
!> \param[in] eps             (optional) on-the-fly filtering epsilon
! *****************************************************************************
  SUBROUTINE dbcsr_multrec_finalize(carrier, product, error)
    TYPE(carrier_type), INTENT(inout)        :: carrier
    TYPE(dbcsr_type), INTENT(inout)          :: product
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multrec_finalize', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handler, i, ithread

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM()

    IF (dbg) THEN
       multrec_calls = multrec_calls - 1
       WRITE(*,*)routineN//" multrec calls=", multrec_calls
    ENDIF
    ! Release the carrier
    DEALLOCATE (carrier%m_sizes)
    DEALLOCATE (carrier%n_sizes)
    DEALLOCATE (carrier%row_max_epss)
    ! Clear hash tables
    DO i = 1, SIZE(carrier%c_hashes)
       CALL hash_table_release (carrier%c_hashes (i))
    ENDDO
    DEALLOCATE (carrier%c_hashes)
    ! Clear MM stack
    DEALLOCATE (carrier%mult_params)
    CALL dbcsr_data_release (carrier%c_data)
    !
    IF (carrier%local_indexing) THEN
       CALL remap_local2global(carrier%my_wm%row_i, carrier%my_wm%col_i,&
            carrier%c_local_rows, carrier%c_local_cols,&
            carrier%original_lastblk+1, carrier%lastblk)
    ENDIF
    ! Reinstate WM
    carrier%my_wm%lastblk = carrier%lastblk
    carrier%my_wm%datasize = carrier%datasize
    product%wms(ithread+1) = carrier%my_wm
    !
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_multrec_finalize

  !> \brief Packs a globally-indexed array into a locally-indexed array.
  PURE SUBROUTINE remap_local2global(row_i, col_i, local_rows, local_cols,&
       first, last)
    INTEGER, INTENT(in)                      :: last, first
    INTEGER, DIMENSION(:), INTENT(in)        :: local_cols, local_rows
    INTEGER, DIMENSION(1:last), &
      INTENT(inout)                          :: col_i, row_i

    INTEGER                                  :: i

    FORALL (i = first : last)
       row_i(i) = local_rows(row_i(i))
       col_i(i) = local_cols(col_i(i))
    END FORALL
  END SUBROUTINE remap_local2global

  !> \brief Maps between locally-indexed arrays with different local indexing.
  PURE SUBROUTINE remap_local2local(n_local_src, local_remap, local_src, &
       n_globals, global_dst)
    INTEGER, INTENT(IN)                      :: n_local_src
    INTEGER, DIMENSION(1:n_local_src), &
      INTENT(OUT)                            :: local_remap
    INTEGER, DIMENSION(1:n_local_src), &
      INTENT(IN)                             :: local_src
    INTEGER, INTENT(IN)                      :: n_globals
    INTEGER, DIMENSION(1:n_globals), &
      INTENT(IN)                             :: global_dst

    INTEGER                                  :: i

    FORALL (i = 1 : n_local_src)
       local_remap(i) = global_dst (local_src (i))
    END FORALL
  END SUBROUTINE remap_local2local



  !> \author JV
  PURE FUNCTION find_cut_row(ai,af,a,val) RESULT(res)
    INTEGER, INTENT(IN)                      :: ai, af
    INTEGER, DIMENSION(3, 1:af), INTENT(IN)  :: a
    INTEGER, INTENT(IN)                      :: val
    INTEGER                                  :: res

    INTEGER                                  :: i, ihigh, ilow

! do a log(N) search along the ordered index

    ilow = ai
    IF (a(1,ilow) > val)  THEN
       res = ilow
       RETURN
    ENDIF

    ihigh = af
    IF (a(1,ihigh) <= val)  THEN
       res = ihigh+1
       RETURN
    ENDIF

    DO
       IF (ihigh-ilow == 1) EXIT
       i = (ilow + ihigh)/2
       IF (a(1,i)>val) THEN
          ihigh=i
       ELSE
          ilow=i
       ENDIF
    ENDDO
    res=ihigh

    ! the linear search version
    ! DO i=ai,af
    !    IF (a(i)%r>val) EXIT
    !ENDDO
    !res=i
  END FUNCTION find_cut_row

  !> \author JV
  PURE FUNCTION find_cut_col(ai,af,a,val) RESULT(res)
    INTEGER, INTENT(IN)                      :: ai, af
    INTEGER, DIMENSION(3, 1:af), INTENT(IN)  :: a
    INTEGER, INTENT(IN)                      :: val
    INTEGER                                  :: res

    INTEGER                                  :: i, ihigh, ilow

! do a log(N) search along the ordered index

    ilow = ai
    IF (a(2,ilow) > val)  THEN
       res = ilow
       RETURN
    ENDIF

    ihigh = af
    IF (a(2,ihigh) <= val)  THEN
       res = ihigh+1
       RETURN
    ENDIF

    DO
       IF (ihigh-ilow == 1) EXIT
       i = (ilow + ihigh)/2
       IF (a(2,i) > val) THEN
          ihigh = i
       ELSE
          ilow = i
       ENDIF
    ENDDO
    res = ihigh

    ! the linear search version
    ! DO i=ai,af
    !    IF (a(i)%c>val) EXIT
    !ENDDO
    !res=i
  END FUNCTION find_cut_col

! *****************************************************************************
!> \brief Performs recursive multiplication
!> \author Joost VandeVondele
! *****************************************************************************
  RECURSIVE SUBROUTINE sparse_multrec(mi, mf, ni, nf, ki, kf,&
       ai, af, a_index, bi, bf, b_index, &
       carrier, error, d)
    INTEGER, INTENT(IN)                      :: mi, mf, ni, nf, ki, kf, ai, af
    INTEGER, DIMENSION(3, 1:af), INTENT(IN)  :: a_index
    INTEGER, INTENT(IN)                      :: bi, bf
    INTEGER, DIMENSION(3, 1:bf), INTENT(IN)  :: b_index
    TYPE(carrier_type), INTENT(INOUT)        :: carrier
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER, INTENT(IN)                      :: d

    CHARACTER(len=*), PARAMETER :: routineN = 'sparse_multrec', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: norec = 512
    LOGICAL, PARAMETER                       :: careful = careful_mod, &
                                                dbg = .FALSE.

    INTEGER                                  :: acut, bcut, cut, K, M, N, s1

!   ---------------------------------------------------------------------------

    IF (dbg) THEN
       WRITE(*,'(I7,1X,5(A,2(1X,I7)))')d," rm", mi, mf,",",ni,nf,",",ki,kf,"/",ai,af,",",bi,bf
    ENDIF
    IF (.TRUE.) THEN
       IF (af .LT. ai .OR. bf .LT. bi .OR. mf .LT. mi .OR. nf .LT. ni .OR. kf .LT. ki) THEN
          IF (dbg) WRITE(*,*)"Empty"
          RETURN
       ENDIF
    ENDIF
    IF (af-ai+1 <= norec .AND. bf-bi+1 <= norec) THEN
       CALL csr_multiply_unwrap(&
            mi, mf, ni, nf, ki, kf,&
            ai, af, a_index,&
            bi, bf, b_index,&
            carrier, error)
       RETURN
    ENDIF

    M = mf-mi + 1
    N = nf-ni + 1
    K = kf-ki + 1
    IF (dbg) THEN
       WRITE(*,*)'m,k,n',M,K,N
    ENDIF
    IF (M >= MAX(N, K)) cut = 1
    IF (K >= MAX(N, M)) cut = 2
    IF (N >= MAX(M, K)) cut = 3
    SELECT CASE(cut)
    CASE(1)
       s1=M/2
       acut = find_cut_row(ai,af,a_index,mi+s1-1)
       CALL sparse_multrec(mi,mi+s1-1, ni,nf, ki,kf,&
            ai,acut-1,a_index, bi,bf,b_index, carrier, error,d+1)
       CALL sparse_multrec(mi+s1,mf, ni,nf, ki,kf,&
            acut,af,a_index, bi,bf,b_index, carrier, error,d+1)
    CASE(2)
       s1=K/2
       acut = find_cut_col(ai,af,a_index,ki+s1-1)
       IF (dbg) THEN
          WRITE(*,*)N,s1,ni+s1-1,"/",ai,af,acut
          WRITE(*,'(3(I7))')a_index
       ENDIF
       bcut = find_cut_row(bi,bf,b_index,ki+s1-1)
       IF (dbg) THEN
          WRITE(*,*)N,s1,ni+s1-1,"/",bi,bf,bcut
          WRITE(*,'(3(I7))')b_index
       ENDIF
       CALL sparse_multrec(mi,mf, ni,nf, ki,ki+s1-1,&
            ai,acut-1,a_index, bi,bcut-1,b_index, carrier, error,d+1)
       CALL sparse_multrec(mi,mf, ni,nf, ki+s1,kf,&
            acut,af,a_index, bcut,bf,b_index, carrier, error,d+1)
    CASE(3)
       s1=N/2
       bcut = find_cut_col(bi,bf,b_index,ni+s1-1)
       IF (dbg) THEN
          WRITE(*,*)N,s1,ni+s1-1,"/",bi,bf,bcut
          WRITE(*,'(3(I7))')b_index
       ENDIF
       CALL sparse_multrec(mi,mf, ni,ni+s1-1, ki,kf,&
            ai,af,a_index, bi,bcut-1,b_index, carrier, error,d+1)
       CALL sparse_multrec(mi,mf, ni+s1,nf, ki,kf,&
            ai,af,a_index, bcut,bf,b_index, carrier, error,d+1)
    END SELECT
  END SUBROUTINE sparse_multrec

! *****************************************************************************
!> \brief Performs multiplication of smaller submatrices.
!>
!> This routine is used for setting up and calling csr_multiply with
!> as few argument indirections as possible.
! *****************************************************************************
  SUBROUTINE csr_multiply_unwrap(mi, mf, ni, nf, ki, kf,&
       ai, af, a_index, bi, bf, b_index, &
       carrier, error)
    INTEGER, INTENT(IN)                      :: mi, mf, ni, nf, ki, kf, ai, af
    INTEGER, DIMENSION(1:3, 1:af), &
      INTENT(IN)                             :: a_index
    INTEGER, INTENT(IN)                      :: bi, bf
    INTEGER, DIMENSION(1:3, 1:bf), &
      INTENT(IN)                             :: b_index
    TYPE(carrier_type), INTENT(INOUT)        :: carrier
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'csr_multiply_unwrap', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: max_new_nblks, n_a_norms, &
                                                n_b_norms, nblks_new
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: a_row_p, b_row_p
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: a_blk_info, b_blk_info
    REAL(KIND=sp), DIMENSION(1:af-ai+1)      :: csr_a_norms
    REAL(KIND=sp), DIMENSION(1:bf-bi+1)      :: csr_b_norms

!   ---------------------------------------------------------------------------

    IF (dbg) THEN
       WRITE(*,'(I7,1X,5(A,2(1X,I7)))')0,"uwr", mi, mf,",",ni,nf,",",ki,kf,"/",ai,af,",",bi,bf
    ENDIF

    IF (af-ai+1 .GT. 0 .AND. bf-bi+1 .GT. 0) THEN
      
       ! the maximum number of blocks can be safely estimated by considering both the rowxcol,
       ! but also the blocks the latter can never be larger than norec**2, which is a 'small' constant
       max_new_nblks = MIN(INT(mf-mi+1,int_8) * INT(nf-ni+1,int_8), &
                           INT(af-ai+1,int_8) * INT(bf-bi+1,int_8))
       nblks_new = carrier%lastblk + max_new_nblks
       
       CALL ensure_array_size(carrier%my_wm%row_i, ub=nblks_new,&
            factor=default_resize_factor, error=error)
       CALL ensure_array_size(carrier%my_wm%col_i, ub=nblks_new,&
            factor=default_resize_factor, error=error)
       CALL ensure_array_size(carrier%my_wm%blk_p, ub=nblks_new,&
            factor=default_resize_factor, error=error)

       ALLOCATE (a_row_p(mi:mf+1))
       ALLOCATE (b_row_p(ki:kf+1))
       ALLOCATE (a_blk_info(2,af-ai+1))
       ALLOCATE (b_blk_info(2,bf-bi+1))
       !
       IF (carrier%use_eps) THEN
          n_a_norms = af-ai+1
          n_b_norms = bf-bi+1
       ELSE
          n_a_norms = 0
          n_b_norms = 0
       ENDIF
       !
       ! Build the indices
       CALL build_csr_index (mi,mf,ai,af,a_row_p, a_blk_info, a_index,&
            n_a_norms, csr_a_norms, carrier%a_norms)
       CALL build_csr_index (ki,kf,bi,bf,b_row_p, b_blk_info, b_index,&
            n_b_norms, csr_b_norms, carrier%b_norms)
       !
       CALL csr_multiply (mi, mf, ni, nf, ki, kf,&
            ai, af, a_row_p, a_blk_info, bi, bf, b_row_p, b_blk_info,&
            carrier%c_hashes,&
            carrier%my_wm%row_i, carrier%my_wm%col_i, carrier%my_wm%blk_p,&
            carrier%lastblk, carrier%datasize,&
            carrier%m_sizes, carrier%n_sizes, carrier%k_sizes,&
            carrier%c_local_rows, carrier%c_local_cols,&
            carrier%c_has_symmetry, carrier%keep_sparsity, carrier%use_eps,&
            carrier%mult_params, carrier%mm_stack_p,&
            carrier%zero_first, carrier%zero_last,&
            carrier%row_max_epss, csr_a_norms, csr_b_norms,&
            carrier%a_data, carrier%b_data, carrier%c_data,&
            carrier%flop, carrier%t_gemm, error)
       !
       DEALLOCATE (a_row_p)
       DEALLOCATE (b_row_p)
       DEALLOCATE (a_blk_info)
       DEALLOCATE (b_blk_info)
    ENDIF
  END SUBROUTINE csr_multiply_unwrap

! *****************************************************************************
!> \brief Performs multiplication of smaller submatrices.
! *****************************************************************************
  SUBROUTINE csr_multiply(mi, mf, ni, nf, ki, kf,&
       ai, af, a_row_p, a_blk_info, bi, bf, b_row_p, b_blk_info,&
       c_hashes, c_row_i, c_col_i, c_blk_p, lastblk, datasize,&
       m_sizes, n_sizes, k_sizes,&
       c_local_rows, c_local_cols,&
       c_has_symmetry, keep_sparsity, use_eps,&
       mult_params, stack_p, zero_first, zero_last,&
       row_max_epss, left_norms, right_norms,&
       a_data, b_data, c_data, flop, t_gemm, error)
    INTEGER, INTENT(IN)                      :: mi, mf, ni, nf, ki, kf, ai, af
    INTEGER, DIMENSION(mi:mf+1), INTENT(IN)  :: a_row_p
    INTEGER, DIMENSION(2, 1:af-ai+1), &
      INTENT(IN)                             :: a_blk_info
    INTEGER, INTENT(IN)                      :: bi, bf
    INTEGER, DIMENSION(ki:kf+1), INTENT(IN)  :: b_row_p
    INTEGER, DIMENSION(2, 1:bf-bi+1), &
      INTENT(IN)                             :: b_blk_info
    TYPE(hash_table_type), DIMENSION(:), &
      INTENT(INOUT)                          :: c_hashes
    INTEGER, DIMENSION(*), INTENT(INOUT)     :: c_row_i, c_col_i, c_blk_p
    INTEGER, INTENT(INOUT)                   :: lastblk, datasize
    INTEGER, DIMENSION(*), INTENT(IN)        :: m_sizes, n_sizes, k_sizes, &
                                                c_local_rows, c_local_cols
    LOGICAL, INTENT(IN)                      :: c_has_symmetry, &
                                                keep_sparsity, use_eps
    INTEGER, DIMENSION(n_mult_params, &
      mm_stack_size), INTENT(INOUT)          :: mult_params
    INTEGER, INTENT(INOUT)                   :: stack_p, zero_first, zero_last
    REAL(kind=sp), DIMENSION(*)              :: row_max_epss, left_norms, &
                                                right_norms
    TYPE(dbcsr_data_obj), INTENT(IN)         :: a_data, b_data
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: c_data
    INTEGER(KIND=int_8), INTENT(INOUT)       :: flop
    REAL(KIND=dp), INTENT(INOUT)             :: t_gemm
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'csr_multiply', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE., &
                                                local_timing = .FALSE.

    INTEGER :: a_blk, a_col_l, a_row_l, b_blk, b_col_l, c_blk_id, c_blk_pt, &
      c_col_logical, c_nze, c_row_logical, ithread, k_size, m_size, n_size, &
      new_blk
    INTEGER(KIND=int_4)                      :: offset
    LOGICAL                                  :: block_exists, flush_stack
    REAL(kind=dp)                            :: t_gemm_me
    REAL(kind=sp)                            :: a_norm, a_row_eps, b_norm

!   ---------------------------------------------------------------------------

    ithread = 0
    !$ ithread = omp_get_thread_num()
    t_gemm_me = REAL(0, KIND(t_gemm_me))
    !
    ! New data blocks will be put into the data area starting at
    ! the c_blk_pt position.
    new_blk = lastblk
    c_blk_pt = datasize + 1
    !
    a_row_cycle: DO a_row_l = mi, mf
       m_size = m_sizes(a_row_l)

       a_row_eps = row_max_epss (a_row_l)
       
       a_blk_cycle: DO a_blk = a_row_p(a_row_l)+1, a_row_p(a_row_l+1)
          a_col_l = a_blk_info(1, a_blk)
          IF (debug_mod) WRITE(*,*)ithread,routineN//" A col", a_col_l,";",a_row_l
          k_size = k_sizes (a_col_l)
          a_norm = left_norms(a_blk)
          b_blk_cycle: DO b_blk = b_row_p(a_col_l)+1, b_row_p(a_col_l+1)
             IF (dbg) THEN
                WRITE(*,'(1X,A,3(1X,I7),1X,A,1X,I16)')routineN//" trying B",&
                     a_row_l, b_blk_info(1,b_blk), a_col_l, "at", b_blk_info(2,b_blk)
             ENDIF
             b_norm = right_norms(b_blk)
             IF (a_norm * b_norm .LT. a_row_eps) THEN
                CYCLE
             ENDIF
             b_col_l = b_blk_info(1,b_blk)
             ! Don't calculate symmetric blocks.
             symmetric_product: IF (c_has_symmetry) THEN
                c_row_logical = c_local_rows (a_row_l)
                c_col_logical = c_local_cols (b_col_l)
                IF (c_row_logical .NE. c_col_logical&
                     .AND. my_checker_tr (c_row_logical, c_col_logical)) THEN
                   IF (dbg) THEN
                      WRITE(*,*)"Skipping symmetric block!", c_row_logical,&
                           c_col_logical
                   ENDIF
                   CYCLE
                ENDIF
             ENDIF symmetric_product

             n_size = n_sizes(b_col_l)
             c_nze = m_size * n_size
             !
             c_blk_id = hash_table_get (c_hashes(a_row_l), b_col_l)
             IF (.FALSE.) THEN
                WRITE(*,'(1X,A,3(1X,I7),1X,A,1X,I16)')routineN//" coor",&
                     a_row_l, a_col_l, b_col_l,"c blk", c_blk_id
             ENDIF
             block_exists = c_blk_id .GT. 0

             new_block_case: IF (.NOT. block_exists) THEN
                sparsity_enforcement: IF (keep_sparsity) THEN
                   CYCLE
                ENDIF sparsity_enforcement
                offset = c_blk_pt
                new_blk = new_blk+1
                c_blk_id = new_blk
                IF (dbg) WRITE(*,*)routineN//" new block offset, nze", offset, c_nze
                CALL hash_table_add(c_hashes(a_row_l),&
                     b_col_l, c_blk_id, error=error)
                !
                ! We still keep the linear index because it's
                ! easier than getting the values out of the
                ! hashtable in the end.
                c_row_i(new_blk) = a_row_l
                c_col_i(new_blk) = b_col_l
                c_blk_p(new_blk) = offset
                !
                c_blk_pt = c_blk_pt + c_nze
                datasize = datasize + c_nze
                zero_last = c_blk_pt - 1
             ELSE
                offset = c_blk_p(c_blk_id)
             ENDIF new_block_case
             !
             ! We should not call certain MM routines (netlib BLAS)
             ! with zero LDs; however, we still need to get to here
             ! to get new blocks.
             IF (careful_mod) THEN
                IF (c_nze .EQ. 0 .OR. k_size .EQ. 0) THEN
                   CALL dbcsr_assert (.FALSE.,&
                        dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                        "Can not call MM with LDx=0.", __LINE__, error=error)
                   CYCLE
                ENDIF
             ENDIF
             !
             stack_p = stack_p + 1
             !
             mult_params(p_a_first,stack_p) = a_blk_info(2, a_blk)
             !
             mult_params(p_b_first,stack_p) = b_blk_info(2, b_blk)
             !
             mult_params(p_c_first,stack_p) = offset
             mult_params(p_c_blk,stack_p) = c_blk_id
             !
             mult_params(p_m,stack_p) = m_size
             mult_params(p_n,stack_p) = n_size
             mult_params(p_k,stack_p) = k_size
             !
             flop = flop + INT(2, int_8) * &
                  INT(m_size, int_8) * &
                  INT(n_size, int_8) * &
                  INT(k_size, int_8)
             !
             flush_stack = (stack_p .GE. mm_stack_size-1)
             IF (flush_stack) THEN
                IF (local_timing) THEN
                   !$ IF (.TRUE.) THEN
                   !$    t_gemm_me = t_gemm_me - OMP_GET_WTIME()
                   !$ ELSE
                   t_gemm_me = t_gemm_me - m_walltime()
                   !$ ENDIF
                ENDIF
                CALL process_mm_stack(mult_params,&
                     stack_p,&
                     a_data, b_data,&
                     c_data,&
                     zero_first, zero_last,&
                     new_blk,&
                     error=error)
                IF (local_timing) THEN
                   !$ IF (.TRUE.) THEN
                   !$    t_gemm_me = t_gemm_me + OMP_GET_WTIME()
                   !$ ELSE
                   t_gemm_me = t_gemm_me + m_walltime()
                   !$ ENDIF
                ENDIF
                stack_p = 0
             ENDIF
          ENDDO b_blk_cycle ! b
       ENDDO a_blk_cycle ! a_col
    ENDDO a_row_cycle ! a_row
    !
    lastblk = new_blk
    datasize = c_blk_pt -1
    !
    t_gemm = t_gemm + t_gemm_me
    !
  END SUBROUTINE csr_multiply

! *****************************************************************************
!> \brief  Builds and sorts a CSR index from a list index.
!> \author JV
!> <b>Modification history:</b>
!> - 2011-02-15 [UB] Adapted to use DBCSR-type CSR indexing
! *****************************************************************************
!  PURE SUBROUTINE build_csr_index(mi,mf,ai,af, row_p, blk_info, list_index)
  SUBROUTINE build_csr_index(mi,mf,ai,af, row_p, blk_info, list_index,&
       nnorms, csr_norms, list_norms)
    INTEGER, INTENT(IN)                      :: mi, mf, ai, af
    INTEGER, DIMENSION(mi:mf+1), INTENT(OUT) :: row_p
    INTEGER, DIMENSION(2, 1:af-ai+1), &
      INTENT(OUT)                            :: blk_info
    INTEGER, DIMENSION(3, 1:af), INTENT(IN)  :: list_index
    INTEGER, INTENT(IN)                      :: nnorms
    REAL(KIND=sp), DIMENSION(1:af-ai+1), &
      INTENT(OUT)                            :: csr_norms
    REAL(KIND=sp), DIMENSION(:), INTENT(IN)  :: list_norms

    CHARACTER(len=*), PARAMETER :: routineN = 'build_csr_index', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: careful = .FALSE., &
                                                dbg = .FALSE.

    INTEGER                                  :: i, row
    INTEGER, DIMENSION(mi:mf)                :: counts
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------
! Counts blocks per row and calculates the offsets.

    IF (dbg) THEN
       WRITE(*,'(I7,1X,5(A,2(1X,I7)))')0,"bci", mi, mf,";",ai,af
       !write(*,'(3(I7))')list_index(:,ai:af)
    ENDIF

    counts(:) = 0
    DO i = ai, af
       IF (careful) THEN
          CALL dbcsr_assert (list_index(1,i), "GE", mi,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Out of range", __LINE__, error=error)
          CALL dbcsr_assert (list_index(1,i), "LE", mf,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Out of range", __LINE__, error=error)
       ENDIF
       counts(list_index(1,i)) = counts(list_index(1,i))+1
    ENDDO
    row_p(mi) = 0
    DO i = mi+1, mf+1
       row_p(i) = row_p(i-1) + counts(i-1)
    ENDDO
    ! Adds every block to its corresponding row.
    counts(:) = 0
    DO i = ai, af
       row = list_index(1,i)
       counts(row) = counts(row)+1
       IF (careful) THEN
          CALL dbcsr_assert (row_p(row) + counts(row), "LE", af-ai+1,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Out of range", __LINE__, error=error)
          CALL dbcsr_assert (row_p(row) + counts(row), "GE", 1,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Out of range", __LINE__, error=error)
       ENDIF
       blk_info(1, row_p(row) + counts(row)) = list_index(2,i)
       blk_info(2, row_p(row) + counts(row)) = list_index(3,i)
       IF (nnorms .GT. 0) THEN
          csr_norms(row_p(row) + counts(row)) = list_norms(i)
       ENDIF
    ENDDO
    IF (nnorms .EQ. 0) THEN
       csr_norms(:) = 0.0_sp
    ENDIF
  END SUBROUTINE build_csr_index


! *****************************************************************************
!> \brief Used to thunk a call to rec_sort_index
! *****************************************************************************
  SUBROUTINE call_rec_sort_index (m,n,nblks,idx, error)
    INTEGER, INTENT(IN)                      :: m, n, nblks
    INTEGER, DIMENSION(3, 1:nblks), &
      INTENT(INOUT)                          :: idx
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'call_rec_sort_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    IF (.FALSE.) WRITE(*,*)" Calling rec_sort_index, size", nblks
    CALL rec_sort_index(1, m, 1, n, nblks, idx, 0)
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE call_rec_sort_index


! *****************************************************************************
!> \brief Sorts index for recursing.
!> \author JV
!> \note Always cut longest first. On a tie cut N
!> \par History
!> - 2011-02-17 [UB] modified for use in DBCSR; reduced memory usage.
! *****************************************************************************
  RECURSIVE SUBROUTINE rec_sort_index(mi,mf,ni,nf,nele,a,d)
    INTEGER, INTENT(IN)                      :: mi, mf, ni, nf, nele
    INTEGER, DIMENSION(3, 1:nele), &
      INTENT(inout)                          :: a
    INTEGER, INTENT(IN)                      :: d

    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: half, M, N, nlow
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: tmp

!   ---------------------------------------------------------------------------

    IF (dbg) THEN
       WRITE(*,*)" rs", mi, mf,"/",ni,nf,"=>",nele, d
       !WRITE(*,'(3(1X,I7))')a(:,1:nele)
    ENDIF
    IF (dbg) THEN
       IF (d .GT. 20) THEN
          WRITE(*,*)a(1,-d*1000)
       ENDIF
    ENDIF
    IF (nele .LE. 1) RETURN
    ALLOCATE(tmp(3,nele))
    M = mf-mi+1
    N = nf-ni+1
    IF (M > N) THEN
       half = M/2
       CALL rec_split (nele, a, tmp, 1, nlow, mi, mf, half)
       a = tmp
       DEALLOCATE (tmp)
       CALL rec_sort_index(mi,mi+half-1,ni,nf, nlow, a(:,1:nlow), d+1)
       CALL rec_sort_index(mi+half,mf,ni,nf, nele-nlow, a(:,nlow+1:nele), d+1)
    ELSE
       half = N/2
       CALL rec_split (nele, a, tmp, 2, nlow, ni, nf, half)
       a = tmp
       DEALLOCATE (tmp)
       CALL rec_sort_index(mi,mf,ni,ni+half-1, nlow, a(:,1:nlow), d+1)
       CALL rec_sort_index(mi,mf,ni+half,nf, nele-nlow, a(:,nlow+1:nele), d+1)
    ENDIF
  END SUBROUTINE rec_sort_index


  SUBROUTINE rec_split (nele, a, split, row_or_col, nlow, mi, mf, half)
    INTEGER, INTENT(IN)                      :: nele
    INTEGER, DIMENSION(3, nele), INTENT(IN)  :: a
    INTEGER, DIMENSION(3, nele), INTENT(OUT) :: split
    INTEGER, INTENT(IN)                      :: row_or_col
    INTEGER, INTENT(OUT)                     :: nlow
    INTEGER, INTENT(IN)                      :: mi, mf, half

    INTEGER                                  :: el, half_m, p_high, p_low

    half_m = mi+half-1
    p_low = 1
    p_high = nele
    DO el = 1, nele
       IF (a(row_or_col,el) <= half_m) THEN
          split(1:3, p_low) = a(1:3, el)
          p_low = p_low + 1
       ELSE
          split(1:3, p_high) = a(1:3, el)
          p_high = p_high - 1
       ENDIF
    ENDDO
    nlow = p_low - 1
    IF (p_high .NE. nlow) STOP
  END SUBROUTINE rec_split

! *****************************************************************************
!> \brief Ensures a real variable is in the range [0, 1].
! *****************************************************************************
  ELEMENTAL SUBROUTINE crop_0_1(v)
    REAL(KIND=real_8), INTENT(INOUT)         :: v

    v = MIN (v, 1.0_real_8)
    v = MAX (v, 0.0_real_8)
  END SUBROUTINE crop_0_1


! *****************************************************************************
!> \brief Issues actual GEMM calls.
!>
!> \param[in] params           Stack of GEMM parameters
!> \param[in] n                Number of parameters
!> \param[in] left_data_area   Left-matrix data
!> \param[in] right_data_area  Right-matrix data
!> \param[in,out] zero_first   Zero product data area starting from this
!>                             element
!> \param[in,out] zero_last    Zero product data area up to this element
!> \param[in] lastblk          Number of blocks in product
!> \param[in,out] product_data_area  Data for results
! *****************************************************************************
  SUBROUTINE process_mm_stack(params,&
       stack_size, &
       left_data_area, right_data_area, product_data_area,&
       zero_first, zero_last, nblks,&
       error)
    INTEGER, INTENT(INOUT)                   :: stack_size
    INTEGER, &
      DIMENSION(n_mult_params, stack_size), &
      INTENT(IN)                             :: params
    TYPE(dbcsr_data_obj), INTENT(IN)         :: left_data_area, &
                                                right_data_area
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: product_data_area
    INTEGER, INTENT(INOUT)                   :: zero_first, zero_last
    INTEGER, INTENT(IN)                      :: nblks
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'process_mm_stack', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: careful = .FALSE., &
                                                remote_memory = .FALSE.
    REAL, PARAMETER                          :: resize_factor = 1.618034

    INTEGER :: istat, max_size_a, max_size_b, max_size_c, maxs, min_size_a, &
      min_size_b, min_size_c, sp, tmp_i
    INTEGER, DIMENSION(6)                    :: bounds

!   ---------------------------------------------------------------------------

    t_process_stack = t_process_stack - m_walltime()
    IF (remote_memory) THEN
       min_size_a = dbcsr_data_get_size(left_data_area)
       min_size_b = dbcsr_data_get_size(right_data_area)
       min_size_c = dbcsr_data_get_size(product_data_area)
       max_size_a = 1
       max_size_b = 1
       max_size_c = 1
       DO sp = 1, stack_size
          min_size_a = MIN (min_size_a, params(p_a_first,sp))
          max_size_a = MAX (max_size_a, params(p_a_first,sp) &
               + params(p_m,sp) * params(p_k,sp))
          min_size_b = MIN (min_size_b, params(p_b_first,sp))
          max_size_b = MAX (max_size_b, params(p_b_first,sp) &
               + params(p_k,sp) * params(p_n,sp))
          min_size_c = MIN (min_size_c, params(p_c_first,sp))
          max_size_c = MAX (max_size_c, params(p_c_first,sp) &
               + params(p_m,sp) * params(p_n,sp))
       ENDDO
       bounds = ((/ min_size_a, max_size_a, min_size_b, max_size_b, min_size_c, max_size_c /))
    ELSE
       bounds = 0
    ENDIF
    !
    ! Resize target data area if necessary.
    maxs = dbcsr_data_get_size(product_data_area)
    IF (zero_last .GT. maxs) THEN
       maxs = zero_last
       IF (debug_mod) &
            WRITE(*,*)routineN//" Resizing to", LOG(REAL(maxs))/LOG(10.0)
       CALL dbcsr_data_ensure_size (product_data_area,&
            maxs, factor=resize_factor, error=error)
       cuda_realloc_c: IF (mm_driver .EQ. mm_driver_cuda) THEN
          tmp_i = maxs
          CALL dbcsr_cuda_thread_sync(error=error)
          IF (verbose_acc) WRITE(*,*)routineN//" reallocating c_dev"
          CALL dbcsr_cuda_dev_mem_realloc (c_dev, tmp_i, stat=istat)
          IF (istat /= 0) THEN
             IF (verbose_acc) WRITE(*,*)routineN//" Running out of memory"
             CALL dbcsr_cuda_dev_mem_realloc (c_dev, maxs, stat=istat)
          ENDIF
       ENDIF cuda_realloc_c
    ENDIF
    CALL dbcsr_data_set_size_referenced (product_data_area, zero_last)
    !
    ! Zero new blocks
    IF (zero_last .GE. zero_first) THEN
       CALL dbcsr_data_clear (product_data_area, lb=zero_first, ub=zero_last)
       cuda_clear_c: IF (mm_driver .EQ. mm_driver_cuda) THEN
          CALL dbcsr_cuda_dev_mem_zero(c_dev,&
               first=zero_first, last=zero_last, error=error)
       ENDIF cuda_clear_c
    ENDIF
    zero_first = zero_last + 1
    !
    ! Resize block count.
    cuda_realloc_locks: IF (mm_driver .EQ. mm_driver_cuda) THEN
       IF (nblks .GT. dbcsr_cuda_dev_mem_get_alloc(c_locks_dev)) THEN
          CALL dbcsr_cuda_thread_sync(error=error)
          CALL dbcsr_cuda_dev_mem_dealloc(c_locks_dev, error=error)
          tmp_i = INT(REAL(nblks,kind=dp)*default_resize_factor)
          IF (verbose_acc) WRITE(*,*)routineN//" reallocating locks", tmp_i, &
               dbcsr_cuda_dev_mem_get_alloc(c_locks_dev)
          CALL dbcsr_cuda_dev_mem_alloc(c_locks_dev, tmp_i, stat=istat)
          IF (istat /= 0) THEN
             IF (verbose_acc) WRITE(*,*)routineN//" Trying smaller allocation"
             tmp_i = nblks
             CALL dbcsr_cuda_dev_mem_alloc(c_locks_dev, tmp_i, error=error)
          ENDIF
          CALL dbcsr_cuda_dev_mem_zero(c_locks_dev,&
               first=1, last=dbcsr_cuda_dev_mem_get_alloc(c_locks_dev),&
               error=error)
       ENDIF
    ENDIF cuda_realloc_locks
    !
    ! Verify stack consistency.  Only the upper bound is verified.
    IF (careful) THEN
       DO sp = 1, stack_size
          CALL dbcsr_assert (params(p_a_first,sp)&
               + params(p_m,sp) * params(p_k,sp) - 1,&
               "LE", dbcsr_data_get_size (left_data_area),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "A data out of bounds.", __LINE__, error=error)
          CALL dbcsr_assert (params(p_b_first,sp)&
               + params(p_k,sp) * params(p_n,sp) - 1,&
               "LE", dbcsr_data_get_size (right_data_area),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "B data out of bounds.", __LINE__, error=error)
          CALL dbcsr_assert (params(p_c_first,sp)&
               + params(p_m,sp) * params(p_n,sp) - 1,&
               "LE", dbcsr_data_get_size (product_data_area),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "C data out of bounds.", __LINE__, error=error)
       ENDDO
    ENDIF
    !
    SELECT CASE (mm_driver)
    CASE (mm_driver_matmul)
       SELECT CASE (product_data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL internal_process_mm_stack_s (params, &
               stack_size, &
               left_data_area%d%r_sp, right_data_area%d%r_sp, product_data_area%d%r_sp,&
               error=error)
       CASE (dbcsr_type_real_8)
          CALL internal_process_mm_stack_d (params,&
               stack_size,&
               left_data_area%d%r_dp, right_data_area%d%r_dp, product_data_area%d%r_dp,&
               error=error)
       CASE (dbcsr_type_complex_4)
          CALL internal_process_mm_stack_c (params,&
               stack_size,&
               left_data_area%d%c_sp, right_data_area%d%c_sp, product_data_area%d%c_sp,&
               error=error)
       CASE (dbcsr_type_complex_8)
          CALL internal_process_mm_stack_z (params,&
               stack_size,&
               left_data_area%d%c_dp, right_data_area%d%c_dp, product_data_area%d%c_dp,&
               error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    CASE (mm_driver_smm)
       SELECT CASE (product_data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL smm_process_mm_stack_s (params, &
               stack_size, &
               left_data_area%d%r_sp, right_data_area%d%r_sp, product_data_area%d%r_sp,&
               error=error)
       CASE (dbcsr_type_real_8)
          CALL smm_process_mm_stack_d (params,&
               stack_size,&
               left_data_area%d%r_dp, right_data_area%d%r_dp, product_data_area%d%r_dp,&
               error=error)
       CASE (dbcsr_type_complex_4)
          CALL smm_process_mm_stack_c (params,&
               stack_size,&
               left_data_area%d%c_sp, right_data_area%d%c_sp, product_data_area%d%c_sp,&
               error=error)
       CASE (dbcsr_type_complex_8)
          CALL smm_process_mm_stack_z (params,&
               stack_size,&
               left_data_area%d%c_dp, right_data_area%d%c_dp, product_data_area%d%c_dp,&
               error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    CASE (mm_driver_plasma)
       SELECT CASE (product_data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL plasma_process_mm_stack_s (params,&
               stack_size,&
               left_data_area%d%r_sp, right_data_area%d%r_sp, product_data_area%d%r_sp,&
               error=error)
       CASE (dbcsr_type_real_8)
          CALL plasma_process_mm_stack_d (params,&
               stack_size,&
               left_data_area%d%r_dp, right_data_area%d%r_dp, product_data_area%d%r_dp,&
               error=error)
       CASE (dbcsr_type_complex_4)
          CALL plasma_process_mm_stack_c (params,&
               stack_size,&
               left_data_area%d%c_sp, right_data_area%d%c_sp, product_data_area%d%c_sp,&
               error=error)
       CASE (dbcsr_type_complex_8)
          CALL plasma_process_mm_stack_z (params,&
               stack_size,&
               left_data_area%d%c_dp, right_data_area%d%c_dp, product_data_area%d%c_dp,&
               error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    CASE (mm_driver_blas)
       SELECT CASE (product_data_area%d%data_type)
       CASE (dbcsr_type_real_4)
          CALL blas_process_mm_stack_s (params,&
               stack_size,&
               left_data_area%d%r_sp, right_data_area%d%r_sp, product_data_area%d%r_sp,&
               error=error)
       CASE (dbcsr_type_real_8)
          CALL blas_process_mm_stack_d (params,&
               stack_size,&
               left_data_area%d%r_dp, right_data_area%d%r_dp, product_data_area%d%r_dp,&
               error=error)
       CASE (dbcsr_type_complex_4)
          CALL blas_process_mm_stack_c (params,&
               stack_size,&
               left_data_area%d%c_sp, right_data_area%d%c_sp, product_data_area%d%c_sp,&
               error=error)
       CASE (dbcsr_type_complex_8)
          CALL blas_process_mm_stack_z (params,&
               stack_size,&
               left_data_area%d%c_dp, right_data_area%d%c_dp, product_data_area%d%c_dp,&
               error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    CASE (mm_driver_cuda)
       !CALL swap (param_dev_calc, param_dev_xfer)
       !CALL swap (data_dev_calc, data_dev_xfer)
       CALL cuda_process_mm_stack (params,&
            stack_size,&
            a_dev(data_dev_calc), b_dev(data_dev_calc), c_dev,&
            c_locks_dev,&
            params_dev(param_dev_calc),&
            error=error)
       !SELECT CASE (product_data_area%d%data_type)
       !CASE (dbcsr_type_real_4)
       !   CALL smm_process_mm_stack_s (params, &
       !        stack_size, &
       !        left_data_area%d%r_sp, right_data_area%d%r_sp, product_data_area%d%r_sp,&
       !        error=error)
       !CASE (dbcsr_type_real_8)
       !   CALL cuda_process_mm_stack_d (params,&
       !        stack_size,&
       !        a_dev(data_dev_calc), b_dev(data_dev_calc), c_dev,&
       !        c_locks_dev,&
       !        params_dev(param_dev_calc),&
       !        error=error)
       !CASE (dbcsr_type_complex_4)
       !   CALL smm_process_mm_stack_c (params,&
       !        stack_size,&
       !        left_data_area%d%c_sp, right_data_area%d%c_sp, product_data_area%d%c_sp,&
       !        error=error)
       !CASE (dbcsr_type_complex_8)
       !   CALL smm_process_mm_stack_z (params,&
       !        stack_size,&
       !        left_data_area%d%c_dp, right_data_area%d%c_dp, product_data_area%d%c_dp,&
       !        error=error)
       !CASE default
       !   CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
       !        routineN, "Invalid data type",__LINE__,error)
       !END SELECT
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Invalid multiplication driver",__LINE__,error)
    END SELECT
    stack_size = 0
    t_process_stack = t_process_stack + m_walltime()
  END SUBROUTINE process_mm_stack

  SUBROUTINE print_gemm_parameters(params)
    INTEGER, DIMENSION(:, :), INTENT(in)     :: params

    INTEGER                                  :: sp

    DO sp = 1, SIZE(params)
       WRITE(*,*)"SP",sp
       WRITE(*,'(1X,A,1X,3(1X,I4,","),".",3(1X,I7,","))')&
            "GEMM PARAMETERS",&
               params(p_m,sp),&
               params(p_k,sp),&
               params(p_n,sp),&
               params(p_a_first,sp),&
               params(p_b_first,sp),&
               params(p_c_first,sp)
    ENDDO
  END SUBROUTINE print_gemm_parameters


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

#ifdef __PLASMA
  ELEMENTAL FUNCTION plasma_mat_type (t)
    LOGICAL, INTENT(IN)                      :: t
    INTEGER                                  :: plasma_mat_type

    INCLUDE 'plasmaf.h'

    IF (t) THEN
       plasma_mat_type = PlasmaTrans
    ELSE
       plasma_mat_type = PlasmaNoTrans
    ENDIF
  END FUNCTION plasma_mat_type
#endif

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
!> \brief Determines whether a transpose must be applied
!> \par Source
!> This function is copied from dbcsr_dist_operations for speed reasons.
!> \param[in] row   The absolute matrix row.
!> \param[in] column          The absolute matrix column.
! *****************************************************************************
  ELEMENTAL FUNCTION my_checker_tr(row, column) RESULT(transpose)
    INTEGER, INTENT(IN)                      :: row, column
    LOGICAL                                  :: transpose

    transpose = BTEST(column+row, 0) .EQV. column.GE.row

  END FUNCTION my_checker_tr

  SUBROUTINE xtime_set(timer)
    REAL(kind=dp), INTENT(out)               :: timer

    timer = -m_walltime()
  END SUBROUTINE xtime_set
  SUBROUTINE xtime_stop(timer)
    REAL(kind=dp), INTENT(inout)             :: timer

    timer = timer + m_walltime()
  END SUBROUTINE xtime_stop


! *****************************************************************************
!> \brief Fills row hashtable from an existing matrix.
!> \param[in] block_estimate guess for the number of blocks in the product matrix, can be zero
! *****************************************************************************
  SUBROUTINE fill_hash_tables(hashes, matrix, block_estimate, row_map, col_map, error)
    TYPE(hash_table_type), DIMENSION(:), &
      INTENT(inout)                          :: hashes
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER                                  :: block_estimate
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: row_map, col_map
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fill_hash_tables', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, error_handler, i, imat, &
                                                n_local_rows, row

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    imat = 1
    !$ imat = OMP_GET_THREAD_NUM() + 1
    n_local_rows = matrix%nblkrows_local
    CALL dbcsr_assert (SIZE(hashes), "EQ", n_local_rows,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Local row count mismatch", __LINE__, error=error)
    DO row = 1, n_local_rows
       ! create the hash table row with a reasonable initial size
       CALL hash_table_create (hashes(row), &
            MAX(8,(3*block_estimate)/MAX(1,n_local_rows)))
    ENDDO
    ! We avoid using the iterator because we will use the existing
    ! work matrix instead of the BCSR index.
    DO i = 1, matrix%wms(imat)%lastblk
       row = matrix%wms(imat)%row_i(i)
       col = matrix%wms(imat)%col_i(i)
       IF (PRESENT (row_map)) row = row_map(row)
       IF (PRESENT (col_map)) col = col_map(col)
       CALL hash_table_add(hashes(row), col, i, error=error)
    ENDDO
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE fill_hash_tables


! *****************************************************************************
!> \brief Adds blocks to a matrix
!>
!>        Existing blocks are replaced (overwritten).
!> \param[in,out] matrix_a   DBCSR matrix into which blocks are added
!> \param[in] matrix_b       DBCSR matrix from which blocks are added
!> \param[in,out] error      error
! *****************************************************************************
  SUBROUTINE dbcsr_insert_blocks(matrix_a, matrix_b, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_insert_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, data_type_b, &
                                                error_handler, nblkrows, &
                                                nblks, row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: b_row_i
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_block
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    ! Checks for validity
    CALL dbcsr_assert (dbcsr_valid_index (matrix_a),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Target matrix A must be valid.", __LINE__, error)
    CALL dbcsr_assert (dbcsr_valid_index (matrix_b),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Source matrix B must be valid.", __LINE__, error)
    ! Reserve the blocks to be added
    nblks = dbcsr_get_num_blocks (matrix_b)
    nblkrows = dbcsr_nblkrows_total (matrix_b)
    ALLOCATE (b_row_i(nblks))
    CALL dbcsr_expand_row_index (matrix_b%m%row_p, b_row_i, nblkrows, nblks)
    CALL dbcsr_reserve_blocks (matrix_a, b_row_i, matrix_b%m%col_i, error=error)
    DEALLOCATE (b_row_i)
    ! Prepare data structures
    data_type_b = dbcsr_get_data_type (matrix_b)
    ! Now add the blocks
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_new (data_block, data_type_b)
    CALL dbcsr_iterator_start(iter, matrix_b)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_block, tr, blk)
       CALL dbcsr_put_block(matrix_a, row, col, data_block, tr,&
            summation=.FALSE.)
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_data_clear_pointer (data_block)
    CALL dbcsr_data_release (data_block)
    !
    CALL dbcsr_finalize (matrix_a, error=error)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_insert_blocks


! *****************************************************************************
!> \brief Calculates per-block norms.
!>
!> Rewritten to be very low-level.
!> \param[in,out] matrix     DBCSR matrix for which to calculate norms
!> \param[in] norms          Block norms
!> \param[in,out] error      error
! *****************************************************************************
  SUBROUTINE calculate_norms(matrix, norms, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    REAL(kind=sp), DIMENSION(:), INTENT(OUT) :: norms
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_norms', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_type, error_handle, &
                                                nblks, nrows, row
    INTEGER, DIMENSION(1), TARGET            :: tmp
    INTEGER, DIMENSION(:), POINTER           :: local_cols, local_rows

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handle, error)
    ! Checks for validity
    CALL dbcsr_assert (dbcsr_valid_index (matrix),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "The matrix must be valid.", __LINE__, error)
    data_type = dbcsr_get_data_type (matrix)
    IF (matrix%m%local_indexing) THEN
       IF (careful_mod) &
            CALL dbcsr_assert (array_exists (matrix%m%local_rows),&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Global row mapping should exist", __LINE__, error=error)
       local_rows => array_data (matrix%m%local_rows)
       nrows = SIZE(local_rows)
       local_cols => array_data (matrix%m%local_cols)
    ELSE
       local_rows => tmp ! Have something valid to point to
       local_cols => tmp
       nrows = matrix%m%nblkrows_total
    ENDIF
    IF (matrix%m%list_indexing) THEN
       nblks = matrix%m%nblks
       SELECT CASE (data_type)
       CASE (dbcsr_type_real_4)
          CALL calc_norms_list_s(norms, nblks,&
               matrix%m%row_p, &
               array_data (matrix%m%row_blk_size),&
               array_data (matrix%m%col_blk_size),&
               dbcsr_get_data_p_s (matrix%m%data_area),&
               local=matrix%m%local_indexing,&
               local2global_rows=local_rows,&
               local2global_cols=local_cols)
       CASE (dbcsr_type_real_8)
          CALL calc_norms_list_d(norms, nblks,&
               matrix%m%row_p, &
               array_data (matrix%m%row_blk_size),&
               array_data (matrix%m%col_blk_size),&
               dbcsr_get_data_p_d (matrix%m%data_area),&
               local=matrix%m%local_indexing,&
               local2global_rows=local_rows,&
               local2global_cols=local_cols)
       CASE (dbcsr_type_complex_4)
          CALL calc_norms_list_c(norms, nblks,&
               matrix%m%row_p, &
               array_data (matrix%m%row_blk_size),&
               array_data (matrix%m%col_blk_size),&
               dbcsr_get_data_p_c (matrix%m%data_area),&
               local=matrix%m%local_indexing,&
               local2global_rows=local_rows,&
               local2global_cols=local_cols)
       CASE (dbcsr_type_complex_8)
          CALL calc_norms_list_z(norms, nblks,&
               matrix%m%row_p, &
               array_data (matrix%m%row_blk_size),&
               array_data (matrix%m%col_blk_size),&
               dbcsr_get_data_p_z (matrix%m%data_area),&
               local=matrix%m%local_indexing,&
               local2global_rows=local_rows,&
               local2global_cols=local_cols)
          CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
               routineN, "Invalid data type.",__LINE__,error)
       END SELECT
    ELSE
       SELECT CASE (data_type)
       CASE (dbcsr_type_real_4)
          CALL calc_norms_s(norms, nrows,&
               matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p,&
               array_data (matrix%m%row_blk_size),&
               array_data (matrix%m%col_blk_size),&
               dbcsr_get_data_p_s (matrix%m%data_area),&
               local=matrix%m%local_indexing,&
               local2global=local_rows)
       CASE (dbcsr_type_real_8)
          CALL calc_norms_d(norms, nrows,&
               matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p,&
               array_data (matrix%m%row_blk_size),&
               array_data (matrix%m%col_blk_size),&
               dbcsr_get_data_p_d (matrix%m%data_area),&
               local=matrix%m%local_indexing,&
               local2global=local_rows)
       CASE (dbcsr_type_complex_4)
          CALL calc_norms_c(norms, nrows,&
               matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p,&
               array_data (matrix%m%row_blk_size),&
               array_data (matrix%m%col_blk_size),&
               dbcsr_get_data_p_c (matrix%m%data_area),&
               local=matrix%m%local_indexing,&
               local2global=local_rows)
       CASE (dbcsr_type_complex_8)
          CALL calc_norms_z(norms, nrows,&
               matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p,&
               array_data (matrix%m%row_blk_size),&
               array_data (matrix%m%col_blk_size),&
               dbcsr_get_data_p_z (matrix%m%data_area),&
               local=matrix%m%local_indexing,&
               local2global=local_rows)
          CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
               routineN, "Invalid data type.",__LINE__,error)
       END SELECT
    ENDIF
       !
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE calculate_norms


! -----------------------------------------------------------------------------
! Beginning of hashtable
  ! finds a prime equal or larger than i
  FUNCTION matching_prime(i) RESULT(res)
    INTEGER, INTENT(IN)                      :: i
    INTEGER                                  :: res

    INTEGER                                  :: j

    res=i  
    j=0 
    DO WHILE (j<res) 
      DO j=2,res-1
         IF (MOD(res,j)==0) THEN
            res=res+1
            EXIT
         ENDIF
      ENDDO
    ENDDO
  END FUNCTION

  SUBROUTINE hash_table_create(hash_table,table_size) 
    TYPE(hash_table_type)                    :: hash_table
    INTEGER, INTENT(IN)                      :: table_size

    INTEGER                                  :: j

! guarantee a minimal hash table size (8), so that expansion works

   j=3
   DO WHILE(2**j-1<table_size)
      j=j+1
   ENDDO
   hash_table%nmax=2**j-1
   hash_table%prime=matching_prime(hash_table%nmax)
   hash_table%nele=0
   ALLOCATE(hash_table%table(0:hash_table%nmax))
  END SUBROUTINE hash_table_create

  SUBROUTINE hash_table_release(hash_table)
    TYPE(hash_table_type)                    :: hash_table

   hash_table%nmax=0
   hash_table%nele=0
   DEALLOCATE(hash_table%table)

  END SUBROUTINE hash_table_release

  RECURSIVE SUBROUTINE hash_table_add(hash_table,c,p, error)
    TYPE(hash_table_type), INTENT(INOUT)     :: hash_table
    INTEGER, INTENT(IN)                      :: c, p
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    REAL(KIND=real_8), PARAMETER :: hash_table_expand = 1.5_real_8, &
      inv_hash_table_fill = 2.5_real_8

    INTEGER                                  :: i, j
    TYPE(ele_type), ALLOCATABLE, &
      DIMENSION(:)                           :: tmp_hash

! if too small, make a copy and rehash in a larger table

    IF (hash_table%nele*inv_hash_table_fill>hash_table%nmax) THEN
       ALLOCATE(tmp_hash(LBOUND(hash_table%table,1):UBOUND(hash_table%table,1)))
       tmp_hash=hash_table%table
       CALL hash_table_release(hash_table) 
       CALL hash_table_create(hash_table,INT((UBOUND(tmp_hash,1)+8)*hash_table_expand))
       DO i=LBOUND(tmp_hash,1),UBOUND(tmp_hash,1)
          IF (tmp_hash(i)%c.NE.0) THEN
             CALL hash_table_add(hash_table,tmp_hash(i)%c,tmp_hash(i)%p,error)
          ENDIF
       ENDDO
       DEALLOCATE(tmp_hash)
    ENDIF

   hash_table%nele=hash_table%nele+1
   i=IAND(c*hash_table%prime,hash_table%nmax)

   DO j=i,hash_table%nmax
      IF (hash_table%table(j)%c==0 .OR. hash_table%table(j)%c==c) THEN
         hash_table%table(j)%c=c
         hash_table%table(j)%p=p
         RETURN
      ENDIF
   ENDDO
   DO j=0,i-1
      IF (hash_table%table(j)%c==0 .OR. hash_table%table(j)%c==c) THEN
         hash_table%table(j)%c=c
         hash_table%table(j)%p=p
         RETURN
      ENDIF
   ENDDO
  END SUBROUTINE hash_table_add

  PURE FUNCTION hash_table_get(hash_table,c) RESULT(p)
    TYPE(hash_table_type), INTENT(IN)        :: hash_table
    INTEGER, INTENT(IN)                      :: c
    INTEGER                                  :: p

    INTEGER                                  :: i, j

   i=IAND(c*hash_table%prime,hash_table%nmax)

   ! catch the likely case first
   IF (hash_table%table(i)%c==c) THEN
      p=hash_table%table(i)%p
      RETURN
   ENDIF

   DO j=i,hash_table%nmax
      IF (hash_table%table(j)%c==0 .OR. hash_table%table(j)%c==c) THEN
         p=hash_table%table(j)%p
         RETURN
      ENDIF
   ENDDO
   DO j=0,i-1
      IF (hash_table%table(j)%c==0 .OR. hash_table%table(j)%c==c) THEN
         p=hash_table%table(j)%p
         RETURN
      ENDIF
   ENDDO
  END FUNCTION hash_table_get

! End of hashtable
! -----------------------------------------------------------------------------

  pure SUBROUTINE memcpy (out, in, n)
    REAL(kind=real_8), DIMENSION(*), &
      INTENT(out)                            :: out
    REAL(kind=real_8), DIMENSION(*), &
      INTENT(in)                             :: in
    INTEGER, INTENT(in)                      :: n

    out(1:n) = in(1:n)
  END SUBROUTINE memcpy


  PURE SUBROUTINE local_filter (full_data, nle, local_elements, local_data)
    INTEGER, DIMENSION(:), INTENT(IN)        :: full_data
    INTEGER, INTENT(IN)                      :: nle
    INTEGER, DIMENSION(1:nle), INTENT(IN)    :: local_elements
    INTEGER, DIMENSION(1:nle), INTENT(OUT)   :: local_data

    INTEGER                                  :: l

    FORALL (l = 1 : nle)
       local_data(l) = full_data(local_elements(l))
    END FORALL
  END SUBROUTINE local_filter

  PURE SUBROUTINE local_filter_sp (full_data, nle, local_elements, local_data)
    REAL(KIND=sp), DIMENSION(:), INTENT(IN)  :: full_data
    INTEGER, INTENT(IN)                      :: nle
    INTEGER, DIMENSION(1:nle), INTENT(IN)    :: local_elements
    REAL(KIND=sp), DIMENSION(1:nle), &
      INTENT(OUT)                            :: local_data

    INTEGER                                  :: l

    FORALL (l = 1 : SIZE(local_data))
       local_data(l) = full_data(local_elements(l))
    END FORALL
  END SUBROUTINE local_filter_sp

#include "dbcsr_internal_operations_d.F"
#include "dbcsr_internal_operations_z.F"
#include "dbcsr_internal_operations_s.F"
#include "dbcsr_internal_operations_c.F"


! *****************************************************************************
!> \brief Processes MM stack using CUDA.
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
!> \param[in,out] error        error
! *****************************************************************************
  SUBROUTINE cuda_process_mm_stack(params,&
       stack_size,&
       data_a_dev, data_b_dev, data_c_dev,&
       c_locks,&
       params_dev,&
       error)
    INTEGER, INTENT(IN)                      :: stack_size
    INTEGER, &
      DIMENSION(n_mult_params*stack_size), &
      INTENT(IN), TARGET                     :: params
    TYPE(dbcsr_cuda_mem_type), INTENT(IN)    :: data_a_dev, data_b_dev
    TYPE(dbcsr_cuda_mem_type), INTENT(INOUT) :: data_c_dev
    TYPE(dbcsr_cuda_mem_type_i4), &
      INTENT(INOUT)                          :: c_locks, params_dev
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cuda_process_mm_stack', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, k_max, m_max, &
                                                n_max, sp
    INTEGER, DIMENSION(:), POINTER           :: params_p
    REAL(kind=dp)                            :: kt
    REAL(kind=dp), SAVE                      :: index_time = 0

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    IF (kernel_timing) THEN
       index_time = index_time + m_walltime()
       WRITE(*,'(1X,A20,1(1X,EN12.3,1X,A))')"Index time", index_time, "s"
    ENDIF
    params_p => params
    IF (kernel_timing) kt = -m_walltime()
    m_max = 0
    n_max = 0
    k_max = 0
    !last_c = 0
    DO sp = 1, stack_size
       m_max = MAX(m_max, params((sp-1)*n_mult_params+p_m))
       n_max = MAX(n_max, params((sp-1)*n_mult_params+p_n))
       k_max = MAX(k_max, params((sp-1)*n_mult_params+p_k))
    ENDDO
    IF (kernel_timing) THEN
       kt = kt + m_walltime()
       WRITE(*,'(1X,A20,1(1X,EN12.3,1X,A))')"Max size time", kt, "s"
    ENDIF
    IF (kernel_timing) kt = -m_walltime()
    CALL dbcsr_cuda_thread_sync (error=error)
    IF (kernel_timing) THEN
       kt = kt + m_walltime()
       WRITE(*,'(1X,A20,1(1X,EN12.3,1X,A))')"Thread sync time", kt, "s"
    ENDIF
    !
    IF (kernel_timing) kt = -m_walltime()
    CALL dbcsr_cuda_cp_host_to_dev (params_p,&
         params_dev, n_mult_params*stack_size, error=error)
    IF (kernel_timing) THEN
       kt = kt + m_walltime()
       WRITE(*,'(1X,A20,1(1X,EN12.3,1X,A))')"Parameter copy time", kt, "s"
    ENDIF
    !
    !WRITE(*,*)routineN//" calling kernel", stack_size,&
    !     dbcsr_cuda_dev_mem_get_type(data_c_dev)
    CALL m_flush(6)
    IF (kernel_timing) kt = -m_walltime()
    CALL dbcsr_cuda_do_mm_stack (params_dev, stack_size, n_mult_params,&
         data_a_dev, data_b_dev, data_c_dev,&
         c_locks, m_max, n_max, k_max,&
         error=error)
    IF (kernel_timing) THEN
       CALL dbcsr_cuda_thread_sync (error=error)
       kt = kt + m_walltime()
       WRITE(*,'(1X,A20,1(1X,EN12.3,1X,A))')"Kernel time", kt, "s"
       t_calc_step = t_calc_step + kt
    ENDIF
    IF (kernel_timing) index_time = -m_walltime()
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE cuda_process_mm_stack


END MODULE dbcsr_internal_operations

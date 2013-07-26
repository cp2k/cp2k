!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   First layer of the dbcsr matrix-matrix multiplication.
!>          It performs the MPI parallelization according to Cannon's algorithm.
!>	
!> \author  Urban Borstnik
!>
!> <b>Modification history:</b>
!>  - 2010-02-23 Moved from dbcsr_operations
!>  - 2011-11    Moved parameter-stack processing routines to
!>               dbcsr_mm_methods.
!>  - 2013-01    reorganized code (Ole Schuett)
! *****************************************************************************
MODULE dbcsr_mm_cannon
  USE array_types,                     ONLY: array_data,&
                                             array_exists,&
                                             array_hold
  USE dbcsr_config,                    ONLY: mm_driver,&
                                             mm_driver_blas,&
                                             mm_driver_cuda,&
                                             mm_driver_plasma,&
                                             use_CUDA_host_pinned_memory,&
                                             use_MPI_memory,&
                                             use_combined_types,&
                                             use_comm_thread
  USE dbcsr_cuda_types,                ONLY: dbcsr_cuda_mem_type
  USE dbcsr_data_methods,              ONLY: &
       dbcsr_data_clear_pointer, dbcsr_data_ensure_size, dbcsr_data_get_size, &
       dbcsr_data_get_size_referenced, dbcsr_data_init, dbcsr_data_new, &
       dbcsr_data_release, dbcsr_data_set_pointer, &
       dbcsr_data_set_size_referenced, dbcsr_get_data_p_c, &
       dbcsr_get_data_p_d, dbcsr_get_data_p_s, dbcsr_get_data_p_z
  USE dbcsr_dist_operations,           ONLY: dbcsr_reset_vlocals,&
                                             image_calculator
  USE dbcsr_error_handling
  USE dbcsr_index_operations,          ONLY: dbcsr_count_row_index,&
                                             dbcsr_has_local_row_index,&
                                             dbcsr_repoint_index
  USE dbcsr_io,                        ONLY: dbcsr_print
  USE dbcsr_kinds,                     ONLY: dp,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_4_size,&
                                             real_8,&
                                             real_8_size,&
                                             sp
  USE dbcsr_machine,                   ONLY: default_output_unit,&
                                             m_flush
  USE dbcsr_message_passing,           ONLY: &
       mp_allgather, mp_irecv, mp_isend, mp_request_null, mp_sum, mp_testany, &
       mp_type_descriptor_type, mp_type_free, mp_type_make, mp_waitall
  USE dbcsr_methods,                   ONLY: &
       dbcsr_destroy_array, dbcsr_distribution_mp, dbcsr_get_data_type, &
       dbcsr_get_index_memory_type, dbcsr_image_dist_init, dbcsr_init, &
       dbcsr_mp_grid_setup, dbcsr_mp_group, dbcsr_mp_has_subgroups, &
       dbcsr_mp_my_col_group, dbcsr_mp_my_row_group, dbcsr_mp_mynode, &
       dbcsr_mp_mypcol, dbcsr_mp_myprow, dbcsr_mp_npcols, dbcsr_mp_nprows, &
       dbcsr_mp_numnodes, dbcsr_mp_pgrid, dbcsr_nblkcols_local, &
       dbcsr_nblkcols_total, dbcsr_nblkrows_local, dbcsr_nblkrows_total, &
       dbcsr_nfullcols_local, dbcsr_nfullrows_local, dbcsr_valid_index
  USE dbcsr_mm_multrec,                ONLY: dbcsr_mm_multrec_finalize,&
                                             dbcsr_mm_multrec_init,&
                                             dbcsr_mm_multrec_lib_finalize,&
                                             dbcsr_mm_multrec_lib_init,&
                                             dbcsr_mm_multrec_multiply
  USE dbcsr_mm_types,                  ONLY: multrec_carrier_type
  USE dbcsr_mp_operations,             ONLY: dbcsr_irecv_any,&
                                             dbcsr_isend_any,&
                                             dbcsr_mp_type_from_anytype
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size
  USE dbcsr_types,                     ONLY: &
       dbcsr_2d_array_type, dbcsr_data_obj, dbcsr_memory_CUDA_host_pinned, &
       dbcsr_memory_MPI, dbcsr_memory_default, dbcsr_mp_obj, dbcsr_obj, &
       dbcsr_slot_size, dbcsr_type, dbcsr_type_complex_4, &
       dbcsr_type_complex_8, dbcsr_type_real_4, dbcsr_type_real_8
  USE dbcsr_util,                      ONLY: count_bins,&
                                             dbcsr_checksum
  USE dbcsr_work_operations,           ONLY: dbcsr_create

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_cannon'
  CHARACTER(len=*), PARAMETER, PRIVATE :: int_print = "(10(1X,I7))"
  LOGICAL, PARAMETER :: debug_mod  = .FALSE.
  LOGICAL, PARAMETER :: careful_mod = .FALSE.


  PUBLIC :: dbcsr_mm_cannon_lib_init, dbcsr_mm_cannon_lib_finalize
  PUBLIC :: dbcsr_mm_cannon_multiply


  CONTAINS

! *****************************************************************************
! \brief Initialize the library
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_cannon_lib_init(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

     CALL dbcsr_mm_multrec_lib_init(error)
  END SUBROUTINE


! *****************************************************************************
! \brief Finalize the library
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_cannon_lib_finalize(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

     CALL dbcsr_mm_multrec_lib_finalize(error)
  END SUBROUTINE


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
  SUBROUTINE dbcsr_mm_cannon_multiply(left_set, right_set, product_matrix,&
       error, retain_sparsity, &
       filter_eps, flop)
    TYPE(dbcsr_2d_array_type), POINTER       :: left_set, right_set
    TYPE(dbcsr_obj), INTENT(INOUT)           :: product_matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(kind=real_8), INTENT(in), OPTIONAL  :: filter_eps
    INTEGER(KIND=int_8), INTENT(OUT)         :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_cannon_multiply', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=80), PARAMETER :: &
      fdata = '(A,1X,I4,"(",2(I3),"x",2(I3),")","(",I3,"x",I3,")")', fxfer = &
      '(A,1X,I4,"->",I4,2(1X,"(",I3,"x",I3,")"),1X,"IM (",I3,"x",I3,")")'
    INTEGER, PARAMETER :: id_bytes = 3, id_recv = 2, id_send = 1, &
      id_time = 1, id_waittime = 2, idata = 1, ileft = 0, imeta = 2, &
      iright = 2, M_L = 2, M_P = 1, M_R = 3, RC_C = 2, RC_R = 1
    LOGICAL, PARAMETER                       :: excessive_output = .FALSE.

    INTEGER :: data_size, data_type, error_handler, error_handler2, &
      flop_metronome, grp, i, ithread, left_col_image, left_col_mult, &
      left_col_nimages, left_data_recv_size, left_data_send_size, &
      left_dst_icol, left_dst_irow, left_dst_p, left_dst_pcol, left_dst_prow, &
      left_dst_vcol, left_dst_vrow, left_index_recv_size, &
      left_index_send_size, left_max_nblks, left_max_nze, left_myfirstvcol, &
      left_myfirstvrow, left_mypcol, left_myprow, left_npcols, left_nprows, &
      left_recv_icol, left_recv_irow, left_recv_p, left_recv_pcol, &
      left_recv_prow, left_recv_vcol, left_recv_vrow, left_row_image, &
      left_row_mult, left_row_nimages, left_send_icol
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
      threads_finished, v_k, v_ki
    INTEGER(KIND=int_8)                      :: flop_single, flop_total, size8
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: row_counts, total_row_counts
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: left_sizes, my_sizes, &
                                                right_sizes
    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: all_sizes
    INTEGER, DIMENSION(3, 2)                 :: mp_rc_groups
    INTEGER, DIMENSION(:), POINTER           :: left_index_rp, left_index_sp, &
                                                local_rows, right_index_rp, &
                                                right_index_sp
    INTEGER, DIMENSION(:, :), POINTER :: left_data_rr, left_data_sr, &
      left_index_rr, left_index_sr, left_pgrid, product_pgrid, right_data_rr, &
      right_data_sr, right_index_rr, right_index_sr, right_pgrid
    INTEGER, SAVE                            :: mult_id = 0
    LOGICAL                                  :: keep_sparsity, list_indexing, &
                                                my_use_plasma, otf_filtering
    REAL(KIND=dp)                            :: checksum
    REAL(KIND=real_8)                        :: fill_guess, &
                                                fill_guess_global, left_fill, &
                                                right_fill
    REAL(kind=sp), ALLOCATABLE, DIMENSION(:) :: left_norms, right_norms, &
                                                row_max_epss
    TYPE(dbcsr_2d_array_type), POINTER :: left_buffer_1, left_buffer_2, &
      left_buffer_calc, left_buffer_comm, right_buffer_1, right_buffer_2, &
      right_buffer_calc, right_buffer_comm
    TYPE(dbcsr_cuda_mem_type), POINTER       :: a_dev, b_dev
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
    TYPE(multrec_carrier_type), SAVE         :: multrec_carrier

!$OMP THREADPRIVATE (multrec_carrier)


!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    !
    ALLOCATE (left_buffer_1, left_buffer_2, right_buffer_1, right_buffer_2)
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
                        SIZE(left_set%mats(left_row_image, left_col_image)%m%coo_l),&
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
                     left_set%mats(left_row_image, left_col_image)%m%coo_l(1::3),&
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
!$omp          shared (left_set, right_set,fill_guess_global, &
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
       !$OMP MASTER
       fill_guess_global = fill_guess
       !$OMP END MASTER
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
    left_data_sr = mp_request_null
    left_data_rr = mp_request_null
    left_index_sr = mp_request_null
    left_index_rr = mp_request_null
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
    right_data_sr = mp_request_null
    right_data_rr = mp_request_null
    right_index_sr = mp_request_null
    right_index_rr = mp_request_null
    !
!$omp parallel &
!$omp default (none) &
!$omp firstprivate (fill_guess_global) &
!$omp shared (left_buffer_comm, right_buffer_comm, product_matrix,&
!$omp         keep_sparsity, filter_eps, row_max_epss, &
!$omp         right_data_sr, right_data_rr, left_data_sr, left_data_rr,&
!$omp         right_index_sr, right_index_rr, left_index_sr, left_index_rr), &
!$omp shared (a_dev, b_dev) firstprivate (error)
    CALL dbcsr_mm_multrec_init(&
         left_buffer_comm%mats(1, 1)%m,&
         right_buffer_comm%mats(1, 1)%m,&
         product_matrix%m,&
         multrec_carrier,&
         keep_sparsity=keep_sparsity,&
         eps=filter_eps,&
         row_max_epss = row_max_epss,&
         fill_guess=fill_guess_global,error=error)
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
                CALL mp_irecv (right_recv_type(v_ki+1, 1), right_recv_p,&
                     grp, right_data_rr(v_ki+1, 1), tag=right_src_vrow)
                CALL mp_isend (right_send_type(v_ki+1, 1), right_send_p,&
                     grp, right_data_sr(v_ki+1, 1), tag=right_dst_vrow)
             ELSE
                CALL dbcsr_irecv_any (right_data_rp, right_recv_p,&
                     grp, right_data_rr(v_ki+1, 1), tag=right_src_vrow,&
                     error=error)
                CALL mp_irecv (right_index_rp, right_recv_p,&
                     grp, right_index_rr(v_ki+1, 1), tag=right_src_vrow)
                CALL dbcsr_isend_any (right_data_sp, right_send_p,&
                     grp, right_data_sr(v_ki+1, 1), tag=right_dst_vrow,&
                     error=error)
                CALL mp_isend (right_index_sp, right_send_p,&
                     grp, right_index_sr(v_ki+1, 1), tag=right_dst_vrow)
             ENDIF
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
                CALL mp_irecv (left_recv_type(1, v_ki+1), left_recv_p,&
                     grp, left_data_rr(1, v_ki+1), tag=left_src_vcol)
                CALL mp_isend (left_send_type(1, v_ki+1), left_send_p,&
                     grp, left_data_sr(1, v_ki+1), tag=left_dst_vcol)
             ELSE
                CALL dbcsr_irecv_any (left_data_rp, left_recv_p,&
                     grp, left_data_rr(1, v_ki+1), tag=left_src_vcol,&
                     error=error)
                CALL mp_irecv (left_index_rp, left_recv_p,&
                     grp, left_index_rr(1, v_ki+1), tag=left_src_vcol)
                CALL dbcsr_isend_any (left_data_sp, left_send_p,&
                     grp, left_data_sr(1, v_ki+1), tag=left_dst_vcol,&
                     error=error)
                CALL mp_isend (left_index_sp, left_send_p,&
                     grp, left_index_sr(1, v_ki+1), tag=left_dst_vcol)
             ENDIF
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
             ! form here the code for dbcsr_mm_driver_inner_init was taken 
             !
             IF (.FALSE.) WRITE(*,*)routineN//" TICK", v_ki
             IF (.TRUE. .OR. right_buffer_calc%mats(v_ki+1, 1)%m%local_indexing) THEN
                ! Since the right matrix is shifted vertically, the
                ! received data always has different notions of "local
                ! rows".  Thus the local_rows and global_rows must be
                ! recalculated.
                CALL dbcsr_reset_vlocals (right_buffer_calc%mats(v_ki+1, 1),&
                     right_set%image_dist, error=error)
             ENDIF
             IF (.TRUE. .OR. left_buffer_calc%mats(1, v_ki+1)%m%local_indexing) THEN
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
             flop_single = 0
             threads_finished = 0
!$omp parallel if( .NOT. my_use_plasma ) &
!$omp default (none) &
!$omp shared (left_buffer_calc, right_buffer_calc, &
!$omp         v_ki, &
!$omp         product_matrix, &
!$omp         filter_eps, right_norms, left_norms, row_max_epss, &
!$omp         keep_sparsity, error, threads_finished, &
!$omp         right_data_sr, right_data_rr, right_index_sr, right_index_rr, &
!$omp         left_data_sr, left_data_rr, left_index_sr, left_index_rr, use_comm_thread,error_handler2) &
!$omp private (ithread,nthreads, t_error) &
!$omp reduction (+: flop_single)
             t_error=error
             CALL dbcsr_error_set(routineN//"_multrec_sparse", error_handler2, t_error)
             CALL dbcsr_mm_multrec_multiply(&
                  left_buffer_calc%mats(1, v_ki+1)%m,&
                  right_buffer_calc%mats(v_ki+1, 1)%m,&
                  product_matrix%m,&
                  multrec_carrier,&
                  keep_sparsity=keep_sparsity,&             
                  flop=flop_single,&
                  a_norms=left_norms, b_norms=right_norms,&
                  error=t_error)
             ! once multiplication is complete, thread zero polls MPI until
             ! all other threads have also completed
!$omp atomic
             threads_finished = threads_finished + 1
             ithread = 0
             nthreads = 1
!$           ithread = omp_get_thread_num()
!$           nthreads = omp_get_num_threads()
             IF (use_comm_thread .AND. (ithread .EQ. 0)) THEN
               DO WHILE (threads_finished .NE. nthreads)
                  CALL mp_testany(right_data_sr)
                  CALL mp_testany(right_data_rr)
                  CALL mp_testany(left_data_sr)
                  CALL mp_testany(left_data_rr)
                  CALL mp_testany(right_index_sr)
                  CALL mp_testany(right_index_rr)
                  CALL mp_testany(left_index_sr)
                  CALL mp_testany(left_index_rr)
!$omp flush (threads_finished)
               END DO
             END IF
             !$OMP BARRIER
             CALL dbcsr_error_stop(error_handler2, t_error)
!$omp end parallel
             flop_total = flop_total + flop_single
             flop_metronome=flop_metronome+flop_single
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
    IF (excessive_output) CALL m_flush(1000000+mynode)
  
  
    !
!$omp parallel &                                                      
!$omp default (none) &
!$omp private (t_error) &
!$omp shared (product_matrix, error)
    !
    t_error = error
    ! Finalize multiplication
    CALL dbcsr_mm_multrec_finalize(multrec_carrier, product_matrix%m, t_error)
!$omp end parallel
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
    DEALLOCATE(left_data_rr, left_data_sr, left_index_rr, left_index_sr, &
               right_data_rr, right_data_sr, right_index_rr, right_index_sr)
    DEALLOCATE(left_send_type, left_recv_type, right_send_type, right_recv_type)
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
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_mm_cannon_multiply


! ******************************************************************************
!>
! ******************************************************************************
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
                  SIZE(matrix_set%mats(i_row, i_col)%m%coo_l),&
                  dbcsr_nblkrows_local(matrix_set%mats(i_row, i_col)),&
                  dbcsr_nblkcols_local(matrix_set%mats(i_row, i_col))
             WRITE(*,'(3(1X,I7))')matrix_set%mats(i_row, i_col)%m%coo_l
          ENDIF
          IF (careful_mod) THEN
             CALL dbcsr_assert (SIZE(matrix_set%mats(i_row, i_col)%m%coo_l),&
                  "EQ", matrix_set%mats(i_row, i_col)%m%nblks*3,&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Block count mismatch.", __LINE__, error=error)
          ENDIF
          thread_redist = ASSOCIATED (matrix_set%mats(i_row, i_col)%m%thr_c)
          t_size = SIZE(matrix_set%mats(i_row, i_col)%m%coo_l)/3
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
             IF (matrix_set%mats(i_row, i_col)%m%local_indexing) THEN
                CALL call_rec_sort_index (&
                     dbcsr_nblkrows_local(matrix_set%mats(i_row, i_col)),&
                     dbcsr_nblkcols_local(matrix_set%mats(i_row, i_col)),&
                     t_size,&
                     matrix_set%mats(i_row, i_col)%m%coo_l((t_f*3-2):(t_l*3)), error=error)
             ELSE
                CALL call_rec_sort_index (&
                     dbcsr_nblkrows_total(matrix_set%mats(i_row, i_col)),&
                     dbcsr_nblkcols_total(matrix_set%mats(i_row, i_col)),&
                     t_size,&
                     matrix_set%mats(i_row, i_col)%m%coo_l((t_f*3-2):(t_l*3)), error=error)
             ENDIF
          ENDIF
          !$OMP END PARALLEL
       ENDDO
    ENDDO
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE setup_rec_index



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
       WRITE(*,'(3(1X,I7))')a(:,1:nele)
    ENDIF
    IF (dbg) THEN
       IF (d .GT. 20) THEN
          WRITE(*,*)a(1,-d*1000)
       ENDIF
    ENDIF
    ALLOCATE(tmp(3,nele))
    M = mf-mi+1
    N = nf-ni+1
    IF (M > N) THEN
       half = M/2
       CALL rec_split (nele, a, tmp, 1, nlow, mi, mf, half)
       a = tmp
       DEALLOCATE (tmp)
       IF (nlow .GT. 1) THEN
          CALL rec_sort_index(mi,mi+half-1,ni,nf, nlow, a(:,1:nlow), d+1)
       ENDIF
       IF (nele-nlow .GT. 1) THEN
          CALL rec_sort_index(mi+half,mf,ni,nf, nele-nlow, a(:,nlow+1:nele), d+1)
       ENDIF
    ELSE
       half = N/2
       CALL rec_split (nele, a, tmp, 2, nlow, ni, nf, half)
       a = tmp
       DEALLOCATE (tmp)
       IF (nlow .GT. 1) THEN
         CALL rec_sort_index(mi,mf,ni,ni+half-1, nlow, a(:,1:nlow), d+1)
       ENDIF
       IF (nele-nlow .GT. 1) THEN
         CALL rec_sort_index(mi,mf,ni+half,nf, nele-nlow, a(:,nlow+1:nele), d+1)
       ENDIF
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
               matrix%m%coo_l, &
               array_data (matrix%m%row_blk_size),&
               array_data (matrix%m%col_blk_size),&
               dbcsr_get_data_p_s (matrix%m%data_area),&
               local=matrix%m%local_indexing,&
               local2global_rows=local_rows,&
               local2global_cols=local_cols)
       CASE (dbcsr_type_real_8)
          CALL calc_norms_list_d(norms, nblks,&
               matrix%m%coo_l, &
               array_data (matrix%m%row_blk_size),&
               array_data (matrix%m%col_blk_size),&
               dbcsr_get_data_p_d (matrix%m%data_area),&
               local=matrix%m%local_indexing,&
               local2global_rows=local_rows,&
               local2global_cols=local_cols)
       CASE (dbcsr_type_complex_4)
          CALL calc_norms_list_c(norms, nblks,&
               matrix%m%coo_l, &
               array_data (matrix%m%row_blk_size),&
               array_data (matrix%m%col_blk_size),&
               dbcsr_get_data_p_c (matrix%m%data_area),&
               local=matrix%m%local_indexing,&
               local2global_rows=local_rows,&
               local2global_cols=local_cols)
       CASE (dbcsr_type_complex_8)
          CALL calc_norms_list_z(norms, nblks,&
               matrix%m%coo_l, &
               array_data (matrix%m%row_blk_size),&
               array_data (matrix%m%col_blk_size),&
               dbcsr_get_data_p_z (matrix%m%data_area),&
               local=matrix%m%local_indexing,&
               local2global_rows=local_rows,&
               local2global_cols=local_cols)
       CASE DEFAULT
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
       CASE DEFAULT
          CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
               routineN, "Invalid data type.",__LINE__,error)
       END SELECT
    ENDIF
       !
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE calculate_norms

#include "dbcsr_mm_cannon_d.F"
#include "dbcsr_mm_cannon_z.F"
#include "dbcsr_mm_cannon_s.F"
#include "dbcsr_mm_cannon_c.F"

END MODULE dbcsr_mm_cannon

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  USE acc_event,                       ONLY: acc_event_record,&
                                             acc_event_synchronize,&
                                             acc_stream_wait_event
  USE acc_stream,                      ONLY: acc_stream_associated,&
                                             acc_stream_create,&
                                             acc_stream_destroy,&
                                             acc_stream_type
  USE array_types,                     ONLY: array_data,&
                                             array_equality,&
                                             array_exists,&
                                             array_hold,&
                                             array_i1d_obj,&
                                             array_nullify,&
                                             array_release
  USE dbcsr_acc_operations,            ONLY: dbcsr_acc_transpose
  USE dbcsr_block_operations,          ONLY: dbcsr_block_conjg,&
                                             dbcsr_block_copy_aa,&
                                             dbcsr_block_real_neg,&
                                             dbcsr_block_scale,&
                                             dbcsr_block_transpose_aa
  USE dbcsr_config,                    ONLY: dbcsr_get_conf_avg_elements_images,&
                                             default_resize_factor,&
                                             mm_async,&
                                             mm_driver,&
                                             mm_driver_acc,&
                                             use_comm_thread,&
                                             use_mpi_rma
  USE dbcsr_data_methods,              ONLY: &
       dbcsr_data_clear_pointer, dbcsr_data_copy, dbcsr_data_ensure_size, &
       dbcsr_data_get_size, dbcsr_data_get_size_referenced, dbcsr_data_hold, &
       dbcsr_data_host2dev, dbcsr_data_init, dbcsr_data_new, &
       dbcsr_data_release, dbcsr_data_set_pointer, &
       dbcsr_data_set_size_referenced, dbcsr_get_data_p_c, &
       dbcsr_get_data_p_d, dbcsr_get_data_p_s, dbcsr_get_data_p_z, &
       dbcsr_scalar, dbcsr_scalar_are_equal, dbcsr_scalar_fill_all, &
       dbcsr_scalar_negative, dbcsr_scalar_one, dbcsr_scalar_set_type, &
       dbcsr_scalar_zero
  USE dbcsr_data_types,                ONLY: dbcsr_datatype_sizeof
  USE dbcsr_dist_methods,              ONLY: &
       dbcsr_distribution_col_dist, dbcsr_distribution_data_max_size_allproc, &
       dbcsr_distribution_has_threads, dbcsr_distribution_hold, &
       dbcsr_distribution_make_threads, dbcsr_distribution_mp, &
       dbcsr_distribution_ncols, dbcsr_distribution_no_threads, &
       dbcsr_distribution_nrows, dbcsr_distribution_release, &
       dbcsr_distribution_row_dist
  USE dbcsr_dist_operations,           ONLY: dbcsr_create_image_dist,&
                                             dbcsr_get_local_vcols,&
                                             dbcsr_get_local_vrows,&
                                             dbcsr_make_dists_dense,&
                                             dbcsr_reset_locals,&
                                             dbcsr_reset_vlocals,&
                                             image_calculator,&
                                             make_sizes_dense
  USE dbcsr_error_handling,            ONLY: &
       dbcsr_assert, dbcsr_caller_error, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_failure_level, dbcsr_fatal_level, &
       dbcsr_internal_error, dbcsr_unimplemented_error_nr, &
       dbcsr_warning_level, dbcsr_wrong_args_error
  USE dbcsr_index_operations,          ONLY: dbcsr_count_row_index,&
                                             dbcsr_has_local_row_index,&
                                             dbcsr_index_compact,&
                                             dbcsr_index_prune_deleted,&
                                             dbcsr_make_index_canonical,&
                                             dbcsr_make_index_list,&
                                             dbcsr_make_index_local_row,&
                                             dbcsr_repoint_index
  USE dbcsr_io,                        ONLY: dbcsr_print
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_mem_methods,               ONLY: dbcsr_mempool_clear,&
                                             dbcsr_mempool_destruct,&
                                             dbcsr_mempool_ensure_capacity,&
                                             dbcsr_memtype_setup
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_offsets, dbcsr_col_block_sizes, dbcsr_destroy_array, &
       dbcsr_distribution, dbcsr_get_data_type, dbcsr_get_index_memory_type, &
       dbcsr_get_matrix_type, dbcsr_has_symmetry, dbcsr_image_dist_hold, &
       dbcsr_image_dist_init, dbcsr_image_dist_release, dbcsr_init, &
       dbcsr_nblkcols_local, dbcsr_nblkcols_total, dbcsr_nblkrows_local, &
       dbcsr_nblkrows_total, dbcsr_nfullcols_total, dbcsr_nfullrows_total, &
       dbcsr_release, dbcsr_release_locals, dbcsr_row_block_offsets, &
       dbcsr_valid_index
  USE dbcsr_mm_multrec,                ONLY: dbcsr_mm_multrec_finalize,&
                                             dbcsr_mm_multrec_init,&
                                             dbcsr_mm_multrec_lib_finalize,&
                                             dbcsr_mm_multrec_lib_init,&
                                             dbcsr_mm_multrec_multiply,&
                                             dbcsr_mm_multrec_phaseout,&
                                             dbcsr_mm_multrec_type
  USE dbcsr_mp_methods,                ONLY: &
       dbcsr_mp_grid_setup, dbcsr_mp_group, dbcsr_mp_has_subgroups, &
       dbcsr_mp_my_col_group, dbcsr_mp_my_row_group, dbcsr_mp_mynode, &
       dbcsr_mp_mypcol, dbcsr_mp_myprow, dbcsr_mp_npcols, dbcsr_mp_nprows, &
       dbcsr_mp_numnodes, dbcsr_mp_pgrid
  USE dbcsr_mp_operations,             ONLY: dbcsr_irecv_any,&
                                             dbcsr_isend_any,&
                                             dbcsr_rget_any,&
                                             dbcsr_win_create_any,&
                                             hybrid_alltoall_any,&
                                             hybrid_alltoall_i1
  USE dbcsr_operations,                ONLY: dbcsr_conjg,&
                                             dbcsr_copy,&
                                             dbcsr_crop_matrix,&
                                             dbcsr_filter,&
                                             dbcsr_may_be_dense,&
                                             dbcsr_scale
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size
  USE dbcsr_toollib,                   ONLY: uppercase
  USE dbcsr_transformations,           ONLY: dbcsr_make_dense,&
                                             dbcsr_make_dense_low,&
                                             dbcsr_make_undense,&
                                             dbcsr_make_untransposed_blocks,&
                                             dbcsr_new_transposed
  USE dbcsr_types,                     ONLY: &
       dbcsr_2d_array_type, dbcsr_conjugate_transpose, dbcsr_data_obj, &
       dbcsr_distribution_obj, dbcsr_imagedistribution_obj, dbcsr_iterator, &
       dbcsr_memtype_default, dbcsr_memtype_type, dbcsr_meta_size, &
       dbcsr_mp_obj, dbcsr_mpi_size_limits, dbcsr_mpi_statistics_type, &
       dbcsr_no_transpose, dbcsr_obj, dbcsr_scalar_type, &
       dbcsr_slot_home_coli, dbcsr_slot_home_pcol, dbcsr_slot_home_prow, &
       dbcsr_slot_home_rowi, dbcsr_slot_home_vpcol, dbcsr_slot_home_vprow, &
       dbcsr_slot_size, dbcsr_transpose, dbcsr_type, &
       dbcsr_type_antisymmetric, dbcsr_type_complex_4, dbcsr_type_complex_8, &
       dbcsr_type_int_4, dbcsr_type_no_symmetry, dbcsr_type_real_4, &
       dbcsr_type_real_8
  USE dbcsr_util,                      ONLY: count_bins,&
                                             dbcsr_checksum,&
                                             dbcsr_verify_matrix
  USE dbcsr_work_operations,           ONLY: dbcsr_add_wm_from_matrix,&
                                             dbcsr_create,&
                                             dbcsr_finalize,&
                                             dbcsr_special_finalize,&
                                             dbcsr_work_create
  USE kinds,                           ONLY: dp,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8,&
                                             sp
  USE machine,                         ONLY: default_output_unit
  USE message_passing,                 ONLY: &
       mp_address_kind, mp_allgather, mp_alltoall, mp_irecv, mp_isend, &
       mp_isendrecv, mp_max, mp_min, mp_request_null, mp_rget, mp_sum, &
       mp_testany, mp_waitall, mp_win_free, mp_win_lock_all, mp_win_unlock_all

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_cannon'
  CHARACTER(len=*), PARAMETER, PRIVATE :: int_print = "(10(1X,I7))"
  LOGICAL, PARAMETER :: debug_mod  = .FALSE.
  LOGICAL, PARAMETER :: careful_mod = .FALSE.

  TYPE dbcsr_memtype_type_p
     TYPE(dbcsr_memtype_type), POINTER :: p => Null()
     ! ensure that array-elements are on different cache lines
     INTEGER(kind=int_4), DIMENSION(64)    :: padding
  END TYPE dbcsr_memtype_type_p

  TYPE(dbcsr_memtype_type_p), DIMENSION(:), POINTER, SAVE :: memtype_product_wm => Null()

  TYPE(dbcsr_mpi_statistics_type), PRIVATE, SAVE :: dbcsr_mpi_statistics
  INTEGER(KIND=int_8),             PRIVATE, SAVE :: marketing_flops = 0

  INTEGER, PRIVATE, SAVE :: last_mpi_ranks_used = 0

  TYPE(dbcsr_memtype_type),     PRIVATE, SAVE  :: memtype_abpanel_1, memtype_abpanel_2,&
                                                  memtype_trsbuffer_1, memtype_trsbuffer_2, &
                                                  memtype_mpi_buffer
  TYPE(acc_stream_type), PRIVATE, SAVE         :: stream_1, stream_2
  ! ab-panels and streams are shared between all threads

  TYPE dbcsr_mm_multrec_type_p
    TYPE(dbcsr_mm_multrec_type), POINTER :: p => Null()
    ! ensure that array-elements are on different cache lines
    INTEGER(kind=int_4), DIMENSION(64)       :: padding
  END TYPE dbcsr_mm_multrec_type_p

  PUBLIC :: dbcsr_mm_cannon_lib_init, dbcsr_mm_cannon_lib_finalize
  PUBLIC :: dbcsr_mm_cannon_clear_mempools
  PUBLIC :: dbcsr_mm_cannon_multiply

  CONTAINS

! *****************************************************************************
!> \brief Initialize the library
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_cannon_lib_init(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: ithread, nthreads

    nthreads = 1; ithread = 0
    !$ nthreads = OMP_GET_NUM_THREADS () ; ithread = OMP_GET_THREAD_NUM ()

    CALL dbcsr_mm_multrec_lib_init(error)

    !$OMP MASTER
    dbcsr_mpi_statistics%col_nimages = -1
    dbcsr_mpi_statistics%row_nimages = -1
    dbcsr_mpi_statistics%nexchanged = 0
    dbcsr_mpi_statistics%nfiltered = 0
    dbcsr_mpi_statistics%nprocs = -1
    dbcsr_mpi_statistics%data_size = 0
    dbcsr_mpi_statistics%data_size(:,2) = HUGE(dbcsr_mpi_statistics%data_size(1,2))
    dbcsr_mpi_statistics%index_size = 0
    dbcsr_mpi_statistics%index_size(:,2) = HUGE(dbcsr_mpi_statistics%index_size(1,2))
    dbcsr_mpi_statistics%data_size_breakdown = 0

    marketing_flops = 0
    ALLOCATE(memtype_product_wm(0:nthreads-1))
    !$OMP END MASTER
    !$OMP BARRIER

    ! Each thread has its own working-matrix and its own mempool
    ALLOCATE(memtype_product_wm(ithread)%p)
    CALL dbcsr_memtype_setup(memtype_product_wm(ithread)%p, has_pool=.TRUE., error=error)
    CALL dbcsr_mempool_ensure_capacity(memtype_product_wm(ithread)%p%pool, capacity=1)
  END SUBROUTINE dbcsr_mm_cannon_lib_init


! *****************************************************************************
!> \brief Finalize the library
!> \param group ...
!> \param output_unit ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_cannon_lib_finalize(group, output_unit, error)
    INTEGER, INTENT(IN)                      :: group, output_unit
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_cannon_lib_finalize', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=50)                        :: image_grid_string
    INTEGER                                  :: ilimit, isqrt, isqrt2, ithread
    INTEGER(KIND=int_8)                      :: average, &
                                                total_marketing_flops, &
                                                total_nexchanged, &
                                                total_nfiltered
    INTEGER(KIND=int_8), DIMENSION(2)        :: max_recv_data, min_recv_data, &
                                                total_recv_data
    INTEGER(KIND=int_8), DIMENSION(2, SIZE(&
      dbcsr_mpi_size_limits)+1)              :: total_recv_breakdown

     CALL dbcsr_mm_multrec_lib_finalize(group, output_unit, error)

     ithread = 0
     !$  ithread = omp_get_thread_num()

     ! Each thread has its own working-matrix and its own mempool
     IF (ASSOCIATED(memtype_product_wm(ithread)%p%pool)) &
        CALL dbcsr_mempool_destruct(memtype_product_wm(ithread)%p%pool, error)
     DEALLOCATE(memtype_product_wm(ithread)%p)

     !$OMP BARRIER
     !$omp master
     DEALLOCATE(memtype_product_wm)

     ! this could overflow
     total_marketing_flops=marketing_flops
     CALL mp_sum(total_marketing_flops,group)
     
     total_nexchanged=dbcsr_mpi_statistics%nexchanged
     CALL mp_sum(total_nexchanged,group)

     total_nfiltered=dbcsr_mpi_statistics%nfiltered
     CALL mp_sum(total_nfiltered,group)

     total_recv_data(:)=dbcsr_mpi_statistics%data_size(:,1)
     CALL mp_sum(total_recv_data,group)

     min_recv_data(:)=dbcsr_mpi_statistics%data_size(:,2)
     CALL mp_min(min_recv_data,group)

     max_recv_data(:)=dbcsr_mpi_statistics%data_size(:,3)
     CALL mp_max(max_recv_data,group)

     IF (dbcsr_mpi_statistics%nexchanged.GT.0) THEN
        average = NINT(DBLE(SUM(total_recv_data(:)))/&
                       (dbcsr_mpi_statistics%nprocs*dbcsr_mpi_statistics%nexchanged))
     ELSE
        average=0
        min_recv_data=0
     ENDIF

     total_recv_breakdown(:,:)=dbcsr_mpi_statistics%data_size_breakdown(:,:)
     CALL mp_sum(total_recv_breakdown,group)

     IF (output_unit>0) THEN
       WRITE (output_unit,'(A,T30,I20)') " marketing flops", total_marketing_flops

       WRITE (UNIT=output_unit,FMT="(T2,A)") REPEAT("-",79)       
       WRITE (image_grid_string,'(I4,A,I4)') &
            dbcsr_mpi_statistics%row_nimages, " x ", dbcsr_mpi_statistics%col_nimages 
       WRITE (output_unit,'(A,T39,A)') " image size (#rows x #cols)",image_grid_string
       WRITE (output_unit,'(A,T30,I20)') " # MPI data exchanged",total_nexchanged
       WRITE (output_unit,'(A,T30,I20)') " # MPI data filtered",total_nfiltered
       WRITE (output_unit,'(A)') " MPI data size (elements):"
       WRITE (output_unit,'(A,T30,I20)') "  total",&
            SUM(total_recv_data(:))
       WRITE (output_unit,'(A,T30,I20)') "  min",&
            MINVAL(min_recv_data(:))
       WRITE (output_unit,'(A,T30,I20)') "  max",&
            MAXVAL(max_recv_data(:))
       WRITE (output_unit,'(A,T30,I20)') "  average",average

       WRITE (output_unit,'(A)') " MPI breakdown data size (bytes):"
       WRITE (output_unit,'(A,I8,T40,I10)') "             size <= ",dbcsr_mpi_size_limits(1),&
            SUM(total_recv_breakdown(:,1))
       DO ilimit=2,SIZE(dbcsr_mpi_size_limits)
          WRITE (output_unit,'(A,I8,A,I8,T40,I10)') "  ",dbcsr_mpi_size_limits(ilimit-1),&
               " < size <= ",dbcsr_mpi_size_limits(ilimit),&
               SUM(total_recv_breakdown(:,ilimit))
       ENDDO
       ilimit=SIZE(dbcsr_mpi_size_limits)
       WRITE (output_unit,'(A,I8,A,T40,I10)') "  ",dbcsr_mpi_size_limits(ilimit),&
            " < size    ",SUM(total_recv_breakdown(:,ilimit+1))

       isqrt=NINT(SQRT(REAL(last_mpi_ranks_used,KIND=real_8)))
       isqrt2=NINT(SQRT(REAL(last_mpi_ranks_used*2,KIND=real_8)))
       IF (isqrt*isqrt .NE. last_mpi_ranks_used) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A)") REPEAT("-",79)
          WRITE (UNIT=output_unit,FMT="(T2,A)") &
            "Warning: using a non-square number of MPI ranks might lead to poor performance."
          WRITE (UNIT=output_unit,FMT="(T2,A,I0)") &
            "         used ranks: ",last_mpi_ranks_used
          WRITE (UNIT=output_unit,FMT="(T2,A,2(I0,1X))") &
            "         suggested : ",isqrt**2,isqrt2**2
       ENDIF
     ENDIF
     IF (ASSOCIATED(memtype_trsbuffer_1%pool)) &
        CALL dbcsr_mempool_destruct(memtype_trsbuffer_1%pool, error)
     IF (ASSOCIATED(memtype_trsbuffer_2%pool)) &
        CALL dbcsr_mempool_destruct(memtype_trsbuffer_2%pool, error)
     IF (ASSOCIATED(memtype_abpanel_1%pool)) &
        CALL dbcsr_mempool_destruct(memtype_abpanel_1%pool, error)
     IF (ASSOCIATED(memtype_abpanel_2%pool)) &
        CALL dbcsr_mempool_destruct(memtype_abpanel_2%pool, error)
     IF(acc_stream_associated(stream_1)) &
        CALL acc_stream_destroy(stream_1)
     IF(acc_stream_associated(stream_2)) &
        CALL acc_stream_destroy(stream_2)
     !$omp end master
  END SUBROUTINE dbcsr_mm_cannon_lib_finalize

! *****************************************************************************
!> \brief Deallocate memory contained in mempools
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_cannon_clear_mempools(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: ithread

     ithread = 0
     !$  ithread = omp_get_thread_num()

     ! Each thread has its own working-matrix and its own mempool
     IF (ASSOCIATED(memtype_product_wm(ithread)%p%pool)) &
        CALL dbcsr_mempool_clear(memtype_product_wm(ithread)%p%pool, error)

     !$omp master
     IF (ASSOCIATED(memtype_trsbuffer_1%pool)) &
        CALL dbcsr_mempool_clear(memtype_trsbuffer_1%pool, error)
     IF (ASSOCIATED(memtype_trsbuffer_2%pool)) &
        CALL dbcsr_mempool_clear(memtype_trsbuffer_2%pool, error)
     IF (ASSOCIATED(memtype_abpanel_1%pool)) &
        CALL dbcsr_mempool_clear(memtype_abpanel_1%pool, error)
     IF (ASSOCIATED(memtype_abpanel_2%pool)) &
        CALL dbcsr_mempool_clear(memtype_abpanel_2%pool, error)
     !$omp end master
  END SUBROUTINE dbcsr_mm_cannon_clear_mempools

! *****************************************************************************
!> \brief Performs a multiplication of two dbcsr_type matrices,
!>        as  C := alpha * op( A ) * op( B ) + beta * C.
!> \param[in] transa specifies the form of op( A ) to be used in
!>                            the matrix multiplication
!>                            transa = 'N' or 'n',  op( A ) = A.
!>                            transa = 'T' or 't',  op( A ) = transpose(A).
!>                            transa = 'C' or 'c',  op( A ) = transpose(conjg(A)).
!> \param[in] transb specifies the form of op( B ) to be used in
!>                            the matrix multiplication
!>                            transb = 'N' or 'n',  op( B ) = B.
!>                            transb = 'T' or 't',  op( B ) = transpose(B).
!>                            transb = 'C' or 'c',  op( B ) = transpose(conjg(B)).
!> \param[in] alpha           scaling of product
!> \param[in] matrix_a        left BCSR matrix
!> \param[in] matrix_b        right BCSR matrix
!> \param[in] beta            scaling of existing data
!> \param[out] matrix_c       resulting BCSR product matrix.
!> \param[in] first_row       (optional) first full row of limiting submatrix
!> \param[in] last_row        (optional) last full row of limiting submatrix
!> \param[in] first_column    (optional) first full column of limiting submatrix
!> \param[in] last_column     (optional) last full column of limiting submatrix
!> \param[in] first_k         (optional) first full column of limiting inner
!>                            product
!> \param[in] last_k          (optional) last full column of limiting inner
!>                            product
!> \param[in] retain_sparsity (optional) enforce the sparsity pattern of the
!>                            existing product matrix; default is no
!> \param[in] filter_eps      Filtering of the matrix
!> \param[in,out] error       error
!> \param[out] flop           (optional) effective flop
!> \par Matrices m_a and m_b are multiplied into the m_c product matrix. If the
!>      dist2d parameter is not specified, then a new distribution_2d is
!>      determined for it.
!> \par Non-equal column dimensions of the right and product matrices
!>      The right and product matrix are allowed to have different
!>      (full) column dimensions. If they differ, there are certain
!>      peculiar behaviors, then the last_column is effectively set to
!>      the minimal of the two.
!> \par Beta scaling of the right product matrix
!>      If the effective last_column is less than the full column
!>      dimension of the product matrix, then the scaling of the
!>      product matrix with beta is limited to the submatrix specified
!>      by last_column.
!> \par Filtering
!>      The filter_eps parameter, if present, is used to filter the
!>      resulting matrix.  The filtering criterion is whether the
!>      block-frobenius norm is less than the specified epsilon.
!>      One-the-fly filtering is done such that individual
!>      multiplications are skipped if the product of the frobenius
!>      norms of the left- and right-matrix blocks are less than the
!>      specified epsilon divided by the maximum number of possible
!>      multiplies in each row.  In addition a final filtering is done
!>      as well with the same epsilon value.
! *****************************************************************************
  SUBROUTINE dbcsr_mm_cannon_multiply(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, filter_eps,&
       error, flop)

    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(KIND=real_8), INTENT(IN), OPTIONAL  :: filter_eps
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_cannon_multiply', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE., &
                                                use_list_indexing = .TRUE., &
                                                use_local_indexing = .TRUE.
    REAL(real_8), PARAMETER                  :: make_dense_occ_thresh = 1.0_dp

    CHARACTER                                :: transa_l, transb_l
    INTEGER :: comm, error_handler, error_handler2, f_col, f_k, f_row, &
      ithread, l_col, l_k, l_row, mult_images, numnodes, output_unit
    INTEGER(KIND=int_8)                      :: my_flop
    LOGICAL :: ab_dense, keep_product_data, keep_sparsity, new_left, &
      new_right, product_reindex, release_tdist, use_dense_mult
    REAL(KIND=dp)                            :: cs
    TYPE(array_i1d_obj) :: dense_col_sizes, dense_k_sizes, dense_row_sizes, &
      k_vmap, m_map, n_map, old_product_col_blk_offsets, &
      old_product_col_blk_sizes, old_product_row_blk_offsets, &
      old_product_row_blk_sizes
    TYPE(dbcsr_2d_array_type), POINTER       :: m2s_left, m2s_right
    TYPE(dbcsr_distribution_obj)             :: dense_product_distribution, &
                                                old_product_distribution
    TYPE(dbcsr_imagedistribution_obj)        :: dense_rdist_left, &
                                                dense_rdist_right, &
                                                rdist_left, rdist_right
    TYPE(dbcsr_obj)                          :: matrix_left, matrix_right, &
                                                product_matrix
    TYPE(dbcsr_scalar_type)                  :: eps_any

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL array_nullify(dense_k_sizes)
    CALL array_nullify(dense_col_sizes)
    CALL array_nullify(dense_row_sizes)

    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM ()

    ! setup driver-dependent memory-types and their memory-pools ---------------

    ! the ab_buffers are shared by all threads
    IF (mm_driver==mm_driver_acc) THEN
       IF(.NOT. acc_stream_associated(stream_1)) THEN
          CALL acc_stream_create(stream_1, "MemCpy (odd ticks)")
          CALL acc_stream_create(stream_2, "MemCpy (even ticks)")
       ENDIF

       CALL dbcsr_memtype_setup(memtype_abpanel_1, has_pool=.TRUE.,&
            acc_hostalloc=.TRUE., acc_devalloc=.TRUE., acc_stream=stream_1,&
            mpi=.TRUE., oversize_factor=default_resize_factor, error=error)

       CALL dbcsr_memtype_setup(memtype_abpanel_2, has_pool=.TRUE.,&
            acc_hostalloc=.TRUE., acc_devalloc=.TRUE., acc_stream=stream_2,&
            mpi=.TRUE., oversize_factor=default_resize_factor, error=error)

       !TODO: ensure capacity 2/3?
       CALL dbcsr_memtype_setup(memtype_trsbuffer_1,has_pool=.TRUE.,&
            acc_hostalloc=.TRUE., acc_devalloc=.TRUE., acc_stream=stream_1,error=error)
       CALL dbcsr_mempool_ensure_capacity(memtype_trsbuffer_1%pool, capacity=1)

       CALL dbcsr_memtype_setup(memtype_trsbuffer_2, has_pool=.TRUE.,&
            acc_hostalloc=.TRUE., acc_devalloc=.TRUE.,acc_stream=stream_2, error=error)
       CALL dbcsr_mempool_ensure_capacity(memtype_trsbuffer_2%pool, capacity=1)
    ELSE
       CALL dbcsr_memtype_setup(memtype_abpanel_1, mpi=.TRUE., error=error)
       CALL dbcsr_memtype_setup(memtype_abpanel_2, mpi=.TRUE., error=error)
    ENDIF

    CALL dbcsr_memtype_setup(memtype_mpi_buffer, mpi=.TRUE., error=error)

    ! check parameters ---------------------------------------------------------
    transa_l = transa
    transb_l = transb
    CALL uppercase(transa_l)
    CALL uppercase(transb_l)
    CALL dbcsr_assert(transa_l.EQ.dbcsr_no_transpose.OR.&
                      transa_l.EQ.dbcsr_transpose.OR.&
                      transa_l.EQ.dbcsr_conjugate_transpose,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Invalid transa_l = "//transa_l, __LINE__, error)
    CALL dbcsr_assert(transb_l.EQ.dbcsr_no_transpose.OR.&
                      transb_l.EQ.dbcsr_transpose.OR.&
                      transb_l.EQ.dbcsr_conjugate_transpose,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Invalid transb_l = "//transb_l, __LINE__, error)

    IF (dbg) THEN
       WRITE(*,*)'========== MULTIPLICATION ========================'
       CALL dbcsr_verify_matrix (matrix_a, error=error)
       CALL dbcsr_verify_matrix (matrix_b, error=error)
       CALL dbcsr_verify_matrix (matrix_c, error=error)
       WRITE(*,*)routineN//" ABC checksums",&
            dbcsr_checksum(matrix_a, error=error),&
            dbcsr_checksum(matrix_b, error=error),&
            dbcsr_checksum(matrix_c, error=error)
       IF (dbg) THEN
          CALL dbcsr_print (matrix_a, nodata=.TRUE., error=error)
          CALL dbcsr_print (matrix_b, nodata=.TRUE., error=error)
          CALL dbcsr_print (matrix_c, nodata=.TRUE., error=error)
       ENDIF
    ENDIF

    ! transpose/conjg left and/or right matrices if needed
    SELECT CASE(transa_l)
    CASE(dbcsr_no_transpose)
       matrix_left = matrix_a
       new_left = .FALSE.
    CASE(dbcsr_transpose)
       CALL dbcsr_init(matrix_left)
       IF(dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric) THEN
          !
          ! For antisymmetric matrix, we need to do a hard copy
          ! shallow_data_copy=.TRUE. doesnt handle properly antisymm matrices
          CALL dbcsr_new_transposed (matrix_left, matrix_a,&
               shallow_data_copy=.FALSE., redistribute=.FALSE., &
               transpose_distribution=.FALSE., error=error)
       ELSE
          CALL dbcsr_new_transposed (matrix_left, matrix_a,&
               shallow_data_copy=.TRUE., redistribute=.FALSE.,&
               transpose_distribution=.FALSE., error=error)
       ENDIF
       new_left = .TRUE.
    CASE(dbcsr_conjugate_transpose)
       CALL dbcsr_init(matrix_left)
       CALL dbcsr_new_transposed (matrix_left, matrix_a,&
            shallow_data_copy=.FALSE., redistribute=.FALSE.,&
            transpose_distribution=.FALSE., error=error)
       CALL dbcsr_conjg(matrix_left, error=error)
       new_left = .TRUE.
    CASE DEFAULT
       CALL dbcsr_assert(.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, &
            routineN, "wrong transa_l = "//transa_l, __LINE__, error)
    END SELECT

    SELECT CASE(transb_l)
    CASE(dbcsr_no_transpose)
       matrix_right = matrix_b
       new_right = .FALSE.
    CASE(dbcsr_transpose)
       CALL dbcsr_init(matrix_right)
       IF(dbcsr_get_matrix_type(matrix_b).EQ.dbcsr_type_antisymmetric) THEN
          !
          ! For antisymmetric matrix, we need to do a hard copy
          ! shallow_data_copy=.TRUE. doesnt handle properly antisymm matrices
          CALL dbcsr_new_transposed (matrix_right, matrix_b,&
               shallow_data_copy=.FALSE., redistribute=.FALSE.,&
               transpose_distribution=.FALSE., error=error)
       ELSE
          CALL dbcsr_new_transposed (matrix_right, matrix_b,&
               shallow_data_copy=.TRUE., redistribute=.FALSE.,&
               transpose_distribution=.FALSE., error=error)
       ENDIF
       new_right = .TRUE.
    CASE(dbcsr_conjugate_transpose)
       CALL dbcsr_init(matrix_right)
       CALL dbcsr_new_transposed (matrix_right, matrix_b,&
            shallow_data_copy=.FALSE., redistribute=.FALSE.,&
            transpose_distribution=.FALSE., error=error)
       CALL dbcsr_conjg(matrix_right, error=error)
       new_right = .TRUE.
    CASE DEFAULT
       CALL dbcsr_assert(.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, &
            routineN, "wrong transb_l = "//transb_l, __LINE__, error)
    END SELECT

    !
    ! Ensure matrix compatibility.
    CALL dbcsr_assert (array_equality (dbcsr_row_block_offsets (matrix_c),&
                                       dbcsr_row_block_offsets (matrix_left)),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "C/A rows not equal", __LINE__, error=error)
    CALL dbcsr_assert (array_equality (dbcsr_col_block_offsets (matrix_c),&
                                       dbcsr_col_block_offsets (matrix_right)),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "C/B columns not equal", __LINE__, error=error)
    CALL dbcsr_assert (array_equality (dbcsr_col_block_offsets (matrix_left),&
                                       dbcsr_row_block_offsets (matrix_right)),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "A cols/B rows not equal", __LINE__, error=error)

    !
    ! No dense multiplication when filtering is used.
    use_dense_mult = .NOT. PRESENT (filter_eps)
    IF (mm_async) use_dense_mult = .FALSE.
    ! we skip dense multiply for (anti)symmetric matrices (slowdown for S/H * C)
    IF (use_dense_mult) THEN
       IF(dbcsr_has_symmetry (matrix_left) .OR. &
            dbcsr_has_symmetry(matrix_right)) THEN
          use_dense_mult = .FALSE.
       ELSE
          use_dense_mult = dbcsr_may_be_dense (matrix_left, make_dense_occ_thresh)&
               .AND. dbcsr_may_be_dense (matrix_right, make_dense_occ_thresh)
       ENDIF
    ENDIF
    ab_dense = use_dense_mult
    !
    ! Submatrix selection
    f_row = 1
    l_row = dbcsr_nfullrows_total(matrix_c)
    f_col = 1
    l_col = dbcsr_nfullcols_total(matrix_c)
    f_k = 1
    l_k = dbcsr_nfullcols_total(matrix_left)
    IF (PRESENT (first_row)) THEN
       CALL dbcsr_assert(first_row .GE. 1&
            .AND. first_row .LE. dbcsr_nfullrows_total(matrix_c),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid first row specified", __LINE__, error)
       f_row = first_row
    ENDIF
    IF (PRESENT (last_row)) THEN
       CALL dbcsr_assert(last_row .GE. 1&
            .AND. last_row .LE. dbcsr_nfullrows_total(matrix_c)&
            .OR. last_row .LT. 1 ,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid last row specified", __LINE__, error)
       l_row = last_row
    ENDIF
    IF (PRESENT (first_column)) THEN
       CALL dbcsr_assert(first_column .GE. 1&
            .AND. first_column .LE. dbcsr_nfullcols_total(matrix_c),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid first col specified", __LINE__, error)
       f_col = first_column
    ENDIF
    IF (PRESENT (last_column)) THEN
       CALL dbcsr_assert(last_column .GE. 1&
            .AND. last_column .LE. dbcsr_nfullcols_total(matrix_c)&
            .OR. last_column .LT. 1,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid last column specified (C)", __LINE__, error)
       CALL dbcsr_assert(last_column .GE. 1&
            .AND. last_column .LE. dbcsr_nfullcols_total(matrix_right)&
            .OR. last_column .LT. 1,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid last column specified (B)", __LINE__, error)
       l_col = last_column
    ENDIF
    IF (PRESENT (first_k)) THEN
       CALL dbcsr_assert(first_k .GE. 1&
            .AND. first_k .LE. dbcsr_nfullcols_total(matrix_left),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid first k specified (A)", __LINE__, error)
       f_k = first_k
    ENDIF
    IF (PRESENT (last_k)) THEN
       CALL dbcsr_assert(last_k.GE. 1&
            .AND. last_k .LE. dbcsr_nfullcols_total(matrix_left)&
            .OR. last_k  .LT. 1,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid last k specified (A)", __LINE__, error)
       l_k = last_k
    ENDIF

    ! update statistics (we count marketing flops per MPI rank), 
    ! note that this number can easily overflow
    numnodes = dbcsr_mp_numnodes (dbcsr_distribution_mp (dbcsr_distribution (matrix_c)))
    last_mpi_ranks_used = numnodes
    marketing_flops = marketing_flops + &
                         (2_int_8*INT(l_row-f_row+1, int_8)*INT(l_col-f_col+1, int_8)/numnodes) * &
                         INT(l_k-f_k+1,     int_8)

    !
    ! Now optimize the default submatrix selection values away
    IF (f_row .EQ. 1) f_row = 0
    IF (l_row .EQ. dbcsr_nfullrows_total(matrix_left)) l_row = 0
    IF (f_col .EQ. 1) f_col = 0
    ! The last column must be set if the right and product matrices
    ! differ.
    l_col = MIN (l_col, dbcsr_nfullcols_total(matrix_right))
    l_col = MIN (l_col, dbcsr_nfullcols_total(matrix_c))
    IF (f_col.LE.1.AND.&
        l_col .EQ. dbcsr_nfullcols_total(matrix_right) .AND. &
        dbcsr_nfullcols_total(matrix_right) .EQ.&
        dbcsr_nfullcols_total(matrix_c)) l_col = 0
    IF (f_k .EQ. 1) f_k = 0
    IF (l_k .EQ. dbcsr_nfullcols_total(matrix_left)) l_k = 0
    IF (.NOT. PRESENT(last_column) .AND.&
        .NOT. array_equality (dbcsr_col_block_sizes (matrix_right),&
                              dbcsr_col_block_sizes (matrix_c))) THEN
       l_col = MIN (dbcsr_nfullcols_total(matrix_right),&
                    dbcsr_nfullcols_total(matrix_c))
    ENDIF
    CALL dbcsr_assert (f_row .LE. l_row, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Last row smaller than first row", &
         __LINE__, error)
    CALL dbcsr_assert (f_col .LE. l_col, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Last col smaller than first col", &
         __LINE__, error)

    !
    ! if we have limits we need to turn of make dense for the moment...
    !IF(ANY((/ f_row, l_row, f_col, l_col, f_k, l_k /).NE.0)) use_dense_mult = .FALSE.

    !
    !
    ! Product data needs to be retained when
    ! * beta != 0; or
    ! * there is column limiting (l_col > 0) and the limiting column
    !   is less than the number of full columns in theproduct matrix
    keep_sparsity = .FALSE.
    IF (PRESENT (retain_sparsity)) keep_sparsity=retain_sparsity
    !
    keep_product_data = keep_sparsity&
         .OR. .NOT. dbcsr_scalar_are_equal (beta, dbcsr_scalar_zero(beta%data_type))&
         .OR. (l_col .GT. 0 .AND. l_col .LT. dbcsr_nfullcols_total(matrix_c)) &
         .OR. (l_row .GT. 0 .AND. l_row .LT. dbcsr_nfullrows_total(matrix_c))
    !
    IF (.NOT. dbcsr_scalar_are_equal (beta, dbcsr_scalar_one(beta%data_type)) .AND. keep_product_data) THEN
       CALL dbcsr_scale (matrix_c, alpha_scalar=beta, &
            limits=(/f_row,l_row,f_col,l_col/), error=error)
    ENDIF
    !
    ! The index of the product matrix is twiddled into canonical form
    ! if it is (anti)symmetric (i.e., rows and columns are where the
    ! row/column distributions say they are). Doing this in advance
    ! makes the local multiply more efficient.
    IF (dbcsr_has_symmetry (matrix_c)) THEN
       product_reindex = .TRUE.
    ELSE
       product_reindex = .FALSE.
    ENDIF
    ! Product can not be made dense; however, A & B may still be made
    ! dense unless previously determined otherwise.
    IF (product_reindex.OR.keep_sparsity) THEN
       use_dense_mult = .FALSE.
    ENDIF
    !
    ! The thread distribution must reflect the current (possibly
    ! dense) distribution
    !CALL dbcsr_assert (dbcsr_distribution_has_threads(product_matrix%m%dist),&
    !     dbcsr_fatal_level, dbcsr_internal_error, routineN,&
    !     "Thread distribution not defined.", __LINE__, error=error)
    IF (.NOT. dbcsr_distribution_has_threads(matrix_c%m%dist)) THEN
       release_tdist = .TRUE.
       CALL dbcsr_distribution_make_threads (matrix_c%m%dist)
    ELSE
       release_tdist = .FALSE.
    ENDIF

    mult_images = 1
    IF (dbcsr_get_conf_avg_elements_images().GT.0) THEN
       mult_images = CEILING(DBLE(MAX(dbcsr_distribution_data_max_size_allproc(matrix_left),&
                                      dbcsr_distribution_data_max_size_allproc(matrix_right)))/&
                             dbcsr_get_conf_avg_elements_images())
    ENDIF
    !
    ! Create imaged distributions for the multiply.
    CALL dbcsr_create_image_dist (rdist_right, matrix_right%m%dist,&
         match_row_nbins = dbcsr_mp_npcols (dbcsr_distribution_mp (matrix_left%m%dist)),&
         match_col_nbins = dbcsr_mp_npcols (dbcsr_distribution_mp (matrix_c%m%dist)),&
         match_col_pdist = dbcsr_distribution_col_dist (matrix_c%m%dist),&
         mult_images=mult_images)
    CALL dbcsr_create_image_dist (rdist_left, matrix_left%m%dist,&
         match_row_pdist = dbcsr_distribution_row_dist (matrix_c%m%dist),&
         match_row_nbins = dbcsr_mp_nprows (dbcsr_distribution_mp (matrix_c%m%dist)),&
         match_col_pdist = dbcsr_distribution_row_dist (rdist_right%i%main),&
         match_col_idist = array_data (rdist_right%i%row_image),&
         match_col_nbins = dbcsr_mp_nprows (dbcsr_distribution_mp(matrix_right%m%dist)),&
         mult_images=mult_images)
    IF (ab_dense) THEN
       CALL dbcsr_make_dists_dense (dbcsr_distribution (matrix_c),&
            rdist_left, rdist_right, dense_product_distribution,&
            dense_rdist_left, dense_rdist_right, .not.use_dense_mult,&
            m_map, k_vmap, n_map, matrix_c%m%row_blk_size, error=error)
       CALL make_sizes_dense (matrix_c%m%row_blk_size, m_map,&
            dbcsr_distribution_nrows (dense_product_distribution),&
            dense_row_sizes,&
            error=error)
       CALL make_sizes_dense (matrix_c%m%col_blk_size, n_map, &
            dbcsr_distribution_ncols (dense_product_distribution),&
            dense_col_sizes,&
            error=error)
       CALL make_sizes_dense (matrix_right%m%row_blk_size, k_vmap,&
            dbcsr_distribution_nrows (dense_rdist_right%i%main),&
            dense_k_sizes,&
            error=error)
    ENDIF
    !
    CALL dbcsr_assert (use_dense_mult, "IMP", ab_dense,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Wrong logic when making dense matrices.", __LINE__, error=error)
    IF (use_dense_mult) THEN
       old_product_row_blk_offsets = matrix_c%m%row_blk_offset
       old_product_col_blk_offsets = matrix_c%m%col_blk_offset
       old_product_row_blk_sizes =   matrix_c%m%row_blk_size
       old_product_col_blk_sizes =   matrix_c%m%col_blk_size
       CALL array_hold (old_product_row_blk_offsets)
       CALL array_hold (old_product_col_blk_offsets)
       CALL array_hold (old_product_row_blk_sizes)
       CALL array_hold (old_product_col_blk_sizes)
       old_product_distribution = dbcsr_distribution (matrix_c)
       CALL dbcsr_distribution_hold (old_product_distribution)
       CALL dbcsr_init (product_matrix)
       CALL dbcsr_make_dense (matrix_c, product_matrix,&
            dense_product_distribution,&
            dense_row_sizes, dense_col_sizes,&
            m_map, n_map,&
            error=error)
    ELSE
       CALL dbcsr_init (product_matrix)
       CALL dbcsr_copy(product_matrix, matrix_c, shallow_data=.TRUE., error=error)
    ENDIF
    IF (ab_dense) THEN
       CALL dbcsr_distribution_release (dense_product_distribution)
    ENDIF
    !

    ! This is needed to build the hash tables because they are
    ! locally indexed.
    CALL dbcsr_reset_locals (product_matrix, error=error)

    IF (debug_mod) THEN
       WRITE(*,*)routineN//" Matrices ", dbcsr_get_matrix_type(matrix_a),&
            dbcsr_get_matrix_type(matrix_b), dbcsr_get_matrix_type(matrix_c)
       WRITE(*,*)routineN//" Matrices ", transa_l, transb_l, "keep", keep_product_data
    ENDIF
    IF (keep_product_data) THEN
       IF (product_reindex) THEN
          IF (debug_mod) WRITE(*,*)routineN//" Making canonical index"
          CALL dbcsr_make_index_canonical (product_matrix)
       ENDIF
       CALL dbcsr_assert (.NOT. ASSOCIATED (product_matrix%m%wms),&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Product matrix should be finalized!", __LINE__, error=error)
       CALL dbcsr_make_untransposed_blocks (product_matrix, error=error)
!$omp parallel &
!$omp default (none) shared (product_matrix, error)
       ! For the multiply logic to work correctly, existing data must
       ! be added only after the index has been transformed into the
       ! canonical form.
       CALL dbcsr_add_wm_from_matrix(product_matrix, error=error)
!$omp end parallel
    ELSE
!$omp parallel default(none) private(ithread) &
!$omp shared(product_matrix, memtype_product_wm, error)
       ithread = 0
       !$ ithread = OMP_GET_THREAD_NUM ()
       CALL dbcsr_work_create(product_matrix, work_mutable=.FALSE., &
               memory_type=memtype_product_wm(ithread)%p, error=error)
!$omp end parallel
    ENDIF
    product_matrix%m%nze = 0
    product_matrix%m%row_p(:) = 0
    CALL dbcsr_data_set_size_referenced(product_matrix%m%data_area, 0)
    product_matrix%m%nblks = 0
    product_matrix%m%valid = .FALSE.
    !
    ! Right images
    CALL make_m2s(matrix_right,m2s_right,rdist_right,dense_rdist_right,&
                  use_dense_mult,ab_dense,use_local_indexing,use_list_indexing,&
                  "R",f_k,l_k,f_row,l_row,f_col,l_col,&
                  dense_k_sizes,dense_col_sizes,&
                  k_vmap,m_map,n_map,&
                  alpha,error)
    !
    ! Left images
    CALL make_m2s(matrix_left,m2s_left,rdist_left,dense_rdist_left,&
                  use_dense_mult,ab_dense,use_local_indexing,use_list_indexing,&
                  "L",f_k,l_k,f_row,l_row,f_col,l_col,&
                  dense_row_sizes,dense_k_sizes,&
                  k_vmap,m_map,n_map,&
                  error=error)
    !
    IF (ab_dense) THEN
       CALL array_release (k_vmap)
       CALL array_release (dense_row_sizes)
       CALL array_release (dense_col_sizes)
       CALL array_release (dense_k_sizes)
    ENDIF
    !
    ! The limits were already used.  Reset them.
    f_row = 0 ; l_row = 0
    f_col = 0 ; l_col = 0
    f_k = 0 ; l_k = 0
    !
    my_flop = 0
    IF (use_mpi_rma) THEN
       CALL cannon_multiply_low_rma(m2s_left, m2s_right, product_matrix,&
            retain_sparsity=retain_sparsity,&
            filter_eps=filter_eps, error=error,&
            flop=my_flop)
    ELSE
       CALL cannon_multiply_low(m2s_left, m2s_right, product_matrix,&
            retain_sparsity=retain_sparsity,&
            filter_eps=filter_eps, error=error,&
            flop=my_flop)
    ENDIF

    CALL dbcsr_finalize(product_matrix, error=error)

    IF (PRESENT (flop)) THEN
       ! return the average number of flops per MPI rank. Variance (which is fairly large) could be computed as well.
       CALL dbcsr_error_set(routineN//"_mpsum_flop",error_handler2, error)
       comm = dbcsr_mp_group (dbcsr_distribution_mp (dbcsr_distribution (product_matrix)))
       numnodes = dbcsr_mp_numnodes (dbcsr_distribution_mp (dbcsr_distribution (product_matrix)))
       CALL mp_sum(my_flop,comm)
       flop = (my_flop + numnodes - 1) / numnodes
       CALL dbcsr_error_stop(error_handler2, error)
    ENDIF
    !
    IF (new_left) CALL dbcsr_release (matrix_left)
    IF (new_right) CALL dbcsr_release (matrix_right)
    IF (release_tdist) THEN
       CALL dbcsr_distribution_no_threads (product_matrix%m%dist)
    ENDIF
    !
    IF (.TRUE. .OR. use_local_indexing) &
         CALL dbcsr_release_locals (product_matrix, error=error)
    ! The index of the product matrix is reset to the CP2K form if it
    ! was previously set to the canonical form.
    IF (product_reindex) THEN
       IF (debug_mod) WRITE(*,*)routineN//" Making CP2K index"
       CALL dbcsr_make_index_canonical(product_matrix, cp2k=.TRUE.)
    ENDIF
    IF (use_dense_mult) THEN
       CALL dbcsr_release (matrix_c)
       CALL dbcsr_init (matrix_c)
       CALL dbcsr_make_undense(product_matrix, matrix_c,&
            old_product_distribution,&
            old_product_row_blk_offsets, old_product_col_blk_offsets,&
            old_product_row_blk_sizes, old_product_col_blk_sizes,&
            m_map, n_map, error=error)
       CALL dbcsr_release (product_matrix)
       CALL array_release (old_product_row_blk_offsets)
       CALL array_release (old_product_col_blk_offsets)
       CALL array_release (old_product_row_blk_sizes)
       CALL array_release (old_product_col_blk_sizes)
       CALL dbcsr_distribution_release (old_product_distribution)
    ELSE
       CALL dbcsr_release (matrix_c)
       CALL dbcsr_init (matrix_c)
       CALL dbcsr_copy (matrix_c, product_matrix, shallow_data=.TRUE., error=error)
       CALL dbcsr_release (product_matrix)
    ENDIF
    !

    CALL dbcsr_destroy_array (m2s_left, error=error)
    DEALLOCATE (m2s_left)

    CALL dbcsr_image_dist_release (rdist_left, error=error)
    CALL dbcsr_destroy_array (m2s_right, error=error)
    DEALLOCATE (m2s_right)
    CALL dbcsr_image_dist_release (rdist_right, error=error)
    IF (ab_dense) THEN
       CALL array_release (m_map)
       CALL array_release (n_map)
    ENDIF
    !
    ! if filtering is requested remove small blocks, unless the sparsity needs to be kept.
    !
    IF (PRESENT (filter_eps) .AND. .NOT. keep_sparsity) THEN
       eps_any = dbcsr_scalar(filter_eps)
       CALL dbcsr_scalar_fill_all(eps_any)
       CALL dbcsr_scalar_set_type(eps_any, dbcsr_get_data_type(matrix_c))
       CALL dbcsr_filter (matrix_c, eps_any, quick=.FALSE., error=error)
    ENDIF
    !
    ! To support the canonical multiply (all non-transposed blocks),
    ! blocks may have to be transposed according to the CP2K
    ! triangular index.
    CALL dbcsr_make_untransposed_blocks (matrix_c, error=error)
    !
    IF (debug_mod .OR. careful_mod) THEN
       IF (debug_mod) &
            WRITE(*,*)routineN//" Use dense mult, symm",&
            use_dense_mult, ab_dense, dbcsr_has_symmetry (matrix_c)
       CALL dbcsr_verify_matrix (matrix_c, error=error)
       IF (debug_mod) THEN
          cs = dbcsr_checksum (matrix_c, error=error)
          WRITE(*,*)routineN//" Product checksum", cs
       ENDIF
    ENDIF
    !
    IF (.FALSE.) WRITE(*,*)"Finished with one multiplication."
    output_unit = default_output_unit
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_mm_cannon_multiply



! *****************************************************************************
!> \brief Make images from the matrix (left or right)
!> \param[in] matrix ...
!> \param[out] m2s ...
!> \param[in,out] rdist ...
!> \param[in,out] dense_rdist ...
!> \param[in] use_dense_mult ...
!> \param[in] ab_dense ...
!> \param[in] use_local_indexing ...
!> \param[in] use_list_indexing ...
!> \param[in] predistribute ...
!> \param[in] f_k ...
!> \param[in] l_k ...
!> \param[in] f_row ...
!> \param[in] l_row ...
!> \param[in] f_col ...
!> \param[in] l_col ...
!> \param[in,out] row_blk_size ...
!> \param[in,out] col_blk_size ...
!> \param[in] k_vmap ...
!> \param[in] m_map ...
!> \param[in] n_map ...
!> \param[in] alpha ...
!> \param[in,out] error ...
! *****************************************************************************
  SUBROUTINE make_m2s(matrix,m2s,rdist,dense_rdist,&
                      use_dense_mult,ab_dense,&
                      use_local_indexing,use_list_indexing,&
                      predistribute,f_k,l_k,f_row,l_row,f_col,l_col,&
                      row_blk_size,col_blk_size,&
                      k_vmap,m_map,n_map,&
                      alpha,error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_2d_array_type), INTENT(OUT), &
      POINTER                                :: m2s
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(INOUT)                          :: rdist, dense_rdist
    LOGICAL, INTENT(IN)                      :: use_dense_mult, ab_dense, &
                                                use_local_indexing, &
                                                use_list_indexing
    CHARACTER, INTENT(IN)                    :: predistribute
    INTEGER, INTENT(IN)                      :: f_k, l_k, f_row, l_row, &
                                                f_col, l_col
    TYPE(array_i1d_obj), INTENT(INOUT)       :: row_blk_size, col_blk_size
    TYPE(array_i1d_obj), INTENT(IN)          :: k_vmap, m_map, n_map
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: alpha
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i, im, j, jm
    INTEGER, DIMENSION(4)                    :: f_crop
    LOGICAL                                  :: do_scale, thread_redist
    TYPE(array_i1d_obj)                      :: col_map, row_map
    TYPE(dbcsr_obj)                          :: dense_template, matrix_tmp

    ALLOCATE (m2s)
    do_scale = .FALSE.
    IF (PRESENT(alpha)) THEN
       IF (.NOT. dbcsr_scalar_are_equal (alpha, dbcsr_scalar_one(alpha%data_type))) THEN
          do_scale = .TRUE.
       END IF
    END IF

    IF (do_scale) THEN
       ! Copy and scale matrix if alpha is not 1.
       CALL dbcsr_make_images (matrix, m2s, rdist,&
            predistribute=predistribute, &
            data_memory_type = memtype_abpanel_1,&
            index_memory_type = dbcsr_memtype_default,&
            no_copy_data=use_dense_mult, scale_value=alpha, error=error)
    ELSE
       CALL dbcsr_make_images (matrix, m2s, rdist,&
            predistribute=predistribute, &
            data_memory_type = memtype_abpanel_1,&
            index_memory_type = dbcsr_memtype_default,&
            no_copy_data=use_dense_mult, error=error)
    END IF
    
    SELECT CASE (predistribute)
    CASE('L')
       im = 1
       jm = SIZE (m2s%mats,2)
       f_crop = (/ f_row, l_row, f_k, l_k /)
       row_map = m_map
       col_map = k_vmap
       thread_redist = .TRUE.
    CASE default
       im = SIZE (m2s%mats,1)
       jm = 1
       f_crop = (/ f_k, l_k, f_col, l_col /)
       row_map = k_vmap
       col_map = n_map
       thread_redist = .FALSE.
    END SELECT

    ! Post-processing of images.
    DO i = 1, im
       DO j = 1, jm
          CALL dbcsr_reset_vlocals (m2s%mats(i,j), rdist, error=error)
          ! Crop if necessary
          IF (ANY(f_crop .NE. 0)) THEN
             CALL dbcsr_init (matrix_tmp)
             CALL dbcsr_crop_matrix (matrix_tmp, m2s%mats(i,j),&
                  full_row_bounds=f_crop(1:2),&
                  full_column_bounds=f_crop(3:4),&
                  shallow_data = .FALSE., error=error)
             CALL dbcsr_release (m2s%mats(i,j))
             CALL dbcsr_copy (m2s%mats(i,j), matrix_tmp, shallow_data=.TRUE.,&
                  error=error)
             CALL dbcsr_release (matrix_tmp)
             CALL dbcsr_reset_vlocals (m2s%mats(i,j), rdist, error=error)
          ENDIF
       ENDDO
    ENDDO

    IF (ab_dense) THEN
       CALL dbcsr_init (dense_template)
       CALL dbcsr_create (dense_template, template=matrix,&
            dist=dense_rdist%i%main,&
            row_blk_size=array_data(row_blk_size), col_blk_size=array_data(col_blk_size),&
            error=error)
       CALL dbcsr_make_images_dense (m2s, dense_rdist, &
            row_map = row_map, col_map = col_map,&
            join_cols = use_dense_mult, join_rows=ab_dense, &
            new_template=dense_template, error=error)
       CALL dbcsr_image_dist_release (rdist, error=error)
       rdist = dense_rdist
       CALL dbcsr_image_dist_hold (rdist, error=error)
       DO i = 1, im
          DO j = 1, jm
             CALL dbcsr_reset_vlocals (m2s%mats(i,j), rdist, error=error)
          ENDDO
       ENDDO
    ENDIF

    DO i = 1, im
       DO j = 1, jm
          ! Convert to local-row index
          IF (use_local_indexing) THEN
             CALL dbcsr_make_index_local_row(m2s%mats(i,j), error=error)
          ENDIF
          ! Convert to list index
          IF (use_list_indexing) THEN
             CALL dbcsr_make_index_list(m2s%mats(i,j), thread_redist=thread_redist,&
                  error=error)
          ENDIF
          IF (use_local_indexing .AND. .NOT. use_list_indexing) THEN
             CALL dbcsr_index_compact(m2s%mats(i,j), error=error)
          ENDIF
       ENDDO
    ENDDO

    IF (ab_dense) THEN
       CALL dbcsr_image_dist_release (dense_rdist, error=error)
       CALL dbcsr_release (dense_template)
    ENDIF

  END SUBROUTINE make_m2s


! *****************************************************************************
!> \brief Creates row and column images of a matrix.
!> \param[in] source          input matrix
!> \param[in,out] normalized  image array of the normalized matrix
!> \param[in] target_image_dist          normalize to this image distribution
!> \param[in] predistribute   (optional) predistribute data for multiplication
!> \param[in] data_memory_type     type of memory to use for data
!> \param[in] index_memory_type    type of memory to use for index
!> \param[in] no_copy_data    (optional) try to not merge data at the end
!> \param[in] scale_value     (optional) scale with this value
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_make_images(source, normalized, target_image_dist,&
       predistribute, data_memory_type, index_memory_type,&
       no_copy_data, scale_value, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: source
    TYPE(dbcsr_2d_array_type), INTENT(OUT)   :: normalized
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(IN)                             :: target_image_dist
    CHARACTER, INTENT(IN), OPTIONAL          :: predistribute
    TYPE(dbcsr_memtype_type), INTENT(IN)     :: data_memory_type, &
                                                index_memory_type
    LOGICAL, INTENT(IN), OPTIONAL            :: no_copy_data
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: scale_value
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_images', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_valid_index (source),&
         dbcsr_fatal_level, dbcsr_caller_error, routineN,&
         "Matrix not initialized.",__LINE__,error)
    CALL make_images(source, normalized,&
         target_image_dist, desymmetrize=dbcsr_has_symmetry(source),&
         predistribute=predistribute,&
         data_memory_type = data_memory_type,&
         index_memory_type = index_memory_type,&
         no_copy_data=no_copy_data,&
         scale_value=scale_value,&
         error=error)
    normalized%image_dist = target_image_dist
    CALL dbcsr_image_dist_hold (normalized%image_dist, error=error)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_make_images

! *****************************************************************************
!> \brief Makes column-based and row-based images of a matrix.
!> \param[in] ism input symmetric matrix
!> \param[in,out] ums         normalied matrices
!> \param[in] target_imgdist  image distribution to normalize to
!> \param[in] desymmetrize    (optional) desymmetrize a symmetric matrix
!> \param[in] predistribute   (optional) predistribute data for multiplication
!> \param[in] data_memory_type     type of memory to use for data
!> \param[in] index_memory_type    type of memory to use for index
!> \param[in] no_copy_data    (optional) try to not merge data at the end
!> \param[in] scale_value     (optional) scale with this value
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE make_images(ism, ums, target_imgdist, desymmetrize, predistribute,&
       data_memory_type, index_memory_type, no_copy_data, scale_value, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: ism
    TYPE(dbcsr_2d_array_type), INTENT(OUT)   :: ums
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(IN)                             :: target_imgdist
    LOGICAL, INTENT(IN), OPTIONAL            :: desymmetrize
    CHARACTER, INTENT(IN), OPTIONAL          :: predistribute
    TYPE(dbcsr_memtype_type), INTENT(IN)     :: data_memory_type, &
                                                index_memory_type
    LOGICAL, INTENT(IN), OPTIONAL            :: no_copy_data
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: scale_value
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_images', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: metalen = 5
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    CHARACTER                                :: predist_type, predist_type_fwd
    INTEGER :: blk, blk_l, blk_p, bp, col, col_img, col_size, coli, data_p, &
      data_type, dst_p, error_handler, error_handler2, ithread, mp_group, &
      ncol_images, nrow_images, nsymmetries, nthreads, numproc, nze, pcol, &
      prev_blk_p, prow, row, row_img, row_size, rowi, sd_pos, sm_pos, src_p, &
      stored_blk_p, stored_col, stored_row, symmetry_i, tr_col_size, &
      tr_row_size, vcol, vrow
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: myt_sdp, myt_smp, rd_disp, &
                                                recv_meta, rm_disp, sd_disp, &
                                                sdp, send_meta, sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: all_total_send_offset, blk_ps, &
      blks, myt_total_send_count, total_recv_count, total_send_count
    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: myt_send_count, recv_count, &
                                                send_count
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, col_dist, &
                                                col_img_dist, row_blk_size, &
                                                row_dist, row_img_dist
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    LOGICAL                                  :: nocopy, tr
    TYPE(dbcsr_data_obj)                     :: received_data_area, &
                                                recv_data_area, send_data_area
    TYPE(dbcsr_distribution_obj)             :: old_dist, target_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error, t_error
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_scalar_type)                  :: scale_neg_one
    TYPE(dbcsr_type)                         :: sm

!   ---------------------------------------------------------------------------
! Check input matrix
! Set convenient variables to access input matrix info
!

    CALL dbcsr_error_set(routineN, error_handler, dbcsr_error)
    nocopy = .FALSE.
    IF (PRESENT (no_copy_data)) nocopy = no_copy_data
    sm = ism%m
    nsymmetries = 1
    IF (PRESENT (desymmetrize)) THEN
       IF (desymmetrize .AND. sm%symmetry) THEN
          nsymmetries = 2
       ENDIF
    ENDIF
    SELECT CASE (predistribute)
    CASE('L','l')
       predist_type = 'L'
       predist_type_fwd = 'l'
    CASE('R','r')
       predist_type = 'R'
       predist_type_fwd = 'r'
    CASE default
       CALL dbcsr_assert (.FALSE.,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Incorrect pre-shift specifier.", __LINE__, error=error)
    END SELECT
    data_type = sm%data_type
    CALL dbcsr_assert (data_type .EQ. dbcsr_type_real_8&
         .or.data_type .EQ. dbcsr_type_real_4&
         .or.data_type .EQ. dbcsr_type_complex_8&
         .or.data_type .EQ. dbcsr_type_complex_4, dbcsr_fatal_level,&
         dbcsr_internal_error, routineN, "Invalid data type.",__LINE__,error)
    scale_neg_one = dbcsr_scalar_negative (dbcsr_scalar_one (data_type))
    row_blk_size => array_data (sm%row_blk_size)
    col_blk_size => array_data (sm%col_blk_size)
    old_dist = dbcsr_distribution (ism)
    target_dist = target_imgdist%i%main
    row_dist => dbcsr_distribution_row_dist (target_dist)
    col_dist => dbcsr_distribution_col_dist (target_dist)
    IF (sm%symmetry) THEN
       CALL dbcsr_assert(SIZE(row_dist),'EQ', SIZE(col_dist), dbcsr_warning_level,&
            dbcsr_wrong_args_error, routineN,&
            'Unequal row and column distributions for symmetric matrix.',__LINE__,error)
    ENDIF
    nrow_images = target_imgdist%i%row_decimation
    IF (nrow_images .GT. 1) THEN
       row_img_dist => array_data (target_imgdist%i%row_image)
    ELSE
       NULLIFY (row_img_dist)
    ENDIF
    ncol_images = target_imgdist%i%col_decimation
    IF (ncol_images .GT. 1) THEN
       col_img_dist => array_data (target_imgdist%i%col_image)
    ELSE
       NULLIFY (col_img_dist)
    ENDIF
    mp_obj = dbcsr_distribution_mp (target_dist)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    numproc = dbcsr_mp_numnodes (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    CALL dbcsr_assert(MAXVAL(row_dist),'LE',UBOUND(blacs2mpi,1), dbcsr_failure_level,&
         dbcsr_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows',__LINE__,error)
    IF (dbg) &
         CALL dbcsr_assert(MAXVAL(row_dist),'EQ',UBOUND(blacs2mpi,1), dbcsr_warning_level,&
         dbcsr_wrong_args_error, routineN,&
         'Range of row distribution not equal to processor rows',__LINE__,error)
    CALL dbcsr_assert(MAXVAL(col_dist),'LE',UBOUND(blacs2mpi,2), dbcsr_failure_level,&
         dbcsr_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols',__LINE__,error)
    IF (dbg) &
         CALL dbcsr_assert(MAXVAL(col_dist),'EQ',UBOUND(blacs2mpi,2), dbcsr_warning_level,&
         dbcsr_wrong_args_error, routineN,&
         'Range of col distribution not equal to processor cols',__LINE__,error)

    ! Check threads configuration
!$  CALL dbcsr_assert (dbcsr_distribution_has_threads (old_dist),&
!$       dbcsr_fatal_level, dbcsr_internal_error, routineN,&
!$       "Thread distribution not defined", __LINE__, error=error)

    ! Allocate shared temporary buffers
    !
    ALLOCATE (send_count(2, nrow_images, ncol_images, 0:numproc-1)) ; send_count = 0
    ALLOCATE (recv_count(2, nrow_images, ncol_images, 0:numproc-1))
    ALLOCATE (total_send_count(2, 0:numproc-1)) ; total_send_count = 0
    ALLOCATE (total_recv_count(2, 0:numproc-1))
    ALLOCATE (sdp(0:numproc-1))
    ALLOCATE (smp(0:numproc-1))
    ALLOCATE (sd_disp(0:numproc-1)) ; sd_disp(0) = 1
    ALLOCATE (sm_disp(0:numproc-1)) ; sm_disp(0) = 1
    ALLOCATE (rd_disp(0:numproc-1)) ; rd_disp(0) = 1
    ALLOCATE (rm_disp(0:numproc-1)) ; rm_disp(0) = 1
    ALLOCATE (all_total_send_offset(2, 0:numproc-1))
    ALLOCATE (blk_ps(nrow_images, ncol_images)) ; blk_ps = 1
    ALLOCATE (blks(nrow_images, ncol_images)) ; blks = 1

    ! Allocate and init mats
    !
    ALLOCATE (ums%mats(nrow_images,ncol_images))
    DO row_img = 1, nrow_images
       DO col_img = 1, ncol_images
          CALL dbcsr_init (ums%mats(row_img, col_img))
          CALL dbcsr_create(ums%mats(row_img, col_img), "imaged "//sm%name,&
               target_dist,&
               dbcsr_type_no_symmetry, array_data(sm%row_blk_size), array_data(sm%col_blk_size),&
               0,0, data_type,&
               data_memory_type=data_memory_type,&
               index_memory_type=index_memory_type,&
               error=dbcsr_error)
          ums%mats(row_img, col_img)%m%negate_real = sm%negate_real
          ums%mats(row_img, col_img)%m%negate_imaginary = sm%negate_imaginary
       ENDDO
    ENDDO
    
    ithread = 0
    nthreads = 1
!$omp parallel default (none) &
!$omp private (ithread, symmetry_i, row_img, col_img, &
!$omp          myt_send_count, myt_total_send_count, &
!$omp          iter, row, col, blk, row_size, col_size, stored_row, stored_col, &
!$omp          prow, pcol, rowi, coli, vrow, vcol, dst_p, nze, myt_smp, myt_sdp, &
!$omp          blk_p, bp, sd_pos, sm_pos, t_error, tr, &
!$omp          tr_row_size, tr_col_size) &
!$omp shared (nthreads, nsymmetries, row_img_dist, col_img_dist, &
!$omp         nrow_images, ncol_images, numproc, scale_value, &
!$omp         ums, sm, ism, target_imgdist, row_dist, col_dist, error, &
!$omp         predist_type_fwd, blacs2mpi, row_blk_size, col_blk_size, &
!$omp         send_count, recv_count, error_handler2, dbcsr_error, mp_group, &
!$omp         total_send_count, total_recv_count, recv_data_area, nocopy, &
!$omp         data_type, data_memory_type, recv_meta, send_data_area, send_meta, &
!$omp         sd_disp, sm_disp, rd_disp, rm_disp, all_total_send_offset, blk_ps, blks, &
!$omp         received_data_area, scale_neg_one)
!$  ithread = omp_get_thread_num()
!$omp master
!$  nthreads = omp_get_num_threads()
!$omp end master

    ! Allocate thread private data
    !
    ALLOCATE (myt_send_count(2, nrow_images, ncol_images, 0:numproc-1)) ; myt_send_count(:,:,:,:) = 0
    ALLOCATE (myt_total_send_count(2, 0:numproc-1))
    ! Thread-local pointers of the current adding position into the send buffers
    ALLOCATE (myt_smp(0:numproc-1), myt_sdp(0:numproc-1))

    ! Count sizes for sending
    !
    CALL dbcsr_iterator_start(iter, ism, shared=.TRUE.)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk,&
            row_size=row_size, col_size=col_size)
       nze = row_size * col_size
       IF (nze.EQ.0) CYCLE
       DO symmetry_i = 1, nsymmetries
          IF (symmetry_i .EQ. 1) THEN
             stored_row = row ; stored_col = col
          ELSE
             IF (row .EQ. col) CYCLE
             stored_row = col ; stored_col = row
          ENDIF
          ! Where do we send this block?
          row_img = 1
          IF (nrow_images .GT. 1) row_img = row_img_dist (stored_row)
          col_img = 1
          IF (ncol_images .GT. 1) col_img = col_img_dist (stored_col)
          CALL image_calculator(target_imgdist,&
               prow = prow, rowi = rowi,&
               pcol = pcol, coli = coli,&
               vprow = vrow, vpcol = vcol,&
               myprow = row_dist(stored_row), myrowi = row_img,&
               mypcol = col_dist(stored_col), mycoli = col_img,&
               shifting = predist_type_fwd, error=error)
          dst_p = blacs2mpi(prow, pcol)
          ! These counts are meant for the thread that processes this row.
          myt_send_count(1, rowi, coli, dst_p) =&
               myt_send_count(1, rowi, coli, dst_p) + 1
          myt_send_count(2, rowi, coli, dst_p) =&
               myt_send_count(2, rowi, coli, dst_p) + nze
       END DO  ! symmetry_i
    END DO
    CALL dbcsr_iterator_stop(iter)
    DO dst_p = 0, numproc-1
       myt_total_send_count(1, dst_p) = SUM(myt_send_count(1,:,:,dst_p))
       myt_total_send_count(2, dst_p) = SUM(myt_send_count(2,:,:,dst_p))
    END DO
    ! Merge the send counts
!$omp critical
    send_count(:,:,:,:) = send_count(:,:,:,:) + myt_send_count(:,:,:,:)
    total_send_count(:,:) = total_send_count(:,:) + myt_total_send_count(:,:)
!$omp end critical
!$omp barrier
!$omp master
    CALL dbcsr_error_set(routineN//"_sizes", error_handler2, dbcsr_error)
    CALL mp_alltoall(send_count, recv_count, 2*nrow_images*ncol_images,&
         mp_group)
    CALL dbcsr_error_stop(error_handler2, dbcsr_error)
!$omp end master
!$omp barrier
    ! Fill in the meta data structures and copy the data.
!$omp do
    DO dst_p = 0, numproc-1
       total_recv_count(1, dst_p) = SUM (recv_count (1, :, :, dst_p))
       total_recv_count(2, dst_p) = SUM (recv_count (2, :, :, dst_p))
    ENDDO
!$omp master
    ! Allocate data structures needed for data exchange.
    CALL dbcsr_data_init (recv_data_area)
    IF (nrow_images.EQ.1 .AND. ncol_images.eq.1 .OR. nocopy) THEN
       ! For some cases the faster dbcsr_special_finalize(reshuffle=.FALSE.) can be used.
       ! This basically makes this working matrix the actual data-area.
       ! Hence, for those cases we have to use data_memory_type already here.
       CALL dbcsr_data_new (recv_data_area, data_type, SUM(total_recv_count(2, :)), memory_type=data_memory_type)
    ELSE
       CALL dbcsr_data_new (recv_data_area, data_type, SUM(total_recv_count(2, :)))
    END IF
    ALLOCATE (recv_meta(metalen*SUM(total_recv_count(1, :))))
    CALL dbcsr_data_init (send_data_area)
    CALL dbcsr_data_new (send_data_area, data_type, SUM(total_send_count(2, :)))
    ALLOCATE (send_meta(metalen*SUM(total_send_count(1, :))))
    ! Calculate displacements for processors needed for the exchanges.
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
    myt_smp(:) = sm_disp(:)
    myt_sdp(:) = sd_disp(:)
    IF (nthreads .GT. 1) THEN
       all_total_send_offset(1,:) = myt_smp(:) + metalen*myt_total_send_count(1,:)
       all_total_send_offset(2,:) = myt_sdp(:) + myt_total_send_count(2,:)
    ENDIF
!$omp end master
!$omp barrier
    IF (ithread .GT. 0) THEN
!$omp critical
       myt_smp(:) = all_total_send_offset(1,:)
       myt_sdp(:) = all_total_send_offset(2,:)
       all_total_send_offset(1,:) &
            = all_total_send_offset(1,:) + metalen*myt_total_send_count(1,:)
       all_total_send_offset(2,:) &
            = all_total_send_offset(2,:) + myt_total_send_count(2,:)
!$omp end critical
    ELSE
       CALL dbcsr_data_init (received_data_area)
       received_data_area = recv_data_area
       CALL dbcsr_data_hold(received_data_area)       
       t_error = dbcsr_error
       DO row_img = 1, nrow_images
          DO col_img = 1, ncol_images
             CALL dbcsr_work_create(ums%mats(row_img, col_img),&
                  SUM(recv_count(1,row_img,col_img,:)), n=1, error=t_error)
             CALL dbcsr_data_hold (received_data_area)
             CALL dbcsr_data_release (ums%mats(row_img,col_img)%m%wms(1)%data_area)
             ums%mats(row_img,col_img)%m%wms(1)%data_area = received_data_area
          ENDDO
       ENDDO
    ENDIF
!$omp barrier

    ! Add timing call to the packing of the send buffers 
    !
    CALL dbcsr_error_set(routineN//"_pack", error_handler2, dbcsr_error)
    ! Copies metadata and actual data to be sent into the send buffers.
    CALL dbcsr_iterator_start(iter, ism, shared=.TRUE.)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk, blk_p=blk_p,&
            row_size=row_size, col_size=col_size)
       nze = row_size * col_size
       IF (nze.EQ.0) CYCLE
       bp = ABS(blk_p)
       DO symmetry_i = 1, nsymmetries
          IF (symmetry_i .EQ. 1) THEN
             stored_row = row ; stored_col = col; tr = blk_p.LT.0
             tr_row_size = col_size; tr_col_size = row_size
          ELSE
             IF (row .EQ. col) CYCLE
             stored_row = col ; stored_col = row; tr = blk_p.GT.0
             tr_row_size = row_size; tr_col_size = col_size
          ENDIF
          ! Where do we send this block?
          row_img = 1
          IF (nrow_images .GT. 1) row_img = row_img_dist (stored_row)
          col_img = 1
          IF (ncol_images .GT. 1) col_img = col_img_dist (stored_col)
          CALL image_calculator(target_imgdist,&
               prow = prow, rowi = rowi,&
               pcol = pcol, coli = coli,&
               vprow = vrow, vpcol = vcol,&
               myprow = row_dist(stored_row), myrowi = row_img,&
               mypcol = col_dist(stored_col), mycoli = col_img,&
               shifting = predist_type_fwd, error=error)
          dst_p = blacs2mpi(prow, pcol)
          sm_pos = myt_smp(dst_p)
          myt_smp(dst_p) = myt_smp(dst_p) + metalen
          sd_pos = myt_sdp(dst_p)
          myt_sdp(dst_p) = myt_sdp(dst_p) + nze
          IF (tr) THEN
             CALL dbcsr_block_transpose_aa(send_data_area,sm%data_area,tr_row_size,tr_col_size,&
                                           sd_pos,bp,scale_value,error=dbcsr_error)
             IF (sm%negate_real .AND. sm%negate_imaginary) THEN
                CALL dbcsr_block_scale (send_data_area, scale=scale_neg_one,&
                                        row_size=nze, col_size=1, lb=sd_pos, error=error)
             ELSEIF (sm%negate_real) THEN
                CALL dbcsr_block_real_neg (send_data_area,&
                                           row_size=nze, col_size=1, lb=sd_pos, error=error)
             ELSEIF (sm%negate_imaginary) THEN
                CALL dbcsr_block_conjg (send_data_area,&
                                        row_size=nze, col_size=1, lb=sd_pos, error=error)
             ENDIF
          ELSE
             CALL dbcsr_block_copy_aa(send_data_area,sm%data_area,row_size,col_size,&
                                      sd_pos,bp,scale_value,error=dbcsr_error)
          END IF

          send_meta(sm_pos) = stored_row
          send_meta(sm_pos+1) = stored_col
          send_meta(sm_pos+2) = sd_pos-sd_disp(dst_p)+1
          send_meta(sm_pos+3) = rowi
          send_meta(sm_pos+4) = coli
       ENDDO ! symmetry_i
    ENDDO ! iterator
    CALL dbcsr_iterator_stop(iter)

    ! Deallocate thread private data
    !
    DEALLOCATE (myt_send_count)
    DEALLOCATE (myt_total_send_count)
    DEALLOCATE (myt_smp, myt_sdp)

    CALL dbcsr_error_stop(error_handler2, dbcsr_error)
!$omp end parallel

    ! Exchange the data and metadata structures. In the interesting cases (square grids, row col distribution same),
    ! there are only very few processors that need to exchange data.
    ! The hybrid_alltoall deals with this by doing point to point communication
    CALL dbcsr_error_set(routineN//"_data", error_handler2, dbcsr_error)
    CALL hybrid_alltoall_any(send_data_area,total_send_count(2,:),sd_disp(:)-1,&
                             recv_data_area,total_recv_count(2,:),rd_disp(:)-1,&
                             mp_obj,&
                             most_ptp=.TRUE.,remainder_ptp=.TRUE.,no_hybrid=.FALSE.,&
                             error=dbcsr_error)
    CALL hybrid_alltoall_i1(&
         send_meta(:), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), metalen*total_recv_count(1,:), rm_disp(:)-1,&
         most_ptp=.TRUE., remainder_ptp=.TRUE., no_hybrid=.FALSE.,&
         mp_env = mp_obj)
    CALL dbcsr_error_stop(error_handler2, dbcsr_error)

    ! Now create the work index and/or copy the relevant data from the
    ! receive buffer into the local indices.
    !
    prev_blk_p = 0
    DO src_p = 0, numproc-1
       data_p = 0
       DO blk_l = 1, total_recv_count(1, src_p)
          stored_row = recv_meta(rm_disp(src_p)+metalen*(blk_l-1))
          stored_col = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+1)
          stored_blk_p = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+2)
          row_img = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+3)
          col_img = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+4)
          nze = row_blk_size(ABS(stored_row))&
               * col_blk_size(ABS(stored_col))
          blk = blks(row_img,col_img)
          blks(row_img,col_img) = blks(row_img,col_img) + 1
          IF (stored_col .GT. 0) THEN
             blk_p = data_p
             data_p = data_p + nze
          ELSE
             blk_p = prev_blk_p
          ENDIF
          blk_ps(row_img,col_img) = blk_ps(row_img,col_img) + nze
          ums%mats(row_img,col_img)%m%wms(1)%row_i(blk) = ABS(stored_row)
          ums%mats(row_img,col_img)%m%wms(1)%col_i(blk) = ABS(stored_col)
          ums%mats(row_img,col_img)%m%wms(1)%blk_p(blk) =&
               SIGN(rd_disp(src_p) + ABS(stored_blk_p)-1, stored_blk_p)
          prev_blk_p = blk_p
       ENDDO
    ENDDO

    ! Finalize the actual imaged matrices from the work matrices
    !
    DO row_img = 1, nrow_images
       DO col_img = 1, ncol_images
          ums%mats(row_img,col_img)%m%wms(1)%lastblk = blks(row_img,col_img) - 1
          ums%mats(row_img,col_img)%m%wms(1)%datasize = blk_ps(row_img,col_img) - 1
          CALL dbcsr_data_set_size_referenced (&
               ums%mats(row_img,col_img)%m%wms(1)%data_area,&
               ums%mats(row_img,col_img)%m%wms(1)%datasize)

          IF (nrow_images.EQ.1 .AND. ncol_images.eq.1 .OR. nocopy) THEN
             CALL dbcsr_special_finalize(ums%mats(row_img,col_img),reshuffle=.FALSE.,error=dbcsr_error)
          ELSE
             CALL dbcsr_special_finalize(ums%mats(row_img,col_img),reshuffle=.TRUE.,error=dbcsr_error)
          ENDIF

          ! Save the home process and image row and column
          CALL image_calculator (target_imgdist,&
               ums%mats(row_img,col_img)%m%index(dbcsr_slot_home_prow),&
               ums%mats(row_img,col_img)%m%index(dbcsr_slot_home_rowi),&
               ums%mats(row_img,col_img)%m%index(dbcsr_slot_home_pcol),&
               ums%mats(row_img,col_img)%m%index(dbcsr_slot_home_coli),&
               vprow = ums%mats(row_img,col_img)%m%index(dbcsr_slot_home_vprow),&
               vpcol = ums%mats(row_img,col_img)%m%index(dbcsr_slot_home_vpcol),&
               myrowi=row_img, mycoli=col_img,&
               shifting=predist_type, error=error)
       ENDDO
    ENDDO

    ! Deallocate shared temporary buffers
    !
    DEALLOCATE(send_count,recv_count)
    DEALLOCATE(total_send_count,total_recv_count)
    DEALLOCATE(sdp,smp,sd_disp,sm_disp)
    DEALLOCATE(rd_disp,rm_disp)
    DEALLOCATE(all_total_send_offset)
    DEALLOCATE(blk_ps,blks)
    DEALLOCATE(recv_meta,send_meta)
    
    CALL dbcsr_data_release (send_data_area)
    CALL dbcsr_data_release (received_data_area)
    CALL dbcsr_data_release (recv_data_area)

    CALL dbcsr_error_stop(error_handler, dbcsr_error)
  END SUBROUTINE make_images

! *****************************************************************************
!> \brief Makes dense matrices for the image matrices.
!> \param[in,out] images          current (undense) matrix images, output is
!>                                the dense matrix images
!> \param[in] new_rdist           the new image distribution for dense matrices
!> \param[in] row_map             mapping of current (undense) rows to dense rows
!> \param[in] col_map             mapping of current (undense) columns to
!>                                dense columns
!> \param[in] join_cols           (optional) make columns dense, default is
!>                                yes
!> \param[in] join_rows           (optional) make rows dense, default is yes
!> \param[in] new_template        template dense matrix for creating image
!>                                matrices
!> \param[in,out] error           error
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE dbcsr_make_images_dense (images, new_rdist, &
       row_map, col_map, join_cols, join_rows, new_template, error)
    TYPE(dbcsr_2d_array_type), INTENT(INOUT) :: images
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(INOUT)                          :: new_rdist
    TYPE(array_i1d_obj), INTENT(IN)          :: row_map, col_map
    LOGICAL, INTENT(IN)                      :: join_cols, join_rows
    TYPE(dbcsr_obj), INTENT(IN)              :: new_template
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_images_dense', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handler, mat_col, &
                                                mat_row, mat_vpcol, mat_vprow
    INTEGER, DIMENSION(:), POINTER           :: und_col_blk_offsets, &
                                                und_row_blk_offsets
    INTEGER, DIMENSION(dbcsr_meta_size)      :: old_meta
    REAL(kind=dp)                            :: cs
    TYPE(array_i1d_obj)                      :: dense_local_vcols, &
                                                dense_local_vrows, &
                                                und_local_vcols, &
                                                und_local_vrows
    TYPE(dbcsr_imagedistribution_obj)        :: old_rdist
    TYPE(dbcsr_obj)                          :: tmp_mat

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    old_rdist = images%image_dist
    !
    DO mat_row = 1, images%image_dist%i%row_decimation
       DO mat_col = 1, images%image_dist%i%col_decimation
          IF (dbg) THEN
             cs = dbcsr_checksum (images%mats(mat_row,mat_col), error=error)
             WRITE(*,*)routineN//" cs pre", cs
          ENDIF
          mat_vprow = images%mats(mat_row, mat_col)%m%index(dbcsr_slot_home_vprow)
          mat_vpcol = images%mats(mat_row, mat_col)%m%index(dbcsr_slot_home_vpcol)
          und_row_blk_offsets => array_data (images%mats(mat_row, mat_col)%m%row_blk_offset)
          und_col_blk_offsets => array_data (images%mats(mat_row, mat_col)%m%col_blk_offset)
          CALL dbcsr_get_local_vrows (old_rdist, und_local_vrows, mat_vprow, error=error)
          CALL dbcsr_get_local_vcols (old_rdist, und_local_vcols, mat_vpcol, error=error)

          CALL dbcsr_get_local_vrows (new_rdist, dense_local_vrows, mat_vprow, error=error)
          CALL dbcsr_get_local_vcols (new_rdist, dense_local_vcols, mat_vpcol, error=error)
          ! The old matrix has to be remembered so it is copied to
          ! tmp_mat.
          old_meta(:) = images%mats(mat_row,mat_col)%m%index(1:dbcsr_meta_size)
          CALL dbcsr_init (tmp_mat)
          tmp_mat = images%mats(mat_row,mat_col)
          CALL dbcsr_init (images%mats(mat_row,mat_col))
          CALL dbcsr_create (images%mats(mat_row,mat_col), template=new_template,&
               error=error)
          images%mats(mat_row,mat_col)%m%index(dbcsr_slot_home_prow&
                                              :dbcsr_slot_home_vpcol) =&
               old_meta(dbcsr_slot_home_prow:dbcsr_slot_home_vpcol)
          CALL dbcsr_make_dense_low(tmp_mat, images%mats(mat_row,mat_col),&
               array_data (und_local_vrows), array_data (und_local_vcols),&
               und_row_blk_offsets, und_col_blk_offsets,&
               array_data (dense_local_vrows),&
               array_data (dense_local_vcols),&
               array_data (new_template%m%row_blk_offset),&
               array_data (new_template%m%col_blk_offset),&
               array_data(row_map), array_data(col_map), join_rows, join_cols,&
               error)
          !
          CALL dbcsr_index_prune_deleted (images%mats(mat_row, mat_col), error=error)
          !
          CALL dbcsr_release (tmp_mat)
          IF (dbg) THEN
             cs = dbcsr_checksum (images%mats(mat_row,mat_col), error=error)
             WRITE(*,*)routineN//" cs pst", cs
          ENDIF
       ENDDO
    ENDDO
    CALL dbcsr_image_dist_release (images%image_dist, error=error)
    images%image_dist = new_rdist
    CALL dbcsr_image_dist_hold (images%image_dist, error=error)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_make_images_dense

! *****************************************************************************
!> \brief Multiplies two DBCSR matrices
!>
!> \param[in] left_set             set of imaged left matrices
!> \param[in] right_set            set of imaged right matrices
!> \param[out] product_matrix      DBCSR product matrix 
!> \param[in,out] error            cp2k error
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix; default is no
!> \param filter_eps ...
!> \param[out] flop                (optional) effective flop
! *****************************************************************************
  SUBROUTINE cannon_multiply_low(left_set, right_set, product_matrix,&
       error, retain_sparsity, &
       filter_eps, flop)
    TYPE(dbcsr_2d_array_type), POINTER       :: left_set, right_set
    TYPE(dbcsr_obj), INTENT(INOUT)           :: product_matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(kind=real_8), INTENT(in), OPTIONAL  :: filter_eps
    INTEGER(KIND=int_8), INTENT(OUT)         :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'cannon_multiply_low', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=80), PARAMETER :: &
      fdata = '(A,1X,I4,"(",2(I3),"x",2(I3),")","(",I3,"x",I3,")")', fxfer = &
      '(A,1X,I4,"->",I4,2(1X,"(",I3,"x",I3,")"),1X,"IM (",I3,"x",I3,")")'
    INTEGER, PARAMETER :: id_bytes = 3, id_recv = 2, id_send = 1, &
      id_time = 1, id_waittime = 2, idata = 1, ileft = 0, imeta = 2, &
      iright = 2, M_L = 2, M_P = 1, M_R = 3, RC_C = 2, RC_R = 1

    INTEGER :: data_type, data_type_byte, error_handler, error_handler2, &
      error_handler3, error_handler4, grp, i, ithread, left_col_image, &
      left_col_mult, left_col_nimages, left_dst_icol, left_dst_irow, &
      left_dst_p, left_dst_pcol, left_dst_prow, left_dst_vcol, left_dst_vrow, &
      left_myfirstvcol, left_myfirstvrow, left_mypcol, left_myprow, &
      left_npcols, left_nprows, left_recv_icol, left_recv_irow, left_recv_p, &
      left_recv_pcol, left_recv_prow, left_recv_vcol, left_recv_vrow, &
      left_row_image, left_row_mult, left_row_nimages, left_send_icol, &
      left_send_irow, left_send_p, left_send_pcol, left_send_prow
    INTEGER :: left_send_vcol, left_send_vrow, left_src_icol, left_src_irow, &
      left_src_p, left_src_pcol, left_src_prow, left_src_vcol, left_src_vrow, &
      metronome, min_nimages, mp_group, mynode, nblkrows_total, &
      nblkrows_used, nsteps_k, nthreads, numnodes, nvirt_k, output_unit, &
      right_col_image, right_col_mult, right_col_nimages, right_dst_icol, &
      right_dst_irow, right_dst_p, right_dst_pcol, right_dst_prow, &
      right_dst_vcol, right_dst_vrow, right_myfirstvcol, right_myfirstvrow, &
      right_mypcol, right_myprow, right_npcols, right_nprows, &
      right_recv_icol, right_recv_irow, right_recv_p, right_recv_pcol
    INTEGER :: right_recv_prow, right_recv_vcol, right_recv_vrow, &
      right_row_image, right_row_mult, right_row_nimages, right_send_icol, &
      right_send_irow, right_send_p, right_send_pcol, right_send_prow, &
      right_send_vcol, right_send_vrow, right_src_icol, right_src_irow, &
      right_src_p, right_src_pcol, right_src_prow, right_src_vcol, &
      right_src_vrow, row, size_guess, stat, threads_finished, &
      threads_finished_read, v_k, v_ki
    INTEGER(KIND=int_8)                      :: flop_single, flop_total
    INTEGER(kind=mp_address_kind)            :: left_max_nblks, left_max_nze, &
                                                right_max_nblks, right_max_nze
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: row_counts, total_row_counts
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: left_sizes, my_sizes, &
                                                right_sizes
    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: all_sizes
    INTEGER, DIMENSION(3, 2)                 :: mp_rc_groups
    INTEGER, DIMENSION(:), POINTER :: left_data_rr, left_data_sr, &
      left_index_rp, left_index_rr, left_index_sp, left_index_sr, local_rows, &
      right_data_rr, right_data_sr, right_index_rp, right_index_rr, &
      right_index_sp, right_index_sr
    INTEGER, DIMENSION(:, :), POINTER        :: left_pgrid, product_pgrid, &
                                                right_pgrid
    INTEGER, SAVE                            :: mult_id = 0
    LOGICAL                                  :: keep_sparsity, list_indexing, &
                                                otf_filtering
    REAL(KIND=dp)                            :: checksum

!$  REAL(KIND=real_8)                        :: left_fill, right_fill
    REAL(kind=sp), ALLOCATABLE, DIMENSION(:) :: left_norms, right_norms, &
                                                row_max_epss
    TYPE(dbcsr_2d_array_type), POINTER :: left_buffer_2, left_buffer_calc, &
      left_buffer_comm, right_buffer_2, right_buffer_calc, right_buffer_comm
    TYPE(dbcsr_data_obj)                     :: left_data_rp, left_data_sp, &
                                                right_data_rp, right_data_sp
    TYPE(dbcsr_data_obj), POINTER            :: trs_stackbuf_calc, &
                                                trs_stackbuf_comm
    TYPE(dbcsr_data_obj), TARGET             :: trs_stackbuf_1, trs_stackbuf_2
    TYPE(dbcsr_error_type)                   :: t_error
    TYPE(dbcsr_mm_multrec_type_p), DIMENSION(:), ALLOCATABLE :: multrec
    TYPE(dbcsr_mp_obj)                       :: left_mp_obj, product_mp_obj, &
                                                right_mp_obj

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    NULLIFY(trs_stackbuf_calc, trs_stackbuf_comm)
    !
    ALLOCATE (left_buffer_2, right_buffer_2)
    mult_id=mult_id+1

    IF (PRESENT (retain_sparsity)) THEN
       keep_sparsity = retain_sparsity
    ELSE
       keep_sparsity = .FALSE.
    ENDIF
    otf_filtering = PRESENT (filter_eps)

!$omp parallel default (none) &
!$omp shared (multrec, nthreads, product_matrix, error)
!$omp master
    nthreads = 1
    !$  nthreads = OMP_GET_NUM_THREADS ()
    CALL dbcsr_assert (ASSOCIATED (product_matrix%m%wms),&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Work matrices do not exist",__LINE__,error)
    CALL dbcsr_assert (SIZE (product_matrix%m%wms), "EQ", nthreads,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Work matrices not correctly sized.",__LINE__,error)
    ALLOCATE(multrec(0:nthreads-1))
!$omp end master
!$omp end parallel

    output_unit = default_output_unit
    flop_total = 0
    ! Set up variables
    data_type = dbcsr_get_data_type (product_matrix)
    data_type_byte = dbcsr_datatype_sizeof(data_type)
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
    
    dbcsr_mpi_statistics%col_nimages = left_col_nimages
    dbcsr_mpi_statistics%row_nimages = right_row_nimages
    dbcsr_mpi_statistics%nprocs = numnodes

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
       DO row = 1 , nblkrows_total
          row_max_epss (row) &
               = REAL(filter_eps&
               / REAL(MAX(1, total_row_counts(row)), KIND=KIND(row_max_epss)),&
               KIND=KIND(row_max_epss))
       END DO
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
!$omp          private (i, size_guess, &
!$omp                   left_fill, right_fill, ithread) &
!$omp          shared (product_matrix, left_max_nze, right_max_nze) &
!$omp          shared (left_set, right_set, &
!$omp                 left_col_nimages, right_row_nimages) &
!$omp          shared (error, nthreads, keep_sparsity, mynode)
    !
    !$ ithread = OMP_GET_THREAD_NUM()
    ! The work arrays have to be setup (actually, not quite sure).
    i = ithread + 1
    size_guess = product_matrix%m%wms(i)%datasize ! Should be minimal
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

    ! update capacity of memory-pools
    IF (mm_driver == mm_driver_acc) THEN
      CALL dbcsr_mempool_ensure_capacity(memtype_abpanel_1%pool, &
      capacity=left_row_mult*left_col_nimages + right_row_nimages*right_col_mult)
      CALL dbcsr_mempool_ensure_capacity(memtype_abpanel_2%pool, &
      capacity=left_row_mult*left_col_nimages + right_row_nimages*right_col_mult)
    ENDIF

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

    CALL buffer_matrices_ensure_size(left_set, index_size=left_max_nblks,&
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
    left_buffer_calc => left_set
    left_buffer_comm => left_buffer_2
    ALLOCATE (left_data_sr  (left_col_nimages))
    ALLOCATE (left_index_sr (left_col_nimages))
    ALLOCATE (left_data_rr  (left_col_nimages))
    ALLOCATE (left_index_rr (left_col_nimages))
    left_data_sr = mp_request_null
    left_data_rr = mp_request_null
    left_index_sr = mp_request_null
    left_index_rr = mp_request_null

    ! Setup buffers for right matrix
    CALL buffer_matrices_ensure_size(right_set, index_size=right_max_nblks, &
         data_size=right_max_nze, error=error)
    
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
    right_buffer_calc => right_set
    right_buffer_comm => right_buffer_2
    ALLOCATE (right_data_sr  (right_row_nimages))
    ALLOCATE (right_index_sr (right_row_nimages))
    ALLOCATE (right_data_rr  (right_row_nimages))
    ALLOCATE (right_index_rr (right_row_nimages))
    right_data_sr = mp_request_null
    right_data_rr = mp_request_null
    right_index_sr = mp_request_null
    right_index_rr = mp_request_null
    !
!$omp parallel &
!$omp default (none) &
!$omp shared (left_buffer_comm, right_buffer_comm, product_matrix,&
!$omp         keep_sparsity, filter_eps, row_max_epss, multrec, &
!$omp         right_data_sr, right_data_rr, left_data_sr, left_data_rr,&
!$omp         right_index_sr, right_index_rr, left_index_sr, left_index_rr), &
!$omp private(ithread),&
!$omp firstprivate (error)
    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM ()
    ALLOCATE(multrec(ithread)%p)
    CALL dbcsr_mm_multrec_init(multrec(ithread)%p,&
         left_buffer_comm%mats(1, 1)%m,&
         right_buffer_comm%mats(1, 1)%m,&
         product_matrix%m,&
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

    ! Setup transpose stackbuffers
    IF (mm_driver==mm_driver_acc) THEN
       CALL dbcsr_data_init(trs_stackbuf_1)
       CALL dbcsr_data_init(trs_stackbuf_2)
       CALL dbcsr_data_new(trs_stackbuf_1, data_type=dbcsr_type_int_4, data_size=1000, memory_type=memtype_trsbuffer_1)
       CALL dbcsr_data_new(trs_stackbuf_2, data_type=dbcsr_type_int_4, data_size=1000, memory_type=memtype_trsbuffer_2)
       trs_stackbuf_calc => trs_stackbuf_1
       trs_stackbuf_comm => trs_stackbuf_2
    ENDIF

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
       CALL dbcsr_error_set(routineN//"_metrocomm1", error_handler2, error)
       wait_right: IF (metronome .GT. 1) THEN
          IF (debug_mod) WRITE (*,'(1X,A)')routineN//" waiting for right"
          !
          CALL mp_waitall (right_data_sr)
          CALL mp_waitall (right_data_rr)
          CALL mp_waitall (right_index_sr)
          CALL mp_waitall (right_index_rr)
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
             IF (mm_driver == mm_driver_acc) THEN
                CALL dbcsr_error_set(routineN//"_acc_sync_right", error_handler3, error)
                CALL acc_event_synchronize(right_buffer_comm%mats(v_ki+1, 1)%m%data_area%d%acc_ready)
                CALL dbcsr_error_stop(error_handler3, error)
             ENDIF

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
             CALL dbcsr_error_set(routineN//"_metrocomm2", error_handler2, error)
             CALL dbcsr_irecv_any (right_data_rp, right_recv_p,&
                  grp, right_data_rr(v_ki+1), tag=right_src_vrow,&
                  error=error)
             CALL mp_irecv (right_index_rp, right_recv_p,&
                  grp, right_index_rr(v_ki+1), tag=right_src_vrow)
             CALL dbcsr_isend_any (right_data_sp, right_send_p,&
                  grp, right_data_sr(v_ki+1), tag=right_dst_vrow,&
                  error=error)
             CALL mp_isend (right_index_sp, right_send_p,&
                  grp, right_index_sr(v_ki+1), tag=right_dst_vrow)
             dbcsr_mpi_statistics%nexchanged = dbcsr_mpi_statistics%nexchanged+1
             CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(1,:),&
                  INT(dbcsr_data_get_size(right_data_rp),KIND=mp_address_kind),&
                  dbcsr_mpi_statistics%data_size_breakdown(1,:),data_type_byte)
             CALL count_mpi_statistics(dbcsr_mpi_statistics%index_size(1,:),&
                  INT(SIZE(right_index_rp),KIND=mp_address_kind))
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
       CALL dbcsr_error_set(routineN//"_metrocomm3", error_handler2, error)
       wait_left: IF (metronome .GT. 1) THEN
          IF (debug_mod) WRITE (*,'(1X,A)')routineN//" waiting for left"
          CALL mp_waitall (left_data_sr)
          CALL mp_waitall (left_data_rr)
          CALL mp_waitall (left_index_sr)
          CALL mp_waitall (left_index_rr)
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
             IF (mm_driver == mm_driver_acc) THEN
                CALL dbcsr_error_set(routineN//"_acc_sync_left", error_handler3, error)
                CALL acc_event_synchronize(left_buffer_comm%mats(1, v_ki+1)%m%data_area%d%acc_ready)
                CALL dbcsr_error_stop(error_handler3, error)
             ENDIF

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
             CALL dbcsr_error_set(routineN//"_metrocomm4", error_handler2, error)
             CALL dbcsr_irecv_any (left_data_rp, left_recv_p,&
                  grp, left_data_rr(v_ki+1), tag=left_src_vcol,&
                  error=error)
             CALL mp_irecv (left_index_rp, left_recv_p,&
                  grp, left_index_rr(v_ki+1), tag=left_src_vcol)
             CALL dbcsr_isend_any (left_data_sp, left_send_p,&
                  grp, left_data_sr(v_ki+1), tag=left_dst_vcol,&
                  error=error)
             CALL mp_isend (left_index_sp, left_send_p,&
                  grp, left_index_sr(v_ki+1), tag=left_dst_vcol)
             dbcsr_mpi_statistics%nexchanged = dbcsr_mpi_statistics%nexchanged+1
             CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(2,:),&
                  INT(dbcsr_data_get_size(left_data_rp),KIND=mp_address_kind),&
                  dbcsr_mpi_statistics%data_size_breakdown(2,:),data_type_byte)
             CALL count_mpi_statistics(dbcsr_mpi_statistics%index_size(2,:),&
                  INT(SIZE(left_index_rp),KIND=mp_address_kind))
             CALL dbcsr_error_stop(error_handler2, error)
          ENDDO
       ENDIF xfer_left
       !
       ! Repoint indices of left matrices and do the multiplications.
       calc_case_left:  IF (metronome .GT. 1) THEN
          DO v_ki = 0, left_col_nimages-1
             CALL dbcsr_repoint_index (left_buffer_calc%mats(1,v_ki+1)%m)
             left_buffer_calc%mats(1, v_ki+1)%m%valid=.TRUE.
          ENDDO
       ENDIF calc_case_left
        
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

          IF (mm_driver==mm_driver_acc) THEN
             CALL dbcsr_data_host2dev(left_buffer_calc%mats(1, v_ki+1)%m%data_area, error)
             CALL dbcsr_data_host2dev(right_buffer_calc%mats(v_ki+1, 1)%m%data_area, error)
             CALL acc_transpose_blocks(right_buffer_calc%mats(v_ki+1, 1), trs_stackbuf_calc, error)
          END IF
          
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

!$omp parallel default (none) &
!$omp shared (left_buffer_calc, right_buffer_calc, &
!$omp         v_ki, &
!$omp         product_matrix, multrec,&
!$omp         filter_eps, right_norms, left_norms, row_max_epss, &
!$omp         keep_sparsity, error, threads_finished, &
!$omp         right_data_sr, right_data_rr, right_index_sr, right_index_rr, &
!$omp         left_data_sr, left_data_rr, left_index_sr, left_index_rr, &
!$omp         use_comm_thread,error_handler2, error_handler4) &
!$omp private (ithread,nthreads, t_error, threads_finished_read) &
!$omp firstprivate (metronome, nsteps_k, min_nimages) &
!$omp reduction (+: flop_single)
          ithread = 0; nthreads = 1
!$        ithread = omp_get_thread_num(); nthreads = omp_get_num_threads()

          t_error=error
          CALL dbcsr_error_set(routineN//"_multrec", error_handler2, t_error)
          IF(ithread==0)&
               CALL dbcsr_error_set(routineN//"_multrec_master", error_handler4, t_error)
          
          IF(metronome==nsteps_k .AND. v_ki==min_nimages-1) &
               CALL dbcsr_mm_multrec_phaseout(multrec(ithread)%p, t_error)

          CALL dbcsr_mm_multrec_multiply(multrec(ithread)%p,&
               left=left_buffer_calc%mats(1, v_ki+1)%m,&
               right=right_buffer_calc%mats(v_ki+1, 1)%m,&
               flop=flop_single,&
               a_norms=left_norms, b_norms=right_norms,&
               error=t_error)

          IF(ithread==0) CALL dbcsr_error_stop(error_handler4, t_error)
          
          IF(metronome==nsteps_k .AND. v_ki==min_nimages-1) THEN
             CALL dbcsr_mm_multrec_finalize(multrec(ithread)%p, t_error)
             DEALLOCATE(multrec(ithread)%p)
          ENDIF

!$omp atomic
          threads_finished = threads_finished + 1
          IF (use_comm_thread .AND. (ithread .EQ. 0)) THEN
             DO 
! requires OMP 3.1 (e.g. gcc >=4.7), for correctness, otherwise we keep fingers crossed
#if defined _OPENMP && _OPENMP >= 200711
                !$OMP ATOMIC READ
#endif
                threads_finished_read=threads_finished
                IF (threads_finished_read .EQ. nthreads) EXIT
                CALL mp_testany(right_data_sr)
                CALL mp_testany(right_data_rr)
                CALL mp_testany(left_data_sr)
                CALL mp_testany(left_data_rr)
                CALL mp_testany(right_index_sr)
                CALL mp_testany(right_index_rr)
                CALL mp_testany(left_index_sr)
                CALL mp_testany(left_index_rr)
             END DO
          END IF
          !$OMP BARRIER
          CALL dbcsr_error_stop(error_handler2, t_error)
!$omp end parallel
          flop_total = flop_total + flop_single
       ENDDO

       CALL dbcsr_switch_sets (left_buffer_calc, left_buffer_comm)
       CALL dbcsr_switch_sets (right_buffer_calc, right_buffer_comm)
       CALL dbcsr_switch_d_ptrs(trs_stackbuf_calc, trs_stackbuf_comm)

    ENDDO grouped_k_index

    IF (mm_driver==mm_driver_acc) THEN
       CALL dbcsr_data_release(trs_stackbuf_1)
       CALL dbcsr_data_release(trs_stackbuf_2)
    END IF

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
    CALL dbcsr_destroy_array (right_buffer_2, error=error)
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
    DEALLOCATE (left_buffer_2, right_buffer_2)
    !
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE cannon_multiply_low


! *****************************************************************************
!> \brief Multiplies two DBCSR matrices by means of RMA MPI operations
!>
!> \param[in] left_set             set of imaged left matrices
!> \param[in] right_set            set of imaged right matrices
!> \param[out] product_matrix      DBCSR product matrix 
!> \param[in,out] error            cp2k error
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix; default is no
!> \param filter_eps ...
!> \param[out] flop                (optional) effective flop
! *****************************************************************************
  SUBROUTINE cannon_multiply_low_rma(left_set, right_set, product_matrix,&
       error, retain_sparsity, &
       filter_eps, flop)
    TYPE(dbcsr_2d_array_type), POINTER       :: left_set, right_set
    TYPE(dbcsr_obj), INTENT(INOUT)           :: product_matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(kind=real_8), INTENT(in), OPTIONAL  :: filter_eps
    INTEGER(KIND=int_8), INTENT(OUT)         :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'cannon_multiply_low_rma', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: idata_displ = 3, &
                                                idata_size = 1, &
                                                imeta_displ = 4, &
                                                imeta_size = 2

    INTEGER :: data_type, data_type_byte, dst_icol, dst_irow, dst_p, &
      dst_pcol, dst_prow, error_handler, error_handler2, error_handler3, &
      error_handler4, get_vcol, get_vrow, grp_left, grp_right, i, ithread, &
      left_col_mult, left_col_nimages, left_data_win, left_meta_win, &
      left_myfirstvcol, left_myfirstvrow, left_mypcol, left_myprow, &
      left_npcols, left_nprows, left_recv_p, left_row_mult, left_row_nimages, &
      left_send_p, metronome, min_nimages, mp_group, mynode, nblkrows_total, &
      nblkrows_used, nsteps_k, nthreads, numnodes, nvirt_k, output_unit, &
      recv_pcol, recv_prow, right_col_mult, right_col_nimages
    INTEGER :: right_data_win, right_meta_win, right_myfirstvcol, &
      right_myfirstvrow, right_mypcol, right_myprow, right_npcols, &
      right_nprows, right_recv_p, right_row_mult, right_row_nimages, &
      right_send_p, row, send_pcol, send_prow, size_guess, src_vcol, &
      src_vrow, stat, v_ki
    INTEGER(kind=mp_address_kind) :: left_data_size, left_max_nblks, &
      left_max_nze, left_meta_size, right_data_size, right_max_nblks, &
      right_max_nze, right_meta_size
    INTEGER(kind=mp_address_kind), &
      ALLOCATABLE, DIMENSION(:, :)           :: local_left_gets, &
                                                local_right_gets
    INTEGER(kind=mp_address_kind), &
      ALLOCATABLE, DIMENSION(:, :, :)        :: all_left_gets, &
                                                all_right_gets, left_gets, &
                                                right_gets
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: row_counts, total_row_counts
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: left_get_requests, left_gets_p, &
      precannon_requests, right_get_requests, right_gets_p
    INTEGER, DIMENSION(:), POINTER           :: left_meta_get, local_rows, &
                                                right_meta_get
    INTEGER, DIMENSION(:, :), POINTER        :: left_pgrid, product_pgrid, &
                                                right_pgrid
    LOGICAL                                  :: keep_sparsity, list_indexing, &
                                                otf_filtering
    REAL(KIND=dp)                            :: checksum
    REAL(kind=sp), ALLOCATABLE, DIMENSION(:) :: left_norms, right_norms, &
                                                row_max_epss
    TYPE(dbcsr_2d_array_type), POINTER :: left_buffer, left_buffer_calc, &
      left_buffer_comm, right_buffer, right_buffer_calc, right_buffer_comm
    TYPE(dbcsr_data_obj) :: left_data_get, left_data_orig, left_meta_orig, &
      right_data_get, right_data_orig, right_meta_orig
    TYPE(dbcsr_data_obj), POINTER            :: trs_stackbuf_calc, &
                                                trs_stackbuf_comm
    TYPE(dbcsr_data_obj), TARGET             :: trs_stackbuf_1, trs_stackbuf_2
    TYPE(dbcsr_error_type)                   :: t_error
    TYPE(dbcsr_mm_multrec_type_p), &
      ALLOCATABLE, DIMENSION(:)              :: multrec
    TYPE(dbcsr_mp_obj)                       :: left_mp_obj, product_mp_obj, &
                                                right_mp_obj

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    NULLIFY(trs_stackbuf_calc, trs_stackbuf_comm)
    !
    ALLOCATE (left_buffer, right_buffer)
    !
    IF (PRESENT (retain_sparsity)) THEN
       keep_sparsity = retain_sparsity
    ELSE
       keep_sparsity = .FALSE.
    ENDIF
    otf_filtering = PRESENT (filter_eps)

!$omp parallel default (none) &
!$omp shared (multrec, nthreads, product_matrix, error)
!$omp master
    nthreads = 1
    !$  nthreads = OMP_GET_NUM_THREADS ()
    CALL dbcsr_assert (ASSOCIATED (product_matrix%m%wms),&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Work matrices do not exist",__LINE__,error)
    CALL dbcsr_assert (SIZE (product_matrix%m%wms), "EQ", nthreads,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Work matrices not correctly sized.",__LINE__,error)
    ALLOCATE(multrec(0:nthreads-1))
!$omp end master
!$omp end parallel

    output_unit = default_output_unit
    flop = 0
    ! Set up variables
    data_type = dbcsr_get_data_type (product_matrix)
    data_type_byte = dbcsr_datatype_sizeof(data_type)
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
    left_myfirstvrow  = left_myprow*left_row_nimages
    left_myfirstvcol  = left_mypcol*left_col_nimages
    right_nprows      = dbcsr_mp_nprows(right_mp_obj)
    right_npcols      = dbcsr_mp_npcols(right_mp_obj)
    right_myprow      = dbcsr_mp_myprow(right_mp_obj)
    right_mypcol      = dbcsr_mp_mypcol(right_mp_obj)
    right_myfirstvrow = right_myprow*right_row_nimages
    right_myfirstvcol = right_mypcol*right_col_nimages
    mp_group = dbcsr_mp_group (product_mp_obj)
    left_pgrid => dbcsr_mp_pgrid (left_mp_obj)
    right_pgrid => dbcsr_mp_pgrid (right_mp_obj)
    product_pgrid => dbcsr_mp_pgrid (product_mp_obj)
    CALL dbcsr_mp_grid_setup (product_mp_obj)
    CALL dbcsr_mp_grid_setup (left_mp_obj)
    CALL dbcsr_mp_grid_setup (right_mp_obj)
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
    
    dbcsr_mpi_statistics%col_nimages = left_col_nimages
    dbcsr_mpi_statistics%row_nimages = right_row_nimages
    dbcsr_mpi_statistics%nprocs = numnodes
    !
    ! The main transfer loop goes through the virtual rows/columns.
    ! The number of steps may be smaller if the grid dimension is very
    ! non-optimal (both left column images and right row images are >
    ! 1).
    min_nimages = MIN (left_col_nimages, right_row_nimages)
    nvirt_k = left_npcols * left_col_nimages
    nsteps_k = nvirt_k / min_nimages
    !
    ! Take the dimensions of the local data and meta
    ALLOCATE(local_right_gets(4,right_row_nimages))
    right_data_size = 0 ; right_meta_size = 0
    DO v_ki = 1, right_row_nimages
       local_right_gets(idata_size,v_ki) = dbcsr_data_get_size_referenced(&
            right_set%mats(v_ki, 1)%m%data_area)
       local_right_gets(imeta_size,v_ki) = &
            right_set%mats(v_ki, 1)%m%index(dbcsr_slot_size)
       local_right_gets(idata_displ,v_ki) = right_data_size
       local_right_gets(imeta_displ,v_ki) = right_meta_size
       right_data_size = right_data_size+local_right_gets(idata_size,v_ki)
       right_meta_size = right_meta_size+local_right_gets(imeta_size,v_ki)       
    ENDDO
    !
    ALLOCATE(local_left_gets(4,left_col_nimages))
    left_data_size = 0 ; left_meta_size = 0
    DO v_ki = 1, left_col_nimages
       local_left_gets(idata_size,v_ki) = dbcsr_data_get_size_referenced(&
            left_set%mats(1, v_ki)%m%data_area)
       local_left_gets(imeta_size,v_ki) = &
            left_set%mats(1, v_ki)%m%index(dbcsr_slot_size)
       local_left_gets(idata_displ,v_ki) = left_data_size
       local_left_gets(imeta_displ,v_ki) = left_meta_size
       left_data_size = left_data_size+local_left_gets(idata_size,v_ki)
       left_meta_size = left_meta_size+local_left_gets(imeta_size,v_ki)       
    ENDDO
    !
    ! Exchange size data
    ALLOCATE (all_right_gets(4,right_row_nimages,0:numnodes-1))
    CALL mp_allgather(local_right_gets,all_right_gets,mp_group)
    ALLOCATE (all_left_gets(4,left_col_nimages,0:numnodes-1))
    CALL mp_allgather(local_left_gets,all_left_gets,mp_group)
    !
    ! Translate the all_sizes to account for pre-distribution.  This
    ! is just done to simplify lookups.
    ALLOCATE (right_gets(4, 0:nvirt_k-1, 0:right_npcols-1))
    right_gets = -1
    DO src_vcol = 0, right_npcols-1
       DO src_vrow = 0, right_row_nimages*right_nprows-1
          ! Calculate what was shifted.  The right_src_v{row,col} are
          ! the "source" rows/columns; the right_dst are the shifted
          ! targets where the data was placed in make_images.
          CALL image_calculator(right_set%image_dist,&
               prow = dst_prow, pcol = dst_pcol,&
               rowi = dst_irow, coli = dst_icol,&
               myvprow = src_vrow, myvpcol = src_vcol,&
               shifting = 'r', error=error)
          dst_p = right_pgrid (dst_prow, dst_pcol)
          right_gets(:, src_vrow, src_vcol) =&
               all_right_gets(:, dst_irow, dst_p)
       ENDDO
    ENDDO
    !
    ALLOCATE (left_gets(4, 0:left_nprows-1, 0:nvirt_k-1))
    left_gets = -1
    DO src_vcol = 0, left_col_nimages*left_npcols-1
       DO src_vrow = 0, left_nprows-1
          ! Calculate what was shifted.  The left_src_v{row,col} are
          ! the "source" rows/columns; the left_dst are the shifted
          ! targets where the data was placed in make_images.
          CALL image_calculator(left_set%image_dist,&
               prow = dst_prow, pcol = dst_pcol,&
               rowi = dst_irow, coli = dst_icol,&
               myvprow = src_vrow, myvpcol = src_vcol,&
               shifting = 'l', error=error)
          dst_p = left_pgrid (dst_prow, dst_pcol)
          left_gets(:, src_vrow, src_vcol) =&
               all_left_gets(:, dst_icol, dst_p)
       ENDDO
    ENDDO
    !
    right_max_nze = MAXVAL(all_right_gets(idata_size,:,:))
    right_max_nblks = MAXVAL(all_right_gets(imeta_size,:,:))
    left_max_nze = MAXVAL(all_left_gets(idata_size,:,:))
    left_max_nblks = MAXVAL(all_left_gets(imeta_size,:,:))
    !
    ! Setup product work areas
!$omp parallel default(none) &
!$omp          private (i, size_guess, ithread) &
!$omp          shared (product_matrix, error)
    !
    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM()
    ! The work arrays have to be setup (actually, not quite sure).
    i = ithread + 1
    size_guess = product_matrix%m%wms(i)%datasize ! Should be minimal
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
    ! update capacity of memory-pools
    IF (mm_driver == mm_driver_acc) THEN
       CALL dbcsr_mempool_ensure_capacity(memtype_abpanel_1%pool, &
            capacity=left_col_nimages + right_row_nimages)
       CALL dbcsr_mempool_ensure_capacity(memtype_abpanel_2%pool, &
            capacity=left_col_nimages + right_row_nimages)
    ENDIF
    !
    ! Prepare buffers for computation
    ! Left:
    CALL buffer_matrices_ensure_size(left_set, index_size=left_max_nblks,&
         data_size=left_max_nze, error=error)
    CALL setup_buffer_matrices (left_buffer, 1, left_col_nimages,&
         left_set%mats(1,1), index_size=left_max_nblks, data_size=left_max_nze,&
         error=error)
    left_buffer_calc => left_set
    left_buffer_comm => left_buffer
    ! Right:
    CALL buffer_matrices_ensure_size(right_set, index_size=right_max_nblks, &
         data_size=right_max_nze, error=error)
    CALL setup_buffer_matrices (right_buffer, right_row_nimages, 1,&
         right_set%mats(1,1), index_size=right_max_nblks, data_size=right_max_nze,&
         error=error)
    right_buffer_calc => right_set
    right_buffer_comm => right_buffer
    !
    ! Pre-allocate array for norms
    ! Count the maximum possible multiplies per row for on-the-fly filtering.
    IF (otf_filtering) THEN
       IF (careful_mod) THEN
          CALL dbcsr_assert ("NOT", left_set%mats(1, 1)%m%bcsc,&
               dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
               "Can not do on-the-fly filtering with CSC-indexed matrices.",&
               __LINE__, error=error)
       ENDIF
       nblkrows_total = dbcsr_nblkrows_total (left_set%mats(1, 1))
       ALLOCATE(left_norms(left_max_nblks),right_norms(right_max_nblks),&
                row_max_epss(nblkrows_total),stat=stat)
       CALL dbcsr_assert (stat, "EQ", 0, dbcsr_fatal_level,&
            dbcsr_internal_error,&
            routineN, "Could not allocate memory for left/right norms", __LINE__,&
            error=error)
       IF (dbcsr_has_local_row_index (left_set%mats(1, 1))) THEN
          nblkrows_used = dbcsr_nblkrows_local (left_set%mats(1, 1))
       ELSE
          nblkrows_used = nblkrows_total
       ENDIF
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
       DO v_ki = 1, left_col_nimages
          list_indexing =&
               left_set%mats(1, v_ki)%m%list_indexing
          IF (careful_mod) THEN
             IF (list_indexing) THEN
                CALL dbcsr_assert ((left_set%mats(1, v_ki)%m%nblks)*3, "EQ",&
                     SIZE(left_set%mats(1, v_ki)%m%coo_l),&
                     dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                     "Row count mismatch", __LINE__, error=error)
             ELSE
                CALL dbcsr_assert (nblkrows_used+1, "EQ",&
                     SIZE(left_set%mats(1, v_ki)%m%row_p),&
                     dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                     "Row count mismatch", __LINE__, error=error)
             ENDIF
          ENDIF
          IF (list_indexing) THEN
             CALL count_bins (&
                  left_set%mats(1, v_ki)%m%nblks,&
                  left_set%mats(1, v_ki)%m%coo_l(1::3),&
                  nblkrows_used, row_counts)
          ELSE
             CALL dbcsr_count_row_index (&
                  left_set%mats(1, v_ki)%m%row_p,&
                  row_counts, nblkrows_used)
          ENDIF
          IF (dbcsr_has_local_row_index (left_set%mats(1, v_ki))) THEN
             local_rows => array_data (left_set%mats(1, v_ki)%m%local_rows)
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
       ! The counted blocks are then summed up
       CALL mp_sum(total_row_counts, mp_group)
       ! and used to determine the maximum per-block epsilon.
       DO row = 1 , nblkrows_total
          row_max_epss (row) &
               = REAL(filter_eps&
               / REAL(MAX(1, total_row_counts(row)), KIND=KIND(row_max_epss)),&
               KIND=KIND(row_max_epss))
       END DO
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
    ELSE
       !
       ! These arrays must be valid when passed to called subroutines.
       ALLOCATE(left_norms(0),right_norms(0),&
                row_max_epss(0),stat=stat)
       CALL dbcsr_assert (stat, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Could not allocate memory",&
            __LINE__, error=error)
    ENDIF

!$omp parallel &
!$omp default (none) &
!$omp shared (left_buffer_comm, right_buffer_comm,&
!$omp         product_matrix,&
!$omp         keep_sparsity, filter_eps, row_max_epss, multrec) &
!$omp private(ithread),&
!$omp firstprivate (error)
    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM ()
    ALLOCATE(multrec(ithread)%p)
    CALL dbcsr_mm_multrec_init(multrec(ithread)%p,&
         left_buffer_comm%mats(1, 1)%m,&
         right_buffer_comm%mats(1, 1)%m,&
         product_matrix%m,&
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
    CALL dbcsr_data_init(left_data_get)
    CALL dbcsr_data_init(right_data_get)
    CALL dbcsr_data_new(left_data_get, data_type)
    CALL dbcsr_data_new(right_data_get, data_type)
    
    ! Setup transpose stackbuffers
    IF (mm_driver==mm_driver_acc) THEN
       CALL dbcsr_data_init(trs_stackbuf_1)
       CALL dbcsr_data_init(trs_stackbuf_2)
       CALL dbcsr_data_new(trs_stackbuf_1, data_type=dbcsr_type_int_4, data_size=1000, memory_type=memtype_trsbuffer_1)
       CALL dbcsr_data_new(trs_stackbuf_2, data_type=dbcsr_type_int_4, data_size=1000, memory_type=memtype_trsbuffer_2)
       trs_stackbuf_calc => trs_stackbuf_1
       trs_stackbuf_comm => trs_stackbuf_2
    ENDIF
    !
    ! Pre-CANNON execution:
    ! take communicators for right and left,
    ! set the current position for data.
    ALLOCATE(right_gets_p(0:right_row_nimages-1,0:nsteps_k-1))
    ALLOCATE(left_gets_p(0:left_col_nimages-1,0:nsteps_k-1))
    !
    IF (nsteps_k.GT.1) THEN
       !
       CALL dbcsr_error_set(routineN//"_precannon", error_handler2, error)
       IF (dbcsr_mp_has_subgroups (right_mp_obj)) THEN
          grp_right = dbcsr_mp_my_col_group (right_mp_obj)
          right_gets_p(:,0) = right_myprow
       ELSE
          grp_right = dbcsr_mp_group (right_mp_obj)
          right_gets_p(:,0) = mynode
       ENDIF
       IF (dbcsr_mp_has_subgroups (left_mp_obj)) THEN
          grp_left = dbcsr_mp_my_row_group (left_mp_obj)
          left_gets_p(:,0) = left_mypcol
       ELSE
          grp_left = dbcsr_mp_group (left_mp_obj)
          left_gets_p(:,0) = mynode
       END IF
       !
       ! send/recv requestes for left and right processors and info
       ALLOCATE (precannon_requests(2,0:left_col_nimages+right_row_nimages-1))
       DO metronome = 1, nsteps_k-1
          !
          ! Right matrix transfer
          DO v_ki = 0, right_row_nimages-1
             ! Calculate the process to send to.  It's the virtual
             ! process row -min_nimages up (i.e., smaller row number)
             ! from me.
             CALL image_calculator (right_set%image_dist,&
                  prow=send_prow, pcol=send_pcol,&   ! output
                  ! myvprow goes through all of my (process row) images
                  myvprow=v_ki+right_myfirstvrow,&
                  myvpcol=right_myfirstvcol,& ! nothing happens in the columns
                  ! send to process min_nimages up in the grid
                  vprow_shift=-min_nimages,&
                  shifting='0', error=error)
             ! Calculate the process to receive from, i.e. row +min_nimages down
             CALL image_calculator (right_set%image_dist,&
                  prow=recv_prow, pcol=recv_pcol,&
                  myvprow=v_ki+right_myfirstvrow,&
                  myvpcol=right_myfirstvcol,&
                  vprow_shift=+min_nimages,& ! just the opposite as "send to"
                  shifting='0', error=error)
             !
             ! These are row-communicator relative
             IF (dbcsr_mp_has_subgroups (right_mp_obj)) THEN
                right_send_p = send_prow
                right_recv_p = recv_prow
             ELSE
                right_send_p = right_pgrid (send_prow, send_pcol)
                right_recv_p = right_pgrid (recv_prow, recv_pcol)
             ENDIF
             CALL mp_isendrecv(right_gets_p(v_ki,metronome-1),right_send_p,&
                  right_gets_p(v_ki,metronome),right_recv_p,&
                  grp_right,precannon_requests(1,v_ki),&
                  precannon_requests(2,v_ki))
          END DO
          !
          ! Left matrix transfer
          DO v_ki = 0, left_col_nimages-1
             ! Calculate the process to send to
             CALL image_calculator (left_set%image_dist,&
                  prow=send_prow, pcol=send_pcol,&   ! output
                  myvprow=left_myfirstvrow,& ! nothing happens in the rows
                  ! go through all my column images
                  myvpcol=v_ki+left_myfirstvcol,&
                  ! send to process min_nimages left in the grid
                  vpcol_shift=-min_nimages,&
                  shifting='0', error=error)
             ! Calculate the process to receive from
             CALL image_calculator (left_set%image_dist,&
                  prow=recv_prow, pcol=recv_pcol,&
                  myvprow=left_myfirstvrow,&
                  myvpcol=v_ki+left_myfirstvcol,&
                  vpcol_shift=+min_nimages,& ! just the opposite as "send to"
                  shifting='0', error=error)
             !
             ! These are column-communicator relative
             IF (dbcsr_mp_has_subgroups (left_mp_obj)) THEN
                left_send_p = send_pcol
                left_recv_p = recv_pcol
             ELSE
                left_send_p = left_pgrid (send_prow, send_pcol)
                left_recv_p = left_pgrid (recv_prow, recv_pcol)
             ENDIF
             CALL mp_isendrecv(left_gets_p(v_ki,metronome-1),left_send_p,&
                  left_gets_p(v_ki,metronome),left_recv_p,&
                  grp_left,precannon_requests(1,right_row_nimages+v_ki),&
                  precannon_requests(2,right_row_nimages+v_ki))
          ENDDO
          CALL mp_waitall (precannon_requests)
       ENDDO
       DEALLOCATE(precannon_requests)
       CALL dbcsr_error_stop(error_handler2, error)
       ! 
       ! Set the original data buffers and corresponding windows
       ! Set the right matrix buffers
       !
       CALL dbcsr_data_init(right_data_orig)
       CALL dbcsr_data_new(right_data_orig,data_type=data_type,&
                           data_size=INT(right_data_size),memory_type=memtype_mpi_buffer)
       CALL dbcsr_data_init(right_meta_orig)
       CALL dbcsr_data_new(right_meta_orig,data_type=dbcsr_type_int_4,&
                           data_size=INT(right_meta_size),memory_type=memtype_mpi_buffer)
       DO v_ki = 1, right_row_nimages
          CALL dbcsr_data_copy(right_data_orig,right_set%mats(v_ki, 1)%m%data_area,&
                               INT(local_right_gets(idata_size,v_ki)),&
                               INT(local_right_gets(idata_displ,v_ki))+1)
          CALL dbcsr_data_copy(right_meta_orig,right_set%mats(v_ki, 1)%m%index,&
                               INT(local_right_gets(imeta_size,v_ki)),&
                               INT(local_right_gets(imeta_displ,v_ki))+1)
       ENDDO
       !
       ! Set the left matrix buffers
       !
       CALL dbcsr_data_init(left_data_orig)
       CALL dbcsr_data_new(left_data_orig,data_type=data_type,&
                           data_size=INT(left_data_size),memory_type=memtype_mpi_buffer)       
       CALL dbcsr_data_init(left_meta_orig)
       CALL dbcsr_data_new(left_meta_orig,data_type=dbcsr_type_int_4,&
                           data_size=INT(left_meta_size),memory_type=memtype_mpi_buffer)
       DO v_ki = 1, left_col_nimages
          CALL dbcsr_data_copy(left_data_orig,left_set%mats(1, v_ki)%m%data_area,&
                               INT(local_left_gets(idata_size,v_ki)),&
                               INT(local_left_gets(idata_displ,v_ki))+1)
          CALL dbcsr_data_copy(left_meta_orig,left_set%mats(1, v_ki)%m%index,&
                               INT(local_left_gets(imeta_size,v_ki)),&
                               INT(local_left_gets(imeta_displ,v_ki))+1)
       ENDDO
       !
       ! Create the windows
       ! Right:
       CALL dbcsr_win_create_any(right_data_orig,grp_right,right_data_win,error)
       CALL dbcsr_win_create_any(right_meta_orig,grp_right,right_meta_win,error)
       CALL mp_win_lock_all(right_data_win)
       CALL mp_win_lock_all(right_meta_win)
       ! Left:
       CALL dbcsr_win_create_any(left_data_orig,grp_left,left_data_win,error)
       CALL dbcsr_win_create_any(left_meta_orig,grp_left,left_meta_win,error)
       CALL mp_win_lock_all(left_data_win)
       CALL mp_win_lock_all(left_meta_win)
    ENDIF

    ALLOCATE (right_get_requests(2,0:right_row_nimages-1))
    right_get_requests = mp_request_null
    ALLOCATE (left_get_requests(2,0:left_col_nimages-1))
    left_get_requests = mp_request_null
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
       CALL dbcsr_error_set(routineN//"_metrocomm1", error_handler2, error)
       wait_right: IF (metronome .GT. 1) THEN
          IF (debug_mod) WRITE (*,'(1X,A)')routineN//" waiting for right"
          !
          CALL mp_waitall (right_get_requests)
          right_get_requests = mp_request_null
          ! clean-up wins
          IF (metronome.eq.nsteps_k) THEN
             CALL mp_win_unlock_all(right_data_win)
             CALL mp_win_free(right_data_win)
             CALL mp_win_unlock_all(right_meta_win)
             CALL mp_win_free(right_meta_win)
          ENDIF
       ENDIF wait_right
       CALL dbcsr_error_stop(error_handler2, error)
       ! Right matrix transfer. Transfer in all but the last loop
       ! iteration.
       xfer_right: IF (metronome .LT. nsteps_k) THEN
          DO v_ki = 0, right_row_nimages-1
             ! Calculate which data I receive
             CALL image_calculator (right_set%image_dist,&
                  vprow=get_vrow,vpcol=get_vcol,&
                  myvprow=v_ki+right_myfirstvrow,&
                  myvpcol=right_myfirstvcol,&
                  ! receive window moves with the metronome
                  vprow_shift=metronome*min_nimages,&
                  shifting='R', error=error)
             !
             IF (mm_driver == mm_driver_acc) THEN
                CALL dbcsr_error_set(routineN//"_acc_sync_right", error_handler3, error)
                CALL acc_event_synchronize(right_buffer_comm%mats(v_ki+1, 1)%m%data_area%d%acc_ready)
                CALL dbcsr_error_stop(error_handler3, error)
             ENDIF
             !
             CALL dbcsr_data_set_pointer(&
                  area=right_data_get,&
                  rsize=INT(right_gets(idata_size,get_vrow,get_vcol)),&
                  csize=1,&
                  pointee=right_buffer_comm%mats(v_ki+1, 1)%m%data_area)
             right_meta_get => right_buffer_comm%mats(v_ki+1, 1)%&
                  m%index(1:right_gets(imeta_size,get_vrow,get_vcol))
             !
             CALL dbcsr_error_set(routineN//"_metrocomm2", error_handler2, error)
             IF (right_gets(idata_size,get_vrow,get_vcol).GT.0) THEN
                CALL dbcsr_rget_any(right_data_get,right_gets_p(v_ki,metronome),&
                     right_data_win,right_gets(idata_displ,get_vrow,get_vcol),&
                     right_get_requests(1,v_ki),error=error)
                CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(1,:),right_gets(idata_size,get_vrow,get_vcol),&
                                          dbcsr_mpi_statistics%data_size_breakdown(1,:),data_type_byte)
                dbcsr_mpi_statistics%nexchanged = dbcsr_mpi_statistics%nexchanged+1
             ELSE
                dbcsr_mpi_statistics%nfiltered = dbcsr_mpi_statistics%nfiltered+1
             ENDIF
             CALL mp_rget(right_meta_get,right_gets_p(v_ki,metronome),&
                  right_meta_win,right_gets(imeta_displ,get_vrow,get_vcol),&
                  right_get_requests(2,v_ki))
             CALL count_mpi_statistics(dbcsr_mpi_statistics%index_size(1,:),right_gets(imeta_size,get_vrow,get_vcol))
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
       CALL dbcsr_error_set(routineN//"_metrocomm3", error_handler2, error)
       wait_left: IF (metronome .GT. 1) THEN
          IF (debug_mod) WRITE (*,'(1X,A)')routineN//" waiting for left"
          CALL mp_waitall (left_get_requests)
          left_get_requests = mp_request_null
          ! clean-up wins
          IF (metronome.eq.nsteps_k) THEN
             CALL mp_win_unlock_all(left_data_win)
             CALL mp_win_free(left_data_win)
             CALL mp_win_unlock_all(left_meta_win)
             CALL mp_win_free(left_meta_win)
          ENDIF
       ENDIF wait_left
       CALL dbcsr_error_stop(error_handler2, error)
       ! Left matrix transfer. Transfer in all but the last loop
       ! iteration.
       xfer_left: IF (metronome .LT. nsteps_k) THEN
          DO v_ki = 0, left_col_nimages-1
             ! Calculate which data I receive
             CALL image_calculator (left_set%image_dist,&
                  vprow=get_vrow,vpcol=get_vcol,&
                  myvprow=left_myfirstvrow,&
                  myvpcol=v_ki+left_myfirstvcol,&
                  ! receive window moves with the metronome
                  vpcol_shift=metronome*min_nimages,&
                  shifting='L', error=error)
             !
             IF (mm_driver == mm_driver_acc) THEN
                CALL dbcsr_error_set(routineN//"_acc_sync_left", error_handler3, error)
                CALL acc_event_synchronize(left_buffer_comm%mats(1, v_ki+1)%m%data_area%d%acc_ready)
                CALL dbcsr_error_stop(error_handler3, error)
             ENDIF
             !
             CALL dbcsr_data_set_pointer(&
                  area=left_data_get,&
                  rsize=INT(left_gets(idata_size,get_vrow,get_vcol)),&
                  csize=1,&
                  pointee=left_buffer_comm%mats(1, v_ki+1)%m%data_area)
             left_meta_get => left_buffer_comm%mats(1, v_ki+1)%&
                  m%index(1:left_gets(imeta_size,get_vrow,get_vcol))
             !
             CALL dbcsr_error_set(routineN//"_metrocomm4", error_handler2, error)
             IF (left_gets(idata_size,get_vrow,get_vcol).GT.0) THEN
                CALL dbcsr_rget_any(left_data_get,left_gets_p(v_ki,metronome),&
                     left_data_win,left_gets(idata_displ,get_vrow,get_vcol),&
                     left_get_requests(1,v_ki),error=error)
                CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(2,:),left_gets(idata_size,get_vrow,get_vcol),&
                                          dbcsr_mpi_statistics%data_size_breakdown(2,:),data_type_byte)
                dbcsr_mpi_statistics%nexchanged = dbcsr_mpi_statistics%nexchanged+1
             ELSE
                dbcsr_mpi_statistics%nfiltered = dbcsr_mpi_statistics%nfiltered+1
             ENDIF
             CALL mp_rget(left_meta_get,left_gets_p(v_ki,metronome),&
                  left_meta_win,left_gets(imeta_displ,get_vrow,get_vcol),&
                  left_get_requests(2,v_ki))
             CALL count_mpi_statistics(dbcsr_mpi_statistics%index_size(2,:),left_gets(imeta_size,get_vrow,get_vcol))
             CALL dbcsr_error_stop(error_handler2, error)
          ENDDO
       ENDIF xfer_left
       !
       ! Repoint indices of left matrices
       calc_case_left: IF (metronome .GT. 1) THEN
          DO v_ki = 0, left_col_nimages-1
             CALL dbcsr_repoint_index (left_buffer_calc%mats(1,v_ki+1)%m)
             left_buffer_calc%mats(1, v_ki+1)%m%valid=.TRUE.
          ENDDO
       ENDIF calc_case_left
       !
       ! Do the multiplications.
       !
!$omp parallel default (none) &
!$omp shared (nthreads,min_nimages,error, &
!$omp         left_buffer_calc,right_buffer_calc, &
!$omp         right_set,left_set, &
!$omp         mm_driver,trs_stackbuf_calc, &
!$omp         otf_filtering,left_norms,right_norms, &
!$omp         error_handler2,error_handler4, &           
!$omp         metronome,nsteps_k,multrec) &
!$omp private (ithread,v_ki,t_error) &
!$omp reduction (+: flop)
       ithread = 0
!$     ithread = omp_get_thread_num()
       DO v_ki = 0, min_nimages-1
!$omp master
          IF (debug_mod) THEN
             CALL dbcsr_print(left_buffer_calc%mats(1, v_ki+1), nodata=.TRUE., error=error)
             CALL dbcsr_print(right_buffer_calc%mats(v_ki+1, 1), nodata=.TRUE., error=error)
          ENDIF
          !
          ! form here the code for dbcsr_mm_driver_inner_init was taken 
          ! Since the right matrix is shifted vertically, the
          ! received data always has different notions of "local
          ! rows".  Thus the local_rows and global_rows must be
          ! recalculated.
          CALL dbcsr_reset_vlocals (right_buffer_calc%mats(v_ki+1, 1),&
               right_set%image_dist, error=error)
          ! Since the right matrix is shifted vertically, the
          ! received data always has different notions of "local
          ! rows".  Thus the local_rows and global_rows must be
          ! recalculated.
          CALL dbcsr_reset_vlocals (left_buffer_calc%mats(1, v_ki+1),&
               left_set%image_dist, error=error)

          IF (mm_driver==mm_driver_acc) THEN
             CALL dbcsr_data_host2dev(left_buffer_calc%mats(1, v_ki+1)%m%data_area, error)
             CALL dbcsr_data_host2dev(right_buffer_calc%mats(v_ki+1, 1)%m%data_area, error)
             CALL acc_transpose_blocks(right_buffer_calc%mats(v_ki+1, 1), trs_stackbuf_calc, error)
          END IF

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
!$omp end master
!$omp barrier
          t_error=error
          CALL dbcsr_error_set(routineN//"_multrec", error_handler2, t_error)
          IF(ithread==0)&
               CALL dbcsr_error_set(routineN//"_multrec_master", error_handler4, t_error)

          IF(metronome==nsteps_k .AND. v_ki==min_nimages-1) &
               CALL dbcsr_mm_multrec_phaseout(multrec(ithread)%p, t_error)
          
          CALL dbcsr_mm_multrec_multiply(multrec(ithread)%p,&
               left=left_buffer_calc%mats(1, v_ki+1)%m,&
               right=right_buffer_calc%mats(v_ki+1, 1)%m,&
               flop=flop,&
               a_norms=left_norms, b_norms=right_norms,&
               error=t_error)
          
          IF(ithread==0) CALL dbcsr_error_stop(error_handler4, t_error)
          
          IF(metronome==nsteps_k .AND. v_ki==min_nimages-1) THEN
             CALL dbcsr_mm_multrec_finalize(multrec(ithread)%p, t_error)
             DEALLOCATE(multrec(ithread)%p)
          ENDIF
!$omp barrier
          CALL dbcsr_error_stop(error_handler2, t_error)
       ENDDO
!$omp end parallel

       CALL dbcsr_switch_sets (left_buffer_calc, left_buffer_comm)
       CALL dbcsr_switch_sets (right_buffer_calc, right_buffer_comm)
       CALL dbcsr_switch_d_ptrs(trs_stackbuf_calc, trs_stackbuf_comm)

    ENDDO grouped_k_index

    IF (mm_driver==mm_driver_acc) THEN
       CALL dbcsr_data_release(trs_stackbuf_1)
       CALL dbcsr_data_release(trs_stackbuf_2)
    END IF
    
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
    CALL dbcsr_destroy_array (right_buffer, error=error)
    CALL dbcsr_destroy_array (left_buffer, error=error)
    !
    CALL dbcsr_data_clear_pointer(left_data_get)
    CALL dbcsr_data_clear_pointer(right_data_get)
    CALL dbcsr_data_release(left_data_get)
    CALL dbcsr_data_release(right_data_get)
    !
    IF (nsteps_k.GT.1) THEN
       CALL dbcsr_data_release(left_data_orig)
       CALL dbcsr_data_release(left_meta_orig)
       CALL dbcsr_data_release(right_data_orig)
       CALL dbcsr_data_release(right_meta_orig)
    ENDIF
    !
    DEALLOCATE(left_get_requests, right_get_requests)
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
       checksum = dbcsr_checksum (product_matrix, error=error)
       IF ((output_unit>0)) THEN
          WRITE(output_unit,'(1X,A,1X,F9.4)')"Product Checksum=",checksum
       ENDIF
    ENDIF
    !
    DEALLOCATE(local_right_gets,local_left_gets)
    DEALLOCATE(all_right_gets,all_left_gets)
    DEALLOCATE(right_gets,left_gets)
    DEALLOCATE(right_gets_p,left_gets_p)
    DEALLOCATE(multrec)
    DEALLOCATE(left_buffer, right_buffer)
    !
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE cannon_multiply_low_rma

! *****************************************************************************
!> \brief ...
!> \param mpi_statistics ...
!> \param data_size ...
!> \param size_breakdown ...
!> \param element_size_bytes ...
! *****************************************************************************
  SUBROUTINE count_mpi_statistics(mpi_statistics,data_size,&
                                  size_breakdown,element_size_bytes)
    INTEGER(KIND=int_8), DIMENSION(:), &
      INTENT(INOUT)                          :: mpi_statistics
    INTEGER(KIND=int_8), INTENT(IN)          :: data_size
    INTEGER(KIND=int_8), DIMENSION(:), &
      INTENT(INOUT), OPTIONAL                :: size_breakdown
    INTEGER, INTENT(IN), OPTIONAL            :: element_size_bytes

    INTEGER                                  :: ilimit, nlimits
    INTEGER(KIND=int_8)                      :: data_size_bytes, llimit

    mpi_statistics(1) = mpi_statistics(1) + data_size
    mpi_statistics(2) = MIN(mpi_statistics(2),data_size)
    mpi_statistics(3) = MAX(mpi_statistics(3),data_size)

    IF (PRESENT(size_breakdown)) THEN
       data_size_bytes=data_size
       ! change in bytes
       IF (PRESENT(element_size_bytes)) data_size_bytes=data_size_bytes*element_size_bytes
       nlimits=SIZE(dbcsr_mpi_size_limits)
       ! check for oversize messages
       IF (data_size_bytes.GT.dbcsr_mpi_size_limits(nlimits)) THEN
          size_breakdown(nlimits+1)=size_breakdown(nlimits+1)+1
          RETURN
       ENDIF
       llimit=0
       DO ilimit = 1, nlimits
          IF (data_size_bytes.GE.llimit.AND.data_size_bytes.LE.dbcsr_mpi_size_limits(ilimit)) THEN
             size_breakdown(ilimit)=size_breakdown(ilimit)+1
             RETURN
          ENDIF
          llimit=dbcsr_mpi_size_limits(ilimit)
       ENDDO
    ENDIF

  END SUBROUTINE count_mpi_statistics

! ******************************************************************************
! *****************************************************************************
!> \brief ...
!> \param buffer_set ...
!> \param buff_rows ...
!> \param buff_cols ...
!> \param source_matrix ...
!> \param index_size ...
!> \param data_size ...
!> \param error ...
!>
! *****************************************************************************
  SUBROUTINE setup_buffer_matrices (buffer_set, buff_rows, buff_cols,&
       source_matrix, index_size, data_size, error)
    TYPE(dbcsr_2d_array_type), INTENT(OUT)   :: buffer_set
    INTEGER, INTENT(IN)                      :: buff_rows, buff_cols
    TYPE(dbcsr_obj), INTENT(IN)              :: source_matrix
    INTEGER(kind=mp_address_kind), &
      INTENT(IN)                             :: index_size, data_size
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_buffer_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col_image, error_handler, &
                                                row_image
    INTEGER, DIMENSION(:), POINTER           :: i1

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

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
               nblks = INT(index_size), nze = INT(data_size), &
               data_memory_type = memtype_abpanel_2,&
               index_memory_type = dbcsr_memtype_default,&
               error = error)
          ! Dummy allocation only needed for NAG (at least for 5.1(327))
          DEALLOCATE(i1)
          CALL dbcsr_data_ensure_size (&
               buffer_set%mats(row_image, col_image)%m%data_area,&
               INT(data_size), nocopy=.TRUE.,error=error)
          CALL ensure_array_size (&
               buffer_set%mats(row_image, col_image)%m%index,&
               ub=INT(index_size), nocopy=.TRUE.,&
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


! ******************************************************************************
!> \brief Enlarge left_set and right_set to hold any a/b-panel.
!>  left_set and right_set are created by make_images to hold the a/b-panels
!>  used for the inital cannon-tick. This routine ensures that these buffers
!>  can hold any a/b-panel occuring during a matrix multiply and makes them
!>  therefore suitable as buffers for the entire cannon algorithm.
!>  This saves memory since no seperate buffers for the first cannon-tick
!>  have to be allocated.
!>  
!> \param buffer_set ...
!> \param index_size ...
!> \param data_size ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE buffer_matrices_ensure_size(buffer_set, index_size, data_size, error)
    TYPE(dbcsr_2d_array_type), INTENT(INOUT) :: buffer_set
    INTEGER(kind=mp_address_kind), &
      INTENT(IN)                             :: index_size, data_size
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'buffer_matrices_ensure_size', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col_image, error_handler, &
                                                row_image

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    DO row_image = 1, SIZE(buffer_set%mats, 1)
       DO col_image = 1, SIZE(buffer_set%mats, 2)
          CALL dbcsr_data_ensure_size (&
               buffer_set%mats(row_image, col_image)%m%data_area,&
               INT(data_size), error=error)
          CALL ensure_array_size (&
               buffer_set%mats(row_image, col_image)%m%index,&
               ub=INT(index_size),&
               memory_type=dbcsr_get_index_memory_type(buffer_set%mats(row_image, col_image)),&
               error=error)
          CALL dbcsr_repoint_index(buffer_set%mats(row_image, col_image)%m)
       ENDDO
    ENDDO
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE buffer_matrices_ensure_size


! ******************************************************************************
! *****************************************************************************
!> \brief ...
!> \param matrix_set ...
!> \param error ...
!>
! *****************************************************************************
  SUBROUTINE setup_rec_index (matrix_set, error)
    TYPE(dbcsr_2d_array_type), INTENT(INOUT) :: matrix_set
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_rec_index', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handler, i_col, i_row, &
                                                t_f, t_l, t_size

!$  INTEGER                                  :: ithread
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
!> \param m ...
!> \param n ...
!> \param nblks ...
!> \param idx ...
!> \param error ...
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
!> \param mi ...
!> \param mf ...
!> \param ni ...
!> \param nf ...
!> \param nele ...
!> \param a ...
!> \param d ...
!> \par History
!> - 2011-02-17 [UB] modified for use in DBCSR; reduced memory usage.
!> \author JV
!> \note Always cut longest first. On a tie cut N
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
       CALL rec_split (nele, a, tmp, 1, nlow, mi, half)
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
       CALL rec_split (nele, a, tmp, 2, nlow, ni, half)
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


! *****************************************************************************
!> \brief ...
!> \param nele ...
!> \param a ...
!> \param split ...
!> \param row_or_col ...
!> \param nlow ...
!> \param mi ...
!> \param half ...
! *****************************************************************************
  SUBROUTINE rec_split (nele, a, split, row_or_col, nlow, mi, half)
    INTEGER, INTENT(IN)                      :: nele
    INTEGER, DIMENSION(3, nele), INTENT(IN)  :: a
    INTEGER, DIMENSION(3, nele), INTENT(OUT) :: split
    INTEGER, INTENT(IN)                      :: row_or_col
    INTEGER, INTENT(OUT)                     :: nlow
    INTEGER, INTENT(IN)                      :: mi, half

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
!> \brief Switches pointers between two matrix sets
!> \param[in,out] set1p ...
!> \param[in,out] set2p ...
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
!> \brief Switches pointers between two data areas
!> \param area1p ...
!> \param area2p ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_switch_d_ptrs (area1p, area2p)
    TYPE(dbcsr_data_obj), POINTER            :: area1p, area2p

    TYPE(dbcsr_data_obj), POINTER            :: tmp_p

    tmp_p  => area1p
    area1p => area2p
    area2p => tmp_p
  END SUBROUTINE dbcsr_switch_d_ptrs

! *****************************************************************************
! The following routines are helped here to help the compiler optimize them
! out.
! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param t ...
!> \retval blas_mat_type ...
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

! *****************************************************************************
!> \brief ...
!> \param t ...
!> \retval flip_type ...
! *****************************************************************************
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
                                                nblks, nrows
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


! *****************************************************************************
!> \brief write out a stack for transposing the blocks
!> \param matrix ...
!> \param trs_stackbuf ...
!> \param error ...
!>  \author Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_transpose_blocks(matrix, trs_stackbuf, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: trs_stackbuf
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'acc_transpose_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_p, col, error_handle, i, &
                                                m, mi, mi_max, n, nblks, ni, &
                                                ni_max, offset, row, x
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: counters, filled, offsets
    LOGICAL                                  :: local
    INTEGER, DIMENSION(:), POINTER :: blk_index, col_blk_sizes, &
      col_blk_sizes2enum, enum2col_blk_sizes, enum2row_blk_sizes, &
      local2global_cols, local2global_rows, row_blk_sizes, &
      row_blk_sizes2enum, trs_stack
    INTEGER, DIMENSION(3, matrix%m%nblks)    :: tmp_stack

    CALL dbcsr_error_set(routineN, error_handle, error)

    NULLIFY(row_blk_sizes2enum, enum2row_blk_sizes)
    NULLIFY(col_blk_sizes2enum, enum2col_blk_sizes)
    NULLIFY(local2global_rows, local2global_cols, trs_stack)

    IF (.NOT. matrix%m%list_indexing) STOP "build_trs_stack: only list_indexing supported."
    IF (trs_stackbuf%d%data_type /= dbcsr_type_int_4) STOP "build_trs_stack: stac_buf has wrong datatype"
    blk_index => matrix%m%coo_l
    row_blk_sizes => array_data (matrix%m%row_blk_size)
    col_blk_sizes => array_data (matrix%m%col_blk_size)
    local = matrix%m%local_indexing
    local2global_rows => array_data (matrix%m%local_rows)
    local2global_cols => array_data (matrix%m%local_cols)
    nblks = matrix%m%nblks

    ! enumerate the blocksizes to keep the following 2D-arrays small.
    CALL enumerate_blk_sizes(row_blk_sizes, row_blk_sizes2enum, enum2row_blk_sizes)
    CALL enumerate_blk_sizes(col_blk_sizes, col_blk_sizes2enum, enum2col_blk_sizes)
    mi_max = SIZE(enum2row_blk_sizes); ni_max = SIZE(enum2col_blk_sizes)
    ALLOCATE(counters(mi_max, ni_max), offsets(mi_max, ni_max), filled(mi_max, ni_max))
    counters(:,:)=0;  offsets(:,:)=0; filled(:,:)=0

    ! make sure buffer from previous cannon-tick was uploaded
    CALL acc_event_synchronize(trs_stackbuf%d%acc_ready)

    CALL dbcsr_data_ensure_size(trs_stackbuf, data_size=nblks, nocopy=.TRUE., error=error)
    trs_stack => trs_stackbuf%d%i4

    ! collect block addresses and dimensions in a temporary stack
    ! while doing so, also count number of blocks per block-dimensions
    DO i = 1, nblks
       row   = blk_index(3*(i-1) + 1)
       col   = blk_index(3*(i-1) + 2)
       blk_p = blk_index(3*(i-1) + 3)
       IF (blk_p == 0) CYCLE
       IF (local) THEN
          row = local2global_rows(row)
          col = local2global_cols(col)
       ENDIF
       m = row_blk_sizes(row)
       n = col_blk_sizes(col)
       mi = row_blk_sizes2enum(m)
       ni = col_blk_sizes2enum(n)
       tmp_stack(1, i) = mi
       tmp_stack(2, i) = ni
       tmp_stack(3, i) = blk_p - 1
       counters(mi,ni) = counters(mi,ni) + 1
    ENDDO

    ! calculate offsets for first element of each sub-stack
    offset = 0
    DO mi=1, mi_max
      DO ni=1, ni_max
        offsets(mi, ni) = offset
        offset = offset + counters(mi, ni)
      ENDDO
    ENDDO

    ! write all sub-stacks into the host-pinned buffer
    DO i=1, nblks
      mi     = tmp_stack(1, i)
      ni     = tmp_stack(2, i)
      blk_p  = tmp_stack(3, i)
      x = offsets(mi,ni) + filled(mi,ni) + 1
      trs_stack(x) = blk_p
      filled(mi,ni) = filled(mi,ni) + 1
    ENDDO

    !sanity check
    DO mi=1, mi_max
      DO ni=1, ni_max
         IF(filled(mi,ni) /= counters(mi,ni)) STOP "acc_transpose_blocks: bug"
      END DO
    END DO

    !transfer all stacks
    CALL dbcsr_data_host2dev(trs_stackbuf, error)

    ! make sure block-buffer is uploaded befor running the kernels
    CALL acc_stream_wait_event(trs_stackbuf%d%memory_type%acc_stream, matrix%m%data_area%d%acc_ready)

    ! launch kernels
    DO mi=1, mi_max
      DO ni=1, ni_max
        IF(counters(mi, ni) > 0) THEN
          m = enum2row_blk_sizes(mi)
          n = enum2col_blk_sizes(ni)
          CALL dbcsr_acc_transpose( &
            trs_stack=trs_stackbuf%d%acc_devmem,&
            offset=offsets(mi, ni),&
            nblks=counters(mi, ni),&
            datatype=matrix%m%data_type,&
            buffer=matrix%m%data_area%d%acc_devmem,&
            m=m, n=n,&
            stream=trs_stackbuf%d%memory_type%acc_stream, error=error)
        END IF
      ENDDO
    ENDDO

    ! make sure block-buffer are not used until transpose kernels finnished
    CALL acc_event_record(trs_stackbuf%d%acc_ready, trs_stackbuf%d%memory_type%acc_stream)
    CALL acc_stream_wait_event(matrix%m%data_area%d%memory_type%acc_stream, trs_stackbuf%d%acc_ready)
    CALL acc_event_record(matrix%m%data_area%d%acc_ready, matrix%m%data_area%d%memory_type%acc_stream)

    DEALLOCATE(row_blk_sizes2enum, enum2row_blk_sizes)
    DEALLOCATE(col_blk_sizes2enum, enum2col_blk_sizes)
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE acc_transpose_blocks


! *****************************************************************************
!> \brief Enumerate all occuring blocksizes
!> \param blk_sizes ...
!> \param enum ...
!> \param rev_enum ...
!>  \author Ole Schuett
! *****************************************************************************
  SUBROUTINE enumerate_blk_sizes(blk_sizes, enum, rev_enum)
    INTEGER, DIMENSION(:), POINTER           :: blk_sizes, enum, rev_enum

    INTEGER                                  :: i, n

     n = MAXVAL(blk_sizes)
     ALLOCATE(enum(0:n))
     enum(:) = 0

     DO i=1, SIZE(blk_sizes)
       enum(blk_sizes(i)) = 1
     ENDDO

     n = SUM(enum)
     ALLOCATE(rev_enum(n))

     n = 0
     DO i=0, SIZE(enum)-1
       IF(enum(i) > 0) THEN
           n = n + 1
           enum(i) = n
           rev_enum(n) = i
       END IF
     ENDDO
  END SUBROUTINE enumerate_blk_sizes


#include "dbcsr_mm_cannon_d.f90"
#include "dbcsr_mm_cannon_z.f90"
#include "dbcsr_mm_cannon_s.f90"
#include "dbcsr_mm_cannon_c.f90"

END MODULE dbcsr_mm_cannon

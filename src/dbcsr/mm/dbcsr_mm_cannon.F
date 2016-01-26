!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
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
!>  - 2015-01    Introducting RMA (Alfio Lazzaro)
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
                                             dbcsr_block_transpose_aa,&
                                             dbcsr_data_set
  USE dbcsr_config,                    ONLY: default_resize_factor,&
                                             has_acc,&
                                             use_comm_thread,&
                                             use_mpi_filtering,&
                                             use_mpi_rma
  USE dbcsr_data_methods,              ONLY: &
       dbcsr_data_clear_pointer, dbcsr_data_ensure_size, dbcsr_data_get_size, &
       dbcsr_data_get_size_referenced, dbcsr_data_hold, dbcsr_data_host2dev, &
       dbcsr_data_init, dbcsr_data_new, dbcsr_data_release, &
       dbcsr_data_set_pointer, dbcsr_data_set_size_referenced, &
       dbcsr_get_data_p_c, dbcsr_get_data_p_d, dbcsr_get_data_p_s, &
       dbcsr_get_data_p_z, dbcsr_scalar, dbcsr_scalar_are_equal, &
       dbcsr_scalar_fill_all, dbcsr_scalar_negative, dbcsr_scalar_one, &
       dbcsr_scalar_set_type, dbcsr_scalar_zero
  USE dbcsr_data_types,                ONLY: dbcsr_datatype_sizeof
  USE dbcsr_dist_methods,              ONLY: &
       dbcsr_distribution_col_dist, dbcsr_distribution_get_num_images_1d, &
       dbcsr_distribution_has_col_clusters, &
       dbcsr_distribution_has_row_clusters, dbcsr_distribution_has_threads, &
       dbcsr_distribution_hold, dbcsr_distribution_local_cols, &
       dbcsr_distribution_local_rows, dbcsr_distribution_make_threads, &
       dbcsr_distribution_mp, dbcsr_distribution_ncols, &
       dbcsr_distribution_no_threads, dbcsr_distribution_nrows, &
       dbcsr_distribution_num_local_col_clusters, &
       dbcsr_distribution_num_local_row_clusters, dbcsr_distribution_release, &
       dbcsr_distribution_row_dist, dbcsr_distribution_thread_dist
  USE dbcsr_dist_operations,           ONLY: dbcsr_create_image_dist,&
                                             dbcsr_get_local_vcols,&
                                             dbcsr_get_local_vrows,&
                                             dbcsr_make_dists_dense,&
                                             dbcsr_reset_locals,&
                                             dbcsr_reset_vlocals,&
                                             image_calculator,&
                                             make_sizes_dense
  USE dbcsr_index_operations,          ONLY: &
       dbcsr_count_row_index, dbcsr_has_local_row_index, dbcsr_index_compact, &
       dbcsr_index_prune_deleted, dbcsr_make_index_canonical, &
       dbcsr_make_index_list, dbcsr_make_index_local_row, &
       dbcsr_repoint_index, dbcsr_sort_indices
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
       dbcsr_row_block_sizes, dbcsr_valid_index
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
                                             dbcsr_sendrecv_any,&
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
       dbcsr_1d_array_type, dbcsr_2d_array_type, dbcsr_conjugate_transpose, &
       dbcsr_data_obj, dbcsr_distribution_obj, dbcsr_imagedistribution_obj, &
       dbcsr_iterator, dbcsr_memtype_type, dbcsr_meta_size, dbcsr_mp_obj, &
       dbcsr_mpi_size_limits, dbcsr_mpi_statistics_type, dbcsr_no_transpose, &
       dbcsr_num_slots, dbcsr_obj, dbcsr_scalar_type, dbcsr_slot_blk_p, &
       dbcsr_slot_col_i, dbcsr_slot_coo_l, dbcsr_slot_dense, &
       dbcsr_slot_home_coli, dbcsr_slot_home_pcol, dbcsr_slot_home_prow, &
       dbcsr_slot_home_rowi, dbcsr_slot_home_vpcol, dbcsr_slot_home_vprow, &
       dbcsr_slot_nblkcols_local, dbcsr_slot_nblkcols_total, &
       dbcsr_slot_nblkrows_local, dbcsr_slot_nblkrows_total, &
       dbcsr_slot_nblks, dbcsr_slot_nfullcols_local, &
       dbcsr_slot_nfullcols_total, dbcsr_slot_nfullrows_local, &
       dbcsr_slot_nfullrows_total, dbcsr_slot_nze, dbcsr_slot_row_p, &
       dbcsr_slot_size, dbcsr_slot_thr_c, dbcsr_slot_type, dbcsr_transpose, &
       dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_complex_4, &
       dbcsr_type_complex_8, dbcsr_type_int_4, dbcsr_type_no_symmetry, &
       dbcsr_type_real_4, dbcsr_type_real_8
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
  USE machine,                         ONLY: default_output_unit,&
                                             m_memory
  USE message_passing,                 ONLY: &
       mp_allgather, mp_alltoall, mp_irecv, mp_isend, mp_max, mp_min, &
       mp_request_null, mp_rget, mp_sendrecv, mp_sum, mp_testany, mp_waitall, &
       mp_win_create, mp_win_free, mp_win_lock_all, mp_win_unlock_all
#include "../../base/base_uses.f90"

  !$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads,&
  !$                    omp_set_lock, omp_unset_lock, omp_init_lock, omp_lock_kind, omp_destroy_lock

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_cannon'
  CHARACTER(len=*), PARAMETER, PRIVATE :: int_print = "(10(1X,I7))"
  LOGICAL, PARAMETER :: debug_mod  = .FALSE.
  LOGICAL, PARAMETER :: careful_mod = .FALSE.
  REAL, PARAMETER :: huge_norm = HUGE(1.0)**(1.0/3.0)

  TYPE dbcsr_memtype_type_p
     TYPE(dbcsr_memtype_type), POINTER :: p => Null()
     ! ensure that array-elements are on different cache lines
     INTEGER(kind=int_4), DIMENSION(64)    :: padding
  END TYPE dbcsr_memtype_type_p

  TYPE(dbcsr_memtype_type_p), DIMENSION(:), POINTER, SAVE :: memtype_product_wm => Null()

  TYPE(dbcsr_mpi_statistics_type), PRIVATE, SAVE :: dbcsr_mpi_statistics
  INTEGER,                         PRIVATE, SAVE :: max_nblocks = 0
  REAL,                            PRIVATE, SAVE :: marketing_flops = 0
  REAL,                            PRIVATE, SAVE :: max_memory = 0
  
  INTEGER, PRIVATE, SAVE :: last_mpi_ranks_used = 0

  TYPE(dbcsr_memtype_type),     PRIVATE, SAVE  :: memtype_abpanel_1, memtype_abpanel_2,&
                                                  memtype_trsbuffer_1, memtype_trsbuffer_2, &
                                                  memtype_mpi_buffer
  TYPE(acc_stream_type), PRIVATE, SAVE         :: stream_1, stream_2
  ! ab-panels and streams are shared between all threads

  TYPE dbcsr_rma_buffers
     TYPE(dbcsr_data_obj)               :: left_data, right_data
     INTEGER, DIMENSION(:), POINTER     :: left_meta => Null(), right_meta => Null(), &
                                           left_offset_data => Null(), right_offset_data => Null(), &
                                           left_img_map_rows => Null(), left_img_map_cols => Null(), &
                                           right_img_map_rows => Null(), right_img_map_cols => Null()
     ! count and offset per each image
     INTEGER, DIMENSION(:,:), POINTER   :: left_img_refs_rows => Null(), &
                                           left_img_refs_cols => Null(), &
                                           right_img_refs_rows => Null(), &
                                           right_img_refs_cols => Null()
     INTEGER, DIMENSION(:, :, :), POINTER ::  left_local_refs => Null(), right_local_refs => Null()
     TYPE(dbcsr_1d_array_type)          :: left_buffer, right_buffer
     TYPE(dbcsr_data_obj)               :: trs_stackbuf
     LOGICAL                            :: is_valid_left = .FALSE., is_valid_right = .FALSE.
  END TYPE dbcsr_rma_buffers

  ! RMA buffers
  TYPE(dbcsr_rma_buffers), PRIVATE, TARGET, SAVE :: rma_buffers_orig, &
                                                    rma_buffers_1, rma_buffers_2
  
  INTEGER, PARAMETER , PRIVATE               :: idata_displ = 3, &
                                                idata_size = 1, &
                                                imeta_displ = 4, &
                                                imeta_size = 2, &
                                                ilocal_proc = 1, &
                                                isym_proc = 2

  LOGICAL                                :: my_use_mpi_rma

  TYPE dbcsr_mm_multrec_type_p
    TYPE(dbcsr_mm_multrec_type), POINTER :: p => Null()
    ! ensure that array-elements are on different cache lines
    INTEGER(kind=int_4), DIMENSION(64)       :: padding
  END TYPE dbcsr_mm_multrec_type_p

  INTERFACE dbcsr_switch
     MODULE PROCEDURE dbcsr_switch_sets
     MODULE PROCEDURE dbcsr_switch_d_ptrs
  END INTERFACE

  INTERFACE setup_rec_index
     MODULE PROCEDURE setup_rec_index_1d
     MODULE PROCEDURE setup_rec_index_2d
  END INTERFACE

  PUBLIC :: dbcsr_mm_cannon_lib_init, dbcsr_mm_cannon_lib_finalize
  PUBLIC :: dbcsr_mm_cannon_clear_mempools
  PUBLIC :: dbcsr_mm_cannon_multiply

  CONTAINS

! *****************************************************************************
!> \brief Initialize the library
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_cannon_lib_init()

    INTEGER                                  :: ithread, nthreads

    nthreads = 1; ithread = 0
    !$ nthreads = OMP_GET_NUM_THREADS () ; ithread = OMP_GET_THREAD_NUM ()

    CALL dbcsr_mm_multrec_lib_init()

    !$OMP MASTER
    dbcsr_mpi_statistics%nimages = -1
    dbcsr_mpi_statistics%nexchanged = 0
    dbcsr_mpi_statistics%nfiltered = 0
    dbcsr_mpi_statistics%data_size = 0
    dbcsr_mpi_statistics%data_size(:,2) = HUGE(dbcsr_mpi_statistics%data_size(1,2))
    dbcsr_mpi_statistics%data_size_breakdown = 0

    marketing_flops = 0
    max_memory = 0
    ALLOCATE(memtype_product_wm(0:nthreads-1))
    !$OMP END MASTER
    !$OMP BARRIER

    ! Each thread has its own working-matrix and its own mempool
    ALLOCATE(memtype_product_wm(ithread)%p)
    CALL dbcsr_memtype_setup(memtype_product_wm(ithread)%p, has_pool=.TRUE.)
    CALL dbcsr_mempool_ensure_capacity(memtype_product_wm(ithread)%p%pool, capacity=1)
  END SUBROUTINE dbcsr_mm_cannon_lib_init


! *****************************************************************************
!> \brief Finalize the library
!> \param group ...
!> \param output_unit ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_cannon_lib_finalize(group, output_unit)
    INTEGER, INTENT(IN)                      :: group, output_unit

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_cannon_lib_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilimit, isqrt, isqrt2, ithread
    INTEGER(KIND=int_8)                      :: total_nexchanged, &
                                                total_nfiltered
    INTEGER(KIND=int_8), DIMENSION(SIZE(&
      dbcsr_mpi_size_limits)+1, 2, 2)        :: total_recv_breakdown
    REAL                                     :: average, &
                                                total_marketing_flops, &
                                                total_max_memory
    REAL, DIMENSION(2)                       :: max_recv_data, min_recv_data, &
                                                total_recv_data

     CALL dbcsr_mm_multrec_lib_finalize(group, output_unit)

     ithread = 0
     !$  ithread = omp_get_thread_num()

     ! Each thread has its own working-matrix and its own mempool
     IF (ASSOCIATED(memtype_product_wm(ithread)%p%pool)) &
        CALL dbcsr_mempool_destruct(memtype_product_wm(ithread)%p%pool)
     DEALLOCATE(memtype_product_wm(ithread)%p)

     !$OMP BARRIER
     !$omp master
     DEALLOCATE(memtype_product_wm)

     ! Deallocate RMA buffers
     CALL rma_buffers_release(rma_buffers_orig)
     CALL rma_buffers_release(rma_buffers_1)
     CALL rma_buffers_release(rma_buffers_2)

     total_max_memory=max_memory
     CALL mp_max(total_max_memory,group)

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
        average = SUM(total_recv_data(:))/REAL(total_nexchanged)
     ELSE
        average=0
        min_recv_data=0
     ENDIF

     total_recv_breakdown(:,:,:)=dbcsr_mpi_statistics%data_size_breakdown(:,:,:)
     CALL mp_sum(total_recv_breakdown,group)

     IF (output_unit>0) THEN
       WRITE (output_unit,'(A,T30,EN20.6)') " marketing flops", total_marketing_flops

       IF (dbcsr_mpi_statistics%nimages.GT.0) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A)") REPEAT("-",79)
          WRITE (output_unit,'(A,T30,EN20.6)') " max memory usage/rank",total_max_memory
          WRITE (output_unit,'(A,T30,I20)') " # max total images/rank",dbcsr_mpi_statistics%nimages
          WRITE (output_unit,'(A,T30,I20)') " # MPI messages exchanged",total_nexchanged
          WRITE (output_unit,'(A,T30,I20)') " # MPI messages filtered",total_nfiltered
          WRITE (output_unit,'(A)') " MPI messages size (elements):"
          WRITE (output_unit,'(A,T30,EN20.6)') "  total size",&
               SUM(total_recv_data(:))
          WRITE (output_unit,'(A,T30,EN20.6)') "  min size",&
               MINVAL(min_recv_data(:))
          WRITE (output_unit,'(A,T30,EN20.6)') "  max size",&
               MAXVAL(max_recv_data(:))
          WRITE (output_unit,'(A,T30,EN20.6)') "  average size",average

          WRITE (output_unit,'(A)') " MPI breakdown and total messages size (bytes):"
          WRITE (output_unit,'(A,I8,T40,I10,T55,I20)') "             size <= ",dbcsr_mpi_size_limits(1),&
               SUM(total_recv_breakdown(1,1,:)),SUM(total_recv_breakdown(1,2,:))
          DO ilimit=2,SIZE(dbcsr_mpi_size_limits)
             WRITE (output_unit,'(A,I8,A,I8,T40,I10,T55,I20)') "  ",dbcsr_mpi_size_limits(ilimit-1),&
                  " < size <= ",dbcsr_mpi_size_limits(ilimit),&
                  SUM(total_recv_breakdown(ilimit,1,:)),SUM(total_recv_breakdown(ilimit,2,:))
          ENDDO
          ilimit=SIZE(dbcsr_mpi_size_limits)
          WRITE (output_unit,'(A,I8,A,T40,I10,T55,I20)') "  ",dbcsr_mpi_size_limits(ilimit),&
               " < size    ",SUM(total_recv_breakdown(ilimit+1,1,:)),SUM(total_recv_breakdown(ilimit+1,2,:))
       ENDIF

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
        CALL dbcsr_mempool_destruct(memtype_trsbuffer_1%pool)
     IF (ASSOCIATED(memtype_trsbuffer_2%pool)) &
        CALL dbcsr_mempool_destruct(memtype_trsbuffer_2%pool)
     IF (ASSOCIATED(memtype_abpanel_1%pool)) &
        CALL dbcsr_mempool_destruct(memtype_abpanel_1%pool)
     IF (ASSOCIATED(memtype_abpanel_2%pool)) &
        CALL dbcsr_mempool_destruct(memtype_abpanel_2%pool)
     IF(acc_stream_associated(stream_1)) &
        CALL acc_stream_destroy(stream_1)
     IF(acc_stream_associated(stream_2)) &
        CALL acc_stream_destroy(stream_2)
     !$omp end master
  END SUBROUTINE dbcsr_mm_cannon_lib_finalize

! *****************************************************************************
!> \brief Deallocate memory contained in mempools
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mm_cannon_clear_mempools()

    INTEGER                                  :: ithread

     ithread = 0
     !$  ithread = omp_get_thread_num()

     ! Each thread has its own working-matrix and its own mempool
     IF (ASSOCIATED(memtype_product_wm(ithread)%p%pool)) &
        CALL dbcsr_mempool_clear(memtype_product_wm(ithread)%p%pool)

     !$omp master
     IF (ASSOCIATED(memtype_trsbuffer_1%pool)) &
        CALL dbcsr_mempool_clear(memtype_trsbuffer_1%pool)
     IF (ASSOCIATED(memtype_trsbuffer_2%pool)) &
        CALL dbcsr_mempool_clear(memtype_trsbuffer_2%pool)
     IF (ASSOCIATED(memtype_abpanel_1%pool)) &
        CALL dbcsr_mempool_clear(memtype_abpanel_1%pool)
     IF (ASSOCIATED(memtype_abpanel_2%pool)) &
        CALL dbcsr_mempool_clear(memtype_abpanel_2%pool)
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
       flop)

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
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_cannon_multiply', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE., &
                                                use_list_indexing = .TRUE., &
                                                use_local_indexing = .TRUE.
    REAL(real_8), PARAMETER                  :: make_dense_occ_thresh = 1.0_dp

    CHARACTER                                :: transa_l, transb_l
    INTEGER :: comm, f_col, f_k, f_row, handle, handle2, ithread, l_col, l_k, &
      l_row, nimages_left_rows, nimages_match, nimages_right_cols, numnodes, &
      output_unit
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

    CALL timeset(routineN, handle)
    CALL array_nullify(dense_k_sizes)
    CALL array_nullify(dense_col_sizes)
    CALL array_nullify(dense_row_sizes)

    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM ()

    ! setup driver-dependent memory-types and their memory-pools ---------------

    ! the ab_buffers are shared by all threads
    IF (has_acc) THEN
       IF(.NOT. acc_stream_associated(stream_1)) THEN
          CALL acc_stream_create(stream_1, "MemCpy (odd ticks)")
          CALL acc_stream_create(stream_2, "MemCpy (even ticks)")
       ENDIF

       CALL dbcsr_memtype_setup(memtype_abpanel_1, has_pool=.TRUE.,&
            acc_hostalloc=.TRUE., acc_devalloc=.TRUE., acc_stream=stream_1,&
            mpi=.TRUE., oversize_factor=default_resize_factor)

       CALL dbcsr_memtype_setup(memtype_abpanel_2, has_pool=.TRUE.,&
            acc_hostalloc=.TRUE., acc_devalloc=.TRUE., acc_stream=stream_2,&
            mpi=.TRUE., oversize_factor=default_resize_factor)

       !TODO: ensure capacity 2/3?
       CALL dbcsr_memtype_setup(memtype_trsbuffer_1,has_pool=.TRUE.,&
            acc_hostalloc=.TRUE., acc_devalloc=.TRUE., acc_stream=stream_1)
       CALL dbcsr_mempool_ensure_capacity(memtype_trsbuffer_1%pool, capacity=1)

       CALL dbcsr_memtype_setup(memtype_trsbuffer_2, has_pool=.TRUE.,&
            acc_hostalloc=.TRUE., acc_devalloc=.TRUE.,acc_stream=stream_2)
       CALL dbcsr_mempool_ensure_capacity(memtype_trsbuffer_2%pool, capacity=1)
    ELSE
       CALL dbcsr_memtype_setup(memtype_abpanel_1, mpi=.TRUE.)
       CALL dbcsr_memtype_setup(memtype_abpanel_2, mpi=.TRUE.)
    ENDIF

    CALL dbcsr_memtype_setup(memtype_mpi_buffer, mpi=.TRUE.)

    ! check parameters ---------------------------------------------------------
    transa_l = transa
    transb_l = transb
    CALL uppercase(transa_l)
    CALL uppercase(transb_l)
    IF(transa_l.NE.dbcsr_no_transpose.AND.&
       transa_l.NE.dbcsr_transpose.AND.&
       transa_l.NE.dbcsr_conjugate_transpose)&
       CPABORT("Invalid transa_l = "//transa_l)

    IF(transb_l.NE.dbcsr_no_transpose.AND.&
       transb_l.NE.dbcsr_transpose.AND.&
       transb_l.NE.dbcsr_conjugate_transpose)&
       CPABORT("Invalid transb_l = "//transb_l)

    IF (dbg) THEN
       WRITE(*,*)'========== MULTIPLICATION ========================'
       CALL dbcsr_verify_matrix (matrix_a)
       CALL dbcsr_verify_matrix (matrix_b)
       CALL dbcsr_verify_matrix (matrix_c)
       WRITE(*,*)routineN//" ABC checksums",&
            dbcsr_checksum(matrix_a),&
            dbcsr_checksum(matrix_b),&
            dbcsr_checksum(matrix_c)
       IF (dbg) THEN
          CALL dbcsr_print (matrix_a, nodata=.TRUE.)
          CALL dbcsr_print (matrix_b, nodata=.TRUE.)
          CALL dbcsr_print (matrix_c, nodata=.TRUE.)
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
               transpose_distribution=.FALSE.)
       ELSE
          CALL dbcsr_new_transposed (matrix_left, matrix_a,&
               shallow_data_copy=.TRUE., redistribute=.FALSE.,&
               transpose_distribution=.FALSE.)
       ENDIF
       new_left = .TRUE.
    CASE(dbcsr_conjugate_transpose)
       CALL dbcsr_init(matrix_left)
       CALL dbcsr_new_transposed (matrix_left, matrix_a,&
            shallow_data_copy=.FALSE., redistribute=.FALSE.,&
            transpose_distribution=.FALSE.)
       CALL dbcsr_conjg(matrix_left)
       new_left = .TRUE.
    CASE DEFAULT
       CPABORT("wrong transa_l = "//transa_l)
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
               transpose_distribution=.FALSE.)
       ELSE
          CALL dbcsr_new_transposed (matrix_right, matrix_b,&
               shallow_data_copy=.TRUE., redistribute=.FALSE.,&
               transpose_distribution=.FALSE.)
       ENDIF
       new_right = .TRUE.
    CASE(dbcsr_conjugate_transpose)
       CALL dbcsr_init(matrix_right)
       CALL dbcsr_new_transposed (matrix_right, matrix_b,&
            shallow_data_copy=.FALSE., redistribute=.FALSE.,&
            transpose_distribution=.FALSE.)
       CALL dbcsr_conjg(matrix_right)
       new_right = .TRUE.
    CASE DEFAULT
       CPABORT("wrong transb_l = "//transb_l)
    END SELECT
    !
    ! Ensure matrix compatibility.
    IF(.NOT. array_equality(dbcsr_row_block_offsets (matrix_c), dbcsr_row_block_offsets (matrix_left)))&
       CPABORT("C/A rows not equal")
    IF(.NOT.array_equality(dbcsr_col_block_offsets (matrix_c), dbcsr_col_block_offsets (matrix_right)))&
       CPABORT("C/B columns not equal")
    IF(.NOT.array_equality (dbcsr_col_block_offsets (matrix_left),dbcsr_row_block_offsets (matrix_right)))&
       CPABORT("A cols/B rows not equal")
    !
    ! No dense multiplication when filtering is used.
    use_dense_mult = .NOT. PRESENT (filter_eps)
    IF (has_acc) use_dense_mult = .FALSE.
    !
    ! Temporary check for RMA algorithm
    IF (new_right.OR.new_left.OR.&
        dbcsr_mp_nprows(dbcsr_distribution_mp(matrix_c%m%dist)).NE.&
        dbcsr_mp_npcols(dbcsr_distribution_mp(matrix_c%m%dist))) THEN
       my_use_mpi_rma = .FALSE.
    ELSE
       my_use_mpi_rma = use_mpi_rma
    ENDIF
    !
    ! No dense multiplication and local indexing when 2D clustering is used.
    IF (my_use_mpi_rma) THEN
       use_dense_mult = .FALSE.
    ENDIF
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
       IF(first_row.LT.1 .OR. first_row.GT.dbcsr_nfullrows_total(matrix_c))&
          CPABORT("Invalid first row specified")
       f_row = first_row
    ENDIF
    IF (PRESENT (last_row)) THEN
       IF(last_row .GT. dbcsr_nfullrows_total(matrix_c))&
          CPABORT("Invalid last row specified")
       l_row = last_row
    ENDIF
    IF (PRESENT (first_column)) THEN
       IF(first_column.LT.1 .OR. first_column.GT.dbcsr_nfullcols_total(matrix_c))&
          CPABORT("Invalid first col specified")
       f_col = first_column
    ENDIF
    IF (PRESENT (last_column)) THEN
       IF(last_column .GT. dbcsr_nfullcols_total(matrix_c))&
          CPABORT("Invalid last column specified (C)")
       IF(last_column .GT. dbcsr_nfullcols_total(matrix_right))&
          CPABORT("Invalid last column specified (B)")
       l_col = last_column
    ENDIF
    IF (PRESENT (first_k)) THEN
       IF(first_k.LT.1 .OR. first_k.GT.dbcsr_nfullcols_total(matrix_left))&
          CPABORT("Invalid first k specified (A)")
       f_k = first_k
    ENDIF
    IF (PRESENT (last_k)) THEN
       IF(last_k .GT. dbcsr_nfullcols_total(matrix_left))&
          CPABORT("Invalid last k specified (A)")
       l_k = last_k
    ENDIF
    !
    ! update statistics (we count marketing flops per MPI rank)
    numnodes = dbcsr_mp_numnodes (dbcsr_distribution_mp (dbcsr_distribution (matrix_c)))
    last_mpi_ranks_used = numnodes
    marketing_flops = marketing_flops + &
         (2.0*(l_row-f_row+1.0)*(l_col-f_col+1.0)/numnodes)*(l_k-f_k+1.0)
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
    IF(f_row .GT. l_row)&
       CPABORT("Last row smaller than first row")
    IF(f_col .GT. l_col)&
       CPABORT("Last col smaller than first col")
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
            limits=(/f_row,l_row,f_col,l_col/))
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
    IF (.NOT. dbcsr_distribution_has_threads(matrix_c%m%dist)) THEN
       release_tdist = .TRUE.
       CALL dbcsr_distribution_make_threads (matrix_c%m%dist)
    ELSE
       release_tdist = .FALSE.
    ENDIF
    !
    ! Compute number of clusters
    IF (my_use_mpi_rma) THEN
       IF (dbcsr_distribution_has_row_clusters(matrix_left%m%dist)) THEN
          nimages_left_rows = dbcsr_distribution_num_local_row_clusters(matrix_left%m%dist)*&
               dbcsr_mp_nprows (dbcsr_distribution_mp (matrix_left%m%dist))
       ELSE 
          ! 2D clustering
          nimages_left_rows = dbcsr_distribution_get_num_images_1d(&
               dbcsr_nfullrows_total(matrix_left),&
               dbcsr_nblkrows_total(matrix_left),&
               dbcsr_mp_nprows (dbcsr_distribution_mp (matrix_left%m%dist)),&
               dbcsr_mp_npcols (dbcsr_distribution_mp (matrix_left%m%dist)))
       ENDIF
    ELSE
       ! 1D clustering for the old algorithm (only RMA support 2D clustering)
       nimages_left_rows = dbcsr_mp_nprows(dbcsr_distribution_mp(matrix_left%m%dist))       
    ENDIF
    !
    IF (my_use_mpi_rma.AND.dbcsr_distribution_has_col_clusters(matrix_left%m%dist)) THEN
       nimages_match = dbcsr_distribution_num_local_col_clusters(matrix_left%m%dist)*&
            dbcsr_mp_npcols (dbcsr_distribution_mp (matrix_left%m%dist))
    ELSE IF (my_use_mpi_rma.AND.dbcsr_distribution_has_row_clusters(matrix_right%m%dist)) THEN
       nimages_match = dbcsr_distribution_num_local_row_clusters(matrix_right%m%dist)*&
            dbcsr_mp_nprows (dbcsr_distribution_mp (matrix_right%m%dist))
    ELSE
       ! 2D clustering
       nimages_match = dbcsr_distribution_get_num_images_1d(&
            dbcsr_nfullcols_total(matrix_left),&
            dbcsr_nblkcols_total(matrix_left),&
            dbcsr_mp_nprows (dbcsr_distribution_mp (matrix_left%m%dist)),&
            dbcsr_mp_npcols (dbcsr_distribution_mp (matrix_left%m%dist)))
    ENDIF
    !
    IF (my_use_mpi_rma) THEN 
       IF (dbcsr_distribution_has_col_clusters(matrix_right%m%dist)) THEN
          nimages_right_cols = dbcsr_distribution_num_local_col_clusters(matrix_right%m%dist)*&
               dbcsr_mp_npcols (dbcsr_distribution_mp (matrix_right%m%dist))
       ELSE 
          nimages_right_cols = dbcsr_distribution_get_num_images_1d(&
               dbcsr_nfullcols_total(matrix_right),&
               dbcsr_nblkcols_total(matrix_right),&
               dbcsr_mp_nprows (dbcsr_distribution_mp (matrix_right%m%dist)),&
               dbcsr_mp_npcols (dbcsr_distribution_mp (matrix_right%m%dist)))
       ENDIF
    ELSE
       ! 1D clustering for the old algorithm (only RMA support 2D clustering)
       nimages_right_cols = dbcsr_mp_npcols(dbcsr_distribution_mp(matrix_right%m%dist))       
    ENDIF
    !
    ! Create imaged distributions for the multiply.
    CALL dbcsr_create_image_dist (rdist_right, matrix_right%m%dist,&
         match_row_nbins = dbcsr_mp_npcols (dbcsr_distribution_mp (matrix_left%m%dist)),&
         match_col_nbins = dbcsr_mp_npcols (dbcsr_distribution_mp (matrix_c%m%dist)),&
         match_col_pdist = dbcsr_distribution_col_dist (matrix_c%m%dist),&
         nimages_rows=nimages_match,&
         nimages_cols=nimages_right_cols)
    !
    CALL dbcsr_create_image_dist (rdist_left, matrix_left%m%dist,&
         match_row_pdist = dbcsr_distribution_row_dist (matrix_c%m%dist),&
         match_row_nbins = dbcsr_mp_nprows (dbcsr_distribution_mp (matrix_c%m%dist)),&
         match_col_pdist = dbcsr_distribution_row_dist (rdist_right%i%main),&
         match_col_idist = array_data (rdist_right%i%row_image),&
         match_col_nbins = dbcsr_mp_nprows (dbcsr_distribution_mp(matrix_right%m%dist)),&
         nimages_rows=nimages_left_rows,&
         nimages_cols=nimages_match)
    !
    IF (ab_dense) THEN
       CALL dbcsr_make_dists_dense (dbcsr_distribution (matrix_c),&
            rdist_left, rdist_right, dense_product_distribution,&
            dense_rdist_left, dense_rdist_right, .not.use_dense_mult,&
            m_map, k_vmap, n_map, matrix_c%m%row_blk_size)
       CALL make_sizes_dense (matrix_c%m%row_blk_size, m_map,&
            dbcsr_distribution_nrows (dense_product_distribution),&
            dense_row_sizes)
       CALL make_sizes_dense (matrix_c%m%col_blk_size, n_map, &
            dbcsr_distribution_ncols (dense_product_distribution),&
            dense_col_sizes)
       CALL make_sizes_dense (matrix_right%m%row_blk_size, k_vmap,&
            dbcsr_distribution_nrows (dense_rdist_right%i%main),&
            dense_k_sizes)
    ENDIF
    !
    IF(use_dense_mult .AND. .NOT.ab_dense)&
       CPABORT("Wrong logic when making dense matrices.")
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
            m_map, n_map)
    ELSE
       CALL dbcsr_init (product_matrix)
       CALL dbcsr_copy(product_matrix, matrix_c, shallow_data=.TRUE.)
    ENDIF
    IF (ab_dense) THEN
       CALL dbcsr_distribution_release (dense_product_distribution)
    ENDIF
    !
    ! This is needed to build the hash tables because they are
    ! locally indexed.
    CALL dbcsr_reset_locals (product_matrix)
    !
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
       IF(ASSOCIATED (product_matrix%m%wms))&
          CPABORT("Product matrix should be finalized!")
       CALL dbcsr_make_untransposed_blocks (product_matrix)
!$omp parallel &
!$omp default (none) shared (product_matrix)
       ! For the multiply logic to work correctly, existing data must
       ! be added only after the index has been transformed into the
       ! canonical form.
       CALL dbcsr_add_wm_from_matrix(product_matrix)
!$omp end parallel
    ELSE
!$omp parallel default(none) private(ithread) &
!$omp shared(product_matrix, memtype_product_wm)
       ithread = 0
       !$ ithread = OMP_GET_THREAD_NUM ()
       CALL dbcsr_work_create(product_matrix, work_mutable=.FALSE., &
               memory_type=memtype_product_wm(ithread)%p)
!$omp end parallel
    ENDIF
    product_matrix%m%nze = 0
    product_matrix%m%row_p(:) = 0
    CALL dbcsr_data_set_size_referenced(product_matrix%m%data_area, 0)
    product_matrix%m%nblks = 0
    product_matrix%m%valid = .FALSE.
    !
    IF (my_use_mpi_rma) THEN
       !
       max_nblocks = product_matrix%m%nblks
       !
       ! Right buffer images
       CALL dbcsr_make_buffers_rma(matrix_right,rdist_right,.FALSE.,&
                                   f_k,l_k,f_col,l_col,&
                                   alpha)
       !
       ! Left buffer images
       CALL dbcsr_make_buffers_rma(matrix_left,rdist_left,.TRUE.,&
                                   f_row,l_row,f_k,l_k)
    ELSE
       !
       ! Right images
       CALL make_m2s(matrix_right,m2s_right,rdist_right,dense_rdist_right,&
                     use_dense_mult,ab_dense,use_local_indexing,use_list_indexing,&
                     "R",f_k,l_k,f_row,l_row,f_col,l_col,&
                     dense_k_sizes,dense_col_sizes,&
                     k_vmap,m_map,n_map,&
                     alpha)
       !
       ! Left images
       CALL make_m2s(matrix_left,m2s_left,rdist_left,dense_rdist_left,&
                     use_dense_mult,ab_dense,use_local_indexing,use_list_indexing,&
                     "L",f_k,l_k,f_row,l_row,f_col,l_col,&
                     dense_row_sizes,dense_k_sizes,&
                     k_vmap,m_map,n_map)
    ENDIF
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
    IF (my_use_mpi_rma) THEN
       CALL cannon_multiply_low_rma(rdist_left,rdist_right,&
            matrix_left,matrix_right,product_matrix,&
            retain_sparsity=retain_sparsity,&
            filter_eps=filter_eps,&
            flop=my_flop)
    ELSE
       CALL cannon_multiply_low(m2s_left, m2s_right, product_matrix,&
            retain_sparsity=retain_sparsity,&
            filter_eps=filter_eps,&
            flop=my_flop)
    ENDIF

    CALL dbcsr_finalize(product_matrix)

    IF (PRESENT (flop)) THEN
       ! return the average number of flops per MPI rank. Variance (which is fairly large) could be computed as well.
       CALL timeset(routineN//"_mpsum_flop",handle2)
       comm = dbcsr_mp_group (dbcsr_distribution_mp (dbcsr_distribution (product_matrix)))
       numnodes = dbcsr_mp_numnodes (dbcsr_distribution_mp (dbcsr_distribution (product_matrix)))
       CALL mp_sum(my_flop,comm)
       flop = (my_flop + numnodes - 1) / numnodes
       CALL timestop(handle2)
    ENDIF
    !
    IF (new_left) CALL dbcsr_release (matrix_left)
    IF (new_right) CALL dbcsr_release (matrix_right)
    IF (release_tdist) THEN
       CALL dbcsr_distribution_no_threads (product_matrix%m%dist)
    ENDIF
    !
    CALL dbcsr_release_locals (product_matrix)
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
            m_map, n_map)
       CALL dbcsr_release (product_matrix)
       CALL array_release (old_product_row_blk_offsets)
       CALL array_release (old_product_col_blk_offsets)
       CALL array_release (old_product_row_blk_sizes)
       CALL array_release (old_product_col_blk_sizes)
       CALL dbcsr_distribution_release (old_product_distribution)
    ELSE
       CALL dbcsr_release (matrix_c)
       CALL dbcsr_init (matrix_c)
       CALL dbcsr_copy (matrix_c, product_matrix, shallow_data=.TRUE.)
       CALL dbcsr_release (product_matrix)
    ENDIF
    !
    IF (.NOT.my_use_mpi_rma) THEN
       CALL dbcsr_destroy_array (m2s_left)
       DEALLOCATE (m2s_left)
       CALL dbcsr_destroy_array (m2s_right)
       DEALLOCATE (m2s_right)
    ENDIF
    !
    CALL dbcsr_image_dist_release (rdist_left)
    CALL dbcsr_image_dist_release (rdist_right)
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
       CALL dbcsr_filter (matrix_c, eps_any, quick=.FALSE.)
    ENDIF
    !
    ! To support the canonical multiply (all non-transposed blocks),
    ! blocks may have to be transposed according to the CP2K
    ! triangular index.
    CALL dbcsr_make_untransposed_blocks (matrix_c)
    !
    IF (debug_mod .OR. careful_mod) THEN
       IF (debug_mod) &
            WRITE(*,*)routineN//" Use dense mult, symm",&
            use_dense_mult, ab_dense, dbcsr_has_symmetry (matrix_c)
       CALL dbcsr_verify_matrix (matrix_c)
       IF (debug_mod) THEN
          cs = dbcsr_checksum (matrix_c)
          WRITE(*,*)routineN//" Product checksum", cs
       ENDIF
    ENDIF

    ! This tends to trigger only when all of these conditions are fulfilled:
    !  - transa=="T"
    !  - matrix_c contains already blocks and beta is not zero
    !  - Cuda is enabled
    !  - multiple OpenMP threads are used
    IF(matrix_c%m%nblks > SIZE(array_data(matrix_c%m%row_blk_size)) * SIZE(array_data(matrix_c%m%col_blk_size)))&
       CPABORT("Bug: Matrix contains too many blocks")
    output_unit = default_output_unit
    CALL timestop(handle)
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
! *****************************************************************************
  SUBROUTINE make_m2s(matrix,m2s,rdist,dense_rdist,&
                      use_dense_mult,ab_dense,&
                      use_local_indexing,use_list_indexing,&
                      predistribute,f_k,l_k,f_row,l_row,f_col,l_col,&
                      row_blk_size,col_blk_size,&
                      k_vmap,m_map,n_map,&
                      alpha)
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
            no_copy_data=use_dense_mult, scale_value=alpha)
    ELSE
       CALL dbcsr_make_images (matrix, m2s, rdist,&
            predistribute=predistribute, &
            no_copy_data=use_dense_mult)
    END IF
    
    im = SIZE (m2s%mats,1)
    jm = SIZE (m2s%mats,2)
    SELECT CASE (predistribute)
    CASE('L')
       f_crop = (/ f_row, l_row, f_k, l_k /)
       row_map = m_map
       col_map = k_vmap
       thread_redist = .TRUE.
    CASE default
       f_crop = (/ f_k, l_k, f_col, l_col /)
       row_map = k_vmap
       col_map = n_map
       thread_redist = .FALSE.
    END SELECT

    ! Post-processing of images.
    DO i = 1, im
       DO j = 1, jm
          CALL dbcsr_reset_vlocals (m2s%mats(i,j), rdist)
          ! Crop if necessary
          IF (ANY(f_crop .NE. 0)) THEN
             CALL dbcsr_init (matrix_tmp)
             CALL dbcsr_crop_matrix (matrix_tmp, m2s%mats(i,j),&
                  full_row_bounds=f_crop(1:2),&
                  full_column_bounds=f_crop(3:4),&
                  shallow_data = .FALSE.)
             CALL dbcsr_release (m2s%mats(i,j))
             CALL dbcsr_copy (m2s%mats(i,j), matrix_tmp, shallow_data=.TRUE.)
             CALL dbcsr_release (matrix_tmp)
             CALL dbcsr_reset_vlocals (m2s%mats(i,j), rdist)
          ENDIF
       ENDDO
    ENDDO

    IF (ab_dense) THEN
       CALL dbcsr_init (dense_template)
       CALL dbcsr_create (dense_template, template=matrix,&
            dist=dense_rdist%i%main,&
            row_blk_size=array_data(row_blk_size), col_blk_size=array_data(col_blk_size))
       CALL dbcsr_make_images_dense (m2s, dense_rdist, &
            row_map = row_map, col_map = col_map,&
            join_cols = use_dense_mult, join_rows=ab_dense, &
            new_template=dense_template)
       CALL dbcsr_image_dist_release (rdist)
       rdist = dense_rdist
       CALL dbcsr_image_dist_hold (rdist)
       DO i = 1, im
          DO j = 1, jm
             CALL dbcsr_reset_vlocals (m2s%mats(i,j), rdist)
          ENDDO
       ENDDO
    ENDIF

    DO i = 1, im
       DO j = 1, jm
          ! Convert to local-row index
          IF (use_local_indexing) THEN
             CALL dbcsr_make_index_local_row(m2s%mats(i,j))
          ENDIF
          ! Convert to list index
          IF (use_list_indexing) THEN
             CALL dbcsr_make_index_list(m2s%mats(i,j), thread_redist=thread_redist)
          ENDIF
          IF (use_local_indexing .AND. .NOT. use_list_indexing) THEN
             CALL dbcsr_index_compact(m2s%mats(i,j))
          ENDIF
       ENDDO
    ENDDO

    IF (ab_dense) THEN
       CALL dbcsr_image_dist_release (dense_rdist)
       CALL dbcsr_release (dense_template)
    ENDIF

  END SUBROUTINE make_m2s

! *****************************************************************************
!> \brief Prepare orig images
!> \param matrix ...
!> \param imgdist ...
!> \param is_left ...
!> \param f_row ...
!> \param l_row ...
!> \param f_col ...
!> \param l_col ...
!> \param alpha ...
!> \author Alfio Lazzaro
! *****************************************************************************
  SUBROUTINE dbcsr_make_buffers_rma(matrix,imgdist,is_left,&
                                    f_row,l_row,f_col,l_col,&
                                    alpha)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(IN)                             :: imgdist
    LOGICAL, INTENT(IN)                      :: is_left
    INTEGER, INTENT(IN)                      :: f_row, l_row, f_col, l_col
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: alpha

    LOGICAL                                  :: do_scale

    do_scale = .FALSE.
    IF (PRESENT(alpha)) THEN
       IF (.NOT. dbcsr_scalar_are_equal (alpha, dbcsr_scalar_one(alpha%data_type))) THEN
          do_scale = .TRUE.
       END IF
    END IF
    !
    IF (do_scale) THEN
       CALL make_buffers_rma(matrix,imgdist,is_left,&
            f_row,l_row,f_col,l_col,alpha)
    ELSE
       CALL make_buffers_rma(matrix,imgdist,is_left,&
            f_row,l_row,f_col,l_col)
    ENDIF
  END SUBROUTINE dbcsr_make_buffers_rma

! *****************************************************************************
!> \brief Prepare orig images
!> \param matrix ...
!> \param imgdist ...
!> \param is_left ...
!> \param f_row ...
!> \param l_row ...
!> \param f_col ...
!> \param l_col ...
!> \param scale_value ...
!> \author Alfio Lazzaro
! *****************************************************************************
  SUBROUTINE make_buffers_rma(matrix,imgdist,is_left,&
                              f_row,l_row,f_col,l_col,&
                              scale_value)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(IN)                             :: imgdist
    LOGICAL, INTENT(IN)                      :: is_left
    INTEGER, INTENT(IN)                      :: f_row, l_row, f_col, l_col
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: scale_value

    CHARACTER(len=*), PARAMETER :: routineN = 'make_buffers_rma', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, blk_p, bp, col, col_img, col_size, coli, data_size, &
      data_type, dst_p, grp, handle, ithread, meta_size, mi, mynode, myt, &
      nblkcols_local, nblkrows_local, nblocks, ncols_images, nimages_merged, &
      nimages_unmerged, nrows_images, nsymmetries, nthreads, nze, pcol, prow, &
      row, row_img, row_size, rowi, size_threads_index, stored_col, &
      stored_row, sym_p, symmetry_i, tr_col_size, tr_row_size, ui
    INTEGER, ALLOCATABLE, DIMENSION(:), &
      TARGET                                 :: g2l_map_cols, g2l_map_rows, &
                                                meta_recv, meta_send, &
                                                tmp_img_offset
    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: recv_refs
    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: refs
    INTEGER, DIMENSION(2)                    :: recv_sym_size, sym_size
    INTEGER, DIMENSION(:), POINTER :: col_dist, col_img_dist, img_map_cols, &
      img_map_rows, local_cols, local_rows, meta_buffer, meta_buffer_p, &
      row_dist, row_img_dist, threads_dist
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi, img_refs_cols, &
                                                img_refs_rows, local_refs_p
    INTEGER, DIMENSION(:, :, :), POINTER     :: local_refs
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj), POINTER            :: data_buffer, data_buffer_p
    TYPE(dbcsr_data_obj), TARGET             :: data_recv, data_send
    TYPE(dbcsr_distribution_obj)             :: set_dist
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_scalar_type)                  :: scale_neg_one
    TYPE(dbcsr_type)                         :: sm

!$  INTEGER(kind=omp_lock_kind), ALLOCATABLE, DIMENSION(:) :: locks

    CALL timeset(routineN, handle)
    IF (.NOT.my_use_mpi_rma)&
         CPABORT("Only RMA algo supported here!")
    !
    ! Take input values and check validity
    IF (.NOT.dbcsr_valid_index (matrix))&
         CPABORT("Matrix not initialized.")
    sm = matrix%m
    data_type = sm%data_type
    IF(data_type .NE. dbcsr_type_real_8 .AND.&
       data_type .NE. dbcsr_type_real_4 .AND.&
       data_type .NE. dbcsr_type_complex_8 .AND.&
       data_type .NE. dbcsr_type_complex_4)&
       CPABORT("Invalid data type.")
    scale_neg_one = dbcsr_scalar_negative (dbcsr_scalar_one (data_type))
    set_dist = imgdist%i%main
    row_dist => dbcsr_distribution_row_dist (set_dist)
    col_dist => dbcsr_distribution_col_dist (set_dist)
    local_rows => dbcsr_distribution_local_rows(set_dist)
    local_cols => dbcsr_distribution_local_cols(set_dist)
    nblkrows_local = SIZE(local_rows)
    nblkcols_local = SIZE(local_cols)
    IF (sm%symmetry) THEN
       IF(SIZE(row_dist).NE.SIZE(col_dist))&
            CPWARN('Unequal row and column distributions for symmetric matrix.')
    ENDIF
    nrows_images = imgdist%i%row_decimation
    row_img_dist => array_data (imgdist%i%row_image)
    ncols_images = imgdist%i%col_decimation
    col_img_dist => array_data (imgdist%i%col_image)
    mp_obj = dbcsr_distribution_mp (set_dist)
    grp = dbcsr_mp_group (mp_obj)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    mynode = dbcsr_mp_mynode(mp_obj)
    IF (MAXVAL(row_dist).GT.UBOUND(blacs2mpi,1))&
         CPABORT("Row distribution references unexistent processor rows")
    IF (MAXVAL(col_dist).GT.UBOUND(blacs2mpi,2))&
         CPABORT("Col distribution references unexistent processor cols")
    IF (.NOT.dbcsr_mp_has_subgroups (mp_obj))&
         CPABORT("RMA requires subcommunicators!")
    ! Check threads configuration
    NULLIFY(threads_dist)
!$  IF (.NOT.dbcsr_distribution_has_threads (dbcsr_distribution (matrix)))&
!$       CPABORT("Thread distribution not defined")
!$  threads_dist => array_data(dbcsr_distribution_thread_dist (dbcsr_distribution(matrix)))
    !
    ! Crop matrix
    IF (ANY((/ f_row, l_row, f_col, l_col /).NE.0)) THEN
       IF (f_row.LT.0)&
            CPABORT("Invalid first row bound.")
       IF (l_row.GT.dbcsr_nfullrows_total(matrix))&
            CPABORT("Invalid last row bound.")
       IF (f_col.LT.0)&
            CPABORT("Invalid first column bound.")
       IF (l_col.GT.dbcsr_nfullcols_total(matrix))&
            CPABORT("Invalid last column bound.")
       !
       ! For the moment crop is not supported for RMA
       CPABORT("Crop not supported for RMA algo yet.")
       !
    ENDIF
    !
    IF (dbcsr_has_symmetry(matrix)) THEN
       nsymmetries = 2
       sym_p = blacs2mpi(dbcsr_mp_mypcol(mp_obj),dbcsr_mp_myprow(mp_obj))
    ELSE
       nsymmetries = 1
    ENDIF
    !
    size_threads_index = dbcsr_num_slots
    !
    IF (is_left) THEN
       nimages_merged = nrows_images
       nimages_unmerged = ncols_images
    ELSE
       nimages_merged = ncols_images
       nimages_unmerged = nrows_images
    ENDIF
    !
    ! Evaluate maps for global -> local indexing (g2l_map_rows, g2l_map_cols)
    ! Evaluate maps for local_image -> local indexing (img_map_rows, img_map_cols)
    ! Evaluate refs (count and offset) per each image
    ALLOCATE(g2l_map_rows(sm%nblkrows_total),g2l_map_cols(sm%nblkcols_total))
    ALLOCATE(img_refs_rows(2,nrows_images),img_refs_cols(2,ncols_images))
    ALLOCATE(img_map_rows(nblkrows_local),img_map_cols(nblkcols_local))
    !
    img_refs_rows(:,:) = 0
    g2l_map_rows(:) = 0
    DO row = 1, nblkrows_local
       row_img = row_img_dist (local_rows(row))
       rowi = MOD(row_img-1,nrows_images)+1
       img_refs_rows(1,rowi) = img_refs_rows(1,rowi)+1
       g2l_map_rows(local_rows(row)) = img_refs_rows(1,rowi)
    ENDDO
    DO row_img = 1, nrows_images-1
       img_refs_rows(2,row_img+1) = img_refs_rows(2,row_img) + img_refs_rows(1,row_img)
    ENDDO
    ALLOCATE(tmp_img_offset(nrows_images))
    tmp_img_offset(:) = img_refs_rows(2,:)
    DO row = 1, nblkrows_local
       row_img = row_img_dist (local_rows(row))
       rowi = MOD(row_img-1,nrows_images)+1
       tmp_img_offset(rowi) = tmp_img_offset(rowi)+1
       img_map_rows(tmp_img_offset(rowi)) = row
    ENDDO
    DEALLOCATE(tmp_img_offset)
    !
    img_refs_cols(:,:) = 0
    g2l_map_cols(:) = 0
    DO col = 1, nblkcols_local
       col_img = col_img_dist (local_cols(col))
       coli = MOD(col_img-1,ncols_images)+1
       img_refs_cols(1,coli) = img_refs_cols(1,coli)+1
       g2l_map_cols(local_cols(col)) = img_refs_cols(1,coli)
    ENDDO
    DO col_img = 1, ncols_images-1
       img_refs_cols(2,col_img+1) = img_refs_cols(2,col_img) + img_refs_cols(1,col_img)
    ENDDO
    ALLOCATE(tmp_img_offset(ncols_images))
    tmp_img_offset(:) = img_refs_cols(2,:)
    DO col = 1, nblkcols_local
       col_img = col_img_dist (local_cols(col))
       coli = MOD(col_img-1,ncols_images)+1
       tmp_img_offset(coli) = tmp_img_offset(coli)+1
       img_map_cols(tmp_img_offset(coli)) = col
    ENDDO
    DEALLOCATE(tmp_img_offset)
    !
    ALLOCATE(local_refs(4,nrows_images,ncols_images))
    local_refs(:,:,:) = 0
    data_size = 0
    meta_size = 0
    nblocks = 0
    !
!$omp parallel default (none) &
!$omp private (ithread,myt,iter,row,col,blk,row_size,col_size,&
!$omp          stored_row,stored_col,mi,ui,blk_p,bp,tr,&
!$omp          nze,symmetry_i,row_img,col_img,rowi,coli,&
!$omp          tr_row_size,tr_col_size,prow,pcol,dst_p,&
!$omp          data_buffer_p,meta_buffer_p) &
!$omp shared (nthreads,refs,matrix,nsymmetries,&
!$omp         row_img_dist,col_img_dist,imgdist,row_dist,col_dist,&
!$omp         is_left,nrows_images,ncols_images,&
!$omp         local_refs,nimages_merged,nimages_unmerged,local_refs_p,&
!$omp         data_size,meta_size,data_type,rma_buffers_orig,&
!$omp         memtype_mpi_buffer,data_buffer,meta_buffer,sm,&
!$omp         img_refs_cols,img_refs_rows,img_map_cols,img_map_rows,&
!$omp         g2l_map_cols,g2l_map_rows,recv_refs,grp,&
!$omp         scale_value,scale_neg_one,data_send,data_recv,&
!$omp         size_threads_index,sym_size,recv_sym_size,meta_send,&
!$omp         mp_obj,threads_dist,meta_recv,&
!$omp         locks,mynode,blacs2mpi,sym_p) &
!$omp         reduction(+ : nblocks)
    ithread = 0
!$  ithread = omp_get_thread_num()
    myt = ithread
!$omp master
    nthreads = 1
!$  nthreads = omp_get_num_threads()    
    ALLOCATE(refs(2,nrows_images,ncols_images,nsymmetries,0:nthreads))
    refs(:,:,:,:,1:nthreads) = 0
!$  IF (is_left) THEN
!$     size_threads_index = size_threads_index + nthreads + 1
!$  ENDIF
!$  IF (is_left.AND.nsymmetries.GT.1) THEN
!$     ALLOCATE(locks(0:nthreads-1))
!$  ENDIF
!$omp end master
!$omp barrier
!$  IF (is_left.AND.nsymmetries.GT.1) THEN
!$     call omp_init_lock(locks(ithread))
!$  ENDIF
    !
    ! By default all data are considered already in local proc
    dst_p = ilocal_proc
    !
    ! Take data and meta dimensions per each thread
    CALL dbcsr_iterator_start(iter, matrix, shared=.TRUE.)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk, &
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
          row_img = row_img_dist (stored_row)
          col_img = col_img_dist (stored_col)
          CALL image_calculator(imgdist,&
               prow = prow, pcol = pcol,&
               rowi = rowi, coli = coli,&
               myprow = row_dist(stored_row), myrowi = row_img,&
               mypcol = col_dist(stored_col), mycoli = col_img,&
               shifting = '0')
          IF (nsymmetries.GT.1) THEN
             IF (blacs2mpi(prow, pcol).EQ.mynode) THEN
                dst_p = ilocal_proc
             ELSE
                ! Move data to remove symmetric proc
                dst_p = isym_proc
             ENDIF
          ENDIF
!$        IF (is_left.AND.nsymmetries.GT.1) THEN
!$           myt = threads_dist(stored_row)
!$        ENDIF
!$omp atomic
          refs(imeta_size,rowi,coli,dst_p,myt+1) = refs(imeta_size,rowi,coli,dst_p,myt+1) + 3
!$omp atomic
          refs(idata_size,rowi,coli,dst_p,myt+1) = refs(idata_size,rowi,coli,dst_p,myt+1) + nze
       ENDDO ! symmetry_i
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    IF (nsymmetries.GT.1) THEN
!$omp barrier
!$omp master
       ALLOCATE(recv_refs(2,nrows_images,ncols_images,0:nthreads))
       ! exchange sizes with the symmetric proc
       CALL mp_sendrecv(refs(:,:,:,isym_proc,1:nthreads),sym_p,&
            recv_refs(:,:,:,1:nthreads),sym_p,grp)
!$omp end master
    ENDIF
!$omp barrier
!$omp critical(local_refs_acc)
    !
    ! Store local data and meta dimensions
    local_refs(imeta_size,:,:) = local_refs(imeta_size,:,:) + refs(imeta_size,:,:,ilocal_proc,ithread+1)
    local_refs(idata_size,:,:) = local_refs(idata_size,:,:) + refs(idata_size,:,:,ilocal_proc,ithread+1)
    IF (nsymmetries.GT.1) THEN
       local_refs(imeta_size,:,:) = local_refs(imeta_size,:,:) + recv_refs(imeta_size,:,:,ithread+1)
       local_refs(idata_size,:,:) = local_refs(idata_size,:,:) + recv_refs(idata_size,:,:,ithread+1)
    ENDIF
!$omp end critical(local_refs_acc)
!$omp barrier
!$omp master
    !
    ! Take the displacements, taking in account the images merging
    ! Note that local diplacements refer to position 0
    DO ui = 1, nimages_unmerged
       IF (is_left) THEN
          ! merging over rows
          local_refs_p => local_refs(:,:,ui)
       ELSE
          ! merging over cols
          local_refs_p => local_refs(:,ui,:)
       ENDIF
       DO mi = 1, nimages_merged
          local_refs_p(imeta_displ,mi) = meta_size
          local_refs_p(idata_displ,mi) = data_size
          ! Include stats slots
          IF (local_refs_p(imeta_size,mi).GT.0) THEN
             local_refs_p(imeta_size,mi) = local_refs_p(imeta_size,mi) + size_threads_index
          ENDIF
          meta_size = meta_size + local_refs_p(imeta_size,mi)
          data_size = data_size + local_refs_p(idata_size,mi)
       ENDDO
    ENDDO
    !
    ! Change displacements for the first thread
    ! In this case they refer to position 1
    refs(imeta_size,:,:,ilocal_proc,0) = local_refs(imeta_displ,:,:)+1
    refs(idata_size,:,:,ilocal_proc,0) = local_refs(idata_displ,:,:)+1
    !
    ! Evaluate data and meta diplacements per each thread
    DO mi=1, nthreads-1
       refs(:,:,:,ilocal_proc,mi) = refs(:,:,:,ilocal_proc,mi) + refs(:,:,:,ilocal_proc,mi-1)
    ENDDO
    !
    ! Take offsets for symmetric case
    IF (nsymmetries.GT.1) THEN
       !
       ! For the symmetric case, always goes in rows and columns order w.r.t the the images
       sym_size(:) = 1
       recv_sym_size(:) = 1
       DO coli=1, ncols_images
          DO rowi=1, nrows_images
             ! Carry previous value
             refs(:,rowi,coli,isym_proc,0) = sym_size(:)
             recv_refs(:,rowi,coli,0) = recv_sym_size(:)
             DO mi=1, nthreads
                refs(:,rowi,coli,isym_proc,mi) = refs(:,rowi,coli,isym_proc,mi) + refs(:,rowi,coli,isym_proc,mi-1)
                recv_refs(:,rowi,coli,mi) = recv_refs(:,rowi,coli,mi) + recv_refs(:,rowi,coli,mi-1)
                ! Add symetric events to the diplacements of local data
                refs(:,rowi,coli,ilocal_proc,mi) = refs(:,rowi,coli,ilocal_proc,mi) + &
                     recv_refs(:,rowi,coli,mi) - recv_sym_size(:)
             ENDDO
             sym_size(:) = refs(:,rowi,coli,isym_proc,nthreads)
             recv_sym_size(:) = recv_refs(:,rowi,coli,nthreads)
             ! Use the slot of the last thread to store the displ
             refs(:,rowi,coli,isym_proc,nthreads) = refs(:,rowi,coli,isym_proc,0)-1
          ENDDO
       ENDDO
    ENDIF
    !
    ! Allocate data and meta buffers
    IF (is_left) THEN
       CALL rma_buffers_init(rma_buffers_orig,data_type,&
                             left_data_size=data_size,&
                             left_meta_size=meta_size,&
                             left_local_refs=local_refs,&
                             left_img_map_rows=img_map_rows,&
                             left_img_map_cols=img_map_cols,&
                             left_img_refs_rows=img_refs_rows,&
                             left_img_refs_cols=img_refs_cols,&
                             data_memory_type=memtype_mpi_buffer)
       data_buffer => rma_buffers_orig%left_data
       meta_buffer => rma_buffers_orig%left_meta
    ELSE
       CALL rma_buffers_init(rma_buffers_orig,data_type,&
                             right_data_size=data_size,&
                             right_meta_size=meta_size,&
                             right_local_refs=local_refs,&
                             right_img_map_rows=img_map_rows,&
                             right_img_map_cols=img_map_cols,&
                             right_img_refs_rows=img_refs_rows,&
                             right_img_refs_cols=img_refs_cols,&
                             data_memory_type=memtype_mpi_buffer)
       data_buffer => rma_buffers_orig%right_data
       meta_buffer => rma_buffers_orig%right_meta
    ENDIF
    !
    IF (nsymmetries.GT.1) THEN
       ! Allocate data/meta to send
       CALL dbcsr_data_init (data_send)
       CALL dbcsr_data_new (data_send, data_type, sym_size(idata_size)-1)
       ALLOCATE(meta_send(sym_size(imeta_size)-1))
    ENDIF
!$omp end master
!$omp barrier
    !
    data_buffer_p => data_buffer
    meta_buffer_p => meta_buffer
    !
    ! Copy data and meta in the buffers 
    CALL dbcsr_iterator_start(iter, matrix, shared=.TRUE.)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk, blk_p=blk_p, &
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
          row_img = row_img_dist (stored_row)
          col_img = col_img_dist (stored_col)
          CALL image_calculator(imgdist,&
               prow = prow, pcol = pcol,&
               rowi = rowi, coli = coli,&
               myprow = row_dist(stored_row), myrowi = row_img,&
               mypcol = col_dist(stored_col), mycoli = col_img,&
               shifting = '0')
          IF (nsymmetries.GT.1) THEN
             IF (blacs2mpi(prow, pcol).EQ.mynode) THEN
                dst_p = ilocal_proc
                data_buffer_p => data_buffer
                meta_buffer_p => meta_buffer
             ELSE
                ! Move data to remove symmetric proc
                dst_p = isym_proc
                data_buffer_p => data_send
                meta_buffer_p => meta_send
             ENDIF
          ENDIF
!$        IF (is_left.AND.nsymmetries.GT.1) THEN
!$           myt = threads_dist(stored_row)
!$           call omp_set_lock(locks(myt))
!$        ENDIF
          IF (tr) THEN
             CALL dbcsr_block_transpose_aa(data_buffer_p,sm%data_area,tr_row_size,tr_col_size,&
                                           refs(idata_size,rowi,coli,dst_p,myt),bp,&
                                           scale_value)
             IF (sm%negate_real .AND. sm%negate_imaginary) THEN
                CALL dbcsr_block_scale (data_buffer_p, scale=scale_neg_one,&
                                        row_size=nze, col_size=1, &
                                        lb=refs(idata_size,rowi,coli,dst_p,myt))
             ELSEIF (sm%negate_real) THEN
                CALL dbcsr_block_real_neg (data_buffer_p, row_size=nze, col_size=1,&
                                           lb=refs(idata_size,rowi,coli,dst_p,myt))
             ELSEIF (sm%negate_imaginary) THEN
                CALL dbcsr_block_conjg (data_buffer_p, row_size=nze, col_size=1,&
                                        lb=refs(idata_size,rowi,coli,dst_p,myt))
             ENDIF
          ELSE
             CALL dbcsr_block_copy_aa(data_buffer_p,sm%data_area,row_size,col_size,&
                                      refs(idata_size,rowi,coli,dst_p,myt),bp,&
                                      scale_value)
          ENDIF
          !
          ! Set meta data (local indexing)
          IF (nsymmetries.EQ.1.OR.blacs2mpi(prow, pcol).EQ.mynode) THEN
             meta_buffer_p(refs(imeta_size,rowi,coli,dst_p,myt)+size_threads_index) = g2l_map_rows(stored_row)
             meta_buffer_p(refs(imeta_size,rowi,coli,dst_p,myt)+size_threads_index+1) = g2l_map_cols(stored_col)
             meta_buffer_p(refs(imeta_size,rowi,coli,dst_p,myt)+size_threads_index+2) = &
                  refs(idata_size,rowi,coli,dst_p,myt) - local_refs(idata_displ,rowi,coli)
          ELSE
             ! For the symmetric case we need to swap rows and cols maps
             meta_buffer_p(refs(imeta_size,rowi,coli,dst_p,myt)) = g2l_map_cols(stored_row)
             meta_buffer_p(refs(imeta_size,rowi,coli,dst_p,myt)+1) = g2l_map_rows(stored_col)
             ! The offset will be adjusted on the symmetric proc
             meta_buffer_p(refs(imeta_size,rowi,coli,dst_p,myt)+2) = &
                  refs(idata_size,rowi,coli,dst_p,myt)
          ENDIF
          refs(imeta_size,rowi,coli,dst_p,myt) = refs(imeta_size,rowi,coli,dst_p,myt) + 3
          refs(idata_size,rowi,coli,dst_p,myt) = refs(idata_size,rowi,coli,dst_p,myt) + nze
!$        IF (is_left.AND.nsymmetries.GT.1) THEN
!$           call omp_unset_lock(locks(myt))
!$        ENDIF
       ENDDO
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    !
    IF (nsymmetries.GT.1) THEN
!$omp barrier
!$omp master
       ! Send/recv data to the symmetric proc
       CALL dbcsr_data_init (data_recv)
       CALL dbcsr_data_new (data_recv, data_type, recv_sym_size(idata_size)-1)
       CALL dbcsr_sendrecv_any(data_send,sym_p,data_recv,sym_p,grp)       
       ALLOCATE(meta_recv(recv_sym_size(imeta_size)-1))
       CALL mp_sendrecv(meta_send,sym_p,meta_recv,sym_p,grp)       
       ! Delete send data
       CALL dbcsr_data_release (data_send)
       DEALLOCATE(meta_send)
!$omp end master
    ENDIF
!$omp barrier
!$  IF (is_left.AND.nsymmetries.GT.1) THEN
!$     call omp_destroy_lock(locks(ithread))    
!$  ENDIF
    !
    ! Set stats slots
!$omp do schedule(dynamic) collapse(2)
    DO coli = 1, ncols_images
       DO rowi = 1, nrows_images
          IF (local_refs(imeta_size,rowi,coli).EQ.0) CYCLE
          !
          ! Insert symmetric data in local buffers
          IF (nsymmetries.GT.1) THEN
             IF (recv_refs(imeta_size,rowi,coli,nthreads).GT.1) THEN
                DO ithread=0, nthreads-1
                   ! Copy meta, block by block
                   DO blk = recv_refs(imeta_size,rowi,coli,ithread),recv_refs(imeta_size,rowi,coli,ithread+1)-1,3
                      meta_buffer(refs(imeta_size,rowi,coli,ilocal_proc,ithread)+size_threads_index) = &
                           meta_recv(blk)
                      meta_buffer(refs(imeta_size,rowi,coli,ilocal_proc,ithread)+size_threads_index+1) = &
                           meta_recv(blk+1)
                      meta_buffer(refs(imeta_size,rowi,coli,ilocal_proc,ithread)+size_threads_index+2) = &
                           meta_recv(blk+2)-recv_refs(idata_size,rowi,coli,ithread)+&
                           refs(idata_size,rowi,coli,ilocal_proc,ithread)-local_refs(idata_displ,rowi,coli)
                      refs(imeta_size,rowi,coli,ilocal_proc,ithread) = refs(imeta_size,rowi,coli,ilocal_proc,ithread) + 3
                   ENDDO
                   ! Copy data
                   CALL dbcsr_data_set(data_buffer,refs(idata_size,rowi,coli,ilocal_proc,ithread),&
                        recv_refs(idata_size,rowi,coli,ithread+1)-recv_refs(idata_size,rowi,coli,ithread),&
                        data_recv,recv_refs(idata_size,rowi,coli,ithread))
                ENDDO
             ENDIF
          ENDIF
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_size) = &
               local_refs(imeta_size,rowi,coli)
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_nblks) = &
               (local_refs(imeta_size,rowi,coli)-size_threads_index)/3
          nblocks = nblocks+meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_nblks)
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_nze) = &
               local_refs(idata_size,rowi,coli)
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_dense) = 0
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_nblkrows_total) = &
               sm%nblkrows_total
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_nblkcols_total) = &
               sm%nblkcols_total
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_nfullrows_total) = &
               sm%nfullrows_total
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_nfullcols_total) = &
               sm%nfullcols_total
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_nblkrows_local) = &
               img_refs_rows(1,rowi)
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_nblkcols_local) = &
               img_refs_cols(1,coli)
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_nfullrows_local) = &
               sm%nfullrows_local
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_nfullcols_local) = &
               sm%nfullcols_local
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_type) = 0
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_meta_size+1) = 0
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_row_p) = 1
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_row_p+1) = 0
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_col_i) = 1
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_col_i+1) = 0
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_blk_p) = 1
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_blk_p+1) = 0
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_thr_c) = 0
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_thr_c+1) = 0
          IF (is_left) THEN
             meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_thr_c) = dbcsr_num_slots+1
             meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_thr_c+1) = size_threads_index
             ! first thread
             meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_num_slots+1) = 0
             DO ithread=0, nthreads-1
                meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_num_slots+ithread+2) = &
                     (refs(imeta_size,rowi,coli,ilocal_proc,ithread)-local_refs(imeta_displ,rowi,coli)-1)/3
             ENDDO
          ENDIF
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_coo_l) = size_threads_index+1
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_num_slots) = &
               local_refs(imeta_size,rowi,coli)
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_home_prow) = &
               dbcsr_mp_myprow(mp_obj)
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_home_pcol) = &
               dbcsr_mp_mypcol(mp_obj)
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_home_rowi) = rowi
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_home_coli) = coli
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_home_vprow) = &
               dbcsr_mp_myprow(mp_obj)*nrows_images+rowi-1
          meta_buffer(local_refs(imeta_displ,rowi,coli)+dbcsr_slot_home_vpcol) = &
               dbcsr_mp_mypcol(mp_obj)*ncols_images+coli-1
       ENDDO
    ENDDO
!$omp end do
!$omp end parallel
    DEALLOCATE(refs)
    IF (nsymmetries.GT.1) THEN
       CALL dbcsr_data_release (data_recv)
       DEALLOCATE(recv_refs,meta_recv)
    ENDIF
    DEALLOCATE(g2l_map_rows,g2l_map_cols)
!$  IF (is_left.AND.nsymmetries.GT.1) THEN    
!$     DEALLOCATE(locks)
!$  ENDIF
    !
    max_nblocks = MAX(max_nblocks,nblocks)
    !
    CALL timestop(handle)
  END SUBROUTINE make_buffers_rma

! *****************************************************************************
!> \brief Creates row and column images of a matrix.
!> \param[in] source          input matrix
!> \param[in,out] normalized  image array of the normalized matrix
!> \param[in] target_image_dist          normalize to this image distribution
!> \param[in] predistribute   (optional) predistribute data for multiplication
!> \param[in] no_copy_data    (optional) try to not merge data at the end
!> \param[in] scale_value     (optional) scale with this value
! *****************************************************************************
  SUBROUTINE dbcsr_make_images(source, normalized, target_image_dist,&
       predistribute, no_copy_data, scale_value)
    TYPE(dbcsr_obj), INTENT(IN)              :: source
    TYPE(dbcsr_2d_array_type), INTENT(OUT)   :: normalized
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(IN)                             :: target_image_dist
    CHARACTER, INTENT(IN), OPTIONAL          :: predistribute
    LOGICAL, INTENT(IN), OPTIONAL            :: no_copy_data
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: scale_value

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_images', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, handle)
    IF(my_use_mpi_rma)&
         CPABORT("RMA algo not supported here!")
    IF(.NOT.dbcsr_valid_index (source))&
         CPABORT("Matrix not initialized.")
    CALL make_images(source, normalized,&
         target_image_dist, desymmetrize=dbcsr_has_symmetry(source),&
         predistribute=predistribute,&
         no_copy_data=no_copy_data,&
         scale_value=scale_value)
    normalized%image_dist = target_image_dist
    CALL dbcsr_image_dist_hold (normalized%image_dist)
    CALL timestop(handle)
  END SUBROUTINE dbcsr_make_images

! *****************************************************************************
!> \brief Makes column-based and row-based images of a matrix.
!> \param[in] ism input symmetric matrix
!> \param[in,out] ums         normalied matrices
!> \param[in] target_imgdist  image distribution to normalize to
!> \param[in] desymmetrize    (optional) desymmetrize a symmetric matrix
!> \param[in] predistribute   (optional) predistribute data for multiplication
!> \param[in] no_copy_data    (optional) try to not merge data at the end
!> \param[in] scale_value     (optional) scale with this value
! *****************************************************************************
  SUBROUTINE make_images(ism, ums, target_imgdist, desymmetrize, predistribute,&
       no_copy_data, scale_value)
    TYPE(dbcsr_obj), INTENT(IN)              :: ism
    TYPE(dbcsr_2d_array_type), INTENT(OUT)   :: ums
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(IN)                             :: target_imgdist
    LOGICAL, INTENT(IN), OPTIONAL            :: desymmetrize
    CHARACTER, INTENT(IN), OPTIONAL          :: predistribute
    LOGICAL, INTENT(IN), OPTIONAL            :: no_copy_data
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: scale_value

    CHARACTER(len=*), PARAMETER :: routineN = 'make_images', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: metalen = 5
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    CHARACTER                                :: predist_type, predist_type_fwd
    INTEGER :: blk, blk_l, blk_p, bp, col, col_img, col_size, coli, data_p, &
      data_type, dst_p, handle, handle2, ithread, mp_group, ncol_images, &
      nrow_images, nsymmetries, nthreads, numproc, nze, pcol, prev_blk_p, &
      prow, row, row_img, row_size, rowi, sd_pos, sm_pos, src_p, &
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
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_memtype_type)                 :: data_memory_type
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_scalar_type)                  :: scale_neg_one
    TYPE(dbcsr_type)                         :: sm

!   ---------------------------------------------------------------------------
! Check input matrix
! Set convenient variables to access input matrix info
!

    CALL timeset(routineN, handle)
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
       CPABORT("Incorrect pre-shift specifier.")
    END SELECT
    data_type = sm%data_type
    IF(data_type .NE. dbcsr_type_real_8 .AND.&
       data_type .NE. dbcsr_type_real_4 .AND.&
       data_type .NE. dbcsr_type_complex_8 .AND.&
       data_type .NE. dbcsr_type_complex_4)&
       CPABORT("Invalid data type.")
    scale_neg_one = dbcsr_scalar_negative (dbcsr_scalar_one (data_type))
    row_blk_size => array_data (sm%row_blk_size)
    col_blk_size => array_data (sm%col_blk_size)
    old_dist = dbcsr_distribution (ism)
    target_dist = target_imgdist%i%main
    row_dist => dbcsr_distribution_row_dist (target_dist)
    col_dist => dbcsr_distribution_col_dist (target_dist)
    IF (sm%symmetry) THEN
       IF(SIZE(row_dist).NE.SIZE(col_dist))&
          CPWARN('Unequal row and column distributions for symmetric matrix.')
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
    IF(MAXVAL(row_dist).GT.UBOUND(blacs2mpi,1))&
       CPABORT('Row distribution references unexistent processor rows')
    IF (dbg) THEN
         IF(MAXVAL(row_dist).NE.UBOUND(blacs2mpi,1))&
            CPWARN('Range of row distribution not equal to processor rows')
    ENDIF
    IF(MAXVAL(col_dist).GT.UBOUND(blacs2mpi,2))&
       CPABORT('Col distribution references unexistent processor cols')
    IF (dbg) THEN
         IF(MAXVAL(col_dist).NE.UBOUND(blacs2mpi,2))&
            CPWARN('Range of col distribution not equal to processor cols')
    ENDIF

    ! Check threads configuration
!$  IF(.NOT.dbcsr_distribution_has_threads (old_dist))&
!$     CPABORT("Thread distribution not defined")

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
    !
    ! Allocate and init mats
    !
    ALLOCATE (ums%mats(nrow_images,ncol_images))
    data_memory_type = memtype_abpanel_1
    DO row_img = 1, nrow_images
       DO col_img = 1, ncol_images
          CALL dbcsr_init (ums%mats(row_img, col_img))
          CALL dbcsr_create(ums%mats(row_img, col_img), "imaged "//sm%name,&
               target_dist,&
               dbcsr_type_no_symmetry, array_data(sm%row_blk_size), array_data(sm%col_blk_size),&
               sm%row_blk_size, sm%col_blk_size,&
               nze=0, data_type=data_type,&
               max_rbs=sm%max_rbs,max_cbs=sm%max_cbs,&
               row_blk_offset=sm%row_blk_offset,col_blk_offset=sm%col_blk_offset,&
               thread_dist=sm%dist,&
               data_memory_type=data_memory_type,&
               index_memory_type=memtype_mpi_buffer,&
               make_index=.FALSE.)
          ums%mats(row_img, col_img)%m%negate_real = sm%negate_real
          ums%mats(row_img, col_img)%m%negate_imaginary = sm%negate_imaginary
       ENDDO
    ENDDO
    
    nthreads = 1
!$omp parallel default (none) &
!$omp private (ithread, symmetry_i, row_img, col_img, &
!$omp          myt_send_count, myt_total_send_count, &
!$omp          iter, row, col, blk, row_size, col_size, stored_row, stored_col, &
!$omp          prow, pcol, rowi, coli, vrow, vcol, dst_p, nze, myt_smp, myt_sdp, &
!$omp          blk_p, bp, sd_pos, sm_pos,tr, &
!$omp          tr_row_size, tr_col_size) &
!$omp shared (nthreads, nsymmetries, row_img_dist, col_img_dist, &
!$omp         nrow_images, ncol_images, numproc, scale_value, &
!$omp         ums, sm, ism, target_imgdist, row_dist, col_dist,&
!$omp         predist_type_fwd, blacs2mpi, row_blk_size, col_blk_size, &
!$omp         send_count, recv_count, handle2,mp_group, &
!$omp         total_send_count, total_recv_count, recv_data_area, nocopy, &
!$omp         data_type, recv_meta, send_data_area, send_meta, &
!$omp         sd_disp, sm_disp, rd_disp, rm_disp, all_total_send_offset, blk_ps, blks, &
!$omp         received_data_area, scale_neg_one, memtype_abpanel_1)
    ithread = 0
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
               shifting = predist_type_fwd)
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
    CALL timeset(routineN//"_sizes", handle2)
    CALL mp_alltoall(send_count, recv_count, 2*nrow_images*ncol_images,&
         mp_group)
    CALL timestop(handle2)
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
       CALL dbcsr_data_new (recv_data_area, data_type, SUM(total_recv_count(2, :)), memory_type=memtype_abpanel_1)
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
       DO row_img = 1, nrow_images
          DO col_img = 1, ncol_images
             CALL dbcsr_work_create(ums%mats(row_img, col_img),&
                  SUM(recv_count(1,row_img,col_img,:)), n=1)
             CALL dbcsr_data_hold (received_data_area)
             CALL dbcsr_data_release (ums%mats(row_img,col_img)%m%wms(1)%data_area)
             ums%mats(row_img,col_img)%m%wms(1)%data_area = received_data_area
          ENDDO
       ENDDO
    ENDIF
!$omp barrier
!$omp master
    CALL m_memory()
!$omp end master

    ! Add timing call to the packing of the send buffers 
    !
    CALL timeset(routineN//"_pack", handle2)
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
               shifting = predist_type_fwd)
          dst_p = blacs2mpi(prow, pcol)
          sm_pos = myt_smp(dst_p)
          myt_smp(dst_p) = myt_smp(dst_p) + metalen
          sd_pos = myt_sdp(dst_p)
          myt_sdp(dst_p) = myt_sdp(dst_p) + nze
          IF (tr) THEN
             CALL dbcsr_block_transpose_aa(send_data_area,sm%data_area,tr_row_size,tr_col_size,&
                                           sd_pos,bp,scale_value)
             IF (sm%negate_real .AND. sm%negate_imaginary) THEN
                CALL dbcsr_block_scale (send_data_area, scale=scale_neg_one,&
                                        row_size=nze, col_size=1, lb=sd_pos)
             ELSEIF (sm%negate_real) THEN
                CALL dbcsr_block_real_neg (send_data_area,&
                                           row_size=nze, col_size=1, lb=sd_pos)
             ELSEIF (sm%negate_imaginary) THEN
                CALL dbcsr_block_conjg (send_data_area,&
                                        row_size=nze, col_size=1, lb=sd_pos)
             ENDIF
          ELSE
             CALL dbcsr_block_copy_aa(send_data_area,sm%data_area,row_size,col_size,&
                                      sd_pos,bp,scale_value)
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

    CALL timestop(handle2)
!$omp end parallel

    ! Exchange the data and metadata structures. In the interesting cases (square grids, row col distribution same),
    ! there are only very few processors that need to exchange data.
    ! The hybrid_alltoall deals with this by doing point to point communication
    CALL timeset(routineN//"_data", handle2)
    CALL hybrid_alltoall_any(send_data_area,total_send_count(2,:),sd_disp(:)-1,&
                             recv_data_area,total_recv_count(2,:),rd_disp(:)-1,&
                             mp_obj,&
                             most_ptp=.TRUE.,remainder_ptp=.TRUE.,no_hybrid=.FALSE.)
    CALL hybrid_alltoall_i1(&
         send_meta(:), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), metalen*total_recv_count(1,:), rm_disp(:)-1,&
         most_ptp=.TRUE., remainder_ptp=.TRUE., no_hybrid=.FALSE.,&
         mp_env = mp_obj)
    CALL timestop(handle2)

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
          !
          CALL dbcsr_data_set_size_referenced (&
               ums%mats(row_img,col_img)%m%wms(1)%data_area,&
               ums%mats(row_img,col_img)%m%wms(1)%datasize)
          IF (nrow_images.EQ.1 .AND. ncol_images.eq.1 .OR. nocopy) THEN
             CALL dbcsr_special_finalize(ums%mats(row_img,col_img),reshuffle=.FALSE.)
          ELSE
             CALL dbcsr_special_finalize(ums%mats(row_img,col_img),reshuffle=.TRUE.)
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
               shifting=predist_type)
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

    CALL timestop(handle)
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
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE dbcsr_make_images_dense (images, new_rdist, &
       row_map, col_map, join_cols, join_rows, new_template)
    TYPE(dbcsr_2d_array_type), INTENT(INOUT) :: images
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(INOUT)                          :: new_rdist
    TYPE(array_i1d_obj), INTENT(IN)          :: row_map, col_map
    LOGICAL, INTENT(IN)                      :: join_cols, join_rows
    TYPE(dbcsr_obj), INTENT(IN)              :: new_template

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_images_dense', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: handle, mat_col, mat_row, &
                                                mat_vpcol, mat_vprow
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

    CALL timeset(routineN, handle)
    old_rdist = images%image_dist
    !
    DO mat_row = 1, images%image_dist%i%row_decimation
       DO mat_col = 1, images%image_dist%i%col_decimation
          IF (dbg) THEN
             cs = dbcsr_checksum (images%mats(mat_row,mat_col))
             WRITE(*,*)routineN//" cs pre", cs
          ENDIF
          mat_vprow = images%mats(mat_row, mat_col)%m%index(dbcsr_slot_home_vprow)
          mat_vpcol = images%mats(mat_row, mat_col)%m%index(dbcsr_slot_home_vpcol)
          und_row_blk_offsets => array_data (images%mats(mat_row, mat_col)%m%row_blk_offset)
          und_col_blk_offsets => array_data (images%mats(mat_row, mat_col)%m%col_blk_offset)
          CALL dbcsr_get_local_vrows (old_rdist, und_local_vrows, mat_vprow)
          CALL dbcsr_get_local_vcols (old_rdist, und_local_vcols, mat_vpcol)

          CALL dbcsr_get_local_vrows (new_rdist, dense_local_vrows, mat_vprow)
          CALL dbcsr_get_local_vcols (new_rdist, dense_local_vcols, mat_vpcol)
          ! The old matrix has to be remembered so it is copied to
          ! tmp_mat.
          old_meta(:) = images%mats(mat_row,mat_col)%m%index(1:dbcsr_meta_size)
          CALL dbcsr_init (tmp_mat)
          tmp_mat = images%mats(mat_row,mat_col)
          CALL dbcsr_init (images%mats(mat_row,mat_col))
          CALL dbcsr_create (images%mats(mat_row,mat_col), template=new_template)
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
               array_data(row_map), array_data(col_map), join_rows, join_cols)
          !
          CALL dbcsr_index_prune_deleted (images%mats(mat_row, mat_col))
          !
          CALL dbcsr_release (tmp_mat)
          IF (dbg) THEN
             cs = dbcsr_checksum (images%mats(mat_row,mat_col))
             WRITE(*,*)routineN//" cs pst", cs
          ENDIF
       ENDDO
    ENDDO
    CALL dbcsr_image_dist_release (images%image_dist)
    images%image_dist = new_rdist
    CALL dbcsr_image_dist_hold (images%image_dist)
    CALL timestop(handle)
  END SUBROUTINE dbcsr_make_images_dense

! *****************************************************************************
!> \brief Multiplies two DBCSR matrices
!>
!> \param[in] left_set             set of imaged left matrices
!> \param[in] right_set            set of imaged right matrices
!> \param[out] product_matrix      DBCSR product matrix 
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix; default is no
!> \param filter_eps ...
!> \param[out] flop                (optional) effective flop
! *****************************************************************************
  SUBROUTINE cannon_multiply_low(left_set, right_set, product_matrix,&
       retain_sparsity, &
       filter_eps, flop)
    TYPE(dbcsr_2d_array_type), POINTER       :: left_set, right_set
    TYPE(dbcsr_obj), INTENT(INOUT)           :: product_matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(kind=real_8), INTENT(in), OPTIONAL  :: filter_eps
    INTEGER(KIND=int_8), INTENT(OUT)         :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'cannon_multiply_low', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: idata = 1, ileft = 0, &
                                                imeta = 2, iright = 2, &
                                                M_L = 2, M_P = 1, M_R = 3

    INTEGER :: data_type, data_type_byte, grp, handle, handle2, handle3, i, &
      ithread, left_col_image, left_col_mult, left_col_nimages, &
      left_dst_icol, left_dst_irow, left_dst_p, left_dst_pcol, left_dst_prow, &
      left_dst_vcol, left_dst_vrow, left_max_nblks, left_max_nze, &
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
      right_dst_vcol, right_dst_vrow, right_max_nblks, right_max_nze, &
      right_myfirstvcol, right_myfirstvrow, right_mypcol, right_myprow, &
      right_npcols, right_nprows, right_recv_icol, right_recv_irow
    INTEGER :: right_recv_p, right_recv_pcol, right_recv_prow, &
      right_recv_vcol, right_recv_vrow, right_row_image, right_row_mult, &
      right_row_nimages, right_send_icol, right_send_irow, right_send_p, &
      right_send_pcol, right_send_prow, right_send_vcol, right_send_vrow, &
      right_src_icol, right_src_irow, right_src_p, right_src_pcol, &
      right_src_prow, right_src_vcol, right_src_vrow, row, size_guess, stat, &
      threads_finished, threads_finished_read, v_k, v_ki
    INTEGER(KIND=int_8)                      :: flop_single, flop_total, mem
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
    TYPE(dbcsr_mm_multrec_type_p), DIMENSION(:), ALLOCATABLE :: multrec
    TYPE(dbcsr_mp_obj)                       :: left_mp_obj, product_mp_obj, &
                                                right_mp_obj

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, handle)
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
!$omp shared (multrec, nthreads, product_matrix)
!$omp master
    nthreads = 1
    !$  nthreads = OMP_GET_NUM_THREADS ()
    IF(.NOT.ASSOCIATED (product_matrix%m%wms))&
       CPABORT("Work matrices do not exist")
    IF(SIZE (product_matrix%m%wms).NE.nthreads)&
       CPABORT("Work matrices not correctly sized.")
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
    IF(left_col_nimages.NE.right_row_mult)&
       CPABORT("Left/Right image mismatch")
    IF(left_col_mult.NE.right_row_nimages)&
       CPABORT("Left/Right image mismatch")
    IF(left_col_nimages * left_npcols .NE. right_row_nimages * right_nprows)&
       CPABORT("Left/Right total mismatch")
    ! product/left matching
    IF(left_row_mult * dbcsr_mp_nprows (product_mp_obj) .NE. left_row_nimages * left_nprows)&
       CPABORT("Product/Left total mismatch")
    ! product/left matching
    IF(right_col_mult * dbcsr_mp_npcols (product_mp_obj) .NE. right_col_nimages * right_npcols)&
       CPABORT("Product/Right total mismatch")
    ! Limitations
    IF(left_row_nimages .NE. 1)&
       CPABORT("Product/Left matrix process grid mismatch")
    IF(left_row_mult .NE. 1)&
       CPABORT("Product/Left matrix process grid mismatch")
    IF(right_col_nimages .NE. 1)&
       CPABORT("Product/Right matrix process grid mismatch")
    IF(right_col_mult .NE. 1)&
       CPABORT("Product/Right matrix process grid mismatch")

    dbcsr_mpi_statistics%nimages = MAX(dbcsr_mpi_statistics%nimages,left_row_nimages*left_col_nimages)
    dbcsr_mpi_statistics%nimages = MAX(dbcsr_mpi_statistics%nimages,right_row_nimages*right_col_nimages)
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
       IF(stat .NE. 0)&
          CPABORT("Could not allocate memory")
    ELSE
       IF (careful_mod) THEN
          IF(left_set%mats(1, 1)%m%bcsc)&
             CPABORT("Can not do on-the-fly filtering with CSC-indexed matrices.")
       ENDIF
       IF (dbcsr_has_local_row_index (left_set%mats(1, 1))) THEN
          nblkrows_used = dbcsr_nblkrows_local (left_set%mats(1, 1))
       ELSE
          nblkrows_used = dbcsr_nblkrows_total (left_set%mats(1, 1))
       ENDIF
       nblkrows_total = dbcsr_nblkrows_total (left_set%mats(1, 1))
       ALLOCATE (row_max_epss (nblkrows_total), stat=stat)
       IF(stat .NE. 0)&
          CPABORT("Could not allocate memory for left epsilons")
       ALLOCATE (row_counts (nblkrows_used), stat=stat)
       IF(stat .NE. 0)&
          CPABORT("Could not allocate memory for left row counts")
       ! The summation could be done prow-locally but it would
       ! complicate the pre-row eps calculation.
       ALLOCATE (total_row_counts (nblkrows_total), stat=stat)
       IF(stat .NE. 0)&
          CPABORT("Could not allocate memory for left row counts")
       ! Each prow member matrix (npcols * row_images) counts the
       ! blocks present in each of its rows.
       total_row_counts(:) = 0
       DO left_row_image = 1, left_row_nimages
          DO left_col_image = 1, left_col_nimages
             list_indexing =&
                  left_set%mats(left_row_image, left_col_image)%m%list_indexing
             IF (careful_mod) THEN
                IF (list_indexing) THEN
                   IF((left_set%mats(left_row_image, left_col_image)%m%nblks)*3 .NE. &
                      SIZE(left_set%mats(left_row_image, left_col_image)%m%coo_l))&
                      CPABORT("Row count mismatch")
                ELSE
                   IF(nblkrows_used+1 .NE. SIZE(left_set%mats(left_row_image, left_col_image)%m%row_p))&
                      CPABORT("Row count mismatch")
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
                IF(SIZE(local_rows) .NE. SIZE(row_counts))&
                   CPABORT("Mismatch in number of local rows.")
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
       DEALLOCATE (row_counts)
       DEALLOCATE (total_row_counts)
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
               shifting = 'l')
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
               shifting = 'r')
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
!$omp          shared (nthreads, keep_sparsity, mynode)
    !
    !$ ithread = OMP_GET_THREAD_NUM()
    ! The work arrays have to be setup (actually, not quite sure).
    i = ithread + 1
    size_guess = product_matrix%m%wms(i)%datasize ! Should be minimal
    CALL dbcsr_data_ensure_size(product_matrix%m%wms(i)%data_area,&
         size_guess)
    CALL dbcsr_data_set_size_referenced (product_matrix%m%wms(i)%data_area,&
         product_matrix%m%wms(i)%datasize)
    ! XXXXXXX a quick fix right now, allocation with size 1 might actually not be needed at all,
    !         but something expects this to be associated
    CALL ensure_array_size(product_matrix%m%wms(i)%row_i, ub=1)
    CALL ensure_array_size(product_matrix%m%wms(i)%col_i, ub=1)
    CALL ensure_array_size(product_matrix%m%wms(i)%blk_p, ub=1)
!$omp end parallel

    ! update capacity of memory-pools
    IF (has_acc) THEN
      CALL dbcsr_mempool_ensure_capacity(memtype_abpanel_1%pool, &
      capacity=left_row_mult*left_col_nimages + right_row_nimages*right_col_mult)
      CALL dbcsr_mempool_ensure_capacity(memtype_abpanel_2%pool, &
      capacity=left_row_mult*left_col_nimages + right_row_nimages*right_col_mult)
    ENDIF

    !
    ! Setup the left buffer matrices
    !
    CALL buffer_matrices_ensure_size(left_set, index_size=left_max_nblks,&
         data_size=left_max_nze)

    CALL setup_buffer_matrices (left_buffer_2, left_row_mult, left_col_nimages,&
         left_set%mats(1,1), index_size=left_max_nblks,&
         data_size=left_max_nze)
    IF (otf_filtering) THEN
       ALLOCATE (left_norms (left_max_nblks), stat=stat)
       IF(stat .NE. 0)&
          CPABORT("Could not allocate memory for left norms")
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
         data_size=right_max_nze)
    
    CALL setup_buffer_matrices (right_buffer_2, right_row_nimages, right_col_mult,&
         right_set%mats(1,1), index_size=right_max_nblks, data_size=right_max_nze)
    IF (otf_filtering) THEN
       ALLOCATE (right_norms (right_max_nblks), stat=stat)
       IF(stat.NE.0)&
          CPWARN("Could not allocate memory for right norms")
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
!$omp         keep_sparsity, filter_eps, row_max_epss, multrec, nthreads, &
!$omp         right_data_sr, right_data_rr, left_data_sr, left_data_rr,&
!$omp         right_index_sr, right_index_rr, left_index_sr, left_index_rr), &
!$omp private(ithread)
    ithread = 0
    !$ ithread = OMP_GET_THREAD_NUM ()
    ALLOCATE(multrec(ithread)%p)
    CALL dbcsr_mm_multrec_init(multrec(ithread)%p,&
         left_buffer_comm%mats(1, 1)%m,&
         right_buffer_comm%mats(1, 1)%m,&
         product_matrix%m,&
         keep_sparsity=keep_sparsity,&
         eps=filter_eps,&
         row_max_epss=row_max_epss,&
         block_estimate=MAX(product_matrix%m%nblks,&
                            left_buffer_comm%mats(1, 1)%m%nblks,&
                            right_buffer_comm%mats(1, 1)%m%nblks)/nthreads,&
         right_row_blk_size=array_data(right_buffer_comm%mats(1, 1)%m%row_blk_size))
!$omp end parallel
    !
    ! Setup indexing
    CALL setup_rec_index (left_set, left_row_nimages, left_col_nimages)
    CALL setup_rec_index (right_set, right_row_nimages, right_col_nimages)
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
    IF (has_acc) THEN
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
       CALL timeset(routineN//"_metrocomm1", handle2)
       wait_right: IF (metronome .GT. 1) THEN
          IF (debug_mod) WRITE (*,'(1X,A)')routineN//" waiting for right"
          !
          CALL mp_waitall (right_data_sr)
          CALL mp_waitall (right_data_rr)
          CALL mp_waitall (right_index_sr)
          CALL mp_waitall (right_index_rr)
       ENDIF wait_right
       CALL timestop(handle2)
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
                  shifting='0')
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
                  shifting='R')
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
                  shifting='0')
             ! Calculate which data I receive
             CALL image_calculator (right_set%image_dist,&
                  prow=right_src_prow, rowi=right_src_irow,&
                  pcol=right_src_pcol, coli=right_src_icol,&
                  vprow=right_src_vrow, vpcol=right_src_vcol,&
                  myvprow=v_ki+right_myfirstvrow,&
                  myvpcol=right_myfirstvcol,&
                  ! receive window moves with the metronome
                  vprow_shift=metronome*min_nimages,&
                  shifting='R')
             !
             IF (has_acc) THEN
                CALL timeset(routineN//"_acc_sync_right", handle3)
                CALL acc_event_synchronize(right_buffer_comm%mats(v_ki+1, 1)%m%data_area%d%acc_ready)
                CALL timestop(handle3)
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
             CALL timeset(routineN//"_metrocomm2", handle2)
             CALL dbcsr_irecv_any (right_data_rp, right_recv_p,&
                  grp, right_data_rr(v_ki+1), tag=right_src_vrow)
             CALL mp_irecv (right_index_rp, right_recv_p,&
                  grp, right_index_rr(v_ki+1), tag=right_src_vrow)
             CALL dbcsr_isend_any (right_data_sp, right_send_p,&
                  grp, right_data_sr(v_ki+1), tag=right_dst_vrow)
             CALL mp_isend (right_index_sp, right_send_p,&
                  grp, right_index_sr(v_ki+1), tag=right_dst_vrow)
             dbcsr_mpi_statistics%nexchanged = dbcsr_mpi_statistics%nexchanged+1
             CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(1,:),&
                  dbcsr_data_get_size(right_data_rp),&
                  dbcsr_mpi_statistics%data_size_breakdown(:,:,1),data_type_byte)
             CALL timestop(handle2)
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
       CALL timeset(routineN//"_metrocomm3", handle2)
       wait_left: IF (metronome .GT. 1) THEN
          IF (debug_mod) WRITE (*,'(1X,A)')routineN//" waiting for left"
          CALL mp_waitall (left_data_sr)
          CALL mp_waitall (left_data_rr)
          CALL mp_waitall (left_index_sr)
          CALL mp_waitall (left_index_rr)
       ENDIF wait_left
       CALL timestop(handle2)
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
                  shifting='0')
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
                  shifting='L')
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
                  shifting='0')
             ! Calculate which data I receive
             CALL image_calculator (left_set%image_dist,&
                  prow=left_src_prow, rowi=left_src_irow,&
                  pcol=left_src_pcol, coli=left_src_icol,&
                  vprow=left_src_vrow, vpcol=left_src_vcol,&
                  myvprow=left_myfirstvrow,&
                  myvpcol=v_ki+left_myfirstvcol,&
                  ! receive window moves with the metronome
                  vpcol_shift=metronome*min_nimages,&
                  shifting='L')
             !
             IF (has_acc) THEN
                CALL timeset(routineN//"_acc_sync_left", handle3)
                CALL acc_event_synchronize(left_buffer_comm%mats(1, v_ki+1)%m%data_area%d%acc_ready)
                CALL timestop(handle3)
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
             CALL timeset(routineN//"_metrocomm4", handle2)
             CALL dbcsr_irecv_any (left_data_rp, left_recv_p,&
                  grp, left_data_rr(v_ki+1), tag=left_src_vcol)
             CALL mp_irecv (left_index_rp, left_recv_p,&
                  grp, left_index_rr(v_ki+1), tag=left_src_vcol)
             CALL dbcsr_isend_any (left_data_sp, left_send_p,&
                  grp, left_data_sr(v_ki+1), tag=left_dst_vcol)
             CALL mp_isend (left_index_sp, left_send_p,&
                  grp, left_index_sr(v_ki+1), tag=left_dst_vcol)
             dbcsr_mpi_statistics%nexchanged = dbcsr_mpi_statistics%nexchanged+1
             CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(2,:),&
                  dbcsr_data_get_size(left_data_rp),&
                  dbcsr_mpi_statistics%data_size_breakdown(:,:,2),data_type_byte)
             CALL timestop(handle2)
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

       CALL timeset(routineN//"_multrec", handle2)
       DO v_ki = 0, min_nimages-1
          IF (debug_mod) THEN
             CALL dbcsr_print(left_buffer_calc%mats(1, v_ki+1), nodata=.TRUE.)
             CALL dbcsr_print(right_buffer_calc%mats(v_ki+1, 1), nodata=.TRUE.)
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
                  right_set%image_dist)
          ENDIF
          IF (.TRUE. .OR. left_buffer_calc%mats(1, v_ki+1)%m%local_indexing) THEN
             ! Since the right matrix is shifted vertically, the
             ! received data always has different notions of "local
             ! rows".  Thus the local_rows and global_rows must be
             ! recalculated.
             CALL dbcsr_reset_vlocals (left_buffer_calc%mats(1, v_ki+1),&
                  left_set%image_dist)
          ENDIF

          IF (has_acc) THEN
             CALL dbcsr_data_host2dev(left_buffer_calc%mats(1, v_ki+1)%m%data_area)
             CALL dbcsr_data_host2dev(right_buffer_calc%mats(v_ki+1, 1)%m%data_area)
             CALL acc_transpose_blocks(right_buffer_calc%mats(v_ki+1, 1), trs_stackbuf_calc)
          END IF
          
          ! Sets the local right-matrix columns
          IF (otf_filtering) THEN
             left_norms(:) = huge_norm
             right_norms(:) = huge_norm
             CALL calculate_norms(right_buffer_calc%mats(v_ki+1, 1),&
                  right_norms)
             CALL calculate_norms(left_buffer_calc%mats(1, v_ki+1),&
                  left_norms)
          ENDIF
          !
          flop_single = 0
          threads_finished = 0

!$omp parallel default (none) &
!$omp shared (left_buffer_calc, right_buffer_calc, &
!$omp         v_ki, &
!$omp         product_matrix, multrec,&
!$omp         filter_eps, right_norms, left_norms, row_max_epss, &
!$omp         keep_sparsity,threads_finished, &
!$omp         right_data_sr, right_data_rr, right_index_sr, right_index_rr, &
!$omp         left_data_sr, left_data_rr, left_index_sr, left_index_rr, &
!$omp         use_comm_thread) &
!$omp private (ithread,nthreads,threads_finished_read) &
!$omp firstprivate (metronome, nsteps_k, min_nimages) &
!$omp reduction (+: flop_single)
          ithread = 0; nthreads = 1
!$        ithread = omp_get_thread_num(); nthreads = omp_get_num_threads()

          IF(metronome==nsteps_k .AND. v_ki==min_nimages-1) &
               CALL dbcsr_mm_multrec_phaseout(multrec(ithread)%p)

          CALL dbcsr_mm_multrec_multiply(multrec(ithread)%p,&
               left=left_buffer_calc%mats(1, v_ki+1)%m,&
               right=right_buffer_calc%mats(v_ki+1, 1)%m,&
               flop=flop_single,&
               a_norms=left_norms, b_norms=right_norms)

          IF(metronome==nsteps_k .AND. v_ki==min_nimages-1) THEN
             CALL dbcsr_mm_multrec_finalize(multrec(ithread)%p)
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
!$omp end parallel
          flop_total = flop_total + flop_single
       ENDDO

       CALL timestop(handle2)

       CALL dbcsr_switch (left_buffer_calc, left_buffer_comm)
       CALL dbcsr_switch (right_buffer_calc, right_buffer_comm)
       CALL dbcsr_switch (trs_stackbuf_calc, trs_stackbuf_comm)

    ENDDO grouped_k_index
    CALL m_memory(mem)
    max_memory = MAX(max_memory,REAL(mem))

    IF (has_acc) THEN
       CALL dbcsr_data_release(trs_stackbuf_1)
       CALL dbcsr_data_release(trs_stackbuf_2)
    END IF

    IF (ALLOCATED (right_norms)) THEN
       DEALLOCATE (right_norms)
    ENDIF
    IF (ALLOCATED (left_norms)) THEN
       DEALLOCATE (left_norms)
    ENDIF
    IF (ALLOCATED (row_max_epss)) THEN
       DEALLOCATE (row_max_epss)
    ENDIF
    !
    CALL dbcsr_destroy_array (right_buffer_2)
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
    DEALLOCATE(left_data_rr, left_data_sr, left_index_rr, left_index_sr, &
               right_data_rr, right_data_sr, right_index_rr, right_index_sr)
    !
    !
    IF (debug_mod) THEN
       v_ki = 0
       DO i = 1, SIZE(product_matrix%m%blk_p)
          v_ki = MAX(v_ki, ABS(product_matrix%m%blk_p(i)))
       ENDDO
       WRITE(*,*)routineN//" Actual final size",&
            LOG(REAL(dbcsr_data_get_size(product_matrix%m%data_area)))/LOG(10.0),&
            LOG(REAL(v_ki))/LOG(10.0)
    ENDIF
    !
    flop = flop_total
    DEALLOCATE (left_buffer_2, right_buffer_2)
    !
    CALL timestop(handle)
  END SUBROUTINE cannon_multiply_low


! *****************************************************************************
!> \brief Multiplies two DBCSR matrices by means of RMA MPI operations
!>
!> \param imgdist_left ...
!> \param imgdist_right ...
!> \param matrix_left ...
!> \param matrix_right ...
!> \param[out] product_matrix      DBCSR product matrix 
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix; default is no
!> \param filter_eps ...
!> \param[out] flop                (optional) effective flop
!> \author Alfio Lazzaro
! *****************************************************************************
  SUBROUTINE cannon_multiply_low_rma(imgdist_left, imgdist_right,&
                                     matrix_left, matrix_right,&
                                     product_matrix,&
                                     retain_sparsity,&
                                     filter_eps, flop)
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(INOUT)                          :: imgdist_left, imgdist_right
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_left, matrix_right
    TYPE(dbcsr_obj), INTENT(INOUT)           :: product_matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    REAL(kind=real_8), INTENT(IN), OPTIONAL  :: filter_eps
    INTEGER(KIND=int_8), INTENT(OUT)         :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'cannon_multiply_low_rma', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, data_type, data_type_byte, dst_icol, dst_irow, dst_p, &
      grp_left, grp_right, handle, handle2, handle3, handle4, &
      index_row_max_epss, ithread, left_col_mult, left_col_nimages, &
      left_data_size, left_data_win, left_max_data_size, left_max_meta_size, &
      left_meta_offset, left_meta_size, left_meta_win, left_myfirstvcol, &
      left_myfirstvrow, left_mypcol, left_myprow, left_npcols, left_nprows, &
      left_nrequests, left_row_mult, left_row_nimages, left_vcol, left_vrow, &
      meta_size_image, metronome, min_nimages, mynode, nblkrows_local, nblks, &
      nsteps_k, nthreads, numnodes, nvirt_k, recv_coli, recv_pcol
    INTEGER :: recv_prow, recv_rowi, recv_vcol, recv_vrow, right_col_mult, &
      right_col_nimages, right_data_size, right_data_win, &
      right_max_data_size, right_max_meta_size, right_meta_offset, &
      right_meta_size, right_meta_win, right_myfirstvcol, right_myfirstvrow, &
      right_mypcol, right_myprow, right_npcols, right_nprows, &
      right_nrequests, right_row_mult, right_row_nimages, right_vcol, &
      right_vrow, row, size_guess, size_threads_index, src_vcol, src_vrow, &
      threads_finished, threads_finished_read, v_ci, v_i, v_ki, v_ri
    INTEGER(KIND=int_8)                      :: mem
    INTEGER, ALLOCATABLE, DIMENSION(:) :: left_max_meta_size_merged, &
      right_max_meta_size_merged, row_counts, total_row_counts
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: left_get_requests, &
                                                left_recv_filtered, &
                                                right_get_requests, &
                                                right_recv_filtered
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: left_refs, right_refs
    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: left_all_refs, right_all_refs
    INTEGER, DIMENSION(:), POINTER           :: meta_get, rowi
    INTEGER, DIMENSION(&
      dbcsr_slot_nblkrows_total:&
      dbcsr_slot_nfullcols_local)            :: left_first_indices, &
                                                right_first_indices
    LOGICAL                                  :: keep_sparsity, otf_filtering
    REAL(kind=sp)                            :: min_row_max_epss
    REAL(kind=sp), ALLOCATABLE, DIMENSION(:) :: row_max_epss
    REAL(kind=sp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: left_max_norms, left_norms, &
                                                right_max_norms, right_norms
    REAL(kind=sp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: left_all_max_norms, &
                                                right_all_max_norms
    TYPE(dbcsr_data_obj)                     :: data_get
    TYPE(dbcsr_mm_multrec_type_p), &
      ALLOCATABLE, DIMENSION(:)              :: multrec
    TYPE(dbcsr_mp_obj)                       :: left_mp_obj, product_mp_obj, &
                                                right_mp_obj
    TYPE(dbcsr_rma_buffers), POINTER         :: rma_buffers, &
                                                rma_buffers_calc, &
                                                rma_buffers_comm

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, handle)
    NULLIFY(rowi)
    NULLIFY(meta_get)
    NULLIFY(rma_buffers_calc, rma_buffers_comm, rma_buffers)
    !
    IF (PRESENT (retain_sparsity)) THEN
       keep_sparsity = retain_sparsity
    ELSE
       keep_sparsity = .FALSE.
    ENDIF
    otf_filtering = PRESENT (filter_eps)
    !
!$omp parallel default (none) &
!$omp shared (nthreads)
!$omp master
    nthreads = 1
!$  nthreads = OMP_GET_NUM_THREADS ()
!$omp end master
!$omp end parallel
    !
    size_threads_index = 0
!$  size_threads_index = nthreads+1
    !
    ! Dummy checks
    IF(.NOT.ASSOCIATED (product_matrix%m%wms))&
       CPABORT("Work matrices do not exist")
    IF(SIZE (product_matrix%m%wms) .NE. nthreads)&
       CPABORT("Work matrices not correctly sized.")
    IF(.NOT.rma_buffers_orig%is_valid_left.OR.&
       .NOT.rma_buffers_orig%is_valid_right.OR.&
       .NOT.ASSOCIATED(rma_buffers_orig%left_meta).OR.&
       .NOT.ASSOCIATED(rma_buffers_orig%right_meta).OR.&
       .NOT.ASSOCIATED(rma_buffers_orig%left_local_refs).OR.&
       .NOT.ASSOCIATED(rma_buffers_orig%right_local_refs).OR.&
       .NOT.ASSOCIATED(rma_buffers_orig%left_img_refs_rows).OR.&
       .NOT.ASSOCIATED(rma_buffers_orig%left_img_refs_cols).OR.&
       .NOT.ASSOCIATED(rma_buffers_orig%right_img_refs_rows).OR.&
       .NOT.ASSOCIATED(rma_buffers_orig%right_img_refs_cols))&
       CPABORT("No buffers associated for the RMA algo!")
    !
    ! Set up variables
    flop = 0
    data_type = dbcsr_get_data_type (product_matrix)
    data_type_byte = dbcsr_datatype_sizeof(data_type)
    left_row_nimages =  imgdist_left%i%row_decimation
    left_row_mult =     imgdist_left%i%row_multiplicity
    left_col_nimages =  imgdist_left%i%col_decimation
    left_col_mult =     imgdist_left%i%col_multiplicity
    right_row_nimages = imgdist_right%i%row_decimation
    right_row_mult =    imgdist_right%i%row_multiplicity
    right_col_nimages = imgdist_right%i%col_decimation
    right_col_mult =    imgdist_right%i%col_multiplicity
    left_mp_obj    = dbcsr_distribution_mp (imgdist_left%i%main)
    right_mp_obj   = dbcsr_distribution_mp (imgdist_right%i%main)
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
    grp_right = dbcsr_mp_my_col_group (right_mp_obj)
    grp_left = dbcsr_mp_my_row_group (left_mp_obj)
    CALL dbcsr_mp_grid_setup (product_mp_obj)
    CALL dbcsr_mp_grid_setup (left_mp_obj)
    CALL dbcsr_mp_grid_setup (right_mp_obj)
    !
    ! Dummy checks
    IF(SIZE(rma_buffers_orig%left_local_refs,2).NE.left_row_nimages.OR.&
       SIZE(rma_buffers_orig%left_local_refs,3).NE.left_col_nimages.OR.&
       SIZE(rma_buffers_orig%right_local_refs,2).NE.right_row_nimages.OR.&
       SIZE(rma_buffers_orig%right_local_refs,3).NE.right_col_nimages)&
       CPABORT("Number of images mismatch for the RMA algo!")
    IF(SIZE(rma_buffers_orig%left_img_refs_rows,2).NE.left_row_nimages.OR.&
       SIZE(rma_buffers_orig%left_img_refs_cols,2).NE.left_col_nimages.OR.&
       SIZE(rma_buffers_orig%right_img_refs_rows,2).NE.right_row_nimages.OR.&
       SIZE(rma_buffers_orig%right_img_refs_cols,2).NE.right_col_nimages)&
       CPABORT("Number of images mismatch for the RMA algo!")
    ! subcommunicators
    IF(.NOT.dbcsr_mp_has_subgroups (right_mp_obj))&
         CPABORT("RMA requires rows subcommunicators for right matrix!")
    IF(.NOT.dbcsr_mp_has_subgroups (left_mp_obj))&
         CPABORT("RMA requires columns subcommunicators for left matrix!")
    ! left/right matching
    IF(left_col_nimages .NE. right_row_mult)&
         CPABORT("Left/Right image mismatch")
    IF(left_col_mult .NE. right_row_nimages)&
         CPABORT("Left/Right image mismatch")
    IF(left_col_nimages * left_npcols .NE. right_row_nimages * right_nprows)&
         CPABORT("Left/Right total mismatch")
    ! product/left matching
    IF(left_row_mult * dbcsr_mp_nprows (product_mp_obj) .NE. left_row_nimages * left_nprows)&
         CPABORT("Product/Left total mismatch")
    ! product/left matching
    IF(right_col_mult * dbcsr_mp_npcols (product_mp_obj) .NE. right_col_nimages * right_npcols)&
         CPABORT("Product/Right total mismatch")
    !
    dbcsr_mpi_statistics%nimages = MAX(dbcsr_mpi_statistics%nimages,&
                                       left_row_nimages*left_col_nimages)
    dbcsr_mpi_statistics%nimages = MAX(dbcsr_mpi_statistics%nimages,&
                                       right_row_nimages*right_col_nimages)
    !
    ! The main transfer loop goes through the virtual rows/columns.
    ! The number of steps may be smaller if the grid dimension is very
    ! non-optimal (both left column images and right row images are >
    ! 1).
    min_nimages = MIN (left_col_nimages, right_row_nimages)
    nvirt_k = left_npcols * left_col_nimages
    nsteps_k = nvirt_k / min_nimages
    !
    ! Exchange size data
    ALLOCATE (left_all_refs(4,left_row_nimages,left_col_nimages,0:left_npcols-1))
    CALL mp_allgather(rma_buffers_orig%left_local_refs,left_all_refs,grp_left)
    ALLOCATE (left_refs(4,left_myfirstvrow:left_myfirstvrow+left_row_nimages-1,0:nvirt_k-1))
    !
    ALLOCATE (right_all_refs(4,right_row_nimages,right_col_nimages,0:right_nprows-1))
    CALL mp_allgather(rma_buffers_orig%right_local_refs,right_all_refs,grp_right)
    ALLOCATE (right_refs(4,right_myfirstvcol:right_myfirstvcol+right_col_nimages-1,0:nvirt_k-1))
    !
    ! Evaluate and exchange max norms for MPI filtering
    IF (otf_filtering.AND.use_mpi_filtering) THEN
       ALLOCATE(left_max_norms(left_row_nimages,left_col_nimages))
       CALL calculate_max_image_norms(rma_buffers_orig%left_meta,rma_buffers_orig%left_data,&
                                      rma_buffers_orig%left_local_refs,&
                                      rma_buffers_orig%left_img_map_rows,&
                                      rma_buffers_orig%left_img_map_cols,&
                                      rma_buffers_orig%left_img_refs_rows(2,:),&
                                      rma_buffers_orig%left_img_refs_cols(2,:),&
                                      dbcsr_row_block_sizes(matrix_left),&
                                      dbcsr_col_block_sizes(matrix_left),&
                                      dbcsr_distribution_local_rows(imgdist_left%i%main),&
                                      dbcsr_distribution_local_cols(imgdist_left%i%main),&
                                      left_max_norms)
       ALLOCATE (left_all_max_norms(left_row_nimages,left_col_nimages,0:left_npcols-1))
       CALL mp_allgather(left_max_norms,left_all_max_norms,grp_left)
       !
       ALLOCATE(right_max_norms(right_row_nimages,right_col_nimages))
       CALL calculate_max_image_norms(rma_buffers_orig%right_meta,rma_buffers_orig%right_data,&
                                      rma_buffers_orig%right_local_refs,&
                                      rma_buffers_orig%right_img_map_rows,&
                                      rma_buffers_orig%right_img_map_cols,&
                                      rma_buffers_orig%right_img_refs_rows(2,:),&
                                      rma_buffers_orig%right_img_refs_cols(2,:),&
                                      dbcsr_row_block_sizes(matrix_right),&
                                      dbcsr_col_block_sizes(matrix_right),&
                                      dbcsr_distribution_local_rows(imgdist_right%i%main),&
                                      dbcsr_distribution_local_cols(imgdist_right%i%main),&
                                      right_max_norms)
       !
       ALLOCATE (right_all_max_norms(right_row_nimages,right_col_nimages,0:right_nprows-1))
       CALL mp_allgather(right_max_norms,right_all_max_norms,grp_right)
       DEALLOCATE(right_max_norms,left_max_norms)
       ALLOCATE(right_max_norms(right_myfirstvcol:right_myfirstvcol+right_col_nimages-1,0:nvirt_k-1))
       ALLOCATE(left_max_norms(left_myfirstvrow:left_myfirstvrow+left_row_nimages-1,0:nvirt_k-1))
    ENDIF
    !
    ! Translate the all_sizes to account for pre-distribution.  This
    ! is just done to simplify lookups.
    DO src_vrow = 0, nvirt_k-1
       dst_p = src_vrow / right_row_nimages
       dst_irow = MODULO (src_vrow, right_row_nimages)+1
       right_refs(:, :, src_vrow) =&
            right_all_refs(:, dst_irow, :, dst_p)
       IF (otf_filtering.AND.use_mpi_filtering) THEN
          right_max_norms(:, src_vrow) = right_all_max_norms(dst_irow, :, dst_p)
       ENDIF
    ENDDO
    !
    DO src_vcol = 0, nvirt_k-1
       dst_p = src_vcol / left_col_nimages
       dst_icol = MODULO (src_vcol, left_col_nimages)+1
       left_refs(:, :, src_vcol) =&
            left_all_refs(:, :, dst_icol, dst_p)
       IF (otf_filtering.AND.use_mpi_filtering) THEN
          left_max_norms(:, src_vcol) = left_all_max_norms(:, dst_icol, dst_p)
       ENDIF
    ENDDO
    DEALLOCATE(right_all_refs,left_all_refs)
    IF (otf_filtering.AND.use_mpi_filtering) THEN
       DEALLOCATE(right_all_max_norms,left_all_max_norms)
    ENDIF
    !
    ! Count the maximum possible multiplies per row for on-the-fly filtering.
    IF (otf_filtering) THEN
       !
       nblkrows_local = SUM(rma_buffers_orig%left_img_refs_rows(1,:))
       ALLOCATE(row_max_epss(nblkrows_local),&
                row_counts(nblkrows_local),&
                total_row_counts(nblkrows_local))
       ! Each prow member matrix counts the blocks present in each of its rows
       total_row_counts = 0
       index_row_max_epss = 1
!$omp parallel default(none) &
!$omp          private (v_ri, v_ci, blk, nblks, rowi, row_counts, row) &
!$omp          shared (rma_buffers_orig,total_row_counts,left_row_nimages,&
!$omp                  left_col_nimages,grp_left,index_row_max_epss,&
!$omp                  row_max_epss,nblkrows_local,filter_eps)
       row_counts = 0
!$omp do schedule(dynamic) collapse(2)
       DO v_ci = 1, left_col_nimages
          DO v_ri = 1, left_row_nimages
             IF (rma_buffers_orig%left_local_refs(imeta_size,v_ri,v_ci).NE.0) THEN
                nblks = rma_buffers_orig%left_meta(&
                     rma_buffers_orig%left_local_refs(imeta_displ,v_ri,v_ci)+dbcsr_slot_nblks)
                rowi => rma_buffers_orig%left_meta(&
                     rma_buffers_orig%left_meta(&
                         rma_buffers_orig%left_local_refs(imeta_displ,v_ri,v_ci)+dbcsr_slot_coo_l&
                     )+&
                     rma_buffers_orig%left_local_refs(imeta_displ,v_ri,v_ci):&
                     rma_buffers_orig%left_meta(&
                         rma_buffers_orig%left_local_refs(imeta_displ,v_ri,v_ci)+dbcsr_num_slots&
                     )+&
                     rma_buffers_orig%left_local_refs(imeta_displ,v_ri,v_ci):3)
                DO blk = 1, nblks
                   row_counts(rma_buffers_orig%left_img_map_rows(&
                        rowi(blk)+rma_buffers_orig%left_img_refs_rows(2,v_ri))) = &
                        row_counts(rma_buffers_orig%left_img_map_rows(&
                        rowi(blk)+rma_buffers_orig%left_img_refs_rows(2,v_ri))) + 1 
                ENDDO
             ENDIF
          ENDDO
       ENDDO
!$omp end do
!$omp critical(total_row_counts_acc)
       total_row_counts(:) = total_row_counts(:) + row_counts(:)
!$omp end critical(total_row_counts_acc)
!$omp barrier
!$omp master
       ! The counted blocks are then summed up
       CALL mp_sum(total_row_counts, grp_left)
!$omp end master
!$omp barrier
       ! Determine the maximum per-block epsilon
!$omp do reduction(max:index_row_max_epss)
       DO row = 1, nblkrows_local
          row_max_epss(row) = &
               REAL(filter_eps/REAL(MAX(1,total_row_counts(row)),KIND=KIND(row_max_epss)),&
               KIND=KIND(row_max_epss))
          ! Use integers for a fast comparison
          index_row_max_epss = MAX(index_row_max_epss,total_row_counts(row))
       ENDDO
!$omp end do
!$omp end parallel
       min_row_max_epss = REAL(filter_eps/REAL(index_row_max_epss,KIND=KIND(row_max_epss)),&
            KIND=KIND(row_max_epss))
       DEALLOCATE (row_counts, total_row_counts)
    ELSE
       !
       ! The array must be valid when passed to called subroutines.
       ALLOCATE(left_norms(0,min_nimages),right_norms(0,min_nimages),row_max_epss(0))
    ENDIF
    !
    ALLOCATE(multrec(0:nthreads-1))
    !
!$omp parallel default(none) &
!$omp          private (size_guess, ithread) &
!$omp          shared (product_matrix, multrec, &
!$omp                  keep_sparsity, filter_eps, &
!$omp                  row_max_epss, max_nblocks, &
!$omp                  matrix_right, nthreads)
    !
    ! Setup product work areas
    !
    ithread = 0
!$  ithread = OMP_GET_THREAD_NUM()
    !
    ! The work arrays have to be setup
    size_guess = product_matrix%m%wms(ithread+1)%datasize ! Should be minimal
    CALL dbcsr_data_ensure_size(product_matrix%m%wms(ithread+1)%data_area,&
         size_guess)
    CALL dbcsr_data_set_size_referenced (product_matrix%m%wms(ithread+1)%data_area,&
         product_matrix%m%wms(ithread+1)%datasize)
    ! XXXXXXX a quick fix right now, allocation with size 1 might actually not be needed at all,
    !         but something expects this to be associated
    CALL ensure_array_size(product_matrix%m%wms(ithread+1)%row_i, ub=1)
    CALL ensure_array_size(product_matrix%m%wms(ithread+1)%col_i, ub=1)
    CALL ensure_array_size(product_matrix%m%wms(ithread+1)%blk_p, ub=1)
    !
    ALLOCATE(multrec(ithread)%p)
    CALL dbcsr_mm_multrec_init(multrec(ithread)%p,&
         product=product_matrix%m,&
         keep_sparsity=keep_sparsity,&
         eps=filter_eps,&
         row_max_epss=row_max_epss,&
         block_estimate=max_nblocks/nthreads,&
         right_row_blk_size=dbcsr_row_block_sizes(matrix_right))
!$omp end parallel
    !
    ! update capacity of memory-pools
    IF (has_acc) THEN
       CALL dbcsr_mempool_ensure_capacity(memtype_abpanel_1%pool, &
            capacity=left_col_nimages + right_row_nimages)
       CALL dbcsr_mempool_ensure_capacity(memtype_abpanel_2%pool, &
            capacity=left_col_nimages + right_row_nimages)
    ENDIF
    !
    ! Pre-CANNON execution: filtering and resize buffers for communications
    CALL timeset(routineN//"_precannon", handle2)
    !
    ALLOCATE(right_recv_filtered(0:right_col_nimages-1,0:min_nimages-1))
    ALLOCATE(left_recv_filtered(0:left_row_nimages-1,0:min_nimages-1))
    ! By default all clusters are communicated
    right_recv_filtered = left_row_nimages ; left_recv_filtered = right_col_nimages
    !
    right_max_data_size = 0 ; right_max_meta_size = 0
    left_max_data_size = 0 ; left_max_meta_size = 0
    !
    ALLOCATE(left_max_meta_size_merged(0:min_nimages-1))
    left_max_meta_size_merged = dbcsr_num_slots+size_threads_index
    ALLOCATE(right_max_meta_size_merged(0:min_nimages-1))
    right_max_meta_size_merged = dbcsr_num_slots
    !
!$omp parallel &
!$omp default (none) &
!$omp shared (nsteps_k,min_nimages,right_col_nimages,imgdist_right,&
!$omp         otf_filtering,use_mpi_filtering,&
!$omp         right_myfirstvrow,right_myfirstvcol,right_max_norms,&
!$omp         left_row_nimages,imgdist_left,left_myfirstvrow,left_myfirstvcol,&
!$omp         left_max_norms,min_row_max_epss,&
!$omp         left_recv_filtered,right_recv_filtered,right_refs,left_refs,&
!$omp         right_max_data_size,left_max_data_size,&
!$omp         right_max_meta_size,left_max_meta_size,&
!$omp         right_data_size,left_data_size,right_meta_size,left_meta_size,&
!$omp         size_threads_index,left_max_meta_size_merged,right_max_meta_size_merged) &
!$omp private (metronome,v_ki,v_ci,v_ri,recv_vrow,recv_vcol,right_vrow,&
!$omp          right_vcol,left_vrow,left_vcol,meta_size_image)
    DO metronome = 0, nsteps_k-1
       !
!$omp master
       right_data_size = 0
       left_data_size = 0
       right_meta_size = 0
       left_meta_size = 0
!$omp end master
!$omp barrier
       !
       ! Take first cluster global virtual coordinates
       CALL image_calculator (imgdist_right,&
            vprow=right_vrow, vpcol=right_vcol,&
            myvprow=right_myfirstvrow,&
            myvpcol=right_myfirstvcol,&
            vprow_shift=metronome*min_nimages,&
            shifting='R')
       CALL image_calculator (imgdist_left,&
            vprow=left_vrow,vpcol=left_vcol,&
            myvprow=left_myfirstvrow,&
            myvpcol=left_myfirstvcol,&
            vpcol_shift=metronome*min_nimages,&
            shifting='L')
       !
!$omp do reduction(+:right_data_size,left_data_size,right_meta_size,left_meta_size)
       DO v_ki = 0, min_nimages-1
          !
          recv_vrow = v_ki+right_vrow
          recv_vcol = v_ki+left_vcol
          !
          ! Max norms evaluation for filtering
          IF (otf_filtering.AND.use_mpi_filtering) THEN
             !
             ! Reset clusters counters
             left_recv_filtered(:,v_ki) = 0
             right_recv_filtered(:,v_ki) = 0
             !
             ! multiplication of norms
             DO v_ri = 0, left_row_nimages-1
                IF (left_refs(imeta_size,v_ri+left_vrow,recv_vcol).EQ.0) CYCLE
                DO v_ci = 0, right_col_nimages-1
                   IF (right_refs(imeta_size,v_ci+right_vcol,recv_vrow).EQ.0) CYCLE
                   IF (left_max_norms(v_ri+left_vrow,recv_vcol)*right_max_norms(v_ci+right_vcol,recv_vrow).LT.&
                       min_row_max_epss) CYCLE
                   left_recv_filtered(v_ri,v_ki) = left_recv_filtered(v_ri,v_ki)+1
                   right_recv_filtered(v_ci,v_ki) = right_recv_filtered(v_ci,v_ki)+1
                ENDDO
             ENDDO
          ENDIF
          !
          ! resize buffers taking in account filtering
          meta_size_image = 0
          DO v_ci = 0, right_col_nimages-1
             IF (right_recv_filtered(v_ci,v_ki).EQ.0.OR.&
                 right_refs(imeta_size,v_ci+right_vcol,recv_vrow).EQ.0) THEN
                meta_size_image = meta_size_image + dbcsr_num_slots
                right_refs(idata_size,v_ci+right_vcol,recv_vrow) = 0
                right_refs(imeta_size,v_ci+right_vcol,recv_vrow) = 0
             ELSE
                right_data_size = right_data_size + right_refs(idata_size,v_ci+right_vcol,recv_vrow)
                meta_size_image = meta_size_image + right_refs(imeta_size,v_ci+right_vcol,recv_vrow)
             ENDIF
          ENDDO
          right_meta_size = right_meta_size + meta_size_image
          right_max_meta_size_merged(v_ki) = MAX(right_max_meta_size_merged(v_ki),meta_size_image)
          !
          meta_size_image = 0
          DO v_ri = 0, left_row_nimages-1
             IF (left_recv_filtered(v_ri,v_ki).EQ.0.OR.&
                 left_refs(imeta_size,v_ri+left_vrow,recv_vcol).EQ.0) THEN
                meta_size_image = meta_size_image + dbcsr_num_slots+size_threads_index
                left_refs(idata_size,v_ri+left_vrow,recv_vcol) = 0
                left_refs(imeta_size,v_ri+left_vrow,recv_vcol) = 0
             ELSE
                left_data_size = left_data_size + left_refs(idata_size,v_ri+left_vrow,recv_vcol)
                meta_size_image = meta_size_image + left_refs(imeta_size,v_ri+left_vrow,recv_vcol)
             ENDIF
          ENDDO
          left_meta_size = left_meta_size + meta_size_image
          left_max_meta_size_merged(v_ki) = MAX(left_max_meta_size_merged(v_ki),meta_size_image)
       ENDDO
!$omp end do
!$omp master
       right_max_data_size = MAX(right_max_data_size,right_data_size)
       left_max_data_size = MAX(left_max_data_size,left_data_size)
       right_max_meta_size = MAX(right_max_meta_size,right_meta_size)
       left_max_meta_size = MAX(left_max_meta_size,left_meta_size)
!$omp end master
!$omp barrier
    ENDDO
!$omp end parallel
    IF (otf_filtering.AND.use_mpi_filtering) THEN
       DEALLOCATE(right_max_norms,left_max_norms)
    ENDIF
    DEALLOCATE(left_recv_filtered,right_recv_filtered)
    !
    CALL timestop(handle2)
    !
    ! Set RMA windows for the original data
    ! Right:
    CALL dbcsr_win_create_any(rma_buffers_orig%right_data,grp_right,right_data_win)
    CALL mp_win_create(rma_buffers_orig%right_meta,grp_right,right_meta_win)
    CALL mp_win_lock_all(right_data_win)
    CALL mp_win_lock_all(right_meta_win)
    ! Left:
    CALL dbcsr_win_create_any(rma_buffers_orig%left_data,grp_left,left_data_win)
    CALL mp_win_create(rma_buffers_orig%left_meta,grp_left,left_meta_win)
    CALL mp_win_lock_all(left_data_win)
    CALL mp_win_lock_all(left_meta_win)
    !
    ! Prepare buffers for computation
    IF (nsteps_k.GT.1) THEN
       CALL rma_buffers_init(rma_buffers_1,data_type,&
                             left_max_data_size,right_max_data_size,&
                             left_max_meta_size,right_max_meta_size,&
                             data_memory_type=memtype_abpanel_1)
       ! Right:
       CALL setup_buffer_matrices_images (rma_buffers_1%right_buffer, min_nimages,&
                                          imgdist_right, matrix_right,&
                                          right_max_meta_size_merged,&
                                          rma_buffers_1%right_data)
       ALLOCATE(rma_buffers_1%right_offset_data(min_nimages*right_col_nimages+1))
       rma_buffers_1%right_offset_data(1) = 1
       ! Left:
       CALL setup_buffer_matrices_images (rma_buffers_1%left_buffer, min_nimages,&
                                          imgdist_left, matrix_left,&
                                          left_max_meta_size_merged,&
                                          rma_buffers_1%left_data)
       ALLOCATE(rma_buffers_1%left_offset_data(left_row_nimages*min_nimages+1))
       rma_buffers_1%left_offset_data(1) = 1
       !
       IF (has_acc) THEN
          CALL dbcsr_data_init(rma_buffers_1%trs_stackbuf)
          CALL dbcsr_data_new(rma_buffers_1%trs_stackbuf,&
               data_type=dbcsr_type_int_4, data_size=1000, &
               memory_type=memtype_trsbuffer_1)
       ENDIF
       !
       rma_buffers_calc => rma_buffers_1
    ENDIF
    !
    ! Prepare buffers for communication
    CALL rma_buffers_init(rma_buffers_2,data_type,&
                          left_max_data_size,right_max_data_size,&
                          left_max_meta_size,right_max_meta_size,&
                          data_memory_type=memtype_abpanel_2)
    ! Right:
    CALL setup_buffer_matrices_images (rma_buffers_2%right_buffer, min_nimages,&
                                       imgdist_right, matrix_right,&
                                       right_max_meta_size_merged,&
                                       rma_buffers_2%right_data)
    ALLOCATE(rma_buffers_2%right_offset_data(min_nimages*right_col_nimages+1))
    rma_buffers_2%right_offset_data(1) = 1
    ! Left:
    CALL setup_buffer_matrices_images (rma_buffers_2%left_buffer, min_nimages,&
                                       imgdist_left, matrix_left,&
                                       left_max_meta_size_merged,&
                                       rma_buffers_2%left_data)
    ALLOCATE(rma_buffers_2%left_offset_data(left_row_nimages*min_nimages+1))
    rma_buffers_2%left_offset_data(1) = 1
    !
    IF (has_acc) THEN
       CALL dbcsr_data_init(rma_buffers_2%trs_stackbuf)
       CALL dbcsr_data_new(rma_buffers_2%trs_stackbuf,&
            data_type=dbcsr_type_int_4, data_size=1000, &
            memory_type=memtype_trsbuffer_2)
    ENDIF
    !
    rma_buffers_comm => rma_buffers_2
    !
    ! Preallocate norms arrays
    IF (otf_filtering) THEN
       ALLOCATE(left_norms((MAXVAL(left_max_meta_size_merged)-&
                            dbcsr_num_slots-size_threads_index)/3,min_nimages))
       ALLOCATE(right_norms((MAXVAL(right_max_meta_size_merged)-&
                             dbcsr_num_slots)/3,min_nimages))
    ENDIF
    DEALLOCATE(right_max_meta_size_merged,left_max_meta_size_merged)
    !
    ALLOCATE (right_get_requests(2,min_nimages*right_col_nimages))
    right_get_requests = mp_request_null
    ALLOCATE (left_get_requests(2,left_row_nimages*min_nimages))
    left_get_requests = mp_request_null           
    !
    ! Setup the receive data pointers
    CALL dbcsr_data_init(data_get)
    CALL dbcsr_data_new(data_get, data_type)
    !
    ! These values for meta data are used for empty images
    right_first_indices(dbcsr_slot_nblkrows_total:dbcsr_slot_nfullcols_local) = &
         (/ &
         matrix_right%m%nblkrows_total, &
         matrix_right%m%nblkcols_total, &
         matrix_right%m%nfullrows_total, &
         matrix_right%m%nfullcols_total, &
         0, 0, &
         matrix_right%m%nfullrows_local, &
         matrix_right%m%nfullcols_local /)
    left_first_indices(dbcsr_slot_nblkrows_total:dbcsr_slot_nfullcols_local) = &
         (/ &
         matrix_left%m%nblkrows_total, &
         matrix_left%m%nblkcols_total, &
         matrix_left%m%nfullrows_total, &
         matrix_left%m%nfullcols_total, &
         0, 0, &
         matrix_left%m%nfullrows_local, &
         matrix_left%m%nfullcols_local /)
    !
    ! Here is the main loop
    !
    CALL timeset(routineN//"_loop", handle2)
    !
    grouped_k_index: DO metronome = 0, nsteps_k
       IF (debug_mod) WRITE(*,'(1X,A,3(1X,A,1X,I5))')routineN,&
            "step",metronome,&
            "first k",metronome*min_nimages,&
            "last k",(metronome+1)*min_nimages-1
       ! Wait for right matrix transfer completion. Wait in all but
       ! the first loop iteration.
       wait_right: IF (metronome .GT. 0) THEN
          IF (debug_mod) WRITE (*,'(1X,A)')routineN//" waiting for right"
          !
          CALL mp_waitall (right_get_requests(:,1:right_nrequests))
          right_get_requests(:,1:right_nrequests) = mp_request_null
          ! clean-up wins
          IF (metronome.eq.nsteps_k) THEN
             CALL mp_win_unlock_all(right_data_win)
             CALL mp_win_free(right_data_win)
             CALL mp_win_unlock_all(right_meta_win)
             CALL mp_win_free(right_meta_win)
          ENDIF
       ENDIF wait_right
       ! Right matrix transfer. Transfer in all but the last loop
       ! iteration.
       xfer_right: IF (metronome .LT. nsteps_k) THEN
          !
          IF (has_acc) THEN
             CALL timeset(routineN//"_acc_sync_right", handle3)
             CALL acc_event_synchronize(rma_buffers_comm%right_data%d%acc_ready)
             CALL timestop(handle3)
          ENDIF
          !
          right_meta_offset = 0
          right_nrequests = 0
          !
          ! Take first cluster global virtual coordinates
          CALL image_calculator (imgdist_right,&
               pcol=recv_pcol,&
               vprow=right_vrow,vpcol=right_vcol,&
               myvprow=right_myfirstvrow,&
               myvpcol=right_myfirstvcol,&
               vprow_shift=metronome*min_nimages,&
               shifting='R')
          !
          DO v_ki = 0, min_nimages-1
             v_i = v_ki*right_col_nimages
             recv_vrow = v_ki+right_vrow
             recv_prow = recv_vrow / right_row_nimages
             recv_rowi = MODULO (recv_vrow, right_row_nimages)+1
             DO v_ci = 0, right_col_nimages-1
                recv_vcol = v_ci+right_vcol
                !
                IF (right_refs(imeta_size,recv_vcol,recv_vrow).GT.0) THEN
                   !
                   rma_buffers_comm%right_offset_data(v_ci+2+v_i) = &
                        rma_buffers_comm%right_offset_data(v_ci+1+v_i) + &
                        right_refs(idata_size,recv_vcol,recv_vrow)
                   !
                   meta_get => rma_buffers_comm%right_meta(right_meta_offset+1:&
                        right_meta_offset+right_refs(imeta_size,recv_vcol,recv_vrow))
                   right_meta_offset = right_meta_offset + right_refs(imeta_size,recv_vcol,recv_vrow)
                   !
                   right_nrequests = right_nrequests + 1
                   !
                   CALL dbcsr_data_set_pointer(&
                        area=data_get,&
                        rsize=right_refs(idata_size,recv_vcol,recv_vrow),&
                        csize=1,&
                        pointee=rma_buffers_comm%right_data,&
                        source_lb=rma_buffers_comm%right_offset_data(v_ci+1+v_i))
                   CALL dbcsr_rget_any(data_get,recv_prow,&
                        right_data_win,&
                        rma_buffers_orig%right_data,&
                        right_refs(idata_displ,recv_vcol,recv_vrow),&
                        right_get_requests(1,right_nrequests))
                   CALL mp_rget(meta_get,recv_prow,&
                        right_meta_win,&
                        rma_buffers_orig%right_meta,&
                        right_refs(imeta_displ,recv_vcol,recv_vrow),&
                        right_get_requests(2,right_nrequests))
                   IF (metronome.GT.0) THEN
                      CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(1,:),&
                           right_refs(idata_size,recv_vcol,recv_vrow),&
                           dbcsr_mpi_statistics%data_size_breakdown(:,:,1),data_type_byte)
                      dbcsr_mpi_statistics%nexchanged = dbcsr_mpi_statistics%nexchanged+1
                   ENDIF
                ELSE
                   recv_coli = v_ci+1
                   rma_buffers_comm%right_offset_data(v_ci+2+v_i) = &
                        rma_buffers_comm%right_offset_data(v_ci+1+v_i)
                   right_first_indices(dbcsr_slot_nblkrows_local) = &
                        rma_buffers_orig%right_img_refs_rows(1,recv_rowi)
                   right_first_indices(dbcsr_slot_nblkcols_local) = &
                        rma_buffers_orig%right_img_refs_cols(1,recv_coli)
                   CALL set_empty_meta_index(rma_buffers_comm%right_meta(right_meta_offset+1:&
                                             right_meta_offset+dbcsr_num_slots),&
                                             right_first_indices,&
                                             recv_prow,recv_pcol,recv_rowi,recv_coli,&
                                             recv_vrow,recv_vcol)
                   right_meta_offset = right_meta_offset + dbcsr_num_slots
                   IF (metronome.GT.0) THEN
                      dbcsr_mpi_statistics%nfiltered = dbcsr_mpi_statistics%nfiltered+1
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF xfer_right
       !
       ! Wait for left matrix transfer completion. Wait in all but
       ! the first loop iteration.
       wait_left: IF (metronome .GT. 0) THEN
          IF (debug_mod) WRITE (*,'(1X,A)')routineN//" waiting for left"
          CALL mp_waitall (left_get_requests(:,1:left_nrequests))
          left_get_requests(:,1:left_nrequests) = mp_request_null
          ! clean-up wins
          IF (metronome.eq.nsteps_k) THEN
             CALL mp_win_unlock_all(left_data_win)
             CALL mp_win_free(left_data_win)
             CALL mp_win_unlock_all(left_meta_win)
             CALL mp_win_free(left_meta_win)
          ENDIF
       ENDIF wait_left
       ! Left matrix transfer. Transfer in all but the last loop
       ! iteration.
       xfer_left: IF (metronome .LT. nsteps_k) THEN
          !
          IF (has_acc) THEN
             CALL timeset(routineN//"_acc_sync_left", handle3)
             CALL acc_event_synchronize(rma_buffers_comm%left_data%d%acc_ready)
             CALL timestop(handle3)
          ENDIF
          !
          left_meta_offset = 0
          left_nrequests = 0
          !
          ! Take first cluster global virtual coordinates
          CALL image_calculator (imgdist_left,&
               prow=recv_prow,&
               vprow=left_vrow,vpcol=left_vcol,&
               myvprow=left_myfirstvrow,&
               myvpcol=left_myfirstvcol,&
               vpcol_shift=metronome*min_nimages,&
               shifting='L')
          DO v_ki = 0, min_nimages-1
             v_i = v_ki*left_row_nimages
             recv_vcol = v_ki+left_vcol
             recv_pcol = recv_vcol / left_col_nimages
             recv_coli = MODULO (recv_vcol, left_col_nimages)+1
             DO v_ri = 0, left_row_nimages-1
                recv_vrow = v_ri+left_vrow
                !
                IF (left_refs(imeta_size,recv_vrow,recv_vcol).GT.0) THEN
                   !
                   rma_buffers_comm%left_offset_data(v_ri+2+v_i) = &
                        rma_buffers_comm%left_offset_data(v_ri+1+v_i) + &
                        left_refs(idata_size,recv_vrow,recv_vcol)
                   !
                   meta_get => rma_buffers_comm%left_meta(left_meta_offset+1:&
                        left_meta_offset+left_refs(imeta_size,recv_vrow,recv_vcol))
                   left_meta_offset = left_meta_offset + left_refs(imeta_size,recv_vrow,recv_vcol)
                   !
                   left_nrequests = left_nrequests + 1
                   !                   
                   CALL dbcsr_data_set_pointer(&
                        area=data_get,&
                        rsize=left_refs(idata_size,recv_vrow,recv_vcol),&
                        csize=1,&
                        pointee=rma_buffers_comm%left_data,&
                        source_lb=rma_buffers_comm%left_offset_data(v_ri+1+v_i))
                   CALL dbcsr_rget_any(data_get,recv_pcol,&
                        left_data_win,&
                        rma_buffers_orig%left_data,&
                        left_refs(idata_displ,recv_vrow,recv_vcol),&
                        left_get_requests(1,left_nrequests))
                   CALL mp_rget(meta_get,recv_pcol,&
                        left_meta_win,&
                        rma_buffers_orig%left_meta,&
                        left_refs(imeta_displ,recv_vrow,recv_vcol),&
                        left_get_requests(2,left_nrequests))
                   IF (metronome.GT.0) THEN
                      CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(2,:),&
                           left_refs(idata_size,recv_vrow,recv_vcol),&
                           dbcsr_mpi_statistics%data_size_breakdown(:,:,2),data_type_byte)
                      dbcsr_mpi_statistics%nexchanged = dbcsr_mpi_statistics%nexchanged+1
                   ENDIF
                ELSE
                   recv_rowi = v_ri+1
                   rma_buffers_comm%left_offset_data(v_ri+2+v_i) = &
                        rma_buffers_comm%left_offset_data(v_ri+1+v_i)
                   left_first_indices(dbcsr_slot_nblkrows_local) = &
                        rma_buffers_orig%left_img_refs_rows(1,recv_rowi)
                   left_first_indices(dbcsr_slot_nblkcols_local) = &
                        rma_buffers_orig%left_img_refs_cols(1,recv_coli)
                   CALL set_empty_meta_index(rma_buffers_comm%left_meta(left_meta_offset+1:&
                                             left_meta_offset+dbcsr_num_slots+size_threads_index),&
                                             left_first_indices,&
                                             recv_prow,recv_pcol,recv_rowi,&
                                             recv_coli,recv_vrow,recv_vcol,nthreads)
                   left_meta_offset = left_meta_offset + dbcsr_num_slots + size_threads_index
                   IF (metronome.GT.0) THEN
                      dbcsr_mpi_statistics%nfiltered = dbcsr_mpi_statistics%nfiltered+1
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF xfer_left
       !
       ! Merge and do the multiplications. 
       ! Exclude the first interation
       merge_calc: IF (metronome .GT. 0) THEN
          !
          ! Merge col-images of right matrices
          ! Repoint indices of right matrices
          CALL merge_images(rma_buffers_calc%right_buffer,min_nimages,right_col_nimages,&
               rma_buffers_calc%right_offset_data,rma_buffers_calc%right_meta,&
               imgdist=imgdist_right,do_merge_rows=.FALSE.,&
               vmap_indices=rma_buffers_orig%right_img_map_cols)
          !
          ! Merge col-images of left matrices
          ! Repoint indices of left matrices
          CALL merge_images(rma_buffers_calc%left_buffer,left_row_nimages,min_nimages,&
               rma_buffers_calc%left_offset_data,rma_buffers_calc%left_meta,&
               imgdist=imgdist_left,do_merge_rows=.TRUE.,&
               vmap_indices=rma_buffers_orig%left_img_map_rows,nthreads=nthreads)
          !
          CALL timeset(routineN//"_multrec", handle4)
          !
          IF (has_acc) THEN
             CALL dbcsr_data_host2dev(rma_buffers_calc%right_data)
             CALL dbcsr_data_host2dev(rma_buffers_calc%left_data)
          ENDIF
          !
          DO v_ki = 1, min_nimages
             IF (rma_buffers_calc%right_buffer%mats(v_ki)%m%nblks.EQ.0 .OR. &
                 rma_buffers_calc%left_buffer%mats(v_ki)%m%nblks.EQ.0) CYCLE
             IF (debug_mod) THEN
                CALL dbcsr_print(rma_buffers_calc%right_buffer%mats(v_ki),&
                     nodata=.TRUE.)
                CALL dbcsr_print(rma_buffers_calc%left_buffer%mats(v_ki),&
                     nodata=.TRUE.)
             ENDIF
             !
             IF (has_acc) THEN
                CALL acc_transpose_blocks(rma_buffers_calc%right_buffer%mats(v_ki), &
                     rma_buffers_calc%trs_stackbuf)
             ENDIF
             !
             IF (otf_filtering) THEN
                CALL calculate_norms(rma_buffers_calc%left_buffer%mats(v_ki),&
                     left_norms(:,v_ki))
                CALL calculate_norms(rma_buffers_calc%right_buffer%mats(v_ki),&
                     right_norms(:,v_ki))
             ENDIF
          ENDDO
          !
          threads_finished = 0
          !
!$omp parallel default (none) &
!$omp shared (min_nimages,rma_buffers_calc, &
!$omp         left_norms,right_norms,use_comm_thread, &
!$omp         metronome,nsteps_k,multrec,nthreads,threads_finished, &
!$omp         right_get_requests,left_get_requests,right_nrequests, &
!$omp         left_nrequests) &
!$omp private (ithread,v_ki,threads_finished_read) &
!$omp reduction (+: flop)
          ithread = 0
!$        ithread = omp_get_thread_num()
          !
          IF(metronome==nsteps_k) &
               CALL dbcsr_mm_multrec_phaseout(multrec(ithread)%p)
          !
          DO v_ki = 1, min_nimages
             IF (rma_buffers_calc%right_buffer%mats(v_ki)%m%nblks.EQ.0 .OR. &
                 rma_buffers_calc%left_buffer%mats(v_ki)%m%nblks.EQ.0) CYCLE
             CALL dbcsr_mm_multrec_multiply(multrec(ithread)%p,&
                  left=rma_buffers_calc%left_buffer%mats(v_ki)%m,&
                  right=rma_buffers_calc%right_buffer%mats(v_ki)%m,&
                  flop=flop,&
                  a_norms=left_norms(:,v_ki),&
                  b_norms=right_norms(:,v_ki))
          ENDDO
          !
          IF(metronome==nsteps_k) THEN
             CALL dbcsr_mm_multrec_finalize(multrec(ithread)%p)
             DEALLOCATE(multrec(ithread)%p)
          ENDIF
          !
!$omp atomic
          threads_finished = threads_finished + 1
          IF (use_comm_thread .AND. (nsteps_k.GT.1) .AND. (ithread .EQ. 0)) THEN
             DO
! requires OMP 3.1 (e.g. gcc >=4.7), for correctness, otherwise we keep fingers crossed
#if defined _OPENMP && _OPENMP >= 200711
!$omp atomic read
#endif
                threads_finished_read=threads_finished
                IF (threads_finished_read .EQ. nthreads) EXIT
                CALL mp_testany(right_get_requests(:,1:right_nrequests))
                CALL mp_testany(left_get_requests(:,1:left_nrequests))
             ENDDO
          ENDIF
!$omp end parallel
          !
          CALL timestop(handle4)
          !
       ENDIF merge_calc
       !
       IF (metronome.LT.nsteps_k) THEN
          rma_buffers => rma_buffers_comm
          rma_buffers_comm => rma_buffers_calc
          rma_buffers_calc => rma_buffers
       ENDIF
    ENDDO grouped_k_index
    !
    CALL timestop(handle2)
    !
    CALL m_memory(mem)
    max_memory = MAX(max_memory,REAL(mem))
    !
    IF (debug_mod) THEN
       v_ki = 0
       DO v_i = 1, SIZE(product_matrix%m%blk_p)
          v_ki = MAX(v_ki, ABS(product_matrix%m%blk_p(v_i)))
       ENDDO
       WRITE(*,*)routineN//" Actual final size",&
            LOG(REAL(dbcsr_data_get_size(product_matrix%m%data_area)))/LOG(10.0),&
            LOG(REAL(v_ki))/LOG(10.0)
    ENDIF
    !
    IF (has_acc) THEN
       CALL dbcsr_data_release(rma_buffers_1%trs_stackbuf)
       CALL dbcsr_data_release(rma_buffers_2%trs_stackbuf)
    ENDIF
    !
    DEALLOCATE (right_norms,left_norms,row_max_epss)
    !
    CALL dbcsr_data_clear_pointer(data_get)
    CALL dbcsr_data_release(data_get)
    !
    ! clean-up of communication buffers
    IF (nsteps_k.GT.1) THEN
       DEALLOCATE(rma_buffers_1%right_offset_data, rma_buffers_1%left_offset_data)       
       CALL dbcsr_destroy_array (rma_buffers_1%right_buffer)
       CALL dbcsr_destroy_array (rma_buffers_1%left_buffer)
    ENDIF
    !
    DEALLOCATE(right_refs,left_refs)
    DEALLOCATE(multrec)
    DEALLOCATE(rma_buffers_2%right_offset_data, rma_buffers_2%left_offset_data)
    DEALLOCATE(left_get_requests, right_get_requests)
    CALL dbcsr_destroy_array (rma_buffers_2%right_buffer)
    CALL dbcsr_destroy_array (rma_buffers_2%left_buffer)
    !
    CALL timestop(handle)
  END SUBROUTINE cannon_multiply_low_rma

! *****************************************************************************
!> \brief ...
!> \param mpi_statistics ...
!> \param data_size ...
!> \param size_breakdown ...
!> \param element_size_bytes ...
!> \author Alfio Lazzaro
! *****************************************************************************
  SUBROUTINE count_mpi_statistics(mpi_statistics,data_size,&
                                  size_breakdown,element_size_bytes)
    REAL, DIMENSION(:), INTENT(INOUT)        :: mpi_statistics
    INTEGER, INTENT(IN)                      :: data_size
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: size_breakdown
    INTEGER, INTENT(IN), OPTIONAL            :: element_size_bytes

    INTEGER                                  :: ilimit, nlimits
    INTEGER(KIND=int_8)                      :: data_size_bytes, llimit

    mpi_statistics(1) = mpi_statistics(1) + REAL(data_size)
    mpi_statistics(2) = MIN(mpi_statistics(2),REAL(data_size))
    mpi_statistics(3) = MAX(mpi_statistics(3),REAL(data_size))

    IF (PRESENT(size_breakdown)) THEN
       data_size_bytes=data_size
       ! change in bytes
       IF (PRESENT(element_size_bytes)) data_size_bytes=data_size_bytes*element_size_bytes
       nlimits=SIZE(dbcsr_mpi_size_limits)
       ! check for oversize messages
       IF (data_size_bytes.GT.dbcsr_mpi_size_limits(nlimits)) THEN
          size_breakdown(nlimits+1,1)=size_breakdown(nlimits+1,1)+1
          size_breakdown(nlimits+1,2)=size_breakdown(nlimits+1,2)+data_size_bytes
          RETURN
       ENDIF
       llimit=0
       DO ilimit = 1, nlimits
          IF (data_size_bytes.GE.llimit.AND.data_size_bytes.LE.dbcsr_mpi_size_limits(ilimit)) THEN
             size_breakdown(ilimit,1)=size_breakdown(ilimit,1)+1
             size_breakdown(ilimit,2)=size_breakdown(ilimit,2)+data_size_bytes
             RETURN
          ENDIF
          llimit=dbcsr_mpi_size_limits(ilimit)
       ENDDO
    ENDIF

  END SUBROUTINE count_mpi_statistics

! *****************************************************************************
!> \brief ...
!> \param meta_index ...
!> \param template_index ...
!> \param prow ...
!> \param pcol ...
!> \param rowi ...
!> \param coli ...
!> \param vprow ...
!> \param vpcol ...
!> \param nthreads ...
!> \author Alfio Lazzaro
! *****************************************************************************
  SUBROUTINE set_empty_meta_index(meta_index, template_index, &
                                  prow, pcol, rowi, coli, vprow, vpcol, &
                                  nthreads)
    INTEGER, DIMENSION(:), INTENT(OUT)       :: meta_index
    INTEGER, DIMENSION(:), INTENT(IN)        :: template_index
    INTEGER, INTENT(IN)                      :: prow, pcol, rowi, coli, &
                                                vprow, vpcol
    INTEGER, INTENT(IN), OPTIONAL            :: nthreads

!   ---------------------------------------------------------------------------

    meta_index(dbcsr_slot_size) = dbcsr_num_slots
    IF (PRESENT(nthreads)) THEN
!$     meta_index(dbcsr_slot_size) = meta_index(dbcsr_slot_size)+nthreads+1
    ENDIF
    meta_index(dbcsr_slot_size+1:meta_index(dbcsr_slot_size)) = 0
    meta_index(dbcsr_slot_nblkrows_total:dbcsr_slot_nfullcols_local) = &
         template_index(:)
    meta_index(dbcsr_slot_home_prow) = prow
    meta_index(dbcsr_slot_home_pcol) = pcol
    meta_index(dbcsr_slot_home_rowi) = rowi
    meta_index(dbcsr_slot_home_coli) = coli
    meta_index(dbcsr_slot_home_vprow) = vprow
    meta_index(dbcsr_slot_home_vpcol) = vpcol
    meta_index(dbcsr_slot_row_p) = 1
    meta_index(dbcsr_slot_col_i) = 1
    meta_index(dbcsr_slot_blk_p) = 1
!$ IF (PRESENT(nthreads)) THEN
!$     meta_index(dbcsr_slot_thr_c) = dbcsr_num_slots+1
!$     meta_index(dbcsr_slot_thr_c+1) = meta_index(dbcsr_slot_thr_c)+nthreads
!$  ENDIF
    meta_index(dbcsr_slot_coo_l) = meta_index(dbcsr_slot_size)+1
    meta_index(dbcsr_num_slots) = meta_index(dbcsr_slot_size)
  END SUBROUTINE set_empty_meta_index

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param source_matrix ...
!> \param index_size ...
!> \param data_size ...
!> \param data_buffer ...
!> \param data_memory_type ...
! *****************************************************************************
  SUBROUTINE setup_buffer_matrix (matrix, source_matrix, &
       index_size, data_size, data_buffer, data_memory_type)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_obj), INTENT(IN)              :: source_matrix
    INTEGER, INTENT(IN)                      :: index_size
    INTEGER, INTENT(IN), OPTIONAL            :: data_size
    TYPE(dbcsr_data_obj), INTENT(IN), &
      OPTIONAL                               :: data_buffer
    TYPE(dbcsr_memtype_type), INTENT(IN), &
      OPTIONAL                               :: data_memory_type

    CALL dbcsr_init(matrix)
    CALL dbcsr_create(matrix,&
         template = source_matrix,&
         name = TRIM("Buffer of "//TRIM(source_matrix%m%name)),&
         nze = data_size, &
         data_buffer = data_buffer, &
         data_memory_type = data_memory_type,&
         index_memory_type = memtype_mpi_buffer,&
         make_index=.FALSE.)
    IF (PRESENT(data_size)) THEN
       CALL dbcsr_data_ensure_size (&
            matrix%m%data_area,&
            data_size, nocopy=.TRUE.)
    ENDIF
    CALL ensure_array_size (&
         matrix%m%index,&
         ub=index_size, nocopy=.TRUE.,&
         memory_type=dbcsr_get_index_memory_type(matrix))
    matrix%m%negate_real = source_matrix%m%negate_real
    matrix%m%negate_imaginary = source_matrix%m%negate_imaginary
    matrix%m%local_indexing = source_matrix%m%local_indexing
    matrix%m%list_indexing = source_matrix%m%list_indexing
    !
    IF (source_matrix%m%has_local_rows) THEN
       matrix%m%local_rows = source_matrix%m%local_rows
       CALL array_hold (matrix%m%local_rows)
       matrix%m%has_local_rows = .TRUE.
    ENDIF
    IF (source_matrix%m%has_global_rows) THEN
       matrix%m%global_rows = source_matrix%m%global_rows
       CALL array_hold (matrix%m%global_rows)
       matrix%m%has_global_rows = .TRUE.
    ENDIF
    IF (source_matrix%m%has_local_cols) THEN
       matrix%m%local_cols = source_matrix%m%local_cols
       CALL array_hold (matrix%m%local_cols)
       matrix%m%has_local_cols = .TRUE.
    ENDIF
    IF (source_matrix%m%has_global_cols) THEN
       matrix%m%global_cols = source_matrix%m%global_cols
       CALL array_hold (matrix%m%global_cols)
       matrix%m%has_global_cols = .TRUE.
    ENDIF
    
  END SUBROUTINE setup_buffer_matrix

! *****************************************************************************
!> \brief ...
!> \param buffer_set ...
!> \param nimages ...
!> \param imgdist ...
!> \param template_matrix ...
!> \param index_size_merged ...
!> \param data_buffer ...
!> \author Alfio Lazzaro
! *****************************************************************************
  SUBROUTINE setup_buffer_matrices_images (buffer_set, nimages, imgdist, &
       template_matrix, index_size_merged, data_buffer)
    TYPE(dbcsr_1d_array_type), INTENT(INOUT) :: buffer_set
    INTEGER, INTENT(IN)                      :: nimages
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(INOUT)                          :: imgdist
    TYPE(dbcsr_obj), INTENT(IN)              :: template_matrix
    INTEGER, DIMENSION(:), INTENT(IN)        :: index_size_merged
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: data_buffer

    INTEGER                                  :: v_ki

!   ---------------------------------------------------------------------------

    buffer_set%image_dist = imgdist
    CALL dbcsr_image_dist_hold(imgdist)
    !
    ALLOCATE (buffer_set%mats(nimages))
    DO v_ki = 1, nimages
       CALL dbcsr_init(buffer_set%mats(v_ki))
       CALL dbcsr_create(buffer_set%mats(v_ki),&
            "Buffer image of "//template_matrix%m%name,&
            imgdist%i%main,&
            dbcsr_type_no_symmetry,&
            array_data(template_matrix%m%row_blk_size), array_data(template_matrix%m%col_blk_size),&
            template_matrix%m%row_blk_size, template_matrix%m%col_blk_size,&
            data_buffer = data_buffer,&
            max_rbs=template_matrix%m%max_rbs,max_cbs=template_matrix%m%max_cbs,&
            row_blk_offset=template_matrix%m%row_blk_offset,col_blk_offset=template_matrix%m%col_blk_offset,&
            thread_dist=template_matrix%m%dist,&
            index_memory_type=memtype_mpi_buffer,&
            make_index=.FALSE.)
       CALL ensure_array_size (buffer_set%mats(v_ki)%m%index,&
            ub=index_size_merged(v_ki), nocopy=.TRUE.,&
            memory_type=dbcsr_get_index_memory_type(buffer_set%mats(v_ki)))
       buffer_set%mats(v_ki)%m%negate_real = template_matrix%m%negate_real
       buffer_set%mats(v_ki)%m%negate_imaginary = template_matrix%m%negate_imaginary
       buffer_set%mats(v_ki)%m%local_indexing = .TRUE.
       buffer_set%mats(v_ki)%m%list_indexing = .TRUE.
    ENDDO
  END SUBROUTINE setup_buffer_matrices_images

! *****************************************************************************
!> \brief ...
!> \param buffer_set ...
!> \param buff_nrows ...
!> \param buff_ncols ...
!> \param source_matrix ...
!> \param index_size ...
!> \param data_size ...
! *****************************************************************************
  SUBROUTINE setup_buffer_matrices (buffer_set, buff_nrows, buff_ncols,&
       source_matrix, index_size, data_size)
    TYPE(dbcsr_2d_array_type), INTENT(OUT)   :: buffer_set
    INTEGER, INTENT(IN)                      :: buff_nrows, buff_ncols
    TYPE(dbcsr_obj), INTENT(IN)              :: source_matrix
    INTEGER, INTENT(IN)                      :: index_size, data_size

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_buffer_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col_image, handle, row_image

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, handle)

    CALL dbcsr_image_dist_init (buffer_set%image_dist)
    ALLOCATE (buffer_set%mats(buff_nrows, buff_ncols))
    DO row_image = 1, buff_nrows
       DO col_image = 1, buff_ncols
          CALL setup_buffer_matrix(buffer_set%mats(row_image, col_image),&
               source_matrix, index_size, data_size=data_size,&
               data_memory_type=memtype_abpanel_2)
       ENDDO
    ENDDO
    IF (source_matrix%m%local_indexing .AND. careful_mod) THEN
       IF(.NOT.array_exists (source_matrix%m%local_rows))&
          CPABORT("Local rows should exist.")
       IF(.NOT.array_exists (source_matrix%m%global_rows))&
          CPABORT("Global rows should exist.")
       !
       IF(.NOT.array_exists (source_matrix%m%local_cols))&
          CPABORT("Local cols should exist.")
       IF(.NOT.array_exists (source_matrix%m%global_cols))&
          CPABORT("Global cols should exist.")
    ENDIF
    CALL timestop(handle)
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
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE buffer_matrices_ensure_size(buffer_set, index_size, data_size)
    TYPE(dbcsr_2d_array_type), INTENT(INOUT) :: buffer_set
    INTEGER, INTENT(IN)                      :: index_size, data_size

    CHARACTER(len=*), PARAMETER :: routineN = 'buffer_matrices_ensure_size', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col_image, handle, row_image

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, handle)

    DO row_image = 1, SIZE(buffer_set%mats, 1)
       DO col_image = 1, SIZE(buffer_set%mats, 2)
          CALL dbcsr_data_ensure_size (buffer_set%mats(row_image, col_image)%m%data_area,&
               data_size)
          CALL ensure_array_size (buffer_set%mats(row_image, col_image)%m%index,&
               ub=index_size,&
               memory_type=dbcsr_get_index_memory_type(buffer_set%mats(row_image, col_image)))
          CALL dbcsr_repoint_index(buffer_set%mats(row_image, col_image)%m)    
       ENDDO
    ENDDO
    CALL timestop(handle)
  END SUBROUTINE buffer_matrices_ensure_size

! *****************************************************************************
!> \brief ...
!> \param matrix_set ...
!> \param ns ...
!> \author Alfio Lazzaro
! *****************************************************************************
  SUBROUTINE setup_rec_index_1d (matrix_set, ns)
    TYPE(dbcsr_1d_array_type), INTENT(INOUT) :: matrix_set
    INTEGER, INTENT(IN)                      :: ns

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_rec_index_1d', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: handle, in, t_f, t_l, t_size

!$  INTEGER                                  :: ithread
    LOGICAL                                  :: thread_redist

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, handle)
    DO in = 1, ns
       IF (dbg) THEN
          WRITE(*,*)routineN//" m, n, size",&
               SIZE(matrix_set%mats(in)%m%coo_l),&
               dbcsr_nblkrows_local(matrix_set%mats(in)),&
               dbcsr_nblkcols_local(matrix_set%mats(in)),&
               matrix_set%mats(in)%m%nblks
       ENDIF
       IF (careful_mod) THEN
          IF(SIZE(matrix_set%mats(in)%m%coo_l) .NE. matrix_set%mats(in)%m%nblks*3)&
               CPABORT("Block count mismatch.")
       ENDIF
       thread_redist = ASSOCIATED (matrix_set%mats(in)%m%thr_c)
       t_size = SIZE(matrix_set%mats(in)%m%coo_l)/3
       t_f = 1
       t_l = t_size
       !$OMP PARALLEL IF (thread_redist) DEFAULT (none) &
       !$OMP PRIVATE (ithread) &
       !$OMP FIRSTPRIVATE (t_f, t_l, t_size) &
       !$OMP SHARED (matrix_set, in, thread_redist)
       !$ ithread = OMP_GET_THREAD_NUM()
       !$ IF (thread_redist) THEN
       !$    t_f = matrix_set%mats(in)%m%thr_c(ithread+1)+1
       !$    t_l = matrix_set%mats(in)%m%thr_c(ithread+2)
       !$ ENDIF
       t_size =  t_l - t_f + 1
       !$OMP BARRIER
       IF (t_size .GT. 0) THEN
          CALL rec_sort_index(1,dbcsr_nblkrows_local(matrix_set%mats(in)),&
               1, dbcsr_nblkcols_local(matrix_set%mats(in)),&
               t_size, matrix_set%mats(in)%m%coo_l((t_f*3-2):(t_l*3)),0)
       ENDIF
       !$OMP END PARALLEL
    ENDDO
    CALL timestop(handle)
  END SUBROUTINE setup_rec_index_1d

! *****************************************************************************
!> \brief ...
!> \param matrix_set ...
!> \param n_rows ...
!> \param n_cols ...
! *****************************************************************************
  SUBROUTINE setup_rec_index_2d (matrix_set, n_rows, n_cols)
    TYPE(dbcsr_2d_array_type), INTENT(INOUT) :: matrix_set
    INTEGER, INTENT(IN)                      :: n_rows, n_cols

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_rec_index_2d', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: handle, i_col, i_row, t_f, &
                                                t_l, t_size

!$  INTEGER                                  :: ithread
    LOGICAL                                  :: thread_redist

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, handle)
    DO i_row = 1, n_rows
       DO i_col = 1, n_cols
          IF (.FALSE.) &
               CALL dbcsr_reset_vlocals (matrix_set%mats(i_row, i_col),&
               matrix_set%image_dist)
          IF (dbg) THEN
             WRITE(*,*)routineN//" m, n, size",&
                  SIZE(matrix_set%mats(i_row, i_col)%m%coo_l),&
                  dbcsr_nblkrows_local(matrix_set%mats(i_row, i_col)),&
                  dbcsr_nblkcols_local(matrix_set%mats(i_row, i_col))
             WRITE(*,'(3(1X,I7))')matrix_set%mats(i_row, i_col)%m%coo_l
          ENDIF
          IF (careful_mod) THEN
             IF(SIZE(matrix_set%mats(i_row, i_col)%m%coo_l) .NE. matrix_set%mats(i_row, i_col)%m%nblks*3)&
                CPABORT("Block count mismatch.")
          ENDIF
          thread_redist = ASSOCIATED (matrix_set%mats(i_row, i_col)%m%thr_c)
          t_size = SIZE(matrix_set%mats(i_row, i_col)%m%coo_l)/3
          t_f = 1
          t_l = t_size
          !$OMP PARALLEL IF (thread_redist) DEFAULT (none) &
          !$OMP PRIVATE (ithread) &
          !$OMP FIRSTPRIVATE (t_f, t_l, t_size) &
          !$OMP SHARED (matrix_set, i_row, i_col, thread_redist)
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
                     matrix_set%mats(i_row, i_col)%m%coo_l((t_f*3-2):(t_l*3)))
             ELSE
                CALL call_rec_sort_index (&
                     dbcsr_nblkrows_total(matrix_set%mats(i_row, i_col)),&
                     dbcsr_nblkcols_total(matrix_set%mats(i_row, i_col)),&
                     t_size,&
                     matrix_set%mats(i_row, i_col)%m%coo_l((t_f*3-2):(t_l*3)))
             ENDIF
          ENDIF
          !$OMP END PARALLEL
       ENDDO
    ENDDO
    CALL timestop(handle)
  END SUBROUTINE setup_rec_index_2d

! *****************************************************************************
!> \brief Used to thunk a call to rec_sort_index
!> \param m ...
!> \param n ...
!> \param nblks ...
!> \param idx ...
! *****************************************************************************
  SUBROUTINE call_rec_sort_index (m,n,nblks,idx)
    INTEGER, INTENT(IN)                      :: m, n, nblks
    INTEGER, DIMENSION(3, 1:nblks), &
      INTENT(INOUT)                          :: idx

    CHARACTER(len=*), PARAMETER :: routineN = 'call_rec_sort_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

!   ---------------------------------------------------------------------------

    CALL timeset (routineN, handle)
    IF (.FALSE.) WRITE(*,*)" Calling rec_sort_index, size", nblks
    CALL rec_sort_index(1, m, 1, n, nblks, idx, 0)
    CALL timestop (handle)
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

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rec_split', &
      routineP = moduleN//':'//routineN

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
    CPASSERT(p_high .EQ. nlow)

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
!> \param[out] norms         Block norms
! *****************************************************************************
  SUBROUTINE calculate_norms(matrix, norms)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    REAL(kind=sp), DIMENSION(:), INTENT(OUT) :: norms

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_norms', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_type, handle, nblks, &
                                                nrows
    INTEGER, DIMENSION(1), TARGET            :: tmp
    INTEGER, DIMENSION(:), POINTER           :: local_cols, local_rows

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, handle)
    ! Checks for validity
    IF(.NOT.dbcsr_valid_index (matrix))&
       CPABORT("The matrix must be valid.")
    data_type = dbcsr_get_data_type (matrix)
    IF (matrix%m%local_indexing) THEN
       IF (careful_mod .AND. .NOT.array_exists (matrix%m%local_rows))&
           CPABORT("Global row mapping should exist")
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
          CPABORT("Invalid data type.")
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
          CPABORT("Invalid data type.")
       END SELECT
    ENDIF
    !
    CALL timestop(handle)
  END SUBROUTINE calculate_norms

! *****************************************************************************
!> \brief Calculates max norms per each image
!> \param[in] meta           Meta data array
!> \param DATA ...
!> \param[in] refs           Data/meta offsets/sizes
!> \param img_map_rows ...
!> \param img_map_cols ...
!> \param img_offset_rows ...
!> \param img_offset_cols ...
!> \param row_blk_size ...
!> \param col_blk_size ...
!> \param local_rows ...
!> \param local_cols ...
!> \param[out] max_norms     Max norms per image array
! *****************************************************************************
  SUBROUTINE calculate_max_image_norms(meta,DATA,refs,&
                                       img_map_rows, img_map_cols,&
                                       img_offset_rows,img_offset_cols,&
                                       row_blk_size,col_blk_size,&
                                       local_rows,local_cols,&
                                       max_norms)
    INTEGER, DIMENSION(:), INTENT(IN)        :: meta
    TYPE(dbcsr_data_obj), INTENT(IN)         :: DATA
    INTEGER, DIMENSION(:, :, :), INTENT(IN)  :: refs
    INTEGER, DIMENSION(:), INTENT(IN) :: img_map_rows, img_map_cols, &
      img_offset_rows, img_offset_cols, row_blk_size, col_blk_size, &
      local_rows, local_cols
    REAL(kind=sp), ALLOCATABLE, &
      DIMENSION(:, :), INTENT(INOUT)         :: max_norms

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_max_image_norms', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN, handle)
    !
    ! Check for row/columns image sizes
    IF(.NOT.ALLOCATED(max_norms))&
       CPABORT("Norms array not allocated!")
    IF(SIZE(refs,2).NE.SIZE(max_norms,1).OR.&
       SIZE(refs,3).NE.SIZE(max_norms,2))&
       CPABORT("Wrong number of images!")
    !
    SELECT CASE (data%d%data_type)
    CASE (dbcsr_type_real_4)
       CALL calc_max_image_norms_s(meta,dbcsr_get_data_p_s(DATA),&
            refs,&
            img_map_rows, img_map_cols,&
            img_offset_rows,img_offset_cols,&
            row_blk_size,col_blk_size,&
            local_rows,local_cols,&
            max_norms)
    CASE (dbcsr_type_real_8)
       CALL calc_max_image_norms_d(meta,dbcsr_get_data_p_d(DATA),&
            refs,&
            img_map_rows, img_map_cols,&
            img_offset_rows,img_offset_cols,&
            row_blk_size,col_blk_size,&
            local_rows,local_cols,&
            max_norms)
    CASE (dbcsr_type_complex_4)
       CALL calc_max_image_norms_c(meta,dbcsr_get_data_p_c(DATA),&
            refs,&
            img_map_rows, img_map_cols,&
            img_offset_rows,img_offset_cols,&
            row_blk_size,col_blk_size,&
            local_rows,local_cols,&
            max_norms)
    CASE (dbcsr_type_complex_8)
       CALL calc_max_image_norms_z(meta,dbcsr_get_data_p_z(DATA),&
            refs,&
            img_map_rows, img_map_cols,&
            img_offset_rows,img_offset_cols,&
            row_blk_size,col_blk_size,&
            local_rows,local_cols,&
            max_norms)
    CASE DEFAULT
       CPABORT("Invalid data type.")
    END SELECT
    !
    CALL timestop(handle)
  END SUBROUTINE calculate_max_image_norms

! *****************************************************************************
!> \brief write out a stack for transposing the blocks
!> \param matrix ...
!> \param trs_stackbuf ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_transpose_blocks(matrix, trs_stackbuf)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: trs_stackbuf

    CHARACTER(len=*), PARAMETER :: routineN = 'acc_transpose_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_p, col, handle, i, m, mi, &
                                                mi_max, n, nblks, ni, ni_max, &
                                                offset, row, x
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: counters, filled, offsets
    LOGICAL                                  :: local
    INTEGER, DIMENSION(:), POINTER :: blk_index, col_blk_sizes, &
      col_blk_sizes2enum, enum2col_blk_sizes, enum2row_blk_sizes, &
      local2global_cols, local2global_rows, row_blk_sizes, &
      row_blk_sizes2enum, trs_stack
    INTEGER, DIMENSION(3, matrix%m%nblks)    :: tmp_stack

    CALL timeset(routineN, handle)

    NULLIFY(row_blk_sizes2enum, enum2row_blk_sizes)
    NULLIFY(col_blk_sizes2enum, enum2col_blk_sizes)
    NULLIFY(local2global_rows, local2global_cols, trs_stack)

    IF (.NOT. matrix%m%list_indexing)&
       CPABORT("build_trs_stack: only list_indexing supported.")
    IF (trs_stackbuf%d%data_type /= dbcsr_type_int_4)&
       CPABORT("build_trs_stack: stac_buf has wrong datatype")
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

    CALL dbcsr_data_ensure_size(trs_stackbuf, data_size=nblks, nocopy=.TRUE.)
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
         IF(filled(mi,ni) /= counters(mi,ni))&
            CPABORT("acc_transpose_blocks: bug")
      END DO
    END DO

    !transfer all stacks
    CALL dbcsr_data_host2dev(trs_stackbuf)

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
            stream=trs_stackbuf%d%memory_type%acc_stream)
        END IF
      ENDDO
    ENDDO

    ! make sure block-buffer are not used until transpose kernels finnished
    CALL acc_event_record(trs_stackbuf%d%acc_ready, trs_stackbuf%d%memory_type%acc_stream)
    CALL acc_stream_wait_event(matrix%m%data_area%d%memory_type%acc_stream, trs_stackbuf%d%acc_ready)
    CALL acc_event_record(matrix%m%data_area%d%acc_ready, matrix%m%data_area%d%memory_type%acc_stream)

    DEALLOCATE(row_blk_sizes2enum, enum2row_blk_sizes)
    DEALLOCATE(col_blk_sizes2enum, enum2col_blk_sizes)
    CALL timestop(handle)
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

! *****************************************************************************
!> \brief Init RMA buffer
!>
!> \param buffer ...
!> \param data_type ...
!> \param left_data_size ...
!> \param right_data_size ...
!> \param left_meta_size ...
!> \param right_meta_size ...
!> \param left_local_refs ...
!> \param right_local_refs ...
!> \param left_img_map_rows ...
!> \param left_img_map_cols ...
!> \param right_img_map_rows ...
!> \param right_img_map_cols ...
!> \param left_img_refs_rows ...
!> \param left_img_refs_cols ...
!> \param right_img_refs_rows ...
!> \param right_img_refs_cols ...
!> \param data_memory_type ...
!> \author Alfio Lazzaro
! *****************************************************************************
  SUBROUTINE rma_buffers_init(buffer, data_type, &
                              left_data_size, right_data_size, &
                              left_meta_size, right_meta_size, &
                              left_local_refs, right_local_refs, &
                              left_img_map_rows, left_img_map_cols, &
                              right_img_map_rows, right_img_map_cols, &
                              left_img_refs_rows, left_img_refs_cols, &
                              right_img_refs_rows, right_img_refs_cols, &
                              data_memory_type)
    TYPE(dbcsr_rma_buffers), INTENT(INOUT)   :: buffer
    INTEGER, INTENT(IN)                      :: data_type
    INTEGER, INTENT(IN), OPTIONAL            :: left_data_size, &
                                                right_data_size, &
                                                left_meta_size, &
                                                right_meta_size
    INTEGER, DIMENSION(:, :, :), &
      INTENT(IN), OPTIONAL, POINTER          :: left_local_refs, &
                                                right_local_refs
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL, POINTER                      :: left_img_map_rows, &
                                                left_img_map_cols, &
                                                right_img_map_rows, &
                                                right_img_map_cols
    INTEGER, DIMENSION(:, :), INTENT(IN), &
      OPTIONAL, POINTER                      :: left_img_refs_rows, &
                                                left_img_refs_cols, &
                                                right_img_refs_rows, &
                                                right_img_refs_cols
    TYPE(dbcsr_memtype_type), INTENT(IN)     :: data_memory_type

    IF (PRESENT(left_data_size)) THEN
       IF (.NOT.buffer%is_valid_left) THEN
          ! First initialization
          CALL dbcsr_data_init(buffer%left_data)
          CALL dbcsr_data_new(buffer%left_data,data_type=data_type,&
                              data_size=left_data_size,memory_type=data_memory_type)
          CALL dbcsr_data_set_size_referenced(buffer%left_data, left_data_size)
          buffer%is_valid_left = .TRUE.
       ELSE
          CALL dbcsr_data_ensure_size(buffer%left_data,left_data_size,nocopy=.TRUE.)
       ENDIF
    ENDIF

    IF (PRESENT(right_data_size)) THEN
       IF (.NOT.buffer%is_valid_right) THEN
          ! First initialization
          CALL dbcsr_data_init(buffer%right_data)
          CALL dbcsr_data_new(buffer%right_data,data_type=data_type,&
                              data_size=right_data_size,memory_type=data_memory_type)
          CALL dbcsr_data_set_size_referenced(buffer%right_data, right_data_size)
          buffer%is_valid_right = .TRUE.
       ELSE
          CALL dbcsr_data_ensure_size(buffer%right_data,right_data_size,nocopy=.TRUE.)
       ENDIF
    ENDIF
    
    IF (PRESENT(left_meta_size)) THEN
       CALL ensure_array_size(buffer%left_meta,ub=left_meta_size,nocopy=.TRUE.,&
                              memory_type=memtype_mpi_buffer)
    ENDIF

    IF (PRESENT(right_meta_size)) THEN
       CALL ensure_array_size(buffer%right_meta,ub=right_meta_size,nocopy=.TRUE.,&
                              memory_type=memtype_mpi_buffer)
    ENDIF

    IF (PRESENT(left_local_refs)) THEN
       IF (ASSOCIATED(buffer%left_local_refs)) DEALLOCATE(buffer%left_local_refs)
       buffer%left_local_refs => left_local_refs
    ENDIF
    IF (PRESENT(right_local_refs)) THEN
       IF (ASSOCIATED(buffer%right_local_refs)) DEALLOCATE(buffer%right_local_refs)
       buffer%right_local_refs => right_local_refs
    ENDIF

    IF (PRESENT(left_img_map_rows)) THEN
       IF (ASSOCIATED(buffer%left_img_map_rows)) DEALLOCATE(buffer%left_img_map_rows)
       buffer%left_img_map_rows => left_img_map_rows
    ENDIF
    IF (PRESENT(left_img_map_cols)) THEN
       IF (ASSOCIATED(buffer%left_img_map_cols)) DEALLOCATE(buffer%left_img_map_cols)
       buffer%left_img_map_cols => left_img_map_cols
    ENDIF
    IF (PRESENT(right_img_map_rows)) THEN
       IF (ASSOCIATED(buffer%right_img_map_rows)) DEALLOCATE(buffer%right_img_map_rows)
       buffer%right_img_map_rows => right_img_map_rows
    ENDIF
    IF (PRESENT(right_img_map_cols)) THEN
       IF (ASSOCIATED(buffer%right_img_map_cols)) DEALLOCATE(buffer%right_img_map_cols)
       buffer%right_img_map_cols => right_img_map_cols
    ENDIF

    IF (PRESENT(left_img_refs_rows)) THEN
       IF (ASSOCIATED(buffer%left_img_refs_rows)) DEALLOCATE(buffer%left_img_refs_rows)
       buffer%left_img_refs_rows => left_img_refs_rows
    ENDIF
    IF (PRESENT(left_img_refs_cols)) THEN
       IF (ASSOCIATED(buffer%left_img_refs_cols)) DEALLOCATE(buffer%left_img_refs_cols)
       buffer%left_img_refs_cols => left_img_refs_cols
    ENDIF
    IF (PRESENT(right_img_refs_rows)) THEN
       IF (ASSOCIATED(buffer%right_img_refs_rows)) DEALLOCATE(buffer%right_img_refs_rows)
       buffer%right_img_refs_rows => right_img_refs_rows
    ENDIF
    IF (PRESENT(right_img_refs_cols)) THEN
       IF (ASSOCIATED(buffer%right_img_refs_cols)) DEALLOCATE(buffer%right_img_refs_cols)
       buffer%right_img_refs_cols => right_img_refs_cols
    ENDIF

    NULLIFY(buffer%left_offset_data)
    NULLIFY(buffer%right_offset_data)

  END SUBROUTINE rma_buffers_init

! *****************************************************************************
!> \brief Release RMA buffers
!>
!> \param buffer ...
!> \author Alfio Lazzaro
! *****************************************************************************
  SUBROUTINE rma_buffers_release(buffer)
    TYPE(dbcsr_rma_buffers), INTENT(INOUT)   :: buffer

    IF (buffer%is_valid_left) THEN
       CALL dbcsr_data_release(buffer%left_data)
       buffer%is_valid_left = .FALSE.
    ENDIF
    IF (ASSOCIATED(buffer%left_meta)) THEN
       DEALLOCATE(buffer%left_meta)
       NULLIFY(buffer%left_meta)
    ENDIF
    IF (buffer%is_valid_right) THEN
       CALL dbcsr_data_release(buffer%right_data)
       buffer%is_valid_right = .FALSE.
    ENDIF
    IF (ASSOCIATED(buffer%right_meta)) THEN
       DEALLOCATE(buffer%right_meta)
       NULLIFY(buffer%right_meta)
    ENDIF       
    
    IF (ASSOCIATED(buffer%left_local_refs)) THEN
       DEALLOCATE(buffer%left_local_refs)
       NULLIFY(buffer%left_local_refs)
    ENDIF
    IF (ASSOCIATED(buffer%right_local_refs)) THEN
       DEALLOCATE(buffer%right_local_refs)
       NULLIFY(buffer%right_local_refs)
    ENDIF

    IF (ASSOCIATED(buffer%left_img_map_rows)) THEN
       DEALLOCATE(buffer%left_img_map_rows)
       NULLIFY(buffer%left_img_map_rows)
    ENDIF
    IF (ASSOCIATED(buffer%left_img_map_cols)) THEN
       DEALLOCATE(buffer%left_img_map_cols)
       NULLIFY(buffer%left_img_map_cols)
    ENDIF
    IF (ASSOCIATED(buffer%right_img_map_rows)) THEN
       DEALLOCATE(buffer%right_img_map_rows)
       NULLIFY(buffer%right_img_map_rows)
    ENDIF
    IF (ASSOCIATED(buffer%right_img_map_cols)) THEN
       DEALLOCATE(buffer%right_img_map_cols)
       NULLIFY(buffer%right_img_map_cols)
    ENDIF
    
    IF (ASSOCIATED(buffer%left_img_refs_rows)) THEN
       DEALLOCATE(buffer%left_img_refs_rows)
       NULLIFY(buffer%left_img_refs_rows)
    ENDIF
    IF (ASSOCIATED(buffer%left_img_refs_cols)) THEN
       DEALLOCATE(buffer%left_img_refs_cols)
       NULLIFY(buffer%left_img_refs_cols)
    ENDIF
    IF (ASSOCIATED(buffer%right_img_refs_rows)) THEN
       DEALLOCATE(buffer%right_img_refs_rows)
       NULLIFY(buffer%right_img_refs_rows)
    ENDIF
    IF (ASSOCIATED(buffer%right_img_refs_cols)) THEN
       DEALLOCATE(buffer%right_img_refs_cols)
       NULLIFY(buffer%right_img_refs_cols)
    ENDIF

    NULLIFY(buffer%left_offset_data,buffer%right_offset_data)

  END SUBROUTINE rma_buffers_release

! *****************************************************************************
!> \brief Merge images
!>
!> \param images ...
!> \param nrows_images ...
!> \param ncols_images ...
!> \param offset_images ...
!> \param meta_buffer ...
!> \param imgdist ...
!> \param do_merge_rows ...
!> \param vmap_indices ...
!> \param nthreads ...
!> \author Alfio Lazzaro
! *****************************************************************************
  SUBROUTINE merge_images(images, nrows_images, ncols_images, offset_images,&
                          meta_buffer, imgdist, do_merge_rows,&
                          vmap_indices, nthreads)
    TYPE(dbcsr_1d_array_type), INTENT(INOUT) :: images
    INTEGER, INTENT(IN)                      :: nrows_images, ncols_images
    INTEGER, DIMENSION(:), INTENT(IN)        :: offset_images
    INTEGER, DIMENSION(:), INTENT(IN), &
      TARGET                                 :: meta_buffer
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(INOUT)                          :: imgdist
    LOGICAL, INTENT(IN)                      :: do_merge_rows
    INTEGER, DIMENSION(:), INTENT(IN), &
      POINTER                                :: vmap_indices
    INTEGER, INTENT(IN), OPTIONAL            :: nthreads

    CHARACTER(len=*), PARAMETER :: routineN = 'merge_images', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, index_offset, ithread, ki, kl, ks, mi, my_nthreads, &
      ncols_merged, nimages_merged, nimages_unmerged, nrows_merged, &
      offset_data, offset_meta, offset_meta_first, size_prefix_meta, &
      slot_nblk_local, ui
    INTEGER, DIMENSION(:), POINTER           :: unmerged_indices
    TYPE(dbcsr_obj), POINTER                 :: image_obj
    TYPE(dbcsr_type), POINTER                :: image

    CALL timeset(routineN, handle)
    !
    IF (do_merge_rows) THEN
       ! merge over rows (left)
       nimages_unmerged = ncols_images
       nimages_merged = nrows_images
       slot_nblk_local = dbcsr_slot_nblkrows_local
       index_offset = 1
       nrows_merged = 1
       ncols_merged = nimages_unmerged
    ELSE
       ! merge over cols (right)
       nimages_unmerged = nrows_images
       nimages_merged = ncols_images
       slot_nblk_local = dbcsr_slot_nblkcols_local
       index_offset = 2
       nrows_merged = nimages_unmerged
       ncols_merged = 1
    ENDIF
    !
    offset_meta_first = 0
    merge: IF (nimages_merged.GT.1) THEN
       size_prefix_meta = dbcsr_num_slots
       my_nthreads = 1
       IF (PRESENT(nthreads)) THEN
!$        size_prefix_meta = size_prefix_meta+nthreads+1
          my_nthreads = nthreads
       ENDIF
       !
       ! Merge images
       images_unmerge: DO ui = 1, nimages_unmerged
          image => images%mats(ui)%m
          !
          image%index(dbcsr_slot_size) = size_prefix_meta
          image%index(dbcsr_slot_size+1:size_prefix_meta) = 0
          image%index(dbcsr_num_slots) = size_prefix_meta
          !
          ! Merge corresponding thread blocks for all images
          threads: DO ithread = 1, my_nthreads
             !
             ! Merge images per a given thread accross all images
             image%index(slot_nblk_local) = 0
             offset_meta = offset_meta_first
             images_merge: DO mi = 1, nimages_merged
                unmerged_indices => meta_buffer(offset_meta+dbcsr_slot_size:&
                     offset_meta+meta_buffer(offset_meta+dbcsr_slot_size))
                offset_meta = offset_meta + unmerged_indices(dbcsr_slot_size)
                offset_data = offset_images(mi+(ui-1)*nimages_merged)
                IF (unmerged_indices(dbcsr_slot_nblks).NE.0) THEN
                   IF (ithread.EQ.1) THEN
                      ! nblks, nze
                      image%index(dbcsr_slot_nblks) = &
                           image%index(dbcsr_slot_nblks)+&
                           unmerged_indices(dbcsr_slot_nblks)
                      image%index(dbcsr_slot_nze) = &
                           image%index(dbcsr_slot_nze)+&
                           unmerged_indices(dbcsr_slot_nze)
                      ! threads distribution
!$                    IF (PRESENT(nthreads)) THEN
!$                       image%index(dbcsr_num_slots+1:size_prefix_meta) = &
!$                            image%index(dbcsr_num_slots+1:size_prefix_meta)+&
!$                            unmerged_indices(dbcsr_num_slots+1:size_prefix_meta)
!$                    ENDIF
                   ENDIF
                   ! Copy CSR indices (3 values per block: row-index, col-index, offset)
                   ! remap indices taking in account threads and images
                   ks = size_prefix_meta+1
                   kl = unmerged_indices(dbcsr_slot_size)
!$                 IF (PRESENT(nthreads)) THEN
!$                    ks = ks+unmerged_indices(dbcsr_num_slots+ithread)*3
!$                    kl = size_prefix_meta+unmerged_indices(dbcsr_num_slots+ithread+1)*3
!$                 ENDIF
                   DO ki = ks, kl, 3
                      ! Indices remapping
                      image%index(image%index(dbcsr_slot_size)+index_offset) = &
                           vmap_indices(unmerged_indices(ki+index_offset-1)+&
                           image%index(slot_nblk_local))
                      ! Leave as it is
                      image%index(image%index(dbcsr_slot_size)+3-index_offset) = &
                           unmerged_indices(ki+2-index_offset)
                      ! Data offset
                      image%index(image%index(dbcsr_slot_size)+3) = &
                           unmerged_indices(ki+2)+offset_data-1
                      ! increase meta size
                      image%index(dbcsr_slot_size) = image%index(dbcsr_slot_size)+3
                   ENDDO
                ENDIF
                image%index(slot_nblk_local) = &
                     image%index(slot_nblk_local)+&
                     unmerged_indices(slot_nblk_local)
             ENDDO images_merge
             !
             CALL dbcsr_sort_indices((image%index(dbcsr_slot_size)-image%index(dbcsr_num_slots))/3,&
                  image%index(image%index(dbcsr_num_slots)+1:image%index(dbcsr_slot_size):3),&
                  image%index(image%index(dbcsr_num_slots)+2:image%index(dbcsr_slot_size):3),&
                  blk_d=image%index(image%index(dbcsr_num_slots)+3:image%index(dbcsr_slot_size):3))
             image%index(dbcsr_num_slots) = image%index(dbcsr_slot_size)
          ENDDO threads
          !
          ! Total stats
          image%index(dbcsr_slot_nblkrows_total) = &
               meta_buffer(offset_meta_first+dbcsr_slot_nblkrows_total)
          image%index(dbcsr_slot_nblkcols_total) = &
               meta_buffer(offset_meta_first+dbcsr_slot_nblkcols_total)
          image%index(dbcsr_slot_nfullrows_total) = &
               meta_buffer(offset_meta_first+dbcsr_slot_nfullrows_total)
          image%index(dbcsr_slot_nfullcols_total) = &
               meta_buffer(offset_meta_first+dbcsr_slot_nfullcols_total)
          !
          IF (do_merge_rows) THEN
             ! Local cols
             image%index(dbcsr_slot_nblkcols_local) = &
                  meta_buffer(offset_meta_first+dbcsr_slot_nblkcols_local)
             ! Virtual coords
             image%index(dbcsr_slot_home_vprow) = &
                  meta_buffer(offset_meta_first+dbcsr_slot_home_vprow)/nimages_unmerged
             image%index(dbcsr_slot_home_vpcol) = &
                  meta_buffer(offset_meta_first+dbcsr_slot_home_vpcol)
          ELSE
             ! Local rows
             image%index(dbcsr_slot_nblkrows_local) = &
                  meta_buffer(offset_meta_first+dbcsr_slot_nblkrows_local)
             ! Virtual coords
             image%index(dbcsr_slot_home_vprow) = &
                  meta_buffer(offset_meta_first+dbcsr_slot_home_vprow)
             image%index(dbcsr_slot_home_vpcol) = &
                  meta_buffer(offset_meta_first+dbcsr_slot_home_vpcol)/nimages_unmerged
          ENDIF
          !
          ! Full coords local
          image%index(dbcsr_slot_nfullrows_local) = &
               meta_buffer(offset_meta_first+dbcsr_slot_nfullrows_local)
          image%index(dbcsr_slot_nfullcols_local) = &
               meta_buffer(offset_meta_first+dbcsr_slot_nfullcols_local)
          ! Coords
          image%index(dbcsr_slot_home_prow) = &
               meta_buffer(offset_meta_first+dbcsr_slot_home_prow)
          image%index(dbcsr_slot_home_pcol) = &
               meta_buffer(offset_meta_first+dbcsr_slot_home_pcol)
          image%index(dbcsr_slot_home_rowi) = &
               meta_buffer(offset_meta_first+dbcsr_slot_home_rowi)
          image%index(dbcsr_slot_home_coli) = &
               meta_buffer(offset_meta_first+dbcsr_slot_home_coli)
          ! row_p, col_i, blk_p
          image%index(dbcsr_slot_row_p:dbcsr_slot_row_p+1) = &
               meta_buffer(offset_meta_first+dbcsr_slot_row_p:offset_meta_first+dbcsr_slot_row_p+1)
          image%index(dbcsr_slot_col_i:dbcsr_slot_col_i+1) = &
               meta_buffer(offset_meta_first+dbcsr_slot_col_i:offset_meta_first+dbcsr_slot_col_i+1)
          image%index(dbcsr_slot_blk_p:dbcsr_slot_blk_p+1) = &
               meta_buffer(offset_meta_first+dbcsr_slot_blk_p:offset_meta_first+dbcsr_slot_blk_p+1)
          ! thr_c
!$        IF (PRESENT(nthreads)) THEN
!$           image%index(dbcsr_slot_thr_c) = &
!$                meta_buffer(offset_meta_first+dbcsr_slot_thr_c)
!$           image%index(dbcsr_slot_thr_c+1) = &
!$                meta_buffer(offset_meta_first+dbcsr_slot_thr_c+1)
!$        ENDIF
          ! coo_l
          image%index(dbcsr_slot_coo_l) = &
               meta_buffer(offset_meta_first+dbcsr_slot_coo_l)
          !
          offset_meta_first = offset_meta
       ENDDO images_unmerge
    ELSE
       ! no merging, just copy
       nomerge: DO ui = 1, nimages_unmerged
          image => images%mats(ui)%m
          offset_meta = offset_meta_first + meta_buffer(offset_meta_first+dbcsr_slot_size)
          image%index(dbcsr_slot_size:offset_meta) = &
               meta_buffer(offset_meta_first+dbcsr_slot_size:offset_meta)
          offset_meta_first = offset_meta
       ENDDO nomerge
    ENDIF merge
    !
    ! repoint index and reset
    DO ui = 1, nimages_unmerged
       image_obj => images%mats(ui)
!$     IF (PRESENT(nthreads).AND.image_obj%m%index(dbcsr_slot_nblks).EQ.0.AND.nimages_merged.EQ.1) THEN
!$       image_obj%m%index(image_obj%m%index(dbcsr_slot_thr_c):image_obj%m%index(dbcsr_slot_thr_c+1)) = 0
!$     ENDIF
       CALL dbcsr_repoint_index (image_obj%m)
       image_obj%m%valid = .TRUE.
       CALL dbcsr_reset_vlocals (image_obj,imgdist,&
            do_rows=.NOT.do_merge_rows)
    ENDDO
    !
    ! setup indexing
    CALL setup_rec_index (images, nimages_unmerged)
    !
    CALL timestop(handle)
  END SUBROUTINE merge_images

#include "dbcsr_mm_cannon_d.f90"
#include "dbcsr_mm_cannon_z.f90"
#include "dbcsr_mm_cannon_s.f90"
#include "dbcsr_mm_cannon_c.f90"

END MODULE dbcsr_mm_cannon

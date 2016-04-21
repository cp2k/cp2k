!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
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
! **************************************************************************************************
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
                                              use_mpi_exp,&
                                              use_mpi_filtering
   USE dbcsr_data_methods,              ONLY: &
        dbcsr_data_clear_pointer, dbcsr_data_ensure_size, dbcsr_data_get_size, &
        dbcsr_data_get_size_referenced, dbcsr_data_hold, dbcsr_data_host2dev, dbcsr_data_init, &
        dbcsr_data_new, dbcsr_data_release, dbcsr_data_set_pointer, &
        dbcsr_data_set_size_referenced, dbcsr_data_valid, dbcsr_get_data_p_c, dbcsr_get_data_p_d, &
        dbcsr_get_data_p_s, dbcsr_get_data_p_z, dbcsr_scalar, dbcsr_scalar_are_equal, &
        dbcsr_scalar_fill_all, dbcsr_scalar_negative, dbcsr_scalar_one, dbcsr_scalar_set_type, &
        dbcsr_scalar_zero
   USE dbcsr_data_types,                ONLY: dbcsr_datatype_sizeof
   USE dbcsr_dist_methods,              ONLY: &
        dbcsr_distribution_col_dist, dbcsr_distribution_get_num_images_1d, &
        dbcsr_distribution_has_col_clusters, dbcsr_distribution_has_row_clusters, &
        dbcsr_distribution_has_threads, dbcsr_distribution_hold, dbcsr_distribution_local_cols, &
        dbcsr_distribution_local_rows, dbcsr_distribution_make_threads, dbcsr_distribution_mp, &
        dbcsr_distribution_ncols, dbcsr_distribution_no_threads, dbcsr_distribution_nrows, &
        dbcsr_distribution_num_local_col_clusters, dbcsr_distribution_num_local_row_clusters, &
        dbcsr_distribution_release, dbcsr_distribution_row_dist, dbcsr_distribution_thread_dist
   USE dbcsr_dist_operations,           ONLY: dbcsr_create_image_dist,&
                                              dbcsr_get_local_vcols,&
                                              dbcsr_get_local_vrows,&
                                              dbcsr_make_dists_dense,&
                                              dbcsr_reset_locals,&
                                              dbcsr_reset_vlocals,&
                                              image_calculator,&
                                              make_sizes_dense
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
        dbcsr_col_block_offsets, dbcsr_col_block_sizes, dbcsr_destroy_array, dbcsr_distribution, &
        dbcsr_get_data_type, dbcsr_get_index_memory_type, dbcsr_get_matrix_type, &
        dbcsr_has_symmetry, dbcsr_image_dist_hold, dbcsr_image_dist_init, &
        dbcsr_image_dist_release, dbcsr_init, dbcsr_nblkcols_local, dbcsr_nblkcols_total, &
        dbcsr_nblkrows_local, dbcsr_nblkrows_total, dbcsr_nfullcols_total, dbcsr_nfullrows_total, &
        dbcsr_release, dbcsr_release_locals, dbcsr_row_block_offsets, dbcsr_row_block_sizes, &
        dbcsr_valid_index
   USE dbcsr_mm_multrec,                ONLY: dbcsr_mm_multrec_finalize,&
                                              dbcsr_mm_multrec_init,&
                                              dbcsr_mm_multrec_lib_finalize,&
                                              dbcsr_mm_multrec_lib_init,&
                                              dbcsr_mm_multrec_multiply,&
                                              dbcsr_mm_multrec_phaseout,&
                                              dbcsr_mm_multrec_type
   USE dbcsr_mp_methods,                ONLY: &
        dbcsr_mp_grid_setup, dbcsr_mp_group, dbcsr_mp_has_subgroups, dbcsr_mp_my_col_group, &
        dbcsr_mp_my_row_group, dbcsr_mp_mynode, dbcsr_mp_mypcol, dbcsr_mp_myprow, dbcsr_mp_npcols, &
        dbcsr_mp_nprows, dbcsr_mp_numnodes, dbcsr_mp_pgrid
   USE dbcsr_mp_operations,             ONLY: &
        dbcsr_gatherv_any, dbcsr_ibcast_any, dbcsr_irecv_any, dbcsr_iscatterv_any, &
        dbcsr_isend_any, dbcsr_rget_any, dbcsr_sendrecv_any, dbcsr_win_create_any, &
        hybrid_alltoall_any, hybrid_alltoall_i1
   USE dbcsr_operations,                ONLY: dbcsr_add,&
                                              dbcsr_conjg,&
                                              dbcsr_copy,&
                                              dbcsr_crop_matrix,&
                                              dbcsr_filter,&
                                              dbcsr_may_be_dense,&
                                              dbcsr_scale
   USE dbcsr_ptr_util,                  ONLY: ensure_array_size,&
                                              memory_copy
   USE dbcsr_toollib,                   ONLY: ordered_search,&
                                              uppercase
   USE dbcsr_transformations,           ONLY: dbcsr_make_dense,&
                                              dbcsr_make_dense_low,&
                                              dbcsr_make_undense,&
                                              dbcsr_make_untransposed_blocks,&
                                              dbcsr_new_transposed
   USE dbcsr_types,                     ONLY: &
        dbcsr_1d_array_type, dbcsr_2d_array_type, dbcsr_conjugate_transpose, dbcsr_data_obj, &
        dbcsr_distribution_obj, dbcsr_imagedistribution_obj, dbcsr_iterator, dbcsr_memtype_type, &
        dbcsr_meta_size, dbcsr_mp_obj, dbcsr_mpi_size_limits, dbcsr_mpi_statistics_type, &
        dbcsr_no_transpose, dbcsr_num_slots, dbcsr_obj, dbcsr_scalar_type, dbcsr_slot_blk_p, &
        dbcsr_slot_col_i, dbcsr_slot_coo_l, dbcsr_slot_home_coli, dbcsr_slot_home_pcol, &
        dbcsr_slot_home_prow, dbcsr_slot_home_rowi, dbcsr_slot_home_vpcol, dbcsr_slot_home_vprow, &
        dbcsr_slot_nblkcols_local, dbcsr_slot_nblkrows_total, dbcsr_slot_nblks, &
        dbcsr_slot_nfullcols_local, dbcsr_slot_nze, dbcsr_slot_row_p, dbcsr_slot_size, &
        dbcsr_slot_thr_c, dbcsr_transpose, dbcsr_type, dbcsr_type_antisymmetric, &
        dbcsr_type_complex_4, dbcsr_type_complex_8, dbcsr_type_int_4, dbcsr_type_no_symmetry, &
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
        mp_allgather, mp_alltoall, mp_bcast, mp_gather, mp_gatherv, mp_iallgather, mp_ibcast, &
        mp_irecv, mp_iscatter, mp_isend, mp_isum, mp_max, mp_min, mp_request_null, mp_rget, &
        mp_sendrecv, mp_sum, mp_testall, mp_testany, mp_wait, mp_waitall, mp_win_create, &
        mp_win_free, mp_win_lock_all, mp_win_unlock_all
#include "../../base/base_uses.f90"

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads,&
!$                    omp_set_lock, omp_unset_lock, omp_init_lock, omp_lock_kind, omp_destroy_lock

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_cannon'
   CHARACTER(len=*), PARAMETER, PRIVATE :: int_print = "(10(1X,I7))"
   LOGICAL, PARAMETER :: debug_mod = .FALSE.
   LOGICAL, PARAMETER :: careful_mod = .FALSE.
   REAL, PARAMETER :: huge_norm = HUGE(1.0)**(1.0/3.0)

   TYPE dbcsr_memtype_type_p
      TYPE(dbcsr_memtype_type), POINTER :: p => Null()
      ! ensure that array-elements are on different cache lines
      INTEGER(kind=int_4), DIMENSION(64)    :: padding
   END TYPE dbcsr_memtype_type_p

   TYPE(dbcsr_memtype_type_p), DIMENSION(:), POINTER, SAVE :: memtype_product_wm => Null(), &
                                                              memtype_product_wm_scatter => Null()

   TYPE(dbcsr_mpi_statistics_type), PRIVATE, SAVE :: dbcsr_mpi_statistics
   INTEGER, PRIVATE, SAVE :: max_nblocks = 0
   REAL, PRIVATE, SAVE :: marketing_flops = 0
   REAL, PRIVATE, SAVE :: max_memory = 0

   INTEGER, PRIVATE, SAVE :: last_mpi_ranks_used = 0

   TYPE(dbcsr_memtype_type), PRIVATE, SAVE  :: memtype_abpanel_1, memtype_abpanel_2, &
                                               memtype_trsbuffer_1, memtype_trsbuffer_2, &
                                               memtype_mpi_buffer
   TYPE(acc_stream_type), PRIVATE, SAVE         :: stream_1, stream_2
   ! ab-panels and streams are shared between all threads

   TYPE dbcsr_rma_buffer
      TYPE(dbcsr_data_obj)               :: DATA, data_diag
      INTEGER                            :: vprow, vpcol
      INTEGER, DIMENSION(:), POINTER     :: meta => Null(), &
                                            meta_diag => Null()
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: get_requests, &
                                               get_requests_map, &
                                               offset
      INTEGER                            :: nrequests
      TYPE(dbcsr_1d_array_type)          :: buffer
      LOGICAL                            :: is_valid = .FALSE.
   END TYPE dbcsr_rma_buffer

   TYPE dbcsr_rma_buffers
      TYPE(dbcsr_rma_buffer)             :: left, right
      TYPE(dbcsr_data_obj)               :: trs_stackbuf
   END TYPE dbcsr_rma_buffers

   ! RMA buffers
   TYPE(dbcsr_rma_buffers), PRIVATE, TARGET, SAVE :: rma_buffers_orig, &
                                                     rma_buffers_1, rma_buffers_2

   INTEGER, PARAMETER, PRIVATE               :: idata_displ = 3, &
                                                idata_size = 1, &
                                                imeta_displ = 4, &
                                                imeta_size = 2, &
                                                ilocal_proc = 1, &
                                                isym_proc = 2

   TYPE(dbcsr_data_obj), SAVE, PRIVATE :: local_data_scatter, local_data_product_scatter
   INTEGER, PRIVATE                    :: local_meta_size_scatter, local_data_size_scatter
   INTEGER, DIMENSION(:), SAVE, POINTER :: meta_scatter => Null(), local_meta_scatter => Null(), &
                                           local_meta_product_scatter => Null()
   INTEGER, ALLOCATABLE, DIMENSION(:), PRIVATE :: left_total_row_counts
   INTEGER, ALLOCATABLE, DIMENSION(:), TARGET, PRIVATE :: left_no_empty_images, &
                                                          left_no_empty_images_displ, &
                                                          right_no_empty_images, &
                                                          right_no_empty_images_displ, &
                                                          left_refs_data_size, &
                                                          right_refs_data_size, &
                                                          left_local_refs_data_size, &
                                                          right_local_refs_data_size
   INTEGER, ALLOCATABLE, DIMENSION(:, :), TARGET, PRIVATE  :: left_refs_displ_unmerged, &
                                                              right_refs_displ_unmerged, &
                                                              left_local_refs_displ_unmerged, &
                                                              right_local_refs_displ_unmerged, &
                                                              left_refs_meta_size, &
                                                              right_refs_meta_size, &
                                                              left_local_refs_meta_size, &
                                                              right_local_refs_meta_size
   REAL(kind=sp), ALLOCATABLE, DIMENSION(:), TARGET, PRIVATE :: left_max_norms, right_max_norms, &
                                                                left_local_max_norms, right_local_max_norms

   INTEGER, PRIVATE                :: request_count_rows
   INTEGER, DIMENSION(8), PRIVATE  :: requests
   INTEGER, DIMENSION(4), PRIVATE  :: requests_diag
   INTEGER, DIMENSION(2), PRIVATE  :: requests_scatter
   LOGICAL, PRIVATE                :: my_use_mpi_exp

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
      MODULE PROCEDURE setup_rec_index_images
      MODULE PROCEDURE setup_rec_index_2d
   END INTERFACE

   PUBLIC :: dbcsr_mm_cannon_lib_init, dbcsr_mm_cannon_lib_finalize
   PUBLIC :: dbcsr_mm_cannon_clear_mempools
   PUBLIC :: dbcsr_mm_cannon_multiply

CONTAINS

! **************************************************************************************************
!> \brief Initialize the library
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbcsr_mm_cannon_lib_init()

      INTEGER                                            :: ithread, nthreads

      nthreads = 1; ithread = 0
!$    nthreads = OMP_GET_NUM_THREADS(); ithread = OMP_GET_THREAD_NUM()

      CALL dbcsr_mm_multrec_lib_init()

!$OMP     MASTER
      dbcsr_mpi_statistics%nimages = -1
      dbcsr_mpi_statistics%nexchanged = 0
      dbcsr_mpi_statistics%nfiltered = 0
      dbcsr_mpi_statistics%data_size = 0
      dbcsr_mpi_statistics%data_size(:, 2) = HUGE(dbcsr_mpi_statistics%data_size(1, 2))
      dbcsr_mpi_statistics%data_size_breakdown = 0

      marketing_flops = 0
      max_memory = 0
      ALLOCATE (memtype_product_wm(0:nthreads-1))
      ALLOCATE (memtype_product_wm_scatter(0:nthreads-1))
!$OMP     END MASTER
!$OMP     BARRIER

      ! Each thread has its own working-matrix and its own mempool
      ALLOCATE (memtype_product_wm(ithread)%p)
      CALL dbcsr_memtype_setup(memtype_product_wm(ithread)%p, has_pool=.TRUE.)
      CALL dbcsr_mempool_ensure_capacity(memtype_product_wm(ithread)%p%pool, capacity=1)
      ALLOCATE (memtype_product_wm_scatter(ithread)%p)
      CALL dbcsr_memtype_setup(memtype_product_wm_scatter(ithread)%p, has_pool=.TRUE.)
      CALL dbcsr_mempool_ensure_capacity(memtype_product_wm_scatter(ithread)%p%pool, capacity=1)
   END SUBROUTINE dbcsr_mm_cannon_lib_init

! **************************************************************************************************
!> \brief Finalize the library
!> \param group ...
!> \param output_unit ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbcsr_mm_cannon_lib_finalize(group, output_unit)
      INTEGER, INTENT(IN)                                :: group, output_unit

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_cannon_lib_finalize', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ilimit, isqrt, isqrt2, ithread
      INTEGER(KIND=int_8)                                :: total_nexchanged, total_nfiltered
      INTEGER(KIND=int_8), &
         DIMENSION(SIZE(dbcsr_mpi_size_limits)+1, 2, 2)  :: total_recv_breakdown
      REAL                                               :: average, total_marketing_flops, &
                                                            total_max_memory
      REAL, DIMENSION(2)                                 :: max_recv_data, min_recv_data, &
                                                            total_recv_data

      CALL dbcsr_mm_multrec_lib_finalize(group, output_unit)

      ithread = 0
!$    ithread = omp_get_thread_num()

      ! Each thread has its own working-matrix and its own mempool
      IF (ASSOCIATED(memtype_product_wm(ithread)%p%pool)) &
         CALL dbcsr_mempool_destruct(memtype_product_wm(ithread)%p%pool)
      DEALLOCATE (memtype_product_wm(ithread)%p)
      IF (ASSOCIATED(memtype_product_wm_scatter(ithread)%p%pool)) &
         CALL dbcsr_mempool_destruct(memtype_product_wm_scatter(ithread)%p%pool)
      DEALLOCATE (memtype_product_wm_scatter(ithread)%p)

!$OMP      BARRIER
!$OMP      MASTER
      DEALLOCATE (memtype_product_wm)
      DEALLOCATE (memtype_product_wm_scatter)

      ! Deallocate RMA buffers
      CALL rma_buffers_release(rma_buffers_orig)
      CALL rma_buffers_release(rma_buffers_1)
      CALL rma_buffers_release(rma_buffers_2)

      total_max_memory = max_memory
      CALL mp_max(total_max_memory, group)

      total_marketing_flops = marketing_flops
      CALL mp_sum(total_marketing_flops, group)

      total_nexchanged = dbcsr_mpi_statistics%nexchanged
      CALL mp_sum(total_nexchanged, group)

      total_nfiltered = dbcsr_mpi_statistics%nfiltered
      CALL mp_sum(total_nfiltered, group)

      total_recv_data(:) = dbcsr_mpi_statistics%data_size(:, 1)
      CALL mp_sum(total_recv_data, group)

      min_recv_data(:) = dbcsr_mpi_statistics%data_size(:, 2)
      CALL mp_min(min_recv_data, group)

      max_recv_data(:) = dbcsr_mpi_statistics%data_size(:, 3)
      CALL mp_max(max_recv_data, group)

      IF (dbcsr_mpi_statistics%nexchanged .GT. 0) THEN
         average = SUM(total_recv_data(:))/REAL(total_nexchanged)
      ELSE
         average = 0
         min_recv_data = 0
      ENDIF

      total_recv_breakdown(:, :, :) = dbcsr_mpi_statistics%data_size_breakdown(:, :, :)
      CALL mp_sum(total_recv_breakdown, group)

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(A,T30,EN20.6)') " marketing flops", total_marketing_flops

         IF (dbcsr_mpi_statistics%nimages .GT. 0) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("-", 79)
            WRITE (output_unit, '(A,T30,EN20.6)') " max memory usage/rank", total_max_memory
            WRITE (output_unit, '(A,T30,I20)') " # max total images/rank", dbcsr_mpi_statistics%nimages
            WRITE (output_unit, '(A,T30,I20)') " # MPI messages exchanged", total_nexchanged
            WRITE (output_unit, '(A,T30,I20)') " # MPI messages filtered", total_nfiltered
            WRITE (output_unit, '(A)') " MPI messages size (elements):"
            WRITE (output_unit, '(A,T30,EN20.6)') "  total size", &
               SUM(total_recv_data(:))
            WRITE (output_unit, '(A,T30,EN20.6)') "  min size", &
               MINVAL(min_recv_data(:))
            WRITE (output_unit, '(A,T30,EN20.6)') "  max size", &
               MAXVAL(max_recv_data(:))
            WRITE (output_unit, '(A,T30,EN20.6)') "  average size", average

            WRITE (output_unit, '(A)') " MPI breakdown and total messages size (bytes):"
            WRITE (output_unit, '(A,I8,T40,I10,T55,I20)') "             size <= ", dbcsr_mpi_size_limits(1), &
               SUM(total_recv_breakdown(1, 1, :)), SUM(total_recv_breakdown(1, 2, :))
            DO ilimit = 2, SIZE(dbcsr_mpi_size_limits)
               WRITE (output_unit, '(A,I8,A,I8,T40,I10,T55,I20)') "  ", dbcsr_mpi_size_limits(ilimit-1), &
                  " < size <= ", dbcsr_mpi_size_limits(ilimit), &
                  SUM(total_recv_breakdown(ilimit, 1, :)), SUM(total_recv_breakdown(ilimit, 2, :))
            ENDDO
            ilimit = SIZE(dbcsr_mpi_size_limits)
            WRITE (output_unit, '(A,I8,A,T40,I10,T55,I20)') "  ", dbcsr_mpi_size_limits(ilimit), &
               " < size    ", SUM(total_recv_breakdown(ilimit+1, 1, :)), SUM(total_recv_breakdown(ilimit+1, 2, :))
         ENDIF

         isqrt = NINT(SQRT(REAL(last_mpi_ranks_used, KIND=real_8)))
         isqrt2 = NINT(SQRT(REAL(last_mpi_ranks_used*2, KIND=real_8)))
         IF (isqrt*isqrt .NE. last_mpi_ranks_used) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("-", 79)
            WRITE (UNIT=output_unit, FMT="(T2,A)") &
               "Warning: using a non-square number of MPI ranks might lead to poor performance."
            WRITE (UNIT=output_unit, FMT="(T2,A,I0)") &
               "         used ranks: ", last_mpi_ranks_used
            WRITE (UNIT=output_unit, FMT="(T2,A,2(I0,1X))") &
               "         suggested : ", isqrt**2, isqrt2**2
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
      IF (acc_stream_associated(stream_1)) &
         CALL acc_stream_destroy(stream_1)
      IF (acc_stream_associated(stream_2)) &
         CALL acc_stream_destroy(stream_2)
!$OMP      END MASTER
   END SUBROUTINE dbcsr_mm_cannon_lib_finalize

! **************************************************************************************************
!> \brief Deallocate memory contained in mempools
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbcsr_mm_cannon_clear_mempools()

      INTEGER                                            :: ithread

      ithread = 0
!$    ithread = omp_get_thread_num()

      ! Each thread has its own working-matrix and its own mempool
      IF (ASSOCIATED(memtype_product_wm(ithread)%p%pool)) &
         CALL dbcsr_mempool_clear(memtype_product_wm(ithread)%p%pool)
      IF (ASSOCIATED(memtype_product_wm_scatter(ithread)%p%pool)) &
         CALL dbcsr_mempool_clear(memtype_product_wm_scatter(ithread)%p%pool)

!$OMP      MASTER
      IF (ASSOCIATED(memtype_trsbuffer_1%pool)) &
         CALL dbcsr_mempool_clear(memtype_trsbuffer_1%pool)
      IF (ASSOCIATED(memtype_trsbuffer_2%pool)) &
         CALL dbcsr_mempool_clear(memtype_trsbuffer_2%pool)
      IF (ASSOCIATED(memtype_abpanel_1%pool)) &
         CALL dbcsr_mempool_clear(memtype_abpanel_1%pool)
      IF (ASSOCIATED(memtype_abpanel_2%pool)) &
         CALL dbcsr_mempool_clear(memtype_abpanel_2%pool)
!$OMP      END MASTER
   END SUBROUTINE dbcsr_mm_cannon_clear_mempools

! **************************************************************************************************
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
! **************************************************************************************************
   SUBROUTINE dbcsr_mm_cannon_multiply(transa, transb, &
                                       alpha, matrix_a, matrix_b, beta, matrix_c, &
                                       first_row, last_row, first_column, last_column, first_k, last_k, &
                                       retain_sparsity, filter_eps, &
                                       flop)

      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb
      TYPE(dbcsr_scalar_type), INTENT(IN)                :: alpha
      TYPE(dbcsr_obj), INTENT(IN)                        :: matrix_a, matrix_b
      TYPE(dbcsr_scalar_type), INTENT(IN)                :: beta
      TYPE(dbcsr_obj), INTENT(INOUT)                     :: matrix_c
      INTEGER, INTENT(IN), OPTIONAL                      :: first_row, last_row, first_column, &
                                                            last_column, first_k, last_k
      LOGICAL, INTENT(IN), OPTIONAL                      :: retain_sparsity
      REAL(KIND=real_8), INTENT(IN), OPTIONAL            :: filter_eps
      INTEGER(KIND=int_8), INTENT(OUT), OPTIONAL         :: flop

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_cannon_multiply', &
         routineP = moduleN//':'//routineN
      LOGICAL, PARAMETER                                 :: dbg = .FALSE., &
                                                            use_list_indexing = .TRUE., &
                                                            use_local_indexing = .TRUE.
      REAL(real_8), PARAMETER                            :: make_dense_occ_thresh = 1.0_dp

      CHARACTER                                          :: transa_l, transb_l
      INTEGER :: comm, f_col, f_k, f_row, handle, handle2, ithread, l_col, l_k, l_row, &
         nimages_left_rows, nimages_match, nimages_right_cols, numnodes, output_unit
      INTEGER(KIND=int_8)                                :: my_flop
      LOGICAL :: ab_dense, has_clusters, keep_product_data, keep_sparsity, new_left, new_right, &
         product_reindex, release_tdist, use_dense_mult
      REAL(KIND=dp)                                      :: cs
      TYPE(array_i1d_obj) :: dense_col_sizes, dense_k_sizes, dense_row_sizes, k_vmap, m_map, &
         n_map, old_product_col_blk_offsets, old_product_col_blk_sizes, &
         old_product_row_blk_offsets, old_product_row_blk_sizes
      TYPE(dbcsr_2d_array_type), POINTER                 :: m2s_left, m2s_right
      TYPE(dbcsr_distribution_obj)                       :: dense_product_distribution, &
                                                            old_product_distribution
      TYPE(dbcsr_imagedistribution_obj)                  :: dense_rdist_left, dense_rdist_right, &
                                                            rdist_left, rdist_right
      TYPE(dbcsr_obj)                                    :: matrix_left, matrix_right, product_matrix
      TYPE(dbcsr_scalar_type)                            :: eps_any

      CALL timeset(routineN, handle)
      CALL array_nullify(dense_k_sizes)
      CALL array_nullify(dense_col_sizes)
      CALL array_nullify(dense_row_sizes)

      ithread = 0
!$    ithread = OMP_GET_THREAD_NUM()

      ! setup driver-dependent memory-types and their memory-pools ---------------

      ! the ab_buffers are shared by all threads
      IF (has_acc) THEN
         IF (.NOT. acc_stream_associated(stream_1)) THEN
            CALL acc_stream_create(stream_1, "MemCpy (odd ticks)")
            CALL acc_stream_create(stream_2, "MemCpy (even ticks)")
         ENDIF

         CALL dbcsr_memtype_setup(memtype_abpanel_1, has_pool=.TRUE., &
                                  acc_hostalloc=.TRUE., acc_devalloc=.TRUE., acc_stream=stream_1, &
                                  mpi=.TRUE., oversize_factor=default_resize_factor)

         CALL dbcsr_memtype_setup(memtype_abpanel_2, has_pool=.TRUE., &
                                  acc_hostalloc=.TRUE., acc_devalloc=.TRUE., acc_stream=stream_2, &
                                  mpi=.TRUE., oversize_factor=default_resize_factor)

         !TODO: ensure capacity 2/3?
         CALL dbcsr_memtype_setup(memtype_trsbuffer_1, has_pool=.TRUE., &
                                  acc_hostalloc=.TRUE., acc_devalloc=.TRUE., acc_stream=stream_1)
         CALL dbcsr_mempool_ensure_capacity(memtype_trsbuffer_1%pool, capacity=1)

         CALL dbcsr_memtype_setup(memtype_trsbuffer_2, has_pool=.TRUE., &
                                  acc_hostalloc=.TRUE., acc_devalloc=.TRUE., acc_stream=stream_2)
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
      IF (transa_l .NE. dbcsr_no_transpose .AND. &
          transa_l .NE. dbcsr_transpose .AND. &
          transa_l .NE. dbcsr_conjugate_transpose) &
         CPABORT("Invalid transa_l = "//transa_l)

      IF (transb_l .NE. dbcsr_no_transpose .AND. &
          transb_l .NE. dbcsr_transpose .AND. &
          transb_l .NE. dbcsr_conjugate_transpose) &
         CPABORT("Invalid transb_l = "//transb_l)

      IF (dbg) THEN
         WRITE (*, *) '========== MULTIPLICATION ========================'
         CALL dbcsr_verify_matrix(matrix_a)
         CALL dbcsr_verify_matrix(matrix_b)
         CALL dbcsr_verify_matrix(matrix_c)
         WRITE (*, *) routineN//" ABC checksums", &
            dbcsr_checksum(matrix_a), &
            dbcsr_checksum(matrix_b), &
            dbcsr_checksum(matrix_c)
         IF (dbg) THEN
            CALL dbcsr_print(matrix_a, nodata=.TRUE.)
            CALL dbcsr_print(matrix_b, nodata=.TRUE.)
            CALL dbcsr_print(matrix_c, nodata=.TRUE.)
         ENDIF
      ENDIF

      ! transpose/conjg left and/or right matrices if needed
      SELECT CASE (transa_l)
      CASE (dbcsr_no_transpose)
         matrix_left = matrix_a
         new_left = .FALSE.
      CASE (dbcsr_transpose)
         CALL dbcsr_init(matrix_left)
         IF (dbcsr_get_matrix_type(matrix_a) .EQ. dbcsr_type_antisymmetric) THEN
            !
            ! For antisymmetric matrix, we need to do a hard copy
            ! shallow_data_copy=.TRUE. doesnt handle properly antisymm matrices
            CALL dbcsr_new_transposed(matrix_left, matrix_a, &
                                      shallow_data_copy=.FALSE., redistribute=.FALSE., &
                                      transpose_distribution=.FALSE.)
         ELSE
            CALL dbcsr_new_transposed(matrix_left, matrix_a, &
                                      shallow_data_copy=.TRUE., redistribute=.FALSE., &
                                      transpose_distribution=.FALSE.)
         ENDIF
         new_left = .TRUE.
      CASE (dbcsr_conjugate_transpose)
         CALL dbcsr_init(matrix_left)
         CALL dbcsr_new_transposed(matrix_left, matrix_a, &
                                   shallow_data_copy=.FALSE., redistribute=.FALSE., &
                                   transpose_distribution=.FALSE.)
         CALL dbcsr_conjg(matrix_left)
         new_left = .TRUE.
      CASE DEFAULT
         CPABORT("wrong transa_l = "//transa_l)
      END SELECT

      SELECT CASE (transb_l)
      CASE (dbcsr_no_transpose)
         matrix_right = matrix_b
         new_right = .FALSE.
      CASE (dbcsr_transpose)
         CALL dbcsr_init(matrix_right)
         IF (dbcsr_get_matrix_type(matrix_b) .EQ. dbcsr_type_antisymmetric) THEN
            !
            ! For antisymmetric matrix, we need to do a hard copy
            ! shallow_data_copy=.TRUE. doesnt handle properly antisymm matrices
            CALL dbcsr_new_transposed(matrix_right, matrix_b, &
                                      shallow_data_copy=.FALSE., redistribute=.FALSE., &
                                      transpose_distribution=.FALSE.)
         ELSE
            CALL dbcsr_new_transposed(matrix_right, matrix_b, &
                                      shallow_data_copy=.TRUE., redistribute=.FALSE., &
                                      transpose_distribution=.FALSE.)
         ENDIF
         new_right = .TRUE.
      CASE (dbcsr_conjugate_transpose)
         CALL dbcsr_init(matrix_right)
         CALL dbcsr_new_transposed(matrix_right, matrix_b, &
                                   shallow_data_copy=.FALSE., redistribute=.FALSE., &
                                   transpose_distribution=.FALSE.)
         CALL dbcsr_conjg(matrix_right)
         new_right = .TRUE.
      CASE DEFAULT
         CPABORT("wrong transb_l = "//transb_l)
      END SELECT
      !
      ! Ensure matrix compatibility.
      IF (.NOT. array_equality(dbcsr_row_block_offsets(matrix_c), dbcsr_row_block_offsets(matrix_left))) &
         CPABORT("C/A rows not equal")
      IF (.NOT. array_equality(dbcsr_col_block_offsets(matrix_c), dbcsr_col_block_offsets(matrix_right))) &
         CPABORT("C/B columns not equal")
      IF (.NOT. array_equality(dbcsr_col_block_offsets(matrix_left), dbcsr_row_block_offsets(matrix_right))) &
         CPABORT("A cols/B rows not equal")
      !
      ! No dense multiplication when filtering is used.
      use_dense_mult = .NOT. PRESENT(filter_eps)
      IF (has_acc) use_dense_mult = .FALSE.
      !
      ! Temporary check for RMA algorithm
      IF (new_right .OR. new_left .OR. &
          dbcsr_mp_nprows(dbcsr_distribution_mp(matrix_c%m%dist)) .NE. &
          dbcsr_mp_npcols(dbcsr_distribution_mp(matrix_c%m%dist))) THEN
         my_use_mpi_exp = .FALSE.
      ELSE
         my_use_mpi_exp = use_mpi_exp
      ENDIF
      !
      ! No dense multiplication and local indexing when 2D clustering is used.
      IF (my_use_mpi_exp) THEN
         use_dense_mult = .FALSE.
      ENDIF
      ! we skip dense multiply for (anti)symmetric matrices (slowdown for S/H * C)
      IF (use_dense_mult) THEN
         IF (dbcsr_has_symmetry(matrix_left) .OR. &
             dbcsr_has_symmetry(matrix_right)) THEN
            use_dense_mult = .FALSE.
         ELSE
            use_dense_mult = dbcsr_may_be_dense(matrix_left, make_dense_occ_thresh) &
                             .AND. dbcsr_may_be_dense(matrix_right, make_dense_occ_thresh)
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
      IF (PRESENT(first_row)) THEN
         IF (first_row .LT. 1 .OR. first_row .GT. dbcsr_nfullrows_total(matrix_c)) &
            CPABORT("Invalid first row specified")
         f_row = first_row
      ENDIF
      IF (PRESENT(last_row)) THEN
         IF (last_row .GT. dbcsr_nfullrows_total(matrix_c)) &
            CPABORT("Invalid last row specified")
         l_row = last_row
      ENDIF
      IF (PRESENT(first_column)) THEN
         IF (first_column .LT. 1 .OR. first_column .GT. dbcsr_nfullcols_total(matrix_c)) &
            CPABORT("Invalid first col specified")
         f_col = first_column
      ENDIF
      IF (PRESENT(last_column)) THEN
         IF (last_column .GT. dbcsr_nfullcols_total(matrix_c)) &
            CPABORT("Invalid last column specified (C)")
         IF (last_column .GT. dbcsr_nfullcols_total(matrix_right)) &
            CPABORT("Invalid last column specified (B)")
         l_col = last_column
      ENDIF
      IF (PRESENT(first_k)) THEN
         IF (first_k .LT. 1 .OR. first_k .GT. dbcsr_nfullcols_total(matrix_left)) &
            CPABORT("Invalid first k specified (A)")
         f_k = first_k
      ENDIF
      IF (PRESENT(last_k)) THEN
         IF (last_k .GT. dbcsr_nfullcols_total(matrix_left)) &
            CPABORT("Invalid last k specified (A)")
         l_k = last_k
      ENDIF
      !
      ! update statistics (we count marketing flops per MPI rank)
      numnodes = dbcsr_mp_numnodes(dbcsr_distribution_mp(dbcsr_distribution(matrix_c)))
      last_mpi_ranks_used = numnodes
      marketing_flops = marketing_flops+ &
                        (2.0*(l_row-f_row+1.0)*(l_col-f_col+1.0)/numnodes)*(l_k-f_k+1.0)
      !
      ! Now optimize the default submatrix selection values away
      IF (f_row .EQ. 1) f_row = 0
      IF (l_row .EQ. dbcsr_nfullrows_total(matrix_left)) l_row = 0
      IF (f_col .EQ. 1) f_col = 0
      ! The last column must be set if the right and product matrices
      ! differ.
      l_col = MIN(l_col, dbcsr_nfullcols_total(matrix_right))
      l_col = MIN(l_col, dbcsr_nfullcols_total(matrix_c))
      IF (f_col .LE. 1 .AND. &
          l_col .EQ. dbcsr_nfullcols_total(matrix_right) .AND. &
          dbcsr_nfullcols_total(matrix_right) .EQ. &
          dbcsr_nfullcols_total(matrix_c)) l_col = 0
      IF (f_k .EQ. 1) f_k = 0
      IF (l_k .EQ. dbcsr_nfullcols_total(matrix_left)) l_k = 0
      IF (.NOT. PRESENT(last_column) .AND. &
          .NOT. array_equality(dbcsr_col_block_sizes(matrix_right), &
                               dbcsr_col_block_sizes(matrix_c))) THEN
         l_col = MIN(dbcsr_nfullcols_total(matrix_right), &
                     dbcsr_nfullcols_total(matrix_c))
      ENDIF
      IF (f_row .GT. l_row) &
         CPABORT("Last row smaller than first row")
      IF (f_col .GT. l_col) &
         CPABORT("Last col smaller than first col")
      !
      ! Product data needs to be retained when
      ! * beta != 0; or
      ! * there is column limiting (l_col > 0) and the limiting column
      !   is less than the number of full columns in theproduct matrix
      keep_sparsity = .FALSE.
      IF (PRESENT(retain_sparsity)) keep_sparsity = retain_sparsity
      !
      keep_product_data = keep_sparsity &
                          .OR. .NOT. dbcsr_scalar_are_equal(beta, dbcsr_scalar_zero(beta%data_type)) &
                          .OR. (l_col .GT. 0 .AND. l_col .LT. dbcsr_nfullcols_total(matrix_c)) &
                          .OR. (l_row .GT. 0 .AND. l_row .LT. dbcsr_nfullrows_total(matrix_c))
      !
      IF (.NOT. dbcsr_scalar_are_equal(beta, dbcsr_scalar_one(beta%data_type)) .AND. keep_product_data) THEN
         CALL dbcsr_scale(matrix_c, alpha_scalar=beta, &
                          limits=(/f_row, l_row, f_col, l_col/))
      ENDIF
      !
      ! The index of the product matrix is twiddled into canonical form
      ! if it is (anti)symmetric (i.e., rows and columns are where the
      ! row/column distributions say they are). Doing this in advance
      ! makes the local multiply more efficient.
      IF (dbcsr_has_symmetry(matrix_c)) THEN
         product_reindex = .TRUE.
      ELSE
         product_reindex = .FALSE.
      ENDIF
      ! Product can not be made dense; however, A & B may still be made
      ! dense unless previously determined otherwise.
      IF (product_reindex .OR. keep_sparsity) THEN
         use_dense_mult = .FALSE.
      ENDIF
      !
      ! The thread distribution must reflect the current (possibly
      ! dense) distribution
      IF (.NOT. dbcsr_distribution_has_threads(matrix_c%m%dist)) THEN
         release_tdist = .TRUE.
         CALL dbcsr_distribution_make_threads(matrix_c%m%dist)
      ELSE
         release_tdist = .FALSE.
      ENDIF
      !
      ! Compute number of images (rows and columns)
      ! By construction clusters require RMA algo
      has_clusters = .FALSE.
      IF (my_use_mpi_exp .AND. dbcsr_distribution_has_row_clusters(matrix_left%m%dist)) THEN
         has_clusters = .TRUE.
         nimages_left_rows = dbcsr_distribution_num_local_row_clusters(matrix_left%m%dist)* &
                             dbcsr_mp_nprows(dbcsr_distribution_mp(matrix_left%m%dist))
      ELSE
         ! 1D clustering (only RMA support 2D clustering)
         nimages_left_rows = dbcsr_mp_nprows(dbcsr_distribution_mp(matrix_left%m%dist))
      ENDIF
      !
      IF (my_use_mpi_exp .AND. dbcsr_distribution_has_col_clusters(matrix_left%m%dist)) THEN
         has_clusters = .TRUE.
         nimages_match = dbcsr_distribution_num_local_col_clusters(matrix_left%m%dist)* &
                         dbcsr_mp_npcols(dbcsr_distribution_mp(matrix_left%m%dist))
      ELSE IF (my_use_mpi_exp .AND. dbcsr_distribution_has_row_clusters(matrix_right%m%dist)) THEN
         has_clusters = .TRUE.
         nimages_match = dbcsr_distribution_num_local_row_clusters(matrix_right%m%dist)* &
                         dbcsr_mp_nprows(dbcsr_distribution_mp(matrix_right%m%dist))
      ELSE
         ! 1D clustering
         nimages_match = dbcsr_distribution_get_num_images_1d( &
                         dbcsr_nfullcols_total(matrix_left), &
                         dbcsr_nblkcols_total(matrix_left), &
                         dbcsr_mp_nprows(dbcsr_distribution_mp(matrix_left%m%dist)), &
                         dbcsr_mp_npcols(dbcsr_distribution_mp(matrix_left%m%dist)))
      ENDIF
      !
      IF (my_use_mpi_exp .AND. dbcsr_distribution_has_col_clusters(matrix_right%m%dist)) THEN
         has_clusters = .TRUE.
         nimages_right_cols = dbcsr_distribution_num_local_col_clusters(matrix_right%m%dist)* &
                              dbcsr_mp_npcols(dbcsr_distribution_mp(matrix_right%m%dist))
      ELSE
         ! 1D clustering for the old algorithm (only RMA support 2D clustering)
         nimages_right_cols = dbcsr_mp_npcols(dbcsr_distribution_mp(matrix_right%m%dist))
      ENDIF
      !
      ! Create imaged distributions for the multiply.
      CALL dbcsr_create_image_dist(rdist_right, matrix_right%m%dist, &
                                   match_row_nbins=dbcsr_mp_npcols(dbcsr_distribution_mp(matrix_left%m%dist)), &
                                   match_col_nbins=dbcsr_mp_npcols(dbcsr_distribution_mp(matrix_c%m%dist)), &
                                   match_col_pdist=dbcsr_distribution_col_dist(matrix_c%m%dist), &
                                   nimages_rows=nimages_match, &
                                   nimages_cols=nimages_right_cols)
      !
      CALL dbcsr_create_image_dist(rdist_left, matrix_left%m%dist, &
                                   match_row_pdist=dbcsr_distribution_row_dist(matrix_c%m%dist), &
                                   match_row_nbins=dbcsr_mp_nprows(dbcsr_distribution_mp(matrix_c%m%dist)), &
                                   match_col_pdist=dbcsr_distribution_row_dist(rdist_right%i%main), &
                                   match_col_idist=array_data(rdist_right%i%row_image), &
                                   match_col_nbins=dbcsr_mp_nprows(dbcsr_distribution_mp(matrix_right%m%dist)), &
                                   nimages_rows=nimages_left_rows, &
                                   nimages_cols=nimages_match)
      !
      IF (ab_dense) THEN
         CALL dbcsr_make_dists_dense(dbcsr_distribution(matrix_c), &
                                     rdist_left, rdist_right, dense_product_distribution, &
                                     dense_rdist_left, dense_rdist_right,.NOT. use_dense_mult, &
                                     m_map, k_vmap, n_map, matrix_c%m%row_blk_size)
         CALL make_sizes_dense(matrix_c%m%row_blk_size, m_map, &
                               dbcsr_distribution_nrows(dense_product_distribution), &
                               dense_row_sizes)
         CALL make_sizes_dense(matrix_c%m%col_blk_size, n_map, &
                               dbcsr_distribution_ncols(dense_product_distribution), &
                               dense_col_sizes)
         CALL make_sizes_dense(matrix_right%m%row_blk_size, k_vmap, &
                               dbcsr_distribution_nrows(dense_rdist_right%i%main), &
                               dense_k_sizes)
      ENDIF
      !
      IF (use_dense_mult .AND. .NOT. ab_dense) &
         CPABORT("Wrong logic when making dense matrices.")
      IF (use_dense_mult) THEN
         old_product_row_blk_offsets = matrix_c%m%row_blk_offset
         old_product_col_blk_offsets = matrix_c%m%col_blk_offset
         old_product_row_blk_sizes = matrix_c%m%row_blk_size
         old_product_col_blk_sizes = matrix_c%m%col_blk_size
         CALL array_hold(old_product_row_blk_offsets)
         CALL array_hold(old_product_col_blk_offsets)
         CALL array_hold(old_product_row_blk_sizes)
         CALL array_hold(old_product_col_blk_sizes)
         old_product_distribution = dbcsr_distribution(matrix_c)
         CALL dbcsr_distribution_hold(old_product_distribution)
         CALL dbcsr_init(product_matrix)
         CALL dbcsr_make_dense(matrix_c, product_matrix, &
                               dense_product_distribution, &
                               dense_row_sizes, dense_col_sizes, &
                               m_map, n_map)
      ELSE
         CALL dbcsr_init(product_matrix)
         CALL dbcsr_copy(product_matrix, matrix_c, shallow_data=.TRUE.)
      ENDIF
      IF (ab_dense) THEN
         CALL dbcsr_distribution_release(dense_product_distribution)
      ENDIF
      !
      ! This is needed to build the hash tables because they are
      ! locally indexed.
      CALL dbcsr_reset_locals(product_matrix)
      !
      IF (debug_mod) THEN
         WRITE (*, *) routineN//" Matrices ", dbcsr_get_matrix_type(matrix_a), &
            dbcsr_get_matrix_type(matrix_b), dbcsr_get_matrix_type(matrix_c)
         WRITE (*, *) routineN//" Matrices ", transa_l, transb_l, "keep", keep_product_data
      ENDIF
      IF (keep_product_data) THEN
         IF (product_reindex) THEN
            IF (debug_mod) WRITE (*, *) routineN//" Making canonical index"
            CALL dbcsr_make_index_canonical(product_matrix)
         ENDIF
         IF (ASSOCIATED(product_matrix%m%wms)) &
            CPABORT("Product matrix should be finalized!")
         CALL dbcsr_make_untransposed_blocks(product_matrix)
!$OMP PARALLEL &
!$OMP DEFAULT (NONE) SHARED (product_matrix)
         ! For the multiply logic to work correctly, existing data must
         ! be added only after the index has been transformed into the
         ! canonical form.
         CALL dbcsr_add_wm_from_matrix(product_matrix)
!$OMP END PARALLEL
      ELSE
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ithread) &
!$OMP SHARED(product_matrix, memtype_product_wm)
         ithread = 0
!$       ithread = OMP_GET_THREAD_NUM()
         CALL dbcsr_work_create(product_matrix, work_mutable=.FALSE., &
                                memory_type=memtype_product_wm(ithread)%p)
!$OMP END PARALLEL
      ENDIF
      CALL dbcsr_data_set_size_referenced(product_matrix%m%data_area, 0)
      product_matrix%m%valid = .FALSE.
      !
      IF (my_use_mpi_exp) THEN
         !
         max_nblocks = product_matrix%m%nblks
         !
         ! Right buffer images
         CALL dbcsr_make_buffers(matrix_right, rdist_right, .FALSE., &
                                 f_k, l_k, f_col, l_col, &
                                 PRESENT(filter_eps), &
                                 has_acc, &
                                 alpha)
         !
         ! Left buffer images
         CALL dbcsr_make_buffers(matrix_left, rdist_left, .TRUE., &
                                 f_row, l_row, f_k, l_k, &
                                 PRESENT(filter_eps), &
                                 has_acc)
      ELSE
         product_matrix%m%nblks = 0
         product_matrix%m%nze = 0
         product_matrix%m%row_p(:) = 0
         !
         ! Right images
         CALL make_m2s(matrix_right, m2s_right, rdist_right, dense_rdist_right, &
                       use_dense_mult, ab_dense, use_local_indexing, use_list_indexing, &
                       "R", f_k, l_k, f_row, l_row, f_col, l_col, &
                       dense_k_sizes, dense_col_sizes, &
                       k_vmap, m_map, n_map, &
                       alpha)
         !
         ! Left images
         CALL make_m2s(matrix_left, m2s_left, rdist_left, dense_rdist_left, &
                       use_dense_mult, ab_dense, use_local_indexing, use_list_indexing, &
                       "L", f_k, l_k, f_row, l_row, f_col, l_col, &
                       dense_row_sizes, dense_k_sizes, &
                       k_vmap, m_map, n_map)
      ENDIF
      !
      IF (ab_dense) THEN
         CALL array_release(k_vmap)
         CALL array_release(dense_row_sizes)
         CALL array_release(dense_col_sizes)
         CALL array_release(dense_k_sizes)
      ENDIF
      !
      ! The limits were already used.  Reset them.
      f_row = 0; l_row = 0
      f_col = 0; l_col = 0
      f_k = 0; l_k = 0
      !
      my_flop = 0
      IF (my_use_mpi_exp) THEN
         CALL cannon_multiply_low_rma(rdist_left, rdist_right, &
                                      matrix_left, matrix_right, product_matrix, &
                                      retain_sparsity=retain_sparsity, &
                                      filter_eps=filter_eps, &
                                      flop=my_flop, &
                                      keep_product_data=keep_product_data)
      ELSE
         CALL cannon_multiply_low(m2s_left, m2s_right, product_matrix, &
                                  retain_sparsity=retain_sparsity, &
                                  filter_eps=filter_eps, &
                                  flop=my_flop)
         CALL dbcsr_finalize(product_matrix)
      ENDIF
      !
      IF (PRESENT(flop)) THEN
         ! return the average number of flops per MPI rank.
         ! Variance (which is fairly large) could be computed as well.
         CALL timeset(routineN//"_mpsum_flop", handle2)
         comm = dbcsr_mp_group(dbcsr_distribution_mp(dbcsr_distribution(matrix_c)))
         numnodes = dbcsr_mp_numnodes(dbcsr_distribution_mp(dbcsr_distribution(matrix_c)))
         CALL mp_sum(my_flop, comm)
         flop = (my_flop+numnodes-1)/numnodes
         CALL timestop(handle2)
      ENDIF
      !
      IF (new_left) CALL dbcsr_release(matrix_left)
      IF (new_right) CALL dbcsr_release(matrix_right)
      IF (release_tdist) THEN
         CALL dbcsr_distribution_no_threads(product_matrix%m%dist)
      ENDIF
      !
      CALL dbcsr_release_locals(product_matrix)
      ! The index of the product matrix is reset to the CP2K form if it
      ! was previously set to the canonical form.
      IF (product_reindex) THEN
         IF (debug_mod) WRITE (*, *) routineN//" Making CP2K index"
         CALL dbcsr_make_index_canonical(product_matrix, cp2k=.TRUE.)
      ENDIF
      IF (use_dense_mult) THEN
         CALL dbcsr_release(matrix_c)
         CALL dbcsr_init(matrix_c)
         CALL dbcsr_make_undense(product_matrix, matrix_c, &
                                 old_product_distribution, &
                                 old_product_row_blk_offsets, old_product_col_blk_offsets, &
                                 old_product_row_blk_sizes, old_product_col_blk_sizes, &
                                 m_map, n_map)
         CALL dbcsr_release(product_matrix)
         CALL array_release(old_product_row_blk_offsets)
         CALL array_release(old_product_col_blk_offsets)
         CALL array_release(old_product_row_blk_sizes)
         CALL array_release(old_product_col_blk_sizes)
         CALL dbcsr_distribution_release(old_product_distribution)
      ELSE
         CALL dbcsr_release(matrix_c)
         CALL dbcsr_init(matrix_c)
         CALL dbcsr_copy(matrix_c, product_matrix, shallow_data=.TRUE.)
         CALL dbcsr_release(product_matrix)
      ENDIF
      !
      IF (.NOT. my_use_mpi_exp) THEN
         CALL dbcsr_destroy_array(m2s_left)
         DEALLOCATE (m2s_left)
         CALL dbcsr_destroy_array(m2s_right)
         DEALLOCATE (m2s_right)
      ENDIF
      !
      CALL dbcsr_image_dist_release(rdist_left)
      CALL dbcsr_image_dist_release(rdist_right)
      IF (ab_dense) THEN
         CALL array_release(m_map)
         CALL array_release(n_map)
      ENDIF
      !
      ! if filtering is requested remove small blocks, unless the sparsity needs to be kept.
      !
      IF (PRESENT(filter_eps) .AND. .NOT. keep_sparsity) THEN
         eps_any = dbcsr_scalar(filter_eps)
         CALL dbcsr_scalar_fill_all(eps_any)
         CALL dbcsr_scalar_set_type(eps_any, dbcsr_get_data_type(matrix_c))
         CALL dbcsr_filter(matrix_c, eps_any, quick=.FALSE.)
      ENDIF
      !
      ! To support the canonical multiply (all non-transposed blocks),
      ! blocks may have to be transposed according to the CP2K
      ! triangular index.
      CALL dbcsr_make_untransposed_blocks(matrix_c)
      !
      IF (debug_mod .OR. careful_mod) THEN
         IF (debug_mod) &
            WRITE (*, *) routineN//" Use dense mult, symm", &
            use_dense_mult, ab_dense, dbcsr_has_symmetry(matrix_c)
         CALL dbcsr_verify_matrix(matrix_c)
         IF (debug_mod) THEN
            cs = dbcsr_checksum(matrix_c)
            WRITE (*, *) routineN//" Product checksum", cs
         ENDIF
      ENDIF

      ! This tends to trigger only when all of these conditions are fulfilled:
      !  - transa=="T"
      !  - matrix_c contains already blocks and beta is not zero
      !  - Cuda is enabled
      !  - multiple OpenMP threads are used
      IF (INT(matrix_c%m%nblks, KIND=int_8) > &
          INT(SIZE(array_data(matrix_c%m%row_blk_size)), KIND=int_8)* &
          INT(SIZE(array_data(matrix_c%m%col_blk_size)), KIND=int_8)) &
         CPABORT("Bug: Matrix contains too many blocks")
      output_unit = default_output_unit
      CALL timestop(handle)
   END SUBROUTINE dbcsr_mm_cannon_multiply

! **************************************************************************************************
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
! **************************************************************************************************
   SUBROUTINE make_m2s(matrix, m2s, rdist, dense_rdist, &
                       use_dense_mult, ab_dense, &
                       use_local_indexing, use_list_indexing, &
                       predistribute, f_k, l_k, f_row, l_row, f_col, l_col, &
                       row_blk_size, col_blk_size, &
                       k_vmap, m_map, n_map, &
                       alpha)
      TYPE(dbcsr_obj), INTENT(IN)                        :: matrix
      TYPE(dbcsr_2d_array_type), INTENT(OUT), POINTER    :: m2s
      TYPE(dbcsr_imagedistribution_obj), INTENT(INOUT)   :: rdist, dense_rdist
      LOGICAL, INTENT(IN)                                :: use_dense_mult, ab_dense, &
                                                            use_local_indexing, use_list_indexing
      CHARACTER, INTENT(IN)                              :: predistribute
      INTEGER, INTENT(IN)                                :: f_k, l_k, f_row, l_row, f_col, l_col
      TYPE(array_i1d_obj), INTENT(INOUT)                 :: row_blk_size, col_blk_size
      TYPE(array_i1d_obj), INTENT(IN)                    :: k_vmap, m_map, n_map
      TYPE(dbcsr_scalar_type), INTENT(IN), OPTIONAL      :: alpha

      INTEGER                                            :: i, im, j, jm
      INTEGER, DIMENSION(4)                              :: f_crop
      LOGICAL                                            :: do_scale, thread_redist
      TYPE(array_i1d_obj)                                :: col_map, row_map
      TYPE(dbcsr_obj)                                    :: dense_template, matrix_tmp

      ALLOCATE (m2s)
      do_scale = .FALSE.
      IF (PRESENT(alpha)) THEN
         IF (.NOT. dbcsr_scalar_are_equal(alpha, dbcsr_scalar_one(alpha%data_type))) THEN
            do_scale = .TRUE.
         END IF
      END IF

      IF (do_scale) THEN
         ! Copy and scale matrix if alpha is not 1.
         CALL dbcsr_make_images(matrix, m2s, rdist, &
                                predistribute=predistribute, &
                                no_copy_data=use_dense_mult, scale_value=alpha)
      ELSE
         CALL dbcsr_make_images(matrix, m2s, rdist, &
                                predistribute=predistribute, &
                                no_copy_data=use_dense_mult)
      END IF

      im = SIZE(m2s%mats, 1)
      jm = SIZE(m2s%mats, 2)
      SELECT CASE (predistribute)
      CASE ('L')
         f_crop = (/f_row, l_row, f_k, l_k/)
         row_map = m_map
         col_map = k_vmap
         thread_redist = .TRUE.
      CASE default
         f_crop = (/f_k, l_k, f_col, l_col/)
         row_map = k_vmap
         col_map = n_map
         thread_redist = .FALSE.
      END SELECT

      ! Post-processing of images.
      DO i = 1, im
         DO j = 1, jm
            CALL dbcsr_reset_vlocals(m2s%mats(i, j), rdist)
            ! Crop if necessary
            IF (ANY(f_crop .NE. 0)) THEN
               CALL dbcsr_init(matrix_tmp)
               CALL dbcsr_crop_matrix(matrix_tmp, m2s%mats(i, j), &
                                      full_row_bounds=f_crop(1:2), &
                                      full_column_bounds=f_crop(3:4), &
                                      shallow_data=.FALSE.)
               CALL dbcsr_release(m2s%mats(i, j))
               CALL dbcsr_copy(m2s%mats(i, j), matrix_tmp, shallow_data=.TRUE.)
               CALL dbcsr_release(matrix_tmp)
               CALL dbcsr_reset_vlocals(m2s%mats(i, j), rdist)
            ENDIF
         ENDDO
      ENDDO

      IF (ab_dense) THEN
         CALL dbcsr_init(dense_template)
         CALL dbcsr_create(dense_template, template=matrix, &
                           dist=dense_rdist%i%main, &
                           row_blk_size=array_data(row_blk_size), col_blk_size=array_data(col_blk_size))
         CALL dbcsr_make_images_dense(m2s, dense_rdist, &
                                      row_map=row_map, col_map=col_map, &
                                      join_cols=use_dense_mult, join_rows=ab_dense, &
                                      new_template=dense_template)
         CALL dbcsr_image_dist_release(rdist)
         rdist = dense_rdist
         CALL dbcsr_image_dist_hold(rdist)
         DO i = 1, im
            DO j = 1, jm
               CALL dbcsr_reset_vlocals(m2s%mats(i, j), rdist)
            ENDDO
         ENDDO
      ENDIF

      DO i = 1, im
         DO j = 1, jm
            ! Convert to local-row index
            IF (use_local_indexing) THEN
               CALL dbcsr_make_index_local_row(m2s%mats(i, j))
            ENDIF
            ! Convert to list index
            IF (use_list_indexing) THEN
               CALL dbcsr_make_index_list(m2s%mats(i, j), thread_redist=thread_redist)
            ENDIF
            IF (use_local_indexing .AND. .NOT. use_list_indexing) THEN
               CALL dbcsr_index_compact(m2s%mats(i, j))
            ENDIF
         ENDDO
      ENDDO

      IF (ab_dense) THEN
         CALL dbcsr_image_dist_release(dense_rdist)
         CALL dbcsr_release(dense_template)
      ENDIF

   END SUBROUTINE make_m2s

! **************************************************************************************************
!> \brief Prepare orig images
!> \param matrix ...
!> \param imgdist ...
!> \param is_left ...
!> \param f_row ...
!> \param l_row ...
!> \param f_col ...
!> \param l_col ...
!> \param otf_filtering ...
!> \param has_acc ...
!> \param alpha ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE dbcsr_make_buffers(matrix, imgdist, is_left, &
                                 f_row, l_row, f_col, l_col, &
                                 otf_filtering, has_acc, &
                                 alpha)
      TYPE(dbcsr_obj), INTENT(IN)                        :: matrix
      TYPE(dbcsr_imagedistribution_obj), INTENT(IN)      :: imgdist
      LOGICAL, INTENT(IN)                                :: is_left
      INTEGER, INTENT(IN)                                :: f_row, l_row, f_col, l_col
      LOGICAL, INTENT(IN)                                :: otf_filtering, has_acc
      TYPE(dbcsr_scalar_type), INTENT(IN), OPTIONAL      :: alpha

      LOGICAL                                            :: do_scale

      do_scale = .FALSE.
      IF (PRESENT(alpha)) THEN
         IF (.NOT. dbcsr_scalar_are_equal(alpha, dbcsr_scalar_one(alpha%data_type))) THEN
            do_scale = .TRUE.
         END IF
      END IF
      !
      IF (do_scale) THEN
         CALL make_buffers(matrix, imgdist, is_left, &
                           f_row, l_row, f_col, l_col, otf_filtering, has_acc, alpha)
      ELSE
         CALL make_buffers(matrix, imgdist, is_left, &
                           f_row, l_row, f_col, l_col, otf_filtering, has_acc)
      ENDIF
   END SUBROUTINE dbcsr_make_buffers

! **************************************************************************************************
!> \brief Prepare orig images
!> \param matrix ...
!> \param imgdist ...
!> \param is_left ...
!> \param f_row ...
!> \param l_row ...
!> \param f_col ...
!> \param l_col ...
!> \param otf_filtering ...
!> \param has_acc ...
!> \param scale_value ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE make_buffers(matrix, imgdist, is_left, &
                           f_row, l_row, f_col, l_col, &
                           otf_filtering, has_acc, &
                           scale_value)
      TYPE(dbcsr_obj), INTENT(IN)                        :: matrix
      TYPE(dbcsr_imagedistribution_obj), INTENT(IN)      :: imgdist
      LOGICAL, INTENT(IN)                                :: is_left
      INTEGER, INTENT(IN)                                :: f_row, l_row, f_col, l_col
      LOGICAL, INTENT(IN)                                :: otf_filtering, has_acc
      TYPE(dbcsr_scalar_type), INTENT(IN), OPTIONAL      :: scale_value

      CHARACTER(len=*), PARAMETER :: routineN = 'make_buffers', routineP = moduleN//':'//routineN

      INTEGER :: blk, blk_p, bp, col, col_img, col_size, data_size, data_size_diag, data_type, &
         dst_proc, g2l_img_ref, grp, handle, iproc, irequests, it, ithread, local_no_empty_images, &
         meta_size, meta_size_diag, meta_size_diag_scatter, mygrp, mygrp_scatter, mynode, mypcol, &
         myprow, myt, nblkcols_local, nblkrows_local, nblocks, ncols_images, nimages_merged, &
         nimages_unmerged, nprocs, nprocs_scatter, nrows_images, nsymmetries, nthreads, nze, pcol, &
         pdiag, prow, request_data_size_diag, request_data_size_scatter, request_meta_size_diag, &
         request_meta_size_scatter, request_no_empty_images, row, row_img
      INTEGER :: row_size, scatter_images, scatter_images_proc, size_index_unmerged, stored_col, &
         stored_row, sym_p, symmetry_i, total_no_empty_images, tr_col_size, tr_row_size
      INTEGER, ALLOCATABLE, DIMENSION(:) :: data_displ_scatter, data_size_scatter, g2l_map_cols, &
         g2l_map_cols_merge, g2l_map_rows, g2l_map_rows_merge, img_map, img_refs, img_refs_cols, &
         img_refs_rows, meta_displ_cluster_scatter, meta_displ_scatter, meta_size_cluster_scatter, &
         meta_size_scatter, tmp_img_offset
      INTEGER, ALLOCATABLE, DIMENSION(:), TARGET         :: meta_recv, meta_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: local_refs_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)        :: recv_refs
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :, :)     :: refs
      INTEGER, DIMENSION(2)                              :: recv_sym_size, sym_size
      INTEGER, DIMENSION(:), POINTER :: col_dist, col_img_dist, local_cols, local_refs_data_size, &
         local_rows, meta_buffer_p, no_empty_images, no_empty_images_displ, row_dist, &
         row_img_dist, threads_dist
      INTEGER, DIMENSION(:, :), POINTER                  :: blacs2mpi, local_refs_displ_unmerged, &
                                                            local_refs_meta_size
      INTEGER, DIMENSION(:, :, :), POINTER               :: local_refs_displ
      INTEGER, DIMENSION(dbcsr_slot_nblkrows_total:&
         dbcsr_slot_nfullcols_local)                     :: global_indices
      INTEGER, POINTER                                   :: coli, rowi
      INTEGER, TARGET                                    :: mi, ui
      LOGICAL                                            :: do_scatter, do_symmetry, is_diagonal, tr
      REAL(kind=sp), DIMENSION(:), POINTER               :: local_max_norms, max_norms
      TYPE(dbcsr_data_obj), POINTER                      :: data_buffer_p
      TYPE(dbcsr_data_obj), TARGET                       :: data_recv, data_send
      TYPE(dbcsr_distribution_obj)                       :: set_dist
      TYPE(dbcsr_iterator)                               :: iter
      TYPE(dbcsr_mp_obj)                                 :: mp_obj
      TYPE(dbcsr_rma_buffer), POINTER                    :: buffer
      TYPE(dbcsr_scalar_type)                            :: scale_neg_one
      TYPE(dbcsr_type)                                   :: sm

!$    INTEGER(kind=omp_lock_kind), ALLOCATABLE, DIMENSION(:) :: locks

      CALL timeset(routineN, handle)
      !
      ! Take input values and check validity
      IF (.NOT. dbcsr_valid_index(matrix)) &
         CPABORT("Matrix not initialized.")
      sm = matrix%m
      data_type = sm%data_type
      IF (data_type .NE. dbcsr_type_real_8 .AND. &
          data_type .NE. dbcsr_type_real_4 .AND. &
          data_type .NE. dbcsr_type_complex_8 .AND. &
          data_type .NE. dbcsr_type_complex_4) &
         CPABORT("Invalid data type.")
      scale_neg_one = dbcsr_scalar_negative(dbcsr_scalar_one(data_type))
      set_dist = imgdist%i%main
      row_dist => dbcsr_distribution_row_dist(set_dist)
      col_dist => dbcsr_distribution_col_dist(set_dist)
      local_rows => dbcsr_distribution_local_rows(set_dist)
      local_cols => dbcsr_distribution_local_cols(set_dist)
      nblkrows_local = SIZE(local_rows)
      nblkcols_local = SIZE(local_cols)
      IF (sm%symmetry) THEN
         IF (SIZE(row_dist) .NE. SIZE(col_dist)) &
            CPWARN('Unequal row and column distributions for symmetric matrix.')
      ENDIF
      nrows_images = imgdist%i%row_decimation
      row_img_dist => array_data(imgdist%i%row_image)
      ncols_images = imgdist%i%col_decimation
      col_img_dist => array_data(imgdist%i%col_image)
      mp_obj = dbcsr_distribution_mp(set_dist)
      grp = dbcsr_mp_group(mp_obj)
      blacs2mpi => dbcsr_mp_pgrid(mp_obj)
      mynode = dbcsr_mp_mynode(mp_obj)
      myprow = dbcsr_mp_myprow(mp_obj)
      mypcol = dbcsr_mp_mypcol(mp_obj)
      IF (MAXVAL(row_dist) .GT. UBOUND(blacs2mpi, 1)) &
         CPABORT("Row distribution references unexistent processor rows")
      IF (MAXVAL(col_dist) .GT. UBOUND(blacs2mpi, 2)) &
         CPABORT("Col distribution references unexistent processor cols")
      IF (.NOT. dbcsr_mp_has_subgroups(mp_obj)) &
         CPABORT("RMA requires subcommunicators!")
      ! Check threads configuration
      NULLIFY (threads_dist)
!$    IF (.NOT. dbcsr_distribution_has_threads(dbcsr_distribution(matrix))) &
!$       CPABORT("Thread distribution not defined")
!$    threads_dist => array_data(dbcsr_distribution_thread_dist(dbcsr_distribution(matrix)))
      !
      ! Crop matrix
      IF (ANY((/f_row, l_row, f_col, l_col/) .NE. 0)) THEN
         IF (f_row .LT. 0) &
            CPABORT("Invalid first row bound.")
         IF (l_row .GT. dbcsr_nfullrows_total(matrix)) &
            CPABORT("Invalid last row bound.")
         IF (f_col .LT. 0) &
            CPABORT("Invalid first column bound.")
         IF (l_col .GT. dbcsr_nfullcols_total(matrix)) &
            CPABORT("Invalid last column bound.")
         !
         ! For the moment crop is not supported for RMA
         CPABORT("Crop not supported for RMA algo yet.")
         !
      ENDIF
      !
      IF (dbcsr_has_symmetry(matrix)) THEN
         nsymmetries = 2
         do_symmetry = .TRUE.
      ELSE
         nsymmetries = 1
         do_symmetry = .FALSE.
      ENDIF
      sym_p = blacs2mpi(mypcol, myprow)
      is_diagonal = sym_p .EQ. mynode
      !
      size_index_unmerged = dbcsr_slot_nblks
      !
      IF (is_left) THEN
         ! merging over rows
         nimages_merged = nrows_images
         nimages_unmerged = ncols_images
         nprocs = dbcsr_mp_npcols(mp_obj)
         ALLOCATE (left_refs_meta_size(myprow*nimages_merged:(myprow+1)*nimages_merged-1, &
                                       0:nprocs*nimages_unmerged-1))
         ALLOCATE (left_local_refs_meta_size(nimages_merged, nimages_unmerged))
         left_local_refs_meta_size(:, :) = 0
         local_refs_meta_size => left_local_refs_meta_size
         ALLOCATE (left_local_refs_data_size(nimages_merged*nimages_unmerged))
         local_refs_data_size => left_local_refs_data_size
         ALLOCATE (left_refs_displ_unmerged(idata_displ:imeta_displ, &
                                            0:nprocs*nimages_unmerged-1))
         ALLOCATE (left_local_refs_displ_unmerged(idata_displ:imeta_displ, nimages_unmerged))
         local_refs_displ_unmerged => left_local_refs_displ_unmerged
         mygrp = dbcsr_mp_my_row_group(mp_obj)
         requests(7) = mp_request_null
         pdiag = myprow
         irequests = 1
         !
         ! Count the maximum possible multiplies per row for on-the-fly filtering
         IF (otf_filtering) THEN
            ALLOCATE (left_total_row_counts(nblkrows_local))
            left_total_row_counts = 0
         ENDIF
         buffer => rma_buffers_orig%left
         ALLOCATE (left_no_empty_images(nprocs))
         no_empty_images => left_no_empty_images
         ALLOCATE (left_no_empty_images_displ(nprocs))
         no_empty_images_displ => left_no_empty_images_displ
         do_scatter = .FALSE.
      ELSE
         ! merging over cols
         nimages_merged = ncols_images
         nimages_unmerged = nrows_images
         nprocs = dbcsr_mp_nprows(mp_obj)
         nprocs_scatter = dbcsr_mp_npcols(mp_obj)
         ALLOCATE (right_refs_meta_size(mypcol*nimages_merged:(mypcol+1)*nimages_merged-1, &
                                        0:nprocs*nimages_unmerged-1))
         ALLOCATE (right_local_refs_meta_size(nimages_merged, nimages_unmerged))
         right_local_refs_meta_size(:, :) = 0
         local_refs_meta_size => right_local_refs_meta_size
         ALLOCATE (right_local_refs_data_size(nimages_merged*nimages_unmerged))
         local_refs_data_size => right_local_refs_data_size
         ALLOCATE (right_refs_displ_unmerged(idata_displ:imeta_displ, &
                                             0:nprocs*nimages_unmerged-1))
         ALLOCATE (right_local_refs_displ_unmerged(idata_displ:imeta_displ, nimages_unmerged))
         local_refs_displ_unmerged => right_local_refs_displ_unmerged
         mygrp = dbcsr_mp_my_col_group(mp_obj)
         mygrp_scatter = dbcsr_mp_my_row_group(mp_obj)
         requests(8) = mp_request_null
         pdiag = mypcol
         irequests = 2
         buffer => rma_buffers_orig%right
         ALLOCATE (right_no_empty_images(nprocs))
         no_empty_images => right_no_empty_images
         ALLOCATE (right_no_empty_images_displ(nprocs))
         no_empty_images_displ => right_no_empty_images_displ
         global_indices(dbcsr_slot_nblkrows_total:dbcsr_slot_nfullcols_local) = &
            (/ &
            sm%nblkrows_total, &
            sm%nblkcols_total, &
            sm%nfullrows_total, &
            sm%nfullcols_total, &
            0, 0, &
            sm%nfullrows_local, &
            sm%nfullcols_local/)
         ! Do not scatter data when there is a single cluster or single column proc
         do_scatter = nprocs_scatter .GT. 1 .AND. nimages_merged .GT. 1 .AND. .NOT. has_acc
         IF (do_scatter .AND. is_diagonal) THEN
            ALLOCATE (meta_displ_cluster_scatter(nimages_unmerged))
            ALLOCATE (meta_size_cluster_scatter(nimages_unmerged))
            meta_displ_cluster_scatter(1) = 0
         ENDIF
         local_meta_size_scatter = 0
         local_data_size_scatter = 0
      ENDIF
      !
      ALLOCATE (local_refs_size(idata_size:imeta_size, nimages_merged, nimages_unmerged))
      local_refs_size(:, :, :) = 0
      !
      ! Evaluate maps for global -> local indexing (g2l_map_rows, g2l_map_cols)
      ! Evaluate maps for local_image -> local indexing (img_map)
      ! Evaluate refs (count and offset) per each image
      ALLOCATE (g2l_map_rows(sm%nblkrows_total), g2l_map_cols(sm%nblkcols_total))
      ALLOCATE (g2l_map_rows_merge(sm%nblkrows_total), g2l_map_cols_merge(sm%nblkcols_total))
      ALLOCATE (img_refs_rows(nimages_unmerged), img_refs_cols(nimages_unmerged))
      !
      img_refs_rows(:) = 0
      g2l_map_rows(:) = 0
      g2l_map_rows_merge(:) = 0
      IF (nrows_images .GT. 1) THEN
         g2l_img_ref = 0
         DO row = 1, nblkrows_local
            row_img = row_img_dist(local_rows(row))
            mi = MOD(row_img-1, nrows_images)+1
            img_refs_rows(mi) = img_refs_rows(mi)+1
            g2l_map_rows(local_rows(row)) = img_refs_rows(mi)
            g2l_img_ref = g2l_img_ref+1
            g2l_map_rows_merge(local_rows(row)) = g2l_img_ref
         ENDDO
      ELSE
         DO row = 1, nblkrows_local
            g2l_map_rows_merge(local_rows(row)) = row
         ENDDO
         ! Symmetric case, swap row with col images
         IF (do_symmetry) THEN
            DO row = 1, nblkrows_local
               col_img = col_img_dist(local_rows(row))
               mi = MOD(col_img-1, ncols_images)+1
               img_refs_rows(mi) = img_refs_rows(mi)+1
               g2l_map_rows(local_rows(row)) = img_refs_rows(mi)
            ENDDO
         ELSE
            g2l_map_rows(:) = g2l_map_rows_merge(:)
         ENDIF
         img_refs_rows(1) = nblkrows_local
      ENDIF
      !
      img_refs_cols(:) = 0
      g2l_map_cols(:) = 0
      g2l_map_cols_merge(:) = 0
      IF (ncols_images .GT. 1) THEN
         g2l_img_ref = 0
         DO col = 1, nblkcols_local
            col_img = col_img_dist(local_cols(col))
            mi = MOD(col_img-1, ncols_images)+1
            img_refs_cols(mi) = img_refs_cols(mi)+1
            g2l_map_cols(local_cols(col)) = img_refs_cols(mi)
            g2l_img_ref = g2l_img_ref+1
            g2l_map_cols_merge(local_cols(col)) = g2l_img_ref
         ENDDO
      ELSE
         DO col = 1, nblkcols_local
            g2l_map_cols_merge(local_cols(col)) = col
         ENDDO
         ! Symmetric case, swap col with row images
         IF (do_symmetry) THEN
            DO col = 1, nblkcols_local
               row_img = row_img_dist(local_cols(col))
               mi = MOD(row_img-1, nrows_images)+1
               img_refs_cols(mi) = img_refs_cols(mi)+1
               g2l_map_cols(local_cols(col)) = img_refs_cols(mi)
            ENDDO
         ELSE
            g2l_map_cols(:) = g2l_map_cols_merge(:)
         ENDIF
         img_refs_cols(1) = nblkcols_local
      ENDIF
      !
      IF (is_left) THEN
         ALLOCATE (img_refs(ncols_images))
         img_refs(1) = 0
         DO col_img = 1, ncols_images-1
            img_refs(col_img+1) = img_refs(col_img)+img_refs_cols(col_img)
         ENDDO
         ALLOCATE (img_map(nblkcols_local))
         IF (ncols_images .GT. 1) THEN
            ALLOCATE (tmp_img_offset(ncols_images))
            tmp_img_offset(:) = img_refs(:)
            DO col = 1, nblkcols_local
               col_img = col_img_dist(local_cols(col))
               ui = MOD(col_img-1, ncols_images)+1
               tmp_img_offset(ui) = tmp_img_offset(ui)+1
               img_map(tmp_img_offset(ui)) = col
            ENDDO
         ELSE
            DO col = 1, nblkcols_local
               img_map(col) = col
            ENDDO
         ENDIF
      ELSE
         ALLOCATE (img_refs(nrows_images))
         img_refs(1) = 0
         DO row_img = 1, nrows_images-1
            img_refs(row_img+1) = img_refs(row_img)+img_refs_rows(row_img)
         ENDDO
         ALLOCATE (img_map(nblkrows_local))
         IF (nrows_images .GT. 1) THEN
            ALLOCATE (tmp_img_offset(nrows_images))
            tmp_img_offset(:) = img_refs(:)
            DO row = 1, nblkrows_local
               row_img = row_img_dist(local_rows(row))
               mi = MOD(row_img-1, nrows_images)+1
               tmp_img_offset(mi) = tmp_img_offset(mi)+1
               img_map(tmp_img_offset(mi)) = row
            ENDDO
         ELSE
            DO row = 1, nblkrows_local
               img_map(row) = row
            ENDDO
         ENDIF
      ENDIF
      IF (ALLOCATED(tmp_img_offset)) DEALLOCATE (tmp_img_offset)
      !
      ALLOCATE (local_refs_displ(idata_displ:imeta_displ, nimages_merged, nimages_unmerged))
      data_size = 0
      meta_size = 0
      data_size_diag = 0
      meta_size_diag = 0
      meta_size_diag_scatter = 0
      local_no_empty_images = 0
      !
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP PRIVATE (ithread,myt,iter,row,col,blk,row_size,col_size,&
!$OMP          stored_row,stored_col,blk_p,bp,tr,&
!$OMP          nze,symmetry_i,row_img,col_img,rowi,coli,&
!$OMP          tr_row_size,tr_col_size,prow,pcol,dst_proc,&
!$OMP          data_buffer_p,meta_buffer_p,&
!$OMP          mi,ui,it,nblocks) &
!$OMP SHARED (nthreads,refs,matrix,nsymmetries,do_symmetry,&
!$OMP         row_img_dist,col_img_dist,imgdist,row_dist,col_dist,&
!$OMP         is_left,local_refs_size,nimages_merged,nimages_unmerged,&
!$OMP         data_size,meta_size,data_type,memtype_mpi_buffer,sm,&
!$OMP         img_refs_cols,img_refs_rows,img_refs,img_map,&
!$OMP         g2l_map_cols,g2l_map_rows,recv_refs,grp,&
!$OMP         scale_value,scale_neg_one,data_send,data_recv,&
!$OMP         size_index_unmerged,sym_size,recv_sym_size,meta_send,&
!$OMP         mp_obj,threads_dist,meta_recv,nrows_images,ncols_images,&
!$OMP         locks,mynode,blacs2mpi,myprow,mypcol,sym_p,mygrp,&
!$OMP         left_refs_meta_size,right_refs_meta_size,&
!$OMP         left_refs_data_size,right_refs_data_size,&
!$OMP         local_refs_meta_size,local_refs_data_size,&
!$OMP         local_refs_displ,nprocs_scatter,meta_scatter,&
!$OMP         left_refs_displ_unmerged,right_refs_displ_unmerged,&
!$OMP         local_refs_displ_unmerged,requests,is_diagonal,&
!$OMP         g2l_map_cols_merge,g2l_map_rows_merge,&
!$OMP         meta_size_diag,data_size_diag,pdiag,buffer,&
!$OMP         left_total_row_counts,otf_filtering,meta_size_diag_scatter,&
!$OMP         request_meta_size_diag,request_data_size_diag,&
!$OMP         no_empty_images,local_no_empty_images,&
!$OMP         request_no_empty_images,irequests,&
!$OMP         meta_size_scatter,data_size_scatter,iproc,&
!$OMP         scatter_images_proc,scatter_images,&
!$OMP         request_meta_size_scatter,mygrp_scatter,do_scatter,&
!$OMP         request_data_size_scatter,meta_displ_scatter,data_displ_scatter,&
!$OMP         local_meta_size_scatter,local_data_size_scatter,&
!$OMP         meta_displ_cluster_scatter,meta_size_cluster_scatter,global_indices) &
!$OMP         REDUCTION(MAX : max_nblocks)
      ithread = 0
!$    ithread = omp_get_thread_num()
      myt = ithread
      IF (is_left) THEN
         rowi => mi
         coli => ui
      ELSE
         rowi => ui
         coli => mi
      ENDIF
!$OMP MASTER
      nthreads = 1
!$    nthreads = omp_get_num_threads()
      ALLOCATE (refs(idata_size:imeta_size, &
                     nimages_merged, nimages_unmerged, &
                     0:nthreads, nsymmetries))
      refs(:, :, :, :, :) = 0
!$    IF (is_left) THEN
!$       size_index_unmerged = size_index_unmerged+nthreads+1
!$    ENDIF
!$    IF (is_left .AND. do_symmetry) THEN
!$       ALLOCATE (locks(0:nthreads-1))
!$    ENDIF
!$OMP END MASTER
!$OMP BARRIER
!$    IF (is_left .AND. do_symmetry) THEN
!$       call omp_init_lock(locks(ithread))
!$    ENDIF
      !
      ! By default all data are considered already in local proc
      dst_proc = ilocal_proc
      !
      ! Take data and meta dimensions per each thread
      CALL dbcsr_iterator_start(iter, matrix, shared=.TRUE.)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, blk, &
                                        row_size=row_size, col_size=col_size)
         nze = row_size*col_size
         IF (nze .EQ. 0) CYCLE
         DO symmetry_i = 1, nsymmetries
            IF (symmetry_i .EQ. 1) THEN
               stored_row = row; stored_col = col
            ELSE
               IF (row .EQ. col) CYCLE
               stored_row = col; stored_col = row
            ENDIF
            row_img = row_img_dist(stored_row)
            col_img = col_img_dist(stored_col)
            CALL image_calculator(imgdist, &
                                  prow=prow, pcol=pcol, &
                                  rowi=rowi, coli=coli, &
                                  myprow=row_dist(stored_row), myrowi=row_img, &
                                  mypcol=col_dist(stored_col), mycoli=col_img, &
                                  shifting='0')
            IF (do_symmetry) THEN
               ! Avoid unnecessary copy of data for diagonal procs
               IF (blacs2mpi(prow, pcol) .EQ. mynode .OR. is_diagonal) THEN
                  dst_proc = ilocal_proc
               ELSE
                  ! Move data to remove symmetric proc
                  dst_proc = isym_proc
               ENDIF
            ENDIF
!$          IF (is_left .AND. do_symmetry) THEN
!$             myt = threads_dist(stored_row)
!$          ENDIF
!$OMP ATOMIC
            refs(imeta_size, mi, ui, myt+1, dst_proc) = &
               refs(imeta_size, mi, ui, myt+1, dst_proc)+3
!$OMP ATOMIC
            refs(idata_size, mi, ui, myt+1, dst_proc) = &
               refs(idata_size, mi, ui, myt+1, dst_proc)+nze
         ENDDO ! symmetry_i
      ENDDO
      CALL dbcsr_iterator_stop(iter)
      ! Avoid unnecessary copy if data is already local
      IF (do_symmetry .AND. .NOT. is_diagonal) THEN
!$OMP BARRIER
!$OMP MASTER
         ALLOCATE (recv_refs(idata_size:imeta_size, &
                             nimages_merged, nimages_unmerged, &
                             0:nthreads))
         ! exchange sizes with the symmetric proc
         CALL mp_sendrecv(refs(:, :, :, 1:nthreads, isym_proc), sym_p, &
                          recv_refs(:, :, :, 1:nthreads), sym_p, grp)
!$OMP END MASTER
      ENDIF
!$OMP BARRIER
!$OMP CRITICAL(local_refs_acc)
      !
      ! Store local data and meta dimensions
      local_refs_size(:, :, :) = local_refs_size(:, :, :)+refs(:, :, :, ithread+1, ilocal_proc)
      IF (do_symmetry .AND. .NOT. is_diagonal) THEN
         local_refs_size(:, :, :) = local_refs_size(:, :, :)+recv_refs(:, :, :, ithread+1)
      ENDIF
!$OMP END CRITICAL(local_refs_acc)
!$OMP BARRIER
!$OMP MASTER
      !
      ! Take the displacements, taking in account the images merging
      ! Note that local diplacements refer to position 0
      DO ui = 1, nimages_unmerged
         local_refs_displ_unmerged(imeta_displ, ui) = meta_size
         local_refs_displ_unmerged(idata_displ, ui) = data_size
         DO mi = 1, nimages_merged
            IF (is_diagonal .AND. (ui .EQ. mi .OR. nimages_merged .EQ. 1)) THEN
               local_refs_displ(imeta_displ, mi, ui) = meta_size_diag
               local_refs_displ(idata_displ, mi, ui) = data_size_diag
            ELSE
               local_refs_displ(imeta_displ, mi, ui) = meta_size
               local_refs_displ(idata_displ, mi, ui) = data_size
            ENDIF
            IF (local_refs_size(imeta_size, mi, ui) .EQ. 0) CYCLE
            local_no_empty_images = local_no_empty_images+1
            local_refs_data_size(local_no_empty_images) = local_refs_size(idata_size, mi, ui)
            ! Include stats slots
            local_refs_size(imeta_size, mi, ui) = local_refs_size(imeta_size, mi, ui)+ &
                                                  size_index_unmerged
            local_refs_meta_size(mi, ui) = local_refs_size(imeta_size, mi, ui)
            IF (is_diagonal .AND. (ui .EQ. mi .OR. nimages_merged .EQ. 1)) THEN
               meta_size_diag = meta_size_diag+local_refs_size(imeta_size, mi, ui)
               data_size_diag = data_size_diag+local_refs_size(idata_size, mi, ui)
            ELSE
               meta_size = meta_size+local_refs_size(imeta_size, mi, ui)
               data_size = data_size+local_refs_size(idata_size, mi, ui)
            ENDIF
         ENDDO
      ENDDO
      !
      ! Evaluate diagonal images to send for the right matrix
      IF (do_scatter) THEN
         IF (is_diagonal) THEN
            ALLOCATE (meta_size_scatter(0:nprocs_scatter-1), &
                      data_size_scatter(0:nprocs_scatter-1), &
                      meta_displ_scatter(0:nprocs_scatter-1), &
                      data_displ_scatter(0:nprocs_scatter-1))
            meta_size_scatter(:) = 0; data_size_scatter(:) = 0
            meta_displ_scatter(0) = 0; data_displ_scatter(0) = 0
            scatter_images_proc = nimages_unmerged/nprocs_scatter+1
            ui = MOD(nimages_unmerged, nprocs_scatter)
            mi = 0
            DO iproc = 0, nprocs_scatter-1
               IF (iproc .GT. 0) THEN
                  meta_displ_scatter(iproc) = meta_displ_scatter(iproc-1)+ &
                                              meta_size_scatter(iproc-1)
                  data_displ_scatter(iproc) = data_displ_scatter(iproc-1)+ &
                                              data_size_scatter(iproc-1)
               ENDIF
               IF (iproc .EQ. ui) THEN
                  scatter_images_proc = scatter_images_proc-1
               ENDIF
               scatter_images = 0
               DO WHILE (scatter_images .LT. scatter_images_proc)
                  mi = mi+1
                  scatter_images = scatter_images+1
                  meta_displ_cluster_scatter(mi) = meta_size_diag_scatter
                  IF (local_refs_size(imeta_size, mi, mi) .EQ. 0) CYCLE
                  IF (iproc .NE. pdiag) THEN
                     meta_size_cluster_scatter(mi) = local_refs_size(imeta_size, mi, mi)+ &
                                                     dbcsr_num_slots- &
                                                     size_index_unmerged
                  ELSE
                     ! For diagonal procs the empty data to scatter are considered empty
                     meta_size_cluster_scatter(mi) = dbcsr_num_slots
                  ENDIF
                  meta_size_diag_scatter = meta_size_diag_scatter+ &
                                           meta_size_cluster_scatter(mi)
                  meta_size_scatter(iproc) = &
                     meta_size_scatter(iproc)+meta_size_cluster_scatter(mi)
                  data_size_scatter(iproc) = &
                     data_size_scatter(iproc)+local_refs_size(idata_size, mi, mi)
               ENDDO
            ENDDO
            ! Skip diagonal proc
            data_size_scatter(pdiag) = 0
            CALL ensure_array_size(meta_scatter, ub=meta_size_diag_scatter, &
                                   nocopy=.TRUE., memory_type=memtype_mpi_buffer)
         ENDIF
         !
         ! Scatter the diagonal sizes
         CALL mp_iscatter(meta_size_scatter, local_meta_size_scatter, &
                          myprow, mygrp_scatter, request_meta_size_scatter)
         CALL mp_iscatter(data_size_scatter, local_data_size_scatter, &
                          myprow, mygrp_scatter, request_data_size_scatter)
      ENDIF
      !
      ! Bcast the diagonal sizes
      CALL mp_ibcast(meta_size_diag, pdiag, mygrp, request_meta_size_diag)
      CALL mp_ibcast(data_size_diag, pdiag, mygrp, request_data_size_diag)
      !
      ! Exchange refs
      IF (is_left) THEN
         CALL mp_iallgather(local_refs_meta_size, left_refs_meta_size, mygrp, requests(irequests))
         CALL mp_iallgather(local_refs_displ_unmerged, left_refs_displ_unmerged, mygrp, requests(2+irequests))
      ELSE
         CALL mp_iallgather(local_refs_meta_size, right_refs_meta_size, mygrp, requests(irequests))
         CALL mp_iallgather(local_refs_displ_unmerged, right_refs_displ_unmerged, mygrp, requests(2+irequests))
      ENDIF
      CALL mp_iallgather(local_no_empty_images, no_empty_images, mygrp, request_no_empty_images)
      !
      ! Allocate data and meta buffers
      CALL rma_buffer_init(buffer, data_type, &
                           data_size, meta_size, &
                           data_memory_type=memtype_mpi_buffer)
      IF (is_diagonal) THEN
         CALL resize_buffer_diag(buffer, data_type, &
                                 data_size_diag, meta_size_diag, &
                                 data_memory_type=memtype_mpi_buffer)
      ENDIF
      !
      ! Change displacements for the first thread
      ! In this case they refer to position 1
      refs(imeta_size, :, :, 0, ilocal_proc) = local_refs_displ(imeta_displ, :, :)+1
      refs(idata_size, :, :, 0, ilocal_proc) = local_refs_displ(idata_displ, :, :)+1
      !
      ! Evaluate data and meta diplacements per each thread
      DO it = 1, nthreads-1
         refs(:, :, :, it, ilocal_proc) = refs(:, :, :, it, ilocal_proc)+refs(:, :, :, it-1, ilocal_proc)
      ENDDO
      !
      ! Take offsets for symmetric case
      IF (do_symmetry .AND. .NOT. is_diagonal) THEN
         sym_size(:) = 1
         recv_sym_size(:) = 1
         DO ui = 1, nimages_unmerged
            DO mi = 1, nimages_merged
               ! Empty cluster
               IF (local_refs_size(imeta_size, mi, ui) .EQ. 0) THEN
                  DO it = 0, nthreads
                     recv_refs(:, mi, ui, it) = recv_sym_size(:)
                  ENDDO
               ELSE
                  ! Carry previous value
                  recv_refs(:, mi, ui, 0) = recv_sym_size(:)
                  DO it = 1, nthreads
                     recv_refs(:, mi, ui, it) = &
                        recv_refs(:, mi, ui, it)+recv_refs(:, mi, ui, it-1)
                     ! Add symetric events to the diplacements of local data
                     refs(:, mi, ui, it, ilocal_proc) = refs(:, mi, ui, it, ilocal_proc)+ &
                                                        recv_refs(:, mi, ui, it)-recv_sym_size(:)
                  ENDDO
                  recv_sym_size(:) = recv_refs(:, mi, ui, nthreads)
               ENDIF
               refs(:, mi, ui, 0, isym_proc) = sym_size(:)
               DO it = 1, nthreads
                  refs(:, mi, ui, it, isym_proc) = &
                     refs(:, mi, ui, it, isym_proc)+refs(:, mi, ui, it-1, isym_proc)
               ENDDO
               sym_size(:) = refs(:, mi, ui, nthreads, isym_proc)
            ENDDO
         ENDDO
         !
         ! Allocate data/meta to send
         CALL dbcsr_data_init(data_send)
         CALL dbcsr_data_new(data_send, data_type, sym_size(idata_size)-1)
         ALLOCATE (meta_send(sym_size(imeta_size)-1))
      ENDIF
!$OMP END MASTER
!$OMP BARRIER
      !
      data_buffer_p => buffer%data
      meta_buffer_p => buffer%meta
      !
      ! Copy data and meta in the buffers
      CALL dbcsr_iterator_start(iter, matrix, shared=.TRUE.)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, blk, blk_p=blk_p, &
                                        row_size=row_size, col_size=col_size)
         nze = row_size*col_size
         IF (nze .EQ. 0) CYCLE
         bp = ABS(blk_p)
         DO symmetry_i = 1, nsymmetries
            IF (symmetry_i .EQ. 1) THEN
               stored_row = row; stored_col = col; tr = blk_p .LT. 0
               tr_row_size = col_size; tr_col_size = row_size
            ELSE
               IF (row .EQ. col) CYCLE
               stored_row = col; stored_col = row; tr = blk_p .GT. 0
               tr_row_size = row_size; tr_col_size = col_size
            ENDIF
            row_img = row_img_dist(stored_row)
            col_img = col_img_dist(stored_col)
            CALL image_calculator(imgdist, &
                                  prow=prow, pcol=pcol, &
                                  rowi=rowi, coli=coli, &
                                  myprow=row_dist(stored_row), myrowi=row_img, &
                                  mypcol=col_dist(stored_col), mycoli=col_img, &
                                  shifting='0')
            IF (do_symmetry .OR. is_diagonal) THEN
               IF (is_diagonal) THEN
                  IF (ui .EQ. mi .OR. nimages_merged .EQ. 1) THEN
                     data_buffer_p => buffer%data_diag
                     meta_buffer_p => buffer%meta_diag
                  ELSE
                     data_buffer_p => buffer%data
                     meta_buffer_p => buffer%meta
                  ENDIF
               ELSEIF (blacs2mpi(prow, pcol) .EQ. mynode) THEN
                  dst_proc = ilocal_proc
                  data_buffer_p => buffer%data
                  meta_buffer_p => buffer%meta
               ELSE
                  ! Move data to remove symmetric proc
                  dst_proc = isym_proc
                  data_buffer_p => data_send
                  meta_buffer_p => meta_send
               ENDIF
            ENDIF
!$          IF (is_left .AND. do_symmetry) THEN
!$             myt = threads_dist(stored_row)
!$             call omp_set_lock(locks(myt))
!$          ENDIF
            IF (tr) THEN
               CALL dbcsr_block_transpose_aa(data_buffer_p, sm%data_area, tr_row_size, tr_col_size, &
                                             refs(idata_size, mi, ui, myt, dst_proc), bp, &
                                             scale_value)
               IF (sm%negate_real .AND. sm%negate_imaginary) THEN
                  CALL dbcsr_block_scale(data_buffer_p, scale=scale_neg_one, &
                                         row_size=nze, col_size=1, &
                                         lb=refs(idata_size, mi, ui, myt, dst_proc))
               ELSEIF (sm%negate_real) THEN
                  CALL dbcsr_block_real_neg(data_buffer_p, row_size=nze, col_size=1, &
                                            lb=refs(idata_size, mi, ui, myt, dst_proc))
               ELSEIF (sm%negate_imaginary) THEN
                  CALL dbcsr_block_conjg(data_buffer_p, row_size=nze, col_size=1, &
                                         lb=refs(idata_size, mi, ui, myt, dst_proc))
               ENDIF
            ELSE
               CALL dbcsr_block_copy_aa(data_buffer_p, sm%data_area, row_size, col_size, &
                                        refs(idata_size, mi, ui, myt, dst_proc), bp, &
                                        scale_value)
            ENDIF
            !
            ! Set meta data (local indexing)
            IF (do_symmetry .AND. blacs2mpi(prow, pcol) .NE. mynode) THEN
               ! For the symmetric case we need to swap rows and cols maps
               IF (is_left) THEN
                  meta_buffer_p(refs(imeta_size, mi, ui, myt, dst_proc)) = g2l_map_cols_merge(stored_row)
                  meta_buffer_p(refs(imeta_size, mi, ui, myt, dst_proc)+1) = g2l_map_rows(stored_col)
               ELSE
                  meta_buffer_p(refs(imeta_size, mi, ui, myt, dst_proc)) = g2l_map_cols(stored_row)
                  meta_buffer_p(refs(imeta_size, mi, ui, myt, dst_proc)+1) = g2l_map_rows_merge(stored_col)
               ENDIF
               ! The offset will be adjusted on the symmetric proc
               meta_buffer_p(refs(imeta_size, mi, ui, myt, dst_proc)+2) = &
                  refs(idata_size, mi, ui, myt, dst_proc)
            ELSE
               IF (is_left) THEN
                  meta_buffer_p(refs(imeta_size, mi, ui, myt, dst_proc)+size_index_unmerged) = &
                     g2l_map_rows_merge(stored_row)
                  meta_buffer_p(refs(imeta_size, mi, ui, myt, dst_proc)+size_index_unmerged+1) = &
                     g2l_map_cols(stored_col)
               ELSE
                  meta_buffer_p(refs(imeta_size, mi, ui, myt, dst_proc)+size_index_unmerged) = &
                     g2l_map_rows(stored_row)
                  meta_buffer_p(refs(imeta_size, mi, ui, myt, dst_proc)+size_index_unmerged+1) = &
                     g2l_map_cols_merge(stored_col)
               ENDIF
               meta_buffer_p(refs(imeta_size, mi, ui, myt, dst_proc)+size_index_unmerged+2) = &
                  refs(idata_size, mi, ui, myt, dst_proc)-local_refs_displ(idata_displ, mi, ui)
            ENDIF
            refs(imeta_size, mi, ui, myt, dst_proc) = refs(imeta_size, mi, ui, myt, dst_proc)+3
            refs(idata_size, mi, ui, myt, dst_proc) = refs(idata_size, mi, ui, myt, dst_proc)+nze
!$          IF (is_left .AND. do_symmetry) THEN
!$             call omp_unset_lock(locks(myt))
!$          ENDIF
         ENDDO
      ENDDO
      CALL dbcsr_iterator_stop(iter)
      !
      IF (do_symmetry .AND. .NOT. is_diagonal) THEN
!$OMP BARRIER
!$OMP MASTER
         ! Send/recv data to the symmetric proc
         CALL dbcsr_data_init(data_recv)
         CALL dbcsr_data_new(data_recv, data_type, recv_sym_size(idata_size)-1)
         CALL dbcsr_sendrecv_any(data_send, sym_p, data_recv, sym_p, grp)
         ALLOCATE (meta_recv(recv_sym_size(imeta_size)-1))
         CALL mp_sendrecv(meta_send, sym_p, meta_recv, sym_p, grp)
         ! Delete send data
         CALL dbcsr_data_release(data_send)
         DEALLOCATE (meta_send)
!$OMP END MASTER
      ENDIF
!$OMP BARRIER
!$    IF (is_left .AND. do_symmetry) THEN
!$       call omp_destroy_lock(locks(ithread))
!$    ENDIF
      !
      meta_buffer_p => buffer%meta
      !
      ! Set stats slots
!$OMP DO schedule(dynamic)
      DO ui = 1, nimages_unmerged
         IF (is_left) THEN
            coli => ui
         ELSE
            rowi => ui
         ENDIF
         nblocks = 0
         DO mi = 1, nimages_merged
            IF (local_refs_size(imeta_size, mi, ui) .EQ. 0) CYCLE
            !
            ! Insert symmetric data in local buffers
            IF (do_symmetry .AND. .NOT. is_diagonal) THEN
               IF (recv_refs(imeta_size, mi, ui, nthreads) .GT. 1) THEN
                  DO it = 0, nthreads-1
                     ! Copy meta, block by block
                     DO blk = recv_refs(imeta_size, mi, ui, it), recv_refs(imeta_size, mi, ui, it+1)-1, 3
                        buffer%meta(refs(imeta_size, mi, ui, it, ilocal_proc)+size_index_unmerged) = &
                           meta_recv(blk)
                        buffer%meta(refs(imeta_size, mi, ui, it, ilocal_proc)+size_index_unmerged+1) = &
                           meta_recv(blk+1)
                        buffer%meta(refs(imeta_size, mi, ui, it, ilocal_proc)+size_index_unmerged+2) = &
                           meta_recv(blk+2)-recv_refs(idata_size, mi, ui, it)+ &
                           refs(idata_size, mi, ui, it, ilocal_proc)-local_refs_displ(idata_displ, mi, ui)
                        refs(imeta_size, mi, ui, it, ilocal_proc) = &
                           refs(imeta_size, mi, ui, it, ilocal_proc)+3
                     ENDDO
                     ! Copy data
                     CALL dbcsr_data_set(buffer%data, refs(idata_size, mi, ui, it, ilocal_proc), &
                                         recv_refs(idata_size, mi, ui, it+1)-recv_refs(idata_size, mi, ui, it), &
                                         data_recv, recv_refs(idata_size, mi, ui, it))
                  ENDDO
               ENDIF
            ENDIF
            IF (is_diagonal) THEN
               IF (ui .EQ. mi .OR. nimages_merged .EQ. 1) THEN
                  meta_buffer_p => buffer%meta_diag
               ELSE
                  meta_buffer_p => buffer%meta
               ENDIF
            ENDIF
            meta_buffer_p(local_refs_displ(imeta_displ, mi, ui)+dbcsr_slot_size) = &
               local_refs_size(imeta_size, mi, ui)
            meta_buffer_p(local_refs_displ(imeta_displ, mi, ui)+dbcsr_slot_nblks) = &
               (local_refs_size(imeta_size, mi, ui)-size_index_unmerged)/3
            nblocks = nblocks+meta_buffer_p(local_refs_displ(imeta_displ, mi, ui)+dbcsr_slot_nblks)
            IF (is_left) THEN
               ! first thread
!$             meta_buffer_p(local_refs_displ(imeta_displ, mi, ui)+dbcsr_slot_nblks+1) = 0
!$             DO it = 0, nthreads-1
!$                meta_buffer_p(local_refs_displ(imeta_displ, mi, ui)+dbcsr_slot_nblks+it+2) = &
!$                   (refs(imeta_size, mi, ui, it, ilocal_proc)-local_refs_displ(imeta_displ, mi, ui)-1)/3
!$             ENDDO
               ! Count the maximum possible multiplies per row for on-the-fly filtering
               IF (otf_filtering) THEN
                  DO row = local_refs_displ(imeta_displ, mi, ui)+size_index_unmerged+1, &
                           local_refs_displ(imeta_displ, mi, ui)+ &
                           meta_buffer_p(local_refs_displ(imeta_displ, mi, ui)+dbcsr_slot_size), 3
!$OMP ATOMIC
                     left_total_row_counts(meta_buffer_p(row)) = &
                        left_total_row_counts(meta_buffer_p(row))+1
                  ENDDO
               ENDIF
            ELSEIF (do_scatter .AND. is_diagonal .AND. (ui .EQ. mi .OR. nimages_merged .EQ. 1)) THEN
               ! Prepare meta buffer to scatter
               CALL set_empty_meta_index( &
                  meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_slot_size: &
                               meta_displ_cluster_scatter(mi)+dbcsr_num_slots), &
                  global_indices, dbcsr_num_slots)
               meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_slot_size) = &
                  meta_size_cluster_scatter(mi)
               meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_slot_nblks) = &
                  meta_buffer_p(local_refs_displ(imeta_displ, mi, ui)+dbcsr_slot_nblks)
               meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_slot_nze) = &
                  local_refs_size(idata_size, mi, ui)
               CALL image_calculator(imgdist, &
                                     prow=meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_slot_home_prow), &
                                     pcol=meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_slot_home_pcol), &
                                     rowi=meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_slot_home_rowi), &
                                     coli=meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_slot_home_coli), &
                                     vprow=meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_slot_home_vprow), &
                                     myrowi=rowi, mycoli=coli, &
                                     shifting='0')
               ! Considering merging of the rows
               meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_slot_home_vpcol) = &
                  meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_slot_home_pcol)
               meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_num_slots) = &
                  meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_slot_size)
               IF (meta_size_cluster_scatter(mi) .GT. dbcsr_num_slots) THEN
                  meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_num_slots+1: &
                               meta_displ_cluster_scatter(mi)+ &
                               meta_scatter(meta_displ_cluster_scatter(mi)+dbcsr_slot_size)) = &
                     meta_buffer_p(local_refs_displ(imeta_displ, mi, ui)+size_index_unmerged+1: &
                                   local_refs_displ(imeta_displ, mi, ui)+ &
                                   meta_buffer_p(local_refs_displ(imeta_displ, mi, ui)+dbcsr_slot_size))
               ENDIF
            ENDIF
         ENDDO
         !
         max_nblocks = MAX(max_nblocks, nblocks)
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      DEALLOCATE (refs)
      !
      IF (is_left .AND. otf_filtering) THEN
         CALL mp_isum(left_total_row_counts, mygrp, request_count_rows)
      ENDIF
      !
      IF (nimages_merged .EQ. 1) THEN
         IF (is_diagonal) THEN
            meta_buffer_p => buffer%meta_diag
         ELSE
            meta_buffer_p => buffer%meta
         ENDIF
         CALL setup_rec_index(meta_buffer_p, img_refs_rows, img_refs_cols, &
                              local_refs_size(imeta_size, 1, :), local_refs_displ(imeta_displ, 1, :), &
                              size_index_unmerged, is_left)
      ENDIF
      !
      IF (do_symmetry .AND. .NOT. is_diagonal) THEN
         CALL dbcsr_data_release(data_recv)
         DEALLOCATE (meta_recv, recv_refs)
      ENDIF
      DEALLOCATE (g2l_map_rows, g2l_map_cols)
      DEALLOCATE (g2l_map_rows_merge, g2l_map_cols_merge)
!$    IF (is_left .AND. do_symmetry) THEN
!$       DEALLOCATE (locks)
!$    ENDIF
      !
      CALL mp_wait(request_no_empty_images)
      no_empty_images_displ(1) = 0
      DO iproc = 1, nprocs-1
         no_empty_images_displ(iproc+1) = no_empty_images_displ(iproc)+ &
                                          no_empty_images(iproc)
      ENDDO
      total_no_empty_images = no_empty_images_displ(nprocs)+ &
                              no_empty_images(nprocs)
      !
      IF (is_left) THEN
         ALLOCATE (left_refs_data_size(total_no_empty_images))
         CALL mp_iallgather(local_refs_data_size(1:local_no_empty_images), &
                            left_refs_data_size, no_empty_images, no_empty_images_displ, &
                            mygrp, requests(4+irequests))
      ELSE
         ALLOCATE (right_refs_data_size(total_no_empty_images))
         CALL mp_iallgather(local_refs_data_size(1:local_no_empty_images), &
                            right_refs_data_size, no_empty_images, no_empty_images_displ, &
                            mygrp, requests(4+irequests))
      ENDIF
      !
      ! Evaluate and exchange max norms for MPI filtering
      IF (otf_filtering .AND. use_mpi_filtering) THEN
         IF (is_left) THEN
            ALLOCATE (left_local_max_norms(local_no_empty_images))
            ALLOCATE (left_max_norms(total_no_empty_images))
            local_max_norms => left_local_max_norms
            max_norms => left_max_norms
         ELSE
            ALLOCATE (right_local_max_norms(local_no_empty_images))
            ALLOCATE (right_max_norms(total_no_empty_images))
            local_max_norms => right_local_max_norms
            max_norms => right_max_norms
         ENDIF
         CALL calculate_max_image_norms(buffer, &
                                        local_refs_size(imeta_size, :, :), &
                                        local_refs_displ, &
                                        img_map, &
                                        img_refs, &
                                        dbcsr_row_block_sizes(matrix), &
                                        dbcsr_col_block_sizes(matrix), &
                                        local_rows, &
                                        local_cols, &
                                        size_index_unmerged+1, &
                                        local_max_norms, is_left, &
                                        .NOT. is_diagonal)
         CALL mp_iallgather(local_max_norms, max_norms, &
                            no_empty_images, no_empty_images_displ, &
                            mygrp, requests(6+irequests))
      ENDIF
      !
      DEALLOCATE (img_refs_cols)
      DEALLOCATE (img_refs)
      DEALLOCATE (img_map)
      DEALLOCATE (local_refs_displ)
      DEALLOCATE (img_refs_rows)
      DEALLOCATE (local_refs_size)
      !
      CALL mp_wait(request_meta_size_diag)
      CALL mp_wait(request_data_size_diag)
      IF (.NOT. is_diagonal) THEN
         CALL resize_buffer_diag(buffer, data_type, &
                                 data_size_diag, meta_size_diag, &
                                 data_memory_type=memtype_mpi_buffer)
      ENDIF
      !
      ! Bcast the diagonal data
      CALL mp_ibcast(buffer%meta_diag, pdiag, mygrp, requests_diag(irequests))
      CALL dbcsr_ibcast_any(buffer%data_diag, pdiag, mygrp, requests_diag(2+irequests))
      !
      ! Scatter the diagonal data
      IF (do_scatter) THEN
         CALL mp_wait(request_meta_size_scatter)
         CALL mp_wait(request_data_size_scatter)
         CALL ensure_array_size(local_meta_scatter, ub=local_meta_size_scatter, &
                                nocopy=.TRUE., memory_type=memtype_mpi_buffer)
         IF (dbcsr_data_valid(local_data_scatter)) THEN
            CALL dbcsr_data_ensure_size(local_data_scatter, local_data_size_scatter, &
                                        nocopy=.TRUE.)
         ELSE
            CALL dbcsr_data_init(local_data_scatter)
            CALL dbcsr_data_new(local_data_scatter, data_type=data_type, &
                                data_size=local_data_size_scatter, memory_type=memtype_mpi_buffer)
            CALL dbcsr_data_set_size_referenced(local_data_scatter, local_data_size_scatter)
         ENDIF
         !
         CALL mp_iscatter(meta_scatter, meta_size_scatter, meta_displ_scatter, &
                          local_meta_scatter, local_meta_size_scatter, &
                          myprow, mygrp_scatter, requests_scatter(1))
         CALL dbcsr_iscatterv_any(buffer%data_diag, data_size_scatter, data_displ_scatter, &
                                  local_data_scatter, local_data_size_scatter, &
                                  myprow, mygrp_scatter, requests_scatter(2))
         !
         IF (is_diagonal) THEN
            DEALLOCATE (meta_displ_cluster_scatter, meta_size_cluster_scatter)
            DEALLOCATE (meta_size_scatter, data_size_scatter, &
                        meta_displ_scatter, data_displ_scatter)
         ENDIF
      ENDIF
      !
      CALL timestop(handle)
   END SUBROUTINE make_buffers

! **************************************************************************************************
!> \brief Creates row and column images of a matrix.
!> \param[in] source          input matrix
!> \param[in,out] normalized  image array of the normalized matrix
!> \param[in] target_image_dist          normalize to this image distribution
!> \param[in] predistribute   (optional) predistribute data for multiplication
!> \param[in] no_copy_data    (optional) try to not merge data at the end
!> \param[in] scale_value     (optional) scale with this value
! **************************************************************************************************
   SUBROUTINE dbcsr_make_images(source, normalized, target_image_dist, &
                                predistribute, no_copy_data, scale_value)
      TYPE(dbcsr_obj), INTENT(IN)                        :: source
      TYPE(dbcsr_2d_array_type), INTENT(OUT)             :: normalized
      TYPE(dbcsr_imagedistribution_obj), INTENT(IN)      :: target_image_dist
      CHARACTER, INTENT(IN), OPTIONAL                    :: predistribute
      LOGICAL, INTENT(IN), OPTIONAL                      :: no_copy_data
      TYPE(dbcsr_scalar_type), INTENT(IN), OPTIONAL      :: scale_value

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_images', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

!   ---------------------------------------------------------------------------

      CALL timeset(routineN, handle)
      IF (my_use_mpi_exp) &
         CPABORT("RMA algo not supported here!")
      IF (.NOT. dbcsr_valid_index(source)) &
         CPABORT("Matrix not initialized.")
      CALL make_images(source, normalized, &
                       target_image_dist, desymmetrize=dbcsr_has_symmetry(source), &
                       predistribute=predistribute, &
                       no_copy_data=no_copy_data, &
                       scale_value=scale_value)
      normalized%image_dist = target_image_dist
      CALL dbcsr_image_dist_hold(normalized%image_dist)
      CALL timestop(handle)
   END SUBROUTINE dbcsr_make_images

! **************************************************************************************************
!> \brief Makes column-based and row-based images of a matrix.
!> \param[in] ism input symmetric matrix
!> \param[in,out] ums         normalied matrices
!> \param[in] target_imgdist  image distribution to normalize to
!> \param[in] desymmetrize    (optional) desymmetrize a symmetric matrix
!> \param[in] predistribute   (optional) predistribute data for multiplication
!> \param[in] no_copy_data    (optional) try to not merge data at the end
!> \param[in] scale_value     (optional) scale with this value
! **************************************************************************************************
   SUBROUTINE make_images(ism, ums, target_imgdist, desymmetrize, predistribute, &
                          no_copy_data, scale_value)
      TYPE(dbcsr_obj), INTENT(IN)                        :: ism
      TYPE(dbcsr_2d_array_type), INTENT(OUT)             :: ums
      TYPE(dbcsr_imagedistribution_obj), INTENT(IN)      :: target_imgdist
      LOGICAL, INTENT(IN), OPTIONAL                      :: desymmetrize
      CHARACTER, INTENT(IN), OPTIONAL                    :: predistribute
      LOGICAL, INTENT(IN), OPTIONAL                      :: no_copy_data
      TYPE(dbcsr_scalar_type), INTENT(IN), OPTIONAL      :: scale_value

      CHARACTER(len=*), PARAMETER :: routineN = 'make_images', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: metalen = 5
      LOGICAL, PARAMETER                                 :: dbg = .FALSE.

      CHARACTER                                          :: predist_type, predist_type_fwd
      INTEGER :: blk, blk_l, blk_p, bp, col, col_img, col_size, coli, data_p, data_type, dst_p, &
         handle, handle2, ithread, mp_group, ncol_images, nrow_images, nsymmetries, nthreads, &
         numproc, nze, pcol, prev_blk_p, prow, row, row_img, row_size, rowi, sd_pos, sm_pos, &
         src_p, stored_blk_p, stored_col, stored_row, symmetry_i, tr_col_size, tr_row_size, vcol, &
         vrow
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: myt_sdp, myt_smp, rd_disp, recv_meta, &
                                                            rm_disp, sd_disp, sdp, send_meta, &
                                                            sm_disp, smp
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: all_total_send_offset, blk_ps, blks, &
                                                            myt_total_send_count, &
                                                            total_recv_count, total_send_count
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)        :: myt_send_count, recv_count, send_count
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, col_dist, col_img_dist, &
                                                            row_blk_size, row_dist, row_img_dist
      INTEGER, DIMENSION(:, :), POINTER                  :: blacs2mpi
      LOGICAL                                            :: nocopy, tr
      TYPE(dbcsr_data_obj)                               :: received_data_area, recv_data_area, &
                                                            send_data_area
      TYPE(dbcsr_distribution_obj)                       :: old_dist, target_dist
      TYPE(dbcsr_iterator)                               :: iter
      TYPE(dbcsr_memtype_type)                           :: data_memory_type
      TYPE(dbcsr_mp_obj)                                 :: mp_obj
      TYPE(dbcsr_scalar_type)                            :: scale_neg_one
      TYPE(dbcsr_type)                                   :: sm

!   ---------------------------------------------------------------------------
! Check input matrix
! Set convenient variables to access input matrix info
!

      CALL timeset(routineN, handle)
      nocopy = .FALSE.
      IF (PRESENT(no_copy_data)) nocopy = no_copy_data
      sm = ism%m
      nsymmetries = 1
      IF (PRESENT(desymmetrize)) THEN
         IF (desymmetrize .AND. sm%symmetry) THEN
            nsymmetries = 2
         ENDIF
      ENDIF
      SELECT CASE (predistribute)
      CASE ('L', 'l')
         predist_type = 'L'
         predist_type_fwd = 'l'
      CASE ('R', 'r')
         predist_type = 'R'
         predist_type_fwd = 'r'
      CASE default
         CPABORT("Incorrect pre-shift specifier.")
      END SELECT
      data_type = sm%data_type
      IF (data_type .NE. dbcsr_type_real_8 .AND. &
          data_type .NE. dbcsr_type_real_4 .AND. &
          data_type .NE. dbcsr_type_complex_8 .AND. &
          data_type .NE. dbcsr_type_complex_4) &
         CPABORT("Invalid data type.")
      scale_neg_one = dbcsr_scalar_negative(dbcsr_scalar_one(data_type))
      row_blk_size => array_data(sm%row_blk_size)
      col_blk_size => array_data(sm%col_blk_size)
      old_dist = dbcsr_distribution(ism)
      target_dist = target_imgdist%i%main
      row_dist => dbcsr_distribution_row_dist(target_dist)
      col_dist => dbcsr_distribution_col_dist(target_dist)
      IF (sm%symmetry) THEN
         IF (SIZE(row_dist) .NE. SIZE(col_dist)) &
            CPWARN('Unequal row and column distributions for symmetric matrix.')
      ENDIF
      nrow_images = target_imgdist%i%row_decimation
      IF (nrow_images .GT. 1) THEN
         row_img_dist => array_data(target_imgdist%i%row_image)
      ELSE
         NULLIFY (row_img_dist)
      ENDIF
      ncol_images = target_imgdist%i%col_decimation
      IF (ncol_images .GT. 1) THEN
         col_img_dist => array_data(target_imgdist%i%col_image)
      ELSE
         NULLIFY (col_img_dist)
      ENDIF
      mp_obj = dbcsr_distribution_mp(target_dist)
      blacs2mpi => dbcsr_mp_pgrid(mp_obj)
      numproc = dbcsr_mp_numnodes(mp_obj)
      mp_group = dbcsr_mp_group(mp_obj)
      IF (MAXVAL(row_dist) .GT. UBOUND(blacs2mpi, 1)) &
         CPABORT('Row distribution references unexistent processor rows')
      IF (dbg) THEN
         IF (MAXVAL(row_dist) .NE. UBOUND(blacs2mpi, 1)) &
            CPWARN('Range of row distribution not equal to processor rows')
      ENDIF
      IF (MAXVAL(col_dist) .GT. UBOUND(blacs2mpi, 2)) &
         CPABORT('Col distribution references unexistent processor cols')
      IF (dbg) THEN
         IF (MAXVAL(col_dist) .NE. UBOUND(blacs2mpi, 2)) &
            CPWARN('Range of col distribution not equal to processor cols')
      ENDIF

      ! Check threads configuration
!$    IF (.NOT. dbcsr_distribution_has_threads(old_dist)) &
!$       CPABORT("Thread distribution not defined")

      ! Allocate shared temporary buffers
      !
      ALLOCATE (send_count(2, nrow_images, ncol_images, 0:numproc-1)); send_count = 0
      ALLOCATE (recv_count(2, nrow_images, ncol_images, 0:numproc-1))
      ALLOCATE (total_send_count(2, 0:numproc-1)); total_send_count = 0
      ALLOCATE (total_recv_count(2, 0:numproc-1))
      ALLOCATE (sdp(0:numproc-1))
      ALLOCATE (smp(0:numproc-1))
      ALLOCATE (sd_disp(0:numproc-1)); sd_disp(0) = 1
      ALLOCATE (sm_disp(0:numproc-1)); sm_disp(0) = 1
      ALLOCATE (rd_disp(0:numproc-1)); rd_disp(0) = 1
      ALLOCATE (rm_disp(0:numproc-1)); rm_disp(0) = 1
      ALLOCATE (all_total_send_offset(2, 0:numproc-1))
      ALLOCATE (blk_ps(nrow_images, ncol_images)); blk_ps = 1
      ALLOCATE (blks(nrow_images, ncol_images)); blks = 1
      !
      ! Allocate and init mats
      !
      ALLOCATE (ums%mats(nrow_images, ncol_images))
      data_memory_type = memtype_abpanel_1
      DO row_img = 1, nrow_images
         DO col_img = 1, ncol_images
            CALL dbcsr_init(ums%mats(row_img, col_img))
            CALL dbcsr_create(ums%mats(row_img, col_img), "imaged "//sm%name, &
                              target_dist, &
                              dbcsr_type_no_symmetry, array_data(sm%row_blk_size), array_data(sm%col_blk_size), &
                              sm%row_blk_size, sm%col_blk_size, &
                              nze=0, data_type=data_type, &
                              max_rbs=sm%max_rbs, max_cbs=sm%max_cbs, &
                              row_blk_offset=sm%row_blk_offset, col_blk_offset=sm%col_blk_offset, &
                              thread_dist=sm%dist, &
                              data_memory_type=data_memory_type, &
                              index_memory_type=memtype_mpi_buffer)
            ums%mats(row_img, col_img)%m%negate_real = sm%negate_real
            ums%mats(row_img, col_img)%m%negate_imaginary = sm%negate_imaginary
         ENDDO
      ENDDO

      nthreads = 1
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP PRIVATE (ithread, symmetry_i, row_img, col_img, &
!$OMP          myt_send_count, myt_total_send_count, &
!$OMP          iter, row, col, blk, row_size, col_size, stored_row, stored_col, &
!$OMP          prow, pcol, rowi, coli, vrow, vcol, dst_p, nze, myt_smp, myt_sdp, &
!$OMP          blk_p, bp, sd_pos, sm_pos,tr, &
!$OMP          tr_row_size, tr_col_size) &
!$OMP SHARED (nthreads, nsymmetries, row_img_dist, col_img_dist, &
!$OMP         nrow_images, ncol_images, numproc, scale_value, &
!$OMP         ums, sm, ism, target_imgdist, row_dist, col_dist,&
!$OMP         predist_type_fwd, blacs2mpi, row_blk_size, col_blk_size, &
!$OMP         send_count, recv_count, handle2,mp_group, &
!$OMP         total_send_count, total_recv_count, recv_data_area, nocopy, &
!$OMP         data_type, recv_meta, send_data_area, send_meta, &
!$OMP         sd_disp, sm_disp, rd_disp, rm_disp, all_total_send_offset, blk_ps, blks, &
!$OMP         received_data_area, scale_neg_one, memtype_abpanel_1)
      ithread = 0
!$    ithread = omp_get_thread_num()
!$OMP MASTER
!$    nthreads = omp_get_num_threads()
!$OMP END MASTER

      ! Allocate thread private data
      !
      ALLOCATE (myt_send_count(2, nrow_images, ncol_images, 0:numproc-1)); myt_send_count(:, :, :, :) = 0
      ALLOCATE (myt_total_send_count(2, 0:numproc-1))
      ! Thread-local pointers of the current adding position into the send buffers
      ALLOCATE (myt_smp(0:numproc-1), myt_sdp(0:numproc-1))

      ! Count sizes for sending
      !
      CALL dbcsr_iterator_start(iter, ism, shared=.TRUE.)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, blk, &
                                        row_size=row_size, col_size=col_size)
         nze = row_size*col_size
         IF (nze .EQ. 0) CYCLE
         DO symmetry_i = 1, nsymmetries
            IF (symmetry_i .EQ. 1) THEN
               stored_row = row; stored_col = col
            ELSE
               IF (row .EQ. col) CYCLE
               stored_row = col; stored_col = row
            ENDIF
            ! Where do we send this block?
            row_img = 1
            IF (nrow_images .GT. 1) row_img = row_img_dist(stored_row)
            col_img = 1
            IF (ncol_images .GT. 1) col_img = col_img_dist(stored_col)
            CALL image_calculator(target_imgdist, &
                                  prow=prow, rowi=rowi, &
                                  pcol=pcol, coli=coli, &
                                  vprow=vrow, vpcol=vcol, &
                                  myprow=row_dist(stored_row), myrowi=row_img, &
                                  mypcol=col_dist(stored_col), mycoli=col_img, &
                                  shifting=predist_type_fwd)
            dst_p = blacs2mpi(prow, pcol)
            ! These counts are meant for the thread that processes this row.
            myt_send_count(1, rowi, coli, dst_p) = &
               myt_send_count(1, rowi, coli, dst_p)+1
            myt_send_count(2, rowi, coli, dst_p) = &
               myt_send_count(2, rowi, coli, dst_p)+nze
         END DO ! symmetry_i
      END DO
      CALL dbcsr_iterator_stop(iter)
      DO dst_p = 0, numproc-1
         myt_total_send_count(1, dst_p) = SUM(myt_send_count(1, :, :, dst_p))
         myt_total_send_count(2, dst_p) = SUM(myt_send_count(2, :, :, dst_p))
      END DO
      ! Merge the send counts
!$OMP CRITICAL
      send_count(:, :, :, :) = send_count(:, :, :, :)+myt_send_count(:, :, :, :)
      total_send_count(:, :) = total_send_count(:, :)+myt_total_send_count(:, :)
!$OMP END CRITICAL
!$OMP BARRIER
!$OMP MASTER
      CALL timeset(routineN//"_sizes", handle2)
      CALL mp_alltoall(send_count, recv_count, 2*nrow_images*ncol_images, &
                       mp_group)
      CALL timestop(handle2)
!$OMP END MASTER
!$OMP BARRIER
      ! Fill in the meta data structures and copy the data.
!$OMP DO
      DO dst_p = 0, numproc-1
         total_recv_count(1, dst_p) = SUM(recv_count(1, :, :, dst_p))
         total_recv_count(2, dst_p) = SUM(recv_count(2, :, :, dst_p))
      ENDDO
!$OMP MASTER
      ! Allocate data structures needed for data exchange.
      CALL dbcsr_data_init(recv_data_area)
      IF (nrow_images .EQ. 1 .AND. ncol_images .EQ. 1 .OR. nocopy) THEN
         ! For some cases the faster dbcsr_special_finalize(reshuffle=.FALSE.) can be used.
         ! This basically makes this working matrix the actual data-area.
         ! Hence, for those cases we have to use data_memory_type already here.
         CALL dbcsr_data_new(recv_data_area, data_type, SUM(total_recv_count(2, :)), memory_type=memtype_abpanel_1)
      ELSE
         CALL dbcsr_data_new(recv_data_area, data_type, SUM(total_recv_count(2, :)))
      END IF
      ALLOCATE (recv_meta(metalen*SUM(total_recv_count(1, :))))
      CALL dbcsr_data_init(send_data_area)
      CALL dbcsr_data_new(send_data_area, data_type, SUM(total_send_count(2, :)))
      ALLOCATE (send_meta(metalen*SUM(total_send_count(1, :))))
      ! Calculate displacements for processors needed for the exchanges.
      DO dst_p = 1, numproc-1
         sm_disp(dst_p) = sm_disp(dst_p-1) &
                          +metalen*total_send_count(1, dst_p-1)
         sd_disp(dst_p) = sd_disp(dst_p-1) &
                          +total_send_count(2, dst_p-1)
         rm_disp(dst_p) = rm_disp(dst_p-1) &
                          +metalen*total_recv_count(1, dst_p-1)
         rd_disp(dst_p) = rd_disp(dst_p-1) &
                          +total_recv_count(2, dst_p-1)
      ENDDO
      myt_smp(:) = sm_disp(:)
      myt_sdp(:) = sd_disp(:)
      IF (nthreads .GT. 1) THEN
         all_total_send_offset(1, :) = myt_smp(:)+metalen*myt_total_send_count(1, :)
         all_total_send_offset(2, :) = myt_sdp(:)+myt_total_send_count(2, :)
      ENDIF
!$OMP END MASTER
!$OMP BARRIER
      IF (ithread .GT. 0) THEN
!$OMP CRITICAL
         myt_smp(:) = all_total_send_offset(1, :)
         myt_sdp(:) = all_total_send_offset(2, :)
         all_total_send_offset(1, :) &
            = all_total_send_offset(1, :)+metalen*myt_total_send_count(1, :)
         all_total_send_offset(2, :) &
            = all_total_send_offset(2, :)+myt_total_send_count(2, :)
!$OMP END CRITICAL
      ELSE
         CALL dbcsr_data_init(received_data_area)
         received_data_area = recv_data_area
         CALL dbcsr_data_hold(received_data_area)
         DO row_img = 1, nrow_images
            DO col_img = 1, ncol_images
               CALL dbcsr_work_create(ums%mats(row_img, col_img), &
                                      SUM(recv_count(1, row_img, col_img, :)), n=1)
               CALL dbcsr_data_hold(received_data_area)
               CALL dbcsr_data_release(ums%mats(row_img, col_img)%m%wms(1)%data_area)
               ums%mats(row_img, col_img)%m%wms(1)%data_area = received_data_area
            ENDDO
         ENDDO
      ENDIF
!$OMP BARRIER
!$OMP MASTER
      CALL m_memory()
!$OMP END MASTER

      ! Add timing call to the packing of the send buffers
      !
      CALL timeset(routineN//"_pack", handle2)
      ! Copies metadata and actual data to be sent into the send buffers.
      CALL dbcsr_iterator_start(iter, ism, shared=.TRUE.)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, blk, blk_p=blk_p, &
                                        row_size=row_size, col_size=col_size)
         nze = row_size*col_size
         IF (nze .EQ. 0) CYCLE
         bp = ABS(blk_p)
         DO symmetry_i = 1, nsymmetries
            IF (symmetry_i .EQ. 1) THEN
               stored_row = row; stored_col = col; tr = blk_p .LT. 0
               tr_row_size = col_size; tr_col_size = row_size
            ELSE
               IF (row .EQ. col) CYCLE
               stored_row = col; stored_col = row; tr = blk_p .GT. 0
               tr_row_size = row_size; tr_col_size = col_size
            ENDIF
            ! Where do we send this block?
            row_img = 1
            IF (nrow_images .GT. 1) row_img = row_img_dist(stored_row)
            col_img = 1
            IF (ncol_images .GT. 1) col_img = col_img_dist(stored_col)
            CALL image_calculator(target_imgdist, &
                                  prow=prow, rowi=rowi, &
                                  pcol=pcol, coli=coli, &
                                  vprow=vrow, vpcol=vcol, &
                                  myprow=row_dist(stored_row), myrowi=row_img, &
                                  mypcol=col_dist(stored_col), mycoli=col_img, &
                                  shifting=predist_type_fwd)
            dst_p = blacs2mpi(prow, pcol)
            sm_pos = myt_smp(dst_p)
            myt_smp(dst_p) = myt_smp(dst_p)+metalen
            sd_pos = myt_sdp(dst_p)
            myt_sdp(dst_p) = myt_sdp(dst_p)+nze
            IF (tr) THEN
               CALL dbcsr_block_transpose_aa(send_data_area, sm%data_area, tr_row_size, tr_col_size, &
                                             sd_pos, bp, scale_value)
               IF (sm%negate_real .AND. sm%negate_imaginary) THEN
                  CALL dbcsr_block_scale(send_data_area, scale=scale_neg_one, &
                                         row_size=nze, col_size=1, lb=sd_pos)
               ELSEIF (sm%negate_real) THEN
                  CALL dbcsr_block_real_neg(send_data_area, &
                                            row_size=nze, col_size=1, lb=sd_pos)
               ELSEIF (sm%negate_imaginary) THEN
                  CALL dbcsr_block_conjg(send_data_area, &
                                         row_size=nze, col_size=1, lb=sd_pos)
               ENDIF
            ELSE
               CALL dbcsr_block_copy_aa(send_data_area, sm%data_area, row_size, col_size, &
                                        sd_pos, bp, scale_value)
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
!$OMP END PARALLEL

      ! Exchange the data and metadata structures. In the interesting cases (square grids, row col distribution same),
      ! there are only very few processors that need to exchange data.
      ! The hybrid_alltoall deals with this by doing point to point communication
      CALL timeset(routineN//"_data", handle2)
      CALL hybrid_alltoall_any(send_data_area, total_send_count(2, :), sd_disp(:)-1, &
                               recv_data_area, total_recv_count(2, :), rd_disp(:)-1, &
                               mp_obj, &
                               most_ptp=.TRUE., remainder_ptp=.TRUE., no_hybrid=.FALSE.)
      CALL hybrid_alltoall_i1( &
         send_meta(:), metalen*total_send_count(1, :), sm_disp(:)-1, &
         recv_meta(:), metalen*total_recv_count(1, :), rm_disp(:)-1, &
         most_ptp=.TRUE., remainder_ptp=.TRUE., no_hybrid=.FALSE., &
         mp_env=mp_obj)
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
            nze = row_blk_size(ABS(stored_row)) &
                  *col_blk_size(ABS(stored_col))
            blk = blks(row_img, col_img)
            blks(row_img, col_img) = blks(row_img, col_img)+1
            IF (stored_col .GT. 0) THEN
               blk_p = data_p
               data_p = data_p+nze
            ELSE
               blk_p = prev_blk_p
            ENDIF
            blk_ps(row_img, col_img) = blk_ps(row_img, col_img)+nze
            ums%mats(row_img, col_img)%m%wms(1)%row_i(blk) = ABS(stored_row)
            ums%mats(row_img, col_img)%m%wms(1)%col_i(blk) = ABS(stored_col)
            ums%mats(row_img, col_img)%m%wms(1)%blk_p(blk) = &
               SIGN(rd_disp(src_p)+ABS(stored_blk_p)-1, stored_blk_p)
            prev_blk_p = blk_p
         ENDDO
      ENDDO

      ! Finalize the actual imaged matrices from the work matrices
      !
      DO row_img = 1, nrow_images
         DO col_img = 1, ncol_images
            ums%mats(row_img, col_img)%m%wms(1)%lastblk = blks(row_img, col_img)-1
            ums%mats(row_img, col_img)%m%wms(1)%datasize = blk_ps(row_img, col_img)-1
            !
            CALL dbcsr_data_set_size_referenced( &
               ums%mats(row_img, col_img)%m%wms(1)%data_area, &
               ums%mats(row_img, col_img)%m%wms(1)%datasize)
            IF (nrow_images .EQ. 1 .AND. ncol_images .EQ. 1 .OR. nocopy) THEN
               CALL dbcsr_special_finalize(ums%mats(row_img, col_img), reshuffle=.FALSE.)
            ELSE
               CALL dbcsr_special_finalize(ums%mats(row_img, col_img), reshuffle=.TRUE.)
            ENDIF

            ! Save the home process and image row and column
            CALL image_calculator(target_imgdist, &
                                  ums%mats(row_img, col_img)%m%index(dbcsr_slot_home_prow), &
                                  ums%mats(row_img, col_img)%m%index(dbcsr_slot_home_rowi), &
                                  ums%mats(row_img, col_img)%m%index(dbcsr_slot_home_pcol), &
                                  ums%mats(row_img, col_img)%m%index(dbcsr_slot_home_coli), &
                                  vprow=ums%mats(row_img, col_img)%m%index(dbcsr_slot_home_vprow), &
                                  vpcol=ums%mats(row_img, col_img)%m%index(dbcsr_slot_home_vpcol), &
                                  myrowi=row_img, mycoli=col_img, &
                                  shifting=predist_type)
         ENDDO
      ENDDO

      ! Deallocate shared temporary buffers
      !
      DEALLOCATE (send_count, recv_count)
      DEALLOCATE (total_send_count, total_recv_count)
      DEALLOCATE (sdp, smp, sd_disp, sm_disp)
      DEALLOCATE (rd_disp, rm_disp)
      DEALLOCATE (all_total_send_offset)
      DEALLOCATE (blk_ps, blks)
      DEALLOCATE (recv_meta, send_meta)

      CALL dbcsr_data_release(send_data_area)
      CALL dbcsr_data_release(received_data_area)
      CALL dbcsr_data_release(recv_data_area)

      CALL timestop(handle)
   END SUBROUTINE make_images

! **************************************************************************************************
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
! **************************************************************************************************
   SUBROUTINE dbcsr_make_images_dense(images, new_rdist, &
                                      row_map, col_map, join_cols, join_rows, new_template)
      TYPE(dbcsr_2d_array_type), INTENT(INOUT)           :: images
      TYPE(dbcsr_imagedistribution_obj), INTENT(INOUT)   :: new_rdist
      TYPE(array_i1d_obj), INTENT(IN)                    :: row_map, col_map
      LOGICAL, INTENT(IN)                                :: join_cols, join_rows
      TYPE(dbcsr_obj), INTENT(IN)                        :: new_template

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_images_dense', &
         routineP = moduleN//':'//routineN
      LOGICAL, PARAMETER                                 :: dbg = .FALSE.

      INTEGER                                            :: handle, mat_col, mat_row, mat_vpcol, &
                                                            mat_vprow
      INTEGER, DIMENSION(:), POINTER                     :: und_col_blk_offsets, und_row_blk_offsets
      INTEGER, DIMENSION(dbcsr_meta_size)                :: old_meta
      REAL(kind=dp)                                      :: cs
      TYPE(array_i1d_obj)                                :: dense_local_vcols, dense_local_vrows, &
                                                            und_local_vcols, und_local_vrows
      TYPE(dbcsr_imagedistribution_obj)                  :: old_rdist
      TYPE(dbcsr_obj)                                    :: tmp_mat

!   ---------------------------------------------------------------------------

      CALL timeset(routineN, handle)
      old_rdist = images%image_dist
      !
      DO mat_row = 1, images%image_dist%i%row_decimation
         DO mat_col = 1, images%image_dist%i%col_decimation
            IF (dbg) THEN
               cs = dbcsr_checksum(images%mats(mat_row, mat_col))
               WRITE (*, *) routineN//" cs pre", cs
            ENDIF
            mat_vprow = images%mats(mat_row, mat_col)%m%index(dbcsr_slot_home_vprow)
            mat_vpcol = images%mats(mat_row, mat_col)%m%index(dbcsr_slot_home_vpcol)
            und_row_blk_offsets => array_data(images%mats(mat_row, mat_col)%m%row_blk_offset)
            und_col_blk_offsets => array_data(images%mats(mat_row, mat_col)%m%col_blk_offset)
            CALL dbcsr_get_local_vrows(old_rdist, und_local_vrows, mat_vprow)
            CALL dbcsr_get_local_vcols(old_rdist, und_local_vcols, mat_vpcol)

            CALL dbcsr_get_local_vrows(new_rdist, dense_local_vrows, mat_vprow)
            CALL dbcsr_get_local_vcols(new_rdist, dense_local_vcols, mat_vpcol)
            ! The old matrix has to be remembered so it is copied to
            ! tmp_mat.
            old_meta(:) = images%mats(mat_row, mat_col)%m%index(1:dbcsr_meta_size)
            CALL dbcsr_init(tmp_mat)
            tmp_mat = images%mats(mat_row, mat_col)
            CALL dbcsr_init(images%mats(mat_row, mat_col))
            CALL dbcsr_create(images%mats(mat_row, mat_col), template=new_template)
            images%mats(mat_row, mat_col)%m%index(dbcsr_slot_home_prow &
                                                  :dbcsr_slot_home_vpcol) = &
               old_meta(dbcsr_slot_home_prow:dbcsr_slot_home_vpcol)
            CALL dbcsr_make_dense_low(tmp_mat, images%mats(mat_row, mat_col), &
                                      array_data(und_local_vrows), array_data(und_local_vcols), &
                                      und_row_blk_offsets, und_col_blk_offsets, &
                                      array_data(dense_local_vrows), &
                                      array_data(dense_local_vcols), &
                                      array_data(new_template%m%row_blk_offset), &
                                      array_data(new_template%m%col_blk_offset), &
                                      array_data(row_map), array_data(col_map), join_rows, join_cols)
            !
            CALL dbcsr_index_prune_deleted(images%mats(mat_row, mat_col))
            !
            CALL dbcsr_release(tmp_mat)
            IF (dbg) THEN
               cs = dbcsr_checksum(images%mats(mat_row, mat_col))
               WRITE (*, *) routineN//" cs pst", cs
            ENDIF
         ENDDO
      ENDDO
      CALL dbcsr_image_dist_release(images%image_dist)
      images%image_dist = new_rdist
      CALL dbcsr_image_dist_hold(images%image_dist)
      CALL timestop(handle)
   END SUBROUTINE dbcsr_make_images_dense

! **************************************************************************************************
!> \brief Multiplies two DBCSR matrices
!>
!> \param[in] left_set             set of imaged left matrices
!> \param[in] right_set            set of imaged right matrices
!> \param[out] product_matrix      DBCSR product matrix
!> \param[in] retain_sparsity      (optional) retain the sparsity of the
!>                                 existing product matrix; default is no
!> \param filter_eps ...
!> \param[out] flop                (optional) effective flop
! **************************************************************************************************
   SUBROUTINE cannon_multiply_low(left_set, right_set, product_matrix, &
                                  retain_sparsity, &
                                  filter_eps, flop)
      TYPE(dbcsr_2d_array_type), POINTER                 :: left_set, right_set
      TYPE(dbcsr_obj), INTENT(INOUT)                     :: product_matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: retain_sparsity
      REAL(kind=real_8), INTENT(in), OPTIONAL            :: filter_eps
      INTEGER(KIND=int_8), INTENT(OUT)                   :: flop

      CHARACTER(len=*), PARAMETER :: routineN = 'cannon_multiply_low', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: idata = 1, ileft = 0, imeta = 2, &
                                                            iright = 2, M_L = 2, M_P = 1, M_R = 3

      INTEGER :: data_type, data_type_byte, grp, handle, handle2, handle3, i, ithread, &
         left_col_image, left_col_mult, left_col_nimages, left_dst_icol, left_dst_irow, &
         left_dst_p, left_dst_pcol, left_dst_prow, left_dst_vcol, left_dst_vrow, left_max_nblks, &
         left_max_nze, left_myfirstvcol, left_myfirstvrow, left_mypcol, left_myprow, left_npcols, &
         left_nprows, left_recv_icol, left_recv_irow, left_recv_p, left_recv_pcol, left_recv_prow, &
         left_recv_vcol, left_recv_vrow, left_row_image, left_row_mult, left_row_nimages, &
         left_send_icol, left_send_irow, left_send_p, left_send_pcol, left_send_prow
      INTEGER :: left_send_vcol, left_send_vrow, left_src_icol, left_src_irow, left_src_p, &
         left_src_pcol, left_src_prow, left_src_vcol, left_src_vrow, metronome, min_nimages, &
         mp_group, mynode, nblkrows_total, nblkrows_used, nsteps_k, nthreads, numnodes, nvirt_k, &
         output_unit, right_col_image, right_col_mult, right_col_nimages, right_dst_icol, &
         right_dst_irow, right_dst_p, right_dst_pcol, right_dst_prow, right_dst_vcol, &
         right_dst_vrow, right_max_nblks, right_max_nze, right_myfirstvcol, right_myfirstvrow, &
         right_mypcol, right_myprow, right_npcols, right_nprows, right_recv_icol, right_recv_irow
      INTEGER :: right_recv_p, right_recv_pcol, right_recv_prow, right_recv_vcol, right_recv_vrow, &
         right_row_image, right_row_mult, right_row_nimages, right_send_icol, right_send_irow, &
         right_send_p, right_send_pcol, right_send_prow, right_send_vcol, right_send_vrow, &
         right_src_icol, right_src_irow, right_src_p, right_src_pcol, right_src_prow, &
         right_src_vcol, right_src_vrow, row, size_guess, stat, threads_finished, &
         threads_finished_read, v_k, v_ki
      INTEGER(KIND=int_8)                                :: flop_single, flop_total, mem
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: row_counts, total_row_counts
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: left_sizes, my_sizes, right_sizes
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)        :: all_sizes
      INTEGER, DIMENSION(3, 2)                           :: mp_rc_groups
      INTEGER, DIMENSION(:), POINTER :: left_data_rr, left_data_sr, left_index_rp, left_index_rr, &
         left_index_sp, left_index_sr, local_rows, right_data_rr, right_data_sr, right_index_rp, &
         right_index_rr, right_index_sp, right_index_sr
      INTEGER, DIMENSION(:, :), POINTER                  :: left_pgrid, product_pgrid, right_pgrid
      INTEGER, SAVE                                      :: mult_id = 0
      LOGICAL                                            :: keep_sparsity, list_indexing, &
                                                            otf_filtering

!$    REAL(KIND=real_8)                        :: left_fill, right_fill
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
      NULLIFY (trs_stackbuf_calc, trs_stackbuf_comm)
      !
      ALLOCATE (left_buffer_2, right_buffer_2)
      mult_id = mult_id+1

      IF (PRESENT(retain_sparsity)) THEN
         keep_sparsity = retain_sparsity
      ELSE
         keep_sparsity = .FALSE.
      ENDIF
      otf_filtering = PRESENT(filter_eps)

!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (multrec, nthreads, product_matrix)
!$OMP MASTER
      nthreads = 1
!$    nthreads = OMP_GET_NUM_THREADS()
      IF (.NOT. ASSOCIATED(product_matrix%m%wms)) &
         CPABORT("Work matrices do not exist")
      IF (SIZE(product_matrix%m%wms) .NE. nthreads) &
         CPABORT("Work matrices not correctly sized.")
      ALLOCATE (multrec(0:nthreads-1))
!$OMP END MASTER
!$OMP END PARALLEL

      output_unit = default_output_unit
      flop_total = 0
      ! Set up variables
      data_type = dbcsr_get_data_type(product_matrix)
      data_type_byte = dbcsr_datatype_sizeof(data_type)
      left_row_nimages = left_set%image_dist%i%row_decimation
      left_row_mult = left_set%image_dist%i%row_multiplicity
      left_col_nimages = left_set%image_dist%i%col_decimation
      left_col_mult = left_set%image_dist%i%col_multiplicity
      right_row_nimages = right_set%image_dist%i%row_decimation
      right_row_mult = right_set%image_dist%i%row_multiplicity
      right_col_nimages = right_set%image_dist%i%col_decimation
      right_col_mult = right_set%image_dist%i%col_multiplicity
      left_mp_obj = dbcsr_distribution_mp(left_set%image_dist%i%main)
      right_mp_obj = dbcsr_distribution_mp(right_set%image_dist%i%main)
      product_mp_obj = dbcsr_distribution_mp(product_matrix%m%dist)
      numnodes = dbcsr_mp_numnodes(product_mp_obj)
      mynode = dbcsr_mp_mynode(product_mp_obj)
      left_nprows = dbcsr_mp_nprows(left_mp_obj)
      left_npcols = dbcsr_mp_npcols(left_mp_obj)
      left_myprow = dbcsr_mp_myprow(left_mp_obj)
      left_mypcol = dbcsr_mp_mypcol(left_mp_obj)
      left_myfirstvrow = dbcsr_mp_myprow(left_mp_obj)*left_row_nimages
      left_myfirstvcol = dbcsr_mp_mypcol(left_mp_obj)*left_col_nimages
      right_nprows = dbcsr_mp_nprows(right_mp_obj)
      right_npcols = dbcsr_mp_npcols(right_mp_obj)
      right_myprow = dbcsr_mp_myprow(right_mp_obj)
      right_mypcol = dbcsr_mp_mypcol(right_mp_obj)
      right_myfirstvrow = dbcsr_mp_myprow(right_mp_obj)*right_row_nimages
      right_myfirstvcol = dbcsr_mp_mypcol(right_mp_obj)*right_col_nimages
      mp_group = dbcsr_mp_group(product_mp_obj)
      left_pgrid => dbcsr_mp_pgrid(left_mp_obj)
      right_pgrid => dbcsr_mp_pgrid(right_mp_obj)
      product_pgrid => dbcsr_mp_pgrid(product_mp_obj)
      CALL dbcsr_mp_grid_setup(product_mp_obj)
      CALL dbcsr_mp_grid_setup(left_mp_obj)
      CALL dbcsr_mp_grid_setup(right_mp_obj)
      IF (dbcsr_mp_has_subgroups(product_mp_obj)) THEN
         mp_rc_groups(M_P, 1:2) = (/dbcsr_mp_my_row_group(product_mp_obj), &
                                    dbcsr_mp_my_col_group(product_mp_obj)/)
      ENDIF
      IF (dbcsr_mp_has_subgroups(left_mp_obj)) THEN
         mp_rc_groups(M_L, 1:2) = (/dbcsr_mp_my_row_group(left_mp_obj), &
                                    dbcsr_mp_my_col_group(left_mp_obj)/)
      ENDIF
      IF (dbcsr_mp_has_subgroups(right_mp_obj)) THEN
         mp_rc_groups(M_R, 1:2) = (/dbcsr_mp_my_row_group(right_mp_obj), &
                                    dbcsr_mp_my_col_group(right_mp_obj)/)
      ENDIF
      !
      ! Dummy checks
      ! left/right matching
      IF (left_col_nimages .NE. right_row_mult) &
         CPABORT("Left/Right image mismatch")
      IF (left_col_mult .NE. right_row_nimages) &
         CPABORT("Left/Right image mismatch")
      IF (left_col_nimages*left_npcols .NE. right_row_nimages*right_nprows) &
         CPABORT("Left/Right total mismatch")
      ! product/left matching
      IF (left_row_mult*dbcsr_mp_nprows(product_mp_obj) .NE. left_row_nimages*left_nprows) &
         CPABORT("Product/Left total mismatch")
      ! product/left matching
      IF (right_col_mult*dbcsr_mp_npcols(product_mp_obj) .NE. right_col_nimages*right_npcols) &
         CPABORT("Product/Right total mismatch")
      ! Limitations
      IF (left_row_nimages .NE. 1) &
         CPABORT("Product/Left matrix process grid mismatch")
      IF (left_row_mult .NE. 1) &
         CPABORT("Product/Left matrix process grid mismatch")
      IF (right_col_nimages .NE. 1) &
         CPABORT("Product/Right matrix process grid mismatch")
      IF (right_col_mult .NE. 1) &
         CPABORT("Product/Right matrix process grid mismatch")

      dbcsr_mpi_statistics%nimages = MAX(dbcsr_mpi_statistics%nimages, left_row_nimages*left_col_nimages)
      dbcsr_mpi_statistics%nimages = MAX(dbcsr_mpi_statistics%nimages, right_row_nimages*right_col_nimages)
      !
      ! Exchange size data
      ALLOCATE (my_sizes(4, MAX(left_row_nimages, right_row_nimages), &
                         MAX(left_col_nimages, right_col_nimages)))
      my_sizes(:, :, :) = 0
      DO left_row_image = 1, left_row_nimages
         DO left_col_image = 1, left_col_nimages
            my_sizes(idata+ileft, left_row_image, left_col_image) &
               = dbcsr_data_get_size_referenced( &
                 left_set%mats(left_row_image, left_col_image)%m%data_area)
            my_sizes(imeta+ileft, left_row_image, left_col_image) = &
               left_set%mats(left_row_image, left_col_image)%m%index &
               (dbcsr_slot_size)
         ENDDO
      ENDDO

      DO right_row_image = 1, right_row_nimages
         DO right_col_image = 1, right_col_nimages
            my_sizes(idata+iright, right_row_image, right_col_image) &
               = dbcsr_data_get_size_referenced( &
                 right_set%mats(right_row_image, right_col_image)%m%data_area)
            my_sizes(imeta+iright, right_row_image, right_col_image) = &
               right_set%mats(right_row_image, right_col_image)%m%index &
               (dbcsr_slot_size)
         ENDDO
      ENDDO

      ALLOCATE (all_sizes(4, LBOUND(my_sizes, 2):UBOUND(my_sizes, 2), &
                          LBOUND(my_sizes, 3):UBOUND(my_sizes, 3), 0:numnodes-1))
      CALL mp_allgather(my_sizes, all_sizes, mp_group)
      !
      ! Count the maximum possible multiplies per row for on-the-fly
      ! filtering.
      per_row_eps: IF (.NOT. otf_filtering) THEN
         ! These arrays must be valid when passed to called subroutines.
         ALLOCATE (left_norms(0), right_norms(0), row_max_epss(0), stat=stat)
         IF (stat .NE. 0) &
            CPABORT("Could not allocate memory")
      ELSE
         IF (careful_mod) THEN
            IF (left_set%mats(1, 1)%m%bcsc) &
               CPABORT("Can not do on-the-fly filtering with CSC-indexed matrices.")
         ENDIF
         IF (dbcsr_has_local_row_index(left_set%mats(1, 1))) THEN
            nblkrows_used = dbcsr_nblkrows_local(left_set%mats(1, 1))
         ELSE
            nblkrows_used = dbcsr_nblkrows_total(left_set%mats(1, 1))
         ENDIF
         nblkrows_total = dbcsr_nblkrows_total(left_set%mats(1, 1))
         ALLOCATE (row_max_epss(nblkrows_total), stat=stat)
         IF (stat .NE. 0) &
            CPABORT("Could not allocate memory for left epsilons")
         ALLOCATE (row_counts(nblkrows_used), stat=stat)
         IF (stat .NE. 0) &
            CPABORT("Could not allocate memory for left row counts")
         ! The summation could be done prow-locally but it would
         ! complicate the pre-row eps calculation.
         ALLOCATE (total_row_counts(nblkrows_total), stat=stat)
         IF (stat .NE. 0) &
            CPABORT("Could not allocate memory for left row counts")
         ! Each prow member matrix (npcols * row_images) counts the
         ! blocks present in each of its rows.
         total_row_counts(:) = 0
         DO left_row_image = 1, left_row_nimages
            DO left_col_image = 1, left_col_nimages
               list_indexing = &
                  left_set%mats(left_row_image, left_col_image)%m%list_indexing
               IF (careful_mod) THEN
                  IF (list_indexing) THEN
                     IF ((left_set%mats(left_row_image, left_col_image)%m%nblks)*3 .NE. &
                         SIZE(left_set%mats(left_row_image, left_col_image)%m%coo_l)) &
                        CPABORT("Row count mismatch")
                  ELSE
                     IF (nblkrows_used+1 .NE. SIZE(left_set%mats(left_row_image, left_col_image)%m%row_p)) &
                        CPABORT("Row count mismatch")
                  ENDIF
               ENDIF
               IF (list_indexing) THEN
                  CALL count_bins( &
                     left_set%mats(left_row_image, left_col_image)%m%nblks, &
                     left_set%mats(left_row_image, left_col_image)%m%coo_l(1::3), &
                     nblkrows_used, row_counts)
               ELSE
                  CALL dbcsr_count_row_index( &
                     left_set%mats(left_row_image, left_col_image)%m%row_p, &
                     row_counts, nblkrows_used)
               ENDIF
               IF (dbcsr_has_local_row_index(left_set%mats(left_row_image, left_col_image))) THEN
                  local_rows => array_data(left_set%mats(left_row_image, left_col_image)%m%local_rows)
                  IF (SIZE(local_rows) .NE. SIZE(row_counts)) &
                     CPABORT("Mismatch in number of local rows.")
                  total_row_counts(local_rows) = total_row_counts(local_rows) &
                                                 +row_counts(1:nblkrows_used)
               ELSE
                  total_row_counts(:) = total_row_counts(:) &
                                        +row_counts(:)
               ENDIF
            ENDDO
         ENDDO
         ! The counted blocks are then summed up
         CALL mp_sum(total_row_counts, mp_group)
         ! and used to determine the maximum per-block epsilon.
         DO row = 1, nblkrows_total
            row_max_epss(row) &
               = REAL(filter_eps &
                      /REAL(MAX(1, total_row_counts(row)), KIND=KIND(row_max_epss)), &
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
      min_nimages = MIN(left_col_nimages, right_row_nimages)
      nvirt_k = left_npcols*left_col_nimages
      nsteps_k = nvirt_k/min_nimages
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
            CALL image_calculator(left_set%image_dist, &
                                  prow=left_dst_prow, pcol=left_dst_pcol, &
                                  rowi=left_dst_irow, coli=left_dst_icol, &
                                  myvprow=left_src_vrow, myvpcol=left_src_vcol, &
                                  shifting='l')
            left_dst_p = left_pgrid(left_dst_prow, left_dst_pcol)
            left_sizes(idata, left_src_vrow, left_src_vcol) = &
               all_sizes( &
               idata+ileft, left_dst_irow, left_dst_icol, left_dst_p)
            left_sizes(imeta, left_src_vrow, left_src_vcol) = &
               all_sizes( &
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
            CALL image_calculator(right_set%image_dist, &
                                  prow=right_dst_prow, pcol=right_dst_pcol, &
                                  rowi=right_dst_irow, coli=right_dst_icol, &
                                  myvprow=right_src_vrow, myvpcol=right_src_vcol, &
                                  shifting='r')
            right_dst_p = right_pgrid(right_dst_prow, right_dst_pcol)
            right_sizes(idata, right_src_vrow, right_src_vcol) = &
               all_sizes( &
               idata+iright, right_dst_irow, right_dst_icol, right_dst_p)
            right_sizes(imeta, right_src_vrow, right_src_vcol) = &
               all_sizes( &
               imeta+iright, right_dst_irow, right_dst_icol, right_dst_p)
         ENDDO
      ENDDO
      !
      ! Setup product work areas
      left_max_nze = MAXVAL(all_sizes(idata+ileft, :, :, :))
      left_max_nblks = MAXVAL(all_sizes(imeta+ileft, :, :, :))
      right_max_nze = MAXVAL(all_sizes(idata+iright, :, :, :))
      right_max_nblks = MAXVAL(all_sizes(imeta+iright, :, :, :))
      !!
      ithread = 0
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          PRIVATE (i, size_guess, &
!$OMP                   left_fill, right_fill, ithread) &
!$OMP          SHARED (product_matrix, left_max_nze, right_max_nze) &
!$OMP          SHARED (left_set, right_set, &
!$OMP                 left_col_nimages, right_row_nimages) &
!$OMP          SHARED (nthreads, keep_sparsity, mynode)
      !
!$    ithread = OMP_GET_THREAD_NUM()
      ! The work arrays have to be setup (actually, not quite sure).
      i = ithread+1
      size_guess = product_matrix%m%wms(i)%datasize ! Should be minimal
      CALL dbcsr_data_ensure_size(product_matrix%m%wms(i)%data_area, &
                                  size_guess)
      CALL dbcsr_data_set_size_referenced(product_matrix%m%wms(i)%data_area, &
                                          product_matrix%m%wms(i)%datasize)
      ! XXXXXXX a quick fix right now, allocation with size 1 might actually not be needed at all,
      !         but something expects this to be associated
      CALL ensure_array_size(product_matrix%m%wms(i)%row_i, ub=1)
      CALL ensure_array_size(product_matrix%m%wms(i)%col_i, ub=1)
      CALL ensure_array_size(product_matrix%m%wms(i)%blk_p, ub=1)
!$OMP END PARALLEL

      ! update capacity of memory-pools
      IF (has_acc) THEN
         CALL dbcsr_mempool_ensure_capacity(memtype_abpanel_1%pool, &
                                            capacity=left_row_mult*left_col_nimages+right_row_nimages*right_col_mult)
         CALL dbcsr_mempool_ensure_capacity(memtype_abpanel_2%pool, &
                                            capacity=left_row_mult*left_col_nimages+right_row_nimages*right_col_mult)
      ENDIF

      !
      ! Setup the left buffer matrices
      !
      CALL buffer_matrices_ensure_size(left_set, index_size=left_max_nblks, &
                                       data_size=left_max_nze)

      CALL setup_buffer_matrices(left_buffer_2, left_row_mult, left_col_nimages, &
                                 left_set%mats(1, 1), index_size=left_max_nblks, &
                                 data_size=left_max_nze)
      IF (otf_filtering) THEN
         ALLOCATE (left_norms(left_max_nblks), stat=stat)
         IF (stat .NE. 0) &
            CPABORT("Could not allocate memory for left norms")
         IF (stat .NE. 0) otf_filtering = .FALSE.
      ENDIF
      left_buffer_calc => left_set
      left_buffer_comm => left_buffer_2
      ALLOCATE (left_data_sr(left_col_nimages))
      ALLOCATE (left_index_sr(left_col_nimages))
      ALLOCATE (left_data_rr(left_col_nimages))
      ALLOCATE (left_index_rr(left_col_nimages))
      left_data_sr = mp_request_null
      left_data_rr = mp_request_null
      left_index_sr = mp_request_null
      left_index_rr = mp_request_null

      ! Setup buffers for right matrix
      CALL buffer_matrices_ensure_size(right_set, index_size=right_max_nblks, &
                                       data_size=right_max_nze)

      CALL setup_buffer_matrices(right_buffer_2, right_row_nimages, right_col_mult, &
                                 right_set%mats(1, 1), index_size=right_max_nblks, data_size=right_max_nze)
      IF (otf_filtering) THEN
         ALLOCATE (right_norms(right_max_nblks), stat=stat)
         IF (stat .NE. 0) &
            CPWARN("Could not allocate memory for right norms")
         IF (stat .NE. 0) otf_filtering = .FALSE.
      ENDIF
      right_buffer_calc => right_set
      right_buffer_comm => right_buffer_2
      ALLOCATE (right_data_sr(right_row_nimages))
      ALLOCATE (right_index_sr(right_row_nimages))
      ALLOCATE (right_data_rr(right_row_nimages))
      ALLOCATE (right_index_rr(right_row_nimages))
      right_data_sr = mp_request_null
      right_data_rr = mp_request_null
      right_index_sr = mp_request_null
      right_index_rr = mp_request_null
      !
!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED (left_buffer_comm, right_buffer_comm, product_matrix,&
!$OMP         keep_sparsity, filter_eps, row_max_epss, multrec, nthreads, &
!$OMP         right_data_sr, right_data_rr, left_data_sr, left_data_rr,&
!$OMP         right_index_sr, right_index_rr, left_index_sr, left_index_rr), &
!$OMP PRIVATE(ithread)
      ithread = 0
!$    ithread = OMP_GET_THREAD_NUM()
      ALLOCATE (multrec(ithread)%p)
      CALL dbcsr_mm_multrec_init(multrec(ithread)%p, &
                                 left_buffer_comm%mats(1, 1)%m, &
                                 right_buffer_comm%mats(1, 1)%m, &
                                 product_matrix%m, &
                                 keep_sparsity=keep_sparsity, &
                                 eps=filter_eps, &
                                 row_max_epss=row_max_epss, &
                                 block_estimate=MAX(product_matrix%m%nblks, &
                                                    left_buffer_comm%mats(1, 1)%m%nblks, &
                                                    right_buffer_comm%mats(1, 1)%m%nblks)/nthreads, &
                                 right_row_blk_size=array_data(right_buffer_comm%mats(1, 1)%m%row_blk_size))
!$OMP END PARALLEL
      !
      ! Setup indexing
      CALL setup_rec_index(left_set, left_row_nimages, left_col_nimages)
      CALL setup_rec_index(right_set, right_row_nimages, right_col_nimages)
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
         IF (debug_mod) WRITE (*, '(1X,A,3(1X,A,1X,I5))') routineN, &
            "step", metronome, &
            "first k", metronome*min_nimages, &
            "last k", (metronome+1)*min_nimages-1
         ! Wait for right matrix transfer completion. Wait in all but
         ! the first loop iteration.
         CALL timeset(routineN//"_metrocomm1", handle2)
         wait_right: IF (metronome .GT. 1) THEN
            IF (debug_mod) WRITE (*, '(1X,A)') routineN//" waiting for right"
            !
            CALL mp_waitall(right_data_sr)
            CALL mp_waitall(right_data_rr)
            CALL mp_waitall(right_index_sr)
            CALL mp_waitall(right_index_rr)
         ENDIF wait_right
         CALL timestop(handle2)
         ! Right matrix transfer. Transfer in all but the last loop
         ! iteration.
         xfer_right: IF (metronome .LT. nsteps_k) THEN
            DO v_ki = 0, right_row_nimages-1
               v_k = metronome*min_nimages+v_ki
               ! Calculate the process to send to.  It's the virtual
               ! process row -min_nimages up (i.e., smaller row number)
               ! from me.
               CALL image_calculator(right_set%image_dist, &
                                     prow=right_send_prow, rowi=right_send_irow, & ! output
                                     pcol=right_send_pcol, coli=right_send_icol, & ! output
                                     vprow=right_send_vrow, vpcol=right_send_vcol, & ! output
                                     ! myvprow goes through all of my (process row) images
                                     myvprow=v_ki+right_myfirstvrow, &
                                     myvpcol=right_myfirstvcol, & ! nothing happens in the columns
                                     ! send to process min_nimages up in the grid
                                     vprow_shift=-min_nimages, &
                                     shifting='0')
               ! Calculate which data I send.
               CALL image_calculator(right_set%image_dist, &
                                     prow=right_dst_prow, rowi=right_dst_irow, &
                                     pcol=right_dst_pcol, coli=right_dst_icol, &
                                     vprow=right_dst_vrow, vpcol=right_dst_vcol, &
                                     ! myvprows goes through all of my (process row) images
                                     myvprow=v_ki+right_myfirstvrow, &
                                     myvpcol=right_myfirstvcol, & ! nothing happens in the columns
                                     ! send what I got from min_nimages down, appropriate
                                     ! to the metronome tick
                                     vprow_shift=-min_nimages+metronome*min_nimages, &
                                     ! This is with relative shifting.
                                     shifting='R')
               right_dst_p = right_pgrid(right_dst_prow, right_dst_pcol)
               CALL dbcsr_data_set_pointer( &
                  area=right_data_sp, &
                  rsize=right_sizes(idata, right_dst_vrow, right_dst_vcol), &
                  csize=1, &
                  pointee=right_buffer_calc%mats(v_ki+1, 1)%m%data_area)
               right_index_sp => right_buffer_calc%mats( &
                                 v_ki+1, 1 &
                                 )%m%index(1: &
                                           right_sizes(imeta, right_dst_vrow, right_dst_vcol))
               !
               ! Calculate the process to receive from
               CALL image_calculator(right_set%image_dist, &
                                     prow=right_recv_prow, rowi=right_recv_irow, &
                                     pcol=right_recv_pcol, coli=right_recv_icol, &
                                     vprow=right_recv_vrow, vpcol=right_recv_vcol, &
                                     myvprow=v_ki+right_myfirstvrow, &
                                     myvpcol=right_myfirstvcol, &
                                     vprow_shift=+min_nimages, & ! just the opposite as "send to"
                                     shifting='0')
               ! Calculate which data I receive
               CALL image_calculator(right_set%image_dist, &
                                     prow=right_src_prow, rowi=right_src_irow, &
                                     pcol=right_src_pcol, coli=right_src_icol, &
                                     vprow=right_src_vrow, vpcol=right_src_vcol, &
                                     myvprow=v_ki+right_myfirstvrow, &
                                     myvpcol=right_myfirstvcol, &
                                     ! receive window moves with the metronome
                                     vprow_shift=metronome*min_nimages, &
                                     shifting='R')
               !
               IF (has_acc) THEN
                  CALL timeset(routineN//"_acc_sync_right", handle3)
                  CALL acc_event_synchronize(right_buffer_comm%mats(v_ki+1, 1)%m%data_area%d%acc_ready)
                  CALL timestop(handle3)
               ENDIF

               right_src_p = right_pgrid(right_src_prow, right_src_pcol)
               CALL dbcsr_data_set_pointer( &
                  area=right_data_rp, &
                  rsize=right_sizes(idata, right_src_vrow, right_src_vcol), &
                  csize=1, &
                  pointee=right_buffer_comm%mats(v_ki+1, 1)%m%data_area)
               right_index_rp => right_buffer_comm%mats( &
                                 v_ki+1, 1 &
                                 )%m%index(1: &
                                           right_sizes(imeta, right_src_vrow, right_src_vcol))
               !
               right_send_p = right_pgrid(right_send_prow, right_send_pcol)
               right_recv_p = right_pgrid(right_recv_prow, right_recv_pcol)
               ! These are column-communicator relative
               IF (dbcsr_mp_has_subgroups(right_mp_obj)) THEN
                  right_send_p = right_send_prow
                  right_recv_p = right_recv_prow
                  grp = dbcsr_mp_my_col_group(right_mp_obj)
               ELSE
                  grp = dbcsr_mp_group(right_mp_obj)
               ENDIF
               !
               CALL timeset(routineN//"_metrocomm2", handle2)
               CALL dbcsr_irecv_any(right_data_rp, right_recv_p, &
                                    grp, right_data_rr(v_ki+1), tag=right_src_vrow)
               CALL mp_irecv(right_index_rp, right_recv_p, &
                             grp, right_index_rr(v_ki+1), tag=right_src_vrow)
               CALL dbcsr_isend_any(right_data_sp, right_send_p, &
                                    grp, right_data_sr(v_ki+1), tag=right_dst_vrow)
               CALL mp_isend(right_index_sp, right_send_p, &
                             grp, right_index_sr(v_ki+1), tag=right_dst_vrow)
               dbcsr_mpi_statistics%nexchanged = dbcsr_mpi_statistics%nexchanged+1
               CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(1, :), &
                                         dbcsr_data_get_size(right_data_rp), &
                                         dbcsr_mpi_statistics%data_size_breakdown(:, :, 1), data_type_byte)
               CALL timestop(handle2)
            ENDDO
         ENDIF xfer_right
         !
         ! Repoint indices of right matrices
         calc_case_right: IF (metronome .GT. 1) THEN
            DO v_ki = 0, right_row_nimages-1
               CALL dbcsr_repoint_index(right_buffer_calc%mats(v_ki+1, 1)%m)
               right_buffer_calc%mats(v_ki+1, 1)%m%valid = .TRUE.
            ENDDO
         ENDIF calc_case_right
         !
         ! Wait for left matrix transfer completion. Wait in all but
         ! the first loop iteration.
         CALL timeset(routineN//"_metrocomm3", handle2)
         wait_left: IF (metronome .GT. 1) THEN
            IF (debug_mod) WRITE (*, '(1X,A)') routineN//" waiting for left"
            CALL mp_waitall(left_data_sr)
            CALL mp_waitall(left_data_rr)
            CALL mp_waitall(left_index_sr)
            CALL mp_waitall(left_index_rr)
         ENDIF wait_left
         CALL timestop(handle2)
         ! Left matrix transfer. Transfer in all but the last loop
         ! iteration.
         xfer_left: IF (metronome .LT. nsteps_k) THEN
            DO v_ki = 0, left_col_nimages-1
               v_k = metronome*min_nimages+v_ki
               ! Calculate the process to send to.
               CALL image_calculator(left_set%image_dist, &
                                     prow=left_send_prow, rowi=left_send_irow, & ! output
                                     pcol=left_send_pcol, coli=left_send_icol, & ! output
                                     vprow=left_send_vrow, vpcol=left_send_vcol, & ! output
                                     myvprow=left_myfirstvrow, & ! nothing happens in the rows
                                     ! go through all my column images
                                     myvpcol=v_ki+left_myfirstvcol, &
                                     ! send to process min_nimages left in the grid
                                     vpcol_shift=-min_nimages, &
                                     shifting='0')
               ! Calculate which data I send.
               CALL image_calculator(left_set%image_dist, &
                                     prow=left_dst_prow, rowi=left_dst_irow, &
                                     pcol=left_dst_pcol, coli=left_dst_icol, &
                                     vprow=left_dst_vrow, vpcol=left_dst_vcol, &
                                     myvprow=left_myfirstvrow, &
                                     ! go through all my column images
                                     myvpcol=v_ki+left_myfirstvcol, &
                                     ! send what I got from min_nimages left, appropriate
                                     ! to the metronome tick
                                     vpcol_shift=-min_nimages+metronome*min_nimages, &
                                     ! This is with relative shifting.
                                     shifting='L')
               !
               left_dst_p = left_pgrid(left_dst_prow, left_dst_pcol)
               CALL dbcsr_data_set_pointer( &
                  area=left_data_sp, &
                  rsize=left_sizes(idata, left_dst_vrow, left_dst_vcol), &
                  csize=1, &
                  pointee=left_buffer_calc%mats(1, v_ki+1)%m%data_area)
               left_index_sp => left_buffer_calc%mats( &
                                1, v_ki+1 &
                                )%m%index(1: &
                                          left_sizes(imeta, left_dst_vrow, left_dst_vcol))
               !
               ! Calculate the process to receive from
               CALL image_calculator(left_set%image_dist, &
                                     prow=left_recv_prow, rowi=left_recv_irow, &
                                     pcol=left_recv_pcol, coli=left_recv_icol, &
                                     vprow=left_recv_vrow, vpcol=left_recv_vcol, &
                                     myvprow=left_myfirstvrow, &
                                     myvpcol=v_ki+left_myfirstvcol, &
                                     vpcol_shift=+min_nimages, & ! just the opposite as "send to"
                                     shifting='0')
               ! Calculate which data I receive
               CALL image_calculator(left_set%image_dist, &
                                     prow=left_src_prow, rowi=left_src_irow, &
                                     pcol=left_src_pcol, coli=left_src_icol, &
                                     vprow=left_src_vrow, vpcol=left_src_vcol, &
                                     myvprow=left_myfirstvrow, &
                                     myvpcol=v_ki+left_myfirstvcol, &
                                     ! receive window moves with the metronome
                                     vpcol_shift=metronome*min_nimages, &
                                     shifting='L')
               !
               IF (has_acc) THEN
                  CALL timeset(routineN//"_acc_sync_left", handle3)
                  CALL acc_event_synchronize(left_buffer_comm%mats(1, v_ki+1)%m%data_area%d%acc_ready)
                  CALL timestop(handle3)
               ENDIF

               left_src_p = left_pgrid(left_src_prow, left_src_pcol)
               CALL dbcsr_data_set_pointer( &
                  area=left_data_rp, &
                  rsize=left_sizes(idata, left_src_vrow, left_src_vcol), &
                  csize=1, &
                  pointee=left_buffer_comm%mats(1, v_ki+1)%m%data_area)
               left_index_rp => left_buffer_comm%mats( &
                                1, v_ki+1 &
                                )%m%index(1: &
                                          left_sizes(imeta, left_src_vrow, left_src_vcol))
               !
               left_send_p = left_pgrid(left_send_prow, left_send_pcol)
               left_recv_p = left_pgrid(left_recv_prow, left_recv_pcol)
               ! These are column-communicator relative
               IF (dbcsr_mp_has_subgroups(left_mp_obj)) THEN
                  left_send_p = left_send_pcol
                  left_recv_p = left_recv_pcol
                  grp = dbcsr_mp_my_row_group(left_mp_obj)
               ELSE
                  grp = dbcsr_mp_group(left_mp_obj)
               ENDIF
               !
               CALL timeset(routineN//"_metrocomm4", handle2)
               CALL dbcsr_irecv_any(left_data_rp, left_recv_p, &
                                    grp, left_data_rr(v_ki+1), tag=left_src_vcol)
               CALL mp_irecv(left_index_rp, left_recv_p, &
                             grp, left_index_rr(v_ki+1), tag=left_src_vcol)
               CALL dbcsr_isend_any(left_data_sp, left_send_p, &
                                    grp, left_data_sr(v_ki+1), tag=left_dst_vcol)
               CALL mp_isend(left_index_sp, left_send_p, &
                             grp, left_index_sr(v_ki+1), tag=left_dst_vcol)
               dbcsr_mpi_statistics%nexchanged = dbcsr_mpi_statistics%nexchanged+1
               CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(2, :), &
                                         dbcsr_data_get_size(left_data_rp), &
                                         dbcsr_mpi_statistics%data_size_breakdown(:, :, 2), data_type_byte)
               CALL timestop(handle2)
            ENDDO
         ENDIF xfer_left
         !
         ! Repoint indices of left matrices and do the multiplications.
         calc_case_left: IF (metronome .GT. 1) THEN
            DO v_ki = 0, left_col_nimages-1
               CALL dbcsr_repoint_index(left_buffer_calc%mats(1, v_ki+1)%m)
               left_buffer_calc%mats(1, v_ki+1)%m%valid = .TRUE.
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
            IF (.FALSE.) WRITE (*, *) routineN//" TICK", v_ki
            IF (.TRUE. .OR. right_buffer_calc%mats(v_ki+1, 1)%m%local_indexing) THEN
               ! Since the right matrix is shifted vertically, the
               ! received data always has different notions of "local
               ! rows".  Thus the local_rows and global_rows must be
               ! recalculated.
               CALL dbcsr_reset_vlocals(right_buffer_calc%mats(v_ki+1, 1), &
                                        right_set%image_dist)
            ENDIF
            IF (.TRUE. .OR. left_buffer_calc%mats(1, v_ki+1)%m%local_indexing) THEN
               ! Since the right matrix is shifted vertically, the
               ! received data always has different notions of "local
               ! rows".  Thus the local_rows and global_rows must be
               ! recalculated.
               CALL dbcsr_reset_vlocals(left_buffer_calc%mats(1, v_ki+1), &
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
               CALL calculate_norms(right_buffer_calc%mats(v_ki+1, 1), &
                                    right_norms)
               CALL calculate_norms(left_buffer_calc%mats(1, v_ki+1), &
                                    left_norms)
            ENDIF
            !
            flop_single = 0
            threads_finished = 0

!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (left_buffer_calc, right_buffer_calc, &
!$OMP         v_ki, &
!$OMP         product_matrix, multrec,&
!$OMP         filter_eps, right_norms, left_norms, row_max_epss, &
!$OMP         keep_sparsity,threads_finished, &
!$OMP         right_data_sr, right_data_rr, right_index_sr, right_index_rr, &
!$OMP         left_data_sr, left_data_rr, left_index_sr, left_index_rr, &
!$OMP         use_comm_thread) &
!$OMP PRIVATE (ithread,nthreads,threads_finished_read) &
!$OMP FIRSTPRIVATE (metronome, nsteps_k, min_nimages) &
!$OMP REDUCTION (+: flop_single)
            ithread = 0; nthreads = 1
!$          ithread = omp_get_thread_num(); nthreads = omp_get_num_threads()

            IF (metronome == nsteps_k .AND. v_ki == min_nimages-1) &
               CALL dbcsr_mm_multrec_phaseout(multrec(ithread)%p)

            CALL dbcsr_mm_multrec_multiply(multrec(ithread)%p, &
                                           left=left_buffer_calc%mats(1, v_ki+1)%m, &
                                           right=right_buffer_calc%mats(v_ki+1, 1)%m, &
                                           flop=flop_single, &
                                           a_norms=left_norms, b_norms=right_norms)

            IF (metronome == nsteps_k .AND. v_ki == min_nimages-1) THEN
               CALL dbcsr_mm_multrec_finalize(multrec(ithread)%p)
               DEALLOCATE (multrec(ithread)%p)
            ENDIF

!$OMP ATOMIC
            threads_finished = threads_finished+1
            IF (use_comm_thread .AND. (ithread .EQ. 0)) THEN
               DO
! requires OMP 3.1 (e.g. gcc >=4.7), for correctness, otherwise we keep fingers crossed
#if defined _OPENMP && _OPENMP >= 200711
!$OMP                 ATOMIC READ
#endif
                  threads_finished_read = threads_finished
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
!$OMP END PARALLEL
            flop_total = flop_total+flop_single
         ENDDO

         CALL timestop(handle2)

         CALL dbcsr_switch(left_buffer_calc, left_buffer_comm)
         CALL dbcsr_switch(right_buffer_calc, right_buffer_comm)
         CALL dbcsr_switch(trs_stackbuf_calc, trs_stackbuf_comm)

      ENDDO grouped_k_index
      CALL m_memory(mem)
      max_memory = MAX(max_memory, REAL(mem))

      IF (has_acc) THEN
         CALL dbcsr_data_release(trs_stackbuf_1)
         CALL dbcsr_data_release(trs_stackbuf_2)
      END IF

      IF (ALLOCATED(right_norms)) THEN
         DEALLOCATE (right_norms)
      ENDIF
      IF (ALLOCATED(left_norms)) THEN
         DEALLOCATE (left_norms)
      ENDIF
      IF (ALLOCATED(row_max_epss)) THEN
         DEALLOCATE (row_max_epss)
      ENDIF
      !
      CALL dbcsr_destroy_array(right_buffer_2)
      CALL dbcsr_destroy_array(left_buffer_2)
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
      DEALLOCATE (left_data_rr, left_data_sr, left_index_rr, left_index_sr, &
                  right_data_rr, right_data_sr, right_index_rr, right_index_sr)
      !
      !
      IF (debug_mod) THEN
         v_ki = 0
         DO i = 1, SIZE(product_matrix%m%blk_p)
            v_ki = MAX(v_ki, ABS(product_matrix%m%blk_p(i)))
         ENDDO
         WRITE (*, *) routineN//" Actual final size", &
            LOG(REAL(dbcsr_data_get_size(product_matrix%m%data_area)))/LOG(10.0), &
            LOG(REAL(v_ki))/LOG(10.0)
      ENDIF
      !
      flop = flop_total
      DEALLOCATE (left_buffer_2, right_buffer_2)
      !
      CALL timestop(handle)
   END SUBROUTINE cannon_multiply_low

! **************************************************************************************************
!> \brief Multiplies two DBCSR matrices by means of RMA MPI operations.
!>        This algorithm is experimental and it should be not used in
!>        production runs.
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
!> \param keep_product_data ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE cannon_multiply_low_rma(imgdist_left, imgdist_right, &
                                      matrix_left, matrix_right, &
                                      product_matrix, &
                                      retain_sparsity, &
                                      filter_eps, flop, keep_product_data)
      TYPE(dbcsr_imagedistribution_obj), INTENT(INOUT)   :: imgdist_left, imgdist_right
      TYPE(dbcsr_obj), INTENT(IN)                        :: matrix_left, matrix_right
      TYPE(dbcsr_obj), INTENT(INOUT)                     :: product_matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: retain_sparsity
      REAL(kind=real_8), INTENT(IN), OPTIONAL            :: filter_eps
      INTEGER(KIND=int_8), INTENT(OUT)                   :: flop
      LOGICAL, INTENT(IN)                                :: keep_product_data

      CHARACTER(len=*), PARAMETER :: routineN = 'cannon_multiply_low_rma', &
         routineP = moduleN//':'//routineN

      INTEGER :: blk, col, data_type, data_type_byte, grp_gather, grp_left, grp_right, handle, &
         handle2, handle3, iget, index_row_max_epss, iproc, ithread, left_col_mult, &
         left_col_nimages, left_data_offset, left_data_offset_diag, left_data_size, left_data_win, &
         left_max_data_size, left_max_meta_size, left_meta_offset, left_meta_offset_diag, &
         left_meta_size, left_meta_win, left_myfirstvcol, left_myfirstvrow, left_mypcol, &
         left_myprow, left_npcols, left_nprows, left_nrequests, left_row_mult, left_row_nimages, &
         left_v_i, left_vcol, left_vrow, local_data_size_gather, local_meta_size_gather
      INTEGER :: meta_size_image, meta_size_scatter, metronome, min_nimages, mynode, mypcol, &
         myprow, nblkrows_local, nimages_scatter, nsteps_k, nthreads, numnodes, nvirt_k, &
         recv_vcol, recv_vcol_displ, recv_vrow, recv_vrow_displ, request_product_matrix, &
         right_col_mult, right_col_nimages, right_data_offset, right_data_offset_diag, &
         right_data_size, right_data_win, right_max_data_size, right_max_meta_size, &
         right_meta_offset, right_meta_offset_diag, right_meta_size, right_meta_win, &
         right_myfirstvcol, right_myfirstvrow, right_mypcol, right_myprow, right_npcols, &
         right_nprows, right_nrequests
      INTEGER :: right_row_mult, right_row_nimages, right_v_i, right_vcol, right_vrow, row, &
         size_guess, size_index_merged, size_index_unmerged, v_ci, v_k_loc, v_ki, v_ri
      INTEGER(KIND=int_8)                                :: mem
      INTEGER, ALLOCATABLE, DIMENSION(:) :: data_displ_gather, data_size_gather, &
         indices_common_unmerged, left_displ_vunmerged, left_get_requests_images, &
         left_max_meta_size_merged, meta_displ_gather, meta_size_gather, right_displ_vunmerged, &
         right_get_requests_images, right_max_meta_size_merged, v_k_scatter
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: left_recv_filtered, left_vunmerged, &
                                                            right_recv_filtered, right_vunmerged
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: left_vmerged, right_vmerged
      INTEGER, DIMENSION(:), POINTER                     :: blk_p, col_i
      INTEGER, DIMENSION(dbcsr_slot_nblkrows_total:&
         dbcsr_slot_nfullcols_local)                     :: left_global_indices, right_global_indices
      LOGICAL :: do_diagonal, do_k_scatter, do_rec, do_scatter, is_diagonal, keep_sparsity, &
         left_set_displ, otf_filtering, right_set_displ
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: common_unmerged
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :)           :: left_nofiltered, right_nofiltered
      REAL(kind=sp)                                      :: filter_eps_sp, min_row_max_epss
      REAL(kind=sp), ALLOCATABLE, DIMENSION(:)           :: row_max_epss
      REAL(kind=sp), ALLOCATABLE, DIMENSION(:, :)        :: left_norms, right_norms
      TYPE(dbcsr_1d_array_type)                          :: buffer_scatter
      TYPE(dbcsr_data_obj)                               :: data_get
      TYPE(dbcsr_mm_multrec_type_p), ALLOCATABLE, &
         DIMENSION(:)                                    :: multrec, multrec_scatter
      TYPE(dbcsr_mp_obj)                                 :: left_mp_obj, product_mp_obj, right_mp_obj
      TYPE(dbcsr_obj)                                    :: product_matrix_scatter
      TYPE(dbcsr_rma_buffers), POINTER                   :: rma_buffers, rma_buffers_calc, &
                                                            rma_buffers_comm

!   ---------------------------------------------------------------------------

      CALL timeset(routineN, handle)
      NULLIFY (rma_buffers_calc, rma_buffers_comm, rma_buffers)
      !
      IF (PRESENT(retain_sparsity)) THEN
         keep_sparsity = retain_sparsity
      ELSE
         keep_sparsity = .FALSE.
      ENDIF
      otf_filtering = PRESENT(filter_eps)
      !
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (nthreads)
!$OMP MASTER
      nthreads = 1
!$    nthreads = OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL
      !
      size_index_unmerged = dbcsr_slot_nblks
!$    size_index_unmerged = size_index_unmerged+nthreads+1
      size_index_merged = dbcsr_num_slots
!$    size_index_merged = size_index_merged+nthreads+1
      !
      ! Dummy checks
      IF (.NOT. ASSOCIATED(product_matrix%m%wms)) &
         CPABORT("Work matrices do not exist")
      IF (SIZE(product_matrix%m%wms) .NE. nthreads) &
         CPABORT("Work matrices not correctly sized.")
      IF (.NOT. rma_buffers_orig%left%is_valid .OR. &
          .NOT. rma_buffers_orig%right%is_valid .OR. &
          .NOT. ASSOCIATED(rma_buffers_orig%left%meta) .OR. &
          .NOT. ASSOCIATED(rma_buffers_orig%right%meta) .OR. &
          .NOT. ALLOCATED(left_refs_meta_size) .OR. &
          .NOT. ALLOCATED(left_refs_data_size) .OR. &
          .NOT. ALLOCATED(right_refs_meta_size) .OR. &
          .NOT. ALLOCATED(right_refs_data_size) .OR. &
          .NOT. ALLOCATED(left_local_refs_meta_size) .OR. &
          .NOT. ALLOCATED(left_local_refs_data_size) .OR. &
          .NOT. ALLOCATED(right_local_refs_meta_size) .OR. &
          .NOT. ALLOCATED(right_local_refs_data_size) .OR. &
          .NOT. ALLOCATED(left_refs_displ_unmerged) .OR. &
          .NOT. ALLOCATED(right_refs_displ_unmerged) .OR. &
          .NOT. ALLOCATED(left_local_refs_displ_unmerged) .OR. &
          .NOT. ALLOCATED(right_local_refs_displ_unmerged) .OR. &
          (otf_filtering .AND. use_mpi_filtering .AND. &
           (.NOT. (ALLOCATED(left_max_norms) .AND. ALLOCATED(right_max_norms) .AND. &
                   ALLOCATED(left_local_max_norms) .AND. ALLOCATED(right_local_max_norms))))) &
         CPABORT("No buffers associated for the RMA algo!")
      !
      ! Set up variables
      flop = 0
      data_type = dbcsr_get_data_type(product_matrix)
      data_type_byte = dbcsr_datatype_sizeof(data_type)
      left_row_nimages = imgdist_left%i%row_decimation
      left_row_mult = imgdist_left%i%row_multiplicity
      left_col_nimages = imgdist_left%i%col_decimation
      left_col_mult = imgdist_left%i%col_multiplicity
      right_row_nimages = imgdist_right%i%row_decimation
      right_row_mult = imgdist_right%i%row_multiplicity
      right_col_nimages = imgdist_right%i%col_decimation
      right_col_mult = imgdist_right%i%col_multiplicity
      left_mp_obj = dbcsr_distribution_mp(imgdist_left%i%main)
      right_mp_obj = dbcsr_distribution_mp(imgdist_right%i%main)
      product_mp_obj = dbcsr_distribution_mp(product_matrix%m%dist)
      numnodes = dbcsr_mp_numnodes(product_mp_obj)
      mynode = dbcsr_mp_mynode(product_mp_obj)
      myprow = dbcsr_mp_myprow(product_mp_obj)
      mypcol = dbcsr_mp_mypcol(product_mp_obj)
      left_nprows = dbcsr_mp_nprows(left_mp_obj)
      left_npcols = dbcsr_mp_npcols(left_mp_obj)
      left_myprow = dbcsr_mp_myprow(left_mp_obj)
      left_mypcol = dbcsr_mp_mypcol(left_mp_obj)
      left_myfirstvrow = left_myprow*left_row_nimages
      left_myfirstvcol = left_mypcol*left_col_nimages
      right_nprows = dbcsr_mp_nprows(right_mp_obj)
      right_npcols = dbcsr_mp_npcols(right_mp_obj)
      right_myprow = dbcsr_mp_myprow(right_mp_obj)
      right_mypcol = dbcsr_mp_mypcol(right_mp_obj)
      right_myfirstvrow = right_myprow*right_row_nimages
      right_myfirstvcol = right_mypcol*right_col_nimages
      grp_right = dbcsr_mp_my_col_group(right_mp_obj)
      grp_gather = dbcsr_mp_my_row_group(right_mp_obj)
      grp_left = dbcsr_mp_my_row_group(left_mp_obj)
      CALL dbcsr_mp_grid_setup(product_mp_obj)
      CALL dbcsr_mp_grid_setup(left_mp_obj)
      CALL dbcsr_mp_grid_setup(right_mp_obj)
      is_diagonal = myprow .EQ. mypcol
      do_scatter = right_npcols .GT. 1 .AND. .NOT. is_diagonal .AND. &
                   right_col_nimages .GT. 1 .AND. .NOT. has_acc
      do_diagonal = right_npcols .GT. 1 .AND. is_diagonal .AND. &
                    right_col_nimages .GT. 1 .AND. .NOT. has_acc
      !
      ! Dummy checks
      ! subcommunicators
      IF (.NOT. dbcsr_mp_has_subgroups(right_mp_obj)) &
         CPABORT("RMA requires rows subcommunicators for right matrix!")
      IF (.NOT. dbcsr_mp_has_subgroups(left_mp_obj)) &
         CPABORT("RMA requires columns subcommunicators for left matrix!")
      ! left/right matching
      IF (left_col_nimages .NE. right_row_mult) &
         CPABORT("Left/Right image mismatch")
      IF (left_col_mult .NE. right_row_nimages) &
         CPABORT("Left/Right image mismatch")
      IF (left_col_nimages*left_npcols .NE. right_row_nimages*right_nprows) &
         CPABORT("Left/Right total mismatch")
      ! product/left matching
      IF (left_row_mult*dbcsr_mp_nprows(product_mp_obj) .NE. left_row_nimages*left_nprows) &
         CPABORT("Product/Left total mismatch")
      ! product/left matching
      IF (right_col_mult*dbcsr_mp_npcols(product_mp_obj) .NE. right_col_nimages*right_npcols) &
         CPABORT("Product/Right total mismatch")
      !
      dbcsr_mpi_statistics%nimages = MAX(dbcsr_mpi_statistics%nimages, &
                                         left_row_nimages*left_col_nimages)
      dbcsr_mpi_statistics%nimages = MAX(dbcsr_mpi_statistics%nimages, &
                                         right_row_nimages*right_col_nimages)
      !
      ! The main transfer loop goes through the virtual rows/columns.
      ! The number of steps may be smaller if the grid dimension is very
      ! non-optimal (both left column images and right row images are >
      ! 1).
      min_nimages = MIN(left_col_nimages, right_row_nimages)
      nvirt_k = left_npcols*left_col_nimages
      nsteps_k = nvirt_k/min_nimages
      !
      do_rec = right_col_nimages .EQ. 1 .AND. left_row_nimages .EQ. 1
      !
      ! Exchange product matrix blocks coodinates
      IF (keep_sparsity .AND. right_npcols .GT. 1 .AND. &
          right_col_nimages .GT. 1 .AND. .NOT. has_acc) THEN
         IF (do_diagonal) THEN
            meta_size_scatter = product_matrix%m%index(dbcsr_slot_size)
         ENDIF
         CALL mp_bcast(meta_size_scatter, right_myprow, grp_gather)
         CALL ensure_array_size(local_meta_product_scatter, ub=meta_size_scatter, &
                                nocopy=.TRUE., memory_type=memtype_mpi_buffer)
         IF (do_diagonal) THEN
            local_meta_product_scatter(1:meta_size_scatter) = &
               product_matrix%m%index(1:meta_size_scatter)
         ENDIF
         CALL mp_ibcast(local_meta_product_scatter, right_myprow, &
                        grp_gather, request_product_matrix)
      ENDIF
      !
      product_matrix%m%nblks = 0
      product_matrix%m%nze = 0
      product_matrix%m%row_p(:) = 0
      !
      ! Count the maximum possible multiplies per row for on-the-fly filtering
      IF (otf_filtering) THEN
         ! Wait for counts (sent in make_buffers)
         CALL mp_wait(request_count_rows)
         nblkrows_local = SIZE(left_total_row_counts)
         ALLOCATE (row_max_epss(nblkrows_local))
         index_row_max_epss = 1
         filter_eps_sp = REAL(filter_eps, KIND=KIND(row_max_epss))
!$OMP PARALLEL DO DEFAULT (NONE) &
!$OMP SHARED(nblkrows_local,row_max_epss,filter_eps_sp,&
!$OMP        left_total_row_counts) &
!$OMP REDUCTION(MAX:index_row_max_epss)
         ! Determine the maximum per-block epsilon
         DO row = 1, nblkrows_local
            row_max_epss(row) = &
               filter_eps_sp/REAL(MAX(1, left_total_row_counts(row)), KIND=KIND(row_max_epss))
            ! Use integers for a fast comparison
            index_row_max_epss = MAX(index_row_max_epss, left_total_row_counts(row))
         ENDDO
!$OMP END PARALLEL DO
         min_row_max_epss = filter_eps_sp/REAL(index_row_max_epss, KIND=KIND(row_max_epss))
         DEALLOCATE (left_total_row_counts)
      ELSE
         !
         ! The array must be valid when passed to called subroutines.
         ALLOCATE (left_norms(0, min_nimages), right_norms(0, min_nimages), row_max_epss(0))
      ENDIF
      !
      CALL timeset(routineN//"_scatter", handle2)
      IF (right_npcols .GT. 1 .AND. right_col_nimages .GT. 1 .AND. .NOT. has_acc) THEN
         CALL mp_waitall(requests_scatter)
         IF (keep_sparsity) THEN
            CALL mp_wait(request_product_matrix)
         ENDIF
      ENDIF
      !
      ! Prepare result matrix for scatter data
      nimages_scatter = 0
      IF (do_scatter .OR. do_diagonal) THEN
         CALL dbcsr_init(product_matrix_scatter)
         !
         ! Check scattered meta
         IF (local_meta_size_scatter .GT. 1) THEN
            ALLOCATE (v_k_scatter((right_row_nimages+right_npcols-1)/right_npcols))
            meta_size_scatter = 0
            DO WHILE (meta_size_scatter+dbcsr_num_slots .LE. local_meta_size_scatter)
               nimages_scatter = nimages_scatter+1
               v_k_scatter(nimages_scatter) = local_meta_scatter(meta_size_scatter+dbcsr_slot_home_rowi)
               meta_size_scatter = meta_size_scatter+local_meta_scatter(meta_size_scatter+dbcsr_slot_size)
            ENDDO
         ENDIF
      ENDIF
      IF (do_scatter) THEN
         IF (.NOT. dbcsr_data_valid(local_data_product_scatter)) THEN
            CALL dbcsr_data_init(local_data_product_scatter)
            CALL dbcsr_data_new(local_data_product_scatter, data_type=data_type, &
                                data_size=1, memory_type=memtype_mpi_buffer)
         ENDIF
         CALL dbcsr_create(product_matrix_scatter, &
                           template=product_matrix, &
                           name=TRIM("Scatter buffer of "//TRIM(product_matrix%m%name)), &
                           make_index=.NOT. keep_sparsity .OR. local_meta_size_scatter .LE. 1, &
                           data_buffer=local_data_product_scatter)
         IF (local_meta_size_scatter .GT. 1) THEN
            max_nblocks = MAX(max_nblocks, meta_size_scatter/3)
            !
            IF (keep_sparsity) THEN
               product_matrix_scatter%m%index => local_meta_product_scatter
               CALL dbcsr_data_ensure_size(product_matrix_scatter%m%data_area, &
                                           product_matrix_scatter%m%index(dbcsr_slot_nze), &
                                           zero_pad=.TRUE.)
            ELSE
               CALL dbcsr_data_set_size_referenced(product_matrix_scatter%m%data_area, 0)
            ENDIF
            !
            product_matrix_scatter%m%index(dbcsr_slot_home_pcol) = myprow
            CALL dbcsr_reset_locals(product_matrix_scatter)
            CALL dbcsr_repoint_index(product_matrix_scatter%m)
            product_matrix_scatter%m%valid = .FALSE.
            ALLOCATE (multrec_scatter(0:nthreads-1))
            !
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          PRIVATE (size_guess, ithread) &
!$OMP          SHARED (product_matrix_scatter, multrec_scatter, &
!$OMP                  keep_sparsity, filter_eps, &
!$OMP                  row_max_epss, max_nblocks, &
!$OMP                  matrix_right, nthreads, &
!$OMP                  memtype_product_wm_scatter, &
!$OMP                  local_meta_size_scatter, &
!$OMP                  keep_product_data)
            !
            ! Setup product work areas
            !
            ithread = 0
!$          ithread = OMP_GET_THREAD_NUM()
            !
            IF (keep_product_data) THEN
               CALL dbcsr_add_wm_from_matrix(product_matrix_scatter)
            ELSE
               CALL dbcsr_work_create(product_matrix_scatter, work_mutable=.FALSE., &
                                      memory_type=memtype_product_wm_scatter(ithread)%p)
            ENDIF
!$OMP BARRIER
            ! The work arrays have to be setup
            size_guess = product_matrix_scatter%m%wms(ithread+1)%datasize ! Should be minimal
            CALL dbcsr_data_ensure_size(product_matrix_scatter%m%wms(ithread+1)%data_area, &
                                        size_guess)
            CALL dbcsr_data_set_size_referenced(product_matrix_scatter%m%wms(ithread+1)%data_area, &
                                                product_matrix_scatter%m%wms(ithread+1)%datasize)
            CALL ensure_array_size(product_matrix_scatter%m%wms(ithread+1)%row_i, ub=1)
            CALL ensure_array_size(product_matrix_scatter%m%wms(ithread+1)%col_i, ub=1)
            CALL ensure_array_size(product_matrix_scatter%m%wms(ithread+1)%blk_p, ub=1)
            ALLOCATE (multrec_scatter(ithread)%p)
            CALL dbcsr_mm_multrec_init(multrec_scatter(ithread)%p, &
                                       product=product_matrix_scatter%m, &
                                       keep_sparsity=keep_sparsity, &
                                       eps=filter_eps, &
                                       row_max_epss=row_max_epss, &
                                       block_estimate=max_nblocks/nthreads, &
                                       right_row_blk_size=dbcsr_row_block_sizes(matrix_right))
!$OMP END PARALLEL
            product_matrix_scatter%m%nblks = 0
            product_matrix_scatter%m%nze = 0
            product_matrix_scatter%m%row_p(:) = 0
            !
            ! Set images with scatter data
            CALL setup_buffer_matrices_images(buffer_scatter, nimages_scatter, &
                                              imgdist_right, matrix_right, &
                                              data_buffer=local_data_scatter)
            meta_size_scatter = 0
            right_data_size = 0
            DO v_ki = 1, nimages_scatter
               buffer_scatter%mats(v_ki)%m%index => local_meta_scatter(meta_size_scatter+dbcsr_slot_size: &
                                                                       meta_size_scatter+ &
                                                                       local_meta_scatter(meta_size_scatter+ &
                                                                                          dbcsr_slot_size))
               meta_size_scatter = meta_size_scatter+local_meta_scatter(meta_size_scatter+dbcsr_slot_size)
               CALL dbcsr_reset_vlocals(buffer_scatter%mats(v_ki), imgdist_right, do_rows=.TRUE.)
               ! Set local cols
               CALL array_release(buffer_scatter%mats(v_ki)%m%local_cols)
               buffer_scatter%mats(v_ki)%m%local_cols = product_matrix_scatter%m%local_cols
               CALL array_hold(buffer_scatter%mats(v_ki)%m%local_cols)
               buffer_scatter%mats(v_ki)%m%has_local_cols = .TRUE.
               buffer_scatter%mats(v_ki)%m%nblkcols_local = product_matrix_scatter%m%nblkcols_local
               buffer_scatter%mats(v_ki)%m%index(dbcsr_slot_nblkcols_local) = &
                  product_matrix_scatter%m%nblkcols_local
               CALL dbcsr_repoint_index(buffer_scatter%mats(v_ki)%m)
               buffer_scatter%mats(v_ki)%m%index(dbcsr_num_slots+3: &
                                                 buffer_scatter%mats(v_ki)%m%index(dbcsr_slot_size):3) = &
                  buffer_scatter%mats(v_ki)%m%index(dbcsr_num_slots+3: &
                                                    buffer_scatter%mats(v_ki)%m%index(dbcsr_slot_size):3)+ &
                  right_data_size
               right_data_size = right_data_size+buffer_scatter%mats(v_ki)%m%index(dbcsr_slot_nze)
               buffer_scatter%mats(v_ki)%m%valid = .TRUE.
            ENDDO
         ENDIF
      ENDIF
      CALL timestop(handle2)
      !
      ! Set RMA windows for the original data
      CALL dbcsr_win_create_any(rma_buffers_orig%left%data, grp_left, left_data_win)
      CALL mp_win_create(rma_buffers_orig%left%meta, grp_left, left_meta_win)
      CALL dbcsr_win_create_any(rma_buffers_orig%right%data, grp_right, right_data_win)
      CALL mp_win_create(rma_buffers_orig%right%meta, grp_right, right_meta_win)
      CALL mp_win_lock_all(left_data_win)
      CALL mp_win_lock_all(left_meta_win)
      CALL mp_win_lock_all(right_data_win)
      CALL mp_win_lock_all(right_meta_win)
      !
      ALLOCATE (multrec(0:nthreads-1))
      !
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          PRIVATE (size_guess, ithread) &
!$OMP          SHARED (product_matrix, multrec, &
!$OMP                  keep_sparsity, filter_eps, &
!$OMP                  row_max_epss, max_nblocks, &
!$OMP                  matrix_right, nthreads)
      !
      ! Setup product work areas
      !
      ithread = 0
!$    ithread = OMP_GET_THREAD_NUM()
      !
      ! The work arrays have to be setup
      size_guess = product_matrix%m%wms(ithread+1)%datasize ! Should be minimal
      CALL dbcsr_data_ensure_size(product_matrix%m%wms(ithread+1)%data_area, &
                                  size_guess)
      CALL dbcsr_data_set_size_referenced(product_matrix%m%wms(ithread+1)%data_area, &
                                          product_matrix%m%wms(ithread+1)%datasize)
      CALL ensure_array_size(product_matrix%m%wms(ithread+1)%row_i, ub=1)
      CALL ensure_array_size(product_matrix%m%wms(ithread+1)%col_i, ub=1)
      CALL ensure_array_size(product_matrix%m%wms(ithread+1)%blk_p, ub=1)
      ALLOCATE (multrec(ithread)%p)
      CALL dbcsr_mm_multrec_init(multrec(ithread)%p, &
                                 product=product_matrix%m, &
                                 keep_sparsity=keep_sparsity, &
                                 eps=filter_eps, &
                                 row_max_epss=row_max_epss, &
                                 block_estimate=max_nblocks/nthreads, &
                                 right_row_blk_size=dbcsr_row_block_sizes(matrix_right))
!$OMP END PARALLEL
      !
      ! update capacity of memory-pools
      IF (has_acc) THEN
         CALL dbcsr_mempool_ensure_capacity(memtype_abpanel_1%pool, &
                                            capacity=left_col_nimages+right_row_nimages)
         CALL dbcsr_mempool_ensure_capacity(memtype_abpanel_2%pool, &
                                            capacity=left_col_nimages+right_row_nimages)
      ENDIF
      !
      ! Wait refs and max norms (sent in make_buffers)
      CALL mp_waitall(requests)
      ! Delete buffers used for collectives
      DEALLOCATE (right_local_refs_displ_unmerged, left_local_refs_displ_unmerged)
      IF (otf_filtering .AND. use_mpi_filtering) THEN
         DEALLOCATE (right_local_max_norms, left_local_max_norms)
      ENDIF
      DEALLOCATE (right_local_refs_meta_size, left_local_refs_meta_size)
      DEALLOCATE (right_local_refs_data_size, left_local_refs_data_size)
      DEALLOCATE (right_no_empty_images, left_no_empty_images)
      !
      ! Pre-CANNON execution: filtering and resize buffers for communications
      CALL timeset(routineN//"_precannon", handle2)
      !
      ALLOCATE (right_vunmerged(0:min_nimages-1, 0:nsteps_k-1))
      ALLOCATE (left_vunmerged(0:min_nimages-1, 0:nsteps_k-1))
      right_vunmerged = 0; left_vunmerged = 0
      !
      ALLOCATE (right_displ_vunmerged(0:nvirt_k-1))
      ALLOCATE (left_displ_vunmerged(0:nvirt_k-1))
      !
      ALLOCATE (right_nofiltered(right_col_nimages, 0:min_nimages-1, 0:nsteps_k-1))
      ALLOCATE (left_nofiltered(left_row_nimages, 0:min_nimages-1, 0:nsteps_k-1))
      right_nofiltered = .FALSE.; left_nofiltered = .FALSE.
      ALLOCATE (right_vmerged(right_col_nimages, 0:min_nimages-1, 0:nsteps_k-1))
      ALLOCATE (left_vmerged(left_row_nimages, 0:min_nimages-1, 0:nsteps_k-1))
      !
      ALLOCATE (right_recv_filtered(right_col_nimages, 0:min_nimages-1))
      ALLOCATE (left_recv_filtered(left_row_nimages, 0:min_nimages-1))
      !
      right_max_data_size = 0; right_max_meta_size = 0
      left_max_data_size = 0; left_max_meta_size = 0
      !
      ALLOCATE (left_max_meta_size_merged(0:min_nimages-1))
      left_max_meta_size_merged = 0
      ALLOCATE (right_max_meta_size_merged(0:min_nimages-1))
      right_max_meta_size_merged = 0
      !
!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED (nsteps_k,min_nimages,right_col_nimages,imgdist_right,&
!$OMP         otf_filtering,use_mpi_filtering,right_nprows,left_npcols,&
!$OMP         right_myfirstvrow,right_myfirstvcol,right_max_norms,&
!$OMP         left_row_nimages,imgdist_left,left_myfirstvrow,left_myfirstvcol,&
!$OMP         left_max_norms,min_row_max_epss,right_row_nimages,&
!$OMP         left_recv_filtered,right_recv_filtered,right_refs_meta_size,&
!$OMP         left_refs_meta_size,right_refs_data_size,left_refs_data_size,&
!$OMP         right_max_data_size,left_max_data_size,left_col_nimages,&
!$OMP         right_max_meta_size,left_max_meta_size,do_diagonal,&
!$OMP         right_data_size,left_data_size,right_meta_size,left_meta_size,&
!$OMP         size_index_unmerged,left_max_meta_size_merged,right_max_meta_size_merged,&
!$OMP         right_vunmerged,left_vunmerged,right_nofiltered,left_nofiltered,myprow,&
!$OMP         right_vmerged,left_vmerged,left_displ_vunmerged,right_displ_vunmerged,&
!$OMP         right_no_empty_images_displ,left_no_empty_images_displ,do_scatter,&
!$OMP         local_meta_size_scatter,v_k_scatter,nimages_scatter) &
!$OMP PRIVATE (metronome,v_ki,v_ci,v_ri,recv_vrow,recv_vcol,right_vrow,&
!$OMP          right_vcol,left_vrow,left_vcol,meta_size_image,iproc,v_k_loc,&
!$OMP          do_k_scatter)
      !
!$OMP DO
      DO iproc = 0, right_nprows-1
         recv_vrow = iproc*right_row_nimages
         right_displ_vunmerged(recv_vrow) = right_no_empty_images_displ(iproc+1)
         DO v_ki = recv_vrow+1, recv_vrow+right_row_nimages-1
            right_displ_vunmerged(v_ki) = right_displ_vunmerged(v_ki-1)
            DO v_ci = LBOUND(right_refs_meta_size, 1), UBOUND(right_refs_meta_size, 1)
               IF (right_refs_meta_size(v_ci, v_ki-1) .NE. 0) THEN
                  right_displ_vunmerged(v_ki) = right_displ_vunmerged(v_ki)+1
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!$OMP END DO
!$OMP DO
      DO iproc = 0, left_npcols-1
         recv_vcol = iproc*left_col_nimages
         left_displ_vunmerged(recv_vcol) = left_no_empty_images_displ(iproc+1)
         DO v_ki = recv_vcol+1, recv_vcol+left_col_nimages-1
            left_displ_vunmerged(v_ki) = left_displ_vunmerged(v_ki-1)
            DO v_ri = LBOUND(left_refs_meta_size, 1), UBOUND(left_refs_meta_size, 1)
               IF (left_refs_meta_size(v_ri, v_ki-1) .NE. 0) THEN
                  left_displ_vunmerged(v_ki) = left_displ_vunmerged(v_ki)+1
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!$OMP END DO
      !
      DO metronome = 0, nsteps_k-1
         !
!$OMP MASTER
         right_data_size = 0
         left_data_size = 0
         right_meta_size = 0
         left_meta_size = 0
!$OMP END MASTER
!$OMP BARRIER
         !
         ! Take first cluster global virtual coordinates
         CALL image_calculator(imgdist_right, &
                               vprow=right_vrow, vpcol=right_vcol, &
                               myvprow=right_myfirstvrow, &
                               myvpcol=right_myfirstvcol, &
                               vprow_shift=metronome*min_nimages, &
                               shifting='R')
         CALL image_calculator(imgdist_left, &
                               vprow=left_vrow, vpcol=left_vcol, &
                               myvprow=left_myfirstvrow, &
                               myvpcol=left_myfirstvcol, &
                               vpcol_shift=metronome*min_nimages, &
                               shifting='L')
         !
         v_k_loc = 0
!$OMP DO REDUCTION(+:right_data_size,left_data_size,right_meta_size,left_meta_size) &
!$OMP    schedule(dynamic)
         DO v_ki = 0, min_nimages-1
            !
            recv_vrow = v_ki+right_vrow
            recv_vcol = v_ki+left_vcol
            !
            DO v_ci = right_vcol, right_col_nimages+right_vcol-1
               IF (right_refs_meta_size(v_ci, recv_vrow) .NE. 0) THEN
                  right_vunmerged(v_ki, metronome) = right_vunmerged(v_ki, metronome)+1
                  right_vmerged(right_vunmerged(v_ki, metronome), v_ki, metronome) = v_ci
               ENDIF
            ENDDO
            DO v_ri = left_vrow, left_row_nimages+left_vrow-1
               IF (left_refs_meta_size(v_ri, recv_vcol) .NE. 0) THEN
                  left_vunmerged(v_ki, metronome) = left_vunmerged(v_ki, metronome)+1
                  left_vmerged(left_vunmerged(v_ki, metronome), v_ki, metronome) = v_ri
               ENDIF
            ENDDO
            !
            do_k_scatter = .FALSE.
            IF ((do_diagonal .OR. do_scatter) .AND. &
                local_meta_size_scatter .GT. 1 .AND. &
                recv_vcol/min_nimages .EQ. myprow) THEN
               CALL ordered_search(v_k_scatter, v_ki+1, v_k_loc, do_k_scatter, v_k_loc+1, nimages_scatter)
            ENDIF
            !
            ! Max norms evaluation for filtering
            !
            ! By default all images are not communicated
            left_recv_filtered(:, v_ki) = 0
            right_recv_filtered(:, v_ki) = 0
            !
            ! multiplication of norms
            DO v_ri = 1, left_vunmerged(v_ki, metronome)
               DO v_ci = 1, right_vunmerged(v_ki, metronome)
                  IF (do_diagonal .AND. &
                      right_vmerged(v_ci, v_ki, metronome) .EQ. recv_vrow .AND. &
                      .NOT. do_k_scatter) CYCLE
                  IF (otf_filtering .AND. use_mpi_filtering) THEN
                     IF (left_max_norms(left_displ_vunmerged(recv_vcol)+v_ri)* &
                         right_max_norms(right_displ_vunmerged(recv_vrow)+v_ci) .LT. &
                         min_row_max_epss) CYCLE
                  ENDIF
                  left_recv_filtered(v_ri, v_ki) = left_recv_filtered(v_ri, v_ki)+1
                  right_recv_filtered(v_ci, v_ki) = right_recv_filtered(v_ci, v_ki)+1
               ENDDO
            ENDDO
            !
            ! resize buffers taking in account filtering
            meta_size_image = 0
            DO v_ci = 1, right_vunmerged(v_ki, metronome)
               IF (right_recv_filtered(v_ci, v_ki) .NE. 0) THEN
                  right_nofiltered(v_ci, v_ki, metronome) = .TRUE.
                  right_data_size = right_data_size+ &
                                    right_refs_data_size(right_displ_vunmerged(recv_vrow)+v_ci)
                  right_meta_size = right_meta_size+ &
                                    right_refs_meta_size(right_vmerged(v_ci, v_ki, metronome), recv_vrow)
                  meta_size_image = meta_size_image+ &
                                    right_refs_meta_size(right_vmerged(v_ci, v_ki, metronome), recv_vrow)- &
                                    dbcsr_slot_nblks
               ENDIF
            ENDDO
            right_max_meta_size_merged(v_ki) = MAX(right_max_meta_size_merged(v_ki), meta_size_image)
            !
            meta_size_image = 0
            DO v_ri = 1, left_vunmerged(v_ki, metronome)
               IF (left_recv_filtered(v_ri, v_ki) .NE. 0 .OR. do_k_scatter) THEN
                  left_nofiltered(v_ri, v_ki, metronome) = .TRUE.
                  left_data_size = left_data_size+ &
                                   left_refs_data_size(left_displ_vunmerged(recv_vcol)+v_ri)
                  left_meta_size = left_meta_size+ &
                                   left_refs_meta_size(left_vmerged(v_ri, v_ki, metronome), recv_vcol)
                  meta_size_image = meta_size_image+ &
                                    left_refs_meta_size(left_vmerged(v_ri, v_ki, metronome), recv_vcol)- &
                                    size_index_unmerged
               ENDIF
            ENDDO
            left_max_meta_size_merged(v_ki) = MAX(left_max_meta_size_merged(v_ki), meta_size_image)
         ENDDO
!$OMP END DO
!$OMP MASTER
         right_max_data_size = MAX(right_max_data_size, right_data_size)
         left_max_data_size = MAX(left_max_data_size, left_data_size)
         right_max_meta_size = MAX(right_max_meta_size, right_meta_size)
         left_max_meta_size = MAX(left_max_meta_size, left_meta_size)
!$OMP END MASTER
!$OMP BARRIER
      ENDDO
!$OMP END PARALLEL
      IF (otf_filtering .AND. use_mpi_filtering) THEN
         DEALLOCATE (right_max_norms, left_max_norms)
      ENDIF
      DEALLOCATE (left_recv_filtered, right_recv_filtered)
      !
      CALL timestop(handle2)
      !
      ! Preallocate norms arrays
      IF (otf_filtering) THEN
         ALLOCATE (left_norms(MAXVAL(left_max_meta_size_merged)/3, min_nimages))
         ALLOCATE (right_norms(MAX(MAXVAL(right_max_meta_size_merged), &
                                   local_meta_size_scatter)/3, min_nimages))
      ENDIF
      left_max_meta_size_merged = left_max_meta_size_merged+size_index_merged
      right_max_meta_size_merged = right_max_meta_size_merged+dbcsr_num_slots
      !
      ! Prepare buffers for computation
      IF (nsteps_k .GT. 1) THEN
         ! Right:
         CALL rma_buffer_init(rma_buffers_1%right, data_type, &
                              right_max_data_size, &
                              right_max_meta_size, &
                              requests_size=(min_nimages*right_col_nimages), &
                              data_memory_type=memtype_abpanel_1)
         CALL setup_buffer_matrices_images(rma_buffers_1%right%buffer, min_nimages, &
                                           imgdist_right, matrix_right, &
                                           right_max_meta_size_merged, &
                                           rma_buffers_1%right%data)
         ! Left:
         CALL rma_buffer_init(rma_buffers_1%left, data_type, &
                              left_max_data_size, &
                              left_max_meta_size, &
                              requests_size=(left_row_nimages*min_nimages), &
                              data_memory_type=memtype_abpanel_1)
         CALL setup_buffer_matrices_images(rma_buffers_1%left%buffer, min_nimages, &
                                           imgdist_left, matrix_left, &
                                           left_max_meta_size_merged, &
                                           rma_buffers_1%left%data)
         !
         IF (has_acc) THEN
            CALL dbcsr_data_init(rma_buffers_1%trs_stackbuf)
            CALL dbcsr_data_new(rma_buffers_1%trs_stackbuf, &
                                data_type=dbcsr_type_int_4, data_size=1000, &
                                memory_type=memtype_trsbuffer_1)
         ENDIF
         !
         rma_buffers_calc => rma_buffers_1
      ENDIF
      !
      ! Prepare buffers for communication
      ! Right:
      CALL rma_buffer_init(rma_buffers_2%right, data_type, &
                           right_max_data_size, &
                           right_max_meta_size, &
                           requests_size=(min_nimages*right_col_nimages), &
                           data_memory_type=memtype_abpanel_2)
      CALL setup_buffer_matrices_images(rma_buffers_2%right%buffer, min_nimages, &
                                        imgdist_right, matrix_right, &
                                        right_max_meta_size_merged, &
                                        rma_buffers_2%right%data)
      ! Left:
      CALL rma_buffer_init(rma_buffers_2%left, data_type, &
                           left_max_data_size, &
                           left_max_meta_size, &
                           requests_size=(left_row_nimages*min_nimages), &
                           data_memory_type=memtype_abpanel_2)
      CALL setup_buffer_matrices_images(rma_buffers_2%left%buffer, min_nimages, &
                                        imgdist_left, matrix_left, &
                                        left_max_meta_size_merged, &
                                        rma_buffers_2%left%data)
      !
      IF (has_acc) THEN
         CALL dbcsr_data_init(rma_buffers_2%trs_stackbuf)
         CALL dbcsr_data_new(rma_buffers_2%trs_stackbuf, &
                             data_type=dbcsr_type_int_4, data_size=1000, &
                             memory_type=memtype_trsbuffer_2)
      ENDIF
      !
      rma_buffers_comm => rma_buffers_2
      !
      DEALLOCATE (right_max_meta_size_merged, left_max_meta_size_merged)
      !
      ! Setup the receive data pointers
      CALL dbcsr_data_init(data_get)
      CALL dbcsr_data_new(data_get, data_type)
      !
      ! These values for meta data are used for global values
      right_global_indices(dbcsr_slot_nblkrows_total:dbcsr_slot_nfullcols_local) = &
         (/ &
         matrix_right%m%nblkrows_total, &
         matrix_right%m%nblkcols_total, &
         matrix_right%m%nfullrows_total, &
         matrix_right%m%nfullcols_total, &
         0, 0, &
         matrix_right%m%nfullrows_local, &
         matrix_right%m%nfullcols_local/)
      left_global_indices(dbcsr_slot_nblkrows_total:dbcsr_slot_nfullcols_local) = &
         (/ &
         matrix_left%m%nblkrows_total, &
         matrix_left%m%nblkcols_total, &
         matrix_left%m%nfullrows_total, &
         matrix_left%m%nfullcols_total, &
         0, 0, &
         matrix_left%m%nfullrows_local, &
         matrix_left%m%nfullcols_local/)
      !
      ALLOCATE (right_get_requests_images(min_nimages))
      ALLOCATE (left_get_requests_images(min_nimages))
      ALLOCATE (common_unmerged(min_nimages))
      ALLOCATE (indices_common_unmerged(min_nimages))
      !
      ! Wait for diagonal images data
      CALL timeset(routineN//"_diag", handle2)
      CALL mp_waitall(requests_diag)
      CALL timestop(handle2)
      !
      ! Here is the main loop
      !
      CALL timeset(routineN//"_loop", handle2)
      !
      right_meta_offset_diag = 0
      right_data_offset_diag = 0
      left_meta_offset_diag = 0
      left_data_offset_diag = 0
      !
      grouped_k_index: DO metronome = 0, nsteps_k
         !
         IF (debug_mod) WRITE (*, '(1X,A,3(1X,A,1X,I5))') routineN, &
            "step", metronome, &
            "first k", metronome*min_nimages, &
            "last k", (metronome+1)*min_nimages-1
         !
         ! Matrix transfer. Transfer in all but the last loop
         ! iteration.
         xfer: IF (metronome .LT. nsteps_k) THEN
            !
            IF (has_acc) THEN
               CALL timeset(routineN//"_acc_sync", handle3)
               CALL acc_event_synchronize(rma_buffers_comm%right%data%d%acc_ready)
               CALL acc_event_synchronize(rma_buffers_comm%left%data%d%acc_ready)
               CALL timestop(handle3)
            ENDIF
            !
            right_meta_offset = 0
            right_meta_size = 0
            right_data_offset = 0
            right_data_size = 0
            rma_buffers_comm%right%nrequests = 0
            !
            ! Take first cluster global virtual coordinates
            CALL image_calculator(imgdist_right, &
                                  vprow=rma_buffers_comm%right%vprow, &
                                  vpcol=rma_buffers_comm%right%vpcol, &
                                  myvprow=right_myfirstvrow, &
                                  myvpcol=right_myfirstvcol, &
                                  vprow_shift=metronome*min_nimages, &
                                  shifting='R')
            !
            left_meta_offset = 0
            left_meta_size = 0
            left_data_offset = 0
            left_data_size = 0
            rma_buffers_comm%left%nrequests = 0
            !
            ! Take first cluster global virtual coordinates
            CALL image_calculator(imgdist_left, &
                                  vprow=rma_buffers_comm%left%vprow, &
                                  vpcol=rma_buffers_comm%left%vpcol, &
                                  myvprow=left_myfirstvrow, &
                                  myvpcol=left_myfirstvcol, &
                                  vpcol_shift=metronome*min_nimages, &
                                  shifting='L')
            !
            right_set_displ = .TRUE.
            left_set_displ = .TRUE.
            right_v_i = 0
            left_v_i = 0
            rma_buffers_comm%right%offset(:, :) = 0
            rma_buffers_comm%left%offset(:, :) = 0
            rma_buffers_comm%right%get_requests(:, :) = mp_request_null
            rma_buffers_comm%left%get_requests(:, :) = mp_request_null
            !
            DO v_ki = 0, min_nimages-1
               ! Right
               CALL cannon_transfer_rma(v_ki, rma_buffers_comm%right%vprow, rma_buffers_comm%right%vpcol, &
                                        right_row_nimages, right_col_nimages, &
                                        right_vunmerged(:, metronome), right_vmerged(:, :, metronome), &
                                        right_displ_vunmerged, &
                                        right_v_i, right_nofiltered(:, :, metronome), right_set_displ, &
                                        recv_vrow_displ, right_refs_displ_unmerged, &
                                        rma_buffers_comm%right, right_refs_meta_size, right_refs_data_size, &
                                        right_meta_size, right_data_size, &
                                        min_nimages, right_meta_offset, right_data_offset, &
                                        right_meta_win, right_data_win, right_myprow, metronome .GT. 0, &
                                        data_get, data_type_byte, rma_buffers_orig%right, &
                                        right_meta_offset_diag, right_data_offset_diag)
               ! Left
               CALL cannon_transfer_rma(v_ki, rma_buffers_comm%left%vpcol, rma_buffers_comm%left%vprow, &
                                        left_col_nimages, left_row_nimages, &
                                        left_vunmerged(:, metronome), left_vmerged(:, :, metronome), &
                                        left_displ_vunmerged, &
                                        left_v_i, left_nofiltered(:, :, metronome), left_set_displ, &
                                        recv_vcol_displ, left_refs_displ_unmerged, &
                                        rma_buffers_comm%left, left_refs_meta_size, left_refs_data_size, &
                                        left_meta_size, left_data_size, &
                                        min_nimages, left_meta_offset, left_data_offset, &
                                        left_meta_win, left_data_win, left_mypcol, metronome .GT. 0, &
                                        data_get, data_type_byte, rma_buffers_orig%left, &
                                        left_meta_offset_diag, left_data_offset_diag)
            ENDDO
            ! Set the referenced sizes to the actual data moved via MPI
            CALL dbcsr_data_set_size_referenced(rma_buffers_comm%right%data, &
                                                rma_buffers_comm%right%offset(idata_size, right_v_i+1))
            CALL dbcsr_data_set_size_referenced(rma_buffers_comm%left%data, &
                                                rma_buffers_comm%left%offset(idata_size, left_v_i+1))
         ENDIF xfer
         !
         ! Wait data, merge and do the multiplications.
         ! Exclude the first interation
         wait_merge_calc: IF (metronome .GT. 0) THEN
            IF (debug_mod) WRITE (*, '(1X,A)') routineN//" waiting for right and left"
            !
            rma_buffers_calc%right%buffer%mats(:)%m%nblks = 0
            rma_buffers_calc%left%buffer%mats(:)%m%nblks = 0
            !
            common_unmerged(:) = .FALSE.
            !
            IF (metronome == nsteps_k) THEN
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (multrec) &
!$OMP PRIVATE (ithread)
               ithread = 0
!$             ithread = omp_get_thread_num()
               !
               CALL dbcsr_mm_multrec_phaseout(multrec(ithread)%p)
!$OMP END PARALLEL
            ENDIF
            !
            right_nrequests = rma_buffers_calc%right%nrequests
            left_nrequests = rma_buffers_calc%left%nrequests
            IF (min_nimages .EQ. 1 .OR. has_acc) THEN
               IF (right_nrequests .GT. 0) THEN
                  CALL mp_waitall(rma_buffers_calc%right%get_requests(:, 1:right_nrequests))
                  IF (has_acc) CALL dbcsr_data_host2dev(rma_buffers_calc%right%data)
                  CALL merge_calc(rma_buffers_calc%right, &
                                  rma_buffers_calc%left, &
                                  multrec, do_rec, flop, &
                                  imgdist_right, .FALSE., &
                                  1, min_nimages, &
                                  min_nimages, right_col_nimages, &
                                  right_global_indices, right_vunmerged(:, metronome-1), &
                                  right_vmerged(:, :, metronome-1), otf_filtering, data_type, &
                                  right_norms, left_norms, &
                                  common_unmerged, indices_common_unmerged, has_acc)
               ENDIF
               IF (left_nrequests .GT. 0) THEN
                  CALL mp_waitall(rma_buffers_calc%left%get_requests(:, 1:left_nrequests))
                  IF (has_acc) CALL dbcsr_data_host2dev(rma_buffers_calc%left%data)
                  CALL merge_calc(rma_buffers_calc%right, &
                                  rma_buffers_calc%left, &
                                  multrec, do_rec, flop, &
                                  imgdist_left, .TRUE., &
                                  1, min_nimages, &
                                  left_row_nimages, min_nimages, &
                                  left_global_indices, left_vunmerged(:, metronome-1), &
                                  left_vmerged(:, :, metronome-1), otf_filtering, data_type, &
                                  right_norms, left_norms, &
                                  common_unmerged, indices_common_unmerged, &
                                  has_acc, nthreads, rma_buffers_calc%trs_stackbuf)
               ENDIF
            ELSE
               !
               ! Use async communication-computation over images
               CALL timeset("mp_waitall_2", handle3)
               !
               ! Counting the requests per each unmerged image
               right_get_requests_images(:) = 0
               DO iget = 1, rma_buffers_calc%right%nrequests
                  right_get_requests_images(rma_buffers_calc%right%get_requests_map(:, iget)) = &
                     right_get_requests_images(rma_buffers_calc%right%get_requests_map(:, iget))+1
               ENDDO
               left_get_requests_images(:) = 0
               DO iget = 1, rma_buffers_calc%left%nrequests
                  left_get_requests_images(rma_buffers_calc%left%get_requests_map(:, iget)) = &
                     left_get_requests_images(rma_buffers_calc%left%get_requests_map(:, iget))+1
               ENDDO
               !
               DO WHILE (right_nrequests .GT. 0 .OR. &
                         left_nrequests .GT. 0)
                  IF (right_nrequests .GT. 0) THEN
                     CALL wait_merge_calc_low(rma_buffers_calc%right, &
                                              rma_buffers_calc%left, &
                                              right_nrequests, &
                                              multrec, do_rec, flop, &
                                              right_get_requests_images, &
                                              imgdist_right, .FALSE., min_nimages, right_col_nimages, &
                                              right_global_indices, right_vunmerged(:, metronome-1), &
                                              right_vmerged(:, :, metronome-1), otf_filtering, data_type, &
                                              right_norms, left_norms, &
                                              common_unmerged, indices_common_unmerged)
                  ENDIF
                  IF (left_nrequests .GT. 0) THEN
                     CALL wait_merge_calc_low(rma_buffers_calc%right, &
                                              rma_buffers_calc%left, &
                                              left_nrequests, &
                                              multrec, do_rec, flop, &
                                              left_get_requests_images, &
                                              imgdist_left, .TRUE., left_row_nimages, min_nimages, &
                                              left_global_indices, left_vunmerged(:, metronome-1), &
                                              left_vmerged(:, :, metronome-1), otf_filtering, data_type, &
                                              right_norms, left_norms, &
                                              common_unmerged, indices_common_unmerged, &
                                              nthreads)
                  ENDIF
               ENDDO
               !
               CALL timestop(handle3)
            ENDIF
            !
            IF (debug_mod) THEN
               DO v_ki = 1, min_nimages
                  IF (rma_buffers_calc%right%buffer%mats(v_ki)%m%nblks .NE. 0) &
                     CALL dbcsr_print(rma_buffers_calc%right%buffer%mats(v_ki), &
                                      nodata=.TRUE.)
                  IF (rma_buffers_calc%left%buffer%mats(v_ki)%m%nblks .NE. 0) &
                     CALL dbcsr_print(rma_buffers_calc%left%buffer%mats(v_ki), &
                                      nodata=.TRUE.)
               ENDDO
            ENDIF
         ENDIF wait_merge_calc
         !
         ! Multiply scatter data
         IF (do_scatter .AND. local_meta_size_scatter .GT. 1 .AND. metronome .GT. 0 .AND. &
             rma_buffers_calc%left%nrequests .GT. 0 .AND. &
             rma_buffers_calc%left%vpcol/min_nimages .EQ. myprow) THEN
            IF (otf_filtering) THEN
               CALL calculate_image_norms(buffer_scatter, &
                                          data_type, right_norms, &
                                          uf=1, ul=nimages_scatter)
            ENDIF
            CALL timeset("cannon_multiply_low_rma_multrec", handle3)
            !
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (buffer_scatter,rma_buffers_calc, &
!$OMP         left_norms,right_norms,nimages_scatter,&
!$OMP         multrec_scatter,v_k_scatter) &
!$OMP PRIVATE (ithread,v_ki) &
!$OMP REDUCTION (+: flop)
            ithread = 0
!$          ithread = omp_get_thread_num()
            !
            DO v_ki = 1, nimages_scatter
               IF (rma_buffers_calc%left%buffer%mats(v_k_scatter(v_ki))%m%nblks .EQ. 0) CYCLE
               CALL dbcsr_mm_multrec_multiply(multrec_scatter(ithread)%p, &
                                              left=rma_buffers_calc%left%buffer%mats(v_k_scatter(v_ki))%m, &
                                              right=buffer_scatter%mats(v_ki)%m, &
                                              flop=flop, &
                                              a_norms=left_norms(:, v_k_scatter(v_ki)), &
                                              b_norms=right_norms(:, v_ki), &
                                              do_rec=.FALSE.)
            ENDDO
!$OMP END PARALLEL
            CALL timestop(handle3)
         ENDIF
         !
         IF (metronome .LT. nsteps_k) THEN
            rma_buffers => rma_buffers_comm
            rma_buffers_comm => rma_buffers_calc
            rma_buffers_calc => rma_buffers
         ENDIF
         !
      ENDDO grouped_k_index
      !
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (multrec, multrec_scatter, do_scatter, &
!$OMP         local_meta_size_scatter) &
!$OMP PRIVATE (ithread)
      ithread = 0
!$    ithread = omp_get_thread_num()
      CALL dbcsr_mm_multrec_finalize(multrec(ithread)%p)
      DEALLOCATE (multrec(ithread)%p)
      IF (do_scatter .AND. local_meta_size_scatter .GT. 1) THEN
         CALL dbcsr_mm_multrec_finalize(multrec_scatter(ithread)%p)
         DEALLOCATE (multrec_scatter(ithread)%p)
      ENDIF
!$OMP END PARALLEL
      !
      CALL timestop(handle2)
      !
      CALL m_memory(mem)
      max_memory = MAX(max_memory, REAL(mem))
      !
      IF (debug_mod) THEN
         v_ki = 0
         DO blk = 1, SIZE(product_matrix%m%blk_p)
            v_ki = MAX(v_ki, ABS(product_matrix%m%blk_p(blk)))
         ENDDO
         WRITE (*, *) routineN//" Actual final size", &
            LOG(REAL(dbcsr_data_get_size(product_matrix%m%data_area)))/LOG(10.0), &
            LOG(REAL(v_ki))/LOG(10.0)
      ENDIF
      !
      IF (has_acc) THEN
         CALL dbcsr_data_release(rma_buffers_1%trs_stackbuf)
         CALL dbcsr_data_release(rma_buffers_2%trs_stackbuf)
      ENDIF
      !
      DEALLOCATE (right_norms, left_norms, row_max_epss)
      !
      CALL dbcsr_data_clear_pointer(data_get)
      CALL dbcsr_data_release(data_get)
      !
      ! clean-up of communication buffers
      IF (nsteps_k .GT. 1) THEN
         CALL dbcsr_destroy_array(rma_buffers_1%right%buffer)
         CALL dbcsr_destroy_array(rma_buffers_1%left%buffer)
      ENDIF
      !
      DEALLOCATE (right_get_requests_images, left_get_requests_images)
      DEALLOCATE (common_unmerged, indices_common_unmerged)
      DEALLOCATE (right_refs_meta_size, left_refs_meta_size)
      DEALLOCATE (right_refs_data_size, left_refs_data_size)
      DEALLOCATE (right_refs_displ_unmerged, left_refs_displ_unmerged)
      DEALLOCATE (multrec)
      DEALLOCATE (right_vunmerged, left_vunmerged)
      DEALLOCATE (right_nofiltered, left_nofiltered)
      DEALLOCATE (right_vmerged, left_vmerged)
      DEALLOCATE (right_displ_vunmerged, left_displ_vunmerged)
      DEALLOCATE (right_no_empty_images_displ, left_no_empty_images_displ)
      CALL dbcsr_destroy_array(rma_buffers_2%right%buffer)
      CALL dbcsr_destroy_array(rma_buffers_2%left%buffer)
      !
      CALL dbcsr_finalize(product_matrix)
      !
      ! clean-up wins
      CALL mp_win_unlock_all(right_data_win)
      CALL mp_win_unlock_all(right_meta_win)
      CALL mp_win_unlock_all(left_data_win)
      CALL mp_win_unlock_all(left_meta_win)
      CALL mp_win_free(right_data_win)
      CALL mp_win_free(right_meta_win)
      CALL mp_win_free(left_data_win)
      CALL mp_win_free(left_meta_win)
      !
      ! Gather data
      IF (do_scatter .OR. do_diagonal) THEN
         CALL timeset("cannon_multiply_low_rma_gather", handle2)
         local_meta_size_gather = 1
         local_data_size_gather = 1
         IF (local_meta_size_scatter .GT. 1) THEN
            DEALLOCATE (v_k_scatter)
            IF (do_scatter) THEN
               DO v_ki = 1, nimages_scatter
                  NULLIFY (buffer_scatter%mats(v_ki)%m%index)
               ENDDO
               CALL dbcsr_destroy_array(buffer_scatter)
               DEALLOCATE (multrec_scatter)
               CALL dbcsr_finalize(product_matrix_scatter)
               IF (product_matrix_scatter%m%nblks .NE. 0) THEN
                  local_meta_size_gather = product_matrix_scatter%m%index(dbcsr_slot_size)
                  local_data_size_gather = product_matrix_scatter%m%index(dbcsr_slot_nze)
               ENDIF
            ENDIF
         ENDIF
         IF (do_diagonal) THEN
            ! Gather sizes
            ALLOCATE (meta_size_gather(0:right_npcols-1), &
                      data_size_gather(0:right_npcols-1), &
                      meta_displ_gather(0:right_npcols-1), &
                      data_displ_gather(0:right_npcols-1))
         ENDIF
         ! Gather sizes
         CALL mp_gather(local_meta_size_gather, meta_size_gather, &
                        right_myprow, grp_gather)
         CALL mp_gather(local_data_size_gather, data_size_gather, &
                        right_myprow, grp_gather)
         IF (do_diagonal) THEN
            meta_displ_gather(0) = 0; data_displ_gather(0) = 0
            DO iproc = 0, right_npcols-2
               meta_displ_gather(iproc+1) = meta_displ_gather(iproc)+ &
                                            meta_size_gather(iproc)
               data_displ_gather(iproc+1) = data_displ_gather(iproc)+ &
                                            data_size_gather(iproc)
            ENDDO
            ! Prepare result matrix
            ! Reuse scatter buffers
            CALL ensure_array_size(local_meta_scatter, &
                                   ub=meta_displ_gather(right_npcols-1)+meta_size_gather(right_npcols-1), &
                                   nocopy=.TRUE., memory_type=memtype_mpi_buffer)
            CALL dbcsr_data_ensure_size(local_data_scatter, &
                                        data_displ_gather(right_npcols-1)+data_size_gather(right_npcols-1), &
                                        nocopy=.TRUE.)
            CALL dbcsr_data_set_size_referenced(local_data_scatter, &
                                                data_displ_gather(right_npcols-1)+data_size_gather(right_npcols-1))
            CALL dbcsr_create(product_matrix_scatter, &
                              template=product_matrix, &
                              name=TRIM("Gather buffer of "//TRIM(product_matrix%m%name)), &
                              make_index=.NOT. keep_sparsity, &
                              data_buffer=local_data_scatter)
            IF (keep_sparsity) THEN
               product_matrix_scatter%m%index => local_meta_product_scatter
               CALL dbcsr_repoint_index(product_matrix_scatter%m)
            ENDIF
            CALL dbcsr_reset_locals(product_matrix_scatter)
            product_matrix_scatter%m%nblks = 0
            product_matrix_scatter%m%nze = 0
            product_matrix_scatter%m%row_p(:) = 0
            product_matrix_scatter%m%valid = .FALSE.
         ENDIF
         CALL mp_gatherv(product_matrix_scatter%m%index(1:local_meta_size_gather), &
                         local_meta_scatter, meta_size_gather, meta_displ_gather, &
                         right_myprow, grp_gather)
         IF (do_diagonal) THEN
            ! Receiver
            ! Use a dummy buffer for the sender of the root proc to avoid overlap
            ! between send and receive buffers
            CALL dbcsr_gatherv_any(rma_buffers_orig%right%data_diag, &
                                   local_data_size_gather, &
                                   product_matrix_scatter%m%data_area, &
                                   data_size_gather, data_displ_gather, &
                                   right_myprow, grp_gather)
         ENDIF
         IF (do_scatter) THEN
            ! Sender
            CALL dbcsr_gatherv_any(product_matrix_scatter%m%data_area, &
                                   local_data_size_gather, &
                                   product_matrix_scatter%m%data_area, &
                                   data_size_gather, data_displ_gather, &
                                   right_myprow, grp_gather)
         ENDIF
         IF (do_diagonal) THEN
            ! Merge gathered images
            DO iproc = 0, right_npcols-1
               IF (meta_size_gather(iproc) .GT. 1) THEN
                  product_matrix_scatter%m%index(dbcsr_slot_nblks) = &
                     product_matrix_scatter%m%index(dbcsr_slot_nblks)+ &
                     local_meta_scatter(meta_displ_gather(iproc)+dbcsr_slot_nblks)
                  product_matrix_scatter%m%index(dbcsr_slot_nze) = &
                     product_matrix_scatter%m%index(dbcsr_slot_nze)+ &
                     local_meta_scatter(meta_displ_gather(iproc)+dbcsr_slot_nze)
                  product_matrix_scatter%m%row_p(:) = &
                     product_matrix_scatter%m%row_p(:)+ &
                     local_meta_scatter(meta_displ_gather(iproc)+ &
                                        local_meta_scatter(meta_displ_gather(iproc)+ &
                                                           dbcsr_slot_row_p): &
                                        meta_displ_gather(iproc)+ &
                                        local_meta_scatter(meta_displ_gather(iproc)+ &
                                                           dbcsr_slot_row_p+1))
               ENDIF
            ENDDO
            product_matrix_scatter%m%index(dbcsr_slot_size) = &
               product_matrix_scatter%m%index(dbcsr_slot_size)+ &
               product_matrix_scatter%m%index(dbcsr_slot_nblks)*2
            product_matrix_scatter%m%index(dbcsr_slot_col_i+1) = &
               product_matrix_scatter%m%index(dbcsr_slot_col_i)+ &
               product_matrix_scatter%m%index(dbcsr_slot_nblks)-1
            product_matrix_scatter%m%index(dbcsr_slot_blk_p) = &
               product_matrix_scatter%m%index(dbcsr_slot_col_i+1)+1
            product_matrix_scatter%m%index(dbcsr_slot_blk_p+1) = &
               product_matrix_scatter%m%index(dbcsr_slot_blk_p)+ &
               product_matrix_scatter%m%index(dbcsr_slot_nblks)-1
            CALL ensure_array_size(product_matrix_scatter%m%index, &
                                   ub=product_matrix_scatter%m%index(dbcsr_slot_size), &
                                   nocopy=.FALSE., memory_type=memtype_mpi_buffer)
            !
            CALL dbcsr_repoint_index(product_matrix_scatter%m)
            !
            ! Collect col_i and blk_p, row by row
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (product_matrix_scatter,meta_size_gather,&
!$OMP         local_meta_scatter,right_npcols,meta_displ_gather,&
!$OMP         data_displ_gather) &
!$OMP PRIVATE (local_meta_size_gather,iproc,col,v_ki,&
!$OMP          col_i,blk_p,row)
!$OMP DO
            DO row = 0, SIZE(product_matrix_scatter%m%row_p)-2
               local_meta_size_gather = product_matrix_scatter%m%row_p(row+1)+1
               ! Skip empty rows
               IF (local_meta_size_gather .GT. product_matrix_scatter%m%row_p(row+2)) CYCLE
               DO iproc = 0, right_npcols-1
                  IF (meta_size_gather(iproc) .GT. 1) THEN
                     col = meta_displ_gather(iproc)+ &
                           local_meta_scatter(meta_displ_gather(iproc)+dbcsr_slot_row_p)+row
                     v_ki = local_meta_scatter(col+1)-local_meta_scatter(col)
                     ! Skip empty rows
                     IF (v_ki .LE. 0) CYCLE
                     col_i => local_meta_scatter(meta_displ_gather(iproc)+ &
                                                 local_meta_scatter(meta_displ_gather(iproc)+dbcsr_slot_col_i): &
                                                 meta_displ_gather(iproc)+ &
                                                 local_meta_scatter(meta_displ_gather(iproc)+dbcsr_slot_col_i+1))
                     blk_p => local_meta_scatter(meta_displ_gather(iproc)+ &
                                                 local_meta_scatter(meta_displ_gather(iproc)+dbcsr_slot_blk_p): &
                                                 meta_displ_gather(iproc)+ &
                                                 local_meta_scatter(meta_displ_gather(iproc)+dbcsr_slot_blk_p+1))
                     product_matrix_scatter%m%col_i(local_meta_size_gather:local_meta_size_gather+v_ki-1) = &
                        col_i(local_meta_scatter(col)+1:local_meta_scatter(col+1))
                     product_matrix_scatter%m%blk_p(local_meta_size_gather:local_meta_size_gather+v_ki-1) = &
                        blk_p(local_meta_scatter(col)+1:local_meta_scatter(col+1))+ &
                        data_displ_gather(iproc)
                     local_meta_size_gather = local_meta_size_gather+v_ki
                  ENDIF
               ENDDO
            ENDDO
!$OMP END DO
            CALL dbcsr_finalize(product_matrix_scatter)
!$OMP END PARALLEL
            CALL dbcsr_add(product_matrix, product_matrix_scatter, flop=flop)
            DEALLOCATE (meta_size_gather, data_size_gather, &
                        meta_displ_gather, data_displ_gather)
         ENDIF
         CALL dbcsr_release_locals(product_matrix_scatter)
         IF (keep_sparsity .AND. local_meta_size_scatter .GT. 1) NULLIFY (product_matrix_scatter%m%index)
         CALL dbcsr_release(product_matrix_scatter)
         CALL timestop(handle2)
      ENDIF
      !
      CALL timestop(handle)
   END SUBROUTINE cannon_multiply_low_rma

! **************************************************************************************************
!> \brief ...
!> \param v_ki ...
!> \param vpunmerged ...
!> \param vpmerged ...
!> \param unmerged_nimages ...
!> \param merged_nimages ...
!> \param vunmerged ...
!> \param vmerged ...
!> \param displ_vunmerged ...
!> \param v_i ...
!> \param nofiltered ...
!> \param set_displ ...
!> \param recv_vunmerged_displ ...
!> \param refs_displ_unmerged ...
!> \param buffer ...
!> \param refs_meta_size ...
!> \param refs_data_size ...
!> \param meta_size ...
!> \param data_size ...
!> \param min_nimages ...
!> \param meta_offset ...
!> \param data_offset ...
!> \param meta_win ...
!> \param data_win ...
!> \param mypunmerged ...
!> \param save_statistics ...
!> \param data_get ...
!> \param data_type_byte ...
!> \param buffer_orig ...
!> \param meta_offset_diag ...
!> \param data_offset_diag ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE cannon_transfer_rma(v_ki, vpunmerged, vpmerged, &
                                  unmerged_nimages, merged_nimages, &
                                  vunmerged, vmerged, &
                                  displ_vunmerged, &
                                  v_i, nofiltered, set_displ, &
                                  recv_vunmerged_displ, refs_displ_unmerged, &
                                  buffer, refs_meta_size, refs_data_size, &
                                  meta_size, data_size, &
                                  min_nimages, meta_offset, data_offset, &
                                  meta_win, data_win, mypunmerged, save_statistics, &
                                  data_get, data_type_byte, &
                                  buffer_orig, meta_offset_diag, data_offset_diag)
      INTEGER, INTENT(IN)                                :: v_ki, vpunmerged, vpmerged, &
                                                            unmerged_nimages, merged_nimages
      INTEGER, DIMENSION(0:), INTENT(IN)                 :: vunmerged
      INTEGER, DIMENSION(:, 0:), INTENT(IN)              :: vmerged
      INTEGER, DIMENSION(0:), INTENT(IN)                 :: displ_vunmerged
      INTEGER, INTENT(INOUT)                             :: v_i
      LOGICAL, DIMENSION(:, 0:), INTENT(IN)              :: nofiltered
      LOGICAL, INTENT(INOUT)                             :: set_displ
      INTEGER, INTENT(INOUT)                             :: recv_vunmerged_displ
      INTEGER, ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: refs_displ_unmerged
      TYPE(dbcsr_rma_buffer), INTENT(INOUT)              :: buffer
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(IN)  :: refs_meta_size
      INTEGER, DIMENSION(:), INTENT(IN)                  :: refs_data_size
      INTEGER, INTENT(INOUT)                             :: meta_size, data_size
      INTEGER, INTENT(IN)                                :: min_nimages
      INTEGER, INTENT(INOUT)                             :: meta_offset, data_offset
      INTEGER, INTENT(IN)                                :: meta_win, data_win, mypunmerged
      LOGICAL, INTENT(IN)                                :: save_statistics
      TYPE(dbcsr_data_obj)                               :: data_get
      INTEGER, INTENT(IN)                                :: data_type_byte
      TYPE(dbcsr_rma_buffer), INTENT(IN)                 :: buffer_orig
      INTEGER, INTENT(INOUT)                             :: meta_offset_diag, data_offset_diag

      INTEGER                                            :: recv_punmerged, recv_vunmerged, ui, &
                                                            ui_next, v_ki_next
      INTEGER, DIMENSION(:), POINTER                     :: meta_get
      INTEGER, DIMENSION(idata_displ:imeta_displ)        :: refs_displ_unmerged_first
      LOGICAL                                            :: is_contiguous

! Left  : merged => row, unmerged => col
! Right : merged => col, unmerged => row

      recv_vunmerged = v_ki+vpunmerged
      recv_punmerged = recv_vunmerged/unmerged_nimages
      refs_displ_unmerged_first(:) = 0
      DO ui = 1, vunmerged(v_ki)
         v_i = v_i+1
         IF (nofiltered(ui, v_ki)) THEN
            buffer%offset(idata_size, v_i+1) = &
               buffer%offset(idata_size, v_i)+ &
               refs_data_size(displ_vunmerged(recv_vunmerged)+ui)
            buffer%offset(imeta_size, v_i+1) = &
               buffer%offset(imeta_size, v_i)+ &
               refs_meta_size(vmerged(ui, v_ki), recv_vunmerged)
            meta_size = meta_size+refs_meta_size(vmerged(ui, v_ki), recv_vunmerged)
            data_size = data_size+refs_data_size(displ_vunmerged(recv_vunmerged)+ui)
            !
            ! Check for diagonal images and do the local copy
            IF ((merged_nimages .GT. 1 .AND. vmerged(ui, v_ki) .EQ. recv_vunmerged) .OR. &
                (merged_nimages .EQ. 1 .AND. recv_punmerged .EQ. vpmerged)) THEN
               buffer%nrequests = buffer%nrequests+1
               buffer%get_requests_map(1, buffer%nrequests) = v_ki+1
               buffer%get_requests_map(2, buffer%nrequests) = v_ki+1
               CALL memory_copy(buffer%meta(meta_offset+1:meta_size), &
                                buffer_orig%meta_diag(meta_offset_diag+1: &
                                                      meta_offset_diag+refs_meta_size(vmerged(ui, v_ki), &
                                                                                      recv_vunmerged)), &
                                refs_meta_size(vmerged(ui, v_ki), recv_vunmerged))
               meta_offset = meta_size
               CALL dbcsr_data_set(buffer%data, data_offset+1, &
                                   refs_data_size(displ_vunmerged(recv_vunmerged)+ui), &
                                   buffer_orig%data_diag, data_offset_diag+1)
               data_offset = data_size
            ELSE
               IF (set_displ) THEN
                  recv_vunmerged_displ = recv_vunmerged
                  refs_displ_unmerged(:, recv_vunmerged_displ) = &
                     refs_displ_unmerged(:, recv_vunmerged_displ)+ &
                     refs_displ_unmerged_first(:)
                  refs_displ_unmerged_first(:) = 0
                  set_displ = .FALSE.
                  buffer%get_requests_map(1, buffer%nrequests+1) = v_ki+1
               ENDIF
               !
               ! Check for contiguous data
               is_contiguous = .FALSE.
               ui_next = ui+1
               v_ki_next = v_ki
               DO WHILE (.TRUE.)
                  IF (ui_next .LE. vunmerged(v_ki_next)) THEN
                     IF (nofiltered(ui_next, v_ki_next) .AND. &
                         vmerged(ui_next, v_ki_next) .NE. v_ki_next+vpunmerged) THEN
                        is_contiguous = .TRUE.
                     ENDIF
                     EXIT
                  ELSEIF (v_ki_next .LT. min_nimages-1) THEN
                     v_ki_next = v_ki_next+1
                     IF (recv_punmerged .NE. ((v_ki_next+vpunmerged)/unmerged_nimages)) EXIT
                     ui_next = 1
                  ELSE
                     EXIT
                  ENDIF
               ENDDO
               !
               IF (.NOT. is_contiguous) THEN
                  buffer%nrequests = buffer%nrequests+1
                  buffer%get_requests_map(2, buffer%nrequests) = v_ki+1
                  meta_get => buffer%meta(meta_offset+1:meta_size)
                  CALL mp_rget(meta_get, recv_punmerged, &
                               meta_win, &
                               buffer_orig%meta, &
                               mypunmerged, &
                               refs_displ_unmerged(imeta_displ, recv_vunmerged_displ), &
                               buffer%get_requests(1, buffer%nrequests))
                  meta_offset = meta_size
                  CALL dbcsr_data_set_pointer( &
                     area=data_get, &
                     rsize=data_size-data_offset, &
                     csize=1, &
                     pointee=buffer%data, &
                     source_lb=data_offset+1)
                  CALL dbcsr_rget_any(data_get, recv_punmerged, &
                                      data_win, &
                                      buffer_orig%data, &
                                      mypunmerged, &
                                      refs_displ_unmerged(idata_displ, recv_vunmerged_displ), &
                                      buffer%get_requests(2, buffer%nrequests))
                  IF (save_statistics) THEN
                     CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(1, :), &
                                               data_size-data_offset, &
                                               dbcsr_mpi_statistics%data_size_breakdown(:, :, 1), data_type_byte)
                     dbcsr_mpi_statistics%nexchanged = dbcsr_mpi_statistics%nexchanged+1
                  ENDIF
                  data_offset = data_size
                  set_displ = .TRUE.
               ENDIF
            ENDIF
         ELSE
            buffer%offset(:, v_i+1) = buffer%offset(:, v_i)
            dbcsr_mpi_statistics%nfiltered = dbcsr_mpi_statistics%nfiltered+1
         ENDIF
         IF ((merged_nimages .GT. 1 .AND. vmerged(ui, v_ki) .EQ. recv_vunmerged) .OR. &
             (merged_nimages .EQ. 1 .AND. recv_punmerged .EQ. vpmerged)) THEN
            meta_offset_diag = meta_offset_diag+ &
                               refs_meta_size(vmerged(ui, v_ki), recv_vunmerged)
            data_offset_diag = data_offset_diag+ &
                               refs_data_size(displ_vunmerged(recv_vunmerged)+ui)
         ELSE
            refs_displ_unmerged_first(imeta_displ) = &
               refs_displ_unmerged_first(imeta_displ)+ &
               refs_meta_size(vmerged(ui, v_ki), recv_vunmerged)
            refs_displ_unmerged_first(idata_displ) = &
               refs_displ_unmerged_first(idata_displ)+ &
               refs_data_size(displ_vunmerged(recv_vunmerged)+ui)
         ENDIF
      ENDDO
   END SUBROUTINE cannon_transfer_rma

! **************************************************************************************************
!> \brief ...
!> \param right_buffer ...
!> \param left_buffer ...
!> \param in_nrequests ...
!> \param multrec ...
!> \param do_rec ...
!> \param flop ...
!> \param get_requests_images ...
!> \param imgdist ...
!> \param is_left ...
!> \param nrows_images ...
!> \param ncols_images ...
!> \param global_indices ...
!> \param vunmerged ...
!> \param vmerged ...
!> \param otf_filtering ...
!> \param data_type ...
!> \param right_norms ...
!> \param left_norms ...
!> \param common_unmerged ...
!> \param indices_common_unmerged ...
!> \param nthreads ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE wait_merge_calc_low(right_buffer, left_buffer, &
                                  in_nrequests, &
                                  multrec, do_rec, flop, &
                                  get_requests_images, &
                                  imgdist, is_left, &
                                  nrows_images, ncols_images, &
                                  global_indices, vunmerged, vmerged, &
                                  otf_filtering, data_type, &
                                  right_norms, left_norms, &
                                  common_unmerged, indices_common_unmerged, &
                                  nthreads)
      TYPE(dbcsr_rma_buffer), INTENT(INOUT), TARGET      :: right_buffer, left_buffer
      INTEGER, INTENT(INOUT)                             :: in_nrequests
      TYPE(dbcsr_mm_multrec_type_p), DIMENSION(0:), &
         INTENT(INOUT)                                   :: multrec
      LOGICAL, INTENT(IN)                                :: do_rec
      INTEGER(KIND=int_8), INTENT(INOUT)                 :: flop
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: get_requests_images
      TYPE(dbcsr_imagedistribution_obj), INTENT(INOUT)   :: imgdist
      LOGICAL, INTENT(IN)                                :: is_left
      INTEGER, INTENT(IN)                                :: nrows_images, ncols_images
      INTEGER, DIMENSION(:), INTENT(IN)                  :: global_indices, vunmerged
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: vmerged
      LOGICAL, INTENT(IN)                                :: otf_filtering
      INTEGER, INTENT(IN)                                :: data_type
      REAL(kind=sp), DIMENSION(:, :), INTENT(INOUT)      :: right_norms, left_norms
      LOGICAL, DIMENSION(:), INTENT(INOUT)               :: common_unmerged
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: indices_common_unmerged
      INTEGER, INTENT(IN), OPTIONAL                      :: nthreads

      INTEGER                                            :: iget, nrequests, uf, ul
      TYPE(dbcsr_rma_buffer), POINTER                    :: buffer

      IF (is_left) THEN
         buffer => left_buffer
      ELSE
         buffer => right_buffer
      ENDIF
      nrequests = 0
      DO iget = 1, in_nrequests
         IF (mp_testall(buffer%get_requests(:, iget))) THEN
            get_requests_images(buffer%get_requests_map(:, iget)) = &
               get_requests_images(buffer%get_requests_map(:, iget))-1
            uf = buffer%get_requests_map(1, iget)
            IF (get_requests_images(uf) .NE. 0) uf = uf+1
            ul = buffer%get_requests_map(2, iget)
            IF (get_requests_images(ul) .NE. 0) ul = ul-1
            IF (uf .LE. ul) THEN
               CALL merge_calc(right_buffer, left_buffer, &
                               multrec, do_rec, flop, &
                               imgdist, is_left, &
                               uf, ul, &
                               nrows_images, ncols_images, &
                               global_indices, vunmerged, vmerged, &
                               otf_filtering, data_type, &
                               right_norms, left_norms, &
                               common_unmerged, indices_common_unmerged, &
                               .FALSE., nthreads)
            ENDIF
         ELSE
            nrequests = nrequests+1
            IF (nrequests .NE. iget) THEN
               buffer%get_requests(:, nrequests) = &
                  buffer%get_requests(:, iget)
               buffer%get_requests_map(:, nrequests) = &
                  buffer%get_requests_map(:, iget)
            ENDIF
         ENDIF
      ENDDO
      !
      in_nrequests = nrequests
   END SUBROUTINE wait_merge_calc_low

! **************************************************************************************************
!> \brief ...
!> \param right_buffer ...
!> \param left_buffer ...
!> \param multrec ...
!> \param do_rec ...
!> \param flop ...
!> \param imgdist ...
!> \param is_left ...
!> \param uf ...
!> \param ul ...
!> \param nrows_images ...
!> \param ncols_images ...
!> \param global_indices ...
!> \param vunmerged ...
!> \param vmerged ...
!> \param otf_filtering ...
!> \param data_type ...
!> \param right_norms ...
!> \param left_norms ...
!> \param common_unmerged ...
!> \param indices_common_unmerged ...
!> \param has_acc ...
!> \param nthreads ...
!> \param trs_stackbuf ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE merge_calc(right_buffer, left_buffer, &
                         multrec, do_rec, flop, &
                         imgdist, is_left, &
                         uf, ul, &
                         nrows_images, ncols_images, &
                         global_indices, vunmerged, vmerged, &
                         otf_filtering, data_type, &
                         right_norms, left_norms, &
                         common_unmerged, indices_common_unmerged, &
                         has_acc, nthreads, trs_stackbuf)
      TYPE(dbcsr_rma_buffer), INTENT(INOUT), TARGET      :: right_buffer, left_buffer
      TYPE(dbcsr_mm_multrec_type_p), DIMENSION(0:), &
         INTENT(INOUT)                                   :: multrec
      LOGICAL, INTENT(IN)                                :: do_rec
      INTEGER(KIND=int_8), INTENT(INOUT)                 :: flop
      TYPE(dbcsr_imagedistribution_obj), INTENT(INOUT)   :: imgdist
      LOGICAL, INTENT(IN)                                :: is_left
      INTEGER, INTENT(IN)                                :: uf, ul, nrows_images, ncols_images
      INTEGER, DIMENSION(:), INTENT(IN)                  :: global_indices, vunmerged
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: vmerged
      LOGICAL, INTENT(IN)                                :: otf_filtering
      INTEGER, INTENT(IN)                                :: data_type
      REAL(kind=sp), DIMENSION(:, :), INTENT(INOUT), &
         TARGET                                          :: right_norms, left_norms
      LOGICAL, DIMENSION(:), INTENT(INOUT)               :: common_unmerged
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: indices_common_unmerged
      LOGICAL, INTENT(IN)                                :: has_acc
      INTEGER, INTENT(IN), OPTIONAL                      :: nthreads
      TYPE(dbcsr_data_obj), INTENT(INOUT), OPTIONAL      :: trs_stackbuf

      INTEGER                                            :: handle, ithread, ncommons_unmerged, ui
      REAL(kind=sp), DIMENSION(:, :), POINTER            :: norms
      TYPE(dbcsr_rma_buffer), POINTER                    :: buffer

      IF (is_left) THEN
         buffer => left_buffer
         norms => left_norms
      ELSE
         buffer => right_buffer
         norms => right_norms
      ENDIF
      !
      ! Merge row-images/col-images of right/left matrices
      ! Repoint indices of matrices
      CALL merge_images(buffer%buffer, nrows_images, ncols_images, &
                        buffer%vprow, buffer%vpcol, &
                        buffer%offset, buffer%meta, &
                        uf=uf, ul=ul, &
                        imgdist=imgdist, do_merge_rows=is_left, &
                        global_indices=global_indices, &
                        vunmerged=vunmerged, &
                        vmerged=vmerged, &
                        common_unmerged=common_unmerged, &
                        indices_common_unmerged=indices_common_unmerged, &
                        ncommons_unmerged=ncommons_unmerged, &
                        nthreads=nthreads)
      IF (otf_filtering) THEN
         CALL calculate_image_norms(buffer%buffer, &
                                    data_type, norms, uf=uf, ul=ul)
      ENDIF
      !
      IF (ncommons_unmerged .GT. 0) THEN
         !
         IF (has_acc .AND. PRESENT(trs_stackbuf)) THEN
            CALL acc_transpose_blocks_images(right_buffer%buffer, &
                                             right_buffer%data, &
                                             trs_stackbuf, &
                                             indices_common_unmerged(1:ncommons_unmerged))
         ENDIF
         !
         CALL timeset("cannon_multiply_low_rma_multrec", handle)
         !
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (right_buffer,left_buffer, &
!$OMP         left_norms,right_norms, &
!$OMP         multrec,do_rec, &
!$OMP         ncommons_unmerged,indices_common_unmerged) &
!$OMP PRIVATE (ithread,ui) &
!$OMP REDUCTION (+: flop)
         ithread = 0
!$       ithread = omp_get_thread_num()
         !
         DO ui = 1, ncommons_unmerged
            CALL dbcsr_mm_multrec_multiply(multrec(ithread)%p, &
                                           left=left_buffer%buffer%mats(indices_common_unmerged(ui))%m, &
                                           right=right_buffer%buffer%mats(indices_common_unmerged(ui))%m, &
                                           flop=flop, &
                                           a_norms=left_norms(:, indices_common_unmerged(ui)), &
                                           b_norms=right_norms(:, indices_common_unmerged(ui)), &
                                           do_rec=do_rec)
         ENDDO
!$OMP END PARALLEL
         !
         CALL timestop(handle)
      ENDIF
   END SUBROUTINE merge_calc

! **************************************************************************************************
!> \brief ...
!> \param mpi_statistics ...
!> \param data_size ...
!> \param size_breakdown ...
!> \param element_size_bytes ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE count_mpi_statistics(mpi_statistics, data_size, &
                                   size_breakdown, element_size_bytes)
      REAL, DIMENSION(:), INTENT(INOUT)                  :: mpi_statistics
      INTEGER, INTENT(IN)                                :: data_size
      INTEGER(KIND=int_8), DIMENSION(:, :), &
         INTENT(INOUT), OPTIONAL                         :: size_breakdown
      INTEGER, INTENT(IN), OPTIONAL                      :: element_size_bytes

      INTEGER                                            :: ilimit, nlimits
      INTEGER(KIND=int_8)                                :: data_size_bytes, llimit

      mpi_statistics(1) = mpi_statistics(1)+REAL(data_size)
      mpi_statistics(2) = MIN(mpi_statistics(2), REAL(data_size))
      mpi_statistics(3) = MAX(mpi_statistics(3), REAL(data_size))
      IF (PRESENT(size_breakdown)) THEN
         data_size_bytes = data_size
         ! change in bytes
         IF (PRESENT(element_size_bytes)) data_size_bytes = data_size_bytes*element_size_bytes
         nlimits = SIZE(dbcsr_mpi_size_limits)
         ! check for oversize messages
         IF (data_size_bytes .GT. dbcsr_mpi_size_limits(nlimits)) THEN
            size_breakdown(nlimits+1, 1) = size_breakdown(nlimits+1, 1)+1
            size_breakdown(nlimits+1, 2) = size_breakdown(nlimits+1, 2)+data_size_bytes
            RETURN
         ENDIF
         llimit = 0
         DO ilimit = 1, nlimits
            IF (data_size_bytes .GE. llimit .AND. data_size_bytes .LE. dbcsr_mpi_size_limits(ilimit)) THEN
               size_breakdown(ilimit, 1) = size_breakdown(ilimit, 1)+1
               size_breakdown(ilimit, 2) = size_breakdown(ilimit, 2)+data_size_bytes
               RETURN
            ENDIF
            llimit = dbcsr_mpi_size_limits(ilimit)
         ENDDO
      ENDIF
   END SUBROUTINE count_mpi_statistics

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param source_matrix ...
!> \param index_size ...
!> \param data_size ...
!> \param data_buffer ...
!> \param data_memory_type ...
! **************************************************************************************************
   SUBROUTINE setup_buffer_matrix(matrix, source_matrix, &
                                  index_size, data_size, data_buffer, data_memory_type)
      TYPE(dbcsr_obj), INTENT(INOUT)                     :: matrix
      TYPE(dbcsr_obj), INTENT(IN)                        :: source_matrix
      INTEGER, INTENT(IN)                                :: index_size
      INTEGER, INTENT(IN), OPTIONAL                      :: data_size
      TYPE(dbcsr_data_obj), INTENT(IN), OPTIONAL         :: data_buffer
      TYPE(dbcsr_memtype_type), INTENT(IN), OPTIONAL     :: data_memory_type

      CALL dbcsr_init(matrix)
      CALL dbcsr_create(matrix, &
                        template=source_matrix, &
                        name=TRIM("Buffer of "//TRIM(source_matrix%m%name)), &
                        nze=data_size, &
                        data_buffer=data_buffer, &
                        data_memory_type=data_memory_type, &
                        index_memory_type=memtype_mpi_buffer)
      IF (PRESENT(data_size)) THEN
         CALL dbcsr_data_ensure_size( &
            matrix%m%data_area, &
            data_size, nocopy=.TRUE.)
      ENDIF
      CALL ensure_array_size( &
         matrix%m%index, &
         ub=index_size, nocopy=.TRUE., &
         memory_type=dbcsr_get_index_memory_type(matrix))
      matrix%m%negate_real = source_matrix%m%negate_real
      matrix%m%negate_imaginary = source_matrix%m%negate_imaginary
      matrix%m%local_indexing = source_matrix%m%local_indexing
      matrix%m%list_indexing = source_matrix%m%list_indexing
      !
      IF (source_matrix%m%has_local_rows) THEN
         matrix%m%local_rows = source_matrix%m%local_rows
         CALL array_hold(matrix%m%local_rows)
         matrix%m%has_local_rows = .TRUE.
      ENDIF
      IF (source_matrix%m%has_global_rows) THEN
         matrix%m%global_rows = source_matrix%m%global_rows
         CALL array_hold(matrix%m%global_rows)
         matrix%m%has_global_rows = .TRUE.
      ENDIF
      IF (source_matrix%m%has_local_cols) THEN
         matrix%m%local_cols = source_matrix%m%local_cols
         CALL array_hold(matrix%m%local_cols)
         matrix%m%has_local_cols = .TRUE.
      ENDIF
      IF (source_matrix%m%has_global_cols) THEN
         matrix%m%global_cols = source_matrix%m%global_cols
         CALL array_hold(matrix%m%global_cols)
         matrix%m%has_global_cols = .TRUE.
      ENDIF

   END SUBROUTINE setup_buffer_matrix

! **************************************************************************************************
!> \brief ...
!> \param buffer_set ...
!> \param nimages ...
!> \param imgdist ...
!> \param template_matrix ...
!> \param index_size_merged ...
!> \param data_buffer ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE setup_buffer_matrices_images(buffer_set, nimages, imgdist, &
                                           template_matrix, index_size_merged, data_buffer)
      TYPE(dbcsr_1d_array_type), INTENT(INOUT)           :: buffer_set
      INTEGER, INTENT(IN)                                :: nimages
      TYPE(dbcsr_imagedistribution_obj), INTENT(INOUT)   :: imgdist
      TYPE(dbcsr_obj), INTENT(IN)                        :: template_matrix
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: index_size_merged
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: data_buffer

      INTEGER                                            :: v_ki

!   ---------------------------------------------------------------------------

      buffer_set%image_dist = imgdist
      CALL dbcsr_image_dist_hold(imgdist)
      !
      ALLOCATE (buffer_set%mats(nimages))
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(nimages,buffer_set,template_matrix,&
!$OMP        imgdist,data_buffer,memtype_mpi_buffer,&
!$OMP        index_size_merged)
      DO v_ki = 1, nimages
         CALL dbcsr_init(buffer_set%mats(v_ki))
         CALL dbcsr_create(buffer_set%mats(v_ki), &
                           "Buffer image of "//template_matrix%m%name, &
                           imgdist%i%main, &
                           dbcsr_type_no_symmetry, &
                           array_data(template_matrix%m%row_blk_size), array_data(template_matrix%m%col_blk_size), &
                           template_matrix%m%row_blk_size, template_matrix%m%col_blk_size, &
                           data_buffer=data_buffer, &
                           max_rbs=template_matrix%m%max_rbs, max_cbs=template_matrix%m%max_cbs, &
                           row_blk_offset=template_matrix%m%row_blk_offset, col_blk_offset=template_matrix%m%col_blk_offset, &
                           index_memory_type=memtype_mpi_buffer, &
                           make_index=.FALSE.)
         IF (PRESENT(index_size_merged)) THEN
            CALL ensure_array_size(buffer_set%mats(v_ki)%m%index, &
                                   ub=index_size_merged(v_ki), nocopy=.TRUE., &
                                   memory_type=dbcsr_get_index_memory_type(buffer_set%mats(v_ki)))
         ENDIF
         buffer_set%mats(v_ki)%m%negate_real = template_matrix%m%negate_real
         buffer_set%mats(v_ki)%m%negate_imaginary = template_matrix%m%negate_imaginary
         buffer_set%mats(v_ki)%m%local_indexing = .TRUE.
         buffer_set%mats(v_ki)%m%list_indexing = .TRUE.
      ENDDO
!$OMP END PARALLEL DO
   END SUBROUTINE setup_buffer_matrices_images

! **************************************************************************************************
!> \brief ...
!> \param buffer_set ...
!> \param buff_nrows ...
!> \param buff_ncols ...
!> \param source_matrix ...
!> \param index_size ...
!> \param data_size ...
! **************************************************************************************************
   SUBROUTINE setup_buffer_matrices(buffer_set, buff_nrows, buff_ncols, &
                                    source_matrix, index_size, data_size)
      TYPE(dbcsr_2d_array_type), INTENT(OUT)             :: buffer_set
      INTEGER, INTENT(IN)                                :: buff_nrows, buff_ncols
      TYPE(dbcsr_obj), INTENT(IN)                        :: source_matrix
      INTEGER, INTENT(IN)                                :: index_size, data_size

      CHARACTER(len=*), PARAMETER :: routineN = 'setup_buffer_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: col_image, handle, row_image

!   ---------------------------------------------------------------------------

      CALL timeset(routineN, handle)

      CALL dbcsr_image_dist_init(buffer_set%image_dist)
      ALLOCATE (buffer_set%mats(buff_nrows, buff_ncols))
      DO row_image = 1, buff_nrows
         DO col_image = 1, buff_ncols
            CALL setup_buffer_matrix(buffer_set%mats(row_image, col_image), &
                                     source_matrix, index_size, data_size=data_size, &
                                     data_memory_type=memtype_abpanel_2)
         ENDDO
      ENDDO
      IF (source_matrix%m%local_indexing .AND. careful_mod) THEN
         IF (.NOT. array_exists(source_matrix%m%local_rows)) &
            CPABORT("Local rows should exist.")
         IF (.NOT. array_exists(source_matrix%m%global_rows)) &
            CPABORT("Global rows should exist.")
         !
         IF (.NOT. array_exists(source_matrix%m%local_cols)) &
            CPABORT("Local cols should exist.")
         IF (.NOT. array_exists(source_matrix%m%global_cols)) &
            CPABORT("Global cols should exist.")
      ENDIF
      CALL timestop(handle)
   END SUBROUTINE setup_buffer_matrices

! ***************************************************************************************************
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
! **************************************************************************************************
   SUBROUTINE buffer_matrices_ensure_size(buffer_set, index_size, data_size)
      TYPE(dbcsr_2d_array_type), INTENT(INOUT)           :: buffer_set
      INTEGER, INTENT(IN)                                :: index_size, data_size

      CHARACTER(len=*), PARAMETER :: routineN = 'buffer_matrices_ensure_size', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: col_image, handle, row_image

!   ---------------------------------------------------------------------------

      CALL timeset(routineN, handle)

      DO row_image = 1, SIZE(buffer_set%mats, 1)
         DO col_image = 1, SIZE(buffer_set%mats, 2)
            CALL dbcsr_data_ensure_size(buffer_set%mats(row_image, col_image)%m%data_area, &
                                        data_size)
            CALL ensure_array_size(buffer_set%mats(row_image, col_image)%m%index, &
                                   ub=index_size, &
                                   memory_type=dbcsr_get_index_memory_type(buffer_set%mats(row_image, col_image)))
            CALL dbcsr_repoint_index(buffer_set%mats(row_image, col_image)%m)
         ENDDO
      ENDDO
      CALL timestop(handle)
   END SUBROUTINE buffer_matrices_ensure_size

! **************************************************************************************************
!> \brief ...
!> \param meta_buffer ...
!> \param img_refs_rows ...
!> \param img_refs_cols ...
!> \param refs_size ...
!> \param refs_displ ...
!> \param size_index_unmerged ...
!> \param has_threads ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE setup_rec_index_images(meta_buffer, img_refs_rows, img_refs_cols, &
                                     refs_size, refs_displ, size_index_unmerged, has_threads)
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: meta_buffer
      INTEGER, DIMENSION(:), INTENT(IN)                  :: img_refs_rows, img_refs_cols, refs_size, &
                                                            refs_displ
      INTEGER, INTENT(IN)                                :: size_index_unmerged
      LOGICAL, INTENT(IN)                                :: has_threads

      CHARACTER(len=*), PARAMETER :: routineN = 'setup_rec_index_images', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, in, nblkcols_local, &
                                                            nblkrows_local, t_f, t_l, t_size

!$    INTEGER                           :: ithread

!   ---------------------------------------------------------------------------

      CALL timeset(routineN, handle)
      IF (has_threads) THEN
         nblkrows_local = img_refs_rows(1)
      ELSE
         nblkcols_local = img_refs_cols(1)
      ENDIF
      !
      DO in = 1, SIZE(refs_size)
         IF (refs_size(in) .EQ. 0) CYCLE
         t_size = meta_buffer(refs_displ(in)+dbcsr_slot_nblks)
         IF (careful_mod) THEN
            IF (refs_size(in)-size_index_unmerged .NE. t_size*3) &
               CPABORT("Block count mismatch.")
         ENDIF
         IF (has_threads) THEN
            nblkcols_local = img_refs_cols(in)
         ELSE
            nblkrows_local = img_refs_rows(in)
         ENDIF
         t_f = 1
         t_l = t_size
!$OMP        PARALLEL IF (has_threads) DEFAULT (NONE) &
!$OMP        PRIVATE (ithread) &
!$OMP        FIRSTPRIVATE (t_f, t_l, t_size) &
!$OMP        SHARED (meta_buffer, in, has_threads, refs_displ, &
!$OMP                size_index_unmerged, nblkrows_local, nblkcols_local)
!$       ithread = OMP_GET_THREAD_NUM()+dbcsr_slot_nblks+1
!$       IF (has_threads) THEN
!$          t_f = meta_buffer(refs_displ(in)+ithread)+1
!$          t_l = meta_buffer(refs_displ(in)+ithread+1)
!$       ENDIF
         t_size = t_l-t_f+1
!$OMP        BARRIER
         IF (t_size .GT. 0) THEN
            CALL rec_sort_index(1, nblkrows_local, &
                                1, nblkcols_local, &
                                t_size, &
                                meta_buffer(refs_displ(in)+size_index_unmerged+t_f*3-2: &
                                            refs_displ(in)+size_index_unmerged+t_l*3), &
                                0)
         ENDIF
!$OMP        END PARALLEL
      ENDDO
      CALL timestop(handle)
   END SUBROUTINE setup_rec_index_images

! **************************************************************************************************
!> \brief ...
!> \param matrix_set ...
!> \param n_rows ...
!> \param n_cols ...
! **************************************************************************************************
   SUBROUTINE setup_rec_index_2d(matrix_set, n_rows, n_cols)
      TYPE(dbcsr_2d_array_type), INTENT(INOUT)           :: matrix_set
      INTEGER, INTENT(IN)                                :: n_rows, n_cols

      CHARACTER(len=*), PARAMETER :: routineN = 'setup_rec_index_2d', &
         routineP = moduleN//':'//routineN
      LOGICAL, PARAMETER                                 :: dbg = .FALSE.

      INTEGER                                            :: handle, i_col, i_row, t_f, t_l, t_size

!$    INTEGER                                  :: ithread
      LOGICAL                                  :: thread_redist

!   ---------------------------------------------------------------------------

      CALL timeset(routineN, handle)
      DO i_row = 1, n_rows
         DO i_col = 1, n_cols
            IF (.FALSE.) &
               CALL dbcsr_reset_vlocals(matrix_set%mats(i_row, i_col), &
                                        matrix_set%image_dist)
            IF (dbg) THEN
               WRITE (*, *) routineN//" m, n, size", &
                  SIZE(matrix_set%mats(i_row, i_col)%m%coo_l), &
                  dbcsr_nblkrows_local(matrix_set%mats(i_row, i_col)), &
                  dbcsr_nblkcols_local(matrix_set%mats(i_row, i_col))
               WRITE (*, '(3(1X,I7))') matrix_set%mats(i_row, i_col)%m%coo_l
            ENDIF
            IF (careful_mod) THEN
               IF (SIZE(matrix_set%mats(i_row, i_col)%m%coo_l) .NE. matrix_set%mats(i_row, i_col)%m%nblks*3) &
                  CPABORT("Block count mismatch.")
            ENDIF
            thread_redist = ASSOCIATED(matrix_set%mats(i_row, i_col)%m%thr_c)
            t_size = SIZE(matrix_set%mats(i_row, i_col)%m%coo_l)/3
            t_f = 1
            t_l = t_size
!$OMP           PARALLEL IF (thread_redist) DEFAULT (NONE) &
!$OMP           PRIVATE (ithread) &
!$OMP           FIRSTPRIVATE (t_f, t_l, t_size) &
!$OMP           SHARED (matrix_set, i_row, i_col, thread_redist)
!$          ithread = OMP_GET_THREAD_NUM()
!$          IF (thread_redist) THEN
!$             t_f = matrix_set%mats(i_row, i_col)%m%thr_c(ithread+1)+1
!$             t_l = matrix_set%mats(i_row, i_col)%m%thr_c(ithread+2)
!$          ENDIF
            t_size = t_l-t_f+1
!$OMP           BARRIER
            IF (t_size .GT. 0) THEN
               IF (matrix_set%mats(i_row, i_col)%m%local_indexing) THEN
                  CALL call_rec_sort_index( &
                     dbcsr_nblkrows_local(matrix_set%mats(i_row, i_col)), &
                     dbcsr_nblkcols_local(matrix_set%mats(i_row, i_col)), &
                     t_size, &
                     matrix_set%mats(i_row, i_col)%m%coo_l((t_f*3-2):(t_l*3)))
               ELSE
                  CALL call_rec_sort_index( &
                     dbcsr_nblkrows_total(matrix_set%mats(i_row, i_col)), &
                     dbcsr_nblkcols_total(matrix_set%mats(i_row, i_col)), &
                     t_size, &
                     matrix_set%mats(i_row, i_col)%m%coo_l((t_f*3-2):(t_l*3)))
               ENDIF
            ENDIF
!$OMP           END PARALLEL
         ENDDO
      ENDDO
      CALL timestop(handle)
   END SUBROUTINE setup_rec_index_2d

! **************************************************************************************************
!> \brief Used to thunk a call to rec_sort_index
!> \param m ...
!> \param n ...
!> \param nblks ...
!> \param idx ...
! **************************************************************************************************
   SUBROUTINE call_rec_sort_index(m, n, nblks, idx)
      INTEGER, INTENT(IN)                                :: m, n, nblks
      INTEGER, DIMENSION(3, 1:nblks), INTENT(INOUT)      :: idx

      CHARACTER(len=*), PARAMETER :: routineN = 'call_rec_sort_index', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

!   ---------------------------------------------------------------------------

      CALL timeset(routineN, handle)
      IF (.FALSE.) WRITE (*, *) " Calling rec_sort_index, size", nblks
      CALL rec_sort_index(1, m, 1, n, nblks, idx, 0)
      CALL timestop(handle)
   END SUBROUTINE call_rec_sort_index

! **************************************************************************************************
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
! **************************************************************************************************
   RECURSIVE SUBROUTINE rec_sort_index(mi, mf, ni, nf, nele, a, d)
      INTEGER, INTENT(IN)                                :: mi, mf, ni, nf, nele
      INTEGER, DIMENSION(3, 1:nele), INTENT(inout)       :: a
      INTEGER, INTENT(IN)                                :: d

      LOGICAL, PARAMETER                                 :: dbg = .FALSE.

      INTEGER                                            :: half, M, N, nlow
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: tmp

!   ---------------------------------------------------------------------------

      IF (dbg) THEN
         WRITE (*, *) " rs", mi, mf, "/", ni, nf, "=>", nele, d
         WRITE (*, '(3(1X,I7))') a(:, 1:nele)
      ENDIF
      IF (dbg) THEN
         IF (d .GT. 20) THEN
            WRITE (*, *) a(1, -d*1000)
         ENDIF
      ENDIF
      ALLOCATE (tmp(3, nele))
      M = mf-mi+1
      N = nf-ni+1
      IF (M > N) THEN
         half = M/2
         CALL rec_split(nele, a, tmp, 1, nlow, mi, half)
         a = tmp
         DEALLOCATE (tmp)
         IF (nlow .GT. 1) THEN
            CALL rec_sort_index(mi, mi+half-1, ni, nf, nlow, a(:, 1:nlow), d+1)
         ENDIF
         IF (nele-nlow .GT. 1) THEN
            CALL rec_sort_index(mi+half, mf, ni, nf, nele-nlow, a(:, nlow+1:nele), d+1)
         ENDIF
      ELSE
         half = N/2
         CALL rec_split(nele, a, tmp, 2, nlow, ni, half)
         a = tmp
         DEALLOCATE (tmp)
         IF (nlow .GT. 1) THEN
            CALL rec_sort_index(mi, mf, ni, ni+half-1, nlow, a(:, 1:nlow), d+1)
         ENDIF
         IF (nele-nlow .GT. 1) THEN
            CALL rec_sort_index(mi, mf, ni+half, nf, nele-nlow, a(:, nlow+1:nele), d+1)
         ENDIF
      ENDIF
   END SUBROUTINE rec_sort_index

! **************************************************************************************************
!> \brief ...
!> \param nele ...
!> \param a ...
!> \param split ...
!> \param row_or_col ...
!> \param nlow ...
!> \param mi ...
!> \param half ...
! **************************************************************************************************
   SUBROUTINE rec_split(nele, a, split, row_or_col, nlow, mi, half)
      INTEGER, INTENT(IN)                                :: nele
      INTEGER, DIMENSION(3, nele), INTENT(IN)            :: a
      INTEGER, DIMENSION(3, nele), INTENT(OUT)           :: split
      INTEGER, INTENT(IN)                                :: row_or_col
      INTEGER, INTENT(OUT)                               :: nlow
      INTEGER, INTENT(IN)                                :: mi, half

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rec_split', routineP = moduleN//':'//routineN

      INTEGER                                            :: el, half_m, p_high, p_low

      half_m = mi+half-1
      p_low = 1
      p_high = nele
      DO el = 1, nele
         IF (a(row_or_col, el) <= half_m) THEN
            split(1:3, p_low) = a(1:3, el)
            p_low = p_low+1
         ELSE
            split(1:3, p_high) = a(1:3, el)
            p_high = p_high-1
         ENDIF
      ENDDO
      nlow = p_low-1
      CPASSERT(p_high .EQ. nlow)

   END SUBROUTINE rec_split

! **************************************************************************************************
!> \brief Switches pointers between two matrix sets
!> \param[in,out] set1p ...
!> \param[in,out] set2p ...
! **************************************************************************************************
   SUBROUTINE dbcsr_switch_sets(set1p, set2p)
      TYPE(dbcsr_2d_array_type), POINTER                 :: set1p, set2p

      TYPE(dbcsr_2d_array_type), POINTER                 :: tmp_set

!   ---------------------------------------------------------------------------

      tmp_set => set1p
      set1p => set2p
      set2p => tmp_set
   END SUBROUTINE dbcsr_switch_sets

! **************************************************************************************************
!> \brief Switches pointers between two data areas
!> \param area1p ...
!> \param area2p ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbcsr_switch_d_ptrs(area1p, area2p)
      TYPE(dbcsr_data_obj), POINTER                      :: area1p, area2p

      TYPE(dbcsr_data_obj), POINTER                      :: tmp_p

      tmp_p => area1p
      area1p => area2p
      area2p => tmp_p
   END SUBROUTINE dbcsr_switch_d_ptrs

! **************************************************************************************************
! The following routines are helped here to help the compiler optimize them
! out.
! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param t ...
!> \retval blas_mat_type ...
! **************************************************************************************************
   ELEMENTAL FUNCTION blas_mat_type(t)
      LOGICAL, INTENT(IN)                                :: t
      CHARACTER                                          :: blas_mat_type

      IF (t) THEN
         blas_mat_type = 'T'
      ELSE
         blas_mat_type = 'N'
      ENDIF
   END FUNCTION blas_mat_type

! **************************************************************************************************
!> \brief ...
!> \param t ...
!> \retval flip_type ...
! **************************************************************************************************
   ELEMENTAL FUNCTION flip_type(t)
      CHARACTER, INTENT(IN)                              :: t
      CHARACTER                                          :: flip_type

      SELECT CASE (t)
      CASE ('N')
         flip_type = 'T'
      CASE ('T')
         flip_type = 'N'
      CASE DEFAULT
         flip_type = '@'
      END SELECT
   END FUNCTION flip_type

! **************************************************************************************************
!> \brief Calculates per-block norms.
!>
!> Rewritten to be very low-level.
!> \param[in,out] matrix     DBCSR matrix for which to calculate norms
!> \param[out] norms         Block norms
! **************************************************************************************************
   SUBROUTINE calculate_norms(matrix, norms)
      TYPE(dbcsr_obj), INTENT(IN)                        :: matrix
      REAL(kind=sp), DIMENSION(:), INTENT(OUT)           :: norms

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_norms', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: data_type, handle, nblks, nrows
      INTEGER, DIMENSION(1), TARGET                      :: tmp
      INTEGER, DIMENSION(:), POINTER                     :: local_cols, local_rows

!   ---------------------------------------------------------------------------

      CALL timeset(routineN, handle)
      ! Checks for validity
      IF (.NOT. dbcsr_valid_index(matrix)) &
         CPABORT("The matrix must be valid.")
      data_type = dbcsr_get_data_type(matrix)
      IF (matrix%m%local_indexing) THEN
         IF (careful_mod .AND. .NOT. array_exists(matrix%m%local_rows)) &
            CPABORT("Global row mapping should exist")
         local_rows => array_data(matrix%m%local_rows)
         nrows = SIZE(local_rows)
         local_cols => array_data(matrix%m%local_cols)
      ELSE
         local_rows => tmp ! Have something valid to point to
         local_cols => tmp
         nrows = matrix%m%nblkrows_total
      ENDIF
      IF (matrix%m%list_indexing) THEN
         nblks = matrix%m%nblks
         SELECT CASE (data_type)
         CASE (dbcsr_type_real_4)
            CALL calc_norms_list_s(norms, nblks, &
                                   matrix%m%coo_l, &
                                   array_data(matrix%m%row_blk_size), &
                                   array_data(matrix%m%col_blk_size), &
                                   dbcsr_get_data_p_s(matrix%m%data_area), &
                                   local=matrix%m%local_indexing, &
                                   local2global_rows=local_rows, &
                                   local2global_cols=local_cols)
         CASE (dbcsr_type_real_8)
            CALL calc_norms_list_d(norms, nblks, &
                                   matrix%m%coo_l, &
                                   array_data(matrix%m%row_blk_size), &
                                   array_data(matrix%m%col_blk_size), &
                                   dbcsr_get_data_p_d(matrix%m%data_area), &
                                   local=matrix%m%local_indexing, &
                                   local2global_rows=local_rows, &
                                   local2global_cols=local_cols)
         CASE (dbcsr_type_complex_4)
            CALL calc_norms_list_c(norms, nblks, &
                                   matrix%m%coo_l, &
                                   array_data(matrix%m%row_blk_size), &
                                   array_data(matrix%m%col_blk_size), &
                                   dbcsr_get_data_p_c(matrix%m%data_area), &
                                   local=matrix%m%local_indexing, &
                                   local2global_rows=local_rows, &
                                   local2global_cols=local_cols)
         CASE (dbcsr_type_complex_8)
            CALL calc_norms_list_z(norms, nblks, &
                                   matrix%m%coo_l, &
                                   array_data(matrix%m%row_blk_size), &
                                   array_data(matrix%m%col_blk_size), &
                                   dbcsr_get_data_p_z(matrix%m%data_area), &
                                   local=matrix%m%local_indexing, &
                                   local2global_rows=local_rows, &
                                   local2global_cols=local_cols)
         CASE DEFAULT
            CPABORT("Invalid data type.")
         END SELECT
      ELSE
         SELECT CASE (data_type)
         CASE (dbcsr_type_real_4)
            CALL calc_norms_s(norms, nrows, &
                              matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p, &
                              array_data(matrix%m%row_blk_size), &
                              array_data(matrix%m%col_blk_size), &
                              dbcsr_get_data_p_s(matrix%m%data_area), &
                              local=matrix%m%local_indexing, &
                              local2global=local_rows)
         CASE (dbcsr_type_real_8)
            CALL calc_norms_d(norms, nrows, &
                              matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p, &
                              array_data(matrix%m%row_blk_size), &
                              array_data(matrix%m%col_blk_size), &
                              dbcsr_get_data_p_d(matrix%m%data_area), &
                              local=matrix%m%local_indexing, &
                              local2global=local_rows)
         CASE (dbcsr_type_complex_4)
            CALL calc_norms_c(norms, nrows, &
                              matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p, &
                              array_data(matrix%m%row_blk_size), &
                              array_data(matrix%m%col_blk_size), &
                              dbcsr_get_data_p_c(matrix%m%data_area), &
                              local=matrix%m%local_indexing, &
                              local2global=local_rows)
         CASE (dbcsr_type_complex_8)
            CALL calc_norms_z(norms, nrows, &
                              matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p, &
                              array_data(matrix%m%row_blk_size), &
                              array_data(matrix%m%col_blk_size), &
                              dbcsr_get_data_p_z(matrix%m%data_area), &
                              local=matrix%m%local_indexing, &
                              local2global=local_rows)
         CASE DEFAULT
            CPABORT("Invalid data type.")
         END SELECT
      ENDIF
      !
      CALL timestop(handle)
   END SUBROUTINE calculate_norms

! **************************************************************************************************
!> \brief Calculates max norms per each image
!> \param buffer ...
!> \param refs_size ...
!> \param refs_displ ...
!> \param img_map ...
!> \param img_offset ...
!> \param row_blk_size ...
!> \param col_blk_size ...
!> \param local_rows ...
!> \param local_cols ...
!> \param slot_coo_l ...
!> \param[out] max_norms     Max norms per image array
!> \param is_left ...
!> \param off_diagonal ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE calculate_max_image_norms(buffer, refs_size, refs_displ, &
                                        img_map, img_offset, &
                                        row_blk_size, col_blk_size, &
                                        local_rows, local_cols, &
                                        slot_coo_l, &
                                        max_norms, is_left, off_diagonal)
      TYPE(dbcsr_rma_buffer), INTENT(IN)                 :: buffer
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: refs_size
      INTEGER, DIMENSION(:, :, :), INTENT(IN), POINTER   :: refs_displ
      INTEGER, DIMENSION(:), INTENT(IN)                  :: img_map, img_offset, row_blk_size, &
                                                            col_blk_size, local_rows, local_cols
      INTEGER, INTENT(IN)                                :: slot_coo_l
      REAL(kind=sp), DIMENSION(:), INTENT(INOUT)         :: max_norms
      LOGICAL                                            :: is_left, off_diagonal

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_max_image_norms', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      !
      SELECT CASE (buffer%data%d%data_type)
      CASE (dbcsr_type_real_4)
         CALL calc_max_image_norms_s(buffer%meta, dbcsr_get_data_p_s(buffer%data), &
                                     buffer%meta_diag, dbcsr_get_data_p_s(buffer%data_diag), &
                                     refs_size, refs_displ, &
                                     img_map, &
                                     img_offset, &
                                     row_blk_size, col_blk_size, &
                                     local_rows, local_cols, &
                                     slot_coo_l, &
                                     max_norms, is_left, off_diagonal)
      CASE (dbcsr_type_real_8)
         CALL calc_max_image_norms_d(buffer%meta, dbcsr_get_data_p_d(buffer%data), &
                                     buffer%meta_diag, dbcsr_get_data_p_d(buffer%data_diag), &
                                     refs_size, refs_displ, &
                                     img_map, &
                                     img_offset, &
                                     row_blk_size, col_blk_size, &
                                     local_rows, local_cols, &
                                     slot_coo_l, &
                                     max_norms, is_left, off_diagonal)
      CASE (dbcsr_type_complex_4)
         CALL calc_max_image_norms_c(buffer%meta, dbcsr_get_data_p_c(buffer%data), &
                                     buffer%meta_diag, dbcsr_get_data_p_c(buffer%data_diag), &
                                     refs_size, refs_displ, &
                                     img_map, &
                                     img_offset, &
                                     row_blk_size, col_blk_size, &
                                     local_rows, local_cols, &
                                     slot_coo_l, &
                                     max_norms, is_left, off_diagonal)
      CASE (dbcsr_type_complex_8)
         CALL calc_max_image_norms_z(buffer%meta, dbcsr_get_data_p_z(buffer%data), &
                                     buffer%meta_diag, dbcsr_get_data_p_z(buffer%data_diag), &
                                     refs_size, refs_displ, &
                                     img_map, &
                                     img_offset, &
                                     row_blk_size, col_blk_size, &
                                     local_rows, local_cols, &
                                     slot_coo_l, &
                                     max_norms, is_left, off_diagonal)
      CASE DEFAULT
         CPABORT("Invalid data type.")
      END SELECT
      !
      CALL timestop(handle)
   END SUBROUTINE calculate_max_image_norms

! **************************************************************************************************
!> \brief Calculates norms per each image
!> \param images ...
!> \param data_type ...
!> \param norms ...
!> \param uf ...
!> \param ul ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE calculate_image_norms(images, data_type, norms, uf, ul)
      TYPE(dbcsr_1d_array_type), INTENT(IN)              :: images
      INTEGER, INTENT(IN)                                :: data_type
      REAL(kind=sp), DIMENSION(:, :), INTENT(INOUT)      :: norms
      INTEGER, INTENT(IN)                                :: uf, ul

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_image_norms', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      !
      ! Check for row/columns image sizes
      IF (SIZE(images%mats, 1) .GT. SIZE(norms, 2) .OR. SIZE(images%mats, 1) .LT. (ul-uf+1)) &
         CPABORT("Wrong number of images!")
      !
      SELECT CASE (data_type)
      CASE (dbcsr_type_real_4)
         CALL calc_image_norms_s(images, norms, uf, ul)
      CASE (dbcsr_type_real_8)
         CALL calc_image_norms_d(images, norms, uf, ul)
      CASE (dbcsr_type_complex_4)
         CALL calc_image_norms_c(images, norms, uf, ul)
      CASE (dbcsr_type_complex_8)
         CALL calc_image_norms_z(images, norms, uf, ul)
      CASE DEFAULT
         CPABORT("Invalid data type.")
      END SELECT
      !
      CALL timestop(handle)
   END SUBROUTINE calculate_image_norms

! **************************************************************************************************
!> \brief write out a stack for transposing the blocks
!> \param matrix ...
!> \param trs_stackbuf ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE acc_transpose_blocks(matrix, trs_stackbuf)
      TYPE(dbcsr_obj), INTENT(IN)                        :: matrix
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: trs_stackbuf

      CHARACTER(len=*), PARAMETER :: routineN = 'acc_transpose_blocks', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk_p, col, handle, i, m, mi, mi_max, n, &
                                                            nblks, ni, ni_max, offset, row, x
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: counters, filled, offsets
      LOGICAL                                            :: local
      INTEGER, DIMENSION(:), POINTER :: blk_index, col_blk_sizes, col_blk_sizes2enum, &
         enum2col_blk_sizes, enum2row_blk_sizes, local2global_cols, local2global_rows, &
         row_blk_sizes, row_blk_sizes2enum, trs_stack
      INTEGER, DIMENSION(3, matrix%m%nblks)              :: tmp_stack

      CALL timeset(routineN, handle)

      NULLIFY (row_blk_sizes2enum, enum2row_blk_sizes)
      NULLIFY (col_blk_sizes2enum, enum2col_blk_sizes)
      NULLIFY (local2global_rows, local2global_cols, trs_stack)

      IF (.NOT. matrix%m%list_indexing) &
         CPABORT("build_trs_stack: only list_indexing supported.")
      IF (trs_stackbuf%d%data_type /= dbcsr_type_int_4) &
         CPABORT("build_trs_stack: stac_buf has wrong datatype")
      blk_index => matrix%m%coo_l
      row_blk_sizes => array_data(matrix%m%row_blk_size)
      col_blk_sizes => array_data(matrix%m%col_blk_size)
      local = matrix%m%local_indexing
      local2global_rows => array_data(matrix%m%local_rows)
      local2global_cols => array_data(matrix%m%local_cols)
      nblks = matrix%m%nblks

      ! enumerate the blocksizes to keep the following 2D-arrays small.
      CALL enumerate_blk_sizes(row_blk_sizes, row_blk_sizes2enum, enum2row_blk_sizes)
      CALL enumerate_blk_sizes(col_blk_sizes, col_blk_sizes2enum, enum2col_blk_sizes)
      mi_max = SIZE(enum2row_blk_sizes); ni_max = SIZE(enum2col_blk_sizes)
      ALLOCATE (counters(mi_max, ni_max), offsets(mi_max, ni_max), filled(mi_max, ni_max))
      counters(:, :) = 0; offsets(:, :) = 0; filled(:, :) = 0

      ! make sure buffer from previous cannon-tick was uploaded
      CALL acc_event_synchronize(trs_stackbuf%d%acc_ready)

      CALL dbcsr_data_ensure_size(trs_stackbuf, data_size=nblks, nocopy=.TRUE.)
      trs_stack => trs_stackbuf%d%i4

      ! collect block addresses and dimensions in a temporary stack
      ! while doing so, also count number of blocks per block-dimensions
      DO i = 1, nblks
         row = blk_index(3*(i-1)+1)
         col = blk_index(3*(i-1)+2)
         blk_p = blk_index(3*(i-1)+3)
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
         tmp_stack(3, i) = blk_p-1
         counters(mi, ni) = counters(mi, ni)+1
      ENDDO

      ! calculate offsets for first element of each sub-stack
      offset = 0
      DO mi = 1, mi_max
         DO ni = 1, ni_max
            offsets(mi, ni) = offset
            offset = offset+counters(mi, ni)
         ENDDO
      ENDDO

      ! write all sub-stacks into the host-pinned buffer
      DO i = 1, nblks
         mi = tmp_stack(1, i)
         ni = tmp_stack(2, i)
         blk_p = tmp_stack(3, i)
         x = offsets(mi, ni)+filled(mi, ni)+1
         trs_stack(x) = blk_p
         filled(mi, ni) = filled(mi, ni)+1
      ENDDO

      !sanity check
      DO mi = 1, mi_max
         DO ni = 1, ni_max
            IF (filled(mi, ni) /= counters(mi, ni)) &
               CPABORT("acc_transpose_blocks: bug")
         END DO
      END DO

      !transfer all stacks
      CALL dbcsr_data_host2dev(trs_stackbuf)

      ! make sure block-buffer is uploaded befor running the kernels
      CALL acc_stream_wait_event(trs_stackbuf%d%memory_type%acc_stream, matrix%m%data_area%d%acc_ready)

      ! launch kernels
      DO mi = 1, mi_max
         DO ni = 1, ni_max
            IF (counters(mi, ni) > 0) THEN
               m = enum2row_blk_sizes(mi)
               n = enum2col_blk_sizes(ni)
               CALL dbcsr_acc_transpose( &
                  trs_stack=trs_stackbuf%d%acc_devmem, &
                  offset=offsets(mi, ni), &
                  nblks=counters(mi, ni), &
                  datatype=matrix%m%data_type, &
                  buffer=matrix%m%data_area%d%acc_devmem, &
                  m=m, n=n, &
                  stream=trs_stackbuf%d%memory_type%acc_stream)
            END IF
         ENDDO
      ENDDO

      ! make sure block-buffer are not used until transpose kernels finnished
      CALL acc_event_record(trs_stackbuf%d%acc_ready, trs_stackbuf%d%memory_type%acc_stream)
      CALL acc_stream_wait_event(matrix%m%data_area%d%memory_type%acc_stream, trs_stackbuf%d%acc_ready)
      CALL acc_event_record(matrix%m%data_area%d%acc_ready, matrix%m%data_area%d%memory_type%acc_stream)

      DEALLOCATE (row_blk_sizes2enum, enum2row_blk_sizes)
      DEALLOCATE (col_blk_sizes2enum, enum2col_blk_sizes)
      CALL timestop(handle)
   END SUBROUTINE acc_transpose_blocks

! **************************************************************************************************
!> \brief write out a stack for transposing the blocks
!> \param matrices ...
!> \param data ...
!> \param trs_stackbuf ...
!> \param indices ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE acc_transpose_blocks_images(matrices, DATA, trs_stackbuf, indices)
      TYPE(dbcsr_1d_array_type), INTENT(IN)              :: matrices
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: DATA, trs_stackbuf
      INTEGER, DIMENSION(:), INTENT(IN)                  :: indices

      CHARACTER(len=*), PARAMETER :: routineN = 'acc_transpose_blocks_images', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk_p, col, handle, i, imat, m, mi, &
                                                            mi_max, n, nblks, nblks_total, ni, &
                                                            ni_max, offset, row, x
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: counters, filled, offsets, tmp_stack
      INTEGER, DIMENSION(:), POINTER :: blk_index, col_blk_sizes, col_blk_sizes2enum, &
         enum2col_blk_sizes, enum2row_blk_sizes, local2global_cols, local2global_rows, &
         row_blk_sizes, row_blk_sizes2enum, trs_stack

      CALL timeset(routineN, handle)

      NULLIFY (row_blk_sizes2enum, enum2row_blk_sizes)
      NULLIFY (col_blk_sizes2enum, enum2col_blk_sizes)
      NULLIFY (local2global_rows, local2global_cols, trs_stack)

      IF (trs_stackbuf%d%data_type /= dbcsr_type_int_4) &
         CPABORT("build_trs_stack: stac_buf has wrong datatype")

      row_blk_sizes => array_data(matrices%mats(indices(1))%m%row_blk_size)
      col_blk_sizes => array_data(matrices%mats(indices(1))%m%col_blk_size)

      ! enumerate the blocksizes to keep the following 2D-arrays small.
      CALL enumerate_blk_sizes(row_blk_sizes, row_blk_sizes2enum, enum2row_blk_sizes)
      CALL enumerate_blk_sizes(col_blk_sizes, col_blk_sizes2enum, enum2col_blk_sizes)
      mi_max = SIZE(enum2row_blk_sizes); ni_max = SIZE(enum2col_blk_sizes)
      ALLOCATE (counters(mi_max, ni_max), offsets(mi_max, ni_max), filled(mi_max, ni_max))
      counters(:, :) = 0; offsets(:, :) = 0; filled(:, :) = 0

      ! make sure block-buffer is uploaded befor running the kernels
      CALL acc_stream_wait_event(trs_stackbuf%d%memory_type%acc_stream, data%d%acc_ready)

      nblks_total = 0
      DO imat = 1, SIZE(indices)
         nblks_total = nblks_total+matrices%mats(indices(imat))%m%nblks
      ENDDO
      ALLOCATE (tmp_stack(3, nblks_total))

      nblks_total = 0
      DO imat = 1, SIZE(indices)
         blk_index => matrices%mats(indices(imat))%m%coo_l
         local2global_rows => array_data(matrices%mats(indices(imat))%m%local_rows)
         local2global_cols => array_data(matrices%mats(indices(imat))%m%local_cols)
         nblks = matrices%mats(indices(imat))%m%nblks

         ! collect block addresses and dimensions in a temporary stack
         ! while doing so, also count number of blocks per block-dimensions
         DO i = 1, nblks
            row = blk_index(3*(i-1)+1)
            col = blk_index(3*(i-1)+2)
            blk_p = blk_index(3*(i-1)+3)
            IF (blk_p == 0) CYCLE
            row = local2global_rows(row)
            col = local2global_cols(col)
            m = row_blk_sizes(row)
            n = col_blk_sizes(col)
            mi = row_blk_sizes2enum(m)
            ni = col_blk_sizes2enum(n)
            tmp_stack(1, nblks_total+i) = mi
            tmp_stack(2, nblks_total+i) = ni
            tmp_stack(3, nblks_total+i) = blk_p-1
            counters(mi, ni) = counters(mi, ni)+1
         ENDDO
         nblks_total = nblks_total+nblks
      ENDDO

      ! calculate offsets for first element of each sub-stack
      offset = 0
      DO mi = 1, mi_max
         DO ni = 1, ni_max
            offsets(mi, ni) = offset
            offset = offset+counters(mi, ni)
         ENDDO
      ENDDO

      ! make sure buffer from previous cannon-tick was uploaded
      CALL acc_event_synchronize(trs_stackbuf%d%acc_ready)
      CALL dbcsr_data_ensure_size(trs_stackbuf, data_size=nblks_total, nocopy=.TRUE.)
      trs_stack => trs_stackbuf%d%i4

      ! write all sub-stacks into the host-pinned buffer
      DO i = 1, nblks_total
         mi = tmp_stack(1, i)
         ni = tmp_stack(2, i)
         blk_p = tmp_stack(3, i)
         x = offsets(mi, ni)+filled(mi, ni)+1
         trs_stack(x) = blk_p
         filled(mi, ni) = filled(mi, ni)+1
      ENDDO

      !sanity check
      DO mi = 1, mi_max
         DO ni = 1, ni_max
            IF (filled(mi, ni) /= counters(mi, ni)) &
               CPABORT("acc_transpose_blocks: bug")
         END DO
      END DO

      !transfer all stacks
      CALL dbcsr_data_host2dev(trs_stackbuf)

      ! launch kernels
      DO mi = 1, mi_max
         DO ni = 1, ni_max
            IF (counters(mi, ni) > 0) THEN
               m = enum2row_blk_sizes(mi)
               n = enum2col_blk_sizes(ni)
               CALL dbcsr_acc_transpose( &
                  trs_stack=trs_stackbuf%d%acc_devmem, &
                  offset=offsets(mi, ni), &
                  nblks=counters(mi, ni), &
                  datatype=data%d%data_type, &
                  buffer=data%d%acc_devmem, &
                  m=m, n=n, &
                  stream=trs_stackbuf%d%memory_type%acc_stream)
            END IF
         ENDDO
      ENDDO

      DEALLOCATE (tmp_stack)

      ! make sure block-buffer are not used until transpose kernels finished
      CALL acc_event_record(trs_stackbuf%d%acc_ready, trs_stackbuf%d%memory_type%acc_stream)
      CALL acc_stream_wait_event(data%d%memory_type%acc_stream, trs_stackbuf%d%acc_ready)
      CALL acc_event_record(data%d%acc_ready, data%d%memory_type%acc_stream)

      DEALLOCATE (row_blk_sizes2enum, enum2row_blk_sizes)
      DEALLOCATE (col_blk_sizes2enum, enum2col_blk_sizes)
      CALL timestop(handle)
   END SUBROUTINE acc_transpose_blocks_images

! **************************************************************************************************
!> \brief Enumerate all occuring blocksizes
!> \param blk_sizes ...
!> \param enum ...
!> \param rev_enum ...
!>  \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE enumerate_blk_sizes(blk_sizes, enum, rev_enum)
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes, enum, rev_enum

      INTEGER                                            :: i, n

      n = MAXVAL(blk_sizes)
      ALLOCATE (enum(0:n))
      enum(:) = 0

      DO i = 1, SIZE(blk_sizes)
         enum(blk_sizes(i)) = 1
      ENDDO

      n = SUM(enum)
      ALLOCATE (rev_enum(n))

      n = 0
      DO i = 0, SIZE(enum)-1
         IF (enum(i) > 0) THEN
            n = n+1
            enum(i) = n
            rev_enum(n) = i
         END IF
      ENDDO
   END SUBROUTINE enumerate_blk_sizes

! **************************************************************************************************
!> \brief Init RMA buffer
!>
!> \param buffer ...
!> \param data_type ...
!> \param data_size ...
!> \param meta_size ...
!> \param requests_size ...
!> \param data_memory_type ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE rma_buffer_init(buffer, data_type, &
                              data_size, meta_size, &
                              requests_size, &
                              data_memory_type)
      TYPE(dbcsr_rma_buffer), INTENT(INOUT)              :: buffer
      INTEGER, INTENT(IN)                                :: data_type, data_size, meta_size
      INTEGER, INTENT(IN), OPTIONAL                      :: requests_size
      TYPE(dbcsr_memtype_type), INTENT(IN)               :: data_memory_type

      IF (.NOT. buffer%is_valid) THEN
         ! First initialization
         CALL dbcsr_data_init(buffer%data)
         CALL dbcsr_data_new(buffer%data, data_type=data_type, &
                             data_size=data_size, memory_type=data_memory_type)
         CALL dbcsr_data_set_size_referenced(buffer%data, data_size)
         CALL dbcsr_data_init(buffer%data_diag)
         buffer%is_valid = .TRUE.
      ELSE
         CALL dbcsr_data_ensure_size(buffer%data, data_size, nocopy=.TRUE.)
      ENDIF
      !
      CALL ensure_array_size(buffer%meta, ub=meta_size, nocopy=.TRUE., &
                             memory_type=memtype_mpi_buffer)
      !
      IF (PRESENT(requests_size)) THEN
         IF (ALLOCATED(buffer%offset)) THEN
            IF (SIZE(buffer%offset, 2) .LT. requests_size*2) DEALLOCATE (buffer%offset)
         ENDIF
         IF (.NOT. ALLOCATED(buffer%offset)) ALLOCATE (buffer%offset(2, requests_size*2))
         buffer%offset(:, 1) = 0
         IF (ALLOCATED(buffer%get_requests)) THEN
            IF (SIZE(buffer%get_requests, 2) .LT. requests_size) DEALLOCATE (buffer%get_requests)
         ENDIF
         IF (.NOT. ALLOCATED(buffer%get_requests)) ALLOCATE (buffer%get_requests(2, requests_size))
         IF (ALLOCATED(buffer%get_requests_map)) THEN
            IF (SIZE(buffer%get_requests_map, 2) .LT. requests_size) DEALLOCATE (buffer%get_requests_map)
         ENDIF
         IF (.NOT. ALLOCATED(buffer%get_requests_map)) ALLOCATE (buffer%get_requests_map(2, requests_size))
      ENDIF
   END SUBROUTINE rma_buffer_init

! **************************************************************************************************
!> \brief Resize diagonal buffer
!>
!> \param buffer ...
!> \param data_type ...
!> \param data_size_diag ...
!> \param meta_size_diag ...
!> \param data_memory_type ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE resize_buffer_diag(buffer, data_type, &
                                 data_size_diag, meta_size_diag, &
                                 data_memory_type)
      TYPE(dbcsr_rma_buffer), INTENT(INOUT)              :: buffer
      INTEGER, INTENT(IN)                                :: data_type, data_size_diag, meta_size_diag
      TYPE(dbcsr_memtype_type), INTENT(IN)               :: data_memory_type

      IF (.NOT. buffer%is_valid) CPABORT("Buffer not initialized!")
      !
      IF (dbcsr_data_valid(buffer%data_diag)) THEN
         CALL dbcsr_data_ensure_size(buffer%data_diag, data_size_diag, nocopy=.TRUE.)
      ELSE
         CALL dbcsr_data_init(buffer%data_diag)
         CALL dbcsr_data_new(buffer%data_diag, data_type=data_type, &
                             data_size=data_size_diag, memory_type=data_memory_type)
         CALL dbcsr_data_set_size_referenced(buffer%data_diag, data_size_diag)
      ENDIF
      !
      CALL ensure_array_size(buffer%meta_diag, ub=meta_size_diag, nocopy=.TRUE., &
                             memory_type=memtype_mpi_buffer)
   END SUBROUTINE resize_buffer_diag

! **************************************************************************************************
!> \brief Release RMA buffers
!>
!> \param buffers ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE rma_buffers_release(buffers)
      TYPE(dbcsr_rma_buffers), INTENT(INOUT)             :: buffers

      IF (buffers%left%is_valid) THEN
         CALL dbcsr_data_release(buffers%left%data)
         IF (dbcsr_data_valid(buffers%left%data_diag)) &
            CALL dbcsr_data_release(buffers%left%data_diag)
         buffers%left%is_valid = .FALSE.
      ENDIF
      IF (ASSOCIATED(buffers%left%meta)) THEN
         DEALLOCATE (buffers%left%meta)
         NULLIFY (buffers%left%meta)
      ENDIF
      IF (ASSOCIATED(buffers%left%meta_diag)) THEN
         DEALLOCATE (buffers%left%meta_diag)
         NULLIFY (buffers%left%meta_diag)
      ENDIF
      IF (buffers%right%is_valid) THEN
         CALL dbcsr_data_release(buffers%right%data)
         IF (dbcsr_data_valid(buffers%right%data_diag)) &
            CALL dbcsr_data_release(buffers%right%data_diag)
         buffers%right%is_valid = .FALSE.
      ENDIF
      IF (ASSOCIATED(buffers%right%meta)) THEN
         DEALLOCATE (buffers%right%meta)
         NULLIFY (buffers%right%meta)
      ENDIF
      IF (ASSOCIATED(buffers%right%meta_diag)) THEN
         DEALLOCATE (buffers%right%meta_diag)
         NULLIFY (buffers%right%meta_diag)
      ENDIF
      !
      IF (ALLOCATED(buffers%left%offset)) DEALLOCATE (buffers%left%offset)
      IF (ALLOCATED(buffers%right%offset)) DEALLOCATE (buffers%right%offset)
      IF (ALLOCATED(buffers%left%get_requests)) DEALLOCATE (buffers%left%get_requests)
      IF (ALLOCATED(buffers%right%get_requests)) DEALLOCATE (buffers%right%get_requests)
      IF (ALLOCATED(buffers%left%get_requests_map)) DEALLOCATE (buffers%left%get_requests_map)
      IF (ALLOCATED(buffers%right%get_requests_map)) DEALLOCATE (buffers%right%get_requests_map)
      IF (ASSOCIATED(meta_scatter)) THEN
         DEALLOCATE (meta_scatter)
         NULLIFY (meta_scatter)
      ENDIF
      IF (ASSOCIATED(local_meta_scatter)) THEN
         DEALLOCATE (local_meta_scatter)
         NULLIFY (local_meta_scatter)
         CALL dbcsr_data_release(local_data_scatter)
      ENDIF
      IF (ASSOCIATED(local_meta_product_scatter)) THEN
         DEALLOCATE (local_meta_product_scatter)
         NULLIFY (local_meta_product_scatter)
      ENDIF
      IF (dbcsr_data_valid(local_data_product_scatter)) THEN
         CALL dbcsr_data_release(local_data_product_scatter)
      ENDIF
   END SUBROUTINE rma_buffers_release

! **************************************************************************************************
!> \brief ...
!> \param meta_index ...
!> \param global_indices ...
!> \param num_slots ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE set_empty_meta_index(meta_index, global_indices, num_slots)
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: meta_index
      INTEGER, DIMENSION(:), INTENT(IN)                  :: global_indices
      INTEGER, INTENT(IN)                                :: num_slots

!   ---------------------------------------------------------------------------

      meta_index(dbcsr_slot_size) = num_slots
      meta_index(dbcsr_slot_size+1:num_slots) = 0
      meta_index(dbcsr_slot_nblkrows_total:dbcsr_slot_nfullcols_local) = &
         global_indices(:)
      meta_index(dbcsr_slot_row_p) = 1
      meta_index(dbcsr_slot_col_i) = 1
      meta_index(dbcsr_slot_blk_p) = 1
      meta_index(dbcsr_slot_coo_l) = num_slots+1
      meta_index(dbcsr_num_slots) = num_slots
   END SUBROUTINE set_empty_meta_index

! **************************************************************************************************
!> \brief Merge images
!>
!> \param images ...
!> \param nrows_images ...
!> \param ncols_images ...
!> \param vprow ...
!> \param vpcol ...
!> \param offset_images ...
!> \param meta_buffer ...
!> \param uf ...
!> \param ul ...
!> \param imgdist ...
!> \param do_merge_rows ...
!> \param global_indices ...
!> \param vunmerged ...
!> \param vmerged ...
!> \param common_unmerged ...
!> \param indices_common_unmerged ...
!> \param ncommons_unmerged ...
!> \param nthreads ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE merge_images(images, nrows_images, ncols_images, &
                           vprow, vpcol, offset_images, &
                           meta_buffer, &
                           uf, ul, &
                           imgdist, do_merge_rows, &
                           global_indices, &
                           vunmerged, vmerged, &
                           common_unmerged, &
                           indices_common_unmerged, &
                           ncommons_unmerged, &
                           nthreads)
      TYPE(dbcsr_1d_array_type), INTENT(INOUT)           :: images
      INTEGER, INTENT(IN)                                :: nrows_images, ncols_images, vprow, vpcol
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: offset_images
      INTEGER, DIMENSION(:), INTENT(IN), TARGET          :: meta_buffer
      INTEGER, INTENT(IN)                                :: uf, ul
      TYPE(dbcsr_imagedistribution_obj), INTENT(INOUT)   :: imgdist
      LOGICAL, INTENT(IN)                                :: do_merge_rows
      INTEGER, DIMENSION(:), INTENT(IN)                  :: global_indices, vunmerged
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: vmerged
      LOGICAL, DIMENSION(:), INTENT(INOUT)               :: common_unmerged
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: indices_common_unmerged
      INTEGER, INTENT(INOUT)                             :: ncommons_unmerged
      INTEGER, INTENT(IN), OPTIONAL                      :: nthreads

      CHARACTER(len=*), PARAMETER :: routineN = 'merge_images', routineP = moduleN//':'//routineN

      INTEGER :: handle, ioffset, ithread, kf, ki, kl, mi, my_nthreads, my_vpcol, my_vprow, &
         nimages_merged, nimages_unmerged, offset_data, offset_meta, size_index_merged, &
         size_index_unmerged, ui, v_merged
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ioffset_first
      INTEGER, DIMENSION(:), POINTER                     :: unmerged_indices
      TYPE(dbcsr_obj), POINTER                           :: image_obj
      TYPE(dbcsr_type), POINTER                          :: image

      CALL timeset(routineN, handle)
      !
      IF (do_merge_rows) THEN
         ! merge over rows (left)
         nimages_unmerged = ncols_images
         nimages_merged = nrows_images
         my_vprow = vprow/nimages_merged
         my_vpcol = vpcol-1
      ELSE
         ! merge over cols (right)
         nimages_unmerged = nrows_images
         nimages_merged = ncols_images
         my_vprow = vprow-1
         my_vpcol = vpcol/nimages_merged
      ENDIF
      !
      size_index_unmerged = dbcsr_slot_nblks
      size_index_merged = dbcsr_num_slots
      my_nthreads = 1
      IF (PRESENT(nthreads)) THEN
!$       size_index_unmerged = size_index_unmerged+nthreads+1
!$       size_index_merged = size_index_merged+nthreads+1
         my_nthreads = nthreads
      ENDIF
      !
      ALLOCATE (ioffset_first(nimages_unmerged))
      ioffset_first(1) = 0
      DO ui = 1, nimages_unmerged-1
         ioffset_first(ui+1) = ioffset_first(ui)+vunmerged(ui)
      ENDDO
      !
      ncommons_unmerged = 0
      !
!$OMP PARALLEL DO DEFAULT (NONE) &
!$OMP PRIVATE(image_obj,image,ioffset,ui,mi,v_merged, &
!$OMP         offset_meta,offset_data,unmerged_indices, &
!$OMP         ithread,kf,kl,ki) &
!$OMP SHARED(nimages_unmerged,images,global_indices, &
!$OMP        size_index_merged,vunmerged,do_merge_rows, &
!$OMP        my_vprow,my_vpcol,nimages_merged,nthreads, &
!$OMP        my_nthreads,ioffset_first,vmerged,offset_images, &
!$OMP        meta_buffer,size_index_unmerged,imgdist,uf,ul, &
!$OMP        ncommons_unmerged,common_unmerged,indices_common_unmerged) &
!$OMP schedule(dynamic)
      images_unmerge: DO ui = uf, ul
         image_obj => images%mats(ui)
         image => image_obj%m
         CALL set_empty_meta_index(image%index, global_indices, size_index_merged)
         IF (vunmerged(ui) .NE. 0) THEN
            !
!$OMP CRITICAL(merge_images_common)
            IF (common_unmerged(ui)) THEN
               ncommons_unmerged = ncommons_unmerged+1
               indices_common_unmerged(ncommons_unmerged) = ui
            ELSE
               common_unmerged(ui) = .TRUE.
            ENDIF
!$OMP END CRITICAL(merge_images_common)
            !
            ! Virtual coords
            IF (do_merge_rows) THEN
               image%index(dbcsr_slot_home_vprow) = my_vprow
               image%index(dbcsr_slot_home_vpcol) = my_vpcol+ui
            ELSE
               image%index(dbcsr_slot_home_vprow) = my_vprow+ui
               image%index(dbcsr_slot_home_vpcol) = my_vpcol
            ENDIF
            ! thr_c
!$          IF (PRESENT(nthreads)) THEN
!$             image%index(dbcsr_slot_thr_c) = dbcsr_num_slots+1
!$             image%index(dbcsr_slot_thr_c+1) = size_index_merged
!$          ENDIF
            !
            ! Merge images
            merge: IF (nimages_merged .GT. 1) THEN
               !
               ! Merge corresponding thread blocks for all images
               threads: DO ithread = 1, my_nthreads
                  !
                  ! Merge images per a given thread accross all images
                  ioffset = ioffset_first(ui)
                  images_merge: DO v_merged = 1, vunmerged(ui)
                     mi = vmerged(v_merged, ui)+1
                     ioffset = ioffset+1
                     offset_meta = offset_images(imeta_size, ioffset)
                     offset_data = offset_images(idata_size, ioffset)
                     ! Check for filtered images
                     IF (offset_data .LT. offset_images(idata_size, ioffset+1)) THEN
                        unmerged_indices => meta_buffer(offset_meta+dbcsr_slot_size: &
                                                        offset_meta+meta_buffer(offset_meta+dbcsr_slot_size))
                        IF (ithread .EQ. 1) THEN
                           ! nblks, nze
                           image%index(dbcsr_slot_nblks) = &
                              image%index(dbcsr_slot_nblks)+ &
                              unmerged_indices(dbcsr_slot_nblks)
                           image%index(dbcsr_slot_nze) = &
                              image%index(dbcsr_slot_nze)+ &
                              offset_images(idata_size, ioffset+1)- &
                              offset_data
                           ! threads distribution
!$                         IF (PRESENT(nthreads)) THEN
!$                            image%index(dbcsr_num_slots+1:size_index_merged) = &
!$                               image%index(dbcsr_num_slots+1:size_index_merged)+ &
!$                               unmerged_indices(dbcsr_slot_nblks+1:size_index_unmerged)
!$                         ENDIF
                        ENDIF
                        ! Copy CSR indices (3 values per block: row-index, col-index, offset)
                        ! remap indices taking in account threads and images
                        kf = size_index_unmerged+1
                        kl = unmerged_indices(dbcsr_slot_size)
!$                      IF (PRESENT(nthreads)) THEN
!$                         kf = kf+unmerged_indices(dbcsr_slot_nblks+ithread)*3
!$                         kl = size_index_unmerged+unmerged_indices(dbcsr_slot_nblks+ithread+1)*3
!$                      ENDIF
                        DO ki = kf, kl, 3
                           image%index(image%index(dbcsr_slot_size)+1) = &
                              unmerged_indices(ki)
                           image%index(image%index(dbcsr_slot_size)+2) = &
                              unmerged_indices(ki+1)
                           ! Data offset
                           image%index(image%index(dbcsr_slot_size)+3) = &
                              unmerged_indices(ki+2)+offset_data
                           ! increase meta size
                           image%index(dbcsr_slot_size) = image%index(dbcsr_slot_size)+3
                        ENDDO
                     ENDIF
                  ENDDO images_merge
               ENDDO threads
               !
               image%index(dbcsr_num_slots) = image%index(dbcsr_slot_size)
            ELSE
               ! No merging, just copy
               offset_meta = offset_images(imeta_size, ui)
               ! Check for filtered images
               offset_data = offset_images(idata_size, ui)
               IF (offset_data .LT. offset_images(idata_size, ui+1)) THEN
                  image%index(dbcsr_slot_size) = &
                     meta_buffer(offset_meta+dbcsr_slot_size)- &
                     dbcsr_slot_nblks+dbcsr_num_slots
                  image%index(dbcsr_slot_nblks) = &
                     meta_buffer(offset_meta+dbcsr_slot_nblks)
                  image%index(dbcsr_slot_nze) = &
                     offset_images(idata_size, ui+1)-offset_data
                  image%index(dbcsr_num_slots) = image%index(dbcsr_slot_size)
                  !
                  image%index(dbcsr_num_slots+1:image%index(dbcsr_num_slots)) = &
                     meta_buffer(offset_meta+dbcsr_slot_nblks+1: &
                                 offset_meta+meta_buffer(offset_meta+dbcsr_slot_size))
                  image%index(size_index_merged+3:image%index(dbcsr_num_slots):3) = &
                     image%index(size_index_merged+3:image%index(dbcsr_num_slots):3)+ &
                     offset_data
               ENDIF
            ENDIF merge
            !
            ! Reset
!$OMP CRITICAL(merge_images_reset_vlocals)
            CALL dbcsr_reset_vlocals(image_obj, imgdist, do_rows=.NOT. do_merge_rows)
!$OMP END CRITICAL(merge_images_reset_vlocals)
         ENDIF
         !
         ! Repoint index
         CALL dbcsr_repoint_index(image)
         image%valid = .TRUE.
      ENDDO images_unmerge
!$OMP END PARALLEL DO
      !
      DEALLOCATE (ioffset_first)
      CALL timestop(handle)
   END SUBROUTINE merge_images

#include "dbcsr_mm_cannon_d.f90"
#include "dbcsr_mm_cannon_z.f90"
#include "dbcsr_mm_cannon_s.f90"
#include "dbcsr_mm_cannon_c.f90"

END MODULE dbcsr_mm_cannon

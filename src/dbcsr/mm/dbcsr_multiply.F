!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Entry point of the dbcsr matrix-matrix multiplication.
!>
!> \author  Alfio Lazzaro
!>
!> <b>Modification history:</b>
!>  - 2016-08    Code organization (Alfio Lazzaro).
! **************************************************************************************************
MODULE dbcsr_multiply
   USE acc_stream,                      ONLY: acc_stream_associated,&
                                              acc_stream_create,&
                                              acc_stream_destroy
   USE array_types,                     ONLY: array_data,&
                                              array_equality,&
                                              array_hold,&
                                              array_i1d_obj,&
                                              array_nullify,&
                                              array_release
   USE dbcsr_config,                    ONLY: default_resize_factor,&
                                              has_acc,&
                                              num_layers_3D,&
                                              use_mpi_exp
   USE dbcsr_data_methods,              ONLY: dbcsr_data_set_size_referenced,&
                                              dbcsr_scalar,&
                                              dbcsr_scalar_are_equal,&
                                              dbcsr_scalar_fill_all,&
                                              dbcsr_scalar_one,&
                                              dbcsr_scalar_set_type,&
                                              dbcsr_scalar_zero
   USE dbcsr_dist_methods,              ONLY: &
        dbcsr_distribution_col_dist, dbcsr_distribution_get_num_images_1d, &
        dbcsr_distribution_has_col_clusters, dbcsr_distribution_has_row_clusters, &
        dbcsr_distribution_has_threads, dbcsr_distribution_hold, dbcsr_distribution_make_threads, &
        dbcsr_distribution_mp, dbcsr_distribution_ncols, dbcsr_distribution_no_threads, &
        dbcsr_distribution_nrows, dbcsr_distribution_num_local_col_clusters, &
        dbcsr_distribution_num_local_row_clusters, dbcsr_distribution_release, &
        dbcsr_distribution_row_dist
   USE dbcsr_dist_operations,           ONLY: dbcsr_create_image_dist,&
                                              dbcsr_make_dists_dense,&
                                              dbcsr_reset_locals,&
                                              make_sizes_dense
   USE dbcsr_index_operations,          ONLY: dbcsr_make_index_canonical
   USE dbcsr_io,                        ONLY: dbcsr_print
   USE dbcsr_mem_methods,               ONLY: dbcsr_mempool_clear,&
                                              dbcsr_mempool_destruct,&
                                              dbcsr_mempool_ensure_capacity,&
                                              dbcsr_memtype_setup
   USE dbcsr_methods,                   ONLY: &
        dbcsr_col_block_offsets, dbcsr_col_block_sizes, dbcsr_destroy_array, dbcsr_distribution, &
        dbcsr_get_data_type, dbcsr_get_matrix_type, dbcsr_has_symmetry, dbcsr_image_dist_release, &
        dbcsr_init, dbcsr_nblkcols_total, dbcsr_nfullcols_total, dbcsr_nfullrows_total, &
        dbcsr_release, dbcsr_release_locals, dbcsr_row_block_offsets
   USE dbcsr_mm_3D,                     ONLY: &
        buffers_1, buffers_2, buffers_orig, buffers_release, dbcsr_make_buffers, &
        make_layers_3D_C_reduction, multiply_3D, multiply_clusters, release_layers_3D_C_reduction
   USE dbcsr_mm_cannon,                 ONLY: make_m2s,&
                                              multiply_cannon
   USE dbcsr_mm_common,                 ONLY: &
        dbcsr_mpi_statistics, max_memory, max_nblocks, memtype_abpanel_1, memtype_abpanel_2, &
        memtype_mpi_buffer, memtype_product_wm, memtype_trsbuffer_1, memtype_trsbuffer_2, &
        num_multiplications, stream_1, stream_2
   USE dbcsr_mm_multrec,                ONLY: dbcsr_mm_multrec_lib_finalize,&
                                              dbcsr_mm_multrec_lib_init
   USE dbcsr_mp_methods,                ONLY: dbcsr_mp_group,&
                                              dbcsr_mp_npcols,&
                                              dbcsr_mp_nprows,&
                                              dbcsr_mp_numnodes
   USE dbcsr_operations,                ONLY: dbcsr_conjg,&
                                              dbcsr_copy,&
                                              dbcsr_filter,&
                                              dbcsr_may_be_dense,&
                                              dbcsr_scale
   USE dbcsr_toollib,                   ONLY: uppercase
   USE dbcsr_transformations,           ONLY: dbcsr_make_dense,&
                                              dbcsr_make_undense,&
                                              dbcsr_make_untransposed_blocks,&
                                              dbcsr_new_transposed
   USE dbcsr_types,                     ONLY: &
        dbcsr_2d_array_type, dbcsr_conjugate_transpose, dbcsr_distribution_obj, &
        dbcsr_imagedistribution_obj, dbcsr_mp_obj, dbcsr_mpi_size_limits, dbcsr_no_transpose, &
        dbcsr_obj, dbcsr_scalar_type, dbcsr_transpose, dbcsr_type_antisymmetric
   USE dbcsr_util,                      ONLY: dbcsr_checksum,&
                                              dbcsr_verify_matrix
   USE dbcsr_work_operations,           ONLY: dbcsr_add_wm_from_matrix,&
                                              dbcsr_finalize,&
                                              dbcsr_work_create
   USE kinds,                           ONLY: dp,&
                                              int_8,&
                                              real_8
   USE machine,                         ONLY: default_output_unit
   USE message_passing,                 ONLY: mp_max,&
                                              mp_min,&
                                              mp_sum
#include "../../base/base_uses.f90"

!$ USE OMP_LIB, ONLY: omp_get_thread_num, omp_get_num_threads

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_multiply'
   LOGICAL, PARAMETER :: debug_mod = .FALSE.
   LOGICAL, PARAMETER :: careful_mod = .FALSE.

   REAL, PRIVATE, SAVE :: marketing_flops = 0

   PUBLIC :: dbcsr_multiply_lib_init, dbcsr_multiply_lib_finalize
   PUBLIC :: dbcsr_multiply_clear_mempools
   PUBLIC :: dbcsr_multiply_internal

CONTAINS

! **************************************************************************************************
!> \brief Initialize the library
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbcsr_multiply_lib_init()

      INTEGER                                            :: ithread, nthreads

      nthreads = 1; ithread = 0
!$    nthreads = OMP_GET_NUM_THREADS(); ithread = OMP_GET_THREAD_NUM()

      CALL dbcsr_mm_multrec_lib_init()

!$OMP     MASTER
      dbcsr_mpi_statistics%last_mpi_ranks_used = 0
      dbcsr_mpi_statistics%nimages = -1
      dbcsr_mpi_statistics%nexchanged = 0
      dbcsr_mpi_statistics%nfiltered = 0
      dbcsr_mpi_statistics%data_size = 0
      dbcsr_mpi_statistics%data_size(:, 2) = HUGE(dbcsr_mpi_statistics%data_size(1, 2))
      dbcsr_mpi_statistics%data_size_breakdown = 0

      marketing_flops = 0
      max_memory = 0
      ALLOCATE (memtype_product_wm(0:nthreads-1))
!$OMP     END MASTER
!$OMP     BARRIER

      ! Each thread has its own working-matrix and its own mempool
      ALLOCATE (memtype_product_wm(ithread)%p)
      CALL dbcsr_memtype_setup(memtype_product_wm(ithread)%p, has_pool=.TRUE.)
      CALL dbcsr_mempool_ensure_capacity(memtype_product_wm(ithread)%p%pool, capacity=1)
   END SUBROUTINE dbcsr_multiply_lib_init

! **************************************************************************************************
!> \brief Finalize the library
!> \param group ...
!> \param output_unit ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbcsr_multiply_lib_finalize(group, output_unit)
      INTEGER, INTENT(IN)                                :: group, output_unit

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_lib_finalize', &
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

!$OMP      BARRIER
!$OMP      MASTER
      DEALLOCATE (memtype_product_wm)

      ! Deallocate buffers
      CALL buffers_release(buffers_orig)
      CALL buffers_release(buffers_1)
      CALL buffers_release(buffers_2)

      ! Release 3D communicators
      CALL release_layers_3D_C_reduction()

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

         isqrt = NINT(SQRT(REAL(dbcsr_mpi_statistics%last_mpi_ranks_used, KIND=real_8)))
         isqrt2 = NINT(SQRT(REAL(dbcsr_mpi_statistics%last_mpi_ranks_used*2, KIND=real_8)))
         IF (isqrt*isqrt .NE. dbcsr_mpi_statistics%last_mpi_ranks_used) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("-", 79)
            WRITE (UNIT=output_unit, FMT="(T2,A)") &
               "Warning: using a non-square number of MPI ranks might lead to poor performance."
            WRITE (UNIT=output_unit, FMT="(T2,A,I0)") &
               "         used ranks: ", dbcsr_mpi_statistics%last_mpi_ranks_used
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
   END SUBROUTINE dbcsr_multiply_lib_finalize

! **************************************************************************************************
!> \brief Deallocate memory contained in mempools
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbcsr_multiply_clear_mempools()

      INTEGER                                            :: ithread

      ithread = 0
!$    ithread = omp_get_thread_num()

      ! Each thread has its own working-matrix and its own mempool
      IF (ASSOCIATED(memtype_product_wm(ithread)%p%pool)) &
         CALL dbcsr_mempool_clear(memtype_product_wm(ithread)%p%pool)

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
   END SUBROUTINE dbcsr_multiply_clear_mempools

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
   SUBROUTINE dbcsr_multiply_internal(transa, transb, &
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

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_internal', &
         routineP = moduleN//':'//routineN
      LOGICAL, PARAMETER                                 :: dbg = .FALSE., &
                                                            use_list_indexing = .TRUE., &
                                                            use_local_indexing = .TRUE.
      REAL(real_8), PARAMETER                            :: make_dense_occ_thresh = 1.0_dp

      CHARACTER                                          :: transa_l, transb_l
      INTEGER :: comm, f_col, f_k, f_row, handle, handle2, ithread, l_col, l_k, l_row, &
         nimages_left_rows, nimages_match, nimages_right_cols, npcols, nprows, numnodes, &
         output_unit
      INTEGER(KIND=int_8)                                :: my_flop
      LOGICAL :: ab_dense, has_clusters, keep_product_data, keep_sparsity, product_reindex, &
         release_tdist, transpose_left, transpose_right, use_dense_mult
      REAL(KIND=dp)                                      :: cs
      TYPE(array_i1d_obj) :: dense_col_sizes, dense_k_sizes, dense_row_sizes, k_vmap, m_map, &
         n_map, old_product_col_blk_offsets, old_product_col_blk_sizes, &
         old_product_row_blk_offsets, old_product_row_blk_sizes
      TYPE(dbcsr_2d_array_type), POINTER                 :: m2s_left, m2s_right
      TYPE(dbcsr_distribution_obj)                       :: dense_product_distribution, &
                                                            old_product_distribution
      TYPE(dbcsr_imagedistribution_obj)                  :: dense_rdist_left, dense_rdist_right, &
                                                            rdist_left, rdist_right
      TYPE(dbcsr_mp_obj)                                 :: mp_obj
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
         transpose_left = .FALSE.
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
         transpose_left = .TRUE.
      CASE (dbcsr_conjugate_transpose)
         CALL dbcsr_init(matrix_left)
         CALL dbcsr_new_transposed(matrix_left, matrix_a, &
                                   shallow_data_copy=.FALSE., redistribute=.FALSE., &
                                   transpose_distribution=.FALSE.)
         CALL dbcsr_conjg(matrix_left)
         transpose_left = .TRUE.
      CASE DEFAULT
         CPABORT("wrong transa_l = "//transa_l)
      END SELECT

      SELECT CASE (transb_l)
      CASE (dbcsr_no_transpose)
         matrix_right = matrix_b
         transpose_right = .FALSE.
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
         transpose_right = .TRUE.
      CASE (dbcsr_conjugate_transpose)
         CALL dbcsr_init(matrix_right)
         CALL dbcsr_new_transposed(matrix_right, matrix_b, &
                                   shallow_data_copy=.FALSE., redistribute=.FALSE., &
                                   transpose_distribution=.FALSE.)
         CALL dbcsr_conjg(matrix_right)
         transpose_right = .TRUE.
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
      mp_obj = dbcsr_distribution_mp(matrix_c%m%dist)
      numnodes = dbcsr_mp_numnodes(mp_obj)
      nprows = dbcsr_mp_nprows(mp_obj)
      npcols = dbcsr_mp_npcols(mp_obj)
      !
      has_clusters = dbcsr_distribution_has_row_clusters(matrix_left%m%dist) .OR. &
                     dbcsr_distribution_has_col_clusters(matrix_left%m%dist) .OR. &
                     dbcsr_distribution_has_row_clusters(matrix_right%m%dist) .OR. &
                     dbcsr_distribution_has_col_clusters(matrix_right%m%dist)
      !
      ! 3D layers
      CALL make_layers_3D_C_reduction(num_layers_3D, mp_obj, has_clusters)
      !
      ! No dense multiplication and local indexing when 2D clustering is used.
      IF (use_mpi_exp) THEN
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
      dbcsr_mpi_statistics%last_mpi_ranks_used = numnodes
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
      IF (f_row .GT. l_row .AND. l_row .GT. 0) &
         CPABORT("Last row smaller than first row")
      IF (f_col .GT. l_col .AND. l_col .GT. 0) &
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
      IF (use_mpi_exp .AND. dbcsr_distribution_has_row_clusters(matrix_left%m%dist)) THEN
         nimages_left_rows = dbcsr_distribution_num_local_row_clusters(matrix_left%m%dist)* &
                             dbcsr_mp_nprows(dbcsr_distribution_mp(matrix_left%m%dist))
      ELSE
         ! 1D clustering (only RMA support 2D clustering)
         nimages_left_rows = dbcsr_mp_nprows(dbcsr_distribution_mp(matrix_left%m%dist))
      ENDIF
      !
      IF (use_mpi_exp .AND. dbcsr_distribution_has_col_clusters(matrix_left%m%dist)) THEN
         nimages_match = dbcsr_distribution_num_local_col_clusters(matrix_left%m%dist)* &
                         dbcsr_mp_npcols(dbcsr_distribution_mp(matrix_left%m%dist))
      ELSE IF (use_mpi_exp .AND. dbcsr_distribution_has_row_clusters(matrix_right%m%dist)) THEN
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
      IF (use_mpi_exp .AND. dbcsr_distribution_has_col_clusters(matrix_right%m%dist)) THEN
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
                                   match_col_nbins=npcols, &
                                   match_col_pdist=dbcsr_distribution_col_dist(matrix_c%m%dist), &
                                   nimages_rows=nimages_match, &
                                   nimages_cols=nimages_right_cols)
      !
      CALL dbcsr_create_image_dist(rdist_left, matrix_left%m%dist, &
                                   match_row_pdist=dbcsr_distribution_row_dist(matrix_c%m%dist), &
                                   match_row_nbins=nprows, &
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
      IF (use_mpi_exp) THEN
         !
         max_nblocks = product_matrix%m%nblks
         !
         ! Left buffer images
         CALL dbcsr_make_buffers(matrix_left, rdist_left, .TRUE., &
                                 f_row, l_row, f_k, l_k, &
                                 PRESENT(filter_eps), &
                                 transpose_left, &
                                 has_acc)
         !
         ! Right buffer images
         CALL dbcsr_make_buffers(matrix_right, rdist_right, .FALSE., &
                                 f_k, l_k, f_col, l_col, &
                                 PRESENT(filter_eps), &
                                 transpose_right, &
                                 has_acc, &
                                 alpha)
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
      IF (use_mpi_exp) THEN
         IF (has_clusters) THEN
            CALL multiply_clusters(rdist_left, rdist_right, &
                                   matrix_left, matrix_right, product_matrix, &
                                   retain_sparsity=retain_sparsity, &
                                   filter_eps=filter_eps, &
                                   flop=my_flop, &
                                   keep_product_data=keep_product_data)
         ELSE
            CALL multiply_3D(rdist_left, rdist_right, &
                             matrix_left, matrix_right, product_matrix, &
                             retain_sparsity=retain_sparsity, &
                             filter_eps=filter_eps, &
                             flop=my_flop, keep_product_data=keep_product_data)
         ENDIF
      ELSE
         CALL multiply_cannon(m2s_left, m2s_right, product_matrix, &
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
      IF (transpose_left) CALL dbcsr_release(matrix_left)
      IF (transpose_right) CALL dbcsr_release(matrix_right)
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
      IF (.NOT. use_mpi_exp) THEN
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
      num_multiplications = num_multiplications+1
      CALL timestop(handle)
   END SUBROUTINE dbcsr_multiply_internal

END MODULE dbcsr_multiply

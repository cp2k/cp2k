!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   3D matrix-matrix multiplication.
!>
!> \author  Alfio Lazzaro
!>
!> <b>Modification history:</b>
!>  - 2016-08    Code organization (Alfio Lazzaro).
!>  - 2017-02    Remove clusters (Alfio Lazzaro).
! **************************************************************************************************
MODULE dbcsr_mm_3d
   USE acc_event,                       ONLY: acc_event_create,&
                                              acc_event_destroy,&
                                              acc_event_synchronize
   USE array_types,                     ONLY: array_data,&
                                              array_get,&
                                              array_size
   USE dbcsr_block_operations,          ONLY: dbcsr_block_conjg,&
                                              dbcsr_block_copy_aa,&
                                              dbcsr_block_real_neg,&
                                              dbcsr_block_scale,&
                                              dbcsr_block_transpose_aa,&
                                              dbcsr_data_clear,&
                                              dbcsr_data_set
   USE dbcsr_config,                    ONLY: dbcsr_cfg,&
                                              has_acc
   USE dbcsr_data_methods,              ONLY: &
        dbcsr_data_clear_pointer, dbcsr_data_ensure_size, dbcsr_data_exists, &
        dbcsr_data_get_memory_type, dbcsr_data_get_size, dbcsr_data_get_size_referenced, &
        dbcsr_data_get_type, dbcsr_data_host2dev, dbcsr_data_init, dbcsr_data_new, &
        dbcsr_data_release, dbcsr_data_set_pointer, dbcsr_data_set_size_referenced, &
        dbcsr_data_valid, dbcsr_scalar_are_equal, dbcsr_scalar_negative, dbcsr_scalar_one, &
        dbcsr_type_1d_to_2d
   USE dbcsr_data_types,                ONLY: dbcsr_datatype_sizeof
   USE dbcsr_dist_methods,              ONLY: &
        dbcsr_distribution_col_dist, dbcsr_distribution_has_threads, &
        dbcsr_distribution_local_cols, dbcsr_distribution_local_rows, &
        dbcsr_distribution_max_col_dist, dbcsr_distribution_max_row_dist, dbcsr_distribution_mp, &
        dbcsr_distribution_row_dist, dbcsr_distribution_thread_dist
   USE dbcsr_index_operations,          ONLY: dbcsr_repoint_index
   USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop
   USE dbcsr_mem_methods,               ONLY: dbcsr_mempool_limit_capacity
   USE dbcsr_methods,                   ONLY: &
        dbcsr_col_block_offsets, dbcsr_col_block_sizes, dbcsr_distribution, dbcsr_get_data_type, &
        dbcsr_has_symmetry, dbcsr_max_col_size, dbcsr_max_row_size, dbcsr_nblkcols_local, &
        dbcsr_nblkcols_total, dbcsr_nblkrows_local, dbcsr_nblkrows_total, dbcsr_nfullcols_local, &
        dbcsr_nfullcols_total, dbcsr_nfullrows_local, dbcsr_nfullrows_total, dbcsr_release, &
        dbcsr_row_block_offsets, dbcsr_row_block_sizes, dbcsr_valid_index
   USE dbcsr_mm_common,                 ONLY: &
        acc_transpose_blocks, calculate_norms, count_mpi_statistics, dbcsr_mm_multrec_type_p, &
        dbcsr_mpi_statistics, enumerate_blk_sizes, local_filter, max_memory, memtype_abpanel_1, &
        memtype_abpanel_2, memtype_mpi_buffer, memtype_mpi_product, memtype_product_wm, &
        memtype_trsbuffer_1, memtype_trsbuffer_2, product_matrix_size_guess, rec_sort_index, &
        setup_buffer_matrix
   USE dbcsr_mm_dist_operations,        ONLY: dbcsr_reset_locals,&
                                              dbcsr_reset_vlocals,&
                                              image_calculator
   USE dbcsr_mm_multrec,                ONLY: dbcsr_mm_multrec_dev2host_init,&
                                              dbcsr_mm_multrec_finalize,&
                                              dbcsr_mm_multrec_get_nblks,&
                                              dbcsr_mm_multrec_get_nze,&
                                              dbcsr_mm_multrec_init,&
                                              dbcsr_mm_multrec_multiply,&
                                              dbcsr_mm_multrec_phaseout,&
                                              dbcsr_mm_multrec_red3D
   USE dbcsr_mp_methods,                ONLY: &
        dbcsr_mp_grid_setup, dbcsr_mp_group, dbcsr_mp_has_subgroups, dbcsr_mp_my_col_group, &
        dbcsr_mp_my_row_group, dbcsr_mp_mynode, dbcsr_mp_mypcol, dbcsr_mp_myprow, dbcsr_mp_npcols, &
        dbcsr_mp_nprows, dbcsr_mp_numnodes, dbcsr_mp_pgrid
   USE dbcsr_mp_operations,             ONLY: dbcsr_isendrecv_any,&
                                              dbcsr_rget_any,&
                                              dbcsr_win_create_any,&
                                              hybrid_alltoall_any,&
                                              hybrid_alltoall_i1
   USE dbcsr_ptr_util,                  ONLY: ensure_array_size,&
                                              memory_deallocate
   USE dbcsr_types,                     ONLY: &
        dbcsr_2d_array_obj, dbcsr_data_obj, dbcsr_distribution_obj, dbcsr_imagedistribution_obj, &
        dbcsr_iterator, dbcsr_memtype_type, dbcsr_mp_obj, dbcsr_num_slots, dbcsr_scalar_type, &
        dbcsr_slot_blk_p, dbcsr_slot_col_i, dbcsr_slot_coo_l, dbcsr_slot_dense, &
        dbcsr_slot_home_pcol, dbcsr_slot_home_prow, dbcsr_slot_home_vpcol, dbcsr_slot_home_vprow, &
        dbcsr_slot_nblkrows_total, dbcsr_slot_nblks, dbcsr_slot_nfullcols_local, dbcsr_slot_nze, &
        dbcsr_slot_row_p, dbcsr_slot_size, dbcsr_slot_thr_c, dbcsr_slot_type, dbcsr_type, &
        dbcsr_type_complex_4, dbcsr_type_complex_8, dbcsr_type_int_4, dbcsr_type_no_symmetry, &
        dbcsr_type_real_4, dbcsr_type_real_8
   USE dbcsr_util,                      ONLY: find_block_of_element
   USE dbcsr_work_operations,           ONLY: dbcsr_add_wm_from_matrix,&
                                              dbcsr_create,&
                                              dbcsr_finalize,&
                                              dbcsr_work_create,&
                                              dbcsr_work_destroy
   USE kinds,                           ONLY: int_8,&
                                              real_8,&
                                              sp
   USE machine,                         ONLY: m_memory
   USE message_passing,                 ONLY: &
        mp_allgather, mp_alltoall, mp_comm_free, mp_comm_null, mp_comm_split_direct, &
        mp_iallgather, mp_isendrecv, mp_isum, mp_request_null, mp_rget, mp_wait, mp_waitall, &
        mp_win_create, mp_win_free, mp_win_lock_all, mp_win_unlock_all
#include "../../base/base_uses.f90"

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads,&
!$                    omp_set_lock, omp_unset_lock, omp_init_lock, omp_lock_kind, omp_destroy_lock

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_3d'
   LOGICAL, PARAMETER :: debug_mod = .FALSE.
   LOGICAL, PARAMETER :: careful_mod = .FALSE.

   TYPE dbcsr_buffer
      TYPE(dbcsr_data_obj)                        :: DATA, &
                                                     data_before_resize, &
                                                     trs_stackbuf
      INTEGER                                     :: vprow, vpcol
      INTEGER                                     :: grp = mp_comm_null, & ! Global communicator
                                                     subgrp = mp_comm_null ! Communicator for A and B
      INTEGER                                     :: data_win, meta_win
      INTEGER, DIMENSION(:), POINTER              :: meta => Null(), &
                                                     meta_before_resize => Null(), &
                                                     meta_red3D => Null()
      INTEGER, DIMENSION(2)                       :: get_requests
      INTEGER                                     :: meta_size
      INTEGER                                     :: num_layers_3D = 1
      INTEGER                                     :: coord3D = 1
      TYPE(dbcsr_type)                            :: matrix
      LOGICAL                                     :: is_valid = .FALSE., &
                                                     is_comm = .FALSE., &
                                                     has_rma_win = .FALSE.
   END TYPE dbcsr_buffer

   TYPE dbcsr_buffers
      TYPE(dbcsr_buffer)                 :: left, right
   END TYPE dbcsr_buffers

   TYPE mn_local_sizes
      INTEGER, DIMENSION(:), POINTER  :: sizes => Null()
   END TYPE mn_local_sizes

   TYPE dbcsr_layers_3D_C_reduction
      INTEGER                            :: grp = mp_comm_null, &
                                            grp3D = mp_comm_null, &
                                            rowgrp3D = mp_comm_null
      INTEGER                            :: num_layers_3D = 1, &
                                            max_num_layers_3D = 1, &
                                            side3D = HUGE(1)
      ! Use a buffer per each thread
      TYPE(dbcsr_data_obj), &
         DIMENSION(:), ALLOCATABLE       :: data_red3D
      INTEGER                            :: data_type
   END TYPE dbcsr_layers_3D_C_reduction

   ! Buffers
   TYPE(dbcsr_buffers), TARGET, SAVE :: buffers_win, &
                                        buffers_1, buffers_2

   INTEGER, PARAMETER, PRIVATE               :: idata = 1, &
                                                imeta = 2, &
                                                ilocal_proc = 1

   TYPE(dbcsr_layers_3D_C_reduction), SAVE :: layers_3D_C_reduction

   LOGICAL, DIMENSION(2), TARGET, SAVE, PRIVATE       :: do_win_create_left, &
                                                         do_win_create_right

   INTEGER, ALLOCATABLE, DIMENSION(:), PRIVATE :: left_total_row_counts

   INTEGER, ALLOCATABLE, DIMENSION(:, :), TARGET, PRIVATE  :: left_local_images_size, &
                                                              right_local_images_size
   INTEGER, DIMENSION(:, :, :, :), POINTER, PRIVATE        :: left_images_size => Null(), &
                                                              right_images_size => Null()

   INTEGER, ALLOCATABLE, DIMENSION(:), TARGET, PRIVATE :: g2l_map_cols, g2l_map_rows

   INTEGER, PRIVATE                :: request_count_rows
   INTEGER, DIMENSION(2), PRIVATE  :: requests
   INTEGER, DIMENSION(2), PRIVATE  :: requests_win_create

   INTEGER :: request_sync_mult = mp_request_null

   ! Buffers used in make_buffers
   TYPE(dbcsr_data_obj), TARGET, SAVE :: make_buffers_data_recv, make_buffers_data_send
   INTEGER, DIMENSION(:), POINTER     :: make_buffers_meta_recv => Null(), &
                                         make_buffers_meta_send => Null()

   PUBLIC :: multiply_3D
   PUBLIC :: release_layers_3d_C_reduction, buffers_release
   PUBLIC :: dbcsr_make_buffers, make_layers_3d_C_reduction
   PUBLIC :: request_sync_mult
   PUBLIC :: get_max_layers_3D

CONTAINS

! **************************************************************************************************
!> \brief Prepare orig images for MPI windows
!> \param matrix ...
!> \param imgdist ...
!> \param is_left ...
!> \param f_row ...
!> \param l_row ...
!> \param f_col ...
!> \param l_col ...
!> \param otf_filtering ...
!> \param alpha ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE dbcsr_make_buffers(matrix, imgdist, is_left, &
                                 f_row, l_row, f_col, l_col, &
                                 otf_filtering, alpha)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(dbcsr_imagedistribution_obj), INTENT(IN)      :: imgdist
      LOGICAL, INTENT(IN)                                :: is_left
      INTEGER, INTENT(IN)                                :: f_row, l_row, f_col, l_col
      LOGICAL, INTENT(IN)                                :: otf_filtering
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
                           f_row, l_row, f_col, l_col, &
                           otf_filtering, alpha)
      ELSE
         CALL make_buffers(matrix, imgdist, is_left, &
                           f_row, l_row, f_col, l_col, &
                           otf_filtering)
      ENDIF
   END SUBROUTINE dbcsr_make_buffers

! **************************************************************************************************
!> \brief Prepare orig images for MPI windows
!> \param matrix ...
!> \param imgdist ...
!> \param is_left ...
!> \param f_row ...
!> \param l_row ...
!> \param f_col ...
!> \param l_col ...
!> \param otf_filtering ...
!> \param scale_value ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE make_buffers(matrix, imgdist, is_left, &
                           f_row, l_row, f_col, l_col, &
                           otf_filtering, scale_value)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(dbcsr_imagedistribution_obj), INTENT(IN)      :: imgdist
      LOGICAL, INTENT(IN)                                :: is_left
      INTEGER, INTENT(IN)                                :: f_row, l_row, f_col, l_col
      LOGICAL, INTENT(IN)                                :: otf_filtering
      TYPE(dbcsr_scalar_type), INTENT(IN), OPTIONAL      :: scale_value

      CHARACTER(len=*), PARAMETER :: routineN = 'make_buffers', routineP = moduleN//':'//routineN

      INTEGER :: blk, blk_p, bp, col, col_img, col_size, data_type, dst_proc, f_col_f, f_row_f, &
         grp, handle, handle2, irequests, it, ithread, l_col_l, l_row_l, mynode, myt, &
         nblkcols_local, nblkrows_local, ncols_images, nimages, nprocs, nprocs_total, &
         nrows_images, nsymmetries, nthreads, nze, pcol, prow, row, row_img, row_size, size_index, &
         stored_col, stored_row, symmetry_i, tr_col_size, tr_row_size
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: img_nblks_cols, img_nblks_rows
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: local_images_displ, recv_displ_proc, &
                                                            recv_size_proc, send_displ_proc, &
                                                            send_size_proc
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: offset_data, offset_threads
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)        :: recv_displs, recv_sizes, send_displs, &
                                                            send_sizes
      INTEGER, DIMENSION(2)                              :: block_col_bounds, block_row_bounds
      INTEGER, DIMENSION(:), POINTER :: col_dist, col_img_dist, local_cols, local_g2l_map_cols, &
         local_g2l_map_rows, local_rows, meta_buffer_p, row_dist, row_img_dist, threads_dist
      INTEGER, DIMENSION(:, :), POINTER                  :: blacs2mpi, local_images_size
      INTEGER, DIMENSION(idata:imeta)                    :: my_size_recv, my_size_send
      INTEGER, POINTER                                   :: coli, rowi
      INTEGER, TARGET                                    :: mi, ui
      LOGICAL :: do_crop, do_part_crop_col, do_part_crop_f_col, do_part_crop_f_row, &
         do_part_crop_l_col, do_part_crop_l_row, do_part_crop_row, do_symmetry, tr
      LOGICAL, DIMENSION(:), POINTER                     :: do_win_create
      TYPE(dbcsr_buffer), POINTER                        :: buffer
      TYPE(dbcsr_data_obj)                               :: data_block
      TYPE(dbcsr_data_obj), POINTER                      :: data_buffer_p
      TYPE(dbcsr_distribution_obj)                       :: set_dist
      TYPE(dbcsr_iterator)                               :: iter
      TYPE(dbcsr_mp_obj)                                 :: mp_obj
      TYPE(dbcsr_scalar_type)                            :: scale_neg_one
      TYPE(dbcsr_type)                                   :: sm

!$    INTEGER(kind=omp_lock_kind), ALLOCATABLE, DIMENSION(:) :: locks

      CALL timeset(routineN, handle)
      ! Sync with previous multiplication
      IF (request_sync_mult .NE. mp_request_null) &
         CPABORT("Multiplications are not in sync!")
      !
      ! Take input values and check validity
      IF (.NOT. dbcsr_valid_index(matrix)) &
         CPABORT("Matrix not initialized.")
      sm = matrix
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
      mp_obj = dbcsr_distribution_mp(imgdist%i%main)
      CALL dbcsr_mp_grid_setup(mp_obj)
      nprocs_total = dbcsr_mp_numnodes(mp_obj)
      mynode = dbcsr_mp_mynode(mp_obj)
      grp = dbcsr_mp_group(mp_obj)
      blacs2mpi => dbcsr_mp_pgrid(mp_obj)
      IF (dbcsr_distribution_max_row_dist(set_dist) .GT. UBOUND(blacs2mpi, 1)) &
         CPABORT("Row distribution references unexistent processor rows")
      IF (dbcsr_distribution_max_col_dist(set_dist) .GT. UBOUND(blacs2mpi, 2)) &
         CPABORT("Col distribution references unexistent processor cols")
      ! Check threads configuration
      NULLIFY (threads_dist)
!$    IF (.NOT. dbcsr_distribution_has_threads(dbcsr_distribution(matrix))) &
!$       CPABORT("Thread distribution not defined")
!$    threads_dist => array_data(dbcsr_distribution_thread_dist(dbcsr_distribution(matrix)))
      IF (is_left) THEN
         IF (nrows_images .GT. 1) &
            CPABORT("Row nimages for left matrix is not 1!")
      ELSE
         IF (ncols_images .GT. 1) &
            CPABORT("Col nimages for right matrix is not 1!")
      ENDIF
      !
      ! Crop matrix
      do_crop = .FALSE.
      do_part_crop_row = .FALSE.
      do_part_crop_col = .FALSE.
      ! Set no limits
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
         do_crop = .TRUE.
         !
         ! Convert bounds to block addressing
         do_part_crop_f_row = .FALSE.
         IF (f_row .EQ. 0) THEN
            block_row_bounds(1) = 1
         ELSE
            CALL find_block_of_element(f_row, block_row_bounds(1), &
                                       dbcsr_nblkrows_total(matrix), &
                                       dbcsr_row_block_offsets(matrix), &
                                       hint=0)
            do_part_crop_f_row = array_get(dbcsr_row_block_offsets(matrix), block_row_bounds(1)) .NE. f_row
            IF (do_part_crop_f_row) THEN
               ! Block offset of last cleared row
               f_row_f = f_row-array_get(dbcsr_row_block_offsets(matrix), block_row_bounds(1))
            ENDIF
         ENDIF
         !
         do_part_crop_l_row = .FALSE.
         IF (l_row .EQ. 0) THEN
            block_row_bounds(2) = dbcsr_nblkrows_total(matrix)
         ELSE
            CALL find_block_of_element(l_row, block_row_bounds(2), &
                                       dbcsr_nblkrows_total(matrix), &
                                       dbcsr_row_block_offsets(matrix), &
                                       hint=0)
            do_part_crop_l_row = (array_get(dbcsr_row_block_offsets(matrix), block_row_bounds(2)+1)-1) .NE. l_row
            IF (do_part_crop_l_row) THEN
               ! Block offset of first cleared row
               l_row_l = 2+l_row-array_get(dbcsr_row_block_offsets(matrix), block_row_bounds(2))
            ENDIF
         ENDIF
         do_part_crop_row = do_part_crop_f_row .OR. do_part_crop_l_row
         !
         do_part_crop_f_col = .FALSE.
         IF (f_col .EQ. 0) THEN
            block_col_bounds(1) = 1
         ELSE
            CALL find_block_of_element(f_col, block_col_bounds(1), &
                                       dbcsr_nblkcols_total(matrix), &
                                       dbcsr_col_block_offsets(matrix), &
                                       hint=0)
            do_part_crop_f_col = array_get(dbcsr_col_block_offsets(matrix), block_col_bounds(1)) .NE. f_col
            IF (do_part_crop_f_col) THEN
               ! Block offset of last cleared col
               f_col_f = f_col-array_get(dbcsr_col_block_offsets(matrix), block_col_bounds(1))
            ENDIF
         ENDIF
         !
         do_part_crop_l_col = .FALSE.
         IF (l_col .EQ. 0) THEN
            block_col_bounds(2) = dbcsr_nblkcols_total(matrix)
         ELSE
            CALL find_block_of_element(l_col, block_col_bounds(2), &
                                       dbcsr_nblkcols_total(matrix), &
                                       dbcsr_col_block_offsets(matrix), &
                                       hint=0)
            do_part_crop_l_col = (array_get(dbcsr_col_block_offsets(matrix), block_col_bounds(2)+1)-1) .NE. l_col
            IF (do_part_crop_l_col) THEN
               ! Block offset of first cleared col
               l_col_l = 2+l_col-array_get(dbcsr_col_block_offsets(matrix), block_col_bounds(2))
            ENDIF
         ENDIF
         do_part_crop_col = do_part_crop_f_col .OR. do_part_crop_l_col
      ENDIF
      !
      IF (dbcsr_has_symmetry(matrix)) THEN
         nsymmetries = 2
         do_symmetry = .TRUE.
      ELSE
         nsymmetries = 1
         do_symmetry = .FALSE.
      ENDIF
      !
      IF (is_left) THEN
         nimages = ncols_images
         buffer => buffers_win%left
         nprocs = dbcsr_mp_npcols(mp_obj)
         ALLOCATE (left_images_size(idata:imeta, &
                                    nimages, &
                                    MAX(1, dbcsr_mp_nprows(mp_obj)/layers_3D_C_reduction%side3D), &
                                    0:nprocs-1))
         ALLOCATE (left_local_images_size(idata:imeta, nimages))
         local_images_size => left_local_images_size
         irequests = 1
         !
         ! Count the maximum possible multiplies per row for on-the-fly filtering
         IF (otf_filtering) THEN
            ALLOCATE (left_total_row_counts(nblkrows_local))
            left_total_row_counts = 0
         ENDIF
         do_win_create => do_win_create_left
      ELSE
         nimages = nrows_images
         buffer => buffers_win%right
         nprocs = dbcsr_mp_nprows(mp_obj)
         ALLOCATE (right_images_size(idata:imeta, &
                                     nimages, &
                                     MAX(1, dbcsr_mp_npcols(mp_obj)/layers_3D_C_reduction%side3D), &
                                     0:nprocs-1))
         ALLOCATE (right_local_images_size(idata:imeta, nimages))
         local_images_size => right_local_images_size
         irequests = 2
         do_win_create => do_win_create_right
      ENDIF
      !
      ! 3D communicator
      CALL make_layers_3D_AB(layers_3D_C_reduction%num_layers_3D, &
                             layers_3D_C_reduction%side3D, &
                             mp_obj, is_left, buffer)
      !
      ! Evaluate maps for global -> local indexing (g2l_map_rows, g2l_map_cols)
      ! Count the number of blocks per row/column (img_nblks_rows, img_nblks_cols)
      IF (is_left) THEN
         ALLOCATE (g2l_map_rows(sm%nblkrows_total))
         local_g2l_map_rows => g2l_map_rows
         ALLOCATE (local_g2l_map_cols(sm%nblkcols_total))
         ALLOCATE (img_nblks_rows(1), img_nblks_cols(nimages))
      ELSE
         ALLOCATE (g2l_map_cols(sm%nblkcols_total))
         local_g2l_map_cols => g2l_map_cols
         ALLOCATE (local_g2l_map_rows(sm%nblkrows_total))
         ALLOCATE (img_nblks_rows(nimages), img_nblks_cols(1))
      ENDIF
      !
      local_g2l_map_rows(:) = 0
      IF (nrows_images .EQ. 1) THEN
         img_nblks_rows(1) = nblkrows_local
         DO row = 1, nblkrows_local
            local_g2l_map_rows(local_rows(row)) = row
         ENDDO
      ELSE
         img_nblks_rows(:) = 0
         DO row = 1, nblkrows_local
            row_img = row_img_dist(local_rows(row))
            ui = MOD(row_img-1, nrows_images)+1
            img_nblks_rows(ui) = img_nblks_rows(ui)+1
            local_g2l_map_rows(local_rows(row)) = img_nblks_rows(ui)
         ENDDO
      ENDIF
      !
      local_g2l_map_cols(:) = 0
      IF (ncols_images .EQ. 1) THEN
         img_nblks_cols(1) = nblkcols_local
         DO col = 1, nblkcols_local
            local_g2l_map_cols(local_cols(col)) = col
         ENDDO
      ELSE
         img_nblks_cols(:) = 0
         DO col = 1, nblkcols_local
            col_img = col_img_dist(local_cols(col))
            ui = MOD(col_img-1, ncols_images)+1
            img_nblks_cols(ui) = img_nblks_cols(ui)+1
            local_g2l_map_cols(local_cols(col)) = img_nblks_cols(ui)
         ENDDO
      ENDIF
      !
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP PRIVATE (ithread,myt,iter,row,col,blk,row_size,col_size,&
!$OMP          stored_row,stored_col,blk_p,bp,tr,&
!$OMP          nze,symmetry_i,row_img,col_img,rowi,coli,&
!$OMP          tr_row_size,tr_col_size,prow,pcol,dst_proc,&
!$OMP          data_buffer_p,meta_buffer_p,&
!$OMP          mi,ui,it,data_block) &
!$OMP SHARED (nthreads,send_sizes,offset_data,matrix,nsymmetries,do_symmetry,&
!$OMP         row_img_dist,col_img_dist,imgdist,row_dist,col_dist,&
!$OMP         is_left,my_size_send,my_size_recv,nimages,&
!$OMP         local_images_size,data_type,memtype_mpi_buffer,sm,&
!$OMP         img_nblks_cols,img_nblks_rows,mynode,offset_threads,&
!$OMP         local_g2l_map_cols,local_g2l_map_rows,recv_sizes,grp,make_buffers_meta_send,&
!$OMP         scale_value,scale_neg_one,make_buffers_data_send,make_buffers_data_recv,&
!$OMP         size_index,recv_displ_proc,recv_size_proc,send_size_proc,&
!$OMP         mp_obj,threads_dist,make_buffers_meta_recv,nrows_images,ncols_images,&
!$OMP         locks,blacs2mpi,send_displ_proc,recv_displs,send_displs,&
!$OMP         left_images_size,right_images_size,local_images_displ,requests,&
!$OMP         buffer,left_total_row_counts,otf_filtering,&
!$OMP         irequests,do_win_create,handle2,nprocs_total,&
!$OMP         do_crop,do_part_crop_row,do_part_crop_col,block_row_bounds,block_col_bounds,&
!$OMP         do_part_crop_f_row,do_part_crop_l_row,do_part_crop_f_col,do_part_crop_l_col,&
!$OMP         f_row_f,l_row_l,f_col_f,l_col_l,requests_win_create)
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
      ALLOCATE (send_sizes(idata:imeta, 0:nthreads-1, &
                           nimages, 0:nprocs_total-1))
      send_sizes(:, :, :, :) = 0
      !
      size_index = 0
!$    IF (is_left) THEN
!$       size_index = nthreads+1
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
      ! Take data and meta dimensions per each thread, image, proc
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
            ! Apply cropping
            IF (do_crop) THEN
               IF (stored_row .LT. block_row_bounds(1)) CYCLE
               IF (stored_row .GT. block_row_bounds(2)) CYCLE
               IF (stored_col .LT. block_col_bounds(1)) CYCLE
               IF (stored_col .GT. block_col_bounds(2)) CYCLE
            ENDIF
            row_img = row_img_dist(stored_row)
            col_img = col_img_dist(stored_col)
            CALL image_calculator(imgdist, &
                                  prow=prow, pcol=pcol, &
                                  rowi=rowi, coli=coli, &
                                  myprow=row_dist(stored_row), myrowi=row_img, &
                                  mypcol=col_dist(stored_col), mycoli=col_img, &
                                  shifting='0')
            dst_proc = blacs2mpi(prow, pcol)
!$          IF (is_left .AND. do_symmetry) THEN
!$             myt = threads_dist(stored_row)
!$          ENDIF
!$OMP ATOMIC
            send_sizes(imeta, myt, ui, dst_proc) = &
               send_sizes(imeta, myt, ui, dst_proc)+3
!$OMP ATOMIC
            send_sizes(idata, myt, ui, dst_proc) = &
               send_sizes(idata, myt, ui, dst_proc)+nze
         ENDDO ! symmetry_i
      ENDDO
      CALL dbcsr_iterator_stop(iter)
!$OMP BARRIER
!$OMP MASTER
      ! Exchange refs
      ALLOCATE (recv_sizes(idata:imeta, 0:nthreads-1, &
                           nimages, 0:nprocs_total-1))
      CALL timeset(routineN//"_sizes", handle2)
      CALL mp_alltoall(send_sizes(:, :, :, :), &
                       recv_sizes(:, :, :, :), &
                       2*nimages*nthreads, grp)
      CALL timestop(handle2)
      !
      ! Evaluate the local size for each image, accumulating over threads and procs.
      ! Take the local displacement for each image.
      ! Note that diplacement starts at zero.
      my_size_recv(:) = 0
      local_images_size(:, :) = 0
      ALLOCATE (local_images_displ(idata:imeta, nimages))
      DO ui = 1, nimages
         local_images_displ(:, ui) = my_size_recv(:)
         DO dst_proc = 0, nprocs_total-1
            DO it = 0, nthreads-1
               local_images_size(:, ui) = local_images_size(:, ui)+ &
                                          recv_sizes(:, it, ui, dst_proc)
            ENDDO
         ENDDO
         IF (local_images_size(imeta, ui) .EQ. 0) CYCLE
         ! Include stats slots for threads indices
         local_images_size(imeta, ui) = local_images_size(imeta, ui)+size_index
         my_size_recv(:) = my_size_recv(:)+local_images_size(:, ui)
      ENDDO
      !
      ! Exchange sizes
      IF (is_left) THEN
         CALL mp_iallgather(local_images_size, left_images_size, buffer%subgrp, requests(irequests))
      ELSE
         CALL mp_iallgather(local_images_size, right_images_size, buffer%subgrp, requests(irequests))
      ENDIF
      !
      ! Allocate data and meta buffers
      do_win_create(:) = .NOT. buffer%has_rma_win
      IF (buffer%has_rma_win) THEN
         IF (buffer%grp .NE. grp .OR. dbcsr_data_get_type(buffer%data) .NE. data_type) THEN
            do_win_create(:) = .TRUE.
         ENDIF
      ENDIF
      CALL buffer_init(buffer, data_type, &
                       my_size_recv(idata), my_size_recv(imeta), &
                       data_memory_type=memtype_mpi_buffer)
      buffer%grp = grp
      !
      ! Set send and recv buffers sizes and displacements for each proc.
      ! Accumulate over images and threads.
      ! Here diplacement starts at one.
      ALLOCATE (send_displs(idata:imeta, 0:nthreads-1, &
                            nimages, 0:nprocs_total-1))
      ! Displs for local data arrangement, starting at one.
      ALLOCATE (recv_displs(idata:imeta, 0:nthreads-1, &
                            nimages, 0:nprocs_total-1))
      ! Here displacement starts at zero.
      ALLOCATE (send_size_proc(idata:imeta, 0:nprocs_total-1))
      ALLOCATE (recv_size_proc(idata:imeta, 0:nprocs_total-1))
      ALLOCATE (send_displ_proc(idata:imeta, 0:nprocs_total-1))
      ALLOCATE (recv_displ_proc(idata:imeta, 0:nprocs_total-1))
      my_size_send(:) = 1
      my_size_recv(:) = 1
      DO dst_proc = 0, nprocs_total-1
         send_displ_proc(:, dst_proc) = my_size_send(:)-1
         recv_displ_proc(:, dst_proc) = my_size_recv(:)-1
         ! Avoid communication of local data
         IF (dst_proc .NE. mynode) THEN
            DO ui = 1, nimages
               DO it = 0, nthreads-1
                  send_displs(:, it, ui, dst_proc) = my_size_send(:)
                  recv_displs(:, it, ui, dst_proc) = my_size_recv(:)
                  my_size_send(:) = my_size_send(:)+send_sizes(:, it, ui, dst_proc)
                  my_size_recv(:) = my_size_recv(:)+recv_sizes(:, it, ui, dst_proc)
               ENDDO
            ENDDO
         ELSE
            ! Reset all
            send_displs(:, :, :, dst_proc) = 0
            recv_displs(:, :, :, dst_proc) = 0
         ENDIF
         send_size_proc(:, dst_proc) = my_size_send(:)-send_displ_proc(:, dst_proc)-1
         recv_size_proc(:, dst_proc) = my_size_recv(:)-recv_displ_proc(:, dst_proc)-1
      ENDDO
      !
      ! Allocate data/meta to send
      IF (dbcsr_data_valid(make_buffers_data_send)) THEN
         IF (dbcsr_data_get_type(make_buffers_data_send) .NE. data_type) THEN
            CALL dbcsr_data_release(make_buffers_data_send)
         ENDIF
      ENDIF
      IF (dbcsr_data_valid(make_buffers_data_send)) THEN
         CALL dbcsr_data_ensure_size(make_buffers_data_send, my_size_send(idata)-1, nocopy=.TRUE.)
      ELSE
         CALL dbcsr_data_init(make_buffers_data_send)
         CALL dbcsr_data_new(make_buffers_data_send, data_type, my_size_send(idata)-1, &
                             memory_type=memtype_mpi_buffer)
      ENDIF
      CALL ensure_array_size(make_buffers_meta_send, ub=my_size_send(imeta)-1, &
                             nocopy=.TRUE., memory_type=memtype_mpi_buffer)
      ! Displs for data offset
      ALLOCATE (offset_threads(idata:imeta, 0:nthreads-1, nimages))
      offset_threads(:, :, :) = 0
      ! Set offset for local data
      ALLOCATE (offset_data(0:nthreads-1, nimages, 0:nprocs_total-1))
      offset_data(:, :, :) = 1
      ! Evaluate local displs
      DO ui = 1, nimages
         IF (local_images_size(imeta, ui) .EQ. 0) CYCLE
         offset_threads(:, 0, ui) = 0
         DO it = 1, nthreads-1
            offset_threads(:, it, ui) = offset_threads(:, it-1, ui)
            DO dst_proc = 0, nprocs_total-1
               offset_threads(:, it, ui) = offset_threads(:, it, ui)+ &
                                           recv_sizes(:, it-1, ui, dst_proc)
            ENDDO
         ENDDO
         ! Fill meta indices for threads
!$       IF (is_left) THEN
!$          buffer%meta(local_images_displ(imeta, ui)+1:local_images_displ(imeta, ui)+nthreads) = &
!$             offset_threads(imeta, :, ui)/3
!$          buffer%meta(local_images_displ(imeta, ui)+size_index) = &
!$             (local_images_size(imeta, ui)-size_index)/3
!$       ENDIF
         offset_threads(imeta, :, ui) = offset_threads(imeta, :, ui)+local_images_displ(imeta, ui)+size_index+1
         send_displs(:, :, ui, mynode) = offset_threads(:, :, ui)
         !
         ! Allow ordering by proc for insertion
         DO dst_proc = 0, mynode-1
            DO it = 0, nthreads-1
               send_displs(:, it, ui, mynode) = send_displs(:, it, ui, mynode)+ &
                                                recv_sizes(:, it, ui, dst_proc)
            ENDDO
         ENDDO
         offset_data(:, ui, mynode) = send_displs(idata, :, ui, mynode)+1
         send_displs(idata, :, ui, mynode) = send_displs(idata, :, ui, mynode)+local_images_displ(idata, ui)+1
      ENDDO
!$OMP END MASTER
!$OMP BARRIER
      !
      IF (do_part_crop_row .OR. do_part_crop_col) THEN
         CALL dbcsr_data_init(data_block)
         CALL dbcsr_data_new(data_block, dbcsr_type_1d_to_2d(data_type))
      ENDIF
      !
      ! Copy data and meta in the buffers
      CALL timeset(routineN//"_pack", handle2)
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
            ! Apply cropping
            IF (do_crop) THEN
               IF (stored_row .LT. block_row_bounds(1)) CYCLE
               IF (stored_row .GT. block_row_bounds(2)) CYCLE
               IF (stored_col .LT. block_col_bounds(1)) CYCLE
               IF (stored_col .GT. block_col_bounds(2)) CYCLE
            ENDIF
            row_img = row_img_dist(stored_row)
            col_img = col_img_dist(stored_col)
            CALL image_calculator(imgdist, &
                                  prow=prow, pcol=pcol, &
                                  rowi=rowi, coli=coli, &
                                  myprow=row_dist(stored_row), myrowi=row_img, &
                                  mypcol=col_dist(stored_col), mycoli=col_img, &
                                  shifting='0')
            dst_proc = blacs2mpi(prow, pcol)
            IF (dst_proc .EQ. mynode) THEN
               data_buffer_p => buffer%data
               meta_buffer_p => buffer%meta
            ELSE
               data_buffer_p => make_buffers_data_send
               meta_buffer_p => make_buffers_meta_send
            ENDIF
!$          IF (is_left .AND. do_symmetry) THEN
!$             myt = threads_dist(stored_row)
!$             call omp_set_lock(locks(myt))
!$          ENDIF
            IF (tr) THEN
               CALL dbcsr_block_transpose_aa(data_buffer_p, sm%data_area, tr_row_size, tr_col_size, &
                                             send_displs(idata, myt, ui, dst_proc), bp, &
                                             scale_value)
               IF (sm%negate_real .AND. sm%negate_imaginary) THEN
                  CALL dbcsr_block_scale(data_buffer_p, scale=scale_neg_one, &
                                         row_size=nze, col_size=1, &
                                         lb=send_displs(idata, myt, ui, dst_proc))
               ELSEIF (sm%negate_real) THEN
                  CALL dbcsr_block_real_neg(data_buffer_p, row_size=nze, col_size=1, &
                                            lb=send_displs(idata, myt, ui, dst_proc))
               ELSEIF (sm%negate_imaginary) THEN
                  CALL dbcsr_block_conjg(data_buffer_p, row_size=nze, col_size=1, &
                                         lb=send_displs(idata, myt, ui, dst_proc))
               ENDIF
            ELSE
               CALL dbcsr_block_copy_aa(data_buffer_p, sm%data_area, row_size, col_size, &
                                        send_displs(idata, myt, ui, dst_proc), bp, &
                                        scale_value)
            ENDIF
            !
            ! Apply cropping for partial blocks
            IF (do_part_crop_row .OR. do_part_crop_col) THEN
               CALL dbcsr_data_set_pointer( &
                  area=data_block, &
                  rsize=row_size, &
                  csize=col_size, &
                  pointee=data_buffer_p, &
                  source_lb=send_displs(idata, myt, ui, dst_proc))
               IF (do_part_crop_row) THEN
                  IF (do_part_crop_f_row .AND. stored_row .EQ. block_row_bounds(1)) THEN
                     CALL dbcsr_data_clear(data_block, ub=f_row_f)
                  ENDIF
                  IF (do_part_crop_l_row .AND. stored_row .EQ. block_row_bounds(2)) THEN
                     CALL dbcsr_data_clear(data_block, lb=l_row_l)
                  ENDIF
               ENDIF
               IF (do_part_crop_col) THEN
                  IF (do_part_crop_f_col .AND. stored_col .EQ. block_col_bounds(1)) THEN
                     CALL dbcsr_data_clear(data_block, ub2=f_col_f)
                  ENDIF
                  IF (do_part_crop_l_col .AND. stored_col .EQ. block_col_bounds(2)) THEN
                     CALL dbcsr_data_clear(data_block, lb2=l_col_l)
                  ENDIF
               ENDIF
            ENDIF
            !
            ! Set meta data (global or local indexing)
            IF (dst_proc .EQ. mynode) THEN
               stored_row = local_g2l_map_rows(stored_row)
               stored_col = local_g2l_map_cols(stored_col)
               ! Count the maximum possible multiplies per row for on-the-fly filtering
               IF (is_left .AND. otf_filtering) THEN
                  left_total_row_counts(stored_row) = &
                     left_total_row_counts(stored_row)+1
               ENDIF
            ENDIF
            meta_buffer_p(send_displs(imeta, myt, ui, dst_proc)) = stored_row
            meta_buffer_p(send_displs(imeta, myt, ui, dst_proc)+1) = stored_col
            meta_buffer_p(send_displs(imeta, myt, ui, dst_proc)+2) = offset_data(myt, ui, dst_proc)
            !
            send_displs(imeta, myt, ui, dst_proc) = send_displs(imeta, myt, ui, dst_proc)+3
            send_displs(idata, myt, ui, dst_proc) = send_displs(idata, myt, ui, dst_proc)+nze
            offset_data(myt, ui, dst_proc) = offset_data(myt, ui, dst_proc)+nze
!$          IF (is_left .AND. do_symmetry) THEN
!$             call omp_unset_lock(locks(myt))
!$          ENDIF
         ENDDO
      ENDDO
      CALL dbcsr_iterator_stop(iter)
      CALL timestop(handle2)
      !
      IF (do_part_crop_row .OR. do_part_crop_col) THEN
         CALL dbcsr_data_clear_pointer(data_block)
         CALL dbcsr_data_release(data_block)
      ENDIF
      !
!$OMP BARRIER
!$OMP MASTER
      !
      ! Allocate data/meta to recv
      IF (dbcsr_data_valid(make_buffers_data_recv)) THEN
         IF (dbcsr_data_get_type(make_buffers_data_recv) .NE. data_type) THEN
            CALL dbcsr_data_release(make_buffers_data_recv)
         ENDIF
      ENDIF
      IF (dbcsr_data_valid(make_buffers_data_recv)) THEN
         CALL dbcsr_data_ensure_size(make_buffers_data_recv, my_size_recv(idata)-1, nocopy=.TRUE.)
      ELSE
         CALL dbcsr_data_init(make_buffers_data_recv)
         CALL dbcsr_data_new(make_buffers_data_recv, data_type, my_size_recv(idata)-1, &
                             memory_type=memtype_mpi_buffer)
      ENDIF
      CALL ensure_array_size(make_buffers_meta_recv, ub=my_size_recv(imeta)-1, &
                             nocopy=.TRUE., memory_type=memtype_mpi_buffer)
      ! Exchange data
      CALL timeset(routineN//"_data", handle2)
      CALL hybrid_alltoall_any(make_buffers_data_send, send_size_proc(idata, :), send_displ_proc(idata, :), &
                               make_buffers_data_recv, recv_size_proc(idata, :), recv_displ_proc(idata, :), &
                               mp_obj, &
                               most_ptp=.TRUE., remainder_ptp=.TRUE., no_hybrid=.FALSE.)
      CALL hybrid_alltoall_i1(make_buffers_meta_send, send_size_proc(imeta, :), send_displ_proc(imeta, :), &
                              make_buffers_meta_recv, recv_size_proc(imeta, :), recv_displ_proc(imeta, :), &
                              mp_obj, &
                              most_ptp=.TRUE., remainder_ptp=.TRUE., no_hybrid=.FALSE.)
      CALL timestop(handle2)
!$OMP END MASTER
!$OMP BARRIER
!$    IF (is_left .AND. do_symmetry) THEN
!$       call omp_destroy_lock(locks(ithread))
!$    ENDIF
      !
      ! Arrange data in the local buffers in images
      data_buffer_p => buffer%data
      meta_buffer_p => buffer%meta
      DO ui = 1, nimages
         ! Check for empty images
         IF (local_images_size(imeta, ui) .EQ. 0) CYCLE
         DO dst_proc = 0, nprocs_total-1
            IF (recv_sizes(imeta, ithread, ui, dst_proc) .EQ. 0) CYCLE
            ! Skip local data
            IF (dst_proc .EQ. mynode) THEN
               offset_threads(:, ithread, ui) = offset_threads(:, ithread, ui)+ &
                                                recv_sizes(:, ithread, ui, dst_proc)
            ELSE
               ! Copy meta, block by block
               DO blk = recv_displs(imeta, ithread, ui, dst_proc), &
                        recv_displs(imeta, ithread, ui, dst_proc)+recv_sizes(imeta, ithread, ui, dst_proc)-1, 3
                  stored_row = local_g2l_map_rows(make_buffers_meta_recv(blk))
                  stored_col = local_g2l_map_cols(make_buffers_meta_recv(blk+1))
                  meta_buffer_p(offset_threads(imeta, ithread, ui)) = stored_row
                  meta_buffer_p(offset_threads(imeta, ithread, ui)+1) = stored_col
                  meta_buffer_p(offset_threads(imeta, ithread, ui)+2) = make_buffers_meta_recv(blk+2)+ &
                                                                        offset_threads(idata, ithread, ui)
                  offset_threads(imeta, ithread, ui) = offset_threads(imeta, ithread, ui)+3
                  ! Count the maximum possible multiplies per row for on-the-fly filtering
                  IF (is_left .AND. otf_filtering) THEN
!$OMP ATOMIC
                     left_total_row_counts(stored_row) = &
                        left_total_row_counts(stored_row)+1
                  ENDIF
               ENDDO
               ! Copy data
               CALL dbcsr_data_set(data_buffer_p, &
                                   offset_threads(idata, ithread, ui)+local_images_displ(idata, ui)+1, &
                                   recv_sizes(idata, ithread, ui, dst_proc), &
                                   make_buffers_data_recv, recv_displs(idata, ithread, ui, dst_proc))
               offset_threads(idata, ithread, ui) = offset_threads(idata, ithread, ui)+ &
                                                    recv_sizes(idata, ithread, ui, dst_proc)
            ENDIF
         ENDDO
      ENDDO
!$OMP END PARALLEL
      DEALLOCATE (send_sizes, recv_sizes)
      DEALLOCATE (send_displs, recv_displs, offset_data, offset_threads)
      DEALLOCATE (send_size_proc, send_displ_proc, recv_size_proc, recv_displ_proc)
      !
      IF (is_left .AND. otf_filtering) THEN
         CALL mp_isum(left_total_row_counts, dbcsr_mp_my_row_group(mp_obj), request_count_rows)
      ENDIF
      !
      CALL setup_rec_index_images(buffer%meta, img_nblks_rows, img_nblks_cols, &
                                  local_images_size(imeta, :), local_images_displ(imeta, :), &
                                  size_index, is_left)
      IF (buffer%has_rma_win) THEN
         do_win_create(1) = do_win_create(1) .OR. dbcsr_data_exists(buffer%data_before_resize)
         do_win_create(2) = do_win_create(2) .OR. ASSOCIATED(buffer%meta_before_resize)
         CALL mp_isum(do_win_create, buffer%subgrp, requests_win_create(irequests))
      ENDIF
      !
      IF (is_left) THEN
         NULLIFY (local_g2l_map_rows)
         DEALLOCATE (local_g2l_map_cols)
      ELSE
         DEALLOCATE (local_g2l_map_rows)
         NULLIFY (local_g2l_map_cols)
      ENDIF
!$    IF (is_left .AND. do_symmetry) THEN
!$       DEALLOCATE (locks)
!$    ENDIF
      !
      DEALLOCATE (img_nblks_rows, img_nblks_cols)
      DEALLOCATE (local_images_displ)
      !
      CALL timestop(handle)
   END SUBROUTINE make_buffers

! **************************************************************************************************
!> \brief Make communicators for A and B matrices
!> \param my_num_layers_3D ...
!> \param side3D ...
!> \param mp_obj ...
!> \param is_left ...
!> \param buffer ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE make_layers_3D_AB(my_num_layers_3D, side3D, mp_obj, is_left, buffer)
      INTEGER, INTENT(IN)                                :: my_num_layers_3D, side3D
      TYPE(dbcsr_mp_obj), INTENT(IN)                     :: mp_obj
      LOGICAL, INTENT(IN)                                :: is_left
      TYPE(dbcsr_buffer), INTENT(INOUT)                  :: buffer

      INTEGER                                            :: color, key, mygrp, mypcol, myprow

      ! Switch to single layer communicator
      IF (my_num_layers_3D .LE. 1) THEN
         IF (buffer%num_layers_3D .GT. 1 .AND. buffer%subgrp .NE. mp_comm_null) &
            CALL mp_comm_free(buffer%subgrp)
         buffer%num_layers_3D = 1
         IF (is_left) THEN
            buffer%subgrp = dbcsr_mp_my_row_group(mp_obj)
         ELSE
            buffer%subgrp = dbcsr_mp_my_col_group(mp_obj)
         ENDIF
         RETURN
      ENDIF
      !
      ! Check if any existing 3D communicator can be reused
      mygrp = dbcsr_mp_group(mp_obj)
      IF (buffer%grp .EQ. mygrp .AND. buffer%num_layers_3D .EQ. my_num_layers_3D) RETURN
      !
      ! Reset previous 3D communicator
      IF (buffer%num_layers_3D .GT. 1 .AND. buffer%subgrp .NE. mp_comm_null) &
         CALL mp_comm_free(buffer%subgrp)
      !
      myprow = dbcsr_mp_myprow(mp_obj)
      mypcol = dbcsr_mp_mypcol(mp_obj)
      IF (is_left) THEN
         color = MOD(myprow, side3D)
         ! Column-major order
         key = mypcol*(dbcsr_mp_nprows(mp_obj)/side3D)+myprow/side3D
      ELSE
         color = MOD(mypcol, side3D)
         ! Row-major order
         key = myprow*(dbcsr_mp_npcols(mp_obj)/side3D)+mypcol/side3D
      ENDIF
      CALL mp_comm_split_direct(mygrp, buffer%subgrp, color, key)
      buffer%num_layers_3D = my_num_layers_3D
   END SUBROUTINE make_layers_3D_AB

! **************************************************************************************************
!> \brief Return the rank of the 3D layer (3D communicator for C), Column-major order
!> \param myprow ...
!> \param mypcol ...
!> \param nprows ...
!> \param side3D ...
!> \retval get_rank3D ...
! **************************************************************************************************
   PURE FUNCTION get_rank3D(myprow, mypcol, nprows, side3D)
      INTEGER, INTENT(IN)                                :: myprow, mypcol, nprows, side3D
      INTEGER                                            :: get_rank3D

      get_rank3D = myprow/side3D+(nprows/side3D)*(mypcol/side3D)
   END FUNCTION get_rank3D

! **************************************************************************************************
!> \brief Make communicators for 3D layers for C-reduction
!> \param my_num_layers_3D ...
!> \param mp_obj ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE make_layers_3D_C_reduction(my_num_layers_3D, mp_obj)
      INTEGER, INTENT(IN)                                :: my_num_layers_3D
      TYPE(dbcsr_mp_obj), INTENT(INOUT)                  :: mp_obj

      CHARACTER(len=100)                                 :: msg
      INTEGER                                            :: color, key, mygrp, mypcol, myprow, &
                                                            npcols, nprows, numnodes
      LOGICAL                                            :: do_layers_3D
      LOGICAL, SAVE                                      :: warning = .TRUE.

      CALL dbcsr_mp_grid_setup(mp_obj)
      IF (my_num_layers_3D .LE. 1) THEN
         ! Reset 3D communicator if it was previously declared
         IF (layers_3D_C_reduction%num_layers_3D .GT. 1) CALL release_layers_3D_C_reduction()
         RETURN
      ENDIF
      !
      ! Check if any existing 3D communicator can be reused
      mygrp = dbcsr_mp_group(mp_obj)
      IF (layers_3D_C_reduction%grp .EQ. mygrp .AND. &
          layers_3D_C_reduction%num_layers_3D .EQ. my_num_layers_3D) RETURN
      !
      ! Reset 3D communicator
      CALL release_layers_3D_C_reduction()
      !
      ! Checks for 3D algorithm
      numnodes = dbcsr_mp_numnodes(mp_obj)
      nprows = dbcsr_mp_nprows(mp_obj)
      npcols = dbcsr_mp_npcols(mp_obj)
      IF (dbcsr_cfg%use_mpi_rma) THEN
         IF (nprows .NE. npcols) THEN
            ! No square topology, scale the maximum coordinate
            do_layers_3D = MAX(nprows, npcols) .EQ. (my_num_layers_3D*MIN(nprows, npcols)) .AND. &
                           my_num_layers_3D .LE. MIN(nprows, npcols)
         ELSE
            ! Square topology, scale both coordinates
            do_layers_3D = ((nprows/NINT(SQRT(REAL(MAX(1, my_num_layers_3D), KIND=real_8))))**2)* &
                           my_num_layers_3D .EQ. (nprows*npcols)
         ENDIF
         IF (.NOT. do_layers_3D .AND. warning) THEN
            WRITE (UNIT=msg, FMT='(A,I3,A,I3,A,I3,A)') "Cannot make 3D layers with ", my_num_layers_3D, &
               " layers and (", nprows, "x", npcols, ") ranks! Run with a single layer."
            CPWARN(msg)
            warning = .FALSE.
         ENDIF
         IF (do_layers_3D) THEN
            layers_3D_C_reduction%grp = mygrp
            layers_3D_C_reduction%num_layers_3D = my_num_layers_3D
            layers_3D_C_reduction%max_num_layers_3D = &
               MAX(layers_3D_C_reduction%max_num_layers_3D, &
                   my_num_layers_3D)
            layers_3D_C_reduction%side3D = NINT(SQRT(REAL(numnodes/my_num_layers_3D, KIND=real_8)))
            !
            ! Create a new 3D communicator
            myprow = dbcsr_mp_myprow(mp_obj)
            mypcol = dbcsr_mp_mypcol(mp_obj)
            ! Row-wise order for color
            color = MOD(myprow, layers_3D_C_reduction%side3D)* &
                    layers_3D_C_reduction%side3D+MOD(mypcol, layers_3D_C_reduction%side3D)
            ! Column-major order
            key = get_rank3D(myprow, mypcol, nprows, layers_3D_C_reduction%side3D)
            CALL mp_comm_split_direct(mygrp, layers_3D_C_reduction%grp3D, color, key)
            !
            ! Create a 3D-row communicator based on the 3D communicator
            color = key/(nprows/layers_3D_C_reduction%side3D)
            CALL mp_comm_split_direct(layers_3D_C_reduction%grp3D, &
                                      layers_3D_C_reduction%rowgrp3D, color, key)
         ENDIF
      ELSE
         CPWARN('Cannot make 3D layers without experimental MPI algorithm enabled!')
      ENDIF
   END SUBROUTINE make_layers_3D_C_reduction

! **************************************************************************************************
!> \brief Release communicators for 3D layers for C-reduction
!> \param release_buffers ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE release_layers_3D_C_reduction(release_buffers)
      LOGICAL, OPTIONAL                                  :: release_buffers

      INTEGER                                            :: ibuff

      layers_3D_C_reduction%grp = mp_comm_null
      IF (layers_3D_C_reduction%rowgrp3D .NE. mp_comm_null) CALL mp_comm_free(layers_3D_C_reduction%rowgrp3D)
      IF (layers_3D_C_reduction%grp3D .NE. mp_comm_null) CALL mp_comm_free(layers_3D_C_reduction%grp3D)
      layers_3D_C_reduction%rowgrp3D = mp_comm_null
      layers_3D_C_reduction%grp3D = mp_comm_null
      layers_3D_C_reduction%num_layers_3D = 1
      layers_3D_C_reduction%side3D = HUGE(1)

      IF (PRESENT(release_buffers)) THEN
         IF (release_buffers .AND. ALLOCATED(layers_3D_C_reduction%data_red3D)) THEN
            DO ibuff = 1, SIZE(layers_3D_C_reduction%data_red3D)
               CALL dbcsr_data_release(layers_3D_C_reduction%data_red3D(ibuff))
            ENDDO
            DEALLOCATE (layers_3D_C_reduction%data_red3D)
         ENDIF
      ENDIF
   END SUBROUTINE release_layers_3D_C_reduction

! **************************************************************************************************
!> \brief Multiplies two DBCSR matrices (experimental MPI algorithm).
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
   SUBROUTINE multiply_3D(imgdist_left, imgdist_right, &
                          matrix_left, matrix_right, &
                          product_matrix, &
                          retain_sparsity, &
                          filter_eps, flop, keep_product_data)
      TYPE(dbcsr_imagedistribution_obj), INTENT(INOUT)   :: imgdist_left, imgdist_right
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_left, matrix_right
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: product_matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: retain_sparsity
      REAL(kind=real_8), INTENT(IN), OPTIONAL            :: filter_eps
      INTEGER(KIND=int_8), INTENT(OUT)                   :: flop
      LOGICAL, INTENT(IN)                                :: keep_product_data

      CHARACTER(len=*), PARAMETER :: routineN = 'multiply_3D', routineP = moduleN//':'//routineN

      INTEGER :: blk, data_type, data_type_byte, final_step_k, grp_left, grp_right, handle, &
         handle1, handle2, icol3D, icol3D_send, ileft_buffer_calc, ileft_buffer_comm, &
         iright_buffer_calc, iright_buffer_comm, irow3D, irow3D_send, istep_k_ordered, ithread, &
         ivirt_k, last_step_k, left_col_mult, left_col_nimages, left_col_total_nimages, &
         left_max_data_size, left_max_meta_size, left_myfirstvcol, left_myfirstvrow, left_mypcol, &
         left_myprow, left_npcols, left_nprows, left_row_mult, left_row_nimages, &
         leftovers_first_k, leftovers_k, leftovers_shift_k, leftovers_start_k, min_nimages, &
         mycol3D, mypcol, myprow
      INTEGER :: myrank3D, myrow3D, myt, nblkrows_local, nbuffers, nbuffers_norms, ncols3D, &
         nranks3D, nrows3D, nthreads, numnodes, nvirt_k, proc3D_recv, proc3D_send, recv_vcol, &
         recv_vrow, request_epss, request_keep_sparsity, right_col_mult, right_col_nimages, &
         right_max_data_size, right_max_meta_size, right_myfirstvcol, right_myfirstvrow, &
         right_mypcol, right_myprow, right_npcols, right_nprows, right_row_mult, &
         right_row_nimages, right_row_total_nimages, row, shift3D, shift3D_comm, shift3D_recv, &
         size_guess, size_guess_init, start_k, start_k_ordered, v_ki
      INTEGER(KIND=int_8)                                :: mem
      INTEGER, ALLOCATABLE, DIMENSION(:) :: left_vrow, product_matrix_epss_displ, &
         product_matrix_epss_size, product_matrix_meta, product_matrix_size_recv, &
         product_matrix_size_send, requests_reduction, right_vcol
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: product_matrix_meta_displ, &
                                                            product_matrix_meta_size
      INTEGER, DIMENSION(2)                              :: requests_reduction_size
      INTEGER, DIMENSION(:), POINTER :: col_blk_sizes2enum, enum2col_blk_sizes, &
         enum2row_blk_sizes, k_sizes, product_matrix_meta_recv, product_matrix_meta_send, &
         row_blk_sizes2enum
      INTEGER, DIMENSION(:, :, :), POINTER               :: left_displ_layers3D, &
                                                            left_images_size_layers3D, &
                                                            right_displ_layers3D, &
                                                            right_images_size_layers3D
      INTEGER, DIMENSION(dbcsr_slot_nblkrows_total:&
         dbcsr_slot_nfullcols_local)                     :: left_global_indices, right_global_indices
      INTEGER, POINTER                                   :: istep_k_comm
      INTEGER, TARGET                                    :: istep_k, istep_k_comm_curr
      LOGICAL                                            :: do_layers3D, do_square_layers3D, &
                                                            first_k, first_v_k, is_not_comm, &
                                                            keep_sparsity, otf_filtering
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: do_comm_left, do_comm_right
      REAL(kind=sp)                                      :: filter_eps_sp
      REAL(kind=sp), ALLOCATABLE, DIMENSION(:), TARGET   :: row_max_epss
      REAL(kind=sp), ALLOCATABLE, DIMENSION(:, :)        :: left_norms, right_norms
      REAL(kind=sp), DIMENSION(:), POINTER               :: product_matrix_epss
      TYPE(dbcsr_2d_array_obj)                           :: product_matrix3D
      TYPE(dbcsr_buffer), ALLOCATABLE, DIMENSION(:), &
         TARGET                                          :: left_buffers, right_buffers
      TYPE(dbcsr_buffer), POINTER                        :: left_buffer_p, right_buffer_p
      TYPE(dbcsr_data_obj)                               :: data_get, data_send
      TYPE(dbcsr_mm_multrec_type_p), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: multrec
      TYPE(dbcsr_mp_obj)                                 :: left_mp_obj, product_mp_obj, right_mp_obj
      TYPE(mn_local_sizes), ALLOCATABLE, DIMENSION(:)    :: m_sizes, n_sizes

      CALL timeset(routineN, handle)
      !
      NULLIFY (row_blk_sizes2enum, enum2row_blk_sizes)
      NULLIFY (col_blk_sizes2enum, enum2col_blk_sizes)
      NULLIFY (k_sizes)
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
      ! Dummy checks
      IF (.NOT. ASSOCIATED(product_matrix%wms)) &
         CPABORT("Work matrices do not exist")
      IF (SIZE(product_matrix%wms) .NE. nthreads) &
         CPABORT("Work matrices not correctly sized.")
      IF (.NOT. buffers_win%left%is_valid .OR. &
          .NOT. buffers_win%right%is_valid .OR. &
          .NOT. ASSOCIATED(buffers_win%left%meta) .OR. &
          .NOT. ASSOCIATED(buffers_win%right%meta) .OR. &
          .NOT. ASSOCIATED(left_images_size) .OR. &
          .NOT. ASSOCIATED(right_images_size) .OR. &
          .NOT. ALLOCATED(left_local_images_size) .OR. &
          .NOT. ALLOCATED(right_local_images_size)) &
         CPABORT("No buffers associated for the experimental algo!")
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
      product_mp_obj = dbcsr_distribution_mp(product_matrix%dist)
      numnodes = dbcsr_mp_numnodes(product_mp_obj)
      myprow = dbcsr_mp_myprow(product_mp_obj)
      mypcol = dbcsr_mp_mypcol(product_mp_obj)
      left_nprows = dbcsr_mp_nprows(left_mp_obj)
      left_npcols = dbcsr_mp_npcols(left_mp_obj)
      left_myprow = dbcsr_mp_myprow(left_mp_obj)
      left_mypcol = dbcsr_mp_mypcol(left_mp_obj)
      left_myfirstvrow = MOD(left_myprow, layers_3D_C_reduction%side3D)*left_row_nimages
      left_myfirstvcol = MOD(left_mypcol, layers_3D_C_reduction%side3D)*left_col_nimages
      right_nprows = dbcsr_mp_nprows(right_mp_obj)
      right_npcols = dbcsr_mp_npcols(right_mp_obj)
      right_myprow = dbcsr_mp_myprow(right_mp_obj)
      right_mypcol = dbcsr_mp_mypcol(right_mp_obj)
      right_myfirstvrow = MOD(right_myprow, layers_3D_C_reduction%side3D)*right_row_nimages
      right_myfirstvcol = MOD(right_mypcol, layers_3D_C_reduction%side3D)*right_col_nimages
      left_col_total_nimages = left_npcols*left_col_nimages
      right_row_total_nimages = right_nprows*right_row_nimages
      grp_right = buffers_win%right%subgrp
      grp_left = buffers_win%left%subgrp
      !
      do_layers3D = layers_3D_C_reduction%num_layers_3D .GT. 1
      myrow3D = myprow/layers_3D_C_reduction%side3D+1
      mycol3D = mypcol/layers_3D_C_reduction%side3D+1
      nrows3D = SIZE(left_images_size, 3)
      ncols3D = SIZE(right_images_size, 3)
      myrank3D = get_rank3D(myprow, mypcol, dbcsr_mp_nprows(product_mp_obj), layers_3D_C_reduction%side3D)
      nranks3D = layers_3D_C_reduction%num_layers_3D
      myprow = MOD(myprow, layers_3D_C_reduction%side3D)
      mypcol = MOD(mypcol, layers_3D_C_reduction%side3D)
      !
      ! Dummy checks
      ! subcommunicators
      IF (.NOT. dbcsr_mp_has_subgroups(right_mp_obj)) &
         CPABORT("Experimental algorithm requires rows subcommunicators for right matrix!")
      IF (.NOT. dbcsr_mp_has_subgroups(left_mp_obj)) &
         CPABORT("Experimental algorithm requires columns subcommunicators for left matrix!")
      ! Right col nimages
      IF (right_col_nimages .NE. 1) &
         CPABORT("Col nimages for right matrix is not 1!")
      ! Left row nimages
      IF (left_row_nimages .NE. 1) &
         CPABORT("Row nimages for left matrix is not 1!")
      ! left/right matching
      IF (left_col_nimages .NE. right_row_mult) &
         CPABORT("Left/Right image mismatch")
      IF (left_col_mult .NE. right_row_nimages) &
         CPABORT("Left/Right image mismatch")
      IF (left_col_nimages*left_npcols .NE. right_row_nimages*right_nprows) &
         CPABORT("Left/Right total mismatch")
      ! product/left matching
      IF (left_row_mult*dbcsr_mp_nprows(product_mp_obj) .NE. left_nprows) &
         CPABORT("Product/Left total mismatch")
      ! product/left matching
      IF (right_col_mult*dbcsr_mp_npcols(product_mp_obj) .NE. right_npcols) &
         CPABORT("Product/Right total mismatch")
      ! Check sizes from make_buffers
      IF (SIZE(left_images_size, 2) .NE. left_col_nimages .OR. &
          SIZE(right_images_size, 2) .NE. right_row_nimages) &
         CPABORT("Mismatch in the sizes")
      !
      dbcsr_mpi_statistics%nimages = MAX(dbcsr_mpi_statistics%nimages, left_col_nimages)
      dbcsr_mpi_statistics%nimages = MAX(dbcsr_mpi_statistics%nimages, right_row_nimages)
      !
      ! The main transfer loop goes through the virtual rows/columns.
      ! The number of steps may be smaller if the grid dimension is very
      ! non-optimal (both left column images and right row images are >
      ! 1).
      min_nimages = MIN(left_col_nimages, right_row_nimages)
      nvirt_k = left_npcols*left_col_nimages
      !
      ! Check RMA windows creation for original data
      CALL win_setup(buffers_win%left, do_win_create_left, requests_win_create(1))
      CALL win_setup(buffers_win%right, do_win_create_right, requests_win_create(2))
      !
      ! Count the maximum possible multiplies per row for on-the-fly filtering
      ALLOCATE (product_matrix_epss_size(nrows3D), product_matrix_epss_displ(nrows3D))
      IF (otf_filtering) THEN
         ! Wait for counts (sent in make_buffers)
         CALL timeset(routineN//"_count_rows", handle1)
         CALL mp_wait(request_count_rows)
         !
         nblkrows_local = SIZE(left_total_row_counts)
         ALLOCATE (row_max_epss(nblkrows_local))
         filter_eps_sp = REAL(filter_eps, KIND=KIND(row_max_epss))
!$OMP PARALLEL DO DEFAULT (NONE) &
!$OMP SHARED(nblkrows_local,row_max_epss,filter_eps_sp,&
!$OMP        left_total_row_counts)
         ! Determine the maximum per-block epsilon
         DO row = 1, nblkrows_local
            row_max_epss(row) = &
               filter_eps_sp/REAL(MAX(1, left_total_row_counts(row)), KIND=KIND(row_max_epss))
         ENDDO
!$OMP END PARALLEL DO
         DEALLOCATE (left_total_row_counts)
         !
         IF (do_layers3D .AND. nrows3D .GT. 1) THEN
            CALL mp_allgather(SIZE(row_max_epss), &
                              product_matrix_epss_size, &
                              layers_3D_C_reduction%rowgrp3D)
            size_guess = 0
            DO irow3D = 1, nrows3D
               product_matrix_epss_displ(irow3D) = size_guess
               size_guess = size_guess+product_matrix_epss_size(irow3D)
            ENDDO
            ALLOCATE (product_matrix_epss(size_guess))
            CALL mp_iallgather(row_max_epss, &
                               product_matrix_epss, product_matrix_epss_size, product_matrix_epss_displ, &
                               layers_3D_C_reduction%rowgrp3D, request_epss)
         ELSE
            product_matrix_epss_size(nrows3D) = SIZE(row_max_epss)
            product_matrix_epss_displ(nrows3D) = 0
            product_matrix_epss => row_max_epss
         ENDIF
         CALL timestop(handle1)
      ELSE
         product_matrix_epss_size(:) = 0
         product_matrix_epss_displ(:) = 0
         ALLOCATE (product_matrix_epss(0))
      ENDIF
      !
      ! Exchange 3D meta for C matrix
      IF (do_layers3D .AND. keep_sparsity) THEN
         ALLOCATE (product_matrix_meta_size(nrows3D, ncols3D))
         CALL mp_allgather(product_matrix%index(dbcsr_slot_size), &
                           product_matrix_meta_size, layers_3D_C_reduction%grp3D)
         ALLOCATE (product_matrix_meta_displ(nrows3D, ncols3D))
         size_guess = 0
         DO icol3D = 1, ncols3D
            DO irow3D = 1, nrows3D
               product_matrix_meta_displ(irow3D, icol3D) = size_guess
               size_guess = size_guess+product_matrix_meta_size(irow3D, icol3D)
            ENDDO
         ENDDO
         ALLOCATE (product_matrix_meta(size_guess))
         product_matrix%index(dbcsr_slot_nblks) = product_matrix%nblks
         product_matrix%index(dbcsr_slot_nze) = product_matrix%nze
         CALL mp_iallgather(product_matrix%index(1:product_matrix%index(dbcsr_slot_size)), &
                            product_matrix_meta, product_matrix_meta_size, product_matrix_meta_displ, &
                            layers_3D_C_reduction%grp3D, request_keep_sparsity)
      ENDIF
      !
      ! Wait refs and max norms (sent in make_buffers)
      CALL timeset(routineN//"_sizes", handle1)
      CALL mp_waitall(requests)
      CALL timestop(handle1)
      DEALLOCATE (right_local_images_size, left_local_images_size)
      !
      ! Needs to remap refs for virtual coordinates 3D
      CALL remap_layers3D(left_images_size, left_images_size_layers3D, left_displ_layers3D, &
                          left_max_data_size, left_max_meta_size)
      CALL remap_layers3D(right_images_size, right_images_size_layers3D, right_displ_layers3D, &
                          right_max_data_size, right_max_meta_size)
      left_max_meta_size = left_max_meta_size+dbcsr_num_slots
      right_max_meta_size = right_max_meta_size+dbcsr_num_slots
      !
      do_square_layers3D = .FALSE.
      nbuffers_norms = 1
      IF (nvirt_k .EQ. 1) THEN
         nbuffers = 1
      ELSEIF (nrows3D .NE. ncols3D .OR. nranks3D .EQ. 1) THEN
         nbuffers = 2
      ELSE
         ! Note that nrows3D==ncols3D >= 2
         ! Last buffer is used as temporary for communications
         nbuffers = nrows3D+1
         nbuffers_norms = nrows3D
         do_square_layers3D = .TRUE.
      ENDIF
      !
      ! update capacity of memory-pools
      IF (ASSOCIATED(memtype_abpanel_1%pool)) &
         CALL dbcsr_mempool_limit_capacity(memtype_abpanel_1%pool, &
                                           capacity=2)
      IF (ASSOCIATED(memtype_abpanel_2%pool)) &
         CALL dbcsr_mempool_limit_capacity(memtype_abpanel_2%pool, &
                                           capacity=2)
      IF (has_acc) THEN
         ! enumerate the blocksizes to keep the following 2D-arrays small.
         CALL enumerate_blk_sizes(dbcsr_row_block_sizes(matrix_right), &
                                  dbcsr_max_row_size(matrix_right), &
                                  row_blk_sizes2enum, enum2row_blk_sizes)
         CALL enumerate_blk_sizes(dbcsr_col_block_sizes(matrix_right), &
                                  dbcsr_max_col_size(matrix_right), &
                                  col_blk_sizes2enum, enum2col_blk_sizes)
      ENDIF
      IF (nranks3D .GT. 1) THEN
         CALL dbcsr_mempool_limit_capacity(memtype_mpi_product%pool, &
                                           capacity=nranks3D-1)
      ENDIF
      !
      ! Prepare buffers for computation
      IF (nvirt_k .GT. 1) THEN
         ! Right
         CALL buffer_init(buffers_2%right, data_type, &
                          right_max_data_size, &
                          right_max_meta_size, &
                          num_data=(nbuffers/2), &
                          data_memory_type=memtype_abpanel_2, &
                          trs_memory_type=memtype_trsbuffer_2)
         ! Left
         CALL buffer_init(buffers_2%left, data_type, &
                          left_max_data_size, &
                          left_max_meta_size, &
                          num_data=(nbuffers/2), &
                          data_memory_type=memtype_abpanel_2)
      ENDIF
      !
      ! Prepare buffers for communication
      ! Right
      CALL buffer_init(buffers_1%right, data_type, &
                       right_max_data_size, &
                       right_max_meta_size, &
                       num_data=(nbuffers-nbuffers/2), &
                       data_memory_type=memtype_abpanel_1, &
                       trs_memory_type=memtype_trsbuffer_1)
      ! Left
      CALL buffer_init(buffers_1%left, data_type, &
                       left_max_data_size, &
                       left_max_meta_size, &
                       num_data=(nbuffers-nbuffers/2), &
                       data_memory_type=memtype_abpanel_1)
      !
      CALL setup_buffers(buffers_1%right, buffers_2%right, &
                         right_buffers, nbuffers, &
                         right_max_data_size, &
                         right_max_meta_size, &
                         matrix_right, imgdist_right)
      CALL setup_buffers(buffers_1%left, buffers_2%left, &
                         left_buffers, nbuffers, &
                         left_max_data_size, &
                         left_max_meta_size, &
                         matrix_left, imgdist_left)
      !
      ! Setup the receive data pointers
      CALL dbcsr_data_init(data_get)
      CALL dbcsr_data_new(data_get, data_type)
      IF (do_layers3D) THEN
         CALL dbcsr_data_init(data_send)
         CALL dbcsr_data_new(data_send, data_type)
         ! Prepare buffers for 3D reduction
         IF (ALLOCATED(layers_3D_C_reduction%data_red3D)) THEN
            IF (SIZE(layers_3D_C_reduction%data_red3D) .LT. nthreads .OR. &
                layers_3D_C_reduction%data_type .NE. data_type) THEN
               DO myt = 1, SIZE(layers_3D_C_reduction%data_red3D)
                  CALL dbcsr_data_release(layers_3D_C_reduction%data_red3D(myt))
               ENDDO
               DEALLOCATE (layers_3D_C_reduction%data_red3D)
               layers_3D_C_reduction%data_type = 0
            ENDIF
         ENDIF
         IF (.NOT. ALLOCATED(layers_3D_C_reduction%data_red3D)) THEN
            ALLOCATE (layers_3D_C_reduction%data_red3D(nthreads))
            DO myt = 1, nthreads
               CALL dbcsr_data_init(layers_3D_C_reduction%data_red3D(myt))
               CALL dbcsr_data_new(layers_3D_C_reduction%data_red3D(myt), data_type)
            ENDDO
            layers_3D_C_reduction%data_type = data_type
         ENDIF
         ALLOCATE (product_matrix_size_send(nthreads+1), product_matrix_size_recv(nthreads+1))
         ALLOCATE (requests_reduction((nthreads+1)*2))
      ENDIF
      !
      ! These values for meta data are used for global values
      right_global_indices(dbcsr_slot_nblkrows_total:dbcsr_slot_nfullcols_local) = &
         (/ &
         dbcsr_nblkrows_total(matrix_right), &
         dbcsr_nblkcols_total(matrix_right), &
         dbcsr_nfullrows_total(matrix_right), &
         dbcsr_nfullcols_total(matrix_right), &
         0, 0, &
         dbcsr_nfullrows_local(matrix_right), &
         dbcsr_nfullcols_local(matrix_right)/)
      left_global_indices(dbcsr_slot_nblkrows_total:dbcsr_slot_nfullcols_local) = &
         (/ &
         dbcsr_nblkrows_total(matrix_left), &
         dbcsr_nblkcols_total(matrix_left), &
         dbcsr_nfullrows_total(matrix_left), &
         dbcsr_nfullcols_total(matrix_left), &
         0, 0, &
         dbcsr_nfullrows_local(matrix_left), &
         dbcsr_nfullcols_local(matrix_left)/)
      !
      ! Evaluate sizes for workspaces
      size_guess_init = 1
      IF (.NOT. keep_sparsity .AND. has_acc) THEN
         size_guess_init = product_matrix_size_guess(matrix_left, matrix_right, product_matrix, &
                                                     left_max_data_size, right_max_data_size, &
                                                     left_col_nimages, right_row_nimages, &
                                                     nthreads)
      ENDIF
      !
      ! Preallocate norms arrays
      IF (otf_filtering) THEN
         ALLOCATE (right_norms(right_max_meta_size/3, nbuffers_norms))
         ALLOCATE (left_norms(left_max_meta_size/3, nbuffers_norms))
         IF (do_layers3D .AND. nrows3D .GT. 1) THEN
            CALL mp_wait(request_epss)
            DEALLOCATE (row_max_epss)
         ENDIF
      ELSE
         ! The array must be valid when passed to called subroutines.
         ALLOCATE (right_norms(0, nbuffers_norms))
         ALLOCATE (left_norms(0, nbuffers_norms))
      ENDIF
      !
      IF (do_layers3D .AND. keep_sparsity) CALL mp_wait(request_keep_sparsity)
      !
      ALLOCATE (product_matrix3D%mats(nrows3D, ncols3D))
      DO icol3D = 1, ncols3D
         DO irow3D = 1, nrows3D
            NULLIFY (product_matrix3D%mats(irow3D, icol3D)%matrix)
         ENDDO
      ENDDO
      ALLOCATE (multrec(0:nthreads-1, nrows3D, ncols3D))
      !
      ! Here is the main loop
      ! 3D multiplication
      !
      CALL timeset(routineN//"_loop", handle1)
      ! Take into account when ticks are not multiple of 3D layers
      leftovers_k = MOD(nvirt_k, nranks3D)
      leftovers_first_k = leftovers_k*myrank3D
      leftovers_start_k = 0
      leftovers_shift_k = 0
      IF (leftovers_k .GT. 0) THEN
         ! This is only for nrows3D==ncols3D
         leftovers_start_k = (nvirt_k/nrows3D-1)*(myrank3D/nrows3D)- &
                             (leftovers_k/nrows3D-1)*(myrank3D/nrows3D)
         leftovers_shift_k = nranks3D*(leftovers_k/nrows3D)-leftovers_k*(MOD(myrank3D, nrows3D)+1)
      ENDIF
      ! Ticks bounds
      start_k = (nvirt_k/nranks3D)*myrank3D
      last_step_k = nvirt_k+leftovers_first_k
      final_step_k = last_step_k-nranks3D
      ! Shift layers to keep local layer as the last one in computation
      shift3D = (mycol3D-1)*nrows3D+ &
                (nrows3D-myrow3D+1)*(1-MOD(mycol3D, 2))+myrow3D*MOD(mycol3D, 2)
      iright_buffer_comm = 0
      ileft_buffer_comm = 0
      ALLOCATE (do_comm_right(ncols3D), do_comm_left(nrows3D))
      ALLOCATE (right_vcol(ncols3D), left_vrow(nrows3D))
      ALLOCATE (m_sizes(nrows3D), n_sizes(ncols3D))
      irow3D_send = 0
      icol3D_send = 0
      first_k = .TRUE.
      first_v_k = .TRUE.
      istep_k_comm_curr = leftovers_first_k
      istep_k_comm => istep_k_comm_curr
      grouped_steps_index: DO istep_k = leftovers_first_k, last_step_k
         !
         ! Wait data.  Exclude the first interation.
         wait: IF (istep_k .GT. leftovers_first_k) THEN
            IF (debug_mod) WRITE (*, '(1X,A)') routineN//" waiting for right and left"
            right_buffer_p => right_buffers(iright_buffer_calc)
            left_buffer_p => left_buffers(ileft_buffer_calc)
            IF (right_buffer_p%is_comm .AND. left_buffer_p%is_comm) THEN
               ! check if right matrix was already initialized
               IF (.NOT. right_buffer_p%matrix%valid) THEN
                  CALL timeset(routineN//"_comm_right", handle2)
                  CALL mp_waitall(right_buffer_p%get_requests(:))
                  CALL timestop(handle2)
               ENDIF
               ! check if left matrix was already initialized
               IF (.NOT. left_buffer_p%matrix%valid) THEN
                  CALL timeset(routineN//"_comm_left", handle2)
                  CALL mp_waitall(left_buffer_p%get_requests(:))
                  CALL timestop(handle2)
               ENDIF
            ENDIF
         ENDIF wait
         !
         ! Matrix transfer. Transfer in all but the last loop iteration.
         shift3D_comm = shift3D
         xfer: DO WHILE (istep_k_comm .LT. last_step_k)
            start_k_ordered = start_k
            istep_k_ordered = istep_k_comm
            ! Put leftovers ticks always first
            IF (leftovers_k .GT. 0) THEN
               IF (istep_k_comm .LT. leftovers_first_k+leftovers_k) THEN
                  start_k_ordered = leftovers_start_k
               ELSE
                  istep_k_ordered = istep_k_comm+leftovers_shift_k
               ENDIF
            ENDIF
            first_k = MOD(istep_k_ordered, nranks3D) .EQ. 0
            ivirt_k = istep_k_ordered/nranks3D
            IF (istep_k_comm .LT. leftovers_first_k+leftovers_k) THEN
               CALL row_col_3D_reflected(irow3D, icol3D, nrows3D, ncols3D, istep_k_ordered)
            ELSE
               CALL row_col_3D_reflected(irow3D, icol3D, nrows3D, ncols3D, shift3D)
               shift3D = shift3D+1
            ENDIF
            !
            v_ki = MOD(ivirt_k, min_nimages)
            ! Reset communication flags at the first layer
            IF ((first_k .OR. istep_k_comm .EQ. leftovers_first_k) .AND. &
                istep_k_comm .EQ. istep_k_comm_curr) THEN
               do_comm_right(:) = .TRUE.
               do_comm_left(:) = .TRUE.
            ENDIF
            ! Take first image global virtual coordinates
            IF (v_ki .EQ. 0) THEN
               IF (istep_k_comm .GE. leftovers_first_k+leftovers_k) first_v_k = .FALSE.
               start_k_ordered = start_k_ordered+ivirt_k
            ENDIF
            IF (v_ki .EQ. 0 .OR. (first_v_k .AND. min_nimages .GT. 1)) THEN
               CALL image_calculator(imgdist_right, &
                                     vprow=recv_vrow, &
                                     vpcol=right_vcol(icol3D), &
                                     mypcol=mypcol, &
                                     myvprow=right_myfirstvrow, &
                                     myvpcol=right_myfirstvcol+(icol3D-1)*layers_3D_C_reduction%side3D, &
                                     vprow_shift=start_k_ordered, &
                                     shifting='R')
               CALL image_calculator(imgdist_left, &
                                     vprow=left_vrow(irow3D), &
                                     vpcol=recv_vcol, &
                                     myprow=myprow, &
                                     myvprow=left_myfirstvrow+(irow3D-1)*layers_3D_C_reduction%side3D, &
                                     myvpcol=left_myfirstvcol, &
                                     vpcol_shift=start_k_ordered, &
                                     shifting='L')
            ENDIF
            !
            ! Set coordinates
            IF (do_square_layers3D) THEN
               ! Use the temporary buffers for the communication of the first tick
               IF (first_k) THEN
                  iright_buffer_comm = nbuffers
                  ileft_buffer_comm = nbuffers
               ELSE
                  iright_buffer_comm = icol3D
                  ileft_buffer_comm = irow3D
               ENDIF
            ELSE
               IF (do_comm_right(icol3D)) THEN
                  iright_buffer_comm = MOD(iright_buffer_comm, nbuffers)+1
               ENDIF
               IF (do_comm_left(irow3D)) THEN
                  ileft_buffer_comm = MOD(ileft_buffer_comm, nbuffers)+1
               ENDIF
            ENDIF
            !
            ! Exit if data are already communicated
            IF (istep_k_comm .NE. istep_k_comm_curr) EXIT
            !
            right_buffer_p => right_buffers(iright_buffer_comm)
            left_buffer_p => left_buffers(ileft_buffer_comm)
            right_buffer_p%coord3D = icol3D
            left_buffer_p%coord3D = irow3D
            !
            ! First row, communicate right matrix
            IF (do_comm_right(icol3D)) THEN
               right_buffer_p%vprow = MOD(recv_vrow+v_ki, right_row_total_nimages)
               right_buffer_p%vpcol = right_vcol(icol3D)
               right_buffer_p%is_comm = .FALSE.
            ENDIF
            !
            is_not_comm = .TRUE.
            IF (right_images_size_layers3D(imeta, icol3D, right_buffer_p%vprow) .NE. 0) THEN
               ! First col, communicate left matrix
               IF (do_comm_left(irow3D)) THEN
                  left_buffer_p%vprow = left_vrow(irow3D)
                  left_buffer_p%vpcol = MOD(recv_vcol+v_ki, left_col_total_nimages)
                  left_buffer_p%is_comm = .FALSE.
               ENDIF
               !
               IF (left_images_size_layers3D(imeta, irow3D, left_buffer_p%vpcol) .NE. 0) THEN
                  ! Check if data is already communicated
                  is_not_comm = do_comm_right(icol3D) .OR. do_comm_left(irow3D)
                  IF (is_not_comm) THEN
                     ! Right
                     IF (do_comm_right(icol3D)) THEN
                        IF (has_acc) THEN
                           CALL timeset(routineN//"_acc_sync_right", handle2)
                           CALL acc_event_synchronize(right_buffer_p%data%d%acc_ready)
                           CALL timestop(handle2)
                        ENDIF
                        !
                        do_comm_right(icol3D) = .FALSE.
                        CALL rma_transfer(right_buffer_p%vprow, right_row_nimages, &
                                          right_images_size_layers3D(:, icol3D, right_buffer_p%vprow), &
                                          right_displ_layers3D(:, icol3D, right_buffer_p%vprow), &
                                          right_buffer_p, &
                                          buffers_win%right%meta_win, buffers_win%right%data_win, &
                                          data_get, data_type_byte, buffers_win%right, icol3D, ncols3D)
                     ENDIF
                     ! Left
                     IF (do_comm_left(irow3D)) THEN
                        IF (has_acc) THEN
                           CALL timeset(routineN//"_acc_sync_left", handle2)
                           CALL acc_event_synchronize(left_buffer_p%data%d%acc_ready)
                           CALL timestop(handle2)
                        ENDIF
                        !
                        do_comm_left(irow3D) = .FALSE.
                        CALL rma_transfer(left_buffer_p%vpcol, left_col_nimages, &
                                          left_images_size_layers3D(:, irow3D, left_buffer_p%vpcol), &
                                          left_displ_layers3D(:, irow3D, left_buffer_p%vpcol), &
                                          left_buffer_p, &
                                          buffers_win%left%meta_win, buffers_win%left%data_win, &
                                          data_get, data_type_byte, buffers_win%left, irow3D, nrows3D)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            !
            istep_k_comm_curr = istep_k_comm_curr+1
            ! Stop looping when data is communicated
            ! Only works for 4 layers
            IF (is_not_comm .OR. nranks3D .NE. 4) THEN
               istep_k_comm => istep_k
               IF ((istep_k_comm+1) .EQ. istep_k_comm_curr) EXIT
               ! Restore coordinates by looping once again
               shift3D = shift3D_comm
               CYCLE
            ENDIF
            ! Keep looping until it starts a new communication (only for 4 layers)
            istep_k_comm => istep_k_comm_curr
         ENDDO xfer
         !
         ! Create matrices and multrec's, only the first occurrence
         IF (.NOT. ASSOCIATED(product_matrix3D%mats(irow3D, icol3D)%matrix)) THEN
            IF (irow3D .EQ. myrow3D .AND. icol3D .EQ. mycol3D) THEN
               product_matrix3D%mats(irow3D, icol3D)%matrix => product_matrix
            ELSE
               ALLOCATE (product_matrix3D%mats(irow3D, icol3D)%matrix)
               IF (keep_sparsity) THEN
                  size_guess = product_matrix_meta(product_matrix_meta_displ(irow3D, icol3D)+ &
                                                   dbcsr_slot_nze)
                  CALL setup_buffer_matrix(product_matrix3D%mats(irow3D, icol3D)%matrix, &
                                           product_matrix, product_matrix_meta_size(irow3D, icol3D), &
                                           data_size=size_guess, &
                                           data_memory_type=memtype_mpi_product)
                  product_matrix3D%mats(irow3D, icol3D)% &
                     matrix%index(1:product_matrix_meta_size(irow3D, icol3D)) = &
                     product_matrix_meta(product_matrix_meta_displ(irow3D, icol3D)+1: &
                                         product_matrix_meta_displ(irow3D, icol3D)+ &
                                         product_matrix_meta_size(irow3D, icol3D))
                  CALL dbcsr_data_clear(product_matrix3D%mats(irow3D, icol3D)%matrix%data_area, &
                                        ub=size_guess)
               ELSE
                  CALL setup_buffer_matrix(product_matrix3D%mats(irow3D, icol3D)%matrix, &
                                           product_matrix, data_memory_type=memtype_mpi_product)
               ENDIF
               product_matrix3D%mats(irow3D, icol3D)%matrix%index(dbcsr_slot_home_prow) = &
                  (irow3D-1)*layers_3D_C_reduction%side3D+myprow
               product_matrix3D%mats(irow3D, icol3D)%matrix%index(dbcsr_slot_home_pcol) = &
                  (icol3D-1)*layers_3D_C_reduction%side3D+mypcol
               CALL dbcsr_reset_locals(product_matrix3D%mats(irow3D, icol3D)%matrix)
               product_matrix3D%mats(irow3D, icol3D)%matrix%nblks = 0
               CALL dbcsr_repoint_index(product_matrix3D%mats(irow3D, icol3D)%matrix)
            ENDIF
            !
            IF (.NOT. ASSOCIATED(m_sizes(irow3D)%sizes)) THEN
               ALLOCATE (m_sizes(irow3D)%sizes(dbcsr_nblkrows_local(product_matrix3D%mats(irow3D, icol3D)%matrix)))
               CALL local_filter(array_data(product_matrix3D%mats(irow3D, icol3D)%matrix%row_blk_size), &
                                 array_size(product_matrix3D%mats(irow3D, icol3D)%matrix%local_rows), &
                                 array_data(product_matrix3D%mats(irow3D, icol3D)%matrix%local_rows), &
                                 m_sizes(irow3D)%sizes)
            ENDIF
            IF (.NOT. ASSOCIATED(n_sizes(icol3D)%sizes)) THEN
               ALLOCATE (n_sizes(icol3D)%sizes(dbcsr_nblkcols_local(product_matrix3D%mats(irow3D, icol3D)%matrix)))
               CALL local_filter(array_data(product_matrix3D%mats(irow3D, icol3D)%matrix%col_blk_size), &
                                 array_size(product_matrix3D%mats(irow3D, icol3D)%matrix%local_cols), &
                                 array_data(product_matrix3D%mats(irow3D, icol3D)%matrix%local_cols), &
                                 n_sizes(icol3D)%sizes)
            ENDIF
            !
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          PRIVATE (size_guess, ithread) &
!$OMP          SHARED (product_matrix3D, multrec, &
!$OMP                  keep_sparsity, filter_eps, &
!$OMP                  product_matrix_epss, &
!$OMP                  matrix_right, matrix_left, nthreads, &
!$OMP                  irow3D, icol3D, myrow3D, mycol3D, keep_product_data, &
!$OMP                  product_matrix_epss_displ, product_matrix_epss_size, &
!$OMP                  memtype_product_wm, size_guess_init, nranks3D, m_sizes, n_sizes)
            !
            ! Setup product work areas
            !
            ithread = 0
!$          ithread = OMP_GET_THREAD_NUM()
            !
            IF (irow3D .NE. myrow3D .OR. icol3D .NE. mycol3D) THEN
               IF (keep_product_data) THEN
                  CALL dbcsr_add_wm_from_matrix(product_matrix3D%mats(irow3D, icol3D)%matrix)
               ELSE
                  CALL dbcsr_work_create(product_matrix3D%mats(irow3D, icol3D)%matrix, &
                                         work_mutable=.FALSE., memory_type=memtype_product_wm(ithread)%p)
               ENDIF
!$OMP BARRIER
            ENDIF
            ! The work arrays have to be setup
            size_guess = product_matrix3D%mats(irow3D, icol3D)%matrix%wms(ithread+1)%datasize ! Should be minimal
            IF (.NOT. keep_sparsity) THEN
               size_guess = MAX(size_guess, size_guess_init)
            ENDIF
            CALL dbcsr_data_ensure_size(product_matrix3D%mats(irow3D, icol3D)% &
                                        matrix%wms(ithread+1)%data_area, &
                                        size_guess)
            CALL dbcsr_data_set_size_referenced(product_matrix3D%mats(irow3D, icol3D)% &
                                                matrix%wms(ithread+1)%data_area, &
                                                product_matrix3D%mats(irow3D, icol3D)% &
                                                matrix%wms(ithread+1)%datasize)
            CALL ensure_array_size(product_matrix3D%mats(irow3D, icol3D)% &
                                   matrix%wms(ithread+1)%row_i, ub=1)
            CALL ensure_array_size(product_matrix3D%mats(irow3D, icol3D)% &
                                   matrix%wms(ithread+1)%col_i, ub=1)
            CALL ensure_array_size(product_matrix3D%mats(irow3D, icol3D)% &
                                   matrix%wms(ithread+1)%blk_p, ub=1)
            ALLOCATE (multrec(ithread, irow3D, icol3D)%p)
            CALL dbcsr_mm_multrec_init(multrec(ithread, irow3D, icol3D)%p, &
                                       product=product_matrix3D%mats(irow3D, icol3D)%matrix, &
                                       keep_sparsity=keep_sparsity, &
                                       eps=filter_eps, &
                                       row_max_epss=product_matrix_epss(product_matrix_epss_displ(irow3D)+1: &
                                                                        product_matrix_epss_displ(irow3D)+ &
                                                                        product_matrix_epss_size(irow3D)), &
                                       block_estimate=0, &
                                       right_row_blk_size=dbcsr_row_block_sizes(matrix_right), &
                                       m_sizes=m_sizes(irow3D)%sizes, n_sizes=n_sizes(icol3D)%sizes, &
                                       nlayers=nranks3D)
!$OMP END PARALLEL
            !
            product_matrix3D%mats(irow3D, icol3D)%matrix%nblks = 0
            product_matrix3D%mats(irow3D, icol3D)%matrix%nze = 0
            product_matrix3D%mats(irow3D, icol3D)%matrix%row_p(:) = 0
            CALL dbcsr_data_set_size_referenced(product_matrix3D%mats(irow3D, icol3D)%matrix%data_area, 0)
            product_matrix3D%mats(irow3D, icol3D)%matrix%valid = .FALSE.
         ENDIF
         !
         ! Do the multiplication.  Exclude the first interation.
         calc: IF (istep_k .GT. leftovers_first_k) THEN
            right_buffer_p => right_buffers(iright_buffer_calc)
            left_buffer_p => left_buffers(ileft_buffer_calc)
            irow3D = left_buffer_p%coord3D
            icol3D = right_buffer_p%coord3D
            IF (istep_k .GT. final_step_k) THEN
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (multrec, irow3D, icol3D, irow3D_send, icol3D_send, &
!$OMP         istep_k, final_step_k, product_matrix3D, &
!$OMP         handle2, requests_reduction_size, nthreads, &
!$OMP         product_matrix_meta_send, product_matrix_meta_recv, &
!$OMP         product_matrix_size_send, product_matrix_size_recv, &
!$OMP         buffers_win, data_send, data_get, proc3D_send, proc3D_recv, &
!$OMP         layers_3D_C_reduction, requests_reduction, &
!$OMP         dbcsr_mpi_statistics, data_type_byte) &
!$OMP PRIVATE (ithread)
               ithread = 0
!$             ithread = omp_get_thread_num()
               CALL dbcsr_mm_multrec_phaseout(multrec(ithread, irow3D, icol3D)%p)
               ! Prepare data to send for 3D layer
               IF (istep_k .GT. final_step_k+1) THEN
                  CALL dbcsr_mm_multrec_finalize( &
                     multrec(ithread, irow3D_send, icol3D_send)%p, &
                     buffers_win%left%meta_red3D)
!$OMP BARRIER
!$OMP MASTER
                  CALL timeset(routineN//"_red3D_size", handle2)
                  CALL mp_waitall(requests_reduction_size)
                  CALL timestop(handle2)
                  CALL ensure_array_size(buffers_win%right%meta_red3D, &
                                         ub=product_matrix_size_recv(1), &
                                         nocopy=.TRUE.)
                  product_matrix_meta_send => &
                     buffers_win%left%meta_red3D(1:product_matrix_size_send(1))
                  product_matrix_meta_recv => &
                     buffers_win%right%meta_red3D(1:product_matrix_size_recv(1))
                  CALL mp_isendrecv(product_matrix_meta_send, proc3D_send, &
                                    product_matrix_meta_recv, proc3D_recv, &
                                    layers_3D_C_reduction%grp3D, &
                                    requests_reduction(1), requests_reduction(2))
                  DO myt = 1, nthreads
                     CALL dbcsr_data_ensure_size(layers_3D_C_reduction%data_red3D(myt), &
                                                 product_matrix_size_recv(myt+1), &
                                                 nocopy=.TRUE.)
                     CALL dbcsr_data_set_pointer( &
                        area=data_send, &
                        rsize=product_matrix_size_send(myt+1), &
                        csize=1, &
                        pointee=product_matrix3D%mats(irow3D_send, icol3D_send)%matrix%wms(myt)%data_area)
                     CALL dbcsr_data_set_pointer( &
                        area=data_get, &
                        rsize=product_matrix_size_recv(myt+1), &
                        csize=1, &
                        pointee=layers_3D_C_reduction%data_red3D(myt))
                     CALL dbcsr_isendrecv_any(data_send, proc3D_send, &
                                              data_get, proc3D_recv, &
                                              layers_3D_C_reduction%grp3D, &
                                              requests_reduction(3+(myt-1)*2), &
                                              requests_reduction(4+(myt-1)*2))
                     CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(1, :), &
                                               product_matrix_size_send(myt+1), &
                                               data_type_byte, &
                                               dbcsr_mpi_statistics%data_size_breakdown(:, :, 1))
                  ENDDO
!$OMP END MASTER
               ENDIF
!$OMP END PARALLEL
            ENDIF
            !
            IF (right_buffer_p%is_comm .AND. left_buffer_p%is_comm) THEN
               iright_buffer_calc = MIN(iright_buffer_calc, nbuffers_norms)
               ileft_buffer_calc = MIN(ileft_buffer_calc, nbuffers_norms)
               ! check if right matrix was already initialized
               IF (.NOT. right_buffer_p%matrix%valid) THEN
                  IF (has_acc) CALL dbcsr_data_host2dev(right_buffer_p%data)
                  ! Repoint indices of matrices
                  CALL make_meta(right_buffer_p, &
                                 right_row_total_nimages, &
                                 right_buffer_p%vprow, &
                                 right_buffer_p%vpcol, &
                                 imgdist=imgdist_right, do_merge_rows=.FALSE., &
                                 global_indices=right_global_indices)
                  CALL ensure_array_size(k_sizes, ub=array_size(right_buffer_p%matrix%local_rows))
                  CALL local_filter(array_data(right_buffer_p%matrix%row_blk_size), &
                                    array_size(right_buffer_p%matrix%local_rows), &
                                    array_data(right_buffer_p%matrix%local_rows), &
                                    k_sizes)
                  IF (otf_filtering) THEN
                     CALL calculate_norms(right_buffer_p%matrix, &
                                          right_norms(:, iright_buffer_calc), &
                                          k_sizes, n_sizes(icol3D)%sizes)
                  ENDIF
                  IF (has_acc) THEN
                     CALL acc_transpose_blocks(right_buffer_p%matrix, &
                                               right_buffer_p%trs_stackbuf, &
                                               k_sizes, n_sizes(icol3D)%sizes, &
                                               row_blk_sizes2enum, enum2row_blk_sizes, &
                                               col_blk_sizes2enum, enum2col_blk_sizes, &
                                               noresize=.TRUE.)
                  ENDIF
               ENDIF
               ! check if left matrix was already initialized
               IF (.NOT. left_buffer_p%matrix%valid) THEN
                  IF (has_acc) CALL dbcsr_data_host2dev(left_buffer_p%data)
                  ! Repoint indices of matrices
                  CALL make_meta(left_buffer_p, &
                                 left_col_total_nimages, &
                                 left_buffer_p%vprow, &
                                 left_buffer_p%vpcol, &
                                 imgdist=imgdist_left, do_merge_rows=.TRUE., &
                                 global_indices=left_global_indices, &
                                 nthreads=nthreads)
                  IF (otf_filtering) THEN
                     CALL calculate_norms(left_buffer_p%matrix, &
                                          left_norms(:, ileft_buffer_calc), &
                                          m_sizes(irow3D)%sizes, k_sizes)
                  ENDIF
               ENDIF
               !
               CALL timeset(routineN//"_multrec", handle2)
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (left_buffer_p, ileft_buffer_calc, &
!$OMP         right_buffer_p, iright_buffer_calc, &
!$OMP         left_norms,right_norms, &
!$OMP         multrec, irow3D, icol3D, handle2, k_sizes) &
!$OMP PRIVATE (ithread) &
!$OMP REDUCTION (+: flop)
               ithread = 0
!$             ithread = omp_get_thread_num()
               CALL dbcsr_mm_multrec_multiply(multrec(ithread, irow3D, icol3D)%p, &
                                              left=left_buffer_p%matrix, &
                                              right=right_buffer_p%matrix, &
                                              flop=flop, &
                                              a_norms=left_norms(:, ileft_buffer_calc), &
                                              b_norms=right_norms(:, iright_buffer_calc), &
                                              k_sizes=k_sizes)
!$OMP END PARALLEL
               CALL timestop(handle2)
            ENDIF
            ! Reduce 3D layers and finalize the local layer
            IF (istep_k .GT. final_step_k) THEN
               ! Wait for the other 3D layers to reduce
               IF (istep_k .GT. final_step_k+1) THEN
                  CALL timeset(routineN//"_red3D_data", handle2)
                  CALL mp_waitall(requests_reduction)
                  CALL timestop(handle2)
                  DO myt = 0, nthreads-1
                     DEALLOCATE (multrec(myt, irow3D_send, icol3D_send)%p)
                     CALL dbcsr_work_destroy( &
                        product_matrix3D%mats(irow3D_send, icol3D_send)%matrix%wms(myt+1))
                  ENDDO
                  DEALLOCATE (product_matrix3D%mats(irow3D_send, icol3D_send)%matrix%wms)
                  CALL dbcsr_release(product_matrix3D%mats(irow3D_send, icol3D_send)%matrix)
               ENDIF
               irow3D_send = irow3D
               icol3D_send = icol3D
               ! Store the initial shift for the recv node
               IF (istep_k .EQ. final_step_k+1) THEN
                  shift3D_recv = shift3D-4
               ENDIF
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (multrec, irow3D, icol3D, product_matrix3D, &
!$OMP         memtype_mpi_buffer, nthreads, myt, istep_k, &
!$OMP         irow3D_send, icol3D_send, myrow3D, mycol3D, &
!$OMP         last_step_k, proc3D_send, proc3D_recv, &
!$OMP         product_matrix_size_send, product_matrix_size_recv, &
!$OMP         nrows3D, ncols3D, shift3D_recv, myrank3D, &
!$OMP         layers_3D_C_reduction, requests_reduction_size, &
!$OMP         final_step_k, handle2, buffers_win, g2l_map_rows, g2l_map_cols) &
!$OMP PRIVATE (ithread) &
!$OMP REDUCTION (+: flop)
               ithread = 0
!$             ithread = omp_get_thread_num()
               !
               ! Evaluate the size of layers to send and set the buffers
               IF (irow3D .NE. myrow3D .OR. &
                   icol3D .NE. mycol3D) THEN
                  CALL dbcsr_mm_multrec_dev2host_init(multrec(ithread, irow3D, icol3D)%p)
!$OMP ATOMIC
                  product_matrix3D%mats(irow3D_send, icol3D_send)%matrix%nblks = &
                     product_matrix3D%mats(irow3D_send, icol3D_send)%matrix%nblks+ &
                     dbcsr_mm_multrec_get_nblks(multrec(ithread, irow3D_send, icol3D_send)%p)
!$OMP BARRIER
!$OMP MASTER
                  ! First (nthreads+1) positions are reserved for
                  ! the offset sizes of each thread for meta
                  CALL ensure_array_size(buffers_win%left%meta_red3D, &
                                         ub=product_matrix3D%mats(irow3D_send, icol3D_send)% &
                                         matrix%nblks*3+(nthreads+1), &
                                         nocopy=.TRUE.)
                  ! Set the offsets
                  buffers_win%left%meta_red3D(1) = nthreads+1
                  DO myt = 0, nthreads-1
                     buffers_win%left%meta_red3D(myt+2) = &
                        buffers_win%left%meta_red3D(myt+1)+ &
                        dbcsr_mm_multrec_get_nblks(multrec(myt, irow3D_send, icol3D_send)%p)*3
                     product_matrix_size_send(myt+2) = &
                        dbcsr_mm_multrec_get_nze(multrec(myt, irow3D_send, icol3D_send)%p)
                  ENDDO
                  ! Send/recv data and meta sizes
                  product_matrix_size_send(1) = &
                     buffers_win%left%meta_red3D(nthreads+1)
                  proc3D_send = (icol3D_send-1)*nrows3D+irow3D_send-1
                  !
                  CALL row_col_3D_reflected(irow3D, icol3D, nrows3D, ncols3D, shift3D_recv)
                  shift3D_recv = shift3D_recv-1
                  proc3D_recv = (icol3D-1)*nrows3D+irow3D-1
                  CALL mp_isendrecv(product_matrix_size_send, proc3D_send, &
                                    product_matrix_size_recv, proc3D_recv, &
                                    layers_3D_C_reduction%grp3D, &
                                    requests_reduction_size(1), &
                                    requests_reduction_size(2))
!$OMP END MASTER
               ELSE
                  IF (istep_k .NE. last_step_k) &
                     CPABORT("Last layer does not correspond to local layer")
               ENDIF
               ! Reduce to the local layer
               IF (istep_k .GT. final_step_k+1) THEN
                  IF (dbcsr_data_get_size_referenced(layers_3D_C_reduction%data_red3D(ithread+1)) .GT. 0) THEN
                     CALL timeset(routineN//"_red3D", handle2)
                     CALL dbcsr_mm_multrec_red3D(multrec(ithread, myrow3D, mycol3D)%p, &
                                                 buffers_win%right%meta_red3D, &
                                                 layers_3D_C_reduction%data_red3D(ithread+1), &
                                                 flop, g2l_map_rows, g2l_map_cols)
                     CALL timestop(handle2)
                  ENDIF
               ENDIF
!$OMP END PARALLEL
            ENDIF
         ENDIF calc
         !
         ! Swap temporary buffers for the first tick
         IF (do_square_layers3D .AND. first_k .AND. &
             istep_k .LT. last_step_k) THEN
            iright_buffer_comm = right_buffers(iright_buffer_comm)%coord3D
            ileft_buffer_comm = left_buffers(ileft_buffer_comm)%coord3D
            CALL swap_buffers(right_buffers(iright_buffer_comm), right_buffers(nbuffers))
            CALL swap_buffers(left_buffers(ileft_buffer_comm), left_buffers(nbuffers))
         ENDIF
         !
         iright_buffer_calc = iright_buffer_comm
         ileft_buffer_calc = ileft_buffer_comm
      ENDDO grouped_steps_index
      !
      CALL timestop(handle1)
      !
      CALL m_memory(mem)
      max_memory = MAX(max_memory, REAL(mem))
      !
      IF (do_layers3D .AND. keep_sparsity) THEN
         DEALLOCATE (product_matrix_meta_size, product_matrix_meta_displ)
         DEALLOCATE (product_matrix_meta)
      ENDIF
      DEALLOCATE (right_norms, left_norms)
      DEALLOCATE (product_matrix_epss_size, product_matrix_epss_displ)
      IF (.NOT. otf_filtering .OR. (do_layers3D .AND. nrows3D .GT. 1)) THEN
         DEALLOCATE (product_matrix_epss)
      ELSE
         DEALLOCATE (row_max_epss)
      ENDIF
      !
      DEALLOCATE (left_images_size, right_images_size)
      NULLIFY (left_images_size, right_images_size)
      DEALLOCATE (left_images_size_layers3D, left_displ_layers3D)
      DEALLOCATE (right_images_size_layers3D, right_displ_layers3D)
      !
      ! Deallocate 3D layers
      IF (do_layers3D) THEN
         DEALLOCATE (product_matrix_size_send, product_matrix_size_recv)
         DEALLOCATE (requests_reduction)
         DO icol3D = 1, ncols3D
            DO irow3D = 1, nrows3D
               IF (irow3D .NE. myrow3D .OR. icol3D .NE. mycol3D) THEN
                  DEALLOCATE (product_matrix3D%mats(irow3D, icol3D)%matrix)
               ENDIF
            ENDDO
         ENDDO
         CALL dbcsr_data_clear_pointer(data_send)
         CALL dbcsr_data_release(data_send)
      ENDIF
      DEALLOCATE (product_matrix3D%mats)
      ! Finalize local layer
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (multrec, myrow3D, mycol3D) &
!$OMP PRIVATE (ithread)
      ithread = 0
!$    ithread = omp_get_thread_num()
      CALL dbcsr_mm_multrec_finalize(multrec(ithread, myrow3D, mycol3D)%p)
      DEALLOCATE (multrec(ithread, myrow3D, mycol3D)%p)
!$OMP END PARALLEL
      DEALLOCATE (multrec)
      DEALLOCATE (g2l_map_rows, g2l_map_cols)
      CALL dbcsr_finalize(product_matrix, reshuffle=PRESENT(filter_eps) .AND. .NOT. keep_sparsity)
      !
      DO irow3D = 1, nrows3D
         DEALLOCATE (m_sizes(irow3D)%sizes)
      ENDDO
      DEALLOCATE (m_sizes)
      DO icol3D = 1, ncols3D
         DEALLOCATE (n_sizes(icol3D)%sizes)
      ENDDO
      DEALLOCATE (n_sizes)
      IF (ASSOCIATED(k_sizes)) DEALLOCATE (k_sizes)
      !
      CALL dbcsr_data_clear_pointer(data_get)
      CALL dbcsr_data_release(data_get)
      !
      ! clean-up of communication buffers
      DO v_ki = 1, nbuffers
         CALL dbcsr_data_clear_pointer(left_buffers(v_ki)%data)
         IF (left_buffers(v_ki)%data%d%memory_type%acc_devalloc) THEN
            CALL acc_event_destroy(left_buffers(v_ki)%data%d%acc_ready)
         ENDIF
         CALL dbcsr_data_release(left_buffers(v_ki)%data)
         NULLIFY (left_buffers(v_ki)%matrix%index)
         CALL dbcsr_release(left_buffers(v_ki)%matrix)
         !
         CALL dbcsr_data_clear_pointer(right_buffers(v_ki)%data)
         IF (right_buffers(v_ki)%data%d%memory_type%acc_devalloc) THEN
            CALL acc_event_destroy(right_buffers(v_ki)%data%d%acc_ready)
         ENDIF
         CALL dbcsr_data_release(right_buffers(v_ki)%data)
         NULLIFY (right_buffers(v_ki)%matrix%index)
         CALL dbcsr_release(right_buffers(v_ki)%matrix)
         IF (has_acc) THEN
            CALL dbcsr_data_clear_pointer(right_buffers(v_ki)%trs_stackbuf)
            IF (right_buffers(v_ki)%trs_stackbuf%d%memory_type%acc_devalloc) THEN
               CALL acc_event_destroy(right_buffers(v_ki)%trs_stackbuf%d%acc_ready)
            ENDIF
            CALL dbcsr_data_release(right_buffers(v_ki)%trs_stackbuf)
         ENDIF
      ENDDO
      DEALLOCATE (left_buffers, right_buffers)
      DEALLOCATE (do_comm_left, do_comm_right)
      DEALLOCATE (right_vcol, left_vrow)
      IF (has_acc) THEN
         DEALLOCATE (row_blk_sizes2enum, enum2row_blk_sizes)
         DEALLOCATE (col_blk_sizes2enum, enum2col_blk_sizes)
      ENDIF
      !
      IF (debug_mod) THEN
         v_ki = 0
         DO blk = 1, SIZE(product_matrix%blk_p)
            v_ki = MAX(v_ki, ABS(product_matrix%blk_p(blk)))
         ENDDO
         WRITE (*, *) routineN//" Actual final size", &
            LOG(REAL(dbcsr_data_get_size(product_matrix%data_area)))/LOG(10.0), &
            LOG(REAL(v_ki))/LOG(10.0)
      ENDIF
      !
      CALL timestop(handle)
   END SUBROUTINE multiply_3D

! **************************************************************************************************
!> \brief ...
!> \param buffer ...
!> \param do_win_create ...
!> \param request ...
! **************************************************************************************************
   SUBROUTINE win_setup(buffer, do_win_create, request)
      TYPE(dbcsr_buffer), INTENT(INOUT)                  :: buffer
      LOGICAL, DIMENSION(:), INTENT(INOUT)               :: do_win_create
      INTEGER, INTENT(INOUT)                             :: request

      CHARACTER(len=*), PARAMETER :: routineN = 'win_setup', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, handle1

      CALL timeset(routineN, handle)

      IF (buffer%has_rma_win) THEN
         CALL timeset(routineN//"_win_check", handle1)
         CALL mp_wait(request)
         CALL timestop(handle1)
         IF (do_win_create(1)) THEN
            CALL mp_win_unlock_all(buffer%data_win)
            CALL mp_win_free(buffer%data_win)
         ENDIF
         IF (do_win_create(2)) THEN
            CALL mp_win_unlock_all(buffer%meta_win)
            CALL mp_win_free(buffer%meta_win)
         ENDIF
      ENDIF
      CALL dbcsr_data_release(buffer%data_before_resize)
      IF (ASSOCIATED(buffer%meta_before_resize)) THEN
         CALL memory_deallocate(buffer%meta_before_resize, memtype_mpi_buffer)
         NULLIFY (buffer%meta_before_resize)
      ENDIF
      !
      IF (do_win_create(1)) THEN
         CALL dbcsr_win_create_any(buffer%data, buffer%subgrp, buffer%data_win)
         CALL mp_win_lock_all(buffer%data_win)
      ENDIF
      IF (do_win_create(2)) THEN
         CALL mp_win_create(buffer%meta, buffer%subgrp, buffer%meta_win)
         CALL mp_win_lock_all(buffer%meta_win)
      ENDIF
      !
      buffer%has_rma_win = .TRUE.
      CALL timestop(handle)
   END SUBROUTINE win_setup

! **************************************************************************************************
!> \brief Apply reflected order, i.e. row increasing value for odd col value,
!>                                    row decreasing value for even col value
!> \param irow3D ...
!> \param icol3D ...
!> \param nrows3D ...
!> \param ncols3D ...
!> \param shift3D ...
! **************************************************************************************************
   SUBROUTINE row_col_3D_reflected(irow3D, icol3D, nrows3D, ncols3D, shift3D)
      INTEGER, INTENT(INOUT)                             :: irow3D, icol3D
      INTEGER, INTENT(IN)                                :: nrows3D, ncols3D, shift3D

      INTEGER                                            :: odd_or_even

      icol3D = MOD(shift3D/nrows3D, ncols3D)+1
      irow3D = MOD(shift3D, nrows3D)
      odd_or_even = MOD(icol3D, 2)
      irow3D = (nrows3D-irow3D)*(1-odd_or_even)+(irow3D+1)*odd_or_even
   END SUBROUTINE row_col_3D_reflected

! **************************************************************************************************
!> \brief ...
!> \param buffer_1 ...
!> \param buffer_2 ...
!> \param buffers ...
!> \param nbuffers ...
!> \param data_size ...
!> \param meta_size ...
!> \param matrix ...
!> \param imgdist ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE setup_buffers(buffer_1, buffer_2, buffers, nbuffers, data_size, meta_size, matrix, imgdist)
      TYPE(dbcsr_buffer), INTENT(INOUT), TARGET          :: buffer_1, buffer_2
      TYPE(dbcsr_buffer), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: buffers
      INTEGER, INTENT(IN)                                :: nbuffers, data_size, meta_size
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(dbcsr_imagedistribution_obj), INTENT(INOUT)   :: imgdist

      INTEGER                                            :: ibuffer, jbuffer
      INTEGER, DIMENSION(:), POINTER                     :: meta_p
      LOGICAL                                            :: has_trs_stackbuf
      TYPE(dbcsr_buffer), POINTER                        :: buffer_p

      ALLOCATE (buffers(nbuffers))
      has_trs_stackbuf = dbcsr_data_valid(buffer_1%trs_stackbuf) .OR. dbcsr_data_valid(buffer_2%trs_stackbuf)
      DO ibuffer = 1, nbuffers
         jbuffer = (ibuffer-1)/2
         IF (MOD(ibuffer, 2) .EQ. 1) THEN
            buffer_p => buffer_1
         ELSE
            buffer_p => buffer_2
         ENDIF
         ! Use slices for the 3D buffers
         CALL dbcsr_data_init(buffers(ibuffer)%data)
         CALL dbcsr_data_new(buffers(ibuffer)%data, dbcsr_data_get_type(buffer_p%data), &
                             memory_type=dbcsr_data_get_memory_type(buffer_p%data))
         IF (buffers(ibuffer)%data%d%memory_type%acc_devalloc) THEN
            CALL acc_event_create(buffers(ibuffer)%data%d%acc_ready)
         ENDIF
         CALL dbcsr_data_set_pointer( &
            area=buffers(ibuffer)%data, &
            rsize=data_size, &
            csize=1, &
            pointee=buffer_p%data, &
            source_lb=data_size*jbuffer+1)
         ! Use meta_p pointer to avoid warning target-lifetime
         meta_p => buffer_p%meta(meta_size*jbuffer+1: &
                                 meta_size*(jbuffer+1))
         buffers(ibuffer)%meta => meta_p
         IF (has_trs_stackbuf) THEN
            CALL dbcsr_data_init(buffers(ibuffer)%trs_stackbuf)
            CALL dbcsr_data_new(buffers(ibuffer)%trs_stackbuf, dbcsr_data_get_type(buffer_p%trs_stackbuf), &
                                memory_type=dbcsr_data_get_memory_type(buffer_p%trs_stackbuf))
            IF (buffers(ibuffer)%trs_stackbuf%d%memory_type%acc_devalloc) THEN
               CALL acc_event_create(buffers(ibuffer)%trs_stackbuf%d%acc_ready)
            ENDIF
            CALL dbcsr_data_set_pointer( &
               area=buffers(ibuffer)%trs_stackbuf, &
               rsize=meta_size/3, &
               csize=1, &
               pointee=buffer_p%trs_stackbuf, &
               source_lb=(meta_size/3)*jbuffer+1)
         ENDIF
         CALL setup_buffer_matrix_image(buffers(ibuffer)%matrix, imgdist, matrix, &
                                        buffers(ibuffer)%data, &
                                        buffers(ibuffer)%meta)
      ENDDO
   END SUBROUTINE setup_buffers

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param imgdist ...
!> \param template_matrix ...
!> \param data_buffer ...
!> \param meta_buffer ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE setup_buffer_matrix_image(matrix, imgdist, &
                                        template_matrix, data_buffer, &
                                        meta_buffer)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      TYPE(dbcsr_imagedistribution_obj), INTENT(INOUT)   :: imgdist
      TYPE(dbcsr_type), INTENT(IN)                       :: template_matrix
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: data_buffer
      INTEGER, DIMENSION(:), INTENT(IN), TARGET          :: meta_buffer

      matrix = dbcsr_type()
      CALL dbcsr_create(matrix, &
                        "Buffer image of "//template_matrix%name, &
                        imgdist%i%main, &
                        dbcsr_type_no_symmetry, &
                        array_data(template_matrix%row_blk_size), array_data(template_matrix%col_blk_size), &
                        template_matrix%row_blk_size, template_matrix%col_blk_size, &
                        data_type=dbcsr_data_get_type(data_buffer), &
                        data_buffer=data_buffer, &
                        max_rbs=template_matrix%max_rbs, max_cbs=template_matrix%max_cbs, &
                        row_blk_offset=template_matrix%row_blk_offset, &
                        col_blk_offset=template_matrix%col_blk_offset, &
                        index_memory_type=memtype_mpi_buffer, &
                        make_index=.FALSE.)
      matrix%index => meta_buffer
      matrix%negate_real = template_matrix%negate_real
      matrix%negate_imaginary = template_matrix%negate_imaginary
      matrix%local_indexing = .TRUE.
      matrix%list_indexing = .TRUE.
   END SUBROUTINE setup_buffer_matrix_image

! **************************************************************************************************
!> \brief ...
!> \param buffers_1 ...
!> \param buffers_2 ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE swap_buffers(buffers_1, buffers_2)
      TYPE(dbcsr_buffer), INTENT(INOUT)                  :: buffers_1, buffers_2

      TYPE(dbcsr_buffer)                                 :: tmp

      tmp = buffers_1
      buffers_1 = buffers_2
      buffers_2 = tmp
   END SUBROUTINE swap_buffers

! **************************************************************************************************
!> \brief ...
!> \param recv_vproc ...
!> \param nimages ...
!> \param size_layers3D ...
!> \param displ_layers3D ...
!> \param buffer ...
!> \param meta_win ...
!> \param data_win ...
!> \param data_get ...
!> \param data_type_byte ...
!> \param buffer_win ...
!> \param layer3D ...
!> \param nlayers3D ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE rma_transfer(recv_vproc, nimages, &
                           size_layers3D, displ_layers3D, &
                           buffer, &
                           meta_win, data_win, &
                           data_get, data_type_byte, &
                           buffer_win, layer3D, nlayers3D)
      INTEGER, INTENT(IN)                                :: recv_vproc, nimages
      INTEGER, DIMENSION(:), INTENT(IN)                  :: size_layers3D, displ_layers3D
      TYPE(dbcsr_buffer), INTENT(INOUT)                  :: buffer
      INTEGER, INTENT(IN)                                :: meta_win, data_win
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: data_get
      INTEGER, INTENT(IN)                                :: data_type_byte
      TYPE(dbcsr_buffer), INTENT(IN)                     :: buffer_win
      INTEGER, INTENT(IN)                                :: layer3D, nlayers3D

      INTEGER                                            :: recv_proc
      INTEGER, DIMENSION(:), POINTER                     :: meta_get

      buffer%is_comm = .TRUE.
      buffer%get_requests(:) = mp_request_null
      recv_proc = (recv_vproc/nimages)*nlayers3D+layer3D-1
      !
      meta_get => buffer%meta(dbcsr_num_slots+1:dbcsr_num_slots+size_layers3D(imeta))
      buffer%meta_size = size_layers3D(imeta)
      CALL mp_rget(meta_get, recv_proc, &
                   meta_win, &
                   buffer_win%meta, &
                   disp=displ_layers3D(imeta), &
                   request=buffer%get_requests(1))
      CALL dbcsr_data_set_pointer( &
         area=data_get, &
         rsize=size_layers3D(idata), &
         csize=1, &
         pointee=buffer%data, &
         source_lb=1)
      CALL dbcsr_rget_any(data_get, recv_proc, &
                          data_win, &
                          buffer_win%data, &
                          disp=displ_layers3D(idata), &
                          request=buffer%get_requests(2))
      CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(1, :), &
                                size_layers3D(idata), &
                                data_type_byte, &
                                dbcsr_mpi_statistics%data_size_breakdown(:, :, 1))
      dbcsr_mpi_statistics%nexchanged = dbcsr_mpi_statistics%nexchanged+1
      !
      ! Set the referenced sizes to the actual data moved via MPI
      CALL dbcsr_data_set_size_referenced(buffer%data, size_layers3D(idata))
      buffer%matrix%valid = .FALSE.
   END SUBROUTINE rma_transfer

! **************************************************************************************************
!> \brief ...
!> \param meta_buffer ...
!> \param img_nblks_rows ...
!> \param img_nblks_cols ...
!> \param refs_size ...
!> \param refs_displ ...
!> \param size_index ...
!> \param has_threads ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE setup_rec_index_images(meta_buffer, img_nblks_rows, img_nblks_cols, &
                                     refs_size, refs_displ, size_index, has_threads)
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: meta_buffer
      INTEGER, DIMENSION(:), INTENT(IN)                  :: img_nblks_rows, img_nblks_cols, &
                                                            refs_size, refs_displ
      INTEGER, INTENT(IN)                                :: size_index
      LOGICAL, INTENT(IN)                                :: has_threads

      CHARACTER(len=*), PARAMETER :: routineN = 'setup_rec_index_images', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, in, nblkcols_local, &
                                                            nblkrows_local, t_f, t_l, t_size

!$    INTEGER                           :: ithread

      CALL timeset(routineN, handle)
      IF (has_threads) THEN
         nblkrows_local = img_nblks_rows(1)
      ELSE
         nblkcols_local = img_nblks_cols(1)
      ENDIF
      !
      DO in = 1, SIZE(refs_size)
         IF (refs_size(in) .EQ. 0) CYCLE
         ! Number of blocks
         t_size = (refs_size(in)-size_index)/3
         IF (has_threads) THEN
            nblkcols_local = img_nblks_cols(in)
         ELSE
            nblkrows_local = img_nblks_rows(in)
         ENDIF
         t_f = 1
         t_l = t_size
!$OMP    PARALLEL IF (has_threads) DEFAULT (NONE) &
!$OMP    PRIVATE (ithread) &
!$OMP    FIRSTPRIVATE (t_f, t_l, t_size) &
!$OMP    SHARED (meta_buffer, in, has_threads, refs_displ, &
!$OMP            size_index, nblkrows_local, nblkcols_local)
!$       ithread = OMP_GET_THREAD_NUM()+1
!$       IF (has_threads) THEN
!$          t_f = meta_buffer(refs_displ(in)+ithread)+1
!$          t_l = meta_buffer(refs_displ(in)+ithread+1)
!$       ENDIF
         t_size = t_l-t_f+1
         IF (t_size .GT. 0) THEN
            CALL rec_sort_index(1, nblkrows_local, &
                                1, nblkcols_local, &
                                t_size, &
                                meta_buffer(refs_displ(in)+size_index+t_f*3-2: &
                                            refs_displ(in)+size_index+t_l*3), &
                                0)
         ENDIF
!$OMP    END PARALLEL
      ENDDO
      CALL timestop(handle)
   END SUBROUTINE setup_rec_index_images

! **************************************************************************************************
!> \brief Init buffer
!>
!> \param buffer ...
!> \param data_type ...
!> \param data_size ...
!> \param meta_size ...
!> \param num_data ...
!> \param data_memory_type ...
!> \param trs_memory_type ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE buffer_init(buffer, data_type, &
                          data_size, meta_size, &
                          num_data, &
                          data_memory_type, trs_memory_type)
      TYPE(dbcsr_buffer), INTENT(INOUT)                  :: buffer
      INTEGER, INTENT(IN)                                :: data_type, data_size, meta_size
      INTEGER, INTENT(IN), OPTIONAL                      :: num_data
      TYPE(dbcsr_memtype_type), INTENT(IN)               :: data_memory_type
      TYPE(dbcsr_memtype_type), INTENT(IN), OPTIONAL     :: trs_memory_type

      INTEGER                                            :: my_num_data
      LOGICAL                                            :: new_trs_stackbuf

      my_num_data = 1
      IF (PRESENT(num_data)) THEN
         my_num_data = num_data
      ELSE
         IF (dbcsr_data_valid(buffer%data_before_resize) .OR. ASSOCIATED(buffer%meta_before_resize)) &
            CPABORT("Previous data area already initialized.")
         CALL dbcsr_data_init(buffer%data_before_resize)
         CALL dbcsr_data_new(buffer%data_before_resize, data_type, memory_type=data_memory_type)
      ENDIF
      new_trs_stackbuf = PRESENT(trs_memory_type) .AND. has_acc
      !
      IF (buffer%is_valid) THEN
         ! Invalid buffers if data_type is different
         IF (dbcsr_data_get_type(buffer%data) .NE. data_type) THEN
            CALL dbcsr_data_release(buffer%data)
            IF (new_trs_stackbuf) THEN
               CALL dbcsr_data_release(buffer%trs_stackbuf)
            ENDIF
            buffer%is_valid = .FALSE.
         ENDIF
      ENDIF
      !
      IF (.NOT. buffer%is_valid) THEN
         ! First initialization
         CALL dbcsr_data_init(buffer%data)
         CALL dbcsr_data_new(buffer%data, data_type=data_type, &
                             data_size=data_size*my_num_data, memory_type=data_memory_type)
         CALL dbcsr_data_set_size_referenced(buffer%data, data_size*my_num_data)
         IF (new_trs_stackbuf) THEN
            CALL dbcsr_data_init(buffer%trs_stackbuf)
            CALL dbcsr_data_new(buffer%trs_stackbuf, &
                                data_type=dbcsr_type_int_4, data_size=(meta_size/3)*my_num_data, &
                                memory_type=trs_memory_type)
         ENDIF
         buffer%is_valid = .TRUE.
      ELSE
         IF (PRESENT(num_data)) THEN
            CALL dbcsr_data_ensure_size(buffer%data, data_size*my_num_data, nocopy=.TRUE.)
            IF (new_trs_stackbuf) THEN
               CALL dbcsr_data_ensure_size(buffer%trs_stackbuf, (meta_size/3)*my_num_data, nocopy=.TRUE.)
            ENDIF
         ELSE
            ! Case for MPI windows
            ! data_before_resize keeps the pointer to previous data in the case of reallocation
            CALL dbcsr_data_ensure_size(buffer%data, data_size, nocopy=.TRUE., &
                                        area_resize=buffer%data_before_resize)
         ENDIF
      ENDIF
      !
      IF (PRESENT(num_data)) THEN
         CALL ensure_array_size(buffer%meta, ub=meta_size*my_num_data, nocopy=.TRUE., &
                                memory_type=memtype_mpi_buffer)
      ELSE
         ! Case for MPI windows
         ! meta_before_resize keeps the pointer to previous meta in the case of reallocation
         CALL ensure_array_size(buffer%meta, array_resize=buffer%meta_before_resize, &
                                ub=meta_size, nocopy=.TRUE., &
                                memory_type=memtype_mpi_buffer)
      ENDIF
      !
      buffer%is_comm = .FALSE.
   END SUBROUTINE buffer_init

! **************************************************************************************************
!> \brief Release all buffers
!>
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE buffers_release()
      IF (request_sync_mult .NE. mp_request_null) CALL mp_wait(request_sync_mult)
      request_sync_mult = mp_request_null
      CALL buffer_release(buffers_1%right)
      CALL buffer_release(buffers_1%left)
      CALL buffer_release(buffers_2%right)
      CALL buffer_release(buffers_2%left)
      CALL buffer_release(buffers_win%right)
      CALL buffer_release(buffers_win%left)
      !
      IF (dbcsr_data_valid(make_buffers_data_send)) CALL dbcsr_data_release(make_buffers_data_send)
      IF (dbcsr_data_valid(make_buffers_data_recv)) CALL dbcsr_data_release(make_buffers_data_recv)
      IF (ASSOCIATED(make_buffers_meta_send)) CALL memory_deallocate(make_buffers_meta_send, memtype_mpi_buffer)
      IF (ASSOCIATED(make_buffers_meta_recv)) CALL memory_deallocate(make_buffers_meta_recv, memtype_mpi_buffer)
   END SUBROUTINE buffers_release

! **************************************************************************************************
!> \brief Release buffer
!>
!> \param buffer ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE buffer_release(buffer)
      TYPE(dbcsr_buffer), INTENT(INOUT)                  :: buffer

      IF (buffer%has_rma_win) THEN
         CALL mp_win_unlock_all(buffer%data_win)
         CALL mp_win_free(buffer%data_win)
         CALL mp_win_unlock_all(buffer%meta_win)
         CALL mp_win_free(buffer%meta_win)
         buffer%has_rma_win = .FALSE.
         buffer%grp = mp_comm_null
         IF (buffer%subgrp .NE. mp_comm_null .AND. buffer%num_layers_3D .GT. 1) &
            CALL mp_comm_free(buffer%subgrp)
         buffer%subgrp = mp_comm_null
         buffer%num_layers_3D = 1
      ENDIF
      !
      IF (buffer%is_valid) THEN
         CALL dbcsr_data_release(buffer%data)
         IF (dbcsr_data_valid(buffer%trs_stackbuf)) THEN
            CALL dbcsr_data_release(buffer%trs_stackbuf)
         ENDIF
         IF (dbcsr_data_valid(buffer%data_before_resize)) THEN
            CALL dbcsr_data_release(buffer%data_before_resize)
         ENDIF
         buffer%is_valid = .FALSE.
      ENDIF
      IF (ASSOCIATED(buffer%meta)) THEN
         CALL memory_deallocate(buffer%meta, memtype_mpi_buffer)
         NULLIFY (buffer%meta)
      ENDIF
      IF (ASSOCIATED(buffer%meta_before_resize)) THEN
         CALL memory_deallocate(buffer%meta_before_resize, memtype_mpi_buffer)
         NULLIFY (buffer%meta_before_resize)
      ENDIF
      IF (ASSOCIATED(buffer%meta_red3D)) THEN
         CALL memory_deallocate(buffer%meta_red3D, memtype_mpi_buffer)
         NULLIFY (buffer%meta_red3D)
      ENDIF
   END SUBROUTINE buffer_release

! **************************************************************************************************
!> \brief Create meta indices
!>
!> \param buffer ...
!> \param ntotal_images ...
!> \param vprow ...
!> \param vpcol ...
!> \param imgdist ...
!> \param do_merge_rows ...
!> \param global_indices ...
!> \param nthreads ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE make_meta(buffer, ntotal_images, &
                        vprow, vpcol, &
                        imgdist, do_merge_rows, &
                        global_indices, &
                        nthreads)
      TYPE(dbcsr_buffer), INTENT(INOUT)                  :: buffer
      INTEGER, INTENT(IN)                                :: ntotal_images, vprow, vpcol
      TYPE(dbcsr_imagedistribution_obj), INTENT(INOUT)   :: imgdist
      LOGICAL, INTENT(IN)                                :: do_merge_rows
      INTEGER, DIMENSION(:), INTENT(IN)                  :: global_indices
      INTEGER, INTENT(IN), OPTIONAL                      :: nthreads

      CHARACTER(len=*), PARAMETER :: routineN = 'make_meta', routineP = moduleN//':'//routineN

      INTEGER                                            :: size_index

      buffer%matrix%index(dbcsr_slot_size) = &
         buffer%meta_size+dbcsr_num_slots
      size_index = 0
      IF (PRESENT(nthreads)) THEN
!$       size_index = nthreads+1
      ENDIF
      buffer%matrix%index(dbcsr_slot_nblks) = &
         (buffer%meta_size-size_index)/3
      buffer%matrix%index(dbcsr_slot_nze) = &
         dbcsr_data_get_size_referenced(buffer%data)
      buffer%matrix%index(dbcsr_slot_dense) = 0
      buffer%matrix%index(dbcsr_slot_nblkrows_total:dbcsr_slot_nfullcols_local) = &
         global_indices(:)
      buffer%matrix%index(dbcsr_slot_type:dbcsr_num_slots) = 0
      ! Virtual coords
      IF (do_merge_rows) THEN
         buffer%matrix%index(dbcsr_slot_home_vprow) = vprow
         buffer%matrix%index(dbcsr_slot_home_vpcol) = MOD(vpcol, ntotal_images)
      ELSE
         buffer%matrix%index(dbcsr_slot_home_vprow) = MOD(vprow, ntotal_images)
         buffer%matrix%index(dbcsr_slot_home_vpcol) = vpcol
      ENDIF
      buffer%matrix%index(dbcsr_slot_row_p) = 1
      buffer%matrix%index(dbcsr_slot_col_i) = 1
      buffer%matrix%index(dbcsr_slot_blk_p) = 1
      ! thr_c
      size_index = dbcsr_num_slots
!$    IF (PRESENT(nthreads)) THEN
!$       size_index = size_index+nthreads+1
!$       buffer%matrix%index(dbcsr_slot_thr_c) = dbcsr_num_slots+1
!$       buffer%matrix%index(dbcsr_slot_thr_c+1) = size_index
!$    ENDIF
      buffer%matrix%index(dbcsr_slot_coo_l) = size_index+1
      buffer%matrix%index(dbcsr_num_slots) = buffer%matrix%index(dbcsr_slot_size)
      !
      ! Reset
      CALL dbcsr_reset_vlocals(buffer%matrix, imgdist)
      !
      ! Repoint index
      buffer%matrix%nblks = 0
      buffer%matrix%nze = 0
      CALL dbcsr_repoint_index(buffer%matrix)
      buffer%matrix%valid = .TRUE.
   END SUBROUTINE make_meta

! **************************************************************************************************
!> \brief Remap the 4-rank array in a 3-rank array by introducing the virtual coordinate
!> \param refs_size ...
!> \param refs_size_layers3D ...
!> \param refs_displ_layers3D ...
!> \param data_size ...
!> \param meta_size ...
! **************************************************************************************************
   SUBROUTINE remap_layers3D(refs_size, refs_size_layers3D, refs_displ_layers3D, &
                             data_size, meta_size)
      INTEGER, DIMENSION(:, :, :, :), INTENT(IN), &
         POINTER                                         :: refs_size
      INTEGER, DIMENSION(:, :, :), INTENT(OUT), POINTER  :: refs_size_layers3D, refs_displ_layers3D
      INTEGER, INTENT(OUT)                               :: data_size, meta_size

      INTEGER                                            :: ilayer, image, iproc, nimages, &
                                                            nlayers3D, nprocs

      nimages = SIZE(refs_size, 2)
      nlayers3D = SIZE(refs_size, 3)
      nprocs = SIZE(refs_size, 4)
      !
      ALLOCATE (refs_size_layers3D(idata:imeta, nlayers3D, 0:nimages*nprocs-1))
      ALLOCATE (refs_displ_layers3D(idata:imeta, nlayers3D, 0:nimages*nprocs-1))
      data_size = 0; meta_size = 0
      !
!$OMP PARALLEL DO DEFAULT (NONE) &
!$OMP SHARED (nprocs, nimages, nlayers3D, &
!$OMP         refs_size_layers3D, refs_displ_layers3D, refs_size) &
!$OMP PRIVATE (iproc,image,ilayer) &
!$OMP REDUCTION (MAX : data_size, meta_size)
      DO iproc = 0, nprocs-1
         DO ilayer = 1, nlayers3D
            DO image = 1, nimages
               refs_size_layers3D(:, ilayer, image+iproc*nimages-1) = refs_size(:, image, ilayer, iproc)
               data_size = MAX(data_size, refs_size(idata, image, ilayer, iproc))
               meta_size = MAX(meta_size, refs_size(imeta, image, ilayer, iproc))
            ENDDO
            refs_displ_layers3D(:, ilayer, iproc*nimages) = 0
            DO image = 1, nimages-1
               refs_displ_layers3D(:, ilayer, image+iproc*nimages) = &
                  refs_displ_layers3D(:, ilayer, image+iproc*nimages-1)+refs_size(:, image, ilayer, iproc)
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
   END SUBROUTINE remap_layers3D

! **************************************************************************************************
!> \brief ...
!> \retval get_max_layers_3D ...
! **************************************************************************************************
   PURE FUNCTION get_max_layers_3D()
      INTEGER                                            :: get_max_layers_3D

      get_max_layers_3D = layers_3D_C_reduction%max_num_layers_3D

   END FUNCTION get_max_layers_3D

END MODULE dbcsr_mm_3d

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   3D matrix-matrix multiplication.
!>
!> \author  Alfio Lazzaro
!>
!> <b>Modification history:</b>
!>  - 2016-08    Code organization (Alfio Lazzaro).
! **************************************************************************************************
MODULE dbcsr_mm_3d
   USE acc_event,                       ONLY: acc_event_record,&
                                              acc_event_synchronize,&
                                              acc_stream_wait_event
   USE array_types,                     ONLY: array_data,&
                                              array_get
   USE dbcsr_acc_operations,            ONLY: dbcsr_acc_transpose
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
        dbcsr_data_clear_pointer, dbcsr_data_ensure_size, dbcsr_data_exists, dbcsr_data_get_size, &
        dbcsr_data_get_size_referenced, dbcsr_data_get_type, dbcsr_data_host2dev, dbcsr_data_init, &
        dbcsr_data_new, dbcsr_data_release, dbcsr_data_set_pointer, &
        dbcsr_data_set_size_referenced, dbcsr_data_valid, dbcsr_get_data_p_c, dbcsr_get_data_p_d, &
        dbcsr_get_data_p_s, dbcsr_get_data_p_z, dbcsr_scalar_are_equal, dbcsr_scalar_negative, &
        dbcsr_scalar_one, dbcsr_type_1d_to_2d
   USE dbcsr_data_types,                ONLY: dbcsr_datatype_sizeof
   USE dbcsr_dist_methods,              ONLY: dbcsr_distribution_col_dist,&
                                              dbcsr_distribution_has_threads,&
                                              dbcsr_distribution_local_cols,&
                                              dbcsr_distribution_local_rows,&
                                              dbcsr_distribution_mp,&
                                              dbcsr_distribution_row_dist,&
                                              dbcsr_distribution_thread_dist
   USE dbcsr_index_operations,          ONLY: dbcsr_repoint_index
   USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop
   USE dbcsr_mem_methods,               ONLY: dbcsr_mempool_limit_capacity
   USE dbcsr_methods,                   ONLY: &
        dbcsr_col_block_offsets, dbcsr_col_block_sizes, dbcsr_destroy_array, dbcsr_distribution, &
        dbcsr_get_data_type, dbcsr_get_index_memory_type, dbcsr_has_symmetry, &
        dbcsr_image_dist_hold, dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_nfullcols_local, &
        dbcsr_nfullcols_total, dbcsr_nfullrows_local, dbcsr_nfullrows_total, dbcsr_release, &
        dbcsr_row_block_offsets, dbcsr_row_block_sizes, dbcsr_valid_index
   USE dbcsr_mm_common,                 ONLY: &
        count_mpi_statistics, dbcsr_mm_multrec_type_p, dbcsr_mpi_statistics, enumerate_blk_sizes, &
        max_memory, memtype_abpanel_1, memtype_abpanel_2, memtype_mpi_buffer, memtype_mpi_product, &
        memtype_product_wm, memtype_trsbuffer_1, memtype_trsbuffer_2, product_matrix_size_guess, &
        rec_sort_index, setup_buffer_matrix
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
                                              dbcsr_sendrecv_any,&
                                              dbcsr_win_create_any,&
                                              hybrid_alltoall_any,&
                                              hybrid_alltoall_i1
   USE dbcsr_ptr_util,                  ONLY: ensure_array_size,&
                                              memory_deallocate
   USE dbcsr_types,                     ONLY: &
        dbcsr_1d_array_type, dbcsr_2d_array_obj, dbcsr_data_obj, dbcsr_distribution_obj, &
        dbcsr_imagedistribution_obj, dbcsr_iterator, dbcsr_memtype_type, dbcsr_mp_obj, &
        dbcsr_num_slots, dbcsr_scalar_type, dbcsr_slot_blk_p, dbcsr_slot_col_i, dbcsr_slot_coo_l, &
        dbcsr_slot_home_pcol, dbcsr_slot_home_prow, dbcsr_slot_home_vpcol, dbcsr_slot_home_vprow, &
        dbcsr_slot_nblkrows_total, dbcsr_slot_nblks, dbcsr_slot_nfullcols_local, dbcsr_slot_nze, &
        dbcsr_slot_row_p, dbcsr_slot_size, dbcsr_slot_thr_c, dbcsr_type, dbcsr_type_complex_4, &
        dbcsr_type_complex_8, dbcsr_type_int_4, dbcsr_type_no_symmetry, dbcsr_type_real_4, &
        dbcsr_type_real_8
   USE dbcsr_util,                      ONLY: find_block_of_element
   USE dbcsr_work_operations,           ONLY: dbcsr_add_wm_from_matrix,&
                                              dbcsr_create,&
                                              dbcsr_finalize,&
                                              dbcsr_work_create,&
                                              dbcsr_work_destroy
   USE kinds,                           ONLY: int_8,&
                                              real_4,&
                                              real_8,&
                                              sp
   USE machine,                         ONLY: m_memory
   USE message_passing,                 ONLY: &
        mp_allgather, mp_alltoall, mp_comm_free, mp_comm_null, mp_comm_split_direct, &
        mp_iallgather, mp_isendrecv, mp_isum, mp_request_null, mp_rget, mp_sendrecv, mp_wait, &
        mp_waitall, mp_win_create, mp_win_free, mp_win_lock_all, mp_win_unlock_all
#include "../../base/base_uses.f90"

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads,&
!$                    omp_set_lock, omp_unset_lock, omp_init_lock, omp_lock_kind, omp_destroy_lock

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_3d'
   LOGICAL, PARAMETER :: debug_mod = .FALSE.
   LOGICAL, PARAMETER :: careful_mod = .FALSE.

   TYPE dbcsr_buffer
      TYPE(dbcsr_data_obj), DIMENSION(:), ALLOCATABLE :: DATA
      TYPE(dbcsr_data_obj)                            :: data_red3D, data_resize
      INTEGER                                         :: vprow, vpcol
      INTEGER                                         :: grp = mp_comm_null, & ! Global communicator
                                                         subgrp = mp_comm_null ! Communicator for A and B
      INTEGER                                         :: data_win, meta_win
      INTEGER, DIMENSION(:), POINTER                  :: meta => Null(), &
                                                         meta_resize => Null(), &
                                                         meta_red3D => Null()
      INTEGER, DIMENSION(:, :), ALLOCATABLE           :: get_requests, &
                                                         offset
      TYPE(dbcsr_data_obj), DIMENSION(:), ALLOCATABLE :: trs_stackbuf
      INTEGER                            :: nrequests = 0
      INTEGER                            :: num_layers_3D = 1
      INTEGER                            :: coord3D = 1
      TYPE(dbcsr_1d_array_type)          :: buffer
      LOGICAL                            :: is_valid = .FALSE., &
                                            has_rma_win = .FALSE.
   END TYPE dbcsr_buffer

   TYPE dbcsr_buffer_p
      TYPE(dbcsr_buffer), POINTER        :: b => Null()
   END TYPE dbcsr_buffer_p

   TYPE dbcsr_buffers
      TYPE(dbcsr_buffer)             :: left, right
   END TYPE dbcsr_buffers

   TYPE dbcsr_layers_3D_C_reduction
      INTEGER                            :: grp = mp_comm_null, &
                                            grp3D = mp_comm_null, &
                                            rowgrp3D = mp_comm_null
      INTEGER                            :: num_layers_3D = 1, &
                                            side3D = HUGE(1)
   END TYPE dbcsr_layers_3D_C_reduction

   ! Buffers
   TYPE(dbcsr_buffers), TARGET, SAVE :: buffers_orig, &
                                        buffers_1, buffers_2

   INTEGER, PARAMETER, PRIVATE               :: idata = 1, &
                                                imeta = 2, &
                                                ilocal_proc = 1, &
                                                isym_proc = 2

   TYPE(dbcsr_layers_3D_C_reduction), SAVE :: layers_3D_C_reduction

   LOGICAL, DIMENSION(2), TARGET, SAVE, PRIVATE       :: do_win_create_left, &
                                                         do_win_create_right

   INTEGER, ALLOCATABLE, DIMENSION(:), PRIVATE :: left_total_row_counts
   INTEGER, ALLOCATABLE, DIMENSION(:), TARGET, PRIVATE :: left_refs_data_size, &
                                                          right_refs_data_size, &
                                                          left_local_refs_data_size, &
                                                          right_local_refs_data_size
   INTEGER, ALLOCATABLE, DIMENSION(:, :), TARGET, PRIVATE  :: left_no_empty_images, &
                                                              left_no_empty_images_displ, &
                                                              right_no_empty_images, &
                                                              right_no_empty_images_displ, &
                                                              left_local_refs_displ_unmerged, &
                                                              right_local_refs_displ_unmerged, &
                                                              left_local_refs_meta_size, &
                                                              right_local_refs_meta_size
   INTEGER, ALLOCATABLE, DIMENSION(:, :, :), PRIVATE  :: left_refs_displ_unmerged, &
                                                         right_refs_displ_unmerged, &
                                                         left_refs_meta_size, &
                                                         right_refs_meta_size
   REAL(kind=sp), ALLOCATABLE, DIMENSION(:), TARGET, PRIVATE :: left_max_norms, right_max_norms, &
                                                                left_local_max_norms, right_local_max_norms

   INTEGER, ALLOCATABLE, DIMENSION(:), TARGET, PRIVATE :: g2l_map_cols, g2l_map_rows

   INTEGER, PRIVATE                :: request_count_rows
   INTEGER, DIMENSION(8), PRIVATE  :: requests
   INTEGER, DIMENSION(2), PRIVATE  :: request_win_create

   PUBLIC :: multiply_3D
   PUBLIC :: release_layers_3d_C_reduction, buffers_release
   PUBLIC :: dbcsr_make_buffers, make_layers_3d_C_reduction

CONTAINS

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
!> \param transpose ...
!> \param alpha ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE dbcsr_make_buffers(matrix, imgdist, is_left, &
                                 f_row, l_row, f_col, l_col, &
                                 otf_filtering, transpose, &
                                 alpha)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(dbcsr_imagedistribution_obj), INTENT(IN)      :: imgdist
      LOGICAL, INTENT(IN)                                :: is_left
      INTEGER, INTENT(IN)                                :: f_row, l_row, f_col, l_col
      LOGICAL, INTENT(IN)                                :: otf_filtering, transpose
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
                           otf_filtering, transpose, &
                           alpha)
      ELSE
         CALL make_buffers(matrix, imgdist, is_left, &
                           f_row, l_row, f_col, l_col, &
                           otf_filtering, transpose)
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
!> \param transpose ...
!> \param scale_value ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE make_buffers(matrix, imgdist, is_left, &
                           f_row, l_row, f_col, l_col, &
                           otf_filtering, transpose, scale_value)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(dbcsr_imagedistribution_obj), INTENT(IN)      :: imgdist
      LOGICAL, INTENT(IN)                                :: is_left
      INTEGER, INTENT(IN)                                :: f_row, l_row, f_col, l_col
      LOGICAL, INTENT(IN)                                :: otf_filtering, transpose
      TYPE(dbcsr_scalar_type), INTENT(IN), OPTIONAL      :: scale_value

      CHARACTER(len=*), PARAMETER :: routineN = 'make_buffers', routineP = moduleN//':'//routineN

      INTEGER :: blk, blk_p, bp, col, col_img, col_mult, col_size, data_type, dst_proc, f_col_f, &
         f_row_f, grp, handle, handle2, ilayer, iproc, irequests, it, ithread, l_col_l, l_row_l, &
         local_no_empty_images, mynode, mypcol, myprow, myt, nblkcols_local, nblkrows_local, &
         nblocks, ncols_images, nimages_merged, nimages_unmerged, nprocs, nprocs_sym, &
         nrows_images, nsymmetries, nthreads, nze, pcol, prow, request_no_empty_images, row, &
         row_img, row_mult, row_size, size_index_unmerged, stored_col, stored_row, sym_p, &
         symmetry_i, total_no_empty_images, tr_col_size, tr_row_size
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: img_map, img_refs, img_refs_cols, &
                                                            img_refs_rows, tmp_img_offset
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: recv_displ, recv_size, send_displ, &
                                                            send_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: local_refs_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :, :)     :: recv_refs, refs_displ, refs_size
      INTEGER, DIMENSION(2)                              :: block_col_bounds, block_row_bounds
      INTEGER, DIMENSION(:), POINTER :: col_dist, col_img_dist, local_cols, local_g2l_map_cols, &
         local_g2l_map_rows, local_refs_data_size, local_rows, meta_buffer_p, meta_recv, &
         meta_send, refs_data_size, row_dist, row_img_dist, threads_dist
      INTEGER, DIMENSION(:, :), POINTER                  :: blacs2mpi, local_refs_displ_unmerged, &
                                                            local_refs_meta_size, no_empty_images, &
                                                            no_empty_images_displ
      INTEGER, DIMENSION(:, :, :), POINTER               :: local_refs_displ
      INTEGER, DIMENSION(dbcsr_slot_nblkrows_total:&
         dbcsr_slot_nfullcols_local)                     :: global_indices
      INTEGER, DIMENSION(idata:imeta)                    :: local_size
      INTEGER, POINTER                                   :: coli, rowi
      INTEGER, TARGET                                    :: mi, ui
      LOGICAL :: do_crop, do_part_crop_col, do_part_crop_f_col, do_part_crop_f_row, &
         do_part_crop_l_col, do_part_crop_l_row, do_part_crop_row, do_symmetry, do_virt, &
         is_diagonal, tr
      LOGICAL, DIMENSION(:), POINTER                     :: do_win_create
      REAL(kind=sp), DIMENSION(:), POINTER               :: local_max_norms, max_norms
      TYPE(dbcsr_buffer), POINTER                        :: buffer
      TYPE(dbcsr_data_obj)                               :: data_block
      TYPE(dbcsr_data_obj), POINTER                      :: data_buffer_p
      TYPE(dbcsr_data_obj), TARGET                       :: data_recv, data_send
      TYPE(dbcsr_distribution_obj)                       :: set_dist
      TYPE(dbcsr_iterator)                               :: iter
      TYPE(dbcsr_mp_obj)                                 :: mp_obj
      TYPE(dbcsr_scalar_type)                            :: scale_neg_one
      TYPE(dbcsr_type)                                   :: sm

!$    INTEGER(kind=omp_lock_kind), ALLOCATABLE, DIMENSION(:) :: locks

      CALL timeset(routineN, handle)
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
      row_mult = imgdist%i%row_multiplicity
      row_img_dist => array_data(imgdist%i%row_image)
      ncols_images = imgdist%i%col_decimation
      col_mult = imgdist%i%col_multiplicity
      col_img_dist => array_data(imgdist%i%col_image)
      mp_obj = dbcsr_distribution_mp(imgdist%i%main)
      CALL dbcsr_mp_grid_setup(mp_obj)
      grp = dbcsr_mp_group(mp_obj)
      blacs2mpi => dbcsr_mp_pgrid(mp_obj)
      mynode = dbcsr_mp_mynode(mp_obj)
      myprow = dbcsr_mp_myprow(mp_obj)
      mypcol = dbcsr_mp_mypcol(mp_obj)
      IF (MAXVAL(row_dist) .GT. UBOUND(blacs2mpi, 1)) &
         CPABORT("Row distribution references unexistent processor rows")
      IF (MAXVAL(col_dist) .GT. UBOUND(blacs2mpi, 2)) &
         CPABORT("Col distribution references unexistent processor cols")
      ! Check threads configuration
      NULLIFY (threads_dist)
!$    IF (.NOT. dbcsr_distribution_has_threads(dbcsr_distribution(matrix))) &
!$       CPABORT("Thread distribution not defined")
!$    threads_dist => array_data(dbcsr_distribution_thread_dist(dbcsr_distribution(matrix)))
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
      ! Check for virtual topology
      do_virt = row_mult .NE. nrows_images .OR. col_mult .NE. ncols_images
      !
      IF (do_virt) THEN
         ! For virtual topology we just ignore the symmetric proc
         ! Instead an alltoall communication is done between
         ! all processors.
         sym_p = -1
         nprocs_sym = dbcsr_mp_numnodes(mp_obj)
      ELSE
         sym_p = blacs2mpi(mypcol, myprow)
         IF (do_symmetry .OR. transpose) THEN
            ! Symmetric case requires a communication with the
            ! symmetric processor
            nprocs_sym = 2
         ELSE
            nprocs_sym = 1
         ENDIF
      ENDIF
      !
      is_diagonal = sym_p .EQ. mynode
      !
      size_index_unmerged = dbcsr_slot_nblks
      !
      IF (is_left) THEN
         ! merging over rows
         nimages_merged = nrows_images
         nimages_unmerged = ncols_images
         buffer => buffers_orig%left
         nprocs = dbcsr_mp_npcols(mp_obj)
         ALLOCATE (left_refs_meta_size(myprow*nimages_merged:(myprow+1)*nimages_merged-1, &
                                       MAX(1, dbcsr_mp_nprows(mp_obj)/layers_3D_C_reduction%side3D), &
                                       0:nprocs*nimages_unmerged-1))
         ALLOCATE (left_local_refs_meta_size(nimages_merged, nimages_unmerged))
         left_local_refs_meta_size(:, :) = 0
         local_refs_meta_size => left_local_refs_meta_size
         ALLOCATE (left_local_refs_data_size(nimages_merged*nimages_unmerged))
         local_refs_data_size => left_local_refs_data_size
         ALLOCATE (left_refs_displ_unmerged(idata:imeta, &
                                            UBOUND(left_refs_meta_size, 2), &
                                            0:nprocs*nimages_unmerged-1))
         ALLOCATE (left_local_refs_displ_unmerged(idata:imeta, nimages_unmerged))
         local_refs_displ_unmerged => left_local_refs_displ_unmerged
         requests(7) = mp_request_null
         irequests = 1
         !
         ! Count the maximum possible multiplies per row for on-the-fly filtering
         IF (otf_filtering) THEN
            ALLOCATE (left_total_row_counts(nblkrows_local))
            left_total_row_counts = 0
         ENDIF
         ALLOCATE (left_no_empty_images(UBOUND(left_refs_meta_size, 2), nprocs))
         no_empty_images => left_no_empty_images
         ALLOCATE (left_no_empty_images_displ(UBOUND(left_refs_meta_size, 2), nprocs))
         no_empty_images_displ => left_no_empty_images_displ
         do_win_create => do_win_create_left
      ELSE
         ! merging over cols
         nimages_merged = ncols_images
         nimages_unmerged = nrows_images
         buffer => buffers_orig%right
         nprocs = dbcsr_mp_nprows(mp_obj)
         ALLOCATE (right_refs_meta_size(mypcol*nimages_merged:(mypcol+1)*nimages_merged-1, &
                                        MAX(1, dbcsr_mp_npcols(mp_obj)/layers_3D_C_reduction%side3D), &
                                        0:nprocs*nimages_unmerged-1))
         ALLOCATE (right_local_refs_meta_size(nimages_merged, nimages_unmerged))
         right_local_refs_meta_size(:, :) = 0
         local_refs_meta_size => right_local_refs_meta_size
         ALLOCATE (right_local_refs_data_size(nimages_merged*nimages_unmerged))
         local_refs_data_size => right_local_refs_data_size
         ALLOCATE (right_refs_displ_unmerged(idata:imeta, &
                                             UBOUND(right_refs_meta_size, 2), &
                                             0:nprocs*nimages_unmerged-1))
         ALLOCATE (right_local_refs_displ_unmerged(idata:imeta, nimages_unmerged))
         local_refs_displ_unmerged => right_local_refs_displ_unmerged
         requests(8) = mp_request_null
         irequests = 2
         ALLOCATE (right_no_empty_images(UBOUND(right_refs_meta_size, 2), nprocs))
         no_empty_images => right_no_empty_images
         ALLOCATE (right_no_empty_images_displ(UBOUND(right_refs_meta_size, 2), nprocs))
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
         do_win_create => do_win_create_right
      ENDIF
      !
      ! 3D communicator
      IF (layers_3D_C_reduction%num_layers_3D .GT. 1 .AND. nimages_merged .GT. 1) &
         CPABORT("Cannot make 3D layers with 2D images distribution")
      CALL make_layers_3D_AB(layers_3D_C_reduction%num_layers_3D, &
                             layers_3D_C_reduction%side3D, &
                             mp_obj, is_left, buffer)
      !
      ALLOCATE (local_refs_size(idata:imeta, nimages_merged, nimages_unmerged))
      local_refs_size(:, :, :) = 0
      !
      ! Evaluate maps for global -> local indexing (g2l_map_rows, g2l_map_cols)
      ! Evaluate maps for local_image -> local indexing (img_map)
      ! Evaluate refs (count and offset) per each image
      IF (is_left) THEN
         ALLOCATE (g2l_map_rows(sm%nblkrows_total))
         local_g2l_map_rows => g2l_map_rows
         ALLOCATE (local_g2l_map_cols(sm%nblkcols_total))
      ELSE
         ALLOCATE (local_g2l_map_rows(sm%nblkrows_total))
         ALLOCATE (g2l_map_cols(sm%nblkcols_total))
         local_g2l_map_cols => g2l_map_cols
      ENDIF
      ALLOCATE (img_refs_rows(nimages_unmerged), img_refs_cols(nimages_unmerged))
      !
      local_g2l_map_rows(:) = 0
      IF (is_left .OR. nrows_images .EQ. 1) THEN
         ! Use merge global -> local indexing
         DO row = 1, nblkrows_local
            local_g2l_map_rows(local_rows(row)) = row
         ENDDO
         img_refs_rows(1) = nblkrows_local
      ELSE
         img_refs_rows(:) = 0
         DO row = 1, nblkrows_local
            row_img = row_img_dist(local_rows(row))
            mi = MOD(row_img-1, nrows_images)+1
            img_refs_rows(mi) = img_refs_rows(mi)+1
            local_g2l_map_rows(local_rows(row)) = img_refs_rows(mi)
         ENDDO
      ENDIF
      !
      local_g2l_map_cols(:) = 0
      IF (.NOT. is_left .OR. ncols_images .EQ. 1) THEN
         ! Use merge global -> local indexing
         DO col = 1, nblkcols_local
            local_g2l_map_cols(local_cols(col)) = col
         ENDDO
         img_refs_cols(1) = nblkcols_local
      ELSE
         img_refs_cols(:) = 0
         DO col = 1, nblkcols_local
            col_img = col_img_dist(local_cols(col))
            mi = MOD(col_img-1, ncols_images)+1
            img_refs_cols(mi) = img_refs_cols(mi)+1
            local_g2l_map_cols(local_cols(col)) = img_refs_cols(mi)
         ENDDO
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
      ALLOCATE (local_refs_displ(idata:imeta, nimages_merged, nimages_unmerged))
      local_size(:) = 0
      local_no_empty_images = 0
      nblocks = 0
      !
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP PRIVATE (ithread,myt,iter,row,col,blk,row_size,col_size,&
!$OMP          stored_row,stored_col,blk_p,bp,tr,&
!$OMP          nze,symmetry_i,row_img,col_img,rowi,coli,&
!$OMP          tr_row_size,tr_col_size,prow,pcol,dst_proc,&
!$OMP          data_buffer_p,meta_buffer_p,&
!$OMP          mi,ui,it,data_block) &
!$OMP SHARED (nthreads,refs_size,refs_displ,matrix,nsymmetries,do_symmetry,&
!$OMP         row_img_dist,col_img_dist,imgdist,row_dist,col_dist,&
!$OMP         is_left,local_refs_size,nimages_merged,nimages_unmerged,&
!$OMP         local_size,data_type,memtype_mpi_buffer,sm,&
!$OMP         img_refs_cols,img_refs_rows,img_refs,img_map,&
!$OMP         local_g2l_map_cols,local_g2l_map_rows,recv_refs,grp,meta_send,&
!$OMP         scale_value,scale_neg_one,data_send,data_recv,&
!$OMP         size_index_unmerged,recv_displ,send_displ,recv_size,send_size,&
!$OMP         mp_obj,threads_dist,meta_recv,nrows_images,ncols_images,&
!$OMP         locks,mynode,blacs2mpi,myprow,mypcol,sym_p,&
!$OMP         left_refs_meta_size,right_refs_meta_size,&
!$OMP         left_refs_data_size,right_refs_data_size,&
!$OMP         local_refs_meta_size,local_refs_data_size,&
!$OMP         local_refs_displ,&
!$OMP         left_refs_displ_unmerged,right_refs_displ_unmerged,&
!$OMP         local_refs_displ_unmerged,requests,is_diagonal,&
!$OMP         buffer,left_total_row_counts,otf_filtering,&
!$OMP         iproc,ilayer,&
!$OMP         no_empty_images,local_no_empty_images,transpose,&
!$OMP         request_no_empty_images,irequests,do_win_create,&
!$OMP         do_virt,handle2,nprocs_sym,global_indices,&
!$OMP         do_crop,do_part_crop_row,do_part_crop_col,block_row_bounds,block_col_bounds,&
!$OMP         do_part_crop_f_row,do_part_crop_l_row,do_part_crop_f_col,do_part_crop_l_col,&
!$OMP         f_row_f,l_row_l,f_col_f,l_col_l,nblocks,request_win_create,refs_data_size)
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
      ALLOCATE (refs_size(idata:imeta, &
                          0:nthreads, &
                          nimages_merged, &
                          nimages_unmerged, &
                          nprocs_sym))
      ALLOCATE (refs_displ(idata, &
                           0:nthreads, &
                           nimages_merged, &
                           nimages_unmerged, &
                           nprocs_sym))
      IF (.NOT. do_virt) nprocs_sym = 1
      ALLOCATE (send_size(idata:imeta, nprocs_sym))
      ALLOCATE (recv_size(idata:imeta, nprocs_sym))
      IF (do_virt) THEN
         ALLOCATE (send_displ(idata:imeta, nprocs_sym))
         ALLOCATE (recv_displ(idata:imeta, nprocs_sym))
      ENDIF
      refs_size(:, :, :, :, :) = 0
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
            IF (do_virt) THEN
               dst_proc = blacs2mpi(prow, pcol)+1
            ELSEIF (do_symmetry .OR. transpose) THEN
               ! Avoid unnecessary copy of data for diagonal procs
               IF (blacs2mpi(prow, pcol) .EQ. mynode .OR. is_diagonal) THEN
                  dst_proc = ilocal_proc
               ELSE
                  ! Move data to remote symmetric proc
                  dst_proc = isym_proc
               ENDIF
            ENDIF
!$          IF (is_left .AND. do_symmetry) THEN
!$             myt = threads_dist(stored_row)
!$          ENDIF
!$OMP ATOMIC
            refs_size(imeta, myt+1, mi, ui, dst_proc) = &
               refs_size(imeta, myt+1, mi, ui, dst_proc)+3
!$OMP ATOMIC
            refs_size(idata, myt+1, mi, ui, dst_proc) = &
               refs_size(idata, myt+1, mi, ui, dst_proc)+nze
         ENDDO ! symmetry_i
      ENDDO
      CALL dbcsr_iterator_stop(iter)
      ! Avoid unnecessary copy if data is already local
      IF (do_virt .OR. ((do_symmetry .OR. transpose) .AND. .NOT. is_diagonal)) THEN
!$OMP BARRIER
!$OMP MASTER
         ALLOCATE (recv_refs(idata:imeta, &
                             0:nthreads, &
                             nimages_merged, &
                             nimages_unmerged, &
                             nprocs_sym))
         IF (do_virt) THEN
            CALL mp_alltoall(refs_size(:, :, :, :, :), &
                             recv_refs(:, :, :, :, :), &
                             2*nimages_merged*nimages_unmerged*(nthreads+1), grp)
         ELSE
            ! exchange sizes with the symmetric proc
            CALL mp_sendrecv(refs_size(:, :, :, :, isym_proc), sym_p, &
                             recv_refs(:, :, :, :, ilocal_proc), sym_p, grp)
         ENDIF
!$OMP END MASTER
      ENDIF
!$OMP BARRIER
      !
      ! Store local data and meta dimensions
      IF (.NOT. do_virt) THEN
!$OMP DO COLLAPSE(3)
         DO ui = 1, nimages_unmerged
            DO mi = 1, nimages_merged
               DO it = 1, nthreads
!$OMP ATOMIC
                  local_refs_size(idata, mi, ui) = &
                     local_refs_size(idata, mi, ui)+ &
                     refs_size(idata, it, mi, ui, ilocal_proc)
!$OMP ATOMIC
                  local_refs_size(imeta, mi, ui) = &
                     local_refs_size(imeta, mi, ui)+ &
                     refs_size(imeta, it, mi, ui, ilocal_proc)
               ENDDO
            ENDDO
         ENDDO
!$OMP END DO NOWAIT
      ENDIF
      IF (do_virt .OR. ((do_symmetry .OR. transpose) .AND. (.NOT. is_diagonal))) THEN
!$OMP DO COLLAPSE(4)
         DO dst_proc = 1, nprocs_sym
            DO ui = 1, nimages_unmerged
               DO mi = 1, nimages_merged
                  DO it = 1, nthreads
!$OMP ATOMIC
                     local_refs_size(idata, mi, ui) = &
                        local_refs_size(idata, mi, ui)+ &
                        recv_refs(idata, it, mi, ui, dst_proc)
!$OMP ATOMIC
                     local_refs_size(imeta, mi, ui) = &
                        local_refs_size(imeta, mi, ui)+ &
                        recv_refs(imeta, it, mi, ui, dst_proc)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
!$OMP END DO NOWAIT
      ENDIF
!$OMP BARRIER
!$OMP MASTER
      !
      ! Take the displacements, taking in account the images merging
      ! Note that local diplacements refer to position 0
      DO ui = 1, nimages_unmerged
         local_refs_displ_unmerged(:, ui) = local_size(:)
         DO mi = 1, nimages_merged
            local_refs_displ(:, mi, ui) = local_size(:)
            ! Exclude communications of empty data size only when nimages_merged>1
            IF (local_refs_size(imeta, mi, ui) .EQ. 0 .AND. nimages_merged .GT. 1) CYCLE
            local_no_empty_images = local_no_empty_images+1
            local_refs_data_size(local_no_empty_images) = local_refs_size(idata, mi, ui)
            IF (local_refs_size(imeta, mi, ui) .EQ. 0) CYCLE
            ! Include stats slots
            local_refs_size(imeta, mi, ui) = local_refs_size(imeta, mi, ui)+ &
                                             size_index_unmerged
            local_refs_meta_size(mi, ui) = local_refs_size(imeta, mi, ui)
            local_size(:) = local_size(:)+local_refs_size(:, mi, ui)
         ENDDO
      ENDDO
      !
      ! Exchange refs
      IF (is_left) THEN
         CALL mp_iallgather(local_refs_meta_size, left_refs_meta_size, buffer%subgrp, requests(irequests))
         CALL mp_iallgather(local_refs_displ_unmerged, left_refs_displ_unmerged, &
                            buffer%subgrp, requests(2+irequests))
      ELSE
         CALL mp_iallgather(local_refs_meta_size, right_refs_meta_size, buffer%subgrp, requests(irequests))
         CALL mp_iallgather(local_refs_displ_unmerged, right_refs_displ_unmerged, &
                            buffer%subgrp, requests(2+irequests))
      ENDIF
      IF (nimages_merged .GT. 1) THEN
         CALL mp_iallgather(local_no_empty_images, no_empty_images, buffer%subgrp, request_no_empty_images)
      ENDIF
      !
      ! Allocate data and meta buffers
      do_win_create(:) = .NOT. buffer%has_rma_win
      IF (buffer%has_rma_win) THEN
         IF (buffer%grp .NE. grp .OR. dbcsr_data_get_type(buffer%data(1)) .NE. data_type) THEN
            do_win_create(:) = .TRUE.
         ENDIF
      ENDIF
      CALL buffer_init(buffer, data_type, &
                       local_size(idata), local_size(imeta), &
                       data_memory_type=memtype_mpi_buffer)
      IF (buffer%has_rma_win) THEN
         do_win_create(1) = do_win_create(1) .OR. dbcsr_data_exists(buffer%data_resize)
         do_win_create(2) = do_win_create(2) .OR. ASSOCIATED(buffer%meta_resize)
         CALL mp_isum(do_win_create, buffer%subgrp, request_win_create(irequests))
      ENDIF
      buffer%grp = grp
      !
      IF (.NOT. do_virt) THEN
         !
         ! Change displacements for the first thread
         ! In this case they refer to position 1
         refs_size(:, 0, :, :, ilocal_proc) = local_refs_displ(:, :, :)+1
         !
         ! Evaluate data and meta diplacements per each thread
         DO it = 1, nthreads-1
            refs_size(:, it, :, :, ilocal_proc) = &
               refs_size(:, it, :, :, ilocal_proc)+refs_size(:, it-1, :, :, ilocal_proc)
         ENDDO
      ENDIF
      !
      ! Take offsets for virtual and symmetric case
      IF (do_virt .OR. ((do_symmetry .OR. transpose) .AND. (.NOT. is_diagonal))) THEN
         ! Set recv buffer sizes and offsets
         recv_size(:, 1) = 1
         IF (do_virt) recv_displ(:, 1) = 0
         DO dst_proc = 1, nprocs_sym
            DO ui = 1, nimages_unmerged
               DO mi = 1, nimages_merged
                  ! Empty image
                  IF (local_refs_size(imeta, mi, ui) .EQ. 0) THEN
                     DO it = 0, nthreads
                        recv_refs(:, it, mi, ui, dst_proc) = recv_size(:, dst_proc)
                     ENDDO
                  ELSE
                     ! Carry previous value
                     recv_refs(:, 0, mi, ui, dst_proc) = recv_size(:, dst_proc)
                     DO it = 1, nthreads
                        recv_refs(:, it, mi, ui, dst_proc) = &
                           recv_refs(:, it, mi, ui, dst_proc)+ &
                           recv_refs(:, it-1, mi, ui, dst_proc)
                        IF (.NOT. do_virt) THEN
                           ! Add symetric events to the diplacements of local data
                           refs_size(:, it, mi, ui, ilocal_proc) = &
                              refs_size(:, it, mi, ui, ilocal_proc)+ &
                              recv_refs(:, it, mi, ui, ilocal_proc)-recv_size(:, ilocal_proc)
                        ENDIF
                     ENDDO
                     recv_size(:, dst_proc) = recv_refs(:, nthreads, mi, ui, dst_proc)
                  ENDIF
               ENDDO
            ENDDO
            IF (do_virt) THEN
               ! Carry previous value
               IF (dst_proc .LT. nprocs_sym) THEN
                  recv_size(:, dst_proc+1) = recv_size(:, dst_proc)
                  recv_displ(:, dst_proc+1) = recv_size(:, dst_proc)-1
               ENDIF
               ! Set size for the current proc
               recv_size(:, dst_proc) = recv_size(:, dst_proc)-recv_displ(:, dst_proc)-1
            ENDIF
         ENDDO
         !
         send_size(:, 1) = 1
         IF (do_virt) THEN
            DO dst_proc = 1, nprocs_sym
               send_displ(:, dst_proc) = send_size(:, dst_proc)-1
               DO ui = 1, nimages_unmerged
                  DO mi = 1, nimages_merged
                     ! Carry previous value
                     refs_size(:, 0, mi, ui, dst_proc) = send_size(:, dst_proc)
                     DO it = 1, nthreads
                        refs_size(:, it, mi, ui, dst_proc) = &
                           refs_size(:, it, mi, ui, dst_proc)+ &
                           refs_size(:, it-1, mi, ui, dst_proc)
                     ENDDO
                     send_size(:, dst_proc) = refs_size(:, nthreads, mi, ui, dst_proc)
                  ENDDO
               ENDDO
               IF (dst_proc .LT. nprocs_sym) THEN
                  send_size(:, dst_proc+1) = send_size(:, dst_proc)
               ENDIF
               send_size(:, dst_proc) = send_size(:, dst_proc)-send_displ(:, dst_proc)-1
            ENDDO
         ELSE
            DO ui = 1, nimages_unmerged
               DO mi = 1, nimages_merged
                  ! Carry previous value
                  refs_size(:, 0, mi, ui, isym_proc) = send_size(:, 1)
                  DO it = 1, nthreads
                     refs_size(:, it, mi, ui, isym_proc) = &
                        refs_size(:, it, mi, ui, isym_proc)+ &
                        refs_size(:, it-1, mi, ui, isym_proc)
                  ENDDO
                  send_size(:, 1) = refs_size(:, nthreads, mi, ui, isym_proc)
               ENDDO
            ENDDO
         ENDIF
         !
         ! Allocate data/meta to send
         CALL dbcsr_data_init(data_send)
         CALL dbcsr_data_new(data_send, data_type, SUM(send_size(idata, :)), &
                             memory_type=memtype_mpi_buffer)
         CALL dbcsr_data_clear(data_send)
         NULLIFY (meta_send)
         CALL ensure_array_size(meta_send, ub=SUM(send_size(imeta, :)), &
                                nocopy=.TRUE., memory_type=memtype_mpi_buffer, zero_pad=.TRUE.)
      ENDIF
      !
      refs_displ(idata, :, :, :, :) = refs_size(idata, :, :, :, :)-1
!$OMP END MASTER
!$OMP BARRIER
      !
      IF (do_part_crop_row .OR. do_part_crop_col) THEN
         CALL dbcsr_data_init(data_block)
         CALL dbcsr_data_new(data_block, dbcsr_type_1d_to_2d(data_type))
      ENDIF
      !
      IF (do_virt) THEN
         data_buffer_p => data_send
         meta_buffer_p => meta_send
      ELSE
         data_buffer_p => buffer%data(1)
         meta_buffer_p => buffer%meta
      ENDIF
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
            IF (do_virt) THEN
               dst_proc = blacs2mpi(prow, pcol)+1
            ELSEIF (.NOT. is_diagonal .AND. (do_symmetry .OR. transpose)) THEN
               IF (blacs2mpi(prow, pcol) .EQ. mynode) THEN
                  dst_proc = ilocal_proc
                  data_buffer_p => buffer%data(1)
                  meta_buffer_p => buffer%meta
               ELSE
                  ! Move data to remote symmetric proc
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
                                             refs_size(idata, myt, mi, ui, dst_proc), bp, &
                                             scale_value)
               IF (sm%negate_real .AND. sm%negate_imaginary) THEN
                  CALL dbcsr_block_scale(data_buffer_p, scale=scale_neg_one, &
                                         row_size=nze, col_size=1, &
                                         lb=refs_size(idata, myt, mi, ui, dst_proc))
               ELSEIF (sm%negate_real) THEN
                  CALL dbcsr_block_real_neg(data_buffer_p, row_size=nze, col_size=1, &
                                            lb=refs_size(idata, myt, mi, ui, dst_proc))
               ELSEIF (sm%negate_imaginary) THEN
                  CALL dbcsr_block_conjg(data_buffer_p, row_size=nze, col_size=1, &
                                         lb=refs_size(idata, myt, mi, ui, dst_proc))
               ENDIF
            ELSE
               CALL dbcsr_block_copy_aa(data_buffer_p, sm%data_area, row_size, col_size, &
                                        refs_size(idata, myt, mi, ui, dst_proc), bp, &
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
                  source_lb=refs_size(idata, myt, mi, ui, dst_proc))
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
            ! Set meta data (global indexing)
            IF (do_virt .OR. ((do_symmetry .OR. transpose) .AND. blacs2mpi(prow, pcol) .NE. mynode)) THEN
               meta_buffer_p(refs_size(imeta, myt, mi, ui, dst_proc)) = stored_row
               meta_buffer_p(refs_size(imeta, myt, mi, ui, dst_proc)+1) = stored_col
               ! Reset data position per each proc
               meta_buffer_p(refs_size(imeta, myt, mi, ui, dst_proc)+2) = &
                  refs_size(idata, myt, mi, ui, dst_proc)- &
                  refs_displ(idata, myt, mi, ui, dst_proc)
            ELSE
               meta_buffer_p(refs_size(imeta, myt, mi, ui, dst_proc)+size_index_unmerged) = &
                  stored_row
               meta_buffer_p(refs_size(imeta, myt, mi, ui, dst_proc)+size_index_unmerged+1) = &
                  stored_col
               meta_buffer_p(refs_size(imeta, myt, mi, ui, dst_proc)+size_index_unmerged+2) = &
                  refs_size(idata, myt, mi, ui, dst_proc)-local_refs_displ(idata, mi, ui)
            ENDIF
            refs_size(imeta, myt, mi, ui, dst_proc) = refs_size(imeta, myt, mi, ui, dst_proc)+3
            refs_size(idata, myt, mi, ui, dst_proc) = refs_size(idata, myt, mi, ui, dst_proc)+nze
!$          IF (is_left .AND. do_symmetry) THEN
!$             call omp_unset_lock(locks(myt))
!$          ENDIF
         ENDDO
      ENDDO
      CALL dbcsr_iterator_stop(iter)
      !
      IF (do_part_crop_row .OR. do_part_crop_col) THEN
         CALL dbcsr_data_clear_pointer(data_block)
         CALL dbcsr_data_release(data_block)
      ENDIF
      !
      IF (do_virt .OR. ((do_symmetry .OR. transpose) .AND. .NOT. is_diagonal)) THEN
!$OMP BARRIER
!$OMP MASTER
         CALL dbcsr_data_init(data_recv)
         CALL dbcsr_data_new(data_recv, data_type, SUM(recv_size(idata, :)), &
                             memory_type=memtype_mpi_buffer)
         NULLIFY (meta_recv)
         CALL ensure_array_size(meta_recv, ub=SUM(recv_size(imeta, :)), &
                                nocopy=.TRUE., memory_type=memtype_mpi_buffer, zero_pad=.TRUE.)
         CALL timeset(routineN//"_data", handle2)
         IF (do_virt) THEN
            ! Exchange data
            CALL hybrid_alltoall_any(data_send, send_size(idata, :), send_displ(idata, :), &
                                     data_recv, recv_size(idata, :), recv_displ(idata, :), &
                                     mp_obj, &
                                     most_ptp=.TRUE., remainder_ptp=.TRUE., no_hybrid=.FALSE.)
            CALL hybrid_alltoall_i1(meta_send, send_size(imeta, :), send_displ(imeta, :), &
                                    meta_recv, recv_size(imeta, :), recv_displ(imeta, :), &
                                    mp_obj, &
                                    most_ptp=.TRUE., remainder_ptp=.TRUE., no_hybrid=.FALSE.)
         ELSE
            ! Send/recv data to the symmetric proc
            CALL dbcsr_sendrecv_any(data_send, sym_p, data_recv, sym_p, grp)
            CALL mp_sendrecv(meta_send, sym_p, meta_recv, sym_p, grp)
         ENDIF
         CALL timestop(handle2)
         ! Delete send data
         CALL dbcsr_data_release(data_send)
         CALL memory_deallocate(meta_send, memtype_mpi_buffer)
!$OMP END MASTER
      ENDIF
!$OMP BARRIER
!$    IF (is_left .AND. do_symmetry) THEN
!$       call omp_destroy_lock(locks(ithread))
!$    ENDIF
      !
      meta_buffer_p => buffer%meta
      DO ui = 1, nimages_unmerged
         DO mi = 1, nimages_merged
            ! Check for empty images
            IF (local_refs_size(imeta, mi, ui) .EQ. 0) CYCLE
            !
            ! Insert symmetric/remote data in local buffers
            IF (do_virt .OR. ((do_symmetry .OR. transpose) .AND. .NOT. is_diagonal)) THEN
               IF (do_virt) THEN
                  IF (ithread .EQ. 0) THEN
                     refs_size(:, 0, mi, ui, ilocal_proc) = local_refs_displ(:, mi, ui)+1
                     refs_size(:, nthreads, mi, ui, ilocal_proc) = 0
                  ELSE
                     ! Reset refs displacement
                     refs_size(:, ithread, mi, ui, ilocal_proc) = 0
                     ! Accumulate over threads and procs
                     DO dst_proc = 1, nprocs_sym
                        refs_size(:, ithread, mi, ui, ilocal_proc) = &
                           refs_size(:, ithread, mi, ui, ilocal_proc)+ &
                           recv_refs(:, ithread, mi, ui, dst_proc)-recv_refs(:, ithread-1, mi, ui, dst_proc)
                     ENDDO
                  ENDIF
!$OMP BARRIER
!$OMP MASTER
                  DO it = 1, nthreads-1
                     refs_size(:, it, mi, ui, ilocal_proc) = &
                        refs_size(:, it, mi, ui, ilocal_proc)+ &
                        refs_size(:, it-1, mi, ui, ilocal_proc)
                  ENDDO
!$OMP END MASTER
!$OMP BARRIER
               ENDIF
               ! Temporary shift by index prefix
               refs_size(imeta, ithread, mi, ui, ilocal_proc) = &
                  refs_size(imeta, ithread, mi, ui, ilocal_proc)+size_index_unmerged
               !
               DO dst_proc = 1, nprocs_sym
                  IF (recv_refs(imeta, nthreads, mi, ui, dst_proc) .GT. &
                      recv_refs(imeta, 0, mi, ui, dst_proc)) THEN
                     ! Copy meta, block by block
                     IF (recv_refs(imeta, ithread+1, mi, ui, dst_proc) .GT. &
                         recv_refs(imeta, ithread, mi, ui, dst_proc)) THEN
                        DO blk = recv_refs(imeta, ithread, mi, ui, dst_proc), &
                                 recv_refs(imeta, ithread+1, mi, ui, dst_proc)-1, 3
                           buffer%meta(refs_size(imeta, ithread, mi, ui, ilocal_proc)) = &
                              meta_recv(blk)
                           buffer%meta(refs_size(imeta, ithread, mi, ui, ilocal_proc)+1) = &
                              meta_recv(blk+1)
                           buffer%meta(refs_size(imeta, ithread, mi, ui, ilocal_proc)+2) = &
                              meta_recv(blk+2)+refs_size(idata, ithread, mi, ui, ilocal_proc)- &
                              local_refs_displ(idata, mi, ui)-1
                           refs_size(imeta, ithread, mi, ui, ilocal_proc) = &
                              refs_size(imeta, ithread, mi, ui, ilocal_proc)+3
                        ENDDO
                        ! Copy data
                        CALL dbcsr_data_set(buffer%data(1), refs_size(idata, ithread, mi, ui, ilocal_proc), &
                                            recv_refs(idata, ithread+1, mi, ui, dst_proc)- &
                                            recv_refs(idata, ithread, mi, ui, dst_proc), &
                                            data_recv, recv_refs(idata, ithread, mi, ui, dst_proc))
                        refs_size(idata, ithread, mi, ui, ilocal_proc) = &
                           refs_size(idata, ithread, mi, ui, ilocal_proc)+ &
                           recv_refs(idata, ithread+1, mi, ui, dst_proc)- &
                           recv_refs(idata, ithread, mi, ui, dst_proc)
                     ENDIF
                  ENDIF
               ENDDO
               ! Remove index prefix
               refs_size(imeta, ithread, mi, ui, ilocal_proc) = &
                  refs_size(imeta, ithread, mi, ui, ilocal_proc)-size_index_unmerged
!$OMP BARRIER
            ENDIF
            !
            ! Make local indexing
!$OMP DO
            DO blk = local_refs_displ(imeta, mi, ui)+size_index_unmerged+1, &
                     local_refs_displ(imeta, mi, ui)+local_refs_size(imeta, mi, ui), 3
               meta_buffer_p(blk) = local_g2l_map_rows(meta_buffer_p(blk))
               meta_buffer_p(blk+1) = local_g2l_map_cols(meta_buffer_p(blk+1))
            ENDDO
!$OMP END DO
!$OMP MASTER
            ! Set stats slots
            meta_buffer_p(local_refs_displ(imeta, mi, ui)+dbcsr_slot_size) = &
               local_refs_size(imeta, mi, ui)
            meta_buffer_p(local_refs_displ(imeta, mi, ui)+dbcsr_slot_nblks) = &
               (local_refs_size(imeta, mi, ui)-size_index_unmerged)/3
            nblocks = nblocks+meta_buffer_p(local_refs_displ(imeta, mi, ui)+dbcsr_slot_nblks)
!$          IF (is_left) THEN
!$             meta_buffer_p(local_refs_displ(imeta, mi, ui)+dbcsr_slot_nblks+1) = 0
!$          ENDIF
!$OMP END MASTER
!$OMP BARRIER
            IF (is_left) THEN
!$             meta_buffer_p(local_refs_displ(imeta, mi, ui)+dbcsr_slot_nblks+ithread+2) = &
!$                (refs_size(imeta, ithread, mi, ui, ilocal_proc)-local_refs_displ(imeta, mi, ui)-1)/3
               ! Count the maximum possible multiplies per row for on-the-fly filtering
               IF (otf_filtering) THEN
!$OMP DO
                  DO row = local_refs_displ(imeta, mi, ui)+size_index_unmerged+1, &
                           local_refs_displ(imeta, mi, ui)+ &
                           meta_buffer_p(local_refs_displ(imeta, mi, ui)+dbcsr_slot_size), 3
!$OMP ATOMIC
                     left_total_row_counts(meta_buffer_p(row)) = &
                        left_total_row_counts(meta_buffer_p(row))+1
                  ENDDO
!$OMP END DO
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!$OMP END PARALLEL
      DEALLOCATE (refs_size, refs_displ)
      DEALLOCATE (send_size, recv_size)
      IF (do_virt) THEN
         DEALLOCATE (send_displ, recv_displ)
      ENDIF
      !
      IF (is_left .AND. otf_filtering) THEN
         CALL mp_isum(left_total_row_counts, dbcsr_mp_my_row_group(mp_obj), request_count_rows)
      ENDIF
      !
      IF (nimages_merged .EQ. 1) THEN
         CALL setup_rec_index_images(buffer%meta, img_refs_rows, img_refs_cols, &
                                     local_refs_size(imeta, 1, :), local_refs_displ(imeta, 1, :), &
                                     size_index_unmerged, is_left)
      ENDIF
      !
      IF (do_virt .OR. ((do_symmetry .OR. transpose) .AND. .NOT. is_diagonal)) THEN
         CALL dbcsr_data_release(data_recv)
         CALL memory_deallocate(meta_recv, memtype_mpi_buffer)
         DEALLOCATE (recv_refs)
      ENDIF
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
      IF (nimages_merged .GT. 1) THEN
         CALL mp_wait(request_no_empty_images)
      ELSE
         no_empty_images(:, :) = local_no_empty_images
      ENDIF
      !
      total_no_empty_images = 0
      DO iproc = 1, nprocs
         DO ilayer = 1, UBOUND(no_empty_images_displ, 1)
            no_empty_images_displ(ilayer, iproc) = total_no_empty_images
            total_no_empty_images = total_no_empty_images+ &
                                    no_empty_images(ilayer, iproc)
         ENDDO
      ENDDO
      !
      IF (is_left) THEN
         ALLOCATE (left_refs_data_size(total_no_empty_images))
         refs_data_size => left_refs_data_size
      ELSE
         ALLOCATE (right_refs_data_size(total_no_empty_images))
         refs_data_size => right_refs_data_size
      ENDIF
      IF (nimages_merged .GT. 1) THEN
         CALL mp_iallgather(local_refs_data_size(1:local_no_empty_images), &
                            refs_data_size, no_empty_images, no_empty_images_displ, &
                            buffer%subgrp, requests(4+irequests))
      ELSE
         CALL mp_iallgather(local_refs_data_size(1:local_no_empty_images), &
                            refs_data_size, &
                            buffer%subgrp, requests(4+irequests))
      ENDIF
      !
      ! Evaluate and exchange max norms for MPI filtering
      IF (otf_filtering .AND. dbcsr_cfg%use_mpi_filtering) THEN
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
                                        local_refs_size(imeta, :, :), &
                                        local_refs_displ, &
                                        img_map, &
                                        img_refs, &
                                        dbcsr_row_block_sizes(matrix), &
                                        dbcsr_col_block_sizes(matrix), &
                                        local_rows, &
                                        local_cols, &
                                        size_index_unmerged+1, &
                                        local_max_norms, is_left)
         IF (nimages_merged .GT. 1) THEN
            CALL mp_iallgather(local_max_norms, max_norms, &
                               no_empty_images, no_empty_images_displ, &
                               buffer%subgrp, requests(6+irequests))
         ELSE
            CALL mp_iallgather(local_max_norms, max_norms, &
                               buffer%subgrp, requests(6+irequests))
         ENDIF
      ENDIF
      !
      DEALLOCATE (img_refs_cols)
      DEALLOCATE (img_refs)
      DEALLOCATE (img_map)
      DEALLOCATE (local_refs_displ)
      DEALLOCATE (img_refs_rows)
      DEALLOCATE (local_refs_size)
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
      IF (dbcsr_cfg%use_mpi_exp) THEN
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
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE release_layers_3D_C_reduction()
      layers_3D_C_reduction%grp = mp_comm_null
      IF (layers_3D_C_reduction%rowgrp3D .NE. mp_comm_null) CALL mp_comm_free(layers_3D_C_reduction%rowgrp3D)
      IF (layers_3D_C_reduction%grp3D .NE. mp_comm_null) CALL mp_comm_free(layers_3D_C_reduction%grp3D)
      layers_3D_C_reduction%rowgrp3D = mp_comm_null
      layers_3D_C_reduction%grp3D = mp_comm_null
      layers_3D_C_reduction%num_layers_3D = 1
      layers_3D_C_reduction%side3D = HUGE(1)
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
         handle2, handle3, icol3D, icol3D_send, ileft_buffer_calc, ileft_buffer_comm, &
         index_row_max_epss, iproc, iright_buffer_calc, iright_buffer_comm, irow3D, irow3D_send, &
         istep_k, istep_k_ordered, ithread, ivirt_k, last_step_k, left_col_mult, left_col_nimages, &
         left_col_total_nimages, left_data_size, left_meta_size, left_myfirstvcol, &
         left_myfirstvrow, left_mypcol, left_myprow, left_npcols, left_nprows, left_nrequests, &
         left_row_mult, left_row_nimages, leftovers_first_k, leftovers_k, leftovers_shift_k
      INTEGER :: leftovers_start_k, min_nimages, mycol3D, mynode, mypcol, myprow, myrank3D, &
         myrow3D, myt, nblkrows_local, nbuffers, nbuffers_norms, ncols3D, nranks3D, nrows3D, &
         nthreads, numnodes, nvirt_k, proc3D_recv, proc3D_send, recv_vcol, recv_vrow, &
         request_epss, request_keep_sparsity, right_col_mult, right_col_nimages, right_data_size, &
         right_meta_size, right_myfirstvcol, right_myfirstvrow, right_mypcol, right_myprow, &
         right_npcols, right_nprows, right_nrequests, right_row_mult, right_row_nimages, &
         right_row_total_nimages, row, shift3D, shift3D_recv, size_guess, size_guess_init, start_k
      INTEGER :: start_k_ordered, v_ki
      INTEGER(KIND=int_8)                                :: mem
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: left_vrow, product_matrix_epss_displ, &
                                                            product_matrix_epss_size, &
                                                            product_matrix_meta, right_vcol
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: left_displ_vunmerged, &
                                                            product_matrix_meta_displ, &
                                                            product_matrix_meta_size, &
                                                            right_displ_vunmerged
      INTEGER, DIMENSION(2)                              :: requests_reduction_size
      INTEGER, DIMENSION(4)                              :: requests_reduction
      INTEGER, DIMENSION(:), POINTER :: col_blk_sizes2enum, enum2col_blk_sizes, &
         enum2row_blk_sizes, product_matrix_meta_recv, product_matrix_meta_send, row_blk_sizes2enum
      INTEGER, DIMENSION(:, :), POINTER                  :: left_refs_meta_size_layers3D, &
                                                            right_refs_meta_size_layers3D
      INTEGER, DIMENSION(:, :, :), POINTER :: left_refs_displ_unmerged_layers3D, &
         right_refs_displ_unmerged_layers3D
      INTEGER, DIMENSION(dbcsr_slot_nblkrows_total:&
         dbcsr_slot_nfullcols_local)                     :: left_global_indices, right_global_indices
      INTEGER, DIMENSION(idata:imeta)                    :: product_matrix_size_recv, &
                                                            product_matrix_size_send
      LOGICAL                                            :: do_comm, do_layers3D, &
                                                            do_square_layers3D, first_k, &
                                                            first_v_k, keep_sparsity, otf_filtering
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: do_comm_left, do_comm_right
      REAL(kind=sp)                                      :: filter_eps_sp
      REAL(kind=sp), ALLOCATABLE, DIMENSION(:), TARGET   :: row_max_epss
      REAL(kind=sp), ALLOCATABLE, DIMENSION(:, :, :)     :: left_norms, right_norms
      REAL(kind=sp), DIMENSION(:), POINTER               :: product_matrix_epss
      TYPE(dbcsr_2d_array_obj)                           :: product_matrix3D
      TYPE(dbcsr_buffer), POINTER                        :: left_buffer_p, right_buffer_p
      TYPE(dbcsr_buffer_p), ALLOCATABLE, DIMENSION(:)    :: left_buffers, right_buffers
      TYPE(dbcsr_data_obj)                               :: data_get, data_send
      TYPE(dbcsr_mm_multrec_type_p), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: multrec
      TYPE(dbcsr_mp_obj)                                 :: left_mp_obj, product_mp_obj, right_mp_obj

!   ---------------------------------------------------------------------------

      CALL timeset(routineN, handle)
      !
      NULLIFY (row_blk_sizes2enum, enum2row_blk_sizes)
      NULLIFY (col_blk_sizes2enum, enum2col_blk_sizes)
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
      IF (.NOT. buffers_orig%left%is_valid .OR. &
          .NOT. buffers_orig%right%is_valid .OR. &
          .NOT. ASSOCIATED(buffers_orig%left%meta) .OR. &
          .NOT. ASSOCIATED(buffers_orig%right%meta) .OR. &
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
          (otf_filtering .AND. dbcsr_cfg%use_mpi_filtering .AND. &
           (.NOT. (ALLOCATED(left_max_norms) .AND. ALLOCATED(right_max_norms) .AND. &
                   ALLOCATED(left_local_max_norms) .AND. ALLOCATED(right_local_max_norms))))) &
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
      mynode = dbcsr_mp_mynode(product_mp_obj)
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
      grp_right = buffers_orig%right%subgrp
      grp_left = buffers_orig%left%subgrp
      !
      do_layers3D = layers_3D_C_reduction%num_layers_3D .GT. 1
      myrow3D = myprow/layers_3D_C_reduction%side3D+1
      mycol3D = mypcol/layers_3D_C_reduction%side3D+1
      nrows3D = UBOUND(left_refs_meta_size, 2)
      ncols3D = UBOUND(right_refs_meta_size, 2)
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
      CALL win_setup(buffers_orig%left, do_win_create_left, request_win_create(1))
      CALL win_setup(buffers_orig%right, do_win_create_right, request_win_create(2))
      !
      ! Count the maximum possible multiplies per row for on-the-fly filtering
      ALLOCATE (product_matrix_epss_size(nrows3D), product_matrix_epss_displ(nrows3D))
      IF (otf_filtering) THEN
         ! Wait for counts (sent in make_buffers)
         CALL mp_wait(request_count_rows)
         !
         nblkrows_local = SIZE(left_total_row_counts)
         ALLOCATE (row_max_epss(0:nblkrows_local))
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
         row_max_epss(0) = filter_eps_sp/REAL(index_row_max_epss, KIND=KIND(row_max_epss))
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
            ALLOCATE (product_matrix_epss(0:size_guess))
            CALL mp_iallgather(row_max_epss, &
                               product_matrix_epss, product_matrix_epss_size, product_matrix_epss_displ, &
                               layers_3D_C_reduction%rowgrp3D, request_epss)
         ELSE
            product_matrix_epss_size(nrows3D) = SIZE(row_max_epss)
            product_matrix_epss_displ(nrows3D) = 0
            product_matrix_epss => row_max_epss
         ENDIF
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
      CALL mp_waitall(requests)
      !
      ! Delete buffers used for collectives
      DEALLOCATE (right_local_refs_displ_unmerged, left_local_refs_displ_unmerged)
      IF (otf_filtering .AND. dbcsr_cfg%use_mpi_filtering) THEN
         DEALLOCATE (right_local_max_norms, left_local_max_norms)
      ENDIF
      DEALLOCATE (right_local_refs_meta_size, left_local_refs_meta_size)
      DEALLOCATE (right_local_refs_data_size, left_local_refs_data_size)
      DEALLOCATE (right_no_empty_images, left_no_empty_images)
      !
      ! Needs to remap refs for 3D
      CALL remap_size_layers3D(left_col_nimages, nrows3D, left_npcols, &
                               left_refs_meta_size, &
                               left_refs_meta_size_layers3D)
      CALL remap_size_layers3D(right_row_nimages, ncols3D, right_nprows, &
                               right_refs_meta_size, &
                               right_refs_meta_size_layers3D)
      CALL remap_displ_layers3D(left_col_nimages, nrows3D, left_npcols, &
                                left_refs_displ_unmerged, &
                                left_refs_displ_unmerged_layers3D)
      CALL remap_displ_layers3D(right_row_nimages, ncols3D, right_nprows, &
                                right_refs_displ_unmerged, &
                                right_refs_displ_unmerged_layers3D)
      !
      ! Pre-execution: resize buffers for communications and take displacements
      right_data_size = 0; right_meta_size = 0
      left_data_size = 0; left_meta_size = 0
      !
      ALLOCATE (right_displ_vunmerged(ncols3D, 0:nvirt_k-1))
      ALLOCATE (left_displ_vunmerged(nrows3D, 0:nvirt_k-1))
      !
!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED (right_nprows,left_npcols,right_row_nimages,left_col_nimages,&
!$OMP         right_refs_meta_size_layers3D,left_refs_meta_size_layers3D,&
!$OMP         right_refs_data_size,left_refs_data_size,&
!$OMP         right_displ_vunmerged,left_displ_vunmerged,&
!$OMP         right_no_empty_images_displ,left_no_empty_images_displ,&
!$OMP         ncols3D,nrows3D) &
!$OMP PRIVATE (iproc,v_ki,recv_vrow,recv_vcol,icol3D,irow3D) &
!$OMP REDUCTION (MAX : right_data_size,left_data_size,right_meta_size,left_meta_size)
      !
!$OMP DO
      DO iproc = 0, right_nprows-1
         recv_vrow = iproc*right_row_nimages
         DO icol3D = 1, ncols3D
            right_displ_vunmerged(icol3D, recv_vrow) = right_no_empty_images_displ(icol3D, iproc+1)
            DO v_ki = recv_vrow, recv_vrow+right_row_nimages-2
               right_displ_vunmerged(icol3D, v_ki+1) = right_displ_vunmerged(icol3D, v_ki)+1
               right_data_size = &
                  MAX(right_data_size, right_refs_data_size(right_displ_vunmerged(icol3D, v_ki+1)))
               right_meta_size = &
                  MAX(right_meta_size, right_refs_meta_size_layers3D(icol3D, v_ki))
            ENDDO
            v_ki = recv_vrow+right_row_nimages-1
            IF (right_refs_meta_size_layers3D(icol3D, v_ki) .NE. 0) THEN
               right_data_size = &
                  MAX(right_data_size, right_refs_data_size(right_displ_vunmerged(icol3D, v_ki)+1))
               right_meta_size = &
                  MAX(right_meta_size, right_refs_meta_size_layers3D(icol3D, v_ki))
            ENDIF
         ENDDO
      ENDDO
!$OMP END DO
!$OMP DO
      DO iproc = 0, left_npcols-1
         recv_vcol = iproc*left_col_nimages
         DO irow3D = 1, nrows3D
            left_displ_vunmerged(irow3D, recv_vcol) = left_no_empty_images_displ(irow3D, iproc+1)
            DO v_ki = recv_vcol, recv_vcol+left_col_nimages-2
               left_displ_vunmerged(irow3D, v_ki+1) = left_displ_vunmerged(irow3D, v_ki)+1
               left_data_size = &
                  MAX(left_data_size, left_refs_data_size(left_displ_vunmerged(irow3D, v_ki+1)))
               left_meta_size = &
                  MAX(left_meta_size, left_refs_meta_size_layers3D(irow3D, v_ki))
            ENDDO
            v_ki = recv_vcol+left_col_nimages-1
            IF (left_refs_meta_size_layers3D(irow3D, v_ki) .NE. 0) THEN
               left_data_size = &
                  MAX(left_data_size, left_refs_data_size(left_displ_vunmerged(irow3D, v_ki)+1))
               left_meta_size = &
                  MAX(left_meta_size, left_refs_meta_size_layers3D(irow3D, v_ki))
            ENDIF
         ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      left_meta_size = left_meta_size+dbcsr_num_slots-dbcsr_slot_nblks
      right_meta_size = right_meta_size+dbcsr_num_slots-dbcsr_slot_nblks
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
      IF (has_acc) THEN
         CALL dbcsr_mempool_limit_capacity(memtype_abpanel_1%pool, &
                                           capacity=nbuffers)
         CALL dbcsr_mempool_limit_capacity(memtype_abpanel_2%pool, &
                                           capacity=nbuffers)
         CALL dbcsr_mempool_limit_capacity(memtype_trsbuffer_1%pool, &
                                           capacity=nbuffers/2+MOD(nbuffers, 2))
         CALL dbcsr_mempool_limit_capacity(memtype_trsbuffer_2%pool, &
                                           capacity=nbuffers/2+MOD(nbuffers, 2))
         ! enumerate the blocksizes to keep the following 2D-arrays small.
         CALL enumerate_blk_sizes(dbcsr_row_block_sizes(matrix_right), &
                                  row_blk_sizes2enum, enum2row_blk_sizes)
         CALL enumerate_blk_sizes(dbcsr_col_block_sizes(matrix_right), &
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
                          right_data_size, &
                          right_meta_size, &
                          num_data=(nbuffers/2), &
                          data_memory_type=memtype_abpanel_2, &
                          trs_memory_type=memtype_trsbuffer_2)
         ! Left
         CALL buffer_init(buffers_2%left, data_type, &
                          left_data_size, &
                          left_meta_size, &
                          num_data=(nbuffers/2), &
                          data_memory_type=memtype_abpanel_2)
      ENDIF
      !
      ! Prepare buffers for communication
      ! Right
      CALL buffer_init(buffers_1%right, data_type, &
                       right_data_size, &
                       right_meta_size, &
                       num_data=(nbuffers-nbuffers/2), &
                       data_memory_type=memtype_abpanel_1, &
                       trs_memory_type=memtype_trsbuffer_1)
      ! Left
      CALL buffer_init(buffers_1%left, data_type, &
                       left_data_size, &
                       left_meta_size, &
                       num_data=(nbuffers-nbuffers/2), &
                       data_memory_type=memtype_abpanel_1)
      !
      CALL setup_buffers(buffers_1%right, buffers_2%right, &
                         right_buffers, nbuffers, &
                         right_meta_size, matrix_right, &
                         imgdist_right)
      CALL setup_buffers(buffers_1%left, buffers_2%left, &
                         left_buffers, nbuffers, &
                         left_meta_size, matrix_left, &
                         imgdist_left)
      !
      ! Setup the receive data pointers
      CALL dbcsr_data_init(data_get)
      CALL dbcsr_data_new(data_get, data_type)
      IF (do_layers3D) THEN
         CALL dbcsr_data_init(data_send)
         CALL dbcsr_data_new(data_send, data_type)
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
      IF (.NOT. keep_sparsity) THEN
         size_guess_init = product_matrix_size_guess(matrix_left, matrix_right, product_matrix, &
                                                     left_data_size, right_data_size, &
                                                     left_col_nimages, right_row_nimages, &
                                                     nthreads)
      ENDIF
      !
      ! Preallocate norms arrays
      IF (otf_filtering) THEN
         ALLOCATE (right_norms(right_meta_size/3, 1, nbuffers_norms))
         ALLOCATE (left_norms(left_meta_size/3, 1, nbuffers_norms))
         IF (do_layers3D .AND. nrows3D .GT. 1) THEN
            CALL mp_wait(request_epss)
            DEALLOCATE (row_max_epss)
         ENDIF
      ELSE
         ! The array must be valid when passed to called subroutines.
         ALLOCATE (right_norms(0, 1, nbuffers_norms))
         ALLOCATE (left_norms(0, 1, nbuffers_norms))
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
      CALL timeset(routineN//"_loop", handle2)
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
      irow3D_send = 0
      icol3D_send = 0
      first_k = .TRUE.
      first_v_k = .TRUE.
      grouped_steps_index: DO istep_k = leftovers_first_k, last_step_k
         !
         ! Matrix transfer. Transfer in all but the last loop
         ! iteration.
         xfer: IF (istep_k .LT. last_step_k) THEN
            start_k_ordered = start_k
            istep_k_ordered = istep_k
            ! Put leftovers ticks always first
            IF (leftovers_k .GT. 0) THEN
               IF (istep_k .LT. leftovers_first_k+leftovers_k) THEN
                  start_k_ordered = leftovers_start_k
               ELSE
                  istep_k_ordered = istep_k+leftovers_shift_k
               ENDIF
            ENDIF
            first_k = MOD(istep_k_ordered, nranks3D) .EQ. 0
            ivirt_k = istep_k_ordered/nranks3D
            IF (istep_k .LT. leftovers_first_k+leftovers_k) THEN
               CALL row_col_3D_reflected(irow3D, icol3D, nrows3D, ncols3D, istep_k_ordered)
            ELSE
               CALL row_col_3D_reflected(irow3D, icol3D, nrows3D, ncols3D, shift3D)
               shift3D = shift3D+1
            ENDIF
            !
            v_ki = MOD(ivirt_k, min_nimages)
            ! Reset communication flags at the first layer
            IF (first_k .OR. istep_k .EQ. leftovers_first_k) THEN
               do_comm_right(:) = .TRUE.
               do_comm_left(:) = .TRUE.
            ENDIF
            ! Take first image global virtual coordinates
            IF (v_ki .EQ. 0) THEN
               IF (istep_k .GE. leftovers_first_k+leftovers_k) first_v_k = .FALSE.
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
            right_buffer_p => right_buffers(iright_buffer_comm)%b
            left_buffer_p => left_buffers(ileft_buffer_comm)%b
            right_buffer_p%coord3D = icol3D
            left_buffer_p%coord3D = irow3D
            !
            ! First row, communicate right matrix
            IF (do_comm_right(icol3D)) THEN
               IF (has_acc) THEN
                  CALL timeset(routineN//"_acc_sync", handle3)
                  CALL acc_event_synchronize(right_buffer_p%data(1)%d%acc_ready)
                  CALL timestop(handle3)
               ENDIF
               right_buffer_p%vprow = MOD(recv_vrow+v_ki, right_row_total_nimages)
               right_buffer_p%vpcol = right_vcol(icol3D)
               right_buffer_p%nrequests = 0
            ENDIF
            !
            IF (right_refs_meta_size_layers3D(icol3D, right_buffer_p%vprow) .NE. 0) THEN
               ! First col, communicate left matrix
               IF (do_comm_left(irow3D)) THEN
                  IF (has_acc) THEN
                     CALL timeset(routineN//"_acc_sync", handle3)
                     CALL acc_event_synchronize(left_buffer_p%data(1)%d%acc_ready)
                     CALL timestop(handle3)
                  ENDIF
                  !
                  left_buffer_p%vprow = left_vrow(irow3D)
                  left_buffer_p%vpcol = MOD(recv_vcol+v_ki, left_col_total_nimages)
                  left_buffer_p%nrequests = 0
               ENDIF
               !
               IF (left_refs_meta_size_layers3D(irow3D, left_buffer_p%vpcol) .NE. 0) THEN
                  do_comm = .TRUE.
                  IF (otf_filtering .AND. dbcsr_cfg%use_mpi_filtering) THEN
                     IF (left_max_norms(left_displ_vunmerged(irow3D, left_buffer_p%vpcol)+1)* &
                         right_max_norms(right_displ_vunmerged(icol3D, &
                                                               right_buffer_p%vprow)+1) .LT. &
                         product_matrix_epss(product_matrix_epss_displ(irow3D))) THEN
                        dbcsr_mpi_statistics%nfiltered = dbcsr_mpi_statistics%nfiltered+1
                        do_comm = .FALSE.
                     ENDIF
                  ENDIF
                  !
                  IF (do_comm) THEN
                     ! Right
                     IF (do_comm_right(icol3D)) THEN
                        do_comm_right(icol3D) = .FALSE.
                        CALL rma_transfer(right_buffer_p%vprow, &
                                          right_row_nimages, &
                                          right_displ_vunmerged, &
                                          right_refs_displ_unmerged_layers3D, &
                                          right_buffer_p, &
                                          right_refs_meta_size_layers3D, right_refs_data_size, &
                                          buffers_orig%right%meta_win, buffers_orig%right%data_win, &
                                          data_get, data_type_byte, buffers_orig%right, icol3D)
                        ! Set the referenced sizes to the actual data moved via MPI
                        CALL dbcsr_data_set_size_referenced(right_buffer_p%data(1), &
                                                            right_buffer_p%offset(idata, 2))
                        right_buffer_p%buffer%mats(1)%valid = .FALSE.
                     ENDIF
                     ! Left
                     IF (do_comm_left(irow3D)) THEN
                        do_comm_left(irow3D) = .FALSE.
                        CALL rma_transfer(left_buffer_p%vpcol, &
                                          left_col_nimages, &
                                          left_displ_vunmerged, &
                                          left_refs_displ_unmerged_layers3D, &
                                          left_buffer_p, &
                                          left_refs_meta_size_layers3D, left_refs_data_size, &
                                          buffers_orig%left%meta_win, buffers_orig%left%data_win, &
                                          data_get, data_type_byte, buffers_orig%left, irow3D)
                        ! Set the referenced sizes to the actual data moved via MPI
                        CALL dbcsr_data_set_size_referenced(left_buffer_p%data(1), &
                                                            left_buffer_p%offset(idata, 2))
                        left_buffer_p%buffer%mats(1)%valid = .FALSE.
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF xfer
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
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP          PRIVATE (size_guess, ithread) &
!$OMP          SHARED (product_matrix3D, multrec, &
!$OMP                  keep_sparsity, filter_eps, &
!$OMP                  product_matrix_epss, &
!$OMP                  matrix_right, nthreads, &
!$OMP                  irow3D, icol3D, myrow3D, mycol3D, keep_product_data, &
!$OMP                  product_matrix_epss_displ, product_matrix_epss_size, &
!$OMP                  memtype_product_wm, size_guess_init, nranks3D)
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
            size_guess = product_matrix3D%mats(irow3D, icol3D)% &
                         matrix%wms(ithread+1)%datasize ! Should be minimal
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
                                                                        product_matrix_epss_size(irow3D)-1), &
                                       block_estimate=1, &
                                       right_row_blk_size=dbcsr_row_block_sizes(matrix_right), &
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
         ! Wait data and do the multiplications.
         ! Exclude the first interation
         wait_calc: IF (istep_k .GT. leftovers_first_k) THEN
            IF (debug_mod) WRITE (*, '(1X,A)') routineN//" waiting for right and left"
            !
            right_buffer_p => right_buffers(iright_buffer_calc)%b
            left_buffer_p => left_buffers(ileft_buffer_calc)%b
            irow3D = left_buffer_p%coord3D
            icol3D = right_buffer_p%coord3D
            IF (istep_k .GT. final_step_k) THEN
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (multrec, irow3D, icol3D, irow3D_send, icol3D_send, &
!$OMP         istep_k, final_step_k, product_matrix3D, &
!$OMP         handle3, requests_reduction_size, &
!$OMP         product_matrix_meta_send, product_matrix_meta_recv, &
!$OMP         product_matrix_size_send, product_matrix_size_recv, &
!$OMP         buffers_orig, memtype_mpi_buffer, &
!$OMP         data_send, data_get, proc3D_send, proc3D_recv, &
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
                     buffers_orig%left%meta_red3D, &
                     buffers_orig%left%data_red3D)
                  DEALLOCATE (multrec(ithread, irow3D_send, icol3D_send)%p)
                  CALL dbcsr_work_destroy( &
                     product_matrix3D%mats(irow3D_send, icol3D_send)%matrix%wms(ithread+1))
!$OMP BARRIER
!$OMP MASTER
                  DEALLOCATE (product_matrix3D%mats(irow3D_send, icol3D_send)%matrix%wms)
                  CALL timeset(routineN//"_red3D_size", handle3)
                  CALL mp_waitall(requests_reduction_size)
                  CALL timestop(handle3)
                  CALL ensure_array_size(buffers_orig%right%meta_red3D, &
                                         ub=product_matrix_size_recv(imeta), &
                                         nocopy=.TRUE., memory_type=memtype_mpi_buffer)
                  product_matrix_meta_send => &
                     buffers_orig%left%meta_red3D(1:product_matrix_size_send(imeta))
                  product_matrix_meta_recv => &
                     buffers_orig%right%meta_red3D(1:product_matrix_size_recv(imeta))
                  CALL mp_isendrecv(product_matrix_meta_send, proc3D_send, &
                                    product_matrix_meta_recv, proc3D_recv, &
                                    layers_3D_C_reduction%grp3D, &
                                    requests_reduction(1), requests_reduction(2))
                  CALL dbcsr_data_ensure_size(buffers_orig%right%data_red3D, &
                                              product_matrix_size_recv(idata), &
                                              nocopy=.TRUE.)
                  CALL dbcsr_data_set_pointer( &
                     area=data_send, &
                     rsize=product_matrix_size_send(idata), &
                     csize=1, &
                     pointee=buffers_orig%left%data_red3D)
                  CALL dbcsr_data_set_pointer( &
                     area=data_get, &
                     rsize=product_matrix_size_recv(idata), &
                     csize=1, &
                     pointee=buffers_orig%right%data_red3D)
                  CALL dbcsr_isendrecv_any(data_send, proc3D_send, &
                                           data_get, proc3D_recv, &
                                           layers_3D_C_reduction%grp3D, &
                                           requests_reduction(3), requests_reduction(4))
                  CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(1, :), &
                                            product_matrix_size_send(idata), &
                                            dbcsr_mpi_statistics%data_size_breakdown(:, :, 1), data_type_byte)
!$OMP END MASTER
               ENDIF
!$OMP END PARALLEL
            ENDIF
            !
            right_nrequests = right_buffer_p%nrequests
            left_nrequests = left_buffer_p%nrequests
            !
            IF (right_nrequests .GT. 0 .AND. left_nrequests .GT. 0) THEN
               ! check if right matrix was already initialized
               IF (.NOT. right_buffer_p%buffer%mats(1)%valid) THEN
                  CALL mp_waitall(right_buffer_p%get_requests(:, 1:right_nrequests))
                  IF (has_acc) CALL dbcsr_data_host2dev(right_buffer_p%data(1))
                  ! Repoint indices of matrices
                  CALL merge_images(right_buffer_p%buffer, min_nimages, 1, &
                                    right_row_total_nimages, &
                                    right_buffer_p%vprow, &
                                    right_buffer_p%vpcol, &
                                    right_buffer_p%offset, &
                                    right_buffer_p%meta, &
                                    imgdist=imgdist_right, do_merge_rows=.FALSE., &
                                    global_indices=right_global_indices)
                  IF (otf_filtering) THEN
                     CALL calculate_image_norms(right_buffer_p%buffer, &
                                                data_type, &
                                                right_norms(:, :, MIN(iright_buffer_calc, nbuffers_norms)), &
                                                uf=1, ul=1)
                  ENDIF
                  IF (has_acc) THEN
                     CALL acc_transpose_blocks_images(right_buffer_p%buffer, &
                                                      right_buffer_p%data(1), &
                                                      right_buffer_p%trs_stackbuf(1), &
                                                      row_blk_sizes2enum, enum2row_blk_sizes, &
                                                      col_blk_sizes2enum, enum2col_blk_sizes, &
                                                      (/1/))
                  ENDIF
               ENDIF
               ! check if left matrix was already initialized
               IF (.NOT. left_buffer_p%buffer%mats(1)%valid) THEN
                  CALL mp_waitall(left_buffer_p%get_requests(:, 1:left_nrequests))
                  IF (has_acc) CALL dbcsr_data_host2dev(left_buffer_p%data(1))
                  ! Repoint indices of matrices
                  CALL merge_images(left_buffer_p%buffer, 1, min_nimages, &
                                    left_col_total_nimages, &
                                    left_buffer_p%vprow, &
                                    left_buffer_p%vpcol, &
                                    left_buffer_p%offset, &
                                    left_buffer_p%meta, &
                                    imgdist=imgdist_left, do_merge_rows=.TRUE., &
                                    global_indices=left_global_indices, &
                                    nthreads=nthreads)
                  IF (otf_filtering) THEN
                     CALL calculate_image_norms(left_buffer_p%buffer, &
                                                data_type, &
                                                left_norms(:, :, MIN(ileft_buffer_calc, nbuffers_norms)), &
                                                uf=1, ul=1)
                  ENDIF
               ENDIF
               !
               CALL timeset(routineN//"_multrec", handle3)
               !
!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (left_buffer_p, ileft_buffer_calc, &
!$OMP         right_buffer_p, iright_buffer_calc, &
!$OMP         left_norms,right_norms, nbuffers_norms, &
!$OMP         multrec, irow3D, icol3D) &
!$OMP PRIVATE (ithread) &
!$OMP REDUCTION (+: flop)
               ithread = 0
!$             ithread = omp_get_thread_num()
               CALL dbcsr_mm_multrec_multiply(multrec(ithread, irow3D, icol3D)%p, &
                                              left=left_buffer_p%buffer%mats(1), &
                                              right=right_buffer_p%buffer%mats(1), &
                                              flop=flop, &
                                              a_norms=left_norms(:, 1, MIN(ileft_buffer_calc, nbuffers_norms)), &
                                              b_norms=right_norms(:, 1, MIN(iright_buffer_calc, nbuffers_norms)))
!$OMP END PARALLEL
               !
               CALL timestop(handle3)
            ENDIF
            ! Reduce 3D layers and finalize the local layer
            IF (istep_k .GT. final_step_k) THEN
               ! Wait for the other 3D layers to reduce
               IF (istep_k .GT. final_step_k+1) THEN
                  CALL timeset(routineN//"_red3D_data", handle3)
                  CALL mp_waitall(requests_reduction)
                  CALL timestop(handle3)
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
!$OMP         final_step_k, handle3, buffers_orig, g2l_map_rows, g2l_map_cols) &
!$OMP PRIVATE (ithread) &
!$OMP REDUCTION (+: flop)
               ithread = 0
!$             ithread = omp_get_thread_num()
               CALL dbcsr_mm_multrec_dev2host_init(multrec(ithread, irow3D, icol3D)%p)
               !
               ! Evaluate the size of layers to send and set the buffers
               IF (irow3D .NE. myrow3D .OR. &
                   icol3D .NE. mycol3D) THEN
!$OMP ATOMIC
                  product_matrix3D%mats(irow3D_send, icol3D_send)%matrix%nblks = &
                     product_matrix3D%mats(irow3D_send, icol3D_send)%matrix%nblks+ &
                     dbcsr_mm_multrec_get_nblks(multrec(ithread, irow3D_send, icol3D_send)%p)
!$OMP ATOMIC
                  product_matrix3D%mats(irow3D_send, icol3D_send)%matrix%nze = &
                     product_matrix3D%mats(irow3D_send, icol3D_send)%matrix%nze+ &
                     dbcsr_mm_multrec_get_nze(multrec(ithread, irow3D_send, icol3D_send)%p)
!$OMP BARRIER
!$OMP MASTER
                  ! First (nthreads+1)*2 positions are reserved for
                  ! the offset sizes of each thread for meta and data
                  CALL ensure_array_size(buffers_orig%left%meta_red3D, &
                                         ub=product_matrix3D%mats(irow3D_send, icol3D_send)% &
                                         matrix%nblks*3+(nthreads+1)*2, &
                                         nocopy=.TRUE., memory_type=memtype_mpi_buffer)
                  CALL dbcsr_data_ensure_size(buffers_orig%left%data_red3D, &
                                              product_matrix3D%mats(irow3D_send, icol3D_send)%matrix%nze, &
                                              nocopy=.TRUE.)
                  ! Set the offsets
                  buffers_orig%left%meta_red3D(1) = (nthreads+1)*2
                  buffers_orig%left%meta_red3D(nthreads+2) = 0
                  DO myt = 1, nthreads
                     buffers_orig%left%meta_red3D(myt+1) = &
                        buffers_orig%left%meta_red3D(myt)+ &
                        dbcsr_mm_multrec_get_nblks(multrec(myt-1, irow3D_send, icol3D_send)%p)*3
                     buffers_orig%left%meta_red3D(myt+nthreads+2) = &
                        buffers_orig%left%meta_red3D(myt+nthreads+1)+ &
                        dbcsr_mm_multrec_get_nze(multrec(myt-1, irow3D_send, icol3D_send)%p)
                  ENDDO
                  ! Send/recv data and meta sizes
                  product_matrix_size_send(idata) = &
                     buffers_orig%left%meta_red3D((nthreads+1)*2)
                  product_matrix_size_send(imeta) = &
                     buffers_orig%left%meta_red3D(nthreads+1)
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
                  IF (dbcsr_data_get_size_referenced(buffers_orig%right%data_red3D) .GT. 0) THEN
                     CALL timeset(routineN//"_red3D", handle3)
                     CALL dbcsr_mm_multrec_red3D(multrec(ithread, myrow3D, mycol3D)%p, &
                                                 buffers_orig%right%meta_red3D, &
                                                 buffers_orig%right%data_red3D, flop, &
                                                 g2l_map_rows, g2l_map_cols)
                     CALL timestop(handle3)
                  ENDIF
               ENDIF
!$OMP END PARALLEL
            ENDIF
         ENDIF wait_calc
         !
         ! Swap temporary buffers for the first tick
         IF (do_square_layers3D .AND. first_k .AND. &
             istep_k .LT. last_step_k) THEN
            iright_buffer_comm = right_buffers(iright_buffer_comm)%b%coord3D
            ileft_buffer_comm = left_buffers(ileft_buffer_comm)%b%coord3D
            CALL swap_buffers(right_buffers(iright_buffer_comm), right_buffers(nbuffers))
            CALL swap_buffers(left_buffers(ileft_buffer_comm), left_buffers(nbuffers))
         ENDIF
         !
         iright_buffer_calc = iright_buffer_comm
         ileft_buffer_calc = ileft_buffer_comm
      ENDDO grouped_steps_index
      !
      CALL timestop(handle2)
      !
      CALL m_memory(mem)
      max_memory = MAX(max_memory, REAL(mem))
      !
      IF (do_layers3D .AND. keep_sparsity) THEN
         DEALLOCATE (product_matrix_meta_size, product_matrix_meta_displ)
         DEALLOCATE (product_matrix_meta)
      ENDIF
      DEALLOCATE (right_norms, left_norms)
      IF (otf_filtering .AND. dbcsr_cfg%use_mpi_filtering) THEN
         DEALLOCATE (right_max_norms, left_max_norms)
      ENDIF
      DEALLOCATE (product_matrix_epss_size, product_matrix_epss_displ)
      IF (.NOT. otf_filtering .OR. (do_layers3D .AND. nrows3D .GT. 1)) THEN
         DEALLOCATE (product_matrix_epss)
      ELSE
         DEALLOCATE (row_max_epss)
      ENDIF
      !
      DEALLOCATE (left_refs_meta_size, left_refs_displ_unmerged)
      DEALLOCATE (right_refs_meta_size, right_refs_displ_unmerged)
      IF (left_col_nimages .GT. 1) &
         DEALLOCATE (left_refs_meta_size_layers3D, left_refs_displ_unmerged_layers3D)
      IF (right_row_nimages .GT. 1) &
         DEALLOCATE (right_refs_meta_size_layers3D, right_refs_displ_unmerged_layers3D)
      DEALLOCATE (right_refs_data_size, left_refs_data_size)
      DEALLOCATE (right_displ_vunmerged, left_displ_vunmerged)
      DEALLOCATE (right_no_empty_images_displ, left_no_empty_images_displ)
      !
      ! Deallocate 3D layers
      IF (do_layers3D) THEN
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
      CALL dbcsr_finalize(product_matrix)
      !
      CALL dbcsr_data_clear_pointer(data_get)
      CALL dbcsr_data_release(data_get)
      !
      ! clean-up of communication buffers
      DO ileft_buffer_comm = 1, nbuffers
         DEALLOCATE (left_buffers(ileft_buffer_comm)%b%data)
         CALL dbcsr_destroy_array(left_buffers(ileft_buffer_comm)%b%buffer)
         DEALLOCATE (left_buffers(ileft_buffer_comm)%b%offset)
         DEALLOCATE (left_buffers(ileft_buffer_comm)%b%get_requests)
         DEALLOCATE (left_buffers(ileft_buffer_comm)%b)
      ENDDO
      DO iright_buffer_comm = 1, nbuffers
         DEALLOCATE (right_buffers(iright_buffer_comm)%b%data)
         CALL dbcsr_destroy_array(right_buffers(iright_buffer_comm)%b%buffer)
         IF (has_acc) THEN
            DEALLOCATE (right_buffers(iright_buffer_comm)%b%trs_stackbuf)
         ENDIF
         DEALLOCATE (right_buffers(iright_buffer_comm)%b%offset)
         DEALLOCATE (right_buffers(iright_buffer_comm)%b%get_requests)
         DEALLOCATE (right_buffers(iright_buffer_comm)%b)
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

      INTEGER                                            :: handle

      IF (buffer%has_rma_win) THEN
         CALL timeset(routineN//"_win_check", handle)
         CALL mp_wait(request)
         CALL timestop(handle)
         IF (do_win_create(1)) THEN
            CALL mp_win_unlock_all(buffer%data_win)
            CALL mp_win_free(buffer%data_win)
         ENDIF
         IF (do_win_create(2)) THEN
            CALL mp_win_unlock_all(buffer%meta_win)
            CALL mp_win_free(buffer%meta_win)
         ENDIF
      ENDIF
      CALL dbcsr_data_release(buffer%data_resize)
      IF (ASSOCIATED(buffer%meta_resize)) THEN
         CALL memory_deallocate(buffer%meta_resize, memtype_mpi_buffer)
         NULLIFY (buffer%meta_resize)
      ENDIF
      !
      IF (do_win_create(1)) THEN
         CALL dbcsr_win_create_any(buffer%data(1), buffer%subgrp, buffer%data_win)
         CALL mp_win_lock_all(buffer%data_win)
      ENDIF
      IF (do_win_create(2)) THEN
         CALL mp_win_create(buffer%meta, buffer%subgrp, buffer%meta_win)
         CALL mp_win_lock_all(buffer%meta_win)
      ENDIF
      !
      buffer%has_rma_win = .TRUE.
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
!> \param meta_size ...
!> \param matrix ...
!> \param imgdist ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE setup_buffers(buffer_1, buffer_2, buffers, nbuffers, meta_size, matrix, imgdist)
      TYPE(dbcsr_buffer), INTENT(INOUT)                  :: buffer_1, buffer_2
      TYPE(dbcsr_buffer_p), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: buffers
      INTEGER, INTENT(IN)                                :: nbuffers, meta_size
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(dbcsr_imagedistribution_obj), INTENT(INOUT)   :: imgdist

      INTEGER                                            :: ibuffer, jbuffer
      LOGICAL                                            :: has_trs_stackbuff

      ALLOCATE (buffers(nbuffers))
      has_trs_stackbuff = ALLOCATED(buffer_1%trs_stackbuf) .OR. ALLOCATED(buffer_2%trs_stackbuf)
      DO ibuffer = 1, nbuffers
         ALLOCATE (buffers(ibuffer)%b)
         ALLOCATE (buffers(ibuffer)%b%data(1))
         IF (has_trs_stackbuff) THEN
            ALLOCATE (buffers(ibuffer)%b%trs_stackbuf(1))
         ENDIF
         jbuffer = (ibuffer-1)/2
         IF (MOD(ibuffer, 2) .EQ. 1) THEN
            buffers(ibuffer)%b%data(1) = buffer_1%data(jbuffer+1)
            buffers(ibuffer)%b%meta => &
               buffer_1%meta(meta_size*jbuffer+1: &
                             meta_size*(jbuffer+1))
            IF (has_trs_stackbuff) THEN
               buffers(ibuffer)%b%trs_stackbuf(1) = buffer_1%trs_stackbuf(jbuffer+1)
            ENDIF
         ELSE
            buffers(ibuffer)%b%data(1) = buffer_2%data(jbuffer+1)
            buffers(ibuffer)%b%meta => &
               buffer_2%meta(meta_size*jbuffer+1: &
                             meta_size*(jbuffer+1))
            IF (has_trs_stackbuff) THEN
               buffers(ibuffer)%b%trs_stackbuf(1) = buffer_2%trs_stackbuf(jbuffer+1)
            ENDIF
         ENDIF
         ALLOCATE (buffers(ibuffer)%b%offset(2, 2))
         buffers(ibuffer)%b%offset(:, 1) = 0
         ALLOCATE (buffers(ibuffer)%b%get_requests(2, 1))
         CALL setup_buffer_matrices_images(buffers(ibuffer)%b%buffer, 1, &
                                           imgdist, matrix, &
                                           buffers(ibuffer)%b%data(1), &
                                           index_size=meta_size)
      ENDDO
   END SUBROUTINE setup_buffers

! **************************************************************************************************
!> \brief ...
!> \param buffers_1 ...
!> \param buffers_2 ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE swap_buffers(buffers_1, buffers_2)
      TYPE(dbcsr_buffer_p), INTENT(INOUT)                :: buffers_1, buffers_2

      TYPE(dbcsr_buffer_p)                               :: tmp

      tmp = buffers_1
      buffers_1 = buffers_2
      buffers_2 = tmp
   END SUBROUTINE swap_buffers

! **************************************************************************************************
!> \brief ...
!> \param recv_vunmerged ...
!> \param unmerged_nimages ...
!> \param displ_vunmerged ...
!> \param refs_displ_unmerged ...
!> \param buffer ...
!> \param refs_meta_size ...
!> \param refs_data_size ...
!> \param meta_win ...
!> \param data_win ...
!> \param data_get ...
!> \param data_type_byte ...
!> \param buffer_orig ...
!> \param layer3D ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE rma_transfer(recv_vunmerged, &
                           unmerged_nimages, &
                           displ_vunmerged, &
                           refs_displ_unmerged, &
                           buffer, refs_meta_size, refs_data_size, &
                           meta_win, data_win, &
                           data_get, data_type_byte, &
                           buffer_orig, layer3D)
      INTEGER, INTENT(IN)                                :: recv_vunmerged, unmerged_nimages
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(IN)  :: displ_vunmerged
      INTEGER, DIMENSION(:, :, :), INTENT(IN), POINTER   :: refs_displ_unmerged
      TYPE(dbcsr_buffer), INTENT(INOUT)                  :: buffer
      INTEGER, DIMENSION(:, :), INTENT(IN), POINTER      :: refs_meta_size
      INTEGER, DIMENSION(:), INTENT(IN)                  :: refs_data_size
      INTEGER, INTENT(IN)                                :: meta_win, data_win
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: data_get
      INTEGER, INTENT(IN)                                :: data_type_byte
      TYPE(dbcsr_buffer), INTENT(IN)                     :: buffer_orig
      INTEGER, INTENT(IN)                                :: layer3D

      INTEGER                                            :: recv_punmerged
      INTEGER, DIMENSION(:), POINTER                     :: meta_get

! Left  : merged => row, unmerged => col
! Right : merged => col, unmerged => row

      buffer%nrequests = 1
      buffer%get_requests(:, :) = mp_request_null
      !
      recv_punmerged = (recv_vunmerged/unmerged_nimages)*SIZE(refs_meta_size, 1)+layer3D-1
      buffer%offset(idata, 2) = refs_data_size(displ_vunmerged(layer3D, recv_vunmerged)+1)
      buffer%offset(imeta, 2) = refs_meta_size(layer3D, recv_vunmerged)
      !
      meta_get => buffer%meta(1:refs_meta_size(layer3D, recv_vunmerged))
      CALL mp_rget(meta_get, recv_punmerged, &
                   meta_win, &
                   buffer_orig%meta, &
                   disp=refs_displ_unmerged(imeta, layer3D, recv_vunmerged), &
                   request=buffer%get_requests(1, buffer%nrequests))
      CALL dbcsr_data_set_pointer( &
         area=data_get, &
         rsize=refs_data_size(displ_vunmerged(layer3D, recv_vunmerged)+1), &
         csize=1, &
         pointee=buffer%data(1), &
         source_lb=1)
      CALL dbcsr_rget_any(data_get, recv_punmerged, &
                          data_win, &
                          buffer_orig%data(1), &
                          disp=refs_displ_unmerged(idata, layer3D, recv_vunmerged), &
                          request=buffer%get_requests(2, buffer%nrequests))
      CALL count_mpi_statistics(dbcsr_mpi_statistics%data_size(1, :), &
                                refs_data_size(displ_vunmerged(layer3D, recv_vunmerged)+1), &
                                dbcsr_mpi_statistics%data_size_breakdown(:, :, 1), data_type_byte)
      dbcsr_mpi_statistics%nexchanged = dbcsr_mpi_statistics%nexchanged+1
   END SUBROUTINE rma_transfer

! **************************************************************************************************
!> \brief ...
!> \param buffer_set ...
!> \param nimages ...
!> \param imgdist ...
!> \param template_matrix ...
!> \param data_buffer ...
!> \param index_size_images ...
!> \param index_size ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE setup_buffer_matrices_images(buffer_set, nimages, imgdist, &
                                           template_matrix, data_buffer, &
                                           index_size_images, index_size)
      TYPE(dbcsr_1d_array_type), INTENT(INOUT)           :: buffer_set
      INTEGER, INTENT(IN)                                :: nimages
      TYPE(dbcsr_imagedistribution_obj), INTENT(INOUT)   :: imgdist
      TYPE(dbcsr_type), INTENT(IN)                       :: template_matrix
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: data_buffer
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: index_size_images
      INTEGER, INTENT(IN), OPTIONAL                      :: index_size

      INTEGER                                            :: v_ki

!   ---------------------------------------------------------------------------

      buffer_set%image_dist = imgdist
      CALL dbcsr_image_dist_hold(imgdist)
      !
      ALLOCATE (buffer_set%mats(nimages))
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(nimages,buffer_set,template_matrix,&
!$OMP        imgdist,data_buffer,memtype_mpi_buffer,&
!$OMP        index_size_images,index_size)
      DO v_ki = 1, nimages
         buffer_set%mats(v_ki) = dbcsr_type()
         CALL dbcsr_create(buffer_set%mats(v_ki), &
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
         IF (PRESENT(index_size_images)) THEN
            CALL ensure_array_size(buffer_set%mats(v_ki)%index, &
                                   ub=index_size_images(v_ki), nocopy=.TRUE., &
                                   memory_type=dbcsr_get_index_memory_type(buffer_set%mats(v_ki)))
         ELSEIF (PRESENT(index_size)) THEN
            CALL ensure_array_size(buffer_set%mats(v_ki)%index, &
                                   ub=index_size, nocopy=.TRUE., &
                                   memory_type=dbcsr_get_index_memory_type(buffer_set%mats(v_ki)))
         ENDIF
         buffer_set%mats(v_ki)%negate_real = template_matrix%negate_real
         buffer_set%mats(v_ki)%negate_imaginary = template_matrix%negate_imaginary
         buffer_set%mats(v_ki)%local_indexing = .TRUE.
         buffer_set%mats(v_ki)%list_indexing = .TRUE.
      ENDDO
!$OMP END PARALLEL DO
   END SUBROUTINE setup_buffer_matrices_images

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
!$OMP    PARALLEL IF (has_threads) DEFAULT (NONE) &
!$OMP    PRIVATE (ithread) &
!$OMP    FIRSTPRIVATE (t_f, t_l, t_size) &
!$OMP    SHARED (meta_buffer, in, has_threads, refs_displ, &
!$OMP            size_index_unmerged, nblkrows_local, nblkcols_local)
!$       ithread = OMP_GET_THREAD_NUM()+dbcsr_slot_nblks+1
!$       IF (has_threads) THEN
!$          t_f = meta_buffer(refs_displ(in)+ithread)+1
!$          t_l = meta_buffer(refs_displ(in)+ithread+1)
!$       ENDIF
         t_size = t_l-t_f+1
         IF (t_size .GT. 0) THEN
            CALL rec_sort_index(1, nblkrows_local, &
                                1, nblkcols_local, &
                                t_size, &
                                meta_buffer(refs_displ(in)+size_index_unmerged+t_f*3-2: &
                                            refs_displ(in)+size_index_unmerged+t_l*3), &
                                0)
         ENDIF
!$OMP    END PARALLEL
      ENDDO
      CALL timestop(handle)
   END SUBROUTINE setup_rec_index_images

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
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE calculate_max_image_norms(buffer, refs_size, refs_displ, &
                                        img_map, img_offset, &
                                        row_blk_size, col_blk_size, &
                                        local_rows, local_cols, &
                                        slot_coo_l, &
                                        max_norms, is_left)
      TYPE(dbcsr_buffer), INTENT(IN)                     :: buffer
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: refs_size
      INTEGER, DIMENSION(:, :, :), INTENT(IN), POINTER   :: refs_displ
      INTEGER, DIMENSION(:), INTENT(IN)                  :: img_map, img_offset, row_blk_size, &
                                                            col_blk_size, local_rows, local_cols
      INTEGER, INTENT(IN)                                :: slot_coo_l
      REAL(kind=sp), DIMENSION(:), INTENT(INOUT)         :: max_norms
      LOGICAL                                            :: is_left

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_max_image_norms', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      !
      SELECT CASE (buffer%data (1)%d%data_type)
      CASE (dbcsr_type_real_4)
         CALL calc_max_image_norms_s(buffer%meta, dbcsr_get_data_p_s(buffer%data(1)), &
                                     refs_size, refs_displ, &
                                     img_map, &
                                     img_offset, &
                                     row_blk_size, col_blk_size, &
                                     local_rows, local_cols, &
                                     slot_coo_l, &
                                     max_norms, is_left)
      CASE (dbcsr_type_real_8)
         CALL calc_max_image_norms_d(buffer%meta, dbcsr_get_data_p_d(buffer%data(1)), &
                                     refs_size, refs_displ, &
                                     img_map, &
                                     img_offset, &
                                     row_blk_size, col_blk_size, &
                                     local_rows, local_cols, &
                                     slot_coo_l, &
                                     max_norms, is_left)
      CASE (dbcsr_type_complex_4)
         CALL calc_max_image_norms_c(buffer%meta, dbcsr_get_data_p_c(buffer%data(1)), &
                                     refs_size, refs_displ, &
                                     img_map, &
                                     img_offset, &
                                     row_blk_size, col_blk_size, &
                                     local_rows, local_cols, &
                                     slot_coo_l, &
                                     max_norms, is_left)
      CASE (dbcsr_type_complex_8)
         CALL calc_max_image_norms_z(buffer%meta, dbcsr_get_data_p_z(buffer%data(1)), &
                                     refs_size, refs_displ, &
                                     img_map, &
                                     img_offset, &
                                     row_blk_size, col_blk_size, &
                                     local_rows, local_cols, &
                                     slot_coo_l, &
                                     max_norms, is_left)
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
!> \param matrices ...
!> \param DATA ...
!> \param trs_stackbuf ...
!> \param row_blk_sizes2enum ...
!> \param enum2row_blk_sizes ...
!> \param col_blk_sizes2enum ...
!> \param enum2col_blk_sizes ...
!> \param indices ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE acc_transpose_blocks_images(matrices, DATA, trs_stackbuf, &
                                          row_blk_sizes2enum, enum2row_blk_sizes, &
                                          col_blk_sizes2enum, enum2col_blk_sizes, &
                                          indices)
      TYPE(dbcsr_1d_array_type), INTENT(IN)              :: matrices
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: DATA, trs_stackbuf
      INTEGER, DIMENSION(:), INTENT(IN), POINTER         :: row_blk_sizes2enum, enum2row_blk_sizes, &
                                                            col_blk_sizes2enum, enum2col_blk_sizes
      INTEGER, DIMENSION(:), INTENT(IN)                  :: indices

      CHARACTER(len=*), PARAMETER :: routineN = 'acc_transpose_blocks_images', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk_p, col, handle, handle1, i, imat, m, &
                                                            mi, mi_max, n, nblks, nblks_total, ni, &
                                                            ni_max, offset, row, x
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: counters, filled, offsets, tmp_stack
      INTEGER, DIMENSION(:), POINTER                     :: blk_index, col_blk_sizes, &
                                                            local2global_cols, local2global_rows, &
                                                            row_blk_sizes, trs_stack

      CALL timeset(routineN, handle)
      !
      NULLIFY (local2global_rows, local2global_cols, trs_stack)
      !
      IF (trs_stackbuf%d%data_type /= dbcsr_type_int_4) &
         CPABORT("build_trs_stack: stac_buf has wrong datatype")
      !
      nblks_total = 0
      DO imat = 1, SIZE(indices)
         nblks_total = nblks_total+matrices%mats(indices(imat))%nblks
      ENDDO
      !
      ! make sure buffer from previous cannon-tick was uploaded
      CALL timeset(routineN//"_sync", handle1)
      CALL acc_event_synchronize(trs_stackbuf%d%acc_ready)
      CALL timestop(handle1)
      !
      CALL timeset(routineN//"_ensure", handle1)
      CALL dbcsr_data_ensure_size(trs_stackbuf, data_size=nblks_total, nocopy=.TRUE.)
      CALL dbcsr_data_set_size_referenced(trs_stackbuf, nblks_total)
      trs_stack => trs_stackbuf%d%i4
      CALL timestop(handle1)
      !
      mi_max = SIZE(enum2row_blk_sizes); ni_max = SIZE(enum2col_blk_sizes)
      ALLOCATE (counters(mi_max, ni_max), offsets(mi_max, ni_max))
      counters(:, :) = 0; offsets(:, :) = 0
      !
      CALL timeset(routineN//"_comp", handle1)
      IF (mi_max .NE. 1 .OR. ni_max .NE. 1) THEN
         row_blk_sizes => array_data(matrices%mats(indices(1))%row_blk_size)
         col_blk_sizes => array_data(matrices%mats(indices(1))%col_blk_size)
         ALLOCATE (tmp_stack(3, nblks_total))
         nblks_total = 0
      ENDIF
      !
      DO imat = 1, SIZE(indices)
         blk_index => matrices%mats(indices(imat))%coo_l
         nblks = matrices%mats(indices(imat))%nblks
         !
         ! collect block addresses and dimensions in a temporary stack
         ! while doing so, also count number of blocks per block-dimensions
         !
         ! Simplified algorithm for single size blocks
         IF (mi_max .EQ. 1 .AND. ni_max .EQ. 1) THEN
            DO i = 1, nblks
               blk_p = blk_index(3*(i-1)+3)
               IF (blk_p == 0) CYCLE
               counters(1, 1) = counters(1, 1)+1
               trs_stack(counters(1, 1)) = blk_p-1
            ENDDO
         ELSE
            local2global_rows => array_data(matrices%mats(indices(imat))%local_rows)
            local2global_cols => array_data(matrices%mats(indices(imat))%local_cols)
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
         ENDIF
      ENDDO
      !
      IF (mi_max .NE. 1 .OR. ni_max .NE. 1) THEN
         ! calculate offsets for first element of each sub-stack
         offset = 0
         DO mi = 1, mi_max
            DO ni = 1, ni_max
               offsets(mi, ni) = offset
               offset = offset+counters(mi, ni)
            ENDDO
         ENDDO
         !
         ! write all sub-stacks into the host-pinned buffer
         ALLOCATE (filled(mi_max, ni_max))
         filled(:, :) = 0
         DO i = 1, nblks_total
            mi = tmp_stack(1, i)
            ni = tmp_stack(2, i)
            blk_p = tmp_stack(3, i)
            x = offsets(mi, ni)+filled(mi, ni)+1
            trs_stack(x) = blk_p
            filled(mi, ni) = filled(mi, ni)+1
         ENDDO
         !sanity check
         DO ni = 1, ni_max
            DO mi = 1, mi_max
               IF (filled(mi, ni) /= counters(mi, ni)) &
                  CPABORT("acc_transpose_blocks: bug")
            END DO
         END DO
         !
         DEALLOCATE (tmp_stack)
      ENDIF
      CALL timestop(handle1)

      CALL timeset(routineN//"_sync", handle1)
      !transfer all stacks
      CALL dbcsr_data_host2dev(trs_stackbuf)
      ! make sure block-buffer is uploaded befor running the kernels
      CALL acc_stream_wait_event(trs_stackbuf%d%memory_type%acc_stream, data%d%acc_ready)
      CALL timestop(handle1)
      !
      CALL timeset(routineN//"_kernels", handle1)
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
      CALL timestop(handle1)
      !
      CALL timeset(routineN//"_sync", handle1)
      ! make sure block-buffer are not used until transpose kernels finished
      CALL acc_event_record(trs_stackbuf%d%acc_ready, trs_stackbuf%d%memory_type%acc_stream)
      CALL acc_stream_wait_event(data%d%memory_type%acc_stream, trs_stackbuf%d%acc_ready)
      CALL acc_event_record(data%d%acc_ready, data%d%memory_type%acc_stream)
      CALL timestop(handle1)

      CALL timestop(handle)
   END SUBROUTINE acc_transpose_blocks_images

! **************************************************************************************************
!> \brief Init buffer
!>
!> \param buffer ...
!> \param data_type ...
!> \param data_size ...
!> \param meta_size ...
!> \param num_data ...
!> \param requests_size ...
!> \param data_memory_type ...
!> \param trs_memory_type ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE buffer_init(buffer, data_type, &
                          data_size, meta_size, &
                          num_data, requests_size, &
                          data_memory_type, trs_memory_type)
      TYPE(dbcsr_buffer), INTENT(INOUT)                  :: buffer
      INTEGER, INTENT(IN)                                :: data_type, data_size, meta_size
      INTEGER, INTENT(IN), OPTIONAL                      :: num_data, requests_size
      TYPE(dbcsr_memtype_type), INTENT(IN)               :: data_memory_type
      TYPE(dbcsr_memtype_type), INTENT(IN), OPTIONAL     :: trs_memory_type

      INTEGER                                            :: idata, my_num_data
      LOGICAL                                            :: new_trs_stackbuf

      my_num_data = 1
      IF (PRESENT(num_data)) THEN
         my_num_data = num_data
      ELSE
         IF (dbcsr_data_valid(buffer%data_resize) .OR. ASSOCIATED(buffer%meta_resize)) &
            CPABORT("Previous data area already initialized.")
         CALL dbcsr_data_init(buffer%data_resize)
         CALL dbcsr_data_new(buffer%data_resize, data_type, memory_type=data_memory_type)
      ENDIF
      new_trs_stackbuf = PRESENT(trs_memory_type) .AND. has_acc
      !
      IF (buffer%is_valid) THEN
         ! Invalid buffers if data_type is different
         IF (dbcsr_data_get_type(buffer%data(1)) .NE. data_type .OR. &
             SIZE(buffer%data) .LT. my_num_data .OR. &
             (new_trs_stackbuf .AND. SIZE(buffer%trs_stackbuf) .LT. my_num_data)) THEN
            DO idata = 1, SIZE(buffer%data)
               CALL dbcsr_data_release(buffer%data(idata))
            ENDDO
            DEALLOCATE (buffer%data)
            CALL dbcsr_data_release(buffer%data_red3D)
            IF (new_trs_stackbuf) THEN
               DO idata = 1, SIZE(buffer%trs_stackbuf)
                  CALL dbcsr_data_release(buffer%trs_stackbuf(idata))
               ENDDO
               DEALLOCATE (buffer%trs_stackbuf)
            ENDIF
            buffer%is_valid = .FALSE.
         ENDIF
      ENDIF
      !
      IF (.NOT. buffer%is_valid) THEN
         ! First initialization
         ALLOCATE (buffer%data(my_num_data))
         IF (new_trs_stackbuf) ALLOCATE (buffer%trs_stackbuf(my_num_data))
         DO idata = 1, my_num_data
            CALL dbcsr_data_init(buffer%data(idata))
            CALL dbcsr_data_new(buffer%data(idata), data_type=data_type, &
                                data_size=data_size, memory_type=data_memory_type)
            CALL dbcsr_data_set_size_referenced(buffer%data(idata), data_size)
            IF (new_trs_stackbuf) THEN
               CALL dbcsr_data_init(buffer%trs_stackbuf(idata))
               CALL dbcsr_data_new(buffer%trs_stackbuf(idata), &
                                   data_type=dbcsr_type_int_4, data_size=2*meta_size, &
                                   memory_type=trs_memory_type)
            ENDIF
         ENDDO
         CALL dbcsr_data_init(buffer%data_red3D)
         CALL dbcsr_data_new(buffer%data_red3D, data_type=data_type, &
                             data_size=1, memory_type=memtype_mpi_buffer)
         buffer%is_valid = .TRUE.
      ELSE
         IF (PRESENT(num_data)) THEN
            DO idata = 1, my_num_data
               CALL dbcsr_data_ensure_size(buffer%data(idata), data_size, nocopy=.TRUE.)
            ENDDO
         ELSE
            CALL dbcsr_data_ensure_size(buffer%data(1), data_size, nocopy=.TRUE., &
                                        area_resize=buffer%data_resize)
         ENDIF
      ENDIF
      !
      IF (PRESENT(num_data)) THEN
         CALL ensure_array_size(buffer%meta, ub=meta_size*my_num_data, nocopy=.TRUE., &
                                memory_type=memtype_mpi_buffer)
      ELSE
         CALL ensure_array_size(buffer%meta, array_resize=buffer%meta_resize, &
                                ub=meta_size, nocopy=.TRUE., &
                                memory_type=memtype_mpi_buffer)
      ENDIF
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
      ENDIF
      !
      buffer%nrequests = 0
   END SUBROUTINE buffer_init

! **************************************************************************************************
!> \brief Release all buffers
!>
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE buffers_release()
      CALL buffer_release(buffers_1%right)
      CALL buffer_release(buffers_1%left)
      CALL buffer_release(buffers_2%right)
      CALL buffer_release(buffers_2%left)
      CALL buffer_release(buffers_orig%right)
      CALL buffer_release(buffers_orig%left)
   END SUBROUTINE buffers_release

! **************************************************************************************************
!> \brief Release buffer
!>
!> \param buffer ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE buffer_release(buffer)
      TYPE(dbcsr_buffer), INTENT(INOUT)                  :: buffer

      INTEGER                                            :: idata

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
         DO idata = 1, SIZE(buffer%data)
            CALL dbcsr_data_release(buffer%data(idata))
            IF (ALLOCATED(buffer%trs_stackbuf)) THEN
               CALL dbcsr_data_release(buffer%trs_stackbuf(idata))
            ENDIF
         ENDDO
         DEALLOCATE (buffer%data)
         IF (dbcsr_data_valid(buffer%data_resize)) THEN
            CALL dbcsr_data_release(buffer%data_resize)
         ENDIF
         IF (ALLOCATED(buffer%trs_stackbuf)) THEN
            DEALLOCATE (buffer%trs_stackbuf)
         ENDIF
         IF (dbcsr_data_valid(buffer%data_red3D)) &
            CALL dbcsr_data_release(buffer%data_red3D)
         buffer%is_valid = .FALSE.
      ENDIF
      IF (ASSOCIATED(buffer%meta)) THEN
         CALL memory_deallocate(buffer%meta, memtype_mpi_buffer)
         NULLIFY (buffer%meta)
      ENDIF
      IF (ASSOCIATED(buffer%meta_resize)) THEN
         CALL memory_deallocate(buffer%meta_resize, memtype_mpi_buffer)
         NULLIFY (buffer%meta_resize)
      ENDIF
      IF (ASSOCIATED(buffer%meta_red3D)) THEN
         CALL memory_deallocate(buffer%meta_red3D, memtype_mpi_buffer)
         NULLIFY (buffer%meta_red3D)
      ENDIF
      IF (ALLOCATED(buffer%offset)) DEALLOCATE (buffer%offset)
      IF (ALLOCATED(buffer%get_requests)) DEALLOCATE (buffer%get_requests)
   END SUBROUTINE buffer_release

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
!> \param ntotal_images ...
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
                           ntotal_images, &
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
      INTEGER, INTENT(IN)                                :: nrows_images, ncols_images, &
                                                            ntotal_images, vprow, vpcol
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: offset_images
      INTEGER, DIMENSION(:), INTENT(IN), TARGET          :: meta_buffer
      INTEGER, INTENT(IN), OPTIONAL                      :: uf, ul
      TYPE(dbcsr_imagedistribution_obj), INTENT(INOUT)   :: imgdist
      LOGICAL, INTENT(IN)                                :: do_merge_rows
      INTEGER, DIMENSION(:), INTENT(IN)                  :: global_indices
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: vunmerged
      INTEGER, DIMENSION(:, :), INTENT(IN), OPTIONAL     :: vmerged
      LOGICAL, DIMENSION(:), INTENT(INOUT), OPTIONAL     :: common_unmerged
      INTEGER, DIMENSION(:), INTENT(INOUT), OPTIONAL     :: indices_common_unmerged
      INTEGER, INTENT(INOUT), OPTIONAL                   :: ncommons_unmerged
      INTEGER, INTENT(IN), OPTIONAL                      :: nthreads

      CHARACTER(len=*), PARAMETER :: routineN = 'merge_images', routineP = moduleN//':'//routineN

      INTEGER :: handle, ioffset, ithread, kf, ki, kl, mi, my_nthreads, my_uf, my_ul, my_vpcol, &
         my_vprow, nimages_merged, nimages_unmerged, offset_data, offset_meta, size_index_merged, &
         size_index_unmerged, ui, v_merged
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ioffset_first
      INTEGER, DIMENSION(:), POINTER                     :: unmerged_indices
      LOGICAL                                            :: noempty
      TYPE(dbcsr_type), POINTER                          :: image, image_obj

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
      my_uf = 1
      IF (PRESENT(uf)) my_uf = uf
      my_ul = 1
      IF (PRESENT(ul)) my_ul = ul
      !
      ALLOCATE (ioffset_first(my_ul))
      ioffset_first(1) = 0
      DO ui = 1, my_ul-1
         IF (PRESENT(vunmerged)) THEN
            ioffset_first(ui+1) = ioffset_first(ui)+vunmerged(ui)
         ELSE
            ioffset_first(ui+1) = ioffset_first(ui)+1
         ENDIF
      ENDDO
      !
      IF (PRESENT(ncommons_unmerged)) ncommons_unmerged = 0
      !
!$OMP PARALLEL DO DEFAULT (NONE) &
!$OMP PRIVATE(image_obj,image,ioffset,ui,mi,v_merged, &
!$OMP         offset_meta,offset_data,unmerged_indices, &
!$OMP         ithread,kf,kl,ki,noempty) &
!$OMP SHARED(nimages_unmerged,images,global_indices, &
!$OMP        size_index_merged,vunmerged,do_merge_rows, &
!$OMP        my_vprow,my_vpcol,nimages_merged,nthreads, &
!$OMP        my_nthreads,ioffset_first,vmerged,offset_images, &
!$OMP        meta_buffer,size_index_unmerged,imgdist,my_uf,my_ul, &
!$OMP        ncommons_unmerged,common_unmerged,indices_common_unmerged, &
!$OMP        ntotal_images) &
!$OMP SCHEDULE(guided)
      images_unmerge: DO ui = my_uf, my_ul
         image_obj => images%mats(ui)
         image => image_obj
         CALL set_empty_meta_index(image%index, global_indices, size_index_merged)
         noempty = .TRUE.
         IF (PRESENT(vunmerged)) THEN
            noempty = vunmerged(ui) .NE. 0
         ENDIF
         IF (noempty) THEN
            !
            IF (PRESENT(common_unmerged) .AND. PRESENT(ncommons_unmerged) .AND. PRESENT(indices_common_unmerged)) THEN
!$OMP CRITICAL(merge_images_common)
               IF (common_unmerged(ui)) THEN
                  ncommons_unmerged = ncommons_unmerged+1
                  indices_common_unmerged(ncommons_unmerged) = ui
               ELSE
                  common_unmerged(ui) = .TRUE.
               ENDIF
!$OMP END CRITICAL(merge_images_common)
            ENDIF
            !
            ! Virtual coords
            IF (do_merge_rows) THEN
               image%index(dbcsr_slot_home_vprow) = my_vprow
               image%index(dbcsr_slot_home_vpcol) = MOD(my_vpcol+ui, ntotal_images)
            ELSE
               image%index(dbcsr_slot_home_vprow) = MOD(my_vprow+ui, ntotal_images)
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
                     offset_meta = offset_images(imeta, ioffset)
                     offset_data = offset_images(idata, ioffset)
                     ! Check for filtered images
                     IF (offset_data .LT. offset_images(idata, ioffset+1)) THEN
                        unmerged_indices => meta_buffer(offset_meta+dbcsr_slot_size: &
                                                        offset_meta+meta_buffer(offset_meta+dbcsr_slot_size))
                        IF (ithread .EQ. 1) THEN
                           ! nblks, nze
                           image%index(dbcsr_slot_nblks) = &
                              image%index(dbcsr_slot_nblks)+ &
                              unmerged_indices(dbcsr_slot_nblks)
                           image%index(dbcsr_slot_nze) = &
                              image%index(dbcsr_slot_nze)+ &
                              offset_images(idata, ioffset+1)- &
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
               ioffset = ioffset_first(ui)+1
               offset_meta = offset_images(imeta, ioffset)
               ! Check for filtered images
               offset_data = offset_images(idata, ioffset)
               IF (offset_data .LT. offset_images(idata, ioffset+1)) THEN
                  image%index(dbcsr_slot_size) = &
                     meta_buffer(offset_meta+dbcsr_slot_size)- &
                     dbcsr_slot_nblks+dbcsr_num_slots
                  image%index(dbcsr_slot_nblks) = &
                     meta_buffer(offset_meta+dbcsr_slot_nblks)
                  image%index(dbcsr_slot_nze) = &
                     offset_images(idata, ioffset+1)-offset_data
                  image%index(dbcsr_num_slots) = image%index(dbcsr_slot_size)
                  !
                  image%index(dbcsr_num_slots+1:image%index(dbcsr_num_slots)) = &
                     meta_buffer(offset_meta+dbcsr_slot_nblks+1: &
                                 offset_meta+meta_buffer(offset_meta+dbcsr_slot_size))
                  IF (offset_data .GT. 0) THEN
                     image%index(size_index_merged+3:image%index(dbcsr_num_slots):3) = &
                        image%index(size_index_merged+3:image%index(dbcsr_num_slots):3)+ &
                        offset_data
                  ENDIF
               ENDIF
            ENDIF merge
            !
            ! Reset
!$OMP CRITICAL(merge_images_reset_vlocals)
            IF (nimages_merged .GT. 1) THEN
               CALL dbcsr_reset_vlocals(image_obj, imgdist, do_rows=.NOT. do_merge_rows)
            ELSE
               CALL dbcsr_reset_vlocals(image_obj, imgdist)
            ENDIF
!$OMP END CRITICAL(merge_images_reset_vlocals)
         ENDIF
         !
         ! Repoint index
         image%nblks = 0
         image%nze = 0
         CALL dbcsr_repoint_index(image)
         image%valid = .TRUE.
      ENDDO images_unmerge
!$OMP END PARALLEL DO
      !
      DEALLOCATE (ioffset_first)
      CALL timestop(handle)
   END SUBROUTINE merge_images

! **************************************************************************************************
!> \brief ...
!> \param nimages ...
!> \param nlayers ...
!> \param nproc ...
!> \param refs ...
!> \param refs_layer3D ...
! **************************************************************************************************
   SUBROUTINE remap_size_layers3D(nimages, nlayers, nproc, &
                                  refs, refs_layer3D)
      INTEGER, INTENT(IN)                                :: nimages, nlayers, nproc
      INTEGER, DIMENSION(nimages*nlayers, 0:nproc-1), &
         INTENT(IN), TARGET                              :: refs
      INTEGER, DIMENSION(:, :), INTENT(OUT), POINTER     :: refs_layer3D

      INTEGER                                            :: ilayer, image, iproc

      IF (nimages .EQ. 1) THEN
         refs_layer3D => refs
         RETURN
      ENDIF
      ! Remap
      ALLOCATE (refs_layer3D(nlayers, 0:nimages*nproc-1))
!$OMP PARALLEL DO DEFAULT (NONE) &
!$OMP SHARED (nproc, nimages, nlayers, &
!$OMP         refs_layer3D, refs) &
!$OMP PRIVATE (iproc,image,ilayer)
      DO iproc = 0, nproc-1
         DO image = 0, nimages-1
            DO ilayer = 1, nlayers
               refs_layer3D(ilayer, image+iproc*nimages) = refs((ilayer-1)*nimages+image+1, iproc)
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
   END SUBROUTINE remap_size_layers3D

! **************************************************************************************************
!> \brief ...
!> \param nimages ...
!> \param nlayers ...
!> \param nproc ...
!> \param refs ...
!> \param refs_layer3D ...
! **************************************************************************************************
   SUBROUTINE remap_displ_layers3D(nimages, nlayers, nproc, &
                                   refs, refs_layer3D)
      INTEGER, INTENT(IN)                                :: nimages, nlayers, nproc
      INTEGER, DIMENSION(idata:imeta, nimages*nlayers, 0&
         :nproc-1), INTENT(IN), TARGET                   :: refs
      INTEGER, DIMENSION(:, :, :), INTENT(OUT), POINTER  :: refs_layer3D

      INTEGER                                            :: ilayer, image, iproc

      IF (nimages .EQ. 1) THEN
         refs_layer3D => refs
         RETURN
      ENDIF
      ! Remap
      ALLOCATE (refs_layer3D(idata:imeta, nlayers, 0:nimages*nproc-1))
!$OMP PARALLEL DO DEFAULT (NONE) &
!$OMP SHARED (nproc, nimages, nlayers, &
!$OMP         refs_layer3D, refs) &
!$OMP PRIVATE (iproc,image,ilayer)
      DO iproc = 0, nproc-1
         DO image = 0, nimages-1
            DO ilayer = 1, nlayers
               refs_layer3D(:, ilayer, image+iproc*nimages) = refs(:, (ilayer-1)*nimages+image+1, iproc)
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
   END SUBROUTINE remap_displ_layers3D

#include "dbcsr_mm_3d_d.f90"
#include "dbcsr_mm_3d_z.f90"
#include "dbcsr_mm_3d_s.f90"
#include "dbcsr_mm_3d_c.f90"

END MODULE dbcsr_mm_3d

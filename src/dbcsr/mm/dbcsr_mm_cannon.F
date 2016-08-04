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
! **************************************************************************************************
MODULE dbcsr_mm_cannon
   USE acc_event,                       ONLY: acc_event_synchronize
   USE array_types,                     ONLY: array_data,&
                                              array_exists,&
                                              array_i1d_obj
   USE dbcsr_block_operations,          ONLY: dbcsr_block_conjg,&
                                              dbcsr_block_copy_aa,&
                                              dbcsr_block_real_neg,&
                                              dbcsr_block_scale,&
                                              dbcsr_block_transpose_aa
   USE dbcsr_config,                    ONLY: has_acc,&
                                              use_comm_thread,&
                                              use_mpi_exp
   USE dbcsr_data_methods,              ONLY: &
        dbcsr_data_clear_pointer, dbcsr_data_ensure_size, dbcsr_data_get_size, &
        dbcsr_data_get_size_referenced, dbcsr_data_hold, dbcsr_data_host2dev, dbcsr_data_init, &
        dbcsr_data_new, dbcsr_data_release, dbcsr_data_set_pointer, &
        dbcsr_data_set_size_referenced, dbcsr_get_data_p_c, dbcsr_get_data_p_d, &
        dbcsr_get_data_p_s, dbcsr_get_data_p_z, dbcsr_scalar_are_equal, dbcsr_scalar_negative, &
        dbcsr_scalar_one
   USE dbcsr_data_types,                ONLY: dbcsr_datatype_sizeof
   USE dbcsr_dist_methods,              ONLY: dbcsr_distribution_col_dist,&
                                              dbcsr_distribution_has_threads,&
                                              dbcsr_distribution_mp,&
                                              dbcsr_distribution_row_dist
   USE dbcsr_dist_operations,           ONLY: dbcsr_get_local_vcols,&
                                              dbcsr_get_local_vrows,&
                                              dbcsr_reset_vlocals,&
                                              image_calculator
   USE dbcsr_index_operations,          ONLY: dbcsr_count_row_index,&
                                              dbcsr_has_local_row_index,&
                                              dbcsr_index_compact,&
                                              dbcsr_index_prune_deleted,&
                                              dbcsr_make_index_list,&
                                              dbcsr_make_index_local_row,&
                                              dbcsr_repoint_index
   USE dbcsr_io,                        ONLY: dbcsr_print
   USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop
   USE dbcsr_mem_methods,               ONLY: dbcsr_mempool_ensure_capacity
   USE dbcsr_methods,                   ONLY: &
        dbcsr_destroy_array, dbcsr_distribution, dbcsr_get_data_type, dbcsr_get_index_memory_type, &
        dbcsr_has_symmetry, dbcsr_image_dist_hold, dbcsr_image_dist_init, &
        dbcsr_image_dist_release, dbcsr_init, dbcsr_nblkcols_local, dbcsr_nblkcols_total, &
        dbcsr_nblkrows_local, dbcsr_nblkrows_total, dbcsr_release, dbcsr_valid_index
   USE dbcsr_mm_common,                 ONLY: &
        acc_transpose_blocks, count_mpi_statistics, dbcsr_mm_multrec_type_p, dbcsr_mpi_statistics, &
        max_memory, memtype_abpanel_1, memtype_abpanel_2, memtype_mpi_buffer, memtype_trsbuffer_1, &
        memtype_trsbuffer_2, rec_sort_index, setup_buffer_matrix
   USE dbcsr_mm_multrec,                ONLY: dbcsr_mm_multrec_finalize,&
                                              dbcsr_mm_multrec_init,&
                                              dbcsr_mm_multrec_multiply,&
                                              dbcsr_mm_multrec_phaseout
   USE dbcsr_mp_methods,                ONLY: &
        dbcsr_mp_grid_setup, dbcsr_mp_group, dbcsr_mp_has_subgroups, dbcsr_mp_my_col_group, &
        dbcsr_mp_my_row_group, dbcsr_mp_mynode, dbcsr_mp_mypcol, dbcsr_mp_myprow, dbcsr_mp_npcols, &
        dbcsr_mp_nprows, dbcsr_mp_numnodes, dbcsr_mp_pgrid
   USE dbcsr_mp_operations,             ONLY: dbcsr_irecv_any,&
                                              dbcsr_isend_any,&
                                              hybrid_alltoall_any,&
                                              hybrid_alltoall_i1
   USE dbcsr_operations,                ONLY: dbcsr_copy,&
                                              dbcsr_crop_matrix
   USE dbcsr_ptr_util,                  ONLY: ensure_array_size
   USE dbcsr_transformations,           ONLY: dbcsr_make_dense_low
   USE dbcsr_types,                     ONLY: &
        dbcsr_2d_array_type, dbcsr_data_obj, dbcsr_distribution_obj, dbcsr_imagedistribution_obj, &
        dbcsr_iterator, dbcsr_memtype_type, dbcsr_meta_size, dbcsr_mp_obj, dbcsr_obj, &
        dbcsr_scalar_type, dbcsr_slot_home_coli, dbcsr_slot_home_pcol, dbcsr_slot_home_prow, &
        dbcsr_slot_home_rowi, dbcsr_slot_home_vpcol, dbcsr_slot_home_vprow, dbcsr_slot_size, &
        dbcsr_type, dbcsr_type_complex_4, dbcsr_type_complex_8, dbcsr_type_int_4, &
        dbcsr_type_no_symmetry, dbcsr_type_real_4, dbcsr_type_real_8
   USE dbcsr_util,                      ONLY: count_bins,&
                                              dbcsr_checksum
   USE dbcsr_work_operations,           ONLY: dbcsr_create,&
                                              dbcsr_special_finalize,&
                                              dbcsr_work_create
   USE kinds,                           ONLY: dp,&
                                              int_8,&
                                              real_4,&
                                              real_8,&
                                              sp
   USE machine,                         ONLY: default_output_unit,&
                                              m_memory
   USE message_passing,                 ONLY: mp_allgather,&
                                              mp_alltoall,&
                                              mp_irecv,&
                                              mp_isend,&
                                              mp_request_null,&
                                              mp_sum,&
                                              mp_testany,&
                                              mp_waitall
#include "../../base/base_uses.f90"

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_cannon'
   LOGICAL, PARAMETER :: debug_mod = .FALSE.
   LOGICAL, PARAMETER :: careful_mod = .FALSE.
   REAL, PARAMETER :: huge_norm = HUGE(1.0)**(1.0/3.0)

   INTERFACE dbcsr_switch
      MODULE PROCEDURE dbcsr_switch_sets
      MODULE PROCEDURE dbcsr_switch_d_ptrs
   END INTERFACE

   PUBLIC :: multiply_cannon, make_m2s

CONTAINS

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
      IF (use_mpi_exp) &
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
   SUBROUTINE multiply_cannon(left_set, right_set, product_matrix, &
                              retain_sparsity, &
                              filter_eps, flop)
      TYPE(dbcsr_2d_array_type), POINTER                 :: left_set, right_set
      TYPE(dbcsr_obj), INTENT(INOUT)                     :: product_matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: retain_sparsity
      REAL(kind=real_8), INTENT(in), OPTIONAL            :: filter_eps
      INTEGER(KIND=int_8), INTENT(OUT)                   :: flop

      CHARACTER(len=*), PARAMETER :: routineN = 'multiply_cannon', &
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
      CALL setup_rec_index_2d(left_set, left_row_nimages, left_col_nimages)
      CALL setup_rec_index_2d(right_set, right_row_nimages, right_col_nimages)
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
   END SUBROUTINE multiply_cannon

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
      INTEGER, INTENT(IN)                                :: index_size
      INTEGER, INTENT(IN), OPTIONAL                      :: data_size

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

#include "dbcsr_mm_cannon_d.f90"
#include "dbcsr_mm_cannon_z.f90"
#include "dbcsr_mm_cannon_s.f90"
#include "dbcsr_mm_cannon_c.f90"

END MODULE dbcsr_mm_cannon

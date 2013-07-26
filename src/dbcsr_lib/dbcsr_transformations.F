!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR transformations
!> \author  Urban Borstnik
!> \date    2010-02-18
!> \version 0.9
!>
!> <b>Modification history:</b>
!> - 2010-02-18 Moved from dbcsr_util and dbcsr_operations
! *****************************************************************************
MODULE dbcsr_transformations

  USE array_types,                     ONLY: array_data,&
                                             array_hold,&
                                             array_i1d_obj,&
                                             array_release
  USE dbcsr_block_access,              ONLY: dbcsr_access_flush,&
                                             dbcsr_access_start,&
                                             dbcsr_access_stop,&
                                             dbcsr_get_block_p,&
                                             dbcsr_put_block
  USE dbcsr_block_operations,          ONLY: dbcsr_block_conjg,&
                                             dbcsr_block_partial_copy,&
                                             dbcsr_block_real_neg,&
                                             dbcsr_block_scale,&
                                             dbcsr_block_transpose,&
                                             dbcsr_data_clear,&
                                             dbcsr_data_copy,&
                                             dbcsr_data_set
  USE dbcsr_data_methods,              ONLY: &
       dbcsr_data_clear_pointer, dbcsr_data_ensure_size, &
       dbcsr_data_get_memory_type, dbcsr_data_get_size, dbcsr_data_get_type, &
       dbcsr_data_hold, dbcsr_data_init, dbcsr_data_new, dbcsr_data_release, &
       dbcsr_data_set_pointer, dbcsr_data_set_size_referenced, &
       dbcsr_get_data, dbcsr_scalar_negative, dbcsr_scalar_one, &
       dbcsr_type_1d_to_2d
  USE dbcsr_data_operations,           ONLY: dbcsr_copy_sort_data
  USE dbcsr_dist_operations,           ONLY: &
       dbcsr_get_local_cols, dbcsr_get_local_rows, dbcsr_get_local_vcols, &
       dbcsr_get_local_vrows, dbcsr_get_stored_coordinates, &
       dbcsr_reblocking_targets, dbcsr_transpose_dims, &
       dbcsr_transpose_distribution, image_calculator
  USE dbcsr_error_handling
  USE dbcsr_index_operations,          ONLY: dbcsr_addto_index_array,&
                                             dbcsr_clearfrom_index_array,&
                                             dbcsr_index_prune_deleted,&
                                             dbcsr_repoint_index,&
                                             make_dense_index,&
                                             make_undense_index,&
                                             transpose_index_local
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_kinds,                     ONLY: dp,&
                                             int_8,&
                                             sp
  USE dbcsr_message_passing,           ONLY: mp_allgather,&
                                             mp_alltoall
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_sizes, dbcsr_distribution, &
       dbcsr_distribution_col_dist, dbcsr_distribution_has_threads, &
       dbcsr_distribution_hold, dbcsr_distribution_local_cols, &
       dbcsr_distribution_local_rows, dbcsr_distribution_make_threads, &
       dbcsr_distribution_mp, dbcsr_distribution_ncols, &
       dbcsr_distribution_nlocal_cols, dbcsr_distribution_nlocal_rows, &
       dbcsr_distribution_no_threads, dbcsr_distribution_nrows, &
       dbcsr_distribution_release, dbcsr_distribution_row_dist, &
       dbcsr_get_data_size_referenced, dbcsr_get_data_type, &
       dbcsr_get_matrix_type, dbcsr_has_symmetry, dbcsr_image_dist_hold, &
       dbcsr_image_dist_release, dbcsr_init, dbcsr_is_initialized, &
       dbcsr_mp_grid_remove, dbcsr_mp_grid_setup, dbcsr_mp_group, &
       dbcsr_mp_has_subgroups, dbcsr_mp_my_col_group, dbcsr_mp_my_row_group, &
       dbcsr_mp_mynode, dbcsr_mp_mypcol, dbcsr_mp_myprow, dbcsr_mp_npcols, &
       dbcsr_mp_nprows, dbcsr_mp_numnodes, dbcsr_mp_pgrid, &
       dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_nfullcols_total, &
       dbcsr_nfullrows_total, dbcsr_release, dbcsr_row_block_sizes, &
       dbcsr_switch_data_area, dbcsr_valid_index
  USE dbcsr_mp_operations,             ONLY: dbcsr_allgatherv,&
                                             hybrid_alltoall_any,&
                                             hybrid_alltoall_c1,&
                                             hybrid_alltoall_d1,&
                                             hybrid_alltoall_i1,&
                                             hybrid_alltoall_s1,&
                                             hybrid_alltoall_z1
  USE dbcsr_ptr_util,                  ONLY: pointer_view
  USE dbcsr_types,                     ONLY: &
       dbcsr_2d_array_type, dbcsr_data_obj, dbcsr_distribution_obj, &
       dbcsr_imagedistribution_obj, dbcsr_iterator, dbcsr_meta_size, &
       dbcsr_mp_obj, dbcsr_obj, dbcsr_repl_col, dbcsr_repl_full, &
       dbcsr_repl_none, dbcsr_repl_row, dbcsr_scalar_type, dbcsr_slot_blk_p, &
       dbcsr_slot_col_i, dbcsr_slot_home_coli, dbcsr_slot_home_pcol, &
       dbcsr_slot_home_prow, dbcsr_slot_home_rowi, dbcsr_slot_home_vpcol, &
       dbcsr_slot_home_vprow, dbcsr_slot_nblks, dbcsr_slot_nze, &
       dbcsr_slot_row_p, dbcsr_type, dbcsr_type_complex_4, &
       dbcsr_type_complex_8, dbcsr_type_no_symmetry, dbcsr_type_real_4, &
       dbcsr_type_real_8
  USE dbcsr_util,                      ONLY: convert_sizes_to_offsets,&
                                             dbcsr_checksum,&
                                             dbcsr_pack_meta,&
                                             dbcsr_unpack_meta,&
                                             dbcsr_verify_matrix,&
                                             get_internal_offsets,&
                                             global_offsets_to_local,&
                                             nfull_elements
  USE dbcsr_work_operations,           ONLY: dbcsr_create,&
                                             dbcsr_finalize,&
                                             dbcsr_special_finalize,&
                                             dbcsr_work_create

  !$ USE OMP_LIB

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_transformations'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034

  LOGICAL, PARAMETER, PRIVATE          :: careful_mod = .FALSE.

  PUBLIC :: dbcsr_desymmetrize_deep,&
            dbcsr_new_transposed,&
            dbcsr_complete_redistribute,&
            dbcsr_redistribute, dbcsr_make_untransposed_blocks
  PUBLIC :: dbcsr_replicate_all, dbcsr_distribute, dbcsr_datablock_redistribute, dbcsr_replicate
  PUBLIC :: dbcsr_make_dense, dbcsr_make_undense, dbcsr_make_images_dense
  PUBLIC :: dbcsr_make_images

CONTAINS


! *****************************************************************************
!> \brief Transposes a DBCSR matrix.
!> \par Distribution options
!>      By default the distribution is transposed. If transpose_distribution
!>      is false, then an undetermined distribution is created that is
!>      compatible with the same process grid.
!> \param[out] transposed     transposed DBCSR matrix
!> \param[in] normal          input DBCSR matrix
!> \param[in] shallow_data_copy         (optional) only shallow data_copy;
!>                                      default is no; if set, the
!>                                      transpose_data option is ignored
!> \param[in] transpose_data  (optional) transpose data blocks, default is True
!> \param[in] transpose_distribution    (optional) transpose the distribution
!>                                      from the input matrix, default is True
!> \param[in] transpose_index           (optional) transpose the index
!>                                      (default=yes) or turn it into BCSC
!> \param[in] use_distribution          (optional) use this distribution
!> \param[in] redistribute              (optional) redistributes the matrix;
!>                                      default is .TRUE. unless shallow or
!>                                      transpose_distribution are set.
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_new_transposed (transposed, normal, shallow_data_copy,&
       transpose_data, transpose_distribution, transpose_index,&
       use_distribution, redistribute, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: transposed
    TYPE(dbcsr_obj), INTENT(IN)              :: normal
    LOGICAL, INTENT(IN), OPTIONAL            :: shallow_data_copy, &
                                                transpose_data, &
                                                transpose_distribution, &
                                                transpose_index
    TYPE(dbcsr_distribution_obj), &
      INTENT(IN), OPTIONAL                   :: use_distribution
    LOGICAL, INTENT(IN), OPTIONAL            :: redistribute
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_new_transposed', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handler, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: blk_p
    LOGICAL                                  :: redist, shallow, tr_blocks, &
                                                tr_dist, tr_index
    TYPE(dbcsr_distribution_obj)             :: new_dist
    TYPE(dbcsr_obj)                          :: t2

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (dbcsr_valid_index(normal), dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Matrix does not exist.",__LINE__,error)
    CALL dbcsr_access_flush (normal, error=error)
    ! Internalize options
    shallow = .FALSE.
    IF (PRESENT (shallow_data_copy)) shallow = shallow_data_copy
    tr_blocks = .TRUE.
    IF (PRESENT (transpose_data)) tr_blocks = transpose_data
    tr_dist = .TRUE.
    IF (PRESENT (transpose_distribution)) tr_dist = transpose_distribution
    tr_index = .TRUE.
    IF (PRESENT (transpose_index)) tr_index = transpose_index
    ! Prepare the distribution for the transposed matrix
    IF (PRESENT (use_distribution)) THEN
       CALL dbcsr_assert (dbcsr_distribution_nrows (use_distribution),&
            "EQ", dbcsr_distribution_ncols (normal%m%dist), dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Given distribution must be compatible"&
            //" with the current distribution",__LINE__,error)
       CALL dbcsr_assert (dbcsr_distribution_ncols (use_distribution),&
            "EQ", dbcsr_distribution_nrows (normal%m%dist), dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Given distribution must be compatible"&
            //" with the current distribution",__LINE__,error)
       new_dist = use_distribution
       CALL dbcsr_distribution_hold (new_dist)
    ELSE
       IF (tr_dist) THEN
          CALL dbcsr_transpose_distribution (new_dist, normal%m%dist)
       ELSE
          CALL dbcsr_transpose_dims (new_dist, normal%m%dist)
       ENDIF
    ENDIF
    ! Create the transposed matrix
    CALL dbcsr_create (transposed, name="transposed "//TRIM(normal%m%name),&
         template = normal, &
         dist = new_dist,&
         matrix_type = dbcsr_get_matrix_type(normal),&
         row_blk_size = normal%m%col_blk_size,&
         col_blk_size = normal%m%row_blk_size,&
         error=error)
    CALL dbcsr_distribution_release (new_dist)
    ! Reserve the space for the new indices.
    IF (tr_index) THEN
       CALL dbcsr_addto_index_array (transposed%m, dbcsr_slot_row_p,&
            reservation=transposed%m%nblkrows_total+1, extra=transposed%m%nblks*2,&
            error=error)
    ELSE
       CALL dbcsr_addto_index_array (transposed%m, dbcsr_slot_row_p,&
            reservation=normal%m%nblkrows_total+1, extra=transposed%m%nblks*2,&
            error=error)
    ENDIF
    CALL dbcsr_addto_index_array (transposed%m, dbcsr_slot_col_i,&
         reservation=normal%m%nblks, error=error)
    CALL dbcsr_addto_index_array (transposed%m, dbcsr_slot_blk_p,&
         reservation=normal%m%nblks, error=error)
    CALL dbcsr_repoint_index (transposed%m)
    IF (.NOT. shallow) THEN
       CALL dbcsr_data_ensure_size (transposed%m%data_area,&
            dbcsr_get_data_size_referenced(normal),&
            nocopy=.TRUE., error=error)
    ENDIF
    !
    transposed%m%nblks = normal%m%nblks
    transposed%m%nze = normal%m%nze
    transposed%m%index(dbcsr_slot_nblks) = normal%m%nblks
    transposed%m%index(dbcsr_slot_nze) = normal%m%nze
    ! Transpose the local index.
    ALLOCATE (blk_p (normal%m%nblks), stat=stat)
    CALL dbcsr_assert(stat == 0, dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "blk_p",__LINE__,error)
    IF (tr_index) THEN
       CALL transpose_index_local (transposed%m%row_p, transposed%m%col_i,&
            normal%m%row_p, normal%m%col_i, blk_p, normal%m%blk_p)
       IF (dbg) THEN
          WRITE(*,*)'orig. row_p',normal%m%row_p
          WRITE(*,*)'orig. col_i',normal%m%col_i
          WRITE(*,*)'orig. blk_p',normal%m%blk_p
          WRITE(*,*)'new . row_p',transposed%m%row_p
          WRITE(*,*)'new . col_i',transposed%m%col_i
          WRITE(*,*)'new . blk_p',blk_p!transposed%m%blk_p
       ENDIF
    ELSE
       transposed%m%row_p(:) = normal%m%row_p(:)
       transposed%m%col_i(:) = normal%m%col_i(:)
       blk_p(:) = normal%m%blk_p(:)
       !transposed%m%transpose = .TRUE.
    ENDIF
    ! Copy the data
    IF (shallow) THEN
       CALL dbcsr_switch_data_area (transposed, normal%m%data_area, error=error)
       transposed%m%blk_p(1:transposed%m%nblks) =&
            -blk_p(1:transposed%m%nblks)
    ELSE
       CALL dbcsr_copy_sort_data (transposed%m%blk_p, blk_p, transposed%m%row_p,&
            transposed%m%col_i, array_data (transposed%m%row_blk_size),&
            array_data (transposed%m%col_blk_size),&
            transposed%m%data_area, normal%m%data_area,&
            mark_transposed=.not.tr_blocks,&
            transpose_blocks=tr_blocks)
    ENDIF
    transposed%m%valid = .TRUE.
    !CALL dbcsr_copy_sort_data (transposed%m%blk_p, blk_p, transposed%m%row_p,&
    !     transposed%m%col_i, array_data (transposed%m%row_blk_size),&
    !     array_data (transposed%m%col_blk_size),&
    !     transposed%m%data_area, normal%m%data_area,&
    !     transpose_blocks=.TRUE.)
    !
1315 FORMAT (I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5,1X,I5)
    IF (dbg) THEN
       WRITE(*,*)'new FINAL index'
       WRITE(*,1315)transposed%m%row_p
       WRITE(*,1315)transposed%m%col_i
       WRITE(*,1315)transposed%m%blk_p
    ENDIF
    !
    IF (tr_index) DEALLOCATE (blk_p)
    !
    IF (PRESENT (redistribute)) THEN
       redist = redistribute
    ELSE
       redist = .NOT. tr_dist .AND. .NOT. shallow
    ENDIF
    IF (redist) THEN
       !write (*,*)routineN//" redistributing"
       CALL dbcsr_init (t2)
       CALL dbcsr_create (t2, template=transposed, error=error)
       CALL dbcsr_redistribute (transposed, t2, error=error)
       CALL dbcsr_release (transposed)
       transposed = t2
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_new_transposed

! *****************************************************************************
!> \brief Duplicates data in symmetric matrix to make it normal (w.r.t. data
!>        structure
!> \param[in] sm              input symmetric matrix
!> \param[out] desm           desymmetrized matrix
!> \param[in] transpose_data  (optional) make all data blocks be non-transposed
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_desymmetrize_deep(sm, desm, untransposed_data, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: sm
    TYPE(dbcsr_obj), INTENT(INOUT)           :: desm
    LOGICAL, INTENT(IN), OPTIONAL            :: untransposed_data
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_desymmetrize_deep', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, idata = 2, &
                                                imeta = 1, metalen = 3

    INTEGER :: blk, blk_l, blk_p, blk_ps, blks, col, col_size, dst_p, &
      error_handler, mp_group, nsymmetries, numproc, nze, pcol, prow, row, &
      row_size, src_p, stored_col, stored_row, symmetry_i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rd_disp, recv_meta, rm_disp, &
                                                sd_disp, sdp, send_meta, &
                                                sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_count, send_count, &
                                                total_recv_count, &
                                                total_send_count
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, col_dist, &
                                                row_blk_size, row_dist
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    LOGICAL                                  :: make_untr, tr
    TYPE(dbcsr_data_obj)                     :: recv_data, send_data
    TYPE(dbcsr_distribution_obj)             :: target_dist
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_valid_index (sm),&
         dbcsr_fatal_level, dbcsr_caller_error, routineN,&
         "Input matrix is invalid.",__LINE__,error)
    CALL dbcsr_assert (dbcsr_valid_index (desm),&
         dbcsr_fatal_level, dbcsr_caller_error, routineN,&
         "Target matrix must be created.",__LINE__,error)
    CALL dbcsr_assert ("NOT", dbcsr_has_symmetry (desm),&
         dbcsr_fatal_level, dbcsr_caller_error, routineN,&
         "Target matrix must not be symmetric.",__LINE__,error)
    CALL dbcsr_access_flush (sm, error=error)
    CALL dbcsr_access_flush (desm, error=error)
    IF (PRESENT (untransposed_data)) THEN
       make_untr = untransposed_data
    ELSE
       make_untr = .FALSE.
    ENDIF
    nsymmetries = 1
    IF (sm%m%symmetry) THEN
       nsymmetries = 2
    ENDIF
    row_blk_size => array_data (sm%m%row_blk_size)
    col_blk_size => array_data (sm%m%col_blk_size)
    target_dist = sm%m%dist
    row_dist => array_data (dbcsr_distribution_row_dist (target_dist))
    col_dist => array_data (dbcsr_distribution_col_dist (target_dist))
    mp_obj = dbcsr_distribution_mp (target_dist)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    numproc = dbcsr_mp_numnodes (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    !
    ALLOCATE (send_count(2, 0:numproc-1))
    ALLOCATE (recv_count(2, 0:numproc-1))
    ALLOCATE (total_send_count(2, 0:numproc-1))
    ALLOCATE (total_recv_count(2, 0:numproc-1))
    ALLOCATE (sdp(0:numproc-1))
    ALLOCATE (sd_disp(0:numproc-1))
    ALLOCATE (smp(0:numproc-1))
    ALLOCATE (sm_disp(0:numproc-1))
    ALLOCATE (rd_disp(0:numproc-1))
    ALLOCATE (rm_disp(0:numproc-1))
    !
    desm%m%negate_real = sm%m%negate_real
    desm%m%negate_imaginary = sm%m%negate_imaginary
    ! Count initial sizes for sending.
    send_count(:,:) = 0
    CALL dbcsr_iterator_start (iter, sm)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk,&
            row_size=row_size, col_size=col_size)
       DO symmetry_i = 1, nsymmetries
          IF (symmetry_i .EQ. 1) THEN
             stored_row = row ; stored_col = col
          ELSE
             IF (row .EQ. col) CYCLE
             stored_row = col ; stored_col = row
          ENDIF
          ! Where do we send this block?
          prow = row_dist(stored_row)
          pcol = col_dist(stored_col)
          dst_p = blacs2mpi(prow, pcol)
          nze = row_size * col_size
          send_count(imeta, dst_p) = send_count(imeta, dst_p) + 1
          send_count(idata, dst_p) = send_count(idata, dst_p) + nze
       ENDDO ! symmetry_i
    ENDDO ! iter
    CALL dbcsr_iterator_stop (iter)
    !
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    !
    ! Allocate data structures needed for data exchange.
    CALL dbcsr_data_init (recv_data)
    CALL dbcsr_data_new (recv_data, dbcsr_get_data_type (sm),&
         SUM(recv_count(idata, :)))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(imeta, :))))
    CALL dbcsr_data_init (send_data)
    CALL dbcsr_data_new (send_data, dbcsr_get_data_type (sm),&
         SUM(send_count(idata, :)))
    ALLOCATE (send_meta(metalen*SUM(send_count(imeta, :))))
    !
    ! Fill in the meta data structures and copy the data.
    DO dst_p = 0, numproc-1
       total_send_count(imeta, dst_p) = send_count (imeta, dst_p)
       total_send_count(idata, dst_p) = send_count (idata, dst_p)
       total_recv_count(imeta, dst_p) = recv_count (imeta, dst_p)
       total_recv_count(idata, dst_p) = recv_count (idata, dst_p)
    ENDDO
    sd_disp = -1 ; sm_disp = -1
    rd_disp = -1 ; rm_disp = -1
    sd_disp(0) = 1 ; sm_disp(0) = 1
    rd_disp(0) = 1 ; rm_disp(0) = 1
    DO dst_p = 1, numproc-1
       sm_disp(dst_p) = sm_disp(dst_p-1)&
            + metalen*total_send_count(imeta, dst_p-1)
       sd_disp(dst_p) = sd_disp(dst_p-1)&
            + total_send_count(idata, dst_p-1)
       rm_disp(dst_p) = rm_disp(dst_p-1)&
            + metalen*total_recv_count(imeta, dst_p-1)
       rd_disp(dst_p) = rd_disp(dst_p-1)&
            + total_recv_count(idata, dst_p-1)
    ENDDO
    sdp(:) = sd_disp
    smp(:) = sm_disp
    !
    CALL dbcsr_iterator_start (iter, sm)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk, blk_p=blk_p,&
            row_size=row_size, col_size=col_size)
       DO symmetry_i = 1, nsymmetries
          IF (symmetry_i .EQ. 1) THEN
             stored_row = row ; stored_col = col; tr = .FALSE.
          ELSE
             IF (row .EQ. col) CYCLE
             stored_row = col ; stored_col = row; tr = .TRUE.
          ENDIF
          ! Where do we send this block?
          prow = row_dist(stored_row)
          pcol = col_dist(stored_col)
          dst_p = blacs2mpi(prow, pcol)
          nze = row_size * col_size
          send_meta(smp(dst_p)) = stored_row
          send_meta(smp(dst_p)+1) = stored_col
          send_meta(smp(dst_p)+2) = sdp(dst_p) - sd_disp(dst_p) + 1
          IF (make_untr) THEN
             CALL dbcsr_block_partial_copy(dst=send_data,&
                  dst_rs=row_size, dst_cs=col_size,&
                  dst_tr=symmetry_i .GT. 1,&
                  dst_r_lb=1, dst_c_lb=1, dst_offset=sdp(dst_p)-1,&
                  nrow = row_size, ncol = col_size,&
                  src=sm%m%data_area,&
                  src_rs=row_size, src_cs=col_size,&
                  src_tr=blk_p .LT. 0,&
                  src_r_lb=1, src_c_lb=1,&
                  src_offset = ABS(blk_p)-1)
          ELSE
             CALL dbcsr_data_copy (send_data, (/ sdp(dst_p) /), (/ nze /),&
                  sm%m%data_area, (/ ABS(blk_p) /), (/ nze /), error=error)
             IF (tr) &
                  send_meta(smp(dst_p)+2) = -send_meta(smp(dst_p)+2)
          ENDIF
          smp(dst_p) = smp(dst_p) + metalen
          sdp(dst_p) = sdp(dst_p) + nze
       ENDDO ! symmetry_i
    ENDDO ! iter
    CALL dbcsr_iterator_stop (iter)
    ! Exchange the data and metadata structures.
    CALL hybrid_alltoall_any(send_data,&
         total_send_count(idata,:), sd_disp(:)-1,&
         recv_data,&
         total_recv_count(idata,:), rd_disp(:)-1, mp_obj, error=error)
    CALL mp_alltoall(send_meta(:), metalen*total_send_count(imeta,:), sm_disp(:)-1,&
         recv_meta(:), metalen*total_recv_count(imeta,:), rm_disp(:)-1,&
         mp_group)
    ! Now fill in the data.
    CALL dbcsr_work_create(desm,&
         SUM(recv_count(imeta,:)),&
         SUM(recv_count(idata,:)), n=1,&
         work_mutable=.FALSE.,&
         error=error)
    ! Switch data data area of the work matrix with the receieved data
    ! (avoids copying).
    CALL dbcsr_data_hold (recv_data)
    CALL dbcsr_data_release (desm%m%wms(1)%data_area)
    desm%m%wms(1)%data_area = recv_data
    !
    blk_ps = 1
    blks = 1
    ! WRITE(*,*)rm_disp
    ! WRITE(*,*)recv_count
    DO src_p = 0, numproc-1
       IF (careful_mod) THEN
          CALL dbcsr_assert ((blks-1)*3, "EQ", rm_disp(src_p)-1,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Count mismatch", __LINE__, error=error)
       ENDIF
       blks = (rm_disp(src_p)-1)/metalen+1
       DO blk_l = 1, recv_count(imeta, src_p)
          stored_row = recv_meta(rm_disp(src_p)+metalen*(blk_l-1))
          stored_col = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+1)
          blk_p =      recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+2)
          desm%m%wms(1)%row_i(blks) = stored_row
          desm%m%wms(1)%col_i(blks) = stored_col
          desm%m%wms(1)%blk_p(blks) = SIGN(ABS(blk_p) + (rd_disp(src_p)-1), blk_p)
          nze = row_blk_size(ABS(stored_row))&
               * col_blk_size(stored_col)
          blk_ps = blk_ps + nze
          blks = blks + 1
       ENDDO
    ENDDO
    !
    desm%m%wms(1)%lastblk = blks - 1
    desm%m%wms(1)%datasize = blk_ps - 1
    CALL dbcsr_finalize(desm, error=error)
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    CALL dbcsr_data_release (recv_data)
    DEALLOCATE(recv_meta)
    CALL dbcsr_data_release (send_data)
    DEALLOCATE(send_meta)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_desymmetrize_deep


! *****************************************************************************
!> \brief Distributes a matrix that is currently replicated.
!> \param[in,out] matrix      matrix to replicate
!> \param[in] fast            change just the index, don't touch the data
! *****************************************************************************
  SUBROUTINE dbcsr_distribute(matrix, fast, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    LOGICAL, INTENT(in), OPTIONAL            :: fast
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_distribute', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 2

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: c_dp
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: c_sp
    INTEGER                                  :: blk, col, error_handler, &
                                                mynode, nblks, nze, p, row
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size, &
                                                tmp_index
    LOGICAL                                  :: mini, tr
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_dp
    REAL(KIND=sp), DIMENSION(:), POINTER     :: r_sp
    TYPE(dbcsr_data_obj)                     :: tmp_data
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_obj)                          :: distributed

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_valid_index(matrix%m),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Matrix not initialized.",__LINE__,error)
    CALL dbcsr_assert (matrix%m%replication_type .NE. dbcsr_repl_none,&
         dbcsr_warning_level, dbcsr_caller_error, routineN,&
         "Distributing a non-replicated matrix makes no sense.",__LINE__,error)
    CALL dbcsr_access_flush (matrix, error=error)
    IF (PRESENT (fast)) THEN
       mini=fast
    ELSE
       mini=.FALSE.
    ENDIF
    SELECT CASE (matrix%m%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (matrix%m%data_area, r_dp)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (matrix%m%data_area, r_sp)
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
            routineN,"Only real double precision",__LINE__,error)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (matrix%m%data_area, c_dp)
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
            routineN,"Only real double precision",__LINE__,error)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (matrix%m%data_area, c_sp)
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
            routineN,"Only real double precision",__LINE__,error)
    END SELECT
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    dist = dbcsr_distribution (matrix)
    mp_obj = dbcsr_distribution_mp (dist)
    mynode = dbcsr_mp_mynode (dbcsr_distribution_mp (dist))
    !
    IF (mini) THEN
       ! We just mark the blocks as deleted.
       CALL dbcsr_iterator_start (iter, matrix)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, row, col, r_dp, tr, blk)
          tr = .FALSE.
          CALL dbcsr_get_stored_coordinates (matrix%m, row, col, tr, p)
          IF (mynode .EQ. p) THEN
             matrix%m%blk_p(blk) = 0
          ENDIF
       ENDDO
       CALL dbcsr_iterator_stop (iter)
       matrix%m%replication_type = dbcsr_repl_none
    ELSE
       CALL dbcsr_init(distributed)
       CALL dbcsr_create(distributed, name='Distributed '//TRIM(matrix%m%name),&
            template = matrix,&
            matrix_type = dbcsr_type_no_symmetry,&
            replication_type = dbcsr_repl_none,&
            error=error)
       distributed%m%replication_type = dbcsr_repl_none
       ! First count how many blocks are local.
       nze = 0
       nblks = 0
       CALL dbcsr_iterator_start (iter, matrix)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, row, col, r_dp, tr, blk)
          tr = .FALSE.
          CALL dbcsr_get_stored_coordinates (matrix%m, row, col, tr, p)
          IF (mynode .EQ. p) THEN
             nze = nze + row_blk_size(row) * col_blk_size(col)
             nblks = nblks + 1
          ENDIF
       ENDDO
       CALL dbcsr_iterator_stop (iter)
       ! Preallocate the erray
       CALL dbcsr_work_create(distributed, nblks_guess=nblks,&
            sizedata_guess=nze, work_mutable=.FALSE.,error=error)
       ! Now actually do the work
       CALL dbcsr_iterator_start (iter, matrix)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, row, col, r_dp, tr, blk)
          tr = .FALSE.
          CALL dbcsr_get_stored_coordinates (matrix%m, row, col, tr, p)
          IF (mynode .EQ. p) THEN
             CALL dbcsr_put_block (distributed, row, col, r_dp, tr)
          ENDIF
       ENDDO
       CALL dbcsr_iterator_stop (iter)
       CALL dbcsr_finalize (distributed, error=error)
       ! Now replace the data and index
       CALL dbcsr_switch_data_area (matrix, distributed%m%data_area,&
            previous_data_area=tmp_data, error=error)
       CALL dbcsr_switch_data_area (distributed, tmp_data, error=error)
       CALL dbcsr_data_release (tmp_data)
       tmp_index => matrix%m%index
       matrix%m%index => distributed%m%index
       distributed%m%index => tmp_index
       CALL dbcsr_repoint_index (matrix%m)
       matrix%m%nze = distributed%m%nze
       matrix%m%nblks = distributed%m%nblks
       CALL dbcsr_release (distributed)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_distribute


! *****************************************************************************
!> \brief Detransposes all blocks in a matrix
!> \param[in,out] matrix      DBCSR matrix
! *****************************************************************************
  SUBROUTINE dbcsr_make_untransposed_blocks(matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_make_untransposed_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_size, &
                                                error_handler, row, row_size
    INTEGER, DIMENSION(:), POINTER           :: cbs, rbs
    LOGICAL                                  :: sym_negation, tr
    TYPE(dbcsr_data_obj)                     :: block_data
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    rbs => array_data (dbcsr_row_block_sizes (matrix))
    cbs => array_data (dbcsr_col_block_sizes (matrix))
    sym_negation = matrix%m%negate_real

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(block_data,iter,row,col,tr,blk,row_size,col_size) &
!$OMP          SHARED(matrix,rbs,cbs,sym_negation,error)
    CALL dbcsr_data_init (block_data)
    CALL dbcsr_data_new (block_data, dbcsr_get_data_type(matrix))
    CALL dbcsr_iterator_start (iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, block_data,&
            transposed=tr,&
            block_number=blk)
       IF (tr) THEN
          row_size = rbs(row)
          col_size = cbs(col)
          CALL dbcsr_block_transpose(block_data, col_size, row_size, &
               error=error)
          IF (sym_negation) THEN
             SELECT CASE (block_data%d%data_type)
             CASE (dbcsr_type_real_4)
                block_data%d%r_sp(:) = -block_data%d%r_sp(:)
             CASE (dbcsr_type_real_8)
                block_data%d%r_dp(:) = -block_data%d%r_dp(:)
             CASE (dbcsr_type_complex_4)
                block_data%d%c_sp(:) = -block_data%d%c_sp(:)
             CASE (dbcsr_type_complex_8)
                block_data%d%c_dp(:) = -block_data%d%c_dp(:)
             END SELECT
          ENDIF
          matrix%m%blk_p(blk) = -matrix%m%blk_p(blk)
       ENDIF
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    CALL dbcsr_data_clear_pointer (block_data)
    CALL dbcsr_data_release (block_data)
!$OMP END PARALLEL

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_make_untransposed_blocks



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
    INTEGER, INTENT(IN)                      :: data_memory_type, &
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
    CALL dbcsr_access_flush (source, error=error)
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
!> \param[in] sm    input symmetric matrix
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
    INTEGER, INTENT(IN)                      :: data_memory_type, &
                                                index_memory_type
    LOGICAL, INTENT(IN), OPTIONAL            :: no_copy_data
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: scale_value
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_images', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 5
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    CHARACTER                                :: predist_type, predist_type_fwd
    INTEGER :: blk, blk_l, blk_p, bp, col, col_img, col_size, data_p, &
      data_type, dst_p, error_handler, error_handler2, i, icol, irow, &
      ithread, j, mp_group, ncol_images, nrow_images, nsymmetries, nthreads, &
      numproc, nze, pcol, prev_blk_p, prev_dst_p, prow, row, row_img, &
      row_size, sd_pos, sm_pos, src_p, stored_blk_p, stored_col, stored_row, &
      symmetry_i, t, vcol, vrow
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: lsdp, lsmp, rd_disp, &
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
    LOGICAL                                  :: nocopy, release_td, &
                                                same_dst_p, tr
    TYPE(dbcsr_data_obj)                     :: received_data_area, &
                                                recv_data_area, send_data_area
    TYPE(dbcsr_distribution_obj)             :: old_dist, target_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error, t_error
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_scalar_type)                  :: scale_neg_one
    TYPE(dbcsr_type)                         :: sm

!   ---------------------------------------------------------------------------

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
    row_blk_size => array_data (sm%row_blk_size)
    col_blk_size => array_data (sm%col_blk_size)
    target_dist = target_imgdist%i%main
    old_dist = dbcsr_distribution (ism)
    row_dist => array_data (dbcsr_distribution_row_dist (target_dist))
    col_dist => array_data (dbcsr_distribution_col_dist (target_dist))
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
    ALLOCATE (ums%mats(nrow_images,ncol_images))
    IF (sm%symmetry) THEN
       CALL dbcsr_assert(SIZE(row_dist),'EQ', SIZE(col_dist), dbcsr_warning_level,&
            dbcsr_wrong_args_error, routineN,&
            'Unequal row and column distributions for symmetric matrix.',__LINE__,error)
    ENDIF
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
    ALLOCATE (send_count(2, nrow_images, ncol_images, 0:numproc-1))
    ALLOCATE (recv_count(2, nrow_images, ncol_images, 0:numproc-1))
    ALLOCATE (total_send_count(2, 0:numproc-1))
    ALLOCATE (total_recv_count(2, 0:numproc-1))
    ALLOCATE (sdp(0:numproc-1))
    ALLOCATE (sd_disp(0:numproc-1))
    ALLOCATE (smp(0:numproc-1))
    ALLOCATE (sm_disp(0:numproc-1))
    ALLOCATE (rd_disp(0:numproc-1))
    ALLOCATE (rm_disp(0:numproc-1))
    ithread = 0
    nthreads = 1
    release_td = .FALSE.
!$  IF (.NOT. dbcsr_distribution_has_threads (old_dist)) THEN
!$     CALL dbcsr_distribution_make_threads (old_dist)
!$     release_td = .TRUE.
!$  ENDIF
!$  CALL dbcsr_assert (dbcsr_distribution_has_threads (old_dist),&
!$       dbcsr_fatal_level, dbcsr_internal_error, routineN,&
!$       "Thread distribution not defined", __LINE__, error=error)
    DO row_img = 1, nrow_images
       DO col_img = 1, ncol_images
          CALL dbcsr_init (ums%mats(row_img, col_img))
          CALL dbcsr_create(ums%mats(row_img, col_img), "imaged "//sm%name,&
               target_dist,&
               dbcsr_type_no_symmetry, sm%row_blk_size, sm%col_blk_size,&
               0,0, sm%data_type,&
               data_memory_type=data_memory_type,&
               index_memory_type=index_memory_type,&
               error=dbcsr_error)
          ums%mats(row_img, col_img)%m%negate_real = sm%negate_real
          ums%mats(row_img, col_img)%m%negate_imaginary = sm%negate_imaginary
          !ums%mats(row_img, col_img)%m%transpose = sm%transpose
       ENDDO
    ENDDO
    scale_neg_one = dbcsr_scalar_negative (dbcsr_scalar_one (sm%data_type))
!$omp parallel default (none) &
!$omp private (ithread, t_error, &
!$omp          row_img, col_img, iter,&
!$omp          myt_send_count, myt_total_send_count, &
!$omp          prev_dst_p, dst_p, same_dst_p, &
!$omp          row, col, blk, symmetry_i, stored_row, stored_col, &
!$omp          row_size, col_size, &
!$omp          prev_blk_p, blk_p, tr, data_p, stored_blk_p, &
!$omp          prow, pcol, vcol, vrow, i, j, nze, bp, sm_pos, sd_pos,&
!$omp          lsmp, lsdp, t) &
!$omp shared (nthreads, dbcsr_error, nocopy, release_td, &
!$omp         nrow_images, ncol_images, scale_value, error,&
!$omp         ums, sm, ism, target_dist, predistribute, predist_type, &
!$omp         predist_type_fwd,&
!$omp         data_memory_type, index_memory_type,&
!$omp         old_dist, &
!$omp         mp_obj, target_imgdist, mp_group, numproc, row_dist, col_dist,&
!$omp         row_img_dist, col_img_dist, blacs2mpi, row_blk_size,&
!$omp         col_blk_size, data_type, &
!$omp         send_count, recv_count, all_total_send_offset, total_send_count, &
!$omp         total_recv_count, &
!$omp         sd_disp, sm_disp, rd_disp, rm_disp, &
!$omp         send_meta, recv_meta, send_data_area, &
!$omp         recv_data_area, received_data_area, &
!$omp         blk_ps, blks, nsymmetries, scale_neg_one, error_handler2)
!$  ithread = omp_get_thread_num()
!$  IF (release_td) THEN !@@@
!$     CALL dbcsr_assert (.NOT.release_td, dbcsr_fatal_level, dbcsr_internal_error,&
!$     routineN, "No thread distribution defined", __LINE__, error=dbcsr_error)
!$     CALL dbcsr_distribution_make_threads (old_dist)
!$  ENDIF
    ! Create and allocate the imaged matrices.
!$omp master
!$  nthreads = omp_get_num_threads()
    ! C
!$omp end master
    ALLOCATE (myt_send_count(2, nrow_images, ncol_images, 0:numproc-1))
    ALLOCATE (myt_total_send_count(2, 0:numproc-1))
    myt_send_count(:,:,:,:) = 0
!$omp master
    ALLOCATE (all_total_send_offset(2, 0:numproc-1))
!$omp end master
    prev_dst_p = -1
    ! Count sizes for sending.
    CALL dbcsr_iterator_start(iter, ism, shared=.TRUE.)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk,&
            row_size=row_size, col_size=col_size)
       IF (row_size .EQ. 0 .OR. col_size .EQ. 0) CYCLE
       DO symmetry_i = 1, nsymmetries
          IF (symmetry_i .EQ. 1) THEN
             stored_row = row ; stored_col = col
          ELSE
             IF (row .EQ. col) CYCLE
             stored_row = col ; stored_col = row
          ENDIF
          ! Where do we send this block?
          row_img = 1
          col_img = 1
          IF (nrow_images .GT. 1) row_img = row_img_dist (stored_row)
          IF (ncol_images .GT. 1) col_img = col_img_dist (stored_col)
          CALL image_calculator(target_imgdist,&
               prow = prow, rowi = i,&
               pcol = pcol, coli = j,&
               vprow = vrow, vpcol = vcol,&
               myprow = row_dist(stored_row), myrowi = row_img,&
               mypcol = col_dist(stored_col), mycoli = col_img,&
               shifting = predist_type_fwd, error=error)
          row_img = i
          col_img = j
          dst_p = blacs2mpi(prow, pcol)
          same_dst_p = prev_dst_p .EQ. dst_p
          ! To allow normalization to non-transposed blocks, every
          ! index entry must have its own block.
          same_dst_p = .FALSE.
          prev_dst_p = dst_p
          ! These counts are meant for the thread that processes this row.
          myt_send_count(1, row_img, col_img, dst_p) =&
               myt_send_count(1, row_img, col_img, dst_p) + 1
          ! Data can be duplicated if the transpose is destined to the same
          ! process.
          IF (.NOT. same_dst_p .OR. symmetry_i .EQ. 1) THEN
             nze = row_blk_size(stored_row) * col_blk_size(stored_col)
             myt_send_count(2, row_img, col_img, dst_p) =&
                  myt_send_count(2, row_img, col_img, dst_p) + nze
          ENDIF
       ENDDO ! symmetry_i
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    FORALL (dst_p = 0:numproc-1)
       myt_total_send_count(1, dst_p) = SUM(myt_send_count(1,:,:,dst_p))
       myt_total_send_count(2, dst_p) = SUM(myt_send_count(2,:,:,dst_p))
    END FORALL
    ! Merge the send counts
!$omp master
    send_count(:,:,:,:) = 0
!$omp end master
!$omp barrier
!$omp critical
    send_count(:,:,:,:) = send_count(:,:,:,:) + myt_send_count(:,:,:,:)
!$omp end critical
    DEALLOCATE (myt_send_count)
!$omp barrier
!$omp master
    CALL dbcsr_error_set(routineN//"_sizes", error_handler2, dbcsr_error)
    CALL mp_alltoall(send_count, recv_count, 2*nrow_images*ncol_images,&
         mp_group)
    CALL dbcsr_error_stop(error_handler2, dbcsr_error)
    ! Fill in the meta data structures and copy the data.
    DO dst_p = 0, numproc-1
       total_send_count(1, dst_p) = SUM (send_count (1, :, :, dst_p))
       total_send_count(2, dst_p) = SUM (send_count (2, :, :, dst_p))
       total_recv_count(1, dst_p) = SUM (recv_count (1, :, :, dst_p))
       total_recv_count(2, dst_p) = SUM (recv_count (2, :, :, dst_p))
    ENDDO
    ! Allocate data structures needed for data exchange.
    CALL dbcsr_data_init (recv_data_area)
    CALL dbcsr_data_new (recv_data_area, data_type, SUM(recv_count(2, :, :, :)))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :, :, :))))
    CALL dbcsr_data_init (send_data_area)
    CALL dbcsr_data_new (send_data_area, data_type, SUM(send_count(2, :, :, :)))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :, :, :))))
    ! Calculate displacements for processors needed for the exchanges.
    sd_disp = -1 ; sm_disp = -1
    rd_disp = -1 ; rm_disp = -1
    sd_disp(0) = 1 ; sm_disp(0) = 1
    rd_disp(0) = 1 ; rm_disp(0) = 1
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
!$omp end master
!$omp barrier
    ! Thread-local pointers of the current adding position into the
    ! send buffers
    ALLOCATE (lsmp(0:numproc-1), lsdp(0:numproc-1))
    ! Calculate thread-local displacemnts
    IF (ithread .EQ. 0) THEN
       lsmp(:) = sm_disp(:)
       lsdp(:) = sd_disp(:)
       IF (nthreads .GT. 1) THEN
          all_total_send_offset(1,:) = sm_disp(:) + metalen*myt_total_send_count(1,:)
          all_total_send_offset(2,:) = sd_disp(:) + myt_total_send_count(2,:)
       ENDIF
    ENDIF
!$omp barrier
    IF (ithread .GT. 0) THEN
!$omp critical
       lsmp(:) = all_total_send_offset(1,:)
       lsdp(:) = all_total_send_offset(2,:)
       all_total_send_offset(1,:) &
            = all_total_send_offset(1,:) + metalen*myt_total_send_count(1,:)
       all_total_send_offset(2,:) &
            = all_total_send_offset(2,:) + myt_total_send_count(2,:)
!$omp end critical
    ENDIF
    DEALLOCATE (myt_total_send_count)
    ! Prepares some indices needed for the last DO loop that copies
    ! from buffer to local space. Placed early to take advantage of
    ! the SECTIONS.
!$omp master
    ALLOCATE (blk_ps(nrow_images, ncol_images))
    ALLOCATE (blks (nrow_images, ncol_images))
    blk_ps(:,:) = 1
    blks(:,:) = 1
    ! Prepares the work matrices used in the last DO loop. Placed
    ! early.
    CALL dbcsr_data_init (received_data_area)
    received_data_area = recv_data_area
    CALL dbcsr_data_hold(received_data_area)
!$omp end master
!$omp barrier
!$omp master
    DEALLOCATE (all_total_send_offset)
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
!$omp end master
!$omp barrier
    ! add timing call to the packing of the send buffers 
    CALL dbcsr_error_set(routineN//"_pack", error_handler2, dbcsr_error)
    prev_dst_p = -1
    ! Copies metadata and actual data to be sent into the send buffers.
    CALL dbcsr_iterator_start(iter, ism, shared=.TRUE.)
    prev_blk_p = 0
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk, blk_p=blk_p,&
            row_size=row_size, col_size=col_size)
       IF (row_size .EQ. 0 .OR. col_size .EQ. 0) CYCLE
       bp = ABS(blk_p)
       DO symmetry_i = 1, nsymmetries
          IF (symmetry_i .EQ. 1) THEN
             stored_row = row ; stored_col = col; tr = blk_p .LT. 0
          ELSE
             IF (row .EQ. col) CYCLE
             stored_row = col ; stored_col = row; tr = blk_p .GT. 0
          ENDIF
          ! Where do we send this block?
          row_img = 1
          col_img = 1
          IF (nrow_images .GT. 1) row_img = row_img_dist (stored_row)
          IF (ncol_images .GT. 1) col_img = col_img_dist (stored_col)
          CALL image_calculator(target_imgdist,&
               prow = prow, rowi = i,&
               pcol = pcol, coli = j,&
               vprow = vrow, vpcol = vcol,&
               myprow = row_dist(stored_row), myrowi = row_img,&
               mypcol = col_dist(stored_col), mycoli = col_img,&
               shifting = predist_type_fwd, error=error)
          row_img = i
          col_img = j
          dst_p = blacs2mpi(prow, pcol)
          same_dst_p = dst_p .EQ. prev_dst_p
          ! To allow normalization to non-transposed blocks, every
          ! index entry must have its own block.
          same_dst_p = .FALSE.
          prev_dst_p = dst_p
          sm_pos = lsmp(dst_p)
          lsmp(dst_p) = lsmp(dst_p) + metalen
          send_meta(sm_pos) = stored_row
          IF (.NOT. same_dst_p ) THEN
           nze = row_blk_size(stored_row) * col_blk_size(stored_col)
           sd_pos = lsdp(dst_p)
           lsdp(dst_p) = lsdp(dst_p) + nze
           IF ( symmetry_i .GT. 1 ) THEN
              ! transpose variant
             SELECT CASE (send_data_area%d%data_type)
              CASE (dbcsr_type_real_4)
                CALL dbcsr_block_transpose(send_data_area%d%r_sp(sd_pos:lsdp(dst_p)-1),&
                sm%data_area%d%r_sp(ABS(blk_p):ABS(blk_p)+nze-1), row_size, col_size)
              CASE (dbcsr_type_real_8)
                CALL dbcsr_block_transpose(send_data_area%d%r_dp(sd_pos:lsdp(dst_p)-1),&
                sm%data_area%d%r_dp(ABS(blk_p):ABS(blk_p)+nze-1), row_size, col_size)
              CASE (dbcsr_type_complex_4)
                CALL dbcsr_block_transpose(send_data_area%d%c_sp(sd_pos:lsdp(dst_p)-1),&
                sm%data_area%d%c_sp(ABS(blk_p):ABS(blk_p)+nze-1), row_size, col_size)
              CASE (dbcsr_type_complex_8)
                CALL dbcsr_block_transpose(send_data_area%d%c_dp(sd_pos:lsdp(dst_p)-1),&
                sm%data_area%d%c_dp(ABS(blk_p):ABS(blk_p)+nze-1), row_size, col_size)
              CASE default
                CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
                   routineN, "Invalid data type.",__LINE__,error)
             END SELECT
           ELSE   ! normal variant
             SELECT CASE (send_data_area%d%data_type)
              CASE (dbcsr_type_real_4)
                IF ( blk_p .LT. 0 ) THEN
                  DO icol=0,col_size -1 ; DO irow=0,row_size-1
                   send_data_area%d%r_sp(sd_pos+irow + icol*row_size) &
                   = sm%data_area%d%r_sp(bp+icol+(irow*col_size))
                  ENDDO ; ENDDO
                ELSE
                  DO icol=0,col_size -1 ; DO irow=0,row_size-1
                   send_data_area%d%r_sp(sd_pos+irow + icol*row_size) &
                    = sm%data_area%d%r_sp(bp+irow+icol*row_size)
                  ENDDO ; ENDDO
                END IF
              CASE (dbcsr_type_real_8)
                IF ( blk_p .LT. 0 ) THEN
                  DO icol=0,col_size -1 ; DO irow=0,row_size-1
                   send_data_area%d%r_dp(sd_pos+irow + icol*row_size) &
                   = sm%data_area%d%r_dp(bp+icol+(irow*col_size))
                  ENDDO ; ENDDO
                ELSE
                  DO icol=0,col_size -1 ; DO irow=0,row_size-1
                   send_data_area%d%r_dp(sd_pos+irow + icol*row_size) &
                    = sm%data_area%d%r_dp(bp+irow+icol*row_size)
                  ENDDO ; ENDDO
                END IF
              CASE (dbcsr_type_complex_4)
                IF ( blk_p .LT. 0 ) THEN
                  DO icol=0,col_size -1 ; DO irow=0,row_size-1
                   send_data_area%d%c_sp(sd_pos+irow + icol*row_size) &
                   = sm%data_area%d%c_sp(bp+icol+(irow*col_size))
                  ENDDO ; ENDDO
                ELSE
                  DO icol=0,col_size -1 ; DO irow=0,row_size-1
                   send_data_area%d%c_sp(sd_pos+irow + icol*row_size) &
                    = sm%data_area%d%c_sp(bp+irow+icol*row_size)
                  ENDDO ; ENDDO
                END IF
              CASE (dbcsr_type_complex_8)
                IF ( blk_p .LT. 0 ) THEN
                  DO icol=0,col_size -1 ; DO irow=0,row_size-1
                   send_data_area%d%c_dp(sd_pos+irow + icol*row_size) &
                   = sm%data_area%d%c_dp(bp+icol+(irow*col_size))
                  ENDDO ; ENDDO
                ELSE
                  DO icol=0,col_size -1 ; DO irow=0,row_size-1
                   send_data_area%d%c_dp(sd_pos+irow + icol*row_size) &
                    = sm%data_area%d%c_dp(bp+irow+icol*row_size)
                  ENDDO ; ENDDO
                END IF
               CASE DEFAULT
                CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
                   routineN, "Invalid data type.",__LINE__,error)
              END SELECT
           ENDIF
             IF (PRESENT (scale_value)) THEN
                CALL dbcsr_block_scale (send_data_area, scale=scale_value,&
                     row_size=nze, col_size=1, lb=sd_pos, error=error)
             ENDIF
             IF (tr) THEN
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
             ENDIF
             send_meta(sm_pos+1) = stored_col
             send_meta(sm_pos+2) = sd_pos-sd_disp(dst_p)+1
             send_meta(sm_pos+3) = row_img
             send_meta(sm_pos+4) = col_img
             prev_blk_p = send_meta(sm_pos+2)
           ELSE
             send_meta(sm_pos+1) = -stored_col
             send_meta(sm_pos+2) = -prev_blk_p
             send_meta(sm_pos+3) = row_img
             send_meta(sm_pos+4) = col_img
          ENDIF
       ENDDO ! symmetry_i
    ENDDO ! iterator
    CALL dbcsr_iterator_stop(iter)
    DEALLOCATE (lsmp, lsdp)
    CALL dbcsr_error_stop(error_handler2, dbcsr_error)

!$omp end parallel
    ! Exchange the data and metadata structures. In the interesting cases (square grids, row col distribution same),
    ! there are only very few processors that need to exchange data.
    ! The hybrid_alltoall deals with this by doing point to point communication
    CALL dbcsr_error_set(routineN//"_data", error_handler2, dbcsr_error)
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       CALL hybrid_alltoall_s1(&
            send_data_area%d%r_sp(:), total_send_count(2,:), sd_disp(:)-1,&
            recv_data_area%d%r_sp(:), total_recv_count(2,:), rd_disp(:)-1,&
            most_ptp=.TRUE., remainder_ptp=.TRUE., no_hybrid=.FALSE.,&
            mp_env = mp_obj)
    CASE (dbcsr_type_real_8)
       !CALL mp_alltoall(&
       !     send_data_area%d%r_dp(:), total_send_count(2,:), sd_disp(:)-1,&
       !     recv_data_area%d%r_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
       !     mp_group)
       CALL hybrid_alltoall_d1 (&
            send_data_area%d%r_dp(:), total_send_count(2,:), sd_disp(:)-1,&
            recv_data_area%d%r_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
            most_ptp=.TRUE., remainder_ptp=.TRUE., no_hybrid=.FALSE.,&
            mp_env = mp_obj)
    CASE (dbcsr_type_complex_4)
       CALL hybrid_alltoall_c1(&
            send_data_area%d%c_sp(:), total_send_count(2,:), sd_disp(:)-1,&
            recv_data_area%d%c_sp(:), total_recv_count(2,:), rd_disp(:)-1,&
            most_ptp=.TRUE., remainder_ptp=.TRUE., no_hybrid=.FALSE.,&
            mp_env = mp_obj)
    CASE (dbcsr_type_complex_8)
       CALL hybrid_alltoall_z1(&
            send_data_area%d%c_dp(:), total_send_count(2,:), sd_disp(:)-1,&
            recv_data_area%d%c_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
            most_ptp=.TRUE., remainder_ptp=.TRUE., no_hybrid=.FALSE.,&
            mp_env = mp_obj)
    END SELECT
    CALL hybrid_alltoall_i1(&
         send_meta(:), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), metalen*total_recv_count(1,:), rm_disp(:)-1,&
         most_ptp=.TRUE., remainder_ptp=.TRUE., no_hybrid=.FALSE.,&
         mp_env = mp_obj)
    CALL dbcsr_error_stop(error_handler2, dbcsr_error)

    ! Now create the work index and/or copy the relevant data from the
    ! receive buffer into the local indices.
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


    ! Finalize the actual imaged matrices from the work matrices.
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

    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_meta)
    CALL dbcsr_data_release (send_data_area)
    !DEALLOCATE(send_data)
    DEALLOCATE(send_meta)
    ! Get rid of the cuckoo.
!$  IF (release_td) THEN
!$     CALL dbcsr_distribution_no_threads (old_dist)
!$  ENDIF
    CALL dbcsr_data_release (received_data_area)
    !
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
          CALL make_dense (tmp_mat, images%mats(mat_row,mat_col),&
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
!> \brief Makes a dense matrix, inplace.
!> \param[in,out] matrix      matrix to make dense
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE dbcsr_make_dense (matrix, dense_matrix, dense_dist,&
       dense_row_sizes, dense_col_sizes, row_map, col_map, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_obj), INTENT(OUT)             :: dense_matrix
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dense_dist
    TYPE(array_i1d_obj), INTENT(IN)          :: dense_row_sizes, &
                                                dense_col_sizes, row_map, &
                                                col_map
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_dense', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handler
    REAL(kind=dp)                            :: cs
    TYPE(array_i1d_obj)                      :: dense_local_cols, &
                                                dense_local_rows
    TYPE(dbcsr_distribution_obj)             :: old_distribution

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_init (dense_matrix)
    CALL dbcsr_create (dense_matrix, template=matrix,&
         dist = dense_dist,&
         row_blk_size=dense_row_sizes, col_blk_size=dense_col_sizes,&
         error=error)
    !
    IF (dbg) THEN
       cs = dbcsr_checksum (matrix, error=error)
       WRITE(*,*)routineN//" prod cs pre", cs
    ENDIF
    CALL dbcsr_access_flush (matrix, error=error)
    old_distribution = dbcsr_distribution (matrix)
    ! Conversion of global to local offsets for the dense blocks
    CALL dbcsr_get_local_rows (dense_dist, dense_local_rows,&
         dense_matrix%m%index(dbcsr_slot_home_prow), error=error)
    CALL dbcsr_get_local_cols (dense_dist, dense_local_cols,&
         dense_matrix%m%index(dbcsr_slot_home_pcol), error=error)
    !
    CALL make_dense (matrix, dense_matrix,&
         array_data(dbcsr_distribution_local_rows (old_distribution)),&
         array_data(dbcsr_distribution_local_cols (old_distribution)),&
         array_data(matrix%m%row_blk_offset),&
         array_data(matrix%m%col_blk_offset),&
         array_data(dense_local_rows), array_data(dense_local_cols),&
         array_data(dense_matrix%m%row_blk_offset),&
         array_data(dense_matrix%m%col_blk_offset),&
         array_data(row_map), array_data(col_map), .TRUE., .TRUE., error=error)
    IF (dbg) THEN
       cs = dbcsr_checksum (dense_matrix, error=error)
       WRITE(*,*)routineN//" prod cs pst", cs
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_make_dense

! *****************************************************************************
!> \brief Copies a matrix and makes its data dense.
!> \note the dense_matrix must have been created
!> \param[in] und_matrix         Original non-dense matrix
!> \param[in,out] dense_matrix   Dense copy of und_matrix
!> \param[in] und_local_rows     The process-grid local rows of the non-dense
!>                               und_matrix
!> \param[in] und_local_cols     The process-grid local columns of the
!>                               non-dense und_matrix
!> \param[in] und_row_blk_offsets  The block offsets of the rows of the
!>                                 non-dense matrix
!> \param[in] und_col_blk_offsets  The block offsets of the columns of the
!>                                non-dense matrix
!> \param[in] dense_local_rows    The process-grid local rows of the dense
!>                                matrix
!> \param[in] dense_local_cols    The process-grid local columns of the
!>                                dense matrix
!> \param[in] dense_row_blk_offsets  The block offsets of the rows of the
!>                                   dense matrix
!> \param[in] dense_col_blk_offsets  The block offsets of the columns of the
!>                                   dense matrix
!> \param[in] row_map                Mapping of non-dense rows to dense rows
!> \param[in] col_map                Mapping of non-dense columns to dense
!>                                   columns
!> \param[in] join_rows              Make rows dense
!> \param[in] join_cols              Make columns dense
!> \param[in] error                  Error
! *****************************************************************************
  SUBROUTINE make_dense (und_matrix, dense_matrix,&
       und_local_rows, und_local_cols,&
       und_row_blk_offsets, und_col_blk_offsets,&
       dense_local_rows, dense_local_cols,&
       dense_row_blk_offsets, dense_col_blk_offsets,&
       row_map, col_map, join_rows, join_cols,&
       error)
    TYPE(dbcsr_obj), INTENT(IN)              :: und_matrix
    TYPE(dbcsr_obj), INTENT(INOUT)           :: dense_matrix
    INTEGER, DIMENSION(:), INTENT(IN) :: und_local_rows, und_local_cols, &
      und_row_blk_offsets, und_col_blk_offsets, dense_local_rows, &
      dense_local_cols, dense_row_blk_offsets, dense_col_blk_offsets, &
      row_map, col_map
    LOGICAL, INTENT(IN)                      :: join_rows, join_cols
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_dense', &
      routineP = moduleN//':'//routineN

    INTEGER :: dense_nblkcols_local, dense_nblkcols_total, &
      dense_nblkrows_local, dense_nblkrows_total, dense_nlocal_blocks, &
      error_handle, nfullcols, nfullrows, und_nblkcols_total, &
      und_nblkrows_total
    INTEGER, ALLOCATABLE, DIMENSION(:) :: col_internal_offsets, &
      dense_local_col_blk_offsets, dense_local_row_blk_offsets, &
      row_internal_offsets, und_local_col_blk_offsets, &
      und_local_row_blk_offsets
    INTEGER, DIMENSION(dbcsr_meta_size)      :: meta
    TYPE(dbcsr_data_obj)                     :: dense_data, und_data

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handle, error)
    !
    nfullrows = nfull_elements (und_row_blk_offsets, und_local_rows)
    nfullcols = nfull_elements (und_col_blk_offsets, und_local_cols)
    !
    und_nblkrows_total = SIZE (und_row_blk_offsets)-1
    und_nblkcols_total = SIZE (und_col_blk_offsets)-1
    !
    ! Find the local data offsets (but indexed by the global
    ! rows/columns) for the undense data.
    ALLOCATE (und_local_row_blk_offsets(und_nblkrows_total+1))
    ALLOCATE (und_local_col_blk_offsets(und_nblkcols_total+1))
    CALL global_offsets_to_local (&
         und_row_blk_offsets,&
         und_local_rows,&
         und_local_row_blk_offsets)
    CALL global_offsets_to_local (&
         und_col_blk_offsets,&
         und_local_cols,&
         und_local_col_blk_offsets)
    !
    dense_nblkrows_total = SIZE (dense_row_blk_offsets)-1
    dense_nblkcols_total = SIZE (dense_col_blk_offsets)-1
    dense_nblkrows_local = SIZE (dense_local_rows)
    dense_nblkcols_local = SIZE (dense_local_cols)
    !
    ! Find the local data offsets (but indexed by the (dense) global
    ! rows/columns) for the dense data.
    ALLOCATE (dense_local_row_blk_offsets (dense_nblkrows_total+1))
    ALLOCATE (dense_local_col_blk_offsets (dense_nblkcols_total+1))
    CALL global_offsets_to_local (&
         dense_row_blk_offsets,&
         dense_local_rows,&
         dense_local_row_blk_offsets)
    CALL global_offsets_to_local (&
         dense_col_blk_offsets,&
         dense_local_cols,&
         dense_local_col_blk_offsets)
    ! Find the offset of blocks within dense rows/columns.  This is needed
    ! since the blocked rows/columns are not necessarily in the same order.
    ALLOCATE (row_internal_offsets (und_nblkrows_total))
    ALLOCATE (col_internal_offsets (und_nblkcols_total))
    CALL get_internal_offsets (&
         und_local_rows, row_map,&
         und_local_row_blk_offsets,&
         dense_local_row_blk_offsets,&
         row_internal_offsets)
    CALL get_internal_offsets (&
         und_local_cols, col_map,&
         und_local_col_blk_offsets,&
         dense_local_col_blk_offsets,&
         col_internal_offsets)
    !
    und_data = und_matrix%m%data_area
    CALL dbcsr_data_hold (und_data)
    CALL dbcsr_data_init (dense_data)
    CALL dbcsr_data_new (dense_data, dbcsr_data_get_type (und_data),&
         data_size = nfullrows * nfullcols,&
         memory_type=dbcsr_data_get_memory_type (und_data))
    !
    ! Reshuffle the data
    CALL make_dense_data (und_matrix,&
         dense_data, nfullrows, nfullcols,&
         und_local_row_blk_offsets, und_local_col_blk_offsets,&
         dense_local_row_blk_offsets, dense_local_col_blk_offsets,&
         row_map = row_map, col_map = col_map,&
         row_internal_offsets = row_internal_offsets,&
         col_internal_offsets = col_internal_offsets,&
         join_rows=join_rows, join_cols=join_cols, make_tr=.FALSE.)
    CALL dbcsr_switch_data_area (dense_matrix, dense_data, error=error)
    CALL dbcsr_data_release (dense_data)
    CALL dbcsr_data_release (und_data)
    !
    ! Create the new dense index.
    dense_nlocal_blocks = dense_nblkrows_local * dense_nblkcols_local
    CALL dbcsr_addto_index_array (dense_matrix%m,&
         dbcsr_slot_row_p,&
         reservation=dense_nblkrows_total+1, extra=2*dense_nlocal_blocks, error=error)
    CALL dbcsr_addto_index_array (dense_matrix%m,&
         dbcsr_slot_col_i,&
         reservation=dense_nlocal_blocks, error=error)
    CALL dbcsr_addto_index_array (dense_matrix%m,&
         dbcsr_slot_blk_p,&
         reservation=dense_nlocal_blocks, error=error)
    !
    meta = dense_matrix%m%index(1:dbcsr_meta_size)
    CALL dbcsr_pack_meta (dense_matrix%m, meta)
    meta(dbcsr_slot_nze) = nfullrows*nfullcols
    meta(dbcsr_slot_nblks) = dense_nlocal_blocks
    CALL make_dense_index (dense_matrix%m%row_p,&
         dense_matrix%m%col_i,&
         dense_matrix%m%blk_p,&
         dense_nblkrows_total, dense_nblkcols_total,&
         dense_local_rows,&
         dense_local_cols,&
         dense_local_row_blk_offsets,&
         dense_local_col_blk_offsets,&
         make_tr=.FALSE.,&
         meta=meta, error=error)
    CALL dbcsr_unpack_meta (dense_matrix%m, meta)
    !
    DEALLOCATE (und_local_row_blk_offsets)
    DEALLOCATE (und_local_col_blk_offsets)
    DEALLOCATE (dense_local_row_blk_offsets)
    DEALLOCATE (dense_local_col_blk_offsets)
    DEALLOCATE (row_internal_offsets)
    DEALLOCATE (col_internal_offsets)
    !
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE make_dense


! *****************************************************************************
!> \brief Makes a blocked matrix from a dense matrix, inplace
!> \param[in] matrix                 dense matrix
!> \param[in,out] undense_matrix     matrix to make undense
!> \param[in] distribution    distribution of non-dense rows and columns
!> \param[in] row_blk_offsets, col_blk_offsets   non-dense row and column block offsets
!> \param[in] row_blk_sizes, col_blk_sizes       non-dense row and column block sizes
!> \param[in] row_map, col_map   mapping from non-dense rows/columns
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE dbcsr_make_undense (matrix, undense_matrix, distribution,&
       row_blk_offsets, col_blk_offsets, row_blk_sizes, col_blk_sizes,&
       row_map, col_map, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_obj), INTENT(INOUT)           :: undense_matrix
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: distribution
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_offsets, &
                                                col_blk_offsets, &
                                                row_blk_sizes, col_blk_sizes, &
                                                row_map, col_map
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_undense', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER :: error_handler, nblkcols_local, nblkcols_total, nblkrows_local, &
      nblkrows_total, nfullcols_local, nfullrows_local
    INTEGER, ALLOCATABLE, DIMENSION(:) :: col_internal_offsets, &
      dense_local_col_blk_offsets, dense_local_row_blk_offsets, &
      local_col_blk_offsets, local_row_blk_offsets, row_internal_offsets
    INTEGER, DIMENSION(:), POINTER           :: local_cols, local_rows, meta
    REAL(kind=dp)                            :: cs
    TYPE(dbcsr_data_obj)                     :: blocked_data, dense_data
    TYPE(dbcsr_distribution_obj)             :: dense_distribution

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_access_flush (matrix, error=error)
    IF (dbg) THEN
       cs = dbcsr_checksum (matrix, error=error)
       WRITE(*,*)routineN//" prod cs pre", cs
    ENDIF
    dense_distribution = dbcsr_distribution (matrix)
    nfullrows_local = matrix%m%nfullrows_local
    nfullcols_local = matrix%m%nfullcols_local
    nblkrows_local = dbcsr_distribution_nlocal_rows (distribution)
    nblkcols_local = dbcsr_distribution_nlocal_cols (distribution)
    nblkrows_total = dbcsr_distribution_nrows (distribution)
    nblkcols_total = dbcsr_distribution_ncols (distribution)
    local_rows => array_data (dbcsr_distribution_local_rows (distribution))
    local_cols => array_data (dbcsr_distribution_local_cols (distribution))
    CALL dbcsr_create (undense_matrix, template=matrix,&
         dist = distribution,&
         row_blk_size = row_blk_sizes,&
         col_blk_size = col_blk_sizes,&
         error=error)
    ! Restore previous offsets, just to try to keep the same memory.
    CALL array_release (undense_matrix%m%row_blk_offset)
    CALL array_release (undense_matrix%m%col_blk_offset)
    undense_matrix%m%row_blk_offset = row_blk_offsets
    undense_matrix%m%col_blk_offset = col_blk_offsets
    CALL array_hold (undense_matrix%m%row_blk_offset)
    CALL array_hold (undense_matrix%m%col_blk_offset)
    !
    ALLOCATE (local_row_blk_offsets(nblkrows_total+1))
    ALLOCATE (local_col_blk_offsets(nblkcols_total+1))
    CALL dbcsr_clearfrom_index_array (undense_matrix%m, dbcsr_slot_row_p)
    CALL dbcsr_clearfrom_index_array (undense_matrix%m, dbcsr_slot_col_i)
    CALL dbcsr_clearfrom_index_array (undense_matrix%m, dbcsr_slot_blk_p)
    CALL dbcsr_addto_index_array (undense_matrix%m, dbcsr_slot_row_p,&
         reservation=nblkrows_total+1, extra=nblkrows_local*nblkcols_local*2,&
         error=error)
    CALL dbcsr_addto_index_array (undense_matrix%m, dbcsr_slot_col_i,&
         reservation=nblkrows_local*nblkcols_local, error=error)
    CALL dbcsr_addto_index_array (undense_matrix%m, dbcsr_slot_blk_p,&
         reservation=nblkrows_local*nblkcols_local, error=error)
    meta => undense_matrix%m%index(1:dbcsr_meta_size)
    CALL dbcsr_pack_meta (undense_matrix%m, meta)
    meta(dbcsr_slot_nblks) = nblkrows_local*nblkcols_local
    meta(dbcsr_slot_nze) = nfullrows_local * nfullcols_local
    CALL global_offsets_to_local (array_data (row_blk_offsets),&
         local_rows, local_row_blk_offsets(1:nblkrows_local+1))
    CALL global_offsets_to_local (array_data (col_blk_offsets),&
         local_cols, local_col_blk_offsets(1:nblkcols_local+1))
    CALL make_undense_index (undense_matrix%m%row_p, undense_matrix%m%col_i, undense_matrix%m%blk_p,&
         distribution,&
         local_row_blk_offsets(1:nblkrows_local+1),&
         local_col_blk_offsets(1:nblkcols_local+1),&
         meta)
    CALL dbcsr_unpack_meta (undense_matrix%m, meta)
    !
    CALL global_offsets_to_local (array_data (row_blk_offsets),&
         local_rows, local_row_blk_offsets)
    CALL global_offsets_to_local (array_data (col_blk_offsets),&
         local_cols, local_col_blk_offsets)
    !
    ALLOCATE (dense_local_row_blk_offsets (1+dbcsr_distribution_nrows (dense_distribution)))
    ALLOCATE (dense_local_col_blk_offsets (1+dbcsr_distribution_ncols (dense_distribution)))
    CALL global_offsets_to_local (array_data(matrix%m%row_blk_offset),&
         array_data(dbcsr_distribution_local_rows (dense_distribution)),&
         dense_local_row_blk_offsets)
    CALL global_offsets_to_local (array_data(matrix%m%col_blk_offset),&
         array_data(dbcsr_distribution_local_cols (dense_distribution)),&
         dense_local_col_blk_offsets)
    ! Find the offset of blocks within dense rows/columns.  This is needed
    ! since the blocked rows/columns are not necessarily in the same order.
    ALLOCATE (row_internal_offsets (nblkrows_total))
    ALLOCATE (col_internal_offsets (nblkcols_total))
    CALL get_internal_offsets (&
         local_rows, array_data(row_map),&
         local_row_blk_offsets,&
         dense_local_row_blk_offsets,&
         row_internal_offsets)
    CALL get_internal_offsets (&
         local_cols, array_data(col_map),&
         local_col_blk_offsets,&
         dense_local_col_blk_offsets,&
         col_internal_offsets)
    !
    dense_data = matrix%m%data_area
    CALL dbcsr_data_hold (dense_data)
    CALL dbcsr_data_init (blocked_data)
    CALL dbcsr_data_new (blocked_data, dbcsr_data_get_type (dense_data),&
         data_size = nfullrows_local*nfullcols_local,&
         memory_type=dbcsr_data_get_memory_type (dense_data))
    CALL dbcsr_switch_data_area (undense_matrix, blocked_data, error=error)
    CALL dbcsr_data_release (blocked_data)
    ! Reshuffle the data
    CALL make_undense_data (undense_matrix, dense_data,&
         nfullrows_local, nfullcols_local,&
         local_row_blk_offsets, local_col_blk_offsets,&
         dense_local_row_blk_offsets, dense_local_col_blk_offsets,&
         array_data(row_map), array_data(col_map),&
         row_internal_offsets, col_internal_offsets, error)
    CALL dbcsr_data_release (dense_data)
    IF (dbg) THEN
       cs = dbcsr_checksum (matrix, error=error)
       WRITE(*,*)routineN//" prod cs pst", cs
    ENDIF
    DEALLOCATE (local_row_blk_offsets)
    DEALLOCATE (local_col_blk_offsets)
    DEALLOCATE (dense_local_row_blk_offsets)
    DEALLOCATE (dense_local_col_blk_offsets)
    DEALLOCATE (row_internal_offsets)
    DEALLOCATE (col_internal_offsets)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_make_undense

! *****************************************************************************
!> \brief Shuffles the data from blocked to standard dense form
!> \param[in] matrix                 Existing blocked matrix
!> \param[in,out] dense_data         Dense data
!> \param[in] nfullrows, nfullcols    size of new data
!> \param[in] join_cols           (optional) make columns dense, default is
!>                                yes
!> \param[in] join_rows           (optional) make rows dense, default is yes
!> \param[in] make_tr             (optional) make the dense blocks transposed
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE make_dense_data (matrix, dense_data, nfullrows, nfullcols,&
       und_row_blk_offsets, und_col_blk_offsets,&
       dense_row_blk_offsets, dense_col_blk_offsets,&
       row_map, col_map,&
       row_internal_offsets, col_internal_offsets,&
       join_rows, join_cols, make_tr)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dense_data
    INTEGER, INTENT(IN)                      :: nfullrows, nfullcols
    INTEGER, DIMENSION(:), INTENT(IN) :: und_row_blk_offsets, &
      und_col_blk_offsets, dense_row_blk_offsets, dense_col_blk_offsets, &
      row_map, col_map, row_internal_offsets, col_internal_offsets
    LOGICAL, INTENT(IN)                      :: join_rows, join_cols, make_tr

    CHARACTER(len=*), PARAMETER :: routineN = 'make_dense_data', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk_col, blk_col_offset, blk_col_size, blk_row, &
      blk_row_offset, blk_row_size, dense_col, dense_row, error_handle, &
      target_col_offset, target_cs, target_offset, target_row_offset, &
      target_rs, tco, tro
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: block
    TYPE(dbcsr_error_type)                   :: error
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    CALL dbcsr_assert (dbcsr_data_get_size (dense_data), 'GE',&
         nfullrows*nfullcols,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Dense data too small",__LINE__,error)
    CALL dbcsr_assert (join_cols, "OR", join_rows, dbcsr_warning_level,&
         dbcsr_unimplemented_error_nr, routineN,&
         "Joining neither rows nor columns is untested", __LINE__, error=error)
    !
    CALL dbcsr_data_clear (dense_data)
    IF (dbcsr_data_get_size(matrix%m%data_area) .GT. 0&
         .AND. nfullrows .GT. 0 .AND. nfullcols .GT. 0) THEN
!$OMP PARALLEL DEFAULT(none) &
!$OMP PRIVATE (block, iter, &
!$OMP         target_rs, target_cs, blk_row, blk_col, tr, blk_row_size, blk_col_size,&
!$OMP         blk_row_offset, blk_col_offset, tro, tco, target_offset,&
!$OMP         target_row_offset, target_col_offset,&
!$OMP         dense_row, dense_col) &
!$OMP SHARED (error,&
!$OMP         dense_data, matrix, &
!$OMP         make_tr, join_rows, join_cols, &
!$OMP         und_row_blk_offsets, und_col_blk_offsets,&
!$OMP         dense_row_blk_offsets, dense_col_blk_offsets,&
!$OMP         row_internal_offsets, col_internal_offsets,&
!$OMP         row_map, col_map,&
!$OMP         nfullrows, nfullcols)
       CALL dbcsr_data_init (block)
       CALL dbcsr_data_new (block,&
         dbcsr_type_1d_to_2d(dbcsr_data_get_type(dense_data)))
       CALL dbcsr_iterator_start (iter, matrix, dynamic=.TRUE., shared=.TRUE.,&
            contiguous_pointers=.FALSE., read_only=.TRUE.)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, blk_row, blk_col, block, tr,&
               row_size=blk_row_size, col_size=blk_col_size)
          dense_row = row_map (blk_row)
          dense_col = col_map (blk_col)
          !
          ! Calculate the target block row/column size and the offset
          ! within the target block where the undense block is placed.
          IF (join_rows) THEN
             target_row_offset = dense_row_blk_offsets(dense_row)
             target_rs = dense_row_blk_offsets(dense_row+1) -&
                         dense_row_blk_offsets(dense_row)
             tro = 1 + row_internal_offsets(blk_row)
          ELSE
             target_row_offset = und_row_blk_offsets(blk_row)
             target_rs = blk_row_size
             tro = 1
          ENDIF
          IF (join_cols) THEN
             target_col_offset = dense_col_blk_offsets(dense_col)
             target_cs =  dense_col_blk_offsets(dense_col+1) -&
                          dense_col_blk_offsets(dense_col)
             tco = 1 + col_internal_offsets(blk_col)
          ELSE
             target_col_offset = und_col_blk_offsets(blk_col)
             target_cs = und_col_blk_offsets(blk_col)
             tco = 1
          ENDIF
          target_offset = (target_row_offset-1) * nfullcols &
                        + (target_col_offset-1) * (&
                           dense_row_blk_offsets(dense_row+1) -&
                           dense_row_blk_offsets(dense_row))
          CALL dbcsr_block_partial_copy(dst=dense_data,&
               dst_offset=target_offset,&
               dst_rs=target_rs, dst_cs=target_cs, dst_tr=make_tr,&
               dst_r_lb = tro, dst_c_lb = tco,&
               src=block, src_rs=blk_row_size, src_cs=blk_col_size, src_tr=tr,&
               src_r_lb=1, src_c_lb=1, nrow=blk_row_size, ncol=blk_col_size)
       ENDDO
       CALL dbcsr_iterator_stop (iter)
       CALL dbcsr_data_clear_pointer (block)
       CALL dbcsr_data_release (block)
!$OMP END PARALLEL
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE make_dense_data

! *****************************************************************************
!> \brief Shuffles the data from standard dense to blocked form
!> \param[in,out] matrix             Matrix with data to fill
!> \param[in] dense_data             Dense data
!> \param[in] nfullrows              number of full rows in local submatrix
!> \param[in] nfullcols              number of full columns in local submatrix
!> \param[in] und_row_blk_offsets    row block offsets for undense data
!> \param[in] und_col_blk_offsets    column block offsets for undense data
!> \param[in] dense_row_blk_offsets  row block offsets for dense data
!> \param[in] dense_col_blk_offsets  column block offsets for dense data
!> \param[in] row_map                mapping from undense to dense rows
!> \param[in] col_map                mapping from undense to dense rows
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE make_undense_data (matrix, dense_data, nfullrows, nfullcols,&
       und_row_blk_offsets, und_col_blk_offsets,&
       dense_row_blk_offsets, dense_col_blk_offsets,&
       row_map, col_map, row_internal_offsets, col_internal_offsets, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_data_obj), INTENT(IN)         :: dense_data
    INTEGER, INTENT(IN)                      :: nfullrows, nfullcols
    INTEGER, DIMENSION(:), INTENT(IN) :: und_row_blk_offsets, &
      und_col_blk_offsets, dense_row_blk_offsets, dense_col_blk_offsets, &
      row_map, col_map, row_internal_offsets, col_internal_offsets
    TYPE(dbcsr_error_type)                   :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_undense_data', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk_col, blk_col_size, blk_row, blk_row_size, dense_col, &
      dense_col_offset, dense_cs, dense_offset, dense_row, dense_row_offset, &
      dense_rs, sco, sro
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: block
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (dbcsr_data_get_size (dense_data), 'GE', nfullrows*nfullcols,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Dense data too small",__LINE__,error)
    IF (dbcsr_data_get_size(matrix%m%data_area) .GT. 0) THEN
       CALL dbcsr_data_clear (matrix%m%data_area)
!$OMP PARALLEL DEFAULT(none) &
!$OMP PRIVATE (block, iter,&
!$OMP          blk_row, blk_col, tr,&
!$OMP          blk_row_size, blk_col_size, sro, sco,&
!$OMP          dense_row_offset, dense_col_offset, dense_row, dense_col,&
!$OMP          dense_cs, dense_rs,&
!$OMP          dense_offset) &
!$OMP SHARED (error,&
!$OMP         matrix, dense_data, &
!$OMP         nfullrows, nfullcols, und_row_blk_offsets, und_col_blk_offsets,&
!$OMP         dense_row_blk_offsets, dense_col_blk_offsets,&
!$OMP         row_map, col_map,&
!$OMP         row_internal_offsets, col_internal_offsets)
       CALL dbcsr_data_init (block)
       CALL dbcsr_data_new (block,&
            dbcsr_type_1d_to_2d(dbcsr_data_get_type(dense_data)))
       CALL dbcsr_iterator_start (iter, matrix, dynamic=.TRUE., shared=.TRUE.,&
            contiguous_pointers=.FALSE.)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, blk_row, blk_col, block, tr,&
               row_size=blk_row_size, col_size=blk_col_size)
          dense_row = row_map (blk_row)
          dense_col = col_map (blk_col)
          dense_row_offset = dense_row_blk_offsets(dense_row)
          dense_col_offset = dense_col_blk_offsets(dense_col)
          dense_rs = dense_row_blk_offsets(dense_row+1) -&
                     dense_row_blk_offsets(dense_row)
          dense_cs = dense_col_blk_offsets(dense_col+1) -&
                     dense_col_blk_offsets(dense_col)
          sro = 1 + row_internal_offsets(blk_row)
          sco = 1 + col_internal_offsets(blk_col)
          dense_offset = (dense_row_offset-1) * nfullcols &
                       + (dense_col_offset-1) * dense_rs
          CALL dbcsr_block_partial_copy (&
               dst=block, dst_rs=blk_row_size, dst_cs=blk_col_size, dst_tr=tr,&
               dst_r_lb=1, dst_c_lb=1,&
               src=dense_data, src_offset=dense_offset, &
               src_rs=dense_rs, src_cs=dense_cs, src_tr=.FALSE.,&
               src_r_lb=sro, src_c_lb=sco,&
               nrow=blk_row_size, ncol=blk_col_size)
       ENDDO
       CALL dbcsr_iterator_stop (iter)
       CALL dbcsr_data_clear_pointer (block)
       CALL dbcsr_data_release (block)
!$OMP END PARALLEL
    ENDIF
  END SUBROUTINE make_undense_data


! *****************************************************************************
!> \brief Replicates a DBCSR on all processors.
!> \param[in,out] matrix      matrix to replicate
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE dbcsr_replicate_all (matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CALL dbcsr_replicate (matrix, replicate_rows = .TRUE.,&
         replicate_columns=.TRUE., error=error)
  END SUBROUTINE dbcsr_replicate_all


! *****************************************************************************
!> \brief Replicates a DBCSR matrix among process rows and columns
!> \param[in,out] matrix      matrix to replicate
!> \param[in] replicate_rows  Row should be replicated among all
!>                            processors
!> \param[in] replicate_columns  Column should be replicated among
!>                               all processors
!> \param[in] restrict_source    (optional) Send only from this node
!>                               (ignores blocks on other nodes)
!> \param[in,out] error       error
!> \par Direction definition
!>      Row replication means that all processors in a process grid sharing
!>      the same row get the data of the entire row. (In a 1-column grid the
!>      operation has no effect.) Similar logic applies to column replication.
! *****************************************************************************
  SUBROUTINE dbcsr_replicate (matrix, replicate_rows, replicate_columns,&
       restrict_source, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    LOGICAL, INTENT(IN)                      :: replicate_rows, &
                                                replicate_columns
    INTEGER, INTENT(IN), OPTIONAL            :: restrict_source
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_replicate', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 3
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    CHARACTER                                :: rep_type
    INTEGER :: blk, blk_l, blk_p, blk_ps, blks, col, col_size, data_type, &
      dst_p, error_handler, mp_group, mynode, mypcol, myprow, nblks, &
      numnodes, nze, offset, row, row_size, smp, src_p, stored_blk_p, &
      stored_col, stored_row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rd_disp, recv_meta, rm_disp, &
                                                send_meta
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_count
    INTEGER, DIMENSION(2)                    :: send_count
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, col_dist, &
                                                row_blk_size, row_dist, &
                                                tmp_index
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    LOGICAL                                  :: had_subcomms, &
                                                i_am_restricted, rest_src, tr
    TYPE(dbcsr_data_obj)                     :: data_block, recv_data, &
                                                send_data, tmp_data
    TYPE(dbcsr_distribution_obj)             :: target_dist
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_obj)                          :: replicated

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_access_flush (matrix, error=error)
    CALL dbcsr_assert (dbcsr_valid_index(matrix%m),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Matrix not initialized.",__LINE__,error)
    !CALL dbcsr_assert (matrix%m%replication_type .EQ. dbcsr_repl_none,&
    !     dbcsr_warning_level, dbcsr_caller_error, routineN,&
    !     "Replicating a non-distributed matrix makes no sense.",__LINE__,error)
    IF (replicate_rows .AND. replicate_columns) THEN
       rep_type = dbcsr_repl_full
    ELSEIF (replicate_rows .AND. .NOT. replicate_columns) THEN
       rep_type = dbcsr_repl_row
    ELSEIF (replicate_columns .AND. .NOT. replicate_rows) THEN
       rep_type = dbcsr_repl_col
    ELSE
       rep_type = dbcsr_repl_none
       CALL dbcsr_assert(replicate_rows, "OR", replicate_columns, &
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Some replication must be specified", __LINE__, error=error)
    ENDIF
    data_type = dbcsr_get_data_type (matrix)
    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    target_dist = matrix%m%dist
    row_dist => array_data (dbcsr_distribution_row_dist (target_dist))
    col_dist => array_data (dbcsr_distribution_col_dist (target_dist))
    mp_obj = dbcsr_distribution_mp (target_dist)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    CALL dbcsr_assert(MAXVAL(row_dist).LE.UBOUND(blacs2mpi,1), dbcsr_failure_level,&
         dbcsr_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows',__LINE__,error)
    IF (dbg) THEN
       CALL dbcsr_assert(MAXVAL(row_dist).EQ.UBOUND(blacs2mpi,1),&
            dbcsr_warning_level, dbcsr_wrong_args_error, routineN,&
            'Range of row distribution not equal to processor rows',__LINE__,error)
    ENDIF
    CALL dbcsr_assert(MAXVAL(col_dist).LE.UBOUND(blacs2mpi,2), dbcsr_failure_level,&
         dbcsr_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols',__LINE__,error)
    IF (dbg) THEN
       CALL dbcsr_assert(MAXVAL(col_dist).EQ.UBOUND(blacs2mpi,2),&
            dbcsr_warning_level, dbcsr_wrong_args_error, routineN,&
            'Range of col distribution not equal to processor cols',__LINE__,error)
    ENDIF
    ! Define the number of nodes with which I will communicate. Also
    ! setup row and column communicators.
    had_subcomms = dbcsr_mp_has_subgroups (mp_obj)
    SELECT CASE (rep_type)
    CASE (dbcsr_repl_full)
       numnodes = dbcsr_mp_numnodes (mp_obj)
       mp_group = dbcsr_mp_group (mp_obj)
       mynode = dbcsr_mp_mynode (mp_obj)
    CASE (dbcsr_repl_row)
       numnodes = dbcsr_mp_npcols (mp_obj)
       CALL dbcsr_mp_grid_setup (mp_obj, force=.TRUE.)
       mp_group = dbcsr_mp_my_row_group (mp_obj)
       mynode = dbcsr_mp_mypcol (mp_obj)
    CASE (dbcsr_repl_col)
       numnodes = dbcsr_mp_nprows (mp_obj)
       CALL dbcsr_mp_grid_setup (mp_obj, force=.TRUE.)
       mp_group = dbcsr_mp_my_col_group (mp_obj)
       mynode = dbcsr_mp_myprow (mp_obj)
    CASE (dbcsr_repl_none)
       numnodes = 1
       mp_group = dbcsr_mp_group (mp_obj)
       mynode = 0
    END SELECT
    CALL dbcsr_assert (rep_type.EQ.dbcsr_repl_row &
         .OR. rep_type.EQ.dbcsr_repl_col,&
         "IMP", dbcsr_mp_has_subgroups (mp_obj), dbcsr_fatal_level,&
         dbcsr_unimplemented_error_nr, routineN,&
         "Only full replication supported when subcommunicators are turned off.",&
         __LINE__, error=error)
    !
    IF (PRESENT (restrict_source)) THEN
       rest_src = .TRUE.
       i_am_restricted = mynode .NE. restrict_source
    ELSE
       rest_src = .FALSE.
       i_am_restricted = .FALSE.
    ENDIF
    !
    ALLOCATE (recv_count(2, 0:numnodes-1))
    ALLOCATE (rd_disp(0:numnodes-1))
    ALLOCATE (rm_disp(0:numnodes-1))
    CALL dbcsr_init (replicated)
    CALL dbcsr_create(replicated, name='Replicated '//TRIM(matrix%m%name),&
         template = matrix,&
         matrix_type = dbcsr_type_no_symmetry,&
         replication_type=rep_type, error=error)
    !
    ! Count initial sizes for sending. Also, ensure that blocks are on their
    ! home processors.
    send_count(1:2) = 0
    CALL dbcsr_iterator_start (iter, matrix)
    IF (.NOT. i_am_restricted) THEN
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, row, col, &
               row_size=row_size, col_size=col_size, blk=blk)
          !tr = .FALSE.
          !CALL dbcsr_get_stored_coordinates (matrix, row, col, tr, dst_p)
          !CALL dbcsr_assert (dst_p, "EQ", mynode, dbcsr_fatal_level,&
          !     dbcsr_wrong_args_error, routineN, &
          !     "Matrix is not correctly distributed. Call dbcsr_redistribute.",&
          !     __LINE__, error=error)
          nze = row_size * col_size
          send_count(1) = send_count(1) + 1
          send_count(2) = send_count(2) + nze
       ENDDO
       send_count(2) = dbcsr_get_data_size_referenced (matrix)
    ENDIF
    CALL dbcsr_iterator_stop (iter)
    ! Exchange how much data others have.
    CALL mp_allgather (send_count(1:2), recv_count, mp_group)
    CALL dbcsr_data_init (recv_data)
    nze = SUM (recv_count(2, :))
    nblks = SUM (recv_count(1, :))
    CALL dbcsr_data_new (recv_data, data_type=data_type, data_size=nze)
    ! send_data should have the correct size
    CALL dbcsr_data_init (send_data)
    IF (send_count(2) .EQ. 0) THEN
       CALL dbcsr_data_new (send_data, data_type=data_type, data_size=0)
    ELSE
       CALL dbcsr_data_new (send_data, data_type=data_type)
       send_data = pointer_view (send_data, matrix%m%data_area, 1, send_count(2))
    ENDIF
    ALLOCATE (recv_meta(metalen * nblks))
    ALLOCATE (send_meta(metalen * send_count(1)))
    recv_meta(:) = 0
    ! Fill in the meta data structures and copy the data.
    rd_disp = -1 ; rm_disp = -1
    rd_disp(0) = 1 ; rm_disp(0) = 1
    DO dst_p = 1, numnodes-1
       rm_disp(dst_p) = rm_disp(dst_p-1)&
                        + metalen*recv_count(1, dst_p-1)
       rd_disp(dst_p) = rd_disp(dst_p-1)&
                        + recv_count(2, dst_p-1)
    ENDDO
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_new (data_block, data_type=data_type)
    CALL dbcsr_iterator_start (iter, matrix)
    smp = 1
    IF (.NOT. i_am_restricted) THEN
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, row, col, blk,&
               transposed=tr, blk_p=blk_p)
          send_meta(smp+0) = row
          send_meta(smp+1) = col
          send_meta(smp+2) = blk_p
          smp = smp + metalen
       ENDDO
    ENDIF
    CALL dbcsr_iterator_stop (iter)
    CALL dbcsr_data_clear_pointer (data_block)
    CALL dbcsr_data_release (data_block)
    ! Exchange the data and metadata structures.
    CALL mp_allgather (send_meta, recv_meta, metalen*recv_count(1,:),&
         rm_disp-1, mp_group)
    CALL dbcsr_allgatherv (send_data, recv_data, recv_count(2,:),&
         rd_disp-1, mp_group)
    ! Release the send buffer. If it had a non-zero size then it was a
    ! pointer into the regular matrix and the data pointer should be
    ! cleared and not deallocated.
    IF (send_count(2) .NE. 0) THEN
       CALL dbcsr_data_clear_pointer (send_data)
    ENDIF
    CALL dbcsr_data_release (send_data)
    !
    ! Now fill in the data.
    CALL dbcsr_work_create(replicated,&
         SUM(recv_count(1,:)),&
         SUM(recv_count(2,:)), n=1, &
         work_mutable=.FALSE.,&
         error=error)
    CALL dbcsr_data_hold (recv_data)
    CALL dbcsr_data_release (replicated%m%wms(1)%data_area)
    replicated%m%wms(1)%data_area = recv_data
    blk_ps = 1
    blks = 1
    DO src_p = 0, numnodes-1
       nze = recv_count(2, src_p)
       !CALL dbcsr_data_set (replicated%m%wms(1)%data_area, blk_ps, nze,&
       !     recv_data, rd_disp(src_p))
       offset = rd_disp(src_p) - 1
       DO blk_l = 1, recv_count(1, src_p)
          IF (dbg) WRITE(*,*)"src_p, blk_l", src_p, blk_l
          stored_row = recv_meta(rm_disp(src_p)+metalen*(blk_l-1))
          stored_col = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+1)
          stored_blk_p = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+2)
          replicated%m%wms(1)%row_i(blks) = stored_row
          replicated%m%wms(1)%col_i(blks) = stored_col
          replicated%m%wms(1)%blk_p(blks) = SIGN(ABS(stored_blk_p) + offset,&
                                                 stored_blk_p)
          nze = row_blk_size(stored_row)&
               * col_blk_size(stored_col)
          blk_ps = MAX (blk_ps, ABS(stored_blk_p) + nze + offset)
          blks = blks + 1
       ENDDO
    ENDDO
    CALL dbcsr_data_set_size_referenced (replicated%m%wms(1)%data_area, blk_ps - 1)
    !
    replicated%m%wms(1)%lastblk = blks - 1
    replicated%m%wms(1)%datasize = blk_ps - 1
    CALL dbcsr_finalize(replicated, reshuffle=.TRUE., error=error)
    !
    ! Remove communicators if they were forcibly created.
    IF (had_subcomms .AND.&
         (rep_type .EQ. dbcsr_repl_row .OR. rep_type .EQ. dbcsr_repl_col)) THEN
       CALL dbcsr_mp_grid_remove (mp_obj)
    ENDIF
    DEALLOCATE(recv_count)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    CALL dbcsr_data_release (recv_data)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_meta)
    matrix%m%replication_type = replicated%m%replication_type
    ! Now replace the data and index
    CALL dbcsr_switch_data_area (matrix, replicated%m%data_area,&
         previous_data_area=tmp_data, error=error)
    CALL dbcsr_switch_data_area (replicated, tmp_data, error=error)
    CALL dbcsr_data_release (tmp_data)
    tmp_index => matrix%m%index
    matrix%m%index => replicated%m%index
    replicated%m%index => tmp_index
    CALL dbcsr_repoint_index (matrix%m)
    matrix%m%nze = replicated%m%nze
    matrix%m%nblks = replicated%m%nblks
    CALL dbcsr_release (replicated)
    CALL dbcsr_verify_matrix (matrix, error=error)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_replicate


! *****************************************************************************
!> \brief Fully redistributes a DBCSR matrix.
!>
!>        The new distribution may be arbitrary as long as the total
!>        number full rows and columns matches that of the existing
!>        matrix.
!> \param[in] matrix          matrix to redistribute
!> \param[in,out] redist      redistributed matrix
!> \param[in] keep_sparsity   (optional) retains the sparsity of the redist
!>                            matrix
! *****************************************************************************
  SUBROUTINE dbcsr_complete_redistribute(matrix, redist, keep_sparsity, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_obj), INTENT(INOUT)           :: redist
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_complete_redistribute', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 7
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER :: blk, blk_col_new, blk_ps, blk_row_new, blks, cnt_fnd, cnt_new, &
      cnt_skip, col, col_int, col_offset_new, col_offset_old, col_rle, &
      col_size, col_size_new, data_offset_l, data_type, dst_p, error_handler, &
      i, meta_l, mp_group, numnodes, nze_rle, row, row_int, row_offset_new, &
      row_offset_old, row_rle, row_size, row_size_new, src_p, stored_col_new, &
      stored_row_new
    INTEGER, ALLOCATABLE, DIMENSION(:) :: col_end_new, col_end_old, &
      col_start_new, col_start_old, rd_disp, recv_meta, rm_disp, row_end_new, &
      row_end_old, row_start_new, row_start_old, sd_disp, sdp, send_meta, &
      sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: col_reblocks, n_col_reblocks, &
      n_row_reblocks, recv_count, row_reblocks, send_count, total_recv_count, &
      total_send_count
    INTEGER, DIMENSION(:), POINTER :: col_blk_size_new, col_blk_size_old, &
      col_dist_new, row_blk_size_new, row_blk_size_old, row_dist_new
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    LOGICAL                                  :: found, my_keep_sparsity, sym, &
                                                tr, valid_block
    REAL(kind=dp)                            :: cs1, cs2
    TYPE(dbcsr_data_obj)                     :: buff_data, data_block, &
                                                recv_data, send_data
    TYPE(dbcsr_distribution_obj)             :: dist_new
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj_new

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_valid_index(matrix%m),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Input not valid.",__LINE__,error)
    CALL dbcsr_assert (dbcsr_is_initialized(redist%m),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Output matrix not initialized.",__LINE__,error)
    CALL dbcsr_assert (matrix%m%replication_type .EQ. dbcsr_repl_none,&
         dbcsr_warning_level, dbcsr_caller_error, routineN,&
         "Can not redistribute replicated matrix.",__LINE__,error)
    CALL dbcsr_access_flush (matrix, error=error)
    !
    my_keep_sparsity = .FALSE.
    IF (PRESENT (keep_sparsity)) my_keep_sparsity = keep_sparsity
    !
    sym = dbcsr_has_symmetry (redist)
    data_type = matrix%m%data_type
    ! Get row and column start and end positions
    ! Old matrix
    row_blk_size_old => array_data (matrix%m%row_blk_size)
    col_blk_size_old => array_data (matrix%m%col_blk_size)
    ALLOCATE (row_start_old (dbcsr_nblkrows_total (matrix)),&
         row_end_old (dbcsr_nblkrows_total (matrix)),&
         col_start_old (dbcsr_nblkcols_total (matrix)),&
         col_end_old (dbcsr_nblkcols_total (matrix)))
    CALL convert_sizes_to_offsets (row_blk_size_old,&
         row_start_old, row_end_old)
    CALL convert_sizes_to_offsets (col_blk_size_old,&
         col_start_old, col_end_old)
    ! New matrix
    dist_new = dbcsr_distribution (redist)
    row_blk_size_new => array_data (redist%m%row_blk_size)
    col_blk_size_new => array_data (redist%m%col_blk_size)
    ALLOCATE (row_start_new (dbcsr_nblkrows_total (redist)),&
         row_end_new (dbcsr_nblkrows_total (redist)),&
         col_start_new (dbcsr_nblkcols_total (redist)),&
         col_end_new (dbcsr_nblkcols_total (redist)))
    CALL convert_sizes_to_offsets (row_blk_size_new,&
         row_start_new, row_end_new)
    CALL convert_sizes_to_offsets (col_blk_size_new,&
         col_start_new, col_end_new)
    row_dist_new => array_data (dbcsr_distribution_row_dist (dist_new))
    col_dist_new => array_data (dbcsr_distribution_col_dist (dist_new))
    ! Create mappings
    i = dbcsr_nfullrows_total (redist)
    ALLOCATE (row_reblocks (4,i))
    ALLOCATE (n_row_reblocks (2, dbcsr_nblkrows_total (matrix)))
    CALL dbcsr_reblocking_targets (row_reblocks, i, n_row_reblocks,&
         row_blk_size_old, row_blk_size_new, error=error)
    i = dbcsr_nfullcols_total (redist)
    ALLOCATE (col_reblocks (4,i))
    ALLOCATE (n_col_reblocks (2, dbcsr_nblkcols_total (matrix)))
    CALL dbcsr_reblocking_targets (col_reblocks, i, n_col_reblocks,&
         col_blk_size_old, col_blk_size_new, error=error)
    !
    mp_obj_new = dbcsr_distribution_mp (dist_new)
    pgrid => dbcsr_mp_pgrid (mp_obj_new)
    numnodes = dbcsr_mp_numnodes (mp_obj_new)
    mp_group = dbcsr_mp_group (mp_obj_new)
    !
    CALL dbcsr_assert(MAXVAL(row_dist_new),"LE",UBOUND(pgrid,1), dbcsr_failure_level,&
         dbcsr_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows',__LINE__,error)
    IF (dbg) THEN
       CALL dbcsr_assert(MAXVAL(row_dist_new),"EQ",UBOUND(pgrid,1),&
            dbcsr_warning_level, dbcsr_wrong_args_error, routineN,&
            'Range of row distribution not equal to processor rows',__LINE__,error)
    ENDIF
    CALL dbcsr_assert(MAXVAL(col_dist_new),"LE",UBOUND(pgrid,2), dbcsr_failure_level,&
         dbcsr_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols',__LINE__,error)
    IF (dbg) THEN
       CALL dbcsr_assert(MAXVAL(col_dist_new),"EQ",UBOUND(pgrid,2),&
            dbcsr_warning_level, dbcsr_wrong_args_error, routineN,&
            'Range of col distribution not equal to processor cols',__LINE__,error)
    ENDIF
    ALLOCATE (send_count(2, 0:numnodes-1))
    ALLOCATE (recv_count(2, 0:numnodes-1))
    ALLOCATE (total_send_count(2, 0:numnodes-1))
    ALLOCATE (total_recv_count(2, 0:numnodes-1))
    ALLOCATE (sdp(0:numnodes-1))
    ALLOCATE (sd_disp(0:numnodes-1))
    ALLOCATE (smp(0:numnodes-1))
    ALLOCATE (sm_disp(0:numnodes-1))
    ALLOCATE (rd_disp(0:numnodes-1))
    ALLOCATE (rm_disp(0:numnodes-1))
    IF (dbg) THEN
       cs1 = dbcsr_checksum (matrix, error=error)
    ENDIF
    !cs1 = dbcsr_checksum (matrix)
    !call dbcsr_print(matrix)
    !
    !
    ! Count initial sizes for sending.
    !
    ! We go through every element of every local block and determine
    ! to which processor it must be sent. It could be more efficient,
    ! but at least the index data are run-length encoded.
    send_count(:,:) = 0
    CALL dbcsr_iterator_start (iter, matrix)
    dst_p = -1
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk)
       DO col_int = n_col_reblocks (1,col),&
            n_col_reblocks(1,col) + n_col_reblocks(2,col)-1
          blk_col_new = col_reblocks(1, col_int)
          DO row_int = n_row_reblocks (1, row),&
               n_row_reblocks(1,row) + n_row_reblocks(2,row)-1
             blk_row_new = row_reblocks(1, row_int)
             IF (.NOT. sym .OR. blk_col_new .GE. blk_row_new) THEN
                tr = .FALSE.
                CALL dbcsr_get_stored_coordinates (redist%m,&
                     blk_row_new, blk_col_new, tr, dst_p)
                send_count(1, dst_p) = send_count(1, dst_p) + 1
                send_count(2, dst_p) = send_count(2, dst_p) +&
                     col_reblocks(2,col_int)*row_reblocks(2,row_int)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    !
    !
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    ! Allocate data structures needed for data exchange.
    CALL dbcsr_data_init (recv_data)
    CALL dbcsr_data_new (recv_data, data_type, SUM(recv_count(2, :)))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :))))
    CALL dbcsr_data_init (send_data)
    CALL dbcsr_data_new (send_data, data_type, SUM(send_count(2, :)))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :))))
    ! Fill in the meta data structures and copy the data.
    DO dst_p = 0, numnodes-1
       total_send_count(1, dst_p) = send_count (1, dst_p)
       total_send_count(2, dst_p) = send_count (2, dst_p)
       total_recv_count(1, dst_p) = recv_count (1, dst_p)
       total_recv_count(2, dst_p) = recv_count (2, dst_p)
    ENDDO
    sd_disp = -1 ; sm_disp = -1
    rd_disp = -1 ; rm_disp = -1
    sd_disp(0) = 1 ; sm_disp(0) = 1
    rd_disp(0) = 1 ; rm_disp(0) = 1
    DO dst_p = 1, numnodes-1
       sm_disp(dst_p) = sm_disp(dst_p-1)&
                        + metalen*total_send_count(1, dst_p-1)
       sd_disp(dst_p) = sd_disp(dst_p-1)&
                        + total_send_count(2, dst_p-1)
       rm_disp(dst_p) = rm_disp(dst_p-1)&
                        + metalen*total_recv_count(1, dst_p-1)
       rd_disp(dst_p) = rd_disp(dst_p-1)&
                        + total_recv_count(2, dst_p-1)
    ENDDO
    sdp(:) = sd_disp     ! sdp points to the the next place to store
                         ! data. It is postincremented.
    smp(:) = sm_disp - metalen  ! But smp points to the "working" data, not
                                ! the next. It is pre-incremented, so we must
                                ! first rewind it.
    !
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_new (data_block, data_type)
    CALL dbcsr_iterator_start (iter, matrix)
    dst_p = -1
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, data_block, tr, blk,&
            row_size=row_size, col_size=col_size)
       !IF (tr) WRITE(*,*)"block at",row,col," is transposed"
       DO col_int = n_col_reblocks (1,col),&
            n_col_reblocks(1,col) + n_col_reblocks(2,col)-1
          blk_col_new = col_reblocks(1, col_int)
          DO row_int = n_row_reblocks (1, row),&
               n_row_reblocks(1,row) + n_row_reblocks(2,row)-1
             blk_row_new = row_reblocks(1, row_int)
             loc_ok: IF (.NOT. sym .OR. blk_col_new .GE. blk_row_new) THEN
                IF (dbg) &
                     WRITE(*,*)'using block',blk_row_new,'x',blk_col_new
                ! Start a new RLE run
                tr = .FALSE.
                CALL dbcsr_get_stored_coordinates (redist%m,&
                     blk_row_new, blk_col_new, tr, dst_p)
                row_offset_old = row_reblocks(3, row_int)
                col_offset_old = col_reblocks(3, col_int)
                row_offset_new = row_reblocks(4, row_int)
                col_offset_new = col_reblocks(4, col_int)
                row_rle = row_reblocks(2,row_int)
                col_rle = col_reblocks(2,col_int)
                smp(dst_p) = smp(dst_p) + metalen
                send_meta(smp(dst_p)) = blk_row_new   ! new blocked row
                send_meta(smp(dst_p)+1) = blk_col_new ! new blocked column
                send_meta(smp(dst_p)+2) = row_offset_new  ! row in new block
                send_meta(smp(dst_p)+3) = col_offset_new  ! col in new block
                send_meta(smp(dst_p)+4) = row_rle ! RLE rows
                send_meta(smp(dst_p)+5) = col_rle ! RLE columns
                send_meta(smp(dst_p)+6) = sdp(dst_p)-sd_disp(dst_p) ! Offset in data
                nze_rle = row_rle * col_rle
                ! Copy current block into the send buffer
                CALL dbcsr_block_partial_copy(&
                     send_data, dst_offset=sdp(dst_p)-1,&
                     dst_rs=row_rle, dst_cs=col_rle, dst_tr=.FALSE.,&
                     dst_r_lb=1, dst_c_lb=1,&
                     src=data_block,&
                     src_rs=row_size, src_cs=col_size, src_tr=tr,&
                     src_r_lb=row_offset_old, src_c_lb=col_offset_old,&
                     nrow=row_rle, ncol=col_rle)
                sdp(dst_p) = sdp(dst_p)+nze_rle
             ENDIF loc_ok
          ENDDO ! row_int
       ENDDO ! col_int
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    CALL dbcsr_data_clear_pointer (data_block)
    CALL dbcsr_data_release (data_block)

    ! Exchange the data and metadata structures.
    !
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       CALL hybrid_alltoall_s1(&
            send_data%d%r_sp(:), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%r_sp(:), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_real_8)
       !CALL mp_alltoall(&
       !     send_data%d%r_dp(:), total_send_count(2,:), sd_disp(:)-1,&
       !     recv_data%d%r_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
       !     mp_group)
       CALL hybrid_alltoall_d1 (&
            send_data%d%r_dp(:), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%r_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_complex_4)
       CALL hybrid_alltoall_c1(&
            send_data%d%c_sp(:), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%c_sp(:), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_complex_8)
       CALL hybrid_alltoall_z1(&
            send_data%d%c_dp(:), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%c_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid matrix type",__LINE__,error)
    END SELECT
    CALL hybrid_alltoall_i1(send_meta(:), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), metalen*total_recv_count(1,:), rm_disp(:)-1, mp_obj_new)
    !
    ! Now fill in the data.
    CALL dbcsr_work_create(redist,&
         nblks_guess=SUM(recv_count(1,:)),&
         sizedata_guess=SUM(recv_count(2,:)), work_mutable=.TRUE.,error=error)
    CALL dbcsr_data_init (buff_data)
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_new (buff_data, dbcsr_type_1d_to_2d (data_type), &
         redist%m%max_rbs, redist%m%max_cbs)
    CALL dbcsr_data_new (data_block, dbcsr_type_1d_to_2d (data_type))

    CALL dbcsr_access_start (redist, error=error)
    !blk_p = 1
    !blk = 1
    blk_ps = 0
    blks = 0
    cnt_fnd = 0 ; cnt_new = 0 ; cnt_skip = 0
    DO src_p = 0, numnodes-1
       data_offset_l = rd_disp(src_p)
       DO meta_l = 1, recv_count(1, src_p)
          stored_row_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1))
          stored_col_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+1)
          row_offset_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+2)
          col_offset_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+3)
          row_rle = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+4)
          col_rle = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+5)
          data_offset_l = rd_disp(src_p)&
               + recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+6)

          CALL dbcsr_data_clear_pointer (data_block)
          CALL dbcsr_get_block_p(redist, stored_row_new, stored_col_new,&
               data_block, tr, found)
          valid_block = found

          IF (found) cnt_fnd = cnt_fnd + 1
          IF (.NOT. found .AND. .NOT. my_keep_sparsity) THEN
             ! We have to set up a buffer block
             CALL dbcsr_data_set_pointer (data_block,&
                  rsize=row_blk_size_new (stored_row_new),&
                  csize=col_blk_size_new (stored_col_new),&
                  pointee=buff_data)
             CALL dbcsr_data_clear (data_block)
             !r2_dp => r2_dp_buff(1:row_blk_size_new (stored_row_new),&
             !     1:col_blk_size_new (stored_col_new))
             !r2_dp(:,:) = 0.0_dp
             tr = .FALSE.
             blks = blks + 1
             blk_ps = blk_ps + row_blk_size_new (stored_row_new) * &
                  col_blk_size_new (stored_col_new)
             valid_block = .TRUE.
             cnt_new = cnt_new + 1
          ENDIF
          nze_rle = row_rle*col_rle

          IF (valid_block) THEN
             row_size_new = row_blk_size_new (stored_row_new)
             col_size_new = col_blk_size_new (stored_col_new)
             CALL dbcsr_block_partial_copy (&
                  dst=data_block, dst_tr=tr,&
                  dst_rs=row_size_new, dst_cs=col_size_new, &
                  dst_r_lb=row_offset_new, dst_c_lb=col_offset_new,&
                  src=recv_data, src_offset=data_offset_l-1,&
                  src_rs=row_rle, src_cs=col_rle, src_tr=.FALSE.,&
                  src_r_lb=1, src_c_lb=1,&
                  nrow=row_rle, ncol=col_rle)
          ELSE
             cnt_skip = cnt_skip+1
          ENDIF

          data_offset_l = data_offset_l + nze_rle
          IF (.NOT. found .AND. valid_block) THEN
             IF (dbg) WRITE(*,*)routineN//" Adding new block at",&
                  stored_row_new, stored_col_new
             CALL dbcsr_put_block(redist, stored_row_new, stored_col_new,&
                  data_block, tr)
             !DEALLOCATE (r2_dp)
          ELSE
             IF (.NOT. my_keep_sparsity .AND. dbg) &
                  WRITE(*,*)routineN//" Reusing block at",&
                  stored_row_new, stored_col_new
          ENDIF
       ENDDO
    ENDDO

    CALL dbcsr_data_clear_pointer(data_block)
    CALL dbcsr_data_release (buff_data)
    CALL dbcsr_data_release (data_block)
    CALL dbcsr_access_stop (redist, error=error)
    !
    IF (dbg) THEN
       WRITE(*,*)routineN//" Declared blocks=",redist%m%wms(1)%lastblk,&
            "actual=",blks
       WRITE(*,*)routineN//" Declared data size=",redist%m%wms(1)%datasize,&
            "actual=",blk_ps
    ENDIF

    CALL dbcsr_finalize(redist, error=error)

    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)

    CALL dbcsr_data_release (recv_data)
    CALL dbcsr_data_release (send_data)

    DEALLOCATE(recv_meta)
    DEALLOCATE(send_meta)
   
    !if (dbg) call dbcsr_print(redist)
    IF (dbg) THEN
       cs2 = dbcsr_checksum (redist, error=error)
       WRITE(*,*)routineN//" Checksums=",cs1, cs2, cs1-cs2
    ENDIF
    !CALL dbcsr_assert (cs1-cs2 .LT. 0.00001, dbcsr_fatal_level, dbcsr_internal_error,&
    !     routineN, "Mangled data!")
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_complete_redistribute

! *****************************************************************************
!> \brief Redistributes a DBCSR matrix.
!>
!>        The new distribution should have compatible row and column blocks.
!> \param[in] matrix          matrix to redistribute
!> \param[in,out] redist      redistributed matrix, which should already be
!>                            created
! *****************************************************************************
  SUBROUTINE dbcsr_redistribute(matrix, redist, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_obj), INTENT(INOUT)           :: redist
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_redistribute', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 2
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER :: blk, blk_ps, blks, col, col_size, data_type, dst_p, &
      error_handler, meta_l, mp_group, numnodes, nze, row, row_size, src_p, &
      stored_col_new, stored_row_new
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rd_disp, recv_meta, rm_disp, &
                                                sd_disp, sdp, send_meta, &
                                                sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_count, send_count, &
                                                total_recv_count, &
                                                total_send_count
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size_new, &
                                                col_dist_new, &
                                                row_blk_size_new, row_dist_new
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    LOGICAL                                  :: sym_tr, tr
    TYPE(dbcsr_data_obj)                     :: data_block, recv_data, &
                                                send_data
    TYPE(dbcsr_distribution_obj)             :: dist_new
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj_new

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_access_flush (matrix, error=error)
    CALL dbcsr_access_flush (redist, error=error)
    !call dbcsr_print_dist (matrix%m%dist)
    !call dbcsr_print_dist (redist%m%dist)
    CALL dbcsr_assert (dbcsr_valid_index(matrix%m),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Input not valid.",__LINE__,error)
    CALL dbcsr_assert (dbcsr_is_initialized(redist%m),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Output matrix not initialized.",__LINE__,error)
    CALL dbcsr_assert (matrix%m%replication_type .EQ. dbcsr_repl_none,&
         dbcsr_warning_level, dbcsr_caller_error, routineN,&
         "Can not redistribute replicated matrix.",__LINE__,error)
    data_type = matrix%m%data_type
    ! Get row and column start and end positions
    ! Old matrix
    ! New matrix
    dist_new = dbcsr_distribution (redist)
    row_blk_size_new => array_data (redist%m%row_blk_size)
    col_blk_size_new => array_data (redist%m%col_blk_size)
    row_dist_new => array_data (dbcsr_distribution_row_dist (dist_new))
    col_dist_new => array_data (dbcsr_distribution_col_dist (dist_new))
    !
    mp_obj_new = dbcsr_distribution_mp (dist_new)
    pgrid => dbcsr_mp_pgrid (mp_obj_new)
    numnodes = dbcsr_mp_numnodes (mp_obj_new)
    mp_group = dbcsr_mp_group (mp_obj_new)
    !
    CALL dbcsr_assert(MAXVAL(row_dist_new).LE.UBOUND(pgrid,1), dbcsr_failure_level,&
         dbcsr_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows',__LINE__,error)
    IF (dbg) THEN
       CALL dbcsr_assert(MAXVAL(row_dist_new).EQ.UBOUND(pgrid,1),&
            dbcsr_warning_level, dbcsr_wrong_args_error, routineN,&
            'Range of row distribution not equal to processor rows',__LINE__,error)
    ENDIF
    CALL dbcsr_assert(MAXVAL(col_dist_new).LE.UBOUND(pgrid,2), dbcsr_failure_level,&
         dbcsr_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols',__LINE__,error)
    IF (dbg) THEN
       CALL dbcsr_assert(MAXVAL(col_dist_new).EQ.UBOUND(pgrid,2),&
            dbcsr_warning_level, dbcsr_wrong_args_error, routineN,&
            'Range of col distribution not equal to processor cols',__LINE__,error)
    ENDIF
    ALLOCATE (send_count(2, 0:numnodes-1))
    ALLOCATE (recv_count(2, 0:numnodes-1))
    ALLOCATE (total_send_count(2, 0:numnodes-1))
    ALLOCATE (total_recv_count(2, 0:numnodes-1))
    ALLOCATE (sdp(0:numnodes-1))
    ALLOCATE (sd_disp(0:numnodes-1))
    ALLOCATE (smp(0:numnodes-1))
    ALLOCATE (sm_disp(0:numnodes-1))
    ALLOCATE (rd_disp(0:numnodes-1))
    ALLOCATE (rm_disp(0:numnodes-1))
    ! Count initial sizes for sending.
    !
    send_count(:,:) = 0
    CALL dbcsr_iterator_start (iter, matrix)
    dst_p = -1
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk, tr,&
            row_size=row_size, col_size=col_size)
       sym_tr = .FALSE.
       CALL dbcsr_get_stored_coordinates (redist%m,&
            row, col, sym_tr, dst_p)
       nze = row_size*col_size
       send_count(1, dst_p) = send_count(1, dst_p) + 1
       send_count(2, dst_p) = send_count(2, dst_p) + nze
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    ! Allocate data structures needed for data exchange.
    CALL dbcsr_data_init (recv_data)
    CALL dbcsr_data_new (recv_data, data_type, SUM(recv_count(2, :)))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :))))
    CALL dbcsr_data_init (send_data)
    CALL dbcsr_data_new (send_data, data_type, SUM(send_count(2, :)))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :))))
    ! Fill in the meta data structures and copy the data.
    DO dst_p = 0, numnodes-1
       total_send_count(1, dst_p) = send_count (1, dst_p)
       total_send_count(2, dst_p) = send_count (2, dst_p)
       total_recv_count(1, dst_p) = recv_count (1, dst_p)
       total_recv_count(2, dst_p) = recv_count (2, dst_p)
    ENDDO
    sd_disp = -1 ; sm_disp = -1
    rd_disp = -1 ; rm_disp = -1
    sd_disp(0) = 1 ; sm_disp(0) = 1
    rd_disp(0) = 1 ; rm_disp(0) = 1
    DO dst_p = 1, numnodes-1
       sm_disp(dst_p) = sm_disp(dst_p-1)&
                        + metalen*total_send_count(1, dst_p-1)
       sd_disp(dst_p) = sd_disp(dst_p-1)&
                        + total_send_count(2, dst_p-1)
       rm_disp(dst_p) = rm_disp(dst_p-1)&
                        + metalen*total_recv_count(1, dst_p-1)
       rd_disp(dst_p) = rd_disp(dst_p-1)&
                        + total_recv_count(2, dst_p-1)
    ENDDO
    sdp(:) = sd_disp     ! sdp points to the the next place to store
                         ! data. It is postincremented.
    smp(:) = sm_disp - metalen  ! But smp points to the "working" data, not
                                ! the next. It is pre-incremented, so we must
                                ! first rewind it.
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_new (data_block, data_type)
    CALL dbcsr_iterator_start (iter, matrix)
    dst_p = -1
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, data_block, tr, blk)
       !IF (tr) WRITE(*,*)"block at",row,col," is transposed"
       sym_tr = .FALSE.
       CALL dbcsr_get_stored_coordinates (redist%m,&
            row, col, sym_tr, dst_p)
       smp(dst_p) = smp(dst_p) + metalen
       IF (tr) THEN
          send_meta(smp(dst_p)) = -row
       ELSE
          send_meta(smp(dst_p)) = row
       ENDIF
       send_meta(smp(dst_p)+1) = col ! new blocked column
       nze = dbcsr_data_get_size (data_block)
       CALL dbcsr_data_set (send_data, lb=sdp(dst_p), data_size=nze,&
            src=data_block, source_lb=1)
       !send_data(sdp(dst_p):sdp(dst_p)+SIZE(r_dp)-1) &
       !     = r_dp(:)
       sdp(dst_p) = sdp(dst_p) + nze
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    CALL dbcsr_data_clear_pointer (data_block)
    ! Exchange the data and metadata structures.
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       CALL hybrid_alltoall_s1(&
            send_data%d%r_sp(:), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%r_sp(:), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_real_8)
       !CALL mp_alltoall(&
       !     send_data%d%r_dp(:), total_send_count(2,:), sd_disp(:)-1,&
       !     recv_data%d%r_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
       !     mp_group)
       CALL hybrid_alltoall_d1 (&
            send_data%d%r_dp(:), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%r_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_complex_4)
       CALL hybrid_alltoall_c1(&
            send_data%d%c_sp(:), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%c_sp(:), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_complex_8)
       CALL hybrid_alltoall_z1(&
            send_data%d%c_dp(:), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%c_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    END SELECT
    !CALL mp_alltoall(send_data(:), total_send_count(2,:), sd_disp(:)-1,&
    !     recv_data(:), total_recv_count(2,:), rd_disp(:)-1, mp_group)
    CALL hybrid_alltoall_i1(send_meta(:), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), metalen*total_recv_count(1,:), rm_disp(:)-1, mp_obj_new)
    ! Now fill in the data.
    CALL dbcsr_work_create(redist,&
            SUM(recv_count(1,:)),&
            SUM(recv_count(2,:)), work_mutable=.FALSE., n=1, error=error)
    !
    blk_ps = 1
    blks = 0
    DO src_p = 0, numnodes-1
       !data_offset_l = rd_disp(src_p)
       DO meta_l = 1, recv_count(1, src_p)
          row = recv_meta(rm_disp(src_p)+metalen*(meta_l-1))
          tr = row .LT. 0
          stored_row_new = ABS(row)
          stored_col_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+1)
          nze = row_blk_size_new(stored_row_new) * col_blk_size_new(stored_col_new)
          !r_dp => recv_data(blk_ps:blk_ps+nze-1)
          !CALL dbcsr_put_block(redist, stored_row_new, stored_col_new, r_dp, tr)
          !### this should be changed to be like the make images (i.e., copy data in finalize, not here & now)
          data_block = pointer_view (data_block, recv_data, blk_ps, nze)
          CALL dbcsr_put_block(redist, stored_row_new, stored_col_new, data_block, tr)
          blk_ps = blk_ps + nze
          blks = blks + 1
       ENDDO
    ENDDO
    CALL dbcsr_data_clear_pointer (data_block)
    CALL dbcsr_data_release (data_block)
    !
    IF (dbg) THEN
       WRITE(*,*)routineN//" Declared blocks=",redist%m%wms(1)%lastblk,&
            "actual=",blks
       WRITE(*,*)routineN//" Declared data size=",redist%m%wms(1)%datasize,&
            "actual=",blk_ps
    ENDIF
    CALL dbcsr_finalize(redist, error=error)
    CALL dbcsr_data_release (recv_data)
    CALL dbcsr_data_release (send_data)
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_meta)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_redistribute

! *****************************************************************************
!> \brief Redistributes data blocks of a DBCSR matrix read from a file. 
!>        This routine should be used with dbcsr_binary_read in the module 
!>        dbcsr_io.F 
!>
!> \param[in] dblk, row_p, col_i, blk_p   data blocks, row_p, col_i and blk_p
!>    of the DBCSR matrix that the current node possesses after reading the 
!>    data file 
!> \param[in] proc_nblks, proc_darea_sizes   1D arrays holding nblks and data_
!>    area_size of those nodes of the mp environment, that created the file, 
!>    whose contents have been read by the current node of the present mp environment
!> \param[in,out] new_matrix   redistributed matrix 
!> 
!> \par History
!>      11.2012 created [Hossein Bani-Hashemian]
! *****************************************************************************
  SUBROUTINE dbcsr_datablock_redistribute (dblk, row_p, col_i, blk_p,&
                                           proc_nblks, proc_darea_sizes, new_matrix,error)

    TYPE(dbcsr_data_obj), INTENT(IN)         :: dblk
    INTEGER, DIMENSION(:), INTENT(IN), &
      POINTER                                :: row_p, col_i, blk_p, &
                                                proc_nblks, proc_darea_sizes
    TYPE(dbcsr_obj), INTENT(INOUT)           :: new_matrix
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_datablock_redistribute', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: metalen = 2

    COMPLEX(kind=dp), DIMENSION(:), POINTER  :: c_dp
    COMPLEX(kind=sp), DIMENSION(:), POINTER  :: c_sp
    INTEGER :: bcol, blk, blk_ps, blk_size, blks, brow, col_size, data_type, &
      dst_p, error_handler, i, ind, job_count, meta_l, mp_group, &
      nblkrows_total, numnodes, row_size, src_p, stored_col_new, &
      stored_row_new
    INTEGER(kind=int_8)                      :: actual_blk, blkp, end_ind, &
                                                start_ind
    INTEGER(kind=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: extra_darea_size, extra_nblks
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rd_disp, recv_meta, rm_disp, &
                                                sd_disp, sdp, send_meta, &
                                                sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_count, send_count, &
                                                total_recv_count, &
                                                total_send_count
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: sym_tr, tr
    REAL(kind=dp), DIMENSION(:), POINTER     :: r_dp
    REAL(kind=sp), DIMENSION(:), POINTER     :: r_sp
    TYPE(dbcsr_data_obj)                     :: data_block, recv_data, &
                                                send_data
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(dbcsr_mp_obj)                       :: mp_env

    CALL dbcsr_error_set(routineN, error_handler, error)


     dist = dbcsr_distribution (new_matrix)
     mp_env = dbcsr_distribution_mp (dist)
     numnodes = dbcsr_mp_numnodes (mp_env)
     mp_group = dbcsr_mp_group (mp_env)

     data_type = dbcsr_get_data_type(new_matrix)
     nblkrows_total = dbcsr_nblkrows_total(new_matrix)
     row_blk_size => array_data (new_matrix%m%row_blk_size)
     col_blk_size => array_data (new_matrix%m%col_blk_size)

     ALLOCATE (send_count(2, 0:numnodes-1))
     ALLOCATE (recv_count(2, 0:numnodes-1))
     ALLOCATE (total_send_count(2, 0:numnodes-1))
     ALLOCATE (total_recv_count(2, 0:numnodes-1))
     ALLOCATE (sdp(0:numnodes-1))
     ALLOCATE (sd_disp(0:numnodes-1))
     ALLOCATE (smp(0:numnodes-1))
     ALLOCATE (sm_disp(0:numnodes-1))
     ALLOCATE (rd_disp(0:numnodes-1))
     ALLOCATE (rm_disp(0:numnodes-1))
     send_count(:,:) = 0
     dst_p = -1

     job_count = COUNT(proc_nblks .NE. 0)

     ALLOCATE(extra_nblks(job_count))
     ALLOCATE(extra_darea_size(job_count))
     CALL cumsum_l(INT((/0,proc_nblks(1:job_count-1)/),kind=int_8),extra_nblks)
     CALL cumsum_l(INT((/0,proc_darea_sizes(1:job_count-1)/),kind=int_8),extra_darea_size)

     i = 0
     DO ind = 1, job_count*nblkrows_total
        brow = MOD(ind-1,nblkrows_total)+1
        IF (brow .EQ. 1) i = i+1
        row_size = row_blk_size(brow)
        DO blk = row_p(ind+i-1)+1, row_p(ind+i)
           actual_blk = INT(blk,kind=int_8)+extra_nblks(i)
           bcol = col_i(actual_blk)
           col_size = col_blk_size(bcol)
           blk_size = row_size*col_size
           sym_tr = .FALSE.
           CALL dbcsr_get_stored_coordinates (new_matrix%m, brow, bcol, sym_tr, dst_p)
           send_count(1, dst_p) = send_count(1, dst_p) + 1
           send_count(2, dst_p) = send_count(2, dst_p) + blk_size
        ENDDO
     ENDDO
     CALL mp_alltoall(send_count, recv_count, 2, mp_group)
     CALL dbcsr_data_init (recv_data)
     CALL dbcsr_data_new (recv_data, data_type, SUM(recv_count(2, :)))
     ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :))))
     CALL dbcsr_data_init (send_data)
     CALL dbcsr_data_new (send_data, data_type, SUM(send_count(2, :)))
     ALLOCATE (send_meta(metalen*SUM(send_count(1, :))))
     DO dst_p = 0, numnodes-1
        total_send_count(1, dst_p) = send_count (1, dst_p)
        total_send_count(2, dst_p) = send_count (2, dst_p)
        total_recv_count(1, dst_p) = recv_count (1, dst_p)
        total_recv_count(2, dst_p) = recv_count (2, dst_p)
     ENDDO
     sd_disp = -1 ; sm_disp = -1 ; rd_disp = -1 ; rm_disp = -1                                                             
     sd_disp(0) = 1 ; sm_disp(0) = 1 ; rd_disp(0) = 1 ; rm_disp(0) = 1                        
     DO dst_p = 1, numnodes-1                               
        sm_disp(dst_p) = sm_disp(dst_p-1) + metalen*total_send_count(1, dst_p-1)
        sd_disp(dst_p) = sd_disp(dst_p-1) + total_send_count(2, dst_p-1)
        rm_disp(dst_p) = rm_disp(dst_p-1) + metalen*total_recv_count(1, dst_p-1)
        rd_disp(dst_p) = rd_disp(dst_p-1) + total_recv_count(2, dst_p-1)
     ENDDO
     sdp(:) = sd_disp     
     smp(:) = sm_disp - metalen  
    
     SELECT CASE (data_type)
     CASE (dbcsr_type_real_4)
          r_sp => dblk%d%r_sp
     CASE (dbcsr_type_real_8)
          r_dp => dblk%d%r_dp
     CASE (dbcsr_type_complex_4)
          c_sp => dblk%d%c_sp
     CASE (dbcsr_type_complex_8)
          c_dp => dblk%d%c_dp
     END SELECT 

     CALL dbcsr_data_init (data_block)
     CALL dbcsr_data_new (data_block, data_type)
     i = 0
     dst_p = -1
     DO ind = 1, job_count*nblkrows_total
        brow = MOD(ind-1,nblkrows_total)+1
        IF (brow .EQ. 1) i = i+1
        row_size = row_blk_size(brow)
        DO blk = row_p(ind+i-1)+1, row_p(ind+i)
           actual_blk = INT(blk,kind=int_8)+extra_nblks(i)
           bcol = col_i(actual_blk)
           col_size = col_blk_size(bcol)
           blk_size = row_size*col_size
           blkp = INT(blk_p(actual_blk),kind=int_8)
           start_ind = blkp+extra_darea_size(i)
           end_ind = blkp+extra_darea_size(i)+blk_size-1

           SELECT CASE (data_type)
           CASE (dbcsr_type_real_4)
               data_block%d%r_sp => r_sp(start_ind:end_ind)
           CASE (dbcsr_type_real_8)
               data_block%d%r_dp => r_dp(start_ind:end_ind)
           CASE (dbcsr_type_complex_4)
               data_block%d%c_sp => c_sp(start_ind:end_ind)
           CASE (dbcsr_type_complex_8)
               data_block%d%c_dp => c_dp(start_ind:end_ind)
           END SELECT
           sym_tr = .FALSE.
           CALL dbcsr_get_stored_coordinates(new_matrix%m, brow, bcol, sym_tr, dst_p)
           smp(dst_p) = smp(dst_p) + metalen
           tr=.FALSE.
!           IF (tr) THEN
!              send_meta(smp(dst_p)) = -brow
!           ELSE
              send_meta(smp(dst_p)) = brow
!           ENDIF
           send_meta(smp(dst_p)+1) = bcol 
           blk_size = dbcsr_data_get_size (data_block)

           CALL dbcsr_data_set (send_data, lb=sdp(dst_p),&
                                data_size=blk_size, src=data_block, source_lb=1)
           sdp(dst_p) = sdp(dst_p) + blk_size
        ENDDO
     ENDDO
     CALL dbcsr_data_clear_pointer (data_block)

     SELECT CASE (data_type)
     CASE (dbcsr_type_real_4)
        CALL hybrid_alltoall_s1(&
            send_data%d%r_sp(:), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%r_sp(:), total_recv_count(2,:), rd_disp(:)-1,&
            mp_env)
     CASE (dbcsr_type_real_8)
        CALL hybrid_alltoall_d1 (&
           send_data%d%r_dp(:), total_send_count(2,:), sd_disp(:)-1,&
           recv_data%d%r_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
           mp_env)
     CASE (dbcsr_type_complex_4)
        CALL hybrid_alltoall_c1(&
           send_data%d%c_sp(:), total_send_count(2,:), sd_disp(:)-1,&
           recv_data%d%c_sp(:), total_recv_count(2,:), rd_disp(:)-1,&
           mp_env)
     CASE (dbcsr_type_complex_8)
        CALL hybrid_alltoall_z1(&
          send_data%d%c_dp(:), total_send_count(2,:), sd_disp(:)-1,&
          recv_data%d%c_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
          mp_env)
     END SELECT         
     CALL hybrid_alltoall_i1(send_meta(:), metalen*total_send_count(1,:), sm_disp(:)-1,&
          recv_meta(:), metalen*total_recv_count(1,:), rm_disp(:)-1, mp_env)
     CALL dbcsr_work_create(new_matrix,SUM(recv_count(1,:)),&
                            SUM(recv_count(2,:)), work_mutable=.FALSE., n=1, error=error)

     blk_ps = 1
     blks = 0
     DO src_p = 0, numnodes-1
        DO meta_l = 1, recv_count(1, src_p)
           brow = recv_meta(rm_disp(src_p)+metalen*(meta_l-1))
           tr = brow .LT. 0
           stored_row_new = ABS(brow)
           stored_col_new = recv_meta(rm_disp(src_p)+metalen*(meta_l-1)+1)
           blk_size = row_blk_size(stored_row_new) * col_blk_size(stored_col_new)

           data_block = pointer_view (data_block, recv_data, blk_ps, blk_size)
           CALL dbcsr_put_block(new_matrix, stored_row_new, stored_col_new, data_block, tr)

           blk_ps = blk_ps + blk_size
           blks = blks + 1
        ENDDO
     ENDDO
     CALL dbcsr_data_clear_pointer (data_block)
     DEALLOCATE (data_block%d)
     CALL dbcsr_finalize(new_matrix, reshuffle=.TRUE., error=error)
     CALL dbcsr_data_release (recv_data)
     CALL dbcsr_data_release (send_data)

     DEALLOCATE(send_count)
     DEALLOCATE(recv_count)
     DEALLOCATE(sdp); DEALLOCATE(sd_disp)
     DEALLOCATE(smp); DEALLOCATE(sm_disp)
     DEALLOCATE(rd_disp)
     DEALLOCATE(rm_disp)
     DEALLOCATE(recv_meta)
     DEALLOCATE(send_meta)
     DEALLOCATE(extra_nblks); DEALLOCATE(extra_darea_size)

    CALL dbcsr_error_stop(error_handler, error)
    CONTAINS
    
    SUBROUTINE cumsum_l(arr,cumsum)
    INTEGER(kind=int_8), DIMENSION(:), &
      INTENT(IN)                             :: arr
    INTEGER(kind=int_8), &
      DIMENSION(SIZE(arr)), INTENT(OUT)      :: cumsum

    INTEGER                                  :: i

        cumsum(1)=arr(1)
        DO i=2,SIZE(arr)
          cumsum(i)=cumsum(i-1)+arr(i)
        END DO
    END SUBROUTINE cumsum_l
  END SUBROUTINE dbcsr_datablock_redistribute

END MODULE dbcsr_transformations

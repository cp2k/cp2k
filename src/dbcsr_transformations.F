!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
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
                                             array_new,&
                                             array_nullify,&
                                             array_release
  USE dbcsr_block_access,              ONLY: dbcsr_get_block_p,&
                                             dbcsr_put_block
  USE dbcsr_block_operations,          ONLY: dbcsr_block_partial_copy,&
                                             dbcsr_block_transpose,&
                                             dbcsr_data_clear,&
                                             dbcsr_data_set
  USE dbcsr_data_operations,           ONLY: dbcsr_copy_sort_data
  USE dbcsr_dist_operations,           ONLY: dbcsr_get_stored_coordinates,&
                                             dbcsr_reblocking_targets,&
                                             dbcsr_transpose_dims,&
                                             dbcsr_transpose_distribution,&
                                             find_local_virtual_elements,&
                                             make_image_distribution_dense,&
                                             match_sizes_to_dist
  USE dbcsr_error_handling
  USE dbcsr_index_operations,          ONLY: dbcsr_addto_index_array,&
                                             dbcsr_clearfrom_index_array,&
                                             dbcsr_repoint_index,&
                                             make_dense_index,&
                                             make_undense_index,&
                                             transpose_index_local
  USE dbcsr_io,                        ONLY: dbcsr_print
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_kinds,                     ONLY: dp,&
                                             sp
  USE dbcsr_machine,                   ONLY: m_walltime
  USE dbcsr_message_passing,           ONLY: mp_alltoall
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_sizes, dbcsr_data_clear_pointer, &
       dbcsr_data_ensure_size, dbcsr_data_hold, dbcsr_data_init, &
       dbcsr_data_new, dbcsr_data_release, dbcsr_data_set_pointer, &
       dbcsr_distribution, dbcsr_distribution_col_dist, &
       dbcsr_distribution_hold, dbcsr_distribution_init, &
       dbcsr_distribution_local_cols, dbcsr_distribution_local_rows, &
       dbcsr_distribution_mp, dbcsr_distribution_ncols, &
       dbcsr_distribution_new, dbcsr_distribution_nlocal_cols, &
       dbcsr_distribution_nlocal_rows, dbcsr_distribution_nrows, &
       dbcsr_distribution_release, dbcsr_distribution_row_dist, &
       dbcsr_get_data, dbcsr_get_data_size, dbcsr_get_data_type, &
       dbcsr_get_matrix_type, dbcsr_init, dbcsr_is_initialized, &
       dbcsr_mp_group, dbcsr_mp_mynode, dbcsr_mp_mypcol, dbcsr_mp_myprow, &
       dbcsr_mp_npcols, dbcsr_mp_nprows, dbcsr_mp_numnodes, dbcsr_mp_pgrid, &
       dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_nfullcols_total, &
       dbcsr_nfullrows_total, dbcsr_release, dbcsr_row_block_sizes, &
       dbcsr_type_1d_to_2d, dbcsr_uses_special_memory, dbcsr_valid_index
  USE dbcsr_mp_operations,             ONLY: hybrid_alltoall_c1,&
                                             hybrid_alltoall_d1,&
                                             hybrid_alltoall_i1,&
                                             hybrid_alltoall_s1,&
                                             hybrid_alltoall_z1
  USE dbcsr_ptr_util,                  ONLY: pointer_view
  USE dbcsr_types,                     ONLY: &
       dbcsr_2d_array_type, dbcsr_data_obj, dbcsr_distribution_obj, &
       dbcsr_imagedistribution_type, dbcsr_iterator, dbcsr_meta_size, &
       dbcsr_mp_obj, dbcsr_obj, dbcsr_repl_full, dbcsr_repl_none, &
       dbcsr_slot_blk_p, dbcsr_slot_col_i, dbcsr_slot_nblks, dbcsr_slot_nze, &
       dbcsr_slot_row_p, dbcsr_type, dbcsr_type_antisymmetric, &
       dbcsr_type_complex_4, dbcsr_type_complex_8, dbcsr_type_normal, &
       dbcsr_type_real_4, dbcsr_type_real_8, dbcsr_type_symmetric
  USE dbcsr_util,                      ONLY: absmod,&
                                             convert_sizes_to_offsets,&
                                             dbcsr_pack_meta,&
                                             dbcsr_set_debug,&
                                             dbcsr_unpack_meta,&
                                             global_offsets_to_local,&
                                             nfull_elements,&
                                             sgn
  USE dbcsr_work_operations,           ONLY: dbcsr_create,&
                                             dbcsr_finalize,&
                                             dbcsr_work_create
  USE f77_blas
  USE f77_blas_generic

  !$ USE OMP_LIB

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_transformations'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034


  PUBLIC :: dbcsr_desymmetrize_deep,&
            dbcsr_new_transposed,&
            dbcsr_complete_redistribute,&
            dbcsr_redistribute, dbcsr_make_untransposed_blocks
  PUBLIC :: dbcsr_replicate_all, dbcsr_distribute
  PUBLIC :: dbcsr_make_dense, dbcsr_make_undense, dbcsr_make_images_dense
  PUBLIC :: dbcsr_make_images



#define DBG IF (dbg) WRITE(*,*)routineN//" ",
#define DBGV IF (bcsr_verbose) WRITE(*,*)routineN//" ",
#define DBGI IF (info) WRITE(*,*)routineN//" ",
#define DEBUG_HEADER        LOGICAL :: dbg, info
#define DEBUG_BODY        dbg = .FALSE. ; CALL dbcsr_set_debug(dbg, info=info)
#define DEBUG_BODYY       dbg = .TRUE. ; CALL dbcsr_set_debug(dbg, info=info)

  LOGICAL, PARAMETER :: bcsr_debug =   .TRUE.
  LOGICAL, PARAMETER :: bcsr_info =    .FALSE.
  LOGICAL, PARAMETER :: bcsr_verbose = .FALSE.

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
!> \param[in] use_distribution          (optional) use this distribution
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_new_transposed (transposed, normal, shallow_data_copy,&
       transpose_data, transpose_distribution, use_distribution, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: transposed
    TYPE(dbcsr_obj), INTENT(IN)              :: normal
    LOGICAL, INTENT(IN), OPTIONAL            :: shallow_data_copy, &
                                                transpose_data, &
                                                transpose_distribution
    TYPE(dbcsr_distribution_obj), &
      INTENT(IN), OPTIONAL                   :: use_distribution
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_new_transposed', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: new_type
    INTEGER                                  :: error_handler, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: blk_p
    LOGICAL                                  :: shallow, tr_blocks, tr_dist
    TYPE(dbcsr_distribution_obj)             :: new_dist
    TYPE(dbcsr_obj)                          :: t2

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (dbcsr_valid_index(normal), dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Matrix does not exist.",__LINE__,error)
    ! Internalize options
    shallow = .FALSE.
    IF (PRESENT (shallow_data_copy)) shallow = shallow_data_copy
    tr_blocks = .TRUE.
    IF (PRESENT (transpose_data)) tr_blocks = transpose_data
    tr_dist = .TRUE.
    IF (PRESENT (transpose_distribution)) tr_dist = transpose_distribution
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
    new_type = normal%m%matrix_type
    CALL dbcsr_create (transposed, "transposed "//normal%m%name, new_dist,&
         new_type, normal%m%col_blk_size, normal%m%row_blk_size,&
         data_type=normal%m%data_type, special=normal%m%special_memory, error=error)
    CALL dbcsr_distribution_release (new_dist)
    ! Reserve the space for the new indices.
    CALL dbcsr_addto_index_array (transposed%m, dbcsr_slot_row_p,&
         reservation=transposed%m%nblkrows_total+1, extra=transposed%m%nblks*2,&
         error=error)
    CALL dbcsr_addto_index_array (transposed%m, dbcsr_slot_col_i,&
         reservation=normal%m%nblks, error=error)
    CALL dbcsr_addto_index_array (transposed%m, dbcsr_slot_blk_p,&
         reservation=normal%m%nblks, error=error)
    CALL dbcsr_repoint_index (transposed%m)
    IF (.NOT. shallow) THEN
       CALL dbcsr_data_ensure_size (transposed%m%data_area, normal%m%nze,&
            nocopy=.TRUE., error=error)
    ENDIF
    !
    transposed%m%nblks = normal%m%nblks
    transposed%m%nze = normal%m%nze
    transposed%m%index(dbcsr_slot_nblks) = normal%m%nblks
    transposed%m%index(dbcsr_slot_nze) = normal%m%nze
    ! Transpose the local index.
    DBG 'Original size', normal%m%nblks, normal%m%nze
    DBG 'meta', transposed%m%nblkrows_total+1, normal%m%nblks
    DBG 'size', SIZE (transposed%m%row_p), SIZE (transposed%m%col_i)
    ALLOCATE (blk_p (normal%m%nblks), stat=stat)
    CALL dbcsr_assert(stat == 0, dbcsr_fatal_level, dbcsr_caller_error, routineN,&
         "blk_p",__LINE__,error)
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
    ! Copy the data
    IF (shallow) THEN
       CALL dbcsr_data_release (transposed%m%data_area)
       transposed%m%data_area = normal%m%data_area
       CALL dbcsr_data_hold (transposed%m%data_area)
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
    DEALLOCATE (blk_p)
    !
    ! Redistribute if mpgrid was not transposed
    IF (.NOT. tr_dist) THEN
       !write (*,*)routineN//" redistributing"
       CALL dbcsr_init (t2)
       CALL dbcsr_create (t2, transposed%m%name, transposed%m%dist,&
            dbcsr_type_normal, transposed%m%row_blk_size, transposed%m%col_blk_size,&
            0, 0, transposed%m%data_type, error=error)
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
!> \param[in,out] error       cp2k error
!>
!>  THIS ROUTINE IS NOT USED ANYWHERE !
!>  THIS ROUTINE IS NOT USED ANYWHERE !
!>  THIS ROUTINE IS NOT USED ANYWHERE !
!>  THIS ROUTINE IS NOT USED ANYWHERE !
!>
! *****************************************************************************
  SUBROUTINE dbcsr_desymmetrize_deep(sm, desm, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: sm
    TYPE(dbcsr_obj), INTENT(OUT)             :: desm
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_desymmetrize_deep', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 2

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: c_dp
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: c_sp
    INTEGER :: blk, blk_l, blk_p, blk_ps, blks, col, dst_p, error_handler, &
      mp_group, nsymmetries, numproc, nze, pcol, prow, row, src_p, &
      stored_col, stored_row, symmetry_i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rd_disp, recv_meta, rm_disp, &
                                                sd_disp, sdp, send_meta, &
                                                sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_count, send_count, &
                                                total_recv_count, &
                                                total_send_count
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, col_dist, &
                                                row_blk_size, row_dist
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    LOGICAL                                  :: tr
    REAL(KIND=dp)                            :: tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: send_data
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: recv_data
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_dp
    REAL(KIND=sp), DIMENSION(:), POINTER     :: r_sp
    TYPE(dbcsr_distribution_obj)             :: target_dist
    TYPE(dbcsr_mp_obj)                       :: mp_obj

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    tstart = m_walltime ()
    CALL dbcsr_assert (dbcsr_valid_index (sm),&
         dbcsr_fatal_level, dbcsr_caller_error, routineN,&
         "Matrix not initialized.",__LINE__,error)
    nsymmetries = 1
    IF (sm%m%matrix_type .EQ. 'S' .OR. sm%m%matrix_type .EQ. 'A') THEN
       nsymmetries = 2
    ENDIF
    SELECT CASE (sm%m%data_type)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (sm%m%data_area, r_dp)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (sm%m%data_area, r_sp)
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
            routineN,"Only real double precision",__LINE__,error)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (sm%m%data_area, c_dp)
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
            routineN,"Only real double precision",__LINE__,error)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (sm%m%data_area, c_sp)
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
            routineN,"Only real double precision",__LINE__,error)
    END SELECT
    row_blk_size => array_data (sm%m%row_blk_size)
    col_blk_size => array_data (sm%m%col_blk_size)
    target_dist = sm%m%dist
    row_dist => array_data (dbcsr_distribution_row_dist (target_dist))
    col_dist => array_data (dbcsr_distribution_col_dist (target_dist))
    mp_obj = dbcsr_distribution_mp (target_dist)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    numproc = dbcsr_mp_numnodes (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    DBGV 'row_dist', row_dist
    DBGV 'col_dist', col_dist
    DBG 'row, col maxval:', MAXVAL(row_dist), MAXVAL(col_dist)
    DBG 'blacs2mpi',blacs2mpi
    DBG 'blacs2mpi LB',LBOUND(blacs2mpi)
    DBG 'blacs2mpi UB',UBOUND(blacs2mpi)
    IF (sm%m%matrix_type .EQ. 'S' .OR. sm%m%matrix_type .EQ. 'A') THEN
       CALL dbcsr_assert(SIZE(row_dist) .EQ. SIZE(col_dist), dbcsr_warning_level,&
            dbcsr_wrong_args_error, routineN,&
            'Unequal row and column distributions for symmetric matrix.',__LINE__,error)
    ENDIF
    CALL dbcsr_assert(MAXVAL(row_dist).LE.UBOUND(blacs2mpi,1), dbcsr_failure_level,&
         dbcsr_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows',__LINE__,error)
    IF (dbg) &
         CALL dbcsr_assert(MAXVAL(row_dist).EQ.UBOUND(blacs2mpi,1), dbcsr_warning_level,&
         dbcsr_wrong_args_error, routineN,&
         'Range of row distribution not equal to processor rows',__LINE__,error)
    CALL dbcsr_assert(MAXVAL(col_dist).LE.UBOUND(blacs2mpi,2), dbcsr_failure_level,&
         dbcsr_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols',__LINE__,error)
    IF (dbg) &
         CALL dbcsr_assert(MAXVAL(col_dist).EQ.UBOUND(blacs2mpi,2), dbcsr_warning_level,&
         dbcsr_wrong_args_error, routineN,&
         'Range of col distribution not equal to processor cols',__LINE__,error)
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
    CALL dbcsr_create(desm, 'desym '//sm%m%name, sm%m%dist,&
         dbcsr_type_normal, sm%m%row_blk_size, sm%m%col_blk_size,&
         0,0, dbcsr_get_data_type(sm), special=sm%m%special_memory,&
         error=error)
    IF (sm%m%matrix_type .EQ. dbcsr_type_antisymmetric) THEN
       desm%m%sym_negation = .TRUE.
    ENDIF
    ! Count initial sizes for sending.
    DBG 'A'
    send_count(:,:) = 0
    DO row = 1, sm%m%nblkrows_total
       DO blk = sm%m%row_p(row)+1, sm%m%row_p(row+1)
          col = sm%m%col_i(blk)
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
             nze = row_blk_size(stored_row) * col_blk_size(stored_col)
             send_count(1, dst_p) = send_count(1, dst_p) + 1
             send_count(2, dst_p) = send_count(2, dst_p) + nze
          ENDDO ! symmetry_i
       ENDDO ! col_i
    ENDDO ! row
    DBG 'B'
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    DBG 'C'
    DBG 'send counts',send_count
    DBG 'recv counts',recv_count
    ! Allocate data structures needed for data exchange.
    ALLOCATE (recv_data(SUM(recv_count(2, :))))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :))))
    ALLOCATE (send_data(SUM(send_count(2, :))))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :))))
    DBG 'send data len',SIZE (send_data)
    DBG 'send meta len',SIZE (send_meta)
    DBG 'recv data len',SIZE (recv_data)
    DBG 'recv meta len',SIZE (recv_meta)
    ! Fill in the meta data structures and copy the data.
    DO dst_p = 0, numproc-1
       total_send_count(1, dst_p) = send_count (1, dst_p)
       total_send_count(2, dst_p) = send_count (2, dst_p)
       total_recv_count(1, dst_p) = recv_count (1, dst_p)
       total_recv_count(2, dst_p) = recv_count (2, dst_p)
    ENDDO
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
    !DO dst_p = 0, numproc-1
    !      sm_disp (dst_p) = sm_disp (dst_p) + metalen*send_count(1, dst_p)
    !      sd_disp (dst_p) = sd_disp (dst_p) + send_count(2, dst_p)
    !      rm_disp (dst_p) = rm_disp (dst_p) + metalen*recv_count(1, dst_p)
    !      rd_disp (dst_p) = rd_disp (dst_p) + recv_count(2, dst_p)
    !ENDDO
    sdp(:) = sd_disp
    smp(:) = sm_disp
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    DBG 'D'
    DO row = 1, sm%m%nblkrows_total
       DO blk = sm%m%row_p(row)+1, sm%m%row_p(row+1)
          col = sm%m%col_i(blk)
          blk_p = sm%m%blk_p(blk)
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
             nze = row_blk_size(stored_row) * col_blk_size(stored_col)
             send_meta(smp(dst_p)) = stored_row
             send_meta(smp(dst_p)+1) = stored_col
             DBG stored_row, stored_col, blk_p, tr, nze
             IF (.NOT. tr) THEN
                send_data(sdp(dst_p):sdp(dst_p)+nze-1) =&
                     r_dp(blk_p:blk_p+nze-1)
                !send_data(sdp(dst_p):sdp(dst_p)+nze-1) =&
                !     sm%m%data(blk_p:blk_p+nze-1)
             ELSE
                send_meta(smp(dst_p)) = -stored_row
                send_data(sdp(dst_p):sdp(dst_p)+nze-1) =&
                     r_dp(blk_p:blk_p+nze-1)
                !send_data(sdp(dst_p):sdp(dst_p)+nze-1) =&
                !     sm%m%data(blk_p:blk_p+nze-1)
             ENDIF
             smp(dst_p) = smp(dst_p) + metalen
             sdp(dst_p) = sdp(dst_p) + nze
          ENDDO ! symmetry_i
       ENDDO ! col_i
    ENDDO ! row
    ! Exchange the data and metadata structures.
    DBG 'E'
    CALL mp_alltoall(send_data(:), total_send_count(2,:), sd_disp(:)-1,&
         recv_data(:), total_recv_count(2,:), rd_disp(:)-1, mp_group)
    DBG 'F'
    DBG 'send_meta',send_meta
    CALL mp_alltoall(send_meta(:), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), metalen*total_recv_count(1,:), rm_disp(:)-1, mp_group)
    DBG 'recv_meta',recv_meta
    DBG 'G'
    ! Now fill in the data.
    CALL dbcsr_work_create(desm,&
            SUM(recv_count(1,:)),&
            SUM(recv_count(2,:)), n=1, error=error)
    !DBG 'send disps data',sd_disp
    !DBG 'send disps meta',sm_disp
    !DBG 'recv disps data',rd_disp
    !DBG 'recv disps meta',rm_disp
    !blk_p = 1
    !blk = 1
    blk_ps = 1
    blks = 1
    DO src_p = 0, numproc-1
       !wm%data(blk_p : blk_p+recv_count(2,src_p)-1) =&
       !     recv_data(rd_disp(src_p) :&
       !               rd_disp(src_p)+recv_count(2,src_p)-1)
       nze = recv_count(2, src_p)
       CALL dbcsr_data_set (desm%m%wms(1)%data_area, blk_ps, nze,&
            recv_data, rd_disp(src_p))
       !desm%m%w%data(blk_ps:blk_ps+nze-1) =&
       !     recv_data(rd_disp(src_p):rd_disp(src_p)+nze-1)
       DBG 'processing from',src_p,'counts',send_count(1, src_p)
       DO blk_l = 1, recv_count(1, src_p)
          stored_row = recv_meta(rm_disp(src_p)+metalen*(blk_l-1))
          stored_col = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+1)
          !DBG 'blk: p',blk,src_p,'; row, col',stored_row,stored_col,&
          !     '; blk_p',blk_p,'; row_rep, imgdist',row_rep,&
          !     row_rep_dist (ABS (stored_row))
          DBG 'blk: p',blks,src_p,'; row, col',&
               stored_row,stored_col,&
               '; blk_ps',blk_ps
          !row_rep = row_rep_dist (ABS(stored_row))
          desm%m%wms(1)%row_i(blks) = ABS(stored_row)
          desm%m%wms(1)%col_i(blks) = stored_col
          desm%m%wms(1)%blk_p(blks) = SIGN(blk_ps, stored_row)
          nze = row_blk_size(ABS(stored_row))&
               * col_blk_size(stored_col)
          !blk_p = blk_p + nze
          !blk = blk + 1
          blk_ps = blk_ps + nze
          blks = blks + 1
       ENDDO
    ENDDO
    !
    desm%m%wms(1)%lastblk = blks - 1
    desm%m%wms(1)%datasize = blk_ps - 1
    DBG 'Finalizing normalization'
    CALL dbcsr_finalize(desm, error=error)
    tstop = m_walltime ()
    DBG 'time:',tstop-tstart
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_data)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_data)
    DEALLOCATE(send_meta)
    DBG 'desym name',desm%m%name
    DBG 'Done desymmetrizing!'
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
                                                mp_group, mynode, nblks, &
                                                numnodes, nze, p, row
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

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_valid_index(matrix%m),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Matrix not initialized.",__LINE__,error)
    CALL dbcsr_assert (matrix%m%replication_type .NE. dbcsr_repl_none,&
         dbcsr_warning_level, dbcsr_caller_error, routineN,&
         "Distributing a non-replicated matrix makes no sense.",__LINE__,error)
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
    numnodes = dbcsr_mp_numnodes (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
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
       CALL dbcsr_create(distributed, 'Distributed '//matrix%m%name,&
            matrix%m%dist,&
            dbcsr_type_normal, matrix%m%row_blk_size, matrix%m%col_blk_size,&
            0, 0, matrix%m%data_type, special=matrix%m%special_memory, &
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
       tmp_data = matrix%m%data_area
       matrix%m%data_area = distributed%m%data_area
       distributed%m%data_area = tmp_data
       tmp_index => matrix%m%index
       matrix%m%index => distributed%m%index
       distributed%m%index => tmp_index
       CALL dbcsr_repoint_index (matrix%m)
       matrix%m%nze = distributed%m%nze
       matrix%m%nblks = distributed%m%nblks
       CALL dbcsr_release (distributed)
    ENDIF
    DBG 'Done replicating'
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_distribute


! *****************************************************************************
!> \brief Detransposes all blocks in a matrix
!> \param[in,out] matrix      DBCSR matrix
! *****************************************************************************
  SUBROUTINE dbcsr_make_untransposed_blocks(matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_make_untransposed_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_size, row, &
                                                row_size
    INTEGER, DIMENSION(:), POINTER           :: cbs, rbs
    LOGICAL                                  :: sym_negation, tr
    TYPE(dbcsr_data_obj)                     :: block_data
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    rbs => array_data (dbcsr_row_block_sizes (matrix))
    cbs => array_data (dbcsr_col_block_sizes (matrix))
    CALL dbcsr_data_init (block_data)
    CALL dbcsr_data_new (block_data, dbcsr_get_data_type(matrix))
    CALL dbcsr_iterator_start (iter, matrix)
    sym_negation = matrix%m%sym_negation
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, block_data,&
            transposed=tr,&
            block_number=blk)
       IF (tr) THEN
          row_size = rbs(row)
          col_size = cbs(col)
          CALL dbcsr_block_transpose(block_data, col_size, row_size)
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
  END SUBROUTINE dbcsr_make_untransposed_blocks



! *****************************************************************************
!> \brief Creates row and column images of a matrix.
!> \param[in] source          input matrix
!> \param[in,out] normalized  image array of the normalized matrix
!> \param[in] target_image_dist          normalize to this image distribution
!> \param[in] predistribute   (optional) predistribute data for multiplication
!> \param[in] special_memory  (optional) whether to use special memory;
!>                            default is yes
!> \param[in] no_copy_data    (optional) try to not merge data at the end
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_make_images(source, normalized, target_image_dist,&
       predistribute, special_memory, no_copy_data, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: source
    TYPE(dbcsr_2d_array_type), INTENT(OUT)   :: normalized
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(IN)                             :: target_image_dist
    CHARACTER, INTENT(IN), OPTIONAL          :: predistribute
    LOGICAL, INTENT(IN), OPTIONAL            :: special_memory, no_copy_data
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_images', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_valid_index (source),&
         dbcsr_fatal_level, dbcsr_caller_error, routineN,&
         "Matrix not initialized.",__LINE__,error)
    DBG 'Make images for matrix', source%m%matrix_type, source%m%name
    normalized%image_dist = target_image_dist
    SELECT CASE (source%m%matrix_type)
    CASE ('N', 'S', 'A')
       CALL make_images(source, normalized,&
            target_image_dist, desymmetrize=.TRUE.,&
            predistribute=predistribute, special_memory=special_memory,&
            no_copy_data=no_copy_data,&
            error=error)
    CASE ('T')
       CALL dbcsr_assert(.FALSE., dbcsr_failure_level,&
            dbcsr_unimplemented_error_nr, routineP,&
            "Making images of transposed matrices not "//&
            "supported yet.",__LINE__,error)
    CASE DEFAULT
       CALL dbcsr_assert(.FALSE., dbcsr_failure_level,&
            dbcsr_wrong_args_error, routineP, "Invalid matrix type.",__LINE__,error)
    END SELECT
    normalized%image_dist = target_image_dist
    DBG 'parent ums(1) name', normalized%mats(1,1)%m%name
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_make_images


! *****************************************************************************
!> \brief Makes column-based and row-based images of a matrix.
!> \param[in] sm    input symmetric matrix
!> \param[in,out] ums         normalied matrices
!> \param[in] target_imgdist  image distribution to normalize to
!> \param[in] desymmetrize    (optional) desymmetrize a symmetric matrix
!> \param[in] predistribute   (optional) predistribute data for multiplication
!> \param[in] special_memory  (optional) whether to use special memory;
!>                            default is yes
!> \param[in] no_copy_data    (optional) try to not merge data at the end
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE make_images(ism, ums, target_imgdist, desymmetrize, predistribute,&
       special_memory, no_copy_data, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: ism
    TYPE(dbcsr_2d_array_type), INTENT(OUT)   :: ums
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(IN)                             :: target_imgdist
    LOGICAL, INTENT(IN), OPTIONAL            :: desymmetrize
    CHARACTER, INTENT(IN), OPTIONAL          :: predistribute
    LOGICAL, INTENT(IN), OPTIONAL            :: special_memory, no_copy_data
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_images', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 5

    CHARACTER                                :: data_type
    INTEGER :: blk, blk_l, blk_p, bp, col, col_img, data_p, dst_p, &
      error_handler, i, ithread, j, mp_group, ncol_images, nrow_images, &
      nsymmetries, nthreads, numproc, nze, pcol, prev_blk_p, prev_dst_p, &
      prev_row, prow, row, row_img, sd_pos, sm_pos, src_p, stored_blk_p, &
      stored_col, stored_row, symmetry_i, vcol, vrow
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: lsdp, lsmp, rd_disp, &
                                                recv_meta, rm_disp, sd_disp, &
                                                sdp, send_meta, sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: blk_ps, blks, &
                                                total_recv_count, &
                                                total_send_count
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: row_dst_counts, &
                                                row_dst_offsets
    INTEGER, ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: myt_send_count, recv_count, &
                                                send_count
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, col_dist, &
                                                col_img_dist, row_blk_size, &
                                                row_dist, row_img_dist
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    LOGICAL                                  :: nocopy, predist, same_dst_p, &
                                                special, tr
    REAL(KIND=dp)                            :: tstart, tstop
    TYPE(dbcsr_data_obj)                     :: received_data_area, &
                                                recv_data_area, send_data_area
    TYPE(dbcsr_distribution_obj)             :: target_dist
    TYPE(dbcsr_error_type)                   :: dbcsr_error
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_type)                         :: sm

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, dbcsr_error)
    tstart = m_walltime ()
    special = .TRUE.
    IF (PRESENT (special_memory)) special=special_memory
    nocopy = .FALSE.
    IF (PRESENT (no_copy_data)) nocopy = no_copy_data
    sm = ism%m
    nsymmetries = 1
    IF (PRESENT (desymmetrize)) THEN
       IF (desymmetrize .AND. (&
            sm%matrix_type .EQ. dbcsr_type_symmetric&
            .OR. sm%matrix_type .EQ. dbcsr_type_antisymmetric)) THEN
          nsymmetries = 2
       ENDIF
    ENDIF
    predist = .FALSE.
    IF (PRESENT (predistribute)) predist = .TRUE.
    data_type = sm%data_type
    CALL dbcsr_assert (data_type .EQ. dbcsr_type_real_8&
         .or.data_type .EQ. dbcsr_type_real_4&
         .or.data_type .EQ. dbcsr_type_complex_8&
         .or.data_type .EQ. dbcsr_type_complex_4, dbcsr_fatal_level,&
         dbcsr_internal_error, routineN, "Invalid data type.",__LINE__,error)
    row_blk_size => array_data (sm%row_blk_size)
    col_blk_size => array_data (sm%col_blk_size)
    target_dist = target_imgdist%main
    row_dist => array_data (dbcsr_distribution_row_dist (target_dist))
    col_dist => array_data (dbcsr_distribution_col_dist (target_dist))
    nrow_images = target_imgdist%row_decimation
    IF (nrow_images .GT. 1) THEN
       row_img_dist => array_data (target_imgdist%row_image)
    ELSE
       NULLIFY (row_img_dist)
    ENDIF
    ncol_images = target_imgdist%col_decimation
    IF (ncol_images .GT. 1) THEN
       col_img_dist => array_data (target_imgdist%col_image)
    ELSE
       NULLIFY (col_img_dist)
    ENDIF
    mp_obj = dbcsr_distribution_mp (target_dist)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    numproc = dbcsr_mp_numnodes (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    DBG 'row_dist', row_dist
    DBG 'col_dist', col_dist
    DBG 'row, col maxval:', MAXVAL(row_dist), MAXVAL(col_dist)
    DBG 'blacs2mpi',blacs2mpi
    DBG 'blacs2mpi LB',LBOUND(blacs2mpi)
    DBG 'blacs2mpi UB',UBOUND(blacs2mpi)
    DBG 'row,col images',nrow_images,ncol_images
    ALLOCATE (ums%mats(nrow_images,ncol_images))
    IF (sm%matrix_type .EQ. 'S' .OR. sm%matrix_type .EQ. 'A') THEN
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
    ALLOCATE (row_dst_counts(2, 0:numproc-1, sm%nblkrows_total))
    ALLOCATE (row_dst_offsets(2, 0:numproc-1, sm%nblkrows_total))
    ithread = 0
    nthreads = 1
    ! Create and allocate the imaged matrices.
    DO row_img = 1, nrow_images
       DO col_img = 1, ncol_images
          CALL dbcsr_init (ums%mats(row_img, col_img))
          CALL dbcsr_create(ums%mats(row_img, col_img), "imaged "//sm%name,&
               target_dist,&
               dbcsr_type_normal, sm%row_blk_size, sm%col_blk_size,&
               0,0, sm%data_type, special=special, error=dbcsr_error)
          IF (sm%matrix_type .EQ. dbcsr_type_antisymmetric) THEN
             ums%mats(row_img, col_img)%m%sym_negation = .TRUE.
          ENDIF
       ENDDO
    ENDDO
    send_count(:,:,:,:) = 0
    row_dst_counts(:,:,:) = 0
    ! C
    ALLOCATE (myt_send_count(2, nrow_images, ncol_images, 0:numproc-1))
    myt_send_count(:,:,:,:) = 0
    prev_dst_p = -1
    ! Count sizes for sending.
    CALL dbcsr_iterator_start(iter, ism)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk)
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
          row_img = 1
          IF (nrow_images .GT. 1) row_img = row_img_dist (stored_row)
          col_img = 1
          IF (ncol_images .GT. 1) col_img = col_img_dist (stored_col)
          IF (predist) THEN
             IF (dbg) THEN
                WRITE(*,'(1X,A,2(1X,I5))')routineN//" Block", stored_row, stored_col
             ENDIF
             vcol = pcol*target_imgdist%col_decimation + col_img - 1
             vrow = prow*target_imgdist%row_decimation + row_img - 1
             i = vrow ; j = vcol
             IF (predistribute.eq."R") THEN
                ! shift = vpcol
                ! vrow = vrow - shift
                vrow = vrow - pcol*target_imgdist%row_multiplicity
                vrow = absmod(vrow, dbcsr_mp_nprows(mp_obj)*target_imgdist%row_decimation)
                prow = vrow / target_imgdist%row_decimation
                row_img = 1+absmod (vrow, target_imgdist%row_decimation)
                IF (dbg) &
                     WRITE(*,'(2(1X,A,2(1X,I5)))')routineN//" R shift from",&
                     i,j, "to", vrow, vcol
             ELSEIF (predistribute.eq."L") THEN
                ! shift = vprow
                ! vcol = vcol - shift
                vcol = vcol - prow*target_imgdist%col_multiplicity
                vcol = absmod(vcol, dbcsr_mp_npcols(mp_obj)*target_imgdist%col_decimation)
                pcol = vcol / target_imgdist%col_decimation
                col_img = 1+absmod (vcol, target_imgdist%col_decimation)
                IF (dbg) &
                     WRITE(*,'(2(1X,A,2(1X,I5)))')routineN//" L shift from",&
                     i,j, "to", vrow, vcol
             ENDIF
          ENDIF
          dst_p = blacs2mpi(prow, pcol)
          same_dst_p = prev_dst_p .EQ. dst_p
          prev_dst_p = dst_p
          ! These counts are meant for the thread that processes this row.
          row_dst_counts(1, dst_p, row) =&
               row_dst_counts(1, dst_p, row) + metalen
          myt_send_count(1, row_img, col_img, dst_p) =&
               myt_send_count(1, row_img, col_img, dst_p) + 1
          ! Data can be duplicated if the transpose is destined to the same
          ! process.
          IF (.NOT. same_dst_p .OR. symmetry_i .EQ. 1) THEN
             nze = row_blk_size(stored_row) * col_blk_size(stored_col)
             row_dst_counts(2, dst_p, row) =&
                  row_dst_counts(2, dst_p, row) + nze
             myt_send_count(2, row_img, col_img, dst_p) =&
                  myt_send_count(2, row_img, col_img, dst_p) + nze
          ENDIF
       ENDDO ! symmetry_i
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    ! Merge the send counts.
    send_count(:,:,:,:) = send_count(:,:,:,:) + myt_send_count(:,:,:,:)
    CALL mp_alltoall(send_count, recv_count, 2*nrow_images*ncol_images,&
         mp_group)
    ! Fill in the meta data structures and copy the data.
    DO dst_p = 0, numproc-1
       total_send_count(1, dst_p) = SUM (send_count (1, :, :, dst_p))
       total_send_count(2, dst_p) = SUM (send_count (2, :, :, dst_p))
       total_recv_count(1, dst_p) = SUM (recv_count (1, :, :, dst_p))
       total_recv_count(2, dst_p) = SUM (recv_count (2, :, :, dst_p))
    ENDDO
    ! Allocate data structures needed for data exchange.
    !ALLOCATE (recv_data(SUM(recv_count(2, :, :, :))))
    CALL dbcsr_data_init (recv_data_area)
    CALL dbcsr_data_new (recv_data_area, data_type, SUM(recv_count(2, :, :, :)))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :, :, :))))
    !ALLOCATE (send_data(SUM(send_count(2, :, :, :))))
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
    ! Calculate the displacements for rows on processors needed for
    ! the exchanges. (Used to avoid critical sections in the following
    ! DO loop that copies data into the send buffer.)
    IF(sm%nblkrows_total.GT.0) THEN
       row_dst_offsets(:,0,1) = 1
       DO dst_p = 1, numproc-1
          row_dst_offsets(1, dst_p, 1) = row_dst_offsets(1, dst_p-1, 1)&
               + metalen*total_send_count(1, dst_p-1)
          row_dst_offsets(2, dst_p, 1) = row_dst_offsets(2, dst_p-1, 1)&
               + total_send_count(2, dst_p-1)
       ENDDO
    ENDIF
    DO row = 2, sm%nblkrows_total
       row_dst_offsets(:, :, row) =&
            + row_dst_offsets(:, :, row-1)&
            + row_dst_counts(:, :, row-1)
    ENDDO
    ! Prepares some indices needed for the last DO loop that copies
    ! from buffer to local space. Placed early to take advantage of
    ! the SECTIONS.
    ALLOCATE (blk_ps(nrow_images, ncol_images))
    ALLOCATE (blks (nrow_images, ncol_images))
    blk_ps(:,:) = 1
    blks(:,:) = 1
    ! Prepares the work matrices used in the last DO loop. Placed
    ! early.
    CALL dbcsr_data_init (received_data_area)
    !CALL dbcsr_data_new (received_data_area, sm%data_type,&
    !     special_memory=.FALSE.)
    ! This is an ugly cuckoo.
    received_data_area = recv_data_area
    CALL dbcsr_data_hold(received_data_area)
    !received_data_area%d%r_dp => recv_data
    DO row_img = 1, nrow_images
       DO col_img = 1, ncol_images
          CALL dbcsr_work_create(ums%mats(row_img, col_img),&
               SUM(recv_count(1,row_img,col_img,:)), error=dbcsr_error)
          CALL dbcsr_data_hold (received_data_area)
          CALL dbcsr_data_release (ums%mats(row_img,col_img)%m%wms(1)%data_area)
          ums%mats(row_img,col_img)%m%wms(1)%data_area = received_data_area
       ENDDO
    ENDDO
    ALLOCATE (lsmp(0:numproc-1), lsdp(0:numproc-1))
    prev_dst_p = -1
    ! Copies metadata and actual data to be sent into the send buffers.
    CALL dbcsr_iterator_start(iter, ism)
    prev_row = 0
    prev_blk_p = 0
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk, blk_p=blk_p)
       ! New row
       IF (row .NE. prev_row) THEN
          lsmp = row_dst_offsets(1, :, row)
          lsdp = row_dst_offsets(2, :, row)
          prev_row = row
       ENDIF
       bp = ABS(blk_p)
       DO symmetry_i = 1, nsymmetries
          IF (symmetry_i .EQ. 1) THEN
             stored_row = row ; stored_col = col; tr = blk_p .LT. 0
          ELSE                                                
             IF (row .EQ. col) CYCLE
             stored_row = col ; stored_col = row; tr = blk_p .GT. 0
          ENDIF
          ! Where do we send this block?
          prow = row_dist(stored_row)
          pcol = col_dist(stored_col)
          row_img = 1
          IF (nrow_images .GT. 1) row_img = row_img_dist (stored_row)
          col_img = 1
          IF (ncol_images .GT. 1) col_img = col_img_dist (stored_col)
          IF (predist) THEN
             vcol = pcol*target_imgdist%col_decimation + col_img - 1
             vrow = prow*target_imgdist%row_decimation + row_img - 1
             IF (predistribute.eq."R") THEN
                ! shift = vpcol
                ! vrow = vrow + shift
                vrow = vrow - pcol*target_imgdist%row_multiplicity
                vrow = absmod(vrow, dbcsr_mp_nprows(mp_obj)*target_imgdist%row_decimation)
                prow = vrow / target_imgdist%row_decimation
                row_img = 1+absmod (vrow, target_imgdist%row_decimation)
             ELSEIF (predistribute.eq."L") THEN
                ! shift = vprow
                ! vcol = vcol + shift
                vcol = vcol - prow*target_imgdist%col_multiplicity
                vcol = absmod(vcol, dbcsr_mp_npcols(mp_obj)*target_imgdist%col_decimation)
                pcol = vcol / target_imgdist%col_decimation
                col_img = 1+absmod (vcol, target_imgdist%col_decimation)
             ENDIF
          ENDIF
          dst_p = blacs2mpi(prow, pcol)
          same_dst_p = dst_p .EQ. prev_dst_p
          prev_dst_p = dst_p
          sm_pos = lsmp(dst_p)
          lsmp(dst_p) = lsmp(dst_p) + metalen
          send_meta(sm_pos) = stored_row
          IF (.NOT. same_dst_p .OR. symmetry_i .EQ. 1) THEN
             nze = row_blk_size(stored_row) * col_blk_size(stored_col)
             sd_pos = lsdp(dst_p)
             lsdp(dst_p) = lsdp(dst_p) + nze
             CALL dbcsr_data_set (send_data_area, sd_pos, nze,&
                  sm%data_area, bp)
             !send_data(sd_pos:sd_pos+nze-1) = r_dp(bp:bp+nze-1)
             send_meta(sm_pos+1) = stored_col
             send_meta(sm_pos+2) = SGN (&
                  sd_pos-sd_disp(dst_p)+1, blk_p, tr)
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
    ! Exchange the data and metadata structures.
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       CALL hybrid_alltoall_s1(&
            send_data_area%d%r_sp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data_area%d%r_sp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            no_hybrid = (nsymmetries .GT. 1),&
            mp_env = mp_obj)
    CASE (dbcsr_type_real_8)
       !CALL mp_alltoall(&
       !     send_data_area%d%r_dp(:), total_send_count(2,:), sd_disp(:)-1,&
       !     recv_data_area%d%r_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
       !     mp_group)
       CALL hybrid_alltoall_d1 (&
            send_data_area%d%r_dp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data_area%d%r_dp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            no_hybrid = (nsymmetries .GT. 1),&
            mp_env = mp_obj)
    CASE (dbcsr_type_complex_4)
       CALL hybrid_alltoall_c1(&
            send_data_area%d%c_sp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data_area%d%c_sp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            no_hybrid = (nsymmetries .GT. 1),&
            mp_env = mp_obj)
    CASE (dbcsr_type_complex_8)
       CALL hybrid_alltoall_z1(&
            send_data_area%d%c_dp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data_area%d%c_dp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            no_hybrid = (nsymmetries .GT. 1),&
            mp_env = mp_obj)
    END SELECT
    !CALL mp_alltoall(send_data(:), total_send_count(2,:), sd_disp(:)-1,&
    !     recv_data(:), total_recv_count(2,:), rd_disp(:)-1, mp_group)
    CALL hybrid_alltoall_i1(&
         send_meta(:), SUM(metalen*total_send_count(1,:)), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), SUM(metalen*total_recv_count(1,:)), metalen*total_recv_count(1,:), rm_disp(:)-1,&
         mp_env = mp_obj,&
         no_hybrid = nsymmetries .GT. 1)
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
          blk_ps(row_img,col_img) = blk_ps(row_img,col_img) + nze
          IF (stored_col .GT. 0) THEN
             blk_p = data_p
             data_p = data_p + nze
          ELSE
             blk_p = prev_blk_p
          ENDIF
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
          IF (nrow_images.EQ.1 .AND. ncol_images.eq.1 .OR. nocopy) THEN
             CALL dbcsr_finalize(ums%mats(row_img,col_img),reshuffle=.FALSE.,error=dbcsr_error)
          ELSE
             CALL dbcsr_finalize(ums%mats(row_img,col_img),reshuffle=.TRUE.,error=dbcsr_error)
          ENDIF
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
    CALL dbcsr_data_release (received_data_area)
    !
    CALL dbcsr_data_release (recv_data_area)
    !
    tstop = m_walltime ()
    !WRITE(*,'(1X,A,1X,EN12.4)')'make_column_images time:',tstop-tstart
    DBG 'ums(1) name',ums%mats(1,1)%m%name
    DBG 'Done making images'
    CALL dbcsr_error_stop(error_handler, dbcsr_error)
  END SUBROUTINE make_images




! *****************************************************************************
!> \brief Makes dense matrices for the image matrices.
!> \param[in,out] images          matrix images
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE dbcsr_make_images_dense (images, predistributed, error)
    TYPE(dbcsr_2d_array_type), INTENT(INOUT) :: images
    CHARACTER, INTENT(IN), OPTIONAL          :: predistributed
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_images_dense', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: predist
    INTEGER :: error_handler, i, mat_col, mat_row, mypcol, myprow, &
      ncol_images, nfullcols, nfullrows, nlocal_blocks, nlvc, nlvr, npcols, &
      nprows, nrow_images, nvpcols, nvprows, which_icol, which_irow, &
      which_pcol, which_prow, which_vcol, which_vrow
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: local_col_blk_offsets, &
                                                local_row_blk_offsets, &
                                                local_virt_cols, &
                                                local_virt_rows
    INTEGER, DIMENSION(:), POINTER           :: col_blk_offsets, meta, &
                                                old_col_images, &
                                                old_row_images, &
                                                row_blk_offsets
    TYPE(array_i1d_obj)                      :: new_col_dist, new_row_dist, &
                                                old_cbo, old_rbo
    TYPE(dbcsr_data_obj)                     :: dense_data, old_data
    TYPE(dbcsr_distribution_obj)             :: new_distribution, &
                                                old_distribution
    TYPE(dbcsr_mp_obj)                       :: mp

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    predist = '0'
    IF (PRESENT (predistributed)) predist = predistributed
    old_distribution = images%image_dist%main
    mp = dbcsr_distribution_mp (old_distribution)
    nprows = dbcsr_mp_nprows (mp)
    npcols = dbcsr_mp_npcols (mp)
    myprow = dbcsr_mp_myprow (mp)
    mypcol = dbcsr_mp_mypcol (mp)
    !
    nrow_images = images%image_dist%row_decimation
    ncol_images = images%image_dist%col_decimation
    nvprows = nprows * nrow_images
    nvpcols = npcols * ncol_images
    !
    ! Max allocation
    ALLOCATE (local_virt_rows (dbcsr_distribution_nrows (old_distribution)))
    ALLOCATE (local_virt_cols (dbcsr_distribution_ncols (old_distribution)))
    ALLOCATE (local_row_blk_offsets(MAX (&
         dbcsr_distribution_nrows(old_distribution)+1, 1)))
    ALLOCATE (local_col_blk_offsets(MAX (&
         dbcsr_distribution_ncols(old_distribution)+1, 1)))
    old_row_images => array_data (images%image_dist%row_image)
    old_col_images => array_data (images%image_dist%col_image)
    old_rbo = images%mats(1,1)%m%row_blk_offset
    old_cbo = images%mats(1,1)%m%col_blk_offset
    CALL array_hold (old_rbo)
    CALL array_hold (old_cbo)
    row_blk_offsets => array_data (old_rbo)
    col_blk_offsets => array_data (old_cbo)
    !
    ! New distribution
    CALL array_nullify (new_row_dist)
    CALL array_nullify (new_col_dist)
    CALL array_new (new_row_dist, &
         (/ (INT(i/nrow_images), i=0, nvprows-1) /), 1)
    CALL array_new (new_col_dist,&
         (/ (INT(i/ncol_images), i=0, nvpcols-1) /), 1)
    CALL dbcsr_distribution_init (new_distribution)
    CALL dbcsr_distribution_new (new_distribution, mp,&
         new_row_dist, new_col_dist)
    CALL array_release (new_row_dist)
    CALL array_release (new_col_dist)
    !
    DO mat_row = 1, nrow_images
       DO mat_col = 1, ncol_images
          which_vrow = myprow*nrow_images+mat_row-1
          which_vcol = mypcol*ncol_images+mat_col-1
          SELECT CASE (predist)
          CASE ('R')
             which_vrow = MODULO (&
                  which_vrow + mypcol*images%image_dist%row_multiplicity,&
                  nprows*nrow_images)
          CASE ('L')
             which_vcol = MODULO (&
                  which_vcol + myprow*images%image_dist%col_multiplicity,&
                  npcols*ncol_images)
          END SELECT
          which_prow = which_vrow / nrow_images
          which_irow = MODULO (which_vrow, nrow_images)+1
          which_pcol = which_vcol / ncol_images
          which_icol = MODULO (which_vcol, ncol_images)+1
          CALL find_local_virtual_elements (local_virt_rows, nlvr,&
               array_data(dbcsr_distribution_row_dist (&
               old_distribution)),&
               old_row_images, which_prow, which_irow)
          CALL find_local_virtual_elements (local_virt_cols, nlvc,&
               array_data(dbcsr_distribution_col_dist (&
               old_distribution)),&
               old_col_images, which_pcol, which_icol)
          ! Conversion of global to local offsets (old blocks)
          CALL global_offsets_to_local (&
               row_blk_offsets,&
               local_virt_rows(1:nlvr),&
               local_row_blk_offsets)
          CALL global_offsets_to_local (&
               col_blk_offsets,&
               local_virt_cols(1:nlvc),&
               local_col_blk_offsets)
          ! New data area
          old_data = images%mats(mat_row, mat_col)%m%data_area
          CALL dbcsr_data_init (dense_data)
          nfullrows = nfull_elements (row_blk_offsets, local_virt_rows(1:nlvr))
          nfullcols = nfull_elements (col_blk_offsets, local_virt_cols(1:nlvc))
          CALL dbcsr_data_new (dense_data, dbcsr_get_data_type (old_data),&
               data_size = nfullrows * nfullcols,&
               special_memory=dbcsr_uses_special_memory (old_data))
          ! Reshuffle the data
          CALL make_dense_data (images%mats(mat_row,mat_col),&
               dense_data, nfullrows, nfullcols,&
               local_row_blk_offsets, local_col_blk_offsets)
          CALL dbcsr_data_release (old_data)
          images%mats(mat_row, mat_col)%m%data_area = dense_data
          ! Switch disttributions
          CALL dbcsr_distribution_release (images%mats(mat_row, mat_col)%m%dist)
          CALL dbcsr_distribution_hold (new_distribution)
          images%mats(mat_row, mat_col)%m%dist = new_distribution
          ! New sizes
          CALL match_sizes_to_dist (new_distribution,&
               old_distribution,&
               images%mats(mat_row, mat_col)%m%row_blk_size,&
               images%mats(mat_row, mat_col)%m%col_blk_size,&
               images%mats(mat_row, mat_col)%m%row_blk_offset,&
               images%mats(mat_row, mat_col)%m%col_blk_offset,&
               old_image_dist=images%image_dist)
          nlocal_blocks = MIN (nfullrows*nfullcols, 1)
          CALL dbcsr_assert (nlocal_blocks, "LE", 1,&
               dbcsr_warning_level,&
               dbcsr_internal_error, routineN,&
               "Too many local dense blocks.",__LINE__,error)
          ! New index structure
          CALL dbcsr_addto_index_array (images%mats(mat_row,mat_col)%m,&
               dbcsr_slot_row_p,&
               reservation=nvprows+1, error=error)
          CALL dbcsr_addto_index_array (images%mats(mat_row,mat_col)%m,&
               dbcsr_slot_col_i,&
               reservation=nlocal_blocks, error=error)
          CALL dbcsr_addto_index_array (images%mats(mat_row,mat_col)%m,&
               dbcsr_slot_blk_p,&
               reservation=nlocal_blocks, error=error)
          meta => images%mats(mat_row,mat_col)%m%index(1:dbcsr_meta_size)
          CALL dbcsr_pack_meta (images%mats(mat_row,mat_col)%m, meta)
          meta(dbcsr_slot_nze) = nfullrows*nfullcols
          meta(dbcsr_slot_nblks) = nlocal_blocks
          CALL make_dense_index (images%mats(mat_row,mat_col)%m%row_p,&
               images%mats(mat_row,mat_col)%m%col_i,&
               images%mats(mat_row,mat_col)%m%blk_p,&
               nvprows, nvpcols,&
               which_vrow,&
               which_vcol, meta, nlocal_blocks)
          CALL dbcsr_unpack_meta (images%mats(mat_row, mat_col)%m, meta)
       ENDDO
    ENDDO
    CALL array_release (old_rbo)
    CALL array_release (old_cbo)
    CALL dbcsr_distribution_release (images%image_dist%main)
    images%image_dist%main = new_distribution
    CALL make_image_distribution_dense (images%image_dist)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_make_images_dense


! *****************************************************************************
!> \brief Makes a dense matrix, inplace.
!> \param[in,out] matrix      matrix to make dense
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE dbcsr_make_dense (matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_dense', &
      routineP = moduleN//':'//routineN

    INTEGER :: error_handler, i, mypcol, myprow, nfullcols_local, &
      nfullrows_local, nlocal_blocks, npcols, nprows
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: local_col_blk_offsets, &
                                                local_row_blk_offsets
    INTEGER, DIMENSION(:), POINTER           :: meta
    TYPE(array_i1d_obj)                      :: new_col_dist, new_row_dist
    TYPE(dbcsr_data_obj)                     :: dense_data, old_data
    TYPE(dbcsr_distribution_obj)             :: distribution, old_distribution
    TYPE(dbcsr_mp_obj)                       :: mp

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    distribution = dbcsr_distribution (matrix)
    old_distribution = distribution
    mp = dbcsr_distribution_mp (distribution)
    nprows = dbcsr_mp_nprows (mp)
    npcols = dbcsr_mp_npcols (mp)
    myprow = dbcsr_mp_myprow (mp)
    mypcol = dbcsr_mp_mypcol (mp)
    nfullrows_local = matrix%m%nfullrows_local
    nfullcols_local = matrix%m%nfullcols_local
    ! Conversion of global to local offsets
    ALLOCATE (local_row_blk_offsets(matrix%m%nblkrows_total+1))
    ALLOCATE (local_col_blk_offsets(matrix%m%nblkcols_total+1))
    CALL global_offsets_to_local (array_data(matrix%m%row_blk_offset),&
         array_data(dbcsr_distribution_local_rows (distribution)),&
         local_row_blk_offsets)
    CALL global_offsets_to_local (array_data(matrix%m%col_blk_offset),&
         array_data(dbcsr_distribution_local_cols (distribution)),&
         local_col_blk_offsets)
    ! New data area
    old_data = matrix%m%data_area
    CALL dbcsr_data_init (dense_data)
    CALL dbcsr_data_new (dense_data, dbcsr_get_data_type (old_data),&
         data_size = nfullrows_local*nfullcols_local,&
         special_memory=dbcsr_uses_special_memory (old_data))
    ! Reshuffle the data
    CALL make_dense_data (matrix, dense_data, nfullrows_local, nfullcols_local,&
         local_row_blk_offsets, local_col_blk_offsets)
    CALL dbcsr_data_release (old_data)
    matrix%m%data_area = dense_data
    ! New distribution
    CALL array_nullify (new_row_dist)
    CALL array_nullify (new_col_dist)
    CALL array_new (new_row_dist, (/ (i, i=0, nprows-1) /), 1)
    CALL array_new (new_col_dist, (/ (i, i=0, npcols-1) /), 1)
    CALL dbcsr_distribution_new (distribution, mp,&
         new_row_dist, new_col_dist)
    CALL array_release (new_row_dist)
    CALL array_release (new_col_dist)
    nlocal_blocks = dbcsr_distribution_nlocal_rows(distribution)
    CALL dbcsr_assert (nlocal_blocks, "LE", 1, dbcsr_warning_level,&
         dbcsr_internal_error, routineN,&
         "There should be 0 or 1 local dense blocks.",__LINE__,error)
    matrix%m%dist = distribution
    ! New sizes
    CALL match_sizes_to_dist (matrix%m%dist, old_distribution,&
         matrix%m%row_blk_size, matrix%m%col_blk_size,&
         matrix%m%row_blk_offset, matrix%m%col_blk_offset)
    CALL dbcsr_distribution_release (old_distribution)
    ! New index structure
    CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_row_p,&
         reservation=nprows+1, error=error)
    CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_col_i,&
         reservation=nlocal_blocks, error=error)
    CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_blk_p,&
         reservation=nlocal_blocks, error=error)
    meta => matrix%m%index(1:dbcsr_meta_size)
    CALL dbcsr_pack_meta (matrix%m, meta)
    meta(dbcsr_slot_nze) = nfullrows_local*nfullcols_local
    meta(dbcsr_slot_nblks) = nlocal_blocks
    CALL make_dense_index (matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p,&
         nprows, npcols, myprow, mypcol, meta, nlocal_blocks)
    CALL dbcsr_unpack_meta (matrix%m, meta)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_make_dense


! *****************************************************************************
!> \brief Makes a blocked matrix from a dense matrix, inplace
!> \param[in,out] matrix      matrix to make blocked
!> \param[in] distribution    distribution of blocked rows and columns
!> \param[in] row_blk_offsets, col_blk_offsets   row and column block offsets
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE dbcsr_make_undense (matrix, distribution,&
       row_blk_offsets, col_blk_offsets, row_blk_sizes, col_blk_sizes, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: distribution
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_offsets, &
                                                col_blk_offsets, &
                                                row_blk_sizes, col_blk_sizes
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_make_undense', &
      routineP = moduleN//':'//routineN

    INTEGER :: error_handler, nblkcols_local, nblkcols_total, nblkrows_local, &
      nblkrows_total, nfullcols_local, nfullrows_local
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: local_col_offsets, &
                                                local_row_offsets
    INTEGER, DIMENSION(:), POINTER           :: local_cols, local_rows, meta
    TYPE(dbcsr_data_obj)                     :: blocked_data, dense_data
    TYPE(dbcsr_distribution_obj)             :: dense_distribution

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    dense_distribution = dbcsr_distribution (matrix)
    nfullrows_local = matrix%m%nfullrows_local
    nfullcols_local = matrix%m%nfullcols_local
    nblkrows_local = dbcsr_distribution_nlocal_rows (distribution)
    nblkcols_local = dbcsr_distribution_nlocal_cols (distribution)
    nblkrows_total = dbcsr_distribution_nrows (distribution)
    nblkcols_total = dbcsr_distribution_ncols (distribution)
    local_rows => array_data (dbcsr_distribution_local_rows (distribution))
    local_cols => array_data (dbcsr_distribution_local_cols (distribution))
    !
    CALL array_release (matrix%m%row_blk_offset)
    CALL array_release (matrix%m%col_blk_offset)
    CALL array_release (matrix%m%row_blk_size)
    CALL array_release (matrix%m%col_blk_size)
    matrix%m%row_blk_offset = row_blk_offsets
    matrix%m%col_blk_offset = col_blk_offsets
    matrix%m%row_blk_size = row_blk_sizes
    matrix%m%col_blk_size = col_blk_sizes
    CALL array_hold (matrix%m%row_blk_offset)
    CALL array_hold (matrix%m%col_blk_offset)
    CALL array_hold (matrix%m%row_blk_size)
    CALL array_hold (matrix%m%col_blk_size)
    ! New distribution
    matrix%m%dist = distribution
    CALL dbcsr_distribution_hold (matrix%m%dist)
    CALL dbcsr_distribution_release (dense_distribution)
    !
    ALLOCATE (local_row_offsets(nblkrows_total+1))
    ALLOCATE (local_col_offsets(nblkcols_total+1))
    CALL dbcsr_clearfrom_index_array (matrix%m, dbcsr_slot_row_p)
    CALL dbcsr_clearfrom_index_array (matrix%m, dbcsr_slot_col_i)
    CALL dbcsr_clearfrom_index_array (matrix%m, dbcsr_slot_blk_p)
    CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_row_p,&
         reservation=nblkrows_total+1, error=error)
    CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_col_i,&
         reservation=nblkrows_local*nblkcols_local, error=error)
    CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_blk_p,&
         reservation=nblkrows_local*nblkcols_local, error=error)
    meta => matrix%m%index(1:dbcsr_meta_size)
    CALL dbcsr_pack_meta (matrix%m, meta)
    meta(dbcsr_slot_nblks) = nblkrows_local*nblkcols_local
    meta(dbcsr_slot_nze) = nfullrows_local * nfullcols_local
    CALL global_offsets_to_local (array_data (row_blk_offsets),&
         local_rows, local_row_offsets(1:nblkrows_local+1))
    CALL global_offsets_to_local (array_data (col_blk_offsets),&
         local_cols, local_col_offsets(1:nblkcols_local+1))
    CALL make_undense_index (matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p,&
         distribution,&
         local_row_offsets(1:nblkrows_local+1),&
         local_col_offsets(1:nblkcols_local+1),&
         meta)
    CALL dbcsr_unpack_meta (matrix%m, meta)
    !
    ! New data area
    CALL global_offsets_to_local (array_data (row_blk_offsets),&
         local_rows, local_row_offsets)
    CALL global_offsets_to_local (array_data (col_blk_offsets),&
         local_cols, local_col_offsets)
    dense_data = matrix%m%data_area
    CALL dbcsr_data_init (blocked_data)
    CALL dbcsr_data_new (blocked_data, dbcsr_get_data_type (dense_data),&
         data_size = nfullrows_local*nfullcols_local,&
         special_memory=dbcsr_uses_special_memory (dense_data))
    matrix%m%data_area = blocked_data
    ! Reshuffle the data
    CALL make_undense_data (matrix, dense_data,&
         nfullrows_local, nfullcols_local,&
         local_row_offsets, local_col_offsets)
    CALL dbcsr_data_release (dense_data)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_make_undense

! *****************************************************************************
!> \brief Shuffles the data from blocked to standard dense form
!> \param[in] matrix                 Existing blocked matrix
!> \param[in,out] dense_data         Dense data
!> \param[in] nfullrows, nfullcols    size of new data
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE make_dense_data (matrix, dense_data, nfullrows, nfullcols,&
       row_blk_offsets, col_blk_offsets)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: dense_data
    INTEGER, INTENT(IN)                      :: nfullrows, nfullcols
    INTEGER, DIMENSION(:), INTENT(IN)        :: row_blk_offsets, &
                                                col_blk_offsets

    CHARACTER(len=*), PARAMETER :: routineN = 'make_dense_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_col, blk_row, col_offset, &
                                                col_size, row_offset, row_size
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: block, dense2d
    TYPE(dbcsr_error_type)                   :: error
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (dbcsr_get_data_size (dense_data), 'GE', nfullrows*nfullcols,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Dense data too small",__LINE__,error)
    CALL dbcsr_data_clear (dense_data)
    IF (dbcsr_get_data_size(matrix%m%data_area) .GT. 0&
         .AND. nfullrows .GT. 0 .AND. nfullcols .GT. 0) THEN
       CALL dbcsr_data_init (block)
       CALL dbcsr_data_new (block,&
         dbcsr_type_1d_to_2d(dbcsr_get_data_type(dense_data)))
       CALL dbcsr_data_init (dense2d)
       CALL dbcsr_data_new (dense2d, dbcsr_type_1d_to_2d(dbcsr_get_data_type(dense_data)))
       CALL dbcsr_data_set_pointer (dense2d, rsize=nfullrows, csize=nfullcols,&
            pointee=dense_data)
       CALL dbcsr_iterator_start (iter, matrix)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, blk_row, blk_col, block, tr,&
               row_size=row_size, col_size=col_size)
          row_offset = row_blk_offsets(blk_row)
          col_offset = col_blk_offsets(blk_col)
          IF (tr) THEN
             CALL dbcsr_block_transpose (&
                  dst=dense2d, lb=row_offset, lb2=col_offset,&
                  row_size=row_size, col_size=col_size,&
                  src=block, source_lb=1, source_lb2=1)
          ELSE
             CALL dbcsr_data_set (&
                  dst=dense2d, lb=row_offset, lb2=col_offset,&
                  data_size=row_size, data_size2=col_size,&
                  src=block, source_lb=1, source_lb2=1)
          ENDIF
       ENDDO
       CALL dbcsr_iterator_stop (iter)
       CALL dbcsr_data_clear_pointer (block)
       CALL dbcsr_data_release (block)
       CALL dbcsr_data_clear_pointer (dense2d)
       CALL dbcsr_data_release (dense2d)
    ENDIF
  END SUBROUTINE make_dense_data

! *****************************************************************************
!> \brief Shuffles the data from standard dense to blocked form
!> \param[in,out] matrix             Matrix with data to fill
!> \param[in] dense_data             Dense data
!> \param[in] nfullrows, nfullcols   size of dense data
!> \note Used for making matrices dense/undense
! *****************************************************************************
  SUBROUTINE make_undense_data (matrix, dense_data, nfullrows, nfullcols,&
       row_blk_offsets, col_blk_offsets)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_data_obj), INTENT(IN)         :: dense_data
    INTEGER, INTENT(IN)                      :: nfullrows, nfullcols
    INTEGER, DIMENSION(:), INTENT(IN)        :: row_blk_offsets, &
                                                col_blk_offsets

    CHARACTER(len=*), PARAMETER :: routineN = 'make_undense_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_col, blk_row, col_offset, &
                                                col_size, row_offset, row_size
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: block, dense2d
    TYPE(dbcsr_error_type)                   :: error
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (dbcsr_get_data_size (dense_data), 'GE', nfullrows*nfullcols,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Dense data too small",__LINE__,error)
    IF (dbcsr_get_data_size(matrix%m%data_area) .GT. 0) THEN
       CALL dbcsr_data_clear (matrix%m%data_area)
       CALL dbcsr_data_init (block)
       CALL dbcsr_data_new (block,&
            dbcsr_type_1d_to_2d(dbcsr_get_data_type(dense_data)))
       CALL dbcsr_data_init (dense2d)
       CALL dbcsr_data_new (dense2d,&
            dbcsr_type_1d_to_2d (dbcsr_get_data_type(dense_data)))
       CALL dbcsr_data_set_pointer (dense2d, rsize=nfullrows, csize=nfullcols,&
            pointee=dense_data)
       CALL dbcsr_iterator_start (iter, matrix)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, blk_row, blk_col, block, tr,&
               row_size=row_size, col_size=col_size)
          row_offset = row_blk_offsets(blk_row)
          col_offset = col_blk_offsets(blk_col)
          IF (tr) THEN
             CALL dbcsr_block_transpose (&
                  dst=block, lb=1, lb2=1,&
                  row_size=row_size, col_size=col_size,&
                  src=dense2d, source_lb=row_offset, source_lb2=col_offset)
          ELSE
             CALL dbcsr_data_set (&
                  dst=block, lb=1, lb2=1,&
                  data_size=row_size, data_size2=col_size,&
                  src=dense2d, source_lb=row_offset, source_lb2=col_offset)
          ENDIF
       ENDDO
       CALL dbcsr_iterator_stop (iter)
       CALL dbcsr_data_clear_pointer (block)
       CALL dbcsr_data_release (block)
       CALL dbcsr_data_clear_pointer (dense2d)
       CALL dbcsr_data_release (dense2d)
    ENDIF
  END SUBROUTINE make_undense_data





! *****************************************************************************
!> \brief Replicates a DBCSR on all processors.
!> \param[in,out] matrix      matrix to replicate
! *****************************************************************************
  SUBROUTINE dbcsr_replicate_all(matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_replicate_all', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: debug_level = 0, metalen = 2

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: c_dp
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: c_sp
    INTEGER :: blk, blk_l, blk_ps, blks, col, dst_p, error_handler, mp_group, &
      numnodes, nze, row, src_p, stored_col, stored_row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rd_disp, recv_meta, rm_disp, &
                                                sd_disp, sdp, send_meta, &
                                                sm_disp, smp
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: recv_count, send_count, &
                                                total_recv_count, &
                                                total_send_count
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, col_dist, &
                                                row_blk_size, row_dist, &
                                                tmp_index
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    LOGICAL                                  :: tr
    REAL(KIND=dp)                            :: tstart, tstop
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: send_data
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: recv_data
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_dp
    REAL(KIND=sp), DIMENSION(:), POINTER     :: r_sp
    TYPE(dbcsr_data_obj)                     :: tmp_data
    TYPE(dbcsr_distribution_obj)             :: target_dist
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_obj)                          :: replicated

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    tstart = m_walltime ()
    CALL dbcsr_assert (dbcsr_valid_index(matrix%m),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Matrix not initialized.",__LINE__,error)
    CALL dbcsr_assert (matrix%m%replication_type .EQ. dbcsr_repl_none,&
         dbcsr_warning_level, dbcsr_caller_error, routineN,&
         "Replicating a non-distributed matrix makes no sense.",__LINE__,error)
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
    target_dist = matrix%m%dist
    row_dist => array_data (dbcsr_distribution_row_dist (target_dist))
    col_dist => array_data (dbcsr_distribution_col_dist (target_dist))
    mp_obj = dbcsr_distribution_mp (target_dist)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    numnodes = dbcsr_mp_numnodes (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)
    DBGV 'row_dist', row_dist
    DBGV 'col_dist', col_dist
    DBG 'row, col maxval:', MAXVAL(row_dist), MAXVAL(col_dist)
    DBG 'blacs2mpi',blacs2mpi
    DBG 'blacs2mpi LB',LBOUND(blacs2mpi)
    DBG 'blacs2mpi UB',UBOUND(blacs2mpi)
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
    CALL dbcsr_init (replicated)
    CALL dbcsr_create(replicated, 'Replicated '//matrix%m%name, matrix%m%dist,&
         dbcsr_type_normal, matrix%m%row_blk_size, matrix%m%col_blk_size,&
         0, 0, matrix%m%data_type, special=matrix%m%special_memory,&
         replication_type=dbcsr_repl_full, error=error)
    ! Count initial sizes for sending.
    DBG 'A'
    send_count(:,:) = 0
    CALL dbcsr_iterator_start (iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, r_dp, tr, blk)
       nze = SIZE(r_dp)
       DO dst_p = 0, numnodes-1
          send_count(1, dst_p) = send_count(1, dst_p) + 1
          send_count(2, dst_p) = send_count(2, dst_p) + nze
       END DO
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    DBG 'B'
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    DBG 'C'
    DBG 'send counts',send_count
    DBG 'recv counts',recv_count
    ! Allocate data structures needed for data exchange.
    ALLOCATE (recv_data(SUM(recv_count(2, :))))
    ALLOCATE (recv_meta(metalen*SUM(recv_count(1, :))))
    ALLOCATE (send_data(SUM(send_count(2, :))))
    ALLOCATE (send_meta(metalen*SUM(send_count(1, :))))
    DBG 'send data len',SIZE (send_data)
    DBG 'send meta len',SIZE (send_meta)
    DBG 'recv data len',SIZE (recv_data)
    DBG 'recv meta len',SIZE (recv_meta)
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
    !DO dst_p = 0, numnodes-1
    !      sm_disp (dst_p) = sm_disp (dst_p) + metalen*send_count(1, dst_p)
    !      sd_disp (dst_p) = sd_disp (dst_p) + send_count(2, dst_p)
    !      rm_disp (dst_p) = rm_disp (dst_p) + metalen*recv_count(1, dst_p)
    !      rd_disp (dst_p) = rd_disp (dst_p) + recv_count(2, dst_p)
    !ENDDO
    sdp(:) = sd_disp
    smp(:) = sm_disp
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    DBG 'D'
    CALL dbcsr_iterator_start (iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, r_dp, tr, blk)
       nze = SIZE(r_dp)
       DO dst_p = 0, numnodes-1
          send_meta(smp(dst_p)+1) = col
          IF (.NOT. tr) THEN
             send_meta(smp(dst_p)) = row
             send_data(sdp(dst_p):sdp(dst_p)+nze-1) =&
                  r_dp(:)
          ELSE
             send_meta(smp(dst_p)) = -row ! Mark block transposed
             send_data(sdp(dst_p):sdp(dst_p)+nze-1) =&
                  r_dp(:)
          ENDIF
          smp(dst_p) = smp(dst_p) + metalen
          sdp(dst_p) = sdp(dst_p) + nze
       END DO
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    ! Exchange the data and metadata structures.
    DBG 'E'
    CALL mp_alltoall(send_data(:), total_send_count(2,:), sd_disp(:)-1,&
         recv_data(:), total_recv_count(2,:), rd_disp(:)-1, mp_group)
    DBG 'F'
    DBG 'send_meta',send_meta
    CALL mp_alltoall(send_meta(:), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), metalen*total_recv_count(1,:), rm_disp(:)-1, mp_group)
    DBG 'recv_meta',recv_meta
    DBG 'G'
    ! Now fill in the data.
    CALL dbcsr_work_create(replicated,&
            SUM(recv_count(1,:)),&
            SUM(recv_count(2,:)), n=1, error=error)
    !DBG 'send disps data',sd_disp
    !DBG 'send disps meta',sm_disp
    !DBG 'recv disps data',rd_disp
    !DBG 'recv disps meta',rm_disp
    !blk_p = 1
    !blk = 1
    blk_ps = 1
    blks = 1
    DO src_p = 0, numnodes-1
       nze = recv_count(2, src_p)
       CALL dbcsr_data_set (replicated%m%wms(1)%data_area, blk_ps, nze,&
            recv_data, rd_disp(src_p))
       DO blk_l = 1, recv_count(1, src_p)
          stored_row = recv_meta(rm_disp(src_p)+metalen*(blk_l-1))
          stored_col = recv_meta(rm_disp(src_p)+metalen*(blk_l-1)+1)
          DBG 'blk: p',blks,src_p,'; row, col',&
               stored_row,stored_col,&
               '; blk_ps',blk_ps
          replicated%m%wms(1)%row_i(blks) = ABS(stored_row)
          replicated%m%wms(1)%col_i(blks) = stored_col
          replicated%m%wms(1)%blk_p(blks) = SIGN(blk_ps, stored_row)
          nze = row_blk_size(ABS(stored_row))&
               * col_blk_size(stored_col)
          blk_ps = blk_ps + nze
          blks = blks + 1
       ENDDO
    ENDDO
    !
    replicated%m%wms(1)%lastblk = blks - 1
    replicated%m%wms(1)%datasize = blk_ps - 1
    DBG 'Finalizing normalization'
    CALL dbcsr_finalize(replicated, resort=.TRUE., reshuffle=.TRUE., error=error)
    tstop = m_walltime ()
    DBG 'time:',tstop-tstart
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_data)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_data)
    DEALLOCATE(send_meta)
    matrix%m%replication_type = dbcsr_repl_full
    ! Now replace the data and index
    tmp_data = matrix%m%data_area
    matrix%m%data_area = replicated%m%data_area
    replicated%m%data_area = tmp_data
    tmp_index => matrix%m%index
    matrix%m%index => replicated%m%index
    replicated%m%index => tmp_index
    CALL dbcsr_repoint_index (matrix%m)
    matrix%m%nze = replicated%m%nze
    matrix%m%nblks = replicated%m%nblks
    CALL dbcsr_release (replicated)
    DBG 'Done replicating'
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_replicate_all


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

    CHARACTER                                :: data_type
    INTEGER :: blk, blk_col_new, blk_col_prev_new, blk_ps, blk_row_new, blks, &
      cnt_fnd, cnt_new, cnt_skip, col, col_int, col_offset_new, &
      col_offset_old, col_rle, col_size, data_offset_l, dst_p, error_handler, &
      i, meta_l, mp_group, nblkcols_total_new, nblkrows_total_new, numnodes, &
      nze_rle, row, row_int, row_offset_new, row_offset_old, row_prev, &
      row_rle, row_size, src_p, stored_col_new, stored_row_new
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
    LOGICAL                                  :: found, my_keep_sparsity, &
                                                new_meta, sym, tr, valid_block
    REAL(kind=dp)                            :: t_all, t_blk, t_blk2, &
                                                t_count, t_fill, t_prolog, &
                                                t_unpack, t_xfer1, t_xfer2
    TYPE(dbcsr_data_obj)                     :: buff_data, data_block, &
                                                recv_data, send_data, src_2d
    TYPE(dbcsr_distribution_obj)             :: dist_new, dist_old
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj_new

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    t_all = m_walltime()
    t_prolog = t_all
    CALL dbcsr_assert (dbcsr_valid_index(matrix%m),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Input not valid.",__LINE__,error)
    CALL dbcsr_assert (dbcsr_is_initialized(redist%m),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Output matrix not initialized.",__LINE__,error)
    CALL dbcsr_assert (matrix%m%replication_type .EQ. dbcsr_repl_none,&
         dbcsr_warning_level, dbcsr_caller_error, routineN,&
         "Can not redistribute replicated matrix.",__LINE__,error)
    !
    my_keep_sparsity = .FALSE.
    IF (PRESENT (keep_sparsity)) my_keep_sparsity = keep_sparsity
    !
    sym = dbcsr_get_matrix_type(redist) .EQ. dbcsr_type_symmetric .OR.&
         dbcsr_get_matrix_type(redist) .EQ. dbcsr_type_antisymmetric
    data_type = matrix%m%data_type
    ! Get row and column start and end positions
    ! Old matrix
    dist_old = dbcsr_distribution (matrix)
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
    nblkrows_total_new = dbcsr_nblkrows_total (redist)
    nblkcols_total_new = dbcsr_nblkcols_total (redist)
    ! Create mappings
    i = dbcsr_nfullrows_total (redist)
    ALLOCATE (row_reblocks (2,i))
    ALLOCATE (n_row_reblocks (2, dbcsr_nblkrows_total (matrix)))
    CALL dbcsr_reblocking_targets (row_reblocks, i, n_row_reblocks,&
         row_blk_size_old, row_blk_size_new)
    i = dbcsr_nfullcols_total (redist)
    ALLOCATE (col_reblocks (2,i))
    ALLOCATE (n_col_reblocks (2, dbcsr_nblkcols_total (matrix)))
    CALL dbcsr_reblocking_targets (col_reblocks, i, n_col_reblocks,&
         col_blk_size_old, col_blk_size_new)
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
    IF (dbg) CALL dbcsr_print(matrix,error=error)
    !cs1 = dbcsr_checksum (matrix)
    !cs1 = dbcsr_checksum (matrix)
    !call dbcsr_print(matrix)
    !
    t_count = m_walltime()
    t_prolog = t_count - t_prolog
    !
    ! Count initial sizes for sending.
    !
    ! We go through every element of every local block and determine
    ! to which processor it must be sent. It could be more efficient,
    ! but at least the index data are run-length encoded.
    DBG 'A'
    send_count(:,:) = 0
    CALL dbcsr_iterator_start (iter, matrix)
    dst_p = -1
    new_meta = .TRUE.
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
    t_xfer1 = m_walltime()
    t_count = t_xfer1 - t_count
    !
    DBG 'B'
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    DBG 'C'
    DBG 'send counts',send_count
    DBG 'recv counts',recv_count
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
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    DBG 'D'
    IF (dbg) THEN
       WRITE(*,*)'row_start_old',row_start_old
       WRITE(*,*)'row_start_new',row_start_new
       WRITE(*,*)'row_blk_size_new',row_blk_size_new
    ENDIF
    !
    t_fill = m_walltime()
    t_xfer1 = t_fill - t_xfer1
    !
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_new (data_block, data_type)
    CALL dbcsr_iterator_start (iter, matrix)
    row_prev = 0
    dst_p = -1
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, data_block, tr, blk,&
            row_size=row_size, col_size=col_size)
       IF (row .NE. row_prev) blk_col_prev_new = 0 ! Reset the "hint"
       !IF (tr) WRITE(*,*)"block at",row,col," is transposed"
       col_offset_old = 1
       DO col_int = n_col_reblocks (1,col),&
            n_col_reblocks(1,col) + n_col_reblocks(2,col)-1
          blk_col_new = col_reblocks(1, col_int)
          row_offset_old = 1
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
                row_rle = row_reblocks(2,row_int)
                col_rle = col_reblocks(2,col_int)
                row_offset_new = row_offset_old + row_start_old(row) - &
                     row_start_new(blk_row_new)
                col_offset_new = col_offset_old + col_start_old(col) - &
                     col_start_new(blk_col_new)
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
                !IF (.NOT. tr) THEN
                !      CALL dbcsr_block_copy(&
                !           send_data(sdp(dst_p):sdp(dst_p)+nze_rle-1),&
                !           r2_dp(&
                !           row_offset_old:row_offset_old+row_rle-1,&
                !           col_offset_old:col_offset_old+col_rle-1),&
                !           row_rle, col_rle)
                !ELSE
                !   CALL dbcsr_block_transpose(&
                !        send_data(sdp(dst_p):sdp(dst_p)+nze_rle-1),&
                !        r2_dp(&
                !        col_offset_old:col_offset_old+col_rle-1,&
                !        row_offset_old:row_offset_old+row_rle-1),&
                !        row_rle, col_rle)
                !ENDIF
                sdp(dst_p) = sdp(dst_p)+nze_rle
             ENDIF loc_ok
             row_offset_old = row_offset_old + row_reblocks(2, row_int)
          ENDDO ! row_int
          col_offset_old = col_offset_old + col_reblocks(2, col_int)
       ENDDO ! col_int
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    CALL dbcsr_data_clear_pointer (data_block)
    CALL dbcsr_data_release (data_block)

    ! Exchange the data and metadata structures.
    !
    t_xfer2 = m_walltime()
    t_fill = t_xfer2 - t_fill
    !
    DBG 'E'
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       CALL hybrid_alltoall_s1(&
            send_data%d%r_sp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%r_sp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_real_8)
       !CALL mp_alltoall(&
       !     send_data%d%r_dp(:), total_send_count(2,:), sd_disp(:)-1,&
       !     recv_data%d%r_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
       !     mp_group)
       CALL hybrid_alltoall_d1 (&
            send_data%d%r_dp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%r_dp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_complex_4)
       CALL hybrid_alltoall_c1(&
            send_data%d%c_sp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%c_sp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_complex_8)
       CALL hybrid_alltoall_z1(&
            send_data%d%c_dp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%c_dp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid matrix type",__LINE__,error)
    END SELECT
    DBG 'F'
    DBG 'send_data'
    CALL hybrid_alltoall_i1(send_meta(:), SUM(metalen*total_send_count(1,:)), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), SUM(metalen*total_recv_count(1,:)), metalen*total_recv_count(1,:), rm_disp(:)-1, mp_obj_new)
    DBG 'G'
    !
    t_unpack = m_walltime()
    t_xfer2 = t_unpack - t_xfer2
    !
    ! Now fill in the data.
    CALL dbcsr_work_create(redist,&
         nblks_guess=SUM(recv_count(1,:)),&
         sizedata_guess=SUM(recv_count(2,:)), work_mutable=.TRUE.,error=error)
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    CALL dbcsr_data_init (buff_data)
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_init (src_2d)
    CALL dbcsr_data_new (buff_data, dbcsr_type_1d_to_2d (data_type), &
         redist%m%max_rbs, redist%m%max_cbs)
    CALL dbcsr_data_new (data_block, dbcsr_type_1d_to_2d (data_type))
    CALL dbcsr_data_new (src_2d, dbcsr_type_1d_to_2d (data_type))

    !blk_p = 1
    !blk = 1
    blk_ps = 0
    blks = 0
    t_blk2 = 0.0_dp
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
          DBG 'blk: p',blks,src_p,'; row, col',&
               stored_row_new,stored_col_new,&
               '; blk_ps',blk_ps

          CALL dbcsr_data_clear_pointer (data_block)
          t_blk = m_walltime()
          CALL dbcsr_get_block_p(redist, stored_row_new, stored_col_new,&
               data_block, tr, found)
          t_blk2 = t_blk2 + (m_walltime() - t_blk)
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
             CALL dbcsr_data_set_pointer (src_2d,&
                  rsize=row_rle, csize=col_rle,&
                  pointee=recv_data, source_lb=data_offset_l)

             IF (.NOT. tr) THEN
                CALL dbcsr_data_set (dst=data_block,&
                     lb=row_offset_new, data_size=row_rle,&
                     lb2=col_offset_new, data_size2=col_rle,&
                     src=src_2d,&
                     source_lb=1, source_lb2=1)
                !CALL dbcsr_block_copy (&
                !     r2_dp(&
                !     row_offset_new:row_offset_new+row_rle-1,&
                !     col_offset_new:col_offset_new+col_rle-1),&
                !     recv_data(data_offset_l:data_offset_l+nze_rle-1),&
                !     row_rle, col_rle)
             ELSE
                CALL dbcsr_block_transpose(&
                     dst=data_block,&
                     lb=col_offset_new, lb2=row_offset_new,&
                     src=src_2d,&
                     row_size=row_rle, col_size=col_rle)
                !CALL dbcsr_block_transpose (&
                !     r2_dp(&
                !     col_offset_new:col_offset_new+col_rle-1,&
                !     row_offset_new:row_offset_new+row_rle-1),&
                !     recv_data(data_offset_l:data_offset_l+nze_rle-1),&
                !     row_rle, col_rle)
             ENDIF
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
    CALL dbcsr_data_clear_pointer(src_2d)
    CALL dbcsr_data_release (buff_data)
    CALL dbcsr_data_release (data_block)
    CALL dbcsr_data_release (src_2d)
    !
    IF (dbg) THEN
       WRITE(*,*)routineN//" Declared blocks=",redist%m%wms(1)%lastblk,&
            "actual=",blks
       WRITE(*,*)routineN//" Declared data size=",redist%m%wms(1)%datasize,&
            "actual=",blk_ps
    ENDIF
    DBG 'Finalizing redistribution'

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
    DBG 'Done redistributing'
    t_unpack = m_walltime() - t_unpack
    t_all = MAX(m_walltime() - t_all, 0.001_dp)
    !write(*,'(A,2(1X,F9.3),F5.1,"%",4(1X,I7))')&
    !     " Times total, lookup", t_all, t_blk2,&
    !     100.0_dp*t_blk2/t_all,&
    !     cnt_fnd, cnt_new, cnt_skip, redist%m%nblks
    !write(*,'(A,6(1X,F9.3))')" Times              ",&
    !     t_prolog, t_count, t_xfer1, t_fill, t_xfer2, t_unpack
    !write(*,*)" "
    !if (dbg) call dbcsr_print(redist)
    !call dbcsr_print(redist)
    !cs2 = dbcsr_checksum (redist)
    !WRITE(*,*)routineN//" Checksums=",cs1, cs2, cs1-cs2
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

    CHARACTER                                :: data_type
    CHARACTER(len=80)                        :: ifmt
    INTEGER :: blk, blk_ps, blks, col, col_size, dst_p, error_handler, &
      meta_l, mp_group, numnodes, nze, row, row_prev, row_size, src_p, &
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
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_block, recv_data, &
                                                send_data
    TYPE(dbcsr_distribution_obj)             :: dist_new, dist_old
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_mp_obj)                       :: mp_obj_new

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    IF (dbg) CALL dbcsr_print (matrix,error=error)
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
    dist_old = dbcsr_distribution (matrix)
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
    DBG 'A'
    send_count(:,:) = 0
    CALL dbcsr_iterator_start (iter, matrix)
    row_prev = 0
    dst_p = -1
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, blk, tr,&
            row_size=row_size, col_size=col_size)
       CALL dbcsr_get_stored_coordinates (redist%m,&
            row, col, tr, dst_p)
       IF (dbg) WRITE(*,*)routinen//" Sending",row,col,"to",dst_p
       nze = row_size*col_size
       send_count(1, dst_p) = send_count(1, dst_p) + 1
       send_count(2, dst_p) = send_count(2, dst_p) + nze
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    DBG 'B'
    CALL mp_alltoall(send_count, recv_count, 2, mp_group)
    DBG 'C'
    DBG 'send counts',send_count
    DBG 'recv counts',recv_count
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
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
    DBG 'D'
    IF (dbg) THEN
       WRITE(*,*)'row_blk_size_new',row_blk_size_new
    ENDIF
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_new (data_block, data_type)
    CALL dbcsr_iterator_start (iter, matrix)
    dst_p = -1
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, data_block, tr, blk)
       !IF (tr) WRITE(*,*)"block at",row,col," is transposed"
       CALL dbcsr_get_stored_coordinates (redist%m,&
            row, col, tr, dst_p)
       smp(dst_p) = smp(dst_p) + metalen
       IF (tr) THEN
          send_meta(smp(dst_p)) = -row
       ELSE
          send_meta(smp(dst_p)) = row
       ENDIF
       send_meta(smp(dst_p)+1) = col ! new blocked column
       nze = dbcsr_get_data_size (data_block)
       CALL dbcsr_data_set (send_data, lb=sdp(dst_p), data_size=nze,&
            src=data_block, source_lb=1)
       !send_data(sdp(dst_p):sdp(dst_p)+SIZE(r_dp)-1) &
       !     = r_dp(:)
       sdp(dst_p) = sdp(dst_p) + nze
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    CALL dbcsr_data_clear_pointer (data_block)
    ! Exchange the data and metadata structures.
    DBG 'E'
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       CALL hybrid_alltoall_s1(&
            send_data%d%r_sp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%r_sp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_real_8)
       !CALL mp_alltoall(&
       !     send_data%d%r_dp(:), total_send_count(2,:), sd_disp(:)-1,&
       !     recv_data%d%r_dp(:), total_recv_count(2,:), rd_disp(:)-1,&
       !     mp_group)
       CALL hybrid_alltoall_d1 (&
            send_data%d%r_dp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%r_dp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_complex_4)
       CALL hybrid_alltoall_c1(&
            send_data%d%c_sp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%c_sp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    CASE (dbcsr_type_complex_8)
       CALL hybrid_alltoall_z1(&
            send_data%d%c_dp(:), SUM(total_send_count(2,:)), total_send_count(2,:), sd_disp(:)-1,&
            recv_data%d%c_dp(:), SUM(total_recv_count(2,:)), total_recv_count(2,:), rd_disp(:)-1,&
            mp_obj_new)
    END SELECT
    !CALL mp_alltoall(send_data(:), total_send_count(2,:), sd_disp(:)-1,&
    !     recv_data(:), total_recv_count(2,:), rd_disp(:)-1, mp_group)
    DBG 'F'
    DBG 'send_meta'
    IF (dbg) WRITE(ifmt,'("(",I2,"(1X,I5))")') metalen
    IF (dbg) WRITE(*,ifmt) send_meta
    DBG 'send_data'
    !IF (dbg) WRITE(*,'(10(1X,F5.2))') send_data
    CALL hybrid_alltoall_i1(send_meta(:), SIZE(send_meta), metalen*total_send_count(1,:), sm_disp(:)-1,&
         recv_meta(:), SIZE(recv_meta), metalen*total_recv_count(1,:), rm_disp(:)-1, mp_obj_new)
    DBG 'recv_meta'
    IF (dbg) WRITE(*,ifmt) recv_meta
    DBG 'recv_data'
    !IF (dbg) WRITE(*,'(10(1X,F5.2))') recv_data
    DBG 'G'
    ! Now fill in the data.
    CALL dbcsr_work_create(redist,&
            SUM(recv_count(1,:)),&
            SUM(recv_count(2,:)), work_mutable=.FALSE., n=1, error=error)
    !
    DBG 'send disps data',sd_disp
    DBG 'send disps meta',sm_disp
    DBG 'recv disps data',rd_disp
    DBG 'recv disps meta',rm_disp
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
    IF (dbg) WRITE(*,*)routineN//" blk_p", redist%m%wms(1)%blk_p
    CALL dbcsr_data_clear_pointer (data_block)
    CALL dbcsr_data_release (data_block)
    !
    IF (dbg) THEN
       WRITE(*,*)routineN//" Declared blocks=",redist%m%wms(1)%lastblk,&
            "actual=",blks
       WRITE(*,*)routineN//" Declared data size=",redist%m%wms(1)%datasize,&
            "actual=",blk_ps
    ENDIF
    DBG 'Finalizing redistribution'
    CALL dbcsr_finalize(redist, error=error)
    CALL dbcsr_data_release (recv_data)
    CALL dbcsr_data_release (send_data)
    IF (dbg) CALL dbcsr_print (redist,error=error)
    DEALLOCATE(send_count)
    DEALLOCATE(recv_count)
    DEALLOCATE(sdp); DEALLOCATE(sd_disp)
    DEALLOCATE(smp); DEALLOCATE(sm_disp)
    DEALLOCATE(rd_disp)
    DEALLOCATE(rm_disp)
    DEALLOCATE(recv_meta)
    DEALLOCATE(send_meta)
    DBG 'Done redistributing'
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_redistribute





END MODULE dbcsr_transformations

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR work matrix utilities
!> \author  Urban Borstnik
!> \date    2010-02-18
!> \version 0.9
!>
!> <b>Modification history:</b>
!> - 2010-02-18 Copied from dbcsr_util
! *****************************************************************************
MODULE dbcsr_work_operations

  USE array_types,                     ONLY: array_data,&
                                             array_hold,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_release,&
                                             array_size
  USE btree_I8_k_cp2d_v,               ONLY: btree_2d_data_c => cp2d,&
                                             btree_add_c => btree_add,&
                                             btree_c => btree,&
                                             btree_get_c => btree_find
  USE btree_I8_k_dp2d_v,               ONLY: btree_2d_data_d => dp2d,&
                                             btree_add_d => btree_add,&
                                             btree_d => btree,&
                                             btree_get_d => btree_find
  USE btree_I8_k_sp2d_v,               ONLY: btree_2d_data_s => sp2d,&
                                             btree_add_s => btree_add,&
                                             btree_get_s => btree_find,&
                                             btree_s => btree
  USE btree_I8_k_zp2d_v,               ONLY: btree_2d_data_z => zp2d,&
                                             btree_add_z => btree_add,&
                                             btree_get_z => btree_find,&
                                             btree_z => btree
  USE dbcsr_block_operations,          ONLY: dbcsr_data_set
  USE dbcsr_data_operations,           ONLY: dbcsr_data_copyall
  USE dbcsr_error_handling
  USE dbcsr_index_operations,          ONLY: dbcsr_addto_index_array,&
                                             dbcsr_repoint_index,&
                                             dbcsr_sort_indices,&
                                             dbcsr_sort_many_indices
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_kinds,                     ONLY: dp,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE dbcsr_message_passing,           ONLY: mp_allocate
  USE dbcsr_methods,                   ONLY: &
       dbcsr_data_ensure_size, dbcsr_data_hold, dbcsr_data_init, &
       dbcsr_data_new, dbcsr_data_release, dbcsr_destroy_array, &
       dbcsr_destroy_image_dist, dbcsr_distribution_hold, &
       dbcsr_distribution_release, dbcsr_get_data_p, dbcsr_get_data_size, &
       dbcsr_get_data_type, dbcsr_get_info, dbcsr_get_num_blocks, dbcsr_init, &
       dbcsr_is_initialized, dbcsr_modify_lock, dbcsr_modify_unlock, &
       dbcsr_mutable_destroy, dbcsr_mutable_init, dbcsr_mutable_new, &
       dbcsr_mutable_release, dbcsr_uses_special_memory, dbcsr_valid_index, &
       dbcsr_wm_use_mutable
  USE dbcsr_ptr_util,                  ONLY: addto_array,&
                                             ensure_array_size
  USE dbcsr_types,                     ONLY: &
       dbcsr_data_obj, dbcsr_distribution_obj, dbcsr_iterator, &
       dbcsr_magic_number, dbcsr_meta_size, dbcsr_num_slots, dbcsr_obj, &
       dbcsr_repl_col, dbcsr_repl_full, dbcsr_repl_none, dbcsr_repl_row, &
       dbcsr_slot_blk_p, dbcsr_slot_col_i, dbcsr_slot_nblkcols_local, &
       dbcsr_slot_nblkcols_total, dbcsr_slot_nblkrows_local, &
       dbcsr_slot_nblkrows_total, dbcsr_slot_nblks, &
       dbcsr_slot_nfullcols_local, dbcsr_slot_nfullcols_total, &
       dbcsr_slot_nfullrows_local, dbcsr_slot_nfullrows_total, &
       dbcsr_slot_nze, dbcsr_slot_row_p, dbcsr_slot_size, dbcsr_type, &
       dbcsr_type_complex_4, dbcsr_type_complex_8, dbcsr_type_real_4, &
       dbcsr_type_real_8, dbcsr_type_real_default, dbcsr_work_type
  USE dbcsr_util,                      ONLY: convert_sizes_to_offsets,&
                                             dbcsr_set_debug,&
                                             dbcsr_unpack_i8_2i4,&
                                             dbcsr_verify_matrix,&
                                             meta_from_dist

  !$ USE OMP_LIB

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_work_operations'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034


  PUBLIC :: dbcsr_create, dbcsr_work_create, dbcsr_finalize
  PUBLIC :: dbcsr_work_destroy, dbcsr_add_wm_from_matrix,&
            add_work_coordinate


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
!> \brief Creates a matrix, allocating the essentials.
!> \par The matrix itself is allocated, as well as the essential parts of
!>      the index. When passed the nze argument, the data is also allocated
!>      to that size.
!> \param[in,out] matrix      new matrix
!> \param[in] dist            distribution_2d distribution
!> \param[in] matrix_type     'N' for normal, 'T' for transposed, 'S' for
!>                            symmetric, and 'A' for antisymmetric
!> \param[in] nblks           (optional) number of blocks
!> \param[in] nze             (optional) number of elements
!> \param[in] data_type       type of data from [rRcC] for single/double
!>                            precision real/complex, default is 'R'
!> \param[in] special         (optional) allocate indices and data using
!>                            special memory
!> \param[in] reuse           (optional) reuses an existing matrix, default
!>                            is to create a fresh one
!> \param[in] mutable_work    uses the mutable data for working and not the
!>                            append-only data; default is append-only
!> \param[in] replication_type     replication to be used for this matrix;
!>                                 default is dbcsr_repl_none
!>                                 \see dbcsr_types.F
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_create(matrix, name, dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, data_type, special, reuse,&
       mutable_work, replication_type, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    CHARACTER(len=*), INTENT(IN)             :: name
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    CHARACTER, INTENT(IN)                    :: matrix_type
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_size, col_blk_size
    INTEGER, INTENT(IN), OPTIONAL            :: nblks, nze
    CHARACTER, INTENT(IN), OPTIONAL          :: data_type
    LOGICAL, INTENT(IN), OPTIONAL            :: special, reuse, mutable_work
    CHARACTER, INTENT(IN), OPTIONAL          :: replication_type
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, my_nblks, &
                                                my_nze, stat
    INTEGER, DIMENSION(:), POINTER           :: col_blk_offset, row_blk_offset
    INTEGER, DIMENSION(dbcsr_meta_size)      :: new_meta
    LOGICAL                                  :: hijack

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    hijack = .FALSE.
    IF (PRESENT (reuse)) THEN
       hijack = reuse
    ELSE
       IF (matrix%m%initialized.EQ.dbcsr_magic_number) THEN
          ! Reuse matrix only if has actually been allocated.
          IF (ASSOCIATED (matrix%m%index)) THEN
             hijack = .TRUE.
          ELSE
             hijack = .FALSE.
          ENDIF
       ELSE
          CALL dbcsr_assert (matrix%m%initialized, 'EQ', 0,&
               dbcsr_fatal_level, dbcsr_caller_error, routineN,&
               "Matrix may not have been initialized with dbcsr_init",__LINE__,error)
          hijack = .FALSE.
       ENDIF
    ENDIF
    IF (.NOT.hijack) THEN
       CALL dbcsr_init (matrix%m)
       matrix%m%refcount = 1
    ENDIF
    ! Mark matrix index as having an invalid index.
    matrix%m%valid = .FALSE.
    matrix%m%name = name
    ! Sets the type of matrix building/modifying work structures.
    IF (PRESENT (mutable_work)) THEN
       matrix%m%work_mutable = mutable_work
    ELSE
       matrix%m%work_mutable = .FALSE.
    ENDIF
    ! Sets the correct data type.
    IF (PRESENT (data_type)) THEN
       SELECT CASE (data_type)
       CASE (dbcsr_type_real_4)
          matrix%m%data_type = dbcsr_type_real_4
       CASE (dbcsr_type_real_8)
          matrix%m%data_type = dbcsr_type_real_8
       CASE (dbcsr_type_complex_4)
          matrix%m%data_type = dbcsr_type_complex_4
       CASE (dbcsr_type_complex_8)
          matrix%m%data_type = dbcsr_type_complex_8
       CASE DEFAULT
         CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
              routineN, "Invalid matrix type",__LINE__,error)
      END SELECT
    ELSE
       matrix%m%data_type = dbcsr_type_real_default
    ENDIF
    IF (hijack) THEN
       ! Release/deallocate elements that are replaced or not needed
       ! by the new matrix. This is similar to what dbcsr_destroy
       ! does, except that it keeps the index and data.
       CALL array_release (matrix%m%row_blk_size)
       CALL array_release (matrix%m%col_blk_size)
       CALL array_release (matrix%m%row_blk_offset)
       CALL array_release (matrix%m%col_blk_offset)
       CALL dbcsr_distribution_release (matrix%m%dist)
       IF (ASSOCIATED (matrix%m%wms)) THEN
          CALL dbcsr_work_destroy_all(matrix%m)
       ENDIF
    ELSE
       ! Invalidate index
       NULLIFY(matrix%m%index)
       ! Invalidate data
       CALL dbcsr_data_init (matrix%m%data_area)
    ENDIF
    ! These are always invalidated.
    NULLIFY(matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p)
    matrix%m%row_blk_size = row_blk_size
    CALL array_hold (matrix%m%row_blk_size)
    IF (array_size (matrix%m%row_blk_size) .GT. 0) THEN
       matrix%m%max_rbs = MAXVAL (array_data (matrix%m%row_blk_size))
    ELSE
       matrix%m%max_rbs = 0
    ENDIF
    matrix%m%col_blk_size = col_blk_size
    CALL array_hold (matrix%m%col_blk_size)
    IF (array_size (matrix%m%col_blk_size) .GT. 0) THEN
       matrix%m%max_cbs = MAXVAL (array_data (matrix%m%col_blk_size))
    ELSE
       matrix%m%max_cbs = 0
    ENDIF

    ! initialize row/col offsets
    ALLOCATE(row_blk_offset(array_size (matrix%m%row_blk_size)+1))
    ALLOCATE(col_blk_offset(array_size (matrix%m%col_blk_size)+1))
    CALL convert_sizes_to_offsets(array_data(matrix%m%col_blk_size), col_blk_offset)
    CALL convert_sizes_to_offsets(array_data(matrix%m%row_blk_size), row_blk_offset)
    CALL array_new(matrix%m%col_blk_offset, col_blk_offset, gift=.TRUE.)
    CALL array_new(matrix%m%row_blk_offset, row_blk_offset, gift=.TRUE.)

    matrix%m%dist = dist
    CALL dbcsr_distribution_hold (matrix%m%dist)
    ! Set up some data.
    CALL meta_from_dist (new_meta, dist, row_blk_size, col_blk_size)
    matrix%m%nblkrows_total  = new_meta(dbcsr_slot_nblkrows_total )
    matrix%m%nblkcols_total  = new_meta(dbcsr_slot_nblkcols_total )
    matrix%m%nfullrows_total = new_meta(dbcsr_slot_nfullrows_total)
    matrix%m%nfullcols_total = new_meta(dbcsr_slot_nfullcols_total)
    matrix%m%nblkrows_local  = new_meta(dbcsr_slot_nblkrows_local )
    matrix%m%nblkcols_local  = new_meta(dbcsr_slot_nblkcols_local )
    matrix%m%nfullrows_local = new_meta(dbcsr_slot_nfullrows_local)
    matrix%m%nfullcols_local = new_meta(dbcsr_slot_nfullcols_local)
    my_nze = 0; IF (PRESENT (nze)) my_nze = nze
    my_nblks = 0; IF (PRESENT (nblks)) my_nblks = nblks
    matrix%m%nblks = 0
    matrix%m%nze = 0
    IF (PRESENT (special)) THEN
       matrix%m%special_memory = special
    ELSE
       matrix%m%special_memory = .FALSE.
    ENDIF
    IF (PRESENT (replication_type)) THEN
       CALL dbcsr_assert (replication_type .EQ. dbcsr_repl_none&
            .OR. replication_type .EQ. dbcsr_repl_full&
            .OR. replication_type .EQ. dbcsr_repl_row&
            .OR. replication_type .EQ. dbcsr_repl_col,&
            dbcsr_failure_level, dbcsr_wrong_args_error, routineN,&
            "Invalid replication type '"//replication_type//"'",__LINE__,error)
       CALL dbcsr_assert (replication_type .NE. dbcsr_repl_row&
            .AND. replication_type .NE. dbcsr_repl_col,&
            dbcsr_warning_level, dbcsr_unimplemented_error_nr, routineN,&
            "Row and column replication not fully supported",__LINE__,error)
       matrix%m%replication_type = replication_type
    ELSE
       matrix%m%replication_type = dbcsr_repl_none
    ENDIF
    !
    ! Setup the data
    IF (.NOT. hijack) THEN
       CALL dbcsr_data_new (matrix%m%data_area, matrix%m%data_type, my_nze,&
            special_memory=matrix%m%special_memory)
       IF (matrix%m%special_memory) THEN
          CALL mp_allocate(matrix%m%index, dbcsr_num_slots, stat=stat)
          CALL dbcsr_assert (stat == 0,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "matrix%m%index",__LINE__,error)
       ELSE
          ALLOCATE(matrix%m%index(dbcsr_num_slots), stat=stat)
          CALL dbcsr_assert (stat == 0,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "matrix%m%index",__LINE__,error)
       ENDIF
    ENDIF
    CALL dbcsr_assert (LBOUND (matrix%m%index, 1) .LE. 1&
         .AND. UBOUND (matrix%m%index, 1) .GE. dbcsr_num_slots,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Index is not large enough",__LINE__,error)
    matrix%m%index(1:dbcsr_num_slots) = 0
    matrix%m%index(dbcsr_slot_size) = dbcsr_num_slots
    !
    matrix%m%symmetry = .FALSE.
    matrix%m%negate_real = .FALSE.
    matrix%m%negate_imaginary = .FALSE.
    matrix%m%transpose = .FALSE.
    SELECT CASE (matrix_type)
    CASE ('N', 'n')
    CASE ('T', 't')
       matrix%m%transpose = .TRUE.
    CASE ('S', 's')
       matrix%m%symmetry = .TRUE.
    CASE ('A', 'a')
       matrix%m%symmetry = .TRUE.
       matrix%m%negate_real = .TRUE.
    CASE DEFAULT
       CALL dbcsr_assert(.FALSE., dbcsr_failure_level,&
            dbcsr_wrong_args_error, routineP, "Invalid matrix type.",__LINE__,error)
    END SELECT
    NULLIFY (matrix%m%predistributed)
    matrix%m%bcsc = .FALSE.
!$  CALL OMP_INIT_LOCK (matrix%m%modification_lock)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_create


! *****************************************************************************
!> \brief Initializes one work matrix
!> \param[out] wm             initialized work matrix
!> \param[in] nblks_guess     (optional) estimated number of blocks
!> \param[in] sizedata_guess  (optional) estimated size of data
! *****************************************************************************
  SUBROUTINE dbcsr_init_wm (wm, data_type, nblks_guess, sizedata_guess, error)
    TYPE(dbcsr_work_type), INTENT(OUT)       :: wm
    CHARACTER, INTENT(IN)                    :: data_type
    INTEGER, INTENT(IN), OPTIONAL            :: nblks_guess, sizedata_guess
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_wm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, nblks, stat

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    wm%lastblk = 0
    wm%datasize = 0
    ! Index
    IF(PRESENT(nblks_guess)) THEN
       nblks = nblks_guess
       ALLOCATE(wm%row_i(nblks), stat=stat)
       CALL dbcsr_assert (stat == 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "wm%row_i",__LINE__,error)
       ALLOCATE(wm%col_i(nblks), stat=stat)
       CALL dbcsr_assert (stat == 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "wm%col_i",__LINE__,error)
       ALLOCATE(wm%blk_p(nblks), stat=stat)
       CALL dbcsr_assert (stat == 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "wm%blk_p",__LINE__,error)
    ELSE
       NULLIFY (wm%row_i, wm%col_i, wm%blk_p)
       !nblks = CEILING (REAL (matrix%nblkrows_local * matrix%nblkcols_local)&
       !     / REAL (dbcsr_mp_numnodes (dbcsr_distribution_mp (matrix%dist))))
    ENDIF
    ! Data
    CALL dbcsr_data_init (wm%data_area)
    IF(PRESENT(sizedata_guess)) THEN
       CALL dbcsr_data_new (wm%data_area, data_type,&
            data_size=sizedata_guess)
    ELSE
       CALL dbcsr_data_new (wm%data_area, data_type)
    ENDIF
    CALL dbcsr_mutable_init (wm%mutable)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_init_wm


! *****************************************************************************
!> \brief Creates a the working matrix(es) for a DBCSR matrix.
!> \param[out] matrix         new matrix
!> \param[in] nblks_guess     (optional) estimated number of blocks
!> \param[in] sizedata_guess  (optional) estimated size of data
!> \param[in] n               (optional) number work matrices to create,
!>                            default is 1
!> \param[in] work_mutable    (optional) use mutable work type, default is
!>                            what was specified in create
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_work_create(matrix, nblks_guess, sizedata_guess, n,&
       work_mutable, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN), OPTIONAL            :: nblks_guess, sizedata_guess, n
    LOGICAL, INTENT(in), OPTIONAL            :: work_mutable
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_work_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, iw, nw, ow
    LOGICAL                                  :: wms_new, wms_realloc
    TYPE(dbcsr_work_type), DIMENSION(:), &
      POINTER                                :: wms

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (matrix%m%initialized,'EQ',dbcsr_magic_number,&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Can not create work matrix for nonexisting matrix object.",__LINE__,error)
    IF (PRESENT (n)) THEN
       nw = n
    ELSE
       nw = 1
!$     nw = omp_get_max_threads()
    ENDIF
!$  CALL dbcsr_modify_lock (matrix)
    wms_new = .NOT. ASSOCIATED (matrix%m%wms)
    wms_realloc = .FALSE.
    nw = 1
!$    nw = omp_get_max_threads()
    IF (ASSOCIATED (matrix%m%wms)) THEN
       ow = SIZE(matrix%m%wms)
       CALL dbcsr_assert (ow, 'GE', nw, dbcsr_warning_level,&
            dbcsr_internal_error, routineN,&
            "Number of work matrices less than threads.",__LINE__,error)
       IF (ow .LT. nw) wms_realloc = .TRUE.
    ENDIF
    IF (PRESENT (work_mutable)) THEN
       matrix%m%work_mutable = work_mutable
    ENDIF
    IF (wms_realloc) THEN
       ALLOCATE (wms(nw))
       wms(1:ow) = matrix%m%wms(1:ow)
       DEALLOCATE (matrix%m%wms)
       matrix%m%wms => wms
       DO iw = ow+1, nw
          CALL dbcsr_init_wm (matrix%m%wms(iw), matrix%m%data_type,&
               nblks_guess=nblks_guess, sizedata_guess=sizedata_guess, error=error)
          IF (matrix%m%work_mutable) &
               CALL dbcsr_mutable_new (matrix%m%wms(iw)%mutable,&
               dbcsr_get_data_type (matrix))
       END DO
    ENDIF
    IF (wms_new) THEN
       ALLOCATE (matrix%m%wms(nw))
       DO iw = 1, nw
          CALL dbcsr_init_wm (matrix%m%wms(iw), matrix%m%data_type,&
               nblks_guess=nblks_guess, sizedata_guess=sizedata_guess, error=error)
          IF (matrix%m%work_mutable) &
               CALL dbcsr_mutable_new (matrix%m%wms(iw)%mutable,&
               dbcsr_get_data_type (matrix))
       END DO
    ENDIF
!$OMP FLUSH (matrix)
!$  CALL dbcsr_modify_unlock (matrix)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_work_create


! *****************************************************************************
!> \brief Creates the final dbcsr_type matrix from the working matrix.
!>
!> If a matrix has a valid index, then nothing is done.
!> \note The columns and rows must be ordered! We don't support appending
!>       data to an existing matrix.
!> \param[in,out] matrix      final matrix
!> \param[in] resort          whether the indices should be sorted, default
!>                            is true
!> \param[in] reshuffle       whether the data should be reshuffled,
!>                            default is false
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_finalize(matrix, resort, reshuffle, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: resort, reshuffle
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: all_nblks, all_nze, &
                                                error_handler, i, j
    LOGICAL                                  :: multidata, sort_data, &
                                                sort_index

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (.NOT. ASSOCIATED (matrix%m%predistributed),&
         dbcsr_warning_level,&
         dbcsr_internal_error, routineN, "Submatrices exist",__LINE__,error)
    IF (ASSOCIATED (matrix%m%predistributed)) THEN
          CALL dbcsr_destroy_array (matrix%m%predistributed,error)
          CALL dbcsr_destroy_image_dist (matrix%m%predistributed%image_dist)
          DEALLOCATE (matrix%m%predistributed)
    ENDIF
    NULLIFY (matrix%m%predistributed)
    CALL dbcsr_assert (dbcsr_is_initialized (matrix),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Can not finalize uninitialized matrix.",__LINE__,error)
    IF (dbcsr_valid_index(matrix)) THEN
       IF (dbg) WRITE(*,'(2(1X,A))')&
            "No need to finalize a valid matrix, skipping."
!vw       CALL dbcsr_assert (.NOT. ASSOCIATED (matrix%m%wms), dbcsr_warning_level,&
!vw            dbcsr_internal_error, routineN,&
!vw            "Valid matrix should not contain work data")
       IF (ASSOCIATED (matrix%m%wms)) &
            CALL dbcsr_work_destroy_all(matrix%m)
       CALL dbcsr_error_stop(error_handler, error)
       RETURN
    ENDIF
    IF (dbg) WRITE(*,'(A,1X,A,2(1X,I9))')routineN//' finalizing matrix ',&
         matrix%m%name, matrix%m%nze, matrix%m%nblks
    !CALL dbcsr_assert (matrix%m%nze.EQ.0 .EQV. matrix%m%nblks.EQ.0,&
    !     dbcsr_warning_level, dbcsr_unimplemented_error_nr, routineN,&
    !     "Matrix data size not consistent with number of blocks.")
    ! If there is data in the matrix proper, create a work matrix from it.
    IF (matrix%m%nze .GT. 0) THEN
       CALL dbcsr_add_wm_from_matrix (matrix, 1, error=error)
    ENDIF
    IF (PRESENT (reshuffle)) THEN
       sort_data = reshuffle
    ELSE
       sort_data = .FALSE.
    ENDIF
    sort_index = .TRUE.
    IF (PRESENT (resort)) sort_index = resort
    IF (ASSOCIATED (matrix%m%wms)) THEN
       multidata = SIZE (matrix%m%wms) .GT. 1
    ELSE
       multidata = .FALSE.
    ENDIF
    IF (.NOT. (matrix%m%nze .GT. 0 .OR. ASSOCIATED (matrix%m%wms))) THEN
       IF (dbg) THEN
          CALL dbcsr_assert (matrix%m%nze .GT. 0 .OR. ASSOCIATED (matrix%m%wms),&
               dbcsr_warning_level, dbcsr_caller_error, routineN,&
               "Nothing to be done.",__LINE__,error)
       ENDIF
       ! Should short-circut and just create an empty index.
       CALL dbcsr_work_create (matrix, error=error)
    ENDIF
    CALL dbcsr_assert (ASSOCIATED (matrix%m%wms), dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Can not finalize matrix containing no work matrices.",__LINE__,error)
    CALL dbcsr_assert (SIZE (matrix%m%wms),'GT', 0, dbcsr_failure_level,&
         dbcsr_caller_error, routineN,&
         "Can not finalize matrix containing no work matrices.",__LINE__,error)
    ! If there are multiple work matrices, the indices have to be sorted to
    ! enable merging the index.
    IF (sort_index .OR. multidata) THEN
       IF (multidata) THEN
          CALL dbcsr_sort_many_indices (matrix%m)
          sort_data = .TRUE.
       ELSE
          IF (.NOT.dbcsr_wm_use_mutable (matrix%m%wms(1))&
               .AND. matrix%m%wms(1)%lastblk .GT. 0) THEN
             CALL dbcsr_sort_indices(matrix%m%wms(1)%lastblk,&
                  matrix%m%wms(1)%row_i,&
                  matrix%m%wms(1)%col_i, matrix%m%wms(1)%blk_p)
          ENDIF
       ENDIF
    ENDIF
    ! Count how large the new matrix will have to be.
    all_nblks = 0
    all_nze = 0
    DO i = 1, SIZE (matrix%m%wms)
       all_nblks = all_nblks + matrix%m%wms(i)%lastblk
       all_nze = all_nze + matrix%m%wms(i)%datasize
    ENDDO
    IF (dbg) THEN
       DO i = 1, SIZE (matrix%m%wms)
          WRITE(*,'(A,4(1X,I4))')"Sorted work matrix",i, matrix%m%wms(i)%lastblk,&
               matrix%m%wms(i)%datasize,&
               dbcsr_get_data_size (matrix%m%wms(i)%data_area)
          IF (.NOT.dbcsr_wm_use_mutable (matrix%m%wms(i))) THEN
             DO j = 1, matrix%m%wms(i)%lastblk
                WRITE(*,'(A,I4,":",3(1X,I5))')"Sorted",j,&
                     matrix%m%wms(i)%row_i(j),&
                     matrix%m%wms(i)%col_i(j), matrix%m%wms(i)%blk_p(j)
             ENDDO
          ENDIF
       ENDDO
    ENDIF
    matrix%m%nblks = all_nblks
    matrix%m%nze = all_nze
    matrix%m%index(dbcsr_slot_nblks) = matrix%m%nblks
    matrix%m%index(dbcsr_slot_nze) = matrix%m%nze
    DBG routineP//' index',matrix%m%index(:dbcsr_num_slots)
    DBG routineP//' # blocks=',matrix%m%nblks, all_nblks
    DBG routinep//' # nze=',matrix%m%nze, all_nze
    ! We first figure out how large the indices should be and allocate them.
    CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_row_p,&
         reservation=matrix%m%nblkrows_total+1, extra=all_nblks*2,&
         error=error)
    CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_col_i,&
         reservation=all_nblks, error=error)
    CALL dbcsr_addto_index_array(matrix%m, dbcsr_slot_blk_p,&
         reservation=all_nblks, error=error)
    CALL dbcsr_repoint_index(matrix%m)
    matrix%m%row_p(1) = 0
    !CALL dbcsr_data_ensure_size (matrix%m%data_area, all_nze)
    ! Now merge the index and data.
    matrix%m%nze = all_nze ! Specifies size of data array to be allocated
    CALL dbcsr_merge_data (matrix%m, sort_data, error=error)
    !
    ! Clean up.
    CALL dbcsr_work_destroy_all(matrix%m)
    matrix%m%valid = .TRUE.
    IF (dbg) CALL dbcsr_verify_matrix (matrix, error=error)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_finalize


! *****************************************************************************
!> \brief Creates a work matrix from the data present in a finalized matrix.
!> \param[in,out] matrix      DBCSR matrix
!> \param[in] nwms            Number of work matrices (i.e., threads)
! *****************************************************************************
  SUBROUTINE dbcsr_add_wm_from_matrix(matrix, nwms, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: nwms
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_wm_from_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, i, n
    TYPE(dbcsr_work_type), DIMENSION(:), &
      POINTER                                :: new_wms

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    DBG "Adding to an existing matrix, wee!"
    IF (ASSOCIATED (matrix%m%wms)) THEN
       n = SIZE (matrix%m%wms)
       IF (n .EQ. 0) DEALLOCATE (matrix%m%wms)
    ELSE
       n = 0
    ENDIF
    IF (n .EQ. 0) THEN
       ! Just create it
       CALL dbcsr_work_create (matrix, matrix%m%nblks, matrix%m%nze, n=nwms,&
            work_mutable=.FALSE., error=error)
       CALL dbcsr_fill_wm_from_matrix (matrix%m%wms, matrix%m, nwms=nwms,&
            offset=0, error=error)
    ELSE
       ALLOCATE (new_wms (1:n+nwms))
       new_wms(1:n) = matrix%m%wms(1:n)
       DEALLOCATE (matrix%m%wms)
       matrix%m%wms => new_wms
       DO i = 1, nwms
          CALL dbcsr_init_wm (matrix%m%wms(n+i), dbcsr_get_data_type (matrix),&
               nblks_guess=dbcsr_get_num_blocks (matrix), error=error)
       ENDDO
       CALL dbcsr_fill_wm_from_matrix (matrix%m%wms, matrix%m,&
            nwms=nwms, offset=n, error=error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_add_wm_from_matrix


! *****************************************************************************
!> \brief Fills index and data of the work matrix from the
!>        previously-finalized one.
!> \param[out] wm             the work matrix to fill
!> \param[in,out] matrix      DBCSR matrix
!> \param[in] nwms            number of work matrices
!> \param[in] offset          first work matrix to use
! *****************************************************************************
  SUBROUTINE dbcsr_fill_wm_from_matrix(wm, matrix, nwms, offset,error)
    TYPE(dbcsr_work_type), DIMENSION(:), &
      INTENT(INOUT)                          :: wm
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, INTENT(IN)                      :: nwms, offset
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_fill_wm_from_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_p, col, &
                                                error_handler, i, nblks, row, &
                                                which_wm
    INTEGER, DIMENSION(:), POINTER           :: tdist
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_obj)                          :: m

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    nblks = matrix%nblks
    m%m = matrix
    IF (nwms .GT. 1) THEN
       tdist => array_data (matrix%dist%d%thread_dist)
       CALL dbcsr_assert (ASSOCIATED (tdist), dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "thread distribution must be defined in threaded env.",__LINE__,error)
       WRITE(*,*)routineN//" Threading fill", tdist
    ENDIF
    CALL dbcsr_iterator_start (iter, m)
    IF (dbg) WRITE(*,*)routineN//" There are",nblks,"blocks."
    i = 1
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       IF (dbg) WRITE(*,*)routineN//" Getting block", i,"/",nblks
       CALL dbcsr_iterator_next_block (iter, row, col, blk=blk, blk_p=blk_p)
       IF (dbg) WRITE(*,*)routineN//" With coor", row, col, "at", blk_p
       !which_wm = MOD (row, nwms) + 1 + offset
       IF (nwms .GT. 1) THEN
          which_wm = tdist (row) + offset
       ELSE
          which_wm = offset+1
       ENDIF
       CALL add_work_coordinate (wm(which_wm), row, col, blk_p, error=error)
       i = i+1
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    !wm%data_area = matrix%data_area
    !CALL dbcsr_data_hold (wm%data_area)
    DO which_wm = 1+offset, nwms+offset
       CALL dbcsr_data_ensure_size (wm(which_wm)%data_area,&
            dbcsr_get_data_size (matrix%data_area), error=error)
       CALL dbcsr_data_copyall (wm(which_wm)%data_area, matrix%data_area)
       wm(which_wm)%datasize = dbcsr_get_data_size (wm(which_wm)%data_area)
    ENDDO
    IF (dbg) WRITE(*,*)routineN//" Size is",wm%datasize, wm%lastblk
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_fill_wm_from_matrix


! *****************************************************************************
!> \brief Deallocates and destroys a work matrix.
!> \param[in,out] wm          work matrix
! *****************************************************************************
  SUBROUTINE dbcsr_work_destroy(wm)
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: wm

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_work_destroy', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (wm%row_i)) THEN
       DEALLOCATE(wm%row_i)
       NULLIFY (wm%row_i)
    ENDIF
    IF (ASSOCIATED (wm%col_i)) THEN
       DEALLOCATE(wm%col_i)
       NULLIFY (wm%row_i)
    ENDIF
    IF (ASSOCIATED (wm%blk_p)) THEN
       DEALLOCATE(wm%blk_p)
       NULLIFY (wm%blk_p)
    ENDIF
    CALL dbcsr_data_release (wm%data_area)
    CALL dbcsr_mutable_destroy (wm%mutable)
  END SUBROUTINE dbcsr_work_destroy


! *****************************************************************************
!> \brief Deallocates and destroys a work matrix.
!> \param[in,out] wm          work matrix
!> \param[in,out] error       cp2k error
!> \param keepdata    do not deallocate data
!> \param keepfinal   do not destroy the final, non-work matrix
!> \param keepfinaldata       do not destroy the data in the final,
!>                            non-work matrix
! *****************************************************************************
  SUBROUTINE dbcsr_work_destroy_all(m)
    TYPE(dbcsr_type), INTENT(INOUT)          :: m

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_work_destroy_all', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (m%wms), dbcsr_warning_level, dbcsr_caller_error,&
         routineN, "Want to destroy nonexisting work matrices.",__LINE__,error)
    IF (ASSOCIATED (m%wms)) THEN
       DO i = 1, SIZE (m%wms)
          CALL dbcsr_work_destroy (m%wms(i))
       ENDDO
       DEALLOCATE (m%wms)
       NULLIFY (m%wms)
    ENDIF
  END SUBROUTINE dbcsr_work_destroy_all



! *****************************************************************************
!> \brief Adds a coordinate (or other data) into a work matrix's row_i and
!>        col_i arrays and returns its position.
!> \note  Uses the matrix%lastblk to keep track of the current position.
!> \param[in,out] matrix      work matrix
!> \param[in] row,col         row, col data to add
!> \param[in] blk   (optional) block pointer to add
!> \param[out] index          (optional) saved position
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE add_work_coordinate(matrix, row, col, blk, index, error)
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: matrix
    INTEGER, INTENT(IN)                      :: row, col
    INTEGER, INTENT(IN), OPTIONAL            :: blk
    INTEGER, INTENT(OUT), OPTIONAL           :: index
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_work_coordinate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    dbg = .FALSE.
    matrix%lastblk = matrix%lastblk+1
    CALL ensure_array_size(matrix%row_i, ub=matrix%lastblk,&
         factor=default_resize_factor,error=error)
    CALL ensure_array_size(matrix%col_i, ub=matrix%lastblk,&
         factor=default_resize_factor,error=error)
    matrix%row_i(matrix%lastblk) = row
    matrix%col_i(matrix%lastblk) = col
    IF (PRESENT(blk)) THEN
       CALL ensure_array_size(matrix%blk_p, ub=matrix%lastblk,&
            factor=default_resize_factor,error=error)
       matrix%blk_p(matrix%lastblk) = blk
    ENDIF
    IF(dbg.AND.PRESENT(blk))&
         WRITE(*,*)routineP//' Adding',row,col,blk,'at',matrix%lastblk
    IF (dbg.AND.bcsr_verbose) THEN
       WRITE(*,*)routineP//' row_i=',matrix%row_i(1:matrix%lastblk)
       WRITE(*,*)routineP//' col_i=',matrix%col_i(1:matrix%lastblk)
       WRITE(*,*)routineP//' blk_p=',matrix%blk_p(1:matrix%lastblk)
    ENDIF
    IF (PRESENT (index)) index = matrix%lastblk
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE add_work_coordinate



! *****************************************************************************
!> \brief Merge data from sorted work matrices indices into the final matrix.
!>
!> \param[in,out] matrix      matrix to work on
!> \param[in] sort_data       (optional) whether to sort data or copy/move
! *****************************************************************************
  SUBROUTINE dbcsr_merge_data(matrix, sort_data, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: sort_data
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_merge_data', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, blk_p, error_handler, i, nwms, nze, o_blk_p, prev_row, &
      row, smallest, smallest_col, smallest_row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: wm_blk_p
    INTEGER, DIMENSION(:), POINTER           :: cbs, rbs
    LOGICAL                                  :: srt, wspecial
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: finished, tree_based

!### outfactor tree-to-work conversion 

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    srt = .FALSE.
    IF (PRESENT (sort_data)) srt = sort_data
    rbs => array_data (matrix%row_blk_size)
    cbs => array_data (matrix%col_blk_size)
    nwms = SIZE (matrix%wms)
    ALLOCATE (wm_blk_p(nwms), finished(nwms), tree_based(nwms))
    wm_blk_p(:) = 1
    ! Convert mutable data type into standard work arrays
    DO i = 1, nwms
       tree_based(i) = dbcsr_wm_use_mutable (matrix%wms(i))
       IF (tree_based(i)) THEN
          SELECT CASE (matrix%data_type)
          CASE (dbcsr_type_real_4)
             CALL tree_to_linear_s (matrix%wms(i), error=error)
          CASE (dbcsr_type_real_8)
             CALL tree_to_linear_d (matrix%wms(i), error=error)
          CASE (dbcsr_type_complex_4)
             CALL tree_to_linear_c (matrix%wms(i), error=error)
          CASE (dbcsr_type_complex_8)
             CALL tree_to_linear_z (matrix%wms(i), error=error)
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
                  routineN, "Invalid data type",__LINE__,error)
          END SELECT
       ENDIF
    ENDDO
    FORALL (i = 1:nwms)
       finished(i) = wm_blk_p(i) .GT. matrix%wms(i)%lastblk
    END FORALL
    IF (nwms .GT. 1) srt = .TRUE.
    !IF (matrix%special_memory) srt = .TRUE.
    DO i = 1, nwms
       wspecial = matrix%special_memory
       IF (dbcsr_get_data_size (matrix%wms(i)%data_area) .GT. 0) THEN
          wspecial = dbcsr_uses_special_memory (matrix%wms(i)%data_area)
          IF (matrix%special_memory .NEQV. wspecial) THEN
             srt = .TRUE.
          ENDIF
       ENDIF
    ENDDO
    IF (.NOT. srt) THEN
       CALL dbcsr_data_release (matrix%data_area)
       matrix%data_area = matrix%wms(1)%data_area
       CALL dbcsr_data_hold (matrix%data_area)
    ELSE
       CALL dbcsr_data_ensure_size (matrix%data_area, matrix%nze, error=error)
    ENDIF
    IF (dbg) THEN
       DO i = 1, nwms
          WRITE(*,'(1X,A,1X,I4,1X,A,3(1X,I7))')routineN//" size of wm",i,"is",&
               matrix%wms(i)%lastblk,&
               matrix%wms(i)%datasize,&
               dbcsr_get_data_size (matrix%wms(i)%data_area)
          WRITE(*,'(1X,A)')routineN//' row_i='
          IF (ASSOCIATED (matrix%wms(i)%row_i))&
               WRITE(*,'(10(1X,I7))')matrix%wms(i)%row_i(:)
          WRITE(*,'(1X,A)')routineN//' col_i='
          IF (ASSOCIATED (matrix%wms(i)%col_i))&
               WRITE(*,'(10(1X,I7))')matrix%wms(i)%col_i(:)
          WRITE(*,'(1X,A)')routineN//' blk_p='
          IF (ASSOCIATED (matrix%wms(i)%blk_p))&
               WRITE(*,'(10(1X,I7))')matrix%wms(i)%blk_p(:)
       ENDDO
    ENDIF
    matrix%row_p(:) = -1
    matrix%row_p(1) = 0
    blk_p = 1
    blk = 1
    row = 1
    prev_row = 1
    smallest_row = 1
    DO WHILE (.NOT. ALL (finished))
       ! Find the (still-valid) source block that has the smallest column in
       ! the smallest row. Start with any valid source block.
       smallest = MINLOC (wm_blk_p, 1, .NOT.finished)
       smallest_row = matrix%wms(smallest)%row_i(wm_blk_p(smallest))
       smallest_col = matrix%wms(smallest)%col_i(wm_blk_p(smallest))
       DO i = 1, nwms
          IF (.NOT. finished(i)) THEN
             IF (matrix%wms(i)%row_i(wm_blk_p(i)) .LE. smallest_row) THEN
                IF (matrix%wms(i)%row_i(wm_blk_p(i)) .LT. smallest_row) THEN
                   smallest_row = matrix%wms(i)%row_i(wm_blk_p(i))
                   smallest_col = matrix%wms(i)%col_i(wm_blk_p(i))
                   smallest = i
                ELSE
                   IF (matrix%wms(i)%col_i(wm_blk_p(i)) .LT. smallest_col) THEN
                      smallest_col = matrix%wms(i)%col_i(wm_blk_p(i))
                      smallest = i
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
       ENDDO
       matrix%col_i(blk) = smallest_col
       o_blk_p = matrix%wms(smallest)%blk_p(wm_blk_p(smallest))
       nze = rbs(smallest_row) * cbs(smallest_col)
       IF (dbg) THEN
          WRITE(*,'(1X,A,"(",I5,",",I5,") at ",I7,1X,A,I5,1X,I7,1X,A,1X,I7)')&
               routineN//" Adding",&
               smallest_row, smallest_col, blk_p, "from",smallest, o_blk_p,&
               "size", nze
       END IF
       IF (srt) THEN
          matrix%blk_p(blk) = SIGN (blk_p, o_blk_p)
          CALL dbcsr_data_set (matrix%data_area, ABS(blk_p), nze,&
               matrix%wms(smallest)%data_area, ABS(o_blk_p))
       ENDIF
       IF (smallest_row .GT. prev_row) THEN
          matrix%row_p(prev_row+1:smallest_row) = blk-1
          IF (dbg) WRITE(*,'(15(1X,I4))')matrix%row_p(:)
       ENDIF
       prev_row = smallest_row
       blk = blk + 1
       blk_p = blk_p + nze
       wm_blk_p(smallest) = wm_blk_p(smallest)+1
       finished(smallest) = wm_blk_p(smallest) .GT. matrix%wms(smallest)%lastblk
    ENDDO
    matrix%row_p(smallest_row+1:matrix%nblkrows_total+1) = blk-1
    matrix%nze = blk_p - 1
    IF (.NOT. srt .AND. ASSOCIATED (matrix%wms(1)%blk_p)) THEN
       matrix%blk_p(1:matrix%wms(1)%lastblk)&
            = matrix%wms(1)%blk_p(1:matrix%wms(1)%lastblk)
    ENDIF
    IF (dbg) WRITE(*,*)routineP//'row_p=',matrix%row_p(:)
    IF (dbg) WRITE(*,*)routineP//'col_i=',matrix%col_i(:)
    IF (dbg) WRITE(*,*)routineP//'blk_p=',matrix%blk_p(:)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_merge_data



! *****************************************************************************
!> \brief Calls the functions that are called in the included files but not
!>        in thin one. Needed to make prettify think the functions ary really
!>        used.
!> \param[in] matrix          data_area
!> \param[out] temp           pointer to the data
! *****************************************************************************
  SUBROUTINE prettify_fooler (data_area, temp)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: data_area
    REAL(KIND=dp), DIMENSION(:), POINTER     :: temp

    INTEGER                                  :: i
    INTEGER(KIND=int_8)                      :: key
    INTEGER, DIMENSION(:), POINTER           :: p_i
    REAL(KIND=real_4), DIMENSION(:), POINTER :: r4
    REAL(KIND=real_8), DIMENSION(:), POINTER :: r8
    TYPE(array_i1d_obj)                      :: array_i
    TYPE(btree_2d_data_c)                    :: bt_dc
    TYPE(btree_2d_data_d)                    :: bt_dd
    TYPE(btree_2d_data_s)                    :: bt_ds
    TYPE(btree_2d_data_z)                    :: bt_dz
    TYPE(btree_c)                            :: bt_c
    TYPE(btree_d)                            :: bt_d
    TYPE(btree_s)                            :: bt_s
    TYPE(btree_z)                            :: bt_z
    TYPE(dbcsr_error_type)                   :: error
    TYPE(dbcsr_obj)                          :: matrix

    IF (.TRUE.) THEN
       r4 => dbcsr_get_data_p (data_area, 0.0_real_4)
       r8 => dbcsr_get_data_p (data_area, 0.0_real_8)
       CALL addto_array (data_area, r8, error=error)
       CALL btree_add_s (bt_s, key, bt_ds)
       CALL btree_add_d (bt_d, key, bt_dd)
       CALL btree_add_c (bt_c, key, bt_dc)
       CALL btree_add_z (bt_z, key, bt_dz)
       CALL btree_get_s (bt_s, key, bt_ds)
       CALL btree_get_d (bt_d, key, bt_dd)
       CALL btree_get_c (bt_c, key, bt_dc)
       CALL btree_get_z (bt_z, key, bt_dz)
       p_i => array_data (array_i)
    ENDIF
    IF (.TRUE.) THEN
       temp => dbcsr_get_data_p (data_area, 0.0_dp)
       CALL dbcsr_get_info (matrix)
       CALL dbcsr_mutable_release(matrix%m%wms(0)%mutable)
       CALL dbcsr_unpack_i8_2i4 ( (/ key /), p_i, p_i)
    ENDIF
  END SUBROUTINE prettify_fooler

#include "dbcsr_work_operations_d.F"
#include "dbcsr_work_operations_z.F"
#include "dbcsr_work_operations_s.F"
#include "dbcsr_work_operations_c.F"


END MODULE dbcsr_work_operations

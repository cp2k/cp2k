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
  USE btree_I8_k_cp2d_v,               ONLY: btree_2d_data_c => cp2d
  USE btree_I8_k_dp2d_v,               ONLY: btree_2d_data_d => dp2d
  USE btree_I8_k_sp2d_v,               ONLY: btree_2d_data_s => sp2d
  USE btree_I8_k_zp2d_v,               ONLY: btree_2d_data_z => zp2d
  USE dbcsr_block_buffers,             ONLY: dbcsr_buffers_2d_needed,&
                                             dbcsr_buffers_flush,&
                                             dbcsr_buffers_init,&
                                             dbcsr_buffers_new,&
                                             dbcsr_buffers_release
  USE dbcsr_block_operations,          ONLY: dbcsr_data_set
  USE dbcsr_data_methods,              ONLY: &
       dbcsr_data_clear_pointer, dbcsr_data_ensure_size, dbcsr_data_get_size, &
       dbcsr_data_get_size_referenced, dbcsr_data_get_type, dbcsr_data_hold, &
       dbcsr_data_init, dbcsr_data_new, dbcsr_data_release
  USE dbcsr_data_operations,           ONLY: dbcsr_data_copyall,&
                                             dbcsr_sort_data
  USE dbcsr_error_handling
  USE dbcsr_index_operations,          ONLY: dbcsr_addto_index_array,&
                                             dbcsr_build_row_index,&
                                             dbcsr_clearfrom_index_array,&
                                             dbcsr_count_row_index,&
                                             dbcsr_make_dbcsr_index,&
                                             dbcsr_repoint_index,&
                                             dbcsr_sort_indices,&
                                             sort
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_kinds,                     ONLY: default_string_length,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE dbcsr_message_passing,           ONLY: mp_allocate
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_sizes, dbcsr_destroy_array, dbcsr_destroy_image_dist, &
       dbcsr_distribution, dbcsr_distribution_has_threads, &
       dbcsr_distribution_hold, dbcsr_distribution_make_threads, &
       dbcsr_distribution_release, dbcsr_get_data_type, &
       dbcsr_get_matrix_type, dbcsr_get_replication_type, dbcsr_init, &
       dbcsr_is_initialized, dbcsr_modify_lock, dbcsr_modify_unlock, &
       dbcsr_mutable_destroy, dbcsr_mutable_init, dbcsr_mutable_instantiated, &
       dbcsr_mutable_new, dbcsr_mutable_release, dbcsr_name, &
       dbcsr_row_block_sizes, dbcsr_switch_data_area, dbcsr_use_mutable, &
       dbcsr_uses_special_memory, dbcsr_valid_index, dbcsr_wm_use_mutable
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size
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
                                             dbcsr_calc_block_sizes,&
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

  INTERFACE dbcsr_create
     MODULE PROCEDURE dbcsr_create_new, dbcsr_create_template
  END INTERFACE

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
  SUBROUTINE dbcsr_create_new(matrix, name, dist, matrix_type,&
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

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_new', &
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
    ELSEIF (dbcsr_buffers_2d_needed) THEN
       CALL dbcsr_buffers_flush (matrix%m%buffers, error=error)
       CALL dbcsr_buffers_release (matrix%m%buffers, error=error)
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
!$  IF (.NOT. dbcsr_distribution_has_threads (matrix%m%dist)) THEN
!$     CALL dbcsr_distribution_make_threads (matrix%m%dist,&
!$          array_data(matrix%m%row_blk_size))
!$  ENDIF
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
    CALL dbcsr_buffers_init (matrix%m%buffers)
    IF (dbcsr_buffers_2d_needed) THEN
       CALL dbcsr_buffers_new (matrix%m%buffers, matrix%m%data_area,&
            error=error)
    ENDIF
!$  CALL OMP_INIT_LOCK (matrix%m%modification_lock)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_create_new

  SUBROUTINE dbcsr_create_template(matrix, template, name, dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, data_type, special, reuse,&
       mutable_work, replication_type, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_obj), INTENT(IN)              :: template
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: name
    TYPE(dbcsr_distribution_obj), &
      INTENT(IN), OPTIONAL                   :: dist
    CHARACTER, INTENT(IN), OPTIONAL          :: matrix_type
    TYPE(array_i1d_obj), INTENT(IN), &
      OPTIONAL                               :: row_blk_size, col_blk_size
    INTEGER, INTENT(IN), OPTIONAL            :: nblks, nze
    CHARACTER, INTENT(IN), OPTIONAL          :: data_type
    LOGICAL, INTENT(IN), OPTIONAL            :: special, reuse, mutable_work
    CHARACTER, INTENT(IN), OPTIONAL          :: replication_type
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_template', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: new_data_type, &
                                                new_matrix_type, &
                                                new_replication_type
    CHARACTER(len=default_string_length)     :: new_name
    LOGICAL                                  :: new_mutable_work, new_special
    TYPE(array_i1d_obj)                      :: new_col_blk_size, &
                                                new_row_blk_size
    TYPE(dbcsr_distribution_obj)             :: new_dist

!   ---------------------------------------------------------------------------

    IF (PRESENT (name)) THEN
       new_name = TRIM(name)
    ELSE
       new_name = TRIM(dbcsr_name (template))
    ENDIF
    IF (PRESENT (dist)) THEN
       new_dist = dist
    ELSE
       new_dist = dbcsr_distribution (template)
    ENDIF
    IF (PRESENT (matrix_type)) THEN
       new_matrix_type = matrix_type
    ELSE
       new_matrix_type = dbcsr_get_matrix_type (template)
    ENDIF
    IF (PRESENT (row_blk_size)) THEN
       new_row_blk_size = row_blk_size
    ELSE
       new_row_blk_size = dbcsr_row_block_sizes (template)
    ENDIF
    IF (PRESENT (col_blk_size)) THEN
       new_col_blk_size = col_blk_size
    ELSE
       new_col_blk_size = dbcsr_col_block_sizes (template)
    ENDIF
    IF (PRESENT (data_type)) THEN
       new_data_type = data_type
    ELSE
       new_data_type = dbcsr_get_data_type (template)
    ENDIF
    IF (PRESENT (special)) THEN
       new_special = special
    ELSE
       new_special = dbcsr_uses_special_memory (template)
    ENDIF
    IF (PRESENT (replication_type)) THEN
       new_replication_type = replication_type
    ELSE
       new_replication_type = dbcsr_get_replication_type (template)
    ENDIF
    IF (PRESENT (mutable_work)) THEN
       new_mutable_work = mutable_work
    ELSE
       new_mutable_work = dbcsr_use_mutable (template%m)
    ENDIF
    CALL dbcsr_create (matrix, name=new_name, dist=new_dist,&
         matrix_type = new_matrix_type,&
         row_blk_size = new_row_blk_size,&
         col_blk_size = new_col_blk_size,&
         nblks = nblks, nze = nze,&
         data_type = new_data_type,&
         special = new_special,&
         mutable_work = new_mutable_work,&
         replication_type = new_replication_type,&
         error=error)
  END SUBROUTINE dbcsr_create_template


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
!$OMP FLUSH
!$  CALL dbcsr_modify_unlock (matrix)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_work_create



! *****************************************************************************
!> \brief Creates the final dbcsr_type matrix from the working matrix.
!>
!> Work matrices (array or tree-based) are merged into the base DBCSR matrix.
!>
!> If a matrix has a valid index, then nothing is done.
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
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handler, nwms, &
                                                start_offset
    INTEGER, ALLOCATABLE, DIMENSION(:), &
      SAVE, TARGET                           :: empty_row_p
    INTEGER, DIMENSION(:), POINTER, SAVE     :: old_blk_p, old_col_i, &
                                                old_row_p
    LOGICAL                                  :: fake_row_p, sort_data, spawn

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    !$OMP BARRIER
    CALL dbcsr_assert (dbcsr_is_initialized (matrix),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Can not finalize uninitialized matrix.",__LINE__,error)
    IF (dbcsr_buffers_2d_needed) THEN
       CALL dbcsr_buffers_flush (matrix%m%buffers, error=error)
    ENDIF
    ! Remove the cached pre-distributed images.
    !$OMP SINGLE
    IF (ASSOCIATED (matrix%m%predistributed)) THEN
          CALL dbcsr_destroy_array (matrix%m%predistributed,error)
          CALL dbcsr_destroy_image_dist (matrix%m%predistributed%image_dist)
          DEALLOCATE (matrix%m%predistributed)
    ENDIF
    NULLIFY (matrix%m%predistributed)
    !$OMP END SINGLE
    ! If the matrix is not marked as dirty then skip the work.
    IF (dbcsr_valid_index(matrix)) THEN
       !"No need to finalize a valid matrix, skipping."
       !
       ! A matrix with a valid index should not have associated work
       ! arrays.  This may happen when this routine is called on a
       ! matrix that was not changed.
       !$OMP BARRIER
       !$OMP SINGLE
       IF (ASSOCIATED (matrix%m%wms)) &
            CALL dbcsr_work_destroy_all(matrix%m)
       matrix%m%valid = .TRUE.
       !$OMP END SINGLE
       CALL dbcsr_error_stop(error_handler, error)
       RETURN
    ENDIF
    !
    ! If possible, data copying is avoided.
    IF (PRESENT (reshuffle)) THEN
       sort_data = reshuffle
    ELSE
       sort_data = .FALSE.
    ENDIF
    !
    ! Now make sure that a valid row_p exists. Also clear the row_p if
    ! the matrix is declared to have 0 blocks.
    !$OMP MASTER
    fake_row_p = .NOT. ASSOCIATED (matrix%m%row_p)
    IF (ASSOCIATED (matrix%m%row_p)) THEN
       fake_row_p = SIZE (matrix%m%row_p) .LE. 1
    ENDIF
    fake_row_p = fake_row_p .OR. matrix%m%nblks .EQ. 0
    !$OMP END MASTER
    !
    ! See where data will be appended in the main data
    ! area. Alternatively, set to the start if the matrix is declared
    ! to have no data. (This value is ignored if reshuffle is true
    ! because the main data area is always new.)
    start_offset = dbcsr_data_get_size_referenced (matrix%m%data_area)+1
    IF (matrix%m%nze .EQ. 0) start_offset=1
    !$OMP MASTER
    matrix%m%index(dbcsr_slot_nze) = matrix%m%nze
    IF (fake_row_p) THEN
       ALLOCATE (empty_row_p (matrix%m%nblkrows_total+1))
       empty_row_p(:) = 0
       CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_row_p,&
            DATA=empty_row_p, extra=0, error=error)
       CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_col_i,&
            reservation=0, error=error)
       CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_blk_p,&
            reservation=0, error=error)
       CALL dbcsr_repoint_index (matrix%m)
    ENDIF
    old_row_p => matrix%m%row_p
    old_col_i => matrix%m%col_i
    old_blk_p => matrix%m%blk_p
    !$OMP END MASTER
!$OMP BARRIER
    !
    ! If the matrix, work matricies, and environment fit several
    ! criteria, then a quick O(1) finalization is performed.
    IF (can_quickly_finalize (matrix) .AND. .NOT. sort_data) THEN
       !$OMP BARRIER
       CALL quick_finalize (matrix, error)
    ELSEIF (ASSOCIATED (matrix%m%wms)) THEN
       ! Check to see if we will need to create a parallel environment
       ! (needed when there are multiple work matrices but we are not
       ! in an OpenMP parallel section.)
       nwms = SIZE (matrix%m%wms)
       spawn = .FALSE.
       !$ IF (.NOT. OMP_IN_PARALLEL ()) THEN
       !$    IF (nwms .GT. 1) spawn = .TRUE.
       !$ ENDIF
       !$ IF (.NOT. spawn) THEN
          CALL dbcsr_merge_all (matrix%m,&
               old_row_p, old_col_i, old_blk_p,&
               data_starting_offset = start_offset,&
               sort_data=sort_data, &
               error=error)
       !$ ELSE
          !
          ! This here is used when the matrix has more work
          ! matrices. It's a shortcut when the finalize is called from
          ! a non-parallel environment whereas the matrix was
          ! built/modified in a parallel environment
          !
          !$OMP PARALLEL DEFAULT (none) &
          !$OMP          SHARED (matrix, old_row_p, old_col_i, old_blk_p,&
          !$OMP                  start_offset, sort_data, error)
          !$   CALL dbcsr_merge_all (matrix%m,&
          !$        old_row_p, old_col_i, old_blk_p,&
          !$        data_starting_offset &
          !$        = start_offset,&
          !$        sort_data=sort_data, &
          !$        error=error)
          !$
          !$OMP END PARALLEL
       !$ END IF
    ENDIF
!$OMP BARRIER
!$OMP SINGLE
    ! Clean up.
    IF (ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_work_destroy_all(matrix%m)
    ENDIF
    matrix%m%valid = .TRUE.
!$OMP END SINGLE
    IF (dbg) THEN
       !$OMP SINGLE
       CALL dbcsr_verify_matrix (matrix, error=error)
       !$OMP END SINGLE
    ENDIF
!$OMP MASTER
    IF (fake_row_p) THEN
       DEALLOCATE (empty_row_p)
    ENDIF
!$OMP END MASTER
    IF (dbcsr_buffers_2d_needed) THEN
       CALL dbcsr_buffers_new (matrix%m%buffers, matrix%m%data_area,&
            error=error)
    ENDIF
!$OMP BARRIER
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_finalize

! *****************************************************************************
!> \brief Checks whether the matrix can be finalized with minimal copying.
!> \param[in] matrix          matrix to check
!> \result quick              whether matrix can be quickly finalized
! *****************************************************************************
  FUNCTION can_quickly_finalize (matrix) RESULT (quick)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL                                  :: quick

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (matrix%m%wms)) THEN
       quick = matrix%m%nblks .EQ. 0
       quick = quick .AND. SIZE(matrix%m%wms) .EQ. 1 .AND.&
            .NOT. dbcsr_wm_use_mutable(matrix%m%wms(1))
       IF (quick) THEN
          quick = quick .AND.&
               dbcsr_uses_special_memory(matrix%m%wms(1)%data_area) .EQV. &
               dbcsr_uses_special_memory(matrix%m%data_area)
       ENDIF
    ELSE
       quick = .FALSE.
    ENDIF
  END FUNCTION can_quickly_finalize

! *****************************************************************************
!> \brief Performs quick finalization of matrix
!>
!> The data area from the work matrix is accepted as the new matrix's data
!> area and the index is built from the work matrix.
!> \param[in,out] matrix      matrix to finalize
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE quick_finalize (matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'quick_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, nblks, nrows

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handle, error)
    !$OMP SECTIONS
    !$OMP SECTION
    nblks = matrix%m%wms(1)%lastblk
    nrows = matrix%m%nblkrows_total
    CALL dbcsr_sort_indices (nblks,&
         matrix%m%wms(1)%row_i,&
         matrix%m%wms(1)%col_i,&
         matrix%m%wms(1)%blk_p)
    CALL dbcsr_clearfrom_index_array (matrix%m, dbcsr_slot_col_i)
    CALL dbcsr_clearfrom_index_array (matrix%m, dbcsr_slot_blk_p)
    CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_row_p,&
         reservation=nrows+1, extra=2*nblks, error=error)
    CALL dbcsr_make_dbcsr_index (matrix%m%row_p, matrix%m%wms(1)%row_i,&
         nrows, nblks)
    CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_col_i,&
         DATA=matrix%m%wms(1)%col_i(1:nblks), error=error)
    CALL dbcsr_addto_index_array (matrix%m, dbcsr_slot_blk_p,&
         DATA=matrix%m%wms(1)%blk_p(1:nblks), error=error)
    matrix%m%nblks = nblks
    matrix%m%nze = matrix%m%wms(1)%datasize
    matrix%m%index(dbcsr_slot_nblks) = nblks
    matrix%m%index(dbcsr_slot_nze) = matrix%m%wms(1)%datasize
    CALL dbcsr_repoint_index (matrix%m)
    !$OMP SECTION
    CALL dbcsr_switch_data_area (matrix, matrix%m%wms(1)%data_area, error=error)
    !$OMP END SECTIONS
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE quick_finalize



! *****************************************************************************
!> \brief Creates a work matrix from the data present in a finalized matrix.
!> \param[in,out] matrix      DBCSR matrix
!> \param[in] limits          (optional) the limits to use for copying
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE dbcsr_add_wm_from_matrix(matrix, limits, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, DIMENSION(4), INTENT(IN), &
      OPTIONAL                               :: limits
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_wm_from_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, ithread, &
                                                nthreads, nwms, old_nwms
    LOGICAL                                  :: preexists

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    !$OMP BARRIER
    preexists = ASSOCIATED (matrix%m%wms)
    IF (preexists) THEN
       old_nwms = SIZE (matrix%m%wms)
       CALL dbcsr_assert (old_nwms, "NE", 0, dbcsr_warning_level,&
            dbcsr_internal_error, routineN, "Nonexisting work matrices?!",&
            __LINE__, error=error)
    ELSE
       old_nwms = 0
    ENDIF
    nthreads = 1; ithread = 0
    !$ nthreads = OMP_GET_NUM_THREADS() ; ithread = OMP_GET_THREAD_NUM()
    IF (nthreads .GT. 1) THEN
       CALL dbcsr_assert (old_nwms .EQ. nthreads, "OR", old_nwms .EQ. 0,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Number of work matrices and threads do not match",&
            __LINE__, error=error)
    ENDIF
    nwms = MAX(1, old_nwms)
    !$ nwms = MAX(nwms, nthreads)
    !$OMP BARRIER
    CALL dbcsr_work_create (matrix, matrix%m%nblks, matrix%m%nze, n=nwms,&
         work_mutable=.FALSE., error=error)
    !$OMP BARRIER
    CALL dbcsr_fill_wm_from_matrix (matrix%m%wms, matrix%m,&
         limits=limits, error=error)
    !$OMP BARRIER
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_add_wm_from_matrix


! *****************************************************************************
!> \brief Fills index and data of the work matrix from the
!>        previously-finalized one.
!> \param[in,out] wm          the work matrix to fill
!> \param[in,out] matrix      DBCSR matrix
!> \param[in] limits          (optional) only fills blocks within this range
!> \param[in,out] error       error
!> \par limits
!> The limits is a 4-tuple
!> (lower_row, higher_row, lower_column, higher_column).
! *****************************************************************************
  SUBROUTINE dbcsr_fill_wm_from_matrix(wm, matrix, limits, error)
    TYPE(dbcsr_work_type), DIMENSION(:), &
      INTENT(INOUT)                          :: wm
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, DIMENSION(4), INTENT(IN), &
      OPTIONAL                               :: limits
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_fill_wm_from_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, blk_p, col, error_handler, ithread, nthreads, nwms, nze, &
      row, wblk_p, which_wm, wm_first, wm_last
    LOGICAL                                  :: limit, mt, tr
    TYPE(dbcsr_data_obj)                     :: data_block
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_obj)                          :: m

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    nwms = SIZE (matrix%wms)
    m%m = matrix
    mt = .FALSE.
    !$ IF (nwms .GT. 1) mt = omp_get_num_threads() .GT. 1
    ithread = 0 ; nthreads = 1
    !$ ithread = omp_get_thread_num ()
    !$ nthreads = omp_get_num_threads()
    limit = PRESENT (limits)
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_new (data_block, dbcsr_data_get_type (matrix%data_area))
    IF (mt) THEN
       wm_first = ithread+1
       wm_last = ithread+1
    ELSE
       wm_first = 1
       wm_last = nwms
    ENDIF
    ! Prepares the work matrices to accept the main data.
    DO which_wm = wm_first, wm_last
       CALL dbcsr_assert ("NOT", dbcsr_wm_use_mutable (wm(which_wm)),&
            dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
            "Adding main matrix into mutable not supported.", __LINE__,&
            error=error)
       IF (dbcsr_wm_use_mutable (wm(which_wm))) THEN
          IF (.NOT. dbcsr_mutable_instantiated (wm(which_wm)%mutable)) THEN
             CALL dbcsr_mutable_new (wm(which_wm)%mutable, matrix%data_type)
          ENDIF
       ELSE
          ! We don't know how much data we'll get so we have to be generous.
          !$OMP CRITICAL
          CALL dbcsr_data_ensure_size (wm(which_wm)%data_area,&
               dbcsr_data_get_size (matrix%data_area), error=error)
          !$OMP END CRITICAL
       ENDIF
    ENDDO
    ! Now copy the data
    CALL dbcsr_iterator_start (iter, m, shared=mt)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, data_block,&
            transposed=tr, block_number=blk)
       IF (limit) THEN
          IF (.NOT. within_limits (row, col, limits)) CYCLE
       ENDIF
       blk_p = matrix%blk_p(blk)
       which_wm = ithread+1
       wblk_p = SIGN (wm(which_wm)%datasize+1, blk_p)
       nze = dbcsr_data_get_size (data_block)
       IF (mt .OR. limit) THEN
          ! The data gets copied block by block so the block pointers
          ! are ordered accordingly.
          CALL add_work_coordinate (wm(which_wm), row, col, wblk_p, error=error)
          CALL dbcsr_data_set (wm(which_wm)%data_area,&
               lb=ABS(wblk_p),&
               data_size=nze,&
               src=data_block, source_lb=1)
       ELSE
          ! The data gets copied all at once so the block pointers
          ! should remain the same as they were.
          CALL add_work_coordinate (wm(which_wm), row, col, blk_p, error=error)
       ENDIF
       wm(which_wm)%datasize = wm(which_wm)%datasize+nze
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    CALL dbcsr_data_clear_pointer (data_block)
    CALL dbcsr_data_release (data_block)
    ! Copy all blocks at once
    IF (.NOT. mt .AND. .NOT. limit) THEN
       DO which_wm = 1, nwms
          CALL dbcsr_data_ensure_size (wm(which_wm)%data_area,&
               dbcsr_data_get_size_referenced (matrix%data_area), error=error)
          CALL dbcsr_data_copyall (wm(which_wm)%data_area, matrix%data_area,&
               error=error)
          wm(which_wm)%datasize = dbcsr_data_get_size_referenced (wm(which_wm)%data_area)
       ENDDO
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_fill_wm_from_matrix

! *****************************************************************************
!> \brief Checks whether a point is within bounds
!> \param[in] row, column    point to check
!> \param[in] limits         limits (low_row, high_row, low_col, high_col)
!> \result                   whether the point is within the bounds
! *****************************************************************************

  PURE FUNCTION within_limits (row, column, limits)
    INTEGER, INTENT(IN)                      :: row, column
    INTEGER, DIMENSION(4), INTENT(IN)        :: limits
    LOGICAL                                  :: within_limits

    within_limits =  row .GE. limits(1) .AND. row .LE. limits(2) .AND.&
         column .GE. limits(3) .AND. column .LE. limits(4)
  END FUNCTION within_limits


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
!> \brief Count actual size taken by work matrix
!> \param[in] row_i           row indices
!> \param[in] col_i           column indices
!> \param[in] row_block_sizes  offsets of rows
!> \param[in] col_block_sizes  offsets of columns
!> \param[in] nblks              number of blocks
!> \param[out] nze               counted data size
! *****************************************************************************
  PURE SUBROUTINE dbcsr_count_wm(row_i, col_i,&
       row_block_sizes, col_block_sizes, nblks, nze)
    INTEGER, INTENT(IN)                      :: nblks
    INTEGER, DIMENSION(:), INTENT(IN)        :: col_block_sizes, &
                                                row_block_sizes
    INTEGER, DIMENSION(nblks), INTENT(IN)    :: col_i, row_i
    INTEGER, INTENT(OUT)                     :: nze

    INTEGER                                  :: blk

!   ---------------------------------------------------------------------------

    nze = 0
    DO blk = 1, nblks
       nze = nze + row_block_sizes(row_i(blk)) * col_block_sizes(col_i(blk))
    ENDDO
  END SUBROUTINE dbcsr_count_wm


! *****************************************************************************
!> \brief Merge data from matrix and work matrices into the final matrix.
!>
!> \param[in,out] matrix      matrix to work on
!> \param[in] data_starting_offset   Where to add data
!> \param[in] sort_data       whether data will be fully sorted
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE dbcsr_merge_all (matrix, old_row_p, old_col_i, old_blk_p,&
       data_starting_offset, sort_data, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    INTEGER, DIMENSION(:), POINTER           :: old_row_p, old_col_i, &
                                                old_blk_p
    INTEGER, INTENT(IN)                      :: data_starting_offset
    LOGICAL, INTENT(IN)                      :: sort_data
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_merge_all', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: index_blocks = 1, &
                                                index_data = 2
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER :: data_current_offset, data_to_add = 0, error_handler, i, &
      ithread, my_num_new_blocks, nblks_to_add = 0, new_nblks = 0, nrows, &
      nthreads, nwms, which_wm, wm_first, wm_last
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: additions_per_row
    INTEGER, DIMENSION(2)                    :: additions_total
    INTEGER, DIMENSION(:), POINTER, SAVE :: mt_data_offsets, new_blk_d, &
      new_blk_p, new_blk_sizes, new_col_i, new_row_p, row_size_counter
    INTEGER, DIMENSION(:, :), POINTER        :: my_blk_p
    LOGICAL                                  :: mt
    TYPE(dbcsr_data_obj), DIMENSION(:), &
      POINTER, SAVE                          :: all_data_areas

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
!$OMP BARRIER
    nwms = SIZE(matrix%wms)
!$OMP MASTER
    nblks_to_add = 0
    data_to_add = 0
    !
    ! If the data will be sorted then the source data areas are all
    ! put into an array. The new_blk_d elements select from which data
    ! area a block will be copied.
    NULLIFY (all_data_areas)
    IF (sort_data) THEN
       ALLOCATE (all_data_areas(1:nwms+1))
       DO i = 1, nwms+1
          CALL dbcsr_data_init (all_data_areas(i))
       ENDDO
       all_data_areas(1) = matrix%data_area
       CALL dbcsr_data_hold (all_data_areas(1))
    ENDIF
!$OMP END MASTER
!$OMP BARRIER
    ithread = 0 ; nthreads = 1
    !$ ithread = OMP_GET_THREAD_NUM() ; nthreads = OMP_GET_NUM_THREADS()
    mt = .FALSE.
    !$ mt = nthreads .GT. 1
    !$ IF (mt) THEN
    !$    CALL dbcsr_assert (nwms, "EQ", nthreads, dbcsr_fatal_level,&
    !$         dbcsr_caller_error, routineN, &
    !$         "Number of work matrices does not match number of threads.",&
    !$         __LINE__, error=error)
    !$ ENDIF
    ! One thread has to handle many work matrices
    nrows = matrix%nblkrows_total
    !
    ! The first step is to count how many blocks will be added per row.
    data_current_offset = data_starting_offset
    IF (mt) THEN
       wm_first = ithread + 1
       wm_last = ithread + 1
    ELSE
       wm_first = 1
       wm_last = nwms
       CALL dbcsr_assert (nwms, "EQ", 1,&
            dbcsr_warning_level, dbcsr_caller_error, routineN,&
            "Finalize should be called from parallel.",&
            __LINE__, error=error)
    ENDIF
    DO which_wm = wm_first, wm_last
       old_row_p => matrix%row_p
       old_col_i => matrix%col_i
       old_blk_p => matrix%blk_p
       ALLOCATE (additions_per_row(2, nrows))
       additions_per_row(:,:) = 0
       my_num_new_blocks = matrix%wms(which_wm)%lastblk
       ALLOCATE (my_blk_p(2, my_num_new_blocks))
       CALL count_blocks_to_add (matrix%wms(which_wm), additions_per_row,&
            additions_total, my_blk_p,&
            array_data (matrix%row_blk_size), array_data (matrix%col_blk_size),&
            error=error)
       !$OMP ATOMIC
       nblks_to_add = nblks_to_add + additions_total(index_blocks)
       !$OMP ATOMIC
       data_to_add = data_to_add + additions_total(index_data)
!$OMP MASTER
       ALLOCATE (mt_data_offsets (nwms))
       mt_data_offsets(1) = data_starting_offset
       ALLOCATE (new_row_p (nrows+1))
       ! Counts new number of blocks per row (to be later made into
       ! the new row_p index).
       CALL dbcsr_count_row_index (old_row_p, new_row_p, nrows)
!$OMP END MASTER
!$OMP BARRIER
       !
       ! Now add the blocks to be added to the row counts
       CALL add_to_row_count (new_row_p, additions_per_row(index_blocks,:))
       IF (sort_data) THEN
          ! When sort_data is true, then the following arry holds the
          ! source data array, not the true data offset.
          mt_data_offsets(which_wm) = which_wm + 1
       ELSE
          ! Calculate the data offsets at which the work data is
          ! copied into the main data area.
          CALL tally_offsets (mt_data_offsets, additions_total(index_data),&
               data_starting_offset)
       ENDIF
!$OMP BARRIER
!$OMP MASTER
       !!$OMP SECTION
       matrix%nblks = old_row_p(nrows+1)
       new_nblks = matrix%nblks + nblks_to_add
       ALLOCATE (new_col_i (new_nblks))
       ALLOCATE (new_blk_p (new_nblks))
       IF (sort_data) THEN
          ALLOCATE (new_blk_d(new_nblks))
          ALLOCATE (new_blk_sizes(new_nblks))
       ELSE
          NULLIFY (new_blk_d)
       ENDIF
       !!$OMP SECTION
       CALL dbcsr_build_row_index (new_row_p, nrows)
       !!$OMP SECTION
       ALLOCATE (row_size_counter (nrows))
       row_size_counter(:) = 0
       !!$OMP SECTION
       IF (.NOT. sort_data) THEN
          ! Resize the data area to accept the new data to be appended.
          CALL dbcsr_data_ensure_size (matrix%data_area, &
               data_size=data_starting_offset+data_to_add-1, error=error)
       ELSE
          ! Makes a new data area for the main matrix (previous data
          ! is already saved in all_data_areas(1)).
          CALL dbcsr_data_release (matrix%data_area)
          CALL dbcsr_data_init (matrix%data_area)
          CALL dbcsr_data_new (matrix%data_area,&
               dbcsr_data_get_type(all_data_areas(1)),&
               data_size = data_starting_offset+data_to_add-1,&
               special_memory = dbcsr_uses_special_memory(all_data_areas(1)))
       ENDIF
!$OMP END MASTER
!$OMP BARRIER
       ! Inserts new blocks into the index.
       CALL index_insertion (old_row_p, old_col_i, old_blk_p,&
            new_row_p, new_col_i, new_blk_p, new_blk_d, &
            additions_per_row(index_blocks,:),&
            matrix%wms(which_wm)%col_i, matrix%wms(which_wm)%blk_p, my_blk_p,&
            mt_data_offsets(which_wm),&
            row_size_counter, error=error)
!$OMP BARRIER
       IF (sort_data) THEN
          all_data_areas(mt_data_offsets(which_wm)) &
               = matrix%wms(which_wm)%data_area
          !$OMP CRITICAL
          CALL dbcsr_data_hold (all_data_areas(mt_data_offsets(which_wm)))
          !$OMP END CRITICAL
       ENDIF
       !$OMP MASTER
       ! Creates new index array.
       CALL dbcsr_clearfrom_index_array (matrix, dbcsr_slot_col_i)
       CALL dbcsr_clearfrom_index_array (matrix, dbcsr_slot_blk_p)
       CALL dbcsr_addto_index_array (matrix, dbcsr_slot_row_p,&
            DATA=new_row_p(1:nrows+1), error=error)
       CALL dbcsr_addto_index_array (matrix, dbcsr_slot_col_i,&
            reservation=new_nblks, extra=SIZE(new_blk_p), error=error)
       CALL dbcsr_addto_index_array (matrix, dbcsr_slot_blk_p,&
            reservation=new_nblks, error=error)
       matrix%nblks = new_nblks
       matrix%nze = matrix%nze + data_to_add
       matrix%index(dbcsr_slot_nblks) = matrix%nblks
       matrix%index(dbcsr_slot_nze) = matrix%nze
       CALL dbcsr_repoint_index(matrix)
       !$OMP END MASTER
!$OMP BARRIER
       IF (.NOT. sort_data) THEN
          ! Copies data into the main data area.
          CALL merge_my_wm (matrix%data_area,&
               matrix%wms(which_wm), my_blk_p, error=error)
          !$OMP WORKSHARE
          matrix%col_i(1:new_nblks) = new_col_i(1:new_nblks)
          matrix%blk_p(1:new_nblks) = new_blk_p(1:new_nblks)
          !$OMP END WORKSHARE
       ELSE
          CALL dbcsr_calc_block_sizes (new_blk_sizes, new_row_p, new_col_i,&
               array_data (matrix%row_blk_size),&
               array_data (matrix%col_blk_size))
          ! Merges all data into the main data area.
          CALL dbcsr_sort_data (matrix%blk_p, new_blk_p, new_blk_sizes,&
               matrix%data_area, all_data_areas(1), all_data_areas,&
               new_blk_d, error=error)
          !$OMP WORKSHARE
          matrix%col_i(1:new_nblks) = new_col_i(1:new_nblks)
          !$OMP END WORKSHARE
          !$OMP MASTER
          DEALLOCATE (new_blk_d)
          DEALLOCATE (new_blk_sizes)
          !$OMP END MASTER
       ENDIF
       !$OMP BARRIER
       !$OMP MASTER
       DEALLOCATE (row_size_counter)
       DEALLOCATE (mt_data_offsets)
       DEALLOCATE (new_row_p)
       DEALLOCATE (new_col_i)
       DEALLOCATE (new_blk_p)
       CALL dbcsr_repoint_index (matrix)
       !$OMP END MASTER
       DEALLOCATE (my_blk_p)
       DEALLOCATE (additions_per_row)
    ENDDO
    !$OMP MASTER
    IF (sort_data) THEN
       DO i = 1, nwms+1
          CALL dbcsr_data_release (all_data_areas(i))
       ENDDO
       DEALLOCATE (all_data_areas)
    ENDIF
    !$OMP END MASTER
    !$OMP BARRIER
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_merge_all


! *****************************************************************************
!> \brief Merges the indices of the work and actual matrix.
!> \param[in] old_row_p       old row pointers
!> \param[in] old_col_i       old column numbers
!> \param[in] old_blk_p       old block pointers
!> \param[in] new_row_p       new row pointers
!> \param[out] new_col_i      new column numbers
!> \param[out] new_blk_p      new block pointers
!> \param[out] new_blk_d      which data area a block resides
!> \param[in] additions_per_row   the number of blocks added by this thread
!>                                for each row
!> \param[in] my_col_i        column number for this thread's blocks
!> \param[in] my_blk_meta     this thread's block meta (for transposes)
!> \param[out] my_blk_info    new block pointers for this thread's blocks
!> \param[in] data_offset     data offset for my thread
!> \param[in,out] row_size_counter    current offsets of blocks in row_p
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE index_insertion (old_row_p, old_col_i, old_blk_p,&
       new_row_p, new_col_i, new_blk_p, new_blk_d,&
       additions_per_row, my_col_i, my_blk_meta, my_blk_info, &
       data_offset, row_size_counter,&
       error)
    INTEGER, DIMENSION(:), INTENT(IN)        :: old_row_p, old_col_i, &
                                                old_blk_p, new_row_p
    INTEGER, DIMENSION(:), POINTER           :: new_col_i, new_blk_p, &
                                                new_blk_d
    INTEGER, DIMENSION(:), INTENT(IN)        :: additions_per_row, my_col_i, &
                                                my_blk_meta
    INTEGER, DIMENSION(:, :), INTENT(INOUT)  :: my_blk_info
    INTEGER, INTENT(IN)                      :: data_offset
    INTEGER, DIMENSION(:), POINTER           :: row_size_counter
    TYPE(dbcsr_error_type)                   :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'index_insertion', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, curr_blk_p, error_handle, i, ithread, my_block_counter, &
      nblks, ncols, new_first_blk, new_last_blk, new_nblks, nrows, nthreads, &
      nze, old_nblks, row, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: buff, buff2, permutation
    LOGICAL                                  :: copy_data

!
!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    copy_data = .NOT. ASSOCIATED (new_blk_d)
    nrows = SIZE(old_row_p)-1
    ithread = 0 ; nthreads = 1
    !$ ithread = OMP_GET_THREAD_NUM() ; nthreads = OMP_GET_NUM_THREADS()
    my_block_counter = 1
    IF (copy_data) curr_blk_p = data_offset
    DO row = 1, nrows
       IF (additions_per_row(row) .GT. 0) THEN
          new_nblks = additions_per_row(row)
          ! Find where in the col_i and blk_p arrays to add this
          ! thread's blocks.
          !$OMP CRITICAL
          blk = row_size_counter(row)
          !$OMP ATOMIC
          row_size_counter(row) = row_size_counter(row) + new_nblks
          !$OMP END CRITICAL
          old_nblks = old_row_p(row+1) - old_row_p(row)
          new_first_blk = new_row_p(row)+1 + old_nblks + blk
          new_last_blk = new_first_blk + new_nblks - 1
          ! Now add the new blocks to the global column index
          new_col_i(new_first_blk : new_last_blk)&
               = my_col_i(my_block_counter:my_block_counter+new_nblks-1)
          ! and calculate the new block pointers from the running offset and the
          ! stored sizes
          IF (copy_data) THEN
             DO blk = 1, new_nblks
                nze = my_blk_info(2, my_block_counter+blk-1)
                my_blk_info(1, my_block_counter+blk-1)&
                     = SIGN(curr_blk_p, my_blk_meta(my_block_counter+blk-1))
                curr_blk_p = curr_blk_p + nze
             ENDDO
          ELSE
             new_blk_d(new_first_blk : new_last_blk) = data_offset
          ENDIF
          new_blk_p(new_first_blk : new_last_blk)&
               = my_blk_info(1, my_block_counter:my_block_counter+new_nblks-1)
          my_block_counter = my_block_counter + new_nblks
       ENDIF
    ENDDO
    ! Copy existing blocks and sort everything
    nblks = SIZE (new_col_i)
    ALLOCATE (permutation(nblks), stat=stat)
    ALLOCATE (buff(nblks), stat=stat)
    IF (.NOT. copy_data) ALLOCATE (buff2(nblks), stat=stat)
    !$OMP BARRIER
    !$OMP DO
    DO row = 1, nrows
       old_nblks = old_row_p(row+1) - old_row_p(row)
       new_first_blk = new_row_p(row)+1
       new_last_blk = new_row_p(row+1)
       new_col_i(new_first_blk : new_first_blk + old_nblks - 1) &
            = old_col_i (old_row_p(row)+1 : old_row_p(row+1))
       new_blk_p(new_first_blk : new_first_blk + old_nblks - 1) &
            = old_blk_p (old_row_p(row)+1 : old_row_p(row+1))
       ncols = new_last_blk - new_first_blk + 1
       ! Sort
       IF (ncols .GT. 0) THEN
          IF (.NOT.copy_data) THEN
             new_blk_d(new_first_blk : new_first_blk + old_nblks - 1) = 1
          ENDIF
          permutation(1:ncols) = (/ (i, i=1,ncols) /)
          buff(1:ncols) = new_blk_p(new_first_blk:new_last_blk)
          IF (.NOT. copy_data) THEN
             buff2(1:ncols) = new_blk_d(new_first_blk:new_last_blk)
          ENDIF
          CALL sort (new_col_i(new_first_blk:new_last_blk), ncols, permutation)
          IF (copy_data) THEN
             DO i = 1, ncols
                new_blk_p(new_first_blk+i-1) = buff(permutation(i))
             ENDDO
          ELSE
             DO i = 1, ncols
                new_blk_p(new_first_blk+i-1) = buff(permutation(i))
                new_blk_d(new_first_blk+i-1) = buff2(permutation(i))
             ENDDO
          ENDIF
       ENDIF
    ENDDO
    !$OMP ENDDO
    DEALLOCATE (permutation)
    DEALLOCATE (buff)
    IF (.NOT. copy_data) DEALLOCATE (buff2)
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE index_insertion


! *****************************************************************************
!> \brief Merges a work matrix into the main matrix
!> \param[in,out] target_data_area   where to copy data
!> \param[in] wm                     work matrix to merge
!> \param[in] new_blk_p              block offsets and sizes for the final
!>                                   matrix
!> \param[in,out] error              error
! *****************************************************************************
  SUBROUTINE merge_my_wm (target_data_area, wm, new_blk_p, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: target_data_area
    TYPE(dbcsr_work_type), INTENT(IN)        :: wm
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: new_blk_p
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'merge_my_wm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    CALL merge_my_wm_linear (target_data_area,&
         wm%data_area, wm%lastblk, wm%blk_p, new_blk_p,&
         error=error)
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE merge_my_wm


! *****************************************************************************
!> \brief Merges an array-based work matrix into the main matrix
!> \param[in,out] target_data_area   where to copy data
!> \param[in] source_data_area       from where to copy data
!> \param[in] old_blk_p              block offsets and sizes in the work
!>                                   matrix
!> \param[in] new_blk_p              block offsets and sizes for the final
!>                                   matrix
!> \param[in,out] error              error
! *****************************************************************************
  SUBROUTINE merge_my_wm_linear (target_data_area,&
       source_data_area, source_nblocks, old_blk_p, new_blk_p, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: target_data_area
    TYPE(dbcsr_data_obj), INTENT(IN)         :: source_data_area
    INTEGER, INTENT(IN)                      :: source_nblocks
    INTEGER, DIMENSION(:), INTENT(IN)        :: old_blk_p
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: new_blk_p
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'merge_my_wm_linear', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, error_handle
    INTEGER, DIMENSION(1)                    :: nze_array

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    DO blk = 1, source_nblocks
       nze_array(1) = new_blk_p(2,blk)
       CALL dbcsr_data_set (target_data_area,&
            lb=ABS(new_blk_p(1,blk)), data_size=new_blk_p(2,blk),&
            src=source_data_area, source_lb=ABS(old_blk_p(blk)))
    ENDDO
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE merge_my_wm_linear

! *****************************************************************************
!> \brief Counts the number of blocks and data to be added
!>
!> \param[in,out] wm          work matrix
!> \param[in,out] additions_per_row     number of blocks and data to add for
!>                                      every row
!> \param[out] additions_total          total number of blocks and data to add
!> \param[out] blk_p          block sizes of data to add
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE count_blocks_to_add (wm, additions_per_row, additions_total,&
       blk_p, row_sizes, col_sizes, error)
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: wm
    INTEGER, DIMENSION(:, :)                 :: additions_per_row
    INTEGER, DIMENSION(2)                    :: additions_total
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: blk_p
    INTEGER, DIMENSION(:), INTENT(IN)        :: row_sizes, col_sizes
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'count_blocks_to_add', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: index_blocks = 1, &
                                                index_data = 2

    INTEGER                                  :: blk, col, error_handle, &
                                                nblks, nze, row

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    IF (dbcsr_wm_use_mutable (wm)) THEN
       SELECT CASE (wm%mutable%m%data_type)
       CASE (dbcsr_type_real_4)
          CALL tree_to_linear_s (wm, error=error)
       CASE (dbcsr_type_real_8)
          CALL tree_to_linear_d (wm, error=error)
       CASE (dbcsr_type_complex_4)
          CALL tree_to_linear_c (wm, error=error)
       CASE (dbcsr_type_complex_8)
          CALL tree_to_linear_z (wm, error=error)
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
               routineN, "Invalid data type",__LINE__,error)
       END SELECT
    ENDIF
    nblks = wm%lastblk
    ! Sort
    CALL dbcsr_sort_indices (nblks, wm%row_i, wm%col_i, wm%blk_p)
    !
    additions_total(:) = 0
    DO blk = 1, wm%lastblk
       row = wm%row_i(blk)
       col = wm%col_i(blk)
       nze = row_sizes(row) * col_sizes(col)
       blk_p(1,blk) = wm%blk_p(blk) ! offset will go here
       blk_p(2,blk) = nze
       additions_per_row(index_blocks, row)&
            = additions_per_row(index_blocks, row) + 1
       additions_per_row(index_data, row)&
            = additions_per_row(index_data, row) + nze
       additions_total(index_blocks) = additions_total(index_blocks) + 1
       additions_total(index_data) = additions_total(index_data) + nze
    ENDDO
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE count_blocks_to_add


! *****************************************************************************
!> \brief Adds (a thread's) count of blocks to add to total count
!>
!> \param[in,out] row_counts            total counts of blocks to add per row
!> \param[in]additions_per_row          number of blocks and data to add for
!>                                      every row (for current thread/work
!>                                      matrix)
! *****************************************************************************
  SUBROUTINE add_to_row_count (row_counts, additions_per_row)
    INTEGER, DIMENSION(:), POINTER           :: row_counts
    INTEGER, DIMENSION(:), INTENT(IN)        :: additions_per_row

    INTEGER                                  :: nrows, row

!   ---------------------------------------------------------------------------

    nrows = SIZE (additions_per_row)
    DO row = 1, nrows
       IF (additions_per_row(row) .GT. 0) THEN
          !$OMP ATOMIC
          row_counts(row) = row_counts(row) + additions_per_row(row)
       ENDIF
    ENDDO
  END SUBROUTINE add_to_row_count


! *****************************************************************************
!> \brief Computes offsets for each work matrix
!> \param[in,out] all_offsets    array of all offsets
!> \param[in] my_size            my thead's size
!> \param[in] first_offset       the offset for the first thread
! *****************************************************************************
  SUBROUTINE tally_offsets (all_offsets, my_size, first_offset)
    INTEGER, DIMENSION(:), POINTER           :: all_offsets
    INTEGER, INTENT(IN)                      :: my_size, first_offset

    INTEGER                                  :: me, o, old_size

!   ---------------------------------------------------------------------------

    me = 1
    !$ me = omp_get_thread_num() + 1
    all_offsets(me) = my_size
!$OMP BARRIER
    !$OMP MASTER
    old_size = all_offsets(1)
    all_offsets(1) = first_offset
    DO me = 2, SIZE (all_offsets)
       o = all_offsets(me)
       all_offsets(me) = all_offsets(me-1) + old_size
       old_size = o
    ENDDO
    !$OMP END MASTER
    !$OMP BARRIER
  END SUBROUTINE tally_offsets
  


#include "dbcsr_work_operations_d.F"
#include "dbcsr_work_operations_z.F"
#include "dbcsr_work_operations_s.F"
#include "dbcsr_work_operations_c.F"


END MODULE dbcsr_work_operations

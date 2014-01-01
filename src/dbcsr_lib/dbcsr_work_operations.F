!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
                                             array_nullify,&
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
  USE dbcsr_block_operations,          ONLY: block_copy_c,&
                                             block_copy_d,&
                                             block_copy_s,&
                                             block_copy_z,&
                                             dbcsr_data_copy,&
                                             dbcsr_data_set
  USE dbcsr_config,                    ONLY: default_resize_factor
  USE dbcsr_data_methods,              ONLY: &
       dbcsr_data_clear_pointer, dbcsr_data_ensure_size, &
       dbcsr_data_get_memory_type, dbcsr_data_get_size, &
       dbcsr_data_get_size_referenced, dbcsr_data_get_type, dbcsr_data_hold, &
       dbcsr_data_init, dbcsr_data_new, dbcsr_data_release, &
       dbcsr_data_set_size_referenced, dbcsr_get_data_p_c, &
       dbcsr_get_data_p_d, dbcsr_get_data_p_s, dbcsr_get_data_p_z
  USE dbcsr_data_operations,           ONLY: dbcsr_data_copyall,&
                                             dbcsr_sort_data
  USE dbcsr_error_handling
  USE dbcsr_index_operations,          ONLY: dbcsr_addto_index_array,&
                                             dbcsr_build_row_index,&
                                             dbcsr_clearfrom_index_array,&
                                             dbcsr_index_prune_deleted,&
                                             dbcsr_make_dbcsr_index,&
                                             dbcsr_make_index_exist,&
                                             dbcsr_repoint_index,&
                                             dbcsr_sort_indices
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_kinds,                     ONLY: default_string_length,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE dbcsr_mem_methods,               ONLY: dbcsr_memtype_equal
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_sizes, dbcsr_distribution, &
       dbcsr_distribution_has_threads, dbcsr_distribution_hold, &
       dbcsr_distribution_make_threads, dbcsr_distribution_ncols, &
       dbcsr_distribution_nrows, dbcsr_distribution_release, &
       dbcsr_get_data_memory_type, dbcsr_get_data_size_used, &
       dbcsr_get_data_type, dbcsr_get_index_memory_type, &
       dbcsr_get_matrix_type, dbcsr_get_replication_type, dbcsr_init, &
       dbcsr_is_initialized, dbcsr_matrix_counter, dbcsr_mutable_destroy, &
       dbcsr_mutable_init, dbcsr_mutable_instantiated, dbcsr_mutable_new, &
       dbcsr_mutable_release, dbcsr_name, dbcsr_row_block_sizes, &
       dbcsr_switch_data_area, dbcsr_use_mutable, dbcsr_valid_index, &
       dbcsr_wm_use_mutable
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size
  USE dbcsr_toollib,                   ONLY: dbcsr_unpack_i8_2i4,&
                                             sort,&
                                             uppercase
  USE dbcsr_types,                     ONLY: &
       dbcsr_data_obj, dbcsr_distribution_obj, dbcsr_iterator, &
       dbcsr_memtype_default, dbcsr_memtype_type, dbcsr_meta_size, &
       dbcsr_num_slots, dbcsr_obj, dbcsr_repl_col, dbcsr_repl_full, &
       dbcsr_repl_none, dbcsr_repl_row, dbcsr_slot_blk_p, dbcsr_slot_col_i, &
       dbcsr_slot_home_coli, dbcsr_slot_home_pcol, dbcsr_slot_home_prow, &
       dbcsr_slot_home_rowi, dbcsr_slot_home_vpcol, dbcsr_slot_home_vprow, &
       dbcsr_slot_nblkcols_local, dbcsr_slot_nblkcols_total, &
       dbcsr_slot_nblkrows_local, dbcsr_slot_nblkrows_total, &
       dbcsr_slot_nblks, dbcsr_slot_nfullcols_local, &
       dbcsr_slot_nfullcols_total, dbcsr_slot_nfullrows_local, &
       dbcsr_slot_nfullrows_total, dbcsr_slot_nze, dbcsr_slot_row_p, &
       dbcsr_slot_size, dbcsr_type, dbcsr_type_antihermitian, &
       dbcsr_type_antisymmetric, dbcsr_type_complex_4, dbcsr_type_complex_8, &
       dbcsr_type_hermitian, dbcsr_type_no_symmetry, dbcsr_type_real_4, &
       dbcsr_type_real_8, dbcsr_type_real_default, dbcsr_type_symmetric, &
       dbcsr_work_type
  USE dbcsr_util,                      ONLY: convert_sizes_to_offsets,&
                                             dbcsr_calc_block_sizes,&
                                             dbcsr_verify_matrix,&
                                             meta_from_dist

  !$ USE OMP_LIB

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_work_operations'

  LOGICAL, PARAMETER :: careful_mod = .FALSE.



  PUBLIC :: dbcsr_create, dbcsr_work_create, dbcsr_finalize, dbcsr_special_finalize
  PUBLIC :: dbcsr_work_destroy, dbcsr_add_wm_from_matrix,&
            add_work_coordinate

  INTERFACE dbcsr_create
     MODULE PROCEDURE dbcsr_create_new, dbcsr_create_template
  END INTERFACE


  TYPE i_array_p
     INTEGER, DIMENSION(:), POINTER :: p
  END TYPE i_array_p


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
!> \param[in] memory_type     (optional) allocate indices and data using
!>                            special memory
!> \param[in] index_memory_type     (optional) allocate indices using
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
       row_blk_size, col_blk_size, nblks, nze, data_type,&
       data_memory_type, index_memory_type, reuse,&
       mutable_work, replication_type, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    CHARACTER(len=*), INTENT(IN)             :: name
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    CHARACTER, INTENT(IN)                    :: matrix_type
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_size, col_blk_size
    INTEGER, INTENT(IN), OPTIONAL            :: nblks, nze, data_type
    TYPE(dbcsr_memtype_type), INTENT(IN), &
      OPTIONAL                               :: data_memory_type, &
                                                index_memory_type
    LOGICAL, INTENT(IN), OPTIONAL            :: reuse, mutable_work
    CHARACTER, INTENT(IN), OPTIONAL          :: replication_type
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_new', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: matrix_type_l
    INTEGER                                  :: error_handler, my_nblks, &
                                                my_nze
    INTEGER, DIMENSION(:), POINTER           :: col_blk_offset, row_blk_offset
    INTEGER, DIMENSION(dbcsr_meta_size)      :: new_meta
    LOGICAL                                  :: hijack

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert(dbcsr_is_initialized(matrix), &
           dbcsr_fatal_level, dbcsr_caller_error, routineN,&
           "Matrix may not have been initialized with dbcsr_init",__LINE__,error)

    ! Reuse matrix only if has actually been allocated.
    hijack = ASSOCIATED(matrix%m%index)
    IF (PRESENT (reuse)) hijack = reuse

    IF (.NOT.hijack) THEN
       CALL dbcsr_init (matrix%m)
       matrix%m%refcount = 1
    ELSEIF (dbcsr_buffers_2d_needed) THEN
       CALL dbcsr_buffers_flush (matrix%m%buffers, error=error)
       CALL dbcsr_buffers_release (matrix%m%buffers, error=error)
    ENDIF
    !$OMP CRITICAL
    matrix%m%serial_number = dbcsr_matrix_counter
    dbcsr_matrix_counter = dbcsr_matrix_counter + 1
    !$OMP END CRITICAL
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
       IF (matrix%m%has_local_rows) &
            CALL array_release (matrix%m%local_rows)
       IF (matrix%m%has_global_rows) &
            CALL array_release (matrix%m%global_rows)
       IF (matrix%m%has_local_cols) &
            CALL array_release (matrix%m%local_cols)
       IF (matrix%m%has_global_cols) &
            CALL array_release (matrix%m%global_cols)
       CALL dbcsr_distribution_release (matrix%m%dist)
       IF (ASSOCIATED (matrix%m%wms)) THEN
          CALL dbcsr_work_destroy_all(matrix%m)
       ENDIF
       CALL array_nullify (matrix%m%local_rows)
       CALL array_nullify (matrix%m%global_rows)
       CALL array_nullify (matrix%m%local_cols)
       CALL array_nullify (matrix%m%global_cols)
    ELSE
       ! Invalidate index
       NULLIFY(matrix%m%index)
       ! Invalidate data
       CALL dbcsr_data_init (matrix%m%data_area)
    ENDIF
    ! These are always invalidated.
    NULLIFY(matrix%m%row_p, matrix%m%col_i, matrix%m%blk_p, matrix%m%thr_c,&
         matrix%m%coo_l)
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
    !
    CALL dbcsr_assert (array_size (row_blk_size), "EQ",&
         dbcsr_distribution_nrows (dist),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Number of blocked rows does match blocked row distribution.",&
         __LINE__, error=error)
    CALL dbcsr_assert (array_size (col_blk_size), "EQ",&
         dbcsr_distribution_ncols (dist),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Number of blocked columns does match blocked column distribution.",&
         __LINE__, error=error)
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

    matrix%m%data_memory_type = dbcsr_memtype_default
    IF (PRESENT (data_memory_type)) &
       matrix%m%data_memory_type = data_memory_type

    matrix%m%index_memory_type = dbcsr_memtype_default
    IF (PRESENT (index_memory_type)) &
       matrix%m%index_memory_type = index_memory_type

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
    ! Setup a matrix from scratch
    IF (.NOT. hijack) THEN
       CALL dbcsr_data_new (matrix%m%data_area, matrix%m%data_type, my_nze,&
            memory_type=matrix%m%data_memory_type)
       CALL dbcsr_data_set_size_referenced (matrix%m%data_area, 0)
       !
       NULLIFY (matrix%m%index)
       CALL ensure_array_size (matrix%m%index, lb=1, ub=dbcsr_num_slots,&
            zero_pad = .TRUE., memory_type = matrix%m%index_memory_type,&
            error=error)
    ENDIF
    CALL dbcsr_assert (LBOUND (matrix%m%index, 1) .LE. 1&
         .AND. UBOUND (matrix%m%index, 1) .GE. dbcsr_num_slots,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Index is not large enough",__LINE__,error)
    matrix%m%index(1:dbcsr_num_slots) = 0
    matrix%m%index(1:dbcsr_meta_size) = new_meta(1:dbcsr_meta_size)
    matrix%m%index(dbcsr_slot_size) = dbcsr_num_slots
    !
    matrix%m%symmetry = .FALSE.
    matrix%m%negate_real = .FALSE.
    matrix%m%negate_imaginary = .FALSE.
    !matrix%m%transpose = .FALSE.
    matrix_type_l = matrix_type
    CALL uppercase(matrix_type_l)
    SELECT CASE (matrix_type_l)
    CASE (dbcsr_type_no_symmetry)
    CASE (dbcsr_type_symmetric)
       matrix%m%symmetry = .TRUE.
    CASE (dbcsr_type_antisymmetric)
       matrix%m%symmetry = .TRUE.
       matrix%m%negate_real = .TRUE.
       matrix%m%negate_imaginary = .TRUE.
    CASE (dbcsr_type_hermitian)
       matrix%m%symmetry = .TRUE.
       matrix%m%negate_imaginary = .TRUE.
    CASE (dbcsr_type_antihermitian)
       matrix%m%symmetry = .TRUE.
       matrix%m%negate_real = .TRUE.
    CASE DEFAULT
       CALL dbcsr_assert(.FALSE., dbcsr_failure_level,&
            dbcsr_wrong_args_error, routineP, "Invalid matrix type.",__LINE__,error)
    END SELECT
    matrix%m%bcsc = .FALSE.
    matrix%m%local_indexing = .FALSE.
    matrix%m%list_indexing = .FALSE.
    CALL array_nullify (matrix%m%local_rows)
    CALL array_nullify (matrix%m%global_rows)
    CALL array_nullify (matrix%m%local_cols)
    CALL array_nullify (matrix%m%global_cols)
    matrix%m%has_local_rows  = .FALSE.
    matrix%m%has_global_rows = .FALSE.
    matrix%m%has_local_cols  = .FALSE.
    matrix%m%has_global_cols = .FALSE.
    CALL dbcsr_buffers_init (matrix%m%buffers)
    IF (dbcsr_buffers_2d_needed) THEN
       CALL dbcsr_buffers_new (matrix%m%buffers, matrix%m%data_area,&
            error=error)
    ENDIF
!$  CALL OMP_INIT_LOCK (matrix%m%modification_lock)
    CALL dbcsr_make_index_exist (matrix%m, error=error)
    matrix%m%valid = .TRUE.
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_create_new


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE dbcsr_create_template(matrix, template, name, dist, matrix_type,&
       row_blk_size, col_blk_size, nblks, nze, data_type,&
       data_memory_type, index_memory_type, reuse,&
       mutable_work, replication_type, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_obj), INTENT(IN)              :: template
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: name
    TYPE(dbcsr_distribution_obj), &
      INTENT(IN), OPTIONAL                   :: dist
    CHARACTER, INTENT(IN), OPTIONAL          :: matrix_type
    TYPE(array_i1d_obj), INTENT(IN), &
      OPTIONAL                               :: row_blk_size, col_blk_size
    INTEGER, INTENT(IN), OPTIONAL            :: nblks, nze, data_type
    TYPE(dbcsr_memtype_type), INTENT(IN), &
      OPTIONAL                               :: data_memory_type, &
                                                index_memory_type
    LOGICAL, INTENT(IN), OPTIONAL            :: reuse, mutable_work
    CHARACTER, INTENT(IN), OPTIONAL          :: replication_type
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_create_template', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: new_matrix_type, &
                                                new_replication_type
    CHARACTER(len=default_string_length)     :: new_name
    INTEGER                                  :: new_data_type
    LOGICAL                                  :: new_mutable_work
    TYPE(array_i1d_obj)                      :: new_col_blk_size, &
                                                new_row_blk_size
    TYPE(dbcsr_distribution_obj)             :: new_dist
    TYPE(dbcsr_memtype_type)                 :: new_data_memory_type, &
                                                new_index_memory_type

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
    IF (PRESENT (data_memory_type)) THEN
       new_data_memory_type = data_memory_type
    ELSE
       new_data_memory_type = dbcsr_get_data_memory_type (template)
    ENDIF
    IF (PRESENT (index_memory_type)) THEN
       new_index_memory_type = index_memory_type
    ELSE
       new_index_memory_type = dbcsr_get_index_memory_type (template)
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
         data_memory_type = new_data_memory_type,&
         index_memory_type = new_index_memory_type,&
         mutable_work = new_mutable_work,&
         replication_type = new_replication_type,&
         error=error)
    ! Copy stuff from the meta-array.  These are not normally needed,
    ! but have to be here for creating matrices from "image" matrices.
    matrix%m%index(dbcsr_slot_home_prow) = template%m%index(dbcsr_slot_home_prow)
    matrix%m%index(dbcsr_slot_home_rowi) = template%m%index(dbcsr_slot_home_rowi)
    matrix%m%index(dbcsr_slot_home_pcol) = template%m%index(dbcsr_slot_home_pcol)
    matrix%m%index(dbcsr_slot_home_coli) = template%m%index(dbcsr_slot_home_coli)
    matrix%m%index(dbcsr_slot_home_vprow) = template%m%index(dbcsr_slot_home_vprow)
    matrix%m%index(dbcsr_slot_home_vpcol) = template%m%index(dbcsr_slot_home_vpcol)
  END SUBROUTINE dbcsr_create_template


! *****************************************************************************
!> \brief Initializes one work matrix
!> \param[out] wm             initialized work matrix
!> \param[in] nblks_guess     (optional) estimated number of blocks
!> \param[in] sizedata_guess  (optional) estimated size of data
! *****************************************************************************
  SUBROUTINE dbcsr_init_wm (wm, data_type, nblks_guess, sizedata_guess, memory_type, error)
    TYPE(dbcsr_work_type), INTENT(OUT)       :: wm
    INTEGER, INTENT(IN)                      :: data_type
    INTEGER, INTENT(IN), OPTIONAL            :: nblks_guess, sizedata_guess
    TYPE(dbcsr_memtype_type), INTENT(IN), &
      OPTIONAL                               :: memory_type
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_wm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, nblks, stat

!   ---------------------------------------------------------------------------

    IF(careful_mod) CALL dbcsr_error_set(routineN, error_handler, error)
    wm%lastblk = 0
    wm%datasize = 0
    ! Index
    IF(PRESENT(nblks_guess)) THEN
       nblks = nblks_guess
       CALL dbcsr_assert (nblks_guess, "GE", 0, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN, &
            "Can not have negative block count.", __LINE__, error=error)
       ALLOCATE(wm%row_i(nblks), stat=stat)
       CALL dbcsr_assert (stat == 0, dbcsr_fatal_level,&
            dbcsr_internal_error, routineN,&
            "wm%row_i",__LINE__,error)
       ALLOCATE(wm%col_i(nblks), stat=stat)
       CALL dbcsr_assert (stat == 0, dbcsr_fatal_level,&
            dbcsr_internal_error, routineN,&
            "wm%col_i",__LINE__,error)
       ALLOCATE(wm%blk_p(nblks), stat=stat)
       CALL dbcsr_assert (stat == 0, dbcsr_fatal_level,&
            dbcsr_internal_error, routineN,&
            "wm%blk_p",__LINE__,error)
    ELSE
       NULLIFY (wm%row_i, wm%col_i, wm%blk_p)
       !nblks = CEILING (REAL (matrix%nblkrows_local * matrix%nblkcols_local)&
       !     / REAL (dbcsr_mp_numnodes (dbcsr_distribution_mp (matrix%dist))))
    ENDIF
    ! Data
    CALL dbcsr_data_init (wm%data_area)
    IF(PRESENT(sizedata_guess)) THEN
       CALL dbcsr_assert (sizedata_guess, "GE", 0, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Can not have negative data size.", __LINE__, error=error)
       CALL dbcsr_data_new (wm%data_area, data_type,&
            data_size=sizedata_guess, memory_type=memory_type)
    ELSE
       CALL dbcsr_data_new (wm%data_area, data_type, memory_type=memory_type)
    ENDIF
    CALL dbcsr_mutable_init (wm%mutable)
    IF(careful_mod) CALL dbcsr_error_stop(error_handler, error)
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
       work_mutable, memory_type, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN), OPTIONAL            :: nblks_guess, sizedata_guess, n
    LOGICAL, INTENT(in), OPTIONAL            :: work_mutable
    TYPE(dbcsr_memtype_type), INTENT(IN), &
      OPTIONAL                               :: memory_type
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_work_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, iw, nw, ow
    LOGICAL                                  :: wms_new, wms_realloc
    TYPE(dbcsr_work_type), DIMENSION(:), &
      POINTER                                :: wms

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_is_initialized(matrix),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Can not create work matrix for nonexisting matrix object.",__LINE__,error)
    IF (PRESENT (n)) THEN
       nw = n
    ELSE
       nw = 1
       !$ IF (omp_in_parallel()) THEN
       !$    nw = omp_get_num_threads()
       !$ ELSE
       !$    nw = omp_get_max_threads()
       !$ ENDIF
    ENDIF
    !$OMP MASTER
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
               nblks_guess=nblks_guess, sizedata_guess=sizedata_guess,&
               memory_type=memory_type, error=error)
          IF (matrix%m%work_mutable) &
               CALL dbcsr_mutable_new (matrix%m%wms(iw)%mutable,&
               dbcsr_get_data_type (matrix))
       END DO
    ENDIF
    IF (wms_new) THEN
       ALLOCATE (matrix%m%wms(nw))
       DO iw = 1, nw
          CALL dbcsr_init_wm (matrix%m%wms(iw), matrix%m%data_type,&
               nblks_guess=nblks_guess, sizedata_guess=sizedata_guess,&
               memory_type=memory_type, error=error)
          IF (matrix%m%work_mutable) &
               CALL dbcsr_mutable_new (matrix%m%wms(iw)%mutable,&
               dbcsr_get_data_type (matrix))
       END DO
    ENDIF
    matrix%m%valid = .FALSE.
    !$OMP END MASTER
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_work_create



! *****************************************************************************
!> \brief Creates the final dbcsr_type matrix from the working matrix.
!>
!> Work matrices (array or tree-based) are merged into the base DBCSR matrix.
!>
!> If a matrix is marked as having a valid index, then nothing is done.
!>
!> Deleted blocks are pruned from the index.
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

    INTEGER                                  :: error_handler, i, nblks, &
                                                nwms, start_offset
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

    ! If the matrix is not marked as dirty then skip the work.
    IF (dbcsr_valid_index(matrix)) THEN
       !"No need to finalize a valid matrix, skipping."
       !
       ! A matrix with a valid index should not have associated work
       ! arrays.  This may happen when this routine is called on a
       ! matrix that was not changed.
       !$OMP BARRIER
       !$OMP MASTER
       IF (ASSOCIATED (matrix%m%wms)) &
            CALL dbcsr_work_destroy_all(matrix%m)
       matrix%m%valid = .TRUE.
       !$OMP END MASTER
       !$OMP BARRIER
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
    start_offset = matrix%m%nze
    i = dbcsr_get_data_size_used (matrix, error)
    !$OMP MASTER
    matrix%m%nze = 0
    !$OMP END MASTER
    !$OMP BARRIER
    !$OMP ATOMIC
    matrix%m%nze = matrix%m%nze + i
    !$OMP BARRIER
    IF (dbg) THEN
       WRITE(*,*)routineN//" sizes", matrix%m%nze, i,&
            dbcsr_data_get_size_referenced (matrix%m%data_area), &
            dbcsr_data_get_size (matrix%m%data_area)
    ENDIF
    IF (.FALSE. .AND. dbcsr_data_get_size_referenced (matrix%m%data_area) .NE. &
         matrix%m%nze) THEN
       CALL dbcsr_assert (matrix%m%nze, "EQ", &
            dbcsr_data_get_size_referenced (matrix%m%data_area),&
            dbcsr_warning_level, dbcsr_caller_error, routineN,&
            "Should reshuffle.", __LINE__, error=error)
       IF (ASSOCIATED (matrix%m%wms)) THEN
          sort_data = .NOT. dbcsr_wm_use_mutable (matrix%m%wms(1))
       ENDIF
    ENDIF
    IF (sort_data .AND. matrix%m%nze .GT. 0) THEN
       CALL dbcsr_add_wm_from_matrix (matrix, error=error)
       matrix%m%nze = 0
       !$OMP MASTER
       fake_row_p  = .TRUE.
       !$OMP END MASTER
    ENDIF
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
    !$OMP END MASTER
    !
    ! If the matrix, work matricies, and environment fit several
    ! criteria, then a quick O(1) finalization is performed.
    IF (can_quickly_finalize (matrix) .AND. .NOT. sort_data) THEN
       !$OMP BARRIER
       CALL quick_finalize (matrix, error=error)
    ELSE
       !
       !$OMP MASTER
       !
       ! Create work matrices if not yet existing
       IF (.NOT. ASSOCIATED (matrix%m%wms)) THEN
          nwms = 1
          !$ nwms = omp_get_num_threads()
          CALL dbcsr_work_create (matrix, n=nwms, error=error)
       ENDIF
       !$OMP END MASTER
       !$OMP BARRIER
       !
       ! Ensure index arrays at least exist.
       !$OMP DO SCHEDULE (STATIC, 1)
       DO i = 1, SIZE(matrix%m%wms)
          IF (.NOT. ASSOCIATED (matrix%m%wms(i)%row_i)) THEN
             CALL ensure_array_size(matrix%m%wms(i)%row_i, ub=0, error=error)
          ENDIF
          IF (.NOT. ASSOCIATED (matrix%m%wms(i)%col_i)) THEN
             CALL ensure_array_size(matrix%m%wms(i)%col_i, ub=0, error=error)
          ENDIF
          IF (.NOT. ASSOCIATED (matrix%m%wms(i)%blk_p)) THEN
             CALL ensure_array_size(matrix%m%wms(i)%blk_p, ub=0, error=error)
          ENDIF
       ENDDO
       !$OMP ENDDO
       !
       ! Check for deleted blocks
       !$OMP MASTER
       nblks = matrix%m%row_p(matrix%m%nblkrows_total+1)
       IF (ANY(matrix%m%blk_p(1:nblks) .EQ. 0)) THEN
          CALL dbcsr_index_prune_deleted (matrix, error)
       ENDIF
       old_row_p => matrix%m%row_p
       old_col_i => matrix%m%col_i
       old_blk_p => matrix%m%blk_p
       !$OMP END MASTER
       !
       !$OMP BARRIER
       ! Check to see if we will need to create a parallel environment
       ! (needed when there are multiple work matrices but we are not
       ! in an OpenMP parallel section.)
       !
       ! A parallel section is created is used when the matrix has
       ! more work matrices. It's a shortcut when the finalize is
       ! called from a non-parallel environment whereas the matrix was
       ! built/modified in a parallel environment
       nwms = SIZE (matrix%m%wms)
       spawn = .FALSE.
       !$ IF (.NOT. OMP_IN_PARALLEL ()) THEN
       !$    IF (nwms .GT. 1) spawn = .TRUE.
       !$ ENDIF
       IF (spawn) THEN
          !$OMP PARALLEL IF (spawn) &
          !$OMP          DEFAULT (none) &
          !$OMP          SHARED (matrix, old_row_p, old_col_i, old_blk_p,&
          !$OMP                  start_offset, sort_data, error)
          CALL dbcsr_merge_all (matrix%m,&
               old_row_p, old_col_i, old_blk_p,&
               data_starting_offset = start_offset,&
               sort_data=sort_data, &
               error=error)
          !$OMP END PARALLEL
       ELSE
          CALL dbcsr_merge_all (matrix%m,&
               old_row_p, old_col_i, old_blk_p,&
               data_starting_offset = start_offset,&
               sort_data=sort_data, &
               error=error)
       ENDIF
    ENDIF
!$OMP BARRIER
!$OMP MASTER
    ! Clean up.
    IF (ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_work_destroy_all(matrix%m)
    ENDIF
    matrix%m%valid = .TRUE.
!$OMP END MASTER
!$OMP BARRIER
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
! *****************************************************************************
  SUBROUTINE dbcsr_special_finalize(matrix, reshuffle, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: reshuffle
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_special_finalize', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handler
    LOGICAL                                  :: sort_data

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (dbcsr_is_initialized (matrix),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Can not finalize uninitialized matrix.",__LINE__,error)
    CALL dbcsr_assert (matrix%m%nblks .EQ. 0,&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Optimized finalize requires empty matrix.",__LINE__,error)
    CALL dbcsr_assert (.NOT. dbcsr_valid_index(matrix),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Optimized finalize requires invalid matrix.",__LINE__,error)
    CALL dbcsr_assert (ASSOCIATED (matrix%m%wms),&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Optimized finalize requires work matrices exist.",__LINE__,error)
    CALL dbcsr_assert (SIZE(matrix%m%wms).EQ. 1,&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Optimized finalize requires single work matrix",__LINE__,error)

    IF (dbcsr_buffers_2d_needed) THEN
       CALL dbcsr_buffers_flush (matrix%m%buffers, error=error)
    ENDIF
    !
    ! If possible, data copying is avoided.
    IF (PRESENT (reshuffle)) THEN
       sort_data = reshuffle
    ELSE
       sort_data = .FALSE.
    ENDIF

    !$OMP BARRIER

    ! If the matrix, work matricies, and environment fit several
    ! criteria, then a quick O(1) finalization is performed.
    IF (can_quickly_finalize (matrix) .AND. .NOT. sort_data) THEN
       CALL quick_finalize (matrix, error=error)
    ELSE
       CALL dbcsr_assert (sort_data,&
         dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "merge_single_wm only supports data sorting",__LINE__,error)
       !
       ! Ensure needed index arrays at least exist.
       !$OMP MASTER
       !
       IF (.NOT. ASSOCIATED (matrix%m%wms(1)%row_i)) THEN
          CALL ensure_array_size(matrix%m%wms(1)%row_i, ub=0, error=error)
       ENDIF
       IF (.NOT. ASSOCIATED (matrix%m%wms(1)%col_i)) THEN
          CALL ensure_array_size(matrix%m%wms(1)%col_i, ub=0, error=error)
       ENDIF
       IF (.NOT. ASSOCIATED (matrix%m%wms(1)%blk_p)) THEN
          CALL ensure_array_size(matrix%m%wms(1)%blk_p, ub=0, error=error)
       ENDIF
       !
       !$OMP END MASTER
       !$OMP BARRIER
       !
       !$OMP PARALLEL DEFAULT (none), SHARED(matrix, error)
         CALL dbcsr_merge_single_wm (matrix%m,&
            error=error)
       !$OMP END PARALLEL

    ENDIF

    !$OMP MASTER
    ! Clean up.
    IF (ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_work_destroy_all(matrix%m)
    ENDIF
    matrix%m%valid = .TRUE.
    !$OMP END MASTER
    !$OMP BARRIER
    IF (dbg) THEN
       !$OMP SINGLE
       CALL dbcsr_verify_matrix (matrix, error=error)
       !$OMP END SINGLE
    ENDIF
    IF (dbcsr_buffers_2d_needed) THEN
       CALL dbcsr_buffers_new (matrix%m%buffers, matrix%m%data_area,&
            error=error)
    ENDIF
    !$OMP BARRIER
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_special_finalize

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
          quick = quick .AND. &
               dbcsr_memtype_equal( &
                  dbcsr_data_get_memory_type(matrix%m%wms(1)%data_area), &
                 dbcsr_data_get_memory_type(matrix%m%data_area))
          quick = quick .AND. &
               ASSOCIATED (matrix%m%wms(1)%row_i)
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
    CALL dbcsr_clearfrom_index_array (matrix%m, dbcsr_slot_row_p)
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
                                                nthreads, nwms, old_nwms, &
                                                size_used
    LOGICAL                                  :: preexists

!   ---------------------------------------------------------------------------

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
    !$OMP MASTER
    IF (.NOT. ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_work_create (matrix, &
              INT(matrix%m%nblks * default_resize_factor / nwms), &
              INT(matrix%m%nze * default_resize_factor / nwms), &
              n=nwms, work_mutable=.FALSE., error=error)
    ENDIF
    !$OMP END MASTER
    !$OMP BARRIER
    size_used = matrix%m%nze
    CALL dbcsr_fill_wm_from_matrix (matrix%m%wms, matrix, size_used,&
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
  SUBROUTINE dbcsr_fill_wm_from_matrix(wm, matrix, size_used, limits, error)
    TYPE(dbcsr_work_type), DIMENSION(:), &
      INTENT(INOUT)                          :: wm
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: size_used
    INTEGER, DIMENSION(4), INTENT(IN), &
      OPTIONAL                               :: limits
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_fill_wm_from_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, blk_p, col, error_handler, ithread, nthreads, nwms, nze, &
      row, wblk_p, which_wm, wm_first, wm_last
    LOGICAL                                  :: careful, limit, mt, tr
    LOGICAL, SAVE                            :: mutable
    TYPE(dbcsr_data_obj)                     :: data_block
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    nwms = SIZE (matrix%m%wms)
    mt = .FALSE.
    !$ IF (nwms .GT. 1) mt = omp_get_num_threads() .GT. 1
    ithread = 0 ; nthreads = 1
    !$ ithread = omp_get_thread_num ()
    !$ nthreads = omp_get_num_threads()
    limit = PRESENT (limits)
    careful = size_used + size_used/8 &
         .LT. dbcsr_data_get_size_referenced (matrix%m%data_area)
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_new (data_block, dbcsr_data_get_type (matrix%m%data_area))
    IF (mt) THEN
       wm_first = ithread+1
       wm_last = ithread+1
    ELSE
       wm_first = 1
       wm_last = nwms
    ENDIF
    ! Prepares the work matrices to accept the main data.
    !$OMP MASTER
    mutable = .FALSE.
    DO which_wm = 1, nwms
       mutable = mutable .OR. dbcsr_wm_use_mutable (wm(which_wm))
    ENDDO
    !$OMP END MASTER
    !$OMP BARRIER
    DO which_wm = wm_first, wm_last
       CALL dbcsr_assert ("NOT", dbcsr_wm_use_mutable (wm(which_wm)),&
            dbcsr_fatal_level, dbcsr_unimplemented_error_nr, routineN,&
            "Adding main matrix into mutable not supported.", __LINE__,&
            error=error)
       IF (mutable) THEN
          IF (.NOT. dbcsr_mutable_instantiated (wm(which_wm)%mutable)) THEN
             CALL dbcsr_mutable_new (wm(which_wm)%mutable, matrix%m%data_type)
          ENDIF
       ELSE
          ! We don't know how much data we'll get so we have to be generous.
          CALL dbcsr_data_ensure_size (wm(which_wm)%data_area,&
               size_used/nwms, error=error)
          CALL dbcsr_data_set_size_referenced (wm(which_wm)%data_area, 0)
       ENDIF
    ENDDO
    ! Now copy the data
    CALL dbcsr_iterator_start (iter, matrix, shared=mt,&
         contiguous_pointers=.FALSE., read_only=.TRUE.)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, data_block,&
            transposed=tr, block_number=blk)
       IF (limit) THEN
          IF (.NOT. within_limits (row, col, limits)) CYCLE
       ENDIF
       blk_p = matrix%m%blk_p(blk)
       which_wm = ithread+1
       wblk_p = SIGN (wm(which_wm)%datasize+1, blk_p)
       nze = dbcsr_data_get_size (data_block)
       IF (mt .OR. limit .OR. careful .OR. mutable) THEN
          ! The data gets copied block by block so the block pointers
          ! are ordered accordingly.
          IF (.not.mutable) THEN
             CALL add_work_coordinate (wm(which_wm), row, col, wblk_p, error=error)
             CALL dbcsr_data_ensure_size (wm(which_wm)%data_area,&
                  ABS(wblk_p)+nze-1,factor=default_resize_factor, error=error)
             CALL dbcsr_data_set_size_referenced (wm(which_wm)%data_area,&
                  ABS(wblk_p) + nze - 1)
             CALL dbcsr_data_set (wm(which_wm)%data_area,&
                  lb=ABS(wblk_p),&
                  data_size=nze,&
                  src=data_block, source_lb=1)
          ENDIF
       ELSE
          ! The data gets copied all at once so the block pointers
          ! should remain the same as they were.
          CALL add_work_coordinate (wm(which_wm), row, col, blk_p, error=error)
       ENDIF
       IF (.NOT. mutable) &
            wm(which_wm)%datasize = wm(which_wm)%datasize+nze
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    CALL dbcsr_data_clear_pointer (data_block)
    CALL dbcsr_data_release (data_block)
    ! Copy all blocks at once
    IF (.NOT. mt .AND. .NOT. limit .AND. .NOT. careful .AND. .NOT. mutable) THEN
       DO which_wm = 1, nwms
          CALL dbcsr_data_ensure_size (wm(which_wm)%data_area,&
               dbcsr_data_get_size_referenced (matrix%m%data_area), error=error)
          CALL dbcsr_data_copyall (wm(which_wm)%data_area, matrix%m%data_area,&
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

    INTEGER                                  :: error_handler
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handler, error)

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

    IF (careful_mod) CALL dbcsr_error_stop(error_handler, error)

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

    INTEGER                                  :: error_handler, i
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (ASSOCIATED (m%wms), dbcsr_warning_level, dbcsr_caller_error,&
         routineN, "Want to destroy nonexisting work matrices.",__LINE__,error)
    IF (ASSOCIATED (m%wms)) THEN
       DO i = 1, SIZE (m%wms)
          CALL dbcsr_work_destroy (m%wms(i))
       ENDDO
       DEALLOCATE (m%wms)
       NULLIFY (m%wms)
    ENDIF

    CALL dbcsr_error_stop(error_handler, error)

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

!   ---------------------------------------------------------------------------

    IF (careful_mod) CALL dbcsr_error_set(routineN, error_handler, error)
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
    IF (PRESENT (index)) index = matrix%lastblk
    IF (careful_mod) CALL dbcsr_error_stop(error_handler, error)
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
    INTEGER, DIMENSION(*), INTENT(IN)        :: old_row_p, old_col_i, &
                                                old_blk_p
    INTEGER, INTENT(IN)                      :: data_starting_offset
    LOGICAL, INTENT(IN)                      :: sort_data
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_merge_all', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: index_blocks = 1, &
                                                index_data = 2
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handler, my_row_count, &
                                                nblks, new_data_size, nrows, &
                                                nwms, row, t, which_wm, &
                                                wm_datasize
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: all_data_offsets, &
                                                all_data_sizes, &
                                                new_blk_p_sorted, &
                                                new_blk_sizes, new_row_p
    INTEGER, ALLOCATABLE, DIMENSION(:), &
      SAVE, TARGET                           :: blk_d, new_blk_p, new_col_i
    INTEGER, DIMENSION(:), POINTER           :: my_row_p
    INTEGER, DIMENSION(:), POINTER, SAVE     :: cbs, rbs
    INTEGER, SAVE                            :: max_row_count, &
                                                nblks_to_add = 0, &
                                                new_nblks = 0, new_nze = 0
    TYPE(dbcsr_data_obj), ALLOCATABLE, &
      DIMENSION(:), SAVE                     :: all_data_areas
    TYPE(dbcsr_work_type), POINTER           :: wm
    TYPE(i_array_p), DIMENSION(:), POINTER, &
      SAVE                                   :: all_blk_p, all_col_i, &
                                                all_row_p

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    !> \par Outline

    !> Each thread has a work matrix.  These must be merged and made
    !> into a new index.  If sort_data is False, then the data areas
    !> are simply appended.  This is probably quicker but the data
    !> blocks are not in order and the data size may expand beyond
    !> what is needed.  If sort_data is True, then data blocks are
    !> sorted in order.

    IF (dbg) WRITE(*,*)routineN//" starting, sort?", sort_data
    !$OMP BARRIER
    nrows = matrix%nblkrows_total
    nwms = SIZE(matrix%wms)
    !$ call dbcsr_assert (nwms, "EQ", OMP_GET_NUM_THREADS (),&
    !$      dbcsr_fatal_level, dbcsr_caller_error, routineN,&
    !$      "Number of threads does not match number of work matrices.",&
    !$      __LINE__, error=error)
    which_wm = 1
    !$ which_wm = OMP_GET_THREAD_NUM () + 1
    wm => matrix%wms(which_wm)
    !
    ! Currently B-tree-based work matrices are converted to array form.
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
    !
    ! Initializations and some counts from the threads are summed.
    !
    !$OMP MASTER
    nblks_to_add = 0
    new_nblks = old_row_p(nrows+1)
    new_nze = matrix%nze
    ALLOCATE (all_row_p (nwms))
    ALLOCATE (all_col_i (nwms))
    ALLOCATE (all_blk_p (nwms))
    ALLOCATE (all_data_sizes (0:nwms))
    ALLOCATE (all_data_offsets (nwms))
    IF (sort_data) ALLOCATE (all_data_areas(0:nwms))
    IF (sort_data) THEN
       CALL dbcsr_data_init (all_data_areas(0))
       all_data_areas(0) = matrix%data_area
       !$OMP CRITICAL (crit_data)
       CALL dbcsr_data_hold (all_data_areas(0))
       !$OMP END CRITICAL (crit_data)
       !
       ! The following is valid because the data actually referenced
       ! by blocks is explicitly calculated in dbcsr_finalize.
       all_data_sizes(0) = matrix%nze
    ELSE
       all_data_sizes(0) = dbcsr_data_get_size_referenced(matrix%data_area)
    ENDIF
    rbs => array_data(matrix%row_blk_size)
    cbs => array_data(matrix%col_blk_size)
    !$OMP END MASTER
    !
    !$OMP BARRIER
    wm_datasize = wm%datasize
    nblks = wm%lastblk
    !$OMP ATOMIC
    nblks_to_add = nblks_to_add + nblks
    !$OMP ATOMIC
    new_nblks = new_nblks + nblks
    !$OMP ATOMIC
    new_nze = new_nze + wm_datasize
    !$OMP BARRIER
    !
    !$OMP MASTER
    ALLOCATE (new_row_p (nrows+1))
    ALLOCATE (new_col_i (new_nblks))
    ALLOCATE (new_blk_p (new_nblks))
    IF (sort_data) THEN
       ALLOCATE (blk_d(new_nblks))
    ELSE
       ALLOCATE (blk_d(1))
    ENDIF
    !$OMP END MASTER
    !
    ! Each thread creates a row_p index for its new blocks.  This gives the
    ! number of new blocks per thread per row.
    CALL dbcsr_sort_indices (nblks, wm%row_i, wm%col_i, wm%blk_p)
    ALLOCATE (my_row_p (nrows+1))
    CALL dbcsr_make_dbcsr_index (my_row_p, wm%row_i, nrows, nblks)
    !
    ! The threads must share their row_p, col_i, blk_p, and data areas
    ! among themselves.
    all_row_p(which_wm)%p => my_row_p
    IF (sort_data) THEN
       all_data_sizes(which_wm) = wm%datasize
    ELSE
       all_data_sizes(which_wm) = wm%datasize
       CALL dbcsr_assert (wm%datasize, "LE", dbcsr_data_get_size_referenced (wm%data_area),&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Data size mismatch.", __LINE__, error=error)
    ENDIF
    all_col_i(which_wm)%p => wm%col_i
    all_blk_p(which_wm)%p => wm%blk_p
    !$OMP BARRIER
    IF (sort_data) THEN
       !$OMP MASTER
       all_data_offsets(:) = 0
       !$OMP END MASTER
       CALL dbcsr_data_init (all_data_areas(which_wm))
       all_data_areas(which_wm) = wm%data_area
       !$OMP CRITICAL (crit_data)
       CALL dbcsr_data_hold (all_data_areas(which_wm))
       !$OMP END CRITICAL (crit_data)
    ELSE
       !$OMP MASTER
       all_data_offsets(1) = all_data_sizes(0)
       DO t = 2, nwms
          all_data_offsets(t) = all_data_offsets(t-1) + all_data_sizes(t-1)
       ENDDO
       !$OMP END MASTER
    ENDIF
    !
    ! The new row counts are created, then the new row_p index is created.
    !
    !$OMP DO
    DO row = 1, nrows
       my_row_count = old_row_p(row+1)-old_row_p(row)
       DO t = 1, nwms
          my_row_count = my_row_count &
               &       + all_row_p(t)%p(row+1)-all_row_p(t)%p(row)
       ENDDO
       new_row_p(row) = my_row_count
    ENDDO
    !$OMP END DO
    !$OMP MASTER
    max_row_count = MAXVAL (new_row_p(1:nrows))
    CALL dbcsr_build_row_index(new_row_p, nrows)
    !$OMP END MASTER
    !$OMP BARRIER
    !
    ! The new index is built
    CALL merge_index (new_row_p, new_col_i, new_blk_p,&
         blk_d,&
         old_row_p, old_col_i, old_blk_p,&
         all_row_p, all_col_i, all_blk_p,&
         all_data_offsets, nwms, nrows, max_row_count, sort_data, &
         error)
    !
    ! The data is then merged in one of two ways.
    !
    !$OMP MASTER
    new_data_size = SUM (all_data_sizes)
    IF (sort_data) THEN
       ! The matrix gets a fresh data area.  Blocks from the work
       ! matrices will be copied into it in order.
       !
       !$OMP CRITICAL (crit_data)
       CALL dbcsr_data_release (matrix%data_area)
       CALL dbcsr_data_init (matrix%data_area)
       !$OMP END CRITICAL (crit_data)
       CALL dbcsr_data_new (matrix%data_area,&
            data_size = new_data_size,&
            data_type=dbcsr_data_get_type (all_data_areas(0)),&
            memory_type = dbcsr_data_get_memory_type (all_data_areas(0)))
       ALLOCATE (new_blk_p_sorted (new_nblks))
       ALLOCATE (new_blk_sizes (new_nblks))
    ELSE
       ! Data stored in the work matrices will be just appended to the
       ! current data area.
       CALL dbcsr_data_ensure_size (matrix%data_area, new_data_size, error=error)
    ENDIF
    !$OMP END MASTER
    !$OMP BARRIER
    IF (sort_data) THEN
       CALL dbcsr_calc_block_sizes(new_blk_sizes,&
            new_row_p, new_col_i, rbs, cbs)
       CALL dbcsr_sort_data (new_blk_p_sorted, new_blk_p,&
            new_blk_sizes, dsts=matrix%data_area,&
            src=all_data_areas(0),&
            srcs=all_data_areas, old_blk_d=blk_d, error=error)
    ELSE
       IF (all_data_sizes(which_wm) .GT. 0) THEN
          CALL dbcsr_data_copy(dst=matrix%data_area,&
               dst_lb=(/ all_data_offsets(which_wm)+1 /),&
               dst_sizes=(/ all_data_sizes(which_wm) /),&
               src=wm%data_area,&
               src_lb=(/ 1 /),&
               src_sizes=(/ all_data_sizes(which_wm) /),&
               error=error)
       ENDIF
    ENDIF
    !
    ! Creates a new index array.
    !
    !$OMP MASTER
    CALL dbcsr_clearfrom_index_array (matrix, dbcsr_slot_row_p)
    CALL dbcsr_clearfrom_index_array (matrix, dbcsr_slot_col_i)
    CALL dbcsr_clearfrom_index_array (matrix, dbcsr_slot_blk_p)
    CALL dbcsr_addto_index_array (matrix, dbcsr_slot_row_p,&
         DATA=new_row_p(1:nrows+1), extra=new_nblks*2, error=error)
    CALL dbcsr_addto_index_array (matrix, dbcsr_slot_col_i,&
         DATA=new_col_i(1:new_nblks), error=error)
    IF (sort_data) THEN
       CALL dbcsr_addto_index_array (matrix, dbcsr_slot_blk_p,&
            DATA=new_blk_p_sorted, error=error)
    ELSE
       CALL dbcsr_addto_index_array (matrix, dbcsr_slot_blk_p,&
            DATA=new_blk_p, error=error)
    ENDIF
    matrix%nblks = new_nblks
    matrix%nze = new_data_size
    matrix%index(dbcsr_slot_nblks) = matrix%nblks
    matrix%index(dbcsr_slot_nze) = matrix%nze
    CALL dbcsr_repoint_index(matrix)
    !$OMP END MASTER
    !
    !$OMP BARRIER
    !
    DEALLOCATE (my_row_p)
    IF (sort_data) THEN
       !$OMP MASTER
       !$OMP CRITICAL (crit_data)
       CALL dbcsr_data_release (all_data_areas(0))
       !$OMP END CRITICAL (crit_data)
       DO which_wm = 1, nwms
          !$OMP CRITICAL (crit_data)
          CALL dbcsr_data_release (all_data_areas(which_wm))
          !$OMP END CRITICAL (crit_data)
       ENDDO
       !$OMP END MASTER
    ENDIF
    !$OMP BARRIER
    !$OMP MASTER
    DEALLOCATE (all_row_p)
    DEALLOCATE (all_col_i)
    DEALLOCATE (all_blk_p)
    DEALLOCATE (new_row_p)
    DEALLOCATE (new_col_i)
    DEALLOCATE (new_blk_p)
    DEALLOCATE (all_data_sizes)
    DEALLOCATE (all_data_offsets)
    IF (sort_data) THEN
       DEALLOCATE (all_data_areas)
       DEALLOCATE (new_blk_sizes)
       DEALLOCATE (new_blk_p_sorted)
    ENDIF
    DEALLOCATE (blk_d)
    !$OMP END MASTER
    !$OMP BARRIER
    IF (dbg) WRITE(*,*)routineN//" stopped"
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_merge_all

! *****************************************************************************
!> \brief Sort data from the WM into the final matrix, based closely on dbcsr_merge_all
!>
!> \param[in,out] matrix      matrix to work on
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE dbcsr_merge_single_wm (matrix, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_merge_single_wm', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: index_blocks = 1, &
                                                index_data = 2
    LOGICAL, PARAMETER                       :: dbg = .FALSE.

    INTEGER                                  :: error_handler, nblks, nrows
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: new_blk_p_sorted, &
                                                new_blk_sizes, new_row_p
    INTEGER, ALLOCATABLE, DIMENSION(:), &
      SAVE, TARGET                           :: blk_d
    INTEGER, DIMENSION(:), POINTER, SAVE     :: cbs, rbs
    TYPE(dbcsr_work_type), POINTER           :: wm

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    !> \par Outline

    !> There is a single work matrix.  Data blocks are sorted and copied
    !> into the matrix data area (which is empty).  The index is made consistent

    nrows = matrix%nblkrows_total
    wm => matrix%wms(1)
    nblks = wm%lastblk
    CALL dbcsr_assert (.NOT.dbcsr_wm_use_mutable (wm),&
          dbcsr_fatal_level, dbcsr_caller_error, routineN,&
          "Number of threads does not match number of work matrices.",&
          __LINE__, error=error)
    !
    ! The following is valid because the data actually referenced
    ! by blocks is explicitly calculated in dbcsr_finalize.
    !$OMP MASTER
    rbs => array_data(matrix%row_blk_size)
    cbs => array_data(matrix%col_blk_size)
    !
    ! Initializations
    !
    ALLOCATE (new_row_p(nrows+1))
    ALLOCATE (blk_d(nblks))
    ALLOCATE (new_blk_p_sorted (nblks))
    ALLOCATE (new_blk_sizes (nblks))
    !
    ! Master thread creates a row_p index for the (sorted) blocks.
    CALL dbcsr_sort_indices (nblks, wm%row_i, wm%col_i, wm%blk_p)
    CALL dbcsr_make_dbcsr_index (new_row_p, wm%row_i, nrows, nblks)
    !
    !
    ! The matrix data area is resized.  Blocks from the work
    ! matrices will be copied into it in order.
    !
    CALL dbcsr_data_ensure_size(matrix%data_area, wm%datasize, nocopy=.TRUE., error=error)
    !$OMP END MASTER
    !$OMP BARRIER
    CALL dbcsr_calc_block_sizes(new_blk_sizes,&
         new_row_p, wm%col_i, rbs, cbs)
    CALL dbcsr_sort_data (new_blk_p_sorted, wm%blk_p,&
         new_blk_sizes, dsts=matrix%data_area,&
         src=wm%data_area,&
         error=error)
    !
    ! Creates a new index array.
    !
    !$OMP BARRIER
    !$OMP MASTER
    CALL dbcsr_clearfrom_index_array (matrix, dbcsr_slot_row_p)
    CALL dbcsr_clearfrom_index_array (matrix, dbcsr_slot_col_i)
    CALL dbcsr_clearfrom_index_array (matrix, dbcsr_slot_blk_p)
    CALL dbcsr_addto_index_array (matrix, dbcsr_slot_row_p,&
         DATA=new_row_p(1:nrows+1), extra=nblks*2, error=error)
    CALL dbcsr_addto_index_array (matrix, dbcsr_slot_col_i,&
         DATA=wm%col_i(1:nblks), error=error)
    CALL dbcsr_addto_index_array (matrix, dbcsr_slot_blk_p,&
            DATA=new_blk_p_sorted, error=error)
    matrix%nblks = nblks
    matrix%nze = wm%datasize
    matrix%index(dbcsr_slot_nblks) = matrix%nblks
    matrix%index(dbcsr_slot_nze) = matrix%nze
    CALL dbcsr_repoint_index(matrix)
    DEALLOCATE (new_row_p)
    DEALLOCATE (new_blk_sizes)
    DEALLOCATE (new_blk_p_sorted)
    DEALLOCATE (blk_d)
    !$OMP END MASTER
    IF (dbg) WRITE(*,*)routineN//" stopped"
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_merge_single_wm



! *****************************************************************************
!> \brief Builds a new index from several work matrices.
!>
!> \param[in] wm          work matrix
!> \param[in,out] additions_per_row     number of blocks and data to add for
!>                                      every row
!> \param[out] additions_total          total number of blocks and data to add
!> \param[out] blk_p          block sizes of data to add
!> \param[in,out] error       error
! *****************************************************************************
  SUBROUTINE merge_index (new_row_p, new_col_i, new_blk_p,&
       blk_d, old_row_p, old_col_i, old_blk_p,&
       all_row_p, all_col_i, all_blk_p,&
       all_data_offsets, nwms, nrows, max_row_count, sort_data, error)
    INTEGER, DIMENSION(*), INTENT(IN)        :: new_row_p
    INTEGER, DIMENSION(*), INTENT(OUT), &
      TARGET                                 :: new_col_i, new_blk_p
    INTEGER, DIMENSION(*), INTENT(IN), &
      TARGET                                 :: blk_d
    INTEGER, DIMENSION(*), INTENT(IN)        :: old_row_p, old_col_i, &
                                                old_blk_p
    TYPE(i_array_p), DIMENSION(*), &
      INTENT(IN)                             :: all_row_p, all_col_i, &
                                                all_blk_p
    INTEGER, DIMENSION(*), INTENT(IN)        :: all_data_offsets
    INTEGER, INTENT(IN)                      :: nwms, nrows, max_row_count
    LOGICAL, INTENT(IN)                      :: sort_data
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'merge_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk1, blk2, error_handle, &
                                                first_row_blk, last_row_blk, &
                                                row, src_blk_1, src_blk_2, t
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: blk_p_buff, tmp_arr
    INTEGER, DIMENSION(:), POINTER           :: blk_span, col_span, d_span

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    !
    ALLOCATE (tmp_arr(max_row_count))
    ALLOCATE (blk_p_buff(max_row_count))
    !
    !$OMP DO
    DO row = 1, nrows
       first_row_blk = new_row_p(row)+1
       last_row_blk = new_row_p(row+1)
       col_span => new_col_i (first_row_blk : last_row_blk)
       blk_span => new_blk_p (first_row_blk : last_row_blk)
       IF (sort_data) d_span => blk_d (first_row_blk : last_row_blk)
       src_blk_1 = old_row_p(row)+1
       src_blk_2 = old_row_p(row+1)
       blk1 = 1
       blk2 = blk1 + (src_blk_2 - src_blk_1+1) - 1
       col_span(blk1:blk2) = old_col_i(src_blk_1:src_blk_2)
       blk_span(blk1:blk2) = old_blk_p(src_blk_1:src_blk_2)
       IF (sort_data) THEN
          d_span(blk1:blk2) = 1
          DO t = 1, nwms
             src_blk_1 = all_row_p(t)%p(row)+1
             src_blk_2 = all_row_p(t)%p(row+1)
             blk1 = blk2+1
             blk2 = blk1 + (src_blk_2 - src_blk_1 + 1) - 1
             col_span(blk1:blk2) = all_col_i(t)%p(src_blk_1:src_blk_2)
             blk_span(blk1:blk2) = all_blk_p(t)%p(src_blk_1:src_blk_2)
             d_span(blk1:blk2) = t+1
          ENDDO
       ELSE
          DO t = 1, nwms
             src_blk_1 = all_row_p(t)%p(row)+1
             src_blk_2 = all_row_p(t)%p(row+1)
             blk1 = blk2+1
             blk2 = blk1 + (src_blk_2 - src_blk_1 + 1) - 1
             col_span(blk1:blk2) = all_col_i(t)%p(src_blk_1:src_blk_2)
             blk_span(blk1:blk2) = all_blk_p(t)%p(src_blk_1:src_blk_2) &
                  &              + all_data_offsets(t)
          ENDDO
       ENDIF
       CALL sort (col_span, SIZE(col_span), tmp_arr)
       blk_p_buff(1:SIZE(blk_span)) = blk_span(:)
       blk_span(1:SIZE(blk_span)) = blk_p_buff(tmp_arr(1:SIZE(blk_span)))
    ENDDO
    !$OMP END DO
    !
    DEALLOCATE (tmp_arr)
    DEALLOCATE (blk_p_buff)
    !
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE merge_index


#include "dbcsr_work_operations_d.F"
#include "dbcsr_work_operations_z.F"
#include "dbcsr_work_operations_s.F"
#include "dbcsr_work_operations_c.F"


END MODULE dbcsr_work_operations

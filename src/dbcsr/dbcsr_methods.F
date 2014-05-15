!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Base methods on DBCSR data structures
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.95
!>
!> <b>Modification history:</b>
!> - Created
! *****************************************************************************
MODULE dbcsr_methods
  USE array_types,                     ONLY: array_data,&
                                             array_exists,&
                                             array_hold,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release,&
                                             array_size
  USE btree_I8_k_cp2d_v,               ONLY: btree_destroy_c => btree_delete,&
                                             btree_new_c => btree_new
  USE btree_I8_k_dp2d_v,               ONLY: btree_destroy_d => btree_delete,&
                                             btree_new_d => btree_new
  USE btree_I8_k_sp2d_v,               ONLY: btree_destroy_s => btree_delete,&
                                             btree_new_s => btree_new
  USE btree_I8_k_zp2d_v,               ONLY: btree_destroy_z => btree_delete,&
                                             btree_new_z => btree_new
  USE dbcsr_block_buffers,             ONLY: dbcsr_buffers_2d_needed,&
                                             dbcsr_buffers_flush,&
                                             dbcsr_buffers_init,&
                                             dbcsr_buffers_new,&
                                             dbcsr_buffers_release
  USE dbcsr_config,                    ONLY: comm_thread_load,&
                                             use_subcommunicators
  USE dbcsr_data_methods,              ONLY: dbcsr_data_get_size,&
                                             dbcsr_data_get_size_referenced,&
                                             dbcsr_data_get_type,&
                                             dbcsr_data_hold,&
                                             dbcsr_data_init,&
                                             dbcsr_data_release,&
                                             dbcsr_data_set_size_referenced
  USE dbcsr_error_handling,            ONLY: &
       dbcsr_assert, dbcsr_caller_error, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_failure_level, dbcsr_fatal_level, &
       dbcsr_internal_error, dbcsr_warning_level, dbcsr_wrong_args_error
  USE dbcsr_ptr_util,                  ONLY: memory_deallocate
  USE dbcsr_toollib,                   ONLY: sort
  USE dbcsr_types,                     ONLY: &
       dbcsr_1d_array_obj, dbcsr_2d_array_obj, dbcsr_2d_array_type, &
       dbcsr_array_type, dbcsr_block_buffer_obj, dbcsr_data_obj, &
       dbcsr_distribution_obj, dbcsr_imagedistribution_obj, &
       dbcsr_imagedistribution_type, dbcsr_memtype_type, dbcsr_mp_obj, &
       dbcsr_mutable_obj, dbcsr_obj, dbcsr_type, dbcsr_type_antihermitian, &
       dbcsr_type_antisymmetric, dbcsr_type_complex_4, dbcsr_type_complex_8, &
       dbcsr_type_hermitian, dbcsr_type_invalid, dbcsr_type_no_symmetry, &
       dbcsr_type_real_4, dbcsr_type_real_8, dbcsr_type_symmetric, &
       dbcsr_work_type
  USE kinds,                           ONLY: default_string_length,&
                                             int_8,&
                                             real_8,&
                                             sp
  USE message_passing,                 ONLY: mp_cart_create,&
                                             mp_cart_sub,&
                                             mp_comm_free,&
                                             mp_sum
  USE min_heap,                        ONLY: heap_fill,&
                                             heap_get_first,&
                                             heap_new,&
                                             heap_release,&
                                             heap_reset_first,&
                                             heap_t

  !$ USE OMP_LIB
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_methods'

  INTEGER, PUBLIC, SAVE :: dbcsr_matrix_counter = 111111

  PUBLIC :: dbcsr_init, dbcsr_init_p
  PUBLIC :: dbcsr_hold, dbcsr_release, dbcsr_release_p, dbcsr_destroy
  PUBLIC :: dbcsr_valid_index, dbcsr_is_initialized
  PUBLIC :: dbcsr_mp_new, dbcsr_mp_hold, dbcsr_mp_release
  PUBLIC :: dbcsr_mp_pgrid, dbcsr_mp_numnodes, dbcsr_mp_mynode, dbcsr_mp_group,&
            dbcsr_mp_new_transposed, dbcsr_mp_nprows, dbcsr_mp_npcols,&
            dbcsr_mp_myprow, dbcsr_mp_mypcol, dbcsr_mp_pgrid_equal,&
            dbcsr_mp_get_coordinates,&
            dbcsr_mp_my_row_group, dbcsr_mp_my_col_group,&
            dbcsr_mp_has_subgroups, dbcsr_mp_get_process,&
            dbcsr_mp_grid_setup, dbcsr_mp_grid_remove,&
            dbcsr_mp_init, dbcsr_mp_exists, dbcsr_mp_active
  PUBLIC :: dbcsr_distribution_new, dbcsr_distribution_hold,&
            dbcsr_distribution_release, dbcsr_distribution_init
  PUBLIC :: dbcsr_distribution_mp, dbcsr_distribution_processor,&
            dbcsr_distribution_nrows, dbcsr_distribution_ncols,&
            dbcsr_distribution_row_dist, dbcsr_distribution_col_dist,&
            dbcsr_distribution_nlocal_rows, dbcsr_distribution_nlocal_cols,&
            dbcsr_distribution_local_rows, dbcsr_distribution_local_cols
  PUBLIC :: dbcsr_distribution_thread_dist, dbcsr_distribution_has_threads,&
            dbcsr_distribution_make_threads, dbcsr_distribution_no_threads,&
            dbcsr_distribution_num_threads
  PUBLIC :: dbcsr_distribution_add_row_map, dbcsr_distribution_add_col_map,&
            dbcsr_distribution_del_row_map, dbcsr_distribution_del_col_map
  PUBLIC :: dbcsr_release_locals, dbcsr_dist_release_locals
  PUBLIC :: dbcsr_get_info, dbcsr_distribution,&
            dbcsr_get_matrix_type, dbcsr_get_data_type, dbcsr_get_replication_type,&
            dbcsr_row_block_sizes, dbcsr_col_block_sizes,&
            dbcsr_nblkrows_total, dbcsr_nblkcols_total, dbcsr_nfullrows_total,&
            dbcsr_nfullcols_total, dbcsr_nblkcols_local, dbcsr_nblkrows_local,&
            dbcsr_nfullrows_local, dbcsr_nfullcols_local,&
            dbcsr_max_row_size, dbcsr_max_col_size,&
            dbcsr_get_occupation,&
            dbcsr_get_index_memory_type, dbcsr_get_data_memory_type, &
            dbcsr_is_real, dbcsr_is_complex,&
            dbcsr_name, dbcsr_may_be_dense,&
            dbcsr_use_mutable, dbcsr_wm_use_mutable, dbcsr_has_symmetry,&
            dbcsr_get_data_size
  PUBLIC :: dbcsr_get_data_size_used, dbcsr_get_data_size_referenced
  PUBLIC :: dbcsr_col_block_offsets, dbcsr_row_block_offsets
  PUBLIC :: dbcsr_data_area, dbcsr_switch_data_area
  PUBLIC :: dbcsr_get_num_blocks, dbcsr_get_nze


  PUBLIC :: dbcsr_blk_row_size, dbcsr_blk_column_size,&
            dbcsr_blk_row_offset, dbcsr_blk_col_offset

  PUBLIC :: dbcsr_allocate_matrix_array, dbcsr_array_new, dbcsr_array_put,&
            dbcsr_array_get, dbcsr_array_destroy, dbcsr_array_hold,&
            dbcsr_array_release
  PUBLIC :: dbcsr_destroy_array
  PUBLIC :: dbcsr_image_dist_init, dbcsr_image_dist_hold, dbcsr_image_dist_release

  PUBLIC :: dbcsr_mutable_init, dbcsr_mutable_new, dbcsr_mutable_destroy,&
            dbcsr_mutable_release, dbcsr_mutable_hold,&
            dbcsr_mutable_instantiated

  PUBLIC :: dbcsr_modify_lock, dbcsr_modify_unlock



  INTERFACE dbcsr_init
     MODULE PROCEDURE dbcsr_init_type, dbcsr_init_obj
  END INTERFACE

  INTERFACE dbcsr_init_p
     MODULE PROCEDURE dbcsr_init_obj_p
  END INTERFACE

  INTERFACE dbcsr_valid_index
     MODULE PROCEDURE dbcsr_valid_index_type, dbcsr_valid_index_obj
  END INTERFACE

  INTERFACE dbcsr_is_initialized
     MODULE PROCEDURE dbcsr_is_initialized_type, dbcsr_is_initialized_obj
  END INTERFACE

  !INTERFACE dbcsr_uses_special_memory
  !   MODULE PROCEDURE uses_special_memory_matrix, uses_special_memory_area
  !END INTERFACE

  INTERFACE dbcsr_get_data_size
     !MODULE PROCEDURE get_data_size_area, get_data_size_matrix
     MODULE PROCEDURE get_data_size_matrix
  END INTERFACE

  INTERFACE dbcsr_get_data_size_used
     !MODULE PROCEDURE get_data_size_ref_area
     MODULE PROCEDURE dbcsr_get_data_size_used
  END INTERFACE

  INTERFACE dbcsr_get_data_size_referenced
     !MODULE PROCEDURE get_data_size_ref_area, get_data_size_ref_matrix
     MODULE PROCEDURE get_data_size_ref_matrix
  END INTERFACE

  INTERFACE dbcsr_get_data_type
     !MODULE PROCEDURE dbcsr_get_data_type_obj, data_get_data_type
     MODULE PROCEDURE dbcsr_get_data_type_obj
  END INTERFACE

  INTERFACE dbcsr_is_real
     MODULE PROCEDURE dbcsr_is_real_obj
  END INTERFACE

  INTERFACE dbcsr_is_complex
     MODULE PROCEDURE dbcsr_is_complex_obj
  END INTERFACE

  INTERFACE dbcsr_blk_row_size
     MODULE PROCEDURE dbcsr_blk_row_size_type, dbcsr_blk_row_size_obj
  END INTERFACE

  INTERFACE dbcsr_blk_column_size
     MODULE PROCEDURE dbcsr_blk_column_size_type, dbcsr_blk_column_size_obj
  END INTERFACE

  ! For the 1-D and 2-D arrays

  INTERFACE dbcsr_allocate_matrix_array
     MODULE PROCEDURE array_init_1d, array_init_2d
  END INTERFACE

  INTERFACE dbcsr_array_new
     MODULE PROCEDURE array_new_1d, array_new_2d
  END INTERFACE

  INTERFACE dbcsr_array_put
     MODULE PROCEDURE array_put_1d, array_put_2d
  END INTERFACE

  INTERFACE dbcsr_array_get
     MODULE PROCEDURE array_get_1d, array_get_2d
  END INTERFACE

  INTERFACE dbcsr_array_destroy
     MODULE PROCEDURE array_destroy_1d, array_destroy_2d
  END INTERFACE

  INTERFACE dbcsr_array_hold
     MODULE PROCEDURE array_hold_1d, array_hold_2d
  END INTERFACE

  INTERFACE dbcsr_array_release
     MODULE PROCEDURE array_release_1d, array_release_2d
  END INTERFACE

  INTERFACE dbcsr_destroy_array
     MODULE PROCEDURE dbcsr_destroy_2d_array, dbcsr_destroy_1d_array
  END INTERFACE



CONTAINS


! *****************************************************************************
!> \brief Initializes a DBCSR matrix but does not create it.
!> \param[out] matrix         uninitialized matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_init_type (matrix)
    TYPE(dbcsr_type), INTENT(OUT)            :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_type', &
      routineP = moduleN//':'//routineN

!TYPE(dbcsr_error_type), INTENT(INOUT)   :: error
!   ---------------------------------------------------------------------------

    matrix%initialized = .TRUE.
    matrix%valid = .FALSE.
    matrix%refcount = 0
    ! Nullifies all pointers.
    NULLIFY (matrix%index, matrix%row_p, matrix%col_i,&
         matrix%blk_p, matrix%thr_c, matrix%coo_l)
    CALL dbcsr_data_init (matrix%data_area)
    CALL dbcsr_distribution_init (matrix%dist)
    CALL array_nullify (matrix%row_blk_size)
    CALL array_nullify (matrix%col_blk_size)
    CALL array_nullify (matrix%row_blk_offset)
    CALL array_nullify (matrix%col_blk_offset)
    CALL array_nullify (matrix%local_rows)
    CALL array_nullify (matrix%global_rows)
    CALL array_nullify (matrix%local_cols)
    CALL array_nullify (matrix%global_cols)
    NULLIFY (matrix%wms)
  END SUBROUTINE dbcsr_init_type

! *****************************************************************************
!> \brief Initializes a DBCSR matrix but does not allocate any memory.
!> \param[out] matrix         uninitialized matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_init_obj (matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_obj', &
      routineP = moduleN//':'//routineN

!TYPE(dbcsr_error_type), INTENT(INOUT)   :: error
!   ---------------------------------------------------------------------------

    CALL dbcsr_init(matrix%m)
  END SUBROUTINE dbcsr_init_obj


! *****************************************************************************
!> \brief Initializes a DBCSR matrix but does not allocate any memory.
!> \param[out] matrix         uninitialized matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_init_obj_p (matrix, error)
    TYPE(dbcsr_obj), POINTER                 :: matrix
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_obj_p', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF(ASSOCIATED(matrix)) CALL dbcsr_release_p(matrix, error)
    ALLOCATE(matrix)
    CALL dbcsr_init(matrix%m)
  END SUBROUTINE dbcsr_init_obj_p

! *****************************************************************************
!> \brief Returns whether the index structure of the matrix is valid.
!> \param[in] matrix          verify index validity of this matrix
!> \retval valid_index        index validity
! *****************************************************************************
  PURE FUNCTION dbcsr_valid_index_type (matrix) RESULT (valid_index)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    LOGICAL                                  :: valid_index

!   ---------------------------------------------------------------------------
!valid_index = .FALSE.
!IF (ASSOCIATED (matrix%row_p)) THEN
!   valid_index = SIZE (matrix%row_p) .GT. 0
!ENDIF

    valid_index = matrix%valid
  END FUNCTION dbcsr_valid_index_type

! *****************************************************************************
!> \brief Returns whether the index structure of the matrix is valid.
!> \param[in] matrix          verify index validity of this matrix
!> \retval valid_index        index validity
! *****************************************************************************
  PURE FUNCTION dbcsr_valid_index_obj (matrix) RESULT (valid_index)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL                                  :: valid_index

!   ---------------------------------------------------------------------------

    valid_index = dbcsr_valid_index_type (matrix%m)
  END FUNCTION dbcsr_valid_index_obj

! *****************************************************************************
!> \brief Returns whether the index structure of the matrix is valid.
!> \param[in] matrix          verify index validity of this matrix
!> \retval initialized        the matrix is initialized
! *****************************************************************************
  PURE FUNCTION dbcsr_is_initialized_type (matrix) RESULT (initialized)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    LOGICAL                                  :: initialized

!   ---------------------------------------------------------------------------

    initialized = matrix%initialized
  END FUNCTION dbcsr_is_initialized_type

! *****************************************************************************
!> \brief Returns whether the index structure of the matrix is valid.
!> \param[in] matrix          verify index validity of this matrix
!> \retval initialized        the matrix is initialized
! *****************************************************************************
  PURE FUNCTION dbcsr_is_initialized_obj (matrix) RESULT (initialized)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL                                  :: initialized

!   ---------------------------------------------------------------------------

    initialized = dbcsr_is_initialized_type (matrix%m)
  END FUNCTION dbcsr_is_initialized_obj


! *****************************************************************************
!> \brief Registers another reference for a DBCSR matrix
!> \param[in,out] matrix    DBCSR matrix
! *****************************************************************************
  SUBROUTINE dbcsr_hold (matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_hold', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (matrix%m%initialized,&
         dbcsr_failure_level, dbcsr_caller_error, routineN,&
         "Matrix not initialized",__LINE__,error)
    matrix%m%refcount = matrix%m%refcount + 1
  END SUBROUTINE dbcsr_hold


! *****************************************************************************
!> \brief Releases a reference for a DBCSR matrix
!>
!> If there are no references left, the matrix is destroyed.
!> \param[in,out] set    DBCSR matrix
! *****************************************************************************
  RECURSIVE SUBROUTINE dbcsr_release (matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_release', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (dbcsr_is_initialized (matrix),&
         dbcsr_warning_level, dbcsr_caller_error, routineN,&
         "Matrix not initialized",__LINE__,error)
    IF (matrix%m%initialized) THEN
       matrix%m%refcount = matrix%m%refcount - 1
       IF (matrix%m%refcount .EQ. 0) THEN
          CALL dbcsr_destroy (matrix,error)
       ENDIF
    ENDIF
  END SUBROUTINE dbcsr_release


! *****************************************************************************
!> \brief Releases a reference for a DBCSR matrix
!>
!> If there are no references left, the matrix is destroyed.
!> \param[in,out] set    DBCSR matrix
! *****************************************************************************
  SUBROUTINE dbcsr_release_p (matrix,error)
    TYPE(dbcsr_obj), POINTER                 :: matrix
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_release_p', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF(ASSOCIATED(matrix)) THEN
       CALL dbcsr_assert (dbcsr_is_initialized (matrix),&
            dbcsr_warning_level, dbcsr_caller_error, routineN,&
            "Matrix not initialized",__LINE__,error)
       IF (matrix%m%initialized) THEN
          matrix%m%refcount = matrix%m%refcount - 1
          IF (matrix%m%refcount .EQ. 0) THEN
             CALL dbcsr_destroy (matrix, error)
             DEALLOCATE(matrix)
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE dbcsr_release_p

! *****************************************************************************
!> \brief Deallocates and destroys a matrix.
!> \param[in,out] matrix      matrix
!> \param[in,out] error       cp2k error
!> \param[in] force           (optional) force deallocation
! *****************************************************************************
  RECURSIVE SUBROUTINE dbcsr_destroy(matrix, error, force)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: force

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_destroy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle
    LOGICAL                                  :: force_all

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error=error)
    CALL dbcsr_assert (dbcsr_is_initialized (matrix),&
         dbcsr_warning_level, dbcsr_caller_error, routineN,&
         "Can not destroy uninitialized matrix object.",__LINE__,error)
    force_all = .FALSE.
    IF (PRESENT (force)) force_all = force
    IF (dbcsr_is_initialized (matrix)) THEN
       CALL dbcsr_assert (force_all .OR. matrix%m%refcount.EQ.0,&
            dbcsr_warning_level, dbcsr_caller_error,&
            routineN, "You should not destroy referenced matrix.",__LINE__,error)
       CALL dbcsr_assert (.NOT.force_all .OR. matrix%m%refcount.LE.1,&
            dbcsr_warning_level, dbcsr_caller_error,&
            routineN, "You should not destroy referenced matrix.",__LINE__,error)
       IF (force_all .OR. matrix%m%refcount .EQ. 0) THEN
          CALL dbcsr_assert (.NOT. ASSOCIATED (matrix%m%wms), dbcsr_warning_level,&
               dbcsr_caller_error, routineN, "Destroying unfinalized matrix",__LINE__,error)
          IF (dbcsr_buffers_2d_needed) THEN
             CALL dbcsr_buffers_flush (matrix%m%buffers, error=error)
             CALL dbcsr_buffers_release (matrix%m%buffers, error=error)
          ENDIF
          CALL memory_deallocate(matrix%m%index, matrix%m%index_memory_type,&
               error=error)
          CALL dbcsr_data_release (matrix%m%data_area)
          CALL array_release (matrix%m%row_blk_size)
          CALL array_release (matrix%m%col_blk_size)
          CALL array_release (matrix%m%row_blk_offset)
          CALL array_release (matrix%m%col_blk_offset)
          CALL dbcsr_distribution_release(matrix%m%dist)
          IF (matrix%m%local_indexing) THEN
             CALL dbcsr_assert (array_exists (matrix%m%local_rows),&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Global row mapping should exist", __LINE__, error=error)
             CALL dbcsr_assert (array_exists (matrix%m%global_rows),&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Global row mapping should exist", __LINE__, error=error)
             CALL dbcsr_assert (array_exists (matrix%m%local_cols),&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Global column mapping should exist", __LINE__, error=error)
             CALL dbcsr_assert (array_exists (matrix%m%global_cols),&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Global column mapping should exist", __LINE__, error=error)
          ENDIF
          CALL dbcsr_release_locals (matrix, error)
          matrix%m%valid = .FALSE.
!$        CALL OMP_DESTROY_LOCK (matrix%m%modification_lock)
          CALL dbcsr_init (matrix%m)
       ENDIF
    ENDIF
    CALL dbcsr_error_stop (error_handle, error=error)
  END SUBROUTINE dbcsr_destroy

  SUBROUTINE dbcsr_release_locals (matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (matrix%m%has_local_rows) &
         CALL array_release (matrix%m%local_rows)
    IF (matrix%m%has_global_rows) &
         CALL array_release (matrix%m%global_rows)
    IF (matrix%m%has_local_cols) &
         CALL array_release (matrix%m%local_cols)
    IF (matrix%m%has_global_cols) &
         CALL array_release (matrix%m%global_cols)
    matrix%m%has_local_rows  = .FALSE.
    matrix%m%has_global_rows = .FALSE.
    matrix%m%has_local_cols  = .FALSE.
    matrix%m%has_global_cols = .FALSE.
  END SUBROUTINE dbcsr_release_locals

!  SUBROUTINE dbcsr_release_vlocals (matrix, error)
!    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
!    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
!
!    IF (matrix%m%has_local_vrows) &
!         CALL array_release (matrix%m%local_vrows)
!    IF (matrix%m%has_global_vrows) &
!         CALL array_release (matrix%m%global_vrows)
!    IF (matrix%m%has_local_vcols) &
!         CALL array_release (matrix%m%local_vcols)
!    IF (matrix%m%has_global_vcols) &
!         CALL array_release (matrix%m%global_vcols)
!    matrix%m%has_local_vrows  = .FALSE.
!    matrix%m%has_global_vrows = .FALSE.
!    matrix%m%has_local_vcols  = .FALSE.
!    matrix%m%has_global_vcols = .FALSE.
!  END SUBROUTINE dbcsr_release_vlocals
!


! *****************************************************************************
!> \brief Initializes a new process grid
! *****************************************************************************
  SUBROUTINE dbcsr_mp_init (mp_env)
    TYPE(dbcsr_mp_obj), INTENT(OUT)          :: mp_env

    NULLIFY (mp_env%mp)
  END SUBROUTINE dbcsr_mp_init

! *****************************************************************************
!> \brief Checks whether the message passing environment exists.
! *****************************************************************************
  FUNCTION dbcsr_mp_exists (mp_env) RESULT (exists)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    LOGICAL                                  :: exists

    exists = ASSOCIATED (mp_env%mp)
  END FUNCTION dbcsr_mp_exists

! *****************************************************************************
!> \brief Checks whether this process is part of the message passing environment
! *****************************************************************************
  FUNCTION dbcsr_mp_active (mp_env) RESULT (active)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    LOGICAL                                  :: active

    active = ASSOCIATED (mp_env%mp)
  END FUNCTION dbcsr_mp_active


! *****************************************************************************
!> \brief Creates new process grid
!> \param[out] mp_env         multiprocessor environment
!> \param[in] pgrid           process grid
!> \param[in] mynode          my processor number
!> \param[in] numnodes        total number of processors (processes)
!> \par History
!>      UB (2010-02-04)       Duplicates own communicator and sets up
!>                            cartesian grid
! *****************************************************************************
  SUBROUTINE dbcsr_mp_new(mp_env, pgrid, mp_group, mynode, numnodes, myprow,&
       mypcol)
    TYPE(dbcsr_mp_obj), INTENT(OUT)          :: mp_env
    INTEGER, DIMENSION(0:, 0:), INTENT(IN)   :: pgrid
    INTEGER, INTENT(IN)                      :: mp_group, mynode
    INTEGER, INTENT(IN), OPTIONAL            :: numnodes, myprow, mypcol

    INTEGER                                  :: pcol, prow

!   ---------------------------------------------------------------------------

     ALLOCATE(mp_env%mp)
     mp_env%mp%refcount = 1
     ALLOCATE (mp_env%mp%pgrid (0:SIZE(pgrid, 1)-1, 0:SIZE(pgrid, 2)-1 ))
     mp_env%mp%pgrid(:,:) = pgrid(:,:)
     mp_env%mp%mynode = mynode
     mp_env%mp%mp_group = mp_group
     IF (PRESENT (numnodes)) THEN
        mp_env%mp%numnodes = numnodes
     ELSE
        mp_env%mp%numnodes = SIZE (pgrid)
     ENDIF
     IF (PRESENT (myprow) .AND. PRESENT (mypcol)) THEN
        mp_env%mp%myprow = myprow
        mp_env%mp%mypcol = mypcol
     ELSE
        mp_env%mp%myprow = -33777
        mp_env%mp%mypcol = -33777
        column_loop: DO pcol = LBOUND (pgrid, 2), UBOUND (pgrid, 2)
           row_loop: DO prow = LBOUND (pgrid, 1), UBOUND (pgrid, 1)
              test_position: IF (pgrid (prow, pcol) .EQ. mynode) THEN
                 mp_env%mp%myprow = prow
                 mp_env%mp%mypcol = pcol
                 EXIT column_loop
              ENDIF test_position
           ENDDO row_loop
        ENDDO column_loop
     ENDIF
     mp_env%mp%subgroups_defined = .FALSE.
!!!! KG workaround in place.
!!!! This will be the replacement:
!    ALLOCATE(mp_env%mp)
!    mp_env%mp%refcount = 1
!    ndims = 2
!    dims(1:2) = (/ SIZE (pgrid, 1), SIZE (pgrid, 2) /)
!    CALL mp_cart_create (mp_group, ndims,&
!         dims, my_pos,&
!         mp_env%mp%mp_group)
!    CALL mp_environ (mp_env%mp%numnodes, mp_env%mp%mynode, mp_env%mp%mp_group)
!    mp_env%mp%myprow = my_pos(1)
!    mp_env%mp%mypcol = my_pos(2)
!    ALLOCATE (mp_env%mp%pgrid (0:SIZE(pgrid, 1)-1, 0:SIZE(pgrid, 2)-1 ))
!    column_loop: DO pcol = 0, SIZE (mp_env%mp%pgrid, 2)-1
!       row_loop: DO prow = 0, SIZE (mp_env%mp%pgrid, 1)-1
!          CALL mp_cart_rank (mp_env%mp%mp_group, (/ prow, pcol /),&
!               mp_env%mp%pgrid (prow, pcol))
!       ENDDO row_loop
!    ENDDO column_loop
!    mp_env%mp%subgroups_defined = .FALSE.
  END SUBROUTINE dbcsr_mp_new

! *****************************************************************************
!> \brief Sets up MPI cartesian process grid
!> \param[in,out] mp_env      multiprocessor environment
!> \param[in] force           (optional) force creation of subcommunicators
! *****************************************************************************
  SUBROUTINE dbcsr_mp_grid_setup(mp_env, force)
    TYPE(dbcsr_mp_obj), INTENT(INOUT)        :: mp_env
    LOGICAL, INTENT(IN), OPTIONAL            :: force

    INTEGER                                  :: ndims, tmp_group
    INTEGER, DIMENSION(2)                    :: dims, my_pos
    LOGICAL                                  :: should_make
    LOGICAL, DIMENSION(2)                    :: remain
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------
! Intel MPI Workaround
!CALL dbcsr_assert (..NOT. mp_env%mp%subgroups_defined,&
!     dbcsr_fatal_level, dbcsr_internal_error, "mp_grid_setup",&
!     "Subcommunicators should be turned off")

    IF (PRESENT (force)) THEN
       should_make = use_subcommunicators .OR. force
    ELSE
       should_make = use_subcommunicators
    ENDIF

    IF (.NOT. mp_env%mp%subgroups_defined .AND. should_make) THEN
       ! KG workaround.
       ! This will be deleted (replaced by code in mp_new).
       ndims = 2
       dims(1:2) = (/ SIZE (mp_env%mp%pgrid, 1), SIZE (mp_env%mp%pgrid, 2) /)
       CALL mp_cart_create (mp_env%mp%mp_group, ndims,&
            dims, my_pos,&
            tmp_group)
       CALL dbcsr_assert (my_pos(1) .EQ. mp_env%mp%myprow,&
            dbcsr_fatal_level, dbcsr_internal_error, "mp_grid_setup",&
            "Got different MPI process grid",__LINE__,error)
       CALL dbcsr_assert (my_pos(2) .EQ. mp_env%mp%mypcol,&
            dbcsr_fatal_level, dbcsr_internal_error, "mp_grid_setup",&
            "Got different MPI process grid",__LINE__,error)
       !
       remain = (/.FALSE., .TRUE./)
       CALL mp_cart_sub (tmp_group, remain, mp_env%mp%prow_group)
       remain = (/.TRUE., .FALSE./)
       CALL mp_cart_sub (tmp_group, remain, mp_env%mp%pcol_group)
       CALL mp_comm_free (tmp_group)
       mp_env%mp%subgroups_defined = .TRUE.
    ENDIF
  END SUBROUTINE dbcsr_mp_grid_setup

! *****************************************************************************
!> \brief Removes an MPI cartesian process grid
!> \param[in,out] mp_env      multiprocessor environment
! *****************************************************************************
  SUBROUTINE dbcsr_mp_grid_remove (mp_env)
    TYPE(dbcsr_mp_obj), INTENT(INOUT)        :: mp_env

    IF (mp_env%mp%subgroups_defined) THEN
       CALL mp_comm_free (mp_env%mp%prow_group)
       CALL mp_comm_free (mp_env%mp%pcol_group)
    ENDIF
  END SUBROUTINE dbcsr_mp_grid_remove

! *****************************************************************************
!> \brief Locks the DBCSR data structure for modification.
!> \param[in,out] matrix      the matrix to lock
! *****************************************************************************
  SUBROUTINE dbcsr_modify_lock (matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    !$ CALL OMP_SET_LOCK (matrix%m%modification_lock)
  END SUBROUTINE dbcsr_modify_lock


! *****************************************************************************
!> \brief Releases the modification lock on a DBCSR data structure.
!> \param[in,out] matrix      the matrix to lock
! *****************************************************************************
  SUBROUTINE dbcsr_modify_unlock (matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    !$ CALL OMP_UNSET_LOCK (matrix%m%modification_lock)
  END SUBROUTINE dbcsr_modify_unlock

! *****************************************************************************
!> \brief Locks the DBCSR data structure for reading.
!> \param[in,out] matrix      the matrix to lock
! *****************************************************************************
  SUBROUTINE dbcsr_read_lock (matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    !$ CALL OMP_SET_LOCK (matrix%m%modification_lock)
  END SUBROUTINE dbcsr_read_lock

! *****************************************************************************
!> \brief Releases the reading lock on a DBCSR data structure.
!> \param[in,out] matrix      the matrix to lock
! *****************************************************************************
  SUBROUTINE dbcsr_read_unlock (matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    !$ CALL OMP_UNSET_LOCK (matrix%m%modification_lock)
  END SUBROUTINE dbcsr_read_unlock


! *****************************************************************************
!> \brief Marks another use of the mp_env
!> \param[in,out] mp_env      multiprocessor environment
! *****************************************************************************
  PURE SUBROUTINE dbcsr_mp_hold(mp_env)
    TYPE(dbcsr_mp_obj), INTENT(INOUT)        :: mp_env

!   ---------------------------------------------------------------------------

    mp_env%mp%refcount = mp_env%mp%refcount+1
  END SUBROUTINE dbcsr_mp_hold

! *****************************************************************************
!> \brief Releases and potentially destrops an mp_env
!> \param[in,out] mp_env         multiprocessor environment
! *****************************************************************************
  SUBROUTINE dbcsr_mp_release(mp_env)
    TYPE(dbcsr_mp_obj), INTENT(INOUT)        :: mp_env

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (mp_env%mp)) THEN
       mp_env%mp%refcount = mp_env%mp%refcount - 1
       IF (mp_env%mp%refcount .LE. 0) THEN
          CALL dbcsr_mp_grid_remove (mp_env)
          ! KG workaround
          !CALL mp_comm_free (mp_env%mp%mp_group)
          DEALLOCATE (mp_env%mp%pgrid)
          DEALLOCATE (mp_env%mp)
          NULLIFY (mp_env%mp)
       ENDIF
    ENDIF
  END SUBROUTINE dbcsr_mp_release

  PURE FUNCTION dbcsr_mp_get_process(mp_env, prow, pcol) RESULT (process)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, INTENT(IN)                      :: prow, pcol
    INTEGER                                  :: process

    process = mp_env%mp%pgrid(prow, pcol)
  END FUNCTION dbcsr_mp_get_process

  PURE SUBROUTINE dbcsr_mp_get_coordinates(mp_env, node, prow, pcol)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, INTENT(IN)                      :: node
    INTEGER, INTENT(OUT)                     :: prow, pcol

    INTEGER                                  :: tmp

    pcol = MOD (node, SIZE(mp_env%mp%pgrid,2))
    prow = node / SIZE(mp_env%mp%pgrid,2)
    IF (mp_env%mp%pgrid(prow, pcol) .NE. node) THEN
       tmp = pcol ; pcol = prow ; prow = tmp
    ENDIF
  END SUBROUTINE dbcsr_mp_get_coordinates

  FUNCTION dbcsr_mp_pgrid(mp_env) RESULT (pgrid)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid

    pgrid => mp_env%mp%pgrid
  END FUNCTION dbcsr_mp_pgrid
  PURE FUNCTION dbcsr_mp_numnodes(mp_env) RESULT (numnodes)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: numnodes

    numnodes = mp_env%mp%numnodes
  END FUNCTION dbcsr_mp_numnodes
  PURE FUNCTION dbcsr_mp_mynode(mp_env) RESULT (mynode)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: mynode

    mynode = mp_env%mp%mynode
  END FUNCTION dbcsr_mp_mynode
  PURE FUNCTION dbcsr_mp_group(mp_env) RESULT (mp_group)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: mp_group

    mp_group = mp_env%mp%mp_group
  END FUNCTION dbcsr_mp_group
  PURE FUNCTION dbcsr_mp_nprows(mp_env) RESULT (nprows)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: nprows

    nprows = SIZE (mp_env%mp%pgrid, 1)
  END FUNCTION dbcsr_mp_nprows
  PURE FUNCTION dbcsr_mp_npcols(mp_env) RESULT (npcols)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: npcols

    npcols = SIZE (mp_env%mp%pgrid, 2)
  END FUNCTION dbcsr_mp_npcols
  PURE FUNCTION dbcsr_mp_myprow(mp_env) RESULT (myprow)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: myprow

    myprow = mp_env%mp%myprow
  END FUNCTION dbcsr_mp_myprow
  PURE FUNCTION dbcsr_mp_mypcol(mp_env) RESULT (mypcol)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: mypcol

    mypcol = mp_env%mp%mypcol
  END FUNCTION dbcsr_mp_mypcol
  PURE FUNCTION dbcsr_mp_has_subgroups(mp_env) RESULT (has_subgroups)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    LOGICAL                                  :: has_subgroups

    has_subgroups = mp_env%mp%subgroups_defined
  END FUNCTION dbcsr_mp_has_subgroups
  PURE FUNCTION dbcsr_mp_my_row_group(mp_env) RESULT (row_group)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: row_group

    row_group = mp_env%mp%prow_group
  END FUNCTION dbcsr_mp_my_row_group
  PURE FUNCTION dbcsr_mp_my_col_group(mp_env) RESULT (col_group)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: col_group

    col_group = mp_env%mp%pcol_group
  END FUNCTION dbcsr_mp_my_col_group
  FUNCTION dbcsr_mp_pgrid_equal(mp_env1, mp_env2) RESULT (equal_pgrid)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env1, mp_env2
    LOGICAL                                  :: equal_pgrid

    INTEGER, DIMENSION(:, :), POINTER        :: pgrid1, pgrid2

    IF (dbcsr_mp_nprows (mp_env1) .EQ. dbcsr_mp_nprows (mp_env2)&
         .AND. dbcsr_mp_nprows (mp_env1) .EQ. dbcsr_mp_nprows (mp_env2)) THEN
       pgrid1 => dbcsr_mp_pgrid (mp_env1)
       pgrid2 => dbcsr_mp_pgrid (mp_env2)
       IF (ALL (pgrid1 .EQ. pgrid2)) THEN
          equal_pgrid = .TRUE.
       ELSE
          equal_pgrid = .FALSE.
       ENDIF
    ELSE
       equal_pgrid = .FALSE.
    ENDIF
  END FUNCTION dbcsr_mp_pgrid_equal

! *****************************************************************************
!> \brief Transposes a multiprocessor environment
!> \param[out] mp_t           transposed multiprocessor environment
!> \param[in] mp              original multiprocessor environment
! *****************************************************************************
  SUBROUTINE dbcsr_mp_new_transposed(mp_t, mp)
    TYPE(dbcsr_mp_obj), INTENT(OUT)          :: mp_t
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp

!   ---------------------------------------------------------------------------

    CALL dbcsr_mp_new (mp_t, TRANSPOSE (dbcsr_mp_pgrid (mp)),&
         dbcsr_mp_group (mp),&
         dbcsr_mp_mynode (mp), dbcsr_mp_numnodes (mp),&
         dbcsr_mp_mypcol (mp), dbcsr_mp_myprow (mp))
  END SUBROUTINE dbcsr_mp_new_transposed


! *****************************************************************************
!> \brief Creates new distribution
!> \param[out] dist           distribution
!> \param[in] mp_env          multiprocessing environment
!> \param[in] row_dist, col_dist        row and column distributions
! *****************************************************************************
  SUBROUTINE dbcsr_distribution_new(dist, mp_env, row_dist, col_dist,&
       local_rows, local_cols)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    TYPE(array_i1d_obj), INTENT(IN)          :: row_dist, col_dist
    TYPE(array_i1d_obj), INTENT(IN), &
      OPTIONAL                               :: local_rows, local_cols

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_distribution_new', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, mypcoor, seq
    INTEGER, DIMENSION(:), POINTER           :: dd, ld
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    ALLOCATE (dist%d)
    dist%d%refcount = 1
    dist%d%row_dist = row_dist
    CALL array_hold (dist%d%row_dist)
    dist%d%col_dist = col_dist
    CALL array_hold (dist%d%col_dist)
    dist%d%mp_env = mp_env
    CALL dbcsr_mp_hold (dist%d%mp_env)
    ! Verify given process row distribution.
    dd => array_data (row_dist)
    CALL dbcsr_assert (MAXVAL (dd), "LT", dbcsr_mp_nprows (mp_env),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "A process row is too big for process grid",&
         __LINE__, error=error)
    ! Verify given process column distribution.
    dd => array_data (col_dist)
    CALL dbcsr_assert (MAXVAL (dd), "LT", dbcsr_mp_npcols (mp_env),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "A process column is too big for process grid",&
         __LINE__, error=error)
    IF (PRESENT (local_rows)) THEN
       dist%d%local_rows = local_rows
       CALL array_hold (dist%d%local_rows)
    ELSE
       dd => array_data (row_dist)
       mypcoor = dbcsr_mp_myprow (mp_env)
       i = COUNT (dd.EQ.mypcoor)
       NULLIFY (ld)
       ALLOCATE (ld(i))
       seq = 1
       DO i = 1, array_size (row_dist)
          IF (dd(i) .EQ. mypcoor) THEN
             ld(seq) = i
             seq = seq+1
          ENDIF
       ENDDO
       CALL array_new (dist%d%local_rows, ld, gift=.TRUE.)
    ENDIF
    IF (PRESENT (local_cols)) THEN
       dist%d%local_cols = local_cols
       CALL array_hold (dist%d%local_cols)
    ELSE
       dd => array_data (col_dist)
       mypcoor = dbcsr_mp_mypcol (mp_env)
       i = COUNT (dd.EQ.mypcoor)
       NULLIFY (ld)
       ALLOCATE (ld(i))
       seq = 1
       DO i = 1, array_size (col_dist)
          IF (dd(i) .EQ. mypcoor) THEN
             ld(seq) = i
             seq = seq+1
          ENDIF
       ENDDO
       CALL array_new (dist%d%local_cols, ld, gift=.TRUE.)
    ENDIF
    dist%d%num_threads = 1
!$  dist%d%num_threads = OMP_GET_MAX_THREADS()
    dist%d%has_thread_dist = .FALSE.
    CALL array_nullify (dist%d%thread_dist)
    CALL array_nullify (dist%d%row_map)
    CALL array_nullify (dist%d%col_map)
    NULLIFY (dist%d%other_l_rows)
    NULLIFY (dist%d%other_l_cols)
    dist%d%has_other_l_rows = .FALSE.
    dist%d%has_other_l_cols = .FALSE.
    CALL array_nullify (dist%d%global_row_map)
    CALL array_nullify (dist%d%global_col_map)
    dist%d%has_global_row_map = .FALSE.
    dist%d%has_global_col_map = .FALSE.
  END SUBROUTINE dbcsr_distribution_new

! *****************************************************************************
!> \brief Marks another use of the distribution
!> \param[in,out] mp_env      multiprocessor environment
! *****************************************************************************
  SUBROUTINE dbcsr_distribution_hold(dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist

!   ---------------------------------------------------------------------------

    dist%d%refcount = dist%d%refcount + 1
  END SUBROUTINE dbcsr_distribution_hold

! *****************************************************************************
!> \brief Releases and potentially destrops a distribution
!> \param[in,out] mp_env         multiprocessor environment
! *****************************************************************************
  SUBROUTINE dbcsr_distribution_release(dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (dist%d)) THEN
       dist%d%refcount = dist%d%refcount - 1
       IF (dist%d%refcount .EQ. 0) THEN
          CALL array_release (dist%d%row_dist)
          CALL array_release (dist%d%col_dist)
          CALL array_release (dist%d%local_rows)
          CALL array_release (dist%d%local_cols)
          CALL dbcsr_mp_release (dist%d%mp_env)
          IF (dist%d%has_thread_dist) &
               CALL array_release (dist%d%thread_dist)
          CALL array_release (dist%d%row_map)
          CALL array_release (dist%d%col_map)
          CALL dbcsr_dist_release_locals (dist, error)
          DEALLOCATE (dist%d)
          CALL dbcsr_distribution_init (dist)
       ENDIF
    ENDIF
  END SUBROUTINE dbcsr_distribution_release

  SUBROUTINE dbcsr_dist_release_locals (dist, error)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: i

    IF (dist%d%has_other_l_rows) THEN
       DO i = LBOUND (dist%d%other_l_rows,1), UBOUND (dist%d%other_l_rows,1)
          CALL array_release (dist%d%other_l_rows(i))
       END DO
       DEALLOCATE (dist%d%other_l_rows)
       NULLIFY (dist%d%other_l_rows)
    ENDIF
    IF (dist%d%has_other_l_cols) THEN
       DO i = LBOUND (dist%d%other_l_cols,1), UBOUND (dist%d%other_l_cols,1)
          CALL array_release (dist%d%other_l_cols(i))
       END DO
       DEALLOCATE (dist%d%other_l_cols)
       NULLIFY (dist%d%other_l_cols)
    ENDIF
    IF (dist%d%has_global_row_map) THEN
       CALL array_release (dist%d%global_row_map)
    ENDIF
    IF (dist%d%has_global_col_map) THEN
       CALL array_release (dist%d%global_col_map)
    ENDIF
    dist%d%has_other_l_rows = .FALSE.
    dist%d%has_other_l_cols = .FALSE.
    dist%d%has_global_row_map = .FALSE.
    dist%d%has_global_col_map = .FALSE.
  END SUBROUTINE dbcsr_dist_release_locals



  SUBROUTINE dbcsr_distribution_init (dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist

    NULLIFY (dist%d)
  END SUBROUTINE dbcsr_distribution_init

  SUBROUTINE dbcsr_distribution_add_row_map (dist, row_map)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist
    TYPE(array_i1d_obj), INTENT(IN)          :: row_map

    dist%d%row_map = row_map
    CALL array_hold (dist%d%row_map)
  END SUBROUTINE dbcsr_distribution_add_row_map
  SUBROUTINE dbcsr_distribution_add_col_map (dist, col_map)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist
    TYPE(array_i1d_obj), INTENT(IN)          :: col_map

    dist%d%col_map = col_map
    CALL array_hold (dist%d%col_map)
  END SUBROUTINE dbcsr_distribution_add_col_map
  SUBROUTINE dbcsr_distribution_del_row_map (dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist

    CALL array_release (dist%d%row_map)
    CALL array_nullify (dist%d%row_map)
  END SUBROUTINE dbcsr_distribution_del_row_map
  SUBROUTINE dbcsr_distribution_del_col_map (dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist

    CALL array_release (dist%d%col_map)
    CALL array_nullify (dist%d%col_map)
  END SUBROUTINE dbcsr_distribution_del_col_map

  FUNCTION dbcsr_distribution_mp(dist) RESULT (mp_env)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(dbcsr_mp_obj)                       :: mp_env

!   ---------------------------------------------------------------------------

    mp_env = dist%d%mp_env
  END FUNCTION dbcsr_distribution_mp
  PURE FUNCTION dbcsr_distribution_nrows(dist) RESULT (nrows)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER                                  :: nrows

    nrows = array_size (dist%d%row_dist)
  END FUNCTION dbcsr_distribution_nrows
  PURE FUNCTION dbcsr_distribution_ncols(dist) RESULT (ncols)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER                                  :: ncols

    ncols = array_size (dist%d%col_dist)
  END FUNCTION dbcsr_distribution_ncols
  FUNCTION dbcsr_distribution_row_dist(dist) RESULT (row_dist)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: row_dist

!   ---------------------------------------------------------------------------

    row_dist = dist%d%row_dist
  END FUNCTION dbcsr_distribution_row_dist
  FUNCTION dbcsr_distribution_col_dist(dist) RESULT (col_dist)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: col_dist

!   ---------------------------------------------------------------------------

    col_dist = dist%d%col_dist
  END FUNCTION dbcsr_distribution_col_dist

  PURE FUNCTION dbcsr_distribution_nlocal_rows(dist) RESULT (nlocalrows)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER                                  :: nlocalrows

    nlocalrows = array_size (dist%d%local_rows)
  END FUNCTION dbcsr_distribution_nlocal_rows
  PURE FUNCTION dbcsr_distribution_nlocal_cols(dist) RESULT (nlocalcols)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER                                  :: nlocalcols

    nlocalcols = array_size (dist%d%local_cols)
  END FUNCTION dbcsr_distribution_nlocal_cols
  FUNCTION dbcsr_distribution_local_rows(dist) RESULT (local_rows)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: local_rows

    local_rows = dist%d%local_rows
  END FUNCTION dbcsr_distribution_local_rows
  FUNCTION dbcsr_distribution_local_cols(dist) RESULT (local_cols)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: local_cols

    local_cols = dist%d%local_cols
  END FUNCTION dbcsr_distribution_local_cols
  !
  PURE FUNCTION dbcsr_distribution_processor(dist, row, col)&
       RESULT (processor)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER, INTENT(IN)                      :: row, col
    INTEGER                                  :: processor

    INTEGER                                  :: c, r

    IF (ASSOCIATED (dist%d%row_map%low)) THEN ! instead of array_exists
       r = dist%d%row_map%low%data(row)
    ELSE
       r = row
    ENDIF
    IF (ASSOCIATED (dist%d%col_map%low)) THEN ! instead of array_exists
       c = dist%d%col_map%low%data(col)
    ELSE
       c = col
    ENDIF
    processor = dist%d%mp_env%mp%pgrid(dist%d%row_dist%low%data(r),&
         dist%d%col_dist%low%data(c))
  END FUNCTION dbcsr_distribution_processor

  FUNCTION dbcsr_distribution_thread_dist(dist) RESULT (thread_dist)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: thread_dist

!   ---------------------------------------------------------------------------

    thread_dist = dist%d%thread_dist
  END FUNCTION dbcsr_distribution_thread_dist


  PURE FUNCTION dbcsr_distribution_has_threads(dist) RESULT (has_thread_dist)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    LOGICAL                                  :: has_thread_dist

!   ---------------------------------------------------------------------------

    has_thread_dist = dist%d%has_thread_dist
  END FUNCTION dbcsr_distribution_has_threads

  PURE FUNCTION dbcsr_distribution_num_threads(dist) RESULT (num_threads)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER                                  :: num_threads

!   ---------------------------------------------------------------------------

    num_threads = dist%d%num_threads
  END FUNCTION dbcsr_distribution_num_threads

! *****************************************************************************
!> \brief Creates a distribution for threads
!> \param[in,out] dist   Add thread distribution to this distribution
!> \param[in] row_sizes  (optional) row block sizes
! *****************************************************************************
  SUBROUTINE dbcsr_distribution_make_threads(dist, row_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT), TARGET                  :: dist
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: row_sizes

    TYPE(dbcsr_distribution_obj), POINTER    :: dist_p

!   ---------------------------------------------------------------------------

    dist_p => dist
    !$ IF (.NOT. OMP_IN_PARALLEL ()) THEN
       !$OMP PARALLEL DEFAULT(NONE) &
       !$OMP          SHARED(dist_p,row_sizes)
       !$    CALL make_threads (dist_p, row_sizes=row_sizes)
       !$OMP END PARALLEL
    !$ ELSE
       CALL make_threads (dist_p, row_sizes=row_sizes)
       !$OMP BARRIER
    !$ ENDIF
  END SUBROUTINE dbcsr_distribution_make_threads

! *****************************************************************************
!> \brief Creates a distribution for threads
!> \par Presence of row_sizes
!>      When row_sizes is present then the thread distribution
!>      attempts to distribute rows to threads such that the sum of
!>      delegated row sizes is approximately matched for all rows.
!>
!>      When row_sizes is not present then a random distribution is chosen.
!> \param[in,out] dist   Add thread distribution to this distribution
!> \param[in] row_sizes  (optional) row block sizes
! *****************************************************************************
  SUBROUTINE make_threads(dist, row_sizes)
    TYPE(dbcsr_distribution_obj), POINTER    :: dist
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: row_sizes

    CHARACTER(len=*), PARAMETER :: routineN = 'make_threads', &
      routineP = moduleN//':'//routineN

    INTEGER :: block_size, block_size0, cur_block, group_size, i, last_row, &
      nlrows, nrows, nthreads, row, t, t_cnt
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: itemp1, itemp2, reorder, &
                                                sorted_row_sizes
    INTEGER, DIMENSION(:), POINTER           :: lrows, td
    LOGICAL                                  :: assigned, found, heap_error
    REAL(kind=sp)                            :: load_fraction, rn, soft_thr

!$  TYPE(dbcsr_error_type)                   :: error
    TYPE(heap_t)                             :: t_heap

!   ---------------------------------------------------------------------------

    nthreads = 1
!$  nthreads = OMP_GET_NUM_THREADS () ;

!$  CALL dbcsr_assert (dist%d%num_threads, "EQ", nthreads,&
!$                     dbcsr_fatal_level, dbcsr_internal_error, routineN,&
!$                     "Thread number has changed", __LINE__, error=error)
    nrows = dbcsr_distribution_nrows (dist)
    nlrows = dbcsr_distribution_nlocal_rows (dist)
    lrows => array_data (dbcsr_distribution_local_rows (dist))

    !$OMP BARRIER
    !$OMP MASTER


    load_fraction = REAL(comm_thread_load)/100.0
    IF(nthreads==1) load_fraction = 1.0

    IF (.NOT. dist%d%has_thread_dist) THEN
       dist%d%num_threads = nthreads
       group_size = 0 ; cur_block = 0

       ALLOCATE (td(nrows))
       dist%d%has_thread_dist = .TRUE.
       CALL array_new (dist%d%thread_dist, td, gift=.TRUE.)
       td => array_data (dist%d%thread_dist)

       IF (PRESENT (row_sizes)) THEN
          ! The goal is to distribute rows to threads as equally as
          ! possible. The row sizes are first sorted. Each group of
          ! equally sized rows (group_size rows of size cur_block) is
          ! distributed to threads (keeping consecutive rows
          ! together). The group is divided into equally-sized blocks
          ! (block_size0, block_size).  Leftover rows (those that can
          ! not be equally distributed to threads) are then assigned
          ! to threads so that each thread's commulative load attempts
          ! to be equal. This distribution is achieved using a heap.
          !
          ! The heap is used to distribute "leftover"rows to threads.
          ! Lefotver rows are those of the same size that can not be
          ! evenly distributed among all threads.
          CALL heap_new (t_heap, nthreads-1, heap_error)
          ! We do not want thread 0 to be in the heap.
          ALLOCATE(itemp1(1:nthreads-1))
          ALLOCATE(itemp2(1:nthreads-1))
          DO i=1,nthreads-1
             itemp1(i)=i
             itemp2(i)=0
          ENDDO
          CALL heap_fill (t_heap,itemp1,itemp2,heap_error)
          DEALLOCATE(itemp1,itemp2)
          ALLOCATE (sorted_row_sizes (nrows))
          ALLOCATE (reorder (nrows))
          sorted_row_sizes(:) = row_sizes(:)
          CALL sort (sorted_row_sizes, nrows, reorder)

          row = 1
          DO WHILE ( row .LE. nrows)
             cur_block = sorted_row_sizes(nrows-row+1)
             assigned = .FALSE.
             group_size = 0

             last_row = nrows-row+1
             DO i = last_row, 1, -1
                IF ( cur_block == sorted_row_sizes(i) ) THEN
                   group_size = group_size + 1
                   row = row + 1
                ELSE
                   EXIT
                END IF
             END DO

             soft_thr = load_fraction + nthreads - 1
             block_size0 = load_fraction*(group_size/soft_thr)
             block_size = group_size/soft_thr

             !blocks for master thread
             IF (block_size0 > 0) &
                td(reorder(last_row:last_row-block_size0+1:-1)) = 0

             !Other threads
             IF (block_size > 0) THEN
                DO t=1, nthreads-1
                   td(reorder(last_row-block_size0-(t-1)*block_size:&
                        last_row-block_size0-(t)*block_size+1:-1)) = t
                END DO
             END IF

             !Leftover bocks
             DO i=last_row-block_size0-(nthreads-1)*block_size, last_row+1-group_size, -1
                CALL heap_get_first (t_heap, t, t_cnt, found,heap_error)
                t_cnt = t_cnt + cur_block
                CALL heap_reset_first (t_heap, t_cnt, heap_error)
                td(reorder(i)) = t
             END DO

          END DO
          CALL heap_release (t_heap, heap_error)
          DEALLOCATE (sorted_row_sizes)
          DEALLOCATE (reorder)
       ELSE
          DO t = 1, nrows
             IF (.FALSE.) THEN
                td(t) = MOD(t-1, nthreads)
             ELSE
                CALL RANDOM_NUMBER (rn)
                ! Makes sure the numbers are in the proper integer range.
                td(t) = MOD (INT (rn*REAL(nthreads)), nthreads)
             ENDIF
          END DO
       ENDIF
    ENDIF
    !$OMP END MASTER
  END SUBROUTINE make_threads

!> \brief Removes the thread distribution from a distribution
  SUBROUTINE dbcsr_distribution_no_threads(dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist

!$OMP MASTER
    CALL array_release (dist%d%thread_dist)
    dist%d%has_thread_dist = .FALSE.
!$OMP END MASTER
  END SUBROUTINE dbcsr_distribution_no_threads


! Pertaining to the dbcsr matrix.

  FUNCTION dbcsr_nblkrows_total(matrix) RESULT (nblkrows_total)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nblkrows_total

    nblkrows_total = matrix%m%nblkrows_total
  END FUNCTION dbcsr_nblkrows_total

  FUNCTION dbcsr_nblkcols_total(matrix) RESULT (nblkcols_total)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nblkcols_total

    nblkcols_total = matrix%m%nblkcols_total
  END FUNCTION dbcsr_nblkcols_total
  FUNCTION dbcsr_nfullrows_total(matrix) RESULT (nfullrows_total)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nfullrows_total

    nfullrows_total = matrix%m%nfullrows_total
  END FUNCTION dbcsr_nfullrows_total
  FUNCTION dbcsr_nfullcols_total(matrix) RESULT (nfullcols_total)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nfullcols_total

    nfullcols_total = matrix%m%nfullcols_total
  END FUNCTION dbcsr_nfullcols_total
  FUNCTION dbcsr_nblkrows_local(matrix) RESULT (nblkrows_local)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nblkrows_local

    nblkrows_local = matrix%m%nblkrows_local
  END FUNCTION dbcsr_nblkrows_local
  FUNCTION dbcsr_nblkcols_local(matrix) RESULT (nblkcols_local)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nblkcols_local

    nblkcols_local = matrix%m%nblkcols_local
  END FUNCTION dbcsr_nblkcols_local
  FUNCTION dbcsr_nfullrows_local(matrix) RESULT (nfullrows_local)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nfullrows_local

    nfullrows_local = matrix%m%nfullrows_local
  END FUNCTION dbcsr_nfullrows_local
  FUNCTION dbcsr_nfullcols_local(matrix) RESULT (nfullcols_local)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nfullcols_local

    nfullcols_local = matrix%m%nfullcols_local
  END FUNCTION dbcsr_nfullcols_local
  FUNCTION dbcsr_max_row_size(matrix) RESULT (max_row_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: max_row_size

    max_row_size = matrix%m%max_rbs
  END FUNCTION dbcsr_max_row_size
  FUNCTION dbcsr_max_col_size(matrix) RESULT (max_col_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: max_col_size

    max_col_size = matrix%m%max_cbs
  END FUNCTION dbcsr_max_col_size

  FUNCTION dbcsr_distribution (matrix) RESULT (distribution)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_distribution_obj)             :: distribution

    distribution = matrix%m%dist
  END FUNCTION dbcsr_distribution

  FUNCTION dbcsr_name (matrix) RESULT (name)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    CHARACTER(len=default_string_length)     :: name

    name = matrix%m%name
  END FUNCTION dbcsr_name

! *****************************************************************************
!> \brief Returns whether this work matrix uses the mutable type
!> \param[in] wm              work matrix
!> \param[out] use_mutable    use the mutable and not append-only working
!>                            structures
! *****************************************************************************
  PURE FUNCTION dbcsr_wm_use_mutable (wm) RESULT (use_mutable)
    TYPE(dbcsr_work_type), INTENT(IN)        :: wm
    LOGICAL                                  :: use_mutable

!   ---------------------------------------------------------------------------

    use_mutable = dbcsr_mutable_instantiated (wm%mutable)
  END FUNCTION dbcsr_wm_use_mutable

! *****************************************************************************
!> \brief Returns whether work matrices should use the mutable data type
!> \param[in] matrix          matrix
!> \param[out] use_mutable    use the mutable and not append-only working
!>                            structures
! *****************************************************************************
  PURE FUNCTION dbcsr_use_mutable (matrix) RESULT (use_mutable)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    LOGICAL                                  :: use_mutable

!   ---------------------------------------------------------------------------

    use_mutable = matrix%work_mutable
  END FUNCTION dbcsr_use_mutable


  FUNCTION dbcsr_row_block_sizes (matrix) RESULT (row_blk_sizes)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(array_i1d_obj)                      :: row_blk_sizes

    row_blk_sizes = matrix%m%row_blk_size
  END FUNCTION dbcsr_row_block_sizes

  FUNCTION dbcsr_col_block_sizes (matrix) RESULT (col_blk_sizes)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(array_i1d_obj)                      :: col_blk_sizes

    col_blk_sizes = matrix%m%col_blk_size
  END FUNCTION dbcsr_col_block_sizes

  FUNCTION dbcsr_col_block_offsets (matrix) RESULT (col_blk_offsets)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(array_i1d_obj)                      :: col_blk_offsets

    col_blk_offsets = matrix%m%col_blk_offset
  END FUNCTION dbcsr_col_block_offsets

  FUNCTION dbcsr_row_block_offsets (matrix) RESULT (row_blk_offsets)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(array_i1d_obj)                      :: row_blk_offsets

    row_blk_offsets = matrix%m%row_blk_offset
  END FUNCTION dbcsr_row_block_offsets

! *****************************************************************************
!> \brief Gets information about a matrix
!> \param[in] matrix          matrix to query
!> \param[out] nblkrows_total
!> \param[out] nblkcols_total
!> \param[out] nfullrows_total
!> \param[out] nfullcols_total
!> \param[out] nblkrows_local
!> \param[out] nblkrows_local
!> \param[out] nfullrows_local
!> \param[out] nfullrows_local
!> \param[out] my_prow, my_pcol
!> \param[out] local_rows, local_cols
!> \param[out] proc_row_dist, proc_col_dist
!> \param[out] row_blk_size
!> \param[out] col_blk_size
!> \param[out] row_blk_offset
!> \param[out] col_blk_offset
!> \param[out] distribution   the data distribution of the matrix
!> \param[out] name           matrix name
!> \param[out] data_area      data_area
!> \param[out] matrix_type    matrix type (regular, symmetric, see
!>                            dbcsr_types.F for values)
!> \param[out] data_type      data type (single/double precision real/complex)
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_get_info(matrix, nblkrows_total, nblkcols_total,&
       nfullrows_total, nfullcols_total,&
       nblkrows_local, nblkcols_local,&
       nfullrows_local, nfullcols_local,&
       my_prow, my_pcol,&
       local_rows, local_cols, proc_row_dist, proc_col_dist,&
       row_blk_size, col_blk_size, row_blk_offset, col_blk_offset, distribution, name, data_area,&
       matrix_type, data_type)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(OUT), OPTIONAL :: nblkrows_total, nblkcols_total, &
      nfullrows_total, nfullcols_total, nblkrows_local, nblkcols_local, &
      nfullrows_local, nfullcols_local, my_prow, my_pcol
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: local_rows, local_cols, &
                                                proc_row_dist, proc_col_dist
    TYPE(array_i1d_obj), INTENT(OUT), &
      OPTIONAL                               :: row_blk_size, col_blk_size, &
                                                row_blk_offset, col_blk_offset
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT), OPTIONAL                  :: distribution
    CHARACTER(len=*), INTENT(OUT), OPTIONAL  :: name
    TYPE(dbcsr_data_obj), INTENT(OUT), &
      OPTIONAL                               :: data_area
    CHARACTER, OPTIONAL                      :: matrix_type
    INTEGER, OPTIONAL                        :: data_type

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_info', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (matrix%m%initialized,&
         dbcsr_fatal_level, dbcsr_caller_error, routineP,&
         "Matrix not initialized",__LINE__,error)
    !vw avoid massive priting of warnings
    !CALL cp_assert (matrix%m%valid, cp_warning_level, cp_caller_error,&
    !     routineP,"Invalid matrix")
    IF (PRESENT (nblkrows_total)) nblkrows_total = matrix%m%nblkrows_total
    IF (PRESENT (nblkcols_total)) nblkcols_total = matrix%m%nblkcols_total
    IF (PRESENT (nfullrows_total)) nfullrows_total = matrix%m%nfullrows_total
    IF (PRESENT (nfullcols_total)) nfullcols_total = matrix%m%nfullcols_total
    IF (PRESENT (nblkrows_local)) nblkrows_local = matrix%m%nblkrows_local
    IF (PRESENT (nblkcols_local)) nblkcols_local = matrix%m%nblkcols_local
    IF (PRESENT (nfullrows_local)) nfullrows_local = matrix%m%nfullrows_local
    IF (PRESENT (nfullcols_local)) nfullcols_local = matrix%m%nfullcols_local
    IF (PRESENT (row_blk_size)) row_blk_size = matrix%m%row_blk_size
    IF (PRESENT (col_blk_size)) col_blk_size = matrix%m%col_blk_size
    IF (PRESENT (row_blk_offset)) row_blk_offset = matrix%m%row_blk_offset
    IF (PRESENT (col_blk_offset)) col_blk_offset = matrix%m%col_blk_offset
    IF (PRESENT (distribution)) distribution = matrix%m%dist
    IF (PRESENT (name)) name = matrix%m%name
    IF (PRESENT (data_area)) data_area = matrix%m%data_area
    IF (PRESENT (matrix_type)) THEN
       matrix_type = dbcsr_get_matrix_type (matrix)
       CALL dbcsr_assert (matrix_type.NE.dbcsr_type_invalid, dbcsr_failure_level, dbcsr_caller_error,&
            routineN, "Incorrect symmetry",__LINE__,error)
    ENDIF
    IF (PRESENT (data_type)) data_type = matrix%m%data_type
    IF (PRESENT (local_rows)) &
         local_rows => array_data (dbcsr_distribution_local_rows (matrix%m%dist))
    IF (PRESENT (local_cols)) &
         local_cols => array_data (dbcsr_distribution_local_cols (matrix%m%dist))
    IF (PRESENT (proc_row_dist)) &
         proc_row_dist => array_data (dbcsr_distribution_row_dist (matrix%m%dist))
    IF (PRESENT (proc_col_dist)) &
         proc_col_dist => array_data (dbcsr_distribution_col_dist (matrix%m%dist))
    IF (PRESENT (my_prow)) &
       my_prow = dbcsr_mp_myprow (dbcsr_distribution_mp (matrix%m%dist))
    IF (PRESENT (my_pcol)) &
       my_pcol = dbcsr_mp_mypcol (dbcsr_distribution_mp (matrix%m%dist))
  END SUBROUTINE dbcsr_get_info

! *****************************************************************************
!> \brief Returns whether the matrix could be represeneted in a dense form
!> \param[in] matrix          matrix
!> \result may_be_dense    use the mutable and not append-only working
!>                            structures
! *****************************************************************************
  FUNCTION dbcsr_may_be_dense (matrix, occ_thresh) RESULT (may_be_dense)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    REAL(real_8), INTENT(in)                 :: occ_thresh
    LOGICAL                                  :: may_be_dense

    REAL(real_8)                             :: occ

!   ---------------------------------------------------------------------------

    occ = dbcsr_get_occupation(matrix)
    may_be_dense =  .NOT.( occ .LT. occ_thresh )
    ! make sure every proc sees the same
    CALL mp_sum(may_be_dense, dbcsr_mp_group (dbcsr_distribution_mp (matrix%m%dist)))
  END FUNCTION dbcsr_may_be_dense

! *****************************************************************************
!> \brief Returns the blocked row size of a row
!>
!> This routine is optimized for speed and no checks are performed.
!> \param[in] matrix          DBCSR matrix
!> \param[in] row             row number
!> \retval row_size           blocked row size
! *****************************************************************************
  PURE FUNCTION dbcsr_blk_row_size_type (matrix, row) RESULT (row_size)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(IN)                      :: row
    INTEGER                                  :: row_size

    row_size = matrix%row_blk_size%low%data(row)
  END FUNCTION dbcsr_blk_row_size_type
  PURE FUNCTION dbcsr_blk_row_size_obj (matrix, row) RESULT (row_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row
    INTEGER                                  :: row_size

    row_size = matrix%m%row_blk_size%low%data(row)
  END FUNCTION dbcsr_blk_row_size_obj
  PURE FUNCTION dbcsr_blk_row_offset (matrix, row) RESULT (row_offset)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row
    INTEGER                                  :: row_offset

    row_offset = matrix%m%row_blk_offset%low%data(row)
  END FUNCTION dbcsr_blk_row_offset


! *****************************************************************************
!> \brief Returns the blocked column size of a column
!>
!> This routine is optimized for speed and no checks are performed.
!> \param[in] matrix          DBCSR matrix
!> \param[in] column          column number
!> \retval column_size        blocked row size
! *****************************************************************************
  PURE FUNCTION dbcsr_blk_column_size_type (matrix, column) RESULT (column_size)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    INTEGER, INTENT(IN)                      :: column
    INTEGER                                  :: column_size

    column_size = matrix%col_blk_size%low%data(column)
  END FUNCTION dbcsr_blk_column_size_type
  PURE FUNCTION dbcsr_blk_column_size_obj (matrix, column) RESULT (column_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: column
    INTEGER                                  :: column_size

    column_size = matrix%m%col_blk_size%low%data(column)
  END FUNCTION dbcsr_blk_column_size_obj
  PURE FUNCTION dbcsr_blk_col_offset (matrix, col) RESULT (col_offset)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: col
    INTEGER                                  :: col_offset

    col_offset = matrix%m%col_blk_offset%low%data(col)
  END FUNCTION dbcsr_blk_col_offset

! *****************************************************************************
!> \brief Returns the data area
!> \param matrix     matrix from which to get data
!> \param area       data area
! *****************************************************************************
  FUNCTION dbcsr_data_area (matrix) RESULT (data_area)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    TYPE(dbcsr_data_obj)                     :: data_area

    data_area = matrix%data_area
  END FUNCTION dbcsr_data_area

! *****************************************************************************
!> \brief Sets the data area of a matrix
!> \param[in,out] matrix     matrix for which to set the data area
!> \param[in] data_area      data area to set
!> \param[out] previous_data_area  (optional) previous data area
!> \param[in,out] error      error
! *****************************************************************************
  SUBROUTINE dbcsr_switch_data_area (matrix, data_area, previous_data_area, &
       error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_data_obj), INTENT(IN)         :: data_area
    TYPE(dbcsr_data_obj), INTENT(OUT), &
      OPTIONAL                               :: previous_data_area
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_switch_data_area', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    TYPE(dbcsr_block_buffer_obj)             :: buffers

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    IF (dbcsr_buffers_2d_needed) THEN
       buffers = matrix%m%buffers
       CALL dbcsr_buffers_flush (buffers, error=error)
       CALL dbcsr_buffers_release (buffers, error=error)
       CALL dbcsr_buffers_init (buffers)
       CALL dbcsr_buffers_new (buffers, data_area, error=error)
       matrix%m%buffers = buffers
    ENDIF
    IF (PRESENT (previous_data_area)) THEN
       previous_data_area = matrix%m%data_area
    ELSE
       CALL dbcsr_data_release (matrix%m%data_area)
    ENDIF
    matrix%m%data_area = data_area
    CALL dbcsr_data_hold (matrix%m%data_area)
    CALL dbcsr_error_stop (error_handler, error)
  END SUBROUTINE dbcsr_switch_data_area


! *****************************************************************************
!> \brief Returns the matrix type
!> \param matrix              query this matrix
!> \param matrix_type         matrix_type (see dbcsr_types.F for possible
!>                            values)
! *****************************************************************************
  PURE FUNCTION dbcsr_get_matrix_type (matrix) RESULT (matrix_type)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    CHARACTER                                :: matrix_type

    matrix_type = dbcsr_type_invalid
    IF (matrix%m%symmetry) THEN
       IF ( (.NOT.matrix%m%negate_real).AND. matrix%m%negate_imaginary ) THEN
          matrix_type = dbcsr_type_hermitian
       ELSEIF ( matrix%m%negate_real .AND. (.NOT.matrix%m%negate_imaginary) ) THEN
          matrix_type = dbcsr_type_antihermitian
       ELSEIF (matrix%m%negate_real .AND. matrix%m%negate_imaginary) THEN
          matrix_type = dbcsr_type_antisymmetric
       ELSEIF ( (.NOT.matrix%m%negate_real) .AND. (.NOT.matrix%m%negate_imaginary)) THEN
          matrix_type = dbcsr_type_symmetric
       ENDIF
    ELSE
       matrix_type = dbcsr_type_no_symmetry
    ENDIF
  END FUNCTION dbcsr_get_matrix_type

! *****************************************************************************
!> \brief Whether matrix has symmetry
!> \param matrix              query this matrix
!> \result has_symmetry       matrix has symmetry
! *****************************************************************************
  PURE FUNCTION dbcsr_has_symmetry (matrix) RESULT (has_symmetry)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL                                  :: has_symmetry

    has_symmetry = matrix%m%symmetry
  END FUNCTION dbcsr_has_symmetry

! *****************************************************************************
!> \brief Whether matrix is (logically) transposed
!> \param matrix          query this matrix
!> \result transposed     matrix is transposed
! *****************************************************************************
  !PURE FUNCTION dbcsr_is_transposed (matrix) RESULT (transposed)
  !  TYPE(dbcsr_obj), INTENT(IN)              :: matrix
  !  LOGICAL                                  :: transposed
  !
  !  transposed = matrix%m%transpose
  !END FUNCTION dbcsr_is_transposed

! *****************************************************************************
!> \brief Returns the data type stored in the matrix
!> \param matrix              query this matrix
!> \param data_type           repl_type (see dbcsr_types.F for possible
!>                            values)
! *****************************************************************************
  PURE FUNCTION dbcsr_get_replication_type (matrix) RESULT (repl_type)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    CHARACTER                                :: repl_type

    repl_type = matrix%m%replication_type
  END FUNCTION dbcsr_get_replication_type

! *****************************************************************************
!> \brief Returns the data type stored in the matrix
!> \param matrix              query this matrix
!> \param data_type           data_type (see dbcsr_types.F for possible
!>                            values)
! *****************************************************************************
  PURE FUNCTION dbcsr_get_data_type_obj (matrix) RESULT (data_type)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: data_type

    data_type = matrix%m%data_type
  END FUNCTION dbcsr_get_data_type_obj

! *****************************************************************************
!> \brief Returns data type of a data area
!> \param[in] area         data area
!> \result data_type       data type of the data area
! *****************************************************************************
  PURE FUNCTION data_get_data_type (area) RESULT (data_type)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    INTEGER                                  :: data_type

    data_type = dbcsr_data_get_type (area)
  END FUNCTION data_get_data_type


! *****************************************************************************
!> \brief Returns true if the data type of the matrix is real
!> \param[in]  matrix              query this matrix
!> \param[out] is_real
! *****************************************************************************
  PURE FUNCTION dbcsr_is_real_obj (matrix) RESULT (is_real)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL                                  :: is_real

    IF(matrix%m%data_type.EQ.dbcsr_type_real_4.OR.&
       matrix%m%data_type.EQ.dbcsr_type_real_8) THEN
       is_real = .TRUE.
    ELSE
       is_real = .FALSE.
    ENDIF
  END FUNCTION dbcsr_is_real_obj

! *****************************************************************************
!> \brief Returns true if the data type of the matrix is complex
!> \param[in]  matrix              query this matrix
!> \param[out] is_complex
! *****************************************************************************
  PURE FUNCTION dbcsr_is_complex_obj (matrix) RESULT (is_complex)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL                                  :: is_complex

    IF(matrix%m%data_type.EQ.dbcsr_type_complex_4.OR.&
       matrix%m%data_type.EQ.dbcsr_type_complex_8) THEN
       is_complex = .TRUE.
    ELSE
       is_complex = .FALSE.
    ENDIF
  END FUNCTION dbcsr_is_complex_obj


! *****************************************************************************
!> \brief Returns the type of memory used for data in the matrix
!> \note It returns the declared data type, not the actually used type
!> \param[in] matrix           query this matrix
!> \param[out] same_as_actual  (optional) whether the data actually uses the
!>                             declared memory type
!> \result memory_type         memory type used for data
! *****************************************************************************
  FUNCTION dbcsr_get_data_memory_type (matrix) &
       RESULT (memory_type)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_memtype_type)                 :: memory_type

    memory_type = matrix%m%data_memory_type
  END FUNCTION dbcsr_get_data_memory_type

! *****************************************************************************
!> \brief Returns the type of memory used for the index in the matrix
!> \param[in] matrix          query this matrix
!> \result memory_type        memory type used for the index
! *****************************************************************************
  FUNCTION dbcsr_get_index_memory_type (matrix) RESULT (memory_type)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_memtype_type)                 :: memory_type

    memory_type = matrix%m%index_memory_type
  END FUNCTION dbcsr_get_index_memory_type


!! *****************************************************************************
!!> \brief Returns whether the matrix uses specially-allocated memory
!!> \param[in] matrix          query this matrix
!!> \param[out] uses_special   whether the matrix uses specially allocated
!!>                            memory
!! *****************************************************************************
!  PURE FUNCTION uses_special_memory_matrix (matrix) RESULT (uses_special)
!    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
!    LOGICAL                                  :: uses_special
!
!    uses_special = matrix%m%data_memory_type .NE. dbcsr_memory_default
!  END FUNCTION uses_special_memory_matrix
!
!
!! *****************************************************************************
!!> \brief Returns whether the data area uses speciall-allocated memory
!!> \param[in] data_area       query this data area
!!> \param[out] uses_special   whether the data area uses specially allocated
!!>                            memory
!! *****************************************************************************
!  PURE FUNCTION uses_special_memory_area (area) RESULT (uses_special)
!    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
!    LOGICAL                                  :: uses_special
!
!    IF (ASSOCIATED (area%d)) THEN
!       uses_special = area%d%memory_type .NE. dbcsr_memory_default
!    ELSE
!       uses_special = .FALSE.
!    ENDIF
!  END FUNCTION uses_special_memory_area
!

! *****************************************************************************
!> \brief Returns the allocated data size
!> \param[in] area       data area
!> \retval data_size      size of data
! *****************************************************************************
  FUNCTION get_data_size_area (area) RESULT (data_size)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    INTEGER                                  :: data_size

    CHARACTER(len=*), PARAMETER :: routineN = 'get_data_size_area', &
      routineP = moduleN//':'//routineN

    data_size = dbcsr_data_get_size (area)
  END FUNCTION get_data_size_area

! *****************************************************************************
!> \brief Returns data type of a data area
!> \param[in] area         data area
!> \result data_type       data type of the data area
! *****************************************************************************
  PURE FUNCTION data_data_get_type (area) RESULT (data_type)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    INTEGER                                  :: data_type

    data_type = dbcsr_data_get_type (area)
  END FUNCTION data_data_get_type


! *****************************************************************************
!> \brief Returns the allocated data size of a DBCSR matrix
!> \param[in] matrix      matrix
!> \retval data_size      size of data
! *****************************************************************************
  FUNCTION get_data_size_matrix (matrix) RESULT (data_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: data_size

    CHARACTER(len=*), PARAMETER :: routineN = 'get_data_size_matrix', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

    data_size = 0
    IF (ASSOCIATED (matrix%m%data_area%d)) THEN
       SELECT CASE (matrix%m%data_area%d%data_type)
          CASE (dbcsr_type_real_8)
             IF (ASSOCIATED (matrix%m%data_area%d%r_dp))&
                  data_size = SIZE (matrix%m%data_area%d%r_dp)
          CASE (dbcsr_type_real_4)
             IF (ASSOCIATED (matrix%m%data_area%d%r_sp))&
                  data_size = SIZE (matrix%m%data_area%d%r_sp)
          CASE (dbcsr_type_complex_8)
             IF (ASSOCIATED (matrix%m%data_area%d%c_dp))&
                  data_size = SIZE (matrix%m%data_area%d%c_dp)
          CASE (dbcsr_type_complex_4)
             IF (ASSOCIATED (matrix%m%data_area%d%c_sp))&
                  data_size = SIZE (matrix%m%data_area%d%c_sp)
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
                  routineN, "Incorrect data type",__LINE__,error)
          END SELECT
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_warning_level, dbcsr_caller_error, routineN,&
            "Uninitialized data area",__LINE__,error)
       data_size = 0
    ENDIF
  END FUNCTION get_data_size_matrix

! *****************************************************************************
!> \brief Get actual data storage used for matrix
!> \param[in] matrix          Count data of this matrix
!> \result data_size          Data size used by matrix
! *****************************************************************************
  PURE FUNCTION get_data_size_ref_matrix (matrix) RESULT (data_size_referenced)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: data_size_referenced

    data_size_referenced = dbcsr_data_get_size_referenced (matrix%m%data_area)
  END FUNCTION get_data_size_ref_matrix

! *****************************************************************************
!> \brief Get actual data storage used for matrix
!> \param[in] area            Count data of this matrix
!> \result data_size          Data size used by matrix
! *****************************************************************************
  PURE FUNCTION get_data_size_ref_area (area) RESULT (data_size_referenced)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    INTEGER                                  :: data_size_referenced

    data_size_referenced = dbcsr_data_get_size_referenced (area)
  END FUNCTION get_data_size_ref_area

! *****************************************************************************
!> \brief Sets the referenced size of the data area
!> \param[in,out] data_area  area for which to set referenced data size
!> \param[in]                set referenced data size to this value
! *****************************************************************************
  PURE SUBROUTINE set_data_size_ref_area (data_area, referenced_size)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: data_area
    INTEGER, INTENT(IN)                      :: referenced_size

    CALL dbcsr_data_set_size_referenced (data_area, referenced_size)
  END SUBROUTINE set_data_size_ref_area


! *****************************************************************************
!> \brief Count actual data storage used for matrix data.
!> \param[in] matrix          Count data of this matrix
!> \param[in,out] error       Error
!> \result data_size          Data size used by matrix
! *****************************************************************************
  FUNCTION dbcsr_get_data_size_used (matrix, error) RESULT (data_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER                                  :: data_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_size_used', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_p, col, &
                                                error_handle, nze, row
    INTEGER, DIMENSION(:), POINTER           :: col_blk_sizes, row_blk_sizes

!type(dbcsr_iterator_type) :: iter
!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    row_blk_sizes => array_data (dbcsr_row_block_sizes (matrix))
    col_blk_sizes => array_data (dbcsr_col_block_sizes (matrix))
    data_size = 0
    !$OMP DO
    DO row = 1, matrix%m%nblkrows_total
       DO blk = matrix%m%row_p(row)+1, matrix%m%row_p(row+1)
          col = matrix%m%col_i(blk)
          blk_p = ABS (matrix%m%blk_p(blk))
          IF (matrix%m%blk_p(blk) .NE. 0) THEN
             nze = row_blk_sizes(row) * col_blk_sizes(col)
             data_size = data_size + nze
          ENDIF
       ENDDO
    ENDDO
    !$OMP END DO
    CALL dbcsr_error_stop (error_handle, error)
  END FUNCTION dbcsr_get_data_size_used



! *****************************************************************************
!> \brief Returns the number of blocks in the matrix
!> \param matrix     matrix from which to get data
!> \param area       data area
! *****************************************************************************
  PURE FUNCTION dbcsr_get_num_blocks (matrix) RESULT (num_blocks)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: num_blocks

    num_blocks = matrix%m%nblks
  END FUNCTION dbcsr_get_num_blocks

! *****************************************************************************
!> \brief Returns the number of non-zero elements in the matrix
!> \param matrix     matrix from which to get data
!> \param area       data area
! *****************************************************************************
  PURE FUNCTION dbcsr_get_nze (matrix) RESULT (num_nze)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: num_nze

    num_nze = matrix%m%nze
  END FUNCTION dbcsr_get_nze

! *****************************************************************************
!> \brief Returns the occupation of the matrix
!> \param matrix     matrix from which to get the occupation
! *****************************************************************************
  FUNCTION dbcsr_get_occupation (matrix) RESULT (occupation)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    REAL(KIND=real_8)                        :: occupation

    INTEGER                                  :: nfullcols, nfullrows
    INTEGER(KIND=int_8)                      :: nze_global
    INTEGER, DIMENSION(:), POINTER           :: row_blk_size

    nze_global = matrix%m%nze
    CALL mp_sum(nze_global,dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))

    nfullrows = dbcsr_nfullrows_total(matrix)
    nfullcols = dbcsr_nfullcols_total(matrix)

    row_blk_size => array_data (matrix%m%row_blk_size)

    IF(nfullrows.NE.0.AND.nfullcols.NE.0) THEN
       IF(dbcsr_has_symmetry (matrix)) THEN
          IF( 2*nze_global .EQ. &
              ( INT(nfullrows,KIND=int_8)*INT(nfullrows+1,KIND=int_8) + SUM( row_blk_size*(row_blk_size-1) ) ) ) THEN
             occupation = 1.0_real_8
          ELSE
             occupation = 2.0_real_8 * REAL(nze_global,real_8) / &
               ( REAL(nfullrows,real_8) * REAL(nfullrows+1,real_8) + &
               SUM( REAL(row_blk_size,real_8) * REAL(row_blk_size-1,real_8) ) )
          ENDIF
       ELSE
          IF (nze_global .EQ. INT(nfullrows,KIND=int_8) * INT(nfullcols,KIND=int_8)) THEN
             occupation = 1.0_real_8
          ELSE
             occupation = REAL(nze_global,real_8) / ( REAL(nfullrows,real_8) * REAL(nfullcols,real_8) )
          ENDIF
       ENDIF
    ELSE
       occupation = 0.0_real_8
    ENDIF
  END FUNCTION dbcsr_get_occupation


! *****************************************************************************
! Arrays
! *****************************************************************************


! *****************************************************************************
!> \brief Initializes a 1-d array of DBCSR matrices
!> \param[out] array    the array of matrices
!> \param[in] n       (optional) number of matrices to be included
! *****************************************************************************
  SUBROUTINE array_init_1d (array, n)
    TYPE(dbcsr_1d_array_obj), INTENT(OUT)    :: array
    INTEGER, INTENT(IN), OPTIONAL            :: n

    NULLIFY (array%s, array%refcount)
    IF (PRESENT (n)) CALL dbcsr_array_new (array, n)
  END SUBROUTINE array_init_1d

! *****************************************************************************
!> \brief Initializes a 2-d array of DBCSR matrices
!> \param[out] array    the array of matrices
!> \param[in] n       (optional) number of matrices to be included
! *****************************************************************************
  SUBROUTINE array_init_2d (array, n)
    TYPE(dbcsr_2d_array_obj), INTENT(OUT)    :: array
    INTEGER, DIMENSION(2), INTENT(IN), &
      OPTIONAL                               :: n

    NULLIFY (array%s, array%refcount)
    IF (PRESENT (n)) CALL dbcsr_array_new (array, n)
  END SUBROUTINE array_init_2d

! *****************************************************************************
!> \brief Allocates a 1-d array of DBCSR matrices
!> \param[in,out] array  the array of matrices
!> \param[in] n        number of matrices to be included
! *****************************************************************************
  SUBROUTINE array_new_1d (array, n)
    TYPE(dbcsr_1d_array_obj), INTENT(INOUT)  :: array
    INTEGER, INTENT(IN)                      :: n

    CHARACTER(len=*), PARAMETER :: routineN = 'array_new_1d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (.NOT. ASSOCIATED (array%s), dbcsr_warning_level,&
         dbcsr_caller_error, routineN,&
         "Array was already allocated, memory leakage will occur",__LINE__,error)
    ALLOCATE (array%s(n))
    ALLOCATE (array%refcount)
    array%refcount = 1
  END SUBROUTINE array_new_1d

! *****************************************************************************
!> \brief Allocates a 2-d array of DBCSR matrices
!> \param[in,out] array  the array of matrices
!> \param[in] n        number of matrices to be included
! *****************************************************************************
  SUBROUTINE array_new_2d (array, n)
    TYPE(dbcsr_2d_array_obj), INTENT(INOUT)  :: array
    INTEGER, DIMENSION(2), INTENT(IN)        :: n

    CHARACTER(len=*), PARAMETER :: routineN = 'array_new_2d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (.NOT. ASSOCIATED (array%s), dbcsr_warning_level,&
         dbcsr_caller_error, routineN,&
         "Array was already allocated, memory leakage will occur",__LINE__,error)
    ALLOCATE (array%s(n(1), n(2)))
    ALLOCATE (array%refcount)
    array%refcount = 1
  END SUBROUTINE array_new_2d

! *****************************************************************************
!> \brief Inserts a DBCSR matrix into the array
!>
!> The matrix is inserted into the array and its reference is registered.
!> \param[in,out] array    the array of matrices
!> \param[in] position   array position into which matrix is inserted
!> \param[in] matrix     matrices to be included
! *****************************************************************************
  SUBROUTINE array_put_1d (array, position, matrix)
    TYPE(dbcsr_1d_array_obj), INTENT(INOUT)  :: array
    INTEGER, INTENT(IN)                      :: position
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'array_put_1d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (array%s), dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Can not insert into unallocated array",__LINE__,error)
    CALL dbcsr_assert(position .GE. LBOUND (array%s, 1)&
         .AND. position .LE. UBOUND (array%s, 1),&
         dbcsr_failure_level, dbcsr_caller_error, routineN,&
         "Insertion position out of bounds.",__LINE__,error)
    array%s(position) = matrix
    CALL dbcsr_hold (array%s(position))
  END SUBROUTINE array_put_1d

! *****************************************************************************
!> \brief Inserts a DBCSR matrix into the array
!>
!> The matrix is inserted into the array and its reference is registered.
!> \param[in,out] array the array of matrices
!> \param[in] position   array position into which matrix is inserted
!> \param[in] matrix  matrices to be included
! *****************************************************************************
  SUBROUTINE array_put_2d (array, position, matrix)
    TYPE(dbcsr_2d_array_obj), INTENT(INOUT)  :: array
    INTEGER, DIMENSION(2), INTENT(IN)        :: position
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'array_put_2d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (array%s), dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Can not insert into unallocated array",__LINE__,error)
    CALL dbcsr_assert(position(1) .GE. LBOUND (array%s, 1)&
         .AND. position(1) .LE. UBOUND (array%s, 1)&
         .AND. position(2) .GE. LBOUND (array%s, 2)&
         .AND. position(2) .LE. UBOUND (array%s, 2),&
         dbcsr_failure_level, dbcsr_caller_error, routineN,&
         "Insertion position out of bounds.",__LINE__,error)
    array%s(position(1), position(2)) = matrix
    CALL dbcsr_hold (array%s(position(1), position(2)))
  END SUBROUTINE array_put_2d

! *****************************************************************************
!> \brief Gets a DBCSR matrix from the array
!> \param[in] array     the array of matrices
!> \param[in] position   array position into which matrix is inserted
!> \param matrix      matrix (no new reference is registered)
! *****************************************************************************
  FUNCTION array_get_1d (array, position) RESULT (matrix)
    TYPE(dbcsr_1d_array_obj), INTENT(IN)     :: array
    INTEGER, INTENT(IN)                      :: position
    TYPE(dbcsr_obj)                          :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'array_get_1d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (array%s), dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Can not read from an unallocated array",__LINE__,error)
    CALL dbcsr_assert (position .GE. LBOUND (array%s,1)&
         .AND. position .LE. UBOUND (array%s,1),&
         dbcsr_failure_level, dbcsr_caller_error, routineN,&
         "Read position out of bounds.",__LINE__,error)
    matrix = array%s(position)
  END FUNCTION array_get_1d

! *****************************************************************************
!> \brief Gets a DBCSR matrix from the array
!> \param[in] array     the array of matrices
!> \param[in] position   array position into which matrix is inserted
!> \param matrix      matrix (no new reference is registered)
! *****************************************************************************
  FUNCTION array_get_2d (array, position) RESULT (matrix)
    TYPE(dbcsr_2d_array_obj), INTENT(IN)     :: array
    INTEGER, DIMENSION(2), INTENT(IN)        :: position
    TYPE(dbcsr_obj)                          :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'array_get_2d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (array%s), dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Can not read from an unallocated array",__LINE__,error)
    CALL dbcsr_assert(position(1) .GE. LBOUND (array%s, 1)&
         .AND. position(1) .LE. UBOUND (array%s, 1)&
         .AND. position(2) .GE. LBOUND (array%s, 2)&
         .AND. position(2) .LE. UBOUND (array%s, 2),&
         dbcsr_failure_level, dbcsr_caller_error, routineN,&
         "Read position out of bounds.",__LINE__,error)
    matrix = array%s(position(1), position(2))
  END FUNCTION array_get_2d

! *****************************************************************************
!> \brief Destroys a 1-d array and releases all held resources
!> \param[in,out] array    the array of matrices
! *****************************************************************************
  SUBROUTINE array_destroy_1d (array)
    TYPE(dbcsr_1d_array_obj), INTENT(INOUT)  :: array

    CHARACTER(len=*), PARAMETER :: routineN = 'array_destroy_1d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (array%s), dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Can not destroy an unallocated array",__LINE__,error)
    CALL dbcsr_assert (ASSOCIATED (array%refcount),&
         dbcsr_failure_level, dbcsr_caller_error, routineN,&
         "Reference count does not exist.",__LINE__,error)
    CALL dbcsr_assert (array%refcount .EQ. 0, dbcsr_warning_level, dbcsr_caller_error,&
         routineN, "References are still held",__LINE__,error)
    DO i = LBOUND (array%s, 1), UBOUND (array%s, 1)
       CALL dbcsr_release (array%s(i))
    ENDDO
    DEALLOCATE (array%s)
    DEALLOCATE (array%refcount)
    CALL dbcsr_allocate_matrix_array (array)
  END SUBROUTINE array_destroy_1d

! *****************************************************************************
!> \brief Destroys a 2-d array and releases all held resources
!> \param[in,out] array    the array of matrices
! *****************************************************************************
  SUBROUTINE array_destroy_2d (array)
    TYPE(dbcsr_2d_array_obj), INTENT(INOUT)  :: array

    CHARACTER(len=*), PARAMETER :: routineN = 'array_destroy_2d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (array%s), dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Can not destroy an unallocated array",__LINE__,error)
    CALL dbcsr_assert (ASSOCIATED (array%refcount),&
         dbcsr_failure_level, dbcsr_caller_error, routineN,&
         "Reference count does not exist.",__LINE__,error)
    CALL dbcsr_assert (array%refcount .EQ. 0, dbcsr_warning_level, dbcsr_caller_error,&
         routineN, "References are still held",__LINE__,error)
    DO i = LBOUND (array%s, 1), UBOUND (array%s, 1)
       DO j = LBOUND (array%s, 2), UBOUND (array%s, 2)
          CALL dbcsr_release (array%s(i, j))
       ENDDO
    ENDDO
    DEALLOCATE (array%s)
    DEALLOCATE (array%refcount)
    CALL dbcsr_allocate_matrix_array (array)
  END SUBROUTINE array_destroy_2d


! *****************************************************************************
!> \brief Registers another reference for a 1-d array
!> \param[in,out] array    the array of matrices
! *****************************************************************************
  SUBROUTINE array_hold_1d (array)
    TYPE(dbcsr_1d_array_obj), INTENT(INOUT)  :: array

    CHARACTER(len=*), PARAMETER :: routineN = 'array_hold_1d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (array%refcount),&
         dbcsr_failure_level, dbcsr_caller_error, routineN,&
         "Reference count does not exist.",__LINE__,error)
    array%refcount = array%refcount + 1
  END SUBROUTINE array_hold_1d

! *****************************************************************************
!> \brief Registers another reference for a 2-d array
!> \param[in,out] array    the array of matrices
! *****************************************************************************
  SUBROUTINE array_hold_2d (array)
    TYPE(dbcsr_2d_array_obj), INTENT(INOUT)  :: array

    CHARACTER(len=*), PARAMETER :: routineN = 'array_hold_2d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (array%refcount),&
         dbcsr_failure_level, dbcsr_caller_error, routineN,&
         "Reference count does not exist.",__LINE__,error)
    array%refcount = array%refcount + 1
  END SUBROUTINE array_hold_2d


! *****************************************************************************
!> \brief Releases a reference for a 1-d array
!>
!> If there are no references left, the array is destroyed.
!> \param[in,out] array    the array of matrices
! *****************************************************************************
  SUBROUTINE array_release_1d (array)
    TYPE(dbcsr_1d_array_obj), INTENT(INOUT)  :: array

    CHARACTER(len=*), PARAMETER :: routineN = 'array_release_1d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (array%refcount),&
         dbcsr_failure_level, dbcsr_caller_error, routineN,&
         "Reference count does not exist.",__LINE__,error)
    array%refcount = array%refcount - 1
    IF (array%refcount .EQ. 0) THEN
       CALL dbcsr_array_destroy (array)
    ENDIF
  END SUBROUTINE array_release_1d

! *****************************************************************************
!> \brief Releases a reference for a 2-d array
!>
!> If there are no references left, the array is destroyed.
!> \param[in,out] array    the array of matrices
! *****************************************************************************
  SUBROUTINE array_release_2d (array)
    TYPE(dbcsr_2d_array_obj), INTENT(INOUT)  :: array

    CHARACTER(len=*), PARAMETER :: routineN = 'array_release_2d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (array%refcount),&
         dbcsr_failure_level, dbcsr_caller_error, routineN,&
         "Reference count does not exist.",__LINE__,error)
    array%refcount = array%refcount - 1
    IF (array%refcount .EQ. 0) THEN
       CALL dbcsr_array_destroy (array)
    ENDIF
  END SUBROUTINE array_release_2d

! *****************************************************************************
!> \brief Releases all matrices in a 1-d arrray.
!> \param[in] source          input matrix
!> \param[in,out] marray      matrix array
! *****************************************************************************
  SUBROUTINE dbcsr_destroy_1d_array(marray, error)
    TYPE(dbcsr_array_type), INTENT(INOUT)    :: marray
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_destroy_1d_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    DO i = LBOUND(marray%mats,1), UBOUND (marray%mats,1)
       CALL dbcsr_destroy (marray%mats(i), error=error, force=.TRUE.)
    ENDDO
    DEALLOCATE (marray%mats)
    !CALL dbcsr_destroy_image_dist(marray%image_dist)
  END SUBROUTINE dbcsr_destroy_1d_array


! *****************************************************************************
!> \brief Releases all matrices in 2-d arrray.
!> \param[in] source          input matrix
!> \param[in,out] marray      matrix array
! *****************************************************************************
  SUBROUTINE dbcsr_destroy_2d_array(marray,error)
    TYPE(dbcsr_2d_array_type), INTENT(INOUT) :: marray
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_destroy_2d_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, row

!   ---------------------------------------------------------------------------

    DO row = LBOUND(marray%mats,1), UBOUND (marray%mats,1)
       DO col = LBOUND(marray%mats,2), UBOUND (marray%mats,2)
          CALL dbcsr_destroy (marray%mats(row, col), error=error, force=.TRUE.)
       ENDDO
    ENDDO
    CALL dbcsr_image_dist_release (marray%image_dist, error)
    DEALLOCATE (marray%mats)
    NULLIFY (marray%mats)
  END SUBROUTINE dbcsr_destroy_2d_array

! *****************************************************************************
!> \brief Releases a reference to and possible deallocates an image
!         distribution
!*****************************************************************************
  SUBROUTINE dbcsr_image_dist_release (imgdist, error)
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(INOUT)                          :: imgdist
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF (ASSOCIATED (imgdist%i)) THEN
       imgdist%i%refcount = imgdist%i%refcount - 1
       IF (imgdist%i%refcount .EQ. 0) THEN
          CALL dbcsr_destroy_image_dist (imgdist%i)
          DEALLOCATE (imgdist%i)
          NULLIFY (imgdist%i)
       ENDIF
    ENDIF
  END SUBROUTINE dbcsr_image_dist_release

! *****************************************************************************
!> \brief Retains a reference to an image distribution
! *****************************************************************************
  SUBROUTINE dbcsr_image_dist_hold (imgdist, error)
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(INOUT)                          :: imgdist
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    imgdist%i%refcount = imgdist%i%refcount + 1
  END SUBROUTINE dbcsr_image_dist_hold

! *****************************************************************************
!> \brief Initialized an image distribution
!> \par Akin to nullify.
! *****************************************************************************
  SUBROUTINE dbcsr_image_dist_init (imgdist, error)
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(OUT)                            :: imgdist
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    NULLIFY (imgdist%i)
  END SUBROUTINE dbcsr_image_dist_init

! *****************************************************************************
!> \brief Destroys a DBCSR distribution for a matrix multiplication based on
!>        the right matrix
!> \param[inout] imgdist_product        product distribution repetition
! *****************************************************************************
  SUBROUTINE dbcsr_destroy_image_dist(imgdist)
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(INOUT)                          :: imgdist

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    CALL array_release (imgdist%row_image)
    CALL array_release (imgdist%col_image)
    CALL dbcsr_distribution_release (imgdist%main)
    !
    CALL array_release (imgdist%vrow_dist)
    CALL array_release (imgdist%vcol_dist)
    !
    IF (imgdist%has_other_vl_rows) THEN
       DO i = LBOUND (imgdist%other_vl_rows,1), UBOUND (imgdist%other_vl_rows,1)
          CALL array_release (imgdist%other_vl_rows(i))
       END DO
       DEALLOCATE (imgdist%other_vl_rows)
       NULLIFY (imgdist%other_vl_rows)
       imgdist%has_other_vl_rows = .FALSE.
    ENDIF
    !
    IF (imgdist%has_other_vl_cols) THEN
       DO i = LBOUND (imgdist%other_vl_cols,1), UBOUND (imgdist%other_vl_cols,1)
          CALL array_release (imgdist%other_vl_cols(i))
       END DO
       DEALLOCATE (imgdist%other_vl_cols)
       NULLIFY (imgdist%other_vl_cols)
       imgdist%has_other_vl_cols = .FALSE.
    ENDIF
    !
    IF (imgdist%has_global_vrow_map) THEN
       CALL array_release (imgdist%global_vrow_map)
    ENDIF
    IF (imgdist%has_global_vcol_map) THEN
       CALL array_release (imgdist%global_vcol_map)
    ENDIF
  END SUBROUTINE dbcsr_destroy_image_dist

! *****************************************************************************
! Mutable data
! *****************************************************************************

! *****************************************************************************
!> \brief Initializes a mutable data type
!> \param[out] mutable        mutable data
! *****************************************************************************
  SUBROUTINE dbcsr_mutable_init (mutable)
    TYPE(dbcsr_mutable_obj), INTENT(OUT)     :: mutable

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mutable_init', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    NULLIFY (mutable%m)
  END SUBROUTINE dbcsr_mutable_init

! *****************************************************************************
!> \brief Destroys a mutable data type
!> \param[in,out] mutable     mutable data
! *****************************************************************************
  SUBROUTINE dbcsr_mutable_destroy (mutable)
    TYPE(dbcsr_mutable_obj), INTENT(INOUT)   :: mutable

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mutable_destroy', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (mutable%m)) THEN
       !CALL dbcsr_assert (mutable%m%refcount .EQ. 0, dbcsr_warning_level,&
       !     dbcsr_caller_error, routineN, "Destroying with non-0 reference count")
       CALL btree_destroy_s (mutable%m%btree_s)
       CALL btree_destroy_d (mutable%m%btree_d)
       CALL btree_destroy_c (mutable%m%btree_c)
       CALL btree_destroy_z (mutable%m%btree_z)
       DEALLOCATE (mutable%m)
    ENDIF
    NULLIFY (mutable%m)
  END SUBROUTINE dbcsr_mutable_destroy


! *****************************************************************************
!> \brief Registers another reference to the mutable data type
!> \param[in,out] mutable     mutable data
! *****************************************************************************
  SUBROUTINE dbcsr_mutable_hold (mutable)
    TYPE(dbcsr_mutable_obj), INTENT(INOUT)   :: mutable

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mutable_hold', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (mutable%m), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Mutable data area not instantiated",__LINE__,error)
    mutable%m%refcount = mutable%m%refcount + 1
  END SUBROUTINE dbcsr_mutable_hold

! *****************************************************************************
!> \brief Deregisters a reference to the mutable data type
!>
!> The object is destroy when there is no reference to it left.
!> \param[in,out] mutable     mutable data
! *****************************************************************************
  SUBROUTINE dbcsr_mutable_release (mutable)
    TYPE(dbcsr_mutable_obj), INTENT(INOUT)   :: mutable

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mutable_release', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (mutable%m), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Mutable data area not instantiated",__LINE__,error)
    mutable%m%refcount = mutable%m%refcount - 1
    IF (mutable%m%refcount .EQ. 0) THEN
       CALL dbcsr_mutable_destroy (mutable)
    ENDIF
  END SUBROUTINE dbcsr_mutable_release

! *****************************************************************************
!> \brief Creates a new mutable instance.
!>
!> \param[in,out] mutable     mutable data
!> \param[in] data_type       data type to be stored here (see dbcsr_types for
!>                            possibilities)
! *****************************************************************************
  SUBROUTINE dbcsr_mutable_new (mutable, data_type)
    TYPE(dbcsr_mutable_obj), INTENT(INOUT)   :: mutable
    INTEGER, INTENT(IN)                      :: data_type

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mutable_new', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (.NOT.ASSOCIATED (mutable%m), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Mutable data area already instantiated",__LINE__,error)
    CALL dbcsr_assert (data_type.EQ.dbcsr_type_real_4&
         .OR. data_type.EQ.dbcsr_type_real_8&
         .OR. data_type.EQ.dbcsr_type_complex_4&
         .OR. data_type.EQ.dbcsr_type_complex_8, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Invalid data type",__LINE__,error)
    ALLOCATE (mutable%m)
    mutable%m%refcount = 1
    mutable%m%data_type = data_type
    CALL btree_new_s (mutable%m%btree_s)
    CALL btree_new_d (mutable%m%btree_d)
    CALL btree_new_c (mutable%m%btree_c)
    CALL btree_new_z (mutable%m%btree_z)
  END SUBROUTINE dbcsr_mutable_new

! *****************************************************************************
!> \brief Deregisters a reference to the mutable data type
!>
!> The object is destroy when there is no reference to it left.
!> \param[in] mutable         mutable data
!> \param[out] instantiated   whether the object is instantiated
! *****************************************************************************
  PURE FUNCTION dbcsr_mutable_instantiated (mutable) RESULT (instantiated)
    TYPE(dbcsr_mutable_obj), INTENT(IN)      :: mutable
    LOGICAL                                  :: instantiated

!   ---------------------------------------------------------------------------

    instantiated = ASSOCIATED (mutable%m)
  END FUNCTION dbcsr_mutable_instantiated


!! *****************************************************************************
!!> \brief Returns the entire data for a matrix.
!!> \par Warning
!!>      This routine should only be used within DBCSR code.
!!> \param[in] area       data area
!!> \param[in] coersion   force datatype
!!> \param[out] data      pointer to data
!! *****************************************************************************
!  SUBROUTINE get_data_m_[nametype1] (matrix, DATA)
!    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
!    [type1], DIMENSION(:), POINTER :: DATA
!
!    CALL get_data_[nametype1] (matrix%m%data_area, DATA)
!  END SUBROUTINE get_data_m_[nametype1]


END MODULE dbcsr_methods

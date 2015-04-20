!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
                                             array_nullify,&
                                             array_release
  USE btree_I8_k_cp2d_v,               ONLY: btree_destroy_c => btree_delete,&
                                             btree_new_c => btree_new
  USE btree_I8_k_dp2d_v,               ONLY: btree_destroy_d => btree_delete,&
                                             btree_new_d => btree_new
  USE btree_I8_k_sp2d_v,               ONLY: btree_destroy_s => btree_delete,&
                                             btree_new_s => btree_new
  USE btree_I8_k_zp2d_v,               ONLY: btree_destroy_z => btree_delete,&
                                             btree_new_z => btree_new
  USE dbcsr_data_methods,              ONLY: dbcsr_data_init,&
                                             dbcsr_data_release
  USE dbcsr_error_handling,            ONLY: &
       dbcsr_assert, dbcsr_caller_error, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_failure_level, dbcsr_fatal_level, &
       dbcsr_internal_error, dbcsr_warning_level, dbcsr_wrong_args_error
  USE dbcsr_ptr_util,                  ONLY: memory_deallocate
  USE dbcsr_types,                     ONLY: &
       dbcsr_2d_array_type, dbcsr_array_type, dbcsr_data_obj, &
       dbcsr_distribution_obj, dbcsr_imagedistribution_obj, &
       dbcsr_imagedistribution_type, dbcsr_memtype_type, dbcsr_mp_obj, &
       dbcsr_mutable_obj, dbcsr_obj, dbcsr_type, dbcsr_type_antihermitian, &
       dbcsr_type_antisymmetric, dbcsr_type_complex_4, dbcsr_type_complex_8, &
       dbcsr_type_hermitian, dbcsr_type_invalid, dbcsr_type_no_symmetry, &
       dbcsr_type_real_4, dbcsr_type_real_8, dbcsr_type_symmetric, &
       dbcsr_work_type
  USE kinds,                           ONLY: default_string_length
  USE message_passing,                 ONLY: mp_comm_free

  !$ USE OMP_LIB
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_methods'

  INTEGER, PUBLIC, SAVE :: dbcsr_matrix_counter = 111111

  PUBLIC :: dbcsr_init
  PUBLIC :: dbcsr_release
  PUBLIC :: dbcsr_valid_index, dbcsr_is_initialized
  PUBLIC :: dbcsr_release_locals
  PUBLIC :: dbcsr_distribution,&
            dbcsr_get_matrix_type, dbcsr_get_data_type, dbcsr_get_replication_type,&
            dbcsr_row_block_sizes, dbcsr_col_block_sizes,&
            dbcsr_nblkrows_total, dbcsr_nblkcols_total, dbcsr_nfullrows_total,&
            dbcsr_nfullcols_total, dbcsr_nblkcols_local, dbcsr_nblkrows_local,&
            dbcsr_max_row_size, dbcsr_max_col_size,&
            dbcsr_get_index_memory_type, dbcsr_get_data_memory_type,&
            dbcsr_name, dbcsr_setname, dbcsr_get_data_size,&
            dbcsr_use_mutable, dbcsr_wm_use_mutable, dbcsr_has_symmetry,&
            dbcsr_get_nze, dbcsr_nfullrows_local
  PUBLIC :: dbcsr_get_data_size_used
  PUBLIC :: dbcsr_col_block_offsets, dbcsr_row_block_offsets
  PUBLIC :: dbcsr_data_area
  PUBLIC :: dbcsr_get_num_blocks


  PUBLIC :: dbcsr_blk_row_size, dbcsr_blk_column_size,&
            dbcsr_blk_row_offset, dbcsr_blk_col_offset

  PUBLIC :: dbcsr_destroy_array
  PUBLIC :: dbcsr_image_dist_init, dbcsr_image_dist_hold, dbcsr_image_dist_release

  PUBLIC :: dbcsr_mutable_init, dbcsr_mutable_new, dbcsr_mutable_destroy,&
            dbcsr_mutable_release,&
            dbcsr_mutable_instantiated

  PUBLIC :: dbcsr_distribution_release, dbcsr_distribution_init
  PUBLIC :: dbcsr_mp_release, dbcsr_mp_grid_remove

  INTERFACE dbcsr_init
     MODULE PROCEDURE dbcsr_init_type, dbcsr_init_obj
  END INTERFACE

  INTERFACE dbcsr_valid_index
     MODULE PROCEDURE dbcsr_valid_index_type, dbcsr_valid_index_obj
  END INTERFACE

  INTERFACE dbcsr_is_initialized
     MODULE PROCEDURE dbcsr_is_initialized_type, dbcsr_is_initialized_obj
  END INTERFACE

  INTERFACE dbcsr_blk_row_size
     MODULE PROCEDURE dbcsr_blk_row_size_type, dbcsr_blk_row_size_obj
  END INTERFACE

  INTERFACE dbcsr_blk_column_size
     MODULE PROCEDURE dbcsr_blk_column_size_type, dbcsr_blk_column_size_obj
  END INTERFACE

  ! For the 1-D and 2-D arrays

  INTERFACE dbcsr_destroy_array
     MODULE PROCEDURE dbcsr_destroy_2d_array, dbcsr_destroy_1d_array
  END INTERFACE



CONTAINS


! *****************************************************************************
!> \brief Initializes a DBCSR matrix but does not create it.
!> \param[out] matrix         uninitialized matrix
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
!> \brief ...
!> \param dist ...
! *****************************************************************************
  SUBROUTINE dbcsr_distribution_init (dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist

    NULLIFY (dist%d)
  END SUBROUTINE dbcsr_distribution_init

! *****************************************************************************
!> \brief Initializes a DBCSR matrix but does not allocate any memory.
!> \param[out] matrix         uninitialized matrix
! *****************************************************************************
  SUBROUTINE dbcsr_init_obj (matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_obj', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_init(matrix%m)
  END SUBROUTINE dbcsr_init_obj


! *****************************************************************************
!> \brief Returns whether the index structure of the matrix is valid.
!> \param[in] matrix          verify index validity of this matrix
!> \retval valid_index        index validity
! *****************************************************************************
  PURE FUNCTION dbcsr_valid_index_type (matrix) RESULT (valid_index)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    LOGICAL                                  :: valid_index

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

    initialized = dbcsr_is_initialized_type (matrix%m)
  END FUNCTION dbcsr_is_initialized_obj


! *****************************************************************************
!> \brief Releases a reference for a DBCSR matrix
!>
!> If there are no references left, the matrix is destroyed.
!> \param[in,out] matrix DBCSR matrix
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
          CALL dbcsr_init (matrix%m)
       ENDIF
    ENDIF
    CALL dbcsr_error_stop (error_handle, error=error)
  END SUBROUTINE dbcsr_destroy

! *****************************************************************************
!> \brief Releases and potentially destrops a distribution
!> \param dist ...
! *****************************************************************************
  SUBROUTINE dbcsr_distribution_release(dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (dist%d)) THEN
       dist%d%refcount = dist%d%refcount - 1
       IF (dist%d%refcount .EQ. 0) THEN
          CALL array_release (dist%d%row_dist_block)
          CALL array_release (dist%d%col_dist_block)
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

! *****************************************************************************
!> \brief ...
!> \param dist ...
!> \param error ...
! *****************************************************************************
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
!> \brief ...
!> \param matrix ...
!> \param error ...
! *****************************************************************************
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

! Pertaining to the dbcsr matrix.

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nblkrows_total ...
! *****************************************************************************
  FUNCTION dbcsr_nblkrows_total(matrix) RESULT (nblkrows_total)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nblkrows_total

    nblkrows_total = matrix%m%nblkrows_total
  END FUNCTION dbcsr_nblkrows_total

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nblkcols_total ...
! *****************************************************************************
  FUNCTION dbcsr_nblkcols_total(matrix) RESULT (nblkcols_total)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nblkcols_total

    nblkcols_total = matrix%m%nblkcols_total
  END FUNCTION dbcsr_nblkcols_total
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nfullrows_total ...
! *****************************************************************************
  FUNCTION dbcsr_nfullrows_total(matrix) RESULT (nfullrows_total)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nfullrows_total

    nfullrows_total = matrix%m%nfullrows_total
  END FUNCTION dbcsr_nfullrows_total
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nfullcols_total ...
! *****************************************************************************
  FUNCTION dbcsr_nfullcols_total(matrix) RESULT (nfullcols_total)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nfullcols_total

    nfullcols_total = matrix%m%nfullcols_total
  END FUNCTION dbcsr_nfullcols_total
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nblkrows_local ...
! *****************************************************************************
  FUNCTION dbcsr_nblkrows_local(matrix) RESULT (nblkrows_local)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nblkrows_local

    nblkrows_local = matrix%m%nblkrows_local
  END FUNCTION dbcsr_nblkrows_local
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nblkcols_local ...
! *****************************************************************************
  FUNCTION dbcsr_nblkcols_local(matrix) RESULT (nblkcols_local)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nblkcols_local

    nblkcols_local = matrix%m%nblkcols_local
  END FUNCTION dbcsr_nblkcols_local
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nfullrows_local ...
! *****************************************************************************
  FUNCTION dbcsr_nfullrows_local(matrix) RESULT (nfullrows_local)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nfullrows_local

    nfullrows_local = matrix%m%nfullrows_local
  END FUNCTION dbcsr_nfullrows_local
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval nfullcols_local ...
! *****************************************************************************
  FUNCTION dbcsr_nfullcols_local(matrix) RESULT (nfullcols_local)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nfullcols_local

    nfullcols_local = matrix%m%nfullcols_local
  END FUNCTION dbcsr_nfullcols_local
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval max_row_size ...
! *****************************************************************************
  FUNCTION dbcsr_max_row_size(matrix) RESULT (max_row_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: max_row_size

    max_row_size = matrix%m%max_rbs
  END FUNCTION dbcsr_max_row_size
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval max_col_size ...
! *****************************************************************************
  FUNCTION dbcsr_max_col_size(matrix) RESULT (max_col_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: max_col_size

    max_col_size = matrix%m%max_cbs
  END FUNCTION dbcsr_max_col_size

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval distribution ...
! *****************************************************************************
  FUNCTION dbcsr_distribution (matrix) RESULT (distribution)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_distribution_obj)             :: distribution

    distribution = matrix%m%dist
  END FUNCTION dbcsr_distribution

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval name ...
! *****************************************************************************
  FUNCTION dbcsr_name (matrix) RESULT (name)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    CHARACTER(len=default_string_length)     :: name

    name = matrix%m%name
  END FUNCTION dbcsr_name

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param newname ...
!> \retval name ...
! *****************************************************************************

  SUBROUTINE dbcsr_setname (matrix,newname)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    CHARACTER(len=*), INTENT(IN)             :: newname

    matrix%m%name = newname
  END SUBROUTINE dbcsr_setname

! *****************************************************************************
!> \brief Returns whether this work matrix uses the mutable type
!> \param[in] wm              work matrix
!> \retval use_mutable        use the mutable and not append-only working
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
!> \retval use_mutable        use the mutable and not append-only working
!>                            structures
! *****************************************************************************
  PURE FUNCTION dbcsr_use_mutable (matrix) RESULT (use_mutable)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    LOGICAL                                  :: use_mutable

!   ---------------------------------------------------------------------------

    use_mutable = matrix%work_mutable
  END FUNCTION dbcsr_use_mutable


! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval row_blk_sizes ...
! *****************************************************************************
  FUNCTION dbcsr_row_block_sizes (matrix) RESULT (row_blk_sizes)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, DIMENSION(:), POINTER           :: row_blk_sizes

    row_blk_sizes => array_data(matrix%m%row_blk_size)
  END FUNCTION dbcsr_row_block_sizes

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval col_blk_sizes ...
! *****************************************************************************
  FUNCTION dbcsr_col_block_sizes (matrix) RESULT (col_blk_sizes)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, DIMENSION(:), POINTER           :: col_blk_sizes

    col_blk_sizes => array_data(matrix%m%col_blk_size)
  END FUNCTION dbcsr_col_block_sizes

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval col_blk_offsets ...
! *****************************************************************************
  FUNCTION dbcsr_col_block_offsets (matrix) RESULT (col_blk_offsets)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, DIMENSION(:), POINTER           :: col_blk_offsets

    col_blk_offsets => array_data(matrix%m%col_blk_offset)
  END FUNCTION dbcsr_col_block_offsets

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \retval row_blk_offsets ...
! *****************************************************************************
  FUNCTION dbcsr_row_block_offsets (matrix) RESULT (row_blk_offsets)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, DIMENSION(:), POINTER           :: row_blk_offsets

    row_blk_offsets => array_data(matrix%m%row_blk_offset)
  END FUNCTION dbcsr_row_block_offsets

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
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param row ...
!> \retval row_size ...
! *****************************************************************************
  PURE FUNCTION dbcsr_blk_row_size_obj (matrix, row) RESULT (row_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: row
    INTEGER                                  :: row_size

    row_size = matrix%m%row_blk_size%low%data(row)
  END FUNCTION dbcsr_blk_row_size_obj
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param row ...
!> \retval row_offset ...
! *****************************************************************************
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
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param column ...
!> \retval column_size ...
! *****************************************************************************
  PURE FUNCTION dbcsr_blk_column_size_obj (matrix, column) RESULT (column_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: column
    INTEGER                                  :: column_size

    column_size = matrix%m%col_blk_size%low%data(column)
  END FUNCTION dbcsr_blk_column_size_obj
! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param col ...
!> \retval col_offset ...
! *****************************************************************************
  PURE FUNCTION dbcsr_blk_col_offset (matrix, col) RESULT (col_offset)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: col
    INTEGER                                  :: col_offset

    col_offset = matrix%m%col_blk_offset%low%data(col)
  END FUNCTION dbcsr_blk_col_offset

! *****************************************************************************
!> \brief Returns the data area
!> \param matrix     matrix from which to get data
!> \retval data_area data area
! *****************************************************************************
  FUNCTION dbcsr_data_area (matrix) RESULT (data_area)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    TYPE(dbcsr_data_obj)                     :: data_area

    data_area = matrix%data_area
  END FUNCTION dbcsr_data_area


! *****************************************************************************
!> \brief Returns the matrix type
!> \param matrix              query this matrix
!> \retval matrix_type        matrix_type (see dbcsr_types.F for possible
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
!> \retval has_symmetry       matrix has symmetry
! *****************************************************************************
  PURE FUNCTION dbcsr_has_symmetry (matrix) RESULT (has_symmetry)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL                                  :: has_symmetry

    has_symmetry = matrix%m%symmetry
  END FUNCTION dbcsr_has_symmetry

! *****************************************************************************
!> \brief Returns the data type stored in the matrix
!> \param matrix              query this matrix
!> \retval repl_type repl_type (see dbcsr_types.F for possible
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
!> \retval data_type          data_type (see dbcsr_types.F for possible
!>                            values)
! *****************************************************************************
  PURE FUNCTION dbcsr_get_data_type (matrix) RESULT (data_type)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: data_type

    data_type = matrix%m%data_type
  END FUNCTION dbcsr_get_data_type

! *****************************************************************************
!> \brief Returns the type of memory used for data in the matrix
!> \param[in] matrix           query this matrix
!> \retval memory_type         memory type used for data
!> \note It returns the declared data type, not the actually used type
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
!> \retval memory_type        memory type used for the index 
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
!> \brief Returns the allocated data size of a DBCSR matrix
!> \param[in] matrix      matrix
!> \retval data_size      size of data
! *****************************************************************************
  FUNCTION dbcsr_get_data_size (matrix) RESULT (data_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: data_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_size', &
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
  END FUNCTION dbcsr_get_data_size

! *****************************************************************************
!> \brief Count actual data storage used for matrix data.
!> \param[in] matrix          Count data of this matrix
!> \param[in,out] error       Error
!> \retval data_size          Data size used by matrix 
! *****************************************************************************
  FUNCTION dbcsr_get_data_size_used (matrix, error) RESULT (data_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER                                  :: data_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_size_used', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, error_handle, nze, &
                                                row
    INTEGER, DIMENSION(:), POINTER           :: col_blk_sizes, row_blk_sizes

!type(dbcsr_iterator_type) :: iter
!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    row_blk_sizes => dbcsr_row_block_sizes (matrix)
    col_blk_sizes => dbcsr_col_block_sizes (matrix)
    data_size = 0
    !$OMP DO
    DO row = 1, matrix%m%nblkrows_total
       DO blk = matrix%m%row_p(row)+1, matrix%m%row_p(row+1)
          col = matrix%m%col_i(blk)
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
!> \retval num_blocks ...
! *****************************************************************************
  PURE FUNCTION dbcsr_get_num_blocks (matrix) RESULT (num_blocks)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: num_blocks

    num_blocks = matrix%m%nblks
  END FUNCTION dbcsr_get_num_blocks


! *****************************************************************************
!> \brief Returns the number of non-zero elements in the matrix
!> \param matrix     matrix from which to get data
!> \retval num_nze ...
! *****************************************************************************
  PURE FUNCTION dbcsr_get_nze (matrix) RESULT (num_nze)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: num_nze

    num_nze = matrix%m%nze
  END FUNCTION dbcsr_get_nze


! *****************************************************************************
! Arrays
! *****************************************************************************

! *****************************************************************************
!> \brief Releases all matrices in a 1-d arrray.
!> \param[in,out] marray      matrix array
!> \param error ...
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
!> \param[in,out] marray      matrix array
!> \param error ...
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
!>         distribution
!> \param imgdist ...
!> \param error ...
! *****************************************************************************
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
!> \param imgdist ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_image_dist_hold (imgdist, error)
    TYPE(dbcsr_imagedistribution_obj), &
      INTENT(INOUT)                          :: imgdist
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    imgdist%i%refcount = imgdist%i%refcount + 1
  END SUBROUTINE dbcsr_image_dist_hold

! *****************************************************************************
!> \brief Initialized an image distribution
!> \param imgdist ...
!> \param error ...
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
!> \param imgdist product distribution repetition
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
!> \retval instantiated       whether the object is instantiated
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

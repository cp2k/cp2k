!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
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
  USE dbcsr_error_handling,            ONLY: &
       dbcsr_assert, dbcsr_caller_error, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_failure_level, dbcsr_fatal_level, &
       dbcsr_internal_error, dbcsr_unimplemented_error_nr, &
       dbcsr_warning_level, dbcsr_wrong_args_error
  USE dbcsr_kinds,                     ONLY: default_string_length,&
                                             real_4,&
                                             real_8
  USE dbcsr_message_passing,           ONLY: mp_allocate,&
                                             mp_cart_create,&
                                             mp_cart_sub,&
                                             mp_comm_free,&
                                             mp_deallocate,&
                                             mp_sum
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size
  USE dbcsr_types,                     ONLY: &
       dbcsr_1d_array_obj, dbcsr_2d_array_obj, dbcsr_2d_array_type, &
       dbcsr_array_type, dbcsr_data_obj, dbcsr_distribution_obj, &
       dbcsr_imagedistribution_type, dbcsr_magic_number, dbcsr_mp_obj, &
       dbcsr_mutable_obj, dbcsr_obj, dbcsr_scalar_type, dbcsr_type, &
       dbcsr_type_antisymmetric, dbcsr_type_complex_4, &
       dbcsr_type_complex_4_2d, dbcsr_type_complex_8, &
       dbcsr_type_complex_8_2d, dbcsr_type_conjugate, dbcsr_type_normal, &
       dbcsr_type_real_4, dbcsr_type_real_4_2d, dbcsr_type_real_8, &
       dbcsr_type_real_8_2d, dbcsr_type_symmetric, dbcsr_type_transposed, &
       dbcsr_work_type

  !$ USE OMP_LIB
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_methods'

  ! Workaround for Intel MPI.
  LOGICAL, PARAMETER, PRIVATE :: use_subcommunicators = .FALSE.


  PUBLIC :: dbcsr_init, dbcsr_init_p
  PUBLIC :: dbcsr_hold, dbcsr_release, dbcsr_release_p, dbcsr_destroy
  PUBLIC :: dbcsr_valid_index, dbcsr_is_initialized
  PUBLIC :: dbcsr_scalar, dbcsr_scalar_one, dbcsr_scalar_zero,&
            dbcsr_scalar_are_equal, dbcsr_scalar_negative
  PUBLIC :: dbcsr_mp_new, dbcsr_mp_hold, dbcsr_mp_release
  PUBLIC :: dbcsr_mp_pgrid, dbcsr_mp_numnodes, dbcsr_mp_mynode, dbcsr_mp_group,&
            dbcsr_mp_new_transposed, dbcsr_mp_nprows, dbcsr_mp_npcols,&
            dbcsr_mp_myprow, dbcsr_mp_mypcol, dbcsr_mp_pgrid_equal,&
            dbcsr_mp_get_coordinates,&
            dbcsr_mp_my_row_group, dbcsr_mp_my_col_group,&
            dbcsr_mp_has_subgroups, dbcsr_mp_get_process,&
            dbcsr_mp_grid_setup, dbcsr_mp_grid_remove
  PUBLIC :: dbcsr_distribution_new, dbcsr_distribution_hold,&
            dbcsr_distribution_release, dbcsr_distribution_init
  PUBLIC :: dbcsr_distribution_mp, dbcsr_distribution_processor,&
            dbcsr_distribution_nrows, dbcsr_distribution_ncols,&
            dbcsr_distribution_row_dist, dbcsr_distribution_col_dist,&
            dbcsr_distribution_nlocal_rows, dbcsr_distribution_nlocal_cols,&
            dbcsr_distribution_local_rows, dbcsr_distribution_local_cols,&
            dbcsr_distribution_thread_dist, dbcsr_distribution_has_threads,&
            dbcsr_distribution_make_threads, dbcsr_distribution_no_threads
  PUBLIC :: dbcsr_get_info, dbcsr_distribution,&
            dbcsr_get_matrix_type, dbcsr_get_data_type, dbcsr_get_replication_type,&
            dbcsr_row_block_sizes, dbcsr_col_block_sizes,&
            dbcsr_nblkrows_total, dbcsr_nblkcols_total, dbcsr_nfullrows_total,&
            dbcsr_nfullcols_total, dbcsr_nblkcols_local, dbcsr_nblkrows_local,&
            dbcsr_nfullrows_local, dbcsr_nfullcols_local,&
            dbcsr_max_row_size, dbcsr_max_col_size,&
            dbcsr_uses_special_memory, dbcsr_get_occupation,&
            dbcsr_is_real, dbcsr_is_complex,&
            dbcsr_name, dbcsr_may_be_dense,&
            dbcsr_use_mutable, dbcsr_wm_use_mutable, dbcsr_has_symmetry,&
            dbcsr_is_transposed
  PUBLIC :: dbcsr_get_data_size_used
  PUBLIC :: dbcsr_data_init, dbcsr_get_data, dbcsr_data_new, dbcsr_data_hold,&
            dbcsr_data_release, dbcsr_get_data_size,&
            dbcsr_data_set_pointer, dbcsr_data_area,&
            dbcsr_data_clear_pointer, dbcsr_data_set_2d_pointer,&
            dbcsr_data_clear_2d_pointer, dbcsr_data_ensure_size,&
            dbcsr_col_block_offsets, dbcsr_row_block_offsets
  PUBLIC :: dbcsr_get_num_blocks
  PUBLIC :: dbcsr_get_data_p, dbcsr_get_data_p_s, dbcsr_get_data_p_c,&
            dbcsr_get_data_p_d, dbcsr_get_data_p_z,&
            dbcsr_get_data_p_2d_s, dbcsr_get_data_p_2d_d,&
            dbcsr_get_data_p_2d_c, dbcsr_get_data_p_2d_z


  PUBLIC :: dbcsr_blk_row_size, dbcsr_blk_column_size

  PUBLIC :: dbcsr_allocate_matrix_array, dbcsr_array_new, dbcsr_array_put,&
            dbcsr_array_get, dbcsr_array_destroy, dbcsr_array_hold,&
            dbcsr_array_release
  PUBLIC :: dbcsr_destroy_array, dbcsr_destroy_image_dist

  PUBLIC :: dbcsr_mutable_init, dbcsr_mutable_new, dbcsr_mutable_destroy,&
            dbcsr_mutable_release, dbcsr_mutable_hold,&
            dbcsr_mutable_instantiated

  PUBLIC :: dbcsr_modify_lock, dbcsr_modify_unlock

  PUBLIC :: dbcsr_type_is_2d, dbcsr_type_2d_to_1d, dbcsr_type_1d_to_2d



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

  INTERFACE dbcsr_uses_special_memory
     MODULE PROCEDURE uses_special_memory_matrix, uses_special_memory_area
  END INTERFACE

  INTERFACE dbcsr_get_data
     MODULE PROCEDURE get_data_s, get_data_d, get_data_c, get_data_z
     MODULE PROCEDURE get_data_m_s, get_data_m_d, get_data_m_c, get_data_m_z
     MODULE PROCEDURE get_data_2d_s, get_data_2d_d, get_data_2d_c, get_data_2d_z
  END INTERFACE

  INTERFACE dbcsr_get_data_size
     MODULE PROCEDURE get_data_size_area, get_data_size_matrix
  END INTERFACE

  INTERFACE dbcsr_data_set_pointer
     MODULE PROCEDURE set_data_p_s, set_data_p_d, set_data_p_c, set_data_p_z
     MODULE PROCEDURE set_data_p_2d_s, set_data_p_2d_d,&
                      set_data_p_2d_c, set_data_p_2d_z
     MODULE PROCEDURE set_data_area_area
  END INTERFACE

  INTERFACE dbcsr_get_data_type
     MODULE PROCEDURE dbcsr_get_data_type_obj, data_get_data_type
  END INTERFACE

  INTERFACE dbcsr_is_real
     MODULE PROCEDURE dbcsr_is_real_obj
  END INTERFACE

  INTERFACE dbcsr_is_complex
     MODULE PROCEDURE dbcsr_is_complex_obj
  END INTERFACE

  !> \brief Encapsulates a scalar.
  INTERFACE dbcsr_scalar
     MODULE PROCEDURE dbcsr_scalar_s, dbcsr_scalar_d,&
                      dbcsr_scalar_c, dbcsr_scalar_z
  END INTERFACE

  INTERFACE dbcsr_get_data_p
     MODULE PROCEDURE dbcsr_get_data_c_s, dbcsr_get_data_c_c,&
                      dbcsr_get_data_c_d, dbcsr_get_data_c_z
  END INTERFACE

  INTERFACE dbcsr_blk_row_size
     MODULE PROCEDURE dbcsr_blk_row_size_type, dbcsr_blk_row_size_obj
  END INTERFACE

  INTERFACE dbcsr_blk_column_size
     MODULE PROCEDURE dbcsr_blk_column_size_type, dbcsr_blk_column_size_obj
  END INTERFACE

  INTERFACE pointer_rank_remap
     MODULE PROCEDURE pointer_s_rank_remap2, pointer_d_rank_remap2,&
                      pointer_c_rank_remap2, pointer_z_rank_remap2
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


#define DBG IF (dbg) WRITE(*,*)routineN//" ",
#define DBGV IF (bcsr_verbose) WRITE(*,*)routineN//" ",
#define DBGI IF (info) WRITE(*,*)routineN//" ",
#define DEBUG_HEADER        LOGICAL :: dbg, info
#define DEBUG_BODY        dbg = .FALSE. ; CALL cp_dbcsr_set_debug(dbg, info=info)
#define DEBUG_BODYY       dbg = .TRUE. ; CALL cp_dbcsr_set_debug(dbg, info=info)



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

    matrix%initialized = dbcsr_magic_number
    matrix%valid = .FALSE.
    matrix%refcount = 0
    ! Nullifies all pointers.
    NULLIFY (matrix%index, matrix%row_p, matrix%col_i,&
         matrix%blk_p)
    CALL dbcsr_data_init (matrix%data_area)
    CALL dbcsr_distribution_init (matrix%dist)
    CALL array_nullify (matrix%row_blk_size)
    CALL array_nullify (matrix%col_blk_size)
    CALL array_nullify (matrix%row_blk_offset)
    CALL array_nullify (matrix%col_blk_offset)
    NULLIFY (matrix%wms)
    NULLIFY (matrix%predistributed)
  END SUBROUTINE dbcsr_init_type

! *****************************************************************************
!> \brief Initializes a DBCSR matrix but does not allocate any memory.
!> \param[out] matrix         uninitialized matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_init_obj (matrix)
    TYPE(dbcsr_obj), INTENT(OUT)             :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_obj', &
      routineP = moduleN//':'//routineN

!TYPE(dbcsr_error_type), INTENT(INOUT)   :: error
!   ---------------------------------------------------------------------------

    CALL dbcsr_init(matrix%m)
    matrix%m%initialized = 0
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
    matrix%m%initialized = 0
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
    !valid_index = .FALSE.
    !IF (matrix%m%initialized .EQ. dbcsr_magic_number) THEN
    !   IF (ASSOCIATED (matrix%m%row_p)) THEN
    !      valid_index = SIZE (matrix%m%row_p) .GT. 0
    !   ENDIF
    !ELSE
    !   valid_index = .FALSE.
    !ENDIF
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

    initialized = matrix%initialized .EQ. dbcsr_magic_number
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

    CALL dbcsr_assert (matrix%m%initialized .EQ. dbcsr_magic_number,&
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
    IF (matrix%m%initialized .EQ. dbcsr_magic_number) THEN
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
       IF (matrix%m%initialized .EQ. dbcsr_magic_number) THEN
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
    TYPE(dbcsr_error_type)                   :: error
    LOGICAL, INTENT(IN), OPTIONAL            :: force

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_destroy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: force_all

!   ---------------------------------------------------------------------------

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
          IF (matrix%m%special_memory) THEN
             CALL MP_DEALLOCATE(matrix%m%index)
          ELSE
             DEALLOCATE(matrix%m%index)
          ENDIF
          !CALL cp_error_init(my_error)
          IF (ASSOCIATED (matrix%m%predistributed)) THEN
             CALL dbcsr_assert (ASSOCIATED (matrix%m%predistributed%mats),&
                  dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "Predistributed data should exist if allocated",__LINE__,error)
             CALL dbcsr_destroy_array (matrix%m%predistributed,error)
             CALL dbcsr_destroy_image_dist (matrix%m%predistributed%image_dist)
             DEALLOCATE (matrix%m%predistributed)
          ENDIF
          NULLIFY (matrix%m%predistributed)
          CALL dbcsr_data_release (matrix%m%data_area)
          CALL array_release (matrix%m%row_blk_size)
          CALL array_release (matrix%m%col_blk_size)
          CALL array_release (matrix%m%row_blk_offset)
          CALL array_release (matrix%m%col_blk_offset)
          CALL dbcsr_distribution_release(matrix%m%dist)
          matrix%m%valid = .FALSE.
!$        CALL OMP_DESTROY_LOCK (matrix%m%modification_lock)
          CALL dbcsr_init (matrix%m)
       ENDIF
    ENDIF
  END SUBROUTINE dbcsr_destroy


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
       CALL mp_cart_sub (tmp_group, (/.FALSE., .TRUE./), mp_env%mp%prow_group)
       CALL mp_cart_sub (tmp_group, (/.TRUE., .FALSE./), mp_env%mp%pcol_group)
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

    INTEGER                                  :: i, mypcoor, seq
    INTEGER, DIMENSION(:), POINTER           :: dd, ld

!   ---------------------------------------------------------------------------

    ALLOCATE (dist%d)
    dist%d%refcount = 1
    dist%d%row_dist = row_dist
    CALL array_hold (dist%d%row_dist)
    dist%d%col_dist = col_dist
    CALL array_hold (dist%d%col_dist)
    dist%d%mp_env = mp_env
    CALL dbcsr_mp_hold (dist%d%mp_env)
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
    CALL array_nullify (dist%d%thread_dist)
    !WRITE(*,*)' row_dist', array_data(dist%d%row_dist)
    !WRITE(*,*)' col_dist', array_data(dist%d%col_dist)
    !WRITE(*,*)' loc_rows', array_data(dist%d%local_rows)
    !WRITE(*,*)' loc_cols', array_data(dist%d%local_cols)
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

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (dist%d)) THEN
       dist%d%refcount = dist%d%refcount - 1
       IF (dist%d%refcount .EQ. 0) THEN
          CALL array_release (dist%d%row_dist)
          CALL array_release (dist%d%col_dist)
          CALL array_release (dist%d%local_rows)
          CALL array_release (dist%d%local_cols)
          CALL dbcsr_mp_release (dist%d%mp_env)
          CALL array_release (dist%d%thread_dist)
          DEALLOCATE (dist%d)
          CALL dbcsr_distribution_init (dist)
       ENDIF
    ENDIF
  END SUBROUTINE dbcsr_distribution_release

  SUBROUTINE dbcsr_distribution_init (dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist

    NULLIFY (dist%d)
  END SUBROUTINE dbcsr_distribution_init

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

    processor = dist%d%mp_env%mp%pgrid(dist%d%row_dist%low%data(row),&
         dist%d%col_dist%low%data(col))
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

    has_thread_dist = array_exists (dist%d%thread_dist)
  END FUNCTION dbcsr_distribution_has_threads

! *****************************************************************************
!> \brief Creates a distribution for threads
!> \param[in,out] dist   Add thread distribution to this distribution
!> \param[in] row_sizes  (optional) row block sizes
! *****************************************************************************
  SUBROUTINE dbcsr_distribution_make_threads(dist, row_sizes)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT), VOLATILE                :: dist
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: row_sizes

    INTEGER                                  :: first, ithread, last, nlrows, &
                                                nrows, nthreads, rpt, &
                                                size_local, size_sum, t
    INTEGER, DIMENSION(:), POINTER           :: lrows, td

!   ---------------------------------------------------------------------------

    ithread = 0 ; nthreads = 1
    !$ nthreads = OMP_GET_NUM_THREADS () ; ithread = OMP_GET_THREAD_NUM()
    nrows = dbcsr_distribution_nrows (dist)
    nlrows = dbcsr_distribution_nlocal_rows (dist)
    lrows => array_data (dbcsr_distribution_local_rows (dist))
    !
!$OMP SINGLE
    CALL array_release (dist%d%thread_dist)
    ALLOCATE (td(nrows))
    CALL array_new (dist%d%thread_dist, td, gift=.TRUE.)
!$OMP END SINGLE
    td => array_data (dist%d%thread_dist)
    rpt = CEILING (REAL(nrows) / REAL(nthreads))
    t = ithread
    first = MIN (t*rpt+1, nrows)
    last = MIN ((t+1)*rpt, nrows)
    IF(first.GE.1.AND.last.LE.nrows) THEN
       td(first:last) = t
    ENDIF
!$OMP BARRIER
    ! Calculate row distribution
    IF (PRESENT (row_sizes)) THEN
       size_local = 0
       DO t = 1, nlrows
          size_local = size_local + row_sizes(lrows(t))
       END DO
       rpt = CEILING (REAL(size_local) / REAL(nthreads))
       first = MIN (ithread*rpt, nlrows)
       last = MIN ((ithread+1)*rpt, nlrows)
       size_sum = 0
       DO t = 1, nlrows
          IF (size_sum .GE. first .AND. size_sum .LT. last) THEN
             td(lrows(t)) = ithread
          ELSEIF (size_sum .GE. size_local) THEN
             td(lrows(t)) = nthreads-1
          ENDIF
          size_sum = size_sum + row_sizes(lrows(t))
          IF (size_sum .LT. last) EXIT
       ENDDO
    ELSE
       size_local = nlrows
       rpt = CEILING (REAL(nlrows) / REAL(nthreads))
       first = MIN (ithread*rpt+1, nlrows)
       last = MIN ((ithread+1)*rpt, nlrows)
       IF(first.GE.1.AND.last.LE.nlrows) THEN
          FORALL (t = first:last)
             td(lrows(t)) = ithread
          END FORALL
       ENDIF
    ENDIF
!$OMP BARRIER
  END SUBROUTINE dbcsr_distribution_make_threads

  SUBROUTINE dbcsr_distribution_no_threads(dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT), VOLATILE                :: dist

!   ---------------------------------------------------------------------------

!$OMP MASTER
    CALL array_release (dist%d%thread_dist)
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
    CHARACTER, OPTIONAL                      :: matrix_type, data_type

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_info', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (matrix%m%initialized.EQ.dbcsr_magic_number,&
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
       ! Doesn't check all valid types
       matrix_type = '0' ! Invalid
       IF (matrix%m%transpose) THEN
          matrix_type = dbcsr_type_transposed
       ELSEIF (matrix%m%symmetry) THEN
          IF (.NOT. matrix%m%negate_real .AND. matrix%m%negate_imaginary) THEN
             matrix_type = dbcsr_type_conjugate
          ELSEIF (matrix%m%negate_real .AND. .NOT. matrix%m%negate_imaginary) THEN
             matrix_type = dbcsr_type_antisymmetric
          ELSEIF (.NOT. matrix%m%negate_real .AND. .NOT. matrix%m%negate_imaginary) THEN
             matrix_type = dbcsr_type_symmetric
          ENDIF
       ELSE
          matrix_type = dbcsr_type_normal
       ENDIF
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



! Data type transformations
  FUNCTION data_type_2d_from_1d (type_1d) RESULT (type_2d)
    CHARACTER, INTENT(in)                    :: type_1d
    CHARACTER                                :: type_2d

    CHARACTER(len=*), PARAMETER :: routineN = 'data_type_2d_from_1d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

    SELECT CASE (type_1d)
    CASE (dbcsr_type_real_4)
       type_2d = dbcsr_type_real_4_2d
    CASE (dbcsr_type_real_8)
       type_2d = dbcsr_type_real_8_2d
    CASE (dbcsr_type_complex_4)
       type_2d = dbcsr_type_complex_4_2d
    CASE (dbcsr_type_complex_8)
       type_2d = dbcsr_type_complex_8_2d
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
            routineN, "Invalid data type.",__LINE__,error)
    END SELECT
  END FUNCTION data_type_2d_from_1d

! Data type transformations
  FUNCTION data_type_1d_from_2d (type_2d) RESULT (type_1d)
    CHARACTER, INTENT(IN)                    :: type_2d
    CHARACTER                                :: type_1d

    CHARACTER(len=*), PARAMETER :: routineN = 'data_type_1d_from_2d', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

    SELECT CASE (type_2d)
    CASE (dbcsr_type_real_4_2d)
       type_1d = dbcsr_type_real_4
    CASE (dbcsr_type_real_8_2d)
       type_1d = dbcsr_type_real_8
    CASE (dbcsr_type_complex_4_2d)
       type_1d = dbcsr_type_complex_4
    CASE (dbcsr_type_complex_8_2d)
       type_1d = dbcsr_type_complex_8
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
            routineN, "Invalid data type.",__LINE__,error)
    END SELECT
  END FUNCTION data_type_1d_from_2d


! *****************************************************************************
!> \brief Initializes a data area
!> \param[inout] area         data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_init (area)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area

    NULLIFY (area%d)
  END SUBROUTINE dbcsr_data_init
! *****************************************************************************
!> \brief Initializes a data area and all the actual data pointers
!> \param[inout] area         data area
!> \param[in] data_type       select data type to use
!> \param[in] data_size       (optional) allocate this much data
!> \param[in] data_size2      (optional) second dimension data size
!> \param[in] special_memory  (optional) whether to use specially-allocated
!>                            memory, default is False
! *****************************************************************************
  SUBROUTINE dbcsr_data_new (area, data_type, data_size, data_size2,&
       special_memory)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    CHARACTER, INTENT(IN)                    :: data_type
    INTEGER, INTENT(IN), OPTIONAL            :: data_size, data_size2
    LOGICAL, INTENT(IN), OPTIONAL            :: special_memory

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_new', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (.NOT. ASSOCIATED (area%d)) THEN
       ALLOCATE (area%d)
    ENDIF
    NULLIFY (area%d%r_sp)
    NULLIFY (area%d%r_dp)
    NULLIFY (area%d%c_sp)
    NULLIFY (area%d%c_dp)
    NULLIFY (area%d%r2_sp)
    NULLIFY (area%d%r2_dp)
    NULLIFY (area%d%c2_sp)
    NULLIFY (area%d%c2_dp)
    area%d%refcount = 1
    IF (PRESENT (special_memory)) THEN
       area%d%special = special_memory
    ELSE
       area%d%special = .FALSE.
    ENDIF
    IF (PRESENT (data_size)) THEN
       IF (area%d%special) THEN
          SELECT CASE (data_type)
          CASE (dbcsr_type_real_4)
             CALL mp_allocate (area%d%r_sp, data_size)
          CASE (dbcsr_type_real_8)
             CALL mp_allocate (area%d%r_dp, data_size)
          CASE (dbcsr_type_complex_4)
             CALL mp_allocate (area%d%c_sp, data_size)
          CASE (dbcsr_type_complex_8)
             CALL mp_allocate (area%d%c_dp, data_size)
          CASE (dbcsr_type_real_8_2d, dbcsr_type_real_4_2d,&
               dbcsr_type_complex_8_2d, dbcsr_type_complex_4_2d)
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
                  routineN, "Can not use special memory with 2D data areas.",__LINE__,error)
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
                  routineN, "Invalid data type.",__LINE__,error)
          END SELECT
       ELSE
          SELECT CASE (data_type)
          CASE (dbcsr_type_real_8)
             ALLOCATE (area%d%r_dp(data_size))
          CASE (dbcsr_type_real_4)
             ALLOCATE (area%d%r_sp(data_size))
          CASE (dbcsr_type_complex_8)
             ALLOCATE (area%d%c_dp(data_size))
          CASE (dbcsr_type_complex_4)
             ALLOCATE (area%d%c_sp(data_size))
          CASE (dbcsr_type_real_8_2d)
             CALL dbcsr_assert (PRESENT(data_size2), dbcsr_fatal_level,&
                  dbcsr_caller_error, routineN,&
                  "Must specify either both or no sizes for 2D data",__LINE__,error)
             ALLOCATE (area%d%r2_dp(data_size, data_size2))
          CASE (dbcsr_type_real_4_2d)
             CALL dbcsr_assert (PRESENT(data_size2), dbcsr_fatal_level,&
                  dbcsr_caller_error, routineN,&
                  "Must specify either both or no sizes for 2D data",__LINE__,error)
             ALLOCATE (area%d%r2_sp(data_size, data_size2))
          CASE (dbcsr_type_complex_8_2d)
             CALL dbcsr_assert (PRESENT(data_size2), dbcsr_fatal_level,&
                  dbcsr_caller_error, routineN,&
                  "Must specify either both or no sizes for 2D data",__LINE__,error)
             ALLOCATE (area%d%c2_dp(data_size, data_size2))
          CASE (dbcsr_type_complex_4_2d)
             CALL dbcsr_assert (PRESENT(data_size2), dbcsr_fatal_level,&
                  dbcsr_caller_error, routineN,&
                  "Must specify either both or no sizes for 2D data",__LINE__,error)
             ALLOCATE (area%d%c2_sp(data_size, data_size2))
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
                  routineN, "Invalid data type.",__LINE__,error)
          END SELECT
       ENDIF
    ENDIF
    area%d%data_type = data_type
  END SUBROUTINE dbcsr_data_new

! *****************************************************************************
!> \brief Returns data type of a data area
!> \param[in] area         data area
!> \result data_type       data type of the data area
! *****************************************************************************
  PURE FUNCTION data_get_data_type (area) RESULT (data_type)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    CHARACTER                                :: data_type

    data_type = area%d%data_type
  END FUNCTION data_get_data_type

! *****************************************************************************
!> \brief Removes a reference and/or clears the data area.
!> \param[inout] area         data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_release (area)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_release', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (.NOT. ASSOCIATED (area%d)) THEN
       RETURN
    ENDIF
    CALL dbcsr_assert (area%d%refcount .GT. 0, dbcsr_warning_level, dbcsr_caller_error,&
         routineN, "Data seems to be unreferenced.",__LINE__,error)
    IF (area%d%refcount .GT. 0) THEN
       area%d%refcount = area%d%refcount - 1
    ELSE
       RETURN
    ENDIF
    ! If we're releasing the last reference, then free the memory.
    IF (area%d%refcount .EQ. 0) THEN
       IF (area%d%special) THEN
          SELECT CASE (area%d%data_type)
          CASE (dbcsr_type_real_4)
             CALL mp_deallocate (area%d%r_sp)
             NULLIFY (area%d%r_sp)
          CASE (dbcsr_type_real_8)
             CALL mp_deallocate (area%d%r_dp)
             NULLIFY (area%d%r_dp)
          CASE (dbcsr_type_complex_4)
             CALL mp_deallocate (area%d%c_sp)
             NULLIFY (area%d%c_sp)
          CASE (dbcsr_type_complex_8)
             CALL mp_deallocate (area%d%c_dp)
             NULLIFY (area%d%c_dp)
          CASE (dbcsr_type_real_8_2d, dbcsr_type_real_4_2d,&
               dbcsr_type_complex_8_2d, dbcsr_type_complex_4_2d)
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_internal_error,&
                  routineN, "How could 2D data areas use special memory?",__LINE__,error)
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
                  routineN, "Invalid data type.",__LINE__,error)
          END SELECT
       ELSE
          SELECT CASE (area%d%data_type)
          CASE (dbcsr_type_real_8)
             IF (ASSOCIATED (area%d%r_dp)) &
                  DEALLOCATE (area%d%r_dp)
             NULLIFY (area%d%r_dp)
          CASE (dbcsr_type_real_4)
             IF (ASSOCIATED (area%d%r_sp)) &
                  DEALLOCATE (area%d%r_sp)
             NULLIFY (area%d%r_sp)
          CASE (dbcsr_type_complex_8)
             IF (ASSOCIATED (area%d%c_dp)) &
                  DEALLOCATE (area%d%c_dp)
             NULLIFY (area%d%c_dp)
          CASE (dbcsr_type_complex_4)
             IF (ASSOCIATED (area%d%c_sp)) &
                  DEALLOCATE (area%d%c_sp)
             NULLIFY (area%d%c_sp)
          CASE (dbcsr_type_real_8_2d)
             IF (ASSOCIATED (area%d%r2_dp)) &
                  DEALLOCATE (area%d%r2_dp)
             NULLIFY (area%d%r2_dp)
          CASE (dbcsr_type_real_4_2d)
             IF (ASSOCIATED (area%d%r2_sp)) &
                  DEALLOCATE (area%d%r2_sp)
             NULLIFY (area%d%r2_sp)
          CASE (dbcsr_type_complex_8_2d)
             IF (ASSOCIATED (area%d%c2_dp)) &
                  DEALLOCATE (area%d%c2_dp)
             NULLIFY (area%d%c2_dp)
          CASE (dbcsr_type_complex_4_2d)
             IF (ASSOCIATED (area%d%c2_sp)) &
                  DEALLOCATE (area%d%c2_sp)
             NULLIFY (area%d%c2_sp)
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
                  routineN, "Invalid data type.",__LINE__,error)
          END SELECT
       ENDIF
       DEALLOCATE (area%d)
       NULLIFY (area%d)
    ENDIF
  END SUBROUTINE dbcsr_data_release

! *****************************************************************************
!> \brief Clears pointers from the data area.
!> \param[inout] area         data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_clear_pointer (area)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_clear_pointer', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    IF (.NOT. ASSOCIATED (area%d)) THEN
       RETURN
    ENDIF
    CALL dbcsr_assert (area%d%refcount .GT. 0, dbcsr_warning_level, dbcsr_caller_error,&
         routineN, "Data seems to be unreferenced.",__LINE__,error)
    SELECT CASE (area%d%data_type)
    CASE (dbcsr_type_real_4)
       NULLIFY (area%d%r_sp)
    CASE (dbcsr_type_real_8)
       NULLIFY (area%d%r_dp)
    CASE (dbcsr_type_complex_4)
       NULLIFY (area%d%c_sp)
    CASE (dbcsr_type_complex_8)
       NULLIFY (area%d%c_dp)
    CASE (dbcsr_type_real_8_2d)
       NULLIFY (area%d%r2_dp)
    CASE (dbcsr_type_real_4_2d)
       NULLIFY (area%d%r2_sp)
    CASE (dbcsr_type_complex_8_2d)
       NULLIFY (area%d%c2_dp)
    CASE (dbcsr_type_complex_4_2d)
       NULLIFY (area%d%c2_sp)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
            routineN, "Invalid data type.",__LINE__,error)
    END SELECT
  END SUBROUTINE dbcsr_data_clear_pointer



! *****************************************************************************
!> \brief Registers another use of the data area
!> \param[inout] area         data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_hold (area)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_hold', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    CALL dbcsr_assert (ASSOCIATED (area%d), dbcsr_caller_error, dbcsr_warning_level,&
         routineN, "Can't hold an empty data area.",__LINE__,error)
    CALL dbcsr_assert (area%d%refcount .GT. 0, dbcsr_caller_error, dbcsr_warning_level,&
         routineN, "Should not hold an area with zero references.",__LINE__,error)
    IF (.NOT. ASSOCIATED (area%d)) THEN
       RETURN
    ENDIF
    area%d%refcount = area%d%refcount + 1
  END SUBROUTINE dbcsr_data_hold

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
!> \brief Returns the matrix type
!> \param matrix              query this matrix
!> \param matrix_type         matrix_type (see dbcsr_types.F for possible
!>                            values)
! *****************************************************************************
  PURE FUNCTION dbcsr_get_matrix_type (matrix) RESULT (matrix_type)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    CHARACTER                                :: matrix_type

! Doesn't check all valid types

    matrix_type = '0' ! Invalid
    IF (matrix%m%transpose) THEN
       matrix_type = dbcsr_type_transposed
    ELSEIF (matrix%m%symmetry) THEN
       IF (.NOT. matrix%m%negate_real .AND. matrix%m%negate_imaginary) THEN
          matrix_type = dbcsr_type_conjugate
       ELSEIF (matrix%m%negate_real .AND. .NOT. matrix%m%negate_imaginary) THEN
          matrix_type = dbcsr_type_antisymmetric
       ELSEIF (.NOT. matrix%m%negate_real .AND. .NOT. matrix%m%negate_imaginary) THEN
          matrix_type = dbcsr_type_symmetric
       ENDIF
    ELSE
       matrix_type = dbcsr_type_normal
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
  PURE FUNCTION dbcsr_is_transposed (matrix) RESULT (transposed)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL                                  :: transposed

    transposed = matrix%m%transpose
  END FUNCTION dbcsr_is_transposed

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
    CHARACTER                                :: data_type

    data_type = matrix%m%data_type
  END FUNCTION dbcsr_get_data_type_obj

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
!> \brief Returns whether the matrix uses specially-allocated memory
!> \param[in] matrix          query this matrix
!> \param[out] uses_special   whether the matrix uses specially allocated
!>                            memory
! *****************************************************************************
  PURE FUNCTION uses_special_memory_matrix (matrix) RESULT (uses_special)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL                                  :: uses_special

    uses_special = matrix%m%special_memory
  END FUNCTION uses_special_memory_matrix


! *****************************************************************************
!> \brief Returns whether the data area uses speciall-allocated memory
!> \param[in] data_area       query this data area
!> \param[out] uses_special   whether the data area uses specially allocated
!>                            memory
! *****************************************************************************
  PURE FUNCTION uses_special_memory_area (area) RESULT (uses_special)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    LOGICAL                                  :: uses_special

    IF (ASSOCIATED (area%d)) THEN
       uses_special = area%d%special
    ELSE
       uses_special = .FALSE.
    ENDIF
  END FUNCTION uses_special_memory_area


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

    TYPE(dbcsr_error_type)                   :: error

    data_size = 0
    IF (ASSOCIATED (area%d)) THEN
       SELECT CASE (area%d%data_type)
          CASE (dbcsr_type_real_8)
             IF (ASSOCIATED (area%d%r_dp))&
                  data_size = SIZE (area%d%r_dp)
          CASE (dbcsr_type_real_4)
             IF (ASSOCIATED (area%d%r_sp))&
                  data_size = SIZE (area%d%r_sp)
          CASE (dbcsr_type_complex_8)
             IF (ASSOCIATED (area%d%c_dp))&
                  data_size = SIZE (area%d%c_dp)
          CASE (dbcsr_type_complex_4)
             IF (ASSOCIATED (area%d%c_sp))&
                  data_size = SIZE (area%d%c_sp)
          CASE (dbcsr_type_real_8_2d)
             IF (ASSOCIATED (area%d%r2_dp))&
                  data_size = SIZE (area%d%r2_dp)
          CASE (dbcsr_type_real_4_2d)
             IF (ASSOCIATED (area%d%r2_sp))&
                  data_size = SIZE (area%d%r2_sp)
          CASE (dbcsr_type_complex_8_2d)
             IF (ASSOCIATED (area%d%c2_dp))&
                  data_size = SIZE (area%d%c2_dp)
          CASE (dbcsr_type_complex_4_2d)
             IF (ASSOCIATED (area%d%c2_sp))&
                  data_size = SIZE (area%d%c2_sp)
          CASE default
             CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
                  routineN, "Incorrect data type",__LINE__,error)
          END SELECT
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_warning_level, dbcsr_caller_error, routineN,&
            "Uninitialized data area",__LINE__,error)
       data_size = 0
    ENDIF
  END FUNCTION get_data_size_area


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
!> \brief Count actual data storage used for matrix
!> \param[in] matrix          Count data of this matrix
!> \param[in,out] error       Error
!> \result data_size          Data size used by matrix
! *****************************************************************************
  FUNCTION dbcsr_get_data_size_used (matrix) RESULT (data_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: data_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_size_used', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_p, col, nze, row
    INTEGER, DIMENSION(:), POINTER           :: col_blk_sizes, row_blk_sizes

!type(dbcsr_iterator_type) :: iter
!   ---------------------------------------------------------------------------

    row_blk_sizes => array_data (dbcsr_row_block_sizes (matrix))
    col_blk_sizes => array_data (dbcsr_col_block_sizes (matrix))
    data_size = 0
    DO row = 1, matrix%m%nblkrows_total
       DO blk = matrix%m%row_p(row)+1, matrix%m%row_p(row+1)
          col = matrix%m%col_i(blk)
          blk_p = ABS (matrix%m%blk_p(blk))
          nze = row_blk_sizes(row) * col_blk_sizes(col)
          data_size = MAX (data_size, blk_p + nze-1)
       ENDDO
    ENDDO
    !call dbcsr_iterator_start (iter, matrix)
    !data_size = 0
    !do while (dbcsr_iterator_blocks_left (iter))
    !   call dbcsr_iterator_next_block (iter, row, col, blk, &
    !        blk_p=blk_p, row_size=row_size, col_size=col_size)
    !   data_size = MAX (data_size, blk_p + row_size*col_size-1)
    !enddo
    !call dbcsr_iterator_stop (iter)
  END FUNCTION dbcsr_get_data_size_used


! *****************************************************************************
!> \brief Returns an encapsulated scalar "1"
!> \param[in] data_type       use the data type
!> \result one                enpsulated value of one
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_one (data_type) RESULT (one)
    CHARACTER, INTENT(IN)                    :: data_type
    TYPE(dbcsr_scalar_type)                  :: one

    one%data_type = data_type
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       one%r_sp = 1.0_real_4
    CASE (dbcsr_type_real_8)
       one%r_dp = 1.0_real_8
    CASE (dbcsr_type_complex_4)
       one%c_sp = CMPLX(1.0, 0.0, real_4)
    CASE (dbcsr_type_complex_8)
       one%c_dp = CMPLX(1.0, 0.0, real_8)
    END SELECT
  END FUNCTION dbcsr_scalar_one

! *****************************************************************************
!> \brief Returns an encapsulated scalar "0"
!> \param[in] data_type       use the data type
!> \result zero               enpsulated value of zero
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_zero (data_type) RESULT (zero)
    CHARACTER, INTENT(IN)                    :: data_type
    TYPE(dbcsr_scalar_type)                  :: zero

    zero%data_type = data_type
    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       zero%r_sp = 0.0_real_4
    CASE (dbcsr_type_real_8)
       zero%r_dp = 0.0_real_8
    CASE (dbcsr_type_complex_4)
       zero%c_sp = CMPLX(0.0, 0.0, real_4)
    CASE (dbcsr_type_complex_8)
       zero%c_dp = CMPLX(0.0, 0.0, real_8)
    END SELECT
  END FUNCTION dbcsr_scalar_zero


! *****************************************************************************
!> \brief Returns whether an encapsulated scalar is equal to another value
!> \param[in] s1                    one value
!> \param[in] s2                    another value
!> \result are_equal                whether values are equal
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_are_equal (s1, s2) RESULT (are_equal)
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: s1, s2
    LOGICAL                                  :: are_equal

    IF (s1%data_type .NE. s2%data_type) THEN
       are_equal = .FALSE.
    ELSE
       SELECT CASE (s1%data_type)
       CASE (dbcsr_type_real_4)
          are_equal = s1%r_sp .EQ. s2%r_sp
       CASE (dbcsr_type_real_8)
          are_equal = s1%r_dp .EQ. s2%r_dp
       CASE (dbcsr_type_complex_4)
          are_equal = s1%c_sp .EQ. s2%c_sp
       CASE (dbcsr_type_complex_8)
          are_equal = s1%c_dp .EQ. s2%c_dp
       CASE default
          are_equal = .FALSE.
       END SELECT
    ENDIF
  END FUNCTION dbcsr_scalar_are_equal

! *****************************************************************************
!> \brief Returns an encapsulated scalar as a negation of the given
!> \param[in] v                    given value
!> \result negated                 negated value
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_negative (s) RESULT (negated)
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: s
    TYPE(dbcsr_scalar_type)                  :: negated

    negated%data_type = s%data_type
    SELECT CASE (s%data_type)
    CASE (dbcsr_type_real_4)
       negated%r_sp = -s%r_sp
    CASE (dbcsr_type_real_8)
       negated%r_dp = -s%r_dp
    CASE (dbcsr_type_complex_4)
       negated%c_sp = -s%c_sp
    CASE (dbcsr_type_complex_8)
       negated%c_dp = -s%c_dp
    CASE default
       negated = dbcsr_scalar_zero (s%data_type)
    END SELECT
  END FUNCTION dbcsr_scalar_negative



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
!> \brief Returns the occupation of the matrix
!> \param matrix     matrix from which to get the occupation
! *****************************************************************************
  FUNCTION dbcsr_get_occupation (matrix) RESULT (occupation)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    REAL(KIND=real_8)                        :: occupation

    INTEGER                                  :: nfullcols, nfullrows
    INTEGER, DIMENSION(:), POINTER           :: row_blk_size

    occupation = REAL(matrix%m%nze,real_8)
    CALL mp_sum(occupation,dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))

    nfullrows = dbcsr_nfullrows_total(matrix)
    nfullcols = dbcsr_nfullcols_total(matrix)

    row_blk_size => array_data (matrix%m%row_blk_size)

    IF(nfullrows.NE.0.AND.nfullcols.NE.0) THEN
       IF(dbcsr_has_symmetry (matrix)) THEN
          occupation = 2.0_real_8 * occupation / &
               ( REAL(nfullrows,real_8) * REAL(nfullrows+1,real_8) + &
               SUM( REAL(row_blk_size,real_8) * REAL(row_blk_size-1,real_8) ) )
       ELSE
          occupation = occupation / ( REAL(nfullrows,real_8) * REAL(nfullcols,real_8) )
       ENDIF
    ELSE
       occupation = 0.0_real_8
    ENDIF
  END FUNCTION dbcsr_get_occupation



! *****************************************************************************
!> \brief Points data area data pointers to another data area
!>
!> Assumes that no memory will be lost when repointing the pointer in the data
!> area and that the area is initialized.
!> \param[in,out] area         data area to repoint
!> \param[in] rsize, csize     size of data area to point to
!> \param[in] pointee          data area to point to
!> \param[in] source_lb        (optional) point to this offset in pointee
! *****************************************************************************
  SUBROUTINE set_data_area_area (area, rsize, csize, pointee, source_lb)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, INTENT(IN)                      :: rsize, csize
    TYPE(dbcsr_data_obj), INTENT(IN)         :: pointee
    INTEGER, INTENT(IN), OPTIONAL            :: source_lb

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_area_area', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: dt1, dt2
    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: c_sp
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: c_dp
    INTEGER                                  :: bp, nze
    LOGICAL                                  :: compatible, pointee_is_2d, rmp
    REAL(KIND=real_4), DIMENSION(:), POINTER :: r_sp
    REAL(KIND=real_8), DIMENSION(:), POINTER :: r_dp
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    bp = 1 ; IF (PRESENT (source_lb)) bp = source_lb
    nze = rsize*csize
    dt1 = area%d%data_type
    dt2 = pointee%d%data_type
    compatible = dt1 .EQ. dt2&
         .OR. dt1 .EQ. dbcsr_type_real_4_2d .AND. dt2 .EQ. dbcsr_type_real_4&
         .OR. dt1 .EQ. dbcsr_type_real_8_2d .AND. dt2 .EQ. dbcsr_type_real_8&
         .OR. dt1 .EQ. dbcsr_type_complex_4_2d .AND. dt2 .EQ. dbcsr_type_complex_4&
         .OR. dt1 .EQ. dbcsr_type_complex_8_2d .AND. dt2 .EQ. dbcsr_type_complex_8
    CALL dbcsr_assert (compatible, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Can not point 1-d pointer to 2-d data",__LINE__,error)
    pointee_is_2d = dt2 .EQ. dbcsr_type_real_4_2d &
         .OR. dt2 .EQ. dbcsr_type_real_8_2d &
         .OR. dt2 .EQ. dbcsr_type_complex_4_2d &
         .OR. dt2 .EQ. dbcsr_type_complex_8_2d
    CALL dbcsr_assert (.NOT. PRESENT (source_lb) .OR. .NOT. pointee_is_2d, &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Lower bound specification not possible with 2-d data",__LINE__,error)
    ! Check if size is OK.
    CALL dbcsr_assert (bp, "GE", 1,&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Attempt to point out of bounds",__LINE__,error)
    CALL dbcsr_assert (bp + nze - 1, "LE", dbcsr_get_data_size (pointee),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Attempt to point out of bounds",__LINE__,error)
    ! There's a remap if the ranks are compatible but not equal.
    rmp = dt1 .NE. dt2
    SELECT CASE (dt2)
    CASE (dbcsr_type_real_4_2d)
       area%d%r2_sp => pointee%d%r2_sp(1:rsize,1:csize)
    CASE (dbcsr_type_real_4)
       IF (rmp) THEN
          r_sp => dbcsr_get_data_p_s (pointee, bp, bp+nze-1)
          CALL pointer_s_rank_remap2(area%d%r2_sp, rsize, csize,&
               r_sp)
       ELSE
          area%d%r_sp => dbcsr_get_data_p_s (pointee, bp, bp+nze-1)
       ENDIF
    CASE (dbcsr_type_real_8_2d)
       area%d%r2_dp => pointee%d%r2_dp(1:rsize,1:csize)
    CASE (dbcsr_type_real_8)
       IF (rmp) THEN
          r_dp => dbcsr_get_data_p_d (pointee, bp, bp+nze-1)
          CALL pointer_d_rank_remap2(area%d%r2_dp, rsize, csize,&
               r_dp)
       ELSE
          area%d%r_dp => dbcsr_get_data_p_d (pointee, bp, bp+nze-1)
       ENDIF
    CASE (dbcsr_type_complex_4_2d)
       area%d%c2_sp => pointee%d%c2_sp(1:rsize,1:csize)
    CASE (dbcsr_type_complex_4)
       IF (rmp) THEN
          c_sp => dbcsr_get_data_p_c (pointee, bp, bp+nze-1)
          CALL pointer_c_rank_remap2(area%d%c2_sp, rsize, csize,&
               c_sp)
       ELSE
          area%d%c_sp => dbcsr_get_data_p_c (pointee, bp, bp+nze-1)
       ENDIF
    CASE (dbcsr_type_complex_8_2d)
       area%d%c2_dp => pointee%d%c2_dp(1:rsize,1:csize)
    CASE (dbcsr_type_complex_8)
       IF (rmp) THEN
          c_dp => dbcsr_get_data_p_z (pointee, bp, bp+nze-1)
          CALL pointer_z_rank_remap2(area%d%c2_dp, rsize, csize,&
               c_dp)
       ELSE
          area%d%c_dp => dbcsr_get_data_p_z (pointee, bp, bp+nze-1)
       ENDIF
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
  END SUBROUTINE set_data_area_area


! *****************************************************************************
!> \brief Ensures a minimum size of a previously-setup data area.
!>
!> The data area must have been previously setup with dbcsr_data_new.
!> \param[inout] area         data area
!> \param[in] data_size       (optional) allocate this much data
!> \param[in] nocopy          (optional) do not keep potentially existing data,
!>                            default is to keep it
!> \param[in] zero_pad        (optional) pad new data with zeros
!> \param[in] factor          (optional) increase size by this factor
! *****************************************************************************
  SUBROUTINE dbcsr_data_ensure_size (area, data_size, nocopy, zero_pad, factor, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, INTENT(IN), OPTIONAL            :: data_size
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, zero_pad
    REAL, INTENT(IN), OPTIONAL               :: factor
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_ensure_size', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: current_size, error_handler
    LOGICAL                                  :: nocp, pad

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert(ASSOCIATED (area%d), dbcsr_fatal_level, dbcsr_caller_error,&
         routineN, "Data area must be setup.",__LINE__,error)
    current_size = dbcsr_get_data_size (area)
    IF (current_size .GT. 1 .AND. current_size .GE. data_size) THEN
       CALL dbcsr_error_stop(error_handler, error)
       RETURN
    ENDIF
    !
    nocp = .FALSE.
    IF (PRESENT (nocopy)) nocp = nocopy
    pad = .FALSE.
    IF (PRESENT (zero_pad)) pad = zero_pad
    !
    SELECT CASE (area%d%data_type)
    CASE (dbcsr_type_real_8)
       IF (.NOT. ASSOCIATED (area%d%r_dp)) THEN
          IF (area%d%special) THEN
             CALL mp_allocate (area%d%r_dp, data_size)
          ELSE
             ALLOCATE (area%d%r_dp(data_size))
          ENDIF
          IF (pad) area%d%r_dp(:) = 0.0_real_8
       ELSE
          CALL ensure_array_size (area%d%r_dp, ub=data_size,&
               special=area%d%special, nocopy=nocp, zero_pad=zero_pad,&
               factor=factor,error=error)
       ENDIF
    CASE (dbcsr_type_real_4)
       IF (.NOT. ASSOCIATED (area%d%r_sp)) THEN
          IF (area%d%special) THEN
             CALL mp_allocate (area%d%r_sp, data_size)
          ELSE
             ALLOCATE (area%d%r_sp(data_size))
          ENDIF
          IF (pad) area%d%r_sp(:) = 0.0_real_4
       ELSE
          CALL ensure_array_size (area%d%r_sp, ub=data_size,&
               special=area%d%special, nocopy=nocp, zero_pad=zero_pad,&
               factor=factor,error=error)
       ENDIF
    CASE (dbcsr_type_complex_8)
       IF (.NOT. ASSOCIATED (area%d%c_dp)) THEN
          IF (area%d%special) THEN
             CALL mp_allocate (area%d%c_dp, data_size)
          ELSE
             ALLOCATE (area%d%c_dp(data_size))
          ENDIF
          IF (pad) area%d%c_dp(:) = 0.0_real_8
       ELSE
          CALL ensure_array_size (area%d%c_dp, ub=data_size,&
               special=area%d%special, nocopy=nocp, zero_pad=zero_pad,&
               factor=factor,error=error)
       ENDIF
    CASE (dbcsr_type_complex_4)
       IF (.NOT. ASSOCIATED (area%d%c_sp)) THEN
          IF (area%d%special) THEN
             CALL mp_allocate (area%d%c_sp, data_size)
          ELSE
             ALLOCATE (area%d%c_sp(data_size))
          ENDIF
          IF (pad) area%d%c_sp(:) = 0.0_real_4
       ELSE
          CALL ensure_array_size (area%d%c_sp, ub=data_size,&
               special=area%d%special, nocopy=nocp, zero_pad=zero_pad,&
               factor=factor,error=error)
       ENDIF
    CASE default
       CALL dbcsr_assert(.FALSE., dbcsr_failure_level,&
            dbcsr_unimplemented_error_nr, routineN,&
            "Invalid data type are supported",__LINE__,error)
    END SELECT
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_data_ensure_size



! *****************************************************************************
!> \brief Points a 2-d pointer to current 1-d data
!>
!> \param[in,out] area         data area to work on
!> \param[in] rowsize          row size of 2-d data
!> \param[in] colsize          (optional) column size of 2-d data
!> \param[in] offset           (optional) offset
! *****************************************************************************
  SUBROUTINE dbcsr_data_set_2d_pointer (area, rowsize, colsize, offset)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    INTEGER, INTENT(IN)                      :: rowsize
    INTEGER, INTENT(IN), OPTIONAL            :: colsize, offset

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_set_2d_pointer', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: dt1
    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: c_sp
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: c_dp
    INTEGER                                  :: bp, csize, nze
    LOGICAL                                  :: compatible
    REAL(KIND=real_4), DIMENSION(:), POINTER :: r_sp
    REAL(KIND=real_8), DIMENSION(:), POINTER :: r_dp
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    dt1 = area%d%data_type
    compatible = dt1 .EQ. dbcsr_type_real_4&
            .OR. dt1 .EQ. dbcsr_type_real_8&
            .OR. dt1 .EQ. dbcsr_type_complex_4&
            .OR. dt1 .EQ. dbcsr_type_complex_8
    CALL dbcsr_assert (compatible, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Must target 1-D data",__LINE__,error)
    IF (PRESENT (offset)) THEN
       bp = offset
       nze = dbcsr_get_data_size (area) - offset + 1
    ELSE
       bp = 1
       nze = dbcsr_get_data_size (area)
    ENDIF
    IF (PRESENT (colsize)) THEN
       csize = colsize
    ELSE
       csize = nze / rowsize
    ENDIF
    SELECT CASE (dt1)
    CASE (dbcsr_type_real_4)
       r_sp => dbcsr_get_data_p_s (area, bp, bp+nze-1)
       CALL pointer_s_rank_remap2(area%d%r2_sp, rowsize, csize,&
            r_sp)
    CASE (dbcsr_type_real_8)
       r_dp => dbcsr_get_data_p_d (area, bp, bp+nze-1)
       CALL pointer_d_rank_remap2(area%d%r2_dp, rowsize, csize,&
            r_dp)
    CASE (dbcsr_type_complex_4)
       c_sp => dbcsr_get_data_p_c (area, bp, bp+nze-1)
       CALL pointer_c_rank_remap2(area%d%c2_sp, rowsize, csize,&
            c_sp)
    CASE (dbcsr_type_complex_8)
       c_dp => dbcsr_get_data_p_z (area, bp, bp+nze-1)
       CALL pointer_z_rank_remap2(area%d%c2_dp, rowsize, csize,&
            c_dp)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
  END SUBROUTINE dbcsr_data_set_2d_pointer

! *****************************************************************************
!> \brief Nullifies a 2-d pointer to current 1-d data
!>
!> \param[in,out] area         data area to work on
! *****************************************************************************
  SUBROUTINE dbcsr_data_clear_2d_pointer (area)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_clear_2d_pointer', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: dt1
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------

    dt1 = area%d%data_type
    SELECT CASE (dt1)
    CASE (dbcsr_type_real_4)
       NULLIFY (area%d%r2_sp)
    CASE (dbcsr_type_real_8)
       NULLIFY (area%d%r2_dp)
    CASE (dbcsr_type_complex_4)
       NULLIFY (area%d%c2_sp)
    CASE (dbcsr_type_complex_8)
       NULLIFY (area%d%c2_dp)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT
  END SUBROUTINE dbcsr_data_clear_2d_pointer




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
    DEALLOCATE (marray%mats)
    NULLIFY (marray%mats)
  END SUBROUTINE dbcsr_destroy_2d_array

! *****************************************************************************
!> \brief Destroys a DBCSR distribution for a matrix multiplication based on
!>        the right matrix
!> \param[inout] imgdist_product        product distribution repetition
! *****************************************************************************
  SUBROUTINE dbcsr_destroy_image_dist(imgdist)
    TYPE(dbcsr_imagedistribution_type), &
      INTENT(INOUT)                          :: imgdist

!
!   ---------------------------------------------------------------------------

    CALL array_release (imgdist%row_image)
    CALL array_release (imgdist%col_image)
    CALL dbcsr_distribution_release (imgdist%main)
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
    CHARACTER, INTENT(IN)                    :: data_type

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

! *****************************************************************************
!> \brief Checks whether the data type is 2-D.
!> \result Data type is 2-D.
! *****************************************************************************
  PURE FUNCTION dbcsr_type_is_2d (data_type)
    CHARACTER, INTENT(IN)                    :: data_type
    LOGICAL                                  :: dbcsr_type_is_2d

    dbcsr_type_is_2d = data_type .EQ. dbcsr_type_real_4_2d .OR.&
         data_type .EQ. dbcsr_type_real_8_2d .OR.&
         data_type .EQ. dbcsr_type_complex_4_2d .OR.&
         data_type .EQ. dbcsr_type_complex_8_2d
  END FUNCTION dbcsr_type_is_2d

! *****************************************************************************
!> \brief Returns 1-d data type corresponding to the given 2-D one.
!> \result 1-D data type
! *****************************************************************************
  PURE FUNCTION dbcsr_type_2d_to_1d (data_type)
    CHARACTER, INTENT(IN)                    :: data_type
    CHARACTER                                :: dbcsr_type_2d_to_1d

    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4_2d)
       dbcsr_type_2d_to_1d = dbcsr_type_real_4
    CASE (dbcsr_type_real_8_2d)
       dbcsr_type_2d_to_1d = dbcsr_type_real_8
    CASE (dbcsr_type_complex_4_2d)
       dbcsr_type_2d_to_1d = dbcsr_type_complex_4
    CASE (dbcsr_type_complex_8_2d)
       dbcsr_type_2d_to_1d = dbcsr_type_complex_8
    CASE (dbcsr_type_real_4)
       dbcsr_type_2d_to_1d = dbcsr_type_real_4
    CASE (dbcsr_type_real_8)
       dbcsr_type_2d_to_1d = dbcsr_type_real_8
    CASE (dbcsr_type_complex_4)
       dbcsr_type_2d_to_1d = dbcsr_type_complex_4
    CASE (dbcsr_type_complex_8)
       dbcsr_type_2d_to_1d = dbcsr_type_complex_8
    CASE default
       dbcsr_type_2d_to_1d = '0'
    END SELECT
  END FUNCTION dbcsr_type_2d_to_1d

! *****************************************************************************
!> \brief Returns 2-D data type corresponding to the given 1-D one.
!> \result 2-D data type
! *****************************************************************************
  PURE FUNCTION dbcsr_type_1d_to_2d (data_type)
    CHARACTER, INTENT(IN)                    :: data_type
    CHARACTER                                :: dbcsr_type_1d_to_2d

    SELECT CASE (data_type)
    CASE (dbcsr_type_real_4)
       dbcsr_type_1d_to_2d = dbcsr_type_real_4_2d
    CASE (dbcsr_type_real_8)
       dbcsr_type_1d_to_2d = dbcsr_type_real_8_2d
    CASE (dbcsr_type_complex_4)
       dbcsr_type_1d_to_2d = dbcsr_type_complex_4_2d
    CASE (dbcsr_type_complex_8)
       dbcsr_type_1d_to_2d = dbcsr_type_complex_8_2d
    CASE (dbcsr_type_real_4_2d)
       dbcsr_type_1d_to_2d = dbcsr_type_real_4_2d
    CASE (dbcsr_type_real_8_2d)
       dbcsr_type_1d_to_2d = dbcsr_type_real_8_2d
    CASE (dbcsr_type_complex_4_2d)
       dbcsr_type_1d_to_2d = dbcsr_type_complex_4_2d
    CASE (dbcsr_type_complex_8_2d)
       dbcsr_type_1d_to_2d = dbcsr_type_complex_8_2d
    CASE default
       dbcsr_type_1d_to_2d = '0'
    END SELECT
  END FUNCTION dbcsr_type_1d_to_2d


#include "dbcsr_methods_d_.F"
#include "dbcsr_methods_z_.F"
#include "dbcsr_methods_s_.F"
#include "dbcsr_methods_c_.F"

END MODULE dbcsr_methods

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Higher-level operations on DBCSR matrices.
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.9
!>
!> <b>Modification history:</b>
!  - Created 2009-05-12
! *****************************************************************************

MODULE dbcsr_operations
  USE array_types,                     ONLY: array_data,&
                                             array_get,&
                                             array_i1d_obj
  USE dbcsr_block_access,              ONLY: dbcsr_get_block_p,&
                                             dbcsr_put_block,&
                                             dbcsr_remove_block,&
                                             dbcsr_reserve_blocks
  USE dbcsr_block_operations,          ONLY: block_add_on_diag,&
                                             dbcsr_block_conjg,&
                                             dbcsr_block_partial_copy,&
                                             dbcsr_block_real_neg,&
                                             dbcsr_block_scale,&
                                             dbcsr_data_clear,&
                                             get_block2d_diagonal,&
                                             set_block2d_diagonal
  USE dbcsr_data_methods,              ONLY: &
       dbcsr_data_clear_pointer, dbcsr_data_ensure_size, dbcsr_data_get_size, &
       dbcsr_data_get_size_referenced, dbcsr_data_get_type, dbcsr_data_init, &
       dbcsr_data_new, dbcsr_data_release, dbcsr_data_set_pointer, &
       dbcsr_get_data, dbcsr_scalar, dbcsr_scalar_are_equal, &
       dbcsr_scalar_fill_all, dbcsr_scalar_get_type, dbcsr_scalar_one, &
       dbcsr_scalar_zero, dbcsr_type_1d_to_2d
  USE dbcsr_data_operations,           ONLY: dbcsr_data_convert,&
                                             dbcsr_data_copyall,&
                                             dbcsr_switch_data_area
  USE dbcsr_dist_methods,              ONLY: dbcsr_distribution_col_dist,&
                                             dbcsr_distribution_local_cols,&
                                             dbcsr_distribution_local_rows,&
                                             dbcsr_distribution_mp,&
                                             dbcsr_distribution_row_dist
  USE dbcsr_dist_operations,           ONLY: checker_square_proc,&
                                             checker_tr,&
                                             dbcsr_find_column,&
                                             dbcsr_get_stored_coordinates
  USE dbcsr_error_handling,            ONLY: &
       dbcsr_assert, dbcsr_caller_error, dbcsr_error_set, dbcsr_error_stop, &
       dbcsr_error_type, dbcsr_failure_level, dbcsr_fatal_level, &
       dbcsr_internal_error, dbcsr_unimplemented_error_nr, &
       dbcsr_warning_level, dbcsr_wrong_args_error
  USE dbcsr_index_operations,          ONLY: dbcsr_index_checksum,&
                                             dbcsr_repoint_index
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_offsets, dbcsr_distribution, dbcsr_get_data_size, &
       dbcsr_get_data_type, dbcsr_get_index_memory_type, &
       dbcsr_get_matrix_type, dbcsr_get_num_blocks, &
       dbcsr_get_replication_type, dbcsr_has_symmetry, dbcsr_is_initialized, &
       dbcsr_max_col_size, dbcsr_max_row_size, dbcsr_name, &
       dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_nfullcols_total, &
       dbcsr_nfullrows_total, dbcsr_row_block_offsets, dbcsr_valid_index
  USE dbcsr_mp_methods,                ONLY: dbcsr_mp_group,&
                                             dbcsr_mp_mynode,&
                                             dbcsr_mp_mypcol,&
                                             dbcsr_mp_myprow,&
                                             dbcsr_mp_numnodes,&
                                             dbcsr_mp_pgrid
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size,&
                                             pointer_view
  USE dbcsr_toollib,                   ONLY: swap
  USE dbcsr_types,                     ONLY: &
       dbcsr_data_obj, dbcsr_distribution_obj, dbcsr_filter_frobenius, &
       dbcsr_func_artanh, dbcsr_func_asin, dbcsr_func_cos, dbcsr_func_ddsin, &
       dbcsr_func_ddtanh, dbcsr_func_dsin, dbcsr_func_dtanh, &
       dbcsr_func_inverse, dbcsr_func_inverse_special, dbcsr_func_sin, &
       dbcsr_func_spread_from_zero, dbcsr_func_tanh, dbcsr_func_truncate, &
       dbcsr_iterator, dbcsr_mp_obj, dbcsr_norm_column, dbcsr_norm_frobenius, &
       dbcsr_norm_gershgorin, dbcsr_norm_maxabsnorm, dbcsr_obj, &
       dbcsr_repl_full, dbcsr_repl_none, dbcsr_scalar_type, &
       dbcsr_type_antihermitian, dbcsr_type_antisymmetric, &
       dbcsr_type_complex_4, dbcsr_type_complex_8, dbcsr_type_hermitian, &
       dbcsr_type_invalid, dbcsr_type_no_symmetry, dbcsr_type_real_4, &
       dbcsr_type_real_8, dbcsr_type_symmetric
  USE dbcsr_util,                      ONLY: find_block_of_element
  USE dbcsr_work_operations,           ONLY: dbcsr_create,&
                                             dbcsr_finalize,&
                                             dbcsr_work_create
  USE kinds,                           ONLY: dp,&
                                             int_8,&
                                             real_4,&
                                             real_8,&
                                             sp
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_max,&
                                             mp_sum

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_operations'

  ! prettify protection
  CHARACTER, PARAMETER :: xa=dbcsr_type_hermitian, xb=dbcsr_type_antihermitian,&
                          xc=dbcsr_type_no_symmetry

  PUBLIC :: dbcsr_trace, dbcsr_add_on_diag,&
            dbcsr_set, dbcsr_scale, dbcsr_add, dbcsr_copy,&
            dbcsr_copy_into_existing,&
            dbcsr_get_diag, dbcsr_set_diag,&
            dbcsr_get_block_diag, dbcsr_hadamard_product,&
            dbcsr_filter, dbcsr_scale_by_vector,&
            dbcsr_function_of_elements,&
            dbcsr_triu,&
            dbcsr_init_random
  PUBLIC :: dbcsr_sum_replicated
  PUBLIC :: dbcsr_norm, dbcsr_conjg,&
            dbcsr_gershgorin_norm, dbcsr_maxabs, dbcsr_frobenius_norm

  PUBLIC :: dbcsr_crop_matrix
  PUBLIC :: dbcsr_get_info, dbcsr_may_be_dense, dbcsr_get_occupation

! The interfaces for the generic routines found in the generated
! generic files.


  INTERFACE dbcsr_trace
     MODULE PROCEDURE dbcsr_trace_a_any
     MODULE PROCEDURE dbcsr_trace_a_s, dbcsr_trace_a_d,&
                      dbcsr_trace_a_c, dbcsr_trace_a_z
     MODULE PROCEDURE dbcsr_trace_ab_s, dbcsr_trace_a_b_d,&
                      dbcsr_trace_ab_c, dbcsr_trace_ab_z
  END INTERFACE

  INTERFACE dbcsr_scale
     MODULE PROCEDURE dbcsr_scale_anytype
     MODULE PROCEDURE dbcsr_scale_s, dbcsr_scale_d,&
                      dbcsr_scale_c, dbcsr_scale_z
  END INTERFACE

  INTERFACE dbcsr_scale_by_vector
     MODULE PROCEDURE dbcsr_scale_by_vector_anytype
     MODULE PROCEDURE dbcsr_scale_by_vector_s, dbcsr_scale_by_vector_d,&
                      dbcsr_scale_by_vector_c, dbcsr_scale_by_vector_z
  END INTERFACE

  INTERFACE dbcsr_set
     MODULE PROCEDURE dbcsr_set_anytype
     MODULE PROCEDURE dbcsr_set_s, dbcsr_set_d, dbcsr_set_c, dbcsr_set_z
  END INTERFACE

  INTERFACE dbcsr_add
     MODULE PROCEDURE dbcsr_add_anytype
     MODULE PROCEDURE dbcsr_add_s, dbcsr_add_d,&
                      dbcsr_add_c, dbcsr_add_z
  END INTERFACE

  INTERFACE dbcsr_add_on_diag
     MODULE PROCEDURE dbcsr_add_on_diag_anytype
     MODULE PROCEDURE dbcsr_add_on_diag_s, dbcsr_add_on_diag_d,&
                      dbcsr_add_on_diag_c, dbcsr_add_on_diag_z
  END INTERFACE

  INTERFACE dbcsr_filter
     MODULE PROCEDURE dbcsr_filter_anytype
     MODULE PROCEDURE dbcsr_filter_s, dbcsr_filter_d,&
                      dbcsr_filter_c, dbcsr_filter_z
  END INTERFACE

  INTERFACE dbcsr_get_diag
     MODULE PROCEDURE dbcsr_get_diag_anytype
     MODULE PROCEDURE dbcsr_get_diag_s, dbcsr_get_diag_d,&
                      dbcsr_get_diag_c, dbcsr_get_diag_z
  END INTERFACE

  INTERFACE dbcsr_set_diag
     MODULE PROCEDURE dbcsr_set_diag_anytype
     MODULE PROCEDURE dbcsr_set_diag_s, dbcsr_set_diag_d,&
                      dbcsr_set_diag_c, dbcsr_set_diag_z
  END INTERFACE

  INTERFACE dbcsr_norm
     MODULE PROCEDURE dbcsr_norm_anytype
     MODULE PROCEDURE dbcsr_norm_r4_scal
     MODULE PROCEDURE dbcsr_norm_r4_vec, dbcsr_norm_r8_vec
  END INTERFACE
  INTERFACE dbcsr_gershgorin_norm
     MODULE PROCEDURE dbcsr_gershgorin_norm_r8
  END INTERFACE
  INTERFACE dbcsr_maxabs
     MODULE PROCEDURE dbcsr_maxabs_r8
  END INTERFACE
  INTERFACE dbcsr_frobenius_norm
     MODULE PROCEDURE dbcsr_frobenius_norm_r8
  END INTERFACE

  LOGICAL, PARAMETER :: debug_mod = .FALSE.
  LOGICAL, PARAMETER :: careful_mod = .FALSE.

#define temp_transpose(v, r, c) RESHAPE(TRANSPOSE(RESHAPE(v,(/r,c/))),(/r*c/))

  INTEGER, PARAMETER, PRIVATE :: rpslot_owner = 1
  INTEGER, PARAMETER, PRIVATE :: rpslot_addblks = 2
  INTEGER, PARAMETER, PRIVATE :: rpslot_addoffset = 3
  INTEGER, PARAMETER, PRIVATE :: rpslot_oldblks = 4
  INTEGER, PARAMETER, PRIVATE :: rpslot_oldoffset = 5
  INTEGER, PARAMETER, PRIVATE :: rpslot_totaloffset = 6
  INTEGER, PARAMETER, PRIVATE :: rpnslots = 6


CONTAINS


! *****************************************************************************
!> \brief Conjugate a DBCSR matrix
!> \param[inout] matrix       DBCSR matrix
!>
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_conjg(matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_conjg', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, data_type, &
                                                error_handler, row
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_any
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------
!

    CALL dbcsr_error_set(routineN, error_handler, error)
    data_type = dbcsr_get_data_type(matrix)
    CALL dbcsr_data_init (data_any)
    CALL dbcsr_data_new (data_any, data_type)
    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_any, tr, blk)
       SELECT CASE (data_type)
       CASE (dbcsr_type_complex_4)
          data_any%d%c_sp = CONJG(data_any%d%c_sp)
       CASE (dbcsr_type_complex_8)
          data_any%d%c_dp = CONJG(data_any%d%c_dp)
       CASE DEFAULT
          ! needed for g95
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_data_clear_pointer(data_any)
    CALL dbcsr_data_release(data_any)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_conjg


! *****************************************************************************
!> \brief Scales a DBCSR matrix by alpha
!> \param[inout] matrix_a       DBCSR matrix
!> \param[in] alpha_scalar      (optional) a scalar
!> \param[in] limits            (optional) Scale only a subbox
!> \param[in,out] error         error
!> \par Limits
!> A 4-tuple describing (first_row, last_row, first_column, last_column). Set
!> to 0 to avoid limiting.
! *****************************************************************************
  SUBROUTINE dbcsr_scale_anytype(matrix_a, alpha_scalar, limits, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: alpha_scalar
    INTEGER, DIMENSION(4), INTENT(IN), &
      OPTIONAL                               :: limits
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_anytype', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: first_col_i = 3, &
                                                first_row_i = 1, &
                                                last_col_i = 4, last_row_i = 2

    INTEGER :: a_col, a_col_size, a_row, a_row_size, col_offset, &
      error_handler, row_offset, scale_col_offset, scale_col_size, &
      scale_row_offset, scale_row_size
    INTEGER, DIMENSION(4)                    :: my_limits
    LOGICAL                                  :: do_scale, has_limits, tr
    TYPE(dbcsr_data_obj)                     :: data_any
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_scalar_type)                  :: one

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    ! Limits are only honored if the argument is present and any are
    ! non-zero.
    IF (PRESENT (limits)) THEN
       has_limits = ANY (limits(:) .NE. 0)
    ELSE
       has_limits = .FALSE.
    ENDIF
    my_limits(first_row_i) = 1
    my_limits(last_row_i)  = dbcsr_nfullrows_total (matrix_a)
    my_limits(first_col_i) = 1
    my_limits(last_col_i)  = dbcsr_nfullcols_total (matrix_a)
    IF (has_limits) THEN
       IF (limits(last_col_i) .NE. 0) THEN
          IF (debug_mod) THEN
             CALL dbcsr_assert (limits(last_col_i) .GE. 0, "AND",&
                  limits(last_col_i) .LE. dbcsr_nfullcols_total (matrix_a),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
               "Specified last column is out of bounds.",__LINE__,error)
          ENDIF
          my_limits(last_col_i) = limits(last_col_i)
       ENDIF
       IF (limits(first_col_i) .NE. 0) THEN
          IF (debug_mod) THEN
             CALL dbcsr_assert (limits(first_col_i) .GE. 0, "AND",&
                  limits(first_col_i) .LE. dbcsr_nfullcols_total (matrix_a),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
               "Specified first column is out of bounds.",__LINE__,error)
          ENDIF
          my_limits(first_col_i) = limits(first_col_i)
       ENDIF
       IF (limits(last_row_i) .NE. 0) THEN
          IF (debug_mod) THEN
             CALL dbcsr_assert (limits(last_row_i) .GE. 0, "AND",&
                  limits(last_row_i) .LE. dbcsr_nfullrows_total (matrix_a),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
               "Specified last row is out of bounds.",__LINE__,error)
          ENDIF
          my_limits(last_row_i) = limits(last_row_i)
       ENDIF
       IF (limits(first_row_i) .NE. 0) THEN
          IF (debug_mod) THEN
             CALL dbcsr_assert (limits(first_row_i) .GE. 0, "AND",&
                  limits(first_row_i) .LE. dbcsr_nfullrows_total (matrix_a),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
               "Specified first row is out of bounds.",__LINE__,error)
          ENDIF
          my_limits(first_row_i) = limits(first_row_i)
       ENDIF
    ENDIF
    !
    ! quick return if possible
    one = dbcsr_scalar_one (dbcsr_scalar_get_type (alpha_scalar))
    do_scale = .NOT. dbcsr_scalar_are_equal (alpha_scalar, one)
    !
    ! let's go
    IF(do_scale) THEN
       !$OMP PARALLEL DEFAULT (none) &
       !$OMP          PRIVATE (iter, data_any) &
       !$OMP          PRIVATE (a_row, a_col, tr, a_row_size, a_col_size, &
       !$OMP                   row_offset, col_offset) &
       !$OMP          PRIVATE (scale_row_size, scale_col_size,&
       !$OMP                   scale_row_offset, scale_col_offset) &
       !$OMP          SHARED (matrix_a, my_limits, error, alpha_scalar)
       CALL dbcsr_data_init (data_any)
       CALL dbcsr_data_new (data_any, dbcsr_type_1d_to_2d (dbcsr_get_data_type (matrix_a)))
       CALL dbcsr_iterator_start(iter, matrix_a, read_only=.FALSE.,&
            contiguous_pointers = .FALSE., dynamic = .TRUE.,&
            dynamic_byrows = .TRUE., shared=.TRUE.)
       iterations: DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block(iter, a_row, a_col, data_any, tr,&
               row_size=a_row_size, col_size=a_col_size, &
               row_offset=row_offset, col_offset=col_offset)
          IF (a_row_size .GT. 0 .AND. a_col_size .GT. 0) THEN
             CALL frame_block_limit (a_row_size, row_offset,&
                  my_limits(first_row_i), my_limits(last_row_i),&
                  scale_row_size, scale_row_offset)
             CALL frame_block_limit (a_col_size, col_offset,&
                  my_limits(first_col_i), my_limits(last_col_i),&
                  scale_col_size, scale_col_offset)
             IF (tr) THEN
                CALL swap (scale_row_size, scale_col_size)
                CALL swap (scale_row_offset, scale_col_offset)
             ENDIF
             CALL dbcsr_block_scale (data_any, scale=alpha_scalar,&
                  row_size=scale_row_size, col_size=scale_col_size,&
                  lb=scale_row_offset, lb2=scale_col_offset, error=error)
          ENDIF
       ENDDO iterations
       CALL dbcsr_iterator_stop(iter)
       CALL dbcsr_data_clear_pointer (data_any)
       CALL dbcsr_data_release (data_any)
       !$OMP END PARALLEL
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_scale_anytype


! *****************************************************************************
!> \brief Determines the effect of limits on a block
!> \param[in] block_size     size of block
!> \param[in] block_offset   global offset of block
!> \param[in] first_limit    lower limit
!> \param[in] last_limit     upper limit
!> \param[out] frame_size    size of block region within the limits
!> \param[out] frame_offset  starting position of the block region
!>                           that is within the limits
! *****************************************************************************
  ELEMENTAL SUBROUTINE frame_block_limit (block_size, block_offset,&
       first_limit, last_limit,&
       frame_size, frame_offset)
    INTEGER, INTENT(IN)                      :: block_size, block_offset, &
                                                first_limit, last_limit
    INTEGER, INTENT(OUT)                     :: frame_size, frame_offset

    INTEGER                                  :: f, l

    f = MAX (block_offset, first_limit)
    l = MIN (block_offset + block_size - 1, last_limit)
    frame_size = MAX(l - f + 1, 0)
    frame_offset = MIN(f - block_offset + 1, block_size)
  END SUBROUTINE frame_block_limit


! *****************************************************************************
!> \brief Scales a DBCSR matrix by alpha
!> \param[inout] matrix_a       DBCSR matrix
!> \param[in] alpha             the scaling vector
!> \param[in] side              apply the scaling from the side
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_scale_by_vector_anytype(matrix_a, alpha, side, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_data_obj), INTENT(IN), &
      OPTIONAL                               :: alpha
    CHARACTER(LEN=*), INTENT(IN)             :: side
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_scale_by_vector_anytype', &
      routineP = moduleN//':'//routineN

    INTEGER :: a_blk, a_col, a_col_size, a_nze, a_row, a_row_size, &
      col_offset, data_type, error_handler, i, icol, irow, row_offset
    LOGICAL                                  :: right, tr
    TYPE(dbcsr_data_obj)                     :: data_any
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

! check that alpha and matrix have the same data type
    CALL dbcsr_assert (dbcsr_get_data_type (matrix_a).EQ.alpha%d%data_type, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "wrong data type matrix_a",__LINE__,error)
    !
    ! set vars
    right = .TRUE.
    SELECT CASE(side)
    CASE('right');right = .TRUE.
    CASE('left' );right = .FALSE.
    CASE DEFAULT
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, &
            routineN, "wrong side="//side,__LINE__,error)
    END SELECT
    !
    ! let's go
    data_type = dbcsr_get_data_type (matrix_a)
    CALL dbcsr_data_init(data_any)
    CALL dbcsr_data_new(data_any, dbcsr_get_data_type (matrix_a))
    CALL dbcsr_iterator_start(iter, matrix_a)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block(iter, a_row, a_col, data_any, tr,&
            block_number=a_blk,&
            row_size=a_row_size, col_size=a_col_size, &
            row_offset=row_offset, col_offset=col_offset)
       a_nze = a_row_size * a_col_size
       IF (a_nze .EQ. 0) CYCLE ! Skip empty blocks
       !
       ! let's scale
       IF(right) THEN
          SELECT CASE (data_type)
          CASE (dbcsr_type_real_4)
             DO i = 1,a_col_size
                DO icol=(i-1)*a_row_size+1,(i-1)*a_row_size+a_row_size
                   data_any%d%r_sp(icol)=data_any%d%r_sp(icol)*alpha%d%r_sp(col_offset+i-1)
                END DO
             END DO
          CASE (dbcsr_type_real_8)
             DO i = 1,a_col_size
                DO icol=(i-1)*a_row_size+1,(i-1)*a_row_size+a_row_size
                   data_any%d%r_dp(icol)=data_any%d%r_dp(icol)*alpha%d%r_dp(col_offset+i-1)
                END DO
             END DO
          CASE (dbcsr_type_complex_4)
             DO i = 1,a_col_size
                DO icol=(i-1)*a_row_size+1,(i-1)*a_row_size+a_row_size
                   data_any%d%c_sp(icol)=data_any%d%c_sp(icol)*alpha%d%c_sp(col_offset+i-1)
                END DO
             END DO
          CASE (dbcsr_type_complex_8)
             DO i = 1,a_col_size
                DO icol=(i-1)*a_row_size+1,(i-1)*a_row_size+a_row_size
                   data_any%d%c_dp(icol)=data_any%d%c_dp(icol)*alpha%d%c_dp(col_offset+i-1)
                END DO
             END DO
          END SELECT
       ELSE
          SELECT CASE (data_type)
          CASE (dbcsr_type_real_4)
             DO i = 1,a_row_size
                DO irow=i,i+a_col_size**2-1,a_col_size
                   data_any%d%r_sp(irow)=data_any%d%r_sp(irow)*alpha%d%r_sp(row_offset+i-1)
                END DO
             ENDDO
          CASE (dbcsr_type_real_8)
             DO i = 1,a_row_size
                DO irow=i,i+a_col_size**2-1,a_col_size
                   data_any%d%r_dp(irow)=data_any%d%r_dp(irow)*alpha%d%r_dp(row_offset+i-1)
                END DO
             ENDDO
          CASE (dbcsr_type_complex_4)
             DO i = 1,a_row_size
                DO irow=i,i+a_col_size**2-1,a_col_size
                   data_any%d%c_sp(irow)=data_any%d%c_sp(irow)*alpha%d%c_sp(row_offset+i-1)
                END DO
             ENDDO
          CASE (dbcsr_type_complex_8)
             DO i = 1,a_row_size
                DO irow=i,i+a_col_size**2-1,a_col_size
                   data_any%d%c_dp(irow)=data_any%d%c_dp(irow)*alpha%d%c_dp(row_offset+i-1)
                END DO
             ENDDO
          END SELECT
       ENDIF
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_data_clear_pointer (data_any)
    CALL dbcsr_data_release (data_any)
    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_scale_by_vector_anytype


! *****************************************************************************
!> \brief Set a DBCSR matrix by alpha
!> \param[inout] matrix       DBCSR matrix
!> \param[in] alpha           a scalar
!>
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_set_anytype(matrix, alpha, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: alpha
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_anytype', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, error_handler, row
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_block
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------
!

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_data_init (data_block)
    CALL dbcsr_data_new (data_block, dbcsr_get_data_type (matrix))
    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_block, tr, blk)
       CALL dbcsr_data_clear (data_block, value=alpha)
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_data_clear_pointer (data_block)
    CALL dbcsr_data_release (data_block)
    CALL dbcsr_error_stop(error_handler, error)
    !
  END SUBROUTINE dbcsr_set_anytype


! *****************************************************************************
!> \brief add and scale matrices
!>    A = alpha*A + beta*B or
!> \param[in,out] matrix_a   DBCSR matrix
!> \param[in] matrix_b       DBCSR matrix
!> \param[in] alpha_scalar   (optional)
!> \param[in] beta_scalar    (optional)
!>
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_add_anytype(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: alpha_scalar, beta_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_anytype', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, data_type_a, &
                                                data_type_b, error_handler, &
                                                row, size_a, size_b
    LOGICAL                                  :: do_scale, tr
    TYPE(dbcsr_data_obj)                     :: data_block
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_scalar_type)                  :: my_alpha_scalar, &
                                                my_beta_scalar

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_valid_index (matrix_a), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Invalid matrix", __LINE__, error=error)

    IF((dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric.OR.&
        dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric).NEQV. &
       (dbcsr_get_matrix_type(matrix_b).EQ.dbcsr_type_symmetric.OR.&
        dbcsr_get_matrix_type(matrix_b).EQ.dbcsr_type_antisymmetric)) THEN
       CALL dbcsr_assert (.FALSE.,dbcsr_fatal_level,&
            dbcsr_unimplemented_error_nr, routineN, "Summing general with symmetric matrix NYI",&
            __LINE__,error)
    ENDIF

    data_type_a = dbcsr_get_data_type(matrix_a)
    data_type_b = dbcsr_get_data_type(matrix_b)
    !
    my_alpha_scalar = dbcsr_scalar_one (data_type_a)
    IF(PRESENT(alpha_scalar)) my_alpha_scalar = alpha_scalar
    my_beta_scalar = dbcsr_scalar_one (data_type_b)
    IF(PRESENT(beta_scalar)) my_beta_scalar = beta_scalar
    !
    ! let's go
    CALL dbcsr_assert (dbcsr_nblkrows_total(matrix_a).EQ.dbcsr_nblkrows_total(matrix_b), &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, "matrices not consistent",__LINE__,error)

    do_scale = dbcsr_scalar_are_equal (&
         my_beta_scalar, dbcsr_scalar_one (data_type_b))

    CALL dbcsr_scale(matrix_a, alpha_scalar=my_alpha_scalar, error=error)

    ! Pre-size work arrays of matrix_a to avoid continuous reallocation.
    size_a = dbcsr_data_get_size_referenced (matrix_a%m%data_area)
    size_b = dbcsr_data_get_size_referenced (matrix_b%m%data_area)
    IF(.NOT.dbcsr_scalar_are_equal(my_beta_scalar,&
         dbcsr_scalar_zero(data_type_b))) THEN
       !$OMP PARALLEL DEFAULT (none) &
       !$OMP          PRIVATE (iter, data_block) &
       !$OMP          PRIVATE (row, col, tr, blk) &
       !$OMP          SHARED (matrix_a, matrix_b, data_type_b, size_b, size_a) &
       !$OMP          SHARED (do_scale, my_beta_scalar) &
       !$OMP          SHARED (error)
       IF (size_b .GT. size_a .AND. matrix_b%m%nblks .GT. matrix_a%m%nblks) THEN
          CALL dbcsr_work_create (matrix_a,&
               nblks_guess = matrix_b%m%nblks - matrix_a%m%nblks,&
               sizedata_guess = size_b - size_a,&
               work_mutable = .FALSE., error=error)
       ELSE
          CALL dbcsr_work_create (matrix_a,&
               work_mutable = .FALSE., error=error)
       ENDIF
       !$OMP BARRIER
       CALL dbcsr_data_init (data_block)
       CALL dbcsr_data_new (data_block, data_type_b)
       CALL dbcsr_iterator_start(iter, matrix_b,&
            shared = .TRUE., read_only = .TRUE., contiguous_pointers=.FALSE.,&
            dynamic = .TRUE., dynamic_byrows = .TRUE.)

       DO WHILE (dbcsr_iterator_blocks_left(iter))

          CALL dbcsr_iterator_next_block(iter, row, col, data_block, tr, blk)

             IF (do_scale) THEN
                CALL dbcsr_put_block(matrix_a, row, col, data_block, tr,&
                     summation=.TRUE.)
             ELSE
                CALL dbcsr_put_block(matrix_a, row, col, data_block, tr,&
                     summation=.TRUE., scale=my_beta_scalar)
             ENDIF

       ENDDO

       CALL dbcsr_iterator_stop(iter)
       CALL dbcsr_finalize (matrix_a, error=error)
       CALL dbcsr_data_clear_pointer (data_block)
       CALL dbcsr_data_release (data_block)
       !$OMP END PARALLEL
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_add_anytype

!> \brief Interface for dbcsr_add
!> \param matrix_a ...
!> \param matrix_b ...
!> \param alpha_scalar ...
!> \param beta_scalar ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_add_d(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    REAL(real_8), INTENT(IN)                 :: alpha_scalar, beta_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
    IF(    dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_8 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_8) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar),&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(REAL(alpha_scalar,real_4)),&
            beta_scalar=dbcsr_scalar(REAL(beta_scalar,real_4)), error=error)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_4 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_complex_4) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(CMPLX(alpha_scalar,0,real_4)),&
            beta_scalar=dbcsr_scalar(CMPLX(beta_scalar,0,real_4)), error=error)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_8 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_complex_8) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(CMPLX(alpha_scalar,0,real_8)),&
            beta_scalar=dbcsr_scalar(CMPLX(beta_scalar,0,real_8)), error=error)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_8 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar),&
            beta_scalar=dbcsr_scalar(REAL(beta_scalar,real_4)), error=error)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_8) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(REAL(alpha_scalar,real_4)),&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_add_d

! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param alpha_scalar ...
!> \param beta_scalar ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_add_s(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    REAL(real_4), INTENT(IN)                 :: alpha_scalar, beta_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4 .AND. &
       dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar),&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_add_s

! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param alpha_scalar ...
!> \param beta_scalar ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_add_z(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    COMPLEX(real_8), INTENT(IN)              :: alpha_scalar, beta_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_8 .AND. &
       dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_complex_8) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar),&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_4 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_complex_4) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(CMPLX(alpha_scalar,KIND=real_4)),&
            beta_scalar=dbcsr_scalar(CMPLX(beta_scalar,KIND=real_4)), error=error)
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_add_z

! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param alpha_scalar ...
!> \param beta_scalar ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_add_c(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    COMPLEX(real_4), INTENT(IN)              :: alpha_scalar, beta_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_4 .AND. &
       dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_complex_4) THEN
       CALL dbcsr_add_anytype(matrix_a, matrix_b,&
            alpha_scalar=dbcsr_scalar(alpha_scalar),&
            beta_scalar=dbcsr_scalar(beta_scalar), error=error)
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_add_c

! *****************************************************************************
!> \brief Computes varios functions (defined by func) of matrix elements
!>        NOTE: sign(A,B) returns the value of A with the sign of B
!>        dbcsr_func_inverse:   1/(a1*x+a0)
!>           fails if the inversion produces infinite numbers 
!>        dbcsr_func_inverse_special: 1/(x+sign(a0,x))
!>           safe inverse: if a0>0 then the denominator is never zero
!>        dbcsr_func_tanh:    tanh(a1*x+a0)
!>        dbcsr_func_dtanh:   d(tanh(a1*x+a0)) / dx
!>        dbcsr_func_ddtanh:  d2(tanh(a1*x+a0)) / dx2
!>        dbcsr_func_artanh:  artanh(a1*x+a0)=ln[(1+(a1*x+a0))/(1-(a1*x+a0))]/2
!>           fails if |a1*x+a0| >= 1
!>        dbcsr_func_sread_from_zero:  if |x|<|a0| then x=sign(a0,x) 
!>        dbcsr_func_truncate:  if |x|>|a0| then x=sign(a0,x)
!>        dbcsr_func_sin:     sin(a1*x+a0)
!>        dbcsr_func_cos:     cos(a1*x+a0)
!>        dbcsr_func_dsin:    d(sin(a1*x+a0)) / dx = a1*cos(a1*x+a0)
!>        dbcsr_func_ddsin:   d2(sin(a1*x+a0)) / dx2 = -a1*a1*sin(a1*x+a0)
!>        dbcsr_func_asin:    asin(a1*x+a0)
!>           fails if |a1*x+a0| > 1
!> \param[in] matrix_a        DBCSR matrix
!> \param func ...
!> \param a0 ...
!> \param a1 ...
!> \param a2 ...
!> \param a3 ...
!> \param a4 ...
!> \param a5 ...
!> \param a6 ...
!> \param a7 ...
!> \param a8 ...
!> \param a9 ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_function_of_elements(matrix_a, func, a0, a1, a2, a3, a4,&
    a5, a6, a7, a8, a9, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    INTEGER, INTENT(IN)                      :: func
    REAL(kind=dp), INTENT(IN), OPTIONAL      :: a0, a1, a2, a3, a4, a5, a6, &
                                                a7, a8, a9
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_function_of_elements', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_size, &
                                                data_type, error_handler, &
                                                ielem, nze, row, row_size
    LOGICAL                                  :: tr_a
    REAL(kind=dp)                            :: p0, p1, p2
    TYPE(dbcsr_data_obj)                     :: a_data
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    
    IF (PRESENT(a0)) THEN
       p0=a0
    ELSE
       p0=0.0_dp
    ENDIF
    IF (PRESENT(a1)) THEN
       p1=a1
    ELSE
       p1=1.0_dp
    ENDIF
    IF (PRESENT(a2)) THEN
       p2=a2
    ELSE
       p2=0.0_dp
    ENDIF

    data_type = dbcsr_get_data_type(matrix_a)
    CALL dbcsr_data_init (a_data)
    CALL dbcsr_data_new (a_data, data_type)
    CALL dbcsr_iterator_start(iter, matrix_a)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, a_data, tr_a, blk, &
            row_size=row_size, col_size=col_size)
       nze = row_size * col_size
       SELECT CASE (data_type)
       !CASE (dbcsr_type_real_4)
       !   a_data%d%r_sp(1:nze) = 1.0_real_4/a_data%d%r_sp(1:nze)
       !   CALL dbcsr_assert(MAXVAL(ABS(a_data%d%r_sp)).LT.HUGE(0.0_real_4),&
       !           dbcsr_fatal_level,dbcsr_internal_error, routineN,&
       !           "Division by zero", __LINE__, error=error)
       CASE (dbcsr_type_real_8)
          SELECT CASE (func)
          CASE (dbcsr_func_spread_from_zero)
             ! if |x|<|a0| then x=|a0|*sign(x)
             DO ielem=1,nze
               IF (ABS(a_data%d%r_dp(ielem)).lt.ABS(p0)) THEN
                  a_data%d%r_dp(ielem) = SIGN(p0,a_data%d%r_dp(ielem))
               ENDIF
             ENDDO
          CASE (dbcsr_func_truncate)
             ! if |x|>|a0| then x=|a0|*sign(x)
             DO ielem=1,nze
               IF (ABS(a_data%d%r_dp(ielem)).gt.ABS(p0)) THEN
                  a_data%d%r_dp(ielem) = SIGN(p0,a_data%d%r_dp(ielem))
               ENDIF
             ENDDO
          CASE (dbcsr_func_inverse_special)
             !IF (MINVAL(ABS(a_data%d%r_dp)).le.ABS(p2)) THEN
             !   ! there is at least one near-zero element,
             !   ! invert element-by-element
             !   DO ielem=1,nze
             !     IF (a_data%d%r_dp(ielem).le.ABS(p2)) THEN
             !        a_data%d%r_dp(ielem) = 0.0_real_8
             !     ELSE
             !        a_data%d%r_dp(ielem) = &
             !           1.0_real_8/(p1*a_data%d%r_dp(ielem)+p0)
             !     ENDIF
             !   ENDDO
             !ELSE
             !   a_data%d%r_dp(1:nze) = 1.0_real_8/(p1*a_data%d%r_dp(1:nze)+p0)
             !ENDIF
             a_data%d%r_dp(1:nze) = 1.0_real_8/(a_data%d%r_dp(1:nze)+SIGN(p0,a_data%d%r_dp(1:nze)))
          CASE (dbcsr_func_inverse)
             a_data%d%r_dp(1:nze) = 1.0_real_8/(p1*a_data%d%r_dp(1:nze)+p0)
             CALL dbcsr_assert(MAXVAL(ABS(a_data%d%r_dp)).LT.HUGE(0.0_real_8),&
                     dbcsr_fatal_level,dbcsr_internal_error, routineN,&
                     "Division by zero", __LINE__, error=error)
          CASE (dbcsr_func_tanh)
             a_data%d%r_dp(1:nze) = TANH(p1*a_data%d%r_dp(1:nze)+p0)
          CASE (dbcsr_func_dtanh)
             a_data%d%r_dp(1:nze) = TANH(p1*a_data%d%r_dp(1:nze)+p0)
             a_data%d%r_dp(1:nze) = a_data%d%r_dp(1:nze)**2
             a_data%d%r_dp(1:nze) = p1*(1.0_real_8-a_data%d%r_dp(1:nze))
          CASE (dbcsr_func_ddtanh)
             a_data%d%r_dp(1:nze) = TANH(p1*a_data%d%r_dp(1:nze)+p0)
             a_data%d%r_dp(1:nze) = a_data%d%r_dp(1:nze)**3-a_data%d%r_dp(1:nze)
             a_data%d%r_dp(1:nze) = 2.0_real_8*(p1**2)*a_data%d%r_dp(1:nze)
          CASE (dbcsr_func_artanh)
             a_data%d%r_dp(1:nze) = p1*a_data%d%r_dp(1:nze)+p0
             CALL dbcsr_assert(MAXVAL(ABS(a_data%d%r_dp)).LT.1.0_real_8,&
                     dbcsr_fatal_level,dbcsr_internal_error, routineN,&
                     "ARTANH is undefined for |x|>=1", __LINE__, error=error)
             a_data%d%r_dp(1:nze) = (1.0_real_8+a_data%d%r_dp(1:nze))&
                                  / (1.0_real_8-a_data%d%r_dp(1:nze))
             a_data%d%r_dp(1:nze) = 0.5_real_8*LOG(a_data%d%r_dp(1:nze))
          CASE (dbcsr_func_sin)
             a_data%d%r_dp(1:nze) = SIN(p1*a_data%d%r_dp(1:nze)+p0)
          CASE (dbcsr_func_cos)
             a_data%d%r_dp(1:nze) = COS(p1*a_data%d%r_dp(1:nze)+p0)
          CASE (dbcsr_func_dsin)
             a_data%d%r_dp(1:nze) = p1*COS(p1*a_data%d%r_dp(1:nze)+p0)
          CASE (dbcsr_func_ddsin)
             a_data%d%r_dp(1:nze) = -p1*p1*SIN(p1*a_data%d%r_dp(1:nze)+p0)
          CASE (dbcsr_func_asin)
             a_data%d%r_dp(1:nze) = p1*a_data%d%r_dp(1:nze)+p0
             CALL dbcsr_assert(MAXVAL(ABS(a_data%d%r_dp)).LE.1.0_real_8,&
                     dbcsr_fatal_level,dbcsr_internal_error, routineN,&
                     "ASIN is undefined for |x|>1", __LINE__, error=error)
             a_data%d%r_dp(1:nze) = ASIN(a_data%d%r_dp(1:nze))
          CASE DEFAULT
             CALL dbcsr_assert(.FALSE.,&
                     dbcsr_fatal_level,dbcsr_internal_error, routineN,&
                     "Unknown function of matrix elements", __LINE__,&
                     error=error)
          END SELECT
       !CASE (dbcsr_type_complex_4)
       !CASE (dbcsr_type_complex_8)
       CASE DEFAULT
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level,&
                  dbcsr_unimplemented_error_nr, routineN,&
                  "Operation is implemented only for dp real values",&
                  __LINE__,error)
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_data_clear_pointer (a_data)
    CALL dbcsr_data_release (a_data)
    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_function_of_elements

! *****************************************************************************
!> \brief Hadamard product
!>    C = A . B (C needs to be different from A and B)
!> \param[in] matrix_a DBCSR matrix
!> \param[in] matrix_b DBCSR matrix
!> \param[in,out] matrix_c        DBCSR matrix
!> \param b_assume_value ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_hadamard_product(matrix_a, matrix_b, matrix_c,&
    b_assume_value, error)

    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: b_assume_value
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_hadamard_product', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_size, &
                                                data_type, error_handler, &
                                                nze, row, row_size
    LOGICAL                                  :: assume_blocks_in_b, found, &
                                                tr_a, tr_b
    REAL(KIND=dp)                            :: assumed_b_value
    TYPE(dbcsr_data_obj)                     :: a_data, b_data, c_data
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    IF (PRESENT(b_assume_value)) THEN
       assume_blocks_in_b=.TRUE.
       assumed_b_value=b_assume_value
    ELSE
       assume_blocks_in_b=.FALSE.
       assumed_b_value=0.0_dp
    ENDIF

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_get_data_type(matrix_a).EQ.dbcsr_get_data_type(matrix_b).AND.&
         dbcsr_get_data_type(matrix_a).EQ.dbcsr_get_data_type(matrix_c), &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, &
         "data types not consistent, need to fix that",__LINE__,error)

    CALL dbcsr_assert (dbcsr_nblkrows_total(matrix_a).EQ.dbcsr_nblkrows_total(matrix_b).AND.&
         dbcsr_nblkrows_total(matrix_c).EQ.dbcsr_nblkrows_total(matrix_a), &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, &
         "matrices not consistent",__LINE__,error)

    data_type = dbcsr_get_data_type(matrix_a)
    CALL dbcsr_data_init (c_data)
    CALL dbcsr_data_new (c_data, data_type,&
         data_size=dbcsr_max_row_size(matrix_a)*dbcsr_max_col_size(matrix_a))
    CALL dbcsr_set(matrix_c, dbcsr_scalar_zero(data_type), error=error)
    CALL dbcsr_data_init (a_data)
    CALL dbcsr_data_new (a_data, data_type)
    CALL dbcsr_data_init (b_data)
    CALL dbcsr_data_new (b_data, data_type)
   CALL dbcsr_iterator_start(iter, matrix_a)
   DO WHILE (dbcsr_iterator_blocks_left(iter))
       SELECT CASE (dbcsr_get_data_type(matrix_a))
          !CASE (dbcsr_type_real_4)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, a_data, tr_a, blk, &
               row_size=row_size, col_size=col_size)
          nze = row_size * col_size
          CALL dbcsr_get_block_p(matrix_b, row, col, b_data, tr_b, found)
          CALL dbcsr_assert (tr_a.EQV.tr_b, dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, &
               "tr not consistent, need to fix that",__LINE__,error)
          IF(found) THEN
             SELECT CASE (data_type)
             CASE (dbcsr_type_real_4)
                c_data%d%r_sp(1:nze) = a_data%d%r_sp(1:nze) * b_data%d%r_sp(1:nze)
             CASE (dbcsr_type_real_8)
                c_data%d%r_dp(1:nze) = a_data%d%r_dp(1:nze) * b_data%d%r_dp(1:nze)
             CASE (dbcsr_type_complex_4)
                c_data%d%c_sp(1:nze) = a_data%d%c_sp(1:nze) * b_data%d%c_sp(1:nze)
             CASE (dbcsr_type_complex_8)
                c_data%d%c_dp(1:nze) = a_data%d%c_dp(1:nze) * b_data%d%c_dp(1:nze)
             END SELECT
             CALL dbcsr_put_block(matrix_c, row, col, c_data, tr_a, &
                  summation=.FALSE.)
          ELSE
             IF (assume_blocks_in_b) THEN ! this makes not too much sense, to delete ?
                SELECT CASE (data_type)
                CASE (dbcsr_type_real_4)
                   c_data%d%r_sp(1:nze) = a_data%d%r_sp(1:nze) * REAL(assumed_b_value,KIND=sp)
                CASE (dbcsr_type_real_8)
                   c_data%d%r_dp(1:nze) = a_data%d%r_dp(1:nze) * assumed_b_value 
                CASE (dbcsr_type_complex_4)
                   c_data%d%c_sp(1:nze) = a_data%d%c_sp(1:nze) * REAL(assumed_b_value,KIND=sp)
                CASE (dbcsr_type_complex_8)
                   c_data%d%c_dp(1:nze) = a_data%d%c_dp(1:nze) * assumed_b_value
                END SELECT
                CALL dbcsr_put_block(matrix_c, row, col, c_data, tr_a, &
                     summation=.FALSE.)
             ENDIF
          ENDIF
          !CASE (dbcsr_type_complex_4)
          !CASE (dbcsr_type_complex_8)
       CASE DEFAULT
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Only real double precision",__LINE__,error)
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_finalize (matrix_c, error=error)
    CALL dbcsr_data_clear_pointer (a_data)
    CALL dbcsr_data_clear_pointer (b_data)
    CALL dbcsr_data_release (c_data)
    CALL dbcsr_data_release (a_data)
    CALL dbcsr_data_release (b_data)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_hadamard_product

! *****************************************************************************
!> \brief add a constant to the diagonal of a matrix
!> \param[inout] matrix       DBCSR matrix
!> \param[in]    alpha_scalar scalar
!> \param first_row ...
!> \param last_row ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_add_on_diag_anytype(matrix, alpha_scalar, first_row, last_row, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: alpha_scalar
    INTEGER, INTENT(in), OPTIONAL            :: first_row, last_row
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_on_diag_anytype', &
      routineP = moduleN//':'//routineN

    INTEGER :: error_handler, hold, imax, imin, my_first_row, my_last_row, &
      mynode, offset_beg, offset_end, row, row_size, stored_row
    INTEGER, DIMENSION(:), POINTER           :: row_blk_offsets
    LOGICAL                                  :: found, tr
    TYPE(dbcsr_data_obj)                     :: buff, data_a, small_buff

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_nblkrows_total(matrix).EQ.dbcsr_nblkcols_total(matrix).AND.&
         dbcsr_nfullrows_total(matrix).EQ.dbcsr_nfullrows_total(matrix), &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, "matrices not consistent",__LINE__,error)

    my_first_row = 1
    my_last_row = dbcsr_nfullrows_total(matrix)
    IF(PRESENT(first_row)) my_first_row = first_row
    IF(PRESENT(last_row)) my_last_row = last_row

    mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(dbcsr_distribution(matrix)))
    row_blk_offsets => array_data(dbcsr_col_block_offsets(matrix))

    CALL dbcsr_work_create(matrix, work_mutable=.TRUE., error=error)
    CALL dbcsr_data_init (buff)
    CALL dbcsr_data_init (data_a)
    CALL dbcsr_data_init (small_buff)
    CALL dbcsr_data_new (data_a,&
         dbcsr_type_1d_to_2d(dbcsr_get_data_type(matrix)))
    CALL dbcsr_data_new (buff,&
         dbcsr_type_1d_to_2d(dbcsr_get_data_type(matrix)),&
         dbcsr_max_row_size(matrix), dbcsr_max_col_size(matrix))
    CALL dbcsr_data_new (small_buff,&
         dbcsr_type_1d_to_2d(dbcsr_get_data_type(matrix)))

    DO row = 1,dbcsr_nblkrows_total(matrix)
       tr = .FALSE.
       stored_row = row
       CALL dbcsr_get_stored_coordinates (matrix, stored_row, stored_row, tr,&
            hold)
       IF(hold.EQ.mynode) THEN
          CALL dbcsr_get_block_p(matrix, stored_row, stored_row, data_a, tr,&
               found, row_size=row_size)
          offset_beg = row_blk_offsets(row)
          offset_end = row_blk_offsets(row+1)-1
          IF(my_first_row.GT.offset_end.OR.my_last_row.LT.offset_beg) CYCLE
          imin = 1
          IF(my_first_row.gt.offset_beg) THEN
             imin = my_first_row - offset_beg + 1
          ENDIF
          imax = row_size
          IF(my_last_row.lt.offset_end) THEN
             imax = my_last_row - offset_end + row_size
          ENDIF
          IF(found) THEN
             CALL block_add_on_diag(data_a, alpha_scalar, row_size, &
                  imin=imin, imax=imax, error=error)
          ELSE
             CALL dbcsr_data_set_pointer (small_buff, row_size, row_size,&
                  buff)
             CALL dbcsr_data_clear (small_buff)
             CALL block_add_on_diag (small_buff, alpha_scalar, row_size,&
                  imin=imin, imax=imax, error=error)
             CALL dbcsr_put_block (matrix, stored_row, stored_row, small_buff)
          ENDIF
       ENDIF
    ENDDO

    CALL dbcsr_data_clear_pointer (data_a)
    CALL dbcsr_data_clear_pointer (small_buff)
    CALL dbcsr_data_release (small_buff)
    CALL dbcsr_data_release (buff)
    CALL dbcsr_data_release (data_a)
    CALL dbcsr_finalize(matrix, error=error)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_add_on_diag_anytype

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param keep_sparsity ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_init_random(matrix, keep_sparsity, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    LOGICAL, OPTIONAL                        :: keep_sparsity
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_random', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, col_size, error_handler, &
                                                hold, iseed(4), mynode, row, &
                                                row_size, stored_col, &
                                                stored_row
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found, my_keep_sparsity, tr
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: rnd
    REAL(real_8), DIMENSION(:, :), POINTER   :: buff, data_d

!   ---------------------------------------------------------------------------

    my_keep_sparsity=.FALSE.
    IF (PRESENT(keep_sparsity)) my_keep_sparsity=keep_sparsity

    CALL dbcsr_error_set(routineN, error_handler, error)

    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(dbcsr_distribution(matrix)))
    CALL dbcsr_work_create(matrix, work_mutable=.TRUE., error=error)

    ALLOCATE(rnd(MAXVAL(row_blk_size)*MAXVAL(col_blk_size)))
    DO row = 1,dbcsr_nblkrows_total(matrix)
    DO col = 1,dbcsr_nblkcols_total(matrix)
       row_size = row_blk_size(row)
       col_size = col_blk_size(col)
       tr = .FALSE.
       stored_row = row
       stored_col = col
       CALL dbcsr_get_stored_coordinates(matrix, stored_row, stored_col, tr, hold)
       IF(hold.EQ.mynode) THEN
          CALL dbcsr_get_block_p(matrix, stored_row, stored_col, data_d, tr, found)
          IF(found.OR.(.NOT.my_keep_sparsity)) THEN
             ! set the seed for dlarnv, is here to guarantee same value of the random numbers
             ! for all layouts (and block distributions)
             iseed(1)=4;iseed(2)=row;iseed(3)=col;iseed(4)=1
             CALL dlarnv(1,iseed,row_size*col_size,rnd(1))
          ENDIF
          IF(found) THEN
             CALL dcopy(row_size*col_size,rnd,1,data_d,1)
          ELSE
             IF (.NOT. my_keep_sparsity) THEN
                ALLOCATE(buff(row_size,col_size))
                CALL dcopy(row_size*col_size,rnd,1,buff,1)
                CALL dbcsr_put_block (matrix, stored_row, stored_col, buff)
                DEALLOCATE(buff)
             ENDIF
          ENDIF
       ENDIF
    ENDDO
    ENDDO
    DEALLOCATE(rnd)

    CALL dbcsr_finalize(matrix, error=error)
    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_init_random

! *****************************************************************************
!> \brief get the diagonal of a dbcsr matrix
!> \param[in] matrix    the matrix
!> \param[inout] diag   the diagonal
!>
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_get_block_diag(matrix, diag, error)

    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_obj), INTENT(INOUT)           :: diag
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_block_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, error_handler, row
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_a
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_create(diag, name='diag of '//TRIM(matrix%m%name), &
         template=matrix, error=error)

    CALL dbcsr_data_init (data_a)
    CALL dbcsr_data_new (data_a, dbcsr_get_data_type(matrix))
    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_a, tr, blk)
       IF(row.EQ.col) CALL dbcsr_put_block(diag, row, col, data_a, tr)
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_data_clear_pointer (data_a)
    CALL dbcsr_data_release (data_a)
    CALL dbcsr_finalize(diag, error=error)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_get_block_diag


! *****************************************************************************
!> \brief get the diagonal of a dbcsr matrix
!> \param[in] matrix    the matrix
!> \param[inout] diag   diagonal
!>
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_get_diag_anytype(matrix, diag, error)

    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: diag
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_diag_anytype', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, data_type, &
                                                error_handler, row, &
                                                row_offset, row_size
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_a, diag_a
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (dbcsr_nfullrows_total(matrix), "LE", dbcsr_data_get_size(diag), &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Diagonal is too small",__LINE__,error)

    CALL dbcsr_data_clear (diag)
    data_type = dbcsr_get_data_type (matrix)

    CALL dbcsr_data_init (data_a)
    CALL dbcsr_data_new (data_a,&
         dbcsr_type_1d_to_2d (data_type))
    CALL dbcsr_data_init (diag_a)
    CALL dbcsr_data_new (diag_a, data_type)

    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_a, tr, blk, &
            row_offset=row_offset, row_size=row_size)
       IF(row.EQ.col) THEN
          diag_a = pointer_view (diag_a, diag, offset=row_offset, len=row_size)
          CALL get_block2d_diagonal (data_a, diag_a, row_size, error=error)
       ENDIF
    ENDDO
    CALL dbcsr_iterator_stop(iter)

    CALL dbcsr_data_clear_pointer (diag_a)
    CALL dbcsr_data_release (diag_a)
    CALL dbcsr_data_clear_pointer (data_a)
    CALL dbcsr_data_release (data_a)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_get_diag_anytype

! *****************************************************************************
!> \brief set the diagonal of a dbcsr matrix
!> \param[in] matrix    the matrix
!> \param[inout] diag   diagonal
!>
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_set_diag_anytype(matrix, diag, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_data_obj), INTENT(IN)         :: diag
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_diag_anytype', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, dt, error_handler, &
                                                row, row_offset, row_size
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_a, diag_a
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (dbcsr_nfullrows_total(matrix), "LE", dbcsr_data_get_size(diag), &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Diagonal too small",__LINE__,error)

    dt = dbcsr_get_data_type (matrix)
    CALL dbcsr_data_init (data_a)
    CALL dbcsr_data_new (data_a, dbcsr_type_1d_to_2d (dt))
    CALL dbcsr_data_init (diag_a)
    CALL dbcsr_data_new (diag_a, dt)

    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, data_a, tr, blk, &
            row_offset=row_offset, row_size=row_size)
       IF(row.EQ.col) THEN
          diag_a = pointer_view (diag_a, diag, offset=row_offset, len=row_size)
          CALL set_block2d_diagonal (data_a, diag_a, row_size, error=error)
       ENDIF
    ENDDO
    CALL dbcsr_iterator_stop(iter)

    CALL dbcsr_data_clear_pointer (diag_a)
    CALL dbcsr_data_release (diag_a)
    CALL dbcsr_data_clear_pointer (data_a)
    CALL dbcsr_data_release (data_a)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_set_diag_anytype


! *****************************************************************************
!> \brief checks if matrix symmetry and data_type are consistent
!> \brief note: does not check the symmetry of the data itself
!> \param matrix_type ...
!> \param data_type ...
!> \param error ...
!> \retval symmetry_consistent ...
! *****************************************************************************
  LOGICAL FUNCTION symmetry_consistent(matrix_type,data_type,error)
    CHARACTER, INTENT(IN)                    :: matrix_type
    INTEGER, INTENT(IN)                      :: data_type
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'symmetry_consistent', &
      routineP = moduleN//':'//routineN

    symmetry_consistent = .FALSE.

    SELECT CASE (data_type)
    CASE(dbcsr_type_real_4,dbcsr_type_real_8)
      SELECT CASE (matrix_type)
        CASE (dbcsr_type_no_symmetry, dbcsr_type_symmetric, dbcsr_type_antisymmetric)
          symmetry_consistent = .TRUE.
      END SELECT
    CASE(dbcsr_type_complex_4,dbcsr_type_complex_8)
      SELECT CASE (matrix_type)
        CASE (dbcsr_type_no_symmetry, dbcsr_type_hermitian, dbcsr_type_antihermitian)
          symmetry_consistent = .TRUE.
      END SELECT
    CASE DEFAULT
      CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
           routineN, "Invalid data type.",__LINE__,error)
    END SELECT 

  END FUNCTION symmetry_consistent

! *****************************************************************************
!> \brief checks if symmetries of two matrices are compatible for copying
!> \brief data from matrix_a(source) to matrix_b(target)
!> \param matrix_type_a ...
!> \param matrix_type_b ...
!> \param error ...
!> \retval symmetry_compatible ...
! *****************************************************************************
  LOGICAL FUNCTION symmetry_compatible(matrix_type_a,matrix_type_b,error)
    CHARACTER, INTENT(IN)                    :: matrix_type_a, matrix_type_b
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'symmetry_compatible', &
      routineP = moduleN//':'//routineN

    symmetry_compatible = .FALSE.

    SELECT CASE (matrix_type_a)
    CASE (dbcsr_type_no_symmetry)
      SELECT CASE(matrix_type_b)
      CASE(dbcsr_type_no_symmetry)
        symmetry_compatible = .TRUE.
      END SELECT
    CASE(dbcsr_type_symmetric, dbcsr_type_hermitian)
      SELECT CASE(matrix_type_b)
      CASE(dbcsr_type_symmetric, dbcsr_type_hermitian)
        symmetry_compatible = .TRUE.
      END SELECT
    CASE(dbcsr_type_antisymmetric, dbcsr_type_antihermitian)
      SELECT CASE(matrix_type_b)
      CASE(dbcsr_type_antisymmetric, dbcsr_type_antihermitian)
        symmetry_compatible = .TRUE.
      END SELECT
    CASE DEFAULT
      CALL dbcsr_assert(.FALSE., dbcsr_failure_level,&
           dbcsr_wrong_args_error, routineP, "Invalid matrix type.",__LINE__,error)
    END SELECT

  END FUNCTION symmetry_compatible


! *****************************************************************************
!> \brief copy a matrix
!> \param[inout] matrix_b       target DBCSR matrix
!> \param[in]    matrix_a       source DBCSR matrix
!> \param[in]    name           (optional) name of the new matrix
!> \param[in]    keep_sparsity  (optional) keep the target matrix sparsity;
!>                              default is False.
!> \param[in]    shallow_data   (optional) shallow data copy
!> \param[in]    keep_imaginary  (optional) when copy from complex to real,&
!>                               the default is to keep only the real part; if
!>                               this flag is set, the imaginary part is used
!> \param[in]    matrix_type     'N' for normal, 'T' for transposed, 'S' for
!>                               symmetric, and 'A' for antisymmetric
!> \param[in,out] error         cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_copy(matrix_b, matrix_a, name, keep_sparsity,&
       shallow_data, keep_imaginary, matrix_type, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity, shallow_data, &
                                                keep_imaginary
    CHARACTER, INTENT(IN), OPTIONAL          :: matrix_type
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: new_matrix_type, repl_type
    INTEGER                                  :: error_handler, new_type
    LOGICAL                                  :: keep_sparse, shallow

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (symmetry_consistent(dbcsr_get_matrix_type(matrix_a),&
                                           dbcsr_get_data_type(matrix_a), error),&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Source matrix symmetry not consistent with its data type.",__LINE__,error)
    shallow = .FALSE. ; IF (PRESENT (shallow_data)) shallow = shallow_data
    keep_sparse = .FALSE.
    IF (PRESENT (keep_sparsity)) keep_sparse = keep_sparsity
    CALL dbcsr_assert (dbcsr_is_initialized (matrix_b),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Target matrix must be initialized",__LINE__,error)
    CALL dbcsr_assert (.not.keep_sparse.or.dbcsr_valid_index(matrix_b),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Target matrix must be valid to keep its sparsity",__LINE__,error)
    CALL dbcsr_assert (.not.keep_sparse.or..not.shallow, dbcsr_warning_level,&
         dbcsr_wrong_args_error, routineN,&
         "Shallow copy not compatibly with sparsity retainment",__LINE__,error)
    IF (keep_sparse) THEN
       IF (PRESENT (name)) matrix_b%m%name = name
       CALL dbcsr_copy_into_existing(matrix_b, matrix_a, error)
    ELSE
       IF (dbcsr_valid_index(matrix_b)) THEN
          new_type = dbcsr_get_data_type (matrix_b)
          repl_type = dbcsr_get_replication_type(matrix_b)
       ELSE
          new_type = dbcsr_get_data_type (matrix_a)
          repl_type = dbcsr_get_replication_type(matrix_a)
       ENDIF
       new_matrix_type = dbcsr_get_matrix_type (matrix_a)
       IF (PRESENT (matrix_type)) THEN
         CALL dbcsr_assert (symmetry_compatible(dbcsr_get_matrix_type(matrix_a),&
                                                matrix_type, error),&
              dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
              "Specified target matrix symmetry "//matrix_type//&
               " not compatible with source matrix type "//dbcsr_get_matrix_type(matrix_a),&
              __LINE__,error)
         new_matrix_type = matrix_type
       END IF
       CALL dbcsr_assert (symmetry_consistent(new_matrix_type,new_type,error),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Target matrix symmetry "//new_matrix_type//" not consistent with its data type.",&
            __LINE__,error)
       IF(PRESENT(name)) THEN
          CALL dbcsr_create(matrix_b, name=TRIM(name), &
               template = matrix_a,&
               matrix_type = new_matrix_type,&
               data_type = new_type,&
               error=error)
       ELSE
          CALL dbcsr_create(matrix_b,&
               name='copy of '//TRIM(dbcsr_name(matrix_a)),&
               data_type = new_type,&
               matrix_type = new_matrix_type,&
               template = matrix_a, error=error)
       ENDIF
       CALL ensure_array_size(matrix_b%m%index, ub=SIZE(matrix_a%m%index),&
            memory_type = dbcsr_get_index_memory_type(matrix_b), error=error)
       !
       ! copy index and data
       matrix_b%m%index(1:SIZE(matrix_a%m%index)) = matrix_a%m%index(:)
       IF (.NOT. shallow) THEN
          CALL dbcsr_assert (matrix_a%m%nze, "LE", dbcsr_get_data_size(matrix_a),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Source matrix sizes not consistent!",__LINE__,error)
          CALL dbcsr_data_ensure_size (matrix_b%m%data_area,&
               dbcsr_get_data_size(matrix_a), error=error)
          IF (dbcsr_get_data_type(matrix_a) .EQ. dbcsr_get_data_type(matrix_b))&
               THEN
             CALL dbcsr_data_copyall (matrix_b%m%data_area,&
                  matrix_a%m%data_area, error=error)
          ELSE
             CALL dbcsr_data_convert (matrix_b%m%data_area,&
                  matrix_a%m%data_area, drop_real=keep_imaginary)
          ENDIF
       ELSE
          CALL dbcsr_assert (dbcsr_get_data_type(matrix_a) &
               .EQ. dbcsr_get_data_type(matrix_b), dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "Shallow copy only possible when retaining data type.", __LINE__, error)
          CALL dbcsr_switch_data_area (matrix_b, matrix_a%m%data_area,&
               error=error)
       ENDIF
       !
       ! the row_p, col_i and blk_p ...
       CALL dbcsr_repoint_index(matrix_b%m)
       matrix_b%m%nze = matrix_a%m%nze
       matrix_b%m%nblks = matrix_b%m%nblks
       matrix_b%m%valid = .TRUE.

       matrix_b%m%sparsity_id = matrix_a%m%sparsity_id
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_copy

! *****************************************************************************
!> \brief copy a matrix, retaining current sparsity
!> \param[inout] matrix_b       target DBCSR matrix
!> \param[in]    matrix_a       source DBCSR matrix
!> \param[in,out] error         cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_copy_into_existing(matrix_b, matrix_a, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_into_existing', &
      routineP = moduleN//':'//routineN

    INTEGER :: col_size, data_type, dst_col, dst_row, error_handler, rel, &
      row_size, src_col, src_cs, src_row, src_rs
    LOGICAL                                  :: dst_tr, making_symmetric, &
                                                neg_imag, neg_real, src_tr
    TYPE(dbcsr_data_obj)                     :: dst_data, src_data
    TYPE(dbcsr_iterator)                     :: dst_iter, src_iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_get_data_type(matrix_b)&
         .EQ. dbcsr_get_data_type(matrix_a), dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Matrices have different data types.",__LINE__,error)
    data_type = dbcsr_get_data_type (matrix_b)
    neg_real = matrix_b%m%negate_real
    neg_imag = matrix_b%m%negate_imaginary
    making_symmetric = dbcsr_has_symmetry (matrix_b)&
         .AND. .NOT. dbcsr_has_symmetry (matrix_a)
    IF (making_symmetric) THEN
       CALL dbcsr_copy_into_existing_sym(matrix_b, matrix_a, error)
       CALL dbcsr_error_stop(error_handler, error)
       RETURN
    ENDIF
    CALL dbcsr_data_init (src_data)
    CALL dbcsr_data_init (dst_data)
    CALL dbcsr_data_new (src_data, data_type)
    CALL dbcsr_data_new (dst_data, data_type)
    CALL dbcsr_iterator_start (src_iter, matrix_a)
    CALL dbcsr_iterator_start (dst_iter, matrix_b)
    ! Iterate through the blocks of the source and destination
    ! matrix. There are three possibilites: 1. copy the data for
    ! blocks present in both; 2 skip source blocks not present in the
    ! target; 3 zero blocks not present in the source.
    IF (dbcsr_iterator_blocks_left (src_iter)) THEN
       CALL dbcsr_iterator_next_block (src_iter, src_row, src_col, src_data,&
            src_tr)
    ELSE
       src_row = 0 ; src_col = 0
    ENDIF
    DO WHILE (dbcsr_iterator_blocks_left (dst_iter))
       CALL dbcsr_iterator_next_block (dst_iter, dst_row, dst_col, dst_data,&
            dst_tr, row_size=row_size, col_size=col_size)
       ! Now find the source position that is greater or equal to the
       ! target one. I.e, skip blocks that the target doesn't have.
       rel = pos_relation (dst_row, dst_col, src_row, src_col)
       DO WHILE (rel .EQ. 1 .AND. dbcsr_iterator_blocks_left (src_iter))
          CALL dbcsr_iterator_next_block (src_iter, src_row, src_col,&
               src_data, src_tr, row_size=src_rs, col_size=src_cs)
          rel = pos_relation (dst_row, dst_col, src_row, src_col)
       ENDDO
       SELECT CASE (rel)
       CASE (-1, 1)
           ! Target lags source or ran out of source
          CALL dbcsr_data_clear (dst_data)
       CASE (0)
          ! Copy the data
          CALL dbcsr_assert (dbcsr_data_get_size (src_data)&
               .EQ. dbcsr_data_get_size (dst_data), dbcsr_fatal_level,&
               dbcsr_internal_error, routineN, "Block sizes not equal!",__LINE__,error)
          IF (src_tr .EQV. dst_tr) THEN
             CALL dbcsr_data_copyall (dst_data, src_data, error=error)
          ELSE
             CALL dbcsr_block_partial_copy(dst=dst_data, dst_tr = dst_tr,&
                  dst_rs = row_size, dst_cs = col_size,&
                  dst_r_lb = 1, dst_c_lb = 1,&
                  src = src_data, src_tr = src_tr,&
                  src_rs = src_rs, src_cs = src_cs,&
                  src_r_lb = 1, src_c_lb = 1,&
                  nrow = row_size, ncol = col_size)
             IF (neg_real) THEN
                CALL dbcsr_block_real_neg (dst_data, row_size, col_size, error=error)
             ENDIF
             IF (neg_imag) THEN
                CALL dbcsr_block_conjg (dst_data, row_size, col_size, error=error)
             ENDIF
          ENDIF
       CASE default
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
               routineN, "Trouble syncing iterators",__LINE__,error)
       END SELECT
    END DO
    CALL dbcsr_iterator_stop (src_iter)
    CALL dbcsr_iterator_stop (dst_iter)
    CALL dbcsr_data_clear_pointer (src_data)
    CALL dbcsr_data_clear_pointer (dst_data)
    CALL dbcsr_data_release (src_data)
    CALL dbcsr_data_release (dst_data)
    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_copy_into_existing

! *****************************************************************************
!> \brief copy a matrix, retaining current sparsity
!> \param[inout] matrix_b       target DBCSR matrix
!> \param[in]    matrix_a       source DBCSR matrix
!> \param[in,out] error         cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_copy_into_existing_sym(matrix_b, matrix_a, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_into_existing_sym', &
      routineP = moduleN//':'//routineN

    INTEGER :: col_size, data_type, dst_col, dst_row, error_handler, &
      row_size, src_col, src_cs, src_row, src_rs
    LOGICAL                                  :: dst_tr, found, neg_imag, &
                                                neg_real, src_tr
    TYPE(dbcsr_data_obj)                     :: dst_data, src_data
    TYPE(dbcsr_iterator)                     :: dst_iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_get_data_type(matrix_b)&
         .EQ. dbcsr_get_data_type(matrix_a), dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Matrices have different data types.",__LINE__,error)
    data_type = dbcsr_get_data_type (matrix_b)
    CALL dbcsr_assert (dbcsr_has_symmetry (matrix_b)&
         .AND. .NOT. dbcsr_has_symmetry (matrix_a),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Must copy from non-symmetric to symmetric matrix.",&
         __LINE__, error=error)
    neg_real = matrix_b%m%negate_real
    neg_imag = matrix_b%m%negate_imaginary

    CALL dbcsr_data_init (src_data)
    CALL dbcsr_data_init (dst_data)
    CALL dbcsr_data_new (src_data, data_type)
    CALL dbcsr_data_new (dst_data, data_type)
    CALL dbcsr_iterator_start (dst_iter, matrix_b)
    ! Iterate through the blocks of the destination matrix.  For each
    ! one, try to find an appropriate source matrix block and copy it
    ! into the destination matrix.
    DO WHILE (dbcsr_iterator_blocks_left (dst_iter))
       CALL dbcsr_iterator_next_block (dst_iter, dst_row, dst_col, dst_data,&
            dst_tr, row_size=row_size, col_size=col_size)
       src_row = dst_row
       src_col = dst_col
       IF (checker_tr (dst_row, dst_col))&
            CALL swap (src_row, src_col)
       CALL dbcsr_get_block_p (matrix_a, src_row, src_col, src_data, src_tr,&
            found=found, row_size=src_rs, col_size=src_cs)
       IF (.NOT. found) THEN
          CALL dbcsr_data_clear (dst_data)
       ELSE
          IF (checker_tr (dst_row, dst_col)) THEN
             src_tr = .NOT. src_tr
             CALL swap (src_rs, src_cs)
          ENDIF
          CALL dbcsr_block_partial_copy(dst=dst_data, dst_tr = dst_tr,&
               dst_rs = row_size, dst_cs = col_size,&
               dst_r_lb = 1, dst_c_lb = 1,&
               src = src_data, src_tr = src_tr,&
               src_rs = src_rs, src_cs = src_cs,&
               src_r_lb = 1, src_c_lb = 1,&
               nrow = row_size, ncol = col_size)
          IF (neg_real) THEN
             CALL dbcsr_block_real_neg (dst_data, row_size, col_size, error=error)
          ENDIF
          IF (neg_imag) THEN
             CALL dbcsr_block_conjg (dst_data, row_size, col_size, error=error)
          ENDIF
       ENDIF
    END DO
    CALL dbcsr_iterator_stop (dst_iter)
    CALL dbcsr_data_clear_pointer (src_data)
    CALL dbcsr_data_clear_pointer (dst_data)
    CALL dbcsr_data_release (src_data)
    CALL dbcsr_data_release (dst_data)
    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_copy_into_existing_sym

! *****************************************************************************
!> \brief Determines the relation between two matrix positions.
!> \param row1 ...
!> \param col1 ...
!> \param row2 ...
!> \param col2 ...
!> \retval relation  Relation between positions 1 and 2.
!>                    0: same
!>                   -1: pos1 < pos2
!>                    1: pos1 > pos2
! *****************************************************************************
  ELEMENTAL FUNCTION pos_relation (row1, col1, row2, col2) RESULT (relation)
    INTEGER, INTENT(IN)                      :: row1, col1, row2, col2
    INTEGER                                  :: relation

    IF (row1 .LT. row2) THEN
       relation = -1
    ELSEIF (row1 .GT. row2) THEN
       relation = 1
    ELSE ! rows are equal, check column
       IF (col1 .LT. col2) THEN
          relation = -1
       ELSEIF (col1 .GT. col2) THEN
          relation = 1
       ELSE
          relation = 0
       ENDIF
    ENDIF
  END FUNCTION pos_relation


! *****************************************************************************
!> \brief Copy a submatrix.
!> \param[inout] matrix_b       target DBCSR matrix
!> \param[in]    matrix_a       source DBCSR matrix
!> \param[in]    name           (optional) name of the new matrix
!> \param[in]    block_row_bounds     (optional) rows to extract (array of
!>                                    size 2 holding the lower and upper
!>                                    inclusive bounds)
!> \param[in]    block_column_bounds  (optional) columns to extract (array of
!>                                    size 2 holding the lower and upper
!>                                    inclusive bounds)
!> \param[in]    shallow_data   (optional) shallow data copy
!> \param[in,out] error         cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_copy_submatrix(matrix_b, matrix_a, name,&
       block_row_bounds, block_column_bounds, &
       shallow_data, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name
    INTEGER, DIMENSION(2), INTENT(IN), &
      OPTIONAL                               :: block_row_bounds, &
                                                block_column_bounds
    LOGICAL, INTENT(IN), OPTIONAL            :: shallow_data
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_submatrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_p, col, error_handler, &
                                                nblocks, new_blk, old_blk, row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: blkp_list, col_list, row_list
    LOGICAL                                  :: shallow, tr
    TYPE(dbcsr_data_obj)                     :: data_block
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    IF (PRESENT (shallow_data)) THEN
       shallow = shallow_data
    ELSE
       shallow = .FALSE.
    ENDIF
    ! Verify assumptions.
    IF (PRESENT (block_row_bounds)) THEN
       CALL dbcsr_assert (SIZE(block_row_bounds), "EQ", 2, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Size of bounds specifier must be 2", __LINE__, error=error)
    ENDIF
    IF (PRESENT (block_column_bounds)) THEN
       CALL dbcsr_assert (SIZE(block_column_bounds), "EQ", 2,&
            dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Size of bounds specifier must be 2", __LINE__, error=error)
    ENDIF
    ! Setup target matrix
    CALL dbcsr_create (matrix_b, name=name, template=matrix_a, error=error)
    CALL dbcsr_finalize(matrix_b, error=error)
    IF (.NOT. shallow) THEN
       ! Non-shallow copy uses the standard iterator on the source and
       ! block put on the target.
!
!$OMP PARALLEL DEFAULT (none) &
!$OMP          PRIVATE (data_block, iter, row, col, tr) &
!$OMP          SHARED (matrix_a, matrix_b, error,&
!$OMP                  block_row_bounds, block_column_bounds)
       CALL dbcsr_work_create (matrix_b, work_mutable=.FALSE., error=error)
       CALL dbcsr_data_init (data_block)
       CALL dbcsr_data_new (data_block, dbcsr_get_data_type (matrix_a))
       CALL dbcsr_iterator_start (iter, matrix_a, dynamic=.TRUE.,&
            dynamic_byrows = .TRUE.)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, row, col, data_block, tr)
          ! Only keep the block if they are within the specified bounds.
          IF (PRESENT (block_row_bounds)) THEN
             IF (row .LT. block_row_bounds(1)) CYCLE
             IF (row .GT. block_row_bounds(2)) CYCLE
          ENDIF
          IF (PRESENT (block_column_bounds)) THEN
             IF (col .LT. block_column_bounds(1)) CYCLE
             IF (col .GT. block_column_bounds(2)) CYCLE
          ENDIF
          CALL dbcsr_put_block (matrix_b, row, col, data_block, transposed=tr)
       END DO
       CALL dbcsr_iterator_stop (iter)
       CALL dbcsr_data_clear_pointer (data_block)
       CALL dbcsr_data_release (data_block)
       CALL dbcsr_finalize (matrix_b, error=error)
!$OMP END PARALLEL
    ELSE
       ! For the shallow copy the source matrix data is referenced.
       CALL dbcsr_switch_data_area (matrix_b, matrix_a%m%data_area, error=error)
       nblocks = dbcsr_get_num_blocks (matrix_a) ! High estimate.
       ! Shallow copy goes through source's data blocks and inserts
       ! the only the ones corresponding to the submatrix specifier
       ! into the target. Block pointers must remain the same as in
       ! the source.
       ALLOCATE (row_list(nblocks), col_list(nblocks), blkp_list(nblocks))
       !
       CALL dbcsr_iterator_start (iter, matrix_a)
       new_blk = 1
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, row, col,&
               blk=old_blk, blk_p=blk_p)
          ! Only keep the block if they are within the specified bounds.
          IF (PRESENT (block_row_bounds)) THEN
             IF (row .LT. block_row_bounds(1)) CYCLE
             IF (row .GT. block_row_bounds(2)) CYCLE
          ENDIF
          IF (PRESENT (block_column_bounds)) THEN
             IF (col .LT. block_column_bounds(1)) CYCLE
             IF (col .GT. block_column_bounds(2)) CYCLE
          ENDIF
          row_list(new_blk) = row
          col_list(new_blk) = col
          blkp_list(new_blk) = blk_p
          new_blk = new_blk + 1
       ENDDO
       new_blk = new_blk - 1
       CALL dbcsr_iterator_stop (iter)
       CALL dbcsr_reserve_blocks (matrix_b, row_list(1:new_blk),&
            col_list(1:new_blk), blkp_list(1:new_blk), error=error)
    ENDIF
    !
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_copy_submatrix


! *****************************************************************************
!> \brief Crop and copies a matrix.
!> \param[inout] matrix_b       target DBCSR matrix
!> \param[in]    matrix_a       source DBCSR matrix
!> \param[in]    name           (optional) name of the new matrix
!> \param[in]    full_row_bounds      (optional) rows to extract (array of
!>                                    size 2 holding the lower and upper
!>                                    inclusive bounds)
!> \param[in]    full_column_bounds   (optional) columns to extract (array of
!>                                    size 2 holding the lower and upper
!>                                    inclusive bounds)
!> \param shallow_data ...
!> \param[in,out] error         cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_crop_matrix(matrix_b, matrix_a, name,&
       full_row_bounds, full_column_bounds, &
       shallow_data, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name
    INTEGER, DIMENSION(2), INTENT(IN), &
      OPTIONAL                               :: full_row_bounds, &
                                                full_column_bounds
    LOGICAL, INTENT(IN), OPTIONAL            :: shallow_data
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_crop_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, error_handler, f_col_f, &
                                                f_row_f, l_col_l, l_row_l, row
    INTEGER, DIMENSION(2)                    :: block_col_bounds, &
                                                block_row_bounds
    LOGICAL                                  :: part_col, part_f_col, &
                                                part_f_row, part_l_col, &
                                                part_l_row, part_row, &
                                                shallow, tr
    TYPE(dbcsr_data_obj)                     :: data_block
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    IF (PRESENT (shallow_data)) THEN
       shallow = shallow_data
    ELSE
       shallow = .FALSE.
    ENDIF
    block_row_bounds = 0
    block_col_bounds = 0
    part_col = .FALSE.
    part_row = .FALSE.
    !
    ! If row bounds are present, they must be converted to block
    ! addressing.
    IF (PRESENT (full_row_bounds)) THEN
       CALL dbcsr_assert (SIZE(full_row_bounds), "EQ", 2, dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Size of bounds specifier must be 2", __LINE__, error=error)
       CALL dbcsr_assert (full_row_bounds(1), "GE", 0,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid first row bound.", __LINE__, error=error)
       CALL dbcsr_assert (full_row_bounds(2), "LE", dbcsr_nfullrows_total(matrix_a),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid last row bound.", __LINE__, error=error)
       part_f_row = .FALSE.
       IF (full_row_bounds(1) .EQ. 0) THEN
          block_row_bounds(1) = 1
       ELSE
          CALL find_block_of_element (full_row_bounds(1), block_row_bounds(1),&
               dbcsr_nblkrows_total (matrix_a),&
               array_data (dbcsr_row_block_offsets (matrix_a)),&
               hint = 0, error=error)
          part_f_row = array_get (dbcsr_row_block_offsets (matrix_a), block_row_bounds(1))&
               .NE. full_row_bounds(1)
       ENDIF
       f_row_f = -7
       IF (part_f_row) THEN
          ! Block offset of last cleared row
          f_row_f = full_row_bounds(1) -&
               array_get (dbcsr_row_block_offsets (matrix_a), block_row_bounds(1))
       ENDIF
       part_l_row = .FALSE.
       IF (full_row_bounds(2) .EQ. 0) THEN
          block_row_bounds(2) = dbcsr_nblkrows_total (matrix_a)
       ELSE
          CALL find_block_of_element (full_row_bounds(2), block_row_bounds(2),&
               dbcsr_nblkrows_total (matrix_a),&
               array_data (dbcsr_row_block_offsets (matrix_a)),&
               hint = 0, error=error)
          part_l_row = array_get (dbcsr_row_block_offsets (matrix_a), block_row_bounds(2)+1)-1&
               .NE. full_row_bounds(2)
       ENDIF
       ! Block offset of first cleared row
       l_row_l = -7
       IF (part_l_row) THEN
          l_row_l = 2 + full_row_bounds(2) - &
               array_get (dbcsr_row_block_offsets (matrix_a), block_row_bounds(2))
       ENDIF
       part_row = part_f_row .OR. part_l_row
    ENDIF
    !
    ! If column bounds are present, they must be converted to block
    ! addressing.
    IF (PRESENT (full_column_bounds)) THEN
       CALL dbcsr_assert (SIZE(full_column_bounds), "EQ", 2,&
            dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Size of bounds specifier must be 2", __LINE__, error=error)
       CALL dbcsr_assert (full_column_bounds(1), "GE", 0,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid first column bound.", __LINE__, error=error)
       CALL dbcsr_assert (full_column_bounds(2), "LE", dbcsr_nfullcols_total(matrix_a),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid last column bound.", __LINE__, error=error)
       part_f_col = .FALSE.
       IF (full_column_bounds(1) .EQ. 0) THEN
          block_col_bounds(1) = 1
       ELSE
          CALL find_block_of_element (full_column_bounds(1), block_col_bounds(1),&
               dbcsr_nblkcols_total (matrix_a),&
               array_data (dbcsr_col_block_offsets (matrix_a)),&
               hint=0, error=error)
          part_f_col = array_get (dbcsr_col_block_offsets (matrix_a), block_col_bounds(1))&
               .NE. full_column_bounds(1)
       ENDIF
       f_col_f = -7
       IF (part_f_col) THEN
          ! Block offset of last cleared column
          f_col_f = full_column_bounds(1) -&
               array_get (dbcsr_col_block_offsets (matrix_a), block_col_bounds(1))
       ENDIF
       part_l_col = .FALSE.
       IF (full_column_bounds(2) .EQ. 0) THEN
          block_col_bounds(2) = dbcsr_nblkcols_total (matrix_a)
       ELSE
          CALL find_block_of_element (full_column_bounds(2), block_col_bounds(2),&
               dbcsr_nblkcols_total (matrix_a),&
               array_data (dbcsr_col_block_offsets (matrix_a)),&
               hint=0, error=error)
          part_l_col = array_get (dbcsr_col_block_offsets (matrix_a), block_col_bounds(2)+1)-1&
               .NE. full_column_bounds(2)
       ENDIF
       l_col_l = -7
       IF (part_l_col) THEN
          ! Block offset of first cleared column
          l_col_l = 2 + full_column_bounds(2) - &
               array_get (dbcsr_col_block_offsets (matrix_a), block_col_bounds(2))
       ENDIF
       part_col = part_f_col .OR. part_l_col
    ENDIF
    !
    ! First copy the blocks then perform the intra-block zeroing.
    CALL dbcsr_copy_submatrix (matrix_b, matrix_a,&
         block_row_bounds=block_row_bounds,&
         block_column_bounds=block_col_bounds,&
         shallow_data=shallow, error=error)
    IF (part_row .OR. part_col) THEN
!$OMP PARALLEL DEFAULT (none) &
!$OMP          PRIVATE (data_block, iter, row, col, tr) &
!$OMP          SHARED (matrix_b, error,&
!$OMP                  part_row, part_f_row, part_l_row, f_row_f, l_row_l, &
!$OMP                  part_col, part_f_col, part_l_col, f_col_f, l_col_l,&
!$OMP                  block_row_bounds, block_col_bounds)
       CALL dbcsr_data_init (data_block)
       CALL dbcsr_data_new (data_block, dbcsr_type_1d_to_2d (dbcsr_get_data_type (matrix_b)))
       CALL dbcsr_iterator_start (iter, matrix_b, &
            dynamic=.TRUE., dynamic_byrows = .TRUE.)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, row, col, data_block, tr)
          IF (part_row) THEN
             IF (row .LT. block_row_bounds(1)) CYCLE
             IF (row .GT. block_row_bounds(2)) CYCLE
          ENDIF
          IF (part_col) THEN
             IF (col .LT. block_col_bounds(1)) CYCLE
             IF (col .GT. block_col_bounds(2)) CYCLE
          ENDIF
          IF (part_row) THEN
             IF (part_f_row .AND. row .EQ. block_row_bounds(1)) THEN
                CALL dbcsr_data_clear (data_block, ub=f_row_f, tr=tr)
             ENDIF
             IF (part_l_row .AND. row .EQ. block_row_bounds(2)) THEN
                CALL dbcsr_data_clear (data_block, lb=l_row_l, tr=tr)
             ENDIF
          ENDIF
          IF (part_col) THEN
             IF (part_f_col .AND. col .EQ. block_col_bounds(1)) THEN
                CALL dbcsr_data_clear (data_block, ub2=f_col_f, tr=tr)
             ENDIF
             IF (part_l_col .AND. col .EQ. block_col_bounds(2)) THEN
                CALL dbcsr_data_clear (data_block, lb2=l_col_l, tr=tr)
             ENDIF
          ENDIF
       END DO
       CALL dbcsr_iterator_stop (iter)
       CALL dbcsr_data_clear_pointer (data_block)
       CALL dbcsr_data_release (data_block)
       CALL dbcsr_finalize (matrix_b, error=error)
!$OMP END PARALLEL
    ENDIF
    !
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_crop_matrix




! *****************************************************************************
!> \brief  triu of a dbcsr matrix
!> \param[inout] matrix_a  the matrix
!>
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_triu(matrix_a, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_triu', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_nze, col, col_size, &
                                                error_handler, i, j, row, &
                                                row_size
    LOGICAL                                  :: tr
    REAL(dp), DIMENSION(:, :), POINTER       :: DATA
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_iterator_start(iter, matrix_a)

    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, DATA, tr, &
            block_number=blk, row_size=row_size, col_size=col_size)
       blk_nze = row_size * col_size
       IF(row.GT.col) CALL dbcsr_remove_block(matrix_a, row, col, blk_nze, blk)
       IF(row.EQ.col) THEN
          DO j = 1,col_size
          DO i = j+1,row_size
             DATA(i,j) = 0.0_dp
          ENDDO
          ENDDO
       ENDIF
    ENDDO

    CALL dbcsr_iterator_stop(iter)

    CALL dbcsr_finalize (matrix_a, error=error)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_triu

! *****************************************************************************
!> \brief filter a dbcsr matrix
!> \param[inout] matrix the matrix
!> \param[in] eps          the threshold
!> \param[in] method       how the matrix is filtered
!> \param[in] use_absolute ?
!> \param[in] filter_diag  NYI
!> \param[in] quick        (optional) filter just the index (no data copying)
!> \param[error] error     error
! *****************************************************************************
  SUBROUTINE dbcsr_filter_anytype(matrix, eps, method,&
       use_absolute, filter_diag, quick, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: eps
    INTEGER, INTENT(IN), OPTIONAL            :: method
    LOGICAL, INTENT(in), OPTIONAL            :: use_absolute, filter_diag, &
                                                quick
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_filter_anytype', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: data_c
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: data_z
    INTEGER                                  :: blk, blk_nze, col, col_size, &
                                                error_handler, my_method, &
                                                row, row_size
    LOGICAL                                  :: gt0, my_filter_diag, tr
    REAL(KIND=real_4)                        :: nrm_s, sdot
    REAL(KIND=real_4), DIMENSION(:), POINTER :: data_s
    REAL(KIND=real_8)                        :: ddot, my_absolute, nrm_d
    REAL(KIND=real_8), DIMENSION(:), POINTER :: data_d
    REAL(KIND=real_8), EXTERNAL              :: DZNRM2, SCNRM2
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    my_method = dbcsr_filter_frobenius
    IF(PRESENT(method)) my_method = method
    my_absolute = 1.0_dp
    IF(PRESENT(use_absolute)) my_absolute = dbcsr_maxabs(matrix)
    my_filter_diag = .TRUE.
    IF(PRESENT(filter_diag)) my_filter_diag = filter_diag

    SELECT CASE (eps%data_type)
    CASE (dbcsr_type_real_4)
       gt0 = eps%r_sp .GT. 0.0_real_4
    CASE (dbcsr_type_real_8)
       gt0 = eps%r_dp .GT. 0.0_real_8
    CASE (dbcsr_type_complex_4)
       gt0 = ABS(eps%c_sp) .GT. 0.0_real_4
    CASE (dbcsr_type_complex_8)
       gt0 = ABS(eps%c_dp) .GT. 0.0_real_8
    CASE default
       gt0 = .FALSE.
    END SELECT

    IF(gt0) THEN

       !$OMP PARALLEL DEFAULT(NONE) PRIVATE(iter,row,col,data_s,data_d,data_c,data_z,tr, &
       !$OMP                                blk,row_size,col_size,blk_nze,nrm_d,nrm_s) &
       !$OMP                        SHARED(my_method,quick,error,my_absolute,eps,matrix)

       CALL dbcsr_iterator_start(iter, matrix, contiguous_pointers=.TRUE.)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          SELECT CASE (dbcsr_get_data_type(matrix))

          CASE (dbcsr_type_real_4)
             CALL dbcsr_iterator_next_block(iter, row, col, data_s, tr, blk, &
                  row_size, col_size)
             blk_nze = row_size * col_size
             IF (blk_nze .EQ. 0) CYCLE ! Skip empty blocks
             SELECT CASE(my_method)
             CASE(dbcsr_filter_frobenius)
                !
                ! Frobenius based
                nrm_s = SQRT(SDOT(SIZE(data_s), data_s(1), 1, data_s(1), 1))
                IF(nrm_s.LT.my_absolute*eps%r_sp) &
                     CALL dbcsr_remove_block(matrix, row, col, blk_nze, blk)
             CASE DEFAULT
                CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
                     routineN,"Only Frobenius based filtering",__LINE__,error)
             END SELECT

          CASE (dbcsr_type_real_8)
             CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk, &
                  row_size, col_size)
             blk_nze = row_size * col_size
             IF (blk_nze .EQ. 0) CYCLE ! Skip empty blocks
             SELECT CASE(my_method)
             CASE(dbcsr_filter_frobenius)
                !
                ! Frobenius based
                nrm_d = SQRT(DDOT(SIZE(data_d), data_d(1), 1, data_d(1), 1))
                IF(nrm_d.LT.my_absolute*eps%r_dp) &
                     CALL dbcsr_remove_block(matrix, row, col, blk_nze, blk)
             CASE DEFAULT
                CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
                     routineN,"Only Frobenius based filtering",__LINE__,error)
             END SELECT

          CASE (dbcsr_type_complex_4)
               CALL dbcsr_iterator_next_block(iter, row, col, data_c, tr, blk, &
                  row_size, col_size)
             blk_nze = row_size * col_size
             IF (blk_nze .EQ. 0) CYCLE ! Skip empty blocks
             SELECT CASE(my_method)
             CASE(dbcsr_filter_frobenius)
                !
                ! Frobenius based
                nrm_d = SCNRM2(SIZE(data_c), data_c(1), 1)
                IF(nrm_d.LT.my_absolute*eps%r_dp) &
                     CALL dbcsr_remove_block(matrix, row, col, blk_nze, blk)
             CASE DEFAULT
                CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
                     routineN,"Only Frobenius based filtering",__LINE__,error)
             END SELECT

          CASE (dbcsr_type_complex_8)
             CALL dbcsr_iterator_next_block(iter, row, col, data_z, tr, blk, &
                  row_size, col_size)
             blk_nze = row_size * col_size
             IF (blk_nze .EQ. 0) CYCLE ! Skip empty blocks
             SELECT CASE(my_method)
             CASE(dbcsr_filter_frobenius)
                !
                ! Frobenius based
                nrm_d = DZNRM2(SIZE(data_z), data_z(1), 1)
                IF(nrm_d.LT.my_absolute*eps%r_dp) &
                     CALL dbcsr_remove_block(matrix, row, col, blk_nze, blk)
             CASE DEFAULT
                CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
                     routineN,"Only Frobenius based filtering",__LINE__,error)
             END SELECT

          CASE DEFAULT
             CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error, &
                  routineN,"Wrong data type",__LINE__,error)
          END SELECT
       ENDDO
       CALL dbcsr_iterator_stop(iter)
       IF (PRESENT (quick)) THEN
          CALL dbcsr_finalize (matrix, reshuffle=.NOT. quick, error=error)
       ELSE
          CALL dbcsr_finalize (matrix, reshuffle=.TRUE., error=error)
       ENDIF

       !$OMP END PARALLEL

    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_filter_anytype

! *****************************************************************************
!> \brief compute a norm of a dbcsr matrix
!> \param[in] matrix  the matrix
!> \param[in] which_norm ...
!>
!> \param norm_scalar ...
!> \param norm_vector ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_norm_anytype(matrix, which_norm, norm_scalar, norm_vector, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: which_norm
    REAL(KIND=real_8), INTENT(OUT), OPTIONAL :: norm_scalar
    TYPE(dbcsr_data_obj), INTENT(INOUT), &
      OPTIONAL                               :: norm_vector
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_norm_anytype', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_offset, &
                                                error_handler, i, j, row, &
                                                row_offset
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_a
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    SELECT CASE(which_norm)

    CASE(dbcsr_norm_frobenius)

       IF(PRESENT(norm_scalar)) norm_scalar = dbcsr_frobenius_norm(matrix)

    CASE(dbcsr_norm_maxabsnorm)

       IF(PRESENT(norm_scalar)) norm_scalar = dbcsr_maxabs(matrix)

    CASE(dbcsr_norm_gershgorin)

       IF(PRESENT(norm_scalar)) norm_scalar = dbcsr_gershgorin_norm(matrix)

    CASE(dbcsr_norm_column)

       IF(PRESENT(norm_vector)) THEN
          CALL dbcsr_assert (dbcsr_data_get_type(norm_vector), "EQ",&
               dbcsr_get_data_type(matrix), dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "Mismatched vector/matrix data types", __LINE__, error=error)
          IF(dbcsr_has_symmetry(matrix)) THEN
            CALL dbcsr_assert (dbcsr_data_get_size(norm_vector), "GE",&
                 dbcsr_nfullrows_total(matrix), dbcsr_fatal_level,&
                 dbcsr_wrong_args_error, routineN, "Passed vector too small",&
                 __LINE__, error=error)
          END IF
          CALL dbcsr_assert (dbcsr_data_get_size(norm_vector), "GE",&
               dbcsr_nfullcols_total(matrix), dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN, "Passed vector too small",&
               __LINE__, error=error)
          CALL dbcsr_data_init (data_a)
          CALL dbcsr_data_new (data_a, dbcsr_type_1d_to_2d(dbcsr_get_data_type(matrix)))
          CALL dbcsr_data_clear (norm_vector)
          CALL dbcsr_iterator_start(iter, matrix)
          DO WHILE (dbcsr_iterator_blocks_left(iter))
             CALL dbcsr_iterator_next_block(iter, row, col, data_a, tr,&
                  blk, row_offset=row_offset, col_offset=col_offset)
             SELECT CASE (dbcsr_get_data_type(matrix))
             CASE (dbcsr_type_real_4)
                IF(dbcsr_has_symmetry(matrix) .AND. row.NE.col) THEN
                   DO j=1,SIZE(data_a%d%r2_sp,2)
                   DO i=1,SIZE(data_a%d%r2_sp,1)
                      norm_vector%d%r_sp(col_offset+j-1) &
                           = norm_vector%d%r_sp(col_offset+j-1) &
                           + data_a%d%r2_sp(i,j)**2
                      norm_vector%d%r_sp(row_offset+i-1) &
                           = norm_vector%d%r_sp(row_offset+i-1) &
                           + data_a%d%r2_sp(i,j)**2
                   ENDDO
                   ENDDO
                ELSE
                   DO j=1,SIZE(data_a%d%r2_sp,2)
                   DO i=1,SIZE(data_a%d%r2_sp,1)
                      norm_vector%d%r_sp(col_offset+j-1) &
                           = norm_vector%d%r_sp(col_offset+j-1) &
                           + data_a%d%r2_sp(i,j) * data_a%d%r2_sp(i,j)
                   ENDDO
                   ENDDO
                ENDIF
             CASE (dbcsr_type_real_8)
                IF(dbcsr_has_symmetry(matrix) .AND. row.NE.col) THEN
                   DO j=1,SIZE(data_a%d%r2_dp,2)
                   DO i=1,SIZE(data_a%d%r2_dp,1)
                      norm_vector%d%r_dp(col_offset+j-1) &
                           = norm_vector%d%r_dp(col_offset+j-1) &
                           + data_a%d%r2_dp(i,j)**2
                      norm_vector%d%r_dp(row_offset+i-1) &
                           = norm_vector%d%r_dp(row_offset+i-1) &
                           + data_a%d%r2_dp(i,j)**2
                   ENDDO
                   ENDDO
                ELSE
                   DO j=1,SIZE(data_a%d%r2_dp,2)
                   DO i=1,SIZE(data_a%d%r2_dp,1)
                      norm_vector%d%r_dp(col_offset+j-1) &
                           = norm_vector%d%r_dp(col_offset+j-1) &
                           + data_a%d%r2_dp(i,j) * data_a%d%r2_dp(i,j)
                   ENDDO
                   ENDDO
                ENDIF
             !CASE (dbcsr_type_complex_4)
             !CASE (dbcsr_type_complex_8)
             CASE DEFAULT
                CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
                     routineN,"Only real values",__LINE__,error)
             END SELECT
          ENDDO
          CALL dbcsr_iterator_stop(iter)
          CALL dbcsr_data_clear_pointer (data_a)
          CALL dbcsr_data_release (data_a)
          SELECT CASE (dbcsr_get_data_type(matrix))
          CASE (dbcsr_type_real_4)
             CALL mp_sum(norm_vector%d%r_sp,&
                  dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))
             norm_vector%d%r_sp = SQRT(norm_vector%d%r_sp)
          CASE (dbcsr_type_real_8)
             CALL mp_sum(norm_vector%d%r_dp,&
                  dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))
             norm_vector%d%r_dp = SQRT(norm_vector%d%r_dp)
          END SELECT
       ENDIF

    CASE DEFAULT

       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, &
            "this norm is NYI",__LINE__,error)

    END SELECT
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_norm_anytype

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param which_norm ...
!> \param norm_vector ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_norm_r8_vec(matrix, which_norm, norm_vector, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: which_norm
    REAL(KIND=real_8), DIMENSION(:), &
      INTENT(OUT), TARGET                    :: norm_vector
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    REAL(KIND=real_8), DIMENSION(:), POINTER :: v_p
    TYPE(dbcsr_data_obj)                     :: norm_vector_a

    CALL dbcsr_data_init (norm_vector_a)
    CALL dbcsr_data_new (norm_vector_a, dbcsr_type_real_8)
    v_p => norm_vector
    CALL dbcsr_data_set_pointer (norm_vector_a, v_p)
    CALL dbcsr_norm_anytype (matrix, which_norm, norm_vector=norm_vector_a,&
         error=error)
    CALL dbcsr_data_clear_pointer (norm_vector_a)
    CALL dbcsr_data_release (norm_vector_a)
  END SUBROUTINE dbcsr_norm_r8_vec

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param which_norm ...
!> \param norm_scalar ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_norm_r4_scal(matrix, which_norm, norm_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: which_norm
    REAL(KIND=real_4), INTENT(OUT)           :: norm_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    REAL(KIND=real_8)                        :: norm_r8

    CALL dbcsr_norm_anytype (matrix, which_norm, norm_scalar=norm_r8,&
         error=error)
    norm_scalar = REAL(norm_r8, KIND=real_4)
  END SUBROUTINE dbcsr_norm_r4_scal

! *****************************************************************************
!> \brief ...
!> \param matrix ...
!> \param which_norm ...
!> \param norm_vector ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_norm_r4_vec(matrix, which_norm, norm_vector, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: which_norm
    REAL(KIND=real_4), DIMENSION(:), &
      INTENT(OUT), TARGET                    :: norm_vector
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    REAL(KIND=real_4), DIMENSION(:), POINTER :: v_p
    TYPE(dbcsr_data_obj)                     :: norm_vector_a

    CALL dbcsr_data_init (norm_vector_a)
    CALL dbcsr_data_new (norm_vector_a, dbcsr_type_real_4)
    v_p => norm_vector
    CALL dbcsr_data_set_pointer (norm_vector_a, v_p)
    CALL dbcsr_norm_anytype (matrix, which_norm, norm_vector=norm_vector_a,&
         error=error)
    CALL dbcsr_data_clear_pointer (norm_vector_a)
    CALL dbcsr_data_release (norm_vector_a)
  END SUBROUTINE dbcsr_norm_r4_vec

! *****************************************************************************
!> \brief compute a norm of a dbcsr matrix
!> \param[in] matrix  the matrix
!> \retval norm ...
! *****************************************************************************
  FUNCTION dbcsr_gershgorin_norm_r8(matrix) RESULT (norm)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(KIND=real_8)                        :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_gershgorin_norm_r8', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=real_4), DIMENSION(:, :), &
      POINTER                                :: data_c
    COMPLEX(KIND=real_8), DIMENSION(:, :), &
      POINTER                                :: data_z
    INTEGER                                  :: blk, col, col_offset, i, j, &
                                                nc, nr, row, row_offset
    LOGICAL                                  :: any_sym, tr
    REAL(KIND=real_4), DIMENSION(:, :), &
      POINTER                                :: data_r
    REAL(KIND=real_8), DIMENSION(:, :), &
      POINTER                                :: data_d
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: buff_d
    TYPE(dbcsr_error_type)                   :: error
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    nr = dbcsr_nfullrows_total(matrix)
    nc = dbcsr_nfullcols_total(matrix)

    any_sym = dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_symmetric.OR.&
              dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_antisymmetric

    CALL dbcsr_assert (nr.EQ.nc, dbcsr_fatal_level, dbcsr_wrong_args_error, &
         routineN, "not a square matrix",__LINE__,error)

    norm = 0.0_dp
    ALLOCATE(buff_d(nr))
    buff_d = 0.0_dp
    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       SELECT CASE (dbcsr_get_data_type(matrix))
       CASE (dbcsr_type_real_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_r, tr, blk, &
               row_offset=row_offset, col_offset=col_offset)
          DO j=1,SIZE(data_r,2)
          DO i=1,SIZE(data_r,1)
             buff_d(row_offset+i-1) = buff_d(row_offset+i-1) + ABS(data_r(i,j))
             IF(any_sym.AND.row.NE.col) &
                  buff_d(col_offset+j-1) = buff_d(col_offset+j-1) + ABS(data_r(i,j))
          ENDDO
          ENDDO
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk, &
               row_offset=row_offset, col_offset=col_offset)
          DO j=1,SIZE(data_d,2)
          DO i=1,SIZE(data_d,1)
             buff_d(row_offset+i-1) = buff_d(row_offset+i-1) + ABS(data_d(i,j))
             IF(any_sym.AND.row.NE.col) &
                  buff_d(col_offset+j-1) = buff_d(col_offset+j-1) + ABS(data_d(i,j))
          ENDDO
          ENDDO
       CASE (dbcsr_type_complex_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_c, tr, blk, &
               row_offset=row_offset, col_offset=col_offset)
          DO j=1,SIZE(data_c,2)
          DO i=1,SIZE(data_c,1)
             buff_d(row_offset+i-1) = buff_d(row_offset+i-1) + ABS(data_c(i,j))
             IF(any_sym.AND.row.NE.col) &
                  CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
                  routineN,"Only nonsymmetric matrix so far",__LINE__,error)
             !     buff_d(col_offset+j-1) = buff_d(col_offset+j-1) + ABS(data_c(i,j))
          ENDDO
          ENDDO
       CASE (dbcsr_type_complex_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_z, tr, blk, &
               row_offset=row_offset, col_offset=col_offset)
          DO j=1,SIZE(data_z,2)
          DO i=1,SIZE(data_z,1)
             buff_d(row_offset+i-1) = buff_d(row_offset+i-1) + ABS(data_z(i,j))
             IF(any_sym.AND.row.NE.col) &
                  CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
                  routineN,"Only nonsymmetric matrix so far",__LINE__,error)
             !     buff_d(col_offset+j-1) = buff_d(col_offset+j-1) + ABS(data_z(i,j))
          ENDDO
          ENDDO
       CASE DEFAULT
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Wrong data type",__LINE__,error)
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL mp_sum(buff_d,dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))
    norm = MAXVAL(buff_d)
    DEALLOCATE(buff_d)

  END FUNCTION dbcsr_gershgorin_norm_r8

! *****************************************************************************
!> \brief compute a norm of a dbcsr matrix
!> \param[in] matrix  the matrix
!> \retval norm ...
! *****************************************************************************
  FUNCTION dbcsr_maxabs_r8(matrix) RESULT (norm)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    REAL(real_8)                             :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_maxabs_r8', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=real_4), DIMENSION(:, :), &
      POINTER                                :: data_c
    COMPLEX(KIND=real_8), DIMENSION(:, :), &
      POINTER                                :: data_z
    INTEGER                                  :: blk, col, row
    LOGICAL                                  :: tr
    REAL(KIND=real_4), DIMENSION(:, :), &
      POINTER                                :: data_r
    REAL(KIND=real_8), DIMENSION(:, :), &
      POINTER                                :: data_d
    TYPE(dbcsr_error_type)                   :: error
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    norm = 0.0_dp
    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       SELECT CASE (dbcsr_get_data_type(matrix))
       CASE (dbcsr_type_real_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_r, tr, blk)
          norm = MAX(norm,REAL(MAXVAL(ABS(data_r)),dp))
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
          norm = MAX(norm,MAXVAL(ABS(data_d)))
       CASE (dbcsr_type_complex_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_c, tr, blk)
          norm = MAX(norm,REAL(MAXVAL(ABS(data_c)),dp))
       CASE (dbcsr_type_complex_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_z, tr, blk)
          norm = MAX(norm,MAXVAL(ABS(data_z)))
       CASE DEFAULT
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Wrong data type",__LINE__,error)
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)

    CALL mp_max(norm, dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))

  END FUNCTION dbcsr_maxabs_r8

! *****************************************************************************
!> \brief compute a norm of a dbcsr matrix
!> \param[in] matrix  the matrix
!> \param local ...
!> \retval norm ...
! *****************************************************************************
  FUNCTION dbcsr_frobenius_norm_r8(matrix, local) RESULT (norm)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    LOGICAL, INTENT(in), OPTIONAL            :: local
    REAL(KIND=real_8)                        :: norm

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_frobenius_norm_r8', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=real_4), DIMENSION(:, :), &
      POINTER                                :: data_c
    COMPLEX(KIND=real_8), DIMENSION(:, :), &
      POINTER                                :: data_z
    INTEGER                                  :: blk, col, error_handler, row
    LOGICAL                                  :: any_sym, my_local, tr
    REAL(KIND=real_4), DIMENSION(:, :), &
      POINTER                                :: data_r
    REAL(KIND=real_8), DIMENSION(:, :), &
      POINTER                                :: data_d
    REAL(real_8)                             :: fac
    TYPE(dbcsr_error_type)                   :: error
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    my_local = .FALSE.
    IF(PRESENT(local)) my_local = local

    any_sym = dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_symmetric.OR.&
              dbcsr_get_matrix_type(matrix).EQ.dbcsr_type_antisymmetric

    norm = 0.0_dp
    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       SELECT CASE (dbcsr_get_data_type(matrix))
       CASE (dbcsr_type_real_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_r, tr, blk)
          fac = 1.0_dp
          IF(any_sym.AND.row.NE.col) fac = 2.0_dp
          norm = norm + fac * SUM(data_r**2)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
          fac = 1.0_dp
          IF(any_sym.AND.row.NE.col) fac = 2.0_dp
          norm = norm + fac * SUM(data_d**2)
       CASE (dbcsr_type_complex_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_c, tr, blk)
          fac = 1.0_dp
          IF(any_sym.AND.row.NE.col) &
               CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Only nonsymmetric matrix so far",__LINE__,error)
          norm = norm + fac * REAL( SUM(CONJG(data_c)*data_c), KIND=real_8)
       CASE (dbcsr_type_complex_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_z, tr, blk)
          fac = 1.0_dp
          IF(any_sym.AND.row.NE.col) &
               CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Only nonsymmetric matrix so far",__LINE__,error)
          norm = norm + fac * REAL(SUM(CONJG(data_z)*data_z),KIND=real_8)
       CASE DEFAULT
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Wrong data type",__LINE__,error)
       END SELECT
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    IF(.NOT.my_local) CALL mp_sum(norm,dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))
    norm = SQRT(norm)

    CALL dbcsr_error_stop(error_handler, error)

  END FUNCTION dbcsr_frobenius_norm_r8


! *****************************************************************************
!> \brief Sums blocks in a replicated dbcsr matrix, which has the same structure on all ranks.
!> \param[in,out] matrix      dbcsr matrix to operate on
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_sum_replicated (matrix, error)
    TYPE(dbcsr_obj), INTENT(inout)           :: matrix
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sum_replicated', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: comm, error_handler, &
                                                index_checksum, mynode, &
                                                numnodes
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: all_checksums
    TYPE(dbcsr_mp_obj)                       :: mp

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    mp = dbcsr_distribution_mp (dbcsr_distribution (matrix))
    comm = dbcsr_mp_group (mp)
    numnodes = dbcsr_mp_numnodes (mp)
    mynode = dbcsr_mp_mynode (mp)
    !
    ALLOCATE(all_checksums(numnodes))
    CALL dbcsr_index_checksum (matrix, index_checksum, error)
    CALL mp_allgather (index_checksum, all_checksums, comm)
    !
    CALL dbcsr_assert (ALL (all_checksums .EQ. index_checksum),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Replicated matrices do not all have the same index structure.",&
         __LINE__, error=error)
    !
    SELECT CASE (dbcsr_data_get_type(matrix%m%data_area))
    CASE (dbcsr_type_real_4)
       CALL mp_sum (matrix%m%data_area%d%r_sp, comm)
    CASE (dbcsr_type_real_8)
       CALL mp_sum (matrix%m%data_area%d%r_dp, comm)
    CASE (dbcsr_type_complex_4)
       CALL mp_sum (matrix%m%data_area%d%c_sp, comm)
    CASE (dbcsr_type_complex_8)
       CALL mp_sum (matrix%m%data_area%d%c_dp, comm)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Incorrect data type", __LINE__, error=error)
    END SELECT
    !
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_sum_replicated


! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param matrix_b ...
!> \param trace ...
!> \param trans_a ...
!> \param trans_b ...
!> \param local_sum ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_trace_a_b_d(matrix_a, matrix_b, trace, trans_a, trans_b, local_sum, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a, matrix_b
    REAL(kind=real_8), INTENT(INOUT)         :: trace
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: trans_a, trans_b
    LOGICAL, INTENT(IN), OPTIONAL            :: local_sum
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_a_b_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler
    REAL(kind=real_4)                        :: trace_4

    CALL dbcsr_error_set(routineN, error_handler, error)
    IF(    dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_8 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_8 .OR. &
           dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_8 .OR. &
           dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_8 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_trace_ab_d(matrix_a, matrix_b, trace, trans_a, trans_b, local_sum, error)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4 .AND. &
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_4) THEN
       trace_4 = 0.0_real_4
       CALL dbcsr_trace_ab_s(matrix_a, matrix_b, trace_4, trans_a, trans_b, local_sum, error)
       trace = REAL(trace_4,real_8)
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid combination of data type, NYI",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_trace_a_b_d



! *****************************************************************************
!> \brief ...
!> \param matrix_a ...
!> \param trace ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_trace_a_any(matrix_a, trace, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_scalar_type), INTENT(INOUT)   :: trace
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_trace_a_any', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handle, error)
    SELECT CASE (dbcsr_scalar_get_type (trace))
    CASE (dbcsr_type_real_4)
       CALL dbcsr_trace (matrix_a, trace%r_sp, error)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_trace (matrix_a, trace%r_dp, error)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_trace (matrix_a, trace%c_sp, error)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_trace (matrix_a, trace%c_dp, error)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_caller_error,&
            routineN, "Invalid data type.",__LINE__,error)
    END SELECT

    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_trace_a_any


! *****************************************************************************
!> \brief check if a block is not in the limits
!> \param[in] row ...
!> \param[in] col ...
!> \param[in] block_row_limits ...
!> \param[in] block_column_limits ...
!> \retval dbcsr_block_in_limits ...
! *****************************************************************************
  FUNCTION dbcsr_block_in_limits(row, col, block_row_limits, block_column_limits)
    INTEGER, INTENT(in)                      :: row, col
    INTEGER, DIMENSION(2), INTENT(in), &
      OPTIONAL                               :: block_row_limits, &
                                                block_column_limits
    LOGICAL                                  :: dbcsr_block_in_limits

    dbcsr_block_in_limits = .TRUE.
    IF (PRESENT (block_row_limits)) THEN
       IF (row .LT. block_row_limits(1)) dbcsr_block_in_limits = .FALSE.
       IF (row .GT. block_row_limits(2)) dbcsr_block_in_limits = .FALSE.
    ENDIF
    IF (PRESENT (block_column_limits)) THEN
       IF (col .LT. block_column_limits(1)) dbcsr_block_in_limits = .FALSE.
       IF (col .GT. block_column_limits(2)) dbcsr_block_in_limits = .FALSE.
    ENDIF
  END FUNCTION dbcsr_block_in_limits


! *****************************************************************************
!> \brief Gets information about a matrix
!> \param[in] matrix          matrix to query
!> \param[out] nblkrows_total ...
!> \param[out] nblkcols_total ...
!> \param[out] nfullrows_total ...
!> \param[out] nfullcols_total ...
!> \param[out] nblkrows_local ...
!> \param[out] nblkcols_local ...
!> \param[out] nfullrows_local ...
!> \param[out] nfullcols_local ...
!> \param[out] my_prow ...
!> \param[out] my_pcol ...
!> \param[out] local_rows ...
!> \param[out] local_cols ...
!> \param[out] proc_row_dist ...
!> \param[out] proc_col_dist ...
!> \param[out] row_blk_size ...
!> \param[out] col_blk_size ...
!> \param[out] row_blk_offset ...
!> \param[out] col_blk_offset ...
!> \param[out] distribution   the data distribution of the matrix
!> \param[out] name           matrix name
!> \param[out] data_area      data_area
!> \param[out] matrix_type    matrix type (regular, symmetric, see
!>                            dbcsr_types.F for values)
!> \param[out] data_type      data type (single/double precision real/complex)
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
!> \param occ_thresh ...
!> \retval may_be_dense    use the mutable and not append-only working
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
!> \brief Returns the occupation of the matrix
!> \param matrix     matrix from which to get the occupation
!> \retval occupation ...
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


#include "dbcsr_operations_d.f90"
#include "dbcsr_operations_z.f90"
#include "dbcsr_operations_s.f90"
#include "dbcsr_operations_c.f90"

END MODULE dbcsr_operations

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
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
                                             array_equality,&
                                             array_hold,&
                                             array_i1d_obj,&
                                             array_release
  USE dbcsr_blas_operations,           ONLY: dbcsr_blas_copy,&
                                             dbcsr_blas_gemm,&
                                             dbcsr_blas_scal
  USE dbcsr_block_access,              ONLY: dbcsr_get_block_p,&
                                             dbcsr_put_block,&
                                             dbcsr_remove_block,&
                                             dbcsr_reserve_blocks
  USE dbcsr_block_operations,          ONLY: block_add_on_diag,&
                                             dbcsr_block_partial_copy,&
                                             dbcsr_block_transpose,&
                                             dbcsr_data_clear,&
                                             get_block2d_diagonal,&
                                             set_block2d_diagonal
  USE dbcsr_data_operations,           ONLY: dbcsr_data_convert,&
                                             dbcsr_data_copyall
  USE dbcsr_dist_operations,           ONLY: checker_square_proc,&
                                             create_bl_distribution,&
                                             dbcsr_create_image_dist,&
                                             dbcsr_find_column,&
                                             dbcsr_get_stored_coordinates
  USE dbcsr_error_handling
  USE dbcsr_index_operations,          ONLY: dbcsr_expand_row_index,&
                                             dbcsr_make_index_canonical,&
                                             dbcsr_repoint_index
  USE dbcsr_internal_operations,       ONLY: dbcsr_insert_blocks,&
                                             dbcsr_mult_m_e_e
  USE dbcsr_io,                        ONLY: dbcsr_print
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_kinds,                     ONLY: dp,&
                                             int_8,&
                                             real_4,&
                                             real_8,&
                                             sp
  USE dbcsr_message_passing,           ONLY: dmp_max,&
                                             mp_allgather,&
                                             mp_bcast,&
                                             mp_recv,&
                                             mp_send,&
                                             mp_sum
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_offsets, dbcsr_col_block_sizes, &
       dbcsr_data_clear_pointer, dbcsr_data_ensure_size, dbcsr_data_hold, &
       dbcsr_data_init, dbcsr_data_new, dbcsr_data_release, &
       dbcsr_data_set_pointer, dbcsr_destroy, dbcsr_destroy_array, &
       dbcsr_destroy_image_dist, dbcsr_distribution, &
       dbcsr_distribution_col_dist, dbcsr_distribution_hold, &
       dbcsr_distribution_make_threads, dbcsr_distribution_mp, &
       dbcsr_distribution_new, dbcsr_distribution_no_threads, &
       dbcsr_distribution_release, dbcsr_distribution_row_dist, &
       dbcsr_get_data, dbcsr_get_data_p, dbcsr_get_data_size, &
       dbcsr_get_data_type, dbcsr_get_info, dbcsr_get_matrix_type, &
       dbcsr_get_num_blocks, dbcsr_get_replication_type, dbcsr_has_symmetry, &
       dbcsr_init, dbcsr_is_initialized, dbcsr_max_col_size, &
       dbcsr_max_row_size, dbcsr_may_be_dense, dbcsr_mp_grid_setup, &
       dbcsr_mp_group, dbcsr_mp_has_subgroups, dbcsr_mp_my_col_group, &
       dbcsr_mp_my_row_group, dbcsr_mp_mynode, dbcsr_mp_mypcol, &
       dbcsr_mp_myprow, dbcsr_mp_npcols, dbcsr_mp_nprows, dbcsr_mp_numnodes, &
       dbcsr_mp_pgrid, dbcsr_name, dbcsr_nblkcols_total, &
       dbcsr_nblkrows_total, dbcsr_nfullcols_total, dbcsr_nfullrows_total, &
       dbcsr_release, dbcsr_row_block_sizes, dbcsr_scalar, &
       dbcsr_scalar_are_equal, dbcsr_scalar_one, dbcsr_scalar_zero, &
       dbcsr_type_1d_to_2d, dbcsr_valid_index
  USE dbcsr_mp_operations,             ONLY: dbcsr_recv_any,&
                                             dbcsr_send_any
  USE dbcsr_ptr_util,                  ONLY: ensure_array_size,&
                                             pointer_view
  USE dbcsr_transformations,           ONLY: dbcsr_make_dense,&
                                             dbcsr_make_images,&
                                             dbcsr_make_images_dense,&
                                             dbcsr_make_undense,&
                                             dbcsr_new_transposed
  USE dbcsr_types,                     ONLY: &
       dbcsr_2d_array_type, dbcsr_data_obj, dbcsr_distribution_obj, &
       dbcsr_filter_frobenius, dbcsr_imagedistribution_type, dbcsr_iterator, &
       dbcsr_mp_obj, dbcsr_norm_column, dbcsr_norm_frobenius, &
       dbcsr_norm_gershgorin, dbcsr_norm_maxabsnorm, dbcsr_obj, &
       dbcsr_repl_col, dbcsr_repl_full, dbcsr_repl_none, dbcsr_repl_row, &
       dbcsr_scalar_type, dbcsr_type, dbcsr_type_antisymmetric, &
       dbcsr_type_complex_4, dbcsr_type_complex_8, dbcsr_type_normal, &
       dbcsr_type_real_4, dbcsr_type_real_8, dbcsr_type_symmetric
  USE dbcsr_util,                      ONLY: ceil_log2,&
                                             convert_sizes_to_offsets,&
                                             dbcsr_checksum,&
                                             dbcsr_set_debug,&
                                             dbcsr_verify_matrix,&
                                             swap
  USE dbcsr_work_operations,           ONLY: dbcsr_add_wm_from_matrix,&
                                             dbcsr_create,&
                                             dbcsr_finalize,&
                                             dbcsr_work_create

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_operations'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034


  PUBLIC :: dbcsr_multiply,&
            dbcsr_mult_m_e_e, dbcsr_trace, dbcsr_add_on_diag,&
            dbcsr_set, dbcsr_scale, dbcsr_add, dbcsr_copy,&
            dbcsr_copy_submatrix, &
            dbcsr_get_diag, dbcsr_set_diag, &
            dbcsr_get_block_diag, dbcsr_hadamard_product, &
            dbcsr_filter, dbcsr_scale_by_vector,&
            dbcsr_replace_blocks, dbcsr_conjg,&
            dbcsr_btriu, dbcsr_triu, dbcsr_tril,&
            dbcsr_symmetrize_block_diag, dbcsr_copy_columns,&
            dbcsr_init_random, dbcsr_lanczos
  PUBLIC :: dbcsr_sum_replicated
  PUBLIC :: dbcsr_norm, &
            dbcsr_gershgorin_norm, dbcsr_maxabs, dbcsr_frobenius_norm

! The interfaces for the generic routines found in the generated
! generic files.

  INTERFACE dbcsr_multiply
     MODULE PROCEDURE dbcsr_multiply_anytype
     MODULE PROCEDURE dbcsr_multiply_s, dbcsr_multiply_d,&
                      dbcsr_multiply_c, dbcsr_multiply_z
  END INTERFACE

  INTERFACE dbcsr_trace
     MODULE PROCEDURE dbcsr_trace_a_s, dbcsr_trace_a_d,&
          dbcsr_trace_a_c, dbcsr_trace_a_z,&
          dbcsr_trace_ab_s, dbcsr_trace_a_b_d,&!dbcsr_trace_ab_d,&
          dbcsr_trace_ab_c, dbcsr_trace_ab_z
  END INTERFACE

  INTERFACE dbcsr_scale
     MODULE PROCEDURE dbcsr_scale_anytype
     MODULE PROCEDURE dbcsr_scale_s, dbcsr_scale_d,&
                      dbcsr_scale_c, dbcsr_scale_z
     MODULE PROCEDURE dbcsr_scale_s_m, dbcsr_scale_d_m,&
                      dbcsr_scale_c_m, dbcsr_scale_z_m
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

  INTERFACE dbcsr_lanczos
     MODULE PROCEDURE dbcsr_lanczos_extremal_eig
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


#define temp_transpose(v, r, c) RESHAPE(TRANSPOSE(RESHAPE(v,(/r,c/))),(/r*c/))

  INTEGER, PARAMETER, PRIVATE :: rpslot_owner = 1
  INTEGER, PARAMETER, PRIVATE :: rpslot_addblks = 2
  INTEGER, PARAMETER, PRIVATE :: rpslot_addoffset = 3
  INTEGER, PARAMETER, PRIVATE :: rpslot_oldblks = 4
  INTEGER, PARAMETER, PRIVATE :: rpslot_oldoffset = 5
  INTEGER, PARAMETER, PRIVATE :: rpslot_totaloffset = 6
  INTEGER, PARAMETER, PRIVATE :: rpnslots = 6


  LOGICAL, PARAMETER, PRIVATE :: detailed_timing = .FALSE.

  INTEGER, PARAMETER :: dgemm_stack_size = 32

  TYPE block_parameters
     LOGICAL :: tr
     INTEGER :: logical_rows, logical_cols
     INTEGER :: offset, nze
  END TYPE block_parameters

  TYPE dgemm_join
     INTEGER :: p_a, p_b, p_c
     INTEGER :: last_k, last_n
     TYPE(dbcsr_scalar_type) :: alpha, beta
  END TYPE dgemm_join

CONTAINS


! *****************************************************************************
!> \brief Performs a multiplication of two dbcsr_type matrices,
!>        as  C := alpha * op( A ) * op( B ) + beta * C.
!> \par Matrices m_a and m_b are multiplied into the m_c product matrix. If the
!>      dist2d parameter is not specified, then a new distribution_2d is
!>      determined for it.
!> \par Non-equal column dimensions of the right and product matrices
!>      The right and product matrix are allowed to have different
!>      (full) column dimensions. If they differ, there are certain
!>      peculiar behaviors, then the last_column is effectively set to
!>      the minimal of the two.
!> \par Beta scaling of the right product matrix
!>      If the effective last_column is less than the full column
!>      dimension of the product matrix, then the scaling of the
!>      product matrix with beta is limited to the submatrix specified
!>      by last_column.
!> \par Sets
!>      The left and right sets may be passed in (if the left_set or
!>      right_set argument is present and non-NULL), or they may be
!>      requeseted to be saved (present and NULL). No verification is done
!>      that the set is consistent with the given matrix.
!> \param[in] transa, transb  specifies the form of op( A or B ) to be used in 
!>                            the matrix multiplication 
!>                            transx = 'N' or 'n',  op( X ) = X.
!>                            transx = 'T' or 't',  op( X ) = transpose(X).
!>                            transx = 'C' or 'c',  op( X ) = transpose(conjg(X)).
!> \param[in] alpha           scaling of product
!> \param[in] matrix_a        left BCSR matrix
!> \param[in] matrix_b        right BCSR matrix
!> \param[in] beta            scaling of existing data
!> \param[out] matrix_c       resulting BCSR product matrix.
!> \param[in] first_row       (optional) first full row of limiting submatrix
!> \param[in] last_row        (optional) last full row of limiting submatrix
!> \param[in] first_column    (optional) first full column of limiting submatrix
!> \param[in] last_column     (optional) last full column of limiting submatrix
!> \param[in] first_k         (optional) first full column of limiting inner
!>                            product
!> \param[in] last_k          (optional) last full column of limiting inner
!>                            product
!> \param[in] retain_sparsity (optional) enforce the sparsity pattern of the
!>                            existing product matrix; default is no
!> \param[in,out] left_set    (optional) multiplication set for left matrix
!> \param[in,out] right_set   (optional) multiplication set for right matrix
!> \param[in,out] error       error
!> \param[out] flop           (optional) effective flop
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_anytype(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, left_set, right_set, error, flop)
    
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_2d_array_type), OPTIONAL, &
      POINTER                                :: left_set, right_set
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_anytype', &
      routineP = moduleN//':'//routineN
    REAL(real_8), PARAMETER                  :: make_dense_occ_thresh = 1.0_dp

    INTEGER                                  :: error_handler, f_col, f_k, &
                                                f_row, i, l_col, l_k, l_row
    INTEGER(KIND=int_8)                      :: my_flop
    LOGICAL :: ab_dense, keep_product_data, keep_sparsity, left_set_given, &
      new_left, new_right, product_reindex, right_set_given, use_dense_mult
    REAL(KIND=dp)                            :: cs
    TYPE(array_i1d_obj) :: old_product_col_blk_offsets, &
      old_product_col_blk_sizes, old_product_row_blk_offsets, &
      old_product_row_blk_sizes
    TYPE(dbcsr_2d_array_type), POINTER       :: m2s_left, m2s_right
    TYPE(dbcsr_distribution_obj)             :: old_product_distribution
    TYPE(dbcsr_imagedistribution_type)       :: rdist_left, rdist_right
    TYPE(dbcsr_obj)                          :: matrix_left, matrix_right, &
                                                product_matrix

! this will become the default

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL dbcsr_error_set(routineN, error_handler, error)
    !
    ! check parameters
    CALL dbcsr_assert(transa.EQ.'N'.OR.transa.EQ.'T'.OR.transa.EQ.'C'.OR.&
                   transa.EQ.'n'.OR.transa.EQ.'t'.OR.transa.EQ.'c',&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Invalid transa = "//transa, __LINE__, error)
    CALL dbcsr_assert(transb.EQ.'N'.OR.transb.EQ.'T'.OR.transb.EQ.'C'.OR.&
                   transb.EQ.'n'.OR.transb.EQ.'t'.OR.transb.EQ.'c',&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Invalid transb = "//transb, __LINE__, error)

    IF (dbg) THEN
       WRITE(*,*)'========== MULTIPLICATION ========================'
       CALL dbcsr_verify_matrix (matrix_a, error=error)
       CALL dbcsr_verify_matrix (matrix_b, error=error)
       CALL dbcsr_verify_matrix (matrix_c, error=error)
       WRITE(*,*)routineN//" ABC checksums",&
            dbcsr_checksum(matrix_a, error=error),&
            dbcsr_checksum(matrix_b, error=error),&
            dbcsr_checksum(matrix_c, error=error)
       CALL dbcsr_print (matrix_a, nodata=.TRUE., error=error)
       CALL dbcsr_print (matrix_b, nodata=.TRUE., error=error)
       CALL dbcsr_print (matrix_c, nodata=.TRUE., error=error)
    ENDIF
    !
    ! transpose/conjg left and/or right matrices if needed 
    SELECT CASE(transa)
    CASE('N','n')
       matrix_left = matrix_a
       new_left = .FALSE.
    CASE('T','t')
       CALL dbcsr_init(matrix_left)
       CALL dbcsr_new_transposed (matrix_left, matrix_a,&
            shallow_data_copy=.TRUE.,&
            transpose_distribution=.FALSE., error=error)
       new_left = .TRUE.
    CASE('C','c')
       CALL dbcsr_init(matrix_left)
       CALL dbcsr_new_transposed (matrix_left, matrix_a,&
            shallow_data_copy=.FALSE.,&
            transpose_distribution=.FALSE., error=error)
       CALL dbcsr_conjg(matrix_left, error=error)
       new_left = .TRUE.
    CASE DEFAULT
       CALL dbcsr_assert(.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, &
            routineN, "wrong transa = "//transa, __LINE__, error)
    END SELECT

    SELECT CASE(transb)
    CASE('N','n')
       matrix_right = matrix_b
       new_right = .FALSE.
    CASE('T','t')
       CALL dbcsr_init(matrix_right)
       CALL dbcsr_new_transposed (matrix_right, matrix_b,&
            shallow_data_copy=.TRUE.,&
            transpose_distribution=.FALSE., error=error)
       new_right = .TRUE.
    CASE('C','c')
       CALL dbcsr_init(matrix_right)
       CALL dbcsr_new_transposed (matrix_right, matrix_b,&
            shallow_data_copy=.FALSE.,&
            transpose_distribution=.FALSE., error=error)
       CALL dbcsr_conjg(matrix_right, error=error)
       new_right = .TRUE.
    CASE DEFAULT
       CALL dbcsr_assert(.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, &
            routineN, "wrong transb = "//transb, __LINE__, error)
    END SELECT
    !
    ! we skip dense multiply for (anti)symmetric matrices (slowdown for S/H * C)
    IF(dbcsr_has_symmetry (matrix_left) .OR. &
         dbcsr_has_symmetry(matrix_right)) THEN
       use_dense_mult = .FALSE.
    ELSE
       use_dense_mult = dbcsr_may_be_dense (matrix_left, make_dense_occ_thresh)&
            .AND. dbcsr_may_be_dense (matrix_right, make_dense_occ_thresh)
    ENDIF
    ab_dense = use_dense_mult
    ! Submatrix selection
    f_row = 1
    l_row = dbcsr_nfullrows_total(matrix_c)
    f_col = 1
    l_col = dbcsr_nfullcols_total(matrix_c)
    f_k = 0
    l_k = 0
    IF (PRESENT (first_row)) THEN
       CALL dbcsr_assert(first_row .GE. 1&
            .AND. first_row .LE. dbcsr_nfullrows_total(matrix_c),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid first row specified", __LINE__, error)
       f_row = first_row
    ENDIF
    IF (PRESENT (last_row)) THEN
       CALL dbcsr_assert(last_row .GE. 1&
            .AND. last_row .LE. dbcsr_nfullrows_total(matrix_c)&
            .OR. last_row .LT. 1 ,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid last row specified", __LINE__, error)
       l_row = last_row
    ENDIF
    IF (PRESENT (first_column)) THEN
       CALL dbcsr_assert(first_column .GE. 1&
            .AND. first_column .LE. dbcsr_nfullcols_total(matrix_c),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid first col specified", __LINE__, error)
       f_col = first_column
    ENDIF
    IF (PRESENT (last_column)) THEN
       CALL dbcsr_assert(last_column .GE. 1&
            .AND. last_column .LE. dbcsr_nfullcols_total(matrix_c)&
            .OR. last_column .LT. 1,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid last column specified (C)", __LINE__, error)
       CALL dbcsr_assert(last_column .GE. 1&
            .AND. last_column .LE. dbcsr_nfullcols_total(matrix_right)&
            .OR. last_column .LT. 1,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid last column specified (B)", __LINE__, error)
       l_col = last_column
    ENDIF
    IF (PRESENT (first_k)) THEN
       CALL dbcsr_assert(first_k .GE. 1&
            .AND. first_k .LE. dbcsr_nfullcols_total(matrix_left),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid first k specified (A)", __LINE__, error)
       f_k = first_k
    ENDIF
    IF (PRESENT (last_k)) THEN
       CALL dbcsr_assert(last_k.GE. 1&
            .AND. last_k .LE. dbcsr_nfullcols_total(matrix_left)&
            .OR. last_k  .LT. 1,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid last k specified (A)", __LINE__, error)
       l_k = last_k
    ENDIF
    !
    ! Now optimize the default submatrix selection values away
    IF (f_row .EQ. 1) f_row = 0
    IF (l_row .EQ. dbcsr_nfullrows_total(matrix_left)) l_row = 0
    IF (f_col .EQ. 1) f_col = 0
    ! The last column must be set if the right and product matrices
    ! differ.
    l_col = MIN (l_col, dbcsr_nfullcols_total(matrix_right))
    l_col = MIN (l_col, dbcsr_nfullcols_total(matrix_c))
    IF (f_col.LE.1.AND.&
        l_col .EQ. dbcsr_nfullcols_total(matrix_right) .AND. &
        dbcsr_nfullcols_total(matrix_right) .EQ.&
        dbcsr_nfullcols_total(matrix_c)) l_col = 0
    IF (f_k .EQ. 1) f_k = 0
    IF (l_k .EQ. dbcsr_nfullcols_total(matrix_left)) l_k = 0
    IF (.NOT. PRESENT(last_column) .AND.&
        .NOT. array_equality (dbcsr_col_block_sizes (matrix_right),&
                              dbcsr_col_block_sizes (matrix_c))) THEN
       l_col = MIN (dbcsr_nfullcols_total(matrix_right),&
                    dbcsr_nfullcols_total(matrix_c))
    ENDIF
    CALL dbcsr_assert (f_row .LE. l_row, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Last row smaller than first row", &
         __LINE__, error)
    CALL dbcsr_assert (f_col .LE. l_col, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Last col smaller than first col", &
         __LINE__, error)
    !
    !
    ! Product data needs to be retained when
    ! * beta != 0; or
    ! * there is column limiting (l_col > 0) and the limiting column
    !   is less than the number of full columns in theproduct matrix
    keep_sparsity = .FALSE.
    IF (PRESENT (retain_sparsity)) keep_sparsity=retain_sparsity
    !
    keep_product_data = keep_sparsity&
         .OR. .NOT. dbcsr_scalar_are_equal (beta, dbcsr_scalar_zero(beta%data_type))&
         .OR. (l_col .GT. 0 .AND. l_col .LT. dbcsr_nfullcols_total(matrix_c))
    !
    IF (.NOT. dbcsr_scalar_are_equal (beta, dbcsr_scalar_one(beta%data_type)) .AND. keep_product_data) THEN
       CALL dbcsr_scale (matrix_c, scale_type='S',&
            alpha_scalar=beta, &
            last_column=last_column, error=error)
    ENDIF
    !
    ! The index of the product matrix is twiddled into canonical form
    ! if it is (anti)symmetric (i.e., rows and columns are where the
    ! row/column distributions say they are). Doing this in advance
    ! makes the local multiply more efficient.
    IF (dbcsr_has_symmetry (matrix_c)) THEN
       product_reindex = .TRUE.
    ELSE
       product_reindex = .FALSE.
    ENDIF
    ! Product can not be made dense; however, A & B may still be made
    ! dense unless previously determined otherwise.
    IF (product_reindex) THEN
       use_dense_mult = .FALSE.
    ENDIF
    product_matrix = matrix_c
    IF (use_dense_mult) THEN
       old_product_row_blk_offsets = product_matrix%m%row_blk_offset
       old_product_col_blk_offsets = product_matrix%m%col_blk_offset
       old_product_row_blk_sizes = product_matrix%m%row_blk_size
       old_product_col_blk_sizes = product_matrix%m%col_blk_size
       CALL array_hold (old_product_row_blk_offsets)
       CALL array_hold (old_product_col_blk_offsets)
       CALL array_hold (old_product_row_blk_sizes)
       CALL array_hold (old_product_col_blk_sizes)
       old_product_distribution = dbcsr_distribution (product_matrix)
       CALL dbcsr_distribution_hold (old_product_distribution)
       CALL dbcsr_make_dense (product_matrix, error=error)
    ENDIF
!$omp parallel default (none) shared (product_matrix, error, i)
    ! The thread distribution must reflect the current (possibly
    ! dense) distribution
    CALL dbcsr_distribution_make_threads (product_matrix%m%dist)
!$omp single
    i = 1
!$  i = omp_get_num_threads()
!$omp end single nowait
!$omp end parallel
    IF (dbg) THEN
       WRITE(*,*)routineN//" Matrices ", dbcsr_get_matrix_type(matrix_a),&
            dbcsr_get_matrix_type(matrix_b), dbcsr_get_matrix_type(matrix_c)
       WRITE(*,*)routineN//" Matrices ", transA, transB, "keep", keep_product_data
    ENDIF
    IF (keep_product_data) THEN
       IF (product_reindex) THEN
          IF (dbg) WRITE(*,*)routineN//" Making canonical index"
          CALL dbcsr_make_index_canonical (product_matrix)
       ENDIF
       CALL dbcsr_assert (.NOT. ASSOCIATED (product_matrix%m%wms),&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Product matrix should be finalized!", __LINE__, error=error)
!$omp parallel default (none) shared (product_matrix, i, error)
       ! For the multiply logic to work correctly, existing data must
       ! be added only after the index has been transformed into the
       ! canonical form.
!@@@
       !vwWRITE(*,*)routineN//' nwms', ASSOCIATED (product_matrix%m%wms)
       !vwIF (ASSOCIATED (product_matrix%m%wms)) THEN
       !vw   WRITE(*,*)routineN//' nwms', SIZE(product_matrix%m%wms)
       !vwENDIF

!$omp barrier
       CALL dbcsr_add_wm_from_matrix(product_matrix, error=error)
!$omp end parallel
    ENDIF
    product_matrix%m%nze = 0
    product_matrix%m%nblks = 0
    product_matrix%m%valid = .FALSE.
    !
    ! Create imaged distributions for the multiply.
    !
    ! matrix_c is used instead of product_matrix so that the original
    ! blocked (non-dense) form is used.
    CALL dbcsr_create_image_dist (rdist_right, matrix_right%m%dist,&
         match_row_nbins = dbcsr_mp_npcols (dbcsr_distribution_mp (matrix_left%m%dist)),&
         match_col_nbins = dbcsr_mp_npcols (dbcsr_distribution_mp (matrix_c%m%dist)),&
         match_col_pdist = array_data (dbcsr_distribution_col_dist (matrix_c%m%dist)))
    CALL dbcsr_create_image_dist (rdist_left, matrix_left%m%dist,&
         match_row_pdist = array_data (dbcsr_distribution_row_dist (matrix_c%m%dist)),&
         match_row_nbins = dbcsr_mp_nprows (dbcsr_distribution_mp (matrix_c%m%dist)),&
         match_col_pdist = array_data (dbcsr_distribution_row_dist (rdist_right%main)),&
         match_col_idist = array_data (rdist_right%row_image),&
         match_col_nbins = dbcsr_mp_nprows (dbcsr_distribution_mp(matrix_right%m%dist)))
    NULLIFY (m2s_right)
    right_set_given = .FALSE.
    IF (PRESENT (right_set)) THEN
       IF (ASSOCIATED (right_set)) THEN
          right_set_given = .TRUE.
          m2s_right => right_set
       ENDIF
    ENDIF
    NULLIFY (m2s_left)
    left_set_given = .FALSE.
    IF (PRESENT (left_set)) THEN
       IF (ASSOCIATED (left_set)) THEN
          left_set_given = .TRUE.
          m2s_left => left_set
       ENDIF
    ENDIF
    IF (.NOT. right_set_given) THEN
       ALLOCATE (m2s_right)
       CALL dbcsr_make_images (matrix_right, m2s_right, rdist_right,&
            predistribute="R", special_memory=.FALSE.,&
            no_copy_data=use_dense_mult, error=error)
       IF (ab_dense) THEN
          CALL dbcsr_make_images_dense (m2s_right, "R", &
               join_cols = use_dense_mult, error=error)
          rdist_right = m2s_right%image_dist
       ENDIF
       IF (PRESENT (right_set)) THEN
          right_set => m2s_right
          right_set%image_dist = rdist_right
       ENDIF
    ELSE
       m2s_right => right_set
    ENDIF
    IF (.NOT. left_set_given) THEN
       ALLOCATE (m2s_left)
       CALL dbcsr_make_images (matrix_left, m2s_left, rdist_left,&
            predistribute="L", special_memory=.FALSE.,&
            no_copy_data=use_dense_mult, error=error)
       IF (ab_dense) THEN
          CALL dbcsr_make_images_dense (m2s_left, "L",&
               join_rows = use_dense_mult, error=error)
          rdist_left = m2s_left%image_dist
       ENDIF
       IF (PRESENT (left_set)) THEN
          left_set =>  m2s_left
          left_set%image_dist = rdist_left
       ENDIF
    ELSE
       m2s_left => left_set
    ENDIF
    IF (dbg) THEN
       DO i = 1, SIZE (m2s_right%mats,1)
          CALL dbcsr_verify_matrix(m2s_right%mats(i,1), error=error)
          CALL dbcsr_print(m2s_right%mats(i,1),nodata=.TRUE., error=error)
       ENDDO
       DO i = 1, SIZE (m2s_left%mats,2)
          CALL dbcsr_verify_matrix(m2s_left%mats(1,i), error=error)
          CALL dbcsr_print(m2s_left%mats(1,i),nodata=.TRUE., error=error)
       ENDDO
    ENDIF
    !
    IF (PRESENT (flop)) flop = 0
    my_flop = 0
    CALL dbcsr_mult_m_e_e (m2s_left, m2s_right, product_matrix,&
         (/ f_row, l_row, f_col, l_col, f_k, l_k /),&
         alpha=alpha, beta=beta, retain_sparsity=retain_sparsity,&
         error=error, flop=my_flop)
    IF (PRESENT (flop)) flop = my_flop
    !
    IF (new_left) CALL dbcsr_release (matrix_left)
    IF (new_right) CALL dbcsr_release (matrix_right)
    CALL dbcsr_distribution_no_threads (product_matrix%m%dist)
    ! The index of the product matrix is reset to the CP2K form if it
    ! was previously set to the canonical form.
    IF (product_reindex) THEN
       IF (dbg) WRITE(*,*)routineN//" Making CP2K index"
       CALL dbcsr_make_index_canonical(product_matrix, cp2k=.TRUE.)
    ELSEIF (use_dense_mult) THEN
       CALL dbcsr_make_undense(product_matrix, old_product_distribution,&
            old_product_row_blk_offsets, old_product_col_blk_offsets,&
            old_product_row_blk_sizes, old_product_col_blk_sizes, error=error)
       CALL array_release (old_product_row_blk_offsets)
       CALL array_release (old_product_col_blk_offsets)
       CALL array_release (old_product_row_blk_sizes)
       CALL array_release (old_product_col_blk_sizes)
       CALL dbcsr_distribution_release (old_product_distribution)
    ENDIF
    matrix_c = product_matrix
    IF (dbg) THEN
       CALL dbcsr_print (matrix_c, nodata=.TRUE., error=error)
    ENDIF
    IF (.NOT. PRESENT (left_set)) THEN
       CALL dbcsr_destroy_array (m2s_left, error=error)
       CALL dbcsr_destroy_image_dist (rdist_left)
       DEALLOCATE (m2s_left)
    ELSE
       IF (.NOT. left_set_given) THEN
          left_set%image_dist = m2s_left%image_dist
       ELSE
          CALL dbcsr_destroy_image_dist (rdist_left)
       ENDIF
    ENDIF
    IF (.NOT. PRESENT (right_set)) THEN
       CALL dbcsr_destroy_array (m2s_right, error=error)
       CALL dbcsr_destroy_image_dist (rdist_right)
       DEALLOCATE (m2s_right)
    ELSE
       IF (.not.right_set_given) THEN
          right_set%image_dist = m2s_right%image_dist
       ELSE
          CALL dbcsr_destroy_image_dist (rdist_right)
       ENDIF
    ENDIF
    IF (dbg) THEN
       WRITE(*,*)routineN//" Use dense mult, symm", use_dense_mult, ab_dense, &
            dbcsr_has_symmetry (matrix_c)
       CALL dbcsr_verify_matrix (matrix_c, error=error)
       CALL dbcsr_print (matrix_c, nodata=.TRUE., error=error)
       cs = dbcsr_checksum (matrix_c, error=error)
       WRITE(*,*)routineN//" Product checksum", cs
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_multiply_anytype
  


  SUBROUTINE dbcsr_multiply_s(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, left_set, right_set, error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    REAL(KIND=real_4), INTENT(IN)            :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    REAL(KIND=real_4), INTENT(IN)            :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_2d_array_type), OPTIONAL, &
      POINTER                                :: left_set, right_set
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CALL dbcsr_multiply_anytype(transa, transb,&
         dbcsr_scalar(alpha), matrix_a, matrix_b, dbcsr_scalar(beta), matrix_c,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, left_set=left_set, right_set=right_set, error=error, flop=flop)
  END SUBROUTINE dbcsr_multiply_s
  SUBROUTINE dbcsr_multiply_d(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, left_set, right_set, error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    REAL(KIND=real_8), INTENT(IN)            :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    REAL(KIND=real_8), INTENT(IN)            :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_2d_array_type), OPTIONAL, &
      POINTER                                :: left_set, right_set
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_d', &
      routineP = moduleN//':'//routineN

    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4 .AND.&
       dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_4 .AND.&
       dbcsr_get_data_type(matrix_c) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_multiply_anytype(transa, transb,&
            dbcsr_scalar(REAL(alpha,real_4)), matrix_a, matrix_b, &
            dbcsr_scalar(REAL(beta,real_4)), matrix_c,&
            first_row, last_row, first_column, last_column, first_k, last_k,&
            retain_sparsity, left_set=left_set, right_set=right_set, error=error, flop=flop)
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_8 .AND.&
           dbcsr_get_data_type(matrix_b) .EQ. dbcsr_type_real_8 .AND.&
           dbcsr_get_data_type(matrix_c) .EQ. dbcsr_type_real_8) THEN
       CALL dbcsr_multiply_anytype(transa, transb,&
            dbcsr_scalar(alpha), matrix_a, matrix_b, dbcsr_scalar(beta), matrix_c,&
            first_row, last_row, first_column, last_column, first_k, last_k,&
            retain_sparsity, left_set=left_set, right_set=right_set, error=error, flop=flop)
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_failure_level, dbcsr_internal_error,&
            routineP, "This combination of data types NYI",__LINE__, error)       
    ENDIF
  END SUBROUTINE dbcsr_multiply_d
  SUBROUTINE dbcsr_multiply_c(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, left_set, right_set, error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    COMPLEX(KIND=real_4), INTENT(IN)         :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    COMPLEX(KIND=real_4), INTENT(IN)         :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_2d_array_type), OPTIONAL, &
      POINTER                                :: left_set, right_set
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CALL dbcsr_multiply_anytype(transa, transb,&
         dbcsr_scalar(alpha), matrix_a, matrix_b, dbcsr_scalar(beta), matrix_c,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, left_set=left_set, right_set=right_set, error=error, flop=flop)
  END SUBROUTINE dbcsr_multiply_c
  SUBROUTINE dbcsr_multiply_z(transa, transb,&
       alpha, matrix_a, matrix_b, beta, matrix_c,&
       first_row, last_row, first_column, last_column, first_k, last_k,&
       retain_sparsity, left_set, right_set, error, flop)
    CHARACTER(LEN=1), INTENT(IN)             :: transa, transb
    COMPLEX(KIND=real_8), INTENT(IN)         :: alpha
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    COMPLEX(KIND=real_8), INTENT(IN)         :: beta
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_column, last_column, &
                                                first_k, last_k
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_sparsity
    TYPE(dbcsr_2d_array_type), OPTIONAL, &
      POINTER                                :: left_set, right_set
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CALL dbcsr_multiply_anytype(transa, transb,&
         dbcsr_scalar(alpha), matrix_a, matrix_b, dbcsr_scalar(beta), matrix_c,&
         first_row, last_row, first_column, last_column, first_k, last_k,&
         retain_sparsity, left_set=left_set, right_set=right_set, error=error, flop=flop)
  END SUBROUTINE dbcsr_multiply_z

! *****************************************************************************
!> \brief Scales a DBCSR matrix by alpha
!> \param[inout] matrix_a       DBCSR matrix
!> \param[in] alpha_scalar      (optional) a scalar
!> \param[in] alpha_matrix      (optional) a list of matrices
!> \param[in] side              (optional) ?
!> \param[in] last_column       (optional) only scale up to this full column
!> \param[in] scale type        corresponds to scalar, vector, matrix
! *****************************************************************************
  SUBROUTINE dbcsr_scale_anytype(matrix_a, last_column, error, scale_type,&
         side,&
         alpha_scalar, alpha_matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    CHARACTER, INTENT(IN)                    :: scale_type
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: side
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: alpha_scalar
    TYPE(dbcsr_data_obj), INTENT(IN), &
      OPTIONAL                               :: alpha_matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_anytype', &
      routineP = moduleN//':'//routineN

    INTEGER :: a_blk, a_col, a_col_size, a_nze, a_row, a_row_size, &
      col_offset, error_handler, i, last_col, row_offset
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: col_end, col_start, m_offset
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: do_scale, limit_col, right, tr
    TYPE(dbcsr_data_obj)                     :: buffer, data_any, mat_tmp
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (.NOT. PRESENT (last_column) .OR. scale_type .EQ. 'S', &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Column limiting supported only with scalar scaling",__LINE__,error)
    last_col = 0
    IF (PRESENT (last_column)) THEN
       CALL dbcsr_assert (last_column .GE. 0 .AND. last_column .LE.&
            dbcsr_nfullcols_total (matrix_a), dbcsr_fatal_level,&
            dbcsr_wrong_args_error, routineN,&
            "Specified last column is out of bounds.",__LINE__,error)
       last_col = last_column
       ALLOCATE (col_start(dbcsr_nblkcols_total(matrix_a)))
       ALLOCATE (col_end(dbcsr_nblkcols_total(matrix_a)))
       CALL convert_sizes_to_offsets(array_data (matrix_a%m%col_blk_size), &
            col_start, col_end)
    ENDIF
    !
    ! set vars
    right = .TRUE.
    IF(PRESENT(side)) THEN
       SELECT CASE(side)
       CASE('right');right = .TRUE.
       CASE('left' );right = .FALSE.
       CASE DEFAULT
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "wrong side="//side,__LINE__,error)
       END SELECT
    ENDIF
    IF(scale_type .EQ. 'M') THEN     ! Matrix scaling
       CALL dbcsr_data_init (buffer)
       ! Clean up these maxvalues
       CALL dbcsr_data_new (buffer, dbcsr_get_data_type(matrix_a),&
            data_size=MAXVAL(array_data (matrix_a%m%row_blk_size)) &
            *MAXVAL(array_data (matrix_a%m%col_blk_size)))
       CALL dbcsr_data_clear (buffer) ! Zero the buffer
       !
       ! TODO: Clean up these offsets using the matrix-supplied ones
       row_blk_size => array_data (matrix_a%m%row_blk_size)
       col_blk_size => array_data (matrix_a%m%col_blk_size)
       ALLOCATE(m_offset(dbcsr_nblkrows_total(matrix_a)))
       m_offset(1) = 1
       IF(right) THEN
          DO a_col = 2,dbcsr_nblkcols_total(matrix_a)
             a_col_size = col_blk_size(a_col-1)
             m_offset(a_col) = m_offset(a_col-1) + a_col_size**2
          ENDDO
       ELSE
          DO a_row = 2,dbcsr_nblkrows_total(matrix_a)
             a_row_size = row_blk_size(a_row-1)
             m_offset(a_row) = m_offset(a_row-1) + a_row_size**2
          ENDDO
       ENDIF
       !
       ! This data area holds the current column
       CALL dbcsr_data_init (mat_tmp)
       CALL dbcsr_data_new (mat_tmp, dbcsr_get_data_type (alpha_matrix))
    ENDIF
    !
    ! quick return if possible
    do_scale = .TRUE.
    IF (PRESENT (alpha_scalar)) THEN
       SELECT CASE (alpha_scalar%data_type)
       CASE (dbcsr_type_real_4)
          IF (alpha_scalar%r_sp .EQ. 1.0_real_4) do_scale = .FALSE.
       CASE (dbcsr_type_real_8)
          IF (alpha_scalar%r_dp .EQ. 1.0_real_8) do_scale = .FALSE.
       CASE (dbcsr_type_complex_4)
          IF (alpha_scalar%c_sp .EQ. CMPLX(1.0, 0.0, real_4)) do_scale = .FALSE.
       CASE (dbcsr_type_complex_8)
          IF (alpha_scalar%c_dp .EQ. CMPLX(1.0, 0.0, real_8)) do_scale = .FALSE.
       END SELECT
    END IF
    !
    ! let's go
    IF(do_scale) THEN
       CALL dbcsr_data_init (data_any)
       CALL dbcsr_data_new (data_any, dbcsr_get_data_type (matrix_a))
       CALL dbcsr_iterator_start(iter, matrix_a)
       iterations: DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block(iter, a_row, a_col, data_any, tr,&
               block_number=a_blk,&
               row_size=a_row_size, col_size=a_col_size, &
               row_offset=row_offset, col_offset=col_offset)
          ! ### UB: Find where the transposed blocks come from.
          !CALL cp_assert (.NOT. tr, dbcsr_fatal_level, cp_unimplemented_error_nr,&
          !     routineN, "No support for transposed blocks.")
          a_nze = a_row_size * a_col_size
          IF (a_nze .EQ. 0) CYCLE ! Skip empty blocks
          !
          ! let's scale
          IF(scale_type .EQ. 'S') THEN
             limit_col = last_col .GT. 0
             IF (limit_col) limit_col = col_end(a_col).gt.last_col
             IF (limit_col) THEN
                IF (col_start(a_col) .LE. last_col) THEN
                   i = MIN(last_col, col_end(a_col)) - col_start(a_col) + 1
                   a_nze = a_row_size*i
                ELSE
                   ! Nothing to do, the block is beyond the boundry.
                   a_nze = 0
                ENDIF
             ENDIF
             IF (a_nze .GT. 0) THEN
                CALL dbcsr_blas_scal (a_nze, alpha_scalar, data_any, 1)
             ENDIF
          ELSEIF(scale_type .EQ. 'M') THEN
             IF(right) THEN
                !A = A * alpha
                mat_tmp = pointer_view (mat_tmp, alpha_matrix, m_offset(a_col))
                CALL dbcsr_blas_gemm ('N','N',&
                     a_row_size, a_col_size, a_col_size,&
                     dbcsr_scalar_one (alpha_matrix%d%data_type),&
                     data_any, a_row_size,&
                     mat_tmp, a_col_size,&
                     dbcsr_scalar_zero (alpha_matrix%d%data_type),&
                     buffer, a_row_size)
                !CALL dgemm('N','N',a_row_size,a_col_size,a_col_size,&
                !     1.0_dp,data_any%d%r_dp(1),a_row_size,&
                !     alpha_matrix%d%r_dp(m_offset(a_col)),a_col_size,&
                !     0.0_dp,buffer%d%r_dp(1),a_row_size)
                CALL dbcsr_blas_copy (a_nze, buffer, 1, data_any, 1)
             ELSE
                !A = alpha * A
                mat_tmp = pointer_view (mat_tmp, alpha_matrix, m_offset(a_row))
                CALL dbcsr_blas_gemm ('N','N',&
                     a_row_size, a_col_size, a_row_size,&
                     dbcsr_scalar_one (alpha_matrix%d%data_type),&
                     mat_tmp, a_row_size,&
                     data_any, a_row_size,&
                     dbcsr_scalar_zero (alpha_matrix%d%data_type),&
                     buffer, a_row_size)
                !CALL dgemm('N','N',a_row_size,a_col_size,a_row_size,&
                !        &     1.0_dp,alpha_matrix%d%r_dp(m_offset(a_col)),a_row_size,&
                !        &            data_any%d%r_dp(1),a_row_size,&
                !        &     0.0_dp,buffer%d%r_dp(1),a_row_size)
                CALL dbcsr_blas_copy(a_nze, buffer, 1, data_any, 1)
             ENDIF
          ELSE
             CALL dbcsr_assert (.FALSE., dbcsr_fatal_level,&
                  dbcsr_wrong_args_error, routineN,&
                  "wrong logic",__LINE__,error)
          ENDIF
       ENDDO iterations
       CALL dbcsr_iterator_stop(iter)
       CALL dbcsr_data_clear_pointer (data_any)
       CALL dbcsr_data_release (data_any)
    ENDIF

    IF(scale_type .EQ. 'M') THEN
       DEALLOCATE(m_offset)
       !@@@
       WRITE(*,*)"refs:", buffer%d%refcount, buffer%d%data_type
       CALL dbcsr_data_release (buffer)
       CALL dbcsr_data_clear_pointer (mat_tmp)
       CALL dbcsr_data_release (mat_tmp)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_scale_anytype

  !> \brief Interface for matrix scaling by a scalar
  SUBROUTINE dbcsr_scale_d(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(real_8), INTENT(IN)                 :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_d', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: scale_type
    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
    scale_type = 'S'
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_8) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (alpha_scalar))
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (REAL(alpha_scalar,real_4)))
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid data type",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_scale_d

  SUBROUTINE dbcsr_scale_s(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    REAL(real_4), INTENT(IN)                 :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_s', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: scale_type
    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
    scale_type = 'S'
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_real_4) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (alpha_scalar))
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid data type",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_scale_s
  
  SUBROUTINE dbcsr_scale_z(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    COMPLEX(real_8), INTENT(IN)              :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_z', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: scale_type
    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
    scale_type = 'S'
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_8) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (alpha_scalar))
    ELSEIF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_4) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (CMPLX(alpha_scalar,KIND=real_4)))
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid data type",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_scale_z

  SUBROUTINE dbcsr_scale_c(matrix_a, alpha_scalar, last_column, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    COMPLEX(real_4), INTENT(IN)              :: alpha_scalar
    INTEGER, INTENT(IN), OPTIONAL            :: last_column
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scale_c', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: scale_type
    INTEGER                                  :: error_handler

    CALL dbcsr_error_set(routineN, error_handler, error)
    scale_type = 'S'
    IF(dbcsr_get_data_type(matrix_a) .EQ. dbcsr_type_complex_4) THEN
       CALL dbcsr_scale_anytype(matrix_a, last_column, error,&
            scale_type=scale_type, alpha_scalar=dbcsr_scalar (alpha_scalar))
    ELSE
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid data type",__LINE__,error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_scale_c

! *****************************************************************************
!> \brief Scales a DBCSR matrix by alpha
!> \param[inout] matrix_a       DBCSR matrix
!> \param[in] alpha             the scaling vector
!> \param[in] side              apply the scaling from the side
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

    CHARACTER                                :: data_type
    INTEGER                                  :: a_blk, a_col, a_col_size, &
                                                a_nze, a_row, a_row_size, &
                                                col_offset, error_handler, i, &
                                                row_offset
    LOGICAL                                  :: right, tr
    TYPE(dbcsr_data_obj)                     :: data_any
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

! check that alpha and matrix have the same data type
    CALL dbcsr_assert (dbcsr_get_data_type (matrix_a).EQ.alpha%d%data_type, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "wrong data type matrix_a="//&
         dbcsr_get_data_type (matrix_a)//' alpha='//alpha%d%data_type,__LINE__,error)
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
          DO i = 1,a_col_size
             SELECT CASE (data_type)
             CASE (dbcsr_type_real_4)
                CALL sscal(a_row_size,alpha%d%r_sp(col_offset+i-1),&
                     data_any%d%r_sp((i-1)*a_row_size+1),1)
             CASE (dbcsr_type_real_8)
                CALL dscal(a_row_size,alpha%d%r_dp(col_offset+i-1),&
                     data_any%d%r_dp((i-1)*a_row_size+1),1)
             CASE (dbcsr_type_complex_4)
                CALL cscal(a_row_size,alpha%d%c_sp(col_offset+i-1),&
                     data_any%d%c_sp((i-1)*a_row_size+1),1)
             CASE (dbcsr_type_complex_8)
                CALL zscal(a_row_size,alpha%d%c_dp(col_offset+i-1),&
                     data_any%d%c_dp((i-1)*a_row_size+1),1)
             END SELECT
          ENDDO
       ELSE
          DO i = 1,a_row_size
             SELECT CASE (data_type)
             CASE (dbcsr_type_real_4)
                CALL sscal(a_col_size,alpha%d%r_sp(row_offset+i-1),&
                     data_any%d%r_sp(i),a_col_size)
             CASE (dbcsr_type_real_8)
                CALL dscal(a_col_size,alpha%d%r_dp(row_offset+i-1),&
                     data_any%d%r_dp(i),a_col_size)
             CASE (dbcsr_type_complex_4)
                CALL cscal(a_col_size,alpha%d%c_sp(row_offset+i-1),&
                     data_any%d%c_sp(i),a_col_size)
             CASE (dbcsr_type_complex_8)
                CALL zscal(a_col_size,alpha%d%c_dp(row_offset+i-1),&
                     data_any%d%c_dp(i),a_col_size)
             END SELECT
          ENDDO
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
!> \brief Conjugate a DBCSR matrix
!> \param[inout] matrix       DBCSR matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_conjg(matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_conjg', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER                                  :: blk, col, error_handler, row
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
!> \brief add and scale matrices 
!>    A = alpha*A + beta*B or 
!> \param[in,out] matrix_a   DBCSR matrix
!> \param[in] matrix_b       DBCSR matrix
!> \param[in] alpha_scalar   (optional) 
!> \param[in] beta_scalar    (optional)
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_add_anytype(matrix_a, matrix_b, alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: alpha_scalar, beta_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_anytype', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type_a, data_type_b
    INTEGER                                  :: blk, col, error_handler, row, &
                                                size_a, size_b
    LOGICAL                                  :: do_scale, tr
    TYPE(dbcsr_data_obj)                     :: data_block
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_scalar_type)                  :: my_alpha_scalar, &
                                                my_beta_scalar

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

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

    CALL dbcsr_scale(matrix_a,scale_type='S',&
         alpha_scalar=my_alpha_scalar,error=error)

    ! Pre-size work arrays of matrix_a to avoid continuous reallocation.
    size_a = dbcsr_get_data_size (matrix_a)
    size_b = dbcsr_get_data_size (matrix_b)
    IF (size_b .GT. size_a) THEN
       CALL dbcsr_work_create (matrix_a,&
            nblks_guess=matrix_b%m%nblks-matrix_a%m%nblks,&
            sizedata_guess=size_b-size_a,error=error)
    ENDIF

    IF(.NOT.dbcsr_scalar_are_equal(my_beta_scalar,&
         dbcsr_scalar_zero(data_type_b))) THEN
       CALL dbcsr_data_init (data_block)
       CALL dbcsr_data_new (data_block, data_type_b)
       CALL dbcsr_iterator_start(iter, matrix_b)

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

    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_add_anytype

! *****************************************************************************
!> \brief add and scale matrices
!>    A = alpha*A + beta*B or 
!> \par Internal implementation
!>      This routine avoids using work matrices. Instead, it pre-sizes the
!>      matrix so that the add functions
!> \param[in,out] matrix_a   DBCSR matrix
!> \param[in] matrix_b       DBCSR matrix
!> \param[in] alpha_scalar   (optional) scale matrix_a by this factor
!> \param[in] beta_scalar    (optional) scale matrix_b by this factor
!>                           before adding
!> \param[in,out] error      error
! *****************************************************************************
  SUBROUTINE dbcsr_add_reserved(matrix_a, matrix_b,&
       alpha_scalar, beta_scalar, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    TYPE(dbcsr_scalar_type), INTENT(IN), &
      OPTIONAL                               :: alpha_scalar, beta_scalar
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_reserved', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type_a, data_type_b
    INTEGER                                  :: blk, col, error_handler, &
                                                nblkrows, nblks, row
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: b_row_i
    LOGICAL                                  :: do_scale, tr
    TYPE(dbcsr_data_obj)                     :: data_block
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_scalar_type)                  :: my_alpha_scalar, &
                                                my_beta_scalar

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    ! Checks for validity
    CALL dbcsr_assert (dbcsr_valid_index (matrix_a),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Target matrix A must be valid.", __LINE__, error)
    CALL dbcsr_assert (dbcsr_valid_index (matrix_b),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Source matrix B must be valid.", __LINE__, error)
    ! Shortcuts
    data_type_a = dbcsr_get_data_type(matrix_a)
    data_type_b = dbcsr_get_data_type(matrix_b)
    ! Process arguments
    my_alpha_scalar = dbcsr_scalar_one (data_type_a)
    IF(PRESENT(alpha_scalar)) my_alpha_scalar = alpha_scalar
    my_beta_scalar = dbcsr_scalar_one (data_type_b)
    IF(PRESENT(beta_scalar)) my_beta_scalar = beta_scalar
    ! Scale target matrix
    do_scale = dbcsr_scalar_are_equal (&
         my_beta_scalar, dbcsr_scalar_one (data_type_b))
    CALL dbcsr_scale(matrix_a,scale_type='S',&
         alpha_scalar=my_alpha_scalar,error=error)
    !
    IF(.NOT.dbcsr_scalar_are_equal(my_beta_scalar,&
         dbcsr_scalar_zero(data_type_b))) THEN
       ! Reserve the blocks to be added
       nblks = dbcsr_get_num_blocks (matrix_b)
       nblkrows = dbcsr_nblkrows_total (matrix_b)
       ALLOCATE (b_row_i(nblks))
       CALL dbcsr_expand_row_index (matrix_b%m%row_p, b_row_i, nblkrows, nblks)
       CALL dbcsr_reserve_blocks (matrix_a, b_row_i, matrix_b%m%col_i,&
            error=error)
       DEALLOCATE (b_row_i)
       ! Now add the blocks
       CALL dbcsr_data_init (data_block)
       CALL dbcsr_data_new (data_block, data_type_b)
       CALL dbcsr_iterator_start(iter, matrix_b)
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
       CALL dbcsr_assert (dbcsr_valid_index (matrix_a),&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Block reservations seem to be incompletely done.", __LINE__,&
            error=error)
       CALL dbcsr_iterator_stop(iter)
       CALL dbcsr_data_clear_pointer (data_block)
       CALL dbcsr_data_release (data_block)
       !CALL dbcsr_finalize (matrix_a, error=error)
    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_add_reserved


  !> \brief Interface for dbcsr_add
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
!> \brief Performs a non-scattering reduction on a DBCSR matrix.
!> \par Meaning
!>      The blocks on all processes are incorporated according to the
!>      reduction operation into a block located in the target
!>      process, process row, or process column. E.g., if all
!>      processes have block (1,1) and the target is process 0, then
!>      the sum of all (1,1) blocks is placed onto process 0.
!> \param[in] matrix             matrix to reduce
!> \param[in,out] reduced        reduced version of input matrix
!> \param[in] reduction_ target  target process, process row, or process column
!> \param[in] reduce_rows        reduce by rows
!> \param[in] reduce_columns     reduce by columns
!> \param[in] operation          reduction operation
!> \param[in,out] error          error
!> \par Operation types
!>      Specified as a string; currently supported are <ul>
!>      <li> "+" or "SUM" for summation</li></ul>
! *****************************************************************************
  SUBROUTINE dbcsr_reduce (matrix, reduced,&
       reduction_target, reduce_rows, reduce_columns, &
       operation, error)
    TYPE(dbcsr_obj), INTENT(IN), TARGET      :: matrix
    TYPE(dbcsr_obj), INTENT(INOUT), TARGET   :: reduced
    INTEGER, INTENT(IN)                      :: reduction_target
    LOGICAL, INTENT(IN)                      :: reduce_rows, reduce_columns
    CHARACTER(LEN=*), INTENT(IN)             :: operation
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_reduce', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: i_d = 2, i_i = 1
    LOGICAL, PARAMETER                       :: dbg = .TRUE.
    REAL(kind=dp), PARAMETER                 :: comp_eps = 1.0E-06_dp

    CHARACTER                                :: data_type, red_type
    INTEGER :: dst_p, error_handler, half_spread, most_spread, mp_group, &
      my_norm_node, mynode, nsteps, numnodes, spread, src_p, step, tag
    INTEGER, DIMENSION(2)                    :: get_sizes, my_sizes
    LOGICAL                                  :: i_recv, i_send
    REAL(kind=dp)                            :: cs1, cs2
    TYPE(dbcsr_data_obj)                     :: data_ptr
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(dbcsr_obj), POINTER                 :: mat_comm
    TYPE(dbcsr_obj), TARGET                  :: mat_recv, mat_send

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_valid_index (matrix%m), dbcsr_fatal_level,&
         dbcsr_caller_error, routineN, "Source matrix not initialized.",&
         __LINE__, error=error)
    CALL dbcsr_assert (operation, "EQ", "+", dbcsr_fatal_level,&
         dbcsr_unimplemented_error_nr, routineN,&
         "Reduction operation not yet implemented.", __LINE__, error=error)
    IF (reduce_rows .AND. reduce_columns) THEN
       red_type = dbcsr_repl_full
    ELSEIF (reduce_rows .AND. .NOT. reduce_columns) THEN
       red_type = dbcsr_repl_row
    ELSEIF (reduce_columns .AND. .NOT. reduce_rows) THEN
       red_type = dbcsr_repl_col
    ELSE
       red_type = dbcsr_repl_none
       CALL dbcsr_assert(reduce_rows, "OR", reduce_columns, &
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "At least one of row or column reduction must be specified",&
            __LINE__, error=error)
    ENDIF
    ! Remember stuff for easier handling.
    mp_obj = dbcsr_distribution_mp (dbcsr_distribution (matrix))
    data_type = dbcsr_get_data_type (matrix)
    ! Set communication variables.
    SELECT CASE (red_type)
    CASE (dbcsr_repl_full)
       numnodes = dbcsr_mp_numnodes (mp_obj)
       mp_group = dbcsr_mp_group (mp_obj)
       mynode = dbcsr_mp_mynode (mp_obj)
    CASE (dbcsr_repl_row)
       numnodes = dbcsr_mp_npcols (mp_obj)
       CALL dbcsr_mp_grid_setup (mp_obj, force=.TRUE.)
       mp_group = dbcsr_mp_my_row_group (mp_obj)
       mynode = dbcsr_mp_mynode (mp_obj)
    CASE (dbcsr_repl_col)
       numnodes = dbcsr_mp_nprows (mp_obj)
       CALL dbcsr_mp_grid_setup (mp_obj, force=.TRUE.)
       mp_group = dbcsr_mp_my_col_group (mp_obj)
       mynode = dbcsr_mp_mynode (mp_obj)
    CASE (dbcsr_repl_none)
       numnodes = 1
       mp_group = dbcsr_mp_group (mp_obj)
       mynode = dbcsr_mp_mynode (mp_obj)
    END SELECT
    ! Check that we actually have the row/column communicators when
    ! doing row or column-limited reduction.
    CALL dbcsr_assert (red_type.EQ.dbcsr_repl_row &
         .OR. red_type.EQ.dbcsr_repl_col,&
         "IMP", dbcsr_mp_has_subgroups (mp_obj), dbcsr_fatal_level,&
         dbcsr_unimplemented_error_nr, routineN,&
         "Only full reduction supported when subcommunicators are turned off.",&
         __LINE__, error=error)
    ! Check the reduction_target is sane
    CALL dbcsr_assert (reduction_target, "GE", 0, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Can not have non-negative target.",&
         __LINE__, error=error)
    CALL dbcsr_assert (reduction_target, "LT", numnodes, dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Target must fit process grid.", &
         __LINE__, error=error)
    !
    ! Send/receive buffer setup
    CALL dbcsr_init (mat_send)
    CALL dbcsr_init (mat_recv)
    CALL dbcsr_create (reduced, "Reduced "//matrix%m%name,&
         dbcsr_distribution (matrix), dbcsr_get_matrix_type (matrix),&
         dbcsr_row_block_sizes (matrix), dbcsr_col_block_sizes (matrix),&
         data_type=data_type, error=error)
    CALL dbcsr_finalize (reduced, error=error)
    CALL dbcsr_create (mat_recv, "Recv buffer for "//matrix%m%name,&
         dbcsr_distribution (matrix), dbcsr_get_matrix_type (matrix),&
         dbcsr_row_block_sizes (matrix), dbcsr_col_block_sizes (matrix),&
         data_type=data_type, special=.TRUE., error=error)
    CALL dbcsr_data_init (data_ptr)
    CALL dbcsr_data_new (data_ptr, data_type)
    !
    IF (dbg) THEN
       cs1 = dbcsr_checksum (matrix, error=error)
    ENDIF
    CALL dbcsr_insert_blocks (reduced, matrix, error=error)
    ! Now do the iterative folding in ceil(log_2(numnodes)) steps.
    nsteps = ceil_log2 (numnodes)
    ! Write folding info.
    !write(*,*)routineN//" ceil_log2", numnodes, nsteps,&
    !     32, ceil_log2(32), 7, ceil_log2(7), 2, ceil_log2(2)
    !write(*,*)routineN//" Reduction target:", reduction_target
    most_spread = ISHFT(1, nsteps)
    my_norm_node = MODULO (mynode-reduction_target, numnodes)
    step = 1
    NULLIFY (mat_comm)
    log_steps: DO WHILE (step .LE. nsteps)
       !
       half_spread = ISHFT (1, step-1)
       spread = 2 * half_spread
       ! Write step info
       !WRITE(*,*)routineN//" Step", step, half_spread, spread
       ! Check whether we participate in a send/recv according to the
       ! rules and according to the number of nodes.
       !---
       i_recv = MODULO (my_norm_node, spread) .EQ. 0
       src_p = MODULO (mynode+half_spread, numnodes)
       i_recv = i_recv .AND. (my_norm_node+half_spread) .LT. numnodes
       i_send = MODULO (my_norm_node, spread) .EQ. half_spread
       dst_p = MODULO (mynode-half_spread, numnodes)
       i_send = i_send .AND. MODULO(my_norm_node,spread)-half_spread .GE. 0
       ! Write comm info
       !WRITE(*,*)routineN//" i_recv,send", i_recv, src_p, i_send, dst_p
       CALL dbcsr_assert ("NOT", i_recv .AND. i_send, dbcsr_fatal_level,&
            dbcsr_internal_error, routineN, &
            "Can not send and receive in the same step.", __LINE__, error=error)
       !
       if_recv: IF (i_recv) THEN
          mat_comm => mat_recv
          CALL mp_recv (get_sizes, src_p, step, mp_group)
          CALL dbcsr_data_ensure_size (mat_comm%m%data_area, get_sizes(i_d),&
               nocopy=.TRUE., zero_pad=.FALSE., error=error)
          CALL ensure_array_size (mat_comm%m%index, lb=1, ub=get_sizes(i_i),&
               nocopy=.TRUE., special=mat_comm%m%special_memory,&
               zero_pad=.TRUE., error=error)
          CALL mp_recv (mat_comm%m%index(1:get_sizes(i_i)), &
               src_p, step, mp_group)
          IF (get_sizes(i_d) .GT. 0) THEN
             CALL dbcsr_data_set_pointer (area=data_ptr,&
                  rsize=get_sizes(i_d), csize=1,&
                  pointee=mat_comm%m%data_area)
             CALL dbcsr_recv_any (data_ptr, src_p, tag=step, comm=mp_group,&
                  error=error)
          ENDIF
          CALL dbcsr_repoint_index (mat_comm%m)
          mat_comm%m%valid = .TRUE.
          ! Now add the new data into the receiving matrix.
          CALL dbcsr_add_reserved (reduced, mat_comm, error=error)
       ENDIF if_recv
       if_send: IF (i_send) THEN
          mat_comm => reduced
          my_sizes(i_i) = SIZE (mat_comm%m%index)
          my_sizes(i_d) = dbcsr_get_data_size (mat_comm)
          CALL mp_send (my_sizes, dst_p, step, mp_group)
          CALL mp_send (mat_comm%m%index(1:my_sizes(i_i)), dst_p,&
               step, mp_group)
          IF (my_sizes(i_d) .GT. 0) THEN
             CALL dbcsr_data_set_pointer (area=data_ptr,&
                  rsize=my_sizes(i_d), csize=1,&
                  pointee=mat_comm%m%data_area)
             CALL dbcsr_send_any (data_ptr, dst_p, tag=step, comm=mp_group,&
                  error=error)
          ENDIF
       ENDIF if_send
       step = step + 1
    ENDDO log_steps
    CALL dbcsr_data_clear_pointer (data_ptr)
    CALL dbcsr_data_release (data_ptr)
    IF (dbg .AND. (mynode .EQ. reduction_target)) THEN
       cs2 = dbcsr_checksum (reduced, local=.TRUE., error=error)
       IF (ABS(cs1) .GT. comp_eps) THEN
          CALL dbcsr_assert (ABS((cs2 - cs1) / cs1) .LT. comp_eps,&
               dbcsr_warning_level, dbcsr_internal_error, routineN,&
               "Reduced checksum differs", __LINE__, error=error)
       ENDIF
    ENDIF
    CALL dbcsr_error_stop (error_handler, error)
  END SUBROUTINE dbcsr_reduce


! *****************************************************************************
!> \brief Hadamard product
!>    C = A . B (C needs to be different from A and B)
!> \param[in,out] matrix_c        DBCSR matrix
!> \param[in] matrix_a, matrix_b  DBCSR matrix
! *****************************************************************************
  SUBROUTINE dbcsr_hadamard_product(matrix_a, matrix_b, matrix_c, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_hadamard_product', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER                                  :: blk, col, col_size, &
                                                error_handler, nze, row, &
                                                row_size
    LOGICAL                                  :: found, tr_a, tr_b
    TYPE(dbcsr_data_obj)                     :: a_data, b_data, c_data
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

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
!> \brief Replaces blocks with blocks from another matrix.
!>
!> Every block in matrix B is copied over to matrix A, replacing any possibly
!> pre-existing blocks.
!> \param[in,out] matrix_a   replace blocks in this DBCSR matrix
!> \param[in] matrix_b       replace blocks taken from the DBCSR matrix
! *****************************************************************************
  SUBROUTINE dbcsr_replace_blocks(matrix_a, matrix_b, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_replace_blocks', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: data_c
    COMPLEX(KIND=sp), DIMENSION(:, :), &
      POINTER                                :: data_z
    INTEGER                                  :: blk, col, error_handler, row
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    REAL(KIND=sp), DIMENSION(:, :), POINTER  :: data_r
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (dbcsr_get_data_type(matrix_a).EQ.dbcsr_get_data_type(matrix_b), &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, "data types not consistent",__LINE__,error)
    !
    ! let's go
    CALL dbcsr_assert (matrix_a%m%nblkrows_total.EQ.matrix_b%m%nblkrows_total, &
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN, "matrices not consistent",__LINE__,error)

    CALL dbcsr_iterator_start(iter, matrix_b)

    DO WHILE (dbcsr_iterator_blocks_left(iter))

       SELECT CASE (dbcsr_get_data_type(matrix_a))
       CASE (dbcsr_type_real_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_r, tr, blk)
          CALL dbcsr_put_block(matrix_a, row, col, data_r)
       CASE (dbcsr_type_real_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)
          CALL dbcsr_put_block(matrix_a, row, col, data_d, tr,&
               summation=.FALSE.)
       CASE (dbcsr_type_complex_4)
          CALL dbcsr_iterator_next_block(iter, row, col, data_c, tr, blk)
          CALL dbcsr_put_block(matrix_a, row, col, data_c)
       CASE (dbcsr_type_complex_8)
          CALL dbcsr_iterator_next_block(iter, row, col, data_z, tr, blk)
          CALL dbcsr_put_block(matrix_a, row, col, data_z)
       CASE DEFAULT
          CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
               routineN,"Unkown precision",__LINE__,error)
       END SELECT

    ENDDO

    CALL dbcsr_iterator_stop(iter)
    CALL dbcsr_finalize (matrix_a, error=error)
    !
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_replace_blocks

! *****************************************************************************
!> \brief add a constant to the diagonal of a matrix
!> \param[inout] matrix       DBCSR matrix
!> \param[in]    alpha_scalar scalar
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

  SUBROUTINE dbcsr_init_random(matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_random', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, col_size, error_handler, &
                                                hold, iseed(4), mynode, row, &
                                                row_size, stored_col, &
                                                stored_row
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found, tr
    REAL(real_8), ALLOCATABLE, DIMENSION(:)  :: rnd
    REAL(real_8), DIMENSION(:, :), POINTER   :: buff, data_d

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    row_blk_size => array_data (matrix%m%row_blk_size)
    col_blk_size => array_data (matrix%m%col_blk_size)
    mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(dbcsr_distribution(matrix)))
    CALL dbcsr_work_create(matrix, work_mutable=.TRUE., error=error)

    iseed(1)=4;iseed(2)=3;iseed(3)=2;iseed(4)=1! set the seed for dlarnv
    ALLOCATE(rnd(MAXVAL(row_blk_size)*MAXVAL(col_blk_size)))
    DO row = 1,dbcsr_nblkrows_total(matrix)
    DO col = 1,dbcsr_nblkcols_total(matrix)
       row_size = row_blk_size(row)
       col_size = col_blk_size(col)
       CALL dlarnv(1,iseed,row_size*col_size,rnd(1))
       tr = .FALSE.
       stored_row = row
       stored_col = col
       CALL dbcsr_get_stored_coordinates(matrix, stored_row, stored_col, tr, hold)
       IF(hold.EQ.mynode) THEN
          CALL dbcsr_get_block_p(matrix, stored_row, stored_col, data_d, tr, found)
          IF(found) THEN
             CALL dcopy(row_size*col_size,rnd,1,data_d,1)
          ELSE
             ALLOCATE(buff(row_size,col_size))
             CALL dcopy(row_size*col_size,rnd,1,buff,1)
             CALL dbcsr_put_block (matrix, stored_row, stored_col, buff)
             DEALLOCATE(buff)
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
    CALL dbcsr_create(diag, 'diag of '//matrix%m%name, &
         matrix%m%dist, dbcsr_get_matrix_type (matrix),&
         matrix%m%row_blk_size,&
         matrix%m%col_blk_size, matrix%m%nblks, matrix%m%nze,&
         matrix%m%data_type, matrix%m%special_memory, error=error)

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
! *****************************************************************************
  SUBROUTINE dbcsr_get_diag_anytype(matrix, diag, error)

    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: diag
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_diag_anytype', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER                                  :: blk, col, error_handler, row, &
                                                row_offset, row_size
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_a, diag_a
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (dbcsr_nfullrows_total(matrix), "LE", dbcsr_get_data_size(diag), &
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
! *****************************************************************************
  SUBROUTINE dbcsr_set_diag_anytype(matrix, diag, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_data_obj), INTENT(IN)         :: diag
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_diag_anytype', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: dt
    INTEGER                                  :: blk, col, error_handler, row, &
                                                row_offset, row_size
    LOGICAL                                  :: tr
    TYPE(dbcsr_data_obj)                     :: data_a, diag_a
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_assert (dbcsr_nfullrows_total(matrix), "LE", dbcsr_get_data_size(diag), &
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
!> \brief copy a matrix 
!> \param[inout] matrix_b       target DBCSR matrix
!> \param[in]    matrix_a       source DBCSR matrix
!> \param[in]    name           (optional) name of the new matrix
!> \param[in,out] error         cp2k error
!> \param[in]    keep_sparsity  (optional) keep the target matrix sparsity;
!>                              default is False.
!> \param[in]    shallow_data   (optional) shallow data copy
!> \param[in]    keep_imaginary  (optional) when copy from complex to real,&
!>                               the default is to keep only the real part; if
!>                               this flag is set, the imaginary part is used
! *****************************************************************************
  SUBROUTINE dbcsr_copy(matrix_b, matrix_a, name, keep_sparsity,&
       shallow_data, keep_imaginary, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: name
    LOGICAL, INTENT(IN), OPTIONAL            :: keep_sparsity, shallow_data, &
                                                keep_imaginary
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: new_type, repl_type
    INTEGER                                  :: error_handler
    LOGICAL                                  :: keep_sparse, shallow

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    shallow = .FALSE. ; IF (PRESENT (shallow_data)) shallow = shallow_data
    keep_sparse = .FALSE.
    IF (PRESENT (keep_sparsity)) keep_sparse = keep_sparsity

    !CALL dbcsr_assert (dbcsr_is_initialized (matrix_b),&
    !     dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
    !     "Target matrix must be initialized", error=error)
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
       IF (dbcsr_is_initialized (matrix_b)) THEN
          new_type = dbcsr_get_data_type (matrix_b)
          repl_type = dbcsr_get_replication_type(matrix_b)
       ELSE
          new_type = dbcsr_get_data_type (matrix_a)
          repl_type = dbcsr_get_replication_type(matrix_a)
       ENDIF
       IF(PRESENT(name)) THEN
          CALL dbcsr_create(matrix_b, name, &
               dbcsr_distribution (matrix_a), dbcsr_get_matrix_type (matrix_a),&
               matrix_a%m%row_blk_size, matrix_a%m%col_blk_size,&
               nblks=matrix_a%m%nblks, data_type=new_type,&
               special=matrix_a%m%special_memory, &
               replication_type=repl_type, error=error)
       ELSE
          CALL dbcsr_create(matrix_b,'copy of '//dbcsr_name(matrix_a),&
               dbcsr_distribution (matrix_a), dbcsr_get_matrix_type (matrix_a),&
               matrix_a%m%row_blk_size, matrix_a%m%col_blk_size,&
               nblks=matrix_a%m%nblks, data_type=new_type,&
               special=matrix_a%m%special_memory,&
               replication_type=repl_type, error=error)
       ENDIF
       CALL ensure_array_size(matrix_b%m%index, ub=SIZE(matrix_a%m%index), error=error)
       !
       ! copy index and data
       matrix_b%m%index(1:SIZE(matrix_a%m%index)) = matrix_a%m%index(:)
       IF (.NOT. shallow) THEN
          CALL dbcsr_assert (matrix_a%m%nze .LE. dbcsr_get_data_size(matrix_a),&
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
          CALL dbcsr_data_release (matrix_b%m%data_area)
          matrix_b%m%data_area = matrix_a%m%data_area
          CALL dbcsr_data_hold (matrix_b%m%data_area)
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

    CHARACTER                                :: data_type
    INTEGER                                  :: col_size, cs, dst_col, &
                                                dst_row, error_handler, rel, &
                                                row_size, rs, src_col, src_row
    LOGICAL                                  :: dst_tr, src_tr
    TYPE(dbcsr_data_obj)                     :: dst_data, src_data
    TYPE(dbcsr_iterator)                     :: dst_iter, src_iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_get_data_type(matrix_b)&
         .EQ. dbcsr_get_data_type(matrix_a), dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Matrices have different data types.",__LINE__,error)
    data_type = dbcsr_get_data_type (matrix_b)

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
               src_data, src_tr)
          rel = pos_relation (dst_row, dst_col, src_row, src_col)
       ENDDO
       SELECT CASE (rel)
       CASE (-1, 1)
           ! Target lags source or ran out of source
          CALL dbcsr_data_clear (dst_data)
          !dst_r_dp(:) = 0.0_dp
       CASE (0)
          ! Copy the data
          CALL dbcsr_assert (dbcsr_get_data_size (src_data)&
               .EQ. dbcsr_get_data_size (dst_data), dbcsr_fatal_level,&
               dbcsr_internal_error, routineN, "Block sizes not equal!",__LINE__,error)
          IF (src_tr .EQV. dst_tr) THEN
             CALL dbcsr_data_copyall (dst_data, src_data, error=error)
             !dst_r_dp(:) = src_r_dp(:)
          ELSE
             rs = row_size
             cs = col_size
             IF (src_tr) CALL swap (rs, cs)
             CALL dbcsr_block_transpose (dst_data, src_data, rs, cs, &
                  error=error)
             !CALL dbcsr_block_transpose (dst_r_dp, src_r_dp, rs, cs)
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
!> \brief Copy only a subset of matrix columns (element-wise, not block-size)
!> \par Column block compatibility
!>      The process column distributions of the copied columns must be
!>      equal. Therefore, the source and target offsets must also be
!>      the same.
!> \param[in,out] matrix_b      target DBCSR matrix
!> \param[in]    matrix_a       source DBCSR matrix
!> \param[in]    ncol           number of elemental column to copy
!> \param[in]    source_start   starting elemental column in source
!> \param[in]    target_start   starting elemental column in target
!> \param[in,out] error         cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_copy_columns(matrix_b, matrix_a,&
       ncol, source_start, target_start, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    INTEGER, INTENT(IN)                      :: ncol, source_start, &
                                                target_start
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_columns', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: data_type
    INTEGER :: col_size, dst_col, dst_row, error_handler, first_col, &
      last_col, num_col, rel, row_size, src_col, src_row
    INTEGER, DIMENSION(:), POINTER           :: dst_col_offsets, &
                                                dst_col_sizes, &
                                                src_col_offsets, src_col_sizes
    LOGICAL                                  :: dst_tr, src_tr
    TYPE(dbcsr_data_obj)                     :: dst_data, src_data, tmp_buffer
    TYPE(dbcsr_iterator)                     :: dst_iter, src_iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_assert (dbcsr_get_data_type(matrix_b)&
         .EQ. dbcsr_get_data_type(matrix_a), dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Matrices have different data types.",__LINE__,error)
    data_type = dbcsr_get_data_type (matrix_b)

    CALL dbcsr_assert (source_start .EQ. target_start, dbcsr_fatal_level,&
         dbcsr_unimplemented_error_nr, routineN, "Column shifting not supported",__LINE__,error)

    src_col_offsets => array_data (matrix_a%m%col_blk_offset)
    dst_col_offsets => array_data (matrix_b%m%col_blk_offset)
    src_col_sizes => array_data (matrix_a%m%col_blk_size)
    dst_col_sizes => array_data (matrix_b%m%col_blk_size)

    CALL dbcsr_data_init (src_data)
    CALL dbcsr_data_init (dst_data)
    CALL dbcsr_data_init (tmp_buffer)
    CALL dbcsr_data_new (src_data, data_type)
    CALL dbcsr_data_new (dst_data, data_type)
    CALL dbcsr_data_new (tmp_buffer, data_type,&
         matrix_a%m%max_rbs*matrix_a%m%max_cbs)
    CALL dbcsr_iterator_start (src_iter, matrix_a)
    CALL dbcsr_iterator_start (dst_iter, matrix_b)
    ! Iterate through the blocks of the source and destination
    ! matrix. There are three possibilites: 1 copy data for blocks
    ! present in both that are fully within the specified range; 2
    ! copy partial data for blacks partially within the specified
    ! range; 3 add (partial or full) data for blocks present in source
    ! but not in the target
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
               src_data, src_tr)
          rel = pos_relation (dst_row, dst_col, src_row, src_col)
       ENDDO
       SELECT CASE (rel)
       CASE (0) ! Equal block coordinates; block exists in both src and dst
          ! Check if the full block is being copied
          CALL dbcsr_assert (src_col_sizes(src_col) .EQ. dst_col_sizes(dst_col),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Unequal column sizes",__LINE__,error)
          !CALL dbcsr_unimplemented_error (routineN, "Partial block copy not tested",&
          !     error=error, error_level=dbcsr_warning_level)
          CALL dbcsr_assert (src_col_offsets(src_col).eq.dst_col_offsets(dst_col),&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Block offsets must be equal",__LINE__,error)
          first_col = 1 + MAX (src_col_offsets(src_col), source_start)&
               - src_col_offsets(src_col)
          last_col = MIN (src_col_offsets(src_col)+ncol-1, source_start+ncol-1) - &
               src_col_offsets(src_col+1)-1 + 1
          num_col = last_col - first_col + 1
          IF (num_col .GT. 0) THEN
             CALL dbcsr_block_partial_copy(&
                  dst_data, row_size, col_size, dst_tr,&
                  src_data, row_size, col_size, src_tr,&
                  1, first_col, 1, first_col,& ! offsets
                  row_size, num_col) ! sizes
          ENDIF
       CASE (-1) ! Block exists in dst but must be added to dst
          IF (source_start .GE. src_col_offsets(src_col) .AND. &
              src_col_offsets(src_col)+ncol-1 .LE. source_start+ncol-1) THEN
             CALL dbcsr_put_block (matrix_b, src_row, src_col, src_data,&
                  transposed=src_tr)
          ELSE
             CALL dbcsr_assert (.FALSE., dbcsr_fatal_level,  dbcsr_internal_error, routineN,&
                  "Partial block addition not yet tested",__LINE__,error)
             CALL dbcsr_data_clear (tmp_buffer)
             CALL dbcsr_block_partial_copy (&
                  tmp_buffer, row_size, col_size, src_tr,&
                  src_data, row_size, col_size, src_tr,&
                  1, first_col, 1, first_col,&
                  row_size, num_col)
             CALL dbcsr_put_block (matrix_b, src_row, src_col, tmp_buffer,&
                  src_tr)
          ENDIF
       CASE (1)
          ! Ran out of source. Do nothing.
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
  END SUBROUTINE dbcsr_copy_columns


! *****************************************************************************
!> \brief Determines the relation between two matrix positions.
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

    CHARACTER                                :: new_type, repl_type
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
    IF (dbcsr_is_initialized (matrix_b)) THEN
       new_type = dbcsr_get_data_type (matrix_b)
       repl_type = dbcsr_get_replication_type(matrix_b)
    ELSE
       new_type = dbcsr_get_data_type (matrix_a)
       repl_type = dbcsr_get_replication_type(matrix_a)
    ENDIF
    ! Setup target matrix
    IF(PRESENT(name)) THEN
       CALL dbcsr_create(matrix_b, name, &
            dbcsr_distribution (matrix_a), dbcsr_get_matrix_type (matrix_a),&
            matrix_a%m%row_blk_size, matrix_a%m%col_blk_size,&
            nblks=matrix_a%m%nblks, data_type=new_type,&
            special=matrix_a%m%special_memory, &
            replication_type=repl_type, error=error)
    ELSE
       CALL dbcsr_create(matrix_b,'submatrix copy of '//TRIM(dbcsr_name(matrix_a)),&
            dbcsr_distribution (matrix_a), dbcsr_get_matrix_type (matrix_a),&
            matrix_a%m%row_blk_size, matrix_a%m%col_blk_size,&
            nblks=matrix_a%m%nblks, data_type=new_type,&
            special=matrix_a%m%special_memory,&
            replication_type=repl_type, error=error)
    ENDIF
    CALL dbcsr_finalize(matrix_b, error=error)
    IF (.NOT. shallow) THEN
       ! Non-shallow copy uses the standard iterator on the source and
       ! block put on the target.
       CALL dbcsr_data_init (data_block)
       CALL dbcsr_data_new (data_block, dbcsr_get_data_type (matrix_a))
       CALL dbcsr_iterator_start (iter, matrix_a)
       DO WHILE (dbcsr_iterator_blocks_left (iter))
          CALL dbcsr_iterator_next_block (iter, row, col, data_block, tr)
          ! Only keep the block if they are within the specified bounds.
          IF (PRESENT (block_row_bounds)) THEN
             IF (row .LT. block_row_bounds(1)) CYCLE
             IF (row .GT. block_row_bounds(2)) CYCLE
          ENDIF
          IF (PRESENT (block_column_bounds)) THEN
             IF (col .LT. block_column_bounds(1)) CYCLE
             IF (col .GT. block_row_bounds(2)) CYCLE
          ENDIF
          CALL dbcsr_put_block (matrix_b, row, col, data_block, transposed=tr)
       END DO
       CALL dbcsr_iterator_stop (iter)
       CALL dbcsr_data_clear_pointer (data_block)
       CALL dbcsr_data_release (data_block)
       CALL dbcsr_finalize (matrix_b, error=error)
    ELSE
       ! For the shallow copy the source matrix data is referenced.
       CALL dbcsr_data_release (matrix_b%m%data_area)
       matrix_b%m%data_area = matrix_a%m%data_area
       CALL dbcsr_data_hold (matrix_b%m%data_area)
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
             IF (col .GT. block_row_bounds(2)) CYCLE
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
!> \brief block triu of a dbcsr matrix
!> \param[in] matrix_b  triu of the matrix
!> \param[in] matrix_a  the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_btriu(matrix_b, matrix_a, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_b
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_btriu', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_nze, col, col_size, &
                                                error_handler, row, row_size
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_copy(matrix_b, matrix_a, name="triu of "//matrix_a%m%name, error=error)

    CALL dbcsr_iterator_start(iter, matrix_b)

    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, blk=blk,&
            row_size=row_size, col_size=col_size)
       blk_nze = row_size * col_size
       IF(row.GT.col) CALL dbcsr_remove_block(matrix_b, row, col, blk_nze, blk)
    ENDDO

    CALL dbcsr_iterator_stop(iter)

    CALL dbcsr_finalize (matrix_b, error=error)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_btriu

! *****************************************************************************
!> \brief  triu of a dbcsr matrix
!> \param[inout] matrix_a  the matrix
!> 
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
!> \brief  symmetrize the block diagonal of a dbcsr matrix
!> \param[inout] matrix_a  the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_symmetrize_block_diag(matrix_a, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_symmetrize_block_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_size, &
                                                error_handler, i, j, row, &
                                                row_size
    LOGICAL                                  :: tr
    REAL(dp)                                 :: dum
    REAL(dp), DIMENSION(:, :), POINTER       :: DATA
    TYPE(dbcsr_iterator)                     :: iter

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    CALL dbcsr_iterator_start(iter, matrix_a)

    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, DATA, tr, &
            block_number=blk, row_size=row_size, col_size=col_size)
       IF(col.NE.row) CYCLE
       DO j = 1,col_size
       DO i = j+1,row_size
          dum = (DATA(i,j) + DATA(j,i)) / 2.0_dp
          DATA(i,j) = dum
          DATA(j,i) = dum
       ENDDO
       ENDDO
    ENDDO

    CALL dbcsr_iterator_stop(iter)

    CALL dbcsr_finalize(matrix_a, error=error)
    CALL dbcsr_error_stop(error_handler, error)

  END SUBROUTINE dbcsr_symmetrize_block_diag

! *****************************************************************************
!> \brief  tril of a dbcsr matrix
!> \param[inout] matrix_a  the matrix
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_tril(matrix_a, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_tril', &
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
       CALL dbcsr_iterator_next_block(iter, row, col, block=DATA, transposed=tr, &
            block_number=blk, row_size=row_size, col_size=col_size)
       blk_nze = row_size * col_size
       IF(row.GT.col) CALL dbcsr_remove_block(matrix_a, row, col, blk_nze, blk)
       IF(row.EQ.col) THEN
          DO j = 1,col_size
          DO i = 1,j-1
             DATA(i,j) = 0.0_dp
          ENDDO
          ENDDO
       ENDIF
    ENDDO

    CALL dbcsr_iterator_stop(iter)

    CALL dbcsr_finalize (matrix_a, error=error)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_tril

! *****************************************************************************
!> \brief filter a dbcsr matrix
!> \param[inout] matrix_a  the matrix
!> \param[in] eps          the threshold
!> \param[in] method       how the matrix is filtered 
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_filter_anytype(matrix, eps, method, use_absolute, filter_diag, error)

    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(dbcsr_scalar_type), INTENT(IN)      :: eps
    INTEGER, INTENT(IN), OPTIONAL            :: method
    LOGICAL, INTENT(in), OPTIONAL            :: use_absolute, filter_diag
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_filter_anytype', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_nze, col, col_size, &
                                                error_handler, my_method, &
                                                row, row_size
    LOGICAL                                  :: gt0, my_filter_diag, tr
    REAL(KIND=real_4)                        :: nrm_s, sdot
    REAL(KIND=real_4), DIMENSION(:), POINTER :: data_s
    REAL(KIND=real_8)                        :: ddot, my_absolute, nrm_d
    REAL(KIND=real_8), DIMENSION(:), POINTER :: data_d
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
             !CASE (dbcsr_type_complex_4)
             !CASE (dbcsr_type_complex_8)
          CASE DEFAULT
             CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_unimplemented_error_nr, &
                  routineN,"Only real filtering",__LINE__,error)
          END SELECT
       ENDDO
       CALL dbcsr_iterator_stop(iter)
       CALL dbcsr_finalize (matrix, error=error)

    ENDIF
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_filter_anytype

! *****************************************************************************
!> \brief compute a norm of a dbcsr matrix
!> \param[in] matrix  the matrix
!> \param[out] norm
!> \param[in] which_norm
!> 
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
          CALL dbcsr_assert (dbcsr_get_data_type(norm_vector), "EQ",&
               dbcsr_get_data_type(matrix), dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN,&
               "Mismatched vector/matrix data types", __LINE__, error=error)
          CALL dbcsr_assert (dbcsr_get_data_size(norm_vector), "GE",&
               dbcsr_nfullrows_total(matrix), dbcsr_fatal_level,&
               dbcsr_wrong_args_error, routineN, "Passed vector too small",&
               __LINE__, error=error)
          CALL dbcsr_assert (dbcsr_get_data_size(norm_vector), "GE",&
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
!> \param[out] norm
!> 
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
!> \param[out] norm
!> 
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

    !dmp_max, this fixes a bug in g95
    ! -> ambigous interface for mp_max
    !    in module dbcsr_message_passing and wrongly
    !    exported from the module timings
    CALL dmp_max(norm, dbcsr_mp_group(dbcsr_distribution_mp(matrix%m%dist)))

  END FUNCTION dbcsr_maxabs_r8

! *****************************************************************************
!> \brief compute a norm of a dbcsr matrix
!> \param[in] matrix  the matrix
!> \param[out] norm
!> 
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
          norm = norm + fac * SUM(CONJG(data_z)*data_z)
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
!> \brief Sums blocks in a replicated dbcsr matrix.
!> \param[in,out] matrix      dbcsr matrix to operate on
! *****************************************************************************
  SUBROUTINE dbcsr_sum_replicated (matrix, error)
    TYPE(dbcsr_obj), INTENT(inout)           :: matrix
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sum_replicated', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: comm, error_handler, &
                                                max_data_size, &
                                                max_index_size, mynode, &
                                                numnodes, p
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: send_index
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: all_sizes
    INTEGER, DIMENSION(2)                    :: my_sizes
    INTEGER, DIMENSION(:), POINTER           :: recv_index
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: send_dp
    REAL(kind=dp), DIMENSION(:), POINTER     :: recv_dp
    TYPE(dbcsr_mp_obj)                       :: mp
    TYPE(dbcsr_obj)                          :: recv_matrix

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handler, error)
    mp = dbcsr_distribution_mp (dbcsr_distribution (matrix))
    comm = dbcsr_mp_group (mp)
    numnodes = dbcsr_mp_numnodes (mp)
    mynode = dbcsr_mp_mynode (mp)
    ! Count data sizes
    my_sizes(1) = SIZE (matrix%m%index)
    recv_dp => dbcsr_get_data_p (matrix%m%data_area, 0.0_dp)
    my_sizes(2) = SIZE (recv_dp)
    ! Allocate receive buffers
    ALLOCATE (all_sizes (2, 0:numnodes-1))
    CALL mp_allgather (my_sizes, all_sizes, comm)
    max_index_size = MAXVAL (all_sizes(1,:))
    max_data_size = MAXVAL (all_sizes(2,:))
    CALL dbcsr_init (recv_matrix)
    CALL dbcsr_create(recv_matrix, "Temporary receiving buffer",&
         dbcsr_distribution(matrix), dbcsr_get_matrix_type(matrix),&
         dbcsr_row_block_sizes(matrix), dbcsr_col_block_sizes(matrix),&
         0, 0, dbcsr_get_data_type(matrix), error=error)
    CALL ensure_array_size (recv_matrix%m%index, ub=max_index_size, error=error)
    CALL dbcsr_data_ensure_size (recv_matrix%m%data_area, max_data_size, error=error)
    ! Setup data to be sent.
    ALLOCATE (send_dp(my_sizes(2)))
    send_dp(:) = recv_dp(:)
    ALLOCATE (send_index(my_sizes(1)))
    send_index(:) = matrix%m%index(:)
    ! Broadcasts a matrix and every processor sums it. Inefficient and
    ! slow, but it gets the job done with little memory use.
    DO p = 0, numnodes-1
       IF (p .EQ. mynode) THEN
          CALL mp_bcast (send_index, p, comm)
          CALL mp_bcast (send_dp, p, comm)
       ELSE
          recv_index => recv_matrix%m%index(1:all_sizes(1,p))
          recv_dp => pointer_view ( dbcsr_get_data_p (&
               recv_matrix%m%data_area, 0.0_dp), 1, all_sizes(2,p))
          CALL mp_bcast (recv_index, p, comm)
          CALL mp_bcast (recv_dp, p, comm)
          CALL dbcsr_repoint_index(recv_matrix%m)
          CALL dbcsr_add (matrix, recv_matrix, 1.0_dp, 1.0_dp, error=error)
        ENDIF
    ENDDO
    CALL dbcsr_destroy (recv_matrix, error=error,force=.TRUE.)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_sum_replicated


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
!> \brief compute the extremal eigenvalues of a symmetric real matrix 
!>        with a (simple) Lanczos approach
!> \param[in] matrix
!> \param[in] max_iter            maximum iteration
!> \param[in] eps                 convergence parameter
!> \param[out] min_eig, max_eig   approximation to the extremal eigenvalues
!> \param[out] approx_norm_2      approximation of the 2 norm (1 < max_iter < 10)
!> \param[out] converged          true if the iteration converged
!> \param[inout] error
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_lanczos_extremal_eig(matrix, max_iter, eps, min_eig, max_eig, &
       approx_norm_2, converged, error)

    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: max_iter
    REAL(dp), INTENT(in), OPTIONAL           :: eps
    REAL(dp), INTENT(out), OPTIONAL          :: min_eig, max_eig, &
                                                approx_norm_2
    LOGICAL, INTENT(OUT), OPTIONAL           :: converged
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_lanczos_extremal_eig', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, i, info, &
                                                iwork, lwork, n
    LOGICAL                                  :: my_converged
    REAL(dp)                                 :: alpha, beta, max_eig_old, &
                                                min_eig_old, my_eps, &
                                                my_max_eig, my_min_eig, &
                                                nrm_f, nrm_v
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: work
    REAL(dp), DIMENSION(max_iter)            :: evals
    REAL(dp), DIMENSION(max_iter, max_iter)  :: T, T0
    TYPE(array_i1d_obj)                      :: col_blk_size, col_dist
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(dbcsr_obj)                          :: f, v, v0

    CALL dbcsr_error_set(routineN, error_handler, error)

    CALL dbcsr_get_info(matrix,nfullrows_total=n)

    CALL create_bl_distribution (col_dist, col_blk_size, 1, &
         dbcsr_mp_npcols(dbcsr_distribution_mp(dbcsr_distribution(matrix))))
    CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp (dbcsr_distribution(matrix)),&
         dbcsr_distribution_row_dist(dbcsr_distribution(matrix)), col_dist)

    CALL dbcsr_init(v)
    CALL dbcsr_create(v, 'v', dist, dbcsr_type_normal, matrix%m%row_blk_size,&
         col_blk_size, 0, 0, error=error)
    CALL dbcsr_finalize(v, error=error)

    CALL dbcsr_init(v0)
    CALL dbcsr_create(v0, 'v0', dist, dbcsr_type_normal, matrix%m%row_blk_size,&
         col_blk_size, 0, 0, error=error)
    CALL dbcsr_finalize(v0, error=error)

    CALL dbcsr_init(f)
    CALL dbcsr_create(f, 'f', dist, dbcsr_type_normal, matrix%m%row_blk_size,&
         col_blk_size, 0, 0, error=error)
    CALL dbcsr_finalize(f, error=error)

    CALL dbcsr_distribution_release (dist)
    CALL array_release (col_blk_size)
    CALL array_release (col_dist)

    lwork = 1+2*max_iter+100
    ALLOCATE(work(lwork))

    my_eps = 1.0e-1_dp
    IF(PRESENT(eps)) my_eps = eps

    min_eig_old = 0.0_dp
    max_eig_old = 0.0_dp
    T(:,:) = 0.0_dp
    ! v = rand(n,1)
    CALL dbcsr_init_random(v,error=error)
    ! v = v / norm(v)
    nrm_v = dbcsr_frobenius_norm(v)
    CALL dbcsr_scale(v,1.0_dp/nrm_v,error=error)
    ! f = A * v
    CALL dbcsr_multiply('N','N',1.0_dp,matrix,v,0.0_dp,f,error=error)
    ! alpha = f' * v
    CALL dbcsr_trace(f,v,alpha,error=error)
    ! f = f - alpha * v
    CALL dbcsr_add(f,v,1.0_dp,-alpha,error=error)
    T(1,1) = alpha
    my_min_eig = alpha; my_max_eig = alpha
    my_converged = .FALSE.
    DO i = 2,max_iter
       ! beta = norm(f)
       beta = dbcsr_frobenius_norm(f)
       ! v0 = v
       CALL dbcsr_copy(v0,v,error=error)
       ! v = f / beta
       CALL dbcsr_add(v,f,0.0_dp,1.0_dp/beta,error=error)
       ! f = A * v
       CALL dbcsr_multiply('N','N',1.0_dp,matrix,v,0.0_dp,f,error=error)
       ! f = f - beta * v0
       CALL dbcsr_add(f,v0,1.0_dp,-beta,error=error)
       ! alpha = f' * v
       CALL dbcsr_trace(f,v,alpha,error=error)
       ! f = f - alpha * v
       CALL dbcsr_add(f,v,1.0_dp,-alpha,error=error)
       T(i  ,i-1) = beta
       T(i-1,i  ) = beta
       T(i  ,i  ) = alpha
       !
       max_eig_old = my_max_eig; min_eig_old = my_min_eig
       T0(:,:) = T(:,:)
       CALL DSYEVD('N','U',i,T0(1,1),max_iter,evals(1),work(1),lwork,iwork,1,info)
       CALL dbcsr_assert(info.EQ.0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "DSYEVD", __LINE__, error)
       my_max_eig = MAXVAL(evals(1:i)); my_min_eig = MINVAL(evals(1:i))
       !write(*,*) routinen//' i',i,'max_eig',my_max_eig,' min_eig',my_min_eig
       IF(ABS(my_max_eig-max_eig_old).LT.my_eps.AND.ABS(my_min_eig-min_eig_old).LT.my_eps) THEN
          my_converged = .TRUE.
          EXIT
       ENDIF
    ENDDO

    IF(PRESENT(approx_norm_2)) THEN
       ! norm(f,2)
       nrm_f = dbcsr_frobenius_norm(f)
       ! norm(T,2)
       T0(:,:) = T(:,:)
       CALL DSYEVD('N','U',max_iter,T0(1,1),max_iter,evals(1),work(1),lwork,iwork,1,info)
       CALL dbcsr_assert(info.EQ.0, dbcsr_fatal_level, dbcsr_internal_error, &
            routineN, "DSYEVD", __LINE__, error)
       ! norm(T,2)+norm(f,2)
       approx_norm_2 = MAXVAL(ABS(evals)) + nrm_f
    ENDIF

    IF(PRESENT(min_eig)) min_eig = my_min_eig
    IF(PRESENT(max_eig)) max_eig = my_max_eig
    IF(PRESENT(converged)) converged = my_converged

    DEALLOCATE(work)
    CALL dbcsr_release(v)
    CALL dbcsr_release(v0)
    CALL dbcsr_release(f)

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_lanczos_extremal_eig



#include "dbcsr_operations_d_.F"
#include "dbcsr_operations_z_.F"
#include "dbcsr_operations_s_.F"
#include "dbcsr_operations_c_.F"

END MODULE dbcsr_operations

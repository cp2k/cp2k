!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR operations in CP2K
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.8
!>
!> <b>Modification history:</b>
!> - Created 2009-05-12
! *****************************************************************************
MODULE cp_dbcsr_operations
  USE array_types,                     ONLY: array_data,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release
  USE cp_blacs_calls,                  ONLY: cp_blacs_pnum
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type,&
                                             cp_blacs_env_write,&
                                             get_blacs_info
  USE cp_fm_basic_linalg,              ONLY: cp_fm_frobenius_norm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_sizes, dbcsr_distribution, &
       dbcsr_distribution_col_dist, dbcsr_distribution_init, &
       dbcsr_distribution_local_cols, dbcsr_distribution_local_rows, &
       dbcsr_distribution_mp, dbcsr_distribution_ncols, &
       dbcsr_distribution_new, dbcsr_distribution_nlocal_cols, &
       dbcsr_distribution_nlocal_rows, dbcsr_distribution_nrows, &
       dbcsr_distribution_release, dbcsr_distribution_row_dist, &
       dbcsr_get_data, dbcsr_get_data_p, dbcsr_get_data_type, dbcsr_get_info, &
       dbcsr_get_matrix_type, dbcsr_init, dbcsr_iterator_blocks_left, &
       dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, &
       dbcsr_mp_group, dbcsr_mp_mynode, dbcsr_mp_mypcol, dbcsr_mp_myprow, &
       dbcsr_mp_npcols, dbcsr_mp_nprows, dbcsr_mp_numnodes, dbcsr_mp_pgrid, &
       dbcsr_mutable_new, dbcsr_nblkcols_total, dbcsr_nblkrows_total, &
       dbcsr_nfullcols_local, dbcsr_nfullcols_total, dbcsr_nfullrows_local, &
       dbcsr_nfullrows_total, dbcsr_release, dbcsr_row_block_sizes, &
       dbcsr_uses_special_memory, dbcsr_valid_index
  USE dbcsr_operations,                ONLY: &
       block_set_d, dbcsr_add, dbcsr_complete_redistribute, dbcsr_copy, &
       dbcsr_distribute, dbcsr_filter, dbcsr_get_block_p, &
       dbcsr_get_tmp_block2d, dbcsr_mult_e_ec_m, dbcsr_norm, dbcsr_put_block, &
       dbcsr_remove_block, dbcsr_replace_blocks, dbcsr_replicate_all, &
       dbcsr_reserve_block2d, dbcsr_scale, dbcsr_set, dbcsr_sum_replicated, &
       dbcsr_trace
  USE dbcsr_types,                     ONLY: &
       dbcsr_array_type, dbcsr_block_node_type, dbcsr_distribution_obj, &
       dbcsr_imagedistribution_type, dbcsr_iterator, dbcsr_mp_obj, dbcsr_obj, &
       dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_complex_4, &
       dbcsr_type_complex_8, dbcsr_type_normal, dbcsr_type_real_4, &
       dbcsr_type_real_8, dbcsr_type_real_default, dbcsr_type_symmetric
  USE dbcsr_util,                      ONLY: &
       addto_array, convert_offsets_to_sizes, convert_sizes_to_offsets, &
       dbcsr_checksum, dbcsr_create, dbcsr_create_col_image_dist, &
       dbcsr_create_dist_elemental, dbcsr_create_dist_fullest, &
       dbcsr_create_dist_r_rot, dbcsr_create_dist_r_unrot, &
       dbcsr_destroy_array, dbcsr_destroy_image_dist, dbcsr_dist2d_to_dist, &
       dbcsr_finalize, dbcsr_find_column, dbcsr_get_block_index, &
       dbcsr_make_col_images, dbcsr_new_transposed, dbcsr_set_debug, &
       dbcsr_verify_matrix, dbcsr_work_create, ensure_array_size, &
       pointer_set, pointer_view, swap
  USE distribution_2d_types,           ONLY: distribution_2d_get,&
                                             distribution_2d_release,&
                                             distribution_2d_retain,&
                                             distribution_2d_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8,&
                                             sp
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum_i8,&
                                             mp_sync
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE sparse_matrix_output,            ONLY: write_matrix_distribution,&
                                             write_sparse_matrix
  USE sparse_matrix_types,             ONLY: &
       add_block_node, allocate_matrix, checksum_matrix, deallocate_matrix, &
       first_block_node, get_block_node, next_block_node, &
       real_block_node_type, real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_operations'

  PUBLIC :: dbcsr_sm_sm_multiply, dbcsr_check, dbcsr_multiply_local
  PUBLIC :: cp_dbcsr_add_dist2d, cp_dbcsr_get_dist2d, cp_dbcsr_remove_dist2d
  PUBLIC :: dbcsr_multiply_vec, dbcsr_pack_vec,&
       packed_vec_scale, dbcsr_unpack_vec,&
       add_packed_vecs, dbcsr_mult_pack_vec, packed_vec_bif_tech, &
       packed_vec_bif_tech2, dbcsr_copy_vec, packed_vec_build_u
  ! CP2K API emulation
  PUBLIC :: dbcsr_add_block_node,dbcsr_get_block_node,&
            dbcsr_replicate_matrix_struct,&
            dbcsr_allocate_matrix, dbcsr_deallocate_matrix,&
            dbcsr_allocate_matrix_set, dbcsr_deallocate_matrix_set,&
            dbcsr_obj_p_type, dbcsr_block_node_type, test_sm,&
            dbcsr_from_fm, dbcsr_from_sm, &!copy_fm_to_dbcsr, copy_dbcsr_to_fm,&
            sm_from_dbcsr, dbcsr_type_from_sm_symmetry, dbcsr_write_sparse_matrix, &
            dbcsr_write_matrix_distribution, dbcsr_multiply_sm_fm

  INTERFACE dbcsr_allocate_matrix_set
     MODULE PROCEDURE allocate_dbcsr_matrix_set, allocate_dbcsr_matrix_set_2d
  END INTERFACE

  INTERFACE dbcsr_deallocate_matrix_set
     MODULE PROCEDURE deallocate_dbcsr_matrix_set,&
                      deallocate_dbcsr_matrix_set_2d
  END INTERFACE

  INTERFACE dbcsr_get_block_node
     MODULE PROCEDURE get_dbcsr_block_node, get_dbcsr_matrix_block
  END INTERFACE

  !INTERFACE dbcsr_put_block_node
  !   MODULE PROCEDURE put_block_node, put_matrix_block
  !END INTERFACE
  !
  !INTERFACE dbcsr_next_block_node
  !   MODULE PROCEDURE next_block_node
  !END INTERFACE

!> \brief Emulation of real_matrix_p_type
  TYPE dbcsr_obj_p_type
     TYPE(dbcsr_obj), POINTER :: matrix
  END TYPE dbcsr_obj_p_type

  PRIVATE

#define DBG IF (dbg) WRITE(*,*)routineN//" ",
#define DBGV IF (bcsr_verbose) WRITE(*,*)routineN//" ",
#define DBGI IF (info) WRITE(*,*)routineN//" ",
#define DEBUG_HEADER        LOGICAL :: dbg, info
#define DEBUG_BODY        dbg = .FALSE. ; CALL dbcsr_set_debug(dbg, info=info)
#define DEBUG_BODYY       dbg = .TRUE. ; CALL dbcsr_set_debug(dbg, info=info)


CONTAINS

! *****************************************************************************
!> \brief Performs a multiplication of two dbcsr_type matrices.
!> \par Matrices m_a and m_b are multiplied into the m_c product matrix. If the
!>      dist2d parameter is not specified, then a new distribution_2d is
!>      determined for it.
!> \param[in] matrix_a        left BCSR matrix
!> \param[in] matrix_b        right BCSR matrix
!> \param[out] matrix_c       resulting BCSR product matrix.
!> \param[in] dist            distribution for the product matrix.
!> \param[in] hint            code path selection
!> \param error               cp2k error
!> \param[out] flop           (optional) effective flop
! *****************************************************************************
  SUBROUTINE dbcsr_sm_sm_multiply(matrix_a, matrix_b, matrix_c, distribution, &
       hint, error, flop)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_c
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: distribution
    INTEGER, INTENT(IN)                      :: hint
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sm_sm_multiply', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, timing_handle
    INTEGER(KIND=int_8)                      :: my_flop
    REAL(KIND=dp)                            :: flops
    TYPE(dbcsr_array_type)                   :: ms_left
    TYPE(dbcsr_imagedistribution_type)       :: rdist_left
    TYPE(dbcsr_obj)                          :: matrix_d

!TYPE(dbcsr_distribution_obj)             :: dist_left, dist_right

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    !> pre The matrices must be valid and they must have compatible shapes.
    ! Assign new dist2d
    ! Create & initialize the final & working matrix.
    IF (dbg) WRITE(*,*)routineP//'Creating 2-d distribution...'
    IF (dbg) CALL mp_sync(dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%m%dist)))
    CALL dbcsr_create(matrix_c, "C", distribution, 'N',&
         matrix_a%m%row_blk_size, matrix_b%m%col_blk_size,&
         0, 0, error=error)
    IF (dbg) WRITE(*,*)routineP//'Updating new matrix...'
    IF (dbg) CALL mp_sync(dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%m%dist)))
    !IF (dbg) WRITE(*,*)routineP//' C initial index is', matrix_c%m%index
    IF (dbg) WRITE(*,*)routineP//'Creating work matrix...'
    IF (dbg) CALL mp_sync(dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%m%dist)))
    i = dbcsr_nfullrows_local(matrix_a) * dbcsr_nfullcols_local(matrix_b)
    !CALL dbcsr_work_create(matrix_c, sizedata_guess=i, error=error)
    IF (dbg) WRITE(*,*)routineP//'Multiplying...'
    IF (dbg) CALL mp_sync(dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%m%dist)))
    !
    DBG 'Creating left image distribution'
    !CALL dbcsr_create_image_dist_left(rdist_left, matrix_a%m%dist, matrix_b%m%dist)
    CALL dbcsr_create_col_image_dist (rdist_left,&
         dbcsr_distribution_row_dist (matrix_a%m%dist), matrix_b%m%dist)
    !
    DBG 'Making A images'
    CALL dbcsr_make_col_images(matrix_a, ms_left, rdist_left, error)
    !ms_left%image_dist = rdist_left
    matrix_d = ms_left%mats(1)
    DBG 'master name', ms_left%mats(1)%m%name
    DBG "Verifying A images"
    DO i = 1, SIZE (ms_left%mats)
       CALL dbcsr_verify_matrix(ms_left%mats(i), error)
    ENDDO
    IF (PRESENT (flop)) flop = 0
    flops = dbcsr_checksum(matrix_d, error)
    DBG "DBCSR norm. S checksum", flops
    !CALL dbcsr_print(matrix_e, error)
    DBG 'Multiplying matrices'
    !
    my_flop = 0
    IF (hint .GE. 0) THEN
       !CALL dbcsr_multiply_from_left(ms_left, matrix_b, matrix_c, error,&
       !     flop=my_flop)
       CALL dbcsr_mult_e_ec_m (ms_left, matrix_b, matrix_c, error,&
            flop=my_flop)
       IF (PRESENT (flop)) flop = my_flop
    ELSE
       IF (PRESENT (flop)) flop = 0
    ENDIF
    !
    !IF (dbg) THEN
    !   CALL dbcsr_verify_matrix(matrix_c, error)
    !   checksum = dbcsr_checksum(matrix_c, error)
    !   WRITE(*,*)routineP//' BCSR Sx Matrix checksum is', checksum
    !ENDIF
    !CALL mp_sum(flop, dbcsr_mp_group (dbcsr_distribution_mp (matrix_a%m%dist)))
    !WRITE(*,*)'FLOP=',flop
    !
    !CALL dbcsr_distribution_release (dist_left)
    !CALL dbcsr_distribution_release (dist_right)
    CALL dbcsr_destroy_image_dist (rdist_left)
    CALL dbcsr_destroy_array (ms_left)
    ! Perform multiplication
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_sm_sm_multiply



! *****************************************************************************
!> \brief   Checks the new bcsr routines, to find bugs early
!> \author  Urban
!> \date    2008-10-07
! *****************************************************************************
  SUBROUTINE dbcsr_check(matrix_s, matrix_x, hint, error)

    TYPE(real_matrix_type), POINTER          :: matrix_s
    TYPE(cp_fm_type), POINTER                :: matrix_x
    INTEGER, INTENT(IN)                      :: hint
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col, colsize, handle, k, &
                                                output_unit, row, rowsize
    INTEGER(KIND=int_8)                      :: flop
    INTEGER, DIMENSION(:), POINTER           :: cbs, rbs
    LOGICAL                                  :: failure, found, ionode, ldbg, &
                                                tr
    REAL(kind=dp)                            :: checksum, newtime, oldtime, &
                                                trace, walltime1, walltime2
    REAL(KIND=dp), DIMENSION(:), POINTER     :: block1d
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block2d, block2d_2
    TYPE(array_i1d_obj)                      :: col_blk_size, &
                                                col_blk_size_right, &
                                                row_blk_size
    TYPE(cp_fm_type), POINTER                :: fmout
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_distribution_obj)             :: dist_right, product_dist, &
                                                weird_dist
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_obj)                          :: matrix_a, matrix_b, &
                                                matrix_bt, matrix_c
    TYPE(dbcsr_obj), POINTER                 :: matrix_a2
    TYPE(distribution_2d_type), POINTER      :: dist2d
    TYPE(real_matrix_type), POINTER          :: matrix_s2

    CALL timeset(routineN,handle)

    failure=.FALSE.
    ldbg = .TRUE.
    logger => cp_error_get_logger (error)
    ionode = logger%para_env%mepos == logger%para_env%source
    output_unit = 0
    IF (ionode) output_unit = cp_logger_get_default_unit_nr (logger)
    !Do the new multiplication
    CALL dbcsr_init (matrix_a, error)
    CALL dbcsr_from_sm(matrix_a, matrix_s, error)
    CALL dbcsr_verify_matrix(matrix_a, error)
    !WRITE(*,*)'S checksum=',&
    !     checksum_matrix(matrix_s, matrix_s%distribution_2d%blacs_env%para_env)
    IF (ldbg) THEN
       checksum = checksum_matrix(&
            matrix_s, matrix_s%distribution_2d%blacs_env%para_env)
       IF (ldbg.AND.ionode) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"   SM S  checksum", checksum
       ENDIF
    ENDIF
    checksum = 0.0_dp
    checksum = dbcsr_checksum(matrix_a, error)
    !WRITE(*,*)'BCSR S checksum',checksum
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR S  checksum", checksum
    ENDIF
    !
    ! Convert back to see if we get the same thing.
    dist2d => matrix_s%distribution_2d
    NULLIFY (matrix_s2)
    CALL sm_from_dbcsr (matrix_s2, matrix_a, dist2d, error)
    checksum = checksum_matrix(matrix_s2,&
                               matrix_s2%distribution_2d%blacs_env%para_env)
    !WRITE(*,*)'BCSR S checksum',checksum
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"  SM2 S  checksum", checksum
    ENDIF
    ! Check that pointer reuse works.
    CALL sm_from_dbcsr (matrix_s2, matrix_a, dist2d, error)
    checksum = checksum_matrix(matrix_s2,&
                               matrix_s2%distribution_2d%blacs_env%para_env)
    !WRITE(*,*)'BCSR S checksum',checksum
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"  SM2 S2 checksum", checksum
    ENDIF
    
    CALL deallocate_matrix (matrix_s2, error)
    !
    ! test the scale and trace routines on a square matrix
    CALL dbcsr_scale(matrix_a,alpha_scalar=2.0_dp,error=error)
    checksum = 0.0_dp
    checksum = dbcsr_checksum(matrix_a, error)
    CALL dbcsr_trace(matrix_a,trace,error=error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,2(A,F12.4))')"DBCSR 2*S  checksum", checksum, " trace(2*S)",trace
    ENDIF
    CALL dbcsr_scale(matrix_a,alpha_scalar=0.5_dp,error=error)
    !
    ! test the copy and add
    ALLOCATE (matrix_a2)
    CALL dbcsr_init (matrix_a2)
    CALL dbcsr_copy(matrix_a2,matrix_a,error=error)
    checksum = 0.0_dp
    checksum = dbcsr_checksum(matrix_a2, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR Sbis  checksum", checksum
    ENDIF
    CALL dbcsr_add(matrix_a2,matrix_a,alpha_scalar=-2.0_dp,beta_scalar=2.0_dp,error=error)
    checksum = 0.0_dp
    checksum = dbcsr_checksum(matrix_a2, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR copy and add: 2*S-2*S  checksum", checksum
    ENDIF
    !CALL dbcsr_release(matrix_a2)
    !
    ! Replace blocks
    CALL dbcsr_replace_blocks(matrix_a2, matrix_a)
    CALL dbcsr_verify_matrix(matrix_a2, error)
    checksum = dbcsr_checksum(matrix_a2, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR replace=S checksum", checksum
    ENDIF
    CALL dbcsr_release(matrix_a2)
    !
    ! filter
    CALL dbcsr_copy(matrix_a2,matrix_a,error=error)
    CALL dbcsr_filter(matrix_a2,1.0_dp,error=error)
    checksum = dbcsr_checksum(matrix_a2, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR filter(S,1.0)  checksum", checksum
    ENDIF
    CALL dbcsr_release(matrix_a2)
    !
    ! some norms here
    CALL dbcsr_norm(matrix_a, 'MaxAbsNorm', checksum, error=error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR maxabsval(S)", checksum
    ENDIF
    CALL dbcsr_norm(matrix_a, 'Frobenius', checksum, error=error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR Frobenius(S)", checksum
    ENDIF
    CALL dbcsr_norm(matrix_a, 'Gershgorin', checksum, error=error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR Gershgorin(S)", checksum
    ENDIF
    !
    ! Fill a matrix with all zeros.
    CALL dbcsr_init (matrix_a2, error)
    CALL dbcsr_from_sm(matrix_a2, matrix_s, error)
    ! Now delete all the blocks.
    !WRITE(*,*)'Deleting all blocks!'
    CALL dbcsr_iterator_start (iter, matrix_a2)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col, k)
       CALL dbcsr_remove_block (matrix_a2, row, col, k)
    ENDDO
    CALL dbcsr_iterator_stop (iter)
    !WRITE(*,*)"blk_p=", matrix_a2%blk_p
    CALL dbcsr_finalize (matrix_a2)
    !WRITE(*,*)"blk_p=", matrix_a2%nblks
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR 0 checksum ", checksum
    ENDIF
    !write(*,*)'zeroing...'
    !call dbcsr_work_create(matrix_a2, dbcsr_get_num_blocks (matrix_a2),&
    !     dbcsr_nfullrows_total(matrix_a2)*dbcsr_nfullcols_total(matrix_a2))
    rbs => array_data (dbcsr_row_block_sizes (matrix_a2))
    cbs => array_data (dbcsr_col_block_sizes (matrix_a2))
    ALLOCATE (block2d(MAXVAL (rbs), MAXVAL (cbs)))
    block2d(:,:) = 0
    ALLOCATE (block1d(MAXVAL(rbs)*MAXVAL(cbs)))
    block1d(:) = 0
    k = 0
    DO row = 1, MIN(10, dbcsr_nblkrows_total (matrix_a2))
       rowsize = rbs (row)
       DO col = 1, MIN(10, dbcsr_nblkcols_total (matrix_a2))
          colsize = cbs (col)
          CALL dbcsr_put_block (matrix_a2, row, col,&
               block2d(1:rowsize, 1:colsize))
       ENDDO
    ENDDO
    CALL dbcsr_finalize (matrix_a2)
    checksum = dbcsr_checksum(matrix_a2, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR 0 checksum2", checksum
    ENDIF
    DEALLOCATE (block2d)
    DEALLOCATE (block1d)
    CALL dbcsr_release (matrix_a2)
    !
    ! Tests (undefined) block addition then subsequently setting the
    ! data area.
    DEALLOCATE (matrix_a2)
    CALL dbcsr_allocate_matrix (matrix_a2, matrix_s%nrow, matrix_s%ncol,&
         matrix_s%nblock_row, matrix_s%nblock_col,&
         matrix_s%first_row, matrix_s%last_row,&
         matrix_s%first_col, matrix_s%last_col,&
         matrix_s%name, matrix_s%symmetry,&
         0, matrix_s%distribution_2d, error=error)
    !write(*,*)'zeroing..., second time around'
    rbs => array_data (dbcsr_row_block_sizes (matrix_a2))
    cbs => array_data (dbcsr_col_block_sizes (matrix_a2))
    DO row = 1, dbcsr_nblkrows_total (matrix_a2)
       rowsize = rbs (row)
       DO col = row, dbcsr_nblkcols_total (matrix_a2)
          colsize = cbs (col)
          NULLIFY (block2d)
          !CALL dbcsr_reserve_block2d (matrix_a2, row, col,&
          !     block2d)
          CALL dbcsr_add_block_node (matrix_a2, row, col, block2d, error)
!write(*,*)routineN//" setting block at virt. coor", row, col, rowsize, colsize          
          !block2d(1:rowsize, 1:colsize) = 1.0_dp
          block2d(:,:) = 1.0_dp
       ENDDO
    ENDDO
    DO row = 1, dbcsr_nblkrows_total (matrix_a2)
       DO col = 1, dbcsr_nblkcols_total (matrix_a2)
          CALL dbcsr_get_tmp_block2d (matrix_a2, row, col, block2d_2, tr, found)
          IF (found) THEN
             block2d_2(:,:) = 0.0_dp
          ENDIF
          !write(*,*)'size1:',size(block2d_2), rowsize, colsize
          !write(*,'(10(1X,F5.1))')block2d_2
       ENDDO
    ENDDO
    DO row = 1, dbcsr_nblkrows_total (matrix_a2)
       DO col = 1, dbcsr_nblkcols_total (matrix_a2)
          CALL dbcsr_get_tmp_block2d (matrix_a2, row, col, block2d, tr, found)
          !write(*,*)'size3:',size(block2d), rowsize, colsize
          !write(*,'(10(1X,F5.1))')block2d
          !IF (.NOT. found) WRITE(*,*)"Oops, missing block."
       ENDDO
    ENDDO
    CALL dbcsr_finalize (matrix_a2, error=error)
    CALL dbcsr_verify_matrix(matrix_a2, error)
    checksum = dbcsr_checksum(matrix_a2, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR 0 checksum4", checksum
    ENDIF
    !
    ! Tests block getting
    DO row = 1, dbcsr_nblkrows_total (matrix_a2)
       DO col = 1, dbcsr_nblkcols_total (matrix_a2)
          CALL dbcsr_get_block_node (matrix_a2, row, col, block=block2d)
          !WRITE(*,*)'getting block_node:',row, col, ASSOCIATED(block2d)
       ENDDO
    ENDDO
    !
    !CALL dbcsr_release (matrix_a2)
    CALL dbcsr_deallocate_matrix (matrix_a2, error)
    !
    ! Tests matrix structure replication
    CALL dbcsr_replicate_matrix_struct(matrix_a, matrix_a2, "replication",&
         error=error)
    CALL dbcsr_finalize (matrix_a2, error=error)
    CALL dbcsr_verify_matrix(matrix_a2, error)
    checksum = dbcsr_checksum(matrix_a2, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR 0 checksum5", checksum
    ENDIF
    !
    CALL dbcsr_deallocate_matrix (matrix_a2, error)
    

    checksum = 0.0_dp
    CALL cp_fm_frobenius_norm(matrix_x, checksum, error)
    !WRITE(*,*)'x checksum',checksum*checksum
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"   FM x  checksum", checksum*checksum
    ENDIF
    CALL cp_fm_get_info(matrix_x, ncol_global=k, error=error)
    
    IF (MOD(hint,2) .EQ. 0) THEN
       CALL dbcsr_create_dist_r_unrot (dist_right, matrix_a%m%dist, k,&
            col_blk_size_right)
    ELSE
       CALL dbcsr_create_dist_r_rot (dist_right, matrix_a%m%dist, k,&
            col_blk_size_right)
    ENDIF
    CALL dbcsr_init (matrix_b, error)
    CALL dbcsr_from_fm(matrix_b, matrix_x, 0.0_dp, dist_right,&
         row_blk_size=matrix_a%m%row_blk_size,&
         col_blk_size=col_blk_size_right,&
         error=error)
    CALL dbcsr_verify_matrix(matrix_b, error)
    !CALL dbcsr_print (matrix_b, error)

    checksum = dbcsr_checksum(matrix_b, error)
    !WRITE(*,*)'DBCSR x checksum',checksum
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x  checksum", checksum
    ENDIF
    !
    ! Test replacing a block.
    !CALL dbcsr_print (matrix_b, error)
    CALL dbcsr_iterator_start (iter, matrix_b)
    IF (dbcsr_iterator_blocks_left (iter)) THEN
       checksum = dbcsr_checksum(matrix_b, error, local=.TRUE.)
       IF (ldbg.AND.ionode) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x checksumL", checksum
       ENDIF
       CALL dbcsr_iterator_next_block (iter, row, col, block2d, tr)
       !write(*,*)'got block',row,col,tr
       CALL dbcsr_put_block (matrix_b, row, col, block2d, tr)
       CALL dbcsr_verify_matrix (matrix_b, error, local=.TRUE.)
       !CALL dbcsr_print (matrix_b, error)
       checksum = dbcsr_checksum(matrix_b, error, local=.TRUE.)
       IF (ldbg.AND.ionode) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x checksum2", checksum
       ENDIF
       ALLOCATE (block2d_2(SIZE(block2d,1),SIZE(block2d,2)))
       block2d_2(:,:) = block2d(:,:)
       !write(*,'(8(1X,F7.2))')block2d_2
       block2d(:,:) = 0.0_dp
       CALL dbcsr_put_block (matrix_b, row, col, block2d, tr)
       CALL dbcsr_finalize (matrix_b, error=error)
       CALL dbcsr_verify_matrix (matrix_b, error, local=.TRUE.)
       checksum = dbcsr_checksum(matrix_b, error, local=.TRUE.)
       IF (ldbg.AND.ionode) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x wrng cksm", checksum
       ENDIF
       !CALL dbcsr_print (matrix_b, error)
       CALL dbcsr_put_block (matrix_b, row, col, block2d_2, tr)
       CALL dbcsr_finalize (matrix_b, error=error)
       CALL dbcsr_verify_matrix (matrix_b, error, local=.TRUE.)
       checksum = dbcsr_checksum(matrix_b, error, local=.TRUE.)
       IF (ldbg.AND.ionode) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x checksum3", checksum
       ENDIF
       !write(*,*)'removing block'
       CALL dbcsr_remove_block (matrix_b, row, col)
       CALL dbcsr_finalize (matrix_b, error=error)
       CALL dbcsr_verify_matrix (matrix_b, error=error, local=.TRUE.)
       checksum = dbcsr_checksum(matrix_b, error, local=.TRUE.)
       IF (ldbg.AND.ionode) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x wrng cks4", checksum
       ENDIF
       !CALL dbcsr_print (matrix_b, error)
       !write(*,*)'use block',row,col,tr
       !write(*,'(8(1X,F7.2))')block2d_2
       CALL dbcsr_put_block (matrix_b, row, col, block2d_2, tr)
       CALL dbcsr_finalize (matrix_b, error=error)
       CALL dbcsr_verify_matrix (matrix_b, error=error, local=.TRUE.)
       !CALL dbcsr_print (matrix_b, error)
       checksum = dbcsr_checksum(matrix_b, error, local=.TRUE.)
       IF (ldbg.AND.ionode) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x checksm4L", checksum
       ENDIF
       DEALLOCATE (block2d_2)
    ENDIF
    CALL dbcsr_iterator_stop (iter)
    CALL dbcsr_verify_matrix (matrix_b, error=error)
    checksum = dbcsr_checksum(matrix_b, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x checksm4g", checksum
    ENDIF
    CALL dbcsr_finalize (matrix_b)
    CALL dbcsr_verify_matrix (matrix_b, error=error)
    checksum = dbcsr_checksum(matrix_b, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x checksum4", checksum
    ENDIF
    !
    ! test the scale routine on a general matrix
    CALL dbcsr_scale(matrix_b,alpha_scalar=3.0_dp,error=error)
    checksum = dbcsr_checksum(matrix_b, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR 3*x  checksum", checksum
    ENDIF
    CALL dbcsr_scale(matrix_b,alpha_scalar=1.0_dp/3.0_dp,error=error)
    checksum = dbcsr_checksum(matrix_b, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR 3*x/3  checksum", checksum
    ENDIF
    !
    ! Replication
    CALL dbcsr_init (matrix_bt)
    CALL dbcsr_copy (matrix_bt, matrix_b, error=error)
    CALL dbcsr_replicate_all (matrix_bt)
    CALL dbcsr_verify_matrix (matrix_bt)
    checksum = dbcsr_checksum(matrix_bt, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x.*p checksum", checksum
    ENDIF
    ! Sum the replicates.
    CALL dbcsr_sum_replicated (matrix_bt)
    CALL dbcsr_verify_matrix (matrix_bt)
    checksum = dbcsr_checksum(matrix_bt, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSRpx.*p checksum", checksum
    ENDIF
    ! Again, to test unity of distribution of replication
    !CALL dbcsr_release(matrix_bt)
    CALL dbcsr_copy (matrix_bt, matrix_b, error=error)
    CALL dbcsr_replicate_all (matrix_bt)
    CALL dbcsr_verify_matrix (matrix_bt)
    CALL dbcsr_distribute (matrix_bt)
    CALL dbcsr_verify_matrix (matrix_bt)
    checksum = dbcsr_checksum(matrix_bt, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x s  checksum", checksum
    ENDIF
    !CALL dbcsr_release(matrix_bt)
    !
    CALL dbcsr_copy (matrix_bt, matrix_b, error=error)
    CALL dbcsr_replicate_all (matrix_bt)
    CALL dbcsr_verify_matrix (matrix_bt)
    CALL dbcsr_distribute (matrix_bt, fast=.TRUE.)
    CALL dbcsr_verify_matrix (matrix_bt)
    !checksum = dbcsr_checksum(matrix_bt, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x d  checksum", checksum
    ENDIF

    !CALL dbcsr_release (matrix_bt)
    ! Transpose
    CALL dbcsr_new_transposed (matrix_bt, matrix_b, error=error)
    CALL dbcsr_verify_matrix (matrix_bt)
    !
    CALL dbcsr_distribution_new (product_dist,&
         dbcsr_distribution_mp (matrix_a%m%dist),&
         dbcsr_distribution_row_dist (matrix_a%m%dist),&
         dbcsr_distribution_col_dist (matrix_b%m%dist))
    !
    walltime1 = m_walltime()
    flop = 0
    CALL dbcsr_init (matrix_c, error=error)
    CALL dbcsr_sm_sm_multiply(matrix_a, matrix_b, matrix_c,&
         product_dist, hint, error, flop=flop)
    !CALL dbcsr_sm_sm_multiply(matrix_bt, matrix_a, matrix_c,&
    !     dist_result, hint, error)
    !CALL dbcsr_distribution_release (dist_result)
    walltime2 = m_walltime()
    newtime = walltime2-walltime1
    IF (newtime .LT. 0.1_dp) newtime = 0.1_dp
    !
    CALL dbcsr_verify_matrix(matrix_c, error)
    !
    checksum = dbcsr_checksum(matrix_c, error)
    !WRITE(*,*)'BCSR Sx checksum',checksum,'time=',newtime
    CALL mp_sum_i8 (flop, dbcsr_mp_group(dbcsr_distribution_mp(&
         dbcsr_distribution(matrix_c))))
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4,A,F12.4,";",F12.4,A)')&
            "DBCSR Sx checksum", checksum,&
            " time ", newtime,&
            REAL(flop, dp)/1000000000.0_dp/newtime, " GFLOP/s"
       !WRITE(output_unit,'(A,I30,5X,EN12.4)')'FLOP=',flop, REAL (flop, dp)
    ENDIF
    !
    ! Verify copy_dbcsr_to_fm
    CALL dbcsr_copy(matrix_bt, matrix_b, error=error)
    CALL dbcsr_iterator_start(iter,matrix_bt)
    DO WHILE(dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter,row,col,k)
       !write(*,*)'i have',row,col
       IF (MOD(row+col,3) .EQ. 0) &
            CALL dbcsr_remove_block(matrix_bt,row,col,k)
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    !call dbcsr_remove_block(matrix_bt, 1, 1)
    !call dbcsr_finalize(matrix_bt)
    CALL dbcsr_finalize(matrix_bt)
    CALL dbcsr_verify_matrix(matrix_bt)
    CALL copy_fm_to_dbcsr(matrix_x, matrix_bt)
    CALL dbcsr_verify_matrix (matrix_bt)
    checksum = dbcsr_checksum(matrix_bt, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"copy fm2dbcsr chksm", checksum
    ENDIF

    CALL dbcsr_release (matrix_bt)
    ! check trace(c'*s*c)
    !call dbcsr_replicate_all (matrix_b)
    !call dbcsr_replicate_all (matrix_c)
    CALL dbcsr_trace(matrix_b, matrix_c, trace, error=error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4,A,I4)')"DBCSR trace(C'*S*C)", trace," should be equal to ",k
    ENDIF

    CALL dbcsr_release (matrix_c)

    ! Test redistribution
    CALL dbcsr_create_dist_elemental (weird_dist,&
         dbcsr_nfullrows_total(matrix_b),&
         dbcsr_nfullcols_total(matrix_b),&
         dbcsr_distribution_mp(dist_right),&
         row_blk_size, col_blk_size)
    CALL dbcsr_create (matrix_bt, "Elemental", weird_dist, 'N',&
         row_blk_size, col_blk_size, 0, 0)
    CALL dbcsr_distribution_release (weird_dist)
    CALL array_release (row_blk_size)
    CALL array_release (col_blk_size)

    CALL dbcsr_complete_redistribute (matrix_b, matrix_bt)
    CALL dbcsr_verify_matrix (matrix_bt)
    checksum = dbcsr_checksum(matrix_bt, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"Elemental B chcksum", checksum
    ENDIF
    !
    CALL dbcsr_create_dist_fullest (weird_dist,&
         dbcsr_nfullrows_total(matrix_b),&
         dbcsr_nfullcols_total(matrix_b),&
         dbcsr_distribution_mp(dist_right),&
         row_blk_size, col_blk_size)
    CALL dbcsr_create (matrix_bt, "Fullest", weird_dist, 'N',&
         row_blk_size, col_blk_size, 0, 0)
    CALL dbcsr_distribution_release (weird_dist)
    CALL array_release (row_blk_size)
    CALL array_release (col_blk_size)
    CALL dbcsr_complete_redistribute (matrix_b, matrix_bt)
    CALL dbcsr_verify_matrix (matrix_bt)
    checksum = dbcsr_checksum(matrix_bt, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')" Fullest B checksum", checksum
    ENDIF
    CALL dbcsr_release (matrix_b)
    CALL dbcsr_init (matrix_b)

    CALL dbcsr_release (matrix_bt)
    CALL dbcsr_init (matrix_bt)

    CALL timestop(handle)

    CALL timeset(routineN//'std', handle)

    CALL cp_fm_get_info(matrix_x, ncol_global=k, error=error)
    CALL cp_fm_create(fmout, matrix_x%matrix_struct, "FMout", error=error)
    walltime1 = m_walltime()
    CALL cp_sm_fm_multiply(matrix_s, matrix_x, fmout, k, error=error)
    walltime2 = m_walltime()
    oldtime = walltime2-walltime1
    CALL cp_fm_frobenius_norm(fmout, checksum, error)
    !WRITE(*,*)'Sx checksum',checksum*checksum,'time=',oldtime
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4,A,F12.4)')" SMFM Sx checksum",&
            checksum*checksum, " time ", oldtime
    ENDIF
    CALL dbcsr_from_fm(matrix_b, fmout, 0.0_dp, dist_right,&
         row_blk_size=matrix_a%m%row_blk_size,&
         col_blk_size=col_blk_size_right,&
         error=error)
    
    CALL dbcsr_multiply_sm_fm(matrix_a, matrix_x, fmout, k, error=error)
    CALL cp_fm_frobenius_norm(fmout, checksum, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')" dbcsr_SM_FM Sx checksum",&
            checksum*checksum
    ENDIF

    CALL array_release (col_blk_size_right)
    CALL dbcsr_release (matrix_a)
    CALL dbcsr_release (matrix_b)
    CALL cp_fm_release (fmout, error)
    CALL dbcsr_distribution_release (dist_right)
    CALL dbcsr_distribution_release (product_dist)

    CALL timestop(handle)


  END SUBROUTINE dbcsr_check

! *****************************************************************************
!> \brief   Adds distribution_2d_type to a DBCSR matrix
!> \author  Urban
!> \date    2009-05-27
! *****************************************************************************
  SUBROUTINE cp_dbcsr_add_dist2d (matrix, distribution_2d, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(cp_error_type), INTENT(INOUT)       :: error

!   ---------------------------------------------------------------------------

    CALL distribution_2d_retain (distribution_2d, error)
    matrix%dist2d => distribution_2d
  END SUBROUTINE cp_dbcsr_add_dist2d

! *****************************************************************************
!> \brief   Gets the distibution_2d of a DBCSR matrix.
!> \author  Urban
!> \date    2008-05-27
! *****************************************************************************
  SUBROUTINE cp_dbcsr_get_dist2d (matrix, distribution_2d, error)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(cp_error_type), INTENT(INOUT)       :: error

!   ---------------------------------------------------------------------------

    distribution_2d => matrix%dist2d
  END SUBROUTINE cp_dbcsr_get_dist2d

! *****************************************************************************
!> \brief   Removes distibution_2d from a DBCSR matrix.
!> \author  Urban
!> \date    2008-05-27
! *****************************************************************************
  SUBROUTINE cp_dbcsr_remove_dist2d (matrix, error)
    TYPE(dbcsr_type), INTENT(INOUT)          :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

!   ---------------------------------------------------------------------------

    CALL distribution_2d_release (matrix%dist2d, error)
    NULLIFY (matrix%dist2d)
  END SUBROUTINE cp_dbcsr_remove_dist2d

! *****************************************************************************
!> \brief multiply a dbcsr with a column vector
!>        c (packed column) = A (dbscr) * B_b (dbcsr)
!> \param[in]  matrix_a  DBCSR matrix
!> \param[in]  matrix_b  DBCSR matrix
!> \param[in]  b_col     the column of the matrix B
!> \param[out] pkd       the result in a packed format
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_vec(matrix_a, matrix_b, a_row_beg, a_row_end, b_col, pkd, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a, matrix_b
    INTEGER, INTENT(in)                      :: a_row_beg, a_row_end, b_col
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pkd
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_vec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: timing_handle
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: pkd_b

    CALL timeset(routineN, timing_handle)

    ALLOCATE(pkd_b(SIZE(pkd)))
    !
    ! packed the b_col
    CALL dbcsr_pack_vec(matrix_b, b_col, pkd_b, 'column', error)
    !
    ! send the packed col to the right guy
    !CALL packed_vec_bcast(pkd_b, source, , mp_obj, error)
    !
    ! local multiply
    CALL dbcsr_mult_pack_vec(matrix_a, pkd_b, a_row_beg, a_row_end, b_col, pkd, error)
    !
    ! redistribute the result
    !CALL ...

    DEALLOCATE(pkd_b)

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_multiply_vec

! *****************************************************************************
!> \brief multiply a dbcsr with a packed vector as
!>        c (packed column) = A (dbscr) * b (packed column)
!> \param[in]  matrix      DBCSR matrix
!> \param[in]  pkd_b       packed vector b
!> \param[in]  b_col       the column of the packed b (this should be changed in b_col_size)
!> \param[out] pkd_c       packed vector c
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_mult_pack_vec(matrix_a, pkd_b, a_row_beg, a_row_end, b_col, pkd_c, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    REAL(dp), DIMENSION(:), INTENT(in)       :: pkd_b
    INTEGER, INTENT(in)                      :: a_row_beg, a_row_end, b_col
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pkd_c
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mult_pack_vec', &
      routineP = moduleN//':'//routineN

    INTEGER :: a_blk, a_col, a_col_size, a_offset, a_row, a_row_size, &
      b_col_size, b_offset, c_offset, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: a_col_blk_size, a_row_blk_size
    LOGICAL                                  :: new_c
    REAL(dp), DIMENSION(:), POINTER          :: data_p

    CALL timeset(routineN, timing_handle)

    a_row_blk_size => array_data (matrix_a%m%row_blk_size)
    a_col_blk_size => array_data (matrix_a%m%col_blk_size)


    IF(pkd_b(matrix_a%m%nblkrows_total+1).LE.0) CALL stop_program(routineN,"1pkd_b(n+1).LE.0")

    pkd_c(:) = 0.0_dp ! if not done
    c_offset = matrix_a%m%nblkrows_total + 2
    b_col_size = a_col_blk_size(b_col)

    CALL dbcsr_get_data (matrix_a, data_p)
    DO a_row = a_row_beg,a_row_end
       a_row_size = a_row_blk_size(a_row)
       new_c = .FALSE.
       DO a_blk = matrix_a%m%row_p(a_row)+1,matrix_a%m%row_p(a_row+1)
          a_col = matrix_a%m%col_i(a_blk)
          a_offset = ABS(matrix_a%m%blk_p(a_blk))
          a_col_size = a_col_blk_size(a_col)
          b_offset = INT(pkd_b(a_col))
          IF(b_offset.LE.0) CYCLE
          !
          pkd_c(a_row) = REAL(c_offset,dp)
          CALL dgemm('N','N',a_row_size,b_col_size,a_col_size,&
               &     1.0_dp,data_p(a_offset),a_row_size,&
               &            pkd_b(b_offset),a_col_size,&
               &     1.0_dp,pkd_c(c_offset),a_row_size)
          new_c = .TRUE.
       ENDDO ! a_col
       IF(new_c) c_offset = c_offset + a_row_size*b_col_size
    ENDDO ! a_row
    pkd_c(matrix_a%m%nblkrows_total+1) = REAL(c_offset,dp)

    IF(pkd_b(matrix_a%m%nblkrows_total+1).LE.0) CALL stop_program(routineN,"2pkd_b(n+1).LE.0")
    IF(pkd_c(matrix_a%m%nblkrows_total+1).LE.0) CALL stop_program(routineN,"2pkd_c(n+1).LE.0")

    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_mult_pack_vec





! *****************************************************************************
!> \brief multiply a dbcsr with a replicated array
!>        c = alpha_scalar * A (dbscr) * b + c 
!> \param[in]  matrix       DBCSR matrix
!> \param[in]  vec_b        vectors b
!> \param[inout] vec_c      vectors c
!> \param[in]  ncol         nbr of columns
!> \param[in]  alpha        alpha
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_local(matrix_a, vec_b, vec_c, ncol, alpha, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_a
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: vec_b
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: vec_c
    INTEGER, INTENT(in), OPTIONAL            :: ncol
    REAL(dp), INTENT(IN), OPTIONAL           :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, my_ncol, row, &
                                                timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: offset
    INTEGER, DIMENSION(:), POINTER           :: row_blk_size
    LOGICAL                                  :: tr
    REAL(dp)                                 :: my_alpha, my_alpha2
    REAL(dp), DIMENSION(:, :), POINTER       :: data_d
    TYPE(dbcsr_iterator)                     :: iter

    CALL timeset(routineN, timing_handle)
    
    my_alpha = 1.0_dp
    IF (PRESENT(alpha)) my_alpha = alpha

    my_ncol = SIZE(vec_b,2)
    IF(PRESENT(ncol)) my_ncol = ncol

    my_alpha2 = 0.0_dp
    IF(dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric) my_alpha2 = my_alpha
    IF(dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric) my_alpha2 = -my_alpha

    !> this should be precomputed somewhere else
    row_blk_size => array_data (matrix_a%m%row_blk_size)
    ALLOCATE(offset(matrix_a%m%nblkrows_total))
    offset(1) = 1
    DO row = 2,matrix_a%m%nblkrows_total
       offset(row) = offset(row-1) + row_blk_size(row-1)
    ENDDO
    !<

    CALL dbcsr_iterator_start(iter, matrix_a)

    DO WHILE (dbcsr_iterator_blocks_left(iter))

       CALL dbcsr_iterator_next_block(iter, row, col, data_d, tr, blk)

       CALL dgemm('N','N',SIZE(data_d,1),my_ncol,SIZE(data_d,2),&
                  my_alpha,data_d(1,1),SIZE(data_d,1),&
                  vec_b(offset(col),1),SIZE(vec_b,1),&
                  1.0_dp,vec_c(offset(row),1),SIZE(vec_c,1))

       IF((dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_symmetric.OR.&
          dbcsr_get_matrix_type(matrix_a).EQ.dbcsr_type_antisymmetric)) THEN
          IF(row.NE.col) THEN
             CALL dgemm('T','N',SIZE(data_d,1),my_ncol,SIZE(data_d,2),&
                  my_alpha2,data_d(1,1),SIZE(data_d,1),&
                  vec_b(offset(row),1),SIZE(vec_b,1),&
                  1.0_dp,vec_c(offset(col),1),SIZE(vec_c,1))
          ENDIF
       ENDIF
    ENDDO

    CALL dbcsr_iterator_stop(iter)

    DEALLOCATE(offset)

    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_multiply_local

! *****************************************************************************
!> \brief multiply a dbcsr with a fm matrix
!> \param[in]  matrix         DBCSR matrix
!> \param[]    fm_in, fm_out  fm matrices
!> \param[in]  ncol           nbr of columns
!> \param[in]  alpha          alpha
!> \param[in]  beta           beta
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_multiply_sm_fm(matrix, fm_in, fm_out, ncol, alpha, beta, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(cp_fm_type), POINTER                :: fm_in, fm_out
    INTEGER, INTENT(IN)                      :: ncol
    REAL(dp), INTENT(IN), OPTIONAL           :: alpha, beta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_multiply_sm_fm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k, timing_handle
    TYPE(dbcsr_obj)                          :: in, out
    TYPE(dbcsr_distribution_obj)             :: dist_right, product_dist
    TYPE(array_i1d_obj)                      :: col_blk_size_right
    REAL(dp)                                 :: my_alpha, my_beta

    CALL timeset(routineN, timing_handle)
    
    my_alpha=1.0_dp
    my_beta=0.0_dp
    IF (PRESENT(alpha)) my_alpha=alpha
    IF (PRESENT(beta)) my_beta=beta
    IF(my_alpha.NE.1.0_dp.OR.my_beta.NE.0.0_dp) &
         CALL stop_program(routineN,"NYI alpha, beta") 

    CALL cp_fm_get_info(fm_in, ncol_global=k, error=error)
    IF(k.NE.ncol) CALL stop_program(routineN,"NYI ncol")

    CALL dbcsr_create_dist_r_unrot (dist_right, matrix%m%dist, ncol, col_blk_size_right)
    
    CALL dbcsr_init(in,error)
    CALL dbcsr_create(in, "D", dist_right, 'N', &
         matrix%m%row_blk_size, col_blk_size_right,&
         0, 0, error=error)
    CALL copy_fm_to_dbcsr(fm_in, in)

    CALL dbcsr_init(out,error)
    CALL dbcsr_create(out, "D", dist_right, 'N', &
         matrix%m%row_blk_size, col_blk_size_right,&
         0, 0, error=error)
    CALL copy_fm_to_dbcsr(fm_out, out)

    CALL dbcsr_set(out,0.0_dp,error=error)

    CALL dbcsr_distribution_new (product_dist,&
         dbcsr_distribution_mp (matrix%m%dist),&
         dbcsr_distribution_row_dist (matrix%m%dist),&
         dbcsr_distribution_col_dist (in%m%dist))

    !this scaling should be done in dbcsr_sm_sm_multiply
    CALL dbcsr_scale(in,alpha_scalar=my_alpha,error=error)
    CALL dbcsr_sm_sm_multiply(matrix, in, out, product_dist, 1, error)

    !CALL copy_dbcsr_to_fm(out, fm_out)

    CALL dbcsr_release(in)
    CALL dbcsr_release(out)
    CALL array_release(col_blk_size_right)
    CALL dbcsr_distribution_release(dist_right)
    CALL dbcsr_distribution_release(product_dist)

    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_multiply_sm_fm

! *****************************************************************************
!> \brief pack a row or column of a dbcsr
!> \param[in]  matrix      DBCSR matrix
!> \param[in]  ivec        the vector to pack
!> \param[out] pkd_vec     packed vector
!> \param[in]  what        what to pack (row or column)
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_pack_vec(matrix, ivec, pkd_vec, what, error)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(IN)                      :: ivec
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_vec
    CHARACTER(LEN=*), INTENT(IN)             :: what
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_pack_vec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_size, frst_blk, &
                                                last_blk, nze, row, row_size, &
                                                timing_handle, v_offset
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_p

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    row_blk_size => array_data (dbcsr_row_block_sizes (matrix))
    col_blk_size => array_data (dbcsr_col_block_sizes (matrix))

    pkd_vec(:) = 0.0_dp ! should be big enough to hold all the data
    !
    ! let's go
    CALL dbcsr_get_data (matrix, data_p)
    SELECT CASE(what)
    CASE('column')
       v_offset = matrix%m%nblkrows_total+2 ! we need nblkrows_total+1 to store rows 
       !                                    and 1 extra for adding new data
       col = ivec
       col_size = col_blk_size(col)
       DO row = 1, matrix%m%nblkrows_total
          row_size = row_blk_size(row)
          blk = matrix%m%row_p(row)+1
          frst_blk = matrix%m%row_p(row)+1
          last_blk = matrix%m%row_p(row+1)
          CALL dbcsr_find_column(col,frst_blk,last_blk,matrix%m%col_i,&
               matrix%m%blk_p,blk,found)
          IF(found) THEN
             nze = row_size*col_size
             !
             ! let's copy the block
             pkd_vec(row) = REAL(v_offset,dp)
             CALL dcopy(nze,data_p(ABS(matrix%m%blk_p(blk))),1,pkd_vec(v_offset),1)
             v_offset = v_offset + nze
          ENDIF
       ENDDO ! row
       pkd_vec(matrix%m%nblkrows_total+1) = REAL(v_offset,dp)
    CASE('row')
       v_offset = matrix%m%nblkcols_total+2 ! we need nblkcols_total+1 to store rows 
       !                                    and 1 extra for the size of the data
       row = ivec
       row_size = row_blk_size(row)
       DO blk = matrix%m%row_p(row)+1,matrix%m%row_p(row+1)
          col = matrix%m%col_i(blk)
          col_size = col_blk_size(col)
          nze = row_size*col_size
          !
          ! let's copy the block
          pkd_vec(col) = REAL(v_offset,dp)
          CALL dcopy(nze,data_p(ABS(matrix%m%blk_p(blk))),1,pkd_vec(v_offset),1)
          v_offset = v_offset + nze
       ENDDO ! col
       pkd_vec(matrix%m%nblkcols_total+1) = REAL(v_offset,dp)
    CASE DEFAULT
       CALL stop_program(routineN,"pack what?")       
    END SELECT

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_pack_vec

! *****************************************************************************
!> \brief unpack a row or column of a dbcsr
!> \param[inout]  matrix      DBCSR matrix
!> \param[in]     ivec        the vector to pack
!> \param[in]     pkd_vec     packed vector
!> \param[in]     what        what to unpack (row or column)
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_unpack_vec(matrix, ivec, pkd_vec, what, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: ivec
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_vec
    CHARACTER(LEN=*), INTENT(IN)             :: what
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_unpack_vec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, col, col_size, frst_blk, &
                                                last_blk, nze, offset, row, &
                                                row_size, timing_handle, &
                                                v_offset
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: found
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_p

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    row_blk_size => array_data (dbcsr_row_block_sizes (matrix))
    col_blk_size => array_data (dbcsr_col_block_sizes (matrix))
    CALL dbcsr_get_data (matrix, data_p)
    !
    ! let's go
    SELECT CASE(what)
    CASE('column')
       IF(pkd_vec(matrix%m%nblkrows_total+1).LE.0) CALL stop_program(routineN,"1pkd_vec(n+1).LE.0")
       col = ivec
       col_size = col_blk_size(col)
       DO row = 1, matrix%m%nblkrows_total
          v_offset = INT(pkd_vec(row))
          IF(v_offset.LE.0) CYCLE
          row_size = row_blk_size(row)
          blk = matrix%m%row_p(row)+1
          frst_blk = matrix%m%row_p(row)+1
          last_blk = matrix%m%row_p(row+1)
          CALL dbcsr_find_column(col,frst_blk,last_blk,matrix%m%col_i,&
               matrix%m%blk_p,blk,found)
          IF(found) THEN
             nze = row_size*col_size
             !
             ! let's scale the block
             offset = ABS(matrix%m%blk_p(blk))
             ! we just dcopy here, could aslo do a daxpy
             !CALL dcopy(nze,pkd_vec(v_offset),1,matrix%m%data(offset),1)
            data_p(offset:offset+nze-1)=pkd_vec(v_offset:v_offset+nze-1)
          ELSE
             !
             ! need to add it
             WRITE(*,*) routineN//' row',row,' col',col
             CALL stop_program(routineN,"some more work here")
          ENDIF
       ENDDO ! row
    CASE('row')
       IF(pkd_vec(matrix%m%nblkcols_total+1).LE.0) CALL stop_program(routineN,"2pkd_vec(n+1).LE.0")
       row = ivec
       row_size = row_blk_size(row)
       blk = matrix%m%row_p(row)+1
       frst_blk = matrix%m%row_p(row)+1
       last_blk = matrix%m%row_p(row+1)
       DO col = 1,matrix%m%nblkcols_total
          v_offset = INT(pkd_vec(col))
          IF(v_offset.LE.0) CYCLE
          col_size = col_blk_size(col)
          ! we assume that the columns are ordered !
          CALL dbcsr_find_column(col,frst_blk,last_blk,matrix%m%col_i,&
               matrix%m%blk_p,blk,found)
          IF(found) THEN
             nze = row_size*col_size
             offset = ABS(matrix%m%blk_p(blk))
             ! we just dcopy here, could aslo do a daxpy
             CALL dcopy(nze,pkd_vec(v_offset),1,data_p(offset),1)
          ELSE
             !
             ! need to add it
             CALL stop_program(routineN,"some more work here")
          ENDIF
       ENDDO ! row
    CASE DEFAULT
       CALL stop_program(routineN,"unpack what?")
    END SELECT

    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_unpack_vec


! *****************************************************************************
!> \brief Does one of the operations:
!>      left)  pkd_ji = alpha * block_jj (block) * pkd_ji for j = 1,..., nvec
!>      right) pkd_ij = alpha * pkd_ij * block_jj (block) for j = 1,..., nvec
!>
! *****************************************************************************
  SUBROUTINE packed_vec_scale(alpha, block, pkd, nvec, ivec, vec_blk_size, side, error)
    REAL(dp), INTENT(IN)                     :: alpha
    REAL(dp), DIMENSION(:), INTENT(in)       :: block
    REAL(dp), DIMENSION(:), INTENT(inout)    :: pkd
    INTEGER, INTENT(in)                      :: nvec, ivec
    INTEGER, DIMENSION(:), INTENT(in)        :: vec_blk_size
    CHARACTER(LEN=*), INTENT(IN)             :: side
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_scale', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iblk, k, m, n, offset
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff

!   ---------------------------------------------------------------------------

    ALLOCATE(buff(SIZE(pkd)))!too long buff allocated here...
    buff=0.0_dp
    SELECT CASE(side)
    CASE('left')
       iblk = 1
       DO i = 1,nvec
          offset = INT(pkd(i))
          IF(offset.GT.0) THEN
             m = vec_blk_size(i)
             n = vec_blk_size(ivec)
             k = vec_blk_size(i)
             CALL dgemm('N','N',m,n,k,alpha,block(iblk),m,pkd(offset),k,0.0_dp,buff(1),m)
             CALL dcopy(m*n,buff(1),1,pkd(offset),1)
          ENDIF
          iblk = iblk + vec_blk_size(i)**2
       ENDDO
    CASE('right')
       iblk = 1
       DO i = 1,nvec
          offset = INT(pkd(i))
          IF(offset.GT.0) THEN
             m = vec_blk_size(ivec)
             n = vec_blk_size(i)
             k = vec_blk_size(i)
             CALL dgemm('N','N',m,n,k,alpha,pkd(offset),m,block(iblk),k,0.0_dp,buff(1),m)
             CALL dcopy(m*n,buff(1),1,pkd(offset),1)
          ENDIF
          iblk = iblk + vec_blk_size(i)**2
       ENDDO
    CASE DEFAULT
       CALL stop_program(routineN,"side?")
    END SELECT
    DEALLOCATE(buff)

  END SUBROUTINE packed_vec_scale

! *****************************************************************************
!> \brief Does the operations:
!>      v(k,:) = v(k,:) - uAk(i,k)' * inv(d(i,i)) * v(i,:) / s
!>      u(k,:) = u(k,:) -   v(i,k)' * inv(d(i,i)) * u(i,:) / s
!>
! *****************************************************************************
  SUBROUTINE packed_vec_bif_tech(mat_v, mat_u, pkd_v_fac, pkd_u_fac, ivec, pkd_v, pkd_u, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: mat_v, mat_u
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_v_fac, pkd_u_fac
    INTEGER, INTENT(IN)                      :: ivec
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_v, pkd_u
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_bif_tech', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, col, col_size, k_row, k_row_size, offset, row, row_size, &
      timing_handle, u_offset, u_offset_last, ufac_offset, v_offset, &
      v_offset_last, vfac_offset
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_u, data_v

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    row_blk_size => array_data (dbcsr_row_block_sizes (mat_u))
    col_blk_size => array_data (dbcsr_col_block_sizes (mat_u))

    IF(pkd_u(mat_u%m%nblkrows_total+1).LE.0) CALL stop_program(routineN,"1pkd_u(n+1).LE.0")
    IF(pkd_v(mat_v%m%nblkrows_total+1).LE.0) CALL stop_program(routineN,"1pkd_v(n+1).LE.0")
    CALL dbcsr_get_data (mat_v, data_v)
    CALL dbcsr_get_data (mat_u, data_u)
    !
    !v(k,:) (pkd) = v(k,:) (pkd) - (inv(d(i,i)) * uAk(i,k))' / s (pkd) * v(i,:) (bcsr) i = 1,ivec-1
    k_row = ivec
    k_row_size = row_blk_size(k_row)
    v_offset_last = INT(pkd_v(mat_v%m%nblkrows_total+1))
    DO row = 1, mat_v%m%nblkrows_total
       IF(row.GE.k_row) EXIT !i = 1,ivec-1
       row_size = row_blk_size(row)
       ufac_offset = INT(pkd_v_fac(row))
       IF(ufac_offset.LE.0) CYCLE
       DO blk = mat_v%m%row_p(row)+1,mat_v%m%row_p(row+1)
          col = mat_v%m%col_i(blk)
          col_size = col_blk_size(col)
          offset = ABS(mat_v%m%blk_p(blk))
          v_offset = INT(pkd_v(col))
          IF(v_offset.LE.0) THEN
             v_offset = v_offset_last
             v_offset_last = v_offset_last + k_row_size * col_size
             pkd_v(v_offset:v_offset+k_row_size * col_size-1) = 0.0_dp
             pkd_v(mat_v%m%nblkrows_total+1) = REAL(v_offset_last,dp)
          ENDIF
          !
          ! let's multiply and add
          CALL dgemm('T', 'N', k_row_size, col_size, row_size,&
               &    -1.0_dp, pkd_v_fac(ufac_offset), row_size,&
               &             data_v(offset), row_size,&
                  &     1.0_dp, pkd_v(v_offset), k_row_size)
       ENDDO
    ENDDO ! row
    !
    !u(k,:) (pkd) = u(k,:) (pkd) - (inv(d(i,i)) * v(i,k))' / s (pkd) * u(i,:) (bcsr) i = 1,ivec-1
    u_offset_last = INT(pkd_u(mat_u%m%nblkrows_total+1))
    DO row = 1, mat_u%m%nblkrows_total
       IF(row.GE.k_row) EXIT !i = 1,ivec-1
       row_size = row_blk_size(row)
       vfac_offset = INT(pkd_u_fac(row))
       IF(vfac_offset.LE.0) CYCLE
       DO blk = mat_u%m%row_p(row)+1,mat_u%m%row_p(row+1)
          col = mat_u%m%col_i(blk)
          col_size = col_blk_size(col)
          offset = ABS(mat_u%m%blk_p(blk))
          u_offset = INT(pkd_u(col))
          IF(u_offset.LE.0) THEN
             u_offset = u_offset_last
             u_offset_last = u_offset_last + k_row_size * col_size
             pkd_u(u_offset:u_offset+k_row_size * col_size-1) = 0.0_dp
             pkd_u(mat_u%m%nblkrows_total+1) = REAL(u_offset_last,dp)
          ENDIF
          !
          ! let's multiply and add
          CALL dgemm('T', 'N', k_row_size, col_size, row_size,&
               &    -1.0_dp, pkd_u_fac(vfac_offset), row_size,&
               &             data_u(offset), row_size,&
               &     1.0_dp, pkd_u(u_offset), k_row_size)
       ENDDO
    ENDDO ! row

    IF(pkd_u(mat_u%m%nblkrows_total+1).LE.0) CALL stop_program(routineN,"2pkd_u(n+1).LE.0")
    IF(pkd_v(mat_v%m%nblkrows_total+1).LE.0) CALL stop_program(routineN,"2pkd_v(n+1).LE.0")

    CALL timestop(timing_handle)

  END SUBROUTINE packed_vec_bif_tech


  SUBROUTINE packed_vec_build_u(pkd_u, pkd_v, k, n, s, vec_blk_size, error)
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: pkd_u
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_v
    INTEGER, INTENT(IN)                      :: k, n
    REAL(dp), INTENT(IN)                     :: s
    INTEGER, DIMENSION(:), INTENT(IN)        :: vec_blk_size
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_build_u', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i_blk_size, k_blk_size, &
                                                u_offset, v_offset

    pkd_u(1:n+1) = 0.0_dp
    u_offset = n+2
    !
    ! u(k,:) = [ v(k,1:k-1)/s -1 0 ... 0]
    k_blk_size = vec_blk_size(k)
    DO i = 1,k-1
       v_offset = INT(pkd_v(i))
       IF(v_offset.GT.0) THEN
          i_blk_size = vec_blk_size(i)
          CALL dcopy(k_blk_size*i_blk_size, pkd_v(v_offset),1,pkd_u(u_offset),1)
          CALL dscal(k_blk_size*i_blk_size, -1.0_dp/s, pkd_u(u_offset),1)
          pkd_u(i) = u_offset
          u_offset = u_offset + k_blk_size*i_blk_size
       ENDIF
    ENDDO
    !
    CALL block_set_d(k_blk_size, k_blk_size, pkd_u(u_offset:u_offset+ k_blk_size**2-1), &
         &         1.0_dp, 0.0_dp)
    pkd_u(k) = u_offset
    u_offset = u_offset + k_blk_size**2
    pkd_u(n+1) = u_offset

  END SUBROUTINE packed_vec_build_u

! *****************************************************************************
!> \brief Does the operations:
!>      v(k,:) = v(k,:) - uAk(i,k)' * inv(d(i,i)) * v(i,:) / s
!>
! *****************************************************************************
  SUBROUTINE packed_vec_bif_tech2(mat_v, pkd_v_fac, ivec, pkd_v, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: mat_v
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_v_fac
    INTEGER, INTENT(IN)                      :: ivec
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_v
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_bif_tech2', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, col, col_size, k_row, k_row_size, offset, row, row_size, &
      timing_handle, ufac_offset, v_offset, v_offset_last
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_v

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    row_blk_size => array_data (mat_v%m%row_blk_size)
    col_blk_size => array_data (mat_v%m%col_blk_size)

    IF(pkd_v(mat_v%m%nblkrows_total+1).LE.0) CALL stop_program(routineN,"1pkd_v(n+1).LE.0")
    CALL dbcsr_get_data (mat_v, data_v)

    !
    !v(k,:) (pkd) = v(k,:) (pkd) - (inv(d(i,i)) * uAk(i,k))' / s (pkd) * v(i,:) (bcsr) i = 1,ivec-1
    k_row = ivec
    k_row_size = row_blk_size(k_row)
    v_offset_last = INT(pkd_v(mat_v%m%nblkrows_total+1))
    DO row = 1, mat_v%m%nblkrows_total
       IF(row.GE.k_row) EXIT !i = 1,ivec-1
       row_size = row_blk_size(row)
       ufac_offset = INT(pkd_v_fac(row))
       IF(ufac_offset.LE.0) CYCLE
       DO blk = mat_v%m%row_p(row)+1,mat_v%m%row_p(row+1)
          col = mat_v%m%col_i(blk)
          col_size = col_blk_size(col)
          offset = ABS(mat_v%m%blk_p(blk))
          v_offset = INT(pkd_v(col))
          IF(v_offset.LE.0) THEN
             v_offset = v_offset_last
             v_offset_last = v_offset_last + k_row_size * col_size
             pkd_v(v_offset:v_offset+k_row_size * col_size-1) = 0.0_dp
             pkd_v(mat_v%m%nblkrows_total+1) = REAL(v_offset_last,dp)
          ENDIF
          !
          ! let's multiply and add
          CALL dgemm('T', 'N', k_row_size, col_size, row_size,&
               &    -1.0_dp, pkd_v_fac(ufac_offset), row_size,&
               &             data_v(offset), row_size,&
                  &     1.0_dp, pkd_v(v_offset), k_row_size)
       ENDDO
    ENDDO ! row
    !

    IF(pkd_v(mat_v%m%nblkrows_total+1).LE.0) CALL stop_program(routineN,"2pkd_v(n+1).LE.0")

    CALL timestop(timing_handle)

  END SUBROUTINE packed_vec_bif_tech2



  SUBROUTINE packed_vec_bcast(pkd_vec, source, scope, mp_obj, error)
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_vec
    INTEGER, INTENT(IN)                      :: source
    CHARACTER(LEN=*), INTENT(IN)             :: scope
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_bcast', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icol, irow, mp_group, mynode, &
                                                mypcol, myprow, npcols, &
                                                nprows, numnodes, src, &
                                                timing_handle
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    npcols = dbcsr_mp_npcols (mp_obj)
    nprows = dbcsr_mp_nprows (mp_obj)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)

    ALLOCATE(buff(SIZE(pkd_vec)))
    SELECT CASE(scope)
    CASE('rowise')
       !
       ! simple hack
       DO icol = 0,npcols-1
          buff = pkd_vec
          src = blacs2mpi(source,icol)
          CALL mp_bcast(buff,src,mp_group)
          IF(mypcol.EQ.icol) THEN
             pkd_vec = buff
          ENDIF
       ENDDO
    CASE('columnwise')
       !
       ! simple hack
       DO irow = 0,nprows-1
          buff = pkd_vec
          src = blacs2mpi(irow,source)
          CALL mp_bcast(buff,src,mp_group)
          IF(myprow.EQ.irow) THEN
             pkd_vec = buff
          ENDIF
       ENDDO
    CASE('all')
       CALL mp_bcast(pkd_vec,source,mp_group)
    CASE DEFAULT
       CALL stop_program(routineN,"how do you wanna bcast")
    END SELECT
    DEALLOCATE(buff)

    CALL timestop(timing_handle)

  END SUBROUTINE packed_vec_bcast

  SUBROUTINE packed_vec_reduce(pkd_vec, to, scope, vec_blk_size, n, mp_obj, error)
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_vec
    INTEGER, INTENT(IN)                      :: to
    CHARACTER(LEN=*), INTENT(IN)             :: scope
    INTEGER, DIMENSION(:), INTENT(IN)        :: vec_blk_size
    INTEGER, INTENT(IN)                      :: n
    TYPE(dbcsr_mp_obj)                       :: mp_obj
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'packed_vec_reduce', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icol, irow, mp_group, mynode, &
                                                mypcol, myprow, npcols, &
                                                nprows, numnodes, src, &
                                                timing_handle
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: buff, sums

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)

    numnodes = dbcsr_mp_numnodes (mp_obj)
    mynode = dbcsr_mp_mynode (mp_obj)
    myprow = dbcsr_mp_myprow (mp_obj)
    mypcol = dbcsr_mp_mypcol (mp_obj)
    npcols = dbcsr_mp_npcols (mp_obj)
    nprows = dbcsr_mp_nprows (mp_obj)
    blacs2mpi => dbcsr_mp_pgrid (mp_obj)
    mp_group = dbcsr_mp_group (mp_obj)

    IF(pkd_vec(n+1).LE.0) CALL stop_program(routineN,"a_offset_last.LE.0!!!!!!!")

    ALLOCATE(buff(SIZE(pkd_vec)),sums(SIZE(pkd_vec)))
    sums = 0.0_dp
    SELECT CASE(scope)
    CASE('rowise')
       CALL stop_program(routineN,"more work here")
       !
       ! simple hack
       DO icol = 0,npcols-1
          buff = pkd_vec
          src = blacs2mpi(to,icol)
          CALL mp_bcast(buff,src,mp_group)
          IF(mypcol.EQ.icol) THEN
             pkd_vec = buff
             CALL add_packed_vecs(sums,buff,vec_blk_size,n,error)
          ENDIF
       ENDDO
    CASE('columnwise')
       !CALL stop_program(routineN,"more work here")
       !
       ! simple hack
       DO irow = 0,nprows-1
          buff = pkd_vec
          src = blacs2mpi(irow,to)
          CALL mp_bcast(buff,src,mp_group)
          IF(myprow.EQ.irow) THEN
             CALL add_packed_vecs(sums,buff,vec_blk_size,n,error)
          ENDIF
       ENDDO
    CASE('all')
       CALL stop_program(routineN,"more work here")
    CASE DEFAULT
       CALL stop_program(routineN,"how do you wanna reduce")
    END SELECT
    pkd_vec = sums
    DEALLOCATE(buff,sums)

    IF(pkd_vec(n+1).LE.0) CALL stop_program(routineN,"a_offset_last.LE.0!!!!!!!")

    CALL timestop(timing_handle)

  END SUBROUTINE packed_vec_reduce

  SUBROUTINE add_packed_vecs(pkd_vec_a,pkd_vec_b,vec_blk_size,n,error)
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: pkd_vec_a
    REAL(dp), DIMENSION(:), INTENT(IN)       :: pkd_vec_b
    INTEGER, DIMENSION(:), INTENT(IN)        :: vec_blk_size
    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_packed_vecs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a_offset, a_offset_last, &
                                                b_offset, i

!   ---------------------------------------------------------------------------
! this points to the last empty entry

    a_offset_last = INT(pkd_vec_a(n+1)) 
    IF(a_offset_last.LE.0) CALL stop_program(routineN,"a_offset_last.LE.0!!!!!!!")
    DO i = 1,n
       b_offset = INT(pkd_vec_b(i))
       a_offset = INT(pkd_vec_a(i))
       IF(b_offset.LE.0) CYCLE
       IF(a_offset.GT.0) THEN
          ! the block exsits in a, just add
          CALL daxpy(vec_blk_size(i),1.0_dp,pkd_vec_b(b_offset),1,pkd_vec_a(a_offset),1)
       ELSE
          ! the block doesnt exsit in a, copy b at the end
          CALL dcopy(vec_blk_size(i),pkd_vec_b(b_offset),1,pkd_vec_a(a_offset_last),1)
          a_offset_last = a_offset_last + vec_blk_size(i)
       ENDIF
    ENDDO
    ! reset the last empty entry if needed
    pkd_vec_a(n+1) = REAL(a_offset_last,dp)
  END SUBROUTINE add_packed_vecs

! *****************************************************************************
!> \brief copy a vector (row or column) from a matrix to another matrix
!> \param[inout] matrix_b  matrix
!> \param[in]    matrix_a  matrix
!> \param[in]    what      row or column
!> \param[in]    ivec      the column/row index
!> 
! *****************************************************************************
  SUBROUTINE dbcsr_copy_vec(matrix_a, matrix_b, what, ivec, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix_a
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix_b
    CHARACTER(len=*), INTENT(in)             :: what
    INTEGER, INTENT(IN)                      :: ivec
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_copy_vec', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: a_data_z, b_data_z
    COMPLEX(KIND=sp), DIMENSION(:), POINTER  :: a_data_c, b_data_c
    INTEGER :: a_blk, a_col_size, a_frst_blk, a_last_blk, a_offset, &
      a_row_size, b_blk, b_col, b_col_size, b_offset, b_row_size, col, nze, &
      row, timing_handle
    INTEGER, DIMENSION(:), POINTER           :: a_col_blk_size, &
                                                a_row_blk_size, &
                                                b_col_blk_size, b_row_blk_size
    LOGICAL                                  :: found
    REAL(KIND=dp), DIMENSION(:), POINTER     :: a_data_d, b_data_d
    REAL(KIND=sp), DIMENSION(:), POINTER     :: a_data_r, b_data_r

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, timing_handle)
    
    a_row_blk_size => array_data (matrix_a%m%row_blk_size)
    a_col_blk_size => array_data (matrix_a%m%col_blk_size)
    b_row_blk_size => array_data (matrix_b%m%row_blk_size)
    b_col_blk_size => array_data (matrix_b%m%col_blk_size)
    SELECT CASE (matrix_a%m%data_type)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (matrix_a%m%data_area, a_data_r)
       CALL dbcsr_get_data (matrix_b%m%data_area, b_data_r)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (matrix_a%m%data_area, a_data_d)
       CALL dbcsr_get_data (matrix_b%m%data_area, b_data_d)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (matrix_a%m%data_area, a_data_c)
       CALL dbcsr_get_data (matrix_b%m%data_area, b_data_c)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (matrix_a%m%data_area, a_data_z)
       CALL dbcsr_get_data (matrix_b%m%data_area, b_data_z)
    END SELECT
    !
    ! let's go
    SELECT CASE(what)
    CASE('column')
       col = ivec
       CALL stop_program(routineN, "NYI")
    CASE('row')
       row = ivec
       a_row_size = a_row_blk_size(row)
       b_row_size = b_row_blk_size(row)
       IF(a_row_size.NE.b_row_size) CALL stop_program(routineN, "matrices not consistent")
       a_frst_blk = matrix_a%m%row_p(row)+1
       a_last_blk = matrix_a%m%row_p(row+1)
       DO b_blk = matrix_b%m%row_p(row)+1,matrix_b%m%row_p(row+1)
          b_col = matrix_b%m%col_i(b_blk)
          CALL dbcsr_find_column(b_col, a_frst_blk, a_last_blk,&
               matrix_a%m%col_i, matrix_a%m%blk_p, a_blk, found)
          IF(found) THEN
             a_col_size = a_col_blk_size(b_col)
             b_col_size = b_col_blk_size(b_col)
             IF(a_col_size.NE.b_col_size) CALL stop_program(routineN, "matrices not consistent")
             nze = a_row_size * a_col_size
             !
             ! let's copy the block
             a_offset = ABS(matrix_a%m%blk_p(a_blk))
             b_offset = ABS(matrix_b%m%blk_p(b_blk))
             SELECT CASE (matrix_a%m%data_type)
             CASE (dbcsr_type_real_8)
                CALL dcopy (nze, b_data_d(b_offset), 1, a_data_d(a_offset), 1)
             !CALL dcopy(nze, matrix_b%m%data(b_offset), 1, matrix_a%m%data(a_offset), 1)
             CASE default
                CALL cp_unimplemented_error (routineN,&
                     "Only real double precision")
             END SELECT
          ELSE
             CALL stop_program(routineN, "need ot create the block, NYI")
          ENDIF
       ENDDO ! a_col
    CASE DEFAULT
       CALL stop_program(routineN,"copy what?")
    END SELECT
    !
    CALL timestop(timing_handle)

  END SUBROUTINE dbcsr_copy_vec


! *****************************************************************************
!> \brief Emulation of sparse_matrix_types/add_block_node mapped
!>        to add_real_matrix_block.
!>
!> It adds a block to the dbcsr matrix and returns a rank-2 pointer to the
!> block. Currently it only and always uses the mutable data.
!> \param[in,out] matrix      DBCSR matrix
!> \param[in]  row            the row
!> \param[in]  col            the column
!> \param[in]  block          the block to put
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_add_block_node (matrix, block_row, block_col, block, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_add_block_node', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: existed, is_there

!   ---------------------------------------------------------------------------
!    CALL cp_assert (matrix%m%work_mutable, cp_warning_level, cp_caller_error,&
!         routineN, "Mutable not defined upon DBCSR creation, forcing its use.")

    IF (.NOT. ASSOCIATED (matrix%m%wms)) THEN
       CALL dbcsr_work_create (matrix)
       matrix%m%valid = .FALSE.
       CALL dbcsr_mutable_new (matrix%m%wms(1)%mutable,&
            dbcsr_get_data_type (matrix))
    ENDIF
    CALL cp_assert (.NOT. ASSOCIATED (matrix%m%wms(1)%row_i), cp_warning_level,&
         cp_unimplemented_error_nr, routineN,&
         "Data loss due to no conversion of appendable to mutable data")
    is_there = ASSOCIATED(block)
    !r = row ; c = col ; tr = .FALSE.
    !CALL dbcsr_get_stored_coordinates (matrix, r, c, tr)
    !CALL dbcsr_reserve_block2d (matrix, row, col, block)
    !write(*,*) 'add_block_node: block_row',block_row,' block_col',block_col
    CALL dbcsr_reserve_block2d (matrix, block_row, block_col, block,&
         existed=existed)
    CALL cp_assert (.NOT.existed, cp_warning_level, cp_wrong_args_error,&
         routineN, "You should not add existing blocks according to old API.")
    IF(.NOT.is_there) block(:,:) = 0.0_dp
  END SUBROUTINE dbcsr_add_block_node


! *****************************************************************************
!> \brief Emulation of sparse_matrix_types/get_block_node mapped to
!>        get_real_block_node
! *****************************************************************************
  SUBROUTINE get_dbcsr_block_node(block_node,block_col,BLOCK)

    TYPE(dbcsr_block_node_type), INTENT(IN)  :: block_node
    INTEGER, INTENT(OUT), OPTIONAL           :: block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK

!   ---------------------------------------------------------------------------

    IF (PRESENT(block_col)) block_col = block_node%block_col
    IF (PRESENT(BLOCK)) BLOCK => block_node%block
  END SUBROUTINE get_dbcsr_block_node


! *****************************************************************************
!> \brief Emulation of sparse_matrix_types/get_block_node mapped to
!>        get_matrix_block
! *****************************************************************************
  SUBROUTINE get_dbcsr_matrix_block(matrix,block_row,block_col,&
                                   first_row,last_row,first_col,last_col,&
                                   block_node,BLOCK)
    TYPE(dbcsr_obj)                          :: matrix
    INTEGER, INTENT(IN)                      :: block_row, block_col
    INTEGER, INTENT(OUT), OPTIONAL           :: first_row, last_row, &
                                                first_col, last_col
    TYPE(dbcsr_block_node_type), OPTIONAL, &
      POINTER                                :: block_node
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: BLOCK

    INTEGER                                  :: block_number, cbs, first_pos, &
                                                offset, rbs, stored_col, &
                                                stored_row
    INTEGER, DIMENSION(:), POINTER           :: bs
    LOGICAL                                  :: found, stored_tr

!   ---------------------------------------------------------------------------

    IF (PRESENT(first_row) .OR. PRESENT(last_row)) THEN
       bs = array_data (dbcsr_row_block_sizes (matrix))
       first_pos = 1 + SUM (bs(1:block_row-1))
       IF (PRESENT(first_row)) first_row = first_pos
       IF (PRESENT(last_row)) last_row = first_pos - 1 + bs(block_row)
    ENDIF

    IF (PRESENT(first_col) .OR. PRESENT(last_col)) THEN
       bs = array_data (dbcsr_col_block_sizes (matrix))
       first_pos = 1 + SUM (bs(1:block_col-1))
       IF (PRESENT(first_col)) first_col = first_pos
       IF (PRESENT(last_col)) last_col = first_pos - 1 + bs(block_col)
    ENDIF
    CALL make_block_node_exist (matrix)
    CALL dbcsr_get_block_index (matrix, block_row, block_col, stored_row,&
         stored_col, stored_tr, found, block_number, offset)
    IF (found) THEN
       matrix%m%block_node%block_col = block_col
       rbs = matrix%m%block_node%iter%rbs(block_row)
       cbs = matrix%m%block_node%iter%cbs(block_col)
       IF (PRESENT (block_node)) block_node => matrix%m%block_node
       IF (PRESENT (block)) THEN
          CALL cp_assert (.NOT. stored_tr, cp_warning_level, cp_internal_error,&
               "get_matrix_block", "Block is stored transposed")
          CALL pointer_set (block,&
               (/ rbs, cbs /), pointer_view (&
               dbcsr_get_data_p (matrix%m%data_area, 0.0_dp), offset,&
               offset+rbs*cbs-1),&
               matrix%m%block_node%iter%buffer_r_dp)
       ENDIF
    ELSE
       IF (PRESENT (block_node)) NULLIFY (block_node)
       IF (PRESENT (block)) NULLIFY (block)
    ENDIF
  END SUBROUTINE get_dbcsr_matrix_block

  SUBROUTINE make_block_node_exist(matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    IF (.NOT. ASSOCIATED (matrix%m%block_node)) THEN
       ALLOCATE (matrix%m%block_node)
       CALL dbcsr_iterator_start (matrix%m%block_node%iter, matrix)
    ENDIF
  END SUBROUTINE make_block_node_exist

!  INTERFACE dbcsr_put_block_node
!     MODULE PROCEDURE put_block_node, put_matrix_block
!  END INTERFACE
!
!  INTERFACE dbcsr_next_block_node
!     MODULE PROCEDURE next_block_node
!  END INTERFACE


! *****************************************************************************
!> \brief Emulation of sparse_matrix_types/allocate_matrix mapped to
!>        allocate_real_matrix
!> \param[in,out] matrix      DBCSR matrix
!> \param[in] nrow            number of full rows
!> \param[in] ncol            number of full columns
!> \param[in] nblock_row      number of blocked rows
!> \param[in] nblock_col      number of blocked columns
!> \param[in] first_row       starting offsets of blocked rows
!> \param[in] last_row        ending offsets of blocked rows
!> \param[in] first_col       starting offsets of blocked columns
!> \param[in] last_col        ending offsets of blocked columns
!> \param[in] matrix_name     name of matrix
!> \param[in] matrix_symmetry matrix symmetry type (mapped do DBCSR symmetries)
!> \param[in] sparsity_id     sparsity id (not used)
!> \param[in] distribution_2d use this distribution
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_allocate_matrix (matrix,nrow,ncol,nblock_row,nblock_col,&
       first_row,last_row,first_col,last_col,matrix_name,matrix_symmetry,&
       sparsity_id,distribution_2d,replication_type,error)
    TYPE(dbcsr_obj), POINTER                 :: matrix
    INTEGER, INTENT(IN)                      :: nrow, ncol, nblock_row, &
                                                nblock_col
    INTEGER, DIMENSION(nblock_row), &
      INTENT(IN)                             :: first_row, last_row
    INTEGER, DIMENSION(nblock_col), &
      INTENT(IN)                             :: first_col, last_col
    CHARACTER(LEN=*), INTENT(IN)             :: matrix_name, matrix_symmetry
    INTEGER, INTENT(IN)                      :: sparsity_id
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    CHARACTER, INTENT(IN), OPTIONAL          :: replication_type
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_allocate_matrix', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: dbcsr_matrix_type
    INTEGER                                  :: stat
    INTEGER, DIMENSION(:), POINTER           :: cbs, rbs
    TYPE(array_i1d_obj)                      :: col_blk_sizes, row_blk_sizes
    TYPE(dbcsr_distribution_obj)             :: dbcsr_dist

!   ---------------------------------------------------------------------------
! Setup the distribution

    IF (ASSOCIATED (matrix)) CALL dbcsr_deallocate_matrix (matrix, error)
    ALLOCATE (matrix)
    CALL dbcsr_init (matrix, error)
    CALL cp_assert (ASSOCIATED(distribution_2d),&
         cp_fatal_level, cp_caller_error, routineN,&
         "distribution_2d not associated", error=error)
    CALL dbcsr_dist2d_to_dist (distribution_2d, dbcsr_dist, error)
    CALL cp_assert (&
         dbcsr_distribution_nrows (dbcsr_dist) .EQ. nblock_row,&
         cp_warning_level, cp_caller_error, routineN,&
         "Inconsistent number of blocked rows", error=error)
    CALL cp_assert (&
         dbcsr_distribution_ncols (dbcsr_dist) .EQ. nblock_row,&
         cp_warning_level, cp_caller_error, routineN,&
         "Inconsistent number of blocked columns", error=error)
    ! Create row and column block sizes arrays.
    ALLOCATE (rbs(nblock_row), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineP, "rbs", nblock_row)
    ALLOCATE (cbs(nblock_col), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineP, "cbs", nblock_col)
    CALL convert_offsets_to_sizes (first_row, rbs, last_row)
    CALL convert_offsets_to_sizes (first_col, cbs, last_col)
    CALL array_nullify (row_blk_sizes)
    CALL array_nullify (col_blk_sizes)
    CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
    CALL array_new (col_blk_sizes, cbs, gift=.TRUE.)
    ! Transform matrix type
    CALL dbcsr_type_from_sm_symmetry (dbcsr_matrix_type, matrix_symmetry)
    !SELECT CASE (matrix_symmetry)
    !CASE ("symmetric")
    !   dbcsr_matrix_type = dbcsr_type_symmetric
    !CASE ("anisymmetric")
    !   dbcsr_matrix_type = dbcsr_type_antisymmetric
    !CASE ("none", "no symmetry")
    !   dbcsr_matrix_type = dbcsr_type_normal
    !CASE default
    !   CALL cp_assert (.FALSE., cp_warning_level, cp_caller_error, routineN,&
    !        "Unknown matrix type "//"matrix_symmetry", failure=failure,&
    !        error=error)
    !END SELECT
    ! Create the matrix
    CALL dbcsr_create (matrix, matrix_name, dbcsr_dist, dbcsr_matrix_type,&
         row_blk_sizes, col_blk_sizes, 0, 0, dbcsr_type_real_default,&
         mutable_work=.TRUE., replication_type=replication_type, error=error)
    CALL cp_dbcsr_add_dist2d (matrix%m, distribution_2d, error)
    matrix%m%sparsity_id = sparsity_id
    ! Clean up
    CALL dbcsr_distribution_release (dbcsr_dist)
    CALL array_release (row_blk_sizes)
    CALL array_release (col_blk_sizes)
    CALL dbcsr_work_create (matrix)
    matrix%m%valid = .FALSE.
    CALL dbcsr_mutable_new (matrix%m%wms(1)%mutable,&
         dbcsr_get_data_type (matrix))
  END SUBROUTINE dbcsr_allocate_matrix


! *****************************************************************************
!> \brief Deallocates a DBCSR matrix for compatibility with CP2K
!> \param[in,out] matrix      DBCSR matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_deallocate_matrix(matrix, error)
    TYPE(dbcsr_obj), POINTER                 :: matrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_deallocate_matrix', &
      routineP = moduleN//':'//routineN

    CALL dbcsr_release (matrix)
    CALL cp_assert (.NOT. dbcsr_valid_index(matrix), cp_warning_level,&
         cp_caller_error, routineN,&
         'You should not "deallocate" a referenced matrix. '//&
         'Avoid pointers to DBCSR matrices.')
    DEALLOCATE (matrix)
    NULLIFY (matrix)
  END SUBROUTINE dbcsr_deallocate_matrix


! *****************************************************************************
!> \brief Emulation of sparse_matrix_types/replicate_matrix_structure mapped to
!>        replicate_real_matrix_structure
!> \note This just sets up the metadata; it does not acutally replicate the
!>       block structure.
!> \param[in] matrix          source matrix to replicate
!> \param[out] nrow           replicated matrix
!> \param[in] target_name     name of the replicated matrix
!> \param[in] target_symmetry symmetry of the target matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_replicate_matrix_struct(source,TARGET,target_name,&
                                             target_symmetry,error)
    TYPE(dbcsr_obj), INTENT(IN)              :: source
    TYPE(dbcsr_obj), POINTER                 :: TARGET
    CHARACTER(LEN=*), INTENT(IN)             :: target_name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: target_symmetry
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'dbcsr_replicate_matrix_struct', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: dbcsr_matrix_type

!   ---------------------------------------------------------------------------

    IF (PRESENT (target_symmetry)) THEN
       CALL dbcsr_type_from_sm_symmetry (dbcsr_matrix_type, target_symmetry)
    ELSE
       dbcsr_matrix_type = dbcsr_get_matrix_type (source)
    ENDIF
    ALLOCATE (TARGET)
    CALL dbcsr_init (TARGET)
    CALL dbcsr_create (TARGET, target_name, dbcsr_distribution (source),&
         dbcsr_matrix_type,&
         dbcsr_row_block_sizes (source), dbcsr_col_block_sizes (source),&
         0, 0, dbcsr_get_data_type (source),&
         dbcsr_uses_special_memory (source), mutable_work=source%m%work_mutable,&
         error=error)
    CALL cp_dbcsr_add_dist2d (TARGET%m, source%m%dist2d, error)
    TARGET%m%sparsity_id=source%m%sparsity_id
  END SUBROUTINE dbcsr_replicate_matrix_struct


! *****************************************************************************
! CP2k-compatible matrix sets
! *****************************************************************************


! *****************************************************************************
!> \brief   Allocate and initialize a real matrix 1-dimensional set.
!> \param[in,out] matrix_set  Set containing the DBCSR matrices
!> \param[in] nmatrix         Size of set
!> \param[in,out] error       cp2k error
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! *****************************************************************************
  SUBROUTINE allocate_dbcsr_matrix_set(matrix_set, nmatrix, error)
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: matrix_set
    INTEGER, INTENT(IN)                      :: nmatrix
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_dbcsr_matrix_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix, istat

    IF (ASSOCIATED(matrix_set)) CALL dbcsr_deallocate_matrix_set(matrix_set,error=error)
    ALLOCATE (matrix_set(nmatrix),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"matrix_set",0)
    DO imatrix=1,nmatrix
       NULLIFY (matrix_set(imatrix)%matrix)
       !ALLOCATE (matrix_set(imatrix)%matrix)
       !CALL dbcsr_init (matrix_set(imatrix)%matrix)
    END DO
  END SUBROUTINE allocate_dbcsr_matrix_set

! *****************************************************************************
!> \brief   Allocate and initialize a real matrix 2-dimensional set.
!> \param[in,out] matrix_set  Set containing the DBCSR matrix pointer type
!> \param[in] nmatrix         Size of set
!> \param[in,out] error       cp2k error
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! *****************************************************************************
  SUBROUTINE allocate_dbcsr_matrix_set_2d(matrix_set,nmatrix,mmatrix,error)
    TYPE(dbcsr_obj_p_type), &
      DIMENSION(:, :), POINTER               :: matrix_set
    INTEGER, INTENT(IN)                      :: nmatrix, mmatrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_dbcsr_matrix_set_2d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix, istat, jmatrix

    IF (ASSOCIATED(matrix_set)) CALL dbcsr_deallocate_matrix_set(matrix_set,error=error)
    ALLOCATE (matrix_set(nmatrix,mmatrix),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"matrix_set",0)
    DO jmatrix=1,mmatrix
      DO imatrix=1,nmatrix
         NULLIFY (matrix_set(imatrix,jmatrix)%matrix)
         !ALLOCATE (matrix_set(imatrix,jmatrix)%matrix)
         !CALL dbcsr_init (matrix_set(imatrix,jmatrix)%matrix)
      END DO
    END DO
  END SUBROUTINE allocate_dbcsr_matrix_set_2d


! *****************************************************************************
!> \brief Deallocate a real matrix set and release all of the member matrices.
!> \param[in,out] matrix_set  Set containing the DBCSR matrix pointer type
!> \param[in,out] error       cp2k error
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! *****************************************************************************
  SUBROUTINE deallocate_dbcsr_matrix_set(matrix_set,error)

    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: matrix_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_dbcsr_matrix_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix, istat

    IF (ASSOCIATED(matrix_set)) THEN
      DO imatrix=1,SIZE(matrix_set)
        CALL dbcsr_deallocate_matrix(matrix_set(imatrix)%matrix,error=error)
      END DO
      DEALLOCATE (matrix_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"matrix_set")
    END IF
  END SUBROUTINE deallocate_dbcsr_matrix_set

! *****************************************************************************
!> \brief Deallocate a real matrix set and release all of the member matrices.
!> \param[in,out] matrix_set  Set containing the DBCSR matrix pointer type
!> \param[in,out] error       cp2k error
!> \par History
!>      2009-08-17 Adapted from sparse_matrix_type for DBCSR
! *****************************************************************************
  SUBROUTINE deallocate_dbcsr_matrix_set_2d(matrix_set,error)

    TYPE(dbcsr_obj_p_type), &
      DIMENSION(:, :), POINTER               :: matrix_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'deallocate_dbcsr_matrix_set_2d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imatrix, istat, jmatrix

    IF (ASSOCIATED(matrix_set)) THEN
      DO jmatrix=1,SIZE(matrix_set,2)
        DO imatrix=1,SIZE(matrix_set,1)
          CALL dbcsr_deallocate_matrix(matrix_set(imatrix,jmatrix)%matrix,error=error)
        END DO
      END DO
      DEALLOCATE (matrix_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"matrix_set")
    END IF
  END SUBROUTINE deallocate_dbcsr_matrix_set_2d



! *****************************************************************************
!> \brief Maps sparse_matrix_type symmetry to DBCSR matrix type
!> \param[out] dbcsr_matrix_type        matrix_type of the DBCSR matrix
!> \param[in] sm_symmetry               sparse_matrix_symmetry
! *****************************************************************************
  SUBROUTINE dbcsr_type_from_sm_symmetry (dbcsr_matrix_type, sm_symmetry)
    CHARACTER, INTENT(OUT)                   :: dbcsr_matrix_type
    CHARACTER(LEN=*), INTENT(IN)             :: sm_symmetry

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_type_from_sm_symmetry', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    SELECT CASE (sm_symmetry)
    CASE ("symmetric")
       dbcsr_matrix_type = dbcsr_type_symmetric
    CASE ("antisymmetric")
       dbcsr_matrix_type = dbcsr_type_antisymmetric
    CASE ("none", "no symmetry")
       dbcsr_matrix_type = dbcsr_type_normal
    CASE default
       CALL cp_assert (.FALSE., cp_warning_level, cp_caller_error, routineN,&
            "Unknown matrix type "//sm_symmetry)
    END SELECT
  END SUBROUTINE dbcsr_type_from_sm_symmetry



! *****************************************************************************
!> \brief Converts a cp2k full matrix into a DBCSR matrix.
!> \param[out] matrix         the created BCSR matrix
!> \param[in] fm              the cpk full matrix
!> \param[in] threshold       the threshold for determining sparsity
!> \param[in] row_blk_sizes   sizes of row blocks
!> \param[in] col_blk_sizes   sizes of column blocks
!> \param[in] distribution    the distribution to use for the new matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_from_fm(matrix, fm, threshold, distribution, row_blk_size,&
       col_blk_size, error)
    TYPE(dbcsr_obj), INTENT(OUT)             :: matrix
    TYPE(cp_fm_type), POINTER                :: fm
    REAL(KIND=dp), INTENT(IN)                :: threshold
    TYPE(dbcsr_distribution_obj)             :: distribution
    TYPE(array_i1d_obj), INTENT(IN)          :: row_blk_size, col_blk_size
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_from_fm', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=80)                        :: fm_name
    INTEGER :: blkcol, blkcol_l, blkrow, blkrow_l, fullrow, i, j, m, &
      mp_group, mynode, mypcol, myprow, n, nblkcols_local, nblkcols_total, &
      nblkrows_local, nblkrows_total, nfullcols_local, nfullcols_total, &
      nfullrows_local, nfullrows_total, npcols, nprows, numnodes, &
      prev_blkrow, prow, stat, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: blk_ps, col_blk_begin, &
                                                row_blk_begin
    INTEGER, DIMENSION(:), POINTER :: blkcols_local, blkrows_local, cbs, &
      fullcols_local, fullrows_local, pcol_dist, prow_dist, rbs
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    LOGICAL                                  :: dbg, failure = .FALSE.
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: thefmrow
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(dbcsr_mp_obj)                       :: mp_env

!   ---------------------------------------------------------------------------

    dbg = .FALSE.
    CALL dbcsr_set_debug(dbg)
    CALL timeset(routineN, timing_handle)
    CALL cp_fm_get_info(fm, name=fm_name, matrix_struct=fm_struct,&
         row_indices=fullrows_local, col_indices=fullcols_local,&
         nrow_global=nfullrows_total, ncol_global=nfullcols_total,&
         error=error)
    prow_dist => array_data (dbcsr_distribution_row_dist (distribution))
    pcol_dist => array_data (dbcsr_distribution_col_dist (distribution))
    nblkrows_local = dbcsr_distribution_nlocal_rows (distribution)
    nblkcols_local = dbcsr_distribution_nlocal_cols (distribution)
    nblkrows_total = dbcsr_distribution_nrows (distribution)
    nblkcols_total = dbcsr_distribution_ncols (distribution)
    mp_env = dbcsr_distribution_mp (distribution)
    blkrows_local => array_data (dbcsr_distribution_local_rows (distribution))
    blkcols_local => array_data (dbcsr_distribution_local_cols (distribution))
    pgrid => dbcsr_mp_pgrid (mp_env)
    mynode = dbcsr_mp_mynode (mp_env)
    numnodes = dbcsr_mp_numnodes (mp_env)
    myprow = dbcsr_mp_myprow (mp_env)
    mypcol = dbcsr_mp_mypcol (mp_env)
    nprows = dbcsr_mp_nprows (mp_env)
    npcols = dbcsr_mp_npcols (mp_env)
    mp_group = dbcsr_mp_group (mp_env)
    rbs => array_data (row_blk_size)
    cbs => array_data (col_blk_size)
    nfullrows_local = 0
    nfullcols_local = 0
    IF (dbg) WRITE(*,*)routineP//' myprow, mypcol',myprow,'/',nprows,';',&
         mypcol,'/',npcols
    IF (dbg) WRITE(*,*)routineP//' row_blk_size=',rbs
    IF (dbg) WRITE(*,*)routineP//' blkrows_local=',blkrows_local
    IF (dbg) WRITE(*,*)routineP//' blkcols_local=',blkcols_local
    ALLOCATE(row_blk_begin(nblkrows_total))
    row_blk_begin(1) = 1
    DO blkrow = 2, nblkrows_total
       row_blk_begin(blkrow) = row_blk_begin(blkrow-1)+rbs(blkrow-1)
    ENDDO
    ALLOCATE(col_blk_begin(nblkcols_total))
    col_blk_begin(1) = 1
    DO blkcol = 2, nblkcols_total
       col_blk_begin(blkcol) = col_blk_begin(blkcol-1)+cbs(blkcol-1)
    ENDDO
    DO blkrow_l = 1, nblkrows_local
       blkrow = blkrows_local(blkrow_l)
       nfullrows_local = nfullrows_local + rbs(blkrow)
    ENDDO
    DO blkcol_l = 1, nblkcols_local
       blkcol = blkcols_local(blkcol_l)
       nfullcols_local = nfullcols_local + cbs(blkcol)
    ENDDO
    ! Create local work matrix.
    !IF (ASSOCIATED (matrix)) CALL dbcsr_deallocate (matrix, error)
    !ALLOCATE (matrix)
    CALL dbcsr_init (matrix, error)
    CALL dbcsr_create(matrix, fm%name, distribution, 'N',&
         row_blk_size, col_blk_size,&
         0, 0, dbcsr_type_real_8, error=error)
    ! Somehow assign row & column block sizes.
    CALL dbcsr_work_create(matrix, error=error)
    ALLOCATE(blk_ps(nblkcols_local), stat=stat)
    CALL ensure_array_size (matrix%m%wms(1)%row_i, ub=1)
    CALL ensure_array_size (matrix%m%wms(1)%col_i, ub=1)
    CALL ensure_array_size (matrix%m%wms(1)%blk_p, ub=1)
    IF (stat /= 0) CALL stop_memory(routineP, 'blk_ps', nblkcols_local+1)
    ! Figure out the maximum row & column block sizes
    m = MAXVAL(rbs)
    n = MAXVAL(cbs)
    IF(dbg) WRITE(*,*)routineP//' max size',m,'x',n
    ALLOCATE(thefmrow(m,nfullcols_total), stat=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    !> This is what happens here: I go through the rows and obtain the FM block.
    !> How do I update the indexes of the wm? I just add to the array in
    !> the wm matrix and add some data to the data array.
    prev_blkrow = 0
    fullrow = 1
    matrix%m%wms(1)%lastblk = 0
    matrix%m%wms(1)%datasize = 0
    !matrix%m%wms(1)%data(:) = 999
    IF(dbg) WRITE(*,*)routineP//' my local blocks',blkrows_local
    DO blkrow = 1, nblkrows_total
       prow = prow_dist(blkrow)
       ! This processor must participate in forming the thefmrow, but I don't
       ! need it until it's one of my local rows.
       IF(dbg) WRITE(*,*)routineP//' blocked row', blkrow
       m = rbs(blkrow)
       IF(dbg) WRITE(*,*)routineP//' Submatrix',rbs(blkrow),' x 1:',&
            nfullcols_total,'size',m,'x',&
            nfullcols_total
       DEALLOCATE(thefmrow)
       ALLOCATE(thefmrow(m,nfullcols_total), stat=stat)
       IF (stat /= 0) CALL stop_memory(routineP, 'thefmrow', m*n)
       thefmrow(:,:) = 888
       CALL cp_fm_get_submatrix(fm, thefmrow,&
            start_row=row_blk_begin(blkrow), n_rows=m,&
            start_col=1, n_cols=nfullcols_total, error=error)
       !IF (dbg) WRITE(*,*)thefmrow
       IF(dbg) WRITE(*,*)routineP//' Chksum=',&
            DOT_PRODUCT(RESHAPE(thefmrow,(/m*nfullcols_total/)),&
            &           RESHAPE(thefmrow,(/m*nfullcols_total/)))
       IF (myprow .EQ. prow) THEN
          ! The array-adding should be done in one pass
          IF(dbg) WRITE(*,*)routineP//' adding row',blkrow,'sized',m,&
               'for local columns', blkcols_local,'fullrows',fullrow
          CALL addto_array(matrix%m%wms(1)%row_i,&
               (/ (blkrow, i=1,nblkcols_local) /),&
               matrix%m%wms(1)%lastblk+1, error=error)
          CALL addto_array(matrix%m%wms(1)%col_i, blkcols_local,&
               matrix%m%wms(1)%lastblk+1, error=error)
          DO blkcol_l = 1, nblkcols_local
             blkcol = blkcols_local(blkcol_l)
             i = col_blk_begin(blkcol)
             j = i + cbs(blkcol)-1
             n = cbs(blkcol)
             blk_ps(blkcol_l) = matrix%m%wms(1)%datasize + 1
             IF(dbg) WRITE(*,*)routineP//' adding block at column',blkcol,&
                  'sized',n,'(',i,':',j,')'
             IF(dbg) WRITE(*,*)routineP//' with chksum',&
                  DOT_PRODUCT(RESHAPE(thefmrow(1:m,i:j),(/m*n/)),&
                  &           RESHAPE(thefmrow(1:m,i:j),(/m*n/)))
             CALL addto_array(matrix%m%wms(1)%data_area,&
                  RESHAPE(thefmrow(1:m,i:j),(/m*n/)),&
                  offset=blk_ps(blkcol_l))
             !CALL addto_array(matrix%m%wms(1)%data,&
             !     RESHAPE(thefmrow(1:m,i:j),(/m*n/)),&
             !     offset=blk_ps(blkcol_l), error=error)
             matrix%m%wms(1)%datasize = matrix%m%wms(1)%datasize + m*n
          ENDDO
          IF(dbg) WRITE(*,*)routineP//' with blk_ps',blk_ps
          CALL addto_array(matrix%m%wms(1)%blk_p, blk_ps(1:nblkcols_local),&
               matrix%m%wms(1)%lastblk+1, error=error)
          matrix%m%wms(1)%lastblk = matrix%m%wms(1)%lastblk + nblkcols_local
       ENDIF ! myprow = prow
    ENDDO !blkrow
    DEALLOCATE(blk_ps, thefmrow)
    ! Construct the final matrix.
    CALL dbcsr_finalize(matrix, error=error)
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_from_fm


! *****************************************************************************
!> \brief   Copy a BLACS matrix to a dbcsr matrix.
!>
!>          real_matrix=beta*real_matrix+alpha*fm
!>          beta defaults to 0, alpha to 1
!> \author  Matthias Krack
!> \date    06.06.2001
!> \par History
!>          08.2002 adapted to local_data, could be optimized
!>          08.2003 imported form qs_blacs
!>          2009-09-01 adapted for DBCSR
!> \version 1.0
!> \param[in] fm              full matrix
!> \param[out] matrix         DBCSR matrix
!> \param[in] alpha           (optional) scaling of FM
!> \param[in] beta            (optional) scalinig of existing SM
! *****************************************************************************
  SUBROUTINE copy_fm_to_dbcsr(fm,matrix,alpha,beta)
    TYPE(cp_fm_type), POINTER           :: fm
    TYPE(dbcsr_obj), INTENT(INOUT)      :: matrix
    REAL(kind=dp), INTENT(in), OPTIONAL :: alpha,beta

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_fm_to_dbcsr', &
      routineP = moduleN//':'//routineN

    TYPE(real_block_node_type), POINTER :: block_node

    INTEGER :: group,handle,iblock_col,iblock_row,icol,icol_global,&
         icol_local,ipcol,ipe,iprow,irow,irow_global,irow_local,istat,&
         jpcol,jprow,mypcol,mype,myprow,nblock_row,nblkcols_total,&
         ncol_local,npcol,npe,nprow,nblkrows_total,nrow_local, nblock_col,&
         nblkrows_local, nblkcols_local, row_l, col_l, nrow_block, ncol_block

    INTEGER, DIMENSION(:), POINTER    :: row_blk_sizes, col_blk_sizes, local_rows, local_cols
    INTEGER, DIMENSION(:), ALLOCATABLE :: first_col,first_row,last_col,last_row
         
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: fm_block, real_matrix_block
    TYPE(cp_blacs_env_type), POINTER  :: context
    REAL(kind=dp) :: my_alpha,my_beta
    LOGICAL :: tr, found
    TYPE(array_i1d_obj) :: rbs, cbs

#if defined(__SCALAPACK)
    INTEGER, EXTERNAL :: indxg2l,indxg2p

#endif

    CALL timeset(routineN,handle)

    my_beta=0._dp
    my_alpha=1._dp
    IF (PRESENT(alpha)) my_alpha=alpha
    IF (PRESENT(beta)) my_beta=beta
    group = fm%matrix_struct%para_env%group
    context => fm%matrix_struct%context

    CALL dbcsr_get_info(matrix,&
         nblkrows_total=nblkrows_total,&
         nblkcols_total=nblkcols_total,&
         nblkrows_local=nblkrows_local,&
         nblkcols_local=nblkcols_local,&
         local_rows=local_rows,&
         local_cols=local_cols,&
         row_blk_size=rbs,&
         col_blk_size=cbs)
    !write(*,*)'blk size',nblkrows_total, nblkcols_total
    !write(*,*)'blk local size',nblkrows_local, nblkcols_local
    row_blk_sizes => array_data (rbs)
    col_blk_sizes => array_data (cbs)
    ALLOCATE (first_row(nblkrows_total),last_row(nblkrows_total))
    ALLOCATE (first_col(nblkcols_total),last_col(nblkcols_total))
    CALL convert_sizes_to_offsets (row_blk_sizes, first_row, last_row)
    CALL convert_sizes_to_offsets (col_blk_sizes, first_col, last_col)
    mype=context%my_pid
    npe=context%n_pid
    myprow=context%mepos(1)
    mypcol=context%mepos(2)
    nprow=context%num_pe(1)
    npcol=context%num_pe(2)
    nrow_block = fm%matrix_struct%nrow_block
    ncol_block = fm%matrix_struct%ncol_block

    CALL dbcsr_work_create(matrix, nblks_guess=nblkrows_local*nblkcols_local,&
         work_mutable=.TRUE.)

#if defined(__SCALAPACK)
    prow_loop: DO iprow=0,nprow-1
       pcol_loop: DO ipcol=0,npcol-1
          ipe = cp_blacs_pnum(context%group,iprow,ipcol)
          nrow_local = fm%matrix_struct%nrow_locals(iprow)
          ncol_local = fm%matrix_struct%ncol_locals(ipcol)
          IF (ipe /= mype) THEN
             ALLOCATE (fm_block(nrow_local,ncol_local), STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineP,"fm_block",nrow_local*ncol_local)
          ELSE
             fm_block => fm%local_data
          END IF
          IF (nrow_local*ncol_local.NE.0) CALL mp_bcast(fm_block,ipe,group)
          mat_rows: DO row_l = 1, nblkrows_local
             iblock_row = local_rows(row_l)
             mat_cols: DO col_l = 1, nblkcols_local
                iblock_col = local_cols(col_l)
                IF (dbcsr_get_matrix_type(matrix)&
                           .EQ. dbcsr_type_symmetric&
                     .OR. dbcsr_get_matrix_type(matrix)&
                           .EQ. dbcsr_type_antisymmetric&
                     .AND. iblock_row .GT. iblock_col) CYCLE
                CALL dbcsr_get_block_p(matrix, iblock_row, iblock_col,&
                     real_matrix_block, tr, found)
                IF (.NOT. found) THEN
                   ALLOCATE (real_matrix_block(row_blk_sizes(iblock_row),&
                        col_blk_sizes(iblock_col)))
                   real_matrix_block(:,:) = 0.0_dp
                ENDIF
                icol = 1
                bcol_loop: DO icol_global=first_col(iblock_col),last_col(iblock_col)
                   jpcol = indxg2p(icol_global,ncol_block,mypcol,&
                        fm%matrix_struct%first_p_pos(2),npcol)
                   IF (jpcol == ipcol) THEN
                      icol_local = indxg2l(icol_global,ncol_block,mypcol,&
                           fm%matrix_struct%first_p_pos(2),npcol)
                      irow = 1
                      DO irow_global=first_row(iblock_row),last_row(iblock_row)
                         jprow = indxg2p(irow_global,nrow_block,myprow,&
                              fm%matrix_struct%first_p_pos(1),nprow)
                         IF (jprow == iprow) THEN
                            irow_local = indxg2l(irow_global,nrow_block,myprow,&
                                 fm%matrix_struct%first_p_pos(1),nprow)
                            !write(*,*)irow,icol
                            IF (found) THEN
                               real_matrix_block(irow,icol)&
                                    = my_beta*real_matrix_block(irow,icol)+&
                                    my_alpha*fm_block(irow_local,icol_local)
                            ELSE
                               real_matrix_block(irow,icol)&
                                    = my_alpha*fm_block(irow_local,icol_local)
                            ENDIF
                         END IF
                         irow = irow + 1
                      END DO
                   END IF
                   icol = icol + 1
                END DO bcol_loop
                IF (.NOT. found) THEN
                   CALL dbcsr_put_block(matrix,&
                        iblock_row, iblock_col,&
                        real_matrix_block, tr)
                   DEALLOCATE (real_matrix_block)
                ENDIF
             END DO mat_cols
          END DO mat_rows
          IF (ipe /= mype) THEN
             DEALLOCATE (fm_block,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineP,"fm_block")
          END IF
       END DO pcol_loop
    END DO prow_loop

#else
    fm_block => fm%local_data
    DO iblock_row = 1, nblkrows_total
       icol = 1
       DO iblock_col = 1, nblkcols_total
          CALL dbcsr_get_block_p(matrix, iblock_row, iblock_col,&
               real_matrix_block, tr, found)
          IF (.NOT. found) THEN
             ALLOCATE (real_matrix_block(row_blk_sizes(iblock_row),&
                  col_blk_sizes(iblock_col)))
          ENDIF
          DO icol_global=first_col(iblock_col),last_col(iblock_col)
             irow = 1
             DO irow_global=first_row(iblock_row),last_row(iblock_row)
                IF (found) THEN
                   real_matrix_block(irow,icol)&
                        = my_beta*real_matrix_block(irow,icol)+&
                        my_alpha*fm_block(irow_global,icol_global)
                ELSE
                   real_matrix_block(irow,icol)&
                        = my_alpha*fm_block(irow_global,icol_global)
                ENDIF
                irow = irow + 1
             END DO
             icol = icol + 1
          END DO
          IF (.NOT. found) THEN
             CALL dbcsr_put_block(matrix,&
                  iblock_row, iblock_col,&
                  real_matrix_block, tr)
             DEALLOCATE (real_matrix_block)
          ENDIF
       ENDDO
    END DO

#endif

    CALL dbcsr_finalize (matrix)
    CALL timestop(handle)
  END SUBROUTINE copy_fm_to_dbcsr



! *****************************************************************************
!> \brief Create a DBCSR matrix from a real_matrix_type
!> \par Limitations
!>      Data is never moved among processes, so it is more efficient
!>      and gracefully handles any misuse of distribution_2d.
!> \param[out] matrix         bcsr matrix
!> \param[in] sm    cp2k sparse matrix
!> \param[in,out] error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_from_sm(matrix, sm, error)
    TYPE(dbcsr_obj), INTENT(OUT)             :: matrix
    TYPE(real_matrix_type), POINTER          :: sm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_from_sm', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: data_exchange_info_cnt = 3, &
                                                metalen = 2

    INTEGER :: absrow, blk, block_col, mypcol, myprow, nblkcols_local, &
      nblkcols_total, nblkrows_local, nblkrows_total, nfullcols_local, &
      nfullcols_total, nfullrows_local, nfullrows_total, npcols, nprows, nze, &
      row_l, stat, timing_handle
    INTEGER, DIMENSION(:), POINTER :: blkcols_local, blkrows_local, &
      col_blk_begin, col_blk_end, col_blk_size, pcol_dist, prow_dist, &
      row_blk_begin, row_blk_end, row_blk_size
    INTEGER, DIMENSION(:, :), POINTER        :: blacs2mpi
    LOGICAL                                  :: symmetric
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blkdata_2d
    TYPE(array_i1d_obj)                      :: cbs_array, rbs_array
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(distribution_2d_type), POINTER      :: dist2d
    TYPE(real_block_node_type), POINTER      :: sm_block_node

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    dist2d => sm%distribution_2d
    CALL dbcsr_distribution_init (dist)
    CALL dbcsr_dist2d_to_dist (dist2d, dist, error)
    symmetric = (sm%symmetry .EQ. "symmetric")
    nblkrows_total = sm%nblock_row
    nblkcols_total = sm%nblock_col
    row_blk_begin => sm%first_row
    row_blk_end => sm%last_row
    col_blk_begin => sm%first_col
    col_blk_end => sm%last_col
    CALL distribution_2d_get(dist2d, row_distribution=prow_dist,&
         col_distribution=pcol_dist, error=error, blacs_env=blacs_env,&
         n_flat_local_rows=nblkrows_local, n_flat_local_cols=nblkcols_local,&
         n_row_distribution=nblkrows_total, n_col_distribution=nblkcols_total,&
         flat_local_rows=blkrows_local, flat_local_cols=blkcols_local)
    nfullrows_total = MAXVAL(row_blk_end)
    nfullcols_total = MAXVAL(col_blk_end)
    nfullrows_local = 0
    nfullcols_local = 0
    DO row_l = 1, nblkrows_local
       absrow = blkrows_local(row_l)
       nfullrows_local = nfullrows_local +&
            row_blk_end(absrow) - row_blk_begin(absrow) + 1
    ENDDO
    DO row_l = 1, nblkcols_local
       absrow = blkcols_local(row_l)
       nfullcols_local = nfullcols_local +&
            col_blk_end(absrow) - col_blk_begin(absrow) + 1
    ENDDO
    CALL get_blacs_info(blacs_env, para_env=para_env, blacs2mpi=blacs2mpi,&
         my_process_row=myprow, my_process_column=mypcol,&
         number_of_process_rows=nprows, number_of_process_columns=npcols)
    IF (dbg) CALL cp_blacs_env_write(blacs_env, 6, error)
    IF (dbg) WRITE(*,*)routineP//' myprow, mypcol',myprow,'/',nprows,';',&
         mypcol,'/',npcols
    CALL cp_assert(SIZE(prow_dist) .EQ. SIZE(pcol_dist), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Unequal row and column distributions for symmetric matrix.', error)
    CALL cp_assert(MAXVAL(prow_dist).LE.UBOUND(blacs2mpi,1), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Row distribution references unexistent processor rows', error)
    CALL cp_assert(MAXVAL(prow_dist).EQ.UBOUND(blacs2mpi,1), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Range of row distribution not equal to processor rows', error)
    CALL cp_assert(MAXVAL(pcol_dist).LE.UBOUND(blacs2mpi,2), cp_failure_level,&
         cp_wrong_args_error, routineN,&
         'Col distribution references unexistent processor cols', error)
    CALL cp_assert(MAXVAL(pcol_dist).EQ.UBOUND(blacs2mpi,2), cp_warning_level,&
         cp_wrong_args_error, routineN,&
         'Range of col distribution not equal to processor cols', error)
    ALLOCATE (row_blk_size (nblkrows_total), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineP, 'row_blk_size', nblkrows_total)
    ALLOCATE (col_blk_size (nblkcols_total), stat=stat)
    IF (stat /= 0) CALL stop_memory (routineP, 'col_blk_size', nblkcols_total)
    row_blk_size = row_blk_end(:)-row_blk_begin(:)+1
    col_blk_size = col_blk_end(:)-col_blk_begin(:)+1
    !
    CALL array_new (rbs_array, row_blk_size, gift=.TRUE.)
    CALL array_new (cbs_array, col_blk_size, gift=.TRUE.)
    !IF (ASSOCIATED (matrix)) CALL dbcsr_deallocate (matrix, error)
    !ALLOCATE (matrix)
    CALL dbcsr_init (matrix, error)
    CALL dbcsr_create(matrix, sm%name, dist, sm%symmetry,&
         rbs_array, cbs_array, nblks=0, nze=0, data_type=dbcsr_type_real_8,&
         error=error)
    CALL distribution_2d_retain (dist2d, error)
    matrix%m%dist2d => dist2d
    ! Count
    nze = 0
    blk = 0
    DO absrow = 1, nblkrows_total
       sm_block_node => first_block_node(sm, absrow)
       DO WHILE (ASSOCIATED(sm_block_node))
          CALL get_block_node(block_node=sm_block_node, block_col=block_col,&
               block=blkdata_2d)
          ! Only add data if the block is associated
          IF (ASSOCIATED (blkdata_2d)) THEN
             nze = nze + SIZE(blkdata_2d)
             blk = blk + 1
          ENDIF
          sm_block_node => next_block_node(sm_block_node)
       END DO
    ENDDO
    ! Allocate
    CALL dbcsr_work_create(matrix,&
            nblks_guess=blk,&
            sizedata_guess=nze, error=error)
    ! Now copy the data
    DO absrow = 1, nblkrows_total
       sm_block_node => first_block_node(sm, absrow)
       DO WHILE (ASSOCIATED(sm_block_node))
          CALL get_block_node(block_node=sm_block_node, block_col=block_col,&
               block=blkdata_2d)
          ! Only add data if the block is associated
          IF (ASSOCIATED (blkdata_2d)) THEN
             CALL dbcsr_put_block(matrix, absrow, block_col, blkdata_2d)
          ENDIF
          sm_block_node => next_block_node(sm_block_node)
       END DO
    ENDDO
    matrix%m%sparsity_id = sm%sparsity_id
    !
    CALL dbcsr_finalize (matrix)
    CALL dbcsr_distribution_release (dist)
    CALL array_release (rbs_array)
    CALL array_release (cbs_array)
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_from_sm


! *****************************************************************************
!> \brief Create a DBCSR matrix from a real_matrix_type
!> \par Limitations
!>      Data is never moved among processes, so it is more efficient
!>      and gracefully handles any misuse of distribution_2d.
!> \param[out] sm             cp2k sparse matrix
!> \param[in] matrix          DBCSR matrix
!> \param[dist] distribution_2d  distribution_2d for the sparse matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE sm_from_dbcsr (sm, matrix, distribution_2d, error)
    TYPE(real_matrix_type), POINTER          :: sm
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: distribution_2d
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_from_dbcsr', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: data_exchange_info_cnt = 3, &
                                                metalen = 2

    CHARACTER                                :: dbcsr_matrix_type
    CHARACTER(LEN=40)                        :: sm_symmetry
    CHARACTER(LEN=default_string_length)     :: name
    INTEGER :: blk, block_col, nblkcols_total, nblkrows_total, &
      nfullcols_total, nfullrows_total, stored_col, stored_row, timing_handle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sm_first_col, sm_first_row, &
                                                sm_last_col, sm_last_row
    INTEGER, DIMENSION(4)                    :: ver
    LOGICAL                                  :: symmetric, tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blkdata_2d, buff2d
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(distribution_2d_type), POINTER      :: dist2d

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    ! Sanity checks.
    CALL cp_assert (dbcsr_valid_index (matrix), cp_warning_level,&
         cp_caller_error, routineN,&
         "DBCSR Matrix is declared invalid; expect problems.", error=error)
    CALL cp_assert (dbcsr_get_data_type (matrix) .EQ. dbcsr_type_real_8,&
         cp_failure_level, cp_unimplemented_error_nr, routineN,&
         "Non double-precision conversions not possible.", error=error)
    ! Get the dbcsr info and dimensions
    CALL dbcsr_get_info (matrix, nblkrows_total=nblkrows_total,&
         nblkcols_total=nblkcols_total, nfullrows_total=nfullrows_total,&
         nfullcols_total=nfullcols_total, name=name,&
         matrix_type=dbcsr_matrix_type)
    ! Handle dist2d
    IF (PRESENT (distribution_2d)) THEN
       dist2d => distribution_2d
    ELSE
       CALL cp_assert (ASSOCIATED (matrix%m%dist2d), cp_fatal_level,&
            cp_caller_error, routineN,&
            "distribution_2d must be specified or included in DBCSR matrix",&
            error=error)
       dist2d => matrix%m%dist2d
    ENDIF
    ! Get info and verify.
    CALL distribution_2d_get(dist2d, error=error, blacs_env=blacs_env,&
         n_row_distribution=ver(3), n_col_distribution=ver(4))
    CALL cp_assert (nblkrows_total.EQ.ver(3).AND.nblkcols_total.EQ.ver(4),&
         cp_warning_level, cp_internal_error, routineN,&
         "Specified matrix dimensions/distributions seem incompatible.")
    ! Convert type
    symmetric = .FALSE.

    SELECT CASE (dbcsr_matrix_type)
    CASE (dbcsr_type_normal)
       sm_symmetry = "none"
    CASE (dbcsr_type_symmetric)
       sm_symmetry = "symmetric"
       symmetric = .TRUE.
    CASE (dbcsr_type_antisymmetric)
       sm_symmetry = "antisymmetric"
       symmetric = .TRUE.
    CASE default
       CALL cp_unimplemented_error (routineN,&
            "Can not convert this matrix type or invalid matrix type")
    END SELECT
    ! Get row/column offsets.
    ALLOCATE (sm_first_row(nblkrows_total), sm_last_row(nblkrows_total),&
         sm_first_col(nblkcols_total), sm_last_col(nblkcols_total))
    CALL convert_sizes_to_offsets (array_data (dbcsr_row_block_sizes (matrix)),&
         sm_first_row, sm_last_row)
    CALL convert_sizes_to_offsets (array_data (dbcsr_col_block_sizes (matrix)),&
         sm_first_col, sm_last_col)
    CALL allocate_matrix(matrix=sm,&
                         nblock_row=nblkrows_total,&
                         nblock_col=nblkcols_total,&
                         nrow=nfullrows_total,&
                         ncol=nfullcols_total,&
                         first_row=sm_first_row,&
                         last_row=sm_last_row,&
                         first_col=sm_first_col,&
                         last_col=sm_last_col,&
                         matrix_name=TRIM(name),&
                         matrix_symmetry=sm_symmetry,&
                         sparsity_id=matrix%m%sparsity_id,&
                         distribution_2d=dist2d, error=error,&
                         keep_pointer=.TRUE.)
    DEALLOCATE (sm_first_row, sm_last_row, sm_first_col, sm_last_col)
    !
    NULLIFY (buff2d)
    CALL dbcsr_iterator_start(iter, matrix)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       NULLIFY(blkdata_2d)
       CALL dbcsr_iterator_next_block(iter, stored_row, stored_col,&
            blkdata_2d, tr, blk)
       ! We need to get the upper triangular form for the
       ! sparse_matrix_type.
       IF (symmetric .AND. stored_row .GT. stored_col) THEN
          CALL swap (stored_row, stored_col)
          tr = .NOT. tr
       ENDIF
       IF (.NOT. tr) THEN
          CALL add_block_node(matrix=sm,&
               block_row=stored_row, block_col=stored_col,&
               block=blkdata_2d, error=error)
       ELSE
          ! This should never be the case when converting a
          ! sparse_matrix_type-compatible matrix back to the
          ! sparse_matrix_type, but can be the case with other
          ! matrices.
          !CALL cp_assert (.FALSE., cp_warning_level, cp_internal_error,&
          !     routineN, "Should not need to transpose")
          ALLOCATE (buff2d(SIZE(blkdata_2d,1),SIZE(blkdata_2d,2)))
          buff2d = TRANSPOSE(blkdata_2d)
          CALL add_block_node(matrix=sm,&
               block_row=stored_row, block_col=stored_col,&
               block=buff2d, error=error)
          DEALLOCATE (buff2d)
       ENDIF
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    CALL timestop(timing_handle)
  END SUBROUTINE sm_from_dbcsr


  SUBROUTINE dbcsr_write_sparse_matrix(sparse_matrix,before,after,qs_env,para_env,&
                                      first_row,last_row,first_col,last_col,scale,&
                                       output_unit,error)

    TYPE(dbcsr_obj), POINTER                 :: sparse_matrix
    INTEGER, INTENT(IN)                      :: before, after
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_col, last_col
    REAL(dp), INTENT(IN), OPTIONAL           :: scale
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_write_sparse_matrix', &
      routineP = moduleN//':'//routineN

    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(real_matrix_type), POINTER          :: sparse_matrix_sm

    NULLIFY(sparse_matrix_sm,distribution_2d)

    CALL get_qs_env(qs_env=qs_env,distribution_2d=distribution_2d,error=error)

    CALL sm_from_dbcsr ( sparse_matrix_sm, sparse_matrix, distribution_2d, error=error )

    CALL write_sparse_matrix(sparse_matrix_sm,before,after,qs_env,para_env,&
                             first_row,last_row,first_col,last_col,scale,&
                             output_unit,error)

    CALL deallocate_matrix(sparse_matrix_sm, error)

  END SUBROUTINE dbcsr_write_sparse_matrix

  SUBROUTINE dbcsr_write_matrix_distribution(matrix,force_env_section,para_env,error)
    TYPE(dbcsr_obj), POINTER                 :: matrix
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'dbcsr_write_matrix_distribution', &
      routineP = moduleN//':'//routineN

    TYPE(real_matrix_type), POINTER          :: matrix_sm

    NULLIFY(matrix_sm)

    CALL sm_from_dbcsr ( matrix_sm, matrix, matrix%m%dist2d, error=error )

    CALL write_matrix_distribution(matrix_sm,force_env_section,para_env,error)

    CALL deallocate_matrix(matrix_sm, error)

  END SUBROUTINE dbcsr_write_matrix_distribution

  SUBROUTINE test_sm(sm, error)
    TYPE(real_matrix_type), POINTER          :: sm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'test_sm', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: data_exchange_info_cnt = 3, &
                                                metalen = 2

    INTEGER                                  :: blk, block_col, icol, irow, &
                                                proc, timing_handle
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blkdata_2d
    TYPE(real_block_node_type), POINTER      :: block_node

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)

    proc = sm%distribution_2d%blacs_env%para_env%mepos

    DO irow=1,sm%nblock_row
       WRITE(*,*) '1proc',proc,' irow',irow
       block_node => first_block_node(sm,irow)
       DO WHILE (ASSOCIATED(block_node))
          CALL get_block_node(block_node=block_node,block_col=icol,BLOCK=blkdata_2d)

          WRITE(*,*) '1sm_test: proc',proc,' irow',irow,' icol',icol,'blk',blkdata_2d(1,1),&
               SIZE(blkdata_2d,1),SIZE(blkdata_2d,2)
          block_node => next_block_node(block_node)
       END DO
    END DO

    DO irow = 1,sm%nblock_row
       DO icol = irow,sm%nblock_row
          WRITE(*,*) '2proc',proc,' irow',irow,' icol',icol
          CALL get_block_node(sm,block_row=irow,block_col=icol,BLOCK=blkdata_2d)
          IF(ASSOCIATED(blkdata_2d)) THEN
             WRITE(*,*) '2sm_test: proc',proc,' irow',irow,' icol',icol,'blk',blkdata_2d(1,1)
          ENDIF
       ENDDO
    ENDDO

    CALL timestop(timing_handle)
  END SUBROUTINE test_sm


END MODULE cp_dbcsr_operations

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
                                             array_release
  USE cp_fm_basic_linalg,              ONLY: cp_fm_frobenius_norm
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE dbcsr_methods,                   ONLY: &
       dbcsr_col_block_sizes, dbcsr_destroy, dbcsr_distribution, &
       dbcsr_distribution_col_dist, dbcsr_distribution_mp, &
       dbcsr_distribution_new, dbcsr_distribution_release, &
       dbcsr_distribution_row_dist, dbcsr_get_data, &
       dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
       dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_mp_group, &
       dbcsr_mp_mynode, dbcsr_mp_mypcol, dbcsr_mp_myprow, dbcsr_mp_npcols, &
       dbcsr_mp_nprows, dbcsr_mp_numnodes, dbcsr_mp_pgrid, &
       dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_nfullcols_local, &
       dbcsr_nfullrows_local, dbcsr_row_block_sizes
  USE dbcsr_operations,                ONLY: &
       block_set_d, dbcsr_add, dbcsr_copy, dbcsr_find_column, &
       dbcsr_mult_e_ec_m, dbcsr_put_block, dbcsr_remove_block, dbcsr_scale, &
       dbcsr_trace
  USE dbcsr_types,                     ONLY: dbcsr_array_type,&
                                             dbcsr_distribution_obj,&
                                             dbcsr_imagedistribution_type,&
                                             dbcsr_iterator,&
                                             dbcsr_mp_obj,&
                                             dbcsr_obj,&
                                             dbcsr_type
  USE dbcsr_util,                      ONLY: &
       dbcsr_checksum, dbcsr_create, dbcsr_create_col_image_dist, &
       dbcsr_create_dist_r_rot, dbcsr_create_dist_r_unrot, &
       dbcsr_destroy_array, dbcsr_destroy_image_dist, dbcsr_finalize, &
       dbcsr_from_fm, dbcsr_from_sm, dbcsr_init, dbcsr_make_col_images, &
       dbcsr_new_transposed, dbcsr_set_debug, dbcsr_verify_matrix, &
       sm_from_dbcsr, dbcsr_print
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_retain,&
                                             distribution_2d_type
  USE kinds,                           ONLY: dp,&
                                             int_8,&
                                             sp
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum_i8,&
                                             mp_sync
  USE sparse_matrix_types,             ONLY: checksum_matrix,&
                                             deallocate_matrix,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_operations'

  PUBLIC :: dbcsr_sm_sm_multiply, dbcsr_check
  PUBLIC :: cp_dbcsr_add_dist2d, cp_dbcsr_get_dist2d, cp_dbcsr_remove_dist2d
  PUBLIC :: dbcsr_multiply_vec, dbcsr_pack_vec,&
       packed_vec_scale, dbcsr_unpack_vec,&
       add_packed_vecs, dbcsr_mult_pack_vec, packed_vec_bif_tech, packed_vec_bif_tech2, dbcsr_copy_vec, packed_vec_build_u

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
    ! We don't recycle yet.
    IF (ASSOCIATED (matrix_c%m)) THEN
       CALL dbcsr_destroy (matrix_c, error=error)
    ENDIF
    CALL dbcsr_init (matrix_c, error)
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
    CALL dbcsr_destroy_array (ms_left, error)
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
    LOGICAL                                  :: failure, ionode, ldbg
    REAL(kind=dp)                            :: checksum, newtime, oldtime, &
                                                trace, walltime1, walltime2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block2d, block2d_2
    TYPE(array_i1d_obj)                      :: col_blk_size_right
    TYPE(cp_fm_type), POINTER                :: fmout
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_distribution_obj)             :: dist_right, product_dist
    TYPE(dbcsr_iterator)                     :: iter
    TYPE(dbcsr_obj)                          :: matrix_a, matrix_a2, &
                                                matrix_b, matrix_bt, matrix_c
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
    CALL dbcsr_destroy(matrix_a2,error)

    ! Fill a matrix with all zeros.
    CALL dbcsr_init (matrix_a2, error)
    CALL dbcsr_from_sm(matrix_a2, matrix_s, error)
    !write(*,*)'zeroing...'
    rbs => array_data (dbcsr_row_block_sizes (matrix_a2))
    cbs => array_data (dbcsr_col_block_sizes (matrix_a2))
    ALLOCATE (block2d(MAXVAL (rbs), MAXVAL (cbs)))
    block2d(:,:) = 0
    DO row = 1, dbcsr_nblkrows_total (matrix_a2)
       rowsize = rbs (row)
       DO col = 1, dbcsr_nblkcols_total (matrix_a2)
          colsize = cbs (col)
          CALL dbcsr_put_block (matrix_a2,&
               block2d(1:rowsize, 1:colsize), row, col)
       ENDDO
    ENDDO
    DEALLOCATE (block2d)
    CALL dbcsr_finalize (matrix_a2, error=error)
    CALL dbcsr_verify_matrix(matrix_a2, error)
    checksum = dbcsr_checksum(matrix_a2, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR 0 checksum ", checksum
    ENDIF
    ! Now delete all the blocks.
    !WRITE(*,*)'Deleting all blocks!'
    CALL dbcsr_iterator_start (iter, matrix_a2)
    DO WHILE (dbcsr_iterator_blocks_left (iter))
       CALL dbcsr_iterator_next_block (iter, row, col)
       CALL dbcsr_remove_block (matrix_a2, row, col)
    ENDDO
    !WRITE(*,*)"blk_p=", matrix_a2%blk_p
    CALL dbcsr_finalize (matrix_a2, error=error)
    !WRITE(*,*)"blk_p=", matrix_a2%nblks
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR 0 checksum2", checksum
    ENDIF
    CALL dbcsr_iterator_stop (iter)
    CALL dbcsr_destroy (matrix_a2, error)
    !
    

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
       CALL dbcsr_iterator_next_block (iter, row, col, block2d)
    !write(*,*)'got block',row,col
       CALL dbcsr_put_block (matrix_b, block2d, row, col)
       CALL dbcsr_verify_matrix (matrix_b, error, local=.TRUE.)
       !CALL dbcsr_print (matrix_b, error)
       checksum = dbcsr_checksum(matrix_b, error, local=.TRUE.)
       IF (ldbg.AND.ionode) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x checksum2", checksum
       ENDIF
       ALLOCATE (block2d_2(SIZE(block2d,1),SIZE(block2d,2)))
       block2d_2(:,:) = block2d(:,:)
       block2d(:,:) = 0.0_dp
       CALL dbcsr_put_block (matrix_b, block2d, row, col)
       CALL dbcsr_finalize (matrix_b, error=error)
       CALL dbcsr_verify_matrix (matrix_b, error, local=.TRUE.)
       checksum = dbcsr_checksum(matrix_b, error, local=.TRUE.)
       IF (ldbg.AND.ionode) THEN
          WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x wrng cksm", checksum
       ENDIF
       !CALL dbcsr_print (matrix_b, error)
       CALL dbcsr_put_block (matrix_b, block2d_2, row, col)
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
       CALL dbcsr_put_block (matrix_b, block2d_2, row, col)
       CALL dbcsr_finalize (matrix_b, error=error)
       CALL dbcsr_verify_matrix (matrix_b, error=error, local=.TRUE.)
       !CALL dbcsr_print (matrix_b, error)
       DEALLOCATE (block2d_2)
    ENDIF
    CALL dbcsr_iterator_stop (iter)
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
    !
    CALL dbcsr_new_transposed (matrix_bt, matrix_b, error=error)
    CALL dbcsr_verify_matrix (matrix_bt, error)
    !CALL dbcsr_destroy (matrix_bt, error)
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
    newtime = 0.0_dp
    newtime = MAX (walltime2-walltime1, 1.0_dp)
    CALL dbcsr_destroy (matrix_bt, error)
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
    ! check trace(c'*s*c)
    CALL dbcsr_trace(matrix_b, matrix_c, trace, error=error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4,A,I4)')"DBCSR trace(C'*S*C)", trace," should be equal to ",k
    ENDIF

    CALL dbcsr_destroy (matrix_b, error)
    CALL dbcsr_init (matrix_b, error)
    CALL dbcsr_destroy (matrix_c, error)

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
    CALL array_release (col_blk_size_right)
    CALL dbcsr_destroy (matrix_a,error)
    CALL dbcsr_destroy (matrix_b, error)
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
          CALL dbcsr_find_column(col,frst_blk,last_blk,matrix%m%col_i,blk,found)
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
          CALL dbcsr_find_column(col,frst_blk,last_blk,matrix%m%col_i,blk,found)
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
          CALL dbcsr_find_column(col,frst_blk,last_blk,matrix%m%col_i,blk,found)
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
    CASE ('r')
       CALL dbcsr_get_data (matrix_a%m%data_area, a_data_r)
       CALL dbcsr_get_data (matrix_b%m%data_area, b_data_r)
    CASE ('R')
       CALL dbcsr_get_data (matrix_a%m%data_area, a_data_d)
       CALL dbcsr_get_data (matrix_b%m%data_area, b_data_d)
    CASE ('c')
       CALL dbcsr_get_data (matrix_a%m%data_area, a_data_c)
       CALL dbcsr_get_data (matrix_b%m%data_area, b_data_c)
    CASE ('C')
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
          CALL dbcsr_find_column(b_col, a_frst_blk, a_last_blk, matrix_a%m%col_i, a_blk, found)
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
             CASE ('R')
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

END MODULE cp_dbcsr_operations

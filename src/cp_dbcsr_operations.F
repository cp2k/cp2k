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
  USE array_types,                     ONLY: array_i1d_obj,&
                                             array_release
  USE cp_fm_basic_linalg,              ONLY: cp_fm_frobenius_norm
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE dbcsr_methods,                   ONLY: dbcsr_distribution_col_dist,&
                                             dbcsr_distribution_mp,&
                                             dbcsr_distribution_new,&
                                             dbcsr_distribution_release,&
                                             dbcsr_distribution_row_dist,&
                                             dbcsr_mp_group,&
                                             dbcsr_nfullcols_local,&
                                             dbcsr_nfullrows_local
  USE dbcsr_operations,                ONLY: dbcsr_add,&
                                             dbcsr_copy,&
                                             dbcsr_mult_e_ec_m,&
                                             dbcsr_scale,&
                                             dbcsr_trace
  USE dbcsr_types,                     ONLY: dbcsr_array_type,&
                                             dbcsr_distribution_obj,&
                                             dbcsr_imagedistribution_type,&
                                             dbcsr_type
  USE dbcsr_util,                      ONLY: &
       dbcsr_checksum, dbcsr_create, dbcsr_create_col_image_dist, &
       dbcsr_create_dist_r_rot, dbcsr_create_dist_r_unrot, dbcsr_destroy, &
       dbcsr_destroy_array, dbcsr_destroy_image_dist, dbcsr_from_fm, &
       dbcsr_from_sm, dbcsr_init, dbcsr_make_col_images, &
       dbcsr_new_transposed, dbcsr_set_debug, dbcsr_verify_matrix, &
       dbcsr_work_create
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_retain,&
                                             distribution_2d_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: mp_sum_i8,&
                                             mp_sync
  USE sparse_matrix_types,             ONLY: checksum_matrix,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_operations'

  PUBLIC :: dbcsr_sm_sm_multiply, dbcsr_check
  PUBLIC :: cp_dbcsr_add_dist2d, cp_dbcsr_get_dist2d, cp_dbcsr_remove_dist2d

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
    TYPE(dbcsr_type), INTENT(IN)             :: matrix_a, matrix_b
    TYPE(dbcsr_type), INTENT(OUT)            :: matrix_c
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: distribution
    INTEGER, INTENT(IN)                      :: hint
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                               :: flop

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sm_sm_multiply', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, timing_handle
    INTEGER(KIND=int_8)                      :: my_flop
    REAL(KIND=dp)                            :: checksum, flops
    TYPE(dbcsr_array_type)                   :: ms_left
    TYPE(dbcsr_imagedistribution_type)       :: rdist_left
    TYPE(dbcsr_type)                         :: matrix_d

!TYPE(dbcsr_distribution_obj)             :: dist_left, dist_right

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY
    CALL timeset(routineN, timing_handle)
    !> pre The matrices must be valid and they must have compatible shapes.
    ! Assign new dist2d
    ! Create & initialize the final & working matrix.
    IF (dbg) WRITE(*,*)routineP//'Creating 2-d distribution...'
    IF (dbg) CALL mp_sync(dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))
    CALL dbcsr_init (matrix_c, error)
    CALL dbcsr_create(matrix_c, "C", distribution, 'N',&
         matrix_a%row_blk_size, matrix_b%col_blk_size,&
         0, 0, error=error)
    IF (dbg) WRITE(*,*)routineP//'Updating new matrix...'
    IF (dbg) CALL mp_sync(dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))
    !IF (dbg) WRITE(*,*)routineP//' C initial index is', matrix_c%index
    IF (dbg) WRITE(*,*)routineP//'Creating work matrix...'
    IF (dbg) CALL mp_sync(dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))
    i = dbcsr_nfullrows_local(matrix_a) * dbcsr_nfullcols_local(matrix_b)
    CALL dbcsr_work_create(matrix_c, sizedata_guess=i, error=error)
    IF (dbg) WRITE(*,*)routineP//'Multiplying...'
    IF (dbg) CALL mp_sync(dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))
    !
    DBG 'Creating left image distribution'
    !CALL dbcsr_create_image_dist_left(rdist_left, matrix_a%dist, matrix_b%dist)
    CALL dbcsr_create_col_image_dist (rdist_left,&
         dbcsr_distribution_row_dist (matrix_a%dist), matrix_b%dist)
    !
    DBG 'Making A images'
    CALL dbcsr_make_col_images(matrix_a, ms_left, rdist_left, error)
    !ms_left%image_dist = rdist_left
    matrix_d = ms_left%mats(1)
    DBG 'master name', ms_left%mats(1)%name
    DBG "Verifying A images"
    DO i = 1, SIZE (ms_left%mats)
       CALL dbcsr_verify_matrix(ms_left%mats(i), error)
    ENDDO
    IF (PRESENT (flop)) flop = 0.0
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
    !CALL mp_sum(flop, dbcsr_mp_group (dbcsr_distribution_mp (matrix_a%dist)))
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

    INTEGER                                  :: handle, k, output_unit
    INTEGER(KIND=int_8)                      :: flop
    LOGICAL                                  :: failure, ionode, ldbg
    REAL(kind=dp)                            :: checksum, newtime, oldtime, &
                                                trace, walltime1, walltime2
    TYPE(array_i1d_obj)                      :: col_blk_size_right
    TYPE(cp_fm_type), POINTER                :: fmout
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_distribution_obj)             :: dist_right, product_dist
    TYPE(dbcsr_type)                         :: matrix_a, matrix_a2, &
                                                matrix_b, matrix_bt, matrix_c

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

    checksum = 0.0_dp
    CALL cp_fm_frobenius_norm(matrix_x, checksum, error)
    !WRITE(*,*)'x checksum',checksum*checksum
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"   FM x  checksum", checksum*checksum
    ENDIF
    CALL cp_fm_get_info(matrix_x, ncol_global=k, error=error)
    
    IF (MOD(hint,2) .EQ. 0) THEN
       CALL dbcsr_create_dist_r_unrot (dist_right, matrix_a%dist, k,&
            col_blk_size_right)
    ELSE
       CALL dbcsr_create_dist_r_rot (dist_right, matrix_a%dist, k,&
            col_blk_size_right)
    ENDIF
    CALL dbcsr_init (matrix_b, error)
    CALL dbcsr_from_fm(matrix_b, matrix_x, 0.0_dp, dist_right,&
         row_blk_size=matrix_a%row_blk_size,&
         col_blk_size=col_blk_size_right,&
         error=error)
    CALL dbcsr_verify_matrix(matrix_b, error)

    checksum = dbcsr_checksum(matrix_b, error)
    !WRITE(*,*)'DBCSR x checksum',checksum
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR x  checksum", checksum
    ENDIF
    !
    ! test the scale routine on a general matrix
    CALL dbcsr_scale(matrix_b,alpha_scalar=3.0_dp,error=error)
    checksum = 0.0_dp
    checksum = dbcsr_checksum(matrix_b, error)
    IF (ldbg.AND.ionode) THEN
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR 3*x  checksum", checksum
    ENDIF
    CALL dbcsr_scale(matrix_b,alpha_scalar=1.0_dp/3.0_dp,error=error)

   !
    CALL dbcsr_new_transposed (matrix_bt, matrix_b, error=error)
    CALL dbcsr_verify_matrix (matrix_bt, error)
    !CALL dbcsr_destroy (matrix_bt, error)
    !
    CALL dbcsr_distribution_new (product_dist,&
         dbcsr_distribution_mp (matrix_a%dist),&
         dbcsr_distribution_row_dist (matrix_a%dist),&
         dbcsr_distribution_col_dist (matrix_b%dist))
    !
    walltime1 = m_walltime()
    flop = 0
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
    CALL mp_sum_i8 (flop, dbcsr_mp_group(dbcsr_distribution_mp(matrix_c%dist)))
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
       WRITE(output_unit, '(T2,A,F12.4)')"DBCSR trace(C'*S*C)", trace
    ENDIF

    CALL dbcsr_destroy (matrix_b, error)
    CALL dbcsr_init (matrix_b, error)
    CALL dbcsr_destroy (matrix_c, error)

    CALL timestop(handle)

    CALL timeset(routineN//'std', handle)

    CALL cp_fm_get_info(matrix_x, ncol_global=k, error=error)
    CALL cp_fm_create(fmout, matrix_x%matrix_struct, "FMout", error)
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
         row_blk_size=matrix_a%row_blk_size,&
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

END MODULE cp_dbcsr_operations

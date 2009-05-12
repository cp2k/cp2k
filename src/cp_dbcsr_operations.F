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
!> - Created
! *****************************************************************************
MODULE cp_dbcsr_operations
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_frobenius_norm
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE dbcsr_methods,                   ONLY: dbcsr_distribution_release
  USE dbcsr_operations,                ONLY: dbcsr_multiply_sets
  USE dbcsr_types,                     ONLY: dbcsr_array_type,&
                                             dbcsr_distribution_obj,&
                                             dbcsr_rptdistribution_type,&
                                             dbcsr_type,&
                                             dbcsr_work_type
  USE dbcsr_util,                      ONLY: &
       dbcsr_checksum, dbcsr_create, dbcsr_create_dist2d_r_unrot, &
       dbcsr_create_dist_left, dbcsr_create_rdist, dbcsr_destroy, &
       dbcsr_destroy_array, dbcsr_destroy_rdist, dbcsr_dist2d_to_dist, &
       dbcsr_from_fm, dbcsr_from_sm, dbcsr_normalize, dbcsr_print, &
       dbcsr_set_debug, dbcsr_verify_matrix, dbcsr_work_create
  USE distribution_2d_types,           ONLY: distribution_2d_get,&
                                             distribution_2d_release,&
                                             distribution_2d_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE sparse_matrix_types,             ONLY: checksum_matrix,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_dbcsr_operations'

  PUBLIC :: dbcsr_sm_sm_multiply, dbcsr_check

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
!> \param[in] dist2d          optional distribution_2d for the product matrix.
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_sm_sm_multiply(matrix_a, matrix_b, matrix_c, dist2d,&
       error)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix_a, matrix_b
    TYPE(dbcsr_type), INTENT(OUT)            :: matrix_c
    TYPE(distribution_2d_type), POINTER      :: dist2d
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sm_sm_multiply', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: flop, i, nblkcols_local, &
                                                nblkrows_local, timing_handle
    REAL(KIND=dp)                            :: checksum, flops
    TYPE(dbcsr_array_type)                   :: ms_left, ms_right
    TYPE(dbcsr_distribution_obj)             :: dist_left, dist_right
    TYPE(dbcsr_rptdistribution_type)         :: rdist_left, rdist_right
    TYPE(dbcsr_type)                         :: matrix_d, matrix_e
    TYPE(dbcsr_work_type)                    :: wm

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODYY
    CALL timeset(routineN, timing_handle)
    !> pre The matrices must be valid and they must have compatible shapes.
    ! Assign new dist2d
    ! Create & initialize the final & working matrix.
    IF (dbg) WRITE(*,*)routineP//'Creating 2-d distribution...'
    CALL mp_sync(matrix_a%dist2d%blacs_env%para_env%group)
    CALL distribution_2d_get(matrix_b%dist2d, error=error,&
         n_flat_local_rows=nblkrows_local, n_flat_local_cols=nblkcols_local)
    CALL dbcsr_create(matrix_c, "C", matrix_b%dist2d, 'N',&
         (/0, 0,&
         matrix_b%nblkrows_total, matrix_b%nblkcols_total,&
         matrix_b%nfullrows_total, matrix_b%nfullrows_total,&
         matrix_b%nblkrows_local, matrix_b%nblkcols_local,&
         matrix_b%nfullrows_local, matrix_b%nfullcols_local/),&
         0, 0, error=error)
    IF (dbg) WRITE(*,*)routineP//'Updating new matrix...'
    CALL mp_sync(matrix_a%dist2d%blacs_env%para_env%group)
    ALLOCATE(matrix_c%row_blk_size(matrix_a%nblkrows_total),&
         matrix_c%col_blk_size(matrix_b%nblkcols_total))
    matrix_c%row_blk_size(:) =&
         matrix_a%row_blk_size(:)
    matrix_c%col_blk_size(:) =&
         matrix_b%col_blk_size(:)
    !IF (dbg) WRITE(*,*)routineP//' C initial index is', matrix_c%index
    IF (dbg) WRITE(*,*)routineP//'Creating work matrix...'
    CALL mp_sync(matrix_a%dist2d%blacs_env%para_env%group)
    CALL dbcsr_work_create(wm, matrix_c,&
         sizedata_guess=INT(1.5*SIZE(matrix_a%data)), error=error)
    IF (dbg) WRITE(*,*)routineP//'Multiplying...'
    CALL mp_sync(matrix_a%dist2d%blacs_env%para_env%group)
    !CALL dbcsr_sm2_mult_overlap(matrix_a, matrix_b, wm, error)
    !CALL dbcsr_sm2_mult_general(matrix_a, matrix_b, wm, error)
    !CALL dbcsr_sm2_mult_ls_rn(matrix_a, matrix_b, wm, error)
    CALL dbcsr_print(matrix_a, error)
    DBG 'Creating repetition distribution'
    CALL dbcsr_dist2d_to_dist(matrix_a%dist2d, dist_left, error)
    CALL dbcsr_dist2d_to_dist(matrix_b%dist2d, dist_right, error)
    CALL dbcsr_create_rdist (rdist_right, dist_right)
    CALL dbcsr_create_dist_left(rdist_left, dist_left, dist_right)
    DBG 'Normalizing A'
    CALL dbcsr_normalize(matrix_a, ms_left, rdist_left, error)
    ms_left%rdist = rdist_left
    matrix_d = ms_left%mats(1)
    DBG 'master name', ms_left%mats(1)%name
    DBG "Verifying A->D"
    CALL dbcsr_verify_matrix(matrix_d, error)
    DO i = 1, SIZE (ms_left%mats)
       CALL dbcsr_print (ms_left%mats(i), error)
    ENDDO
    flops = dbcsr_checksum(matrix_d, error)
    DBG "DBCSR norm. S checksum", flops
    CALL dbcsr_print (matrix_b, error)
    DBG 'Normalizing B'
    CALL dbcsr_normalize(matrix_b, ms_right, rdist_right, error)
    ms_right%mats(1)%dist = dist_right
    matrix_e = ms_right%mats(1)
    DBG "Verifying B->E"
    CALL dbcsr_verify_matrix(matrix_e, error)
    CALL dbcsr_print(matrix_e, error)
    flops = dbcsr_checksum(matrix_e, error)
    DBG "DBCSR norm. x checksum", flops
    !DBG "Matrix A"
    !CALL dbcsr_print(matrix_a, error)
    !DBG "Matrix D"
    !CALL dbcsr_print(matrix_d, error)
    !DBG "Matrix B"
    !CALL dbcsr_print(matrix_b, error)
    !DBG "Matrix E"
    !CALL dbcsr_print(matrix_e, error)
    DBG 'Multiplying normalized matrices'
    !DO i = 1, SIZE (ms_left%mats)
    !   CALL dbcsr_nn_mult(ms_left%mats(i), matrix_e, matrix_c, wm, error)
    !ENDDO
    CALL dbcsr_dist2d_to_dist (matrix_c%dist2d, matrix_c%dist, error)
    !
    CALL dbcsr_multiply_sets(ms_left, ms_right, matrix_c, wm, error)
    !CALL dbcsr_multiply_sets_sqrt(ms_left, ms_right, matrix_c, wm, error)
    !
    !IF (dbg) WRITE(*,*)routineP//' C prefinal index is', matrix_c%index
    !CALL dbcsr_finalize(wm, matrix_c, .TRUE., .FALSE., error)
    !IF (dbg) WRITE(*,*)routineP//' C final   index is', matrix_c%index
    IF (dbg) THEN
       CALL dbcsr_verify_matrix(matrix_c, error)
       checksum = dbcsr_checksum(matrix_c, error)
       WRITE(*,*)routineP//' BCSR Sx Matrix checksum is', checksum
    ENDIF
    flop = 0
    CALL mp_sum(flop, matrix_a%dist2d%blacs_env%para_env%group)
    WRITE(*,*)'FLOP=',flop
    DBG "Result matrix"
    CALL dbcsr_print(matrix_c, error)
    !
    CALL dbcsr_distribution_release (matrix_c%dist)
    CALL dbcsr_distribution_release (dist_left)
    CALL dbcsr_distribution_release (dist_right)
    CALL dbcsr_destroy_rdist (rdist_left)
    CALL dbcsr_destroy_rdist (rdist_right)
    CALL dbcsr_destroy_array (ms_left, error)
    CALL dbcsr_destroy_array (ms_right, error)
    ! Perform multiplication
    CALL timestop(timing_handle)
  END SUBROUTINE dbcsr_sm_sm_multiply



! *****************************************************************************
!> \brief   Checks the new bcsr routines, to find bugs early
!> \author  Urban
!> \date    2008-10-07
! *****************************************************************************
  SUBROUTINE dbcsr_check(matrix_s, matrix_x, error)

    TYPE(real_matrix_type), POINTER          :: matrix_s
    TYPE(cp_fm_type), POINTER                :: matrix_x
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, k
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: dbg = .FALSE., failure
    REAL(kind=dp)                            :: checksum, newtime, oldtime, &
                                                walltime1, walltime2
    TYPE(cp_blacs_env_type), POINTER         :: tmp_blacs_env
    TYPE(cp_fm_type), POINTER                :: fmout
    TYPE(dbcsr_type)                         :: matrix_a, matrix_b, matrix_c
    TYPE(distribution_2d_type), POINTER      :: dist2d, dist2d_left

    CALL timeset(routineN,handle)

    failure=.FALSE.

    !Do the new multiplication
    dist2d_left => matrix_s%distribution_2d
    IF (dbg) WRITE(*,*)'left-1',dist2d_left%ref_count,dist2d_left%blacs_env%ref_count,&
         dist2d_left%blacs_env%para_env%ref_count
    CALL dbcsr_from_sm(matrix_a, matrix_s, error)
    CALL dbcsr_verify_matrix(matrix_a, error)
    WRITE(*,*)'S checksum=',&
         checksum_matrix(matrix_s, matrix_s%distribution_2d%blacs_env%para_env)
    checksum = 0.0
    checksum = dbcsr_checksum(matrix_a, error)
    WRITE(*,*)'BCSR S checksum',checksum
    !row_blk_begin => matrix_a%row_blk_begin
    !row_blk_end => matrix_a%row_blk_end
    CALL cp_fm_frobenius_norm(matrix_x, checksum, error)
    CALL cp_fm_get_info(matrix_x, ncol_global=k, error=error)
    CALL dbcsr_create_dist2d_r_unrot(dist2d, dist2d_left, k,&
         col_blk_size, error=error)
    IF (dbg) WRITE(*,*)'right0=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
         dist2d%blacs_env%para_env%ref_count
    CALL dbcsr_from_fm(matrix_b, matrix_x, 0.0d0, dist2d,&
         row_blk_size=matrix_a%row_blk_size,&
         col_blk_size=col_blk_size,&
         error=error)
    IF (dbg) WRITE(*,*)'right1=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
         dist2d%blacs_env%para_env%ref_count
    CALL dbcsr_verify_matrix(matrix_b, error)
    WRITE(*,*)'x checksum',checksum*checksum
    checksum = dbcsr_checksum(matrix_b, error)
    WRITE(*,*)'DBCSR x checksum',checksum
    !
    walltime1 = m_walltime()
    !!CALL dbcsr_desymmetrize(matrix_a, matrix_afull, error)
    CALL dbcsr_sm_sm_multiply(matrix_a, matrix_b, matrix_c, dist2d, error)
    IF (dbg) WRITE(*,*)'right2=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
         dist2d%blacs_env%para_env%ref_count
    walltime2 = m_walltime()
    newtime = walltime2-walltime1
    !
    CALL dbcsr_verify_matrix(matrix_c, error)
    !
    checksum = dbcsr_checksum(matrix_c, error)
    WRITE(*,*)'BCSR Sx checksum',checksum,'time=',newtime
    CALL dbcsr_destroy(matrix_b,error)
    CALL dbcsr_destroy(matrix_c,error)
    IF (dbg) WRITE(*,*)'right3=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
         dist2d%blacs_env%para_env%ref_count
    tmp_blacs_env => dist2d%blacs_env
    IF (dbg) WRITE(*,*)'right4=','"',0,'"',tmp_blacs_env%ref_count,&
         tmp_blacs_env%para_env%ref_count
    CALL cp_blacs_env_release(tmp_blacs_env, error)
    IF (dbg) WRITE(*,*)'left99',dist2d_left%ref_count,dist2d_left%blacs_env%ref_count,&
         dist2d_left%blacs_env%para_env%ref_count
    !WRITE(*,*)'right4=',dist2d%ref_count,dist2d%blacs_env%ref_count,&
    !     dist2d%blacs_env%para_env%ref_count

    CALL timestop(handle)

    CALL timeset(routineN//'std', handle)

    CALL cp_fm_get_info(matrix_x, ncol_global=k, error=error)
    CALL cp_fm_create(fmout, matrix_x%matrix_struct, "FMout", error)
    walltime1 = m_walltime()
    CALL cp_sm_fm_multiply(matrix_s, matrix_x, fmout, k, error=error)
    walltime2 = m_walltime()
    oldtime = walltime2-walltime1
    CALL cp_fm_frobenius_norm(fmout, checksum, error)
    WRITE(*,*)'Sx checksum',checksum*checksum,'time=',oldtime
    CALL dbcsr_from_fm(matrix_b, fmout, 0.0d0, dist2d,&
         row_blk_size=matrix_a%row_blk_size,&
         col_blk_size=col_blk_size,&
         error=error)
    DEALLOCATE(col_blk_size)
    CALL distribution_2d_release(dist2d, error)
    CALL dbcsr_destroy(matrix_a,error)
    CALL cp_fm_release(fmout, error)
    CALL dbcsr_destroy (matrix_b, error)

    CALL timestop(handle)


  END SUBROUTINE dbcsr_check


END MODULE cp_dbcsr_operations

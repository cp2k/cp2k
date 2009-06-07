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
  USE array_types,                     ONLY: array_data
  USE cp_fm_basic_linalg,              ONLY: cp_fm_frobenius_norm
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE dbcsr_methods,                   ONLY: dbcsr_distribution_mp,&
                                             dbcsr_distribution_release,&
                                             dbcsr_mp_group
  USE dbcsr_operations,                ONLY: dbcsr_multiply_sets,&
                                             dbcsr_multiply_sets_sqrt
  USE dbcsr_types,                     ONLY: dbcsr_array_type,&
                                             dbcsr_distribution_obj,&
                                             dbcsr_imagedistribution_type,&
                                             dbcsr_type
  USE dbcsr_util,                      ONLY: &
       dbcsr_checksum, dbcsr_create, dbcsr_create_dist2d_r_unrot, &
       dbcsr_create_image_dist, dbcsr_create_image_dist_left, dbcsr_destroy, &
       dbcsr_destroy_array, dbcsr_destroy_image_dist, dbcsr_from_fm, &
       dbcsr_from_sm, dbcsr_init, dbcsr_new_transposed, dbcsr_normalize, &
       dbcsr_set_debug, dbcsr_transpose_distribution, dbcsr_verify_matrix, &
       dbcsr_work_create
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_retain,&
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
!> \param[in] dist          distribution for the product matrix.
!> \param[in] hint            code path selection
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_sm_sm_multiply(matrix_a, matrix_b, matrix_c, distribution, &
       hint, error)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix_a, matrix_b
    TYPE(dbcsr_type), INTENT(OUT)            :: matrix_c
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: distribution
    INTEGER, INTENT(IN)                      :: hint
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_sm_sm_multiply', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: flop, timing_handle
    REAL(KIND=dp)                            :: checksum, flops
    TYPE(dbcsr_array_type)                   :: ms_left, ms_right
    TYPE(dbcsr_imagedistribution_type)       :: rdist_left, rdist_right
    TYPE(dbcsr_type)                         :: matrix_d, matrix_e

!TYPE(dbcsr_distribution_obj)             :: dist_left, dist_right

    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODYY
    CALL timeset(routineN, timing_handle)
    !> pre The matrices must be valid and they must have compatible shapes.
    ! Assign new dist2d
    ! Create & initialize the final & working matrix.
    IF (dbg) WRITE(*,*)routineP//'Creating 2-d distribution...'
    CALL mp_sync(dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))
    CALL dbcsr_init (matrix_c, error)
    CALL dbcsr_create(matrix_c, "C", distribution, 'N',&
         matrix_a%row_blk_size, matrix_b%col_blk_size,&
         0, 0, error=error)
    IF (dbg) WRITE(*,*)routineP//'Updating new matrix...'
    CALL mp_sync(dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))
    !IF (dbg) WRITE(*,*)routineP//' C initial index is', matrix_c%index
    IF (dbg) WRITE(*,*)routineP//'Creating work matrix...'
    CALL mp_sync(dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))
    CALL dbcsr_work_create(matrix_c,&
         sizedata_guess=INT(1.5*SIZE(matrix_a%data)), error=error)
    IF (dbg) WRITE(*,*)routineP//'Multiplying...'
    CALL mp_sync(dbcsr_mp_group(dbcsr_distribution_mp(matrix_a%dist)))
    !
    DBG 'Creating image distribution'
    CALL dbcsr_create_image_dist (rdist_right, matrix_b%dist)
    CALL dbcsr_create_image_dist_left(rdist_left, matrix_a%dist, matrix_b%dist)
    !
    DBG 'Normalizing A'
    CALL dbcsr_normalize(matrix_a, ms_left, rdist_left, error)
    !ms_left%image_dist = rdist_left
    matrix_d = ms_left%mats(1)
    DBG 'master name', ms_left%mats(1)%name
    DBG "Verifying A->D"
    CALL dbcsr_verify_matrix(matrix_d, error)
    flops = dbcsr_checksum(matrix_d, error)
    DBG "DBCSR norm. S checksum", flops
    DBG 'Normalizing B'
    CALL dbcsr_normalize(matrix_b, ms_right, rdist_right, error)
    matrix_e = ms_right%mats(1)
    DBG "Verifying B->E"
    CALL dbcsr_verify_matrix(matrix_e, error)
    !CALL dbcsr_print(matrix_e, error)
    flops = dbcsr_checksum(matrix_e, error)
    DBG "DBCSR norm. x checksum", flops
    DBG 'Multiplying normalized matrices'
    !
    IF (hint .EQ. 1) THEN
       CALL dbcsr_multiply_sets_sqrt(ms_left, ms_right, matrix_c, error)
    ELSE
       CALL dbcsr_multiply_sets(ms_left, ms_right, matrix_c, error)
    ENDIF
    !
    IF (dbg) THEN
       CALL dbcsr_verify_matrix(matrix_c, error)
       checksum = dbcsr_checksum(matrix_c, error)
       WRITE(*,*)routineP//' BCSR Sx Matrix checksum is', checksum
    ENDIF
    flop = 0
    CALL mp_sum(flop, dbcsr_mp_group (dbcsr_distribution_mp (matrix_a%dist)))
    WRITE(*,*)'FLOP=',flop
    DBG "Result matrix"
    !
    !CALL dbcsr_distribution_release (dist_left)
    !CALL dbcsr_distribution_release (dist_right)
    CALL dbcsr_destroy_image_dist (rdist_left)
    CALL dbcsr_destroy_image_dist (rdist_right)
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
  SUBROUTINE dbcsr_check(matrix_s, matrix_x, hint, error)

    TYPE(real_matrix_type), POINTER          :: matrix_s
    TYPE(cp_fm_type), POINTER                :: matrix_x
    INTEGER, INTENT(IN)                      :: hint
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, k
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size, row_blk_size
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: checksum, newtime, oldtime, &
                                                walltime1, walltime2
    TYPE(cp_fm_type), POINTER                :: fmout
    TYPE(dbcsr_distribution_obj)             :: dist_result
    TYPE(dbcsr_type)                         :: matrix_a, matrix_b, &
                                                matrix_bt, matrix_c
    TYPE(distribution_2d_type), POINTER      :: dist2d, dist2d_left

    CALL timeset(routineN,handle)

    failure=.FALSE.

    !Do the new multiplication
    dist2d_left => matrix_s%distribution_2d
    CALL dbcsr_init (matrix_a, error)
    CALL dbcsr_from_sm(matrix_a, matrix_s, error)
    CALL dbcsr_verify_matrix(matrix_a, error)
    WRITE(*,*)'S checksum=',&
         checksum_matrix(matrix_s, matrix_s%distribution_2d%blacs_env%para_env)
    checksum = 0.0_dp
    checksum = dbcsr_checksum(matrix_a, error)
    WRITE(*,*)'BCSR S checksum',checksum
    CALL cp_fm_frobenius_norm(matrix_x, checksum, error)
    CALL cp_fm_get_info(matrix_x, ncol_global=k, error=error)
    
    CALL dbcsr_create_dist2d_r_unrot(dist2d, dist2d_left, k,&
         col_blk_size, error=error)
    !CALL dbcsr_create_dist2d_r(dist2d, dist2d_left, k,&
    !     col_blk_size, error=error)

    CALL dbcsr_init (matrix_b, error)
    CALL dbcsr_from_fm(matrix_b, matrix_x, 0.0_dp, dist2d,&
         row_blk_size=array_data(matrix_a%row_blk_size),&
         col_blk_size=col_blk_size,&
         error=error)
    CALL dbcsr_verify_matrix(matrix_b, error)

    WRITE(*,*)'x checksum',checksum*checksum
    checksum = dbcsr_checksum(matrix_b, error)
    WRITE(*,*)'DBCSR x checksum',checksum
    !
    ! Transpose of the row (transposed column) matrix
    !CALL dbcsr_create_dist2d_transposed (dist2d_result, dist2d, error)
    CALL dbcsr_transpose_distribution (dist_result, matrix_b%dist)
    !CALL dbcsr_create_dist_r_rot(dist_result, matrix_x%dist, k,&
    !     col_blk_size_array, error=error)
    !CALL array_release (col_blk_size_array)
    !CALL dbcsr_distribution_release (dist_right) ! Not used now

   !
    CALL dbcsr_new_transposed (matrix_bt, matrix_b, error=error)
    CALL dbcsr_verify_matrix (matrix_bt, error)
    !CALL dbcsr_destroy (matrix_bt, error)
    !
    !
    walltime1 = m_walltime()
    CALL dbcsr_sm_sm_multiply(matrix_a, matrix_b, matrix_c,&
         matrix_b%dist, hint, error)
    !CALL dbcsr_sm_sm_multiply(matrix_bt, matrix_a, matrix_c,&
    !     dist_result, hint, error)
    CALL dbcsr_distribution_release (dist_result)
    walltime2 = m_walltime()
    newtime = walltime2-walltime1
    CALL dbcsr_destroy (matrix_bt, error)
    !
    CALL dbcsr_verify_matrix(matrix_c, error)
    !
    checksum = dbcsr_checksum(matrix_c, error)
    WRITE(*,*)'BCSR Sx checksum',checksum,'time=',newtime
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
    WRITE(*,*)'Sx checksum',checksum*checksum,'time=',oldtime
    CALL dbcsr_from_fm(matrix_b, fmout, 0.0_dp, dist2d,&
         row_blk_size=array_data(matrix_a%row_blk_size),&
         col_blk_size=col_blk_size,&
         error=error)
    DEALLOCATE(col_blk_size)
    CALL distribution_2d_release (dist2d, error)
    CALL dbcsr_destroy (matrix_a,error)
    CALL dbcsr_destroy (matrix_b, error)
    CALL cp_fm_release (fmout, error)

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

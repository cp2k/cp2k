!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!
MODULE optbas_opt_utils
  USE admm_types,                      ONLY: admm_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_trace
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_types,                     ONLY: cp_fm_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: evaluate_fval

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optbas_opt_utils'

CONTAINS

  SUBROUTINE evaluate_fval(mos,mos_aux_fit,Q,Snew,admm_env,fval,S_cond_number,error)
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(cp_dbcsr_type), POINTER             :: Q, Snew
    TYPE(admm_type), POINTER                 :: admm_env
    REAL(KIND=dp)                            :: fval, S_cond_number
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'evaluate_fval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, nao_aux_fit, nao_orb, &
                                                nmo, nspins
    REAL(KIND=dp)                            :: trace
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)

    CALL copy_dbcsr_to_fm(Q,admm_env%Q,error)
    fval=0.0_dp
    DO ispin=1,nspins
      nmo = admm_env%nmo(ispin)
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
      CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)

      CALL cp_fm_gemm('N','N',nao_aux_fit,nmo,nao_orb,-2.0_dp,admm_env%Q,mo_coeff,&
                      0.0_dp,admm_env%work_aux_nmo(ispin)%matrix,error)
      CALL cp_fm_trace(mo_coeff_aux_fit,admm_env%work_aux_nmo(ispin)%matrix,trace,error)
      fval=fval+trace+2.0_dp*nmo
    END DO

    ALLOCATE(eigenvalues(nao_aux_fit))
    CALL copy_dbcsr_to_fm(Snew,admm_env%work_aux_aux,error)
    CALL cp_fm_syevd(admm_env%work_aux_aux,admm_env%work_aux_aux2,eigenvalues,error=error)
    S_cond_number=MAXVAL(ABS(eigenvalues))/MAX(MINVAL(ABS(eigenvalues)),EPSILON(0.0_dp))
    DEALLOCATE(eigenvalues)

  END SUBROUTINE evaluate_fval

END MODULE optbas_opt_utils

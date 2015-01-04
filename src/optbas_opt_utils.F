!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!
MODULE optbas_opt_utils
  USE admm_types,                      ONLY: admm_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_fm_basic_linalg,              ONLY: cp_fm_trace,&
                                             cp_fm_upper_to_full
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE kinds,                           ONLY: dp
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: evaluate_fval, evaluate_energy

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optbas_opt_utils'

CONTAINS

   
! *****************************************************************************
!> \brief ...
!> \param mos ...
!> \param matrix_ks ...
!> \param S_inv_orb ...
!> \param Q ...
!> \param tmp1 ...
!> \param energy ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE evaluate_energy(mos,matrix_ks,S_inv_orb,Q,tmp1,energy,error)
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_fm_type), POINTER                :: S_inv_orb, Q, tmp1
    REAL(KIND=dp)                            :: energy
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'evaluate_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, naux, nmo, norb
    REAL(KIND=dp)                            :: tmp_energy
    TYPE(cp_fm_type), POINTER                :: mo_coeff, QS_inv, tmp, tmp2, &
                                                work, work_orb

    CALL cp_fm_create(QS_inv,matrix_struct=Q%matrix_struct,error=error)
    CALL cp_fm_create(tmp,matrix_struct=Q%matrix_struct,error=error) 
    CALL cp_fm_create(tmp2,matrix_struct=tmp1%matrix_struct,error=error)
    CALL cp_fm_create(work,matrix_struct=S_inv_orb%matrix_struct,error=error)
    CALL cp_fm_create(work_orb,matrix_struct=S_inv_orb%matrix_struct,error=error)
    CALL cp_fm_get_info(Q,nrow_global=naux,ncol_global=norb,error=error)
    CALL cp_gemm('N','N',naux,norb,norb,1.0_dp,Q,S_inv_orb,0.0_dp,QS_inv,error)
    energy=0.0_dp 
    DO ispin=1,SIZE(matrix_ks)
       CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix,work,error=error)
       CALL cp_fm_upper_to_full(work,work_orb,error=error)

       CALL get_mo_set(mos(ispin)%mo_set,nmo=nmo,mo_coeff=mo_coeff)
       CALL cp_gemm('N','N',naux,norb,norb,1.0_dp,QS_inv,work,0.0_dp,tmp,error)
       CALL cp_gemm('N','T',naux,naux,norb,1.0_dp,tmp,QS_inv,0.0_dp,tmp1,error)
       CALL cp_gemm('N','T',naux,naux,nmo,1.0_dp,mo_coeff,mo_coeff,0.0_dp,tmp2,error)
       CALL cp_fm_trace(tmp1,tmp2,tmp_energy,error)
       energy=energy+tmp_energy*(3.0_dp-REAL(SIZE(matrix_ks),dp))
       
    END DO
    
    CALL cp_fm_release(work_orb,error)
    CALL cp_fm_release(QS_inv,error)
    CALL cp_fm_release(tmp,error)
    CALL cp_fm_release(tmp2,error)
    CALL cp_fm_release(work,error)
    CALL cp_fm_release(work_orb,error)

  END SUBROUTINE evaluate_energy

! *****************************************************************************
!> \brief ...
!> \param mos ...
!> \param mos_aux_fit ...
!> \param Q ...
!> \param Snew ...
!> \param admm_env ...
!> \param fval ...
!> \param S_cond_number ...
!> \param error ...
! *****************************************************************************
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

      CALL cp_gemm('N','N',nao_aux_fit,nmo,nao_orb,-2.0_dp,admm_env%Q,mo_coeff,&
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

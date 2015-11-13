!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief PAO pre-optimization based on McWeeny purification
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_preopt_mcweeny
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_frobenius_norm, &
       cp_dbcsr_init, cp_dbcsr_multiply, cp_dbcsr_p_type, cp_dbcsr_release, &
       cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_type
  USE cp_log_handling,                 ONLY: cp_logger_get_default_unit_nr
  USE dm_ls_scf_qs,                    ONLY: matrix_qs_to_ls
  USE dm_ls_scf_types,                 ONLY: ls_mstruct_type,&
                                             ls_scf_env_type
  USE kinds,                           ONLY: dp
  USE pao_types,                       ONLY: mcweeny_history_p_type,&
                                             mcweeny_history_type,&
                                             pao_env_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_preopt_mcweeny'

  PUBLIC :: pao_preopt_dm_mcweeny, pao_calc_outer_grad_mcweeny, pao_preopt_finalize_mcweeny

CONTAINS

! *****************************************************************************
!> \brief Calculate density matrix using McWeeny purification
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_preopt_dm_mcweeny(pao, qs_env, ls_scf_env)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_preopt_dm_mcweeny', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s, rho_ao
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)
    CALL get_qs_env(qs_env,&
                    matrix_ks=matrix_ks,&
                    matrix_s=matrix_s,&
                    rho=rho)

    IF(.NOT.ASSOCIATED(pao%mcweeny_history))&
       ALLOCATE(pao%mcweeny_history(ls_scf_env%nspins))

    DO ispin=1, ls_scf_env%nspins
       CALL matrix_qs_to_ls(ls_scf_env%matrix_ks(ispin),matrix_ks(ispin)%matrix,&
                            ls_scf_env%ls_mstruct,covariant=.TRUE.)

       CALL qs_rho_get(rho, rho_ao=rho_ao)
       CALL matrix_qs_to_ls(ls_scf_env%matrix_p(ispin),rho_ao(ispin)%matrix,&
                         ls_scf_env%ls_mstruct,covariant=.FALSE.)

       CALL matrix_qs_to_ls(ls_scf_env%matrix_s, matrix_s(1)%matrix,&
                            ls_scf_env%ls_mstruct, covariant=.TRUE.)

       IF(ls_scf_env%nspins==1) CALL cp_dbcsr_scale(ls_scf_env%matrix_p(ispin), 0.5_dp)

       CALL clear_mcweeny_history(pao%mcweeny_history(ispin))
       CALL pao_purify_mcweeny(ls_scf_env%matrix_p(ispin),&
                           ls_scf_env%matrix_s,&
                           pao%mcweeny_history(ispin),&
                           threshold=1e-8_dp, max_steps=100)

       IF(ls_scf_env%nspins==1) CALL cp_dbcsr_scale(ls_scf_env%matrix_p(ispin), 2.0_dp)
    ENDDO
    CALL timestop(handle)
  END SUBROUTINE pao_preopt_dm_mcweeny


! *****************************************************************************
!> \brief Finalize McWeeny data-structures
!> \param pao ...
! *****************************************************************************
  SUBROUTINE pao_preopt_finalize_mcweeny(pao)
    TYPE(pao_env_type), POINTER              :: pao

    INTEGER                                  :: ispin

    DO ispin=1, SIZE(pao%mcweeny_history)
       CALL clear_mcweeny_history(pao%mcweeny_history(ispin))
    ENDDO

    DEALLOCATE(pao%mcweeny_history)

  END SUBROUTINE pao_preopt_finalize_mcweeny


! *****************************************************************************
!> \brief Helper routine to free history
!> \param history ...
! *****************************************************************************
  SUBROUTINE clear_mcweeny_history(history)
    TYPE(mcweeny_history_p_type)             :: history

    TYPE(mcweeny_history_type), POINTER      :: history_curr, history_next

    NULLIFY(history_next, history_curr)
    history_curr => history%p
    NULLIFY(history%p)

    DO WHILE (ASSOCIATED(history_curr))
       CALL cp_dbcsr_release(history_curr%m)
       history_next => history_curr%next
       DEALLOCATE(history_curr)
       history_curr => history_next
       NULLIFY(history_next)
    END DO
  END SUBROUTINE clear_mcweeny_history


! *****************************************************************************
!> \brief Apply McWeeny purification to auxilary density matrix
!> \param matrix_p ...
!> \param matrix_s ...
!> \param history ...
!> \param threshold ...
!> \param max_steps ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE pao_purify_mcweeny(matrix_p, matrix_s, history, threshold, max_steps)
    TYPE(cp_dbcsr_type)                      :: matrix_p, matrix_s
    TYPE(mcweeny_history_p_type)             :: history
    REAL(KIND=dp)                            :: threshold
    INTEGER                                  :: max_steps

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pao_purify_mcweeny', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istep, unit_nr
    REAL(KIND=dp)                            :: frob_norm
    TYPE(cp_dbcsr_type)                      :: matrix_ps, matrix_psp, &
                                                matrix_test
    TYPE(mcweeny_history_type), POINTER      :: new_hist_entry

    CALL timeset(routineN,handle)
    NULLIFY(new_hist_entry)

    unit_nr = cp_logger_get_default_unit_nr()

    IF(ASSOCIATED(history%p)) CPABORT("purify_mcweeny: history already associated")

    CALL cp_dbcsr_init(matrix_PS)
    CALL cp_dbcsr_create(matrix_PS,template=matrix_p,matrix_type="N")
    CALL cp_dbcsr_init(matrix_PSP)
    CALL cp_dbcsr_create(matrix_PSP,template=matrix_p,matrix_type="N")
    CALL cp_dbcsr_init(matrix_test)

    DO istep=1, max_steps
       ! allocate new element in linked list
        ALLOCATE(new_hist_entry)
        new_hist_entry%next => history%p
        history%p => new_hist_entry
        history%p%count = istep
        NULLIFY(new_hist_entry)
        CALL cp_dbcsr_init(history%p%m)
        CALL cp_dbcsr_create(history%p%m, template=matrix_p ,matrix_type="N")
        CALL cp_dbcsr_copy(history%p%m, matrix_p, name="P from McWeeny")

        ! calc PS and PSP
        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_p, matrix_s, 0.0_dp, matrix_ps)
        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_ps, matrix_p, 0.0_dp, matrix_psp)

        !test convergence
        CALL cp_dbcsr_copy(matrix_test, matrix_psp)
        CALL cp_dbcsr_add(matrix_test, matrix_p, 1.0_dp,-1.0_dp)
        frob_norm=cp_dbcsr_frobenius_norm(matrix_test)
        IF(unit_nr>0) WRITE(unit_nr,'(t3,a,i5,a,f16.8)')"McWeeny-Step",istep,&
                            ": Deviation of idempotency", frob_norm
        IF(frob_norm < threshold .AND. istep>1) EXIT

        ! build next P matrix
        CALL cp_dbcsr_copy(matrix_p, matrix_PSP, name="P from McWeeny")
        CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_PS, matrix_PSP, 3.0_dp, matrix_p)
    END DO

    ! clean up
    CALL cp_dbcsr_release(matrix_PS)
    CALL cp_dbcsr_release(matrix_PSP)
    CALL cp_dbcsr_release(matrix_test)
    CALL timestop(handle)
  END SUBROUTINE pao_purify_mcweeny


! *****************************************************************************
!> \brief Calculate partial derivative dE/dU of McWeeny
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param matrix_M ...
! *****************************************************************************
   SUBROUTINE pao_calc_outer_grad_mcweeny(qs_env,ls_scf_env,matrix_M)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env
    TYPE(cp_dbcsr_type)                      :: matrix_M

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_outer_grad_mcweeny', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nspin
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s, rho_ao
    TYPE(cp_dbcsr_type) :: matrix_HB, matrix_HBP, matrix_k, matrix_l, &
      matrix_NHBP, matrix_NPAK, matrix_NSBL, matrix_PA, matrix_PAK, &
      matrix_SB, matrix_SBL
    TYPE(ls_mstruct_type), POINTER           :: ls_mstruct
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)

    ls_mstruct => ls_scf_env%ls_mstruct
    pao => ls_scf_env%pao_env
    nspin = ls_scf_env%nspins

    ! matrix_K is called M^1 in pao_math4.pdf
    CALL cp_dbcsr_init(matrix_k)
    CALL cp_dbcsr_copy(matrix_k, ls_scf_env%matrix_ks(1), name="K")

    ! matrix_L is called M^4 in pao_math4.pdf
    CALL cp_dbcsr_init(matrix_L)
    CALL cp_dbcsr_create(matrix_L, template=ls_scf_env%matrix_s, matrix_type="N")
    CALL cp_dbcsr_set(matrix_L, 0.0_dp)

    CALL pao_revert_purify_mcweeny(matrix_k=matrix_k,&
                               matrix_s=ls_scf_env%matrix_s,&
                               matrix_s_deriv=matrix_l,&
                               history=pao%mcweeny_history(1))

    CALL get_qs_env(qs_env,&
                    rho=rho,&
                    matrix_ks=matrix_ks,&
                    matrix_s=matrix_s)
    CALL qs_rho_get(rho, rho_ao=rho_ao)

    CALL cp_dbcsr_init(matrix_M)
    CALL cp_dbcsr_create(matrix_M, template=matrix_s(1)%matrix, matrix_type="N")

    !------------------------------------------------------------------
    ! calculate M^2 = N^T P A K^T Y^T
    CALL cp_dbcsr_init(matrix_PA)
    CALL cp_dbcsr_create(matrix_PA, template=ls_mstruct%matrix_A, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, rho_ao(1)%matrix, ls_mstruct%matrix_A,&
                           0.0_dp, matrix_PA)

    CALL cp_dbcsr_init(matrix_PAK)
    CALL cp_dbcsr_create(matrix_PAK, template=ls_mstruct%matrix_A, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_PA, matrix_K,&
                           0.0_dp, matrix_PAK)

    CALL cp_dbcsr_init(matrix_NPAK)
    CALL cp_dbcsr_create(matrix_NPAK, template=ls_mstruct%matrix_A, matrix_type="N")
    CALL cp_dbcsr_multiply("T", "N", 1.0_dp, pao%matrix_N_inv, matrix_PAK,&
                           0.0_dp, matrix_NPAK)

    CALL cp_dbcsr_multiply("N", "T", 2.0_dp, matrix_NPAK, pao%matrix_Y,&
                           1.0_dp, matrix_M)

    !------------------------------------------------------------------
    ! calculate M^3 = N^T H  B \tilde P Y^T
    CALL cp_dbcsr_init(matrix_HB)
    CALL cp_dbcsr_create(matrix_HB, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_ks(1)%matrix, ls_mstruct%matrix_B,&
                           0.0_dp, matrix_HB)

    CALL cp_dbcsr_init(matrix_HBP)
    CALL cp_dbcsr_create(matrix_HBP, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_HB, ls_scf_env%matrix_p(1),&
                           0.0_dp, matrix_HBP)

    CALL cp_dbcsr_init(matrix_NHBP)
    CALL cp_dbcsr_create(matrix_NHBP, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("T", "N", 1.0_dp, pao%matrix_N, matrix_HBP,&
                           0.0_dp, matrix_NHBP)

    CALL cp_dbcsr_multiply("N", "T", 2.0_dp, matrix_NHBP, pao%matrix_Y,&
                           1.0_dp, matrix_M)

    !------------------------------------------------------------------
    ! calculate M^5 = N^T S B L^T Y^T
    CALL cp_dbcsr_init(matrix_SB)
    CALL cp_dbcsr_create(matrix_SB, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s(1)%matrix, ls_mstruct%matrix_B,&
                           0.0_dp, matrix_SB)

    CALL cp_dbcsr_init(matrix_SBL)
    CALL cp_dbcsr_create(matrix_SBL, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_SB, matrix_L,&
                           0.0_dp, matrix_SBL)

    CALL cp_dbcsr_init(matrix_NSBL)
    CALL cp_dbcsr_create(matrix_NSBL, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("T", "N", 1.0_dp, pao%matrix_N, matrix_SBL,&
                           0.0_dp, matrix_NSBL)

    CALL cp_dbcsr_multiply("N", "T", 2.0_dp, matrix_NSBL, pao%matrix_Y,&
                           1.0_dp, matrix_M)

    CALL cp_dbcsr_release(matrix_K)
    CALL cp_dbcsr_release(matrix_PA)
    CALL cp_dbcsr_release(matrix_PAK)
    CALL cp_dbcsr_release(matrix_NPAK)
    CALL cp_dbcsr_release(matrix_HB)
    CALL cp_dbcsr_release(matrix_HBP)
    CALL cp_dbcsr_release(matrix_NHBP)
    CALL cp_dbcsr_release(matrix_L)
    CALL cp_dbcsr_release(matrix_SB)
    CALL cp_dbcsr_release(matrix_SBL)
    CALL cp_dbcsr_release(matrix_NSBL)
    CALL timestop(handle)
  END SUBROUTINE pao_calc_outer_grad_mcweeny


! *****************************************************************************
!> \brief Prepare auxilary KS-matrix for merge using reverse McWeeny
!> \param matrix_k ...
!> \param matrix_s ...
!> \param matrix_s_deriv ...
!> \param history ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE pao_revert_purify_mcweeny(matrix_k, matrix_s, matrix_s_deriv, history)
    TYPE(cp_dbcsr_type)                      :: matrix_k, matrix_s, &
                                                matrix_s_deriv
    TYPE(mcweeny_history_p_type)             :: history

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pao_revert_purify_mcweeny', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, unit_nr
    TYPE(mcweeny_history_type), POINTER      :: history_curr, history_next

    CALL timeset(routineN,handle)
    unit_nr = cp_logger_get_default_unit_nr()
    NULLIFY(history_next, history_curr)

    IF(.NOT.ASSOCIATED(history%p))&
       CPABORT("mcweeny history not associated")

    history_curr => history%p
    NULLIFY(history%p)

    CALL cp_dbcsr_set(matrix_s_deriv, 0.0_dp)

    ! reverse McWeeny iteration
    DO WHILE (ASSOCIATED(history_curr))
       IF(unit_nr>0) WRITE(unit_nr,'(t3,a,i5)') "Reverse McWeeny-Step ", history_curr%count
       CALL reverse_mcweeny_add_S_deriv(matrix_k=matrix_k,&
                                        matrix_s=matrix_s,&
                                        matrix_p=history_curr%m,&
                                        matrix_s_deriv=matrix_s_deriv)
       CALL reverse_mcweeny_step(matrix_k=matrix_k,&
                                 matrix_s=matrix_s,&
                                 matrix_p=history_curr%m)
       CALL cp_dbcsr_release(history_curr%m)
       history_next => history_curr%next
       DEALLOCATE(history_curr)
       history_curr => history_next
       NULLIFY(history_next)
    END DO

    !TODO: where is this factor coming from? Symmetry or nspin==1 ???
    CALL cp_dbcsr_scale(matrix_s_deriv, 2.0_dp)

    CALL timestop(handle)

  END SUBROUTINE pao_revert_purify_mcweeny


! *****************************************************************************
!> \brief Helper routine to calculate partial derivative dE/dS
!> \param matrix_k ...
!> \param matrix_s ...
!> \param matrix_p ...
!> \param matrix_s_deriv ...
! *****************************************************************************
  SUBROUTINE reverse_mcweeny_add_S_deriv(matrix_k, matrix_s, matrix_p, matrix_s_deriv)
    TYPE(cp_dbcsr_type)                      :: matrix_k, matrix_s, matrix_p, &
                                                matrix_s_deriv

    TYPE(cp_dbcsr_type)                      :: matrix_PK, matrix_PKP, &
                                                matrix_PS

    CALL cp_dbcsr_init(matrix_PK)
    CALL cp_dbcsr_create(matrix_PK, template=matrix_p, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_P, matrix_K,&
                           0.0_dp, matrix_PK)

    CALL cp_dbcsr_init(matrix_PKP)
    CALL cp_dbcsr_create(matrix_PKP, template=matrix_p, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_PK, matrix_P,&
                           0.0_dp, matrix_PKP)

    CALL cp_dbcsr_init(matrix_PS)
    CALL cp_dbcsr_create(matrix_PS, template=matrix_p, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_P, matrix_S,&
                           0.0_dp, matrix_PS)

    CALL cp_dbcsr_add(matrix_s_deriv, matrix_PKP, 1.0_dp, 3.0_dp)

    CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_PS, matrix_PKP,&
                           1.0_dp, matrix_s_deriv)

    CALL cp_dbcsr_multiply("N", "T", -2.0_dp, matrix_PKP, matrix_PS,&
                           1.0_dp, matrix_s_deriv)

    CALL cp_dbcsr_release(matrix_PK)
    CALL cp_dbcsr_release(matrix_PKP)
    CALL cp_dbcsr_release(matrix_PS)
  END SUBROUTINE reverse_mcweeny_add_S_deriv


! *****************************************************************************
!> \brief Multiply matrix_k with partial derivative of McWeeny by reversing it.
!> \param matrix_k ...
!> \param matrix_s ...
!> \param matrix_p ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE reverse_mcweeny_step(matrix_k, matrix_s, matrix_p)
    TYPE(cp_dbcsr_type)                      :: matrix_k, matrix_s, matrix_p

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reverse_mcweeny_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_dbcsr_type)                      :: matrix_ps, matrix_sp, &
                                                matrix_sum, matrix_tmp

    CALL timeset(routineN,handle)
    CALL cp_dbcsr_init(matrix_ps)
    CALL cp_dbcsr_create(matrix_ps,template=matrix_p,matrix_type="N")
    CALL cp_dbcsr_init(matrix_sp)
    CALL cp_dbcsr_create(matrix_sp,template=matrix_p,matrix_type="N")
    CALL cp_dbcsr_init(matrix_tmp)
    CALL cp_dbcsr_create(matrix_tmp,template=matrix_p,matrix_type="N")
    CALL cp_dbcsr_init(matrix_sum)
    CALL cp_dbcsr_create(matrix_sum,template=matrix_p,matrix_type="N")

    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_p, matrix_s,&
             0.0_dp, matrix_ps)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s, matrix_p,&
             0.0_dp, matrix_sp)

    !TODO: can we exploid more symmetry?
    CALL cp_dbcsr_multiply("N", "N", 3.0_dp, matrix_k, matrix_ps,&
             0.0_dp, matrix_sum)
    CALL cp_dbcsr_multiply("N", "N", 3.0_dp, matrix_sp, matrix_k,&
             1.0_dp, matrix_sum)

    !matrix_tmp = KPS
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_k, matrix_ps,&
             0.0_dp, matrix_tmp)
    CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_tmp, matrix_ps,&
            1.0_dp, matrix_sum)
    CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_sp, matrix_tmp,&
            1.0_dp, matrix_sum)

    !matrix_tmp = SPK
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_sp, matrix_k,&
             0.0_dp, matrix_tmp)
    CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_sp, matrix_tmp,&
            1.0_dp, matrix_sum)

    ! overwrite matrix_k
    CALL cp_dbcsr_copy(matrix_k, matrix_sum, name="K from reverse McWeeny")

    ! clean up
    CALL cp_dbcsr_release(matrix_sum)
    CALL cp_dbcsr_release(matrix_tmp)
    CALL cp_dbcsr_release(matrix_ps)
    CALL cp_dbcsr_release(matrix_sp)
    CALL timestop(handle)
  END SUBROUTINE reverse_mcweeny_step

END MODULE pao_preopt_mcweeny

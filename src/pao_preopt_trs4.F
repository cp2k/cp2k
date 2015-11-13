!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief PAO pre-optimization based on TRS4 purification
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_preopt_trs4
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_create,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_multiply,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_release,&
                                             cp_dbcsr_type
  USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                             cp_logger_type
  USE dm_ls_scf_methods,               ONLY: density_matrix_trs4,&
                                             ls_scf_init_matrix_S
  USE dm_ls_scf_qs,                    ONLY: matrix_qs_to_ls
  USE dm_ls_scf_types,                 ONLY: ls_mstruct_type,&
                                             ls_scf_env_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE pao_types,                       ONLY: pao_env_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_preopt_trs4'

  PUBLIC :: pao_preopt_dm_trs4, pao_calc_outer_grad_lnv, pao_preopt_finalize_trs4

CONTAINS

! *****************************************************************************
!> \brief Finalize TRS4 data-structures
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_preopt_finalize_trs4(ls_scf_env)
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    CALL release_matrix_s_relatives(ls_scf_env)
  END SUBROUTINE pao_preopt_finalize_trs4

! *****************************************************************************
!> \brief Helper routines to release matrices derived from S
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE release_matrix_s_relatives(ls_scf_env)
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    CALL cp_dbcsr_release(ls_scf_env%matrix_s_inv)
    CALL cp_dbcsr_release(ls_scf_env%matrix_s_sqrt)
    CALL cp_dbcsr_release(ls_scf_env%matrix_s_sqrt_inv)
    ! we keep ls_scf_env%matrix_s
  END SUBROUTINE release_matrix_s_relatives

! *****************************************************************************
!> \brief Calculate density matrix using TRS4 purification
!> \param qs_env ...
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE pao_preopt_dm_trs4(qs_env, ls_scf_env)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_preopt_dm_trs4', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: project_name
    INTEGER                                  :: handle, ispin, &
                                                nelectron_spin_real, nspin
    REAL(KIND=dp)                            :: homo_spin, lumo_spin, mu_spin
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    logger => cp_get_default_logger()
    project_name = logger%iter_info%project_name
    nspin = ls_scf_env%nspins

    CALL get_qs_env(qs_env,&
                    matrix_ks=matrix_ks,&
                    matrix_s=matrix_s)

    ! rebuild S, S_inv, S_sqrt, and S_sqrt_inv in pao basis
    CALL release_matrix_s_relatives(ls_scf_env)
    CALL ls_scf_init_matrix_s(matrix_s(1)%matrix, ls_scf_env)

    DO ispin=1, nspin
       CALL matrix_qs_to_ls(ls_scf_env%matrix_ks(ispin),matrix_ks(ispin)%matrix,&
                            ls_scf_env%ls_mstruct,covariant=.TRUE.)

       nelectron_spin_real = ls_scf_env%nelectron_spin(ispin)
       CALL density_matrix_trs4(ls_scf_env%matrix_p(ispin), ls_scf_env%matrix_ks(ispin),&
                                ls_scf_env%matrix_s_sqrt_inv,&
                                nelectron_spin_real, 1.0E-06_dp, homo_spin, lumo_spin, mu_spin, &
                                dynamic_threshold=.FALSE., eps_lanczos=1.0E-03_dp, max_iter_lanczos=128)
    ENDDO
    CALL timestop(handle)
  END SUBROUTINE pao_preopt_dm_trs4


! *****************************************************************************
!> \brief Helper routine, calculates partial derivative dE/dU
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param matrix_M ...
! *****************************************************************************
   SUBROUTINE pao_calc_outer_grad_lnv(qs_env,ls_scf_env,matrix_M)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type), TARGET            :: ls_scf_env
    TYPE(cp_dbcsr_type)                      :: matrix_M

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_outer_grad_lnv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nspin
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s, rho_ao
    TYPE(cp_dbcsr_type) :: matrix_HB, matrix_HPS, matrix_M1, matrix_M2, &
      matrix_M3, matrix_NHB, matrix_NHBM2, matrix_NPA, matrix_NPAM1, &
      matrix_NSB, matrix_NSBM3, matrix_PA, matrix_PH, matrix_PHP, matrix_PS, &
      matrix_PSP, matrix_SB
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ls_mstruct_type), POINTER           :: ls_mstruct
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)

    ls_mstruct => ls_scf_env%ls_mstruct
    pao => ls_scf_env%pao_env

    CALL get_qs_env(qs_env,&
                    rho=rho,&
                    matrix_ks=matrix_ks,&
                    matrix_s=matrix_s,&
                    dft_control=dft_control)
    CALL qs_rho_get(rho, rho_ao=rho_ao)
    nspin = dft_control%nspins

    IF(nspin /= 1) CPABORT("open shell not yet implemented")
    !TODO: handle openshell case properly

    ! notation according to pao_math_lnv.pdf

    CALL cp_dbcsr_init(matrix_M)
    CALL cp_dbcsr_create(matrix_M, template=matrix_s(1)%matrix, matrix_type="N")

    !---------------------------------------------------------------------------
    ! calculate need products in pao basis
    CALL cp_dbcsr_init(matrix_PH)
    CALL cp_dbcsr_create(matrix_PH, template=ls_scf_env%matrix_s, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_p(1), ls_scf_env%matrix_ks(1),&
                          0.0_dp, matrix_PH)

    CALL cp_dbcsr_init(matrix_PHP)
    CALL cp_dbcsr_create(matrix_PHP, template=ls_scf_env%matrix_s, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_PH, ls_scf_env%matrix_p(1),&
                           0.0_dp, matrix_PHP)

    CALL cp_dbcsr_init(matrix_PS)
    CALL cp_dbcsr_create(matrix_PS, template=ls_scf_env%matrix_s, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_p(1), ls_scf_env%matrix_s, &
                          0.0_dp, matrix_PS)

    CALL cp_dbcsr_init(matrix_HPS)
    CALL cp_dbcsr_create(matrix_HPS, template=ls_scf_env%matrix_s, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_ks(1), matrix_PS, &
                          0.0_dp, matrix_HPS)

    CALL cp_dbcsr_init(matrix_PSP)
    CALL cp_dbcsr_create(matrix_PSP, template=ls_scf_env%matrix_s, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_PS, ls_scf_env%matrix_p(1),&
                           0.0_dp, matrix_PSP)

    !---------------------------------------------------------------------------
    ! M1 = dE_lnv / dP_pao
    CALL cp_dbcsr_init(matrix_M1)
    CALL cp_dbcsr_create(matrix_M1, template=ls_scf_env%matrix_s, matrix_type="N")

    CALL cp_dbcsr_multiply("N", "N", 3.0_dp, ls_scf_env%matrix_ks(1), matrix_PS,&
                           1.0_dp, matrix_M1)

    CALL cp_dbcsr_multiply("T", "N", 3.0_dp, matrix_PS, ls_scf_env%matrix_ks(1),&
                           1.0_dp, matrix_M1)

    CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_HPS, matrix_PS,&
                           1.0_dp, matrix_M1)

    CALL cp_dbcsr_multiply("T", "N", -2.0_dp, matrix_PS, matrix_HPS, &
                           1.0_dp, matrix_M1)

    CALL cp_dbcsr_multiply("T", "T", -2.0_dp, matrix_PS, matrix_HPS,&
                           1.0_dp, matrix_M1)


    !---------------------------------------------------------------------------
    ! M2 = dE_lnv / dH
    CALL cp_dbcsr_init(matrix_M2)
    CALL cp_dbcsr_create(matrix_M2, template=ls_scf_env%matrix_s, matrix_type="N")

    CALL cp_dbcsr_add(matrix_M2, matrix_PSP, 1.0_dp, 3.0_dp)

    CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_PS, matrix_PSP,&
                           1.0_dp, matrix_M2)

    !---------------------------------------------------------------------------
    ! M3 = dE_lnv / dS
    CALL cp_dbcsr_init(matrix_M3)
    CALL cp_dbcsr_create(matrix_M3, template=ls_scf_env%matrix_s, matrix_type="N")

    CALL cp_dbcsr_add(matrix_M3, matrix_PHP, 1.0_dp, 3.0_dp)

    CALL cp_dbcsr_multiply("N", "T", -2.0_dp, matrix_PHP, matrix_PS,&
                           1.0_dp, matrix_M3)

    CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_PS, matrix_PHP,&
                           1.0_dp, matrix_M3)


    !---------------------------------------------------------------------------
    ! combine M1 with matrices from primary basis
    CALL cp_dbcsr_init(matrix_PA)
    CALL cp_dbcsr_create(matrix_PA, template=ls_mstruct%matrix_A, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, rho_ao(1)%matrix, ls_mstruct%matrix_A ,&
                           0.0_dp, matrix_PA)

    CALL cp_dbcsr_init(matrix_NPA)
    CALL cp_dbcsr_create(matrix_NPA, template=ls_mstruct%matrix_A, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N_inv, matrix_PA,&
                           0.0_dp, matrix_NPA)

    CALL cp_dbcsr_init(matrix_NPAM1)
    CALL cp_dbcsr_create(matrix_NPAM1, template=ls_mstruct%matrix_A, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_NPA, matrix_M1,&
                           0.0_dp, matrix_NPAM1)

    CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_NPAM1, pao%matrix_Y,&
                           1.0_dp, matrix_M)


    !---------------------------------------------------------------------------
    ! combine M2 with matrices from primary basis
    CALL cp_dbcsr_init(matrix_HB)
    CALL cp_dbcsr_create(matrix_HB, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_ks(1)%matrix, ls_mstruct%matrix_B,&
                           0.0_dp, matrix_HB)

    CALL cp_dbcsr_init(matrix_NHB)
    CALL cp_dbcsr_create(matrix_NHB, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N, matrix_HB,&
                           0.0_dp, matrix_NHB)

    CALL cp_dbcsr_init(matrix_NHBM2)
    CALL cp_dbcsr_create(matrix_NHBM2, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_NHB, matrix_M2,&
                           0.0_dp, matrix_NHBM2)

    CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_NHBM2, pao%matrix_Y,&
                           1.0_dp, matrix_M)


    !---------------------------------------------------------------------------
    ! combine M3 with matrices from primary basis
    CALL cp_dbcsr_init(matrix_SB)
    CALL cp_dbcsr_create(matrix_SB, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s(1)%matrix, ls_mstruct%matrix_B,&
                           0.0_dp, matrix_SB)

    CALL cp_dbcsr_init(matrix_NSB)
    CALL cp_dbcsr_create(matrix_NSB, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N, matrix_SB,&
                           0.0_dp, matrix_NSB)

    CALL cp_dbcsr_init(matrix_NSBM3)
    CALL cp_dbcsr_create(matrix_NSBM3, template=ls_mstruct%matrix_B, matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_NSB, matrix_M3,&
                           0.0_dp, matrix_NSBM3)

    CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_NSBM3, pao%matrix_Y,&
                           1.0_dp, matrix_M)

    !---------------------------------------------------------------------------
    ! cleanup
    CALL cp_dbcsr_release(matrix_PH)
    CALL cp_dbcsr_release(matrix_PHP)
    CALL cp_dbcsr_release(matrix_PS)
    CALL cp_dbcsr_release(matrix_HPS)
    CALL cp_dbcsr_release(matrix_PSP)
    CALL cp_dbcsr_release(matrix_M1)
    CALL cp_dbcsr_release(matrix_M2)
    CALL cp_dbcsr_release(matrix_M3)
    CALL cp_dbcsr_release(matrix_PA)
    CALL cp_dbcsr_release(matrix_NPA)
    CALL cp_dbcsr_release(matrix_NPAM1)
    CALL cp_dbcsr_release(matrix_HB)
    CALL cp_dbcsr_release(matrix_NHB)
    CALL cp_dbcsr_release(matrix_NHBM2)
    CALL cp_dbcsr_release(matrix_SB)
    CALL cp_dbcsr_release(matrix_NSB)
    CALL cp_dbcsr_release(matrix_NSBM3)

    CALL timestop(handle)
  END SUBROUTINE pao_calc_outer_grad_lnv

END MODULE pao_preopt_trs4

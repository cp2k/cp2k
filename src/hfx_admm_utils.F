!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for hfx and admm methods
!>
!>
!> \par History
!>     refactoring 03-2011 [MI]
!> \author MI
! *****************************************************************************
MODULE hfx_admm_utils

  USE admm_methods,                    ONLY: admm_calculate_density_matrix,&
                                             admm_fit_mo_coeffs
  USE admm_types,                      ONLY: admm_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_fm_types,                     ONLY: cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type                   
  USE input_constants,                 ONLY: use_aux_fit_basis_set 
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: mo_set_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: admm_density_matrix

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_admm_utils'

CONTAINS

! *****************************************************************************
!> \brief 
!>
!> \par History
!>     refactoring 03-2011 [MI]
! *****************************************************************************
  SUBROUTINE  admm_density_matrix(qs_env,ks_env,dft_control,input,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_density_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ispin, n_rep_hf, nspins
    LOGICAL :: failure, do_hfx

    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p, matrix_p_aux_fit
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(section_vals_type), POINTER         :: hfx_sections

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(admm_env, hfx_sections,  matrix_p, matrix_p_aux_fit, mos, mos_aux_fit, para_env)    

    hfx_sections => section_vals_get_subs_vals(input,"DFT%XC%HF",error=error)
    CALL section_vals_get(hfx_sections,explicit=do_hfx,error=error)
    !! ** ADMM can only be used with HFX
    IF ( .NOT. do_hfx ) THEN
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                      "Wavefunction fitting requested without Hartree-Fock."//&
                       CPSourceFileRef,&
                       only_ionode=.TRUE.)
    END IF

    ! ** Method only implemented for GPW
    IF( dft_control%qs_control%gapw ) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
                               message="ADMM only implemented for GPW", &
                               error=error, error_level=cp_failure_level)
    END IF

    nspins=dft_control%nspins
    CALL get_qs_env(qs_env,mos_aux_fit=mos_aux_fit,&
                    mos=mos,admm_env=admm_env,&
                    para_env=para_env,error=error)

    CALL section_vals_get(hfx_sections,n_repetition=n_rep_hf,error=error)
    IF ( n_rep_hf > 1 ) THEN
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                      "ADMM can handle only one HF section."//&
                       CPSourceFileRef,&
                       only_ionode=.TRUE.)
    END IF

    DO ispin=1,nspins!fm->dbcsr
       IF(mos(ispin)%mo_set%use_mo_coeff_b) THEN!fm->dbcsr
          CALL copy_dbcsr_to_fm(mos(ispin)%mo_set%mo_coeff_b,mos(ispin)%mo_set%mo_coeff,error=error)!fm->dbcsr            
       ENDIF!fm->dbcsr
    ENDDO!fm->dbcsr

    CALL admm_fit_mo_coeffs(qs_env, admm_env, dft_control%admm_control, para_env, &
                            qs_env%matrix_s_aux_fit, qs_env%matrix_s_aux_fit_vs_orb,&
                            mos, mos_aux_fit, ks_env%s_mstruct_changed, error=error)

    CALL set_qs_env(qs_env,admm_env=admm_env,error=error)

    DO ispin=1,nspins
       matrix_p_aux_fit =>qs_env%rho_aux_fit%rho_ao
       matrix_p => qs_env%rho%rho_ao
       CALL admm_calculate_density_matrix(admm_env, mos_aux_fit(ispin)%mo_set, matrix_p(ispin)%matrix, &
                                          matrix_p_aux_fit(ispin)%matrix,ispin,nspins,error=error)
       CALL calculate_rho_elec(matrix_p=qs_env%rho_aux_fit%rho_ao(ispin)%matrix,&
            rho=qs_env%rho_aux_fit%rho_r(ispin),&
            rho_gspace=qs_env%rho_aux_fit%rho_g(ispin),&
            total_rho=qs_env%rho_aux_fit%tot_rho_r(ispin),&
            qs_env=qs_env,soft_valid=.FALSE.,&
            basis_set_id=use_aux_fit_basis_set, error=error)
    END DO

    qs_env%rho_aux_fit%rho_r_valid=.TRUE.
    qs_env%rho_aux_fit%rho_g_valid=.TRUE.

    CALL timestop(handle)

  END SUBROUTINE admm_density_matrix

END MODULE hfx_admm_utils



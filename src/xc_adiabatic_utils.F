!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief
!>
!>
!> \par History
!>     refactoring 03-2011 [MI]
!> \author MI
! *****************************************************************************
MODULE xc_adiabatic_utils

  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE hfx_communication,               ONLY: scale_and_add_fock_to_ks_matrix
  USE hfx_derivatives,                 ONLY: derivatives_four_center
  USE input_constants,                 ONLY: do_adiabatic_hybrid_mcy3,&
                                             do_adiabatic_model_pade
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE pw_types,                        ONLY: pw_p_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_vxc,                          ONLY: qs_vxc_create
  USE qs_vxc_atom,                     ONLY: calculate_vxc_atom
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_adiabatic_methods,            ONLY: rescale_MCY3_pade
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: rescale_xc_potential

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_adiabatic_utils'

CONTAINS

! *****************************************************************************
!> \brief
!>
!> \param qs_env ...
!> \param ks_matrix ...
!> \param rho ...
!> \param energy ...
!> \param v_rspace_new ...
!> \param v_tau_rspace ...
!> \param hf_energy ...
!> \param just_energy ...
!> \param calculate_forces ...
!> \param use_virial ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE  rescale_xc_potential(qs_env,ks_matrix,rho,energy,v_rspace_new,v_tau_rspace,&
              hf_energy, just_energy,calculate_forces,use_virial,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: v_rspace_new, v_tau_rspace
    REAL(dp), DIMENSION(:)                   :: hf_energy
    LOGICAL, INTENT(in)                      :: just_energy, &
                                                calculate_forces, use_virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rescale_xc_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: adiabatic_functional, &
                                                adiabatic_model, handle, &
                                                n_rep_hf
    LOGICAL                                  :: do_adiabatic_rescaling, &
                                                do_hfx, failure, gapw, gapw_xc
    REAL(dp) :: adiabatic_lambda, adiabatic_omega, scale_dDFA, scale_ddW0, &
      scale_dEx1, scale_dEx2, total_energy_xc
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(section_vals_type), POINTER         :: adiabatic_rescaling_section, &
                                                hfx_sections, input, &
                                                xc_section

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(para_env, dft_control, adiabatic_rescaling_section, hfx_sections, input, xc_section)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         para_env=para_env,&
         input=input,error=error)

    adiabatic_rescaling_section => section_vals_get_subs_vals(input,"DFT%XC%ADIABATIC_RESCALING",error=error)
    CALL section_vals_get(adiabatic_rescaling_section,explicit=do_adiabatic_rescaling,error=error)
    hfx_sections => section_vals_get_subs_vals(input,"DFT%XC%HF",error=error)
    CALL section_vals_get(hfx_sections,explicit=do_hfx,error=error)
    CALL section_vals_get(hfx_sections,n_repetition=n_rep_hf,error=error)

    gapw=dft_control%qs_control%gapw
    gapw_xc=dft_control%qs_control%gapw_xc

    CALL section_vals_val_get(adiabatic_rescaling_section, "FUNCTIONAL_TYPE",&
         i_val=adiabatic_functional,error=error)
    CALL section_vals_val_get(adiabatic_rescaling_section, "FUNCTIONAL_MODEL",&
         i_val=adiabatic_model,error=error)
    CALL section_vals_val_get(adiabatic_rescaling_section, "LAMBDA",&
         r_val=adiabatic_lambda,error=error)
    CALL section_vals_val_get(adiabatic_rescaling_section, "OMEGA",&
         r_val=adiabatic_omega,error=error)
    SELECT CASE(adiabatic_functional)
    CASE (do_adiabatic_hybrid_mcy3)
       SELECT CASE(adiabatic_model)
       CASE(do_adiabatic_model_pade)
          CALL cp_assert( n_rep_hf == 2 , cp_failure_level,cp_assertion_failed,routineP,&
               " For this kind of adiababatic hybrid functional 2 HF sections have to be provided. "//&
               " Please check your input file.",&
               error,failure)
          CALL rescale_MCY3_pade(qs_env, hf_energy, energy, adiabatic_lambda, &
               adiabatic_omega, scale_dEx1, scale_ddW0, scale_dDFA,&
               scale_dEx2, total_energy_xc, error)

          !! Scale and add Fock matrix to KS matrix
          IF(do_hfx) THEN
             CALL scale_and_add_fock_to_ks_matrix(para_env, qs_env, ks_matrix, 1 ,&
                                                  scale_dEx1, error)
             CALL scale_and_add_fock_to_ks_matrix(para_env, qs_env, ks_matrix, 2 ,&
                                                  scale_dEx2, error)
          END IF
          IF( calculate_forces ) THEN
             CPPrecondition(.NOT.use_virial,cp_failure_level,routineP,error,failure)
             !! we also have to scale the forces!!!!
             CALL derivatives_four_center(qs_env, rho%rho_ao, hfx_sections, para_env, 1, use_virial, &
                                          adiabatic_rescale_factor=scale_dEx1, error=error)
             CALL derivatives_four_center(qs_env, rho%rho_ao, hfx_sections, para_env, 2, use_virial, &
                                          adiabatic_rescale_factor=scale_dEx2, error=error)
          END IF

          ! Calculate vxc and rescale it
          xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
          IF (gapw_xc) THEN
            CALL qs_vxc_create(qs_env=qs_env, rho_struct=qs_env%rho_xc, xc_section=xc_section, &
                 vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=energy%exc, &
                 just_energy=just_energy ,adiabatic_rescale_factor=scale_dDFA, error=error)
          ELSE
            CALL qs_vxc_create(qs_env=qs_env, rho_struct=qs_env%rho, xc_section=xc_section, &
                 vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=energy%exc, &
                 just_energy=just_energy, adiabatic_rescale_factor=scale_dDFA, error=error)
          END IF

          ! Calculate vxc and rescale it
          IF (gapw .OR. gapw_xc) THEN
             CALL calculate_vxc_atom(qs_env,just_energy,adiabatic_rescale_factor=scale_dDFA,&
                  error=error)
          END IF
          !! Hack for the total energy expression
          energy%ex = 0.0_dp
          energy%exc1 = 0.0_dp
          energy%exc = total_energy_xc

       END SELECT
    END SELECT
    CALL timestop(handle)

  END SUBROUTINE rescale_xc_potential

END MODULE xc_adiabatic_utils



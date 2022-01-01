!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Initialize a qs_env for kpoint calculations starting from a gamma point qs_env
!> \par History
!>      11.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE qs_gamma2kp
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_p_type
   USE input_constants,                 ONLY: atomic_guess,&
                                              xc_none
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_type
   USE pw_types,                        ONLY: pw_p_type
   USE qs_energy_init,                  ONLY: qs_energies_init
   USE qs_environment,                  ONLY: qs_init
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_create,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_initialization,           ONLY: qs_scf_env_init_basic
   USE qs_scf_types,                    ONLY: qs_scf_env_type,&
                                              scf_env_release
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_gamma2kp'

   PUBLIC :: create_kp_from_gamma

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param qs_env_kp ...
!> \param with_xc_terms ...
! **************************************************************************************************
   SUBROUTINE create_kp_from_gamma(qs_env, qs_env_kp, with_xc_terms)
      TYPE(qs_environment_type), POINTER                 :: qs_env, qs_env_kp
      LOGICAL, OPTIONAL                                  :: with_xc_terms

      INTEGER                                            :: ispin, xc_func
      LOGICAL                                            :: without_xc_terms
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp, rho_ao_kp_gamma
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g_gamma, rho_g_kp, rho_r_gamma, &
                                                            rho_r_kp
      TYPE(qs_rho_type), POINTER                         :: rho, rho_gamma
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: force_env_section, subsys_section, &
                                                            xc_fun_section

      !--------------------------------------------------------------------------------------------!

      CPASSERT(.NOT. ASSOCIATED(qs_env_kp))

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      input=force_env_section, &
                      cp_subsys=cp_subsys)

      NULLIFY (subsys_section)

      NULLIFY (kpoint)
      CALL kpoint_create(kpoint)
      kpoint%kp_scheme = "GAMMA"
      kpoint%symmetry = .FALSE.
      kpoint%verbose = .FALSE.
      kpoint%full_grid = .TRUE.
      kpoint%eps_geo = 1.0e-6_dp
      kpoint%use_real_wfn = .TRUE.
      kpoint%parallel_group_size = 0

      without_xc_terms = .FALSE.
      IF (PRESENT(with_xc_terms)) without_xc_terms = .NOT. with_xc_terms

      IF (without_xc_terms) THEN
         xc_fun_section => section_vals_get_subs_vals(force_env_section, "DFT%XC%XC_FUNCTIONAL")
         CALL section_vals_val_get(xc_fun_section, "_SECTION_PARAMETERS_", i_val=xc_func)
         CALL section_vals_val_set(xc_fun_section, "_SECTION_PARAMETERS_", i_val=xc_none)
      END IF

      CALL qs_env_create(qs_env_kp)
      CALL qs_init(qs_env_kp, para_env, cp_subsys=cp_subsys, kpoint_env=kpoint, &
                   force_env_section=force_env_section, subsys_section=subsys_section, &
                   use_motion_section=.FALSE.)

      CALL get_qs_env(qs_env_kp, scf_control=scf_control)
      scf_control%density_guess = atomic_guess
      CALL qs_energies_init(qs_env_kp, calc_forces=.FALSE.)

      NULLIFY (scf_env)
      CALL qs_scf_env_init_basic(qs_env_kp, scf_env)

      CALL set_qs_env(qs_env_kp, scf_env=scf_env)
      CALL scf_env_release(scf_env)

      ! copy density matrix, n(r) and n(G) from Gamma-only to kpoint calculation
      CALL get_qs_env(qs_env, rho=rho_gamma)
      CALL qs_rho_get(rho_gamma, rho_ao_kp=rho_ao_kp_gamma, rho_r=rho_r_gamma, rho_g=rho_g_gamma)
      CALL get_qs_env(qs_env_kp, rho=rho)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp, rho_r=rho_r_kp, rho_g=rho_g_kp)

      DO ispin = 1, SIZE(rho_r_gamma)
         rho_r_kp(ispin)%pw%cr3d(:, :, :) = rho_r_gamma(ispin)%pw%cr3d(:, :, :)
         rho_g_kp(ispin)%pw%cc(:) = rho_g_gamma(ispin)%pw%cc(:)
         CALL dbcsr_add(matrix_a=rho_ao_kp(ispin, 1)%matrix, alpha_scalar=0.0_dp, &
                        matrix_b=rho_ao_kp_gamma(ispin, 1)%matrix, beta_scalar=1.0_dp)
      END DO

      CALL qs_ks_update_qs_env(qs_env_kp, print_active=.FALSE.)

      IF (without_xc_terms) THEN
         ! set back the functional
         CALL section_vals_val_set(xc_fun_section, "_SECTION_PARAMETERS_", i_val=xc_func)
      END IF

   END SUBROUTINE create_kp_from_gamma

END MODULE qs_gamma2kp

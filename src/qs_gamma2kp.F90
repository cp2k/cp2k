!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
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
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_p_type
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE input_constants,                 ONLY: atomic_guess,&
                                              xc_none
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_methods,                      ONLY: pw_copy
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_energy_init,                  ONLY: qs_energies_init
   USE qs_environment,                  ONLY: qs_init
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_create,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_rho_atom_types,               ONLY: rho_atom_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_scf_initialization,           ONLY: qs_scf_env_init_basic
   USE qs_scf_types,                    ONLY: qs_scf_env_type
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

      INTEGER                                            :: iatom, ispin, xc_func
      LOGICAL                                            :: without_xc_terms
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_gamma
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(rho_atom_type), POINTER                       :: atoms(:), atoms_gamma(:)
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: force_env_section, subsys_section, &
                                                            xc_fun_section

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

      ALLOCATE (qs_env_kp)
      CALL qs_env_create(qs_env_kp)
      CALL qs_init(qs_env_kp, para_env, cp_subsys=cp_subsys, kpoint_env=kpoint, &
                   force_env_section=force_env_section, subsys_section=subsys_section, &
                   use_motion_section=.FALSE., silent=.TRUE.)

      CALL get_qs_env(qs_env_kp, scf_control=scf_control)
      scf_control%density_guess = atomic_guess
      CALL qs_energies_init(qs_env_kp, calc_forces=.FALSE.)

      NULLIFY (scf_env)
      CALL qs_scf_env_init_basic(qs_env_kp, scf_env)

      CALL set_qs_env(qs_env_kp, scf_env=scf_env)

      ! Preserve the converged density fields, not the new environment's atomic guess.
      CALL get_qs_env(qs_env, rho=rho_gamma)
      CALL get_qs_env(qs_env_kp, rho=rho)
      CALL copy_frozen_density(rho_gamma, rho)

      CALL get_qs_env(qs_env, dft_control=dft_control)
      IF (dft_control%qs_control%gapw_xc) THEN
         CALL get_qs_env(qs_env, rho_xc=rho_gamma)
         CALL get_qs_env(qs_env_kp, rho_xc=rho)
         CALL copy_frozen_density(rho_gamma, rho)
      END IF
      IF (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc) THEN
         ! The frozen GAPW potential also depends on these hard/soft onsite densities.
         CALL get_qs_env(qs_env, rho_atom_set=atoms_gamma)
         CALL get_qs_env(qs_env_kp, rho_atom_set=atoms)
         CPASSERT(SIZE(atoms) == SIZE(atoms_gamma))
         DO iatom = 1, SIZE(atoms_gamma)
            IF (.NOT. ASSOCIATED(atoms_gamma(iatom)%cpc_h)) CYCLE
            CPASSERT(ASSOCIATED(atoms(iatom)%cpc_h))
            DO ispin = 1, SIZE(atoms_gamma(iatom)%cpc_h)
               IF (.NOT. ASSOCIATED(atoms_gamma(iatom)%cpc_h(ispin)%r_coef)) CYCLE
               CPASSERT(ASSOCIATED(atoms(iatom)%cpc_h(ispin)%r_coef))
               CPASSERT(ASSOCIATED(atoms(iatom)%cpc_s(ispin)%r_coef))
               atoms(iatom)%cpc_h(ispin)%r_coef = atoms_gamma(iatom)%cpc_h(ispin)%r_coef
               atoms(iatom)%cpc_s(ispin)%r_coef = atoms_gamma(iatom)%cpc_s(ispin)%r_coef
            END DO
         END DO
      END IF

      CALL qs_ks_update_qs_env(qs_env_kp, print_active=.FALSE.)

      IF (without_xc_terms) THEN
         ! set back the functional
         CALL section_vals_val_set(xc_fun_section, "_SECTION_PARAMETERS_", i_val=xc_func)
      END IF

   END SUBROUTINE create_kp_from_gamma

! **************************************************************************************************
!> \brief Copy converged Gamma density fields into the initialized k-point environment.
!> \param source Converged Gamma density
!> \param destination Initialized density with the same real/reciprocal grids
!> \author Thomas D. Kuehne (tkuehne@cp2k.org)
! **************************************************************************************************
   SUBROUTINE copy_frozen_density(source, destination)
      TYPE(qs_rho_type), INTENT(IN)                      :: source
      TYPE(qs_rho_type), INTENT(INOUT)                   :: destination

      INTEGER                                            :: idir, ispin
      LOGICAL                                            :: drho_g_valid, drho_r_valid, rho_g_valid, &
                                                            rho_r_valid, soft_valid, tau_g_valid, &
                                                            tau_r_valid
      REAL(dp), DIMENSION(:), POINTER                    :: total_g_in, total_g_out, total_r_in, &
                                                            total_r_out
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ao_in, ao_out
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g_in, rho_g_out, tau_g_in, tau_g_out
      TYPE(pw_c1d_gs_type), DIMENSION(:, :), POINTER     :: drho_g_in, drho_g_out
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r_in, rho_r_out, tau_r_in, tau_r_out
      TYPE(pw_r3d_rs_type), DIMENSION(:, :), POINTER     :: drho_r_in, drho_r_out

      CALL qs_rho_get(source, rho_ao_kp=ao_in, rho_r=rho_r_in, rho_g=rho_g_in, &
                      tau_r=tau_r_in, tau_g=tau_g_in, drho_r=drho_r_in, drho_g=drho_g_in, &
                      tot_rho_r=total_r_in, tot_rho_g=total_g_in, &
                      rho_r_valid=rho_r_valid, rho_g_valid=rho_g_valid, &
                      tau_r_valid=tau_r_valid, tau_g_valid=tau_g_valid, &
                      drho_r_valid=drho_r_valid, drho_g_valid=drho_g_valid, soft_valid=soft_valid)
      CALL qs_rho_get(destination, rho_ao_kp=ao_out, rho_r=rho_r_out, rho_g=rho_g_out, &
                      tau_r=tau_r_out, tau_g=tau_g_out, drho_r=drho_r_out, drho_g=drho_g_out, &
                      tot_rho_r=total_r_out, tot_rho_g=total_g_out)

      CPASSERT(SIZE(rho_r_in) == SIZE(rho_r_out))
      DO ispin = 1, SIZE(rho_r_in)
         CALL pw_copy(rho_r_in(ispin), rho_r_out(ispin))
         CALL pw_copy(rho_g_in(ispin), rho_g_out(ispin))
         CALL dbcsr_add(matrix_a=ao_out(ispin, 1)%matrix, alpha_scalar=0.0_dp, &
                        matrix_b=ao_in(ispin, 1)%matrix, beta_scalar=1.0_dp)
         ! XC fields need not be allocated when with_xc_terms is false.
         IF (ASSOCIATED(tau_r_out) .AND. tau_r_valid) THEN
            CPASSERT(ASSOCIATED(tau_r_in))
            CALL pw_copy(tau_r_in(ispin), tau_r_out(ispin))
         END IF
         IF (ASSOCIATED(tau_g_out) .AND. tau_g_valid) THEN
            CPASSERT(ASSOCIATED(tau_g_in))
            CALL pw_copy(tau_g_in(ispin), tau_g_out(ispin))
         END IF
         IF (ASSOCIATED(drho_r_out) .AND. drho_r_valid) THEN
            CPASSERT(ASSOCIATED(drho_r_in))
            DO idir = 1, SIZE(drho_r_out, 1)
               CALL pw_copy(drho_r_in(idir, ispin), drho_r_out(idir, ispin))
            END DO
         END IF
         IF (ASSOCIATED(drho_g_out) .AND. drho_g_valid) THEN
            CPASSERT(ASSOCIATED(drho_g_in))
            DO idir = 1, SIZE(drho_g_out, 1)
               CALL pw_copy(drho_g_in(idir, ispin), drho_g_out(idir, ispin))
            END DO
         END IF
      END DO
      IF (ASSOCIATED(total_r_in) .AND. ASSOCIATED(total_r_out)) total_r_out = total_r_in
      IF (ASSOCIATED(total_g_in) .AND. ASSOCIATED(total_g_out)) total_g_out = total_g_in
      CALL qs_rho_set(destination, rho_r_valid=rho_r_valid, rho_g_valid=rho_g_valid, &
                      tau_r_valid=tau_r_valid .AND. ASSOCIATED(tau_r_out), &
                      tau_g_valid=tau_g_valid .AND. ASSOCIATED(tau_g_out), &
                      drho_r_valid=drho_r_valid .AND. ASSOCIATED(drho_r_out), &
                      drho_g_valid=drho_g_valid .AND. ASSOCIATED(drho_g_out), soft_valid=soft_valid)

   END SUBROUTINE copy_frozen_density

END MODULE qs_gamma2kp

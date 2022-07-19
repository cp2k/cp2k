!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate the KS reference potentials
!> \par History
!>       07.2022 created
!> \author JGH
! **************************************************************************************************
MODULE qs_ks_reference
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_create,&
                                              pw_p_type,&
                                              pw_release
   USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                              calculate_ecore_self
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_reference'

   LOGICAL, PARAMETER                   :: debug_stress = .TRUE.
   LOGICAL, PARAMETER                   :: debug_forces = .TRUE.

   PUBLIC :: ks_ref_potential

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief calculate the Kohn-Sham reference potential
!> \param qs_env ...
!> \param vh_rspace ...
!> \param vxc_rspace ...
!> \param vtau_rspace ...
!> \param vadmm_rspace ...
!> \param ehartree ...
!> \param exc ...
!> \param h_stress container for the stress tensor of the Hartree term
!> \par History
!>      10.2019 created [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ks_ref_potential(qs_env, vh_rspace, vxc_rspace, vtau_rspace, vadmm_rspace, ehartree, exc, h_stress)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), POINTER                           :: vh_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rspace, vtau_rspace, vadmm_rspace
      REAL(KIND=dp), INTENT(OUT)                         :: ehartree, exc
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: h_stress

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ks_ref_potential'

      INTEGER                                            :: handle, iab, ispin, nspins
      REAL(dp)                                           :: eadmm, eovrl, eself
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_xc
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_admm_rspace, v_admm_tau_rspace, &
                                                            v_rspace, v_tau_rspace
      TYPE(pw_p_type), POINTER                           :: rho_core
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      ! get all information on the electronic density
      NULLIFY (rho, ks_env)
      CALL get_qs_env(qs_env=qs_env, rho=rho, dft_control=dft_control, &
                      para_env=para_env, ks_env=ks_env, rho_core=rho_core)

      nspins = dft_control%nspins

      NULLIFY (pw_env)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CPASSERT(ASSOCIATED(pw_env))

      NULLIFY (auxbas_pw_pool, poisson_env)
      ! gets the tmp grids
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)

      ! Calculate the Hartree potential
      NULLIFY (v_hartree_gspace%pw, rho_tot_gspace%pw, v_hartree_rspace%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      ! Get the total density in g-space [ions + electrons]
      CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, ehartree, &
                            v_hartree_gspace%pw, h_stress=h_stress, rho_core=rho_core)
      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
      !
      CALL calculate_ecore_self(qs_env, E_self_core=eself)
      CALL calculate_ecore_overlap(qs_env, para_env, PRESENT(h_stress), E_overlap_core=eovrl)
      ehartree = ehartree + eovrl + eself

      ! v_rspace and v_tau_rspace are generated from the auxbas pool
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_primary
      ELSE
         xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      END IF
      NULLIFY (v_rspace, v_tau_rspace)
      CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=xc_section, &
                         vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=exc, just_energy=.FALSE.)

      NULLIFY (v_admm_rspace, v_admm_tau_rspace)
      IF (dft_control%do_admm) THEN
         IF (dft_control%admm_control%aux_exch_func /= do_admm_aux_exch_func_none) THEN
            ! For the virial, we have to save the pv_xc component because it will be reset in qs_vxc_create
            IF (PRESENT(h_stress)) THEN
               CALL get_qs_env(qs_env, virial=virial)
               virial_xc = virial%pv_xc
            END IF
            CALL get_admm_env(admm_env, rho_aux_fit=rho)
            xc_section => admm_env%xc_section_aux
            CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=xc_section, &
                               vxc_rho=v_admm_rspace, vxc_tau=v_admm_tau_rspace, exc=eadmm, just_energy=.FALSE.)
            IF (PRESENT(h_stress)) virial%pv_xc = virial%pv_xc + virial_xc
         END IF
      END IF

      ! allocate potentials
      IF (ASSOCIATED(vh_rspace)) THEN
         CALL pw_release(vh_rspace%pw)
      ELSE
         ALLOCATE (vh_rspace)
         NULLIFY (vh_rspace%pw)
      END IF
      IF (ASSOCIATED(vxc_rspace)) THEN
         DO iab = 1, SIZE(vxc_rspace)
            CALL pw_release(vxc_rspace(iab)%pw)
         END DO
      ELSE
         ALLOCATE (vxc_rspace(nspins))
         DO iab = 1, nspins
            NULLIFY (vxc_rspace(iab)%pw)
         END DO
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
         IF (ASSOCIATED(vtau_rspace)) THEN
            DO iab = 1, SIZE(vtau_rspace)
               CALL pw_release(vtau_rspace(iab)%pw)
            END DO
         ELSE
            ALLOCATE (vtau_rspace(nspins))
            DO iab = 1, nspins
               NULLIFY (vtau_rspace(iab)%pw)
            END DO
         END IF
      ELSE
         NULLIFY (vtau_rspace)
      END IF
      IF (ASSOCIATED(v_admm_rspace)) THEN
         IF (ASSOCIATED(vadmm_rspace)) THEN
            DO iab = 1, SIZE(vadmm_rspace)
               CALL pw_release(vadmm_rspace(iab)%pw)
            END DO
         ELSE
            ALLOCATE (vadmm_rspace(nspins))
            DO iab = 1, nspins
               NULLIFY (vadmm_rspace(iab)%pw)
            END DO
         END IF
      ELSE
         NULLIFY (vadmm_rspace)
      END IF

      pw_grid => v_hartree_rspace%pw%pw_grid
      CALL pw_create(vh_rspace%pw, pw_grid, use_data=REALDATA3D, in_space=REALSPACE)
      DO ispin = 1, nspins
         NULLIFY (vxc_rspace(ispin)%pw)
         CALL pw_create(vxc_rspace(ispin)%pw, pw_grid, &
                        use_data=REALDATA3D, in_space=REALSPACE)
         IF (ASSOCIATED(vtau_rspace)) THEN
            NULLIFY (vtau_rspace(ispin)%pw)
            CALL pw_create(vtau_rspace(ispin)%pw, pw_grid, &
                           use_data=REALDATA3D, in_space=REALSPACE)
         END IF
         IF (ASSOCIATED(vadmm_rspace)) THEN
            NULLIFY (vadmm_rspace(ispin)%pw)
            CALL pw_create(vadmm_rspace(ispin)%pw, pw_grid, &
                           use_data=REALDATA3D, in_space=REALSPACE)
         END IF
      END DO
      !
      CALL pw_transfer(v_hartree_rspace%pw, vh_rspace%pw)
      IF (ASSOCIATED(v_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_transfer(v_rspace(ispin)%pw, vxc_rspace(ispin)%pw)
            CALL pw_scale(vxc_rspace(ispin)%pw, v_rspace(ispin)%pw%pw_grid%dvol)
            IF (ASSOCIATED(v_tau_rspace)) THEN
               CALL pw_transfer(v_tau_rspace(ispin)%pw, vtau_rspace(ispin)%pw)
               CALL pw_scale(vtau_rspace(ispin)%pw, v_tau_rspace(ispin)%pw%pw_grid%dvol)
            END IF
         END DO
      ELSE
         DO ispin = 1, nspins
            CALL pw_zero(vxc_rspace(ispin)%pw)
         END DO
      END IF
      IF (ASSOCIATED(v_admm_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_transfer(v_admm_rspace(ispin)%pw, vadmm_rspace(ispin)%pw)
            CALL pw_scale(vadmm_rspace(ispin)%pw, vadmm_rspace(ispin)%pw%pw_grid%dvol)
         END DO
      END IF

      ! return pw grids
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      IF (ASSOCIATED(v_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace(ispin)%pw)
            IF (ASSOCIATED(v_tau_rspace)) THEN
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw)
            END IF
         END DO
         DEALLOCATE (v_rspace)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) DEALLOCATE (v_tau_rspace)
      IF (ASSOCIATED(v_admm_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_admm_rspace(ispin)%pw)
         END DO
         DEALLOCATE (v_admm_rspace)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ks_ref_potential

! **************************************************************************************************

END MODULE qs_ks_reference

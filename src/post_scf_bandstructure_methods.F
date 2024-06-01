!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE post_scf_bandstructure_methods
   USE gw_main,                         ONLY: gw
   USE input_section_types,             ONLY: section_vals_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: create_and_init_bs_env,&
                                              dos_pdos_ldos
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_scf,                          ONLY: scf
   USE soc_pseudopotential_methods,     ONLY: H_KS_spinor,&
                                              V_SOC_xyz_from_pseudopotential
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'post_scf_bandstructure_methods'

   PUBLIC :: post_scf_bandstructure

CONTAINS

! **************************************************************************************************
!> \brief Perform post-SCF band structure calculations from higher level methods
!> \param qs_env  Quickstep environment
!> \param post_scf_bandstructure_section ...
!> \par History
!>    * 07.2023 created [Jan Wilhelm]
! **************************************************************************************************
   SUBROUTINE post_scf_bandstructure(qs_env, post_scf_bandstructure_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: post_scf_bandstructure_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'post_scf_bandstructure'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! general setup of post SCF bandstructure calculation
      CALL create_and_init_bs_env(qs_env, qs_env%bs_env, post_scf_bandstructure_section)

      ! shifts of eigenvalues/bandstructure due to spin-orbit coupling from pseudopotentials
      IF (qs_env%bs_env%do_soc) THEN
         CALL soc(qs_env, qs_env%bs_env)
      END IF

      ! GW calculation for eigenvalues/bandstructure for molecules and periodic systems
      IF (qs_env%bs_env%do_gw) THEN
         CALL gw(qs_env, qs_env%bs_env, post_scf_bandstructure_section)
      END IF

      ! density of states (DOS), projected DOS, local DOS for DFT, DFT+SOC, G0W0, G0W0+SOC
      CALL dos_pdos_ldos(qs_env, qs_env%bs_env)

      CALL timestop(handle)

   END SUBROUTINE post_scf_bandstructure

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE soc(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'soc'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Compute V^SOC_µν^(α) = ħ/2 < ϕ_µ | sum_ℓ ΔV_ℓ^SO(r,r') L^(α) | ϕ_ν >, α = x, y, z, see
      ! Hartwigsen, Goedecker, Hutter, Eq.(18), (19) (doi.org/10.1103/PhysRevB.58.3641)
      CALL V_SOC_xyz_from_pseudopotential(qs_env, bs_env%mat_V_SOC_xyz)

      ! Spinor KS-matrix H_µν,σσ' = h^SCF_µν*δ_σσ' + sum_α V^SOC_µν^(α)*Pauli-matrix^(α)_σσ', see
      ! Hartwigsen, Goedecker, Hutter, Eq.(18) (doi.org/10.1103/PhysRevB.58.3641)
      CALL H_KS_spinor(bs_env%cfm_ks_spinor_ao_Gamma, bs_env%fm_ks_Gamma(1:2), bs_env%n_spin, &
                       bs_env%mat_V_SOC_xyz(:, 1), bs_env%cfm_s_spinor_Gamma, bs_env%fm_s_Gamma, &
                       bs_env%cfm_SOC_spinor_ao_Gamma)

      CALL timestop(handle)

   END SUBROUTINE soc

END MODULE post_scf_bandstructure_methods

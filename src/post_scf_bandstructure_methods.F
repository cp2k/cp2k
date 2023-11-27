!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE post_scf_bandstructure_methods
   USE gw_methods,                      ONLY: gw
   USE input_section_types,             ONLY: section_vals_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: bandstructure_primitive_cell,&
                                              bandstructure_primitive_cell_spinor,&
                                              create_and_init_bs_env,&
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

      ! GW calculation for eigenvalues/bandstructure
      IF (qs_env%bs_env%do_gw) THEN
         CALL gw(qs_env, qs_env%bs_env, post_scf_bandstructure_section)
      END IF

      ! density of states (DOS) and projected DOS for DFT, DFT+SOC, G0W0, G0W0+SOC, also
      ! quantities from local DOS (LDOS) as local valence band maximum (VBM), local conduction
      ! band minimum (CBM), and local gap are calculated (local: as function of space r)
      CALL dos_pdos_ldos(qs_env, qs_env%bs_env)

      ! band structure of primitive unit cell for DFT, DFT+SOC, G0W0, G0W0+SOC
      IF (qs_env%bs_env%do_bs_primitive_cell) THEN
         CALL bandstructure_primitive_cell_all_methods(qs_env, qs_env%bs_env)
      END IF

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

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE bandstructure_primitive_cell_all_methods(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bandstructure_primitive_cell_all_methods'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL bandstructure_primitive_cell(qs_env, bs_env, &
                                        bs_env%eigenval_prim_cell_scf, &
                                        "bandstructure_SCF.bs", &
                                        bs_env%fm_ks_Gamma(1))
      IF (bs_env%do_gw) THEN
         CALL bandstructure_primitive_cell(qs_env, bs_env, &
                                           bs_env%eigenval_prim_cell_G0W0, &
                                           "bandstructure_G0W0.bs", &
                                           bs_env%fm_h_G0W0_Gamma)
      END IF

      IF (bs_env%do_soc) THEN
         CALL bandstructure_primitive_cell_spinor(qs_env, bs_env, &
                                                  bs_env%eigenval_prim_cell_scf_soc, &
                                                  "bandstructure_SCF_SOC.bs", &
                                                  bs_env%cfm_ks_spinor_ao_Gamma)
      END IF

      CALL timestop(handle)

   END SUBROUTINE bandstructure_primitive_cell_all_methods

END MODULE post_scf_bandstructure_methods

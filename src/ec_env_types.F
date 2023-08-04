!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types needed for a for a Energy Correction
!> \par History
!>       2019.09 created
!> \author JGH
! **************************************************************************************************
MODULE ec_env_types
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type,&
                                              ls_scf_release
   USE hfx_types,                       ONLY: hfx_release,&
                                              hfx_type
   USE input_section_types,             ONLY: section_vals_release,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_types,                        ONLY: pw_release,&
                                              pw_type
   USE qs_dispersion_types,             ONLY: qs_dispersion_release,&
                                              qs_dispersion_type
   USE qs_force_types,                  ONLY: deallocate_qs_force,&
                                              qs_force_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE qs_period_efield_types,          ONLY: efield_berry_release,&
                                              efield_berry_type
   USE task_list_types,                 ONLY: deallocate_task_list,&
                                              task_list_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_env_types'

   PUBLIC :: energy_correction_type, ec_env_release

! *****************************************************************************
!> \brief Contains information on the energy correction functional for KG
!> \par History
!>       03.2014 created
!> \author JGH
! *****************************************************************************
   TYPE energy_correction_type
      CHARACTER(len=20)                                :: ec_name = ""
      INTEGER                                          :: energy_functional = 0
      INTEGER                                          :: ks_solver = 0
      INTEGER                                          :: factorization = 0
      INTEGER                                          :: ec_initial_guess = 0
      REAL(KIND=dp)                                    :: eps_default = 0.0_dp
      LOGICAL                                          :: do_ec_admm = .FALSE.
      LOGICAL                                          :: do_ec_hfx = .FALSE.
      LOGICAL                                          :: should_update = .FALSE.
      LOGICAL                                          :: use_ls_solver = .FALSE.
      LOGICAL                                          :: reuse_hfx = .FALSE.
      LOGICAL                                          :: basis_inconsistent = .FALSE.
      ! debug
      LOGICAL                                          :: debug_forces = .FALSE.
      LOGICAL                                          :: debug_stress = .FALSE.
      ! basis set
      CHARACTER(len=20)                                :: basis = ""
      LOGICAL                                          :: mao = .FALSE.
      ! Skip EC calculation if ground-state didnt converge
      LOGICAL                                          :: do_skip = .FALSE., skip_ec = .FALSE.
      INTEGER                                          :: mao_max_iter = 0
      REAL(KIND=dp)                                    :: mao_eps_grad = 0.0_dp
      REAL(KIND=dp)                                    :: mao_eps1 = 0.0_dp
      INTEGER                                          :: mao_iolevel = 0
      ! energy components
      REAL(KIND=dp)                                    :: etotal = 0.0_dp, old_etotal = 0.0_dp
      REAL(KIND=dp)                                    :: eband = 0.0_dp, ecore = 0.0_dp, exc = 0.0_dp, &
                                                          ehartree = 0.0_dp, vhxc = 0.0_dp
      REAL(KIND=dp)                                    :: edispersion = 0.0_dp, efield_elec = 0.0_dp, &
                                                          efield_nuclear = 0.0_dp, ex = 0.0_dp, exc_aux_fit = 0.0_dp
      ! forces
      TYPE(qs_force_type), DIMENSION(:), POINTER       :: force => Null()
      ! full neighbor lists and corresponding task list
      TYPE(neighbor_list_set_p_type), &
         DIMENSION(:), POINTER                         :: sab_orb => Null(), sac_ppl => Null(), sap_ppnl => Null()
      TYPE(task_list_type), POINTER                    :: task_list => Null()
      ! the XC function to be used for the correction, dispersion info
      TYPE(section_vals_type), POINTER                 :: xc_section => Null()
      TYPE(qs_dispersion_type), POINTER                :: dispersion_env => Null()
      ! matrices in complete basis
      ! KS: Kohn-Sham; H: Core; S: overlap; T: kinetic energy;
      ! P: Harris density, W: Harris energy weighted density
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER     :: matrix_ks => Null(), &
                                                          matrix_h => Null(), &
                                                          matrix_s => Null(), &
                                                          matrix_t => Null(), &
                                                          matrix_p => Null(), &
                                                          matrix_w => Null()
      ! reduce basis
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER        :: mao_coef => Null()
      ! CP equations
      TYPE(qs_p_env_type), POINTER                     :: p_env => Null()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER        :: matrix_hz => Null(), matrix_wz => Null(), &
                                                          matrix_z => Null(), z_admm => Null()
      ! Harris (rhoout), and response density (rhoz) on grid
      TYPE(pw_type), DIMENSION(:), POINTER             :: rhoout_r => Null(), rhoz_r => Null()
      ! potentials from input density
      TYPE(pw_type)                                    :: vh_rspace = pw_type()
      TYPE(pw_type), DIMENSION(:), POINTER             :: vxc_rspace => Null(), vtau_rspace => Null(), vadmm_rspace => Null()
      ! efield
      TYPE(efield_berry_type), POINTER                 :: efield => NULL()
      ! LS matrices and types
      TYPE(ls_scf_env_type), POINTER                   :: ls_env => Null()
      ! Environment for Hartree-Fock exchange
      TYPE(hfx_type), DIMENSION(:, :), POINTER         :: x_data => Null()
      ! ADMM XC environments
      TYPE(section_vals_type), POINTER         :: xc_section_primary => Null(), &
                                                  xc_section_aux => Null()
   END TYPE energy_correction_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param ec_env ...
! **************************************************************************************************
   SUBROUTINE ec_env_release(ec_env)
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ec_env_release'

      INTEGER                                            :: handle, iab

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(ec_env)) THEN
         ! neighbor lists
         CALL release_neighbor_list_sets(ec_env%sab_orb)
         CALL release_neighbor_list_sets(ec_env%sac_ppl)
         CALL release_neighbor_list_sets(ec_env%sap_ppnl)
         ! forces
         IF (ASSOCIATED(ec_env%force)) CALL deallocate_qs_force(ec_env%force)
         ! operator matrices
         IF (ASSOCIATED(ec_env%matrix_ks)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_ks)
         IF (ASSOCIATED(ec_env%matrix_h)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_h)
         IF (ASSOCIATED(ec_env%matrix_s)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_s)
         IF (ASSOCIATED(ec_env%matrix_t)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_t)
         IF (ASSOCIATED(ec_env%matrix_p)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_p)
         IF (ASSOCIATED(ec_env%matrix_w)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_w)
         IF (ASSOCIATED(ec_env%task_list)) THEN
            CALL deallocate_task_list(ec_env%task_list)
         END IF
         ! reduced basis
         IF (ASSOCIATED(ec_env%mao_coef)) CALL dbcsr_deallocate_matrix_set(ec_env%mao_coef)
         ! dispersion environment
         IF (ASSOCIATED(ec_env%dispersion_env)) THEN
            CALL qs_dispersion_release(ec_env%dispersion_env)
         END IF

         IF (ASSOCIATED(ec_env%matrix_hz)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_hz)
         IF (ASSOCIATED(ec_env%matrix_wz)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_wz)
         IF (ASSOCIATED(ec_env%matrix_z)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_z)
         IF (ASSOCIATED(ec_env%z_admm)) CALL dbcsr_deallocate_matrix_set(ec_env%z_admm)
         NULLIFY (ec_env%matrix_z, ec_env%matrix_hz, ec_env%matrix_wz)
         NULLIFY (ec_env%z_admm)

         IF (ASSOCIATED(ec_env%p_env)) THEN
            CALL p_env_release(ec_env%p_env)
            DEALLOCATE (ec_env%p_env)
         END IF
         ! potential
         IF (ASSOCIATED(ec_env%vh_rspace%pw_grid)) THEN
            CALL pw_release(ec_env%vh_rspace)
         END IF
         IF (ASSOCIATED(ec_env%vxc_rspace)) THEN
            DO iab = 1, SIZE(ec_env%vxc_rspace)
               CALL pw_release(ec_env%vxc_rspace(iab))
            END DO
            DEALLOCATE (ec_env%vxc_rspace)
         END IF
         IF (ASSOCIATED(ec_env%vtau_rspace)) THEN
            DO iab = 1, SIZE(ec_env%vtau_rspace)
               CALL pw_release(ec_env%vtau_rspace(iab))
            END DO
            DEALLOCATE (ec_env%vtau_rspace)
         END IF
         IF (ASSOCIATED(ec_env%vadmm_rspace)) THEN
            DO iab = 1, SIZE(ec_env%vadmm_rspace)
               CALL pw_release(ec_env%vadmm_rspace(iab))
            END DO
            DEALLOCATE (ec_env%vadmm_rspace)
         END IF
         CALL efield_berry_release(ec_env%efield)

         IF (ASSOCIATED(ec_env%ls_env)) THEN
            CALL ls_scf_release(ec_env%ls_env)
         END IF

         IF (.NOT. ec_env%reuse_hfx) THEN
            IF (ASSOCIATED(ec_env%x_data)) CALL hfx_release(ec_env%x_data)
         END IF

         IF (ASSOCIATED(ec_env%xc_section_aux)) CALL section_vals_release(ec_env%xc_section_aux)
         IF (ASSOCIATED(ec_env%xc_section_primary)) CALL section_vals_release(ec_env%xc_section_primary)

         DEALLOCATE (ec_env)

      END IF

      NULLIFY (ec_env)

      CALL timestop(handle)

   END SUBROUTINE ec_env_release

END MODULE ec_env_types

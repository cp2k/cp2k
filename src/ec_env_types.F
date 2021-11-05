!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
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
   USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                              cp_fm_release
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type,&
                                              ls_scf_release
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_types,                        ONLY: pw_p_type,&
                                              pw_release
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
      CHARACTER(len=20)                                :: ec_name
      INTEGER                                          :: energy_functional
      INTEGER                                          :: ks_solver
      INTEGER                                          :: factorization
      INTEGER                                          :: ec_initial_guess
      REAL(KIND=dp)                                    :: eps_default
      LOGICAL                                          :: should_update
      LOGICAL                                          :: use_ls_solver
      ! basis set
      CHARACTER(len=20)                                :: basis
      LOGICAL                                          :: mao
      ! Skip EC calculation if ground-state didnt converge
      LOGICAL                                          :: do_skip, skip_ec
      INTEGER                                          :: mao_max_iter
      REAL(KIND=dp)                                    :: mao_eps_grad
      REAL(KIND=dp)                                    :: mao_eps1
      INTEGER                                          :: mao_iolevel
      ! energy components
      REAL(KIND=dp)                                    :: etotal
      REAL(KIND=dp)                                    :: eband, exc, ehartree, vhxc
      REAL(KIND=dp)                                    :: edispersion, efield_nuclear
      ! forces
      TYPE(qs_force_type), DIMENSION(:), POINTER       :: force => Null()
      ! full neighbor lists and corresponding task list
      TYPE(neighbor_list_set_p_type), &
         DIMENSION(:), POINTER                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(task_list_type), POINTER                    :: task_list
      ! the XC function to be used for the correction, dispersion info
      TYPE(section_vals_type), POINTER                 :: xc_section
      TYPE(qs_dispersion_type), POINTER                :: dispersion_env
      ! matrices in complete basis
      ! KS: Kohn-Sham; H: Core; S: overlap; T: kinetic energy;
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER     :: matrix_ks, &
                                                          matrix_h, &
                                                          matrix_s, &
                                                          matrix_t, &
                                                          matrix_p, &
                                                          matrix_w
      ! reduce basis
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER        :: mao_coef
      ! CP equations
      TYPE(qs_p_env_type), POINTER                     :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER        :: cpmos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER        :: matrix_hz, matrix_z, matrix_wz, z_admm
      ! potentials from input density
      TYPE(pw_p_type), POINTER                         :: vh_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER           :: vxc_rspace, vtau_rspace, vadmm_rspace
      ! efield
      TYPE(efield_berry_type), POINTER                 :: efield => NULL()
      ! LS matrices and types
      TYPE(ls_scf_env_type), POINTER                   :: ls_env
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
         ! CP env
         IF (ASSOCIATED(ec_env%cpmos)) THEN
            DO iab = 1, SIZE(ec_env%cpmos)
               CALL cp_fm_release(ec_env%cpmos(iab)%matrix)
            END DO
            DEALLOCATE (ec_env%cpmos)
            NULLIFY (ec_env%cpmos)
         END IF

         IF (ASSOCIATED(ec_env%matrix_z)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_z)
         IF (ASSOCIATED(ec_env%matrix_hz)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_hz)
         IF (ASSOCIATED(ec_env%matrix_wz)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_wz)
         IF (ASSOCIATED(ec_env%z_admm)) CALL dbcsr_deallocate_matrix_set(ec_env%z_admm)
         NULLIFY (ec_env%matrix_z, ec_env%matrix_hz, ec_env%matrix_wz)
         NULLIFY (ec_env%z_admm)

         IF (ASSOCIATED(ec_env%p_env)) THEN
            CALL p_env_release(ec_env%p_env)
         END IF
         ! potential
         IF (ASSOCIATED(ec_env%vh_rspace)) THEN
            CALL pw_release(ec_env%vh_rspace%pw)
            DEALLOCATE (ec_env%vh_rspace)
         END IF
         IF (ASSOCIATED(ec_env%vxc_rspace)) THEN
            DO iab = 1, SIZE(ec_env%vxc_rspace)
               CALL pw_release(ec_env%vxc_rspace(iab)%pw)
            END DO
            DEALLOCATE (ec_env%vxc_rspace)
         END IF
         IF (ASSOCIATED(ec_env%vtau_rspace)) THEN
            DO iab = 1, SIZE(ec_env%vtau_rspace)
               CALL pw_release(ec_env%vtau_rspace(iab)%pw)
            END DO
            DEALLOCATE (ec_env%vtau_rspace)
         END IF
         IF (ASSOCIATED(ec_env%vadmm_rspace)) THEN
            DO iab = 1, SIZE(ec_env%vadmm_rspace)
               CALL pw_release(ec_env%vadmm_rspace(iab)%pw)
            END DO
            DEALLOCATE (ec_env%vadmm_rspace)
         END IF
         CALL efield_berry_release(ec_env%efield)

         IF (ASSOCIATED(ec_env%ls_env)) THEN
            CALL ls_scf_release(ec_env%ls_env)
         ENDIF

         DEALLOCATE (ec_env)

      END IF

      NULLIFY (ec_env)

      CALL timestop(handle)

   END SUBROUTINE ec_env_release

END MODULE ec_env_types

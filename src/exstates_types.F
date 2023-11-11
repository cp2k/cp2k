!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types for excited states potential energies
!> \par History
!>       2020.01 created
!> \author JGH
! **************************************************************************************************
MODULE exstates_types
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_constants,                 ONLY: xc_kernel_method_best
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE pw_types,                        ONLY: pw_type
   USE qs_local_rho_types,              ONLY: local_rho_set_release,&
                                              local_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'exstates_types'

   PUBLIC :: excited_energy_type, exstate_release, exstate_create
   PUBLIC :: exstate_matrix_release, exstate_potential_release
   PUBLIC :: wfn_history_type

! *****************************************************************************

   TYPE wfn_history_type
      INTEGER                                            :: state
      REAL(KIND=dp)                                      :: evalue
      TYPE(cp_fm_type), POINTER, DIMENSION(:)            :: evect => NULL()
      TYPE(cp_fm_type), POINTER, DIMENSION(:)            :: cpmos => NULL()
      REAL(KIND=dp)                                      :: xsval, gsval, gsmin
   END TYPE wfn_history_type

! *****************************************************************************
!> \brief Contains information on the excited states energy
!> \par History
!>       01.2020 created
!> \author JGH
! *****************************************************************************
   TYPE excited_energy_type
      INTEGER                                            :: state
      REAL(KIND=dp)                                      :: evalue
      INTEGER                                            :: xc_kernel_method
      LOGICAL                                            :: debug_forces = .FALSE.
      TYPE(cp_fm_type), POINTER, DIMENSION(:)            :: evect => NULL()
      TYPE(cp_fm_type), POINTER, DIMENSION(:)            :: cpmos => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pe => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hz => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pe_admm => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_px1 => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_px1_admm => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_px1_asymm => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_px1_admm_asymm => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_wx1 => NULL()
      TYPE(pw_type)                                      :: vh_rspace
      TYPE(pw_type), DIMENSION(:), POINTER               :: vxc_rspace => NULL()
      TYPE(pw_type), DIMENSION(:), POINTER               :: vtau_rspace => NULL()
      TYPE(pw_type), DIMENSION(:), POINTER               :: vadmm_rspace => NULL()
      TYPE(local_rho_type), POINTER                      :: local_rho_set => NULL()
      TYPE(local_rho_type), POINTER                      :: local_rho_set_admm => NULL()
      TYPE(wfn_history_type)                             :: wfn_history
   END TYPE excited_energy_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param ex_env ...
! **************************************************************************************************
   SUBROUTINE exstate_release(ex_env)
      TYPE(excited_energy_type), POINTER                 :: ex_env

      IF (ASSOCIATED(ex_env)) THEN

         CALL cp_fm_release(ex_env%evect)
         CALL cp_fm_release(ex_env%cpmos)
         !
         CALL exstate_matrix_release(ex_env)
         !
         CALL exstate_potential_release(ex_env)
         !
         CALL cp_fm_release(ex_env%wfn_history%evect)
         CALL cp_fm_release(ex_env%wfn_history%cpmos)

         DEALLOCATE (ex_env)

      END IF

   END SUBROUTINE exstate_release

! **************************************************************************************************
!> \brief ...
!> \param ex_env ...
! **************************************************************************************************
   SUBROUTINE exstate_matrix_release(ex_env)
      TYPE(excited_energy_type), POINTER                 :: ex_env

      IF (ASSOCIATED(ex_env)) THEN
         IF (ASSOCIATED(ex_env%matrix_pe)) CALL dbcsr_deallocate_matrix_set(ex_env%matrix_pe)
         NULLIFY (ex_env%matrix_pe)
         IF (ASSOCIATED(ex_env%matrix_hz)) CALL dbcsr_deallocate_matrix_set(ex_env%matrix_hz)
         NULLIFY (ex_env%matrix_hz)
         IF (ASSOCIATED(ex_env%matrix_pe_admm)) CALL dbcsr_deallocate_matrix_set(ex_env%matrix_pe_admm)
         NULLIFY (ex_env%matrix_pe_admm)
         IF (ASSOCIATED(ex_env%matrix_px1)) CALL dbcsr_deallocate_matrix_set(ex_env%matrix_px1)
         NULLIFY (ex_env%matrix_px1)
         IF (ASSOCIATED(ex_env%matrix_px1_admm)) CALL dbcsr_deallocate_matrix_set(ex_env%matrix_px1_admm)
         NULLIFY (ex_env%matrix_px1_admm)
         IF (ASSOCIATED(ex_env%matrix_px1_asymm)) CALL dbcsr_deallocate_matrix_set(ex_env%matrix_px1_asymm)
         NULLIFY (ex_env%matrix_px1_asymm)
         IF (ASSOCIATED(ex_env%matrix_px1_admm_asymm)) CALL dbcsr_deallocate_matrix_set(ex_env%matrix_px1_admm_asymm)
         NULLIFY (ex_env%matrix_px1_admm_asymm)
         IF (ASSOCIATED(ex_env%matrix_wx1)) CALL dbcsr_deallocate_matrix_set(ex_env%matrix_wx1)
         NULLIFY (ex_env%matrix_wx1)
      END IF

   END SUBROUTINE exstate_matrix_release

! **************************************************************************************************
!> \brief ...
!> \param ex_env ...
! **************************************************************************************************
   SUBROUTINE exstate_potential_release(ex_env)
      TYPE(excited_energy_type), POINTER                 :: ex_env

      INTEGER                                            :: iab

      IF (ASSOCIATED(ex_env)) THEN
         IF (ASSOCIATED(ex_env%vh_rspace%pw_grid)) THEN
            CALL ex_env%vh_rspace%release()
         END IF
         IF (ASSOCIATED(ex_env%vxc_rspace)) THEN
            DO iab = 1, SIZE(ex_env%vxc_rspace)
               CALL ex_env%vxc_rspace(iab)%release()
            END DO
            DEALLOCATE (ex_env%vxc_rspace)
            NULLIFY (ex_env%vxc_rspace)
         END IF
         IF (ASSOCIATED(ex_env%vtau_rspace)) THEN
            DO iab = 1, SIZE(ex_env%vtau_rspace)
               CALL ex_env%vtau_rspace(iab)%release()
            END DO
            DEALLOCATE (ex_env%vtau_rspace)
            NULLIFY (ex_env%vtau_rspace)
         END IF
         IF (ASSOCIATED(ex_env%vadmm_rspace)) THEN
            DO iab = 1, SIZE(ex_env%vadmm_rspace)
               CALL ex_env%vadmm_rspace(iab)%release()
            END DO
            DEALLOCATE (ex_env%vadmm_rspace)
            NULLIFY (ex_env%vadmm_rspace)
         END IF
         IF (ASSOCIATED(ex_env%local_rho_set)) THEN
            CALL local_rho_set_release(ex_env%local_rho_set)
            NULLIFY (ex_env%local_rho_set)
         END IF
         IF (ASSOCIATED(ex_env%local_rho_set_admm)) THEN
            CALL local_rho_set_release(ex_env%local_rho_set_admm)
            NULLIFY (ex_env%local_rho_set_admm)
         END IF
      END IF

   END SUBROUTINE exstate_potential_release

! **************************************************************************************************
!> \brief Allocates and intitializes exstate_env
!> \param ex_env the object to create
!> \param excited_state ...
!> \param dft_section ...
!> \par History
!>       2020.01 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE exstate_create(ex_env, excited_state, dft_section)
      TYPE(excited_energy_type), POINTER                 :: ex_env
      LOGICAL, INTENT(IN)                                :: excited_state
      TYPE(section_vals_type), POINTER                   :: dft_section

      CPASSERT(.NOT. ASSOCIATED(ex_env))
      ALLOCATE (ex_env)
      ex_env%evalue = 0.0_dp
      NULLIFY (ex_env%evect)
      IF (excited_state) THEN
         CALL section_vals_val_get(dft_section, "EXCITED_STATES%STATE", i_val=ex_env%state)
         CALL section_vals_val_get(dft_section, "EXCITED_STATES%XC_KERNEL_METHOD", &
                                   i_val=ex_env%xc_kernel_method)
         CALL section_vals_val_get(dft_section, "EXCITED_STATES%DEBUG_FORCES", &
                                   l_val=ex_env%debug_forces)
      ELSE
         ex_env%state = 0
         ex_env%xc_kernel_method = xc_kernel_method_best
      END IF
      ex_env%wfn_history%evalue = 0.0_dp
      ex_env%wfn_history%state = ex_env%state
      ex_env%wfn_history%xsval = 1.0_dp
      ex_env%wfn_history%gsval = 1.0_dp
      ex_env%wfn_history%gsmin = 1.0_dp
      NULLIFY (ex_env%wfn_history%cpmos)
      NULLIFY (ex_env%wfn_history%evect)

   END SUBROUTINE exstate_create

END MODULE exstates_types

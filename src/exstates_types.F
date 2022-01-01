!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
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
   USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                              cp_fm_release
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE pw_types,                        ONLY: pw_p_type,&
                                              pw_release
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'exstates_types'

   PUBLIC :: excited_energy_type, exstate_release, exstate_create

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
      TYPE(cp_fm_p_type), POINTER, DIMENSION(:)          :: evect => NULL()
      TYPE(cp_fm_p_type), POINTER, DIMENSION(:)          :: cpmos => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pe => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hz => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pe_admm => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_px1 => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_px1_admm => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_wx1 => NULL()
      TYPE(pw_p_type), POINTER                           :: vh_rspace => NULL()
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rspace => NULL()
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vtau_rspace => NULL()
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vadmm_rspace => NULL()
   END TYPE excited_energy_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param ex_env ...
! **************************************************************************************************
   SUBROUTINE exstate_release(ex_env)
      TYPE(excited_energy_type), POINTER                 :: ex_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'exstate_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iab, is

      IF (ASSOCIATED(ex_env)) THEN

         IF (ASSOCIATED(ex_env%evect)) THEN
            DO is = 1, SIZE(ex_env%evect)
               CALL cp_fm_release(ex_env%evect(is)%matrix)
            END DO
            DEALLOCATE (ex_env%evect)
         END IF
         IF (ASSOCIATED(ex_env%cpmos)) THEN
            DO is = 1, SIZE(ex_env%cpmos)
               CALL cp_fm_release(ex_env%cpmos(is)%matrix)
            END DO
            DEALLOCATE (ex_env%cpmos)
         END IF
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
         IF (ASSOCIATED(ex_env%matrix_wx1)) CALL dbcsr_deallocate_matrix_set(ex_env%matrix_wx1)
         NULLIFY (ex_env%matrix_wx1)
         !
         IF (ASSOCIATED(ex_env%vh_rspace)) THEN
            CALL pw_release(ex_env%vh_rspace%pw)
            DEALLOCATE (ex_env%vh_rspace)
         END IF
         IF (ASSOCIATED(ex_env%vxc_rspace)) THEN
            DO iab = 1, SIZE(ex_env%vxc_rspace)
               CALL pw_release(ex_env%vxc_rspace(iab)%pw)
            END DO
            DEALLOCATE (ex_env%vxc_rspace)
         END IF
         IF (ASSOCIATED(ex_env%vtau_rspace)) THEN
            DO iab = 1, SIZE(ex_env%vtau_rspace)
               CALL pw_release(ex_env%vtau_rspace(iab)%pw)
            END DO
            DEALLOCATE (ex_env%vtau_rspace)
         END IF
         IF (ASSOCIATED(ex_env%vadmm_rspace)) THEN
            DO iab = 1, SIZE(ex_env%vadmm_rspace)
               CALL pw_release(ex_env%vadmm_rspace(iab)%pw)
            END DO
            DEALLOCATE (ex_env%vadmm_rspace)
         END IF

         DEALLOCATE (ex_env)

      END IF

   END SUBROUTINE exstate_release

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

      CHARACTER(len=*), PARAMETER :: routineN = 'exstate_create', routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(ex_env))
      ALLOCATE (ex_env)
      ex_env%evalue = 0.0_dp
      NULLIFY (ex_env%evect)
      IF (excited_state) THEN
         CALL section_vals_val_get(dft_section, "EXCITED_STATES%STATE", i_val=ex_env%state)
         CALL section_vals_val_get(dft_section, "EXCITED_STATES%XC_KERNEL_METHOD", &
                                   i_val=ex_env%xc_kernel_method)
      ELSE
         ex_env%state = -1
      END IF

   END SUBROUTINE exstate_create

END MODULE exstates_types

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
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
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
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
      TYPE(cp_fm_p_type), POINTER, DIMENSION(:)          :: evect
   END TYPE excited_energy_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param ec_env ...
! **************************************************************************************************
   SUBROUTINE exstate_release(ex_env)
      TYPE(excited_energy_type), POINTER                 :: ex_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'exstate_release', routineP = moduleN//':'//routineN

      INTEGER                                            :: is

      IF (ASSOCIATED(ex_env)) THEN

         IF(ASSOCIATED(ex_env%evect)) THEN
            DO is = 1, SIZE(ex_env%evect)
               CALL cp_fm_release(ex_env%evect(is)%matrix)
            END DO
            DEALLOCATE (ex_env%evect)
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
      NULLIFY(ex_env%evect)
      IF (excited_state) THEN
         CALL section_vals_val_get(dft_section, "EXCITED_STATES%STATE", i_val=ex_env%state)
      ELSE
         ex_env%state = -1
      END IF

   END SUBROUTINE exstate_create

END MODULE exstates_types

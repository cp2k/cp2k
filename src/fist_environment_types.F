!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>    gt sept-23-02 added atomic_kind_set to replica_environment_type
!>    to allow use of kind_based neighbor list
!>    CJM rewrite
!> \author CJM SEPT-01-02
! **************************************************************************************************
MODULE fist_environment_types
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                              atomic_kind_list_release,&
                                              atomic_kind_list_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_release,&
                                              cell_retain,&
                                              cell_type
   USE cp_result_types,                 ONLY: cp_result_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_release,&
                                              cp_subsys_set,&
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE ewald_environment_types,         ONLY: ewald_env_release,&
                                              ewald_environment_type
   USE ewald_pw_types,                  ONLY: ewald_pw_release,&
                                              ewald_pw_type
   USE exclusion_types,                 ONLY: exclusion_release,&
                                              exclusion_type
   USE fist_efield_types,               ONLY: fist_efield_type
   USE fist_energy_types,               ONLY: deallocate_fist_energy,&
                                              fist_energy_type
   USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_release,&
                                              fist_nonbond_env_type
   USE input_section_types,             ONLY: section_vals_release,&
                                              section_vals_retain,&
                                              section_vals_type
   USE message_passing,                 ONLY: mp_para_env_release,&
                                              mp_para_env_type
   USE molecule_kind_list_types,        ONLY: molecule_kind_list_create,&
                                              molecule_kind_list_release,&
                                              molecule_kind_list_type
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_list_types,             ONLY: molecule_list_create,&
                                              molecule_list_release,&
                                              molecule_list_type
   USE molecule_types,                  ONLY: molecule_type
   USE multipole_types,                 ONLY: multipole_type
   USE particle_list_types,             ONLY: particle_list_create,&
                                              particle_list_release,&
                                              particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE qmmm_types_low,                  ONLY: qmmm_env_mm_create,&
                                              qmmm_env_mm_release,&
                                              qmmm_env_mm_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

! **************************************************************************************************
!> \par History
!>      11/03
!> \author CJM
! **************************************************************************************************
   TYPE fist_environment_type
      PRIVATE
      LOGICAL                                      :: qmmm = .FALSE.
      LOGICAL                                      :: shell_model = .FALSE., shell_model_ad = .FALSE.
      TYPE(qmmm_env_mm_type), POINTER              :: qmmm_env => NULL()
      TYPE(cell_type), POINTER                     :: cell_ref => NULL()
      TYPE(ewald_environment_type), POINTER        :: ewald_env => NULL()
      TYPE(ewald_pw_type), POINTER                 :: ewald_pw => NULL()
      TYPE(fist_energy_type), POINTER              :: thermo => NULL()
      TYPE(mp_para_env_type), POINTER              :: para_env => NULL()
      TYPE(cp_subsys_type), POINTER                :: subsys => NULL()
      TYPE(fist_nonbond_env_type), POINTER         :: fist_nonbond_env => NULL()
      TYPE(section_vals_type), POINTER             :: input => NULL()
      TYPE(exclusion_type), DIMENSION(:), POINTER  :: exclusions => NULL()
      TYPE(fist_efield_type), POINTER              :: efield => NULL()
   END TYPE fist_environment_type

! *** Public data types ***
   PUBLIC :: fist_environment_type

! *** Public subroutines ***
   PUBLIC :: fist_env_get, &
             fist_env_set, &
             fist_env_create, &
             fist_env_release

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_environment_types'

!***

CONTAINS

! **************************************************************************************************
!> \brief Purpose: Get the FIST environment.
!> \param fist_env the pointer to the fist_env
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param ewald_pw ...
!> \param local_particles ...
!> \param local_molecules ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param cell ...
!> \param cell_ref ...
!> \param ewald_env ...
!> \param fist_nonbond_env ...
!> \param thermo ...
!> \param para_env ...
!> \param subsys ...
!> \param qmmm ...
!> \param qmmm_env ...
!> \param input ...
!> \param shell_model ...
!> \param shell_model_ad ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param multipoles ...
!> \param results ...
!> \param exclusions ...
!> \param efield ...
!> \par History
!>      11/03
!> \author CJM
! **************************************************************************************************
   SUBROUTINE fist_env_get(fist_env, atomic_kind_set, particle_set, ewald_pw, &
                           local_particles, local_molecules, molecule_kind_set, molecule_set, cell, &
                           cell_ref, ewald_env, fist_nonbond_env, thermo, para_env, subsys, qmmm, &
                           qmmm_env, input, shell_model, shell_model_ad, shell_particle_set, &
                           core_particle_set, multipoles, results, exclusions, efield)

      TYPE(fist_environment_type), INTENT(IN)            :: fist_env
      TYPE(atomic_kind_type), OPTIONAL, POINTER          :: atomic_kind_set(:)
      TYPE(particle_type), OPTIONAL, POINTER             :: particle_set(:)
      TYPE(ewald_pw_type), OPTIONAL, POINTER             :: ewald_pw
      TYPE(distribution_1d_type), OPTIONAL, POINTER      :: local_particles, local_molecules
      TYPE(molecule_kind_type), OPTIONAL, POINTER        :: molecule_kind_set(:)
      TYPE(molecule_type), OPTIONAL, POINTER             :: molecule_set(:)
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell, cell_ref
      TYPE(ewald_environment_type), OPTIONAL, POINTER    :: ewald_env
      TYPE(fist_nonbond_env_type), OPTIONAL, POINTER     :: fist_nonbond_env
      TYPE(fist_energy_type), OPTIONAL, POINTER          :: thermo
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: subsys
      LOGICAL, OPTIONAL                                  :: qmmm
      TYPE(qmmm_env_mm_type), OPTIONAL, POINTER          :: qmmm_env
      TYPE(section_vals_type), OPTIONAL, POINTER         :: input
      LOGICAL, OPTIONAL                                  :: shell_model, shell_model_ad
      TYPE(particle_type), OPTIONAL, POINTER             :: shell_particle_set(:), &
                                                            core_particle_set(:)
      TYPE(multipole_type), OPTIONAL, POINTER            :: multipoles
      TYPE(cp_result_type), OPTIONAL, POINTER            :: results
      TYPE(exclusion_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: exclusions
      TYPE(fist_efield_type), OPTIONAL, POINTER          :: efield

      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(multipole_type), POINTER                      :: fist_multipoles
      TYPE(particle_list_type), POINTER                  :: core_particles, particles, &
                                                            shell_particles

      NULLIFY (atomic_kinds, particles, molecules, molecule_kinds, fist_multipoles)

      IF (PRESENT(input)) input => fist_env%input
      IF (PRESENT(qmmm)) qmmm = fist_env%qmmm
      IF (PRESENT(qmmm_env)) qmmm_env => fist_env%qmmm_env
      IF (PRESENT(cell_ref)) cell_ref => fist_env%cell_ref
      IF (PRESENT(ewald_env)) ewald_env => fist_env%ewald_env
      IF (PRESENT(thermo)) thermo => fist_env%thermo
      IF (PRESENT(exclusions)) exclusions => fist_env%exclusions
      IF (PRESENT(para_env)) para_env => fist_env%para_env
      IF (PRESENT(ewald_pw)) ewald_pw => fist_env%ewald_pw
      IF (PRESENT(fist_nonbond_env)) fist_nonbond_env => fist_env%fist_nonbond_env
      IF (PRESENT(shell_model)) shell_model = fist_env%shell_model
      IF (PRESENT(shell_model_ad)) shell_model_ad = fist_env%shell_model_ad
      IF (PRESENT(subsys)) subsys => fist_env%subsys
      IF (PRESENT(efield)) efield => fist_env%efield

      IF (ASSOCIATED(fist_env%subsys)) &
         CALL cp_subsys_get(fist_env%subsys, &
                            atomic_kinds=atomic_kinds, &
                            local_molecules=local_molecules, &
                            local_particles=local_particles, &
                            particles=particles, &
                            molecule_kinds=molecule_kinds, &
                            molecules=molecules, &
                            shell_particles=shell_particles, &
                            core_particles=core_particles, &
                            multipoles=fist_multipoles, &
                            results=results, &
                            cell=cell)
      IF (PRESENT(atomic_kind_set)) atomic_kind_set => atomic_kinds%els
      IF (PRESENT(particle_set)) particle_set => particles%els
      IF (PRESENT(molecule_kind_set)) molecule_kind_set => molecule_kinds%els
      IF (PRESENT(molecule_set)) molecule_set => molecules%els
      IF (PRESENT(shell_particle_set)) shell_particle_set => shell_particles%els
      IF (PRESENT(core_particle_set)) core_particle_set => core_particles%els
      IF (PRESENT(multipoles)) multipoles => fist_multipoles
   END SUBROUTINE fist_env_get

! **************************************************************************************************
!> \brief Initialise the FIST environment.
!> \param fist_env the pointer to the fist_env
!> \param para_env ...
!> \par History
!>      11/03
!> \author CJM
! **************************************************************************************************
   SUBROUTINE init_fist_env(fist_env, para_env)

      TYPE(fist_environment_type), INTENT(OUT)           :: fist_env
      TYPE(mp_para_env_type), POINTER                    :: para_env

      NULLIFY (fist_env%input)
      NULLIFY (fist_env%qmmm_env)
      NULLIFY (fist_env%cell_ref)
      NULLIFY (fist_env%ewald_env)
      NULLIFY (fist_env%ewald_pw)
      NULLIFY (fist_env%thermo)
      NULLIFY (fist_env%fist_nonbond_env)
      NULLIFY (fist_env%subsys)
      NULLIFY (fist_env%exclusions)
      NULLIFY (fist_env%efield)
      fist_env%qmmm = .FALSE.
      fist_env%shell_model = .FALSE.
      fist_env%shell_model_ad = .FALSE.
      ALLOCATE (fist_env%qmmm_env)
      CALL qmmm_env_mm_create(fist_env%qmmm_env)
      NULLIFY (fist_env%subsys)
      CALL para_env%retain()
      fist_env%para_env => para_env

   END SUBROUTINE init_fist_env

! **************************************************************************************************
!> \brief Set the FIST environment.
!> \param fist_env the pointer to the fist_env
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param ewald_pw ...
!> \param local_particles ...
!> \param local_molecules ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param cell_ref ...
!> \param ewald_env ...
!> \param fist_nonbond_env ...
!> \param thermo ...
!> \param subsys ...
!> \param qmmm ...
!> \param qmmm_env ...
!> \param input ...
!> \param shell_model ...
!> \param shell_model_ad ...
!> \param exclusions ...
!> \param efield ...
!> \par History
!>      11/03
!> \author CJM
! **************************************************************************************************
   SUBROUTINE fist_env_set(fist_env, atomic_kind_set, particle_set, ewald_pw, &
                           local_particles, local_molecules, molecule_kind_set, &
                           molecule_set, cell_ref, ewald_env, &
                           fist_nonbond_env, thermo, subsys, qmmm, qmmm_env, &
                           input, shell_model, shell_model_ad, exclusions, efield)

      TYPE(fist_environment_type), INTENT(INOUT)         :: fist_env
      TYPE(atomic_kind_type), OPTIONAL, POINTER          :: atomic_kind_set(:)
      TYPE(particle_type), OPTIONAL, POINTER             :: particle_set(:)
      TYPE(ewald_pw_type), OPTIONAL, POINTER             :: ewald_pw
      TYPE(distribution_1d_type), OPTIONAL, POINTER      :: local_particles, local_molecules
      TYPE(molecule_kind_type), OPTIONAL, POINTER        :: molecule_kind_set(:)
      TYPE(molecule_type), OPTIONAL, POINTER             :: molecule_set(:)
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell_ref
      TYPE(ewald_environment_type), OPTIONAL, POINTER    :: ewald_env
      TYPE(fist_nonbond_env_type), OPTIONAL, POINTER     :: fist_nonbond_env
      TYPE(fist_energy_type), OPTIONAL, POINTER          :: thermo
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: subsys
      LOGICAL, OPTIONAL                                  :: qmmm
      TYPE(qmmm_env_mm_type), OPTIONAL, POINTER          :: qmmm_env
      TYPE(section_vals_type), OPTIONAL, POINTER         :: input
      LOGICAL, OPTIONAL                                  :: shell_model, shell_model_ad
      TYPE(exclusion_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: exclusions
      TYPE(fist_efield_type), OPTIONAL, POINTER          :: efield

      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(particle_list_type), POINTER                  :: particles

      IF (PRESENT(qmmm)) fist_env%qmmm = qmmm
      IF (PRESENT(qmmm_env)) THEN
         IF (ASSOCIATED(fist_env%qmmm_env)) THEN
            CALL qmmm_env_mm_release(fist_env%qmmm_env)
            DEALLOCATE (fist_env%qmmm_env)
         END IF
         fist_env%qmmm_env => qmmm_env
      END IF
      IF (PRESENT(ewald_env)) THEN
         IF (ASSOCIATED(fist_env%ewald_env)) THEN
         IF (.NOT. ASSOCIATED(fist_env%ewald_env, ewald_env)) THEN
            CALL ewald_env_release(fist_env%ewald_env)
            DEALLOCATE (fist_env%ewald_env)
         END IF
         END IF
         fist_env%ewald_env => ewald_env
      END IF
      IF (PRESENT(ewald_pw)) THEN
         IF (ASSOCIATED(fist_env%ewald_pw)) THEN
            IF (.NOT. ASSOCIATED(fist_env%ewald_pw, ewald_pw)) THEN
               CALL ewald_pw_release(fist_env%ewald_pw)
               DEALLOCATE (fist_env%ewald_pw)
            END IF
         END IF
         fist_env%ewald_pw => ewald_pw
      END IF
      IF (PRESENT(cell_ref)) THEN
         CALL cell_retain(cell_ref)
         CALL cell_release(fist_env%cell_ref)
         fist_env%cell_ref => cell_ref
      END IF
      IF (PRESENT(fist_nonbond_env)) THEN
         IF (ASSOCIATED(fist_env%fist_nonbond_env)) THEN
         IF (.NOT. ASSOCIATED(fist_env%fist_nonbond_env, fist_nonbond_env)) THEN
            CALL fist_nonbond_env_release(fist_env%fist_nonbond_env)
            DEALLOCATE (fist_env%fist_nonbond_env)
         END IF
         END IF
         fist_env%fist_nonbond_env => fist_nonbond_env
      END IF
      IF (PRESENT(input)) THEN
         CALL section_vals_retain(input)
         CALL section_vals_release(fist_env%input)
         fist_env%input => input
      END IF
      IF (PRESENT(thermo)) fist_env%thermo => thermo
      IF (PRESENT(subsys)) THEN
         IF (ASSOCIATED(fist_env%subsys)) THEN
         IF (.NOT. ASSOCIATED(fist_env%subsys, subsys)) THEN
            CALL cp_subsys_release(fist_env%subsys)
         END IF
         END IF
         fist_env%subsys => subsys
      END IF
      IF (PRESENT(atomic_kind_set)) THEN
         CALL atomic_kind_list_create(atomic_kinds, &
                                      els_ptr=atomic_kind_set)
         CALL cp_subsys_set(fist_env%subsys, &
                            atomic_kinds=atomic_kinds)
         CALL atomic_kind_list_release(atomic_kinds)
      END IF
      IF (PRESENT(particle_set)) THEN
         CALL particle_list_create(particles, &
                                   els_ptr=particle_set)
         CALL cp_subsys_set(fist_env%subsys, &
                            particles=particles)
         CALL particle_list_release(particles)
      END IF
      IF (PRESENT(local_particles)) THEN
         CALL cp_subsys_set(fist_env%subsys, &
                            local_particles=local_particles)
      END IF
      IF (PRESENT(local_molecules)) THEN
         CALL cp_subsys_set(fist_env%subsys, &
                            local_molecules=local_molecules)
      END IF
      IF (PRESENT(molecule_kind_set)) THEN
         CALL molecule_kind_list_create(molecule_kinds, &
                                        els_ptr=molecule_kind_set)
         CALL cp_subsys_set(fist_env%subsys, &
                            molecule_kinds=molecule_kinds)
         CALL molecule_kind_list_release(molecule_kinds)
      END IF
      IF (PRESENT(molecule_set)) THEN
         CALL molecule_list_create(molecules, &
                                   els_ptr=molecule_set)
         CALL cp_subsys_set(fist_env%subsys, &
                            molecules=molecules)
         CALL molecule_list_release(molecules)
      END IF
      IF (PRESENT(exclusions)) fist_env%exclusions => exclusions
      IF (PRESENT(shell_model)) THEN
         fist_env%shell_model = shell_model
      END IF
      IF (PRESENT(shell_model_ad)) THEN
         fist_env%shell_model_ad = shell_model_ad
      END IF
      IF (PRESENT(efield)) fist_env%efield => efield

   END SUBROUTINE fist_env_set

! **************************************************************************************************
!> \brief allocates and intitializes a fist_env
!> \param fist_env the object to create
!> \param para_env the parallel environment for the qs_env
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE fist_env_create(fist_env, para_env)
      TYPE(fist_environment_type), INTENT(OUT)           :: fist_env
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env

      CALL init_fist_env(fist_env, para_env=para_env)
   END SUBROUTINE fist_env_create

! **************************************************************************************************
!> \brief releases the given fist_env (see doc/ReferenceCounting.html)
!> \param fist_env the object to release
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE fist_env_release(fist_env)
      TYPE(fist_environment_type), INTENT(INOUT)         :: fist_env

      IF (ASSOCIATED(fist_env%qmmm_env)) THEN
         CALL qmmm_env_mm_release(fist_env%qmmm_env)
         DEALLOCATE (fist_env%qmmm_env)
      END IF
      CALL cell_release(fist_env%cell_ref)
      IF (ASSOCIATED(fist_env%ewald_pw)) THEN
         CALL ewald_pw_release(fist_env%ewald_pw)
         DEALLOCATE (fist_env%ewald_pw)
      END IF
      IF (ASSOCIATED(fist_env%ewald_env)) THEN
         CALL ewald_env_release(fist_env%ewald_env)
         DEALLOCATE (fist_env%ewald_env)
      END IF
      CALL mp_para_env_release(fist_env%para_env)
      CALL deallocate_fist_energy(fist_env%thermo)

      IF (ASSOCIATED(fist_env%fist_nonbond_env)) THEN
         CALL fist_nonbond_env_release(fist_env%fist_nonbond_env)
         DEALLOCATE (fist_env%fist_nonbond_env)
      END IF
      CALL cp_subsys_release(fist_env%subsys)
      CALL section_vals_release(fist_env%input)
      CALL exclusion_release(fist_env%exclusions)

      IF (ASSOCIATED(fist_env%efield)) THEN
         DEALLOCATE (fist_env%efield)
      END IF

   END SUBROUTINE fist_env_release

END MODULE fist_environment_types

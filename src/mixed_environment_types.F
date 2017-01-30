!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \author fschiff SEPT-11-06
! **************************************************************************************************
MODULE mixed_environment_types
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                              atomic_kind_list_release,&
                                              atomic_kind_list_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_release,&
                                              cell_retain,&
                                              cell_type
   USE cp_log_handling,                 ONLY: cp_logger_p_type,&
                                              cp_logger_release
   USE cp_para_env,                     ONLY: cp_para_env_release,&
                                              cp_para_env_retain
   USE cp_para_types,                   ONLY: cp_para_env_p_type,&
                                              cp_para_env_type
   USE cp_result_types,                 ONLY: cp_result_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_release,&
                                              cp_subsys_retain,&
                                              cp_subsys_set,&
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE input_section_types,             ONLY: section_vals_release,&
                                              section_vals_retain,&
                                              section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE mixed_cdft_types,                ONLY: mixed_cdft_type,&
                                              mixed_cdft_type_release
   USE mixed_energy_types,              ONLY: deallocate_mixed_energy,&
                                              mixed_energy_type
   USE molecule_kind_list_types,        ONLY: molecule_kind_list_create,&
                                              molecule_kind_list_release,&
                                              molecule_kind_list_type
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_list_types,             ONLY: molecule_list_create,&
                                              molecule_list_release,&
                                              molecule_list_type
   USE molecule_types,                  ONLY: molecule_type
   USE particle_list_types,             ONLY: particle_list_create,&
                                              particle_list_release,&
                                              particle_list_type
   USE particle_types,                  ONLY: particle_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

! **************************************************************************************************
!> \param mixed_env the pointer to the mixed_env
!> \par History
!>      11/06 Created [fschiff]
!>      12/15-12/16 Mixed CDFT [Nico Holmberg]
! **************************************************************************************************
   TYPE mixed_environment_type
      INTEGER                                          :: id_nr, ref_count
      TYPE(cell_type), POINTER                         :: cell_ref
      TYPE(mixed_energy_type), POINTER                 :: mixed_energy
      TYPE(cp_para_env_type), POINTER                  :: para_env
      TYPE(cp_subsys_type), POINTER                    :: subsys
      TYPE(section_vals_type), POINTER                 :: input
      REAL(KIND=dp), DIMENSION(:), POINTER             :: energies
      ! Parallelization of multiple force_eval
      INTEGER                                          :: new_group, ngroups
      INTEGER, DIMENSION(:), POINTER                   :: group_distribution
      TYPE(cp_para_env_p_type), DIMENSION(:), POINTER  :: sub_para_env
      TYPE(cp_logger_p_type), DIMENSION(:), POINTER    :: sub_logger
      REAL(KIND=dp), POINTER, DIMENSION(:)             :: val
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                         :: par
      REAL(KIND=dp)                                    :: dx, lerr
      CHARACTER(default_path_length)                   :: coupling_function
      ! Mixed CDFT control parameters
      LOGICAL                                          :: do_mixed_cdft, do_mixed_et, &
                                                          do_mixed_qmmm_cdft
      INTEGER                                          :: et_freq
      REAL(KIND=dp), DIMENSION(:), POINTER             :: strength
      TYPE(mixed_cdft_type), POINTER                   :: cdft_control
   END TYPE mixed_environment_type

! **************************************************************************************************
!> \brief to build arrays of pointers
!> \param mixed_env the pointer to the mixed_env
!> \par History
!>      11/06
!> \author fschiff
! **************************************************************************************************
   TYPE mixed_environment_p_type
      TYPE(mixed_environment_type), POINTER :: mixed_env
   END TYPE mixed_environment_p_type

! *** Public data types ***

   PUBLIC :: mixed_environment_type

! *** Public subroutines ***

   PUBLIC :: get_mixed_env, &
             set_mixed_env, &
             mixed_env_create, &
             mixed_env_release, &
             mixed_env_retain

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mixed_environment_types'
   INTEGER, PRIVATE, SAVE :: last_mixed_env_id_nr = 0

CONTAINS

! **************************************************************************************************
!> \brief Get the MIXED environment.
!> \param mixed_env the pointer to the mixed_env
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param local_molecules ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param cell ...
!> \param cell_ref ...
!> \param mixed_energy ...
!> \param para_env ...
!> \param sub_para_env ...
!> \param subsys ...
!> \param input ...
!> \param results ...
!> \param cdft_control ...
! **************************************************************************************************
   SUBROUTINE get_mixed_env(mixed_env, atomic_kind_set, particle_set, &
                            local_particles, local_molecules, molecule_kind_set, &
                            molecule_set, cell, cell_ref, &
                            mixed_energy, para_env, sub_para_env, subsys, &
                            input, results, cdft_control)

      TYPE(mixed_environment_type), INTENT(IN)           :: mixed_env
      TYPE(atomic_kind_type), OPTIONAL, POINTER          :: atomic_kind_set(:)
      TYPE(particle_type), OPTIONAL, POINTER             :: particle_set(:)
      TYPE(distribution_1d_type), OPTIONAL, POINTER      :: local_particles, local_molecules
      TYPE(molecule_kind_type), OPTIONAL, POINTER        :: molecule_kind_set(:)
      TYPE(molecule_type), OPTIONAL, POINTER             :: molecule_set(:)
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell, cell_ref
      TYPE(mixed_energy_type), OPTIONAL, POINTER         :: mixed_energy
      TYPE(cp_para_env_type), OPTIONAL, POINTER          :: para_env
      TYPE(cp_para_env_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: sub_para_env
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: subsys
      TYPE(section_vals_type), OPTIONAL, POINTER         :: input
      TYPE(cp_result_type), OPTIONAL, POINTER            :: results
      TYPE(mixed_cdft_type), OPTIONAL, POINTER           :: cdft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'get_mixed_env', routineP = moduleN//':'//routineN

      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(particle_list_type), POINTER                  :: particles

      NULLIFY (atomic_kinds, particles, molecules, molecule_kinds)
      CPASSERT(ASSOCIATED(mixed_env%subsys))

      IF (PRESENT(input)) input => mixed_env%input
      IF (PRESENT(cell_ref)) cell_ref => mixed_env%cell_ref
      IF (PRESENT(mixed_energy)) mixed_energy => mixed_env%mixed_energy
      IF (PRESENT(para_env)) para_env => mixed_env%para_env
      IF (PRESENT(sub_para_env)) sub_para_env => mixed_env%sub_para_env
      IF (PRESENT(cdft_control)) cdft_control => mixed_env%cdft_control
      IF (PRESENT(subsys)) subsys => mixed_env%subsys
      CALL cp_subsys_get(mixed_env%subsys, &
                         atomic_kinds=atomic_kinds, &
                         local_molecules=local_molecules, &
                         local_particles=local_particles, &
                         particles=particles, &
                         molecule_kinds=molecule_kinds, &
                         molecules=molecules, &
                         results=results, &
                         cell=cell)
      IF (PRESENT(atomic_kind_set)) atomic_kind_set => atomic_kinds%els
      IF (PRESENT(particle_set)) particle_set => particles%els
      IF (PRESENT(molecule_kind_set)) molecule_kind_set => molecule_kinds%els
      IF (PRESENT(molecule_set)) molecule_set => molecules%els

   END SUBROUTINE get_mixed_env

! **************************************************************************************************
!> \brief Initialise the MIXED environment.
!> \param mixed_env the pointer to the mixed_env
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE init_mixed_env(mixed_env, para_env)

      TYPE(mixed_environment_type), INTENT(OUT)          :: mixed_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      NULLIFY (mixed_env%input)
      NULLIFY (mixed_env%cell_ref)
      NULLIFY (mixed_env%mixed_energy)
      NULLIFY (mixed_env%subsys)
      NULLIFY (mixed_env%sub_para_env)
      NULLIFY (mixed_env%sub_logger)
      NULLIFY (mixed_env%energies)
      NULLIFY (mixed_env%par)
      NULLIFY (mixed_env%val)
      NULLIFY (mixed_env%subsys)
      NULLIFY (mixed_env%cdft_control)
      NULLIFY (mixed_env%strength)
      mixed_env%do_mixed_cdft = .FALSE.
      mixed_env%do_mixed_et = .FALSE.
      mixed_env%do_mixed_qmmm_cdft = .FALSE.
      mixed_env%et_freq = -1
      CALL cp_para_env_retain(para_env)
      mixed_env%para_env => para_env
      mixed_env%ref_count = 1
      last_mixed_env_id_nr = last_mixed_env_id_nr+1
      mixed_env%id_nr = last_mixed_env_id_nr

   END SUBROUTINE init_mixed_env

! **************************************************************************************************
!> \brief Set the MIXED environment.
!> \param mixed_env the pointer to the mixed_env
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param local_molecules ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param cell_ref ...
!> \param mixed_energy ...
!> \param subsys ...
!> \param input ...
!> \param sub_para_env ...
!> \param cdft_control ...
! **************************************************************************************************
   SUBROUTINE set_mixed_env(mixed_env, atomic_kind_set, particle_set, &
                            local_particles, local_molecules, molecule_kind_set, &
                            molecule_set, cell_ref, mixed_energy, subsys, &
                            input, sub_para_env, cdft_control)

      TYPE(mixed_environment_type), POINTER              :: mixed_env
      TYPE(atomic_kind_type), OPTIONAL, POINTER          :: atomic_kind_set(:)
      TYPE(particle_type), OPTIONAL, POINTER             :: particle_set(:)
      TYPE(distribution_1d_type), OPTIONAL, POINTER      :: local_particles, local_molecules
      TYPE(molecule_kind_type), OPTIONAL, POINTER        :: molecule_kind_set(:)
      TYPE(molecule_type), OPTIONAL, POINTER             :: molecule_set(:)
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell_ref
      TYPE(mixed_energy_type), OPTIONAL, POINTER         :: mixed_energy
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: subsys
      TYPE(section_vals_type), OPTIONAL, POINTER         :: input
      TYPE(cp_para_env_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: sub_para_env
      TYPE(mixed_cdft_type), OPTIONAL, POINTER           :: cdft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'set_mixed_env', routineP = moduleN//':'//routineN

      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(particle_list_type), POINTER                  :: particles

      CPASSERT(mixed_env%ref_count > 0)

      IF (PRESENT(cell_ref)) THEN
         CALL cell_retain(cell_ref)
         CALL cell_release(mixed_env%cell_ref)
         mixed_env%cell_ref => cell_ref
      END IF
      IF (PRESENT(input)) THEN
         CALL section_vals_retain(input)
         CALL section_vals_release(mixed_env%input)
         mixed_env%input => input
      END IF
      IF (PRESENT(mixed_energy)) mixed_env%mixed_energy => mixed_energy
      IF (PRESENT(subsys)) THEN
         CALL cp_subsys_retain(subsys)
         CALL cp_subsys_release(mixed_env%subsys)
         mixed_env%subsys => subsys
      END IF
      IF (PRESENT(sub_para_env)) THEN
         mixed_env%sub_para_env => sub_para_env
      END IF
      IF (PRESENT(cdft_control)) mixed_env%cdft_control => cdft_control
      IF (PRESENT(atomic_kind_set)) THEN
         CALL atomic_kind_list_create(atomic_kinds, &
                                      els_ptr=atomic_kind_set)
         CALL cp_subsys_set(mixed_env%subsys, &
                            atomic_kinds=atomic_kinds)
         CALL atomic_kind_list_release(atomic_kinds)
      END IF
      IF (PRESENT(particle_set)) THEN
         CALL particle_list_create(particles, &
                                   els_ptr=particle_set)
         CALL cp_subsys_set(mixed_env%subsys, &
                            particles=particles)
         CALL particle_list_release(particles)
      END IF
      IF (PRESENT(local_particles)) THEN
         CALL cp_subsys_set(mixed_env%subsys, &
                            local_particles=local_particles)
      END IF
      IF (PRESENT(local_molecules)) THEN
         CALL cp_subsys_set(mixed_env%subsys, &
                            local_molecules=local_molecules)
      END IF
      IF (PRESENT(molecule_kind_set)) THEN
         CALL molecule_kind_list_create(molecule_kinds, els_ptr=molecule_kind_set)
         CALL cp_subsys_set(mixed_env%subsys, molecule_kinds=molecule_kinds)
         CALL molecule_kind_list_release(molecule_kinds)
      END IF
      IF (PRESENT(molecule_set)) THEN
         CALL molecule_list_create(molecules, els_ptr=molecule_set)
         CALL cp_subsys_set(mixed_env%subsys, molecules=molecules)
         CALL molecule_list_release(molecules)
      END IF

   END SUBROUTINE set_mixed_env

! **************************************************************************************************
!> \brief allocates and intitializes a mixed_env
!> \param mixed_env the object to create
!> \param para_env the parallel environement for the qs_env
!> \author fschiff 11.06
! **************************************************************************************************
   SUBROUTINE mixed_env_create(mixed_env, para_env)
      TYPE(mixed_environment_type), POINTER              :: mixed_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_env_create', &
         routineP = moduleN//':'//routineN

      ALLOCATE (mixed_env)
      CALL init_mixed_env(mixed_env, para_env=para_env)
   END SUBROUTINE mixed_env_create

! **************************************************************************************************
!> \brief retains the given mixed_env (see doc/ReferenceCounting.html)
!> \param mixed_env the object to retain
!> \author fschiff 11.06
! **************************************************************************************************
   SUBROUTINE mixed_env_retain(mixed_env)
      TYPE(mixed_environment_type), POINTER              :: mixed_env

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_env_retain', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(mixed_env))
      CPASSERT(mixed_env%ref_count > 0)
      mixed_env%ref_count = mixed_env%ref_count+1
   END SUBROUTINE mixed_env_retain

! **************************************************************************************************
!> \brief releases the given mixed_env (see doc/ReferenceCounting.html)
!> \param mixed_env the object to release
!> \author fschiff 11.06
! **************************************************************************************************
   SUBROUTINE mixed_env_release(mixed_env)
      TYPE(mixed_environment_type), POINTER              :: mixed_env

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_env_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ngroups

      IF (ASSOCIATED(mixed_env)) THEN
         CPASSERT(mixed_env%ref_count > 0)
         mixed_env%ref_count = mixed_env%ref_count-1
         IF (mixed_env%ref_count < 1) THEN
            ngroups = SIZE(mixed_env%sub_para_env)
            DO i = 1, ngroups
               IF (ASSOCIATED(mixed_env%sub_para_env(i)%para_env)) THEN
                  CALL cp_logger_release(mixed_env%sub_logger(i)%p)
                  CALL cp_para_env_release(mixed_env%sub_para_env(i)%para_env)
               END IF
            END DO
            DEALLOCATE (mixed_env%sub_para_env)
            DEALLOCATE (mixed_env%sub_logger)
            DEALLOCATE (mixed_env%energies)
            IF (ASSOCIATED(mixed_env%par)) THEN
               DEALLOCATE (mixed_env%par)
            END IF
            IF (ASSOCIATED(mixed_env%val)) THEN
               DEALLOCATE (mixed_env%val)
            END IF
            CALL cell_release(mixed_env%cell_ref)
            CALL cp_para_env_release(mixed_env%para_env)
            CALL deallocate_mixed_energy(mixed_env%mixed_energy)
            CALL cp_subsys_release(mixed_env%subsys)
            CALL section_vals_release(mixed_env%input)
            IF (ASSOCIATED(mixed_env%group_distribution)) THEN
               DEALLOCATE (mixed_env%group_distribution)
            END IF
            IF (ASSOCIATED(mixed_env%cdft_control)) &
               CALL mixed_cdft_type_release(mixed_env%cdft_control)
            IF (ASSOCIATED(mixed_env%strength)) &
               DEALLOCATE (mixed_env%strength)
            DEALLOCATE (mixed_env)
         END IF
      END IF
      NULLIFY (mixed_env)

   END SUBROUTINE mixed_env_release

END MODULE mixed_environment_types

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author fschiff SEPT-11-06
! *****************************************************************************
MODULE mixed_environment_types
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_methods,                    ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_p_type,&
                                             cp_para_env_type
  USE cp_result_types,                 ONLY: cp_result_create,&
                                             cp_result_release,&
                                             cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_create,&
                                             cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_retain,&
                                             cp_subsys_set,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mixed_energy_types,              ONLY: deallocate_mixed_energy,&
                                             mixed_energy_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_create,&
                                             mol_kind_new_list_release,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_create,&
                                             mol_new_list_release,&
                                             mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

! *****************************************************************************
!> \param mixed_env the pointer to the mixed_env
!> \par History
!>      11/06
!> \author fschiff
! *****************************************************************************
  TYPE mixed_environment_type
     INTEGER :: id_nr, ref_count
     TYPE ( cell_type ), POINTER                      :: cell
     TYPE ( cell_type ), POINTER                      :: cell_ref
     TYPE ( mixed_energy_type ), POINTER              :: mixed_energy
     TYPE ( cp_para_env_type ), POINTER               :: para_env
     TYPE ( cp_subsys_type ), POINTER                 :: subsys
     TYPE(section_vals_type), POINTER                 :: input
     REAL(KIND=dp), DIMENSION(:), POINTER             :: energies
     TYPE(cp_result_type),POINTER                     :: results
     ! Parallelization of multiple force_eval
     INTEGER                                          :: new_group, ngroups
     INTEGER, DIMENSION(:), POINTER                   :: group_distribution
     TYPE (cp_para_env_p_type),DIMENSION(:), POINTER  :: sub_para_env
     TYPE (cp_error_type), DIMENSION(:), POINTER      :: sub_error
     REAL(KIND=dp), POINTER, DIMENSION(:)             :: val
     CHARACTER(LEN=default_string_length),&
          DIMENSION(:), POINTER                       :: par
     REAL(KIND=dp)                                    :: dx, lerr
     CHARACTER(default_path_length)                   :: coupling_function
  END TYPE mixed_environment_type

! *****************************************************************************
!> \brief to build arrays of pointers
!> \param mixed_env the pointer to the mixed_env
!> \par History
!>      11/06
!> \author fschiff
! *****************************************************************************
  TYPE mixed_environment_p_type
     TYPE(mixed_environment_type), POINTER :: mixed_env
  END TYPE mixed_environment_p_type

! *** Public data types ***

  PUBLIC :: mixed_environment_type

! *** Public subroutines ***

  PUBLIC :: get_mixed_env,&
            set_mixed_env,&
            mixed_env_create,&
            mixed_env_release,&
            mixed_env_retain

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mixed_environment_types'
  INTEGER, PRIVATE, SAVE :: last_mixed_env_id_nr=0

CONTAINS

! *****************************************************************************
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
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_mixed_env( mixed_env, atomic_kind_set, particle_set,  &
                           local_particles, local_molecules, molecule_kind_set, &
                           molecule_set, cell, cell_ref, &
                           mixed_energy, para_env, sub_para_env, subsys,  &
                           input, results, error )

    TYPE(mixed_environment_type), INTENT(IN) :: mixed_env
    TYPE(atomic_kind_type), OPTIONAL, &
      POINTER                                :: atomic_kind_set( : )
    TYPE(particle_type), OPTIONAL, POINTER   :: particle_set( : )
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles, &
                                                local_molecules
    TYPE(molecule_kind_type), OPTIONAL, &
      POINTER                                :: molecule_kind_set( : )
    TYPE(molecule_type), OPTIONAL, POINTER   :: molecule_set( : )
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    TYPE(mixed_energy_type), OPTIONAL, &
      POINTER                                :: mixed_energy
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_para_env_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: sub_para_env
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    TYPE(cp_result_type), OPTIONAL, POINTER  :: results
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_mixed_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.
    NULLIFY( atomic_kinds, particles, molecules_new, molecule_kinds_new )
    CPPrecondition(ASSOCIATED(mixed_env%subsys),cp_failure_level,routineP,error,failure)

    IF ( PRESENT ( input ) ) input => mixed_env%input
    IF ( PRESENT ( cell ) ) cell => mixed_env % cell
    IF ( PRESENT ( cell_ref ) ) cell_ref => mixed_env % cell_ref
    IF ( PRESENT ( mixed_energy ) ) mixed_energy => mixed_env % mixed_energy
    IF ( PRESENT ( para_env ) ) para_env => mixed_env % para_env
    IF ( PRESENT ( sub_para_env ) ) sub_para_env => mixed_env % sub_para_env
    IF (PRESENT(subsys)) subsys => mixed_env%subsys
    CALL cp_subsys_get(mixed_env%subsys,&
                      atomic_kinds=atomic_kinds,&
                      local_molecules_new=local_molecules,&
                      local_particles=local_particles,&
                      particles=particles,&
                      molecule_kinds_new=molecule_kinds_new,&
                      molecules_new=molecules_new,&
                      error=error)
    IF (PRESENT(atomic_kind_set)) atomic_kind_set => atomic_kinds%els
    IF (PRESENT(particle_set)) particle_set => particles%els
    IF (PRESENT(molecule_kind_set)) molecule_kind_set => molecule_kinds_new%els
    IF (PRESENT(molecule_set)) molecule_set => molecules_new%els
    IF (PRESENT(results)) results => mixed_env%results

  END SUBROUTINE get_mixed_env

! *****************************************************************************
!> \brief Initialise the MIXED environment.
!> \param mixed_env the pointer to the mixed_env
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE init_mixed_env ( mixed_env, para_env, error )

    TYPE(mixed_environment_type), &
      INTENT(OUT)                            :: mixed_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    NULLIFY (mixed_env%input)
    NULLIFY (mixed_env%cell)
    NULLIFY (mixed_env%cell_ref)
    NULLIFY (mixed_env%mixed_energy)
    NULLIFY (mixed_env%subsys)
    NULLIFY (mixed_env%sub_para_env)
    NULLIFY (mixed_env%sub_error)
    NULLIFY (mixed_env%energies)
    NULLIFY (mixed_env%par)
    NULLIFY (mixed_env%val)
    CALL cp_subsys_create(mixed_env%subsys, para_env=para_env, error=error)
    CALL cp_result_create(results=mixed_env%results,error=error)
    CALL cp_para_env_retain(para_env,error=error)
    mixed_env%para_env => para_env
    mixed_env%ref_count=1
    last_mixed_env_id_nr=last_mixed_env_id_nr+1
    mixed_env%id_nr=last_mixed_env_id_nr

  END SUBROUTINE init_mixed_env

! *****************************************************************************
!> \brief Set the MIXED environment.
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
!> \param subsys ...
!> \param input ...
!> \param sub_para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE set_mixed_env( mixed_env, atomic_kind_set, particle_set, &
                           local_particles, local_molecules, molecule_kind_set, &
                           molecule_set, cell, cell_ref, mixed_energy, subsys, &
                           input, sub_para_env, error )

    TYPE(mixed_environment_type), POINTER    :: mixed_env
    TYPE(atomic_kind_type), OPTIONAL, &
      POINTER                                :: atomic_kind_set( : )
    TYPE(particle_type), OPTIONAL, POINTER   :: particle_set( : )
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles, &
                                                local_molecules
    TYPE(molecule_kind_type), OPTIONAL, &
      POINTER                                :: molecule_kind_set( : )
    TYPE(molecule_type), OPTIONAL, POINTER   :: molecule_set( : )
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    TYPE(mixed_energy_type), OPTIONAL, &
      POINTER                                :: mixed_energy
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    TYPE(cp_para_env_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: sub_para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_mixed_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.
    CPPrecondition(mixed_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(cell)) THEN
       CALL cell_retain(cell, error=error)
       CALL cell_release(mixed_env%cell,error=error)
       mixed_env%cell => cell
    END IF

    IF (PRESENT(cell_ref)) THEN
       CALL cell_retain(cell_ref, error=error)
       CALL cell_release(mixed_env%cell_ref,error=error)
       mixed_env%cell_ref => cell_ref
    END IF
    IF (PRESENT(input)) THEN
       CALL section_vals_retain(input,error=error)
       CALL section_vals_release(mixed_env%input,error=error)
       mixed_env%input => input
    END IF
    IF ( PRESENT ( mixed_energy ) ) mixed_env % mixed_energy => mixed_energy
    IF (PRESENT(subsys)) THEN
       CALL cp_subsys_retain(subsys,error=error)
       CALL cp_subsys_release(mixed_env%subsys,error=error)
       mixed_env%subsys => subsys
    END IF
    IF (PRESENT(sub_para_env)) THEN
       mixed_env%sub_para_env => sub_para_env
    END IF
    IF (PRESENT(atomic_kind_set)) THEN
      CALL atomic_kind_list_create(atomic_kinds,&
                                   els_ptr=atomic_kind_set,&
                                   error=error)
      CALL cp_subsys_set(mixed_env%subsys,&
                        atomic_kinds=atomic_kinds,&
                        error=error)
      CALL atomic_kind_list_release(atomic_kinds,error=error)
    END IF
    IF (PRESENT(particle_set)) THEN
      CALL particle_list_create(particles,&
                                els_ptr=particle_set,&
                                error=error)
      CALL cp_subsys_set(mixed_env%subsys,&
                        particles=particles,&
                        error=error)
      CALL particle_list_release(particles,error=error)
    END IF
    IF (PRESENT(local_particles)) THEN
      CALL cp_subsys_set(mixed_env%subsys,&
                        local_particles=local_particles,&
                        error=error)
    END IF
    IF (PRESENT(local_molecules)) THEN
      CALL cp_subsys_set(mixed_env%subsys,&
                        local_molecules_new=local_molecules,&
                        error=error)
    END IF
    IF (PRESENT(molecule_kind_set)) THEN
      CALL mol_kind_new_list_create(molecule_kinds_new,&
                                    els_ptr=molecule_kind_set,&
                                    error=error)
      CALL cp_subsys_set(mixed_env%subsys,&
                        molecule_kinds_new=molecule_kinds_new,&
                        error=error)
      CALL mol_kind_new_list_release(molecule_kinds_new,error=error)
    END IF
    IF (PRESENT(molecule_set)) THEN
      CALL mol_new_list_create(molecules_new,&
                               els_ptr=molecule_set,&
                               error=error)
      CALL cp_subsys_set(mixed_env%subsys,&
                        molecules_new=molecules_new,&
                        error=error)
      CALL mol_new_list_release(molecules_new,error=error)
    END IF

  END SUBROUTINE set_mixed_env

! *****************************************************************************
!> \brief allocates and intitializes a mixed_env
!> \param mixed_env the object to create
!> \param para_env the parallel environement for the qs_env
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fschiff 11.06
! *****************************************************************************
  SUBROUTINE mixed_env_create(mixed_env,para_env,error)
    TYPE(mixed_environment_type), POINTER    :: mixed_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mixed_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(mixed_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL init_mixed_env(mixed_env,para_env=para_env, error=error)
    END IF
  END SUBROUTINE mixed_env_create

! *****************************************************************************
!> \brief retains the given mixed_env (see doc/ReferenceCounting.html)
!> \param mixed_env the object to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fschiff 11.06
! *****************************************************************************
  SUBROUTINE mixed_env_retain(mixed_env,error)
    TYPE(mixed_environment_type), POINTER    :: mixed_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mixed_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(mixed_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(mixed_env%ref_count>0,cp_failure_level,routineP,error,failure)
       mixed_env%ref_count=mixed_env%ref_count+1
    END IF
  END SUBROUTINE mixed_env_retain

! *****************************************************************************
!> \brief releases the given mixed_env (see doc/ReferenceCounting.html)
!> \param mixed_env the object to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fschiff 11.06
! *****************************************************************************
  SUBROUTINE mixed_env_release(mixed_env,error)
    TYPE(mixed_environment_type), POINTER    :: mixed_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mixed_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ngroups, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(mixed_env)) THEN
       CPPrecondition(mixed_env%ref_count>0,cp_failure_level,routineP,error,failure)
       mixed_env%ref_count=mixed_env%ref_count-1
       IF (mixed_env%ref_count<1) THEN
          CALL cell_release(mixed_env%cell,error=error)
          ngroups = SIZE(mixed_env%sub_para_env)
          DO i = 1, ngroups
             IF (ASSOCIATED(mixed_env%sub_para_env(i)%para_env)) THEN
                CALL cp_error_dealloc_ref(mixed_env%sub_error(i), error)
                CALL cp_para_env_release(mixed_env%sub_para_env(i)%para_env,error=error)
             END IF
          END DO
          DEALLOCATE(mixed_env%sub_para_env, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(mixed_env%sub_error, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(mixed_env%energies, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (ASSOCIATED(mixed_env%par)) THEN
             DEALLOCATE(mixed_env%par, stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(mixed_env%val)) THEN
             DEALLOCATE(mixed_env%val, stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          CALL cell_release(mixed_env%cell_ref,error=error)
          CALL cp_result_release(mixed_env%results,error=error)
          CALL cp_para_env_release(mixed_env%para_env,error=error)
          CALL deallocate_mixed_energy(mixed_env%mixed_energy,error=error)
          CALL cp_subsys_release(mixed_env%subsys,error=error)
          CALL section_vals_release(mixed_env%input,error=error)
          IF (ASSOCIATED(mixed_env%group_distribution)) THEN
             DEALLOCATE(mixed_env%group_distribution, stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          DEALLOCATE(mixed_env, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(mixed_env)

  END SUBROUTINE mixed_env_release

END MODULE mixed_environment_types

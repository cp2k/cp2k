!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>    gt sept-23-02 added atomic_kind_set to replica_environment_type
!>    to allow use of kind_based neighbor list
!>    CJM rewrite
!> \author CJM SEPT-01-02
! *****************************************************************************
MODULE fist_environment_types
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_methods,                    ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_types,                 ONLY: cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_retain,&
                                             cp_subsys_set,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_env_release,&
                                             ewald_env_retain,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_release,&
                                             ewald_pw_retain,&
                                             ewald_pw_type
  USE exclusion_types,                 ONLY: exclusion_release,&
                                             exclusion_type
  USE f77_blas
  USE fist_energy_types,               ONLY: deallocate_fist_energy,&
                                             fist_energy_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_release,&
                                             fist_nonbond_env_retain,&
                                             fist_nonbond_env_type
  USE input_section_types,             ONLY: section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_create,&
                                             mol_kind_new_list_release,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_create,&
                                             mol_new_list_release,&
                                             mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE multipole_types,                 ONLY: multipole_type
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: qmmm_env_mm_create,&
                                             qmmm_env_mm_release,&
                                             qmmm_env_mm_retain,&
                                             qmmm_env_mm_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

! *****************************************************************************
!> \par History
!>      11/03
!> \author CJM
! *****************************************************************************
  TYPE fist_environment_type
    PRIVATE
    INTEGER                                      :: id_nr, ref_count
    LOGICAL                                      :: qmmm
    LOGICAL                                      :: shell_model, shell_model_ad
    TYPE (qmmm_env_mm_type), POINTER             :: qmmm_env
    TYPE (cell_type), POINTER                    :: cell_ref
    TYPE (ewald_environment_type), POINTER       :: ewald_env
    TYPE (ewald_pw_type),          POINTER       :: ewald_pw
    TYPE (fist_energy_type), POINTER             :: thermo
    TYPE (cp_para_env_type), POINTER             :: para_env
    TYPE (cp_subsys_type), POINTER               :: subsys
    TYPE (fist_nonbond_env_type), POINTER        :: fist_nonbond_env
    TYPE(section_vals_type), POINTER             :: input
    TYPE(exclusion_type), DIMENSION(:), POINTER  :: exclusions
 END TYPE fist_environment_type

! *****************************************************************************
!> \brief to build arrays of pointers
!> \param fist_env the pointer to the fist_env
!> \par History
!>      11/03
!> \author CJM
! *****************************************************************************
  TYPE fist_environment_p_type
     TYPE(fist_environment_type), POINTER :: fist_env
  END TYPE fist_environment_p_type

! *** Public data types ***
  PUBLIC :: fist_environment_type

! *** Public subroutines ***
  PUBLIC :: fist_env_get,&
            fist_env_set,&
            fist_env_create,&
            fist_env_release,&
            fist_env_retain

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_environment_types'
  INTEGER, PRIVATE, SAVE :: last_fist_env_id_nr=0

!***

CONTAINS

! *****************************************************************************
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
!> \param error ...
!> \par History
!>      11/03
!> \author CJM
! *****************************************************************************
  SUBROUTINE fist_env_get( fist_env, atomic_kind_set, particle_set, ewald_pw,  &
       local_particles, local_molecules, molecule_kind_set, molecule_set, cell,&
       cell_ref, ewald_env, fist_nonbond_env, thermo, para_env, subsys, qmmm,&
       qmmm_env, input, shell_model, shell_model_ad, shell_particle_set,&
       core_particle_set, multipoles, results, exclusions, error )

    TYPE(fist_environment_type), INTENT(IN)  :: fist_env
    TYPE(atomic_kind_type), OPTIONAL, &
      POINTER                                :: atomic_kind_set( : )
    TYPE(particle_type), OPTIONAL, POINTER   :: particle_set( : )
    TYPE(ewald_pw_type), OPTIONAL, POINTER   :: ewald_pw
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles, &
                                                local_molecules
    TYPE(molecule_kind_type), OPTIONAL, &
      POINTER                                :: molecule_kind_set( : )
    TYPE(molecule_type), OPTIONAL, POINTER   :: molecule_set( : )
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    TYPE(ewald_environment_type), OPTIONAL, &
      POINTER                                :: ewald_env
    TYPE(fist_nonbond_env_type), OPTIONAL, &
      POINTER                                :: fist_nonbond_env
    TYPE(fist_energy_type), OPTIONAL, &
      POINTER                                :: thermo
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    LOGICAL, OPTIONAL                        :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    LOGICAL, OPTIONAL                        :: shell_model, shell_model_ad
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set( : ), &
                                                core_particle_set( : )
    TYPE(multipole_type), OPTIONAL, POINTER  :: multipoles
    TYPE(cp_result_type), OPTIONAL, POINTER  :: results
    TYPE(exclusion_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: exclusions
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(multipole_type), POINTER            :: fist_multipoles
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles

    failure=.FALSE.
    NULLIFY( atomic_kinds, particles, molecules_new, molecule_kinds_new, fist_multipoles)

    IF ( PRESENT ( input ) ) input => fist_env%input
    IF ( PRESENT ( qmmm ) ) qmmm = fist_env % qmmm
    IF ( PRESENT ( qmmm_env )) qmmm_env => fist_env % qmmm_env
    IF ( PRESENT ( cell_ref ) ) cell_ref => fist_env % cell_ref
    IF ( PRESENT ( ewald_env ) ) ewald_env => fist_env % ewald_env
    IF ( PRESENT ( thermo ) ) thermo => fist_env % thermo
    IF ( PRESENT ( exclusions ) ) exclusions => fist_env % exclusions
    IF ( PRESENT ( para_env ) ) para_env => fist_env % para_env
    IF ( PRESENT ( ewald_pw ) ) ewald_pw => fist_env % ewald_pw
    IF ( PRESENT ( fist_nonbond_env ) ) fist_nonbond_env => fist_env % fist_nonbond_env
    IF ( PRESENT ( shell_model ) ) shell_model = fist_env % shell_model
    IF ( PRESENT ( shell_model_ad ) ) shell_model_ad = fist_env % shell_model_ad
    IF (PRESENT(subsys)) subsys => fist_env%subsys

    IF(ASSOCIATED(fist_env%subsys)) &
       CALL cp_subsys_get(fist_env%subsys,&
                          atomic_kinds=atomic_kinds,&
                          local_molecules_new=local_molecules,&
                          local_particles=local_particles,&
                          particles=particles,&
                          molecule_kinds_new=molecule_kinds_new,&
                          molecules_new=molecules_new,&
                          shell_particles=shell_particles,&
                          core_particles=core_particles,&
                          multipoles=fist_multipoles,&
                          results=results,&
                          cell=cell,&
                          error=error)
    IF (PRESENT(atomic_kind_set)) atomic_kind_set => atomic_kinds%els
    IF (PRESENT(particle_set)) particle_set => particles%els
    IF (PRESENT(molecule_kind_set)) molecule_kind_set => molecule_kinds_new%els
    IF (PRESENT(molecule_set)) molecule_set => molecules_new%els
    IF (PRESENT(shell_particle_set)) shell_particle_set => shell_particles%els
    IF (PRESENT(core_particle_set)) core_particle_set => core_particles%els
    IF (PRESENT(multipoles)) multipoles => fist_multipoles
  END SUBROUTINE fist_env_get

! *****************************************************************************
!> \brief Initialise the FIST environment.
!> \param fist_env the pointer to the fist_env
!> \param para_env ...
!> \param error ...
!> \par History
!>      11/03
!> \author CJM
! *****************************************************************************
  SUBROUTINE init_fist_env ( fist_env, para_env, error )

    TYPE(fist_environment_type), INTENT(OUT) :: fist_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    NULLIFY (fist_env%input)
    NULLIFY (fist_env%qmmm_env)
    NULLIFY (fist_env%cell_ref)
    NULLIFY (fist_env%ewald_env)
    NULLIFY (fist_env%ewald_pw)
    NULLIFY (fist_env%thermo)
    NULLIFY (fist_env%fist_nonbond_env)
    NULLIFY (fist_env%subsys)
    NULLIFY (fist_env%exclusions)
    fist_env % qmmm = .FALSE.
    fist_env % shell_model = .FALSE.
    fist_env % shell_model_ad = .FALSE.
    CALL qmmm_env_mm_create(fist_env%qmmm_env, error)
    NULLIFY (fist_env%subsys)
    CALL cp_para_env_retain(para_env,error=error)
    fist_env%para_env => para_env
    fist_env%ref_count=1
    last_fist_env_id_nr=last_fist_env_id_nr+1
    fist_env%id_nr=last_fist_env_id_nr

  END SUBROUTINE init_fist_env

! *****************************************************************************
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
!> \param error ...
!> \par History
!>      11/03
!> \author CJM
! *****************************************************************************
  SUBROUTINE fist_env_set( fist_env, atomic_kind_set, particle_set, ewald_pw, &
                           local_particles, local_molecules, molecule_kind_set, &
                           molecule_set, cell_ref, ewald_env,  &
                           fist_nonbond_env, thermo, subsys, qmmm, qmmm_env, &
                           input, shell_model, shell_model_ad, exclusions, error )

    TYPE(fist_environment_type), POINTER     :: fist_env
    TYPE(atomic_kind_type), OPTIONAL, &
      POINTER                                :: atomic_kind_set( : )
    TYPE(particle_type), OPTIONAL, POINTER   :: particle_set( : )
    TYPE(ewald_pw_type), OPTIONAL, POINTER   :: ewald_pw
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles, &
                                                local_molecules
    TYPE(molecule_kind_type), OPTIONAL, &
      POINTER                                :: molecule_kind_set( : )
    TYPE(molecule_type), OPTIONAL, POINTER   :: molecule_set( : )
    TYPE(cell_type), OPTIONAL, POINTER       :: cell_ref
    TYPE(ewald_environment_type), OPTIONAL, &
      POINTER                                :: ewald_env
    TYPE(fist_nonbond_env_type), OPTIONAL, &
      POINTER                                :: fist_nonbond_env
    TYPE(fist_energy_type), OPTIONAL, &
      POINTER                                :: thermo
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    LOGICAL, OPTIONAL                        :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    LOGICAL, OPTIONAL                        :: shell_model, shell_model_ad
    TYPE(exclusion_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: exclusions
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_env_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.
    CPPrecondition(fist_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(qmmm)) fist_env%qmmm = qmmm
    IF (PRESENT(qmmm_env)) THEN
       CALL qmmm_env_mm_retain(qmmm_env, error)
       CALL qmmm_env_mm_release(fist_env%qmmm_env, error=error)
       fist_env%qmmm_env => qmmm_env
    END IF
    IF ( PRESENT ( ewald_env ) ) THEN
       CALL ewald_env_retain ( ewald_env, error = error )
       CALL ewald_env_release ( fist_env % ewald_env, error = error )
       fist_env % ewald_env => ewald_env
    ENDIF
    IF ( PRESENT ( ewald_pw ) ) THEN
       CALL ewald_pw_retain ( ewald_pw, error = error )
       CALL ewald_pw_release ( fist_env % ewald_pw, error = error )
       fist_env % ewald_pw => ewald_pw
    ENDIF
    IF (PRESENT(cell_ref)) THEN
       CALL cell_retain(cell_ref, error=error)
       CALL cell_release(fist_env%cell_ref,error=error)
       fist_env%cell_ref => cell_ref
    END IF
    IF ( PRESENT ( fist_nonbond_env ) ) THEN
       CALL fist_nonbond_env_retain ( fist_nonbond_env, error )
       CALL fist_nonbond_env_release ( fist_env % fist_nonbond_env, error )
       fist_env % fist_nonbond_env => fist_nonbond_env
    ENDIF
    IF (PRESENT(input)) THEN
       CALL section_vals_retain(input,error=error)
       CALL section_vals_release(fist_env%input,error=error)
       fist_env%input => input
    END IF
    IF ( PRESENT ( thermo ) ) fist_env % thermo => thermo
    IF (PRESENT(subsys)) THEN
       CALL cp_subsys_retain(subsys,error=error)
       CALL cp_subsys_release(fist_env%subsys,error=error)
       fist_env%subsys => subsys
    END IF
    IF (PRESENT(atomic_kind_set)) THEN
      CALL atomic_kind_list_create(atomic_kinds,&
                                   els_ptr=atomic_kind_set,&
                                   error=error)
      CALL cp_subsys_set(fist_env%subsys,&
                        atomic_kinds=atomic_kinds,&
                        error=error)
      CALL atomic_kind_list_release(atomic_kinds,error=error)
    END IF
    IF (PRESENT(particle_set)) THEN
      CALL particle_list_create(particles,&
                                els_ptr=particle_set,&
                                error=error)
      CALL cp_subsys_set(fist_env%subsys,&
                        particles=particles,&
                        error=error)
      CALL particle_list_release(particles,error=error)
    END IF
    IF (PRESENT(local_particles)) THEN
      CALL cp_subsys_set(fist_env%subsys,&
                        local_particles=local_particles,&
                        error=error)
    END IF
    IF (PRESENT(local_molecules)) THEN
      CALL cp_subsys_set(fist_env%subsys,&
                        local_molecules_new=local_molecules,&
                        error=error)
    END IF
    IF (PRESENT(molecule_kind_set)) THEN
      CALL mol_kind_new_list_create(molecule_kinds_new,&
                                    els_ptr=molecule_kind_set,&
                                    error=error)
      CALL cp_subsys_set(fist_env%subsys,&
                        molecule_kinds_new=molecule_kinds_new,&
                        error=error)
      CALL mol_kind_new_list_release(molecule_kinds_new,error=error)
    END IF
    IF (PRESENT(molecule_set)) THEN
      CALL mol_new_list_create(molecules_new,&
                               els_ptr=molecule_set,&
                               error=error)
      CALL cp_subsys_set(fist_env%subsys,&
                        molecules_new=molecules_new,&
                        error=error)
      CALL mol_new_list_release(molecules_new,error=error)
    END IF
    IF (PRESENT(exclusions)) fist_env%exclusions=>exclusions
    IF (PRESENT(shell_model)) THEN
       fist_env%shell_model = shell_model
    END IF
    IF (PRESENT(shell_model_ad)) THEN
       fist_env%shell_model_ad = shell_model_ad
    END IF

  END SUBROUTINE fist_env_set

! *****************************************************************************
!> \brief allocates and intitializes a fist_env
!> \param fist_env the object to create
!> \param para_env the parallel environement for the qs_env
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fist_env_create(fist_env,para_env,error)
    TYPE(fist_environment_type), POINTER     :: fist_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(fist_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL init_fist_env(fist_env,para_env=para_env, error=error)
    END IF
  END SUBROUTINE fist_env_create

! *****************************************************************************
!> \brief retains the given fist_env (see doc/ReferenceCounting.html)
!> \param fist_env the object to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fist_env_retain(fist_env,error)
    TYPE(fist_environment_type), POINTER     :: fist_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(fist_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(fist_env%ref_count>0,cp_failure_level,routineP,error,failure)
       fist_env%ref_count=fist_env%ref_count+1
    END IF
  END SUBROUTINE fist_env_retain

! *****************************************************************************
!> \brief releases the given fist_env (see doc/ReferenceCounting.html)
!> \param fist_env the object to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fist_env_release(fist_env,error)
    TYPE(fist_environment_type), POINTER     :: fist_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(fist_env)) THEN
       CPPrecondition(fist_env%ref_count>0,cp_failure_level,routineP,error,failure)
       fist_env%ref_count=fist_env%ref_count-1
       IF (fist_env%ref_count<1) THEN
          CALL qmmm_env_mm_release(fist_env%qmmm_env, error=error)
          CALL cell_release(fist_env%cell_ref,error=error)
          CALL ewald_pw_release(fist_env%ewald_pw,error=error)
          CALL ewald_env_release(fist_env%ewald_env,error=error)
          CALL cp_para_env_release(fist_env%para_env,error=error)
          CALL deallocate_fist_energy(fist_env%thermo)

          CALL fist_nonbond_env_release ( fist_env % fist_nonbond_env, error = error )
          CALL cp_subsys_release(fist_env%subsys,error=error)
          CALL section_vals_release(fist_env%input,error=error)
          CALL exclusion_release(fist_env%exclusions,error=error)

          DEALLOCATE(fist_env, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(fist_env)
  END SUBROUTINE fist_env_release

END MODULE fist_environment_types

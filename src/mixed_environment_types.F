!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/mixed_environment_types [1.0] *
!!
!!   NAME
!!     mixed_environment_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     fschiff SEPT-11-06
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE mixed_environment_types
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_create,&
                                             cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_retain,&
                                             cp_subsys_set,&
                                             cp_subsystem_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE mixed_energy_types,               ONLY: deallocate_mixed_energy,&
                                             mixed_energy_type
  USE input_section_types,             ONLY: section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
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
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
!


  TYPE mixed_environment_type
    PRIVATE
    INTEGER :: id_nr, ref_count
    TYPE ( cell_type ), POINTER                    :: cell
    TYPE ( cell_type ), POINTER                    :: cell_ref
    TYPE ( mixed_energy_type ), POINTER            :: mixed_energy
    TYPE ( cp_para_env_type ), POINTER             :: para_env
    TYPE ( cp_subsystem_type ), POINTER            :: subsys
    TYPE(section_vals_type), POINTER               :: input
 END TYPE mixed_environment_type

! *****************************************************************************

!!****s* mixed_environment_types/mixed_environment_p_type [1.0] *
!!
!!   NAME
!!     mixed_environment_p_type
!!
!!   FUNCTION
!!     to build arrays of pointers
!!
!!   ATTRIBUTES
!!     - mixed_env: the pointer to the mixed_env
!!
!!   AUTHOR
!!     fschiff
!!
!!   MODIFICATION HISTORY
!!     11/06
!!
!!   SOURCE
  !***************************************************************************
  TYPE mixed_environment_p_type
     TYPE(mixed_environment_type), POINTER :: mixed_env
  END TYPE mixed_environment_p_type
!!***
  !****************************************************************************

! *** Public data types ***

  PUBLIC :: mixed_environment_type,        &
            mixed_environment_p_type

! *** Public subroutines ***

  PUBLIC :: get_mixed_env,&
            init_mixed_env,&
            set_mixed_env,&
            mixed_env_create, &
            mixed_env_release, &
            mixed_env_retain

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mixed_environment_types'
  INTEGER, PRIVATE, SAVE :: last_mixed_env_id_nr=0

!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_mixed_env( mixed_env, atomic_kind_set, particle_set,  &
                           local_particles, local_molecules, molecule_kind_set, &
                           molecule_set, cell, cell_ref, &
                           mixed_energy, para_env, subsys,  &
                           input, error )

!   Purpose: Get the MIXED environment.

!   ***************************************************************************

    TYPE(mixed_environment_type), INTENT(IN)  :: mixed_env
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
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsys
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    TYPE(cp_error_type), INTENT(INOUT)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_mixed_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: particles

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY( atomic_kinds, particles, molecules_new, molecule_kinds_new )
!    CPPrecondition(fist_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mixed_env%subsys),cp_failure_level,routineP,error,failure)

    IF ( PRESENT ( input ) ) input => mixed_env%input
    IF ( PRESENT ( cell ) ) cell => mixed_env % cell
    IF ( PRESENT ( cell_ref ) ) cell_ref => mixed_env % cell_ref
    IF ( PRESENT ( mixed_energy ) ) mixed_energy => mixed_env % mixed_energy
    IF ( PRESENT ( para_env ) ) para_env => mixed_env % para_env
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

  END SUBROUTINE get_mixed_env
!
! *****************************************************************************

  SUBROUTINE init_mixed_env ( mixed_env, para_env, error )

!   Purpose: Initialise the MIXED environment.

!   ***************************************************************************

    TYPE(mixed_environment_type), INTENT(OUT) :: mixed_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

!   ---------------------------------------------------------------------------

    NULLIFY (mixed_env%input)
    NULLIFY (mixed_env%cell)
    NULLIFY (mixed_env%cell_ref)
    NULLIFY (mixed_env%mixed_energy)
    NULLIFY (mixed_env%subsys)
    CALL cp_subsys_create(mixed_env%subsys, para_env=para_env, error=error)
    CALL cp_para_env_retain(para_env,error=error)
    mixed_env%para_env => para_env
    mixed_env%ref_count=1
    last_mixed_env_id_nr=last_mixed_env_id_nr+1
    mixed_env%id_nr=last_mixed_env_id_nr

  END SUBROUTINE init_mixed_env

! *****************************************************************************

  SUBROUTINE set_mixed_env( mixed_env, atomic_kind_set, particle_set, &
                           local_particles, local_molecules, molecule_kind_set, &
                           molecule_set, cell, cell_ref,  &
                           mixed_energy, subsys, &
                           input, error )

!   Purpose: Set the MIXED environment.

    TYPE(mixed_environment_type), POINTER     :: mixed_env
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
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsys
     TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    TYPE(cp_error_type), INTENT(INOUT)  &
                                             :: error

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

!***************************************************************************
  END SUBROUTINE set_mixed_env
!***************************************************************************
!
!****f* mixed_environment_types/mixed_env_create [1.0] *
!
!   NAME
!     mixed_env_create
!
!   SYNOPSIS
!     Subroutine mixed_env_create(mixed_env, para_env, error)
!       Type(mixed_environment_type), Pointer:: mixed_env
!       Type(cp_para_env_type), Pointer:: para_env
!       Type(cp_error_type), Optional, Intent (INOUT):: error
!     End Subroutine mixed_env_create
!
!   FUNCTION
!     allocates and intitializes a mixed_env
!
!   ARGUMENTS
!     - mixed_env: the object to create
!     - para_env: the parallel environement for the qs_env
!     - error: variable to control error logging, stopping,...
!       see module cp_error_handling
!
!   AUTHOR
!     fschiff 11.06
!
!
!*** **********************************************************************
SUBROUTINE mixed_env_create(mixed_env,para_env,error)
    TYPE(mixed_environment_type), POINTER     :: mixed_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

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
!***************************************************************************

!****f* mixed_environment_types/mixed_env_retain [1.0] *
!
!   NAME
!     mixed_env_retain
!
!   SYNOPSIS
!     Subroutine mixed_env_retain(mixed_env, error)
!       Type(mixed_environment_type), Pointer:: mixed_env
!       Type(cp_error_type), Optional, Intent (INOUT):: error
!     End Subroutine mixed_env_retain
!
!   FUNCTION
!     retains the given mixed_env (see doc/ReferenceCounting.html)
!
!   NOTES
!     -
!
!   ARGUMENTS
!     - mixed_env: the object to retain
!     - error: variable to control error logging, stopping,...
!       see module cp_error_handling
!
!   AUTHOR
!     fschiff 11.06
!
!
!*** **********************************************************************
SUBROUTINE mixed_env_retain(mixed_env,error)
    TYPE(mixed_environment_type), POINTER     :: mixed_env
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mixed_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
!!
  CPPrecondition(ASSOCIATED(mixed_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(mixed_env%ref_count>0,cp_failure_level,routineP,error,failure)
     mixed_env%ref_count=mixed_env%ref_count+1
  END IF
END SUBROUTINE mixed_env_retain
!***************************************************************************

!****f* mixed_environment_types/mixed_env_release [1.0] *
!
!   NAME
!     mixed_env_release
!
!   SYNOPSIS
!     Subroutine mixed_env_release(mixed_env, error)
!       Type(mixed_environment_type), Pointer:: mixed_env
!       Type(cp_error_type), Optional, Intent (INOUT):: error
!     End Subroutine mixed_env_release
!
!   FUNCTION
!     releases the given mixed_env (see doc/ReferenceCounting.html)
!
!   NOTES
!     -
!
!   ARGUMENTS
!     - mixed_env: the object to release
!     - error: variable to control error logging, stopping,...
!       see module cp_error_handling
!
!   AUTHOR
!     fschiff 11.06
!
!
!*** **********************************************************************
SUBROUTINE mixed_env_release(mixed_env,error)
    TYPE(mixed_environment_type), POINTER     :: mixed_env
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mixed_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(mixed_env)) THEN
     CPPrecondition(mixed_env%ref_count>0,cp_failure_level,routineP,error,failure)
     mixed_env%ref_count=mixed_env%ref_count-1
     IF (mixed_env%ref_count<1) THEN
        CALL cell_release(mixed_env%cell,error=error)
        CALL cell_release(mixed_env%cell_ref,error=error)
        CALL cp_para_env_release(mixed_env%para_env,error=error)
        CALL deallocate_mixed_energy(mixed_env%mixed_energy)

        CALL cp_subsys_release(mixed_env%subsys,error=error)
        CALL section_vals_release(mixed_env%input,error=error)
        DEALLOCATE(mixed_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(mixed_env)
END SUBROUTINE mixed_env_release
!***************************************************************************
END MODULE mixed_environment_types

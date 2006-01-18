!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/tbmd_environment_types [1.0] *
!!
!!   NAME
!!     tbmd_environment_types
!!
!!   FUNCTION
!!     The tight-binding environment type, which contains everything needed 
!!     to actually do tight-binding.
!!
!!   NOTES
!!     -
!!     
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE tbmd_environment_types
  USE atomic_kind_list_types,                   ONLY: atomic_kind_list_type
  USE atomic_kind_types,                        ONLY: atomic_kind_type
  USE cp_subsystem_types,                       ONLY: cp_subsys_get,&
                                                      cp_subsys_release,&
                                                      cp_subsys_retain,&
                                                      cp_subsystem_type
  USE distribution_1d_types,                    ONLY: distribution_1d_type
  USE mol_kind_new_list_types,                  ONLY: mol_kind_new_list_type
  USE mol_new_list_types,                       ONLY: mol_new_list_type
  USE molecule_kind_types,                      ONLY: molecule_kind_type
  USE molecule_types_new,                       ONLY: molecule_type
  USE particle_list_types,                      ONLY: particle_list_type
  USE particle_types,                           ONLY: particle_type
  USE tbmd_energy_types,                        ONLY: tbmd_energy_type, &
                                                      tbmd_energy_release, &
                                                      tbmd_energy_clear, &
                                                      tbmd_energy_create
  USE tbmd_force_types,                         ONLY: tbmd_force_type, &
                                                      !tbmd_force_create, &
                                                      !tbmd_force_clear, &
                                                      tbmd_force_release
  USE tbmd_model_types,                         ONLY: tbmd_model_type, &
                                                      tbmd_model_release, &
                                                      tbmd_model_clear, &
                                                      tbmd_model_create
  USE virial_types,                             ONLY: virial_type, &
                                                      zero_virial, &
                                                      virial_release

#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='tbmd_environment_types'
  INTEGER, PRIVATE, save :: last_tbmd_id=0

  ! *** Public data types ***
  PUBLIC :: tbmd_environment_type, &
            tbmd_environment_p_type

  ! *** Public subroutines ***
  PUBLIC :: tbmd_env_retain, &
            tbmd_env_release, &
            tbmd_env_set,&
            tbmd_env_get, &
            tbmd_env_clear, &
            tbmd_env_create
!!***
!****************************************************************************

!!****s* tbmd_environment_types/tbmd_environment_type [1.0] *
!!
!!   NAME
!!     tbmd_environment_type
!!
!!   FUNCTION
!!     The TBMD Environment
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: unique identifier
!!     - ref_count: reference counter (see doc/ReferenceCounting.html)
!!     - tbmd_energy: The tight-binding energies
!!     - tbmd_force: The tight-binding forces
!!     - tbmd_model: The tight-binding parametrization model
!!     - virial: Dummy virial pointer
!!     - subsys: the particles, molecules,... of this environment
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
  TYPE tbmd_environment_type
    INTEGER                                          :: id_nr, &
                                                        ref_count
    TYPE(tbmd_energy_type), POINTER                  :: tbmd_energy
    TYPE(tbmd_force_type), POINTER                   :: tbmd_force
    TYPE(tbmd_model_type), POINTER                   :: tbmd_model
    TYPE(virial_type), POINTER                       :: virial
    TYPE(cp_subsystem_type), POINTER                 :: subsys

  END TYPE tbmd_environment_type
!!***
!****************************************************************************

!!****s* tbmd_environment_types/tbmd_environment_p_type [1.0] *
!!
!!   NAME
!!     tbmd_environment_p_type
!!
!!   FUNCTION
!!     Type which allows the use of arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - tbmd_env_ptr: Pointer to the actual tbmd environment
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!***************************************************************************
  TYPE tbmd_environment_p_type
    TYPE(tbmd_environment_type), POINTER :: tbmd_env_ptr
  END TYPE tbmd_environment_p_type
!!***
!****************************************************************************

CONTAINS

!!****f* tbmd_environment_types/tbmd_env_retain [1.0] *
!!
!!   NAME
!!     tbmd_env_retain
!!
!!   FUNCTION
!!     Retains a tbmd environment (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_env_retain(tbmd_env, error)

    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------
   
    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_env), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
      CPPrecondition(tbmd_env%ref_count>0, cp_failure_level, routineP, error, failure)
      tbmd_env%ref_count = tbmd_env%ref_count+1
    END IF
  END SUBROUTINE tbmd_env_retain
!***************************************************************************

!!****f* tbmd_environment_types/tbmd_env_release [1.0] *
!!
!!   NAME
!!     tbmd_env_release
!!
!!   FUNCTION
!!     Releases the given tbmd environment (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_env_release(tbmd_env, error)

    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    IF (ASSOCIATED(tbmd_env)) THEN
      CPPrecondition(tbmd_env%ref_count>0, cp_failure_level, routineP, error, failure)
      tbmd_env%ref_count = tbmd_env%ref_count-1
      IF (tbmd_env%ref_count<1) THEN
        IF (ASSOCIATED(tbmd_env%tbmd_energy)) THEN
          CALL tbmd_energy_release(tbmd_env%tbmd_energy, error=error)
        END IF
        IF (ASSOCIATED(tbmd_env%tbmd_force)) THEN
          CALL tbmd_force_release(tbmd_env%tbmd_force, error=error)
        END IF
        IF (ASSOCIATED(tbmd_env%tbmd_model)) THEN
          CALL tbmd_model_release(tbmd_env%tbmd_model, error=error)
        END IF
        IF (ASSOCIATED(tbmd_env%virial)) THEN
          CALL virial_release(tbmd_env%virial, error=error)
        END IF
        IF (ASSOCIATED(tbmd_env%subsys)) THEN
          CALL cp_subsys_release(tbmd_env%subsys, error=error)
        END IF
        DEALLOCATE(tbmd_env, stat=stat)
        CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
      END IF
    END IF
    NULLIFY(tbmd_env)
  END SUBROUTINE tbmd_env_release
!***************************************************************************

!!****f* tbmd_environment_types/tbmd_env_get [1.0] *
!!
!!   NAME
!!     tbmd_env_get
!!
!!   FUNCTION
!!     Returns various attributes of the tbmd environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_env: The enquired tbmd environment
!!     - id_nr: unique identifier
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!     - tbmd_energy: The tight-binding energies
!!     - tbmd_force: The tight-binding forces
!!     - tbmd_model: The tight-binding parametrization model
!!     - virial: Dummy virial pointer
!!     - subsys: the particles, molecules,... of this environment
!!
!!     For possible missing arguments see the attributes of 
!!     tbmd_environment_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_env_get(tbmd_env, tbmd_energy, tbmd_force, tbmd_model, &
                          virial, subsys, atomic_kind_set, particle_set , &
                          local_particles, molecule_kind_set, molecule_set, &
                          local_molecules, id_nr, error)

    TYPE(tbmd_environment_type), POINTER         :: tbmd_env
    TYPE(tbmd_energy_type), POINTER, OPTIONAL    :: tbmd_energy
    TYPE(tbmd_force_type), POINTER, OPTIONAL     :: tbmd_force
    TYPE(tbmd_model_type), POINTER, OPTIONAL     :: tbmd_model
    TYPE(virial_type), POINTER, OPTIONAL         :: virial
    TYPE(cp_subsystem_type), POINTER, OPTIONAL   :: subsys
    TYPE(atomic_kind_type), DIMENSION(:), &
         OPTIONAL, POINTER                       :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
         OPTIONAL, POINTER                       :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
         OPTIONAL, POINTER                       :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
         OPTIONAL, POINTER                       :: molecule_set
    TYPE(distribution_1d_type),OPTIONAL, POINTER :: local_particles,&
                                                    local_molecules
    INTEGER, INTENT(OUT), OPTIONAL               :: id_nr
    TYPE(cp_error_type), INTENT(INOUT), &
         OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                      :: failure
    TYPE (atomic_kind_list_type), POINTER        :: atomic_kinds
    TYPE (particle_list_type), POINTER           :: particles
    TYPE (mol_new_list_type), POINTER            :: molecules_new
    TYPE (mol_kind_new_list_type), POINTER       :: molecule_kinds_new

!   ------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(tbmd_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(tbmd_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) id_nr = tbmd_env%id_nr
      IF (PRESENT(tbmd_energy)) tbmd_energy => tbmd_env%tbmd_energy
      IF (PRESENT(tbmd_force)) tbmd_force => tbmd_env%tbmd_force
      IF (PRESENT(tbmd_model)) tbmd_model => tbmd_env%tbmd_model
      IF (PRESENT(virial)) virial => tbmd_env%virial
      IF (PRESENT(subsys)) subsys => tbmd_env%subsys
      CALL cp_subsys_get(tbmd_env%subsys,&
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
    END IF

  END SUBROUTINE tbmd_env_get
!***************************************************************************

!!****f* tbmd_environment_types/tbmd_env_set [1.0] *
!!
!!   NAME
!!     tbmd_env_set
!!
!!   FUNCTION
!!     Sets various attributes of the tbmd environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_env: The enquired tbmd environment
!!     - id_nr: unique identifier
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!     - tbmd_energy: The tight-binding energies
!!     - tbmd_force: The tight-binding forces
!!     - tbmd_model: The tight-binding parametrization model
!!     - virial: Dummy virial pointer
!!     - subsys: the particles, molecules,... of this environment
!!
!!     For possible missing arguments see the attributes of 
!!     tbmd_environment_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_env_set(tbmd_env, tbmd_energy, tbmd_force, tbmd_model, &
                          virial, subsys, id_nr, error)

    TYPE(tbmd_environment_type), POINTER        :: tbmd_env
    TYPE(tbmd_energy_type), POINTER, OPTIONAL   :: tbmd_energy
    TYPE(tbmd_force_type), POINTER, OPTIONAL    :: tbmd_force
    TYPE(tbmd_model_type), POINTER, OPTIONAL    :: tbmd_model
    TYPE(virial_type), POINTER, OPTIONAL        :: virial
    TYPE(cp_subsystem_type), POINTER, OPTIONAL  :: subsys
    INTEGER, INTENT(IN), OPTIONAL               :: id_nr
    TYPE(cp_error_type), INTENT(INOUT), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_env_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                     :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(tbmd_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) tbmd_env%id_nr = id_nr
      IF (PRESENT(tbmd_energy)) tbmd_env%tbmd_energy => tbmd_energy
      IF (PRESENT(tbmd_force)) tbmd_env%tbmd_force => tbmd_force
      IF (PRESENT(tbmd_model)) tbmd_env%tbmd_model => tbmd_model
      IF (PRESENT(virial)) tbmd_env%virial => virial
      IF (PRESENT(subsys)) THEN
       CALL cp_subsys_retain(subsys, error=error)
       CALL cp_subsys_release(tbmd_env%subsys, error=error)
       tbmd_env%subsys => subsys
      END IF
    END IF
  END SUBROUTINE tbmd_env_set
!***************************************************************************

!!****f* tbmd_environment_types/tbmd_env_clear [1.0] *
!!
!!   NAME
!!     tbmd_env_clear
!!
!!   FUNCTION
!!     Reinitializes the tbmd environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_env: The tbmd environment to be reinitialized
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of 
!!     tbmd_environment_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_env_clear(tbmd_env, error)

    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_env_clear', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN
      IF (ASSOCIATED(tbmd_env%tbmd_energy)) THEN
        CALL tbmd_energy_clear(tbmd_env%tbmd_energy, error=error)
      END IF
      IF (ASSOCIATED(tbmd_env%tbmd_force)) THEN
        !CALL tbmd_force_clear(tbmd_env%tbmd_force, error=error)
      END IF
      IF (ASSOCIATED(tbmd_env%tbmd_model)) THEN
        CALL tbmd_model_clear(tbmd_env%tbmd_model, error=error)
      END IF
      IF (ASSOCIATED(tbmd_env%virial)) THEN
        CALL zero_virial(tbmd_env%virial)
      END IF
      IF (ASSOCIATED(tbmd_env%subsys)) THEN
        CALL cp_subsys_release(tbmd_env%subsys, error=error)
      END IF
     ELSE
      CALL tbmd_env_create(tbmd_env, error=error)
    END IF
  END SUBROUTINE tbmd_env_clear
!***************************************************************************

!!****f* tbmd_environment_types/tbmd_env_create [1.0] *
!!
!!   NAME
!!     tbmd_env_create
!!
!!   FUNCTION
!!     Creates the tbmd environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - tbmd_env: The tbmd environment to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of 
!!     tbmd_environment_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_env_create(tbmd_env, error)

    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    ALLOCATE(tbmd_env, stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)

    CALL tbmd_energy_create(tbmd_env%tbmd_energy, error=error)
    !CALL tbmd_force_create(tbmd_env%tbmd_force, error=error)
    CALL tbmd_model_create(tbmd_env%tbmd_model, error=error)

    !CALL virial_create(tbmd_env%virial, error=error)
    !CALL cp_subsys_create(tbmd_env%subsys, para_env=para_env, error=error)
   
    IF (.NOT. failure) THEN
      tbmd_env%ref_count = 1
      last_tbmd_id=last_tbmd_id+1
      tbmd_env%id_nr = last_tbmd_id
    END IF
  END SUBROUTINE tbmd_env_create
!***************************************************************************

END MODULE tbmd_environment_types

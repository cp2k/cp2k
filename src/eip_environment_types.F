!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/eip_environment_types [1.0] *
!!
!!   NAME
!!     eip_environment_types
!!
!!   FUNCTION
!!     The environment for the empirical interatomic potential methods.
!!
!!   NOTES
!!     -
!!     
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE eip_environment_types
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type
  USE cp_subsystem_types,              ONLY: cp_subsystem_type, &
                                             cp_subsys_release, &
                                             cp_subsys_retain, &
                                             cp_subsys_get, &
                                             cp_subsys_set
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE input_section_types,             ONLY: section_vals_release,&
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
  USE virial_types,                    ONLY: virial_release,&
                                             virial_type,&
                                             zero_virial

#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eip_environment_types'
  INTEGER, PRIVATE, SAVE :: last_eip_id=0

  ! *** Public data types ***
  PUBLIC :: eip_environment_type, &
            eip_environment_p_type

  ! *** Public subroutines ***
  PUBLIC :: eip_env_retain, &
            eip_env_release, &
            eip_env_set,&
            eip_env_get, &
            eip_env_clear, &
            eip_env_create
!!***
!****************************************************************************

!!****s* eip_environment_types/eip_environment_type [1.0] *
!!
!!   NAME
!!     eip_environment_type
!!
!!   FUNCTION
!!     The empirical interatomic potential environment
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: unique identifier
!!     - ref_count: reference counter (see doc/ReferenceCounting.html)
!!     - eip_model: Specifies which EIP model is in use.
!!     - eip_kinetic_energy: The EIP kinetic energy
!!     - eip_potential_energy: The EIP potential energy
!!     - eip_energy: The total eip energy
!!     - eip_energy_var: Variance of the energy/atom
!!     - eip_forces: The final eip forces [eV/A]
!!     - coord_avg: The average coordination number
!!     - coord_var: The variance of the coordination number
!!     - count: Counts how often the empirical interatomic potential function
!!              is called. Don't ask why this is a real!
!!     - subsystem: the particles, molecules,... of this environment
!!     - eip_input: Pointer to the EIP input section
!!     - force_env_input: Pointer to the force_env input section
!!     - cell: The simulation cell
!!     - cell_ref: The reference simulation cell
!!     - use_ref_cell: Logical which indicates if reference 
!!                     simulation cell is used
!!     - virial: Dummy virial pointer
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
  TYPE eip_environment_type
    INTEGER                                          :: id_nr, &
                                                        ref_count, &
                                                        eip_model
    REAL(KIND=dp)                                    :: eip_energy, &
                                                        eip_kinetic_energy, &
                                                        eip_potential_energy, &
                                                        eip_energy_var
    REAL(KIND=dp), DIMENSION(:,:), POINTER           :: eip_forces
    REAL(KIND=dp)                                    :: coord_avg, &
                                                        coord_var, &
                                                        count
    TYPE(cp_subsystem_type), POINTER                 :: subsystem
    TYPE(section_vals_type), POINTER                 :: eip_input, &
                                                        force_env_input
    TYPE(cell_type), POINTER                         :: cell
    TYPE(cell_type), POINTER                         :: cell_ref
    LOGICAL                                          :: use_ref_cell
    TYPE(virial_type), POINTER                       :: virial

  END TYPE eip_environment_type
!!***
!****************************************************************************

!!****s* eip_environment_types/eip_environment_p_type [1.0] *
!!
!!   NAME
!!     eip_environment_p_type
!!
!!   FUNCTION
!!     Type which allows the use of arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - eip_env_ptr: Pointer to the actual eip environment
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!   SOURCE
!***************************************************************************
  TYPE eip_environment_p_type
    TYPE(eip_environment_type), POINTER :: eip_env_ptr
  END TYPE eip_environment_p_type
!!***
!****************************************************************************

CONTAINS

!!****f* eip_environment_types/eip_env_retain [1.0] *
!!
!!   NAME
!!     eip_env_retain
!!
!!   FUNCTION
!!     Retains a eip environment (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - eip_env: The eip environment to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_env_retain(eip_env, error)

    TYPE(eip_environment_type), POINTER      :: eip_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(eip_env), cp_failure_level, routineP, error, failure)
    IF (.NOT. failure) THEN
      CPPrecondition(eip_env%ref_count>0, cp_failure_level, routineP, error, failure)
      eip_env%ref_count = eip_env%ref_count+1
    END IF
  END SUBROUTINE eip_env_retain
!***************************************************************************

!!****f* eip_environment_types/eip_env_release [1.0] *
!!
!!   NAME
!!     eip_env_release
!!
!!   FUNCTION
!!     Releases the given eip environment (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - eip_env: The eip environment to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_env_release(eip_env, error)

    TYPE(eip_environment_type), POINTER      :: eip_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    IF (ASSOCIATED(eip_env)) THEN
      CPPrecondition(eip_env%ref_count>0, cp_failure_level, routineP, error, failure)
      eip_env%ref_count = eip_env%ref_count-1
      IF (eip_env%ref_count<1) THEN
        IF (ASSOCIATED(eip_env%eip_forces)) THEN
          DEALLOCATE(eip_env%eip_forces, stat=stat)
          CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
        END IF
        IF (ASSOCIATED(eip_env%subsystem)) THEN
          CALL cp_subsys_release(eip_env%subsystem, error=error)
        END IF
        IF (ASSOCIATED(eip_env%subsystem)) THEN
          CALL cp_subsys_release(eip_env%subsystem, error=error)
        END IF
        !IF (ASSOCIATED(eip_env%eip_input)) THEN
        !  CALL section_vals_release(eip_env%eip_input, error=error)
        !END IF
        !IF (ASSOCIATED(eip_env%force_env_input)) THEN
        !  CALL section_vals_release(eip_env%force_env_input, error=error)
        !END IF
        IF (ASSOCIATED(eip_env%cell)) THEN
          CALL cell_release(eip_env%cell, error=error)
        END IF
        IF (ASSOCIATED(eip_env%cell_ref)) THEN
          CALL cell_release(eip_env%cell_ref, error=error)
        END IF
        IF (ASSOCIATED(eip_env%virial)) THEN
          CALL virial_release(eip_env%virial, error=error)
        END IF
        DEALLOCATE(eip_env, stat=stat)
        CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
      END IF
    END IF
    NULLIFY(eip_env)
  END SUBROUTINE eip_env_release
!***************************************************************************

!!****f* eip_environment_types/eip_env_get [1.0] *
!!
!!   NAME
!!     eip_env_get
!!
!!   FUNCTION
!!     Returns various attributes of the eip environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - eip_env: The enquired eip environment
!!     - id_nr: unique identifier
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!     - eip_model: Specifies which EIP model is in use.
!!     - eip_kinetic_energy: The EIP kinetic energy
!!     - eip_potential_energy: The EIP potential energy
!!     - eip_energy: The total eip energy
!!     - eip_energy_var: Variance of the energy/atom
!!     - eip_forces: The final eip forces [eV/A]
!!     - coord_avg: The average coordination number
!!     - coord_var: The variance of the coordination number
!!     - count: Counts how often the empirical interatomic potential function
!!              is called. Don't ask why this is a real!
!!     - subsystem: the particles, molecules,... of this environment
!!     - atomic_kind_set: The set of all atomic kinds involved
!!     - particle_set: The set of all particles
!!     - local_particles: All particles on this particular node
!!     - molecule_kind_set: The set of all different molecule kinds involved
!!     - molecule_set: The set of all molecules
!!     - local_molecules: All molecules on this particular node
!!     - eip_input: the pointer to the EIP input section
!!     - force_env_input: Pointer to the force_env input section
!!     - cell: The simulation cell
!!     - cell_ref: The reference simulation cell
!!     - use_ref_cell: Logical which indicates if reference 
!!                     simulation cell is used
!!     - virial: Dummy virial pointer
!!
!!     For possible missing arguments see the attributes of 
!!     eip_environment_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_env_get(eip_env, id_nr, eip_model, eip_energy, eip_energy_var, &
                         eip_forces, coord_avg, coord_var, count, subsystem, &
                         atomic_kind_set, particle_set , local_particles, &
                         molecule_kind_set, molecule_set, local_molecules, &
                         eip_input, force_env_input, cell, cell_ref, &
                         use_ref_cell, eip_kinetic_energy, eip_potential_energy, &
                         virial, error)

    TYPE(eip_environment_type), POINTER      :: eip_env
    INTEGER, INTENT(OUT), OPTIONAL           :: id_nr
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    INTEGER, INTENT(OUT), OPTIONAL           :: eip_model
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: eip_energy, &
                                                eip_kinetic_energy, &
                                                eip_potential_energy, &
                                                eip_energy_var
    REAL(KIND=dp), DIMENSION(:,:), OPTIONAL, &
      POINTER                                :: eip_forces
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: coord_avg, &
                                                coord_var, &
                                                count
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsystem
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_set
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles, &
                                                local_molecules
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: eip_input, &
                                                force_env_input
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    LOGICAL, INTENT(OUT), OPTIONAL           :: use_ref_cell
    TYPE(virial_type), OPTIONAL, POINTER     :: virial

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: particles

!   ------------------------------------------------------------------------

    failure = .FALSE.

    NULLIFY(atomic_kinds, particles, molecules_new, molecule_kinds_new)

    CPPrecondition(ASSOCIATED(eip_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(eip_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) id_nr = eip_env%id_nr
      IF (PRESENT(eip_model)) eip_model = eip_env%eip_model
      IF (PRESENT(eip_kinetic_energy)) THEN
        eip_kinetic_energy = eip_env%eip_kinetic_energy
      END IF
      IF (PRESENT(eip_potential_energy)) THEN
        eip_potential_energy = eip_env%eip_potential_energy
      END IF
      IF (PRESENT(eip_energy)) eip_energy = eip_env%eip_energy
      IF (PRESENT(eip_energy_var)) eip_energy_var = eip_env%eip_energy_var
      IF (PRESENT(eip_forces)) eip_forces = eip_env%eip_forces
      IF (PRESENT(coord_avg)) coord_avg = eip_env%coord_avg
      IF (PRESENT(coord_var)) coord_var = eip_env%coord_var
      IF (PRESENT(count)) count = eip_env%count
      IF (PRESENT(subsystem)) subsystem => eip_env%subsystem
      CALL cp_subsys_get(eip_env%subsystem, &
                         atomic_kinds=atomic_kinds, &
                         particles=particles, &
                         molecule_kinds_new=molecule_kinds_new, &
                         molecules_new=molecules_new, &
                         local_molecules_new=local_molecules, &
                         local_particles=local_particles, &
                         error=error)
      IF (PRESENT(atomic_kind_set)) atomic_kind_set => atomic_kinds%els
      IF (PRESENT(particle_set)) particle_set => particles%els
      IF (PRESENT(molecule_kind_set)) molecule_kind_set => molecule_kinds_new%els
      IF (PRESENT(molecule_set)) molecule_set => molecules_new%els

      IF (PRESENT(eip_input)) eip_input => eip_env%eip_input
      IF (PRESENT(force_env_input)) force_env_input => eip_env%force_env_input
      IF (PRESENT(cell)) cell => eip_env%cell
      IF (PRESENT(cell_ref)) cell_ref => eip_env%cell_ref
      IF (PRESENT(use_ref_cell)) use_ref_cell = eip_env%use_ref_cell
      IF (PRESENT(virial)) virial => eip_env%virial
    END IF

  END SUBROUTINE eip_env_get
!***************************************************************************

!!****f* eip_environment_types/eip_env_set [1.0] *
!!
!!   NAME
!!     eip_env_set
!!
!!   FUNCTION
!!     Sets various attributes of the eip environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - eip_env: The enquired eip environment
!!     - id_nr: unique identifier
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!     - eip_model: Specifies which EIP model is in use
!!     - eip_kinetic_energy: The EIP kinetic energy
!!     - eip_potential_energy: The EIP potential energy
!!     - eip_energy: The total eip energy
!!     - eip_energy_var: Variance of the energy/atom
!!     - eip_forces: The final eip forces [eV/A]
!!     - coord_avg: The average coordination number
!!     - coord_var: The variance of the coordination number
!!     - count: Counts how often the empirical interatomic potential function
!!              is called. Don't ask why this is a real!
!!     - subsystem: the particles, molecules,... of this environment
!!     - atomic_kind_set: The set of all atomic kinds involved
!!     - particle_set: The set of all particles
!!     - local_particles: All particles on this particular node
!!     - molecule_kind_set: The set of all different molecule kinds involved
!!     - molecule_set: The set of all molecules
!!     - local_molecules: All molecules on this particular node
!!     - eip_input: the pointer to the EIP input section
!!     - force_env_input: Pointer to the force_env input section
!!     - cell: The simulation cell
!!     - cell_ref: The reference simulation cell
!!     - use_ref_cell: Logical which indicates if reference 
!!                     simulation cell is used
!!     - virial: Dummy virial pointer
!!
!!     For possible missing arguments see the attributes of 
!!     eip_environment_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_env_set(eip_env, id_nr, eip_model, eip_energy, eip_energy_var, &
                         eip_forces, coord_avg, coord_var, count, subsystem, &
                         atomic_kind_set, particle_set , local_particles, &
                         molecule_kind_set, molecule_set, local_molecules, &
                         eip_input, force_env_input, cell, cell_ref, &
                         use_ref_cell, eip_kinetic_energy, eip_potential_energy, &
                         virial, error)

    TYPE(eip_environment_type), POINTER      :: eip_env
    INTEGER, INTENT(IN), OPTIONAL            :: id_nr
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    INTEGER, INTENT(IN), OPTIONAL            :: eip_model
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: eip_energy, &
                                                eip_kinetic_energy, &
                                                eip_potential_energy, &
                                                eip_energy_var
    REAL(KIND=dp), DIMENSION(:,:), OPTIONAL, &
      POINTER                                :: eip_forces
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: coord_avg, &
                                                coord_var, &
                                                count
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsystem
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_set
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles, &
                                                local_molecules
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: eip_input, &
                                                force_env_input
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    LOGICAL, INTENT(IN), OPTIONAL            :: use_ref_cell
    TYPE(virial_type), OPTIONAL, POINTER     :: virial

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_env_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: particles

!   ------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(eip_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(eip_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      IF (PRESENT(id_nr)) eip_env%id_nr = id_nr
      IF (PRESENT(eip_model)) eip_env%eip_model = eip_model
      IF (PRESENT(eip_kinetic_energy)) THEN
        eip_env%eip_kinetic_energy = eip_kinetic_energy
      END IF
      IF (PRESENT(eip_potential_energy)) THEN
        eip_env%eip_potential_energy = eip_potential_energy
      END IF
      IF (PRESENT(eip_energy)) eip_env%eip_energy = eip_energy
      IF (PRESENT(eip_energy_var)) eip_env%eip_energy_var = eip_energy_var
      IF (PRESENT(eip_forces)) eip_env%eip_forces = eip_forces
      IF (PRESENT(coord_avg)) eip_env%coord_avg = coord_avg
      IF (PRESENT(coord_var)) eip_env%coord_var = coord_var
      IF (PRESENT(count)) eip_env%count = count
      IF (PRESENT(subsystem)) THEN
        CALL cp_subsys_retain(subsystem, error=error)
        CALL cp_subsys_release(eip_env%subsystem, error=error)
        eip_env%subsystem => subsystem
      END IF
      IF (PRESENT(atomic_kind_set)) THEN
        CALL atomic_kind_list_create(atomic_kinds, &
                                     els_ptr=atomic_kind_set, &
                                     error=error)
        CALL cp_subsys_set(eip_env%subsystem, &
                           atomic_kinds=atomic_kinds, &
                           error=error)
        CALL atomic_kind_list_release(atomic_kinds, error=error)
      END IF
      IF (PRESENT(particle_set)) THEN
        CALL particle_list_create(particles, &
                                  els_ptr=particle_set, &
                                  error=error)
        CALL cp_subsys_set(eip_env%subsystem, &
                           particles=particles, &
                           error=error)
        CALL particle_list_release(particles, error=error)
      END IF
      IF (PRESENT(molecule_kind_set)) THEN
        CALL mol_kind_new_list_create(molecule_kinds_new, &
                                      els_ptr=molecule_kind_set, &
                                      error=error)
        CALL cp_subsys_set(eip_env%subsystem, &
                           molecule_kinds_new=molecule_kinds_new, &
                           error=error)
        CALL mol_kind_new_list_release(molecule_kinds_new, error=error)
      END IF
      IF (PRESENT(molecule_set)) THEN
        CALL mol_new_list_create(molecules_new, &
                                 els_ptr=molecule_set, &
                                 error=error)
        CALL cp_subsys_set(eip_env%subsystem, &
                           molecules_new=molecules_new, &
                           error=error)
        CALL mol_new_list_release(molecules_new,error=error)
      END IF
      IF (PRESENT(local_particles)) THEN
        CALL cp_subsys_set(eip_env%subsystem, &
                           local_particles=local_particles, &
                           error=error)
      END IF
      IF (PRESENT(local_molecules)) THEN
        CALL cp_subsys_set(eip_env%subsystem, &
                           local_molecules_new=local_molecules, &
                           error=error)
      END IF

      IF (PRESENT(eip_input)) eip_env%eip_input => eip_input
      IF (PRESENT(force_env_input)) THEN
        eip_env%force_env_input => force_env_input
      END IF
      IF (PRESENT(cell)) THEN
        CALL cell_retain(cell, error=error)
        CALL cell_release(eip_env%cell,error=error)
        eip_env%cell => cell
      END IF
      IF (PRESENT(cell_ref)) THEN
        CALL cell_retain(cell_ref, error=error)
        CALL cell_release(eip_env%cell_ref,error=error)
        eip_env%cell_ref => cell_ref
      END IF
      IF (PRESENT(use_ref_cell)) eip_env%use_ref_cell = use_ref_cell
      IF (PRESENT(virial)) eip_env%virial => virial
    END IF
  END SUBROUTINE eip_env_set
!***************************************************************************

!!****f* eip_environment_types/eip_env_clear [1.0] *
!!
!!   NAME
!!     eip_env_clear
!!
!!   FUNCTION
!!     Reinitializes the eip environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - eip_env: The eip environment to be reinitialized
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of 
!!     eip_environment_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_env_clear(eip_env, error)

    TYPE(eip_environment_type), POINTER      :: eip_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_env_clear', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(eip_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN
      eip_env%eip_model            = 0
      eip_env%eip_kinetic_energy   = 0.0_dp
      eip_env%eip_potential_energy = 0.0_dp
      eip_env%eip_energy           = 0.0_dp
      eip_env%eip_energy_var       = 0.0_dp
      eip_env%eip_forces(:,:)      = 0.0_dp
      eip_env%coord_avg            = 0.0_dp
      eip_env%coord_var            = 0.0_dp
      eip_env%count                = 0.0_dp
      IF (ASSOCIATED(eip_env%subsystem)) THEN
        CALL cp_subsys_release(eip_env%subsystem, error=error)
      END IF
      IF (ASSOCIATED(eip_env%eip_input)) THEN
        CALL section_vals_release(eip_env%eip_input, error=error)
      END IF
      IF (ASSOCIATED(eip_env%force_env_input)) THEN
        CALL section_vals_release(eip_env%force_env_input, error=error)
      END IF
      IF (ASSOCIATED(eip_env%cell)) THEN
        CALL cell_release(eip_env%cell, error=error)
      END IF
      IF (ASSOCIATED(eip_env%cell_ref)) THEN
        CALL cell_release(eip_env%cell_ref, error=error)
      END IF
      IF (ASSOCIATED(eip_env%virial)) THEN
        CALL zero_virial(eip_env%virial)
      END IF
    ELSE
      CALL eip_env_create(eip_env, error=error)
    END IF
  END SUBROUTINE eip_env_clear
!***************************************************************************

!!****f* eip_environment_types/eip_env_create [1.0] *
!!
!!   NAME
!!     eip_env_create
!!
!!   FUNCTION
!!     Creates the eip environment
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - eip_env: The eip environment to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!     For possible missing arguments see the attributes of 
!!     eip_environment_type
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_env_create(eip_env, error)

    TYPE(eip_environment_type), POINTER      :: eip_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    ALLOCATE(eip_env, stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)

    !ALLOCATE(eip_env%eip_forces(3,natom), stat=stat)
    !CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)

    NULLIFY(eip_env%eip_forces)
    NULLIFY(eip_env%subsystem)
    NULLIFY(eip_env%eip_input)
    NULLIFY(eip_env%force_env_input)
    NULLIFY(eip_env%cell)
    NULLIFY(eip_env%cell_ref)
    NULLIFY(eip_env%virial)

    IF (.NOT. failure) THEN
      eip_env%ref_count = 1
      last_eip_id=last_eip_id+1
      eip_env%id_nr = last_eip_id
      eip_env%use_ref_cell = .FALSE.
      CALL eip_env_clear(eip_env, error=error)
    END IF
  END SUBROUTINE eip_env_create
!***************************************************************************

END MODULE eip_environment_types

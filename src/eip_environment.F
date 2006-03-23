!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/eip_environment [1.0] *
!!
!!   NAME
!!     eip_environment
!!
!!   FUNCTION
!!     Methods and functions on the EIP environment
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
MODULE eip_environment
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             write_atomic_kind_set
  USE cell_types,                      ONLY: cell_release,&
                                             cell_type,&
                                             get_cell,&
                                             read_cell,&
                                             write_cell
  USE cp_subsystem_methods,            ONLY: cp_subsys_read_colvar
  USE cp_subsystem_types,              ONLY: cp_subsys_create,&
                                             cp_subsys_release,&
                                             cp_subsys_set,&
                                             cp_subsystem_type
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d
  USE eip_environment_types,           ONLY: eip_env_set,&
                                             eip_environment_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_create,&
                                             mol_kind_new_list_release,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_create,&
                                             mol_new_list_release,&
                                             mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type,&
                                             write_molecule_kind_set
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_distances,&
                                             write_qs_particle_coordinates
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: topology_control
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eip_environment'

! *** Public subroutines ***

  PUBLIC :: eip_init

!!***
! *****************************************************************************

CONTAINS

!!****f* eip_environment/eip_init [1.0] *
!!
!!   NAME
!!     eip_init
!!
!!   FUNCTION
!!     Initialize the eip environment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - eip_env: The eip environment to retain
!!     - globenv: The global environment
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
  SUBROUTINE eip_init(eip_env, globenv, error)
    TYPE(eip_environment_type), POINTER      :: eip_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, use_ref_cell
    REAL(KIND=dp), DIMENSION(3)              :: abc
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_subsystem_type), POINTER         :: subsystem
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: cell_section, colvar_section, &
                                                eip_section, force_env_section

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(eip_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN
      ! nullifying  pointers
      NULLIFY(force_env_section, cell_section, colvar_section, eip_section, &
              cell, cell_ref, subsystem, molecule_kind_set, molecule_set, &
              particle_set, atomic_kind_set, particles, molecules, &
              molecule_kinds, atomic_kinds)

      force_env_section => section_vals_get_subs_vals(globenv%input_file, &
                                                      "FORCE_EVAL", error=error)
      cell_section => section_vals_get_subs_vals(force_env_section,&
                                                 "SUBSYS%CELL", error=error)
      colvar_section => section_vals_get_subs_vals(force_env_section, &
                                                   "SUBSYS%COLVAR", error=error)
      eip_section => section_vals_get_subs_vals(force_env_section, &
                                                 "EIP", error=error)

      CALL eip_env_set(eip_env=eip_env, eip_input=eip_section, &
                       force_env_input=force_env_section, error=error)

      CALL read_cell(cell=cell, cell_ref=cell_ref, globenv=globenv, &
                     use_ref_cell=use_ref_cell, cell_section=cell_section, &
                     error=error)
      CALL get_cell(cell=cell, abc=abc)
      CALL write_cell(cell=cell, force_env_section=force_env_section,error=error)

      CALL cp_subsys_create(subsystem, para_env=globenv%para_env, error=error)
      CALL cp_subsys_read_colvar(subsys=subsystem, colvar_section=colvar_section, &
                                 globenv=globenv, error=error)

      CALL topology_control(atomic_kind_set, particle_set, &
                            molecule_kind_set, molecule_set, &
                            subsystem%colvar_p, globenv, error=error)

      !CALL read_atomic_kind_set(atomic_kind_set, globenv)

      CALL particle_list_create(particles, els_ptr=particle_set, error=error)
      CALL atomic_kind_list_create(atomic_kinds, els_ptr=atomic_kind_set, error=error)
      CALL mol_new_list_create(molecules, els_ptr=molecule_set, error=error)
      CALL mol_kind_new_list_create(molecule_kinds, els_ptr=molecule_kind_set, error=error)
      CALL cp_subsys_set(subsystem, particles=particles, atomic_kinds=atomic_kinds, &
                         molecules_new=molecules, molecule_kinds_new=molecule_kinds, &
                         error=error)
      CALL particle_list_release(particles, error=error)
      CALL atomic_kind_list_release(atomic_kinds, error=error)
      CALL mol_new_list_release(molecules, error=error)
      CALL mol_kind_new_list_release(molecule_kinds, error=error)

      CALL eip_init_subsys(eip_env=eip_env, subsystem=subsystem, cell=cell, &
                           cell_ref=cell_ref, use_ref_cell=use_ref_cell, &
                           error=error)

      CALL cell_release(cell, error=error)
      CALL cell_release(cell_ref, error=error)
      CALL cp_subsys_release(subsystem, error=error)
    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE eip_init
!***************************************************************************

!!****f* eip_environment/eip_init_subsys [1.0] *
!!
!!   NAME
!!     eip_init_subsys
!!
!!   FUNCTION
!!     Initialize the eip environment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - eip_env: The eip environment of matter
!!     - subsystem: the subsystem
!!     - cell: Pointer to the actual simulation cell
!!     - cell_ref: Pointer to the reference cell, used e.g. in NPT simulations
!!     - use_ref_cell: Logical which indicates if cell_ref is in use
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
  SUBROUTINE eip_init_subsys(eip_env, subsystem, cell, cell_ref, &
                             use_ref_cell, error)
    TYPE(eip_environment_type), POINTER      :: eip_env
    TYPE(cp_subsystem_type), POINTER         :: subsystem
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    LOGICAL, INTENT(in)                      :: use_ref_cell
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_init_subsys', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, natom, stat
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(atomic_kind_set, molecule_kind_set, particle_set, molecule_set, &
            local_molecules, local_particles)

    failure = .FALSE.

    particle_set => subsystem%particles%els
    atomic_kind_set => subsystem%atomic_kinds%els
    molecule_kind_set => subsystem%molecule_kinds_new%els
    molecule_set => subsystem%molecules_new%els

!   *** Print the atomic kind set ***
    CALL write_atomic_kind_set(atomic_kind_set, eip_env%force_env_input, &
                               error=error)

!   *** Print the molecule kind set ***
    CALL write_molecule_kind_set(molecule_kind_set, eip_env%force_env_input, &
                                 error=error)

!   *** Print the atomic coordinates 
    CALL write_qs_particle_coordinates(particle_set, cell=cell, &
                                       force_env_section=eip_env%force_env_input, &
                                       error=error)
    CALL write_particle_distances(particle_set, cell=cell, &
                                  force_env_section=eip_env%force_env_input, &
                                  error=error)

!   *** Distribute molecules and atoms using the new data structures ***
    CALL distribute_molecules_1d(particle_kind_set=atomic_kind_set, &
                                 particle_set=particle_set, &
                                 local_particles=local_particles, &
                                 molecule_kind_set=molecule_kind_set, &
                                 molecule_set=molecule_set, &
                                 local_molecules=local_molecules, &
                                 force_env_section=eip_env%force_env_input, &
                                 error=error)

    natom = SIZE(particle_set)

    ALLOCATE(eip_env%eip_forces(3,natom), stat=stat)
    CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)

    eip_env%eip_forces(:,:) = 0.0_dp

    CALL eip_env_set(eip_env=eip_env, subsystem=subsystem, cell=cell, &
                     cell_ref=cell_ref, use_ref_cell=use_ref_cell, &
                     local_molecules=local_molecules, &
                     local_particles=local_particles, error=error)

    CALL distribution_1d_release(local_particles, error=error)
    CALL distribution_1d_release(local_molecules, error=error)

    CALL eip_init_model(eip_env=eip_env, error=error)

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE eip_init_subsys
!***************************************************************************

!!****f* eip_environment/eip_init_model [1.0] *
!!
!!   NAME
!!     eip_init_model
!!
!!   FUNCTION
!!     Initialize the empirical interatomic potnetial (force field) model 
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
  SUBROUTINE eip_init_model(eip_env, error)
    TYPE(eip_environment_type), POINTER      :: eip_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_init_model', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: eip_atomic_kind_name
    INTEGER                                  :: handle, i
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_ptr
    TYPE(section_vals_type), POINTER         :: eip_section

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(atomic_kind_set, atomic_kind_ptr, eip_section)

    eip_section => section_vals_get_subs_vals(eip_env%force_env_input, &
                                              "EIP", error=error)

    atomic_kind_set => eip_env%subsystem%atomic_kinds%els

    ! loop over all kinds
    DO i=1,SIZE(atomic_kind_set)
      atomic_kind_ptr => eip_env%subsystem%atomic_kinds%els(i)
      CALL get_atomic_kind(atomic_kind=atomic_kind_ptr, &
                           name=eip_atomic_kind_name)
      SELECT CASE(eip_atomic_kind_name)
        CASE("SI", "Si")
          CALL section_vals_val_get(section_vals=eip_section, &
                                    keyword_name="EIP-Model", &
                                    i_val=eip_env%eip_model, &
                                    error=error)
        CASE DEFAULT
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="EIP models for other elements" //&
               "than Si isn't implemented yet.", &
               error=error, error_level=cp_failure_level)
      END SELECT
    END DO

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE eip_init_model
!***************************************************************************

END MODULE eip_environment

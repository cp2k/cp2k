!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2003  CP2K developers group
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/trajana_methods [1.0] *
!!
!!   NAME
!!     trajana_methods
!!
!!   FUNCTION
!!     TRAJectory ANAlysis
!!
!!   AUTHOR
!!     Joost VandeVondele (11.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE trajana_methods
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_retain,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind, write_atomic_kind_set, &
                                             read_atomic_kind_set,deallocate_atomic_kind_set
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_check,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error
  USE cp_subsystem_types,               ONLY: cp_subsystem_p_type,&
                                             cp_subsystem_type,&
                                             cp_subsys_get
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d
  USE cntl_input,                      ONLY: read_cntl_section
  USE setup_input,                     ONLY: read_setup_section
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             dp,&
                                             dp_size
  USE input_types,                     ONLY: setup_parameters_type
  USE machine,                         ONLY: m_walltime
  USE md,                              ONLY: initialize_velocities,&
                                             mdio_parameters_type,&
                                             qs_md_parameters_type,&
                                             read_md_section,&
                                             simulation_parameters_type,&
                                             virial_type
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_create,&
                                             mol_kind_new_list_release,&
                                             mol_kind_new_list_retain,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_create,&
                                             mol_new_list_release,&
                                             mol_new_list_retain,&
                                             mol_new_list_type
  USE mol_struct_list_types,           ONLY: mol_struct_list_create,&
                                             mol_struct_list_release,&
                                             mol_struct_list_retain,&
                                             mol_struct_list_type
  USE molecule_types,                  ONLY: molecule_structure_type,&
                                             molecule_type

  USE molecule_types_new,              ONLY: molecule_type_new => molecule_type,deallocate_molecule_set
  USE molecule_kind_types,             ONLY: molecule_kind_type,&
                                             write_molecule_kind_set, &
                                             deallocate_molecule_kind_set
  USE output_utilities,                ONLY: print_warning
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_retain,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type, &
                                             write_qs_particle_coordinates, &
                                             deallocate_particle_set
  USE physcon,                         ONLY: kelvin
  USE simulation_cell,                 ONLY: cell_type, &
                                             read_cell, &
                                             write_cell
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: topology_control
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "trajana"

  PUBLIC :: trajana

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE trajana ( globenv, error )

! ***************************************************************************

    TYPE(global_environment_type), POINTER :: globenv
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN = "trajana",&
         routineP=module_name//':'//routineN

!   *** Local variables ***
    INTEGER                                   :: handle
    LOGICAL                                   :: failure
    TYPE(cell_type), POINTER                  :: cell

    TYPE(atomic_kind_type), DIMENSION(:), POINTER   :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), POINTER :: molecule_kind_set
    TYPE(molecule_type_new), DIMENSION(:), POINTER      :: molecule_set    
    TYPE(distribution_1d_type), POINTER :: local_molecules,&
                                           local_particles
    LOGICAL  :: found_frame

!   ---------------------------------------------------------------------------
    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(cell,atomic_kind_set,particle_set, molecule_kind_set,molecule_set)
    NULLIFY(local_molecules, local_particles)

    CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)

    ! CELL
    CALL read_cell(cell,globenv) 
    CALL write_cell(cell,globenv)

    ! get a list of particles / molecules set up 
    CALL topology_control(atomic_kind_set,particle_set,&
                         molecule_kind_set,molecule_set,&
                         globenv)
    CALL read_atomic_kind_set(atomic_kind_set,globenv)
    CALL write_atomic_kind_set(atomic_kind_set,globenv)
    CALL write_qs_particle_coordinates(particle_set,cell,globenv)

    ! set up a parallel distribution 
    CALL distribute_molecules_1d(particle_kind_set=atomic_kind_set,&
                                 particle_set=particle_set,&
                                 local_particles=local_particles,&
                                 molecule_kind_set=molecule_kind_set,&
                                 molecule_set=molecule_set,&
                                 local_molecules=local_molecules,&
                                 error=error)

    ! loop over the frames
    DO
      CALL get_next_frame(particle_set,found_frame)
      IF (found_frame) THEN
          
      ELSE
          EXIT
      ENDIF
    END DO

    ! release stuff
    ! CALL deallocate_molecule_kind_set(molecule_kind_set) seems buggy right now
    CALL deallocate_molecule_set(molecule_set)
    CALL deallocate_particle_set(particle_set)
    CALL deallocate_atomic_kind_set(atomic_kind_set)

    CALL distribution_1d_release(local_molecules)
    CALL distribution_1d_release(local_particles)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE trajana

  SUBROUTINE get_next_frame(particle_set,found_frame)
    TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
    LOGICAL                                         :: found_frame
    
    ! unhappily no way to read our own trajectory
    found_frame = .FALSE.
  END SUBROUTINE get_next_frame

END MODULE trajana_methods

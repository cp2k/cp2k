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
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             deallocate_atomic_kind_set,&
                                             read_atomic_kind_set,&
                                             write_atomic_kind_set
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE molecule_kind_types,             ONLY: deallocate_molecule_kind_set,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: deallocate_molecule_set,&
                                             molecule_type
  USE particle_types,                  ONLY: deallocate_particle_set,&
                                             particle_type,&
                                             write_qs_particle_coordinates
  USE simulation_cell,                 ONLY: cell_type,&
                                             read_cell,&
                                             write_cell
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

    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "trajana", &
      routineP = module_name//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, found_frame
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

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
    CALL deallocate_molecule_kind_set(molecule_kind_set)
    CALL deallocate_molecule_set(molecule_set)
    CALL deallocate_particle_set(particle_set)
    CALL deallocate_atomic_kind_set(atomic_kind_set)

    CALL distribution_1d_release(local_molecules)
    CALL distribution_1d_release(local_particles)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE trajana

  SUBROUTINE get_next_frame(particle_set,found_frame)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    LOGICAL                                  :: found_frame

! unhappily no way to read our own trajectory

    found_frame = .FALSE.
  END SUBROUTINE get_next_frame

END MODULE trajana_methods

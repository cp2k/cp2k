!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief initialize mixed environment
!> \author fschiff
! *****************************************************************************
MODULE mixed_environment
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_methods,                    ONLY: read_cell,&
                                             write_cell
  USE cell_types,                      ONLY: cell_release,&
                                             cell_type,&
                                             get_cell
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_methods,               ONLY: cp_subsys_create
  USE cp_subsys_types,                 ONLY: cp_subsys_release,&
                                             cp_subsys_set,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mixed_energy_types,              ONLY: allocate_mixed_energy,&
                                             mixed_energy_type
  USE mixed_environment_types,         ONLY: mixed_environment_type,&
                                             set_mixed_env
  USE molecule_kind_types,             ONLY: molecule_kind_type,&
                                             write_molecule_kind_set
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mixed_environment'
  PUBLIC :: mixed_init

CONTAINS

! *****************************************************************************
!> \brief reads the input and database file for mixed
!> \param mixed_env ...
!> \param root_section ...
!> \param para_env ...
!> \param force_env_section ...
!> \param use_motion_section ...
!> \par Used By
!>      mixed_main
!> \author fschiff
! *****************************************************************************
  SUBROUTINE mixed_init ( mixed_env, root_section,para_env, force_env_section,&
       use_motion_section)

    TYPE(mixed_environment_type), POINTER    :: mixed_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: force_env_section
    LOGICAL, INTENT(IN)                      :: use_motion_section

    CHARACTER(len=*), PARAMETER :: routineN = 'mixed_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, use_ref_cell
    REAL(KIND=dp), DIMENSION(3)              :: abc
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(section_vals_type), POINTER         :: cell_section, subsys_section

    CALL timeset(routineN,handle)
    failure = .FALSE.

    NULLIFY (subsys, cell, cell_ref)
    NULLIFY ( cell_section)

    subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS")
    cell_section => section_vals_get_subs_vals(subsys_section,"CELL")

    CALL set_mixed_env(mixed_env,input=force_env_section)
    CALL cp_subsys_create(subsys, para_env, root_section,&
                         force_env_section=force_env_section,&
                         use_motion_section=use_motion_section)

    CALL read_cell( cell, cell_ref, use_ref_cell = use_ref_cell, &
         cell_section=cell_section, para_env=para_env)
    CALL get_cell (cell, abc=abc)

    ! Print the cell parameters ***
    CALL write_cell ( cell, subsys_section, cell_ref)

    CALL mixed_init_subsys ( mixed_env, subsys, cell, cell_ref,  &
         force_env_section, subsys_section)

    CALL cell_release(cell)
    CALL cell_release(cell_ref)
    CALL cp_subsys_release(subsys)

    CALL timestop(handle)

  END SUBROUTINE mixed_init

! *****************************************************************************
!> \brief   Read the input and the database files for the setup of the
!>          mixed environment.
!> \param mixed_env ...
!> \param subsys ...
!> \param cell ...
!> \param cell_ref ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \date    11.06
!> \author  fschiff
!> \version 1.0
! *****************************************************************************
  SUBROUTINE mixed_init_subsys ( mixed_env, subsys, cell, cell_ref,  &
       force_env_section, subsys_section)

    TYPE(mixed_environment_type), POINTER    :: mixed_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section

    CHARACTER(len=*), PARAMETER :: routineN = 'mixed_init_subsys', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(mixed_energy_type), POINTER         :: mixed_energy
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CALL timeset(routineN,handle)
    NULLIFY ( mixed_energy, local_molecules, local_particles)
    particle_set => subsys % particles % els
    atomic_kind_set => subsys % atomic_kinds % els
    molecule_set => subsys % molecules_new % els
    molecule_kind_set => subsys % molecule_kinds_new % els

    ! Create the mixed_energy_type
    CALL allocate_mixed_energy( mixed_energy)

    ! Print the molecule kind set
    CALL write_molecule_kind_set(molecule_kind_set,subsys_section)

    ! Distribute molecules and atoms using the new data structures ***
    CALL distribute_molecules_1d (atomic_kind_set=atomic_kind_set,&
         particle_set=particle_set,&
         local_particles=local_particles,&
         molecule_kind_set=molecule_kind_set,&
         molecule_set=molecule_set,&
         local_molecules=local_molecules,&
         force_env_section=force_env_section)

    CALL cp_subsys_set(subsys, cell=cell)

    ! set the mixed_env
    CALL set_mixed_env ( mixed_env = mixed_env, subsys = subsys)
    CALL set_mixed_env ( mixed_env=mixed_env,&
         cell_ref=cell_ref, &
         local_molecules=local_molecules,&
         local_particles=local_particles,&
         mixed_energy=mixed_energy)

    CALL distribution_1d_release ( local_particles)
    CALL distribution_1d_release ( local_molecules)

    CALL timestop(handle)

  END SUBROUTINE mixed_init_subsys

END MODULE mixed_environment

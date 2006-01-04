
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****** cp2k/fist_environment [1.0] *
!!
!!   NAME
!!     fist_environment
!!
!!   FUNCTION
!!    initialize fist environment             
!!   AUTHOR
!!     CJM
!!   HISTORY 
!!   SOURCE
!******************************************************************************

MODULE fist_environment

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             write_atomic_kind_set
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_subsystem_types,              ONLY: cp_subsys_create,&
                                             cp_subsys_set,&
                                             cp_subsys_release,&
                                             cp_subsystem_type
  USE cp_subsystem_methods,            ONLY: cp_subsys_read_colvar
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d
  USE force_fields,                    ONLY: force_field_control
  USE global_types,                    ONLY: global_environment_type
  USE ewalds,                          ONLY: ewald_self
  USE ewald_environment_types,         ONLY: ewald_environment_type, &
                                             ewald_env_create, &
                                             ewald_env_release, &
                                             read_ewald_section,&
                                             ewald_env_get, ewald_env_set
  USE ewald_pw_types,                  ONLY: ewald_pw_type, &
                                             ewald_pw_create, &
                                             ewald_pw_release
  USE fist_environment_types,          ONLY: init_fist_env,&
                                             get_fist_env,&
                                             fist_environment_type,&
                                             set_fist_env  
  USE fist_energy_types,               ONLY: fist_energy_type,    &
                                             allocate_fist_energy
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type, &
                                             fist_nonbond_env_release
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE machine,                         ONLY: m_flush
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
                                             write_particle_distances, &
                                             write_qs_particle_coordinates,&
                                             write_structure_data
  USE simulation_cell,                 ONLY: cell_type,&
                                             cell_release,&
                                             read_cell,&
                                             write_cell,&
                                             get_cell
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: topology_control
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
#include "cp_common_uses.h"
  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: moduleN = "fist_environment"
  PUBLIC :: fist_init
  
!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
! FIST  FIST  FIST  FIST  FIST  FIST  FIST  FIST  FIST  FIST  FIST  FIST  FIST 
!-----------------------------------------------------------------------------!
!!
!!   NAME
!!     fist_init
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!   reads the input and database file for fist
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   USED BY
!!     fist_main
!!
!!*** *************************************************************************

SUBROUTINE fist_init ( fist_env, globenv, error )

  IMPLICIT NONE

! Arguments

  TYPE(fist_environment_type), POINTER           :: fist_env
  TYPE(global_environment_type), &
      POINTER                                    :: globenv
  TYPE(cp_error_type), INTENT(inout), OPTIONAL   :: error

! Local variables 
  LOGICAL                                  :: use_ref_cell, qmmm
  TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
  TYPE(atomic_kind_type), DIMENSION(:), &
    POINTER                                :: atomic_kind_set
  TYPE(cell_type), POINTER                 :: cell, cell_ref
  TYPE(cp_subsystem_type), POINTER         :: subsystem
  TYPE(mol_kind_new_list_type), POINTER    :: mol_kinds
  TYPE(mol_new_list_type), POINTER         :: mols
  TYPE(molecule_kind_type), DIMENSION(:), &
    POINTER                                :: molecule_kind_set
  TYPE(molecule_type), DIMENSION(:), &
    POINTER                                :: molecule_set
  TYPE(particle_list_type), POINTER        :: particles
  TYPE(particle_type), DIMENSION(:), &
    POINTER                                :: particle_set
  TYPE ( distribution_1d_type ), POINTER :: local_molecules,&
                                            local_particles
  TYPE ( ewald_environment_type ), POINTER :: ewald_env
  TYPE ( fist_nonbond_env_type ), POINTER ::  fist_nonbond_env
  INTEGER                                  :: handle
  TYPE(section_vals_type), POINTER         :: cell_section, &
                                              poisson_section, &
                                              colvar_section
! Local parameters 
  REAL(KIND=dp), DIMENSION(3)              :: abc
  CHARACTER(LEN=*), PARAMETER :: routine_name = "fist_init"
  TYPE (qmmm_env_mm_type), POINTER :: qmmm_env
!   ---------------------------------------------------------------------------
  CALL write_checkpoint_information("entering "//routine_name,globenv)

  CALL timeset(routine_name,"I",'',handle)


  NULLIFY ( subsystem, molecule_kind_set, molecule_set, atomic_kind_set,&
            particle_set,atomic_kinds,particles,mols,mol_kinds,cell,cell_ref)
  NULLIFY ( ewald_env, fist_nonbond_env, qmmm_env, cell_section, poisson_section)

  cell_section => section_vals_get_subs_vals(globenv%input_file,&
       "FORCE_EVAL%SUBSYS%CELL",error=error)
  poisson_section => section_vals_get_subs_vals(globenv%input_file,&
       "FORCE_EVAL%MM%POISSON_MM",error=error)
  CALL cp_subsys_create ( subsystem, para_env = globenv % para_env, error = error )

  ! read colvars and add to subsystem
  colvar_section => section_vals_get_subs_vals(globenv%input_file,&
       "FORCE_EVAL%SUBSYS%COLVAR",error=error)
  CALL cp_subsys_read_colvar( subsystem, colvar_section, error=error)

  CALL read_cell( cell, cell_ref, globenv, use_ref_cell = use_ref_cell, error=error,&
       cell_section=cell_section)
  CALL get_cell (cell, abc=abc)

! Print the cell parameters ***
  CALL write_cell ( cell, globenv )

! Create the ewald environment
  CALL ewald_env_create ( ewald_env, globenv % para_env, error )

! Read the input section and set the ewald_env ***
  CALL read_ewald_section ( ewald_env, globenv, error )
  CALL ewald_env_set(ewald_env, poisson_section=poisson_section, error=error)

! topology
  CALL get_fist_env (fist_env, qmmm=qmmm, qmmm_env=qmmm_env, error=error)

  CALL topology_control ( atomic_kind_set, particle_set,&
                          molecule_kind_set, molecule_set,&
                          globenv, qmmm=qmmm, qmmm_env=qmmm_env,&
                          error=error)

  CALL force_field_control ( atomic_kind_set, particle_set,&
                             molecule_kind_set, molecule_set,&
                             ewald_env, fist_nonbond_env, globenv, qmmm=qmmm,&
                             qmmm_env=qmmm_env, error=error )

  CALL particle_list_create ( particles, els_ptr = particle_set, error = error )
  CALL atomic_kind_list_create ( atomic_kinds, els_ptr = atomic_kind_set, error = error )
  CALL mol_new_list_create ( mols, els_ptr = molecule_set, error = error )
  CALL mol_kind_new_list_create ( mol_kinds, els_ptr = molecule_kind_set, error = error )
  CALL cp_subsys_set ( subsystem, particles = particles, atomic_kinds = atomic_kinds, &
       molecules_new = mols, molecule_kinds_new = mol_kinds, error = error )
  CALL particle_list_release ( particles, error = error )
  CALL atomic_kind_list_release ( atomic_kinds, error = error )
  CALL mol_new_list_release ( mols, error = error )
  CALL mol_kind_new_list_release ( mol_kinds, error = error )

  CALL fist_init_subsys ( fist_env, subsystem, cell, cell_ref, use_ref_cell, &
                          fist_nonbond_env, ewald_env, globenv, error )

  CALL cell_release(cell,error=error)
  CALL cell_release(cell_ref,error=error)
  CALL ewald_env_release ( ewald_env, error )
  CALL fist_nonbond_env_release ( fist_nonbond_env, error )
  CALL cp_subsys_release(subsystem,error=error)

  CALL timestop(0.0_dp,handle)
 
  END SUBROUTINE fist_init

!   ***************************************************************************
  SUBROUTINE fist_init_subsys ( fist_env, subsystem, cell, cell_ref, use_ref_cell, &
                                fist_nonbond_env, ewald_env, globenv, error )

!   Purpose: Read the input and the database files for the setup of the
!            FIST environment.

!   History: - Creation (22.05.2000,MK)

!   ***************************************************************************
  TYPE ( fist_environment_type ), POINTER       :: fist_env
  TYPE ( cp_subsystem_type ), POINTER           :: subsystem
  TYPE ( cell_type ), POINTER                   :: cell, cell_ref
  TYPE ( fist_nonbond_env_type ), POINTER       :: fist_nonbond_env
  TYPE ( ewald_environment_type ), POINTER      :: ewald_env
  LOGICAL, INTENT ( IN )                        :: use_ref_cell
  TYPE ( global_environment_type ), &
    pointer                            :: globenv
  TYPE(cp_error_type), INTENT (INOUT), &
    OPTIONAL                                    :: error

! locals
  CHARACTER(LEN=*), PARAMETER              :: routine_name = "fist_init_subsys"

  INTEGER :: handle, istat, natom_types
  INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom_of_kind
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
  TYPE(fist_energy_type), POINTER          :: thermo
  TYPE ( ewald_pw_type ), POINTER          :: ewald_pw


  CALL write_checkpoint_information("entering "//routine_name,globenv)
  CALL timeset(routine_name,"I",'',handle)
  particle_set => subsystem % particles % els
  atomic_kind_set => subsystem % atomic_kinds % els
  molecule_set => subsystem % molecules_new % els
  molecule_kind_set => subsystem % molecule_kinds_new % els

  NULLIFY ( thermo, ewald_pw, local_molecules, local_particles)

! Create the fist_energy_type
  call allocate_fist_energy( thermo )

! Print the atomic kind set
!  CALL write_atomic_kind_set(atomic_kind_set,globenv)

! Print the molecule kind set
  CALL write_molecule_kind_set(molecule_kind_set,globenv)

! Print the atomic coordinates 
  CALL write_qs_particle_coordinates(particle_set,cell,globenv,error)
  CALL write_particle_distances(particle_set,cell,globenv,error)
  CALL write_structure_data(particle_set,cell,globenv)

! Distribute molecules and atoms using the new data structures ***

  CALL distribute_molecules_1d ( particle_kind_set=atomic_kind_set,&
                                 particle_set=particle_set,&
                                 local_particles=local_particles,&
                                 molecule_kind_set=molecule_kind_set,&
                                 molecule_set=molecule_set,&
                                 local_molecules=local_molecules,&
                                 error=error )

! Initialize ewald grids
  CALL ewald_pw_create ( ewald_pw, cell, cell_ref, ewald_env, error )

! Initialize self-interaction and neutralizing background
  CALL ewald_self ( ewald_env, atomic_kind_set, local_particles, &
                    thermo % e_self, thermo % e_neut )

! set the fist_env
  CALL set_fist_env ( fist_env = fist_env, subsys = subsystem, error = error )
  CALL set_fist_env ( fist_env=fist_env,&
                      cell=cell, cell_ref=cell_ref, &
                      local_molecules=local_molecules,&
                      local_particles=local_particles,&
                      ewald_env=ewald_env, ewald_pw=ewald_pw, &
                      fist_nonbond_env=fist_nonbond_env, thermo=thermo )  
    
  CALL distribution_1d_release ( local_particles, error = error )
  CALL distribution_1d_release ( local_molecules, error = error )
  CALL ewald_pw_release ( ewald_pw, error )

  IF (globenv%ionode) CALL m_flush(globenv%scr)

  CALL timestop(0.0_dp,handle)

  CALL write_checkpoint_information("leaving "//routine_name,globenv)

 END SUBROUTINE fist_init_subsys
!******************************************************************************
END MODULE fist_environment
!******************************************************************************

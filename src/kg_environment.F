!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_environment [1.0] *
!!
!!   NAME
!!     kg_environment
!!
!!   FUNCTION
!!    initialize kg environment             
!!   AUTHOR
!!     gt SEPT-23-2002
!!   HISTORY 
!!     GT 11-13-2003 :            
!!     moved the main driver to cp2k 
!!     performs only the initialization of the environment
!!   SOURCE
!******************************************************************************

MODULE kg_environment

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set,&
                                             init_atomic_kind_set,&
                                             read_atomic_kind_set,&
                                             write_atomic_kind_set,&
                                             write_gto_basis_sets
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_control_types,                ONLY: dft_control_release,&
                                             dft_control_type
  USE cp_control_utils,                ONLY: read_dft_control,&
                                             read_kg_control,&
                                             write_dft_control,&
                                             write_kg_control
  USE cp_error_handling,               ONLY: cp_error_type
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d
! USE force_env_types,                   ONLY: force_env_calculate_force,&

!                                            force_env_create,&

!                                            force_env_release,&

!                                            force_env_type

  USE global_types,                    ONLY: global_environment_type
!  USE kg_energy_utils,                 ONLY: energy_min_type,&

!                                             min_info_type,&

!                                             print_convergence,&

!                                             print_energy_info,&

!                                             read_energy_section

  USE kg_energy_types,                 ONLY: allocate_kg_energy,&
                                             kg_energy_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type,&
                                             set_kg_env
  USE kg_force_types,                  ONLY: allocate_kg_force,&
                                             kg_force_type
  USE kg_rho_methods,                  ONLY: kg_rho_create,&
                                             kg_rho_setup
  USE kg_rho_types,                    ONLY: kg_rho_release,&
                                             kg_rho_type
!  USE kg_force,                        ONLY: kg_force_release

!  USE kg_pol_coefs,                    ONLY: allocate_coefs,&

!                                             dump_coefs,&

!                                             initialize_coefs

  USE kinds,                           ONLY: dp,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE molecule_kind_types,             ONLY: molecule_kind_type,&
                                             write_molecule_kind_set
  USE molecule_types_new,              ONLY: molecule_type
  USE orbital_pointers,                ONLY: init_orbital_pointers
  USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_distances
  USE qs_interactions,                 ONLY: init_interaction_radii
  USE simulation_cell,                 ONLY: cell_release,&
                                             cell_type,&
                                             read_cell,&
                                             write_cell
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: topology_control

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: moduleN = "kg_environment"
  PUBLIC :: kg_init
  
!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
! KG KG KG                                                                    !
!-----------------------------------------------------------------------------!
!!
!!   NAME
!!     kg_init
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!   reads the input and database file for KG
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   USED BY
!!     kg_main
!!
!!*** *************************************************************************

SUBROUTINE kg_init (kg_env,globenv,error)



    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "kg_init"

    INTEGER                                  :: handle, i, istat, maxl, &
                                                maxlgto, maxlppl, maxlppnl, &
                                                nkind, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom_of_kind
    LOGICAL                                  :: ionode, kgpol
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(kg_energy_type), POINTER            :: energy
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!  CHARACTER ( LEN = 40 ) :: set_fn, project_name
!  CHARACTER ( LEN = 5 ) :: tag
!------------------------------------------------------------------------------

  CALL write_checkpoint_information("entering "//routine_name,globenv)

  CALL timeset(routine_name,"I",'',handle)

  ionode = globenv%ionode
  output_unit = globenv%scr


! initialize kg_env
! CALL init_kg_env ( kg_env )
! nullifying  pointers
  NULLIFY ( atomic_kind_set )
  NULLIFY ( cell )
  NULLIFY ( dft_control )
  NULLIFY (energy)
  NULLIFY (force)
  NULLIFY (local_molecules)
  NULLIFY (local_particles)
  NULLIFY (molecule_set)
  NULLIFY ( particle_set )

  kgpol = .FALSE.

! Read the input section with the dft control parameters *
  CALL read_dft_control ( dft_control, globenv )

! Print the DFT  parameters 
  CALL write_dft_control ( dft_control, globenv )

! Read the input section with the KG control parameters 
  CALL read_kg_control ( dft_control % kg_control, globenv )

! Print the kg parameters 
  CALL write_kg_control ( dft_control % kg_control, globenv )

  kgpol = dft_control % kg_control % polarization

! Read the input section with the cell parameters ***
  CALL read_cell ( cell, globenv )

! Print the cell parameters ***
  CALL write_cell ( cell, globenv )

! CALL read_setup_section ( setup, globenv )
 
! topology
  CALL topology_control (atomic_kind_set,particle_set,&
                         molecule_kind_set,molecule_set,&
                         globenv)
! CALL xstring ( topology % coord_file_name, ia, ib )
! project_name = topology % coord_file_name ( ia : ib )
! topology%cell => cell


! Read the atomic kind set
  CALL read_atomic_kind_set(atomic_kind_set,globenv)

! Write the atomic kind set
  CALL write_gto_basis_sets(atomic_kind_set,globenv)

! Initialize the spherical harmonics
! the orbital transformation matrices
  CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxlgto=maxlgto,&
                             maxlppl=maxlppl,&
                             maxlppnl=maxlppnl)

  maxl = MAX(maxlgto,maxlppl,maxlppnl)

  IF (dft_control%forces) maxl = maxl + 1

  CALL init_orbital_pointers(maxl)
  CALL init_spherical_harmonics(maxl,globenv)

! Initialise the atomic kind set
  CALL init_atomic_kind_set(atomic_kind_set,globenv)

! Print the atomic kind set
  CALL write_atomic_kind_set(atomic_kind_set,globenv)

! Print the molecule kind set
  CALL write_molecule_kind_set(molecule_kind_set,globenv)

! Print the total number of kinds, atoms, basis functions etc.
!  CALL write_total_numbers(atomic_kind_set,particle_set,globenv)

! Print the atomic coordinates 
!   IF (globenv%print%atomic_coordinates) THEN
!     CALL write_qs_particle_coordinates(particle_set,cell,globenv)
!   END IF

! Print the interatomic distances ***

  CALL write_particle_distances(particle_set,cell,globenv)


! Initialize the atomic interaction radii 
  CALL init_interaction_radii ( dft_control%kg_control, cell, &
                                atomic_kind_set, globenv)

! Distribute molecules and atoms using the new data structures ***

  CALL distribute_molecules_1d(particle_kind_set=atomic_kind_set,&
                                 particle_set=particle_set,&
                                 local_particles=local_particles,&
                                 molecule_kind_set=molecule_kind_set,&
                                 molecule_set=molecule_set,&
                                 local_molecules=local_molecules,&
                                 error=error)

! If polarization is also present, allocate and initialize coefs pointer  

  IF ( kgpol ) THEN


!    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
!                             ncgf_aux=ncoef)
!    CALL allocate_coefs (ncoef, coef_pos ( 1 ), coef_vel ( 1 ), &
!                             coef_force ( 1 ) )
!    CALL initialize_coefs (simpar, coef_pos ( 1 ), coef_vel ( 1 ),&
!                             coef_force ( 1 ), globenv )
!    IF ( globenv % run_type == 'WAVEFUNCTION OPTIMIZATION') THEN
!      CALL read_energy_section ( energy_min, globenv )
!    END IF

  END IF

! Allocate the data structure for energies ***
  CALL allocate_kg_energy(energy)

!   *** Allocate the force data structure ***

  IF (dft_control%forces) THEN
    nkind = SIZE(atomic_kind_set)
    ALLOCATE (natom_of_kind(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,moduleN,__LINE__,&
                                     "natom_of_kind",nkind*int_size)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom_of_kind=natom_of_kind)
    CALL allocate_kg_force(force,natom_of_kind)
    DEALLOCATE (natom_of_kind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,moduleN,__LINE__,&
                                     "natom_of_kind")
  END IF

! set the kg_env
  CALL set_kg_env ( kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell, &
                    dft_control=dft_control, &
                    energy=energy,&
                    force=force,&
                    local_molecules=local_molecules,&
                    local_particles=local_particles,&
                    molecule_kind_set=molecule_kind_set,&
                    molecule_set=molecule_set,&
                    particle_set=particle_set)  
    
  CALL cell_release(cell,error=error)
  CALL distribution_1d_release(local_particles,error=error)
  CALL distribution_1d_release(local_molecules,error=error)
  CALL dft_control_release(dft_control, error=error)

!Sets up pw, grids, densities              

  CALL kg_env_setup(kg_env,error)

  IF (ionode) CALL m_flush(output_unit)

  CALL timestop(0.0_dp,handle)

  CALL write_checkpoint_information("leaving "//routine_name,globenv)

 END SUBROUTINE kg_init
!***************************************************************************
SUBROUTINE kg_env_setup(kg_env,error)

    TYPE(kg_environment_type), INTENT(inout) :: kg_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_env_setup', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(kg_rho_type), POINTER               :: rho

  failure=.FALSE.
  NULLIFY(rho)

  CALL get_kg_env(kg_env=kg_env,&
                  rho=rho)

  CALL kg_rho_create(kg_rho=rho,kg_env=kg_env,error=error)
  CALL kg_rho_setup(kg_rho=rho,kg_env=kg_env,error=error)
  CALL set_kg_env(kg_env,rho=rho)
  CALL kg_rho_release(kg_rho=rho,error=error)
  END SUBROUTINE kg_env_setup

!******************************************************************************
END MODULE kg_environment
!******************************************************************************


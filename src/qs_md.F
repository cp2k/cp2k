!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_md [1.0] *
!!
!!   NAME
!!     qs_md
!!
!!   FUNCTION
!!     Perform a molecular dynamics (MD) run using QUICKSTEP
!!
!!   AUTHOR
!!     Matthias Krack (07.11.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_md
  USE atoms_input,                     ONLY: set_particle_properties
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cntl_input,                      ONLY: read_cntl_section
  USE convert_units,                   ONLY: convert
  USE extended_system_types,           ONLY: extended_system_type
  USE force_control,                   ONLY: force,&
                                             set_force_env
  USE global_types,                    ONLY: global_environment_type
  USE initialize_extended_types,       ONLY: initialize_nhc_baro,&
                                             initialize_nhc_part,&
                                             initialize_npt_type
  USE initialize_molecule_types,       ONLY: initialize_molecule_type
  USE input_types,                     ONLY: setup_parameters_type
  USE integrator,                      ONLY: set_integrator
  USE kinds,                           ONLY: int_size,&
                                             wp => dp,&
                                             wp_size => dp_size
  USE md,                              ONLY: initialize_velocities,&
                                             mdio_parameters_type,&
                                             read_md_section,&
                                             simulation_parameters_type,&
                                             virial_type
  USE md_environment_types,            ONLY: init_md_env,&
                                             md_environment_type,&
                                             set_md_env,&
                                             zero_virial
  USE molecule_input,                  ONLY: read_molecule_section,&
                                             read_setup_section,&
                                             read_topology_section
  USE molecule_types,                  ONLY: intra_parameters_type,&
                                             molecule_structure_type,&
                                             molecule_type,&
                                             particle_node_type,&
                                             topology_type
  USE particle_types,                  ONLY: particle_prop_type,&
                                             particle_type
  USE qs_energy,                       ONLY: qs_energies
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force,                        ONLY: qs_forces
  USE simulation_cell,                 ONLY: cell_type
  USE structure_types,                 ONLY: init_structure_type,&
                                             set_structure_type,&
                                             structure_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE unit,                            ONLY: set_units,&
                                             unit_convert_type
  USE velocity_verlet_control,         ONLY: velocity_verlet
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_md"

  PUBLIC :: qs_mol_dyn

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE qs_mol_dyn(qs_env,globenv)

!   Purpose: Driver routine for MD run using QUICKSTEP.

!   History: - Creation (07.11.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(INOUT) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)     :: qs_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_mol_dyn"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(extended_system_type)  :: extended_type
    TYPE(intra_parameters_type) :: intra_param
    TYPE(md_environment_type)   :: md_env
    TYPE(mdio_parameters_type)  :: mdio
    TYPE(setup_parameters_type) :: setup
    TYPE(topology_type)         :: topo
    TYPE(unit_convert_type)     :: units
    TYPE(virial_type)           :: virial

    TYPE(cell_type), POINTER                  :: cell
    TYPE(simulation_parameters_type), POINTER :: simpar

    TYPE(molecule_structure_type), DIMENSION(:), POINTER :: molecule
    TYPE(molecule_type), DIMENSION(:), POINTER           :: mol_setup
    TYPE(particle_node_type), DIMENSION(:), POINTER      :: pnode
    TYPE(particle_prop_type), DIMENSION(:), POINTER      :: pstat
    TYPE(particle_type), DIMENSION(:), POINTER           :: particle_set
    TYPE(structure_type), DIMENSION(:), POINTER          :: struc
    REAL(wp), POINTER                                    :: constant
    INTEGER, POINTER                                     :: itimes

    INTEGER :: handle,istat,natom,nmolecule

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    NULLIFY (molecule)
    NULLIFY (mol_setup)
    NULLIFY (pnode)
    NULLIFY (pstat)

    ALLOCATE (simpar,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"simpar",0)
    ALLOCATE (struc(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"struc",0)

    CALL read_cntl_section(setup,globenv)

    CALL read_setup_section(mol_setup,setup,globenv)

    CALL read_molecule_section(mol_setup,setup,globenv)

    CALL read_topology_section(topo,setup,globenv)

    CALL read_md_section(simpar,globenv,mdio)
    simpar%program = globenv%program_name

    CALL set_particle_properties(setup,mol_setup,pstat)

    CALL set_units(setup%unit_type,units)
    CALL convert(units=units,simpar=simpar)

    CALL get_qs_env(qs_env=qs_env,&
                    cell=cell,&
                    particle_set=particle_set)

    natom = SIZE(particle_set)

    ALLOCATE (pnode(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pnode",int_size*natom)

    nmolecule = SIZE(mol_setup)

    ALLOCATE (molecule(nmolecule),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"molecule",int_size*nmolecule)

    CALL initialize_molecule_type(mol_setup,intra_param,pnode,&
                                  particle_set,molecule,globenv)

    ALLOCATE (extended_type%nhc_part(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nhc_part",0)

    CALL initialize_nhc_part(cell,simpar,molecule,mol_setup,globenv,&
                             extended_type%nhc_part(1))

    ALLOCATE (extended_type%nhc_baro(1),STAT=istat)
    IF (istat /= 0 ) CALL stop_memory(routine,"nhc_baro",0)

    CALL initialize_nhc_baro(simpar,globenv,extended_type%nhc_baro(1))

    CALL initialize_npt_type(simpar,globenv,extended_type%npt_info,cell)

    IF ((simpar%read_type == "POS").OR.&
        (simpar%read_type == "INIT")) THEN
      CALL initialize_velocities(simpar,particle_set,globenv)
    END IF

    CALL set_force_env(qs_env_type=qs_env)

    CALL init_md_env(md_env)

    CALL set_md_env(md_env=md_env,&
                    simpar=simpar,&
                    struc=struc,&
                    nhc_part=extended_type%nhc_part,&
                    nhc_baro=extended_type%nhc_baro,&
                    npt=extended_type%npt_info,&
                    cell=cell)

    ALLOCATE (md_env%itimes,STAT=istat)
    IF (istat /= 0) CALL stop_memory (routine,"itimes",int_size)
    ALLOCATE (md_env%constant,STAT=istat)
    IF (istat /= 0) CALL stop_memory (routine,"constant",wp_size)

    CALL set_integrator(globenv)

    itimes => md_env%itimes
    itimes = 0

    CALL qs_forces(qs_env,globenv)

    DO itimes=1,md_env%simpar%nsteps
PRINT*,"Step: :",itimes," of ",md_env%simpar%nsteps
      CALL velocity_verlet(md_env)

    END DO

    CALL timestop(0.0_wp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE qs_mol_dyn

! *****************************************************************************

END MODULE qs_md

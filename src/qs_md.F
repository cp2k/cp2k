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
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cntl_input,                      ONLY: read_cntl_section
  USE convert_units,                   ONLY: convert
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
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE extended_system_types,           ONLY: extended_system_type
  USE force_control,                   ONLY: force_env_create,&
                                             force_env_release,&
                                             force_env_type
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
                                             qs_md_parameters_type,&
                                             read_md_section,&
                                             simulation_parameters_type,&
                                             virial_type
  USE md_conserved_quantities,         ONLY: compute_conserved_quantity
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
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_coordinates
  USE qs_energy,                       ONLY: qs_energies
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force,                        ONLY: qs_forces
  USE qs_wf_history_methods,           ONLY: wfi_create
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                             wfi_release
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

  SUBROUTINE qs_mol_dyn(qs_env,globenv, error)

!   Purpose: Driver routine for MD run using QUICKSTEP.

!   History: - Creation (07.11.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(INOUT) :: globenv
    TYPE(qs_environment_type), POINTER           :: qs_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

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

    TYPE(cell_type), POINTER                  :: cell
    TYPE(simulation_parameters_type), POINTER :: simpar

    TYPE(molecule_structure_type), DIMENSION(:), POINTER :: molecule
    TYPE(molecule_type), DIMENSION(:), POINTER           :: mol_setup
    TYPE(particle_node_type), DIMENSION(:), POINTER      :: pnode
    TYPE(particle_type), DIMENSION(:), POINTER           :: particle_set
    TYPE(structure_type), DIMENSION(:), POINTER          :: struc
    TYPE(virial_type), POINTER                           :: virial
    REAL(wp), POINTER                                    :: constant
    INTEGER, POINTER                                     :: itimes
    REAL(wp)                                             :: ekin
    REAL(wp)                                             :: econs 
    INTEGER                                              :: i,nnodes
    TYPE(qs_energy_type), POINTER                        :: energy
    TYPE(qs_wf_history_type), POINTER                    :: wf_history
    TYPE(force_env_type), POINTER                        :: force_env

    INTEGER :: handle,istat,natom,nmolecule

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    NULLIFY (molecule)
    NULLIFY (mol_setup)
    NULLIFY (pnode)
    NULLIFY (energy, wf_history, force_env)

    ALLOCATE (simpar,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"simpar",0)
    ALLOCATE (struc(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"struc",0)
    ALLOCATE (virial,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"virial",0)
    CALL zero_virial(virial)

    CALL read_cntl_section(setup,globenv)

    CALL read_setup_section(mol_setup,setup,globenv)

    CALL read_molecule_section(mol_setup,setup,globenv)

    CALL read_topology_section(topo,setup,globenv)

    CALL read_md_section(simpar,globenv,mdio)
    simpar%program = globenv%program_name

    CALL set_units(setup%unit_type,units)
    CALL convert(units=units,simpar=simpar)

    ! activate the interpolation
    CALL wfi_create(wf_history, &
         interpolation_method_nr=simpar%qs_md_param%wf_interpolation_method_nr,&
         error=error)
    CALL set_qs_env(qs_env, wf_history=wf_history, error=error)

    CALL get_qs_env(qs_env=qs_env,&
                    cell=cell,&
                    particle_set=particle_set,energy=energy)

    ! is not supposed to be working in parallel case !
    IF (globenv%num_pe .NE. 1) CALL stop_program("qs_md","MD not parallel")
    natom = SIZE(particle_set)

    ALLOCATE (pnode(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pnode",int_size*natom)

    nmolecule = SUM(mol_setup(:)%num_mol)

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

    CALL force_env_create(force_env, qs_env=qs_env, error=error)

    ALLOCATE (struc(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"struc",0)
    CALL init_structure_type(struc(1))
    struc(1)%pnode => pnode
    struc(1)%part =>  particle_set

    CALL init_md_env(md_env)

    CALL set_md_env(md_env=md_env,&
                    simpar=simpar,&
                    struc=struc,&
                    virial=virial,&
                    nhc_part=extended_type%nhc_part,&
                    nhc_baro=extended_type%nhc_baro,&
                    npt=extended_type%npt_info,&
                    cell=cell,&
                    force_env=force_env)

    ALLOCATE (md_env%itimes,STAT=istat)
    IF (istat /= 0) CALL stop_memory (routine,"itimes",int_size)
    ALLOCATE (md_env%constant,STAT=istat)
    IF (istat /= 0) CALL stop_memory (routine,"constant",wp_size)
    CALL force_env_release(force_env,error=error)

    CALL set_integrator(globenv)

    itimes => md_env%itimes
    itimes = 0

!FM    CALL qs_forces(qs_env,globenv)

    DO itimes=1,md_env%simpar%nsteps
      PRINT*,"Step: :",itimes," of ",md_env%simpar%nsteps

      CALL velocity_verlet(md_env)

      CALL compute_conserved_quantity (md_env, energy%total, econs, ekin=ekin)
      WRITE(6,'(A,3F20.12)') "XXXXXXXXXXXXXXXX EKIN,EPOT,TOTAL ",ekin,energy%total,econs
!!      nnodes = SIZE ( pnode)
!!      ekin=0.0_wp
!!      DO i = 1, nnodes
!!           ekin = ekin + 0.5_wp*pnode(i)%p%prop%mass*( &
!!           pnode(i)%p%v(1)*pnode(i)%p%v(1)+pnode(i)%p%v(2)* &
!!           pnode(i)%p%v(2)+pnode(i)%p%v(3)*pnode(i)%p%v(3))
!!      END DO
!!    write(6,'(A,3F20.12)') "XXXXXXXXXXXXXXXX EKIN,EPOT,TOTAL ",ekin,energy%total,energy%total+ekin
      CALL write_particle_coordinates(qs_env%particle_set,qs_env%cell,globenv)
    END DO

    CALL timestop(0.0_wp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE qs_mol_dyn

! *****************************************************************************

END MODULE qs_md

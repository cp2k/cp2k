!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
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
!  USE initialize_molecule_types,       ONLY: initialize_molecule_type
! cjm   USE util,                            ONLY: get_share
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE atoms_input,                     ONLY: read_md_restart
  USE checkpoint_handler,              ONLY: write_checkpoint_information
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
  USE dump,                            ONLY: dump_variables
  USE extended_system_dynamics,        ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE extended_system_types,           ONLY: extended_system_type
  USE force_control,                   ONLY: force_env_calculate_force,&
                                             force_env_create,&
                                             force_env_release,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE initialize_extended_types,       ONLY: initialize_nhc_baro,&
                                             initialize_nhc_forces,&
                                             initialize_nhc_part,&
                                             initialize_npt_type
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
  USE md_environment_types,            ONLY: destroy_md_env,&
                                             init_md_env,&
                                             md_environment_type,&
                                             set_md_env,&
                                             zero_virial
  USE md_qs_energies,                  ONLY: md_energy_qs,&
                                             set_qs_energies
  USE message_passing,                 ONLY: mp_bcast
  USE molecule_types,                  ONLY: intra_parameters_type,&
                                             molecule_structure_type,&
                                             molecule_type,&
                                             particle_node_type
  USE output_utilities,                ONLY: print_warning
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: kelvin
  USE qs_energy,                       ONLY: qs_energies
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force,                        ONLY: qs_forces
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE simulation_cell,                 ONLY: cell_type
  USE structure_types,                 ONLY: init_structure_type,&
                                             set_structure_type,&
                                             structure_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
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

  SUBROUTINE qs_mol_dyn ( qs_env, mol_setup, globenv, error )

!   Purpose: Driver routine for MD run using QUICKSTEP.

!   History: - Creation (07.11.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(INOUT) :: globenv
    TYPE ( molecule_type ), DIMENSION ( : ), POINTER :: mol_setup
    TYPE(qs_environment_type), POINTER           :: qs_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_mol_dyn"

!   *** Local variables ***

    TYPE ( lnhc_parameters_type ), POINTER :: nhc_part, nhc_baro
    TYPE ( npt_info_type ), DIMENSION ( :, : ), POINTER :: npt_info
    TYPE(extended_system_type)  :: extended_type
    TYPE(intra_parameters_type) :: intra_param
    TYPE(md_environment_type)   :: md_env
    TYPE(mdio_parameters_type)  :: mdio
    TYPE(setup_parameters_type) :: setup

    TYPE(cell_type), POINTER                  :: cell
    TYPE(cell_type), POINTER                  :: cell_ref
    TYPE(simulation_parameters_type), POINTER :: simpar

    TYPE(molecule_structure_type), DIMENSION(:), POINTER :: molecule
    TYPE(particle_node_type), DIMENSION(:), POINTER      :: pnode
    TYPE(particle_type), DIMENSION(:), POINTER           :: particle_set
    TYPE(structure_type), DIMENSION(:), POINTER          :: struc
    TYPE(virial_type), POINTER                           :: virial
    REAL(wp), POINTER                                    :: constant
    INTEGER, POINTER                                     :: itimes
    REAL(wp)                                             :: ekin
    REAL(wp)                                             :: econs 
    INTEGER                                              :: i
    TYPE(qs_energy_type), POINTER                        :: energy
    TYPE(force_env_type), POINTER                        :: force_env

    INTEGER :: handle,istat,istep,natom,nmolecule
    INTEGER :: nmoltype, nmol, nnodes, nsh
    LOGICAL :: should_stop

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    NULLIFY ( molecule )
    NULLIFY ( pnode )
    NULLIFY ( energy, force_env )
    NULLIFY ( nhc_part, nhc_baro, npt_info )

    ALLOCATE (simpar,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "simpar",0)
    ALLOCATE (struc(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "struc",0)
    ALLOCATE (virial,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "virial",0)
    CALL zero_virial(virial)

    CALL read_md_section(simpar,globenv,mdio)
    simpar%program = globenv%program_name

    CALL get_qs_env(qs_env=qs_env, cell=cell,  &
                    pnode = pnode, molecule = molecule, &
                    particle_set=particle_set,energy=energy)

! Define a reference cell for compatibility with MD restarts
    cell_ref => cell

!================ THIS NEEDS TO BE MOVED TO QS_ENV===================
    ! just to get a working parallel version.
    ! clearly a full clean-up of fist/qs/kg is needed 
!    natom = SIZE ( particle_set )
!    nmoltype = SIZE ( mol_setup )
!    nmol = 0
!    nnodes = 0
!    DO i = 1, nmoltype
!       nsh = get_share ( mol_setup ( i ) % num_mol, &
!            globenv % num_pe, globenv % mepos )
!       nmol = nmol + nsh
!       nnodes = nnodes + nsh * mol_setup ( i ) % molpar % natom
!    END DO
!    ALLOCATE (pnode(nnodes),STAT=istat)
!    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
!                                     "pnode",int_size*natom)
!    ALLOCATE (molecule(nmol),STAT=istat)
!    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
!                                     "molecule",int_size*nmolecule)
!
!    CALL initialize_molecule_type( mol_setup,intra_param,pnode,&
!                                   particle_set,molecule,globenv )
!    CALL set_qs_env ( qs_env, molecule = molecule, error = error )
!================ THIS NEEDS TO BE MOVED TO QS_ENV===================

    ALLOCATE (extended_type%nhc_part(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nhc_part",0)

    CALL initialize_nhc_part(cell,simpar,molecule,mol_setup,globenv,&
                             extended_type%nhc_part(1))

    ALLOCATE (extended_type%nhc_baro(1),STAT=istat)
    IF (istat /= 0 ) CALL stop_memory(routine_name,module_name,__LINE__,&
                                      "nhc_baro",0)

    CALL initialize_nhc_baro(simpar,globenv,extended_type%nhc_baro(1))

    CALL initialize_npt_type(simpar,globenv,extended_type%npt_info,cell)

! initialize velocities and read restart files if needed
    npt_info => extended_type % npt_info 
    nhc_part => extended_type % nhc_part ( 1 )
    nhc_baro => extended_type % nhc_baro ( 1 ) 
    SELECT CASE ( simpar % read_type )
    CASE ( 'INIT' )
      CALL initialize_velocities ( simpar, particle_set, globenv )
    CASE ( 'POS' )
      CALL read_md_restart ( simpar % read_type, globenv,   &
                             particle_set, cell, cell_ref,  &
                             nhc_part, nhc_baro, npt_info )
      CALL initialize_velocities ( simpar, particle_set, globenv )
    CASE DEFAULT
      CALL read_md_restart ( simpar % read_type, globenv,   &
                             particle_set, cell, cell_ref,  &
                             nhc_part, nhc_baro, npt_info )
    END SELECT
    CALL initialize_nhc_forces ( nhc_part, nhc_baro )


    CALL force_env_create(force_env, qs_env=qs_env, globenv=globenv, error=error)

    ALLOCATE (struc(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "struc",0)
    CALL init_structure_type(struc(1))
    struc(1)%pnode => pnode
    struc(1)%part =>  particle_set
    struc(1)%molecule =>  molecule

    CALL init_md_env(md_env, globenv % para_env)

    CALL set_md_env(md_env=md_env,&
                    simpar=simpar,&
                    struc=struc,&
                    virial=virial,&
                    nhc_part=extended_type%nhc_part,&
                    nhc_baro=extended_type%nhc_baro,&
                    npt=extended_type%npt_info,&
                    cell=cell,&
                    force_env=force_env)

    ! this is in fact pretty ugly !!!
    ALLOCATE (md_env%itimes,STAT=istat)
    IF (istat /= 0) CALL stop_memory (routine_name,module_name,__LINE__,&
                                      "itimes",int_size)
    ALLOCATE (md_env%constant,STAT=istat)
    IF (istat /= 0) CALL stop_memory (routine_name,module_name,__LINE__,&
                                      "constant",wp_size)
    CALL force_env_release(force_env,error=error)

    CALL set_integrator(globenv)

    itimes => md_env%itimes
    itimes = 0

    CALL set_qs_energies ( globenv, mdio )

    ! the verlet integrator needs to be started with the correct forces alreay initialized
    CALL force_env_calculate_force ( md_env%force_env )
    CALL md_energy_qs(md_env) ! inits the print env at itimes == 0

    DO istep=1,md_env%simpar%nsteps
      CALL velocity_verlet(md_env)
      itimes = itimes + 1
      CALL temperature_control(simpar,qs_env,globenv)
      CALL md_energy_qs(md_env)
      IF (MODULO(itimes,mdio%idump) == 0) THEN
        CALL dump_variables(md_env,qs_env,mdio%dump_file_name,globenv)
      END IF
      CALL external_control(qs_env,globenv,should_stop)
      IF (should_stop) EXIT
    END DO

    CALL destroy_md_env(md_env)

    CALL timestop(0.0_wp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE qs_mol_dyn

! *****************************************************************************

  SUBROUTINE temperature_control(simpar,qs_env,globenv)

!   Purpose: Perform all temperature manipulations during a QS MD run.

!   History: - Creation (15.09.2003,MK)

!   ***************************************************************************

    TYPE(simulation_parameters_type), POINTER :: simpar
    TYPE(qs_environment_type), POINTER         :: qs_env
    TYPE(global_environment_type), INTENT(IN)  :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "temperature_control"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind

    REAL(wp) :: e_kinetic,mass,temperature,scale,v2
    INTEGER  :: atom,handle,iatom,ikind,istat,natom,nkind,output_unit
    LOGICAL  :: ionode

    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    INTEGER, DIMENSION(:), POINTER                :: atom_list

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

    IF (simpar%temp_tol > 0.0_wp) THEN

      CALL get_qs_env(qs_env=qs_env,&
                      atomic_kind_set=atomic_kind_set,&
                      particle_set=particle_set)

      nkind = SIZE(atomic_kind_set)
      natom = SIZE(particle_set)

!     *** Calculate the current temperature ***

      e_kinetic = 0.0_wp

      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             atom_list=atom_list,&
                             mass=mass)
        v2 = 0.0_wp
        natom = SIZE(atom_list)
        DO iatom=1,natom
          atom = atom_list(iatom)
          v2 = v2 + particle_set(atom)%v(1)*particle_set(atom)%v(1) +&
                    particle_set(atom)%v(2)*particle_set(atom)%v(2) +&
                    particle_set(atom)%v(3)*particle_set(atom)%v(3)
        END DO
        e_kinetic = e_kinetic + 0.5_wp*mass*v2
      END DO

      temperature = 2.0_wp*e_kinetic/REAL(simpar%nfree,wp)

!     *** Scale temperature to reference temperture, if requested ***

      IF (ABS(simpar%temp_ext - temperature) > simpar%temp_tol) THEN

        scale = SQRT(simpar%temp_ext/temperature)

        DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               atom_list=atom_list,&
                               mass=mass)
          v2 = 0.0_wp
          natom = SIZE(atom_list)
          DO iatom=1,natom
            atom = atom_list(iatom)
            particle_set(atom)%v(:) = scale*particle_set(atom)%v(:)
          END DO
        END DO

        IF (ionode) THEN
          WRITE (UNIT=output_unit,FMT="(/,T2,A,F10.2,A,F10.2,A)")&
            "Temperature scaled to requested temperature:",&
            temperature*kelvin," K  ->",simpar%temp_ext*kelvin," K"
        END IF

      END IF

    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE temperature_control

! *****************************************************************************

  SUBROUTINE external_control(qs_env,globenv,should_stop)

!   Purpose: External manipulations during a QS MD run.

!   History: - Creation (16.09.2003,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), POINTER        :: qs_env
    LOGICAL, INTENT(OUT)                      :: should_stop

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "external_control",&
                                   exit_file_name = "EXIT_MD"

!   *** Local variables ***

    INTEGER :: group,handle,output_unit,source,unit_number
    LOGICAL :: file_exists,ionode

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr
    source = globenv%source

    should_stop = .FALSE.

    INQUIRE (FILE=exit_file_name,EXIST=file_exists)

    IF (file_exists) THEN

      should_stop = .TRUE.

      CALL mp_bcast(should_stop,source,group)

!     *** Delete exit file ***

      IF (ionode) THEN
        CALL open_file(file_name=exit_file_name,unit_number=unit_number)
        CALL close_file(unit_number=unit_number,file_status="DELETE")
        WRITE (UNIT=output_unit,FMT="(/,T2,A,/)")&
          "*** MD run terminated by external request ***"
      END IF

    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE external_control

! *****************************************************************************

END MODULE qs_md

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
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
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_retain,&
                                             atomic_kind_list_type
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
  USE cp_fragment_types,               ONLY: cp_fragment_p_type,&
                                             cp_fragment_type,&
                                             fragment_get
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dump,                            ONLY: dump_variables
  USE extended_system_dynamics,        ONLY: lnhc_parameters_type,&
                                             npt_info_type
  USE extended_system_types,           ONLY: extended_system_type
  USE force_control,                   ONLY: force_env_calculate_force,&
                                             force_env_create,&
                                             force_env_get,&
                                             force_env_release,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE initialize_extended_types,       ONLY: initialize_nhc_baro,&
                                             initialize_nhc_forces,&
                                             initialize_nhc_part_new,&
                                             initialize_npt_type
  USE input_types,                     ONLY: setup_parameters_type
  USE integrator,                      ONLY: set_integrator
  USE kinds,                           ONLY: int_size,&
                                             dp,&
                                             dp_size
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
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE mol_kind_list_types,             ONLY: mol_kind_list_create,&
                                             mol_kind_list_release,&
                                             mol_kind_list_retain,&
                                             mol_kind_list_type
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
  USE molecule_types,                  ONLY: intra_parameters_type,&
                                             molecule_structure_type,&
                                             molecule_type,&
                                             particle_node_type
  USE output_utilities,                ONLY: print_warning
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_retain,&
                                             particle_list_type
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

  SUBROUTINE qs_mol_dyn ( force_env, globenv, error )

!   Purpose: Driver routine for MD run using QUICKSTEP.

!   History: - Creation (07.11.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), POINTER :: globenv
    TYPE(force_env_type), POINTER :: force_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_mol_dyn",&
         routineP=module_name//':'//routine_name

!   *** Local variables ***

    TYPE(extended_system_type)            :: extended_type
    TYPE(intra_parameters_type)           :: intra_param
    TYPE(md_environment_type)             :: md_env
    TYPE(mdio_parameters_type)            :: mdio
    TYPE(setup_parameters_type)           :: setup
    TYPE(mol_kind_list_type), POINTER     :: molecule_kinds
    TYPE(mol_kind_new_list_type), POINTER :: molecule_kinds_new
    TYPE(particle_list_type), POINTER     :: particles
    TYPE(mol_struct_list_type), POINTER   :: local_molecules

    TYPE(cell_type), POINTER                  :: cell
    TYPE(cell_type), POINTER                  :: cell_ref
    TYPE(simulation_parameters_type), POINTER :: simpar

    TYPE(structure_type), DIMENSION(:), POINTER          :: struc
    TYPE(virial_type), POINTER                           :: virial
    REAL(KIND = dp), POINTER                                    :: constant
    INTEGER, POINTER                                     :: itimes
    REAL(KIND = dp)                                             :: ekin
    REAL(KIND = dp)                                             :: econs 
    INTEGER                                              :: i
    TYPE(cp_fragment_p_type), DIMENSION(:), POINTER      :: subsys
    TYPE ( lnhc_parameters_type ), POINTER :: nhc_part, nhc_baro

    INTEGER :: handle,istat,istep,natom,nmolecule
    INTEGER :: nmoltype, nmol, nnodes, nsh, isubsys
    LOGICAL :: should_stop, failure

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)
    failure=.FALSE.
    NULLIFY ( molecule_kinds, particles, local_molecules, cell, cell_ref,&
         simpar,struc,virial,constant,itimes,subsys,nhc_part, nhc_baro )
    NULLIFY (molecule_kinds_new)

    CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)

    ALLOCATE (simpar,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "simpar",0)
    ALLOCATE (virial,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "virial",0)
    CALL zero_virial(virial)

    CALL read_md_section(simpar,globenv,mdio)
    simpar%program = globenv%program_name

    CALL force_env_get(force_env,subsys=subsys, cell=cell, error=error)

    ALLOCATE (struc(SIZE(subsys)),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    !FM probably in some cases it might also be good to connect the same
    !FM termostat to different subsystems (add a way in the input to specify 
    !FM this?)
    ALLOCATE (extended_type%nhc_part(SIZE(subsys)),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (extended_type%nhc_baro(SIZE(subsys)),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! Define a reference cell for compatibility with MD restarts
    cell_ref => cell

    DO isubsys=1,SIZE(subsys)
       CALL fragment_get(subsys(isubsys)%fragment,&
            local_molecules=local_molecules,&
            molecule_kinds=molecule_kinds,&
            molecule_kinds_new=molecule_kinds_new,&
            error=error)
       CALL initialize_nhc_part_new(cell,simpar,molecule_kinds_new%els,&
                                    extended_type%nhc_part(isubsys),globenv)
       CALL initialize_nhc_baro(simpar,globenv,&
                                extended_type%nhc_baro(isubsys))
    END DO
    CALL initialize_npt_type(simpar,globenv,extended_type%npt_info,cell)

! initialize velocities and read restart files if needed
    SELECT CASE ( simpar % read_type )
    CASE ( 'INIT' )
       DO isubsys=1,SIZE(subsys)
          CALL fragment_get(subsys(isubsys)%fragment,&
               particles=particles,&
               error=error)
          CALL initialize_velocities ( simpar, &
               particles%els, globenv )
       END DO
    CASE ( 'POS' )
       !FM should be adapted to the existence of many fragments
       CPPrecondition(SIZE(subsys)==1,cp_failure_level,routineP,error,failure)
       isubsys=1
       CALL fragment_get(subsys(isubsys)%fragment,&
            particles=particles,&
            error=error)
       nhc_part => extended_type%nhc_part(isubsys)
       nhc_baro => extended_type%nhc_baro(isubsys)
       CALL read_md_restart ( simpar % read_type, globenv,   &
            particles%els, cell, cell_ref,  &
            nhc_part, &
            nhc_baro, extended_type%npt_info )
       CALL initialize_velocities ( simpar, particles%els, globenv )
    CASE DEFAULT
       !FM should be adapted to the existence of many fragments
       CPPrecondition(SIZE(subsys)==1,cp_failure_level,routineP,error,failure)
       isubsys=1
       CALL fragment_get(subsys(isubsys)%fragment,&
            particles=particles,&
            error=error)
       nhc_part => extended_type%nhc_part(isubsys)
       nhc_baro => extended_type%nhc_baro(isubsys)
       CALL read_md_restart ( simpar % read_type, globenv,   &
                              particles%els, &
                              cell, cell_ref,  &
                              nhc_part, &
                              nhc_baro, &
                              extended_type%npt_info )
    END SELECT

!MK    DO isubsys=1,SIZE(subsys)
!MK      nhc_part => extended_type%nhc_part(isubsys)
!MK      nhc_baro => extended_type%nhc_baro(isubsys)
!MK      CALL initialize_nhc_forces (nhc_part,nhc_baro)
!MK    END DO

    DO isubsys=1,SIZE(subsys)
       CALL init_structure_type(struc(isubsys))
       CALL fragment_get(subsys(isubsys)%fragment,&
            pnode=struc(isubsys)%pnode, particles=particles,&
            local_molecules=local_molecules, error=error)
       struc(isubsys)%part => particles%els
!MK       struc(isubsys)%molecule => local_molecules%els
    END DO

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
                                      "constant",dp_size)

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
      CALL temperature_control(simpar,force_env,error)
      CALL md_energy_qs(md_env)
      IF (MODULO(itimes,mdio%idump) == 0) THEN
        CALL dump_variables(md_env,force_env,mdio%dump_file_name,error)
      END IF
      CALL external_control(force_env,should_stop,error=error)
      IF (should_stop) EXIT
    END DO

    CALL destroy_md_env(md_env)

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE qs_mol_dyn

! *****************************************************************************

  SUBROUTINE temperature_control(simpar,force_env,error)

    !   Purpose: Perform all temperature manipulations during a QS MD run.

    !   History: - Creation (15.09.2003,MK)
    !            - adapted to force_env (05.10.2003,fawzi)

    !   ***************************************************************************

    TYPE(simulation_parameters_type), POINTER :: simpar
    TYPE(force_env_type), POINTER :: force_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    !   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "temperature_control",&
         routineP=module_name//':'//routine_name

    !   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    TYPE(atomic_kind_list_type), POINTER :: atomic_kinds
    TYPE(particle_list_type), POINTER :: particles
    TYPE(distribution_1d_type), POINTER :: local_particles
    INTEGER, DIMENSION(:), POINTER :: atom_list

    REAL(KIND = dp) :: e_kinetic,mass,temperature,scale,v2
    INTEGER  :: atom,handle,iatom,ikind,istat,natom,nkind, isubsys
    LOGICAL  :: failure

    TYPE(cp_fragment_p_type), DIMENSION(:), POINTER :: subsys
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(cp_para_env_type), POINTER :: para_env

    !   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)
    NULLIFY(atomic_kind, atomic_kinds, particles, local_particles, atom_list,&
         subsys, para_env)
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(simpar),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)

    IF (simpar%temp_tol > 0.0_dp) THEN

       CALL force_env_get(force_env,subsys=subsys, para_env=para_env,error=error)

       DO isubsys=1,SIZE(subsys)

          CALL fragment_get(subsys(isubsys)%fragment,&
               local_particles=local_particles,&
               particles=particles, atomic_kinds=atomic_kinds,&
               error=error)

          !     *** Calculate the current temperature ***
          e_kinetic = 0.0_dp

          DO ikind=1,SIZE(local_particles%list)
             atomic_kind => atomic_kinds%els(ikind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  mass=mass)
             atom_list => local_particles%list(ikind)%array
             v2 = 0.0_dp
             DO iatom=1,SIZE(atom_list)
                atom = atom_list(iatom)
                v2 = v2 + particles%els(atom)%v(1)*particles%els(atom)%v(1) +&
                     particles%els(atom)%v(2)*particles%els(atom)%v(2) +&
                     particles%els(atom)%v(3)*particles%els(atom)%v(3)
             END DO
             e_kinetic = e_kinetic + 0.5_dp*mass*v2
          END DO
       END DO
       CALL mp_sum(e_kinetic,para_env%group) !FM loop on all particles and avoid sum?

       temperature = 2.0_dp*e_kinetic/REAL(simpar%nfree,dp)

       !     *** Scale temperature to reference temperture, if requested ***

       IF (ABS(simpar%temp_ext - temperature) > simpar%temp_tol) THEN

          scale = SQRT(simpar%temp_ext/temperature)

          DO isubsys=1,SIZE(subsys)

             CALL fragment_get(subsys(isubsys)%fragment,&
                  local_particles=local_particles,&
                  particles=particles, atomic_kinds=atomic_kinds,&
                  error=error)

             DO ikind=1,atomic_kinds%n_els
                atomic_kind => atomic_kinds%els(ikind)
                CALL get_atomic_kind(atomic_kind=atomic_kind,&
                     atom_list=atom_list,&
                     mass=mass)
                v2 = 0.0_dp
                natom = SIZE(atom_list)
                DO iatom=1,natom
                   atom = atom_list(iatom)
                   particles%els(atom)%v(:) = scale*particles%els(atom)%v(:)
                END DO
             END DO

          END DO

          IF (logger%para_env%mepos==logger%para_env%source) THEN
             WRITE (UNIT=cp_logger_get_default_unit_nr(logger),&
                  FMT="(/,T2,A,F10.2,A,F10.2,A)")&
                  "Temperature scaled to requested temperature:",&
                  temperature*kelvin," K  ->",simpar%temp_ext*kelvin," K"
          END IF

       END IF

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE temperature_control

! *****************************************************************************

  SUBROUTINE external_control(force_env,should_stop,error)

!   Purpose: External manipulations during an MD run.

!   History: - Creation (16.09.2003,MK)

!   ***************************************************************************

    TYPE(force_env_type), POINTER :: force_env
    LOGICAL, INTENT(OUT)                         :: should_stop
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "external_control",&
                                   exit_file_name = "EXIT_MD",&
                                   routineP=module_name//':'//routine_name

!   *** Local variables ***

    INTEGER :: handle,unit_number
    LOGICAL :: file_exists, failure
    TYPE(cp_logger_type), POINTER :: logger

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)

    logger => cp_error_get_logger(error)
    should_stop = .FALSE.

    IF (logger%para_env%source==logger%para_env%mepos) THEN
       INQUIRE (FILE=exit_file_name,EXIST=should_stop)

       IF (should_stop) THEN
          CALL open_file(file_name=exit_file_name,unit_number=unit_number)
          CALL close_file(unit_number=unit_number,file_status="DELETE")
          WRITE (UNIT=cp_logger_get_default_unit_nr(logger),FMT="(/,T2,A,/)")&
               "*** MD run terminated by external request ***"
       END IF
    END IF
    CALL mp_bcast(should_stop,logger%para_env%source,logger%para_env%group)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE external_control

! *****************************************************************************

END MODULE qs_md

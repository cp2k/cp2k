!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/md_run [1.0] *
!!
!!   NAME
!!     md_run
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

MODULE md_run
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE atoms_input,                     ONLY: read_md_restart
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dump,                            ONLY: dump_variables
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_initialize,&
                                             dyn_coeff_set_p_type
  USE extended_system_types,           ONLY: extended_system_type,&
                                             lnhc_parameters_type
  USE force_env_methods,               ONLY: force_env_calculate_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             force_env_p_type
  USE global_types,                    ONLY: global_environment_type
  USE initialize_extended_types,       ONLY: initialize_nhc_baro,&
                                             initialize_nhc_coef,&
                                             initialize_nhc_forces,&
                                             initialize_nhc_part,&
                                             initialize_npt_type
  USE integrator,                      ONLY: set_integrator
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE md,                              ONLY: initialize_velocities,&
                                             mdio_parameters_type,&
                                             read_md_section,&
                                             simulation_parameters_type
  USE md_energies,                     ONLY: md_energy,&
                                             set_energies
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_env_create,&
                                             md_env_release,&
                                             md_environment_type,&
                                             set_md_env
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: read_fixed_particles
  USE physcon,                         ONLY: kelvin
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE velocity_verlet_control,         ONLY: velocity_verlet

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "md_run"

  PUBLIC :: qs_mol_dyn

  INTERFACE qs_mol_dyn
     MODULE PROCEDURE qs_mol_dyn1, qs_mol_dyn2
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE qs_mol_dyn1 ( force_env, globenv, error )

    TYPE(force_env_type), POINTER  :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    TYPE(force_env_p_type), DIMENSION(:), POINTER  :: force_env_p
    INTEGER :: stat
    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_mol_dyn1", &
      routineP = module_name//':'//routine_name

    ALLOCATE(force_env_p(1), stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    force_env_p(1)%force_env => force_env
    CALL qs_mol_dyn2 ( force_env_p, globenv, error )
  
  END SUBROUTINE qs_mol_dyn1

! *****************************************************************************

  SUBROUTINE qs_mol_dyn2 ( force_env_p, globenv, error )

!   Purpose: Driver routine for MD run using QUICKSTEP.

!   History: - Creation (07.11.2002,MK)

!   ***************************************************************************

    TYPE(force_env_p_type), DIMENSION(:), POINTER  :: force_env_p
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_mol_dyn2", &
      routineP = module_name//':'//routine_name

    INTEGER                                  :: handle, i, iset, istat, &
                                                istep, isubsys, ibead, nbeads,&
                                                nfixed
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: failure, should_stop
    REAL(KIND=dp)                            :: md_run_start_time, &
                                                time_iter_start, &
                                                time_iter_stop
    REAL(KIND=dp), POINTER                   :: used_time
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), POINTER                  :: dyn_coeff_set
    TYPE(extended_system_type)               :: extended_type
    TYPE(lnhc_parameters_type), POINTER      :: nhc_baro, nhc_part
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(mdio_parameters_type)               :: mdio
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    ! will be used to stop if a maximum walltime is given
    md_run_start_time=m_walltime()

    failure=.FALSE.
    NULLIFY ( particles, local_molecules, cell, cell_ref, dyn_coeff_set,&
              simpar,itimes,used_time,subsys,nhc_part,nhc_baro)
    NULLIFY (md_env,molecule_kinds_new, molecules_new)

    CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env_p),cp_failure_level,routineP,error,failure)
    nbeads = SIZE ( force_env_p )
    DO i = 1, nbeads
      CPPrecondition ( ASSOCIATED ( force_env_p ( i ) % force_env ), cp_failure_level,routineP,error,failure)
    END DO

    ALLOCATE (simpar,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "simpar",0)
    CALL read_md_section ( simpar, globenv, mdio )
    simpar%program = globenv%program_name

!    CALL force_env_get(force_env,subsys=subsys, cell=cell,&
!                       dyn_coeff_set=dyn_coeff_set, error=error)
!    !FM probably in some cases it might also be good to connect the same
!    !FM termostat to different subsystems (add a way in the input to specify 
!    !FM this?)
!    ALLOCATE (extended_type%nhc_part(SIZE(subsys)),STAT=istat)
    ALLOCATE (extended_type%nhc_part(nbeads),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (extended_type%nhc_baro(nbeads),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

! loop over number of beads
    DO ibead = 1, nbeads

      CALL force_env_get( force_env_p(ibead)%force_env, subsys=subsys, &
                          cell=cell, error=error )

! Define a reference cell for compatibility with MD restarts
      cell_ref => cell

      DO isubsys=1,SIZE(subsys)
         CALL cp_subsys_get(subsys(isubsys)%subsys,&
              local_molecules_new=local_molecules,&
              molecule_kinds_new=molecule_kinds_new,&
              molecules_new=molecules_new,&
              particles=particles,&
              error=error)
          IF (simpar%fix_atoms) THEN
            CALL read_fixed_particles(particles%els,globenv,nfixed)
          END IF
         CALL initialize_nhc_part(cell,simpar,local_molecules, &
                                  molecules_new%els, molecule_kinds_new%els,&
                                  globenv, extended_type%nhc_part(ibead))
         CALL initialize_nhc_baro(simpar,globenv, extended_type%nhc_baro(ibead))
      END DO
      CALL initialize_npt_type(simpar,globenv,extended_type%npt_info,cell)

! initialize coef thermostats and setup coefficient forces, velocities and effective masses
    IF (ASSOCIATED(dyn_coeff_set)) THEN
       ALLOCATE (extended_type%nhc_coef(SIZE(dyn_coeff_set)),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO iset=1,SIZE(dyn_coeff_set)
          CALL initialize_nhc_coef(simpar,dyn_coeff_set(iset)%dyn_coeff_set,&
                                   globenv,extended_type%nhc_coef(iset))
          CALL dyn_coeff_set_initialize(simpar,dyn_coeff_set(iset)%dyn_coeff_set,globenv) 
          CALL initialize_nhc_forces (nhc_coef=extended_type% nhc_coef(iset))
       END DO
    END IF

! initialize velocities and read restart files if needed
      SELECT CASE ( simpar % read_type )
      CASE ( 'INIT' )
         DO isubsys=1,SIZE(subsys)
            CALL cp_subsys_get(subsys(isubsys)%subsys,&
                 particles=particles,&
                 error=error)
            CALL initialize_velocities ( simpar, &
                 particles%els, globenv )
         END DO
      CASE ( 'POS' )
!FM should be adapted to the existence of many fragments
         CPPrecondition(SIZE(subsys)==1,cp_failure_level,routineP,error,failure)
         isubsys=1
         CALL cp_subsys_get(subsys(isubsys)%subsys,&
              particles=particles,&
              error=error)
         nhc_part => extended_type%nhc_part(ibead)
         nhc_baro => extended_type%nhc_baro(ibead)
         CALL read_md_restart ( simpar % read_type, globenv,   &
              particles%els, cell, cell_ref,  &
              nhc_part, nhc_baro, extended_type%npt_info )
         CALL initialize_velocities ( simpar, particles%els, globenv )
      CASE DEFAULT
!FM should be adapted to the existence of many fragments
         CPPrecondition(SIZE(subsys)==1,cp_failure_level,routineP,error,failure)
          isubsys=1
          CALL cp_subsys_get(subsys(isubsys)%subsys,&
               particles=particles,&
               error=error)
          nhc_part => extended_type%nhc_part(ibead)
          nhc_baro => extended_type%nhc_baro(ibead)
          CALL read_md_restart ( simpar % read_type, globenv,   &
                              particles%els, &
                              cell, cell_ref,  &
                              nhc_part, &
                              nhc_baro, &
                              extended_type%npt_info )
       END SELECT

       nhc_part => extended_type%nhc_part(ibead)
       nhc_baro => extended_type%nhc_baro(ibead)
       CALL initialize_nhc_forces ( nhc_part=nhc_part, nhc_baro=nhc_baro )

       DO isubsys=1,SIZE(subsys)
          CALL cp_subsys_get(subsys(isubsys)%subsys,&
                             particles=particles,&
                             local_molecules_new=local_molecules, &
                             error=error)
       END DO
! end loop over number of beads
    END DO 

    CALL md_env_create ( md_env, globenv % para_env, error = error )

    CALL set_md_env(md_env=md_env,&
                    simpar=simpar,&
                    nhc_part=extended_type%nhc_part,&
                    nhc_baro=extended_type%nhc_baro,&
                    npt=extended_type%npt_info,&
                    cell=cell,&
                    force_env_p=force_env_p)

    IF (ASSOCIATED(dyn_coeff_set)) THEN
       CALL set_md_env( md_env, nhc_coef=extended_type%nhc_coef)
    END IF

    CALL set_integrator ( globenv )

    CALL set_energies ( globenv, mdio )
    ! the verlet integrator needs to be started with the correct forces alreay initialized

    time_iter_start=m_walltime()

    CALL get_md_env ( md_env, force_env_p=force_env_p, itimes = itimes, &
                      used_time = used_time ) 
    DO ibead = 1, nbeads
       CALL force_env_calculate_force ( force_env_p(ibead)%force_env )
    END DO
    CALL md_energy ( md_env ) ! inits the print env at itimes == 0

    DO istep=1,simpar%nsteps
      CALL velocity_verlet(md_env)
      itimes = itimes + 1
      do ibead=1,nbeads
        CALL temperature_control(simpar,force_env_p(ibead)%force_env,error)
      end do
      time_iter_stop=m_walltime()
      used_time = time_iter_stop - time_iter_start
      time_iter_start=time_iter_stop
      CALL md_energy ( md_env )
      IF (MODULO(itimes,mdio%idump) == 0) THEN
        CALL dump_variables(md_env,mdio%dump_file_name,error)
      END IF

      ! test for EXIT_MD to exit
      CALL external_control(force_env_p,should_stop,error=error)
      IF (should_stop) EXIT

      ! test if we have used more walltime than allowed      
      IF (simpar%walltime > 0.0_dp .AND. simpar%walltime < m_walltime()-md_run_start_time) THEN
          should_stop=.TRUE.
      ENDIF
      CALL mp_bcast(should_stop,globenv%para_env%source,globenv%para_env%group) ! sync values of should_stop
      IF (should_stop) EXIT
    END DO 

    CALL md_env_release ( md_env, error )

    ! just deallocate this stuff ?? 
    ! seem to need proper init / destroy subroutines ... to be fixed
    DEALLOCATE(extended_type%nhc_part,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(extended_type%nhc_baro,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(simpar,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE qs_mol_dyn2

! *****************************************************************************

  SUBROUTINE temperature_control(simpar,force_env,error)

!   Purpose: Perform all temperature manipulations during a QS MD run.

!   History: - Creation (15.09.2003,MK)
!            - adapted to force_env (05.10.2003,fawzi)

!   ***************************************************************************

    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routine_name = "temperature_control", &
      routineP = module_name//':'//routine_name

    INTEGER                                  :: atom, handle, iatom, ikind, &
                                                isubsys, natom
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: e_kinetic, mass, scale, &
                                                temperature, v2
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_list_type), POINTER        :: particles

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

          CALL cp_subsys_get(subsys(isubsys)%subsys,&
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

             CALL cp_subsys_get(subsys(isubsys)%subsys,&
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

  SUBROUTINE external_control(force_env_p,should_stop,error)

!   Purpose: External manipulations during an MD run.

!   History: - Creation (16.09.2003,MK)

!   ***************************************************************************

    TYPE(force_env_p_type), DIMENSION(:), POINTER  :: force_env_p
    LOGICAL, INTENT(OUT)                     :: should_stop
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: exit_file_name = "EXIT_MD", &
      routine_name = "external_control", &
      routineP = module_name//':'//routine_name

    INTEGER                                  :: handle, unit_number
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(force_env_p),cp_failure_level,routineP,error,failure)

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

END MODULE md_run

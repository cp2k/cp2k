!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
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
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_initialize,&
                                             dyn_coeff_set_p_type
  USE extended_system_types,           ONLY: extended_system_type,&
                                             extended_type_init,&
                                             lnhc_dealloc,&
                                             lnhc_init,&
                                             npt_dealloc
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE init_extended_system_types,      ONLY: initialize_nhc_baro,&
                                             initialize_nhc_coef,&
                                             initialize_nhc_part,&
                                             initialize_npt
  USE init_extended_system_variables,  ONLY: init_nhc_forces
  USE input_constants,                 ONLY: langevin_ensemble
  USE input_cp2k_restarts,             ONLY: write_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE md,                              ONLY: initialize_velocities,&
                                             read_md_section,&
                                             simulation_parameters_type
  USE md_energies,                     ONLY: md_energy
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_env_create,&
                                             md_env_release,&
                                             md_environment_type,&
                                             set_md_env
  USE message_passing,                 ONLY: mp_sum
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: kelvin
  USE termination,                     ONLY: external_control,&
                                             stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE velocity_verlet_control,         ONLY: velocity_verlet
  USE virial_methods,                  ONLY: virial_evaluate
!JCS
  USE virial_types,                    ONLY: virial_set,&
                                             virial_type
  USE wiener_process,                  ONLY: create_wiener_process
!JCS 
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_run'

  PUBLIC :: qs_mol_dyn, &
            temperature_control


  INTERFACE qs_mol_dyn
     MODULE PROCEDURE qs_mol_dyn1, qs_mol_dyn2
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE qs_mol_dyn1 ( force_env, globenv, error )

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_mol_dyn1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p

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

    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_mol_dyn2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ibead, in_use, &
                                                istat, istep, isubsys, &
                                                iter_level, nbeads
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: failure, should_stop
    REAL(KIND=dp)                            :: time_iter_start, &
                                                time_iter_stop
    REAL(KIND=dp), POINTER                   :: t, used_time
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), POINTER                  :: dyn_coeff_set
    TYPE(extended_system_type)               :: extended_type
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                motion_section, &
                                                subsys_section, work_section
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(virial_type), POINTER               :: virial

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    NULLIFY ( particles, local_molecules, cell, cell_ref, dyn_coeff_set,&
              simpar,itimes,used_time,subsys,motion_section )
    NULLIFY (md_env,molecule_kinds_new, molecules_new, force_env_section,&
         work_section, atomic_kinds, local_particles, t)
    logger => cp_error_get_logger(error)

    force_env_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL",error=error)
    motion_section => section_vals_get_subs_vals(globenv%input_file,"MOTION",error=error)
    subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)
    CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env_p),cp_failure_level,routineP,error,failure)
    nbeads = SIZE ( force_env_p )
    DO i = 1, nbeads
      CPPrecondition ( ASSOCIATED ( force_env_p ( i ) % force_env ), cp_failure_level,routineP,error,failure)
    END DO

    ALLOCATE (simpar,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "simpar",0)

    CALL cp_add_iter_level(logger%iter_info,"MD",n_rlevel_new=iter_level,&
         error=error)
    CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=0,error=error)    
    CALL read_md_section ( simpar, globenv, error )
    
    CALL extended_type_init(extended_type)
    ALLOCATE (extended_type%nhc_part(nbeads),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (extended_type%nhc_baro(nbeads),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (extended_type%nhc_coef(nbeads),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ibead=1,nbeads
       CALL lnhc_init(extended_type%nhc_part(ibead))
       CALL lnhc_init(extended_type%nhc_baro(ibead))
       CALL lnhc_init(extended_type%nhc_coef(ibead))
    ENDDO

    DO ibead = 1, nbeads

       CALL force_env_get( force_env_p(ibead)%force_env, subsys=subsys, &
                          dyn_coeff_set=dyn_coeff_set,in_use=in_use,&
                          cell=cell, cell_ref = cell_ref, virial = virial,&
                          error=error )

! Define a reference cell for compatibility with MD restarts

       DO isubsys=1,SIZE(subsys)
          CALL cp_subsys_get(subsys(isubsys)%subsys,&
              local_molecules_new=local_molecules,&
              molecule_kinds_new=molecule_kinds_new,&
              atomic_kinds=atomic_kinds,&
              local_particles=local_particles,&
              molecules_new=molecules_new,&
              particles=particles,&
              error=error)
          !
          ! Possibly restart Nose particles
          !
          NULLIFY(work_section)
          work_section => section_vals_get_subs_vals(motion_section,"MD%NOSE",error=error)
          CALL initialize_nhc_part(cell,simpar,local_molecules, &
                                   molecules_new%els, molecule_kinds_new%els,&
                                   globenv, extended_type%nhc_part(ibead),work_section=work_section,error=error)
          !
          ! Possibly restart barostat
          !
          NULLIFY(work_section)
          work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%NOSE",error=error)
          CALL initialize_nhc_baro(simpar,globenv, extended_type%nhc_baro(ibead),work_section=work_section,error=error)

          IF (ASSOCIATED(dyn_coeff_set)) THEN
             !
             ! Initialize coef thermostats and setup coefficient forces, velocities and effective masses
             ! Possibly restart KG coefficients
             !
             NULLIFY(work_section)
             work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%NOSE",error=error)
             CALL initialize_nhc_coef(dyn_coeff_set(isubsys)%dyn_coeff_set,&
                                   globenv,extended_type%nhc_coef(ibead),work_section=work_section,error=error)

             NULLIFY(work_section)
             work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF",error=error)
             CALL dyn_coeff_set_initialize(dyn_coeff_set(isubsys)%dyn_coeff_set,globenv,&
                  work_section=work_section,error=error)
             CALL init_nhc_forces (nhc_coef=extended_type% nhc_coef(ibead),error=error)
          END IF
       END DO
       !
       ! Possibly restart Barostat
       !
       NULLIFY(work_section)
       work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT",error=error)
       CALL initialize_npt (simpar,globenv,extended_type%npt_info,cell,work_section=work_section,error=error)

       CALL virial_set ( virial,pv_availability=simpar%pv_availability,&
                         pv_numer=simpar%pv_numerical)
       !
       ! Always initialize velocities 
       ! and possibly restart them
       !
       DO isubsys=1,SIZE(subsys)
          CALL cp_subsys_get(subsys(isubsys)%subsys,&
              particles=particles,&
              error=error)
          NULLIFY(work_section)
          work_section => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)       
          CALL initialize_velocities ( simpar, particles%els, globenv, work_section=work_section,error=error)
       END DO

       ! Set V0 if needed
       IF ( ASSOCIATED ( cell_ref ) ) THEN
          IF ( simpar % v0 == 0._dp ) simpar % v0 = cell_ref % deth
       ENDIF
       
       CALL init_nhc_forces ( nhc_part=extended_type%nhc_part(ibead),error=error)
       CALL init_nhc_forces ( nhc_baro=extended_type%nhc_baro(ibead),error=error)

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

    ! loop over number of beads
    IF (simpar%ensemble == langevin_ensemble) THEN
      CALL create_wiener_process(md_env,subsys_section,error)
    END IF

    time_iter_start=m_walltime()

    CALL get_md_env ( md_env, force_env_p=force_env_p, itimes = itimes, &
                      t = t, used_time = used_time ) 

    DO ibead = 1, nbeads
       CALL force_env_calc_energy_force ( force_env_p(ibead)%force_env,&
            calc_force=.TRUE.,error=error)

       CALL force_env_get( force_env_p(ibead)%force_env, subsys=subsys, &
                           error=error )

       CALL cp_subsys_get(subsys(1)%subsys,&
              atomic_kinds=atomic_kinds,&
              local_particles=local_particles,&
              particles=particles,&
              error=error)

       CALL virial_evaluate(atomic_kinds%els, particles%els, &
            local_particles, force_env_p(ibead)%force_env% virial,&
            force_env_p(ibead)%force_env%para_env%group)

    END DO

    CALL section_vals_val_get(motion_section,"MD%STEP_START_VAL",i_val=itimes,error=error)
    CALL section_vals_val_get(motion_section,"MD%TIME_START_VAL",r_val=t,error=error)
    CALL md_energy ( md_env, globenv, error ) ! inits the print env at itimes == 0

    DO istep=1,simpar%nsteps
      itimes = itimes + 1
      CALL cp_iterate(logger%iter_info,last=(istep==simpar%nsteps),&
           iter_nr=itimes,error=error)
      CALL velocity_verlet(md_env,error)
      ! test for EXIT_MD to exit, if so we don't overwrite the restart or append to the trajectory
      CALL external_control(should_stop,"MD",globenv,error=error)
      IF (should_stop) EXIT
      DO ibead=1,nbeads
        CALL temperature_control(simpar,force_env_p(ibead)%force_env,error)
      END DO
      time_iter_stop=m_walltime()
      used_time = time_iter_stop - time_iter_start
      time_iter_start=time_iter_stop

      CALL md_output(md_env,force_env_section,motion_section,globenv,error)

    END DO 

    CALL cp_rm_iter_level(logger%iter_info,n_rlevel_att=iter_level,&
         error=error)

    CALL md_env_release ( md_env, error )

    ! just deallocate this stuff ?? 
    ! seem to need proper init / destroy subroutines ... to be fixed
    CALL lnhc_dealloc(extended_type%nhc_part)
    CALL lnhc_dealloc(extended_type%nhc_baro)
    CALL lnhc_dealloc(extended_type%nhc_coef)
    CALL npt_dealloc(extended_type%npt_info)
    DEALLOCATE(simpar,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_mol_dyn2

!!****f* md_run/md_output *
!!
!!   NAME
!!     md_output 
!!
!!   FUNCTION
!!     collects the part of the MD that, basically, does the output
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE md_output(md_env,force_env_section,motion_section,globenv,error)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                motion_section
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'md_output', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)


    CALL md_energy ( md_env, globenv, error ) ! also writes trajectories

    CALL write_restart(md_env=md_env,globenv=globenv,force_env_section=force_env_section,&
                 error=error)

    CALL timestop(handle)

END SUBROUTINE md_output

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

    CHARACTER(LEN=*), PARAMETER :: routineN = 'temperature_control', &
      routineP = moduleN//':'//routineN

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

    CALL timeset(routineN,"I","",handle)
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

END MODULE md_run

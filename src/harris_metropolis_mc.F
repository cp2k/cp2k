!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2005  CP2K developers group                            !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/harris_metropolis_mc [1.0] *
!!
!!   NAME
!!     harris_metropolis_mc
!!
!!   FUNCTION
!!     Metropolis MC Langevin dynamics with the Harris functional
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     09.2005 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE harris_metropolis_mc

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE atoms_input,                     ONLY: read_md_restart
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_error_handling,               ONLY: cp_a_l, &
                                             cp_assert, &
                                             cp_assertion_failed, &
                                             cp_debug, &
                                             cp_error_get_logger, &
                                             cp_error_message, &
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level, &
                                             cp_note_level, &
                                             cp_to_string, &
                                             cp_warning_level, &
                                             cp_logger_get_default_unit_nr, &
                                             cp_logger_type
  USE cp_output_handling,              ONLY: cp_add_iter_level, &
                                             cp_rm_iter_level, &
                                             cp_iterate
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type, &
                                             cp_subsys_get
  USE cp2k_units,                      ONLY: get_cp2k_units, &
                                             unit_convert_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dump,                            ONLY: dump_variables
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_p_type, &
                                             dyn_coeff_set_initialize
  USE extended_system_types,           ONLY: extended_system_type, &
                                             extended_type_init, &
                                             lnhc_parameters_type, &
                                             lnhc_init, &
                                             lnhc_dealloc, &
                                             npt_dealloc
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_p_type, &
                                             force_env_type, &
                                             force_env_get
  USE global_types,                    ONLY: global_environment_type
  USE initialize_extended_types,       ONLY: initialize_nhc_forces, &
                                             initialize_nhc_part, &
                                             initialize_nhc_baro, &
                                             initialize_nhc_coef, &
                                             initialize_npt_type
  USE input_constants,                 ONLY: langevin_ensemble, &
                                             do_qs, &
                                             do_kg, &
                                             md_pos_restart, &
                                             md_pos_vel_restart, &
                                             md_all_restart
  USE kinds,                           ONLY: dp, &
                                             dp_size
  USE machine,                         ONLY: m_walltime
  USE md,                              ONLY: mdio_parameters_type, &
                                             !read_md_new_section, &
                                             read_md_section, &
                                             simulation_parameters_type, &
                                             initialize_velocities
  USE md_conserved_quantities,         ONLY: compute_conserved_quantity
  USE md_energies,                     ONLY: md_energy
  USE md_environment_types,            ONLY: md_environment_type, &
                                             md_env_create, &
                                             md_env_release, &
                                             get_md_env, &
                                             set_md_env
  USE md_run,                          ONLY: temperature_control
  USE message_passing,                 ONLY: mp_bcast, &
                                             mp_sum
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE parallel_rng_types,              ONLY: rng_stream_type, &
                                             create_rng_stream, &
                                             next_random_number, &
                                             delete_rng_stream
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: boltzmann
  USE qs_wf_history_methods,           ONLY: wfi_duplicate_history, &
                                             wfi_copy_history
  USE qs_wf_history_types,             ONLY: qs_wf_history_type, &
                                             wfi_release
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory, &
                                             external_control
  USE timings,                         ONLY: timeset, &
                                             timestop
  USE util,                            ONLY: gasdev
  USE velocity_verlet_control,         ONLY: velocity_verlet
  USE wiener_process,                  ONLY: create_wiener_process,&
                                             read_md_restart_file,&
                                             write_md_restart_file

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: moduleN = 'harris_metropolis_mc'

  ! *** Data types ***
  TYPE reduced_particle_type
    REAL(KIND=dp), DIMENSION(3) :: r, &
                                   v
  END TYPE reduced_particle_type


  ! *** Public subroutines ***
  PUBLIC :: harris_mc_ld

!***
!****************************************************************************

CONTAINS

!!****f* harris_metropolis_mc/harris_mc_ld [1.0] *
!!
!!   NAME
!!     harris_mc_ld
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - force_env: The force environment, which is actual a wrapper to abstract
!!                  the various force evaluation methods.
!!     - globenv: The global environment which contains the initially parsed
!!                file and the parallel environment.
!!     - error: variable to control error logging, stopping,... 
!!              see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     09.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_mc_ld(force_env, globenv, error)

    TYPE(force_env_type), POINTER                   :: force_env
    TYPE(global_environment_type), POINTER          :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_mc_ld', &
      routineP = moduleN//':'//routineN

    INTEGER                                         :: stat
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                       :: force_env_p

    ALLOCATE(force_env_p(1), stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    force_env_p(1)%force_env => force_env
    CALL harris_mc_md(force_env_p, globenv, total_energy=.TRUE., &
                      external_temp=.TRUE., qs_accept=.TRUE., &
                      random_velocities=.TRUE., velocity_scaling=.FALSE., &
                      qs_start=.TRUE., restore_wf_history=.TRUE., error=error)

  END SUBROUTINE harris_mc_ld
!***************************************************************************

!!****f* harris_metropolis_mc/harris_mc_md [1.0] *
!!
!!   NAME
!!     harris_mc_md
!!
!!   FUNCTION
!!     Driver routine for the Langevin Dynamics Monte Carlo simulation
!!
!!   NOTES
!!     By using the default values a hybrid Monte Carlo (HMC) simulation is 
!!     performed.
!!     For a high-level quick-start to following references are suggested:
!!
!!     S.Duane, A.D. Kennedy, et.al.: Hybrid Monte Carlo
!!     A.Brass, B.J. Pendleton, et.al.: Hybrid Monte Carlo Simulations Theory 
!!                                      and Initial Comparison with Molecular 
!!                                      Dynamics
!!
!!   INPUTS
!!     - force_env: The force environment, which is actual a wrapper to abstract
!!                  the various force evaluation methods.
!!     - globenv: The global environment which contains the initially parsed
!!                file and the parallel environment.
!!     - total_energy: If TRUE the differences of the total energies are used, 
!!                     instead of the potential energies.
!!     - external_temp: If TRUE the external temperature is used, instead of
!!                      the computed current temperature.
!!     - qs_accept: If TRUE the QS energies are used for the MC acceptance 
!!                  criterion, otherwise Harris or KG energies are used.
!!     - random_velocities: Activates random velocities to perform hybrid MC 
!!                          simulations.
!!     - error: variable to control error logging, stopping,... 
!!              see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     09.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_mc_md(force_env_p, globenv, total_energy, external_temp, &
                          qs_accept, random_velocities, velocity_scaling, &
                          qs_start, restore_wf_history, error)

    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(IN), OPTIONAL            :: total_energy, &
                                                external_temp, &
                                                qs_accept, &
                                                random_velocities, &
                                                velocity_scaling, &
                                                qs_start, &
                                                restore_wf_history
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "harris_mc_md", &
      routineP = moduleN//':'//routineN

    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER                                  :: handle, &
                                                nbeads, &
                                                i, &
                                                istat, &
                                                ibead, &
                                                isubsys, &
                                                iter_level, &
                                                istep, &
                                                iparticles, &
                                                ikind, &
                                                atom, &
                                                iatom, &
                                                iw, &
                                                max_scf_tmp, &
                                                reject, &
                                                accept, &
                                                in_use
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: failure, &
                                                accepted, &
                                                tot_energy, &
                                                ext_temp, &
                                                ionode, &
                                                qs_ener, &
                                                use_harris_tmp, &
                                                rand_vel, &
                                                vel_scaling, &
                                                start_qs, &
                                                restore_wf_past
    REAL(KIND=dp)                            :: md_run_start_time, &
                                                time_iter_start, &
                                                time_iter_stop, &
                                                prev_ldmc_energy, &
                                                new_ldmc_energy, &
                                                random_number, &
                                                prob, &
                                                temperature, &
                                                e_kinetic, &
                                                v2, &
                                                mass, &
                                                Ekin, &
                                                Epot, &
                                                Econs, &
                                                akin, &
                                                temp, &
                                                denom, &
                                                sigma
    REAL(KIND=dp), DIMENSION(3)              :: vcom, &
                                                vcom_new
    REAL(KIND=dp), POINTER                   :: used_time
    TYPE(extended_system_type)               :: extended_type
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(cell_type), POINTER                 :: cell, &
                                                cell_ref
    TYPE(lnhc_parameters_type), POINTER      :: nhc_baro, &
                                                nhc_part
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), POINTER                  :: dyn_coeff_set
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(mdio_parameters_type), &
      POINTER                                :: mdio
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(reduced_particle_type), &
      DIMENSION(:), POINTER                  :: tmp_particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    TYPE(unit_convert_type), POINTER         :: unit
    TYPE(qs_wf_history_type), POINTER        :: tmp_wf_history

!   ------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routineN, globenv)

    CALL timeset(routineN, "I", "", handle)

    ! will be used to stop if a maximum walltime is given
    md_run_start_time=m_walltime()

    IF (PRESENT(total_energy)) THEN
      tot_energy = total_energy
    ELSE
      tot_energy = .TRUE.
    END IF

    IF (PRESENT(external_temp)) THEN
      ext_temp = external_temp
    ELSE
      ext_temp = .TRUE.
    END IF

    IF (PRESENT(qs_accept)) THEN
      qs_ener = qs_accept
    ELSE
      qs_ener = .TRUE.
    END IF

    IF (PRESENT(random_velocities)) THEN
      rand_vel = random_velocities
    ELSE
      rand_vel = .TRUE.
    END IF

    IF (PRESENT(velocity_scaling) .AND. rand_vel) THEN
      vel_scaling = velocity_scaling
    ELSE
      vel_scaling = .FALSE.
    END IF

    IF (PRESENT(qs_start)) THEN
      start_qs = qs_start
    ELSE
      start_qs = .TRUE.
    END IF

    IF (PRESENT(restore_wf_history)) THEN
      restore_wf_past = restore_wf_history
    ELSE
      restore_wf_past = .FALSE.
    END IF

    failure=.FALSE.
    accepted = .TRUE.
    accept = 0
    reject = 0

    NULLIFY(logger, simpar, mdio, subsys, dyn_coeff_set, md_env, nhc_baro, &
            nhc_part, molecule_kinds_new, molecules_new, local_molecules, &
            particles, itimes, used_time, rng_stream, local_particles, atomic_kind, &
            atomic_kinds, atom_list, unit, tmp_wf_history)
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(globenv), cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(force_env_p), cp_failure_level, routineP, error, failure)

    CALL create_rng_stream(rng_stream=rng_stream, name="MonteCarlo acceptance", &
                           extended_precision=.TRUE., error=error)

    nbeads = SIZE(force_env_p)
    DO i = 1, nbeads
      CPPrecondition(ASSOCIATED(force_env_p(i)%force_env), cp_failure_level, routineP, error, failure)
    END DO

    ALLOCATE (simpar,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "simpar", 0)
    ALLOCATE (mdio,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "mdio", 0)

    CALL read_md_section(simpar, globenv, mdio, error)
    !CALL read_md_new_section(simpar, globenv, md_section, error)

    CALL extended_type_init(extended_type)
    ALLOCATE(extended_type%nhc_part(nbeads), STAT=istat)
    CPPostcondition(istat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(extended_type%nhc_baro(nbeads), STAT=istat)
    CPPostcondition(istat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(extended_type%nhc_coef(nbeads), STAT=istat)
    CPPostcondition(istat==0, cp_failure_level, routineP, error, failure)
    DO ibead=1,nbeads
       CALL lnhc_init(extended_type%nhc_part(ibead))
       CALL lnhc_init(extended_type%nhc_baro(ibead))
       CALL lnhc_init(extended_type%nhc_coef(ibead))
    ENDDO

    ! loop over number of beads
    DO ibead = 1,nbeads

      CALL force_env_get(force_env_p(ibead)%force_env, subsys=subsys, &
                         dyn_coeff_set=dyn_coeff_set, in_use=in_use, &
                         cell=cell, cell_ref = cell_ref, error=error )

! Define a reference cell for compatibility with MD restarts

      DO isubsys = 1,SIZE(subsys)
         CALL cp_subsys_get(subsys(isubsys)%subsys, &
              local_molecules_new=local_molecules, &
              molecule_kinds_new=molecule_kinds_new, &
              molecules_new=molecules_new, &
              particles=particles, &
              error=error)

         CALL initialize_nhc_part(cell,simpar,local_molecules, &
                                  molecules_new%els, molecule_kinds_new%els, &
                                  globenv, extended_type%nhc_part(ibead))
         CALL initialize_nhc_baro(simpar,globenv, extended_type%nhc_baro(ibead))
! initialize coef thermostats and setup coefficient forces, velocities and effective masses
         IF (ASSOCIATED(dyn_coeff_set)) THEN
            CALL initialize_nhc_coef(dyn_coeff_set(isubsys)%dyn_coeff_set, &
                                     globenv, extended_type%nhc_coef(ibead))
            CALL dyn_coeff_set_initialize(dyn_coeff_set(isubsys)%dyn_coeff_set, globenv) 
            CALL initialize_nhc_forces(nhc_coef=extended_type%nhc_coef(ibead))
          END IF

      END DO
      CALL initialize_npt_type(simpar, globenv, extended_type%npt_info, cell)

! Always initialize velocities (to be overwritten if needed)
      DO isubsys = 1,SIZE(subsys)
         CALL cp_subsys_get(subsys(isubsys)%subsys, &
              particles=particles, &
              error=error)
         CALL initialize_velocities(simpar, &
              particles%els, globenv)
      END DO

      IF (simpar%ensemble /= langevin_ensemble) THEN
        SELECT CASE(simpar%read_type)
        CASE(md_pos_restart)
!FM should be adapted to the existence of many fragments
          CPPrecondition(SIZE(subsys)==1,cp_failure_level,routineP,error,failure)
          isubsys=1
          CALL cp_subsys_get(subsys(isubsys)%subsys, &
               particles=particles, &
               error=error)
          nhc_part => extended_type%nhc_part(ibead)
          nhc_baro => extended_type%nhc_baro(ibead)
          CALL read_md_restart(in_use, simpar % read_type, globenv, &
               particles%els, cell, cell_ref, &
               nhc_part, nhc_baro, extended_type%npt_info)
          CALL initialize_velocities(simpar, particles%els, globenv)
        CASE(md_pos_vel_restart, md_all_restart)
!FM should be adapted to the existence of many fragments
          CPPrecondition(SIZE(subsys)==1, cp_failure_level, routineP, error, failure)
          isubsys=1
          CALL cp_subsys_get(subsys(isubsys)%subsys, &
               particles=particles, &
               error=error)
          nhc_part => extended_type%nhc_part(ibead)
          nhc_baro => extended_type%nhc_baro(ibead)
          CALL read_md_restart(in_use, simpar % read_type, globenv, &
                              particles%els, &
                              cell, cell_ref, &
                              nhc_part, &
                              nhc_baro, &
                              extended_type%npt_info)
        END SELECT
      END IF

      nhc_part => extended_type%nhc_part(ibead)
      nhc_baro => extended_type%nhc_baro(ibead)
      CALL initialize_nhc_forces (nhc_part=nhc_part, nhc_baro=nhc_baro)

      DO isubsys=1,SIZE(subsys)
        CALL cp_subsys_get(subsys(isubsys)%subsys, &
                           particles=particles, &
                           local_molecules_new=local_molecules, &
                           error=error)
      END DO
! end loop over number of beads
    END DO

    CALL md_env_create(md_env, globenv%para_env, error=error)

    CALL set_md_env(md_env=md_env, &
                    simpar=simpar, &
                    mdio=mdio, &
                    nhc_part=extended_type%nhc_part, &
                    nhc_baro=extended_type%nhc_baro, &
                    npt=extended_type%npt_info, &
                    cell=cell, &
                    force_env_p=force_env_p)

    IF (ASSOCIATED(dyn_coeff_set)) THEN
       CALL set_md_env(md_env, nhc_coef=extended_type%nhc_coef)
    END IF

    IF (simpar%ensemble == langevin_ensemble) THEN
      CALL create_wiener_process(md_env,error)
      CALL read_md_restart_file(md_env,error)
    END IF

    time_iter_start=m_walltime()

    CALL get_md_env(md_env, force_env_p=force_env_p, itimes=itimes, &
                    used_time = used_time)

    DO ibead = 1,nbeads
      IF (start_qs .AND. (globenv%prog_name_id == do_qs)) THEN
        use_harris_tmp = force_env_p(1)%force_env%qs_env%use_harris
        max_scf_tmp = force_env_p(1)%force_env%qs_env%scf_control%max_scf
        force_env_p(1)%force_env%qs_env%use_harris = .FALSE.
        force_env_p(1)%force_env%qs_env%scf_control%max_scf = 100
        CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                                         calc_force=.FALSE., error=error)
        force_env_p(1)%force_env%qs_env%use_harris = use_harris_tmp
        force_env_p(1)%force_env%qs_env%scf_control%max_scf = max_scf_tmp
      ELSE
        CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                                         calc_force=.FALSE., error=error)
      END IF

      ! CALL force_env_calc_energy_force(force_env_p(ibead)%force_env, &
      !                                  calc_force=.TRUE.)
    END DO
    !IF (globenv%prog_name_id == do_kg) THEN
    !  Epot = force_env_p(1)%force_env%kg_env%energy%total
    !  CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot, ekin=Ekin)
    !  IF (.NOT. tot_energy) THEN
    !    prev_ldmc_energy = Epot
    !  ELSE
    !    prev_ldmc_energy = Econs
    !  END IF
    !ELSE
    !  Epot = force_env_p(1)%force_env%qs_env%energy%total
    !  CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot, ekin=Ekin)
    !  IF (.NOT. tot_energy) THEN
    !    prev_ldmc_energy = Epot
    !  ELSE
    !    prev_ldmc_energy = Econs
    !  END IF
    !END IF
    CALL md_energy(md_env, globenv) ! inits the print env at itimes == 0

    CALL cp_add_iter_level(logger%iter_info, "HMC", n_rlevel_new=iter_level, &
                           error=error)

!    itimes = itimes + 1
!    CALL cp_iterate(logger%iter_info, last=(istep==simpar%nsteps), &
!                    iter_nr=itimes,error=error)
!    CALL velocity_verlet(md_env)
!    DO ibead=1,nbeads
!      CALL temperature_control(simpar, force_env_p(ibead)%force_env, error)
!    END DO
!    time_iter_stop = m_walltime()
!    used_time = time_iter_stop - time_iter_start
!    !time_iter_start = time_iter_stop
!    CALL md_energy(md_env, globenv)
!    !IF (MODULO(itimes,mdio%idump) == 0) THEN
!      IF (simpar%ensemble == langevin_ensemble) THEN
!        CALL write_md_restart_file(md_env, error=error)
!      ELSE
!        CALL dump_variables(md_env, mdio%dump_file_name, error)
!      END IF
!    !END IF

    ! Save the wavefunction history to be able to restore it in case the 
    ! first MD approach is rejected
    IF (restore_wf_past) THEN
      CALL wfi_duplicate_history(input_wf_history=force_env_p(1)%force_env%qs_env%wf_history, &
                                 output_wf_history=tmp_wf_history, &
                                 qs_env=force_env_p(1)%force_env%qs_env, error=error)
    END IF

    DO isubsys = 1,SIZE(subsys)
      !CALL cp_subsys_get(subsys=subsys(isubsys)%subsys, &
      !                   particles=particles, &
      !                   error=error)

      ALLOCATE(tmp_particle_set(particles%n_els), stat=istat)
      IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "simpar", particles%n_els*dp_size*6)

      DO iparticles = 1, particles%n_els
        tmp_particle_set(iparticles)%r(1:3) = particles%els(iparticles)%r(1:3)
        tmp_particle_set(iparticles)%v(1:3) = particles%els(iparticles)%v(1:3)
      END DO
    END DO

    istep = 1
    DO WHILE (istep < simpar%nsteps)
    !  IF (.NOT. accepted) THEN
    !    IF (simpar%ensemble == langevin_ensemble) THEN
    !      CALL read_md_restart_file(md_env,error)
    !    ELSE
    !      CPPrecondition(SIZE(subsys)==1,cp_failure_level,routineP,error,failure)
    !      isubsys=1
    !      CALL cp_subsys_get(subsys(isubsys)%subsys, &
    !           particles=particles, error=error)
    !      nhc_part => extended_type%nhc_part(ibead)
    !      nhc_baro => extended_type%nhc_baro(ibead)
    !      CALL read_md_restart(in_use, simpar%read_type, globenv, &
    !                          particles%els, cell, cell_ref, &
    !                          nhc_part, nhc_baro, extended_type%npt_info)
    !    END IF
    !  END IF

      ! previous energy
      IF (accepted) THEN
        IF (qs_ener .AND. (globenv%prog_name_id == do_qs)) THEN
          use_harris_tmp = force_env_p(1)%force_env%qs_env%use_harris
          max_scf_tmp = force_env_p(1)%force_env%qs_env%scf_control%max_scf
          force_env_p(1)%force_env%qs_env%use_harris = .FALSE.
          force_env_p(1)%force_env%qs_env%scf_control%max_scf = 100
          CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                                           calc_force=.FALSE., error=error)
          force_env_p(1)%force_env%qs_env%use_harris = use_harris_tmp
          force_env_p(1)%force_env%qs_env%scf_control%max_scf = max_scf_tmp
        ELSE
          CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                                           calc_force=.FALSE., error=error)
        END IF

        IF (globenv%prog_name_id == do_kg) THEN
          Epot = force_env_p(1)%force_env%kg_env%energy%total
            CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot, ekin=Ekin)
          IF (.NOT. tot_energy) THEN
            prev_ldmc_energy = Epot
          ELSE
            prev_ldmc_energy = Econs
          END IF
        ELSE
          Epot = force_env_p(1)%force_env%qs_env%energy%total
          CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot, ekin=Ekin)
          IF (.NOT. tot_energy) THEN
            prev_ldmc_energy = Epot
          ELSE
            prev_ldmc_energy = Econs
          END IF
        END IF

        IF (globenv%ionode .AND. globenv%print_level >= 1) THEN
          iw = globenv%scr
          WRITE(iw, *) ""
          WRITE(iw, *) "Kinetic Energy =", Ekin, " & Potential Energy = ", Epot
          WRITE(iw, *) "prev_ldmc_energy =", prev_ldmc_energy
          WRITE(iw, *) ""
        END IF
      END IF

      ! MC MD-approach
      CALL harris_mc_md_approach(md_env, globenv, simpar, logger, error)

      ! Compare energies
      IF (qs_ener .AND. (globenv%prog_name_id == do_qs)) THEN
        use_harris_tmp = force_env_p(1)%force_env%qs_env%use_harris
        max_scf_tmp = force_env_p(1)%force_env%qs_env%scf_control%max_scf
        force_env_p(1)%force_env%qs_env%use_harris = .FALSE.
        force_env_p(1)%force_env%qs_env%scf_control%max_scf = 100
        CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                                         calc_force=.FALSE., error=error)
        force_env_p(1)%force_env%qs_env%use_harris = use_harris_tmp
        force_env_p(1)%force_env%qs_env%scf_control%max_scf = max_scf_tmp
      ELSE
        CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                                         calc_force=.FALSE., error=error)
      END IF

      IF (globenv%prog_name_id == do_kg) THEN
        Epot = force_env_p(1)%force_env%kg_env%energy%total
        CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot, ekin=Ekin)
        IF (.NOT. tot_energy) THEN
          new_ldmc_energy = Epot
        ELSE
          new_ldmc_energy = Econs
        END IF
      ELSE
        Epot = force_env_p(1)%force_env%qs_env%energy%total
        CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot, ekin=Ekin)
        IF (.NOT. tot_energy) THEN
          new_ldmc_energy = Epot
        ELSE
          new_ldmc_energy = Econs
        END IF
      END IF
      ! Metropolis Monte Carlo acceptance criterion
      ionode = globenv%ionode
      iw = globenv%scr
      IF (new_ldmc_energy < prev_ldmc_energy) THEN
      !IF (new_ldmc_energy > prev_ldmc_energy) THEN
        accepted = .TRUE.
        accept = accept + 1
        IF (ionode .AND. globenv%print_level >= 1) THEN
          WRITE(iw, *) ""
          WRITE(iw, *) "Kinetic Energy =", Ekin, " & Potential Energy = ", Epot
          WRITE(iw, *) ""
          WRITE(iw, *) "New Energy =", new_ldmc_energy, "<", prev_ldmc_energy, "= Previous Energy"
          WRITE(iw, *) "Energy difference = ", new_ldmc_energy - prev_ldmc_energy
          WRITE(iw, *) ""
          WRITE(iw, *) "This HMC approach is directly accepted"
          WRITE(iw, *) "Current acceptance ration =", (accept*100)/(accept+reject), "%"
          WRITE(iw, *) "accepts =", accept, "  rejects =", reject
          WRITE(iw, *) ""
        END IF
      ELSE
        ! Random Number
        random_number = next_random_number(rng_stream=rng_stream, error=error); 

        IF (.NOT. ext_temp) THEN
          ! *** Calculate the current temperature ***
          DO isubsys=1,SIZE(subsys)

            CALL cp_subsys_get(subsys(isubsys)%subsys,&
                 local_particles=local_particles,&
                 particles=particles, atomic_kinds=atomic_kinds,&
                 error=error)

            e_kinetic = 0.0_dp

            DO ikind=1,SIZE(particles%els)
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
          CALL mp_sum(e_kinetic, globenv%para_env%group) !FM loop on all particles and avoid sum?

          temperature = 2.0_dp*e_kinetic/REAL(simpar%nfree,dp)
        END IF

        !prob = EXP((-1)*(ABS(prev_ldmc_energy-new_ldmc_energy)/(boltzmann* &
        !                 (simpar%temp_beta*(simpar%temp_ext)+simpar%temp_alpha))))
        IF (.NOT. ext_temp) THEN
          prob = EXP((-1)*(ABS(prev_ldmc_energy-new_ldmc_energy)/(simpar%temp_beta* &
                               temperature+simpar%temp_alpha)))
        ELSE
          prob = EXP((-1)*(ABS(prev_ldmc_energy-new_ldmc_energy)/(simpar%temp_beta* &
                               simpar%temp_ext+simpar%temp_alpha)))
        END IF

        IF (prob > random_number) THEN
          accepted = .TRUE.
          accept = accept + 1
          IF (ionode .AND. globenv%print_level >= 1) THEN
            WRITE(iw, *) ""

            WRITE(iw, *) ""
            WRITE(iw, *) "New Energy             =", new_ldmc_energy, ">", prev_ldmc_energy, "= Previous Energy"
            WRITE(iw, *) "Energy difference = ", new_ldmc_energy - prev_ldmc_energy
            WRITE(iw, *) "Transition probability =", prob, ">", random_number, "= Random Number"
            WRITE(iw, *) ""
            WRITE(iw, *) "This HMC approach is accepted"
            WRITE(iw, *) "Current acceptance ratio =", (accept*100)/(accept+reject), "%"
            WRITE(iw, *) "accepts =", accept, "  rejects =", reject
            WRITE(iw, *) ""
          END IF 
        ELSE
          accepted = .FALSE.
          reject = reject + 1
          IF (ionode .AND. globenv%print_level >= 1) THEN
            WRITE(iw, *) ""

            WRITE(iw, *) ""
            WRITE(iw, *) "New Energy             =", new_ldmc_energy, ">", prev_ldmc_energy, "= Previous Energy"
            WRITE(iw, *) "Energy difference = ", new_ldmc_energy - prev_ldmc_energy
            WRITE(iw, *) "Transition probability =", prob, "<", random_number, "= Random Number"
            WRITE(iw, *) ""
            WRITE(iw, *) "This HMC approach is rejected"
            WRITE(iw, *) "Current acceptance ratio =", (accept*100)/(accept+reject), "%"
            WRITE(iw, *) "accepts =", accept, "  rejects =", reject
            WRITE(iw, *) ""
          END IF
        END IF
      END IF

      IF (rand_vel) THEN
        IF (vel_scaling) THEN
          CALL initialize_velocities(simpar, particles%els, globenv)
        ELSE
          iw = globenv%scr
          IF (globenv%ionode .AND. globenv%print_level >= 1) THEN
            WRITE ( iw, '( A, A, A )' ) ' **********************', &
                 ' begin of velocity initialization ', '***********************'
          END IF

          CALL get_cp2k_units(unit)

          ! initializing the velocities
          IF (globenv%ionode) THEN

            denom = 0.0_dp
            vcom = 0.0_dp
            ! computing the COM velocities
            DO i = 1,SIZE(particles%els)
              CALL get_atomic_kind(atomic_kind = particles%els(i)%atomic_kind, mass=mass)
              IF (particles%els(i)%is_fixed) THEN
                particles%els(i)%v = 0.0_dp
              ELSE
                sigma = SQRT(simpar%temp_ext/mass)
                particles%els(i)%v(1) = gasdev(globenv%idum) * sigma
                particles%els(i)%v(2) = gasdev(globenv%idum) * sigma
                particles%els(i)%v(3) = gasdev(globenv%idum) * sigma
                vcom(1) = vcom(1) + mass * particles%els(i)%v(1)
                vcom(2) = vcom(2) + mass * particles%els(i)%v(2)
                vcom(3) = vcom(3) + mass * particles%els(i)%v(3)
                denom = denom + mass
              END IF
            END DO
            vcom = vcom / denom

            akin = 0.0_dp
            vcom_new = 0.0_dp
            DO i = 1,SIZE(particles%els)
              atomic_kind => particles%els(i)%atomic_kind
              CALL get_atomic_kind (atomic_kind = atomic_kind, mass = mass)
              ! first subtract the COM
              IF (.NOT.particles%els(i)%is_fixed) THEN
                particles%els(i)%v(:) = particles%els(i)%v(:) - vcom(:)
                ! now compute the temp
                akin = akin + 0.5_dp * (mass * particles%els(i)%v(1) * particles%els(i)%v(1))
                akin = akin + 0.5_dp * (mass * particles%els(i)%v(2) * particles%els(i)%v(2))
                akin = akin + 0.5_dp * (mass * particles%els(i)%v(3) * particles%els(i)%v(3))
              END IF
              vcom_new(:) = vcom_new(:) + mass * particles%els(i)%v(:)
            END DO
          END IF
        
#if defined(__parallel)
          DO i = 1,SIZE(particles%els)
            CALL mp_bcast(particles%els(i)%v, globenv%source, globenv%group)
          END DO
#endif

          temp = 0.0_dp
          IF (globenv%ionode .AND. globenv%print_level >= 1) THEN
            IF (simpar % nfree == 0) THEN
              temp = 2.0_dp * akin / REAL(1, KIND=dp) ! Should always be zero...
            ELSE
              temp = 2.0_dp * akin / REAL(simpar%nfree, KIND=dp)
            END IF
            WRITE ( iw, '( A, T61, F18.2, A2 )' ) &
                 ' Initial Temperature ', temp / unit%kelvin, " K"
            WRITE ( iw, '( A, T61, F20.12 )' ) &
                 ' Centre of mass velocity in direction x:', vcom_new(1)
            WRITE ( iw, '( A, T61, F20.12 )' ) &
                 ' Centre of mass velocity in direction y:', vcom_new(2)
            WRITE ( iw, '( A, T61, F20.12 )' ) &
                 ' Centre of mass velocity in direction z:', vcom_new(3)
            WRITE ( iw, '( A, A, A, / )' ) ' ***********************', &
                 ' end of velocity initialization ', '************************'
          END IF
        END IF
      END IF


      IF (accepted) THEN
      !  prev_ldmc_energy = new_ldmc_energy

        !itimes = itimes + simpar%ld_steps
        istep = istep + simpar%ld_steps
        !CALL cp_iterate(logger%iter_info, last=(istep>=simpar%nsteps), &
        !                iter_nr=itimes, error=error)

        IF (restore_wf_past) THEN
          !CALL wfi_copy_history(input_wf_history=force_env_p(1)%force_env%qs_env%wf_history, &
          !                      output_wf_history=tmp_wf_history, &
          !                      qs_env=force_env_p(1)%force_env%qs_env, error=error)
          CALL wfi_release(wf_history=tmp_wf_history, error=error)
          CALL wfi_duplicate_history(input_wf_history=force_env_p(1)%force_env%qs_env%wf_history, &
                                     output_wf_history=tmp_wf_history, &
                                     qs_env=force_env_p(1)%force_env%qs_env, error=error)
        END IF

        !IF (MODULO(itimes, mdio%idump) == 0) THEN
          IF (simpar%ensemble == langevin_ensemble) THEN
            CALL write_md_restart_file(md_env, error=error)
          ELSE
            CALL dump_variables(md_env, mdio%dump_file_name, error)
          END IF
        !END IF
      ELSE
        IF (restore_wf_past) THEN
          ! Restore the state of the wavefunction history as it was before the unsuccessful MD approach
          !CALL wfi_copy_history(input_wf_history=tmp_wf_history, &
          !                      output_wf_history=force_env_p(1)%force_env%qs_env%wf_history, &
          !                      qs_env=force_env_p(1)%force_env%qs_env, error=error)
          CALL wfi_release(wf_history=force_env_p(1)%force_env%qs_env%wf_history, error=error)
          CALL wfi_duplicate_history(input_wf_history=tmp_wf_history, &
                                     output_wf_history=force_env_p(1)%force_env%qs_env%wf_history, &
                                     qs_env=force_env_p(1)%force_env%qs_env, error=error)
        END IF

        ! If not accepted, set back velocities and positions
         DO isubsys = 1,SIZE(subsys)
           DO iparticles = 1, particles%n_els
             particles%els(iparticles)%r(1:3) = tmp_particle_set(iparticles)%r(1:3)
             particles%els(iparticles)%v(1:3) = tmp_particle_set(iparticles)%v(1:3)
           END DO
         END DO

        !IF (MODULO(itimes,mdio%idump) == 0) THEN
        !  IF (simpar%ensemble == langevin_ensemble) THEN
        !    CALL write_md_restart_file(md_env, error=error)
        !  ELSE
        !    CALL dump_variables(md_env, mdio%dump_file_name, error)
        !  END IF
        !END IF

        itimes = itimes - simpar%ld_steps
        
        CALL cp_iterate(logger%iter_info, last=.FALSE., &
                        iter_nr=itimes, error=error)
      END IF
    END DO
    ! Last iteration
    time_iter_start = m_walltime()

    itimes = itimes + 1
    CALL cp_iterate(logger%iter_info, last=.TRUE., &
                    iter_nr=itimes, error=error)

    CALL velocity_verlet(md_env)
    DO ibead=1,nbeads
      CALL temperature_control(simpar,force_env_p(ibead)%force_env,error)
    END DO
    time_iter_stop = m_walltime()
    used_time = time_iter_stop - time_iter_start
    CALL md_energy(md_env, globenv)

    IF (ionode .AND. globenv%print_level >= 1) THEN
      WRITE(iw, *) "Final acceptance ratio =", (accept*100)/(accept+reject), "%"
      WRITE(iw, *) "accepts =", accept, "  rejects =", reject
    END IF

    CALL cp_rm_iter_level(logger%iter_info,n_rlevel_att=iter_level,&
         error=error)

    CALL md_env_release (md_env, error)

    ! just deallocate this stuff ?? 
    ! seem to need proper init / destroy subroutines ... to be fixed
    CALL lnhc_dealloc(extended_type%nhc_part)
    CALL lnhc_dealloc(extended_type%nhc_baro)
    CALL lnhc_dealloc(extended_type%nhc_coef)
    CALL npt_dealloc(extended_type%npt_info)
    DEALLOCATE(simpar, STAT=istat)
    DEALLOCATE(mdio, STAT=istat)
    CPPostcondition(istat==0, cp_failure_level, routineP, error, failure)

    CALL delete_rng_stream(rng_stream=rng_stream, error=error)

    CALL wfi_release(wf_history=tmp_wf_history, error=error)
    
    CALL timestop(0.0_dp, handle)

    CALL write_checkpoint_information("leaving "//routineN,globenv)

  END SUBROUTINE harris_mc_md
!***************************************************************************

!!****f* harris_metropolis_mc/harris_mc_md_approach [1.0] *
!!
!!   NAME
!!     harris_mc_md_approach
!!
!!   FUNCTION
!!     Driver routine for the Langevin Dynamics Monte Carlo simulation
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - force_env: The force environment, which is actual a wrapper to abstract
!!                  the various force evaluation methods.
!!     - globenv: The global environment which contains the initially parsed
!!                file and the parallel environment.
!!     - error: variable to control error logging, stopping,... 
!!              see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     09.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_mc_md_approach(md_env, globenv, simpar, logger, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "harris_mc_md", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, &
                                                iter, &
                                                ibead, &
                                                nbeads
    LOGICAL                                  :: failure, &
                                                should_stop
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    INTEGER, POINTER                         :: itimes
    REAL(KIND=dp), POINTER                   :: used_time
    REAL(KIND=dp)                            :: time_iter_start, &
                                                time_iter_stop

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(force_env_p, itimes, used_time)

    CALL get_md_env(md_env, force_env_p=force_env_p, itimes=itimes, &
                    used_time=used_time)

    time_iter_start=m_walltime()
    nbeads = SIZE(force_env_p)

    DO iter = 1, simpar%ld_steps
      itimes = itimes + 1
      CALL cp_iterate(logger%iter_info, last=.FALSE., &
                      iter_nr=itimes, error=error)

      CALL velocity_verlet(md_env)

      DO ibead = 1,nbeads
        CALL temperature_control(simpar, force_env_p(ibead)%force_env, error)
      END DO
      time_iter_stop = m_walltime()
      used_time = time_iter_stop - time_iter_start
      time_iter_start = time_iter_stop
      CALL md_energy(md_env, globenv, error)

      !IF (MODULO(itimes, mdio%idump) == 0) THEN
      !  IF (simpar%ensemble == langevin_ensemble) THEN
      !    CALL write_md_restart_file(md_env, error=error)
      !  ELSE
      !    CALL dump_variables(md_env, mdio%dump_file_name, error)
      !  END IF
      !END IF

      CALL external_control(should_stop,"HMC",error=error)
      IF (should_stop) EXIT

      CALL mp_bcast(should_stop, globenv%para_env%source, globenv%para_env%group) ! sync values of should_stop
      IF (should_stop) EXIT
    END DO

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE harris_mc_md_approach
!***************************************************************************

END MODULE harris_metropolis_mc

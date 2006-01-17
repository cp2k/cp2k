!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2005  CP2K developers group                            !
!-----------------------------------------------------------------------------!


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
  USE cp_control_types,                ONLY: dft_control_type
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
                                             force_env_get, &
                                             use_qs_force, &
                                             use_kg_force
  USE global_types,                    ONLY: global_environment_type
  USE init_extended_system_types,      ONLY: initialize_nhc_part, &
                                             initialize_nhc_baro, &
                                             initialize_nhc_coef, &
                                             initialize_npt
  USE init_extended_system_variables,  ONLY: init_nhc_forces
  USE input_constants,                 ONLY: langevin_ensemble, &
                                             do_qs, &
                                             do_kg, &
                                             md_pos_restart, &
                                             md_pos_vel_restart, &
                                             restart_guess,&
                                             md_all_restart, &
                                             md_init_start, &
                                             FIST_MD_method, &
                                             KG_MD_method, &
                                             QS_MD_method, &
                                             Harris_MD_method, &
                                             KG_MC_method, &
                                             QS_MC_method, &
                                             Harris_MC_method
  USE input_section_types,             ONLY: section_vals_type,section_vals_get_subs_vals
  USE kg_environment_types,            ONLY: kg_environment_type
  USE kinds,                           ONLY: dp, &
                                             dp_size
  USE machine,                         ONLY: m_walltime, &
                                             m_flush
  USE md,                              ONLY: mdio_parameters_type, &
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
  USE md_energies,                     ONLY: open_md_file
  USE md_run,                          ONLY: temperature_control
  USE message_passing,                 ONLY: mp_bcast, &
                                             mp_sum
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE output_utilities,                ONLY: print_warning
  USE parallel_rng_types,              ONLY: rng_stream_type, &
                                             create_rng_stream, &
                                             next_random_number, &
                                             delete_rng_stream
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: boltzmann
  USE qs_environment_types,            ONLY: qs_environment_type, &
                                             get_qs_env
  USE qs_mo_types,                     ONLY: read_mo_set
                                             !write_mo_set
  USE qs_wf_history_methods,           ONLY: wfi_duplicate_history, &
                                             wfi_copy_history
  USE qs_wf_history_types,             ONLY: qs_wf_history_type, &
                                             wfi_release, &
                                             wfi_use_guess_method_nr
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory, &
                                             external_control, &
                                             stop_program
  USE timings,                         ONLY: timeset, &
                                             timestop
  USE util,                            ONLY: gasdev
  USE velocity_verlet_control,         ONLY: velocity_verlet
  USE wiener_process,                  ONLY: create_wiener_process,&
                                             read_md_restart_file,&
                                             write_md_restart_file
#include "cp_common_uses.h"

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

!   ------------------------------------------------------------------------

    ALLOCATE(force_env_p(1), stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    force_env_p(1)%force_env => force_env
    !CALL harris_mc_md(force_env_p, globenv, total_energy=.TRUE., &
    !                  external_temp=.TRUE., qs_accept=.FALSE., &
    !                  random_velocities=.TRUE., velocity_scaling=.FALSE., &
    !                  qs_start=.TRUE., restore_wf_history=.TRUE., error=error)
    CALL harris_mc_md(force_env_p=force_env_p, globenv=globenv, error=error)

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
!!     - velocity_scaling: The velocities are scaled, such that the current 
!!                         temperature is equivalent with the external temperature
!!     - qs_start: Denotes if a fully converged QS energy evaluation is performed 
!!                 at the start-up or not.
!!     - restore_wf_history: If switched on the WF history is fully restored in
!!                           case a MD approach is rejected.
!!     - restore_velocitites: Causes that the momentas are recovered in case a MD
!!                            approach is rejected.
!!     - semi_hybrid: Results that the momenta are choosen at random, only in case 
!!                    a MD approach is rejected.
!!     - md_method: Indicates which method is used for the MD approach.
!!     - mc_method: Indicates according to which method the MC criterion is evaluated.
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
  SUBROUTINE harris_mc_md(force_env_p, globenv, error)

    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(global_environment_type), POINTER   :: globenv
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
                                                in_use, &
                                                force_env_in_use, &
                                                md_method, &
                                                mc_method, &
                                                md_approach_nr, &
                                                accepted_md
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
                                                restore_wf_past, &
                                                restore_velocities, &
                                                semi_hybrid, &
                                                file_exists
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
                                                Ekin_prev, &
                                                Ekin_new, &
                                                Epot_prev, &
                                                Epot_new, &
                                                Econs, &
                                                akin, &
                                                temp, &
                                                denom, &
                                                sigma, &
                                                eps_default_tmp, &
                                                eps_scf_tmp, &
                                                acceptance_rate
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
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type),POINTER          :: dft_section, input
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER, SAVE                            :: hmc, &
                                                md_approach_old, &
                                                accepted_md_old
    CHARACTER(LEN=50)                        :: hmc_file_name

!   ------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routineN, globenv)

    CALL timeset(routineN, "I", "", handle)

    ! will be used to stop if a maximum walltime is given
    md_run_start_time=m_walltime()

    failure=.FALSE.
    accepted = .TRUE.
    accept = 0
    reject = 0

    NULLIFY(logger, simpar, mdio, subsys, dyn_coeff_set, md_env, nhc_baro, &
            nhc_part, molecule_kinds_new, molecules_new, local_molecules, &
            particles, itimes, used_time, rng_stream, tmp_particle_set, &
            local_particles, atomic_kind, atomic_kinds, atom_list, unit, &
            tmp_wf_history, dft_control, qs_env, kg_env, particle_set, atomic_kind_set)
    NULLIFY(dft_section, input)
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(globenv), cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(force_env_p), cp_failure_level, routineP, error, failure)

    CALL create_rng_stream(rng_stream=rng_stream, name="MonteCarlo acceptance", &
                           extended_precision=.TRUE., error=error)

    IF (globenv%prog_name_id == do_qs) THEN
      CALL force_env_get(force_env=force_env_p(1)%force_env, in_use=force_env_in_use, qs_env=qs_env, &
                         error=error)
    ELSE IF (globenv%prog_name_id == do_kg) THEN
      CALL force_env_get(force_env=force_env_p(1)%force_env, in_use=force_env_in_use, kg_env=kg_env, &
                         error=error)
    END IF
    IF (force_env_in_use == use_qs_force) THEN
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set,input=input, error=error)
    END IF

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

    hmc_file_name = mdio%hmc_file_name    

    tot_energy = simpar%tot_energy
    ext_temp = simpar%ext_temp
    qs_ener = simpar%qs_accept
    rand_vel = simpar%rnd_velocities
    vel_scaling = simpar%vel_scaling
    start_qs = simpar%qs_start
    restore_wf_past = simpar%restore_history
    restore_velocities = simpar%restore_momenta
    semi_hybrid = simpar%semi_hybrid
    md_method = simpar%md_method
    mc_method = simpar%mc_method

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
            CALL init_nhc_forces(nhc_coef=extended_type%nhc_coef(ibead))
          END IF

      END DO
      CALL initialize_npt(simpar, globenv, extended_type%npt_info, cell)

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
          CALL read_md_restart(in_use, simpar%read_type, globenv, &
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
          CALL read_md_restart(in_use, simpar%read_type, globenv, &
                              particles%els, &
                              cell, cell_ref, &
                              nhc_part, &
                              nhc_baro, &
                              extended_type%npt_info)
        END SELECT
      END IF

      nhc_part => extended_type%nhc_part(ibead)
      nhc_baro => extended_type%nhc_baro(ibead)
      CALL init_nhc_forces (nhc_part=nhc_part, nhc_baro=nhc_baro)

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
        use_harris_tmp = qs_env%use_harris
        max_scf_tmp = qs_env%scf_control%max_scf

        eps_default_tmp = MAX(dft_control%qs_control%eps_rho_rspace, &
                              dft_control%qs_control%eps_rho_gspace)
        eps_scf_tmp = qs_env%scf_control%eps_scf

        qs_env%use_harris = .FALSE.
        qs_env%scf_control%max_scf = 100

        qs_env%scf_control%eps_scf = 1.0e-8
        dft_control%qs_control%eps_rho_rspace = 1.0e-16
        dft_control%qs_control%eps_rho_gspace = 1.0e-16
        dft_control%qs_control%eps_gvg_rspace = 1.0e-4
        dft_control%qs_control%eps_core_charge = 1.0e-18
        dft_control%qs_control%eps_pgf_orb = 1.0e-4
        dft_control%qs_control%eps_ppl = 1.0e-18
        dft_control%qs_control%eps_ppnl = 1.0e-18
        CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                                         calc_force=.FALSE., error=error)
        qs_env%use_harris = use_harris_tmp
        qs_env%scf_control%max_scf = max_scf_tmp

        qs_env%scf_control%eps_scf = eps_scf_tmp
        dft_control%qs_control%eps_rho_rspace = eps_default_tmp
        dft_control%qs_control%eps_rho_gspace = eps_default_tmp
        dft_control%qs_control%eps_gvg_rspace = SQRT(eps_default_tmp)
        dft_control%qs_control%eps_core_charge = eps_default_tmp/100.0_dp
        dft_control%qs_control%eps_pgf_orb = SQRT(eps_default_tmp)
        dft_control%qs_control%eps_ppl = eps_default_tmp/100.0_dp
        dft_control%qs_control%eps_ppnl = eps_default_tmp/100.0_dp
      ELSE
        CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                                         calc_force=.FALSE., error=error)
      END IF
    END DO

    CALL md_energy(md_env, globenv) ! inits the print env at itimes == 0

    CALL cp_add_iter_level(logger%iter_info, "HMC", n_rlevel_new=iter_level, &
                           error=error)

    IF (simpar%ensemble == langevin_ensemble) THEN
      CALL write_md_restart_file(md_env, error=error)
    ELSE
      CALL dump_variables(md_env, mdio%dump_file_name, error)
    END IF

    ! Save the wavefunction history to be able to restore it in case the 
    ! first MD approach is rejected
    IF (restore_wf_past) THEN
      CALL wfi_duplicate_history(input_wf_history=qs_env%wf_history, &
                                 output_wf_history=tmp_wf_history, &
                                 qs_env=qs_env, error=error)
    END IF
    !IF ((qs_env%scf_control%density_guess == restart_guess) .AND. &
    !    (dft_control%qs_control%wf_interpolation_method_nr == wfi_use_guess_method_nr)) THEN
    !  CALL write_mo_set(qs_env%mos, atomic_kind_set, particle_set, globenv, id_nr=1)
    !END IF

    DO isubsys = 1,SIZE(subsys)
      !CALL cp_subsys_get(subsys=subsys(isubsys)%subsys, &
      !                   particles=particles, &
      !                   error=error)

      ALLOCATE(tmp_particle_set(particles%n_els), stat=istat)
      IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "tmp_particle_set", particles%n_els*dp_size*6)

      IF (restore_velocities) THEN
        DO iparticles = 1, particles%n_els
          tmp_particle_set(iparticles)%r(1:3) = particles%els(iparticles)%r(1:3)
          tmp_particle_set(iparticles)%v(1:3) = particles%els(iparticles)%v(1:3)
        END DO
      ELSE
        DO iparticles = 1, particles%n_els
          tmp_particle_set(iparticles)%r(1:3) = particles%els(iparticles)%r(1:3)
          tmp_particle_set(iparticles)%v(1:3) = 0.0_dp
        END DO
      END IF
    END DO

    IF (globenv%ionode) THEN
      IF (simpar%read_type == md_init_start) THEN
        CALL open_md_file(file_name=hmc_file_name, unit_number=hmc, replace_file=.TRUE.)
        md_approach_old = 0
        accepted_md_old = 0
      ELSE
        CALL open_md_file(file_name=hmc_file_name, unit_number=hmc, replace_file=.FALSE.)
        INQUIRE (FILE=hmc_file_name, EXIST=file_exists)
        IF (file_exists) THEN
          BACKSPACE(UNIT=hmc, IOSTAT=istat)
          IF (istat /= 0) THEN
            CALL stop_program(routineN, moduleN, __LINE__, &
                              "BACKSPACE command failed for file <"//&
                              TRIM(hmc_file_name)//">")
          END IF

          READ (UNIT=hmc, FMT=*, IOSTAT=istat) md_approach_old, accepted_md_old
          IF (istat /= 0) THEN
            CALL print_warning(routineN, &
                               "Could not read the last line of the file <"//&
                               TRIM(hmc_file_name)//">. MD step number "//&
                               "and time are set to zero.")
            md_approach_old = 0
            accepted_md_old = 0
          END IF
        END IF
      END IF
    END IF
    CALL mp_bcast(md_approach_old, globenv%source, globenv%group)
    CALL mp_bcast(accepted_md_old, globenv%source, globenv%group)
    md_approach_nr = md_approach_old
    accepted_md = accepted_md_old

    istep = 1
    DO WHILE (istep < simpar%nsteps)
      ! previous energy
      !IF (accepted) THEN
        IF (qs_ener .AND. (globenv%prog_name_id == do_qs)) THEN
          use_harris_tmp = qs_env%use_harris
          max_scf_tmp = qs_env%scf_control%max_scf
          qs_env%use_harris = .FALSE.
          qs_env%scf_control%max_scf = 100
          CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                                           calc_force=.FALSE., error=error)
          qs_env%use_harris = use_harris_tmp
          qs_env%scf_control%max_scf = max_scf_tmp
        ELSE
          CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                                           calc_force=.FALSE., error=error)
        END IF

        IF (globenv%prog_name_id == do_kg) THEN
          Epot_prev = force_env_p(1)%force_env%kg_env%energy%total
            CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot_prev, ekin=Ekin_prev)
          IF (.NOT. tot_energy) THEN
            prev_ldmc_energy = Epot_prev
          ELSE
            prev_ldmc_energy = Econs
          END IF
        ELSE
          Epot_prev = force_env_p(1)%force_env%qs_env%energy%total
          CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot_prev, ekin=Ekin_prev)
          IF (.NOT. tot_energy) THEN
            prev_ldmc_energy = Epot_prev
          ELSE
            prev_ldmc_energy = Econs
          END IF
        END IF

        IF (globenv%ionode .AND. globenv%print_level >= 1) THEN
          iw = globenv%scr
          WRITE(iw, *) ""
          WRITE(iw, *) "Kinetic Energy =", Ekin_prev, " & Potential Energy = ", Epot_prev
          WRITE(iw, *) "prev_ldmc_energy =", prev_ldmc_energy
          WRITE(iw, *) ""
        END IF
      !END IF

      ! MC MD-approach
      md_approach_nr = md_approach_nr+1
      CALL harris_mc_md_approach(md_env, globenv, simpar, logger, error)

      ! Compare energies
      IF (qs_ener .AND. (globenv%prog_name_id == do_qs)) THEN
        use_harris_tmp = qs_env%use_harris
        max_scf_tmp = qs_env%scf_control%max_scf
        qs_env%use_harris = .FALSE.
        qs_env%scf_control%max_scf = 100
        CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                                         calc_force=.FALSE., error=error)
        qs_env%use_harris = use_harris_tmp
        qs_env%scf_control%max_scf = max_scf_tmp
      ELSE
        CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                                         calc_force=.FALSE., error=error)
      END IF

      ! new energy
      IF (globenv%prog_name_id == do_kg) THEN
        Epot_new = force_env_p(1)%force_env%kg_env%energy%total
        CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot_new, ekin=Ekin_new)
        IF (.NOT. tot_energy) THEN
          new_ldmc_energy = Epot_new
        ELSE
          new_ldmc_energy = Econs
        END IF
      ELSE
        Epot_new = force_env_p(1)%force_env%qs_env%energy%total
        CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot_new, ekin=Ekin_new)
        IF (.NOT. tot_energy) THEN
          new_ldmc_energy = Epot_new
        ELSE
          new_ldmc_energy = Econs
        END IF
      END IF

      ! Metropolis Monte Carlo acceptance criterion
      ionode = globenv%ionode
      iw = globenv%scr
      !direct accept
      IF (new_ldmc_energy < prev_ldmc_energy) THEN
      !IF (new_ldmc_energy > prev_ldmc_energy) THEN
        accepted = .TRUE.
        accepted_md = accepted_md+1
        accept = accept + 1
        acceptance_rate = 100*(DBLE(accept)/DBLE(accept+reject))
        prob = 1.0_dp
        IF (ionode .AND. globenv%print_level >= 1) THEN
          WRITE(iw, *) ""
          WRITE(iw, *) "Kinetic Energy =", Ekin_new, " & Potential Energy = ", Epot_new
          WRITE(iw, *) ""
          WRITE(iw, *) "New Energy =", new_ldmc_energy, "<", prev_ldmc_energy, "= Previous Energy"
          WRITE(iw, *) "Energy difference = ", new_ldmc_energy - prev_ldmc_energy
          WRITE(iw, *) ""
          WRITE(iw, *) "This HMC approach is directly accepted"
          WRITE(iw, *) "Current acceptance ration =", acceptance_rate, "%"
          WRITE(iw, *) "accepts =", accept, "  rejects =", reject
          WRITE(iw, *) ""
        END IF
        random_number = 0.0_dp
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
          accepted_md = accepted_md+1
          accept = accept + 1
          acceptance_rate = 100*(DBLE(accept)/DBLE(accept+reject))
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
          acceptance_rate = 100*(DBLE(accept)/DBLE(accept+reject))
          IF (ionode .AND. globenv%print_level >= 1) THEN
            WRITE(iw, *) ""

            WRITE(iw, *) ""
            WRITE(iw, *) "New Energy             =", new_ldmc_energy, ">", prev_ldmc_energy, "= Previous Energy"
            WRITE(iw, *) "Energy difference = ", new_ldmc_energy - prev_ldmc_energy
            WRITE(iw, *) "Transition probability =", prob, "<", random_number, "= Random Number"
            WRITE(iw, *) ""
            WRITE(iw, *) "This HMC approach is rejected"
            WRITE(iw, *) "Current acceptance ratio =", acceptance_rate, "%"
            WRITE(iw, *) "accepts =", accept, "  rejects =", reject
            WRITE(iw, *) ""
          END IF
        END IF
      END IF

      ! New maxwell distributed momenta
      IF (rand_vel .OR. (semi_hybrid .AND. (.NOT. accepted))) THEN
!      IF (rand_vel) THEN
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
      ELSE
        ! Restore old momenta if the velocities are not choosen at random.
      END IF ! rand_vel

      ! Write .hmc file
      IF (ionode) THEN
        WRITE(hmc,"(I6, I5, F18.9, F18.9, F18.9, F18.9, F18.9, F18.9, F18.9, F18.9, F18.9, I8, I8, F8.4)") &
              md_approach_nr, accepted_md, prev_ldmc_energy, Epot_prev, Ekin_prev, new_ldmc_energy, &
              Epot_new, Ekin_new, (new_ldmc_energy-prev_ldmc_energy), prob, random_number, accept, &
              reject, acceptance_rate
        CALL m_flush(hmc)
      END IF

      IF (accepted) THEN
      ! prev_ldmc_energy = new_ldmc_energy

        !itimes = itimes + simpar%ld_steps
        istep = istep + simpar%ld_steps
        !CALL cp_iterate(logger%iter_info, last=(istep>=simpar%nsteps), &
        !                iter_nr=itimes, error=error)

        IF (restore_wf_past) THEN
          !CALL wfi_copy_history(input_wf_history=force_env_p(1)%force_env%qs_env%wf_history, &
          !                      output_wf_history=tmp_wf_history, &
          !                      qs_env=force_env_p(1)%force_env%qs_env, error=error)
          CALL wfi_release(wf_history=tmp_wf_history, error=error)
          CALL wfi_duplicate_history(input_wf_history=qs_env%wf_history, &
                                     output_wf_history=tmp_wf_history, &
                                     qs_env=qs_env, error=error)
        END IF
        !IF ((qs_env%scf_control%density_guess == restart_guess) .AND. &
        !    (dft_control%qs_control%wf_interpolation_method_nr == wfi_use_guess_method_nr)) THEN
        !  CALL write_mo_set(qs_env%mos, atomic_kind_set, particle_set, globenv, id_nr=1)
        !END IF

        ! If accepted, save velocities and positions
        IF (restore_velocities) THEN
          DO isubsys = 1,SIZE(subsys)
            DO iparticles = 1,particles%n_els
              tmp_particle_set(iparticles)%r(1:3) = particles%els(iparticles)%r(1:3)
              tmp_particle_set(iparticles)%v(1:3) = particles%els(iparticles)%v(1:3)
            END DO
          END DO
        ELSE
          DO isubsys = 1,SIZE(subsys)
            DO iparticles = 1,particles%n_els
              tmp_particle_set(iparticles)%r(1:3) = particles%els(iparticles)%r(1:3)
            END DO
          END DO
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
          CALL wfi_release(wf_history=qs_env%wf_history, error=error)
          CALL wfi_duplicate_history(input_wf_history=tmp_wf_history, &
                                     output_wf_history=qs_env%wf_history, &
                                     qs_env=qs_env, error=error)
        END IF
        IF ((qs_env%scf_control%density_guess == restart_guess) .AND. &
            (dft_control%qs_control%wf_interpolation_method_nr == wfi_use_guess_method_nr)) THEN
          dft_section =>  section_vals_get_subs_vals(input,"DFT",error=error)
          CALL read_mo_set(qs_env%mos, atomic_kind_set, particle_set, globenv,&
               id_nr=1,dft_section=dft_section)
        END IF

        ! If not accepted, set back velocities and positions
        IF (restore_velocities) THEN
          DO isubsys = 1,SIZE(subsys)
            DO iparticles = 1,particles%n_els
              particles%els(iparticles)%r(1:3) = tmp_particle_set(iparticles)%r(1:3)
              particles%els(iparticles)%v(1:3) = tmp_particle_set(iparticles)%v(1:3)
            END DO
          END DO
        ELSE
          DO isubsys = 1,SIZE(subsys)
            DO iparticles = 1,particles%n_els
              particles%els(iparticles)%r(1:3) = tmp_particle_set(iparticles)%r(1:3)
            END DO
          END DO
        END IF

        !IF (MODULO(itimes,mdio%idump) == 0) THEN
        !  IF (simpar%ensemble == langevin_ensemble) THEN
        !    CALL read_md_restart_file(md_env,error)
        !  ELSE
        !    CALL read_md_restart(in_use, simpar%read_type, globenv, &
        !                         particles%els, cell, cell_ref, &
        !                         nhc_part, nhc_baro, extended_type%npt_info)
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

    IF (simpar%ensemble == langevin_ensemble) THEN
      CALL write_md_restart_file(md_env, error=error)
    ELSE
      CALL dump_variables(md_env, mdio%dump_file_name, error)
    END IF
    !CALL write_mo_set(qs_env%mos, atomic_kind_set, particle_set, globenv, id_nr=1)

    IF (ionode .AND. globenv%print_level >= 1) THEN
      accept = accept+1
      acceptance_rate = 100*(DBLE(accept)/DBLE(accept+reject))
      WRITE(iw, *) "Final acceptance ratio =", acceptance_rate, "%"
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
    CPPostcondition(istat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(mdio, STAT=istat)
    CPPostcondition(istat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(tmp_particle_set, STAT=istat)
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

      CALL external_control(should_stop,"HMC",globenv,error)
      IF (should_stop) EXIT

    END DO

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE harris_mc_md_approach
!***************************************************************************

END MODULE harris_metropolis_mc

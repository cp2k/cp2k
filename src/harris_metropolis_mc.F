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
  USE cell_types,                      ONLY: cell_type
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsystem_p_type,&
                                             cp_subsystem_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_initialize,&
                                             dyn_coeff_set_p_type
  USE eip_environment,                 ONLY: eip_init
  USE eip_environment_types,           ONLY: eip_env_create,&
                                             eip_env_get,&
                                             eip_env_release,&
                                             eip_env_retain,&
                                             eip_environment_type
  USE extended_system_types,           ONLY: extended_system_type,&
                                             extended_type_init,&
                                             lnhc_dealloc,&
                                             lnhc_init,&
                                             npt_dealloc
  USE fist_environment,                ONLY: fist_init
  USE fist_environment_types,          ONLY: fist_env_create,&
                                             fist_env_release,&
                                             fist_env_retain,&
                                             fist_environment_type,&
                                             get_fist_env,&
                                             set_fist_env
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_type,&
                                             use_eip_force,&
                                             use_fist_force,&
                                             use_kg_force,&
                                             use_qs_force,&
                                             use_tbmd_force
  USE global_types,                    ONLY: global_environment_type
  USE harris_energy_types,             ONLY: harris_energy_get,&
                                             harris_energy_type
  USE harris_env_types,                ONLY: harris_env_get
  USE init_extended_system_types,      ONLY: initialize_nhc_baro,&
                                             initialize_nhc_coef,&
                                             initialize_nhc_part,&
                                             initialize_npt
  USE init_extended_system_variables,  ONLY: init_nhc_forces
  USE input_constants,                 ONLY: &
       EIP_MC_method, EIP_MD_method, Fist_MD_method, Harris_MC_method, &
       Harris_MD_method, KG_DFT_MC_method, KG_DFT_MD_method, KG_MD_method, &
       QS_MC_method, QS_MD_method, SE_MC_method, SE_MD_method, &
       TBMD_MC_method, do_eip, do_fist, do_kg, do_qs, do_tbmd, &
       langevin_ensemble, restart_guess
  USE input_cp2k_restarts,             ONLY: write_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kg_environment,                  ONLY: kg_init
  USE kg_environment_types,            ONLY: kg_env_create,&
                                             kg_env_release,&
                                             kg_env_retain,&
                                             kg_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE md,                              ONLY: initialize_velocities,&
                                             read_md_section,&
                                             simulation_parameters_type
  USE md_conserved_quantities,         ONLY: compute_conserved_quantity
  USE md_energies,                     ONLY: md_energy
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_env_create,&
                                             md_env_release,&
                                             md_environment_type,&
                                             set_md_env
  USE md_run,                          ONLY: temperature_control
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qs_environment,                  ONLY: qs_init
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_create,&
                                             qs_env_release,&
                                             qs_env_retain,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: read_mo_set
  USE qs_wf_history_methods,           ONLY: wfi_duplicate_history
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                             wfi_release,&
                                             wfi_use_guess_method_nr,&
                                             wfi_use_prev_p_method_nr
  USE tbmd_energy_types,               ONLY: tbmd_energy_get,&
                                             tbmd_energy_type
  USE tbmd_environment,                ONLY: tbmd_init
  USE tbmd_environment_types,          ONLY: tbmd_env_create,&
                                             tbmd_env_get,&
                                             tbmd_env_release,&
                                             tbmd_env_retain,&
                                             tbmd_environment_type
  USE termination,                     ONLY: external_control,&
                                             stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE velocity_verlet_control,         ONLY: velocity_verlet
  USE virial_types,                    ONLY: virial_create,&
                                             virial_retain
  USE wiener_process,                  ONLY: create_wiener_process
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'harris_metropolis_mc'

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

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_mc_ld', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, mc_prog_name_id, &
                                                md_method, stat
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(cp_subsystem_type), POINTER         :: subsystem
    TYPE(eip_environment_type), POINTER      :: eip_env
    TYPE(fist_environment_type), POINTER     :: fist_env
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: hmc_section
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env

!TYPE(cp_error_type)                      :: suberror
!   ------------------------------------------------------------------------

    NULLIFY(force_env_p, hmc_section, kg_env, qs_env, tbmd_env, eip_env, &
         subsys, subsystem, fist_env)

    CALL timeset(routineN, "I", "", handle)

    hmc_section => section_vals_get_subs_vals(globenv%input_file, &
         "MOTION%MD%HMC", error=error)
    CALL section_vals_val_get(hmc_section, "MD_METHOD", &
         i_val=md_method, error=error)

    mc_prog_name_id = globenv%prog_name_id

    SELECT CASE(md_method)
    CASE(FIST_MD_method)
      globenv%prog_name_id = do_fist
    CASE(EIP_MD_method)
      globenv%prog_name_id = do_eip
    CASE(SE_MD_method)
      globenv%prog_name_id = do_qs
    CASE(KG_MD_method)
      globenv%prog_name_id = do_kg
    CASE(KG_DFT_MD_method)
      globenv%prog_name_id = do_kg
    CASE(QS_MD_method)
      globenv%prog_name_id = do_qs
    CASE(Harris_MD_method)
      globenv%prog_name_id = do_qs
    CASE DEFAULT
      PRINT *, "The default DFT/QS-method is used for the MC energy evaluation."
      globenv%prog_name_id = do_qs
    END SELECT
    
    IF (mc_prog_name_id .NE. globenv%prog_name_id) THEN
      SELECT CASE (globenv%prog_name_id)
      CASE (do_fist)
        CALL fist_env_create(fist_env, para_env = globenv%para_env)
        CALL set_fist_env(fist_env, qmmm=.FALSE., error=error)

        ! *** Read the input and the database files and perform further  ***
        ! *** initializations for the setup of the FIST environment ***

        CALL fist_init(fist_env, globenv, error)

        CALL get_fist_env(fist_env=fist_env, subsys=subsystem, error=error)
        CALL cp_subsys_release(subsys=subsystem, error=error)

        CALL force_env_get(force_env=force_env, subsys=subsys, error=error)
        CALL set_fist_env(fist_env=fist_env, subsys=subsys(1)%subsys)

        !force_env%in_use=use_fist_force
        force_env%fist_env => fist_env
        CALL fist_env_retain(fist_env, error = error)
        CALL virial_create(force_env%virial, virial_in_use = .TRUE.)

        CALL fist_env_release (fist_env, error = error)
      CASE (do_kg)
        CALL kg_env_create(kg_env, para_env=globenv%para_env, &
                           error=error)

        CALL kg_init(kg_env, globenv, error)

        CALL cp_subsys_release(subsys=kg_env%subsys, error=error)

        CALL force_env_get(force_env=force_env, subsys=subsys, error=error)
        kg_env%subsys => subsys(1)%subsys

        !force_env%in_use = use_kg_force
        force_env%kg_env => kg_env
        CALL kg_env_retain(kg_env, error=error)
        CALL virial_create(force_env%virial, virial_in_use = .FALSE.)
        kg_env%virial => force_env%virial
        CALL virial_retain(kg_env%virial)

        CALL kg_env_release(kg_env, error=error)
      CASE (do_qs)
        CALL qs_env_create(qs_env, globenv%para_env)
         
        !   *** Read the input and the database files and perform further  ***
        !   *** initializations for the setup of the QUICKSTEP environment ***
        CALL qs_init(qs_env, globenv, error=error)

        CALL cp_subsys_release(subsys=qs_env%subsys, error=error)

        CALL force_env_get(force_env=force_env, subsys=subsys, error=error)
        qs_env%subsys => subsys(1)%subsys

        !force_env%in_use = use_qs_force
        force_env%qs_env => qs_env
        CALL qs_env_retain(qs_env,error=error)
        CALL virial_create(force_env%virial, virial_in_use = .TRUE.)
        qs_env%virial => force_env%virial
        CALL virial_retain(qs_env%virial)

        CALL qs_env_release(qs_env, error=error)
      CASE (do_tbmd)
        CALL tbmd_env_create(tbmd_env, error=error)

        CALL tbmd_init(tbmd_env, globenv, error=error)

        CALL cp_subsys_release(subsys=tbmd_env%subsystem, error=error)

        CALL force_env_get(force_env=force_env, subsys=subsys, error=error)
        tbmd_env%subsystem => subsys(1)%subsys

        !force_env%in_use = use_tbmd_force
        force_env%tbmd_env => tbmd_env
        CALL tbmd_env_retain(tbmd_env, error=error)
        CALL virial_create(force_env%virial, virial_in_use = .FALSE.)
        tbmd_env%virial => force_env%virial
        CALL virial_retain(tbmd_env%virial)

        CALL tbmd_env_release(tbmd_env, error=error)
      CASE (do_eip)
        CALL eip_env_create(eip_env, error=error)

        CALL eip_init(eip_env, globenv, error=error)

        CALL cp_subsys_release(subsys=eip_env%subsystem, error=error)

        CALL force_env_get(force_env=force_env, subsys=subsys, error=error)
        eip_env%subsystem => subsys(1)%subsys

        !force_env%in_use = use_eip_force
        force_env%eip_env => eip_env
        CALL eip_env_retain(eip_env, error=error)
        CALL virial_create(force_env%virial, virial_in_use = .FALSE.)
        eip_env%virial => force_env%virial
        CALL virial_retain(eip_env%virial)

        CALL eip_env_release(eip_env, error=error)
      CASE DEFAULT
        CALL cp_unimplemented_error(fromWhere=routineP, &
              message="This program isn't integrated in HMC yet.", &
              error=error, error_level=cp_failure_level)
      END SELECT
      !IF (ASSOCIATED(force_env)) THEN
      !  NULLIFY(meta_env)
      !  metadyn_section => section_vals_get_subs_vals(globenv%input_file, &
      !                     "MOTION%METADYN", error=suberror)
      !  CALL metadyn_read(meta_env, globenv, metadyn_section, error=suberror)
      !  CALL force_env_set(force_env, meta_env=meta_env, error=suberror)
      !  CALL meta_env_release(meta_env, error=suberror)
      !END IF
    END IF

    ALLOCATE(force_env_p(1), stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    force_env_p(1)%force_env => force_env
    CALL harris_mc_md(force_env_p=force_env_p, globenv=globenv, error=error)

    CALL timestop(0.0_dp, handle)

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
!!     - soft_start: Flushes the wavefunction history after each metropolis 
!!                   criterion evaluation and iterates to convergence afterwards, 
!!                   until it's refilled again.
!!     - fill_wf_history: Iterates to convergence until the wavefunction 
!!                        extrapolation pipeline is filled.
!!     - interMMC_interpolation: Interpolates between the Metropolis energy evaluations.
!!     - always_accepted: Accept all MD approaches
!!     - drld: Density restoring Langevin Dynamics
!!     - multiple_scf: Performs a MD with varying number of SCF cycles.
!!     - hmc_max_scf: Max SCF for DRLD & MULTIPLE_SCF
!!     - lpsc: Langevin Dynamics phase space correction
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

    CHARACTER(LEN=*), PARAMETER :: routineN = 'harris_mc_md', &
      routineP = moduleN//':'//routineN

    INTEGER :: accept, atom, force_env_in_use, handle, hmc_max_scf, i, iatom, &
      ibead, ikind, in_use, iparticles, istat, istep, isubsys, iter_level, &
      iw, max_scf_aspc_tmp, max_scf_tmp, mc_method, md_method, myi, nbeads, &
      reject, wf_interpolation_method_tmp
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    INTEGER, POINTER                         :: accepted_md, itimes, &
                                                md_approach_nr
    INTEGER, SAVE                            :: hmc
    LOGICAL :: accepted, always_accepted, drld, ext_temp, failure, &
      fill_wf_history, interMMC_extrapolation, ionode, lpsc, multiple_scf, &
      qs_ener, rand_vel, restore_velocities, restore_wf_past, semi_hybrid, &
      soft_start, start_qs, tot_energy, use_harris_tmp, vel_scaling
    REAL(KIND=dp) :: acceptance_rate, akin, denom, e_kinetic, Econs, &
      Ekin_new, Ekin_prev, Epot_new, Epot_prev, eps_default_tmp, eps_scf_tmp, &
      mass, md_run_start_time, new_ldmc_energy, prev_ldmc_energy, prob, &
      random_number, sigma, temp, temperature, time_iter_start, &
      time_iter_stop, v2
    REAL(KIND=dp), DIMENSION(3)              :: vcom, vcom_new
    REAL(KIND=dp), POINTER                   :: t, used_time
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), POINTER                  :: dyn_coeff_set
    TYPE(eip_environment_type), POINTER      :: eip_env
    TYPE(extended_system_type)               :: extended_type
    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_wf_history_type), POINTER        :: tmp_wf_history
    TYPE(reduced_particle_type), &
      DIMENSION(:), POINTER                  :: tmp_particle_set
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(section_vals_type), POINTER         :: dft_section, &
                                                force_env_section, input, &
                                                motion_section, &
                                                subsys_section, work_section
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(tbmd_energy_type), POINTER          :: tbmd_energy
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(unit_convert_type), POINTER         :: unit

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)
    NULLIFY(force_env_section, motion_section,  subsys_section, md_approach_nr,&
            accepted_md)
    force_env_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL",error=error)
    motion_section => section_vals_get_subs_vals(globenv%input_file,"MOTION",error=error)
    subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)

    ! will be used to stop if a maximum walltime is given
    md_run_start_time=m_walltime()

    failure=.FALSE.
    accepted = .TRUE.
    accept = 0
    reject = 0

    NULLIFY(logger, simpar, subsys, dyn_coeff_set, md_env, molecule_kinds_new, &
            molecules_new, local_molecules, particles, itimes, used_time, rng_stream, &
            tmp_particle_set, t, local_particles, atomic_kind, atomic_kinds, &
            atom_list, unit, tmp_wf_history, dft_control, qs_env, particle_set, &
            atomic_kind_set, kg_env, tbmd_env, eip_env, tbmd_energy, harris_energy)
    NULLIFY(dft_section, input, motion_section)
    logger => cp_error_get_logger(error)

    CALL cp_add_iter_level(logger%iter_info, "HMC", n_rlevel_new=iter_level, &
         error=error)

    CALL cp_iterate(logger%iter_info, last=.FALSE., &
         iter_nr=1, error=error)

    CPPrecondition(ASSOCIATED(globenv), cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(force_env_p), cp_failure_level, routineP, error, failure)

    CALL create_rng_stream(rng_stream=rng_stream,&
                           name="MonteCarlo acceptance",&
                           distribution_type=UNIFORM,&
                           extended_precision=.TRUE.,&
                           error=error)

    nbeads = SIZE(force_env_p)
    DO i = 1, nbeads
       CPPrecondition(ASSOCIATED(force_env_p(i)%force_env), cp_failure_level, routineP, error, failure)
    END DO

    ALLOCATE (simpar,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
         "simpar", 0)

    CALL read_md_section(simpar, globenv, error)

    tot_energy = simpar%tot_energy
    ext_temp = simpar%ext_temp
    qs_ener = simpar%qs_accept
    rand_vel = simpar%rnd_velocities
    vel_scaling = simpar%vel_scaling
    start_qs = simpar%qs_start
    restore_wf_past = simpar%restore_history
    restore_velocities = simpar%restore_momenta
    semi_hybrid = simpar%semi_hybrid
    soft_start = simpar%soft_start
    fill_wf_history = simpar%fill_wf_history
    interMMC_extrapolation = simpar%MMC_extrapolation
    always_accepted = simpar%always_accept
    drld = simpar%drld
    multiple_scf = simpar%multiple_scf
    hmc_max_scf = simpar%hmc_max_scf
    lpsc = simpar%lpsc
    md_method = simpar%md_method
    mc_method = simpar%mc_method
    motion_section => section_vals_get_subs_vals(globenv%input_file,"MOTION",error=error)

    SELECT CASE(mc_method)
    CASE(EIP_MC_method)
      globenv%prog_name_id = do_eip
    CASE(SE_MC_method)
      globenv%prog_name_id = do_qs
    CASE(TBMD_MC_method)
      globenv%prog_name_id = do_tbmd
    CASE(KG_DFT_MC_method)
      globenv%prog_name_id = do_kg
    CASE(QS_MC_method)
      globenv%prog_name_id = do_qs
    CASE(Harris_MC_method)
      globenv%prog_name_id = do_qs
    CASE DEFAULT
      globenv%prog_name_id = do_qs
    END SELECT

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

    IF (globenv%prog_name_id == do_qs) THEN
       CALL force_env_get(force_env=force_env_p(1)%force_env, in_use=force_env_in_use, qs_env=qs_env, &
            error=error)
    ELSE IF (globenv%prog_name_id == do_kg) THEN
       CALL force_env_get(force_env=force_env_p(1)%force_env, in_use=force_env_in_use, kg_env=kg_env, &
            error=error)
    ELSE IF (globenv%prog_name_id == do_tbmd) THEN
       CALL force_env_get(force_env=force_env_p(1)%force_env, in_use=force_env_in_use, tbmd_env=tbmd_env, &
            error=error)
    ELSE IF (globenv%prog_name_id == do_eip) THEN
       CALL force_env_get(force_env=force_env_p(1)%force_env, in_use=force_env_in_use, eip_env=eip_env, &
            error=error)
    END IF
    IF (force_env_in_use == use_qs_force) THEN
       CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, particle_set=particle_set, &
            atomic_kind_set=atomic_kind_set,input=input, error=error)
    END IF

    DO ibead = 1, nbeads

       CALL force_env_get( force_env_p(ibead)%force_env, subsys=subsys, &
                          dyn_coeff_set=dyn_coeff_set,in_use=in_use,&
                          cell=cell, cell_ref = cell_ref,&
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
          CALL init_nhc_forces ( nhc_part=extended_type%nhc_part(ibead),error=error)
          !
          ! Possibly restart barostat
          !
          NULLIFY(work_section)
          work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%NOSE",error=error)
          CALL initialize_nhc_baro(simpar,globenv, extended_type%nhc_baro(ibead),work_section=work_section,error=error)
          CALL init_nhc_forces ( nhc_baro=extended_type%nhc_baro(ibead),error=error)

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
       
       DO isubsys=1,SIZE(subsys)
          CALL cp_subsys_get(subsys(isubsys)%subsys,&
               particles=particles,&
               local_molecules_new=local_molecules, &
               error=error)
       END DO
       ! end loop over number of beads
    END DO

    CALL md_env_create(md_env, globenv%para_env, error=error)

    CALL set_md_env(md_env=md_env, &
         simpar=simpar, &
         nhc_part=extended_type%nhc_part, &
         nhc_baro=extended_type%nhc_baro, &
         npt=extended_type%npt_info, &
         cell=cell, &
         force_env_p=force_env_p)

    IF (ASSOCIATED(dyn_coeff_set)) THEN
       CALL set_md_env(md_env, nhc_coef=extended_type%nhc_coef)
    END IF

    IF (simpar%ensemble == langevin_ensemble) THEN
       CALL create_wiener_process(md_env,subsys_section,error)
    END IF

    time_iter_start=m_walltime()

    CALL get_md_env(md_env, force_env_p=force_env_p, itimes=itimes, t =t , &
         used_time = used_time, md_approach_nr = md_approach_nr, accepted_md = accepted_md)

    DO ibead = 1,nbeads
       IF (start_qs .AND. (globenv%prog_name_id == do_qs)) THEN
          use_harris_tmp = qs_env%use_harris
          max_scf_tmp = qs_env%scf_control%max_scf

          eps_default_tmp = MAX(dft_control%qs_control%eps_rho_rspace, &
               dft_control%qs_control%eps_rho_gspace)
          eps_scf_tmp = qs_env%scf_control%eps_scf

          qs_env%use_harris = .FALSE.
          qs_env%scf_control%max_scf = 1000

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

    CALL section_vals_val_get(motion_section,"MD%STEP_START_VAL",i_val=itimes,error=error)
    CALL section_vals_val_get(motion_section,"MD%TIME_START_VAL",r_val=t,error=error)

    CALL section_vals_val_get(motion_section,"MD%HMC%MD_APPROACH_START_VAL",&
         i_val=myi,error=error)
    md_approach_nr = myi
    CALL section_vals_val_get(motion_section,"MD%HMC%ACCEPTED_MD_START_VAL",&
         i_val=myi,error=error)    
    accepted_md = myi

    CALL md_energy(md_env, globenv) ! inits the print env at itimes == 0

    CALL write_restart(md_env=md_env,globenv=globenv,force_env_section=force_env_section,&
              error=error)

    IF ((.NOT. drld) .AND. (.NOT. always_accepted) .AND. (.NOT. multiple_scf)) THEN
      ! Save the wavefunction history to be able to restore it in case the 
      ! first MD approach is rejected
      IF (restore_wf_past) THEN
         CALL wfi_duplicate_history(input_wf_history=qs_env%wf_history, &
              output_wf_history=tmp_wf_history, &
              qs_env=qs_env, error=error)
      END IF
    END IF

    DO isubsys = 1,SIZE(subsys)

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

    ! Inter MMC extrapolation
    IF (interMMC_extrapolation) THEN
      IF ((mc_method == QS_MC_method) .OR. (mc_method == Harris_MC_method)) THEN
        wf_interpolation_method_tmp = dft_control%qs_control%wf_interpolation_method_nr
      END IF
    END IF

    istep = 1
    DO WHILE (istep < simpar%nsteps)
       IF ((.NOT. drld) .AND. (.NOT. multiple_scf) .AND. (.NOT. multiple_scf)) THEN
         ! previous energy
         !IF (accepted) THEN
         IF (qs_ener .AND. (globenv%prog_name_id == do_qs)) THEN
            use_harris_tmp = qs_env%use_harris
            max_scf_tmp = qs_env%scf_control%max_scf
            max_scf_aspc_tmp = qs_env%scf_control%max_scf_aspc
            qs_env%use_harris = .FALSE.
            qs_env%scf_control%max_scf = 1000
            qs_env%scf_control%max_scf_aspc = 1000
            CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                 calc_force=.FALSE., error=error)
            qs_env%use_harris = use_harris_tmp
            qs_env%scf_control%max_scf = max_scf_tmp
            qs_env%scf_control%max_scf_aspc = max_scf_aspc_tmp
         ELSE
            CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                 calc_force=.FALSE., error=error)
         END IF
       END IF

       ! Inter MMC extrapolation
       IF (interMMC_extrapolation) THEN
         IF ((mc_method == QS_MC_method) .OR. (mc_method == Harris_MC_method)) THEN
           dft_control%qs_control%wf_interpolation_method_nr = wf_interpolation_method_tmp
         END IF
       END IF

       IF ((.NOT. drld) .AND. (.NOT. always_accepted) .AND. (.NOT. multiple_scf)) THEN
         IF (globenv%prog_name_id == do_kg) THEN
            Epot_prev = force_env_p(1)%force_env%kg_env%energy%total
            CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot_prev, ekin=Ekin_prev)
            IF (.NOT. tot_energy) THEN
               prev_ldmc_energy = Epot_prev
            ELSE
               prev_ldmc_energy = Econs
            END IF
         ELSE IF (globenv%prog_name_id == do_qs) THEN
            IF (mc_method == Harris_MC_method) THEN
              CALL harris_env_get(harris_env=qs_env%harris_env, harris_energy=harris_energy, error=error)
              CALL harris_energy_get(harris_energy=harris_energy, Eharris=Econs, error=error)
              Epot_prev = harris_energy%sum_of_eigenvalues - 2*harris_energy%Ehartree_elec &
                        - harris_energy%integral_vxc
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
         ELSE IF (globenv%prog_name_id == do_tbmd) THEN
            CALL tbmd_env_get(tbmd_env=tbmd_env, tbmd_energy=tbmd_energy, error=error)
            CALL tbmd_energy_get(tbmd_energy=tbmd_energy, E_tb=Econs, E_kinetic=Ekin_prev, error=error)
            Epot_prev = Econs - Ekin_prev
            IF (.NOT. tot_energy) THEN
               prev_ldmc_energy = Epot_prev
            ELSE
               prev_ldmc_energy = Econs
            END IF
         ELSE IF (globenv%prog_name_id == do_eip) THEN
            CALL eip_env_get(eip_env=eip_env, eip_energy=Econs, eip_potential_energy=Epot_prev, &
                             error=error)
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
       END IF

       SELECT CASE(md_method)
       CASE(FIST_MD_method)
         globenv%prog_name_id = do_fist
         force_env_p(1)%force_env%in_use=use_fist_force
       CASE(KG_MD_method)
         globenv%prog_name_id = do_kg
         force_env_p(1)%force_env%in_use=use_kg_force
       CASE(KG_DFT_MD_method)
         globenv%prog_name_id = do_kg
         force_env_p(1)%force_env%in_use=use_kg_force
       CASE(QS_MD_method)
         globenv%prog_name_id = do_qs
         force_env_p(1)%force_env%in_use=use_qs_force
       CASE(Harris_MD_method)
         globenv%prog_name_id = do_qs
         force_env_p(1)%force_env%in_use=use_qs_force
       CASE(EIP_MD_method)
         globenv%prog_name_id = do_eip
         force_env_p(1)%force_env%in_use=use_eip_force
       CASE(SE_MD_method)
         globenv%prog_name_id = do_qs
         force_env_p(1)%force_env%in_use=use_qs_force
       CASE DEFAULT
         globenv%prog_name_id = do_qs
         force_env_p(1)%force_env%in_use=use_qs_force
       END SELECT

       ! MC MD-approach
       md_approach_nr = md_approach_nr+1
       CALL harris_mc_md_approach(md_env, globenv, simpar, logger, error)

       SELECT CASE(mc_method)
       CASE(EIP_MC_method)
         globenv%prog_name_id = do_eip
         force_env_p(1)%force_env%in_use=use_eip_force
       CASE(SE_MC_method)
         globenv%prog_name_id = do_qs
         force_env_p(1)%force_env%in_use=use_qs_force
       CASE(TBMD_MC_method)
         globenv%prog_name_id = do_tbmd
         force_env_p(1)%force_env%in_use=use_tbmd_force
       CASE(KG_DFT_MC_method)
         globenv%prog_name_id = do_kg
         force_env_p(1)%force_env%in_use=use_kg_force
       CASE(QS_MC_method)
         globenv%prog_name_id = do_qs
         force_env_p(1)%force_env%in_use=use_qs_force
       CASE(Harris_MC_method)
         globenv%prog_name_id = do_qs
         force_env_p(1)%force_env%in_use=use_qs_force
       CASE DEFAULT
         globenv%prog_name_id = do_qs
         force_env_p(1)%force_env%in_use=use_qs_force
       END SELECT

       ! Inter MMC extrapolation
       IF (interMMC_extrapolation) THEN
         IF ((mc_method == QS_MC_method) .OR. (mc_method == Harris_MC_method)) THEN
           wf_interpolation_method_tmp = dft_control%qs_control%wf_interpolation_method_nr
           dft_control%qs_control%wf_interpolation_method_nr = wfi_use_prev_p_method_nr
         END IF
       END IF

       IF (.NOT. multiple_scf) THEN
         ! Compare energies
         IF (qs_ener .AND. (globenv%prog_name_id == do_qs)) THEN
            use_harris_tmp = qs_env%use_harris
            max_scf_tmp = qs_env%scf_control%max_scf
            max_scf_aspc_tmp = qs_env%scf_control%max_scf_aspc
            qs_env%use_harris = .FALSE.
            qs_env%scf_control%max_scf = 1000
            qs_env%scf_control%max_scf_aspc = 1000
            CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                 calc_force=.FALSE., error=error)
            qs_env%use_harris = use_harris_tmp
            qs_env%scf_control%max_scf = max_scf_tmp
            qs_env%scf_control%max_scf_aspc = max_scf_aspc_tmp
         ELSE
            !IF (.NOT. multiple_scf) THEN
              IF ((drld) .AND. (globenv%prog_name_id == do_qs)) THEN
                max_scf_tmp = qs_env%scf_control%max_scf
                max_scf_aspc_tmp = qs_env%scf_control%max_scf_aspc
                qs_env%scf_control%max_scf = hmc_max_scf
                qs_env%scf_control%max_scf_aspc = hmc_max_scf
              END IF
              CALL force_env_calc_energy_force(force_env=force_env_p(1)%force_env, &
                   calc_force=.FALSE., error=error)
              IF ((drld) .AND. (globenv%prog_name_id == do_qs)) THEN
                qs_env%scf_control%max_scf = max_scf_tmp
                qs_env%scf_control%max_scf_aspc = max_scf_aspc_tmp
              END IF
            !END IF
         END IF
       END IF

       IF ((.NOT. drld) .AND. (.NOT. always_accepted) .AND. (.NOT. multiple_scf)) THEN
         ! new energy
         IF (globenv%prog_name_id == do_kg) THEN
            Epot_new = force_env_p(1)%force_env%kg_env%energy%total
            CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot_new, ekin=Ekin_new)
            IF (.NOT. tot_energy) THEN
               new_ldmc_energy = Epot_new
            ELSE
               new_ldmc_energy = Econs
            END IF
         ELSE IF (globenv%prog_name_id == do_qs) THEN
            IF (mc_method == Harris_MC_method) THEN
              CALL harris_env_get(harris_env=qs_env%harris_env, harris_energy=harris_energy, error=error)
              CALL harris_energy_get(harris_energy=harris_energy, Eharris=Econs, error=error)
              Epot_new = harris_energy%sum_of_eigenvalues - 2*harris_energy%Ehartree_elec &
                       - harris_energy%integral_vxc
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
         ELSE IF (globenv%prog_name_id == do_tbmd) THEN
            CALL tbmd_env_get(tbmd_env=tbmd_env, tbmd_energy=tbmd_energy, error=error)
            CALL tbmd_energy_get(tbmd_energy=tbmd_energy, E_tb=Econs, E_kinetic=Ekin_new, error=error)
            Epot_new = Econs - Ekin_new
            IF (.NOT. tot_energy) THEN
               new_ldmc_energy = Epot_new
            ELSE
               new_ldmc_energy = Econs
            END IF
         ELSE IF (globenv%prog_name_id == do_eip) THEN
            CALL eip_env_get(eip_env=eip_env, eip_energy=Econs, eip_potential_energy=Epot_new, &
                             error=error)
            IF (.NOT. tot_energy) THEN
               new_ldmc_energy = Epot_new
            ELSE
               new_ldmc_energy = Econs
            END IF
         END IF
       END IF

       IF (always_accepted .OR. drld .OR. multiple_scf) THEN
         accepted = .TRUE.
       ELSE
         ! Metropolis Monte Carlo acceptance criterion
         ionode = globenv%ionode
         iw = globenv%scr
         !direct accept
         IF (new_ldmc_energy < prev_ldmc_energy) THEN
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
            random_number = next_random_number(rng_stream=rng_stream,error=error)

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
       END IF

       IF ((.NOT. drld) .AND. (.NOT. multiple_scf)) THEN
         ! New maxwell distributed momenta
         IF (rand_vel .OR. (semi_hybrid .AND. (.NOT. accepted))) THEN         
            IF (vel_scaling) THEN
               NULLIFY(work_section)
               CALL initialize_velocities ( simpar, particles%els, globenv, work_section=work_section,error=error)
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
                        sigma = simpar%temp_ext/mass
                        particles%els(i)%v(1) = next_random_number(globenv%gaussian_rng_stream,&
                                                                   variance=sigma)
                        particles%els(i)%v(2) = next_random_number(globenv%gaussian_rng_stream,&
                                                                   variance=sigma)
                        particles%els(i)%v(3) = next_random_number(globenv%gaussian_rng_stream,&
                                                                   variance=sigma)
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
       END IF

       ! Write .hmc file
       hmc = cp_print_key_unit_nr(logger,motion_section,"MD%PRINT%HMC",&
            extension=".hmc", error=error)
       IF (hmc>0) THEN
          WRITE(hmc,"(I6, I5, F18.9, F18.9, F18.9, F18.9, F18.9, F18.9, F18.9, F18.9, F18.9, I8, I8, F8.4)") &
               md_approach_nr, accepted_md, prev_ldmc_energy, Epot_prev, Ekin_prev, new_ldmc_energy, &
               Epot_new, Ekin_new, (new_ldmc_energy-prev_ldmc_energy), prob, random_number, accept, &
               reject, acceptance_rate
          CALL m_flush(hmc)
       END IF
       CALL cp_print_key_finished_output(hmc,logger,motion_section,&
            "MD%PRINT%HMC", error=error) 

       IF (accepted) THEN
          istep = istep + simpar%ld_steps

          IF (.NOT. multiple_scf) THEN
            IF ((.NOT. drld) .AND. (.NOT. always_accepted)) THEN
              IF (restore_wf_past) THEN
                 CALL wfi_release(wf_history=tmp_wf_history, error=error)
                 CALL wfi_duplicate_history(input_wf_history=qs_env%wf_history, &
                      output_wf_history=tmp_wf_history, &
                      qs_env=qs_env, error=error)
              END IF
            END IF
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
          END IF
 
          CALL write_restart(md_env=md_env,globenv=globenv,force_env_section=force_env_section,&
              error=error)

       ELSE
          IF (restore_wf_past) THEN
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

    CALL velocity_verlet(md_env,error=error)
    DO ibead=1,nbeads
       CALL temperature_control(simpar,force_env_p(ibead)%force_env,error)
    END DO
    time_iter_stop = m_walltime()
    used_time = time_iter_stop - time_iter_start
    CALL md_energy(md_env, globenv)

    CALL write_restart(md_env=md_env,globenv=globenv,force_env_section=force_env_section,&
              error=error)

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
    DEALLOCATE(tmp_particle_set, STAT=istat)
    CPPostcondition(istat==0, cp_failure_level, routineP, error, failure)

    CALL delete_rng_stream(rng_stream=rng_stream,error=error)

    CALL wfi_release(wf_history=tmp_wf_history, error=error)

    CALL timestop(0.0_dp, handle)

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

    CHARACTER(LEN=*), PARAMETER :: routineN = 'harris_mc_md_approach', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ibead, iter, iw, &
                                                max_scf_aspc_tmp, &
                                                max_scf_tmp, nbeads
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: ionode, should_stop, &
                                                tot_energy
    REAL(KIND=dp) :: Econs, Ekin_new, Ekin_prev, Epot_new, Epot_prev, &
      new_lpsc_energy, prev_lpsc_energy, time_iter_start, time_iter_stop
    REAL(KIND=dp), POINTER                   :: used_time
    TYPE(eip_environment_type), POINTER      :: eip_env
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(tbmd_energy_type), POINTER          :: tbmd_energy
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(force_env_p, itimes, used_time, wf_history, qs_env, tbmd_env, &
            eip_env, tbmd_energy, harris_energy)

    CALL get_md_env(md_env, force_env_p=force_env_p, itimes=itimes, &
                    used_time=used_time)

    time_iter_start=m_walltime()
    nbeads = SIZE(force_env_p)

    qs_env => force_env_p(1)%force_env%qs_env
    wf_history => qs_env%wf_history

    max_scf_tmp = qs_env%scf_control%max_scf
    max_scf_aspc_tmp = qs_env%scf_control%max_scf_aspc

    tot_energy = simpar%tot_energy

    DO iter = 1, simpar%ld_steps
      IF ((globenv%prog_name_id == do_qs) .OR. (globenv%prog_name_id == do_kg)) THEN
        IF (simpar%soft_start) THEN
          IF ((MIN(wf_history%memory_depth,wf_history%snapshot_count) &
              >= iter)) THEN
            qs_env%scf_control%max_scf = 1000
            qs_env%scf_control%max_scf_aspc = 1000
          ELSE
            qs_env%scf_control%max_scf = max_scf_tmp
            qs_env%scf_control%max_scf_aspc = max_scf_aspc_tmp
          END IF
        ELSE
          IF (simpar%fill_wf_history) THEN
            IF ((MIN(wf_history%memory_depth,wf_history%snapshot_count) &
                > itimes)) THEN
              !max_scf_tmp = qs_env%scf_control%max_scf
              !max_scf_aspc_tmp = qs_env%scf_control%max_scf_aspc
              qs_env%scf_control%max_scf = 1000
              qs_env%scf_control%max_scf_aspc = 1000
            ELSE
              qs_env%scf_control%max_scf = max_scf_tmp
              qs_env%scf_control%max_scf_aspc = max_scf_aspc_tmp
            END IF
          END IF
        END IF
        IF (simpar%multiple_scf) THEN
          IF (simpar%ld_steps == iter) THEN
            qs_env%scf_control%max_scf = simpar%hmc_max_scf
            qs_env%scf_control%max_scf_aspc = simpar%hmc_max_scf
          END IF
        END IF
      END IF
      itimes = itimes + 1
      CALL cp_iterate(logger%iter_info, last=.FALSE., &
                      iter_nr=itimes, error=error)

      CALL velocity_verlet(md_env,error=error)

      IF (simpar%multiple_scf) THEN
        IF (simpar%ld_steps == iter) THEN
          qs_env%scf_control%max_scf = max_scf_tmp
          qs_env%scf_control%max_scf_aspc = max_scf_aspc_tmp 
        END IF
      END IF

      IF (simpar%lpsc) THEN
        IF (iter==1) THEN
          IF (globenv%prog_name_id == do_kg) THEN
            Epot_prev = force_env_p(1)%force_env%kg_env%energy%total
            CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot_prev, ekin=Ekin_prev)
            IF (.NOT. tot_energy) THEN
              prev_lpsc_energy = Epot_prev
            ELSE
              prev_lpsc_energy = Econs
            END IF
          ELSE IF (globenv%prog_name_id == do_qs) THEN
            IF (simpar%mc_method == Harris_MC_method) THEN
              CALL harris_env_get(harris_env=qs_env%harris_env, harris_energy=harris_energy, error=error)
              CALL harris_energy_get(harris_energy=harris_energy, Eharris=Econs, error=error)
              Epot_prev = harris_energy%sum_of_eigenvalues - 2*harris_energy%Ehartree_elec &
                        - harris_energy%integral_vxc
              IF (.NOT. tot_energy) THEN
                prev_lpsc_energy = Epot_prev
              ELSE
                prev_lpsc_energy = Econs
              END IF
            ELSE
              Epot_prev = force_env_p(1)%force_env%qs_env%energy%total
              CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot_prev, ekin=Ekin_prev)
              IF (.NOT. tot_energy) THEN
                prev_lpsc_energy = Epot_prev
              ELSE
                prev_lpsc_energy = Econs
              END IF
            END IF
          ELSE IF (globenv%prog_name_id == do_tbmd) THEN
            CALL tbmd_env_get(tbmd_env=tbmd_env, tbmd_energy=tbmd_energy, error=error)
            CALL tbmd_energy_get(tbmd_energy=tbmd_energy, E_tb=Econs, E_kinetic=Ekin_prev, error=error)
            Epot_prev = Econs - Ekin_prev
            IF (.NOT. tot_energy) THEN
              prev_lpsc_energy = Epot_prev
            ELSE
              prev_lpsc_energy = Econs
            END IF
          ELSE IF (globenv%prog_name_id == do_eip) THEN
            CALL eip_env_get(eip_env=eip_env, eip_energy=Econs, eip_potential_energy=Epot_prev, &
                             error=error)
            IF (.NOT. tot_energy) THEN
              prev_lpsc_energy = Epot_prev
            ELSE
              prev_lpsc_energy = Econs
            END IF
          END IF

          IF (globenv%ionode .AND. globenv%print_level >= 1) THEN
            iw = globenv%scr
            WRITE(iw, *) ""
            WRITE(iw, *) "Kinetic MD Energy =", Ekin_prev, " & Potential MD Energy = ", Epot_prev
            WRITE(iw, *) "prev_lpsc_energy =", prev_lpsc_energy
            WRITE(iw, *) ""
          END IF
        ELSE IF (iter==simpar%ld_steps) THEN
          IF (globenv%prog_name_id == do_kg) THEN
            Epot_new = force_env_p(1)%force_env%kg_env%energy%total
            CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot_new, ekin=Ekin_new)
            IF (.NOT. tot_energy) THEN
              new_lpsc_energy = Epot_new
            ELSE
              new_lpsc_energy = Econs
            END IF
          ELSE IF (globenv%prog_name_id == do_qs) THEN
            IF (simpar%mc_method == Harris_MC_method) THEN
              CALL harris_env_get(harris_env=qs_env%harris_env, harris_energy=harris_energy, error=error)
              CALL harris_energy_get(harris_energy=harris_energy, Eharris=Econs, error=error)
              Epot_new = harris_energy%sum_of_eigenvalues - 2*harris_energy%Ehartree_elec &
                       - harris_energy%integral_vxc
              IF (.NOT. tot_energy) THEN
                 new_lpsc_energy = Epot_new
              ELSE
                 new_lpsc_energy = Econs
              END IF
            ELSE
              Epot_new = force_env_p(1)%force_env%qs_env%energy%total
              CALL compute_conserved_quantity(md_env=md_env, cons=Econs, epot=Epot_new, ekin=Ekin_new)
              IF (.NOT. tot_energy) THEN
                new_lpsc_energy = Epot_new
              ELSE
                new_lpsc_energy = Econs
              END IF
            END IF
          ELSE IF (globenv%prog_name_id == do_tbmd) THEN
            CALL tbmd_env_get(tbmd_env=tbmd_env, tbmd_energy=tbmd_energy, error=error)
            CALL tbmd_energy_get(tbmd_energy=tbmd_energy, E_tb=Econs, E_kinetic=Ekin_new, error=error)
            Epot_new = Econs - Ekin_new
            IF (.NOT. tot_energy) THEN
              new_lpsc_energy = Epot_new
            ELSE
              new_lpsc_energy = Econs
            END IF
          ELSE IF (globenv%prog_name_id == do_eip) THEN
            CALL eip_env_get(eip_env=eip_env, eip_energy=Econs, eip_potential_energy=Epot_new, &
                             error=error)
            IF (.NOT. tot_energy) THEN
              new_lpsc_energy = Epot_new
            ELSE
              new_lpsc_energy = Econs
            END IF
          END IF
        END IF
        
        ionode = globenv%ionode
        iw = globenv%scr
        IF (ionode .AND. globenv%print_level >= 1) THEN
          WRITE(iw, *) ""
          WRITE(iw, *) "New MC Energy             =", new_lpsc_energy
          WRITE(iw, *) "MC Energy difference = ", new_lpsc_energy - prev_lpsc_energy
        END IF
      END IF

      

      DO ibead = 1,nbeads
        CALL temperature_control(simpar, force_env_p(ibead)%force_env, error)
      END DO
      time_iter_stop = m_walltime()
      used_time = time_iter_stop - time_iter_start
      time_iter_start = time_iter_stop
      CALL md_energy(md_env, globenv, error)
      CALL external_control(should_stop,"HMC",globenv,error)
      IF (should_stop) EXIT

    END DO

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE harris_mc_md_approach
!***************************************************************************

END MODULE harris_metropolis_mc

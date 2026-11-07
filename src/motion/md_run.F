!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Perform a molecular dynamics (MD) run using QUICKSTEP
!> \par History
!>   - Added support for Langevin regions (2014/02/05, LT)
!> \author Matthias Krack (07.11.2002)
! **************************************************************************************************
MODULE md_run
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE averages_types,                  ONLY: average_quantities_type
   USE barostat_types,                  ONLY: barostat_type,&
                                              create_barostat_type
   USE cell_types,                      ONLY: cell_type
   USE constraint_clv,                  ONLY: shake_update_colv_low
   USE cp_external_control,             ONLY: external_control
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE force_env_methods,               ONLY: force_env_calc_energy_force
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE force_env_utils,                 ONLY: force_env_rattle
   USE free_energy_methods,             ONLY: free_energy_evaluate
   USE free_energy_types,               ONLY: fe_env_create,&
                                              free_energy_type
   USE global_types,                    ONLY: global_environment_type
   USE hfx_ace_methods,                 ONLY: hfx_ace_set_dynamic_mode
   USE input_constants,                 ONLY: &
        ehrenfest, langevin_ensemble, npe_f_ensemble, npe_i_ensemble, &
        nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, npt_i_ensemble, &
        npt_ia_ensemble, nvt_ensemble, reftraj_ensemble
   USE input_cp2k_check,                ONLY: remove_restart_info
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_remove_values,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE md_ener_types,                   ONLY: create_md_ener,&
                                              md_ener_type
   USE md_energies,                     ONLY: initialize_md_ener,&
                                              md_ener_reftraj,&
                                              md_energy,&
                                              md_write_output
   USE md_environment_types,            ONLY: get_md_env,&
                                              md_env_create,&
                                              md_env_release,&
                                              md_environment_type,&
                                              need_per_atom_wiener_process,&
                                              set_md_env
   USE md_util,                         ONLY: md_output,&
                                              update_expected_temperature
   USE md_vel_utils,                    ONLY: angvel_control,&
                                              comvel_control,&
                                              setup_velocities,&
                                              temperature_control
   USE mdctrl_methods,                  ONLY: mdctrl_callback
   USE mdctrl_types,                    ONLY: mdctrl_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE metadynamics,                    ONLY: metadyn_finalise_plumed,&
                                              metadyn_forces,&
                                              metadyn_initialise_plumed,&
                                              metadyn_integrator,&
                                              metadyn_write_colvar
   USE metadynamics_types,              ONLY: set_meta_env
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_types,                  ONLY: global_constraint_type
   USE particle_list_types,             ONLY: particle_list_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt,&
                                              femtoseconds
   USE qs_electrons_charge_mixing,      ONLY: charge_mixing_maybe_rebuild,&
                                              charge_mixing_sync_endpoint_state,&
                                              compute_mixed_workfunction,&
                                              push_lambda_to_mixed_env,&
                                              qs_electrons_charge_mixing_lambda
   USE qs_electrons_ne_electrode,       ONLY: apply_ne_core_correction
   USE qs_electrons_types,              ONLY: qs_electrons_init,&
                                              qs_electrons_type
   USE qs_electrons_workfunction,       ONLY: qs_electrons_compute_workfunction
   USE qs_environment_methods,          ONLY: qs_env_time_update
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE reftraj_types,                   ONLY: create_reftraj,&
                                              reftraj_type
   USE reftraj_util,                    ONLY: initialize_reftraj,&
                                              write_output_reftraj
   USE rt_propagation,                  ONLY: rt_prop_setup
   USE simpar_methods,                  ONLY: read_md_section
   USE simpar_types,                    ONLY: create_simpar_type,&
                                              release_simpar_type,&
                                              simpar_type
   USE thermal_region_types,            ONLY: thermal_regions_type
   USE thermal_region_utils,            ONLY: create_thermal_regions,&
                                              print_thermal_regions_langevin
   USE thermostat_methods,              ONLY: create_thermostats
   USE thermostat_types,                ONLY: thermostats_type
   USE velocity_verlet_control,         ONLY: velocity_verlet
   USE virial_methods,                  ONLY: virial_evaluate
   USE virial_types,                    ONLY: virial_type
   USE wiener_process,                  ONLY: create_wiener_process,&
                                              create_wiener_process_cv
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_run'

   PUBLIC :: qs_mol_dyn

CONTAINS

! **************************************************************************************************
!> \brief Main driver module for Molecular Dynamics
!> \param force_env ...
!> \param globenv ...
!> \param averages ...
!> \param rm_restart_info ...
!> \param hmc_e_initial ...
!> \param hmc_e_final ...
!> \param mdctrl ...
! **************************************************************************************************
   SUBROUTINE qs_mol_dyn(force_env, globenv, averages, rm_restart_info, hmc_e_initial, hmc_e_final, mdctrl)

      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(average_quantities_type), OPTIONAL, POINTER   :: averages
      LOGICAL, INTENT(IN), OPTIONAL                      :: rm_restart_info
      REAL(KIND=dp), OPTIONAL                            :: hmc_e_initial, hmc_e_final
      TYPE(mdctrl_type), OPTIONAL, POINTER               :: mdctrl

      LOGICAL                                            :: my_rm_restart_info
      TYPE(md_environment_type), POINTER                 :: md_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: md_section, motion_section

      my_rm_restart_info = .TRUE.
      IF (PRESENT(rm_restart_info)) my_rm_restart_info = rm_restart_info
      NULLIFY (md_env, para_env)

      ! Tell ACE that this is a dynamic run: Bypass C will use full HFX
      ! for the entire first MD step so wavefunction propagation delivers
      ! a near-converged C_occ to step 1, making the ACE BUILD accurate.
      CALL hfx_ace_set_dynamic_mode(.TRUE.)

      para_env => force_env%para_env
      motion_section => section_vals_get_subs_vals(force_env%root_section, "MOTION")
      md_section => section_vals_get_subs_vals(motion_section, "MD")

      ! Real call to MD driver - Low Level
      ALLOCATE (md_env)
      CALL md_env_create(md_env, md_section, para_env, force_env=force_env)
      CALL set_md_env(md_env, averages=averages)
      IF (PRESENT(hmc_e_initial) .AND. PRESENT(hmc_e_final)) THEN
         CALL qs_mol_dyn_low(md_env, md_section, motion_section, force_env, globenv, &
                             hmc_e_initial=hmc_e_initial, hmc_e_final=hmc_e_final)
      ELSE
         CALL qs_mol_dyn_low(md_env, md_section, motion_section, force_env, globenv, mdctrl=mdctrl)
      END IF
      CALL md_env_release(md_env)
      DEALLOCATE (md_env)

      ! Clean restartable sections..
      IF (my_rm_restart_info) CALL remove_restart_info(force_env%root_section)
   END SUBROUTINE qs_mol_dyn

! **************************************************************************************************
!> \brief Purpose: Driver routine for MD run using QUICKSTEP.
!> \param md_env ...
!> \param md_section ...
!> \param motion_section ...
!> \param force_env ...
!> \param globenv ...
!> \param hmc_e_initial ...
!> \param hmc_e_final ...
!> \param mdctrl ...
!> \par History
!>   - Cleaning (09.2007) Teodoro Laino [tlaino] - University of Zurich
!>   - Added lines to print out langevin regions (2014/02/04, LT)
!> \author Creation (07.11.2002,MK)
! **************************************************************************************************
   SUBROUTINE qs_mol_dyn_low(md_env, md_section, motion_section, force_env, globenv, hmc_e_initial, hmc_e_final, mdctrl)

      TYPE(md_environment_type), POINTER                 :: md_env
      TYPE(section_vals_type), POINTER                   :: md_section, motion_section
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(global_environment_type), POINTER             :: globenv
      REAL(KIND=dp), OPTIONAL                            :: hmc_e_initial, hmc_e_final
      TYPE(mdctrl_type), OPTIONAL, POINTER               :: mdctrl

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'qs_mol_dyn_low'

      CHARACTER(LEN=23)                                  :: cpcont_target_label
      CHARACTER(LEN=6)                                   :: cpcont_pos
      CHARACTER(LEN=80)                                  :: md_step_line
      CHARACTER(LEN=default_string_length)               :: my_act, my_pos
      INTEGER                                            :: charge_span, cpcont_step0, handle, i, &
                                                            istep, iw, md_stride, nelec_report, &
                                                            output_unit, run_type_id, wf_axis
      INTEGER, POINTER                                   :: itimes
      LOGICAL :: charge_mixing, check, cm_inplace, control_fermi, ehrenfest_md, electrons_enabled, &
         electrons_explicit, electrons_report, lambda_is_fixed, save_mem, should_stop, &
         write_binary_restart_file
      REAL(KIND=dp) :: cur_ef, cur_wf, dummy, grand_canonical_eV, lambda_fix, lambda_mix, &
         legendre_eV, mass_electrons, neutral_e, p_init, pot_energy, q_init, target_col, &
         target_ef, target_wf, tau_pot, temp_pot, time_iter_start, time_iter_stop, wf_max_a, &
         wf_min_a
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wf_region
      REAL(KIND=dp), POINTER                             :: constant, time, used_time
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(barostat_type), POINTER                       :: barostat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys, subsys_i
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(free_energy_type), POINTER                    :: fe_env
      TYPE(md_ener_type), POINTER                        :: md_ener
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(qs_electrons_type), POINTER                   :: electrons
      TYPE(qs_environment_type), POINTER                 :: qs_env_report
      TYPE(reftraj_type), POINTER                        :: reftraj
      TYPE(section_vals_type), POINTER :: constraint_section, electrons_section, &
         force_env_section, free_energy_section, global_section, reftraj_section, subsys_section, &
         work_section
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(thermal_regions_type), POINTER                :: thermal_regions
      TYPE(thermostats_type), POINTER                    :: thermostats
      TYPE(virial_type), POINTER                         :: virial

! Restartable state of the electronic extended variable: read from
! &ELECTRONS at start-up, written back by store_electrons_restart.
! Fermi-level control (mutually exclusive with work-function control)
! Frozen mixing weight (diagnostic control experiment)
! Report-only mode, selected by "&ELECTRONS 0": section present, parameter
! false. electrons_explicit is what separates that from an absent section,
! since both leave _SECTION_PARAMETERS_ at .FALSE.
! Column 9 of the .cpcont and its header. The column carries whichever
! target the run is actually controlling, so the label has to follow it.
! .cpcont append-vs-truncate decision
! legendre_eV is the Legendre term Phi_target*(q - N_neutral), used by the
! grand-canonical potential below. It is not added to Conserved[eV], which is
! CP2K's own md_ener%constant, so that column means the same thing as in any
! other MD run. The electronic contributions to an extended-system invariant
! can be reconstructed from the printed columns:
!   p^2/2m           = 0.5_dp*electrons%p**2/electrons%mass * evolt
!   thermostat heat  = electrons%thermostat_energy * evolt
!   Legendre term    = Grand_canonical_potential[eV] - Potential[eV]

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(globenv))
      CPASSERT(ASSOCIATED(force_env))

      NULLIFY (particles, cell, simpar, itimes, used_time, subsys, &
               md_ener, thermostats, barostat, reftraj, force_env_section, &
               reftraj_section, work_section, atomic_kinds, &
               local_particles, time, fe_env, free_energy_section, &
               constraint_section, thermal_regions, virial, subsys_i)
      logger => cp_get_default_logger()
      para_env => force_env%para_env
      output_unit = cp_logger_get_default_io_unit(logger)

      global_section => section_vals_get_subs_vals(force_env%root_section, "GLOBAL")
      free_energy_section => section_vals_get_subs_vals(motion_section, "FREE_ENERGY")
      constraint_section => section_vals_get_subs_vals(motion_section, "CONSTRAINT")
      CALL section_vals_val_get(global_section, "SAVE_MEM", l_val=save_mem)

      CALL section_vals_val_get(global_section, "RUN_TYPE", i_val=run_type_id)
      IF (run_type_id == ehrenfest) CALL set_md_env(md_env, ehrenfest_md=.TRUE.)

      CALL create_simpar_type(simpar)
      force_env_section => force_env%force_env_section
      subsys_section => section_vals_get_subs_vals(force_env_section, "SUBSYS")
      CALL cp_add_iter_level(logger, "MD")
      CALL cp_iterate(logger%iter_info, iter_nr=0)
      ! Read MD section
      CALL read_md_section(simpar, motion_section, md_section)
      ! Setup print_keys
      simpar%info_constraint = cp_print_key_unit_nr(logger, constraint_section, &
                                                    "CONSTRAINT_INFO", extension=".shakeLog", log_filename=.FALSE.)
      simpar%lagrange_multipliers = cp_print_key_unit_nr(logger, constraint_section, &
                                                         "LAGRANGE_MULTIPLIERS", extension=".LagrangeMultLog", log_filename=.FALSE.)
      simpar%dump_lm = BTEST(cp_print_key_should_output(logger%iter_info, constraint_section, &
                                                        "LAGRANGE_MULTIPLIERS"), cp_p_file)

      ! Create the structure for the md energies
      ALLOCATE (md_ener)
      CALL create_md_ener(md_ener)
      CALL set_md_env(md_env, md_ener=md_ener)
      NULLIFY (md_ener)

      ! If requested setup Thermostats
      CALL create_thermostats(thermostats, md_section, force_env, simpar, para_env, &
                              globenv, global_section)

      ! If requested setup Barostat
      CALL create_barostat_type(barostat, md_section, force_env, simpar, globenv)

      ! If requested setup different thermal regions
      CALL create_thermal_regions(thermal_regions, md_section, simpar, force_env)

      NULLIFY (electrons_section, wf_region)
      electrons_section => section_vals_get_subs_vals(md_section, "ELECTRONS")
      CALL section_vals_val_get(electrons_section, "_SECTION_PARAMETERS_", l_val=electrons_enabled)

      ! "&ELECTRONS 0" -- section present with a false parameter -- selects
      ! report-only mode: the .cpcont file is written for an ordinary MD run,
      ! with the Fermi level and work function computed and logged every step,
      ! but q is never propagated and no potentiostat, charge mixing or core
      ! correction is applied. An ABSENT section must keep meaning "do nothing
      ! at all", and both cases leave _SECTION_PARAMETERS_ at .FALSE., so the
      ! two are separated by asking whether the section was written explicitly.
      ! (CP2K's logical parser accepts 0/1 as well as F/T, so "&ELECTRONS 0"
      ! and "&ELECTRONS F" are the same thing.)
      CALL section_vals_get(electrons_section, explicit=electrons_explicit)
      electrons_report = electrons_explicit .AND. (.NOT. electrons_enabled)

      NULLIFY (electrons, qs_env_report)
      IF (electrons_enabled) THEN
         ! The constant-potential hooks (electronic kick/drift, endpoint rebuild,
         ! lambda push, potentiostat) live inside the per-ensemble integrator
         ! routines in integrator.F, and are currently implemented only in nvt()
         ! and langevin().  Under any other ensemble the atoms would propagate
         ! normally while q(t) stayed frozen at Q_INIT and the mixing weight was
         ! never updated -- the run would finish and write entirely plausible
         ! output that means nothing.  Refuse rather than fail silently.
         IF (simpar%ensemble /= nvt_ensemble .AND. &
             simpar%ensemble /= langevin_ensemble) THEN
            CPABORT("&ELECTRONS requires ENSEMBLE NVT or LANGEVIN")
         END IF

         CALL section_vals_val_get(electrons_section, "POTENTIAL_AVERAGE_AXIS", i_val=wf_axis)
         ! r_vals returns a POINTER INTO THE INPUT TREE, not a copy. The previous
         ! line here was
         !
         !    wf_region = wf_region*angstrom
         !
         ! which converted the stored value from bohr to angstrom IN PLACE,
         ! corrupting the input tree. On a fresh run that was invisible: the
         ! converted numbers were what the work-function routine wanted, and
         ! nothing read the section again. But CP2K's restart writer does read it,
         ! and converts bohr -> angstrom a second time on output, so every restart
         ! divided the averaging region by 1.8897 (bohr per angstrom) and the
         ! error COMPOUNDED with each continuation. A 12.5-15.0 angstrom vacuum
         ! window became 6.61-7.94 after one restart -- inside the slab rather
         ! than in the vacuum, making the work function, dPhi/dq and hence the
         ! whole potentiostat meaningless while still looking plausible.
         !
         ! Convert into local scalars and leave the input tree alone.
         CALL section_vals_val_get(electrons_section, "POTENTIAL_AVERAGE_REGION", r_vals=wf_region)
         wf_min_a = wf_region(1)*angstrom
         wf_max_a = wf_region(2)*angstrom
         CALL section_vals_val_get(electrons_section, "NEUTRAL_ELECTRONS", r_val=neutral_e)
         CALL section_vals_val_get(electrons_section, "Q_INIT", r_val=q_init)
         CALL section_vals_val_get(electrons_section, "TARGET_WORKFUNCTION", r_val=target_wf)
         ! TARGET_WORKFUNCTION carries unit_str="eV", so the input machinery
         ! returns it in atomic units (Hartree) whether the deck wrote a bare
         ! number, [eV] or [hartree].  The potentiostat works in eV throughout --
         ! qs_electrons_force forms (current_workfunction - target_workfunction)
         ! in eV before dividing by HARTREE_TO_EV, and the .cpcont writes both
         ! columns as eV -- so convert back here, exactly as
         ! POTENTIAL_AVERAGE_REGION above is converted with `angstrom'.
         target_wf = target_wf*evolt
         ! Fermi-level control is selected simply by setting TARGET_FERMI_LEVEL;
         ! its default is HUGE, the sentinel for "not given".
         CALL section_vals_val_get(electrons_section, "CHARGE_MIXING_INPLACE", l_val=cm_inplace)
         CALL section_vals_val_get(electrons_section, "LAMBDA", r_val=lambda_fix)
         lambda_is_fixed = (lambda_fix < 0.5_dp*HUGE(0.0_dp))
         IF (.NOT. lambda_is_fixed) lambda_fix = 0.0_dp
         IF (lambda_is_fixed) THEN
            CPWARN("LAMBDA is fixed: q(t) frozen, potentiostat inactive, no endpoint rebuilds")
         END IF
         CALL section_vals_val_get(electrons_section, "TARGET_FERMI_LEVEL", r_val=target_ef)
         control_fermi = (target_ef < 0.5_dp*HUGE(0.0_dp))
         IF (control_fermi) THEN
            ! TARGET_FERMI_LEVEL takes PRECEDENCE rather than being mutually
            ! exclusive. It cannot be an abort: TARGET_WORKFUNCTION has a
            ! non-zero default (5.0 eV), so "was it set?" is unanswerable from
            ! its value alone, and an exclusivity check would fire on every
            ! Fermi-mode run. Precedence also means a restart file written in
            ! work-function mode -- which necessarily carries an explicit
            ! TARGET_WORKFUNCTION -- can be switched to Fermi control by adding
            ! one line, without having to delete the old target by hand.
            IF (ABS(target_wf - 5.0_dp) > 1.0E-12_dp) THEN
               CPWARN("TARGET_FERMI_LEVEL is set: TARGET_WORKFUNCTION is ignored")
            END IF
         ELSE
            target_ef = 0.0_dp
         END IF
         CALL section_vals_val_get(electrons_section, "MASS", r_val=mass_electrons)
         CALL section_vals_val_get(electrons_section, "TAU_POTENTIOSTAT", r_val=tau_pot)
         CALL section_vals_val_get(electrons_section, "TEMP_POTENTIOSTAT", r_val=temp_pot)
         CALL section_vals_val_get(electrons_section, "CHARGE_MIXING", l_val=charge_mixing)
         charge_span = 1
         IF (charge_mixing) CALL section_vals_val_get(electrons_section, "CHARGE_MIXING_SPAN", i_val=charge_span)
         ! The mixing code clamps S with MAX(1, S), but a zero or negative span is
         ! an input error, not something to correct silently.
         IF (charge_span < 1) CPABORT("CHARGE_MIXING_SPAN must be a positive integer")
         ! Restartable state. On a fresh run these carry their defaults (p = 0,
         ! no accumulated potentiostat heat) so behaviour is unchanged; on a
         ! restart they carry the live values written by store_electrons_restart
         ! below, and the electronic dynamics continues rather than restarting.
         CALL section_vals_val_get(electrons_section, "P_INIT", r_val=p_init)
         CALL section_vals_val_get(electrons_section, "POTENTIOSTAT_ENERGY", r_val=pot_energy)
         CALL section_vals_val_get(electrons_section, "CURRENT_FERMI_LEVEL", r_val=cur_ef)
         CALL section_vals_val_get(electrons_section, "CURRENT_WORKFUNCTION", r_val=cur_wf)

         IF (charge_mixing) THEN
            CALL qs_electrons_init(electrons, q_init=q_init, p_init=p_init, mass=mass_electrons, &
                                   neutral_electrons=neutral_e, target_workfunction=target_wf, &
                                   tau_potentiostat=tau_pot, temp_potentiostat=temp_pot, &
                                   potential_average_axis=wf_axis, potential_average_min=wf_min_a, &
                                   potential_average_max=wf_max_a, charge_mixing_span=charge_span, &
                                   target_fermi_level=target_ef, control_fermi=control_fermi, &
                                   lambda_fixed=lambda_fix, lambda_is_fixed=lambda_is_fixed, &
                                   charge_mixing_inplace=cm_inplace)
         ELSE
            CALL qs_electrons_init(electrons, q_init=q_init, p_init=p_init, mass=mass_electrons, &
                                   neutral_electrons=neutral_e, target_workfunction=target_wf, &
                                   tau_potentiostat=tau_pot, temp_potentiostat=temp_pot, &
                                   potential_average_axis=wf_axis, potential_average_min=wf_min_a, &
                                   potential_average_max=wf_max_a, &
                                   target_fermi_level=target_ef, control_fermi=control_fermi, &
                                   lambda_fixed=lambda_fix, lambda_is_fixed=lambda_is_fixed, &
                                   charge_mixing_inplace=cm_inplace)
         END IF
         ! qs_electrons_init zeroes these; restore them after it.
         !
         ! p and thermostat_energy are load-bearing: without them a restart
         ! resumes with the electronic variable at rest and its potentiostat
         ! heat accounting reset to zero, mid-trajectory.
         !
         ! current_fermi_level_eV and current_workfunction are NOT load-bearing:
         ! the initial force evaluation below (line ~523) recomputes both before
         ! the first electronic kick, on restarts as well as fresh runs. They are
         ! restored only so the step-0 .cpcont line of a continued run reports the
         ! values the previous run ended on rather than zeros.
         electrons%p = p_init
         electrons%thermostat_energy = pot_energy
         electrons%current_fermi_level_eV = cur_ef
         electrons%current_workfunction = cur_wf

         IF (charge_mixing) THEN
            ! Normalise the endpoints BEFORE the first force evaluation, so that
            ! slot 1 holds N_lo and slot 2 holds N_hi whatever order the input's
            ! two FORCE_EVAL blocks were written in, and whether or not their
            ! CHARGE values bracket Q_INIT.  Without this the step-0 energy,
            ! Fermi level and work function are computed from endpoints that need
            ! not bracket q -- the "startup endpoint mismatch" -- and if the two
            ! blocks happen to be in the wrong order the mixing is inverted for
            ! the whole run.  Reading the counts back off the sub_force_envs
            ! first means the common case -- an input that already agrees with
            ! Q_INIT -- costs nothing at all, while a disagreeing or misordered
            ! one is corrected before it can affect a single step.
            CALL charge_mixing_sync_endpoint_state(force_env, electrons)
            CALL charge_mixing_maybe_rebuild(force_env, electrons)
            CALL qs_electrons_charge_mixing_lambda(electrons, lambda_mix)
            CALL push_lambda_to_mixed_env(force_env, lambda_mix)
            ! Report the bracket explicitly.  The whole scheme rests on slot 1
            ! being the LOW endpoint, and until now nothing in the output let
            ! you check that; a misordered input was invisible.
            iw = cp_logger_get_default_io_unit()
            IF (iw > 0) THEN
               WRITE (iw, '(/,T2,A)') REPEAT("-", 77)
               WRITE (iw, '(T2,A)') "CHARGE MIXING| endpoint bracket at start of MD"
               WRITE (iw, '(T2,A,T50,I20)') "CHARGE MIXING|   FORCE_EVAL 1, N_lo [electrons]", &
                  electrons%charge_mixing_current_n_lo
               WRITE (iw, '(T2,A,T50,I20)') "CHARGE MIXING|   FORCE_EVAL 2, N_hi [electrons]", &
                  electrons%charge_mixing_current_n_hi
               WRITE (iw, '(T2,A,T50,F20.10)') "CHARGE MIXING|   q [electrons]", electrons%q
               WRITE (iw, '(T2,A,T50,F20.10)') "CHARGE MIXING|   lambda (weight on N_hi)", lambda_mix
               WRITE (iw, '(T2,A)') REPEAT("-", 77)
            END IF
         ELSE
            CALL apply_ne_core_correction(force_env, electrons%q, electrons%neutral_electrons)
         END IF

         logger => cp_get_default_logger()
         output_unit = cp_logger_get_default_io_unit()
         IF (output_unit > 0) THEN
            ! ALWAYS append, never truncate. This deliberately does not mirror
            ! CP2K's own .ener file, which rewinds on a fresh run.
            !
            ! The previous version keyed the choice on STEP_START_VAL > 0, which
            ! is set only when a restart file supplies it. Any other route back
            ! into a directory -- resuming from copied coordinates, a restart
            ! that did not carry STEP_START_VAL, or simply relaunching after a
            ! job was killed mid-run -- left it at 0 and silently overwrote the
            ! whole preceding trajectory. The asymmetry matters: an overwritten
            ! trajectory is unrecoverable, whereas a file that accumulates across
            ! runs costs only a little disk and is trivially split on the header
            ! lines. STEP_START_VAL is still read, but now only to label the
            ! boundary in the header below.
            CALL section_vals_val_get(md_section, "STEP_START_VAL", i_val=cpcont_step0)
            cpcont_pos = "APPEND"
            CALL open_file(file_name=TRIM(logger%iter_info%project_name)//"-1.cpcont", &
                           file_status="UNKNOWN", file_form="FORMATTED", file_action="WRITE", &
                           file_position=TRIM(cpcont_pos), unit_number=electrons%cpcont_unit)
            ! Column meanings:
            !
            ! Conserved[eV] -- CP2K's own md_ener%constant: the atomic degrees of
            !   freedom plus, under NVT, the atomic thermostat. It does not include
            !   the electronic kinetic energy, the potentiostat heat or the Legendre
            !   term, so that it keeps its usual meaning in every CP2K MD run.
            !
            ! Grand_canonical_potential[eV] -- Omega = E_pot - mu*(q - N_neutral).
            !   In work-function control mu = -Phi_target, so
            !   Omega = E_pot + Phi_target*(q - N_neutral).
            ! Header on every run, fresh or continued. Since the file is now
            ! always opened APPEND, this lands at the start of each run's block
            ! and marks where that run began -- the boundary is explicit rather
            ! than inferred. It costs nothing downstream: gnuplot skips '#' lines
            ! by default, and the awk idiom used throughout this project already
            ! filters them with !/^#/.
            ! Column 9 holds whichever target this run controls, so the header
            ! names that one. Under Fermi control TARGET_WORKFUNCTION is ignored
            ! everywhere else, and labelling the column after it made a
            ! Fermi-mode file look as though a 5.0 eV work-function target were
            ! in force.
            IF (control_fermi) THEN
               cpcont_target_label = "Fermi_target[eV]"
            ELSE
               cpcont_target_label = "Workfunction_target[eV]"
            END IF
            WRITE (electrons%cpcont_unit, '(A,I0)') &
               "# --- run opened here, STEP_START_VAL = ", cpcont_step0
            WRITE (electrons%cpcont_unit, '(A)') "# Step  Time[fs]  Potential[eV]  Conserved[eV]  "// &
               "Temperature[K]  Electronic_charge  Fermi_level[eV]  Workfunction[eV]  "// &
               TRIM(cpcont_target_label)//"  Grand_canonical_potential[eV]"
            CALL m_flush(electrons%cpcont_unit)
         END IF
      END IF

      ! ------------------------------------------------------------------
      ! Report-only mode: "&ELECTRONS 0"
      ! ------------------------------------------------------------------
      IF (electrons_report) THEN
         ! Plain single-FORCE_EVAL only. A MIXED force_env has no single qs_env
         ! to read the electron count or the Hartree potential from, and routing
         ! through compute_mixed_workfunction would drag in exactly the
         ! charge-mixing machinery this mode exists to stay clear of.
         IF (ASSOCIATED(force_env%sub_force_env)) THEN
            CPABORT("&ELECTRONS 0 (report-only) requires a plain FORCE_EVAL, not METHOD MIXED")
         END IF

         ! Unlike the full mode this works under any ensemble: q is never
         ! propagated, so there is no integrator hook that has to exist.
         CALL section_vals_val_get(electrons_section, "POTENTIAL_AVERAGE_AXIS", i_val=wf_axis)
         CALL section_vals_val_get(electrons_section, "POTENTIAL_AVERAGE_REGION", r_vals=wf_region)
         wf_min_a = wf_region(1)*angstrom
         wf_max_a = wf_region(2)*angstrom

         ! Electronic_charge is the system's real electron count, constant for
         ! the whole run. Read it from the environment rather than from
         ! NEUTRAL_ELECTRONS so the column is right whatever DFT%CHARGE was set
         ! to, and so the keyword is not silently required in this mode.
         CALL force_env_get(force_env, qs_env=qs_env_report)
         CALL get_qs_env(qs_env_report, nelectron_total=nelec_report)

         ! target_workfunction = 0 does double duty: the Workfunction_target
         ! column reads 0 throughout, as intended, and the Legendre term
         ! Phi_target*(q - N) vanishes so Grand_canonical_potential reduces to
         ! E_pot. Every column therefore keeps the same definition it has in the
         ! full mode -- the two file types stay directly comparable, and the same
         ! plotting scripts work on both. mass/tau/temp are never used because
         ! report_only gates every kick, drift and thermostat call in
         ! integrator.F; they are given harmless non-zero values only so that
         ! nothing downstream can divide by zero.
         CALL qs_electrons_init(electrons, q_init=REAL(nelec_report, KIND=dp), p_init=0.0_dp, &
                                mass=1.0_dp, neutral_electrons=REAL(nelec_report, KIND=dp), &
                                target_workfunction=0.0_dp, tau_potentiostat=1.0_dp, &
                                temp_potentiostat=0.0_dp, &
                                potential_average_axis=wf_axis, potential_average_min=wf_min_a, &
                                potential_average_max=wf_max_a, report_only=.TRUE.)

         logger => cp_get_default_logger()
         output_unit = cp_logger_get_default_io_unit()
         IF (output_unit > 0) THEN
            CALL section_vals_val_get(md_section, "STEP_START_VAL", i_val=cpcont_step0)
            CALL open_file(file_name=TRIM(logger%iter_info%project_name)//"-1.cpcont", &
                           file_status="UNKNOWN", file_form="FORMATTED", file_action="WRITE", &
                           file_position="APPEND", unit_number=electrons%cpcont_unit)
            ! Report-only controls nothing, so column 9 has no target to name.
            ! It is written as 0 (target_workfunction is 0 in this mode and
            ! control_fermi is false) and the header says so plainly rather than
            ! borrowing a label that would imply a target is in force.
            cpcont_target_label = "Nothing_important"
            WRITE (electrons%cpcont_unit, '(A,I0,A)') &
               "# --- run opened here, STEP_START_VAL = ", cpcont_step0, &
               "   (report-only: no potentiostat, q fixed at the system electron count)"
            WRITE (electrons%cpcont_unit, '(A)') "# Step  Time[fs]  Potential[eV]  Conserved[eV]  "// &
               "Temperature[K]  Electronic_charge  Fermi_level[eV]  Workfunction[eV]  "// &
               TRIM(cpcont_target_label)//"  Grand_canonical_potential[eV]"
            CALL m_flush(electrons%cpcont_unit)
         END IF
      END IF

      ! If doing langevin_ensemble, then print out langevin_regions information upon request
      IF (simpar%ensemble == langevin_ensemble) THEN
         my_pos = "REWIND"
         my_act = "WRITE"
         CALL print_thermal_regions_langevin(thermal_regions, simpar, &
                                             pos=my_pos, act=my_act)
      END IF

      CALL set_md_env(md_env, thermostats=thermostats, barostat=barostat, thermal_regions=thermal_regions, &
                      electrons=electrons)

      CALL get_md_env(md_env, ehrenfest_md=ehrenfest_md)

      !If requested set up the REFTRAJ run
      IF (simpar%ensemble == reftraj_ensemble .AND. ehrenfest_md) THEN
         CPABORT("Ehrenfest MD does not support reftraj ensemble ")
      END IF
      IF (simpar%ensemble == reftraj_ensemble) THEN
         reftraj_section => section_vals_get_subs_vals(md_section, "REFTRAJ")
         ALLOCATE (reftraj)
         CALL create_reftraj(reftraj, reftraj_section, para_env)
         CALL set_md_env(md_env, reftraj=reftraj)
      END IF

      CALL force_env_get(force_env, subsys=subsys, cell=cell, &
                         force_env_section=force_env_section)
      CALL cp_subsys_get(subsys, virial=virial)

      ! Set V0 if needed
      IF (simpar%ensemble == nph_uniaxial_ensemble .OR. simpar%ensemble == nph_uniaxial_damped_ensemble) THEN
         IF (simpar%v0 == 0._dp) simpar%v0 = cell%deth
      END IF

      ! Initialize velocities possibly applying constraints at the zeroth MD step
      CALL section_vals_val_get(motion_section, "PRINT%RESTART%SPLIT_RESTART_FILE", &
                                l_val=write_binary_restart_file)
      CALL setup_velocities(force_env, simpar, globenv, md_env, md_section, constraint_section, &
                            write_binary_restart_file)

      ! Setup Free Energy Calculation (if required)
      CALL fe_env_create(fe_env, free_energy_section)

      CALL set_md_env(md_env=md_env, simpar=simpar, fe_env=fe_env, cell=cell, &
                      force_env=force_env)
      IF (simpar%constraint) CALL initialize_moving_colvars(force_env, motion_section, simpar)

      ! Possibly initialize Wiener processes
      ![NB] Tested again within create_wiener_process.  Why??
      IF (need_per_atom_wiener_process(md_env)) CALL create_wiener_process(md_env)

      time_iter_start = m_walltime()

      CALL get_md_env(md_env, force_env=force_env, itimes=itimes, constant=constant, &
                      md_ener=md_ener, t=time, used_time=used_time)

      ! Attach the time counter of the meta_env to the one of the MD
      CALL set_meta_env(force_env%meta_env, time=time)

      ! Initialize the md_ener structure
      CALL initialize_md_ener(md_ener, force_env, simpar)

      ! Check for ensembles requiring the stress tensor - takes into account the possibility for
      ! multiple force_evals
      IF ((simpar%ensemble == npt_i_ensemble) .OR. &
          (simpar%ensemble == npt_ia_ensemble) .OR. &
          (simpar%ensemble == npt_f_ensemble) .OR. &
          (simpar%ensemble == npe_f_ensemble) .OR. &
          (simpar%ensemble == npe_i_ensemble) .OR. &
          (simpar%ensemble == nph_uniaxial_ensemble) .OR. &
          (simpar%ensemble == nph_uniaxial_damped_ensemble)) THEN
         check = virial%pv_availability
         IF (.NOT. check) THEN
            CALL cp_abort(__LOCATION__, &
                          "Virial evaluation not requested for this run in the input file!"// &
                          " You may consider to switch on the virial evaluation with the keyword: STRESS_TENSOR."// &
                          " Be sure the method you are using can compute the virial!")
         END IF
         IF (ASSOCIATED(force_env%sub_force_env)) THEN
            DO i = 1, SIZE(force_env%sub_force_env)
               IF (ASSOCIATED(force_env%sub_force_env(i)%force_env)) THEN
                  CALL force_env_get(force_env%sub_force_env(i)%force_env, subsys=subsys_i)
                  CALL cp_subsys_get(subsys_i, virial=virial)
                  check = check .AND. virial%pv_availability
               END IF
            END DO
         END IF
         IF (.NOT. check) THEN
            CALL cp_abort(__LOCATION__, &
                          "Virial evaluation not requested for all the force_eval sections present in"// &
                          " the input file! You have to switch on the virial evaluation with the keyword: STRESS_TENSOR"// &
                          " in each force_eval section. Be sure the method you are using can compute the virial!")
         END IF
      END IF

      ! Computing Forces at zero MD step
      IF (simpar%ensemble /= reftraj_ensemble) THEN
         CALL section_vals_val_get(md_section, "STEP_START_VAL", i_val=itimes)
         CALL section_vals_val_get(md_section, "TIME_START_VAL", r_val=time)
         CALL section_vals_val_get(md_section, "ECONS_START_VAL", r_val=constant)
         CALL cp_iterate(logger%iter_info, iter_nr=itimes)
         IF (save_mem) THEN
            work_section => section_vals_get_subs_vals(subsys_section, "VELOCITY")
            CALL section_vals_remove_values(work_section)
            work_section => section_vals_get_subs_vals(subsys_section, "SHELL_VELOCITY")
            CALL section_vals_remove_values(work_section)
            work_section => section_vals_get_subs_vals(subsys_section, "CORE_VELOCITY")
            CALL section_vals_remove_values(work_section)
         END IF

         IF (ehrenfest_md) THEN
            CALL rt_prop_setup(force_env)
            force_env%qs_env%rtp%dt = simpar%dt
         ELSE
            ![NB] Lets let all methods, even ones without consistent energies, succeed here.
            !     They'll fail in actual integrator if needed
            ! consistent_energies=.FALSE. by default
            CALL force_env_calc_energy_force(force_env, calc_force=.TRUE.)
         END IF
         ! ASSOCIATED rather than electrons_enabled: report-only mode needs this
         ! readback too -- computing the Fermi level and work function every step
         ! is the entire point of it. charge_mixing_enabled is false there, so it
         ! takes the plain single-force_env branch below.
         IF (ASSOCIATED(electrons)) THEN
            IF (electrons%charge_mixing_enabled) THEN
               CALL compute_mixed_workfunction(force_env, electrons%potential_average_axis, &
                                               electrons%potential_average_min, electrons%potential_average_max, &
                                               electrons%q, electrons%charge_mixing_span, &
                                               electrons%current_fermi_level_eV, electrons%current_workfunction, &
                                               lambda_in=lambda_mix)
            ELSE
               CALL qs_electrons_compute_workfunction(force_env, electrons%potential_average_axis, &
                                                      electrons%potential_average_min, electrons%potential_average_max, &
                                                      electrons%current_fermi_level_eV, electrons%current_workfunction)
            END IF
         END IF

         IF (ASSOCIATED(force_env%qs_env)) THEN
            CALL qs_env_time_update(force_env%qs_env, time, itimes)
         END IF
         ! Warm-up engines for metadynamics
         IF (ASSOCIATED(force_env%meta_env)) THEN
            ! Setup stuff for plumed if needed
            IF (force_env%meta_env%use_plumed .EQV. .TRUE.) THEN
               CALL metadyn_initialise_plumed(force_env, simpar, itimes)
               ! Include the bias in the first half kick and initial energy/pressure.
               ! This also initializes PLUMED's first-step bookkeeping on restart.
               CALL metadyn_integrator(force_env, itimes)
            ELSE
               IF (force_env%meta_env%langevin) THEN
                  CALL create_wiener_process_cv(force_env%meta_env)
               END IF
               IF (force_env%meta_env%well_tempered) THEN
                  force_env%meta_env%wttemperature = simpar%temp_ext
                  IF (force_env%meta_env%wtgamma > EPSILON(1._dp)) THEN
                     dummy = force_env%meta_env%wttemperature*(force_env%meta_env%wtgamma - 1._dp)
                     IF (force_env%meta_env%delta_t > EPSILON(1._dp)) THEN
                        check = ABS(force_env%meta_env%delta_t - dummy) < 1.E+3_dp*EPSILON(1._dp)
                        IF (.NOT. check) THEN
                           CALL cp_abort(__LOCATION__, &
                                         "Inconsistency between DELTA_T and WTGAMMA (both specified):"// &
                                         " please, verify that DELTA_T=(WTGAMMA-1)*TEMPERATURE")
                        END IF
                     ELSE
                        force_env%meta_env%delta_t = dummy
                     END IF
                  ELSE
                     force_env%meta_env%wtgamma = 1._dp &
                                                  + force_env%meta_env%delta_t/force_env%meta_env%wttemperature
                  END IF
                  force_env%meta_env%invdt = 1._dp/force_env%meta_env%delta_t
               END IF
               CALL metadyn_forces(force_env)
               CALL metadyn_write_colvar(force_env)
            END IF
         END IF

         IF (simpar%do_respa) THEN
            CALL force_env_calc_energy_force(force_env%sub_force_env(1)%force_env, &
                                             calc_force=.TRUE.)
         END IF

         CALL force_env_get(force_env, subsys=subsys)

         CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, local_particles=local_particles, &
                            particles=particles, virial=virial)

         CALL virial_evaluate(atomic_kinds%els, particles%els, local_particles, &
                              virial, force_env%para_env)

         CALL md_energy(md_env, md_ener)

         ! ASSOCIATED rather than electrons_enabled, so report-only mode writes
         ! its step-0 line too. There target_workfunction is 0, so legendre_eV
         ! vanishes and grand_canonical_eV reduces to E_pot.
         IF (ASSOCIATED(electrons)) THEN
            IF (output_unit > 0 .AND. electrons%cpcont_unit > 0) THEN
               ! Grand-canonical potential, Omega = E - mu*(q - N_neutral), with
               ! mu = -Phi_target in work-function control and mu = e_F_target in
               ! Fermi-level control.
               ! Column 9 carries the control target actually in use (the work function
               ! or the Fermi level); the header written at file-open names which one.
               IF (electrons%control_fermi) THEN
                  legendre_eV = -electrons%target_fermi_level* &
                                (electrons%q - electrons%neutral_electrons)
                  target_col = electrons%target_fermi_level
               ELSE
                  legendre_eV = electrons%target_workfunction* &
                                (electrons%q - electrons%neutral_electrons)
                  target_col = electrons%target_workfunction
               END IF
               grand_canonical_eV = md_ener%epot*evolt + legendre_eV
               WRITE (electrons%cpcont_unit, '(I10,9ES20.10E3)') itimes, time*femtoseconds, &
                  md_ener%epot*evolt, md_ener%constant*evolt, md_ener%temp_part, &
                  electrons%q, electrons%current_fermi_level_eV, electrons%current_workfunction, &
                  target_col, grand_canonical_eV
               CALL m_flush(electrons%cpcont_unit)
            END IF
         END IF

         IF (electrons_enabled) CALL store_electrons_restart(electrons_section, electrons)
         CALL md_write_output(md_env) !inits the print env at itimes == 0 also writes trajectories
         md_stride = 1
      ELSE
         CALL get_md_env(md_env, reftraj=reftraj)
         CALL initialize_reftraj(reftraj, reftraj_section, md_env)
         itimes = reftraj%info%first_snapshot - 1
         md_stride = reftraj%info%stride
         IF (ASSOCIATED(force_env%meta_env)) THEN
            IF (force_env%meta_env%use_plumed .EQV. .TRUE.) THEN
               CALL metadyn_initialise_plumed(force_env, simpar, itimes)
            END IF
         END IF
      END IF

      CALL cp_print_key_finished_output(simpar%info_constraint, logger, &
                                        constraint_section, "CONSTRAINT_INFO")
      CALL cp_print_key_finished_output(simpar%lagrange_multipliers, logger, &
                                        constraint_section, "LAGRANGE_MULTIPLIERS")

! if we need the initial kinetic energy for Hybrid Monte Carlo
      IF (PRESENT(hmc_e_initial)) hmc_e_initial = md_ener%ekin

      IF (itimes >= simpar%max_steps) CALL cp_abort(__LOCATION__, &
                                                    "maximum step number smaller than initial step value")

      ! Real MD Loop
      DO istep = 1, simpar%nsteps, md_stride
         ! Increase counters
         itimes = itimes + 1
         time = time + simpar%dt
         !needed when electric field fields are applied
         IF (ASSOCIATED(force_env%qs_env)) THEN
            CALL qs_env_time_update(force_env%qs_env, time, itimes)
         END IF
         IF (ehrenfest_md) force_env%qs_env%rtp%istep = istep

         IF (.NOT. logger%iter_info%last_iter(logger%iter_info%n_rlevel)) THEN
            CALL cp_iterate(logger%iter_info, last=(istep == simpar%nsteps), iter_nr=itimes)
         ELSE
            CALL cp_iterate(logger%iter_info, last=.TRUE., iter_nr=itimes)
         END IF

         IF (output_unit > 0 .AND. &
             BTEST(cp_print_key_should_output(logger%iter_info, md_section, &
                                              "PRINT%PROGRAM_RUN_INFO"), cp_p_file)) THEN
            WRITE (md_step_line, FMT="(A,I12,A,I0)") "MD STEP: ", istep, " / ", simpar%nsteps

            WRITE (UNIT=output_unit, FMT="(/,T2,A)") REPEAT("-", LEN_TRIM(md_step_line))
            WRITE (UNIT=output_unit, FMT="(T2,A)") TRIM(md_step_line)
            WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("-", LEN_TRIM(md_step_line))
            CALL m_flush(output_unit)
         END IF

         ! Open possible Shake output units
         simpar%info_constraint = cp_print_key_unit_nr(logger, constraint_section, "CONSTRAINT_INFO", &
                                                       extension=".shakeLog", log_filename=.FALSE.)
         simpar%lagrange_multipliers = cp_print_key_unit_nr( &
                                       logger, constraint_section, &
                                       "LAGRANGE_MULTIPLIERS", extension=".LagrangeMultLog", log_filename=.FALSE.)
         simpar%dump_lm = BTEST(cp_print_key_should_output(logger%iter_info, constraint_section, &
                                                           "LAGRANGE_MULTIPLIERS"), cp_p_file)

         ! Update temperature for thermal regions and thermostat regions
         CALL update_expected_temperature(md_env)

         ! Velocity Verlet Integrator
         CALL velocity_verlet(md_env, globenv)

         ! Close Shake output if requested...
         CALL cp_print_key_finished_output(simpar%info_constraint, logger, &
                                           constraint_section, "CONSTRAINT_INFO")
         CALL cp_print_key_finished_output(simpar%lagrange_multipliers, logger, &
                                           constraint_section, "LAGRANGE_MULTIPLIERS")

         ! Free Energy calculation
         CALL free_energy_evaluate(md_env, should_stop, free_energy_section)

         IF (should_stop) EXIT

         ! Test for <PROJECT_NAME>.EXIT_MD or for WALL_TIME to exit
         ! Default:
         ! IF so we don't overwrite the restart or append to the trajectory
         ! because the execution could in principle stop inside the SCF where energy
         ! and forces are not converged.
         ! But:
         ! You can force to print the last step (for example if the method used
         ! to compute energy and forces is not SCF based) activating the print_key
         ! MOTION%MD%PRINT%FORCE_LAST.
         CALL external_control(should_stop, "MD", globenv=globenv)

         !check if upper bound of total steps has been reached
         IF (.NOT. (istep == simpar%nsteps) .AND. logger%iter_info%last_iter(logger%iter_info%n_rlevel)) should_stop = .TRUE.
         IF (itimes >= simpar%max_steps) should_stop = .TRUE.

         ! call external hook e.g. from global optimization
         IF (PRESENT(mdctrl)) THEN
            CALL mdctrl_callback(mdctrl, md_env, should_stop)
         END IF

         IF (should_stop) THEN
            CALL cp_iterate(logger%iter_info, last=.TRUE., iter_nr=itimes)
            !In Ehrenfest molecular dynamics the external control is only checked after a converged propagation
            !The restart needs to be written in order to be consistent with the mos/density matrix for the restart
            IF (run_type_id == ehrenfest) THEN
               IF (electrons_enabled) CALL store_electrons_restart(electrons_section, electrons)
               CALL md_output(md_env, md_section, force_env%root_section, .FALSE.)
            ELSE
               IF (electrons_enabled) CALL store_electrons_restart(electrons_section, electrons)
               CALL md_output(md_env, md_section, force_env%root_section, should_stop)
            END IF
            EXIT
         END IF

         IF (simpar%ensemble /= reftraj_ensemble) THEN
            CALL md_energy(md_env, md_ener)
            CALL temperature_control(simpar, md_env, md_ener, force_env, logger)
            CALL comvel_control(md_ener, force_env, md_section, logger)
            CALL angvel_control(md_ener, force_env, md_section, logger)

            ! ASSOCIATED rather than electrons_enabled, so report-only mode
            ! writes a line every step as well.
            IF (ASSOCIATED(electrons)) THEN
               IF (output_unit > 0 .AND. electrons%cpcont_unit > 0) THEN
                  ! Column 9 is the control target -- see the comment at the
                  ! step-0 write above.
                  IF (electrons%control_fermi) THEN
                     legendre_eV = -electrons%target_fermi_level* &
                                   (electrons%q - electrons%neutral_electrons)
                     target_col = electrons%target_fermi_level
                  ELSE
                     legendre_eV = electrons%target_workfunction* &
                                   (electrons%q - electrons%neutral_electrons)
                     target_col = electrons%target_workfunction
                  END IF
                  grand_canonical_eV = md_ener%epot*evolt + legendre_eV
                  WRITE (electrons%cpcont_unit, '(I10,9ES20.10E3)') itimes, time*femtoseconds, &
                     md_ener%epot*evolt, md_ener%constant*evolt, md_ener%temp_part, &
                     electrons%q, electrons%current_fermi_level_eV, electrons%current_workfunction, &
                     target_col, grand_canonical_eV
                  CALL m_flush(electrons%cpcont_unit)
               END IF
            END IF
         ELSE
            CALL md_ener_reftraj(md_env, md_ener)
         END IF

         time_iter_stop = m_walltime()
         used_time = time_iter_stop - time_iter_start
         time_iter_start = time_iter_stop

         IF (electrons_enabled) CALL store_electrons_restart(electrons_section, electrons)
         CALL md_output(md_env, md_section, force_env%root_section, should_stop)
         IF (simpar%ensemble == reftraj_ensemble) THEN
            CALL write_output_reftraj(md_env)
         END IF
      END DO

      ! Guarded on ASSOCIATED: with no &ELECTRONS section at all, `electrons` is
      ! NULL and this dereferenced a null pointer on the ionode of every ordinary
      ! MD run. Also guarded on cpcont_unit > 0, since the unit keeps its -1
      ! default on every rank that never opened the file.
      IF (ASSOCIATED(electrons)) THEN
         IF (output_unit > 0 .AND. electrons%cpcont_unit > 0) THEN
            CALL close_file(electrons%cpcont_unit)
         END IF
      END IF

! if we need the final kinetic energy for Hybrid Monte Carlo
      IF (PRESENT(hmc_e_final)) hmc_e_final = md_ener%ekin

      ! Remove the iteration level
      CALL cp_rm_iter_level(logger, "MD")

      ! Clean up PLUMED
      IF (ASSOCIATED(force_env%meta_env)) THEN
         IF (force_env%meta_env%use_plumed .EQV. .TRUE.) THEN
            CALL metadyn_finalise_plumed()
         END IF
      END IF

      ! Deallocate Thermostats and Barostats
      CALL release_simpar_type(simpar)
      CALL timestop(handle)

   END SUBROUTINE qs_mol_dyn_low

! **************************************************************************************************
!> \brief Project initial or restarted velocities onto the upcoming target rate.
!> \param force_env ...
!> \param motion_section ...
!> \param simpar ...
! **************************************************************************************************
   SUBROUTINE initialize_moving_colvars(force_env, motion_section, simpar)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(section_vals_type), POINTER                   :: motion_section
      TYPE(simpar_type), POINTER                         :: simpar

      INTEGER                                            :: i
      LOGICAL                                            :: moving
      REAL(KIND=dp)                                      :: dt
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(global_constraint_type), POINTER              :: gci
      TYPE(molecule_kind_type), POINTER                  :: kinds(:)

      CALL force_env_get(force_env, subsys=subsys)
      CALL cp_subsys_get(subsys, molecule_kind_set=kinds, gci=gci)
      dt = simpar%dt
      IF (simpar%do_respa) dt = dt/simpar%n_time_steps
      moving = .FALSE.
      DO i = 1, SIZE(kinds)
         IF (.NOT. ASSOCIATED(kinds(i)%colv_list)) CYCLE
         CALL shake_update_colv_low(kinds(i)%colv_list, dt, motion_section, initialize=.TRUE.)
         moving = moving .OR. ANY(kinds(i)%colv_list%expected_value_growth_speed /= 0.0_dp)
      END DO
      IF (ASSOCIATED(gci%colv_list)) THEN
         CALL shake_update_colv_low(gci%colv_list, dt, motion_section, initialize=.TRUE.)
         moving = moving .OR. ANY(gci%colv_list%expected_value_growth_speed /= 0.0_dp)
      END IF
      gci%joint_colvar_checks = moving
      IF (moving) CALL force_env_rattle(force_env, dt=dt, shake_tol=simpar%shake_tol, reset=.TRUE.)
   END SUBROUTINE initialize_moving_colvars

! **************************************************************************************************
!> \brief Write the live state of the electronic extended variable back into the
!>        &ELECTRONS input section, so that CP2K's ordinary restart dump carries
!>        it and a continued run resumes the electronic dynamics instead of
!>        silently resetting q to Q_INIT with p = 0 mid-trajectory.
!>
!>        Called immediately before every md_output / md_write_output, which are
!>        the only places a restart file is produced. Deliberately NOT guarded by
!>        output_unit > 0: section_vals_val_set mutates the replicated input tree,
!>        so every rank must apply the same update or the trees diverge.
!>
!>        Q_INIT doubles as the restart slot for q, in the same way CP2K reuses
!>        the input coordinates as the restart slot for the positions. The units
!>        match the getters in qs_mol_dyn_low: q and p in atomic units,
!>        thermostat_energy in Hartree, Fermi level and work function in eV.
!> \param electrons_section the &MOTION%MD%ELECTRONS section
!> \param electrons the live electronic state
! **************************************************************************************************
   SUBROUTINE store_electrons_restart(electrons_section, electrons)
      TYPE(section_vals_type), POINTER                   :: electrons_section
      TYPE(qs_electrons_type), POINTER                   :: electrons

      IF (.NOT. ASSOCIATED(electrons_section)) RETURN
      IF (.NOT. ASSOCIATED(electrons)) RETURN

      CALL section_vals_val_set(electrons_section, "Q_INIT", r_val=electrons%q)
      CALL section_vals_val_set(electrons_section, "P_INIT", r_val=electrons%p)
      CALL section_vals_val_set(electrons_section, "POTENTIOSTAT_ENERGY", &
                                r_val=electrons%thermostat_energy)
      CALL section_vals_val_set(electrons_section, "CURRENT_FERMI_LEVEL", &
                                r_val=electrons%current_fermi_level_eV)
      CALL section_vals_val_set(electrons_section, "CURRENT_WORKFUNCTION", &
                                r_val=electrons%current_workfunction)

   END SUBROUTINE store_electrons_restart

END MODULE md_run

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Methods for storing MD parameters type
!> \author CJM
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
!>         reorganization of the original routines/modules
! **************************************************************************************************
MODULE simpar_methods
   USE bibliography,                    ONLY: Evans1983,&
                                              Kuhne2007,&
                                              Minary2003,&
                                              Rengaraj2020,&
                                              Ricci2003,&
                                              cite_reference
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_unit_nr
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE input_constants,                 ONLY: &
        isokin_ensemble, langevin_ensemble, npe_f_ensemble, npe_i_ensemble, &
        nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, npt_i_ensemble, &
        nvt_ensemble, reftraj_ensemble
   USE input_cp2k_md,                   ONLY: create_md_section
   USE input_enumeration_types,         ONLY: enum_i2c,&
                                              enumeration_type
   USE input_keyword_types,             ONLY: keyword_get,&
                                              keyword_type
   USE input_section_types,             ONLY: section_get_keyword,&
                                              section_release,&
                                              section_type,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE simpar_types,                    ONLY: simpar_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'simpar_methods'
   PUBLIC :: read_md_section

CONTAINS

! **************************************************************************************************
!> \brief Reads the MD section and setup the simulation parameters type
!> \param simpar ...
!> \param motion_section ...
!> \param md_section ...
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE read_md_section(simpar, motion_section, md_section)
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(section_vals_type), POINTER                   :: motion_section, md_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_md_section', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: iprint, iw
      REAL(kind=dp)                                      :: tmp_r1, tmp_r2, tmp_r3
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(enumeration_type), POINTER                    :: enum
      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: section
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (logger, print_key, enum, keyword, section)
      logger => cp_get_default_logger()
      iw = cp_print_key_unit_nr(logger, md_section, "PRINT%PROGRAM_RUN_INFO", extension=".log")

      CALL read_md_low(simpar, motion_section, md_section)
      IF (iw > 0) WRITE (iw, *)

      ! Begin setup Langevin dynamics
      IF (simpar%ensemble == langevin_ensemble) THEN
         CALL cite_reference(Ricci2003)
         IF (simpar%noisy_gamma > 0.0_dp) CALL cite_reference(Kuhne2007)
         IF (simpar%shadow_gamma > 0.0_dp) CALL cite_reference(Rengaraj2020)
         ! Normalization factor using a normal Gaussian random number distribution
         simpar%var_w = 2.0_dp*simpar%temp_ext*simpar%dt*(simpar%gamma + simpar%noisy_gamma)
         IF (iw > 0) THEN
            tmp_r1 = cp_unit_from_cp2k(simpar%gamma, "fs^-1")
            tmp_r2 = cp_unit_from_cp2k(simpar%noisy_gamma, "fs^-1")
            tmp_r3 = cp_unit_from_cp2k(simpar%shadow_gamma, "fs^-1")
            WRITE (UNIT=iw, FMT="(T2,A,T71,ES10.3)") &
               "LD| Gamma [1/fs] ", tmp_r1, &
               "LD| Noisy Gamma [1/fs]", tmp_r2, &
               "LD| Shadow Gamma [1/fs]", tmp_r3, &
               "LD| Variance [a.u.]", simpar%var_w
         END IF
      END IF

      ! Create section for output enumeration infos
      CALL create_md_section(section)
      keyword => section_get_keyword(section, "ENSEMBLE")
      CALL keyword_get(keyword, enum=enum)

      !..write some information to output
      IF (iw > 0) THEN
         WRITE (iw, '( A )') ' MD| Molecular Dynamics Protocol '
         WRITE (iw, '( A,T61,A20)') ' MD| Ensemble Type ', &
            ADJUSTR(TRIM(enum_i2c(enum, simpar%ensemble)))
         WRITE (iw, '( A,T71,I10 )') ' MD| Number of Time Steps ', &
            simpar%nsteps
         IF (simpar%variable_dt) THEN
            WRITE (iw, '( A )') ' MD| Variable Time Step is activated '
            tmp_r1 = cp_unit_from_cp2k(simpar%dt, "fs")
            WRITE (iw, '( A,A2,A,T71,F10.2 )') ' MD| Max. Time Step [', 'fs', '] ', tmp_r1
            tmp_r1 = cp_unit_from_cp2k(simpar%dr_tol, "angstrom")
            WRITE (iw, '( A, T71, F10.4 )') ' MD| Max. atomic displacement permitted [A] ', tmp_r1
         ELSE
            tmp_r1 = cp_unit_from_cp2k(simpar%dt, "fs")
            WRITE (iw, '( A,A2,A,T71,F10.2 )') ' MD| Time Step [', 'fs', '] ', tmp_r1
         END IF
         tmp_r1 = cp_unit_from_cp2k(simpar%temp_ext, "K")
         WRITE (iw, '( A,T71,F10.2 )') ' MD| Temperature [K] ', tmp_r1
         tmp_r1 = cp_unit_from_cp2k(simpar%temp_tol, "K")
         WRITE (iw, '( A,T71,F10.2 )') ' MD| Temperature tolerance [K] ', tmp_r1

         IF (simpar%annealing) &
            WRITE (iw, '( A,T71,F10.2 )') ' MD| Annealing ion factor      ', &
            simpar%f_annealing
         IF ((simpar%ensemble == npe_f_ensemble .OR. simpar%ensemble == npe_i_ensemble) &
             .AND. simpar%annealing_cell) &
            WRITE (iw, '( A,T71,F10.2 )') ' MD| Annealing cell factor      ', &
            simpar%f_annealing_cell
         IF (simpar%ensemble == npt_i_ensemble .OR. &
             simpar%ensemble == npt_f_ensemble .OR. &
             simpar%ensemble == npe_i_ensemble .OR. &
             simpar%ensemble == npe_f_ensemble) THEN
            tmp_r1 = cp_unit_from_cp2k(simpar%p_ext, "bar")
            WRITE (iw, '( A,A3,A, T71, F10.2 )') &
               ' MD| Pressure [', 'Bar', '] ', tmp_r1
            tmp_r1 = cp_unit_from_cp2k(simpar%tau_cell, "fs")
            WRITE (iw, '( A,A4,A, T71, F10.2 )') &
               ' MD| Barostat time constant [', 'fs', '] ', tmp_r1
         END IF
         IF (simpar%ensemble == isokin_ensemble) THEN
            CALL cite_reference(Evans1983)
            CALL cite_reference(Minary2003)
            WRITE (iw, '( A )') ' MD| Simulation in the isokinetic ensemble'
         END IF
         IF (simpar%constraint) THEN
            WRITE (iw, '( A )') ' MD| Constraints activated '
            WRITE (iw, '( A,T71,G10.4 )') ' MD| Tolerance for shake ', &
               simpar%shake_tol
         END IF

         print_key => section_vals_get_subs_vals(motion_section, "MD%PRINT%PROGRAM_RUN_INFO")
         CALL section_vals_val_get(print_key, "EACH%MD", i_val=iprint)
         WRITE (iw, '( A,T63,i10,A )') ' MD| Print MD information every', iprint, ' step(s)'
         WRITE (iw, '( A,T20,A,T71,A10 )') ' MD| File type', 'Print frequency[steps]', 'File names'

         print_key => section_vals_get_subs_vals(motion_section, "PRINT%TRAJECTORY")
         CALL section_vals_val_get(print_key, "EACH%MD", i_val=iprint)
         filename = cp_print_key_generate_filename(logger, print_key, &
                                                   extension=".xyz", middle_name="pos", my_local=.FALSE.)
         WRITE (iw, '( A,T20,i10,T31,A50 )') ' MD| Coordinates', iprint, &
            ADJUSTR(TRIM(filename))

         IF (simpar%ensemble == nph_uniaxial_ensemble .OR. &
             simpar%ensemble == nph_uniaxial_damped_ensemble .OR. &
             simpar%ensemble == npt_i_ensemble .OR. &
             simpar%ensemble == npt_f_ensemble .OR. &
             simpar%ensemble == npe_i_ensemble .OR. &
             simpar%ensemble == npe_f_ensemble) THEN

            print_key => section_vals_get_subs_vals(motion_section, "PRINT%CELL")
            CALL section_vals_val_get(print_key, "EACH%MD", i_val=iprint)
            filename = cp_print_key_generate_filename(logger, print_key, &
                                                      extension=".cell", my_local=.FALSE.)
            WRITE (iw, '( A,T20,i10,T31,A50 )') ' MD| Simulation Cell', iprint, &
               ADJUSTR(TRIM(filename))
         END IF

         print_key => section_vals_get_subs_vals(motion_section, "PRINT%VELOCITIES")
         CALL section_vals_val_get(print_key, "EACH%MD", i_val=iprint)
         filename = cp_print_key_generate_filename(logger, print_key, &
                                                   extension=".xyz", middle_name="vel", my_local=.FALSE.)
         WRITE (iw, '( A,T20,i10,T31,A50 )') ' MD| Velocities', iprint, &
            ADJUSTR(TRIM(filename))

         print_key => section_vals_get_subs_vals(motion_section, "MD%PRINT%ENERGY")
         CALL section_vals_val_get(print_key, "EACH%MD", i_val=iprint)
         filename = cp_print_key_generate_filename(logger, print_key, &
                                                   extension=".ener", my_local=.FALSE.)
         WRITE (iw, '( A,T20,i10,T31,A50 )') ' MD| Energies', iprint, &
            ADJUSTR(TRIM(filename))

         print_key => section_vals_get_subs_vals(motion_section, "PRINT%RESTART")
         CALL section_vals_val_get(print_key, "EACH%MD", i_val=iprint)
         filename = cp_print_key_generate_filename(logger, print_key, &
                                                   extension=".restart", my_local=.FALSE.)
         WRITE (iw, '( A,T20,i10,T31,A50 )') ' MD| Dump', iprint, &
            ADJUSTR(TRIM(filename))

         WRITE (iw, *)
         IF (simpar%ensemble == nph_uniaxial_ensemble .OR. &
             simpar%ensemble == nph_uniaxial_damped_ensemble) THEN
            WRITE (iw, '( A )') ' SHOCK| Uniaxial Shock Parameters: '
            tmp_r1 = cp_unit_from_cp2k(simpar%v_shock, "m*s^-1")
            WRITE (iw, '( A,A4,A,T71,F10.4 )') &
               ' SHOCK| Shock Velocity [', 'm/s', '] ', tmp_r1
            tmp_r1 = cp_unit_from_cp2k(simpar%gamma_nph, "fs^-1")
            WRITE (iw, '( A,A4,A,T71,F10.4 )') &
               ' SHOCK| Damping Coefficient [', '1/fs', '] ', tmp_r1
            tmp_r1 = cp_unit_from_cp2k(simpar%p0, "bar")
            WRITE (iw, '( A,A3,A, T71, F10.2 )') &
               ' SHOCK| Pressure [', 'Bar', '] ', tmp_r1
            WRITE (iw, '( A,A4,A, T71, E10.4 )') &
               ' SHOCK| Barostat Mass [', 'au', '] ', simpar%cmass
         END IF
         ! Print warning for temp_tol
         IF (simpar%temp_tol > 0.0_dp) THEN
            CALL cp_warn(__LOCATION__, &
                         "A temperature tolerance (TEMP_TOL) is used during the MD. "// &
                         "Due to the velocity rescaling algorithm jumps may appear in the conserved quantity.")
         END IF
         ! Print warning for annealing
         IF (simpar%annealing) THEN
            IF (simpar%ensemble == nvt_ensemble .OR. &
                simpar%ensemble == npt_i_ensemble .OR. &
                simpar%ensemble == npt_f_ensemble) THEN
               CALL cp_abort(__LOCATION__, &
                             "Annealing of the ions has been required "// &
                             "even if the thermostat is active (nvt or npt_i or npt_f) "// &
                             "These two methods to control the temperature act one against the other.")
            END IF
         END IF
         ! Print warning for variable time step
         IF (simpar%variable_dt) THEN
            IF (simpar%ensemble == langevin_ensemble .OR. &
                simpar%ensemble == reftraj_ensemble .OR. &
                simpar%do_respa) THEN
               CALL cp_warn( &
                  __LOCATION__, &
                  "The variable timestep  has been required, however "// &
                  "this option is not available either with the Langevin ensemble or with the multiple timestep schme. "// &
                  "The run will proceed with constant timestep, as read from input.")
            END IF
         END IF
      END IF
      CALL section_release(section)
      CALL cp_print_key_finished_output(iw, logger, md_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

   END SUBROUTINE read_md_section

! **************************************************************************************************
!> \brief Low Level: Parses the MD input section
!> \param simpar ...
!> \param motion_section ...
!> \param md_section ...
!> \author teo
! **************************************************************************************************
   SUBROUTINE read_md_low(simpar, motion_section, md_section)
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(section_vals_type), POINTER                   :: motion_section, md_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_md_low', routineP = moduleN//':'//routineN

      LOGICAL                                            :: explicit
      TYPE(section_vals_type), POINTER                   :: tmp_section

      NULLIFY (tmp_section)
      CALL section_vals_val_get(md_section, "ENSEMBLE", i_val=simpar%ensemble)
      CALL section_vals_val_get(md_section, "STEPS", i_val=simpar%nsteps)
      CALL section_vals_val_get(md_section, "MAX_STEPS", i_val=simpar%max_steps)
      CALL section_vals_val_get(md_section, "TEMPERATURE", r_val=simpar%temp_ext)
      CALL section_vals_val_get(md_section, "TEMP_TOL", r_val=simpar%temp_tol)
      CALL section_vals_val_get(md_section, "ANGVEL_ZERO", l_val=simpar%angvel_zero)
      CALL section_vals_val_get(md_section, "TEMP_KIND", l_val=simpar%temperature_per_kind)
      CALL section_vals_val_get(md_section, "SCALE_TEMP_KIND", l_val=simpar%scale_temperature_per_kind)
      CALL section_vals_val_get(md_section, "ANNEALING", r_val=simpar%f_annealing, explicit=simpar%annealing)
      CALL section_vals_val_get(md_section, "ANNEALING_CELL", r_val=simpar%f_annealing_cell, &
                                explicit=simpar%annealing_cell)
      CALL section_vals_val_get(md_section, "TEMPERATURE_ANNEALING", r_val=simpar%f_temperature_annealing, &
                                explicit=simpar%temperature_annealing)
      CALL section_vals_val_get(md_section, "DISPLACEMENT_TOL", r_val=simpar%dr_tol, &
                                explicit=simpar%variable_dt)
      CALL section_vals_val_get(md_section, "TIMESTEP", r_val=simpar%dt)
      CALL section_vals_val_get(md_section, "INITIALIZATION_METHOD", &
                                i_val=simpar%initialization_method)
      ! Initialize dt_fact to 1.0
      simpar%dt_fact = 1.0_dp

      IF (simpar%ensemble == langevin_ensemble) THEN
         CALL section_vals_val_get(md_section, "LANGEVIN%GAMMA", r_val=simpar%gamma)
         CALL section_vals_val_get(md_section, "LANGEVIN%NOISY_GAMMA", r_val=simpar%noisy_gamma)
         CALL section_vals_val_get(md_section, "LANGEVIN%SHADOW_GAMMA", r_val=simpar%shadow_gamma)
      END IF

      tmp_section => section_vals_get_subs_vals(motion_section, "CONSTRAINT")
      CALL section_vals_get(tmp_section, explicit=simpar%constraint)
      IF (simpar%constraint) THEN
         CALL section_vals_val_get(tmp_section, "SHAKE_TOLERANCE", r_val=simpar%shake_tol)
         IF (simpar%shake_tol <= EPSILON(0.0_dp)*1000.0_dp) &
            CALL cp_warn(__LOCATION__, &
                         "Shake tolerance lower than 1000*EPSILON, where EPSILON is the machine precision. "// &
                         "This may lead to numerical problems. Setting up shake_tol to 1000*EPSILON!")
         simpar%shake_tol = MAX(EPSILON(0.0_dp)*1000.0_dp, simpar%shake_tol)

         CALL section_vals_val_get(tmp_section, "ROLL_TOLERANCE", r_val=simpar%roll_tol)
         IF (simpar%roll_tol <= EPSILON(0.0_dp)*1000.0_dp) &
            CALL cp_warn(__LOCATION__, &
                         "Roll tolerance lower than 1000*EPSILON, where EPSILON is the machine precision. "// &
                         "This may lead to numerical problems. Setting up roll_tol to 1000*EPSILON!")
         simpar%roll_tol = MAX(EPSILON(0.0_dp)*1000.0_dp, simpar%roll_tol)
      END IF

      IF (simpar%ensemble == nph_uniaxial_ensemble .OR. simpar%ensemble == nph_uniaxial_damped_ensemble) THEN
         tmp_section => section_vals_get_subs_vals(md_section, "MSST")
         CALL section_vals_val_get(tmp_section, "PRESSURE", r_val=simpar%p0)
         CALL section_vals_val_get(tmp_section, "ENERGY", r_val=simpar%e0)
         CALL section_vals_val_get(tmp_section, "VOLUME", r_val=simpar%v0)
         CALL section_vals_val_get(tmp_section, "GAMMA", r_val=simpar%gamma_nph)
         IF (simpar%gamma_nph /= 0.0_dp) simpar%ensemble = nph_uniaxial_damped_ensemble
         CALL section_vals_val_get(tmp_section, "CMASS", r_val=simpar%cmass)
         CALL section_vals_val_get(tmp_section, "VSHOCK", r_val=simpar%v_shock)
      END IF

      SELECT CASE (simpar%ensemble)
      CASE (nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, &
            npt_f_ensemble, npt_i_ensemble, npe_f_ensemble, npe_i_ensemble)
         tmp_section => section_vals_get_subs_vals(md_section, "BAROSTAT")
         CALL section_vals_val_get(tmp_section, "PRESSURE", r_val=simpar%p_ext)
         CALL section_vals_val_get(tmp_section, "TIMECON", r_val=simpar%tau_cell)
      END SELECT

      ! RESPA
      tmp_section => section_vals_get_subs_vals(md_section, "RESPA")
      CALL section_vals_get(tmp_section, explicit=simpar%do_respa)
      CALL section_vals_val_get(tmp_section, "FREQUENCY", i_val=simpar%n_time_steps)
      simpar%multi_time_switch = simpar%do_respa

      ! CORE-SHELL MODEL
      tmp_section => section_vals_get_subs_vals(md_section, "SHELL")
      CALL section_vals_val_get(tmp_section, "TEMPERATURE", r_val=simpar%temp_sh_ext)
      CALL section_vals_val_get(tmp_section, "TEMP_TOL", r_val=simpar%temp_sh_tol)

      CALL section_vals_val_get(tmp_section, "DISPLACEMENT_SHELL_TOL", r_val=simpar%dsc_tol, &
                                explicit=explicit)
      simpar%variable_dt = simpar%variable_dt .OR. explicit
      ! ADIABATIC DYNAMICS
      tmp_section => section_vals_get_subs_vals(md_section, "ADIABATIC_DYNAMICS")
      CALL section_vals_val_get(tmp_section, "TEMP_FAST", r_val=simpar%temp_fast)
      CALL section_vals_val_get(tmp_section, "TEMP_SLOW", r_val=simpar%temp_slow)
      CALL section_vals_val_get(tmp_section, "TEMP_TOL_FAST", r_val=simpar%temp_tol_fast)
      CALL section_vals_val_get(tmp_section, "TEMP_TOL_SLOW", r_val=simpar%temp_tol_slow)
      CALL section_vals_val_get(tmp_section, "N_RESP_FAST", i_val=simpar%n_resp_fast)

      ! VELOCITY SOFTENING
      tmp_section => section_vals_get_subs_vals(md_section, "VELOCITY_SOFTENING")
      CALL section_vals_val_get(tmp_section, "STEPS", i_val=simpar%soften_nsteps)
      CALL section_vals_val_get(tmp_section, "ALPHA", r_val=simpar%soften_alpha)
      CALL section_vals_val_get(tmp_section, "DELTA", r_val=simpar%soften_delta)
   END SUBROUTINE read_md_low

END MODULE simpar_methods

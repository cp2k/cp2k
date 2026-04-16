!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for Molecular Dynamics
!> \author Teodoro Laino [tlaino] - University of Zurich - 09.2007
! **************************************************************************************************
MODULE md_util

   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_generate_filename
   USE cp_units,                        ONLY: cp_unit_from_cp2k,&
                                              cp_unit_to_cp2k
   USE fparser,                         ONLY: EvalErrType,&
                                              evalf,&
                                              finalizef,&
                                              initf,&
                                              parsef
   USE input_constants,                 ONLY: &
        do_region_thermal, do_thermo_al, do_thermo_csvr, do_thermo_nose, npt_f_ensemble, &
        npt_i_ensemble, npt_ia_ensemble, nvt_adiabatic_ensemble, nvt_ensemble
   USE input_cp2k_restarts,             ONLY: write_restart
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE md_energies,                     ONLY: md_write_output
   USE md_environment_types,            ONLY: get_md_env,&
                                              md_environment_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE simpar_types,                    ONLY: simpar_type
   USE string_utilities,                ONLY: integer_to_string
   USE thermal_region_types,            ONLY: thermal_region_type,&
                                              thermal_regions_type
   USE thermostat_types,                ONLY: thermostat_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_util'

   PUBLIC :: md_output, &
             read_vib_eigs_unformatted, &
             update_expected_temperature

CONTAINS

! **************************************************************************************************
!> \brief update_expected_temperature
!> \param md_env ...
!> \par History
!>      03.2026 created with functionality to update thermal regions
!>      04.2026 updated with functionality to update thermostat regions
!> \author HE Zilong
! **************************************************************************************************
   SUBROUTINE update_expected_temperature(md_env)
      TYPE(md_environment_type), POINTER                 :: md_env

      CHARACTER(LEN=default_string_length)               :: my_itimes, my_region
      INTEGER                                            :: i, ireg
      INTEGER, POINTER                                   :: itimes
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(thermal_region_type), POINTER                 :: t_region
      TYPE(thermal_regions_type), POINTER                :: thermal_regions
      TYPE(thermostat_type), POINTER                     :: thermostat_part

      NULLIFY (itimes, simpar)
      CALL get_md_env(md_env, itimes=itimes, simpar=simpar)

      IF (simpar%do_thermal_region) THEN
         ! 1. Update thermal regions
         ! Note on units: t_region%temperature is updated by
         ! scale_velocity_region() in md_vel_utils.F and is
         ! in Kelvin, but t_region%temp_expected is in the
         ! cp2k internal units.
         ! See also: print_thermal_regions_temperature()
         ! from thermal_region_utils.F.
         NULLIFY (thermal_regions)
         CALL get_md_env(md_env, thermal_regions=thermal_regions)
         DO ireg = 1, thermal_regions%nregions
            NULLIFY (t_region)
            t_region => thermal_regions%thermal_region(ireg)
            IF (LEN_TRIM(t_region%temperature_function) > 0) THEN
               CALL initf(1)
               CALL parsef(1, TRIM(t_region%temperature_function), t_region%temperature_function_parameters)
               t_region%temperature_function_values(1) = itimes*1.0_dp
               t_region%temperature_function_values(2) = t_region%temperature
               t_region%temp_expected = cp_unit_to_cp2k(evalf(1, t_region%temperature_function_values), "K")
               CPASSERT(EvalErrType <= 0)
               CALL finalizef()
               IF (t_region%temp_expected < 0.0_dp) THEN
                  CALL integer_to_string(ireg, my_region)
                  CALL integer_to_string(itimes, my_itimes)
                  CALL cp_abort(__LOCATION__, &
                                "At step "//TRIM(my_itimes)//" , the temperature "// &
                                "evaluated for thermal region "//TRIM(my_region)//" is "// &
                                TRIM(ADJUSTL(cp_to_string(cp_unit_from_cp2k(t_region%temp_expected, "K"))))// &
                                " K, which is negative and unphysical!")
               END IF
            END IF
         END DO

         ! 2. Update thermostat regions by thermal regions
         ! Number of thermostat regions is the same as that of thermal regions
         ! if every particle is assigned to a thermal region, or one more than
         ! that of thermal regions if one or more particles are not in any of
         ! the thermal regions (which the extra thermostat adjusts simply to
         ! the constant target temperature simpar%temp_ext).
         ! The key variable is al%nvt(i)%nkt (see al_system_mapping.F) or csvr%nvt(i)%nkt
         ! (see csvr_system_mapping.F) or nhc%nvt(:, i)%nkt (see extended_system_mapping.F)
         ! where i loops over number of thermostats; previously ther are set by a fixed
         ! simpar%temp_ext in cp2k internal units. The thermostat algorithms use it to
         ! determine map_info%v_scale(imap) with imap = map_info%map_index(i).
         ! In thermostat_mapping_region_low(), it is referenced by pointer
         ! map_info%p_scale(ii, jj)%point => map_info%v_scale(imap)
         ! And finally in vel_rescale_particles() the particle velocities are updated.
         ! See also: print_thermostat_status() from thermostat_utils.F.
         IF ((simpar%ensemble == npt_f_ensemble) .OR. &
             (simpar%ensemble == npt_ia_ensemble) .OR. &
             (simpar%ensemble == npt_i_ensemble) .OR. &
             (simpar%ensemble == nvt_ensemble) .AND. &
             (.NOT. simpar%ensemble == nvt_adiabatic_ensemble)) THEN
            NULLIFY (thermostat_part)
            CALL get_md_env(md_env=md_env, thermostat_part=thermostat_part)
            SELECT CASE (thermostat_part%type_of_thermostat)
            CASE (do_thermo_al)
               CPASSERT(ASSOCIATED(thermostat_part%al))
               IF (thermostat_part%al%region == do_region_thermal) THEN
                  DO ireg = 1, thermal_regions%nregions
                     CPASSERT(ireg <= thermostat_part%al%loc_num_al) ! Does not work
                     thermostat_part%al%nvt(ireg)%nkt = &
                        thermal_regions%thermal_region(ireg)%temp_expected*thermostat_part%al%nvt(ireg)%fac
                  END DO
               END IF
            CASE (do_thermo_csvr)
               CPASSERT(ASSOCIATED(thermostat_part%csvr))
               IF (thermostat_part%csvr%region == do_region_thermal) THEN
                  DO ireg = 1, thermal_regions%nregions
                     CPASSERT(ireg <= thermostat_part%csvr%loc_num_csvr) ! Does not work
                     thermostat_part%csvr%nvt(ireg)%nkt = &
                        thermal_regions%thermal_region(ireg)%temp_expected*thermostat_part%csvr%nvt(ireg)%fac
                  END DO
               END IF
            CASE (do_thermo_nose)
               CPASSERT(ASSOCIATED(thermostat_part%nhc))
               IF (thermostat_part%nhc%region == do_region_thermal) THEN
                  DO ireg = 1, thermal_regions%nregions
                     CPASSERT(ireg <= thermostat_part%nhc%loc_num_nhc) ! Does not work
                     thermostat_part%nhc%nvt(1, ireg)%nkt = &
                        thermal_regions%thermal_region(ireg)%temp_expected*thermostat_part%nhc%nvt(1, ireg)%fac
                     DO i = 2, thermostat_part%nhc%nhc_len
                        thermostat_part%nhc%nvt(i, ireg)%nkt = &
                           thermal_regions%thermal_region(ireg)%temp_expected
                     END DO
                  END DO
               END IF
            CASE DEFAULT
               ! Also includes GLE which is supposed to work with REGION MASSIVE
               CALL cp_abort(__LOCATION__, &
                             "The selected type of thermostat does not support "// &
                             "a variable target temperature that is synchronized "// &
                             "between thermal regions and thermostat regions.")
            END SELECT
         END IF
      END IF

   END SUBROUTINE update_expected_temperature

! **************************************************************************************************
!> \brief collects the part of the MD that, basically, does the output
!> \param md_env ...
!> \param md_section ...
!> \param root_section ...
!> \param forced_io ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE md_output(md_env, md_section, root_section, forced_io)
      TYPE(md_environment_type), POINTER                 :: md_env
      TYPE(section_vals_type), POINTER                   :: md_section, root_section
      LOGICAL, INTENT(IN)                                :: forced_io

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'md_output'

      INTEGER                                            :: handle
      LOGICAL                                            :: do_print
      TYPE(section_vals_type), POINTER                   :: print_section

      CALL timeset(routineN, handle)
      do_print = .TRUE.
      IF (forced_io) THEN
         print_section => section_vals_get_subs_vals(md_section, "PRINT")
         CALL section_vals_val_get(print_section, "FORCE_LAST", l_val=do_print)
      END IF
      IF (do_print) THEN
         ! Dumps all files related to the MD run
         CALL md_write_output(md_env)
         CALL write_restart(md_env=md_env, root_section=root_section)
      END IF
      CALL timestop(handle)

   END SUBROUTINE md_output

! **************************************************************************************************
!> \brief read eigenvalues and eigenvectors of Hessian from vibrational analysis results, for use
!>        of initialising MD simulations. Expects to read an unformatted binary file
!> \param md_section : input section object containing MD subsections and keywords. This should
!>                     provide the filename to read vib analysis eigenvalues and eigenvectors.
!>                     If the filename is not explicitly specified by the user in the input, then
!>                     it will use the default CARTESIAN_EIGS print key filename defined in the
!>                     vibrational analysis input section as the filename.
!> \param vib_section : input section object containing vibrational analysis subsections
!>                      and keywords
!> \param para_env : cp2k mpi environment object, needed for IO in parallel computations
!> \param dof : outputs the total number of eigenvalues (no. degrees of freedom) read from the file
!> \param eigenvalues : outputs the eigenvalues (Cartesian frequencies) read from the file
!> \param eigenvectors : outputs the corresponding eigenvectors read from the file
!> \author Lianheng Tong, lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   SUBROUTINE read_vib_eigs_unformatted(md_section, &
                                        vib_section, &
                                        para_env, &
                                        dof, &
                                        eigenvalues, &
                                        eigenvectors)
      TYPE(section_vals_type), POINTER                   :: md_section, vib_section
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(OUT)                               :: dof
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: eigenvectors

      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: jj, n_rep_val, unit_nr
      LOGICAL                                            :: exist
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      logger => cp_get_default_logger()
      dof = 0
      eigenvalues = 0.0_dp
      eigenvectors = 0.0_dp
      ! obtain file name
      CALL section_vals_val_get(md_section, "INITIAL_VIBRATION%VIB_EIGS_FILE_NAME", &
                                n_rep_val=n_rep_val)
      IF (n_rep_val > 0) THEN
         CALL section_vals_val_get(md_section, "INITIAL_VIBRATION%VIB_EIGS_FILE_NAME", c_val=filename)
      ELSE
         print_key => section_vals_get_subs_vals(vib_section, "PRINT%CARTESIAN_EIGS")
         filename = cp_print_key_generate_filename(logger, print_key, extension="eig", &
                                                   my_local=.FALSE.)
      END IF
      ! read file
      IF (para_env%is_source()) THEN
         INQUIRE (FILE=filename, exist=exist)
         IF (.NOT. exist) THEN
            CPABORT("File "//filename//" is not found.")
         END IF
         CALL open_file(file_name=filename, &
                        file_action="READ", &
                        file_form="UNFORMATTED", &
                        file_status="OLD", &
                        unit_number=unit_nr)
         ! the first record contains one integer giving degrees of freedom
         READ (unit_nr) dof
         IF (dof > SIZE(eigenvalues)) THEN
            CPABORT("Too many DoFs found in "//filename)
         END IF
         ! the second record contains the eigenvalues
         READ (unit_nr) eigenvalues(1:dof)
         ! the rest of the records contain the eigenvectors
         DO jj = 1, dof
            READ (unit_nr) eigenvectors(1:dof, jj)
         END DO
      END IF
      ! broadcast to all compulational nodes. note that it is assumed
      ! that source is the ionode
      CALL para_env%bcast(dof)
      CALL para_env%bcast(eigenvalues)
      CALL para_env%bcast(eigenvectors)
      ! close file
      IF (para_env%is_source()) THEN
         CALL close_file(unit_number=unit_nr)
      END IF
   END SUBROUTINE read_vib_eigs_unformatted

END MODULE md_util

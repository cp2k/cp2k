!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Setup of regions with different temperature
!> \par History
!>   - Added support for Langevin regions (2014/01/08, LT)
!>   - Added print subroutine for langevin regions (2014/02/04, LT)
!>   - Changed print_thermal_regions to print_thermal_regions_temperature
!>     (2014/02/04, LT)
!> \author MI
! **************************************************************************************************
MODULE thermal_region_utils

   USE bibliography,                    ONLY: Kantorovich2008,&
                                              Kantorovich2008a,&
                                              cite_reference
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k,&
                                              cp_unit_to_cp2k
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE force_fields_util,               ONLY: get_generic_info
   USE fparser,                         ONLY: EvalErrType,&
                                              evalf,&
                                              finalizef,&
                                              initf,&
                                              parsef
   USE input_constants,                 ONLY: do_region_thermal,&
                                              langevin_ensemble,&
                                              npt_f_ensemble,&
                                              npt_i_ensemble,&
                                              npt_ia_ensemble,&
                                              nvt_ensemble
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE molecule_kind_list_types,        ONLY: molecule_kind_list_type
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_list_types,             ONLY: molecule_list_type
   USE molecule_types,                  ONLY: get_molecule,&
                                              molecule_type
   USE particle_list_types,             ONLY: particle_list_type
   USE physcon,                         ONLY: femtoseconds,&
                                              kelvin
   USE simpar_types,                    ONLY: simpar_type
   USE string_utilities,                ONLY: integer_to_string
   USE thermal_region_types,            ONLY: allocate_thermal_regions,&
                                              release_thermal_regions,&
                                              thermal_region_type,&
                                              thermal_regions_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE :: set_t_region_index
   PUBLIC :: create_thermal_regions, &
             print_thermal_regions_temperature, &
             print_thermal_regions_langevin

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermal_region_utils'

CONTAINS

! **************************************************************************************************
!> \brief create thermal_regions
!> \param thermal_regions ...
!> \param md_section ...
!> \param simpar ...
!> \param force_env ...
!> \par History
!>   - Added support for Langevin regions (2014/01/08, LT)
!> \author
! **************************************************************************************************
   SUBROUTINE create_thermal_regions(thermal_regions, md_section, simpar, force_env)
      TYPE(thermal_regions_type), POINTER                :: thermal_regions
      TYPE(section_vals_type), POINTER                   :: md_section
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(LEN=default_string_length)               :: my_region
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: tmpstringlist
      INTEGER                                            :: first_atom, i, ikind, il, imol, ipart, &
                                                            ireg, itimes, last_atom, nlist, &
                                                            nnpart, nregions, output_unit, region
      INTEGER, DIMENSION(:), POINTER                     :: tmplist
      LOGICAL                                            :: apply_thermostat, do_langevin, &
                                                            do_langevin_default, do_read_ngr, &
                                                            explicit, use_t, use_temp_tol, &
                                                            use_tfunc
      REAL(KIND=dp)                                      :: temp, temp_tol
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind, molecule_kind_set(:)
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(molecule_type), POINTER                       :: molecule
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(section_vals_type), POINTER                   :: region_sections, tfunc_section, &
                                                            thermal_region_section
      TYPE(thermal_region_type), POINTER                 :: t_region

      NULLIFY (logger, molecule, molecule_kind, molecule_kind_set, molecule_kinds, &
               region_sections, t_region, tfunc_section, &
               thermal_region_section, particles, subsys, tmplist)
      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)
      ALLOCATE (thermal_regions)
      CALL allocate_thermal_regions(thermal_regions)
      thermal_region_section => section_vals_get_subs_vals(md_section, "THERMAL_REGION")
      CALL section_vals_get(thermal_region_section, explicit=explicit)
      IF (explicit) THEN
         apply_thermostat = (simpar%ensemble == nvt_ensemble) .OR. &
                            (simpar%ensemble == npt_f_ensemble) .OR. &
                            (simpar%ensemble == npt_ia_ensemble) .OR. &
                            (simpar%ensemble == npt_i_ensemble)
         IF (apply_thermostat) THEN
            CALL section_vals_val_get(md_section, "THERMOSTAT%REGION", i_val=region)
            IF (region /= do_region_thermal) THEN
               CALL cp_warn(__LOCATION__, &
                            "An ensemble has been chosen where one or more thermostats are "// &
                            "used to control temperature, and one or more thermal regions "// &
                            "are also defined. To ensure consistency between thermostat "// &
                            "regions and thermal regions, set THERMOSTAT/REGION to THERMAL.")
            END IF
         END IF
         CALL section_vals_val_get(thermal_region_section, "FORCE_RESCALING", &
                                   l_val=thermal_regions%force_rescaling)
         region_sections => section_vals_get_subs_vals(thermal_region_section, &
                                                       "DEFINE_REGION")
         CALL section_vals_get(region_sections, n_repetition=nregions)
         IF (nregions > 0) THEN
            thermal_regions%nregions = nregions
            thermal_regions%section => thermal_region_section
            ALLOCATE (thermal_regions%thermal_region(nregions))
            simpar%do_thermal_region = .TRUE.
            CALL force_env_get(force_env, subsys=subsys)
            CALL cp_subsys_get(subsys, molecules=molecules, &
                               molecule_kinds=molecule_kinds, &
                               particles=particles)
            molecule_kind_set => molecule_kinds%els
            molecule_set => molecules%els

            DO ireg = 1, nregions
               NULLIFY (t_region)
               CALL integer_to_string(ireg, my_region)
               t_region => thermal_regions%thermal_region(ireg)
               t_region%region_index = ireg
               NULLIFY (t_region%part_index)

               ! TODO: parsing logic of thermal region should be aligned
               ! with thermostat region as in SUBROUTINE get_defined_region_info,
               ! src/motion/thermostat/thermostat_utils.F, so that
               ! MM_SUBSYS and QM_SUBSYS keywords work for thermal regions
               CALL section_vals_val_get(region_sections, "LIST", &
                                         i_rep_section=ireg, n_rep_val=nlist)
               IF (nlist > 0) THEN
                  DO il = 1, nlist
                     CALL section_vals_val_get(region_sections, "LIST", i_rep_section=ireg, &
                                               i_rep_val=il, i_vals=tmplist)
                     DO i = 1, SIZE(tmplist)
                        ipart = tmplist(i)
                        CALL set_t_region_index(particles, ipart, ireg)
                     END DO
                  END DO
               END IF
               CALL section_vals_val_get(region_sections, "MOLNAME", &
                                         i_rep_section=ireg, n_rep_val=nlist)
               IF (nlist > 0) THEN
                  DO il = 1, nlist
                     CALL section_vals_val_get(region_sections, "MOLNAME", i_rep_section=ireg, &
                                               i_rep_val=il, c_vals=tmpstringlist)
                     DO i = 1, SIZE(tmpstringlist)
                        DO ikind = 1, SIZE(molecule_kind_set)
                           molecule_kind => molecule_kind_set(ikind)
                           IF (molecule_kind%name == tmpstringlist(i)) THEN
                              DO imol = 1, SIZE(molecule_kind%molecule_list)
                                 molecule => molecule_set(molecule_kind%molecule_list(imol))
                                 CALL get_molecule(molecule, first_atom=first_atom, last_atom=last_atom)
                                 DO ipart = first_atom, last_atom
                                    CALL set_t_region_index(particles, ipart, ireg)
                                 END DO
                              END DO
                           END IF
                        END DO
                     END DO
                  END DO
               END IF

               ! Set up t_region%part_index after parsing input
               nnpart = 0
               t_region%npart = 0
               DO ipart = 1, particles%n_els
                  IF (particles%els(ipart)%t_region_index == ireg) nnpart = nnpart + 1
               END DO
               ALLOCATE (t_region%part_index(nnpart))
               DO ipart = 1, particles%n_els
                  IF (particles%els(ipart)%t_region_index == ireg) THEN
                     t_region%npart = t_region%npart + 1
                     t_region%part_index(t_region%npart) = ipart
                  END IF
               END DO

               ! Set up temperature function, if any
               tfunc_section => section_vals_get_subs_vals(region_sections, "TFUNC", i_rep_section=ireg)
               CALL section_vals_get(tfunc_section, explicit=use_tfunc)
               IF (use_tfunc) THEN
                  CALL get_generic_info(tfunc_section, "FUNCTION", &
                                        t_region%temperature_function, &
                                        t_region%temperature_function_parameters, &
                                        t_region%temperature_function_values, &
                                        input_variables=["S", "T"], i_rep_sec=1)
               END IF

               ! Set up t_region%temp_expected
               ! Precedence of temperature setting: the first positive value in
               ! 1. MD/THERMAL_REGION/DEFINE_REGION/TEMPERATURE
               ! 2. The temperature function evaluated at S = STEP_START_VAL and T = 0
               ! 3. MD/TEMPERATURE
               ! Otherwise it will be set to 0
               temp = 0.0_dp
               CALL section_vals_val_get(region_sections, "TEMPERATURE", &
                                         i_rep_section=ireg, explicit=use_t)
               IF (use_t) THEN
                  CALL section_vals_val_get(region_sections, "TEMPERATURE", &
                                            i_rep_section=ireg, r_val=temp)
               END IF
               IF (temp > 0.0_dp) THEN
                  t_region%temp_expected = temp
               ELSE
                  IF (use_tfunc) THEN
                     CALL initf(1)
                     CALL parsef(1, TRIM(t_region%temperature_function), &
                                 t_region%temperature_function_parameters)
                     CALL section_vals_val_get(md_section, "STEP_START_VAL", i_val=itimes)
                     t_region%temperature_function_values(1) = itimes*1.0_dp
                     t_region%temperature_function_values(2) = 0.0_dp
                     temp = cp_unit_to_cp2k(evalf(1, t_region%temperature_function_values), "K")
                     CPASSERT(EvalErrType <= 0)
                     CALL finalizef()
                  END IF
                  IF (temp > 0.0_dp) THEN
                     t_region%temp_expected = temp
                  ELSE
                     temp = simpar%temp_ext
                     IF (temp > 0.0_dp) THEN
                        t_region%temp_expected = temp
                     ELSE
                        CALL cp_warn(__LOCATION__, &
                                     "The initial temperature for the region "//TRIM(my_region)// &
                                     "cannot be assigned to a positive value in K from input. "// &
                                     "A value of zero will be used instead.")
                        t_region%temp_expected = 0.0_dp
                     END IF
                  END IF
               END IF

               ! Set up t_region%temp_tol
               ! Precedence of tolerance setting: the first positive value in
               ! 1. MD/THERMAL_REGION/DEFINE_REGION/TEMP_TOL
               ! 2. MD/TEMP_TOL
               ! Otherwise it will be set to 0 (i.e. not used)
               temp_tol = 0.0_dp
               CALL section_vals_val_get(region_sections, "TEMP_TOL", &
                                         i_rep_section=ireg, explicit=use_temp_tol)
               IF (use_temp_tol) THEN
                  CALL section_vals_val_get(region_sections, "TEMP_TOL", &
                                            i_rep_section=ireg, r_val=temp_tol)
               END IF
               IF (temp_tol > 0.0_dp) THEN
                  t_region%temp_tol = temp_tol
               ELSE
                  temp_tol = simpar%temp_tol
                  IF (temp_tol > 0.0_dp) THEN
                     t_region%temp_tol = temp_tol
                  ELSE
                     t_region%temp_tol = 0.0_dp
                  END IF
               END IF
               ! End of one thermal region in the loop
            END DO

            ! For Langevin ensemble, determine thermal_regions%do_langevin
            ! and t_region%noisy_gamma_region
            IF (simpar%ensemble == langevin_ensemble) THEN
               CALL cite_reference(Kantorovich2008)
               CALL cite_reference(Kantorovich2008a)
               CALL section_vals_val_get(thermal_region_section, "DO_LANGEVIN_DEFAULT", &
                                         l_val=do_langevin_default)
               ALLOCATE (thermal_regions%do_langevin(particles%n_els))
               thermal_regions%do_langevin = do_langevin_default
               DO ireg = 1, nregions
                  NULLIFY (t_region)
                  t_region => thermal_regions%thermal_region(ireg)
                  CALL integer_to_string(ireg, my_region)
                  CALL section_vals_val_get(region_sections, "DO_LANGEVIN", &
                                            i_rep_section=ireg, l_val=do_langevin)
                  DO ipart = 1, particles%n_els
                     IF (particles%els(ipart)%t_region_index == ireg) THEN
                        thermal_regions%do_langevin(ipart) = do_langevin
                     END IF
                  END DO
                  CALL section_vals_val_get(region_sections, "NOISY_GAMMA_REGION", i_rep_section=ireg, explicit=do_read_ngr)
                  IF (do_read_ngr) THEN
                     CALL section_vals_val_get(region_sections, "NOISY_GAMMA_REGION", i_rep_section=ireg, &
                                               r_val=t_region%noisy_gamma_region)
                     IF (.NOT. do_langevin) THEN
                        CALL cp_warn(__LOCATION__, &
                                     "You provided NOISY_GAMMA_REGION but atoms in thermal region "//TRIM(my_region)// &
                                     " will not undergo Langevin MD. "// &
                                     "NOISY_GAMMA_REGION will be ignored and its value discarded!")
                     END IF
                  ELSE
                     t_region%noisy_gamma_region = simpar%noisy_gamma
                  END IF
               END DO
            END IF
            ! End of Langevin treatment

            ! Output thermal region mapping to particles
            IF (output_unit > 0) THEN
               WRITE (output_unit, '(/,T2,A)') &
                  "THERMAL_REGION| Mapping of thermal region indices to particles"
               DO ipart = 1, particles%n_els, 16
                  WRITE (output_unit, '(T2,A,T17,16I4)') &
                     "THERMAL_REGION|", particles%els(ipart:MIN(ipart + 15, particles%n_els))%t_region_index
               END DO
            END IF
         ELSE
            CALL release_thermal_regions(thermal_regions)
            DEALLOCATE (thermal_regions)
            simpar%do_thermal_region = .FALSE.
         END IF
      ELSE
         CALL release_thermal_regions(thermal_regions)
         DEALLOCATE (thermal_regions)
         simpar%do_thermal_region = .FALSE.
      END IF

   END SUBROUTINE create_thermal_regions

! **************************************************************************************************
!> \brief set particles%els(ipart)%t_region_index to ireg
!> \param particles ...
!> \param ipart ...
!> \param ireg ...
!> \par History
!>   - Created (04/2026)
!> \author
! **************************************************************************************************
   SUBROUTINE set_t_region_index(particles, ipart, ireg)
      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER, INTENT(IN)                                :: ipart, ireg

      CHARACTER(LEN=default_string_length)               :: my_part, my_reg, my_region
      INTEGER                                            :: iregion

      CALL integer_to_string(ipart, my_part)
      IF ((ipart > 0) .AND. (ipart <= particles%n_els)) THEN
         CALL integer_to_string(ireg, my_reg)
         iregion = particles%els(ipart)%t_region_index
         IF ((iregion == 0) .OR. (iregion == ireg)) THEN
            ! No thermal region has been assigned, or
            ! the same one has already been assigned
            particles%els(ipart)%t_region_index = ireg
         ELSE
            CALL integer_to_string(iregion, my_region)
            CALL cp_abort(__LOCATION__, &
                          "Atom "//TRIM(my_part)//" has been "// &
                          "assigned to two different thermal "// &
                          "regions "//TRIM(my_region)//" and "// &
                          TRIM(my_reg)//" which is not allowed!")
         END IF
      ELSE
         IF (.NOT. ipart > 0) &
            CALL cp_abort(__LOCATION__, &
                          "Input atom index "//TRIM(my_part)//" is non-positive!")
         IF (.NOT. ipart <= particles%n_els) &
            CALL cp_abort(__LOCATION__, &
                          "Input atom index "//TRIM(my_part)//" is out of bounds!")
      END IF

   END SUBROUTINE set_t_region_index

! **************************************************************************************************
!> \brief print_thermal_regions_temperature
!> \param thermal_regions : thermal regions type contains information
!>                          about the regions
!> \param itimes          : iteration number of the time step
!> \param time            : simulation time of the time step
!> \param pos             : file position
!> \param act             : file action
!> \par History
!>   - added doxygen header and changed subroutine name from
!>     print_thermal_regions to print_thermal_regions_temperature
!>     (2014/02/04, LT)
!> \author
! **************************************************************************************************
   SUBROUTINE print_thermal_regions_temperature(thermal_regions, itimes, time, pos, act)
      TYPE(thermal_regions_type), POINTER                :: thermal_regions
      INTEGER, INTENT(IN)                                :: itimes
      REAL(KIND=dp), INTENT(IN)                          :: time
      CHARACTER(LEN=default_string_length)               :: pos, act

      CHARACTER(LEN=default_string_length)               :: fmd
      INTEGER                                            :: ireg, nregions, unit
      LOGICAL                                            :: new_file
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: temp, temp_tfunc
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (logger)
      logger => cp_get_default_logger()

      IF (ASSOCIATED(thermal_regions)) THEN
         print_key => section_vals_get_subs_vals(thermal_regions%section, "PRINT%TEMPERATURE")
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            unit = cp_print_key_unit_nr(logger, thermal_regions%section, "PRINT%TEMPERATURE", &
                                        extension=".tregion", file_position=pos, &
                                        file_action=act, is_new_file=new_file)
            IF (unit > 0) THEN
               IF (new_file) THEN
                  WRITE (unit, '(A)') "# Temperature per Region"
                  WRITE (unit, '("#",3X,A,2X,A,13X,A)') "Step Nr.", "Time[fs]", "Temp.[K] ...."
               END IF
               nregions = thermal_regions%nregions
               ALLOCATE (temp(0:nregions))
               ALLOCATE (temp_tfunc(1:nregions))
               temp = 0.0_dp
               temp(0) = thermal_regions%temp_reg0
               DO ireg = 1, nregions
                  ! Note the unit conversion here
                  temp(ireg) = thermal_regions%thermal_region(ireg)%temperature
                  temp_tfunc(ireg) = cp_unit_from_cp2k(thermal_regions%thermal_region(ireg)%temp_expected, "K")
               END DO
               fmd = "(I10,F20.3,"//TRIM(ADJUSTL(cp_to_string(2*nregions + 1)))//"F20.10)"
               fmd = TRIM(fmd)
               WRITE (UNIT=unit, FMT=fmd) itimes, time, temp(0), (temp_tfunc(ireg), temp(ireg), ireg=1, nregions)
               DEALLOCATE (temp)
               DEALLOCATE (temp_tfunc)
            END IF
            CALL cp_print_key_finished_output(unit, logger, thermal_regions%section, "PRINT%TEMPERATURE")
         END IF
      END IF
   END SUBROUTINE print_thermal_regions_temperature

! **************************************************************************************************
!> \brief print out information regarding to langevin regions defined in
!>        thermal_regions section
!> \param thermal_regions : thermal regions type containing the relevant
!>                          langevin regions information
!> \param simpar          : wrapper for simulation parameters
!> \param pos             : file position
!> \param act             : file action
!> \par History
!>   - created (2014/02/02, LT)
!> \author Lianheng Tong [LT] (tonglianheng@gmail.com)
! **************************************************************************************************
   SUBROUTINE print_thermal_regions_langevin(thermal_regions, simpar, pos, act)
      TYPE(thermal_regions_type), POINTER                :: thermal_regions
      TYPE(simpar_type), POINTER                         :: simpar
      CHARACTER(LEN=default_string_length)               :: pos, act

      INTEGER                                            :: ipart, ipart_reg, ireg, natoms, &
                                                            print_unit
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: region_id
      LOGICAL                                            :: new_file
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: noisy_gamma_region, temperature
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (logger)
      logger => cp_get_default_logger()

      IF (ASSOCIATED(thermal_regions)) THEN
         IF (ASSOCIATED(thermal_regions%do_langevin)) THEN
            print_key => section_vals_get_subs_vals(thermal_regions%section, &
                                                    "PRINT%LANGEVIN_REGIONS")
            IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), &
                      cp_p_file)) THEN
               print_unit = cp_print_key_unit_nr(logger, thermal_regions%section, &
                                                 "PRINT%LANGEVIN_REGIONS", &
                                                 extension=".lgv_regions", &
                                                 file_position=pos, file_action=act, &
                                                 is_new_file=new_file)
               IF (print_unit > 0) THEN
                  IF (new_file) THEN
                     WRITE (print_unit, '(A)') "# Atoms Undergoing Langevin MD"
                     WRITE (print_unit, '(A,3X,A,3X,A,3X,A,3X,A,3X,A)') &
                        "#", "Atom_ID", "Region_ID", "Langevin(L)/NVE(N)", "Expected_T[K]", "[NoisyGamma]"
                  END IF
                  natoms = SIZE(thermal_regions%do_langevin)
                  ALLOCATE (temperature(natoms))
                  ALLOCATE (region_id(natoms))
                  ALLOCATE (noisy_gamma_region(natoms))
                  temperature(:) = simpar%temp_ext
                  region_id(:) = 0
                  noisy_gamma_region(:) = simpar%noisy_gamma
                  DO ireg = 1, thermal_regions%nregions
                     DO ipart_reg = 1, thermal_regions%thermal_region(ireg)%npart
                        ipart = thermal_regions%thermal_region(ireg)%part_index(ipart_reg)
                        temperature(ipart) = thermal_regions%thermal_region(ireg)%temp_expected
                        region_id(ipart) = thermal_regions%thermal_region(ireg)%region_index
                        noisy_gamma_region(ipart) = thermal_regions%thermal_region(ireg)%noisy_gamma_region
                     END DO
                  END DO
                  DO ipart = 1, natoms
                     WRITE (print_unit, '(1X,I10,2X)', advance='no') ipart
                     WRITE (print_unit, '(I10,3X)', advance='no') region_id(ipart)
                     IF (thermal_regions%do_langevin(ipart)) THEN
                        WRITE (print_unit, '(A,3X)', advance='no') "L"
                        IF (noisy_gamma_region(ipart) > 0._dp) THEN
                           WRITE (print_unit, '(10X,F20.3,3X,ES9.3)') temperature(ipart)*kelvin, &
                              noisy_gamma_region(ipart)/femtoseconds
                        ELSE
                           WRITE (print_unit, '(10X,F20.3)') temperature(ipart)*kelvin
                        END IF
                     ELSE
                        WRITE (print_unit, '(A,3X)', advance='no') "N"
                        WRITE (print_unit, '(18X,A)') "--"
                     END IF
                  END DO
                  DEALLOCATE (region_id)
                  DEALLOCATE (temperature)
                  DEALLOCATE (noisy_gamma_region)
               END IF
               CALL cp_print_key_finished_output(print_unit, logger, thermal_regions%section, &
                                                 "PRINT%LANGEVIN_REGIONS")
            END IF
         END IF
      END IF
   END SUBROUTINE print_thermal_regions_langevin

END MODULE thermal_region_utils

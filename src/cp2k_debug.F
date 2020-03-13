!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Debug energy and derivatives w.r.t. finite differences
!> \note
!>      Use INTERPOLATION USE_GUESS, in order to perform force and energy
!>      calculations with the same density. This is not compulsory when iterating
!>      to selfconsistency, but essential in the non-selfconsistent case [08.2005,TdK].
!> \par History
!>      12.2004 created [tlaino]
!>      08.2005 consistent_energies option added, to allow FD calculations
!>              with the correct energies in the non-selfconsistent case, but
!>              keep in mind, that the QS energies and forces are then NOT
!>              consistent to each other [TdK].
!>      08.2005 In case the Harris functional is used, consistent_energies is
!>              et to .FALSE., otherwise the QS energies are spuriously used [TdK].
!>      01.2015 Remove Harris functional option
!>      - Revised (20.11.2013,MK)
!> \author Teodoro Laino
! **************************************************************************************************
MODULE cp2k_debug
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_result_methods,               ONLY: get_results,&
                                              test_for_result
   USE cp_result_types,                 ONLY: cp_result_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE force_env_methods,               ONLY: force_env_calc_energy_force,&
                                              force_env_calc_num_pressure
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type,&
                                              use_qs_force
   USE input_constants,                 ONLY: do_stress_analytical,&
                                              do_stress_diagonal_anal,&
                                              do_stress_diagonal_numer,&
                                              do_stress_numerical
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE particle_methods,                ONLY: write_fist_particle_coordinates,&
                                              write_qs_particle_coordinates
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE virial_types,                    ONLY: cp_virial,&
                                              virial_create,&
                                              virial_release,&
                                              virial_set,&
                                              virial_type,&
                                              zero_virial
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp2k_debug'
   REAL(KIND=dp), PRIVATE, PARAMETER :: maxerr = 5.0_dp

   PUBLIC :: cp2k_debug_energy_and_forces

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param force_env ...
! **************************************************************************************************
   SUBROUTINE cp2k_debug_energy_and_forces(force_env)

      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(LEN=3*default_string_length)             :: message
      CHARACTER(LEN=60)                                  :: line
      CHARACTER(LEN=default_string_length)               :: description
      INTEGER                                            :: i, ip, iw, j, k, np, stress_tensor
      LOGICAL                                            :: check_failed, debug_dipole, &
                                                            debug_forces, debug_polar, &
                                                            debug_stress_tensor, skip, &
                                                            stop_on_mismatch
      REAL(KIND=dp)                                      :: amplitude, dd, de, derr, difference, dx, &
                                                            eps_no_error_check, std_value, &
                                                            sum_of_differences
      REAL(KIND=dp), DIMENSION(2)                        :: numer_energy
      REAL(KIND=dp), DIMENSION(3)                        :: dipole_moment, dipole_numer, err, &
                                                            my_maxerr, poldir
      REAL(KIND=dp), DIMENSION(3, 2)                     :: dipn
      REAL(KIND=dp), DIMENSION(3, 3)                     :: polar_analytic, polar_numeric
      REAL(KIND=dp), DIMENSION(9)                        :: pvals
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: analyt_forces, numer_forces
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: root_section, subsys_section
      TYPE(virial_type), POINTER                         :: virial, virial_analytical, &
                                                            virial_numerical

      NULLIFY (analyt_forces, numer_forces, subsys, particles)

      root_section => force_env%root_section

      CALL force_env_get(force_env, para_env=para_env, subsys=subsys, cell=cell)
      subsys_section => section_vals_get_subs_vals(force_env%force_env_section, &
                                                   "SUBSYS")

      CALL section_vals_val_get(root_section, "DEBUG%DEBUG_STRESS_TENSOR", &
                                l_val=debug_stress_tensor)
      CALL section_vals_val_get(root_section, "DEBUG%DEBUG_FORCES", &
                                l_val=debug_forces)
      CALL section_vals_val_get(root_section, "DEBUG%DEBUG_DIPOLE", &
                                l_val=debug_dipole)
      CALL section_vals_val_get(root_section, "DEBUG%DEBUG_POLARIZABILITY", &
                                l_val=debug_polar)
      CALL section_vals_val_get(root_section, "DEBUG%DX", &
                                r_val=dx)
      CALL section_vals_val_get(root_section, "DEBUG%DE", &
                                r_val=de)
      dx = ABS(dx)
      CALL section_vals_val_get(root_section, "DEBUG%EPS_NO_ERROR_CHECK", &
                                r_val=eps_no_error_check)
      eps_no_error_check = MAX(eps_no_error_check, EPSILON(0.0_dp))
      CALL section_vals_val_get(root_section, "DEBUG%STOP_ON_MISMATCH", &
                                l_val=stop_on_mismatch)

      logger => cp_get_default_logger()
      iw = cp_print_key_unit_nr(logger, root_section, "DEBUG%PROGRAM_RUN_INFO", &
                                extension=".log")
      IF (debug_stress_tensor) THEN
         ! To debug stress tensor the stress tensor calculation must be
         ! first enabled..
         CALL section_vals_val_get(force_env%force_env_section, "STRESS_TENSOR", &
                                   i_val=stress_tensor)
         skip = .FALSE.
         SELECT CASE (stress_tensor)
         CASE (do_stress_analytical, do_stress_diagonal_anal)
            ! OK
         CASE (do_stress_numerical, do_stress_diagonal_numer)
            ! Nothing to check
            CALL cp_warn(__LOCATION__, "Numerical stress tensor was requested in "// &
                         "the FORCE_EVAL section. Nothing to debug!")
            skip = .TRUE.
         CASE DEFAULT
            CALL cp_warn(__LOCATION__, "Stress tensor calculation was not enabled in "// &
                         "the FORCE_EVAL section. Nothing to debug!")
            skip = .TRUE.
         END SELECT

         IF (.NOT. skip) THEN

            ! Compute the analytical stress tensor
            CALL cp_subsys_get(subsys, virial=virial)
            CALL virial_set(virial, pv_numer=.FALSE.)
            CALL force_env_calc_energy_force(force_env, &
                                             calc_force=.TRUE., &
                                             calc_stress_tensor=.TRUE.)

            ! Retrieve the analytical virial
            CALL virial_create(virial_analytical)
            CALL zero_virial(virial_analytical)
            CALL cp_virial(virial, virial_analytical)

            ! Debug stress tensor (numerical vs analytical)
            CALL virial_set(virial, pv_numer=.TRUE.)
            CALL force_env_calc_num_pressure(force_env, dx=dx)

            ! Retrieve the numerical virial
            CALL cp_subsys_get(subsys, virial=virial)
            CALL virial_create(virial_numerical)
            CALL zero_virial(virial_numerical)
            CALL cp_virial(virial, virial_numerical)

            ! Print results
            IF (iw > 0) THEN
               WRITE (UNIT=iw, FMT="((T2,A))") &
                  "DEBUG| Numerical pv_virial [a.u.]"
               WRITE (UNIT=iw, FMT="((T2,A,T21,3(1X,F19.12)))") &
                  ("DEBUG|", virial_numerical%pv_virial(i, 1:3), i=1, 3)
               WRITE (UNIT=iw, FMT="(/,(T2,A))") &
                  "DEBUG| Analytical pv_virial [a.u.]"
               WRITE (UNIT=iw, FMT="((T2,A,T21,3(1X,F19.12)))") &
                  ("DEBUG|", virial_analytical%pv_virial(i, 1:3), i=1, 3)
               WRITE (UNIT=iw, FMT="(/,(T2,A))") &
                  "DEBUG| Difference pv_virial [a.u.]"
               WRITE (UNIT=iw, FMT="((T2,A,T21,3(1X,F19.12)))") &
                  ("DEBUG|", virial_numerical%pv_virial(i, 1:3) - virial_analytical%pv_virial(i, 1:3), i=1, 3)
               WRITE (UNIT=iw, FMT="(T2,A,T61,F20.12)") &
                  "DEBUG| Sum of differences", &
                  SUM(ABS(virial_numerical%pv_virial(:, :) - virial_analytical%pv_virial(:, :)))
            END IF

            ! Check and abort on failure
            check_failed = .FALSE.
            IF (iw > 0) THEN
               WRITE (UNIT=iw, FMT="(/T2,A)") &
                  "DEBUG| Relative error pv_virial"
               WRITE (UNIT=iw, FMT="(T2,A,T61,ES20.8)") &
                  "DEBUG| Threshold value for error check [a.u.]", eps_no_error_check
            END IF
            DO i = 1, 3
               err(:) = 0.0_dp
               DO k = 1, 3
                  IF (ABS(virial_analytical%pv_virial(i, k)) >= eps_no_error_check) THEN
                     err(k) = 100.0_dp*(virial_numerical%pv_virial(i, k) - virial_analytical%pv_virial(i, k))/ &
                              virial_analytical%pv_virial(i, k)
                     WRITE (UNIT=line(20*(k - 1) + 1:20*k), FMT="(1X,F17.2,A2)") err(k), " %"
                  ELSE
                     WRITE (UNIT=line(20*(k - 1) + 1:20*k), FMT="(17X,A3)") "- %"
                  END IF
               END DO
               IF (iw > 0) THEN
                  WRITE (UNIT=iw, FMT="(T2,A,T21,A60)") &
                     "DEBUG|", line
               END IF
               IF (ANY(ABS(err(1:3)) > maxerr)) check_failed = .TRUE.
            END DO
            IF (iw > 0) THEN
               WRITE (UNIT=iw, FMT="(T2,A,T61,F18.2,A2)") &
                  "DEBUG| Maximum accepted error", maxerr, " %"
            END IF
            CALL virial_release(virial_analytical)
            CALL virial_release(virial_numerical)
            IF (check_failed) THEN
               message = "A mismatch between the analytical and the numerical "// &
                         "stress tensor has been detected. Check the implementation "// &
                         "of the stress tensor"
               IF (stop_on_mismatch) THEN
                  CPABORT(TRIM(message))
               ELSE
                  CPWARN(TRIM(message))
               END IF
            END IF
         END IF
      END IF

      IF (debug_forces) THEN
         ! Debug forces (numerical vs analytical)
         particles => subsys%particles%els
         SELECT CASE (force_env%in_use)
         CASE (use_qs_force)
            CALL get_qs_env(force_env%qs_env, qs_kind_set=qs_kind_set)
            CALL write_qs_particle_coordinates(particles, qs_kind_set, subsys_section, "DEBUG")
         CASE DEFAULT
            CALL write_fist_particle_coordinates(particles, subsys_section, charges=Null())
         END SELECT
         ! First evaluate energy and forces
         CALL force_env_calc_energy_force(force_env, &
                                          calc_force=.TRUE., &
                                          calc_stress_tensor=.FALSE.)
         ! Copy forces in array and start the numerical calculation
         IF (ASSOCIATED(analyt_forces)) DEALLOCATE (analyt_forces)
         np = subsys%particles%n_els
         ALLOCATE (analyt_forces(np, 3))
         DO ip = 1, np
            analyt_forces(ip, 1:3) = particles(ip)%f(1:3)
         END DO
         ! Loop on atoms and coordinates
         IF (ASSOCIATED(numer_forces)) DEALLOCATE (numer_forces)
         ALLOCATE (numer_forces(subsys%particles%n_els, 3))
         Atom: DO ip = 1, np
            Coord: DO k = 1, 3
               numer_energy = 0.0_dp
               std_value = particles(ip)%r(k)
               DO j = 1, 2
                  particles(ip)%r(k) = std_value - (-1.0_dp)**j*dx
                  SELECT CASE (force_env%in_use)
                  CASE (use_qs_force)
                     CALL get_qs_env(force_env%qs_env, qs_kind_set=qs_kind_set)
                     CALL write_qs_particle_coordinates(particles, qs_kind_set, subsys_section, "DEBUG")
                  CASE DEFAULT
                     CALL write_fist_particle_coordinates(particles, subsys_section, charges=Null())
                  END SELECT
                  ! Compute energy
                  CALL force_env_calc_energy_force(force_env, &
                                                   calc_force=.FALSE., &
                                                   calc_stress_tensor=.FALSE., &
                                                   consistent_energies=.TRUE.)
                  CALL force_env_get(force_env, potential_energy=numer_energy(j))
               END DO
               particles(ip)%r(k) = std_value
               numer_forces(ip, k) = -0.5_dp*(numer_energy(1) - numer_energy(2))/dx
               IF (iw > 0) THEN
                  WRITE (UNIT=iw, FMT="(/,T2,A,T17,A,F7.4,A,T34,A,F7.4,A,T52,A,T68,A)") &
                     "DEBUG| Atom", "E("//ACHAR(119 + k)//" +", dx, ")", &
                     "E("//ACHAR(119 + k)//" -", dx, ")", &
                     "f(numerical)", "f(analytical)"
                  WRITE (UNIT=iw, FMT="(T2,A,I5,4(1X,F16.8))") &
                     "DEBUG|", ip, numer_energy(1:2), numer_forces(ip, k), analyt_forces(ip, k)
               END IF
            END DO Coord
            ! Check analytical forces using the numerical forces as reference
            my_maxerr = maxerr
            err(1:3) = 0.0_dp
            DO k = 1, 3
               ! Calculate percentage but ignore very small force values
               IF (ABS(analyt_forces(ip, k)) >= eps_no_error_check) THEN
                  err(k) = 100.0_dp*(numer_forces(ip, k) - analyt_forces(ip, k))/analyt_forces(ip, k)
               END IF
               ! Increase error tolerance for small force values
               IF (ABS(analyt_forces(ip, k)) <= 0.0001_dp) my_maxerr(k) = 5.0_dp*my_maxerr(k)
            END DO
            IF (iw > 0) THEN
               WRITE (UNIT=iw, FMT="(/,T2,A)") &
                  "DEBUG| Atom  Coordinate   f(numerical)   f(analytical)   Difference   Error [%]"
               DO k = 1, 3
                  difference = analyt_forces(ip, k) - numer_forces(ip, k)
                  IF (ABS(analyt_forces(ip, k)) >= eps_no_error_check) THEN
                     WRITE (UNIT=iw, FMT="(T2,A,I5,T19,A1,T26,F14.8,T42,F14.8,T57,F12.8,T71,F10.2)") &
                        "DEBUG|", ip, ACHAR(119 + k), numer_forces(ip, k), analyt_forces(ip, k), difference, err(k)
                  ELSE
                     WRITE (UNIT=iw, FMT="(T2,A,I5,T19,A1,T26,F14.8,T42,F14.8,T57,F12.8,T80,A1)") &
                        "DEBUG|", ip, ACHAR(119 + k), numer_forces(ip, k), analyt_forces(ip, k), difference, "-"
                  END IF
               END DO
            END IF
            IF (ANY(ABS(err(1:3)) > my_maxerr(1:3))) THEN
               message = "A mismatch between analytical and numerical forces "// &
                         "has been detected. Check the implementation of the "// &
                         "analytical force calculation"
               IF (stop_on_mismatch) THEN
                  CPABORT(message)
               ELSE
                  CPWARN(message)
               END IF
            END IF
         END DO Atom
         ! Print summary
         IF (iw > 0) THEN
            WRITE (UNIT=iw, FMT="(/,(T2,A))") &
               "DEBUG|======================== BEGIN OF SUMMARY ===============================", &
               "DEBUG| Atom  Coordinate   f(numerical)   f(analytical)   Difference   Error [%]"
            sum_of_differences = 0.0_dp
            DO ip = 1, np
               err(1:3) = 0.0_dp
               DO k = 1, 3
                  difference = analyt_forces(ip, k) - numer_forces(ip, k)
                  IF (ABS(analyt_forces(ip, k)) >= eps_no_error_check) THEN
                     err(k) = 100_dp*(numer_forces(ip, k) - analyt_forces(ip, k))/analyt_forces(ip, k)
                     WRITE (UNIT=iw, FMT="(T2,A,I5,T19,A1,T26,F14.8,T42,F14.8,T57,F12.8,T71,F10.2)") &
                        "DEBUG|", ip, ACHAR(119 + k), numer_forces(ip, k), analyt_forces(ip, k), difference, err(k)
                  ELSE
                     WRITE (UNIT=iw, FMT="(T2,A,I5,T19,A1,T26,F14.8,T42,F14.8,T57,F12.8,T80,A1)") &
                        "DEBUG|", ip, ACHAR(119 + k), numer_forces(ip, k), analyt_forces(ip, k), difference, "-"
                  END IF
                  sum_of_differences = sum_of_differences + ABS(difference)
               END DO
            END DO
            WRITE (UNIT=iw, FMT="(T2,A,T57,F12.8)") &
               "DEBUG| Sum of differences:", sum_of_differences
            WRITE (UNIT=iw, FMT="(T2,A)") &
               "DEBUG|======================== END OF SUMMARY ================================="
         END IF
         ! Release work storage
         IF (ASSOCIATED(analyt_forces)) DEALLOCATE (analyt_forces)
         IF (ASSOCIATED(numer_forces)) DEALLOCATE (numer_forces)
      END IF

      IF (debug_dipole) THEN
         IF (force_env%in_use == use_qs_force) THEN
            CALL get_qs_env(force_env%qs_env, dft_control=dft_control)
            poldir = (/0.0_dp, 0.0_dp, 1.0_dp/)
            amplitude = 0.0_dp
            CALL set_efield(dft_control, amplitude, poldir)
            CALL force_env_calc_energy_force(force_env, calc_force=.TRUE.)
            CALL get_qs_env(force_env%qs_env, results=results)
            description = "[DIPOLE]"
            IF (test_for_result(results, description=description)) THEN
               CALL get_results(results, description=description, values=dipole_moment)
            ELSE
               CALL cp_warn(__LOCATION__, "Debug of dipole moments needs DFT/PRINT/MOMENTS section enabled")
               CPABORT("DEBUG")
            END IF
            amplitude = de
            DO k = 1, 3
               poldir = 0.0_dp
               poldir(k) = 1.0_dp
               DO j = 1, 2
                  poldir = -1.0_dp*poldir
                  CALL set_efield(dft_control, amplitude, poldir)
                  CALL force_env_calc_energy_force(force_env, calc_force=.FALSE.)
                  CALL force_env_get(force_env, potential_energy=numer_energy(j))
               END DO
               dipole_numer(k) = 0.5_dp*(numer_energy(1) - numer_energy(2))/de
            END DO
            IF (iw > 0) THEN
               WRITE (UNIT=iw, FMT="(/,(T2,A))") &
                  "DEBUG|========================= DIPOLE MOMENTS ================================", &
                  "DEBUG| Coordinate      D(numerical)    D(analytical)    Difference    Error [%]"
               DO k = 1, 3
                  dd = dipole_moment(k) - dipole_numer(k)
                  IF (ABS(dipole_moment(k)) > eps_no_error_check) THEN
                     derr = 100._dp*dd/dipole_moment(k)
                     WRITE (UNIT=iw, FMT="(T13,A1,T21,F16.8,T38,F16.8,T56,G12.3,T72,F9.3)") &
                        ACHAR(119 + k), dipole_numer(k), dipole_moment(k), dd, derr
                  ELSE
                     derr = 0.0_dp
                     WRITE (UNIT=iw, FMT="(T13,A1,T21,F16.8,T38,F16.8,T56,G12.3)") &
                        ACHAR(119 + k), dipole_numer(k), dipole_moment(k), dd
                  END IF
                  err(k) = derr
               END DO
               WRITE (UNIT=iw, FMT="((T2,A))") &
                  "DEBUG|========================================================================="
               WRITE (UNIT=iw, FMT="(T2,A,T61,E20.12)") 'DIPOLE : CheckSum  =', SUM(dipole_moment)
               IF (ANY(ABS(err(1:3)) > 0.1_dp)) THEN
                  message = "A mismatch between analytical and numerical dipoles "// &
                            "has been detected. Check the implementation of the "// &
                            "analytical dipole calculation"
                  IF (stop_on_mismatch) THEN
                     CPABORT(message)
                  ELSE
                     CPWARN(message)
                  END IF
               END IF
            END IF

         ELSE
            CALL cp_warn(__LOCATION__, "Debug of dipole moments only for Quickstep code available")
         END IF
      END IF

      IF (debug_polar) THEN
         IF (force_env%in_use == use_qs_force) THEN
            CALL get_qs_env(force_env%qs_env, dft_control=dft_control)
            poldir = (/0.0_dp, 0.0_dp, 1.0_dp/)
            amplitude = 0.0_dp
            CALL set_efield(dft_control, amplitude, poldir)
            CALL force_env_calc_energy_force(force_env, calc_force=.FALSE.)
            CALL get_qs_env(force_env%qs_env, results=results)
            description = "[POLAR]"
            IF (test_for_result(results, description=description)) THEN
               CALL get_results(results, description=description, values=pvals)
               polar_analytic(1:3, 1:3) = RESHAPE(pvals(1:9), (/3, 3/))
            ELSE
               CALL cp_warn(__LOCATION__, "Debug of polarizabilities needs PROPERTIES/LINRES/POLAR section enabled")
               CPABORT("DEBUG")
            END IF
            description = "[DIPOLE]"
            IF (.NOT. test_for_result(results, description=description)) THEN
               CALL cp_warn(__LOCATION__, "Debug of polarizabilities need DFT/PRINT/MOMENTS section enabled")
               CPABORT("DEBUG")
            END IF
            amplitude = de
            DO k = 1, 3
               poldir = 0.0_dp
               poldir(k) = 1.0_dp
               DO j = 1, 2
                  poldir = -1.0_dp*poldir
                  CALL set_efield(dft_control, amplitude, poldir)
                  CALL force_env_calc_energy_force(force_env, calc_force=.FALSE., linres=.TRUE.)
                  CALL get_results(results, description=description, values=dipn(1:3, j))
               END DO
               polar_numeric(1:3, k) = 0.5_dp*(dipn(1:3, 2) - dipn(1:3, 1))/de
            END DO
            IF (iw > 0) THEN
               WRITE (UNIT=iw, FMT="(/,(T2,A))") &
                  "DEBUG|========================= POLARIZABILITY ================================", &
                  "DEBUG| Coordinates     P(numerical)    P(analytical)    Difference    Error [%]"
               DO k = 1, 3
                  DO j = 1, 3
                     dd = polar_analytic(k, j) - polar_numeric(k, j)
                     IF (ABS(polar_analytic(k, j)) > eps_no_error_check) THEN
                        derr = 100._dp*dd/polar_analytic(k, j)
                        WRITE (UNIT=iw, FMT="(T12,A1,A1,T21,F16.8,T38,F16.8,T56,G12.3,T72,F9.3)") &
                           ACHAR(119 + k), ACHAR(119 + j), polar_numeric(k, j), polar_analytic(k, j), dd, derr
                     ELSE
                        WRITE (UNIT=iw, FMT="(T12,A1,A1,T21,F16.8,T38,F16.8,T56,G12.3)") &
                           ACHAR(119 + k), ACHAR(119 + j), polar_numeric(k, j), polar_analytic(k, j), dd
                     END IF
                  END DO
               END DO
               WRITE (UNIT=iw, FMT="((T2,A))") &
                  "DEBUG|========================================================================="
               WRITE (UNIT=iw, FMT="(T2,A,T61,E20.12)") ' POLAR : CheckSum  =', SUM(polar_analytic)
            END IF
         ELSE
            CALL cp_warn(__LOCATION__, "Debug of polarizabilities only for Quickstep code available")
         END IF
      END IF

      CALL cp_print_key_finished_output(iw, logger, root_section, "DEBUG%PROGRAM_RUN_INFO")

   END SUBROUTINE cp2k_debug_energy_and_forces

! **************************************************************************************************
!> \brief ...
!> \param dft_control ...
!> \param amplitude ...
!> \param poldir ...
! **************************************************************************************************
   SUBROUTINE set_efield(dft_control, amplitude, poldir)
      TYPE(dft_control_type), POINTER                    :: dft_control
      REAL(KIND=dp), INTENT(IN)                          :: amplitude
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: poldir

      IF (dft_control%apply_efield) THEN
         dft_control%efield_fields(1)%efield%strength = amplitude
         dft_control%efield_fields(1)%efield%polarisation(1:3) = poldir(1:3)
      ELSEIF (dft_control%apply_period_efield) THEN
         dft_control%period_efield%strength = amplitude
         dft_control%period_efield%polarisation(1:3) = poldir(1:3)
      ELSE
         CPABORT("No EFIELD definition available")
      END IF

   END SUBROUTINE set_efield

END MODULE cp2k_debug

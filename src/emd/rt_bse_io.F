!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Input/output from the propagation via RT-BSE method.
!> \author Stepan Marek (08.24)
! **************************************************************************************************

MODULE rt_bse_io
   USE cp_fm_types, ONLY: cp_fm_type, &
                          cp_fm_read_unformatted, &
                          cp_fm_write_unformatted, &
                          cp_fm_write_formatted
   USE cp_cfm_types, ONLY: cp_cfm_type, &
                           cp_fm_to_cfm, &
                           cp_cfm_to_fm
   USE kinds, ONLY: dp, &
                    default_path_length
   USE cp_fm_basic_linalg, ONLY: cp_fm_trace, &
                                 cp_fm_transpose
   USE cp_control_types, ONLY: dft_control_type
   USE cp_log_handling, ONLY: cp_logger_type, &
                              cp_get_default_logger, &
                              cp_to_string
   USE cp_output_handling, ONLY: cp_print_key_unit_nr, &
                                 cp_print_key_finished_output, &
                                 cp_print_key_generate_filename, &
                                 low_print_level, &
                                 medium_print_level
   USE force_env_types, ONLY: force_env_type
   USE input_section_types, ONLY: section_vals_get_subs_vals, &
                                  section_vals_type, &
                                  section_vals_val_get
   USE rt_bse_types, ONLY: rtbse_env_type, &
                           multiply_cfm_fm, &
                           multiply_fm_cfm
   USE cp_files, ONLY: open_file, &
                       file_exists, &
                       close_file
   USE input_constants, ONLY: do_exact, &
                              do_bch, &
                              evGW0, &
                              rtp_bse_ham_gw, &
                              rtp_bse_ham_ks, &
                              use_rt_restart, &
                              restart_guess
   USE qs_environment_types, ONLY: get_qs_env
   USE scf_control_types, ONLY: scf_control_type
   USE physcon, ONLY: evolt, femtoseconds
   USE rt_propagation_output, ONLY: print_moments, &
                                    print_rt_file, &
                                    rt_file_comp_real

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_bse_io"

   ! RESTART.trace on-disk format version - bump whenever the header or record layout changes
   INTEGER, PARAMETER, PRIVATE :: restart_trace_version = 1

   #:include "rt_bse_macros.fypp"

   PUBLIC :: output_moments, &
             output_field, &
             read_field, &
             output_mos_contravariant, &
             output_mos_covariant, &
             output_restart, &
             read_restart, &
             output_restart_linearized, &
             check_restart_files_before_scf, &
             read_restart_info, &
             read_restart_trace, &
             read_restart_density, &
             read_restart_C, &
             check_restart_eps_consistency, &
             print_etrs_info_header, &
             print_etrs_info, &
             print_timestep_info, &
             print_rtbse_header_info

CONTAINS

! **************************************************************************************************
!> \brief Writes the header and basic info to the standard output
!> \param rtbse_env Entry point - rtbse environment
! **************************************************************************************************
   SUBROUTINE print_rtbse_header_info(rtbse_env)
      TYPE(rtbse_env_type)                               :: rtbse_env
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => cp_get_default_logger()

      IF (rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, *) ''
         WRITE (rtbse_env%unit_nr, '(A)') ' /-----------------------------------------------'// &
            '------------------------------\'
         WRITE (rtbse_env%unit_nr, '(A)') ' |                                               '// &
            '                              |'
         WRITE (rtbse_env%unit_nr, '(A)') ' |                    Real Time Bethe-Salpeter Propagation'// &
            '                     |'
         WRITE (rtbse_env%unit_nr, '(A)') ' |                                               '// &
            '                              |'
         WRITE (rtbse_env%unit_nr, '(A)') ' \-----------------------------------------------'// &
            '------------------------------/'
         WRITE (rtbse_env%unit_nr, *) ''

         ! Methods used
         WRITE (rtbse_env%unit_nr, '(A19)', advance="no") ' Exponential method'
         SELECT CASE (rtbse_env%mat_exp_method)
         CASE (do_bch)
            WRITE (rtbse_env%unit_nr, '(A61)') 'BCH'
         CASE (do_exact)
            WRITE (rtbse_env%unit_nr, '(A61)') 'EXACT'
         END SELECT

         WRITE (rtbse_env%unit_nr, '(A22)', advance="no") ' Reference Hamiltonian'
         SELECT CASE (rtbse_env%ham_reference_type)
         CASE (rtp_bse_ham_gw)
            IF (rtbse_env%bs_env%gw_flavour == evGW0) THEN
               WRITE (rtbse_env%unit_nr, '(A58)') 'evGW0'
            ELSE
               WRITE (rtbse_env%unit_nr, '(A58)') 'G0W0'
            END IF
         CASE (rtp_bse_ham_ks)
            WRITE (rtbse_env%unit_nr, '(A58)') 'Kohn-Sham'
         END SELECT

         WRITE (rtbse_env%unit_nr, '(A18,L62)') ' Apply delta pulse', &
            rtbse_env%dft_control%rtp_control%apply_delta_pulse

         WRITE (rtbse_env%unit_nr, '(A)') ''
      END IF

   END SUBROUTINE print_rtbse_header_info

! **************************************************************************************************
!> \brief Writes the update after single etrs iteration - only for log level > medium
!> \param rtbse_env Entry point - rtbse environment
! **************************************************************************************************
   SUBROUTINE print_etrs_info(rtbse_env, step, metric)
      TYPE(rtbse_env_type)                               :: rtbse_env
      INTEGER                                            :: step
      REAL(kind=dp)                                      :: metric
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => cp_get_default_logger()

      IF (logger%iter_info%print_level > medium_print_level .AND. rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A7,I5, E20.8E3)') ' RTBSE|', step, metric
      END IF

   END SUBROUTINE print_etrs_info
! **************************************************************************************************
!> \brief Writes the header for the etrs iteration updates - only for log level > medium
!> \param rtbse_env Entry point - rtbse environment
! **************************************************************************************************
   SUBROUTINE print_etrs_info_header(rtbse_env)
      TYPE(rtbse_env_type)                               :: rtbse_env
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => cp_get_default_logger()

      IF (logger%iter_info%print_level > medium_print_level .AND. rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A13, A20)') ' RTBSE| Iter.', 'Convergence'
      END IF

   END SUBROUTINE print_etrs_info_header
! **************************************************************************************************
!> \brief Writes the summary line of a completed propagation timestep
!> \param rtbse_env Entry point - rtbse environment
!> \param step Index of the completed timestep
!> \param electron_num_re Real part of the electron number, one entry per spin channel
!> \param convergence Optional convergence metric reached by the ETRS iteration
!> \param etrs_num Optional number of ETRS iterations used in this timestep
!> \param step_walltime Optional wall time of this timestep, in seconds
! **************************************************************************************************
   SUBROUTINE print_timestep_info(rtbse_env, step, electron_num_re, convergence, etrs_num, &
                                  step_walltime)
      TYPE(rtbse_env_type)                               :: rtbse_env
      INTEGER                                            :: step
      REAL(kind=dp), DIMENSION(:), INTENT(IN)            :: electron_num_re
      REAL(kind=dp), OPTIONAL                            :: convergence
      INTEGER, OPTIONAL                                  :: etrs_num
      REAL(kind=dp), OPTIONAL                            :: step_walltime
      TYPE(cp_logger_type), POINTER                      :: logger
      LOGICAL                                          :: flag_lrrtbse
      INTEGER                                            :: nch

      logger => cp_get_default_logger()
      ! one electron-number column per spin channel (1 = closed shell, 2 = open shell alpha/beta)
      nch = SIZE(electron_num_re)

      IF (.NOT. PRESENT(convergence) .OR. .NOT. PRESENT(etrs_num)) THEN
         flag_lrrtbse = .TRUE.
      ELSE
         flag_lrrtbse = .FALSE.
      END IF

      IF (logger%iter_info%print_level > low_print_level .AND. rtbse_env%unit_nr > 0) THEN
         IF (flag_lrrtbse) THEN
            IF (step == 0) THEN
               IF (PRESENT(step_walltime)) THEN
                  WRITE (rtbse_env%unit_nr, '(A45,T70,F11.3)') &
                     " RTBSE| Estimated runtime for propagation [s]", &
                     step_walltime*REAL(rtbse_env%sim_nsteps, dp)
                  WRITE (rtbse_env%unit_nr, '(A)') &
                     " RTBSE|"
                  IF (nch == 1) THEN
                     WRITE (rtbse_env%unit_nr, '(A23,T27,A13,T66,A15)') &
                        " RTBSE| Simulation step", "Step time [s]", "Electron number"
                  ELSE
                     ! T67/T88 (not T66/T86): each α/β is 2 bytes but 1 display col, so the byte
                     ! anchor is +1 per unicode char to right-align ')' under the value's last digit.
                     ! T anchors carry +1 byte per α/β before the ')' (α col → +1, β col → +2),
                     ! since each is 2 bytes but 1 display col; right-aligns ')' on the last digit.
                     WRITE (rtbse_env%unit_nr, '(A23,T27,A13,T47,A15,T68,A15)') &
                        " RTBSE| Simulation step", "Step time [s]", "El. number (α)", "El. number (β)"
                  END IF
               ELSE
                  IF (nch == 1) THEN
                     WRITE (rtbse_env%unit_nr, '(A23,T66,A15)') " RTBSE| Simulation step", "Electron number"
                  ELSE
                     WRITE (rtbse_env%unit_nr, '(A23,T47,A15,T68,A15)') &
                        " RTBSE| Simulation step", "El. number (α)", "El. number (β)"
                  END IF
               END IF
            END IF
            IF (PRESENT(step_walltime)) THEN
               IF (nch == 1) THEN
                  WRITE (rtbse_env%unit_nr, '(A7,I16,T30,F10.3,T69,E12.3E3)') &
                     ' RTBSE|', step, step_walltime, electron_num_re(1)
               ELSE
                  WRITE (rtbse_env%unit_nr, '(A7,I16,T30,F10.3,T49,E12.3E3,T69,E12.3E3)') &
                     ' RTBSE|', step, step_walltime, electron_num_re(1), electron_num_re(2)
               END IF
            ELSE
               IF (nch == 1) THEN
                  WRITE (rtbse_env%unit_nr, '(A7,I16,T61,E20.8E3)') ' RTBSE|', step, electron_num_re(1)
               ELSE
                  WRITE (rtbse_env%unit_nr, '(A7,I16,T49,E12.3E3,T69,E12.3E3)') &
                     ' RTBSE|', step, electron_num_re(1), electron_num_re(2)
               END IF
            END IF
         ELSE
            WRITE (rtbse_env%unit_nr, '(A23,A20,A20,A17)') " RTBSE| Simulation step", "Convergence", &
               "Electron number", "ETRS Iterations"
            WRITE (rtbse_env%unit_nr, '(A7,I16,E20.8E3,E20.8E3,I17)') ' RTBSE|', step, convergence, &
               electron_num_re(1), etrs_num
         END IF
      END IF

   END SUBROUTINE print_timestep_info

! **************************************************************************************************
!> \brief Outputs the matrix in MO basis for matrix coefficients corresponding to contravariant
!>        operator, i.e. density matrix
!> \param rtbse_env Entry point - gwbse environment
!> \param rho Density matrix in AO basis
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE output_mos_contravariant(rtbse_env, rho, print_key_section)
      TYPE(rtbse_env_type)                               :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho
      TYPE(section_vals_type), POINTER                   :: print_key_section
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: j, rho_unit_re, rho_unit_im
      CHARACTER(len=14), DIMENSION(4)                    :: file_labels

      file_labels(1) = "_SPIN_A_RE.dat"
      file_labels(2) = "_SPIN_A_IM.dat"
      file_labels(3) = "_SPIN_B_RE.dat"
      file_labels(4) = "_SPIN_B_IM.dat"
      logger => cp_get_default_logger()

      ! In the linearized RT-BSE active-MO path, rho is already in the MO basis
      ! restricted to the active window (sized mo_active x mo_active).  Dump it
      ! directly without the AO-side C^T S * rho * S C transformation.
      IF (ASSOCIATED(rtbse_env%real_workspace_mo)) THEN
         DO j = 1, rtbse_env%n_spin
            rho_unit_re = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(2*j - 1))
            rho_unit_im = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(2*j))
            CALL cp_cfm_to_fm(rho(j), rtbse_env%real_workspace_mo(1), rtbse_env%real_workspace_mo(2))
            CALL cp_fm_write_formatted(rtbse_env%real_workspace_mo(1), rho_unit_re)
            CALL cp_fm_write_formatted(rtbse_env%real_workspace_mo(2), rho_unit_im)
            CALL cp_print_key_finished_output(rho_unit_re, logger, print_key_section)
            CALL cp_print_key_finished_output(rho_unit_im, logger, print_key_section)
         END DO
         RETURN
      END IF
      ! Start by multiplying the current density by MOS
      DO j = 1, rtbse_env%n_spin
         rho_unit_re = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(2*j - 1))
         rho_unit_im = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(2*j))
         ! Transform the density matrix into molecular orbitals basis and print it out
         ! S * rho
         CALL multiply_fm_cfm("N", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                              1.0_dp, rtbse_env%S_fm, rho(j), &
                              0.0_dp, rtbse_env%rho_workspace(1))
         ! C^T * S * rho
         CALL multiply_fm_cfm("T", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                              1.0_dp, rtbse_env%bs_env%fm_mo_coeff_Gamma(j), rtbse_env%rho_workspace(1), &
                              0.0_dp, rtbse_env%rho_workspace(2))
         ! C^T * S * rho * S
         CALL multiply_cfm_fm("N", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                              1.0_dp, rtbse_env%rho_workspace(2), rtbse_env%S_fm, &
                              0.0_dp, rtbse_env%rho_workspace(1))
         ! C^T * S * rho * S * C
         CALL multiply_cfm_fm("N", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                              1.0_dp, rtbse_env%rho_workspace(1), rtbse_env%bs_env%fm_mo_coeff_Gamma(j), &
                              0.0_dp, rtbse_env%rho_workspace(2))
         ! Print real and imaginary parts separately
         CALL cp_cfm_to_fm(rtbse_env%rho_workspace(2), &
                           rtbse_env%real_workspace(1), rtbse_env%real_workspace(2))
         CALL cp_fm_write_formatted(rtbse_env%real_workspace(1), rho_unit_re)
         CALL cp_fm_write_formatted(rtbse_env%real_workspace(2), rho_unit_im)
         CALL cp_print_key_finished_output(rho_unit_re, logger, print_key_section)
         CALL cp_print_key_finished_output(rho_unit_im, logger, print_key_section)
      END DO
   END SUBROUTINE output_mos_contravariant
! **************************************************************************************************
!> \brief Outputs the matrix in MO basis for matrix components corresponding to covariant representation,
!>        i.e. the Hamiltonian matrix
!> \param rtbse_env Entry point - gwbse environment
!> \param cohsex cohsex matrix in AO basis, covariant representation
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE output_mos_covariant(rtbse_env, ham, print_key_section)
      TYPE(rtbse_env_type)                               :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: ham
      TYPE(section_vals_type), POINTER                   :: print_key_section
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: j, rho_unit_re, rho_unit_im
      CHARACTER(len=21), DIMENSION(4)                    :: file_labels

      file_labels(1) = "_SPIN_A_RE.dat"
      file_labels(2) = "_SPIN_A_IM.dat"
      file_labels(3) = "_SPIN_B_RE.dat"
      file_labels(4) = "_SPIN_B_IM.dat"
      logger => cp_get_default_logger()
      DO j = 1, rtbse_env%n_spin
         rho_unit_re = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(2*j - 1))
         rho_unit_im = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(2*j))
         ! C^T * cohsex
         CALL multiply_fm_cfm("T", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                              1.0_dp, rtbse_env%bs_env%fm_mo_coeff_Gamma(j), ham(j), &
                              0.0_dp, rtbse_env%rho_workspace(1))
         ! C^T * cohsex * C
         CALL multiply_cfm_fm("N", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                              1.0_dp, rtbse_env%rho_workspace(1), rtbse_env%bs_env%fm_mo_coeff_Gamma(j), &
                              0.0_dp, rtbse_env%rho_workspace(2))
         ! Print real and imaginary parts separately
         CALL cp_cfm_to_fm(rtbse_env%rho_workspace(2), &
                           rtbse_env%real_workspace(1), rtbse_env%real_workspace(2))
         CALL cp_fm_write_formatted(rtbse_env%real_workspace(1), rho_unit_re)
         CALL cp_fm_write_formatted(rtbse_env%real_workspace(2), rho_unit_im)
         CALL cp_print_key_finished_output(rho_unit_re, logger, print_key_section)
         CALL cp_print_key_finished_output(rho_unit_im, logger, print_key_section)
      END DO
   END SUBROUTINE output_mos_covariant
! **************************************************************************************************
!> \brief Prints the current field components into a file provided by input
!> \param rtbse_env Entry point - gwbse environment
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE output_field(rtbse_env, append_opt)
      TYPE(rtbse_env_type)                               :: rtbse_env
      LOGICAL, OPTIONAL                                  :: append_opt
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: field_unit, n, i
      LOGICAL                                            :: append

      ! Figure out whether we are appending or not, true by default
      append = .TRUE.
      IF (PRESENT(append_opt)) append = .FALSE.

      ! First, write the current field to memory
      ! Need the absolute index
      n = rtbse_env%sim_step - rtbse_env%sim_start_orig + 1
      DO i = 1, 3
         rtbse_env%field_trace(i, n) = CMPLX(rtbse_env%field(i), 0.0, kind=dp)
      END DO
      rtbse_env%time_trace(n) = rtbse_env%sim_time

      ! Now, continue to file output
      ! Get logger
      logger => cp_get_default_logger()
      ! Get file descriptor
      field_unit = cp_print_key_unit_nr(logger, rtbse_env%field_section, extension=".dat")
      IF (append) THEN
         CALL print_rt_file(field_unit, xvals=rtbse_env%time_trace(n:n), &
                            yvals=rtbse_env%field_trace(:, n:n), &
                            xscale_opt=femtoseconds, comp_opt=rt_file_comp_real)
      ELSE
         CALL print_rt_file(field_unit, [ &
                            "#          Time [fs]", &
                            "     field x [at.u.]", &
                            "     field y [at.u.]", &
                            "     field z [at.u.]"], &
                            rtbse_env%time_trace(n:n), rtbse_env%field_trace(:, n:n), &
                            xscale_opt=femtoseconds, comp_opt=rt_file_comp_real)
      END IF
      CALL cp_print_key_finished_output(field_unit, logger, rtbse_env%field_section)

   END SUBROUTINE output_field
! **************************************************************************************************
!> \brief Reads the field from the files provided by input - useful for the continuation run
!> \param rtbse_env Entry point - gwbse environment
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE read_field(rtbse_env)
      TYPE(rtbse_env_type)                               :: rtbse_env
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(len=default_path_length)                 :: save_name
      INTEGER                                            :: k, n, field_unit
      REAL(kind=dp), DIMENSION(3)                        :: real_field

      ! Get logger
      logger => cp_get_default_logger()
      ! Get file name
      save_name = cp_print_key_generate_filename(logger, rtbse_env%field_section, extension=".dat", my_local=.FALSE.)
      IF (file_exists(save_name)) THEN
         CALL open_file(save_name, file_status="OLD", file_form="FORMATTED", file_action="READ", &
                        unit_number=field_unit)
         ! Skip the first line - it contains headers
         READ (field_unit, '()')
         DO k = rtbse_env%sim_start_orig, rtbse_env%sim_start
            n = k - rtbse_env%sim_start_orig + 1
            READ (field_unit, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') rtbse_env%time_trace(n), &
               real_field(1), real_field(2), real_field(3)
            rtbse_env%field_trace(:, n) = CMPLX(real_field(:), 0.0, kind=dp)
            ! Set the time units back to atomic units
            rtbse_env%time_trace(n) = rtbse_env%time_trace(n)/femtoseconds
         END DO
         CALL close_file(field_unit)
      ELSE IF (.NOT. rtbse_env%dft_control%rtp_control%apply_delta_pulse .AND. &
               rtbse_env%dft_control%rtp_control%initial_wfn == use_rt_restart) THEN
         CPWARN("Restart without RT field file - unknown field trace set to zero.")
      END IF
   END SUBROUTINE read_field

! **************************************************************************************************
!> \brief Outputs the expectation value of moments from a given density matrix
!> \note  Moments matrix is provided by the rtbse_env, uses rho_workspace(1:3)
!> \param rtbse_env Entry point - gwbse environment
!> \param rho Density matrix in AO basis
!> \param rtp_section RTP section of the input parameters, where moments destination may be present
! **************************************************************************************************
   SUBROUTINE output_moments(rtbse_env, rho)
      TYPE(rtbse_env_type)                               :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho
      INTEGER                                            :: i, j, n
      REAL(kind=dp), DIMENSION(3)                        :: moments_re
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: ws

      n = rtbse_env%sim_step - rtbse_env%sim_start_orig + 1

      ! In linearized RT-BSE rho and moments are MO-active sized; otherwise AO sized.
      IF (ASSOCIATED(rtbse_env%real_workspace_mo)) THEN
         ws => rtbse_env%real_workspace_mo
      ELSE
         ws => rtbse_env%real_workspace
      END IF

      DO j = 1, rtbse_env%n_spin
         ! Need to transpose due to the definition of trace function
         CALL cp_cfm_to_fm(msource=rho(j), mtargetr=ws(2))
         DO i = 1, 3
            ! Moments should be symmetric, test without transopose?
            CALL cp_fm_transpose(rtbse_env%moments(i, j), ws(1))
            CALL cp_fm_trace(ws(1), ws(2), moments_re(i))
            ! Scale by spin degeneracy and electron charge
            moments_re(i) = -moments_re(i)*rtbse_env%spin_degeneracy
            rtbse_env%moments_trace(j, i, n) = CMPLX(moments_re(i), 0.0, kind=dp)
         END DO
         ! Same for imaginary part
         CALL cp_cfm_to_fm(msource=rho(j), mtargeti=ws(2))
         DO i = 1, 3
            CALL cp_fm_transpose(rtbse_env%moments(i, j), ws(1))
            CALL cp_fm_trace(ws(1), ws(2), moments_re(i))
            ! Scale by spin degeneracy and electron charge
            moments_re(i) = -moments_re(i)*rtbse_env%spin_degeneracy
            rtbse_env%moments_trace(j, i, n) = rtbse_env%moments_trace(j, i, n) + CMPLX(0.0, moments_re(i), kind=dp)
         END DO
      END DO
      ! Output to the file
      CALL print_moments(rtbse_env%moments_section, rtbse_env%unit_nr, rtbse_env%moments_trace(:, :, n), &
                         rtbse_env%sim_time, .TRUE., append_opt=(rtbse_env%sim_step /= rtbse_env%sim_start_orig))
   END SUBROUTINE output_moments
! **************************************************************************************************
!> \brief Outputs the restart info (last finished iteration step) + restard density matrix
!> \param restart_section Print key section for the restart files
!> \param rho Density matrix in AO basis
!> \param time_index Time index to be written into the info file
! **************************************************************************************************
   SUBROUTINE output_restart(rtbse_env, rho, time_index)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho
      INTEGER                                            :: time_index
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: workspace
      CHARACTER(len=17), DIMENSION(4)                    :: file_labels
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: rho_unit_nr, i

      ! Default labels distinguishing up to two spin species and real/imaginary parts
      file_labels(1) = "_SPIN_A_RE.matrix"
      file_labels(2) = "_SPIN_A_IM.matrix"
      file_labels(3) = "_SPIN_B_RE.matrix"
      file_labels(4) = "_SPIN_B_IM.matrix"

      logger => cp_get_default_logger()

      workspace => rtbse_env%real_workspace

      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_to_fm(rho(i), workspace(1), workspace(2))
         ! Real part
         rho_unit_nr = cp_print_key_unit_nr(logger, rtbse_env%restart_section, extension=file_labels(2*i - 1), &
                                            file_form="UNFORMATTED", file_position="REWIND")
         CALL cp_fm_write_unformatted(workspace(1), rho_unit_nr)
         CALL cp_print_key_finished_output(rho_unit_nr, logger, rtbse_env%restart_section)
         ! Imag part
         rho_unit_nr = cp_print_key_unit_nr(logger, rtbse_env%restart_section, extension=file_labels(2*i), &
                                            file_form="UNFORMATTED", file_position="REWIND")
         CALL cp_fm_write_unformatted(workspace(2), rho_unit_nr)
         CALL cp_print_key_finished_output(rho_unit_nr, logger, rtbse_env%restart_section)
         ! Info
         rho_unit_nr = cp_print_key_unit_nr(logger, rtbse_env%restart_section, extension=".info", &
                                            file_form="UNFORMATTED", file_position="REWIND")
         IF (rho_unit_nr > 0) WRITE (rho_unit_nr) time_index
         CALL cp_print_key_finished_output(rho_unit_nr, logger, rtbse_env%restart_section)
      END DO
   END SUBROUTINE output_restart
! **************************************************************************************************
!> \brief Reads the density matrix from restart files and updates the starting time
!> \param restart_section Print key section for the restart files
!> \param rho Density matrix in AO basis
!> \param time_index Time index to be written into the info file
! **************************************************************************************************
   SUBROUTINE read_restart(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(len=default_path_length)                 :: save_name, save_name_2
      INTEGER                                            :: rho_unit_nr, j
      CHARACTER(len=17), DIMENSION(4)                    :: file_labels

      ! This allows the delta kick and output of moment at time 0 in all cases
      ! except the case when both imaginary and real parts of the density are read
      rtbse_env%restart_extracted = .FALSE.
      logger => cp_get_default_logger()
      ! Start by probing/loading info file
      save_name = cp_print_key_generate_filename(logger, rtbse_env%restart_section, extension=".info", my_local=.FALSE.)
      IF (file_exists(save_name)) THEN
         CALL open_file(save_name, file_status="OLD", file_form="UNFORMATTED", file_action="READ", &
                        unit_number=rho_unit_nr)
         READ (rho_unit_nr) rtbse_env%sim_start
         CALL close_file(rho_unit_nr)
         IF (rtbse_env%unit_nr > 0) WRITE (rtbse_env%unit_nr, '(A31,I25,A24)') " RTBSE| Starting from timestep ", &
            rtbse_env%sim_start, ", delta kick NOT applied"
      ELSE
         CPWARN("Restart required but no info file found - starting from sim_step given in input")
      END IF

      ! Default labels distinguishing up to two spin species and real/imaginary parts
      file_labels(1) = "_SPIN_A_RE.matrix"
      file_labels(2) = "_SPIN_A_IM.matrix"
      file_labels(3) = "_SPIN_B_RE.matrix"
      file_labels(4) = "_SPIN_B_IM.matrix"
      DO j = 1, rtbse_env%n_spin
         save_name = cp_print_key_generate_filename(logger, rtbse_env%restart_section, &
                                                    extension=file_labels(2*j - 1), my_local=.FALSE.)
         save_name_2 = cp_print_key_generate_filename(logger, rtbse_env%restart_section, &
                                                      extension=file_labels(2*j), my_local=.FALSE.)
         IF (file_exists(save_name) .AND. file_exists(save_name_2)) THEN
            CALL open_file(save_name, file_status="OLD", file_form="UNFORMATTED", file_action="READ", &
                           unit_number=rho_unit_nr)
            CALL cp_fm_read_unformatted(rtbse_env%real_workspace(1), rho_unit_nr)
            CALL close_file(rho_unit_nr)
            CALL open_file(save_name_2, file_status="OLD", file_form="UNFORMATTED", file_action="READ", &
                           unit_number=rho_unit_nr)
            CALL cp_fm_read_unformatted(rtbse_env%real_workspace(2), rho_unit_nr)
            CALL close_file(rho_unit_nr)
            CALL cp_fm_to_cfm(rtbse_env%real_workspace(1), rtbse_env%real_workspace(2), &
                              rtbse_env%rho(j))
            rtbse_env%restart_extracted = .TRUE.
         ELSE
            CPWARN("Restart without some restart matrices - starting from SCF density.")
         END IF
      END DO
   END SUBROUTINE read_restart
! **************************************************************************************************
!> \brief Linearized RT-BSE restart writer. Writes the restart set: lab-frame MO-active density
!>        matrices, the .info step index (sim_step = steps completed = the resume step), the
!>        once-per-run C_active gauge reference, and the appended RESTART.trace record feeding the
!>        FT prefix on continuation. All indices derive from sim_step so the .info/.trace
!>        bookkeeping cannot drift apart. (Full RTBSE uses the upstream output_restart above.)
!> \param rtbse_env RT-BSE environment
!> \param rho Density matrix (lab frame) to store
! **************************************************************************************************
   SUBROUTINE output_restart_linearized(rtbse_env, rho)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: workspace
      CHARACTER(len=17), DIMENSION(4)                    :: file_labels
      CHARACTER(len=16), DIMENSION(2)                    :: c_file_labels
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: rho_unit_nr, i

      ! Default labels distinguishing up to two spin species and real/imaginary parts
      file_labels(1) = "_SPIN_A_RE.matrix"
      file_labels(2) = "_SPIN_A_IM.matrix"
      file_labels(3) = "_SPIN_B_RE.matrix"
      file_labels(4) = "_SPIN_B_IM.matrix"
      c_file_labels(1) = "_SPIN_A_C.matrix"
      c_file_labels(2) = "_SPIN_B_C.matrix"

      logger => cp_get_default_logger()

      ! In linearized RT-BSE rho is MO-active sized; otherwise AO sized.
      IF (ASSOCIATED(rtbse_env%real_workspace_mo)) THEN
         workspace => rtbse_env%real_workspace_mo
      ELSE
         workspace => rtbse_env%real_workspace
      END IF

      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_to_fm(rho(i), workspace(1), workspace(2))
         ! Real part
         rho_unit_nr = cp_print_key_unit_nr(logger, rtbse_env%restart_section, extension=file_labels(2*i - 1), &
                                            file_form="UNFORMATTED", file_position="REWIND")
         CALL cp_fm_write_unformatted(workspace(1), rho_unit_nr)
         CALL cp_print_key_finished_output(rho_unit_nr, logger, rtbse_env%restart_section)
         ! Imag part
         rho_unit_nr = cp_print_key_unit_nr(logger, rtbse_env%restart_section, extension=file_labels(2*i), &
                                            file_form="UNFORMATTED", file_position="REWIND")
         CALL cp_fm_write_unformatted(workspace(2), rho_unit_nr)
         CALL cp_print_key_finished_output(rho_unit_nr, logger, rtbse_env%restart_section)
         ! Info
         rho_unit_nr = cp_print_key_unit_nr(logger, rtbse_env%restart_section, extension=".info", &
                                            file_form="UNFORMATTED", file_position="REWIND")
         IF (rho_unit_nr > 0) WRITE (rho_unit_nr) rtbse_env%sim_step
         CALL cp_print_key_finished_output(rho_unit_nr, logger, rtbse_env%restart_section)
      END DO

      ! Once-per-run C_active dump (linearized only): gauge reference for the restart basis bridge
      IF (ASSOCIATED(rtbse_env%real_workspace_mo) .AND. .NOT. rtbse_env%restart_C_written) THEN
         DO i = 1, rtbse_env%n_spin
            rho_unit_nr = cp_print_key_unit_nr(logger, rtbse_env%restart_section, extension=c_file_labels(i), &
                                               file_form="UNFORMATTED", file_position="REWIND")
            CALL cp_fm_write_unformatted(rtbse_env%C_active(i), rho_unit_nr)
            CALL cp_print_key_finished_output(rho_unit_nr, logger, rtbse_env%restart_section)
         END DO
         rtbse_env%restart_C_written = .TRUE.
      END IF

      CALL write_restart_trace(rtbse_env)
   END SUBROUTINE output_restart_linearized
! **************************************************************************************************
!> \brief Appends the current observable-trace record to RESTART.trace. Records are keyed by the
!>        observable slot n_slot = sim_step - sim_start_orig + 1 - the SAME index output_moments/
!>        output_field write (both drivers bump sim_step inside the propagation call before the
!>        output calls run). On the first call of a run the file is rewritten from memory (header +
!>        records 1..n_slot), truncating leftovers from an aborted run; subsequent calls append one
!>        record. Ionode writes; layout matches read_restart_trace verbatim.
!> \param rtbse_env RT-BSE environment
! **************************************************************************************************
   SUBROUTINE write_restart_trace(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(len=default_path_length)                 :: save_name
      INTEGER                                            :: trace_unit, n, n_eps, n_slot

      logger => cp_get_default_logger()
      save_name = cp_print_key_generate_filename(logger, rtbse_env%restart_section, extension=".trace", my_local=.FALSE.)

      IF (rtbse_env%unit_nr > 0) THEN
         n_slot = rtbse_env%sim_step - rtbse_env%sim_start_orig + 1
         IF (ASSOCIATED(rtbse_env%real_workspace_mo)) THEN
            ! scalar row count is safe: determine_active_mo_window takes the union window across
            ! spins, so eps_active is rectangular (mo_active, n_spin) by construction
            n_eps = SIZE(rtbse_env%eps_active, 1)
         ELSE
            n_eps = 0
         END IF
         IF (.NOT. rtbse_env%restart_trace_written) THEN
            CALL open_file(save_name, file_status="UNKNOWN", file_form="UNFORMATTED", file_action="WRITE", &
                           file_position="REWIND", unit_number=trace_unit)
            WRITE (trace_unit) restart_trace_version
            IF (ASSOCIATED(rtbse_env%real_workspace_mo)) THEN
               WRITE (trace_unit) rtbse_env%n_spin, rtbse_env%mo_active, rtbse_env%n_ao
            ELSE
               WRITE (trace_unit) rtbse_env%n_spin, rtbse_env%n_ao, rtbse_env%n_ao
            END IF
            WRITE (trace_unit) rtbse_env%sim_dt
            WRITE (trace_unit) REAL(rtbse_env%dft_control%rtp_control%delta_pulse_direction, dp), &
               rtbse_env%dft_control%rtp_control%delta_pulse_scale
            WRITE (trace_unit) n_eps
            IF (n_eps > 0) WRITE (trace_unit) rtbse_env%eps_active
            DO n = 1, n_slot
               WRITE (trace_unit) n, rtbse_env%time_trace(n), rtbse_env%field_trace(:, n), &
                  rtbse_env%moments_trace(:, :, n)
            END DO
         ELSE
            CALL open_file(save_name, file_status="OLD", file_form="UNFORMATTED", file_action="WRITE", &
                           file_position="APPEND", unit_number=trace_unit)
            WRITE (trace_unit) n_slot, rtbse_env%time_trace(n_slot), rtbse_env%field_trace(:, n_slot), &
               rtbse_env%moments_trace(:, :, n_slot)
         END IF
         CALL close_file(trace_unit)
      END IF
      rtbse_env%restart_trace_written = .TRUE.
   END SUBROUTINE write_restart_trace
! **************************************************************************************************
!> \brief Reads the RESTART.trace header that write_restart_trace writes, through the n_eps record
!> \param trace_unit unit of the trace, opened for reading at its start
!> \param ios nonzero when a record could not be read
!> \param version format version of the file; nothing after it is read when it differs from this binary's
!> \param n_spin_file number of spins of the original run
!> \param n_basis_file active-MO count (linearized) or AO count of the original run's density
!> \param n_ao_file AO count of the original run
!> \param dt_file time step of the original run, in atomic units
!> \param kick_dir_file delta-kick direction of the original run
!> \param kick_scale_file delta-kick scale of the original run
!> \param n_eps active eigenvalues per spin in the record that follows; the caller reads it when positive
!> \author Maximilian Graml
! **************************************************************************************************
   SUBROUTINE read_restart_trace_header(trace_unit, ios, version, n_spin_file, n_basis_file, n_ao_file, &
                                        dt_file, kick_dir_file, kick_scale_file, n_eps)
      INTEGER, INTENT(IN)                                :: trace_unit
      INTEGER, INTENT(OUT)                               :: ios, version, n_spin_file, n_basis_file, &
                                                            n_ao_file
      REAL(kind=dp), INTENT(OUT)                         :: dt_file
      REAL(kind=dp), DIMENSION(3), INTENT(OUT)           :: kick_dir_file
      REAL(kind=dp), INTENT(OUT)                         :: kick_scale_file
      INTEGER, INTENT(OUT)                               :: n_eps

      CHARACTER(len=*), PARAMETER :: routineN = 'read_restart_trace_header'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      n_spin_file = 0
      n_basis_file = 0
      n_ao_file = 0
      dt_file = 0.0_dp
      kick_dir_file(:) = 0.0_dp
      kick_scale_file = 0.0_dp
      n_eps = 0
      READ (trace_unit, IOSTAT=ios) version
      IF (ios == 0 .AND. version == restart_trace_version) THEN
         READ (trace_unit, IOSTAT=ios) n_spin_file, n_basis_file, n_ao_file
         IF (ios == 0) READ (trace_unit, IOSTAT=ios) dt_file
         IF (ios == 0) READ (trace_unit, IOSTAT=ios) kick_dir_file, kick_scale_file
         IF (ios == 0) READ (trace_unit, IOSTAT=ios) n_eps
         ! a failed read leaves n_eps undefined, and the caller reads an eps record when it is positive
         IF (ios /= 0) n_eps = 0
      END IF
      CALL timestop(handle)
   END SUBROUTINE read_restart_trace_header
! **************************************************************************************************
!> \brief Checks the linearized RT-BSE restart files before the SCF: their presence, the step in .info,
!>        the RESTART.trace header against the input, and the trace's record count. All ranks read.
!> \param force_env the force environment, for the input and the DFT control
!> \param restart_found true when INITIAL_WFN is RT_RESTART and the restart files exist
!> \author Maximilian Graml
! **************************************************************************************************
   SUBROUTINE check_restart_files_before_scf(force_env, restart_found)
      TYPE(force_env_type), POINTER                      :: force_env
      LOGICAL, INTENT(OUT)                               :: restart_found

      CHARACTER(len=*), PARAMETER :: routineN = 'check_restart_files_before_scf'
      CHARACTER(len=1), DIMENSION(2), PARAMETER          :: spin_letters = ["A", "B"]
      CHARACTER(len=10), DIMENSION(3), PARAMETER :: &
         matrix_suffixes = ["_RE.matrix", "_IM.matrix", "_C.matrix "]

      CHARACTER(len=default_path_length), DIMENSION(8)   :: file_names
      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:, :)     :: moments_record
      COMPLEX(kind=dp), DIMENSION(3)                     :: field_record
      INTEGER :: handle, i_file, i_missing, i_spin, info_unit, ios, n_ao_file, n_basis_file, n_eps, &
                 n_files, n_missing, n_records, n_spin, n_spin_file, n_steps, n_steps_input, record_index, &
                 sim_start, trace_unit, version
      LOGICAL                                            :: enforce_max_dt
      REAL(kind=dp)                                      :: dt_file, dt_input, kick_scale_file, &
                                                            time_record
      REAL(kind=dp), DIMENSION(3)                        :: kick_dir_file
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input, md_section, restart_section

      CALL timeset(routineN, handle)
      NULLIFY (dft_control, input)
      CALL get_qs_env(force_env%qs_env, dft_control=dft_control, input=input)
      restart_found = .FALSE.
      IF (dft_control%rtp_control%initial_wfn /= use_rt_restart) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! the names the readers open: .info, .trace, then RE, IM and C per spin
      logger => cp_get_default_logger()
      restart_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION%PRINT%RESTART")
      n_spin = dft_control%nspins
      file_names(1) = cp_print_key_generate_filename(logger, restart_section, extension=".info", my_local=.FALSE.)
      file_names(2) = cp_print_key_generate_filename(logger, restart_section, extension=".trace", my_local=.FALSE.)
      n_files = 2
      DO i_spin = 1, n_spin
         DO i_file = 1, 3
            n_files = n_files + 1
            file_names(n_files) = cp_print_key_generate_filename(logger, restart_section, &
                                                                 extension="_SPIN_"//spin_letters(i_spin)// &
                                                                 TRIM(matrix_suffixes(i_file)), my_local=.FALSE.)
         END DO
      END DO
      i_missing = 0
      n_missing = 0
      DO i_file = 1, n_files
         IF (.NOT. file_exists(file_names(i_file))) THEN
            n_missing = n_missing + 1
            IF (i_missing == 0) i_missing = i_file
         END IF
      END DO
      IF (n_missing == n_files) THEN
         CALL cp_warn(__LOCATION__, &
                      "INITIAL_WFN RT_RESTART was requested, but no restart file exists ("// &
                      TRIM(file_names(1))//" not found); the calculation continues as a fresh run from "// &
                      "the SCF density with the delta kick")
         CALL timestop(handle)
         RETURN
      END IF
      IF (n_missing > 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "RT_RESTART: restart file not found: "//TRIM(file_names(i_missing))// &
                       "; restore the first run's restart files beside the input, or remove "// &
                       "DFT%REAL_TIME_PROPAGATION%INITIAL_WFN RT_RESTART")
      END IF
      restart_found = .TRUE.

      CALL open_file(file_names(1), file_status="OLD", file_form="UNFORMATTED", file_action="READ", &
                     unit_number=info_unit)
      READ (info_unit, IOSTAT=ios) sim_start
      CALL close_file(info_unit)
      IF (ios /= 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "RT_RESTART: cannot read the step from "//TRIM(file_names(1))// &
                       "; restore the restart files of one completed restart step of the first run")
      END IF

      CALL open_file(file_names(2), file_status="OLD", file_form="UNFORMATTED", file_action="READ", &
                     unit_number=trace_unit)
      CALL read_restart_trace_header(trace_unit, ios, version, n_spin_file, n_basis_file, n_ao_file, &
                                     dt_file, kick_dir_file, kick_scale_file, n_eps)
      ! the active eigenvalues close the header
      IF (ios == 0 .AND. version == restart_trace_version .AND. n_eps > 0) READ (trace_unit, IOSTAT=ios)
      IF (ios /= 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "RT_RESTART: cannot read the RESTART.trace header of "//TRIM(file_names(2))// &
                       "; restore the restart files of one completed restart step of the first run")
      END IF
      IF (version /= restart_trace_version) THEN
         CALL cp_abort(__LOCATION__, &
                       "RT_RESTART: RESTART.trace format version differs from this binary's ("// &
                       TRIM(cp_to_string(version))//" against "//TRIM(cp_to_string(restart_trace_version))// &
                       "); continue the run with the CP2K version that wrote it")
      END IF
      IF (n_spin_file /= n_spin) THEN
         CALL cp_abort(__LOCATION__, &
                       "RT_RESTART: RESTART.trace was written for another number of spins ("// &
                       TRIM(cp_to_string(n_spin_file))//", this input has "//TRIM(cp_to_string(n_spin))// &
                       "); set DFT%UKS as in the first run")
      END IF
      md_section => section_vals_get_subs_vals(force_env%root_section, "MOTION%MD")
      CALL section_vals_val_get(md_section, "TIMESTEP", r_val=dt_input)
      CALL section_vals_val_get(md_section, "STEPS", i_val=n_steps_input)
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%RTBSE%ENFORCE_MAX_DT", l_val=enforce_max_dt)
      ! ENFORCE_MAX_DT inherits the stored time step instead
      IF (.NOT. enforce_max_dt .AND. ABS(dt_file - dt_input) > 1.0e-12_dp*MAX(1.0_dp, ABS(dt_input))) THEN
         CALL cp_abort(__LOCATION__, &
                       "RT_RESTART: MOTION%MD%TIMESTEP differs from the original run; set it to "// &
                       TRIM(ADJUSTL(cp_to_string(dt_file*femtoseconds, '(ES16.9)')))//" fs, or enable "// &
                       "DFT%REAL_TIME_PROPAGATION%RTBSE%ENFORCE_MAX_DT to inherit it")
      END IF
      IF (MAXVAL(ABS(kick_dir_file - REAL(dft_control%rtp_control%delta_pulse_direction, dp))) > 1.0e-12_dp .OR. &
          ABS(kick_scale_file - dft_control%rtp_control%delta_pulse_scale) > 1.0e-12_dp) THEN
         CALL cp_abort(__LOCATION__, &
                       "RT_RESTART: delta kick differs from the stored kick; set "// &
                       "DFT%REAL_TIME_PROPAGATION%DELTA_PULSE_DIRECTION and DELTA_PULSE_SCALE as in the first run")
      END IF
      ! with ENFORCE_MAX_DT, initialize_maximum_timestep keeps the window STEPS*TIMESTEP at the stored step
      n_steps = n_steps_input
      IF (enforce_max_dt .AND. dt_file > 0.0_dp) THEN
         n_steps = MAX(1, NINT(REAL(MAX(0, n_steps_input), dp)*dt_input/dt_file))
      END IF
      IF (sim_start >= n_steps) THEN
         CALL cp_abort(__LOCATION__, &
                       "RT_RESTART: the restart step is not below MOTION%MD%STEPS: the restart files hold step "// &
                       TRIM(cp_to_string(sim_start))//"; raise MOTION%MD%STEPS")
      END IF

      ! the continuation reloads records 1 to sim_start + 1, each read with the reader's item list
      ALLOCATE (moments_record(n_spin, 3))
      n_records = 0
      DO
         READ (trace_unit, IOSTAT=ios) record_index, time_record, field_record, moments_record
         IF (ios /= 0) EXIT
         n_records = n_records + 1
      END DO
      DEALLOCATE (moments_record)
      CALL close_file(trace_unit)
      IF (n_records < sim_start + 1) THEN
         CALL cp_abort(__LOCATION__, &
                       "RT_RESTART: RESTART.trace holds fewer records than the restart step needs ("// &
                       TRIM(cp_to_string(n_records))//" for step "//TRIM(cp_to_string(sim_start))// &
                       "); a run killed while writing its restart files leaves this state: restore the "// &
                       "restart files of one completed restart step")
      END IF

      CALL timestop(handle)
   END SUBROUTINE check_restart_files_before_scf
! **************************************************************************************************
!> \brief Early phase of the restart read: the starting step index from the .info file, the original
!>        run's dt peeked from the RESTART.trace header (so ENFORCE_MAX_DT can inherit it), and an
!>        SCF_GUESS hygiene check. Runs BEFORE initialize_maximum_timestep; the trace prefix records
!>        are loaded separately by read_restart_trace once the trace arrays are sized.
!> \param rtbse_env RT-BSE environment
! **************************************************************************************************
   SUBROUTINE read_restart_info(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(len=default_path_length)                 :: save_name
      INTEGER                                            :: info_unit, ios, n_ao_file, n_basis_file, &
                                                            n_eps, n_spin_file, trace_unit, version
      REAL(kind=dp)                                      :: dt_file, kick_scale_file
      REAL(kind=dp), DIMENSION(3)                        :: kick_dir_file
      TYPE(scf_control_type), POINTER                    :: scf_control

      logger => cp_get_default_logger()

      ! check_restart_files_before_scf found both files readable before the SCF
      save_name = cp_print_key_generate_filename(logger, rtbse_env%restart_section, extension=".info", my_local=.FALSE.)
      CALL open_file(save_name, file_status="OLD", file_form="UNFORMATTED", file_action="READ", &
                     unit_number=info_unit)
      READ (info_unit) rtbse_env%sim_start
      CALL close_file(info_unit)
      IF (rtbse_env%unit_nr > 0) WRITE (rtbse_env%unit_nr, '(A31,I25,A24)') " RTBSE| Starting from timestep ", &
         rtbse_env%sim_start, ", delta kick NOT applied"

      ! the original run's dt, which ENFORCE_MAX_DT inherits rather than recompute a window-dependent dt
      save_name = cp_print_key_generate_filename(logger, rtbse_env%restart_section, extension=".trace", my_local=.FALSE.)
      CALL open_file(save_name, file_status="OLD", file_form="UNFORMATTED", file_action="READ", &
                     unit_number=trace_unit)
      CALL read_restart_trace_header(trace_unit, ios, version, n_spin_file, n_basis_file, n_ao_file, &
                                     dt_file, kick_dir_file, kick_scale_file, n_eps)
      CALL close_file(trace_unit)
      ! check_restart_files_before_scf read this header whole before the SCF
      CPASSERT(ios == 0)
      rtbse_env%sim_dt_restart = dt_file

      ! Hygiene nudge only - correctness is protected by the restart basis bridge (linearized path)
      NULLIFY (scf_control)
      CALL get_qs_env(rtbse_env%qs_env, scf_control=scf_control)
      IF (scf_control%density_guess /= restart_guess) THEN
         CPWARN("RT_RESTART without SCF_GUESS RESTART - SCF may reconverge to a gauge-rotated MO basis.")
      END IF
   END SUBROUTINE read_restart_info
! **************************************************************************************************
!> \brief Reads the RESTART.trace prefix (records 1..sim_start) into the in-memory moment/field/
!>        time traces so the continuation FT covers the full history. check_restart_files_before_scf
!>        checked the header against the input and the record count before the SCF; this reader
!>        checks what needs the SCF: the AO and active-MO counts, and eps_active > 0.1 meV (warn).
!>        All ranks read (the traces are replicated). Record layout matches write_restart_trace verbatim.
!> \param rtbse_env RT-BSE environment
! **************************************************************************************************
   SUBROUTINE read_restart_trace(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(len=default_path_length)                 :: save_name
      INTEGER                                            :: trace_unit, version, n_spin_file, &
                                                            n_basis_file, n_ao_file, n_eps, n, ios
      REAL(kind=dp)                                      :: dt_file, kick_scale_file, t_rec
      REAL(kind=dp), DIMENSION(3)                        :: kick_dir_file
      REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE        :: eps_file
      COMPLEX(kind=dp), DIMENSION(3)                     :: field_rec
      COMPLEX(kind=dp), DIMENSION(:, :), ALLOCATABLE     :: mom_rec

      logger => cp_get_default_logger()
      save_name = cp_print_key_generate_filename(logger, rtbse_env%restart_section, extension=".trace", my_local=.FALSE.)
      CALL open_file(save_name, file_status="OLD", file_form="UNFORMATTED", file_action="READ", &
                     unit_number=trace_unit)
      CALL read_restart_trace_header(trace_unit, ios, version, n_spin_file, n_basis_file, n_ao_file, &
                                     dt_file, kick_dir_file, kick_scale_file, n_eps)
      ! check_restart_files_before_scf read this header whole before the SCF
      CPASSERT(ios == 0)
      IF (n_ao_file /= rtbse_env%n_ao) CPABORT("RESTART.trace: n_ao mismatch")
      ! v1 is strict same-window: the linearized active-MO count must match (D6)
      IF (ASSOCIATED(rtbse_env%real_workspace_mo) .AND. n_basis_file /= rtbse_env%mo_active) THEN
         CALL cp_abort(__LOCATION__, &
                       "RESTART.trace: active-MO count differs from the original run - same-window continuation only")
      END IF
      IF (n_eps > 0) THEN
         ALLOCATE (eps_file(n_eps, rtbse_env%n_spin))
         READ (trace_unit) eps_file
         IF (ASSOCIATED(rtbse_env%real_workspace_mo)) THEN
            IF (n_eps /= SIZE(rtbse_env%eps_active, 1)) THEN
               CALL cp_abort(__LOCATION__, "RESTART.trace: active-window size mismatch")
            END IF
            ! eps_active is not populated yet (built later in the Hamiltonian init); stash for the
            ! post-Hamiltonian consistency check in check_restart_eps_consistency.
            ALLOCATE (rtbse_env%eps_active_restart(n_eps, rtbse_env%n_spin))
            rtbse_env%eps_active_restart(:, :) = eps_file
         END IF
         DEALLOCATE (eps_file)
      END IF

      ALLOCATE (mom_rec(rtbse_env%n_spin, 3))
      DO
         READ (trace_unit, IOSTAT=ios) n, t_rec, field_rec, mom_rec
         IF (ios /= 0) EXIT
         ! Fill through slot sim_start+1: both drivers bump sim_step inside the propagation call
         ! (etrs_scf_loop / solve_rk4_timestep), so the live loop's first step writes slot
         ! sim_start+2 and slot sim_start+1 must be reloaded here.
         IF (n > rtbse_env%sim_start + 1) EXIT
         IF (n > SIZE(rtbse_env%time_trace)) THEN
            CPABORT("RESTART.trace: record index exceeds trace size - increase MOTION%MD%STEPS")
         END IF
         rtbse_env%time_trace(n) = t_rec
         rtbse_env%field_trace(:, n) = field_rec
         rtbse_env%moments_trace(:, :, n) = mom_rec
      END DO
      DEALLOCATE (mom_rec)
      CALL close_file(trace_unit)
   END SUBROUTINE read_restart_trace
! **************************************************************************************************
!> \brief Compares the original run's active eigenvalues (stashed by read_restart_trace) against
!>        the recomputed eps_active, once the Hamiltonian is built. Prints the max deviation and
!>        warns above 0.1 meV (GW analytic continuation gives run-to-run QP noise above FP, so the
!>        threshold is deliberately loose). Frees the stash. No-op if nothing was stashed.
!> \param rtbse_env RT-BSE environment
! **************************************************************************************************
   SUBROUTINE check_restart_eps_consistency(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      REAL(kind=dp)                                      :: eps_dev

      IF (.NOT. ASSOCIATED(rtbse_env%eps_active_restart)) RETURN

      eps_dev = MAXVAL(ABS(rtbse_env%eps_active_restart - rtbse_env%eps_active))
      IF (rtbse_env%unit_nr > 0) WRITE (rtbse_env%unit_nr, '(A,ES12.3,A)') &
         " RTBSE| Restart eps_active max deviation vs original run ", eps_dev, " Ha"
      ! 0.1 meV = 1.0e-4 eV, converted to Ha via evolt (Ha -> eV factor)
      IF (eps_dev > 1.0e-4_dp/evolt) THEN
         CALL cp_warn(__LOCATION__, &
                      "RESTART.trace: active eigenvalues deviate beyond 0.1 meV - "// &
                      "Hamiltonian changed; continuation is physically inconsistent.")
      END IF

      DEALLOCATE (rtbse_env%eps_active_restart)
      NULLIFY (rtbse_env%eps_active_restart)
   END SUBROUTINE check_restart_eps_consistency
! **************************************************************************************************
!> \brief Late phase of the restart read: overwrites rho from the lab-frame restart matrices and
!>        sets restart_extracted. MO-active for the linearized path, AO otherwise
!>        (cp_fm_read_unformatted aborts on a size mismatch, guarding a changed active window).
!> \param rtbse_env RT-BSE environment
! **************************************************************************************************
   SUBROUTINE read_restart_density(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(len=default_path_length)                 :: save_name, save_name_2
      INTEGER                                            :: rho_unit_nr, j
      CHARACTER(len=17), DIMENSION(4)                    :: file_labels
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: ws

      ! This allows the delta kick and output of moment at time 0 in all cases
      ! except the case when both imaginary and real parts of the density are read
      rtbse_env%restart_extracted = .FALSE.
      logger => cp_get_default_logger()

      IF (ASSOCIATED(rtbse_env%real_workspace_mo)) THEN
         ws => rtbse_env%real_workspace_mo
      ELSE
         ws => rtbse_env%real_workspace
      END IF

      ! Default labels distinguishing up to two spin species and real/imaginary parts
      file_labels(1) = "_SPIN_A_RE.matrix"
      file_labels(2) = "_SPIN_A_IM.matrix"
      file_labels(3) = "_SPIN_B_RE.matrix"
      file_labels(4) = "_SPIN_B_IM.matrix"
      DO j = 1, rtbse_env%n_spin
         save_name = cp_print_key_generate_filename(logger, rtbse_env%restart_section, &
                                                    extension=file_labels(2*j - 1), my_local=.FALSE.)
         save_name_2 = cp_print_key_generate_filename(logger, rtbse_env%restart_section, &
                                                      extension=file_labels(2*j), my_local=.FALSE.)
         ! check_restart_files_before_scf found both files before the SCF
         CALL open_file(save_name, file_status="OLD", file_form="UNFORMATTED", file_action="READ", &
                        unit_number=rho_unit_nr)
         CALL cp_fm_read_unformatted(ws(1), rho_unit_nr)
         CALL close_file(rho_unit_nr)
         CALL open_file(save_name_2, file_status="OLD", file_form="UNFORMATTED", file_action="READ", &
                        unit_number=rho_unit_nr)
         CALL cp_fm_read_unformatted(ws(2), rho_unit_nr)
         CALL close_file(rho_unit_nr)
         CALL cp_fm_to_cfm(ws(1), ws(2), &
                           rtbse_env%rho(j))
         rtbse_env%restart_extracted = .TRUE.
      END DO
   END SUBROUTINE read_restart_density
! **************************************************************************************************
!> \brief Reads the previous run's C_active slabs (gauge reference for the restart basis bridge).
!> \param rtbse_env RT-BSE environment
!> \param C_old Caller-created fm array (n_spin) on fm_struct_ao_mo_active, filled here
! **************************************************************************************************
   SUBROUTINE read_restart_C(rtbse_env, C_old)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: C_old
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(len=default_path_length)                 :: save_name
      CHARACTER(len=16), DIMENSION(2)                    :: c_file_labels
      INTEGER                                            :: c_unit, j

      c_file_labels(1) = "_SPIN_A_C.matrix"
      c_file_labels(2) = "_SPIN_B_C.matrix"
      logger => cp_get_default_logger()
      ! check_restart_files_before_scf found every file before the SCF
      DO j = 1, rtbse_env%n_spin
         save_name = cp_print_key_generate_filename(logger, rtbse_env%restart_section, &
                                                    extension=c_file_labels(j), my_local=.FALSE.)
         CALL open_file(save_name, file_status="OLD", file_form="UNFORMATTED", file_action="READ", &
                        unit_number=c_unit)
         CALL cp_fm_read_unformatted(C_old(j), c_unit)
         CALL close_file(c_unit)
      END DO
   END SUBROUTINE read_restart_C
END MODULE rt_bse_io

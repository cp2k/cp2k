!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Input/output from the propagation via RT-BSE method.
!> \author Stepan Marek (08.24)
! **************************************************************************************************

MODULE rt_bse_io
   USE cp_fm_types, ONLY: cp_fm_type, &
                          cp_fm_p_type, &
                          cp_fm_create, &
                          cp_fm_release, &
                          cp_fm_read_unformatted, &
                          cp_fm_write_unformatted, &
                          cp_fm_write_formatted
   USE cp_cfm_types, ONLY: cp_cfm_type, &
                           cp_cfm_p_type, &
                           cp_fm_to_cfm, &
                           cp_cfm_to_fm
   USE kinds, ONLY: dp, &
                    default_path_length
   USE cp_fm_basic_linalg, ONLY: cp_fm_trace, &
                                 cp_fm_transpose, &
                                 cp_fm_norm
   USE parallel_gemm_api, ONLY: parallel_gemm
   USE cp_log_handling, ONLY: cp_logger_type, &
                              cp_get_default_logger
   USE cp_output_handling, ONLY: cp_print_key_unit_nr, &
                                 cp_print_key_finished_output, &
                                 cp_print_key_generate_filename, &
                                 low_print_level, &
                                 medium_print_level
   USE input_section_types, ONLY: section_vals_type
   USE rt_bse_types, ONLY: rtbse_env_type, &
                           multiply_cfm_fm, &
                           multiply_fm_cfm
   USE cp_files, ONLY: open_file, &
                       file_exists, &
                       close_file
   USE input_constants, ONLY: do_exact, &
                              do_bch, &
                              rtp_bse_ham_g0w0, &
                              rtp_bse_ham_ks, &
                              use_rt_restart
   USE physcon, ONLY: femtoseconds, &
                      evolt
   USE mathconstants, ONLY: twopi

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_bse_io"

   #:include "rt_bse_macros.fypp"

   PUBLIC :: output_moments, &
             read_moments, &
             output_moments_ft, &
             output_polarizability, &
             output_field, &
             read_field, &
             output_mos_contravariant, &
             output_mos_covariant, &
             output_restart, &
             read_restart, &
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
      CASE (rtp_bse_ham_g0w0)
         WRITE (rtbse_env%unit_nr, '(A58)') 'G0W0'
      CASE (rtp_bse_ham_ks)
         WRITE (rtbse_env%unit_nr, '(A58)') 'Kohn-Sham'
      END SELECT

      WRITE (rtbse_env%unit_nr, '(A18,L62)') ' Apply delta pulse', &
         rtbse_env%dft_control%rtp_control%apply_delta_pulse

      WRITE (rtbse_env%unit_nr, '(A)') ''

   END SUBROUTINE

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

   END SUBROUTINE
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

   END SUBROUTINE
! **************************************************************************************************
!> \brief Writes the header for the etrs iteration updates - only for log level > medium
!> \param rtbse_env Entry point - rtbse environment
! **************************************************************************************************
   SUBROUTINE print_timestep_info(rtbse_env, step, convergence, electron_num_re, etrs_num)
      TYPE(rtbse_env_type)                               :: rtbse_env
      INTEGER                                            :: step
      REAL(kind=dp)                                      :: convergence
      REAL(kind=dp)                                      :: electron_num_re
      INTEGER                                            :: etrs_num
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => cp_get_default_logger()

      IF (logger%iter_info%print_level > low_print_level .AND. rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, '(A23,A20,A20,A17)') " RTBSE| Simulation step", "Convergence", &
            "Electron number", "ETRS Iterations"
         WRITE (rtbse_env%unit_nr, '(A7,I16,E20.8E3,E20.8E3,I17)') ' RTBSE|', step, convergence, &
            electron_num_re, etrs_num
      END IF

   END SUBROUTINE

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
   END SUBROUTINE
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
   END SUBROUTINE
! **************************************************************************************************
!> \brief Prints the current field components into a file provided by input
!> \param rtbse_env Entry point - gwbse environment
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE output_field(rtbse_env)
      TYPE(rtbse_env_type)                               :: rtbse_env
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: field_unit, n, i

      ! Get logger
      logger => cp_get_default_logger()
      ! Get file descriptor
      field_unit = cp_print_key_unit_nr(logger, rtbse_env%field_section, extension=".dat")
      ! If the file descriptor is non-zero, output field
      ! TODO : Output also in SI
      IF (field_unit /= -1) THEN
         WRITE (field_unit, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') rtbse_env%sim_time*femtoseconds, &
            rtbse_env%field(1), rtbse_env%field(2), rtbse_env%field(3)
      END IF
      ! Write the output to memory for FT
      ! Need the absolute index
      n = rtbse_env%sim_step - rtbse_env%sim_start_orig + 1
      DO i = 1, 3
         rtbse_env%field_trace(i)%series(n) = rtbse_env%field(i)
      END DO
      rtbse_env%time_trace%series(n) = rtbse_env%sim_time
      CALL cp_print_key_finished_output(field_unit, logger, rtbse_env%field_section)

   END SUBROUTINE
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

      ! Get logger
      logger => cp_get_default_logger()
      ! Get file name
      save_name = cp_print_key_generate_filename(logger, rtbse_env%field_section, extension=".dat", my_local=.FALSE.)
      IF (file_exists(save_name)) THEN
         CALL open_file(save_name, file_status="OLD", file_form="FORMATTED", file_action="READ", &
                        unit_number=field_unit)
         DO k = rtbse_env%sim_start_orig, rtbse_env%sim_start
            n = k - rtbse_env%sim_start_orig + 1
            READ (field_unit, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') rtbse_env%time_trace%series(n), &
               rtbse_env%field_trace(1)%series(n), rtbse_env%field_trace(2)%series(n), rtbse_env%field_trace(3)%series(n)
            ! Set the time units back to atomic units
            rtbse_env%time_trace%series(n) = rtbse_env%time_trace%series(n)/femtoseconds
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
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: i, j, n, moments_unit_re, moments_unit_im
      CHARACTER(len=14), DIMENSION(4)                    :: file_labels
      REAL(kind=dp), DIMENSION(3)                        :: moments

      ! Start by getting the relevant file unit
      file_labels(1) = "_SPIN_A_RE.dat"
      file_labels(2) = "_SPIN_A_IM.dat"
      file_labels(3) = "_SPIN_B_RE.dat"
      file_labels(4) = "_SPIN_B_IM.dat"
      logger => cp_get_default_logger()
      DO j = 1, rtbse_env%n_spin
         moments_unit_re = cp_print_key_unit_nr(logger, rtbse_env%moments_section, extension=file_labels(2*j - 1))
         moments_unit_im = cp_print_key_unit_nr(logger, rtbse_env%moments_section, extension=file_labels(2*j))
         ! If, for any reason, the file unit is not provided, skip to next cycle immediately
         ! TODO : Specify output units in config
         ! Need to transpose due to the definition of trace function
         CALL cp_cfm_to_fm(msource=rho(j), mtargetr=rtbse_env%real_workspace(2))
         DO i = 1, 3
            ! Moments should be symmetric, test without transopose?
            CALL cp_fm_transpose(rtbse_env%moments(i), rtbse_env%real_workspace(1))
            CALL cp_fm_trace(rtbse_env%real_workspace(1), rtbse_env%real_workspace(2), moments(i))
            ! Scale by spin degeneracy and electron charge
            moments(i) = -moments(i)*rtbse_env%spin_degeneracy
         END DO
         ! Output to the file
         IF (moments_unit_re > 0) THEN
            ! If outputting to standard output, also prepend identifying 'MOMENTS_TRACE|' string
            IF (rtbse_env%unit_nr == moments_unit_re) THEN
               WRITE (moments_unit_re, '(A18,E20.8E3,E20.8E3,E20.8E3,E20.8E3)') &
                  " MOMENTS_TRACE_RE|", rtbse_env%sim_time*femtoseconds, moments(1), moments(2), moments(3)
            ELSE
               WRITE (moments_unit_re, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') &
                  rtbse_env%sim_time*femtoseconds, moments(1), moments(2), moments(3)
               CALL cp_print_key_finished_output(moments_unit_re, logger, rtbse_env%moments_section)
            END IF
         END IF
         ! Save to memory for FT - real part
         n = rtbse_env%sim_step - rtbse_env%sim_start_orig + 1
         DO i = 1, 3
            rtbse_env%moments_trace(i)%series(n) = CMPLX(moments(i), 0.0, kind=dp)
         END DO
         ! Same for imaginary part
         CALL cp_cfm_to_fm(msource=rho(j), mtargeti=rtbse_env%real_workspace(2))
         DO i = 1, 3
            CALL cp_fm_transpose(rtbse_env%moments(i), rtbse_env%real_workspace(1))
            CALL cp_fm_trace(rtbse_env%real_workspace(1), rtbse_env%real_workspace(2), moments(i))
            ! Scale by spin degeneracy and electron charge
            moments(i) = -moments(i)*rtbse_env%spin_degeneracy
         END DO
         ! Output to the file
         IF (moments_unit_im > 0) THEN
            ! If outputting to standard output, also prepend identifying 'MOMENTS_TRACE|' string
            IF (rtbse_env%unit_nr == moments_unit_im) THEN
               WRITE (moments_unit_im, '(A18,E20.8E3,E20.8E3,E20.8E3,E20.8E3)') &
                  " MOMENTS_TRACE_IM|", rtbse_env%sim_time*femtoseconds, moments(1), moments(2), moments(3)
            ELSE
               WRITE (moments_unit_im, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') &
                  rtbse_env%sim_time*femtoseconds, moments(1), moments(2), moments(3)
               ! Close the files
               CALL cp_print_key_finished_output(moments_unit_im, logger, rtbse_env%moments_section)
            END IF
         END IF
         ! Save to memory for FT - imag part
         DO i = 1, 3
            rtbse_env%moments_trace(i)%series(n) = rtbse_env%moments_trace(i)%series(n) + CMPLX(0.0, moments(i), kind=dp)
         END DO
      END DO
   END SUBROUTINE
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
!> \brief Reads the moments and time traces from the save files
!> \param rtbse_env GW-BSE environment (assumes consistent setup, i.e. a continuation run).
!>                  Otherwise, the traces are set at zero
! **************************************************************************************************
   SUBROUTINE read_moments(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(len=default_path_length)                 :: save_name, save_name_2
      INTEGER                                            :: i, j, k, moments_unit_re, moments_unit_im, n
      CHARACTER(len=17), DIMENSION(4)                    :: file_labels
      REAL(kind=dp), DIMENSION(3)                        :: moments_re, moments_im
      REAL(kind=dp)                                      :: timestamp

      logger => cp_get_default_logger()

      ! Read moments from the previous run
      ! Default labels distinguishing up to two spin species and real/imaginary parts
      file_labels(1) = "_SPIN_A_RE.dat"
      file_labels(2) = "_SPIN_A_IM.dat"
      file_labels(3) = "_SPIN_B_RE.dat"
      file_labels(4) = "_SPIN_B_IM.dat"
      DO j = 1, rtbse_env%n_spin
         save_name = cp_print_key_generate_filename(logger, rtbse_env%moments_section, &
                                                    extension=file_labels(2*j - 1), my_local=.FALSE.)
         save_name_2 = cp_print_key_generate_filename(logger, rtbse_env%moments_section, &
                                                      extension=file_labels(2*j), my_local=.FALSE.)
         IF (file_exists(save_name) .AND. file_exists(save_name_2)) THEN
            CALL open_file(save_name, file_status="OLD", file_form="FORMATTED", file_action="READ", &
                           unit_number=moments_unit_re)
            CALL open_file(save_name_2, file_status="OLD", file_form="FORMATTED", file_action="READ", &
                           unit_number=moments_unit_im)
            ! Extra time step for the initial one
            DO k = rtbse_env%sim_start_orig, rtbse_env%sim_start
               ! Determine the absolute time step - offset in memory
               n = k - rtbse_env%sim_start_orig + 1
               READ (moments_unit_re, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') timestamp, &
                  moments_re(1), moments_re(2), moments_re(3)
               READ (moments_unit_im, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') timestamp, &
                  moments_im(1), moments_im(2), moments_im(3)
               DO i = 1, 3
                  rtbse_env%moments_trace(i)%series(n) = CMPLX(moments_re(i), moments_im(i), kind=dp)
               END DO
            END DO
            ! Change back to atomic units in the trace
            rtbse_env%time_trace%series(:) = rtbse_env%time_trace%series(:)/femtoseconds
            CALL close_file(moments_unit_re)
            CALL close_file(moments_unit_im)
         ELSE IF (rtbse_env%dft_control%rtp_control%initial_wfn == use_rt_restart) THEN
            CPWARN("Restart without previous moments - unknown moments trace set to zero.")
         END IF
      END DO
   END SUBROUTINE read_moments
! **************************************************************************************************
!> \brief Outputs the Fourier transform of moments stored in the environment memory to the configured file
!> \param rtbse_env GW-BSE environment
! **************************************************************************************************
   SUBROUTINE output_moments_ft(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_logger_type), POINTER                      :: logger
      REAL(kind=dp), DIMENSION(:), POINTER               :: omega_series, &
                                                            ft_real_series, &
                                                            ft_imag_series, &
                                                            value_series
      ! Stores the data in ready for output format
      !  - first dimension is 6 - 1 - real part along x, 2 - imag part along x, 3 - real part along y, ...
      REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE         :: ft_full_series
      INTEGER                                            :: i, n, ft_unit

      logger => cp_get_default_logger()

      n = rtbse_env%sim_nsteps + 2
      NULLIFY (omega_series)
      ALLOCATE (omega_series(n), source=0.0_dp)
      NULLIFY (ft_real_series)
      ALLOCATE (ft_real_series(n), source=0.0_dp)
      NULLIFY (ft_imag_series)
      ALLOCATE (ft_imag_series(n), source=0.0_dp)
      NULLIFY (value_series)
      ALLOCATE (value_series(n), source=0.0_dp)
      ALLOCATE (ft_full_series(6, n))
      ! Carry out for each direction independently and real and imaginary parts also independently
      DO i = 1, 3
         ! Real part of the value first
         value_series(:) = REAL(rtbse_env%moments_trace(i)%series(:))
         CALL ft_simple(rtbse_env%time_trace%series, value_series, omega_series, ft_real_series, ft_imag_series, &
                        rtbse_env%ft_damping, rtbse_env%ft_start)
         ft_full_series(2*i - 1, :) = ft_real_series(:)
         ft_full_series(2*i, :) = ft_imag_series(:)
         ! Now imaginary part
         value_series(:) = AIMAG(rtbse_env%moments_trace(i)%series(:))
         CALL ft_simple(rtbse_env%time_trace%series, value_series, omega_series, ft_real_series, ft_imag_series, &
                        rtbse_env%ft_damping, rtbse_env%ft_start)
         ft_full_series(2*i - 1, :) = ft_full_series(2*i - 1, :) - ft_imag_series
         ft_full_series(2*i, :) = ft_full_series(2*i, :) + ft_real_series
      END DO
      DEALLOCATE (ft_real_series)
      DEALLOCATE (ft_imag_series)
      DEALLOCATE (value_series)
      ! Now, write these to file
      ft_unit = cp_print_key_unit_nr(logger, rtbse_env%ft_section, extension=".dat", &
                                     file_form="FORMATTED", file_position="REWIND")
      IF (ft_unit > 0) THEN
         DO i = 1, n
            WRITE (ft_unit, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3,E20.8E3,E20.8E3,E20.8E3)') &
               omega_series(i), ft_full_series(1, i), ft_full_series(2, i), &
               ft_full_series(3, i), ft_full_series(4, i), &
               ft_full_series(5, i), ft_full_series(6, i)
         END DO
      END IF
      CALL cp_print_key_finished_output(ft_unit, logger, rtbse_env%ft_section)
      DEALLOCATE (omega_series)
      DEALLOCATE (ft_full_series)
   END SUBROUTINE output_moments_ft
! **************************************************************************************************
!> \brief Outputs the isotropic polarizability tensor element alpha _ ij = mu_i(omega)/E_j(omega),
!>        where i and j are provided by the configuration. The tensor element is energy dependent and
!>        has real and imaginary parts
!> \param rtbse_env GW-BSE environment
! **************************************************************************************************
   SUBROUTINE output_polarizability(rtbse_env)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_logger_type), POINTER                      :: logger
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: omega_series, &
                                                            ft_real_series, &
                                                            ft_imag_series, &
                                                            value_series
      COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE        :: moment_series, &
                                                            field_series, &
                                                            polarizability_series
      INTEGER                                            :: pol_unit, &
                                                            i, n

      logger => cp_get_default_logger()

      n = rtbse_env%sim_nsteps + 2
      ! All allocations together, although could save some memory, if required by consequent deallocations
      ALLOCATE (omega_series(n), source=0.0_dp)
      ALLOCATE (ft_real_series(n), source=0.0_dp)
      ALLOCATE (ft_imag_series(n), source=0.0_dp)
      ALLOCATE (value_series(n), source=0.0_dp)
      ALLOCATE (moment_series(n), source=CMPLX(0.0, 0.0, kind=dp))
      ALLOCATE (field_series(n), source=CMPLX(0.0, 0.0, kind=dp))
      ALLOCATE (polarizability_series(n), source=CMPLX(0.0, 0.0, kind=dp))

      ! The moment ft
      ! Real part
      value_series(:) = REAL(rtbse_env%moments_trace(rtbse_env%pol_element(1))%series(:))
      CALL ft_simple(rtbse_env%time_trace%series, value_series, omega_series, ft_real_series, ft_imag_series, &
                     rtbse_env%ft_damping, rtbse_env%ft_start, .TRUE.)
      moment_series(:) = CMPLX(ft_real_series(:), ft_imag_series(:), kind=dp)
      ! Imaginary part
      value_series(:) = AIMAG(rtbse_env%moments_trace(rtbse_env%pol_element(1))%series(:))
      CALL ft_simple(rtbse_env%time_trace%series, value_series, omega_series, ft_real_series, ft_imag_series, &
                     rtbse_env%ft_damping, rtbse_env%ft_start, .TRUE.)
      moment_series(:) = moment_series(:) + CMPLX(-ft_imag_series(:), ft_real_series(:), kind=dp)

      ! Calculate the field transform - store it in ft_real_series
      IF (rtbse_env%dft_control%rtp_control%apply_delta_pulse) THEN
         ! Only divide by constant magnitude in atomic units
         field_series(:) = CMPLX(rtbse_env%dft_control%rtp_control%delta_pulse_scale, 0.0, kind=dp)
      ELSE
         ! Calculate the transform of the field as well and divide by it
         ! The field FT is not damped - assume field is localised in time
         ! The field is strictly real
         CALL ft_simple(rtbse_env%time_trace%series, rtbse_env%field_trace(rtbse_env%pol_element(2))%series, &
                        omega_series, ft_real_series, ft_imag_series, 0.0_dp, rtbse_env%ft_start, .FALSE.)
         ! Regularization for the case when ft_series becomes identically zero - TODO : Set in config
         field_series(:) = CMPLX(ft_real_series(:), ft_imag_series(:) + 1.0e-20, kind=dp)
      END IF
      ! Divide to get the polarizability series
      ! Real part
      polarizability_series(:) = moment_series(:)/field_series(:)

      ! Change units to eV for energy
      ! use value_series for energy and moment_series for polarizability
      value_series(:) = omega_series(:)*evolt
      ! Use below for printing of field FT in case of problems
      ! PRINT *, "=== Field FT"
      ! DO i=1,n
      !    PRINT '(E20.8E3,E20.8E3,E20.8E3)', value_series(i), REAL(field_series(i)), AIMAG(field_series(i))
      ! END DO
      ! PRINT *, "=== Field FT"
      ! Print out the polarizability to a file
      pol_unit = cp_print_key_unit_nr(logger, rtbse_env%pol_section, extension=".dat", &
                                      file_form="FORMATTED", file_position="REWIND")
      IF (pol_unit > 0) THEN
         IF (pol_unit == rtbse_env%unit_nr) THEN
            WRITE (pol_unit, '(A16,A20,A22,A22)') &
               " POLARIZABILITY|", "Energy [eV]", "Real polarizability", "Imag polarizability"
            DO i = 1, n
               WRITE (pol_unit, '(A16,E20.8E3,E22.8E3,E22.8E3)') &
                  " POLARIZABILITY|", value_series(i), REAL(polarizability_series(i)), AIMAG(polarizability_series(i))
            END DO
         ELSE
            WRITE (pol_unit, '(A1,A19,A20,A20,A20)') &
               "#", "omega [a.u.]", "Energy [eV]", "Real polarizability", "Imag polarizability"
            DO i = 1, n
               WRITE (pol_unit, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') &
                  omega_series(i), value_series(i), REAL(polarizability_series(i)), AIMAG(polarizability_series(i))
            END DO
            CALL cp_print_key_finished_output(pol_unit, logger, rtbse_env%ft_section)
         END IF
      END IF

      DEALLOCATE (value_series)
      DEALLOCATE (ft_real_series)
      DEALLOCATE (ft_imag_series)

      DEALLOCATE (field_series)
      DEALLOCATE (moment_series)

      DEALLOCATE (omega_series)
      DEALLOCATE (polarizability_series)
   END SUBROUTINE output_polarizability
! **************************************************************************************************
!> \brief Naively calculates the Fourier transform - it is not the bottleneck of this calculation
!> \param time_series Timestamps in atomic units of time
!> \param value_series Values to be Fourier transformed - moments, field etc. So far only real.
!> \param omega_series Array to be filled by sampled values of frequency
!> \param result_series FT of the value series - real values (cosines)
!> \param iresult_series FT of the value series - imaginary values (sines)
!> \param damping_opt Supply custom exponential damping - default is 4.0/totalTime, i.e. ratio
!>                    of last and first element in windowed value series is reduced by e^(-4)
!> \param t0_opt Carry the FT only starting from certain time - allows for exclusion of trace before
!>               the pulse application etc.
!> \author Stepan Marek
!> \date 09.2024
! **************************************************************************************************
   ! So far only for defined one dimensional series
   SUBROUTINE ft_simple(time_series, value_series, omega_series, result_series, iresult_series, &
                        damping_opt, t0_opt, subtract_initial_opt)
      REAL(kind=dp), DIMENSION(:)                        :: time_series
      REAL(kind=dp), DIMENSION(:)                        :: value_series
      REAL(kind=dp), DIMENSION(:)                        :: omega_series
      REAL(kind=dp), DIMENSION(:)                        :: result_series
      REAL(kind=dp), DIMENSION(:)                        :: iresult_series
      REAL(kind=dp), OPTIONAL                            :: damping_opt
      REAL(kind=dp), OPTIONAL                            :: t0_opt
      LOGICAL, OPTIONAL                                  :: subtract_initial_opt
      CHARACTER(len=*), PARAMETER                        :: routineN = "ft_simple"
      INTEGER                                            :: N, i, j, t0_i, j_wrap, handle
      REAL(kind=dp)                                      :: t0, delta_t, delta_omega, damping, &
                                                            value_subtract
      LOGICAL                                            :: subtract_initial

      CALL timeset(routineN, handle)

      N = SIZE(time_series)

      t0_i = 1
      IF (PRESENT(t0_opt)) THEN
         ! Determine the index at which we start applying the damping
         ! Also the index around which we fold around
         DO i = 1, N
            ! Increase until we break or reach the end of the time series
            t0_i = i
            IF (time_series(i) >= t0_opt) THEN
               EXIT
            END IF
         END DO
      END IF

      t0 = time_series(t0_i)

      ! Default damping so that at the end of the time series, divide value by e^-4
      damping = 4.0_dp/(time_series(N) - time_series(t0_i))
      IF (PRESENT(damping_opt)) THEN
         IF (damping_opt >= 0.0_dp) damping = damping_opt
      END IF

      IF (PRESENT(subtract_initial_opt)) THEN
         subtract_initial = subtract_initial_opt
      ELSE
         subtract_initial = .TRUE.
      END IF

      ! Construct the grid
      ! Assume series equidistant
      delta_t = time_series(2) - time_series(1)
      delta_omega = twopi/(time_series(N) - time_series(1))
      ! Subtract initial value, if requested (default is to subtract the value)
      value_subtract = 0.0_dp
      IF (subtract_initial) value_subtract = value_series(t0_i)
      DO i = 1, N
         result_series(i) = 0.0_dp
         iresult_series(i) = 0.0_dp
         omega_series(i) = delta_omega*(i - 1)
         DO j = 1, N
            j_wrap = MODULO(j + t0_i - 2, N) + 1
            result_series(i) = result_series(i) + COS(twopi*(i - 1)*(j - 1)/N)* &
                               EXP(-damping*delta_t*(j - 1))*(value_series(j_wrap) - value_subtract)
            iresult_series(i) = iresult_series(i) + SIN(twopi*(i - 1)*(j - 1)/N)* &
                                EXP(-damping*delta_t*(j - 1))*(value_series(j_wrap) - value_subtract)
         END DO
      END DO
      delta_omega = twopi/(time_series(N) - time_series(1))
      result_series(:) = delta_t*result_series(:)
      iresult_series(:) = delta_t*iresult_series(:)

      CALL timestop(handle)

   END SUBROUTINE
END MODULE rt_bse_io

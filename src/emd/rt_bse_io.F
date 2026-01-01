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
   USE physcon, ONLY: femtoseconds
   USE rt_propagation_output, ONLY: print_moments

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_bse_io"

   #:include "rt_bse_macros.fypp"

   PUBLIC :: output_moments, &
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
         CASE (rtp_bse_ham_g0w0)
            WRITE (rtbse_env%unit_nr, '(A58)') 'G0W0'
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
!> \brief Writes the header for the etrs iteration updates - only for log level > low
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
         rtbse_env%field_trace(i, n) = CMPLX(rtbse_env%field(i), 0.0, kind=dp)
      END DO
      rtbse_env%time_trace(n) = rtbse_env%sim_time
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

      n = rtbse_env%sim_step - rtbse_env%sim_start_orig + 1

      DO j = 1, rtbse_env%n_spin
         ! Need to transpose due to the definition of trace function
         CALL cp_cfm_to_fm(msource=rho(j), mtargetr=rtbse_env%real_workspace(2))
         DO i = 1, 3
            ! Moments should be symmetric, test without transopose?
            CALL cp_fm_transpose(rtbse_env%moments(i), rtbse_env%real_workspace(1))
            CALL cp_fm_trace(rtbse_env%real_workspace(1), rtbse_env%real_workspace(2), moments_re(i))
            ! Scale by spin degeneracy and electron charge
            moments_re(i) = -moments_re(i)*rtbse_env%spin_degeneracy
            rtbse_env%moments_trace(j, i, n) = CMPLX(moments_re(i), 0.0, kind=dp)
         END DO
         ! Same for imaginary part
         CALL cp_cfm_to_fm(msource=rho(j), mtargeti=rtbse_env%real_workspace(2))
         DO i = 1, 3
            CALL cp_fm_transpose(rtbse_env%moments(i), rtbse_env%real_workspace(1))
            CALL cp_fm_trace(rtbse_env%real_workspace(1), rtbse_env%real_workspace(2), moments_re(i))
            ! Scale by spin degeneracy and electron charge
            moments_re(i) = -moments_re(i)*rtbse_env%spin_degeneracy
            rtbse_env%moments_trace(j, i, n) = rtbse_env%moments_trace(j, i, n) + CMPLX(0.0, moments_re(i), kind=dp)
         END DO
      END DO
      ! Output to the file
      CALL print_moments(rtbse_env%moments_section, rtbse_env%unit_nr, rtbse_env%moments_trace(:, :, n), &
                         rtbse_env%sim_time, .TRUE.)
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
END MODULE rt_bse_io

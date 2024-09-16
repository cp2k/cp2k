!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Input/output from the real time propagation via GWBSE method.
!> \author Stepan Marek (08.24)
! **************************************************************************************************

MODULE rt_gwbse_io
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
                                 cp_print_key_generate_filename
   USE input_section_types, ONLY: section_vals_type
   USE rt_gwbse_types, ONLY: gwbse_env_type, &
                             multiply_cfm_fm, &
                             multiply_fm_cfm
   USE cp_files, ONLY: open_file, &
                       file_exists, &
                       close_file
   USE physcon, ONLY: femtoseconds
   USE mathconstants, ONLY: twopi

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_gwbse_io"

   #:include "rt_gwbse_macros.fypp"

   PUBLIC :: output_moments, &
             read_moments, &
             output_moments_ft, &
             output_polarizability, &
             output_field, &
             read_field, &
             output_mos_contravariant, &
             output_mos_covariant, &
             output_restart, &
             read_restart

CONTAINS

! **************************************************************************************************
!> \brief Outputs the matrix in MO basis for matrix coefficients corresponding to contravariant
!>        operator, i.e. density matrix
!> \param gwbse_env Entry point - gwbse environment
!> \param rho Density matrix in AO basis
!> \param mos Spin dependent molecular orbitals
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE output_mos_contravariant(gwbse_env, rho, mos, print_key_section)
      TYPE(gwbse_env_type)                               :: gwbse_env
      TYPE(cp_cfm_p_type), DIMENSION(:)                  :: rho
      TYPE(cp_fm_type), DIMENSION(:)                     :: mos
      TYPE(section_vals_type), POINTER                   :: print_key_section
      CHARACTER(len=*), PARAMETER                        :: routineN = "output_mos_rho"
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: j, rho_unit_re, rho_unit_im
      CHARACTER(len=14), DIMENSION(4)                    :: file_labels

      file_labels(1) = "_SPIN_A_RE.dat"
      file_labels(2) = "_SPIN_A_IM.dat"
      file_labels(3) = "_SPIN_B_RE.dat"
      file_labels(4) = "_SPIN_B_IM.dat"
      logger => cp_get_default_logger()
      ! Start by multiplying the current density by MOS
      DO j = 1, gwbse_env%n_spin
         rho_unit_re = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(2*j - 1))
         rho_unit_im = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(2*j))
         ! Transform the density matrix into molecular orbitals basis and print it out
         ! S * rho
         CALL multiply_fm_cfm("N", "N", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao, &
                              1.0_dp, gwbse_env%S_fm%matrix, rho(j)%matrix, &
                              0.0_dp, gwbse_env%rho_workspace(1)%matrix)
         ! C^T * S * rho
         CALL multiply_fm_cfm("T", "N", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao, &
                              1.0_dp, mos(j), gwbse_env%rho_workspace(1)%matrix, &
                              0.0_dp, gwbse_env%rho_workspace(2)%matrix)
         ! C^T * S * rho * S
         CALL multiply_cfm_fm("N", "N", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao, &
                              1.0_dp, gwbse_env%rho_workspace(2)%matrix, gwbse_env%S_fm%matrix, &
                              0.0_dp, gwbse_env%rho_workspace(1)%matrix)
         ! C^T * S * rho * S * C
         CALL multiply_cfm_fm("N", "N", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao, &
                              1.0_dp, gwbse_env%rho_workspace(1)%matrix, mos(j), &
                              0.0_dp, gwbse_env%rho_workspace(2)%matrix)
         ! Print real and imaginary parts separately
         CALL cp_cfm_to_fm(gwbse_env%rho_workspace(2)%matrix, &
                           gwbse_env%real_workspace(1)%matrix, gwbse_env%real_workspace(2)%matrix)
         CALL cp_fm_write_formatted(gwbse_env%real_workspace(1)%matrix, rho_unit_re)
         CALL cp_fm_write_formatted(gwbse_env%real_workspace(2)%matrix, rho_unit_im)
         CALL cp_print_key_finished_output(rho_unit_re, logger, print_key_section)
         CALL cp_print_key_finished_output(rho_unit_im, logger, print_key_section)
      END DO
   END SUBROUTINE
! **************************************************************************************************
!> \brief Outputs the matrix in MO basis for matrix components corresponding to covariant representation,
!>        i.e. the Hamiltonian matrix
!> \param gwbse_env Entry point - gwbse environment
!> \param cohsex cohsex matrix in AO basis, covariant representation
!> \param mos Spin dependent molecular orbitals
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE output_mos_covariant(gwbse_env, ham, mos, print_key_section)
      TYPE(gwbse_env_type)                               :: gwbse_env
      TYPE(cp_cfm_p_type), DIMENSION(:)                  :: ham
      TYPE(cp_fm_type), DIMENSION(:)                     :: mos
      TYPE(section_vals_type), POINTER                   :: print_key_section
      CHARACTER(len=*), PARAMETER                        :: routineN = "output_mos_cohsex"
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: j, rho_unit_re, rho_unit_im
      CHARACTER(len=21), DIMENSION(4)                    :: file_labels

      file_labels(1) = "_SPIN_A_RE.dat"
      file_labels(2) = "_SPIN_A_IM.dat"
      file_labels(3) = "_SPIN_B_RE.dat"
      file_labels(4) = "_SPIN_B_IM.dat"
      logger => cp_get_default_logger()
      DO j = 1, gwbse_env%n_spin
         rho_unit_re = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(2*j - 1))
         rho_unit_im = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(2*j))
         ! C^T * cohsex
         CALL multiply_fm_cfm("T", "N", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao, &
                              1.0_dp, mos(j), ham(j)%matrix, &
                              0.0_dp, gwbse_env%rho_workspace(1)%matrix)
         ! C^T * cohsex * C
         CALL multiply_cfm_fm("N", "N", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao, &
                              1.0_dp, gwbse_env%rho_workspace(1)%matrix, mos(j), &
                              0.0_dp, gwbse_env%rho_workspace(2)%matrix)
         ! Print real and imaginary parts separately
         CALL cp_cfm_to_fm(gwbse_env%rho_workspace(2)%matrix, &
                           gwbse_env%real_workspace(1)%matrix, gwbse_env%real_workspace(2)%matrix)
         CALL cp_fm_write_formatted(gwbse_env%real_workspace(1)%matrix, rho_unit_re)
         CALL cp_fm_write_formatted(gwbse_env%real_workspace(2)%matrix, rho_unit_im)
         CALL cp_print_key_finished_output(rho_unit_re, logger, print_key_section)
         CALL cp_print_key_finished_output(rho_unit_im, logger, print_key_section)
      END DO
   END SUBROUTINE
! **************************************************************************************************
!> \brief Prints the current field components into a file provided by input
!> \param gwbse_env Entry point - gwbse environment
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE output_field(gwbse_env)
      TYPE(gwbse_env_type)                               :: gwbse_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "output_field"
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: field_unit, n, i

      ! Get logger
      logger => cp_get_default_logger()
      ! Get file descriptor
      field_unit = cp_print_key_unit_nr(logger, gwbse_env%field_section, extension=".dat")
      ! If the file descriptor is non-zero, output field
      ! TODO : Output also in SI
      IF (field_unit /= -1) THEN
         WRITE (field_unit, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') gwbse_env%sim_time*femtoseconds, &
            gwbse_env%field(1), gwbse_env%field(2), gwbse_env%field(3)
      END IF
      ! Write the output to memory for FT
      ! Need the absolute index
      n = gwbse_env%sim_step - gwbse_env%sim_start_orig + 1
      DO i=1, 3
         gwbse_env%field_trace(i)%series(n) = gwbse_env%field(i)
      END DO
      gwbse_env%time_trace%series(n) = gwbse_env%sim_time
      CALL cp_print_key_finished_output(field_unit, logger, gwbse_env%field_section)

   END SUBROUTINE
! **************************************************************************************************
!> \brief Reads the field from the files provided by input - useful for the continuation run
!> \param gwbse_env Entry point - gwbse environment
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE read_field(gwbse_env)
      TYPE(gwbse_env_type)                               :: gwbse_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "read_field"
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(len=default_path_length)                 :: save_name
      INTEGER                                            :: k, n, field_unit

      ! Get logger
      logger => cp_get_default_logger()
      ! Get file name
      save_name = cp_print_key_generate_filename(logger, gwbse_env%field_section, extension=".dat", my_local=.FALSE.)
      IF (file_exists(save_name)) THEN
         CALL open_file(save_name, file_status="OLD", file_form="FORMATTED", file_action="READ", &
                        unit_number=field_unit)
         DO k=gwbse_env%sim_start_orig,gwbse_env%sim_start
            n = k - gwbse_env%sim_start_orig + 1
            WRITE (field_unit, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') gwbse_env%time_trace%series(n), &
               gwbse_env%field_trace(1)%series(n), gwbse_env%field_trace(2)%series(n), gwbse_env%field_trace(3)%series(n)
            ! Set the time units back to atomic units
            gwbse_env%time_trace%series(n) = gwbse_env%time_trace%series(n) / femtoseconds
         END DO
         CALL close_file(field_unit)
      END IF
   END SUBROUTINE read_field

! **************************************************************************************************
!> \brief Outputs the expectation value of moments from a given density matrix
!> \note  Moments matrix is provided by the gwbse_env, uses rho_workspace(1:3)
!> \param gwbse_env Entry point - gwbse environment
!> \param rho Density matrix in AO basis
!> \param rtp_section RTP section of the input parameters, where moments destination may be present
! **************************************************************************************************
   SUBROUTINE output_moments(gwbse_env, rho)
      TYPE(gwbse_env_type)                               :: gwbse_env
      TYPE(cp_cfm_p_type), DIMENSION(:)                  :: rho
      CHARACTER(len=*), PARAMETER                        :: routineN = "output_moments"
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
      DO j = 1, gwbse_env%n_spin
         moments_unit_re = cp_print_key_unit_nr(logger, gwbse_env%moments_section, extension=file_labels(2*j - 1))
         moments_unit_im = cp_print_key_unit_nr(logger, gwbse_env%moments_section, extension=file_labels(2*j))
         ! If, for any reason, the file unit is not provided, skip to next cycle immediately
         ! TODO : Either centralize or let user define the output format
         ! TODO : Let user define the output units
         ! TODO : Handle case of complex molecular orbital matrix?
         ! TODO : Specify output units in config
         ! Need to transpose due to the definition of trace function
         CALL cp_cfm_to_fm(msource=rho(j)%matrix, mtargetr=gwbse_env%real_workspace(2)%matrix)
         DO i = 1, 3
            ! Moments should be symmetric, test without transopose?
            CALL cp_fm_transpose(gwbse_env%moments(i)%matrix, gwbse_env%real_workspace(1)%matrix)
            CALL cp_fm_trace(gwbse_env%real_workspace(1)%matrix, gwbse_env%real_workspace(2)%matrix, moments(i))
            ! Scale by spin degeneracy and electron charge
            moments(i) = -moments(i)*gwbse_env%spin_degeneracy
         END DO
         ! Output to the file
         IF (moments_unit_re > 0) WRITE (moments_unit_re, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') gwbse_env%sim_time*femtoseconds, &
            moments(1), moments(2), moments(3)
         ! Save to memory for FT - real part
         n = gwbse_env%sim_step - gwbse_env%sim_start_orig + 1
         DO i = 1, 3
            gwbse_env%moments_trace(i)%series(n) = CMPLX(moments(i), 0.0, kind=dp)
         END DO
         CALL cp_print_key_finished_output(moments_unit_re, logger, gwbse_env%moments_section)
         ! Same for imaginary part
         CALL cp_cfm_to_fm(msource=rho(j)%matrix, mtargeti=gwbse_env%real_workspace(2)%matrix)
         DO i = 1, 3
            CALL cp_fm_transpose(gwbse_env%moments(i)%matrix, gwbse_env%real_workspace(1)%matrix)
            CALL cp_fm_trace(gwbse_env%real_workspace(1)%matrix, gwbse_env%real_workspace(2)%matrix, moments(i))
            ! Scale by spin degeneracy and electron charge
            moments(i) = -moments(i)*gwbse_env%spin_degeneracy
         END DO
         ! Output to the file
         IF (moments_unit_im > 0) WRITE (moments_unit_im, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') gwbse_env%sim_time*femtoseconds, &
            moments(1), moments(2), moments(3)
         ! Save to memory for FT - imag part
         DO i = 1, 3
            gwbse_env%moments_trace(i)%series(n) = gwbse_env%moments_trace(i)%series(n) + CMPLX(0.0, moments(i), kind=dp)
         END DO
         ! Close the files
         CALL cp_print_key_finished_output(moments_unit_im, logger, gwbse_env%moments_section)
      END DO
   END SUBROUTINE
! **************************************************************************************************
!> \brief Outputs the restart info (last finished iteration step) + restard density matrix
!> \param restart_section Print key section for the restart files
!> \param rho Density matrix in AO basis
!> \param time_index Time index to be written into the info file
! **************************************************************************************************
   SUBROUTINE output_restart(gwbse_env, rho, time_index)
      TYPE(gwbse_env_type), POINTER                      :: gwbse_env
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER         :: rho
      INTEGER                                            :: time_index
      CHARACTER(len=*), PARAMETER                        :: routineN = "output_restart"
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: workspace
      CHARACTER(len=17), DIMENSION(4)                    :: file_labels
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: rho_unit_nr, i

      ! Default labels distinguishing up to two spin species and real/imaginary parts
      file_labels(1) = "_SPIN_A_RE.matrix"
      file_labels(2) = "_SPIN_A_IM.matrix"
      file_labels(3) = "_SPIN_B_RE.matrix"
      file_labels(4) = "_SPIN_B_IM.matrix"

      logger => cp_get_default_logger()

      workspace => gwbse_env%real_workspace

      DO i = 1, gwbse_env%n_spin
         CALL cp_cfm_to_fm(rho(i)%matrix, workspace(1)%matrix, workspace(2)%matrix)
         ! Real part
         rho_unit_nr = cp_print_key_unit_nr(logger, gwbse_env%restart_section, extension=file_labels(2*i - 1), &
                                            file_form="UNFORMATTED", file_position="REWIND")
         CALL cp_fm_write_unformatted(workspace(1)%matrix, rho_unit_nr)
         CALL cp_print_key_finished_output(rho_unit_nr, logger, gwbse_env%restart_section)
         ! Imag part
         rho_unit_nr = cp_print_key_unit_nr(logger, gwbse_env%restart_section, extension=file_labels(2*i), &
                                            file_form="UNFORMATTED", file_position="REWIND")
         CALL cp_fm_write_unformatted(workspace(2)%matrix, rho_unit_nr)
         CALL cp_print_key_finished_output(rho_unit_nr, logger, gwbse_env%restart_section)
         ! Info
         rho_unit_nr = cp_print_key_unit_nr(logger, gwbse_env%restart_section, extension=".info", &
                                            file_form="UNFORMATTED", file_position="REWIND")
         IF (rho_unit_nr > 0) WRITE (rho_unit_nr) time_index
         CALL cp_print_key_finished_output(rho_unit_nr, logger, gwbse_env%restart_section)
      END DO
   END SUBROUTINE output_restart
! **************************************************************************************************
!> \brief Reads the density matrix from restart files and updates the starting time
!> \param restart_section Print key section for the restart files
!> \param rho Density matrix in AO basis
!> \param time_index Time index to be written into the info file
! **************************************************************************************************
   SUBROUTINE read_restart(gwbse_env)
      TYPE(gwbse_env_type), POINTER                      :: gwbse_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "read_restart"
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(len=default_path_length)                 :: save_name, save_name_2
      INTEGER                                            :: rho_unit_nr, j
      CHARACTER(len=17), DIMENSION(4)                    :: file_labels

      ! This allows the delta kick and output of moment at time 0 in all cases
      ! except the case when both imaginary and real parts of the density are read
      gwbse_env%restart_extracted = .FALSE.
      logger => cp_get_default_logger()
      ! Start by probing/loading info file
      save_name = cp_print_key_generate_filename(logger, gwbse_env%restart_section, extension=".info", my_local=.FALSE.)
      IF (file_exists(save_name)) THEN
         CALL open_file(save_name, file_status="OLD", file_form="UNFORMATTED", file_action="READ", &
                        unit_number=rho_unit_nr)
         ! TODO : Some error catching and handling
         READ (rho_unit_nr) gwbse_env%sim_start
         CALL close_file(rho_unit_nr)
         IF (gwbse_env%unit_nr > 0) WRITE (gwbse_env%unit_nr, *) "Starting from timestep ", gwbse_env%sim_start, &
            ", delta kick NOT applied"
      ELSE
         CPWARN("Restart required but no info file found - starting from sim_step given in input")
      END IF

      ! Default labels distinguishing up to two spin species and real/imaginary parts
      file_labels(1) = "_SPIN_A_RE.matrix"
      file_labels(2) = "_SPIN_A_IM.matrix"
      file_labels(3) = "_SPIN_B_RE.matrix"
      file_labels(4) = "_SPIN_B_IM.matrix"
      DO j = 1, gwbse_env%n_spin
       save_name = cp_print_key_generate_filename(logger, gwbse_env%restart_section, extension=file_labels(2*j - 1), my_local=.FALSE.)
       save_name_2 = cp_print_key_generate_filename(logger, gwbse_env%restart_section, extension=file_labels(2*j), my_local=.FALSE.)
       IF (file_exists(save_name) .AND. file_exists(save_name_2)) THEN
            CALL open_file(save_name, file_status="OLD", file_form="UNFORMATTED", file_action="READ", &
                           unit_number=rho_unit_nr)
            CALL cp_fm_read_unformatted(gwbse_env%real_workspace(1)%matrix, rho_unit_nr)
            CALL close_file(rho_unit_nr)
            CALL open_file(save_name_2, file_status="OLD", file_form="UNFORMATTED", file_action="READ", &
                           unit_number=rho_unit_nr)
            CALL cp_fm_read_unformatted(gwbse_env%real_workspace(2)%matrix, rho_unit_nr)
            CALL close_file(rho_unit_nr)
            CALL cp_fm_to_cfm(gwbse_env%real_workspace(1)%matrix, gwbse_env%real_workspace(2)%matrix, &
                              gwbse_env%rho(j)%matrix)
            gwbse_env%restart_extracted = .TRUE.
         ELSE
            ! TODO : Spin dependent warning - for some spins, SCF may be used while for others restart is used
            CPWARN("Restart required but some restart density matrices not found - starting from SCF density.")
         END IF
      END DO
   END SUBROUTINE read_restart
! **************************************************************************************************
!> \brief Reads the moments and time traces from the save files
!> \param gwbse_env GW-BSE environment (assumes consistent setup, i.e. a continuation run).
!>                  Otherwise, the traces are set at zero
! **************************************************************************************************
   SUBROUTINE read_moments(gwbse_env)
      TYPE(gwbse_env_type), POINTER                      :: gwbse_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "read_moments"
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
      DO j=1,gwbse_env%n_spin
         save_name = cp_print_key_generate_filename(logger, gwbse_env%moments_section, extension=file_labels(2*j - 1), my_local=.FALSE.)
         save_name_2 = cp_print_key_generate_filename(logger, gwbse_env%moments_section, extension=file_labels(2*j), my_local=.FALSE.)
         IF (file_exists(save_name) .AND. file_exists(save_name_2)) THEN
            CALL open_file(save_name, file_status="OLD", file_form="FORMATTED", file_action="READ", &
                           unit_number=moments_unit_re)
            CALL open_file(save_name_2, file_status="OLD", file_form="FORMATTED", file_action="READ", &
                           unit_number=moments_unit_im)
            IF (gwbse_env%unit_nr > 0) WRITE(gwbse_env%unit_nr, *) "Expecting following number of time steps : ", &
                    gwbse_env%sim_start - gwbse_env%sim_start_orig + 1
            ! Extra time step for the initial one
            DO k=gwbse_env%sim_start_orig,gwbse_env%sim_start
               ! Determine the absolute time step - offset in memory
               n = k - gwbse_env%sim_start_orig + 1
               READ (moments_unit_re, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') timestamp, &
                     moments_re(1), moments_re(2), moments_re(3)
               READ (moments_unit_im, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') timestamp, &
                     moments_im(1), moments_im(2), moments_im(3)
               DO i = 1, 3
                  gwbse_env%moments_trace(i)%series(n) =CMPLX(moments_re(i), moments_im(i), kind=dp)
               END DO
            END DO
            ! Change back to atomic units in the trace
            gwbse_env%time_trace%series(:) = gwbse_env%time_trace%series(:) / femtoseconds
            CALL close_file(moments_unit_re)
            CALL close_file(moments_unit_im)
         END IF
      END DO
   END SUBROUTINE read_moments
! **************************************************************************************************
!> \brief Outputs the Fourier transform of moments stored in the environment memory to the configured file
!> \param gwbse_env GW-BSE environment
! **************************************************************************************************
   SUBROUTINE output_moments_ft(gwbse_env)
      TYPE(gwbse_env_type), POINTER                      :: gwbse_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "output_moments_ft"
      TYPE(cp_logger_type), POINTER                      :: logger
      REAL(kind=dp), DIMENSION(:), POINTER               :: omega_series, &
                                                            ft_real_series, &
                                                            ft_imag_series, &
                                                            value_series
      ! Stores the data in ready for output format
      !  - first dimension is 6 - 1 - real part along x, 2 - imag part along x, 3 - real part along y, ...
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE         :: ft_full_series
      INTEGER                                            :: i, n, ft_unit

      logger => cp_get_default_logger()

      n = gwbse_env%sim_nsteps+2
      NULLIFY (omega_series)
      ALLOCATE (omega_series(n), source=0.0_dp)
      NULLIFY (ft_real_series)
      ALLOCATE (ft_real_series(n), source=0.0_dp)
      NULLIFY (ft_imag_series)
      ALLOCATE (ft_imag_series(n), source=0.0_dp)
      NULLIFY (value_series)
      ALLOCATE (value_series(n), source=0.0_dp)
      ALLOCATE (ft_full_series(6,n))
      ! Carry out for each direction independently and real and imaginary parts also independently
      ! TODO : Finite pulse
      DO i = 1, 3
         ! Real part of the value first 
         value_series(:) = REAL(gwbse_env%moments_trace(i)%series(:))
         CALL ft_slow(gwbse_env%time_trace%series, value_series, omega_series, ft_real_series, ft_imag_series, &
                 gwbse_env%ft_damping, gwbse_env%ft_start)
         ft_full_series(2*i-1, :) = ft_real_series(:)
         ft_full_series(2*i, :) = ft_imag_series(:)
         ! Now imaginary part
         value_series(:) = AIMAG(gwbse_env%moments_trace(i)%series(:))
         CALL ft_slow(gwbse_env%time_trace%series, value_series, omega_series, ft_real_series, ft_imag_series, &
                 gwbse_env%ft_damping, gwbse_env%ft_start)
         ft_full_series(2*i-1, :) = ft_full_series(2*i-1,:) - ft_imag_series
         ft_full_series(2*i, :) = ft_full_series(2*i, :) + ft_real_series
      END DO
      DEALLOCATE (ft_real_series)
      DEALLOCATE (ft_imag_series)
      DEALLOCATE (value_series)
      ! Now, write these to file
      ft_unit = cp_print_key_unit_nr(logger, gwbse_env%ft_section, extension=".dat", &
                                         file_form="FORMATTED", file_position="REWIND")
      DO i=1, n
         WRITE (ft_unit, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3,E20.8E3,E20.8E3,E20.8E3)') &
                 omega_series(i), ft_full_series(1,i), ft_full_series(2,i), &
                 ft_full_series(3,i), ft_full_series(4,i), &
                 ft_full_series(5,i), ft_full_series(6,i)
      END DO
      CALL cp_print_key_finished_output(ft_unit, logger, gwbse_env%ft_section)
      DEALLOCATE (omega_series)
      DEALLOCATE (ft_full_series)
   END SUBROUTINE output_moments_ft
! **************************************************************************************************
!> \brief Outputs the isotropic polarizability tensor element alpha _ ij = mu_i(omega)/E_j(omega),
!>        where i and j are provided by the configuration. The tensor element is energy dependent and
!>        has real and imaginary parts
!> \param gwbse_env GW-BSE environment
! **************************************************************************************************
   SUBROUTINE output_polarizability(gwbse_env)
      TYPE(gwbse_env_type), POINTER                      :: gwbse_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "output_moments_ft"
      TYPE(cp_logger_type), POINTER                      :: logger
      REAL(kind=dp), DIMENSION(:), POINTER               :: omega_series, &
                                                            ft_real_series, &
                                                            ft_imag_series, &
                                                            value_series
      COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE        :: moment_series, &
                                                            field_series, &
                                                            polarizability_series
      INTEGER                                            :: pol_unit, &
                                                            i, n

      logger => cp_get_default_logger()

      n = gwbse_env%sim_nsteps+2
      ! All allocations together, although could save some memory, if required by consequent deallocations
      NULLIFY (omega_series)
      ALLOCATE (omega_series(n), source=0.0_dp)
      NULLIFY (ft_real_series)
      ALLOCATE (ft_real_series(n), source=0.0_dp)
      NULLIFY (ft_imag_series)
      ALLOCATE (ft_imag_series(n), source=0.0_dp)
      NULLIFY (value_series)
      ALLOCATE (value_series(n), source=0.0_dp)
      ALLOCATE (moment_series(n), source=CMPLX(0.0, 0.0, kind=dp))
      ALLOCATE (field_series(n), source=CMPLX(0.0, 0.0, kind=dp))
      ALLOCATE (polarizability_series(n), source=CMPLX(0.0, 0.0, kind=dp))

      ! The moment ft
      ! Real part
      value_series(:) = REAL(gwbse_env%moments_trace(gwbse_env%pol_element(1))%series(:))
      CALL ft_slow(gwbse_env%time_trace%series, value_series, omega_series, ft_real_series, ft_imag_series, &
                   gwbse_env%ft_damping, gwbse_env%ft_start)
      moment_series(:) = CMPLX(ft_real_series(:), ft_imag_series(:), kind=dp)
      ! Imaginary part
      value_series(:) = AIMAG(gwbse_env%moments_trace(gwbse_env%pol_element(1))%series(:))
      CALL ft_slow(gwbse_env%time_trace%series, value_series, omega_series, ft_real_series, ft_imag_series, &
                   gwbse_env%ft_damping, gwbse_env%ft_start)
      moment_series(:) = moment_series(:) + CMPLX(-ft_imag_series(:), ft_real_series(:), kind=dp)

      ! Calculate the field transform - store it in ft_real_series
      IF (gwbse_env%dft_control%rtp_control%apply_delta_pulse) THEN
         ! Only divide by constant magnitude in atomic units
         field_series(:) = CMPLX(gwbse_env%dft_control%rtp_control%delta_pulse_scale, 0.0, kind=dp)
      ELSE
         ! Calculate the transform of the field as well and divide by it
         ! The field is strictly real
         CALL ft_slow(gwbse_env%time_trace%series, gwbse_env%field_trace(gwbse_env%pol_element(2))%series, omega_series, ft_real_series, ft_imag_series, &
                      gwbse_env%ft_damping, gwbse_env%ft_start)
         field_series(:) = CMPLX(ft_real_series(:), ft_imag_series(:), kind=dp)
      END IF
      ! Divide to get the polarizability series
      ! Real part
      polarizability_series(:) = moment_series(:) / field_series(:)

      ! Change units to eV for energy
      ! use value_series for energy and moment_series for polarizability
      value_series(:) = omega_series(:) * 27.211386
      ! TODO : Use the physcon for conversion
      ! Print out the polarizability to a file
      ! TODO : Change to custom file
      pol_unit = cp_print_key_unit_nr(logger, gwbse_env%pol_section, extension=".dat", &
                                         file_form="FORMATTED", file_position="REWIND")
      WRITE (pol_unit, '(A1,A19,A20,A20,A20)') "#", "omega [a.u.]", "Energy [eV]", "Real polarizability", "Imag polarizability"
      DO i=1, n
         WRITE (pol_unit, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') &
                 omega_series(i), value_series(i), REAL(polarizability_series(i)), AIMAG(polarizability_series(i))
      END DO
      CALL cp_print_key_finished_output(pol_unit, logger, gwbse_env%ft_section)

      DEALLOCATE (value_series)

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
   ! TODO : Maybe simple instead of slow
   SUBROUTINE ft_slow(time_series, value_series, omega_series, result_series, iresult_series, damping_opt, t0_opt)
      REAL(kind=dp), DIMENSION(:), POINTER               :: time_series
      REAL(kind=dp), DIMENSION(:), POINTER               :: value_series
      REAL(kind=dp), DIMENSION(:), POINTER               :: omega_series
      REAL(kind=dp), DIMENSION(:), POINTER               :: result_series
      REAL(kind=dp), DIMENSION(:), POINTER               :: iresult_series
      REAL(kind=dp), OPTIONAL                            :: damping_opt
      REAL(kind=dp), OPTIONAL                            :: t0_opt
      CHARACTER(len=*), PARAMETER                        :: routineN="ft_slow"
      INTEGER                                            :: M, N, i, j, t0_i
      REAL(kind=dp)                                      :: t0, delta_t, delta_omega, damping

      N = SIZE(time_series)
      ! TODO : Some sanity checks?

      t0_i = 1
      IF (PRESENT(t0_opt)) THEN
         ! Determine the index at which we start applying the damping
         DO i=1,N
            ! Increase until we break or reach the end of the time series
            t0_i = i
            IF (time_series(i) >= t0_opt) THEN
               EXIT
            END IF
         END DO
      END IF

      t0 = time_series(t0_i)

      ! Default damping so that at the end of the time series, divide value by e^-4
      damping = 4.0_dp / (time_series(N) - time_series(t0_i))
      ! damping = 0.02_dp
      IF (PRESENT(damping_opt)) THEN
         IF (damping_opt >= 0.0_dp) damping = damping_opt
      END IF

      ! Construct the grid
      ! Number of points that we actually integrate over
      M = N - t0_i + 1
      ! Time differential
      delta_t = (time_series(N) - time_series(t0_i)) / REAL(M, kind=dp)
      ! Frequency differential
      delta_omega = twopi / (time_series(N) - time_series(t0_i))
      ! TODO : Some parallelism/vectorisation? But this will never be a bottle neck of the program
      ! TODO : Maybe change and allow for different sizes in the timestep? Or is that not a relevant case?
      DO i = 1, M
         ! Initialise the integral
         result_series(i) = 0.0_dp
         iresult_series(i) = 0.0_dp
         omega_series(i) = delta_omega * (i-1)
         DO j = 1, M
            result_series(i) = result_series(i) + COS(twopi * (i-1) * (j-1)/ M) * EXP(- damping * delta_t * (j-1)) * value_series(j + t0_i - 1)
            iresult_series(i) = iresult_series(i) + SIN(twopi * (i-1) * (j-1)/ M) * EXP(- damping * delta_t * (j-1)) * value_series(j + t0_i - 1)
         END DO
      END DO
      result_series(:) = delta_t * result_series(:)
      iresult_series(:) = delta_t * iresult_series(:)

   END SUBROUTINE
END MODULE rt_gwbse_io

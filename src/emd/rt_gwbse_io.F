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

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_gwbse_io"

   #:include "rt_gwbse_macros.fypp"

   PUBLIC :: output_moments, &
             output_field, &
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
         ! TODO : Number of molecular orbitals and atomic orbitals might not be the same
         ! TODO : Spins in S_fm
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
         ! TODO : Number of molecular orbitals and atomic orbitals might not be the same
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
      INTEGER                                            :: field_unit

      ! Get logger
      logger => cp_get_default_logger()
      ! Get file descriptor
      field_unit = cp_print_key_unit_nr(logger, gwbse_env%field_section, extension=".dat")
      ! If the file descriptor is non-zero, output field
      IF (field_unit /= -1) THEN
         WRITE (field_unit, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') gwbse_env%sim_time*femtoseconds, &
            gwbse_env%field(1), gwbse_env%field(2), gwbse_env%field(3)
      END IF
      CALL cp_print_key_finished_output(field_unit, logger, gwbse_env%field_section)

   END SUBROUTINE

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
      INTEGER                                            :: i, j, moments_unit_re, moments_unit_im
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
         ! Same for imaginary part
         CALL cp_cfm_to_fm(msource=rho(j)%matrix, mtargeti=gwbse_env%real_workspace(2)%matrix)
         DO i = 1, 3
            CALL cp_fm_transpose(gwbse_env%moments(i)%matrix, gwbse_env%real_workspace(1)%matrix)
            CALL cp_fm_trace(gwbse_env%real_workspace(1)%matrix, gwbse_env%real_workspace(2)%matrix, moments(i))
            ! Scale by spin degeneracy and electron charge
            moments(i) = -moments(i)*gwbse_env%spin_degeneracy
         END DO
         ! Output to the file
        IF (moments_unit_im > 0) WRITE (moments_unit_im, '(D16.8,D16.8,D16.8,D16.8)') gwbse_env%sim_time*femtoseconds, moments(1), &
            moments(2), moments(3)
         ! Close the files
         CALL cp_print_key_finished_output(moments_unit_re, logger, gwbse_env%moments_section)
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
END MODULE rt_gwbse_io

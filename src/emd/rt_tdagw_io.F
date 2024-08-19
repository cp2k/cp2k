!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Input/output from the real time propagation via TD-aGW method.
!> \author Stepan Marek (08.24)
! **************************************************************************************************

MODULE rt_tdagw_io
   USE cp_fm_types,                     ONLY: cp_fm_type,&
                                              cp_fm_p_type,&
                                              cp_fm_read_unformatted,&
                                              cp_fm_write_unformatted,&
                                              cp_fm_write_formatted
   USE kinds,                           ONLY: dp,&
                                              default_path_length
   USE cp_fm_basic_linalg,              ONLY: cp_fm_trace,&
                                              cp_fm_transpose,&
                                              cp_fm_norm
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE cp_log_handling,                 ONLY: cp_logger_type,&
                                              cp_get_default_logger
   USE cp_output_handling,              ONLY: cp_print_key_unit_nr,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_generate_filename
   USE input_section_types,             ONLY: section_vals_type
   USE rt_tdagw_types,                  ONLY: tdagw_env_type
   USE cp_files,                        ONLY: open_file,&
                                              file_exists,&
                                              close_file
   USE physcon,                         ONLY: femtoseconds

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_tdagw_io"

   #:include "tdagw_macros.fypp"

   PUBLIC :: output_moments,&
             output_field,&
             output_mos_contravariant,&
             output_mos_covariant,&
             output_restart,&
             read_restart

CONTAINS

! **************************************************************************************************
!> \brief Outputs the matrix in MO basis for matrix coefficients corresponding to contravariant
!>        operator, i.e. density matrix
!> \param tdagw_env Entry point - tdagw environment
!> \param rho Density matrix in AO basis
!> \param mos Spin dependent molecular orbitals
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE output_mos_contravariant(tdagw_env, rho, mos, print_key_section)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: rho
      TYPE(cp_fm_type), DIMENSION(:)                     :: mos
      TYPE(section_vals_type), POINTER                   :: print_key_section
      CHARACTER(len=*), PARAMETER                        :: routineN = "output_mos_rho"
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: j, re, im,&
                                                            rho_unit_re, rho_unit_im
      CHARACTER(len=14), DIMENSION(4)                    :: file_labels
      
      file_labels(1) = "_SPIN_A_RE.dat"
      file_labels(2) = "_SPIN_A_IM.dat"
      file_labels(3) = "_SPIN_B_RE.dat"
      file_labels(4) = "_SPIN_B_IM.dat"
      logger => cp_get_default_logger()
      ! Start by multiplying the current density by MOS
      @:SPIN_DO(j,re, im, tdagw_env%n_spin)
         rho_unit_re = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(re))
         rho_unit_im = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(im))
         ! Transform the density matrix into molecular orbitals basis and print it out
         ! TODO : Number of molecular orbitals and atomic orbitals might not be the same
         ! TODO : Spins in S_fm
         ! Real part
         ! Scale by spin degeneracy
         ! S * rho
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            tdagw_env%spin_degeneracy, tdagw_env%S_fm%matrix, rho(re)%matrix,&
                            0.0_dp,tdagw_env%rho_workspace(1)%matrix)
         ! C^T * S * rho
         CALL parallel_gemm("T", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, mos(j), tdagw_env%rho_workspace(1)%matrix, 0.0_dp, tdagw_env%rho_workspace(2)%matrix)
         ! C^T * S * rho * S
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, tdagw_env%rho_workspace(2)%matrix, tdagw_env%S_fm%matrix,&
                            0.0_dp, tdagw_env%rho_workspace(1)%matrix)
         ! C^T * S * rho * S * C
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, tdagw_env%rho_workspace(1)%matrix, mos(j), 0.0_dp, tdagw_env%rho_workspace(2)%matrix)
         CALL cp_fm_write_formatted(tdagw_env%rho_workspace(2)%matrix, rho_unit_re)
         ! Imag part
         ! Scale by spin degeneracy
         ! S * rho
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            tdagw_env%spin_degeneracy, tdagw_env%S_fm%matrix, rho(im)%matrix,&
                            0.0_dp, tdagw_env%rho_workspace(1)%matrix)
         ! C^T * S * rho
         CALL parallel_gemm("T", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, mos(j), tdagw_env%rho_workspace(1)%matrix, 0.0_dp, tdagw_env%rho_workspace(2)%matrix)
         ! C^T * rho * S
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, tdagw_env%rho_workspace(2)%matrix, tdagw_env%S_fm%matrix,&
                            0.0_dp, tdagw_env%rho_workspace(1)%matrix)
         ! C^T * rho * S * C
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, tdagw_env%rho_workspace(1)%matrix, mos(j), 0.0_dp, tdagw_env%rho_workspace(2)%matrix)
         CALL cp_fm_write_formatted(tdagw_env%rho_workspace(2)%matrix, rho_unit_im)
         CALL cp_print_key_finished_output(rho_unit_re, logger, print_key_section)
         CALL cp_print_key_finished_output(rho_unit_im, logger, print_key_section)
      END DO
   END SUBROUTINE
! **************************************************************************************************
!> \brief Outputs the matrix in MO basis for matrix components corresponding to covariant representation,
!>        i.e. the Hamiltonian matrix
!> \param tdagw_env Entry point - tdagw environment
!> \param cohsex cohsex matrix in AO basis, covariant representation
!> \param mos Spin dependent molecular orbitals
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE output_mos_covariant(tdagw_env, ham, mos, print_key_section)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: ham
      TYPE(cp_fm_type), DIMENSION(:)                     :: mos
      TYPE(section_vals_type), POINTER                   :: print_key_section
      CHARACTER(len=*), PARAMETER                        :: routineN = "output_mos_cohsex"
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: j, re, im,&
                                                            rho_unit_re, rho_unit_im
      CHARACTER(len=21), DIMENSION(4)                    :: file_labels
      
      file_labels(1) = "_SPIN_A_RE.dat"
      file_labels(2) = "_SPIN_A_IM.dat"
      file_labels(3) = "_SPIN_B_RE.dat"
      file_labels(4) = "_SPIN_B_IM.dat"
      logger => cp_get_default_logger()
      ! Start by multiplying the current density by MOS
      ! TODO : Avoid multiplications if not required
      @:SPIN_DO(j,re, im, tdagw_env%n_spin)
         rho_unit_re = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(re))
         rho_unit_im = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(im))
         ! Transform the density matrix into molecular orbitals basis and print it out
         ! TODO : Number of molecular orbitals and atomic orbitals might not be the same
         ! Real part
         ! C^T * cohsex
         CALL parallel_gemm("T", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, mos(j), ham(re)%matrix, 0.0_dp, tdagw_env%rho_workspace(1)%matrix)
         ! C^T * cohsex * C
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, tdagw_env%rho_workspace(1)%matrix, mos(j),&
                            0.0_dp, tdagw_env%rho_workspace(2)%matrix)
         CALL cp_fm_write_formatted(tdagw_env%rho_workspace(2)%matrix, rho_unit_re)
         ! Imag part
         ! C^T * cohsex
         CALL parallel_gemm("T", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, mos(j), ham(im)%matrix, 0.0_dp, tdagw_env%rho_workspace(1)%matrix)
         ! C^T * cohsex * C
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, tdagw_env%rho_workspace(1)%matrix, mos(j), 0.0_dp, tdagw_env%rho_workspace(2)%matrix)
         CALL cp_fm_write_formatted(tdagw_env%rho_workspace(2)%matrix, rho_unit_im)
         CALL cp_print_key_finished_output(rho_unit_re, logger, print_key_section)
         CALL cp_print_key_finished_output(rho_unit_im, logger, print_key_section)
      END DO
   END SUBROUTINE
! **************************************************************************************************
!> \brief Prints the current field components into a file provided by input
!> \param tdagw_env Entry point - tdagw environment
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE output_field(tdagw_env, print_key_section)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(section_vals_type) , POINTER                  :: print_key_section
      CHARACTER(len=*), PARAMETER                        :: routineN="output_field"
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: field_unit
      
      ! Get logger
      logger => cp_get_default_logger()
      ! Get file descriptor
      field_unit = cp_print_key_unit_nr(logger, print_key_section, extension=".dat")
      ! If the file descriptor is non-zero, output field
      IF (field_unit /= -1) THEN
         WRITE(field_unit, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') tdagw_env%sim_time*femtoseconds,&
               tdagw_env%field(1), tdagw_env%field(2), tdagw_env%field(3)
      END IF
      CALL cp_print_key_finished_output(field_unit, logger, print_key_section)
      
   END SUBROUTINE

! **************************************************************************************************
!> \brief Outputs the expectation value of moments from a given density matrix
!> \note  Moments matrix is provided by the tdagw_env, uses rho_workspace(1:3) 
!> \param tdagw_env Entry point - tdagw environment
!> \param rho Density matrix in AO basis
!> \param rtp_section RTP section of the input parameters, where moments destination may be present
! **************************************************************************************************
   SUBROUTINE output_moments(tdagw_env, rho, mos, print_key_section)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: rho
      TYPE(cp_fm_type), DIMENSION(:)                     :: mos
      TYPE(section_vals_type), POINTER                   :: print_key_section
      CHARACTER(len=*), PARAMETER                        :: routineN = "output_moments"
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: i,j,re,im,&
                                                            moments_unit_re,&
                                                            moments_unit_im
      CHARACTER(len=14), DIMENSION(4)                    :: file_labels
      REAL(kind=dp), DIMENSION(3)                        :: moments
      
      ! Start by getting the relevant file unit
      MARK_USED(mos)
      file_labels(1) = "_SPIN_A_RE.dat"
      file_labels(2) = "_SPIN_A_IM.dat"
      file_labels(3) = "_SPIN_B_RE.dat"
      file_labels(4) = "_SPIN_B_IM.dat"
      logger => cp_get_default_logger()
      @:SPIN_DO(j, re, im, tdagw_env%n_spin)
         moments_unit_re = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(re))
         moments_unit_im = cp_print_key_unit_nr(logger, print_key_section, extension=file_labels(im))
         ! If, for any reason, the file unit is not provided, skip to next cycle immediately
         ! TODO : Either centralize or let user define the output format
         ! TODO : Let user define the output units
         ! TODO : Handle case of complex molecular orbital matrix?
         ! TODO : Specify output units in config
         ! Need to transpose due to the definition of trace function
         CALL cp_fm_transpose(rho(re)%matrix, tdagw_env%rho_workspace(1)%matrix)
         DO i=1,3
            CALL cp_fm_trace(tdagw_env%rho_workspace(1)%matrix, tdagw_env%moments(i)%matrix, moments(i))
            ! Scale by spin degeneracy and electron charge
            moments(i) = - moments(i) * tdagw_env%spin_degeneracy
         END DO
         ! Output to the file
         IF (moments_unit_re > 0) WRITE(moments_unit_re, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') tdagw_env%sim_time*femtoseconds,&
                 moments(1),moments(2), moments(3)
         ! Need to transpose due to the definition of trace function
         CALL cp_fm_transpose(rho(im)%matrix, tdagw_env%rho_workspace(1)%matrix)
         DO i=1,3
            CALL cp_fm_trace(tdagw_env%rho_workspace(1)%matrix, tdagw_env%moments(i)%matrix, moments(i))
         END DO
         ! Output to the file
         IF(moments_unit_im > 0) WRITE(moments_unit_im, '(D16.8,D16.8,D16.8,D16.8)') tdagw_env%sim_time*femtoseconds, moments(1),&
             moments(2), moments(3)
         ! Close the files
         CALL cp_print_key_finished_output(moments_unit_re, logger, print_key_section)
         CALL cp_print_key_finished_output(moments_unit_im, logger, print_key_section)
      END DO
   END SUBROUTINE
! **************************************************************************************************
!> \brief Outputs the restart info (last finished iteration step) + restard density matrix
!> \param restart_section Print key section for the restart files
!> \param rho Density matrix in AO basis
!> \param time_index Time index to be written into the info file
! **************************************************************************************************
   SUBROUTINE output_restart(restart_section, rho, time_index)
      TYPE(section_vals_type), POINTER                   :: restart_section
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: rho
      CHARACTER(len=*), PARAMETER                        :: routineN="output_restart"
      CHARACTER(len=14), DIMENSION(4)                    :: file_labels
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: rho_unit_nr, j, re, im, n_spin, time_index

      ! Default labels distinguishing up to two spin species and real/imaginary parts
      file_labels(1) = "_SPIN_A_RE.dat"
      file_labels(2) = "_SPIN_A_IM.dat"
      file_labels(3) = "_SPIN_B_RE.dat"
      file_labels(4) = "_SPIN_B_IM.dat"

      logger => cp_get_default_logger()
      n_spin = SIZE(rho) /2

      @:SPIN_DO(j, re, im, n_spin)
         ! Real part
         rho_unit_nr = cp_print_key_unit_nr(logger, restart_section, extension=file_labels(re),&
                                            file_form="UNFORMATTED", file_position="REWIND")
         CALL cp_fm_write_unformatted(rho(re)%matrix, rho_unit_nr)
         CALL cp_print_key_finished_output(rho_unit_nr, logger, restart_section)
         ! Imag part
         rho_unit_nr = cp_print_key_unit_nr(logger, restart_section, extension=file_labels(im),&
                                            file_form="UNFORMATTED", file_position="REWIND")
         CALL cp_fm_write_unformatted(rho(im)%matrix, rho_unit_nr)
         CALL cp_print_key_finished_output(rho_unit_nr, logger, restart_section)
         ! Info
         rho_unit_nr = cp_print_key_unit_nr(logger, restart_section, extension=".info",&
                                            file_form="UNFORMATTED", file_position="REWIND")
         IF (rho_unit_nr > 0) WRITE(rho_unit_nr) time_index
         CALL cp_print_key_finished_output(rho_unit_nr, logger, restart_section)
      END DO
   END SUBROUTINE output_restart
! **************************************************************************************************
!> \brief Reads the density matrix from restart files and updates the starting time
!> \param restart_section Print key section for the restart files
!> \param rho Density matrix in AO basis
!> \param time_index Time index to be written into the info file
! **************************************************************************************************
   SUBROUTINE read_restart(tdagw_env, restart_section, run_restart)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(section_vals_type), POINTER                   :: restart_section
      LOGICAL                                            :: run_restart
      CHARACTER(len=*), PARAMETER                        :: routineN="read_restart"
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(len=default_path_length)                 :: save_name, save_name_2
      INTEGER                                            :: rho_unit_nr, j, re, im
      CHARACTER(len=14), DIMENSION(4)                    :: file_labels

      run_restart = .FALSE.
      logger => cp_get_default_logger()
      ! Start by probing/loading info file
      save_name = cp_print_key_generate_filename(logger, restart_section, extension=".info", my_local=.FALSE.)
      IF (file_exists(save_name)) THEN
         CALL open_file(save_name, file_status="OLD", file_form="UNFORMATTED", file_action="READ",&
                        unit_number=rho_unit_nr)
         ! TODO : Some error catching and handling
         READ(rho_unit_nr) tdagw_env%sim_start
         CALL close_file(rho_unit_nr)
         run_restart = .TRUE.
         IF (tdagw_env%unit_nr > 0) WRITE(tdagw_env%unit_nr, *) "Starting from timestep ", tdagw_env%sim_start,&
                 ", delta kick NOT applied"
      ELSE
         CPWARN("Restart required but no info file found - starting from sim_step given in input")
      END IF

      ! Default labels distinguishing up to two spin species and real/imaginary parts
      file_labels(1) = "_SPIN_A_RE.dat"
      file_labels(2) = "_SPIN_A_IM.dat"
      file_labels(3) = "_SPIN_B_RE.dat"
      file_labels(4) = "_SPIN_B_IM.dat"
      @:SPIN_DO(j, re, im, tdagw_env%n_spin)
         save_name = cp_print_key_generate_filename(logger, restart_section, extension=file_labels(re), my_local=.FALSE.)
         save_name_2 = cp_print_key_generate_filename(logger, restart_section, extension=file_labels(im), my_local=.FALSE.)
         IF (file_exists(save_name) .AND. file_exists(save_name_2)) THEN
            CALL open_file(save_name, file_status="OLD", file_form="UNFORMATTED", file_action="READ",&
                           unit_number=rho_unit_nr)
            CALL cp_fm_read_unformatted(tdagw_env%rho(re)%matrix, rho_unit_nr)
            CALL close_file(rho_unit_nr)
            CALL open_file(save_name_2, file_status="OLD", file_form="UNFORMATTED", file_action="READ",&
                           unit_number=rho_unit_nr)
            CALL cp_fm_read_unformatted(tdagw_env%rho(im)%matrix, rho_unit_nr)
            CALL close_file(rho_unit_nr)
         ELSE
            ! TODO : Spin dependent warning - for some spins, SCF may be used while for others restart is used
            CPWARN("Restart required but some restart density matrices not found - starting from SCF density.")
         END IF
      END DO
   END SUBROUTINE read_restart
END MODULE rt_tdagw_io

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the propagation via RT-BSE method.
!> \note  The control is handed directly from cp2k_runs
!> \author Maximilian Graml (03.26)
! **************************************************************************************************

MODULE rt_bse_linearized
   USE bibliography, ONLY: Marek2025, &
                           cite_reference
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE force_env_types, ONLY: force_env_type
   USE post_scf_bandstructure_types, ONLY: post_scf_bandstructure_type
   USE cp_fm_types, ONLY: cp_fm_type, &
                          cp_fm_to_fm, &
                          cp_fm_create, &
                          cp_fm_set_all, &
                          cp_fm_release
   USE cp_cfm_types, ONLY: cp_cfm_type, &
                           cp_fm_to_cfm, &
                           cp_cfm_to_cfm, &
                           cp_cfm_to_fm, &
                           cp_cfm_create, &
                           cp_cfm_get_info, &
                           cp_cfm_release
   USE kinds, ONLY: dp
   USE cp_dbcsr_api, ONLY: dbcsr_p_type, &
                           dbcsr_type, &
                           dbcsr_create, &
                           dbcsr_release, &
                           dbcsr_copy, &
                           dbcsr_add, &
                           dbcsr_set, &
                           dbcsr_clear, &
                           dbcsr_iterator_type, &
                           dbcsr_iterator_start, &
                           dbcsr_iterator_stop, &
                           dbcsr_iterator_next_block, &
                           dbcsr_reserve_blocks, &
                           dbcsr_get_num_blocks, &
                           dbcsr_get_block_p, &
                           dbcsr_get_info
   USE dbt_api, ONLY: dbt_clear, &
                      dbt_contract, &
                      dbt_copy_matrix_to_tensor, &
                      dbt_copy_tensor_to_matrix, &
                      dbt_type
   USE libint_2c_3c, ONLY: libint_potential_type
   USE qs_tensors, ONLY: build_2c_integrals, &
                         build_2c_neighbor_lists
   USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type, &
                                     release_neighbor_list_sets
   USE cp_dbcsr_operations, ONLY: copy_dbcsr_to_fm, &
                                  copy_fm_to_dbcsr, &
                                  cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg, ONLY: cp_fm_scale, &
                                 cp_fm_invert, &
                                 cp_fm_transpose, &
                                 cp_fm_column_scale, &
                                 cp_fm_scale_and_add
   USE cp_cfm_basic_linalg, ONLY: cp_cfm_scale_and_add, &
                                  cp_cfm_scale, &
                                  cp_cfm_transpose, &
                                  cp_cfm_norm, &
                                  cp_cfm_trace, &
                                  cp_cfm_column_scale
   USE cp_cfm_diag, ONLY: cp_cfm_geeig
   USE parallel_gemm_api, ONLY: parallel_gemm
   USE qs_moments, ONLY: build_local_moment_matrix
   USE moments_utils, ONLY: get_reference_point
   USE force_env_methods, ONLY: force_env_calc_energy_force
   USE efield_utils, ONLY: make_field
   USE message_passing, ONLY: mp_para_env_type
   USE input_constants, ONLY: rtp_bse_ham_g0w0, &
                              use_mom_ref_zero, &
                              do_bch, &
                              do_exact, &
                              use_rt_restart
   USE rt_bse_types, ONLY: rtbse_env_type, &
                           create_rtbse_env, &
                           release_rtbse_env, &
                           multiply_fm_cfm
   USE rt_bse_io, ONLY: output_moments, &
                        read_field, &
                        output_field, &
                        output_mos_contravariant, &
                        read_restart, &
                        output_restart, &
                        print_timestep_info, &
                        print_etrs_info_header, &
                        print_etrs_info, &
                        print_rtbse_header_info
   USE cp_log_handling, ONLY: cp_logger_type, &
                              cp_get_default_logger
   USE cp_output_handling, ONLY: cp_add_iter_level, &
                                 cp_rm_iter_level, &
                                 cp_iterate
   USE rt_propagation_output, ONLY: print_ft
   USE rt_propagation_utils, ONLY: read_moments

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_bse"

   PUBLIC :: run_propagation_linearized_bse


CONTAINS

! **************************************************************************************************
!> \brief Runs the electron-only real time propagation of the linearized BSE
!> \param force_env Force environment data, entry point of the calculation
! **************************************************************************************************
   SUBROUTINE run_propagation_linearized_bse(force_env)
      TYPE(force_env_type), POINTER                      :: force_env
      CHARACTER(len=*), PARAMETER                        :: routineN = 'run_propagation_linearized_bse'
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      INTEGER                                            :: i, j, k, handle
      LOGICAL                                            :: converged
      REAL(kind=dp)                                      :: metric, enum_re, enum_im, &
                                                            idempotence_dev, a_metric_1, a_metric_2
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      ! To Do: Bibliography information

      logger => cp_get_default_logger()

      ! Run the initial SCF calculation / read SCF restart information
      CALL force_env_calc_energy_force(force_env, calc_force=.FALSE., consistent_energies=.FALSE.)

      ! Allocate all persistant storage and read input that does not need further processing
      CALL create_rtbse_env(rtbse_env, force_env)

      CALL print_linrtbse_header_info(rtbse_env)

      ! Initiate iteration level "MD" in order to copy the structure of other RTP codes
      CALL cp_add_iter_level(logger%iter_info, "MD")
      ! Initialize non-trivial values
      !  - calculates the moment operators
      !  - populates the initial density matrix
      !     - reads the restart density if requested
      !  - reads the field and moment trace from previous runs
      !  - populates overlap and inverse overlap matrices
      !  - calculates/populates the G0W0/KS Hamiltonian, respectively
      !  - calculates the Hartree reference potential
      !  - calculates the COHSEX reference self-energy
      !  - prints some info about loaded files into the output
      !CALL initialize_rtbse_env(rtbse_env)

      ! Setup the time based on the starting step
      ! Assumes identical dt between two runs
      rtbse_env%sim_time = REAL(rtbse_env%sim_start, dp)*rtbse_env%sim_dt
      ! Output 0 time moments and field
      IF (.NOT. rtbse_env%restart_extracted) THEN
         CALL output_field(rtbse_env)
         CALL output_moments(rtbse_env, rtbse_env%rho)
      END IF

      ! Do not apply the delta kick if we are doing a restart calculation
    !   IF (rtbse_env%dft_control%rtp_control%apply_delta_pulse .AND. (.NOT. rtbse_env%restart_extracted)) THEN
    !      CALL apply_delta_pulse(rtbse_env)
    !   END IF

      ! Time loop

      CALL cp_rm_iter_level(logger%iter_info, "MD")

      ! Carry out the FT
      CALL print_ft(rtbse_env%rtp_section, &
                    rtbse_env%moments_trace, &
                    rtbse_env%time_trace, &
                    rtbse_env%field_trace, &
                    rtbse_env%dft_control%rtp_control, &
                    info_opt=rtbse_env%unit_nr)

      ! Deallocate everything
      CALL release_rtbse_env(rtbse_env)

      CALL timestop(handle)
   END SUBROUTINE run_propagation_linearized_bse

   ! **************************************************************************************************
!> \brief Writes the header and basic info to the standard output
!> \param rtbse_env Entry point - rtbse environment
! **************************************************************************************************
   SUBROUTINE print_linrtbse_header_info(rtbse_env)
      TYPE(rtbse_env_type)                               :: rtbse_env
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => cp_get_default_logger()

      IF (rtbse_env%unit_nr > 0) THEN
         WRITE (rtbse_env%unit_nr, *) ''
         WRITE (rtbse_env%unit_nr, '(A)') ' /-----------------------------------------------'// &
            '------------------------------\'
         WRITE (rtbse_env%unit_nr, '(A)') ' |                                               '// &
            '                              |'
         WRITE (rtbse_env%unit_nr, '(A)') ' |              Linearized Real Time Bethe-Salpeter Propagation'// &
            '               |'
         WRITE (rtbse_env%unit_nr, '(A)') ' |                                               '// &
            '                              |'
         WRITE (rtbse_env%unit_nr, '(A)') ' \-----------------------------------------------'// &
            '------------------------------/'
         WRITE (rtbse_env%unit_nr, *) ''


         WRITE (rtbse_env%unit_nr, '(A18,L62)') ' Apply delta pulse', &
            rtbse_env%dft_control%rtp_control%apply_delta_pulse

         WRITE (rtbse_env%unit_nr, '(A)') ''
      END IF

   END SUBROUTINE print_linrtbse_header_info

END MODULE
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
                          cp_fm_release, &
                           cp_fm_get_info
   USE cp_cfm_types, ONLY: cp_cfm_type, &
                           cp_fm_to_cfm, &
                           cp_cfm_to_cfm, &
                           cp_cfm_to_fm, &
                           cp_cfm_create, &
                           cp_cfm_get_info, &
                           cp_cfm_release
   USE cp_cfm_basic_linalg, ONLY: cp_cfm_gemm
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
    USE rt_bse, ONLY:             initialize_rtbse_env, &
            get_idempotence_deviation, &
            antiherm_metric, &
            init_hartree, &
            get_hartree, &
            get_sigma
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
      !  - prints some info about loaded files into the output
      CALL initialize_rtbse_env(rtbse_env)

      !  - populates the initial density matrix
      !     - reads the restart density if requested
      CALL initialize_density_matrix(rtbse_env)
      !  - calculates/populates the G0W0/KS Hamiltonian, respectively
      CALL initialize_singleparticle_hamiltonian(rtbse_env)
      !  - calculates the Hartree reference potential
      CALL initialize_hartree_potential(rtbse_env)
      !  - calculates the COHSEX reference self-energy
      CALL initialize_cohsex_selfenergy(rtbse_env)

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


      ! ********************** Start the time loop **********************
      ! NOTE : Time-loop starts at index sim_start = 0, unless restarted or configured otherwise
      DO i = rtbse_env%sim_start, rtbse_env%sim_nsteps - 1

         ! Update the simulation time
         rtbse_env%sim_time = REAL(i, dp)*rtbse_env%sim_dt
         rtbse_env%sim_step = i

         ! Initial setup - calculate the Hamiltonian
         !CALL update_effective_ham_MO(rtbse_env, rtbse_env%rho)
         
         CALL solve_rk4_timestep(rtbse_env, rtbse_env%rho, rtbse_env%rho_new, converged, k, metric)

         !CALL get_electron_number(rtbse_env, rtbse_env%rho_new, enum_re, enum_im)
         IF (.FALSE.) THEN
            ! Not all of these are used, but they are all good metrics to check the convergence in problematic cases
            ! TODO : Allow for conditional warning
            CALL get_idempotence_deviation(rtbse_env, rtbse_env%rho_new, idempotence_dev)
            DO j = 1, rtbse_env%n_spin
               CALL cp_cfm_to_fm(rtbse_env%sigma_SEX(j), rtbse_env%real_workspace(1), rtbse_env%real_workspace(2))
               CALL antiherm_metric(real_fm=rtbse_env%real_workspace(1), imag_fm=rtbse_env%real_workspace(2), &
                                    workspace=rtbse_env%rho_workspace, metric=a_metric_1)
               CALL antiherm_metric(real_fm=rtbse_env%hartree_curr(j), &
                                    workspace=rtbse_env%rho_workspace, metric=a_metric_2)
            END DO
         END IF
         CALL print_timestep_info(rtbse_env, i, metric, enum_re, k)
         IF (.NOT. converged) CPABORT("ETRS did not converge")
         CALL cp_iterate(logger%iter_info, iter_nr=i, last=(i == rtbse_env%sim_nsteps))
         DO j = 1, rtbse_env%n_spin
            CALL cp_cfm_to_cfm(rtbse_env%rho_new(j), rtbse_env%rho(j))
         END DO
         ! Print the updated field
         CALL output_field(rtbse_env)
         ! If needed, print out the density matrix in MO basis
         CALL output_mos_contravariant(rtbse_env, rtbse_env%rho, rtbse_env%rho_section)
         ! Also handles outputting to memory
         CALL output_moments(rtbse_env, rtbse_env%rho)
         ! Output restart files, so that the restart starts at the following time index
         CALL output_restart(rtbse_env, rtbse_env%rho, i + 1)
      END DO
      ! ********************** End the time loop **********************

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
            '                |'
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

   
! **************************************************************************************************
!> \brief Calculates the initial density matrix, based on the SCF density
!> \param rtbse_env RT-BSE environment
!> \author Maximilian Graml (03/26) 
! **************************************************************************************************
   SUBROUTINE initialize_density_matrix(rtbse_env)
      TYPE(rtbse_env_type), POINTER                :: rtbse_env
      CHARACTER(len=*), PARAMETER                  :: routineN = "initialize_density_matrix"
      INTEGER                                      :: i, handle, ii, jj, nrow_global, ncol_global, &
                                                      nrow_local, ncol_local, &
                                                      i_row_global, j_col_global
      INTEGER, DIMENSION(:), POINTER               :: row_indices, col_indices


      CALL timeset(routineN, handle)

      ! Get distribution of real_workspace 
      CALL cp_fm_get_info(rtbse_env%real_workspace(1), &
                          nrow_global=nrow_global, ncol_global=ncol_global, &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)

      !Ensure that workspace is set to 0
      CALL cp_fm_set_all(rtbse_env%real_workspace(1), 0.0_dp)
      ! Iterate over both spins
      DO i = 1, rtbse_env%n_spin
         DO ii = 1, nrow_local
            i_row_global = row_indices(ii)
            DO jj = 1, ncol_local
               j_col_global = col_indices(jj)
               IF (i_row_global == j_col_global .AND. i_row_global <= rtbse_env%n_occ(i)) THEN
                  rtbse_env%real_workspace(1)%local_data(ii, jj) = 1.0_dp
               END IF
            END DO
         END DO
         ! Sets imaginary part to zero
         CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace(1), mtarget=rtbse_env%rho(i))
         ! Save the reference value for the case of delta kick
         CALL cp_cfm_to_cfm(rtbse_env%rho(i), rtbse_env%rho_orig(i))
      END DO
      
      !MG to do: restart reading of density matrix

      CALL timestop(handle)
   END SUBROUTINE initialize_density_matrix

! **************************************************************************************************
!> \brief Calculates the single particle Hamiltonian in the MO basis
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
!> \author Maximilian Graml (03/26) - refactor in prep. of linearized propagation
! **************************************************************************************************
   SUBROUTINE initialize_singleparticle_hamiltonian(rtbse_env)
      TYPE(rtbse_env_type), POINTER                :: rtbse_env
      CHARACTER(len=*), PARAMETER                  :: routineN = "initialize_singleparticle_hamiltonian"
      TYPE(post_scf_bandstructure_type), POINTER   :: bs_env
      INTEGER                                      :: i, ii,jj, handle, nrow_global, ncol_global, &
                                                      nrow_local, ncol_local, &
                                                      i_row_global, j_col_global
      INTEGER, DIMENSION(:), POINTER               :: row_indices, col_indices

      CALL timeset(routineN, handle)
      ! Get pointers to parameters from qs_env
      CALL get_qs_env(rtbse_env%qs_env, bs_env=bs_env)

      ! Get distribution of real_workspace 
      CALL cp_fm_get_info(rtbse_env%real_workspace(1), &
                          nrow_global=nrow_global, ncol_global=ncol_global, &
                          nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)

      !Ensure that workspace is set to 0
      CALL cp_fm_set_all(rtbse_env%real_workspace(1), 0.0_dp)
      ! ****** START SINGLE PARTICLE HAMILTONIAN CALCULATION
      DO i = 1, rtbse_env%n_spin
         DO ii = 1, nrow_local
            i_row_global = row_indices(ii)
            DO jj = 1, ncol_local
               j_col_global = col_indices(jj)
                  IF (i_row_global == j_col_global) THEN
                     IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0) THEN
                        ! G0W0 Hamiltonian
                        rtbse_env%real_workspace(1)%local_data(ii, jj) = bs_env%eigenval_G0W0(i_row_global, 1, i)
                     ELSE
                        ! KS Hamiltonian
                        rtbse_env%real_workspace(1)%local_data(ii, jj) = bs_env%eigenval_scf_Gamma(i_row_global, i)
                     END IF
                  END IF
            END DO
         END DO
         CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace(i), mtarget=rtbse_env%ham_reference(i))
      END DO
      ! ****** END SINGLE PARTICLE HAMILTONIAN CALCULATION

      CALL timestop(handle)
   END SUBROUTINE initialize_singleparticle_hamiltonian

! **************************************************************************************************
!> \brief Calculates the Hartree potential
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
!> \author Maximilian Graml (03/26) - add transform to MO
! **************************************************************************************************
   SUBROUTINE initialize_hartree_potential(rtbse_env)
      TYPE(rtbse_env_type), POINTER                :: rtbse_env
      CHARACTER(len=*), PARAMETER                  :: routineN = "initialize_hartree_potential"
      TYPE(post_scf_bandstructure_type), POINTER   :: bs_env
      INTEGER                                      :: i, handle

      CALL timeset(routineN, handle)
      ! Get pointers to parameters from qs_env
      CALL get_qs_env(rtbse_env%qs_env, bs_env=bs_env)

      ! ****** START HARTREE POTENTIAL REFERENCE CALCULATION
      ! Calculate Coulomb RI elements, necessary for Hartree calculation
      CALL init_hartree(rtbse_env, rtbse_env%v_dbcsr)
      ! Calculate the original Hartree potential
      ! Uses rho_orig - same as rho for initial run but different for continued run
      DO i = 1, rtbse_env%n_spin
         CALL get_hartree(rtbse_env, rtbse_env%rho_orig(i), rtbse_env%hartree_curr(i))
         ! Scaling by spin degeneracy
         CALL cp_fm_scale(rtbse_env%spin_degeneracy, rtbse_env%hartree_curr(i))
         ! Transform to MO basis
         CALL transform_ao_to_mo_fm(rtbse_env, rtbse_env%hartree_curr(i), rtbse_env%hartree_curr(i))
         ! Subtract the reference from the reference Hamiltonian
         ! following H_eff = ... + V_Hartree(delta_rho) - V_Hartree(rho_0), 
         CALL cp_fm_to_cfm(msourcer=rtbse_env%hartree_curr(i), mtarget=rtbse_env%ham_workspace(1))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                   CMPLX(-1.0, 0.0, kind=dp), rtbse_env%ham_workspace(1))
      END DO
      ! ****** END HARTREE POTENTIAL REFERENCE CALCULATION

      CALL timestop(handle)
   END SUBROUTINE initialize_hartree_potential

! **************************************************************************************************
!> \brief Calculates the single particle Hamiltonian
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
!> \author Maximilian Graml (03/26) - add transform to MO
! **************************************************************************************************
   SUBROUTINE initialize_cohsex_selfenergy(rtbse_env)
      TYPE(rtbse_env_type), POINTER                :: rtbse_env
      CHARACTER(len=*), PARAMETER                  :: routineN = "initialize_cohsex_selfenergy"
      INTEGER                                      :: i, handle

      CALL timeset(routineN, handle)

      ! ****** START COHSEX REFERENCE CALCULATION
      ! Calculate the COHSEX starting energies
      IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0) THEN
         ! Subtract the v_xc from COH part of the self-energy, as V_xc is also not updated during the timestepping
         DO i = 1, rtbse_env%n_spin
            ! TODO : Allow no COH calculation for static screening
            CALL get_sigma(rtbse_env, rtbse_env%sigma_COH(i), -0.5_dp, rtbse_env%S_inv_fm)
            ! Transform to MO basis
            CALL transform_ao_to_mo_fm(rtbse_env, rtbse_env%sigma_COH(i), rtbse_env%sigma_COH(i))
            ! Copy and subtract from the complex reference hamiltonian
            CALL cp_fm_to_cfm(msourcer=rtbse_env%sigma_COH(i), mtarget=rtbse_env%ham_workspace(1))
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                      CMPLX(-1.0, 0.0, kind=dp), rtbse_env%ham_workspace(1))
            ! Calculate exchange part - TODO : should this be applied for different spins? - TEST with O2 HF propagation?
            ! So far only closed shell tested
            ! Uses rho_orig - same as rho for initial run but different for continued run
            CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX(i), -1.0_dp, rtbse_env%rho_orig(i))
            ! Transform to MO basis
            CALL transform_ao_to_mo_cfm(rtbse_env, rtbse_env%sigma_SEX(i), rtbse_env%sigma_SEX(i))
            ! Subtract from the complex reference Hamiltonian
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                      CMPLX(-1.0, 0.0, kind=dp), rtbse_env%sigma_SEX(i))
         END DO
      ELSE
         ! KS Hamiltonian - use time-dependent Fock exchange
         DO i = 1, rtbse_env%n_spin
            CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX(i), -1.0_dp, rtbse_env%rho_orig(i))
            ! Transform to MO basis
            CALL transform_ao_to_mo_cfm(rtbse_env, rtbse_env%sigma_SEX(i), rtbse_env%sigma_SEX(i))
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                      CMPLX(-1.0, 0.0, kind=dp), rtbse_env%sigma_SEX(i))
         END DO
      END IF
      ! ****** END COHSEX REFERENCE CALCULATION

      CALL timestop(handle)
   END SUBROUTINE initialize_cohsex_selfenergy


! **************************************************************************************************
!> \brief Propagates the density to the next timestep using Runge-Kutta 4th order (RK4) method
!> \param rho_start Initial density matrix
!> \param rho_end Final density matrix 
!> \param converged Whether the resulting rho_end is self-consistent
!> \param k How many SC iterations were done
!> \param metric The difference metric from the last self-consistent iteration (for printing/evaluation)
! **************************************************************************************************
   SUBROUTINE solve_rk4_timestep(rtbse_env, rho_start, rho_end, converged, k, metric)
      TYPE(rtbse_env_type), POINTER                     :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER          :: rho_start, &
                                                           rho_mid, &
                                                           rho_end
      LOGICAL                                           :: converged
      INTEGER                                           :: k
      REAL(kind=dp)                                     :: metric
      CHARACTER(len=*), PARAMETER                       :: routineN = "solve_rk4_timestep"
      INTEGER                                           :: j, handle

      CALL timeset(routineN, handle)

            

      CALL timestop(handle)
   END SUBROUTINE solve_rk4_timestep

   ! **************************************************************************************************
!> \brief Updates the effective Hamiltonian in the MO basis, given a density matrix rho
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param qs_env QS env
!> \param rho Real and imaginary parts ( + spin) of the density at current time
! **************************************************************************************************
   SUBROUTINE update_effective_ham_MO(rtbse_env, rho)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho
      CHARACTER(len=*), PARAMETER                        :: routineN = "update_effective_ham_MO"
      INTEGER                                            :: k, j, nspin, handle

      CALL timeset(routineN, handle)
      ! Shorthand
      nspin = rtbse_env%n_spin
      ! Reset the effective Hamiltonian to KS Hamiltonian + G0W0 - reference COHSEX - reference Hartree
      DO j = 1, nspin
         ! Sets the imaginary part to zero
         CALL cp_cfm_to_cfm(rtbse_env%ham_reference(j), rtbse_env%ham_effective(j))
      END DO

      CALL timestop(handle)
   END SUBROUTINE update_effective_ham_MO

! **************************************************************************************************
!> \brief Transforms from AO to MO basis
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param qs_env QS env
!> \param rho Real and imaginary parts ( + spin) of the density at current time
! **************************************************************************************************
   SUBROUTINE transform_ao_to_mo_fm(rtbse_env, fm_ao, fm_mo)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_fm_type)                         :: fm_ao, fm_mo
      CHARACTER(len=*), PARAMETER                        :: routineN = "transform_ao_to_mo_fm"
      INTEGER                                            :: k, j, nspin, handle

      CALL timeset(routineN, handle)
      ! Shorthand
      nspin = rtbse_env%n_spin
      ! Transform from AO to MO basis
      ! M_μν * C_νn
      CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                           1.0_dp, fm_ao, rtbse_env%bs_env%fm_mo_coeff_Gamma(1), &
                           0.0_dp, rtbse_env%real_workspace(1))
      ! (C_μm)^T * M_μν * C_νn
      CALL parallel_gemm("T", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                           1.0_dp, rtbse_env%bs_env%fm_mo_coeff_Gamma(1), rtbse_env%real_workspace(1), &
                           0.0_dp, fm_mo)

      CALL timestop(handle)
   END SUBROUTINE transform_ao_to_mo_fm

! **************************************************************************************************
!> \brief Transforms from AO to MO basis for complex matrices
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param qs_env QS env
!> \param rho Real and imaginary parts ( + spin) of the density at current time
! **************************************************************************************************
   SUBROUTINE transform_ao_to_mo_cfm(rtbse_env, fm_ao, fm_mo)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                         :: fm_ao, fm_mo
      CHARACTER(len=*), PARAMETER                        :: routineN = "transform_ao_to_mo_cfm"
      INTEGER                                            :: k, j, nspin, handle
      TYPE(cp_fm_type)                         :: my_real_workspace

      CALL timeset(routineN, handle)
      ! Shorthand
      nspin = rtbse_env%n_spin
      ! MG To Do: Remove/Clean for final code
      ! Additional workspace needed
      CALL cp_fm_create(my_real_workspace, rtbse_env%bs_env%fm_ks_Gamma(1)%matrix_struct)
      CALL cp_fm_set_all(my_real_workspace, 0.0_dp)
      ! Transform from AO to MO basis
      ! Copy parts to workspace(1/2)
      CALL cp_cfm_to_fm(msource=fm_ao, mtargetr=rtbse_env%real_workspace(1), &
                        mtargeti=rtbse_env%real_workspace(2))
      ! REAL
      ! M_μν * C_νn
      CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                           1.0_dp, rtbse_env%real_workspace(1), rtbse_env%bs_env%fm_mo_coeff_Gamma(1), &
                           0.0_dp, my_real_workspace)
      ! (C_μm)^T * M_μν * C_νn
      CALL parallel_gemm("T", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                           1.0_dp, rtbse_env%bs_env%fm_mo_coeff_Gamma(1), my_real_workspace, &
                           0.0_dp, rtbse_env%real_workspace(1))
      ! IMAGINARY
      ! M_μν * C_νn
      CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                           1.0_dp, rtbse_env%real_workspace(2), rtbse_env%bs_env%fm_mo_coeff_Gamma(1), &
                           0.0_dp, my_real_workspace)
      ! (C_μm)^T * M_μν * C_νn
      CALL parallel_gemm("T", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                           1.0_dp, rtbse_env%bs_env%fm_mo_coeff_Gamma(1), my_real_workspace, &
                           0.0_dp, rtbse_env%real_workspace(2))
      ! Copy back to real part of cfm
      CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace(1), &
                        msourcei=rtbse_env%real_workspace(2), &
                        mtarget=fm_mo)


      CALL timestop(handle)
   END SUBROUTINE transform_ao_to_mo_cfm

END MODULE
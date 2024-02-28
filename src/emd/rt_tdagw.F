!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the real time propagation via TD-aGW method.
!> \note  The control is handed to the routine run_propagation_gw from motion/rt_propagation, where
!>        the REAL_TIME_PROPAGATION section is parsed together with MD section of input
!> \author Stepan Marek (12.23)
! **************************************************************************************************

MODULE rt_tdagw
   USE cp_control_types,                ONLY: dft_control_type,&
                                              rtp_control_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: mo_set_type
   USE rt_propagation_types,            ONLY: get_rtp,&
                                              rt_prop_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE cp_fm_types,                     ONLY: cp_fm_type,&
                                              cp_fm_p_type,&
                                              cp_fm_to_fm,&
                                              cp_fm_create,&
                                              cp_fm_set_all,&
                                              cp_fm_release,&
                                              cp_fm_get_element,&
                                              cp_fm_write_formatted
   USE cp_cfm_types,                    ONLY: cp_cfm_type,&
                                              cp_fm_to_cfm,&
                                              cp_cfm_create,&
                                              cp_cfm_release
   USE kinds,                           ONLY: dp
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_type,&
                                              dbcsr_print,&
                                              dbcsr_has_symmetry,&
                                              dbcsr_desymmetrize,&
                                              dbcsr_create,&
                                              dbcsr_release,&
                                              dbcsr_copy,&
                                              dbcsr_scale,&
                                              dbcsr_set,&
                                              dbcsr_clear,&
                                              dbcsr_setname,&
                                              dbcsr_iterator_type,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop,&
                                              dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_put_block,&
                                              dbcsr_reserve_blocks,&
                                              dbcsr_get_num_blocks,&
                                              dbcsr_get_block_p
   USE OMP_LIB,                         ONLY: omp_get_thread_num,&
                                              omp_get_num_threads,&
                                              omp_set_num_threads,&
                                              omp_get_max_threads
   USE dbt_api,                         ONLY: dbt_create,&
                                              dbt_clear,&
                                              dbt_contract,&
                                              dbt_copy_matrix_to_tensor,&
                                              dbt_copy_tensor_to_matrix,&
                                              dbt_copy,&
                                              dbt_destroy,&
                                              dbt_type,&
                                              dbt_pgrid_type,&
                                              dbt_pgrid_create,&
                                              dbt_pgrid_destroy,&
                                              dbt_mp_environ_pgrid,&
                                              dbt_default_distvec,&
                                              dbt_distribution_type,&
                                              dbt_distribution_new,&
                                              dbt_distribution_destroy,&
                                              dbt_iterator_type,&
                                              dbt_iterator_start,&
                                              dbt_iterator_stop,&
                                              dbt_iterator_blocks_left,&
                                              dbt_iterator_next_block,&
                                              dbt_put_block,&
                                              dbt_get_block,&
                                              dbt_reserve_blocks
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat
   USE qs_tensors,                      ONLY: neighbor_list_3c_destroy,&
                                              build_2c_integrals,&
                                              build_2c_neighbor_lists
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_invert,&
                                              cp_fm_trace,&
                                              cp_fm_transpose,&
                                              cp_fm_norm,&
                                              cp_fm_scale_and_add
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale_and_add,&
                                              cp_cfm_norm
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
   USE efield_utils,                    ONLY: make_field
   USE rt_propagator_init,              ONLY: rt_initialize_rho_from_mos
   USE rt_propagation_methods,          ONLY: s_matrices_create
   USE gw_methods,                      ONLY: compute_3c_integrals
   USE matrix_exp,                      ONLY: taylor_full_complex
   USE cp_log_handling,                 ONLY: cp_logger_type,&
                                              cp_get_default_logger,&
                                              cp_logger_get_unit_nr
   USE cp_output_handling,              ONLY: cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE input_constants,                 ONLY: rtp_tdagw_ham_ks,&
                                              rtp_tdagw_ham_g0w0,&
                                              rtp_tdagw_hartree_dbt,&
                                              rtp_tdagw_hartree_pw
   USE rt_tdagw_types,                  ONLY: tdagw_env_type,&
                                              create_tdagw_env,&
                                              release_tdagw_env
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE qs_integrate_potential_product,  ONLY: integrate_v_rspace
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE physcon,                         ONLY: femtoseconds

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_tdagw"

   #:include "tdagw_macros.fypp"

   PUBLIC :: run_propagation_gw

CONTAINS

! **************************************************************************************************
!> \brief Runs the electron-only real time adiabatic GW (TD Bethe-Salpeter) propagation
!> \param qs_env Quickstep environment data, entry point of the calculation
! **************************************************************************************************
   SUBROUTINE run_propagation_gw(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(len=*), PARAMETER                        :: routineN = 'run_propagation_gw'
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(section_vals_type), POINTER                   :: input, rtp_section
      ! TODO : Maybe complex?
      INTEGER                                            :: i,j,re,im, k
      LOGICAL                                            :: converged
      REAL(kind=dp)                                      :: metric

      CALL get_qs_env(qs_env, bs_env=bs_env, input=input, rtp=rtp)
      ! ********************** ALLOCATIONS **********************
      CALL create_tdagw_env(tdagw_env, qs_env)

      ! Calculate Coulomb RI elements, necessary for Hartree calculation
      CALL init_hartree(qs_env, tdagw_env%hartree_dbt)
      PRINT *, "Hartree initialised"
      ! ********************** EXTRACT PREVIOUS DATA **********************

      ! Subtract the v_xc from COH part of the self-energy, as V_xc is also not updated during the timestepping
      @:SPIN_DO(j, re, im, tdagw_env%n_spin)
         CALL get_sigma(tdagw_env, qs_env, tdagw_env%sigma_COH(re), -0.5_dp, greens_dbcsr_opt=tdagw_env%S_inv)
         ! TODO : Maybe enable printing underneath debug option?
         ! PRINT *, "== Sigma COH =="
         ! CALL cp_fm_write_formatted(sigma_COH, 6)
         CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%sigma_COH(re)%matrix, -1.0_dp, bs_env%fm_V_xc_Gamma(j))
         PRINT *, "metric sigma_COH + Vxc, RE spin", j, cp_fm_norm(tdagw_env%sigma_COH(re)%matrix, "M")
         PRINT *, "metric sigma_COH + Vxc, IM spin", j, cp_fm_norm(tdagw_env%sigma_COH(im)%matrix, "M")
      END DO

      ! Calculate the initial Hartree interaction matrix
      ! TODO : k-points and the possible divergence
      @:SPIN_DO(j, re, im, tdagw_env%n_spin)
         CALL get_hartree(tdagw_env, qs_env, tdagw_env%hartree_dbt, tdagw_env%rho(re)%matrix, tdagw_env%hartree_orig(re)%matrix)
         CALL get_hartree(tdagw_env, qs_env, tdagw_env%hartree_dbt, tdagw_env%rho(im)%matrix, tdagw_env%hartree_orig(im)%matrix)
         PRINT *, "metric hartree_orig, RE spin", j, cp_fm_norm(tdagw_env%hartree_orig(re)%matrix, "M")
         PRINT *, "metric hartree_orig, IM spin", j, cp_fm_norm(tdagw_env%hartree_orig(im)%matrix, "M")
      END DO
      ! TODO : Look at Hartree from qs_ks_env and compare these
      ! Print the RI-based Hartree matrix
      ! PRINT *, "=== RI Hartree matrix"
      ! CALL cp_fm_write_formatted(tdagw_env%hartree_orig(1)%matrix, 6)

      ! Setup the files
      rtp_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION")
      ! ********************** Start the time loop **********************
      DO i = tdagw_env%sim_start, tdagw_env%sim_nsteps

         ! Update the simulation time
         tdagw_env%sim_time = REAL(i, dp)*rtp%dt
         tdagw_env%sim_step = i
         PRINT *, "metric rho, RE spin", j, "time", tdagw_env%sim_time, cp_fm_norm(tdagw_env%rho(re)%matrix, "M")
         PRINT *, "metric rho, IM spin", j, "time", tdagw_env%sim_time, cp_fm_norm(tdagw_env%rho(im)%matrix, "M")
         ! Start the enforced time reversal method
         ! This method determines the density matrix at time (t+dt) by guessing the effective Hamiltonian at (t + dt)
         ! and using the Hamiltonian at time (t), it propagates density from time (t) while ensuring that the density
         ! at (t + dt/2) is the same for both forward and backwards propagation. Then, density at (t + dt) is 
         ! used to calculate the new Hamiltonian at (t+dt), which is then used to get the new propagator, and so on
         ! until the density matrix does not change within certain limit
         ! Pseudocode of the algorithm
         !      rho_M = exp(-i H[rho(t)] S_inv dt/2) rho(t) exp(i S_inv H[rho(t)] dt/2)
         !      rho(t+dt, 0) = rho_M
         !      for j in 0,max_self_iter
         !              rho(t+dt,j+1) = exp(- i H[rho(t+dt,j)] S_inv dt/2) rho_M exp(i S_inv H [rho(t+dt,j)] dt/2)
         !              if ||rho(t+dt,j+1) - rho(t+dt,j)|| < epsilon
         !                      break
         ! *************** Determine rho_M ***************
         ! Update the effective Hamiltonian
         CALL update_effective_ham(tdagw_env, qs_env, tdagw_env%rho)
         PRINT *, "metric ham_effective, RE spin", j, "time", tdagw_env%sim_time,&
                cp_fm_norm(tdagw_env%ham_effective(re)%matrix, "M")
         PRINT *, "metric ham_effective, IM spin", j, "time", tdagw_env%sim_time,&
                cp_fm_norm(tdagw_env%ham_effective(im)%matrix, "M")
         ! Print the updated field TODO : Use printing function that is MPI rank sensitive
         PRINT *, "Field : ", tdagw_env%field(1), tdagw_env%field(2), tdagw_env%field(3)
         ! Convert the effective hamiltonian into the exponential
         CALL ham_to_exp(tdagw_env)
         ! Propagate the density - TODO : So far, only one step - immediately overwrite previous cycle density
         CALL propagate_density(tdagw_env, tdagw_env%ham_workspace, tdagw_env%rho, tdagw_env%rho_M)
         PRINT *, "metric rho_M, RE spin", j, "time", tdagw_env%sim_time, cp_fm_norm(tdagw_env%rho_M(re)%matrix, "M")
         PRINT *, "metric rho_M, IM spin", j, "time", tdagw_env%sim_time, cp_fm_norm(tdagw_env%rho_M(im)%matrix, "M")
         ! In initial iteration, copy rho_M to rho_new_last
         @:SPIN_DO(j,re,im, tdagw_env%n_spin)
            CALL cp_fm_to_fm(tdagw_env%rho_M(re)%matrix, tdagw_env%rho_new_last(re)%matrix)
            CALL cp_fm_to_fm(tdagw_env%rho_M(im)%matrix, tdagw_env%rho_new_last(im)%matrix)
         END DO
         ! *********** Start the self-consistent loop ************************
         tdagw_env%sim_time = REAL(i+1, dp)*rtp%dt
         tdagw_env%sim_step = i+1
         ! TODO : configurable self-iteration limit
         DO k=1,20
            CALL update_effective_ham(tdagw_env, qs_env, tdagw_env%rho_new_last)
            CALL ham_to_exp(tdagw_env)
            CALL propagate_density(tdagw_env, tdagw_env%ham_workspace, tdagw_env%rho_M, tdagw_env%rho_new)
            ! *** Self-consistency check ***
            metric = rho_metric(tdagw_env%rho_new, tdagw_env%rho_new_last, tdagw_env%n_spin)
            PRINT *, "ETRS Iteration : ", k, "Rho convergence metric : ", metric
            PRINT *, "metric rho_new, RE spin", j, "time", tdagw_env%sim_time,&
                "iter", k, cp_fm_norm(tdagw_env%rho_new(re)%matrix, "M")
            PRINT *, "metric rho_new, IM spin", j, "time", tdagw_env%sim_time,&
                "iter", k, cp_fm_norm(tdagw_env%rho_new(im)%matrix, "M")
            IF (metric < tdagw_env%etrs_threshold) THEN
               converged = .TRUE.
               EXIT
            ELSE
               ! Copy rho_new to rho_new_last
               @:SPIN_DO(j, re, im, tdagw_env%n_spin)
                  CALL cp_fm_to_fm(tdagw_env%rho_new(re)%matrix, tdagw_env%rho_new_last(re)%matrix)
                  CALL cp_fm_to_fm(tdagw_env%rho_new(im)%matrix, tdagw_env%rho_new_last(im)%matrix)
               END DO
            END IF
         END DO
         CPASSERT(converged)
         @:SPIN_DO(j, re, im, tdagw_env%n_spin)
            CALL cp_fm_to_fm(tdagw_env%rho_new(re)%matrix, tdagw_env%rho(re)%matrix)
            CALL cp_fm_to_fm(tdagw_env%rho_new(im)%matrix, tdagw_env%rho(im)%matrix)
         END DO
         ! TODO : Implement more relevant outputs
         ! TODO : k-points
         CALL output_mos_rho(tdagw_env, tdagw_env%rho, bs_env%fm_mo_coeff_Gamma, bs_env%fm_s_Gamma, rtp_section)
         CALL output_moments(tdagw_env, tdagw_env%rho, bs_env%fm_mo_coeff_Gamma, rtp_section)
      END DO
      ! ********************** End the time loop **********************

      ! ********************** DEALLOCATIONS **********************
      CALL release_tdagw_env(tdagw_env)
      ! Deallocate the neighbour list that is not deallocated in gw anymore
      IF (ASSOCIATED(bs_env%nl_3c%ij_list)) CALL neighbor_list_3c_destroy(bs_env%nl_3c)
   END SUBROUTINE run_propagation_gw

   FUNCTION rho_metric(rho_new, rho_old, nspin) RESULT(metric)
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: rho_new,&
                                                            rho_old
      INTEGER, INTENT(IN)                                :: nspin
      CHARACTER(len=*), PARAMETER                        :: routineN = "rho_metric"
      TYPE(cp_cfm_type), DIMENSION(2)                    :: workspace
      REAL(kind=dp)                                      :: metric
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: partial_metric
      INTEGER                                            :: j, re, im
      COMPLEX(kind=dp)                                   :: scale_factor
      
      ! TODO : find the max element of rho_new - rho_old
      ! difference is easy, absolute value of complex matrix elements and searching for max less so?
      ! ha, nevermind, cp_fm_norm should be pretty much exactly what is needed
      ! TODO : use cp_fm_to_cfm + cp_cfm_norm to get the norm
      ! Complex matrix for each spin
      ALLOCATE(partial_metric(nspin))
      CALL cp_cfm_create(workspace(1), rho_new(1)%matrix%matrix_struct)
      CALL cp_cfm_create(workspace(2), rho_new(1)%matrix%matrix_struct)
      scale_factor = 1.0
      @:SPIN_DO(j, re, im, nspin)
         ! Create the complex matrix from real matrices (start with new)
         CALL cp_fm_to_cfm(rho_new(re)%matrix, rho_new(im)%matrix, workspace(1))
         ! Copy the other matrix to the workspace
         CALL cp_fm_to_cfm(rho_old(re)%matrix, rho_old(im)%matrix, workspace(2))
         ! Get the difference in the resulting matrix
         CALL cp_cfm_scale_and_add(scale_factor, workspace(1), -scale_factor, workspace(2))
         ! Now, get the relevant number
         partial_metric(j) = cp_cfm_norm(workspace(1), 'M')
      END DO
      metric = 0.0_dp
      ! For more than one spin, do Cartesian sum of the different spin norms
      IF (nspin > 1) THEN
         DO j=1,nspin
            metric = metric + partial_metric(j)*partial_metric(j)
         END DO
         metric = SQRT(metric)
      ELSE
         ! For single spin, just use the norm provided
         metric = partial_metric(1)
      END IF
      ! Deallocate workspace
      CALL cp_cfm_release(workspace(1))
      CALL cp_cfm_release(workspace(2))
      DEALLOCATE(partial_metric)
   END FUNCTION

   SUBROUTINE ham_to_exp(tdagw_env)
      TYPE(tdagw_env_type)                               :: tdagw_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "ham_to_exp"
      INTEGER                                            :: j, re, im
      @:SPIN_DO(j, re, im, tdagw_env%n_spin)

         ! In order to produce correct result, need to remultiply by inverse overlap matrix
         CALL cp_dbcsr_sm_fm_multiply(tdagw_env%S_inv, tdagw_env%ham_effective(re)%matrix,&
                                      tdagw_env%rho_workspace(1)%matrix, tdagw_env%n_ao)
         CALL cp_dbcsr_sm_fm_multiply(tdagw_env%S_inv, tdagw_env%ham_effective(im)%matrix,&
                                      tdagw_env%rho_workspace(2)%matrix, tdagw_env%n_ao)
         CALL cp_fm_to_fm(tdagw_env%rho_workspace(1)%matrix,tdagw_env%ham_effective(re)%matrix)
         CALL cp_fm_to_fm(tdagw_env%rho_workspace(2)%matrix,tdagw_env%ham_effective(im)%matrix)

         ! The evolution of density matrix is in the opposite way to the evolution of operators, i.e. multiply by -dt
         ! TODO : Maybe better to do opposite exponentiation
         CALL cp_fm_scale(tdagw_env%sim_dt/2, tdagw_env%ham_effective(re)%matrix)
         ! Also scale by -1 as we are also multiplying by imaginary unit - this is convention used by exponential routines
         CALL cp_fm_scale(-tdagw_env%sim_dt/2, tdagw_env%ham_effective(im)%matrix)
         ! TODO : Understand which input options set the orders
         CALL taylor_full_complex([tdagw_env%ham_workspace(re)%matrix, tdagw_env%ham_workspace(im)%matrix],&
                                  tdagw_env%ham_effective(im)%matrix, tdagw_env%ham_effective(re)%matrix,&
                                  tdagw_env%orders(1,j), tdagw_env%orders(2,j))
      END DO
   END SUBROUTINE
! **************************************************************************************************
!> \brief Updates the effective Hamiltonian, given a density matrix rho
!> \param tdaggw_env Entry point of the calculation - contains current state of variables
!> \param rho Real and imaginary parts ( + spin) of the density at current time 
! **************************************************************************************************
   SUBROUTINE update_effective_ham(tdagw_env, qs_env, rho)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: rho
      INTEGER                                            :: k, j, re, im, nspin
      ! Shorthand
      nspin = tdagw_env%n_spin
      CALL get_qs_env(qs_env, bs_env=bs_env)
      ! Reset the effective Hamiltonian to KS Hamiltonian + G0W0
      @:SPIN_DO(j,re,im,nspin)
         CALL cp_fm_to_fm(tdagw_env%ham_single_particle(j)%matrix, tdagw_env%ham_effective(re)%matrix)
         ! Imaginary part of KS + G0W0 is zero
         CALL cp_fm_set_all(tdagw_env%ham_effective(im)%matrix, 0.0_dp)
      END DO
      ! Determine the field at current time
      ! TODO : Understand and implement rtp%istep - that should mainly be used for restart calculations?
      CALL make_field(tdagw_env%dft_control, tdagw_env%field, tdagw_env%sim_step, tdagw_env%sim_time)
      IF ( tdagw_env%hartree_method == rtp_tdagw_hartree_pw ) THEN
            ! Use the planewave implementation of Hartree
            ! TODO : Real part, imag part etc.
            CALL get_hartree_pw(qs_env, rho, tdagw_env%hartree_curr(1))
      END IF
      @:SPIN_DO(j, re, im, nspin)
         DO k=1,3
            CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(re)%matrix,&
                                     tdagw_env%field(k), tdagw_env%moments(k)%matrix)
         END DO
         ! Add the COH part
         CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(re)%matrix,&
                                  1.0_dp, tdagw_env%sigma_COH(re)%matrix)
         CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(im)%matrix,&
                                  1.0_dp, tdagw_env%sigma_COH(im)%matrix)
         ! Calculate the SEX part - based on provided rho, not on the 
         ! TODO : Can get away with smaller number of matrices? - depends on memory consumption
         ! iGW = - rho W = - Re(rho) W - i Im(rho) W
         CALL get_sigma(tdagw_env, qs_env, tdagw_env%sigma_SEX(re),&
                        -1.0_dp, greens_fm_opt=rho(re))
         CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(re)%matrix,&
                                  1.0_dp, tdagw_env%sigma_SEX(re)%matrix)
         CALL get_sigma(tdagw_env, qs_env, tdagw_env%sigma_SEX(im), -1.0_dp,&
                        greens_fm_opt=rho(im))
         CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(im)%matrix,&
                                  1.0_dp, tdagw_env%sigma_SEX(im)%matrix)
         ! Calculate update Hartree potential
         ! TODO : Spin (in)dependence check
         IF ( tdagw_env%hartree_method == rtp_tdagw_hartree_dbt ) THEN
            CALL get_hartree(tdagw_env, qs_env, tdagw_env%hartree_dbt, rho(re)%matrix,&
                             tdagw_env%hartree_curr(re)%matrix)
            CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(re)%matrix,&
                                     1.0_dp, tdagw_env%hartree_curr(re)%matrix)
            ! TODO : Check signs
            CALL get_hartree(tdagw_env, qs_env, tdagw_env%hartree_dbt, rho(im)%matrix,&
                             tdagw_env%hartree_curr(im)%matrix)
            CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(im)%matrix,&
                                     1.0_dp, tdagw_env%hartree_curr(im)%matrix)
         END IF
         ! Subtract original Hartree, which is purely real
         CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(re)%matrix,&
                                 -1.0_dp, tdagw_env%hartree_orig(re)%matrix)
         CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(im)%matrix,&
                                 -1.0_dp, tdagw_env%hartree_orig(im)%matrix)
      END DO
   END SUBROUTINE

! **************************************************************************************************
!> \brief Updates the density in tdagw_env, using the provided exponential
!>        The new density is saved to a different matrix, which enables for comparison of matrices
!> \param tdaggw_env Entry point of the calculation - contains current state of variables
!> \param exponential Real and imaginary parts ( + spin) of the exponential propagator 
! **************************************************************************************************
   SUBROUTINE propagate_density(tdagw_env, exponential, rho_old, rho_new)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: exponential,&
                                                            rho_old,&
                                                            rho_new
      INTEGER                                            :: j, re, im
      ! TODO : Consult/think about using complex routines and types
      ! (ham_exp(re)^T - i * ham_exp(im)^T) * (rho_full(re) + i * rho_full(im)) * (ham_exp(re) + i * ham_exp(im)) = 
      ! + ham_exp(re)^T * (rho_full(re) * ham_exp(re) - rho_full(im) * ham_exp(im)) | Term 1
      ! + ham_exp(im)^T * (rho_full(re) * ham_exp(im) + rho_full(im) * ham_exp(re)) | Term 2
      ! - i * ham_exp(im)^T * (rho_full(re) * ham_exp(re) - rho_full(im) * ham_exp(im)) | Term 1
      ! + i * ham_exp(re)^T * (rho_full(re) * ham_exp(im) + rho_full(im) * ham_exp(re)) | Term 2
      ! workspace 2 - real part, workspace 3 - imag part
      ! Real part
      ! TODO : Basis dimensions to tdagw_env
      ! Term 1 - right multiplication first
      @:SPIN_DO(j, re, im, tdagw_env%n_spin)
         CALL parallel_gemm("N","N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, rho_old(re)%matrix, exponential(re)%matrix,&
                            0.0_dp, tdagw_env%rho_workspace(1)%matrix)
         CALL parallel_gemm("N","N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                           -1.0_dp, rho_old(im)%matrix, exponential(im)%matrix,&
                            1.0_dp, tdagw_env%rho_workspace(1)%matrix)
         ! Now term 1 left multiplication real part
         CALL parallel_gemm("T", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, exponential(re)%matrix, tdagw_env%rho_workspace(1)%matrix,&
                            0.0_dp, tdagw_env%rho_workspace(2)%matrix)
         ! Now term 1 left multiplication imag part
         CALL parallel_gemm("T", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                           -1.0_dp, exponential(im)%matrix, tdagw_env%rho_workspace(1)%matrix,&
                            0.0_dp, tdagw_env%rho_workspace(3)%matrix)
         ! Repeat the same procedure for the second term
         ! Term 2 - right multiplication first
         CALL parallel_gemm("N","N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, rho_old(re)%matrix, exponential(im)%matrix,&
                            0.0_dp, tdagw_env%rho_workspace(1)%matrix)
         CALL parallel_gemm("N","N",tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                           -1.0_dp, rho_old(im)%matrix, exponential(re)%matrix,&
                            1.0_dp, tdagw_env%rho_workspace(1)%matrix)
         ! Now term 2 left multiplication real part
         CALL parallel_gemm("T", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, exponential(im)%matrix, tdagw_env%rho_workspace(1)%matrix,&
                            1.0_dp, tdagw_env%rho_workspace(2)%matrix)
         ! Now term 2 left multiplication imag part
         CALL parallel_gemm("T","N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, exponential(re)%matrix, tdagw_env%rho_workspace(1)%matrix,&
                            1.0_dp, tdagw_env%rho_workspace(3)%matrix)
         ! Finally, copy the results to rho - not used anymore TODO : Unless there is spin coupling
         CALL cp_fm_to_fm(tdagw_env%rho_workspace(2)%matrix, rho_new(re)%matrix)
         CALL cp_fm_to_fm(tdagw_env%rho_workspace(3)%matrix, rho_new(im)%matrix)
      END DO
   END SUBROUTINE

! **************************************************************************************************
!> \brief Extracts the relevant information from previous run
!> \param tdagw_env Entry point - tdagw environment
!> \param rho Density matrix in AO basis
!> \param mos Spin dependent molecular orbitals
! **************************************************************************************************
   SUBROUTINE output_mos_rho(tdagw_env, rho, mos, overlap, rtp_section)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: rho
      TYPE(cp_fm_type), DIMENSION(2)                     :: mos
      TYPE(cp_fm_type)                                   :: overlap
      TYPE(section_vals_type), POINTER                   :: rtp_section
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
         rho_unit_re = cp_print_key_unit_nr(logger, rtp_section, "PRINT%DENSITY_MATRIX", extension=file_labels(re))
         rho_unit_im = cp_print_key_unit_nr(logger, rtp_section, "PRINT%DENSITY_MATRIX", extension=file_labels(im))
         ! Transform the density matrix into molecular orbitals basis and print it out
         ! TODO : Number of molecular orbitals and atomic orbitals might not be the same
         !  - Store results in rho_work for now, maybe MO x MO matrix for results is not too expensive?
         ! S * rho
         CALL parallel_gemm("N","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
                         1.0_dp, overlap, rho(re)%matrix, 0.0_dp, tdagw_env%rho_workspace(1)%matrix)
         CALL parallel_gemm("N","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
                         1.0_dp, overlap, rho(im)%matrix, 0.0_dp, tdagw_env%rho_workspace(2)%matrix)
         ! C^T * S * rho
         CALL parallel_gemm("T","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
                         1.0_dp,mos(j),tdagw_env%rho_workspace(1)%matrix,0.0_dp,&
                         tdagw_env%rho_workspace(3)%matrix)
         CALL parallel_gemm("T","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
                         1.0_dp,mos(j),tdagw_env%rho_workspace(2)%matrix,0.0_dp,&
                         tdagw_env%rho_workspace(1)%matrix)
         ! C^T * S * rho * S
         CALL parallel_gemm("N","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
                         1.0_dp, tdagw_env%rho_workspace(3)%matrix,overlap, 0.0_dp,&
                         tdagw_env%rho_workspace(2)%matrix)
         CALL parallel_gemm("N","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
                         1.0_dp, tdagw_env%rho_workspace(1)%matrix,overlap, 0.0_dp,&
                         tdagw_env%rho_workspace(3)%matrix)
         ! C^T * S * rho * S * C
         CALL parallel_gemm("N","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
                         1.0_dp,tdagw_env%rho_workspace(2)%matrix,mos(j),0.0_dp,&
                         tdagw_env%rho_workspace(1)%matrix)
         CALL parallel_gemm("N","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
                         1.0_dp,tdagw_env%rho_workspace(3)%matrix,mos(j),0.0_dp,&
                         tdagw_env%rho_workspace(2)%matrix)
         ! TODO : Save these as separate files
         ! WRITE(UNIT=rho_unit, FMT="(A,I6,A)") "#Iter : ", i, " Real"
         CALL cp_fm_write_formatted(tdagw_env%rho_workspace(1)%matrix, rho_unit_re)
         ! Line below outputs AO basis matrix
         ! CALL cp_fm_write_formatted(rho(re)%matrix, rho_unit_re)
         CALL cp_fm_write_formatted(tdagw_env%rho_workspace(2)%matrix, rho_unit_im)
         ! Line below outputs AO basis matrix
         ! CALL cp_fm_write_formatted(rho(im)%matrix, rho_unit_im)
      END DO
   END SUBROUTINE

! **************************************************************************************************
!> \brief Outputs the expectation value of moments from a given density matrix
!> \note  Moments matrix is provided by the tdagw_env, uses rho_workspace(1:3) 
!> \param tdagw_env Entry point - tdagw environment
!> \param rho Density matrix in AO basis
!> \param rtp_section RTP section of the input parameters, where moments destination may be present
! **************************************************************************************************
   SUBROUTINE output_moments(tdagw_env, rho, mos, rtp_section)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: rho
      TYPE(cp_fm_type), DIMENSION(:)                     :: mos
      TYPE(section_vals_type), POINTER                   :: rtp_section
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
         moments_unit_re = cp_print_key_unit_nr(logger, rtp_section, "PRINT%MOMENTS", extension=file_labels(re))
         moments_unit_im = cp_print_key_unit_nr(logger, rtp_section, "PRINT%MOMENTS", extension=file_labels(im))
         ! If, for any reason, the file unit is not provided, skip to next cycle immediately
         ! TODO : Either centralize or let user define the output format
         ! TODO : Let user define the output units
         ! TODO : Handle case of complex molecular orbital matrix?
         ! TODO : Specify output units in config
         IF (.NOT.(moments_unit_re == -1)) THEN
            DO i=1,3
               ! S_inv * r
               CALL cp_dbcsr_sm_fm_multiply(tdagw_env%S_inv, tdagw_env%moments(i)%matrix,&
                                            tdagw_env%rho_workspace(1)%matrix, tdagw_env%n_ao)
               ! S_inv * rho
               CALL cp_dbcsr_sm_fm_multiply(tdagw_env%S_inv, rho(re)%matrix,&
                                            tdagw_env%rho_workspace(2)%matrix, tdagw_env%n_ao)
               ! Need to transpose due to the definition of trace function
               CALL cp_fm_transpose(tdagw_env%rho_workspace(2)%matrix, tdagw_env%rho_workspace(3)%matrix)
               ! Tr(S_inv * rho * S_inv * r)
               CALL cp_fm_trace(tdagw_env%rho_workspace(3)%matrix, tdagw_env%rho_workspace(1)%matrix, moments(i))
            !   ! r * C
            !   CALL parallel_gemm("N","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
            !                   1.0_dp,tdagw_env%moments(i)%matrix,mos(j),0.0_dp,&
            !                   tdagw_env%rho_workspace(1)%matrix)
            !   ! C^T * r * C
            !   CALL parallel_gemm("T","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
            !                   1.0_dp,mos(j),tdagw_env%rho_workspace(1)%matrix,0.0_dp,&
            !                   tdagw_env%rho_workspace(2)%matrix)
            !   ! C * C^T * r * C
            !   CALL parallel_gemm("N","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
            !                   1.0_dp,mos(j),tdagw_env%rho_workspace(2)%matrix,0.0_dp,&
            !                   tdagw_env%rho_workspace(1)%matrix)
            !   ! rho * C * C^T * r * C
            !   CALL parallel_gemm("N","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
            !                   1.0_dp,rho(re)%matrix,tdagw_env%rho_workspace(1)%matrix,0.0_dp,&
            !                   tdagw_env%rho_workspace(2)%matrix)
            !   ! Tr( C^T * rho * C * C^T * r * C)
            !   CALL cp_fm_trace(mos(j), tdagw_env%rho_workspace(2)%matrix, moments(i))
            !   PRINT *, "Cart index", i, "moment", moments(i)
            END DO
            ! Output to the file
            WRITE(moments_unit_re, '(E16.8,E16.8,E16.8,E16.8)') tdagw_env%sim_time*femtoseconds, moments(1),&
                moments(2), moments(3)
         END IF
         IF (.NOT.(moments_unit_im == -1)) THEN
            DO i=1,3
               ! S_inv * r
               CALL cp_dbcsr_sm_fm_multiply(tdagw_env%S_inv, tdagw_env%moments(i)%matrix,&
                                            tdagw_env%rho_workspace(1)%matrix, tdagw_env%n_ao)
               ! S_inv * rho
               CALL cp_dbcsr_sm_fm_multiply(tdagw_env%S_inv, rho(im)%matrix,&
                                            tdagw_env%rho_workspace(2)%matrix, tdagw_env%n_ao)
               ! Need to transpose due to the definition of trace function
               CALL cp_fm_transpose(tdagw_env%rho_workspace(2)%matrix, tdagw_env%rho_workspace(3)%matrix)
               ! Tr(S_inv * rho * S_inv * r)
               CALL cp_fm_trace(tdagw_env%rho_workspace(3)%matrix, tdagw_env%rho_workspace(1)%matrix, moments(i))
            !   ! r * C
            !   CALL parallel_gemm("N","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
            !                   1.0_dp,tdagw_env%moments(i)%matrix,mos(j),0.0_dp,&
            !                   tdagw_env%rho_workspace(1)%matrix)
            !   ! C^T * r * C
            !   CALL parallel_gemm("T","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
            !                   1.0_dp,mos(j),tdagw_env%rho_workspace(1)%matrix,0.0_dp,&
            !                   tdagw_env%rho_workspace(2)%matrix)
            !   ! C * C^T * r * C
            !   CALL parallel_gemm("N","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
            !                   1.0_dp,mos(j),tdagw_env%rho_workspace(2)%matrix,0.0_dp,&
            !                   tdagw_env%rho_workspace(1)%matrix)
            !   ! rho * C * C^T * r * C
            !   CALL parallel_gemm("N","N",tdagw_env%n_ao,tdagw_env%n_ao,tdagw_env%n_ao,&
            !                   1.0_dp,rho(im)%matrix,tdagw_env%rho_workspace(1)%matrix,0.0_dp,&
            !                   tdagw_env%rho_workspace(2)%matrix)
            !   ! Tr( C^T * rho * C * C^T * r * C)
            !   CALL cp_fm_trace(mos(j), tdagw_env%rho_workspace(2)%matrix, moments(i))
            END DO
            ! Output to the file
            WRITE(moments_unit_im, '(E16.8,E16.8,E16.8,E16.8)') tdagw_env%sim_time*femtoseconds, moments(1),&
                moments(2), moments(3)
         END IF
      END DO
   END SUBROUTINE

! **************************************************************************************************
!> \brief Calculates the self-energy by contraction of screened potential
!> \note Can be used for both the Coulomb hole part and screened exchange part 
!> \param qs_env Quickstep environment data, entry point of the calculation
!> \param greens_fm Pointer to the Green's function matrix, which is used as input data
!> \param sigma_fm Pointer to the self-energy full matrix, which is overwritten by this routine
!> \author Stepan Marek
!> \date 01.2024
! **************************************************************************************************
   SUBROUTINE get_sigma(tdagw_env, qs_env, sigma_fm, prefactor_opt, greens_dbcsr_opt, greens_fm_opt)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type)                                 :: sigma_fm ! resulting self energy
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: prefactor_opt
      TYPE(cp_fm_p_type), INTENT(IN), OPTIONAL           :: greens_fm_opt ! matrix to contract with RI_W
      TYPE(dbcsr_type), POINTER, INTENT(IN), OPTIONAL    :: greens_dbcsr_opt ! matrix, but already in dbcsr form
      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_sigma'
      REAL(kind=dp)                                      :: prefactor
      TYPE(dbcsr_type)                                   :: sigma_dbcsr
      TYPE(dbcsr_type), POINTER                          :: greens_dbcsr
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env,&
                      bs_env=bs_env)

      IF(PRESENT(prefactor_opt)) THEN
         prefactor = prefactor_opt
      ELSE
         prefactor = 1.0_dp
      END IF

      IF (PRESENT(greens_dbcsr_opt)) THEN
         greens_dbcsr => greens_dbcsr_opt
      ELSE IF (PRESENT(greens_fm_opt)) THEN
         ALLOCATE(greens_dbcsr)
         CALL dbcsr_create(greens_dbcsr, name="greens", template=bs_env%mat_ao_ao%matrix)
         CALL copy_fm_to_dbcsr(greens_fm_opt%matrix, greens_dbcsr)
      ELSE
         ! TODO : Warn and exit
         PRINT *, "Either dbcsr or fm greens function has to be provided"
         CPASSERT(.FALSE.)
      END IF

      ! Three-centre integrals are obtained from build_3c_integrals, from qs_tensors
      ! These should use sparcity, while W and Sigma can be full matrices
      ! The summation is carried out by dbt library - dbt_contract in dbt_api
      ! The building of the tensors might be a bit hard, because it requires a lot of parallel information
      ! Probably just use the tensors already present in bs_env? They seem to be mostly work tensors
      ! TODO : Starting with unbounded integrals, then probably need to reduce to a subset
      ! Create by template
      CALL compute_3c_integrals(qs_env,&
                         bs_env,&
                         tdagw_env%t_3c_w)
      CALL dbt_contract(alpha=1.0_dp,&
                   tensor_1=tdagw_env%screened_dbt,&
                   tensor_2=tdagw_env%t_3c_w,&
                   beta=0.0_dp,&
                   tensor_3=tdagw_env%t_3c_work_RI_AO__AO,&
                   contract_1=[2],notcontract_1=[1], map_1=[1],&
                   contract_2=[1],notcontract_2=[2,3], map_2=[2,3],&
                   filter_eps=bs_env%eps_filter)
      ! t_work1 now contains sum _ P (alpha beta | P) W _ (PQ) (iomega = 0)
      ! Next step is to convert the greens full matrix to dbcsr matrix 
      CALL dbt_copy_matrix_to_tensor(greens_dbcsr, tdagw_env%greens_dbt) 
      ! Then contract it
      ! no scaling applied - this has to be applied externally TODO : Is this a good approach?
      CALL dbt_contract(alpha=1.0_dp,&
                   tensor_1=tdagw_env%greens_dbt,&
                   tensor_2=tdagw_env%t_3c_w,&
                   beta=0.0_dp,&
                   tensor_3=tdagw_env%t_3c_work2_RI_AO__AO,&
                   contract_1=[1], notcontract_1=[2], map_1=[2],&
                   contract_2=[3], notcontract_2=[1,2], map_2=[1,3],&
                   filter_eps=bs_env%eps_filter)
      ! Final contraction, which introduces the given prefactor
      CALL dbt_contract(alpha=prefactor,&
                   tensor_1=tdagw_env%t_3c_work_RI_AO__AO,&
                   tensor_2=tdagw_env%t_3c_work2_RI_AO__AO,&
                   beta=0.0_dp,&
                   tensor_3=tdagw_env%sigma_dbt,&
                   contract_1=[1,2], notcontract_1=[3], map_1=[2],&
                   contract_2=[1,2], notcontract_2=[3], map_2=[1],&
                   filter_eps=bs_env%eps_filter)
      ! Finally, convert the COH tensor to matrix and then to fm matrix
      CALL dbcsr_create(sigma_dbcsr, name="sigma", template=bs_env%mat_ao_ao%matrix)
      CALL dbt_copy_tensor_to_matrix(tdagw_env%sigma_dbt, sigma_dbcsr)
      CALL copy_dbcsr_to_fm(sigma_dbcsr, sigma_fm%matrix)
      CALL dbcsr_release(sigma_dbcsr)
      IF (PRESENT(greens_fm_opt)) THEN
         CALL dbcsr_release(greens_dbcsr)
         DEALLOCATE(greens_dbcsr)
      END IF
      ! Clear workspaces - saves memory?
      CALL dbt_clear(tdagw_env%t_3c_w)
      CALL dbt_clear(tdagw_env%t_3c_work_RI_AO__AO)
      CALL dbt_clear(tdagw_env%t_3c_work2_RI_AO__AO)
      CALL dbt_clear(tdagw_env%sigma_dbt)
      CALL dbt_clear(tdagw_env%greens_dbt)
      CALL timestop(handle)

   END SUBROUTINE
! TODO : Maybe combine with get_sigma? - this will be updated at every step and requires rebuilding of 3c?
! **************************************************************************************************
!> \brief Creates the V_dbt tensor and populates it with correct values
!> \note Tensor contains Hartree elements in the auxiliary basis
!> \param qs_env Quickstep environment - entry point of calculation
!> \param V_dbt Tensor representation of the Hartree matrix in the auxiliary basis
!> \author Stepan Marek
!> \date 01.2024
! **************************************************************************************************
   SUBROUTINE init_hartree(qs_env, V_dbt)
      TYPE(qs_environment_type), POINTER, INTENT(IN)     :: qs_env
      TYPE(dbt_type), POINTER, INTENT(OUT)               :: V_dbt
      CHARACTER(len=*), PARAMETER                        :: routineN = "init_hartree"
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(libint_potential_type)                        :: coulomb_op
      TYPE(cp_fm_type)                                   :: V_fm
      TYPE(cp_fm_type)                                   :: metric_fm
      ! TYPE(cp_fm_type), DIMENSION(:,:), ALLOCATABLE      :: metric_fm
      TYPE(cp_fm_type)                                   :: metric_inv_fm,&
                                                            work_fm
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: V_dbcsr,&
                                                            metric_dbcsr
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c

      CALL get_qs_env(qs_env, bs_env=bs_env)

      ! Update KS in order to get the reference element
      ! qs_env%test_hartree_density = .TRUE.
      ! CALL qs_ks_build_kohn_sham_matrix(qs_env, calculate_forces=.FALSE., just_energy=.FALSE.)

      ! Allocate for bare Hartree term
      ! TODO : k-points
      ALLOCATE(V_dbcsr(1))
      ALLOCATE(metric_dbcsr(1))
      CALL dbcsr_create(V_dbcsr(1), name="Hartree_dbcsr", template=bs_env%mat_RI_RI%matrix)
      CALL dbcsr_create(metric_dbcsr(1), name="RI_metric_dbcsr", template=bs_env%mat_RI_RI%matrix)

      ! Calculate full coulomb RI basis elements - V _ (PQ) matrix
      NULLIFY(nl_2c)
      CALL build_2c_neighbor_lists(nl_2c, bs_env%basis_set_RI, bs_env%basis_set_RI,&
                                   coulomb_op, "Coulomb_neighbor_2c_list", qs_env,&
                                   sym_ij=.FALSE., molecular=.TRUE.)
      CALL build_2c_integrals(V_dbcsr, bs_env%eps_filter, qs_env, nl_2c,&
                              bs_env%basis_set_RI, bs_env%basis_set_RI, coulomb_op,&
                              do_kpoints=.FALSE., regularization_RI=bs_env%regularization_RI)
      PRINT *, "Hartree RI dbcsr"
      ! CALL dbcsr_print(V_dbcsr(1))
      ! Calculate the RI metric elements
      ! nl_2c is automatically rewritten (even reallocated) in this routine
      CALL build_2c_neighbor_lists(nl_2c, bs_env%basis_set_RI, bs_env%basis_set_RI,&
                                   bs_env%ri_metric, "Metric_neighbor_2c_list", qs_env,&
                                   sym_ij=.FALSE., molecular=.TRUE.)
      CALL build_2c_integrals(metric_dbcsr, bs_env%eps_filter, qs_env, nl_2c,&
                              bs_env%basis_set_RI, bs_env%basis_set_RI, bs_env%ri_metric,&
                              do_kpoints=.FALSE., regularization_RI=bs_env%regularization_RI)
      ! nl_2c no longer needed
      CALL release_neighbor_list_sets(nl_2c)
      ! CALL RI_2c_integral_mat(qs_env, metric_fm, bs_env%fm_RI_RI, &
      !                         bs_env%n_RI, bs_env%trunc_coulomb, do_kpoints=.FALSE., &
      !                         regularization_RI=bs_env%regularization_RI)
      ! Actual data are in the metric_fm matrix
      ! ! TODO : Spin and k-points
      ! CALL cp_fm_to_fm(metric_fm(1,1), V_fm)
      ! ! Since the routine RI_2c_integral_mat allocates the matrices, they have to be freed again
      ! CALL cp_fm_release(metric_fm(1,1))
      ! DEALLOCATE(metric_fm)
      ! ! Now, calculate the matrices in the standard RI metric
      ! CALL RI_2c_integral_mat(qs_env, metric_fm, bs_env%fm_RI_RI,&
      !                         bs_env%n_RI, bs_env%ri_metric, do_kpoints=.FALSE.,&
      !                         regularization_RI=bs_env%regularization_RI)
      ! ! These matrices have yet to be inverted
      ! CALL cp_fm_create(metric_inv_fm, bs_env%fm_RI_RI%matrix_struct)
      ! CALL cp_fm_create(work_fm, bs_env%fm_RI_RI%matrix_struct)
      ! TODO : Spins and k-points
      ! TODO : Better inverse without conversion?
      CALL cp_fm_create(metric_fm, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(metric_fm, 0.0_dp)
      CALL cp_fm_create(metric_inv_fm, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(metric_inv_fm, 0.0_dp)
      CALL cp_fm_create(work_fm, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(work_fm, 0.0_dp)
      CALL copy_dbcsr_to_fm(metric_dbcsr(1), metric_fm)
      ! CALL cp_fm_invert(metric_fm(1,1), metric_inv_fm)
      CALL cp_fm_invert(metric_fm, metric_inv_fm)
      CALL cp_fm_create(V_fm, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(V_fm, 0.0_dp)
      ! Multiply by the inverse from each site (M^-1 is symmetric)
      CALL cp_dbcsr_sm_fm_multiply(V_dbcsr(1), metric_inv_fm,&
                                   work_fm, bs_env%n_RI)
      CALL parallel_gemm("N", "N", bs_env%n_RI, bs_env%n_RI, bs_env%n_RI,&
                         1.0_dp, metric_inv_fm, work_fm, 0.0_dp, V_fm)
      ! Now, create the tensor from the matrix
      ! First, convert full matrix to dbcsr
      ! TODO : No k-points so far
      CALL dbcsr_clear(V_dbcsr(1))
      CALL copy_fm_to_dbcsr(V_fm, V_dbcsr(1))
      ! Create and copy distinctly, so that unnecessary objects can be destroyed
      CALL dbt_create(bs_env%t_W, V_dbt, name="Hartree_dbt")
      CALL dbt_copy_matrix_to_tensor(V_dbcsr(1), V_dbt)
      ! Destroy all unnecessary matrices
      CALL dbcsr_release(V_dbcsr(1))
      CALL dbcsr_release(metric_dbcsr(1))
      DEALLOCATE(V_dbcsr)
      DEALLOCATE(metric_dbcsr)
      ! TODO : No k-points so far
      CALL cp_fm_release(V_fm)
      ! CALL cp_fm_release(metric_fm(1,1))
      CALL cp_fm_release(metric_fm)
      ! DEALLOCATE(metric_fm)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_release(metric_inv_fm)
   END SUBROUTINE
! **************************************************************************************************
!> \brief Copies the data from a matrix to a 3d tensor, to the first index in the 3rd dimension
!> \param matrix DBCSR Matrix from which the data are copied to the tensor
!> \param tensor Tensor representation of the Hartree matrix in the auxiliary basis
!> \author Stepan Marek
!> \date 01.2024
! **************************************************************************************************
   SUBROUTINE copy_matrix_to_3d_tensor(matrix, tensor)
      TYPE(dbt_type), POINTER, INTENT(OUT)               :: tensor
      TYPE(dbcsr_type), POINTER, INTENT(IN)              :: matrix
      CHARACTER(len=*), PARAMETER                        :: routineN = 'copy_matrix_to_3d_tensor'
      TYPE(dbcsr_type), ALLOCATABLE, TARGET              :: matrix_desym
      TYPE(dbcsr_type), POINTER                          :: matrix_desym_p
      TYPE(dbcsr_iterator_type)                          :: iterator_matrix
      LOGICAL                                            :: data_transpose, found
      INTEGER, DIMENSION(2)                              :: block_shape
      ! Pair of indices for current block
      INTEGER, DIMENSION(3)                              :: index_tensor
      REAL(kind=dp), DIMENSION(:,:), POINTER             :: block_matrix
      REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE       :: block_tensor
      INTEGER                                            :: j,number_of_blocks,backup
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: indices_block_matrix_1,&
                                                            indices_block_matrix_2,&
                                                            indices_block_matrix_custom
      INTEGER                                            :: thread_lower,&
                                                            thread_upper,&
                                                            thread_number,&
                                                            per_thread,&
                                                            per_this_thread,&
                                                            extra_threads,&
                                                            thread_id
      
      ! Desymmetrize the matrix
      IF (dbcsr_has_symmetry(matrix)) THEN
         ALLOCATE(matrix_desym)
         CALL dbcsr_desymmetrize(matrix, matrix_desym)
         matrix_desym_p => matrix_desym
      ELSE
         matrix_desym_p => matrix
      END IF
      ! Get the number of blocks to initiate index vectors
      number_of_blocks = dbcsr_get_num_blocks(matrix_desym_p)
      ALLOCATE(indices_block_matrix_1(number_of_blocks))
      ALLOCATE(indices_block_matrix_2(number_of_blocks))
      ALLOCATE(indices_block_matrix_custom(number_of_blocks), source=1)
      ! TODO : Think about parallelisation
      CALL dbcsr_iterator_start(iterator_matrix, matrix_desym_p)
      DO j=1,number_of_blocks
         CALL dbcsr_iterator_next_block(iterator_matrix, indices_block_matrix_1(j), indices_block_matrix_2(j), block_matrix)
      END DO
      CALL dbcsr_iterator_stop(iterator_matrix)
      ! Construct parallelisation bounds
      ! PRINT *, indices_block_matrix_1
      ! PRINT *, indices_block_matrix_2
      ! PRINT *, indices_block_matrix_custom
      ! Proceed as in dbt - determine distribution of blocks among threads - maybe distvec function? Or integer division
      ! TODO : The reservations are a confusing mess of a spaghetti code - tread lightly
      !        Sure, the statement below does nthreads-times more work than necessary, but it works - any attempts at
      !        parallelisations have failed
      backup = omp_get_max_threads()
      IF (backup > number_of_blocks) THEN
         CALL omp_set_num_threads(number_of_blocks)
         per_thread = 1
         extra_threads = 0
      ELSE
         ! At least one job per thread should be safe ...
         per_thread = number_of_blocks / backup
         extra_threads = MODULO(number_of_blocks, backup)
      END IF
      !$OMP PARALLEL DEFAULT(none) SHARED(tensor, indices_block_matrix_1, indices_block_matrix_2, indices_block_matrix_custom,& 
      !$OMP number_of_blocks, per_thread, extra_threads) &
      !$OMP PRIVATE(thread_lower, thread_upper, per_this_thread, thread_id, thread_number)
         ! Get the omp variables
         thread_number = omp_get_num_threads()
         thread_id = omp_get_thread_num()
         ! Default (higher) bound on reservations to be made by every thread
         per_this_thread = per_thread
         IF (thread_id < extra_threads) THEN
            ! This thread is one of the lucky few who gets to do extra work
            per_this_thread = per_this_thread + 1
         END IF
         ! Get the lower bound for this thread indices, accounting for difference in conventions for omp
         thread_lower = MIN(extra_threads, thread_id) * (per_thread + 1) +&
                        MAX(thread_id - extra_threads, 0) * per_thread + 1
         thread_upper = MIN(thread_lower + per_this_thread - 1, number_of_blocks)
         ! If the lower bound is within the range, some reservations will be made
         IF (thread_lower .LE. number_of_blocks) THEN
            ! Determine the upper bound by minimum
            thread_upper = MIN(thread_lower + per_this_thread - 1, number_of_blocks)
            ! All index fields should be identical 
            CALL dbt_reserve_blocks(tensor,&
                                    indices_block_matrix_1(thread_lower:thread_upper),&
                                    indices_block_matrix_2(thread_lower:thread_upper),&
                                    indices_block_matrix_custom(thread_lower:thread_upper))
         END IF
         ! PRINT *, "thread", thread_id, "done"
      !$OMP END PARALLEL
      ! TODO : Only run when necessary?
      CALL omp_set_num_threads(backup)
      ! PRINT *, "Blocks reserved"
      DO j=1,number_of_blocks
         CALL dbcsr_get_block_p(matrix_desym_p, indices_block_matrix_1(j), indices_block_matrix_2(j),&
                             block_matrix, data_transpose,found, block_shape(1), block_shape(2))
         ALLOCATE(block_tensor(block_shape(1),block_shape(2),1))
         block_tensor(:,:,1) = block_matrix(:,:)
         index_tensor = [indices_block_matrix_1(j), indices_block_matrix_2(j),1]
         CALL dbt_put_block(tensor, index_tensor, SHAPE(block_tensor), block_tensor)
         DEALLOCATE(block_tensor)
      END DO
      ! PRINT *, "Blocks copied"
      DEALLOCATE(indices_block_matrix_1)
      DEALLOCATE(indices_block_matrix_2)
      DEALLOCATE(indices_block_matrix_custom)

      IF (dbcsr_has_symmetry(matrix)) THEN
         CALL dbcsr_release(matrix_desym)
         DEALLOCATE(matrix_desym)
      END IF
   END SUBROUTINE
! **************************************************************************************************
!> \brief Copies the data from a 3d tensor to a matrix, only taking into account fist slice
!>        along the third dimension
!> \param matrix DBCSR Matrix from which the data are copied to the tensor
!> \param tensor Tensor representation of the Hartree matrix in the auxiliary basis
!> \author Stepan Marek
!> \date 01.2024
! **************************************************************************************************
   SUBROUTINE copy_3d_tensor_to_matrix(tensor, matrix)
      TYPE(dbt_type), POINTER, INTENT(IN)                :: tensor
      TYPE(dbcsr_type), POINTER, INTENT(OUT)             :: matrix
      CHARACTER(len=*), PARAMETER                        :: routineN = "copy_3d_tensor_to_matrix"
      TYPE(dbt_iterator_type)                            :: iterator_tensor
      INTEGER, DIMENSION(3)                              :: indices_block_tensor,&
                                                            sizes_block_tensor
      REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE       :: block_tensor
      REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE         :: block_matrix
      LOGICAL                                            :: found
      INTEGER                                            :: backup
      
      ! Start the block iterator
      ! Has to be inside an omp block
      backup = omp_get_max_threads()
      CALL omp_set_num_threads(1)
      !$OMP PARALLEL SHARED(iterator_tensor, tensor, indices_block_tensor, sizes_block_tensor, matrix, &
      !$OMP block_tensor, found, block_matrix)
      CALL dbt_iterator_start(iterator_tensor, tensor)
      DO WHILE (dbt_iterator_blocks_left(iterator_tensor))
         ! Get the tensor indices - last should always be one
         CALL dbt_iterator_next_block(iterator_tensor, indices_block_tensor, sizes_block_tensor)
         CPASSERT(indices_block_tensor(3) == 1)
         ! Reserve the blocks in the matrix
         CALL dbcsr_reserve_blocks(matrix, [indices_block_tensor(1)], [indices_block_tensor(2)]) 
         ! Initiate storage for the block
         ALLOCATE(block_tensor(sizes_block_tensor(1), sizes_block_tensor(2), 1))
         ! Get the block
         CALL dbt_get_block(tensor, indices_block_tensor, sizes_block_tensor, block_tensor, found)
         CPASSERT(found)
         ! Now copy data to matrix
         ALLOCATE(block_matrix(sizes_block_tensor(1), sizes_block_tensor(2)))
         block_matrix(:,:) = block_tensor(:,:,1)
         ! TODO : Handle symmetry and transpose
         CALL dbcsr_put_block(matrix, indices_block_tensor(1), indices_block_tensor(2), block_matrix)
         DEALLOCATE(block_tensor)
         DEALLOCATE(block_matrix)
      END DO
      CALL dbt_iterator_stop(iterator_tensor)
      !$OMP END PARALLEL
      CALL omp_set_num_threads(backup)
   END SUBROUTINE
! **************************************************************************************************
!> \brief Calculates the Hartree matrix in the atomic orbital basis, given a density matrix
!> \note Uses precalculated Hartree matrix kernel in the auxiliary basis. Hartree term is one-electron.
!>       Furthermore uses rho_dbcsr, t_two_plus_one and t_one_plus_one workspaces stored in tdagw_env
!> \param qs_env Quickstep environment - entry point of calculation
!> \param V_dbt Tensor representation of the Hartree matrix in the auxiliary basis
!> \param rho Density matrix in ao basis
!> \param V_ao Overwritten by the Hartree matrix in the atomic orbital basis
!> \author Stepan Marek
!> \date 01.2024
! **************************************************************************************************
   SUBROUTINE get_hartree(tdagw_env, qs_env, V_dbt, rho, V_ao)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbt_type), POINTER, INTENT(INOUT)             :: V_dbt
      TYPE(cp_fm_type), POINTER, INTENT(IN)              :: rho
      TYPE(cp_fm_type), POINTER, INTENT(OUT)             :: V_ao
      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_hartree'
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CALL get_qs_env(qs_env, bs_env=bs_env)

      ! build 3c integrals
      CALL compute_3c_integrals(qs_env,&
                         bs_env,&
                         tdagw_env%t_3c_w)

      ! Convert density to tensor
      CALL copy_fm_to_dbcsr(rho, tdagw_env%rho_dbcsr)
      ! Maybe doable via creation
      CALL copy_matrix_to_3d_tensor(tdagw_env%rho_dbcsr, tdagw_env%t_two_plus_one)
      ! Contract T _ (Q 1) = sum _ Q (Q|mu'nu') D _ (mu' nu' 1)
      ! TODO : Copy filter to tdagw options?
      CALL dbt_contract(alpha=1.0_dp,&
                   tensor_1=tdagw_env%t_3c_w,&
                   tensor_2=tdagw_env%t_two_plus_one,&
                   beta=0.0_dp,&
                   tensor_3=tdagw_env%t_one_plus_one,&
                   contract_1=[2,3],notcontract_1=[1],map_1=[1],&
                   contract_2=[1,2],notcontract_2=[3],map_2=[2],&
                   filter_eps=bs_env%eps_filter)
      ! Contract Q _ (mu nu Q) = sum _ P (mu nu | P) V _ (PQ)
      CALL dbt_contract(alpha=1.0_dp,&
                   tensor_1=tdagw_env%t_3c_w,&
                   tensor_2=V_dbt,&
                   beta=0.0_dp,&
                   tensor_3=tdagw_env%t_3c_work_RI__AO_AO,&
                   contract_1=[1],notcontract_1=[2,3],map_1=[2,3],&
                   contract_2=[1],notcontract_2=[2],map_2=[1],&
                   filter_eps=bs_env%eps_filter)
      ! Final contraction V ^ H _ (mu nu 1) = sum _ (Q) Q _ (mu nu Q) T _ (Q 1)
      CALL dbt_clear(tdagw_env%t_two_plus_one)
      CALL dbt_contract(alpha=1.0_dp,&
                   tensor_1=tdagw_env%t_3c_work_RI__AO_AO,&
                   tensor_2=tdagw_env%t_one_plus_one,&
                   beta=0.0_dp,&
                   tensor_3=tdagw_env%t_two_plus_one,&
                   contract_1=[1], notcontract_1=[2,3], map_1=[1,2],&
                   contract_2=[1], notcontract_2=[2], map_2=[3],&
                   filter_eps=bs_env%eps_filter)
      ! Copy result to output matrix
      CALL dbcsr_clear(tdagw_env%rho_dbcsr)
      CALL copy_3d_tensor_to_matrix(tdagw_env%t_two_plus_one, tdagw_env%rho_dbcsr)
      CALL copy_dbcsr_to_fm(tdagw_env%rho_dbcsr, V_ao)
      ! Free memory by clearing the tensors
      CALL dbt_clear(tdagw_env%t_two_plus_one)
      CALL dbt_clear(tdagw_env%t_one_plus_one)
      CALL dbt_clear(tdagw_env%t_3c_w)
      CALL dbt_clear(tdagw_env%t_3c_work_RI__AO_AO)
      CALL dbcsr_clear(tdagw_env%rho_dbcsr)
   END SUBROUTINE
! **************************************************************************************************
!> \brief Calculates the Hartree matrix in the atomic orbital basis, given a density matrix, using PW approach
!> \param qs_env Entry point
!> \param pw_env Planewave environment for FT calculations
!> \param rho_ao Density matrix in ao basis
!> \param v_ao Overwritten by the Hartree matrix in the atomic orbital basis
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE get_hartree_pw(qs_env, rho_ao, v_ao)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: rho_ao
      TYPE(cp_fm_p_type), INTENT(OUT)                    :: v_ao
      CHARACTER(len=*), PARAMETER                        :: routineN = "get_hartree_pw"
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_ks_env_type), POINTER                      :: qs_ks_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: rho_gspace_tot,&
                                                            v_hartree_gspace,&
                                                            v_hartree_rspace
      TYPE(pw_type), DIMENSION(:), POINTER               :: rho_gspace, rho_rspace
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_dbcsr
      TYPE(dbcsr_p_type)                                 :: v_dbcsr
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      INTEGER                                            :: i, re, im, n_spin
      ! TODO : Imaginary contributions from Hartree? Are these physical?

      ! Get environments and parameters from qs_env
      CALL get_qs_env(qs_env,&
                      dft_control=dft_control,&
                      pw_env=pw_env,&
                      ks_env=qs_ks_env,&
                      matrix_s=matrix_s)
      n_spin = dft_control%nspins
      ! Initialise rho_dbcsr
      ALLOCATE(rho_dbcsr(n_spin))
      @:SPIN_DO(i, re, im, n_spin)
         NULLIFY(rho_dbcsr(re)%matrix)
         ALLOCATE(rho_dbcsr(re)%matrix)
         CALL dbcsr_create(rho_dbcsr(re)%matrix, name="Rho dbcsr", template=matrix_s(i)%matrix)
         CALL copy_fm_to_dbcsr(rho_ao(re)%matrix, rho_dbcsr(re)%matrix)
         NULLIFY(rho_dbcsr(re)%matrix)
         ALLOCATE(rho_dbcsr(re)%matrix)
         CALL dbcsr_create(rho_dbcsr(im)%matrix, name="Rho dbcsr", template=matrix_s(i)%matrix)
         CALL copy_fm_to_dbcsr(rho_ao(im)%matrix, rho_dbcsr(im)%matrix)
      END DO
      ! Initialise v_dbcsr
      ALLOCATE(v_dbcsr%matrix)
      CALL dbcsr_create(v_dbcsr%matrix, template=matrix_s(1)%matrix)
      ! TODO : Redo via block reservations?
      ! TODO : Take into account symmetry present in integrate_v_rspace
      CALL dbcsr_copy(v_dbcsr%matrix, rho_dbcsr(1)%matrix)
      CALL dbcsr_scale(v_dbcsr%matrix, 0.0_dp)
      ! Get the task pool and poisson environment to create pw types
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)
      ! Allocate the densities in g-space and r-space, for each spin
      ALLOCATE(rho_gspace(n_spin))
      ALLOCATE(rho_rspace(n_spin))
      CALL auxbas_pw_pool%create_pw(rho_gspace_tot, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_zero(rho_gspace_tot)
      @:SPIN_DO(i,re,im,n_spin)
         ! Prepare the total density
         CALL auxbas_pw_pool%create_pw(rho_gspace(i), use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         CALL pw_zero(rho_gspace(i))
         CALL auxbas_pw_pool%create_pw(rho_rspace(i), use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_zero(rho_rspace(i))
         ! get the density in rspace and gspace - rho_dbcsr has to be target
         ! CALL calculate_rho_elec(rho_dbcsr(i)%matrix, rho=rho_rspace(i), rho_gspace=rho_gspace(i),&
         !                         ks_env=qs_ks_env, total_rho=total_rho)
         ! PRINT *, "Total rho for spin type ", i, total_rho
         CALL calculate_rho_elec(rho_dbcsr(i)%matrix, rho=rho_rspace(i), rho_gspace=rho_gspace(i),&
                                 ks_env=qs_ks_env)
         ! Do not need rho rspace
         CALL auxbas_pw_pool%give_back_pw(rho_rspace(i))
         CALL pw_axpy(rho_gspace(i), rho_gspace_tot)
         CALL auxbas_pw_pool%give_back_pw(rho_gspace(i))
      END DO
      CALL auxbas_pw_pool%create_pw(v_hartree_gspace, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_zero(v_hartree_gspace)
      ! Use the g-space density to determine the Hartree in g-space
      CALL pw_poisson_solve(poisson_env, rho_gspace_tot, vhartree=v_hartree_gspace)
      CALL auxbas_pw_pool%give_back_pw(rho_gspace_tot)
      ! Transfer to real space
      CALL auxbas_pw_pool%create_pw(v_hartree_rspace, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_zero(v_hartree_rspace)
      CALL pw_transfer(v_hartree_gspace, v_hartree_rspace)
      CALL auxbas_pw_pool%give_back_pw(v_hartree_gspace)
      ! Scale by volume element
      CALL pw_scale(v_hartree_rspace, v_hartree_rspace%pw_grid%dvol)
      ! integrate in rspace
      ! Copy high level in order to create relevant blocks (that is all blocks)
      CALL integrate_v_rspace(v_rspace=v_hartree_rspace, hmat=v_dbcsr, pmat=rho_dbcsr(1), qs_env=qs_env, calculate_forces=.FALSE.)
      CALL auxbas_pw_pool%give_back_pw(v_hartree_rspace)
      DEALLOCATE(rho_gspace)
      DEALLOCATE(rho_rspace)
      CALL copy_dbcsr_to_fm(v_dbcsr%matrix, v_ao%matrix)
      CALL dbcsr_release(v_dbcsr%matrix)
      DEALLOCATE(v_dbcsr%matrix)
      @:SPIN_DO(i, re, im, n_spin)
         CALL dbcsr_release(rho_dbcsr(re)%matrix)
         CALL dbcsr_release(rho_dbcsr(im)%matrix)
      END DO
      !CALL dbcsr_release(rho_dbcsr)
      DEALLOCATE(rho_dbcsr)
   END SUBROUTINE
END MODULE rt_tdagw

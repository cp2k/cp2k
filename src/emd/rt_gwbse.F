!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the real time propagation via GW/BSE method.
!> \note  The control is handed to the routine run_propagation_gw from motion/rt_propagation, where
!>        the REAL_TIME_PROPAGATION section is parsed together with MD section of input
!>        and where the initial guess is calculated/loaded
!> \author Stepan Marek (12.23)
! **************************************************************************************************

MODULE rt_gwbse
   USE cp_control_types,                ONLY: dft_control_type,&
                                              rtp_control_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE force_env_types,                 ONLY: force_env_type
   USE qs_mo_types,                     ONLY: mo_set_type,&
                                              init_mo_set
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE particle_types,                  ONLY: particle_type
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
                                              cp_fm_set_element,&
                                              cp_fm_get_info,&
                                              cp_fm_write_unformatted,&
                                              cp_fm_read_unformatted,&
                                              cp_fm_write_formatted
   USE cp_cfm_types,                    ONLY: cp_cfm_type,&
                                              cp_cfm_p_type,&
                                              cp_fm_to_cfm,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release
   USE kinds,                           ONLY: dp,&
                                              default_path_length
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_type,&
                                              dbcsr_print,&
                                              dbcsr_has_symmetry,&
                                              dbcsr_desymmetrize,&
                                              dbcsr_create,&
                                              dbcsr_release,&
                                              dbcsr_copy,&
                                              dbcsr_scale,&
                                              dbcsr_add,&
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
                                              dbcsr_get_block_p,&
                                              dbcsr_get_info
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
                                              dbt_reserve_blocks,&
                                              dbt_get_num_blocks,&
                                              dbt_get_info
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
                                              cp_dbcsr_sm_fm_multiply,&
                                              copy_cfm_to_dbcsr,&
                                              copy_dbcsr_to_cfm
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_invert,&
                                              cp_fm_trace,&
                                              cp_fm_transpose,&
                                              cp_fm_norm,&
                                              cp_fm_column_scale,&
                                              cp_fm_scale_and_add
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale_and_add,&
                                              cp_cfm_scale,&
                                              cp_cfm_transpose,&
                                              cp_cfm_norm,&
                                              cp_cfm_trace,&
                                              cp_cfm_column_scale
   USE cp_fm_diag,                      ONLY: cp_fm_geeig
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_moments,                      ONLY: build_local_moment_matrix,&
                                              build_berry_moment_matrix
   USE moments_utils,                   ONLY: get_reference_point
   USE qs_mo_io,                        ONLY: read_mo_set_from_restart
   USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
   USE force_env_methods,               ONLY: force_env_calc_energy_force
   USE qs_energy_init,                  ONLY: qs_energies_init
   USE qs_energy_utils,                 ONLY: qs_energies_properties
   USE efield_utils,                    ONLY: make_field
   USE rt_propagator_init,              ONLY: rt_initialize_rho_from_mos
   USE rt_propagation_methods,          ONLY: s_matrices_create
   USE gw_large_cell_Gamma,             ONLY: compute_3c_integrals
   USE gw_integrals,                    ONLY: build_3c_integral_block
   USE matrix_exp,                      ONLY: taylor_full_complex
   USE message_passing,                 ONLY: mp_para_env_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE cell_types,                      ONLY: cell_type
   USE input_constants,                 ONLY: rtp_gwbse_ham_ks,&
                                              rtp_gwbse_ham_g0w0,&
                                              rtp_gwbse_hartree_dbt,&
                                              rtp_gwbse_hartree_pw,&
                                              use_mom_ref_coac,&
                                              use_mom_ref_zero,&
                                              do_taylor,&
                                              do_bch,&
                                              do_exact,&
                                              use_rt_restart
   USE rt_gwbse_types,                  ONLY: gwbse_env_type,&
                                              create_gwbse_env,&
                                              release_gwbse_env,&
                                              multiply_cfm_fm,&
                                              multiply_fm_cfm
   USE rt_gwbse_io,                     ONLY: output_moments,&
                                              output_field,&
                                              output_mos_contravariant,&
                                              read_restart,&
                                              output_restart
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_integrate_potential_product,  ONLY: integrate_v_rspace
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE cp_files,                        ONLY: open_file,&
                                              file_exists,&
                                              close_file
   USE physcon,                         ONLY: femtoseconds
   USE mathconstants,                   ONLY: twopi

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_gwbse"

   #:include "rt_gwbse_macros.fypp"

   PUBLIC :: run_propagation_gw

   INTERFACE get_sigma
      MODULE PROCEDURE get_sigma_dbcsr, get_sigma_real, get_sigma_complex
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Runs the electron-only real time adiabatic GW (TD Bethe-Salpeter) propagation
!> \param qs_env Quickstep environment data, entry point of the calculation
! **************************************************************************************************
   SUBROUTINE run_propagation_gw(qs_env, force_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(force_env_type), POINTER                      :: force_env
      CHARACTER(len=*), PARAMETER                        :: routineN = 'run_propagation_gw'
      TYPE(gwbse_env_type), POINTER                      :: gwbse_env
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cell_type), POINTER                           :: cell
      TYPE(section_vals_type), POINTER                   :: input, rtp_section,&
                                                            output_section
      ! TODO : Maybe complex?
      INTEGER                                            :: i,j,k, handle
      LOGICAL                                            :: converged
      REAL(kind=dp)                                      :: metric, enum_re, enum_im,&
                                                            idempotence_dev, a_metric_1, a_metric_2,&
                                                            intensity
      REAL(kind=dp), DIMENSION(3)                        :: kvec
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(rtp_control_type), POINTER                    :: rtp_control

      CALL timeset(routineN, handle)

      ! Run the initial SCF calculation / read restart information
      CALL force_env_calc_energy_force(force_env, calc_force=.FALSE., &
                                       consistent_energies=.FALSE.)

      ! Allocate all persistant storage
      CALL create_gwbse_env(gwbse_env, qs_env, force_env)

      ! Initialize non-trivial values
      CALL initialize_gwbse_env(gwbse_env, qs_env)

      CALL get_qs_env(qs_env, bs_env=bs_env, rtp=rtp, cell=cell, input=input, dft_control=dft_control)
      ! TODO : Can this be removed?
      rtp_control => dft_control%rtp_control
      ! TODO : Print a nicer header
      IF (gwbse_env%unit_nr > 0) WRITE(gwbse_env%unit_nr, *) "=== Starting time-dependent adiabatic GW real time propagation ==="

      ! Setup the files
      rtp_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION")
      output_section => section_vals_get_subs_vals(rtp_section, "PRINT%FIELD")
      gwbse_env%sim_time = REAL(gwbse_env%sim_start, dp)*gwbse_env%sim_dt
      ! Output 0 time moments
      IF (.NOT. gwbse_env%restart_extracted) CALL output_moments(gwbse_env, gwbse_env%rho, output_section)

      ! Do not apply the delta kick if we are doing a restart calculation
      IF (rtp_control%apply_delta_pulse .AND. (.NOT. gwbse_env%restart_extracted)) THEN
         ! Extra minus for the propagation of density
         intensity = -rtp_control%delta_pulse_scale
         metric = 0.0_dp
         ! TODO : Decide whether the following formulation employed by TDDFT code is not reasonable
         ! TODO : Output intensity multiplied by factors relevant for the output
         ! kvec(:) = cell%h_inv(1, :)*rtp_control%delta_pulse_direction(1) + &
         !           cell%h_inv(2, :)*rtp_control%delta_pulse_direction(2) + &
         !           cell%h_inv(3, :)*rtp_control%delta_pulse_direction(3)
         ! kvec = kvec*twopi*rtp_control%delta_pulse_scale
         kvec(:) = rtp_control%delta_pulse_direction(:)
         DO j=1,gwbse_env%n_spin
            ! So far no spin dependence, but can be added by different structure of delta pulse
            CALL cp_fm_set_all(gwbse_env%real_workspace(1)%matrix, 0.0_dp)
            DO k=1,3
               CALL cp_fm_scale_and_add(1.0_dp, gwbse_env%real_workspace(1)%matrix,&
                       kvec(k), gwbse_env%moments(k)%matrix)
            END DO
            ! enforce hermiticity of the effective Hamiltonian - apply to components instead?
            CALL cp_fm_transpose(gwbse_env%real_workspace(1)%matrix, gwbse_env%real_workspace(2)%matrix)
            CALL cp_fm_scale_and_add(0.5_dp, gwbse_env%real_workspace(1)%matrix,&
                                     0.5_dp, gwbse_env%real_workspace(2)%matrix)
            IF (gwbse_env%mat_exp_method == do_bch) THEN
               ! Multiply by the S_inv matrix - in the classic ordering
               CALL parallel_gemm("N", "N", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao,&
                                  intensity, gwbse_env%S_inv_fm%matrix, gwbse_env%real_workspace(1)%matrix,&
                                  0.0_dp, gwbse_env%real_workspace(2)%matrix)
               ! Sets real part to zero
               CALL cp_fm_to_cfm(msourcei=gwbse_env%real_workspace(2)%matrix, mtarget=gwbse_env%ham_effective(j)%matrix)
            ELSE IF (gwbse_env%mat_exp_method == do_exact) THEN
               CALL cp_fm_to_cfm(msourcer=gwbse_env%real_workspace(1)%matrix, mtarget=gwbse_env%ham_workspace(j)%matrix)
               CALL cp_cfm_gexp(gwbse_env%ham_workspace(j)%matrix, gwbse_env%S_cfm%matrix,gwbse_env%ham_effective(j)%matrix,&
                                CMPLX(0.0, intensity, kind=dp), gwbse_env%rho_workspace)
            END IF
         END DO
         IF (gwbse_env%mat_exp_method == do_bch) THEN
            CALL bch_propagate(gwbse_env%ham_effective, gwbse_env%rho, gwbse_env%rho_new, gwbse_env%rho_workspace,&
                    threshold_opt=gwbse_env%exp_accuracy)
         ELSE IF (gwbse_env%mat_exp_method == do_exact) THEN
            ! Exponential is constructed - apply it to the density matrix
            CALL propagate_density(gwbse_env, gwbse_env%ham_effective, gwbse_env%rho, gwbse_env%rho_new)
         ELSE
            CPABORT("Only exact and BCH exponential schemes adopted for delta kick")
         END IF
         metric = rho_metric(gwbse_env%rho_new, gwbse_env%rho, gwbse_env%n_spin)
         IF (gwbse_env%unit_nr > 0) WRITE(gwbse_env%unit_nr, *) "Metric difference after delta kick", metric
         ! Copy the new density to the old density
         DO j=1,gwbse_env%n_spin
            CALL cp_cfm_to_cfm(gwbse_env%rho_new(j)%matrix, gwbse_env%rho(j)%matrix)
         END DO
      END IF

      ! ********************** Start the time loop **********************
      DO i = gwbse_env%sim_start, gwbse_env%sim_nsteps

         ! Update the simulation time
         gwbse_env%sim_time = REAL(i, dp)*gwbse_env%sim_dt
         gwbse_env%sim_step = i
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
         CALL update_effective_ham(gwbse_env, qs_env, gwbse_env%rho)
         ! Print the updated field
         output_section => section_vals_get_subs_vals(rtp_section, "PRINT%FIELD")
         CALL output_field(gwbse_env, output_section)
         ! Convert the effective hamiltonian into the exponential
         ! TODO : Store the result more explicitly - explicit dependence on input/output matrices
         CALL ham_to_exp(gwbse_env)
         ! Propagate the density to mid-point
         CALL propagate_density(gwbse_env, gwbse_env%ham_workspace, gwbse_env%rho, gwbse_env%rho_M)
         CALL get_electron_number(gwbse_env, gwbse_env%rho_M, enum_re, enum_im)
         IF (gwbse_env%unit_nr > 0) WRITE (gwbse_env%unit_nr, *) "Electron number RHO_M re : ", enum_re,&
                                                                 "Electron number RHO_M im : ", enum_im
         ! In initial iteration, copy rho_M to rho_new_last
         DO j=1,gwbse_env%n_spin
            CALL cp_cfm_to_cfm(gwbse_env%rho_M(j)%matrix, gwbse_env%rho_new_last(j)%matrix)
         END DO
         ! *********** Start the self-consistent loop ************************
         gwbse_env%sim_time = REAL(i+1, dp)*gwbse_env%sim_dt
         gwbse_env%sim_step = i+1
         converged = .FALSE.
         DO k=1,gwbse_env%etrs_max_iter
            CALL update_effective_ham(gwbse_env, qs_env, gwbse_env%rho_new_last)
            CALL ham_to_exp(gwbse_env)
            CALL cp_cfm_to_fm(msource=gwbse_env%ham_workspace(1)%matrix, mtargetr=gwbse_env%real_workspace(1)%matrix,&
                 mtargeti=gwbse_env%real_workspace(2)%matrix)
            CALL propagate_density(gwbse_env, gwbse_env%ham_workspace, gwbse_env%rho_M, gwbse_env%rho_new)
            ! *** Self-consistency check ***
            metric = rho_metric(gwbse_env%rho_new, gwbse_env%rho_new_last, gwbse_env%n_spin)
            IF (gwbse_env%unit_nr > 0) WRITE(gwbse_env%unit_nr, *) "Self-consistent iteration : ", k,&
                                                                   "convergence metric : ", metric
            IF (metric < gwbse_env%etrs_threshold) THEN
               converged = .TRUE.
               EXIT
            ELSE
               ! Copy rho_new to rho_new_last
               DO j=1,gwbse_env%n_spin
                  ! Leaving for free convergence
                  CALL cp_cfm_to_cfm(gwbse_env%rho_new(j)%matrix, gwbse_env%rho_new_last(j)%matrix)
               END DO
            END IF
         END DO
         CALL get_electron_number(gwbse_env, gwbse_env%rho_new, enum_re, enum_im)
         CALL get_idempotence_deviation(gwbse_env, gwbse_env%rho_new, idempotence_dev)
         DO j=1,gwbse_env%n_spin
            CALL cp_cfm_to_fm(gwbse_env%sigma_SEX(j)%matrix, gwbse_env%real_workspace(1)%matrix, gwbse_env%real_workspace(2)%matrix)
            CALL antiherm_metric(real_fm=gwbse_env%real_workspace(1), imag_fm=gwbse_env%real_workspace(2), workspace=gwbse_env%rho_workspace, metric=a_metric_1)
            CALL antiherm_metric(real_fm=gwbse_env%hartree_curr(j), workspace=gwbse_env%rho_workspace, metric=a_metric_2)
         END DO
         ! TODO : Move to io
         IF (gwbse_env%unit_nr > 0) WRITE(gwbse_env%unit_nr, *) "Sim step : ", gwbse_env%sim_step, "Convergence : ", metric,&
                                                                "Threshold : ", gwbse_env%etrs_threshold,&
                                                                "Electron number", enum_re, enum_im, "SC iterations : ", k,&
                                                                "Idempotence deviation : ", idempotence_dev,&
                                                                "Antiherm. metric Sigma : ", a_metric_1,& 
                                                                "Antiherm. metric Hartree : ", a_metric_2
         CPASSERT(converged)
         DO j=1,gwbse_env%n_spin
            CALL cp_cfm_to_cfm(gwbse_env%rho_new(j)%matrix, gwbse_env%rho(j)%matrix)
         END DO
         ! TODO : k-points
         output_section => section_vals_get_subs_vals(rtp_section, "PRINT%DENSITY_MATRIX")
         CALL output_mos_contravariant(gwbse_env, gwbse_env%rho, bs_env%fm_mo_coeff_Gamma, output_section)
         output_section => section_vals_get_subs_vals(rtp_section, "PRINT%MOMENTS")
         CALL output_moments(gwbse_env, gwbse_env%rho, output_section)
         ! Output restart files, so that the restart starts at the following time index
         CALL output_restart(gwbse_env, gwbse_env%rho, i+1)
      END DO
      ! ********************** End the time loop **********************

      ! ********************** DEALLOCATIONS **********************
      CALL release_gwbse_env(gwbse_env)
      ! Deallocate the neighbour list that is not deallocated in gw anymore
      IF (ASSOCIATED(bs_env%nl_3c%ij_list)) CALL neighbor_list_3c_destroy(bs_env%nl_3c)

      CALL timestop(handle)
   END SUBROUTINE run_propagation_gw

! **************************************************************************************************
!> \brief Calculates the initial values, based on restart/scf density, and other non-trivial values
!> \param gwbse_env GW-BSE environment
!> \param qs_env Quickstep environment (needed for reference to previous calculations)
!> \author Stepan Marek (09.24)
! **************************************************************************************************
   SUBROUTINE initialize_gwbse_env(gwbse_env, qs_env)
      TYPE(gwbse_env_type), POINTER                :: gwbse_env
      TYPE(qs_environment_type), POINTER           :: qs_env
      CHARACTER(len=*), PARAMETER                  :: routineN="initialize_gwbse_env"
      TYPE(post_scf_bandstructure_type), POINTER   :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER    :: moments_dbcsr_p
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER    :: matrix_s
      REAL(kind=dp), DIMENSION(:), POINTER         :: custom_ref_point,&
                                                      occupations
      REAL(kind=dp), DIMENSION(3)                  :: rpoint
      INTEGER                                      :: i,k

      ! Get pointers to parameters from qs_env
      CALL get_qs_env(qs_env, bs_env=bs_env, matrix_s=matrix_s)

      ! ****** START MOMENTS OPERATOR CALCULATION
      ! Construct moments from dbcsr
      NULLIFY(moments_dbcsr_p)
      ALLOCATE(moments_dbcsr_p(3))
      DO k=1,3
         ! Make sure the pointer is empty
         NULLIFY(moments_dbcsr_p(k)%matrix)
         ! Allocate a new matrix that the pointer points to
         ALLOCATE(moments_dbcsr_p(k)%matrix)
         ! Create the matrix storage - matrix copies the structure of overlap matrix
         ! TODO : Do a proper matrix create?
         ! CALL dbcsr_create(moments_dbcsr_p(i)%matrix, name="moments", template=matrix_s(1)%matrix)
         CALL dbcsr_copy(moments_dbcsr_p(k)%matrix, matrix_s(1)%matrix)
      END DO
      ! Run the moment calculation
      ! TODO : Read custom reference point from input
      ! TODO : Maybe use the moments section that is present in DFT?
      NULLIFY(custom_ref_point)
      ALLOCATE(custom_ref_point(3), source=0.0_dp)
      rpoint(:) = 0.0_dp
      CALL get_reference_point(rpoint, qs_env=qs_env, reference=use_mom_ref_coac, ref_point=custom_ref_point)
      CALL build_local_moment_matrix(qs_env, moments_dbcsr_p, 1, rpoint)
      ! Copy to full matrix
      DO k=1,3
         ! Again, matrices are created from overlap template
         CALL copy_dbcsr_to_fm(moments_dbcsr_p(k)%matrix, gwbse_env%moments(k)%matrix)
      END DO
      ! Now, repeat without reference point to get the moments for field
      CALL get_reference_point(rpoint, qs_env=qs_env, reference=use_mom_ref_zero, ref_point=custom_ref_point)
      DEALLOCATE(custom_ref_point)
      CALL build_local_moment_matrix(qs_env, moments_dbcsr_p, 1, rpoint)
      DO k=1,3
         CALL copy_dbcsr_to_fm(moments_dbcsr_p(k)%matrix, gwbse_env%moments_field(k)%matrix)
      END DO

      ! Now can deallocate dbcsr matrices
      DO k=1,3
         CALL dbcsr_release(moments_dbcsr_p(k)%matrix)
         DEALLOCATE(moments_dbcsr_p(k)%matrix)
      END DO
      DEALLOCATE(moments_dbcsr_p)
      ! ****** END MOMENTS OPERATOR CALCULATION

      ! ****** START INITIAL DENSITY MATRIX CALCULATION
      ! TODO : Handle cases, when initial method uses density only? Do such cases exist?
      ! Get the rho from fm_MOS
      ! Uses real orbitals only - no kpoints
      ALLOCATE(occupations(gwbse_env%n_ao))
      ! Iterate over both spins
      DO i=1,gwbse_env%n_spin
         occupations(:) = 0.0_dp
         occupations(1:gwbse_env%n_occ(i)) = 1.0_dp
         ! Create real part
         CALL cp_fm_to_fm(bs_env%fm_mo_coeff_Gamma(i), gwbse_env%real_workspace(1)%matrix)
         CALL cp_fm_column_scale(gwbse_env%real_workspace(1)%matrix, occupations)
         CALL parallel_gemm("N", "T", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao,&
                            1.0_dp, gwbse_env%real_workspace(1)%matrix, bs_env%fm_mo_coeff_Gamma(i),&
                            0.0_dp, gwbse_env%real_workspace(2)%matrix)
         ! Sets imaginary part to zero
         CALL cp_fm_to_cfm(msourcer=gwbse_env%real_workspace(2)%matrix, mtarget=gwbse_env%rho(i)%matrix)
         ! Save the reference value for the case of delta kick
         CALL cp_cfm_to_cfm(gwbse_env%rho(i)%matrix, gwbse_env%rho_orig(i)%matrix)
      END DO
      DEALLOCATE(occupations)
      ! If the restart field is provided, overwrite rho from restart
      IF (gwbse_env%dft_control%rtp_control%initial_wfn == use_rt_restart ) THEN
         CALL read_restart(gwbse_env)
      END IF
      ! ****** END INITIAL DENSITY MATRIX CALCULATION

      ! ****** START OVERLAP + INVERSE OVERLAP CALCULATION
      ! TODO : Are different overlaps for different spins reasonable for us?
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, gwbse_env%S_fm%matrix)
      CALL cp_fm_to_cfm(msourcer=gwbse_env%S_fm%matrix, mtarget=gwbse_env%S_cfm%matrix)
      CALL cp_fm_invert(gwbse_env%S_fm%matrix, gwbse_env%S_inv_fm%matrix)
      ! ****** END OVERLAP + INVERSE OVERLAP CALCULATION

      ! ****** START SINGLE PARTICLE HAMILTONIAN CALCULATION
      DO i=1,gwbse_env%n_spin
         IF ( gwbse_env%ham_reference_type == rtp_gwbse_ham_g0w0 ) THEN
            ! G0W0 Hamiltonian
            CALL cp_fm_to_fm(bs_env%fm_mo_coeff_Gamma(i), gwbse_env%real_workspace(1)%matrix)
            ! TODO : Check that Gamma point is the zero k-point
            ! C * Lambda
            CALL cp_fm_column_scale(gwbse_env%real_workspace(1)%matrix, bs_env%eigenval_G0W0(:,1,i))
            ! C * Lambda * C^T
            CALL parallel_gemm("N", "T", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao,&
                               1.0_dp, gwbse_env%real_workspace(1)%matrix, bs_env%fm_mo_coeff_Gamma(i),&
                               0.0_dp, gwbse_env%real_workspace(2)%matrix)
            ! S * C * Lambda * C^T
            CALL parallel_gemm("N", "N", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao,&
                               1.0_dp, gwbse_env%S_fm%matrix, gwbse_env%real_workspace(2)%matrix,&
                               0.0_dp, gwbse_env%real_workspace(1)%matrix)
            ! S * C * Lambda * C^T * S = H
            CALL parallel_gemm("N", "N", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao,&
                               1.0_dp, gwbse_env%real_workspace(1)%matrix, gwbse_env%S_fm%matrix,&
                               0.0_dp, gwbse_env%real_workspace(2)%matrix)
            CALL cp_fm_to_cfm(msourcer=gwbse_env%real_workspace(2)%matrix, mtarget=gwbse_env%ham_reference(i)%matrix)
         ELSE
            ! KS Hamiltonian
            CALL cp_fm_to_cfm(msourcer=bs_env%fm_ks_Gamma(i), mtarget=gwbse_env%ham_reference(i)%matrix)
         END IF
      END DO
      ! ****** END SINGLE PARTICLE HAMILTONIAN CALCULATION

      ! ****** START HARTREE POTENTIAL REFERENCE CALCULATION
      ! Calculate Coulomb RI elements, necessary for Hartree calculation
      CALL init_hartree(qs_env, gwbse_env%v_dbcsr%matrix)
      IF (gwbse_env%ham_reference_type == rtp_gwbse_ham_g0w0) THEN
         ! In a non-HF calculation, copy the actual correlation part of the interaction
         CALL copy_fm_to_dbcsr(bs_env%fm_W_MIC_freq_zero, gwbse_env%w_dbcsr%matrix)
      ELSE
         ! In HF, correlation is set to zero
         CALL dbcsr_set(gwbse_env%w_dbcsr%matrix, 0.0_dp)
      END IF
      ! Add the Hartree to the screened_dbt tensor - now W = V + W^c
      CALL dbcsr_add(gwbse_env%w_dbcsr%matrix, gwbse_env%v_dbcsr%matrix, 1.0_dp, 1.0_dp)
      CALL dbt_copy_matrix_to_tensor(gwbse_env%w_dbcsr%matrix, gwbse_env%screened_dbt)
      ! Calculate the original Hartree potential
      ! Uses rho_orig - same as rho for initial run but different for continued run
      DO i=1,gwbse_env%n_spin
         CALL get_hartree_local(qs_env, gwbse_env, gwbse_env%rho_orig(i)%matrix, gwbse_env%hartree_curr(i)%matrix)
         ! Scaling by spin degeneracy
         CALL cp_fm_scale(gwbse_env%spin_degeneracy, gwbse_env%hartree_curr(i)%matrix)
         ! Subtract the reference from the reference Hamiltonian
         CALL cp_fm_to_cfm(msourcer=gwbse_env%hartree_curr(i)%matrix, mtarget=gwbse_env%ham_workspace(1)%matrix)
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), gwbse_env%ham_reference(i)%matrix,&
                                   CMPLX(-1.0, 0.0, kind=dp), gwbse_env%ham_workspace(1)%matrix)
      END DO
      ! ****** END HARTREE POTENTIAL REFERENCE CALCULATION

      ! ****** START COHSEX REFERENCE CALCULATION
      ! Calculate the COHSEX starting energies
      IF (gwbse_env%ham_reference_type == rtp_gwbse_ham_g0w0) THEN
         ! Subtract the v_xc from COH part of the self-energy, as V_xc is also not updated during the timestepping
         DO i=1,gwbse_env%n_spin
            ! TODO : Allow no COH calculation for static screening
            CALL get_sigma(gwbse_env, qs_env, gwbse_env%sigma_COH(i), -0.5_dp, gwbse_env%S_inv_fm)
            ! Copy and subtract from the complex reference hamiltonian
            CALL cp_fm_to_cfm(msourcer=gwbse_env%sigma_COH(i)%matrix, mtarget=gwbse_env%ham_workspace(1)%matrix)
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), gwbse_env%ham_reference(i)%matrix,&
                                      CMPLX(-1.0, 0.0, kind=dp), gwbse_env%ham_workspace(1)%matrix)
            ! Uses rho_orig - same as rho for initial run but different for continued run
            CALL get_sigma(gwbse_env, qs_env, gwbse_env%sigma_SEX(i), -1.0_dp, gwbse_env%rho_orig(i))
            ! Subtract from the complex reference Hamiltonian
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), gwbse_env%ham_reference(i)%matrix,&
                                      CMPLX(-1.0, 0.0, kind=dp), gwbse_env%sigma_SEX(i)%matrix)
         END DO
      ELSE
         ! KS Hamiltonian - use time-dependent Fock exchange
         DO i=1,gwbse_env%n_spin
            CALL get_sigma(gwbse_env, qs_env, gwbse_env%sigma_SEX(i), -1.0_dp, gwbse_env%rho_orig(i))
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), gwbse_env%ham_reference(i)%matrix,&
                                      CMPLX(-1.0, 0.0, kind=dp), gwbse_env%sigma_SEX(i)%matrix)
         END DO
      END IF
      ! ****** END COHSEX REFERENCE CALCULATION
   END SUBROUTINE initialize_gwbse_env

! **************************************************************************************************
!> \brief Determines the metric for the density matrix, used for convergence criterion
!> \param rho_new Array of new density matrices (one for each spin index)
!> \param rho_old Array of old density matrices (one for each spin index)
!> \param nspin Number of spin indices
! **************************************************************************************************
! TODO : Provide workspace by outside body
   FUNCTION rho_metric(rho_new, rho_old, nspin, workspace_opt) RESULT(metric)
      TYPE(cp_cfm_p_type), DIMENSION(:), INTENT(IN)      :: rho_new,&
                                                            rho_old
      INTEGER, INTENT(IN)                                :: nspin
      TYPE(cp_cfm_type), POINTER, OPTIONAL               :: workspace_opt
      CHARACTER(len=*), PARAMETER                        :: routineN = "rho_metric"
      TYPE(cp_cfm_type)                                  :: workspace
      REAL(kind=dp)                                      :: metric
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: partial_metric
      INTEGER                                            :: j
      COMPLEX(kind=dp)                                   :: scale_factor
      
      ALLOCATE(partial_metric(nspin))

      ! Only allocate/deallocate storage if required
      IF (PRESENT(workspace_opt)) THEN
         workspace = workspace_opt
      ELSE
         CALL cp_cfm_create(workspace, rho_new(1)%matrix%matrix_struct)
      END IF
      scale_factor = 1.0
      DO j=1,nspin
         CALL cp_cfm_to_cfm(rho_new(j)%matrix, workspace)
         ! Get the difference in the resulting matrix
         CALL cp_cfm_scale_and_add(scale_factor, workspace, -scale_factor, rho_old(j)%matrix)
         ! Now, get the relevant number
         partial_metric(j) = cp_cfm_norm(workspace, 'M')
      END DO
      metric = 0.0_dp
      ! For more than one spin, do Cartesian sum of the different spin norms
      DO j=1,nspin
         metric = metric + partial_metric(j)*partial_metric(j)
      END DO
      metric = SQRT(metric)
      ! Deallocate workspace
      IF (.NOT. PRESENT(workspace_opt)) CALL cp_cfm_release(workspace)
      DEALLOCATE(partial_metric)
   END FUNCTION

! **************************************************************************************************
!> \brief Determines the metric of the antihermitian part of the matrix
!> \param real_fm Real part of the full matrix
!> \param imag_fm Imaginary part of the full matrix
! **************************************************************************************************
   SUBROUTINE antiherm_metric(real_fm, imag_fm, workspace, metric)
      TYPE(cp_fm_p_type), INTENT(IN)                    :: real_fm
      TYPE(cp_fm_p_type), INTENT(IN), OPTIONAL          :: imag_fm
      REAL(kind=dp), INTENT(OUT)                        :: metric
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER        :: workspace
      COMPLEX(kind=dp)                                  :: complex_one

      ! Get the complex and complex conjugate matrix
      IF (PRESENT(imag_fm)) THEN
         CALL cp_fm_to_cfm(real_fm%matrix, imag_fm%matrix, workspace(1)%matrix)
      ELSE
         CALL cp_fm_to_cfm(msourcer=real_fm%matrix, mtarget=workspace(1)%matrix)
      END IF
      CALL cp_cfm_transpose(workspace(1)%matrix, "C", workspace(2)%matrix)
      ! Subtract these, and get the metric
      complex_one = CMPLX(1.0, 0.0, kind=dp)
      CALL cp_cfm_scale_and_add(complex_one, workspace(1)%matrix, -complex_one, workspace(2)%matrix)
      metric = cp_cfm_norm(workspace(1)%matrix, "M")
   END SUBROUTINE

! **************************************************************************************************
!> \brief For Taylor and Exact exp_method, calculates the matrix exponential of the
!>        effective Hamiltonian. For BCH, calculates just the effective Hamiltonian. For other methods,
!>        aborts the execution, as they are not implemented yet.
!> \param gwbse_env Entry point of the calculation. Uses rho_workspace for Taylor and BCH. For exact,
!>                  uses complex_workspace, complex_ham, complex_s, real_eigvals and exp_eigvals.
!>                  Results are stored in ham_workspace.
! **************************************************************************************************
   SUBROUTINE ham_to_exp(gwbse_env)
      TYPE(gwbse_env_type)                               :: gwbse_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "ham_to_exp"
      INTEGER                                            :: j, handle
      CALL timeset(routineN, handle)
      DO, j=1,gwbse_env%n_spin
         IF (gwbse_env%mat_exp_method == do_taylor .OR. gwbse_env%mat_exp_method == do_bch) THEN
            ! In Taylor and BCH, we first evaluate the entire exponent and then evaluate exponential in series
            ! In order to produce correct result, need to remultiply by inverse overlap matrix
            CALL multiply_fm_cfm("N", "N", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao,&
                                 1.0_dp, gwbse_env%S_inv_fm%matrix, gwbse_env%ham_effective(j)%matrix,&
                                 0.0_dp, gwbse_env%rho_workspace(1)%matrix)

            ! The evolution of density matrix is derived from the right multiplying term
            ! Imaginary part of the exponent = -real part of the matrix
            CALL cp_cfm_scale(CMPLX(0.0, -gwbse_env%sim_dt/2, kind=dp), gwbse_env%rho_workspace(1)%matrix)
            ! TODO : Understand which input options set the orders
            IF (gwbse_env%mat_exp_method == do_taylor) THEN
               CPABORT("Taylor exponentiaton no longer supported")
            ELSE
               ! In BCH, exponential is not calculated explicitly, but the propagation is solved in series
               CALL cp_cfm_to_cfm(gwbse_env%rho_workspace(1)%matrix, gwbse_env%ham_workspace(j)%matrix)
            END IF
         ELSE IF (gwbse_env%mat_exp_method == do_exact) THEN
            CALL cp_cfm_gexp(gwbse_env%ham_effective(j)%matrix, gwbse_env%S_cfm%matrix, gwbse_env%ham_workspace(j)%matrix,&
                             CMPLX(0.0, -gwbse_env%sim_dt/2, kind=dp), gwbse_env%rho_workspace)
         ELSE
            CPABORT("Only BCH and Taylor matrix exponentiation implemented")
         END IF
      END DO

      CALL timestop(handle)
   END SUBROUTINE
! **************************************************************************************************
!> \brief Updates the effective Hamiltonian, given a density matrix rho
!> \param tdaggw_env Entry point of the calculation - contains current state of variables
!> \param rho Real and imaginary parts ( + spin) of the density at current time 
! **************************************************************************************************
   SUBROUTINE update_effective_ham(gwbse_env, qs_env, rho)
      TYPE(gwbse_env_type)                               :: gwbse_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_cfm_p_type), DIMENSION(:)                  :: rho
      CHARACTER(len=*), PARAMETER                        :: routineN = "update_effective_ham"
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER                                            :: k, j, nspin, handle

      CALL timeset(routineN, handle)
      ! Shorthand
      nspin = gwbse_env%n_spin
      CALL get_qs_env(qs_env, bs_env=bs_env)
      ! Reset the effective Hamiltonian to KS Hamiltonian + G0W0 - reference COHSEX - reference Hartree
      DO j=1,nspin
         ! Sets the imaginary part to zero
         CALL cp_cfm_to_cfm(gwbse_env%ham_reference(j)%matrix, gwbse_env%ham_effective(j)%matrix)
      END DO
      ! Determine the field at current time
      IF (gwbse_env%dft_control%apply_efield_field) THEN
         CALL make_field(gwbse_env%dft_control, gwbse_env%field, gwbse_env%sim_step, gwbse_env%sim_time)
      ELSE
         ! No field
         gwbse_env%field(:) = 0.0_dp
      END IF
      IF ( gwbse_env%hartree_method == rtp_gwbse_hartree_pw ) THEN
         ! TODO : Remove hartree_pw
         CPABORT("PW Hartree no longer supported.")
      END IF
      DO j=1,nspin
         DO k=1,3
            ! Minus sign due to charge of electrons
            CALL cp_fm_to_cfm(msourcer=gwbse_env%moments_field(k)%matrix, mtarget=gwbse_env%ham_workspace(1)%matrix)
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), gwbse_env%ham_effective(j)%matrix,&
                                      CMPLX(gwbse_env%field(k), 0.0, kind=dp), gwbse_env%ham_workspace(1)%matrix)
         END DO
         IF (gwbse_env%ham_reference_type == rtp_gwbse_ham_g0w0) THEN
            ! Add the COH part - so far static but can be dynamic in principle throught the W updates
            CALL get_sigma(gwbse_env, qs_env, gwbse_env%sigma_COH(j),-0.5_dp, gwbse_env%S_inv_fm)
            CALL cp_fm_to_cfm(msourcer=gwbse_env%sigma_COH(j)%matrix, mtarget=gwbse_env%ham_workspace(1)%matrix)
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), gwbse_env%ham_effective(j)%matrix,&
                                      CMPLX(1.0, 0.0, kind=dp), gwbse_env%ham_workspace(1)%matrix)
         END IF
         ! Calculate the (S)EX part - based on provided rho
         ! iGW = - rho W
         CALL get_sigma(gwbse_env, qs_env, gwbse_env%sigma_SEX(j),-1.0_dp, rho(j))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), gwbse_env%ham_effective(j)%matrix,&
                                   CMPLX(1.0, 0.0, kind=dp), gwbse_env%sigma_SEX(j)%matrix)
         ! Calculate Hartree potential
         ! Hartree potential is scaled by number of electrons in each MO - spin degeneracy
         IF ( gwbse_env%hartree_method == rtp_gwbse_hartree_dbt ) THEN
            CALL get_hartree_local(qs_env, gwbse_env, rho(j)%matrix,&
                             gwbse_env%hartree_curr(j)%matrix)
            CALL cp_fm_to_cfm(msourcer=gwbse_env%hartree_curr(j)%matrix, mtarget=gwbse_env%ham_workspace(1)%matrix)
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), gwbse_env%ham_effective(j)%matrix,&
                                      CMPLX(gwbse_env%spin_degeneracy, 0.0, kind=dp), gwbse_env%ham_workspace(1)%matrix)
         ELSEIF (gwbse_env%hartree_method == rtp_gwbse_hartree_pw ) THEN
            ! TODO : Remove PW Hartree
            CPABORT("PW Hartree no longer supported")
         END IF
         ! TODO - enforce hermiticity of the effective Hamiltonian - apply to components instead?
         ! Important components without forced Hermiticity - moments matrix, sigma matrices, Hartree matrix
         ! single particle Ham
         CALL cp_cfm_transpose(gwbse_env%ham_effective(j)%matrix, 'C', gwbse_env%ham_workspace(1)%matrix)
         CALL cp_cfm_scale_and_add(CMPLX(0.5, 0.0, kind=dp), gwbse_env%ham_effective(j)%matrix,&
                                  CMPLX(0.5, 0.0, kind=dp), gwbse_env%ham_workspace(1)%matrix)
      END DO
      CALL timestop(handle)
   END SUBROUTINE update_effective_ham

   ! TODO : Integrate into the main workflow, but also into the Delta kick formalism
! **************************************************************************************************
!> \brief Does the BCH iterative determination of the exponential
!> \param propagator_matrix Matrix X which is to be exponentiated
!> \param target_matrix Matrix Y which the exponential acts upon
!> \param result_matrix Propagated matrix
!> \param workspace Matrices dedicated for work, 4 fm matrices with dimensions of X required
!> \param threshold_opt Optionally, a threshold under which the iteration is considered converged (default 1e-10)
!> \param max_iter_opt Optionally, maximum number of BCH iterations (default 20)
! **************************************************************************************************
   SUBROUTINE bch_propagate(propagator_matrix, target_matrix, result_matrix, workspace, threshold_opt, max_iter_opt)
      ! Array of complex propagator matrix X, such that
      ! the propagated matrix will follow Y' = e^X Y e^(-X), for each spin
      ! effect of e^(-X) is calculated - provide the X on the left hand side
      TYPE(cp_cfm_p_type), DIMENSION(:)                 :: propagator_matrix
      ! Matrix Y to be propagated into matrix Y'
      TYPE(cp_cfm_p_type), DIMENSION(:)                 :: target_matrix
      ! Matrix Y' is stored here on exit
      TYPE(cp_cfm_p_type), DIMENSION(:)                 :: result_matrix, workspace
      ! Threshold for the metric which decides when to truncate the BCH expansion 
      REAL(kind=dp), OPTIONAL                           :: threshold_opt
      INTEGER, OPTIONAL                                 :: max_iter_opt
      CHARACTER(len=*), PARAMETER                       :: routineN="bch_propagate"
      REAL(kind=dp)                                     :: threshold, prefactor, metric
      INTEGER                                           :: max_iter, i, n_spin, n_ao, k,&
                                                           w_stride, handle
      LOGICAL                                           :: converged

      CALL timeset(routineN, handle)

      converged = .FALSE.

      IF (PRESENT(threshold_opt)) THEN
         threshold = threshold_opt
      ELSE
         threshold = 1.0e-10
      END IF

      IF (PRESENT(max_iter_opt)) THEN
         max_iter = max_iter_opt
      ELSE
         max_iter = 20
      END IF

      n_spin = SIZE(target_matrix)
      n_ao = 0
      CALL cp_cfm_get_info(target_matrix(1)%matrix, nrow_global=n_ao)
      w_stride = n_spin

      ! Initiate
      DO i=1,n_spin
         CALL cp_cfm_to_cfm(target_matrix(i)%matrix, result_matrix(i)%matrix)
         CALL cp_cfm_to_cfm(target_matrix(i)%matrix, workspace(i)%matrix)
      END DO

      ! Start the BCH iterations
      ! TODO : Spin mixing terms?
      DO k=1,max_iter
         prefactor = 1.0_dp/REAL(k, kind=dp)
         DO i=1,n_spin
            CALL parallel_gemm("N", "N", n_ao, n_ao, n_ao,&
                               CMPLX(prefactor, 0.0, kind=dp), propagator_matrix(i)%matrix, workspace(i)%matrix,&
                               CMPLX(0.0, 0.0, kind=dp), workspace(i+w_stride)%matrix)
            CALL parallel_gemm("N", "C", n_ao, n_ao, n_ao,&
                               CMPLX(prefactor, 0.0, kind=dp), workspace(i)%matrix, propagator_matrix(i)%matrix,&
                               CMPLX(1.0, 0.0, kind=dp), workspace(i+w_stride)%matrix)
            ! Add to the result
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), result_matrix(i)%matrix,&
                                      CMPLX(1.0, 0.0, kind=dp), workspace(i+w_stride)%matrix)
         END DO
         metric = rho_metric(workspace(w_stride+1:), workspace(1:w_stride), n_spin)
         IF (metric <= threshold) THEN
            converged = .TRUE.
            EXIT
         ELSE
            DO i=1,n_spin
               CALL cp_cfm_to_cfm(workspace(i+w_stride)%matrix, workspace(i)%matrix)
            END DO
         END IF
      END DO
      ! TODO : An informative message, maybe about number of terms and convergence
      IF (.NOT. converged) THEN
         PRINT *, "BCH Metric", metric, "BCH Threshold", threshold
         CPABORT("BCH did not converge")
      END IF

      CALL timestop(handle)
   END SUBROUTINE bch_propagate

! **************************************************************************************************
!> \brief Updates the density in gwbse_env, using the provided exponential
!>        The new density is saved to a different matrix, which enables for comparison of matrices
!> \param gwbse_env Entry point of the calculation - contains current state of variables
!> \param exponential Real and imaginary parts ( + spin) of the exponential propagator 
! **************************************************************************************************
   SUBROUTINE propagate_density(gwbse_env, exponential, rho_old, rho_new)
      TYPE(gwbse_env_type)                               :: gwbse_env
      TYPE(cp_cfm_p_type), DIMENSION(:)                  :: exponential,&
                                                            rho_old,&
                                                            rho_new
      CHARACTER(len=*), PARAMETER                        :: routineN = "propagate_density"
      INTEGER                                            :: j, handle

      CALL timeset(routineN, handle)
      IF (gwbse_env%mat_exp_method == do_taylor .OR. gwbse_env%mat_exp_method == do_exact) THEN
         ! For these methods, exponential is explicitly constructed
         DO j=1,gwbse_env%n_spin
            ! rho * (exp^dagger)
            CALL parallel_gemm("N", "C", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao,&
                               CMPLX(1.0, 0.0, kind=dp), rho_old(j)%matrix, exponential(j)%matrix,&
                               CMPLX(0.0, 0.0, kind=dp), gwbse_env%rho_workspace(1)%matrix)
            ! exp * rho * (exp^dagger)
            CALL parallel_gemm("N", "N", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao,&
                               CMPLX(1.0, 0.0, kind=dp), exponential(j)%matrix, gwbse_env%rho_workspace(1)%matrix,&
                               CMPLX(0.0, 0.0, kind=dp), rho_new(j)%matrix)
         END DO
      ELSE IF (gwbse_env%mat_exp_method == do_bch) THEN
         ! TODO : Specific option for BCH max iter?
         CALL bch_propagate(exponential, rho_old, rho_new, gwbse_env%rho_workspace, threshold_opt=gwbse_env%exp_accuracy,&
                            max_iter_opt=gwbse_env%etrs_max_iter)
      ELSE
         CPABORT("Only BCH, exact and Taylor matrix exponentiation implemented.")
      END IF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Outputs the number of electrons in the system from the density matrix
!> \note  Moments matrix is provided by the gwbse_env, uses rho_workspace(1:3) 
!> \param gwbse_env Entry point - gwbse environment
!> \param rho Density matrix in AO basis
!> \param electron_n_re Real number of electrons
!> \param electron_n_im Imaginary number of electrons, which can arise from numerical non-hermiticity
! **************************************************************************************************
   SUBROUTINE get_electron_number(gwbse_env, rho, electron_n_re, electron_n_im)
      TYPE(gwbse_env_type)                               :: gwbse_env
      TYPE(cp_cfm_p_type), DIMENSION(:)                  :: rho
      REAL(kind=dp), INTENT(OUT)                         :: electron_n_re, electron_n_im
      CHARACTER(len=*), PARAMETER                        :: routineN="get_electron_number"
      COMPLEX(kind=dp)                                   :: electron_n_buffer
      INTEGER                                            :: j

      electron_n_re = 0.0_dp
      electron_n_im = 0.0_dp
      CALL cp_fm_to_cfm(msourcer=gwbse_env%S_fm%matrix, mtarget=gwbse_env%rho_workspace(1)%matrix)
      DO j=1,gwbse_env%n_spin
         CALL cp_cfm_trace(gwbse_env%rho_workspace(1)%matrix, rho(j)%matrix, electron_n_buffer)
         electron_n_re = electron_n_re + REAL(electron_n_buffer, kind=dp)
         electron_n_im = electron_n_im + REAL(AIMAG(electron_n_buffer), kind=dp)
      END DO
      ! Scale by spin degeneracy
      electron_n_re = electron_n_re * gwbse_env%spin_degeneracy
      electron_n_im = electron_n_im * gwbse_env%spin_degeneracy
   END SUBROUTINE get_electron_number
! **************************************************************************************************
!> \brief Outputs the deviation from idempotence of density matrix
!> \note  Moments matrix is provided by the gwbse_env, uses rho_workspace(1:3) 
!> \param gwbse_env Entry point - gwbse environment
!> \param rho Density matrix in AO basis
!> \param electron_n_re Real number of electrons
!> \param electron_n_im Imaginary number of electrons, which can arise from numerical non-hermiticity
! **************************************************************************************************
   SUBROUTINE get_idempotence_deviation(gwbse_env, rho, deviation_metric)
      TYPE(gwbse_env_type)                              :: gwbse_env
      TYPE(cp_cfm_p_type),DIMENSION(:)                  :: rho
      REAL(kind=dp), INTENT(OUT)                        :: deviation_metric
      CHARACTER(len=*), PARAMETER                       :: routineN="get_idempotence_deviation"
      COMPLEX(kind=dp)                                  :: buffer_1, buffer_2
      REAL(kind=dp)                                     :: buffer_dev
      INTEGER                                           :: j

      deviation_metric = 0.0_dp
      buffer_dev = 0.0_dp
      ! First, determine Tr(S * rho_re) + i Tr (S * rho_im)
      ! TODO : Handle better operations where real and complex matrix trace are present
      CALL cp_fm_to_cfm(msourcer=gwbse_env%S_fm%matrix, mtarget=gwbse_env%rho_workspace(1)%matrix)
      DO j=1,gwbse_env%n_spin
         CALL cp_cfm_trace(gwbse_env%rho_workspace(1)%matrix, rho(j)%matrix, buffer_1)
         buffer_dev = buffer_dev + REAL(ABS(buffer_1) * ABS(buffer_1))
      END DO
      ! Now, determine Tr(S * rho_re * S * rho_re) - Tr(S * rho_im * S * rho_im) + 2i Tr(S * rho_re * S * rho_im)
      DO j=1,gwbse_env%n_spin
         ! S * rho
         CALL multiply_fm_cfm("N", "N", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao,&
                            1.0_dp, gwbse_env%S_fm%matrix, rho(j)%matrix,&
                            0.0_dp, gwbse_env%rho_workspace(2)%matrix)
         ! rho * S * rho
         CALL parallel_gemm("N", "N", gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_ao,&
                            CMPLX(1.0, 0.0, kind=dp), rho(j)%matrix, gwbse_env%rho_workspace(2)%matrix,&
                            CMPLX(1.0, 0.0, kind=dp), gwbse_env%rho_workspace(3)%matrix)
         ! Tr (S * rho * S * rho)
         CALL cp_cfm_trace(gwbse_env%rho_workspace(1)%matrix, gwbse_env%rho_workspace(3)%matrix, buffer_2)
         deviation_metric = deviation_metric + REAL(ABS(buffer_2) * ABS(buffer_2))
      END DO
      deviation_metric = SQRT(deviation_metric) - SQRT(buffer_dev)
   END SUBROUTINE get_idempotence_deviation

! **************************************************************************************************
!> \brief Calculates the self-energy by contraction of screened potential, for complex density
!> \note Can be used for both the Coulomb hole part and screened exchange part 
!> \param gwbse_env Quickstep environment data, entry point of the calculation
!> \param qs_env Quickstep environment data, entry point of the calculation
!> \param sigma_cfm Pointer to the self-energy full matrix, which is overwritten by this routine
!> \param greens_cfm Pointer to the Green's function matrix, which is used as input data
!> \author Stepan Marek
!> \date 09.2024
! **************************************************************************************************
SUBROUTINE get_sigma_complex(gwbse_env, qs_env, sigma_cfm, prefactor_opt, greens_cfm)
      TYPE(gwbse_env_type)                               :: gwbse_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_cfm_p_type)                                :: sigma_cfm ! resulting self energy
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: prefactor_opt
      TYPE(cp_cfm_p_type), INTENT(IN)                    :: greens_cfm ! matrix to contract with RI_W
      REAL(kind=dp)                                      :: prefactor

      prefactor = 1.0_dp
      IF (PRESENT(prefactor_opt)) prefactor = prefactor_opt

      ! Carry out the sigma part twice
      ! Real part
      CALL cp_cfm_to_fm(msource=greens_cfm%matrix, mtargetr=gwbse_env%real_workspace(1)%matrix)
      CALL get_sigma(gwbse_env, qs_env, gwbse_env%real_workspace(2), prefactor, gwbse_env%real_workspace(1))
      CALL cp_fm_to_cfm(msourcer=gwbse_env%real_workspace(2)%matrix, mtarget=gwbse_env%ham_workspace(1)%matrix)
      ! Imaginary part
      CALL cp_cfm_to_fm(msource=greens_cfm%matrix, mtargeti=gwbse_env%real_workspace(1)%matrix)
      CALL get_sigma(gwbse_env, qs_env, gwbse_env%real_workspace(2), prefactor, gwbse_env%real_workspace(1))
      CALL cp_fm_to_cfm(msourcei=gwbse_env%real_workspace(2)%matrix, mtarget=sigma_cfm%matrix)
      ! Add the real part
      CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), sigma_cfm%matrix,&
                                CMPLX(1.0, 0.0, kind=dp), gwbse_env%ham_workspace(1)%matrix)
    
END SUBROUTINE get_sigma_complex
! **************************************************************************************************
!> \brief Calculates the self-energy by contraction of screened potential, for complex density
!> \note Can be used for both the Coulomb hole part and screened exchange part 
!> \param gwbse_env Quickstep environment data, entry point of the calculation
!> \param qs_env Quickstep environment data, entry point of the calculation
!> \param sigma_fm Pointer to the self-energy full matrix, which is overwritten by this routine
!> \param greens_fm Pointer to the Green's function matrix, which is used as input data
!> \author Stepan Marek
!> \date 09.2024
! **************************************************************************************************
SUBROUTINE get_sigma_real(gwbse_env, qs_env, sigma_fm, prefactor_opt, greens_fm)
      TYPE(gwbse_env_type)                               :: gwbse_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type)                                 :: sigma_fm ! resulting self energy
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: prefactor_opt
      TYPE(cp_fm_p_type), INTENT(IN)                     :: greens_fm ! matrix to contract with RI_W
      REAL(kind=dp)                                      :: prefactor
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      prefactor = 1.0_dp
      IF (PRESENT(prefactor_opt)) prefactor = prefactor_opt
      
      CALL get_qs_env(qs_env, bs_env=bs_env)

      ! Carry out the sigma part twice
      ! Convert to dbcsr
      CALL copy_fm_to_dbcsr(greens_fm%matrix, gwbse_env%rho_dbcsr)
      CALL get_sigma(gwbse_env, qs_env, sigma_fm, prefactor, gwbse_env%rho_dbcsr)
    
END SUBROUTINE get_sigma_real
! **************************************************************************************************
!> \brief Calculates the self-energy by contraction of screened potential
!> \note Can be used for both the Coulomb hole part and screened exchange part 
!> \param qs_env Quickstep environment data, entry point of the calculation
!> \param greens_fm Pointer to the Green's function matrix, which is used as input data
!> \param sigma_fm Pointer to the self-energy full matrix, which is overwritten by this routine
!> \author Stepan Marek
!> \date 01.2024
! **************************************************************************************************
   SUBROUTINE get_sigma_dbcsr(gwbse_env, qs_env, sigma_fm, prefactor_opt, greens_dbcsr)
      TYPE(gwbse_env_type)                               :: gwbse_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type)                                 :: sigma_fm ! resulting self energy
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: prefactor_opt
      TYPE(dbcsr_type), POINTER                          :: greens_dbcsr
      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_sigma'
      REAL(kind=dp)                                      :: prefactor
      TYPE(dbcsr_type)                                   :: sigma_dbcsr
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

      ! Three-centre integrals are obtained from build_3c_integrals, from qs_tensors
      ! These should use sparcity, while W and Sigma can be full matrices
      ! The summation is carried out by dbt library - dbt_contract in dbt_api
      ! The building of the tensors might be a bit hard, because it requires a lot of parallel information
      ! Probably just use the tensors already present in bs_env? They seem to be mostly work tensors
      ! TODO : Starting with unbounded integrals, then probably need to reduce to a subset
      ! Create by template
      CALL dbt_contract(alpha=1.0_dp,&
                   tensor_1=gwbse_env%screened_dbt,&
                   tensor_2=gwbse_env%t_3c_w,&
                   beta=0.0_dp,&
                   tensor_3=gwbse_env%t_3c_work_RI_AO__AO,&
                   contract_1=[2],notcontract_1=[1], map_1=[1],&
                   contract_2=[1],notcontract_2=[2,3], map_2=[2,3])!,&
                   !filter_eps=bs_env%eps_filter)
      ! t_work1 now contains B^P_(nu beta) = sum _ Q W _ (PQ) (iomega = 0) (Q| nu beta)
      ! Next step is to convert the greens full matrix to dbcsr matrix 
      CALL dbt_copy_matrix_to_tensor(greens_dbcsr, gwbse_env%greens_dbt) 
      ! Then contract it
      ! no scaling applied - this has to be applied externally TODO : Is this a good approach?
      CALL dbt_contract(alpha=1.0_dp,&
              tensor_1=gwbse_env%t_3c_work_RI_AO__AO,&
              tensor_2=gwbse_env%greens_dbt,&
              beta=0.0_dp,&
              tensor_3=gwbse_env%t_3c_work2_RI_AO__AO,&
              contract_1=[2],notcontract_1=[1,3],map_1=[1,3],&
              contract_2=[2],notcontract_2=[1],map_2=[2])
      ! workspace 2 now contains C ^ P _ (mu beta) sum _ nu B ^ P _ (nu beta) g _ (mu nu)
      CALL dbt_contract(alpha=prefactor,&
              tensor_1=gwbse_env%t_3c_w,&
              tensor_2=gwbse_env%t_3c_work2_RI_AO__AO,&
              beta=0.0_dp,&
              tensor_3=gwbse_env%sigma_dbt,&
              contract_1=[1,3],notcontract_1=[2],map_1=[1],&
              contract_2=[1,2],notcontract_2=[3],map_2=[2])!,&
              !filter_eps=bs_env%eps_filter)
      ! Finally, convert the COH tensor to matrix and then to fm matrix
      CALL dbcsr_create(sigma_dbcsr, name="sigma", template=bs_env%mat_ao_ao%matrix)
      CALL dbt_copy_tensor_to_matrix(gwbse_env%sigma_dbt, sigma_dbcsr)
      CALL copy_dbcsr_to_fm(sigma_dbcsr, sigma_fm%matrix)
      CALL dbcsr_release(sigma_dbcsr)
      ! Clear workspaces - saves memory? - TODO : Decide, whether 3c integrals can be stored or not
      CALL dbt_clear(gwbse_env%t_3c_work_RI_AO__AO)
      CALL dbt_clear(gwbse_env%t_3c_work2_RI_AO__AO)
      CALL dbt_clear(gwbse_env%sigma_dbt)
      CALL dbt_clear(gwbse_env%greens_dbt)
      CALL timestop(handle)

   END SUBROUTINE get_sigma_dbcsr
! **************************************************************************************************
!> \brief Creates the RI matrix and populates it with correct values
!> \note Tensor contains Hartree elements in the auxiliary basis
!> \param qs_env Quickstep environment - entry point of calculation
!> \author Stepan Marek
!> \date 01.2024
! **************************************************************************************************
   SUBROUTINE init_hartree(qs_env, v_dbcsr)
      TYPE(qs_environment_type), POINTER, INTENT(IN)     :: qs_env
      TYPE(dbcsr_type), POINTER                          :: v_dbcsr
      CHARACTER(len=*), PARAMETER                        :: routineN = "init_hartree"
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(libint_potential_type)                        :: coulomb_op
      TYPE(cp_fm_type)                                   :: V_fm
      TYPE(cp_fm_type)                                   :: metric_fm
      ! TYPE(cp_fm_type), DIMENSION(:,:), ALLOCATABLE      :: metric_fm
      TYPE(cp_fm_type)                                   :: metric_inv_fm,&
                                                            work_fm
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: V_dbcsr_a,&
                                                            metric_dbcsr
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c

      CALL get_qs_env(qs_env, bs_env=bs_env)

      ! Allocate for bare Hartree term
      ! TODO : k-points
      ALLOCATE(V_dbcsr_a(1))
      ALLOCATE(metric_dbcsr(1))
      CALL dbcsr_create(V_dbcsr_a(1), name="Hartree_dbcsr", template=bs_env%mat_RI_RI%matrix)
      CALL dbcsr_create(metric_dbcsr(1), name="RI_metric_dbcsr", template=bs_env%mat_RI_RI%matrix)

      ! Calculate full coulomb RI basis elements - V _ (PQ) matrix
      NULLIFY(nl_2c)
      CALL build_2c_neighbor_lists(nl_2c, bs_env%basis_set_RI, bs_env%basis_set_RI,&
                                   coulomb_op, "Coulomb_neighbor_2c_list", qs_env,&
                                   sym_ij=.FALSE., molecular=.TRUE.)
      CALL build_2c_integrals(V_dbcsr_a, bs_env%eps_filter, qs_env, nl_2c,&
                              bs_env%basis_set_RI, bs_env%basis_set_RI, coulomb_op,&
                              do_kpoints=.FALSE., regularization_RI=bs_env%regularization_RI)
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
      ! TODO : Spins and k-points
      ! TODO : Better inverse without conversion?
      CALL cp_fm_create(metric_fm, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(metric_fm, 0.0_dp)
      CALL cp_fm_create(metric_inv_fm, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(metric_inv_fm, 0.0_dp)
      CALL cp_fm_create(work_fm, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(work_fm, 0.0_dp)
      CALL copy_dbcsr_to_fm(metric_dbcsr(1), metric_fm)
      CALL cp_fm_invert(metric_fm, metric_inv_fm)
      CALL cp_fm_create(V_fm, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(V_fm, 0.0_dp)
      ! Multiply by the inverse from each side (M^-1 is symmetric)
      CALL cp_dbcsr_sm_fm_multiply(V_dbcsr_a(1), metric_inv_fm,&
                                   work_fm, bs_env%n_RI)
      CALL parallel_gemm("N", "N", bs_env%n_RI, bs_env%n_RI, bs_env%n_RI,&
                         1.0_dp, metric_inv_fm, work_fm, 0.0_dp, V_fm)
      ! Now, create the tensor from the matrix
      ! First, convert full matrix to dbcsr
      ! TODO : No k-points so far
      CALL dbcsr_clear(V_dbcsr_a(1))
      CALL copy_fm_to_dbcsr(V_fm, V_dbcsr_a(1))
      CALL dbcsr_create(v_dbcsr, "Hartree ri", V_dbcsr_a(1))
      CALL dbcsr_copy(v_dbcsr, V_dbcsr_a(1))
      ! Create and copy distinctly, so that unnecessary objects can be destroyed
      ! Destroy all unnecessary matrices
      CALL dbcsr_release(V_dbcsr_a(1))
      CALL dbcsr_release(metric_dbcsr(1))
      DEALLOCATE(V_dbcsr_a)
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
!> \brief Calculates the Hartree matrix in the atomic orbital basis, given a density matrix, in local arrays
!>        Calculates the values for single spin species present in given rho
!> \param qs_env Entry point
!> \param gwbse_env Entry point of GWBSE - uses rho_dbcsr and some complex_workspace
!> \param rho_ao Density matrix in ao basis
!> \param v_ao Overwritten by the Hartree matrix in the atomic orbital basis
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE get_hartree_local(qs_env, gwbse_env, rho_ao, v_ao)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gwbse_env_type)                               :: gwbse_env
      TYPE(cp_cfm_type), INTENT(IN), POINTER             :: rho_ao
      TYPE(cp_fm_type), INTENT(OUT), POINTER             :: v_ao
      CHARACTER(len=*), PARAMETER                        :: routineN="get_hartree_local"
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iterator_matrix
      INTEGER                                            :: i, j, k, n, nblocks, ind_1, ind_2, row_offset, col_offset,&
                                                            row_size, col_size, j_n_AO, k_n_AO, i_n_RI, n_RI,&
                                                            ri_offset, ind_i, handle
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: Pvector, Qvector
      REAL(kind=dp), DIMENSION(:,:), POINTER             :: block_matrix
      REAL(kind=dp), DIMENSION(:,:,:), POINTER           :: int_3c
      TYPE(dbcsr_type)                                   :: v_dbcsr_ao

      ! Very ineffective first implementation - calculate all 3cs an all ranks
      ! Importantly - dbcsr blocks are ordered by atoms - i.e. ethene with 6 atoms will have 6x6 block structure
      ! Number of basis states on each basis set is known is post_scf_bandstructure env
      ! TODO : Only calculate integrals which are needed for a given rank, which are given by structure
      !        of the density matrix

      CALL timeset(routineN, handle)
      ! Get the relevant environments
      CALL get_qs_env(qs_env, bs_env=bs_env, para_env=para_env)
      n_RI = gwbse_env%n_RI
      int_3c => gwbse_env%int_3c_array

      ! Allocate the Q and Pvector on each rank
      ALLOCATE(Qvector(gwbse_env%n_RI), source=0.0_dp)
      ALLOCATE(Pvector(gwbse_env%n_RI), source=0.0_dp)

      ! First step - analyze the structure of copied dbcsr matrix on all ranks
      CALL dbcsr_clear(gwbse_env%rho_dbcsr)
      ! Only the real part of the density matrix contributes
      CALL cp_cfm_to_fm(msource=rho_ao, mtargetr=gwbse_env%real_workspace(1)%matrix)
      CALL copy_fm_to_dbcsr(gwbse_env%real_workspace(1)%matrix, gwbse_env%rho_dbcsr)
      nblocks = dbcsr_get_num_blocks(gwbse_env%rho_dbcsr)
      CALL dbcsr_iterator_start(iterator_matrix, gwbse_env%rho_dbcsr)
      DO n=1,nblocks
         CALL dbcsr_iterator_next_block(iterator_matrix, ind_1, ind_2, block_matrix,&
                 row_offset=row_offset, col_offset=col_offset, row_size=row_size, col_size=col_size)
         ! Now we have a block corresponding to a single atom pair
         j_n_AO = bs_env%sizes_AO(ind_1)
         k_n_AO = bs_env%sizes_AO(ind_2)
         ! For each atom pair, we need to get contributions from all RI atoms
         ! The allocations are as follows
         ! Now, get the relevant 3c integrals
         ri_offset = 0
         DO ind_i=1,bs_env%n_atom
            i_n_RI = bs_env%sizes_RI(ind_i)
            !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j,k) &
            !$OMP SHARED(Qvector, int_3c, block_matrix,i_n_RI,j_n_AO,k_n_AO, ri_offset, ind_1, ind_2, ind_i, bs_env)
            DO i=1,i_n_RI
               DO j=1,j_n_AO
                  DO k=1,k_n_AO
                     ! Different OMP threads write to different places in memory - should be safe from data race
                     Qvector(ri_offset + i) = Qvector(ri_offset + i) + int_3c(j+bs_env%i_ao_start_from_atom(ind_1)-1,&
                                                                              k+bs_env%i_ao_start_from_atom(ind_2)-1,&
                                                                              i+bs_env%i_RI_start_from_atom(ind_i)-1)&
                                                                              * block_matrix(j,k)
                  END DO
               END DO
            END DO
            !$OMP END PARALLEL DO
            ri_offset = ri_offset + i_n_RI
         END DO
      END DO
      CALL dbcsr_iterator_stop(iterator_matrix)
      ! Now, each rank has contributions from D_jk within its scope
      ! Need to sum over different ranks to get the total vector on all ranks
      CALL para_env%sum(Qvector)
      ! TODO : Maybe search for matrix_vector multiply/trace scale for the V^PQ B_Q operation
      ! Once this is done, Pvector is current on all ranks
      ! Continue with V_PQ summation
      nblocks = dbcsr_get_num_blocks(gwbse_env%v_dbcsr%matrix)
      CALL dbcsr_iterator_start(iterator_matrix, gwbse_env%v_dbcsr%matrix)
      DO n=1,nblocks
         ! TODO : Try OMP parallelisation over different blocks - expect many more available speedup for large systems
         CALL dbcsr_iterator_next_block(iterator_matrix, ind_1, ind_2, block_matrix, &
                                        row_offset=row_offset, col_offset=col_offset, row_size=row_size, col_size=col_size)
         ! TODO : Better names for RI
         j_n_AO = bs_env%sizes_RI(ind_1)
         k_n_AO = bs_env%sizes_RI(ind_2)
         ! The allocations are as follows
         !$OMP PARALLEL DO DEFAULT(none) PRIVATE(j,k) &
         !$OMP SHARED(block_matrix, Pvector, Qvector,j_n_AO,k_n_AO,row_offset,col_offset)
         DO j=1,j_n_AO
            DO k=1,k_n_AO
               Pvector(j+row_offset-1) = Pvector(j+row_offset-1) + block_matrix(j,k) * Qvector(k+col_offset-1)
            END DO
         END DO
         !$OMP END PARALLEL DO
      END DO
      CALL dbcsr_iterator_stop(iterator_matrix)
      ! Again, make sure that the P vector is present on all ranks
      CALL para_env%sum(Pvector)
      ! Now, for the final trick, iterate over local blocks of v_dbcsr_ao to get the Hartree as dbcsr, then convert to fm
      CALL dbcsr_create(v_dbcsr_ao, "Hartree ao", gwbse_env%rho_dbcsr)
      CALL copy_fm_to_dbcsr(v_ao, v_dbcsr_ao)
      nblocks = dbcsr_get_num_blocks(v_dbcsr_ao)
      CALL dbcsr_iterator_start(iterator_matrix, v_dbcsr_ao)
      DO n=1,nblocks
         CALL dbcsr_iterator_next_block(iterator_matrix, ind_1, ind_2, block_matrix,&
                                        row_offset=row_offset, col_offset=col_offset)
         j_n_AO = bs_env%sizes_AO(ind_1)
         k_n_AO = bs_env%sizes_AO(ind_2)
         ri_offset = 0
         block_matrix(:,:) = 0.0_dp
         DO ind_i=1,bs_env%n_atom
            i_n_RI = bs_env%sizes_RI(ind_i)
            ! TODO : In principle, can parallelize over both directions here
            !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j,k) &
            !$OMP SHARED(block_matrix, Pvector, int_3c, i_n_RI, j_n_AO, k_n_AO, ri_offset, ind_1, ind_2, ind_i, bs_env)
            DO j=1,j_n_AO
               DO k=1,k_n_AO
                  ! Add all the summed up values
                  DO i=1,i_n_RI
                     block_matrix(j, k) = block_matrix(j, k) + int_3c(j+bs_env%i_ao_start_from_atom(ind_1)-1,&
                                                                      k+bs_env%i_ao_start_from_atom(ind_2)-1,&
                                                                      i+bs_env%i_RI_start_from_atom(ind_i)-1) * Pvector(i+ri_offset)
                  END DO
               END DO
            END DO
            !$OMP END PARALLEL DO
            ri_offset = ri_offset + i_n_RI
         END DO
      END DO
      CALL dbcsr_iterator_stop(iterator_matrix)
      ! Since P vector was present on all the ranks, v_dbcsr_ao has the complete Hartree result
      CALL copy_dbcsr_to_fm(v_dbcsr_ao, v_ao)
      CALL dbcsr_release(v_dbcsr_ao)
      DEALLOCATE(Qvector)
      DEALLOCATE(Pvector)

      CALL timestop(handle)
   END SUBROUTINE get_hartree_local
   SUBROUTINE cp_cfm_gexp(amatrix, bmatrix, exponential, eig_scale_opt, work_opt)
      ! TODO : Do interface for real matrices
      TYPE(cp_cfm_type), INTENT(IN)                      :: amatrix
      TYPE(cp_cfm_type), INTENT(IN)                      :: bmatrix
      TYPE(cp_cfm_type)                                  :: exponential
      COMPLEX(kind=dp), INTENT(IN), OPTIONAL             :: eig_scale_opt
      TYPE(cp_cfm_p_type), DIMENSION(:), OPTIONAL        :: work_opt
      CHARACTER(len=*), PARAMETER                        :: routineN="cp_cfm_exp"
      COMPLEX(kind=dp)                                   :: eig_scale
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: eigenvalues
      COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE        :: expvalues
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER         :: work
      LOGICAL                                            :: deallocate_work
      INTEGER                                            :: nrow,i

      ! Argument parsing and sanity checks
      IF (PRESENT(eig_scale_opt)) THEN
         eig_scale = eig_scale_opt
      ELSE
         eig_scale = CMPLX(1.0, 0.0, kind=dp)
      END IF

      NULLIFY(work)
      ALLOCATE(work(4))
      IF (PRESENT(work_opt) .AND. SIZE(work_opt) >= 4) THEN
         deallocate_work = .FALSE.
         DO i=1,4
            work(i)%matrix => work_opt(i)%matrix
         END DO
      ELSE
         deallocate_work = .TRUE.
         ! Allocate the work storage on the fly
         DO i=1,4
            NULLIFY(work(i)%matrix)
            ALLOCATE(work(i)%matrix)
            CALL cp_cfm_create(work(i)%matrix, amatrix%matrix_struct)
         END DO
      END IF

      nrow = amatrix%matrix_struct%nrow_global

      ALLOCATE(eigenvalues(nrow))
      ALLOCATE(expvalues(nrow))

      ! Do not change the amatrix and bmatrix - need to copy them first
      CALL cp_cfm_to_cfm(amatrix, work(1)%matrix)
      CALL cp_cfm_to_cfm(bmatrix, work(2)%matrix)

      ! Solve the generalized eigenvalue equation
      CALL cp_cfm_geeig(work(1)%matrix, work(2)%matrix, work(3)%matrix, eigenvalues, work(4)%matrix)

      ! Scale and exponentiate the eigenvalues
      ! TODO : Check threading/vectorization
      expvalues(:) = EXP(eigenvalues(:) * eig_scale)

      ! Copy eigenvectors to column scale them
      CALL cp_cfm_to_cfm(work(3)%matrix, work(1)%matrix)
      ! C * exp(eig)
      CALL cp_cfm_column_scale(work(1)%matrix, expvalues)

      ! Carry out the remaining operations
      ! C * exp(eig) * C^H
      CALL parallel_gemm("N", "C", nrow, nrow, nrow,&
                         CMPLX(1.0, 0.0, kind=dp), work(1)%matrix, work(3)%matrix,&
                         CMPLX(0.0, 0.0, kind=dp), work(2)%matrix)
      CALL cp_cfm_to_cfm(bmatrix, work(3)%matrix)
      ! C * exp(eig) * C^H * S
      CALL parallel_gemm("N", "N", nrow, nrow, nrow,&
                         CMPLX(1.0, 0.0, kind=dp), work(2)%matrix, work(3)%matrix,&
                         CMPLX(0.0, 0.0, kind=dp), exponential)

      ! Deallocate work storage if necessary
      IF (deallocate_work) THEN
         DO i=1,4
            CALL cp_cfm_release(work(i)%matrix)
            DEALLOCATE(work(i)%matrix)
         END DO
      END IF
      DEALLOCATE(work)

      DEALLOCATE(eigenvalues)
      DEALLOCATE(expvalues)
   END SUBROUTINE cp_cfm_gexp
END MODULE rt_gwbse

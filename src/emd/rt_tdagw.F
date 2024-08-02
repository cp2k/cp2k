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
!>        and where the initial guess is calculated/loaded
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
                                              cp_fm_set_element,&
                                              cp_fm_get_info,&
                                              cp_fm_write_formatted
   USE cp_cfm_types,                    ONLY: cp_cfm_type,&
                                              cp_fm_to_cfm,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_create,&
                                              cp_cfm_release
   USE kinds,                           ONLY: dp
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
                                              cp_cfm_transpose,&
                                              cp_cfm_norm,&
                                              cp_cfm_column_scale
   USE cp_fm_diag,                      ONLY: cp_fm_geeig
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_moments,                      ONLY: build_local_moment_matrix,&
                                              build_berry_moment_matrix
   USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
   USE efield_utils,                    ONLY: make_field
   USE rt_propagator_init,              ONLY: rt_initialize_rho_from_mos
   USE rt_propagation_methods,          ONLY: s_matrices_create
   USE gw_large_cell_Gamma,             ONLY: compute_3c_integrals
   USE gw_integrals,                    ONLY: build_3c_integral_block
   USE matrix_exp,                      ONLY: taylor_full_complex
   USE cp_log_handling,                 ONLY: cp_logger_type,&
                                              cp_get_default_logger,&
                                              cp_logger_get_unit_nr
   USE cp_output_handling,              ONLY: cp_print_key_unit_nr
   USE message_passing,                 ONLY: mp_para_env_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE cell_types,                      ONLY: cell_type
   USE input_constants,                 ONLY: rtp_tdagw_ham_ks,&
                                              rtp_tdagw_ham_g0w0,&
                                              rtp_tdagw_hartree_dbt,&
                                              rtp_tdagw_hartree_pw,&
                                              do_taylor,&
                                              do_bch,&
                                              do_exact
   USE rt_tdagw_types,                  ONLY: tdagw_env_type,&
                                              create_tdagw_env,&
                                              release_tdagw_env
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_integrate_potential_product,  ONLY: integrate_v_rspace
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE physcon,                         ONLY: femtoseconds
   USE mathconstants,                   ONLY: twopi

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
      TYPE(cell_type), POINTER                           :: cell
      TYPE(section_vals_type), POINTER                   :: input, rtp_section
      ! TODO : Maybe complex?
      INTEGER                                            :: i,j,re,im, k, handle, custom_exp_method
      LOGICAL                                            :: converged
      REAL(kind=dp)                                      :: metric, enum_re, enum_im,&
                                                            idempotence_dev, a_metric_1, a_metric_2,&
                                                            intensity
      REAL(kind=dp), DIMENSION(3)                        :: kvec
      REAL(kind=dp), DIMENSION(:), POINTER               :: gon_eigvals
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(rtp_control_type), POINTER                    :: rtp_control
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s

      CALL timeset(routineN, handle)
      ! ********************** ALLOCATIONS **********************
      CALL create_tdagw_env(tdagw_env, qs_env)
      ! Get the logger unit number, only print if we are the source process
      CALL get_qs_env(qs_env, bs_env=bs_env, input=input, rtp=rtp, cell=cell, dft_control=dft_control, matrix_s=matrix_s)
      rtp_control => dft_control%rtp_control
      ! TODO : Print a nicer header
      IF (tdagw_env%unit_nr > 0) WRITE(tdagw_env%unit_nr, *) "=== Starting time-dependent adiabatic GW real time propagation ==="

      ! ********************** EXTRACT PREVIOUS DATA **********************

      ! Calculate the initial Hartree interaction matrix
      ! TODO : k-points and the possible divergence
      IF ( tdagw_env%hartree_method == rtp_tdagw_hartree_dbt ) THEN
         ! Calculate Coulomb RI elements, necessary for Hartree calculation
         CALL init_hartree(qs_env, tdagw_env%hartree_dbt, tdagw_env%v_dbcsr%matrix)
         ! Add the Hartree to the screened_dbt tensor - now W = V + W^c
         CALL dbcsr_add(tdagw_env%w_dbcsr%matrix, tdagw_env%v_dbcsr%matrix, 1.0_dp, 1.0_dp)
         CALL dbt_copy_matrix_to_tensor(tdagw_env%w_dbcsr%matrix, tdagw_env%screened_dbt)
         ! TODO : So far only uncoupled spins
         @:SPIN_DO(j, re, im, tdagw_env%n_spin)
            CALL get_hartree_local(qs_env, tdagw_env, tdagw_env%rho(re)%matrix, tdagw_env%hartree_orig(re)%matrix)
            CALL get_hartree_local(qs_env, tdagw_env, tdagw_env%rho(im)%matrix, tdagw_env%hartree_orig(im)%matrix)
            ! Scaling by spin degeneracy
            CALL cp_fm_scale(tdagw_env%spin_degeneracy, tdagw_env%hartree_orig(re)%matrix)
            CALL cp_fm_scale(tdagw_env%spin_degeneracy, tdagw_env%hartree_orig(im)%matrix)
         END DO
      ELSEIF ( tdagw_env%hartree_method == rtp_tdagw_hartree_pw ) THEN
         ! TODO : Remove pw hartree
         CPABORT("TDAGW Hartree PW no longer supported.")
      ELSE
         CPABORT("Only DBT and PW methods for determination of Hartree term implemented.")
      END IF

      ! Calculate the COHSEX starting energies
      IF (tdagw_env%ham_single_particle_type == rtp_tdagw_ham_g0w0) THEN
         ! Subtract the v_xc from COH part of the self-energy, as V_xc is also not updated during the timestepping
         @:SPIN_DO(j, re, im, tdagw_env%n_spin)
            CALL get_sigma(tdagw_env, qs_env, tdagw_env%sigma_COH(re), -0.5_dp, greens_dbcsr_opt=tdagw_env%S_inv)
            ! TODO : Enable subtraction of COH? So far, COH does not enter, as per the reference paper
            CALL get_sigma(tdagw_env, qs_env, tdagw_env%sigma_SEX(re), -1.0_dp, greens_fm_opt=tdagw_env%rho(re))
            CALL get_sigma(tdagw_env, qs_env, tdagw_env%sigma_SEX(im), -1.0_dp, greens_fm_opt=tdagw_env%rho(im))
            ! TODO : Allow for subtraction of COH
            CALL cp_fm_scale_and_add(0.0_dp, tdagw_env%sigma_COH(re)%matrix, -1.0_dp, tdagw_env%sigma_SEX(re)%matrix)
            CALL cp_fm_scale_and_add(0.0_dp, tdagw_env%sigma_COH(im)%matrix, -1.0_dp, tdagw_env%sigma_SEX(im)%matrix)
         END DO
      ELSE
         ! KS Hamiltonian - use static xc only
         @:SPIN_DO(j, re, im, tdagw_env%n_spin)
            ! CALL cp_fm_scale_and_add(0.0_dp, tdagw_env%sigma_COH(re)%matrix, 1.0_dp, bs_env%fm_V_xc_Gamma(j))
            ! Instead, use time-dependent Fock exchange
            CALL get_sigma(tdagw_env, qs_env, tdagw_env%sigma_SEX(re), -1.0_dp, greens_fm_opt=tdagw_env%rho(re))
            CALL get_sigma(tdagw_env, qs_env, tdagw_env%sigma_SEX(im), -1.0_dp, greens_fm_opt=tdagw_env%rho(im))
            ! CALL get_sigma_local(tdagw_env, qs_env, tdagw_env%sigma_SEX(re), -1.0_dp, greens_fm_opt=tdagw_env%rho(re))
            ! CALL get_sigma_local(tdagw_env, qs_env, tdagw_env%sigma_SEX(im), -1.0_dp, greens_fm_opt=tdagw_env%rho(im))
            CALL cp_fm_scale_and_add(0.0_dp, tdagw_env%sigma_COH(re)%matrix, -1.0_dp, tdagw_env%sigma_SEX(re)%matrix)
            CALL cp_fm_scale_and_add(0.0_dp, tdagw_env%sigma_COH(im)%matrix, -1.0_dp, tdagw_env%sigma_SEX(im)%matrix)
         END DO
      END IF

      ! Output 0 time moments
      ! Setup the files
      rtp_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION")
      tdagw_env%sim_time = REAL(tdagw_env%sim_start, dp)*rtp%dt
      CALL output_moments(tdagw_env, tdagw_env%rho, bs_env%fm_mo_coeff_Gamma, rtp_section)

      ! Debug : Alternative application of delta kick - copying the formalism used in TDDFT
      ! IF (rtp_control%apply_delta_pulse) THEN
      !    ! Get the dbcsr exp matrices
      !    NULLIFY(cosmat)
      !    NULLIFY(sinmat)
      !    ALLOCATE(cosmat)
      !    ALLOCATE(sinmat)
      !    CALL dbcsr_create(cosmat, template=matrix_s(1)%matrix)
      !    CALL dbcsr_create(sinmat, template=matrix_s(1)%matrix)

      !    @:SPIN_DO(j, re, im, tdagw_env%n_spin)
      !       CALL dbcsr_copy(cosmat, matrix_s(j)%matrix)
      !       CALL dbcsr_copy(sinmat, matrix_s(j)%matrix)
      !       kvec(:) = cell%h_inv(1, :)*rtp_control%delta_pulse_direction(1) + &
      !                 cell%h_inv(2, :)*rtp_control%delta_pulse_direction(2) + &
      !                 cell%h_inv(3, :)*rtp_control%delta_pulse_direction(3)
      !       kvec = kvec*twopi*rtp_control%delta_pulse_scale
      !       CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, -kvec)

      !       CALL copy_dbcsr_to_fm(cosmat, tdagw_env%ham_workspace(re)%matrix)
      !       CALL copy_dbcsr_to_fm(sinmat, tdagw_env%ham_workspace(im)%matrix)
      !       CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
      !               1.0_dp, tdagw_env%S_inv_fm%matrix, tdagw_env%ham_workspace(re)%matrix,&
      !               0.0_dp, tdagw_env%ham_effective(re)%matrix)
      !       CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
      !               1.0_dp, tdagw_env%S_inv_fm%matrix, tdagw_env%ham_workspace(im)%matrix,&
      !               0.0_dp, tdagw_env%ham_effective(im)%matrix)
      !    END DO

      !    custom_exp_method = tdagw_env%mat_exp_method
      !    tdagw_env%mat_exp_method = do_exact
      !    CALL propagate_density(tdagw_env, tdagw_env%ham_effective, tdagw_env%rho, tdagw_env%rho_new)
      !    tdagw_env%mat_exp_method = custom_exp_method
      !    metric = rho_metric(tdagw_env%rho_new, tdagw_env%rho, tdagw_env%n_spin)
      !    IF (tdagw_env%unit_nr > 0) WRITE(tdagw_env%unit_nr, *) "Metric difference after delta kick", metric
      !    ! Copy the new density to the old density
      !    @:SPIN_DO(j, re, im, tdagw_env%n_spin)
      !       CALL cp_fm_to_fm(tdagw_env%rho_new(re)%matrix, tdagw_env%rho(re)%matrix)
      !       CALL cp_fm_to_fm(tdagw_env%rho_new(im)%matrix, tdagw_env%rho(im)%matrix)
      !    END DO

      !    CALL dbcsr_release(cosmat)
      !    CALL dbcsr_release(sinmat)
      !    DEALLOCATE(cosmat)
      !    DEALLOCATE(sinmat)
      ! END IF
      ! TODO : Check spectrum of custom delta kick
      IF (rtp_control%apply_delta_pulse) THEN
         ! Extra minus sign for the charge of electrons
         intensity = -rtp_control%delta_pulse_scale
         metric = 0.0_dp
         kvec(:) = cell%h_inv(1, :)*rtp_control%delta_pulse_direction(1) + &
                   cell%h_inv(2, :)*rtp_control%delta_pulse_direction(2) + &
                   cell%h_inv(3, :)*rtp_control%delta_pulse_direction(3)
         kvec = kvec*twopi*rtp_control%delta_pulse_scale
         @:SPIN_DO(j, re, im, tdagw_env%n_spin)
            CALL cp_fm_set_all(tdagw_env%ham_effective(im)%matrix, 0.0_dp)
            DO k=1,3
               ! CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(im)%matrix,&
               !         REAL(rtp_control%delta_pulse_direction(k), kind=dp), tdagw_env%moments_field(k)%matrix)
               CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(im)%matrix,&
                       kvec(k), tdagw_env%moments_field(k)%matrix)
            END DO
            ! enforce hermiticity of the effective Hamiltonian - apply to components instead?
            CALL cp_fm_transpose(tdagw_env%ham_effective(im)%matrix, tdagw_env%rho_workspace(1)%matrix)
            CALL cp_fm_scale_and_add(0.5_dp, tdagw_env%ham_effective(im)%matrix,&
                                    -0.5_dp, tdagw_env%rho_workspace(1)%matrix)
            IF (tdagw_env%mat_exp_method == do_bch) THEN
               ! Multiply by the S_inv matrix - in the classic ordering
               CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                                  intensity, tdagw_env%S_inv_fm%matrix, tdagw_env%ham_effective(im)%matrix,&
                                  0.0_dp, tdagw_env%ham_effective(re)%matrix)
               CALL cp_fm_to_fm(tdagw_env%ham_effective(re)%matrix, tdagw_env%ham_effective(im)%matrix)
               CALL cp_fm_set_all(tdagw_env%ham_effective(re)%matrix, 0.0_dp)
            ELSE IF (tdagw_env%mat_exp_method == do_exact) THEN
               NULLIFY(gon_eigvals)
               ALLOCATE(gon_eigvals(tdagw_env%n_ao))
               CALL cp_fm_to_fm(tdagw_env%S_fm%matrix, tdagw_env%rho_workspace(1)%matrix)
               CALL cp_fm_geeig(tdagw_env%ham_effective(im)%matrix, tdagw_env%rho_workspace(1)%matrix,&
                       tdagw_env%ham_workspace(1)%matrix, tdagw_env%real_eigvals, tdagw_env%ham_workspace(2)%matrix)
               IF (tdagw_env%unit_nr > 0) WRITE(tdagw_env%unit_nr, *) "Real eigvals", tdagw_env%real_eigvals
               ! Determine the exponential eigenvalues
               ! Start with real part
               gon_eigvals(:) = COS(intensity * tdagw_env%real_eigvals(:))
               IF (tdagw_env%unit_nr > 0) WRITE(tdagw_env%unit_nr, *) "Cos eigvals", gon_eigvals
               ! Scale the real eigenvectors
               CALL cp_fm_to_fm(tdagw_env%ham_workspace(1)%matrix, tdagw_env%ham_workspace(2)%matrix)
               CALL cp_fm_column_scale(tdagw_env%ham_workspace(2)%matrix, gon_eigvals)
               CALL parallel_gemm("N", "T", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                                  1.0_dp, tdagw_env%ham_workspace(2)%matrix, tdagw_env%ham_workspace(1)%matrix,&
                                  0.0_dp, tdagw_env%ham_effective(im)%matrix)
               CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                                  1.0_dp, tdagw_env%ham_effective(im)%matrix, tdagw_env%S_fm%matrix,&
                                  0.0_dp, tdagw_env%ham_effective(re)%matrix)
               ! Continue for the imaginary part
               gon_eigvals(:) = -SIN(intensity * tdagw_env%real_eigvals(:))
               IF (tdagw_env%unit_nr > 0) WRITE(tdagw_env%unit_nr, *) "Sin eigvals", gon_eigvals
               CALL cp_fm_to_fm(tdagw_env%ham_workspace(1)%matrix, tdagw_env%ham_workspace(2)%matrix)
               CALL cp_fm_column_scale(tdagw_env%ham_workspace(2)%matrix, gon_eigvals)
               CALL parallel_gemm("N", "T", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                                  1.0_dp, tdagw_env%ham_workspace(2)%matrix, tdagw_env%ham_workspace(1)%matrix,&
                                  0.0_dp, tdagw_env%ham_effective(im)%matrix)
               CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                                  1.0_dp, tdagw_env%ham_effective(im)%matrix, tdagw_env%S_fm%matrix,&
                                  0.0_dp, tdagw_env%ham_workspace(2)%matrix)
               CALL cp_fm_to_fm(tdagw_env%ham_workspace(2)%matrix, tdagw_env%ham_effective(im)%matrix)
               PRINT *, "Exponential from TDaGW"
               CALL cp_fm_write_formatted(tdagw_env%ham_effective(re)%matrix, tdagw_env%unit_nr)
               CALL cp_fm_write_formatted(tdagw_env%ham_effective(im)%matrix, tdagw_env%unit_nr)
               DEALLOCATE(gon_eigvals)
            END IF
         END DO
         IF (tdagw_env%mat_exp_method == do_bch) THEN
            CALL bch_propagate(tdagw_env%ham_effective, tdagw_env%rho, tdagw_env%rho_new, tdagw_env%rho_workspace,&
                    threshold_opt=tdagw_env%exp_accuracy)
         ELSE IF (tdagw_env%mat_exp_method == do_exact) THEN
            ! Exponential is constructed - apply it to the density matrix
            CALL propagate_density(tdagw_env, tdagw_env%ham_effective, tdagw_env%rho, tdagw_env%rho_new)
         ELSE
            CPABORT("Only exact and BCH exponential schemes adopted for delta kick")
         END IF
         metric = rho_metric(tdagw_env%rho_new, tdagw_env%rho, tdagw_env%n_spin)
         IF (tdagw_env%unit_nr > 0) WRITE(tdagw_env%unit_nr, *) "Metric difference after delta kick", metric
         ! Copy the new density to the old density
         @:SPIN_DO(j, re, im, tdagw_env%n_spin)
            CALL cp_fm_to_fm(tdagw_env%rho_new(re)%matrix, tdagw_env%rho(re)%matrix)
            CALL cp_fm_to_fm(tdagw_env%rho_new(im)%matrix, tdagw_env%rho(im)%matrix)
         END DO
      END IF

      ! ********************** Start the time loop **********************
      DO i = tdagw_env%sim_start, tdagw_env%sim_nsteps

         ! Update the simulation time
         tdagw_env%sim_time = REAL(i, dp)*rtp%dt
         tdagw_env%sim_step = i
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
         ! Print the updated field
         CALL output_field(tdagw_env, rtp_section)
         ! Convert the effective hamiltonian into the exponential
         CALL ham_to_exp(tdagw_env)
         ! Propagate the density to mid-point
         CALL propagate_density(tdagw_env, tdagw_env%ham_workspace, tdagw_env%rho, tdagw_env%rho_M)
         CALL get_electron_number(tdagw_env, tdagw_env%rho_M, enum_re, enum_im)
         IF (tdagw_env%unit_nr > 0) WRITE (tdagw_env%unit_nr, *) "Electron number RHO_M re : ", enum_re,&
                                                                 "Electron number RHO_M im : ", enum_im
         ! In initial iteration, copy rho_M to rho_new_last
         @:SPIN_DO(j,re,im, tdagw_env%n_spin)
            CALL cp_fm_to_fm(tdagw_env%rho_M(re)%matrix, tdagw_env%rho_new_last(re)%matrix)
            CALL cp_fm_to_fm(tdagw_env%rho_M(im)%matrix, tdagw_env%rho_new_last(im)%matrix)
         END DO
         ! *********** Start the self-consistent loop ************************
         tdagw_env%sim_time = REAL(i+1, dp)*rtp%dt
         tdagw_env%sim_step = i+1
         converged = .FALSE.
         DO k=1,tdagw_env%etrs_max_iter
            CALL update_effective_ham(tdagw_env, qs_env, tdagw_env%rho_new_last)
            CALL ham_to_exp(tdagw_env)
            CALL propagate_density(tdagw_env, tdagw_env%ham_workspace, tdagw_env%rho_M, tdagw_env%rho_new)
            ! *** Self-consistency check ***
            metric = rho_metric(tdagw_env%rho_new, tdagw_env%rho_new_last, tdagw_env%n_spin)
            IF (tdagw_env%unit_nr > 0) WRITE(tdagw_env%unit_nr, *) "Self-consistent iteration : ", k,&
                                                                   "convergence metric : ", metric
            IF (metric < tdagw_env%etrs_threshold) THEN
               converged = .TRUE.
               EXIT
            ELSE
               ! Copy rho_new to rho_new_last
               @:SPIN_DO(j, re, im, tdagw_env%n_spin)
                  ! Leaving for free convergence
                  CALL cp_fm_to_fm(tdagw_env%rho_new(re)%matrix, tdagw_env%rho_new_last(re)%matrix)
                  CALL cp_fm_to_fm(tdagw_env%rho_new(im)%matrix, tdagw_env%rho_new_last(im)%matrix)
               END DO
            END IF
         END DO
         CALL get_electron_number(tdagw_env, tdagw_env%rho_new, enum_re, enum_im)
         CALL get_idempotence_deviation(tdagw_env, tdagw_env%rho_new, idempotence_dev)
         a_metric_1 = antiherm_metric(tdagw_env%sigma_SEX(re),tdagw_env%sigma_SEX(im))
         a_metric_2 = antiherm_metric(tdagw_env%hartree_curr(re),tdagw_env%hartree_curr(im))
         IF (tdagw_env%unit_nr > 0) WRITE(tdagw_env%unit_nr, *) "Sim step : ", tdagw_env%sim_step, "Convergence : ", metric,&
                                                                "Threshold : ", tdagw_env%etrs_threshold,&
                                                                "Electron number", enum_re, enum_im, "SC iterations : ", k,&
                                                                "Idempotence deviation : ", idempotence_dev,&
                                                                "Antiherm. metric Sigma : ", a_metric_1,& 
                                                                "Antiherm. metric Hartree : ", a_metric_2
         CPASSERT(converged)
         @:SPIN_DO(j, re, im, tdagw_env%n_spin)
            CALL cp_fm_to_fm(tdagw_env%rho_new(re)%matrix, tdagw_env%rho(re)%matrix)
            CALL cp_fm_to_fm(tdagw_env%rho_new(im)%matrix, tdagw_env%rho(im)%matrix)
         END DO
         ! TODO : k-points
         CALL output_mos_rho(tdagw_env, tdagw_env%rho, bs_env%fm_mo_coeff_Gamma, rtp_section)
         CALL output_mos_cohsex(tdagw_env, tdagw_env%sigma_SEX, bs_env%fm_mo_coeff_Gamma, rtp_section)
         CALL output_moments(tdagw_env, tdagw_env%rho, bs_env%fm_mo_coeff_Gamma, rtp_section)
      END DO
      ! ********************** End the time loop **********************

      ! ********************** DEALLOCATIONS **********************
      CALL release_tdagw_env(tdagw_env, qs_env)
      ! Deallocate the neighbour list that is not deallocated in gw anymore
      IF (ASSOCIATED(bs_env%nl_3c%ij_list)) CALL neighbor_list_3c_destroy(bs_env%nl_3c)

      CALL timestop(handle)
   END SUBROUTINE run_propagation_gw

! **************************************************************************************************
!> \brief Determines the metric for the density matrix, used for convergence criterion
!> \param rho_new Array of new density matrices (one for each spin index)
!> \param rho_old Array of old density matrices (one for each spin index)
!> \param nspin Number of spin indices
! **************************************************************************************************
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
      DO j=1,nspin
         metric = metric + partial_metric(j)*partial_metric(j)
      END DO
      metric = SQRT(metric)
      ! Deallocate workspace
      CALL cp_cfm_release(workspace(1))
      CALL cp_cfm_release(workspace(2))
      DEALLOCATE(partial_metric)
   END FUNCTION

! **************************************************************************************************
!> \brief Determines the metric of the antihermitian part of the matrix
!> \param real_fm Real part of the full matrix
!> \param imag_fm Imaginary part of the full matrix
! **************************************************************************************************
   FUNCTION antiherm_metric(real_fm, imag_fm) RESULT (metric)
      TYPE(cp_fm_p_type), INTENT(IN)                    :: real_fm,&
                                                           imag_fm
      REAL(kind=dp)                                     :: metric
      TYPE(cp_cfm_type)                                 :: complex_fm,&
                                                           conjugate_fm
      COMPLEX(kind=dp)                                  :: complex_one

      CALL cp_cfm_create(complex_fm, real_fm%matrix%matrix_struct)
      CALL cp_cfm_create(conjugate_fm, real_fm%matrix%matrix_struct)
      ! Get the complex and complex conjugate matrix
      CALL cp_fm_to_cfm(real_fm%matrix, imag_fm%matrix, complex_fm)
      CALL cp_cfm_transpose(complex_fm, "C", conjugate_fm)
      ! Subtract these, and get the metric
      complex_one = 1.0
      CALL cp_cfm_scale_and_add(complex_one, complex_fm, -complex_one, conjugate_fm)
      metric = cp_cfm_norm(complex_fm, "M")
      CALL cp_cfm_release(complex_fm)
      CALL cp_cfm_release(conjugate_fm)
   END FUNCTION

! **************************************************************************************************
!> \brief For Taylor and Exact exp_method, calculates the matrix exponential of the
!>        effective Hamiltonian. For BCH, calculates just the effective Hamiltonian. For other methods,
!>        aborts the execution, as they are not implemented yet.
!> \param tdagw_env Entry point of the calculation. Uses rho_workspace for Taylor and BCH. For exact,
!>                  uses complex_workspace, complex_ham, complex_s, real_eigvals and exp_eigvals.
!>                  Results are stored in ham_workspace.
! **************************************************************************************************
   SUBROUTINE ham_to_exp(tdagw_env)
      TYPE(tdagw_env_type)                               :: tdagw_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "ham_to_exp"
      INTEGER                                            :: j, re, im, handle
      CALL timeset(routineN, handle)
      @:SPIN_DO(j, re, im, tdagw_env%n_spin)
         IF (tdagw_env%mat_exp_method == do_taylor .OR. tdagw_env%mat_exp_method == do_bch) THEN
            ! In Taylor and BCH, we first evaluate the entire exponent and then evaluate exponential in series
            ! In order to produce correct result, need to remultiply by inverse overlap matrix
            CALL cp_dbcsr_sm_fm_multiply(tdagw_env%S_inv, tdagw_env%ham_effective(re)%matrix,&
                                         tdagw_env%rho_workspace(1)%matrix, tdagw_env%n_ao)
            CALL cp_dbcsr_sm_fm_multiply(tdagw_env%S_inv, tdagw_env%ham_effective(im)%matrix,&
                                         tdagw_env%rho_workspace(2)%matrix, tdagw_env%n_ao)

            ! The evolution of density matrix is derived from the right multiplying term
            ! Imaginary part of the exponent = -real part of the matrix
            CALL cp_fm_scale(-tdagw_env%sim_dt/2, tdagw_env%rho_workspace(1)%matrix)
            ! Real part of the exponent = imag part of the matrix
            CALL cp_fm_scale(tdagw_env%sim_dt/2, tdagw_env%rho_workspace(2)%matrix)
            ! TODO : Understand which input options set the orders
            IF (tdagw_env%mat_exp_method == do_taylor) THEN
               CALL taylor_full_complex([tdagw_env%ham_workspace(re)%matrix, tdagw_env%ham_workspace(im)%matrix],&
                                        tdagw_env%rho_workspace(2)%matrix, tdagw_env%rho_workspace(1)%matrix,&
                                        tdagw_env%orders(1,j), tdagw_env%orders(2,j))
            ELSE
               ! In BCH, exponential is not calculated explicitly, but the propagation is solved in series
               CALL cp_fm_to_fm(tdagw_env%rho_workspace(2)%matrix, tdagw_env%ham_workspace(re)%matrix)
               CALL cp_fm_to_fm(tdagw_env%rho_workspace(1)%matrix, tdagw_env%ham_workspace(im)%matrix)
            END IF
         ELSE IF (tdagw_env%mat_exp_method == do_exact) THEN
            ! In exact exponentiation, we solve the generalized eigenvalue problem by transforming to complex matrices
            CALL cp_fm_to_cfm(tdagw_env%ham_effective(re)%matrix,tdagw_env%ham_effective(im)%matrix,&
                              tdagw_env%complex_workspace(3)%matrix)
            ! Solve the generalized eigenvalue problem, with fm overlap (also complex)
            ! Call the diagonalisation routine
            ! cp_cfm_geeig is done in place - all matrices are overwritten
            CALL cp_cfm_to_cfm(tdagw_env%complex_s(j)%matrix, tdagw_env%complex_workspace(4)%matrix)
            CALL cp_cfm_geeig(tdagw_env%complex_workspace(3)%matrix, tdagw_env%complex_workspace(4)%matrix,&
                              tdagw_env%complex_workspace(1)%matrix, tdagw_env%real_eigvals,&
                              tdagw_env%complex_workspace(2)%matrix)
            ! With these eigenvalues determined, construct the exponential
            tdagw_env%exp_eigvals(:) = EXP(tdagw_env%real_eigvals(:) * CMPLX(0.0, -tdagw_env%sim_dt/2, kind=dp))
            ! Now, get the exponential
            CALL cp_cfm_to_cfm(tdagw_env%complex_workspace(1)%matrix, tdagw_env%complex_workspace(2)%matrix)
            CALL cp_cfm_column_scale(tdagw_env%complex_workspace(2)%matrix, tdagw_env%exp_eigvals)
            CALL parallel_gemm("N", "C", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                               CMPLX(1.0_dp, 0.0_dp, kind=dp), tdagw_env%complex_workspace(2)%matrix,&
                               tdagw_env%complex_workspace(1)%matrix,&
                               CMPLX(0.0_dp, 0.0_dp, kind=dp), tdagw_env%complex_workspace(3)%matrix)
            CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                               CMPLX(1.0_dp, 0.0_dp, kind=dp), tdagw_env%complex_workspace(3)%matrix,&
                               tdagw_env%complex_s(j)%matrix,&
                               CMPLX(0.0_dp, 0.0_dp, kind=dp), tdagw_env%complex_workspace(2)%matrix)
            ! Copy to separated real and imaginary parts
            CALL cp_cfm_to_fm(tdagw_env%complex_workspace(2)%matrix, tdagw_env%ham_workspace(re)%matrix,&
                              tdagw_env%ham_workspace(im)%matrix)
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
   SUBROUTINE update_effective_ham(tdagw_env, qs_env, rho)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: rho
      CHARACTER(len=*), PARAMETER                        :: routineN = "update_effective_ham"
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(kind=dp)                                      :: a_metric
      INTEGER                                            :: k, j, re, im, nspin, handle

      CALL timeset(routineN, handle)
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
      IF (tdagw_env%dft_control%apply_efield_field) THEN
         ! TODO : Understand and implement rtp%istep - that should mainly be used for restart calculations?
         CALL make_field(tdagw_env%dft_control, tdagw_env%field, tdagw_env%sim_step, tdagw_env%sim_time)
      ELSE
         ! No field
         tdagw_env%field(:) = 0.0_dp
      END IF
      IF ( tdagw_env%hartree_method == rtp_tdagw_hartree_pw ) THEN
         ! TODO : Remove hartree_pw
         CPABORT("PW Hartree no longer supported.")
      END IF
      @:SPIN_DO(j, re, im, nspin)
         DO k=1,3
            ! Minus sign due to charge of electrons
            CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(re)%matrix,&
                                    tdagw_env%field(k), tdagw_env%moments_field(k)%matrix)
         END DO
         ! Add the COH part
         ! TODO : Separate the updated and reference part of the Hamiltonian more clearly
         CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(re)%matrix,&
                                  1.0_dp, tdagw_env%sigma_COH(re)%matrix)
         CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(im)%matrix,&
                                  1.0_dp, tdagw_env%sigma_COH(im)%matrix)
         ! Calculate the SEX part - based on provided rho
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
         ! Hartree potential is scaled by number of electrons in each MO - spin degeneracy
         ! TODO : Spin (in)dependence check
         IF ( tdagw_env%hartree_method == rtp_tdagw_hartree_dbt ) THEN
            CALL get_hartree_local(qs_env, tdagw_env, rho(re)%matrix,&
                             tdagw_env%hartree_curr(re)%matrix)
            CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(re)%matrix,&
                                     tdagw_env%spin_degeneracy, tdagw_env%hartree_curr(re)%matrix)
            CALL get_hartree_local(qs_env, tdagw_env, rho(im)%matrix,&
                             tdagw_env%hartree_curr(im)%matrix)
            CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(im)%matrix,&
                                     tdagw_env%spin_degeneracy, tdagw_env%hartree_curr(im)%matrix)
         ELSEIF (tdagw_env%hartree_method == rtp_tdagw_hartree_pw ) THEN
            ! TODO : Remove PW Hartree
            CPABORT("PW Hartree no longer supported")
         END IF
         ! Subtract original Hartree, which is purely real
         CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(re)%matrix,&
                                 -1.0_dp, tdagw_env%hartree_orig(re)%matrix)
         CALL cp_fm_scale_and_add(1.0_dp, tdagw_env%ham_effective(im)%matrix,&
                                 -1.0_dp, tdagw_env%hartree_orig(im)%matrix)
         ! TODO - enforce hermiticity of the effective Hamiltonian - apply to components instead?
         CALL cp_fm_transpose(tdagw_env%ham_effective(re)%matrix, tdagw_env%rho_workspace(1)%matrix)
         CALL cp_fm_scale_and_add(0.5_dp, tdagw_env%ham_effective(re)%matrix,&
                                  0.5_dp, tdagw_env%rho_workspace(1)%matrix)
         CALL cp_fm_transpose(tdagw_env%ham_effective(im)%matrix, tdagw_env%rho_workspace(1)%matrix)
         CALL cp_fm_scale_and_add(0.5_dp, tdagw_env%ham_effective(im)%matrix,&
                                 -0.5_dp, tdagw_env%rho_workspace(1)%matrix)
      END DO
      CALL timestop(handle)
   END SUBROUTINE

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
      ! Array of real and imaginary parts of the propagator matrix X, such that
      ! the propagated matrix will follow Y' = e^X Y e^(-X), for each spin
      ! effect of e^(-X) is calculated - provide the X on the left hand side
      TYPE(cp_fm_p_type), DIMENSION(:)                  :: propagator_matrix
      ! Matrix Y to be propagated into matrix Y'
      TYPE(cp_fm_p_type), DIMENSION(:)                  :: target_matrix
      ! Matrix Y' is stored here on exit
      TYPE(cp_fm_p_type), DIMENSION(:)                  :: result_matrix, workspace
      ! Threshold for the metric which decides when to truncate the BCH expansion 
      REAL(kind=dp), OPTIONAL                           :: threshold_opt
      INTEGER, OPTIONAL                                 :: max_iter_opt
      CHARACTER(len=*), PARAMETER                       :: routineN="bch_propagate"
      REAL(kind=dp)                                     :: threshold, prefactor, metric
      INTEGER                                           :: max_iter, i, re, im, n_spin, n_ao, k,&
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

      n_spin = SIZE(target_matrix) / 2
      n_ao = 0
      CALL cp_fm_get_info(target_matrix(1)%matrix, nrow_global=n_ao)
      w_stride = 2 * n_spin

      ! Initiate
      @:SPIN_DO(i, re, im, n_spin)
         CALL cp_fm_to_fm(target_matrix(re)%matrix, result_matrix(re)%matrix)
         CALL cp_fm_to_fm(target_matrix(im)%matrix, result_matrix(im)%matrix)
         CALL cp_fm_to_fm(target_matrix(re)%matrix, workspace(re)%matrix)
         CALL cp_fm_to_fm(target_matrix(im)%matrix, workspace(im)%matrix)
      END DO

      ! Start the BCH iterations
      ! TODO : Spin mixing terms?
      DO k=1,max_iter
         prefactor = 1.0_dp/REAL(k, kind=dp)
         @:SPIN_DO(i, re, im, n_spin)
            ! Real part
            CALL parallel_gemm("N", "N", n_ao, n_ao, n_ao,&
                               prefactor, propagator_matrix(re)%matrix, workspace(re)%matrix,&
                               0.0_dp, workspace(w_stride+re)%matrix)
            CALL parallel_gemm("N", "N", n_ao, n_ao, n_ao,&
                              -prefactor, propagator_matrix(im)%matrix, workspace(im)%matrix,&
                               1.0_dp, workspace(w_stride+re)%matrix)
            CALL parallel_gemm("N", "T", n_ao, n_ao, n_ao,&
                               prefactor, workspace(re)%matrix, propagator_matrix(re)%matrix,&
                               1.0_dp, workspace(w_stride+re)%matrix)
            CALL parallel_gemm("N", "T", n_ao, n_ao, n_ao,&
                               prefactor, workspace(im)%matrix, propagator_matrix(im)%matrix,&
                               1.0_dp, workspace(w_stride+re)%matrix)
            ! Imag part
            CALL parallel_gemm("N", "N", n_ao, n_ao, n_ao,&
                               prefactor, propagator_matrix(re)%matrix, workspace(im)%matrix,&
                               0.0_dp, workspace(w_stride+im)%matrix)
            CALL parallel_gemm("N", "N", n_ao, n_ao, n_ao,&
                               prefactor, propagator_matrix(im)%matrix, workspace(re)%matrix,&
                               1.0_dp, workspace(w_stride+im)%matrix)
            CALL parallel_gemm("N", "T", n_ao, n_ao, n_ao,&
                              -prefactor, workspace(re)%matrix, propagator_matrix(im)%matrix,&
                               1.0_dp, workspace(w_stride+im)%matrix)
            CALL parallel_gemm("N", "T", n_ao, n_ao, n_ao,&
                               prefactor, workspace(im)%matrix, propagator_matrix(re)%matrix,&
                               1.0_dp, workspace(w_stride+im)%matrix)
            ! Add to the result
            CALL cp_fm_scale_and_add(1.0_dp, result_matrix(re)%matrix, 1.0_dp, workspace(re+w_stride)%matrix)
            CALL cp_fm_scale_and_add(1.0_dp, result_matrix(im)%matrix, 1.0_dp, workspace(im+w_stride)%matrix)
            ! Check for convergence, TODO : Apply to all spins, instead of a per-spin basis?
            CALL cp_fm_set_all(workspace(re)%matrix, 0.0_dp)
            CALL cp_fm_set_all(workspace(im)%matrix, 0.0_dp)
         END DO
         metric = rho_metric(workspace(w_stride+1:), workspace(1:w_stride), n_spin)
         IF (metric <= threshold) THEN
            converged = .TRUE.
            EXIT
         ELSE
            @:SPIN_DO(i, re, im, n_spin)
               CALL cp_fm_to_fm(workspace(re+w_stride)%matrix, workspace(re)%matrix)
               CALL cp_fm_to_fm(workspace(im+w_stride)%matrix, workspace(im)%matrix)
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
!> \brief Updates the density in tdagw_env, using the provided exponential
!>        The new density is saved to a different matrix, which enables for comparison of matrices
!> \param tdagw_env Entry point of the calculation - contains current state of variables
!> \param exponential Real and imaginary parts ( + spin) of the exponential propagator 
! **************************************************************************************************
   SUBROUTINE propagate_density(tdagw_env, exponential, rho_old, rho_new)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: exponential,&
                                                            rho_old,&
                                                            rho_new
      CHARACTER(len=*), PARAMETER                        :: routineN = "propagate_density"
      INTEGER                                            :: j, re, im, k, handle
      REAL(kind=dp)                                      :: convergence,&
                                                            prefactor
      LOGICAL                                            :: converged

      CALL timeset(routineN, handle)
      ! TODO : Consult/think about using complex routines and types
      ! (ham_exp(re)^T - i * ham_exp(im)^T) * (rho_full(re) + i * rho_full(im)) * (ham_exp(re) + i * ham_exp(im)) = 
      ! + ham_exp(re)^T * (rho_full(re) * ham_exp(re) - rho_full(im) * ham_exp(im)) | Term 1
      ! + ham_exp(im)^T * (rho_full(re) * ham_exp(im) + rho_full(im) * ham_exp(re)) | Term 2
      ! - i * ham_exp(im)^T * (rho_full(re) * ham_exp(re) - rho_full(im) * ham_exp(im)) | Term 1
      ! + i * ham_exp(re)^T * (rho_full(re) * ham_exp(im) + rho_full(im) * ham_exp(re)) | Term 2
      ! workspace 2 - real part, workspace 3 - imag part
      ! Real part
      ! Term 1 - right multiplication first
      IF (tdagw_env%mat_exp_method == do_taylor .OR. tdagw_env%mat_exp_method == do_exact) THEN
         @:SPIN_DO(j, re, im, tdagw_env%n_spin)
            CALL parallel_gemm("N","T", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                               1.0_dp, rho_old(re)%matrix, exponential(re)%matrix,&
                               0.0_dp, tdagw_env%rho_workspace(1)%matrix)
            CALL parallel_gemm("N","T", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                               1.0_dp, rho_old(im)%matrix, exponential(im)%matrix,&
                               1.0_dp, tdagw_env%rho_workspace(1)%matrix)
            ! Now term 1 left multiplication real part
            CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                               1.0_dp, exponential(re)%matrix, tdagw_env%rho_workspace(1)%matrix,&
                               0.0_dp, tdagw_env%rho_workspace(2)%matrix)
            ! Now term 1 left multiplication imag part
            CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                               1.0_dp, exponential(im)%matrix, tdagw_env%rho_workspace(1)%matrix,&
                               0.0_dp, tdagw_env%rho_workspace(3)%matrix)
            ! Repeat the same procedure for the second term
            ! Term 2 - right multiplication first
            CALL parallel_gemm("N","T", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                               1.0_dp, rho_old(re)%matrix, exponential(im)%matrix,&
                               0.0_dp, tdagw_env%rho_workspace(1)%matrix)
            CALL parallel_gemm("N","T",tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                              -1.0_dp, rho_old(im)%matrix, exponential(re)%matrix,&
                               1.0_dp, tdagw_env%rho_workspace(1)%matrix)
            ! Now term 2 left multiplication real part
            CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                               1.0_dp, exponential(im)%matrix, tdagw_env%rho_workspace(1)%matrix,&
                               1.0_dp, tdagw_env%rho_workspace(2)%matrix)
            ! Now term 2 left multiplication imag part
            CALL parallel_gemm("N","N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                              -1.0_dp, exponential(re)%matrix, tdagw_env%rho_workspace(1)%matrix,&
                               1.0_dp, tdagw_env%rho_workspace(3)%matrix)
            ! Finally, copy the results to rho - not used anymore TODO : Unless there is spin coupling
            CALL cp_fm_to_fm(tdagw_env%rho_workspace(2)%matrix, rho_new(re)%matrix)
            CALL cp_fm_to_fm(tdagw_env%rho_workspace(3)%matrix, rho_new(im)%matrix)
         END DO
      ELSE IF (tdagw_env%mat_exp_method == do_bch) THEN
         ! TODO : Specific option for BCH max iter?
         CALL bch_propagate(exponential, rho_old, rho_new, tdagw_env%rho_workspace, threshold_opt=tdagw_env%exp_accuracy,&
                            max_iter_opt=tdagw_env%etrs_max_iter)
      ELSE
         CPABORT("Only BCH, exact and Taylor matrix exponentiation implemented.")
      END IF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Outputs rho in basis of MOs to a file
!> \param tdagw_env Entry point - tdagw environment
!> \param rho Density matrix in AO basis
!> \param mos Spin dependent molecular orbitals
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE output_mos_rho(tdagw_env, rho, mos, rtp_section)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: rho
      TYPE(cp_fm_type), DIMENSION(:)                     :: mos
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
      END DO
   END SUBROUTINE
! **************************************************************************************************
!> \brief Prints the current field components into a file provided by input
!> \param tdagw_env Entry point - tdagw environment
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE output_field(tdagw_env, rtp_section)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(section_vals_type) , POINTER                  :: rtp_section
      CHARACTER(len=*), PARAMETER                        :: routineN="output_field"
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: field_unit
      
      ! Get logger
      logger => cp_get_default_logger()
      ! Get file descriptor
      field_unit = cp_print_key_unit_nr(logger, rtp_section, "PRINT%FIELD", extension=".field")
      ! If the file descriptor is non-zero, output field
      IF (field_unit /= -1) THEN
         WRITE(field_unit, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3)') tdagw_env%sim_time*femtoseconds,&
               tdagw_env%field(1), tdagw_env%field(2), tdagw_env%field(3)
      END IF
      
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
      END DO
   END SUBROUTINE
! **************************************************************************************************
!> \brief Outputs the number of electrons in the system from the density matrix
!> \note  Moments matrix is provided by the tdagw_env, uses rho_workspace(1:3) 
!> \param tdagw_env Entry point - tdagw environment
!> \param rho Density matrix in AO basis
!> \param electron_n_re Real number of electrons
!> \param electron_n_im Imaginary number of electrons, which can arise from numerical non-hermiticity
! **************************************************************************************************
   SUBROUTINE get_electron_number(tdagw_env, rho, electron_n_re, electron_n_im)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: rho
      REAL(kind=dp), INTENT(OUT)                         :: electron_n_re, electron_n_im
      CHARACTER(len=*), PARAMETER                        :: routineN="get_electron_number"
      REAL(kind=dp)                                      :: electron_n_buffer
      INTEGER                                            :: i, re, im

      electron_n_re = 0.0_dp
      electron_n_im = 0.0_dp
      @:SPIN_DO(i, re, im, tdagw_env%n_spin)
         CALL cp_fm_trace(tdagw_env%S_fm%matrix, rho(re)%matrix, electron_n_buffer)
         electron_n_re = electron_n_re + electron_n_buffer
         CALL cp_fm_trace(tdagw_env%S_fm%matrix, rho(im)%matrix, electron_n_buffer)
         electron_n_im = electron_n_im + electron_n_buffer
      END DO
      ! Scale by spin degeneracy
      electron_n_re = electron_n_re * tdagw_env%spin_degeneracy
      electron_n_im = electron_n_im * tdagw_env%spin_degeneracy
   END SUBROUTINE get_electron_number
! **************************************************************************************************
!> \brief Outputs the deviation from idempotence of density matrix
!> \note  Moments matrix is provided by the tdagw_env, uses rho_workspace(1:3) 
!> \param tdagw_env Entry point - tdagw environment
!> \param rho Density matrix in AO basis
!> \param electron_n_re Real number of electrons
!> \param electron_n_im Imaginary number of electrons, which can arise from numerical non-hermiticity
! **************************************************************************************************
   SUBROUTINE get_idempotence_deviation(tdagw_env, rho, deviation_metric)
      TYPE(tdagw_env_type)                              :: tdagw_env
      TYPE(cp_fm_p_type),DIMENSION(*)                   :: rho
      REAL(kind=dp), INTENT(OUT)                        :: deviation_metric
      CHARACTER(len=*), PARAMETER                       :: routineN="get_idempotence_deviation"
      REAL(kind=dp)                                     :: buffer_1, buffer_2, buffer_3, buffer_dev
      INTEGER                                           :: re,im,i

      deviation_metric = 0.0_dp
      buffer_dev = 0.0_dp
      ! First, determine Tr(S * rho_re) + i Tr (S * rho_im)
      @:SPIN_DO(i, re, im, tdagw_env%n_spin)
         CALL cp_fm_trace(tdagw_env%S_fm%matrix, rho(re)%matrix, buffer_1)
         CALL cp_fm_trace(tdagw_env%S_fm%matrix, rho(im)%matrix, buffer_2)
         buffer_dev = buffer_dev + buffer_1 * buffer_1 + buffer_2 * buffer_2
      END DO
      ! Now, determine Tr(S * rho_re * S * rho_re) - Tr(S * rho_im * S * rho_im) + 2i Tr(S * rho_re * S * rho_im)
      @:SPIN_DO(i, re, im, tdagw_env%n_spin)
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, tdagw_env%S_fm%matrix, rho(re)%matrix,&
                            0.0_dp, tdagw_env%rho_workspace(1)%matrix)
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, tdagw_env%S_fm%matrix, rho(im)%matrix,&
                            0.0_dp, tdagw_env%rho_workspace(2)%matrix)
         ! rho_re * S * rho_re
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, rho(re)%matrix, tdagw_env%rho_workspace(1)%matrix,&
                            0.0_dp, tdagw_env%rho_workspace(3)%matrix)
         CALL cp_fm_trace(tdagw_env%S_fm%matrix, tdagw_env%rho_workspace(3)%matrix, buffer_1)
         ! rho_im * S * rho_im
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, rho(im)%matrix, tdagw_env%rho_workspace(2)%matrix,&
                            0.0_dp, tdagw_env%rho_workspace(3)%matrix)
         CALL cp_fm_trace(tdagw_env%S_fm%matrix, tdagw_env%rho_workspace(3)%matrix, buffer_2)
         buffer_3 = buffer_1 - buffer_2
         ! rho_im * S * rho_re
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, rho(im)%matrix, tdagw_env%rho_workspace(1)%matrix,&
                            0.0_dp, tdagw_env%rho_workspace(3)%matrix)
         CALL cp_fm_trace(tdagw_env%S_fm%matrix, tdagw_env%rho_workspace(3)%matrix, buffer_1)
         deviation_metric = deviation_metric + buffer_3 * buffer_3 + 4 * buffer_1 * buffer_1
      END DO
      deviation_metric = SQRT(deviation_metric) - SQRT(buffer_dev)
   END SUBROUTINE get_idempotence_deviation
! **************************************************************************************************
!> \brief Outputs COHSEX energy in basis of MOs to a file
!> \note Differs from the density matrix output due to differing representations of the operators
!> \param tdagw_env Entry point - tdagw environment
!> \param cohsex cohsex matrix in AO basis, covariant representation
!> \param mos Spin dependent molecular orbitals
!> \param rtp_section RTP input section
! **************************************************************************************************
   SUBROUTINE output_mos_cohsex(tdagw_env, cohsex, mos, rtp_section, custom_unit)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: cohsex
      TYPE(cp_fm_type), DIMENSION(:)                     :: mos
      TYPE(section_vals_type), POINTER                   :: rtp_section
      INTEGER, DIMENSION(2), OPTIONAL                    :: custom_unit
      CHARACTER(len=*), PARAMETER                        :: routineN = "output_mos_cohsex"
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: j, re, im,&
                                                            rho_unit_re, rho_unit_im
      CHARACTER(len=21), DIMENSION(4)                    :: file_labels
      
      file_labels(1) = "_COHSEX_SPIN_A_RE.dat"
      file_labels(2) = "_COHSEX_SPIN_A_IM.dat"
      file_labels(3) = "_COHSEX_SPIN_B_RE.dat"
      file_labels(4) = "_COHSEX_SPIN_B_IM.dat"
      logger => cp_get_default_logger()
      ! Start by multiplying the current density by MOS
      ! TODO : Avoid multiplications if not required
      @:SPIN_DO(j,re, im, tdagw_env%n_spin)
         IF(.NOT. PRESENT(custom_unit)) THEN
            rho_unit_re = cp_print_key_unit_nr(logger, rtp_section, "PRINT%DENSITY_MATRIX", extension=file_labels(re))
            rho_unit_im = cp_print_key_unit_nr(logger, rtp_section, "PRINT%DENSITY_MATRIX", extension=file_labels(im))
         ELSE
            rho_unit_re = custom_unit(1)
            rho_unit_im = custom_unit(2)
         END IF
         ! Transform the density matrix into molecular orbitals basis and print it out
         ! TODO : Number of molecular orbitals and atomic orbitals might not be the same
         ! TODO : Spins in S_fm
         ! Real part
         ! C^T * cohsex
         CALL parallel_gemm("T", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, mos(j), cohsex(re)%matrix, 0.0_dp, tdagw_env%rho_workspace(1)%matrix)
         ! C^T * cohsex * C
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, tdagw_env%rho_workspace(1)%matrix, mos(j),&
                            0.0_dp, tdagw_env%rho_workspace(2)%matrix)
         CALL cp_fm_write_formatted(tdagw_env%rho_workspace(2)%matrix, rho_unit_re)
         ! Imag part
         ! C^T * cohsex
         CALL parallel_gemm("T", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, mos(j), cohsex(im)%matrix, 0.0_dp, tdagw_env%rho_workspace(1)%matrix)
         ! C^T * cohsex * C
         CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, tdagw_env%rho_workspace(1)%matrix, mos(j), 0.0_dp, tdagw_env%rho_workspace(2)%matrix)
         CALL cp_fm_write_formatted(tdagw_env%rho_workspace(2)%matrix, rho_unit_im)
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
         NULLIFY(greens_dbcsr)
         ALLOCATE(greens_dbcsr)
         CALL dbcsr_create(greens_dbcsr, name="greens", template=bs_env%mat_ao_ao%matrix)
         CALL copy_fm_to_dbcsr(greens_fm_opt%matrix, greens_dbcsr)
      ELSE
         ! TODO : Warn and exit
         IF (tdagw_env%unit_nr > 0) WRITE(tdagw_env%unit_nr, *) "Either dbcsr or fm greens function has to be provided"
         CPASSERT(.FALSE.)
      END IF

      ! Three-centre integrals are obtained from build_3c_integrals, from qs_tensors
      ! These should use sparcity, while W and Sigma can be full matrices
      ! The summation is carried out by dbt library - dbt_contract in dbt_api
      ! The building of the tensors might be a bit hard, because it requires a lot of parallel information
      ! Probably just use the tensors already present in bs_env? They seem to be mostly work tensors
      ! TODO : Starting with unbounded integrals, then probably need to reduce to a subset
      ! Create by template
      CALL dbt_contract(alpha=1.0_dp,&
                   tensor_1=tdagw_env%screened_dbt,&
                   tensor_2=tdagw_env%t_3c_w,&
                   beta=0.0_dp,&
                   tensor_3=tdagw_env%t_3c_work_RI_AO__AO,&
                   contract_1=[2],notcontract_1=[1], map_1=[1],&
                   contract_2=[1],notcontract_2=[2,3], map_2=[2,3])!,&
                   !filter_eps=bs_env%eps_filter)
      ! t_work1 now contains B^P_(nu beta) = sum _ Q W _ (PQ) (iomega = 0) (Q| nu beta)
      ! Next step is to convert the greens full matrix to dbcsr matrix 
      CALL dbt_copy_matrix_to_tensor(greens_dbcsr, tdagw_env%greens_dbt) 
      ! Then contract it
      ! no scaling applied - this has to be applied externally TODO : Is this a good approach?
      CALL dbt_contract(alpha=1.0_dp,&
              tensor_1=tdagw_env%t_3c_work_RI_AO__AO,&
              tensor_2=tdagw_env%greens_dbt,&
              beta=0.0_dp,&
              tensor_3=tdagw_env%t_3c_work2_RI_AO__AO,&
              contract_1=[2],notcontract_1=[1,3],map_1=[1,3],&
              contract_2=[2],notcontract_2=[1],map_2=[2])
      ! workspace 2 now contains C ^ P _ (mu beta) sum _ nu B ^ P _ (nu beta) g _ (mu nu)
      CALL dbt_contract(alpha=prefactor,&
              tensor_1=tdagw_env%t_3c_w,&
              tensor_2=tdagw_env%t_3c_work2_RI_AO__AO,&
              beta=0.0_dp,&
              tensor_3=tdagw_env%sigma_dbt,&
              contract_1=[1,3],notcontract_1=[2],map_1=[1],&
              contract_2=[1,2],notcontract_2=[3],map_2=[2])!,&
              !filter_eps=bs_env%eps_filter)
      ! Finally, convert the COH tensor to matrix and then to fm matrix
      CALL dbcsr_create(sigma_dbcsr, name="sigma", template=bs_env%mat_ao_ao%matrix)
      CALL dbt_copy_tensor_to_matrix(tdagw_env%sigma_dbt, sigma_dbcsr)
      CALL copy_dbcsr_to_fm(sigma_dbcsr, sigma_fm%matrix)
      CALL dbcsr_release(sigma_dbcsr)
      IF (PRESENT(greens_fm_opt)) THEN
         CALL dbcsr_release(greens_dbcsr)
         DEALLOCATE(greens_dbcsr)
      END IF
      ! Clear workspaces - saves memory? - TODO : Decide, whether 3c integrals can be stored or not
      CALL dbt_clear(tdagw_env%t_3c_work_RI_AO__AO)
      CALL dbt_clear(tdagw_env%t_3c_work2_RI_AO__AO)
      CALL dbt_clear(tdagw_env%sigma_dbt)
      CALL dbt_clear(tdagw_env%greens_dbt)
      CALL timestop(handle)

   END SUBROUTINE
! **************************************************************************************************
!> \brief Creates the V_dbt tensor and populates it with correct values
!> \note Tensor contains Hartree elements in the auxiliary basis
!> \param qs_env Quickstep environment - entry point of calculation
!> \param V_dbt Tensor representation of the Hartree matrix in the auxiliary basis
!> \author Stepan Marek
!> \date 01.2024
! **************************************************************************************************
   SUBROUTINE init_hartree(qs_env, V_dbt, v_dbcsr)
      TYPE(qs_environment_type), POINTER, INTENT(IN)     :: qs_env
      TYPE(dbt_type), POINTER, INTENT(OUT)               :: V_dbt
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
      ! CALL dbcsr_print(V_dbcsr_a(1))
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
      ! CALL cp_fm_invert(metric_fm(1,1), metric_inv_fm)
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
      CALL dbt_create(V_dbcsr_a(1), V_dbt, name="Hartree_dbt")
      CALL dbt_copy_matrix_to_tensor(V_dbcsr_a(1), V_dbt)
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
      LOGICAL                                            :: found
      INTEGER, DIMENSION(2)                              :: block_shape
      ! Pair of indices for current block
      INTEGER, DIMENSION(3)                              :: index_tensor
      REAL(kind=dp), DIMENSION(:,:), POINTER             :: block_matrix
      REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE       :: block_tensor
      INTEGER                                            :: j,number_of_blocks,backup
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: indices_block_matrix_1,&
                                                            indices_block_matrix_2,&
                                                            indices_block_matrix_custom,&
                                                            my_ploc
      INTEGER                                            :: thread_lower,&
                                                            thread_upper,&
                                                            thread_number,&
                                                            per_thread,&
                                                            per_this_thread,&
                                                            extra_threads,&
                                                            thread_id,&
                                                            handle

      CALL timeset(routineN, handle)
      
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
      ALLOCATE(my_ploc(3), source=0)
      ! TODO : Think about parallelisation
      ! TODO : This returns the same value for every process - wrong pgrid creation?
      CALL dbt_get_info(tensor, my_ploc)
      CALL dbcsr_iterator_start(iterator_matrix, matrix_desym_p)
      DO j=1,number_of_blocks
         CALL dbcsr_iterator_next_block(iterator_matrix, indices_block_matrix_1(j), indices_block_matrix_2(j), block_matrix)
      END DO
      DEALLOCATE(my_ploc)
      CALL dbcsr_iterator_stop(iterator_matrix)
      ! Construct parallelisation bounds
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
                             block_matrix,found, block_shape(1), block_shape(2))
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

      CALL timestop(handle)

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
      INTEGER                                            :: backup,&
                                                            handle
      
      CALL timeset(routineN, handle)

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
      
      CALL timestop(handle)
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
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, bs_env=bs_env)

      ! build 3c integrals
      CALL compute_3c_integrals(qs_env,&
                         bs_env,&
                         tdagw_env%t_3c_w)

      ! Convert density to tensor
      ! CALL copy_fm_to_dbcsr(tdagw_env%rho_workspace(2)%matrix, tdagw_env%rho_dbcsr)
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
      ! TODO : Alternative - a Pvector approach for debugging
      ! CALL dbt_create(tdagw_env%t_one_plus_one, extra_vector, name="Pvector")
      ! CALL dbt_contract(alpha=1.0_dp,&
      !              tensor_1=V_dbt,&
      !              tensor_2=tdagw_env%t_one_plus_one,&
      !              beta=0.0_dp,&
      !              tensor_3=extra_vector,&
      !              contract_1=[2],notcontract_1=[1],map_1=[1],&
      !              contract_2=[1],notcontract_2=[2],map_2=[2],&
      !              filter_eps=bs_env%eps_filter)
      ! TODO : Get the dbt block and check results
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
      ! TODO : Alternative multiplication
      ! CALL dbt_contract(alpha=1.0_dp,&
      !              tensor_1=tdagw_env%t_3c_w,&
      !              tensor_2=extra_vector,&
      !              beta=0.0_dp,&
      !              tensor_3=tdagw_env%t_two_plus_one,&
      !              contract_1=[1], notcontract_1=[2,3], map_1=[1,2],&
      !              contract_2=[1], notcontract_2=[2], map_2=[3],&
      !              filter_eps=bs_env%eps_filter)
      ! Free memory by clearing the tensors
      CALL dbt_clear(tdagw_env%t_two_plus_one)
      CALL dbt_clear(tdagw_env%t_one_plus_one)
      CALL dbt_clear(tdagw_env%t_3c_w)
      CALL dbt_clear(tdagw_env%t_3c_work_RI__AO_AO)
      CALL dbcsr_clear(tdagw_env%rho_dbcsr)

      CALL timestop(handle)
   END SUBROUTINE
! **************************************************************************************************
!> \brief Calculates the Hartree matrix in the atomic orbital basis, given a density matrix, in local arrays
!> \param qs_env Entry point
!> \param tdagw_env Entry point of TD-aGW - uses rho_dbcsr and some complex_workspace
!> \param rho_ao Density matrix in ao basis
!> \param v_ao Overwritten by the Hartree matrix in the atomic orbital basis
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE get_hartree_local(qs_env, tdagw_env, rho_ao, v_ao)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: rho_ao
      TYPE(cp_fm_type), INTENT(OUT), POINTER             :: v_ao
      CHARACTER(len=*), PARAMETER                        :: routineN="get_hartree_local"
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iterator_matrix
      INTEGER                                            :: my_prow, my_pcol, nblkrows_local, nblkcols_local,&
                                                            i, j, k, n, nblocks, ind_1, ind_2, row_offset, col_offset,&
                                                            row_size, col_size, j_n_AO, k_n_AO, i_n_RI, n_RI,&
                                                            ri_offset, ind_i, handle
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: Pvector, Qvector
      REAL(kind=dp), DIMENSION(:,:), POINTER             :: block_matrix
      ! REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE       :: int_3c
      REAL(kind=dp), DIMENSION(:,:,:), POINTER           :: int_3c
      TYPE(dbcsr_type)                                   :: v_dbcsr_ao
      LOGICAL                                            :: found

      ! Very ineffective first implementation - calculate all 3cs an all ranks
      ! Importantly - dbcsr blocks are ordered by atoms - i.e. ethene with 6 atoms will have 6x6 block structure
      ! Number of basis states on each basis set is known is post_scf_bandstructure env
      ! TODO : Only calculate integrals which are needed for a given rank, which are given by structure
      !        of the density matrix

      CALL timeset(routineN, handle)
      ! Get the relevant environments
      CALL get_qs_env(qs_env, bs_env=bs_env, para_env=para_env)
      n_RI = tdagw_env%n_RI
      int_3c => tdagw_env%int_3c_array

      ! Allocate the Q and Pvector on each rank
      ALLOCATE(Qvector(tdagw_env%n_RI), source=0.0_dp)
      ALLOCATE(Pvector(tdagw_env%n_RI), source=0.0_dp)

      ! First step - analyze the structure of copied dbcsr matrix on both ranks
      CALL dbcsr_clear(tdagw_env%rho_dbcsr)
      CALL copy_fm_to_dbcsr(rho_ao, tdagw_env%rho_dbcsr)
      nblocks = dbcsr_get_num_blocks(tdagw_env%rho_dbcsr)
      CALL dbcsr_iterator_start(iterator_matrix, tdagw_env%rho_dbcsr)
      ! TODO : Test whether these 3c integrals do not give different results
      !         - very similar values observed
      ! CALL compute_3c_integrals(qs_env,&
      !                    bs_env,&
      !                    tdagw_env%t_3c_w)
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
            ! DEALLOCATE(int_3c)
         END DO
      END DO
      CALL dbcsr_iterator_stop(iterator_matrix)
      ! Now, each rank has contributions from D_jk within its scope
      ! Need to sum over different ranks to get the total vector on all ranks
      CALL para_env%sum(Qvector)
      ! TODO : Maybe search for matrix_vector multiply/trace scale for the V^PQ B_Q operation
      ! Once this is done, Pvector is current on all ranks
      ! Continue with V_PQ summation
      nblocks = dbcsr_get_num_blocks(tdagw_env%v_dbcsr%matrix)
      CALL dbcsr_iterator_start(iterator_matrix, tdagw_env%v_dbcsr%matrix)
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
      CALL dbcsr_create(v_dbcsr_ao, "Hartree ao", tdagw_env%rho_dbcsr)
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
            ! ALLOCATE(int_3c(bs_env%sizes_AO(ind_1), bs_env%sizes_AO(ind_2), bs_env%sizes_RI(ind_i)), source=0.0_dp)
            ! CALL build_3c_integrals(int_3c, qs_env, potential_parameter=bs_env%ri_metric,&
            !                         basis_j=bs_env%basis_set_AO, basis_k=bs_env%basis_set_AO, basis_i=bs_env%basis_set_RI,&
            !                         atom_j=ind_1, atom_k=ind_2, atom_i=ind_i)
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
                     ! block_matrix(j, k) = block_matrix(j,k) + int_3c(i, j, k) * Pvector(i+ri_offset)
                  END DO
               END DO
            END DO
            !$OMP END PARALLEL DO
            ! DEALLOCATE(int_3c)
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
END MODULE rt_tdagw

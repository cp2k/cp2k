!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the propagation via RT-BSE method.
!> \note  The control is handed directly from cp2k_runs
!> \author Stepan Marek (12.23)
! **************************************************************************************************

MODULE rt_bse
   USE cp_control_types, ONLY: dft_control_type, &
                               rtp_control_type
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE force_env_types, ONLY: force_env_type
   USE qs_mo_types, ONLY: mo_set_type, &
                          init_mo_set
   USE atomic_kind_types, ONLY: atomic_kind_type
   USE qs_kind_types, ONLY: qs_kind_type
   USE particle_types, ONLY: particle_type
   USE rt_propagation_types, ONLY: get_rtp, &
                                   rt_prop_type
   USE post_scf_bandstructure_types, ONLY: post_scf_bandstructure_type
   USE cp_fm_types, ONLY: cp_fm_type, &
                          cp_fm_p_type, &
                          cp_fm_to_fm, &
                          cp_fm_create, &
                          cp_fm_set_all, &
                          cp_fm_release, &
                          cp_fm_get_element, &
                          cp_fm_set_element, &
                          cp_fm_get_info, &
                          cp_fm_write_unformatted, &
                          cp_fm_read_unformatted, &
                          cp_fm_write_formatted
   USE cp_cfm_types, ONLY: cp_cfm_type, &
                           cp_cfm_p_type, &
                           cp_fm_to_cfm, &
                           cp_cfm_to_cfm, &
                           cp_cfm_to_fm, &
                           cp_cfm_create, &
                           cp_cfm_get_info, &
                           cp_cfm_release
   USE kinds, ONLY: dp, &
                    default_path_length
   USE cp_dbcsr_api, ONLY: dbcsr_p_type, &
                           dbcsr_type, &
                           dbcsr_print, &
                           dbcsr_has_symmetry, &
                           dbcsr_desymmetrize, &
                           dbcsr_create, &
                           dbcsr_release, &
                           dbcsr_copy, &
                           dbcsr_scale, &
                           dbcsr_add, &
                           dbcsr_set, &
                           dbcsr_clear, &
                           dbcsr_setname, &
                           dbcsr_iterator_type, &
                           dbcsr_iterator_start, &
                           dbcsr_iterator_stop, &
                           dbcsr_iterator_blocks_left, &
                           dbcsr_iterator_next_block, &
                           dbcsr_put_block, &
                           dbcsr_reserve_blocks, &
                           dbcsr_get_num_blocks, &
                           dbcsr_get_block_p, &
                           dbcsr_get_info
   USE OMP_LIB, ONLY: omp_get_thread_num, &
                      omp_get_num_threads, &
                      omp_set_num_threads, &
                      omp_get_max_threads
   USE dbt_api, ONLY: dbt_clear, &
                      dbt_contract, &
                      dbt_copy_matrix_to_tensor, &
                      dbt_copy_tensor_to_matrix, &
                      dbt_type
   USE libint_2c_3c, ONLY: libint_potential_type
   USE mp2_ri_2c, ONLY: RI_2c_integral_mat
   USE qs_tensors, ONLY: neighbor_list_3c_destroy, &
                         build_2c_integrals, &
                         build_2c_neighbor_lists
   USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type, &
                                     release_neighbor_list_sets
   USE cp_dbcsr_operations, ONLY: copy_dbcsr_to_fm, &
                                  copy_fm_to_dbcsr, &
                                  dbcsr_allocate_matrix_set, &
                                  dbcsr_deallocate_matrix_set, &
                                  cp_dbcsr_sm_fm_multiply, &
                                  copy_cfm_to_dbcsr, &
                                  copy_dbcsr_to_cfm
   USE cp_fm_basic_linalg, ONLY: cp_fm_scale, &
                                 cp_fm_invert, &
                                 cp_fm_trace, &
                                 cp_fm_transpose, &
                                 cp_fm_norm, &
                                 cp_fm_column_scale, &
                                 cp_fm_scale_and_add
   USE cp_cfm_basic_linalg, ONLY: cp_cfm_scale_and_add, &
                                  cp_cfm_scale, &
                                  cp_cfm_transpose, &
                                  cp_cfm_norm, &
                                  cp_cfm_trace, &
                                  cp_cfm_column_scale
   USE cp_fm_diag, ONLY: cp_fm_geeig
   USE cp_cfm_diag, ONLY: cp_cfm_geeig
   USE parallel_gemm_api, ONLY: parallel_gemm
   USE qs_moments, ONLY: build_local_moment_matrix, &
                         build_berry_moment_matrix
   USE moments_utils, ONLY: get_reference_point
   USE qs_mo_io, ONLY: read_mo_set_from_restart
   USE qs_ks_methods, ONLY: qs_ks_build_kohn_sham_matrix
   USE force_env_methods, ONLY: force_env_calc_energy_force
   USE qs_energy_init, ONLY: qs_energies_init
   USE qs_energy_utils, ONLY: qs_energies_properties
   USE efield_utils, ONLY: make_field
   USE rt_propagator_init, ONLY: rt_initialize_rho_from_mos
   USE rt_propagation_methods, ONLY: s_matrices_create
   USE gw_large_cell_Gamma, ONLY: compute_3c_integrals
   USE gw_integrals, ONLY: build_3c_integral_block
   USE matrix_exp, ONLY: taylor_full_complex
   USE message_passing, ONLY: mp_para_env_type
   USE input_section_types, ONLY: section_vals_get_subs_vals, &
                                  section_vals_type
   USE cell_types, ONLY: cell_type
   USE input_constants, ONLY: rtp_bse_ham_ks, &
                              rtp_bse_ham_g0w0, &
                              use_mom_ref_coac, &
                              use_mom_ref_zero, &
                              do_taylor, &
                              do_bch, &
                              do_exact, &
                              use_rt_restart
   USE rt_bse_types, ONLY: rtbse_env_type, &
                           create_rtbse_env, &
                           release_rtbse_env, &
                           multiply_cfm_fm, &
                           multiply_fm_cfm
   USE rt_bse_io, ONLY: output_moments, &
                        read_moments, &
                        read_field, &
                        output_moments_ft, &
                        output_polarizability, &
                        output_field, &
                        output_mos_contravariant, &
                        read_restart, &
                        output_restart, &
                        print_timestep_info, &
                        print_etrs_info_header, &
                        print_etrs_info, &
                        print_rtbse_header_info
   USE qs_ks_types, ONLY: qs_ks_env_type
   USE qs_integrate_potential_product, ONLY: integrate_v_rspace
   USE qs_collocate_density, ONLY: calculate_rho_elec
   USE cp_files, ONLY: open_file, &
                       file_exists, &
                       close_file
   USE physcon, ONLY: femtoseconds
   USE mathconstants, ONLY: twopi

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_bse"

   #:include "rt_bse_macros.fypp"

   PUBLIC :: run_propagation_bse, &
             get_hartree, &
             get_sigma

   INTERFACE get_sigma
      MODULE PROCEDURE get_sigma_complex, &
         get_sigma_real, &
         get_sigma_dbcsr, &
         get_sigma_noenv
   END INTERFACE
   INTERFACE get_hartree
      MODULE PROCEDURE get_hartree_env, &
         get_hartree_noenv
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Runs the electron-only real time BSE propagation
!> \param qs_env Quickstep environment data, containing the config from input files
!> \param force_env Force environment data, entry point of the calculation
! **************************************************************************************************
   SUBROUTINE run_propagation_bse(qs_env, force_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(force_env_type), POINTER                      :: force_env
      CHARACTER(len=*), PARAMETER                        :: routineN = 'run_propagation_bse'
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      INTEGER                                            :: i, j, k, handle
      LOGICAL                                            :: converged
      REAL(kind=dp)                                      :: metric, enum_re, enum_im, &
                                                            idempotence_dev, a_metric_1, a_metric_2

      CALL timeset(routineN, handle)

      ! Run the initial SCF calculation / read SCF restart information
      CALL force_env_calc_energy_force(force_env, calc_force=.FALSE., consistent_energies=.FALSE.)

      ! Allocate all persistant storage and read input that does not need further processing
      CALL create_rtbse_env(rtbse_env, qs_env, force_env)

      IF (rtbse_env%unit_nr > 0) CALL print_rtbse_header_info(rtbse_env)

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
      CALL initialize_rtbse_env(rtbse_env)

      ! Setup the time based on the starting step
      ! Assumes identical dt between two runs
      rtbse_env%sim_time = REAL(rtbse_env%sim_start, dp)*rtbse_env%sim_dt
      ! Output 0 time moments and field
      IF (.NOT. rtbse_env%restart_extracted) THEN
         CALL output_field(rtbse_env)
         CALL output_moments(rtbse_env, rtbse_env%rho)
      END IF

      ! Do not apply the delta kick if we are doing a restart calculation
      IF (rtbse_env%dft_control%rtp_control%apply_delta_pulse .AND. (.NOT. rtbse_env%restart_extracted)) THEN
         CALL apply_delta_pulse(rtbse_env)
      END IF

      ! ********************** Start the time loop **********************
      DO i = rtbse_env%sim_start, rtbse_env%sim_nsteps

         ! Update the simulation time
         rtbse_env%sim_time = REAL(i, dp)*rtbse_env%sim_dt
         rtbse_env%sim_step = i
         ! Carry out the ETRS self-consistent propagation - propagates rho to rho_new (through rho_M)
         CALL etrs_scf_loop(rtbse_env, rtbse_env%rho, rtbse_env%rho_M, rtbse_env%rho_new, converged, k, metric)
         CALL get_electron_number(rtbse_env, rtbse_env%rho_new, enum_re, enum_im)
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

      ! Carry out the FT
      CALL output_moments_ft(rtbse_env)
      CALL output_polarizability(rtbse_env)

      ! Deallocate everything
      CALL release_rtbse_env(rtbse_env)

      CALL timestop(handle)
   END SUBROUTINE run_propagation_bse

! **************************************************************************************************
!> \brief Calculates the initial values, based on restart/scf density, and other non-trivial values
!> \param rtbse_env RT-BSE environment
!> \param qs_env Quickstep environment (needed for reference to previous calculations)
!> \author Stepan Marek (09.24)
! **************************************************************************************************
   SUBROUTINE initialize_rtbse_env(rtbse_env)
      TYPE(rtbse_env_type), POINTER                :: rtbse_env
      CHARACTER(len=*), PARAMETER                  :: routineN = "initialize_rtbse_env"
      TYPE(post_scf_bandstructure_type), POINTER   :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER    :: moments_dbcsr_p
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER    :: matrix_s
      REAL(kind=dp), DIMENSION(:), POINTER         :: custom_ref_point, &
                                                      occupations
      REAL(kind=dp), DIMENSION(3)                  :: rpoint
      INTEGER                                      :: i, k, handle

      CALL timeset(routineN, handle)

      ! Get pointers to parameters from qs_env
      CALL get_qs_env(rtbse_env%qs_env, bs_env=bs_env, matrix_s=matrix_s)

      ! ****** START MOMENTS OPERATOR CALCULATION
      ! Construct moments from dbcsr
      NULLIFY (moments_dbcsr_p)
      ALLOCATE (moments_dbcsr_p(3))
      DO k = 1, 3
         ! Make sure the pointer is empty
         NULLIFY (moments_dbcsr_p(k)%matrix)
         ! Allocate a new matrix that the pointer points to
         ALLOCATE (moments_dbcsr_p(k)%matrix)
         ! Create the matrix storage - matrix copies the structure of overlap matrix
         CALL dbcsr_copy(moments_dbcsr_p(k)%matrix, matrix_s(1)%matrix)
      END DO
      ! Run the moment calculation
      ! check for presence to prevent memory errors
      NULLIFY (custom_ref_point)
      ALLOCATE (custom_ref_point(3), source=0.0_dp)
      rpoint(:) = 0.0_dp
      CALL get_reference_point(rpoint, qs_env=rtbse_env%qs_env, reference=use_mom_ref_coac, ref_point=custom_ref_point)
      CALL build_local_moment_matrix(rtbse_env%qs_env, moments_dbcsr_p, 1, rpoint)
      ! Copy to full matrix
      DO k = 1, 3
         ! Again, matrices are created from overlap template
         CALL copy_dbcsr_to_fm(moments_dbcsr_p(k)%matrix, rtbse_env%moments(k))
      END DO
      ! Now, repeat without reference point to get the moments for field
      CALL get_reference_point(rpoint, qs_env=rtbse_env%qs_env, reference=use_mom_ref_zero, ref_point=custom_ref_point)
      DEALLOCATE (custom_ref_point)
      CALL build_local_moment_matrix(rtbse_env%qs_env, moments_dbcsr_p, 1, rpoint)
      DO k = 1, 3
         CALL copy_dbcsr_to_fm(moments_dbcsr_p(k)%matrix, rtbse_env%moments_field(k))
      END DO

      ! Now can deallocate dbcsr matrices
      DO k = 1, 3
         CALL dbcsr_release(moments_dbcsr_p(k)%matrix)
         DEALLOCATE (moments_dbcsr_p(k)%matrix)
      END DO
      DEALLOCATE (moments_dbcsr_p)
      ! ****** END MOMENTS OPERATOR CALCULATION

      ! ****** START INITIAL DENSITY MATRIX CALCULATION
      ! Get the rho from fm_MOS
      ! Uses real orbitals only - no kpoints
      ALLOCATE (occupations(rtbse_env%n_ao))
      ! Iterate over both spins
      DO i = 1, rtbse_env%n_spin
         occupations(:) = 0.0_dp
         occupations(1:rtbse_env%n_occ(i)) = 1.0_dp
         ! Create real part
         CALL cp_fm_to_fm(bs_env%fm_mo_coeff_Gamma(i), rtbse_env%real_workspace(1))
         CALL cp_fm_column_scale(rtbse_env%real_workspace(1), occupations)
         CALL parallel_gemm("N", "T", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                            1.0_dp, rtbse_env%real_workspace(1), bs_env%fm_mo_coeff_Gamma(i), &
                            0.0_dp, rtbse_env%real_workspace(2))
         ! Sets imaginary part to zero
         CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace(2), mtarget=rtbse_env%rho(i))
         ! Save the reference value for the case of delta kick
         CALL cp_cfm_to_cfm(rtbse_env%rho(i), rtbse_env%rho_orig(i))
      END DO
      DEALLOCATE (occupations)
      ! If the restart field is provided, overwrite rho from restart
      IF (rtbse_env%dft_control%rtp_control%initial_wfn == use_rt_restart) THEN
         CALL read_restart(rtbse_env)
      END IF
      ! ****** END INITIAL DENSITY MATRIX CALCULATION
      ! ****** START MOMENTS TRACE LOAD
      ! The moments are only loaded if the consistent save files exist
      CALL read_moments(rtbse_env)
      ! ****** END MOMENTS TRACE LOAD

      ! ****** START FIELD TRACE LOAD
      ! The moments are only loaded if the consistent save files exist
      CALL read_field(rtbse_env)
      ! ****** END FIELD TRACE LOAD

      ! ****** START OVERLAP + INVERSE OVERLAP CALCULATION
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, rtbse_env%S_fm)
      CALL cp_fm_to_cfm(msourcer=rtbse_env%S_fm, mtarget=rtbse_env%S_cfm)
      CALL cp_fm_invert(rtbse_env%S_fm, rtbse_env%S_inv_fm)
      ! ****** END OVERLAP + INVERSE OVERLAP CALCULATION

      ! ****** START SINGLE PARTICLE HAMILTONIAN CALCULATION
      DO i = 1, rtbse_env%n_spin
         IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0) THEN
            ! G0W0 Hamiltonian
            CALL cp_fm_to_fm(bs_env%fm_mo_coeff_Gamma(i), rtbse_env%real_workspace(1))
            ! NOTE : Gamma point is not always the zero k-point
            ! C * Lambda
            CALL cp_fm_column_scale(rtbse_env%real_workspace(1), bs_env%eigenval_G0W0(:, 1, i))
            ! C * Lambda * C^T
            CALL parallel_gemm("N", "T", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                               1.0_dp, rtbse_env%real_workspace(1), bs_env%fm_mo_coeff_Gamma(i), &
                               0.0_dp, rtbse_env%real_workspace(2))
            ! S * C * Lambda * C^T
            CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                               1.0_dp, rtbse_env%S_fm, rtbse_env%real_workspace(2), &
                               0.0_dp, rtbse_env%real_workspace(1))
            ! S * C * Lambda * C^T * S = H
            CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                               1.0_dp, rtbse_env%real_workspace(1), rtbse_env%S_fm, &
                               0.0_dp, rtbse_env%real_workspace(2))
            CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace(2), mtarget=rtbse_env%ham_reference(i))
         ELSE
            ! KS Hamiltonian
            CALL cp_fm_to_cfm(msourcer=bs_env%fm_ks_Gamma(i), mtarget=rtbse_env%ham_reference(i))
         END IF
      END DO
      ! ****** END SINGLE PARTICLE HAMILTONIAN CALCULATION

      ! ****** START HARTREE POTENTIAL REFERENCE CALCULATION
      ! Calculate Coulomb RI elements, necessary for Hartree calculation
      CALL init_hartree(rtbse_env, rtbse_env%v_dbcsr)
      IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0) THEN
         ! In a non-HF calculation, copy the actual correlation part of the interaction
         CALL copy_fm_to_dbcsr(bs_env%fm_W_MIC_freq_zero, rtbse_env%w_dbcsr)
      ELSE
         ! In HF, correlation is set to zero
         CALL dbcsr_set(rtbse_env%w_dbcsr, 0.0_dp)
      END IF
      ! Add the Hartree to the screened_dbt tensor - now W = V + W^c
      CALL dbcsr_add(rtbse_env%w_dbcsr, rtbse_env%v_dbcsr, 1.0_dp, 1.0_dp)
      CALL dbt_copy_matrix_to_tensor(rtbse_env%w_dbcsr, rtbse_env%screened_dbt)
      ! Calculate the original Hartree potential
      ! Uses rho_orig - same as rho for initial run but different for continued run
      DO i = 1, rtbse_env%n_spin
         CALL get_hartree(rtbse_env, rtbse_env%rho_orig(i), rtbse_env%hartree_curr(i))
         ! Scaling by spin degeneracy
         CALL cp_fm_scale(rtbse_env%spin_degeneracy, rtbse_env%hartree_curr(i))
         ! Subtract the reference from the reference Hamiltonian
         CALL cp_fm_to_cfm(msourcer=rtbse_env%hartree_curr(i), mtarget=rtbse_env%ham_workspace(1))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                   CMPLX(-1.0, 0.0, kind=dp), rtbse_env%ham_workspace(1))
      END DO
      ! ****** END HARTREE POTENTIAL REFERENCE CALCULATION

      ! ****** START COHSEX REFERENCE CALCULATION
      ! Calculate the COHSEX starting energies
      IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0) THEN
         ! Subtract the v_xc from COH part of the self-energy, as V_xc is also not updated during the timestepping
         DO i = 1, rtbse_env%n_spin
            ! TODO : Allow no COH calculation for static screening
            CALL get_sigma(rtbse_env, rtbse_env%sigma_COH(i), -0.5_dp, rtbse_env%S_inv_fm)
            ! Copy and subtract from the complex reference hamiltonian
            CALL cp_fm_to_cfm(msourcer=rtbse_env%sigma_COH(i), mtarget=rtbse_env%ham_workspace(1))
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                      CMPLX(-1.0, 0.0, kind=dp), rtbse_env%ham_workspace(1))
            ! Calculate exchange part - TODO : should this be applied for different spins? - TEST with O2 HF propagation?
            ! So far only closed shell tested
            ! Uses rho_orig - same as rho for initial run but different for continued run
            CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX(i), -1.0_dp, rtbse_env%rho_orig(i))
            ! Subtract from the complex reference Hamiltonian
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                      CMPLX(-1.0, 0.0, kind=dp), rtbse_env%sigma_SEX(i))
         END DO
      ELSE
         ! KS Hamiltonian - use time-dependent Fock exchange
         DO i = 1, rtbse_env%n_spin
            CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX(i), -1.0_dp, rtbse_env%rho_orig(i))
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_reference(i), &
                                      CMPLX(-1.0, 0.0, kind=dp), rtbse_env%sigma_SEX(i))
         END DO
      END IF
      ! ****** END COHSEX REFERENCE CALCULATION

      CALL timestop(handle)
   END SUBROUTINE initialize_rtbse_env

! **************************************************************************************************
!> \brief Custom reimplementation of the delta pulse routines
!> \param rtbse_env RT-BSE environment
!> \author Stepan Marek (09.24)
! **************************************************************************************************
   SUBROUTINE apply_delta_pulse(rtbse_env)
      TYPE(rtbse_env_type), POINTER                :: rtbse_env
      CHARACTER(len=*), PARAMETER                  :: routineN = "apply_delta_pulse"
      REAL(kind=dp)                                :: intensity, metric
      REAL(kind=dp), DIMENSION(3)                  :: kvec
      INTEGER                                      :: i, k, handle

      CALL timeset(routineN, handle)

      ! Report application
      IF (rtbse_env%unit_nr > 0) WRITE (rtbse_env%unit_nr, '(A28)') ' RTBSE| Applying delta pulse'
      ! Extra minus for the propagation of density
      intensity = -rtbse_env%dft_control%rtp_control%delta_pulse_scale
      metric = 0.0_dp
      kvec(:) = rtbse_env%dft_control%rtp_control%delta_pulse_direction(:)
      IF (rtbse_env%unit_nr > 0) WRITE (rtbse_env%unit_nr, '(A38,E14.4E3,E14.4E3,E14.4E3)') &
         " RTBSE| Delta pulse elements (a.u.) : ", intensity*kvec(:)
      ! So far no spin dependence, but can be added by different structure of delta pulse
      CALL cp_fm_set_all(rtbse_env%real_workspace(1), 0.0_dp)
      DO k = 1, 3
         CALL cp_fm_scale_and_add(1.0_dp, rtbse_env%real_workspace(1), &
                                  kvec(k), rtbse_env%moments_field(k))
      END DO
      ! enforce hermiticity of the effective Hamiltonian
      CALL cp_fm_transpose(rtbse_env%real_workspace(1), rtbse_env%real_workspace(2))
      CALL cp_fm_scale_and_add(0.5_dp, rtbse_env%real_workspace(1), &
                               0.5_dp, rtbse_env%real_workspace(2))
      ! Prepare the exponential/exponent for propagation
      IF (rtbse_env%mat_exp_method == do_bch) THEN
         ! Multiply by the S_inv matrix - in the classic ordering
         CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                            intensity, rtbse_env%S_inv_fm, rtbse_env%real_workspace(1), &
                            0.0_dp, rtbse_env%real_workspace(2))
         DO i = 1, rtbse_env%n_spin
            ! Sets real part to zero
            CALL cp_fm_to_cfm(msourcei=rtbse_env%real_workspace(2), mtarget=rtbse_env%ham_workspace(i))
         END DO
      ELSE IF (rtbse_env%mat_exp_method == do_exact) THEN
         DO i = 1, rtbse_env%n_spin
            CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace(1), mtarget=rtbse_env%ham_effective(i))
            CALL cp_cfm_gexp(rtbse_env%ham_effective(i), rtbse_env%S_cfm, rtbse_env%ham_workspace(i), &
                             CMPLX(0.0, intensity, kind=dp), rtbse_env%rho_workspace)
         END DO
      END IF
      ! Propagate the density by the effect of the delta pulse
      CALL propagate_density(rtbse_env, rtbse_env%ham_workspace, rtbse_env%rho, rtbse_env%rho_new)
      metric = rho_metric(rtbse_env%rho_new, rtbse_env%rho, rtbse_env%n_spin)
      IF (rtbse_env%unit_nr > 0) WRITE (rtbse_env%unit_nr, ('(A42,E38.8E3)')) " RTBSE| Metric difference after delta kick", metric
      ! Copy the new density to the old density
      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_to_cfm(rtbse_env%rho_new(i), rtbse_env%rho(i))
      END DO

      CALL timestop(handle)
   END SUBROUTINE apply_delta_pulse
! **************************************************************************************************
!> \brief Determines the metric for the density matrix, used for convergence criterion
!> \param rho_new Array of new density matrices (one for each spin index)
!> \param rho_old Array of old density matrices (one for each spin index)
!> \param nspin Number of spin indices
!> \param workspace_opt Optionally provide external workspace to save some allocation time
! **************************************************************************************************
   FUNCTION rho_metric(rho_new, rho_old, nspin, workspace_opt) RESULT(metric)
      TYPE(cp_cfm_type), DIMENSION(:), POINTER, INTENT(IN):: rho_new, &
                                                             rho_old
      INTEGER, INTENT(IN)                                :: nspin
      TYPE(cp_cfm_type), POINTER, OPTIONAL               :: workspace_opt
      TYPE(cp_cfm_type)                                  :: workspace
      REAL(kind=dp)                                      :: metric
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: partial_metric
      INTEGER                                            :: j
      COMPLEX(kind=dp)                                   :: scale_factor

      ALLOCATE (partial_metric(nspin))

      ! Only allocate/deallocate storage if required
      IF (PRESENT(workspace_opt)) THEN
         workspace = workspace_opt
      ELSE
         CALL cp_cfm_create(workspace, rho_new(1)%matrix_struct)
      END IF
      scale_factor = 1.0
      DO j = 1, nspin
         CALL cp_cfm_to_cfm(rho_new(j), workspace)
         ! Get the difference in the resulting matrix
         CALL cp_cfm_scale_and_add(scale_factor, workspace, -scale_factor, rho_old(j))
         ! Now, get the relevant number
         partial_metric(j) = cp_cfm_norm(workspace, 'M')
      END DO
      metric = 0.0_dp
      ! For more than one spin, do Cartesian sum of the different spin norms
      DO j = 1, nspin
         metric = metric + partial_metric(j)*partial_metric(j)
      END DO
      metric = SQRT(metric)
      ! Deallocate workspace
      IF (.NOT. PRESENT(workspace_opt)) CALL cp_cfm_release(workspace)
      DEALLOCATE (partial_metric)
   END FUNCTION

! **************************************************************************************************
!> \brief Determines the metric of the antihermitian part of the matrix
!> \param real_fm Real part of the full matrix
!> \param imag_fm Imaginary part of the full matrix
! **************************************************************************************************
   SUBROUTINE antiherm_metric(real_fm, imag_fm, workspace, metric)
      TYPE(cp_fm_type), INTENT(IN)                      :: real_fm
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL            :: imag_fm
      REAL(kind=dp), INTENT(OUT)                        :: metric
      TYPE(cp_cfm_type), DIMENSION(:), POINTER          :: workspace
      COMPLEX(kind=dp)                                  :: complex_one

      ! Get the complex and complex conjugate matrix
      IF (PRESENT(imag_fm)) THEN
         CALL cp_fm_to_cfm(real_fm, imag_fm, workspace(1))
      ELSE
         CALL cp_fm_to_cfm(msourcer=real_fm, mtarget=workspace(1))
      END IF
      CALL cp_cfm_transpose(workspace(1), "C", workspace(2))
      ! Subtract these, and get the metric
      complex_one = CMPLX(1.0, 0.0, kind=dp)
      CALL cp_cfm_scale_and_add(complex_one, workspace(1), -complex_one, workspace(2))
      metric = cp_cfm_norm(workspace(1), "M")
   END SUBROUTINE

! **************************************************************************************************
!> \brief For Taylor and Exact exp_method, calculates the matrix exponential of the
!>        effective Hamiltonian. For BCH, calculates just the effective Hamiltonian. For other methods,
!>        aborts the execution, as they are not implemented yet.
!> \param rtbse_env Entry point of the calculation. Uses rho_workspace for Taylor and BCH. For exact,
!>                  uses complex_workspace, complex_ham, complex_s, real_eigvals and exp_eigvals.
!>                  Results are stored in ham_workspace.
! **************************************************************************************************
   SUBROUTINE ham_to_exp(rtbse_env, ham, ham_exp)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: ham, &
                                                            ham_exp
      CHARACTER(len=*), PARAMETER                        :: routineN = "ham_to_exp"
      INTEGER                                            :: j, handle
      CALL timeset(routineN, handle)
      DO j = 1, rtbse_env%n_spin
         IF (rtbse_env%mat_exp_method == do_bch) THEN
            ! In Taylor and BCH, we first evaluate the entire exponent and then evaluate exponential in series
            ! In order to produce correct result, need to remultiply by inverse overlap matrix
            CALL multiply_fm_cfm("N", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                                 1.0_dp, rtbse_env%S_inv_fm, ham(j), &
                                 0.0_dp, rtbse_env%rho_workspace(1))

            ! The evolution of density matrix is derived from the right multiplying term
            ! Imaginary part of the exponent = -real part of the matrix
            CALL cp_cfm_scale(CMPLX(0.0, -rtbse_env%sim_dt/2, kind=dp), rtbse_env%rho_workspace(1))
            ! In BCH, exponential is not calculated explicitly, but the propagation is solved in series
            CALL cp_cfm_to_cfm(rtbse_env%rho_workspace(1), ham_exp(j))
         ELSE IF (rtbse_env%mat_exp_method == do_exact) THEN
            CALL cp_cfm_gexp(ham(j), rtbse_env%S_cfm, ham_exp(j), &
                             CMPLX(0.0, -rtbse_env%sim_dt/2, kind=dp), rtbse_env%rho_workspace)
         ELSE
            CPABORT("Only BCH and Taylor matrix exponentiation implemented")
         END IF
      END DO

      CALL timestop(handle)
   END SUBROUTINE
! **************************************************************************************************
!> \brief Updates the effective Hamiltonian, given a density matrix rho
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param qs_env QS env
!> \param rho Real and imaginary parts ( + spin) of the density at current time
! **************************************************************************************************
   SUBROUTINE update_effective_ham(rtbse_env, rho)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho
      CHARACTER(len=*), PARAMETER                        :: routineN = "update_effective_ham"
      INTEGER                                            :: k, j, nspin, handle

      CALL timeset(routineN, handle)
      ! Shorthand
      nspin = rtbse_env%n_spin
      ! Reset the effective Hamiltonian to KS Hamiltonian + G0W0 - reference COHSEX - reference Hartree
      DO j = 1, nspin
         ! Sets the imaginary part to zero
         CALL cp_cfm_to_cfm(rtbse_env%ham_reference(j), rtbse_env%ham_effective(j))
      END DO
      ! Determine the field at current time
      IF (rtbse_env%dft_control%apply_efield_field) THEN
         CALL make_field(rtbse_env%dft_control, rtbse_env%field, rtbse_env%sim_step, rtbse_env%sim_time)
      ELSE
         ! No field
         rtbse_env%field(:) = 0.0_dp
      END IF
      DO j = 1, nspin
         DO k = 1, 3
            ! Minus sign due to charge of electrons
            CALL cp_fm_to_cfm(msourcer=rtbse_env%moments_field(k), mtarget=rtbse_env%ham_workspace(1))
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_effective(j), &
                                      CMPLX(rtbse_env%field(k), 0.0, kind=dp), rtbse_env%ham_workspace(1))
         END DO
         IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0) THEN
            ! Add the COH part - so far static but can be dynamic in principle throught the W updates
            CALL get_sigma(rtbse_env, rtbse_env%sigma_COH(j), -0.5_dp, rtbse_env%S_inv_fm)
            CALL cp_fm_to_cfm(msourcer=rtbse_env%sigma_COH(j), mtarget=rtbse_env%ham_workspace(1))
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_effective(j), &
                                      CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_workspace(1))
         END IF
         ! Calculate the (S)EX part - based on provided rho
         ! iGW = - rho W
         CALL get_sigma(rtbse_env, rtbse_env%sigma_SEX(j), -1.0_dp, rho(j))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_effective(j), &
                                   CMPLX(1.0, 0.0, kind=dp), rtbse_env%sigma_SEX(j))
         ! Calculate Hartree potential
         ! Hartree potential is scaled by number of electrons in each MO - spin degeneracy
         CALL get_hartree(rtbse_env, rho(j), &
                          rtbse_env%hartree_curr(j))
         CALL cp_fm_to_cfm(msourcer=rtbse_env%hartree_curr(j), mtarget=rtbse_env%ham_workspace(1))
         CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_effective(j), &
                                   CMPLX(rtbse_env%spin_degeneracy, 0.0, kind=dp), rtbse_env%ham_workspace(1))
         ! Enforce hermiticity of the effective Hamiltonian
         ! Important components without forced Hermiticity - moments matrix, sigma matrices, Hartree matrix
         ! single particle Ham
         CALL cp_cfm_transpose(rtbse_env%ham_effective(j), 'C', rtbse_env%ham_workspace(1))
         CALL cp_cfm_scale_and_add(CMPLX(0.5, 0.0, kind=dp), rtbse_env%ham_effective(j), &
                                   CMPLX(0.5, 0.0, kind=dp), rtbse_env%ham_workspace(1))
      END DO
      CALL timestop(handle)
   END SUBROUTINE update_effective_ham
! **************************************************************************************************
!> \brief Self-consistently (ETRS) propagates the density to the next timestep
!> \note Uses rtbse_env%rho_new_last, assumes correct timestep information is given in rtbse_env
!> \param rho_start Initial density matrix
!> \param rho_mid Midpoint density (propagated to by the initial Hamiltonian)
!> \param rho_end Endpoint density (propagated to by endpoint Hamiltonian)
!> \param converged Whether the resulting rho_end is self-consistent
!> \param k How many SC iterations were done
!> \param metric The difference metric from the last self-consistent iteration (for printing/evaluation)
! **************************************************************************************************
   SUBROUTINE etrs_scf_loop(rtbse_env, rho_start, rho_mid, rho_end, converged, k, metric)
      TYPE(rtbse_env_type), POINTER                     :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER          :: rho_start, &
                                                           rho_mid, &
                                                           rho_end
      LOGICAL                                           :: converged
      INTEGER                                           :: k
      REAL(kind=dp)                                     :: metric
      CHARACTER(len=*), PARAMETER                       :: routineN = "etrs_scf_loop"
      INTEGER                                           :: j, handle

      CALL timeset(routineN, handle)

      ! This method determines the density matrix at time (t+dt) by guessing the effective Hamiltonian at (t + dt)
      ! and using the Hamiltonian at time (t), it propagates density from time (t) while ensuring that the density
      ! at (t + dt/2) is the same for both forward and backwards propagation. Then, density at (t + dt) is
      ! used to calculate the new Hamiltonian at (t+dt), which is then used to get the new propagator, and so on
      ! until the density matrix does not change within certain limit
      ! Pseudocode of the algorithm
      !      rho_M = exp(-i S^(-1) H[rho(t)] dt/2) rho(t) exp(i H[rho(t)] S^(-1) dt/2)
      !      rho(t+dt, 0) = rho_M
      !      for j in 1,max_self_iter
      !              rho(t+dt,j) = exp(- i S^(-1) H[rho(t+dt,j-1)] dt/2) rho_M exp(i H [rho(t+dt,j-1)] S^(-1) dt/2)
      !              if ||rho(t+dt,j) - rho(t+dt,j-1)|| < epsilon
      !                      break

      ! Initial setup - calculate the Hamiltonian
      CALL update_effective_ham(rtbse_env, rho_start)
      ! Create the exponential
      CALL ham_to_exp(rtbse_env, rtbse_env%ham_effective, rtbse_env%ham_workspace)
      ! Propagate to rho_mid
      CALL propagate_density(rtbse_env, rtbse_env%ham_workspace, rho_start, rho_mid)
      ! Propagate to initial guess
      CALL propagate_density(rtbse_env, rtbse_env%ham_workspace, rho_mid, rtbse_env%rho_new_last)
      ! Update bookkeeping to the next timestep - Hamiltonians are now evaluated at the next timestep
      rtbse_env%sim_step = rtbse_env%sim_step + 1
      rtbse_env%sim_time = rtbse_env%sim_time + rtbse_env%sim_dt
      converged = .FALSE.
      CALL print_etrs_info_header(rtbse_env)
      DO k = 1, rtbse_env%etrs_max_iter
         ! Get the Hamiltonian following from the last timestep
         CALL update_effective_ham(rtbse_env, rtbse_env%rho_new_last)
         CALL ham_to_exp(rtbse_env, rtbse_env%ham_effective, rtbse_env%ham_workspace)
         ! Propagate to new guess
         CALL propagate_density(rtbse_env, rtbse_env%ham_workspace, rho_mid, rho_end)
         ! Check for self-consistency
         metric = rho_metric(rho_end, rtbse_env%rho_new_last, rtbse_env%n_spin)
         ! ETRS info - only for log level > medium
         CALL print_etrs_info(rtbse_env, k, metric)
         IF (metric < rtbse_env%etrs_threshold) THEN
            converged = .TRUE.
            EXIT
         ELSE
            ! Copy rho_new to rho_new_last
            DO j = 1, rtbse_env%n_spin
               ! Leaving for free convergence
               CALL cp_cfm_to_cfm(rho_end(j), rtbse_env%rho_new_last(j))
            END DO
         END IF
      END DO
      ! Error handling in the case where the propagation did not converge is left to the main routine
      CALL timestop(handle)
   END SUBROUTINE etrs_scf_loop

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
      TYPE(cp_cfm_type), DIMENSION(:), POINTER          :: propagator_matrix
      ! Matrix Y to be propagated into matrix Y'
      TYPE(cp_cfm_type), DIMENSION(:), POINTER          :: target_matrix
      ! Matrix Y' is stored here on exit
      TYPE(cp_cfm_type), DIMENSION(:), POINTER          :: result_matrix, workspace
      ! Threshold for the metric which decides when to truncate the BCH expansion
      REAL(kind=dp), OPTIONAL                           :: threshold_opt
      INTEGER, OPTIONAL                                 :: max_iter_opt
      CHARACTER(len=*), PARAMETER                       :: routineN = "bch_propagate"
      REAL(kind=dp)                                     :: threshold, prefactor, metric
      INTEGER                                           :: max_iter, i, n_spin, n_ao, k, &
                                                           w_stride, handle
      LOGICAL                                           :: converged
      CHARACTER(len=77)                                 :: error

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
      CALL cp_cfm_get_info(target_matrix(1), nrow_global=n_ao)
      w_stride = n_spin

      ! Initiate
      DO i = 1, n_spin
         CALL cp_cfm_to_cfm(target_matrix(i), result_matrix(i))
         CALL cp_cfm_to_cfm(target_matrix(i), workspace(i))
      END DO

      ! Start the BCH iterations
      ! So far, no spin mixing terms
      DO k = 1, max_iter
         prefactor = 1.0_dp/REAL(k, kind=dp)
         DO i = 1, n_spin
            CALL parallel_gemm("N", "N", n_ao, n_ao, n_ao, &
                               CMPLX(prefactor, 0.0, kind=dp), propagator_matrix(i), workspace(i), &
                               CMPLX(0.0, 0.0, kind=dp), workspace(i + w_stride))
            CALL parallel_gemm("N", "C", n_ao, n_ao, n_ao, &
                               CMPLX(prefactor, 0.0, kind=dp), workspace(i), propagator_matrix(i), &
                               CMPLX(1.0, 0.0, kind=dp), workspace(i + w_stride))
            ! Add to the result
            CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), result_matrix(i), &
                                      CMPLX(1.0, 0.0, kind=dp), workspace(i + w_stride))
         END DO
         metric = rho_metric(workspace(w_stride + 1:), workspace(1:w_stride), n_spin)
         IF (metric <= threshold) THEN
            converged = .TRUE.
            EXIT
         ELSE
            DO i = 1, n_spin
               CALL cp_cfm_to_cfm(workspace(i + w_stride), workspace(i))
            END DO
         END IF
      END DO
      IF (.NOT. converged) THEN
         WRITE (error, '(A35,E13.4E3,A16,E13.4E3)') "BCH did not converge, BCH Metric : ", &
            metric, "BCH Threshold : ", threshold
         CPABORT(error)
      END IF

      CALL timestop(handle)
   END SUBROUTINE bch_propagate

! **************************************************************************************************
!> \brief Updates the density in rtbse_env, using the provided exponential
!>        The new density is saved to a different matrix, which enables for comparison of matrices
!> \param rtbse_env Entry point of the calculation - contains current state of variables
!> \param exponential Real and imaginary parts ( + spin) of the exponential propagator
! **************************************************************************************************
   SUBROUTINE propagate_density(rtbse_env, exponential, rho_old, rho_new)
      TYPE(rtbse_env_type)                               :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: exponential, &
                                                            rho_old, &
                                                            rho_new
      CHARACTER(len=*), PARAMETER                        :: routineN = "propagate_density"
      INTEGER                                            :: j, handle

      CALL timeset(routineN, handle)
      IF (rtbse_env%mat_exp_method == do_exact) THEN
         ! For these methods, exponential is explicitly constructed
         DO j = 1, rtbse_env%n_spin
            ! rho * (exp^dagger)
            CALL parallel_gemm("N", "C", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                               CMPLX(1.0, 0.0, kind=dp), rho_old(j), exponential(j), &
                               CMPLX(0.0, 0.0, kind=dp), rtbse_env%rho_workspace(1))
            ! exp * rho * (exp^dagger)
            CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                               CMPLX(1.0, 0.0, kind=dp), exponential(j), rtbse_env%rho_workspace(1), &
                               CMPLX(0.0, 0.0, kind=dp), rho_new(j))
         END DO
      ELSE IF (rtbse_env%mat_exp_method == do_bch) THEN
         ! Same number of iterations as ETRS
         CALL bch_propagate(exponential, rho_old, rho_new, rtbse_env%rho_workspace, threshold_opt=rtbse_env%exp_accuracy, &
                            max_iter_opt=rtbse_env%etrs_max_iter)
      ELSE
         CPABORT("Only BCH and exact matrix exponentiation implemented.")
      END IF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Outputs the number of electrons in the system from the density matrix
!> \note  Moments matrix is provided by the rtbse_env, uses rho_workspace(1:3)
!> \param rtbse_env Entry point - rtbse environment
!> \param rho Density matrix in AO basis
!> \param electron_n_re Real number of electrons
!> \param electron_n_im Imaginary number of electrons, which can arise from numerical non-hermiticity
! **************************************************************************************************
   SUBROUTINE get_electron_number(rtbse_env, rho, electron_n_re, electron_n_im)
      TYPE(rtbse_env_type)                               :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho
      REAL(kind=dp), INTENT(OUT)                         :: electron_n_re, electron_n_im
      COMPLEX(kind=dp)                                   :: electron_n_buffer
      INTEGER                                            :: j

      electron_n_re = 0.0_dp
      electron_n_im = 0.0_dp
      CALL cp_fm_to_cfm(msourcer=rtbse_env%S_fm, mtarget=rtbse_env%rho_workspace(1))
      DO j = 1, rtbse_env%n_spin
         CALL cp_cfm_trace(rtbse_env%rho_workspace(1), rho(j), electron_n_buffer)
         electron_n_re = electron_n_re + REAL(electron_n_buffer, kind=dp)
         electron_n_im = electron_n_im + REAL(AIMAG(electron_n_buffer), kind=dp)
      END DO
      ! Scale by spin degeneracy
      electron_n_re = electron_n_re*rtbse_env%spin_degeneracy
      electron_n_im = electron_n_im*rtbse_env%spin_degeneracy
   END SUBROUTINE get_electron_number
! **************************************************************************************************
!> \brief Outputs the deviation from idempotence of density matrix
!> \note  Moments matrix is provided by the rtbse_env, uses rho_workspace(1:3)
!> \param rtbse_env Entry point - rtbse environment
!> \param rho Density matrix in AO basis
!> \param electron_n_re Real number of electrons
!> \param electron_n_im Imaginary number of electrons, which can arise from numerical non-hermiticity
! **************************************************************************************************
   SUBROUTINE get_idempotence_deviation(rtbse_env, rho, deviation_metric)
      TYPE(rtbse_env_type)                              :: rtbse_env
      TYPE(cp_cfm_type), DIMENSION(:), POINTER          :: rho
      REAL(kind=dp), INTENT(OUT)                        :: deviation_metric
      COMPLEX(kind=dp)                                  :: buffer_1, buffer_2
      REAL(kind=dp)                                     :: buffer_dev
      INTEGER                                           :: j

      deviation_metric = 0.0_dp
      buffer_dev = 0.0_dp
      ! First, determine Tr(S * rho_re) + i Tr (S * rho_im)
      CALL cp_fm_to_cfm(msourcer=rtbse_env%S_fm, mtarget=rtbse_env%rho_workspace(1))
      DO j = 1, rtbse_env%n_spin
         CALL cp_cfm_trace(rtbse_env%rho_workspace(1), rho(j), buffer_1)
         buffer_dev = buffer_dev + REAL(ABS(buffer_1)*ABS(buffer_1), kind=dp)
      END DO
      ! Now, determine Tr(S * rho_re * S * rho_re) - Tr(S * rho_im * S * rho_im) + 2i Tr(S * rho_re * S * rho_im)
      DO j = 1, rtbse_env%n_spin
         ! S * rho
         CALL multiply_fm_cfm("N", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                              1.0_dp, rtbse_env%S_fm, rho(j), &
                              0.0_dp, rtbse_env%rho_workspace(2))
         ! rho * S * rho
         CALL parallel_gemm("N", "N", rtbse_env%n_ao, rtbse_env%n_ao, rtbse_env%n_ao, &
                            CMPLX(1.0, 0.0, kind=dp), rho(j), rtbse_env%rho_workspace(2), &
                            CMPLX(0.0, 0.0, kind=dp), rtbse_env%rho_workspace(3))
         ! Tr (S * rho * S * rho)
         CALL cp_cfm_trace(rtbse_env%rho_workspace(1), rtbse_env%rho_workspace(3), buffer_2)
         deviation_metric = deviation_metric + REAL(ABS(buffer_2)*ABS(buffer_2), kind=dp)
      END DO
      deviation_metric = SQRT(deviation_metric) - SQRT(buffer_dev)
   END SUBROUTINE get_idempotence_deviation

! **************************************************************************************************
!> \brief Calculates the self-energy by contraction of screened potential, for complex density
!> \note Can be used for both the Coulomb hole part and screened exchange part
!> \param rtbse_env Quickstep environment data, entry point of the calculation
!> \param sigma_cfm Pointer to the self-energy full matrix, which is overwritten by this routine
!> \param greens_cfm Pointer to the Green's function matrix, which is used as input data
!> \author Stepan Marek
!> \date 09.2024
! **************************************************************************************************
   SUBROUTINE get_sigma_complex(rtbse_env, sigma_cfm, prefactor_opt, greens_cfm)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                                  :: sigma_cfm ! resulting self energy
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: prefactor_opt
      TYPE(cp_cfm_type), INTENT(IN)                      :: greens_cfm ! matrix to contract with RI_W
      REAL(kind=dp)                                      :: prefactor

      prefactor = 1.0_dp
      IF (PRESENT(prefactor_opt)) prefactor = prefactor_opt

      ! Carry out the sigma part twice
      ! Real part
      CALL cp_cfm_to_fm(msource=greens_cfm, mtargetr=rtbse_env%real_workspace(1))
      CALL get_sigma(rtbse_env, rtbse_env%real_workspace(2), prefactor, rtbse_env%real_workspace(1))
      CALL cp_fm_to_cfm(msourcer=rtbse_env%real_workspace(2), mtarget=rtbse_env%ham_workspace(1))
      ! Imaginary part
      CALL cp_cfm_to_fm(msource=greens_cfm, mtargeti=rtbse_env%real_workspace(1))
      CALL get_sigma(rtbse_env, rtbse_env%real_workspace(2), prefactor, rtbse_env%real_workspace(1))
      CALL cp_fm_to_cfm(msourcei=rtbse_env%real_workspace(2), mtarget=sigma_cfm)
      ! Add the real part
      CALL cp_cfm_scale_and_add(CMPLX(1.0, 0.0, kind=dp), sigma_cfm, &
                                CMPLX(1.0, 0.0, kind=dp), rtbse_env%ham_workspace(1))

   END SUBROUTINE get_sigma_complex
! **************************************************************************************************
!> \brief Calculates the self-energy by contraction of screened potential, for complex density
!> \note Can be used for both the Coulomb hole part and screened exchange part
!> \param rtbse_env Quickstep environment data, entry point of the calculation
!> \param sigma_fm Pointer to the self-energy full matrix, which is overwritten by this routine
!> \param greens_fm Pointer to the Green's function matrix, which is used as input data
!> \author Stepan Marek
!> \date 09.2024
! **************************************************************************************************
   SUBROUTINE get_sigma_real(rtbse_env, sigma_fm, prefactor_opt, greens_fm)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_fm_type)                                   :: sigma_fm ! resulting self energy
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: prefactor_opt
      TYPE(cp_fm_type), INTENT(IN)                       :: greens_fm ! matrix to contract with RI_W
      REAL(kind=dp)                                      :: prefactor

      prefactor = 1.0_dp
      IF (PRESENT(prefactor_opt)) prefactor = prefactor_opt

      ! Carry out the sigma part twice
      ! Convert to dbcsr
      CALL copy_fm_to_dbcsr(greens_fm, rtbse_env%rho_dbcsr)
      CALL get_sigma_dbcsr(rtbse_env, sigma_fm, prefactor, rtbse_env%rho_dbcsr)

   END SUBROUTINE get_sigma_real
! **************************************************************************************************
!> \brief Calculates the self-energy by contraction of screened potential
!> \note Can be used for both the Coulomb hole part and screened exchange part
!> \param greens_fm Pointer to the Green's function matrix, which is used as input data
!> \param sigma_fm Pointer to the self-energy full matrix, which is overwritten by this routine
!> \author Stepan Marek
!> \date 01.2024
! **************************************************************************************************
   SUBROUTINE get_sigma_dbcsr(rtbse_env, sigma_fm, prefactor_opt, greens_dbcsr)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_fm_type)                                   :: sigma_fm ! resulting self energy
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: prefactor_opt
      TYPE(dbcsr_type)                                   :: greens_dbcsr
      REAL(kind=dp)                                      :: prefactor

      prefactor = 1.0_dp
      IF (PRESENT(prefactor_opt)) prefactor = prefactor_opt

      CALL get_sigma_noenv(sigma_fm, prefactor_opt, greens_dbcsr, &
                           rtbse_env%screened_dbt, rtbse_env%t_3c_w, &
                           rtbse_env%t_3c_work_RI_AO__AO, rtbse_env%t_3c_work2_RI_AO__AO, &
                           rtbse_env%greens_dbt)
   END SUBROUTINE get_sigma_dbcsr
! **************************************************************************************************
!> \brief Calculates the self-energy by contraction of screened potential
!> \note Can be used for both the Coulomb hole part and screened exchange part
!> \note Separated from the rtbse_env - can be in principle called outside of the RTBSE code
!> \param sigma_fm Pointer to the self-energy full matrix, which is overwritten by this routine
!> \param prefactor_opt Optional argument for the prefactor (used for Coulomb hole calculation)
!> \param greens_dbcsr Matrix storing the lesser Green's function elements
!> \param screened_dbt Tensor storing the W_PQ screened Coulomb interaction RI matrix elements
!> \param int_3c_dbt Tensor storing the 3c integrals (RI| ORB ORB )
!> \param work_dbt_3c_1 Tensor workspace optimised for RI_AO__AO contractions
!> \param work_dbt_3c_2 Tensor workspace optimised for RI_AO__AO contractions
!> \param work_dbt_2c Tensor workspace for 2c integrals (Green's function and self-energy)
!> \author Stepan Marek
!> \date 01.2025
! **************************************************************************************************
   SUBROUTINE get_sigma_noenv(sigma_fm, prefactor_opt, greens_dbcsr, screened_dbt, &
                              int_3c_dbt, work_dbt_3c_1, work_dbt_3c_2, work_dbt_2c)
      TYPE(cp_fm_type)                                   :: sigma_fm ! resulting self energy
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: prefactor_opt
      TYPE(dbcsr_type)                                   :: greens_dbcsr
      TYPE(dbt_type)                                     :: screened_dbt, &
                                                            int_3c_dbt, &
                                                            work_dbt_3c_1, &
                                                            work_dbt_3c_2, &
                                                            work_dbt_2c
      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_sigma'
      REAL(kind=dp)                                      :: prefactor
      TYPE(dbcsr_type)                                   :: sigma_dbcsr
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (PRESENT(prefactor_opt)) THEN
         prefactor = prefactor_opt
      ELSE
         prefactor = 1.0_dp
      END IF

      ! Three-centre integrals are obtained from build_3c_integrals, from qs_tensors
      ! These should use sparcity, while W and Sigma can be full matrices
      ! The summation is carried out by dbt library - dbt_contract in dbt_api
      ! The building of the tensors might be a bit hard, because it requires a lot of parallel information
      ! Probably just use the tensors already present in bs_env? They seem to be mostly work tensors
      ! Create by template
      CALL dbt_contract(alpha=1.0_dp, &
                        tensor_1=screened_dbt, &
                        tensor_2=int_3c_dbt, &
                        beta=0.0_dp, &
                        tensor_3=work_dbt_3c_1, &
                        contract_1=[2], notcontract_1=[1], map_1=[1], &
                        contract_2=[1], notcontract_2=[2, 3], map_2=[2, 3])!,&
      !filter_eps=bs_env%eps_filter)
      ! t_work1 now contains B^P_(nu beta) = sum _ Q W _ (PQ) (iomega = 0) (Q| nu beta)
      ! Next step is to convert the greens full matrix to dbcsr matrix
      CALL dbt_copy_matrix_to_tensor(greens_dbcsr, work_dbt_2c)
      ! Then contract it
      ! no scaling applied - this has to be applied externally
      CALL dbt_contract(alpha=1.0_dp, &
                        tensor_1=work_dbt_3c_1, &
                        tensor_2=work_dbt_2c, &
                        beta=0.0_dp, &
                        tensor_3=work_dbt_3c_2, &
                        contract_1=[2], notcontract_1=[1, 3], map_1=[1, 3], &
                        contract_2=[2], notcontract_2=[1], map_2=[2])
      ! workspace 2 now contains C ^ P _ (mu beta) sum _ nu B ^ P _ (nu beta) g _ (mu nu)
      CALL dbt_contract(alpha=prefactor, &
                        tensor_1=int_3c_dbt, &
                        tensor_2=work_dbt_3c_2, &
                        beta=0.0_dp, &
                        tensor_3=work_dbt_2c, &
                        contract_1=[1, 3], notcontract_1=[2], map_1=[1], &
                        contract_2=[1, 2], notcontract_2=[3], map_2=[2])!,&
      !filter_eps=bs_env%eps_filter)
      ! Finally, convert the COH tensor to matrix and then to fm matrix
      ! TODO : extra workspace?
      CALL dbcsr_create(sigma_dbcsr, name="sigma", template=greens_dbcsr)
      CALL dbt_copy_tensor_to_matrix(work_dbt_2c, sigma_dbcsr)
      CALL copy_dbcsr_to_fm(sigma_dbcsr, sigma_fm)
      CALL dbcsr_release(sigma_dbcsr)
      ! Clear workspaces - saves memory?
      CALL dbt_clear(work_dbt_3c_1)
      CALL dbt_clear(work_dbt_3c_2)
      CALL dbt_clear(work_dbt_2c)
      CALL timestop(handle)

   END SUBROUTINE get_sigma_noenv
! **************************************************************************************************
!> \brief Creates the RI matrix and populates it with correct values
!> \note Tensor contains Hartree elements in the auxiliary basis
!> \param qs_env Quickstep environment - entry point of calculation
!> \author Stepan Marek
!> \date 01.2024
! **************************************************************************************************
   SUBROUTINE init_hartree(rtbse_env, v_dbcsr)
      TYPE(rtbse_env_type), POINTER, INTENT(IN)          :: rtbse_env
      TYPE(dbcsr_type)                                   :: v_dbcsr
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(libint_potential_type)                        :: coulomb_op
      TYPE(cp_fm_type)                                   :: V_fm
      TYPE(cp_fm_type)                                   :: metric_fm
      TYPE(cp_fm_type)                                   :: metric_inv_fm, &
                                                            work_fm
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: V_dbcsr_a, &
                                                            metric_dbcsr
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c

      bs_env => rtbse_env%bs_env

      ! Allocate for bare Hartree term
      ALLOCATE (V_dbcsr_a(1))
      ALLOCATE (metric_dbcsr(1))
      CALL dbcsr_create(V_dbcsr_a(1), name="Hartree_dbcsr", template=bs_env%mat_RI_RI%matrix)
      CALL dbcsr_create(metric_dbcsr(1), name="RI_metric_dbcsr", template=bs_env%mat_RI_RI%matrix)

      ! Calculate full coulomb RI basis elements - V _ (PQ) matrix
      NULLIFY (nl_2c)
      CALL build_2c_neighbor_lists(nl_2c, bs_env%basis_set_RI, bs_env%basis_set_RI, &
                                   coulomb_op, "Coulomb_neighbor_2c_list", rtbse_env%qs_env, &
                                   sym_ij=.FALSE., molecular=.TRUE.)
      CALL build_2c_integrals(V_dbcsr_a, bs_env%eps_filter, rtbse_env%qs_env, nl_2c, &
                              bs_env%basis_set_RI, bs_env%basis_set_RI, coulomb_op, &
                              do_kpoints=.FALSE., regularization_RI=bs_env%regularization_RI)
      ! Calculate the RI metric elements
      ! nl_2c is automatically rewritten (even reallocated) in this routine
      CALL build_2c_neighbor_lists(nl_2c, bs_env%basis_set_RI, bs_env%basis_set_RI, &
                                   bs_env%ri_metric, "Metric_neighbor_2c_list", rtbse_env%qs_env, &
                                   sym_ij=.FALSE., molecular=.TRUE.)
      CALL build_2c_integrals(metric_dbcsr, bs_env%eps_filter, rtbse_env%qs_env, nl_2c, &
                              bs_env%basis_set_RI, bs_env%basis_set_RI, bs_env%ri_metric, &
                              do_kpoints=.FALSE., regularization_RI=bs_env%regularization_RI)
      ! nl_2c no longer needed
      CALL release_neighbor_list_sets(nl_2c)
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
      CALL cp_dbcsr_sm_fm_multiply(V_dbcsr_a(1), metric_inv_fm, &
                                   work_fm, bs_env%n_RI)
      CALL parallel_gemm("N", "N", bs_env%n_RI, bs_env%n_RI, bs_env%n_RI, &
                         1.0_dp, metric_inv_fm, work_fm, 0.0_dp, V_fm)
      ! Now, create the tensor from the matrix
      ! First, convert full matrix to dbcsr
      CALL dbcsr_clear(V_dbcsr_a(1))
      CALL copy_fm_to_dbcsr(V_fm, V_dbcsr_a(1))
      CALL dbcsr_create(v_dbcsr, "Hartree ri", V_dbcsr_a(1))
      CALL dbcsr_copy(v_dbcsr, V_dbcsr_a(1))
      ! Create and copy distinctly, so that unnecessary objects can be destroyed
      ! Destroy all unnecessary matrices
      CALL dbcsr_release(V_dbcsr_a(1))
      CALL dbcsr_release(metric_dbcsr(1))
      DEALLOCATE (V_dbcsr_a)
      DEALLOCATE (metric_dbcsr)
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
!> \param rtbse_env Entry point of GWBSE - uses rho_dbcsr and some complex_workspace
!> \param rho_ao Density matrix in ao basis
!> \param v_ao Overwritten by the Hartree matrix in the atomic orbital basis
!> \author Stepan Marek
!> \date 01.2025
! **************************************************************************************************
   SUBROUTINE get_hartree_env(rtbse_env, rho_fm, v_fm)
      TYPE(rtbse_env_type), POINTER                      :: rtbse_env
      TYPE(cp_cfm_type)                                  :: rho_fm
      TYPE(cp_fm_type)                                   :: v_fm
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CALL get_qs_env(rtbse_env%qs_env, para_env=para_env, bs_env=bs_env)

      CALL get_hartree_noenv(v_fm, rho_fm, rtbse_env%int_3c_array, rtbse_env%v_dbcsr, &
                             rtbse_env%n_RI, bs_env%sizes_RI, &
                             para_env, rtbse_env%rho_dbcsr, rtbse_env%v_ao_dbcsr)
   END SUBROUTINE get_hartree_env
! **************************************************************************************************
!> \brief Calculates the Hartree matrix in the atomic orbital basis, given a density matrix, in local arrays
!>        Calculates the values for single spin species present in given rho
!> \param v_fm Hartree potential in atomic orbital basis - is overwritten by the updated potential
!> \param rho_fm Density matrix corresponding to single spin species, in atomic orbital basis
!> \param int_3c Previously allocated array (best to use create_hartree_ri_3c) for 3c integrals
!> \param v_dbcsr Previously calculated 2c Coulomb repulsion between RI orbitals
!> \param n_RI Number of RI basis orbitals
!> \param sizes_RI Number of RI basis orbitals per atom
!> \param para_env MPI Parallel environment (used for summation across ranks)
!> \param rho_dbcsr Previously created dbcsr matrix, used as workspace
!> \param v_ao_dbcsr Previously created dbcsr matrix, used as workspace
!> \author Stepan Marek
!> \date 01.2025
! **************************************************************************************************
   SUBROUTINE get_hartree_noenv(v_fm, rho_fm, int_3c, v_dbcsr, n_RI, sizes_RI, para_env, rho_dbcsr, v_ao_dbcsr)
      TYPE(cp_fm_type)                                   :: v_fm
      TYPE(cp_cfm_type), INTENT(IN)                      :: rho_fm
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: int_3c
      TYPE(dbcsr_type)                                   :: v_dbcsr
      INTEGER                                            :: n_RI
      INTEGER, DIMENSION(:)                              :: sizes_RI
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type)                                   :: rho_dbcsr, v_ao_dbcsr
      CHARACTER(len=*), PARAMETER                        :: routineN = "get_hartree"
      TYPE(dbcsr_iterator_type)                          :: iterator_matrix
      INTEGER                                            :: i, j, k, n, nblocks, ind_1, ind_2, row_offset, col_offset, &
                                                            row_size, col_size, j_n_AO, k_n_AO, i_n_RI, &
                                                            ri_offset, ind_i, handle
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: Pvector, Qvector
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: block_matrix
      INTEGER                                            :: nblkrows_local, nblkcols_local, j_blk, k_blk, j_offset, k_offset
      INTEGER, DIMENSION(:), POINTER                     :: local_blk_rows, local_blk_cols
      LOGICAL                                            :: found

      MARK_USED(i_n_RI)
      MARK_USED(ri_offset)
      MARK_USED(ind_i)

      ! No memory optimisation so far - calculate all 3cs an all ranks
      ! Importantly - dbcsr blocks are ordered by atoms - i.e. ethene with 6 atoms will have 6x6 block structure
      ! Number of basis states on each basis set is known is post_scf_bandstructure env

      CALL timeset(routineN, handle)

      ! Allocate the Q and Pvector on each rank
      ALLOCATE (Qvector(n_RI), source=0.0_dp)
      ALLOCATE (Pvector(n_RI), source=0.0_dp)

      ! First step - analyze the structure of copied dbcsr matrix on all ranks
      CALL dbcsr_clear(rho_dbcsr)
      ! Only the real part of the density matrix contributes
      ! Use v_fm as workspace
      CALL cp_cfm_to_fm(msource=rho_fm, mtargetr=v_fm)
      CALL copy_fm_to_dbcsr(v_fm, rho_dbcsr)
      j_offset = 0
      CALL dbcsr_get_info(rho_dbcsr, nblkrows_local=nblkrows_local, nblkcols_local=nblkcols_local, &
                          local_rows=local_blk_rows, local_cols=local_blk_cols)
      DO j_blk = 1, nblkrows_local
         k_offset = 0
         DO k_blk = 1, nblkcols_local
            ! Check whether we can retrieve the rho block
            ! TODO : Handle transposed case?
            CALL dbcsr_get_block_p(rho_dbcsr, local_blk_rows(j_blk), local_blk_cols(k_blk), &
                                   block=block_matrix, found=found, row_size=row_size, col_size=col_size)
            ! If the block is not found, then the density matrix here has below threshold values
            IF (.NOT. found) CYCLE
            ! With the block retrieved, add its contributions to the Q-vector
            !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j,k) &
            !$OMP SHARED(n_RI, row_size, col_size, Qvector, int_3c, j_offset, k_offset, block_matrix)
            DO i = 1, n_RI
               DO j = 1, row_size
                  DO k = 1, col_size
                     Qvector(i) = Qvector(i) + int_3c(j_offset + j, k_offset + k, i)*block_matrix(j, k)
                  END DO
               END DO
            END DO
            !$OMP END PARALLEL DO
            ! Increment k-offset - setup for the next block
            k_offset = k_offset + col_size
         END DO
         ! Increments the j-offset - row_size is carried over from the last iteration
         j_offset = j_offset + row_size
      END DO
      ! Now, each rank has contributions from D_jk within its scope
      ! Need to sum over different ranks to get the total vector on all ranks
      CALL para_env%sum(Qvector)
      ! Once this is done, Pvector is current on all ranks
      ! Continue with V_PQ summation
      nblocks = dbcsr_get_num_blocks(v_dbcsr)
      CALL dbcsr_iterator_start(iterator_matrix, v_dbcsr)
      DO n = 1, nblocks
         ! TODO : Try OMP parallelisation over different blocks - expect many more available speedup for large systems
         CALL dbcsr_iterator_next_block(iterator_matrix, ind_1, ind_2, block_matrix, &
                                        row_offset=row_offset, col_offset=col_offset, row_size=row_size, col_size=col_size)
         ! TODO : Better names for RI
         j_n_AO = sizes_RI(ind_1)
         k_n_AO = sizes_RI(ind_2)
         ! The allocations are as follows
         !$OMP PARALLEL DO DEFAULT(none) PRIVATE(j,k) &
         !$OMP SHARED(block_matrix, Pvector, Qvector,j_n_AO,k_n_AO,row_offset,col_offset)
         DO j = 1, j_n_AO
            DO k = 1, k_n_AO
               Pvector(j + row_offset - 1) = Pvector(j + row_offset - 1) + block_matrix(j, k)*Qvector(k + col_offset - 1)
            END DO
         END DO
         !$OMP END PARALLEL DO
      END DO
      CALL dbcsr_iterator_stop(iterator_matrix)
      ! Again, make sure that the P vector is present on all ranks
      CALL para_env%sum(Pvector)
      ! Now, for the final trick, iterate over local blocks of v_ao_dbcsr to get the Hartree as dbcsr, then convert to fm
      ! TODO : Clear or set blocks to zero
      ! CALL dbcsr_clear(v_ao_dbcsr)
      j_offset = 0
      DO j_blk = 1, nblkrows_local
         k_offset = 0
         DO k_blk = 1, nblkcols_local
            ! Check whether we can retrieve the rho block
            ! TODO : Handle transposed case?
            CALL dbcsr_get_block_p(v_ao_dbcsr, local_blk_rows(j_blk), local_blk_cols(k_blk), &
                                   block=block_matrix, found=found, row_size=row_size, col_size=col_size)
            ! If the block is not found, reserve it
            IF (.NOT. found) THEN
               ! Reservations
               CALL dbcsr_reserve_blocks(v_ao_dbcsr, local_blk_rows(j_blk:j_blk), local_blk_cols(k_blk:k_blk))
               ! Rerun the getter to get the new block
               CALL dbcsr_get_block_p(v_ao_dbcsr, local_blk_rows(j_blk), local_blk_cols(k_blk), &
                                      block=block_matrix, found=found, row_size=row_size, col_size=col_size)
            END IF
            ! With the block retrieved, contract with the P vector
            !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j,k) &
            !$OMP SHARED(row_size, col_size, n_RI, block_matrix, Pvector, int_3c, j_offset, k_offset)
            DO j = 1, row_size
               DO k = 1, col_size
                  block_matrix(j, k) = 0.0_dp
                  DO i = 1, n_RI
                     block_matrix(j, k) = block_matrix(j, k) + Pvector(i)*int_3c(j_offset + j, k_offset + k, i)
                  END DO
               END DO
            END DO
            !$OMP END PARALLEL DO
            ! Increment k-offset - setup for the next block
            k_offset = k_offset + col_size
         END DO
         ! Increments the j-offset - row_size is carried over from the last iteration
         j_offset = j_offset + row_size
      END DO
      ! Since P vector was present on all the ranks, v_dbcsr_ao has the complete Hartree result
      ! copy_dbcsr_to_fm should set all values in v_fm to zero
      CALL copy_dbcsr_to_fm(v_ao_dbcsr, v_fm)
      DEALLOCATE (Qvector)
      DEALLOCATE (Pvector)

      CALL timestop(handle)
   END SUBROUTINE get_hartree_noenv
! **************************************************************************************************
!> \brief Calculates the exponential of a matrix in a generalized eigenvalue problem. Specifically,
!>        it assumes we have a Hermitian matrix A in the eigenvalue problem AX = BXE, where B is some overlap
!>        matrix and E is a diagonal matrix of real eigenvalues. Then, it calculates
!>        exp(B^(-1) A) = X exp(E) X^C B
!> \param amatrix Matrix to exponentiate
!> \param bmatrix Overlap matrix
!> \param exponential Exponential exp(B^(-1) A) is stored here after the routine is finished
!> \param eig_scale_opt Optionally scale eigenvalues by a complex number before exponentiating them
!> \param work_opt Optionally provide workspace (of size at least 4) that is used in the calculation
!> \author Stepan Marek
!> \date 09.2024
! **************************************************************************************************
   SUBROUTINE cp_cfm_gexp(amatrix, bmatrix, exponential, eig_scale_opt, work_opt)
      ! TODO : Do interface for real matrices
      TYPE(cp_cfm_type), INTENT(IN)                      :: amatrix
      TYPE(cp_cfm_type), INTENT(IN)                      :: bmatrix
      TYPE(cp_cfm_type)                                  :: exponential
      COMPLEX(kind=dp), INTENT(IN), OPTIONAL             :: eig_scale_opt
      TYPE(cp_cfm_type), DIMENSION(:), POINTER, OPTIONAL :: work_opt
      CHARACTER(len=*), PARAMETER                        :: routineN = "cp_cfm_gexp"
      COMPLEX(kind=dp)                                   :: eig_scale
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: eigenvalues
      COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE        :: expvalues
      TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: work
      LOGICAL                                            :: deallocate_work
      INTEGER                                            :: nrow, i, handle

      CALL timeset(routineN, handle)

      ! Argument parsing and sanity checks
      IF (PRESENT(eig_scale_opt)) THEN
         eig_scale = eig_scale_opt
      ELSE
         eig_scale = CMPLX(1.0, 0.0, kind=dp)
      END IF

      NULLIFY (work)
      IF (PRESENT(work_opt) .AND. SIZE(work_opt) >= 4) THEN
         deallocate_work = .FALSE.
         work => work_opt
      ELSE
         deallocate_work = .TRUE.
         ALLOCATE (work(4))
         ! Allocate the work storage on the fly
         DO i = 1, 4
            CALL cp_cfm_create(work(i), amatrix%matrix_struct)
         END DO
      END IF

      nrow = amatrix%matrix_struct%nrow_global

      ALLOCATE (eigenvalues(nrow))
      ALLOCATE (expvalues(nrow))

      ! Do not change the amatrix and bmatrix - need to copy them first
      CALL cp_cfm_to_cfm(amatrix, work(1))
      CALL cp_cfm_to_cfm(bmatrix, work(2))

      ! Solve the generalized eigenvalue equation
      CALL cp_cfm_geeig(work(1), work(2), work(3), eigenvalues, work(4))

      ! Scale and exponentiate the eigenvalues
      expvalues(:) = EXP(eigenvalues(:)*eig_scale)

      ! Copy eigenvectors to column scale them
      CALL cp_cfm_to_cfm(work(3), work(1))
      ! X * exp(E)
      CALL cp_cfm_column_scale(work(1), expvalues)

      ! Carry out the remaining operations
      ! X * exp(E) * X^C
      CALL parallel_gemm("N", "C", nrow, nrow, nrow, &
                         CMPLX(1.0, 0.0, kind=dp), work(1), work(3), &
                         CMPLX(0.0, 0.0, kind=dp), work(2))
      ! X * exp(E) * X^C * B
      CALL parallel_gemm("N", "N", nrow, nrow, nrow, &
                         CMPLX(1.0, 0.0, kind=dp), work(2), bmatrix, &
                         CMPLX(0.0, 0.0, kind=dp), exponential)

      ! Deallocate work storage if necessary
      IF (deallocate_work) THEN
         DO i = 1, 4
            CALL cp_cfm_release(work(i))
         END DO
         DEALLOCATE (work)
      END IF

      DEALLOCATE (eigenvalues)
      DEALLOCATE (expvalues)

      CALL timestop(handle)
   END SUBROUTINE cp_cfm_gexp
END MODULE rt_bse

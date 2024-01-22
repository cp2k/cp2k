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
   USE kinds,                           ONLY: dp
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_type,&
                                              dbcsr_print,&
                                              dbcsr_has_symmetry,&
                                              dbcsr_desymmetrize,&
                                              dbcsr_create,&
                                              dbcsr_release,&
                                              dbcsr_copy,&
                                              dbcsr_set,&
                                              dbcsr_clear,&
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
                                              dbt_contract,&
                                              dbt_copy_matrix_to_tensor,&
                                              dbt_copy_tensor_to_matrix,&
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
   USE qs_tensors,                      ONLY: neighbor_list_3c_destroy
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_invert,&
                                              cp_fm_scale_and_add
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_moments,                      ONLY: build_local_moment_matrix
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
                                              rtp_tdagw_ham_g0w0

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_tdagw"

   PUBLIC :: run_propagation_gw

CONTAINS

! **************************************************************************************************
!> \brief Runs the electron-only real time adiabatic GW (TD Bethe-Salpeter) propagation
!> \param qs_env Quickstep environment data, entry point of the calculation
! **************************************************************************************************
   SUBROUTINE run_propagation_gw(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(len=*), PARAMETER                        :: routineN = 'run_propagation_gw'
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input, rtp_section
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type), POINTER                          :: S_inv
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: start_mos
      TYPE(cp_fm_type), DIMENSION(:), ALLOCATABLE, TARGET:: moments_full,&
                                                            ham_single_particle,&
                                                            rho_full
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: moments_full_p,&
                                                            ham_single_particle_p,&
                                                            rho_full_p
      TYPE(cp_fm_type), DIMENSION(:), ALLOCATABLE        :: ham_effective,&
                                                            ham_exp,&
                                                            rho_work
      ! TODO : Maybe complex?
      TYPE(cp_fm_type), TARGET                           :: sigma_COH,&
                                                            sigma_SEX,&
                                                            S_inv_fm,&
                                                            hartree_fm,&
                                                            hartree_fm_orig
      TYPE(cp_fm_type), POINTER                          :: sigma_p,&
                                                            greens_fm_p,&
                                                            hartree_fm_p
      TYPE(dbt_type), TARGET, ALLOCATABLE                :: hartree_dbt
      TYPE(dbt_type), POINTER                            :: hartree_dbt_p
      REAL(kind=dp)                                      :: field(3)
      INTEGER                                            :: i,j,k,nspin,re,im
      INTEGER                                            :: rho_unit

      ! TODO : Ask Jan whether there is already an identifier present in some of the environments
      CALL get_qs_env(qs_env, mos=start_mos, bs_env=bs_env)
      nspin = SIZE(start_mos)
      ! ********************** ALLOCATIONS **********************

      ! Allocate the single particle Hamiltonian - can be either KS or G0W0 Hamiltonian
      ! Marek TODO : Always uses bs_env - truly only for debug or redo with possibility to refer to qs_env?
      ALLOCATE(ham_single_particle(nspin))
      DO i=1,nspin
         CALL cp_fm_create(ham_single_particle(i), bs_env%fm_h_G0W0_Gamma%matrix_struct)
      END DO
      ham_single_particle_p => ham_single_particle

      ! Allocate the moments matrices
      ALLOCATE(moments_full(3))
      ! Now, for each index, create the matrix
      DO i=1,3
         ! Use cp_fm_struct from Kohn-Sham matrix TODO : Differentiate by spin/kp? - not necessary for moments 
         CALL cp_fm_create(moments_full(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
      END DO
      moments_full_p => moments_full

      ! Allocate the storage for effective Hamiltonian, density and propagator matrix - Hamiltonian exponential
      ALLOCATE(ham_effective(2*nspin))
      ALLOCATE(rho_full(2*nspin))
      ALLOCATE(ham_exp(2*nspin))
      DO i=1,2*nspin
         CALL cp_fm_create(ham_effective(i), bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_create(rho_full(i), bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_create(ham_exp(i), bs_env%fm_h_G0W0_Gamma%matrix_struct)
      END DO
      rho_full_p => rho_full

      ! Allocate and initialise workspace for density propagation
      ALLOCATE(rho_work(3))
      DO i=1,3
         CALL cp_fm_create(rho_work(i), bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_set_all(rho_work(i), 0.0_dp)
      END DO

      ! Allocate sigma matrices
      ! TODO : Spin dependence
      CALL cp_fm_create(sigma_COH, bs_env%fm_h_G0W0_Gamma%matrix_struct)
      CALL cp_fm_create(sigma_SEX, bs_env%fm_h_G0W0_Gamma%matrix_struct)
      CALL cp_fm_set_all(sigma_COH, 0.0_dp)
      CALL cp_fm_set_all(sigma_SEX, 0.0_dp)

      ! Allocate inverse full matrix
      CALL cp_fm_create(S_inv_fm, bs_env%fm_h_G0W0_Gamma%matrix_struct)
      CALL cp_fm_set_all(S_inv_fm, 0.0_dp)

      ! Allocate structures for the Hartree term
      CALL cp_fm_create(hartree_fm, bs_env%fm_h_G0W0_Gamma%matrix_struct)
      CALL cp_fm_set_all(hartree_fm, 0.0_dp)
      CALL cp_fm_create(hartree_fm_orig, bs_env%fm_h_G0W0_Gamma%matrix_struct)
      CALL cp_fm_set_all(hartree_fm_orig, 0.0_dp)
      ALLOCATE(hartree_dbt)
      hartree_dbt_p => hartree_dbt
      CALL init_hartree(qs_env, hartree_dbt_p)
      PRINT *, "Hartree initialised"
      ! ********************** EXTRACT PREVIOUS DATA **********************
      ! Prepare the relevant matrices
      CALL tdagw_prepare(qs_env, ham_single_particle_p, moments_full_p, rho_full_p)

      ! Get data from QS and RTP
      CALL get_qs_env(qs_env,&
                      rtp=rtp,&
                      bs_env=bs_env,&
                      dft_control=dft_control,&
                      input=input)
      CALL get_rtp(rtp, S_inv=S_inv)

      ! TODO : Determine the Coulomb-Hole part of the self-energy - remains constant throughout the time propagation
      ! TODO : Maybe time has come to create a structure that initializes the necessary storage?
      CALL copy_dbcsr_to_fm(S_inv, S_inv_fm)
      ! Point to correct structures
      sigma_p => sigma_COH
      greens_fm_p => S_inv_fm
      CALL get_sigma(qs_env, greens_fm_p, sigma_p, -0.5_dp)
      ! TODO : Maybe enable printing under debug option?
      ! PRINT *, "== Sigma COH =="
      ! CALL cp_fm_write_formatted(sigma_COH, 6)

      ! Subtract the v_xc from COH part of the self-energy, as V_xc is also not updated during the timestepping
      ! TODO : Spin dependence in sigma_COH
      CALL cp_fm_scale_and_add(1.0_dp, sigma_COH, -1.0_dp, bs_env%fm_V_xc_Gamma(1))

      ! Calculate the initial Hartree interaction matrix
      ! TODO : k-points and the possible divergence
      hartree_fm_p => hartree_fm_orig
      ! TODO : Spin dependence
      CALL get_hartree(qs_env, hartree_dbt, rho_full(1), hartree_fm_p)
      PRINT *, "Hartree retrieved"
      ! From now on, hartree_fm_p will point to current hartree potential
      hartree_fm_p => hartree_fm

      ! Setup the files
      logger => cp_get_default_logger()
      rtp_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION")
      rho_unit = cp_print_key_unit_nr(logger, rtp_section, "PRINT%DENSITY_MATRIX", extension=".dat")
      ! ********************** Start the time loop **********************
      DO i = rtp%i_start, rtp%nsteps

         ! Reset the effective Hamiltonian to KS Hamiltonian + G0W0
         DO j=1,nspin
            re = 2*j - 1
            im = 2*j
            CALL cp_fm_to_fm(ham_single_particle(j), ham_effective(re))
            ! Imaginary part of KS + G0W0 is zero
            CALL cp_fm_set_all(ham_effective(im), 0.0_dp)
         END DO

         ! Update the field
         qs_env%sim_time = REAL(i, dp)*rtp%dt
         ! TODO : Understand and implement rtp%istep - that should mainly be used for restart calculations?
         qs_env%sim_step = i
         CALL make_field(dft_control, field, qs_env%sim_step, qs_env%sim_time)
         PRINT *, "Field : ", field(1), field(2), field(3)

         ! TODO : Create the screened-exchange part of COHSEX - updated every loop
         !        Some W * G multiplication, through the auxiliary basis
         ! TODO : Preprocessor macro that sets up the spin loop
         DO j=1,nspin
            re = 2*j - 1
            im = 2*j
            ! Add electric field
            ! The field is added to the real part of both spins
            DO k=1,3
               CALL cp_fm_scale_and_add(1.0_dp, ham_effective(re), field(k), moments_full(k))
            END DO

            ! TODO : Add the Hartree term difference to the Hamiltonian - discuss with Jan
            ! TODO : Add the COHSEX self-energy to the Hamiltonian
            sigma_p => sigma_SEX
            ! Add the COH part
            ! TODO : COH spin dependence
            ! TODO : Proper inverse matric tensor multiplication
            ! TODO : Check the signs with respect to G -> rho transformation
            CALL cp_fm_scale_and_add(1.0_dp, ham_effective(re), 1.0_dp, sigma_COH)
            ! Calculate the SEX part
            ! TODO : SEX spin dependence
            ! Real part of iGW is product of Im G, and G = i rho, so Im G = Re rho, all multiplied by i^2 = -1
            greens_fm_p => rho_full(re)
            CALL get_sigma(qs_env, greens_fm_p, sigma_p, 1.0_dp)
            CALL cp_fm_scale_and_add(1.0_dp, ham_effective(re), -1.0_dp, sigma_SEX)
            ! Imaginary part of iGW is product of Re G = - Im rho
            greens_fm_p => rho_full(im)
            CALL get_sigma(qs_env, greens_fm_p, sigma_p, 1.0_dp)
            CALL cp_fm_scale_and_add(1.0_dp, ham_effective(im), -1.0_dp, sigma_SEX)
            ! Calculate update Hartree potential
            ! TODO : Spin (in)dependence check
            CALL get_hartree(qs_env, hartree_dbt, rho_full(re), hartree_fm_p)
            CALL cp_fm_scale_and_add(1.0_dp, ham_effective(re), 1.0_dp, hartree_fm)
            CALL get_hartree(qs_env, hartree_dbt, rho_full(im), hartree_fm_p)
            CALL cp_fm_scale_and_add(1.0_dp, ham_effective(im), 1.0_dp, hartree_fm)
            ! Subtract original Hartree, which is purely real
            CALL cp_fm_scale_and_add(1.0_dp, ham_effective(re), -1.0_dp, hartree_fm_orig)

            ! In order to produce correct result, need to remultiply by inverse overlap matrix
            ! TODO : Using rho_work(1) - maybe rename?
            CALL cp_dbcsr_sm_fm_multiply(S_inv, ham_effective(re), rho_work(1), bs_env%n_ao)
            CALL cp_dbcsr_sm_fm_multiply(S_inv, ham_effective(im), rho_work(2), bs_env%n_ao)
            CALL cp_fm_to_fm(rho_work(1),ham_effective(re))
            CALL cp_fm_to_fm(rho_work(2),ham_effective(im))
         END DO

         ! TODO : Create the matrix exponential (matrix_exponential.F:taylor_only_imaginary, but might need complex because of
         ! complex W and hence complex self-energy? Furthermore, is the self-energy in COHSEX Hermitian or not?)
         ! TODO : Understand which input options set the orders
         DO j=1,nspin
            re = 2*j - 1
            im = 2*j
            ! The evolution of density matrix is in the opposite way to the evolution of operators, i.e. multiply by -dt
            ! TODO : Maybe better to do opposite exponentiation
            CALL cp_fm_scale(-rtp%dt, ham_effective(re))
            ! Also scale by -1 as we are also multiplying by imaginary unit - this is convention used by exponential routines
            ! TODO : Note that so far, this is only multiplying a zero matrix
            CALL cp_fm_scale(rtp%dt, ham_effective(im))
         END DO
         DO j=1,nspin
            re = 2*j - 1
            im = 2*j
            ! CALL taylor_only_imaginary(ham_exp(re:im), ham_effective(re), rtp%orders(1,j), rtp%orders(2,j))
            CALL taylor_full_complex(ham_exp(re:im), ham_effective(im), ham_effective(re), rtp%orders(1,j), rtp%orders(2,j))
         END DO

         ! TODO : Consult/think about using complex routines and types
         ! (ham_exp(re) + i * ham_exp(im)) * (rho_full(re) + i * rho_full(im)) * (ham_exp(re)^T - i * ham_exp(im)^T) = 
         ! + ham_exp(re) * (rho_full(re) * ham_exp(re)^T + rho_full(im) * ham_exp(im)^T) | Term 1
         ! + ham_exp(im) * (rho_full(re) * ham_exp(im)^T - rho_full(im) * ham_exp(re)^T) | Term 2
         ! + i * ham_exp(im) * (rho_full(re) * ham_exp(re)^T + rho_full(im) * ham_exp(im)^T) | Term 1
         ! - i * ham_exp(re) * (rho_full(re) * ham_exp(im)^T - rho_full(im) * ham_exp(re)^T) | Term 2
         DO j=1,nspin
            re = 2*j - 1
            im = 2*j
            ! Real part
            ! Term 1 - right multiplication first
            CALL parallel_gemm("N","T", bs_env%n_ao, bs_env%n_ao, bs_env%n_ao,&
                               1.0_dp, rho_full(re), ham_exp(re), 0.0_dp, rho_work(1))
            CALL parallel_gemm("N","T", bs_env%n_ao, bs_env%n_ao, bs_env%n_ao,&
                               1.0_dp, rho_full(im), ham_exp(im), 1.0_dp, rho_work(1))
            ! Now term 1 left multiplication real part
            CALL parallel_gemm("N", "N", bs_env%n_ao, bs_env%n_ao, bs_env%n_ao,&
                               1.0_dp, ham_exp(re), rho_work(1), 0.0_dp, rho_work(2))
            ! Now term 1 left multiplication imag part
            CALL parallel_gemm("N", "N", bs_env%n_ao, bs_env%n_ao, bs_env%n_ao,&
                               1.0_dp, ham_exp(im), rho_work(1), 0.0_dp, rho_work(3))
            ! Repeat the same procedure for the second term
            ! Term 2 - right multiplication first
            CALL parallel_gemm("N","T", bs_env%n_ao, bs_env%n_ao, bs_env%n_ao,&
                               1.0_dp, rho_full(re), ham_exp(im), 0.0_dp, rho_work(1))
            CALL parallel_gemm("N","T",bs_env%n_ao, bs_env%n_ao, bs_env%n_ao,&
                               -1.0_dp, rho_full(im), ham_exp(re), 1.0_dp, rho_work(1))
            ! Now term 2 left multiplication real part
            CALL parallel_gemm("N", "N", bs_env%n_ao, bs_env%n_ao, bs_env%n_ao,&
                               1.0_dp, ham_exp(im), rho_work(1), 1.0_dp, rho_work(2))
            ! Now term 2 left multiplication imag part
            CALL parallel_gemm("N","N", bs_env%n_ao, bs_env%n_ao, bs_env%n_ao,&
                               -1.0_dp, ham_exp(re), rho_work(1), 1.0_dp, rho_work(3))
            ! Finally, copy the results to rho_full - not used anymore TODO : Unless there is spin coupling
            CALL cp_fm_to_fm(rho_work(2), rho_full(re))
            CALL cp_fm_to_fm(rho_work(3), rho_full(im))
            ! Transform the density matrix into molecular orbitals basis and print it out
            ! TODO : Number of molecular orbitals and atomic orbitals might not be the same
            !  - Store results in rho_work for now, maybe MO x MO matrix for results is not too expensive?
            ! S * rho
            CALL parallel_gemm("N","N",bs_env%n_ao,bs_env%n_ao,bs_env%n_ao,&
                            1.0_dp, bs_env%fm_s_Gamma, rho_full(re), 0.0_dp, rho_work(1))
            CALL parallel_gemm("N","N",bs_env%n_ao,bs_env%n_ao,bs_env%n_ao,&
                            1.0_dp, bs_env%fm_s_Gamma, rho_full(im), 0.0_dp, rho_work(2))
            ! C^T * S * rho
            CALL parallel_gemm("T","N",bs_env%n_ao,bs_env%n_ao,bs_env%n_ao,&
                            1.0_dp,bs_env%fm_mo_coeff_Gamma(j),rho_work(1),0.0_dp,rho_work(3))
            CALL parallel_gemm("T","N",bs_env%n_ao,bs_env%n_ao,bs_env%n_ao,&
                            1.0_dp,bs_env%fm_mo_coeff_Gamma(j),rho_work(2),0.0_dp,rho_work(1))
            ! C^T * S * rho * S
            CALL parallel_gemm("N","N",bs_env%n_ao,bs_env%n_ao,bs_env%n_ao,&
                            1.0_dp, rho_work(3),bs_env%fm_s_Gamma, 0.0_dp, rho_work(2))
            CALL parallel_gemm("N","N",bs_env%n_ao,bs_env%n_ao,bs_env%n_ao,&
                            1.0_dp, rho_work(1),bs_env%fm_s_Gamma, 0.0_dp, rho_work(3))
            ! C^T * S * rho * S * C
            CALL parallel_gemm("N","N",bs_env%n_ao,bs_env%n_ao,bs_env%n_ao,&
                            1.0_dp,rho_work(2),bs_env%fm_mo_coeff_Gamma(j),0.0_dp,rho_work(1))
            CALL parallel_gemm("N","N",bs_env%n_ao,bs_env%n_ao,bs_env%n_ao,&
                            1.0_dp,rho_work(3),bs_env%fm_mo_coeff_Gamma(j),0.0_dp,rho_work(2))
            ! TODO : Save these as separate files
            WRITE(UNIT=rho_unit, FMT="(A,I6,A)") "#Iter : ", i, " Real"
            CALL cp_fm_write_formatted(rho_work(1), rho_unit)
            WRITE(UNIT=rho_unit, FMT="(A,I6,A)") "#Iter : ", i, " Imag"
            CALL cp_fm_write_formatted(rho_work(2), rho_unit)
         END DO
      END DO
      ! ********************** End the time loop **********************

      ! ********************** DEALLOCATIONS **********************
      DO i=1,3
         CALL cp_fm_release(moments_full(i))
      END DO
      DEALLOCATE(moments_full)
      DO i=1,2*nspin
         ! Deallocate matrix used to store the propagator
         CALL cp_fm_release(ham_exp(i))
         ! Deallocate matrix used to store the effective Hamiltonian
         CALL cp_fm_release(ham_effective(i))
         ! Deallocate the propagated density matrix
         CALL cp_fm_release(rho_full(i))
      END DO
      DO i=1,nspin
         CALL cp_fm_release(ham_single_particle(i))
      END DO
      ! Deallocate the workspace density matrix
      CALL cp_fm_release(rho_work(1))
      CALL cp_fm_release(rho_work(2))
      CALL cp_fm_release(rho_work(3))
      DEALLOCATE(ham_exp)
      DEALLOCATE(ham_effective)
      DEALLOCATE(ham_single_particle)
      DEALLOCATE(rho_full)
      DEALLOCATE(rho_work)
      ! Deallocate sigma matrices
      CALL cp_fm_release(sigma_COH)
      CALL cp_fm_release(sigma_SEX)
      ! Deallocate full overlap matrix
      CALL cp_fm_release(S_inv_fm)
      ! Deallocate bare Hartree matrix and tensor
      CALL cp_fm_release(hartree_fm)
      CALL dbt_destroy(hartree_dbt)
      DEALLOCATE(hartree_dbt)
      ! Deallocate the neighbour list that is not deallocated in gw anymore
      IF (ASSOCIATED(bs_env%nl_3c%ij_list)) CALL neighbor_list_3c_destroy(bs_env%nl_3c)
   END SUBROUTINE run_propagation_gw

! **************************************************************************************************
!> \brief Extracts the relevant information from previous run
!> \note Specifically, extracts the single particle Hamiltonian, sets up inverse overlap matrix,
!>       extracts/calculates the screening potential.
!>       So far, does not allocate any memory - forces allocation and deallocation to be done
!>       in the same place - safety feature, but maybe not necessary. 
!> \param qs_env Quickstep environment data, entry point of the calculation
!> \param ham_single_particle Pointer to array of cp_fm_type, where the single particle H will be
!> \param moments Pointer to array of dbcsr_type, where the moments are stored
! **************************************************************************************************
   SUBROUTINE tdagw_prepare(qs_env, ham_single_particle, moments_full, rho_full)
      CHARACTER(len=*), PARAMETER                        :: routineN = 'tdagw_prepare'
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: start_mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_dft,&
                                                            matrix_s,&
                                                            matrix_ks,&
                                                            moments
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: ham_single_particle,&
                                                            moments_full,&
                                                            rho_full
      INTEGER                                            :: nspin,i,re,im

      ! Get values from QS environment
      CALL get_qs_env(qs_env,&
                 bs_env=bs_env,&
                 rtp=rtp,&
                 dft_control=dft_control,&
                 mos=start_mos,&
                 matrix_s=matrix_s,&
                 matrix_ks=matrix_ks)
      ! Get the number of spins
      nspin = SIZE(start_mos)

      ! Initialize single particle Hamiltonian
      DO i=1,nspin
         IF (dft_control%rtp_control%tdagw_ham == rtp_tdagw_ham_ks) THEN
            ! Setup the KS Hamiltonian in MO basis
            CALL copy_dbcsr_to_fm(matrix_ks(i)%matrix, ham_single_particle(i))
         ELSE
            ! Setup the G0W0 Hamiltonian in MO basis
            ! TODO : Spins ...
            CALL cp_fm_to_fm(bs_env%fm_h_G0W0_Gamma, ham_single_particle(i))
            ! TODO : Redo via pointers? But fm in KS does not exist ...
            ! ham_single_particle(i)%matrix => bs_env%fm_h_G0W0_Gamma
         END IF
      END DO

      ! Initialise the moments matrix
      ! TODO : Set the moments matrix in the create_and_init_bs_env (post_scf_bandstructure_utils.F) 
      NULLIFY(moments)
      CALL dbcsr_allocate_matrix_set(moments, 3)
      ! Initialize individual matrices
      DO i=1,3
         ! Allocate the memory for the moments matrix
         ALLOCATE(moments(i)%matrix)
      END DO
      ! Initialize moments
      DO i=1,3
         ! TODO : Ideally should create by dbcsr_create, but some properties are not copied - but likely this will not be
         ! a performance bottleneck
         CALL dbcsr_copy(moments(i)%matrix, matrix_s(1)%matrix, "Moments matrix")
      END DO
      ! Create the actual moments
      CALL build_local_moment_matrix(qs_env, moments, 1)
      DO i=1,3
         CALL copy_dbcsr_to_fm(moments(i)%matrix, moments_full(i))
      END DO
      ! Once this is done, sparse moment matrices are no longer required
      ! Deallocate dbcsr moments
      CALL dbcsr_deallocate_matrix_set(moments)

      ! Initialize density matrix from mos
      CALL rt_initialize_rho_from_mos(rtp, start_mos)
      ! Create the inverse overlap matrix
      CALL s_matrices_create(matrix_s, rtp)

      ! Get the density matrix created from DFT mos
      ! Rho matrix from rtp has 2*nspin dimensions in the following order
      !  - real of spin1, imag of spin1(, real of spin2, imag of spin2)
      !  - for restricted calculation, second spins are not present
      CALL get_rtp(rtp, rho_new=rho_dft)
      ! Initialise the density matrix
      DO i=1,nspin
         re = 2*i - 1
         im = 2*i
         CALL copy_dbcsr_to_fm(rho_dft(re)%matrix, rho_full(re))
         CALL cp_fm_set_all(rho_full(im), 0.0_dp)
      END DO

   END SUBROUTINE tdagw_prepare

! **************************************************************************************************
!> \brief Calculates the self-energy by contraction of screened potential
!> \note Can be used for both the Coulomb hole part and screened exchange part 
!> \param qs_env Quickstep environment data, entry point of the calculation
!> \param greens_fm Pointer to the Green's function matrix, which is used as input data
!> \param sigma_fm Pointer to the self-energy full matrix, which is overwritten by this routine
!> \author Stepan Marek
!> \date 01.2024
! **************************************************************************************************
   SUBROUTINE get_sigma(qs_env, greens_fm, sigma_fm, prefactor_opt)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), POINTER, INTENT(IN)              :: greens_fm ! matrix to contract with RI_W
      TYPE(cp_fm_type), POINTER, INTENT(OUT)             :: sigma_fm ! resulting self energy
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: prefactor_opt
      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_sigma'
      REAL(kind=dp)                                      :: prefactor
      TYPE(dbt_type)                                     :: t_3c_w, t_3c_g, t_work1, t_work2,&
                                                            w_dbt, greens_dbt, sigma_dbt
      TYPE(dbcsr_type)                                   :: w_dbcsr_mat, greens_dbcsr_mat, sigma_dbcsr
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      IF(PRESENT(prefactor_opt)) THEN
         prefactor = prefactor_opt
      ELSE
         prefactor = 1.0_dp
      END IF

      CALL get_qs_env(qs_env,&
                      bs_env=bs_env)

      ! Three-centre integrals are obtained from build_3c_integrals, from qs_tensors
      ! These should use sparcity, while W and Sigma can be full matrices
      ! The summation is carried out by dbt library - dbt_contract in dbt_api
      ! The building of the tensors might be a bit hard, because it requires a lot of parallel information
      ! Probably just use the tensors already present in bs_env? They seem to be mostly work tensors
      ! TODO : Starting with unbounded integrals, then probably need to reduce to a subset
      ! Create by template
      CALL dbt_create(bs_env%t_RI__AO_AO, t_3c_w)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_work1)
      CALL dbt_create(bs_env%t_W, w_dbt)
      CALL compute_3c_integrals(qs_env,&
                         bs_env,&
                         t_3c_w)
      ! Copy W to dbcsr and then to dbt
      CALL dbcsr_create(w_dbcsr_mat, name="w", template=bs_env%mat_RI_RI%matrix)
      CALL copy_fm_to_dbcsr(bs_env%fm_W_MIC_freq_zero, w_dbcsr_mat)
      ! TODO : Maybe call filter here?
      CALL dbt_copy_matrix_to_tensor(w_dbcsr_mat, w_dbt)
      CALL dbt_contract(alpha=1.0_dp,&
                   tensor_1=w_dbt,&
                   tensor_2=t_3c_w,&
                   beta=0.0_dp,&
                   tensor_3=t_work1,&
                   contract_1=[2],notcontract_1=[1], map_1=[1],&
                   contract_2=[1],notcontract_2=[2,3], map_2=[2,3],&
                   filter_eps=bs_env%eps_filter)
      ! t_work1 now contains sum _ P (alpha beta | P) W _ (PQ) (iomega = 0)
      ! Next step is to convert the greens full matrix to dbcsr matrix 
      CALL dbcsr_create(greens_dbcsr_mat, name="greens", template=bs_env%mat_ao_ao%matrix)
      CALL copy_fm_to_dbcsr(greens_fm, greens_dbcsr_mat)
      CALL dbt_create(bs_env%t_G, greens_dbt)
      CALL dbt_copy_matrix_to_tensor(greens_dbcsr_mat, greens_dbt) 
      ! TODO : Is this the correct way to compute different direction of the ingetrals?
      ! TODO : Do we need different t_3c? Maybe for performance reasons ...
      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c_g)
      CALL compute_3c_integrals(qs_env,&
                         bs_env,&
                         t_3c_g)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_work2)
      ! Then contract it
      ! no scaling applied - this has to be applied externally TODO : Is this a good approach?
      CALL dbt_contract(alpha=1.0_dp,&
                   tensor_1=greens_dbt,&
                   tensor_2=t_3c_g,&
                   beta=0.0_dp,&
                   tensor_3=t_work2,&
                   contract_1=[1], notcontract_1=[2], map_1=[2],&
                   contract_2=[3], notcontract_2=[1,2], map_2=[1,3],&
                   filter_eps=bs_env%eps_filter)
      ! Now, final contraction over remaining indices
      CALL dbt_create(bs_env%t_G, sigma_dbt)
      ! Final contraction, which introduces the given prefactor
      CALL dbt_contract(alpha=prefactor,&
                   tensor_1=t_work1,&
                   tensor_2=t_work2,&
                   beta=0.0_dp,&
                   tensor_3=sigma_dbt,&
                   contract_1=[1,2], notcontract_1=[3], map_1=[2],&
                   contract_2=[1,2], notcontract_2=[3], map_2=[1],&
                   filter_eps=bs_env%eps_filter)
      ! Finally, convert the COH tensor to matrix and then to fm matrix
      CALL dbcsr_create(sigma_dbcsr, name="sigma", template=bs_env%mat_ao_ao%matrix)
      CALL dbt_copy_tensor_to_matrix(sigma_dbt, sigma_dbcsr)
      CALL copy_dbcsr_to_fm(sigma_dbcsr, sigma_fm)
      ! TODO : Deallocations
      CALL dbt_destroy(t_3c_w)
      CALL dbt_destroy(t_work1)
      CALL dbt_destroy(w_dbt)
      CALL dbt_destroy(greens_dbt)
      CALL dbt_destroy(t_3c_g)
      CALL dbt_destroy(t_work2)
      CALL dbt_destroy(sigma_dbt)
      CALL dbcsr_release(sigma_dbcsr)
      CALL dbcsr_release(greens_dbcsr_mat)
      CALL dbcsr_release(w_dbcsr_mat)
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
      TYPE(cp_fm_type), DIMENSION(:,:), ALLOCATABLE      :: metric_fm,&
                                                            metric_fm_2
      TYPE(cp_fm_type)                                   :: metric_inv_fm,&
                                                            work_fm
      TYPE(dbcsr_type)                                   :: V_dbcsr

      CALL get_qs_env(qs_env, bs_env=bs_env)

      ! Allocate for bare Hartree term
      CALL cp_fm_create(V_fm, bs_env%fm_W_MIC_freq_zero%matrix_struct)

      ! Calculate full coulomb RI basis elements - V _ (PQ) matrix
      CALL RI_2c_integral_mat(qs_env, metric_fm, V_fm, &
                              bs_env%n_RI, coulomb_op, do_kpoints=.FALSE., &
                              regularization_RI=bs_env%regularization_RI)
      ! Actual data are in the metric_inv_fm matrix
      ! TODO : Spin and k-points
      CALL cp_fm_to_fm(metric_fm(1,1), V_fm)
      ! Since the routine below allocates the matrices, they have to be freed again
      CALL cp_fm_release(metric_fm(1,1))
      DEALLOCATE(metric_fm)
      ! Now, calculate the matrices in the standard RI metric
      CALL RI_2c_integral_mat(qs_env, metric_fm_2, bs_env%fm_RI_RI,&
                              bs_env%n_RI, bs_env%ri_metric, do_kpoints=.FALSE.,&
                              regularization_RI=bs_env%regularization_RI)
      ! These matrices have yet to be inverted
      CALL cp_fm_create(metric_inv_fm, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(work_fm, bs_env%fm_RI_RI%matrix_struct)
      ! TODO : Spins and k-points
      CALL cp_fm_invert(metric_fm_2(1,1), metric_inv_fm)
      ! Multiply by the inverse from each site (M^-1 is symmetric)
      CALL parallel_gemm("N", "N", bs_env%n_RI, bs_env%n_RI, bs_env%n_RI,&
                         1.0_dp, metric_inv_fm, V_fm, 0.0_dp, work_fm)
      CALL parallel_gemm("N", "N", bs_env%n_RI, bs_env%n_RI, bs_env%n_RI,&
                         1.0_dp, work_fm, metric_inv_fm, 0.0_dp, V_fm)
      ! Now, create the tensor from the matrix
      ! First, convert full matrix to dbcsr
      CALL dbcsr_create(V_dbcsr, name="Hartree_dbcsr", template=bs_env%mat_RI_RI%matrix)
      ! TODO : No k-points so far
      CALL copy_fm_to_dbcsr(V_fm, V_dbcsr)
      ! Create and copy distinctly, so that unnecessary objects can be destroyed
      CALL dbt_create(bs_env%t_W, V_dbt, name="Hartree_dbt")
      CALL dbt_copy_matrix_to_tensor(V_dbcsr, V_dbt)
      ! Destroy all unnecessary matrices
      CALL dbcsr_release(V_dbcsr)
      ! TODO : No k-points so far
      CALL cp_fm_release(V_fm)
      CALL cp_fm_release(metric_fm_2(1,1))
      DEALLOCATE(metric_fm_2)
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
      LOGICAL                                            :: data_transpose
      INTEGER, DIMENSION(2)                              :: block_shape
      ! Pair of indices for current block
      INTEGER, DIMENSION(3)                              :: index_tensor
      REAL(kind=dp), DIMENSION(:,:), POINTER             :: block_matrix
      REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE       :: block_tensor
      INTEGER                                            :: j,number_of_blocks,backup
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: indices_block_matrix_1,&
                                                            indices_block_matrix_2,&
                                                            indices_block_matrix_custom
      LOGICAL                                            :: found
      
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
      ALLOCATE(indices_block_matrix_custom(number_of_blocks))
      CALL dbcsr_iterator_start(iterator_matrix, matrix_desym_p)
      DO j=1,number_of_blocks
         CALL dbcsr_iterator_next_block(iterator_matrix, indices_block_matrix_1(j), indices_block_matrix_2(j), block_matrix)
         indices_block_matrix_custom(j) = 1
      END DO
      CALL dbcsr_iterator_stop(iterator_matrix)
      ! Construct parallelisation bounds
      PRINT *, indices_block_matrix_1
      PRINT *, indices_block_matrix_2
      PRINT *, indices_block_matrix_custom
      ! TODO : OMP DO does not register as OMP parallel - dbm requires reservations to be called within parallel block
      ! Proceed as in dbt - determine distribution of blocks among threads - maybe distvec function? Or integer division
      ! TODO : The reservations are a confusing mess of a spaghetti code - tread lightly
      !        Sure, the statement below does nthreads-times more work than necessary, but it works - any attempts at
      !        parallelisations have failed
      backup = omp_get_max_threads()
      CALL omp_set_num_threads(1)
      !$OMP PARALLEL DEFAULT(none) SHARED(tensor, indices_block_matrix_1, indices_block_matrix_2, indices_block_matrix_custom)
      ! All index fields should be identical 
      CALL dbt_reserve_blocks(tensor, indices_block_matrix_1, indices_block_matrix_2, indices_block_matrix_custom)
      !$OMP END PARALLEL
      CALL omp_set_num_threads(backup)
      PRINT *, "Blocks reserved"
      DO j=1,number_of_blocks
         CALL dbcsr_get_block_p(matrix_desym_p, indices_block_matrix_1(j), indices_block_matrix_2(j),&
                             block_matrix, data_transpose,found, block_shape(1), block_shape(2))
         ALLOCATE(block_tensor(block_shape(1),block_shape(2),1))
         block_tensor(:,:,1) = block_matrix(:,:)
         index_tensor = [indices_block_matrix_1(j), indices_block_matrix_2(j),1]
         CALL dbt_put_block(tensor, index_tensor, SHAPE(block_tensor), block_tensor)
         DEALLOCATE(block_tensor)
      END DO
      PRINT *, "Blocks copied"
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
!> \param qs_env Quickstep environment - entry point of calculation
!> \param V_dbt Tensor representation of the Hartree matrix in the auxiliary basis
!> \param rho Density matrix in ao basis
!> \param V_ao Overwritten by the Hartree matrix in the atomic orbital basis
!> \author Stepan Marek
!> \date 01.2024
! **************************************************************************************************
   SUBROUTINE get_hartree(qs_env, V_dbt, rho, V_ao)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbt_type), INTENT(INOUT)                      :: V_dbt
      TYPE(cp_fm_type), INTENT(IN)                       :: rho
      TYPE(cp_fm_type), POINTER, INTENT(OUT)             :: V_ao
      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_hartree'
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbt_type)                                     :: t_3c_w,&
                                                            t_work1,&
      !                                                     rho_dbt,&
                                                            t_one_plus_one
      TYPE(dbt_type), ALLOCATABLE, TARGET                :: t_two_plus_one
      TYPE(dbt_type), POINTER                            :: t_two_plus_one_p
      INTEGER, ALLOCATABLE, DIMENSION(:,:)               :: distvec
      INTEGER, DIMENSION(1)                              :: distvec_one
      TYPE(dbt_distribution_type)                        :: dist
      TYPE(dbcsr_type), ALLOCATABLE, TARGET              :: rho_dbcsr
      TYPE(dbcsr_type), POINTER                          :: rho_dbcsr_p
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pdims, pcoord
      TYPE(dbt_pgrid_type)                               :: pgrid
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: empty_list
      INTEGER                                            :: j
      
      CALL get_qs_env(qs_env, bs_env=bs_env)

      ! Create the empty list, used for no contraction indices
      ALLOCATE(empty_list(0))

      CALL dbt_create(bs_env%t_RI__AO_AO, t_3c_w)
      CALL dbt_create(bs_env%t_RI__AO_AO, t_work1)
      ! build 3c integrals
      CALL compute_3c_integrals(qs_env,&
                         bs_env,&
                         t_3c_w)
      ! Create a 1D tensor (vector) for the result of contraction with density matrix
      ! TODO : 1D tensor is not available for construction, but maybe can be emulated by
      !        2D tensor with one dimension set to 1. However, extracting 2D matrix from
      !        3D tensor with one dimension set to 1 is not implemented
      ! Taken from gw_utils

      ! Create the 2-index tensor with one dimension of size 1 from tensor primitives
      ! Create 2-index process grid
      ALLOCATE(pdims(3))
      ALLOCATE(pcoord(3))
      ! Set elements to zero - otherwise can get erronous behaviour
      DO j=1,3
         pdims(j) = 0
         pcoord(j) = 0
      END DO
      ! Process grid creation
      CALL dbt_pgrid_create(bs_env%para_env_tensor, pdims, pgrid)
      ! Determine processes on the created grid
      CALL dbt_mp_environ_pgrid(pgrid, pdims, pcoord)
      ! Prepare task distribution on the process grid for each dimension
      ALLOCATE(distvec(2, SIZE(bs_env%sizes_AO)))
      DO j=1,2
         CALL dbt_default_distvec(SIZE(bs_env%sizes_AO), pdims(j), bs_env%sizes_AO, distvec(j,:))
      END DO
      CALL dbt_default_distvec(1, pdims(3), [1], distvec_one)
      ! TODO : Maybe just [1] instead of pdims(2)?
      ! Create the tensor distribution
      CALL dbt_distribution_new(dist, pgrid, distvec(1,:), distvec(2,:), distvec_one)
      ! Create the tensor aligned for the task distribution
      ALLOCATE(t_two_plus_one)
      CALL dbt_create(t_two_plus_one, "four index tensor", dist, [1], [2, 3],&
                      bs_env%sizes_AO,&
                      bs_env%sizes_AO,&
                      [1])
      ! No longer need the process information
      CALL dbt_pgrid_destroy(pgrid)
      CALL dbt_distribution_destroy(dist)
      DEALLOCATE(distvec)
      DEALLOCATE(pdims)
      DEALLOCATE(pcoord)

      ! Create the one + one tensor - true vector
      ALLOCATE(pdims(2))
      ALLOCATE(pcoord(2))
      DO j=1,2
         pdims(j) = 0
         pcoord(j) = 0
      END DO
      ALLOCATE(distvec(1,SIZE(bs_env%sizes_RI)))
      ! Create the process grid
      CALL dbt_pgrid_create(bs_env%para_env_tensor, pdims, pgrid)
      ! Populate the process grid
      CALL dbt_mp_environ_pgrid(pgrid, pdims, pcoord)
      ! Prepare distribution vectors that allocate tasks to tensor blocks
      CALL dbt_default_distvec(SIZE(bs_env%sizes_RI), pdims(1), bs_env%sizes_RI, distvec(1,:))
      CALL dbt_default_distvec(1, pdims(2), [1], distvec_one)
      ! Create tensor distribution
      CALL dbt_distribution_new(dist, pgrid, distvec(1,:), distvec_one)
      CALL dbt_create(t_one_plus_one, "tensor vector", dist, [1], [2], bs_env%sizes_RI, [1])
      ! No longer need the process information
      CALL dbt_pgrid_destroy(pgrid)
      CALL dbt_distribution_destroy(dist)
      DEALLOCATE(distvec)
      DEALLOCATE(pcoord)
      DEALLOCATE(pdims)

      ! Convert density to tensor
      ALLOCATE(rho_dbcsr)
      CALL dbcsr_create(rho_dbcsr, name="Sparse density", template=bs_env%mat_ao_ao%matrix)
      CALL copy_fm_to_dbcsr(rho, rho_dbcsr)
      ! Maybe doable via creation
      ! CALL dbt_create(bs_env%t_G, rho_dbt, "Density tensor")
      ! CALL dbt_copy_matrix_to_tensor(rho_dbcsr, t_two_plus_one)
      rho_dbcsr_p => rho_dbcsr
      t_two_plus_one_p => t_two_plus_one
      CALL copy_matrix_to_3d_tensor(rho_dbcsr_p, t_two_plus_one_p)
      PRINT *, "Matrix copied to tensor"
      ! Contract T _ (Q 1) = sum _ Q (Q|mu'nu') D _ (mu' nu' 1)
      CALL dbt_contract(alpha=1.0_dp,&
                   tensor_1=t_3c_w,&
                   tensor_2=t_two_plus_one,&
                   beta=0.0_dp,&
                   tensor_3=t_one_plus_one,&
                   contract_1=[2,3],notcontract_1=[1],map_1=[1],&
                   contract_2=[1,2],notcontract_2=[3],map_2=[2],&
                   filter_eps=bs_env%eps_filter)
      ! Contract Q _ (mu nu Q) = sum _ P (mu nu | P) V _ (PQ)
      CALL dbt_contract(alpha=1.0_dp,&
                   tensor_1=t_3c_w,&
                   tensor_2=V_dbt,&
                   beta=0.0_dp,&
                   tensor_3=t_work1,&
                   contract_1=[1],notcontract_1=[2,3],map_1=[2,3],&
                   contract_2=[1],notcontract_2=[2],map_2=[1],&
                   filter_eps=bs_env%eps_filter)
      ! Final contraction V ^ H _ (mu nu 1) = sum _ (Q) Q _ (mu nu Q) T _ (Q 1)
      CALL dbt_contract(alpha=1.0_dp,&
                   tensor_1=t_work1,&
                   tensor_2=t_one_plus_one,&
                   beta=0.0_dp,&
                   tensor_3=t_two_plus_one,&
                   contract_1=[1], notcontract_1=[2,3], map_1=[1,2],&
                   contract_2=[1], notcontract_2=[2], map_2=[3],&
                   filter_eps=bs_env%eps_filter)
      ! TODO : Allocations, conversion from t_two_plus_one back to matrix - check that the copying makes sense
      ! Copy result to output matrix
      ! Use rho_dbcsr - TODO : Instead, use some workspace?
      CALL dbcsr_clear(rho_dbcsr)
      CALL copy_3d_tensor_to_matrix(t_two_plus_one_p, rho_dbcsr_p)
      PRINT *, "Tensor copied to matrix"
      CALL copy_dbcsr_to_fm(rho_dbcsr, V_ao)
      ! TODO : Deallocations, initialisations, CALL statements etc.
      CALL dbt_destroy(t_3c_w)
      CALL dbt_destroy(t_work1)
      ! CALL dbt_destroy(rho_dbt)
      CALL dbt_destroy(t_two_plus_one)
      CALL dbt_destroy(t_one_plus_one)
      CALL dbcsr_release(rho_dbcsr)
      DEALLOCATE(rho_dbcsr)
      DEALLOCATE(t_two_plus_one)
      DEALLOCATE(empty_list)
   END SUBROUTINE
END MODULE rt_tdagw

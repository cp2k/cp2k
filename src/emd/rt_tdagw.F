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
                                              dbcsr_create,&
                                              dbcsr_release,&
                                              dbcsr_copy,&
                                              dbcsr_set
   USE dbt_api,                         ONLY: dbt_create,&
                                              dbt_contract,&
                                              dbt_copy_matrix_to_tensor,&
                                              dbt_copy_tensor_to_matrix,&
                                              dbt_destroy,&
                                              dbt_type
   USE qs_tensors,                      ONLY: neighbor_list_3c_destroy
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE efield_utils,                    ONLY: make_field
   USE rt_propagator_init,              ONLY: rt_initialize_rho_from_mos
   USE rt_propagation_methods,          ONLY: s_matrices_create
   USE gw_methods,                      ONLY: compute_3c_integrals
   USE matrix_exp,                      ONLY: taylor_only_imaginary
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
      CHARACTER(len=*), PARAMETER                        :: routineN = 'run_propagation_gw'
      TYPE(qs_environment_type), POINTER                 :: qs_env
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
                                                            S_inv_fm
      TYPE(cp_fm_type), POINTER                          :: sigma_p,&
                                                            greens_fm_p
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
      PRINT *, "== Sigma COH =="
      CALL cp_fm_write_formatted(sigma_COH, 6)

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
            CALL taylor_only_imaginary(ham_exp(re:im), ham_effective(re), rtp%orders(1,j), rtp%orders(2,j))
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
                   beta=1.0_dp,&
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
                   beta=1.0_dp,&
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
                   beta=1.0_dp,&
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
END MODULE rt_tdagw

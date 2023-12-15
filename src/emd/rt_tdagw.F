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
                                              dbcsr_copy,&
                                              dbcsr_set
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
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
      ! TODO : Move to separate file
      CHARACTER(len=*), PARAMETER                        :: routineN = 'run_propagation_gw'
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input, rtp_section
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: start_mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_dft,&
                                                            matrix_s,&
                                                            matrix_ks
      TYPE(dbcsr_type), POINTER                          :: S_inv
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: moments
      TYPE(cp_fm_type), DIMENSION(:), ALLOCATABLE        :: moments_full
      TYPE(cp_fm_type), DIMENSION(:), ALLOCATABLE        :: ham_effective,&
                                                            ham_exp,&
                                                            ham_single_particle,&
                                                            rho_full,&
                                                            rho_work
      REAL(kind=dp)                                      :: field(3)
      INTEGER                                            :: i,j,k,nspin,re,im
      INTEGER                                            :: rho_unit

      ! Marek : Start by getting the DFT results - density matrix
      !         This matrix is identical with the equilibrium time-diagonal lesser Green's function in G0W0
      CALL get_qs_env(qs_env,&
                      rtp=rtp,&
                      bs_env=bs_env,&
                      mos=start_mos,&
                      dft_control=dft_control,&
                      input=input,&
                      matrix_s=matrix_s,&
                      matrix_ks=matrix_ks)
      ! TODO : Ask Jan whether there is already an identifier present in some of the environments
      nspin = SIZE(start_mos)
      CALL rt_initialize_rho_from_mos(rtp, start_mos)
      ! Create the inverse overlap matrix
      CALL s_matrices_create(matrix_s, rtp)
      CALL get_rtp(rtp, rho_new=rho_dft, S_inv=S_inv)
      ! Rho matrix from rtp has 2*nspin dimensions in the following order
      !  - real of spin1, imag of spin1(, real of spin2, imag of spin2)
      !  - for restricted calculation, second spins are not present

      ! Marek : Set up the single particle Hamiltonian - can be either KS or G0W0 Hamiltonian
      ! TODO : Always uses bs_env - truly only for debug or redo with possibility to refer to qs_env?
      ALLOCATE(ham_single_particle(nspin))
      DO i=1,nspin
         CALL cp_fm_create(ham_single_particle(i), bs_env%fm_h_G0W0_Gamma%matrix_struct)
         IF (dft_control%rtp_control%tdagw_ham == rtp_tdagw_ham_ks) THEN
            ! Setup the KS Hamiltonian in MO basis
            CALL copy_dbcsr_to_fm(matrix_ks(i)%matrix, ham_single_particle(i))
         ELSE
            ! Setup the G0W0 Hamiltonian in MO basis
            ! TODO : Redo via pointers, no need to copy this?
            ! TODO : Spins ...
            CALL cp_fm_to_fm(bs_env%fm_h_G0W0_Gamma, ham_single_particle(i))
         END IF
      END DO
      ! Convert the density matrix into full matrix TODO : Spin coupling?
      ALLOCATE(rho_full(2*nspin))
      ALLOCATE(rho_work(3))
      CALL cp_fm_create(rho_work(1), bs_env%fm_h_G0W0_Gamma%matrix_struct)
      CALL cp_fm_create(rho_work(2), bs_env%fm_h_G0W0_Gamma%matrix_struct)
      CALL cp_fm_create(rho_work(3), bs_env%fm_h_G0W0_Gamma%matrix_struct)
      CALL cp_fm_set_all(rho_work(1), 0.0_dp)
      CALL cp_fm_set_all(rho_work(2), 0.0_dp)
      CALL cp_fm_set_all(rho_work(3), 0.0_dp)
      DO i=1,2*nspin
         CALL cp_fm_create(rho_full(i), bs_env%fm_h_G0W0_Gamma%matrix_struct)
      END DO
      DO i=1,nspin
         re = 2*i - 1
         im = 2*i
         CALL copy_dbcsr_to_fm(rho_dft(re)%matrix, rho_full(re))
         CALL cp_fm_set_all(rho_full(im), 0.0_dp)
      END DO

      ! Marek : Now, get the G0W0 results - Hamiltonian + screened interaction matrices
      ! TODO : Ask Jan about correct logger and access to verbosity levels?
      PRINT *, "Writing the full ao G0W0 Hamiltonian"
      CALL cp_fm_write_formatted(bs_env%fm_h_G0W0_Gamma, 6)
      ! Marek : After getting the G0W0 Hamiltonian, first check the no-self-energy time propagation
      ! In order to do that, need to be able to add the electric field to the Hamiltonian at every simulation step
      ! Electric field is added to the core hamiltonian in efield_potential_length_gauge
      ! Here, we reimplement it in order to add it directly to the effective Hamiltonian at each timestep
      ! TODO : Set the moments matrix in the create_and_init_bs_env (post_scf_bandstructure_utils.F) 
      NULLIFY(moments)
      CALL dbcsr_allocate_matrix_set(moments, 3)
      ! Initialize individual matrices
      DO i=1,3
         ! Allocate the memory for the moments matrix
         ALLOCATE(moments(i)%matrix)
         ! Create by template
         CALL dbcsr_copy(moments(i)%matrix, matrix_s(1)%matrix, "Moments matrix")
         ! TODO : Ideally should create by dbcsr_create, but some properties are not copied - but likely this will not be
         ! a performance bottleneck
         ! Fill by zeros - TODO : is this necessary?
         CALL dbcsr_set(moments(i)%matrix, 0.0_dp)
      END DO
      ! Create the actual moments
      CALL build_local_moment_matrix(qs_env, moments, 1)
      ! Transform the moments into full matrices
      ! Start by allocating the array
      ALLOCATE(moments_full(3))
      ! Now, for each index, create the matrix and fill it with dbcsr values
      DO i=1,3
         ! Use cp_fm_struct from Kohn-Sham matrix TODO : Differentiate by spin? - not necessary for moments 
         CALL cp_fm_create(moments_full(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL copy_dbcsr_to_fm(moments(i)%matrix, moments_full(i))
      END DO
      ! Once this is done, sparse moment matrices are no longer required
      ! Deallocate dbcsr moments
      CALL dbcsr_deallocate_matrix_set(moments)
      ! Now the full matrices are populated, try to print some
      PRINT *, "Writing the moments matrices"
      CALL cp_fm_write_formatted(moments_full(1), 6)
      CALL cp_fm_write_formatted(moments_full(2), 6)
      CALL cp_fm_write_formatted(moments_full(3), 6)
      ! storage (available in imaginary time and frequency as fm_W_MIC_freq and fm_W_MIC_time_tensor, respectively)
      ! Create the storage for effective Hamiltonian
      ! Create storage for the propagator matrix - Hamiltonian exponential
      ALLOCATE(ham_effective(2*nspin))
      ALLOCATE(ham_exp(2*nspin))
      DO i=1,2*nspin
         CALL cp_fm_create(ham_effective(i), bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_create(ham_exp(i), bs_env%fm_h_G0W0_Gamma%matrix_struct)
      END DO
      ! Setup the files
      logger => cp_get_default_logger()
      rtp_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION")
      rho_unit = cp_print_key_unit_nr(logger, rtp_section, "PRINT%DENSITY_MATRIX", extension=".dat")
      ! ********************** Start the time loop **********************
      DO i = rtp%i_start, rtp%nsteps
         ! Reset the effective Hamiltonian to KS Hamiltonian + G0W0
         ! TODO : Spin dependence in the KS + G0W0 Hamiltonian
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
         ! TODO : Add the field to the Hamiltonian via the constant moments matrix - might not be constant in 2D bulk via k-points?
         ! The field is added to the real part of both spins
         ! TODO : Preprocessor macro that sets up the spin loop
         DO j=1,nspin
            re = 2*j - 1
            im = 2*j
            DO k=1,3
               CALL cp_fm_scale_and_add(1.0_dp, ham_effective(re), field(k), moments_full(k))
            END DO
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
         ! Now, rho_full should contain density matrix at time t+dt
         ! TODO : Print the density matrix
         ! Marek : In each timestep of the loop
         !          - update the self-energy
         !          - propagate G< via matrix exponential
         !          - if it is required, print out G< or its part
      END DO
      ! ********************** End the time loop **********************
      ! Deallocate fm moments
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
   END SUBROUTINE run_propagation_gw
END MODULE rt_tdagw

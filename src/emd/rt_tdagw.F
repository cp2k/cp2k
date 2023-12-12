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
                                              cp_fm_write_formatted
   USE kinds,                           ONLY: dp
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_print,&
                                              dbcsr_create,&
                                              dbcsr_copy,&
                                              dbcsr_set
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE efield_utils,                    ONLY: make_field
   USE rt_propagator_init,              ONLY: rt_initialize_rho_from_mos

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
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: start_mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_dft,&
                                                            matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: moments
      TYPE(cp_fm_type), DIMENSION(:), ALLOCATABLE        :: moments_full
      TYPE(cp_fm_type), DIMENSION(:), ALLOCATABLE        :: ham_effective,&
                                                            ham_exp,&
                                                            rho_full,&
                                                            rho_work
      REAL(kind=dp)                                      :: field(3)
      INTEGER                                            :: i,nspin,re,im

      ! Marek : Start by getting the DFT results - density matrix
      !         This matrix is identical with the equilibrium time-diagonal lesser Green's function in G0W0
      CALL get_qs_env(qs_env,&
                      rtp=rtp,&
                      bs_env=bs_env,&
                      mos=start_mos,&
                      dft_control=dft_control,&
                      matrix_s=matrix_s)
      ! TODO : Ask Jan whether there is already an identifier present in some of the environments
      nspin = SIZE(start_mos)
      CALL rt_initialize_rho_from_mos(rtp, start_mos)
      CALL get_rtp(rtp, rho_new=rho_dft)
      ! Rho matrix from rtp has 2*nspin dimensions in the following order
      !  - real of spin1, imag of spin1(, real of spin2, imag of spin2)
      !  - for restricted calculation, second spins are not present

      ! CALL dbcsr_print(rho_dft(1)%matrix)

      ! Convert the density matrix into full matrix TODO : Spin coupling?
      ALLOCATE(rho_full(2*nspin))
      ALLOCATE(rho_work(2*nspin))
      DO i=1,2*nspin
         CALL cp_fm_create(rho_full(i), bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_create(rho_work(i), bs_env%fm_h_G0W0_Gamma%matrix_struct)
      END DO
      DO i=1,nspin
         re = 2*i - 1
         im = 2*i
         CALL copy_dbcsr_to_fm(rho_dft(re)%matrix, rho_full(re))
         CALL cp_fm_set_all(rho_full(im), 0.0_dp)
         ! TODO : Initialisation of rho_work might not be necessary
         CALL cp_fm_set_all(rho_work(re), 0.0_dp)
         CALL cp_fm_set_all(rho_work(im), 0.0_dp)
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
      ! Get the field vector
      CALL make_field(dft_control, field, qs_env%sim_step, qs_env%sim_time)
      ! Print the field
      PRINT *, "Field : ", field(1), field(2), field(3)
      ! Add the field to the effective Hamiltonian
      ! Now the full matrices are populated, try to print some
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
      ! Marek : Setup the time loop
      ! DO i = rtp%i_start, rtp%nsteps
      !    ! Reset the effective Hamiltonian to KS Hamiltonian + G0W0
      !    cp_fm_to_fm(bs_env%fm_h_G0W0_Gamma, ham_effective)
      !    ! Update the field
      !    qs_env%sim_time = REAL(i, dp)*rtp%dt
      !    qs_env%sim_step = i
      !    CALL make_field(dft_control, field, qs_env%sim_step, qs_env%sim_time)
      !    ! TODO : Add the field to the Hamiltonian via the constant moments matrix - might not be constant in bulk
      !    DO j=1,3
      !       CALL cp_fm_scale_and_add(1.0_dp, ham_effective, field(j), moments_full(j))
      !    END DO
      !    ! TODO : Create the matrix exponential (matrix_exponential.F:taylor_only_imaginary, but might need complex because of
      !    ! complex W and hence complex self-energy?)
      !    ! TODO : Understand which input options set the orders
      !    ! TODO : Scale the Hamiltonian so that the timestep dt is present in the exponent (with the correct sign)
      !    ! TODO : Set up the ham_exp properly - two components - real and imaginary
      !    CALL taylor_only_imaginary(ham_exp, ham_effective, rtp%orders(1,1), rtp%orders(2,1))
      !    ! TODO : Multiply the density matrix by the propagator from both sides - look at rt_propagation_step for details
      !    ! Use : cp_fm_gemm
      !    CALL cp_fm_gemm("N", "N", bs_env%n_ao, bs_env%n_ao, bs_env%n_ao, 1.0_dp, ham_exp, rho_full, 0.0_dp, rho_work)
      !    CALL cp_fm_gemm("N", "H", bs_env%n_ao, bs_env%n_ao, bs_env%n_ao, 1.0_dp, rho_work, ham_exp, 0.0_dp, rho_full)
      !    ! Now, rho_full should contain density matrix at time t+dt
      !    ! TODO : Print the density matrix
      !    ! Marek : In each timestep of the loop
      !    !          - update the field
      !    !          - update the self-energy
      !    !          - propagate G< via matrix exponential
      !    !          - if it is required, print out G< or its part
      ! END DO
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
         ! Deallocate the workspace density matrix
         CALL cp_fm_release(rho_work(i))
      END DO
      DEALLOCATE(ham_exp)
      DEALLOCATE(ham_effective)
      DEALLOCATE(rho_full)
      DEALLOCATE(rho_work)
   END SUBROUTINE run_propagation_gw
END MODULE rt_tdagw

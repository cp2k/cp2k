!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Data storage and other types for time propagation via TD-aGW method.
!> \author Stepan Marek (01.24)
! **************************************************************************************************

MODULE rt_tdagw_types

   USE kinds,                           ONLY: dp
   USE cp_fm_types,                     ONLY: cp_fm_type,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_create,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm
   USE dbcsr_api,                       ONLY: dbcsr_type,&
                                              dbcsr_p_type,&
                                              dbcsr_create,&
                                              dbcsr_release,&
                                              dbcsr_copy
   USE dbt_api,                         ONLY: dbt_type,&
                                              dbt_destroy
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE qs_mo_types,                     ONLY: mo_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE qs_environment_types,            ONLY: qs_environment_type,&
                                              get_qs_env
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE rt_propagation_types,            ONLY: rt_prop_type,&
                                              get_rtp
   USE rt_propagator_init,              ONLY: rt_initialize_rho_from_mos
   USE rt_propagation_methods,          ONLY: s_matrices_create
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE matrix_exp,                      ONLY: get_nsquare_norder

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_tdagw"

   #:include "tdagw_macros.fypp"
   
   PUBLIC :: tdagw_env_type,&
             create_tdagw_env,&
             release_tdagw_env

! **************************************************************************************************
!> \param nspin Number of spin channels that are present
!> \param field Electric field calculated at the given timestep
!> \param moments Moment operators along cartesian directions - used to coupling to the field and
!>                for derivation of expectation values of moments, which relates to emission
!> \param sim_step Current step of the simulation
!> \param sim_time Current time of the simulation
!> \param sim_dt Timestep of the simulation
!> \param ham_effective Real and imaginary part of the effective Hamiltonian used to propagate
!>                      the density matrix
!> \param ham_single_particle Single particle Hamiltonian, taken from G0W0/DFT TODO : Only real?
!> \param ham_workspace Workspace matrices for use with the Hamiltonian propagation - storage of
!>                      exponential propagators etc.
!> \param rho Density matrix at the current time step
!> \param rho_workspace Matrices for storage of density matrix at different timesteps for
!>                      interpolation and self-consistency checks etc.
!> \param S_inv Inverse overlap matrix, used in propagation of the density and in COH contribution
!>              calculation
!> \author Stepan Marek (01.24)
! **************************************************************************************************
   TYPE tdagw_env_type
      INTEGER                                                   :: n_spin,&
                                                                   n_ao,&
                                                                   n_RI
      REAL(kind=dp), DIMENSION(3)                               :: field
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: moments
      INTEGER                                                   :: sim_step,&
                                                                   sim_start,&
                                                                   sim_nsteps
      REAL(kind=dp)                                             :: sim_time,&
                                                                   sim_dt,&
                                                                   etrs_threshold
      TYPE(dft_control_type), POINTER                           :: dft_control

      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: ham_effective
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: ham_single_particle
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: ham_workspace

      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: rho,&
                                                                   rho_new,&
                                                                   rho_new_last,&
                                                                   rho_M
      ! TODO : Check how much is actually needed and whether it is not better to separate names
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: rho_workspace
      TYPE(dbcsr_type), POINTER                                 :: S_inv
      TYPE(dbt_type), POINTER                                   :: hartree_dbt
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: sigma_SEX, sigma_COH
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: hartree_orig, hartree_curr
      INTEGER, DIMENSION(:,:), POINTER                          :: orders

   END TYPE tdagw_env_type

CONTAINS

   SUBROUTINE create_tdagw_env(tdagw_env, qs_env)
      TYPE(tdagw_env_type)                                      :: tdagw_env
      TYPE(qs_environment_type),POINTER                         :: qs_env
      CHARACTER(len=*), PARAMETER                               :: routineN="create_tdagw_env"
      TYPE(post_scf_bandstructure_type), POINTER                :: bs_env
      TYPE(rt_prop_type), POINTER                               :: rtp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER                 :: moments_dbcsr_p
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER                 :: matrix_s,&
                                                                   rho_dbcsr
      TYPE(mo_set_type), DIMENSION(:), POINTER                  :: mos
      INTEGER                                                   :: i, re, im
      ! Extract the other types first

      CALL get_qs_env(qs_env,&
                      bs_env=bs_env,&
                      rtp=rtp,&
                      matrix_s=matrix_s,&
                      mos=mos,&
                      dft_control=tdagw_env%dft_control)

      ! Number of spins
      tdagw_env%n_spin = bs_env%n_spin
      ! Number of atomic orbitals
      tdagw_env%n_ao = bs_env%n_ao
      ! Number of auxiliary basis orbitals
      tdagw_env%n_RI = bs_env%n_RI
      ! Default field is zero
      tdagw_env%field(:) = 0.0_dp
      ! Default time is zero
      tdagw_env%sim_step = 0
      tdagw_env%sim_time = 0
      ! Time step is taken from rtp
      tdagw_env%sim_dt = rtp%dt
      ! Threshold for etrs is taken from the eps_energy from RT propagation
      tdagw_env%etrs_threshold = tdagw_env%dft_control%rtp_control%eps_ener
      ! Allocate and determine the expansion orders, based on arguments used in TDDFT
      ! TODO : Last two arguments are indicating that we are using Taylor expansion and not doing ion motion
      ALLOCATE(tdagw_env%orders(2, tdagw_env%n_spin))
      CALL get_nsquare_norder(2.0_dp*MAX(ABS(MAXVAL(bs_env%eigenval_G0W0)),&
                              ABS(MINVAL(bs_env%eigenval_G0W0))),&
                              tdagw_env%orders(1,1), tdagw_env%orders(2,1),&
                              tdagw_env%dft_control%rtp_control%eps_exp,&
                              1, .FALSE.)
      IF (tdagw_env%n_spin > 1) THEN
         DO i=2,tdagw_env%n_spin
            tdagw_env%orders(1,i) = tdagw_env%orders(1,1)
            tdagw_env%orders(2,i) = tdagw_env%orders(2,1)
         END DO
      END IF
      ! Sim start index and total number of steps as well
      tdagw_env%sim_start = rtp%i_start
      tdagw_env%sim_nsteps = rtp%nsteps

      ! Construct moments from dbcsr
      ALLOCATE(moments_dbcsr_p(3))
      DO i=1,3
         ! Make sure the pointer is empty
         NULLIFY(moments_dbcsr_p(i)%matrix)
         ! Allocate a new matrix that the pointer points to
         ALLOCATE(moments_dbcsr_p(i)%matrix)
         ! Create the matrix storage - matrix copies the structure of overlap matrix
         ! CALL dbcsr_create(moments_dbcsr_p(i)%matrix, name="moments", template=matrix_s(1)%matrix)
         CALL dbcsr_copy(moments_dbcsr_p(i)%matrix, matrix_s(1)%matrix)
      END DO
      ! Run the moment calculation
      CALL build_local_moment_matrix(qs_env, moments_dbcsr_p, 1)
      ! Copy to full matrix
      ALLOCATE(tdagw_env%moments(3))
      DO i=1,3
         ! Again, matrices are created from overlap template
         ALLOCATE(tdagw_env%moments(i)%matrix)
         CALL cp_fm_create(tdagw_env%moments(i)%matrix, bs_env%fm_s_Gamma%matrix_struct)
         CALL copy_dbcsr_to_fm(moments_dbcsr_p(i)%matrix, tdagw_env%moments(i)%matrix)
      END DO

      ! Now can deallocate dbcsr matrices
      DO i=1,3
         CALL dbcsr_release(moments_dbcsr_p(i)%matrix)
      END DO
      DEALLOCATE(moments_dbcsr_p)

      ! Initiate density through results of DFT
      ! TODO : Handle cases, when initial method uses density only? Do such cases exist?
      CALL rt_initialize_rho_from_mos(rtp, mos)
      CALL get_rtp(rtp, rho_new=rho_dbcsr)
      ! Copy the density
      ALLOCATE(tdagw_env%rho(2*tdagw_env%n_spin))
      DO i=1,2*tdagw_env%n_spin
         ALLOCATE(tdagw_env%rho(i)%matrix)
         CALL cp_fm_create(tdagw_env%rho(i)%matrix, bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL copy_dbcsr_to_fm(rho_dbcsr(i)%matrix, tdagw_env%rho(i)%matrix)
      END DO
      ! Create the inverse overlap matrix, for use in density propagation 
      CALL s_matrices_create(matrix_s, rtp)
      CALL get_rtp(rtp, S_inv=tdagw_env%S_inv)

      ! Copy the single particle hamiltonian
      ALLOCATE(tdagw_env%ham_single_particle(tdagw_env%n_spin))
      @:SPIN_DO(i, re, im, tdagw_env%n_spin)
         ALLOCATE(tdagw_env%ham_single_particle(i)%matrix)
         CALL cp_fm_create(tdagw_env%ham_single_particle(i)%matrix, bs_env%fm_h_G0W0_Gamma%matrix_struct)
         ! TODO : So far no spin resolution in Hamiltonian of GW
         CALL cp_fm_to_fm(bs_env%fm_h_G0W0_Gamma, tdagw_env%ham_single_particle(i)%matrix)
      END DO

      ! Create the matrices and workspaces for effective Hamiltonian
      ALLOCATE(tdagw_env%ham_effective(2*tdagw_env%n_spin))
      ALLOCATE(tdagw_env%ham_workspace(2*tdagw_env%n_spin))
      ALLOCATE(tdagw_env%rho_new(2*tdagw_env%n_spin))
      ALLOCATE(tdagw_env%rho_new_last(2*tdagw_env%n_spin))
      ALLOCATE(tdagw_env%rho_M(2*tdagw_env%n_spin))
      DO i=1,2*tdagw_env%n_spin
         ALLOCATE(tdagw_env%ham_effective(i)%matrix)
         CALL cp_fm_create(tdagw_env%ham_effective(i)%matrix, bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_set_all(tdagw_env%ham_effective(i)%matrix, 0.0_dp)
         ALLOCATE(tdagw_env%ham_workspace(i)%matrix)
         CALL cp_fm_create(tdagw_env%ham_workspace(i)%matrix, bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_set_all(tdagw_env%ham_workspace(i)%matrix, 0.0_dp)
         ALLOCATE(tdagw_env%rho_new(i)%matrix)
         CALL cp_fm_create(tdagw_env%rho_new(i)%matrix, bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_set_all(tdagw_env%rho_new(i)%matrix, 0.0_dp)
         ALLOCATE(tdagw_env%rho_new_last(i)%matrix)
         CALL cp_fm_create(tdagw_env%rho_new_last(i)%matrix, bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_set_all(tdagw_env%rho_new_last(i)%matrix, 0.0_dp)
         ALLOCATE(tdagw_env%rho_M(i)%matrix)
         CALL cp_fm_create(tdagw_env%rho_M(i)%matrix, bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_set_all(tdagw_env%rho_M(i)%matrix, 0.0_dp)
      END DO

      ! Allocate space for density propagation
      ALLOCATE(tdagw_env%rho_workspace(3))
      DO i=1,3
         ALLOCATE(tdagw_env%rho_workspace(i)%matrix)
         CALL cp_fm_create(tdagw_env%rho_workspace(i)%matrix, bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_set_all(tdagw_env%rho_workspace(i)%matrix, 0.0_dp)
      END DO

      ! Self-energy components
      ALLOCATE(tdagw_env%hartree_dbt)

      ! Allocate self-energy parts and dynamic Hartree potential
      ! TODO : Hartree probably has spin independent form
      ALLOCATE(tdagw_env%hartree_orig(2*tdagw_env%n_spin))
      ALLOCATE(tdagw_env%hartree_curr(2*tdagw_env%n_spin))
      ALLOCATE(tdagw_env%sigma_SEX(2*tdagw_env%n_spin))
      ALLOCATE(tdagw_env%sigma_COH(2*tdagw_env%n_spin))
      DO i=1,2*tdagw_env%n_spin
         ALLOCATE(tdagw_env%hartree_orig(i)%matrix)
         ALLOCATE(tdagw_env%hartree_curr(i)%matrix)
         ALLOCATE(tdagw_env%sigma_SEX(i)%matrix)
         ALLOCATE(tdagw_env%sigma_COH(i)%matrix)
         CALL cp_fm_create(tdagw_env%sigma_COH(i)%matrix, bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_create(tdagw_env%sigma_SEX(i)%matrix, bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_create(tdagw_env%hartree_orig(i)%matrix, bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_create(tdagw_env%hartree_curr(i)%matrix, bs_env%fm_h_G0W0_Gamma%matrix_struct)
         CALL cp_fm_set_all(tdagw_env%sigma_COH(i)%matrix, 0.0_dp)
         CALL cp_fm_set_all(tdagw_env%sigma_SEX(i)%matrix, 0.0_dp)
         CALL cp_fm_set_all(tdagw_env%hartree_orig(i)%matrix, 0.0_dp)
         CALL cp_fm_set_all(tdagw_env%hartree_curr(i)%matrix, 0.0_dp)
      END DO
      
   END SUBROUTINE

   SUBROUTINE release_tdagw_env(tdagw_env)
      TYPE(tdagw_env_type)                                      :: tdagw_env
      CHARACTER(len=*), PARAMETER                               :: routineN = "release_tdagw_env"
      ! INTEGER                                                   :: i

      ! DO i=1,3
      !    CALL cp_fm_release(tdagw_env%moments(i))
      ! END DO
      ! DEALLOCATE(tdagw_env%moments)

      ! DO i=1,2*tdagw_env%n_spin
      !    CALL cp_fm_release(tdagw_env%rho(i))
      !    CALL cp_fm_release(tdagw_env%ham_effective(i))
      !    CALL cp_fm_release(tdagw_env%ham_workspace(i))
      !    CALL cp_fm_release(tdagw_env%sigma_COH(i))
      !    CALL cp_fm_release(tdagw_env%sigma_SEX(i))
      !    CALL cp_fm_release(tdagw_env%hartree_orig(i))
      !    CALL cp_fm_release(tdagw_env%hartree_curr(i))
      ! END DO
      ! DEALLOCATE(tdagw_env%rho)
      ! DEALLOCATE(tdagw_env%ham_effective)
      ! DEALLOCATE(tdagw_env%ham_workspace)
      ! DEALLOCATE(tdagw_env%sigma_COH)
      ! DEALLOCATE(tdagw_env%sigma_SEX)
      ! DEALLOCATE(tdagw_env%hartree_orig)
      ! DEALLOCATE(tdagw_env%hartree_curr)
      CALL cp_fm_release(tdagw_env%rho)
      CALL cp_fm_release(tdagw_env%ham_effective)
      CALL cp_fm_release(tdagw_env%ham_workspace)
      CALL cp_fm_release(tdagw_env%sigma_COH)
      CALL cp_fm_release(tdagw_env%sigma_SEX)
      CALL cp_fm_release(tdagw_env%hartree_orig)
      CALL cp_fm_release(tdagw_env%hartree_curr)
      CALL cp_fm_release(tdagw_env%ham_single_particle)
      CALL cp_fm_release(tdagw_env%rho_workspace)
      CALL cp_fm_release(tdagw_env%rho_new)
      CALL cp_fm_release(tdagw_env%rho_new_last)
      CALL cp_fm_release(tdagw_env%rho_M)

      CALL dbt_destroy(tdagw_env%hartree_dbt)

      NULLIFY(tdagw_env%S_inv)
      DEALLOCATE(tdagw_env%orders)
      
   END SUBROUTINE
END MODULE

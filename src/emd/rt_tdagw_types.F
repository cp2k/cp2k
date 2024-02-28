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
                                              cp_fm_write_formatted,&
                                              cp_fm_to_fm
   USE cp_fm_basic_linalg,              ONLY: cp_fm_invert
   USE dbcsr_api,                       ONLY: dbcsr_type,&
                                              dbcsr_p_type,&
                                              dbcsr_create,&
                                              dbcsr_release,&
                                              dbcsr_print,&
                                              dbcsr_copy
   USE dbt_api,                         ONLY: dbt_type,&
                                              dbt_pgrid_type,&
                                              dbt_pgrid_create,&
                                              dbt_pgrid_destroy,&
                                              dbt_mp_environ_pgrid,&
                                              dbt_default_distvec,&
                                              dbt_distribution_type,&
                                              dbt_distribution_new,&
                                              dbt_distribution_destroy,&
                                              dbt_create,&
                                              dbt_copy_matrix_to_tensor,&
                                              dbt_destroy
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
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
   USE moments_utils,                   ONLY: get_reference_point
   USE matrix_exp,                      ONLY: get_nsquare_norder
   USE input_constants,                 ONLY: use_mom_ref_coac,&
                                              rtp_tdagw_ham_g0w0,&
                                              rtp_tdagw_ham_ks,&
                                              rtp_tdagw_hartree_dbt,&
                                              rtp_tdagw_hartree_pw
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get

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
      TYPE(dbcsr_type), POINTER                                 :: S_inv,&
                                                                   rho_dbcsr
      TYPE(dbt_type), POINTER                                   :: hartree_dbt,&
                                                                   screened_dbt,&
                                                                   sigma_dbt,&
                                                                   greens_dbt,&
                                                                   t_two_plus_one,&
                                                                   t_one_plus_one,&
                                                                   t_3c_w,&
                                                                   t_3c_work_RI__AO_AO,&
                                                                   t_3c_work_RI_AO__AO,&
                                                                   t_3c_work2_RI_AO__AO
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: sigma_SEX, sigma_COH
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: hartree_orig, hartree_curr
      INTEGER, DIMENSION(:,:), POINTER                          :: orders
      INTEGER                                                   :: hartree_method
      INTEGER                                                   :: ham_single_particle_type

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
      REAL(kind=dp), DIMENSION(:), POINTER                      :: custom_ref_point
      REAL(kind=dp), DIMENSION(3)                               :: rpoint
      ! TYPE(cp_fm_type)                                          :: S_inv_fm
      TYPE(section_vals_type), POINTER                          :: input
      ! Extract the other types first

      CALL get_qs_env(qs_env,&
                      bs_env=bs_env,&
                      rtp=rtp,&
                      matrix_s=matrix_s,&
                      mos=mos,&
                      dft_control=tdagw_env%dft_control,&
                      input=input)

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
      ! Recover custom options
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%HARTREE_METHOD",&
                                i_val=tdagw_env%hartree_method)
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%TDAGW_HAMILTONIAN",&
                                i_val=tdagw_env%ham_single_particle_type)
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
      ! TODO : Read custom reference point from input
      ! TODO : Maybe use the moments section that is present in DFT?
      ALLOCATE(custom_ref_point(3), source=0.0_dp)
      rpoint(:) = 0.0_dp
      CALL get_reference_point(rpoint, qs_env=qs_env, reference=use_mom_ref_coac, ref_point=custom_ref_point)
      DEALLOCATE(custom_ref_point)
      CALL build_local_moment_matrix(qs_env, moments_dbcsr_p, 1, rpoint)
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
      ! TODO : Check that dbcsr overlap inverse coincides with fm overlap inverse
      ! CALL dbcsr_print(tdagw_env%S_inv)
      ! CALL cp_fm_create(S_inv_fm, bs_env%fm_s_Gamma%matrix_struct)
      ! CALL cp_fm_invert(bs_env%fm_s_Gamma, S_inv_fm)
      ! CALL cp_fm_write_formatted(S_inv_fm, unit=6)
      ! CALL cp_fm_release(S_inv_fm)

      ! Copy the single particle hamiltonian
      ALLOCATE(tdagw_env%ham_single_particle(tdagw_env%n_spin))
      @:SPIN_DO(i, re, im, tdagw_env%n_spin)
         ALLOCATE(tdagw_env%ham_single_particle(i)%matrix)
         CALL cp_fm_create(tdagw_env%ham_single_particle(i)%matrix, bs_env%fm_h_G0W0_Gamma%matrix_struct)
         ! TODO : So far no spin resolution in Hamiltonian of GW
         CALL cp_fm_to_fm(bs_env%fm_h_G0W0_Gamma, tdagw_env%ham_single_particle(i)%matrix)
         ! CALL cp_fm_to_fm(bs_env%fm_ks_Gamma(i), tdagw_env%ham_single_particle(i)%matrix)
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

      ! Allocate workspaces for get_sigma
      CALL create_sigma_workspace(tdagw_env, qs_env)

      ! Depending on the chosen methods, allocate extra workspace
      IF ( tdagw_env%hartree_method == rtp_tdagw_hartree_dbt ) THEN
         CALL create_hartree_dbt_workspace(tdagw_env, qs_env)
      END IF
      
   END SUBROUTINE

! **************************************************************************************************
!> \brief Releases the environment allocated structures
!> \param tdagw_env 
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE release_tdagw_env(tdagw_env)
      TYPE(tdagw_env_type)                                      :: tdagw_env
      CHARACTER(len=*), PARAMETER                               :: routineN = "release_tdagw_env"
      ! INTEGER                                                   :: i

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

      CALL release_sigma_workspace(tdagw_env)

      ! Depending on methods, deallocate additional workspace
      IF ( tdagw_env%hartree_method == rtp_tdagw_hartree_dbt ) THEN
         CALL release_hartree_dbt_workspace(tdagw_env)
      END IF
      
   END SUBROUTINE
! **************************************************************************************************
!> \brief Allocates the workspaces for Hartree dbt method, mainly the density tensors
!> \note DBT method uses 2D and 3D tensors, where one dimension has size 1, in order to emulate vectors
!> \param tdagw_env 
!> \param qs_env Quickstep environment - entry point of calculation
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE create_hartree_dbt_workspace(tdagw_env, qs_env)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "create_hartree_dbt_workspace"
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, ALLOCATABLE, DIMENSION(:,:)               :: distvec
      INTEGER, DIMENSION(1)                              :: distvec_one
      TYPE(dbt_distribution_type)                        :: dist
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pdims, pcoord
      TYPE(dbt_pgrid_type)                               :: pgrid
      INTEGER                                            :: j
      ! Create a 1D tensor (vector) for the result of contraction with density matrix
      ! Taken from gw_utils
      ! Get bs_env
      CALL get_qs_env(qs_env, bs_env=bs_env)

      ! Create the 3-index tensor with one dimension of size 1 from tensor primitives
      ! Create 3-index process grid
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
      NULLIFY(tdagw_env%t_two_plus_one)
      ALLOCATE(tdagw_env%t_two_plus_one)
      CALL dbt_create(tdagw_env%t_two_plus_one, "2+1 index tensor", dist, [1], [2, 3],&
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
      NULLIFY(tdagw_env%t_one_plus_one)
      ALLOCATE(tdagw_env%t_one_plus_one)
      CALL dbt_create(tdagw_env%t_one_plus_one, "tensor vector", dist, [1], [2], bs_env%sizes_RI, [1])
      ! No longer need the process information
      CALL dbt_pgrid_destroy(pgrid)
      CALL dbt_distribution_destroy(dist)
      DEALLOCATE(distvec)
      DEALLOCATE(pcoord)
      DEALLOCATE(pdims)

      ! Rho_dbcsr workspace - a single rho_dbcsr matrix
      NULLIFY(tdagw_env%rho_dbcsr)
      ALLOCATE(tdagw_env%rho_dbcsr)
      CALL dbcsr_create(tdagw_env%rho_dbcsr, name="Sparse density", template=bs_env%mat_ao_ao%matrix)

      ! 3c tensors - do not populate yet, just allocate distributions
      NULLIFY(tdagw_env%t_3c_work_RI__AO_AO)
      ALLOCATE(tdagw_env%t_3c_work_RI__AO_AO)
      CALL dbt_create(bs_env%t_RI__AO_AO, tdagw_env%t_3c_work_RI__AO_AO)
   END SUBROUTINE
! **************************************************************************************************
!> \brief Releases the workspaces for Hartree dbt method
!> \param tdagw_env 
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE release_hartree_dbt_workspace(tdagw_env)
      TYPE(tdagw_env_type)                               :: tdagw_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "release_hartree_dbt_workspace"
      IF (ASSOCIATED(tdagw_env%t_two_plus_one)) THEN
         CALL dbt_destroy(tdagw_env%t_two_plus_one)
         DEALLOCATE(tdagw_env%t_two_plus_one)
      END IF
      IF (ASSOCIATED(tdagw_env%t_one_plus_one)) THEN
         CALL dbt_destroy(tdagw_env%t_one_plus_one)
         DEALLOCATE(tdagw_env%t_one_plus_one)
      END IF
      IF (ASSOCIATED(tdagw_env%rho_dbcsr)) THEN
         CALL dbcsr_release(tdagw_env%rho_dbcsr)
         DEALLOCATE(tdagw_env%rho_dbcsr)
      END IF
      IF (ASSOCIATED(tdagw_env%t_3c_work_RI__AO_AO)) THEN
         CALL dbt_destroy(tdagw_env%t_3c_work_RI__AO_AO)
         DEALLOCATE(tdagw_env%t_3c_work_RI__AO_AO)
      END IF
   END SUBROUTINE
! **************************************************************************************************
!> \brief Allocates the workspaces for self-energy determination routine
!> \param tdagw_env Structure for holding information and workspace structures 
!> \param qs_env Quickstep environment - entry point of calculation
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE create_sigma_workspace(tdagw_env, qs_env)
      TYPE(tdagw_env_type)                               :: tdagw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "create_sigma_workspace"
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_type)                                   :: w_dbcsr

      CALL get_qs_env(qs_env, bs_env=bs_env)

      ! t_3c_w
      NULLIFY(tdagw_env%t_3c_w)
      ALLOCATE(tdagw_env%t_3c_w)
      CALL dbt_create(bs_env%t_RI__AO_AO, tdagw_env%t_3c_w)
      ! t_3c_work_RI_AO__AO
      NULLIFY(tdagw_env%t_3c_work_RI_AO__AO)
      ALLOCATE(tdagw_env%t_3c_work_RI_AO__AO)
      CALL dbt_create(bs_env%t_RI_AO__AO, tdagw_env%t_3c_work_RI_AO__AO)
      ! t_3c_work2_RI_AO__AO
      NULLIFY(tdagw_env%t_3c_work2_RI_AO__AO)
      ALLOCATE(tdagw_env%t_3c_work2_RI_AO__AO)
      CALL dbt_create(bs_env%t_RI_AO__AO, tdagw_env%t_3c_work2_RI_AO__AO)
      ! t_W
      NULLIFY(tdagw_env%screened_dbt)
      ALLOCATE(tdagw_env%screened_dbt)
      CALL dbt_create(bs_env%t_W, tdagw_env%screened_dbt)
      ! Populate t_W from gw run
      CALL dbcsr_create(w_dbcsr, name="W", template=bs_env%mat_RI_RI%matrix)
      CALL copy_fm_to_dbcsr(bs_env%fm_W_MIC_freq_zero, w_dbcsr)
      CALL dbt_copy_matrix_to_tensor(w_dbcsr, tdagw_env%screened_dbt)
      CALL dbcsr_release(w_dbcsr)
      ! sigma_dbt
      NULLIFY(tdagw_env%sigma_dbt)
      ALLOCATE(tdagw_env%sigma_dbt)
      CALL dbt_create(bs_env%t_G, tdagw_env%sigma_dbt)
      ! greens_dbt
      NULLIFY(tdagw_env%greens_dbt)
      ALLOCATE(tdagw_env%greens_dbt)
      CALL dbt_create(bs_env%t_G, tdagw_env%greens_dbt)
   END SUBROUTINE
! **************************************************************************************************
!> \brief Releases the workspaces for self-energy determination
!> \param tdagw_env 
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE release_sigma_workspace(tdagw_env)
      TYPE(tdagw_env_type)                               :: tdagw_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "release_sigma_workspace"
      
      IF (ASSOCIATED(tdagw_env%t_3c_w)) THEN 
         CALL dbt_destroy(tdagw_env%t_3c_w)
         DEALLOCATE(tdagw_env%t_3c_w)
      END IF
      IF (ASSOCIATED(tdagw_env%t_3c_work_RI_AO__AO)) THEN 
         CALL dbt_destroy(tdagw_env%t_3c_work_RI_AO__AO)
         DEALLOCATE(tdagw_env%t_3c_work_RI_AO__AO)
      END IF
      IF (ASSOCIATED(tdagw_env%t_3c_work2_RI_AO__AO)) THEN 
         CALL dbt_destroy(tdagw_env%t_3c_work2_RI_AO__AO)
         DEALLOCATE(tdagw_env%t_3c_work2_RI_AO__AO)
      END IF
      IF (ASSOCIATED(tdagw_env%screened_dbt)) THEN
         CALL dbt_destroy(tdagw_env%screened_dbt)
         DEALLOCATE(tdagw_env%screened_dbt)
      END IF
      IF (ASSOCIATED(tdagw_env%sigma_dbt)) THEN
         CALL dbt_destroy(tdagw_env%sigma_dbt)
         DEALLOCATE(tdagw_env%sigma_dbt)
      END IF
      IF (ASSOCIATED(tdagw_env%greens_dbt)) THEN
         CALL dbt_destroy(tdagw_env%greens_dbt)
         DEALLOCATE(tdagw_env%greens_dbt)
      END IF
   END SUBROUTINE
END MODULE

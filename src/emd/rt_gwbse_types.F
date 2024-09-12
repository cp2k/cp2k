!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Data storage and other types for time propagation via GWBSE method.
!> \author Stepan Marek (01.24)
! **************************************************************************************************

MODULE rt_gwbse_types

   USE kinds,                           ONLY: dp
   USE cp_fm_types,                     ONLY: cp_fm_type,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_create,&
                                              cp_fm_set_all,&
                                              cp_fm_write_formatted,&
                                              cp_fm_to_fm
   USE cp_cfm_types,                    ONLY: cp_cfm_p_type,&
                                              cp_cfm_type,&
                                              cp_cfm_set_all,&
                                              cp_cfm_create,&
                                              cp_fm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_release
   USE cp_fm_basic_linalg,              ONLY: cp_fm_invert,&
                                              cp_fm_transpose,&
                                              cp_fm_column_scale,&
                                              cp_fm_scale_and_add
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale
   USE cp_dbcsr_api,                    ONLY: dbcsr_type,&
                                              dbcsr_p_type,&
                                              dbcsr_create,&
                                              dbcsr_release,&
                                              dbcsr_print,&
                                              dbcsr_get_info,&
                                              dbcsr_copy,&
                                              dbcsr_set,&
                                              dbcsr_scale,&
                                              dbcsr_type_complex_8
   USE parallel_gemm_api,               ONLY: parallel_gemm
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
                                              dbt_get_num_blocks,&
                                              dbt_destroy
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE qs_mo_types,                     ONLY: mo_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE qs_environment_types,            ONLY: qs_environment_type,&
                                              get_qs_env
   USE force_env_types,                 ONLY: force_env_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: create_and_init_bs_env
   USE rt_propagation_types,            ONLY: rt_prop_type,&
                                              get_rtp
   USE rt_propagator_init,              ONLY: rt_initialize_rho_from_mos
   USE rt_propagation_methods,          ONLY: s_matrices_create
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE moments_utils,                   ONLY: get_reference_point
   USE matrix_exp,                      ONLY: get_nsquare_norder
   USE gw_integrals,                    ONLY: build_3c_integral_block
   USE gw_large_cell_Gamma,             ONLY: compute_3c_integrals
!    USE gw_utils,                        ONLY: create_and_init_bs_env_for_gw,&
!                                               setup_AO_and_RI_basis_set,&
!                                               get_RI_basis_and_basis_function_indices,&
!                                               set_heuristic_parameters,&
!                                               set_parallelization_parameters,&
!                                               allocate_and_fill_matrices_and_arrays,&
!                                               create_tensors
   USE libint_wrapper,                  ONLY: cp_libint_static_init
   USE input_constants,                 ONLY: use_mom_ref_coac,&
                                              use_mom_ref_zero,&
                                              use_mom_ref_user,&
                                              use_mom_ref_com,&
                                              rtp_gwbse_ham_g0w0,&
                                              rtp_gwbse_ham_ks,&
                                              rtp_gwbse_hartree_dbt,&
                                              rtp_gwbse_hartree_pw,&
                                              do_taylor,&
                                              do_bch,&
                                              do_exact
   USE physcon,                         ONLY: angstrom
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_get_subs_vals

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_gwbse"

   #:include "rt_gwbse_macros.fypp"
   
   PUBLIC :: gwbse_env_type,&
             create_gwbse_env,&
             release_gwbse_env,&
             multiply_cfm_fm,&
             multiply_fm_cfm

! **************************************************************************************************
!> \param n_spin Number of spin channels that are present
!> \param n_ao Number of atomic orbitals
!> \param n_RI Number of RI orbitals
!> \param n_occ Number of occupied orbitals, spin dependent
!> \param spin_degeneracy Number of electrons per orbital
!> \param field Electric field calculated at the given timestep
!> \param moments Moment operators along cartesian directions - centered at zero charge - used for plotting
!> \param moments_field Moment operators along cartesian directions - used to coupling to the field -
!>                origin bound to unit cell
!> \param sim_step Current step of the simulation
!> \param sim_start Starting step of the simulation
!> \param sim_nsteps Number of steps of the simulation
!> \param sim_time Current time of the simulation
!> \param sim_dt Timestep of the simulation
!> \param etrs_threshold Self-consistency threshold for enforced time reversal symmetry propagation
!> \param exp_accuracy Threshold for matrix exponential calculation
!> \param dft_control DFT control parameters
!> \param ham_effective Real and imaginary part of the effective Hamiltonian used to propagate
!>                      the density matrix
!> \param ham_reference Reference Hamiltonian, which does not change in the propagation = DFT+G0W0 - initial Hartree - initial COHSEX
!> \param ham_workspace Workspace matrices for use with the Hamiltonian propagation - storage of
!>                      exponential propagators etc.
!> \param rho Density matrix at the current time step
!> \param rho_new Density matrix - workspace in ETRS
!> \param rho_last Density matrix - workspace in ETRS
!> \param rho_new_last Density matrix - workspace in ETRS
!> \param rho_M Density matrix - workspace in ETRS
!> \param S_inv_fm Inverse overlap matrix, full matrix
!> \param S_fm Overlap matrix, full matrix
!> \param S_inv Inverse overlap matrix, sparse matrix
!> \param rho_dbcsr Density matrix, sparse matrix
!> \param rho_workspace Matrices for storage of density matrix at different timesteps for
!>                      interpolation and self-consistency checks etc.
!> \param complex_workspace Workspace for complex density (exact diagonalisation)
!> \param complex_s Complex overlap matrix (exact diagonalisation)
!> \param real_eigvals Eigenvalues of hermitian matrix (exact diagonalisation)
!> \param exp_eigvals Exponentiated eigenvalues (exact diagonalisation)
!> \param v_dbcsr Sparse matrix with bare Coulomb in RI basis
!> \param w_dbcsr Sparse matrix with correlation part of dressed Coulomb in RI basis (without bare Coulomb)
!> \param screened_dbt Tensor for screened Coulomb interaction
!> \param sigma_dbt Tensor for self-energy
!> \param greens_dbt Tensor for greens function/density matrix
!> \param t_3c_w Tensor containing 3c integrals
!> \param t_3c_work_RI__AO_AO Tensor sigma contraction
!> \param t_3c_work_RI_AO__AO Tensor sigma contraction
!> \param t_3c_work2_RI_AO__AO Tensor sigma contraction
!> \param sigma_SEX Screened exchange self-energy
!> \param sigma_COH Coulomb hole self-energy
!> \param hartree_curr Current Hartree matrix
!> \param etrs_max_iter Maximum number of ETRS iterations
!> \param hartree_method Which method to use for hartree_calculation
!> \param ham_reference_type Which Hamiltonian to use as single particle basis
!> \param mat_exp_method Which method to use for matrix exponentiation
!> \param unit_nr Number of output unit
!> \param int_3c_array Array containing the local 3c integrals
!> \author Stepan Marek (01.24)
! **************************************************************************************************
   TYPE gwbse_env_type
      INTEGER                                                   :: n_spin,&
                                                                   n_ao,&
                                                                   n_RI
      INTEGER, DIMENSION(2)                                     :: n_occ
      REAL(kind=dp)                                             :: spin_degeneracy
      REAL(kind=dp), DIMENSION(3)                               :: field
      TYPE(cp_fm_p_type), DIMENSION(3)                          :: moments,&
                                                                   moments_field
      INTEGER                                                   :: sim_step,&
                                                                   sim_start,&
                                                                   sim_nsteps
      REAL(kind=dp)                                             :: sim_time,&
                                                                   sim_dt,&
                                                                   etrs_threshold,&
                                                                   exp_accuracy
      TYPE(dft_control_type), POINTER                           :: dft_control
      ! Stores data needed for reading/writing to the restart files
      TYPE(section_vals_type), POINTER                          :: restart_section
      LOGICAL                                                   :: restart_extracted

      ! Different indices signify different spins
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER                :: ham_effective
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER                :: ham_reference
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER                :: ham_workspace
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER                :: sigma_SEX
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: sigma_COH, hartree_curr

      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER                :: rho,&
                                                                   rho_new,&
                                                                   rho_new_last,&
                                                                   rho_M,&
                                                                   rho_orig
      TYPE(cp_fm_p_type)                                        :: S_inv_fm,&
                                                                   S_fm
      ! Many routines require overlap in the complex format
      TYPE(cp_cfm_p_type)                                       :: S_cfm
      TYPE(dbcsr_type), POINTER                                 :: rho_dbcsr
      ! TODO : Check how much is actually needed and whether it is not better to separate names
      ! Indices only correspond to different workspaces
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER                :: rho_workspace
      ! Many methods use real and imaginary parts separately - prevent unnecessary reallocation
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: real_workspace
      ! Workspace required for exact matrix exponentiation
      REAL(kind=dp), DIMENSION(:), POINTER                      :: real_eigvals
      COMPLEX(kind=dp), DIMENSION(:), POINTER                   :: exp_eigvals
      ! Workspace required for hartree_pw
      TYPE(dbcsr_p_type)                                        :: v_dbcsr,&
                                                                   w_dbcsr
      ! TODO : Spins, k-points
      TYPE(dbt_type), POINTER                                   :: screened_dbt,&
                                                                   sigma_dbt,&
                                                                   greens_dbt,&
                                                                   t_3c_w,&
                                                                   t_3c_work_RI__AO_AO,&
                                                                   t_3c_work_RI_AO__AO,&
                                                                   t_3c_work2_RI_AO__AO
      ! These matrices are always real
      INTEGER                                                   :: etrs_max_iter
      INTEGER                                                   :: hartree_method
      INTEGER                                                   :: ham_reference_type
      INTEGER                                                   :: mat_exp_method
      INTEGER                                                   :: unit_nr 
      REAL(kind=dp), DIMENSION(:,:,:), POINTER                  :: int_3c_array

   END TYPE gwbse_env_type

CONTAINS

! **************************************************************************************************
!> \brief Allocates structures and prepares gwbse_env for run
!> \param gwbse_env gwbse_env_type that is initialised
!> \param qs_env Entry point of the calculation
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE create_gwbse_env(gwbse_env, qs_env, force_env)
      TYPE(gwbse_env_type), POINTER                             :: gwbse_env
      TYPE(qs_environment_type),POINTER                         :: qs_env
      TYPE(force_env_type), POINTER                             :: force_env
      CHARACTER(len=*), PARAMETER                               :: routineN="create_gwbse_env"
      TYPE(post_scf_bandstructure_type), POINTER                :: bs_env
      TYPE(rt_prop_type), POINTER                               :: rtp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER                 :: matrix_s
      TYPE(mo_set_type), DIMENSION(:), POINTER                  :: mos
      INTEGER                                                   :: i, k
      TYPE(section_vals_type), POINTER                          :: input, bs_sec, md_sec

      ! Allocate the storage for the gwbse environment
      NULLIFY(gwbse_env)
      ALLOCATE(gwbse_env)
      ! Extract the other types first
      CALL get_qs_env(qs_env,&
                      bs_env=bs_env,&
                      rtp=rtp,&
                      matrix_s=matrix_s,&
                      mos=mos,&
                      dft_control=gwbse_env%dft_control,&
                      input=input)
      bs_sec => section_vals_get_subs_vals(input, "PROPERTIES%BANDSTRUCTURE")
      ! TODO : Ability to run without the GW code?
      IF (.NOT. ASSOCIATED(bs_env)) THEN
         CPABORT("Cannot run RT-aGW without running GW calculation (PROPERTIES) before")
      END IF
      ! Number of spins
      gwbse_env%n_spin = bs_env%n_spin
      ! Number of atomic orbitals
      gwbse_env%n_ao = bs_env%n_ao
      ! Number of auxiliary basis orbitals
      gwbse_env%n_RI = bs_env%n_RI
      ! Number of occupied orbitals - for closed shell equals to half the number of electrons
      gwbse_env%n_occ(:) = bs_env%n_occ(:)
      ! Spin degeneracy - number of spins per orbital
      gwbse_env%spin_degeneracy = bs_env%spin_degeneracy
      ! Default field is zero
      gwbse_env%field(:) = 0.0_dp
      ! Default time is zero
      gwbse_env%sim_step = 0
      gwbse_env%sim_time = 0
      ! Time step is taken from rtp
      md_sec => section_vals_get_subs_vals(force_env%root_section, "MOTION%MD")
      CALL section_vals_val_get(md_sec, "TIMESTEP", r_val=gwbse_env%sim_dt)
      ! gwbse_env%sim_dt = rtp%dt
      ! Threshold for etrs is taken from the eps_energy from RT propagation
      gwbse_env%etrs_threshold = gwbse_env%dft_control%rtp_control%eps_ener
      ! TODO : Check whether dft_control can be dropped from gwbse_env
      gwbse_env%exp_accuracy = gwbse_env%dft_control%rtp_control%eps_exp
      ! Recover custom options
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%HARTREE_METHOD",&
                                i_val=gwbse_env%hartree_method)
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%GWBSE_HAMILTONIAN",&
                                i_val=gwbse_env%ham_reference_type)
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%MAX_ITER",&
                                i_val=gwbse_env%etrs_max_iter)
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%MAT_EXP",&
                                i_val=gwbse_env%mat_exp_method)
      ! Output unit number, recovered from the post_scf_bandstructure_type
      gwbse_env%unit_nr = bs_env%unit_nr
      ! IF (gwbse_env%unit_nr > 0) WRITE (gwbse_env%unit_nr, *) "exp_method", gwbse_env%mat_exp_method
      ! Allocate and determine the expansion orders, based on arguments used in TDDFT
      ! TODO : Last two arguments are indicating that we are using Taylor expansion and not doing ion motion
      ! Sim start index and total number of steps as well
      CALL section_vals_val_get(md_sec, "STEP_START_VAL", i_val=gwbse_env%sim_start)
      CALL section_vals_val_get(md_sec, "STEPS", i_val=gwbse_env%sim_nsteps)
      ! Get the restart section
      gwbse_env%restart_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION%PRINT%RESTART")
      gwbse_env%restart_extracted = .FALSE.

      ! Allocate moments matrices
      DO k=1,3
         ! Matrices are created from overlap template
         ! Values are initialized in initialize_gwbse_env
         NULLIFY(gwbse_env%moments(k)%matrix)
         ALLOCATE(gwbse_env%moments(k)%matrix)
         CALL cp_fm_create(gwbse_env%moments(k)%matrix, bs_env%fm_s_Gamma%matrix_struct)
         NULLIFY(gwbse_env%moments_field(k)%matrix)
         ALLOCATE(gwbse_env%moments_field(k)%matrix)
         CALL cp_fm_create(gwbse_env%moments_field(k)%matrix, bs_env%fm_s_Gamma%matrix_struct)
      END DO

      ! Allocate space for density propagation and other operations
      NULLIFY(gwbse_env%rho_workspace)
      ALLOCATE(gwbse_env%rho_workspace(4))
      DO i=1,SIZE(gwbse_env%rho_workspace)
         ALLOCATE(gwbse_env%rho_workspace(i)%matrix)
         CALL cp_cfm_create(gwbse_env%rho_workspace(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_cfm_set_all(gwbse_env%rho_workspace(i)%matrix, CMPLX(0.0, 0.0, kind=dp))
      END DO
      ! Allocate real workspace
      NULLIFY(gwbse_env%real_workspace)
      SELECT CASE (gwbse_env%mat_exp_method)
         CASE (do_exact)
            ALLOCATE(gwbse_env%real_workspace(4))
         CASE (do_bch)
            ALLOCATE(gwbse_env%real_workspace(2))
         CASE DEFAULT
            CPABORT("Only exact and BCH matrix propagation implemented in rt-GWBSE")
      END SELECT
      DO i=1,SIZE(gwbse_env%real_workspace)
         NULLIFY(gwbse_env%real_workspace(i)%matrix)
         ALLOCATE(gwbse_env%real_workspace(i)%matrix)
         CALL cp_fm_create(gwbse_env%real_workspace(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_set_all(gwbse_env%real_workspace(i)%matrix, 0.0_dp)
      END DO
      ! Allocate density matrix
      ALLOCATE(gwbse_env%rho(gwbse_env%n_spin))
      DO i=1,gwbse_env%n_spin
         NULLIFY(gwbse_env%rho(i)%matrix)
         ALLOCATE(gwbse_env%rho(i)%matrix)
         CALL cp_cfm_create(gwbse_env%rho(i)%matrix, matrix_struct=bs_env%fm_s_Gamma%matrix_struct)
      END DO
      ! Create the inverse overlap matrix, for use in density propagation 
      ! Start by creating the actual overlap matrix
      ! TODO : Handle different spins and k-points
      NULLIFY(gwbse_env%S_fm%matrix)
      ALLOCATE(gwbse_env%S_fm%matrix)
      CALL cp_fm_create(gwbse_env%S_fm%matrix, bs_env%fm_s_Gamma%matrix_struct)
      NULLIFY(gwbse_env%S_inv_fm%matrix)
      ALLOCATE(gwbse_env%S_inv_fm%matrix)
      CALL cp_fm_create(gwbse_env%S_inv_fm%matrix, bs_env%fm_s_Gamma%matrix_struct)
      NULLIFY(gwbse_env%S_cfm%matrix)
      ALLOCATE(gwbse_env%S_cfm%matrix)
      CALL cp_cfm_create(gwbse_env%S_cfm%matrix, bs_env%fm_s_Gamma%matrix_struct)

      ! Create the single particle hamiltonian
      ! Allocate workspace
      NULLIFY(gwbse_env%ham_workspace)
      ALLOCATE(gwbse_env%ham_workspace(gwbse_env%n_spin))
      DO i=1,gwbse_env%n_spin
         ALLOCATE(gwbse_env%ham_workspace(i)%matrix)
         CALL cp_cfm_create(gwbse_env%ham_workspace(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_cfm_set_all(gwbse_env%ham_workspace(i)%matrix, CMPLX(0.0, 0.0, kind=dp))
      END DO
      ! Now onto the Hamiltonian itself
      ALLOCATE(gwbse_env%ham_reference(gwbse_env%n_spin))
      DO i=1,gwbse_env%n_spin
         ALLOCATE(gwbse_env%ham_reference(i)%matrix)
         CALL cp_cfm_create(gwbse_env%ham_reference(i)%matrix, bs_env%fm_ks_Gamma(i)%matrix_struct)
      END DO

      ! Create the matrices and workspaces for ETRS propagation
      ALLOCATE(gwbse_env%ham_effective(gwbse_env%n_spin))
      ALLOCATE(gwbse_env%rho_new(gwbse_env%n_spin))
      ALLOCATE(gwbse_env%rho_new_last(gwbse_env%n_spin))
      ALLOCATE(gwbse_env%rho_M(gwbse_env%n_spin))
      ALLOCATE(gwbse_env%rho_orig(gwbse_env%n_spin))
      DO i=1,gwbse_env%n_spin
         ALLOCATE(gwbse_env%ham_effective(i)%matrix)
         CALL cp_cfm_create(gwbse_env%ham_effective(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_cfm_set_all(gwbse_env%ham_effective(i)%matrix, CMPLX(0.0, 0.0, kind=dp))
         ALLOCATE(gwbse_env%rho_new(i)%matrix)
         CALL cp_cfm_create(gwbse_env%rho_new(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_cfm_set_all(gwbse_env%rho_new(i)%matrix, CMPLX(0.0, 0.0, kind=dp))
         ALLOCATE(gwbse_env%rho_new_last(i)%matrix)
         CALL cp_cfm_create(gwbse_env%rho_new_last(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_cfm_set_all(gwbse_env%rho_new_last(i)%matrix, CMPLX(0.0, 0.0, kind=dp))
         ALLOCATE(gwbse_env%rho_M(i)%matrix)
         CALL cp_cfm_create(gwbse_env%rho_M(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_cfm_set_all(gwbse_env%rho_M(i)%matrix, CMPLX(0.0, 0.0, kind=dp))
         ALLOCATE(gwbse_env%rho_orig(i)%matrix)
         CALL cp_cfm_create(gwbse_env%rho_orig(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
      END DO

      ! Fields for exact diagonalisation
      NULLIFY(gwbse_env%real_eigvals)
      ALLOCATE(gwbse_env%real_eigvals(gwbse_env%n_ao))
      gwbse_env%real_eigvals(:) = 0.0_dp
      NULLIFY(gwbse_env%exp_eigvals)
      ALLOCATE(gwbse_env%exp_eigvals(gwbse_env%n_ao))
      gwbse_env%exp_eigvals(:) = CMPLX(0.0, 0.0, kind=dp)

      ! Allocate self-energy parts and dynamic Hartree potential
      ALLOCATE(gwbse_env%hartree_curr(gwbse_env%n_spin))
      ALLOCATE(gwbse_env%sigma_SEX(gwbse_env%n_spin))
      ALLOCATE(gwbse_env%sigma_COH(gwbse_env%n_spin))
      DO i=1,gwbse_env%n_spin
         ALLOCATE(gwbse_env%hartree_curr(i)%matrix)
         ALLOCATE(gwbse_env%sigma_SEX(i)%matrix)
         ALLOCATE(gwbse_env%sigma_COH(i)%matrix)
         CALL cp_fm_create(gwbse_env%sigma_COH(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_cfm_create(gwbse_env%sigma_SEX(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_create(gwbse_env%hartree_curr(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_set_all(gwbse_env%sigma_COH(i)%matrix, 0.0_dp)
         CALL cp_cfm_set_all(gwbse_env%sigma_SEX(i)%matrix, CMPLX(0.0, 0.0, kind=dp))
         CALL cp_fm_set_all(gwbse_env%hartree_curr(i)%matrix, 0.0_dp)
      END DO

      ! Allocate workspaces for get_sigma
      CALL create_sigma_workspace(gwbse_env, qs_env)

      ! Depending on the chosen methods, allocate extra workspace
      CALL create_hartree_ri_workspace(gwbse_env, qs_env)
      
   END SUBROUTINE

! **************************************************************************************************
!> \brief Simple reimplementation of cp_fm_release_pp1 for complex matrices
!> \param matrices cp_cfm_p_type(:) 
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE cp_cfm_release_pp1(matrices)
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER                :: matrices
      INTEGER                                                   :: i

      DO i=1,SIZE(matrices)
         CALL cp_cfm_release(matrices(i)%matrix)
         DEALLOCATE(matrices(i)%matrix)
      END DO
      DEALLOCATE(matrices)
      NULLIFY(matrices)
   END SUBROUTINE cp_cfm_release_pp1

! **************************************************************************************************
!> \brief Releases the environment allocated structures
!> \param gwbse_env 
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE release_gwbse_env(gwbse_env)
      TYPE(gwbse_env_type), POINTER                             :: gwbse_env
      CHARACTER(len=*), PARAMETER                               :: routineN = "release_gwbse_env"
      INTEGER                                                   :: i

      CALL cp_cfm_release_pp1(gwbse_env%ham_effective)
      CALL cp_cfm_release_pp1(gwbse_env%ham_workspace)
      CALL cp_fm_release(gwbse_env%sigma_COH)
      CALL cp_cfm_release_pp1(gwbse_env%sigma_SEX)
      CALL cp_fm_release(gwbse_env%hartree_curr)
      CALL cp_cfm_release_pp1(gwbse_env%ham_reference)
      CALL cp_cfm_release_pp1(gwbse_env%rho)
      CALL cp_cfm_release_pp1(gwbse_env%rho_workspace)
      CALL cp_cfm_release_pp1(gwbse_env%rho_new)
      CALL cp_cfm_release_pp1(gwbse_env%rho_new_last)
      CALL cp_cfm_release_pp1(gwbse_env%rho_M)
      CALL cp_cfm_release_pp1(gwbse_env%rho_orig)
      CALL cp_fm_release(gwbse_env%real_workspace)
      CALL cp_fm_release(gwbse_env%S_inv_fm%matrix)
      DEALLOCATE(gwbse_env%S_inv_fm%matrix)
      CALL cp_fm_release(gwbse_env%S_fm%matrix)
      DEALLOCATE(gwbse_env%S_fm%matrix)
      CALL cp_cfm_release(gwbse_env%S_cfm%matrix)
      DEALLOCATE(gwbse_env%S_cfm%matrix)

      DO i=1,3
         CALL cp_fm_release(gwbse_env%moments(i)%matrix)
         DEALLOCATE(gwbse_env%moments(i)%matrix)
         NULLIFY(gwbse_env%moments(i)%matrix)
         CALL cp_fm_release(gwbse_env%moments_field(i)%matrix)
         DEALLOCATE(gwbse_env%moments_field(i)%matrix)
         NULLIFY(gwbse_env%moments_field(i)%matrix)
      END DO

      ! TODO : Decide on where and when to release this

      CALL release_sigma_workspace(gwbse_env)

      CALL release_hartree_ri_workspace(gwbse_env)

      DEALLOCATE(gwbse_env%real_eigvals)
      DEALLOCATE(gwbse_env%exp_eigvals)

      ! Deallocate the storage for the environment itself
      DEALLOCATE(gwbse_env)
      ! Nullify to make sure it is not used again
      NULLIFY(gwbse_env)
      
   END SUBROUTINE
! **************************************************************************************************
!> \brief Allocates the workspaces for Hartree RI method
!> \note RI method calculates the Hartree contraction without the use of DBT, as it cannot emulate vectors
!> \param gwbse_env 
!> \param qs_env Quickstep environment - entry point of calculation
!> \author Stepan Marek
!> \date 05.2024
! **************************************************************************************************
   SUBROUTINE create_hartree_ri_workspace(gwbse_env, qs_env)
      TYPE(gwbse_env_type)                              :: gwbse_env
      TYPE(qs_environment_type), POINTER                :: qs_env
      CHARACTER(len=*), PARAMETER                       :: routineN = "create_hartree_ri_workspace"
      TYPE(post_scf_bandstructure_type), POINTER        :: bs_env
      REAL(kind=dp)                                     :: size_mb

      CALL get_qs_env(qs_env, bs_env=bs_env)

      NULLIFY(gwbse_env%rho_dbcsr)
      ALLOCATE(gwbse_env%rho_dbcsr)
      CALL dbcsr_create(gwbse_env%rho_dbcsr, name="Sparse density", template=bs_env%mat_ao_ao%matrix)

      NULLIFY(gwbse_env%v_dbcsr%matrix)
      ALLOCATE(gwbse_env%v_dbcsr%matrix)

      ! TODO : Implement option/decision to not precompute all the 3c integrals
      size_mb = REAL(gwbse_env%n_ao * gwbse_env%n_ao * gwbse_env%n_RI * STORAGE_SIZE(size_mb)) / (1024_dp * 1024_dp)
      IF (gwbse_env%unit_nr > 0) WRITE(gwbse_env%unit_nr, *) "Approximate size of int_3c in MB", size_mb

      ALLOCATE(gwbse_env%int_3c_array(gwbse_env%n_ao, gwbse_env%n_ao, gwbse_env%n_RI))
      CALL build_3c_integral_block(gwbse_env%int_3c_array, qs_env, potential_parameter=bs_env%ri_metric,&
                              basis_j=bs_env%basis_set_AO, basis_k=bs_env%basis_set_AO,&
                              basis_i=bs_env%basis_set_RI,&
                              j_bf_start_from_atom=bs_env%i_ao_start_from_atom,&
                              k_bf_start_from_atom=bs_env%i_ao_start_from_atom,&
                              i_bf_start_from_atom=bs_env%i_RI_start_from_atom)
   END SUBROUTINE create_hartree_ri_workspace
   SUBROUTINE release_hartree_ri_workspace(gwbse_env)
      TYPE(gwbse_env_type)                              :: gwbse_env
      CHARACTER(len=*), PARAMETER                       :: routineN = "release_hartree_ri_workspace"

      DEALLOCATE(gwbse_env%int_3c_array)

      CALL dbcsr_release(gwbse_env%rho_dbcsr)
      DEALLOCATE(gwbse_env%rho_dbcsr)

      CALL dbcsr_release(gwbse_env%v_dbcsr%matrix)
      DEALLOCATE(gwbse_env%v_dbcsr%matrix)

   END SUBROUTINE release_hartree_ri_workspace
! **************************************************************************************************
!> \brief Allocates the workspaces for self-energy determination routine
!> \param gwbse_env Structure for holding information and workspace structures 
!> \param qs_env Quickstep environment - entry point of calculation
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE create_sigma_workspace(gwbse_env, qs_env)
      TYPE(gwbse_env_type)                               :: gwbse_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "create_sigma_workspace"
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CALL get_qs_env(qs_env, bs_env=bs_env)

      ! t_3c_w
      NULLIFY(gwbse_env%t_3c_w)
      ALLOCATE(gwbse_env%t_3c_w)
      CALL dbt_create(bs_env%t_RI__AO_AO, gwbse_env%t_3c_w)
      ! TODO : Provide option/decision whether to store the 3c integrals precomputed
      CALL compute_3c_integrals(qs_env, bs_env, gwbse_env%t_3c_w)
      ! t_3c_work_RI_AO__AO
      NULLIFY(gwbse_env%t_3c_work_RI_AO__AO)
      ALLOCATE(gwbse_env%t_3c_work_RI_AO__AO)
      CALL dbt_create(bs_env%t_RI_AO__AO, gwbse_env%t_3c_work_RI_AO__AO)
      ! t_3c_work2_RI_AO__AO
      NULLIFY(gwbse_env%t_3c_work2_RI_AO__AO)
      ALLOCATE(gwbse_env%t_3c_work2_RI_AO__AO)
      CALL dbt_create(bs_env%t_RI_AO__AO, gwbse_env%t_3c_work2_RI_AO__AO)
      ! t_W
      NULLIFY(gwbse_env%screened_dbt)
      ALLOCATE(gwbse_env%screened_dbt)
      ! Populate screened_dbt from gw run
      ! TODO : Avoid without G0W0 run?
      NULLIFY(gwbse_env%w_dbcsr%matrix)
      ALLOCATE(gwbse_env%w_dbcsr%matrix)
      CALL dbcsr_create(gwbse_env%w_dbcsr%matrix, name="W", template=bs_env%mat_RI_RI%matrix)
      CALL dbt_create(gwbse_env%w_dbcsr%matrix, gwbse_env%screened_dbt)
      ! sigma_dbt
      NULLIFY(gwbse_env%sigma_dbt)
      ALLOCATE(gwbse_env%sigma_dbt)
      CALL dbt_create(bs_env%mat_ao_ao%matrix, gwbse_env%sigma_dbt)
      ! greens_dbt
      NULLIFY(gwbse_env%greens_dbt)
      ALLOCATE(gwbse_env%greens_dbt)
      CALL dbt_create(bs_env%mat_ao_ao%matrix, gwbse_env%greens_dbt)
   END SUBROUTINE
! **************************************************************************************************
!> \brief Releases the workspaces for self-energy determination
!> \param gwbse_env 
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE release_sigma_workspace(gwbse_env)
      TYPE(gwbse_env_type)                               :: gwbse_env
      CHARACTER(len=*), PARAMETER                        :: routineN = "release_sigma_workspace"
      
      IF (ASSOCIATED(gwbse_env%t_3c_w)) THEN 
         CALL dbt_destroy(gwbse_env%t_3c_w)
         DEALLOCATE(gwbse_env%t_3c_w)
      END IF
      IF (ASSOCIATED(gwbse_env%t_3c_work_RI_AO__AO)) THEN 
         CALL dbt_destroy(gwbse_env%t_3c_work_RI_AO__AO)
         DEALLOCATE(gwbse_env%t_3c_work_RI_AO__AO)
      END IF
      IF (ASSOCIATED(gwbse_env%t_3c_work2_RI_AO__AO)) THEN 
         CALL dbt_destroy(gwbse_env%t_3c_work2_RI_AO__AO)
         DEALLOCATE(gwbse_env%t_3c_work2_RI_AO__AO)
      END IF
      IF (ASSOCIATED(gwbse_env%screened_dbt)) THEN
         CALL dbt_destroy(gwbse_env%screened_dbt)
         DEALLOCATE(gwbse_env%screened_dbt)
      END IF
      IF (ASSOCIATED(gwbse_env%sigma_dbt)) THEN
         CALL dbt_destroy(gwbse_env%sigma_dbt)
         DEALLOCATE(gwbse_env%sigma_dbt)
      END IF
      IF (ASSOCIATED(gwbse_env%greens_dbt)) THEN
         CALL dbt_destroy(gwbse_env%greens_dbt)
         DEALLOCATE(gwbse_env%greens_dbt)
      END IF
      IF (ASSOCIATED(gwbse_env%w_dbcsr%matrix)) THEN
         CALL dbcsr_release(gwbse_env%w_dbcsr%matrix)
         DEALLOCATE(gwbse_env%w_dbcsr%matrix)
      END IF
   END SUBROUTINE
   SUBROUTINE multiply_fm_cfm(trans_r, trans_c, na, nb, nc,&
                   alpha, matrix_r, matrix_c, beta, res)
      ! Transposition
      CHARACTER(len=1)                                   :: trans_r, trans_c
      INTEGER                                            :: na, nb, nc
      ! accept real numbers
      ! TODO : complex number support via interface?
      REAL(kind=dp)                                      :: alpha, beta
      TYPE(cp_fm_type)                                   :: matrix_r
      TYPE(cp_cfm_type)                                  :: matrix_c, res
      CHARACTER(len=*), PARAMETER                        :: routineN="multiply_fm_cfm"
      TYPE(cp_fm_type)                                   :: work_re, work_im, res_re, res_im
      REAL(kind=dp)                                      :: i_unit
      CHARACTER(len=1)                                   :: trans_cr

      CALL cp_fm_create(work_re, matrix_c%matrix_struct)
      CALL cp_fm_create(work_im, matrix_c%matrix_struct)
      CALL cp_fm_create(res_re, res%matrix_struct)
      CALL cp_fm_create(res_im, res%matrix_struct)
      CALL cp_cfm_to_fm(matrix_c, work_re, work_im)
      SELECT CASE (trans_c)
         CASE ("C")
            i_unit = -1.0_dp
            trans_cr = "T"
         CASE ("T")
            i_unit = 1.0_dp
            trans_cr = "T"
         CASE default
            i_unit = 1.0_dp
            trans_cr = "N"
      END SELECT
      ! Actual multiplication
      ! TODO : Correct character for trans_c
      CALL parallel_gemm(trans_r, trans_cr, na, nb, nc, &
              alpha, matrix_r, work_re, beta, res_re)
      CALL parallel_gemm(trans_r, trans_cr, na, nb, nc, &
              i_unit*alpha, matrix_r, work_im, beta, res_im)
      CALL cp_fm_to_cfm(res_re, res_im, res)
      CALL cp_fm_release(work_re)
      CALL cp_fm_release(work_im)
      CALL cp_fm_release(res_re)
      CALL cp_fm_release(res_im)

   END SUBROUTINE multiply_fm_cfm
   SUBROUTINE multiply_cfm_fm(trans_c, trans_r, na, nb, nc,&
                   alpha, matrix_c, matrix_r, beta, res)
      ! Transposition
      CHARACTER(len=1)                                   :: trans_c, trans_r
      INTEGER                                            :: na, nb, nc
      ! accept real numbers
      ! TODO : complex number support via interface?
      REAL(kind=dp)                                      :: alpha, beta
      TYPE(cp_cfm_type)                                  :: matrix_c, res
      TYPE(cp_fm_type)                                   :: matrix_r
      CHARACTER(len=*), PARAMETER                        :: routineN="multiply_fm_cfm"
      TYPE(cp_fm_type)                                   :: work_re, work_im, res_re, res_im
      REAL(kind=dp)                                      :: i_unit
      CHARACTER(len=1)                                   :: trans_cr

      CALL cp_fm_create(work_re, matrix_c%matrix_struct)
      CALL cp_fm_create(work_im, matrix_c%matrix_struct)
      CALL cp_fm_create(res_re, res%matrix_struct)
      CALL cp_fm_create(res_im, res%matrix_struct)
      CALL cp_cfm_to_fm(matrix_c, work_re, work_im)
      SELECT CASE (trans_c)
         CASE ("C")
            i_unit = -1.0_dp
            trans_cr = "T"
         CASE ("T")
            i_unit = 1.0_dp
            trans_cr = "T"
         CASE default
            i_unit = 1.0_dp
            trans_cr = "N"
      END SELECT
      ! Actual multiplication
      ! TODO : Correct character for trans_c
      CALL parallel_gemm(trans_cr, trans_r, na, nb, nc, &
              alpha, work_re, matrix_r, beta, res_re)
      CALL parallel_gemm(trans_cr, trans_r, na, nb, nc, &
              i_unit*alpha, work_im, matrix_r, beta, res_im)
      CALL cp_fm_to_cfm(res_re, res_im, res)
      CALL cp_fm_release(work_re)
      CALL cp_fm_release(work_im)
      CALL cp_fm_release(res_re)
      CALL cp_fm_release(res_im)

   END SUBROUTINE multiply_cfm_fm
END MODULE

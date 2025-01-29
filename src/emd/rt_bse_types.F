!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Data storage and other types for propagation via RT-BSE method.
!> \author Stepan Marek (01.24)
! **************************************************************************************************

MODULE rt_bse_types

   USE kinds, ONLY: dp
   USE cp_fm_types, ONLY: cp_fm_type, &
                          cp_fm_p_type, &
                          cp_fm_release, &
                          cp_fm_create, &
                          cp_fm_set_all, &
                          cp_fm_write_formatted, &
                          cp_fm_to_fm
   USE cp_cfm_types, ONLY: cp_cfm_p_type, &
                           cp_cfm_type, &
                           cp_cfm_set_all, &
                           cp_cfm_create, &
                           cp_fm_to_cfm, &
                           cp_cfm_to_fm, &
                           cp_cfm_to_cfm, &
                           cp_cfm_release
   USE cp_fm_basic_linalg, ONLY: cp_fm_invert, &
                                 cp_fm_transpose, &
                                 cp_fm_column_scale, &
                                 cp_fm_scale_and_add
   USE cp_cfm_basic_linalg, ONLY: cp_cfm_column_scale
   USE cp_dbcsr_api, ONLY: dbcsr_type, &
                           dbcsr_p_type, &
                           dbcsr_create, &
                           dbcsr_release, &
                           dbcsr_print, &
                           dbcsr_get_info, &
                           dbcsr_copy, &
                           dbcsr_set, &
                           dbcsr_scale
   USE parallel_gemm_api, ONLY: parallel_gemm
   USE dbt_api, ONLY: dbt_type, &
                      dbt_pgrid_type, &
                      dbt_pgrid_create, &
                      dbt_pgrid_destroy, &
                      dbt_mp_environ_pgrid, &
                      dbt_default_distvec, &
                      dbt_distribution_type, &
                      dbt_distribution_new, &
                      dbt_distribution_destroy, &
                      dbt_create, &
                      dbt_copy_matrix_to_tensor, &
                      dbt_get_num_blocks, &
                      dbt_destroy
   USE cp_dbcsr_operations, ONLY: copy_dbcsr_to_fm, &
                                  copy_fm_to_dbcsr
   USE qs_mo_types, ONLY: mo_set_type
   USE basis_set_types, ONLY: gto_basis_set_p_type
   USE basis_set_types, ONLY: gto_basis_set_p_type
   USE cp_control_types, ONLY: dft_control_type
   USE qs_environment_types, ONLY: qs_environment_type, &
                                   get_qs_env
   USE force_env_types, ONLY: force_env_type
   USE post_scf_bandstructure_types, ONLY: post_scf_bandstructure_type
   USE post_scf_bandstructure_utils, ONLY: create_and_init_bs_env
   USE rt_propagation_types, ONLY: rt_prop_type, &
                                   get_rtp
   USE rt_propagation_utils, ONLY: warn_section_unused
   USE rt_propagator_init, ONLY: rt_initialize_rho_from_mos
   USE rt_propagation_methods, ONLY: s_matrices_create
   USE qs_moments, ONLY: build_local_moment_matrix
   USE moments_utils, ONLY: get_reference_point
   USE matrix_exp, ONLY: get_nsquare_norder
   USE gw_integrals, ONLY: build_3c_integral_block
   USE gw_large_cell_Gamma, ONLY: compute_3c_integrals
   USE qs_tensors, ONLY: neighbor_list_3c_destroy
!    USE gw_utils,                        ONLY: create_and_init_bs_env_for_gw,&
!                                               setup_AO_and_RI_basis_set,&
!                                               get_RI_basis_and_basis_function_indices,&
!                                               set_heuristic_parameters,&
!                                               set_parallelization_parameters,&
!                                               allocate_and_fill_matrices_and_arrays,&
!                                               create_tensors
   USE libint_wrapper, ONLY: cp_libint_static_init
   USE libint_2c_3c, ONLY: libint_potential_type
   USE input_constants, ONLY: use_mom_ref_coac, &
                              use_mom_ref_zero, &
                              use_mom_ref_user, &
                              use_mom_ref_com, &
                              rtp_bse_ham_g0w0, &
                              rtp_bse_ham_ks, &
                              do_taylor, &
                              do_bch, &
                              do_exact
   USE physcon, ONLY: angstrom
   USE mathconstants, ONLY: z_zero
   USE input_section_types, ONLY: section_vals_type, &
                                  section_vals_val_get, &
                                  section_vals_get_subs_vals

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_bse"

   #:include "rt_bse_macros.fypp"

   PUBLIC :: rtbse_env_type, &
             create_rtbse_env, &
             release_rtbse_env, &
             multiply_cfm_fm, &
             multiply_fm_cfm, &
             create_hartree_ri_3c, &
             create_sigma_workspace_qs_only

   ! Created so that we can have an array of pointers to arrays
   TYPE series_real_type
      REAL(kind=dp), DIMENSION(:), POINTER                      :: series => NULL()
   END TYPE series_real_type
   TYPE series_complex_type
      COMPLEX(kind=dp), DIMENSION(:), POINTER                   :: series => NULL()
   END TYPE series_complex_type

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
!> \param greens_dbt Tensor for greens function/density matrix
!> \param t_3c_w Tensor containing 3c integrals
!> \param t_3c_work_RI_AO__AO Tensor sigma contraction
!> \param t_3c_work2_RI_AO__AO Tensor sigma contraction
!> \param sigma_SEX Screened exchange self-energy
!> \param sigma_COH Coulomb hole self-energy
!> \param hartree_curr Current Hartree matrix
!> \param etrs_max_iter Maximum number of ETRS iterations
!> \param ham_reference_type Which Hamiltonian to use as single particle basis
!> \param mat_exp_method Which method to use for matrix exponentiation
!> \param unit_nr Number of output unit
!> \param int_3c_array Array containing the local 3c integrals
!> \author Stepan Marek (01.24)
! **************************************************************************************************
   TYPE rtbse_env_type
      INTEGER                                                   :: n_spin = 1, &
                                                                   n_ao = -1, &
                                                                   n_RI = -1
      INTEGER, DIMENSION(2)                                     :: n_occ = -1
      REAL(kind=dp)                                             :: spin_degeneracy = 2
      REAL(kind=dp), DIMENSION(3)                               :: field = 0.0_dp
      TYPE(cp_fm_type), DIMENSION(:), POINTER                   :: moments => NULL(), &
                                                                   moments_field => NULL()
      INTEGER                                                   :: sim_step = 0, &
                                                                   sim_start = 0, &
                                                                   ! Needed for continuation runs for loading of previous moments trace
                                                                   sim_start_orig = 0, &
                                                                   sim_nsteps = -1
      REAL(kind=dp)                                             :: sim_time = 0.0_dp, &
                                                                   sim_dt = 0.1_dp, &
                                                                   etrs_threshold = 1.0e-7_dp, &
                                                                   exp_accuracy = 1.0e-10_dp, &
                                                                   ft_damping = 0.0_dp, &
                                                                   ft_start = 0.0_dp
      ! Which element of polarizability to print out
      INTEGER, DIMENSION(2)                                     :: pol_element = [1, 1]
      TYPE(dft_control_type), POINTER                           :: dft_control => NULL()
      ! DEBUG : Trying keeping the reference to previous environments inside this one
      TYPE(qs_environment_type), POINTER                        :: qs_env => NULL()
      TYPE(post_scf_bandstructure_type), POINTER                :: bs_env => NULL()
      ! Stores data needed for reading/writing to the restart files
      TYPE(section_vals_type), POINTER                          :: restart_section => NULL(), &
                                                                   field_section => NULL(), &
                                                                   rho_section => NULL(), &
                                                                   ft_section => NULL(), &
                                                                   pol_section => NULL(), &
                                                                   moments_section => NULL()
      LOGICAL                                                   :: restart_extracted = .FALSE.

      ! Different indices signify different spins
      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: ham_effective => NULL()
      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: ham_reference => NULL()
      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: ham_workspace => NULL()
      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: sigma_SEX => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER                   :: sigma_COH => NULL(), &
                                                                   hartree_curr => NULL()

      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: rho => NULL(), &
                                                                   rho_new => NULL(), &
                                                                   rho_new_last => NULL(), &
                                                                   rho_M => NULL(), &
                                                                   rho_orig => NULL()
      TYPE(cp_fm_type)                                          :: S_inv_fm = cp_fm_type(), &
                                                                   S_fm = cp_fm_type()
      ! Many routines require overlap in the complex format
      TYPE(cp_cfm_type)                                         :: S_cfm = cp_cfm_type()
      TYPE(dbcsr_type)                                          :: rho_dbcsr = dbcsr_type(), &
                                                                   v_ao_dbcsr = dbcsr_type()
      ! Indices only correspond to different workspaces
      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: rho_workspace => NULL()
      ! Many methods use real and imaginary parts separately - prevent unnecessary reallocation
      TYPE(cp_fm_type), DIMENSION(:), POINTER                   :: real_workspace => NULL()
      ! Workspace required for exact matrix exponentiation
      REAL(kind=dp), DIMENSION(:), POINTER                      :: real_eigvals => NULL()
      COMPLEX(kind=dp), DIMENSION(:), POINTER                   :: exp_eigvals => NULL()
      ! Workspace for saving the values for FT
      TYPE(series_complex_type), DIMENSION(3)                   :: moments_trace = series_complex_type()
      TYPE(series_real_type)                                    :: time_trace = series_real_type()
      TYPE(series_real_type), DIMENSION(3)                      :: field_trace = series_real_type()
      ! Workspace required for hartree_pw
      TYPE(dbcsr_type)                                          :: v_dbcsr = dbcsr_type(), &
                                                                   w_dbcsr = dbcsr_type()
#if defined(FTN_NO_DEFAULT_INIT)
      TYPE(dbt_type)                                            :: screened_dbt, &
                                                                   greens_dbt, &
                                                                   t_3c_w, &
                                                                   t_3c_work_RI_AO__AO, &
                                                                   t_3c_work2_RI_AO__AO
#else
      TYPE(dbt_type)                                            :: screened_dbt = dbt_type(), &
                                                                   greens_dbt = dbt_type(), &
                                                                   t_3c_w = dbt_type(), &
                                                                   t_3c_work_RI_AO__AO = dbt_type(), &
                                                                   t_3c_work2_RI_AO__AO = dbt_type()
#endif
      ! These matrices are always real
      INTEGER                                                   :: etrs_max_iter = 10
      INTEGER                                                   :: ham_reference_type = 2
      INTEGER                                                   :: mat_exp_method = 4
      INTEGER                                                   :: unit_nr = -1
      REAL(kind=dp), DIMENSION(:, :, :), POINTER                :: int_3c_array => NULL()

   END TYPE rtbse_env_type

CONTAINS

! **************************************************************************************************
!> \brief Allocates structures and prepares rtbse_env for run
!> \param rtbse_env rtbse_env_type that is initialised
!> \param qs_env Entry point of the calculation
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE create_rtbse_env(rtbse_env, qs_env, force_env)
      TYPE(rtbse_env_type), POINTER                             :: rtbse_env
      TYPE(qs_environment_type), POINTER                        :: qs_env
      TYPE(force_env_type), POINTER                             :: force_env
      TYPE(post_scf_bandstructure_type), POINTER                :: bs_env
      TYPE(rt_prop_type), POINTER                               :: rtp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER                 :: matrix_s
      TYPE(mo_set_type), DIMENSION(:), POINTER                  :: mos
      INTEGER                                                   :: i, k
      TYPE(section_vals_type), POINTER                          :: input, bs_sec, md_sec
      INTEGER, DIMENSION(:), POINTER                            :: pol_tmp

      ! Allocate the storage for the gwbse environment
      NULLIFY (rtbse_env)
      ALLOCATE (rtbse_env)
      ! Extract the other types first
      CALL get_qs_env(qs_env, &
                      bs_env=bs_env, &
                      rtp=rtp, &
                      matrix_s=matrix_s, &
                      mos=mos, &
                      dft_control=rtbse_env%dft_control, &
                      input=input)
      bs_sec => section_vals_get_subs_vals(input, "PROPERTIES%BANDSTRUCTURE")
      IF (.NOT. ASSOCIATED(bs_env)) THEN
         CPABORT("Cannot run RT-BSE without running GW calculation (PROPERTIES) before")
      END IF
      ! Number of spins
      rtbse_env%n_spin = bs_env%n_spin
      ! Number of atomic orbitals
      rtbse_env%n_ao = bs_env%n_ao
      ! Number of auxiliary basis orbitals
      rtbse_env%n_RI = bs_env%n_RI
      ! Number of occupied orbitals - for closed shell equals to half the number of electrons
      rtbse_env%n_occ(:) = bs_env%n_occ(:)
      ! Spin degeneracy - number of spins per orbital
      rtbse_env%spin_degeneracy = bs_env%spin_degeneracy
      ! Default field is zero
      rtbse_env%field(:) = 0.0_dp
      ! Default time is zero
      rtbse_env%sim_step = 0
      rtbse_env%sim_time = 0
      ! Time step is taken from rtp
      md_sec => section_vals_get_subs_vals(force_env%root_section, "MOTION%MD")
      CALL section_vals_val_get(md_sec, "TIMESTEP", r_val=rtbse_env%sim_dt)
      ! rtbse_env%sim_dt = rtp%dt
      ! Threshold for etrs is taken from the eps_energy from RT propagation
      rtbse_env%etrs_threshold = rtbse_env%dft_control%rtp_control%eps_ener
      rtbse_env%exp_accuracy = rtbse_env%dft_control%rtp_control%eps_exp
      ! Recover custom options
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%RTBSE%RTBSE_HAMILTONIAN", &
                                i_val=rtbse_env%ham_reference_type)
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%MAX_ITER", &
                                i_val=rtbse_env%etrs_max_iter)
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%MAT_EXP", &
                                i_val=rtbse_env%mat_exp_method)
      ! Output unit number, recovered from the post_scf_bandstructure_type
      rtbse_env%unit_nr = bs_env%unit_nr
      ! Sim start index and total number of steps as well
      CALL section_vals_val_get(md_sec, "STEP_START_VAL", i_val=rtbse_env%sim_start)
      ! Copy this value to sim_start_orig for continuation runs
      rtbse_env%sim_start_orig = rtbse_env%sim_start
      CALL section_vals_val_get(md_sec, "STEPS", i_val=rtbse_env%sim_nsteps)
      ! Get the values for the FT
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%PRINT%MOMENTS_FT%DAMPING", &
                                r_val=rtbse_env%ft_damping)
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%PRINT%MOMENTS_FT%START_TIME", &
                                r_val=rtbse_env%ft_start)
      ! TODO : Units for starting time and damping - so far give atomic units
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%PRINT%POLARIZABILITY%ELEMENT", &
                                i_vals=pol_tmp)
      IF (SIZE(pol_tmp) < 2) CPABORT("Less than two elements provided for polarizability")
      rtbse_env%pol_element(:) = pol_tmp(:)
      ! Do basic sanity checks for pol_element
      DO k = 1, 2
         IF (rtbse_env%pol_element(k) > 3 .OR. rtbse_env%pol_element(k) < 1) &
            CPABORT("Polarisation tensor element not 1,2 or 3 in at least one index")
      END DO

      ! Get the restart section
      rtbse_env%restart_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION%PRINT%RESTART")
      rtbse_env%restart_extracted = .FALSE.
      rtbse_env%field_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION%PRINT%FIELD")
      rtbse_env%moments_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION%PRINT%MOMENTS")
      rtbse_env%rho_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION%PRINT%DENSITY_MATRIX")
      rtbse_env%ft_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION%PRINT%MOMENTS_FT")
      rtbse_env%pol_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION%PRINT%POLARIZABILITY")
      ! Warn the user about print sections which are not yet implemented in the RTBSE run
      CALL warn_section_unused(input, "DFT%REAL_TIME_PROPAGATION%PRINT%CURRENT", &
                               "CURRENT print section not yet implemented for RTBSE.")
      CALL warn_section_unused(input, "DFT%REAL_TIME_PROPAGATION%PRINT%E_CONSTITUENTS", &
                               "E_CONSTITUENTS print section not yet implemented for RTBSE.")
      CALL warn_section_unused(input, "DFT%REAL_TIME_PROPAGATION%PRINT%PROGRAM_RUN_INFO", &
                               "PROGRAM_RUN_INFO print section not yet implemented for RTBSE.")
      CALL warn_section_unused(input, "DFT%REAL_TIME_PROPAGATION%PRINT%PROJECTION_MO", &
                               "PROJECTION_MO print section not yet implemented for RTBSE.")
      CALL warn_section_unused(input, "DFT%REAL_TIME_PROPAGATION%PRINT%RESTART_HISTORY", &
                               "RESTART_HISTORY print section not yet implemented for RTBSE.")
      ! DEBUG : References to previous environments
      rtbse_env%qs_env => qs_env
      rtbse_env%bs_env => bs_env

      ! Allocate moments matrices
      NULLIFY (rtbse_env%moments)
      ALLOCATE (rtbse_env%moments(3))
      NULLIFY (rtbse_env%moments_field)
      ALLOCATE (rtbse_env%moments_field(3))
      DO k = 1, 3
         ! Matrices are created from overlap template
         ! Values are initialized in initialize_rtbse_env
         CALL cp_fm_create(rtbse_env%moments(k), bs_env%fm_s_Gamma%matrix_struct)
         CALL cp_fm_create(rtbse_env%moments_field(k), bs_env%fm_s_Gamma%matrix_struct)
      END DO

      ! Allocate space for density propagation and other operations
      NULLIFY (rtbse_env%rho_workspace)
      ALLOCATE (rtbse_env%rho_workspace(4))
      DO i = 1, SIZE(rtbse_env%rho_workspace)
         CALL cp_cfm_create(rtbse_env%rho_workspace(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_cfm_set_all(rtbse_env%rho_workspace(i), CMPLX(0.0, 0.0, kind=dp))
      END DO
      ! Allocate real workspace
      NULLIFY (rtbse_env%real_workspace)
      SELECT CASE (rtbse_env%mat_exp_method)
      CASE (do_exact)
         ALLOCATE (rtbse_env%real_workspace(4))
      CASE (do_bch)
         ALLOCATE (rtbse_env%real_workspace(2))
      CASE DEFAULT
         CPABORT("Only exact and BCH matrix propagation implemented in RT-BSE")
      END SELECT
      DO i = 1, SIZE(rtbse_env%real_workspace)
         CALL cp_fm_create(rtbse_env%real_workspace(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_set_all(rtbse_env%real_workspace(i), 0.0_dp)
      END DO
      ! Allocate density matrix
      NULLIFY (rtbse_env%rho)
      ALLOCATE (rtbse_env%rho(rtbse_env%n_spin))
      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_create(rtbse_env%rho(i), matrix_struct=bs_env%fm_s_Gamma%matrix_struct)
      END DO
      ! Create the inverse overlap matrix, for use in density propagation
      ! Start by creating the actual overlap matrix
      CALL cp_fm_create(rtbse_env%S_fm, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(rtbse_env%S_inv_fm, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_cfm_create(rtbse_env%S_cfm, bs_env%fm_s_Gamma%matrix_struct)

      ! Create the single particle hamiltonian
      ! Allocate workspace
      NULLIFY (rtbse_env%ham_workspace)
      ALLOCATE (rtbse_env%ham_workspace(rtbse_env%n_spin))
      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_create(rtbse_env%ham_workspace(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_cfm_set_all(rtbse_env%ham_workspace(i), CMPLX(0.0, 0.0, kind=dp))
      END DO
      ! Now onto the Hamiltonian itself
      NULLIFY (rtbse_env%ham_reference)
      ALLOCATE (rtbse_env%ham_reference(rtbse_env%n_spin))
      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_create(rtbse_env%ham_reference(i), bs_env%fm_ks_Gamma(i)%matrix_struct)
      END DO

      ! Create the matrices and workspaces for ETRS propagation
      NULLIFY (rtbse_env%ham_effective)
      NULLIFY (rtbse_env%rho_new)
      NULLIFY (rtbse_env%rho_new_last)
      NULLIFY (rtbse_env%rho_M)
      NULLIFY (rtbse_env%rho_orig)
      ALLOCATE (rtbse_env%ham_effective(rtbse_env%n_spin))
      ALLOCATE (rtbse_env%rho_new(rtbse_env%n_spin))
      ALLOCATE (rtbse_env%rho_new_last(rtbse_env%n_spin))
      ALLOCATE (rtbse_env%rho_M(rtbse_env%n_spin))
      ALLOCATE (rtbse_env%rho_orig(rtbse_env%n_spin))
      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_create(rtbse_env%ham_effective(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_cfm_set_all(rtbse_env%ham_effective(i), CMPLX(0.0, 0.0, kind=dp))
         CALL cp_cfm_create(rtbse_env%rho_new(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_cfm_set_all(rtbse_env%rho_new(i), CMPLX(0.0, 0.0, kind=dp))
         CALL cp_cfm_create(rtbse_env%rho_new_last(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_cfm_set_all(rtbse_env%rho_new_last(i), CMPLX(0.0, 0.0, kind=dp))
         CALL cp_cfm_create(rtbse_env%rho_M(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_cfm_set_all(rtbse_env%rho_M(i), CMPLX(0.0, 0.0, kind=dp))
         CALL cp_cfm_create(rtbse_env%rho_orig(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
      END DO

      ! Fields for exact diagonalisation
      NULLIFY (rtbse_env%real_eigvals)
      ALLOCATE (rtbse_env%real_eigvals(rtbse_env%n_ao))
      rtbse_env%real_eigvals(:) = 0.0_dp
      NULLIFY (rtbse_env%exp_eigvals)
      ALLOCATE (rtbse_env%exp_eigvals(rtbse_env%n_ao))
      rtbse_env%exp_eigvals(:) = CMPLX(0.0, 0.0, kind=dp)

      ! Workspace for FT - includes (in principle) the zeroth step and the extra last step
      DO i = 1, 3
         NULLIFY (rtbse_env%moments_trace(i)%series)
         ALLOCATE (rtbse_env%moments_trace(i)%series(rtbse_env%sim_nsteps + 2), source=z_zero)
         NULLIFY (rtbse_env%field_trace(i)%series)
         ALLOCATE (rtbse_env%field_trace(i)%series(rtbse_env%sim_nsteps + 2), source=0.0_dp)
      END DO
      NULLIFY (rtbse_env%time_trace%series)
      ALLOCATE (rtbse_env%time_trace%series(rtbse_env%sim_nsteps + 2), source=0.0_dp)

      ! Allocate self-energy parts and dynamic Hartree potential
      NULLIFY (rtbse_env%hartree_curr)
      NULLIFY (rtbse_env%sigma_SEX)
      NULLIFY (rtbse_env%sigma_COH)
      ALLOCATE (rtbse_env%hartree_curr(rtbse_env%n_spin))
      ALLOCATE (rtbse_env%sigma_SEX(rtbse_env%n_spin))
      ALLOCATE (rtbse_env%sigma_COH(rtbse_env%n_spin))
      DO i = 1, rtbse_env%n_spin
         CALL cp_fm_create(rtbse_env%sigma_COH(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_cfm_create(rtbse_env%sigma_SEX(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_create(rtbse_env%hartree_curr(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_set_all(rtbse_env%sigma_COH(i), 0.0_dp)
         CALL cp_cfm_set_all(rtbse_env%sigma_SEX(i), CMPLX(0.0, 0.0, kind=dp))
         CALL cp_fm_set_all(rtbse_env%hartree_curr(i), 0.0_dp)
      END DO

      ! Allocate workspaces for get_sigma
      CALL create_sigma_workspace(rtbse_env, qs_env)

      ! Depending on the chosen methods, allocate extra workspace
      CALL create_hartree_ri_workspace(rtbse_env, qs_env)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Simple reimplementation of cp_fm_release_pp1 for complex matrices
!> \param matrices cp_cfm_p_type(:)
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE cp_cfm_release_pa1(matrices)
      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: matrices
      INTEGER                                                   :: i

      DO i = 1, SIZE(matrices)
         CALL cp_cfm_release(matrices(i))
      END DO
      DEALLOCATE (matrices)
      NULLIFY (matrices)
   END SUBROUTINE cp_cfm_release_pa1

! **************************************************************************************************
!> \brief Releases the environment allocated structures
!> \param rtbse_env
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE release_rtbse_env(rtbse_env)
      TYPE(rtbse_env_type), POINTER                             :: rtbse_env
      INTEGER                                                   :: i

      CALL cp_cfm_release_pa1(rtbse_env%ham_effective)
      CALL cp_cfm_release_pa1(rtbse_env%ham_workspace)
      CALL cp_fm_release(rtbse_env%sigma_COH)
      CALL cp_cfm_release_pa1(rtbse_env%sigma_SEX)
      CALL cp_fm_release(rtbse_env%hartree_curr)
      CALL cp_cfm_release_pa1(rtbse_env%ham_reference)
      CALL cp_cfm_release_pa1(rtbse_env%rho)
      CALL cp_cfm_release_pa1(rtbse_env%rho_workspace)
      CALL cp_cfm_release_pa1(rtbse_env%rho_new)
      CALL cp_cfm_release_pa1(rtbse_env%rho_new_last)
      CALL cp_cfm_release_pa1(rtbse_env%rho_M)
      CALL cp_cfm_release_pa1(rtbse_env%rho_orig)
      CALL cp_fm_release(rtbse_env%real_workspace)
      CALL cp_fm_release(rtbse_env%S_inv_fm)
      CALL cp_fm_release(rtbse_env%S_fm)
      CALL cp_cfm_release(rtbse_env%S_cfm)

      ! DO i = 1, 3
      !    CALL cp_fm_release(rtbse_env%moments(i)%matrix)
      !    CALL cp_fm_release(rtbse_env%moments_field(i)%matrix)
      ! END DO
      CALL cp_fm_release(rtbse_env%moments)
      CALL cp_fm_release(rtbse_env%moments_field)

      CALL release_sigma_workspace(rtbse_env)

      CALL release_hartree_ri_workspace(rtbse_env)

      DEALLOCATE (rtbse_env%real_eigvals)
      DEALLOCATE (rtbse_env%exp_eigvals)
      DO i = 1, 3
         DEALLOCATE (rtbse_env%moments_trace(i)%series)
         DEALLOCATE (rtbse_env%field_trace(i)%series)
      END DO
      DEALLOCATE (rtbse_env%time_trace%series)

      ! Deallocate the neighbour list that is not deallocated in gw anymore
      IF (ASSOCIATED(rtbse_env%bs_env%nl_3c%ij_list)) CALL neighbor_list_3c_destroy(rtbse_env%bs_env%nl_3c)
      ! Deallocate the storage for the environment itself
      DEALLOCATE (rtbse_env)
      ! Nullify to make sure it is not used again
      NULLIFY (rtbse_env)

   END SUBROUTINE
! **************************************************************************************************
!> \brief Allocates the workspaces for Hartree RI method
!> \note RI method calculates the Hartree contraction without the use of DBT, as it cannot emulate vectors
!> \param rtbse_env
!> \param qs_env Quickstep environment - entry point of calculation
!> \author Stepan Marek
!> \date 05.2024
! **************************************************************************************************
   SUBROUTINE create_hartree_ri_workspace(rtbse_env, qs_env)
      TYPE(rtbse_env_type)                              :: rtbse_env
      TYPE(qs_environment_type), POINTER                :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER        :: bs_env

      CALL get_qs_env(qs_env, bs_env=bs_env)

      CALL dbcsr_create(rtbse_env%rho_dbcsr, name="Sparse density", template=bs_env%mat_ao_ao%matrix)
      CALL dbcsr_create(rtbse_env%v_ao_dbcsr, name="Sparse Hartree", template=bs_env%mat_ao_ao%matrix)

      CALL create_hartree_ri_3c(rtbse_env%rho_dbcsr, rtbse_env%int_3c_array, rtbse_env%n_ao, rtbse_env%n_RI, &
                                bs_env%basis_set_AO, bs_env%basis_set_RI, bs_env%i_RI_start_from_atom, &
                                bs_env%ri_metric, qs_env, rtbse_env%unit_nr)
   END SUBROUTINE create_hartree_ri_workspace
! **************************************************************************************************
!> \brief Separated method for allocating the 3c integrals for RI Hartree
!> \note RI method calculates the Hartree contraction without the use of DBT, as it cannot emulate vectors
!> \param rho_dbcsr matrix used for the description of shape of 3c array
!> \param int_3c 3-center integral array to be allocated and filled
!> \param n_ao Number of atomic orbitals
!> \param n_RI Number of auxiliary RI orbitals
!> \param basis_set_AO AO basis set
!> \param basis_set_RI RI auxiliary basis set
!> \param i_RI_start_from_atom Array of indices where functions of a given atom in RI basis start
!> \param unit_nr Unit number used for printing information about the size of int_3c
!> \author Stepan Marek
!> \date 01.2025
! **************************************************************************************************
   SUBROUTINE create_hartree_ri_3c(rho_dbcsr, int_3c, n_ao, n_RI, basis_set_AO, basis_set_RI, &
                                   i_RI_start_from_atom, ri_metric, qs_env, unit_nr)
      TYPE(dbcsr_type)                                  :: rho_dbcsr
      REAL(kind=dp), DIMENSION(:, :, :), POINTER          :: int_3c
      INTEGER                                           :: n_ao, n_RI
      TYPE(gto_basis_set_p_type), DIMENSION(:)          :: basis_set_AO, &
                                                           basis_set_RI
      INTEGER, DIMENSION(:)                             :: i_RI_start_from_atom
      TYPE(libint_potential_type)                       :: ri_metric
      TYPE(qs_environment_type), POINTER                :: qs_env
      INTEGER                                           :: unit_nr
      REAL(kind=dp)                                     :: size_mb
      INTEGER                                           :: nblkrows_local, &
                                                           nblkcols_local, &
                                                           i_blk_local, &
                                                           j_blk_local, &
                                                           nrows_local, &
                                                           ncols_local, &
                                                           col_local_offset, &
                                                           row_local_offset, &
                                                           start_col_index, &
                                                           end_col_index, &
                                                           start_row_index, &
                                                           end_row_index
      INTEGER, DIMENSION(:), POINTER                    :: local_blk_rows, &
                                                           local_blk_cols, &
                                                           row_blk_size, &
                                                           col_blk_size
      ! TODO : Implement option/decision to not precompute all the 3c integrals
      size_mb = REAL(n_ao, kind=dp)*REAL(n_ao, kind=dp)*REAL(n_RI, kind=dp)* &
                REAL(STORAGE_SIZE(size_mb), kind=dp)/8.0_dp/1024.0_dp/1024.0_dp
      IF (unit_nr > 0) WRITE (unit_nr, '(A44,E32.2E3,A4)') &
         " RTBSE| Approximate size of the 3c integrals", size_mb, " MiB"

      ! Get the number of block rows and columns
      CALL dbcsr_get_info(rho_dbcsr, nblkrows_local=nblkrows_local, nblkcols_local=nblkcols_local)
      ! Get the global indices of local rows and columns
      CALL dbcsr_get_info(rho_dbcsr, local_rows=local_blk_rows, local_cols=local_blk_cols)
      ! Get the sizes of all blocks
      CALL dbcsr_get_info(rho_dbcsr, row_blk_size=row_blk_size, col_blk_size=col_blk_size)

      ! Get the total required local rows and cols
      nrows_local = 0
      DO i_blk_local = 1, nblkrows_local
         nrows_local = nrows_local + row_blk_size(local_blk_rows(i_blk_local))
      END DO
      ncols_local = 0
      DO j_blk_local = 1, nblkcols_local
         ncols_local = ncols_local + col_blk_size(local_blk_cols(j_blk_local))
      END DO

      ! Allocate the appropriate storage
      ALLOCATE (int_3c(nrows_local, ncols_local, n_RI))

      ! Fill the storage with appropriate values, block by block
      row_local_offset = 1
      DO i_blk_local = 1, nblkrows_local
         col_local_offset = 1
         DO j_blk_local = 1, nblkcols_local
            start_row_index = row_local_offset
            end_row_index = start_row_index + row_blk_size(local_blk_rows(i_blk_local)) - 1
            start_col_index = col_local_offset
            end_col_index = start_col_index + col_blk_size(local_blk_cols(j_blk_local)) - 1
            CALL build_3c_integral_block(int_3c(start_row_index:end_row_index, &
                                                start_col_index:end_col_index, &
                                                1:n_RI), &
                                         qs_env, potential_parameter=ri_metric, &
                                         basis_j=basis_set_AO, basis_k=basis_set_AO, &
                                         basis_i=basis_set_RI, &
                                         atom_j=local_blk_rows(i_blk_local), &
                                         atom_k=local_blk_cols(j_blk_local), &
                                         i_bf_start_from_atom=i_RI_start_from_atom)
            col_local_offset = col_local_offset + col_blk_size(local_blk_cols(j_blk_local))
         END DO
         row_local_offset = row_local_offset + row_blk_size(local_blk_rows(i_blk_local))
      END DO
   END SUBROUTINE create_hartree_ri_3c
! **************************************************************************************************
!> \brief Releases the workspace for the Hartree RI method
!> \param rtbse_env RT-BSE Environment, containing specific RI Hartree storage
!> \author Stepan Marek
!> \date 09.2024
! **************************************************************************************************
   SUBROUTINE release_hartree_ri_workspace(rtbse_env)
      TYPE(rtbse_env_type)                              :: rtbse_env

      DEALLOCATE (rtbse_env%int_3c_array)

      CALL dbcsr_release(rtbse_env%rho_dbcsr)

      CALL dbcsr_release(rtbse_env%v_dbcsr)

      CALL dbcsr_release(rtbse_env%v_ao_dbcsr)

   END SUBROUTINE release_hartree_ri_workspace
! **************************************************************************************************
!> \brief Allocates the workspaces for self-energy determination routine
!> \param rtbse_env Structure for holding information and workspace structures
!> \param qs_env Quickstep environment - entry point of calculation
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE create_sigma_workspace(rtbse_env, qs_env)
      TYPE(rtbse_env_type)                               :: rtbse_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CALL create_sigma_workspace_qs_only(qs_env, rtbse_env%screened_dbt, rtbse_env%w_dbcsr, &
                                          rtbse_env%t_3c_w, rtbse_env%t_3c_work_RI_AO__AO, &
                                          rtbse_env%t_3c_work2_RI_AO__AO, rtbse_env%greens_dbt)
   END SUBROUTINE create_sigma_workspace
! **************************************************************************************************
!> \brief Allocates the workspaces for self-energy determination routine
!> \note Does so without referencing the rtbse_env
!> \note References bs_env
!> \param rtbse_env Structure for holding information and workspace structures
!> \param qs_env Quickstep environment - entry point of calculation
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE create_sigma_workspace_qs_only(qs_env, screened_dbt, screened_dbcsr, int_3c_dbt, &
                                             work_dbt_3c_1, work_dbt_3c_2, work_dbt_2c)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type)                                   :: screened_dbcsr
      TYPE(dbt_type)                                     :: screened_dbt, &
                                                            int_3c_dbt, &
                                                            work_dbt_3c_1, &
                                                            work_dbt_3c_2, &
                                                            work_dbt_2c
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CALL get_qs_env(qs_env, bs_env=bs_env)

      ! t_3c_w
      CALL dbt_create(bs_env%t_RI__AO_AO, int_3c_dbt)
      ! TODO : Provide option/decision whether to store the 3c integrals precomputed
      CALL compute_3c_integrals(qs_env, bs_env, int_3c_dbt)
      ! t_3c_work_RI_AO__AO
      CALL dbt_create(bs_env%t_RI_AO__AO, work_dbt_3c_1)
      ! t_3c_work2_RI_AO__AO
      CALL dbt_create(bs_env%t_RI_AO__AO, work_dbt_3c_2)
      ! t_W
      ! Populate screened_dbt from gw run
      CALL dbcsr_create(screened_dbcsr, name="W", template=bs_env%mat_RI_RI%matrix)
      CALL dbt_create(screened_dbcsr, screened_dbt)
      ! greens_dbt
      CALL dbt_create(bs_env%mat_ao_ao%matrix, work_dbt_2c)
   END SUBROUTINE
! **************************************************************************************************
!> \brief Releases the workspaces for self-energy determination
!> \param rtbse_env
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE release_sigma_workspace(rtbse_env)
      TYPE(rtbse_env_type)                               :: rtbse_env

      CALL dbt_destroy(rtbse_env%t_3c_w)
      CALL dbt_destroy(rtbse_env%t_3c_work_RI_AO__AO)
      CALL dbt_destroy(rtbse_env%t_3c_work2_RI_AO__AO)
      CALL dbt_destroy(rtbse_env%screened_dbt)
      CALL dbt_destroy(rtbse_env%greens_dbt)
      CALL dbcsr_release(rtbse_env%w_dbcsr)
   END SUBROUTINE
! **************************************************************************************************
!> \brief Multiplies real matrix by a complex matrix from the right
!> \note So far only converts the real matrix to complex one, potentially doubling the work
!> \param rtbse_env
!> \author Stepan Marek
!> \date 09.2024
! **************************************************************************************************
   SUBROUTINE multiply_fm_cfm(trans_r, trans_c, na, nb, nc, &
                              alpha, matrix_r, matrix_c, beta, res)
      ! Transposition
      CHARACTER(len=1)                                   :: trans_r, trans_c
      INTEGER                                            :: na, nb, nc
      ! accept real numbers
      ! TODO : Just use complex numbers and import z_one, z_zero etc.
      REAL(kind=dp)                                      :: alpha, beta
      TYPE(cp_fm_type)                                   :: matrix_r
      TYPE(cp_cfm_type)                                  :: matrix_c, res
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
! **************************************************************************************************
!> \brief Multiplies complex matrix by a real matrix from the right
!> \note So far only converts the real matrix to complex one, potentially doubling the work
!> \param rtbse_env
!> \author Stepan Marek
!> \date 09.2024
! **************************************************************************************************
   SUBROUTINE multiply_cfm_fm(trans_c, trans_r, na, nb, nc, &
                              alpha, matrix_c, matrix_r, beta, res)
      ! Transposition
      CHARACTER(len=1)                                   :: trans_c, trans_r
      INTEGER                                            :: na, nb, nc
      ! accept real numbers
      ! TODO : complex number support via interface?
      REAL(kind=dp)                                      :: alpha, beta
      TYPE(cp_cfm_type)                                  :: matrix_c, res
      TYPE(cp_fm_type)                                   :: matrix_r
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

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
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
                          cp_fm_release, &
                          cp_fm_create, &
                          cp_fm_set_all
   USE cp_fm_struct, ONLY: cp_fm_struct_type, &
                           cp_fm_struct_create, &
                           cp_fm_struct_release
   USE cp_cfm_types, ONLY: cp_cfm_type, &
                           cp_cfm_set_all, &
                           cp_cfm_create, &
                           cp_fm_to_cfm, &
                           cp_cfm_to_fm, &
                           cp_cfm_release
   USE cp_dbcsr_api, ONLY: dbcsr_type, &
                           dbcsr_p_type, &
                           dbcsr_create, &
                           dbcsr_release, &
                           dbcsr_get_info
   USE parallel_gemm_api, ONLY: parallel_gemm
   USE dbt_api, ONLY: dbt_type, &
                      dbt_create, &
                      dbt_destroy
   USE qs_mo_types, ONLY: mo_set_type
   USE basis_set_types, ONLY: gto_basis_set_p_type
   USE cp_control_types, ONLY: dft_control_type
   USE qs_environment_types, ONLY: qs_environment_type, &
                                   get_qs_env
   USE force_env_types, ONLY: force_env_type
   USE post_scf_bandstructure_types, ONLY: post_scf_bandstructure_type
   USE rt_propagation_types, ONLY: rt_prop_type
   USE rt_propagation_utils, ONLY: warn_section_unused
   USE gw_integrals, ONLY: build_3c_integral_block
   USE gw_large_cell_Gamma, ONLY: compute_3c_integrals
   USE qs_tensors, ONLY: neighbor_list_3c_destroy
   USE libint_2c_3c, ONLY: libint_potential_type
   USE input_constants, ONLY: use_mom_ref_coac, &
                              do_bch, &
                              do_exact, &
                              rtp_bse_ham_g0w0, &
                              rtp_method_bse_linearized
   USE bse_util, ONLY: determine_cutoff_indices
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

   ! ! Created so that we can have an array of pointers to arrays
   ! TYPE series_real_type
   !    REAL(kind=dp), DIMENSION(:), POINTER                      :: series => NULL()
   ! END TYPE series_real_type
   ! TYPE series_complex_type
   !    COMPLEX(kind=dp), DIMENSION(:), POINTER                   :: series => NULL()
   ! END TYPE series_complex_type

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
!> \param ham_reference Reference Hamiltonian, which does not change in the
!>                      propagation = DFT+G0W0 - initial Hartree - initial COHSEX
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
      ! Active MO window for linearized RT-BSE truncation. When no truncation is requested,
      ! first_active_mo=1, last_active_mo=n_ao, and mo_active=n_ao. The window is the
      ! combined inclusive bound that covers both spin channels.
      INTEGER                                                   :: first_active_mo = 1, &
                                                                   last_active_mo = -1, &
                                                                   mo_active = -1
      REAL(KIND=dp)                                             :: rtbse_energy_cutoff_occ = -1.0_dp, &
                                                                   rtbse_energy_cutoff_empty = -1.0_dp
      LOGICAL                                                   :: active_mo_truncation = .FALSE.
      ! Owned fm structures sized to the active MO window. Equal to the full n_ao x n_ao
      ! when no truncation is active.
      TYPE(cp_fm_struct_type), POINTER                          :: fm_struct_mo_active => NULL()
      ! n_ao x mo_active fm struct used for C_active and AO<->MO rectangular intermediates.
      TYPE(cp_fm_struct_type), POINTER                          :: fm_struct_ao_mo_active => NULL()
      ! Truncated MO coefficient slabs C_active(:,:) of size n_ao x mo_active for each spin
      ! (only allocated for the linearized path).
      TYPE(cp_fm_type), DIMENSION(:), POINTER                   :: C_active => NULL()
      ! Rectangular n_ao x mo_active scratch used by linearized AO<->MO transforms.
      TYPE(cp_fm_type), DIMENSION(:), POINTER                   :: ao_mo_workspace => NULL()
      REAL(kind=dp)                                             :: spin_degeneracy = 2
      REAL(kind=dp), DIMENSION(3)                               :: field = 0.0_dp
      TYPE(cp_fm_type), DIMENSION(:), POINTER                   :: moments => NULL(), &
                                                                   moments_field => NULL()
      INTEGER                                                   :: sim_step = 0, &
                                                                   sim_start = 0, &
                                                                   ! Needed to continue runs by loading previous moments trace
                                                                   sim_start_orig = 0, &
                                                                   sim_nsteps = -1, &
                                                                   ! Default reference point type for output moments
                                                                   ! Field moments always use zero reference
                                                                   moment_ref_type = use_mom_ref_coac
      REAL(kind=dp), DIMENSION(:), POINTER                      :: user_moment_ref_point => NULL()
      REAL(kind=dp)                                             :: sim_time = 0.0_dp, &
                                                                   sim_dt = 0.1_dp, &
                                                                   maximum_timestep = -1.0_dp, &
                                                                   omega_max = -1.0_dp, &
                                                                   etrs_threshold = 1.0e-7_dp, &
                                                                   exp_accuracy = 1.0e-10_dp, &
                                                                   ft_damping = 0.0_dp, &
                                                                   ft_start = 0.0_dp
      ! Which element of polarizability to print out
      INTEGER, DIMENSION(:, :), POINTER                         :: pol_elements => NULL()
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
                                                                   moments_section => NULL(), &
                                                                   rtp_section => NULL()
      LOGICAL                                                   :: restart_extracted = .FALSE.

      ! Different indices signify different spins
      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: ham_effective => NULL()
      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: ham_reference => NULL()
      !Only for linearised RTBSE
      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: ham_reference_singleparticle => NULL()
      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: ham_workspace => NULL()
      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: sigma_SEX => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER                   :: sigma_COH => NULL(), &
                                                                   hartree_curr => NULL()
      ! AO-sized scratch buffers used in the linearized RT-BSE path so that the MO-sized
      ! sigma_COH/sigma_SEX/hartree_curr matrices above can be allocated on
      ! fm_struct_mo_active. Only allocated when my_linearized=.TRUE..
      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: sigma_SEX_ao => NULL()
      TYPE(cp_fm_type), DIMENSION(:), POINTER                   :: sigma_COH_ao => NULL(), &
                                                                   hartree_curr_ao => NULL()
      ! mo_active x mo_active real workspace pair used in MO-side transforms (linearized only).
      TYPE(cp_fm_type), DIMENSION(:), POINTER                   :: real_workspace_mo => NULL()

      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: rho => NULL(), &
                                                                   rho_new => NULL(), &
                                                                   rho_new_last => NULL(), &
                                                                   rho_M => NULL(), &
                                                                   rho_orig => NULL(), &
                                                                   rho_ao_lrrtbse => NULL()
      ! Workspace for rk4 in linearized RTBSE
      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: rk4_coefficients => NULL()

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
      ! AO-sized complex scratch used to stage the real part in get_sigma_complex.
      TYPE(cp_cfm_type), DIMENSION(:), POINTER                  :: sigma_complex_workspace => NULL()
      ! Workspace required for exact matrix exponentiation
      REAL(kind=dp), DIMENSION(:), POINTER                      :: real_eigvals => NULL()
      COMPLEX(kind=dp), DIMENSION(:), POINTER                   :: exp_eigvals => NULL()
      ! Workspace for saving the values for FT
      ! TODO : Change back to multi-dimensional arrays
      ! Index 1 : spin, Index 2 : direction, Index 3 : time point
      COMPLEX(kind=dp), DIMENSION(:, :, :), POINTER               :: moments_trace => NULL()
      REAL(kind=dp), DIMENSION(:), POINTER                      :: time_trace => NULL()
      ! Index 1 : direction, Index 2 : time point
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER                 :: field_trace => NULL()
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
      ! Parameters for Padé refinement
      REAL(kind=dp)                                             :: pade_e_min = 0.0_dp, &
                                                                   pade_e_max = 100.0_dp, &
                                                                   pade_e_step = 0.05_dp, &
                                                                   pade_fit_e_min = 0.0_dp, &
                                                                   pade_fit_e_max = -1.0_dp
      INTEGER                                                   :: pade_npoints = 0
      LOGICAL                                                   :: pade_requested = .FALSE.
      COMPLEX(kind=dp), DIMENSION(:), POINTER                   :: pade_x_eval => NULL()

   END TYPE rtbse_env_type

CONTAINS

! **************************************************************************************************
!> \brief Allocates structures and prepares rtbse_env for run
!> \param rtbse_env rtbse_env_type that is initialised
!> \param qs_env Entry point of the calculation
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE create_rtbse_env(rtbse_env, force_env, linearized)
      TYPE(rtbse_env_type), POINTER                             :: rtbse_env
      TYPE(force_env_type), POINTER                             :: force_env
      LOGICAL, OPTIONAL                                        :: linearized
      TYPE(post_scf_bandstructure_type), POINTER                :: bs_env
      TYPE(rt_prop_type), POINTER                               :: rtp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER                 :: matrix_s
      TYPE(mo_set_type), DIMENSION(:), POINTER                  :: mos
      INTEGER                                                   :: i, k
      TYPE(section_vals_type), POINTER                          :: input, bs_sec, md_sec
      LOGICAL :: my_linearized
      TYPE(cp_fm_struct_type), POINTER                          :: mo_struct => NULL()

      IF (PRESENT(linearized)) THEN
         my_linearized = linearized
      ELSE
         my_linearized = .FALSE.
      END IF

      ! Allocate the storage for the gwbse environment
      NULLIFY (rtbse_env)
      ALLOCATE (rtbse_env)
      ! Extract the other types first
      CALL get_qs_env(force_env%qs_env, &
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
      ! Active MO window for linearized RT-BSE truncation
      rtbse_env%rtbse_energy_cutoff_occ = rtbse_env%dft_control%rtp_control%rtbse_energy_cutoff_occ
      rtbse_env%rtbse_energy_cutoff_empty = rtbse_env%dft_control%rtp_control%rtbse_energy_cutoff_empty
      CALL determine_active_mo_window(rtbse_env, my_linearized, bs_env)
      ! Owned active-MO matrix structure (currently identical to the full n_ao x n_ao struct
      ! when no truncation is active; will be used by the linearized RT-BSE allocation path).
      NULLIFY (rtbse_env%fm_struct_mo_active)
      CALL cp_fm_struct_create(rtbse_env%fm_struct_mo_active, &
                               bs_env%fm_ks_Gamma(1)%matrix_struct%para_env, &
                               bs_env%fm_ks_Gamma(1)%matrix_struct%context, &
                               rtbse_env%mo_active, rtbse_env%mo_active)
      ! Rectangular n_ao x mo_active struct used for C_active and AO<->MO intermediates
      NULLIFY (rtbse_env%fm_struct_ao_mo_active)
      CALL cp_fm_struct_create(rtbse_env%fm_struct_ao_mo_active, &
                               bs_env%fm_ks_Gamma(1)%matrix_struct%para_env, &
                               bs_env%fm_ks_Gamma(1)%matrix_struct%context, &
                               rtbse_env%n_ao, rtbse_env%mo_active)
      ! Choose the matrix struct used for MO-side persistent matrices.
      ! Linearized RT-BSE: mo_active x mo_active. Full RT-BSE: full AO struct (unchanged).
      IF (my_linearized) THEN
         mo_struct => rtbse_env%fm_struct_mo_active
      ELSE
         mo_struct => bs_env%fm_ks_Gamma(1)%matrix_struct
      END IF
      ! Output unit number, recovered from the post_scf_bandstructure_type
      rtbse_env%unit_nr = bs_env%unit_nr
      ! Sim start index and total number of steps as well
      CALL section_vals_val_get(md_sec, "STEP_START_VAL", i_val=rtbse_env%sim_start)
      ! Copy this value to sim_start_orig for continuation runs
      rtbse_env%sim_start_orig = rtbse_env%sim_start
      CALL section_vals_val_get(md_sec, "STEPS", i_val=rtbse_env%sim_nsteps)
      ! Get the values for the FT
      rtbse_env%ft_damping = rtbse_env%dft_control%rtp_control%ft_damping
      rtbse_env%ft_damping = rtbse_env%dft_control%rtp_control%ft_t0
      rtbse_env%pol_elements => rtbse_env%dft_control%rtp_control%print_pol_elements

      rtbse_env%rtp_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION")
      ! Get the restart section
      rtbse_env%restart_section => section_vals_get_subs_vals(rtbse_env%rtp_section, "PRINT%RESTART")
      rtbse_env%restart_extracted = .FALSE.
      rtbse_env%field_section => section_vals_get_subs_vals(rtbse_env%rtp_section, "PRINT%FIELD")
      rtbse_env%moments_section => section_vals_get_subs_vals(rtbse_env%rtp_section, "PRINT%MOMENTS")
      ! Moment specification
      CALL section_vals_val_get(rtbse_env%rtp_section, "PRINT%MOMENTS%REFERENCE", &
                                i_val=rtbse_env%moment_ref_type)
      CALL section_vals_val_get(rtbse_env%rtp_section, "PRINT%MOMENTS%REFERENCE_POINT", &
                                r_vals=rtbse_env%user_moment_ref_point)
      rtbse_env%rho_section => section_vals_get_subs_vals(rtbse_env%rtp_section, "PRINT%DENSITY_MATRIX")
      rtbse_env%ft_section => section_vals_get_subs_vals(rtbse_env%rtp_section, "PRINT%MOMENTS_FT")
      rtbse_env%pol_section => section_vals_get_subs_vals(rtbse_env%rtp_section, "PRINT%POLARIZABILITY")
      ! Warn the user about print sections which are not yet implemented in the RTBSE run
      CALL warn_section_unused(rtbse_env%rtp_section, "PRINT%CURRENT", &
                               "CURRENT print section not yet implemented for RTBSE.")
      CALL warn_section_unused(rtbse_env%rtp_section, "PRINT%E_CONSTITUENTS", &
                               "E_CONSTITUENTS print section not yet implemented for RTBSE.")
      CALL warn_section_unused(rtbse_env%rtp_section, "PRINT%PROGRAM_RUN_INFO", &
                               "PROGRAM_RUN_INFO print section not yet implemented for RTBSE.")
      CALL warn_section_unused(rtbse_env%rtp_section, "PRINT%PROJECTION_MO", &
                               "PROJECTION_MO print section not yet implemented for RTBSE.")
      CALL warn_section_unused(rtbse_env%rtp_section, "PRINT%RESTART_HISTORY", &
                               "RESTART_HISTORY print section not yet implemented for RTBSE.")
      ! DEBUG : References to previous environments
      rtbse_env%qs_env => force_env%qs_env
      rtbse_env%bs_env => bs_env
      ! Padé refinement
      rtbse_env%pade_requested = rtbse_env%dft_control%rtp_control%pade_requested
      rtbse_env%pade_e_min = rtbse_env%dft_control%rtp_control%pade_e_min
      rtbse_env%pade_e_step = rtbse_env%dft_control%rtp_control%pade_e_step
      rtbse_env%pade_e_max = rtbse_env%dft_control%rtp_control%pade_e_max
      rtbse_env%pade_fit_e_min = rtbse_env%dft_control%rtp_control%pade_fit_e_min
      rtbse_env%pade_fit_e_max = rtbse_env%dft_control%rtp_control%pade_fit_e_max
      rtbse_env%pade_npoints = INT((rtbse_env%pade_e_max - rtbse_env%pade_e_min)/rtbse_env%pade_e_step)
      ! Evaluate the evaluation grid
      IF (rtbse_env%pade_requested) THEN
         NULLIFY (rtbse_env%pade_x_eval)
         ALLOCATE (rtbse_env%pade_x_eval(rtbse_env%pade_npoints))
         DO i = 1, rtbse_env%pade_npoints
            rtbse_env%pade_x_eval(i) = CMPLX(rtbse_env%pade_e_step*REAL(i - 1, kind=dp), 0.0, kind=dp)
         END DO
      END IF

      ! Allocate moments matrices.
      ! In linearized RT-BSE these store the MO-active transformed dipole moments;
      ! in full RT-BSE they remain AO-sized (initialized from overlap template).
      NULLIFY (rtbse_env%moments)
      ALLOCATE (rtbse_env%moments(3))
      NULLIFY (rtbse_env%moments_field)
      ALLOCATE (rtbse_env%moments_field(3))
      DO k = 1, 3
         CALL cp_fm_create(rtbse_env%moments(k), mo_struct)
         CALL cp_fm_create(rtbse_env%moments_field(k), mo_struct)
      END DO

      ! Allocate space for density propagation and other operations.
      ! In linearized RT-BSE these workspaces are MO-active sized; in full RT-BSE
      ! they remain at the full AO size.
      NULLIFY (rtbse_env%rho_workspace)
      ALLOCATE (rtbse_env%rho_workspace(4))
      DO i = 1, SIZE(rtbse_env%rho_workspace)
         CALL cp_cfm_create(rtbse_env%rho_workspace(i), mo_struct)
         CALL cp_cfm_set_all(rtbse_env%rho_workspace(i), CMPLX(0.0, 0.0, kind=dp))
      END DO

      !MG To Do: Flags to avoid allocating unnecessary workspaces for methods which do not require them

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
      NULLIFY (rtbse_env%sigma_complex_workspace)
      ALLOCATE (rtbse_env%sigma_complex_workspace(1))
      CALL cp_cfm_create(rtbse_env%sigma_complex_workspace(1), bs_env%fm_ks_Gamma(1)%matrix_struct)
      CALL cp_cfm_set_all(rtbse_env%sigma_complex_workspace(1), CMPLX(0.0_dp, 0.0_dp, kind=dp))
      ! Allocate density matrix (MO-active sized when linearized; AO-sized otherwise)
      NULLIFY (rtbse_env%rho)
      ALLOCATE (rtbse_env%rho(rtbse_env%n_spin))
      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_create(rtbse_env%rho(i), matrix_struct=mo_struct)
      END DO
      ! Allocate additional space for AO density matrix
      ! in linearised RTBSE, where default is MO
      IF (my_linearized) THEN
         NULLIFY (rtbse_env%rho_ao_lrrtbse)
         ALLOCATE (rtbse_env%rho_ao_lrrtbse(rtbse_env%n_spin))
         DO i = 1, rtbse_env%n_spin
            CALL cp_cfm_create(rtbse_env%rho_ao_lrrtbse(i), matrix_struct=bs_env%fm_s_Gamma%matrix_struct)
         END DO
      END IF
      ! Create the inverse overlap matrix, for use in density propagation
      ! Start by creating the actual overlap matrix
      CALL cp_fm_create(rtbse_env%S_fm, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(rtbse_env%S_inv_fm, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_cfm_create(rtbse_env%S_cfm, bs_env%fm_s_Gamma%matrix_struct)

      ! Create the single particle hamiltonian
      ! Allocate workspace (MO-active sized in linearized RT-BSE; AO sized otherwise)
      NULLIFY (rtbse_env%ham_workspace)
      ALLOCATE (rtbse_env%ham_workspace(rtbse_env%n_spin))
      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_create(rtbse_env%ham_workspace(i), mo_struct)
         CALL cp_cfm_set_all(rtbse_env%ham_workspace(i), CMPLX(0.0, 0.0, kind=dp))
      END DO
      ! Now onto the Hamiltonian itself
      ! full RTBSE: Contains energy differences and Hartree/COHSEX ρ_0 parts 
      ! linearised RTBSE: Contains only the Hartree/COHSEX ρ_0 parts as Δε * Δρ(t) need to be updated 
      NULLIFY (rtbse_env%ham_reference)
      ALLOCATE (rtbse_env%ham_reference(rtbse_env%n_spin))
      DO i = 1, rtbse_env%n_spin
         CALL cp_cfm_create(rtbse_env%ham_reference(i), mo_struct)
      END DO
      ! Single particle Hamiltonian (Δε * Δρ(t)) for updates during timesteps in LR-RTBSE
      IF (my_linearized) THEN
         NULLIFY (rtbse_env%ham_reference_singleparticle)
         ALLOCATE (rtbse_env%ham_reference_singleparticle(rtbse_env%n_spin))
         DO i = 1, rtbse_env%n_spin
            CALL cp_cfm_create(rtbse_env%ham_reference_singleparticle(i), mo_struct)
         END DO
      END IF

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
         CALL cp_cfm_create(rtbse_env%ham_effective(i), mo_struct)
         CALL cp_cfm_set_all(rtbse_env%ham_effective(i), CMPLX(0.0, 0.0, kind=dp))
         CALL cp_cfm_create(rtbse_env%rho_new(i), mo_struct)
         CALL cp_cfm_set_all(rtbse_env%rho_new(i), CMPLX(0.0, 0.0, kind=dp))
         CALL cp_cfm_create(rtbse_env%rho_new_last(i), mo_struct)
         CALL cp_cfm_set_all(rtbse_env%rho_new_last(i), CMPLX(0.0, 0.0, kind=dp))
         CALL cp_cfm_create(rtbse_env%rho_M(i), mo_struct)
         CALL cp_cfm_set_all(rtbse_env%rho_M(i), CMPLX(0.0, 0.0, kind=dp))
         CALL cp_cfm_create(rtbse_env%rho_orig(i), mo_struct)
      END DO

      !For LR-RTBSE we need RK4 coefficients - create new workspace
      IF (my_linearized) THEN
         NULLIFY (rtbse_env%rk4_coefficients)
         ALLOCATE (rtbse_env%rk4_coefficients(4))
         DO i = 1, 4
            CALL cp_cfm_create(rtbse_env%rk4_coefficients(i), mo_struct)
            CALL cp_cfm_set_all(rtbse_env%rk4_coefficients(i), CMPLX(0.0, 0.0, kind=dp))
         END DO
      END IF

      ! Fields for exact diagonalisation
      NULLIFY (rtbse_env%real_eigvals)
      ALLOCATE (rtbse_env%real_eigvals(rtbse_env%n_ao))
      rtbse_env%real_eigvals(:) = 0.0_dp
      NULLIFY (rtbse_env%exp_eigvals)
      ALLOCATE (rtbse_env%exp_eigvals(rtbse_env%n_ao))
      rtbse_env%exp_eigvals(:) = CMPLX(0.0, 0.0, kind=dp)

      ! Workspace for FT - includes (in principle) the zeroth step and the extra last step
      NULLIFY (rtbse_env%moments_trace)
      ! TODO : Unite the number of steps with TD-DFT
      ALLOCATE (rtbse_env%moments_trace(rtbse_env%n_spin, 3, rtbse_env%sim_nsteps + 1), source=z_zero)
      NULLIFY (rtbse_env%field_trace)
      ALLOCATE (rtbse_env%field_trace(3, rtbse_env%sim_nsteps + 1), source=z_zero)
      NULLIFY (rtbse_env%time_trace)
      ALLOCATE (rtbse_env%time_trace(rtbse_env%sim_nsteps + 1), source=0.0_dp)

      ! Allocate self-energy parts and dynamic Hartree potential.
      ! In linearized RT-BSE these matrices hold the MO-active-sized result of the
      ! AO->MO transform; the AO-sized buffer is allocated as sigma_*_ao below.
      NULLIFY (rtbse_env%hartree_curr)
      NULLIFY (rtbse_env%sigma_SEX)
      NULLIFY (rtbse_env%sigma_COH)
      ALLOCATE (rtbse_env%hartree_curr(rtbse_env%n_spin))
      ALLOCATE (rtbse_env%sigma_SEX(rtbse_env%n_spin))
      ALLOCATE (rtbse_env%sigma_COH(rtbse_env%n_spin))
      DO i = 1, rtbse_env%n_spin
         CALL cp_fm_create(rtbse_env%sigma_COH(i), mo_struct)
         CALL cp_cfm_create(rtbse_env%sigma_SEX(i), mo_struct)
         CALL cp_fm_create(rtbse_env%hartree_curr(i), mo_struct)
         CALL cp_fm_set_all(rtbse_env%sigma_COH(i), 0.0_dp)
         CALL cp_cfm_set_all(rtbse_env%sigma_SEX(i), CMPLX(0.0, 0.0, kind=dp))
         CALL cp_fm_set_all(rtbse_env%hartree_curr(i), 0.0_dp)
      END DO
      ! AO-sized scratch buffers used by the linearized RT-BSE path
      IF (my_linearized) THEN
         NULLIFY (rtbse_env%hartree_curr_ao)
         NULLIFY (rtbse_env%sigma_SEX_ao)
         NULLIFY (rtbse_env%sigma_COH_ao)
         ALLOCATE (rtbse_env%hartree_curr_ao(rtbse_env%n_spin))
         ALLOCATE (rtbse_env%sigma_SEX_ao(rtbse_env%n_spin))
         ALLOCATE (rtbse_env%sigma_COH_ao(rtbse_env%n_spin))
         DO i = 1, rtbse_env%n_spin
            CALL cp_fm_create(rtbse_env%sigma_COH_ao(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
            CALL cp_cfm_create(rtbse_env%sigma_SEX_ao(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
            CALL cp_fm_create(rtbse_env%hartree_curr_ao(i), bs_env%fm_ks_Gamma(1)%matrix_struct)
            CALL cp_fm_set_all(rtbse_env%sigma_COH_ao(i), 0.0_dp)
            CALL cp_cfm_set_all(rtbse_env%sigma_SEX_ao(i), CMPLX(0.0, 0.0, kind=dp))
            CALL cp_fm_set_all(rtbse_env%hartree_curr_ao(i), 0.0_dp)
         END DO
         ! mo_active x mo_active real workspace pair for MO-side intermediates
         NULLIFY (rtbse_env%real_workspace_mo)
         ALLOCATE (rtbse_env%real_workspace_mo(2))
         DO i = 1, SIZE(rtbse_env%real_workspace_mo)
            CALL cp_fm_create(rtbse_env%real_workspace_mo(i), rtbse_env%fm_struct_mo_active)
            CALL cp_fm_set_all(rtbse_env%real_workspace_mo(i), 0.0_dp)
         END DO
         NULLIFY (rtbse_env%ao_mo_workspace)
         ALLOCATE (rtbse_env%ao_mo_workspace(1))
         CALL cp_fm_create(rtbse_env%ao_mo_workspace(1), rtbse_env%fm_struct_ao_mo_active)
         CALL cp_fm_set_all(rtbse_env%ao_mo_workspace(1), 0.0_dp)
         ! Truncated MO coefficient slabs C_active (n_ao x mo_active) for each spin.
         ! Filled in initialize_rtbse_env from bs_env%fm_mo_coeff_Gamma via submatrix copy.
         NULLIFY (rtbse_env%C_active)
         ALLOCATE (rtbse_env%C_active(rtbse_env%n_spin))
         DO i = 1, rtbse_env%n_spin
            CALL cp_fm_create(rtbse_env%C_active(i), rtbse_env%fm_struct_ao_mo_active)
            CALL cp_fm_set_all(rtbse_env%C_active(i), 0.0_dp)
         END DO
      END IF

      ! Allocate workspaces for get_sigma
      CALL create_sigma_workspace(rtbse_env)

      ! Depending on the chosen methods, allocate extra workspace
      CALL create_hartree_ri_workspace(rtbse_env)

   END SUBROUTINE create_rtbse_env

! **************************************************************************************************
!> \brief Simple reimplementation of cp_fm_release_pp1 for complex matrices
!> \param matrices cp_cfm_type(:)
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

      CALL cp_cfm_release_pa1(rtbse_env%ham_effective)
      CALL cp_cfm_release_pa1(rtbse_env%ham_workspace)
      CALL cp_fm_release(rtbse_env%sigma_COH)
      CALL cp_cfm_release_pa1(rtbse_env%sigma_SEX)
      CALL cp_fm_release(rtbse_env%hartree_curr)
      CALL cp_cfm_release_pa1(rtbse_env%ham_reference)
      IF (ASSOCIATED(rtbse_env%ham_reference_singleparticle)) THEN
         CALL cp_cfm_release_pa1(rtbse_env%ham_reference_singleparticle)
      END IF
      CALL cp_cfm_release_pa1(rtbse_env%rho)
      CALL cp_cfm_release_pa1(rtbse_env%rho_workspace)
      CALL cp_cfm_release_pa1(rtbse_env%rho_new)
      CALL cp_cfm_release_pa1(rtbse_env%rho_new_last)
      CALL cp_cfm_release_pa1(rtbse_env%rho_M)
      CALL cp_cfm_release_pa1(rtbse_env%rho_orig)
      IF (ASSOCIATED(rtbse_env%rk4_coefficients)) THEN
         CALL cp_cfm_release_pa1(rtbse_env%rk4_coefficients)
      END IF
      CALL cp_fm_release(rtbse_env%real_workspace)
      IF (ASSOCIATED(rtbse_env%sigma_complex_workspace)) CALL cp_cfm_release_pa1(rtbse_env%sigma_complex_workspace)
      CALL cp_fm_release(rtbse_env%S_inv_fm)
      CALL cp_fm_release(rtbse_env%S_fm)
      CALL cp_cfm_release(rtbse_env%S_cfm)

      CALL cp_fm_release(rtbse_env%moments)
      CALL cp_fm_release(rtbse_env%moments_field)

      CALL release_sigma_workspace(rtbse_env)

      CALL release_hartree_ri_workspace(rtbse_env)

      DEALLOCATE (rtbse_env%real_eigvals)
      DEALLOCATE (rtbse_env%exp_eigvals)
      DEALLOCATE (rtbse_env%moments_trace)
      DEALLOCATE (rtbse_env%field_trace)
      DEALLOCATE (rtbse_env%time_trace)

      IF (ASSOCIATED(rtbse_env%pol_elements)) DEALLOCATE (rtbse_env%pol_elements)
      IF (ASSOCIATED(rtbse_env%pade_x_eval)) DEALLOCATE (rtbse_env%pade_x_eval)

      ! Deallocate the neighbour list that is not deallocated in gw anymore
      IF (ASSOCIATED(rtbse_env%bs_env%nl_3c%ij_list)) CALL neighbor_list_3c_destroy(rtbse_env%bs_env%nl_3c)
      ! Release linearized-only AO scratches and MO-side workspaces
      IF (ASSOCIATED(rtbse_env%sigma_COH_ao)) CALL cp_fm_release(rtbse_env%sigma_COH_ao)
      IF (ASSOCIATED(rtbse_env%sigma_SEX_ao)) CALL cp_cfm_release_pa1(rtbse_env%sigma_SEX_ao)
      IF (ASSOCIATED(rtbse_env%hartree_curr_ao)) CALL cp_fm_release(rtbse_env%hartree_curr_ao)
      IF (ASSOCIATED(rtbse_env%real_workspace_mo)) CALL cp_fm_release(rtbse_env%real_workspace_mo)
      IF (ASSOCIATED(rtbse_env%ao_mo_workspace)) CALL cp_fm_release(rtbse_env%ao_mo_workspace)
      IF (ASSOCIATED(rtbse_env%C_active)) CALL cp_fm_release(rtbse_env%C_active)
      ! Release owned active-MO matrix structures
      IF (ASSOCIATED(rtbse_env%fm_struct_mo_active)) THEN
         CALL cp_fm_struct_release(rtbse_env%fm_struct_mo_active)
      END IF
      IF (ASSOCIATED(rtbse_env%fm_struct_ao_mo_active)) THEN
         CALL cp_fm_struct_release(rtbse_env%fm_struct_ao_mo_active)
      END IF
      ! Deallocate the storage for the environment itself
      DEALLOCATE (rtbse_env)
      ! Nullify to make sure it is not used again
      NULLIFY (rtbse_env)

   END SUBROUTINE release_rtbse_env

! **************************************************************************************************
!> \brief Determine the combined active MO window for linearized RT-BSE truncation.
!>
!> Evaluates BSE-like cutoff indices per spin from the requested single-particle spectrum
!> (G0W0 or KS Gamma-point eigenvalues) and collapses them into a single combined window
!> covering both spin channels by choosing the most inclusive bounds. Issues a CPWARN if the
!> spin-resolved cutoff candidates differ. When cutoffs are disabled (or the run is not
!> linearized RT-BSE), the window is set to the full MO range.
!> \param rtbse_env RT-BSE environment with cutoff values, n_ao, n_occ, n_spin, ham_reference_type
!>                  already populated.
!> \param my_linearized .TRUE. if linearized propagation is requested.
!> \param bs_env Bandstructure environment providing the eigenvalues.
! **************************************************************************************************
   SUBROUTINE determine_active_mo_window(rtbse_env, my_linearized, bs_env)
      TYPE(rtbse_env_type), POINTER                       :: rtbse_env
      LOGICAL, INTENT(IN)                                 :: my_linearized
      TYPE(post_scf_bandstructure_type), POINTER          :: bs_env
      CHARACTER(LEN=*), PARAMETER                         :: routineN = "determine_active_mo_window"

      INTEGER                                             :: handle, ispin, n_ao_full, n_virt
      INTEGER                                             :: homo_red, virt_red, homo_incl, virt_incl
      INTEGER                                             :: combined_first_occ, combined_last_virt
      INTEGER                                             :: first_occ_prev, last_virt_prev
      LOGICAL                                             :: spins_differ, do_truncation
      REAL(KIND=dp)                                       :: cutoff_occ, cutoff_empty

      CALL timeset(routineN, handle)

      n_ao_full = rtbse_env%n_ao
      cutoff_occ = rtbse_env%rtbse_energy_cutoff_occ
      cutoff_empty = rtbse_env%rtbse_energy_cutoff_empty
      do_truncation = my_linearized .AND. (cutoff_occ > 0.0_dp .OR. cutoff_empty > 0.0_dp)

      ! Default: full MO window
      rtbse_env%first_active_mo = 1
      rtbse_env%last_active_mo = n_ao_full
      rtbse_env%mo_active = n_ao_full
      rtbse_env%active_mo_truncation = .FALSE.

      IF (.NOT. do_truncation) THEN
         CALL timestop(handle)
         RETURN
      END IF

      combined_first_occ = n_ao_full
      combined_last_virt = 1
      first_occ_prev = -1
      last_virt_prev = -1
      spins_differ = .FALSE.

      DO ispin = 1, rtbse_env%n_spin
         n_virt = n_ao_full - rtbse_env%n_occ(ispin)
         IF (rtbse_env%ham_reference_type == rtp_bse_ham_g0w0) THEN
            CALL determine_cutoff_indices(bs_env%eigenval_G0W0(:, 1, ispin), &
                                          rtbse_env%n_occ(ispin), n_virt, &
                                          homo_red, virt_red, homo_incl, virt_incl, &
                                          cutoff_occ, cutoff_empty)
         ELSE
            CALL determine_cutoff_indices(bs_env%eigenval_scf_Gamma(:, ispin), &
                                          rtbse_env%n_occ(ispin), n_virt, &
                                          homo_red, virt_red, homo_incl, virt_incl, &
                                          cutoff_occ, cutoff_empty)
         END IF
         ! Translate the per-spin candidate to global MO indices [homo_incl, homo + virt_incl]
         IF (ispin > 1) THEN
            IF (homo_incl /= first_occ_prev .OR. (rtbse_env%n_occ(ispin) + virt_incl) /= last_virt_prev) THEN
               spins_differ = .TRUE.
            END IF
         END IF
         first_occ_prev = homo_incl
         last_virt_prev = rtbse_env%n_occ(ispin) + virt_incl
         combined_first_occ = MIN(combined_first_occ, homo_incl)
         combined_last_virt = MAX(combined_last_virt, rtbse_env%n_occ(ispin) + virt_incl)
      END DO

      IF (spins_differ) THEN
         CPWARN("RTBSE: spin-resolved active MO cutoff candidates differ; using combined window.")
      END IF

      rtbse_env%first_active_mo = combined_first_occ
      rtbse_env%last_active_mo = combined_last_virt
      rtbse_env%mo_active = combined_last_virt - combined_first_occ + 1
      rtbse_env%active_mo_truncation = (rtbse_env%mo_active < n_ao_full)

      CALL timestop(handle)
   END SUBROUTINE determine_active_mo_window

! **************************************************************************************************
!> \brief Allocates the workspaces for Hartree RI method
!> \note RI method calculates the Hartree contraction without the use of DBT, as it cannot emulate vectors
!> \param rtbse_env
!> \param qs_env Quickstep environment - entry point of calculation
!> \author Stepan Marek
!> \date 05.2024
! **************************************************************************************************
   SUBROUTINE create_hartree_ri_workspace(rtbse_env)
      TYPE(rtbse_env_type)                              :: rtbse_env
      TYPE(post_scf_bandstructure_type), POINTER        :: bs_env

      CALL get_qs_env(rtbse_env%qs_env, bs_env=bs_env)

      CALL dbcsr_create(rtbse_env%rho_dbcsr, name="Sparse density", template=bs_env%mat_ao_ao%matrix)
      CALL dbcsr_create(rtbse_env%v_ao_dbcsr, name="Sparse Hartree", template=bs_env%mat_ao_ao%matrix)

      CALL create_hartree_ri_3c(rtbse_env%rho_dbcsr, rtbse_env%int_3c_array, rtbse_env%n_ao, rtbse_env%n_RI, &
                                bs_env%basis_set_AO, bs_env%basis_set_RI, bs_env%i_RI_start_from_atom, &
                                bs_env%ri_metric, rtbse_env%qs_env, rtbse_env%unit_nr)
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
   SUBROUTINE create_sigma_workspace(rtbse_env)
      TYPE(rtbse_env_type)                               :: rtbse_env

      CALL create_sigma_workspace_qs_only(rtbse_env%qs_env, rtbse_env%screened_dbt, rtbse_env%w_dbcsr, &
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
   END SUBROUTINE create_sigma_workspace_qs_only
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
   END SUBROUTINE release_sigma_workspace
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
END MODULE rt_bse_types

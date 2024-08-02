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
   USE cp_cfm_types,                    ONLY: cp_cfm_p_type,&
                                              cp_cfm_set_all,&
                                              cp_cfm_create,&
                                              cp_fm_to_cfm,&
                                              cp_cfm_release
   USE cp_fm_basic_linalg,              ONLY: cp_fm_invert,&
                                              cp_fm_transpose,&
                                              cp_fm_column_scale,&
                                              cp_fm_scale_and_add
   USE cp_dbcsr_api,                    ONLY: dbcsr_type,&
                                              dbcsr_p_type,&
                                              dbcsr_create,&
                                              dbcsr_release,&
                                              dbcsr_print,&
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
                                              rtp_tdagw_ham_g0w0,&
                                              rtp_tdagw_ham_ks,&
                                              rtp_tdagw_hartree_dbt,&
                                              rtp_tdagw_hartree_pw,&
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_tdagw"

   #:include "tdagw_macros.fypp"
   
   PUBLIC :: tdagw_env_type,&
             create_tdagw_env,&
             release_tdagw_env

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
!> \param ham_single_particle Single particle Hamiltonian, taken from G0W0/DFT TODO : Only real?
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
!> \param hartree_dbt Tensor for hartree
!> \param screened_dbt Tensor for screened Coulomb interaction
!> \param sigma_dbt Tensor for self-energy
!> \param greens_dbt Tensor for greens function/density matrix
!> \param t_two_plus_one Tensor for hartree
!> \param t_one_plus_one Tensor for hartree
!> \param t_3c_w Tensor containing 3c integrals
!> \param t_3c_work_RI__AO_AO Tensor sigma contraction
!> \param t_3c_work_RI_AO__AO Tensor sigma contraction
!> \param t_3c_work2_RI_AO__AO Tensor sigma contraction
!> \param sigma_SEX Screened exchange self-energy
!> \param sigma_COH Coulomb hole self-energy
!> \param hartree_orig Original Hartree matrix
!> \param hartree_curr Current Hartree matrix
!> \param orders Parameter for exponential calculation (Taylor exp.)
!> \param etrs_max_iter Maximum number of ETRS iterations
!> \param hartree_method Which method to use for hartree_calculation
!> \param ham_single_particle_type Which Hamiltonian to use as single particle basis
!> \param mat_exp_method Which method to use for matrix exponentiation
!> \param unit_nr Number of output unit
!> \param int_3c_array Array containing the local 3c integrals
!> \author Stepan Marek (01.24)
! **************************************************************************************************
   TYPE tdagw_env_type
      INTEGER                                                   :: n_spin,&
                                                                   n_ao,&
                                                                   n_RI
      INTEGER, DIMENSION(2)                                     :: n_occ
      REAL(kind=dp)                                             :: spin_degeneracy
      REAL(kind=dp), DIMENSION(3)                               :: field
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: moments,&
                                                                   moments_field
      INTEGER                                                   :: sim_step,&
                                                                   sim_start,&
                                                                   sim_nsteps
      REAL(kind=dp)                                             :: sim_time,&
                                                                   sim_dt,&
                                                                   etrs_threshold,&
                                                                   exp_accuracy
      TYPE(dft_control_type), POINTER                           :: dft_control

      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: ham_effective
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: ham_single_particle
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: ham_workspace

      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: rho,&
                                                                   rho_new,&
                                                                   rho_new_last,&
                                                                   rho_M
      TYPE(cp_fm_p_type)                                        :: S_inv_fm,&
                                                                   S_fm
      TYPE(dbcsr_type), POINTER                                 :: S_inv,&
                                                                   rho_dbcsr
      ! TODO : Check how much is actually needed and whether it is not better to separate names
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER                 :: rho_workspace
      ! Workspace required for exact matrix exponentiation
      ! Complex workspace also used for hartree_pw
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER                :: complex_workspace,&
                                                                   complex_s
      REAL(kind=dp), DIMENSION(:), POINTER                      :: real_eigvals
      COMPLEX(kind=dp), DIMENSION(:), POINTER                   :: exp_eigvals
      ! Workspace required for hartree_pw
      TYPE(dbcsr_p_type)                                        :: v_dbcsr,&
                                                                   w_dbcsr
      ! TODO : Spins, k-points
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
      INTEGER                                                   :: etrs_max_iter
      INTEGER                                                   :: hartree_method
      INTEGER                                                   :: ham_single_particle_type
      INTEGER                                                   :: mat_exp_method
      INTEGER                                                   :: unit_nr 
      REAL(kind=dp), DIMENSION(:,:,:), POINTER                  :: int_3c_array

   END TYPE tdagw_env_type

CONTAINS

! **************************************************************************************************
!> \brief Allocates structures and prepares tdagw_env for run
!> \param tdagw_env tdagw_env_type that is initialised
!> \param qs_env Entry point of the calculation
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
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
      REAL(kind=dp), DIMENSION(:), POINTER                      :: custom_ref_point,&
                                                                   occupations
      REAL(kind=dp), DIMENSION(3)                               :: rpoint
      ! TYPE(cp_fm_type)                                          :: S_inv_fm
      TYPE(section_vals_type), POINTER                          :: input, bs_sec
      ! Extract the other types first

      CALL get_qs_env(qs_env,&
                      bs_env=bs_env,&
                      rtp=rtp,&
                      matrix_s=matrix_s,&
                      mos=mos,&
                      dft_control=tdagw_env%dft_control,&
                      input=input)
      bs_sec => section_vals_get_subs_vals(input, "PROPERTIES%BANDSTRUCTURE")
      ! TODO : Ability to run without the GW code?
      IF (.NOT. ASSOCIATED(bs_env)) THEN
         CPABORT("Cannot run RT-aGW without running GW calculation (PROPERTIES) before")
      END IF
      ! Number of spins
      tdagw_env%n_spin = bs_env%n_spin
      ! Number of atomic orbitals
      tdagw_env%n_ao = bs_env%n_ao
      ! Number of auxiliary basis orbitals
      tdagw_env%n_RI = bs_env%n_RI
      ! Number of occupied orbitals - for closed shell equals to half the number of electrons
      tdagw_env%n_occ(:) = bs_env%n_occ(:)
      ! Spin degeneracy - number of spins per orbital
      tdagw_env%spin_degeneracy = bs_env%spin_degeneracy
      ! Default field is zero
      tdagw_env%field(:) = 0.0_dp
      ! Default time is zero
      tdagw_env%sim_step = 0
      tdagw_env%sim_time = 0
      ! Time step is taken from rtp
      tdagw_env%sim_dt = rtp%dt
      ! Threshold for etrs is taken from the eps_energy from RT propagation
      tdagw_env%etrs_threshold = tdagw_env%dft_control%rtp_control%eps_ener
      ! TODO : Check whether dft_control can be dropped from tdagw_env
      tdagw_env%exp_accuracy = tdagw_env%dft_control%rtp_control%eps_exp
      ! Recover custom options
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%HARTREE_METHOD",&
                                i_val=tdagw_env%hartree_method)
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%TDAGW_HAMILTONIAN",&
                                i_val=tdagw_env%ham_single_particle_type)
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%MAX_ITER",&
                                i_val=tdagw_env%etrs_max_iter)
      CALL section_vals_val_get(input, "DFT%REAL_TIME_PROPAGATION%MAT_EXP",&
                                i_val=tdagw_env%mat_exp_method)
      ! Output unit number, recovered from the post_scf_bandstructure_type
      tdagw_env%unit_nr = bs_env%unit_nr
      IF (tdagw_env%unit_nr > 0) WRITE (tdagw_env%unit_nr, *) "exp_method", tdagw_env%mat_exp_method
      ! Allocate and determine the expansion orders, based on arguments used in TDDFT
      ! TODO : Last two arguments are indicating that we are using Taylor expansion and not doing ion motion
      ALLOCATE(tdagw_env%orders(2, tdagw_env%n_spin))
      IF (tdagw_env%ham_single_particle_type == rtp_tdagw_ham_g0w0) THEN
         @:SPIN_DO(i, re, im, tdagw_env%n_spin)
            CALL get_nsquare_norder(2.0_dp*MAX(ABS(MAXVAL(bs_env%eigenval_G0W0(:,bs_env%nkp_only_DOS,i))),&
                                    ABS(MINVAL(bs_env%eigenval_G0W0(:,bs_env%nkp_only_DOS,i)))),&
                                    tdagw_env%orders(1,i), tdagw_env%orders(2,i),&
                                    tdagw_env%dft_control%rtp_control%eps_exp,&
                                    1, .FALSE.)
         END DO
      ELSE
         ! KS Ham
         @:SPIN_DO(i,re,im,tdagw_env%n_spin)
            CALL get_nsquare_norder(2.0_dp*MAX(ABS(MAXVAL(bs_env%eigenval_scf_Gamma(:,i))),&
                                    ABS(MINVAL(bs_env%eigenval_scf_Gamma(:,i)))),&
                                    tdagw_env%orders(1,i), tdagw_env%orders(2,i),&
                                    tdagw_env%dft_control%rtp_control%eps_exp,&
                                    1, .FALSE.)
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
      ! Now, repeat without reference point to get the moments for field
      CALL build_local_moment_matrix(qs_env, moments_dbcsr_p, 1)
      NULLIFY(tdagw_env%moments_field)
      ALLOCATE(tdagw_env%moments_field(3))
      DO i=1,3
         NULLIFY(tdagw_env%moments_field(i)%matrix)
         ALLOCATE(tdagw_env%moments_field(i)%matrix)
         CALL cp_fm_create(tdagw_env%moments_field(i)%matrix, bs_env%fm_s_Gamma%matrix_struct)
         CALL copy_dbcsr_to_fm(moments_dbcsr_p(i)%matrix, tdagw_env%moments_field(i)%matrix)
      END DO

      ! Now can deallocate dbcsr matrices
      DO i=1,3
         CALL dbcsr_release(moments_dbcsr_p(i)%matrix)
      END DO
      DEALLOCATE(moments_dbcsr_p)

      ! Initiate density through results of DFT
      ! TODO : Handle cases, when initial method uses density only? Do such cases exist?
      ! DEBUG : change to initialisation directly from mo_coeffs
      ! CALL rt_initialize_rho_from_mos(rtp, mos)
      ! CALL get_rtp(rtp, rho_new=rho_dbcsr)
      ! IF (.NOT. tdagw_env%dft_control%rtp_control%apply_delta_pulse ) THEN
      !    CALL rt_initialize_rho_from_mos(rtp, mos)
      ! ELSE
      !    CALL rt_initialize_rho_from_mos(rtp, mos, rtp%mos%old)
      ! END IF
      ! Get the rho from RTP init - includes the RTP implementation of the delta kick
      ! CALL get_rtp(rtp, rho_new=rho_dbcsr)
      ! ! Redo the scaling to 1
      ! DO i=1,SIZE(rho_dbcsr)
      !    CALL dbcsr_scale(rho_dbcsr(i)%matrix, 0.5_dp)
      ! END DO
      ! DO i=1,2*tdagw_env%n_spin
      !    NULLIFY(tdagw_env%rho(i)%matrix)
      !    ALLOCATE(tdagw_env%rho(i)%matrix)
      !    CALL cp_fm_create(tdagw_env%rho(i)%matrix, matrix_struct=bs_env%fm_s_Gamma%matrix_struct)
      !    CALL copy_dbcsr_to_fm(rho_dbcsr(i)%matrix, tdagw_env%rho(i)%matrix)
      ! END DO
      ! Get the rho from fm_MOS - have to use the reimplementation of the delta kick
      ALLOCATE(tdagw_env%rho(2*tdagw_env%n_spin))
      ALLOCATE(occupations(tdagw_env%n_ao))
      ! Iterate over both spins
      @:SPIN_DO(i, re, im, tdagw_env%n_spin)
         occupations(:) = 0.0_dp
         occupations(1:tdagw_env%n_occ(i)) = 1.0_dp
         ! TODO : Handle complex orbitals - then we have non-trivial imaginary part
         ! Create real part
         NULLIFY(tdagw_env%rho(re)%matrix)
         ALLOCATE(tdagw_env%rho(re)%matrix)
         CALL cp_fm_create(tdagw_env%rho(re)%matrix, matrix_struct=bs_env%fm_s_Gamma%matrix_struct)
         ! Create imaginary part
         NULLIFY(tdagw_env%rho(im)%matrix)
         ALLOCATE(tdagw_env%rho(im)%matrix)
         CALL cp_fm_create(tdagw_env%rho(im)%matrix, matrix_struct=bs_env%fm_s_Gamma%matrix_struct)
         ! Populate real part
         CALL cp_fm_to_fm(bs_env%fm_mo_coeff_Gamma(i), tdagw_env%rho(im)%matrix)
         CALL cp_fm_column_scale(tdagw_env%rho(im)%matrix, occupations)
         CALL parallel_gemm("N", "T", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                            1.0_dp, bs_env%fm_mo_coeff_Gamma(i),&
                            tdagw_env%rho(im)%matrix,&
                            0.0_dp, tdagw_env%rho(re)%matrix)
         ! Populate imaginary part
         CALL cp_fm_set_all(tdagw_env%rho(im)%matrix, 0.0_dp)
      END DO
      DEALLOCATE(occupations)
      ! Create the inverse overlap matrix, for use in density propagation 
      CALL s_matrices_create(matrix_s, rtp)
      CALL get_rtp(rtp, S_inv=tdagw_env%S_inv)
      ! TODO : Check that dbcsr overlap inverse coincides with fm overlap inverse
      ! CALL dbcsr_print(tdagw_env%S_inv)
      NULLIFY(tdagw_env%S_inv_fm%matrix)
      ALLOCATE(tdagw_env%S_inv_fm%matrix)
      CALL cp_fm_create(tdagw_env%S_inv_fm%matrix, bs_env%fm_s_Gamma%matrix_struct)
      CALL copy_dbcsr_to_fm(tdagw_env%S_inv, tdagw_env%S_inv_fm%matrix)
      NULLIFY(tdagw_env%S_fm%matrix)
      ALLOCATE(tdagw_env%S_fm%matrix)
      CALL cp_fm_create(tdagw_env%S_fm%matrix, bs_env%fm_s_Gamma%matrix_struct)
      ! TODO : Handle different spins
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, tdagw_env%S_fm%matrix)

      ! Create the single particle hamiltonian
      ! Allocate workspace first
      NULLIFY(tdagw_env%ham_workspace)
      ALLOCATE(tdagw_env%ham_workspace(2*tdagw_env%n_spin))
      DO i=1,2*tdagw_env%n_spin
         ALLOCATE(tdagw_env%ham_workspace(i)%matrix)
         CALL cp_fm_create(tdagw_env%ham_workspace(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_set_all(tdagw_env%ham_workspace(i)%matrix, 0.0_dp)
      END DO
      ! Now onto the Hamiltonian itself
      ALLOCATE(tdagw_env%ham_single_particle(tdagw_env%n_spin))
      @:SPIN_DO(i, re, im, tdagw_env%n_spin)
         ALLOCATE(tdagw_env%ham_single_particle(i)%matrix)
         IF ( tdagw_env%ham_single_particle_type == rtp_tdagw_ham_g0w0 ) THEN
            ! G0W0 Hamiltonian
            CALL cp_fm_create(tdagw_env%ham_single_particle(i)%matrix, bs_env%fm_ks_Gamma(i)%matrix_struct)
            CALL cp_fm_to_fm(bs_env%fm_mo_coeff_Gamma(i), tdagw_env%ham_workspace(1)%matrix)
            ! TODO : Checked that Gamma point is the zero k-point
            ! C * Lambda
            CALL cp_fm_column_scale(tdagw_env%ham_workspace(1)%matrix, bs_env%eigenval_G0W0(:,1,i))
            ! C * Lambda * C^T
            CALL parallel_gemm("N", "T", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                               1.0_dp, tdagw_env%ham_workspace(1)%matrix, bs_env%fm_mo_coeff_Gamma(i),&
                               0.0_dp, tdagw_env%ham_workspace(2)%matrix)
            ! S * C * Lambda * C^T
            CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                               1.0_dp, tdagw_env%S_fm%matrix, tdagw_env%ham_workspace(2)%matrix,&
                               0.0_dp, tdagw_env%ham_workspace(1)%matrix)
            ! S * C * Lambda * C^T * S = H
            CALL parallel_gemm("N", "N", tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_ao,&
                               1.0_dp, tdagw_env%ham_workspace(1)%matrix, tdagw_env%S_fm%matrix,&
                               0.0_dp, tdagw_env%ham_single_particle(i)%matrix)
         ELSE
            ! KS Hamiltonian
            CALL cp_fm_create(tdagw_env%ham_single_particle(i)%matrix, bs_env%fm_ks_Gamma(i)%matrix_struct)
            CALL cp_fm_to_fm(bs_env%fm_ks_Gamma(i), tdagw_env%ham_single_particle(i)%matrix)
         END IF
      END DO

      ! Create the matrices and workspaces for effective Hamiltonian
      ALLOCATE(tdagw_env%ham_effective(2*tdagw_env%n_spin))
      ALLOCATE(tdagw_env%rho_new(2*tdagw_env%n_spin))
      ALLOCATE(tdagw_env%rho_new_last(2*tdagw_env%n_spin))
      ALLOCATE(tdagw_env%rho_M(2*tdagw_env%n_spin))
      DO i=1,2*tdagw_env%n_spin
         ALLOCATE(tdagw_env%ham_effective(i)%matrix)
         CALL cp_fm_create(tdagw_env%ham_effective(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_set_all(tdagw_env%ham_effective(i)%matrix, 0.0_dp)
         ALLOCATE(tdagw_env%rho_new(i)%matrix)
         CALL cp_fm_create(tdagw_env%rho_new(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_set_all(tdagw_env%rho_new(i)%matrix, 0.0_dp)
         ALLOCATE(tdagw_env%rho_new_last(i)%matrix)
         CALL cp_fm_create(tdagw_env%rho_new_last(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_set_all(tdagw_env%rho_new_last(i)%matrix, 0.0_dp)
         ALLOCATE(tdagw_env%rho_M(i)%matrix)
         CALL cp_fm_create(tdagw_env%rho_M(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_set_all(tdagw_env%rho_M(i)%matrix, 0.0_dp)
      END DO

      ! Allocate space for density propagation
      SELECT CASE (tdagw_env%mat_exp_method)
         CASE (do_taylor)
            ALLOCATE(tdagw_env%rho_workspace(3))
         CASE (do_exact)
            ALLOCATE(tdagw_env%rho_workspace(3))
         CASE (do_bch)
            ALLOCATE(tdagw_env%rho_workspace(4))
         CASE DEFAULT
            CPABORT("Only Taylor, exact and BCH matrix propagation implemented")
      END SELECT
      DO i=1,SIZE(tdagw_env%rho_workspace)
         ALLOCATE(tdagw_env%rho_workspace(i)%matrix)
         CALL cp_fm_create(tdagw_env%rho_workspace(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_set_all(tdagw_env%rho_workspace(i)%matrix, 0.0_dp)
      END DO

      ! Exact diagonalisation requires complex workspaces
      ! IF (tdagw_env%mat_exp_method == do_exact) THEN
      IF (.TRUE.) THEN
         NULLIFY(tdagw_env%complex_s)
         ALLOCATE(tdagw_env%complex_s(tdagw_env%n_spin))
         DO i=1,tdagw_env%n_spin
            NULLIFY(tdagw_env%complex_s(i)%matrix)
            ALLOCATE(tdagw_env%complex_s(i)%matrix)
            ! TODO : True spin dependence - from qs_env or present somewhere in bs_env?
            CALL cp_cfm_create(tdagw_env%complex_s(i)%matrix, bs_env%fm_s_Gamma%matrix_struct)
            CALL cp_cfm_set_all(tdagw_env%complex_s(i)%matrix, (0.0_dp, 0.0_dp))
            CALL cp_fm_to_cfm(msourcer=bs_env%fm_s_Gamma, mtarget=tdagw_env%complex_s(i)%matrix)
         END DO
         NULLIFY(tdagw_env%complex_workspace)
         ALLOCATE(tdagw_env%complex_workspace(4))
         DO i=1,SIZE(tdagw_env%complex_workspace)
            NULLIFY(tdagw_env%complex_workspace(i)%matrix)
            ALLOCATE(tdagw_env%complex_workspace(i)%matrix)
            CALL cp_cfm_create(tdagw_env%complex_workspace(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
            CALL cp_cfm_set_all(tdagw_env%complex_workspace(i)%matrix, (0.0_dp, 0.0_dp))
         END DO
         NULLIFY(tdagw_env%real_eigvals)
         ALLOCATE(tdagw_env%real_eigvals(tdagw_env%n_ao))
         tdagw_env%real_eigvals(:) = 0.0_dp
         NULLIFY(tdagw_env%exp_eigvals)
         ALLOCATE(tdagw_env%exp_eigvals(tdagw_env%n_ao))
         tdagw_env%exp_eigvals(:) = (0.0_dp, 0.0_dp)
      END IF

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
         CALL cp_fm_create(tdagw_env%sigma_COH(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_create(tdagw_env%sigma_SEX(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_create(tdagw_env%hartree_orig(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
         CALL cp_fm_create(tdagw_env%hartree_curr(i)%matrix, bs_env%fm_ks_Gamma(1)%matrix_struct)
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
         CALL create_hartree_ri_workspace(tdagw_env, qs_env)
      ELSEIF ( tdagw_env%hartree_method == rtp_tdagw_hartree_pw) THEN
         ! CALL create_hartree_pw_workspace(tdagw_env, qs_env)
         CPABORT("PW Hartree no longer supported")
      END IF
      ! Rho-dbcsr is used by both Hartree methods
      NULLIFY(tdagw_env%rho_dbcsr)
      ALLOCATE(tdagw_env%rho_dbcsr)
      CALL dbcsr_create(tdagw_env%rho_dbcsr, name="Sparse density", template=bs_env%mat_ao_ao%matrix)

      NULLIFY(tdagw_env%v_dbcsr%matrix)
      ALLOCATE(tdagw_env%v_dbcsr%matrix)
      
   END SUBROUTINE

! **************************************************************************************************
!> \brief Releases the environment allocated structures
!> \param tdagw_env 
!> \author Stepan Marek
!> \date 02.2024
! **************************************************************************************************
   SUBROUTINE release_tdagw_env(tdagw_env, qs_env)
      TYPE(tdagw_env_type)                                      :: tdagw_env
      TYPE(qs_environment_type), POINTER                        :: qs_env
      CHARACTER(len=*), PARAMETER                               :: routineN = "release_tdagw_env"
      INTEGER                                                   :: i

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
      CALL cp_fm_release(tdagw_env%S_inv_fm%matrix)

      CALL cp_fm_release(tdagw_env%moments)
      CALL cp_fm_release(tdagw_env%moments_field)

      ! TODO : Decide on where and when to releae this
      CALL dbcsr_release(tdagw_env%v_dbcsr%matrix)
      DEALLOCATE(tdagw_env%v_dbcsr%matrix)

      NULLIFY(tdagw_env%S_inv)
      DEALLOCATE(tdagw_env%orders)

      CALL release_sigma_workspace(tdagw_env)

      ! Depending on methods, deallocate additional workspace
      IF ( tdagw_env%hartree_method == rtp_tdagw_hartree_dbt ) THEN
         CALL release_hartree_dbt_workspace(tdagw_env)
      ELSEIF (tdagw_env%hartree_method == rtp_tdagw_hartree_pw) THEN
         ! CALL release_hartree_pw_workspace(qs_env, tdagw_env)
         CPABORT("PW Hartree no longer supported")
      END IF
      IF (ASSOCIATED(tdagw_env%rho_dbcsr)) THEN
         CALL dbcsr_release(tdagw_env%rho_dbcsr)
         DEALLOCATE(tdagw_env%rho_dbcsr)
      END IF
      ! IF (tdagw_env%mat_exp_method == do_exact) THEN
      IF (.TRUE.) THEN
         DO i=1,SIZE(tdagw_env%complex_s)
            CALL cp_cfm_release(tdagw_env%complex_s(i)%matrix)
         END DO
         DEALLOCATE(tdagw_env%complex_s)
         DO i=1,SIZE(tdagw_env%complex_workspace)
            CALL cp_cfm_release(tdagw_env%complex_workspace(i)%matrix)
         END DO
         DEALLOCATE(tdagw_env%complex_workspace)
         DEALLOCATE(tdagw_env%real_eigvals)
         DEALLOCATE(tdagw_env%exp_eigvals)
      END IF
      IF (tdagw_env%mat_exp_method == do_exact) THEN
         DO i=1,SIZE(tdagw_env%complex_s)
            CALL cp_cfm_release(tdagw_env%complex_s(i)%matrix)
         END DO
         DEALLOCATE(tdagw_env%complex_s)
         DO i=1,SIZE(tdagw_env%complex_workspace)
            CALL cp_cfm_release(tdagw_env%complex_workspace(i)%matrix)
         END DO
         DEALLOCATE(tdagw_env%complex_workspace)
         DEALLOCATE(tdagw_env%real_eigvals)
         DEALLOCATE(tdagw_env%exp_eigvals)
      END IF
      
   END SUBROUTINE
! **************************************************************************************************
!> \brief Allocates the workspaces for Hartree RI method
!> \note RI method calculates the Hartree contraction without the use of DBT, as it cannot emulate vectors
!> \param tdagw_env 
!> \param qs_env Quickstep environment - entry point of calculation
!> \author Stepan Marek
!> \date 05.2024
! **************************************************************************************************
   SUBROUTINE create_hartree_ri_workspace(tdagw_env, qs_env)
      TYPE(tdagw_env_type)                              :: tdagw_env
      TYPE(qs_environment_type), POINTER                :: qs_env
      CHARACTER(len=*), PARAMETER                       :: routineN = "create_hartree_ri_workspace"
      TYPE(post_scf_bandstructure_type), POINTER        :: bs_env
      REAL(kind=dp)                                     :: size_mb

      ! TODO : Implement option/decision to not precompute all the 3c integrals
      size_mb = REAL(tdagw_env%n_ao * tdagw_env%n_ao * tdagw_env%n_RI * STORAGE_SIZE(size_mb)) / (1024_dp * 1024_dp)
      IF (tdagw_env%unit_nr > 0) WRITE(tdagw_env%unit_nr, *) "Approximate size of int_3c in MB", size_mb

      CALL get_qs_env(qs_env, bs_env=bs_env)

      ALLOCATE(tdagw_env%int_3c_array(tdagw_env%n_ao, tdagw_env%n_ao, tdagw_env%n_RI))
      CALL build_3c_integral_block(tdagw_env%int_3c_array, qs_env, potential_parameter=bs_env%ri_metric,&
                              basis_j=bs_env%basis_set_AO, basis_k=bs_env%basis_set_AO,&
                              basis_i=bs_env%basis_set_RI,&
                              j_bf_start_from_atom=bs_env%i_ao_start_from_atom,&
                              k_bf_start_from_atom=bs_env%i_ao_start_from_atom,&
                              i_bf_start_from_atom=bs_env%i_RI_start_from_atom)
   END SUBROUTINE create_hartree_ri_workspace
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
      ALLOCATE(pdims(3), source=0)
      ALLOCATE(pcoord(3), source=0)
      ! Set elements to zero - otherwise can get erronous behaviour
      ! DO j=1,3
      !    pdims(j) = 0
      !    pcoord(j) = 0
      ! END DO
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

      ! 3c tensors - do not populate yet, just allocate distributions
      NULLIFY(tdagw_env%t_3c_work_RI__AO_AO)
      ALLOCATE(tdagw_env%t_3c_work_RI__AO_AO)
      CALL dbt_create(bs_env%t_RI__AO_AO, tdagw_env%t_3c_work_RI__AO_AO)

      NULLIFY(tdagw_env%hartree_dbt)
      ALLOCATE(tdagw_env%hartree_dbt)
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
      IF (ASSOCIATED(tdagw_env%t_3c_work_RI__AO_AO)) THEN
         CALL dbt_destroy(tdagw_env%t_3c_work_RI__AO_AO)
         DEALLOCATE(tdagw_env%t_3c_work_RI__AO_AO)
      END IF

      CALL dbt_destroy(tdagw_env%hartree_dbt)
      DEALLOCATE(tdagw_env%hartree_dbt)
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
      TYPE(cp_fm_type)                                   :: W_fm

      CALL get_qs_env(qs_env, bs_env=bs_env)

      ! t_3c_w
      NULLIFY(tdagw_env%t_3c_w)
      ALLOCATE(tdagw_env%t_3c_w)
      CALL dbt_create(bs_env%t_RI__AO_AO, tdagw_env%t_3c_w)
      ! TODO : Provide option/decision whether to store the 3c integrals precomputed
      CALL compute_3c_integrals(qs_env, bs_env, tdagw_env%t_3c_w)
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
      ! Populate screened_dbt from gw run
      ! TODO : Avoid without G0W0 run?
      NULLIFY(tdagw_env%w_dbcsr%matrix)
      ALLOCATE(tdagw_env%w_dbcsr%matrix)
      ! ! TODO : Debug - subtract the antihermitian part, create by numerical inacuraccy of FFT
      ! CALL cp_fm_create(W_fm, matrix_struct=bs_env%fm_W_MIC_freq_zero%matrix_struct)
      ! CALL cp_fm_transpose(bs_env%fm_W_MIC_freq_zero, W_fm)
      ! CALL cp_fm_scale_and_add(0.5_dp, W_fm, 0.5_dp, bs_env%fm_W_MIC_freq_zero)
      CALL dbcsr_create(tdagw_env%w_dbcsr%matrix, name="W", template=bs_env%mat_RI_RI%matrix)
      ! CALL copy_fm_to_dbcsr(W_fm, tdagw_env%w_dbcsr%matrix)
      IF (tdagw_env%ham_single_particle_type == rtp_tdagw_ham_g0w0) THEN
         CALL copy_fm_to_dbcsr(bs_env%fm_W_MIC_freq_zero, tdagw_env%w_dbcsr%matrix)
      ELSE
         ! Assume only Hartree available - this is added latero
         ! DEBUG : Try w_dbcsr instead - should have just a small effect
         ! CALL dbcsr_set(tdagw_env%w_dbcsr%matrix, 0.0_dp)
         CALL copy_fm_to_dbcsr(bs_env%fm_W_MIC_freq_zero, tdagw_env%w_dbcsr%matrix)
         CALL dbcsr_set(tdagw_env%w_dbcsr%matrix, 0.0_dp)
      END IF
      CALL dbt_create(tdagw_env%w_dbcsr%matrix, tdagw_env%screened_dbt)
      CALL dbt_copy_matrix_to_tensor(tdagw_env%w_dbcsr%matrix, tdagw_env%screened_dbt)
      ! PRINT *, "W_PQ matrix from BS_ENV"
      ! CALL cp_fm_write_formatted(bs_env%fm_W_MIC_freq_zero, 6)
      ! sigma_dbt
      NULLIFY(tdagw_env%sigma_dbt)
      ALLOCATE(tdagw_env%sigma_dbt)
      CALL dbt_create(bs_env%mat_ao_ao%matrix, tdagw_env%sigma_dbt)
      ! greens_dbt
      NULLIFY(tdagw_env%greens_dbt)
      ALLOCATE(tdagw_env%greens_dbt)
      CALL dbt_create(bs_env%mat_ao_ao%matrix, tdagw_env%greens_dbt)
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
      IF (ASSOCIATED(tdagw_env%w_dbcsr%matrix)) THEN
         CALL dbcsr_release(tdagw_env%w_dbcsr%matrix)
         DEALLOCATE(tdagw_env%w_dbcsr%matrix)
      END IF
   END SUBROUTINE
END MODULE

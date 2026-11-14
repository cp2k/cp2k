!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Green functions, screening and self-energy for periodic RI-RS GW.
!> \author Jan Wilhelm (with help of ChatGPT, GPT-6 Sol and GPT-6 Astra)
! **************************************************************************************************
MODULE gw_ri_rs_small_cell_full_kp_utils
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_add_on_diag,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_transpose,&
                                              cp_cfm_uplo_to_full
   USE cp_cfm_cholesky,                 ONLY: cp_cfm_cholesky_decompose,&
                                              cp_cfm_cholesky_invert
   USE cp_cfm_diag,                     ONLY: cp_cfm_power
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_set_submatrix,&
                                              cp_cfm_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE fft_tools,                       ONLY: BWFFT,&
                                              fft_alloc,&
                                              fft_dealloc
   USE gw_ri_rs_small_cell_full_kp_gamma,ONLY: add_gamma_corr_cell_correction,&
                                              add_gamma_corr_pair_cell,&
                                              gamma_point_correction,&
                                              singular_head_cell_factors
   USE gw_ri_rs_small_cell_full_kp_lattice,ONLY: &
        assemble_AA_contraction, block_range, build_full_cell_grid, build_ri_2c_lattice, &
        complex_bytes, gamma_auxiliary_charges, owned_mo_coefficients, real_bytes, &
        ri_2c_lattice_to_kp, ri_metric_head_zero
   USE gw_ri_rs_small_cell_full_kp_types,ONLY: gwkp_type,&
                                              ri_q_matrix_block_type
   USE gw_utils_fm,                     ONLY: local_complex_power
   USE input_constants,                 ONLY: evGW0,&
                                              rskp_backend_fft
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_coulomb_2c,               ONLY: add_ewald_reciprocal
   USE kpoint_lattice_fft,              ONLY: lattice_fft_execute,&
                                              lattice_fft_prepare_geometry,&
                                              lattice_fft_profile_type
   USE kpoint_transform,                ONLY: cell_to_k_transform,&
                                              k_to_real_cell_transform
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathconstants,                   ONLY: fourpi,&
                                              gaussi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_mem_avail_per_rank_GB,&
                                              mp_para_env_release,&
                                              mp_para_env_type,&
                                              mp_request_type,&
                                              mp_waitall
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE ri_rs_grid_setup,                ONLY: setup_ri_rs_grid
   USE ri_rs_kp_fitting,                ONLY: compute_Z_lP_q,&
                                              evaluate_mu_rl_R_values
   USE ri_rs_utils,                     ONLY: precompute_ri_rs_radii
   USE util,                            ONLY: sort

!$ USE omp_lib, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_small_cell_full_kp_utils'

   PUBLIC :: compute_Sigma_x, &
             compute_g_and_chi, &
             compute_W0, &
             compute_sigma_c_from_w_pq, &
             prepare_ri_rs, &
             release_gw_work
   PUBLIC :: product_to_k_pair_centered, gw_pair_fft_work_type

   ! Tolerance for a sum of fractional k-point coordinates to be a reciprocal lattice vector.
   REAL(KIND=dp), PARAMETER :: eps_kpoint = 1.0E-8_dp

   ! Blocking limits and named slots in the six existing CP2K matrix workspaces.
   INTEGER, PARAMETER :: column_block_cols = 32
   INTEGER, PARAMETER :: band_kpoints_per_pass = 64, pair_batch_size = 256
   INTEGER, PARAMETER :: w_chi = 1, w_metric = 2, w_tmp = 3, w_dielectric = 4, &
                         w_output = 5, n_w_matrices = 5

   ! Memory estimates: bytes per MiB and per GB (mp_mem_avail_per_rank_GB).
   REAL(KIND=dp), PARAMETER :: bytes_per_MiB = 1024.0_dp**2, bytes_per_GB = 1.0E9_dp
   ! A single stage may use stage_memory_share of the available memory. Stages that keep earlier
   ! data alive reserve reserve_share of it (at least min_reserve_GB) and add library_bytes.
   REAL(KIND=dp), PARAMETER :: stage_memory_share = 0.5_dp, reserve_share = 0.1_dp, &
                               min_reserve_GB = 1.0_dp, library_bytes = 256.0_dp*bytes_per_MiB
   ! Fixed workspace of the band projection and the Σ(k) of one projection round.
   REAL(KIND=dp), PARAMETER :: projection_workspace = 64.0_dp*bytes_per_MiB, &
                               projection_round_bytes = 512.0_dp*bytes_per_MiB

   ! Optional cache of stencil factors, shared by repeated transforms on one pair block.
   ! Geometry is still prepared by the common lattice FFT interface at each call.
   TYPE gw_pair_fft_batch_type
      LOGICAL :: use_combined = .FALSE.
      INTEGER :: axis_count = 0, combined_count = 0
      INTEGER(KIND=int_8) :: cached_count = 0_int_8
      INTEGER, ALLOCATABLE :: axis_class(:, :), combined_class(:)
      COMPLEX(KIND=dp), ALLOCATABLE :: axis_factor(:, :, :), combined_factor(:, :)
   END TYPE gw_pair_fft_batch_type

   TYPE gw_pair_fft_work_type
      PRIVATE
      LOGICAL :: initialized = .FALSE.
      INTEGER :: first_pair = 0, n_pairs = 0, width = 0, ngrid = 0, nmesh(3) = 0
      INTEGER, ALLOCATABLE :: cells(:, :)
      REAL(KIND=dp), ALLOCATABLE :: xkp(:, :)
      TYPE(gw_pair_fft_batch_type), ALLOCATABLE :: batches(:)
   END TYPE gw_pair_fft_work_type

CONTAINS

! **************************************************************************************************
!> \brief Prepare RI-RS radial data, grid values, and Bloch orbitals.
!> \param bs_env Band-structure settings and physical results.
!> \param gwkp_env Temporary periodic GW work data.
! **************************************************************************************************
   SUBROUTINE prepare_ri_rs(bs_env, gwkp_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT)                     :: gwkp_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'prepare_ri_rs'

      COMPLEX(dp), ALLOCATABLE                           :: phi_lmu_k(:, :, :)
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL precompute_ri_rs_radii(bs_env%ri_rs, bs_env%basis_set_AO, bs_env%basis_set_RI, &
                                  bs_env%eps_filter, bs_env%unit_nr)
      CALL setup_ri_rs_grid(bs_env, bs_env%ri_rs%grid_points)
      CALL compute_Z_lP_q(bs_env, bs_env%ri_rs%grid_points, &
                          Z_lP_q=gwkp_env%Z_lP_q, Z_lP_zero=gwkp_env%Z_lP_zero)
      CALL build_all_phi_k(bs_env, phi_lmu_k)
      CALL build_all_psi_k(bs_env, phi_lmu_k, gwkp_env%psi_lm_k)
      DEALLOCATE (phi_lmu_k)
      IF (bs_env%gw_flavour == evGW0) THEN
         CPASSERT(ALLOCATED(gwkp_env%eigenvalues_G))
         CPASSERT(ALLOCATED(bs_env%eigenval_evGW0))
         CPASSERT(SIZE(bs_env%eigenval_scf, 2) >= bs_env%nkp_scf_desymm)
         bs_env%eigenval_evGW0(:, :, :) = bs_env%eigenval_scf
      END IF

      CALL timestop(handle)

   END SUBROUTINE prepare_ri_rs

! **************************************************************************************************
!> \brief Form the grid Green functions and RI polarizability.
!> \param qs_env Quickstep data required by existing integral/lifecycle interfaces.
!> \param bs_env Band-structure settings and physical results.
!> \param gwkp_env Temporary periodic GW work data.
! **************************************************************************************************
   SUBROUTINE compute_G_and_chi(qs_env, bs_env, gwkp_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT)                     :: gwkp_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_G_and_chi'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(gwkp_env%group_env))
      gwkp_env%group_env => bs_env%para_env
      CALL gwkp_env%group_env%retain()
      gwkp_env%ngroups = 1
      gwkp_env%my_group = 0

      CALL gamma_point_correction(bs_env, qs_env, gwkp_env)
      CPASSERT(.NOT. ALLOCATED(gwkp_env%band_orbitals))
      CALL compute_chi(bs_env, gwkp_env)

      CALL timestop(handle)

   END SUBROUTINE compute_G_and_chi

! **************************************************************************************************
!> \brief Compute ε, Wᶜ and its imaginary-time transform, Eqs. (7)-(8).
!> \param qs_env Required by the existing two-centre integral interfaces.
!> \param bs_env Band-structure settings and physical results.
!> \param gwkp_env Periodic GW work data.
! **************************************************************************************************
   SUBROUTINE compute_W0(qs_env, bs_env, gwkp_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT)                     :: gwkp_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_W0'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: chi_k_w(:, :), MWM_k_w(:, :)
      INTEGER                                            :: handle, i_batch, i_group_size, ikp, &
                                                            ikp_local, j_w, my_qgroup, n_qgroups, &
                                                            n_RI, nkp, nkp_local, q_group, &
                                                            q_owner, q_root
      INTEGER, ALLOCATABLE                               :: group_load(:), q_batch(:, :), &
                                                            q_partner(:), q_roots(:)
      LOGICAL, ALLOCATABLE                               :: solves_q(:)
      REAL(KIND=dp)                                      :: available_GB, required_bytes, step_start
      REAL(KIND=dp), ALLOCATABLE                         :: frequency_elapsed(:)
      TYPE(cp_blacs_env_type), POINTER                   :: w_blacs
      TYPE(cp_cfm_type)                                  :: w_matrices(n_w_matrices)
      TYPE(cp_fm_struct_type), POINTER                   :: w_matrix_struct
      TYPE(mp_para_env_type), POINTER                    :: w_group_env

      CALL timeset(routineN, handle)

      CALL compute_metric_coulomb_factor(bs_env, qs_env, gwkp_env)
      ! Σˣ and Σᶜ use the band-path orbitals, which include the SCF mesh.
      DEALLOCATE (gwkp_env%psi_lm_k)
      n_RI = bs_env%n_RI
      nkp = bs_env%nkp_scf_desymm
      IF (bs_env%nkp_chi_eps_W_orig_plus_extra /= nkp) THEN
         CPABORT("Periodic RI-RS GW requires the complete SCF transfer mesh")
      END IF
      IF (ANY(bs_env%nkp_grid_chi_eps_W_orig /= &
              bs_env%kpoints_scf_desymm%nkp_grid)) THEN
         CPABORT("Periodic RI-RS GW requires the SCF mesh dimensions for W")
      END IF
      nkp_local = SIZE(gwkp_env%metric_coulomb_factor, 3)
      ! Distributed solver matrices plus equal workspace reserve; owner buffers remain local.
      ! Select the smallest feasible divisor, independently of the number of q points.
      CALL mp_mem_avail_per_rank_GB(bs_env%para_env, available_GB)
      DO i_group_size = 1, bs_env%para_env%num_pe
         IF (MODULO(bs_env%para_env%num_pe, i_group_size) /= 0) CYCLE
         required_bytes = REAL(n_RI, dp)**2*(2.0_dp*complex_bytes + &
                                             2.0_dp*complex_bytes*n_w_matrices/REAL(i_group_size, dp)) + &
                          complex_bytes*REAL(n_RI, dp)*REAL(column_block_cols, dp)
         IF (required_bytes <= stage_memory_share*MAX(0.0_dp, available_GB)*bytes_per_GB) EXIT
      END DO
      IF (i_group_size > bs_env%para_env%num_pe) THEN
         CPABORT('Insufficient node memory for the distributed GW W(q) solver')
      END IF
      n_qgroups = bs_env%para_env%num_pe/i_group_size
      CALL grouped_w_partition(bs_env%para_env, n_qgroups, my_qgroup, q_roots)
      ! Time reversal gives W(-q) = W(q)*: compute_chi provides χ for the first q of each ±q pair,
      ! and the owner of -q receives the conjugate block after the last batch.
      CALL minus_q_partners(bs_env%kpoints_chi_eps_W%xkp(1:3, 1:nkp), q_partner)
      ALLOCATE (solves_q(nkp))
      solves_q(:) = [(q_partner(ikp) >= ikp, ikp=1, nkp)]
      ALLOCATE (group_load(0:n_qgroups - 1), SOURCE=0)
      DO ikp = 1, nkp
         IF (solves_q(ikp)) group_load(MODULO(ikp, n_qgroups)) = group_load(MODULO(ikp, n_qgroups)) + 1
      END DO
      ALLOCATE (q_batch(MAXVAL(group_load), 0:n_qgroups - 1), SOURCE=0)
      group_load(:) = 0
      DO ikp = 1, nkp
         IF (.NOT. solves_q(ikp)) CYCLE
         q_group = MODULO(ikp, n_qgroups)
         group_load(q_group) = group_load(q_group) + 1
         q_batch(group_load(q_group), q_group) = ikp
      END DO
      ALLOCATE (gwkp_env%W_PQ_tau(nkp_local))
      IF (nkp_local > 0 .OR. bs_env%para_env%mepos == q_roots(my_qgroup + 1)) THEN
         ALLOCATE (chi_k_w(n_RI, n_RI), MWM_k_w(n_RI, n_RI))
      ELSE
         ALLOCATE (chi_k_w(0, 0), MWM_k_w(0, 0))
      END IF

      ALLOCATE (w_group_env)
      CALL w_group_env%from_split(bs_env%para_env, my_qgroup)
      CALL grouped_w_create(w_matrices, w_blacs, w_matrix_struct, w_group_env, n_RI)
      ALLOCATE (frequency_elapsed(bs_env%num_time_freq_points), SOURCE=0.0_dp)
      ! Each group holds one q at a time, reusing K(q) for all frequencies. Group MODULO(q, n_qgroups)
      ! handles q, so with one rank per group the owner of q solves it without communication.
      ! Ordered world transfers also cover incomplete batches and idle groups.
      DO i_batch = 1, SIZE(q_batch, 1)
         DO q_group = 0, n_qgroups - 1
            ikp = q_batch(i_batch, q_group)
            IF (ikp == 0) CYCLE
            q_root = q_roots(q_group + 1)
            q_owner = MODULO(ikp, bs_env%para_env%num_pe)
            ikp_local = 1 + (ikp - 1)/bs_env%para_env%num_pe
            IF (bs_env%para_env%mepos == q_owner) THEN
               CALL grouped_w_copy_matrix(bs_env%para_env, n_RI, q_owner, q_root, &
                                          gwkp_env%metric_coulomb_factor(:, :, ikp_local), chi_k_w)
               ALLOCATE (gwkp_env%W_PQ_tau(ikp_local)%values( &
                         n_RI, n_RI, bs_env%num_time_freq_points), SOURCE=z_zero)
            ELSE
               CALL grouped_w_copy_matrix(bs_env%para_env, n_RI, q_owner, q_root, MWM_k_w, chi_k_w)
            END IF
            IF (q_group == my_qgroup) THEN
               CALL set_owner_matrix(w_matrices(w_metric), chi_k_w, w_group_env, 0, n_RI)
            END IF
         END DO
         DO j_w = 1, bs_env%num_time_freq_points
            step_start = m_walltime()
            DO q_group = 0, n_qgroups - 1
               ikp = q_batch(i_batch, q_group)
               IF (ikp == 0) CYCLE
               q_owner = MODULO(ikp, bs_env%para_env%num_pe)
               ikp_local = 1 + (ikp - 1)/bs_env%para_env%num_pe
               IF (bs_env%para_env%mepos == q_owner) THEN
                  MWM_k_w(:, :) = gwkp_env%chi_q_frequency(ikp_local)%values(:, :, j_w)
                  ! The RI fit gives occupied-virtual pair densities a spurious charge m_ia that makes
                  ! ε_00 - 1 ~ C(ω)/q and the q sum log-divergent. As in the Γ model (P = 1 - c c^H),
                  ! remove it from the fitted coefficients, d' = d - (c/|c|) m, i.e. Z^H ρ' = Z^H ρ - e m:
                  ! χ' = χ - e X^H - X e^H + L0 e e^H. m and d are taken at q = 0, so the physical
                  ! O(q) dipole charge is kept; the error is O(q).
                  CALL zgerc(n_RI, n_RI, CMPLX(gwkp_env%leak_l0(j_w), 0.0_dp, dp), &
                             gwkp_env%leak_e(:, ikp_local), 1, gwkp_env%leak_e(:, ikp_local), 1, &
                             MWM_k_w, n_RI)
                  CALL zgerc(n_RI, n_RI, -z_one, gwkp_env%leak_e(:, ikp_local), 1, &
                             gwkp_env%leak_x(:, j_w, ikp_local), 1, MWM_k_w, n_RI)
                  CALL zgerc(n_RI, n_RI, -z_one, gwkp_env%leak_x(:, j_w, ikp_local), 1, &
                             gwkp_env%leak_e(:, ikp_local), 1, MWM_k_w, n_RI)
               END IF
               CALL grouped_w_copy_matrix(bs_env%para_env, n_RI, q_owner, q_roots(q_group + 1), MWM_k_w, chi_k_w)
            END DO
            IF (q_batch(i_batch, my_qgroup) > 0) THEN
               CALL grouped_w_solve(w_matrices, w_group_env, n_RI, 0, chi_k_w, &
                                    bs_env%eps_eigval_mat_RI, MWM_k_w)
            END IF
            DO q_group = 0, n_qgroups - 1
               ikp = q_batch(i_batch, q_group)
               IF (ikp == 0) CYCLE
               q_owner = MODULO(ikp, bs_env%para_env%num_pe)
               ikp_local = 1 + (ikp - 1)/bs_env%para_env%num_pe
               CALL grouped_w_copy_matrix(bs_env%para_env, n_RI, q_roots(q_group + 1), q_owner, MWM_k_w, chi_k_w)
               IF (bs_env%para_env%mepos == q_owner) THEN
                  gwkp_env%W_PQ_tau(ikp_local)%values(:, :, j_w) = chi_k_w
                  IF (j_w == bs_env%num_time_freq_points) DEALLOCATE (gwkp_env%chi_q_frequency(ikp_local)%values)
               END IF
            END DO
            frequency_elapsed(j_w) = frequency_elapsed(j_w) + m_walltime() - step_start
         END DO
         DO q_group = 0, n_qgroups - 1
            ikp = q_batch(i_batch, q_group)
            IF (ikp == 0) CYCLE
            IF (bs_env%para_env%mepos /= MODULO(ikp, bs_env%para_env%num_pe)) CYCLE
            CALL cosine_transform_ri_q(gwkp_env%W_PQ_tau(1 + (ikp - 1)/bs_env%para_env%num_pe), &
                                       bs_env%time_frequency_grid%cosine_frequency_to_time)
         END DO
      END DO
      CALL receive_minus_q_blocks(bs_env%para_env, q_partner, solves_q, n_RI, &
                                  bs_env%num_time_freq_points, gwkp_env%W_PQ_tau)
      DO j_w = 1, bs_env%num_time_freq_points
         CALL bs_env%para_env%max(frequency_elapsed(j_w))
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I3,A,I0,T56,A,T70,F10.1,A)') &
               'GW W(iω) point ', j_w, ' / ', bs_env%num_time_freq_points, &
               'Exec. time:', frequency_elapsed(j_w), ' s'
            CALL m_flush(bs_env%unit_nr)
         END IF
      END DO
      DEALLOCATE (frequency_elapsed)

      CALL grouped_w_release(w_matrices, w_matrix_struct, w_blacs)
      CALL mp_para_env_release(w_group_env)

      ! Release frequency-domain matrices after constructing W(q,tau).
      DEALLOCATE (chi_k_w, MWM_k_w)

      CPASSERT(gwkp_env%ngroups == 1 .AND. gwkp_env%my_group == 0)
      CPASSERT(gwkp_env%group_env%num_pe == bs_env%para_env%num_pe)
      CPASSERT(gwkp_env%group_env%mepos == bs_env%para_env%mepos)

      DEALLOCATE (group_load, q_batch, q_partner, solves_q, gwkp_env%leak_e, gwkp_env%leak_x, gwkp_env%leak_l0)
      IF (ALLOCATED(gwkp_env%chi_q_frequency)) DEALLOCATE (gwkp_env%chi_q_frequency)

      CALL timestop(handle)

   END SUBROUTINE compute_W0

! **************************************************************************************************
!> \brief Evaluate exchange and prepare band orbitals for correlation.
!> \param bs_env Band-structure settings and physical results.
!> \param gwkp_env Temporary periodic GW work data.
! **************************************************************************************************
   SUBROUTINE compute_Sigma_x(bs_env, gwkp_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT)                     :: gwkp_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_x'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL collect_band_bloch_orbitals(bs_env, gwkp_env%band_orbitals)
      CALL compute_Sigma_x_from_V(bs_env, gwkp_env)
      CALL prepare_self_energy_groups(bs_env, gwkp_env)
      IF (gwkp_env%ngroups > 1) CALL band_orbitals_to_groups(bs_env, gwkp_env)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x

! **************************************************************************************************
!> \brief Retain Z and fixed W in the smallest memory-feasible time groups.
!>
!> Each group owns τ = group + 1, group + 1 + ngroups, ... . Within a group,
!> Z and W use owner MOD(q,g). W remains in the auxiliary basis. A single q block is broadcast during this
!> one-time redistribution; each auxiliary W(q,τ) goes only to its owning rank.
!> \param bs_env Band-structure settings.
!> \param gwkp_env Work arrays and their MPI distribution.
! **************************************************************************************************
   SUBROUTINE prepare_self_energy_groups(bs_env, gwkp_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT)                     :: gwkp_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_self_energy_groups'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: Z_group(:, :, :)
      INTEGER :: handle, i_dest, i_group_size, i_root, i_source, i_tag, iq, iq_local, it, &
         it_local, n_ao, n_aux, n_grid, n_k, n_local_k, n_local_q, n_local_time, n_spin, &
         n_stored_q, n_threads, n_time, n_world
      REAL(KIND=dp) :: available_GB, block_bytes, budget_bytes, grid_bytes, input_bytes, n_pairs, &
         peak_bytes, persistent_bytes, result_bytes
      TYPE(mp_comm_type)                                 :: column_env
      TYPE(ri_q_matrix_block_type), ALLOCATABLE          :: W_group(:)

      CALL timeset(routineN, handle)

      n_threads = 1
!$    n_threads = omp_get_max_threads()
      n_world = bs_env%para_env%num_pe
      n_grid = bs_env%ri_rs%n_grid_points
      n_ao = bs_env%n_ao
      n_spin = bs_env%n_spin
      n_aux = bs_env%n_RI
      n_k = bs_env%nkp_scf_desymm
      n_time = bs_env%num_time_freq_points
      CALL mp_mem_avail_per_rank_GB(bs_env%para_env, available_GB)
      CALL bs_env%para_env%min(available_GB)
      input_bytes = complex_bytes*(REAL(SIZE(gwkp_env%Z_lP_q, KIND=int_8), dp) + &
                                   REAL(ri_q_storage_elements(gwkp_env%W_PQ_tau), dp))
      CALL bs_env%para_env%min(input_bytes)
      ! World inputs remain alive during redistribution. Count the complete new
      ! group storage in addition to measured usage; retain a library reserve.
      budget_bytes = budget_after_reserve(available_GB)
      block_bytes = complex_bytes*REAL(n_grid, dp)*REAL(n_aux, dp)
      result_bytes = 2.0_dp*complex_bytes*REAL(n_ao, dp)*REAL(SIZE(bs_env%eigenval_scf, 2), dp)* &
                     REAL(n_time*n_spin, dp)
      DO i_group_size = MIN(2, MAX(1, n_world/2)), n_world
         IF (MODULO(n_world, i_group_size) /= 0) CYCLE
         gwkp_env%ngroups = n_world/i_group_size
         IF (gwkp_env%ngroups > n_time) CYCLE
         n_local_k = (n_k + i_group_size - 1)/i_group_size
         n_local_time = (n_time + gwkp_env%ngroups - 1)/gwkp_env%ngroups
         n_stored_q = (n_k + i_group_size - 1)/i_group_size
         persistent_bytes = complex_bytes*(REAL(n_stored_q, dp)*REAL(n_grid, dp)*REAL(n_aux, dp) + &
                                           REAL(n_aux, dp)**2*REAL(n_stored_q*n_local_time, dp))
         n_pairs = REAL(CEILING(REAL(n_grid, dp)**2/REAL(i_group_size, dp), KIND=int_8), dp)
         ! G(k), gwkp_env%W_PQ_tau(R), and Σ(R) coexist: 16 + 8 + 8 bytes per pair and mesh point,
         ! next to the band orbitals of the owned band k points.
         ! q matrices, band projection, and Green blocks belong to separate stages.
         grid_bytes = n_pairs*(complex_bytes + 2.0_dp*real_bytes)*REAL(n_k, dp) + &
                      complex_bytes*REAL(n_grid, dp)*REAL(n_ao, dp)*REAL(n_spin, dp)* &
                      REAL((bs_env%nkp_bs_and_DOS + i_group_size - 1)/i_group_size, dp) + &
                      MAX(complex_bytes*REAL(n_grid, dp)**2*REAL(n_local_k, dp) + &
                          complex_bytes*REAL(n_aux, dp)**2*REAL(n_local_k, dp) + &
                          complex_bytes*REAL(n_threads, dp)*REAL(n_grid, dp)*REAL(n_ao + n_aux, dp), &
                          projection_workspace_bytes(n_grid, i_group_size, &
                                                     PRODUCT(bs_env%kpoints_scf_desymm%nkp_grid + 2), n_ao, n_threads)) + &
                      2.0_dp*complex_bytes*REAL(n_threads*n_k, dp) + library_bytes
         peak_bytes = MAX(persistent_bytes + block_bytes, &
                          persistent_bytes - input_bytes + grid_bytes + result_bytes)
         IF (gwkp_env%ngroups == 1) peak_bytes = grid_bytes + result_bytes
         IF (peak_bytes <= budget_bytes) EXIT
      END DO
      IF (i_group_size > n_world) THEN
         CPABORT('Insufficient memory per MPI rank for GW self-energy time groups')
      END IF
      gwkp_env%ngroups = n_world/i_group_size
      gwkp_env%my_group = bs_env%para_env%mepos/i_group_size
      CALL mp_para_env_release(gwkp_env%group_env)
      ALLOCATE (gwkp_env%group_env)
      CALL gwkp_env%group_env%from_split(bs_env%para_env, gwkp_env%my_group)
      CPASSERT(gwkp_env%group_env%num_pe == i_group_size)
      IF (gwkp_env%ngroups == 1) THEN
         CALL timestop(handle)
         RETURN
      END IF
      n_local_q = COUNT([(MODULO(iq, i_group_size) == gwkp_env%group_env%mepos, iq=1, n_k)])
      n_local_time = (n_time - gwkp_env%my_group + gwkp_env%ngroups - 1)/gwkp_env%ngroups
      ALLOCATE (W_group(n_local_q))
      ! Ranks with the same within-group index need the same k/q slices.
      ! Communicate only along this column, independently of the other columns.
      CALL column_env%from_split(bs_env%para_env, gwkp_env%group_env%mepos)
      CALL gather_group_columns(bs_env, gwkp_env, column_env, i_group_size, gwkp_env%group_env%mepos, Z_group)
      iq_local = 0
      DO iq = MODULO(gwkp_env%group_env%mepos - 1, i_group_size) + 1, n_k, i_group_size
         i_root = MODULO(iq, n_world)/i_group_size
         iq_local = iq_local + 1
         ALLOCATE (W_group(iq_local)%values(n_aux, n_aux, n_local_time))
         DO it = 1, n_time
            i_dest = MODULO(it - 1, gwkp_env%ngroups)
            IF (column_env%mepos == i_root) THEN
               IF (i_dest == i_root) THEN
                  it_local = (it - 1)/gwkp_env%ngroups + 1
                  W_group(iq_local)%values(:, :, it_local) = gwkp_env%W_PQ_tau((iq - 1)/n_world + 1)%values(:, :, it)
               ELSE
                  CALL column_env%send(gwkp_env%W_PQ_tau((iq - 1)/n_world + 1)%values(:, :, it), i_dest, 0)
               END IF
            ELSE IF (column_env%mepos == i_dest) THEN
               i_source = i_root
               i_tag = 0
               it_local = (it - 1)/gwkp_env%ngroups + 1
               CALL column_env%recv(W_group(iq_local)%values(:, :, it_local), i_source, i_tag)
            END IF
         END DO
         IF (column_env%mepos == i_root) DEALLOCATE (gwkp_env%W_PQ_tau((iq - 1)/n_world + 1)%values)
      END DO
      CALL column_env%free()
      CALL MOVE_ALLOC(Z_group, gwkp_env%Z_lP_q)
      CALL MOVE_ALLOC(W_group, gwkp_env%W_PQ_tau)

      CALL timestop(handle)

   END SUBROUTINE prepare_self_energy_groups

! **************************************************************************************************
!> \brief Copy the Z(q), and optionally the orbitals ψ(k), needed by one rank of a time group
!> from their world owners.
!> World rank r owns k with MOD(k-1,n_world)=r and q with MOD(q,n_world)=r. Group rank i
!> receives k with MOD(k-1,g)=i and q with MOD(q,g)=i. column_env connects the ranks with
!> the same group rank; its rank is the time-group index.
!> \param bs_env Band-structure settings.
!> \param gwkp_env World-distributed ψ(k) and Z(q).
!> \param column_env Ranks with the same within-group index.
!> \param group_size Ranks per time group g.
!> \param group_rank Within-group index i.
!> \param Z_group Z(q) of this group rank, indexed (grid, RI, local q).
!> \param psi_group ψ(k) of this group rank, indexed (grid, AO, local k, spin).
! **************************************************************************************************
   SUBROUTINE gather_group_columns(bs_env, gwkp_env, column_env, group_size, group_rank, Z_group, psi_group)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(IN)                        :: gwkp_env
      TYPE(mp_comm_type), INTENT(IN)                     :: column_env
      INTEGER, INTENT(IN)                                :: group_size, group_rank
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: Z_group(:, :, :)
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT), &
         OPTIONAL                                        :: psi_group(:, :, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gather_group_columns'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: psi_block(:, :, :), Z_block(:, :)
      INTEGER                                            :: handle, i_root, ik, ik_local, iq, &
                                                            iq_local, n_grid, n_k, n_world

      CALL timeset(routineN, handle)

      n_world = bs_env%para_env%num_pe
      n_grid = bs_env%ri_rs%n_grid_points
      n_k = bs_env%nkp_scf_desymm
      IF (PRESENT(psi_group)) THEN
         ALLOCATE (psi_group(n_grid, bs_env%n_ao, COUNT([(MODULO(ik - 1, group_size) == group_rank, ik=1, n_k)]), &
                             bs_env%n_spin))
         ALLOCATE (psi_block(n_grid, bs_env%n_ao, bs_env%n_spin))
         ik_local = 0
         DO ik = group_rank + 1, n_k, group_size
            i_root = MODULO(ik - 1, n_world)/group_size
            IF (column_env%mepos == i_root) psi_block(:, :, :) = gwkp_env%psi_lm_k(:, :, (ik - 1)/n_world + 1, :)
            CALL column_env%bcast(psi_block, i_root)
            ik_local = ik_local + 1
            psi_group(:, :, ik_local, :) = psi_block
         END DO
         DEALLOCATE (psi_block)
      END IF
      ALLOCATE (Z_group(n_grid, bs_env%n_RI, COUNT([(MODULO(iq, group_size) == group_rank, iq=1, n_k)])))
      ALLOCATE (Z_block(n_grid, bs_env%n_RI))
      iq_local = 0
      DO iq = MODULO(group_rank - 1, group_size) + 1, n_k, group_size
         i_root = MODULO(iq, n_world)/group_size
         IF (column_env%mepos == i_root) Z_block(:, :) = gwkp_env%Z_lP_q(:, :, (iq - 1)/n_world + 1)
         CALL column_env%bcast(Z_block, i_root)
         iq_local = iq_local + 1
         Z_group(:, :, iq_local) = Z_block
      END DO
      DEALLOCATE (Z_block)

      CALL timestop(handle)

   END SUBROUTINE gather_group_columns

! **************************************************************************************************
!> \brief Compute chi_PQ^R(i tau) through the periodic real-space grid representation.
!> \param bs_env Band-structure settings and physical results.
!> \param gwkp_env Periodic GW work data.
! **************************************************************************************************
   SUBROUTINE compute_chi(bs_env, gwkp_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT), TARGET             :: gwkp_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_chi'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: chi_grid_local_q(:, :, :), &
                                                            chi_grid_q(:, :, :), chi_q(:, :, :), &
                                                            G_grid_k(:, :, :), work(:, :)
      COMPLEX(KIND=dp), ALLOCATABLE, TARGET              :: psi_time(:, :, :, :), Z_group(:, :, :)
      COMPLEX(KIND=dp), POINTER                          :: psi_values(:, :, :, :), Z_time(:, :, :)
      INTEGER :: first_pair, handle, i_group, i_needed, i_q_world, i_rank, i_root, i_size, &
         i_source, i_t, i_tag, i_time_group, i_time_local, iq, iq_local, ispin, n_aux, n_cells, &
         n_grid, n_groups, n_local, n_needed_local, n_pairs_local, n_q_local, n_threads, &
         n_time_local, n_world
      INTEGER, ALLOCATABLE                               :: full_cells(:, :), needed_local(:), &
                                                            q_partner(:)
      LOGICAL, ALLOCATABLE                               :: needed(:)
      REAL(KIND=dp) :: available_GB, budget, input_bytes, live_bytes, n_pair_max, n_q_max, peak, &
         q_bytes, rows_bytes, step_start, tau, time_bytes
      REAL(KIND=dp), ALLOCATABLE                         :: chi_grid_R(:, :, :), G_occ_R(:, :, :), &
                                                            qpoints(:, :), step_times(:)
      TYPE(gw_pair_fft_work_type)                        :: fft_stencil_work
      TYPE(mp_comm_type)                                 :: column_env
      TYPE(mp_para_env_type), POINTER                    :: chi_env
      TYPE(ri_q_matrix_block_type), ALLOCATABLE          :: chi_partial(:)

      CALL timeset(routineN, handle)

      CPASSERT(gwkp_env%ngroups == 1 .AND. gwkp_env%my_group == 0)
      CALL build_full_cell_grid(bs_env, full_cells)
      n_grid = bs_env%ri_rs%n_grid_points
      n_aux = bs_env%n_RI
      n_world = gwkp_env%group_env%num_pe
      ALLOCATE (qpoints(3, bs_env%nkp_scf_desymm))
      qpoints(:, :) = bs_env%kpoints_chi_eps_W%xkp(:, 1:bs_env%nkp_scf_desymm)
      ! χ(-q) = χ(q)*: compute_W0 solves the first q of every ±q pair and conjugates W for -q.
      CALL minus_q_partners(qpoints, q_partner)
      ALLOCATE (needed(bs_env%nkp_scf_desymm))
      needed(:) = [(q_partner(iq) >= iq, iq=1, bs_env%nkp_scf_desymm)]
      CALL mp_mem_avail_per_rank_GB(gwkp_env%group_env, available_GB)
      budget = budget_after_reserve(available_GB)
      n_cells = SIZE(full_cells, 2)
      n_threads = 1
!$    n_threads = omp_get_max_threads()
      DO i_size = MIN(2, MAX(1, n_world/2)), n_world
         IF (MODULO(n_world, i_size) /= 0) CYCLE
         n_groups = n_world/i_size
         IF (n_groups > bs_env%num_time_freq_points) CYCLE
         i_rank = MODULO(gwkp_env%group_env%mepos, i_size)
         i_group = gwkp_env%group_env%mepos/i_size
         n_time_local = (bs_env%num_time_freq_points + n_groups - 1 - i_group)/n_groups
         n_pair_max = REAL(INT(n_grid, int_8)**2*(i_rank + 1)/i_size - INT(n_grid, int_8)**2*i_rank/i_size, dp)
         n_q_max = REAL(COUNT([(MODULO(iq, i_size) == i_rank .AND. needed(iq), &
                                iq=1, bs_env%nkp_scf_desymm)]), dp)
         input_bytes = 0.0_dp
         IF (n_groups > 1) THEN
            input_bytes = complex_bytes*REAL(n_grid, dp)*REAL(n_aux, dp)*n_q_max
            input_bytes = input_bytes + &
                          complex_bytes*REAL(n_grid, dp)*REAL(bs_env%n_ao*bs_env%n_spin, dp)* &
                          REAL((bs_env%nkp_scf_desymm + i_size - 1 - i_rank)/i_size, dp)
         END IF
         time_bytes = complex_bytes*REAL(n_aux, dp)**2*n_q_max*REAL(n_time_local, dp)
         ! Keep only this group's time slices, not a full frequency set per group.
         ! The shared lattice backend uses bounded matrix tiles; library_bytes reserves its scratch.
         rows_bytes = n_pair_max*complex_bytes*REAL(bs_env%nkp_scf_desymm + n_cells, dp)
         q_bytes = complex_bytes*REAL(n_grid, dp)**2*n_q_max
         ! G(k) and the pair transfers hold two full grid matrices.
         live_bytes = MAX(rows_bytes + 2.0_dp*complex_bytes*REAL(n_grid, dp)**2 + &
                          complex_bytes*REAL(n_grid, dp)*REAL(bs_env%n_ao, dp), &
                          complex_bytes*n_pair_max*REAL(bs_env%nkp_scf_desymm, dp) + q_bytes + &
                          2.0_dp*complex_bytes*REAL(n_grid, dp)**2, &
                          q_bytes + complex_bytes*REAL(n_grid, dp)*REAL(n_aux, dp)* &
                          MAX(1.0_dp, MIN(REAL(n_threads, dp), n_q_max)))
         peak = input_bytes + time_bytes + complex_bytes*REAL(n_aux, dp)**2*n_q_max + &
                MAX(live_bytes, complex_bytes*REAL(n_aux, dp)**2*REAL(bs_env%num_time_freq_points, dp)) + &
                2.0_dp*complex_bytes*REAL(n_threads*bs_env%nkp_scf_desymm, dp) + library_bytes
         CALL gwkp_env%group_env%max(peak)
         IF (peak <= budget) EXIT
      END DO
      IF (i_size > n_world) THEN
         CPABORT("Insufficient memory per MPI rank for GW time groups")
      END IF
      n_groups = n_world/i_size
      i_group = gwkp_env%group_env%mepos/i_size
      ALLOCATE (chi_env)
      CALL chi_env%from_split(gwkp_env%group_env, i_group)
      CALL column_env%from_split(gwkp_env%group_env, chi_env%mepos)
      IF (n_groups > 1) THEN
         CALL gather_group_columns(bs_env, gwkp_env, column_env, i_size, chi_env%mepos, Z_group, psi_time)
         psi_values => psi_time
         Z_time => Z_group
      ELSE
         psi_values => gwkp_env%psi_lm_k
         Z_time => gwkp_env%Z_lP_q
      END IF
      n_time_local = (bs_env%num_time_freq_points + n_groups - 1 - i_group)/n_groups
      CALL block_range(n_grid**2, chi_env%mepos, chi_env%num_pe, first_pair, n_pairs_local)
      n_q_local = COUNT([(MODULO(iq, chi_env%num_pe) == chi_env%mepos, iq=1, bs_env%nkp_scf_desymm)])
      ALLOCATE (chi_q(n_aux, n_aux, n_q_local), needed_local(n_q_local))
      ALLOCATE (chi_partial(n_q_local))
      iq_local = 0
      n_needed_local = 0
      DO iq = MODULO(chi_env%mepos - 1, i_size) + 1, SIZE(qpoints, 2), i_size
         iq_local = iq_local + 1
         IF (.NOT. needed(iq)) CYCLE
         n_needed_local = n_needed_local + 1
         needed_local(n_needed_local) = iq_local
         ALLOCATE (chi_partial(iq_local)%values(n_aux, n_aux, n_time_local), SOURCE=z_zero)
      END DO

      ALLOCATE (step_times(bs_env%num_time_freq_points), SOURCE=0.0_dp)
      DO i_t = i_group + 1, bs_env%num_time_freq_points, n_groups
         step_start = m_walltime()
         tau = bs_env%time_frequency_grid%imaginary_time(i_t)
         chi_q = z_zero
         DO ispin = 1, bs_env%n_spin
            ALLOCATE (G_grid_k(n_pairs_local, 1, bs_env%nkp_scf_desymm))
            ALLOCATE (G_occ_R(n_pairs_local, 1, n_cells))
            ALLOCATE (chi_grid_R, MOLD=G_occ_R)
            CALL build_green_pairs(bs_env, psi_values(:, :, :, ispin), chi_env, tau, ispin, .TRUE., G_grid_k)
            CALL k_to_real_cell_transform(G_grid_k, bs_env%kpoints_scf_desymm%xkp, &
                                          bs_env%kpoints_scf_desymm%nkp_grid, full_cells, &
                                          G_occ_R, rskp_backend_fft, profile_fft=.TRUE.)
            CALL build_green_pairs(bs_env, psi_values(:, :, :, ispin), chi_env, tau, ispin, .FALSE., G_grid_k)
            CALL k_to_real_cell_transform(G_grid_k, bs_env%kpoints_scf_desymm%xkp, &
                                          bs_env%kpoints_scf_desymm%nkp_grid, full_cells, &
                                          chi_grid_R, rskp_backend_fft, profile_fft=.TRUE.)
            CALL filter_full_green_cells(bs_env, chi_env, G_occ_R)
            CALL filter_full_green_cells(bs_env, chi_env, chi_grid_R)
            DEALLOCATE (G_grid_k)
            chi_grid_R(:, :, :) = bs_env%spin_degeneracy*G_occ_R*chi_grid_R
            DEALLOCATE (G_occ_R)
            ALLOCATE (chi_grid_q(n_pairs_local, 1, bs_env%nkp_scf_desymm))
            CALL product_to_k_pair_centered(bs_env, full_cells, qpoints, chi_grid_R, &
                                            chi_grid_q, first_pair, fft_work=fft_stencil_work)
            DEALLOCATE (chi_grid_R)

            ALLOCATE (chi_grid_local_q(n_grid, n_grid, n_needed_local))
            CALL fft_pairs_to_q_matrices(chi_grid_q, chi_env, needed, chi_grid_local_q)
            DEALLOCATE (chi_grid_q)
            ! χ_PQ(q) = Σ_ll′ Z_lP(q)* χ_ll′(q) Z_l′Q(q).
            ! Parallelize distinct q matrices; leave singleton GEMM threading to BLAS.
            !$OMP PARALLEL DEFAULT(NONE) IF(n_needed_local > 1) &
            !$OMP SHARED(n_grid, n_aux, n_needed_local, needed_local, chi_grid_local_q, Z_time, chi_q) &
            !$OMP PRIVATE(work, i_needed, iq_local)
            ALLOCATE (work(n_grid, n_aux))
            !$OMP DO SCHEDULE(STATIC)
            DO i_needed = 1, n_needed_local
               iq_local = needed_local(i_needed)
               CALL zgemm('N', 'N', n_grid, n_aux, n_grid, z_one, &
                          chi_grid_local_q(:, :, i_needed), n_grid, &
                          Z_time(:, :, iq_local), n_grid, z_zero, work, n_grid)
               CALL zgemm('C', 'N', n_aux, n_aux, n_grid, z_one, &
                          Z_time(:, :, iq_local), n_grid, work, n_grid, z_one, &
                          chi_q(:, :, iq_local), n_aux)
            END DO
            !$OMP END DO
            DEALLOCATE (work)
            !$OMP END PARALLEL
            DEALLOCATE (chi_grid_local_q)
         END DO
         i_time_local = (i_t - 1)/n_groups + 1
         DO i_needed = 1, n_needed_local
            iq_local = needed_local(i_needed)
            chi_partial(iq_local)%values(:, :, i_time_local) = chi_q(:, :, iq_local)
         END DO
         IF (chi_env%mepos == 0) step_times(i_t) = m_walltime() - step_start
      END DO

      CALL bs_env%para_env%sum(step_times)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(A)') ' '
         DO i_t = 1, bs_env%num_time_freq_points
            WRITE (bs_env%unit_nr, '(T2,A,I3,A,I0,T57,A,T71,F10.1,A)') &
               'GW χ(iτ) point ', i_t, ' / ', bs_env%num_time_freq_points, &
               'Exec. time:', step_times(i_t), ' s'
         END DO
         CALL m_flush(bs_env%unit_nr)
      END IF
      DEALLOCATE (step_times)

      IF (ALLOCATED(psi_time)) DEALLOCATE (psi_time)
      n_local = COUNT([(MODULO(iq, n_world) == gwkp_env%group_env%mepos, iq=1, bs_env%nkp_scf_desymm)])
      ALLOCATE (gwkp_env%chi_q_frequency(n_local))
      iq_local = 0
      DO iq = MODULO(chi_env%mepos - 1, i_size) + 1, bs_env%nkp_scf_desymm, i_size
         iq_local = iq_local + 1
         IF (.NOT. needed(iq)) CYCLE
         i_root = MODULO(iq, n_world)/i_size
         i_q_world = (iq - 1)/n_world + 1
         IF (column_env%mepos == i_root) THEN
            ALLOCATE (gwkp_env%chi_q_frequency(i_q_world)%values(n_aux, n_aux, &
                                                                 bs_env%num_time_freq_points))
         END IF
         ! Send each computed time slice once, point to point, to the world q owner, which alone
         ! forms the complete frequency result; slices and frequency sets are never replicated.
         DO i_t = 1, bs_env%num_time_freq_points
            i_time_group = MODULO(i_t - 1, n_groups)
            i_time_local = (i_t - 1)/n_groups + 1
            IF (column_env%mepos == i_root) THEN
               IF (i_time_group == i_root) THEN
                  gwkp_env%chi_q_frequency(i_q_world)%values(:, :, i_t) = &
                     chi_partial(iq_local)%values(:, :, i_time_local)
               ELSE
                  i_source = i_time_group
                  i_tag = 0
                  CALL column_env%recv(gwkp_env%chi_q_frequency(i_q_world)%values(:, :, i_t), i_source, i_tag)
               END IF
            ELSE IF (column_env%mepos == i_time_group) THEN
               CALL column_env%send(chi_partial(iq_local)%values(:, :, i_time_local), i_root, 0)
            END IF
         END DO
         IF (column_env%mepos == i_root) THEN
            CALL cosine_transform_ri_q(gwkp_env%chi_q_frequency(i_q_world), &
                                       bs_env%time_frequency_grid%cosine_time_to_frequency)
         END IF
         DEALLOCATE (chi_partial(iq_local)%values)
      END DO
      DEALLOCATE (chi_partial)
      NULLIFY (Z_time)
      IF (ALLOCATED(Z_group)) DEALLOCATE (Z_group)
      CALL column_env%free()
      CALL mp_para_env_release(chi_env)
      DEALLOCATE (chi_q, full_cells, needed, needed_local, q_partner, qpoints)

      CALL timestop(handle)

   END SUBROUTINE compute_chi

! **************************************************************************************************
!> \brief Build K(q) = M(q)^-1 V(q)^1/2 from the Gaussian RI metric M and the Ewald Coulomb matrix V,
!>        and the charge-neutralisation vectors e(q) = M(q) c/|c| and X(q,iω) = M(q) M(0)^-1 X0(iω).
!> The Gamma-cell head is corrected separately from these operators.
!> \param bs_env ...
!> \param qs_env ...
!> \param gwkp_env K(q) and the neutralisation vectors for the q points owned by this rank
! **************************************************************************************************
   SUBROUTINE compute_metric_coulomb_factor(bs_env, qs_env, gwkp_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gwkp_type), INTENT(INOUT)                     :: gwkp_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_metric_coulomb_factor'

      COMPLEX(KIND=dp), ALLOCATABLE :: charge(:), charge_hat(:), M_inv(:, :), M_ref(:, :, :), &
         M_zero(:, :), U_Pp(:, :), V_ref(:, :, :), V_sqrt(:, :), X0(:, :), X_coef(:, :)
      INTEGER                                            :: handle, ikp, ilocal, n_ref, n_RI, nkp, &
                                                            nkp_local
      INTEGER, ALLOCATABLE                               :: metric_cells(:, :)
      TYPE(cp_fm_type), ALLOCATABLE                      :: fm_metric(:)

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI
      n_ref = n_RI
      IF (bs_env%auto_ri%enabled) n_ref = SUM(bs_env%auto_ri%sizes_ref_RI)
      nkp = bs_env%nkp_scf_desymm
      CPASSERT(bs_env%nkp_chi_eps_W_orig == nkp)
      nkp_local = COUNT([(MODULO(ikp, bs_env%para_env%num_pe) == &
                          bs_env%para_env%mepos, ikp=1, nkp)])
      ALLOCATE (M_ref(n_ref, n_ref, nkp_local), V_ref(n_ref, n_ref, nkp_local), source=z_zero)
      ALLOCATE (gwkp_env%metric_coulomb_factor(n_RI, n_RI, nkp_local), M_inv(n_RI, n_RI), V_sqrt(n_RI, n_RI))
      CALL build_ri_2c_lattice(fm_metric, metric_cells, bs_env%ri_metric, &
                               0.0_dp, bs_env, qs_env, reference_basis=bs_env%auto_ri%enabled)
      CALL ri_2c_lattice_to_owned_q(fm_metric, metric_cells, &
                                    bs_env%kpoints_chi_eps_W%xkp(1:3, 1:nkp), &
                                    bs_env%kpoints_chi_eps_W%nkp_grid, M_ref, &
                                    bs_env%kpoints_chi_eps_W%rskp_transform)
      ! The erfc(omega r)/r metric is the real-space Ewald part of the Coulomb matrix, which is
      ! finite on the shifted 2D/3D q mesh. No G=0 auxiliary channel is removed at a transfer point.
      V_ref(:, :, :) = M_ref
      CALL add_ewald_reciprocal(V_ref, qs_env, bs_env%kpoints_chi_eps_W, 'RI_AUX', &
                                bs_env%ri_metric%omega, 1, nkp)
      ALLOCATE (M_zero(n_ref, n_ref), charge(n_ref), charge_hat(n_RI))
      CALL ri_2c_lattice_to_kp(fm_metric, metric_cells, [0.0_dp, 0.0_dp, 0.0_dp], M_zero)
      CALL gamma_auxiliary_charges(qs_env, charge)
      CALL cp_fm_release(fm_metric)
      DEALLOCATE (metric_cells)
      IF (bs_env%auto_ri%enabled) CALL assemble_AA_contraction(bs_env, U_Pp)
      IF (bs_env%auto_ri%enabled) THEN
         M_inv(:, :) = MATMUL(CONJG(TRANSPOSE(U_Pp)), MATMUL(M_zero, U_Pp))
         charge_hat(:) = MATMUL(CONJG(TRANSPOSE(U_Pp)), charge)
      ELSE
         M_inv(:, :) = M_zero
         charge_hat(:) = charge
      END IF
      DEALLOCATE (M_zero, charge)
      charge_hat(:) = charge_hat/SQRT(REAL(DOT_PRODUCT(charge_hat, charge_hat), dp))
      CALL local_complex_power(M_inv, -1.0_dp, 0.0_dp, relative_eps=bs_env%regularization_RI)
      ! The pair densities enter χ as Z^H ρ = M d. Neutralise the fitted coefficients d, which are
      ! continuous in q: X0 = Σ (Z(0)^H ρ) m* becomes M(0)^-1 X0 = Σ d m*, and M(q) M(0)^-1 X0 at q.
      ! Carrying X0 itself to q would let M(q)^-1 amplify its O(q) mismatch in ill-conditioned metrics.
      CALL compute_charge_leak_q0(bs_env, gwkp_env, MATMUL(M_inv, charge_hat), X0)
      ALLOCATE (X_coef(n_RI, SIZE(X0, 2)))
      X_coef(:, :) = MATMUL(M_inv, X0)
      DEALLOCATE (X0)
      ALLOCATE (gwkp_env%leak_e(n_RI, nkp_local), gwkp_env%leak_x(n_RI, SIZE(X_coef, 2), nkp_local))
      DO ilocal = 1, nkp_local
         IF (bs_env%auto_ri%enabled) THEN
            M_inv(:, :) = MATMUL(CONJG(TRANSPOSE(U_Pp)), MATMUL(M_ref(:, :, ilocal), U_Pp))
            V_sqrt(:, :) = MATMUL(CONJG(TRANSPOSE(U_Pp)), MATMUL(V_ref(:, :, ilocal), U_Pp))
         ELSE
            M_inv(:, :) = M_ref(:, :, ilocal)
            V_sqrt(:, :) = V_ref(:, :, ilocal)
         END IF
         gwkp_env%leak_e(:, ilocal) = MATMUL(M_inv, charge_hat)
         gwkp_env%leak_x(:, :, ilocal) = MATMUL(M_inv, X_coef)
         CALL local_complex_power(M_inv, -1.0_dp, 0.0_dp, relative_eps=bs_env%regularization_RI)
         CALL local_complex_power(V_sqrt, 0.5_dp, 0.0_dp)
         gwkp_env%metric_coulomb_factor(:, :, ilocal) = MATMUL(M_inv, V_sqrt)
      END DO
      DEALLOCATE (M_ref, V_ref, M_inv, V_sqrt, charge_hat, X_coef)
      IF (ALLOCATED(U_Pp)) DEALLOCATE (U_Pp)

      CALL timestop(handle)

   END SUBROUTINE compute_metric_coulomb_factor

! **************************************************************************************************
!> \brief Spurious charge of the fitted q = 0 occupied-virtual pair densities.
!> The RI fit gives ρ_ia = ψ_ik* ψ_ak (exactly neutral) the charge m_ia(k) = Σ_l g_l* ψ_ik(r_l) ψ_ak(r_l)*
!> with g = Z(0) v0. With the time dependence and cosine transform of compute_chi,
!> L0(iω) = (s/N_k) Σ_k,ia |m_ia|² w_ia(iω) and X0(iω) = Z(0)^H T(iω), T_l = (s/N_k) Σ_k,ia ψ_i ψ_a* m_ia* w_ia.
!> \param bs_env ...
!> \param gwkp_env ...
!> \param v0 M(0)^-1 c/|c|
!> \param X0 X(iω) in the projection space Z(0)^H ρ
! **************************************************************************************************
   SUBROUTINE compute_charge_leak_q0(bs_env, gwkp_env, v0, X0)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT)                     :: gwkp_env
      COMPLEX(KIND=dp), INTENT(IN)                       :: v0(:)
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: X0(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_charge_leak_q0'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: cross(:, :), g0(:), m(:, :), psi(:, :), &
                                                            T_tau(:, :), work(:, :), X_tau(:, :), &
                                                            Y(:, :)
      INTEGER                                            :: a, handle, i, ikp, ispin, j_w, n_ao, &
                                                            n_grid, n_occ, n_t, n_vir, nkp, &
                                                            num_pe, t
      REAL(KIND=dp)                                      :: prefactor, x
      REAL(KIND=dp), ALLOCATABLE                         :: L_tau(:)

      CALL timeset(routineN, handle)

      n_grid = bs_env%ri_rs%n_grid_points
      n_ao = SIZE(gwkp_env%psi_lm_k, 2)
      nkp = bs_env%nkp_scf_desymm
      n_t = bs_env%num_time_freq_points
      num_pe = bs_env%para_env%num_pe
      ALLOCATE (g0(n_grid), psi(n_grid, n_ao), work(n_grid, n_ao), T_tau(n_grid, n_t), L_tau(n_t))
      g0(:) = MATMUL(gwkp_env%Z_lP_zero(:, 2:, 1), v0)
      T_tau(:, :) = z_zero
      L_tau(:) = 0.0_dp
      DO ispin = 1, bs_env%n_spin
         n_occ = bs_env%n_occ(ispin)
         n_vir = n_ao - n_occ
         IF (n_occ == 0 .OR. n_vir == 0) CYCLE
         ALLOCATE (m(n_occ, n_vir), cross(n_occ, n_vir), Y(n_grid, n_occ))
         DO ikp = 1, nkp
            IF (MODULO(ikp - 1, num_pe) /= bs_env%para_env%mepos) CYCLE
            psi(:, :) = gwkp_env%psi_lm_k(:, :, (ikp - 1)/num_pe + 1, ispin)
            DO i = 1, n_occ
               work(:, i) = CONJG(g0)*psi(:, i)
            END DO
            m(:, :) = MATMUL(TRANSPOSE(work(:, 1:n_occ)), CONJG(psi(:, n_occ + 1:n_ao)))
            DO t = 1, n_t
               cross(:, :) = z_zero
               DO a = 1, n_vir
                  DO i = 1, n_occ
                     x = bs_env%time_frequency_grid%imaginary_time(t)* &
                         (gwkp_env%eigenvalues_G(n_occ + a, ikp, ispin) - gwkp_env%eigenvalues_G(i, ikp, ispin))
                     IF (x < 700.0_dp) THEN
                        cross(i, a) = CONJG(m(i, a))*EXP(-x)
                        L_tau(t) = L_tau(t) + ABS(m(i, a))**2*EXP(-x)
                     END IF
                  END DO
               END DO
               Y(:, :) = MATMUL(CONJG(psi(:, n_occ + 1:n_ao)), TRANSPOSE(cross))
               T_tau(:, t) = T_tau(:, t) + SUM(psi(:, 1:n_occ)*Y, DIM=2)
            END DO
         END DO
         DEALLOCATE (m, cross, Y)
      END DO
      CALL bs_env%para_env%sum(T_tau)
      CALL bs_env%para_env%sum(L_tau)
      prefactor = bs_env%spin_degeneracy/REAL(nkp, dp)
      ALLOCATE (X_tau(bs_env%n_RI, n_t))
      X_tau(:, :) = prefactor*MATMUL(CONJG(TRANSPOSE(gwkp_env%Z_lP_zero(:, 2:, 1))), T_tau)
      ALLOCATE (X0(bs_env%n_RI, n_t), SOURCE=z_zero)
      ALLOCATE (gwkp_env%leak_l0(n_t), SOURCE=0.0_dp)
      DO t = 1, n_t
         DO j_w = 1, n_t
            X0(:, j_w) = X0(:, j_w) + bs_env%time_frequency_grid%cosine_time_to_frequency(j_w, t)*X_tau(:, t)
            gwkp_env%leak_l0(j_w) = gwkp_env%leak_l0(j_w) + &
                                    bs_env%time_frequency_grid%cosine_time_to_frequency(j_w, t)*prefactor*L_tau(t)
         END DO
      END DO
      DEALLOCATE (g0, psi, work, L_tau, T_tau, X_tau)

      CALL timestop(handle)

   END SUBROUTINE compute_charge_leak_q0

! **************************************************************************************************
!> \brief Contract the occupied density with the bare RI Coulomb operator.
!> \param bs_env Band-structure settings and physical results.
!> \param gwkp_env Periodic GW work data.
! **************************************************************************************************
   SUBROUTINE compute_Sigma_x_from_V(bs_env, gwkp_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT), TARGET             :: gwkp_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Sigma_x_from_V'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: G_grid_k(:, :, :), head_cell_factor(:), &
                                                            q_matrices(:, :, :), ZK(:, :)
      INTEGER                                            :: first_pair, handle, iq, iq_local, ispin, &
                                                            n_grid, n_pairs_local, n_q_local, n_RI
      INTEGER, ALLOCATABLE                               :: full_cells(:, :)
      REAL(KIND=dp)                                      :: head_metric_zero
      REAL(KIND=dp), ALLOCATABLE                         :: qpoints(:, :), Sigma_grid_R(:, :, :), &
                                                            V_grid_R(:, :, :)

      CALL timeset(routineN, handle)

      ALLOCATE (gwkp_env%sigma_x_bands(bs_env%n_ao, bs_env%nkp_bs_and_DOS, bs_env%n_spin))
      n_grid = bs_env%ri_rs%n_grid_points
      ! Same grid-pair distribution as Σc with q owners MOD(q,P).
      CALL block_range(n_grid**2, bs_env%para_env%mepos, bs_env%para_env%num_pe, first_pair, n_pairs_local)
      n_q_local = SIZE(gwkp_env%Z_lP_q, 3)
      IF (SIZE(gwkp_env%Z_lP_zero, 2) /= bs_env%n_RI + 1 .OR. SIZE(gwkp_env%Z_lP_zero, 3) /= 1) THEN
         CPABORT("Periodic RI-RS GW requires the direct Gamma auxiliary fit")
      END IF
      ALLOCATE (qpoints(3, bs_env%nkp_scf_desymm))
      qpoints(:, :) = bs_env%kpoints_chi_eps_W%xkp(:, 1:bs_env%nkp_scf_desymm)
      CPASSERT(SIZE(gwkp_env%Z_lP_q, 2) == bs_env%n_RI)
      CPASSERT(ALLOCATED(gwkp_env%metric_coulomb_factor))
      n_RI = bs_env%n_RI
      ALLOCATE (q_matrices(n_grid, n_grid, n_q_local))
      ! V_ll'(q) = Z(q) K(q) K(q)† Z(q)† = (ZK)(ZK)† with the metric-Coulomb factor K.
      !$OMP PARALLEL DEFAULT(NONE) SHARED(gwkp_env, q_matrices, n_q_local, n_grid, n_RI) &
      !$OMP PRIVATE(iq_local, ZK)
      ALLOCATE (ZK(n_grid, n_RI))
      !$OMP DO SCHEDULE(STATIC)
      DO iq_local = 1, n_q_local
         CALL zgemm('N', 'N', n_grid, n_RI, n_RI, z_one, gwkp_env%Z_lP_q(:, :, iq_local), n_grid, &
                    gwkp_env%metric_coulomb_factor(:, :, iq_local), n_RI, z_zero, ZK, n_grid)
         CALL zgemm('N', 'C', n_grid, n_grid, n_RI, z_one, ZK, n_grid, ZK, n_grid, &
                    z_zero, q_matrices(:, :, iq_local), n_grid)
      END DO
      !$OMP END DO
      DEALLOCATE (ZK)
      !$OMP END PARALLEL
      DEALLOCATE (gwkp_env%metric_coulomb_factor)
      ALLOCATE (G_grid_k(n_pairs_local, 1, bs_env%nkp_scf_desymm))
      CALL q_matrices_to_fft_pairs(q_matrices, bs_env%para_env, G_grid_k)
      DEALLOCATE (q_matrices)
      CALL build_full_cell_grid(bs_env, full_cells)
      ALLOCATE (V_grid_R(n_pairs_local, 1, SIZE(full_cells, 2)))
      CALL k_to_real_cell_transform(G_grid_k, qpoints, bs_env%kpoints_scf_desymm%nkp_grid, &
                                    full_cells, V_grid_R, rskp_backend_fft, profile_fft=.TRUE.)
      CALL singular_head_cell_factors(bs_env, full_cells, head_cell_factor)
      head_metric_zero = ri_metric_head_zero(bs_env)
      DO iq = 1, SIZE(full_cells, 2)
         IF (ALL(bs_env%periodic == [1, 1, 0])) THEN
            CALL add_gamma_corr_pair_cell(bs_env, full_cells, iq, &
                                          twopi/head_metric_zero**2*head_cell_factor(iq), &
                                          gwkp_env%Z_lP_zero(:, 1, 1), first_pair, V_grid_R)
         ELSE
            CALL add_gamma_corr_pair_cell(bs_env, full_cells, iq, &
                                          fourpi/head_metric_zero**2*head_cell_factor(iq), &
                                          gwkp_env%Z_lP_zero(:, 1, 1), first_pair, V_grid_R)
         END IF
      END DO
      DEALLOCATE (head_cell_factor)

      ALLOCATE (Sigma_grid_R, MOLD=V_grid_R)
      DO ispin = 1, bs_env%n_spin
         CALL build_green_pairs(bs_env, gwkp_env%band_orbitals(:, :, :, ispin), bs_env%para_env, &
                                0.0_dp, ispin, .TRUE., G_grid_k)
         CALL k_to_real_cell_transform(G_grid_k, bs_env%kpoints_scf_desymm%xkp, &
                                       bs_env%kpoints_scf_desymm%nkp_grid, full_cells, Sigma_grid_R, &
                                       rskp_backend_fft, profile_fft=.TRUE.)
         Sigma_grid_R(:, :, :) = -Sigma_grid_R*V_grid_R
         CALL self_energy_on_mesh(bs_env, full_cells, Sigma_grid_R, first_pair, G_grid_k)

         CALL self_energy_pairs_to_Bloch(bs_env, G_grid_k, Sigma_grid_R, first_pair, bs_env%para_env, &
                                         gwkp_env%band_orbitals(:, :, :, ispin), &
                                         gwkp_env%sigma_x_bands(:, :, ispin))
      END DO

      DEALLOCATE (full_cells, G_grid_k, qpoints, Sigma_grid_R, V_grid_R)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x_from_V

! **************************************************************************************************
!> \brief Form Σc with one owner per ordered (l,l') pair through both lattice FFTs.
!> \param bs_env Band-structure settings and physical results.
!> \param gwkp_env Periodic GW work data.
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c_from_W_PQ(bs_env, gwkp_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT), TARGET             :: gwkp_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Sigma_c_from_W_PQ'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: G_grid_k(:, :, :)
      INTEGER                                            :: first_pair, handle, i_group, i_t, ispin, &
                                                            n_grid, n_pairs_local
      INTEGER, ALLOCATABLE                               :: full_cells(:, :)
      REAL(KIND=dp)                                      :: step_elapsed, step_start, tau
      REAL(KIND=dp), ALLOCATABLE :: qpoints(:, :), sigma_band(:, :), Sigma_grid_R(:, :, :), &
         sigma_result_buffer(:), step_times(:), W_grid_R(:, :, :)
      TYPE(gw_pair_fft_work_type)                        :: fft_stencil_work

      CALL timeset(routineN, handle)

      CPASSERT(ALLOCATED(gwkp_env%band_orbitals))
      IF (ALLOCATED(gwkp_env%sigma_c_bands)) DEALLOCATE (gwkp_env%sigma_c_bands)
      CALL build_full_cell_grid(bs_env, full_cells)
      n_grid = bs_env%ri_rs%n_grid_points
      CALL block_range(n_grid**2, gwkp_env%group_env%mepos, gwkp_env%group_env%num_pe, &
                       first_pair, n_pairs_local)
      ALLOCATE (G_grid_k(n_pairs_local, 1, bs_env%nkp_scf_desymm))
      ALLOCATE (W_grid_R(n_pairs_local, 1, SIZE(full_cells, 2)))
      ALLOCATE (Sigma_grid_R, MOLD=W_grid_R)
      ALLOCATE (qpoints(3, bs_env%nkp_scf_desymm))
      qpoints(:, :) = bs_env%kpoints_chi_eps_W%xkp(:, 1:bs_env%nkp_scf_desymm)
      ALLOCATE (sigma_band(bs_env%n_ao, bs_env%nkp_bs_and_DOS))
      ALLOCATE (gwkp_env%sigma_c_bands(bs_env%n_ao, bs_env%nkp_bs_and_DOS, &
                                       bs_env%num_time_freq_points, 2, bs_env%n_spin))
      ALLOCATE (step_times(bs_env%num_time_freq_points), SOURCE=0.0_dp)
      gwkp_env%sigma_c_bands = 0.0_dp

      DO i_t = gwkp_env%my_group + 1, bs_env%num_time_freq_points, gwkp_env%ngroups
         step_start = m_walltime()
         i_group = (i_t - gwkp_env%my_group - 1)/gwkp_env%ngroups + 1
         CALL project_auxiliary_time_q_to_fft_pairs(gwkp_env%Z_lP_q, &
                                                    gwkp_env%W_PQ_tau, i_group, &
                                                    gwkp_env%group_env, G_grid_k)

         CALL k_to_real_cell_transform(G_grid_k, qpoints, &
                                       bs_env%kpoints_scf_desymm%nkp_grid, full_cells, &
                                       W_grid_R, rskp_backend_fft, profile_fft=.TRUE.)
         CALL add_gamma_corr_cell_correction(bs_env, full_cells, gwkp_env%gamma_correction, i_t, &
                                             gwkp_env%Z_lP_zero(:, 1, 1), first_pair, W_grid_R)
         CALL truncate_grid_W_by_pair_distance_pairs(bs_env, full_cells, first_pair, W_grid_R)
         tau = bs_env%time_frequency_grid%imaginary_time(i_t)
         DO ispin = 1, bs_env%n_spin
            CALL build_green_pairs(bs_env, gwkp_env%band_orbitals(:, :, :, ispin), gwkp_env%group_env, &
                                   tau, ispin, .TRUE., G_grid_k)
            CALL k_to_real_cell_transform(G_grid_k, bs_env%kpoints_scf_desymm%xkp, &
                                          bs_env%kpoints_scf_desymm%nkp_grid, full_cells, Sigma_grid_R, &
                                          rskp_backend_fft, profile_fft=.TRUE.)
            CALL filter_full_green_cells(bs_env, gwkp_env%group_env, Sigma_grid_R)
            Sigma_grid_R(:, :, :) = -Sigma_grid_R*W_grid_R
            CALL self_energy_on_mesh(bs_env, full_cells, Sigma_grid_R, first_pair, G_grid_k, &
                                     fft_work=fft_stencil_work)
            CALL self_energy_pairs_to_Bloch(bs_env, G_grid_k, Sigma_grid_R, first_pair, gwkp_env%group_env, &
                                            gwkp_env%band_orbitals(:, :, :, ispin), sigma_band)
            IF (gwkp_env%group_env%mepos == 0) gwkp_env%sigma_c_bands(:, :, i_t, 2, ispin) = sigma_band

            CALL build_green_pairs(bs_env, gwkp_env%band_orbitals(:, :, :, ispin), gwkp_env%group_env, &
                                   tau, ispin, .FALSE., G_grid_k)
            CALL k_to_real_cell_transform(G_grid_k, bs_env%kpoints_scf_desymm%xkp, &
                                          bs_env%kpoints_scf_desymm%nkp_grid, full_cells, Sigma_grid_R, &
                                          rskp_backend_fft, profile_fft=.TRUE.)
            CALL filter_full_green_cells(bs_env, gwkp_env%group_env, Sigma_grid_R)
            Sigma_grid_R(:, :, :) = Sigma_grid_R*W_grid_R
            CALL self_energy_on_mesh(bs_env, full_cells, Sigma_grid_R, first_pair, G_grid_k, &
                                     fft_work=fft_stencil_work)
            CALL self_energy_pairs_to_Bloch(bs_env, G_grid_k, Sigma_grid_R, first_pair, gwkp_env%group_env, &
                                            gwkp_env%band_orbitals(:, :, :, ispin), sigma_band)
            IF (gwkp_env%group_env%mepos == 0) gwkp_env%sigma_c_bands(:, :, i_t, 1, ispin) = sigma_band
         END DO
         step_elapsed = m_walltime() - step_start
         IF (gwkp_env%group_env%mepos == 0) step_times(i_t) = step_elapsed
      END DO

      CALL bs_env%para_env%sum(step_times)
      IF (bs_env%unit_nr > 0) THEN
         DO i_t = 1, bs_env%num_time_freq_points
            WRITE (bs_env%unit_nr, '(T2,A,I3,A,I0,T57,A,T71,F10.1,A)') &
               'GW Σ_c(iτ) point ', i_t, ' / ', bs_env%num_time_freq_points, &
               'Exec. time:', step_times(i_t), ' s'
         END DO
         CALL m_flush(bs_env%unit_nr)
      END IF

      ALLOCATE (sigma_result_buffer(SIZE(gwkp_env%sigma_c_bands)))
      sigma_result_buffer(:) = RESHAPE(gwkp_env%sigma_c_bands, [SIZE(gwkp_env%sigma_c_bands)])
      CALL bs_env%para_env%sum(sigma_result_buffer)
      gwkp_env%sigma_c_bands(:, :, :, :, :) = RESHAPE(sigma_result_buffer, SHAPE(gwkp_env%sigma_c_bands))

      DEALLOCATE (full_cells, G_grid_k, qpoints, &
                  sigma_band, sigma_result_buffer, Sigma_grid_R, &
                  step_times, W_grid_R)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c_from_W_PQ

! **************************************************************************************************
!> \brief Drop an entire real-space Green matrix only when its full Frobenius norm,
!>        divided by the number of atom pairs, is below EPS_FILTER.
!> \param bs_env ...
!> \param para_env ...
!> \param G_R ...
! **************************************************************************************************
   SUBROUTINE filter_full_green_cells(bs_env, para_env, G_R)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), INTENT(INOUT)                       :: G_R(:, :, :)

      INTEGER                                            :: i_cell, n_cells
      REAL(KIND=dp)                                      :: atom_block_count
      REAL(KIND=dp), ALLOCATABLE                         :: norm_squared(:)

      CPASSERT(bs_env%n_atom > 0)
      n_cells = SIZE(G_R, 3)
      atom_block_count = REAL(bs_env%n_atom, dp)**2
      ALLOCATE (norm_squared(n_cells))
      DO i_cell = 1, n_cells
         norm_squared(i_cell) = SUM(G_R(:, :, i_cell)**2)
      END DO
      CALL para_env%sum(norm_squared)
      DO i_cell = 1, n_cells
         IF (SQRT(norm_squared(i_cell))/atom_block_count < bs_env%eps_filter) THEN
            G_R(:, :, i_cell) = 0.0_dp
         END IF
      END DO
      DEALLOCATE (norm_squared)
   END SUBROUTINE filter_full_green_cells

! **************************************************************************************************
!> \brief Apply the same real-space cutoff to ordered pairs owned by this rank.
!> \param bs_env ...
!> \param cells ...
!> \param first_pair First locally owned ordered grid pair.
!> \param W_grid_R ...
! **************************************************************************************************
   SUBROUTINE truncate_grid_W_by_pair_distance_pairs(bs_env, cells, first_pair, W_grid_R)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: cells
      INTEGER, INTENT(IN)                                :: first_pair
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: W_grid_R

      INTEGER                                            :: i_l, i_local, i_lp, i_pair, i_R, n_grid
      REAL(KIND=dp)                                      :: cell_index(3), cutoff_squared, &
                                                            displacement(3), lattice_vector(3)

      IF (bs_env%ri_rs%cutoff_radius_v_w <= 0.0_dp) RETURN
      CPASSERT(ALLOCATED(bs_env%ri_rs%grid_points))
      CPASSERT(SIZE(W_grid_R, 2) == 1)
      CPASSERT(SIZE(cells, 2) == SIZE(W_grid_R, 3))
      cutoff_squared = bs_env%ri_rs%cutoff_radius_v_w**2
      n_grid = bs_env%ri_rs%n_grid_points
      !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
      !$OMP SHARED(bs_env, cells, cutoff_squared, first_pair, n_grid, W_grid_R) &
      !$OMP PRIVATE(cell_index, displacement, i_l, i_local, i_lp, i_pair, i_R, lattice_vector)
      DO i_R = 1, SIZE(cells, 2)
         cell_index(:) = REAL(cells(:, i_R), dp)
         lattice_vector(:) = MATMUL(bs_env%hmat, cell_index)
         DO i_local = 1, SIZE(W_grid_R, 1)
            i_pair = first_pair + i_local - 1
            i_l = MODULO(i_pair - 1, n_grid) + 1
            i_lp = (i_pair - 1)/n_grid + 1
            displacement = bs_env%ri_rs%grid_points(:, i_l) - &
                           bs_env%ri_rs%grid_points(:, i_lp) - lattice_vector
            IF (SUM(displacement**2) > cutoff_squared) W_grid_R(i_local, 1, i_R) = 0.0_dp
         END DO
      END DO
      !$OMP END PARALLEL DO

   END SUBROUTINE truncate_grid_W_by_pair_distance_pairs

! **************************************************************************************************
!> \brief Form G_ll′(k,τ) at the k owners and send each ordered grid pair to its owner.
!> Pair i=l+(l′-1)*n_grid follows Fortran column order, so every rank's pair block is a
!> contiguous section of the full matrix. In each round every rank forms one native k point.
!> \param bs_env ...
!> \param orbitals Band Bloch orbitals of one spin at the k points owned in env, MOD(k-1,P).
!> \param env Ranks that share the grid pairs and the band k points.
!> \param tau ...
!> \param ispin ...
!> \param occupied ...
!> \param G_pairs ...
! **************************************************************************************************
   SUBROUTINE build_green_pairs(bs_env, orbitals, env, tau, ispin, occupied, G_pairs)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), INTENT(IN)                       :: orbitals(:, :, :)
      TYPE(mp_para_env_type), POINTER                    :: env
      REAL(KIND=dp), INTENT(IN)                          :: tau
      INTEGER, INTENT(IN)                                :: ispin
      LOGICAL, INTENT(IN)                                :: occupied
      COMPLEX(KIND=dp), INTENT(OUT)                      :: G_pairs(:, :, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_green_pairs'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: green(:), received(:), scaled_coeff(:, :)
      INTEGER                                            :: first_band, handle, i_band, i_col, &
                                                            i_first, i_rank, ikp, n_bands, n_grid, &
                                                            n_k, n_local, n_rank
      INTEGER, ALLOCATABLE                               :: pair_counts(:), pair_displs(:), &
                                                            recv_counts(:), recv_displs(:), &
                                                            send_counts(:)
      REAL(KIND=dp)                                      :: tau_E
      TYPE(gwkp_type), POINTER                           :: gwkp_env

      CALL timeset(routineN, handle)

      gwkp_env => bs_env%gwkp_env

      n_grid = bs_env%ri_rs%n_grid_points
      n_k = bs_env%nkp_scf_desymm
      n_local = SIZE(G_pairs, 1)
      n_rank = env%num_pe
      first_band = 1
      n_bands = bs_env%n_occ(ispin)
      IF (.NOT. occupied) THEN
         first_band = n_bands + 1
         n_bands = bs_env%n_ao - n_bands
      END IF
      CPASSERT(n_bands >= 0 .AND. n_bands <= bs_env%n_ao)
      IF (n_bands == 0) THEN
         G_pairs = z_zero
         CALL timestop(handle)
         RETURN
      END IF
      CALL grid_pair_blocks(n_grid, n_rank, pair_counts, pair_displs)
      ALLOCATE (green(n_grid**2), received(n_local*n_rank), scaled_coeff(n_grid, n_bands), &
                recv_counts(0:n_rank - 1), recv_displs(0:n_rank - 1), send_counts(0:n_rank - 1))
      recv_displs(:) = [(i_rank*n_local, i_rank=0, n_rank - 1)]
      ! G_ll′(k,τ) = Σ_n [ψ_ln(k) exp(-|τ(ε_nk-ε_F)|/2)]
      !                         [ψ_l′n(k) exp(-|τ(ε_nk-ε_F)|/2)]*.
      DO i_first = 1, n_k, n_rank
         ikp = i_first + env%mepos
         send_counts(:) = 0
         IF (ikp <= n_k) THEN
            scaled_coeff(:, :) = orbitals(:, first_band:first_band + n_bands - 1, (ikp - 1)/n_rank + 1)
            DO i_band = 1, n_bands
               tau_E = ABS(0.5_dp*tau*(gwkp_env%eigenvalues_G(first_band + i_band - 1, ikp, ispin) - &
                                       bs_env%e_fermi(ispin)))
               IF (tau_E < bs_env%stabilize_exp) THEN
                  scaled_coeff(:, i_band) = scaled_coeff(:, i_band)*EXP(-tau_E)
               ELSE
                  scaled_coeff(:, i_band) = z_zero
               END IF
            END DO
            !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
            !$OMP SHARED(green, n_bands, n_grid, scaled_coeff) PRIVATE(i_col)
            DO i_col = 1, n_grid, column_block_cols
               CALL zgemm('N', 'C', n_grid, MIN(column_block_cols, n_grid - i_col + 1), n_bands, z_one, &
                          scaled_coeff, n_grid, scaled_coeff(i_col, 1), n_grid, z_zero, &
                          green(1 + (i_col - 1)*n_grid), n_grid)
            END DO
            !$OMP END PARALLEL DO
            send_counts(:) = pair_counts
         END IF
         DO i_rank = 0, n_rank - 1
            recv_counts(i_rank) = MERGE(n_local, 0, i_first + i_rank <= n_k)
         END DO
         CALL env%alltoall(green, send_counts, pair_displs, received, recv_counts, recv_displs)
         DO i_rank = 0, MIN(n_rank, n_k - i_first + 1) - 1
            G_pairs(:, 1, i_first + i_rank) = received(i_rank*n_local + 1:(i_rank + 1)*n_local)
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE build_green_pairs

! **************************************************************************************************
!> \brief Ordered grid-pair blocks of all ranks, as assigned by block_range.
!> Rank r owns the pairs displs(r)+1, ..., displs(r)+counts(r).
!> \param n_grid ...
!> \param n_rank ...
!> \param counts ...
!> \param displs ...
! **************************************************************************************************
   SUBROUTINE grid_pair_blocks(n_grid, n_rank, counts, displs)
      INTEGER, INTENT(IN)                                :: n_grid, n_rank
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: counts(:), displs(:)

      INTEGER                                            :: i_rank

      ALLOCATE (counts(0:n_rank - 1), displs(0:n_rank - 1))
      DO i_rank = 0, n_rank - 1
         CALL block_range(n_grid**2, i_rank, n_rank, displs(i_rank), counts(i_rank))
      END DO
      displs(:) = displs - 1

   END SUBROUTINE grid_pair_blocks

! **************************************************************************************************
!> \brief Cache grid values of every Bloch molecular orbital without replicating all MO matrices.
!> \param bs_env ...
!> \param phi_lmu_k ...
!> \param psi_lm_k ...
! **************************************************************************************************
   SUBROUTINE build_all_psi_k(bs_env, phi_lmu_k, psi_lm_k)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: phi_lmu_k
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: psi_lm_k(:, :, :, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_all_psi_k'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: mo_coeff(:, :, :, :)
      INTEGER                                            :: handle, ikp_local, ispin

      CALL timeset(routineN, handle)

      CALL owned_mo_coefficients(bs_env, bs_env%cfm_mo_coeff_kp(1:bs_env%nkp_scf_desymm, :), mo_coeff)
      ALLOCATE (psi_lm_k(bs_env%ri_rs%n_grid_points, bs_env%n_ao, SIZE(phi_lmu_k, 3), &
                         bs_env%n_spin), SOURCE=z_zero)
      DO ispin = 1, bs_env%n_spin
         DO ikp_local = 1, SIZE(mo_coeff, 3)
            psi_lm_k(:, :, ikp_local, ispin) = MATMUL(phi_lmu_k(:, :, ikp_local), mo_coeff(:, :, ikp_local, ispin))
         END DO
      END DO
      DEALLOCATE (mo_coeff)

      CALL timestop(handle)

   END SUBROUTINE build_all_psi_k

! **************************************************************************************************
!> \brief Cache Bloch sums of the RI-RS grid values for the k-points owned by this rank.
!> \param bs_env ...
!> \param phi_lmu_k ...
! **************************************************************************************************
   SUBROUTINE build_all_phi_k(bs_env, phi_lmu_k)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: phi_lmu_k(:, :, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_all_phi_k'

      INTEGER                                            :: first, handle, i_R, ikp, image_width, &
                                                            last, n_k_local
      INTEGER, ALLOCATABLE                               :: owned_k(:)
      REAL(dp), ALLOCATABLE                              :: phi_R(:, :, :)

      CALL timeset(routineN, handle)

      n_k_local = (bs_env%nkp_scf_desymm + bs_env%para_env%num_pe - 1 - &
                   bs_env%para_env%mepos)/bs_env%para_env%num_pe
      ALLOCATE (phi_lmu_k(bs_env%ri_rs%n_grid_points, bs_env%n_ao, n_k_local))
      phi_lmu_k = z_zero
      IF (n_k_local > 0) THEN
         ALLOCATE (owned_k(n_k_local))
         owned_k(:) = [(ikp, ikp=bs_env%para_env%mepos + 1, bs_env%nkp_scf_desymm, bs_env%para_env%num_pe)]
         ! Retain a bounded panel of AO images; each image is sampled once.
         image_width = MAX(1, MIN(SIZE(bs_env%index_to_cell_3c, 2), &
                                  INT(projection_workspace/(real_bytes* &
                                                            REAL(bs_env%ri_rs%n_grid_points, dp)*REAL(bs_env%n_ao, dp)))))
         ALLOCATE (phi_R(bs_env%ri_rs%n_grid_points, bs_env%n_ao, image_width))
         DO first = 1, SIZE(bs_env%index_to_cell_3c, 2), image_width
            last = MIN(first + image_width - 1, SIZE(bs_env%index_to_cell_3c, 2))
            DO i_R = first, last
               CALL evaluate_mu_rl_R_values(bs_env, bs_env%ri_rs%grid_points, i_R, phi_R(:, :, i_R - first + 1))
            END DO
            CALL cell_to_k_transform(phi_R(:, :, 1:last - first + 1), &
                                     bs_env%index_to_cell_3c(:, first:last), &
                                     bs_env%kpoints_scf_desymm%xkp(:, owned_k), &
                                     bs_env%kpoints_scf_desymm%nkp_grid, phi_lmu_k, &
                                     bs_env%kpoints_scf_desymm%rskp_transform, accumulate=.TRUE., profile_fft=.TRUE.)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_all_phi_k

! **************************************************************************************************
!> \brief Fourier transform distributed real RI lattice blocks to cyclic q owners in one exchange.
!>
!> I_PQ(q) = Σ_R exp(2πiq·R) I_PQ(R).  Each sender transforms its own matrix block
!> for every q with the selected RSKP backend and sends the blocks to their q owners, iq mod P.
!> \param fm_R ...
!> \param cells ...
!> \param qpoints ...
!> \param nkp_grid regular q mesh of qpoints
!> \param matrices_q ...
!> \param backend GEMM or FFT
! **************************************************************************************************
   SUBROUTINE ri_2c_lattice_to_owned_q(fm_R, cells, qpoints, nkp_grid, matrices_q, backend)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_R(:)
      INTEGER, INTENT(IN)                                :: cells(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: qpoints(:, :)
      INTEGER, INTENT(IN)                                :: nkp_grid(3)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: matrices_q(:, :, :)
      INTEGER, INTENT(IN)                                :: backend

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ri_2c_lattice_to_owned_q'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: recv_buffer(:)
      COMPLEX(KIND=dp), ALLOCATABLE, TARGET              :: send_buffer(:)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: values_q(:, :, :)
      INTEGER                                            :: dest, handle, icell, icol, iq, iq_local, &
                                                            irow, n_q, n_q_local, n_RI, ncols, &
                                                            nprocs, nrecv, nrows, nsend, pos, &
                                                            source
      INTEGER, ALLOCATABLE :: block_metadata(:), q_order(:), recv_counts(:), recv_displs(:), &
         send_counts(:), send_displs(:), source_cols(:, :), source_rows(:, :), source_sizes(:, :)
      INTEGER, POINTER                                   :: cols(:), rows(:)
      REAL(KIND=dp), ALLOCATABLE                         :: values_R(:, :, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(fm_R) > 0)
      CPASSERT(SIZE(cells, 1) == 3 .AND. SIZE(cells, 2) == SIZE(fm_R))
      CPASSERT(SIZE(qpoints, 1) == 3)
      n_q = SIZE(qpoints, 2)
      n_RI = SIZE(matrices_q, 1)
      CPASSERT(SIZE(matrices_q, 2) == n_RI)
      para_env => fm_R(1)%matrix_struct%para_env
      nprocs = para_env%num_pe
      n_q_local = COUNT([(MODULO(iq, nprocs) == para_env%mepos, iq=1, n_q)])
      CPASSERT(SIZE(matrices_q, 3) == n_q_local)
      CALL cp_fm_get_info(fm_R(1), nrow_local=nrows, ncol_local=ncols, &
                          row_indices=rows, col_indices=cols)

      ALLOCATE (block_metadata(2 + 2*n_RI), source_sizes(2, 0:nprocs - 1), &
                source_rows(n_RI, 0:nprocs - 1), source_cols(n_RI, 0:nprocs - 1))
      DO source = 0, nprocs - 1
         block_metadata = 0
         IF (para_env%mepos == source) THEN
            block_metadata(1:2) = [nrows, ncols]
            block_metadata(3:2 + nrows) = rows(1:nrows)
            block_metadata(3 + n_RI:2 + n_RI + ncols) = cols(1:ncols)
         END IF
         CALL para_env%bcast(block_metadata, source)
         source_sizes(:, source) = block_metadata(1:2)
         source_rows(:, source) = block_metadata(3:2 + n_RI)
         source_cols(:, source) = block_metadata(3 + n_RI:2 + 2*n_RI)
      END DO

      ! The q points sorted by owner make the transformed blocks the send buffer itself.
      ALLOCATE (q_order(n_q), send_counts(0:nprocs - 1), recv_counts(0:nprocs - 1), &
                send_displs(0:nprocs - 1), recv_displs(0:nprocs - 1))
      pos = 0
      DO dest = 0, nprocs - 1
         send_displs(dest) = pos*nrows*ncols
         DO iq = 1, n_q
            IF (MODULO(iq, nprocs) /= dest) CYCLE
            pos = pos + 1
            q_order(pos) = iq
         END DO
         send_counts(dest) = pos*nrows*ncols - send_displs(dest)
      END DO
      nsend = n_q*nrows*ncols
      CALL para_env%alltoall(send_counts, recv_counts, 1)
      recv_displs(0) = 0
      DO source = 1, nprocs - 1
         recv_displs(source) = recv_displs(source - 1) + recv_counts(source - 1)
      END DO
      nrecv = recv_displs(nprocs - 1) + recv_counts(nprocs - 1)
      ALLOCATE (send_buffer(MAX(1, nsend)), recv_buffer(MAX(1, nrecv)))
      IF (nsend > 0) THEN
         ALLOCATE (values_R(nrows, ncols, SIZE(fm_R)))
         DO icell = 1, SIZE(fm_R)
            values_R(:, :, icell) = fm_R(icell)%local_data(1:nrows, 1:ncols)
         END DO
         values_q(1:nrows, 1:ncols, 1:n_q) => send_buffer(1:nsend)
         CALL cell_to_k_transform(values_R, cells, qpoints(:, q_order), nkp_grid, values_q, backend, profile_fft=.TRUE.)
         DEALLOCATE (values_R)
      END IF
      CALL para_env%alltoall(send_buffer, send_counts, send_displs, &
                             recv_buffer, recv_counts, recv_displs)
      matrices_q = z_zero
      DO source = 0, nprocs - 1
         nrows = source_sizes(1, source)
         ncols = source_sizes(2, source)
         CPASSERT(recv_counts(source) == nrows*ncols*n_q_local)
         pos = recv_displs(source)
         DO iq_local = 1, n_q_local
            DO icol = 1, ncols
               DO irow = 1, nrows
                  matrices_q(source_rows(irow, source), source_cols(icol, source), iq_local) = &
                     recv_buffer(pos + (icol - 1)*nrows + irow)
               END DO
            END DO
            pos = pos + nrows*ncols
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE ri_2c_lattice_to_owned_q

! **************************************************************************************************
!> \brief Evaluate the band Bloch orbitals of the k points owned by this rank, MOD(k-1,P), one AO
!> cell at a time. No rank stores the orbitals of other k points.
!> \param bs_env ...
!> \param orbitals orbitals(:, :, (k-1)/P+1, spin) on the RI-RS grid
! **************************************************************************************************
   SUBROUTINE collect_band_bloch_orbitals(bs_env, orbitals)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: orbitals(:, :, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'collect_band_bloch_orbitals'

      COMPLEX(dp), ALLOCATABLE                           :: coefficients(:, :, :, :), phi_k(:, :)
      INTEGER                                            :: handle, i_kp, i_local, i_R, ispin, &
                                                            n_grid, n_own
      REAL(dp)                                           :: available_GB, orbital_bytes
      REAL(dp), ALLOCATABLE                              :: phi_R(:, :)
      TYPE(mp_para_env_type), POINTER                    :: env

      CALL timeset(routineN, handle)

      env => bs_env%para_env
      n_grid = bs_env%ri_rs%n_grid_points
      n_own = (bs_env%nkp_bs_and_DOS - env%mepos + env%num_pe - 1)/env%num_pe
      CALL collect_band_mo_coefficients(bs_env, coefficients)
      orbital_bytes = REAL(n_grid, dp)*REAL(bs_env%n_ao, dp)* &
                      (complex_bytes*REAL(n_own*bs_env%n_spin, dp) + complex_bytes + real_bytes)
      CALL mp_mem_avail_per_rank_GB(env, available_GB)
      CALL env%min(available_GB)
      CALL env%max(orbital_bytes)
      IF (orbital_bytes > stage_memory_share*MAX(0.0_dp, available_GB)*bytes_per_GB) THEN
         CPABORT("Insufficient memory for GW band-path orbitals")
      END IF
      ALLOCATE (phi_R(n_grid, bs_env%n_ao))
      ALLOCATE (orbitals(n_grid, bs_env%n_ao, n_own, bs_env%n_spin))
      orbitals = z_zero
      DO i_R = 1, SIZE(bs_env%index_to_cell_3c, 2)
         IF (n_own == 0) EXIT
         CALL evaluate_mu_rl_R_values(bs_env, bs_env%ri_rs%grid_points, i_R, phi_R)
         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
         !$OMP SHARED(bs_env, env, i_R, n_own, orbitals, phi_R) PRIVATE(i_kp, i_local)
         DO i_local = 1, n_own
            i_kp = env%mepos + 1 + (i_local - 1)*env%num_pe
            orbitals(:, :, i_local, 1) = orbitals(:, :, i_local, 1) + phi_R* &
                                         EXP(gaussi*twopi*SUM(bs_env%kpoints_DOS%xkp(:, i_kp)* &
                                                              REAL(bs_env%index_to_cell_3c(:, i_R), dp)))
         END DO
         !$OMP END PARALLEL DO
      END DO
      !$OMP PARALLEL DEFAULT(NONE) SHARED(bs_env, coefficients, n_grid, n_own, orbitals) &
      !$OMP PRIVATE(i_local, ispin, phi_k)
      ALLOCATE (phi_k(n_grid, bs_env%n_ao))
      !$OMP DO SCHEDULE(STATIC)
      DO i_local = 1, n_own
         phi_k(:, :) = orbitals(:, :, i_local, 1)
         DO ispin = 1, bs_env%n_spin
            CALL zgemm('N', 'N', n_grid, bs_env%n_ao, bs_env%n_ao, z_one, &
                       phi_k, n_grid, coefficients(:, :, i_local, ispin), &
                       bs_env%n_ao, z_zero, orbitals(:, :, i_local, ispin), n_grid)
         END DO
      END DO
      !$OMP END DO
      DEALLOCATE (phi_k)
      !$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE collect_band_bloch_orbitals

! **************************************************************************************************
!> \brief Move the band orbitals from their owners in all ranks, MOD(k-1,P), to their owners in
!> the time groups of Σc, MOD(k-1,g). As g divides P, both owners have the same rank within their
!> group, so every k point is broadcast once among the ranks with this group rank.
!> \param bs_env ...
!> \param gwkp_env ...
! **************************************************************************************************
   SUBROUTINE band_orbitals_to_groups(bs_env, gwkp_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT)                     :: gwkp_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'band_orbitals_to_groups'

      COMPLEX(dp), ALLOCATABLE                           :: group_orbitals(:, :, :, :)
      INTEGER                                            :: handle, i_kp, i_local, i_root, ispin, &
                                                            n_group, n_own, n_rank
      TYPE(mp_comm_type)                                 :: column_env

      CALL timeset(routineN, handle)

      n_rank = bs_env%para_env%num_pe
      n_group = gwkp_env%group_env%num_pe
      n_own = (bs_env%nkp_bs_and_DOS - gwkp_env%group_env%mepos + n_group - 1)/n_group
      ALLOCATE (group_orbitals(SIZE(gwkp_env%band_orbitals, 1), SIZE(gwkp_env%band_orbitals, 2), &
                               n_own, bs_env%n_spin))
      CALL column_env%from_split(bs_env%para_env, gwkp_env%group_env%mepos)
      DO i_local = 1, n_own
         i_kp = gwkp_env%group_env%mepos + 1 + (i_local - 1)*n_group
         i_root = MODULO(i_kp - 1, n_rank)/n_group
         DO ispin = 1, bs_env%n_spin
            IF (column_env%mepos == i_root) THEN
               group_orbitals(:, :, i_local, ispin) = gwkp_env%band_orbitals(:, :, (i_kp - 1)/n_rank + 1, ispin)
            END IF
            CALL column_env%bcast(group_orbitals(:, :, i_local, ispin), i_root)
         END DO
      END DO
      CALL column_env%free()
      CALL MOVE_ALLOC(group_orbitals, gwkp_env%band_orbitals)

      CALL timestop(handle)

   END SUBROUTINE band_orbitals_to_groups

! **************************************************************************************************
!> \brief Return the SCF-mesh index of a band k point, or 0 if it lies off the mesh.
!> The first nkp_scf_desymm band points are the SCF mesh in the same order.
!> \param bs_env ...
!> \param i_kp_eval Band k point.
!> \return ...
! **************************************************************************************************
   FUNCTION native_q_index(bs_env, i_kp_eval) RESULT(native_iq)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: i_kp_eval
      INTEGER                                            :: native_iq

      INTEGER                                            :: iq
      REAL(KIND=dp)                                      :: delta_k(3)

      native_iq = 0
      DO iq = 1, bs_env%nkp_scf_desymm
         IF (i_kp_eval <= bs_env%nkp_scf_desymm .AND. iq /= i_kp_eval) CYCLE
         delta_k = bs_env%kpoints_DOS%xkp(:, i_kp_eval) - bs_env%kpoints_scf_desymm%xkp(:, iq)
         IF (MAXVAL(ABS(delta_k - ANINT(delta_k))) <= 64.0_dp*EPSILON(1.0_dp)) THEN
            native_iq = iq
            RETURN
         END IF
      END DO

   END FUNCTION native_q_index

! **************************************************************************************************
!> \brief MO coefficients of the band k points owned by this rank, MOD(k-1,P).
!> \param bs_env ...
!> \param band_coefficients ...
! **************************************************************************************************
   SUBROUTINE collect_band_mo_coefficients(bs_env, band_coefficients)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: band_coefficients(:, :, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'collect_band_mo_coefficients'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: coefficients(:, :)
      INTEGER                                            :: handle, i_kp_eval, ispin
      TYPE(mp_para_env_type), POINTER                    :: env

      CALL timeset(routineN, handle)

      env => bs_env%para_env
      ALLOCATE (coefficients(bs_env%n_ao, bs_env%n_ao))
      ALLOCATE (band_coefficients(bs_env%n_ao, bs_env%n_ao, &
                                  (bs_env%nkp_bs_and_DOS - env%mepos + env%num_pe - 1)/env%num_pe, &
                                  bs_env%n_spin))
      DO ispin = 1, bs_env%n_spin
         DO i_kp_eval = 1, bs_env%nkp_bs_and_DOS
            CALL cp_cfm_get_submatrix(bs_env%cfm_mo_coeff_kp(i_kp_eval, ispin), coefficients)
            IF (MODULO(i_kp_eval - 1, env%num_pe) == env%mepos) THEN
               band_coefficients(:, :, (i_kp_eval - 1)/env%num_pe + 1, ispin) = coefficients
            END IF
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE collect_band_mo_coefficients

! **************************************************************************************************
!> \brief Σ_ll′(k) = Σ_R exp(2πik·R) Σ_ll′(R) on the native k mesh for the ordered grid pairs of
!> this rank. FFT packs the original input without changing Σ(R).
!> \param bs_env ...
!> \param cells ...
!> \param sigma_R ...
!> \param first_pair First locally owned ordered grid pair.
!> \param sigma_mesh ...
!> \param fft_work ...
! **************************************************************************************************
   SUBROUTINE self_energy_on_mesh(bs_env, cells, sigma_R, first_pair, sigma_mesh, fft_work)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: cells(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: sigma_R(:, :, :)
      INTEGER, INTENT(IN)                                :: first_pair
      COMPLEX(KIND=dp), INTENT(OUT)                      :: sigma_mesh(:, :, :)
      TYPE(gw_pair_fft_work_type), INTENT(INOUT), &
         OPTIONAL, TARGET                                :: fft_work

      CALL product_to_k_pair_centered(bs_env, cells, bs_env%kpoints_scf_desymm%xkp, sigma_R, sigma_mesh, &
                                      first_pair, fft_work=fft_work)
   END SUBROUTINE self_energy_on_mesh

! **************************************************************************************************
!> \brief Project Σ_ll′(k) onto the band Bloch orbitals, Σ_n(k) = Re[ψ_n(k)† Σ(k) ψ_n(k)].
!>
!> The pair owners evaluate Σ_ll′(k) = Σ_R exp(2πik·R) Σ_ll′(R) at the band k points. Every k
!> owner, MOD(k-1,P), receives its complete matrices and contracts them with dense kernels. A
!> round covers m k points per rank, so that one pass over Σ(R) serves up to
!> band_kpoints_per_pass k points.
!> \param bs_env ...
!> \param sigma_mesh Σ on the native k mesh for the ordered grid pairs of this rank.
!> \param sigma_R Σ in the cells of the native mesh for the same pairs.
!> \param first_pair First locally owned ordered grid pair.
!> \param env Ranks that share the grid pairs and the band k points.
!> \param orbitals Band Bloch orbitals of one spin at the k points owned in env.
!> \param sigma_band ...
! **************************************************************************************************
   SUBROUTINE self_energy_pairs_to_Bloch(bs_env, sigma_mesh, sigma_R, first_pair, env, orbitals, sigma_band)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), INTENT(IN)                       :: sigma_mesh(:, :, :)
      REAL(KIND=dp), INTENT(IN)                          :: sigma_R(:, :, :)
      INTEGER, INTENT(IN)                                :: first_pair
      TYPE(mp_para_env_type), POINTER                    :: env
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(IN)           :: orbitals(:, :, :)
      REAL(KIND=dp), INTENT(OUT)                         :: sigma_band(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'self_energy_pairs_to_Bloch'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: sigma_k(:, :), work(:, :)
      COMPLEX(KIND=dp), ALLOCATABLE, TARGET              :: send_buffer(:)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: sigma_columns(:, :)
      INTEGER                                            :: handle, i_band, i_batch, i_first, i_kp, &
                                                            i_rank, j, m, n_batch, n_grid, n_kp, &
                                                            n_local, n_rank, n_round
      INTEGER, ALLOCATABLE                               :: pair_counts(:), pair_displs(:), &
                                                            recv_counts(:), send_counts(:), &
                                                            send_displs(:)

      CALL timeset(routineN, handle)

      n_grid = bs_env%ri_rs%n_grid_points
      n_kp = bs_env%nkp_bs_and_DOS
      n_local = SIZE(sigma_mesh, 1)
      n_rank = env%num_pe
      n_batch = MIN(16, bs_env%n_ao)
      m = projection_round_kpoints(n_grid, n_rank)
      CALL grid_pair_blocks(n_grid, n_rank, pair_counts, pair_displs)
      ALLOCATE (send_buffer(n_local*n_rank*m), sigma_k(n_grid**2, m), recv_counts(0:n_rank - 1), &
                send_counts(0:n_rank - 1), send_displs(0:n_rank - 1))
      sigma_columns(1:n_local, 1:n_rank*m) => send_buffer
      send_displs(:) = [(i_rank*n_local, i_rank=0, n_rank - 1)]
      sigma_band = 0.0_dp
      ! Column i_col of a round holds the band k point i_first+i_col-1, owned by MOD(i_col-1,P).
      DO i_first = 1, n_kp, n_rank*m
         n_round = MIN(n_rank*m, n_kp - i_first + 1)
         CALL band_self_energy(bs_env, sigma_mesh, sigma_R, first_pair, i_first, sigma_columns(:, 1:n_round))
         DO j = 1, m
            DO i_rank = 0, n_rank - 1
               send_counts(i_rank) = MERGE(n_local, 0, (j - 1)*n_rank + i_rank < n_round)
            END DO
            i_kp = i_first + (j - 1)*n_rank + env%mepos
            recv_counts(:) = MERGE(pair_counts, 0, i_kp <= n_kp)
            CALL env%alltoall(send_buffer((j - 1)*n_rank*n_local + 1:j*n_rank*n_local), send_counts, &
                              send_displs, sigma_k(:, j), recv_counts, pair_displs)
         END DO
         !$OMP PARALLEL DEFAULT(NONE) &
         !$OMP SHARED(bs_env, env, i_first, m, n_batch, n_grid, n_kp, n_rank, orbitals, sigma_band, sigma_k) &
         !$OMP PRIVATE(i_band, i_batch, i_kp, j, work)
         ALLOCATE (work(n_grid, n_batch))
         !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
         DO j = 1, m
            DO i_batch = 1, bs_env%n_ao, n_batch
               i_kp = i_first + (j - 1)*n_rank + env%mepos
               IF (i_kp > n_kp) CYCLE
               CALL zgemm('N', 'N', n_grid, MIN(n_batch, bs_env%n_ao - i_batch + 1), n_grid, z_one, &
                          sigma_k(1, j), n_grid, &
                          orbitals(:, i_batch:MIN(i_batch + n_batch - 1, bs_env%n_ao), (i_kp - 1)/n_rank + 1), n_grid, &
                          z_zero, work, n_grid)
               DO i_band = i_batch, MIN(i_batch + n_batch - 1, bs_env%n_ao)
                  sigma_band(i_band, i_kp) = REAL(SUM(CONJG(orbitals(:, i_band, (i_kp - 1)/n_rank + 1))* &
                                                      work(:, i_band - i_batch + 1)), dp)
               END DO
            END DO
         END DO
         !$OMP END DO
         DEALLOCATE (work)
         !$OMP END PARALLEL
      END DO
      CALL env%sum(sigma_band)

      CALL timestop(handle)

   END SUBROUTINE self_energy_pairs_to_Bloch

! **************************************************************************************************
!> \brief Weight of the image offset c along one direction of an n_mesh-point k mesh: 1 inside
!> the Born-von Karman supercell, the given boundary weight on its boundary 2|c| = n_mesh, 0
!> outside. A pair on a half-cell tie averages the image sets centred at s and s+1.
!> \param c ...
!> \param n_mesh ...
!> \param tie ...
!> \param boundary Weight of each of the two boundary images.
!> \return ...
! **************************************************************************************************
   PURE REAL(KIND=dp) FUNCTION image_weight(c, n_mesh, tie, boundary) RESULT(weight)
      INTEGER, INTENT(IN)                                :: c, n_mesh
      LOGICAL, INTENT(IN)                                :: tie
      REAL(KIND=dp), INTENT(IN)                          :: boundary

      INTEGER                                            :: i_center, n_centers

      n_centers = MERGE(2, 1, tie)
      weight = 0.0_dp
      DO i_center = 0, n_centers - 1
         IF (2*ABS(c - i_center) < n_mesh) THEN
            weight = weight + 1.0_dp/REAL(n_centers, dp)
         ELSE IF (2*ABS(c - i_center) == n_mesh) THEN
            weight = weight + boundary/REAL(n_centers, dp)
         END IF
      END DO
   END FUNCTION image_weight

! **************************************************************************************************
!> \brief Bound the band-projection and Green-function workspaces of one rank in bytes.
!> The projection holds the Σ(k) columns of one round, the complete Σ(k) of the owned k points,
!> the shared lattice phases, thread-private image, sum and output blocks, and the 64 MiB copy of
!> the mesh transform. The Green function holds one complete G(k), its received pair blocks and
!> the scaled orbitals.
!> \param n_grid ...
!> \param n_rank ...
!> \param n_images Bound of the image cells summed per pair.
!> \param n_ao ...
!> \param n_threads ...
!> \return ...
! **************************************************************************************************
   FUNCTION projection_workspace_bytes(n_grid, n_rank, n_images, n_ao, n_threads) RESULT(bytes)
      INTEGER, INTENT(IN)                                :: n_grid, n_rank, n_images, n_ao, n_threads
      REAL(dp)                                           :: bytes

      REAL(dp)                                           :: m, n_round

      m = REAL(projection_round_kpoints(n_grid, n_rank), dp)
      n_round = m*REAL(n_rank, dp)
      bytes = 2.0_dp*complex_bytes*m*REAL(n_grid, dp)**2 + complex_bytes*n_round*REAL(n_images, dp) + &
              projection_workspace + &
              REAL(n_threads, dp)*(complex_bytes*REAL(n_grid, dp)*REAL(MIN(16, n_ao), dp) + &
                                   real_bytes*REAL(pair_batch_size, dp)*REAL(n_images, dp) + &
                                   complex_bytes*REAL(pair_batch_size, dp)*n_round)
      bytes = MAX(bytes, complex_bytes*REAL(n_grid, dp)*(2.0_dp*REAL(n_grid, dp) + REAL(n_ao, dp)))
   END FUNCTION projection_workspace_bytes

! **************************************************************************************************
!> \brief Band k points per rank and projection round: a round evaluates up to
!> band_kpoints_per_pass k points, while the complete Σ(k) of a round stay below 512 MiB.
!> \param n_grid ...
!> \param n_rank ...
!> \return ...
! **************************************************************************************************
   PURE INTEGER FUNCTION projection_round_kpoints(n_grid, n_rank) RESULT(m)
      INTEGER, INTENT(IN)                                :: n_grid, n_rank

      m = MAX(1, MIN(band_kpoints_per_pass/n_rank, INT(projection_round_bytes/(2.0_dp*complex_bytes*REAL(n_grid, dp)**2))))
   END FUNCTION projection_round_kpoints

! **************************************************************************************************
!> \brief Bytes left for a stage after the reserve for data that stays alive and for libraries.
!> \param available_GB ...
!> \return ...
! **************************************************************************************************
   PURE REAL(KIND=dp) FUNCTION budget_after_reserve(available_GB) RESULT(budget)
      REAL(KIND=dp), INTENT(IN)                          :: available_GB

      budget = MAX(0.0_dp, available_GB - MAX(min_reserve_GB, reserve_share*available_GB))*bytes_per_GB
   END FUNCTION budget_after_reserve

! **************************************************************************************************
!> \brief Σ_ll′(k) = Σ_R exp(2πik·R) Σ_ll′(R) at the band k points of one projection round.
!> Native k points are copied from the mesh transform. Elsewhere every pair sums its images
!> R = s+c around the lattice shift s nearest to H⁻¹(r_l-r_l′), weighted by image_weight in
!> each direction; on the native mesh this is the pair-centered lattice FFT. The two images on
!> the supercell boundary share (1+cos(2πk₀N))/2 for the mesh offset k₀: on a Γ-centered mesh
!> they are equivalent, on a half-shifted mesh G and W are antiperiodic there and the images
!> cancel. All pairs share the phases exp(2πik·c), so batches of pairs are summed by two real
!> GEMMs.
!> \param bs_env ...
!> \param sigma_mesh ...
!> \param sigma_R ...
!> \param first_pair ...
!> \param first_kp_eval Band k point of the first column of sigma_k.
!> \param sigma_k ...
! **************************************************************************************************
   SUBROUTINE band_self_energy(bs_env, sigma_mesh, sigma_R, first_pair, first_kp_eval, sigma_k)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), INTENT(IN)                       :: sigma_mesh(:, :, :)
      REAL(KIND=dp), INTENT(IN)                          :: sigma_R(:, :, :)
      INTEGER, INTENT(IN)                                :: first_pair, first_kp_eval
      COMPLEX(KIND=dp), INTENT(OUT)                      :: sigma_k(:, :)

      COMPLEX(KIND=dp)                                   :: phase
      INTEGER :: c1, c2, c3, first_c(3), i_c, i_col, i_first, i_local, i_off, i_pair, last_c(3), &
         lower(3), n_c, n_grid, n_off, n_pairs, native_iq, nmesh(3), shift(3), wrapped(3)
      INTEGER, ALLOCATABLE                               :: offmesh_col(:), pair_shift(:, :)
      LOGICAL                                            :: any_tie(3), tie(3)
      REAL(KIND=dp)                                      :: angle, boundary(3)
      REAL(KIND=dp), ALLOCATABLE :: cos_c(:, :), fractional_grid(:, :), images(:, :), sin_c(:, :), &
         sum_im(:, :), sum_re(:, :), weight_1(:), weight_2(:), weight_3(:), xk(:, :)

      ALLOCATE (offmesh_col(SIZE(sigma_k, 2)))
      n_off = 0
      DO i_col = 1, SIZE(sigma_k, 2)
         native_iq = native_q_index(bs_env, first_kp_eval + i_col - 1)
         IF (native_iq > 0) THEN
            sigma_k(:, i_col) = sigma_mesh(:, 1, native_iq)
         ELSE
            n_off = n_off + 1
            offmesh_col(n_off) = i_col
         END IF
      END DO
      IF (n_off == 0) RETURN

      n_grid = bs_env%ri_rs%n_grid_points
      nmesh = bs_env%kpoints_scf_desymm%nkp_grid
      lower = -nmesh/2
      ALLOCATE (fractional_grid(3, n_grid))
      fractional_grid(:, :) = MATMUL(bs_env%ri_rs%cell%h_inv, bs_env%ri_rs%grid_points)
      DO i_c = 1, 3
         IF (bs_env%periodic(i_c) == 0) fractional_grid(i_c, :) = 0.0_dp
      END DO
      any_tie = .FALSE.
      DO i_local = 1, SIZE(sigma_R, 1)
         CALL product_pair_shift(first_pair + i_local - 1, n_grid, fractional_grid, shift, tie)
         any_tie = any_tie .OR. tie
      END DO
      boundary = 0.25_dp*(1.0_dp + COS(twopi*bs_env%kpoints_scf_desymm%xkp(:, 1)*nmesh))
      first_c = lower
      last_c = nmesh/2 + MERGE(1, 0, any_tie)
      n_c = PRODUCT(last_c - first_c + 1)
      ALLOCATE (xk(3, n_off), cos_c(n_c, n_off), sin_c(n_c, n_off))
      DO i_off = 1, n_off
         xk(:, i_off) = bs_env%kpoints_DOS%xkp(:, first_kp_eval + offmesh_col(i_off) - 1)
      END DO
      ! Image c of every pair is column i_c of the shared phases, c1 running fastest.
      !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
      !$OMP SHARED(cos_c, first_c, last_c, n_off, sin_c, xk) PRIVATE(c1, c2, c3, i_c, i_off, phase)
      DO i_off = 1, n_off
         i_c = 0
         DO c3 = first_c(3), last_c(3)
            DO c2 = first_c(2), last_c(2)
               DO c1 = first_c(1), last_c(1)
                  i_c = i_c + 1
                  phase = EXP(CMPLX(0.0_dp, twopi*(xk(1, i_off)*c1 + xk(2, i_off)*c2 + xk(3, i_off)*c3), KIND=dp))
                  cos_c(i_c, i_off) = REAL(phase, KIND=dp)
                  sin_c(i_c, i_off) = AIMAG(phase)
               END DO
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(boundary, cos_c, first_c, first_pair, fractional_grid, last_c, lower, n_c, n_grid, n_off, &
      !$OMP        nmesh, offmesh_col, sigma_k, sigma_R, sin_c, xk) &
      !$OMP PRIVATE(angle, c1, c2, c3, i_c, i_first, i_off, i_pair, images, n_pairs, pair_shift, shift, &
      !$OMP         sum_im, sum_re, tie, weight_1, weight_2, weight_3, wrapped)
      ALLOCATE (images(n_c, pair_batch_size), pair_shift(3, pair_batch_size), &
                sum_re(pair_batch_size, n_off), sum_im(pair_batch_size, n_off), &
                weight_1(first_c(1):last_c(1)), weight_2(first_c(2):last_c(2)), weight_3(first_c(3):last_c(3)))
      !$OMP DO SCHEDULE(DYNAMIC)
      DO i_first = 1, SIZE(sigma_R, 1), pair_batch_size
         n_pairs = MIN(pair_batch_size, SIZE(sigma_R, 1) - i_first + 1)
         DO i_pair = 1, n_pairs
            CALL product_pair_shift(first_pair + i_first + i_pair - 2, n_grid, fractional_grid, shift, tie)
            pair_shift(:, i_pair) = shift
            weight_1(:) = [(image_weight(c1, nmesh(1), tie(1), boundary(1)), c1=first_c(1), last_c(1))]
            weight_2(:) = [(image_weight(c2, nmesh(2), tie(2), boundary(2)), c2=first_c(2), last_c(2))]
            weight_3(:) = [(image_weight(c3, nmesh(3), tie(3), boundary(3)), c3=first_c(3), last_c(3))]
            i_c = 0
            DO c3 = first_c(3), last_c(3)
               DO c2 = first_c(2), last_c(2)
                  DO c1 = first_c(1), last_c(1)
                     i_c = i_c + 1
                     IF (weight_1(c1)*weight_2(c2)*weight_3(c3) == 0.0_dp) THEN
                        images(i_c, i_pair) = 0.0_dp
                     ELSE
                        wrapped = MODULO(shift + [c1, c2, c3] - lower, nmesh) + lower
                        images(i_c, i_pair) = weight_1(c1)*weight_2(c2)*weight_3(c3)* &
                                              sigma_R(i_first + i_pair - 1, 1, 1 + wrapped(3) - lower(3) + &
                                                      nmesh(3)*(wrapped(2) - lower(2) + nmesh(2)*(wrapped(1) - lower(1))))
                     END IF
                  END DO
               END DO
            END DO
         END DO
         CALL dgemm('T', 'N', n_pairs, n_off, n_c, 1.0_dp, images, n_c, cos_c, n_c, 0.0_dp, sum_re, pair_batch_size)
         CALL dgemm('T', 'N', n_pairs, n_off, n_c, 1.0_dp, images, n_c, sin_c, n_c, 0.0_dp, sum_im, pair_batch_size)
         DO i_off = 1, n_off
            DO i_pair = 1, n_pairs
               angle = twopi*DOT_PRODUCT(xk(:, i_off), REAL(pair_shift(:, i_pair), dp))
               sigma_k(i_first + i_pair - 1, offmesh_col(i_off)) = &
                  EXP(CMPLX(0.0_dp, angle, KIND=dp))*CMPLX(sum_re(i_pair, i_off), sum_im(i_pair, i_off), KIND=dp)
            END DO
         END DO
      END DO
      !$OMP END DO
      DEALLOCATE (images, pair_shift, sum_re, sum_im, weight_1, weight_2, weight_3)
      !$OMP END PARALLEL
      DEALLOCATE (cos_c, fractional_grid, offmesh_col, sin_c, xk)
   END SUBROUTINE band_self_energy

! **************************************************************************************************
!> \brief Transform a periodic lattice product in a cell centred on each physical grid pair.
!> Products of two half-shifted lattice transforms are periodic. For s nearest to
!> H^{-1}(r_l-r_l'), F_ll'(k)=exp(2πiks) Σ_R exp(2πikR) F_ll'(R+s).
!> The even-mesh Nyquist endpoints and half-cell choices of s are averaged symmetrically.
!> On integer/half-shifted band meshes these averages are real separable endpoint weights,
!> allowing one FFT. Centering is folded into input packing, leaving values_R unchanged.
!> FFT tiles fold the source permutation and output Bloch phase into separable input factors.
!> The rank owns the flattened column-major pairs first_pair, first_pair+1, ...
!> \param bs_env ...
!> \param cells ...
!> \param xkp ...
!> \param values_R ...
!> \param values_k ...
!> \param first_pair ...
!> \param fft_work ...
! **************************************************************************************************
   SUBROUTINE product_to_k_pair_centered(bs_env, cells, xkp, values_R, values_k, first_pair, fft_work)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: cells(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: values_R(:, :, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: values_k(:, :, :)
      INTEGER, INTENT(IN)                                :: first_pair
      TYPE(gw_pair_fft_work_type), INTENT(INOUT), &
         OPTIONAL, TARGET                                :: fft_work

      CHARACTER(LEN=*), PARAMETER :: routineN = 'product_to_k_pair_centered'
      INTEGER(KIND=int_8), PARAMETER :: fft_panel_bytes = 32_int_8*1024_int_8**2
      INTEGER, PARAMETER                                 :: pack_cells = 32, pack_pairs = 16

      COMPLEX(KIND=dp) :: pack_buffer(pack_pairs, pack_cells)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: origin_phase(:)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: fft_in(:), fft_out(:), input(:, :), &
                                                            output(:, :, :, :)
      INTEGER :: cell_first, cell_last, class_x, class_y, class_z, combined_count, combined_id, &
         coordinate, count, first, handle, i, ibatch, icell, idir, ik, ilocal, isource, lower(3), &
         naxis, nfft(3), ngrid, nmesh(3), row_first, row_last, stat, unit_nr, width
      INTEGER(KIND=int_8)                                :: axis_classes, cached_coefficients, &
                                                            combined_classes, nfft_grid, &
                                                            phase_evaluations
      INTEGER, ALLOCATABLE                               :: cell_coordinate(:, :), cell_index(:), &
                                                            grid_index(:, :), keys(:, :), &
                                                            order(:), pair_shift(:, :)
      LOGICAL                                            :: cache_hit, unpadded
      LOGICAL, ALLOCATABLE                               :: pair_tie(:, :)
      REAL(KIND=dp)                                      :: endpoint(3), k_offset(3), setup_seconds, &
                                                            started, total_started, weight
      REAL(KIND=dp), ALLOCATABLE                         :: fractional_grid(:, :)
      TYPE(gw_pair_fft_batch_type), POINTER              :: batch
      TYPE(gw_pair_fft_work_type), POINTER               :: plan
      TYPE(gw_pair_fft_work_type), TARGET                :: local_fft_work
      TYPE(lattice_fft_profile_type)                     :: profile

      CALL timeset(routineN, handle)
      total_started = m_walltime()
      profile = lattice_fft_profile_type()

      ngrid = bs_env%ri_rs%n_grid_points
      nmesh = bs_env%kpoints_scf_desymm%nkp_grid
      lower = -nmesh/2
      endpoint = COS(twopi*xkp(1:3, 1)*REAL(nmesh, dp))
      ALLOCATE (fractional_grid(3, ngrid))
      fractional_grid(:, :) = MATMUL(bs_env%ri_rs%cell%h_inv, bs_env%ri_rs%grid_points)
      DO idir = 1, 3
         IF (bs_env%periodic(idir) == 0) fractional_grid(idir, :) = 0.0_dp
      END DO
      CPASSERT(SIZE(values_R, 2) == 1)
      ! GW's physical centering stays here; the common FFT executor sees only packed grids.
      ALLOCATE (grid_index(3, SIZE(xkp, 2)), cell_coordinate(3, SIZE(cells, 2)), &
                origin_phase(SIZE(xkp, 2)))
      IF (.NOT. lattice_fft_prepare_geometry(xkp, nmesh, cells, grid_index, nfft, k_offset, &
                                             cell_linear=cell_index, cell_origin=lower, &
                                             origin_phase=origin_phase)) THEN
         CPABORT('GW product needs a supported full FFT mesh')
      END IF
      CPASSERT(SIZE(cells, 2, KIND=int_8) == PRODUCT(INT(nmesh, int_8)))
      CPASSERT(SIZE(values_R, 3) == SIZE(cells, 2))
      IF (ANY(SHAPE(values_k) /= [SIZE(values_R, 1), 1, SIZE(xkp, 2)])) THEN
         CPABORT('Inconsistent GW FFT output')
      END IF
      DO icell = 1, SIZE(cells, 2)
         isource = 1 + cells(3, icell) - lower(3) + &
                   nmesh(3)*(cells(2, icell) - lower(2) + nmesh(2)*(cells(1, icell) - lower(1)))
         IF (isource /= icell) CPABORT('GW product needs canonical cells')
         ! The shared geometry maps canonical cells directly into z-fastest FFT storage.
         cell_coordinate(:, icell) = cells(:, icell) - lower + 1
      END DO
      unpadded = ALL(nfft == nmesh([3, 2, 1]))
      profile%seconds(2) = m_walltime() - total_started
      nfft_grid = PRODUCT(INT(nfft, int_8))
      width = INT(MIN(MAX(1_int_8, SIZE(values_R, 1, KIND=int_8)), &
                      MAX(1_int_8, fft_panel_bytes/(48_int_8*nfft_grid))))
      CPASSERT(nfft_grid*width <= HUGE(0))
      plan => local_fft_work
      IF (PRESENT(fft_work)) plan => fft_work
      cache_hit = plan%initialized
      IF (cache_hit) THEN
         CPASSERT(plan%first_pair == first_pair)
         CPASSERT(plan%n_pairs == SIZE(values_R, 1))
         CPASSERT(plan%width == width)
         CPASSERT(plan%ngrid == ngrid)
         CPASSERT(ALL(plan%nmesh == nmesh))
         CPASSERT(ALL(SHAPE(plan%xkp) == [3, SIZE(xkp, 2)]))
         CPASSERT(ALL(SHAPE(plan%cells) == [3, SIZE(cells, 2)]))
         CPASSERT(ALL(plan%xkp == xkp(1:3, :)))
         CPASSERT(ALL(plan%cells == cells(1:3, :)))
      ELSE
         plan%first_pair = first_pair
         plan%n_pairs = SIZE(values_R, 1)
         plan%width = width
         plan%ngrid = ngrid
         plan%nmesh = nmesh
         ALLOCATE (plan%xkp, SOURCE=xkp(1:3, :))
         ALLOCATE (plan%cells, SOURCE=cells(1:3, :))
         ALLOCATE (plan%batches((SIZE(values_R, 1) + width - 1)/width))
      END IF
      started = m_walltime()
      CALL fft_alloc(fft_in, [INT(nfft_grid)*width])
      CALL fft_alloc(fft_out, [INT(nfft_grid)*width])
      IF (.NOT. cache_hit) THEN
         ALLOCATE (pair_shift(3, width), pair_tie(3, width), keys(3, width), order(width))
      END IF
      profile%seconds(3) = m_walltime() - started
      setup_seconds = 0.0_dp
      phase_evaluations = 0_int_8
      axis_classes = 0_int_8
      combined_classes = 0_int_8
      cached_coefficients = 0_int_8
      DO first = 1, SIZE(values_R, 1), width
         count = MIN(width, SIZE(values_R, 1) - first + 1)
         ibatch = (first - 1)/width + 1
         batch => plan%batches(ibatch)
         profile%elements = profile%elements + INT(count, int_8)
         profile%batches = profile%batches + 1_int_8
         profile%tile_min = MIN(profile%tile_min, count)
         profile%tile_max = MAX(profile%tile_max, count)
         started = m_walltime()
         IF (.NOT. cache_hit) THEN
            ALLOCATE (batch%axis_class(3, count), batch%combined_class(count), &
                      batch%axis_factor(MAXVAL(nmesh), count, 3))
            batch%axis_count = 0
            !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
            !$OMP SHARED(first, count, pair_shift, pair_tie, first_pair, ngrid, fractional_grid) &
            !$OMP PRIVATE(i, ilocal)
            DO i = 1, count
               ilocal = first + i - 1
               CALL product_pair_shift(first_pair + ilocal - 1, ngrid, fractional_grid, &
                                       pair_shift(:, i), pair_tie(:, i))
            END DO
            !$OMP END PARALLEL DO
            ! For source cell t, let R = wrap(t-s) and R+s = t+n*h. On k=m/n+offset,
            ! exp(2*pi*i*k*(R+s)) = exp(2*pi*i*m*t/n)*exp(2*pi*i*offset*(R+s)).
            ! Pack F(t) at its original cell with the latter phase and the original weight w(R).
            ! Each axis shares this input factor by exact (shift,tie), including large shifts.
            DO idir = 1, 3
               keys(1, 1:count) = pair_shift(idir, 1:count)
               keys(2, 1:count) = MERGE(1, 0, pair_tie(idir, 1:count))
               CALL sort(keys, 1, count, 1, 2, order)
               naxis = 0
               DO i = 1, count
                  IF (i > 1) THEN
                     IF (ALL(keys(1:2, i) == keys(1:2, i - 1))) THEN
                        batch%axis_class(idir, order(i)) = naxis
                        CYCLE
                     END IF
                  END IF
                  naxis = naxis + 1
                  batch%axis_class(idir, order(i)) = naxis
                  DO coordinate = 1, nmesh(idir)
                     icell = lower(idir) + MODULO(coordinate - 1 - keys(1, i), nmesh(idir))
                     weight = 1.0_dp
                     IF (MODULO(nmesh(idir), 2) == 0) THEN
                        IF (icell == lower(idir)) THEN
                           IF (keys(2, i) == 1) THEN
                              weight = weight*(1.0_dp + 3.0_dp*endpoint(idir))/4.0_dp
                           ELSE
                              weight = weight*(1.0_dp + endpoint(idir))/2.0_dp
                           END IF
                        ELSE IF (keys(2, i) == 1 .AND. icell == lower(idir) + 1) THEN
                           weight = weight*(3.0_dp + endpoint(idir))/4.0_dp
                        END IF
                     ELSE IF (keys(2, i) == 1 .AND. icell == lower(idir)) THEN
                        weight = weight*(1.0_dp + endpoint(idir))/2.0_dp
                     END IF
                     batch%axis_factor(coordinate, naxis, idir) = weight* &
                                                                  EXP(CMPLX(0.0_dp, twopi*k_offset(idir)* &
                                                                            (REAL(icell, dp) + REAL(keys(1, i), dp)), dp))
                  END DO
               END DO
               axis_classes = axis_classes + INT(naxis, int_8)
               phase_evaluations = phase_evaluations + INT(naxis, int_8)*INT(nmesh(idir), int_8)
            END DO
            ! Group the complete (x,y,z) stencil triples, retaining original pair order.
            keys(1:3, 1:count) = batch%axis_class(:, 1:count)
            CALL sort(keys, 1, count, 1, 3, order)
            combined_count = 0
            DO i = 1, count
               IF (i == 1) THEN
                  combined_count = 1
               ELSE IF (ANY(keys(:, i) /= keys(:, i - 1))) THEN
                  combined_count = combined_count + 1
               END IF
               batch%combined_class(order(i)) = combined_count
            END DO
            combined_classes = combined_classes + INT(combined_count, int_8)

            ! Bound cache memory to at most one quarter of the packed input.
            batch%use_combined = combined_count <= count/4
            IF (batch%use_combined) THEN
               ALLOCATE (batch%combined_factor(combined_count, SIZE(cell_index)))
               DO i = 1, count
                  IF (i > 1) THEN
                     IF (ALL(keys(:, i) == keys(:, i - 1))) CYCLE
                  END IF
                  combined_id = batch%combined_class(order(i))
                  DO icell = 1, SIZE(cell_index)
                     batch%combined_factor(combined_id, icell) = &
                        batch%axis_factor(cell_coordinate(1, icell), keys(1, i), 1)* &
                        batch%axis_factor(cell_coordinate(2, icell), keys(2, i), 2)* &
                        batch%axis_factor(cell_coordinate(3, icell), keys(3, i), 3)
                  END DO
               END DO
               cached_coefficients = cached_coefficients + &
                                     SIZE(cell_index, KIND=int_8)*INT(combined_count, int_8)
            END IF
            batch%axis_count = INT(axis_classes, KIND=KIND(batch%axis_count)) - &
                               SUM([(plan%batches(i)%axis_count, i=1, ibatch - 1)])
            batch%combined_count = combined_count
            IF (batch%use_combined) THEN
               batch%cached_count = INT(SIZE(cell_index), int_8)*INT(combined_count, int_8)
               ! The combined path no longer needs its intermediate one-axis tables.
               DEALLOCATE (batch%axis_class, batch%axis_factor)
            END IF
         ELSE
            axis_classes = axis_classes + INT(batch%axis_count, int_8)
            combined_classes = combined_classes + INT(batch%combined_count, int_8)
            cached_coefficients = cached_coefficients + batch%cached_count
         END IF
         IF (.NOT. cache_hit) setup_seconds = setup_seconds + m_walltime() - started
         input(1:INT(nfft_grid), 1:count) => fft_in(1:INT(nfft_grid)*count)
         started = m_walltime()
         ! The canonical cell list covers each mesh slot once. Only padded slots need zeros.
         IF (.NOT. unpadded) input = z_zero
         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
         !$OMP SHARED(input, values_R, cell_index, cell_coordinate, first, count, unpadded) &
         !$OMP SHARED(batch) &
         !$OMP PRIVATE(i, ilocal, icell, isource, class_x, class_y, class_z, row_first, row_last) &
         !$OMP PRIVATE(cell_first, cell_last, pack_buffer)
         DO row_first = 1, count, pack_pairs
            row_last = MIN(count, row_first + pack_pairs - 1)
            DO cell_first = 1, SIZE(cell_index), pack_cells
               cell_last = MIN(SIZE(cell_index), cell_first + pack_cells - 1)
               ! Adjacent pairs share source cache lines. The private tile keeps the subsequent
               ! transpose in cache while the destination is streamed one FFT grid at a time.
               DO icell = cell_first, cell_last
                  !$OMP SIMD PRIVATE(ilocal, class_x, class_y, class_z)
                  DO i = row_first, row_last
                     ilocal = first + i - 1
                     IF (batch%use_combined) THEN
                        pack_buffer(i - row_first + 1, icell - cell_first + 1) = &
                           batch%combined_factor(batch%combined_class(i), icell)*values_R(ilocal, 1, icell)
                     ELSE
                        class_x = batch%axis_class(1, i)
                        class_y = batch%axis_class(2, i)
                        class_z = batch%axis_class(3, i)
                        pack_buffer(i - row_first + 1, icell - cell_first + 1) = &
                           batch%axis_factor(cell_coordinate(1, icell), class_x, 1)* &
                           batch%axis_factor(cell_coordinate(2, icell), class_y, 2)* &
                           batch%axis_factor(cell_coordinate(3, icell), class_z, 3)*values_R(ilocal, 1, icell)
                     END IF
                  END DO
               END DO
               DO i = row_first, row_last
                  IF (unpadded) THEN
                     input(cell_first:cell_last, i) = pack_buffer(i - row_first + 1, 1:cell_last - cell_first + 1)
                  ELSE
                     DO icell = cell_first, cell_last
                        input(cell_index(icell), i) = pack_buffer(i - row_first + 1, icell - cell_first + 1)
                     END DO
                  END IF
               END DO
            END DO
         END DO
         !$OMP END PARALLEL DO
         profile%seconds(4) = profile%seconds(4) + m_walltime() - started
         CALL lattice_fft_execute(BWFFT, nfft, count, fft_in, fft_out, stat, profile)
         IF (stat /= 0) CPABORT('GW product FFT failed')
         output(1:count, 1:nfft(1), 1:nfft(2), 1:nfft(3)) => fft_in(1:INT(nfft_grid)*count)
         started = m_walltime()
         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
         !$OMP SHARED(values_k, output, grid_index, xkp, first, count, origin_phase) PRIVATE(ik, i)
         DO ik = 1, SIZE(xkp, 2)
            DO i = 1, count
               values_k(first + i - 1, 1, ik) = origin_phase(ik)* &
                                                output(i, grid_index(1, ik), grid_index(2, ik), grid_index(3, ik))
            END DO
         END DO
         !$OMP END PARALLEL DO
         profile%seconds(8) = profile%seconds(8) + m_walltime() - started
      END DO
      started = m_walltime()
      CALL fft_dealloc(fft_in)
      CALL fft_dealloc(fft_out)
      IF (.NOT. cache_hit) DEALLOCATE (pair_shift, pair_tie, keys, order)
      DEALLOCATE (cell_index, cell_coordinate, grid_index, origin_phase)
      plan%initialized = .TRUE.
      profile%seconds(2) = profile%seconds(2) + setup_seconds
      profile%seconds(9) = m_walltime() - started
      unit_nr = cp_logger_get_default_io_unit()
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') 'R2K_FFT_PROFILE| GW pair-centered packing/extraction; logger-rank wall seconds'
         WRITE (unit_nr, '(T2,A,3I10)') 'R2K_FFT_PROFILE| source shape: ', SHAPE(values_R)
         WRITE (unit_nr, '(T2,A,6I8)') 'R2K_FFT_PROFILE| grid / padded: ', nmesh, nfft
         WRITE (unit_nr, '(T2,A,2I16,2I8)') 'R2K_FFT_PROFILE| elements/batches/tile min/max: ', &
            profile%elements, profile%batches, profile%tile_min, profile%tile_max
         WRITE (unit_nr, '(T2,A,9F12.6)') 'R2K_FFT_PROFILE| init/geometry/buffers/fill/X/Y/Z/extract/release: ', profile%seconds
         WRITE (unit_nr, '(T2,A,F14.6)') 'R2K_FFT_PROFILE| class-setup (included in geometry): ', setup_seconds
         WRITE (unit_nr, '(T2,A,2I16)') 'R2K_FFT_PROFILE| input-axis phase evaluations / stencil classes: ', &
            phase_evaluations, axis_classes
         WRITE (unit_nr, '(T2,A,2I16)') 'R2K_FFT_PROFILE| 3D stencil classes / cached coefficients: ', &
            combined_classes, cached_coefficients
         WRITE (unit_nr, '(T2,A,L2)') 'R2K_FFT_PROFILE| reused stencil cache: ', cache_hit
         WRITE (unit_nr, '(T2,A,2I8)') 'R2K_FFT_PROFILE| packing tile pairs/cells: ', pack_pairs, pack_cells
         WRITE (unit_nr, '(T2,A,F14.6)') 'R2K_FFT_PROFILE| total         ', m_walltime() - total_started
      END IF
      CALL timestop(handle)

   END SUBROUTINE product_to_k_pair_centered

! **************************************************************************************************
!> \brief Select the lattice image of a grid-point pair and identify half-cell ties.
!> \param ipair flattened column-major pair index l+(l'-1)*ngrid
!> \param ngrid ...
!> \param fractional_grid ...
!> \param pair_shift ...
!> \param pair_tie ...
! **************************************************************************************************
   SUBROUTINE product_pair_shift(ipair, ngrid, fractional_grid, pair_shift, pair_tie)
      INTEGER, INTENT(IN)                                :: ipair, ngrid
      REAL(KIND=dp), INTENT(IN)                          :: fractional_grid(:, :)
      INTEGER, INTENT(OUT)                               :: pair_shift(3)
      LOGICAL, INTENT(OUT)                               :: pair_tie(3)

      INTEGER                                            :: idir, il, ilp
      REAL(KIND=dp)                                      :: delta(3), half_boundary, tolerance

      il = MODULO(ipair - 1, ngrid) + 1
      ilp = (ipair - 1)/ngrid + 1
      delta = fractional_grid(:, il) - fractional_grid(:, ilp)
      pair_shift = NINT(delta)
      pair_tie = .FALSE.
      DO idir = 1, 3
         half_boundary = ANINT(delta(idir) - 0.5_dp) + 0.5_dp
         tolerance = 64.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, ABS(delta(idir)))
         IF (ABS(delta(idir) - half_boundary) <= tolerance) THEN
            pair_shift(idir) = NINT(half_boundary - 0.5_dp)
            pair_tie(idir) = .TRUE.
         END IF
      END DO
   END SUBROUTINE product_pair_shift

! **************************************************************************************************
!> \brief Cosine transform of one q block between imaginary time and frequency,
!>        X_PQ(q,j) = Σ_i T_ji Y_PQ(q,i), as one GEMM over all (P,Q).
!> \param block Y on input, X on output
!> \param transform T, the real minimax cosine transform in either direction
! **************************************************************************************************
   SUBROUTINE cosine_transform_ri_q(block, transform)
      TYPE(ri_q_matrix_block_type), INTENT(INOUT)        :: block
      REAL(KIND=dp), INTENT(IN)                          :: transform(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cosine_transform_ri_q'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: values(:, :, :)
      INTEGER                                            :: handle, n_elements

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(transform, 2) == SIZE(block%values, 3))
      n_elements = SIZE(block%values, 1)*SIZE(block%values, 2)
      ALLOCATE (values(SIZE(block%values, 1), SIZE(block%values, 2), SIZE(transform, 1)))
      CALL zgemm('N', 'T', n_elements, SIZE(transform, 1), SIZE(transform, 2), z_one, block%values, &
                 n_elements, CMPLX(transform, KIND=dp), SIZE(transform, 1), z_zero, values, n_elements)
      CALL MOVE_ALLOC(values, block%values)

      CALL timestop(handle)

   END SUBROUTINE cosine_transform_ri_q

! **************************************************************************************************
!> \brief Count live complex elements across independently allocated q blocks.
!> \param blocks ...
!> \return ...
! **************************************************************************************************
   FUNCTION ri_q_storage_elements(blocks) RESULT(nelements)
      TYPE(ri_q_matrix_block_type), INTENT(IN)           :: blocks(:)
      INTEGER(KIND=int_8)                                :: nelements

      INTEGER                                            :: iq

      nelements = 0_int_8
      DO iq = 1, SIZE(blocks)
         IF (ALLOCATED(blocks(iq)%values)) THEN
            nelements = nelements + SIZE(blocks(iq)%values, KIND=int_8)
         END IF
      END DO
   END FUNCTION ri_q_storage_elements

! **************************************************************************************************
!> \brief Transfer q-local full matrices to owners of ordered grid pairs.
!> Ranks own q points cyclically, MOD(q,P), and ordered grid pairs as contiguous blocks. Round r
!> moves q = (r-1)P+1, ..., rP, the r-th local q of every rank, with one alltoall.
!> \param q_matrices ...
!> \param para_env ...
!> \param fft_pairs ...
! **************************************************************************************************
   SUBROUTINE q_matrices_to_fft_pairs(q_matrices, para_env, fft_pairs)
      COMPLEX(KIND=dp), INTENT(IN)                       :: q_matrices(:, :, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      COMPLEX(KIND=dp), INTENT(OUT)                      :: fft_pairs(:, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'q_matrices_to_fft_pairs'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: received(:), sent(:)
      INTEGER                                            :: handle, i_first, i_rank, i_round, &
                                                            n_grid, n_local, n_q, n_rank, n_round
      INTEGER, ALLOCATABLE                               :: pair_counts(:), pair_displs(:), &
                                                            recv_counts(:), recv_displs(:), &
                                                            send_counts(:)

      CALL timeset(routineN, handle)

      n_rank = para_env%num_pe
      n_grid = SIZE(q_matrices, 1)
      n_local = SIZE(fft_pairs, 1)
      n_q = SIZE(fft_pairs, 3)
      CPASSERT(SIZE(fft_pairs, 2) == 1)
      CALL grid_pair_blocks(n_grid, n_rank, pair_counts, pair_displs)
      ALLOCATE (sent(n_grid**2), received(n_local*n_rank), recv_counts(0:n_rank - 1), &
                recv_displs(0:n_rank - 1), send_counts(0:n_rank - 1))
      recv_displs(:) = [(MODULO(i_rank - 1, n_rank)*n_local, i_rank=0, n_rank - 1)]
      DO i_round = 1, (n_q + n_rank - 1)/n_rank
         i_first = (i_round - 1)*n_rank
         n_round = MIN(n_rank, n_q - i_first)
         DO i_rank = 0, n_rank - 1
            recv_counts(i_rank) = MERGE(n_local, 0, MODULO(i_rank - 1, n_rank) < n_round)
         END DO
         send_counts(:) = 0
         IF (MODULO(para_env%mepos - 1, n_rank) < n_round) THEN
            sent(:) = RESHAPE(q_matrices(:, :, i_round), [n_grid**2])
            send_counts(:) = pair_counts
         END IF
         CALL para_env%alltoall(sent, send_counts, pair_displs, received, recv_counts, recv_displs)
         fft_pairs(:, 1, i_first + 1:i_first + n_round) = RESHAPE(received(1:n_local*n_round), [n_local, n_round])
      END DO

      CALL timestop(handle)

   END SUBROUTINE q_matrices_to_fft_pairs

! **************************************************************************************************
!> \brief Transfer ordered grid pairs to the q owners as full matrices, the reverse of
!>        q_matrices_to_fft_pairs.
!> \param fft_pairs ...
!> \param para_env ...
!> \param needed q points to transfer
!> \param q_matrices needed q points of this rank, MOD(q,P), in ascending order
! **************************************************************************************************
   SUBROUTINE fft_pairs_to_q_matrices(fft_pairs, para_env, needed, q_matrices)
      COMPLEX(KIND=dp), INTENT(IN)                       :: fft_pairs(:, :, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: needed(:)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: q_matrices(:, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fft_pairs_to_q_matrices'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: received(:), sent(:)
      INTEGER                                            :: handle, i_first, i_local, i_rank, &
                                                            i_round, iq, n_grid, n_local, n_q, &
                                                            n_rank, n_round
      INTEGER, ALLOCATABLE                               :: pair_counts(:), pair_displs(:), &
                                                            recv_counts(:), send_counts(:), &
                                                            send_displs(:)
      LOGICAL                                            :: receiving

      CALL timeset(routineN, handle)

      n_rank = para_env%num_pe
      n_grid = SIZE(q_matrices, 1)
      n_local = SIZE(fft_pairs, 1)
      n_q = SIZE(fft_pairs, 3)
      CPASSERT(SIZE(fft_pairs, 2) == 1 .AND. SIZE(needed) == n_q)
      CALL grid_pair_blocks(n_grid, n_rank, pair_counts, pair_displs)
      ALLOCATE (sent(n_local*n_rank), received(n_grid**2), recv_counts(0:n_rank - 1), &
                send_counts(0:n_rank - 1), send_displs(0:n_rank - 1))
      send_displs(:) = [(MODULO(i_rank - 1, n_rank)*n_local, i_rank=0, n_rank - 1)]
      i_local = 0
      DO i_round = 1, (n_q + n_rank - 1)/n_rank
         i_first = (i_round - 1)*n_rank
         n_round = MIN(n_rank, n_q - i_first)
         send_counts(:) = 0
         DO i_rank = 0, n_rank - 1
            iq = i_first + MODULO(i_rank - 1, n_rank) + 1
            IF (iq > n_q) CYCLE
            IF (needed(iq)) send_counts(i_rank) = n_local
         END DO
         sent(1:n_local*n_round) = RESHAPE(fft_pairs(:, 1, i_first + 1:i_first + n_round), [n_local*n_round])
         iq = i_first + MODULO(para_env%mepos - 1, n_rank) + 1
         receiving = iq <= n_q
         IF (receiving) receiving = needed(iq)
         recv_counts(:) = 0
         IF (receiving) recv_counts(:) = pair_counts
         CALL para_env%alltoall(sent, send_counts, send_displs, received, recv_counts, pair_displs)
         IF (receiving) THEN
            i_local = i_local + 1
            q_matrices(:, :, i_local) = RESHAPE(received, [n_grid, n_grid])
         END IF
      END DO
      CPASSERT(i_local == SIZE(q_matrices, 3))

      CALL timestop(handle)

   END SUBROUTINE fft_pairs_to_q_matrices

! **************************************************************************************************
!> \brief Form W_ll'(q,iτ) = Σ_PQ Z_lP(q) W_PQ(q,iτ) Z_l'Q(q)* and transfer FFT pairs.
!> \param Z_q ...
!> \param interaction_q W_PQ(q,iτ) of the locally owned q
!> \param itime ...
!> \param para_env ...
!> \param fft_pairs ...
! **************************************************************************************************
   SUBROUTINE project_auxiliary_time_q_to_fft_pairs(Z_q, interaction_q, itime, para_env, fft_pairs)
      COMPLEX(KIND=dp), INTENT(IN)                       :: Z_q(:, :, :)
      TYPE(ri_q_matrix_block_type), INTENT(IN)           :: interaction_q(:)
      INTEGER, INTENT(IN)                                :: itime
      TYPE(mp_para_env_type), POINTER                    :: para_env
      COMPLEX(KIND=dp), INTENT(OUT)                      :: fft_pairs(:, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'project_auxiliary_time_q_to_fft_pairs'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: q_matrices(:, :, :), ZW(:, :)
      INTEGER                                            :: handle, i_q_local, n_grid, n_q_local, &
                                                            n_RI

      CALL timeset(routineN, handle)

      n_grid = SIZE(Z_q, 1)
      n_RI = SIZE(Z_q, 2)
      n_q_local = SIZE(Z_q, 3)
      CPASSERT(SIZE(interaction_q) == n_q_local)
      ALLOCATE (q_matrices(n_grid, n_grid, n_q_local))
      !$OMP PARALLEL DEFAULT(NONE) SHARED(Z_q, interaction_q, itime, q_matrices, n_q_local, n_grid, n_RI) &
      !$OMP PRIVATE(i_q_local, ZW)
      ALLOCATE (ZW(n_grid, n_RI))
      !$OMP DO SCHEDULE(STATIC)
      DO i_q_local = 1, n_q_local
         CALL zgemm('N', 'N', n_grid, n_RI, n_RI, z_one, Z_q(:, :, i_q_local), n_grid, &
                    interaction_q(i_q_local)%values(:, :, itime), n_RI, z_zero, ZW, n_grid)
         CALL zgemm('N', 'C', n_grid, n_grid, n_RI, z_one, ZW, n_grid, Z_q(:, :, i_q_local), n_grid, &
                    z_zero, q_matrices(:, :, i_q_local), n_grid)
      END DO
      !$OMP END DO
      DEALLOCATE (ZW)
      !$OMP END PARALLEL
      CALL q_matrices_to_fft_pairs(q_matrices, para_env, fft_pairs)

      CALL timestop(handle)

   END SUBROUTINE project_auxiliary_time_q_to_fft_pairs

! **************************************************************************************************
!> \brief ...
!> \param distributed ...
!> \param owner_matrix ...
!> \param group_env ...
!> \param root ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE set_owner_matrix(distributed, owner_matrix, group_env, root, n)
      TYPE(cp_cfm_type), INTENT(IN)                      :: distributed
      COMPLEX(KIND=dp), INTENT(IN)                       :: owner_matrix(:, :)
      TYPE(mp_para_env_type), INTENT(INOUT)              :: group_env
      INTEGER, INTENT(IN)                                :: root, n

      COMPLEX(KIND=dp), ALLOCATABLE                      :: column_block(:, :)
      INTEGER                                            :: first, width

      IF (group_env%mepos == root) THEN
         CPASSERT(ALL(SHAPE(owner_matrix) == [n, n]))
      END IF
      ALLOCATE (column_block(n, column_block_cols))
      DO first = 1, n, column_block_cols
         width = MIN(column_block_cols, n - first + 1)
         IF (group_env%mepos == root) column_block(:, 1:width) = owner_matrix(:, first:first + width - 1)
         CALL group_env%bcast(column_block(:, 1:width), root)
         CALL cp_cfm_set_submatrix(distributed, column_block(:, 1:width), &
                                   start_row=1, start_col=first)
      END DO
      DEALLOCATE (column_block)
   END SUBROUTINE set_owner_matrix

! **************************************************************************************************
!> \brief ...
!> \param distributed ...
!> \param owner_matrix ...
!> \param group_env ...
!> \param root ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE get_owner_matrix(distributed, owner_matrix, group_env, root, n)
      TYPE(cp_cfm_type), INTENT(IN)                      :: distributed
      COMPLEX(KIND=dp), INTENT(OUT)                      :: owner_matrix(:, :)
      TYPE(mp_para_env_type), INTENT(INOUT)              :: group_env
      INTEGER, INTENT(IN)                                :: root, n

      COMPLEX(KIND=dp), ALLOCATABLE                      :: column_block(:, :)
      INTEGER                                            :: first, width

      IF (group_env%mepos == root) THEN
         CPASSERT(ALL(SHAPE(owner_matrix) == [n, n]))
      END IF
      ALLOCATE (column_block(n, column_block_cols))
      DO first = 1, n, column_block_cols
         width = MIN(column_block_cols, n - first + 1)
         CALL cp_cfm_get_submatrix(distributed, column_block(:, 1:width), &
                                   start_row=1, start_col=first)
         IF (group_env%mepos == root) owner_matrix(:, first:first + width - 1) = column_block(:, 1:width)
      END DO
      DEALLOCATE (column_block)
   END SUBROUTINE get_owner_matrix

! **************************************************************************************************
!> \brief Find -q in the same mesh for every q; q is its own partner where -q ≡ q.
!> \param xkp q points in fractional coordinates
!> \param q_partner index of -q for every q
! **************************************************************************************************
   SUBROUTINE minus_q_partners(xkp, q_partner)
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: q_partner(:)

      INTEGER                                            :: ikp, jkp
      REAL(KIND=dp)                                      :: q_sum(3)

      ALLOCATE (q_partner(SIZE(xkp, 2)), SOURCE=0)
      DO ikp = 1, SIZE(xkp, 2)
         DO jkp = 1, SIZE(xkp, 2)
            q_sum(:) = xkp(:, ikp) + xkp(:, jkp)
            IF (ALL(ABS(q_sum - ANINT(q_sum)) < eps_kpoint)) THEN
               q_partner(ikp) = jkp
               EXIT
            END IF
         END DO
         CPASSERT(q_partner(ikp) > 0)
      END DO
   END SUBROUTINE minus_q_partners

! **************************************************************************************************
!> \brief Give the owner of every unsolved q the block of its partner, W(q,τ) = W(-q,τ)*.
!> \param para_env ...
!> \param q_partner index of -q for every q
!> \param solves_q whether W(q) was solved; then the owner of q holds its block
!> \param n_RI ...
!> \param n_time ...
!> \param blocks W(q,τ) of the q points owned by this rank
! **************************************************************************************************
   SUBROUTINE receive_minus_q_blocks(para_env, q_partner, solves_q, n_RI, n_time, blocks)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: q_partner(:)
      LOGICAL, INTENT(IN)                                :: solves_q(:)
      INTEGER, INTENT(IN)                                :: n_RI, n_time
      TYPE(ri_q_matrix_block_type), INTENT(INOUT)        :: blocks(:)

      INTEGER                                            :: destination, ikp, n_requests, &
                                                            partner_local, q_local, source
      TYPE(mp_request_type), ALLOCATABLE                 :: requests(:)

      ! A rank sends at most its solved and receives at most its unsolved blocks.
      ALLOCATE (requests(SIZE(blocks)))
      n_requests = 0
      DO ikp = 1, SIZE(q_partner)
         IF (.NOT. solves_q(ikp) .OR. q_partner(ikp) == ikp) CYCLE
         source = MODULO(ikp, para_env%num_pe)
         destination = MODULO(q_partner(ikp), para_env%num_pe)
         q_local = 1 + (ikp - 1)/para_env%num_pe
         partner_local = 1 + (q_partner(ikp) - 1)/para_env%num_pe
         IF (para_env%mepos == destination) THEN
            ALLOCATE (blocks(partner_local)%values(n_RI, n_RI, n_time))
            IF (source == destination) THEN
               blocks(partner_local)%values(:, :, :) = blocks(q_local)%values
            ELSE
               n_requests = n_requests + 1
               CALL para_env%irecv(blocks(partner_local)%values, source, requests(n_requests), tag=q_local)
            END IF
         ELSE IF (para_env%mepos == source) THEN
            n_requests = n_requests + 1
            CALL para_env%isend(blocks(q_local)%values, destination, requests(n_requests), tag=q_local)
         END IF
      END DO
      CALL mp_waitall(requests(1:n_requests))
      DO ikp = para_env%mepos, SIZE(q_partner), para_env%num_pe
         IF (ikp == 0) CYCLE
         IF (.NOT. solves_q(ikp)) THEN
            q_local = 1 + (ikp - 1)/para_env%num_pe
            blocks(q_local)%values(:, :, :) = CONJG(blocks(q_local)%values)
         END IF
      END DO
   END SUBROUTINE receive_minus_q_blocks

! **************************************************************************************************
!> \brief Partition MPI ranks into balanced contiguous q groups.
!> \param para_env ...
!> \param n_qgroups ...
!> \param my_qgroup ...
!> \param q_roots ...
! **************************************************************************************************
   SUBROUTINE grouped_w_partition(para_env, n_qgroups, my_qgroup, q_roots)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: n_qgroups
      INTEGER, INTENT(OUT)                               :: my_qgroup
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: q_roots(:)

      INTEGER                                            :: i_first, i_group, i_last

      CPASSERT(n_qgroups > 0 .AND. n_qgroups <= para_env%num_pe)
      ALLOCATE (q_roots(n_qgroups))
      DO i_group = 0, n_qgroups - 1
         i_first = (para_env%num_pe*i_group)/n_qgroups
         i_last = (para_env%num_pe*(i_group + 1))/n_qgroups - 1
         q_roots(i_group + 1) = i_first
         IF (para_env%mepos >= i_first .AND. para_env%mepos <= i_last) my_qgroup = i_group
      END DO

   END SUBROUTINE grouped_w_partition

! **************************************************************************************************
!> \brief ...
!> \param matrices ...
!> \param blacs ...
!> \param matrix_struct ...
!> \param group_env ...
!> \param n ...
!> Construct block-cyclic matrices on the MPI ranks assigned to a q group.
! **************************************************************************************************
   SUBROUTINE grouped_w_create(matrices, blacs, matrix_struct, group_env, n)
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: matrices(n_w_matrices)
      TYPE(cp_blacs_env_type), POINTER                   :: blacs
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(mp_para_env_type), INTENT(INOUT), TARGET      :: group_env
      INTEGER, INTENT(IN)                                :: n

      CPASSERT(n > 0)
      CALL cp_blacs_env_create(blacs, group_env)
      CALL cp_fm_struct_create(matrix_struct, context=blacs, &
                               nrow_global=n, ncol_global=n)
      CALL cp_cfm_create(matrices(w_chi), matrix_struct)
      CALL cp_cfm_create(matrices(w_metric), matrix_struct)
      CALL cp_cfm_create(matrices(w_tmp), matrix_struct)
      CALL cp_cfm_create(matrices(w_dielectric), matrix_struct)
      CALL cp_cfm_create(matrices(w_output), matrix_struct)
   END SUBROUTINE grouped_w_create

! **************************************************************************************************
!> \brief ...
!> \param world ...
!> \param n ...
!> \param source_rank ...
!> \param destination_rank ...
!> \param source ...
!> \param destination ...
!> Transfer a matrix between world owners without a full-size communication buffer.
! **************************************************************************************************
   SUBROUTINE grouped_w_copy_matrix(world, n, source_rank, destination_rank, source, destination)
      TYPE(mp_para_env_type), INTENT(INOUT)              :: world
      INTEGER, INTENT(IN)                                :: n, source_rank, destination_rank
      COMPLEX(KIND=dp), INTENT(IN)                       :: source(:, :)
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: destination(:, :)

      INTEGER                                            :: first, last, sender, tag

      IF (source_rank == destination_rank) THEN
         IF (world%mepos == source_rank) destination = source
         RETURN
      END IF
      DO first = 1, n, column_block_cols
         last = MIN(first + column_block_cols - 1, n)
         IF (world%mepos == source_rank) THEN
            CALL world%send(source(:, first:last), destination_rank, 0)
         ELSE IF (world%mepos == destination_rank) THEN
            sender = source_rank
            tag = 0
            CALL world%recv(destination(:, first:last), sender, tag)
         END IF
      END DO
   END SUBROUTINE grouped_w_copy_matrix

! **************************************************************************************************
!> \brief Evaluate screened W with a distributed positive-definite dielectric inverse.
!> For ε=I+BᴴχB after the unchanged positive χ filter, Cholesky gives ε^-1.
!> Reuse K=M^-1 V^1/2 to form W=K(ε^-1-I)Kᴴ.
!> \param matrices ...
!> \param group_env ...
!> \param n ...
!> \param root ...
!> \param chi_owner ...
!> \param cutoff ...
!> \param output_owner ...
! **************************************************************************************************
   SUBROUTINE grouped_w_solve(matrices, group_env, n, root, chi_owner, cutoff, output_owner)
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: matrices(n_w_matrices)
      TYPE(mp_para_env_type), INTENT(INOUT)              :: group_env
      INTEGER, INTENT(IN)                                :: n, root
      COMPLEX(KIND=dp), INTENT(IN)                       :: chi_owner(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      COMPLEX(KIND=dp), INTENT(OUT)                      :: output_owner(:, :)

      INTEGER                                            :: info

      CPASSERT(root >= 0 .AND. root < group_env%num_pe)
      CALL set_owner_matrix(matrices(w_chi), chi_owner, group_env, root, n)

      ! Match the dense reference: symmetrize before each Hermitian eigensolve.
      CALL cp_cfm_transpose(matrices(w_chi), 'C', matrices(w_output))
      CALL cp_cfm_scale_and_add(CMPLX(0.5_dp, 0.0_dp, KIND=dp), matrices(w_chi), &
                                CMPLX(0.5_dp, 0.0_dp, KIND=dp), matrices(w_output))
      CALL cp_cfm_power(matrices(w_chi), exponent=1.0_dp, threshold=cutoff)
      CALL parallel_gemm('C', 'N', n, n, n, z_one, matrices(w_metric), &
                         matrices(w_chi), z_zero, matrices(w_tmp))
      CALL parallel_gemm('N', 'N', n, n, n, z_one, matrices(w_tmp), &
                         matrices(w_metric), z_zero, matrices(w_dielectric))
      CALL cp_cfm_add_on_diag(matrices(w_dielectric), z_one)
      CALL cp_cfm_transpose(matrices(w_dielectric), 'C', matrices(w_output))
      CALL cp_cfm_scale_and_add(CMPLX(0.5_dp, 0.0_dp, KIND=dp), matrices(w_dielectric), &
                                CMPLX(0.5_dp, 0.0_dp, KIND=dp), matrices(w_output))
      ! χ is positive semidefinite after filtering: ε=I+BᴴχB has λ_min>=1.
      CALL cp_cfm_cholesky_decompose(matrices(w_dielectric), info_out=info)
      IF (info /= 0) THEN
         CPABORT("Screened dielectric Cholesky factorization failed")
      END IF
      CALL cp_cfm_cholesky_invert(matrices(w_dielectric), info_out=info)
      IF (info /= 0) CPABORT("Screened dielectric Cholesky inversion failed")
      CALL cp_cfm_uplo_to_full(matrices(w_dielectric), matrices(w_output), uplo='U')
      CALL cp_cfm_add_on_diag(matrices(w_dielectric), -z_one)
      ! With K=M^-1 V^1/2 and Hermitian M^-1, W=K(ε^-1-I)Kᴴ.
      ! The existing factor removes two matrix products at every later frequency.
      CALL parallel_gemm('N', 'N', n, n, n, z_one, matrices(w_metric), &
                         matrices(w_dielectric), z_zero, matrices(w_tmp))
      CALL parallel_gemm('N', 'C', n, n, n, z_one, matrices(w_tmp), &
                         matrices(w_metric), z_zero, matrices(w_output))
      CALL get_owner_matrix(matrices(w_output), output_owner, group_env, root, n)
   END SUBROUTINE grouped_w_solve

! **************************************************************************************************
!> \brief Release the distributed dielectric workspaces for a q-point group.
!> \param matrices Distributed matrix storage.
!> \param matrix_struct ...
!> \param blacs ...
! **************************************************************************************************
   SUBROUTINE grouped_w_release(matrices, matrix_struct, blacs)
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: matrices(n_w_matrices)
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_blacs_env_type), POINTER                   :: blacs

      CALL cp_cfm_release(matrices(w_chi))
      CALL cp_cfm_release(matrices(w_metric))
      CALL cp_cfm_release(matrices(w_tmp))
      CALL cp_cfm_release(matrices(w_dielectric))
      CALL cp_cfm_release(matrices(w_output))
      CALL cp_fm_struct_release(matrix_struct)
      CALL cp_blacs_env_release(blacs)
   END SUBROUTINE grouped_w_release

! **************************************************************************************************
!> \brief Release all periodic RI-RS work data after the last quasiparticle iteration.
!> \param bs_env Band-structure environment owning the work data.
! **************************************************************************************************
   SUBROUTINE release_gw_work(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'release_gw_work'

      INTEGER                                            :: handle
      TYPE(gwkp_type), POINTER                           :: gwkp_env

      CALL timeset(routineN, handle)

      gwkp_env => bs_env%gwkp_env
      CALL mp_para_env_release(gwkp_env%group_env)
      IF (ALLOCATED(gwkp_env%eigenvalues_G)) DEALLOCATE (gwkp_env%eigenvalues_G)
      ! AO image pointers are released immediately after the Gamma contraction.
      CPASSERT(.NOT. ASSOCIATED(gwkp_env%gamma_corr_d))
      CPASSERT(.NOT. ASSOCIATED(gwkp_env%gamma_corr_h))
      CPASSERT(.NOT. ASSOCIATED(gwkp_env%gamma_corr_s))
      CPASSERT(.NOT. ASSOCIATED(gwkp_env%gamma_corr_index_to_cell))
      CPASSERT(.NOT. ASSOCIATED(gwkp_env%gamma_corr_cell_to_index))
      DEALLOCATE (bs_env%gwkp_env)
      NULLIFY (gwkp_env)

      CALL timestop(handle)

   END SUBROUTINE release_gw_work

END MODULE gw_ri_rs_small_cell_full_kp_utils

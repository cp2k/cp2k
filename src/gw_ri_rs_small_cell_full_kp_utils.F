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
                                              cp_cfm_get_diag,&
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
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE fft_tools,                       ONLY: BWFFT,&
                                              fft_allocate,&
                                              fft_deallocate
   USE gw_ri_rs_small_cell_full_kp_lattice,ONLY: &
        assemble_AA_contraction, block_range, build_full_cell_grid, build_ri_2c_lattice, &
        complex_bytes, gamma_auxiliary_charges, grid_pair_gemm_panel, grid_pair_indices, &
        grid_pair_panel_cols, grid_pair_shift, grid_pair_stencil_denominator, &
        grid_pair_stencil_weight, grid_pair_upper_multiply, image_weight, owned_mo_coefficients, &
        real_bytes, ri_2c_lattice_to_kp, ri_metric_head_zero
   USE gw_ri_rs_small_cell_full_kp_q0_corr,ONLY: add_gamma_corr_pairs,&
                                              q_mesh_and_gamma_shell_average,&
                                              singular_head_cell_factors
   USE gw_ri_rs_small_cell_full_kp_types,ONLY: gwkp_type,&
                                              ri_q_matrix_block_type,&
                                              sigma_projection_type
   USE gw_utils_fm,                     ONLY: local_complex_power
   USE input_constants,                 ONLY: G0W0,&
                                              evGW0
   USE kinds,                           ONLY: dp,&
                                              int_1,&
                                              int_1_size,&
                                              int_8
   USE kpoint_coulomb_2c,               ONLY: add_ewald_reciprocal
   USE kpoint_lattice_fft,              ONLY: &
        fft_profile_enabled, k_grid_to_cell_execute, k_grid_to_cell_panel_width, &
        k_grid_to_cell_prepare, k_grid_to_cell_release, k_grid_to_cell_work_type, &
        lattice_fft_execute_pair_major, lattice_fft_profile_type, lattice_fft_real_sink_type
   USE kpoint_transform,                ONLY: cell_to_k_transform
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathconstants,                   ONLY: fourpi,&
                                              gaussi,&
                                              pi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE mathlib,                         ONLY: det_3x3,&
                                              inv_3x3
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
   USE ri_rs_kp_compute_Z_lP,           ONLY: compute_Z_lP_q,&
                                              evaluate_mu_rl_R
   USE ri_rs_utils,                     ONLY: precompute_ri_rs_radii,&
                                              write_ri_rs_progress

!$ USE omp_lib, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_small_cell_full_kp_utils'

   PUBLIC :: compute_Sigma_x, &
             compute_chi, &
             compute_exchange_energy_ri_rs, &
             compute_W0, &
             compute_sigma_c_from_w_pq, &
             prepare_ri_rs, &
             release_gw_work

   ! Tolerance for a sum of fractional k-point coordinates to be a reciprocal lattice vector.
   REAL(KIND=dp), PARAMETER :: eps_kpoint = 1.0E-8_dp
   ! Tolerance (bohr^-1) for |q + G| to be the Γ point.
   REAL(KIND=dp), PARAMETER :: eps_q_gamma = 1.0E-8_dp

   ! Blocking limits and named slots in the six existing CP2K matrix workspaces.
   INTEGER, PARAMETER :: column_block_cols = 32
   INTEGER, PARAMETER :: band_kpoints_per_pass = 512, pair_batch_size = 256, sigma_rhs_cols = 64
   INTEGER, PARAMETER :: w_chi = 1, w_metric = 2, w_tmp = 3, w_dielectric = 4, &
                         w_output = 5, n_w_matrices = 5

   ! Memory estimates: bytes per MiB and per GB (mp_mem_avail_per_rank_GB).
   REAL(KIND=dp), PARAMETER :: bytes_per_MiB = 1024.0_dp**2, bytes_per_GB = 1.0E9_dp
   ! A single stage may use stage_memory_share of the available memory. Stages that keep earlier
   ! data alive reserve reserve_share of it (at least min_reserve_GB) and add library_bytes.
   REAL(KIND=dp), PARAMETER :: stage_memory_share = 0.5_dp, reserve_share = 0.1_dp, &
                               min_reserve_GB = 1.0_dp, library_bytes = 256.0_dp*bytes_per_MiB
   ! Workspace limits for AO-image panels, received chi batches, and Σ(k) exchange rounds.
   REAL(KIND=dp), PARAMETER :: projection_workspace = 64.0_dp*bytes_per_MiB, &
                               projection_round_bytes = 512.0_dp*bytes_per_MiB
   INTEGER(KIND=int_8), PARAMETER :: fft_panel_bytes = 32_int_8*1024_int_8**2

   ! One stencil representation: signed bytes with denominator 64 and a shared cell phase.
   TYPE gw_pair_fft_batch_type
      INTEGER(KIND=int_1), ALLOCATABLE :: weight(:, :)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER :: input(:) => NULL()
   END TYPE gw_pair_fft_batch_type

   TYPE pair_segment_type
      INTEGER :: rank = -1, first_row = 0, last_row = 0, displacement = 0
   END TYPE pair_segment_type

   TYPE gw_pair_fft_work_type
      PRIVATE
      INTEGER :: n_pairs = 0, width = 0, k_batch = 0, green_reuses = 0, nmesh(3) = 0
      INTEGER(KIND=int_8) :: nfft_grid = 0_int_8, spectrum_columns = 0_int_8
      LOGICAL :: cache_hit = .FALSE.
      COMPLEX(KIND=dp), ALLOCATABLE :: product_phase(:)
      INTEGER, ALLOCATABLE :: grid_linear(:), partner_linear(:), extract_q(:), extract_partner(:), mirror_cell(:)
      INTEGER, ALLOCATABLE :: input_to_native(:), green_source(:), pair_counts(:), pair_displs(:), column_segments(:)
      INTEGER, ALLOCATABLE :: q_source(:), q_slot(:)
      TYPE(pair_segment_type), ALLOCATABLE :: segments(:)
      LOGICAL, ALLOCATABLE :: diagonal(:)
      TYPE(k_grid_to_cell_work_type) :: inverse
      TYPE(gw_pair_fft_batch_type), ALLOCATABLE :: batches(:)
      TYPE(lattice_fft_profile_type) :: profile = lattice_fft_profile_type()
      REAL(KIND=dp) :: setup_seconds = 0.0_dp
   END TYPE gw_pair_fft_work_type

   TYPE, EXTENDS(lattice_fft_real_sink_type) :: gw_pair_extraction_type
      LOGICAL :: product = .FALSE.
      REAL(KIND=dp) :: scale = 1.0_dp
      REAL(KIND=dp), ALLOCATABLE :: norms(:, :)
      REAL(KIND=dp), CONTIGUOUS, POINTER :: values(:, :, :) => NULL(), multiply(:, :, :) => NULL()
      TYPE(gw_pair_fft_work_type), POINTER :: fft => NULL()
   CONTAINS
      PROCEDURE :: consume => extract_grid_pairs
   END TYPE gw_pair_extraction_type

CONTAINS

! **************************************************************************************************
!> \brief Release the inverse FFT workspace and the aligned inputs of the product panels.
!> \param work ...
! **************************************************************************************************
   SUBROUTINE release_pair_fft(work)
      TYPE(gw_pair_fft_work_type), INTENT(INOUT)         :: work

      INTEGER                                            :: ibatch

      CALL k_grid_to_cell_release(work%inverse)
      IF (ALLOCATED(work%batches)) THEN
         DO ibatch = 1, SIZE(work%batches)
            IF (ASSOCIATED(work%batches(ibatch)%input)) CALL fft_deallocate(work%batches(ibatch)%input)
         END DO
         DEALLOCATE (work%batches, work%product_phase, work%grid_linear, work%partner_linear, &
                     work%extract_q, work%extract_partner, work%mirror_cell, work%diagonal, work%input_to_native, &
                     work%green_source, work%pair_counts, work%pair_displs, work%column_segments, work%segments)
      END IF
      IF (ALLOCATED(work%q_source)) DEALLOCATE (work%q_source, work%q_slot)
      work%profile = lattice_fft_profile_type()
      work%setup_seconds = 0.0_dp
      work%cache_hit = .FALSE.
      work%green_reuses = 0
   END SUBROUTINE release_pair_fft

! **************************************************************************************************
!> \brief Consume K2R spectra once: original Green norms, physical output and next-FFT input.
!> Chi stores its first Green in the input panel and replaces it in place with the weighted
!> product on the second call. Self-energy retains the physical product for band interpolation.
!> \param sink extraction recipe and result
!> \param first first local pair in this FFT panel
!> \param nbatch physical pairs in this panel
!> \param spectrum pair-major FFT output
!> \param cell_slot FFT slot of each physical cell
!> \param cell_phase common inverse-transform phase and normalization
! **************************************************************************************************
   SUBROUTINE extract_grid_pairs(sink, first, nbatch, spectrum, cell_slot, cell_phase)
      CLASS(gw_pair_extraction_type), INTENT(INOUT) :: sink
      INTEGER(KIND=int_8), INTENT(IN) :: first
      INTEGER, INTENT(IN) :: nbatch, cell_slot(:)
      COMPLEX(KIND=dp), INTENT(IN) :: spectrum(:, :), cell_phase(:)

      COMPLEX(KIND=dp) :: input_phase, phase, raw, saved
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER :: input(:, :)
      INTEGER :: fft_count, i, ia, ib, icell, ipair, jpair, slot
      REAL(KIND=dp) :: a, b, norm, offnorm, product_a, product_b
      LOGICAL :: collect_norm
      TYPE(gw_pair_fft_batch_type), POINTER :: batch
      TYPE(gw_pair_fft_work_type), POINTER :: plan

      fft_count = (nbatch + 1)/2
      IF (.NOT. sink%two_real) THEN
         CPASSERT(nbatch == 1)
      END IF
      NULLIFY (input, batch)
      collect_norm = ALLOCATED(sink%norms)
      plan => sink%fft
      IF (ASSOCIATED(plan)) THEN
         CPASSERT(MODULO(INT(first) - 1, plan%width) == 0)
         batch => plan%batches((INT(first) - 1)/plan%width + 1)
         input(1:fft_count, 1:INT(plan%nfft_grid)) => batch%input
      END IF
      !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
      !$OMP SHARED(sink, first, nbatch, fft_count, spectrum, cell_slot, cell_phase, plan, batch, input) &
      !$OMP SHARED(collect_norm) &
      !$OMP PRIVATE(icell, slot, phase, i, ia, ib, ipair, jpair, raw, saved, A, b, product_a, product_b, norm, offnorm) &
      !$OMP PRIVATE(input_phase)
      DO icell = 1, SIZE(cell_slot)
         slot = cell_slot(icell)
         phase = MERGE(0.5_dp, 1.0_dp, sink%two_real)*cell_phase(icell)
         norm = 0.0_dp
         offnorm = 0.0_dp
         IF (sink%product) input_phase = plan%product_phase(icell)
         !$OMP SIMD REDUCTION(+:norm, offnorm) &
         !$OMP PRIVATE(ia, ib, ipair, jpair, raw, saved, A, b, product_a, product_b)
         DO i = 1, fft_count
            ia = 2*i - 1
            ib = MIN(2*i, nbatch)
            ipair = INT(first) + ia - 1
            jpair = INT(first) + ib - 1
            raw = phase*spectrum(i, slot)
            a = REAL(raw, dp)
            b = 0.0_dp
            IF (ib > ia) b = AIMAG(raw)
            IF (collect_norm) THEN
               norm = norm + a*a + b*b
               IF (.NOT. plan%diagonal(ipair)) offnorm = offnorm + a*a
               IF (ib > ia) THEN
                  IF (.NOT. plan%diagonal(jpair)) offnorm = offnorm + b*b
               END IF
            END IF
            product_a = a
            product_b = b
            IF (ASSOCIATED(sink%multiply)) THEN
               product_a = (sink%scale*a)*sink%multiply(ipair, 1, icell)
               IF (ib > ia) product_b = (sink%scale*b)*sink%multiply(jpair, 1, icell)
            ELSE IF (sink%product) THEN
               saved = input(i, slot)
               product_a = (sink%scale*REAL(saved, dp))*a
               product_b = (sink%scale*AIMAG(saved))*b
            END IF
            IF (ASSOCIATED(sink%values)) THEN
               sink%values(ipair, 1, icell) = product_a
               IF (ib > ia) sink%values(jpair, 1, icell) = product_b
            END IF
            IF (.NOT. ASSOCIATED(plan)) CYCLE
            IF (.NOT. sink%product) THEN
               input(i, slot) = CMPLX(a, b, dp)
               CYCLE
            END IF
            input(i, slot) = input_phase* &
                             CMPLX(REAL(batch%weight(ia, icell), dp)*product_a, REAL(batch%weight(ib, icell), dp)*product_b, dp)
         END DO
         IF (collect_norm) THEN
            sink%norms(1, icell) = sink%norms(1, icell) + norm
            sink%norms(2, icell) = sink%norms(2, icell) + offnorm
         END IF
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE extract_grid_pairs

! **************************************************************************************************
!> \brief Single GW K2R pipeline for unordered pairs and fused product consumers.
!> The squared Frobenius norm counts each upper pair at R and its missing lower partner at -R;
!> a diagonal is counted once. The MP boundary twist drops out of this squared norm.
!> \param values_k unordered pairs on the input mesh
!> \param work shared lattice geometry and next-FFT input
!> \param values_R optional physical field or product required by a later consumer
!> \param multiply_R second physical field for the self-energy product
!> \param product_scale finish the R2K stencil using multiply_R or the saved chi Green
!> \param norm_squared original full-matrix squared norm per cell
! **************************************************************************************************
   SUBROUTINE k_to_grid_pairs(values_k, work, values_R, multiply_R, product_scale, norm_squared)
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(IN)           :: values_k(:, :, :)
      TYPE(gw_pair_fft_work_type), INTENT(INOUT), TARGET :: work
      REAL(KIND=dp), CONTIGUOUS, INTENT(INOUT), &
         OPTIONAL, TARGET                                :: values_R(:, :, :)
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN), OPTIONAL, &
         TARGET                                          :: multiply_R(:, :, :)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: product_scale
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: norm_squared(:)

      CHARACTER(LEN=14), PARAMETER :: stages(9) = [CHARACTER(LEN=14) :: 'initialize', 'geometry', &
         'buffers', 'grid-fill', 'FFT-X/plan', 'FFT-Y/plan', 'FFT-Z/plan', 'extract', 'release']

      INTEGER                                            :: i, n_cells, unit_nr
      REAL(KIND=dp)                                      :: total_started
      TYPE(gw_pair_extraction_type)                      :: sink
      TYPE(lattice_fft_profile_type), ALLOCATABLE        :: profile

      IF (fft_profile_enabled) total_started = m_walltime()
      CPASSERT(ALLOCATED(work%batches))
      CPASSERT(SIZE(values_k, 1) == work%n_pairs .AND. SIZE(values_k, 2) == 1)
      n_cells = SIZE(work%inverse%cell_slot)
      sink%product = PRESENT(product_scale)
      IF (PRESENT(product_scale)) sink%scale = product_scale
      IF (sink%product .OR. .NOT. PRESENT(values_R)) sink%fft => work
      IF (PRESENT(values_R)) THEN
         CPASSERT(ALL(SHAPE(values_R) == [work%n_pairs, 1, n_cells]))
         sink%values => values_R
      END IF
      IF (PRESENT(multiply_R)) THEN
         CPASSERT(sink%product)
         CPASSERT(ALL(SHAPE(multiply_R) == [work%n_pairs, 1, n_cells]))
         sink%multiply => multiply_R
      END IF
      IF (PRESENT(norm_squared)) THEN
         CPASSERT(ASSOCIATED(sink%fft))
         CPASSERT(SIZE(norm_squared) == n_cells)
         ALLOCATE (sink%norms(2, n_cells), SOURCE=0.0_dp)
      END IF
      IF (fft_profile_enabled) ALLOCATE (profile)
      CALL k_grid_to_cell_execute(work%inverse, values_k, require_fft=.TRUE., profile=profile, &
                                  max_workspace_bytes=fft_panel_bytes, real_sink=sink)
      IF (PRESENT(norm_squared)) THEN
         norm_squared = sink%norms(1, :) + sink%norms(2, work%mirror_cell)
      END IF
      unit_nr = 0
      IF (fft_profile_enabled) unit_nr = cp_logger_get_default_io_unit()
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') 'K2R_FFT_PROFILE| Unordered grid pairs; logger-rank wall seconds, no MPI reduction'
         WRITE (unit_nr, '(T2,A,3I10)') 'K2R_FFT_PROFILE| source shape: ', SHAPE(values_k)
         WRITE (unit_nr, '(T2,A,6I8)') 'K2R_FFT_PROFILE| grid / padded: ', work%nmesh, profile%nfft
         WRITE (unit_nr, '(T2,A,L2)') 'K2R_FFT_PROFILE| two real outputs per complex FFT: ', profile%two_real
         WRITE (unit_nr, '(T2,A,L2)') 'K2R_FFT_PROFILE| extraction feeds next FFT input: ', ASSOCIATED(sink%fft)
         WRITE (unit_nr, '(T2,A,2I16)') 'K2R_FFT_PROFILE| FFT transforms / real outputs: ', profile%transforms, profile%elements
         WRITE (unit_nr, '(T2,A,2I16)') 'K2R_FFT_PROFILE| complex FFT elements / full: ', &
            PRODUCT(INT(profile%nfft, int_8))*profile%transforms, PRODUCT(INT(profile%nfft, int_8))*profile%elements
         WRITE (unit_nr, '(T2,A,2I16,2I8)') 'K2R_FFT_PROFILE| elements/batches/tile min/max: ', &
            profile%elements, profile%batches, profile%tile_min, profile%tile_max
         WRITE (unit_nr, '(T2,A,L2)') 'K2R_FFT_PROFILE| pair-major FFT input/output: ', profile%pair_major
         WRITE (unit_nr, '(T2,A)') 'K2R_FFT_PROFILE| FFT-X/plan reports full 3D batched FFT (Y/Z entries zero)'
         DO i = 1, SIZE(stages)
            WRITE (unit_nr, '(T2,A,A14,F14.6)') 'K2R_FFT_PROFILE| ', stages(i), profile%seconds(i)
         END DO
         WRITE (unit_nr, '(T2,A,F14.6)') 'K2R_FFT_PROFILE| total         ', m_walltime() - total_started
      END IF
   END SUBROUTINE k_to_grid_pairs

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
      CALL minus_q_partners(bs_env%kpoints_chi_eps_W%xkp(:, :bs_env%nkp_scf_desymm), gwkp_env%q_partner)
      CALL compact_q_slots(gwkp_env%q_partner, bs_env%para_env%num_pe, gwkp_env%q_slot)
      CALL compute_Z_lP_q(bs_env, bs_env%ri_rs%grid_points, gwkp_env%q_slot, &
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
!> \brief Compute ε, Wᶜ and its imaginary-time transform, Eqs. (7)-(8), and the RPA correlation
!>        energy per unit cell from the same Cholesky factors of ε(q,iω):
!>        E_c^RPA = 1/(2π N_k) Σ_q ∫_0^∞ dω Tr[ln ε(q,iω) - ε(q,iω) + I],
!>        stored in gwkp_env%e_c_rpa. The minimax weights w_j integrate 4 ∫_0^∞ dω (normalization of get_rpa_minimax_coeff),
!>        and ε(-q) = ε(q)* contributes the same trace as ε(q). Wᶜ is only computed for GW.
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
      INTEGER, ALLOCATABLE                               :: group_load(:), q_batch(:, :), q_roots(:)
      REAL(KIND=dp)                                      :: available_GB, required_bytes, &
                                                            rpa_integrand, step_start
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
      nkp_local = SIZE(gwkp_env%metric_coulomb_factor)
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
      ! Store W only at the first q of each inversion pair. Its conjugate is
      ! reconstructed after projecting to grid pairs in the correlation step.
      ALLOCATE (group_load(0:n_qgroups - 1), SOURCE=0)
      DO ikp = 1, nkp
         IF (gwkp_env%q_slot(ikp) > 0) group_load(MODULO(ikp, n_qgroups)) = group_load(MODULO(ikp, n_qgroups)) + 1
      END DO
      ALLOCATE (q_batch(MAXVAL(group_load), 0:n_qgroups - 1), SOURCE=0)
      group_load(:) = 0
      DO ikp = 1, nkp
         IF (gwkp_env%q_slot(ikp) == 0) CYCLE
         q_group = MODULO(ikp, n_qgroups)
         group_load(q_group) = group_load(q_group) + 1
         q_batch(group_load(q_group), q_group) = ikp
      END DO
      ! RPA needs only E_c^RPA from the Cholesky factor of ε, not W
      IF (bs_env%do_gw) ALLOCATE (gwkp_env%W_PQ_tau(nkp_local))
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
            ikp_local = gwkp_env%q_slot(ikp)
            IF (bs_env%para_env%mepos == q_owner) THEN
               CALL grouped_w_copy_matrix(bs_env%para_env, n_RI, q_owner, q_root, &
                                          gwkp_env%metric_coulomb_factor(ikp_local)%values(:, :, 1), chi_k_w)
               IF (bs_env%do_gw) ALLOCATE (gwkp_env%W_PQ_tau(ikp_local)%values( &
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
               ikp_local = gwkp_env%q_slot(ikp)
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
            ikp = q_batch(i_batch, my_qgroup)
            IF (ikp > 0) THEN
               CALL grouped_w_solve(w_matrices, w_group_env, n_RI, 0, chi_k_w, &
                                    bs_env%eps_eigval_mat_RI, MWM_k_w, rpa_integrand, bs_env%do_gw)
               IF (w_group_env%mepos == 0) THEN
                  gwkp_env%e_c_rpa = gwkp_env%e_c_rpa + MERGE(1.0_dp, 2.0_dp, gwkp_env%q_partner(ikp) == ikp)* &
                                     bs_env%time_frequency_grid%frequency_weights(j_w)*rpa_integrand/ &
                                     (8.0_dp*pi*REAL(nkp, dp))
               END IF
            END IF
            DO q_group = 0, n_qgroups - 1
               ikp = q_batch(i_batch, q_group)
               IF (ikp == 0) CYCLE
               q_owner = MODULO(ikp, bs_env%para_env%num_pe)
               ikp_local = gwkp_env%q_slot(ikp)
               IF (bs_env%do_gw) THEN
                  CALL grouped_w_copy_matrix(bs_env%para_env, n_RI, q_roots(q_group + 1), q_owner, MWM_k_w, chi_k_w)
                  IF (bs_env%para_env%mepos == q_owner) gwkp_env%W_PQ_tau(ikp_local)%values(:, :, j_w) = chi_k_w
               END IF
               IF (bs_env%para_env%mepos == q_owner .AND. j_w == bs_env%num_time_freq_points) THEN
                  DEALLOCATE (gwkp_env%chi_q_frequency(ikp_local)%values)
               END IF
            END DO
            frequency_elapsed(j_w) = frequency_elapsed(j_w) + m_walltime() - step_start
         END DO
         DO q_group = 0, n_qgroups - 1
            ikp = q_batch(i_batch, q_group)
            IF (ikp == 0 .OR. .NOT. bs_env%do_gw) CYCLE
            IF (bs_env%para_env%mepos /= MODULO(ikp, bs_env%para_env%num_pe)) CYCLE
            CALL cosine_transform_ri_q(gwkp_env%W_PQ_tau(gwkp_env%q_slot(ikp)), &
                                       bs_env%time_frequency_grid%cosine_frequency_to_time)
         END DO
      END DO
      CALL bs_env%para_env%sum(gwkp_env%e_c_rpa)
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ''
      DO j_w = 1, bs_env%num_time_freq_points
         CALL bs_env%para_env%max(frequency_elapsed(j_w))
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I3,A,I0,T56,A,T70,F10.1,A)') &
               TRIM(MERGE('GW W(iω) point   ', 'RPA ε(iω) point ', bs_env%do_gw))//' ', j_w, ' / ', &
               bs_env%num_time_freq_points, &
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

      DEALLOCATE (group_load, q_batch, gwkp_env%leak_e, gwkp_env%leak_x, gwkp_env%leak_l0)
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
      CALL prepare_sigma_projection(bs_env, bs_env%para_env, gwkp_env%band_orbitals, gwkp_env%sigma_projection)
      CALL compute_Sigma_x_from_V(bs_env, gwkp_env)
      CALL prepare_self_energy_groups(bs_env, gwkp_env)
      IF (gwkp_env%ngroups > 1) THEN
         DEALLOCATE (gwkp_env%sigma_projection%partner_orbitals)
         CALL band_orbitals_to_groups(bs_env, gwkp_env)
         CALL prepare_sigma_projection(bs_env, gwkp_env%group_env, gwkp_env%band_orbitals, gwkp_env%sigma_projection)
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x

! **************************************************************************************************
!> \brief Retain Z and fixed W in the smallest memory-feasible time groups.
!>
!> Each group owns τ = group + 1, group + 1 + ngroups, ... . Within a group,
!> Z and W retain the representative q points, in native order within owner MOD(q,g).
!> \param bs_env Band-structure settings.
!> \param gwkp_env Work arrays and their MPI distribution.
! **************************************************************************************************
   SUBROUTINE prepare_self_energy_groups(bs_env, gwkp_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT)                     :: gwkp_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_self_energy_groups'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: Z_group(:, :, :)
      INTEGER :: handle, i_dest, i_group_size, i_root, i_source, i_tag, iq, iq_local, it, &
         it_local, n_ao, n_aux, n_grid, n_images, n_jobs, n_k, n_local_k, n_local_q, n_local_time, &
         n_offmesh, n_q_exchange, n_spin, n_stored_q, n_threads, n_time, n_world
      REAL(KIND=dp) :: available_GB, block_bytes, budget_bytes, grid_bytes, input_bytes, n_pairs, &
         peak_bytes, persistent_bytes, result_bytes, sigma_geometry_bytes
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
      n_jobs = SUM(gwkp_env%sigma_projection%counts)
      n_offmesh = COUNT(gwkp_env%sigma_projection%mesh_index == 0)
      n_images = PRODUCT(bs_env%kpoints_scf_desymm%nkp_grid + 2)
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
         n_q_exchange = pair_exchange_width(n_grid, n_k, i_group_size)
         n_local_time = (n_time + gwkp_env%ngroups - 1)/gwkp_env%ngroups
         n_stored_q = 0
         DO iq = 1, MIN(n_k, i_group_size)
            n_stored_q = MAX(n_stored_q, COUNT(gwkp_env%q_slot(iq:n_k:i_group_size) > 0))
         END DO
         persistent_bytes = complex_bytes*(REAL(n_stored_q, dp)*REAL(n_grid, dp)*REAL(n_aux, dp) + &
                                           REAL(n_aux, dp)**2*REAL(n_stored_q*n_local_time, dp))
         n_pairs = REAL(CEILING(REAL(n_grid, dp)*REAL(n_grid + 1, dp)/(2.0_dp*REAL(i_group_size, dp)), KIND=int_8), dp)
         sigma_geometry_bytes = 0.0_dp
         IF (n_offmesh > 0) THEN
            sigma_geometry_bytes = 3.0_dp*n_pairs*(real_bytes + REAL(int_1_size + STORAGE_SIZE(0)/8, dp)) + &
                                   complex_bytes*REAL(n_offmesh, dp)*REAL(n_images, dp)
         END IF
         ! G(k), W(R), the retained off-mesh product, and the prepared FFT input coexist.
         ! Count the band orbitals of the owned band k points alongside these panels.
         ! q matrices, band projection, and Green blocks belong to separate stages.
         grid_bytes = (n_pairs*(complex_bytes + 3.0_dp*real_bytes) + real_bytes)*REAL(n_k, dp) + &
                      REAL(int_1_size, dp)*n_pairs*REAL(n_k, dp) + &
                      MAX(REAL(fft_panel_bytes, dp), 3.0_dp*complex_bytes*REAL(n_k, dp)) + &
                      real_bytes*(n_pairs + REAL(n_k, dp)) + &
                      complex_bytes*REAL(n_grid, dp)*REAL(n_ao, dp)*REAL(n_spin, dp)* &
                      REAL((bs_env%nkp_bs_and_DOS + i_group_size - 1)/i_group_size + n_stored_q, dp) + &
                      MAX(complex_bytes*REAL(n_grid, dp)*REAL(n_grid + 1, dp)*REAL(n_local_k, dp)/2.0_dp + &
                          complex_bytes*REAL(n_aux, dp)**2*REAL(n_local_k, dp) + &
                          complex_bytes*REAL(n_threads, dp)*REAL(n_grid, dp)* &
                          REAL(n_ao + n_aux + MIN(n_grid, grid_pair_panel_cols), dp) + &
                          block_bytes*REAL(MAX(0, n_q_exchange - n_threads), dp), &
                          projection_workspace_bytes(n_grid, i_group_size, bs_env%nkp_bs_and_DOS, &
                                                     n_jobs, n_offmesh, n_images, n_ao, n_threads)) + &
                      2.0_dp*complex_bytes*REAL(n_threads*n_k, dp) + sigma_geometry_bytes + library_bytes
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
      CALL column_env%from_split(bs_env%para_env, gwkp_env%group_env%mepos)
      CALL gather_group_columns(bs_env, gwkp_env, column_env, i_group_size, gwkp_env%group_env%mepos, Z_group)
      CALL MOVE_ALLOC(Z_group, gwkp_env%Z_lP_q)
      n_local_q = SIZE(gwkp_env%Z_lP_q, 3)
      n_local_time = (n_time - gwkp_env%my_group + gwkp_env%ngroups - 1)/gwkp_env%ngroups
      ALLOCATE (W_group(n_local_q))
      iq_local = 0
      DO iq = MODULO(gwkp_env%group_env%mepos - 1, i_group_size) + 1, n_k, i_group_size
         IF (gwkp_env%q_slot(iq) == 0) CYCLE
         i_root = MODULO(iq, n_world)/i_group_size
         iq_local = iq_local + 1
         ALLOCATE (W_group(iq_local)%values(n_aux, n_aux, n_local_time))
         DO it = 1, n_time
            i_dest = MODULO(it - 1, gwkp_env%ngroups)
            IF (column_env%mepos == i_root) THEN
               IF (i_dest == i_root) THEN
                  it_local = (it - 1)/gwkp_env%ngroups + 1
                  W_group(iq_local)%values(:, :, it_local) = gwkp_env%W_PQ_tau(gwkp_env%q_slot(iq))%values(:, :, it)
               ELSE
                  CALL column_env%send(gwkp_env%W_PQ_tau(gwkp_env%q_slot(iq))%values(:, :, it), i_dest, 0)
               END IF
            ELSE IF (column_env%mepos == i_dest) THEN
               i_source = i_root
               i_tag = 0
               it_local = (it - 1)/gwkp_env%ngroups + 1
               CALL column_env%recv(W_group(iq_local)%values(:, :, it_local), i_source, i_tag)
            END IF
         END DO
         IF (column_env%mepos == i_root) DEALLOCATE (gwkp_env%W_PQ_tau(gwkp_env%q_slot(iq))%values)
      END DO
      CALL column_env%free()
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

      COMPLEX(KIND=dp), ALLOCATABLE                      :: psi_block(:, :, :)
      INTEGER                                            :: handle, i_root, ik, ik_local, iq, &
                                                            iq_first, iq_local, n_grid, n_k, &
                                                            n_world

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
      iq_first = MODULO(group_rank - 1, group_size) + 1
      ALLOCATE (Z_group(n_grid, bs_env%n_RI, COUNT(gwkp_env%q_slot(iq_first:n_k:group_size) > 0)))
      iq_local = 0
      DO iq = iq_first, n_k, group_size
         IF (gwkp_env%q_slot(iq) == 0) CYCLE
         i_root = MODULO(iq, n_world)/group_size
         iq_local = iq_local + 1
         IF (column_env%mepos == i_root) Z_group(:, :, iq_local) = gwkp_env%Z_lP_q(:, :, gwkp_env%q_slot(iq))
         CALL column_env%bcast(Z_group(:, :, iq_local), i_root)
      END DO

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

      COMPLEX(KIND=dp), ALLOCATABLE                      :: chi_grid_q(:, :, :), chi_q(:, :, :), &
                                                            G_grid_k(:, :, :)
      COMPLEX(KIND=dp), ALLOCATABLE, TARGET              :: psi_time(:, :, :, :), Z_group(:, :, :)
      COMPLEX(KIND=dp), POINTER                          :: psi_values(:, :, :, :), Z_time(:, :, :)
      INTEGER :: first_pair, handle, i_aux, i_group, i_needed, i_q_world, i_rank, i_root, i_size, &
         i_source, i_t, i_tag, i_time_group, i_time_local, iq, iq_local, ispin, n_aux, n_cells, &
         n_grid, n_groups, n_local, n_needed_local, n_pairs_local, n_threads, n_time_local, &
         n_world, q_batch, unit_nr
      INTEGER(KIND=int_8)                                :: n_pairs_total
      INTEGER, ALLOCATABLE                               :: full_cells(:, :), needed_q(:)
      LOGICAL, ALLOCATABLE                               :: cell_active(:), needed(:)
      REAL(KIND=dp) :: available_GB, budget, chi_seconds(7), chi_started, green_bytes, &
         input_bytes, live_bytes, n_pair_max, n_q_max, peak, q_bytes, rows_bytes, stage_started, &
         tau, time_bytes
      REAL(KIND=dp), ALLOCATABLE                         :: green_norms(:, :)
      TYPE(gw_pair_fft_work_type)                        :: fft_stencil_work
      TYPE(mp_comm_type)                                 :: column_env
      TYPE(mp_para_env_type), POINTER                    :: chi_env
      TYPE(ri_q_matrix_block_type), ALLOCATABLE          :: chi_partial(:)

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(gwkp_env%group_env))
      gwkp_env%group_env => bs_env%para_env
      CALL gwkp_env%group_env%retain()
      CALL build_full_cell_grid(bs_env, full_cells)
      n_grid = bs_env%ri_rs%n_grid_points
      n_pairs_total = INT(n_grid, int_8)*(n_grid + 1)/2
      n_aux = bs_env%n_RI
      n_world = gwkp_env%group_env%num_pe
      ALLOCATE (needed(bs_env%nkp_scf_desymm))
      needed(:) = gwkp_env%q_slot > 0
      ALLOCATE (needed_q(COUNT(needed)))
      needed_q(:) = PACK([(iq, iq=1, SIZE(needed))], needed)
      CALL mp_mem_avail_per_rank_GB(gwkp_env%group_env, available_GB)
      budget = budget_after_reserve(available_GB)
      n_cells = SIZE(full_cells, 2)
      n_threads = 1
!$    n_threads = omp_get_max_threads()
      CALL gwkp_env%group_env%max(n_threads)
      ! Keep received q matrices near their consumer; cap a batch at 64 MiB, or one matrix.
      q_batch = MAX(1, MIN(4*n_threads, INT(projection_workspace/(complex_bytes*REAL(n_pairs_total, dp)))))
      DO i_size = MIN(2, MAX(1, n_world/2)), n_world
         IF (MODULO(n_world, i_size) /= 0) CYCLE
         n_groups = n_world/i_size
         IF (n_groups > bs_env%num_time_freq_points) CYCLE
         i_rank = MODULO(gwkp_env%group_env%mepos, i_size)
         i_group = gwkp_env%group_env%mepos/i_size
         n_time_local = (bs_env%num_time_freq_points + n_groups - 1 - i_group)/n_groups
         n_pair_max = REAL(n_pairs_total*(i_rank + 1)/i_size - n_pairs_total*i_rank/i_size, dp)
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
         ! The occupied K2R stores directly in the persistent input; the virtual K2R updates it.
         rows_bytes = n_pair_max*complex_bytes*REAL(bs_env%nkp_scf_desymm, dp) + &
                      real_bytes*(n_pair_max + 1.0_dp)*REAL(n_cells, dp)
         q_bytes = complex_bytes*REAL(n_pairs_total, dp)*MIN(REAL(q_batch, dp), n_q_max)
         ! The pair-local FFT output remains live while a bounded received batch is projected.
         green_bytes = complex_bytes*REAL(n_pairs_total, dp)* &
                       REAL(pair_exchange_width(n_grid, bs_env%nkp_scf_desymm, i_size), dp)
         live_bytes = MAX(rows_bytes + green_bytes + &
                          complex_bytes*REAL(n_threads, dp)*REAL(n_grid, dp)*REAL(bs_env%n_ao + 32, dp), &
                          real_bytes*(n_pair_max + 1.0_dp)*REAL(n_cells, dp) + &
                          complex_bytes*n_pair_max*REAL(SIZE(needed_q), dp) + q_bytes + &
                          complex_bytes*REAL(n_grid, dp)*REAL(n_aux + MIN(n_grid, grid_pair_panel_cols), dp)* &
                          MAX(1.0_dp, MIN(REAL(n_threads, dp), n_q_max, REAL(q_batch, dp)))) + &
                      REAL(int_1_size, dp)*n_pair_max*REAL(n_cells, dp) + &
                      MAX(REAL(fft_panel_bytes, dp), 3.0_dp*complex_bytes*REAL(n_cells, dp)) + &
                      real_bytes*(n_pair_max + REAL(n_cells, dp))
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
      CALL block_range(INT(INT(n_grid, int_8)*(n_grid + 1)/2), chi_env%mepos, chi_env%num_pe, first_pair, n_pairs_local)
      n_needed_local = SIZE(Z_time, 3)
      ALLOCATE (chi_partial(n_needed_local))
      DO iq_local = 1, n_needed_local
         ALLOCATE (chi_partial(iq_local)%values(n_aux, n_aux, n_time_local), SOURCE=z_zero)
      END DO
      ALLOCATE (chi_q(n_aux, n_aux, n_needed_local))

      ALLOCATE (green_norms(n_cells, 2), cell_active(n_cells))
      CALL prepare_pair_fft(bs_env, full_cells, first_pair, n_pairs_local, chi_env%num_pe, &
                            fft_stencil_work, .FALSE., needed_q)
      CALL start_time_loop_progress(bs_env, 'χ(iτ)', 'k points of G and q points of χ', n_groups, i_size, &
                                    INT(n_time_local*bs_env%n_spin, int_8)* &
                                    INT(2*bs_env%nkp_scf_desymm + COUNT(needed), int_8))
      DO i_t = i_group + 1, bs_env%num_time_freq_points, n_groups
         i_time_local = (i_t - 1)/n_groups + 1
         tau = bs_env%time_frequency_grid%imaginary_time(i_t)
         chi_q = z_zero
         DO ispin = 1, bs_env%n_spin
            chi_seconds = 0.0_dp
            IF (fft_profile_enabled) THEN
               chi_started = m_walltime()
               stage_started = chi_started
            END IF
            ALLOCATE (G_grid_k(n_pairs_local, 1, bs_env%nkp_scf_desymm))
            CALL build_green_pairs(bs_env, psi_values(:, :, :, ispin), chi_env, tau, ispin, &
                                   .TRUE., fft_stencil_work, G_grid_k)
            IF (fft_profile_enabled) chi_seconds(1) = m_walltime() - stage_started
            IF (fft_profile_enabled) stage_started = m_walltime()
            CALL k_to_grid_pairs(G_grid_k, fft_stencil_work, norm_squared=green_norms(:, 1))
            IF (fft_profile_enabled) chi_seconds(2) = m_walltime() - stage_started
            IF (fft_profile_enabled) stage_started = m_walltime()
            CALL build_green_pairs(bs_env, psi_values(:, :, :, ispin), chi_env, tau, ispin, &
                                   .FALSE., fft_stencil_work, G_grid_k)
            IF (fft_profile_enabled) chi_seconds(1) = chi_seconds(1) + m_walltime() - stage_started
            IF (fft_profile_enabled) stage_started = m_walltime()
            CALL k_to_grid_pairs(G_grid_k, fft_stencil_work, norm_squared=green_norms(:, 2), &
                                 product_scale=bs_env%spin_degeneracy)
            IF (fft_profile_enabled) chi_seconds(2) = chi_seconds(2) + m_walltime() - stage_started
            IF (fft_profile_enabled) stage_started = m_walltime()
            ! Each original Green is filtered by its own full Frobenius norm.
            ! Scalar cell masks commute with their product: (m_o G_o)(m_v G_v)=m_o m_v G_o G_v.
            CPASSERT(bs_env%n_atom > 0)
            CALL chi_env%sum(green_norms)
            cell_active(:) = .NOT. (SQRT(green_norms(:, 1))/REAL(bs_env%n_atom, dp)**2 < bs_env%eps_filter) .AND. &
                             .NOT. (SQRT(green_norms(:, 2))/REAL(bs_env%n_atom, dp)**2 < bs_env%eps_filter)
            IF (fft_profile_enabled) chi_seconds(3) = m_walltime() - stage_started
            IF (fft_profile_enabled) stage_started = m_walltime()
            DEALLOCATE (G_grid_k)
            ALLOCATE (chi_grid_q(n_pairs_local, 1, SIZE(needed_q)))
            CALL execute_pair_fft(fft_stencil_work, chi_grid_q, cell_active, &
                                  'GW chi (K2R extraction feeds R2K input)')
            IF (fft_profile_enabled) chi_seconds(4) = m_walltime() - stage_started

            IF (fft_profile_enabled) stage_started = m_walltime()
            CALL project_fft_pairs_to_chi(chi_grid_q, chi_env, needed, Z_time, &
                                          chi_q, q_batch, chi_seconds(6))
            DEALLOCATE (chi_grid_q)
            IF (fft_profile_enabled) chi_seconds(5) = MAX(0.0_dp, m_walltime() - stage_started - chi_seconds(6))
            IF (fft_profile_enabled) chi_seconds(7) = m_walltime() - chi_started
            unit_nr = 0
            IF (fft_profile_enabled) unit_nr = cp_logger_get_default_io_unit()
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A)') 'GW_CHI_PROFILE| logger-rank wall seconds, no MPI reduction'
               WRITE (unit_nr, '(T2,A,4I12)') 'GW_CHI_PROFILE| time/spin/pairs/cells: ', &
                  i_t, ispin, n_pairs_local, n_cells
               WRITE (unit_nr, '(T2,A,2I12)') 'GW_CHI_PROFILE| received q batch capacity / full local: ', &
                  MIN(q_batch, n_needed_local), n_needed_local
               ! product+R2K includes allocations, input release and nested FFT profile output.
               WRITE (unit_nr, '(T2,A,7F12.6)') &
                  'GW_CHI_PROFILE| build+alloc/K2R/filter/product+R2K/transfer/projection/total: ', chi_seconds
            END IF
            CALL count_time_loop_k_points(bs_env, INT(COUNT(needed), int_8))
         END DO
         ! The spin sum lives in the upper triangle. Complete each stored time slice once.
         DO i_needed = 1, n_needed_local
            DO i_aux = 1, n_aux
               chi_partial(i_needed)%values(:i_aux, i_aux, i_time_local) = chi_q(:i_aux, i_aux, i_needed)
               chi_partial(i_needed)%values(i_aux + 1:, i_aux, i_time_local) = CONJG(chi_q(i_aux, i_aux + 1:, i_needed))
            END DO
         END DO
      END DO
      CALL release_pair_fft(fft_stencil_work)
      bs_env%gwkp_env%k_points_total = 0

      IF (ALLOCATED(psi_time)) DEALLOCATE (psi_time)
      n_local = SIZE(gwkp_env%Z_lP_q, 3)
      ALLOCATE (gwkp_env%chi_q_frequency(n_local))
      iq_local = 0
      DO iq = MODULO(chi_env%mepos - 1, i_size) + 1, bs_env%nkp_scf_desymm, i_size
         IF (.NOT. needed(iq)) CYCLE
         iq_local = iq_local + 1
         i_root = MODULO(iq, n_world)/i_size
         i_q_world = gwkp_env%q_slot(iq)
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
      DEALLOCATE (chi_q, full_cells, needed, needed_q, green_norms, cell_active)

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
      REAL(KIND=dp), ALLOCATABLE                         :: local_q(:, :)
      TYPE(cp_fm_type), ALLOCATABLE                      :: fm_metric(:)

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI
      n_ref = n_RI
      IF (bs_env%auto_ri%enabled) n_ref = SUM(bs_env%auto_ri%sizes_ref_RI)
      nkp = bs_env%nkp_scf_desymm
      CPASSERT(bs_env%nkp_chi_eps_W_orig == nkp)
      nkp_local = SIZE(gwkp_env%Z_lP_q, 3)
      ALLOCATE (local_q(3, nkp_local))
      DO ikp = 1, nkp
         IF (MODULO(ikp, bs_env%para_env%num_pe) /= bs_env%para_env%mepos) CYCLE
         ilocal = gwkp_env%q_slot(ikp)
         IF (ilocal > 0) local_q(:, ilocal) = bs_env%kpoints_chi_eps_W%xkp(:, ikp)
      END DO
      ALLOCATE (M_ref(n_ref, n_ref, nkp_local), V_ref(n_ref, n_ref, nkp_local), source=z_zero)
      ALLOCATE (gwkp_env%metric_coulomb_factor(nkp_local), M_inv(n_RI, n_RI), V_sqrt(n_RI, n_RI))
      CALL build_ri_2c_lattice(fm_metric, metric_cells, bs_env%ri_metric, &
                               0.0_dp, bs_env, qs_env, reference_basis=bs_env%auto_ri%enabled)
      CALL ri_2c_lattice_to_owned_q(fm_metric, metric_cells, &
                                    bs_env%kpoints_chi_eps_W%xkp(1:3, 1:nkp), &
                                    bs_env%kpoints_chi_eps_W%nkp_grid, gwkp_env%q_slot, M_ref, &
                                    bs_env%kpoints_chi_eps_W%rskp_transform)
      ! The erfc(omega r)/r metric is the real-space Ewald part of the Coulomb matrix. At q = 0 of an
      ! odd mesh the reciprocal part drops the singular G = 0 term, which the Γ correction integrates.
      V_ref(:, :, :) = M_ref
      CALL add_ewald_reciprocal(V_ref, qs_env, local_q, 'RI_AUX', bs_env%ri_metric%omega)
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
         ALLOCATE (gwkp_env%metric_coulomb_factor(ilocal)%values(n_RI, n_RI, 1))
         gwkp_env%metric_coulomb_factor(ilocal)%values(:, :, 1) = MATMUL(M_inv, V_sqrt)
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

      COMPLEX(KIND=dp), ALLOCATABLE                      :: G_grid_k(:, :, :), head_cell_factor(:)
      INTEGER                                            :: first_pair, handle, ispin, n_grid, &
                                                            n_pairs_local
      INTEGER, ALLOCATABLE                               :: full_cells(:, :)
      LOGICAL, ALLOCATABLE                               :: cell_active(:)
      REAL(KIND=dp)                                      :: head_metric_zero
      REAL(KIND=dp), ALLOCATABLE                         :: Sigma_grid_R(:, :, :), V_grid_R(:, :, :)
      TYPE(gw_pair_fft_work_type)                        :: fft_stencil_work

      CALL timeset(routineN, handle)

      ALLOCATE (gwkp_env%sigma_x_bands(bs_env%n_ao, bs_env%nkp_bs_and_DOS, bs_env%n_spin))
      n_grid = bs_env%ri_rs%n_grid_points
      ! Same grid-pair distribution as Σc with q owners MOD(q,P).
      CALL block_range(INT(INT(n_grid, int_8)*(n_grid + 1)/2), bs_env%para_env%mepos, &
                       bs_env%para_env%num_pe, first_pair, n_pairs_local)
      IF (SIZE(gwkp_env%Z_lP_zero, 2) /= bs_env%n_RI + 1 .OR. SIZE(gwkp_env%Z_lP_zero, 3) /= 1) THEN
         CPABORT("Periodic RI-RS GW requires the direct Gamma auxiliary fit")
      END IF
      CPASSERT(SIZE(gwkp_env%Z_lP_q, 2) == bs_env%n_RI)
      CPASSERT(ALLOCATED(gwkp_env%metric_coulomb_factor))
      CALL build_full_cell_grid(bs_env, full_cells)
      CALL prepare_pair_fft(bs_env, full_cells, first_pair, n_pairs_local, bs_env%para_env%num_pe, &
                            fft_stencil_work, .FALSE., q_partner=gwkp_env%q_partner)
      ALLOCATE (G_grid_k(n_pairs_local, 1, bs_env%nkp_scf_desymm))
      CALL project_auxiliary_q_to_fft_pairs(gwkp_env%Z_lP_q, gwkp_env%metric_coulomb_factor, 1, &
                                            bs_env%para_env, fft_stencil_work, G_grid_k, factor=.TRUE.)
      DEALLOCATE (gwkp_env%metric_coulomb_factor)
      CALL q_mesh_and_gamma_shell_average(bs_env, values_q=G_grid_k, q_order=fft_stencil_work%input_to_native)
      ALLOCATE (V_grid_R(n_pairs_local, 1, SIZE(full_cells, 2)))
      CALL k_to_grid_pairs(G_grid_k, fft_stencil_work, values_R=V_grid_R)
      CALL singular_head_cell_factors(bs_env, full_cells, head_cell_factor)
      head_metric_zero = ri_metric_head_zero(bs_env)
      head_cell_factor = (MERGE(twopi, fourpi, ALL(bs_env%periodic == [1, 1, 0]))/head_metric_zero**2)*head_cell_factor
      CALL add_gamma_corr_pairs(bs_env, full_cells, head_cell_factor, &
                                gwkp_env%Z_lP_zero(:, 1, 1), first_pair, V_grid_R)
      DEALLOCATE (head_cell_factor)

      ALLOCATE (Sigma_grid_R, MOLD=V_grid_R)
      ALLOCATE (cell_active(SIZE(full_cells, 2)), SOURCE=.TRUE.)
      DO ispin = 1, bs_env%n_spin
         CALL build_green_pairs(bs_env, gwkp_env%band_orbitals(:, :, :, ispin), bs_env%para_env, &
                                0.0_dp, ispin, .TRUE., fft_stencil_work, G_grid_k)
         CALL k_to_grid_pairs(G_grid_k, fft_stencil_work, values_R=Sigma_grid_R, multiply_R=V_grid_R, product_scale=-1.0_dp)
         CALL execute_pair_fft(fft_stencil_work, G_grid_k, cell_active, &
                               'GW self-energy (K2R extraction feeds R2K input)', Sigma_grid_R)

         CALL self_energy_pairs_to_Bloch(bs_env, G_grid_k, Sigma_grid_R, bs_env%para_env, &
                                         fft_stencil_work, gwkp_env%sigma_projection, gwkp_env%band_orbitals, &
                                         ispin, gwkp_env%sigma_x_bands(:, :, ispin), .FALSE.)
      END DO
      CALL release_pair_fft(fft_stencil_work)

      DEALLOCATE (full_cells, G_grid_k, Sigma_grid_R, V_grid_R, cell_active)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x_from_V

! **************************************************************************************************
!> \brief Exact exchange energy with RI-RS from the real-space grid pairs, the self-energy product
!>        of compute_Sigma_x_from_V contracted with the occupied density matrix,
!>        E^EXX = -s/2 Σ_σ Σ_R Σ_ll' V_ll'(R) D^σ_ll'(R)², s the spin degeneracy,
!>        plus the q → 0 limit of the Coulomb singularity (exchange_energy_singular_correction),
!>        which is missing in V(q = 0) (no G = 0 term).
!> \param bs_env Band-structure settings, orbitals, and RI-RS work data.
!> \param gwkp_env RI-RS coefficients and the metric-Coulomb factor of compute_W0
!> \param e_exx exact exchange energy per unit cell
! **************************************************************************************************
   SUBROUTINE compute_exchange_energy_ri_rs(bs_env, gwkp_env, e_exx)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT), TARGET             :: gwkp_env
      REAL(KIND=dp), INTENT(OUT)                         :: e_exx

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_exchange_energy_ri_rs'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: D_k(:, :, :), V_q(:, :, :)
      INTEGER                                            :: first_in_batch, first_pair, handle, &
                                                            i_batch, i_cell, i_pair, il, ilp, &
                                                            ispin, last_in_batch, n_grid, &
                                                            n_pairs_local
      INTEGER, ALLOCATABLE                               :: full_cells(:, :)
      REAL(KIND=dp), ALLOCATABLE                         :: D_R(:, :, :), pair_weight(:), &
                                                            V_R(:, :, :)
      TYPE(gw_pair_fft_work_type)                        :: fft_stencil_work

      CALL timeset(routineN, handle)

      ! With 2 k-points along a_i the sum over R holds only R = 0 and the boundary cell R = -a_i.
      IF (ANY(MOD(bs_env%nkp_grid_chi_eps_W_orig, 2) == 0 .AND. bs_env%nkp_grid_chi_eps_W_orig < 4)) THEN
         CPABORT("RI_RS_RPA: E^EXX needs >= 4 k-points per direction on even meshes")
      END IF

      CALL collect_band_bloch_orbitals(bs_env, gwkp_env%band_orbitals)

      n_grid = bs_env%ri_rs%n_grid_points
      CALL block_range(INT(INT(n_grid, int_8)*(n_grid + 1)/2), bs_env%para_env%mepos, &
                       bs_env%para_env%num_pe, first_pair, n_pairs_local)
      CPASSERT(ALLOCATED(gwkp_env%metric_coulomb_factor))

      ! V_ll'(q) = (ZK)(ZK)† as in compute_Sigma_x_from_V
      CALL build_full_cell_grid(bs_env, full_cells)
      CALL prepare_pair_fft(bs_env, full_cells, first_pair, n_pairs_local, bs_env%para_env%num_pe, &
                            fft_stencil_work, .FALSE., q_partner=gwkp_env%q_partner)
      ALLOCATE (V_q(n_pairs_local, 1, bs_env%nkp_scf_desymm), D_k(n_pairs_local, 1, bs_env%nkp_scf_desymm))
      CALL project_auxiliary_q_to_fft_pairs(gwkp_env%Z_lP_q, gwkp_env%metric_coulomb_factor, 1, &
                                            bs_env%para_env, fft_stencil_work, V_q, factor=.TRUE.)
      DEALLOCATE (gwkp_env%metric_coulomb_factor)
      ALLOCATE (V_R(n_pairs_local, 1, SIZE(full_cells, 2)), D_R(n_pairs_local, 1, SIZE(full_cells, 2)))
      CALL k_to_grid_pairs(V_q, fft_stencil_work, values_R=V_R)
      DEALLOCATE (V_q)

      ! an unordered pair l < l' stands for (l,l') and (l',l)
      ALLOCATE (pair_weight(n_pairs_local))
      DO i_pair = 1, n_pairs_local
         CALL grid_pair_indices(first_pair + i_pair - 1, il, ilp)
         pair_weight(i_pair) = MERGE(1.0_dp, 2.0_dp, il == ilp)
      END DO

      e_exx = 0.0_dp
      DO ispin = 1, bs_env%n_spin
         CALL build_green_pairs(bs_env, gwkp_env%band_orbitals(:, :, :, ispin), bs_env%para_env, &
                                0.0_dp, ispin, .TRUE., fft_stencil_work, D_k)
         CALL k_to_grid_pairs(D_k, fft_stencil_work, values_R=D_R)
         ! The cells R are centred at the lattice vector nearest to r_l - r_l'. On a k mesh without
         ! the Gamma point (N k-points along a_i), V_ll'(R) changes sign from R to R + N a_i while
         ! D_ll'(R)² does not: the stencil weight carries this sign and removes the cells at
         ! distance N/2 from the centre, where the terms of R and -R cancel.
         DO i_batch = 1, SIZE(fft_stencil_work%batches)
            first_in_batch = (i_batch - 1)*fft_stencil_work%width + 1
            last_in_batch = first_in_batch + SIZE(fft_stencil_work%batches(i_batch)%weight, 1) - 1
            DO i_cell = 1, SIZE(full_cells, 2)
               e_exx = e_exx - 0.5_dp*bs_env%spin_degeneracy/REAL(grid_pair_stencil_denominator, dp)* &
                       SUM(pair_weight(first_in_batch:last_in_batch)* &
                           REAL(fft_stencil_work%batches(i_batch)%weight(:, i_cell), dp)* &
                           V_R(first_in_batch:last_in_batch, 1, i_cell)* &
                           D_R(first_in_batch:last_in_batch, 1, i_cell)**2)
            END DO
         END DO
      END DO
      CALL bs_env%para_env%sum(e_exx)
      CALL release_pair_fft(fft_stencil_work)

      DEALLOCATE (full_cells, D_k, D_R, V_R, pair_weight)

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(T2,A,T55,F26.10)') &
         'RI-RS-RPA| E^EXX without q -> 0 correction:', e_exx
      e_exx = e_exx + exchange_energy_singular_correction(bs_env)

      CALL timestop(handle)

   END SUBROUTINE compute_exchange_energy_ri_rs

! **************************************************************************************************
!> \brief q → 0 correction of the exchange energy for the Coulomb singularity with the periodic
!>        auxiliary function F(q) = Σ_G exp(-α|q+G|²)/|q+G|² of Gygi and Baldereschi,
!>        Phys. Rev. B 34, 4405 (1986). The exchange integrand is a/q² + O(q⁰) with the exact
!>        a = -2π N_e/Ω (N_e electrons per cell, cell volume Ω). The q-mesh sum misses
!>        a [⟨F⟩_BZ - 1/N_q Σ_{q≠0} F(q)] and, if q = 0 is on the mesh, a [α - Σ_{G≠0} F_G(0)]/N_q,
!>        the finite limit of a F(q) - a/q² at q = 0. ⟨F⟩_BZ = 2π^{3/2}/(√α V_BZ).
!>        Not included: the direction-dependent q⁰ term at q = 0 of Γ-centred (odd) meshes.
!> \param bs_env Band-structure settings with cell and q mesh.
!> \return correction to E^EXX per unit cell
! **************************************************************************************************
   FUNCTION exchange_energy_singular_correction(bs_env) RESULT(correction)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp)                                      :: correction

      INTEGER                                            :: i, iq, j, k, n_max(3), n_q
      LOGICAL                                            :: gamma_on_mesh
      REAL(KIND=dp)                                      :: a_singular, alpha, bmat(3, 3), &
                                                            cell_volume, mesh_sum, q(3), qg(3), qg2

      IF (ANY(bs_env%periodic /= 1)) THEN
         CPABORT("RI_RS_RPA exchange energy requires 3D periodicity")
      END IF

      n_q = PRODUCT(bs_env%nkp_grid_chi_eps_W_orig)
      cell_volume = ABS(det_3x3(bs_env%hmat))
      bmat = twopi*TRANSPOSE(inv_3x3(bs_env%hmat))
      a_singular = -twopi*bs_env%spin_degeneracy*SUM(bs_env%n_occ(1:bs_env%n_spin))/cell_volume
      alpha = 1.0_dp/MAXVAL(SUM(bmat**2, DIM=1))
      ! exp(-α|q+G|²) < e^-40 beyond n_max
      DO i = 1, 3
         n_max(i) = CEILING(SQRT(40.0_dp/alpha)*NORM2(bs_env%hmat(:, i))/twopi) + 1
      END DO

      mesh_sum = 0.0_dp
      gamma_on_mesh = .FALSE.
      DO iq = 1, n_q
         q = MATMUL(bmat, bs_env%kpoints_chi_eps_W%xkp(1:3, iq) - ANINT(bs_env%kpoints_chi_eps_W%xkp(1:3, iq)))
         IF (NORM2(q) < eps_q_gamma) THEN
            gamma_on_mesh = .TRUE.
            mesh_sum = mesh_sum - alpha
         END IF
         DO k = -n_max(3), n_max(3)
            DO j = -n_max(2), n_max(2)
               DO i = -n_max(1), n_max(1)
                  qg = q + i*bmat(:, 1) + j*bmat(:, 2) + k*bmat(:, 3)
                  qg2 = SUM(qg**2)
                  IF (qg2 < eps_q_gamma**2) CYCLE
                  mesh_sum = mesh_sum + EXP(-alpha*qg2)/qg2
               END DO
            END DO
         END DO
      END DO

      correction = a_singular*(2.0_dp*pi**1.5_dp/(SQRT(alpha)*twopi**3/cell_volume) - mesh_sum/REAL(n_q, dp))

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T55,F26.10)') 'RI-RS-RPA| E^EXX correction of the 1/q^2 singularity:', correction
         IF (gamma_on_mesh) WRITE (bs_env%unit_nr, '(T2,A)') &
            'RI-RS-RPA| Gamma-centred q mesh: the q^0 term at q = 0 is missing (error ~ 1/N_k)'
      END IF

   END FUNCTION exchange_energy_singular_correction

! **************************************************************************************************
!> \brief Form Σc with one owner per unordered (l,l') pair through both lattice FFTs.
!> \param bs_env Band-structure settings and physical results.
!> \param gwkp_env Periodic GW work data.
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c_from_W_PQ(bs_env, gwkp_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT), TARGET             :: gwkp_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Sigma_c_from_W_PQ'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: G_grid_k(:, :, :)
      INTEGER                                            :: first_pair, handle, i_component, i_t, &
                                                            i_time_local, ispin, n_grid, &
                                                            n_pairs_local, n_time_local
      INTEGER, ALLOCATABLE                               :: full_cells(:, :)
      LOGICAL, ALLOCATABLE                               :: cell_active(:)
      REAL(KIND=dp)                                      :: tau
      REAL(KIND=dp), ALLOCATABLE                         :: green_norms(:), sigma_band(:, :), &
                                                            Sigma_grid_R(:, :, :), &
                                                            sigma_result_buffer(:), &
                                                            W_grid_R(:, :, :)
      TYPE(gw_pair_fft_work_type)                        :: fft_stencil_work

      CALL timeset(routineN, handle)

      CPASSERT(ALLOCATED(gwkp_env%band_orbitals))
      IF (ALLOCATED(gwkp_env%sigma_c_bands)) DEALLOCATE (gwkp_env%sigma_c_bands)
      CALL build_full_cell_grid(bs_env, full_cells)
      n_grid = bs_env%ri_rs%n_grid_points
      CALL block_range(INT(INT(n_grid, int_8)*(n_grid + 1)/2), gwkp_env%group_env%mepos, gwkp_env%group_env%num_pe, &
                       first_pair, n_pairs_local)
      ALLOCATE (G_grid_k(n_pairs_local, 1, bs_env%nkp_scf_desymm))
      ALLOCATE (W_grid_R(n_pairs_local, 1, SIZE(full_cells, 2)))
      ALLOCATE (Sigma_grid_R, MOLD=W_grid_R)
      ALLOCATE (green_norms(SIZE(full_cells, 2)), cell_active(SIZE(full_cells, 2)))
      CPASSERT(bs_env%n_atom > 0)
      CALL prepare_pair_fft(bs_env, full_cells, first_pair, n_pairs_local, gwkp_env%group_env%num_pe, &
                            fft_stencil_work, .TRUE., q_partner=gwkp_env%q_partner)
      ALLOCATE (sigma_band(bs_env%n_ao, bs_env%nkp_bs_and_DOS))
      ALLOCATE (gwkp_env%sigma_c_bands(bs_env%n_ao, bs_env%nkp_bs_and_DOS, &
                                       bs_env%num_time_freq_points, 2, bs_env%n_spin))
      gwkp_env%sigma_c_bands = 0.0_dp

      n_time_local = (bs_env%num_time_freq_points - gwkp_env%my_group + gwkp_env%ngroups - 1)/gwkp_env%ngroups
      CALL start_time_loop_progress(bs_env, 'Σ_c(iτ)', 'k points of G and Σ_c', gwkp_env%ngroups, &
                                    gwkp_env%group_env%num_pe, INT(2*n_time_local*bs_env%n_spin, int_8)* &
                                    INT(bs_env%nkp_scf_desymm + bs_env%nkp_bs_and_DOS, int_8))
      DO i_t = gwkp_env%my_group + 1, bs_env%num_time_freq_points, gwkp_env%ngroups
         i_time_local = (i_t - gwkp_env%my_group - 1)/gwkp_env%ngroups + 1
         CALL project_auxiliary_q_to_fft_pairs(gwkp_env%Z_lP_q, &
                                               gwkp_env%W_PQ_tau, i_time_local, &
                                               gwkp_env%group_env, fft_stencil_work, G_grid_k, factor=.FALSE.)
         CALL q_mesh_and_gamma_shell_average(bs_env, values_q=G_grid_k, q_order=fft_stencil_work%input_to_native)
         CALL k_to_grid_pairs(G_grid_k, fft_stencil_work, values_R=W_grid_R)
         CALL add_gamma_corr_pairs(bs_env, full_cells, gwkp_env%gamma_correction%head_cell_tau(:, i_t), &
                                   gwkp_env%Z_lP_zero(:, 1, 1), first_pair, W_grid_R, gwkp_env%gamma_correction, i_t)
         tau = bs_env%time_frequency_grid%imaginary_time(i_t)
         DO ispin = 1, bs_env%n_spin
            ! Preserve the component order: occupied (2, -G W), then virtual (1, +G W).
            DO i_component = 2, 1, -1
               CALL build_green_pairs(bs_env, gwkp_env%band_orbitals(:, :, :, ispin), gwkp_env%group_env, &
                                      tau, ispin, i_component == 2, fft_stencil_work, G_grid_k)
               CALL k_to_grid_pairs(G_grid_k, fft_stencil_work, values_R=Sigma_grid_R, multiply_R=W_grid_R, &
                                    product_scale=MERGE(-1.0_dp, 1.0_dp, i_component == 2), norm_squared=green_norms)
               CALL gwkp_env%group_env%sum(green_norms)
               cell_active(:) = .NOT. (SQRT(green_norms)/REAL(bs_env%n_atom, dp)**2 < bs_env%eps_filter)
               CALL execute_pair_fft(fft_stencil_work, G_grid_k, cell_active, &
                                     'GW self-energy (K2R extraction feeds R2K input)', Sigma_grid_R)
               CALL self_energy_pairs_to_Bloch(bs_env, G_grid_k, Sigma_grid_R, gwkp_env%group_env, &
                                               fft_stencil_work, gwkp_env%sigma_projection, gwkp_env%band_orbitals, &
                                               ispin, sigma_band, .TRUE.)
               IF (gwkp_env%group_env%mepos == 0) gwkp_env%sigma_c_bands(:, :, i_t, i_component, ispin) = sigma_band
            END DO
         END DO
      END DO
      CALL release_pair_fft(fft_stencil_work)
      gwkp_env%k_points_total = 0

      ALLOCATE (sigma_result_buffer(SIZE(gwkp_env%sigma_c_bands)))
      sigma_result_buffer(:) = RESHAPE(gwkp_env%sigma_c_bands, [SIZE(gwkp_env%sigma_c_bands)])
      CALL bs_env%para_env%sum(sigma_result_buffer)
      gwkp_env%sigma_c_bands(:, :, :, :, :) = RESHAPE(sigma_result_buffer, SHAPE(gwkp_env%sigma_c_bands))

      DEALLOCATE (full_cells, G_grid_k, &
                  sigma_band, sigma_result_buffer, Sigma_grid_R, &
                  W_grid_R, green_norms, cell_active)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c_from_W_PQ

! **************************************************************************************************
!> \brief Print how the imaginary-time points of χ or Σc are split into groups of MPI ranks and
!> start the progress count of the group of rank 0.
!> Each group computes every n_groups-th time point without communicating with the other groups.
!> Only rank 0 writes output, so the progress lines are those of the group of rank 0.
!> \param bs_env Band-structure settings.
!> \param quantity χ(iτ) or Σ_c(iτ).
!> \param counted_in What the progress counts.
!> \param n_groups Number of groups that compute time points in parallel.
!> \param group_size MPI ranks per group.
!> \param k_points_total Number of k and q points the group of rank 0 processes in the time loop.
! **************************************************************************************************
   SUBROUTINE start_time_loop_progress(bs_env, quantity, counted_in, n_groups, group_size, k_points_total)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      CHARACTER(LEN=*), INTENT(IN)                       :: quantity, counted_in
      INTEGER, INTENT(IN)                                :: n_groups, group_size
      INTEGER(KIND=int_8), INTENT(IN)                    :: k_points_total

      INTEGER                                            :: n_time

      bs_env%gwkp_env%k_points_done = 0
      bs_env%gwkp_env%k_points_total = k_points_total
      bs_env%gwkp_env%next_progress = 1
      bs_env%gwkp_env%progress_start = m_walltime()
      bs_env%gwkp_env%progress_label = quantity//' progress of group'
      IF (bs_env%unit_nr <= 0) RETURN
      n_time = bs_env%num_time_freq_points
      IF (n_groups == 1) THEN
         WRITE (bs_env%unit_nr, '(/,T2,3A,I0,A,I0,A)') 'GW ', quantity, ': all ', n_time, &
            ' time points computed by one group of ', group_size, ' MPI ranks'
      ELSE
         WRITE (bs_env%unit_nr, '(/,T2,3A,I0,A,I0,A,I0,A)') 'GW ', quantity, ': ', n_time, &
            ' time points split over ', n_groups, ' groups of ', group_size, ' MPI ranks'
      END IF
      WRITE (bs_env%unit_nr, '(T2,6A)') 'GW ', quantity, ': progress of one group, counted in ', counted_in
      CALL m_flush(bs_env%unit_nr)

   END SUBROUTINE start_time_loop_progress

! **************************************************************************************************
!> \brief Count k points processed in the χ or Σc time loop and print the progress of the group of
!>        rank 0 at 0.1, 0.3, 1, 3, 10, 30, 60, 90 and 100 %. Does nothing outside these loops.
!> \param bs_env Band-structure settings.
!> \param n_k_points Number of k points just processed by the group.
! **************************************************************************************************
   SUBROUTINE count_time_loop_k_points(bs_env, n_k_points)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER(KIND=int_8), INTENT(IN)                    :: n_k_points

      TYPE(gwkp_type), POINTER                           :: gwkp_env

      gwkp_env => bs_env%gwkp_env
      IF (gwkp_env%k_points_total == 0) RETURN
      gwkp_env%k_points_done = gwkp_env%k_points_done + n_k_points
      CALL write_ri_rs_progress(bs_env%unit_nr, TRIM(gwkp_env%progress_label), gwkp_env%k_points_done, &
                                gwkp_env%k_points_total, gwkp_env%next_progress, gwkp_env%progress_start)

   END SUBROUTINE count_time_loop_k_points

! **************************************************************************************************
!> \brief Form the unordered pairs of G(k,τ)=U(k,τ) U(k,τ)^H and exchange bounded k batches.
!> Each upper-triangle Gram tile writes directly to its destination's message. No full
!> grid matrix, conjugate lower triangle, or separate send-packing pass is formed.
!> The exchange sends each reusable inversion pair once; pair owners restore its partner.
!> \param bs_env ...
!> \param orbitals Band Bloch orbitals of one spin at the k points owned in env, MOD(k-1,P).
!> \param env Ranks that share the grid pairs and the band k points.
!> \param tau ...
!> \param ispin ...
!> \param occupied ...
!> \param work shared geometry and exchange layout
!> \param G_pairs unordered pairs in source-rank-major k order within each exchange
! **************************************************************************************************
   SUBROUTINE build_green_pairs(bs_env, orbitals, env, tau, ispin, occupied, work, G_pairs)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), INTENT(IN)                       :: orbitals(:, :, :)
      TYPE(mp_para_env_type), POINTER                    :: env
      REAL(KIND=dp), INTENT(IN)                          :: tau
      INTEGER, INTENT(IN)                                :: ispin
      LOGICAL, INTENT(IN)                                :: occupied
      TYPE(gw_pair_fft_work_type), INTENT(IN)            :: work
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(OUT), TARGET  :: G_pairs(:, :, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_green_pairs'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: gram(:, :), scaled_coeff(:, :, :), &
                                                            sent(:)
      INTEGER :: chunk_first, first_band, handle, i_band, i_batch, i_col, i_first, i_job, ikp, &
         k_batch, n_active, n_bands, n_chunk, n_cols, n_grid, n_jobs, n_k, n_local, n_owned, &
         n_pairs, n_rank, n_round, n_scheduled, unit_nr
      INTEGER(KIND=int_8)                                :: active_bands, possible_bands
      INTEGER, ALLOCATABLE                               :: active_columns(:), k_jobs(:), &
                                                            recv_counts(:), recv_displs(:)
      REAL(KIND=dp)                                      :: seconds(4), started, tau_E, total_started
      TYPE(gwkp_type), POINTER                           :: gwkp_env

      CALL timeset(routineN, handle)
      IF (fft_profile_enabled) total_started = m_walltime()
      seconds = 0.0_dp
      active_bands = 0_int_8
      possible_bands = 0_int_8

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
      k_batch = work%k_batch
      n_pairs = SUM(work%pair_counts)
      CPASSERT(n_local == work%pair_counts(env%mepos))
      n_chunk = 1
!$    n_chunk = MIN(k_batch, omp_get_max_threads())
      ALLOCATE (scaled_coeff(n_grid, n_bands, n_chunk), active_columns(n_chunk))
      ALLOCATE (sent(MAX(1, n_pairs*k_batch)), k_jobs(k_batch), &
                recv_counts(0:n_rank - 1), recv_displs(0:n_rank - 1))
      IF (fft_profile_enabled) seconds(1) = m_walltime() - total_started
      ! G_ll′(k,τ) = Σ_n [ψ_ln(k) exp(-|τ(ε_nk-ε_F)|/2)]
      !                         [ψ_l′n(k) exp(-|τ(ε_nk-ε_F)|/2)]*.
      DO i_first = 1, n_k, n_rank*k_batch
         n_round = MIN(n_rank*k_batch, n_k - i_first + 1)
         CALL pair_exchange_jobs(work, work%green_source, env, i_first, n_round, 0, &
                                 k_jobs, n_owned, n_scheduled, recv_counts, recv_displs)
         IF (n_scheduled == 0) THEN
            CALL count_time_loop_k_points(bs_env, INT(n_round, int_8))
            CYCLE
         END IF
         IF (fft_profile_enabled) started = m_walltime()
         ! Scale at most one reciprocal point per worker, then distribute the Gram
         ! column blocks independently so one owned k point can occupy the whole team.
         !$OMP PARALLEL DEFAULT(NONE) IF(n_owned > 0) &
         !$OMP SHARED(bs_env, first_band, gwkp_env, k_jobs, n_bands, n_chunk, n_grid, n_owned, n_rank) &
         !$OMP SHARED(active_columns, orbitals, scaled_coeff, sent, tau, ispin, work) &
         !$OMP PRIVATE(chunk_first, gram, i_band, i_batch, i_col, i_job, ikp, n_active, n_cols, n_jobs, tau_E) &
         !$OMP REDUCTION(+:active_bands, possible_bands)
         ALLOCATE (gram(n_grid, MIN(n_grid, column_block_cols)))
         DO chunk_first = 1, n_owned, n_chunk
            n_jobs = MIN(n_chunk, n_owned - chunk_first + 1)
            !$OMP DO SCHEDULE(STATIC)
            DO i_job = 1, n_jobs
               ikp = k_jobs(chunk_first + i_job - 1)
               n_active = 0
               DO i_band = 1, n_bands
                  tau_E = ABS(0.5_dp*tau*(gwkp_env%eigenvalues_G(first_band + i_band - 1, ikp, ispin) - &
                                          bs_env%e_fermi(ispin)))
                  IF (tau_E < bs_env%stabilize_exp) THEN
                     n_active = n_active + 1
                     scaled_coeff(:, n_active, i_job) = &
                        orbitals(:, first_band + i_band - 1, (ikp - 1)/n_rank + 1)*EXP(-tau_E)
                  END IF
               END DO
               active_columns(i_job) = n_active
               IF (fft_profile_enabled) THEN
                  active_bands = active_bands + INT(n_active, int_8)
                  possible_bands = possible_bands + INT(n_bands, int_8)
               END IF
            END DO
            !$OMP END DO
            ! Successive column tiles grow in cost. Interleave them among workers.
            !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, 1)
            DO i_job = 1, n_jobs
               DO i_col = 1, n_grid, column_block_cols
                  i_batch = chunk_first + i_job - 1
                  n_cols = MIN(column_block_cols, n_grid - i_col + 1)
                  n_active = active_columns(i_job)
                  IF (n_active == 0) THEN
                     gram(:, :) = z_zero
                  ELSE
                     CALL zherk('U', 'N', n_cols, n_active, 1.0_dp, scaled_coeff(i_col, 1, i_job), n_grid, &
                                0.0_dp, gram(i_col, 1), n_grid)
                     IF (i_col > 1) THEN
                        CALL zgemm('N', 'C', i_col - 1, n_cols, n_active, z_one, scaled_coeff(:, :, i_job), n_grid, &
                                   scaled_coeff(i_col, 1, i_job), n_grid, z_zero, gram, n_grid)
                     END IF
                  END IF
                  CALL pack_grid_pair_panel(work, i_col, n_cols, gram, n_owned, i_batch, sent)
               END DO
            END DO
            !$OMP END DO
         END DO
         DEALLOCATE (gram)
         !$OMP END PARALLEL
         IF (fft_profile_enabled) seconds(2) = seconds(2) + m_walltime() - started
         IF (fft_profile_enabled) started = m_walltime()
         CALL exchange_pair_round(work, work%green_source, env, i_first, n_round, 0, &
                                  n_owned, sent, recv_counts, recv_displs, G_pairs)
         IF (fft_profile_enabled) seconds(3) = seconds(3) + m_walltime() - started
         CALL count_time_loop_k_points(bs_env, INT(n_round, int_8))
      END DO
      IF (work%green_reuses > 0) THEN
         IF (fft_profile_enabled) started = m_walltime()
         CALL restore_pair_conjugates(work%green_source, G_pairs)
         IF (fft_profile_enabled) seconds(3) = seconds(3) + m_walltime() - started
      END IF
      IF (fft_profile_enabled) seconds(4) = m_walltime() - total_started
      unit_nr = 0
      IF (fft_profile_enabled) unit_nr = cp_logger_get_default_io_unit()
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') 'GW_GREEN_PROFILE| logger-rank wall seconds, no MPI reduction'
         WRITE (unit_nr, '(T2,A,L2,2I12)') 'GW_GREEN_PROFILE| occupied / k points / local pairs: ', occupied, n_k, n_local
         WRITE (unit_nr, '(T2,A,2I16)') 'GW_GREEN_PROFILE| active / possible local band columns: ', active_bands, possible_bands
         WRITE (unit_nr, '(T2,A,I12)') 'GW_GREEN_PROFILE| k points per rank / exchange: ', k_batch
         WRITE (unit_nr, '(T2,A,I12)') 'GW_GREEN_PROFILE| conjugate k points reused: ', work%green_reuses
         WRITE (unit_nr, '(T2,A,4F12.6)') 'GW_GREEN_PROFILE| setup/build/transfer/total: ', seconds
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_green_pairs

! **************************************************************************************************
!> \brief Schedule only the representatives in one rank-blocked reciprocal batch.
!> \param work FFT input order and pair ownership
!> \param source representative slot for every FFT input column
!> \param env pair and reciprocal point owners
!> \param first first FFT input column of the batch
!> \param n_round number of FFT input columns of the batch
!> \param owner_shift shift of cyclic source ownership (zero for k, one for q)
!> \param jobs native reciprocal indices to compute on this rank
!> \param n_owned number of local jobs
!> \param n_scheduled number of jobs over all ranks
!> \param recv_counts ...
!> \param recv_displs ...
! **************************************************************************************************
   SUBROUTINE pair_exchange_jobs(work, source, env, first, n_round, owner_shift, &
                                 jobs, n_owned, n_scheduled, recv_counts, recv_displs)
      TYPE(gw_pair_fft_work_type), INTENT(IN)            :: work
      INTEGER, INTENT(IN)                                :: source(:)
      TYPE(mp_para_env_type), INTENT(IN)                 :: env
      INTEGER, INTENT(IN)                                :: first, n_round, owner_shift
      INTEGER, INTENT(OUT)                               :: jobs(:), n_owned, n_scheduled, &
                                                            recv_counts(0:), recv_displs(0:)

      INTEGER                                            :: block_first, block_size, i_input, &
                                                            i_rank, n_local, n_rank, packed, sender

      n_rank = env%num_pe
      n_local = work%pair_counts(env%mepos)
      n_owned = 0
      n_scheduled = 0
      block_first = first
      DO i_rank = 0, n_rank - 1
         sender = MODULO(i_rank + owner_shift, n_rank)
         block_size = (n_round + n_rank - 1 - i_rank)/n_rank
         packed = 0
         DO i_input = block_first, block_first + block_size - 1
            IF (source(i_input) /= i_input) CYCLE
            packed = packed + 1
            IF (sender == env%mepos) jobs(packed) = work%input_to_native(i_input)
         END DO
         recv_counts(sender) = n_local*packed
         recv_displs(sender) = n_local*(block_first - first)
         IF (sender == env%mepos) n_owned = packed
         n_scheduled = n_scheduled + packed
         block_first = block_first + block_size
      END DO
   END SUBROUTINE pair_exchange_jobs

! **************************************************************************************************
!> \brief Receive one batch directly into its FFT input and expand the packed source blocks in place.
!> \param work FFT input order and pair ownership
!> \param source representative slot for every FFT input column
!> \param env pair and reciprocal point owners
!> \param first first FFT input column of the batch
!> \param n_round number of FFT input columns of the batch
!> \param owner_shift shift of cyclic source ownership (zero for k, one for q)
!> \param n_owned number of local jobs
!> \param sent grid pairs packed by destination rank
!> \param recv_counts ...
!> \param recv_displs ...
!> \param pairs partially filled inverse FFT input
! **************************************************************************************************
   SUBROUTINE exchange_pair_round(work, source, env, first, n_round, owner_shift, &
                                  n_owned, sent, recv_counts, recv_displs, pairs)
      TYPE(gw_pair_fft_work_type), INTENT(IN)            :: work
      INTEGER, INTENT(IN)                                :: source(:)
      TYPE(mp_para_env_type), INTENT(IN)                 :: env
      INTEGER, INTENT(IN)                                :: first, n_round, owner_shift, n_owned
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(IN)           :: sent(:)
      INTEGER, CONTIGUOUS, INTENT(IN)                    :: recv_counts(0:), recv_displs(0:)
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(INOUT), &
         TARGET                                          :: pairs(:, :, :)

      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: received(:)
      COMPLEX(KIND=dp), TARGET                           :: unused_receive(1)
      INTEGER :: block_first, block_size, i_input, i_rank, n_local, n_rank, packed, sender, &
         send_displs(0:env%num_pe - 1), send_counts(0:env%num_pe - 1)

      n_rank = env%num_pe
      n_local = SIZE(pairs, 1)
      send_counts(:) = work%pair_counts*n_owned
      send_displs(:) = work%pair_displs*n_owned
      unused_receive = z_zero
      received => unused_receive
      IF (n_local > 0) received(1:n_local*n_round) => pairs(:, :, first:first + n_round - 1)
      CALL env%alltoall(sent, send_counts, send_displs, received, recv_counts, recv_displs)
      IF (n_local == 0) RETURN

      ! Expand from the last column to preserve data still to be read.
      DO i_rank = 0, n_rank - 1
         sender = MODULO(i_rank + owner_shift, n_rank)
         block_first = first + recv_displs(sender)/n_local
         block_size = (n_round + n_rank - 1 - i_rank)/n_rank
         packed = block_first + recv_counts(sender)/n_local - 1
         DO i_input = block_first + block_size - 1, block_first, -1
            IF (source(i_input) /= i_input) CYCLE
            IF (i_input /= packed) pairs(:, 1, i_input) = pairs(:, 1, packed)
            packed = packed - 1
         END DO
      END DO
   END SUBROUTINE exchange_pair_round

! **************************************************************************************************
!> \brief Complete the inverse FFT input after all representatives have been exchanged.
!> \param source representative slot for every FFT input column
!> \param pairs partially filled inverse FFT input
! **************************************************************************************************
   SUBROUTINE restore_pair_conjugates(source, pairs)
      INTEGER, INTENT(IN)                                :: source(:)
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(INOUT)        :: pairs(:, :, :)

      INTEGER                                            :: i_input

      IF (SIZE(pairs, 1) == 0) RETURN
      !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) IF(SIZE(pairs, 1) > 256) &
      !$OMP SHARED(pairs, source) PRIVATE(i_input)
      DO i_input = 1, SIZE(source)
         IF (source(i_input) /= i_input) pairs(:, 1, i_input) = CONJG(pairs(:, 1, source(i_input)))
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE restore_pair_conjugates

! **************************************************************************************************
!> \brief Packed unordered grid-pair blocks of all ranks, as assigned by block_range.
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
         CALL block_range(INT(INT(n_grid, int_8)*(n_grid + 1)/2), i_rank, n_rank, displs(i_rank), counts(i_rank))
      END DO
      displs(:) = displs - 1

   END SUBROUTINE grid_pair_blocks

! **************************************************************************************************
!> \brief Bound one exchange by its send buffer and the MPI wrapper's byte counter.
!> \param n_grid ...
!> \param n_k ...
!> \param n_rank ...
!> \return source k points per rank
! **************************************************************************************************
   PURE INTEGER FUNCTION pair_exchange_width(n_grid, n_k, n_rank) RESULT(width)
      INTEGER, INTENT(IN)                                :: n_grid, n_k, n_rank

      REAL(KIND=dp)                                      :: pair_elements, recv_elements

      pair_elements = REAL(n_grid, dp)*REAL(n_grid + 1, dp)/2.0_dp
      recv_elements = REAL(CEILING(pair_elements/REAL(n_rank, dp)), dp)*REAL(n_rank, dp)
      width = MAX(1, MIN((n_k + n_rank - 1)/n_rank, INT(projection_round_bytes/(complex_bytes*pair_elements)), &
                         INT(REAL(HUGE(0), dp)/(complex_bytes*(pair_elements + recv_elements)))))
   END FUNCTION pair_exchange_width

! **************************************************************************************************
!> \brief Write each upper-triangular column segment directly into its rank's message.
!> \param work ...
!> \param first_col ...
!> \param n_cols ...
!> \param panel ...
!> \param n_owned source points in this batch
!> \param i_batch point within the batch
!> \param sent destination-rank-major send buffer
! **************************************************************************************************
   SUBROUTINE pack_grid_pair_panel(work, first_col, n_cols, panel, n_owned, i_batch, sent)
      TYPE(gw_pair_fft_work_type), INTENT(IN)            :: work
      INTEGER, INTENT(IN)                                :: first_col, n_cols
      COMPLEX(KIND=dp), INTENT(IN)                       :: panel(:, :)
      INTEGER, INTENT(IN)                                :: n_owned, i_batch
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: sent(:)

      INTEGER                                            :: col, i_rank, j_col, offset, row_first, &
                                                            row_last, segment

      DO col = 1, n_cols
         j_col = first_col + col - 1
         DO segment = work%column_segments(j_col), work%column_segments(j_col + 1) - 1
            i_rank = work%segments(segment)%rank
            row_first = work%segments(segment)%first_row
            row_last = work%segments(segment)%last_row
            offset = work%pair_displs(i_rank)*n_owned + (i_batch - 1)*work%pair_counts(i_rank) + &
                     work%segments(segment)%displacement
            sent(offset + 1:offset + row_last - row_first + 1) = panel(row_first:row_last, col)
         END DO
      END DO
   END SUBROUTINE pack_grid_pair_panel

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
               CALL evaluate_mu_rl_R(bs_env, bs_env%ri_rs%grid_points, i_R, phi_R(:, :, i_R - first + 1))
            END DO
            CALL cell_to_k_transform(phi_R(:, :, 1:last - first + 1), &
                                     bs_env%index_to_cell_3c(:, first:last), &
                                     bs_env%kpoints_scf_desymm%xkp(:, owned_k), &
                                     bs_env%kpoints_scf_desymm%nkp_grid, phi_lmu_k, &
                                     bs_env%kpoints_scf_desymm%rskp_transform, accumulate=.TRUE., profile_fft=fft_profile_enabled)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_all_phi_k

! **************************************************************************************************
!> \brief Fourier transform distributed real RI lattice blocks to cyclic q owners in one exchange.
!>
!> I_PQ(q) = Σ_R exp(2πiq·R) I_PQ(R).  Each sender transforms its own matrix block
!> for the retained q points with the selected RSKP backend and sends the blocks to their q owners, iq mod P.
!> \param fm_R ...
!> \param cells ...
!> \param qpoints ...
!> \param nkp_grid regular q mesh of qpoints
!> \param q_slot compact owner-local indices; zero omits a q point
!> \param matrices_q ...
!> \param backend GEMM or FFT
! **************************************************************************************************
   SUBROUTINE ri_2c_lattice_to_owned_q(fm_R, cells, qpoints, nkp_grid, q_slot, matrices_q, backend)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_R(:)
      INTEGER, INTENT(IN)                                :: cells(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: qpoints(:, :)
      INTEGER, INTENT(IN)                                :: nkp_grid(3), q_slot(:)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: matrices_q(:, :, :)
      INTEGER, INTENT(IN)                                :: backend

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ri_2c_lattice_to_owned_q'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: recv_buffer(:)
      COMPLEX(KIND=dp), ALLOCATABLE, TARGET              :: send_buffer(:)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: values_q(:, :, :)
      INTEGER                                            :: dest, handle, icell, icol, iq, iq_local, &
                                                            irow, n_q, n_q_local, n_RI, &
                                                            n_selected, ncols, nprocs, nrecv, &
                                                            nrows, nsend, pos, source
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
      n_q_local = COUNT(q_slot(MODULO(para_env%mepos - 1, nprocs) + 1:n_q:nprocs) > 0)
      n_selected = COUNT(q_slot > 0)
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
      ALLOCATE (q_order(n_selected), send_counts(0:nprocs - 1), recv_counts(0:nprocs - 1), &
                send_displs(0:nprocs - 1), recv_displs(0:nprocs - 1))
      pos = 0
      DO dest = 0, nprocs - 1
         send_displs(dest) = pos*nrows*ncols
         DO iq = 1, n_q
            IF (q_slot(iq) == 0 .OR. MODULO(iq, nprocs) /= dest) CYCLE
            pos = pos + 1
            q_order(pos) = iq
         END DO
         send_counts(dest) = pos*nrows*ncols - send_displs(dest)
      END DO
      nsend = n_selected*nrows*ncols
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
         values_q(1:nrows, 1:ncols, 1:n_selected) => send_buffer(1:nsend)
         CALL cell_to_k_transform(values_R, cells, qpoints(:, q_order), nkp_grid, values_q, backend, &
                                  profile_fft=fft_profile_enabled)
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
                                                            n_grid, n_own, n_partner
      REAL(dp)                                           :: available_GB, orbital_bytes
      REAL(dp), ALLOCATABLE                              :: phi_R(:, :)
      TYPE(mp_para_env_type), POINTER                    :: env

      CALL timeset(routineN, handle)

      env => bs_env%para_env
      n_grid = bs_env%ri_rs%n_grid_points
      CALL owned_mo_coefficients(bs_env, bs_env%cfm_mo_coeff_kp, coefficients)
      n_own = SIZE(coefficients, 3)
      n_partner = 0
      DO i_kp = env%mepos + 1, bs_env%nkp_scf_desymm, env%num_pe
         IF (bs_env%gwkp_env%q_partner(i_kp) > i_kp) n_partner = n_partner + 1
      END DO
      orbital_bytes = REAL(n_grid, dp)*REAL(bs_env%n_ao, dp)* &
                      (complex_bytes*REAL((n_own + n_partner)*bs_env%n_spin, dp) + complex_bytes + real_bytes)
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
         CALL evaluate_mu_rl_R(bs_env, bs_env%ri_rs%grid_points, i_R, phi_R)
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
!> \brief Schedule one Σ matrix per native inversion pair and cache the actual partner orbitals.
!> \param bs_env ...
!> \param env ...
!> \param orbitals Orbitals with cyclic ownership, MOD(k-1,P).
!> \param plan Destination-rank job lists and conjugated partner orbitals.
! **************************************************************************************************
   SUBROUTINE prepare_sigma_projection(bs_env, env, orbitals, plan)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), INTENT(IN)                 :: env
      COMPLEX(KIND=dp), ASYNCHRONOUS, CONTIGUOUS, &
         INTENT(IN)                                      :: orbitals(:, :, :, :)
      TYPE(sigma_projection_type), ASYNCHRONOUS, &
         INTENT(OUT)                                     :: plan

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_sigma_projection'

      COMPLEX(KIND=dp)                                   :: phase
      INTEGER :: c1, c2, c3, first_pair, handle, i_dim, i_image, i_kp, i_local, i_pair, i_partner, &
         i_rank, i_slot, il, ilp, ispin, n_cells, n_grid, n_kp, n_mesh, n_off, n_pairs, n_partner, &
         n_rank, n_requests, nmesh(3), partner_rank, tied
      LOGICAL                                            :: tie(3)
      REAL(KIND=dp)                                      :: angle, boundary(3), xk(3)
      REAL(KIND=dp), ALLOCATABLE                         :: fractional_grid(:, :)
      TYPE(mp_request_type), ALLOCATABLE                 :: requests(:)

      CALL timeset(routineN, handle)
      n_kp = bs_env%nkp_bs_and_DOS
      n_mesh = bs_env%nkp_scf_desymm
      n_rank = env%num_pe
      ALLOCATE (plan%counts(0:n_rank - 1), plan%mesh_index(n_kp), plan%partner(n_kp), plan%partner_slot(n_kp))
      ALLOCATE (plan%jobs((n_kp + n_rank - 1)/n_rank, 0:n_rank - 1))
      plan%counts(:) = 0
      plan%partner(:) = 0
      plan%partner_slot(:) = 0
      DO i_kp = 1, n_kp
         plan%mesh_index(i_kp) = native_q_index(bs_env, i_kp)
      END DO
      n_partner = 0
      DO i_kp = 1, n_kp
         IF (i_kp <= n_mesh) THEN
            i_partner = bs_env%gwkp_env%q_partner(i_kp)
            ! Σ comes from the real-output K2R and a real-space product. Its R2K
            ! conjugacy also holds when evGW0 constructs both Green spectra separately.
            ! Keep independent band-path targets even when they land on the mesh.
            IF (i_kp /= i_partner .AND. plan%mesh_index(i_kp) == i_kp .AND. &
                plan%mesh_index(i_partner) == i_partner) THEN
               IF (i_partner < i_kp) CYCLE
               plan%partner(i_kp) = i_partner
            END IF
         END IF
         i_rank = MODULO(i_kp - 1, n_rank)
         plan%counts(i_rank) = plan%counts(i_rank) + 1
         plan%jobs(plan%counts(i_rank), i_rank) = i_kp
         IF (i_rank /= env%mepos .OR. plan%partner(i_kp) == 0) CYCLE
         n_partner = n_partner + 1
         plan%partner_slot(i_kp) = n_partner
      END DO
      ALLOCATE (plan%partner_orbitals(SIZE(orbitals, 1), SIZE(orbitals, 2), n_partner, SIZE(orbitals, 4)))
      ALLOCATE (requests((n_partner + SIZE(orbitals, 3))*SIZE(orbitals, 4)))
      n_requests = 0
      ! MPI preserves message order between each pair of ranks. Both peers visit
      ! the representatives and spin channels in the same order.
      DO i_kp = 1, n_mesh
         i_partner = plan%partner(i_kp)
         IF (i_partner == 0) CYCLE
         i_rank = MODULO(i_kp - 1, n_rank)
         partner_rank = MODULO(i_partner - 1, n_rank)
         i_local = (i_partner - 1)/n_rank + 1
         i_slot = plan%partner_slot(i_kp)
         DO ispin = 1, SIZE(orbitals, 4)
            IF (env%mepos == i_rank) THEN
               IF (i_rank == partner_rank) THEN
                  plan%partner_orbitals(:, :, i_slot, ispin) = orbitals(:, :, i_local, ispin)
               ELSE
                  n_requests = n_requests + 1
                  CALL env%irecv(plan%partner_orbitals(:, :, i_slot, ispin), partner_rank, requests(n_requests), tag=0)
               END IF
            ELSE IF (env%mepos == partner_rank) THEN
               n_requests = n_requests + 1
               CALL env%isend(orbitals(:, :, i_local, ispin), i_rank, requests(n_requests), tag=0)
            END IF
         END DO
      END DO
      CALL mp_waitall(requests(:n_requests))
      !$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(2) SCHEDULE(STATIC) IF(n_partner > 1) &
      !$OMP SHARED(plan, n_partner, orbitals) PRIVATE(ispin, i_slot)
      DO ispin = 1, SIZE(orbitals, 4)
         DO i_slot = 1, n_partner
            plan%partner_orbitals(:, :, i_slot, ispin) = CONJG(plan%partner_orbitals(:, :, i_slot, ispin))
         END DO
      END DO
      !$OMP END PARALLEL DO
      ! The image geometry and off-mesh phases are common to Σx and every Σc time slice.
      n_off = COUNT(plan%mesh_index == 0)
      n_grid = SIZE(orbitals, 1)
      CALL block_range(INT(INT(n_grid, int_8)*(n_grid + 1)/2), env%mepos, n_rank, first_pair, n_pairs)
      IF (n_off > 0 .AND. n_pairs > 0) THEN
         ALLOCATE (plan%image_shift(3, n_pairs), plan%image_tie(3, n_pairs), plan%displacement(3, n_pairs))
         ALLOCATE (fractional_grid(3, n_grid))
         fractional_grid(:, :) = MATMUL(bs_env%ri_rs%cell%h_inv, bs_env%ri_rs%grid_points)
         DO i_dim = 1, 3
            IF (bs_env%periodic(i_dim) == 0) fractional_grid(i_dim, :) = 0.0_dp
         END DO
         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
         !$OMP SHARED(bs_env, first_pair, fractional_grid, n_pairs, plan) PRIVATE(i_pair, il, ilp, tie)
         DO i_pair = 1, n_pairs
            CALL grid_pair_indices(first_pair + i_pair - 1, il, ilp)
            CALL grid_pair_shift(il, ilp, fractional_grid, plan%image_shift(:, i_pair), tie)
            plan%image_tie(:, i_pair) = MERGE(1_int_1, 0_int_1, tie)
            plan%displacement(:, i_pair) = bs_env%ri_rs%grid_points(:, il) - bs_env%ri_rs%grid_points(:, ilp)
         END DO
         !$OMP END PARALLEL DO
         DEALLOCATE (fractional_grid)

         nmesh = bs_env%kpoints_scf_desymm%nkp_grid
         boundary = 0.25_dp*(1.0_dp + COS(twopi*bs_env%kpoints_scf_desymm%xkp(:, 1)*nmesh))
         ! Omit endpoints with zero weight for both pair-center choices from the GEMM dimension.
         plan%first_image = -nmesh/2
         DO i_dim = 1, 3
            IF (image_weight(plan%first_image(i_dim), nmesh(i_dim), .FALSE., boundary(i_dim)) == 0.0_dp) THEN
               plan%first_image(i_dim) = plan%first_image(i_dim) + 1
            END IF
         END DO
         plan%last_image = -plan%first_image + MERGE(1, 0, ANY(plan%image_tie /= 0_int_1, DIM=2))
         n_cells = PRODUCT(plan%last_image - plan%first_image + 1)
         ALLOCATE (plan%image_weights(MAXVAL(plan%last_image - plan%first_image + 1), 3, 0:1))
         DO tied = 0, 1
            DO i_dim = 1, 3
               DO i_image = plan%first_image(i_dim), plan%last_image(i_dim)
                  plan%image_weights(i_image - plan%first_image(i_dim) + 1, i_dim, tied) = &
                     image_weight(i_image, nmesh(i_dim), tied == 1, boundary(i_dim))
               END DO
            END DO
         END DO
         ALLOCATE (plan%offmesh_slot(n_kp), plan%phase_cos(n_cells, n_off), plan%phase_sin(n_cells, n_off))
         plan%offmesh_slot(:) = 0
         n_off = 0
         DO i_kp = 1, n_kp
            IF (plan%mesh_index(i_kp) /= 0) CYCLE
            n_off = n_off + 1
            plan%offmesh_slot(i_kp) = n_off
         END DO
         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
         !$OMP SHARED(bs_env, n_kp, plan) PRIVATE(angle, c1, c2, c3, i_image, i_kp, i_slot, phase, xk)
         DO i_kp = 1, n_kp
            i_slot = plan%offmesh_slot(i_kp)
            IF (i_slot == 0) CYCLE
            xk = bs_env%kpoints_DOS%xkp(:, i_kp)
            i_image = 0
            DO c3 = plan%first_image(3), plan%last_image(3)
               DO c2 = plan%first_image(2), plan%last_image(2)
                  DO c1 = plan%first_image(1), plan%last_image(1)
                     i_image = i_image + 1
                     angle = twopi*(xk(1)*c1 + xk(2)*c2 + xk(3)*c3)
                     phase = EXP(CMPLX(0.0_dp, angle, KIND=dp))
                     plan%phase_cos(i_image, i_slot) = REAL(phase, KIND=dp)
                     plan%phase_sin(i_image, i_slot) = AIMAG(phase)
                  END DO
               END DO
            END DO
         END DO
         !$OMP END PARALLEL DO
      END IF
      CALL timestop(handle)
   END SUBROUTINE prepare_sigma_projection

! **************************************************************************************************
!> \brief Exchange batched Σ matrices and contract the packed panels with the actual target orbitals.
!> \param bs_env ...
!> \param sigma_mesh Native Σ(k) for the local pairs.
!> \param sigma_R Real-space Σ for the local pairs.
!> \param env Pair and band owners.
!> \param fft FFT pair partition.
!> \param plan Projection jobs and conjugated partner orbitals.
!> \param orbitals Local band orbitals.
!> \param ispin Spin channel.
!> \param sigma_band Diagonals in the original band-point order.
!> \param apply_cutoff Use the W distance cutoff in off-mesh image sums.
! **************************************************************************************************
   SUBROUTINE self_energy_pairs_to_Bloch(bs_env, sigma_mesh, sigma_R, env, fft, plan, orbitals, &
                                         ispin, sigma_band, apply_cutoff)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), INTENT(IN)                       :: sigma_mesh(:, :, :)
      REAL(KIND=dp), INTENT(IN)                          :: sigma_R(:, :, :)
      TYPE(mp_para_env_type), INTENT(IN)                 :: env
      TYPE(gw_pair_fft_work_type), INTENT(IN)            :: fft
      TYPE(sigma_projection_type), INTENT(IN), TARGET    :: plan
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(IN), TARGET   :: orbitals(:, :, :, :)
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp), INTENT(OUT)                         :: sigma_band(:, :)
      LOGICAL, INTENT(IN)                                :: apply_cutoff

      CHARACTER(LEN=*), PARAMETER :: routineN = 'self_energy_pairs_to_Bloch'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: panel(:, :), received(:), work(:, :)
      COMPLEX(KIND=dp), ALLOCATABLE, TARGET              :: send_buffer(:, :)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: psi(:, :), sent(:)
      COMPLEX(KIND=dp), TARGET                           :: unused_send(1)
      INTEGER :: first_col, handle, i_band, i_batch, i_col, i_first, i_job, i_kp, i_local, &
         i_orbital, i_panel, i_rank, i_slot, i_target, last_col, m, n_ao, n_batch, n_cols, n_done, &
         n_grid, n_kp, n_local, n_recv, n_rhs, n_round, n_send, offset, row_first, row_last, &
         segment
      INTEGER, ALLOCATABLE                               :: jobs(:), recv_counts(:), recv_displs(:), &
                                                            send_counts(:), send_displs(:)

      CALL timeset(routineN, handle)

      n_grid = bs_env%ri_rs%n_grid_points
      n_ao = bs_env%n_ao
      n_kp = bs_env%nkp_bs_and_DOS
      n_local = SIZE(sigma_mesh, 1)
      n_batch = MIN(sigma_rhs_cols, n_ao)
      m = projection_round_kpoints(n_grid, n_kp, env%num_pe)
      n_round = SUM(MIN(m, plan%counts))
      n_recv = MIN(m, plan%counts(env%mepos))
      ALLOCATE (received(MAX(1, SUM(fft%pair_counts)*n_recv)), &
                send_buffer(n_local, n_round), jobs(n_round))
      ALLOCATE (recv_counts(0:env%num_pe - 1), recv_displs(0:env%num_pe - 1), &
                send_counts(0:env%num_pe - 1), send_displs(0:env%num_pe - 1))
      unused_send = z_zero
      sigma_band(:, :) = 0.0_dp

      DO i_first = 1, MAXVAL(plan%counts), m
         ! The send buffer is destination-major, with consecutive k columns per destination.
         n_round = 0
         DO i_rank = 0, env%num_pe - 1
            n_send = MAX(0, MIN(m, plan%counts(i_rank) - i_first + 1))
            send_displs(i_rank) = n_local*n_round
            send_counts(i_rank) = n_local*n_send
            IF (n_send > 0) THEN
               jobs(n_round + 1:n_round + n_send) = plan%jobs(i_first:i_first + n_send - 1, i_rank)
            END IF
            n_round = n_round + n_send
         END DO
         sent => unused_send
         IF (n_local > 0) THEN
            CALL band_self_energy(bs_env, sigma_mesh, sigma_R, jobs(:n_round), plan, &
                                  send_buffer(:, :n_round), apply_cutoff)
            sent(1:n_local*n_round) => send_buffer(:, :n_round)
         END IF
         n_recv = MAX(0, MIN(m, plan%counts(env%mepos) - i_first + 1))
         recv_counts(:) = fft%pair_counts*n_recv
         recv_displs(:) = fft%pair_displs*n_recv
         CALL env%alltoall(sent, send_counts, send_displs, received, recv_counts, recv_displs)

         IF (n_recv > 0) THEN
            BLOCK
               REAL(KIND=dp) :: diagonal(n_ao, 2, n_recv)

               diagonal(:, :, :) = 0.0_dp
               !$OMP PARALLEL DEFAULT(NONE) REDUCTION(+:diagonal) &
               !$OMP SHARED(env, fft, i_first, ispin, n_ao, n_batch, n_grid, n_recv, orbitals, plan, received, recv_displs) &
               !$OMP PRIVATE(first_col, i_band, i_batch, i_col, i_job, i_kp, i_local, i_orbital, i_panel, i_rank, i_slot, &
               !$OMP         i_target, last_col, n_cols, n_rhs, offset, panel, psi, row_first, row_last, segment, work)
               ALLOCATE (panel(n_grid, MIN(n_grid, grid_pair_panel_cols)), &
                         work(n_batch, MIN(n_grid, grid_pair_panel_cols)))
               !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, 1)
               DO i_job = 1, n_recv
                  DO first_col = 1, n_grid, grid_pair_panel_cols
                     last_col = MIN(first_col + grid_pair_panel_cols - 1, n_grid)
                     n_cols = last_col - first_col + 1
                     ! Each column may span multiple MPI sources. Read their intersections
                     ! directly from the receive buffer and form U with half the diagonal.
                     DO i_col = first_col, last_col
                        i_panel = i_col - first_col + 1
                        DO segment = fft%column_segments(i_col), fft%column_segments(i_col + 1) - 1
                           i_rank = fft%segments(segment)%rank
                           row_first = fft%segments(segment)%first_row
                           row_last = fft%segments(segment)%last_row
                           offset = recv_displs(i_rank) + (i_job - 1)*fft%pair_counts(i_rank) + &
                                    fft%segments(segment)%displacement
                           panel(row_first:row_last, i_panel) = received(offset + 1:offset + row_last - row_first + 1)
                        END DO
                        panel(i_col, i_panel) = CMPLX(0.5_dp*REAL(panel(i_col, i_panel), dp), 0.0_dp, dp)
                        panel(i_col + 1:last_col, i_panel) = z_zero
                     END DO
                     i_kp = plan%jobs(i_first + i_job - 1, env%mepos)
                     i_local = (i_kp - 1)/env%num_pe + 1
                     i_slot = plan%partner_slot(i_kp)
                     ! Contract Psi^H U_panel first so the result has only one column per panel column.
                     DO i_target = 1, MERGE(2, 1, i_slot > 0)
                        IF (i_target == 1) THEN
                           psi => orbitals(:, :, i_local, ispin)
                        ELSE
                           psi => plan%partner_orbitals(:, :, i_slot, ispin)
                        END IF
                        DO i_batch = 1, n_ao, n_batch
                           n_rhs = MIN(n_batch, n_ao - i_batch + 1)
                           CALL zgemm('C', 'N', n_rhs, n_cols, last_col, z_one, psi(:, i_batch:i_batch + n_rhs - 1), &
                                      n_grid, panel, n_grid, z_zero, work, n_batch)
                           DO i_band = 1, n_rhs
                              i_orbital = i_batch + i_band - 1
                              diagonal(i_orbital, i_target, i_job) = diagonal(i_orbital, i_target, i_job) + &
                                                             REAL(SUM(work(i_band, :n_cols)*psi(first_col:last_col, i_orbital)), dp)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
               !$OMP END DO NOWAIT
               DEALLOCATE (panel, work)
               !$OMP END PARALLEL
               DO i_job = 1, n_recv
                  i_kp = plan%jobs(i_first + i_job - 1, env%mepos)
                  sigma_band(:, i_kp) = 2.0_dp*diagonal(:, 1, i_job)
                  IF (plan%partner(i_kp) > 0) sigma_band(:, plan%partner(i_kp)) = 2.0_dp*diagonal(:, 2, i_job)
               END DO
            END BLOCK
         END IF
         n_done = n_round + COUNT(plan%partner(jobs(:n_round)) > 0)
         CALL count_time_loop_k_points(bs_env, INT(n_done, int_8))
      END DO
      CALL env%sum(sigma_band)

      CALL timestop(handle)

   END SUBROUTINE self_energy_pairs_to_Bloch

! **************************************************************************************************
!> \brief Bound the peak per-rank workspace for the Σ(k) exchange and projection.
!> \param n_grid ...
!> \param n_rank ...
!> \param n_kp Total band targets.
!> \param n_jobs Targets exchanged after pairing conjugate matrices.
!> \param n_offmesh Targets requiring real-space interpolation.
!> \param n_images Bound of the image cells summed per pair.
!> \param n_ao ...
!> \param n_threads ...
!> \return ...
! **************************************************************************************************
   FUNCTION projection_workspace_bytes(n_grid, n_rank, n_kp, n_jobs, n_offmesh, n_images, n_ao, n_threads) RESULT(bytes)
      INTEGER, INTENT(IN)                                :: n_grid, n_rank, n_kp, n_jobs, n_offmesh, &
                                                            n_images, n_ao, n_threads
      REAL(dp)                                           :: bytes

      REAL(dp), PARAMETER :: integer_bytes = REAL(STORAGE_SIZE(0)/8, dp)

      REAL(dp)                                           :: contraction_bytes, exchange_bytes, &
                                                            interpolation_bytes, m, n_c, n_off, &
                                                            n_pairs, n_round

      m = REAL(MIN(n_jobs, projection_round_kpoints(n_grid, n_kp, n_rank)), dp)
      n_round = MIN(m*REAL(n_rank, dp), REAL(n_jobs, dp))
      n_pairs = REAL(n_grid, dp)*REAL(n_grid + 1, dp)/2.0_dp
      exchange_bytes = complex_bytes*(REAL(CEILING(n_pairs/REAL(n_rank, dp), KIND=int_8), dp)*n_round + n_pairs*m + 2.0_dp) + &
                       integer_bytes*(n_round + 4.0_dp*REAL(n_rank, dp))
      n_off = MIN(n_round, REAL(n_offmesh, dp))
      interpolation_bytes = integer_bytes*n_round
      IF (n_off > 0.0_dp) THEN
         n_c = REAL(n_images, dp)
         interpolation_bytes = interpolation_bytes + real_bytes*((2.0_dp*n_c + 3.0_dp)*n_off + &
                                                                 REAL(n_threads*pair_batch_size, dp)*(n_c + 2.0_dp*n_off))
      END IF
      contraction_bytes = REAL(n_threads, dp)*complex_bytes*REAL(MIN(n_grid, grid_pair_panel_cols), dp)* &
                          REAL(n_grid + MIN(sigma_rhs_cols, n_ao), dp) + &
                          2.0_dp*real_bytes*REAL(n_ao, dp)*m*REAL(n_threads + 1, dp)
      ! The interpolation temporaries are released before the panel contraction allocates its work arrays.
      bytes = exchange_bytes + MAX(interpolation_bytes, contraction_bytes)
   END FUNCTION projection_workspace_bytes

! **************************************************************************************************
!> \brief Band k points per rank and projection round: a round evaluates up to
!> band_kpoints_per_pass k points, with the exchange width bounded as for Green and W.
!> \param n_grid ...
!> \param n_kp ...
!> \param n_rank ...
!> \return ...
! **************************************************************************************************
   PURE INTEGER FUNCTION projection_round_kpoints(n_grid, n_kp, n_rank) RESULT(m)
      INTEGER, INTENT(IN)                                :: n_grid, n_kp, n_rank

      m = MAX(1, MIN(band_kpoints_per_pass/n_rank, pair_exchange_width(n_grid, n_kp, n_rank)))
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
!> \param k_targets band k points in send-buffer order
!> \param plan cached native indices, image geometry, and off-mesh phases
!> \param sigma_k ...
!> \param apply_cutoff use the W distance cutoff on the same centered images as the native FFT
! **************************************************************************************************
   SUBROUTINE band_self_energy(bs_env, sigma_mesh, sigma_R, k_targets, plan, sigma_k, apply_cutoff)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), INTENT(IN)                       :: sigma_mesh(:, :, :)
      REAL(KIND=dp), INTENT(IN)                          :: sigma_R(:, :, :)
      INTEGER, INTENT(IN)                                :: k_targets(:)
      TYPE(sigma_projection_type), INTENT(IN)            :: plan
      COMPLEX(KIND=dp), INTENT(OUT)                      :: sigma_k(:, :)
      LOGICAL, INTENT(IN)                                :: apply_cutoff

      INTEGER :: c1, c2, c3, first_c(3), i_c, i_col, i_first, i_local, i_off, i_pair, i_slot, &
         last_c(3), lower(3), n_c, n_off, n_pairs, native_iq, nmesh(3), shift(3), tied(3), &
         wrapped(3)
      INTEGER, ALLOCATABLE                               :: offmesh_col(:)
      REAL(KIND=dp)                                      :: angle, cutoff_squared, distance(3), &
                                                            weight
      REAL(KIND=dp), ALLOCATABLE                         :: cos_c(:, :), images(:, :), sin_c(:, :), &
                                                            sum_im(:, :), sum_re(:, :), xk(:, :)

      ALLOCATE (offmesh_col(SIZE(sigma_k, 2)))
      n_off = 0
      !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
      !$OMP SHARED(k_targets, plan, sigma_k, sigma_mesh) PRIVATE(i_col, native_iq)
      DO i_col = 1, SIZE(sigma_k, 2)
         native_iq = plan%mesh_index(k_targets(i_col))
         IF (native_iq > 0) sigma_k(:, i_col) = sigma_mesh(:, 1, native_iq)
      END DO
      !$OMP END PARALLEL DO
      DO i_col = 1, SIZE(sigma_k, 2)
         IF (plan%mesh_index(k_targets(i_col)) > 0) CYCLE
         n_off = n_off + 1
         offmesh_col(n_off) = i_col
      END DO
      IF (n_off == 0) RETURN

      nmesh = bs_env%kpoints_scf_desymm%nkp_grid
      lower = -nmesh/2
      cutoff_squared = -1.0_dp
      IF (apply_cutoff .AND. bs_env%ri_rs%cutoff_radius_v_w > 0.0_dp) THEN
         cutoff_squared = bs_env%ri_rs%cutoff_radius_v_w**2
      END IF
      first_c = plan%first_image
      last_c = plan%last_image
      n_c = SIZE(plan%phase_cos, 1)
      ALLOCATE (xk(3, n_off), cos_c(n_c, n_off), sin_c(n_c, n_off))
      DO i_off = 1, n_off
         i_col = k_targets(offmesh_col(i_off))
         i_slot = plan%offmesh_slot(i_col)
         xk(:, i_off) = bs_env%kpoints_DOS%xkp(:, i_col)
         cos_c(:, i_off) = plan%phase_cos(:, i_slot)
         sin_c(:, i_off) = plan%phase_sin(:, i_slot)
      END DO

      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(bs_env, cutoff_squared, cos_c, first_c, last_c, lower, n_c, n_off, nmesh) &
      !$OMP SHARED(offmesh_col, plan, sigma_k, sigma_R, sin_c, xk) &
      !$OMP PRIVATE(angle, c1, c2, c3, distance, i_c, i_first, i_local, i_off, i_pair, images, n_pairs, shift, &
      !$OMP         sum_im, sum_re, tied, weight, wrapped)
      ALLOCATE (images(n_c, pair_batch_size), sum_re(pair_batch_size, n_off), sum_im(pair_batch_size, n_off))
      !$OMP DO SCHEDULE(DYNAMIC)
      DO i_first = 1, SIZE(sigma_R, 1), pair_batch_size
         n_pairs = MIN(pair_batch_size, SIZE(sigma_R, 1) - i_first + 1)
         DO i_pair = 1, n_pairs
            i_local = i_first + i_pair - 1
            shift = plan%image_shift(:, i_local)
            tied = INT(plan%image_tie(:, i_local))
            i_c = 0
            DO c3 = first_c(3), last_c(3)
               DO c2 = first_c(2), last_c(2)
                  DO c1 = first_c(1), last_c(1)
                     i_c = i_c + 1
                     images(i_c, i_pair) = 0.0_dp
                     weight = plan%image_weights(c1 - first_c(1) + 1, 1, tied(1))* &
                              plan%image_weights(c2 - first_c(2) + 1, 2, tied(2))* &
                              plan%image_weights(c3 - first_c(3) + 1, 3, tied(3))
                     IF (weight == 0.0_dp) CYCLE
                     IF (cutoff_squared > 0.0_dp) THEN
                        distance = plan%displacement(:, i_local) - MATMUL(bs_env%hmat, REAL(shift + [c1, c2, c3], dp))
                        IF (SUM(distance**2) > cutoff_squared) CYCLE
                     END IF
                     wrapped = MODULO(shift + [c1, c2, c3] - lower, nmesh) + lower
                     images(i_c, i_pair) = weight*sigma_R(i_local, 1, 1 + wrapped(3) - lower(3) + &
                                                          nmesh(3)*(wrapped(2) - lower(2) + nmesh(2)*(wrapped(1) - lower(1))))
                  END DO
               END DO
            END DO
         END DO
         CALL dgemm('T', 'N', n_pairs, n_off, n_c, 1.0_dp, images, n_c, cos_c, n_c, 0.0_dp, sum_re, pair_batch_size)
         CALL dgemm('T', 'N', n_pairs, n_off, n_c, 1.0_dp, images, n_c, sin_c, n_c, 0.0_dp, sum_im, pair_batch_size)
         DO i_off = 1, n_off
            DO i_pair = 1, n_pairs
               angle = twopi*DOT_PRODUCT(xk(:, i_off), REAL(plan%image_shift(:, i_first + i_pair - 1), dp))
               sigma_k(i_first + i_pair - 1, offmesh_col(i_off)) = &
                  EXP(CMPLX(0.0_dp, angle, KIND=dp))*CMPLX(sum_re(i_pair, i_off), sum_im(i_pair, i_off), KIND=dp)
            END DO
         END DO
      END DO
      !$OMP END DO
      DEALLOCATE (images, sum_re, sum_im)
      !$OMP END PARALLEL
      DEALLOCATE (cos_c, offmesh_col, sin_c, xk)
   END SUBROUTINE band_self_energy

! **************************************************************************************************
!> \brief Prepare the common unordered-pair geometry and persistent two-real R2K input panels.
!> The first chi K2R stores its real fields in these panels; the last K2R replaces them with
!> the weighted product. Cell filtering is applied after the global Green norms are known.
!> \param bs_env ...
!> \param cells canonical physical cells
!> \param first_pair first owned unordered pair
!> \param n_pairs local unordered pairs
!> \param n_rank ranks sharing the k points and grid pairs
!> \param plan cached geometry, coefficients, and next-FFT input
!> \param apply_cutoff apply the configured W distance cutoff to the centered product images
!> \param q_indices requested output q points
!> \param q_partner inversion partners in the native W mesh
! **************************************************************************************************
   SUBROUTINE prepare_pair_fft(bs_env, cells, first_pair, n_pairs, n_rank, plan, apply_cutoff, q_indices, q_partner)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: cells(:, :), first_pair, n_pairs, n_rank
      TYPE(gw_pair_fft_work_type), INTENT(INOUT), TARGET :: plan
      LOGICAL, INTENT(IN)                                :: apply_cutoff
      INTEGER, INTENT(IN), OPTIONAL                      :: q_indices(:), q_partner(:)

      INTEGER :: count, fft_count, first, half_shift(3), i, i_col, i_extract, i_first, i_rank, &
         i_time, ibatch, icell, idir, ik, il, ilp, iq, isource, lower(3), n_extract, n_k, n_occ, &
         n_round, n_segments, ngrid, nmesh(3), pair_first, pair_last, partner, partner_native, &
         partner_output, rank_end, reflected(3), row_first, row_last, shift(3), width
      INTEGER(KIND=int_8)                                :: nfft_grid, spectrum_columns
      INTEGER, ALLOCATABLE :: extract_partner(:), extract_q(:), grid_linear(:), &
         input_to_native(:), native_to_input(:), partner_linear(:), q_position(:)
      LOGICAL                                            :: tie(3)
      REAL(KIND=dp)                                      :: cutoff_squared, displacement(3), &
                                                            energy_tolerance, half_tau, started, &
                                                            total_started
      REAL(KIND=dp), ALLOCATABLE                         :: fractional_grid(:, :)
      TYPE(gw_pair_fft_batch_type), POINTER              :: batch

      IF (fft_profile_enabled) total_started = m_walltime()
      CPASSERT(.NOT. ALLOCATED(plan%batches))
      ngrid = bs_env%ri_rs%n_grid_points
      nmesh = bs_env%kpoints_scf_desymm%nkp_grid
      lower = -nmesh/2
      ! GW uses the same native reciprocal mesh for Green functions, chi and W.
      IF (.NOT. ALL(bs_env%kpoints_chi_eps_W%xkp(:, 1:bs_env%nkp_scf_desymm) == bs_env%kpoints_scf_desymm%xkp)) THEN
         CPABORT('GW lattice transforms require one native reciprocal mesh')
      END IF
      n_k = bs_env%nkp_scf_desymm

      plan%k_batch = pair_exchange_width(ngrid, n_k, n_rank)
      CALL grid_pair_blocks(ngrid, n_rank, plan%pair_counts, plan%pair_displs)
      ! Columns and destination blocks partition the same packed triangle. Only a column
      ! or rank boundary can begin a new segment, so the map is linear in ngrid + n_rank.
      ALLOCATE (plan%column_segments(ngrid + 1), plan%segments(ngrid + n_rank))
      i_rank = 0
      n_segments = 0
      pair_first = 1
      DO i_col = 1, ngrid
         plan%column_segments(i_col) = n_segments + 1
         pair_last = pair_first + i_col - 1
         DO WHILE (i_rank < n_rank)
            rank_end = plan%pair_displs(i_rank) + plan%pair_counts(i_rank)
            row_first = MAX(pair_first, plan%pair_displs(i_rank) + 1) - pair_first + 1
            row_last = MIN(pair_last, rank_end) - pair_first + 1
            IF (row_first <= row_last) THEN
               n_segments = n_segments + 1
               plan%segments(n_segments) = pair_segment_type(i_rank, row_first, row_last, &
                                                             pair_first + (row_first - 1) - plan%pair_displs(i_rank) - 1)
            END IF
            IF (rank_end > pair_last) EXIT
            i_rank = i_rank + 1
         END DO
         IF (i_col < ngrid) pair_first = pair_last + 1
      END DO
      plan%column_segments(ngrid + 1) = n_segments + 1
      ALLOCATE (input_to_native(n_k), native_to_input(n_k), plan%green_source(n_k))
      i_extract = 0
      ! Receive contiguous source-rank blocks; the inverse FFT owns the resulting k permutation.
      DO i_first = 1, n_k, n_rank*plan%k_batch
         n_round = MIN(n_rank*plan%k_batch, n_k - i_first + 1)
         DO i_rank = 0, n_rank - 1
            DO ik = i_first + i_rank, i_first + n_round - 1, n_rank
               i_extract = i_extract + 1
               input_to_native(i_extract) = ik
               native_to_input(ik) = i_extract
               plan%green_source(i_extract) = i_extract
            END DO
         END DO
      END DO
      CPASSERT(i_extract == n_k)
      ALLOCATE (plan%product_phase(SIZE(cells, 2)))
      CALL k_grid_to_cell_prepare(plan%inverse, bs_env%kpoints_scf_desymm%xkp(:, input_to_native), nmesh, cells, &
                                  cell_phase_forward=plan%product_phase)
      CPASSERT(ALL(plan%inverse%nfft == nmesh([3, 2, 1])))
      CPASSERT(ALLOCATED(plan%inverse%k_partner))
      IF (PRESENT(q_partner)) THEN
         ALLOCATE (plan%q_source(n_k))
         CALL compact_q_slots(q_partner, n_rank, plan%q_slot)
         DO isource = 1, n_k
            ik = input_to_native(isource)
            plan%q_source(isource) = native_to_input(MIN(ik, q_partner(ik)))
         END DO
      END IF
      ! For the scalar KS problem H(-k)=H(k)* and S(-k)=S(k)*. Its spectral
      ! projectors give G(-k,τ)=G(k,τ)* when the weights and selected subspaces
      ! agree. Eigenvalue updates in evGW0 may give the two points different weights.
      IF (bs_env%gw_flavour == G0W0 .AND. bs_env%n_spin == 1 .AND. .NOT. bs_env%do_soc) THEN
         n_occ = bs_env%n_occ(1)
         ASSOCIATE (energies => bs_env%gwkp_env%eigenvalues_G(:, :, 1))
            DO isource = 1, n_k
               partner = plan%inverse%k_partner(isource)
               IF (partner <= isource) CYCLE
               ik = input_to_native(isource)
               partner_native = input_to_native(partner)
               energy_tolerance = 1024.0_dp*EPSILON(1.0_dp)* &
                                  MAX(1.0_dp, MAXVAL(ABS(energies(:, ik))), MAXVAL(ABS(energies(:, partner_native))))
               IF (MAXVAL(ABS(energies(:, ik) - energies(:, partner_native))) > energy_tolerance) CYCLE
               IF (n_occ > 0 .AND. n_occ < SIZE(energies, 1)) THEN
                  ! A projector cannot select part of a degenerate eigenspace.
                  IF (energies(n_occ + 1, ik) - energies(n_occ, ik) <= 2.0_dp*energy_tolerance) CYCLE
               END IF
               DO i_time = 1, bs_env%num_time_freq_points
                  half_tau = 0.5_dp*ABS(bs_env%time_frequency_grid%imaginary_time(i_time))
                  IF (ANY(ABS(ABS(half_tau*(energies(:, ik) - bs_env%e_fermi(1))) - bs_env%stabilize_exp) <= &
                          half_tau*energy_tolerance)) EXIT
               END DO
               IF (i_time <= bs_env%num_time_freq_points) CYCLE
               plan%green_source(partner) = isource
               plan%green_reuses = plan%green_reuses + 1
            END DO
         END ASSOCIATE
      END IF
      ! The packed forward transform keeps its factor of 1/2 for the two-real unpacking.
      plan%product_phase(:) = plan%product_phase/(2.0_dp*grid_pair_stencil_denominator)
      half_shift = MODULO(NINT(2.0_dp*REAL(nmesh, dp)*bs_env%kpoints_scf_desymm%xkp(:, 1)), 2)
      ! Forward FFT extraction keeps native q order, independently of the inverse input layout.
      n_extract = n_k
      IF (PRESENT(q_indices)) THEN
         IF (ANY(q_indices < 1 .OR. q_indices > n_k)) THEN
            CPABORT('Invalid GW FFT output q indices')
         END IF
         n_extract = SIZE(q_indices)
      END IF
      ! GNU Fortran can give nested vector subscripts a zero lower bound in ALLOCATE(SOURCE=...).
      ALLOCATE (grid_linear(n_extract), partner_linear(n_extract))
      DO iq = 1, n_extract
         ik = iq
         IF (PRESENT(q_indices)) ik = q_indices(iq)
         isource = native_to_input(ik)
         grid_linear(iq) = plan%inverse%k_slot(isource)
         partner_linear(iq) = plan%inverse%k_slot(plan%inverse%k_partner(isource))
      END DO
      CPASSERT(SIZE(cells, 2, KIND=int_8) == PRODUCT(INT(nmesh, int_8)))
      ALLOCATE (plan%mirror_cell(SIZE(cells, 2)), plan%diagonal(n_pairs))
      DO icell = 1, SIZE(cells, 2)
         isource = 1 + cells(3, icell) - lower(3) + &
                   nmesh(3)*(cells(2, icell) - lower(2) + nmesh(2)*(cells(1, icell) - lower(1)))
         IF (isource /= icell) CPABORT('GW product needs canonical cells')
         reflected = MODULO(-cells(:, icell) - lower, nmesh) + lower
         plan%mirror_cell(icell) = 1 + reflected(3) - lower(3) + &
                                   nmesh(3)*(reflected(2) - lower(2) + nmesh(2)*(reflected(1) - lower(1)))
      END DO
      ALLOCATE (fractional_grid(3, ngrid))
      fractional_grid(:, :) = MATMUL(bs_env%ri_rs%cell%h_inv, bs_env%ri_rs%grid_points)
      DO idir = 1, 3
         IF (bs_env%periodic(idir) == 0) fractional_grid(idir, :) = 0.0_dp
      END DO
      nfft_grid = PRODUCT(INT(plan%inverse%nfft, int_8))
      width = MAX(1, k_grid_to_cell_panel_width(plan%inverse, INT(n_pairs, int_8), fft_panel_bytes, real_output=.TRUE.))
      CPASSERT(nfft_grid*width <= HUGE(0))
      ALLOCATE (extract_q(n_extract), extract_partner(n_extract), SOURCE=0)
      spectrum_columns = 0_int_8
      ! When both q and -q are requested, one unpack produces both output columns.
      ALLOCATE (q_position(INT(nfft_grid)), SOURCE=0)
      DO ik = 1, SIZE(grid_linear)
         CPASSERT(q_position(grid_linear(ik)) == 0)
         q_position(grid_linear(ik)) = ik
      END DO
      n_extract = 0
      DO ik = 1, SIZE(grid_linear)
         partner_output = q_position(partner_linear(ik))
         IF (partner_output > 0 .AND. partner_output < ik) CYCLE
         n_extract = n_extract + 1
         extract_q(n_extract) = ik
         IF (partner_output /= ik) extract_partner(n_extract) = partner_output
      END DO
      DEALLOCATE (q_position)
      IF (fft_profile_enabled) THEN
         DO i_extract = 1, n_extract
            ik = extract_q(i_extract)
            spectrum_columns = spectrum_columns + MERGE(1_int_8, 2_int_8, grid_linear(ik) == partner_linear(ik))
         END DO
      END IF

      plan%n_pairs = n_pairs
      plan%width = width
      plan%nmesh = nmesh
      plan%nfft_grid = nfft_grid
      plan%spectrum_columns = spectrum_columns
      cutoff_squared = -1.0_dp
      IF (apply_cutoff .AND. bs_env%ri_rs%cutoff_radius_v_w > 0.0_dp) THEN
         cutoff_squared = bs_env%ri_rs%cutoff_radius_v_w**2
      END IF
      ALLOCATE (plan%batches((n_pairs + width - 1)/width))
      DO first = 1, n_pairs, width
         count = MIN(width, n_pairs - first + 1)
         fft_count = (count + 1)/2
         ibatch = (first - 1)/width + 1
         batch => plan%batches(ibatch)
         IF (fft_profile_enabled) started = m_walltime()
         CALL fft_allocate(batch%input, [INT(nfft_grid)*fft_count])
         IF (fft_profile_enabled) plan%profile%seconds(3) = plan%profile%seconds(3) + m_walltime() - started
         IF (fft_profile_enabled) started = m_walltime()
         ALLOCATE (batch%weight(count, SIZE(cells, 2)))
         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
         !$OMP SHARED(first, count, first_pair, fractional_grid, batch, bs_env, cells, nmesh, half_shift, cutoff_squared, plan) &
         !$OMP PRIVATE(i, il, ilp, shift, tie, displacement, icell)
         DO i = 1, count
            CALL grid_pair_indices(first_pair + first + i - 2, il, ilp)
            plan%diagonal(first + i - 1) = il == ilp
            CALL grid_pair_shift(il, ilp, fractional_grid, shift, tie)
            displacement = bs_env%ri_rs%grid_points(:, il) - bs_env%ri_rs%grid_points(:, ilp)
            DO icell = 1, SIZE(cells, 2)
               batch%weight(i, icell) = INT(grid_pair_stencil_weight(cells(:, icell), shift, tie, nmesh, half_shift, &
                                                                     displacement, bs_env%hmat, cutoff_squared), int_1)
            END DO
         END DO
         !$OMP END PARALLEL DO
         IF (fft_profile_enabled) plan%setup_seconds = plan%setup_seconds + m_walltime() - started
      END DO
      CALL MOVE_ALLOC(grid_linear, plan%grid_linear)
      CALL MOVE_ALLOC(partner_linear, plan%partner_linear)
      CALL MOVE_ALLOC(input_to_native, plan%input_to_native)
      ALLOCATE (plan%extract_q, SOURCE=extract_q(:n_extract))
      ALLOCATE (plan%extract_partner, SOURCE=extract_partner(:n_extract))
      IF (fft_profile_enabled) plan%profile%seconds(2) = m_walltime() - total_started - plan%profile%seconds(3)
   END SUBROUTINE prepare_pair_fft

! **************************************************************************************************
!> \brief Execute prepared pair-product panels after the global physical-cell filtering decision.
!> No real-space gather or packing pass is needed; active input panels are already in FFT order.
!> \param plan prepared geometry and product input
!> \param values_k compact unordered-pair output
!> \param cell_active full-cell mask from the global original-Green norms
!> \param profile_label chi or self-energy
!> \param values_R physical product retained for off-mesh self-energy projection
! **************************************************************************************************
   SUBROUTINE execute_pair_fft(plan, values_k, cell_active, profile_label, values_R)
      TYPE(gw_pair_fft_work_type), INTENT(INOUT), TARGET :: plan
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(OUT)          :: values_k(:, :, :)
      LOGICAL, INTENT(IN)                                :: cell_active(:)
      CHARACTER(LEN=*), INTENT(IN)                       :: profile_label
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: values_R(:, :, :)

      COMPLEX(KIND=dp)                                   :: result_a, result_b, spectrum_a, &
                                                            spectrum_b
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: fft_in(:), fft_out(:), input(:, :), &
                                                            output(:, :), scratch(:)
      INTEGER                                            :: count, fft_count, first, i, i_extract, &
                                                            ibatch, icell, ik, n_extract, nfft(3), &
                                                            nmesh(3), partner_output, unit_nr, &
                                                            width
      INTEGER(KIND=int_8)                                :: extract_reads, fft_transforms, &
                                                            nfft_grid, spectrum_columns
      INTEGER, POINTER                                   :: extract_partner(:), extract_q(:), &
                                                            grid_linear(:), partner_linear(:)
      REAL(KIND=dp)                                      :: elapsed, other_seconds, &
                                                            prepared_seconds, started, &
                                                            total_started
      TYPE(lattice_fft_profile_type), ALLOCATABLE        :: profile

      IF (fft_profile_enabled) total_started = m_walltime()
      CPASSERT(ALLOCATED(plan%batches))
      CPASSERT(SIZE(cell_active) == SIZE(plan%inverse%cell_slot))
      CPASSERT(SIZE(values_k, 1) == plan%n_pairs)
      CPASSERT(SIZE(values_k, 2) == 1)
      CPASSERT(SIZE(values_k, 3) == SIZE(plan%grid_linear))
      IF (fft_profile_enabled) THEN
         ALLOCATE (profile, SOURCE=plan%profile)
         prepared_seconds = SUM(profile%seconds(1:3))
      END IF
      nfft = plan%inverse%nfft
      nmesh = plan%nmesh
      nfft_grid = plan%nfft_grid
      width = plan%width
      n_extract = SIZE(plan%extract_q)
      spectrum_columns = plan%spectrum_columns
      grid_linear => plan%grid_linear
      partner_linear => plan%partner_linear
      extract_q => plan%extract_q
      extract_partner => plan%extract_partner
      fft_transforms = 0_int_8
      extract_reads = 0_int_8
      IF (fft_profile_enabled) started = m_walltime()
      CALL fft_allocate(scratch, [INT(nfft_grid)*((width + 1)/2)])
      IF (fft_profile_enabled) profile%seconds(3) = profile%seconds(3) + m_walltime() - started
      DO first = 1, plan%n_pairs, width
         count = MIN(width, plan%n_pairs - first + 1)
         fft_count = (count + 1)/2
         IF (fft_profile_enabled) fft_transforms = fft_transforms + INT(fft_count, int_8)
         ibatch = (first - 1)/width + 1
         IF (fft_profile_enabled) THEN
            profile%elements = profile%elements + INT(count, int_8)
            profile%batches = profile%batches + 1_int_8
            profile%tile_min = MIN(profile%tile_min, count)
            profile%tile_max = MAX(profile%tile_max, count)
         END IF
         fft_in => plan%batches(ibatch)%input
         fft_out => scratch
         input(1:fft_count, 1:INT(nfft_grid)) => fft_in
         IF (fft_profile_enabled) started = m_walltime()
         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
         !$OMP SHARED(cell_active, input, plan, first, count, values_R) PRIVATE(icell)
         DO icell = 1, SIZE(cell_active)
            IF (cell_active(icell)) CYCLE
            input(:, plan%inverse%cell_slot(icell)) = z_zero
            IF (PRESENT(values_R)) values_R(first:first + count - 1, 1, icell) = 0.0_dp
         END DO
         !$OMP END PARALLEL DO
         IF (fft_profile_enabled) THEN
            profile%seconds(4) = profile%seconds(4) + m_walltime() - started
            profile%pair_major = .TRUE.
         END IF
         CALL lattice_fft_execute_pair_major(BWFFT, nfft, fft_count, fft_in, fft_out, profile)
         output(1:fft_count, 1:INT(nfft_grid)) => fft_in
         IF (fft_profile_enabled) started = m_walltime()
         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
         !$OMP SHARED(values_k, output, grid_linear, partner_linear, first, count, fft_count) &
         !$OMP SHARED(n_extract, extract_q, extract_partner) &
         !$OMP PRIVATE(ik, partner_output, i, spectrum_a, spectrum_b, result_a, result_b)
         DO i_extract = 1, n_extract
            ik = extract_q(i_extract)
            partner_output = extract_partner(i_extract)
            ! The cached stencil packs Z=(F_a+i*F_b)/2. Extract with real
            ! additions/subtractions only; partner_linear includes the MP shift.
            IF (grid_linear(ik) == partner_linear(ik)) THEN
               !$OMP SIMD PRIVATE(spectrum_a)
               DO i = 1, count/2
                  spectrum_a = output(i, grid_linear(ik))
                  values_k(first + 2*i - 2, 1, ik) = CMPLX(2.0_dp*REAL(spectrum_a, dp), 0.0_dp, dp)
                  values_k(first + 2*i - 1, 1, ik) = CMPLX(2.0_dp*AIMAG(spectrum_a), 0.0_dp, dp)
               END DO
            ELSE IF (partner_output > 0) THEN
               !$OMP SIMD PRIVATE(spectrum_a, spectrum_b, result_a, result_b)
               DO i = 1, count/2
                  spectrum_a = output(i, grid_linear(ik))
                  spectrum_b = output(i, partner_linear(ik))
                  result_a = CMPLX(REAL(spectrum_a, dp) + REAL(spectrum_b, dp), &
                                   AIMAG(spectrum_a) - AIMAG(spectrum_b), dp)
                  result_b = CMPLX(AIMAG(spectrum_a) + AIMAG(spectrum_b), &
                                   REAL(spectrum_b, dp) - REAL(spectrum_a, dp), dp)
                  values_k(first + 2*i - 2, 1, ik) = result_a
                  values_k(first + 2*i - 1, 1, ik) = result_b
                  values_k(first + 2*i - 2, 1, partner_output) = CONJG(result_a)
                  values_k(first + 2*i - 1, 1, partner_output) = CONJG(result_b)
               END DO
            ELSE
               !$OMP SIMD PRIVATE(spectrum_a, spectrum_b)
               DO i = 1, count/2
                  spectrum_a = output(i, grid_linear(ik))
                  spectrum_b = output(i, partner_linear(ik))
                  values_k(first + 2*i - 2, 1, ik) = &
                     CMPLX(REAL(spectrum_a, dp) + REAL(spectrum_b, dp), AIMAG(spectrum_a) - AIMAG(spectrum_b), dp)
                  values_k(first + 2*i - 1, 1, ik) = &
                     CMPLX(AIMAG(spectrum_a) + AIMAG(spectrum_b), REAL(spectrum_b, dp) - REAL(spectrum_a, dp), dp)
               END DO
            END IF
            IF (MODULO(count, 2) /= 0) THEN
               ! The lone transform is also half-scaled; preserve its complex roundoff.
               values_k(first + count - 1, 1, ik) = 2.0_dp*output(fft_count, grid_linear(ik))
               IF (partner_output > 0) THEN
                  values_k(first + count - 1, 1, partner_output) = CONJG(values_k(first + count - 1, 1, ik))
               END IF
            END IF
         END DO
         !$OMP END PARALLEL DO
         IF (fft_profile_enabled) THEN
            profile%seconds(8) = profile%seconds(8) + m_walltime() - started
            extract_reads = extract_reads + INT(count/2, int_8)*spectrum_columns + &
                            INT(MODULO(count, 2), int_8)*INT(n_extract, int_8)
         END IF
      END DO

      IF (fft_profile_enabled) started = m_walltime()
      CALL fft_deallocate(scratch)
      unit_nr = 0
      IF (fft_profile_enabled) THEN
         profile%seconds(9) = m_walltime() - started
         elapsed = m_walltime() - total_started + prepared_seconds
         other_seconds = MAX(0.0_dp, elapsed - profile%seconds(4) - SUM(profile%seconds(5:8)))
         unit_nr = cp_logger_get_default_io_unit()
      END IF
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') 'R2K_FFT_PROFILE| GW unordered-pair panels; logger-rank wall seconds'
         WRITE (unit_nr, '(T2,A,A)') 'R2K_FFT_PROFILE| operation: ', profile_label
         WRITE (unit_nr, '(T2,A,L2)') 'R2K_FFT_PROFILE| input packed by K2R extraction: ', .TRUE.
         WRITE (unit_nr, '(T2,A,2I10)') 'R2K_FFT_PROFILE| active / full input cells: ', &
            SUM(MERGE(1, 0, cell_active)), SIZE(cell_active)
         WRITE (unit_nr, '(T2,A,3I10)') 'R2K_FFT_PROFILE| source shape: ', plan%n_pairs, 1, SIZE(cell_active)
         WRITE (unit_nr, '(T2,A,2I10)') 'R2K_FFT_PROFILE| extracted / full q points: ', &
            SIZE(values_k, 3), SIZE(plan%inverse%k_slot)
         WRITE (unit_nr, '(T2,A,2I10)') 'R2K_FFT_PROFILE| unpack q jobs / conjugate reuses: ', &
            n_extract, SIZE(values_k, 3) - n_extract
         WRITE (unit_nr, '(T2,A,2F14.6)') 'R2K_FFT_PROFILE| extract read/write [GB, logical]: ', &
            complex_bytes*REAL(extract_reads, dp)/bytes_per_GB, &
            complex_bytes*REAL(profile%elements, dp)*REAL(SIZE(values_k, 3), dp)/bytes_per_GB
         WRITE (unit_nr, '(T2,A)') 'R2K_FFT_PROFILE| pair-major FFT-X/plan includes the full 3D batched FFT'
         WRITE (unit_nr, '(T2,A,6I8)') 'R2K_FFT_PROFILE| grid / padded: ', nmesh, nfft
         WRITE (unit_nr, '(T2,A,L2)') 'R2K_FFT_PROFILE| two real pairs per complex FFT: ', .TRUE.
         WRITE (unit_nr, '(T2,A,2I16)') 'R2K_FFT_PROFILE| FFT transforms / real pairs: ', &
            fft_transforms, profile%elements
         WRITE (unit_nr, '(T2,A,2I16)') 'R2K_FFT_PROFILE| complex FFT elements / full: ', &
            nfft_grid*fft_transforms, nfft_grid*profile%elements
         WRITE (unit_nr, '(T2,A,2I16,2I8)') 'R2K_FFT_PROFILE| elements/batches/tile min/max: ', &
            profile%elements, profile%batches, profile%tile_min, profile%tile_max
         WRITE (unit_nr, '(T2,A,9F12.6)') 'R2K_FFT_PROFILE| init/geometry/buffers/fill/X/Y/Z/extract/release: ', profile%seconds
         WRITE (unit_nr, '(T2,A,F14.6)') 'R2K_FFT_PROFILE| stencil-setup (included in geometry): ', plan%setup_seconds
         WRITE (unit_nr, '(T2,A,I16)') 'R2K_FFT_PROFILE| cached real stencil [bytes]: ', &
            INT(int_1_size, int_8)*INT(plan%n_pairs, int_8)*SIZE(cell_active, KIND=int_8)
         WRITE (unit_nr, '(T2,A,L2)') 'R2K_FFT_PROFILE| reused stencil cache: ', plan%cache_hit
         WRITE (unit_nr, '(T2,A,5F12.6)') 'R2K_FFT_PROFILE| grid-fill/FFT+plan/extract/other/total: ', &
            profile%seconds(4), SUM(profile%seconds(5:7)), profile%seconds(8), other_seconds, elapsed
         WRITE (unit_nr, '(T2,A,F14.6)') 'R2K_FFT_PROFILE| total         ', elapsed
      END IF
      plan%profile = lattice_fft_profile_type()
      plan%setup_seconds = 0.0_dp
      plan%cache_hit = .TRUE.
   END SUBROUTINE execute_pair_fft

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
!> \brief Transfer unordered grid pairs into bounded q batches and consume the Hermitian
!>        panels as Z(q)^H chi(q) Z(q), without expanding a full grid matrix.
!> \param fft_pairs compact needed-q columns, ordered by ascending original q index
!> \param para_env ...
!> \param needed q points to transfer
!> \param Z_q interpolation matrices at the q points owned by this rank, MOD(q,P)
!> \param chi_q accumulate projected matrices in compact local q order
!> \param batch_rounds common maximum number of communication rounds per projection batch
!> \param projection_seconds logger-rank projection wall time, excluding transfers
! **************************************************************************************************
   SUBROUTINE project_fft_pairs_to_chi(fft_pairs, para_env, needed, Z_q, chi_q, batch_rounds, projection_seconds)
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(IN), TARGET   :: fft_pairs(:, :, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: needed(:)
      COMPLEX(KIND=dp), INTENT(IN)                       :: Z_q(:, :, :)
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: chi_q(:, :, :)
      INTEGER, INTENT(IN)                                :: batch_rounds
      REAL(KIND=dp), INTENT(OUT)                         :: projection_seconds

      CHARACTER(LEN=*), PARAMETER :: routineN = 'project_fft_pairs_to_chi'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: panel(:, :, :), work(:, :, :)
      COMPLEX(KIND=dp), ALLOCATABLE, TARGET              :: q_pairs(:, :)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: received(:), sent(:)
      COMPLEX(KIND=dp), TARGET                           :: unused_receive(1), unused_send(1)
      INTEGER :: capacity, first_round, handle, i_buffer, i_column, i_first, i_local, i_needed, &
         i_rank, i_round, i_thread, iq, n_aux, n_buffer, n_grid, n_local, n_needed_local, n_q, &
         n_rank, n_round, n_selected, n_workers
      INTEGER, ALLOCATABLE                               :: pair_counts(:), pair_displs(:), &
                                                            recv_counts(:), send_counts(:), &
                                                            send_displs(:)
      LOGICAL                                            :: receiving
      REAL(KIND=dp)                                      :: started

      CALL timeset(routineN, handle)

      n_rank = para_env%num_pe
      n_grid = SIZE(Z_q, 1)
      n_aux = SIZE(Z_q, 2)
      n_local = SIZE(fft_pairs, 1)
      n_q = SIZE(needed)
      CPASSERT(SIZE(fft_pairs, 2) == 1)
      CPASSERT(SIZE(fft_pairs, 3) == COUNT(needed))
      CPASSERT(batch_rounds > 0)
      n_needed_local = COUNT([(needed(iq) .AND. MODULO(iq, n_rank) == para_env%mepos, iq=1, n_q)])
      CPASSERT(ALL(SHAPE(chi_q) == [n_aux, n_aux, n_needed_local]))
      capacity = MIN(batch_rounds, n_needed_local)
      n_workers = 1
!$    n_workers = omp_get_max_threads()
      n_workers = MAX(1, MIN(n_workers, capacity))
      CALL grid_pair_blocks(n_grid, n_rank, pair_counts, pair_displs)
      ALLOCATE (recv_counts(0:n_rank - 1), &
                send_counts(0:n_rank - 1), send_displs(0:n_rank - 1), &
                q_pairs(INT(INT(n_grid, int_8)*(n_grid + 1)/2), capacity), work(n_grid, n_aux, n_workers), &
                panel(n_grid, MIN(n_grid, grid_pair_panel_cols), n_workers))
      unused_receive = z_zero
      unused_send = z_zero
      i_column = 0
      i_local = 0
      projection_seconds = 0.0_dp
      DO first_round = 1, (n_q + n_rank - 1)/n_rank, batch_rounds
         n_buffer = 0
         DO i_round = first_round, MIN(first_round + batch_rounds - 1, (n_q + n_rank - 1)/n_rank)
            i_first = (i_round - 1)*n_rank
            n_round = MIN(n_rank, n_q - i_first)
            send_counts(:) = 0
            send_displs(:) = 0
            n_selected = 0
            DO iq = i_first + 1, i_first + n_round
               IF (.NOT. needed(iq)) CYCLE
               i_rank = MODULO(iq, n_rank)
               send_counts(i_rank) = n_local
               send_displs(i_rank) = n_selected*n_local
               n_selected = n_selected + 1
            END DO
            IF (n_selected == 0) CYCLE
            ! Compact q columns are contiguous even when their original q indices have gaps.
            sent => unused_send
            IF (n_local > 0) THEN
               sent(1:n_local*n_selected) => fft_pairs(:, :, i_column + 1:i_column + n_selected)
            END IF
            i_column = i_column + n_selected
            iq = i_first + MODULO(para_env%mepos - 1, n_rank) + 1
            receiving = iq <= n_q
            IF (receiving) receiving = needed(iq)
            recv_counts(:) = 0
            received => unused_receive
            IF (receiving) THEN
               recv_counts(:) = pair_counts
               i_local = i_local + 1
               n_buffer = n_buffer + 1
               CPASSERT(n_buffer <= capacity)
               received => q_pairs(:, n_buffer)
            END IF
            CALL para_env%alltoall(sent, send_counts, send_displs, received, recv_counts, pair_displs)
         END DO
         IF (n_buffer == 0) CYCLE
         IF (fft_profile_enabled) started = m_walltime()
         ! Distinct q matrices remain independent; each thread reuses its GEMM workspace.
         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) IF(n_buffer > 1) NUM_THREADS(n_workers) &
         !$OMP SHARED(n_buffer, n_grid, n_aux, q_pairs, Z_q, chi_q, panel, work, i_local) &
         !$OMP PRIVATE(i_buffer, i_thread, i_needed)
         DO i_buffer = 1, n_buffer
            i_thread = 1
!$          i_thread = omp_get_thread_num() + 1
            i_needed = i_local - n_buffer + i_buffer
            CALL grid_pair_upper_multiply(n_grid, n_aux, q_pairs(:, i_buffer), Z_q(:, :, i_needed), &
                                          work(:, :, i_thread), panel(:, :, i_thread))
            ! Z^H H Z = Z^H (U Z) + (U Z)^H Z; accumulate only the Hermitian upper triangle.
            CALL zher2k('U', 'C', n_aux, n_grid, z_one, Z_q(:, :, i_needed), n_grid, &
                        work(:, :, i_thread), n_grid, 1.0_dp, chi_q(:, :, i_needed), n_aux)
         END DO
         !$OMP END PARALLEL DO
         IF (fft_profile_enabled) projection_seconds = projection_seconds + m_walltime() - started
      END DO
      CPASSERT(i_column == SIZE(fft_pairs, 3))
      CPASSERT(i_local == n_needed_local)

      CALL timestop(handle)

   END SUBROUTINE project_fft_pairs_to_chi

! **************************************************************************************************
!> \brief Project W as Z W Z^H or the bare interaction as (ZK)(ZK)^H, then exchange its panels.
!> D(q) and d(q) come from real-space real coefficients, so the RI fit obeys Z(-q)=Z(q)*.
!> Apply the same conjugation to the projected W after exchanging it with the grid-pair owners.
!> \param Z_q ...
!> \param interaction_q Compact local W(q,iτ) matrices or one K(q) factor per q
!> \param itime ...
!> \param para_env ...
!> \param work shared geometry and exchange layout
!> \param fft_pairs ...
!> \param factor whether interaction_q contains Coulomb factors
! **************************************************************************************************
   SUBROUTINE project_auxiliary_q_to_fft_pairs(Z_q, interaction_q, itime, para_env, work, fft_pairs, factor)
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(IN), TARGET   :: Z_q(:, :, :)
      TYPE(ri_q_matrix_block_type), INTENT(IN)           :: interaction_q(:)
      INTEGER, INTENT(IN)                                :: itime
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(gw_pair_fft_work_type), INTENT(IN)            :: work
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(OUT), TARGET  :: fft_pairs(:, :, :)
      LOGICAL, INTENT(IN)                                :: factor

      CHARACTER(LEN=*), PARAMETER :: routineN = 'project_auxiliary_q_to_fft_pairs'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: panel(:, :), sent(:)
      COMPLEX(KIND=dp), ALLOCATABLE, TARGET              :: ZW(:, :, :)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: right(:, :)
      INTEGER :: first_col, handle, i_first, i_q, i_z, n_aux_cols, n_cols, n_grid, n_owned, &
         n_pairs, n_q, n_q_local, n_rank, n_RI, n_round, n_scheduled, n_threads, ri_block
      INTEGER, ALLOCATABLE                               :: jobs(:), recv_counts(:), recv_displs(:)

      CALL timeset(routineN, handle)

      n_grid = SIZE(Z_q, 1)
      n_RI = SIZE(Z_q, 2)
      n_q_local = SIZE(Z_q, 3)
      n_q = SIZE(fft_pairs, 3)
      n_rank = para_env%num_pe
      n_pairs = SUM(work%pair_counts)
      CPASSERT(SIZE(fft_pairs, 1) == work%pair_counts(para_env%mepos))
      CPASSERT(ALLOCATED(work%q_source))
      CPASSERT(SIZE(interaction_q) == n_q_local)
      ALLOCATE (sent(MAX(1, n_pairs*MIN(n_q_local, work%k_batch))), ZW(n_grid, n_RI, MIN(n_q_local, work%k_batch)), &
                jobs(work%k_batch), recv_counts(0:n_rank - 1), recv_displs(0:n_rank - 1))
      n_threads = 1
!$    n_threads = omp_get_max_threads()
      ri_block = MAX(64, (n_RI + n_threads - 1)/n_threads)

      DO i_first = 1, n_q, n_rank*work%k_batch
         n_round = MIN(n_rank*work%k_batch, n_q - i_first + 1)
         CALL pair_exchange_jobs(work, work%q_source, para_env, i_first, n_round, 1, &
                                 jobs, n_owned, n_scheduled, recv_counts, recv_displs)
         IF (n_scheduled == 0) CYCLE

         !$OMP PARALLEL DEFAULT(NONE) IF(n_owned > 0) &
         !$OMP SHARED(Z_q, ZW, interaction_q, itime, work, sent, jobs, n_owned, n_grid, n_RI, ri_block, factor) &
         !$OMP PRIVATE(first_col, i_q, i_z, n_aux_cols, n_cols, panel, right)
         ALLOCATE (panel(n_grid, MIN(n_grid, grid_pair_panel_cols)))
         !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
         DO i_q = 1, n_owned
            DO first_col = 1, n_RI, ri_block
               i_z = work%q_slot(jobs(i_q))
               n_aux_cols = MIN(ri_block, n_RI - first_col + 1)
               CALL zgemm('N', 'N', n_grid, n_aux_cols, n_RI, z_one, Z_q(:, :, i_z), n_grid, &
                          interaction_q(i_z)%values(1, first_col, itime), n_RI, z_zero, ZW(1, first_col, i_q), n_grid)
            END DO
         END DO
         !$OMP END DO
         !$OMP DO COLLAPSE(2) SCHEDULE(STATIC, 1)
         DO i_q = 1, n_owned
            DO first_col = 1, n_grid, grid_pair_panel_cols
               i_z = work%q_slot(jobs(i_q))
               n_cols = MIN(grid_pair_panel_cols, n_grid - first_col + 1)
               right => Z_q(:, :, i_z)
               IF (factor) right => ZW(:, :, i_q)
               CALL grid_pair_gemm_panel(n_grid, n_RI, first_col, n_cols, ZW(:, :, i_q), right, panel)
               CALL pack_grid_pair_panel(work, first_col, n_cols, panel, n_owned, i_q, sent)
            END DO
         END DO
         !$OMP END DO
         DEALLOCATE (panel)
         !$OMP END PARALLEL

         CALL exchange_pair_round(work, work%q_source, para_env, i_first, n_round, 1, &
                                  n_owned, sent, recv_counts, recv_displs, fft_pairs)
      END DO
      CALL restore_pair_conjugates(work%q_source, fft_pairs)

      CALL timestop(handle)

   END SUBROUTINE project_auxiliary_q_to_fft_pairs

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
!> \brief Number representative q points in native order within each cyclic owner.
!> \param q_partner ...
!> \param n_rank ...
!> \param q_slot local slot for each q, with zero for its omitted partner
! **************************************************************************************************
   SUBROUTINE compact_q_slots(q_partner, n_rank, q_slot)
      INTEGER, INTENT(IN)                                :: q_partner(:), n_rank
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: q_slot(:)

      INTEGER                                            :: iq, owner, counts(0:n_rank - 1)

      ALLOCATE (q_slot(SIZE(q_partner)), SOURCE=0)
      counts(:) = 0
      DO iq = 1, SIZE(q_partner)
         IF (q_partner(iq) < iq) CYCLE
         owner = MODULO(iq, n_rank)
         counts(owner) = counts(owner) + 1
         q_slot(iq) = counts(owner)
      END DO
   END SUBROUTINE compact_q_slots

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
!> The Cholesky factor ε = UᴴU also gives the RPA integrand Tr[ln ε - ε + I] = Σ_P [2 ln U_PP - ε_PP + 1].
!> \param matrices ...
!> \param group_env ...
!> \param n ...
!> \param root ...
!> \param chi_owner ...
!> \param cutoff ...
!> \param output_owner W on the root rank, only computed for GW
!> \param rpa_integrand Tr[ln ε - ε + I] on all ranks of the group
!> \param do_gw W is computed only for GW; RPA needs only the RPA integrand
! **************************************************************************************************
   SUBROUTINE grouped_w_solve(matrices, group_env, n, root, chi_owner, cutoff, output_owner, rpa_integrand, &
                              do_gw)
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: matrices(n_w_matrices)
      TYPE(mp_para_env_type), INTENT(INOUT)              :: group_env
      INTEGER, INTENT(IN)                                :: n, root
      COMPLEX(KIND=dp), INTENT(IN)                       :: chi_owner(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: output_owner(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: rpa_integrand
      LOGICAL, INTENT(IN)                                :: do_gw

      COMPLEX(KIND=dp), ALLOCATABLE                      :: diag_eps(:), diag_U(:)
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
      ALLOCATE (diag_eps(n), diag_U(n))
      CALL cp_cfm_get_diag(matrices(w_dielectric), diag_eps)
      ! χ is positive semidefinite after filtering: ε=I+BᴴχB has λ_min>=1.
      CALL cp_cfm_cholesky_decompose(matrices(w_dielectric), info_out=info)
      IF (info /= 0) THEN
         CPABORT("Screened dielectric Cholesky factorization failed")
      END IF
      CALL cp_cfm_get_diag(matrices(w_dielectric), diag_U)
      rpa_integrand = SUM(2.0_dp*LOG(REAL(diag_U, dp)) - REAL(diag_eps, dp) + 1.0_dp)
      DEALLOCATE (diag_eps, diag_U)
      IF (.NOT. do_gw) RETURN
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
      ! GW releases the AO images in the Gamma contraction of W; RPA never contracts them.
      IF (ASSOCIATED(gwkp_env%gamma_corr_d)) CALL dbcsr_deallocate_matrix_set(gwkp_env%gamma_corr_d)
      IF (ASSOCIATED(gwkp_env%gamma_corr_h)) CALL dbcsr_deallocate_matrix_set(gwkp_env%gamma_corr_h)
      IF (ASSOCIATED(gwkp_env%gamma_corr_s)) CALL dbcsr_deallocate_matrix_set(gwkp_env%gamma_corr_s)
      IF (ASSOCIATED(gwkp_env%gamma_corr_index_to_cell)) DEALLOCATE (gwkp_env%gamma_corr_index_to_cell)
      IF (ASSOCIATED(gwkp_env%gamma_corr_cell_to_index)) DEALLOCATE (gwkp_env%gamma_corr_cell_to_index)
      DEALLOCATE (bs_env%gwkp_env)
      NULLIFY (gwkp_env)

      CALL timestop(handle)

   END SUBROUTINE release_gw_work

END MODULE gw_ri_rs_small_cell_full_kp_utils

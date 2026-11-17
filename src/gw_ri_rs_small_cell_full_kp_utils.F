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
   USE fft_lib,                         ONLY: fft_pair_major_plan_destroy,&
                                              fft_pair_major_plan_type
   USE fft_tools,                       ONLY: BWFFT,&
                                              fft_alloc,&
                                              fft_dealloc
   USE gw_ri_rs_small_cell_full_kp_lattice,ONLY: &
        assemble_AA_contraction, block_range, build_full_cell_grid, build_ri_2c_lattice, &
        complex_bytes, gamma_auxiliary_charges, grid_pair_gemm, grid_pair_indices, &
        grid_pair_panel_cols, grid_pair_shift, grid_pair_stencil_denominator, &
        grid_pair_stencil_weight, grid_pair_upper_multiply, image_weight, owned_mo_coefficients, &
        real_bytes, ri_2c_lattice_to_kp, ri_metric_head_zero
   USE gw_ri_rs_small_cell_full_kp_q0_corr,ONLY: add_gamma_corr_pairs,&
                                              gamma_point_correction,&
                                              q_mesh_and_gamma_shell_average,&
                                              singular_head_cell_factors
   USE gw_ri_rs_small_cell_full_kp_types,ONLY: gwkp_type,&
                                              ri_q_matrix_block_type
   USE gw_utils_fm,                     ONLY: local_complex_power
   USE input_constants,                 ONLY: evGW0
   USE kinds,                           ONLY: dp,&
                                              int_1,&
                                              int_1_size,&
                                              int_8
   USE kpoint_coulomb_2c,               ONLY: add_ewald_reciprocal
   USE kpoint_lattice_fft,              ONLY: k_grid_to_cell_execute,&
                                              k_grid_to_cell_panel_width,&
                                              k_grid_to_cell_prepare,&
                                              k_grid_to_cell_release,&
                                              k_grid_to_cell_work_type,&
                                              lattice_fft_execute_pair_major,&
                                              lattice_fft_profile_type,&
                                              lattice_fft_real_sink_type
   USE kpoint_transform,                ONLY: cell_to_k_transform
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
   USE ri_rs_kp_integrals,              ONLY: write_ri_rs_progress
   USE ri_rs_utils,                     ONLY: precompute_ri_rs_radii

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
   ! Fixed workspace of band projection, received chi batches, and Σ(k) of one projection round.
   REAL(KIND=dp), PARAMETER :: projection_workspace = 64.0_dp*bytes_per_MiB, &
                               projection_round_bytes = 512.0_dp*bytes_per_MiB
   INTEGER(KIND=int_8), PARAMETER :: fft_panel_bytes = 32_int_8*1024_int_8**2

   ! One stencil representation: signed bytes with denominator 64 and a shared cell phase.
   TYPE gw_pair_fft_batch_type
      INTEGER(KIND=int_1), ALLOCATABLE :: weight(:, :)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER :: input(:) => NULL()
   CONTAINS
      FINAL :: release_pair_fft_batch
   END TYPE gw_pair_fft_batch_type

   TYPE gw_pair_fft_work_type
      PRIVATE
      INTEGER :: n_pairs = 0, width = 0, k_batch = 0, nmesh(3) = 0
      INTEGER(KIND=int_8) :: nfft_grid = 0_int_8, spectrum_columns = 0_int_8
      LOGICAL :: cache_hit = .FALSE.
      COMPLEX(KIND=dp), ALLOCATABLE :: product_phase(:)
      INTEGER, ALLOCATABLE :: grid_linear(:), partner_linear(:), extract_q(:), extract_partner(:), mirror_cell(:)
      INTEGER, ALLOCATABLE :: input_to_native(:)
      LOGICAL, ALLOCATABLE :: diagonal(:)
      TYPE(k_grid_to_cell_work_type) :: inverse
      TYPE(gw_pair_fft_batch_type), ALLOCATABLE :: batches(:)
      TYPE(lattice_fft_profile_type) :: profile
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
!> \brief Release the aligned input owned by one persistent FFT panel.
!> \param batch ...
! **************************************************************************************************
   IMPURE ELEMENTAL SUBROUTINE release_pair_fft_batch(batch)
      TYPE(gw_pair_fft_batch_type), INTENT(INOUT)        :: batch

      IF (ASSOCIATED(batch%input)) CALL fft_dealloc(batch%input)
   END SUBROUTINE release_pair_fft_batch

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
!> \param two_real whether real and imaginary channels contain independent real fields
! **************************************************************************************************
   SUBROUTINE extract_grid_pairs(sink, first, nbatch, spectrum, cell_slot, cell_phase, two_real)
      CLASS(gw_pair_extraction_type), INTENT(INOUT) :: sink
      INTEGER(KIND=int_8), INTENT(IN) :: first
      INTEGER, INTENT(IN) :: nbatch, cell_slot(:)
      COMPLEX(KIND=dp), INTENT(IN) :: spectrum(:, :), cell_phase(:)
      LOGICAL, INTENT(IN) :: two_real

      COMPLEX(KIND=dp) :: input_phase, phase, raw, saved
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER :: input(:, :)
      INTEGER :: fft_count, i, ia, ib, icell, ipair, jpair, slot
      REAL(KIND=dp) :: a, b, norm, offnorm, product_a, product_b
      LOGICAL :: collect_norm
      TYPE(gw_pair_fft_batch_type), POINTER :: batch
      TYPE(gw_pair_fft_work_type), POINTER :: plan

      fft_count = (nbatch + 1)/2
      IF (.NOT. two_real) THEN
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
      !$OMP SHARED(sink, first, nbatch, fft_count, spectrum, cell_slot, cell_phase, two_real, plan, batch, input) &
      !$OMP SHARED(collect_norm) &
      !$OMP PRIVATE(icell, slot, phase, i, ia, ib, ipair, jpair, raw, saved, A, b, product_a, product_b, norm, offnorm) &
      !$OMP PRIVATE(input_phase)
      DO icell = 1, SIZE(cell_slot)
         slot = cell_slot(icell)
         phase = MERGE(0.5_dp, 1.0_dp, two_real)*cell_phase(icell)
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
      TYPE(lattice_fft_profile_type)                     :: profile

      total_started = m_walltime()
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
      profile = lattice_fft_profile_type()
      CALL k_grid_to_cell_execute(work%inverse, values_k, require_fft=.TRUE., profile=profile, &
                                  max_workspace_bytes=fft_panel_bytes, real_sink=sink)
      IF (PRESENT(norm_squared)) THEN
         norm_squared = sink%norms(1, :) + sink%norms(2, work%mirror_cell)
      END IF
      unit_nr = cp_logger_get_default_io_unit()
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
      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ''
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
         n_pairs = REAL(CEILING(REAL(n_grid, dp)*REAL(n_grid + 1, dp)/(2.0_dp*REAL(i_group_size, dp)), KIND=int_8), dp)
         ! G(k), W(R), the retained off-mesh product, and the prepared FFT input coexist.
         ! Count the band orbitals of the owned band k points alongside these panels.
         ! q matrices, band projection, and Green blocks belong to separate stages.
         grid_bytes = (n_pairs*(complex_bytes + 3.0_dp*real_bytes) + real_bytes)*REAL(n_k, dp) + &
                      REAL(int_1_size, dp)*n_pairs*REAL(n_k, dp) + &
                      MAX(REAL(fft_panel_bytes, dp), 3.0_dp*complex_bytes*REAL(n_k, dp)) + &
                      real_bytes*(n_pairs + REAL(n_k, dp)) + &
                      complex_bytes*REAL(n_grid, dp)*REAL(n_ao, dp)*REAL(n_spin, dp)* &
                      REAL((bs_env%nkp_bs_and_DOS + i_group_size - 1)/i_group_size, dp) + &
                      MAX(complex_bytes*REAL(n_grid, dp)*REAL(n_grid + 1, dp)*REAL(n_local_k, dp)/2.0_dp + &
                          complex_bytes*REAL(n_aux, dp)**2*REAL(n_local_k, dp) + &
                          complex_bytes*REAL(n_threads, dp)*REAL(n_grid, dp)* &
                          REAL(n_ao + n_aux + MIN(n_grid, grid_pair_panel_cols), dp), &
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
!> \param needed q points whose Z(q) is copied; all q if absent.
! **************************************************************************************************
   SUBROUTINE gather_group_columns(bs_env, gwkp_env, column_env, group_size, group_rank, Z_group, psi_group, &
                                   needed)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(IN)                        :: gwkp_env
      TYPE(mp_comm_type), INTENT(IN)                     :: column_env
      INTEGER, INTENT(IN)                                :: group_size, group_rank
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: Z_group(:, :, :)
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT), &
         OPTIONAL                                        :: psi_group(:, :, :, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: needed(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gather_group_columns'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: psi_block(:, :, :), Z_block(:, :)
      INTEGER                                            :: handle, i_root, ik, ik_local, iq, &
                                                            iq_local, n_grid, n_k, n_world
      LOGICAL, ALLOCATABLE                               :: copy_q(:)

      CALL timeset(routineN, handle)

      n_world = bs_env%para_env%num_pe
      n_grid = bs_env%ri_rs%n_grid_points
      n_k = bs_env%nkp_scf_desymm
      ALLOCATE (copy_q(n_k), SOURCE=.TRUE.)
      IF (PRESENT(needed)) copy_q(:) = needed
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
      ALLOCATE (Z_group(n_grid, bs_env%n_RI, COUNT([(MODULO(iq, group_size) == group_rank .AND. copy_q(iq), &
                                                     iq=1, n_k)])))
      ALLOCATE (Z_block(n_grid, bs_env%n_RI))
      iq_local = 0
      DO iq = MODULO(group_rank - 1, group_size) + 1, n_k, group_size
         IF (.NOT. copy_q(iq)) CYCLE
         i_root = MODULO(iq, n_world)/group_size
         IF (column_env%mepos == i_root) Z_block(:, :) = gwkp_env%Z_lP_q(:, :, (iq - 1)/n_world + 1)
         CALL column_env%bcast(Z_block, i_root)
         iq_local = iq_local + 1
         Z_group(:, :, iq_local) = Z_block
      END DO
      DEALLOCATE (copy_q, Z_block)

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
         n_grid, n_groups, n_local, n_needed_local, n_pairs_local, n_q_local, n_threads, &
         n_time_local, n_world, q_batch, unit_nr
      INTEGER(KIND=int_8)                                :: n_pairs_total
      INTEGER, ALLOCATABLE                               :: full_cells(:, :), needed_local(:), &
                                                            needed_q(:), q_partner(:)
      LOGICAL, ALLOCATABLE                               :: cell_active(:), needed(:)
      REAL(KIND=dp) :: available_GB, budget, chi_seconds(7), chi_started, input_bytes, live_bytes, &
         n_pair_max, n_q_max, peak, q_bytes, rows_bytes, stage_started, tau, time_bytes
      REAL(KIND=dp), ALLOCATABLE                         :: green_norms(:, :), qpoints(:, :)
      TYPE(gw_pair_fft_work_type)                        :: fft_stencil_work
      TYPE(mp_comm_type)                                 :: column_env
      TYPE(mp_para_env_type), POINTER                    :: chi_env
      TYPE(ri_q_matrix_block_type), ALLOCATABLE          :: chi_partial(:)

      CALL timeset(routineN, handle)

      CPASSERT(gwkp_env%ngroups == 1 .AND. gwkp_env%my_group == 0)
      CALL build_full_cell_grid(bs_env, full_cells)
      n_grid = bs_env%ri_rs%n_grid_points
      n_pairs_total = INT(n_grid, int_8)*(n_grid + 1)/2
      n_aux = bs_env%n_RI
      n_world = gwkp_env%group_env%num_pe
      ALLOCATE (qpoints(3, bs_env%nkp_scf_desymm))
      qpoints(:, :) = bs_env%kpoints_chi_eps_W%xkp(:, 1:bs_env%nkp_scf_desymm)
      ! χ(-q) = χ(q)*: compute_W0 solves the first q of every ±q pair and conjugates W for -q.
      CALL minus_q_partners(qpoints, q_partner)
      ALLOCATE (needed(bs_env%nkp_scf_desymm))
      needed(:) = [(q_partner(iq) >= iq, iq=1, bs_env%nkp_scf_desymm)]
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
         live_bytes = MAX(rows_bytes + projection_workspace + &
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
         CALL gather_group_columns(bs_env, gwkp_env, column_env, i_size, chi_env%mepos, Z_group, psi_time, needed)
         psi_values => psi_time
         Z_time => Z_group
      ELSE
         psi_values => gwkp_env%psi_lm_k
         Z_time => gwkp_env%Z_lP_q
      END IF
      n_time_local = (bs_env%num_time_freq_points + n_groups - 1 - i_group)/n_groups
      CALL block_range(INT(INT(n_grid, int_8)*(n_grid + 1)/2), chi_env%mepos, chi_env%num_pe, first_pair, n_pairs_local)
      n_q_local = COUNT([(MODULO(iq, chi_env%num_pe) == chi_env%mepos, iq=1, bs_env%nkp_scf_desymm)])
      ALLOCATE (needed_local(n_q_local))
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
            chi_started = m_walltime()
            stage_started = chi_started
            ALLOCATE (G_grid_k(n_pairs_local, 1, bs_env%nkp_scf_desymm))
            CALL build_green_pairs(bs_env, psi_values(:, :, :, ispin), chi_env, tau, ispin, &
                                   .TRUE., fft_stencil_work%k_batch, G_grid_k)
            chi_seconds(1) = m_walltime() - stage_started
            stage_started = m_walltime()
            CALL k_to_grid_pairs(G_grid_k, fft_stencil_work, norm_squared=green_norms(:, 1))
            chi_seconds(2) = m_walltime() - stage_started
            stage_started = m_walltime()
            CALL build_green_pairs(bs_env, psi_values(:, :, :, ispin), chi_env, tau, ispin, &
                                   .FALSE., fft_stencil_work%k_batch, G_grid_k)
            chi_seconds(1) = chi_seconds(1) + m_walltime() - stage_started
            stage_started = m_walltime()
            CALL k_to_grid_pairs(G_grid_k, fft_stencil_work, norm_squared=green_norms(:, 2), &
                                 product_scale=bs_env%spin_degeneracy)
            chi_seconds(2) = chi_seconds(2) + m_walltime() - stage_started
            stage_started = m_walltime()
            ! Each original Green is filtered by its own full Frobenius norm.
            ! Scalar cell masks commute with their product: (m_o G_o)(m_v G_v)=m_o m_v G_o G_v.
            CPASSERT(bs_env%n_atom > 0)
            CALL chi_env%sum(green_norms)
            cell_active(:) = .NOT. (SQRT(green_norms(:, 1))/REAL(bs_env%n_atom, dp)**2 < bs_env%eps_filter) .AND. &
                             .NOT. (SQRT(green_norms(:, 2))/REAL(bs_env%n_atom, dp)**2 < bs_env%eps_filter)
            chi_seconds(3) = m_walltime() - stage_started
            stage_started = m_walltime()
            DEALLOCATE (G_grid_k)
            ALLOCATE (chi_grid_q(n_pairs_local, 1, SIZE(needed_q)))
            CALL execute_pair_fft(fft_stencil_work, chi_grid_q, cell_active, &
                                  'GW chi (K2R extraction feeds R2K input)')
            chi_seconds(4) = m_walltime() - stage_started

            stage_started = m_walltime()
            CALL project_fft_pairs_to_chi(chi_grid_q, chi_env, needed, Z_time, n_groups > 1, &
                                          chi_q, q_batch, chi_seconds(6))
            DEALLOCATE (chi_grid_q)
            chi_seconds(5) = MAX(0.0_dp, m_walltime() - stage_started - chi_seconds(6))
            chi_seconds(7) = m_walltime() - chi_started
            unit_nr = cp_logger_get_default_io_unit()
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
            iq_local = needed_local(i_needed)
            DO i_aux = 1, n_aux
               chi_partial(iq_local)%values(:i_aux, i_aux, i_time_local) = chi_q(:i_aux, i_aux, i_needed)
               chi_partial(iq_local)%values(i_aux + 1:, i_aux, i_time_local) = CONJG(chi_q(i_aux, i_aux + 1:, i_needed))
            END DO
         END DO
      END DO
      CALL k_grid_to_cell_release(fft_stencil_work%inverse)
      bs_env%gwkp_env%k_points_total = 0

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
      DEALLOCATE (chi_q, full_cells, needed, needed_local, needed_q, q_partner, qpoints, green_norms, cell_active)

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
      ! The erfc(omega r)/r metric is the real-space Ewald part of the Coulomb matrix. At q = 0 of an
      ! odd mesh the reciprocal part drops the singular G = 0 term, which the Γ correction integrates.
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
                                                            panel(:, :), q_pairs(:, :), ZK(:, :)
      INTEGER                                            :: first_pair, handle, iq_local, ispin, &
                                                            n_grid, n_pairs_local, n_q_local, n_RI
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
      n_q_local = SIZE(gwkp_env%Z_lP_q, 3)
      IF (SIZE(gwkp_env%Z_lP_zero, 2) /= bs_env%n_RI + 1 .OR. SIZE(gwkp_env%Z_lP_zero, 3) /= 1) THEN
         CPABORT("Periodic RI-RS GW requires the direct Gamma auxiliary fit")
      END IF
      CPASSERT(SIZE(gwkp_env%Z_lP_q, 2) == bs_env%n_RI)
      CPASSERT(ALLOCATED(gwkp_env%metric_coulomb_factor))
      n_RI = bs_env%n_RI
      ALLOCATE (q_pairs(INT(INT(n_grid, int_8)*(n_grid + 1)/2), n_q_local))
      ! V_ll'(q) = Z(q) K(q) K(q)† Z(q)† = (ZK)(ZK)† with the metric-Coulomb factor K.
      !$OMP PARALLEL DEFAULT(NONE) SHARED(gwkp_env, q_pairs, n_q_local, n_grid, n_RI) &
      !$OMP PRIVATE(iq_local, panel, ZK)
      ALLOCATE (ZK(n_grid, n_RI), panel(n_grid, MIN(n_grid, grid_pair_panel_cols)))
      !$OMP DO SCHEDULE(STATIC)
      DO iq_local = 1, n_q_local
         CALL zgemm('N', 'N', n_grid, n_RI, n_RI, z_one, gwkp_env%Z_lP_q(:, :, iq_local), n_grid, &
                    gwkp_env%metric_coulomb_factor(:, :, iq_local), n_RI, z_zero, ZK, n_grid)
         CALL grid_pair_gemm(n_grid, n_RI, ZK, ZK, q_pairs(:, iq_local), panel)
      END DO
      !$OMP END DO
      DEALLOCATE (ZK, panel)
      !$OMP END PARALLEL
      DEALLOCATE (gwkp_env%metric_coulomb_factor)
      CALL build_full_cell_grid(bs_env, full_cells)
      CALL prepare_pair_fft(bs_env, full_cells, first_pair, n_pairs_local, bs_env%para_env%num_pe, &
                            fft_stencil_work, .FALSE.)
      ALLOCATE (G_grid_k(n_pairs_local, 1, bs_env%nkp_scf_desymm))
      CALL q_matrices_to_fft_pairs(q_pairs, n_grid, bs_env%para_env, fft_stencil_work%k_batch, G_grid_k)
      DEALLOCATE (q_pairs)
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
                                0.0_dp, ispin, .TRUE., fft_stencil_work%k_batch, G_grid_k)
         CALL k_to_grid_pairs(G_grid_k, fft_stencil_work, values_R=Sigma_grid_R, multiply_R=V_grid_R, product_scale=-1.0_dp)
         CALL execute_pair_fft(fft_stencil_work, G_grid_k, cell_active, &
                               'GW self-energy (K2R extraction feeds R2K input)', Sigma_grid_R)

         CALL self_energy_pairs_to_Bloch(bs_env, G_grid_k, Sigma_grid_R, first_pair, bs_env%para_env, &
                                         gwkp_env%band_orbitals(:, :, :, ispin), &
                                         gwkp_env%sigma_x_bands(:, :, ispin), .FALSE.)
      END DO
      CALL k_grid_to_cell_release(fft_stencil_work%inverse)

      DEALLOCATE (full_cells, G_grid_k, Sigma_grid_R, V_grid_R, cell_active)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x_from_V

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
                                                            i_time_local, ispin, n_grid, n_pairs_local, &
                                                            n_time_local
      INTEGER, ALLOCATABLE                               :: full_cells(:, :)
      LOGICAL, ALLOCATABLE                               :: cell_active(:)
      REAL(KIND=dp)                                      :: tau
      REAL(KIND=dp), ALLOCATABLE :: green_norms(:), sigma_band(:, :), Sigma_grid_R(:, :, :), &
         sigma_result_buffer(:), W_grid_R(:, :, :)
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
                            fft_stencil_work, .TRUE.)
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
         CALL project_auxiliary_time_q_to_fft_pairs(gwkp_env%Z_lP_q, &
                                                    gwkp_env%W_PQ_tau, i_time_local, &
                                                    gwkp_env%group_env, fft_stencil_work%k_batch, G_grid_k)
         CALL q_mesh_and_gamma_shell_average(bs_env, values_q=G_grid_k, q_order=fft_stencil_work%input_to_native)
         CALL k_to_grid_pairs(G_grid_k, fft_stencil_work, values_R=W_grid_R)
         CALL add_gamma_corr_pairs(bs_env, full_cells, gwkp_env%gamma_correction%head_cell_tau(:, i_t), &
                                   gwkp_env%Z_lP_zero(:, 1, 1), first_pair, W_grid_R, gwkp_env%gamma_correction, i_t)
         tau = bs_env%time_frequency_grid%imaginary_time(i_t)
         DO ispin = 1, bs_env%n_spin
            ! Preserve the component order: occupied (2, -G W), then virtual (1, +G W).
            DO i_component = 2, 1, -1
               CALL build_green_pairs(bs_env, gwkp_env%band_orbitals(:, :, :, ispin), gwkp_env%group_env, &
                                      tau, ispin, i_component == 2, fft_stencil_work%k_batch, G_grid_k)
               CALL k_to_grid_pairs(G_grid_k, fft_stencil_work, values_R=Sigma_grid_R, multiply_R=W_grid_R, &
                                    product_scale=MERGE(-1.0_dp, 1.0_dp, i_component == 2), norm_squared=green_norms)
               CALL gwkp_env%group_env%sum(green_norms)
               cell_active(:) = .NOT. (SQRT(green_norms)/REAL(bs_env%n_atom, dp)**2 < bs_env%eps_filter)
               CALL execute_pair_fft(fft_stencil_work, G_grid_k, cell_active, &
                                     'GW self-energy (K2R extraction feeds R2K input)', Sigma_grid_R)
               CALL self_energy_pairs_to_Bloch(bs_env, G_grid_k, Sigma_grid_R, first_pair, gwkp_env%group_env, &
                                               gwkp_env%band_orbitals(:, :, :, ispin), sigma_band, .TRUE.)
               IF (gwkp_env%group_env%mepos == 0) gwkp_env%sigma_c_bands(:, :, i_t, i_component, ispin) = sigma_band
            END DO
         END DO
      END DO
      CALL k_grid_to_cell_release(fft_stencil_work%inverse)
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
!> \param bs_env ...
!> \param orbitals Band Bloch orbitals of one spin at the k points owned in env, MOD(k-1,P).
!> \param env Ranks that share the grid pairs and the band k points.
!> \param tau ...
!> \param ispin ...
!> \param occupied ...
!> \param k_batch k points per source rank in one exchange
!> \param G_pairs unordered pairs in source-rank-major k order within each exchange
! **************************************************************************************************
   SUBROUTINE build_green_pairs(bs_env, orbitals, env, tau, ispin, occupied, k_batch, G_pairs)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), INTENT(IN)                       :: orbitals(:, :, :)
      TYPE(mp_para_env_type), POINTER                    :: env
      REAL(KIND=dp), INTENT(IN)                          :: tau
      INTEGER, INTENT(IN)                                :: ispin
      LOGICAL, INTENT(IN)                                :: occupied
      INTEGER, INTENT(IN)                                :: k_batch
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(OUT), TARGET  :: G_pairs(:, :, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_green_pairs'

      TYPE pair_segment_type
         INTEGER :: rank, first_row, last_row, displacement
      END TYPE pair_segment_type

      COMPLEX(KIND=dp), ALLOCATABLE                      :: gram(:, :), scaled_coeff(:, :), sent(:)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: received(:)
      COMPLEX(KIND=dp), TARGET                           :: unused_receive(1)
      INTEGER :: first_band, handle, i_band, i_batch, i_col, i_first, i_rank, ikp, j_col, &
                 n_active, n_bands, n_cols, n_grid, n_k, n_local, n_owned, n_pairs, n_rank, n_round, &
                 n_segments, n_workers, offset, pair_first, pair_last, rank_end, row_first, row_last, segment, unit_nr
      INTEGER(KIND=int_8)                                :: active_bands, possible_bands
      INTEGER, ALLOCATABLE                               :: column_segments(:), pair_counts(:), pair_displs(:), &
                                                            recv_counts(:), recv_displs(:), &
                                                            send_counts(:), send_displs(:)
      REAL(KIND=dp)                                      :: seconds(4), started, tau_E, total_started
      TYPE(gwkp_type), POINTER                           :: gwkp_env
      TYPE(pair_segment_type), ALLOCATABLE                :: segments(:)

      CALL timeset(routineN, handle)
      total_started = m_walltime()
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
      CALL grid_pair_blocks(n_grid, n_rank, pair_counts, pair_displs)
      n_pairs = SUM(pair_counts)
      CPASSERT(k_batch > 0 .AND. n_local == pair_counts(env%mepos))
      ALLOCATE (sent(MAX(1, n_pairs*k_batch)), &
                recv_counts(0:n_rank - 1), recv_displs(0:n_rank - 1), &
                send_counts(0:n_rank - 1), send_displs(0:n_rank - 1))
      ! Columns and destination blocks partition the same packed triangle. Only a column
      ! or rank boundary can begin a new segment, so the map is linear in n_grid + n_rank.
      ALLOCATE (column_segments(n_grid + 1), segments(n_grid + n_rank))
      i_rank = 0
      n_segments = 0
      pair_first = 1
      DO j_col = 1, n_grid
         column_segments(j_col) = n_segments + 1
         pair_last = pair_first + j_col - 1
         DO WHILE (i_rank < n_rank)
            rank_end = pair_displs(i_rank) + pair_counts(i_rank)
            row_first = MAX(pair_first, pair_displs(i_rank) + 1) - pair_first + 1
            row_last = MIN(pair_last, rank_end) - pair_first + 1
            IF (row_first <= row_last) THEN
               n_segments = n_segments + 1
               segments(n_segments) = pair_segment_type(i_rank, row_first, row_last, &
                                                        pair_first + (row_first - 1) - pair_displs(i_rank) - 1)
            END IF
            IF (rank_end > pair_last) EXIT
            i_rank = i_rank + 1
         END DO
         IF (j_col < n_grid) pair_first = pair_last + 1
      END DO
      column_segments(n_grid + 1) = n_segments + 1
      unused_receive = z_zero
      seconds(1) = m_walltime() - total_started
      ! G_ll′(k,τ) = Σ_n [ψ_ln(k) exp(-|τ(ε_nk-ε_F)|/2)]
      !                         [ψ_l′n(k) exp(-|τ(ε_nk-ε_F)|/2)]*.
      DO i_first = 1, n_k, n_rank*k_batch
         n_round = MIN(n_rank*k_batch, n_k - i_first + 1)
         n_owned = MAX(0, (n_round + n_rank - 1 - env%mepos)/n_rank)
         send_counts(:) = pair_counts*n_owned
         send_displs(:) = pair_displs*n_owned
         n_workers = 1
!$       n_workers = MIN(omp_get_max_threads(), MAX(1, n_owned))
         started = m_walltime()
         !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(n_workers) &
         !$OMP SHARED(bs_env, env, first_band, gwkp_env, i_first, n_bands, n_grid, n_owned, n_rank) &
         !$OMP SHARED(orbitals, pair_counts, send_displs, sent, tau, ispin, column_segments, segments) &
         !$OMP PRIVATE(scaled_coeff, gram, i_band, i_batch, ikp, n_active, tau_E) &
         !$OMP PRIVATE(i_col, n_cols, j_col, i_rank, row_first, row_last, offset, segment) &
         !$OMP REDUCTION(+:active_bands, possible_bands)
         ALLOCATE (scaled_coeff(n_grid, n_bands), gram(n_grid, column_block_cols))
         !$OMP DO SCHEDULE(STATIC)
         DO i_batch = 1, n_owned
            ikp = i_first + env%mepos + (i_batch - 1)*n_rank
            n_active = 0
            DO i_band = 1, n_bands
               tau_E = ABS(0.5_dp*tau*(gwkp_env%eigenvalues_G(first_band + i_band - 1, ikp, ispin) - &
                                       bs_env%e_fermi(ispin)))
               IF (tau_E < bs_env%stabilize_exp) THEN
                  n_active = n_active + 1
                  scaled_coeff(:, n_active) = &
                     orbitals(:, first_band + i_band - 1, (ikp - 1)/n_rank + 1)*EXP(-tau_E)
               END IF
            END DO
            active_bands = active_bands + INT(n_active, int_8)
            possible_bands = possible_bands + INT(n_bands, int_8)
            DO i_col = 1, n_grid, column_block_cols
               n_cols = MIN(column_block_cols, n_grid - i_col + 1)
               IF (n_active == 0) THEN
                  gram = z_zero
               ELSE
                  CALL zherk('U', 'N', n_cols, n_active, 1.0_dp, scaled_coeff(i_col, 1), n_grid, &
                             0.0_dp, gram(i_col, 1), n_grid)
                  IF (i_col > 1) THEN
                     CALL zgemm('N', 'C', i_col - 1, n_cols, n_active, z_one, scaled_coeff, n_grid, &
                                scaled_coeff(i_col, 1), n_grid, z_zero, gram, n_grid)
                  END IF
               END IF
               DO j_col = i_col, i_col + n_cols - 1
                  DO segment = column_segments(j_col), column_segments(j_col + 1) - 1
                     i_rank = segments(segment)%rank
                     row_first = segments(segment)%first_row
                     row_last = segments(segment)%last_row
                     offset = send_displs(i_rank) + (i_batch - 1)*pair_counts(i_rank) + segments(segment)%displacement
                     sent(offset + 1:offset + (row_last - row_first + 1)) = gram(row_first:row_last, j_col - i_col + 1)
                  END DO
               END DO
            END DO
         END DO
         !$OMP END DO
         DEALLOCATE (scaled_coeff, gram)
         !$OMP END PARALLEL
         ! Scaling and Gram construction overlap between k jobs; time their shared wall interval.
         seconds(2) = seconds(2) + m_walltime() - started
         offset = 0
         DO i_rank = 0, n_rank - 1
            recv_counts(i_rank) = n_local*MAX(0, (n_round + n_rank - 1 - i_rank)/n_rank)
            recv_displs(i_rank) = offset
            offset = offset + recv_counts(i_rank)
         END DO
         CPASSERT(offset == n_local*n_round)
         received => unused_receive
         IF (n_local > 0) received(1:n_local*n_round) => G_pairs(:, :, i_first:i_first + n_round - 1)
         started = m_walltime()
         CALL env%alltoall(sent, send_counts, send_displs, received, recv_counts, recv_displs)
         seconds(3) = seconds(3) + m_walltime() - started
         CALL count_time_loop_k_points(bs_env, INT(n_round, int_8))
      END DO
      seconds(4) = m_walltime() - total_started
      unit_nr = cp_logger_get_default_io_unit()
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A)') 'GW_GREEN_PROFILE| logger-rank wall seconds, no MPI reduction'
         WRITE (unit_nr, '(T2,A,L2,2I12)') 'GW_GREEN_PROFILE| occupied / k points / local pairs: ', occupied, n_k, n_local
         WRITE (unit_nr, '(T2,A,2I16)') 'GW_GREEN_PROFILE| active / possible local band columns: ', active_bands, possible_bands
         WRITE (unit_nr, '(T2,A,I12)') 'GW_GREEN_PROFILE| k points per rank / exchange: ', k_batch
         WRITE (unit_nr, '(T2,A,4F12.6)') 'GW_GREEN_PROFILE| setup/build/transfer/total: ', seconds
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_green_pairs

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
!> \brief Project Σ_ll′(k) onto the band Bloch orbitals, Σ_n(k) = Re[ψ_n(k)† Σ(k) ψ_n(k)].
!>
!> The pair owners evaluate Σ_ll′(k) = Σ_R exp(2πik·R) Σ_ll′(R) at the band k points. Every k
!> owner, MOD(k-1,P), receives complete packed upper triangles for Hermitian contractions. A
!> round covers m k points per rank, so that one pass over Σ(R) serves up to
!> band_kpoints_per_pass k points.
!> \param bs_env ...
!> \param sigma_mesh Σ on the native k mesh for the unordered grid pairs of this rank.
!> \param sigma_R Σ in the cells of the native mesh for the same pairs.
!> \param first_pair First locally owned unordered grid pair.
!> \param env Ranks that share the grid pairs and the band k points.
!> \param orbitals Band Bloch orbitals of one spin at the k points owned in env.
!> \param sigma_band ...
!> \param apply_cutoff use the W distance cutoff in off-mesh image sums
! **************************************************************************************************
   SUBROUTINE self_energy_pairs_to_Bloch(bs_env, sigma_mesh, sigma_R, first_pair, env, orbitals, sigma_band, apply_cutoff)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(IN), TARGET   :: sigma_mesh(:, :, :)
      REAL(KIND=dp), INTENT(IN)                          :: sigma_R(:, :, :)
      INTEGER, INTENT(IN)                                :: first_pair
      TYPE(mp_para_env_type), POINTER                    :: env
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(IN)           :: orbitals(:, :, :)
      REAL(KIND=dp), INTENT(OUT)                         :: sigma_band(:, :)
      LOGICAL, INTENT(IN)                                :: apply_cutoff

      CHARACTER(LEN=*), PARAMETER :: routineN = 'self_energy_pairs_to_Bloch'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: panel(:, :), sigma_k(:, :), work(:, :)
      COMPLEX(KIND=dp), ALLOCATABLE, TARGET              :: send_buffer(:)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: sent(:), sigma_columns(:, :)
      COMPLEX(KIND=dp), TARGET                           :: unused_send(1)
      INTEGER                                            :: handle, i_band, i_batch, i_col, i_first, &
                                                            i_kp, i_rank, j, m, n_batch, n_grid, &
                                                            n_kp, n_local, n_rank, n_round, n_send
      INTEGER, ALLOCATABLE                               :: pair_counts(:), pair_displs(:), &
                                                            recv_counts(:), send_counts(:), &
                                                            send_displs(:)
      LOGICAL                                            :: native_round

      CALL timeset(routineN, handle)

      n_grid = bs_env%ri_rs%n_grid_points
      n_kp = bs_env%nkp_bs_and_DOS
      n_local = SIZE(sigma_mesh, 1)
      n_rank = env%num_pe
      n_batch = MIN(16, bs_env%n_ao)
      m = projection_round_kpoints(n_grid, n_rank)
      CALL grid_pair_blocks(n_grid, n_rank, pair_counts, pair_displs)
      ALLOCATE (sigma_k(INT(INT(n_grid, int_8)*(n_grid + 1)/2), m), recv_counts(0:n_rank - 1), &
                send_counts(0:n_rank - 1), send_displs(0:n_rank - 1))
      unused_send = z_zero
      send_displs(:) = [(i_rank*n_local, i_rank=0, n_rank - 1)]
      sigma_band = 0.0_dp
      ! Column i_col of a round holds the band k point i_first+i_col-1, owned by MOD(i_col-1,P).
      DO i_first = 1, n_kp, n_rank*m
         n_round = MIN(n_rank*m, n_kp - i_first + 1)
         native_round = i_first + n_round - 1 <= SIZE(sigma_mesh, 3)
         IF (native_round) THEN
            DO i_col = 1, n_round
               IF (native_q_index(bs_env, i_first + i_col - 1) /= i_first + i_col - 1) THEN
                  native_round = .FALSE.
                  EXIT
               END IF
            END DO
         END IF
         IF (.NOT. native_round) THEN
            IF (.NOT. ALLOCATED(send_buffer)) THEN
               ALLOCATE (send_buffer(n_local*n_rank*m))
               sigma_columns(1:n_local, 1:n_rank*m) => send_buffer
            END IF
            CALL band_self_energy(bs_env, sigma_mesh, sigma_R, first_pair, i_first, &
                                  sigma_columns(:, 1:n_round), apply_cutoff)
         END IF
         DO j = 1, m
            DO i_rank = 0, n_rank - 1
               send_counts(i_rank) = MERGE(n_local, 0, (j - 1)*n_rank + i_rank < n_round)
            END DO
            i_kp = i_first + (j - 1)*n_rank + env%mepos
            recv_counts(:) = MERGE(pair_counts, 0, i_kp <= n_kp)
            sent => unused_send
            n_send = MIN(n_rank, MAX(0, n_round - (j - 1)*n_rank))
            IF (n_local > 0 .AND. n_send > 0) THEN
               IF (native_round) THEN
                  sent(1:n_local*n_send) => sigma_mesh(:, :, &
                                                       i_first + (j - 1)*n_rank:i_first + (j - 1)*n_rank + n_send - 1)
               ELSE
                  sent => send_buffer((j - 1)*n_rank*n_local + 1:((j - 1)*n_rank + n_send)*n_local)
               END IF
            END IF
            CALL env%alltoall(sent, send_counts, send_displs, sigma_k(:, j), recv_counts, pair_displs)
         END DO
         !$OMP PARALLEL DEFAULT(NONE) &
         !$OMP SHARED(bs_env, env, i_first, m, n_batch, n_grid, n_kp, n_rank, orbitals, sigma_band, sigma_k) &
         !$OMP PRIVATE(i_band, i_batch, i_kp, j, panel, work)
         ALLOCATE (work(n_grid, n_batch), panel(n_grid, MIN(n_grid, grid_pair_panel_cols)))
         !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
         DO j = 1, m
            DO i_batch = 1, bs_env%n_ao, n_batch
               i_kp = i_first + (j - 1)*n_rank + env%mepos
               IF (i_kp > n_kp) CYCLE
               CALL grid_pair_upper_multiply(n_grid, MIN(n_batch, bs_env%n_ao - i_batch + 1), sigma_k(:, j), &
                                             orbitals(:, i_batch:MIN(i_batch + n_batch - 1, bs_env%n_ao), (i_kp - 1)/n_rank + 1), &
                                             work(:, 1:MIN(n_batch, bs_env%n_ao - i_batch + 1)), panel)
               DO i_band = i_batch, MIN(i_batch + n_batch - 1, bs_env%n_ao)
                  sigma_band(i_band, i_kp) = 2.0_dp*REAL(SUM(CONJG(orbitals(:, i_band, (i_kp - 1)/n_rank + 1))* &
                                                             work(:, i_band - i_batch + 1)), dp)
               END DO
            END DO
         END DO
         !$OMP END DO
         DEALLOCATE (work, panel)
         !$OMP END PARALLEL
         CALL count_time_loop_k_points(bs_env, INT(n_round, int_8))
      END DO
      CALL env%sum(sigma_band)

      CALL timestop(handle)

   END SUBROUTINE self_energy_pairs_to_Bloch

! **************************************************************************************************
!> \brief Bound the band-projection and Green-function workspaces of one rank in bytes.
!> The projection holds the Σ(k) columns of one round, the complete Σ(k) of the owned k points,
!> the shared lattice phases, thread-private image, sum and output blocks, and the 64 MiB copy of
!> the mesh transform. The Green function holds bounded k-batch messages, received pair blocks,
!> thread-private Gram panels and the scaled orbitals.
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
      bytes = complex_bytes*m*REAL(n_grid, dp)*REAL(n_grid + 1, dp) + complex_bytes*n_round*REAL(n_images, dp) + &
              projection_workspace + &
              REAL(n_threads, dp)*(complex_bytes*REAL(n_grid, dp)*REAL(MIN(16, n_ao) + MIN(n_grid, grid_pair_panel_cols), dp) + &
                                   real_bytes*REAL(pair_batch_size, dp)*REAL(n_images, dp) + &
                                   complex_bytes*REAL(pair_batch_size, dp)*n_round)
      bytes = MAX(bytes, projection_workspace + &
                  complex_bytes*REAL(n_threads, dp)*REAL(n_grid, dp)*REAL(n_ao + 32, dp))
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

      m = MAX(1, MIN(band_kpoints_per_pass/n_rank, &
                     INT(projection_round_bytes/(complex_bytes*REAL(n_grid, dp)*REAL(n_grid + 1, dp)))))
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
!> \param apply_cutoff use the W distance cutoff on the same centered images as the native FFT
! **************************************************************************************************
   SUBROUTINE band_self_energy(bs_env, sigma_mesh, sigma_R, first_pair, first_kp_eval, sigma_k, apply_cutoff)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), INTENT(IN)                       :: sigma_mesh(:, :, :)
      REAL(KIND=dp), INTENT(IN)                          :: sigma_R(:, :, :)
      INTEGER, INTENT(IN)                                :: first_pair, first_kp_eval
      COMPLEX(KIND=dp), INTENT(OUT)                      :: sigma_k(:, :)
      LOGICAL, INTENT(IN)                                :: apply_cutoff

      COMPLEX(KIND=dp)                                   :: phase
      INTEGER :: c1, c2, c3, first_c(3), i_c, i_col, i_first, i_local, i_off, i_pair, il, ilp, &
         last_c(3), lower(3), n_c, n_grid, n_off, n_pairs, native_iq, nmesh(3), shift(3), &
         wrapped(3)
      INTEGER, ALLOCATABLE                               :: offmesh_col(:), pair_shift(:, :)
      LOGICAL                                            :: any_tie(3), tie(3)
      REAL(KIND=dp)                                      :: angle, boundary(3), cutoff_squared, &
                                                            displacement(3), distance(3)
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
         CALL grid_pair_indices(first_pair + i_local - 1, il, ilp)
         CALL grid_pair_shift(il, ilp, fractional_grid, shift, tie)
         any_tie = any_tie .OR. tie
      END DO
      boundary = 0.25_dp*(1.0_dp + COS(twopi*bs_env%kpoints_scf_desymm%xkp(:, 1)*nmesh))
      cutoff_squared = -1.0_dp
      IF (apply_cutoff .AND. bs_env%ri_rs%cutoff_radius_v_w > 0.0_dp) THEN
         cutoff_squared = bs_env%ri_rs%cutoff_radius_v_w**2
      END IF
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
      !$OMP SHARED(bs_env, cutoff_squared, boundary, cos_c, first_c, first_pair, fractional_grid, last_c, lower, n_c, n_grid) &
      !$OMP SHARED(n_off, nmesh, offmesh_col, sigma_k, sigma_R, sin_c, xk) &
      !$OMP PRIVATE(angle, c1, c2, c3, i_c, i_first, i_off, i_pair, images, n_pairs, pair_shift, shift, &
      !$OMP         sum_im, sum_re, tie, weight_1, weight_2, weight_3, wrapped, il, ilp, displacement, distance)
      ALLOCATE (images(n_c, pair_batch_size), pair_shift(3, pair_batch_size), &
                sum_re(pair_batch_size, n_off), sum_im(pair_batch_size, n_off), &
                weight_1(first_c(1):last_c(1)), weight_2(first_c(2):last_c(2)), weight_3(first_c(3):last_c(3)))
      !$OMP DO SCHEDULE(DYNAMIC)
      DO i_first = 1, SIZE(sigma_R, 1), pair_batch_size
         n_pairs = MIN(pair_batch_size, SIZE(sigma_R, 1) - i_first + 1)
         DO i_pair = 1, n_pairs
            CALL grid_pair_indices(first_pair + i_first + i_pair - 2, il, ilp)
            CALL grid_pair_shift(il, ilp, fractional_grid, shift, tie)
            displacement = bs_env%ri_rs%grid_points(:, il) - bs_env%ri_rs%grid_points(:, ilp)
            pair_shift(:, i_pair) = shift
            weight_1(:) = [(image_weight(c1, nmesh(1), tie(1), boundary(1)), c1=first_c(1), last_c(1))]
            weight_2(:) = [(image_weight(c2, nmesh(2), tie(2), boundary(2)), c2=first_c(2), last_c(2))]
            weight_3(:) = [(image_weight(c3, nmesh(3), tie(3), boundary(3)), c3=first_c(3), last_c(3))]
            i_c = 0
            DO c3 = first_c(3), last_c(3)
               DO c2 = first_c(2), last_c(2)
                  DO c1 = first_c(1), last_c(1)
                     i_c = i_c + 1
                     images(i_c, i_pair) = 0.0_dp
                     IF (cutoff_squared > 0.0_dp) THEN
                        distance = displacement - MATMUL(bs_env%hmat, REAL(shift + [c1, c2, c3], dp))
                        IF (SUM(distance**2) > cutoff_squared) CYCLE
                     END IF
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
! **************************************************************************************************
   SUBROUTINE prepare_pair_fft(bs_env, cells, first_pair, n_pairs, n_rank, plan, apply_cutoff, q_indices)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: cells(:, :), first_pair, n_pairs, n_rank
      TYPE(gw_pair_fft_work_type), INTENT(INOUT), TARGET :: plan
      LOGICAL, INTENT(IN)                                :: apply_cutoff
      INTEGER, INTENT(IN), OPTIONAL                      :: q_indices(:)

      INTEGER :: count, fft_count, first, half_shift(3), i, i_extract, i_first, i_rank, ibatch, &
         icell, idir, ik, il, ilp, iq, isource, lower(3), n_extract, n_k, n_round, ngrid, &
         nmesh(3), partner_output, reflected(3), shift(3), width
      INTEGER(KIND=int_8)                                :: nfft_grid, spectrum_columns
      INTEGER, ALLOCATABLE :: extract_partner(:), extract_q(:), grid_linear(:), &
         input_to_native(:), native_to_input(:), partner_linear(:), q_position(:)
      LOGICAL                                            :: tie(3)
      REAL(KIND=dp)                                      :: cutoff_squared, displacement(3), &
                                                            started, total_started
      REAL(KIND=dp), ALLOCATABLE                         :: fractional_grid(:, :)
      TYPE(gw_pair_fft_batch_type), POINTER              :: batch

      total_started = m_walltime()
      CPASSERT(.NOT. ALLOCATED(plan%batches))
      ngrid = bs_env%ri_rs%n_grid_points
      nmesh = bs_env%kpoints_scf_desymm%nkp_grid
      lower = -nmesh/2
      ! GW uses the same native reciprocal mesh for Green functions, chi and W.
      IF (.NOT. ALL(bs_env%kpoints_chi_eps_W%xkp(:, 1:bs_env%nkp_scf_desymm) == bs_env%kpoints_scf_desymm%xkp)) THEN
         CPABORT('GW lattice transforms require one native reciprocal mesh')
      END IF
      n_k = bs_env%nkp_scf_desymm
      plan%k_batch = MAX(1, MIN((n_k + n_rank - 1)/n_rank, &
                                INT(projection_workspace/(complex_bytes* &
                                                          (REAL(ngrid, dp)*REAL(ngrid + 1, dp)/2.0_dp + REAL(n_rank, dp))))))
      ALLOCATE (input_to_native(n_k), native_to_input(n_k))
      i_extract = 0
      ! Receive contiguous source-rank blocks; the inverse FFT owns the resulting k permutation.
      DO i_first = 1, n_k, n_rank*plan%k_batch
         n_round = MIN(n_rank*plan%k_batch, n_k - i_first + 1)
         DO i_rank = 0, n_rank - 1
            DO ik = i_first + i_rank, i_first + n_round - 1, n_rank
               i_extract = i_extract + 1
               input_to_native(i_extract) = ik
               native_to_input(ik) = i_extract
            END DO
         END DO
      END DO
      CPASSERT(i_extract == n_k)
      ALLOCATE (plan%product_phase(SIZE(cells, 2)))
      CALL k_grid_to_cell_prepare(plan%inverse, bs_env%kpoints_scf_desymm%xkp(:, input_to_native), nmesh, cells, &
                                  cell_phase_forward=plan%product_phase)
      CPASSERT(ALL(plan%inverse%nfft == nmesh([3, 2, 1])))
      CPASSERT(ALLOCATED(plan%inverse%k_partner))
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
      DO i_extract = 1, n_extract
         ik = extract_q(i_extract)
         spectrum_columns = spectrum_columns + MERGE(1_int_8, 2_int_8, grid_linear(ik) == partner_linear(ik))
      END DO

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
         started = m_walltime()
         CALL fft_alloc(batch%input, [INT(nfft_grid)*fft_count])
         plan%profile%seconds(3) = plan%profile%seconds(3) + m_walltime() - started
         started = m_walltime()
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
         plan%setup_seconds = plan%setup_seconds + m_walltime() - started
      END DO
      CALL MOVE_ALLOC(grid_linear, plan%grid_linear)
      CALL MOVE_ALLOC(partner_linear, plan%partner_linear)
      CALL MOVE_ALLOC(input_to_native, plan%input_to_native)
      ALLOCATE (plan%extract_q, SOURCE=extract_q(:n_extract))
      ALLOCATE (plan%extract_partner, SOURCE=extract_partner(:n_extract))
      plan%profile%seconds(2) = m_walltime() - total_started - plan%profile%seconds(3)
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
                                                            nmesh(3), partner_output, stat, &
                                                            unit_nr, width
      INTEGER(KIND=int_8)                                :: extract_reads, fft_transforms, &
                                                            nfft_grid, spectrum_columns
      INTEGER, POINTER                                   :: extract_partner(:), extract_q(:), &
                                                            grid_linear(:), partner_linear(:)
      REAL(KIND=dp)                                      :: elapsed, other_seconds, &
                                                            prepared_seconds, started, &
                                                            total_started
      TYPE(fft_pair_major_plan_type)                     :: full_fft_plan, tail_fft_plan
      TYPE(lattice_fft_profile_type)                     :: profile

      total_started = m_walltime()
      CPASSERT(ALLOCATED(plan%batches))
      CPASSERT(SIZE(cell_active) == SIZE(plan%inverse%cell_slot))
      CPASSERT(SIZE(values_k, 1) == plan%n_pairs)
      CPASSERT(SIZE(values_k, 2) == 1)
      CPASSERT(SIZE(values_k, 3) == SIZE(plan%grid_linear))
      profile = plan%profile
      prepared_seconds = SUM(profile%seconds(1:3))
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
      started = m_walltime()
      CALL fft_alloc(scratch, [INT(nfft_grid)*((width + 1)/2)])
      profile%seconds(3) = profile%seconds(3) + m_walltime() - started
      DO first = 1, plan%n_pairs, width
         count = MIN(width, plan%n_pairs - first + 1)
         fft_count = (count + 1)/2
         fft_transforms = fft_transforms + INT(fft_count, int_8)
         ibatch = (first - 1)/width + 1
         profile%elements = profile%elements + INT(count, int_8)
         profile%batches = profile%batches + 1_int_8
         profile%tile_min = MIN(profile%tile_min, count)
         profile%tile_max = MAX(profile%tile_max, count)
         fft_in => plan%batches(ibatch)%input
         fft_out => scratch
         input(1:fft_count, 1:INT(nfft_grid)) => fft_in
         started = m_walltime()
         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
         !$OMP SHARED(cell_active, input, plan, first, count, values_R) PRIVATE(icell)
         DO icell = 1, SIZE(cell_active)
            IF (cell_active(icell)) CYCLE
            input(:, plan%inverse%cell_slot(icell)) = z_zero
            IF (PRESENT(values_R)) values_R(first:first + count - 1, 1, icell) = 0.0_dp
         END DO
         !$OMP END PARALLEL DO
         profile%seconds(4) = profile%seconds(4) + m_walltime() - started
         profile%pair_major = .TRUE.
         IF (count == width) THEN
            CALL lattice_fft_execute_pair_major(BWFFT, nfft, fft_count, fft_in, fft_out, stat, profile, full_fft_plan)
         ELSE
            CALL lattice_fft_execute_pair_major(BWFFT, nfft, fft_count, fft_in, fft_out, stat, profile, tail_fft_plan)
         END IF
         IF (stat /= 0) CPABORT('GW product FFT failed')
         output(1:fft_count, 1:INT(nfft_grid)) => fft_in
         started = m_walltime()
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
                  values_k(first + count - 1, 1, partner_output) = 2.0_dp*output(fft_count, partner_linear(ik))
               END IF
            END IF
         END DO
         !$OMP END PARALLEL DO
         profile%seconds(8) = profile%seconds(8) + m_walltime() - started
         extract_reads = extract_reads + INT(count/2, int_8)*spectrum_columns + &
                         INT(MODULO(count, 2), int_8)*SIZE(grid_linear, KIND=int_8)
      END DO

      started = m_walltime()
      CALL fft_pair_major_plan_destroy(full_fft_plan)
      CALL fft_pair_major_plan_destroy(tail_fft_plan)
      CALL fft_dealloc(scratch)
      profile%seconds(9) = m_walltime() - started
      elapsed = m_walltime() - total_started + prepared_seconds
      other_seconds = MAX(0.0_dp, elapsed - profile%seconds(4) - SUM(profile%seconds(5:8)))
      unit_nr = cp_logger_get_default_io_unit()
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
!> \brief Transfer q-local unordered pairs directly to the FFT pair owners.
!> Ranks own q points cyclically, MOD(q,P), and packed grid pairs as contiguous blocks. Round r
!> moves q = (r-1)P+1, ..., rP, the r-th local q of every rank, with one alltoall.
!> Receive displacements place each column directly in the inverse FFT's rank-blocked k batch.
!> \param q_pairs packed upper triangles
!> \param n_grid ...
!> \param para_env ...
!> \param k_batch source-rank block width of the prepared inverse FFT input
!> \param fft_pairs unordered pairs in the prepared inverse FFT input order
! **************************************************************************************************
   SUBROUTINE q_matrices_to_fft_pairs(q_pairs, n_grid, para_env, k_batch, fft_pairs)
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(IN), TARGET   :: q_pairs(:, :)
      INTEGER, INTENT(IN)                                :: n_grid
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: k_batch
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(OUT), TARGET  :: fft_pairs(:, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'q_matrices_to_fft_pairs'

      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: received(:), sent(:)
      COMPLEX(KIND=dp), TARGET                           :: unused_receive(1), unused_send(1)
      INTEGER :: batch_first, batch_round, batch_size, handle, i_first, i_rank, i_round, i_source, &
         n_local, n_owned, n_q, n_rank, n_round, offset
      INTEGER, ALLOCATABLE                               :: pair_counts(:), pair_displs(:), &
                                                            recv_counts(:), recv_displs(:), &
                                                            send_counts(:)

      CALL timeset(routineN, handle)

      n_rank = para_env%num_pe
      n_local = SIZE(fft_pairs, 1)
      n_q = SIZE(fft_pairs, 3)
      CPASSERT(SIZE(fft_pairs, 2) == 1 .AND. k_batch > 0)
      CALL grid_pair_blocks(n_grid, n_rank, pair_counts, pair_displs)
      CPASSERT(SIZE(q_pairs, 1) == SUM(pair_counts))
      ALLOCATE (recv_counts(0:n_rank - 1), &
                recv_displs(0:n_rank - 1), send_counts(0:n_rank - 1))
      unused_send = z_zero
      unused_receive = z_zero
      DO i_round = 1, (n_q + n_rank - 1)/n_rank
         i_first = (i_round - 1)*n_rank
         n_round = MIN(n_rank, n_q - i_first)
         batch_first = ((i_round - 1)/k_batch)*n_rank*k_batch
         batch_round = MODULO(i_round - 1, k_batch)
         batch_size = MIN(n_rank*k_batch, n_q - batch_first)
         offset = 0
         DO i_rank = 0, n_rank - 1
            ! q ownership is rotated by one rank relative to Green k ownership.
            i_source = MODULO(i_rank + 1, n_rank)
            n_owned = MAX(0, (batch_size + n_rank - 1 - i_rank)/n_rank)
            recv_counts(i_source) = MERGE(n_local, 0, batch_round < n_owned)
            recv_displs(i_source) = n_local*(offset + MIN(batch_round, n_owned))
            offset = offset + n_owned
         END DO
         send_counts(:) = 0
         sent => unused_send
         IF (MODULO(para_env%mepos - 1, n_rank) < n_round) THEN
            sent => q_pairs(:, i_round)
            send_counts(:) = pair_counts
         END IF
         received => unused_receive
         IF (n_local > 0) THEN
            received(1:n_local*batch_size) => fft_pairs(:, :, batch_first + 1:batch_first + batch_size)
         END IF
         CALL para_env%alltoall(sent, send_counts, pair_displs, received, recv_counts, recv_displs)
      END DO

      CALL timestop(handle)

   END SUBROUTINE q_matrices_to_fft_pairs

! **************************************************************************************************
!> \brief Transfer unordered grid pairs into bounded q batches and consume the Hermitian
!>        panels as Z(q)^H chi(q) Z(q), without expanding a full grid matrix.
!> \param fft_pairs compact needed-q columns, ordered by ascending original q index
!> \param para_env ...
!> \param needed q points to transfer
!> \param Z_q interpolation matrices at the q points owned by this rank, MOD(q,P)
!> \param compact_z whether Z_q contains only needed q points
!> \param chi_q accumulate projected matrices in compact local q order
!> \param batch_rounds common maximum number of communication rounds per projection batch
!> \param projection_seconds logger-rank projection wall time, excluding transfers
! **************************************************************************************************
   SUBROUTINE project_fft_pairs_to_chi(fft_pairs, para_env, needed, Z_q, compact_z, chi_q, batch_rounds, projection_seconds)
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(IN), TARGET   :: fft_pairs(:, :, :)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: needed(:)
      COMPLEX(KIND=dp), INTENT(IN)                       :: Z_q(:, :, :)
      LOGICAL, INTENT(IN)                                :: compact_z
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: chi_q(:, :, :)
      INTEGER, INTENT(IN)                                :: batch_rounds
      REAL(KIND=dp), INTENT(OUT)                         :: projection_seconds

      CHARACTER(LEN=*), PARAMETER :: routineN = 'project_fft_pairs_to_chi'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: panel(:, :, :), work(:, :, :)
      COMPLEX(KIND=dp), ALLOCATABLE, TARGET              :: q_pairs(:, :)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: received(:), sent(:)
      COMPLEX(KIND=dp), TARGET                           :: unused_receive(1), unused_send(1)
      INTEGER :: capacity, first_round, handle, i_buffer, i_column, i_first, i_local, i_needed, i_rank, &
         i_round, i_thread, i_z, iq, n_aux, n_buffer, n_grid, n_local, n_needed_local, n_q, &
         n_rank, n_round, n_selected, n_workers
      INTEGER, ALLOCATABLE                               :: local_q(:), pair_counts(:), &
                                                            pair_displs(:), recv_counts(:), &
                                                            send_counts(:), send_displs(:)
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
                send_counts(0:n_rank - 1), send_displs(0:n_rank - 1), local_q(capacity), &
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
               local_q(n_buffer) = (iq - 1)/n_rank + 1
               received => q_pairs(:, n_buffer)
            END IF
            CALL para_env%alltoall(sent, send_counts, send_displs, received, recv_counts, pair_displs)
         END DO
         IF (n_buffer == 0) CYCLE
         started = m_walltime()
         ! Distinct q matrices remain independent; each thread reuses its GEMM workspace.
         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) IF(n_buffer > 1) NUM_THREADS(n_workers) &
         !$OMP SHARED(n_buffer, n_grid, n_aux, q_pairs, local_q, Z_q, chi_q, panel, work, i_local, compact_z) &
         !$OMP PRIVATE(i_buffer, i_thread, i_needed, i_z)
         DO i_buffer = 1, n_buffer
            i_thread = 1
!$          i_thread = omp_get_thread_num() + 1
            i_needed = i_local - n_buffer + i_buffer
            i_z = MERGE(i_needed, local_q(i_buffer), compact_z)
            CALL grid_pair_upper_multiply(n_grid, n_aux, q_pairs(:, i_buffer), Z_q(:, :, i_z), &
                                          work(:, :, i_thread), panel(:, :, i_thread))
            ! Z^H H Z = Z^H (U Z) + (U Z)^H Z; accumulate only the Hermitian upper triangle.
            CALL zher2k('U', 'C', n_aux, n_grid, z_one, Z_q(:, :, i_z), n_grid, &
                        work(:, :, i_thread), n_grid, 1.0_dp, chi_q(:, :, i_needed), n_aux)
         END DO
         !$OMP END PARALLEL DO
         projection_seconds = projection_seconds + m_walltime() - started
      END DO
      CPASSERT(i_column == SIZE(fft_pairs, 3))
      CPASSERT(i_local == n_needed_local)

      CALL timestop(handle)

   END SUBROUTINE project_fft_pairs_to_chi

! **************************************************************************************************
!> \brief Form W_ll'(q,iτ) = Σ_PQ Z_lP(q) W_PQ(q,iτ) Z_l'Q(q)* and transfer FFT pairs.
!> \param Z_q ...
!> \param interaction_q W_PQ(q,iτ) of the locally owned q
!> \param itime ...
!> \param para_env ...
!> \param k_batch source-rank block width of the prepared inverse FFT input
!> \param fft_pairs ...
! **************************************************************************************************
   SUBROUTINE project_auxiliary_time_q_to_fft_pairs(Z_q, interaction_q, itime, para_env, k_batch, fft_pairs)
      COMPLEX(KIND=dp), INTENT(IN)                       :: Z_q(:, :, :)
      TYPE(ri_q_matrix_block_type), INTENT(IN)           :: interaction_q(:)
      INTEGER, INTENT(IN)                                :: itime
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: k_batch
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(OUT), TARGET  :: fft_pairs(:, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'project_auxiliary_time_q_to_fft_pairs'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: panel(:, :), q_pairs(:, :), ZW(:, :)
      INTEGER                                            :: handle, i_q_local, n_grid, n_q_local, &
                                                            n_RI

      CALL timeset(routineN, handle)

      n_grid = SIZE(Z_q, 1)
      n_RI = SIZE(Z_q, 2)
      n_q_local = SIZE(Z_q, 3)
      CPASSERT(SIZE(interaction_q) == n_q_local)
      ALLOCATE (q_pairs(INT(INT(n_grid, int_8)*(n_grid + 1)/2), n_q_local))
      !$OMP PARALLEL DEFAULT(NONE) SHARED(Z_q, interaction_q, itime, q_pairs, n_q_local, n_grid, n_RI) &
      !$OMP PRIVATE(i_q_local, panel, ZW)
      ALLOCATE (ZW(n_grid, n_RI), panel(n_grid, MIN(n_grid, grid_pair_panel_cols)))
      !$OMP DO SCHEDULE(STATIC)
      DO i_q_local = 1, n_q_local
         CALL zgemm('N', 'N', n_grid, n_RI, n_RI, z_one, Z_q(:, :, i_q_local), n_grid, &
                    interaction_q(i_q_local)%values(:, :, itime), n_RI, z_zero, ZW, n_grid)
         CALL grid_pair_gemm(n_grid, n_RI, ZW, Z_q(:, :, i_q_local), q_pairs(:, i_q_local), panel)
      END DO
      !$OMP END DO
      DEALLOCATE (ZW, panel)
      !$OMP END PARALLEL
      CALL q_matrices_to_fft_pairs(q_pairs, n_grid, para_env, k_batch, fft_pairs)

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

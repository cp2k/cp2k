!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Gamma-point (q = 0) correction of W and Σ for periodic RI-RS GW.
!> \author Jan Wilhelm (with help of ChatGPT, GPT-6 Sol and GPT-6 Astra)
! **************************************************************************************************
MODULE gw_ri_rs_small_cell_full_kp_gamma
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_transpose
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_set_submatrix,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_get_info,&
                                              dbcsr_get_readonly_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_to_string
   USE gw_auto_ri_types,                ONLY: auto_ri_type
   USE gw_ri_rs_small_cell_full_kp_lattice,ONLY: assemble_AA_contraction,&
                                              block_range,&
                                              build_full_cell_grid,&
                                              build_ri_2c_lattice,&
                                              gamma_auxiliary_charges,&
                                              owned_mo_coefficients,&
                                              ri_2c_lattice_to_kp,&
                                              ri_metric_head_zero
   USE gw_ri_rs_small_cell_full_kp_types,ONLY: gamma_atom_projection_type,&
                                              gamma_correction_type,&
                                              gwkp_type
   USE gw_utils_fm,                     ONLY: local_complex_power
   USE kinds,                           ONLY: dp
   USE kpoint_coulomb_2c,               ONLY: add_ewald_reciprocal,&
                                              build_2c_coulomb_matrix_kp_dense
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index
   USE kpoint_transform,                ONLY: kp_transform_plan_create,&
                                              kp_transform_plan_type
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_create,&
                                              kpoint_release,&
                                              kpoint_type
   USE mathconstants,                   ONLY: fourpi,&
                                              gaussi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE mathlib,                         ONLY: gauss_legendre,&
                                              inv_3x3,&
                                              vector_product
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_moments,                      ONLY: build_local_moment_matrix_rs_img
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_small_cell_full_kp_gamma'

   PUBLIC :: add_gamma_corr_cell_correction, &
             add_gamma_corr_pair_cell, &
             gamma_point_correction, &
             prepare_ri_rs_gamma_correction, &
             singular_head_cell_factors

   ! Below |a*L| = taylor_switch, radial phase moments use their Taylor series.
   REAL(KIND=dp), PARAMETER :: taylor_switch = 1.0E-4_dp

CONTAINS

! **************************************************************************************************
!> \brief Preserve distributed AO position, Hamiltonian and overlap images before QS release.
!> Only the sparse operators and their immutable image traversal survive; SCF grids and
!> neighbor lists can be released by the common GW initialization path.
!> \param qs_env Converged Quickstep state.
!> \param bs_env Band-structure settings.
! **************************************************************************************************
   SUBROUTINE prepare_ri_rs_gamma_correction(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_ri_rs_gamma_correction'

      INTEGER                                            :: handle, i_image, i_spin, n_img
      LOGICAL                                            :: created_scf_image_map
      TYPE(dbcsr_p_type), POINTER                        :: matrix_h(:, :), matrix_s(:, :)
      TYPE(gwkp_type), POINTER                           :: gwkp_env
      TYPE(kpoint_type), POINTER                         :: kpoints_all, kpoints_scf
      TYPE(neighbor_list_set_p_type), POINTER            :: sab_all(:), sab_kp(:)

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ALLOCATED(bs_env%gwkp_env))
      ALLOCATE (bs_env%gwkp_env)
      gwkp_env => bs_env%gwkp_env
      ALLOCATE (gwkp_env%eigenvalues_G(bs_env%n_ao, bs_env%nkp_scf_desymm, bs_env%n_spin))
      gwkp_env%eigenvalues_G(:, :, :) = bs_env%eigenval_scf(:, 1:bs_env%nkp_scf_desymm, :)

      NULLIFY (kpoints_all)

      CALL get_qs_env(qs_env, kpoints=kpoints_scf, matrix_ks_kp=matrix_h, matrix_s_kp=matrix_s, &
                      sab_all=sab_all, sab_kp=sab_kp)
      CPASSERT(ASSOCIATED(sab_all) .AND. ASSOCIATED(sab_kp))
      CPASSERT(ASSOCIATED(matrix_h) .AND. ASSOCIATED(matrix_s))
      CPASSERT(.NOT. ASSOCIATED(gwkp_env%gamma_corr_d))
      created_scf_image_map = .NOT. ASSOCIATED(kpoints_scf)
      IF (created_scf_image_map) THEN
         CALL kpoint_create(kpoints_scf)
         CALL kpoint_init_cell_index(kpoints_scf, sab_kp, bs_env%para_env, n_img)
      END IF
      CALL kpoint_create(kpoints_all)
      CALL kpoint_init_cell_index(kpoints_all, sab_all, bs_env%para_env, n_img)
      CALL gw_gamma_cell_index_union(kpoints_all, kpoints_scf, &
                                     gwkp_env%gamma_corr_index_to_cell, gwkp_env%gamma_corr_cell_to_index)
      CALL dbcsr_allocate_matrix_set(gwkp_env%gamma_corr_d, 3, SIZE(kpoints_all%index_to_cell, 2))
      CALL build_local_moment_matrix_rs_img(qs_env, gwkp_env%gamma_corr_d)
      CALL kp_transform_plan_create(gwkp_env%gamma_corr_plan, kpoints_all%sab_nl, &
                                    kpoints_all%cell_to_index, SIZE(gwkp_env%gamma_corr_d, 2))
      CALL kp_transform_plan_create(gwkp_env%gamma_corr_scf_plan, kpoints_scf%sab_nl, &
                                    kpoints_scf%cell_to_index, SIZE(matrix_h, 2))
      CALL dbcsr_allocate_matrix_set(gwkp_env%gamma_corr_h, SIZE(matrix_h, 1), SIZE(matrix_h, 2))
      CALL dbcsr_allocate_matrix_set(gwkp_env%gamma_corr_s, SIZE(matrix_s, 1), SIZE(matrix_s, 2))
      DO i_image = 1, SIZE(matrix_h, 2)
         DO i_spin = 1, SIZE(matrix_h, 1)
            ALLOCATE (gwkp_env%gamma_corr_h(i_spin, i_image)%matrix)
            CALL dbcsr_copy(gwkp_env%gamma_corr_h(i_spin, i_image)%matrix, matrix_h(i_spin, i_image)%matrix)
         END DO
      END DO
      DO i_image = 1, SIZE(matrix_s, 2)
         DO i_spin = 1, SIZE(matrix_s, 1)
            ALLOCATE (gwkp_env%gamma_corr_s(i_spin, i_image)%matrix)
            CALL dbcsr_copy(gwkp_env%gamma_corr_s(i_spin, i_image)%matrix, matrix_s(i_spin, i_image)%matrix)
         END DO
      END DO
      CALL kpoint_release(kpoints_all)
      IF (created_scf_image_map) CALL kpoint_release(kpoints_scf)

      CALL timestop(handle)

   END SUBROUTINE prepare_ri_rs_gamma_correction

! **************************************************************************************************
!> \brief Gamma-correction head and wing coefficients for the near-Gamma cells.
!> Construct the q->0 dielectric limit along face-cubature and mesh-sample directions
!> from SCF dipole transitions and the directly fitted Gamma Gaussian body.
!> \param bs_env ...
!> \param qs_env ...
!> \param gwkp_env RI-RS work data and Gamma correction.
! **************************************************************************************************
   SUBROUTINE gamma_point_correction(bs_env, qs_env, gwkp_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gwkp_type), INTENT(INOUT)                     :: gwkp_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_point_correction'

      COMPLEX(KIND=dp)                                   :: eps_00, pi_00, pi_q(3), &
                                                            q_direction_complex(3), &
                                                            schur_contribution, schur_denom
      COMPLEX(KIND=dp), ALLOCATABLE :: dipole(:, :, :, :, :), eps_P_cartesian(:, :), &
         F_external(:, :), F_response_cartesian(:, :), F_zero(:, :), metric_ref_zero(:, :), &
         metric_rows(:, :), metric_zero(:, :), pi_cartesian(:, :), response_cartesian(:, :), &
         schur_cartesian(:, :), time_chi(:, :, :), time_head(:, :, :), time_wing(:, :, :), &
         V_sqrt_zero(:, :)
      INTEGER                                            :: a, first_row, handle, i_cartesian, &
                                                            i_direction, ikp, j_cartesian, j_w, &
                                                            n_angular, n_directions, n_near, &
                                                            n_w_kpoints
      INTEGER, ALLOCATABLE                               :: metric_cells(:, :)
      LOGICAL                                            :: is_slab
      REAL(KIND=dp)                                      :: near_q(3, 8), q_cart(3), q_direction(3), &
                                                            reciprocal(3, 3), volume
      REAL(KIND=dp), ALLOCATABLE :: angular_directions(:, :), angular_radii(:), &
         angular_weights(:), directions(:, :), head_response(:, :), radii(:), weights(:)
      TYPE(cp_fm_type), ALLOCATABLE                      :: fm_metric(:)

      CALL timeset(routineN, handle)

      is_slab = ALL(bs_env%periodic == [1, 1, 0])
      volume = bs_env%ri_rs%cell%deth
      n_w_kpoints = PRODUCT(bs_env%nkp_grid_chi_eps_W_orig)
      reciprocal = twopi*TRANSPOSE(inv_3x3(bs_env%hmat))
      IF (is_slab) THEN
         CALL slab_reciprocal_vectors(bs_env, reciprocal(:, 1:2), volume)
         CALL slab_cell_face_directions(bs_env, angular_directions, angular_weights, angular_radii)
      ELSE
         CALL gamma_cell_face_directions(bs_env, angular_directions, angular_weights, angular_radii)
      END IF
      n_angular = SIZE(angular_weights)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T72,I10)') &
            'Number of angular points for the Γ-point correction of W:', n_angular
      END IF
      n_near = PRODUCT(MIN(2, bs_env%nkp_grid_chi_eps_W_orig))
      n_directions = n_angular + n_near
      a = 0
      DO ikp = 1, n_w_kpoints
         q_direction = bs_env%kpoints_chi_eps_W%xkp(1:3, ikp) - &
                       ANINT(bs_env%kpoints_chi_eps_W%xkp(1:3, ikp))
         IF (is_slab) THEN
            IF (.NOT. is_gamma_neighbor(q_direction(1:2), bs_env%nkp_grid_chi_eps_W_orig(1:2))) CYCLE
         ELSE
            IF (.NOT. is_gamma_neighbor(q_direction, bs_env%nkp_grid_chi_eps_W_orig)) CYCLE
         END IF
         a = a + 1
         IF (a > n_near) CPABORT('Too many Gamma-neighbor transfer points')
         near_q(:, a) = q_direction
      END DO
      IF (a /= n_near) CPABORT('Inconsistent Gamma-neighbor transfer points')
      ALLOCATE (directions(3, n_directions), weights(n_directions), radii(n_directions))
      directions(:, 1:n_angular) = angular_directions
      weights(1:n_angular) = angular_weights
      radii(1:n_angular) = angular_radii
      DEALLOCATE (angular_directions, angular_weights, angular_radii)
      DO i_direction = 1, n_near
         q_cart = MATMUL(reciprocal, near_q(:, i_direction))
         radii(n_angular + i_direction) = NORM2(q_cart)
         directions(:, n_angular + i_direction) = q_cart/radii(n_angular + i_direction)
         weights(n_angular + i_direction) = 0.0_dp
      END DO
      ALLOCATE (head_response(n_directions, bs_env%num_time_freq_points))
      ! The finite RI metric at q=0 retains its complete Gaussian channel.
      ! The singular full-Coulomb Gamma-correction head is integrated separately.
      ALLOCATE (metric_zero(bs_env%n_RI, bs_env%n_RI), F_zero(bs_env%n_RI, bs_env%n_RI))
      CALL build_ri_2c_lattice(fm_metric, metric_cells, bs_env%ri_metric, &
                               0.0_dp, bs_env, qs_env, reference_basis=bs_env%auto_ri%enabled)
      IF (bs_env%auto_ri%enabled) THEN
         CALL gamma_metric_contract_atoms(fm_metric, bs_env%auto_ri, metric_zero)
         IF (.NOT. is_slab) THEN
            ALLOCATE (metric_ref_zero(SUM(bs_env%auto_ri%sizes_ref_RI), SUM(bs_env%auto_ri%sizes_ref_RI)))
            CALL ri_2c_lattice_to_kp(fm_metric, metric_cells, [0.0_dp, 0.0_dp, 0.0_dp], metric_ref_zero)
         END IF
      ELSE
         CALL ri_2c_lattice_to_kp(fm_metric, metric_cells, [0.0_dp, 0.0_dp, 0.0_dp], metric_zero)
         IF (.NOT. is_slab) ALLOCATE (metric_ref_zero, SOURCE=metric_zero)
      END IF
      CALL cp_fm_release(fm_metric)
      DEALLOCATE (metric_cells)
      IF (is_slab) THEN
         CALL slab_gamma_factors(qs_env, bs_env, metric_zero, F_zero, F_external)
      ELSE
         CALL coulomb_matrix_gamma(qs_env, bs_env, V_sqrt_zero, metric_ref_zero)
         DEALLOCATE (metric_ref_zero)
         CALL gamma_metric_factor_distributed(bs_env%para_env, bs_env%regularization_RI, &
                                              metric_zero, V_sqrt_zero, F_zero)
         DEALLOCATE (V_sqrt_zero)
      END IF
      DEALLOCATE (metric_zero)
      ! In 3D the external-potential map equals the density map.
      IF (.NOT. is_slab) ALLOCATE (F_external, SOURCE=F_zero)

      ! Gamma-correction matrix elements in the actual SCF MO gauge.  q_x cancels analytically against
      ! the 1/q and 1/q^2 factors, so q_x=1 is used after the limit has been taken.
      CALL gw_gamma_transition_moments(bs_env, gwkp_env, bs_env%kpoints_scf_desymm, &
                                       bs_env%cfm_mo_coeff_kp(1:bs_env%nkp_scf_desymm, :), gwkp_env%eigenvalues_G, &
                                       bs_env%n_occ(1:bs_env%n_spin), dipole)

      ALLOCATE (eps_P_cartesian(bs_env%n_RI, 3), pi_cartesian(3, 3), &
                response_cartesian(bs_env%n_RI, 3), schur_cartesian(3, 3), &
                F_response_cartesian(bs_env%n_RI, 3))
      ALLOCATE (gwkp_env%gamma_correction%wing_cartesian(bs_env%ri_rs%n_grid_points, 3, &
                                                         bs_env%num_time_freq_points))
      CALL accumulate_time_state(bs_env, gwkp_env%psi_lm_k, dipole, volume, first_row, &
                                 time_chi, time_wing, time_head)
      CALL metric_grid_rows(gwkp_env%Z_lP_zero(:, 2:, 1), F_zero, first_row, SIZE(time_chi, 1), metric_rows)
      DO j_w = 1, bs_env%num_time_freq_points
         CALL gamma_ao_time_response(bs_env, time_chi, time_wing, time_head, metric_rows, first_row, &
                                     j_w, eps_P_cartesian, pi_cartesian, response_cartesian)
         CALL zgemm('C', 'N', 3, 3, bs_env%n_RI, CMPLX(fourpi, 0.0_dp, dp), &
                    eps_P_cartesian, bs_env%n_RI, response_cartesian, bs_env%n_RI, &
                    z_zero, schur_cartesian, 3)
         ! The 2D Coulomb kernel is 2π/q, the 3D one 4π/q².
         CALL zgemm('N', 'N', bs_env%n_RI, 3, bs_env%n_RI, z_one, &
                    F_external, bs_env%n_RI, response_cartesian, bs_env%n_RI, &
                    z_zero, F_response_cartesian, bs_env%n_RI)
         CALL zgemm('N', 'N', bs_env%ri_rs%n_grid_points, 3, bs_env%n_RI, &
                    CMPLX(MERGE(twopi, fourpi, is_slab)/ri_metric_head_zero(bs_env), 0.0_dp, dp), &
                    gwkp_env%Z_lP_zero(:, 2:, 1), bs_env%ri_rs%n_grid_points, F_response_cartesian, bs_env%n_RI, &
                    z_zero, gwkp_env%gamma_correction%wing_cartesian(:, :, j_w), bs_env%ri_rs%n_grid_points)
         DO i_direction = 1, n_directions
            q_direction_complex(:) = CMPLX(directions(:, i_direction), 0.0_dp, dp)
            pi_q(:) = MATMUL(pi_cartesian, q_direction_complex)
            pi_00 = DOT_PRODUCT(q_direction_complex, pi_q)
            eps_00 = z_one + fourpi*pi_00
            schur_contribution = z_zero
            DO j_cartesian = 1, 3
               DO i_cartesian = 1, 3
                  schur_contribution = schur_contribution + &
                                       CONJG(q_direction_complex(i_cartesian))* &
                                       schur_cartesian(i_cartesian, j_cartesian)*q_direction_complex(j_cartesian)
               END DO
            END DO
            IF (is_slab) THEN
               ! A(d,ω) = d†[α - E† B⁻¹ E]d is the screened 2D polarizability.
               head_response(i_direction, j_w) = MAX(0.0_dp, REAL(pi_00 - schur_contribution/fourpi, dp))
               CYCLE
            END IF
            schur_denom = eps_00 - schur_contribution
            IF (REAL(schur_denom, dp) <= 0.0_dp) THEN
               CPABORT("Gamma dielectric Schur complement is not positive")
            END IF
            ! The Schur complement is a Hermitian quadratic form in the real d, so ε⁻¹_00(d) is real.
            head_response(i_direction, j_w) = REAL(z_one/schur_denom, dp)
         END DO

      END DO

      CALL gamma_corr_cell_integrals(bs_env, directions, weights, radii, n_angular, head_response, &
                                     gwkp_env%gamma_correction)
      DEALLOCATE (directions, head_response, radii, weights)
      DEALLOCATE (time_chi, time_wing, time_head, metric_rows)
      DEALLOCATE (F_external)

      DEALLOCATE (dipole, &
                  eps_P_cartesian, F_zero, &
                  pi_cartesian, response_cartesian, schur_cartesian, &
                  F_response_cartesian)

      CALL timestop(handle)

   END SUBROUTINE gamma_point_correction

! **************************************************************************************************
!> \brief Coulomb matrix V_PQ(q=0) in the (AUTO_RI contracted) RI basis.
!> In 3D, V(0) is the erfc RI metric M(0) plus the reciprocal Ewald sum, whose q=0 term
!> -pi/(Omega omega^2) F(0) F(0)^H leaves V without its singular G=0 Fourier coefficient.
!> The 2D Ewald sum has no q=0 limit; slabs sum the lattice directly and project the RI charge
!> out afterwards.
!> \param qs_env ...
!> \param bs_env ...
!> \param V_zero ...
!> \param metric_ref_zero M(0) in the reference RI basis, required in 3D
! **************************************************************************************************
   SUBROUTINE coulomb_matrix_gamma(qs_env, bs_env, V_zero, metric_ref_zero)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: V_zero(:, :)
      COMPLEX(KIND=dp), INTENT(IN), OPTIONAL             :: metric_ref_zero(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'coulomb_matrix_gamma'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: U(:, :), V_contracted(:, :), &
                                                            V_owned(:, :, :)
      INTEGER                                            :: handle, n_local, n_ref
      TYPE(kpoint_type), POINTER                         :: gamma_kp

      CALL timeset(routineN, handle)

      n_ref = bs_env%n_RI
      IF (bs_env%auto_ri%enabled) n_ref = SUM(bs_env%auto_ri%sizes_ref_RI)
      n_local = 0
      IF (MODULO(1, bs_env%para_env%num_pe) == bs_env%para_env%mepos) n_local = 1
      ALLOCATE (gamma_kp, V_owned(n_ref, n_ref, n_local), V_zero(n_ref, n_ref))
      ALLOCATE (gamma_kp%xkp(3, 1))
      gamma_kp%xkp = 0.0_dp
      gamma_kp%nkp = 1
      gamma_kp%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig
      V_owned = z_zero
      IF (ALL(bs_env%periodic == [1, 1, 0])) THEN
         CALL build_2c_coulomb_matrix_kp_dense(V_owned, qs_env, gamma_kp, &
                                               bs_env%size_lattice_sum_V, 'RI_AUX', 1, 1)
      ELSE
         CPASSERT(PRESENT(metric_ref_zero))
         IF (n_local == 1) V_owned(:, :, 1) = metric_ref_zero
         CALL add_ewald_reciprocal(V_owned, qs_env, gamma_kp, 'RI_AUX', bs_env%ri_metric%omega, 1, 1)
      END IF
      V_zero(:, :) = z_zero
      IF (n_local == 1) V_zero(:, :) = V_owned(:, :, 1)
      CALL bs_env%para_env%sum(V_zero)
      DEALLOCATE (gamma_kp%xkp)
      DEALLOCATE (gamma_kp, V_owned)
      IF (bs_env%auto_ri%enabled) THEN
         CALL assemble_AA_contraction(bs_env, U)
         ALLOCATE (V_contracted(bs_env%n_RI, bs_env%n_RI))
         V_contracted(:, :) = MATMUL(CONJG(TRANSPOSE(U)), MATMUL(V_zero, U))
         CALL MOVE_ALLOC(V_contracted, V_zero)
         DEALLOCATE (U)
      END IF

      CALL timestop(handle)

   END SUBROUTINE coulomb_matrix_gamma

! **************************************************************************************************
!> \brief Angular cubature over the cells adjacent to Gamma.
!> With mesh dimensions n_i, the fractional half-widths are h_i=min(1/n_i,1/2).
!> The face Jacobian is Omega_BZ*PRODUCT(h_i)/r**3.  The region contains
!> PRODUCT(min(2,n_i)) distinct mesh cells; a one-point direction covers its
!> entire Brillouin-zone interval exactly once, not two periodic copies.
!> \param bs_env ...
!> \param directions ...
!> \param weights ...
!> \param radii ...
! **************************************************************************************************
   SUBROUTINE gamma_cell_face_directions(bs_env, directions, weights, radii)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: directions(:, :), weights(:), radii(:)

      INTEGER, PARAMETER                                 :: n = 6

      INTEGER                                            :: face, i, ia, ib, j, k, n_angular, nk, &
                                                            sign_index
      REAL(KIND=dp)                                      :: bmat(3, 3), jacobian_scale, q(3), &
                                                            qfrac(3), volume
      REAL(KIND=dp), ALLOCATABLE                         :: nodes(:), rule_weights(:)

      ! Six faces with a 6-by-6 Gauss-Legendre rule give 216 angular points.
      n_angular = 6*n*n
      nk = PRODUCT(bs_env%nkp_grid_chi_eps_W_orig)
      jacobian_scale = REAL(PRODUCT(MIN(2, bs_env%nkp_grid_chi_eps_W_orig)), dp)/8.0_dp
      volume = bs_env%ri_rs%cell%deth
      bmat = twopi*TRANSPOSE(inv_3x3(bs_env%hmat))
      ALLOCATE (nodes(n), rule_weights(n))
      CALL gauss_legendre(n, nodes, rule_weights)
      ALLOCATE (directions(3, n_angular), weights(n_angular), radii(n_angular))
      face = 0
      DO i = 1, 3
         j = MODULO(i, 3) + 1
         k = MODULO(i + 1, 3) + 1
         DO sign_index = -1, 1, 2
            DO ia = 1, n
               DO ib = 1, n
                  face = face + 1
                  qfrac(i) = REAL(sign_index, dp)/REAL(MAX(2, bs_env%nkp_grid_chi_eps_W_orig(i)), dp)
                  qfrac(j) = nodes(ia)/REAL(MAX(2, bs_env%nkp_grid_chi_eps_W_orig(j)), dp)
                  qfrac(k) = nodes(ib)/REAL(MAX(2, bs_env%nkp_grid_chi_eps_W_orig(k)), dp)
                  q = MATMUL(bmat, qfrac)
                  radii(face) = NORM2(q)
                  directions(:, face) = q/radii(face)
                  weights(face) = (twopi**3/volume)/(fourpi*REAL(nk, dp)*radii(face)**3)* &
                                  rule_weights(ia)*rule_weights(ib)*jacobian_scale
               END DO
            END DO
         END DO
      END DO
      CPASSERT(face == n_angular)
      DEALLOCATE (nodes, rule_weights)
   END SUBROUTINE gamma_cell_face_directions

! **************************************************************************************************
!> \brief Whether a point of the half-shifted q mesh is one of the transfer points next to Gamma.
!> On a half-shifted mesh, 2*|q_i|*N_i is an odd integer; it equals one next to Gamma.
!> \param xkp fractional q point
!> \param nkp_grid q mesh along the same directions
!> \return ...
! **************************************************************************************************
   PURE FUNCTION is_gamma_neighbor(xkp, nkp_grid) RESULT(neighbor)
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:)
      INTEGER, INTENT(IN)                                :: nkp_grid(:)
      LOGICAL                                            :: neighbor

      neighbor = ALL(NINT(2.0_dp*ABS(xkp - ANINT(xkp))*REAL(nkp_grid, dp)) == 1)

   END FUNCTION is_gamma_neighbor

! **************************************************************************************************
!> \brief Radial moments f0 = int_0^L exp(-i a q) dq and f1 = int_0^L q exp(-i a q) dq.
!> A Taylor series replaces the closed forms where they cancel for small |a L|.
!> \param a ...
!> \param L ...
!> \param f0 ...
!> \param f1 ...
! **************************************************************************************************
   PURE SUBROUTINE radial_phase_moments(a, L, f0, f1)
      REAL(KIND=dp), INTENT(IN)                          :: a, L
      COMPLEX(KIND=dp), INTENT(OUT)                      :: f0
      COMPLEX(KIND=dp), INTENT(OUT), OPTIONAL            :: f1

      COMPLEX(KIND=dp)                                   :: phase

      IF (ABS(a*L) < taylor_switch) THEN
         f0 = CMPLX(L - a*a*L**3/6.0_dp + a**4*L**5/120.0_dp, &
                    -a*L*L/2.0_dp + a**3*L**4/24.0_dp, dp)
         IF (PRESENT(f1)) f1 = CMPLX(L*L/2.0_dp - a*a*L**4/8.0_dp, &
                                     -a*L**3/3.0_dp + a**3*L**5/30.0_dp, dp)
      ELSE
         phase = EXP(-gaussi*a*L)
         f0 = (z_one - phase)/(gaussi*a)
         IF (PRESENT(f1)) f1 = (phase*(z_one + gaussi*a*L) - z_one)/(a*a)
      END IF

   END SUBROUTINE radial_phase_moments

! **************************************************************************************************
!> \brief Integrate 1/q**2 (3D) or 1/|q| (2D slab) over the distinct shifted-mesh cells around Gamma.
!> The regular mesh sum already contains the near-Gamma model samples, which are
!> subtracted here.  In 3D the returned coefficient multiplies the factored head
!> density 4*pi*Z_l0(0)*conjg(Z_l'0(0))/M_0(0)**2; in 2D it has units of inverse
!> in-plane momentum.  The face-to-polar Jacobian resolves the integrable singularity
!> without sampling q=0.
!> \param bs_env ...
!> \param cells ...
!> \param factors ...
! **************************************************************************************************
   SUBROUTINE singular_head_cell_factors(bs_env, cells, factors)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: cells(:, :)
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: factors(:)

      COMPLEX(KIND=dp)                                   :: f0
      INTEGER                                            :: icell, idir, iq, n_dim, nk
      REAL(KIND=dp)                                      :: area, bmat(3, 3), q(3), qnorm, rcart(3), &
                                                            scale
      REAL(KIND=dp), ALLOCATABLE                         :: directions(:, :), radii(:), weights(:)

      nk = PRODUCT(bs_env%nkp_grid_chi_eps_W_orig)
      IF (ALL(bs_env%periodic == [1, 1, 0])) THEN
         n_dim = 2
         CALL slab_cell_face_directions(bs_env, directions, weights, radii)
         bmat(:, 3) = 0.0_dp
         CALL slab_reciprocal_vectors(bs_env, bmat(:, 1:2), area)
         scale = area/twopi**2
      ELSE
         n_dim = 3
         CALL gamma_cell_face_directions(bs_env, directions, weights, radii)
         bmat = twopi*TRANSPOSE(inv_3x3(bs_env%hmat))
         scale = bs_env%ri_rs%cell%deth*fourpi/(twopi**3)
      END IF
      ALLOCATE (factors(SIZE(cells, 2)))
      DO icell = 1, SIZE(cells, 2)
         rcart = MATMUL(bs_env%hmat, REAL(cells(:, icell), dp))
         factors(icell) = z_zero
         DO idir = 1, SIZE(weights)
            CALL radial_phase_moments(DOT_PRODUCT(directions(:, idir), rcart), radii(idir), f0)
            factors(icell) = factors(icell) + weights(idir)*f0
         END DO
         factors(icell) = factors(icell)*scale
         DO iq = 1, nk
            IF (.NOT. is_gamma_neighbor(bs_env%kpoints_chi_eps_W%xkp(1:n_dim, iq), &
                                        bs_env%nkp_grid_chi_eps_W_orig(1:n_dim))) CYCLE
            q = MATMUL(bmat, bs_env%kpoints_chi_eps_W%xkp(1:3, iq) - &
                       ANINT(bs_env%kpoints_chi_eps_W%xkp(1:3, iq)))
            qnorm = NORM2(q)
            factors(icell) = factors(icell) - EXP(-gaussi*DOT_PRODUCT(q, rcart))/ &
                             (REAL(nk, dp)*qnorm**(n_dim - 1))
         END DO
      END DO
      DEALLOCATE (directions, weights, radii)
   END SUBROUTINE singular_head_cell_factors

! **************************************************************************************************
!> \brief Separate the 2D charge mode and factor the neutral Γ Coulomb body.
!> For charge vector f, P=I-ff†/(f†f) and V_b=P V_reg(0) P. The Gamma correction
!> density map is M⁻¹ P V_b¹ᐟ²; the external-potential map is M⁻¹ V_reg P V_b⁻¹ᐟ².
!> Keeping the second map distinct retains the finite charge/body couplings.
!> \param qs_env ...
!> \param bs_env ...
!> \param metric ...
!> \param F_density ...
!> \param F_external ...
! **************************************************************************************************
   SUBROUTINE slab_gamma_factors(qs_env, bs_env, metric, F_density, F_external)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: metric(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: F_density(:, :)
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: F_external(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'slab_gamma_factors'

      COMPLEX(KIND=dp), ALLOCATABLE :: body_inverse(:, :), body_sqrt(:, :), charge(:), &
         charge_contracted(:), projector(:, :), U(:, :), V_regular(:, :)
      INTEGER                                            :: handle, i_col, i_row, n_ref, n_ri

      CALL timeset(routineN, handle)

      n_ri = bs_env%n_RI
      n_ref = n_ri
      IF (bs_env%auto_ri%enabled) n_ref = SUM(bs_env%auto_ri%sizes_ref_RI)
      CALL coulomb_matrix_gamma(qs_env, bs_env, V_regular)
      ALLOCATE (charge(n_ref))
      CALL gamma_auxiliary_charges(qs_env, charge)
      IF (bs_env%auto_ri%enabled) THEN
         CALL assemble_AA_contraction(bs_env, U)
         ALLOCATE (charge_contracted(n_ri))
         charge_contracted(:) = MATMUL(CONJG(TRANSPOSE(U)), charge)
         CALL MOVE_ALLOC(charge_contracted, charge)
         DEALLOCATE (U)
      END IF
      CPASSERT(REAL(DOT_PRODUCT(charge, charge), dp) > 0.0_dp)
      charge(:) = charge/SQRT(REAL(DOT_PRODUCT(charge, charge), dp))
      ALLOCATE (projector(n_ri, n_ri), body_sqrt(n_ri, n_ri), body_inverse(n_ri, n_ri))
      DO i_col = 1, n_ri
         DO i_row = 1, n_ri
            projector(i_row, i_col) = -charge(i_row)*CONJG(charge(i_col))
         END DO
         projector(i_col, i_col) = projector(i_col, i_col) + z_one
      END DO
      body_sqrt(:, :) = MATMUL(projector, MATMUL(V_regular, projector))
      body_inverse(:, :) = body_sqrt
      CALL local_complex_power(body_sqrt, 0.5_dp, 0.0_dp, relative_eps=bs_env%regularization_RI)
      CALL local_complex_power(body_inverse, -0.5_dp, 0.0_dp, relative_eps=bs_env%regularization_RI)
      CALL local_complex_power(metric, -1.0_dp, 0.0_dp, relative_eps=bs_env%regularization_RI)
      F_density = MATMUL(metric, MATMUL(projector, body_sqrt))
      ALLOCATE (F_external(n_ri, n_ri))
      F_external(:, :) = MATMUL(metric, MATMUL(V_regular, MATMUL(projector, body_inverse)))
      DEALLOCATE (V_regular, charge, projector, body_sqrt, body_inverse)

      CALL timestop(handle)

   END SUBROUTINE slab_gamma_factors

! **************************************************************************************************
!> \brief In-plane reciprocal vectors and area for an arbitrary XY cell.
!> \param bs_env ...
!> \param b ...
!> \param area ...
! **************************************************************************************************
   SUBROUTINE slab_reciprocal_vectors(bs_env, b, area)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(OUT)                         :: b(3, 2), area

      REAL(KIND=dp)                                      :: normal(3)

      normal = vector_product(bs_env%hmat(:, 1), bs_env%hmat(:, 2))
      area = NORM2(normal)
      CPASSERT(area > 0.0_dp)
      normal = normal/area
      b(:, 1) = twopi*vector_product(bs_env%hmat(:, 2), normal)/area
      b(:, 2) = twopi*vector_product(normal, bs_env%hmat(:, 1))/area
   END SUBROUTINE slab_reciprocal_vectors

! **************************************************************************************************
!> \brief Polar boundary quadrature of the four shifted Γ-neighbor cells.
!> Directions d, radii L and weights approximate dθ on the bounding parallelogram.
!> As in 3D, the fractional half-widths are min(1/n_i,1/2).
!> \param bs_env ...
!> \param directions ...
!> \param weights ...
!> \param radii ...
! **************************************************************************************************
   SUBROUTINE slab_cell_face_directions(bs_env, directions, weights, radii)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: directions(:, :), weights(:), radii(:)

      INTEGER, PARAMETER                                 :: n_face = 24

      INTEGER                                            :: i_axis, i_direction, i_node, i_sign, &
                                                            j_axis
      REAL(KIND=dp)                                      :: area, b(3, 2), dq(3), q(3)
      REAL(KIND=dp), ALLOCATABLE                         :: nodes(:), rule(:)

      CALL slab_reciprocal_vectors(bs_env, b, area)
      ALLOCATE (nodes(n_face), rule(n_face))
      CALL gauss_legendre(n_face, nodes, rule)
      ALLOCATE (directions(3, 4*n_face), weights(4*n_face), radii(4*n_face))
      i_direction = 0
      DO i_axis = 1, 2
         j_axis = 3 - i_axis
         dq = b(:, j_axis)/REAL(MAX(2, bs_env%nkp_grid_chi_eps_W_orig(j_axis)), dp)
         DO i_sign = -1, 1, 2
            DO i_node = 1, n_face
               i_direction = i_direction + 1
               q = REAL(i_sign, dp)*b(:, i_axis)/REAL(MAX(2, bs_env%nkp_grid_chi_eps_W_orig(i_axis)), dp) &
                   + nodes(i_node)*dq
               radii(i_direction) = NORM2(q)
               directions(:, i_direction) = q/radii(i_direction)
               weights(i_direction) = rule(i_node)*NORM2(vector_product(q, dq))/SUM(q*q)
            END DO
         END DO
      END DO
      DEALLOCATE (nodes, rule)
   END SUBROUTINE slab_cell_face_directions

! **************************************************************************************************
!> \brief Integrate the Gamma-correction head and wings over the near-Gamma cells, in 2D and 3D.
!> The regular mesh sum already contains the near-Gamma mesh samples (directions after
!> n_angular), which are subtracted here; exp(-iqR) is retained per real-space cell R.
!> Per cell only the head and one Cartesian wing coefficient per frequency are stored; the wing
!> on the grid is wing_cartesian(ω) times that coefficient.
!> 3D: ε⁻¹_00(d,ω) does not depend on |q|, so Wc_00 = 4π(ε⁻¹_00-1)/q² and Wc_l0 = -4π ε⁻¹_00 t_l·d/q
!>     have time-independent radial moments.
!> 2D: with A(d,ω)=d†[α-E†B⁻¹E]d and v=2π/q, Wc_00=-4π²A/(1+2πqA) and Wc_l0=-2π t_l·d/(1+2πqA)
!>     are integrated radially by Gauss-Legendre.
!> \param bs_env ...
!> \param directions unit directions d around Gamma
!> \param weights angular weights, zero for the mesh samples
!> \param radii extent of the near-Gamma cell along d, or |q| of a mesh sample
!> \param n_angular number of angular directions, followed by the mesh samples
!> \param head_response ε⁻¹_00(d,ω) in 3D, A(d,ω) in 2D
!> \param model Gamma correction per real-space cell
! **************************************************************************************************
   SUBROUTINE gamma_corr_cell_integrals(bs_env, directions, weights, radii, n_angular, head_response, model)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: directions(:, :), weights(:), radii(:)
      INTEGER, INTENT(IN)                                :: n_angular
      REAL(KIND=dp), INTENT(IN)                          :: head_response(:, :)
      TYPE(gamma_correction_type), INTENT(INOUT)         :: model

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_cell_integrals'
      INTEGER, PARAMETER                                 :: n_radial = 32

      COMPLEX(KIND=dp)                                   :: f0, f1, head, integral, phase, wing
      COMPLEX(KIND=dp), ALLOCATABLE                      :: head_frequency(:, :)
      INTEGER                                            :: handle, i_cell, i_direction, &
                                                            i_frequency, i_node, n_cells, &
                                                            n_frequency, nk
      INTEGER, ALLOCATABLE                               :: cells(:, :)
      LOGICAL                                            :: is_slab
      REAL(KIND=dp)                                      :: area, b(3, 2), head_scale, length, &
                                                            measure, projection, q, rcart(3), &
                                                            response
      REAL(KIND=dp), ALLOCATABLE                         :: nodes(:), rule(:)

      CALL timeset(routineN, handle)

      is_slab = ALL(bs_env%periodic == [1, 1, 0])
      CALL build_full_cell_grid(bs_env, cells)
      nk = PRODUCT(bs_env%nkp_grid_chi_eps_W_orig)
      n_cells = SIZE(cells, 2)
      n_frequency = bs_env%num_time_freq_points
      head_scale = 1.0_dp/ri_metric_head_zero(bs_env)**2
      ! Volume element of the angular rule: Ω/(2π)³ times the 4π of the solid angle in 3D,
      ! A/(2π)² in 2D.
      IF (is_slab) THEN
         CALL slab_reciprocal_vectors(bs_env, b, area)
         measure = area/twopi**2
         ALLOCATE (nodes(n_radial), rule(n_radial))
         CALL gauss_legendre(n_radial, nodes, rule)
      ELSE
         measure = bs_env%ri_rs%cell%deth*fourpi/twopi**3
      END IF
      ALLOCATE (head_frequency(n_cells, n_frequency), SOURCE=z_zero)
      ALLOCATE (model%wing_cell_frequency(3, n_frequency, n_cells), SOURCE=z_zero)
      DO i_cell = bs_env%para_env%mepos + 1, n_cells, bs_env%para_env%num_pe
         rcart = MATMUL(bs_env%hmat, REAL(cells(:, i_cell), dp))
         DO i_direction = 1, SIZE(radii)
            length = radii(i_direction)
            projection = DOT_PRODUCT(directions(:, i_direction), rcart)
            phase = EXP(-gaussi*length*projection)
            IF (.NOT. is_slab) THEN
               IF (i_direction <= n_angular) THEN
                  CALL radial_phase_moments(projection, length, f0, f1)
                  f0 = measure*weights(i_direction)*f0
                  f1 = measure*weights(i_direction)*f1
               ELSE
                  f0 = -phase/(REAL(nk, dp)*length**2)
                  f1 = -phase/(REAL(nk, dp)*length)
               END IF
            END IF
            DO i_frequency = 1, n_frequency
               response = head_response(i_direction, i_frequency)
               IF (is_slab) THEN
                  IF (i_direction <= n_angular) THEN
                     integral = z_zero
                     DO i_node = 1, n_radial
                        q = 0.5_dp*length*(nodes(i_node) + 1.0_dp)
                        integral = integral + 0.5_dp*length*rule(i_node)*q/(1.0_dp + twopi*q*response)* &
                                   EXP(-gaussi*q*projection)
                     END DO
                     integral = measure*weights(i_direction)*integral
                  ELSE
                     integral = -phase/(REAL(nk, dp)*(1.0_dp + twopi*length*response))
                  END IF
                  head = -twopi**2*response*integral
                  wing = -integral
               ELSE
                  head = fourpi*(response - 1.0_dp)*f0
                  wing = -response*f1
               END IF
               head_frequency(i_cell, i_frequency) = head_frequency(i_cell, i_frequency) + head_scale*head
               model%wing_cell_frequency(:, i_frequency, i_cell) = &
                  model%wing_cell_frequency(:, i_frequency, i_cell) + wing*directions(:, i_direction)
            END DO
         END DO
      END DO
      CALL bs_env%para_env%sum(head_frequency)
      CALL bs_env%para_env%sum(model%wing_cell_frequency)
      ALLOCATE (model%head_cell_tau(n_cells, n_frequency))
      model%head_cell_tau(:, :) = MATMUL(head_frequency, &
                                         TRANSPOSE(bs_env%time_frequency_grid%cosine_frequency_to_time))
      DEALLOCATE (cells, head_frequency)
      IF (is_slab) DEALLOCATE (nodes, rule)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_cell_integrals

! **************************************************************************************************
!> \brief Add the Gamma-correction head and wings of one time point to the local grid pairs.
!> The wing of cell R is w_l(R,τ) = Σ_ω cos(τ,ω) Σ_c wing_cartesian_lc(ω) c_c(R,ω); one GEMM per cell
!> forms it together with its conjugate-side partner.
!> \param bs_env ...
!> \param cells ...
!> \param model ...
!> \param i_t ...
!> \param Z_head Grid coefficients Z_l0(0) of the Gamma head function
!> \param first_pair ...
!> \param W_grid_R ...
! **************************************************************************************************
   SUBROUTINE add_gamma_corr_cell_correction(bs_env, cells, model, i_t, Z_head, first_pair, W_grid_R)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: cells(:, :)
      TYPE(gamma_correction_type), INTENT(IN)            :: model
      INTEGER, INTENT(IN)                                :: i_t
      COMPLEX(KIND=dp), INTENT(IN)                       :: Z_head(:)
      INTEGER, INTENT(IN)                                :: first_pair
      REAL(KIND=dp), INTENT(INOUT)                       :: W_grid_R(:, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'add_gamma_corr_cell_correction'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: coefficient(:, :, :), wing(:, :)
      INTEGER                                            :: handle, i_cell, j_w, n_frequency, n_grid

      CALL timeset(routineN, handle)

      n_grid = SIZE(model%wing_cartesian, 1)
      n_frequency = SIZE(model%wing_cartesian, 3)
      ALLOCATE (coefficient(3, n_frequency, 2), wing(n_grid, 2))
      DO i_cell = 1, SIZE(cells, 2)
         DO j_w = 1, n_frequency
            coefficient(:, j_w, 1) = bs_env%time_frequency_grid%cosine_frequency_to_time(i_t, j_w)* &
                                     model%wing_cell_frequency(:, j_w, i_cell)
         END DO
         ! conj(wing_cartesian) c = conj(wing_cartesian conj(c))
         coefficient(:, :, 2) = CONJG(coefficient(:, :, 1))
         CALL zgemm('N', 'N', n_grid, 2, 3*n_frequency, z_one, model%wing_cartesian, n_grid, &
                    coefficient, 3*n_frequency, z_zero, wing, n_grid)
         wing(:, 2) = CONJG(wing(:, 2))
         CALL add_gamma_corr_pair_cell(bs_env, cells, i_cell, model%head_cell_tau(i_cell, i_t), Z_head, &
                                       first_pair, W_grid_R, wing(:, 1), wing(:, 2))
      END DO
      DEALLOCATE (coefficient, wing)

      CALL timestop(handle)

   END SUBROUTINE add_gamma_corr_cell_correction

! **************************************************************************************************
!> \brief Transform cached response and solve (I+F^H Z^H B Z F) X = F^H Z^H E.
!> Z_lP and F are existing metric maps. All three right sides use the same
!> row-distributed frequency matrix. Head and wings use the same time transform.
!> \param bs_env ...
!> \param chi_tau Owned grid rows of B(tau).
!> \param wing_tau Owned grid rows of E(tau).
!> \param head_tau P(tau).
!> \param metric_rows Owned grid rows of T = Z F.
!> \param first_row First owned grid row.
!> \param j_w Frequency index.
!> \param E ...
!> \param P ...
!> \param X ...
! **************************************************************************************************
   SUBROUTINE gamma_ao_time_response(bs_env, chi_tau, wing_tau, head_tau, metric_rows, first_row, &
                                     j_w, E, P, X)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), INTENT(IN)                       :: chi_tau(:, :, :), wing_tau(:, :, :), &
                                                            head_tau(:, :, :), metric_rows(:, :)
      INTEGER, INTENT(IN)                                :: first_row, j_w
      COMPLEX(KIND=dp), INTENT(OUT)                      :: E(:, :), P(:, :), X(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_ao_time_response'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: chi_frequency(:, :), wing_frequency(:, :)
      INTEGER                                            :: handle, nr, t

      CALL timeset(routineN, handle)

      ! Apply the existing cosine coefficients to the local cache, without rebuilding Green matrices.
      nr = SIZE(chi_tau, 1)
      ALLOCATE (chi_frequency(nr, SIZE(chi_tau, 2)), wing_frequency(nr, 3), SOURCE=z_zero)
      P = z_zero
      DO t = 1, bs_env%num_time_freq_points
         chi_frequency(:, :) = chi_frequency + &
                               bs_env%time_frequency_grid%cosine_time_to_frequency(j_w, t)*chi_tau(:, :, t)
         wing_frequency(:, :) = wing_frequency + &
                                bs_env%time_frequency_grid%cosine_time_to_frequency(j_w, t)*wing_tau(:, :, t)
         P = P + bs_env%time_frequency_grid%cosine_time_to_frequency(j_w, t)*head_tau(:, :, t)
      END DO
      ! Project the wing as E_Qc = sum_l conjg(T_lQ) E_lc over disjoint rows.
      E = z_zero
      IF (nr > 0) CALL zgemm('C', 'N', SIZE(metric_rows, 2), 3, nr, z_one, metric_rows, nr, &
                             wing_frequency, nr, z_zero, E, SIZE(E, 1))
      CALL bs_env%para_env%sum(E)
      CALL solve_gamma_corr_body(bs_env, metric_rows, chi_frequency, first_row, E, X)

      CALL timestop(handle)

   END SUBROUTINE gamma_ao_time_response

! **************************************************************************************************
!> \brief Warn if sampled transform errors over the physical transition interval exceed the tolerance.
!>        The interval spans the smallest gap and the largest transition energy of all k points.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE validate_time_transform(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      INTEGER, PARAMETER                                 :: n_sampled_gaps = 32
      REAL(KIND=dp), PARAMETER :: eps_gamma_time_transform = 1.0E-5_dp

      CHARACTER(LEN=256)                                 :: error_message
      INTEGER                                            :: i, no, nv, spin, t, w
      REAL(KIND=dp)                                      :: approximate, delta, exact, hi, lo, &
                                                            omega, worst
      TYPE(gwkp_type), POINTER                           :: gwkp_env

      gwkp_env => bs_env%gwkp_env

      lo = HUGE(1.0_dp)
      hi = 0.0_dp
      DO spin = 1, bs_env%n_spin
         no = bs_env%n_occ(spin)
         nv = bs_env%n_mo_retained - no
         CPASSERT(nv >= 0)
         IF (no == 0 .OR. nv == 0) CYCLE
         lo = MIN(lo, MINVAL(gwkp_env%eigenvalues_G(no + 1, 1:bs_env%nkp_scf_desymm, spin) - &
                             gwkp_env%eigenvalues_G(no, 1:bs_env%nkp_scf_desymm, spin)))
         hi = MAX(hi, MAXVAL(gwkp_env%eigenvalues_G(no + nv, 1:bs_env%nkp_scf_desymm, spin) - &
                             gwkp_env%eigenvalues_G(1, 1:bs_env%nkp_scf_desymm, spin)))
      END DO
      IF (hi == 0.0_dp) RETURN
      IF (lo <= 0.0_dp) THEN
         CPABORT('Gamma response requires a positive transition gap')
      END IF
      worst = 0.0_dp
      DO i = 1, n_sampled_gaps
         delta = lo*(hi/lo)**(REAL(i - 1, dp)/REAL(n_sampled_gaps - 1, dp))
         DO w = 1, bs_env%num_time_freq_points
            omega = bs_env%time_frequency_grid%frequency(w)
            exact = 2.0_dp*delta/(delta*delta + omega*omega)
            approximate = 0.0_dp
            DO t = 1, bs_env%num_time_freq_points
               approximate = approximate + bs_env%time_frequency_grid%cosine_time_to_frequency(w, t)* &
                             EXP(-delta*bs_env%time_frequency_grid%imaginary_time(t))
            END DO
            worst = MAX(worst, ABS(approximate - exact)/exact)
         END DO
      END DO
      IF (worst > eps_gamma_time_transform) THEN
         WRITE (error_message, '(A,ES12.4,A,ES12.4,A)') &
            'Gamma time transform relative error ', worst, &
            ' exceeds tolerance ', eps_gamma_time_transform, &
            '. Increase NUM_TIME_FREQ_POINTS and check convergence.'
         CPWARN(TRIM(error_message))
      END IF
   END SUBROUTINE validate_time_transform

! **************************************************************************************************
!> \brief Cache Gamma-correction response in imaginary time, distributed by RI-RS grid rows.
!>
!> Let psi_li=sum_mu phi_mu(r_l) C_mu,i(k), and let m be the midpoint
!> between the highest occupied and lowest virtual eigenvalue at each k/spin.
!> Gocc_ll'(tau)=sum_i psi_li conjg(psi_l'i) exp[(epsilon_i-m)tau].
!> Gvir_ll'(tau)=sum_a psi_la conjg(psi_l'a) exp[-(epsilon_a-m)tau].
!> These are the projected AO Green matrices, evaluated associatively.
!> With f_k=spin_degeneracy*w_k and d_ia,c=-i*dipole_ia,c/sqrt(volume):
!> B_ll'(tau)=sum_k f_k conjg(Gocc_ll') Gvir_ll'.
!> E_lc(tau)=sum_kia f_k conjg(psi_li) psi_la d_ia,c exp[-Delta_ia*tau].
!> P_cd(tau)=sum_kia f_k conjg(d_ia,c) d_ia,d exp[-Delta_ia*tau].
!> Each rank owns disjoint rows of B and E; only the small P is replicated.
!> A k-owner broadcasts its orbitals once, then all ranks process owned rows.
!> No RI-pair by transition tensor or replicated grid-pair cache is formed.
!> The cache is created once and reused for all frequencies. The existing
!> minimax cosine transform approximates 2*Delta/(Delta**2+omega**2);
!> its relative kernel error is checked independently of the iterative solve.
!> \param bs_env ...
!> \param psi ...
!> \param dipole ...
!> \param volume ...
!> \param first_row First owned grid row.
!> \param chi_tau Owned grid rows of B(tau).
!> \param wing_tau Owned grid rows of E(tau).
!> \param head_tau P(tau).
! **************************************************************************************************
   SUBROUTINE accumulate_time_state(bs_env, psi, dipole, volume, first_row, chi_tau, wing_tau, head_tau)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), INTENT(IN)                       :: psi(:, :, :, :), dipole(:, :, :, :, :)
      REAL(KIND=dp), INTENT(IN)                          :: volume
      INTEGER, INTENT(OUT)                               :: first_row
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: chi_tau(:, :, :), wing_tau(:, :, :), &
                                                            head_tau(:, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'accumulate_time_state'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: gamma_corr(:, :, :), orbitals(:, :)
      INTEGER                                            :: c, handle, k, local_k, n_grid, n_rows, &
                                                            n_time, no, nv, owner, spin

      CALL timeset(routineN, handle)

      CPASSERT(volume > 0.0_dp)
      n_grid = bs_env%ri_rs%n_grid_points
      n_time = bs_env%num_time_freq_points
      CALL block_range(n_grid, bs_env%para_env%mepos, bs_env%para_env%num_pe, first_row, n_rows)
      ALLOCATE (chi_tau(n_rows, n_grid, n_time), wing_tau(n_rows, 3, n_time), &
                head_tau(3, 3, n_time), SOURCE=z_zero)
      CALL validate_time_transform(bs_env)
      DO spin = 1, bs_env%n_spin
         no = bs_env%n_occ(spin)
         nv = bs_env%n_mo_retained - no
         CPASSERT(nv >= 0)
         IF (no == 0 .OR. nv == 0) CYCLE
         ALLOCATE (orbitals(bs_env%ri_rs%n_grid_points, no + nv), gamma_corr(no, nv, 3))
         DO k = 1, bs_env%nkp_scf_desymm
            owner = MODULO(k - 1, bs_env%para_env%num_pe)
            IF (bs_env%para_env%mepos == owner) THEN
               local_k = (k - 1)/bs_env%para_env%num_pe + 1
               CPASSERT(local_k <= SIZE(psi, 3))
               orbitals(:, :) = psi(:, 1:no + nv, local_k, spin)
               DO c = 1, 3
                  gamma_corr(:, :, c) = -gaussi*dipole(spin, local_k, c, 1:no, 1:nv)/SQRT(volume)
               END DO
            END IF
            CALL bs_env%para_env%bcast(orbitals, owner)
            CALL bs_env%para_env%bcast(gamma_corr, owner)
            CALL accumulate_k_time_rows(bs_env, orbitals, gamma_corr, k, spin, owner, &
                                        first_row, chi_tau, wing_tau, head_tau)
         END DO
         DEALLOCATE (orbitals, gamma_corr)
      END DO
      CALL bs_env%para_env%sum(head_tau)

      CALL timestop(handle)

   END SUBROUTINE accumulate_time_state

! **************************************************************************************************
!> \brief Weight only owned orbital rows in Green and Gamma-correction products; accumulate the tiny head once.
!> Gocc_ll'=sum_i [psi_li exp((epsilon_i-m)tau)] conjg(psi_l'i), and analogously for Gvir.
!> The wing uses both fully weighted local factors: conjg(psi_li) psi_la d_ia exp(-Delta_ia*tau).
!> Global orbitals remain unweighted, avoiding a replicated time-weighted orbital array.
!> \param bs_env ...
!> \param orbitals ...
!> \param gamma_corr Transition moments entering the Gamma correction.
!> \param k ...
!> \param spin ...
!> \param owner ...
!> \param first_row ...
!> \param chi_tau ...
!> \param wing_tau ...
!> \param head_tau ...
! **************************************************************************************************
   SUBROUTINE accumulate_k_time_rows(bs_env, orbitals, gamma_corr, k, spin, owner, &
                                     first_row, chi_tau, wing_tau, head_tau)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), INTENT(IN)                       :: orbitals(:, :), gamma_corr(:, :, :)
      INTEGER, INTENT(IN)                                :: k, spin, owner, first_row
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: chi_tau(:, :, :), wing_tau(:, :, :), &
                                                            head_tau(:, :, :)

      COMPLEX(KIND=dp), ALLOCATABLE                      :: block(:, :), go(:, :), gv(:, :), &
                                                            scaled(:, :)
      INTEGER                                            :: a, c, d, first, i, last, ng, no, nr, nv, &
                                                            t
      REAL(KIND=dp)                                      :: head_weight, mid, tau, weight
      REAL(KIND=dp), ALLOCATABLE                         :: occ_factor(:), vir_factor(:)
      TYPE(gwkp_type), POINTER                           :: gwkp_env

      gwkp_env => bs_env%gwkp_env

      no = bs_env%n_occ(spin)
      nv = bs_env%n_mo_retained - no
      CPASSERT(nv >= 0)
      ng = bs_env%ri_rs%n_grid_points
      nr = SIZE(chi_tau, 1)
      first = first_row
      last = first + nr - 1
      mid = 0.5_dp*(gwkp_env%eigenvalues_G(no, k, spin) + gwkp_env%eigenvalues_G(no + 1, k, spin))
      weight = bs_env%spin_degeneracy*bs_env%kpoints_scf_desymm%wkp(k)
      ALLOCATE (occ_factor(no), vir_factor(nv), scaled(nr, no + nv))
      ALLOCATE (go(nr, ng), gv(nr, ng), block(nr, no))
      DO t = 1, bs_env%num_time_freq_points
         tau = bs_env%time_frequency_grid%imaginary_time(t)
         DO i = 1, no
            occ_factor(i) = EXP(0.5_dp*(gwkp_env%eigenvalues_G(i, k, spin) - mid)*tau)
            scaled(:, i) = orbitals(first:last, i)*occ_factor(i)**2
         END DO
         DO a = 1, nv
            vir_factor(a) = EXP(-0.5_dp*(gwkp_env%eigenvalues_G(no + a, k, spin) - mid)*tau)
            scaled(:, no + a) = orbitals(first:last, no + a)*vir_factor(a)**2
         END DO
         IF (nr > 0) THEN
            CALL zgemm('N', 'C', nr, ng, no, z_one, scaled, nr, &
                       orbitals, ng, z_zero, go, nr)
            CALL zgemm('N', 'C', nr, ng, nv, z_one, scaled(1, no + 1), nr, &
                       orbitals(:, no + 1:), ng, z_zero, gv, nr)
            chi_tau(:, :, t) = chi_tau(:, :, t) + weight*CONJG(go)*gv
            DO c = 1, 3
               CALL zgemm('N', 'T', nr, no, nv, z_one, scaled(1, no + 1), nr, &
                          gamma_corr(:, :, c), no, z_zero, block, nr)
               wing_tau(:, c, t) = wing_tau(:, c, t) + &
                                   weight*SUM(CONJG(scaled(:, 1:no))*block, DIM=2)
            END DO
         END IF
         IF (bs_env%para_env%mepos == owner) THEN
            DO a = 1, nv
               DO i = 1, no
                  head_weight = weight*(occ_factor(i)*vir_factor(a))**2
                  DO d = 1, 3
                     DO c = 1, 3
                        head_tau(c, d, t) = head_tau(c, d, t) + &
                                            head_weight*CONJG(gamma_corr(i, a, c))*gamma_corr(i, a, d)
                     END DO
                  END DO
               END DO
            END DO
         END IF
      END DO
      DEALLOCATE (occ_factor, vir_factor, scaled, go, gv, block)
   END SUBROUTINE accumulate_k_time_rows

! **************************************************************************************************
!> \brief Form T_lQ = sum_P Z_lP F_PQ once on the owned grid rows.
!> The same T is reused for all frequencies; no replicated grid-by-RI map is formed.
!> \param Z ...
!> \param F ...
!> \param first_row First owned grid row.
!> \param nr Number of owned grid rows.
!> \param metric_rows ...
! **************************************************************************************************
   SUBROUTINE metric_grid_rows(Z, F, first_row, nr, metric_rows)
      COMPLEX(KIND=dp), INTENT(IN)                       :: Z(:, :), F(:, :)
      INTEGER, INTENT(IN)                                :: first_row, nr
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: metric_rows(:, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'metric_grid_rows'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: Z_rows(:, :)
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ALLOCATE (metric_rows(nr, SIZE(F, 2)))
      IF (nr > 0) THEN
         ALLOCATE (Z_rows, SOURCE=Z(first_row:first_row + nr - 1, :))
         CALL zgemm('N', 'N', nr, SIZE(F, 2), SIZE(F, 1), z_one, &
                    Z_rows, nr, F, SIZE(F, 1), z_zero, metric_rows, nr)
         DEALLOCATE (Z_rows)
      END IF

      CALL timestop(handle)

   END SUBROUTINE metric_grid_rows

! **************************************************************************************************
!> \brief Apply T^H chi T with distributed rows and reductions of three-column vectors.
!> The grid vector is assembled from disjoint row contributions, never from replicated T.
!> \param bs_env ...
!> \param T_rows ...
!> \param chi ...
!> \param first_row ...
!> \param x ...
!> \param y ...
!> \param aux_y ...
!> \param grid_x ...
!> \param grid_y ...
! **************************************************************************************************
   SUBROUTINE apply_gamma_corr_body(bs_env, T_rows, chi, first_row, x, y, aux_y, grid_x, grid_y)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), INTENT(IN)                       :: T_rows(:, :), chi(:, :)
      INTEGER, INTENT(IN)                                :: first_row
      COMPLEX(KIND=dp), INTENT(IN)                       :: x(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: y(:, :)
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: aux_y(:, :), grid_x(:, :), grid_y(:, :)

      INTEGER                                            :: i_last, nc, nr

      nr = SIZE(chi, 1)
      nc = SIZE(x, 2)
      i_last = first_row + nr - 1
      grid_x = z_zero
      IF (nr > 0) THEN
         CALL zgemm('N', 'N', nr, nc, SIZE(T_rows, 2), z_one, T_rows, nr, &
                    x, SIZE(x, 1), z_zero, grid_y, nr)
         grid_x(first_row:i_last, :) = grid_y
      END IF
      CALL bs_env%para_env%sum(grid_x)
      aux_y = z_zero
      IF (nr > 0) THEN
         CALL zgemm('N', 'N', nr, nc, SIZE(chi, 2), z_one, chi, nr, &
                    grid_x, SIZE(grid_x, 1), z_zero, grid_y, nr)
         CALL zgemm('C', 'N', SIZE(T_rows, 2), nc, nr, z_one, T_rows, nr, &
                    grid_y, nr, z_zero, aux_y, SIZE(aux_y, 1))
      END IF
      CALL bs_env%para_env%sum(aux_y)
      y = aux_y
   END SUBROUTINE apply_gamma_corr_body

! **************************************************************************************************
!> \brief Independent conjugate gradients for three directions, sharing matrix actions and reductions.
!> Each column c retains its own alpha_c, beta_c and convergence test; columns are not coupled.
!> \param bs_env ...
!> \param T_rows ...
!> \param chi ...
!> \param first_row ...
!> \param b ...
!> \param x ...
! **************************************************************************************************
   SUBROUTINE solve_gamma_corr_body(bs_env, T_rows, chi, first_row, b, x)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      COMPLEX(KIND=dp), INTENT(IN)                       :: T_rows(:, :), chi(:, :)
      INTEGER, INTENT(IN)                                :: first_row
      COMPLEX(KIND=dp), INTENT(IN)                       :: b(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: x(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_gamma_corr_body'
      REAL(KIND=dp), PARAMETER                           :: tolerance = 1.0E-12_dp

      COMPLEX(KIND=dp)                                   :: denominator
      COMPLEX(KIND=dp), ALLOCATABLE                      :: A_direction(:, :), aux_y(:, :), &
                                                            direction(:, :), grid_x(:, :), &
                                                            grid_y(:, :), r(:, :)
      INTEGER                                            :: handle, i_c, i_iter, nc
      REAL(KIND=dp)                                      :: alpha, rr(SIZE(b, 2)), &
                                                            rr_converged(SIZE(b, 2)), rr_new

      CALL timeset(routineN, handle)

      nc = SIZE(b, 2)
      x = z_zero
      DO i_c = 1, nc
         rr(i_c) = REAL(DOT_PRODUCT(b(:, i_c), b(:, i_c)), dp)
      END DO
      rr_converged = tolerance**2*rr
      ALLOCATE (grid_x(SIZE(chi, 2), nc), grid_y(SIZE(chi, 1), nc), aux_y(SIZE(T_rows, 2), nc))
      ALLOCATE (r, SOURCE=b)
      ALLOCATE (direction, SOURCE=b)
      ALLOCATE (A_direction(SIZE(b, 1), nc))
      DO i_iter = 1, 4*SIZE(b, 1)
         IF (ALL(rr <= rr_converged)) EXIT
         CALL apply_gamma_corr_body(bs_env, T_rows, chi, first_row, direction, A_direction, aux_y, grid_x, grid_y)
         A_direction(:, :) = A_direction + direction
         DO i_c = 1, nc
            IF (rr(i_c) <= rr_converged(i_c)) CYCLE
            denominator = DOT_PRODUCT(direction(:, i_c), A_direction(:, i_c))
            IF (REAL(denominator, dp) <= 0.0_dp) THEN
               CPABORT('Gamma body action is not positive')
            END IF
            alpha = rr(i_c)/REAL(denominator, dp)
            x(:, i_c) = x(:, i_c) + alpha*direction(:, i_c)
            r(:, i_c) = r(:, i_c) - alpha*A_direction(:, i_c)
            rr_new = REAL(DOT_PRODUCT(r(:, i_c), r(:, i_c)), dp)
            direction(:, i_c) = r(:, i_c) + (rr_new/rr(i_c))*direction(:, i_c)
            rr(i_c) = rr_new
         END DO
      END DO
      IF (ANY(rr > rr_converged)) THEN
         CALL cp_warn(__LOCATION__, "Gamma-point body CG not converged after "// &
                      TRIM(cp_to_string(4*SIZE(b, 1)))//" iterations")
      END IF

      CALL timestop(handle)

   END SUBROUTINE solve_gamma_corr_body

! **************************************************************************************************
!> \brief Contract M_pq=Σ_ABPQ U^A_Pp [Σ_R S_AB,PQ(R)] U^B_Qq using atom-local U (1).
!> P,Q enumerate reference functions on atoms A,B; p,q enumerate their optimized
!> functions. Each distributed S element contributes once. No atom pair, image,
!> or coefficient is screened. Sum the resulting optimized columns over MPI (2).
!> Only atom-local transformation rows and one atom-pair intermediate are stored;
!> the final optimized M remains replicated. This is exactly the real q=0 metric.
!> \param images ...
!> \param auto_ri ...
!> \param metric ...
! **************************************************************************************************
   SUBROUTINE gamma_metric_contract_atoms(images, auto_ri, metric)
      TYPE(cp_fm_type), INTENT(IN)                       :: images(:)
      TYPE(auto_ri_type), INTENT(IN)                     :: auto_ri
      COMPLEX(KIND=dp), INTENT(OUT)                      :: metric(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_metric_contract_atoms'

      INTEGER                                            :: handle
      TYPE(gamma_atom_projection_type), ALLOCATABLE      :: left(:), right(:)

      CALL timeset(routineN, handle)

      CALL prepare_gamma_atom_projections(images, auto_ri, left, right)
      CALL accumulate_gamma_atom_metric(images, left, right, metric)
      DEALLOCATE (left, right)

      CALL timestop(handle)

   END SUBROUTINE gamma_metric_contract_atoms

! **************************************************************************************************
!> \brief Extract the distributed reference rows of each atom-local U in equation (1).
!> \param images ...
!> \param auto_ri ...
!> \param left ...
!> \param right ...
! **************************************************************************************************
   SUBROUTINE prepare_gamma_atom_projections(images, auto_ri, left, right)
      TYPE(cp_fm_type), INTENT(IN)                       :: images(:)
      TYPE(auto_ri_type), INTENT(IN)                     :: auto_ri
      TYPE(gamma_atom_projection_type), ALLOCATABLE, &
         INTENT(OUT)                                     :: left(:), right(:)

      INTEGER                                            :: nc, ngc, ngr, nr
      INTEGER, POINTER                                   :: cols(:), rows(:)
      LOGICAL                                            :: valid_layout

      CPASSERT(SIZE(images) > 0)
      CPASSERT(auto_ri%ready)
      CPASSERT(ALLOCATED(auto_ri%sizes_ref_RI))
      CPASSERT(ALLOCATED(auto_ri%sizes_opt_RI))
      CPASSERT(ALLOCATED(auto_ri%U_Pp_atom))
      CPASSERT(ALLOCATED(auto_ri%U_Pp_atom_offset))
      CPASSERT(SIZE(auto_ri%sizes_ref_RI) > 0)
      CPASSERT(SIZE(auto_ri%sizes_ref_RI) == SIZE(auto_ri%sizes_opt_RI))
      CPASSERT(SIZE(auto_ri%U_Pp_atom_offset) == SIZE(auto_ri%sizes_ref_RI))
      valid_layout = ALL(auto_ri%sizes_ref_RI > 0) .AND. ALL(auto_ri%sizes_opt_RI > 0)
      CPASSERT(valid_layout)
      CALL cp_fm_get_info(images(1), nrow_local=nr, ncol_local=nc, &
                          nrow_global=ngr, ncol_global=ngc, row_indices=rows, col_indices=cols)
      CPASSERT(ngr == SUM(auto_ri%sizes_ref_RI) .AND. ngc == ngr)
      CALL extract_gamma_atom_rows(auto_ri, rows(1:nr), left)
      CALL extract_gamma_atom_rows(auto_ri, cols(1:nc), right)
   END SUBROUTINE prepare_gamma_atom_projections

! **************************************************************************************************
!> \brief Map local reference indices onto their atom blocks in equation (1).
!> \param auto_ri ...
!> \param indices ...
!> \param blocks ...
! **************************************************************************************************
   SUBROUTINE extract_gamma_atom_rows(auto_ri, indices, blocks)
      TYPE(auto_ri_type), INTENT(IN)                     :: auto_ri
      INTEGER, INTENT(IN)                                :: indices(:)
      TYPE(gamma_atom_projection_type), ALLOCATABLE, &
         INTENT(OUT)                                     :: blocks(:)

      INTEGER                                            :: first_opt, first_ref, iatom, icol, &
                                                            ilocal, ipos, irow, nopt, nref, offset

      ALLOCATE (blocks(SIZE(auto_ri%sizes_ref_RI)))
      first_ref = 1
      first_opt = 1
      DO iatom = 1, SIZE(blocks)
         nref = auto_ri%sizes_ref_RI(iatom)
         nopt = auto_ri%sizes_opt_RI(iatom)
         offset = auto_ri%U_Pp_atom_offset(iatom)
         CPASSERT(offset > 0 .AND. nref <= HUGE(0)/nopt)
         CPASSERT(offset - 1 <= SIZE(auto_ri%U_Pp_atom) - nref*nopt)
         ALLOCATE (blocks(iatom)%positions(COUNT(indices >= first_ref .AND. indices < first_ref + nref)))
         ALLOCATE (blocks(iatom)%values(SIZE(blocks(iatom)%positions), nopt))
         blocks(iatom)%first = first_opt
         ipos = 0
         DO ilocal = 1, SIZE(indices)
            irow = indices(ilocal) - first_ref + 1
            IF (irow < 1 .OR. irow > nref) CYCLE
            ipos = ipos + 1
            blocks(iatom)%positions(ipos) = ilocal
            DO icol = 1, nopt
               blocks(iatom)%values(ipos, icol) = auto_ri%U_Pp_atom(offset + irow - 1 + (icol - 1)*nref)
            END DO
         END DO
         first_ref = first_ref + nref
         first_opt = first_opt + nopt
      END DO
   END SUBROUTINE extract_gamma_atom_rows

! **************************************************************************************************
!> \brief Evaluate atom-pair products (1) locally and reduce optimized columns (2).
!> \param images ...
!> \param left ...
!> \param right ...
!> \param metric ...
! **************************************************************************************************
   SUBROUTINE accumulate_gamma_atom_metric(images, left, right, metric)
      TYPE(cp_fm_type), INTENT(IN)                       :: images(:)
      TYPE(gamma_atom_projection_type), INTENT(IN)       :: left(:), right(:)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: metric(:, :)

      INTEGER                                            :: ia, ib, icell, n_opt_total, nc, noleft, &
                                                            nopt, noright, nr, nref, pfirst, qfirst
      REAL(KIND=dp), ALLOCATABLE                         :: block(:, :), block_right(:, :), &
                                                            value(:, :)

      nopt = MAXVAL([(SIZE(left(ia)%values, 2), ia=1, SIZE(left))])
      nref = MAX(1, MAXVAL([(SIZE(left(ia)%positions), ia=1, SIZE(left))]), &
                 MAXVAL([(SIZE(right(ib)%positions), ib=1, SIZE(right))]))
      n_opt_total = SUM([(SIZE(left(ia)%values, 2), ia=1, SIZE(left))])
      CPASSERT(ALL(SHAPE(metric) == n_opt_total))
      ALLOCATE (block(nref, nref), block_right(nref, nopt), value(nopt, nopt))
      metric = z_zero
      DO ib = 1, SIZE(right)
         nc = SIZE(right(ib)%positions)
         noright = SIZE(right(ib)%values, 2)
         qfirst = right(ib)%first
         DO ia = 1, SIZE(left)
            nr = SIZE(left(ia)%positions)
            IF (nr == 0 .OR. nc == 0) CYCLE
            noleft = SIZE(left(ia)%values, 2)
            pfirst = left(ia)%first
            block(1:nr, 1:nc) = 0.0_dp
            DO icell = 1, SIZE(images)
               block(1:nr, 1:nc) = block(1:nr, 1:nc) + &
                                   images(icell)%local_data(left(ia)%positions, right(ib)%positions)
            END DO
            CALL dgemm('N', 'N', nr, noright, nc, 1.0_dp, block, nref, &
                       right(ib)%values, nc, 0.0_dp, block_right, nref)
            CALL dgemm('T', 'N', noleft, noright, nr, 1.0_dp, left(ia)%values, nr, &
                       block_right, nref, 0.0_dp, value, nopt)
            metric(pfirst:pfirst + noleft - 1, qfirst:qfirst + noright - 1) = &
               CMPLX(value(1:noleft, 1:noright), 0.0_dp, dp)
         END DO
         CPASSERT(noright <= HUGE(0)/SIZE(metric, 1))
         CALL images(1)%matrix_struct%para_env%sum(metric(:, qfirst:qfirst + noright - 1))
      END DO
      DEALLOCATE (block, block_right, value)
   END SUBROUTINE accumulate_gamma_atom_metric

! **************************************************************************************************
!> \brief Construct F=M^-1 V^1/2 with distributed Hermitian eigensystems.
!> For each input A, ifirst set A=(A+Aᴴ)/2 (1), then diagonalize A=U diag(λ) Uᴴ.
!> Define A^a=U diag(f(λ)) Uᴴ (2), where f(λ)=λ^a only for λ>max(0,η λ_max),
!> and zero otherwise. M uses a=-1 and the existing relative cutoff η; V uses
!> a=1/2 and η=0. Finally F=M^-1 (V^1/2)ᴴ (3), matching the dense Γ path.
!> Inputs and the final F remain replicated; eigenvectors and matrix products
!> are block-cyclic distributed. No spectral threshold or physical input changes.
!> \param para_env ...
!> \param relative_cutoff ...
!> \param metric ...
!> \param coulomb ...
!> \param factor ...
! **************************************************************************************************
   SUBROUTINE gamma_metric_factor_distributed(para_env, relative_cutoff, metric, coulomb, factor)
      TYPE(mp_para_env_type), INTENT(INOUT), TARGET      :: para_env
      REAL(KIND=dp), INTENT(IN)                          :: relative_cutoff
      COMPLEX(KIND=dp), INTENT(IN)                       :: metric(:, :), coulomb(:, :)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: factor(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_metric_factor_distributed'

      INTEGER                                            :: handle, n
      TYPE(cp_blacs_env_type), POINTER                   :: blacs
      TYPE(cp_cfm_type)                                  :: coulomb_d, factor_d, metric_d, work_d
      TYPE(cp_fm_struct_type), POINTER                   :: structure

      CALL timeset(routineN, handle)

      n = SIZE(metric, 1)
      CALL cp_blacs_env_create(blacs, para_env)
      CALL cp_fm_struct_create(structure, context=blacs, nrow_global=n, ncol_global=n)
      CALL cp_cfm_create(metric_d, structure)
      CALL cp_cfm_create(coulomb_d, structure)
      CALL cp_cfm_create(work_d, structure)
      CALL cp_cfm_create(factor_d, structure)
      CALL cp_cfm_set_submatrix(metric_d, metric)
      CALL cp_cfm_set_submatrix(coulomb_d, coulomb)
      CALL retained_hermitian_power(metric_d, work_d, factor_d, n, -1.0_dp, relative_cutoff)
      CALL retained_hermitian_power(coulomb_d, work_d, factor_d, n, 0.5_dp, 0.0_dp)
      CALL parallel_gemm('N', 'C', n, n, n, z_one, metric_d, coulomb_d, z_zero, factor_d)
      CALL collect_factor_columns(factor_d, factor)
      CALL cp_cfm_release(metric_d)
      CALL cp_cfm_release(coulomb_d)
      CALL cp_cfm_release(work_d)
      CALL cp_cfm_release(factor_d)
      CALL cp_fm_struct_release(structure)
      CALL cp_blacs_env_release(blacs)

      CALL timestop(handle)

   END SUBROUTINE gamma_metric_factor_distributed

! **************************************************************************************************
!> \brief Apply Hermitization (1) and the strictly retained spectral power (2).
!> \param matrix ...
!> \param eigenvectors ...
!> \param work ...
!> \param n ...
!> \param exponent ...
!> \param relative_cutoff ...
! **************************************************************************************************
   SUBROUTINE retained_hermitian_power(matrix, eigenvectors, work, n, exponent, relative_cutoff)
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: matrix, eigenvectors, work
      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), INTENT(IN)                          :: exponent, relative_cutoff

      COMPLEX(KIND=dp), ALLOCATABLE                      :: scales(:)
      INTEGER                                            :: i_mode
      REAL(KIND=dp)                                      :: cutoff
      REAL(KIND=dp), ALLOCATABLE                         :: eigenvalues(:)

      CALL cp_cfm_transpose(matrix, 'C', work)
      CALL cp_cfm_scale_and_add(CMPLX(0.5_dp, 0.0_dp, dp), matrix, CMPLX(0.5_dp, 0.0_dp, dp), work)
      ALLOCATE (eigenvalues(n), scales(n))
      CALL cp_cfm_heevd(matrix, eigenvectors, eigenvalues)
      cutoff = MAX(0.0_dp, relative_cutoff*MAXVAL(eigenvalues))
      scales = z_zero
      DO i_mode = 1, n
         IF (eigenvalues(i_mode) > cutoff) scales(i_mode) = CMPLX(eigenvalues(i_mode)**(0.5_dp*exponent), 0.0_dp, dp)
      END DO
      CALL cp_cfm_column_scale(eigenvectors, scales)
      CALL parallel_gemm('N', 'C', n, n, n, z_one, eigenvectors, eigenvectors, z_zero, matrix)
      DEALLOCATE (eigenvalues, scales)
   END SUBROUTINE retained_hermitian_power

! **************************************************************************************************
!> \brief Collect factor (3) by BLACS column blocks with bounded MPI element counts.
!> \param matrix ...
!> \param factor ...
! **************************************************************************************************
   SUBROUTINE collect_factor_columns(matrix, factor)
      TYPE(cp_cfm_type), INTENT(IN)                      :: matrix
      COMPLEX(KIND=dp), INTENT(OUT)                      :: factor(:, :)

      INTEGER                                            :: i_column, i_last, n_columns, n_rows

      n_rows = SIZE(factor, 1)
      CPASSERT(n_rows > 0)
      n_columns = MIN(matrix%matrix_struct%ncol_block, HUGE(0)/n_rows)
      CPASSERT(n_columns > 0)
      i_column = 1
      DO
         i_last = i_column + (MIN(n_columns, SIZE(factor, 2) - i_column + 1) - 1)
         CALL cp_cfm_get_submatrix(matrix, factor(:, i_column:i_last), start_col=i_column)
         IF (i_last == SIZE(factor, 2)) EXIT
         i_column = i_last + 1
      END DO
   END SUBROUTINE collect_factor_columns

! **************************************************************************************************
!> \brief Form Gamma-correction d_ia^α(k) for cyclically owned SCF k points.
!> Reuses the existing AO integral builders and Bloch transformation.
!> \brief Interband Gamma-correction moments in the same SCF MO gauge as the GW body.
!> The SCF Hamiltonian contains the nonlocal pseudopotential; its k derivative
!> therefore includes the nonlocal contribution without a separate commutator.
!> \param bs_env ...
!> \param gwkp_env Gamma-correction AO operators, released after their contraction.
!> \param kpoints ...
!> \param mo_coeff ...
!> \param eigenvalues ...
!> \param nocc_spin ...
!> \param dipole ...
! **************************************************************************************************
   SUBROUTINE gw_gamma_transition_moments(bs_env, gwkp_env, kpoints, mo_coeff, eigenvalues, nocc_spin, &
                                          dipole)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gwkp_type), INTENT(INOUT)                     :: gwkp_env
      TYPE(kpoint_type), INTENT(IN), POINTER             :: kpoints
      TYPE(cp_cfm_type), DIMENSION(:, :), INTENT(IN)     :: mo_coeff
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: eigenvalues
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nocc_spin
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :), INTENT(OUT)           :: dipole

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_gamma_transition_moments'
      REAL(KIND=dp), PARAMETER                           :: eps_degenerate = 1.0E-10_dp

      COMPLEX(KIND=dp), ALLOCATABLE :: bloch(:, :, :, :), coeff_owned(:, :, :, :), &
         operators_k(:, :, :), projected_d(:, :), projected_h(:, :), projected_s(:, :), &
         recv_buffer(:)
      INTEGER :: first_col, first_k, first_rank_col, handle, i_dir, i_k, i_k_local, i_mo_a, &
         i_mo_i, i_rank, i_spin, ib, n_ao, n_batch, n_col, n_img, n_k, n_k_local, n_mat, n_mo, &
         n_mo_matrix, n_rank_col, n_row_matrix, n_spin, n_virt, nocc
      INTEGER, ALLOCATABLE                               :: nmo_spin(:), recv_counts(:), &
                                                            recv_displs(:), send_counts(:), &
                                                            send_displs(:)
      INTEGER, DIMENSION(3)                              :: periodic
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: valid_layout
      REAL(KIND=dp)                                      :: energy_diff, hmat(3, 3)
      REAL(KIND=dp), ALLOCATABLE :: cell_cartesian(:, :), d_cols(:, :, :, :), h_cols(:, :, :, :), &
         phase_angle(:, :), phase_cos(:, :), phase_sin(:, :), s_cols(:, :, :, :), weight_im(:, :), &
         weight_re(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      cell => bs_env%ri_rs%cell
      para_env => bs_env%para_env
      index_to_cell => gwkp_env%gamma_corr_index_to_cell
      cell_to_index => gwkp_env%gamma_corr_cell_to_index
      CPASSERT(ASSOCIATED(index_to_cell) .AND. ASSOCIATED(cell_to_index))
      CALL get_cell(cell=cell, h=hmat, periodic=periodic)
      CALL get_kpoint_info(kpoints, nkp=n_k)
      n_spin = bs_env%n_spin
      CPASSERT(SIZE(mo_coeff, 1) >= n_k .AND. SIZE(mo_coeff, 2) == n_spin)
      valid_layout = SIZE(eigenvalues, 2) >= n_k .AND. SIZE(eigenvalues, 3) == n_spin
      CPASSERT(valid_layout)
      CPASSERT(SIZE(nocc_spin) == n_spin)
      ALLOCATE (nmo_spin(n_spin))
      DO i_spin = 1, n_spin
         CALL cp_cfm_get_info(mo_coeff(1, i_spin), nrow_global=n_row_matrix, ncol_global=n_mo_matrix)
         IF (i_spin == 1) n_ao = n_row_matrix
         valid_layout = n_row_matrix == n_ao .AND. n_mo_matrix == SIZE(eigenvalues, 1)
         CPASSERT(valid_layout)
         nmo_spin(i_spin) = n_mo_matrix
      END DO
      CPASSERT(ALL(nocc_spin >= 0 .AND. nocc_spin <= nmo_spin))

      ! Each rank receives all rows and images of its AO columns; no rank holds all images.
      n_img = SIZE(index_to_cell, 2)
      CALL gamma_corr_image_columns(gwkp_env%gamma_corr_s, gwkp_env%gamma_corr_scf_plan, cell_to_index, &
                                    n_ao, n_img, para_env, s_cols)
      CALL gamma_corr_image_columns(gwkp_env%gamma_corr_h, gwkp_env%gamma_corr_scf_plan, cell_to_index, &
                                    n_ao, n_img, para_env, h_cols)
      CALL gamma_corr_image_columns(gwkp_env%gamma_corr_d, gwkp_env%gamma_corr_plan, cell_to_index, &
                                    n_ao, n_img, para_env, d_cols)
      CALL dbcsr_deallocate_matrix_set(gwkp_env%gamma_corr_s)
      CALL dbcsr_deallocate_matrix_set(gwkp_env%gamma_corr_h)
      CALL dbcsr_deallocate_matrix_set(gwkp_env%gamma_corr_d)
      gwkp_env%gamma_corr_plan = kp_transform_plan_type()
      gwkp_env%gamma_corr_scf_plan = kp_transform_plan_type()
      CALL block_range(n_ao, para_env%mepos, para_env%num_pe, first_col, n_col)
      ALLOCATE (cell_cartesian(3, n_img))
      cell_cartesian(:, :) = MATMUL(hmat, REAL(index_to_cell, dp))
      ALLOCATE (phase_angle(n_img, n_k))
      phase_angle(:, :) = twopi*MATMUL(TRANSPOSE(REAL(index_to_cell, dp)), kpoints%xkp(:, 1:n_k))

      n_k_local = (n_k + para_env%num_pe - 1 - para_env%mepos)/para_env%num_pe
      ALLOCATE (dipole(n_spin, n_k_local, 3, MAXVAL(nocc_spin), MAXVAL(nmo_spin - nocc_spin)), source=z_zero)
      ! Each rank works only on its own k points, with their MO coefficients summed to it.
      CALL owned_mo_coefficients(bs_env, mo_coeff(1:n_k, :), coeff_owned)
      ! Per direction α: ∂S/∂k_α, the moment d_α and ∂H/∂k_α of each spin.
      n_mat = 3*(2 + n_spin)
      ALLOCATE (bloch(n_ao, n_col, n_mat, para_env%num_pe), operators_k(n_ao, n_ao, n_mat), &
                phase_cos(n_img, para_env%num_pe), phase_sin(n_img, para_env%num_pe), &
                weight_re(n_img, para_env%num_pe), weight_im(n_img, para_env%num_pe), &
                recv_buffer(n_ao*n_ao*n_mat), &
                send_counts(0:para_env%num_pe - 1), send_displs(0:para_env%num_pe - 1), &
                recv_counts(0:para_env%num_pe - 1), recv_displs(0:para_env%num_pe - 1))
      ! k points are owned cyclically: batch k_first..k_first+P-1 gives rank r the k point k_first+r.
      DO first_k = 1, n_k, para_env%num_pe
         n_batch = MIN(para_env%num_pe, n_k - first_k + 1)
         DO ib = 1, n_batch
            phase_cos(:, ib) = COS(phase_angle(:, first_k + ib - 1))
            phase_sin(:, ib) = SIN(phase_angle(:, first_k + ib - 1))
         END DO
         bloch(:, :, :, 1:n_batch) = z_zero
         DO i_dir = 1, 3
            IF (periodic(i_dir) == 0) CYCLE
            ! ∂X(k)/∂k_α = Σ_R i R_α e^{ikR} X(R)
            DO ib = 1, n_batch
               weight_re(:, ib) = -cell_cartesian(i_dir, :)*phase_sin(:, ib)
               weight_im(:, ib) = cell_cartesian(i_dir, :)*phase_cos(:, ib)
            END DO
            CALL bloch_sum_columns(s_cols(:, :, :, 1), weight_re, weight_im, n_batch, bloch(:, :, i_dir, :))
            CALL bloch_sum_columns(d_cols(:, :, :, i_dir), phase_cos, phase_sin, n_batch, &
                                   bloch(:, :, 3 + i_dir, :))
            DO i_spin = 1, n_spin
               CALL bloch_sum_columns(h_cols(:, :, :, i_spin), weight_re, weight_im, n_batch, &
                                      bloch(:, :, 3*(1 + i_spin) + i_dir, :))
            END DO
         END DO
         ! Transpose the column slices to the owner of each k point.
         DO i_rank = 0, para_env%num_pe - 1
            send_displs(i_rank) = n_ao*n_col*n_mat*i_rank
            send_counts(i_rank) = MERGE(n_ao*n_col*n_mat, 0, i_rank < n_batch)
            CALL block_range(n_ao, i_rank, para_env%num_pe, first_rank_col, n_rank_col)
            recv_displs(i_rank) = n_ao*(first_rank_col - 1)*n_mat
            recv_counts(i_rank) = MERGE(n_ao*n_rank_col*n_mat, 0, para_env%mepos < n_batch)
         END DO
         CALL para_env%alltoall(RESHAPE(bloch, [SIZE(bloch)]), send_counts, send_displs, &
                                recv_buffer, recv_counts, recv_displs)
         IF (para_env%mepos >= n_batch) CYCLE
         DO i_rank = 0, para_env%num_pe - 1
            CALL block_range(n_ao, i_rank, para_env%num_pe, first_rank_col, n_rank_col)
            operators_k(:, first_rank_col:first_rank_col + n_rank_col - 1, :) = &
               RESHAPE(recv_buffer(recv_displs(i_rank) + 1:recv_displs(i_rank) + recv_counts(i_rank)), &
                       [n_ao, n_rank_col, n_mat])
         END DO
         i_k = first_k + para_env%mepos
         i_k_local = (first_k - 1)/para_env%num_pe + 1
         DO i_spin = 1, n_spin
            n_mo = nmo_spin(i_spin)
            nocc = nocc_spin(i_spin)
            n_virt = n_mo - nocc
            IF (nocc == 0 .OR. n_virt == 0) CYCLE
            ALLOCATE (projected_h(n_virt, nocc), projected_s(n_virt, nocc), projected_d(n_virt, nocc))
            DO i_dir = 1, 3
               IF (periodic(i_dir) == 0) CYCLE
               CALL project_virtual_occupied(operators_k(:, :, 3*(1 + i_spin) + i_dir), &
                                             coeff_owned(:, 1:n_mo, i_k_local, i_spin), nocc, projected_h)
               CALL project_virtual_occupied(operators_k(:, :, i_dir), &
                                             coeff_owned(:, 1:n_mo, i_k_local, i_spin), nocc, projected_s)
               CALL project_virtual_occupied(operators_k(:, :, 3 + i_dir), &
                                             coeff_owned(:, 1:n_mo, i_k_local, i_spin), nocc, projected_d)
               DO i_mo_i = 1, nocc
                  DO i_mo_a = 1, n_virt
                     energy_diff = eigenvalues(nocc + i_mo_a, i_k, i_spin) - eigenvalues(i_mo_i, i_k, i_spin)
                     IF (ABS(energy_diff) <= eps_degenerate) CYCLE
                     dipole(i_spin, i_k_local, i_dir, i_mo_i, i_mo_a) = &
                        (-gaussi*projected_h(i_mo_a, i_mo_i) + &
                         gaussi*eigenvalues(nocc + i_mo_a, i_k, i_spin)* &
                         projected_s(i_mo_a, i_mo_i))/energy_diff + projected_d(i_mo_a, i_mo_i)
                  END DO
               END DO
            END DO
            DEALLOCATE (projected_h, projected_s, projected_d)
         END DO
      END DO
      DEALLOCATE (bloch, cell_cartesian, coeff_owned, d_cols, h_cols, nmo_spin, operators_k, phase_angle, phase_cos, &
                  phase_sin, recv_buffer, recv_counts, recv_displs, s_cols, send_counts, send_displs, &
                  weight_im, weight_re)
      DEALLOCATE (gwkp_env%gamma_corr_index_to_cell, gwkp_env%gamma_corr_cell_to_index)
      NULLIFY (gwkp_env%gamma_corr_index_to_cell, gwkp_env%gamma_corr_cell_to_index)

      CALL timestop(handle)

   END SUBROUTINE gw_gamma_transition_moments

! **************************************************************************************************
!> \brief Scatter one low-rank Gamma correction into the centred grid-pair images.
!> For s=nearest_integer[H^{-1}(r_l-r_l')], C_ll'(R)=C0_ll'(R-s).
!> R-s is wrapped into the finite FFT cell with its transfer-mesh Bloch twist. Half-cell
!> ties average adjacent integer images, retaining Hermiticity and cell-wrapping covariance.
!> The mirrored conjugate-transpose contribution enforces reciprocity also on Nyquist faces.
!> No grid-pair angular tensor is formed. values_R holds flattened column-major pairs
!> beginning at first_pair.
!> \param bs_env ...
!> \param cells ...
!> \param icell ...
!> \param head_scalar ...
!> \param head_vector ...
!> \param first_pair ...
!> \param values_R ...
!> \param wing_left ...
!> \param wing_right ...
! **************************************************************************************************
   SUBROUTINE add_gamma_corr_pair_cell(bs_env, cells, icell, head_scalar, head_vector, &
                                       first_pair, values_R, wing_left, wing_right)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: cells(:, :), icell
      COMPLEX(KIND=dp), INTENT(IN)                       :: head_scalar, head_vector(:)
      INTEGER, INTENT(IN)                                :: first_pair
      REAL(KIND=dp), INTENT(INOUT)                       :: values_R(:, :, :)
      COMPLEX(KIND=dp), INTENT(IN), OPTIONAL             :: wing_left(:), wing_right(:)

      COMPLEX(KIND=dp)                                   :: amplitude(2), twist(27, 2)
      INTEGER :: destination(27, 2), high(3), idir, iimage, il, ilocal, ilp, imirror, ipair, ix, &
         iy, iz, jlocal, low(3), lower(3), ngrid, nmesh(3), shift(3), translated(3), wrapped(3)
      REAL(KIND=dp)                                      :: delta(3), half_boundary, qoffset(3), &
                                                            tolerance, weight
      REAL(KIND=dp), ALLOCATABLE                         :: fractional_grid(:, :)

      ngrid = bs_env%ri_rs%n_grid_points
      nmesh = bs_env%kpoints_scf_desymm%nkp_grid
      lower = -nmesh/2
      qoffset = bs_env%kpoints_chi_eps_W%xkp(1:3, 1)
      ALLOCATE (fractional_grid(3, ngrid))
      fractional_grid(:, :) = MATMUL(bs_env%ri_rs%cell%h_inv, bs_env%ri_rs%grid_points)
      DO idir = 1, 3
         IF (bs_env%periodic(idir) == 0) fractional_grid(idir, :) = 0.0_dp
      END DO
      DO iz = -1, 1
         DO iy = -1, 1
            DO ix = -1, 1
               shift = [ix, iy, iz]
               iimage = 1 + ix + 1 + 3*(iy + 1) + 9*(iz + 1)
               DO imirror = 1, 2
                  translated = (3 - 2*imirror)*cells(:, icell) + shift
                  wrapped = MODULO(translated - lower, nmesh) + lower
                  destination(iimage, imirror) = 1 + (wrapped(3) - lower(3)) + &
                                                 nmesh(3)*(wrapped(2) - lower(2) + nmesh(2)*(wrapped(1) - lower(1)))
                  twist(iimage, imirror) = EXP(CMPLX(0.0_dp, &
                                                     twopi*DOT_PRODUCT(qoffset, REAL(translated - wrapped, dp)), dp))
               END DO
            END DO
         END DO
      END DO
      ! Every local pair row receives its own images, so threads never share an element.
      !$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(2) SCHEDULE(STATIC) &
      !$OMP SHARED(values_R, first_pair, ngrid, fractional_grid, head_scalar, head_vector, &
      !$OMP        wing_left, wing_right, destination, twist) &
      !$OMP PRIVATE(jlocal, ilocal, ipair, il, ilp, delta, low, high, idir, half_boundary, &
      !$OMP         tolerance, weight, amplitude, iz, iy, ix, iimage, imirror)
      DO jlocal = 1, SIZE(values_R, 2)
         DO ilocal = 1, SIZE(values_R, 1)
            ipair = first_pair + ilocal - 1
            il = MODULO(ipair - 1, ngrid) + 1
            ilp = (ipair - 1)/ngrid + 1
            delta = fractional_grid(:, il) - fractional_grid(:, ilp)
            low = NINT(delta)
            high = low
            DO idir = 1, 3
               half_boundary = ANINT(delta(idir) - 0.5_dp) + 0.5_dp
               tolerance = 64.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, ABS(delta(idir)))
               IF (ABS(delta(idir) - half_boundary) <= tolerance) THEN
                  low(idir) = NINT(half_boundary - 0.5_dp)
                  high(idir) = low(idir) + 1
               END IF
            END DO
            weight = 0.5_dp/REAL(PRODUCT(high - low + 1), dp)
            amplitude(1) = head_scalar*head_vector(il)*CONJG(head_vector(ilp))
            amplitude(2) = CONJG(head_scalar)*head_vector(il)*CONJG(head_vector(ilp))
            IF (PRESENT(wing_left)) THEN
               amplitude(1) = amplitude(1) + &
                              head_vector(il)*wing_right(ilp) + wing_left(il)*CONJG(head_vector(ilp))
               amplitude(2) = amplitude(2) + &
                              CONJG(head_vector(ilp)*wing_right(il) + wing_left(ilp)*CONJG(head_vector(il)))
            END IF
            DO iz = low(3), high(3)
               DO iy = low(2), high(2)
                  DO ix = low(1), high(1)
                     iimage = 1 + ix + 1 + 3*(iy + 1) + 9*(iz + 1)
                     DO imirror = 1, 2
                        values_R(ilocal, jlocal, destination(iimage, imirror)) = &
                           values_R(ilocal, jlocal, destination(iimage, imirror)) + &
                           weight*REAL(twist(iimage, imirror)*amplitude(imirror), dp)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE add_gamma_corr_pair_cell

! **************************************************************************************************
!> \brief Gather the sparse Gamma-correction AO images as column slices of their dense form.
!> Rank r receives all rows and images of its AO columns block_range(n_ao, r), so no rank holds
!> all images. Blocks keep the block and mirror conventions of replicate_rs_matrices.
!> \param matrices Distributed sparse AO image matrices, indexed (component, source image).
!> \param plan Prepared source image and atom-pair traversal.
!> \param cell_to_index Common Gamma-correction image map.
!> \param n_ao ...
!> \param n_img ...
!> \param para_env ...
!> \param values values(μ, ν, image, component) for the owned AO columns ν
! **************************************************************************************************
   SUBROUTINE gamma_corr_image_columns(matrices, plan, cell_to_index, n_ao, n_img, para_env, values)
      TYPE(dbcsr_p_type), POINTER                        :: matrices(:, :)
      TYPE(kp_transform_plan_type), INTENT(IN)           :: plan
      INTEGER, POINTER                                   :: cell_to_index(:, :, :)
      INTEGER, INTENT(IN)                                :: n_ao, n_img
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: values(:, :, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gamma_corr_image_columns'

      CHARACTER                                          :: symmetry
      INTEGER :: cell(3), first, first_col, handle, i_col, i_comp, i_entry, i_header, i_image, &
         i_mirror, i_pass, i_rank, i_row, i_source, i_swap, i_value, n_col, n_local, n_rank
      INTEGER, ALLOCATABLE :: col_owner(:), header_counts(:), header_displs(:), header_recv(:), &
         header_recv_counts(:), header_recv_displs(:), header_send(:), value_counts(:), &
         value_displs(:), value_recv_counts(:), value_recv_displs(:)
      INTEGER, POINTER                                   :: col_offsets(:), row_offsets(:)
      LOGICAL                                            :: found
      REAL(dp)                                           :: mirror_sign
      REAL(dp), ALLOCATABLE                              :: value_recv(:), value_send(:)
      REAL(dp), POINTER                                  :: block(:, :)

      CALL timeset(routineN, handle)

      n_rank = para_env%num_pe
      ALLOCATE (col_owner(n_ao))
      DO i_rank = 0, n_rank - 1
         CALL block_range(n_ao, i_rank, n_rank, first, n_local)
         col_owner(first:first + n_local - 1) = i_rank
      END DO
      CALL block_range(n_ao, para_env%mepos, n_rank, first_col, n_col)
      CALL dbcsr_get_info(matrices(1, 1)%matrix, matrix_type=symmetry, &
                          row_blk_offset=row_offsets, col_blk_offset=col_offsets)
      mirror_sign = MERGE(-1.0_dp, 1.0_dp, symmetry == dbcsr_type_antisymmetric)
      ALLOCATE (header_counts(0:n_rank - 1), header_displs(0:n_rank - 1), &
                value_counts(0:n_rank - 1), value_displs(0:n_rank - 1), &
                header_recv_counts(0:n_rank - 1), header_recv_displs(0:n_rank - 1), &
                value_recv_counts(0:n_rank - 1), value_recv_displs(0:n_rank - 1))
      ! The first pass counts the pieces for each column owner, the second packs them.
      DO i_pass = 1, 2
         IF (i_pass == 2) THEN
            header_displs(0) = 0
            value_displs(0) = 0
            DO i_rank = 1, n_rank - 1
               header_displs(i_rank) = header_displs(i_rank - 1) + header_counts(i_rank - 1)
               value_displs(i_rank) = value_displs(i_rank - 1) + value_counts(i_rank - 1)
            END DO
            ALLOCATE (header_send(MAX(1, SUM(header_counts))), value_send(MAX(1, SUM(value_counts))))
         END IF
         header_counts(:) = 0
         value_counts(:) = 0
         DO i_comp = 1, SIZE(matrices, 1)
            DO i_entry = 1, plan%nentry
               i_source = plan%image(i_entry)
               i_row = plan%row(i_entry)
               i_col = plan%col(i_entry)
               cell = plan%cell(:, i_entry)
               IF (plan%symmetry_sign(i_entry) < 0.0_dp) THEN
                  i_swap = i_row
                  i_row = i_col
                  i_col = i_swap
               END IF
               IF (symmetry /= dbcsr_type_no_symmetry .AND. i_row > i_col) THEN
                  i_swap = i_row
                  i_row = i_col
                  i_col = i_swap
                  cell = -cell
               END IF
               CALL dbcsr_get_readonly_block_p(matrices(i_comp, i_source)%matrix, row=i_row, col=i_col, &
                                               block=block, found=found)
               IF (.NOT. found) CYCLE
               i_image = cell_to_index(cell(1), cell(2), cell(3))
               CPASSERT(i_image > 0)
               CALL pack_image_piece(block, [i_comp, i_image, row_offsets(i_row)], col_offsets(i_col), &
                                     col_owner, i_pass == 2, header_counts, value_counts, header_displs, &
                                     value_displs, header_send, value_send)
               IF (symmetry /= dbcsr_type_no_symmetry .AND. i_row /= i_col) THEN
                  i_mirror = cell_to_index(-cell(1), -cell(2), -cell(3))
                  CPASSERT(i_mirror > 0)
                  CALL pack_image_piece(mirror_sign*TRANSPOSE(block), [i_comp, i_mirror, col_offsets(i_col)], &
                                        row_offsets(i_row), col_owner, i_pass == 2, header_counts, &
                                        value_counts, header_displs, value_displs, header_send, value_send)
               END IF
            END DO
         END DO
      END DO
      CALL para_env%alltoall(header_counts, header_recv_counts, 1)
      CALL para_env%alltoall(value_counts, value_recv_counts, 1)
      header_recv_displs(0) = 0
      value_recv_displs(0) = 0
      DO i_rank = 1, n_rank - 1
         header_recv_displs(i_rank) = header_recv_displs(i_rank - 1) + header_recv_counts(i_rank - 1)
         value_recv_displs(i_rank) = value_recv_displs(i_rank - 1) + value_recv_counts(i_rank - 1)
      END DO
      ALLOCATE (header_recv(MAX(1, SUM(header_recv_counts))), value_recv(MAX(1, SUM(value_recv_counts))))
      CALL para_env%alltoall(header_send, header_counts, header_displs, header_recv, header_recv_counts, &
                             header_recv_displs)
      CALL para_env%alltoall(value_send, value_counts, value_displs, value_recv, value_recv_counts, &
                             value_recv_displs)
      DEALLOCATE (header_send, value_send)

      ! Header (component, image, first row, rows, first column, columns); values column-major.
      ALLOCATE (values(n_ao, n_col, n_img, SIZE(matrices, 1)), SOURCE=0.0_dp)
      i_value = 0
      DO i_header = 1, SUM(header_recv_counts), 6
         ASSOCIATE (h => header_recv(i_header:i_header + 5))
            values(h(3):h(3) + h(4) - 1, h(5) - first_col + 1:h(5) - first_col + h(6), h(2), h(1)) = &
               RESHAPE(value_recv(i_value + 1:i_value + h(4)*h(6)), [h(4), h(6)])
            i_value = i_value + h(4)*h(6)
         END ASSOCIATE
      END DO
      DEALLOCATE (col_owner, header_counts, header_displs, header_recv, header_recv_counts, &
                  header_recv_displs, value_counts, value_displs, value_recv, value_recv_counts, &
                  value_recv_displs)

      CALL timestop(handle)

   END SUBROUTINE gamma_corr_image_columns

! **************************************************************************************************
!> \brief Count or pack one dense AO-image piece for the owners of its columns.
!> \param piece block values, rows varying fastest
!> \param label component, image and first AO row of the piece
!> \param first_col first AO column of the piece
!> \param col_owner rank owning each AO column
!> \param fill pack into the buffers; otherwise only count
!> \param header_counts ...
!> \param value_counts ...
!> \param header_displs ...
!> \param value_displs ...
!> \param header_send ...
!> \param value_send ...
! **************************************************************************************************
   SUBROUTINE pack_image_piece(piece, label, first_col, col_owner, fill, header_counts, value_counts, &
                               header_displs, value_displs, header_send, value_send)
      REAL(dp), INTENT(IN)                               :: piece(:, :)
      INTEGER, INTENT(IN)                                :: label(3), first_col, col_owner(:)
      LOGICAL, INTENT(IN)                                :: fill
      INTEGER, INTENT(INOUT)                             :: header_counts(0:), value_counts(0:)
      INTEGER, INTENT(IN)                                :: header_displs(0:), value_displs(0:)
      INTEGER, ALLOCATABLE, INTENT(INOUT)                :: header_send(:)
      REAL(dp), ALLOCATABLE, INTENT(INOUT)               :: value_send(:)

      INTEGER                                            :: i_first, n_cols, n_rows, owner

      n_rows = SIZE(piece, 1)
      i_first = 1
      DO WHILE (i_first <= SIZE(piece, 2))
         owner = col_owner(first_col + i_first - 1)
         n_cols = 1
         DO WHILE (i_first + n_cols <= SIZE(piece, 2))
            IF (col_owner(first_col + i_first + n_cols - 1) /= owner) EXIT
            n_cols = n_cols + 1
         END DO
         IF (fill) THEN
            header_send(header_displs(owner) + header_counts(owner) + 1: &
                        header_displs(owner) + header_counts(owner) + 6) = &
               [label, n_rows, first_col + i_first - 1, n_cols]
            value_send(value_displs(owner) + value_counts(owner) + 1: &
                       value_displs(owner) + value_counts(owner) + n_rows*n_cols) = &
               RESHAPE(piece(:, i_first:i_first + n_cols - 1), [n_rows*n_cols])
         END IF
         header_counts(owner) = header_counts(owner) + 6
         value_counts(owner) = value_counts(owner) + n_rows*n_cols
         i_first = i_first + n_cols
      END DO
   END SUBROUTINE pack_image_piece

! **************************************************************************************************
!> \brief Bloch sums of AO-image column slices for a batch of k points with one GEMM each,
!>        X(k) = Σ_R (w_re(R,k) + i w_im(R,k)) X(R).
!> \param images X(R), indexed (row, owned column, image)
!> \param weight_re ...
!> \param weight_im ...
!> \param n_batch number of k points
!> \param bloch X(k), indexed (row, owned column, k)
! **************************************************************************************************
   SUBROUTINE bloch_sum_columns(images, weight_re, weight_im, n_batch, bloch)
      REAL(dp), INTENT(IN)                               :: images(:, :, :), weight_re(:, :), &
                                                            weight_im(:, :)
      INTEGER, INTENT(IN)                                :: n_batch
      COMPLEX(dp), INTENT(INOUT)                         :: bloch(:, :, :)

      INTEGER                                            :: n_img, n_values
      REAL(dp), ALLOCATABLE                              :: sum_im(:, :), sum_re(:, :)

      n_values = SIZE(images, 1)*SIZE(images, 2)
      n_img = SIZE(images, 3)
      IF (n_values == 0) RETURN
      ALLOCATE (sum_re(n_values, n_batch), sum_im(n_values, n_batch))
      CALL dgemm('N', 'N', n_values, n_batch, n_img, 1.0_dp, images, n_values, weight_re, SIZE(weight_re, 1), &
                 0.0_dp, sum_re, n_values)
      CALL dgemm('N', 'N', n_values, n_batch, n_img, 1.0_dp, images, n_values, weight_im, SIZE(weight_im, 1), &
                 0.0_dp, sum_im, n_values)
      bloch(:, :, 1:n_batch) = RESHAPE(CMPLX(sum_re, sum_im, KIND=dp), [SIZE(images, 1), SIZE(images, 2), n_batch])
      DEALLOCATE (sum_re, sum_im)
   END SUBROUTINE bloch_sum_columns

! **************************************************************************************************
!> \brief Project an AO operator between virtual and occupied MOs, C_virt^H X C_occ.
!> \param operator_ao ...
!> \param coeff MO coefficients, occupied first
!> \param nocc ...
!> \param projected ...
! **************************************************************************************************
   SUBROUTINE project_virtual_occupied(operator_ao, coeff, nocc, projected)
      COMPLEX(dp), INTENT(IN)                            :: operator_ao(:, :), coeff(:, :)
      INTEGER, INTENT(IN)                                :: nocc
      COMPLEX(dp), INTENT(OUT)                           :: projected(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: work(:, :)
      INTEGER                                            :: n_ao, n_virt

      n_ao = SIZE(operator_ao, 1)
      n_virt = SIZE(coeff, 2) - nocc
      ALLOCATE (work(n_ao, nocc))
      CALL zgemm('N', 'N', n_ao, nocc, n_ao, z_one, operator_ao, n_ao, coeff, SIZE(coeff, 1), &
                 z_zero, work, n_ao)
      CALL zgemm('C', 'N', n_virt, nocc, n_ao, z_one, coeff(:, nocc + 1:), SIZE(coeff, 1), work, n_ao, &
                 z_zero, projected, n_virt)
      DEALLOCATE (work)
   END SUBROUTINE project_virtual_occupied

! **************************************************************************************************
!> \brief Collect SCF and position-operator image cells, including mirrored images.
!> \param kp_a ...
!> \param kp_b ...
!> \param index_to_cell ...
!> \param cell_to_index ...
! **************************************************************************************************
   SUBROUTINE gw_gamma_cell_index_union(kp_a, kp_b, index_to_cell, cell_to_index)
      TYPE(kpoint_type), POINTER                         :: kp_a, kp_b
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index

      INTEGER                                            :: i, ic, isgn, n
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: list
      INTEGER, DIMENSION(3)                              :: image, max_image
      INTEGER, DIMENSION(:, :), POINTER                  :: itc_a, itc_b

      NULLIFY (itc_a, itc_b)
      CALL get_kpoint_info(kp_a, index_to_cell=itc_a)
      CALL get_kpoint_info(kp_b, index_to_cell=itc_b)
      CPASSERT(ASSOCIATED(itc_a) .AND. ASSOCIATED(itc_b))
      DO i = 1, 3
         max_image(i) = MAX(MAXVAL(ABS(itc_a(i, :))), MAXVAL(ABS(itc_b(i, :))))
      END DO
      ALLOCATE (cell_to_index(-max_image(1):max_image(1), -max_image(2):max_image(2), &
                              -max_image(3):max_image(3)), source=0)
      ALLOCATE (list(3, 2*(SIZE(itc_a, 2) + SIZE(itc_b, 2))))
      n = 0
      DO isgn = 1, -1, -2
         DO ic = 1, SIZE(itc_a, 2)
            image = isgn*itc_a(:, ic)
            IF (cell_to_index(image(1), image(2), image(3)) /= 0) CYCLE
            n = n + 1
            cell_to_index(image(1), image(2), image(3)) = n
            list(:, n) = image
         END DO
      END DO
      DO isgn = 1, -1, -2
         DO ic = 1, SIZE(itc_b, 2)
            image = isgn*itc_b(:, ic)
            IF (cell_to_index(image(1), image(2), image(3)) /= 0) CYCLE
            n = n + 1
            cell_to_index(image(1), image(2), image(3)) = n
            list(:, n) = image
         END DO
      END DO
      ALLOCATE (index_to_cell(3, n))
      index_to_cell = list(:, 1:n)
      DEALLOCATE (list)
   END SUBROUTINE gw_gamma_cell_index_union

END MODULE gw_ri_rs_small_cell_full_kp_gamma

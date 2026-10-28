!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Complex k-point orbital transformations, response kernels and local batch execution.
!>        Shared real/spectral operations are provided by qs_ot.
! **************************************************************************************************
MODULE qs_ot_kp
   USE cp_cfm_diag,                     ONLY: cp_cfm_diag_workspace_type,&
                                              cp_cfm_heevd_local,&
                                              cp_cfm_local_workspace_prepare
   USE cp_cfm_types,                    ONLY: cp_cfm_p_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_distribution_type, dbcsr_filter, &
        dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_reserve_blocks, dbcsr_scale, dbcsr_set, &
        dbcsr_transposed, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_add_on_diag,&
                                              dbcsr_frobenius_norm,&
                                              dbcsr_hadamard_product,&
                                              dbcsr_reserve_all_blocks,&
                                              dbcsr_scale_by_vector
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_heevd
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_local,&
                                              copy_local_to_dbcsr
   USE cp_fm_diag,                      ONLY: FM_DIAG_TYPE_DLAF,&
                                              diag_lib_explicit,&
                                              diag_type
   USE cp_fm_types,                     ONLY: cp_fm_p_type
   USE kinds,                           ONLY: dp
   USE local_gemm_api,                  ONLY: LOCAL_GEMM_PU_GPU,&
                                              local_gemm_ctxt_type
   USE mathlib,                         ONLY: diag_complex,&
                                              diamat_all,&
                                              gemm_square
   USE preconditioner,                  ONLY: apply_preconditioner,&
                                              apply_preconditioner_complex_batch
   USE preconditioner_types,            ONLY: preconditioner_p_type
   USE qs_ot,                           ONLY: qs_ot_complex_multiply,&
                                              qs_ot_p2m_diag,&
                                              qs_ot_project_gradient_real => qs_ot_project_gradient,&
                                              qs_ot_ref_update,&
                                              qs_ot_spectral_kernel,&
                                              qs_ot_spectral_values,&
                                              qs_ot_square_transpose
   USE qs_ot_types,                     ONLY: qs_ot_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qs_ot_get_p_complex
   PUBLIC :: qs_ot_get_orbitals_complex
   PUBLIC :: qs_ot_get_orbitals_local_batch
   PUBLIC :: qs_ot_get_derivative_complex
   PUBLIC :: qs_ot_get_derivative_local_batch
   PUBLIC :: qs_ot_prepare_complex_tangent_metric
   PUBLIC :: qs_ot_prepare_complex_tangent_metric_batch
   PUBLIC :: qs_ot_get_orbitals_ref_complex
   PUBLIC :: qs_ot_get_derivative_ref_complex
   PUBLIC :: qs_ot_apply_complex_frechet_dbcsr
   PUBLIC :: qs_ot_antihermitian_spectral_norm
   PUBLIC :: qs_ot_complex_exp_frechet_kernel
   PUBLIC :: qs_ot_density_secant_orbital_overlaps
   PUBLIC :: qs_ot_density_secant_projected_hessian
   PUBLIC :: qs_ot_density_tangent
   PUBLIC :: qs_ot_fixed_n_energy_gradient
   PUBLIC :: qs_ot_fixed_n_multigroup_schur_block
   PUBLIC :: qs_ot_fixed_n_response_mu_shift
   PUBLIC :: qs_ot_finite_rotation_response
   PUBLIC :: qs_ot_projected_response_update
   PUBLIC :: qs_ot_symmetric_abs_solve
   PUBLIC :: qs_ot_generate_rotation_complex
   PUBLIC :: qs_ot_rot_mat_derivative_complex
   PUBLIC :: qs_ot_project_gradient

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot_kp'

   INTERFACE qs_ot_complex_multiply_local
      MODULE PROCEDURE qs_ot_complex_multiply_pair_local, qs_ot_complex_multiply_parts_local
   END INTERFACE

   ! Transient per-worker storage. The last index is 1 for real and 2 for imaginary.
   TYPE qs_ot_local_work_type
      INTEGER                                            :: n = 0, k = 0, info = 0
      LOGICAL                                            :: filtered = .FALSE.
      TYPE(local_gemm_ctxt_type)                         :: gemm_ctx
      TYPE(cp_cfm_diag_workspace_type)                   :: solver
      COMPLEX(dp), ALLOCATABLE                           :: vectors(:, :)
      REAL(dp), ALLOCATABLE                              :: values(:)
      REAL(dp), ALLOCATABLE :: hc(:, :, :), x(:, :, :), sx(:, :, :), c0(:, :, :), out(:, :, :), &
                               r(:, :, :), sinp(:, :, :), sc0(:, :, :), TARGET(:, :, :), os(:, :, :), &
                               buf1(:, :, :), buf2(:, :, :), buf3(:, :, :), buf4(:, :, :), buf_nk(:, :, :), &
                               buf1_ortho(:, :, :), buf2_ortho(:, :, :), sinp_b(:, :), cosp_b(:, :), tmp(:, :), tmp_ortho(:, :)
   END TYPE qs_ot_local_work_type

CONTAINS

! **************************************************************************************************
!> \brief spectral norm of a dense anti-Hermitian rotation generator
!> \param rotation_generator anti-Hermitian generator
!> \return largest absolute eigenvalue of i times the generator
! **************************************************************************************************
   FUNCTION qs_ot_antihermitian_spectral_norm(rotation_generator) RESULT(norm)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: rotation_generator
      REAL(KIND=dp)                                      :: norm

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: eigenvectors
      INTEGER                                            :: n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues

      n = SIZE(rotation_generator, 1)
      CPASSERT(SIZE(rotation_generator, 2) == n)
      norm = 0.0_dp
      IF (n == 0) RETURN
      ALLOCATE (eigenvectors(n, n), eigenvalues(n))
      CALL diag_complex(CMPLX(0.0_dp, 1.0_dp, KIND=dp)*rotation_generator, &
                        eigenvectors, eigenvalues)
      norm = MAXVAL(ABS(eigenvalues))
      DEALLOCATE (eigenvalues, eigenvectors)

   END FUNCTION qs_ot_antihermitian_spectral_norm

! **************************************************************************************************
!> \brief chemical-potential response for one fixed-electron-number group
!> \param weighted_energy_response sum_i chi_i de_i over the perturbed local channels
!> \param local_curvature_sum local sum_i chi_i, used as a serial fallback
!> \param fixed_n_curvature_sum global sum_i chi_i for the complete fixed-N group
!> \return first-order chemical-potential shift
! **************************************************************************************************
   PURE FUNCTION qs_ot_fixed_n_response_mu_shift( &
      weighted_energy_response, local_curvature_sum, fixed_n_curvature_sum) RESULT(mu_shift)
      REAL(KIND=dp), INTENT(IN)                          :: weighted_energy_response, &
                                                            local_curvature_sum, &
                                                            fixed_n_curvature_sum
      REAL(KIND=dp)                                      :: mu_shift

      REAL(KIND=dp)                                      :: denominator

      denominator = fixed_n_curvature_sum
      IF (ABS(denominator) <= EPSILON(denominator)) denominator = local_curvature_sum
      mu_shift = 0.0_dp
      IF (ABS(denominator) > EPSILON(denominator)) mu_shift = weighted_energy_response/denominator

   END FUNCTION qs_ot_fixed_n_response_mu_shift

! **************************************************************************************************
!> \brief fixed-N Mermin gradient in auxiliary-energy coordinates
!> \param rayleigh_energy diagonal expectation values of the current Hamiltonian
!> \param energy_coordinate auxiliary band energies controlling the occupations
!> \param response_weight signed weighted occupation responses chi_i
!> \param fixed_n_weight_sum global sum_i chi_i for the fixed-N group
!> \param fixed_n_weighted_residual global sum_i chi_i (h_i-e_i)
!> \param gradient projected fixed-N gradient
! **************************************************************************************************
   PURE SUBROUTINE qs_ot_fixed_n_energy_gradient( &
      rayleigh_energy, energy_coordinate, response_weight, fixed_n_weight_sum, &
      fixed_n_weighted_residual, gradient)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rayleigh_energy, energy_coordinate, &
                                                            response_weight
      REAL(KIND=dp), INTENT(IN)                          :: fixed_n_weight_sum, &
                                                            fixed_n_weighted_residual
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: gradient

      REAL(KIND=dp)                                      :: fixed_n_mean

      fixed_n_mean = 0.0_dp
      IF (ABS(fixed_n_weight_sum) > EPSILON(fixed_n_weight_sum)) THEN
         fixed_n_mean = fixed_n_weighted_residual/fixed_n_weight_sum
      END IF
      gradient(:) = response_weight(:)* &
                    (fixed_n_mean - (rayleigh_energy(:) - energy_coordinate(:)))

   END SUBROUTINE qs_ot_fixed_n_energy_gradient

! **************************************************************************************************
!> \brief eliminate spin-resolved auxiliary energies while retaining every fixed-N constraint
!>
!>        Each energy variable belongs to one particle-number group. The diagonal occupation
!>        response is eliminated locally, while one coupling column per group is retained for the
!>        subsequent global chemical-potential projection. A single group reduces exactly to
!>        the single-group Schur complement A - R^T (D - chi chi^T/sum(chi)) R.
!> \param rotation_hessian finite orbital-rotation Hessian
!> \param rayleigh_response derivative of all spin-resolved Rayleigh energies
!> \param response_weight positive occupation-response weights
!> \param response_group fixed-N group of every energy variable
!> \param rotation_gradient orbital-rotation gradient
!> \param energy_gradient auxiliary-energy gradient
!> \param schur_block energy-eliminated rotation block
!> \param coupling_matrix retained chemical-potential coupling columns
!> \param schur_rhs energy-eliminated rotation right-hand side
! **************************************************************************************************
   SUBROUTINE qs_ot_fixed_n_multigroup_schur_block( &
      rotation_hessian, rayleigh_response, response_weight, response_group, &
      rotation_gradient, energy_gradient, schur_block, coupling_matrix, schur_rhs)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: rotation_hessian, rayleigh_response
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: response_weight
      INTEGER, DIMENSION(:), INTENT(IN)                  :: response_group
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rotation_gradient, energy_gradient
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: schur_block, coupling_matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: schur_rhs

      INTEGER                                            :: group, i, nenergy, ngroups, nrotation
      REAL(KIND=dp)                                      :: weight

      nenergy = SIZE(response_weight)
      nrotation = SIZE(rotation_gradient)
      ngroups = SIZE(coupling_matrix, 2)
      CPASSERT(ngroups > 0)
      CPASSERT(ALL(SHAPE(rotation_hessian) == [nrotation, nrotation]))
      CPASSERT(ALL(SHAPE(rayleigh_response) == [nenergy, nrotation]))
      CPASSERT(SIZE(energy_gradient) == nenergy)
      CPASSERT(SIZE(response_group) == nenergy)
      CPASSERT(ALL(response_group >= 1 .AND. response_group <= ngroups))
      CPASSERT(ALL(SHAPE(schur_block) == [nrotation, nrotation]))
      CPASSERT(SIZE(coupling_matrix, 1) == nrotation)
      CPASSERT(SIZE(schur_rhs) == nrotation)

      schur_block(:, :) = rotation_hessian(:, :)
      coupling_matrix(:, :) = 0.0_dp
      DO i = 1, nenergy
         weight = response_weight(i)
         schur_block(:, :) = schur_block(:, :) - weight* &
                             SPREAD(rayleigh_response(i, :), DIM=2, NCOPIES=nrotation)* &
                             SPREAD(rayleigh_response(i, :), DIM=1, NCOPIES=nrotation)
         group = response_group(i)
         coupling_matrix(:, group) = coupling_matrix(:, group) + &
                                     weight*rayleigh_response(i, :)
      END DO
      schur_block(:, :) = 0.5_dp*(schur_block + TRANSPOSE(schur_block))
      schur_rhs(:) = rotation_gradient + &
                     MATMUL(TRANSPOSE(rayleigh_response), energy_gradient)

   END SUBROUTINE qs_ot_fixed_n_multigroup_schur_block

! **************************************************************************************************
!> \brief apply a positive spectral inverse of a real symmetric response matrix
!>
!>        The magnitude of every resolved eigenmode is retained, including modes with negative
!>        physical curvature. Replacing lambda by abs(lambda) gives a descent metric without the
!>        loss of response information caused by discarding the negative subspace. Unresolved
!>        null modes are projected out instead of being amplified by an artificial eigenvalue floor.
!> \param matrix real symmetric response matrix
!> \param rhs one or more right-hand sides
!> \param solution spectral-absolute inverse applied to rhs
!> \param valid whether finite input and output were obtained
!> \param relative_floor optional relative eigenvalue floor
! **************************************************************************************************
   SUBROUTINE qs_ot_symmetric_abs_solve(matrix, rhs, solution, valid, relative_floor)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: matrix, rhs
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: solution
      LOGICAL, INTENT(OUT)                               :: valid
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: relative_floor

      INTEGER                                            :: i, n, nresolved
      REAL(KIND=dp)                                      :: eigenvalue_floor, relative_floor_eff, &
                                                            scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenvectors, work

      n = SIZE(matrix, 1)
      CPASSERT(ALL(SHAPE(matrix) == [n, n]))
      CPASSERT(SIZE(rhs, 1) == n)
      CPASSERT(ALL(SHAPE(solution) == SHAPE(rhs)))
      valid = ALL(matrix == matrix) .AND. ALL(rhs == rhs)
      IF (.NOT. valid) THEN
         solution(:, :) = 0.0_dp
         RETURN
      END IF

      relative_floor_eff = SQRT(EPSILON(1.0_dp))
      IF (PRESENT(relative_floor)) relative_floor_eff = MAX(relative_floor, relative_floor_eff)
      ALLOCATE (eigenvalues(n), eigenvectors(n, n), work(n, SIZE(rhs, 2)))
      eigenvectors(:, :) = 0.5_dp*(matrix + TRANSPOSE(matrix))
      CALL diamat_all(eigenvectors, eigenvalues)
      scale = MAX(1.0_dp, MAXVAL(ABS(eigenvalues)))
      eigenvalue_floor = relative_floor_eff*scale
      work(:, :) = MATMUL(TRANSPOSE(eigenvectors), rhs)
      nresolved = 0
      DO i = 1, n
         IF (ABS(eigenvalues(i)) > eigenvalue_floor) THEN
            work(i, :) = work(i, :)/ABS(eigenvalues(i))
            nresolved = nresolved + 1
         ELSE
            work(i, :) = 0.0_dp
         END IF
      END DO
      solution(:, :) = MATMUL(eigenvectors, work)
      valid = nresolved > 0 .AND. ALL(solution == solution)
      IF (.NOT. valid) solution(:, :) = 0.0_dp
      DEALLOCATE (eigenvalues, eigenvectors, work)

   END SUBROUTINE qs_ot_symmetric_abs_solve

! **************************************************************************************************
!> \brief update a baseline response direction in a small positive physical-response subspace
!>
!>        The first basis mode is the baseline response direction. With B0 the projected frozen-H
!>        Hessian and K the accepted physical correction, this routine solves
!>
!>          (B0 + K) c = g_Q.
!>
!>        The optional projected physical gradient supplies g_Q. Without it, g_Q=B0*e1, so c=e1
!>        exactly when K=0. An indefinite or unresolved total projected Hessian is rejected instead
!>        of turning its negative modes into an unrelated active direction.
!> \param reference_hessian projected frozen-H Hessian B0
!> \param response_correction projected physical response K
!> \param coefficients response coefficients c in the supplied basis
!> \param valid whether a finite, positive, sufficiently resolved solve was obtained
!> \param projected_gradient optional physical gradient projected onto the supplied basis
!> \param relative_floor optional relative positive-eigenvalue floor
! **************************************************************************************************
   SUBROUTINE qs_ot_projected_response_update( &
      reference_hessian, response_correction, coefficients, valid, projected_gradient, relative_floor)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: reference_hessian, response_correction
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: coefficients
      LOGICAL, INTENT(OUT)                               :: valid
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: projected_gradient
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: relative_floor

      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: eigenvalue_floor, relative_floor_eff, &
                                                            scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, rhs, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenvectors

      n = SIZE(reference_hessian, 1)
      CPASSERT(n > 0)
      CPASSERT(ALL(SHAPE(reference_hessian) == [n, n]))
      CPASSERT(ALL(SHAPE(response_correction) == [n, n]))
      CPASSERT(SIZE(coefficients) == n)
      IF (PRESENT(projected_gradient)) THEN
         CPASSERT(SIZE(projected_gradient) == n)
      END IF
      coefficients(:) = 0.0_dp
      valid = ALL(reference_hessian == reference_hessian) .AND. &
              ALL(response_correction == response_correction)
      IF (PRESENT(projected_gradient)) THEN
         valid = valid .AND. ALL(projected_gradient == projected_gradient)
      END IF
      IF (.NOT. valid) RETURN

      relative_floor_eff = 1.0E-4_dp
      IF (PRESENT(relative_floor)) relative_floor_eff = &
         MAX(relative_floor, SQRT(EPSILON(1.0_dp)))
      ALLOCATE (eigenvalues(n), eigenvectors(n, n), rhs(n), work(n))
      eigenvectors(:, :) = 0.5_dp* &
                           (reference_hessian + response_correction + &
                            TRANSPOSE(reference_hessian + response_correction))
      CALL diamat_all(eigenvectors, eigenvalues)
      scale = MAX(MAXVAL(ABS(eigenvalues)), MAXVAL(ABS(reference_hessian)))
      IF (scale <= TINY(scale)) THEN
         valid = .FALSE.
         coefficients(:) = 0.0_dp
         DEALLOCATE (eigenvalues, eigenvectors, rhs, work)
         RETURN
      END IF
      eigenvalue_floor = relative_floor_eff*scale
      valid = ALL(eigenvalues > eigenvalue_floor)
      IF (valid) THEN
         rhs(:) = reference_hessian(:, 1)
         IF (PRESENT(projected_gradient)) rhs(:) = projected_gradient
         work(:) = MATMUL(TRANSPOSE(eigenvectors), rhs)
         DO i = 1, n
            work(i) = work(i)/eigenvalues(i)
         END DO
         coefficients(:) = MATMUL(eigenvectors, work)
         valid = ALL(coefficients == coefficients) .AND. &
                 ALL(ABS(coefficients) <= HUGE(1.0_dp))
      END IF
      IF (.NOT. valid) coefficients(:) = 0.0_dp
      DEALLOCATE (eigenvalues, eigenvectors, rhs, work)

   END SUBROUTINE qs_ot_projected_response_update

! **************************************************************************************************
!> \brief form a projected self-adjoint Hxc response from distributed density-space overlaps
!> \param density_norm_sq <S,S>
!> \param response_work <S,Y>
!> \param density_overlap <S,B_r>
!> \param response_overlap <Y,B_r>
!> \param correction projected symmetric Hxc response <B_q,K*B_r>
!> \param valid whether finite nonzero secant data were available
!> \param secant_mode optional mode representing the accepted full density secant divided by its
!>        line-search position
!> \param secant_position signed line-search position of the accepted full density secant
! **************************************************************************************************
   SUBROUTINE qs_ot_density_secant_projected_hessian( &
      density_norm_sq, response_work, density_overlap, response_overlap, correction, valid, &
      secant_mode, secant_position)

      REAL(KIND=dp), INTENT(IN)                          :: density_norm_sq, response_work
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: density_overlap, response_overlap
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: correction
      LOGICAL, INTENT(OUT)                               :: valid
      INTEGER, INTENT(IN), OPTIONAL                      :: secant_mode
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: secant_position

      INTEGER                                            :: i, j, nmode
      REAL(KIND=dp)                                      :: inverse_density_norm
      REAL(KIND=dp), DIMENSION(SIZE(density_overlap))    :: projected_density_overlap, &
                                                            projected_response_overlap

      nmode = SIZE(density_overlap)
      CPASSERT(SIZE(response_overlap) == nmode)
      CPASSERT(ALL(SHAPE(correction) == [nmode, nmode]))

      correction(:, :) = 0.0_dp
      valid = .FALSE.
      IF (nmode <= 0 .OR. density_norm_sq <= TINY(1.0_dp)) RETURN
      IF (density_norm_sq /= density_norm_sq .OR. response_work /= response_work) RETURN
      IF (ABS(density_norm_sq) > HUGE(1.0_dp) .OR. ABS(response_work) > HUGE(1.0_dp)) RETURN
      IF (ANY(density_overlap /= density_overlap) .OR. ANY(response_overlap /= response_overlap)) RETURN
      IF (ANY(ABS(density_overlap) > HUGE(1.0_dp)) .OR. &
          ANY(ABS(response_overlap) > HUGE(1.0_dp))) RETURN
      IF (PRESENT(secant_mode) .NEQV. PRESENT(secant_position)) RETURN

      projected_density_overlap(:) = density_overlap
      projected_response_overlap(:) = response_overlap
      IF (PRESENT(secant_mode)) THEN
         IF (secant_mode < 1 .OR. secant_mode > nmode) RETURN
         IF (secant_position /= secant_position .OR. &
             ABS(secant_position) > HUGE(1.0_dp) .OR. &
             ABS(secant_position) <= SQRT(EPSILON(1.0_dp))) RETURN
         ! The complete accepted density change can contain REF-orbital motion that is absent from
         ! a reduced rotation/occupation tangent.  S/alpha supplies its exact projected secant mode.
         projected_density_overlap(secant_mode) = density_norm_sq/secant_position
         projected_response_overlap(secant_mode) = response_work/secant_position
      END IF

      inverse_density_norm = 1.0_dp/density_norm_sq
      DO j = 1, nmode
         DO i = 1, nmode
            correction(i, j) = inverse_density_norm* &
                               (projected_response_overlap(i)*projected_density_overlap(j) + &
                                projected_density_overlap(i)*projected_response_overlap(j) - &
                                response_work*inverse_density_norm* &
                                projected_density_overlap(i)*projected_density_overlap(j))
         END DO
      END DO
      correction(:, :) = 0.5_dp*(correction + TRANSPOSE(correction))
      valid = ALL(correction == correction) .AND. ALL(ABS(correction) <= HUGE(1.0_dp))
      IF (.NOT. valid) correction(:, :) = 0.0_dp

   END SUBROUTINE qs_ot_density_secant_projected_hessian

! **************************************************************************************************
!> \brief finite-chart density tangent for coupled complex rotations and fixed-N occupations
!>
!>        The rotation contribution differentiates
!>
!>          exp(X) diag(w_k f) exp(X)^H
!>
!>        along an anti-Hermitian packed direction.  The supplied weighted occupation response is
!>        added in the same chart, and the result is returned in the current physical orbital basis.
!> \param rotation_generator current anti-Hermitian REF generator X
!> \param occupation current occupations f
!> \param kpoint_weight irreducible K-point weight w_k
!> \param rotation_step interleaved real/imaginary anti-Hermitian direction
!> \param weighted_occupation_step derivative of w_k*f, including the fixed-N mu response
!> \param density_tangent Hermitian tangent in the current physical orbital basis
!> \param difference_step optional finite-chart central-difference step
! **************************************************************************************************
   SUBROUTINE qs_ot_density_tangent( &
      rotation_generator, occupation, kpoint_weight, rotation_step, &
      weighted_occupation_step, density_tangent, difference_step)

      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: rotation_generator
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occupation
      REAL(KIND=dp), INTENT(IN)                          :: kpoint_weight
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rotation_step, weighted_occupation_step
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: density_tangent
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: difference_step

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: direction, generator_minus, &
         generator_plus, mode_ref, rotation, rotation_minus, rotation_plus, weighted_rotation
      INTEGER                                            :: i, j, n, nrotation, r
      REAL(KIND=dp)                                      :: step, step_scale

      n = SIZE(rotation_generator, 1)
      nrotation = n*(n - 1)
      CPASSERT(n > 0)
      CPASSERT(ALL(SHAPE(rotation_generator) == [n, n]))
      CPASSERT(SIZE(occupation) == n)
      CPASSERT(SIZE(rotation_step) == nrotation)
      CPASSERT(SIZE(weighted_occupation_step) == n)
      CPASSERT(ALL(SHAPE(density_tangent) == [n, n]))
      CPASSERT(kpoint_weight > 0.0_dp)

      ALLOCATE (direction(n, n), generator_minus(n, n), generator_plus(n, n), &
                mode_ref(n, n), rotation(n, n), rotation_minus(n, n), &
                rotation_plus(n, n), weighted_rotation(n, n))
      direction(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      r = 0
      DO i = 1, n - 1
         DO j = i + 1, n
            r = r + 1
            direction(i, j) = CMPLX(rotation_step(r), 0.0_dp, KIND=dp)
            direction(j, i) = -direction(i, j)
            r = r + 1
            direction(i, j) = direction(i, j) + &
                              CMPLX(0.0_dp, rotation_step(r), KIND=dp)
            direction(j, i) = direction(j, i) + &
                              CMPLX(0.0_dp, rotation_step(r), KIND=dp)
         END DO
      END DO
      CPASSERT(r == nrotation)

      step_scale = MAX(1.0_dp, MAXVAL(ABS(direction)))
      step = 1.0E-5_dp/step_scale
      IF (PRESENT(difference_step)) step = difference_step/step_scale
      CPASSERT(step > EPSILON(step))
      generator_plus(:, :) = rotation_generator + step*direction
      generator_minus(:, :) = rotation_generator - step*direction
      CALL qs_ot_dense_rotation_state(rotation_generator, rotation)
      CALL qs_ot_dense_rotation_state(generator_plus, rotation_plus)
      CALL qs_ot_dense_rotation_state(generator_minus, rotation_minus)

      weighted_rotation(:, :) = rotation_plus
      DO j = 1, n
         weighted_rotation(:, j) = kpoint_weight*occupation(j)*weighted_rotation(:, j)
      END DO
      mode_ref(:, :) = MATMUL(weighted_rotation, CONJG(TRANSPOSE(rotation_plus)))
      weighted_rotation(:, :) = rotation_minus
      DO j = 1, n
         weighted_rotation(:, j) = kpoint_weight*occupation(j)*weighted_rotation(:, j)
      END DO
      mode_ref(:, :) = (mode_ref - &
                        MATMUL(weighted_rotation, CONJG(TRANSPOSE(rotation_minus))))/(2.0_dp*step)

      weighted_rotation(:, :) = rotation
      DO j = 1, n
         weighted_rotation(:, j) = weighted_occupation_step(j)*weighted_rotation(:, j)
      END DO
      mode_ref(:, :) = mode_ref + MATMUL(weighted_rotation, CONJG(TRANSPOSE(rotation)))
      density_tangent(:, :) = MATMUL(CONJG(TRANSPOSE(rotation)), MATMUL(mode_ref, rotation))
      density_tangent(:, :) = 0.5_dp*(density_tangent + CONJG(TRANSPOSE(density_tangent)))

      DEALLOCATE (direction, generator_minus, generator_plus, mode_ref, rotation, &
                  rotation_minus, rotation_plus, weighted_rotation)

   END SUBROUTINE qs_ot_density_tangent

! **************************************************************************************************
!> \brief project a physical density/Hamiltonian secant between moving orbital subspaces
!>
!>        For separately S-orthonormal endpoint orbitals C0 and C1, O=C0^H*S*C1 retains the
!>        component of the accepted density step that leaves the old subspace.  Density modes are
!>        represented in the current C1 basis and already contain the irreducible K-point weight.
!> \param overlap_start_current cross overlap O
!> \param occupation_start occupations at the accepted start
!> \param occupation_current occupations at the accepted endpoint
!> \param hamiltonian_step_start C0^H*(H1-H0)*C0
!> \param hamiltonian_step_current C1^H*(H1-H0)*C1
!> \param density_modes current-orbital density tangents
!> \param kpoint_weight irreducible K-point weight
!> \param density_norm_sq contribution to <Delta P,Delta P>
!> \param response_work contribution to <Delta P,Delta H>
!> \param density_overlap contributions <Delta P,B_r>
!> \param response_overlap contributions <Delta H,B_r>
!> \param valid whether a finite nonzero secant was available
! **************************************************************************************************
   SUBROUTINE qs_ot_density_secant_orbital_overlaps( &
      overlap_start_current, occupation_start, occupation_current, hamiltonian_step_start, &
      hamiltonian_step_current, density_modes, kpoint_weight, density_norm_sq, response_work, &
      density_overlap, response_overlap, valid)

      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: overlap_start_current
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occupation_start, occupation_current
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: hamiltonian_step_start, &
                                                            hamiltonian_step_current
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: density_modes
      REAL(KIND=dp), INTENT(IN)                          :: kpoint_weight
      REAL(KIND=dp), INTENT(OUT)                         :: density_norm_sq, response_work
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: density_overlap, response_overlap
      LOGICAL, INTENT(OUT)                               :: valid

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: cross_mode
      INTEGER                                            :: i, j, mode, n, nmode
      REAL(KIND=dp)                                      :: cross_density, density_scale

      n = SIZE(occupation_start)
      nmode = SIZE(density_modes, 3)
      CPASSERT(n > 0)
      CPASSERT(SIZE(occupation_current) == n)
      CPASSERT(ALL(SHAPE(overlap_start_current) == [n, n]))
      CPASSERT(ALL(SHAPE(hamiltonian_step_start) == [n, n]))
      CPASSERT(ALL(SHAPE(hamiltonian_step_current) == [n, n]))
      CPASSERT(SIZE(density_modes, 1) == n .AND. SIZE(density_modes, 2) == n)
      CPASSERT(SIZE(density_overlap) == nmode .AND. SIZE(response_overlap) == nmode)

      density_norm_sq = 0.0_dp
      response_work = 0.0_dp
      density_overlap(:) = 0.0_dp
      response_overlap(:) = 0.0_dp
      valid = .FALSE.
      IF (nmode <= 0 .OR. kpoint_weight <= TINY(1.0_dp)) RETURN

      cross_density = 0.0_dp
      DO j = 1, n
         DO i = 1, n
            cross_density = cross_density + occupation_start(i)*occupation_current(j)* &
                            ABS(overlap_start_current(i, j))**2
         END DO
      END DO
      density_norm_sq = kpoint_weight* &
                        (SUM(occupation_start**2) + SUM(occupation_current**2) - &
                         2.0_dp*cross_density)
      response_work = 0.0_dp
      DO i = 1, n
         response_work = response_work + kpoint_weight* &
                         (occupation_current(i)*REAL(hamiltonian_step_current(i, i), KIND=dp) - &
                          occupation_start(i)*REAL(hamiltonian_step_start(i, i), KIND=dp))
      END DO

      ALLOCATE (cross_mode(n, n))
      DO mode = 1, nmode
         ! GCC 15 miscompiles the equivalent complex array expressions with -march=znver3.
         CALL gemm_square(overlap_start_current, "N", density_modes(:, :, mode), "N", &
                          overlap_start_current, "C", cross_mode)
         DO i = 1, n
            density_overlap(mode) = density_overlap(mode) + &
                                    occupation_current(i)* &
                                    REAL(density_modes(i, i, mode), KIND=dp) - &
                                    occupation_start(i)*REAL(cross_mode(i, i), KIND=dp)
         END DO
         response_overlap(mode) = &
            SUM(REAL(hamiltonian_step_current, KIND=dp)* &
                REAL(density_modes(:, :, mode), KIND=dp) + &
                AIMAG(hamiltonian_step_current)*AIMAG(density_modes(:, :, mode)))
      END DO
      DEALLOCATE (cross_mode)

      density_scale = kpoint_weight* &
                      MAX(MAXVAL(ABS(occupation_start)), MAXVAL(ABS(occupation_current)))
      IF (density_scale <= TINY(1.0_dp)) RETURN
      IF (density_norm_sq <= 64.0_dp*EPSILON(1.0_dp)*density_scale**2) RETURN
      valid = density_norm_sq == density_norm_sq .AND. response_work == response_work .AND. &
              ABS(density_norm_sq) <= HUGE(1.0_dp) .AND. ABS(response_work) <= HUGE(1.0_dp) .AND. &
              ALL(density_overlap == density_overlap) .AND. &
              ALL(response_overlap == response_overlap) .AND. &
              ALL(ABS(density_overlap) <= HUGE(1.0_dp)) .AND. &
              ALL(ABS(response_overlap) <= HUGE(1.0_dp))
      IF (.NOT. valid) THEN
         density_norm_sq = 0.0_dp
         response_work = 0.0_dp
         density_overlap(:) = 0.0_dp
         response_overlap(:) = 0.0_dp
      END IF

   END SUBROUTINE qs_ot_density_secant_orbital_overlaps

! **************************************************************************************************
!> \brief finite complex REF rotation Hessian and Rayleigh-energy response
!>
!>        The current projected Hamiltonian is pulled back through the finite rotation and then
!>        differentiated in the independent real-antisymmetric and imaginary-symmetric pair
!>        coordinates. This keeps the response consistent with the exponential chart used by REF
!>        OT instead of replacing it by an infinitesimal commutator away from the chart origin.
!> \param chc current projected Hermitian Hamiltonian U^H H_ref U
!> \param rotation_generator current anti-Hermitian REF generator
!> \param occupation fixed occupations attached to the rotated columns
!> \param kpoint_weight irreducible-k-point weight
!> \param rotation_gradient gradient in interleaved real/imaginary pair coordinates
!> \param rotation_hessian derivative of rotation_gradient in the same coordinates
!> \param rayleigh_response derivative of diag(U^H H_ref U) with respect to the pair coordinates
!> \param difference_step optional central finite-difference step for the Hessian action
! **************************************************************************************************
   SUBROUTINE qs_ot_finite_rotation_response( &
      chc, rotation_generator, occupation, kpoint_weight, rotation_gradient, &
      rotation_hessian, rayleigh_response, difference_step)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: chc, rotation_generator
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occupation
      REAL(KIND=dp), INTENT(IN)                          :: kpoint_weight
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: rotation_gradient
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: rotation_hessian, rayleigh_response
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: difference_step

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: base_hamiltonian, generator_minus, &
                                                            generator_plus, rotation
      INTEGER                                            :: i, j, n, nrotation, r, s
      REAL(KIND=dp)                                      :: step
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: occupation_scale, rayleigh, &
                                                            rayleigh_minus, rayleigh_plus
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: gradient_imag, gradient_minus_imag, &
         gradient_minus_real, gradient_plus_imag, gradient_plus_real, gradient_real

      n = SIZE(chc, 1)
      nrotation = n*(n - 1)
      CPASSERT(n > 0)
      CPASSERT(ALL(SHAPE(chc) == [n, n]))
      CPASSERT(ALL(SHAPE(rotation_generator) == [n, n]))
      CPASSERT(SIZE(occupation) == n)
      CPASSERT(SIZE(rotation_gradient) == nrotation)
      CPASSERT(ALL(SHAPE(rotation_hessian) == [nrotation, nrotation]))
      CPASSERT(ALL(SHAPE(rayleigh_response) == [n, nrotation]))
      CPASSERT(kpoint_weight > 0.0_dp)

      step = 1.0E-4_dp
      IF (PRESENT(difference_step)) step = difference_step
      CPASSERT(step > SQRT(EPSILON(step)))

      ALLOCATE (base_hamiltonian(n, n), generator_minus(n, n), generator_plus(n, n), &
                rotation(n, n), occupation_scale(n), rayleigh(n), rayleigh_minus(n), &
                rayleigh_plus(n), gradient_imag(n, n), gradient_minus_imag(n, n), &
                gradient_minus_real(n, n), gradient_plus_imag(n, n), &
                gradient_plus_real(n, n), gradient_real(n, n))

      CALL qs_ot_dense_rotation_state(rotation_generator, rotation)
      base_hamiltonian(:, :) = MATMUL(rotation, MATMUL(chc, CONJG(TRANSPOSE(rotation))))
      occupation_scale(:) = 2.0_dp*kpoint_weight*occupation(:)
      CALL qs_ot_dense_rotation_gradient(rotation_generator, base_hamiltonian, occupation_scale, &
                                         gradient_real, gradient_imag, rayleigh)

      r = 0
      DO i = 1, n - 1
         DO j = i + 1, n
            r = r + 1
            rotation_gradient(r) = gradient_real(i, j)
            r = r + 1
            rotation_gradient(r) = gradient_imag(i, j)
         END DO
      END DO
      CPASSERT(r == nrotation)

      s = 0
      DO i = 1, n - 1
         DO j = i + 1, n
            generator_plus(:, :) = rotation_generator(:, :)
            generator_minus(:, :) = rotation_generator(:, :)
            generator_plus(i, j) = generator_plus(i, j) + step
            generator_plus(j, i) = generator_plus(j, i) - step
            generator_minus(i, j) = generator_minus(i, j) - step
            generator_minus(j, i) = generator_minus(j, i) + step
            s = s + 1
            CALL qs_ot_dense_rotation_gradient(generator_plus, base_hamiltonian, occupation_scale, &
                                               gradient_plus_real, gradient_plus_imag, rayleigh_plus)
            CALL qs_ot_dense_rotation_gradient(generator_minus, base_hamiltonian, occupation_scale, &
                                               gradient_minus_real, gradient_minus_imag, rayleigh_minus)
            CALL qs_ot_pack_rotation_response( &
               gradient_plus_real, gradient_plus_imag, gradient_minus_real, gradient_minus_imag, &
               rayleigh_plus, rayleigh_minus, step, rotation_hessian(:, s), &
               rayleigh_response(:, s))

            generator_plus(:, :) = rotation_generator(:, :)
            generator_minus(:, :) = rotation_generator(:, :)
            generator_plus(i, j) = generator_plus(i, j) + CMPLX(0.0_dp, step, KIND=dp)
            generator_plus(j, i) = generator_plus(j, i) + CMPLX(0.0_dp, step, KIND=dp)
            generator_minus(i, j) = generator_minus(i, j) - CMPLX(0.0_dp, step, KIND=dp)
            generator_minus(j, i) = generator_minus(j, i) - CMPLX(0.0_dp, step, KIND=dp)
            s = s + 1
            CALL qs_ot_dense_rotation_gradient(generator_plus, base_hamiltonian, occupation_scale, &
                                               gradient_plus_real, gradient_plus_imag, rayleigh_plus)
            CALL qs_ot_dense_rotation_gradient(generator_minus, base_hamiltonian, occupation_scale, &
                                               gradient_minus_real, gradient_minus_imag, rayleigh_minus)
            CALL qs_ot_pack_rotation_response( &
               gradient_plus_real, gradient_plus_imag, gradient_minus_real, gradient_minus_imag, &
               rayleigh_plus, rayleigh_minus, step, rotation_hessian(:, s), &
               rayleigh_response(:, s))
         END DO
      END DO
      CPASSERT(s == nrotation)
      rotation_hessian(:, :) = 0.5_dp*(rotation_hessian + TRANSPOSE(rotation_hessian))

      DEALLOCATE (base_hamiltonian, generator_minus, generator_plus, rotation, occupation_scale, &
                  rayleigh, rayleigh_minus, rayleigh_plus, gradient_imag, gradient_minus_imag, &
                  gradient_minus_real, gradient_plus_imag, gradient_plus_real, gradient_real)

   END SUBROUTINE qs_ot_finite_rotation_response

! **************************************************************************************************
!> \brief pack one finite complex rotation response column
!> \param gradient_plus_real real gradient at the positive endpoint
!> \param gradient_plus_imag imaginary gradient at the positive endpoint
!> \param gradient_minus_real real gradient at the negative endpoint
!> \param gradient_minus_imag imaginary gradient at the negative endpoint
!> \param rayleigh_plus Rayleigh energies at the positive endpoint
!> \param rayleigh_minus Rayleigh energies at the negative endpoint
!> \param step central finite-difference step
!> \param hessian_column packed rotation-Hessian column
!> \param rayleigh_column packed Rayleigh-response column
! **************************************************************************************************
   SUBROUTINE qs_ot_pack_rotation_response( &
      gradient_plus_real, gradient_plus_imag, gradient_minus_real, gradient_minus_imag, &
      rayleigh_plus, rayleigh_minus, step, hessian_column, rayleigh_column)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: gradient_plus_real, gradient_plus_imag, &
                                                            gradient_minus_real, &
                                                            gradient_minus_imag
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rayleigh_plus, rayleigh_minus
      REAL(KIND=dp), INTENT(IN)                          :: step
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: hessian_column, rayleigh_column

      INTEGER                                            :: i, j, n, r

      n = SIZE(gradient_plus_real, 1)
      CPASSERT(ALL(SHAPE(gradient_plus_real) == [n, n]))
      CPASSERT(ALL(SHAPE(gradient_plus_imag) == [n, n]))
      CPASSERT(ALL(SHAPE(gradient_minus_real) == [n, n]))
      CPASSERT(ALL(SHAPE(gradient_minus_imag) == [n, n]))
      CPASSERT(SIZE(hessian_column) == n*(n - 1))
      CPASSERT(SIZE(rayleigh_column) == n)

      r = 0
      DO i = 1, n - 1
         DO j = i + 1, n
            r = r + 1
            hessian_column(r) = &
               (gradient_plus_real(i, j) - gradient_minus_real(i, j))/(2.0_dp*step)
            r = r + 1
            hessian_column(r) = &
               (gradient_plus_imag(i, j) - gradient_minus_imag(i, j))/(2.0_dp*step)
         END DO
      END DO
      rayleigh_column(:) = (rayleigh_plus(:) - rayleigh_minus(:))/(2.0_dp*step)

   END SUBROUTINE qs_ot_pack_rotation_response

! **************************************************************************************************
!> \brief dense complex rotation and projected-Hamiltonian diagonal
!> \param rotation_generator anti-Hermitian REF generator
!> \param base_hamiltonian fixed Hamiltonian in the unrotated REF basis
!> \param occupation_scale twice the weighted occupations
!> \param gradient_real real-antisymmetric gradient component
!> \param gradient_imag imaginary-symmetric gradient component
!> \param rayleigh diagonal of the rotated Hamiltonian
! **************************************************************************************************
   SUBROUTINE qs_ot_dense_rotation_gradient( &
      rotation_generator, base_hamiltonian, occupation_scale, gradient_real, gradient_imag, rayleigh)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: rotation_generator, base_hamiltonian
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occupation_scale
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: gradient_real, gradient_imag
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: rayleigh

      COMPLEX(KIND=dp)                                   :: kernel
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: covector, eigenvectors, frechet, inner, &
                                                            outer, rotation, work
      INTEGER                                            :: i, j, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues

      n = SIZE(rotation_generator, 1)
      ALLOCATE (covector(n, n), eigenvectors(n, n), frechet(n, n), inner(n, n), &
                outer(n, n), rotation(n, n), work(n, n), eigenvalues(n))
      CALL qs_ot_dense_rotation_state(rotation_generator, rotation, &
                                      eigenvectors=eigenvectors, eigenvalues=eigenvalues)

      work(:, :) = MATMUL(base_hamiltonian, rotation)
      covector(:, :) = work(:, :)
      DO j = 1, n
         covector(:, j) = occupation_scale(j)*covector(:, j)
      END DO

      work(:, :) = MATMUL(CONJG(TRANSPOSE(rotation)), &
                          MATMUL(base_hamiltonian, rotation))
      DO i = 1, n
         rayleigh(i) = REAL(work(i, i), KIND=dp)
      END DO

      inner(:, :) = MATMUL(CONJG(TRANSPOSE(eigenvectors)), &
                           MATMUL(covector, eigenvectors))
      DO j = 1, n
         DO i = 1, n
            kernel = CONJG(qs_ot_complex_exp_frechet_kernel(eigenvalues(i), eigenvalues(j)))
            outer(i, j) = inner(i, j)*kernel
         END DO
      END DO
      frechet(:, :) = MATMUL(eigenvectors, &
                             MATMUL(outer, CONJG(TRANSPOSE(eigenvectors))))
      gradient_real(:, :) = REAL(frechet, KIND=dp) - TRANSPOSE(REAL(frechet, KIND=dp))
      gradient_imag(:, :) = AIMAG(frechet) + TRANSPOSE(AIMAG(frechet))

      DEALLOCATE (covector, eigenvectors, frechet, inner, outer, rotation, work, eigenvalues)

   END SUBROUTINE qs_ot_dense_rotation_gradient

! **************************************************************************************************
!> \brief dense exponential of an anti-Hermitian REF generator
!> \param rotation_generator anti-Hermitian generator
!> \param rotation exp(rotation_generator)
!> \param eigenvectors optional eigenvectors of i*rotation_generator
!> \param eigenvalues optional eigenvalues of i*rotation_generator
! **************************************************************************************************
   SUBROUTINE qs_ot_dense_rotation_state(rotation_generator, rotation, eigenvectors, eigenvalues)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: rotation_generator
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: rotation
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT), OPTIONAL :: eigenvalues

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: vectors, weighted_vectors
      INTEGER                                            :: j, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: values

      n = SIZE(rotation_generator, 1)
      ALLOCATE (vectors(n, n), weighted_vectors(n, n), values(n))
      CALL diag_complex(CMPLX(0.0_dp, 1.0_dp, KIND=dp)*rotation_generator, &
                        vectors, values)
      weighted_vectors(:, :) = vectors(:, :)
      DO j = 1, n
         weighted_vectors(:, j) = EXP(CMPLX(0.0_dp, -values(j), KIND=dp))* &
                                  weighted_vectors(:, j)
      END DO
      rotation(:, :) = MATMUL(weighted_vectors, CONJG(TRANSPOSE(vectors)))
      IF (PRESENT(eigenvectors)) eigenvectors(:, :) = vectors(:, :)
      IF (PRESENT(eigenvalues)) eigenvalues(:) = values(:)
      DEALLOCATE (vectors, weighted_vectors, values)

   END SUBROUTINE qs_ot_dense_rotation_state

! **************************************************************************************************
!> \brief Frechet divided-difference kernel for exp(-i*evals)
!> \param e1 ...
!> \param e2 ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION qs_ot_complex_exp_frechet_kernel(e1, e2) RESULT(kernel)
      REAL(KIND=dp), INTENT(IN)                          :: e1, e2
      COMPLEX(KIND=dp)                                   :: kernel

      COMPLEX(KIND=dp)                                   :: l1, l2, x
      INTEGER                                            :: i

      l1 = (0.0_dp, -1.0_dp)*e1
      l2 = (0.0_dp, -1.0_dp)*e2
      IF (ABS(l1 - l2) > 0.5_dp) THEN
         kernel = (EXP(l1) - EXP(l2))/(l1 - l2)
      ELSE
         x = 1.0_dp
         kernel = 0.0_dp
         DO i = 1, 16
            kernel = kernel + x
            x = x*(l1 - l2)/REAL(i + 1, KIND=dp)
         END DO
         kernel = kernel*EXP(l2)
      END IF

   END FUNCTION qs_ot_complex_exp_frechet_kernel

! **************************************************************************************************
!> \brief apply the complex exponential Frechet kernel to sparse DBCSR Re/Im matrices
!> \param evals generator eigenvalues
!> \param inner_deriv_re real part of the matrix in the generator eigenbasis
!> \param inner_deriv_im imaginary part of the matrix in the generator eigenbasis
!> \param outer_deriv_re real part of the mapped matrix
!> \param outer_deriv_im imaginary part of the mapped matrix
!> \param adjoint use the adjoint Frechet kernel for gradients
! **************************************************************************************************
   SUBROUTINE qs_ot_apply_complex_frechet_dbcsr(evals, inner_deriv_re, inner_deriv_im, &
                                                outer_deriv_re, outer_deriv_im, adjoint)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: evals
      TYPE(dbcsr_type)                                   :: inner_deriv_re, inner_deriv_im, &
                                                            outer_deriv_re, outer_deriv_im
      LOGICAL, INTENT(IN), OPTIONAL                      :: adjoint

      COMPLEX(KIND=dp)                                   :: cval_in, kernel
      INTEGER                                            :: col, i, j, max_blocks, nblocks, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cols, rows
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_offset, col_blk_size, &
                                                            row_blk_offset, row_blk_size
      LOGICAL                                            :: found_im, found_out_im, found_out_re, &
                                                            found_re, use_adjoint
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_in_im, block_in_re, block_out_im, &
                                                            block_out_re
      REAL(KIND=dp)                                      :: e1, e2, im_part, re_part
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_iterator_type)                          :: iter

      use_adjoint = .FALSE.
      IF (PRESENT(adjoint)) use_adjoint = adjoint

      ! Re/Im parts can have different sparse block patterns. Build their union
      ! explicitly so a missing partner block is interpreted as zero.
      max_blocks = 0
      CALL dbcsr_iterator_start(iter, inner_deriv_re)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col)
         max_blocks = max_blocks + 1
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_iterator_start(iter, inner_deriv_im)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col)
         max_blocks = max_blocks + 1
      END DO
      CALL dbcsr_iterator_stop(iter)
      ALLOCATE (rows(MAX(max_blocks, 1)), cols(MAX(max_blocks, 1)))
      nblocks = 0

      CALL dbcsr_iterator_start(iter, inner_deriv_re)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col)
         CALL append_union_block(row, col)
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_iterator_start(iter, inner_deriv_im)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col)
         CALL append_union_block(row, col)
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_get_info(inner_deriv_re, distribution=dist, row_blk_size=row_blk_size, &
                          col_blk_size=col_blk_size)
      CALL dbcsr_create(outer_deriv_re, "outer_deriv_re", dist, dbcsr_type_no_symmetry, &
                        row_blk_size, col_blk_size)
      CALL dbcsr_create(outer_deriv_im, "outer_deriv_im", dist, dbcsr_type_no_symmetry, &
                        row_blk_size, col_blk_size)
      IF (nblocks > 0) THEN
         CALL dbcsr_reserve_blocks(outer_deriv_re, rows=rows(1:nblocks), cols=cols(1:nblocks))
         CALL dbcsr_reserve_blocks(outer_deriv_im, rows=rows(1:nblocks), cols=cols(1:nblocks))
      END IF
      CALL dbcsr_finalize(outer_deriv_re)
      CALL dbcsr_finalize(outer_deriv_im)
      CALL dbcsr_set(outer_deriv_re, 0.0_dp)
      CALL dbcsr_set(outer_deriv_im, 0.0_dp)

      CALL dbcsr_get_info(outer_deriv_re, row_blk_offset=row_blk_offset, &
                          col_blk_offset=col_blk_offset)
      CALL dbcsr_iterator_start(iter, outer_deriv_re)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col)
         CALL dbcsr_get_block_p(inner_deriv_re, row, col, block_in_re, found_re)
         CALL dbcsr_get_block_p(inner_deriv_im, row, col, block_in_im, found_im)
         CALL dbcsr_get_block_p(outer_deriv_re, row, col, block_out_re, found_out_re)
         CALL dbcsr_get_block_p(outer_deriv_im, row, col, block_out_im, found_out_im)
         CPASSERT(found_out_re .AND. found_out_im)

         DO i = 1, SIZE(block_out_re, 1)
         DO j = 1, SIZE(block_out_re, 2)
            e1 = evals(row_blk_offset(row) + i - 1)
            e2 = evals(col_blk_offset(col) + j - 1)
            re_part = 0.0_dp
            im_part = 0.0_dp
            IF (found_re) re_part = block_in_re(i, j)
            IF (found_im) im_part = block_in_im(i, j)
            cval_in = CMPLX(re_part, im_part, dp)
            kernel = qs_ot_complex_exp_frechet_kernel(e1, e2)
            IF (use_adjoint) kernel = CONJG(kernel)
            cval_in = cval_in*kernel
            block_out_re(i, j) = REAL(cval_in, KIND=dp)
            block_out_im(i, j) = AIMAG(cval_in)
         END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)
      DEALLOCATE (rows, cols)

   CONTAINS

! **************************************************************************************************
!> \brief append a block coordinate unless it is already present
!> \param row_new block-row index
!> \param col_new block-column index
! **************************************************************************************************
      SUBROUTINE append_union_block(row_new, col_new)
      INTEGER, INTENT(IN)                                :: row_new, col_new

      INTEGER                                            :: iblock

         DO iblock = 1, nblocks
            IF (rows(iblock) == row_new .AND. cols(iblock) == col_new) RETURN
         END DO
         nblocks = nblocks + 1
         rows(nblocks) = row_new
         cols(nblocks) = col_new
      END SUBROUTINE append_union_block

   END SUBROUTINE qs_ot_apply_complex_frechet_dbcsr

! **************************************************************************************************
!> \brief update complex REF k-point orbitals and their S(k)C(k) images
!> \param matrix_c ...
!> \param matrix_c_im ...
!> \param matrix_s ...
!> \param matrix_s_im ...
!> \param qs_ot_env ...
!> \param qs_ot_env1 environment carrying the shared minimizer state
! **************************************************************************************************
   SUBROUTINE qs_ot_get_orbitals_ref_complex(matrix_c, matrix_c_im, matrix_s, matrix_s_im, &
                                             qs_ot_env, qs_ot_env1)

      TYPE(dbcsr_type), POINTER                          :: matrix_c, matrix_c_im, matrix_s, &
                                                            matrix_s_im
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(qs_ot_type), OPTIONAL                         :: qs_ot_env1

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_orbitals_ref_complex'

      INTEGER                                            :: handle, i, k, n
      LOGICAL                                            :: update
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, inverse_sqrt
      TYPE(dbcsr_type)                                   :: rotated_im, rotated_re, rotation_tmp
      TYPE(dbcsr_type), POINTER                          :: c_im, c_re, f_im, f_re, p_im, p_re, &
                                                            sc_im, sc_re, tmp_kk, tmp_nk, v_im, &
                                                            v_re, w_im, w_re

      CALL timeset(routineN, handle)

      CPASSERT(qs_ot_env%has_complex_kpoint_state)
      CPASSERT(ASSOCIATED(matrix_s))
      CPASSERT(ASSOCIATED(matrix_s_im))

      c_re => qs_ot_env%matrix_x
      c_im => qs_ot_env%matrix_x_im
      f_re => qs_ot_env%matrix_ref_inv_sqrt
      f_im => qs_ot_env%matrix_ref_inv_sqrt_im
      sc_re => qs_ot_env%matrix_sx
      sc_im => qs_ot_env%matrix_sx_im
      p_re => qs_ot_env%buf1_k_k_sym
      p_im => qs_ot_env%buf2_k_k_sym
      v_re => qs_ot_env%buf3_k_k_sym
      v_im => qs_ot_env%buf4_k_k_sym
      w_re => qs_ot_env%buf1_k_k_nosym
      w_im => qs_ot_env%buf2_k_k_nosym
      tmp_kk => qs_ot_env%buf3_k_k_nosym
      tmp_nk => qs_ot_env%buf1_n_k

      CALL dbcsr_get_info(c_re, nfullrows_total=n, nfullcols_total=k)
      IF (PRESENT(qs_ot_env1)) THEN
         update = qs_ot_ref_update(qs_ot_env1)
      ELSE
         update = qs_ot_ref_update(qs_ot_env)
      END IF

      ! SC = (S_re + i*S_im) * (C_re + i*C_im)
      CALL qs_ot_complex_multiply('N', 'N', matrix_s, matrix_s_im, c_re, c_im, sc_re, sc_im, tmp_nk)

      ! P = C^H*S*C.  Its imaginary component is real antisymmetric.
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_re, sc_re, 0.0_dp, p_re)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_im, sc_im, 0.0_dp, tmp_kk)
      CALL dbcsr_add(p_re, tmp_kk, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_re, sc_im, 0.0_dp, p_im)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, c_im, sc_re, 0.0_dp, tmp_kk)
      CALL dbcsr_add(p_im, tmp_kk, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      ALLOCATE (eigenvalues(k), inverse_sqrt(k))
      CALL cp_dbcsr_heevd(matrix_re=p_re, matrix_im=p_im, &
                          eigenvectors_re=v_re, eigenvectors_im=v_im, &
                          eigenvalues=eigenvalues, para_env=qs_ot_env%para_env, &
                          blacs_env=qs_ot_env%blacs_env)
      DO i = 1, k
         IF (eigenvalues(i) <= EPSILON(1.0_dp)) THEN
            CPABORT("Complex REF overlap is not positive definite")
         END IF
         inverse_sqrt(i) = 1.0_dp/SQRT(eigenvalues(i))
      END DO

      ! P^(-1/2) = V*diag(lambda^(-1/2))*V^H.
      CALL dbcsr_copy(w_re, v_re)
      CALL dbcsr_copy(w_im, v_im)
      CALL dbcsr_scale_by_vector(w_re, alpha=inverse_sqrt, side='right')
      CALL dbcsr_scale_by_vector(w_im, alpha=inverse_sqrt, side='right')

      CALL dbcsr_multiply('N', 'T', 1.0_dp, w_re, v_re, 0.0_dp, p_re)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, w_im, v_im, 0.0_dp, tmp_kk)
      CALL dbcsr_add(p_re, tmp_kk, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      CALL dbcsr_multiply('N', 'T', 1.0_dp, w_im, v_re, 0.0_dp, p_im)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, w_re, v_im, 0.0_dp, tmp_kk)
      CALL dbcsr_add(p_im, tmp_kk, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      CALL dbcsr_copy(f_re, p_re)
      CALL dbcsr_copy(f_im, p_im)

      ! Return the physical, orthonormal orbitals without changing a rejected
      ! line-search coordinate.
      CALL qs_ot_complex_multiply('N', 'N', c_re, c_im, p_re, p_im, matrix_c, matrix_c_im, tmp_nk)

      IF (update) THEN
         CALL dbcsr_copy(c_re, matrix_c)
         CALL dbcsr_copy(c_im, matrix_c_im)

         CALL dbcsr_set(f_re, 0.0_dp)
         CALL dbcsr_add_on_diag(f_re, alpha=1.0_dp)
         CALL dbcsr_set(f_im, 0.0_dp)

         CALL qs_ot_complex_multiply('N', 'N', matrix_s, matrix_s_im, c_re, c_im, sc_re, sc_im, tmp_nk)
      END IF

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL qs_ot_generate_rotation_complex(qs_ot_env)

         CALL dbcsr_copy(rotated_re, matrix_c, name="rotated_re")
         CALL dbcsr_copy(rotated_im, matrix_c_im, name="rotated_im")
         CALL dbcsr_copy(rotation_tmp, matrix_c, name="rotation_tmp")

         ! C_out = Q*U for complex Q and unitary U.
         CALL qs_ot_complex_multiply('N', 'N', matrix_c, matrix_c_im, qs_ot_env%rot_mat_u, &
                                     qs_ot_env%rot_mat_u_im, rotated_re, rotated_im, rotation_tmp)

         CALL dbcsr_copy(matrix_c, rotated_re)
         CALL dbcsr_copy(matrix_c_im, rotated_im)
         CALL dbcsr_release(rotated_re)
         CALL dbcsr_release(rotated_im)
         CALL dbcsr_release(rotation_tmp)
      END IF

      DEALLOCATE (eigenvalues, inverse_sqrt)

      CALL timestop(handle)
   END SUBROUTINE qs_ot_get_orbitals_ref_complex

! **************************************************************************************************
!> \brief complex k-point REF derivative dE/dX from H(k)C(k), S(k)C(k), and C(k)
!> \param matrix_hc ...
!> \param matrix_hc_im ...
!> \param qs_ot_env ...
!> \param matrix_hc_rotation occupation-weighted H(k)C(k) for the rotation channel
!> \param matrix_hc_rotation_im imaginary component of matrix_hc_rotation
! **************************************************************************************************
   SUBROUTINE qs_ot_get_derivative_ref_complex(matrix_hc, matrix_hc_im, qs_ot_env, &
                                               matrix_hc_rotation, matrix_hc_rotation_im)
      TYPE(dbcsr_type), POINTER                          :: matrix_hc, matrix_hc_im
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: matrix_hc_rotation, matrix_hc_rotation_im

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative_ref_complex'

      INTEGER                                            :: handle, k, n
      TYPE(dbcsr_type)                                   :: tmp_nk
      TYPE(dbcsr_type), POINTER :: b_im, b_re, f_im, f_re, g_im, g_re, hc_im, hc_re, &
         hc_rotation_im, hc_rotation_re, hc_work_im, hc_work_re, q_im, q_re, sc_im, sc_re, tmp_kk
      TYPE(dbcsr_type), TARGET                           :: hc_rot_im, hc_rot_re

      CALL timeset(routineN, handle)

      CPASSERT(qs_ot_env%has_complex_kpoint_state)
      CPASSERT(ASSOCIATED(matrix_hc))
      CPASSERT(ASSOCIATED(matrix_hc_im))

      f_re => qs_ot_env%matrix_ref_inv_sqrt
      f_im => qs_ot_env%matrix_ref_inv_sqrt_im
      sc_re => qs_ot_env%matrix_sx
      sc_im => qs_ot_env%matrix_sx_im
      hc_re => matrix_hc
      hc_im => matrix_hc_im
      hc_rotation_re => hc_re
      hc_rotation_im => hc_im
      IF (PRESENT(matrix_hc_rotation) .OR. PRESENT(matrix_hc_rotation_im)) THEN
         CPASSERT(PRESENT(matrix_hc_rotation) .AND. PRESENT(matrix_hc_rotation_im))
         CPASSERT(ASSOCIATED(matrix_hc_rotation))
         CPASSERT(ASSOCIATED(matrix_hc_rotation_im))
         hc_rotation_re => matrix_hc_rotation
         hc_rotation_im => matrix_hc_rotation_im
      END IF
      hc_work_re => hc_re
      hc_work_im => hc_im
      g_re => qs_ot_env%matrix_gx
      g_im => qs_ot_env%matrix_gx_im
      b_re => qs_ot_env%buf1_k_k_sym
      b_im => qs_ot_env%buf2_k_k_sym
      tmp_kk => qs_ot_env%buf3_k_k_sym
      q_re => qs_ot_env%buf1_n_k
      q_im => qs_ot_env%buf1_n_k_dp

      CALL dbcsr_get_info(sc_re, nfullrows_total=n, nfullcols_total=k)

      ! Q = X*(X^H*S*X)^(-1/2), reconstructed from the current REF coordinate.
      CALL qs_ot_complex_multiply('N', 'N', qs_ot_env%matrix_x, qs_ot_env%matrix_x_im, f_re, f_im, q_re, q_im, &
                                  g_re)

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL qs_ot_generate_rotation_complex(qs_ot_env)

         ! dF/dU = Q^H*G_C for C=Q*U.
         CALL dbcsr_multiply('T', 'N', 1.0_dp, q_re, hc_rotation_re, &
                             0.0_dp, qs_ot_env%rot_mat_dedu)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, q_im, hc_rotation_im, 0.0_dp, tmp_kk)
         CALL dbcsr_add(qs_ot_env%rot_mat_dedu, tmp_kk, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

         CALL dbcsr_multiply('T', 'N', 1.0_dp, q_re, hc_rotation_im, &
                             0.0_dp, qs_ot_env%rot_mat_dedu_im)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, q_im, hc_rotation_re, 0.0_dp, tmp_kk)
         CALL dbcsr_add(qs_ot_env%rot_mat_dedu_im, tmp_kk, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
         CALL qs_ot_rot_mat_derivative_complex(qs_ot_env)

         ! The REF coordinate sees G_Q=G_C*U^H.
         CALL dbcsr_copy(hc_rot_re, hc_re, name="hc_rot_re")
         CALL dbcsr_copy(hc_rot_im, hc_im, name="hc_rot_im")
         CALL dbcsr_copy(tmp_nk, hc_re, name="tmp_nk")
         CALL dbcsr_multiply('N', 'T', 1.0_dp, hc_re, qs_ot_env%rot_mat_u, &
                             0.0_dp, hc_rot_re)
         CALL dbcsr_multiply('N', 'T', 1.0_dp, hc_im, qs_ot_env%rot_mat_u_im, &
                             0.0_dp, tmp_nk)
         CALL dbcsr_add(hc_rot_re, tmp_nk, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

         CALL dbcsr_multiply('N', 'T', 1.0_dp, hc_im, qs_ot_env%rot_mat_u, &
                             0.0_dp, hc_rot_im)
         CALL dbcsr_multiply('N', 'T', 1.0_dp, hc_re, qs_ot_env%rot_mat_u_im, &
                             0.0_dp, tmp_nk)
         CALL dbcsr_add(hc_rot_im, tmp_nk, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
         hc_work_re => hc_rot_re
         hc_work_im => hc_rot_im
      END IF

      ! B = Q^H*G_Q.  For uniform fixed occupations B is Hermitian.
      CALL dbcsr_multiply('T', 'N', 1.0_dp, q_re, hc_work_re, 0.0_dp, b_re)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, q_im, hc_work_im, 0.0_dp, tmp_kk)
      CALL dbcsr_add(b_re, tmp_kk, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      CALL dbcsr_multiply('T', 'N', 1.0_dp, q_re, hc_work_im, 0.0_dp, b_im)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, q_im, hc_work_re, 0.0_dp, tmp_kk)
      CALL dbcsr_add(b_im, tmp_kk, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      IF (qs_ot_env%settings%eps_irac_filter_matrix > 0.0_dp) THEN
         CALL dbcsr_filter(b_re, qs_ot_env%settings%eps_irac_filter_matrix)
         CALL dbcsr_filter(b_im, qs_ot_env%settings%eps_irac_filter_matrix)
      END IF

      ! S*Q = (S*X)*F.  G is temporary storage for this pair.
      CALL qs_ot_complex_multiply('N', 'N', sc_re, sc_im, f_re, f_im, g_re, g_im, q_re)

      ! Form (S*Q)*B.  The Q workspaces are free after B has been formed.
      CALL dbcsr_multiply('N', 'N', 1.0_dp, g_re, b_re, 0.0_dp, q_re)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, g_im, b_im, 0.0_dp, q_im)
      CALL dbcsr_add(q_re, q_im, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      CALL dbcsr_multiply('N', 'N', 1.0_dp, g_re, b_im, 0.0_dp, q_im)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, g_im, b_re, 0.0_dp, g_re)
      CALL dbcsr_add(q_im, g_re, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      CALL dbcsr_add(q_re, hc_work_re, alpha_scalar=-1.0_dp, beta_scalar=1.0_dp)
      CALL dbcsr_add(q_im, hc_work_im, alpha_scalar=-1.0_dp, beta_scalar=1.0_dp)

      ! Pull the physical gradient back to the finite REF coordinate: G_X = (G_C-S*Q*B)*F.
      CALL dbcsr_multiply('N', 'N', 1.0_dp, q_re, f_re, 0.0_dp, g_re)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, q_im, f_im, 0.0_dp, g_im)
      CALL dbcsr_add(g_re, g_im, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      CALL dbcsr_multiply('N', 'N', 1.0_dp, q_re, f_im, 0.0_dp, g_im)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, q_im, f_re, 0.0_dp, q_re)
      CALL dbcsr_add(g_im, q_re, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL qs_ot_add_ref_vertical_response_complex(b_re, b_im, f_re, f_im, &
                                                      sc_re, sc_im, g_re, g_im, qs_ot_env)
      END IF

      IF (qs_ot_env%settings%eps_irac_filter_matrix > 0.0_dp) THEN
         CALL dbcsr_filter(g_re, qs_ot_env%settings%eps_irac_filter_matrix)
         CALL dbcsr_filter(g_im, qs_ot_env%settings%eps_irac_filter_matrix)
      END IF

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL dbcsr_release(hc_rot_re)
         CALL dbcsr_release(hc_rot_im)
         CALL dbcsr_release(tmp_nk)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_ot_get_derivative_ref_complex

! **************************************************************************************************
!> \brief Add the vertical response of the finite complex polar REF map.
!> \param a_re real part of Q^H G_Q
!> \param a_im imaginary part of Q^H G_Q
!> \param f_re real part of (X^H S X)^(-1/2)
!> \param f_im imaginary part of (X^H S X)^(-1/2)
!> \param sx_re real part of S X
!> \param sx_im imaginary part of S X
!> \param gradient_re real REF gradient, updated in place
!> \param gradient_im imaginary REF gradient, updated in place
!> \param qs_ot_env complex OT environment
! **************************************************************************************************
   SUBROUTINE qs_ot_add_ref_vertical_response_complex(a_re, a_im, f_re, f_im, &
                                                      sx_re, sx_im, gradient_re, gradient_im, &
                                                      qs_ot_env)
      TYPE(dbcsr_type)                                   :: a_re, a_im, f_re, f_im, sx_re, sx_im, &
                                                            gradient_re, gradient_im
      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_add_ref_vertical_response_complex'

      INTEGER                                            :: col, handle, i, j, k, row
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_offset, row_blk_offset
      LOGICAL                                            :: found_im, found_re
      REAL(KIND=dp)                                      :: denominator
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: f_evals
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_im, block_re
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type) :: anti_im, anti_re, f_work_im, f_work_re, inner_im, inner_re, sq_im, &
         sq_re, tmp_kk, tmp_nk, v_im, v_re, vertical_im, vertical_re, work_im, work_re, z_im, z_re

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(a_re, nfullrows_total=k)
      IF (k == 0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! Only the anti-Hermitian part of A=Q^H G_Q couples to the unitary
      ! component of the polar differential.
      CALL qs_ot_square_transpose(a_re, tmp_kk, qs_ot_env%rot_mat_u)
      CALL dbcsr_copy(anti_re, a_re, name='anti_re')
      CALL dbcsr_add(anti_re, tmp_kk, alpha_scalar=0.5_dp, beta_scalar=-0.5_dp)
      CALL qs_ot_square_transpose(a_im, tmp_kk, qs_ot_env%rot_mat_u)
      CALL dbcsr_copy(anti_im, a_im, name='anti_im')
      CALL dbcsr_add(anti_im, tmp_kk, alpha_scalar=0.5_dp, beta_scalar=0.5_dp)

      ! P=(X^H S X)^(1/2)=F^(-1).  In the eigenbasis of F, solve
      ! P Z + Z P = antiherm(A) element by element.
      CALL dbcsr_copy(f_work_re, f_re, name='f_work_re')
      CALL dbcsr_copy(f_work_im, f_im, name='f_work_im')
      CALL dbcsr_copy(v_re, f_re, name='v_re')
      CALL dbcsr_copy(v_im, f_im, name='v_im')
      ALLOCATE (f_evals(k))
      CALL cp_dbcsr_heevd(matrix_re=f_work_re, matrix_im=f_work_im, &
                          eigenvectors_re=v_re, eigenvectors_im=v_im, eigenvalues=f_evals, &
                          para_env=qs_ot_env%para_env, blacs_env=qs_ot_env%blacs_env)
      IF (MINVAL(f_evals) <= EPSILON(1.0_dp)) THEN
         CPABORT('Complex REF inverse square root is not positive definite')
      END IF

      ! inner = V^H antiherm(A) V.
      CALL dbcsr_copy(work_re, anti_re, name='work_re')
      CALL dbcsr_copy(work_im, anti_im, name='work_im')
      CALL qs_ot_complex_multiply('N', 'N', anti_re, anti_im, v_re, v_im, work_re, work_im, tmp_kk)

      CALL dbcsr_copy(inner_re, anti_re, name='inner_re')
      CALL dbcsr_copy(inner_im, anti_im, name='inner_im')
      CALL dbcsr_multiply('T', 'N', 1.0_dp, v_re, work_re, 0.0_dp, inner_re)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, v_im, work_im, 1.0_dp, inner_re)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, v_re, work_im, 0.0_dp, inner_im)
      CALL dbcsr_multiply('T', 'N', -1.0_dp, v_im, work_re, 1.0_dp, inner_im)

      CALL dbcsr_get_info(inner_re, row_blk_offset=row_blk_offset, &
                          col_blk_offset=col_blk_offset)
      CALL dbcsr_iterator_start(iter, inner_re)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col)
         CALL dbcsr_get_block_p(inner_re, row, col, block_re, found_re)
         CALL dbcsr_get_block_p(inner_im, row, col, block_im, found_im)
         CPASSERT(found_re .AND. found_im)
         DO i = 1, SIZE(block_re, 1)
         DO j = 1, SIZE(block_re, 2)
            denominator = 1.0_dp/f_evals(row_blk_offset(row) + i - 1) + &
                          1.0_dp/f_evals(col_blk_offset(col) + j - 1)
            block_re(i, j) = block_re(i, j)/denominator
            block_im(i, j) = block_im(i, j)/denominator
         END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      ! Z = V inner V^H.
      CALL qs_ot_complex_multiply('N', 'N', v_re, v_im, inner_re, inner_im, work_re, work_im, tmp_kk)

      CALL dbcsr_copy(z_re, anti_re, name='z_re')
      CALL dbcsr_copy(z_im, anti_im, name='z_im')
      CALL dbcsr_multiply('N', 'T', 1.0_dp, work_re, v_re, 0.0_dp, z_re)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, work_im, v_im, 1.0_dp, z_re)
      CALL dbcsr_multiply('N', 'T', 1.0_dp, work_im, v_re, 0.0_dp, z_im)
      CALL dbcsr_multiply('N', 'T', -1.0_dp, work_re, v_im, 1.0_dp, z_im)

      ! The adjoint vertical contribution is 2 S Q Z, with S Q=(S X)F.
      CALL dbcsr_copy(sq_re, sx_re, name='sq_re')
      CALL dbcsr_copy(sq_im, sx_im, name='sq_im')
      CALL dbcsr_copy(tmp_nk, sx_re, name='tmp_nk')
      CALL qs_ot_complex_multiply('N', 'N', sx_re, sx_im, f_re, f_im, sq_re, sq_im, tmp_nk)

      CALL dbcsr_copy(vertical_re, sx_re, name='vertical_re')
      CALL dbcsr_copy(vertical_im, sx_im, name='vertical_im')
      CALL qs_ot_complex_multiply('N', 'N', sq_re, sq_im, z_re, z_im, vertical_re, vertical_im, tmp_nk)
      CALL dbcsr_add(gradient_re, vertical_re, alpha_scalar=1.0_dp, beta_scalar=2.0_dp)
      CALL dbcsr_add(gradient_im, vertical_im, alpha_scalar=1.0_dp, beta_scalar=2.0_dp)

      DEALLOCATE (f_evals)
      CALL dbcsr_release(anti_im)
      CALL dbcsr_release(anti_re)
      CALL dbcsr_release(f_work_im)
      CALL dbcsr_release(f_work_re)
      CALL dbcsr_release(inner_im)
      CALL dbcsr_release(inner_re)
      CALL dbcsr_release(sq_im)
      CALL dbcsr_release(sq_re)
      CALL dbcsr_release(tmp_kk)
      CALL dbcsr_release(tmp_nk)
      CALL dbcsr_release(v_im)
      CALL dbcsr_release(v_re)
      CALL dbcsr_release(vertical_im)
      CALL dbcsr_release(vertical_re)
      CALL dbcsr_release(work_im)
      CALL dbcsr_release(work_re)
      CALL dbcsr_release(z_im)
      CALL dbcsr_release(z_re)

      CALL timestop(handle)

   END SUBROUTINE qs_ot_add_ref_vertical_response_complex

! **************************************************************************************************
!> \brief computes U=exp(A) for the complex anti-Hermitian generator
!>        A=rot_mat_x+i*rot_mat_x_im
!> \param qs_ot_env a complex k-point OT environment
! **************************************************************************************************
   SUBROUTINE qs_ot_generate_rotation_complex(qs_ot_env)

      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_generate_rotation_complex'

      INTEGER                                            :: handle, k
      REAL(KIND=dp)                                      :: rot_norm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: exp_evals_im, exp_evals_re
      TYPE(dbcsr_type)                                   :: h_re, tmp, w_im, w_re

      CALL timeset(routineN, handle)

      CPASSERT(qs_ot_env%has_complex_kpoint_state)
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_x_im))
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_u_im))

      CALL dbcsr_get_info(qs_ot_env%rot_mat_x, nfullrows_total=k)
      IF (k /= 0) THEN
         rot_norm = dbcsr_frobenius_norm(qs_ot_env%rot_mat_x) + &
                    dbcsr_frobenius_norm(qs_ot_env%rot_mat_x_im)
         IF (rot_norm <= EPSILON(1.0_dp)) THEN
            CALL dbcsr_set(qs_ot_env%rot_mat_u, 0.0_dp)
            CALL dbcsr_add_on_diag(qs_ot_env%rot_mat_u, 1.0_dp)
            CALL dbcsr_set(qs_ot_env%rot_mat_u_im, 0.0_dp)
            CALL timestop(handle)
            RETURN
         END IF

         ! i*A = i*X-Y is Hermitian. Its eigenvectors give
         ! exp(A)=V*diag(exp(-i*lambda))*V^H.
         CALL dbcsr_copy(h_re, qs_ot_env%rot_mat_x_im, name="h_re")
         CALL dbcsr_scale(h_re, -1.0_dp)
         CALL cp_dbcsr_heevd(matrix_re=h_re, matrix_im=qs_ot_env%rot_mat_x, &
                             eigenvectors_re=qs_ot_env%rot_mat_evec_re, &
                             eigenvectors_im=qs_ot_env%rot_mat_evec_im, &
                             eigenvalues=qs_ot_env%rot_mat_evals, &
                             para_env=qs_ot_env%para_env, blacs_env=qs_ot_env%blacs_env)

         ALLOCATE (exp_evals_re(k), exp_evals_im(k))
         exp_evals_re(:) = COS(-qs_ot_env%rot_mat_evals(:))
         exp_evals_im(:) = SIN(-qs_ot_env%rot_mat_evals(:))

         CALL dbcsr_copy(w_re, qs_ot_env%rot_mat_evec_re, name="w_re")
         CALL dbcsr_scale_by_vector(w_re, alpha=exp_evals_re, side='right')
         CALL dbcsr_copy(tmp, qs_ot_env%rot_mat_evec_im, name="tmp")
         CALL dbcsr_scale_by_vector(tmp, alpha=exp_evals_im, side='right')
         CALL dbcsr_add(w_re, tmp, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

         CALL dbcsr_copy(w_im, qs_ot_env%rot_mat_evec_im, name="w_im")
         CALL dbcsr_scale_by_vector(w_im, alpha=exp_evals_re, side='right')
         CALL dbcsr_copy(tmp, qs_ot_env%rot_mat_evec_re)
         CALL dbcsr_scale_by_vector(tmp, alpha=exp_evals_im, side='right')
         CALL dbcsr_add(w_im, tmp, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

         CALL dbcsr_multiply('N', 'T', 1.0_dp, w_re, qs_ot_env%rot_mat_evec_re, &
                             0.0_dp, qs_ot_env%rot_mat_u)
         CALL dbcsr_multiply('N', 'T', 1.0_dp, w_im, qs_ot_env%rot_mat_evec_im, &
                             1.0_dp, qs_ot_env%rot_mat_u)
         CALL dbcsr_multiply('N', 'T', 1.0_dp, w_im, qs_ot_env%rot_mat_evec_re, &
                             0.0_dp, qs_ot_env%rot_mat_u_im)
         CALL dbcsr_multiply('N', 'T', -1.0_dp, w_re, qs_ot_env%rot_mat_evec_im, &
                             1.0_dp, qs_ot_env%rot_mat_u_im)

         CALL dbcsr_release(h_re)
         CALL dbcsr_release(tmp)
         CALL dbcsr_release(w_re)
         CALL dbcsr_release(w_im)
         DEALLOCATE (exp_evals_re, exp_evals_im)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_ot_generate_rotation_complex

! **************************************************************************************************
!> \brief pull the complex dE/dU covector back to the anti-Hermitian generator
!>        using the adjoint Frechet derivative of exp
!> \param qs_ot_env a complex k-point OT environment with an up-to-date U
! **************************************************************************************************
   SUBROUTINE qs_ot_rot_mat_derivative_complex(qs_ot_env)
      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_rot_mat_derivative_complex'

      INTEGER                                            :: handle, k
      REAL(KIND=dp)                                      :: rot_norm
      TYPE(dbcsr_type)                                   :: frechet_im, frechet_re, inner_deriv_im, &
                                                            inner_deriv_re, outer_deriv_im, &
                                                            outer_deriv_re, tmp, work_im, work_re

      CALL timeset(routineN, handle)

      CPASSERT(qs_ot_env%has_complex_kpoint_state)
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_dedu_im))
      CPASSERT(ASSOCIATED(qs_ot_env%rot_mat_gx_im))

      CALL dbcsr_get_info(qs_ot_env%rot_mat_u, nfullrows_total=k)
      IF (k /= 0) THEN
         rot_norm = dbcsr_frobenius_norm(qs_ot_env%rot_mat_x) + &
                    dbcsr_frobenius_norm(qs_ot_env%rot_mat_x_im)
         IF (rot_norm <= EPSILON(1.0_dp)) THEN
            CALL qs_ot_square_transpose(qs_ot_env%rot_mat_dedu, tmp, qs_ot_env%rot_mat_u)
            CALL dbcsr_copy(qs_ot_env%rot_mat_gx, qs_ot_env%rot_mat_dedu)
            CALL dbcsr_add(qs_ot_env%rot_mat_gx, tmp, &
                           alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
            CALL dbcsr_release(tmp)

            CALL qs_ot_square_transpose(qs_ot_env%rot_mat_dedu_im, tmp, qs_ot_env%rot_mat_u)
            CALL dbcsr_copy(qs_ot_env%rot_mat_gx_im, qs_ot_env%rot_mat_dedu_im)
            CALL dbcsr_add(qs_ot_env%rot_mat_gx_im, tmp, &
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            CALL dbcsr_release(tmp)
            CALL timestop(handle)
            RETURN
         END IF

         CALL dbcsr_copy(work_re, qs_ot_env%rot_mat_dedu, name="work_re")
         CALL dbcsr_copy(work_im, qs_ot_env%rot_mat_dedu, name="work_im")
         CALL dbcsr_copy(tmp, qs_ot_env%rot_mat_dedu, name="tmp")
         CALL dbcsr_copy(inner_deriv_re, qs_ot_env%rot_mat_dedu, name="inner_deriv_re")
         CALL dbcsr_copy(inner_deriv_im, qs_ot_env%rot_mat_dedu, name="inner_deriv_im")

         ! V^H*(dE/dU)*V, split into real and imaginary parts.
         CALL qs_ot_complex_multiply('N', 'N', qs_ot_env%rot_mat_dedu, qs_ot_env%rot_mat_dedu_im, &
                                     qs_ot_env%rot_mat_evec_re, qs_ot_env%rot_mat_evec_im, work_re, work_im, &
                                     tmp)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_evec_re, &
                             work_re, 0.0_dp, inner_deriv_re)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_evec_im, &
                             work_im, 1.0_dp, inner_deriv_re)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, qs_ot_env%rot_mat_evec_re, &
                             work_im, 0.0_dp, inner_deriv_im)
         CALL dbcsr_multiply('T', 'N', -1.0_dp, qs_ot_env%rot_mat_evec_im, &
                             work_re, 1.0_dp, inner_deriv_im)

         CALL qs_ot_apply_complex_frechet_dbcsr(qs_ot_env%rot_mat_evals, &
                                                inner_deriv_re, inner_deriv_im, &
                                                outer_deriv_re, outer_deriv_im, adjoint=.TRUE.)

         CALL qs_ot_complex_multiply('N', 'N', qs_ot_env%rot_mat_evec_re, qs_ot_env%rot_mat_evec_im, &
                                     outer_deriv_re, outer_deriv_im, work_re, work_im, tmp)

         CALL dbcsr_copy(frechet_re, qs_ot_env%rot_mat_dedu, name="frechet_re")
         CALL dbcsr_copy(frechet_im, qs_ot_env%rot_mat_dedu, name="frechet_im")
         CALL dbcsr_multiply('N', 'T', 1.0_dp, work_re, qs_ot_env%rot_mat_evec_re, &
                             0.0_dp, frechet_re)
         CALL dbcsr_multiply('N', 'T', 1.0_dp, work_im, qs_ot_env%rot_mat_evec_im, &
                             1.0_dp, frechet_re)
         CALL dbcsr_multiply('N', 'T', 1.0_dp, work_im, qs_ot_env%rot_mat_evec_re, &
                             0.0_dp, frechet_im)
         CALL dbcsr_multiply('N', 'T', -1.0_dp, work_re, qs_ot_env%rot_mat_evec_im, &
                             1.0_dp, frechet_im)

         ! Tangents satisfy X^T=-X and Y^T=Y for A=X+iY.
         CALL qs_ot_square_transpose(frechet_re, tmp, qs_ot_env%rot_mat_u)
         CALL dbcsr_copy(qs_ot_env%rot_mat_gx, frechet_re)
         CALL dbcsr_add(qs_ot_env%rot_mat_gx, tmp, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
         CALL qs_ot_square_transpose(frechet_im, tmp, qs_ot_env%rot_mat_u)
         CALL dbcsr_copy(qs_ot_env%rot_mat_gx_im, frechet_im)
         CALL dbcsr_add(qs_ot_env%rot_mat_gx_im, tmp, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

         CALL dbcsr_release(frechet_re)
         CALL dbcsr_release(frechet_im)
         CALL dbcsr_release(inner_deriv_re)
         CALL dbcsr_release(inner_deriv_im)
         CALL dbcsr_release(outer_deriv_re)
         CALL dbcsr_release(outer_deriv_im)
         CALL dbcsr_release(tmp)
         CALL dbcsr_release(work_re)
         CALL dbcsr_release(work_im)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_ot_rot_mat_derivative_complex

! **************************************************************************************************
!> \brief compute P=X^H*S*X and the STRICT matrix functions for a complex K-point channel
!> \param matrix_x real part of X
!> \param matrix_x_im imaginary part of X
!> \param matrix_sx real part of S*X
!> \param matrix_sx_im imaginary part of S*X
!> \param qs_ot_env OT channel state
!> \param zero_coordinates initialize the known X=0 state without diagonalization
! **************************************************************************************************
   SUBROUTINE qs_ot_get_p_complex(matrix_x, matrix_x_im, matrix_sx, matrix_sx_im, qs_ot_env, zero_coordinates)
      TYPE(dbcsr_type), POINTER                          :: matrix_x, matrix_x_im, matrix_sx, &
                                                            matrix_sx_im
      TYPE(qs_ot_type)                                   :: qs_ot_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: zero_coordinates

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_p_complex'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(qs_ot_env%has_complex_kpoint_state)
      CPASSERT(qs_ot_env%settings%ot_algorithm == 'TOD')

      IF (PRESENT(zero_coordinates)) THEN
         IF (zero_coordinates) THEN
            CALL dbcsr_set(qs_ot_env%matrix_p, 0.0_dp)
            CALL dbcsr_set(qs_ot_env%matrix_p_im, 0.0_dp)
            CALL dbcsr_set(qs_ot_env%matrix_r, 0.0_dp)
            CALL dbcsr_add_on_diag(qs_ot_env%matrix_r, 1.0_dp)
            CALL dbcsr_set(qs_ot_env%matrix_r_im, 0.0_dp)
            CALL dbcsr_copy(qs_ot_env%matrix_cosp, qs_ot_env%matrix_r)
            CALL dbcsr_copy(qs_ot_env%matrix_sinp, qs_ot_env%matrix_r)
            CALL dbcsr_set(qs_ot_env%matrix_cosp_im, 0.0_dp)
            CALL dbcsr_set(qs_ot_env%matrix_sinp_im, 0.0_dp)
            CALL dbcsr_reserve_all_blocks(qs_ot_env%matrix_cosp_b)
            CALL dbcsr_reserve_all_blocks(qs_ot_env%matrix_sinp_b)
            CALL dbcsr_set(qs_ot_env%matrix_cosp_b, -0.5_dp)
            CALL dbcsr_set(qs_ot_env%matrix_sinp_b, -1.0_dp/6.0_dp)
            qs_ot_env%evals(:) = 0.0_dp
            qs_ot_env%dum(:) = 0.0_dp
            qs_ot_env%do_taylor = .FALSE.
            CALL timestop(handle)
            RETURN
         END IF
      END IF

      CALL qs_ot_complex_multiply('C', 'N', matrix_x, matrix_x_im, matrix_sx, matrix_sx_im, &
                                  qs_ot_env%matrix_p, qs_ot_env%matrix_p_im, qs_ot_env%matrix_buf1)
      qs_ot_env%do_taylor = .FALSE.
      CALL qs_ot_p2m_diag(qs_ot_env)

      CALL timestop(handle)

   END SUBROUTINE qs_ot_get_p_complex

! **************************************************************************************************
!> \brief update complex K-point orbitals with the finite STRICT transformation
!> \param matrix_c real output orbitals
!> \param matrix_c_im imaginary output orbitals
!> \param matrix_s real overlap matrix
!> \param matrix_s_im imaginary overlap matrix
!> \param qs_ot_env OT channel state
! **************************************************************************************************
   SUBROUTINE qs_ot_get_orbitals_complex(matrix_c, matrix_c_im, matrix_s, matrix_s_im, qs_ot_env)
      TYPE(dbcsr_type), POINTER                          :: matrix_c, matrix_c_im, matrix_s, &
                                                            matrix_s_im
      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_orbitals_complex'

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: rotated_im, rotated_re, rotation_tmp

      CALL timeset(routineN, handle)

      CPASSERT(qs_ot_env%has_complex_kpoint_state)
      CPASSERT(qs_ot_env%settings%ot_algorithm == 'TOD')

      CALL qs_ot_complex_multiply('N', 'N', matrix_s, matrix_s_im, &
                                  qs_ot_env%matrix_x, qs_ot_env%matrix_x_im, &
                                  qs_ot_env%matrix_sx, qs_ot_env%matrix_sx_im, &
                                  qs_ot_env%matrix_tmp_nk)
      CALL qs_ot_get_p_complex(qs_ot_env%matrix_x, qs_ot_env%matrix_x_im, &
                               qs_ot_env%matrix_sx, qs_ot_env%matrix_sx_im, qs_ot_env)

      CALL qs_ot_complex_multiply('N', 'N', qs_ot_env%matrix_c0, qs_ot_env%matrix_c0_im, &
                                  qs_ot_env%matrix_cosp, qs_ot_env%matrix_cosp_im, &
                                  matrix_c, matrix_c_im, qs_ot_env%matrix_tmp_nk)
      CALL qs_ot_complex_multiply('N', 'N', qs_ot_env%matrix_x, qs_ot_env%matrix_x_im, &
                                  qs_ot_env%matrix_sinp, qs_ot_env%matrix_sinp_im, &
                                  qs_ot_env%matrix_buf_nk, qs_ot_env%matrix_buf_nk_im, &
                                  qs_ot_env%matrix_tmp_nk)
      CALL dbcsr_add(matrix_c, qs_ot_env%matrix_buf_nk, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      CALL dbcsr_add(matrix_c_im, qs_ot_env%matrix_buf_nk_im, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL qs_ot_generate_rotation_complex(qs_ot_env)

         CALL dbcsr_copy(rotated_re, matrix_c, name="strict_rotated_re")
         CALL dbcsr_copy(rotated_im, matrix_c_im, name="strict_rotated_im")
         CALL dbcsr_copy(rotation_tmp, matrix_c, name="strict_rotation_tmp")

         ! C_out = Q(X)*U for the finite STRICT chart Q(X).
         CALL qs_ot_complex_multiply('N', 'N', matrix_c, matrix_c_im, qs_ot_env%rot_mat_u, &
                                     qs_ot_env%rot_mat_u_im, rotated_re, rotated_im, rotation_tmp)

         CALL dbcsr_copy(matrix_c, rotated_re)
         CALL dbcsr_copy(matrix_c_im, rotated_im)
         CALL dbcsr_release(rotated_re)
         CALL dbcsr_release(rotated_im)
         CALL dbcsr_release(rotation_tmp)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_ot_get_orbitals_complex

! **************************************************************************************************
!> \brief Project a complex STRICT gradient; real minimizer fallbacks use the shared real projection.
!> \param qs_ot_env OT channel state
!> \param matrix_gx real component to project in place
!> \param matrix_gx_im imaginary component for a complex channel
!> \param use_preconditioner use the attached preconditioner, default true
! **************************************************************************************************
   SUBROUTINE qs_ot_project_gradient(qs_ot_env, matrix_gx, matrix_gx_im, use_preconditioner)
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type), POINTER                          :: matrix_gx
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: matrix_gx_im
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_preconditioner

      LOGICAL                                            :: preconditioned
      TYPE(dbcsr_type), POINTER                          :: target_im, target_re

      IF (.NOT. qs_ot_env%has_complex_kpoint_state) THEN
         CALL qs_ot_project_gradient_real(qs_ot_env, matrix_gx, use_preconditioner)
         RETURN
      END IF
      preconditioned = .TRUE.
      IF (PRESENT(use_preconditioner)) preconditioned = use_preconditioner
      IF (.NOT. preconditioned) qs_ot_env%os_valid = .FALSE.
      CPASSERT(PRESENT(matrix_gx_im))
      CALL qs_ot_prepare_complex_tangent_metric(qs_ot_env, use_preconditioner=preconditioned)
      IF (preconditioned .AND. ASSOCIATED(qs_ot_env%preconditioner)) THEN
         target_re => qs_ot_env%matrix_psc0
         target_im => qs_ot_env%matrix_psc0_im
      ELSE
         target_re => qs_ot_env%matrix_sc0
         target_im => qs_ot_env%matrix_sc0_im
      END IF

      CALL qs_ot_complex_multiply('C', 'N', target_re, target_im, &
                                  matrix_gx, matrix_gx_im, &
                                  qs_ot_env%matrix_buf1_ortho, qs_ot_env%matrix_buf1_ortho_im, &
                                  qs_ot_env%matrix_tmp_ortho)
      CALL qs_ot_complex_multiply('N', 'N', qs_ot_env%matrix_os, qs_ot_env%matrix_os_im, &
                                  qs_ot_env%matrix_buf1_ortho, qs_ot_env%matrix_buf1_ortho_im, &
                                  qs_ot_env%matrix_buf2_ortho, qs_ot_env%matrix_buf2_ortho_im, &
                                  qs_ot_env%matrix_tmp_ortho)
      CALL qs_ot_complex_multiply('N', 'N', qs_ot_env%matrix_sc0, qs_ot_env%matrix_sc0_im, &
                                  qs_ot_env%matrix_buf2_ortho, qs_ot_env%matrix_buf2_ortho_im, &
                                  qs_ot_env%matrix_buf_nk, qs_ot_env%matrix_buf_nk_im, &
                                  qs_ot_env%matrix_tmp_nk)
      CALL dbcsr_add(matrix_gx, qs_ot_env%matrix_buf_nk, &
                     alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      CALL dbcsr_add(matrix_gx_im, qs_ot_env%matrix_buf_nk_im, &
                     alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

   END SUBROUTINE qs_ot_project_gradient

! **************************************************************************************************
!> \brief Prepare the inverse metric used to project a complex STRICT gradient.
!>        An unusable preconditioner is detached before any minimizer history is updated.
!> \param qs_ot_env OT channel state
!> \param preconditioner_rejected true if the attached preconditioner was not positive definite
!> \param use_preconditioner use the attached preconditioner, default true
! **************************************************************************************************
   SUBROUTINE qs_ot_prepare_complex_tangent_metric(qs_ot_env, preconditioner_rejected, use_preconditioner)
      TYPE(qs_ot_type)                                   :: qs_ot_env
      LOGICAL, INTENT(OUT), OPTIONAL                     :: preconditioner_rejected
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_preconditioner

      INTEGER                                            :: i, k
      LOGICAL                                            :: preconditioned
      REAL(KIND=dp)                                      :: eval_scale, eval_threshold
      TYPE(dbcsr_type), POINTER                          :: target_im, target_re

      preconditioned = .TRUE.
      IF (PRESENT(use_preconditioner)) preconditioned = use_preconditioner
      IF (.NOT. preconditioned) qs_ot_env%os_valid = .FALSE.
      IF (PRESENT(preconditioner_rejected)) preconditioner_rejected = .FALSE.
      IF (qs_ot_env%os_valid) RETURN

      IF (preconditioned .AND. ASSOCIATED(qs_ot_env%preconditioner)) THEN
         target_re => qs_ot_env%matrix_psc0
         target_im => qs_ot_env%matrix_psc0_im
         CALL apply_preconditioner(qs_ot_env%preconditioner, &
                                   qs_ot_env%matrix_sc0, qs_ot_env%matrix_sc0_im, &
                                   target_re, target_im)
      ELSE
         target_re => qs_ot_env%matrix_sc0
         target_im => qs_ot_env%matrix_sc0_im
      END IF
      CALL qs_ot_complex_multiply('C', 'N', qs_ot_env%matrix_sc0, qs_ot_env%matrix_sc0_im, &
                                  target_re, target_im, qs_ot_env%matrix_os, qs_ot_env%matrix_os_im, &
                                  qs_ot_env%matrix_buf1)
      CALL dbcsr_get_info(qs_ot_env%matrix_os, nfullrows_total=k)
      ! Preserve the eigenvectors and eigenvalues of X^H*S*X used by the STRICT derivative.
      CALL cp_dbcsr_heevd(matrix_re=qs_ot_env%matrix_os, matrix_im=qs_ot_env%matrix_os_im, &
                          eigenvectors_re=qs_ot_env%matrix_buf3, eigenvectors_im=qs_ot_env%matrix_buf3_im, &
                          eigenvalues=qs_ot_env%dum, para_env=qs_ot_env%para_env, &
                          blacs_env=qs_ot_env%blacs_env)
      eval_scale = MAX(1.0_dp, MAXVAL(ABS(qs_ot_env%dum(1:k))))
      eval_threshold = 100.0_dp*EPSILON(1.0_dp)*eval_scale
      IF (MINVAL(qs_ot_env%dum(1:k)) <= eval_threshold .AND. &
          preconditioned .AND. ASSOCIATED(qs_ot_env%preconditioner)) THEN
         NULLIFY (qs_ot_env%preconditioner)
         qs_ot_env%preconditioner_changed = .TRUE.
         IF (PRESENT(preconditioner_rejected)) preconditioner_rejected = .TRUE.
         target_re => qs_ot_env%matrix_sc0
         target_im => qs_ot_env%matrix_sc0_im
         CALL qs_ot_complex_multiply('C', 'N', qs_ot_env%matrix_sc0, qs_ot_env%matrix_sc0_im, &
                                     target_re, target_im, qs_ot_env%matrix_os, qs_ot_env%matrix_os_im, &
                                     qs_ot_env%matrix_buf1)
         CALL cp_dbcsr_heevd(matrix_re=qs_ot_env%matrix_os, matrix_im=qs_ot_env%matrix_os_im, &
                             eigenvectors_re=qs_ot_env%matrix_buf3, eigenvectors_im=qs_ot_env%matrix_buf3_im, &
                             eigenvalues=qs_ot_env%dum, para_env=qs_ot_env%para_env, &
                             blacs_env=qs_ot_env%blacs_env)
         eval_scale = MAX(1.0_dp, MAXVAL(ABS(qs_ot_env%dum(1:k))))
         eval_threshold = 100.0_dp*EPSILON(1.0_dp)*eval_scale
      END IF
      CPASSERT(MINVAL(qs_ot_env%dum(1:k)) > eval_threshold)
      DO i = 1, k
         qs_ot_env%dum(i) = 1.0_dp/qs_ot_env%dum(i)
      END DO
      CALL dbcsr_copy(qs_ot_env%matrix_buf1, qs_ot_env%matrix_buf3)
      CALL dbcsr_copy(qs_ot_env%matrix_buf1_im, qs_ot_env%matrix_buf3_im)
      CALL dbcsr_scale_by_vector(qs_ot_env%matrix_buf1, alpha=qs_ot_env%dum, side='right')
      CALL dbcsr_scale_by_vector(qs_ot_env%matrix_buf1_im, alpha=qs_ot_env%dum, side='right')
      CALL qs_ot_complex_multiply('N', 'C', qs_ot_env%matrix_buf1, qs_ot_env%matrix_buf1_im, &
                                  qs_ot_env%matrix_buf3, qs_ot_env%matrix_buf3_im, &
                                  qs_ot_env%matrix_os, qs_ot_env%matrix_os_im, &
                                  qs_ot_env%matrix_buf2)
      ! An identity-metric projection must not replace the cached preconditioned metric.
      qs_ot_env%os_valid = preconditioned

   END SUBROUTINE qs_ot_prepare_complex_tangent_metric

! **************************************************************************************************
!> \brief Prepare uncached single-rank STRICT metrics in bounded OpenMP batches.
!>        Leave unsuccessful channels to the scalar path so preconditioner rejection and history
!>        invalidation retain their existing order. Workers use only private dense storage.
!> \param qs_ot_env local OT channels
! **************************************************************************************************
   SUBROUTINE qs_ot_prepare_complex_tangent_metric_batch(qs_ot_env)
      TYPE(qs_ot_type), INTENT(INOUT)                    :: qs_ot_env(:)

      INTEGER                                            :: first, i, j, k, n, nbatch, nchannels, &
                                                            nworkers
      INTEGER, ALLOCATABLE                               :: channels(:)
      LOGICAL                                            :: handled
      REAL(dp)                                           :: eval_threshold
      TYPE(dbcsr_p_type), ALLOCATABLE                    :: matrix_in(:, :), matrix_out(:, :)
      TYPE(preconditioner_p_type), ALLOCATABLE           :: preconditioners(:)
      TYPE(qs_ot_local_work_type), ALLOCATABLE           :: work(:)

      IF (.NOT. qs_ot_local_channels(qs_ot_env)) RETURN
      IF (.NOT. qs_ot_env(1)%settings%ks .OR. qs_ot_env(1)%settings%ot_algorithm /= 'TOD') RETURN
      IF (diag_lib_explicit .OR. diag_type == FM_DIAG_TYPE_DLAF) RETURN
      nchannels = COUNT(.NOT. qs_ot_env(:)%os_valid)
      IF (nchannels == 0) RETURN
      ALLOCATE (channels(nchannels))
      j = 0
      DO i = 1, SIZE(qs_ot_env)
         IF (qs_ot_env(i)%os_valid) CYCLE
         j = j + 1
         channels(j) = i
      END DO
      nworkers = 1
!$    nworkers = MIN(nchannels, omp_get_max_threads())
      ALLOCATE (work(nworkers), preconditioners(nworkers), &
                matrix_in(2, nworkers), matrix_out(2, nworkers))
      DO j = 1, nworkers
         CALL work(j)%gemm_ctx%create(LOCAL_GEMM_PU_GPU, timing=.FALSE.)
      END DO
      DO first = 1, nchannels, nworkers
         nbatch = MIN(nworkers, nchannels - first + 1)
         DO j = 1, nbatch
            i = channels(first + j - 1)
            preconditioners(j)%preconditioner => qs_ot_env(i)%preconditioner
            matrix_in(1, j)%matrix => qs_ot_env(i)%matrix_sc0
            matrix_in(2, j)%matrix => qs_ot_env(i)%matrix_sc0_im
            matrix_out(:, j) = matrix_in(:, j)
            IF (ASSOCIATED(qs_ot_env(i)%preconditioner)) THEN
               matrix_out(1, j)%matrix => qs_ot_env(i)%matrix_psc0
               matrix_out(2, j)%matrix => qs_ot_env(i)%matrix_psc0_im
            END IF
         END DO
         CALL apply_preconditioner_complex_batch(preconditioners(1:nbatch), &
                                                 matrix_in(:, 1:nbatch), matrix_out(:, 1:nbatch), handled)
         DO j = 1, nbatch
            IF (.NOT. handled .AND. ASSOCIATED(preconditioners(j)%preconditioner)) THEN
               CALL apply_preconditioner(preconditioners(j)%preconditioner, &
                                         matrix_in(1, j)%matrix, matrix_in(2, j)%matrix, &
                                         matrix_out(1, j)%matrix, matrix_out(2, j)%matrix)
            END IF
            CALL dbcsr_get_info(matrix_in(1, j)%matrix, nfullrows_total=n, nfullcols_total=k)
            IF (ALLOCATED(work(j)%sc0)) THEN
               IF (n /= work(j)%n .OR. k /= work(j)%k) THEN
                  DEALLOCATE (work(j)%sc0, work(j)%target, work(j)%os, work(j)%r, &
                              work(j)%buf1, work(j)%buf2, work(j)%vectors, work(j)%values)
               END IF
            END IF
            IF (.NOT. ALLOCATED(work(j)%sc0)) THEN
               ALLOCATE (work(j)%sc0(n, k, 2), work(j)%target(n, k, 2), work(j)%os(k, k, 2), &
                         work(j)%r(k, k, 2), work(j)%buf1(k, k, 2), work(j)%buf2(k, k, 2), &
                         work(j)%vectors(k, k), work(j)%values(k))
               work(j)%n = n
               work(j)%k = k
            END IF
            CALL cp_cfm_local_workspace_prepare(work(j)%solver, k, work(j)%vectors)
            CALL qs_ot_load_local_pair(matrix_in(1, j)%matrix, matrix_in(2, j)%matrix, work(j)%sc0)
            CALL qs_ot_load_local_pair(matrix_out(1, j)%matrix, matrix_out(2, j)%matrix, work(j)%target)
         END DO
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) NUM_THREADS(nbatch) IF(nbatch > 1) &
!$OMP SHARED(work,nbatch) PRIVATE(eval_threshold)
         DO j = 1, nbatch
            CALL qs_ot_complex_multiply_local('C', 'N', work(j)%sc0, work(j)%target, &
                                              work(j)%os, work(j)%buf2(:, :, 1), work(j)%gemm_ctx)
            work(j)%vectors(:, :) = CMPLX(work(j)%os(:, :, 1), work(j)%os(:, :, 2), KIND=dp)
            CALL cp_cfm_heevd_local(work(j)%vectors, info=work(j)%info, ws=work(j)%solver)
            IF (work(j)%info /= 0) CYCLE
            eval_threshold = 100.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(work(j)%solver%evals)))
            IF (MINVAL(work(j)%solver%evals) <= eval_threshold) THEN
               work(j)%info = -1
               CYCLE
            END IF
            work(j)%values(:) = 1.0_dp/work(j)%solver%evals(:)
            work(j)%r(:, :, 1) = REAL(work(j)%vectors, KIND=dp)
            work(j)%r(:, :, 2) = AIMAG(work(j)%vectors)
            CALL qs_ot_spectral_matrix_local(work(j), matrix=work(j)%os)
         END DO
!$OMP END PARALLEL DO
         DO j = 1, nbatch
            IF (work(j)%info /= 0) CYCLE
            i = channels(first + j - 1)
            CALL qs_ot_store_local_pair(qs_ot_env(i)%matrix_os, qs_ot_env(i)%matrix_os_im, work(j)%os)
            qs_ot_env(i)%os_valid = .TRUE.
         END DO
      END DO
      DO j = 1, nworkers
         CALL work(j)%gemm_ctx%destroy()
      END DO
   END SUBROUTINE qs_ot_prepare_complex_tangent_metric_batch

! **************************************************************************************************
!> \brief finite complex STRICT derivative, projected onto C0^H*S*X=0
!> \param matrix_hc real part of H(k)*C(k)
!> \param matrix_hc_im imaginary part of H(k)*C(k)
!> \param qs_ot_env OT channel state
!> \param matrix_hc_rotation occupation-weighted H(k)C(k) for the rotation channel
!> \param matrix_hc_rotation_im imaginary component of matrix_hc_rotation
! **************************************************************************************************
   SUBROUTINE qs_ot_get_derivative_complex(matrix_hc, matrix_hc_im, qs_ot_env, &
                                           matrix_hc_rotation, matrix_hc_rotation_im)
      TYPE(dbcsr_type), POINTER                          :: matrix_hc, matrix_hc_im
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: matrix_hc_rotation, matrix_hc_rotation_im

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative_complex'

      INTEGER                                            :: handle
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_type)                                   :: tmp_nk
      TYPE(dbcsr_type), POINTER                          :: hc_rotation_im, hc_rotation_re, &
                                                            hc_work_im, hc_work_re
      TYPE(dbcsr_type), TARGET                           :: hc_rot_im, hc_rot_re

      CALL timeset(routineN, handle)

      CPASSERT(qs_ot_env%has_complex_kpoint_state)
      CPASSERT(qs_ot_env%settings%ot_algorithm == 'TOD')

      hc_rotation_re => matrix_hc
      hc_rotation_im => matrix_hc_im
      IF (PRESENT(matrix_hc_rotation) .OR. PRESENT(matrix_hc_rotation_im)) THEN
         CPASSERT(PRESENT(matrix_hc_rotation) .AND. PRESENT(matrix_hc_rotation_im))
         CPASSERT(ASSOCIATED(matrix_hc_rotation))
         CPASSERT(ASSOCIATED(matrix_hc_rotation_im))
         hc_rotation_re => matrix_hc_rotation
         hc_rotation_im => matrix_hc_rotation_im
      END IF
      hc_work_re => matrix_hc
      hc_work_im => matrix_hc_im

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL qs_ot_generate_rotation_complex(qs_ot_env)

         ! Reconstruct the unrotated finite STRICT orbitals Q(X).
         CALL qs_ot_complex_multiply('N', 'N', qs_ot_env%matrix_c0, qs_ot_env%matrix_c0_im, &
                                     qs_ot_env%matrix_cosp, qs_ot_env%matrix_cosp_im, &
                                     qs_ot_env%matrix_buf_nk, qs_ot_env%matrix_buf_nk_im, &
                                     qs_ot_env%matrix_tmp_nk)
         CALL qs_ot_complex_multiply('N', 'N', qs_ot_env%matrix_x, qs_ot_env%matrix_x_im, &
                                     qs_ot_env%matrix_sinp, qs_ot_env%matrix_sinp_im, &
                                     qs_ot_env%matrix_gx, qs_ot_env%matrix_gx_im, &
                                     qs_ot_env%matrix_tmp_nk)
         CALL dbcsr_add(qs_ot_env%matrix_buf_nk, qs_ot_env%matrix_gx, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         CALL dbcsr_add(qs_ot_env%matrix_buf_nk_im, qs_ot_env%matrix_gx_im, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

         ! dF/dU = Q(X)^H*G_C for C=Q(X)*U.
         CALL qs_ot_complex_multiply('C', 'N', &
                                     qs_ot_env%matrix_buf_nk, qs_ot_env%matrix_buf_nk_im, &
                                     hc_rotation_re, hc_rotation_im, &
                                     qs_ot_env%rot_mat_dedu, qs_ot_env%rot_mat_dedu_im, &
                                     qs_ot_env%matrix_buf1)
         CALL qs_ot_rot_mat_derivative_complex(qs_ot_env)

         ! The STRICT coordinate sees G_Q=G_C*U^H.
         CALL dbcsr_copy(hc_rot_re, matrix_hc, name="strict_hc_rot_re")
         CALL dbcsr_copy(hc_rot_im, matrix_hc_im, name="strict_hc_rot_im")
         CALL dbcsr_copy(tmp_nk, matrix_hc, name="strict_hc_rot_tmp")
         CALL dbcsr_multiply('N', 'T', 1.0_dp, matrix_hc, qs_ot_env%rot_mat_u, &
                             0.0_dp, hc_rot_re)
         CALL dbcsr_multiply('N', 'T', 1.0_dp, matrix_hc_im, qs_ot_env%rot_mat_u_im, &
                             0.0_dp, tmp_nk)
         CALL dbcsr_add(hc_rot_re, tmp_nk, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

         CALL dbcsr_multiply('N', 'T', 1.0_dp, matrix_hc_im, qs_ot_env%rot_mat_u, &
                             0.0_dp, hc_rot_im)
         CALL dbcsr_multiply('N', 'T', 1.0_dp, matrix_hc, qs_ot_env%rot_mat_u_im, &
                             0.0_dp, tmp_nk)
         CALL dbcsr_add(hc_rot_im, tmp_nk, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
         hc_work_re => hc_rot_re
         hc_work_im => hc_rot_im
      END IF

      ! Direct X contribution, H*C sinc(sqrt(P)).
      CALL qs_ot_complex_multiply('N', 'N', hc_work_re, hc_work_im, &
                                  qs_ot_env%matrix_sinp, qs_ot_env%matrix_sinp_im, &
                                  qs_ot_env%matrix_gx, qs_ot_env%matrix_gx_im, &
                                  qs_ot_env%matrix_tmp_nk)

      CALL qs_ot_strict_frechet(hc_work_re, hc_work_im, qs_ot_env%matrix_x, qs_ot_env%matrix_x_im, &
                                qs_ot_env%matrix_sinp_b, qs_ot_env%matrix_buf3, qs_ot_env%matrix_buf3_im, qs_ot_env)
      CALL qs_ot_strict_frechet(hc_work_re, hc_work_im, qs_ot_env%matrix_c0, qs_ot_env%matrix_c0_im, &
                                qs_ot_env%matrix_cosp_b, qs_ot_env%matrix_buf4, qs_ot_env%matrix_buf4_im, qs_ot_env)
      CALL dbcsr_add(qs_ot_env%matrix_buf3, qs_ot_env%matrix_buf4, &
                     alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      CALL dbcsr_add(qs_ot_env%matrix_buf3_im, qs_ot_env%matrix_buf4_im, &
                     alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      ! Transform back and add the Hermitian adjoint generated by dP.
      CALL qs_ot_complex_multiply('N', 'N', qs_ot_env%matrix_r, qs_ot_env%matrix_r_im, &
                                  qs_ot_env%matrix_buf3, qs_ot_env%matrix_buf3_im, &
                                  qs_ot_env%matrix_buf1, qs_ot_env%matrix_buf1_im, &
                                  qs_ot_env%matrix_buf2)
      CALL qs_ot_complex_multiply('N', 'C', qs_ot_env%matrix_buf1, qs_ot_env%matrix_buf1_im, &
                                  qs_ot_env%matrix_r, qs_ot_env%matrix_r_im, &
                                  qs_ot_env%matrix_buf3, qs_ot_env%matrix_buf3_im, &
                                  qs_ot_env%matrix_buf2)
      CALL dbcsr_get_info(qs_ot_env%matrix_buf3, distribution=dist)
      CALL dbcsr_transposed(qs_ot_env%matrix_buf4, qs_ot_env%matrix_buf3, &
                            shallow_data_copy=.FALSE., use_distribution=dist, &
                            transpose_distribution=.FALSE.)
      CALL dbcsr_transposed(qs_ot_env%matrix_buf4_im, qs_ot_env%matrix_buf3_im, &
                            shallow_data_copy=.FALSE., use_distribution=dist, &
                            transpose_distribution=.FALSE.)
      CALL dbcsr_add(qs_ot_env%matrix_buf3, qs_ot_env%matrix_buf4, &
                     alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      CALL dbcsr_add(qs_ot_env%matrix_buf3_im, qs_ot_env%matrix_buf4_im, &
                     alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      CALL qs_ot_complex_multiply('N', 'N', qs_ot_env%matrix_sx, qs_ot_env%matrix_sx_im, &
                                  qs_ot_env%matrix_buf3, qs_ot_env%matrix_buf3_im, &
                                  qs_ot_env%matrix_buf_nk, qs_ot_env%matrix_buf_nk_im, &
                                  qs_ot_env%matrix_tmp_nk)
      CALL dbcsr_add(qs_ot_env%matrix_gx, qs_ot_env%matrix_buf_nk, &
                     alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      CALL dbcsr_add(qs_ot_env%matrix_gx_im, qs_ot_env%matrix_buf_nk_im, &
                     alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      CALL qs_ot_project_gradient(qs_ot_env, qs_ot_env%matrix_gx, qs_ot_env%matrix_gx_im)

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL dbcsr_release(hc_rot_re)
         CALL dbcsr_release(hc_rot_im)
         CALL dbcsr_release(tmp_nk)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_ot_get_derivative_complex

! **************************************************************************************************
!> \brief Whether all channels can use the single-rank STRICT/IRAC kernels without rotation.
!> \param qs_ot_env local OT channels
!> \return true when every channel is eligible
! **************************************************************************************************
   LOGICAL FUNCTION qs_ot_local_channels(qs_ot_env) RESULT(local)
      TYPE(qs_ot_type), INTENT(IN)                       :: qs_ot_env(:)

      INTEGER                                            :: i

      local = .FALSE.
      IF (SIZE(qs_ot_env) == 0) RETURN
      DO i = 1, SIZE(qs_ot_env)
         IF (.NOT. qs_ot_env(i)%has_complex_kpoint_state) RETURN
         IF (qs_ot_env(i)%settings%do_rotation) RETURN
         SELECT CASE (qs_ot_env(i)%settings%ot_algorithm)
         CASE ('TOD')
            IF (qs_ot_env(i)%settings%do_ener) RETURN
         CASE ('REF')
         CASE DEFAULT
            RETURN
         END SELECT
         IF (qs_ot_env(i)%settings%ot_algorithm /= qs_ot_env(1)%settings%ot_algorithm) RETURN
         IF (qs_ot_env(i)%para_env%num_pe /= 1) RETURN
      END DO
      local = .TRUE.
   END FUNCTION qs_ot_local_channels

! **************************************************************************************************
!> \brief Rebuild single-rank STRICT/IRAC orbitals in bounded batches using cached overlaps.
!>        DBCSR packing and writeback are outside OpenMP; workers only use dense data and private scratch.
!> \param qs_ot_env local OT channels
!> \param matrix_s pointers to cached real/imaginary overlaps, indexed by component and channel
!> \param matrix_c pointers to physical orbitals; only their active columns are updated
!> \param handled false when the original distributed/rotation path is required
!> \param initialize build C0/SC0 from physical orbitals without modifying the physical state
! **************************************************************************************************
   SUBROUTINE qs_ot_get_orbitals_local_batch(qs_ot_env, matrix_s, matrix_c, handled, initialize)
      TYPE(qs_ot_type), INTENT(INOUT)                    :: qs_ot_env(:)
      TYPE(cp_fm_p_type), INTENT(IN)                     :: matrix_s(:, :)
      TYPE(cp_cfm_p_type), INTENT(IN)                    :: matrix_c(:)
      LOGICAL, INTENT(OUT)                               :: handled
      LOGICAL, INTENT(IN), OPTIONAL                      :: initialize

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ot_get_orbitals_local_batch'

      INTEGER                                            :: first, handle, i, j, k, n, nbatch, &
                                                            nchannels, nworkers
      LOGICAL                                            :: initial, ref, update
      TYPE(qs_ot_local_work_type), ALLOCATABLE           :: work(:)

      initial = .FALSE.
      IF (PRESENT(initialize)) initial = initialize
      handled = .FALSE.
      IF (.NOT. qs_ot_local_channels(qs_ot_env)) RETURN
      ref = qs_ot_env(1)%settings%ot_algorithm == 'REF'
      update = .FALSE.
      IF (ref) update = qs_ot_ref_update(qs_ot_env(1))
      CALL timeset(routineN, handle)
      nchannels = SIZE(qs_ot_env)
      nworkers = 1
!$    nworkers = MIN(nchannels, omp_get_max_threads())
      ALLOCATE (work(nworkers))
      DO j = 1, nworkers
         CALL work(j)%gemm_ctx%create(LOCAL_GEMM_PU_GPU, timing=.FALSE.)
      END DO
      DO first = 1, nchannels, nworkers
         nbatch = MIN(nworkers, nchannels - first + 1)
         DO j = 1, nbatch
            i = first + j - 1
            CALL dbcsr_get_info(qs_ot_env(i)%matrix_x, nfullrows_total=n, nfullcols_total=k)
            CALL qs_ot_allocate_local_work(work(j), n, k)
            IF (.NOT. initial .OR. ref) THEN
               IF (ALLOCATED(work(j)%vectors)) THEN
                  IF (SIZE(work(j)%vectors, 1) /= k) DEALLOCATE (work(j)%vectors, work(j)%values)
               END IF
               IF (.NOT. ALLOCATED(work(j)%vectors)) ALLOCATE (work(j)%vectors(k, k), work(j)%values(k))
               CALL cp_cfm_local_workspace_prepare(work(j)%solver, k, work(j)%vectors)
            END IF
            IF (initial) THEN
               work(j)%x(:, :, 1) = REAL(matrix_c(i)%matrix%local_data(1:n, 1:k), KIND=dp)
               work(j)%x(:, :, 2) = AIMAG(matrix_c(i)%matrix%local_data(1:n, 1:k))
            ELSE
               CALL qs_ot_load_local_pair(qs_ot_env(i)%matrix_x, qs_ot_env(i)%matrix_x_im, work(j)%x)
               IF (.NOT. ref) THEN
                  CALL qs_ot_load_local_pair(qs_ot_env(i)%matrix_c0, qs_ot_env(i)%matrix_c0_im, work(j)%c0)
               END IF
            END IF
         END DO
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) NUM_THREADS(nbatch) IF(nbatch > 1) &
!$OMP SHARED(work,matrix_s,matrix_c,nbatch,first,ref,UPDATE,initial) PRIVATE(i,n,k)
         DO j = 1, nbatch
            i = first + j - 1
            CALL qs_ot_orbitals_local(work(j), matrix_s(1, i)%matrix%local_data, &
                                      matrix_s(2, i)%matrix%local_data, ref, update, initial)
            IF (work(j)%info /= 0) CYCLE
            IF (initial) CYCLE
            n = work(j)%n
            k = work(j)%k
            matrix_c(i)%matrix%local_data(1:n, 1:k) = CMPLX(work(j)%out(:, :, 1), work(j)%out(:, :, 2), KIND=dp)
         END DO
!$OMP END PARALLEL DO
         DO j = 1, nbatch
            IF (work(j)%info == -1) CPABORT('Complex REF overlap is not positive definite')
            IF (work(j)%info /= 0) CPABORT('Local complex OT overlap diagonalization failed')
            i = first + j - 1
            IF (initial) THEN
               CALL qs_ot_store_local_pair(qs_ot_env(i)%matrix_c0, qs_ot_env(i)%matrix_c0_im, work(j)%out)
               CALL qs_ot_store_local_pair(qs_ot_env(i)%matrix_sc0, qs_ot_env(i)%matrix_sc0_im, work(j)%buf_nk)
               IF (.NOT. ref) CYCLE
            END IF
            CALL qs_ot_store_local_pair(qs_ot_env(i)%matrix_sx, qs_ot_env(i)%matrix_sx_im, work(j)%sx)
            IF (ref) THEN
               CALL qs_ot_store_local_pair(qs_ot_env(i)%matrix_ref_inv_sqrt, &
                                           qs_ot_env(i)%matrix_ref_inv_sqrt_im, work(j)%sinp)
               IF (update .OR. initial) CALL qs_ot_store_local_pair(qs_ot_env(i)%matrix_x, qs_ot_env(i)%matrix_x_im, work(j)%x)
               CYCLE
            END IF
            CALL qs_ot_store_local_pair(qs_ot_env(i)%matrix_p, qs_ot_env(i)%matrix_p_im, work(j)%buf3)
            CALL qs_ot_store_local_pair(qs_ot_env(i)%matrix_r, qs_ot_env(i)%matrix_r_im, work(j)%r)
            CALL qs_ot_store_local_pair(qs_ot_env(i)%matrix_cosp, qs_ot_env(i)%matrix_cosp_im, work(j)%buf4)
            CALL qs_ot_store_local_pair(qs_ot_env(i)%matrix_sinp, qs_ot_env(i)%matrix_sinp_im, work(j)%sinp)
            CALL copy_local_to_dbcsr(qs_ot_env(i)%matrix_cosp_b, work(j)%cosp_b)
            CALL copy_local_to_dbcsr(qs_ot_env(i)%matrix_sinp_b, work(j)%sinp_b)
            qs_ot_env(i)%evals(:) = work(j)%solver%evals(:)
            qs_ot_env(i)%do_taylor = .FALSE.
         END DO
      END DO
      DO j = 1, nworkers
         CALL work(j)%gemm_ctx%destroy()
      END DO
      handled = .TRUE.
      CALL timestop(handle)
   END SUBROUTINE qs_ot_get_orbitals_local_batch

! **************************************************************************************************
!> \brief Evaluate S*X and reconstruct physical STRICT or IRAC orbitals.
!> \param w private storage and status; buf3 retains P and buf4 retains cos(sqrt(P)) for writeback
!> \param overlap_re cached real overlap, including any leading-dimension padding
!> \param overlap_im cached imaginary overlap
!> \param ref select IRAC instead of STRICT
!> \param update accept the IRAC coordinate at the existing minimizer update point
!> \param initialize retain SC0 and use the analytic zero-coordinate STRICT state
! **************************************************************************************************
   SUBROUTINE qs_ot_orbitals_local(w, overlap_re, overlap_im, ref, update, initialize)
      TYPE(qs_ot_local_work_type), INTENT(INOUT)         :: w
      REAL(dp), CONTIGUOUS, INTENT(IN)                   :: overlap_re(:, :), overlap_im(:, :)
      LOGICAL, INTENT(IN)                                :: ref, update
      LOGICAL, INTENT(IN), OPTIONAL                      :: initialize

      INTEGER                                            :: i

      CALL qs_ot_complex_multiply_local('N', 'N', overlap_re, overlap_im, w%x(:, :, 1), w%x(:, :, 2), &
                                        w%sx(:, :, 1), w%sx(:, :, 2), w%tmp, w%gemm_ctx)
      IF (PRESENT(initialize)) THEN
         IF (initialize) THEN
            ! Keep the incoming SC0, as in the distributed initializer. IRAC may
            ! subsequently normalize C0 and accept a different X/SX pair.
            w%buf_nk(:, :, :) = w%sx(:, :, :)
            IF (.NOT. ref) THEN
               w%out(:, :, :) = w%x(:, :, :)
               w%info = 0
               RETURN
            END IF
         END IF
      END IF
      CALL qs_ot_complex_multiply_local('C', 'N', w%x, w%sx, w%buf3, w%buf2(:, :, 1), w%gemm_ctx)
      w%vectors(:, :) = CMPLX(w%buf3(:, :, 1), w%buf3(:, :, 2), KIND=dp)
      CALL cp_cfm_heevd_local(w%vectors, info=w%info, ws=w%solver)
      IF (w%info /= 0) RETURN
      w%r(:, :, 1) = REAL(w%vectors, KIND=dp)
      w%r(:, :, 2) = AIMAG(w%vectors)
      IF (ref) THEN
         IF (ANY(w%solver%evals <= EPSILON(1.0_dp))) THEN
            w%info = -1
            RETURN
         END IF
         w%values(:) = 1.0_dp/SQRT(w%solver%evals(:))
         CALL qs_ot_spectral_matrix_local(w, matrix=w%sinp)
         CALL qs_ot_complex_multiply_local('N', 'N', w%x, w%sinp, w%out, w%tmp, w%gemm_ctx)
         IF (update) THEN
            w%x(:, :, :) = w%out(:, :, :)
            w%sinp(:, :, :) = 0.0_dp
            DO i = 1, w%k
               w%sinp(i, i, 1) = 1.0_dp
            END DO
            CALL qs_ot_complex_multiply_local('N', 'N', overlap_re, overlap_im, w%x(:, :, 1), w%x(:, :, 2), &
                                              w%sx(:, :, 1), w%sx(:, :, 2), w%tmp, w%gemm_ctx)
         END IF
         RETURN
      END IF
      w%solver%evals(:) = MAX(0.0_dp, w%solver%evals(:))
      CALL qs_ot_spectral_matrix_local(w, .TRUE., w%buf4)
      CALL qs_ot_spectral_matrix_local(w, .FALSE., w%sinp)
      w%values(:) = SQRT(w%solver%evals(:))
      CALL qs_ot_spectral_kernel(w%values, w%values, w%cosp_b, .TRUE.)
      CALL qs_ot_spectral_kernel(w%values, w%values, w%sinp_b, .FALSE.)
      CALL qs_ot_complex_multiply_local('N', 'N', w%c0, w%buf4, w%out, w%tmp, w%gemm_ctx)
      CALL qs_ot_complex_multiply_local('N', 'N', w%x, w%sinp, w%buf_nk, w%tmp, w%gemm_ctx)
      w%out(:, :, :) = w%out(:, :, :) + w%buf_nk(:, :, :)
   END SUBROUTINE qs_ot_orbitals_local

! **************************************************************************************************
!> \brief Reconstruct one spectral matrix function from the local eigenvectors.
!> \param w private eigenvectors, eigenvalues and scratch
!> \param cosine select cos/sinc(sqrt(P)); if absent, use the supplied spectral values
!> \param matrix dense real/imaginary result
! **************************************************************************************************
   SUBROUTINE qs_ot_spectral_matrix_local(w, cosine, matrix)
      TYPE(qs_ot_local_work_type), INTENT(INOUT)         :: w
      LOGICAL, INTENT(IN), OPTIONAL                      :: cosine
      REAL(dp), CONTIGUOUS, INTENT(INOUT)                :: matrix(:, :, :)

      INTEGER                                            :: i

      IF (PRESENT(cosine)) CALL qs_ot_spectral_values(w%solver%evals, w%values, cosine)
      DO i = 1, w%k
         w%buf1(:, i, 1) = w%r(:, i, 1)*w%values(i)
         w%buf1(:, i, 2) = w%r(:, i, 2)*w%values(i)
      END DO
      CALL qs_ot_complex_multiply_local('N', 'C', w%buf1, w%r, matrix, w%buf2(:, :, 1), w%gemm_ctx)
   END SUBROUTINE qs_ot_spectral_matrix_local

! **************************************************************************************************
!> \brief Allocate common dense OT storage only when a worker's matrix dimensions change.
!> \param w worker storage
!> \param n AO count
!> \param k active MO count
! **************************************************************************************************
   SUBROUTINE qs_ot_allocate_local_work(w, n, k)
      TYPE(qs_ot_local_work_type), INTENT(INOUT)         :: w
      INTEGER, INTENT(IN)                                :: n, k

      IF (ALLOCATED(w%x)) THEN
         IF (n == w%n .AND. k == w%k) RETURN
         DEALLOCATE (w%x, w%sx, w%c0, w%out, w%r, w%sinp, w%buf1, w%buf2, w%buf3, w%buf4, &
                     w%buf_nk, w%sinp_b, w%cosp_b, w%tmp)
      END IF
      ALLOCATE (w%x(n, k, 2), w%sx(n, k, 2), w%c0(n, k, 2), w%out(n, k, 2), &
                w%r(k, k, 2), w%sinp(k, k, 2), w%buf1(k, k, 2), w%buf2(k, k, 2), &
                w%buf3(k, k, 2), w%buf4(k, k, 2), w%buf_nk(n, k, 2), &
                w%sinp_b(k, k), w%cosp_b(k, k), w%tmp(n, k))
      w%n = n
      w%k = k
   END SUBROUTINE qs_ot_allocate_local_work

! **************************************************************************************************
!> \brief Write a dense complex pair back to two DBCSR matrices outside OpenMP.
!> \param matrix_re real DBCSR output
!> \param matrix_im imaginary DBCSR output
!> \param dense dense real/imaginary pair
! **************************************************************************************************
   SUBROUTINE qs_ot_store_local_pair(matrix_re, matrix_im, dense)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_re, matrix_im
      REAL(dp), INTENT(IN)                               :: dense(:, :, :)

      CALL copy_local_to_dbcsr(matrix_re, dense(:, :, 1))
      CALL copy_local_to_dbcsr(matrix_im, dense(:, :, 2))
   END SUBROUTINE qs_ot_store_local_pair

! **************************************************************************************************
!> \brief Evaluate independent, single-rank STRICT/IRAC gradients in bounded OpenMP batches.
!>        DBCSR packing, filtering and metric preparation stay outside the worker region.
!> \param qs_ot_env local OT channels
!> \param matrix_hc real orbital gradients
!> \param matrix_hc_im imaginary orbital gradients
!> \param handled false when the original distributed/rotation path is required
!> \param matrix_hc_physical optional physical gradients for occupation preconditioning
!> \param matrix_hc_physical_im imaginary physical gradients
! **************************************************************************************************
   SUBROUTINE qs_ot_get_derivative_local_batch(qs_ot_env, matrix_hc, matrix_hc_im, handled, &
                                               matrix_hc_physical, matrix_hc_physical_im)
      TYPE(qs_ot_type), INTENT(INOUT)                    :: qs_ot_env(:)
      TYPE(dbcsr_p_type), INTENT(IN)                     :: matrix_hc(:), matrix_hc_im(:)
      LOGICAL, INTENT(OUT)                               :: handled
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: matrix_hc_physical(:), &
                                                            matrix_hc_physical_im(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ot_get_derivative_local_batch'

      INTEGER                                            :: first, handle, i, j, nbatch, nchannels, &
                                                            npasses, nworkers, pass
      LOGICAL                                            :: ref
      TYPE(qs_ot_local_work_type), ALLOCATABLE           :: work(:)

      handled = .FALSE.
      IF (.NOT. qs_ot_local_channels(qs_ot_env)) RETURN
      ref = qs_ot_env(1)%settings%ot_algorithm == 'REF'
      npasses = 1
      IF (PRESENT(matrix_hc_physical)) THEN
         IF (.NOT. ref) RETURN
         CPASSERT(PRESENT(matrix_hc_physical_im))
         npasses = 2
      END IF
      nchannels = SIZE(qs_ot_env)

      CALL timeset(routineN, handle)
      nworkers = 1
!$    nworkers = MIN(nchannels, omp_get_max_threads())
      ALLOCATE (work(nworkers))
      DO j = 1, nworkers
         CALL work(j)%gemm_ctx%create(LOCAL_GEMM_PU_GPU, timing=.FALSE.)
      END DO
      DO first = 1, nchannels, nworkers
         nbatch = MIN(nworkers, nchannels - first + 1)
         DO j = 1, nbatch
            i = first + j - 1
            IF (.NOT. ref) CALL qs_ot_prepare_complex_tangent_metric(qs_ot_env(i))
            CALL qs_ot_load_local_work(work(j), qs_ot_env(i), matrix_hc(i)%matrix, matrix_hc_im(i)%matrix)
            work(j)%filtered = ref .AND. qs_ot_env(i)%settings%eps_irac_filter_matrix > 0.0_dp
         END DO
         DO pass = 1, npasses
            IF (pass == 2) THEN
               DO j = 1, nbatch
                  i = first + j - 1
                  CALL qs_ot_load_local_pair(matrix_hc_physical(i)%matrix, matrix_hc_physical_im(i)%matrix, work(j)%hc)
               END DO
            END IF
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) NUM_THREADS(nbatch) IF(nbatch > 1) &
!$OMP SHARED(work,nbatch,ref)
            DO j = 1, nbatch
               IF (ref) THEN
                  CALL qs_ot_complex_multiply_local('N', 'N', work(j)%x, work(j)%sinp, &
                                                    work(j)%c0, work(j)%tmp, work(j)%gemm_ctx)
                  CALL qs_ot_complex_multiply_local('C', 'N', work(j)%c0, work(j)%hc, &
                                                    work(j)%buf3, work(j)%buf2(:, :, 1), work(j)%gemm_ctx)
                  IF (.NOT. work(j)%filtered) CALL qs_ot_ref_gradient_local(work(j))
               ELSE
                  CALL qs_ot_strict_gradient_local(work(j))
               END IF
            END DO
!$OMP END PARALLEL DO
            IF (ANY(work(1:nbatch)%filtered)) THEN
               ! Keep the existing block filter, including independent real/imaginary sparsity.
               DO j = 1, nbatch
                  IF (.NOT. work(j)%filtered) CYCLE
                  i = first + j - 1
                  CALL qs_ot_store_local_pair(qs_ot_env(i)%buf1_k_k_sym, qs_ot_env(i)%buf2_k_k_sym, work(j)%buf3)
                  CALL dbcsr_filter(qs_ot_env(i)%buf1_k_k_sym, qs_ot_env(i)%settings%eps_irac_filter_matrix)
                  CALL dbcsr_filter(qs_ot_env(i)%buf2_k_k_sym, qs_ot_env(i)%settings%eps_irac_filter_matrix)
                  CALL qs_ot_load_local_pair(qs_ot_env(i)%buf1_k_k_sym, qs_ot_env(i)%buf2_k_k_sym, work(j)%buf3)
               END DO
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) NUM_THREADS(nbatch) IF(nbatch > 1) SHARED(work,nbatch)
               DO j = 1, nbatch
                  IF (work(j)%filtered) CALL qs_ot_ref_gradient_local(work(j))
               END DO
!$OMP END PARALLEL DO
            END IF
            DO j = 1, nbatch
               i = first + j - 1
               CALL qs_ot_store_local_pair(qs_ot_env(i)%matrix_gx, qs_ot_env(i)%matrix_gx_im, work(j)%out)
               IF (work(j)%filtered) THEN
                  CALL dbcsr_filter(qs_ot_env(i)%matrix_gx, qs_ot_env(i)%settings%eps_irac_filter_matrix)
                  CALL dbcsr_filter(qs_ot_env(i)%matrix_gx_im, qs_ot_env(i)%settings%eps_irac_filter_matrix)
               END IF
               IF (pass < npasses) THEN
                  CALL dbcsr_copy(qs_ot_env(i)%matrix_preconditioned_gx, qs_ot_env(i)%matrix_gx)
                  CALL dbcsr_copy(qs_ot_env(i)%matrix_preconditioned_gx_im, qs_ot_env(i)%matrix_gx_im)
               END IF
            END DO
         END DO
      END DO
      DO j = 1, nworkers
         CALL work(j)%gemm_ctx%destroy()
      END DO
      handled = .TRUE.
      CALL timestop(handle)
   END SUBROUTINE qs_ot_get_derivative_local_batch

! **************************************************************************************************
!> \brief Complete the IRAC pullback (HC-SX*F*B)*F after forming and optionally filtering B=Q^H*HC.
!> \param w private storage with F in sinp and B in buf3
! **************************************************************************************************
   SUBROUTINE qs_ot_ref_gradient_local(w)
      TYPE(qs_ot_local_work_type), INTENT(INOUT)         :: w

      CALL qs_ot_complex_multiply_local('N', 'N', w%sx, w%sinp, w%buf_nk, w%tmp, w%gemm_ctx)
      CALL qs_ot_complex_multiply_local('N', 'N', w%buf_nk, w%buf3, w%c0, w%tmp, w%gemm_ctx)
      w%c0(:, :, :) = w%hc(:, :, :) - w%c0(:, :, :)
      CALL qs_ot_complex_multiply_local('N', 'N', w%c0, w%sinp, w%out, w%tmp, w%gemm_ctx)
   END SUBROUTINE qs_ot_ref_gradient_local

! **************************************************************************************************
!> \brief Contract one STRICT term with its Frechet kernel in the eigenvector basis.
!> \param matrix_hc real orbital gradient
!> \param matrix_hc_im imaginary orbital gradient
!> \param orbital real X or C0
!> \param orbital_im imaginary X or C0
!> \param kernel real sinc or cos divided differences
!> \param RESULT real contraction
!> \param result_im imaginary contraction
!> \param qs_ot_env channel eigenvectors and scratch
! **************************************************************************************************
   SUBROUTINE qs_ot_strict_frechet(matrix_hc, matrix_hc_im, orbital, orbital_im, kernel, &
                                   RESULT, result_im, qs_ot_env)
      TYPE(dbcsr_type)                                   :: matrix_hc, matrix_hc_im, orbital, &
                                                            orbital_im, kernel, RESULT, result_im
      TYPE(qs_ot_type)                                   :: qs_ot_env

      CALL qs_ot_complex_multiply('C', 'N', matrix_hc, matrix_hc_im, orbital, orbital_im, &
                                  qs_ot_env%matrix_buf2, qs_ot_env%matrix_buf2_im, qs_ot_env%matrix_buf1)
      CALL qs_ot_complex_multiply('N', 'N', qs_ot_env%matrix_buf2, qs_ot_env%matrix_buf2_im, &
                                  qs_ot_env%matrix_r, qs_ot_env%matrix_r_im, &
                                  qs_ot_env%matrix_buf1, qs_ot_env%matrix_buf1_im, qs_ot_env%matrix_buf4)
      CALL qs_ot_complex_multiply('C', 'N', qs_ot_env%matrix_r, qs_ot_env%matrix_r_im, &
                                  qs_ot_env%matrix_buf1, qs_ot_env%matrix_buf1_im, &
                                  qs_ot_env%matrix_buf2, qs_ot_env%matrix_buf2_im, qs_ot_env%matrix_buf4)
      CALL dbcsr_hadamard_product(qs_ot_env%matrix_buf2, kernel, RESULT)
      CALL dbcsr_hadamard_product(qs_ot_env%matrix_buf2_im, kernel, result_im)
   END SUBROUTINE qs_ot_strict_frechet

! **************************************************************************************************
!> \brief Reuse one worker's scratch and pack the local kernel inputs.
!> \param w worker storage
!> \param qs_ot_env channel with an already prepared tangent metric
!> \param matrix_hc real orbital gradient
!> \param matrix_hc_im imaginary orbital gradient
! **************************************************************************************************
   SUBROUTINE qs_ot_load_local_work(w, qs_ot_env, matrix_hc, matrix_hc_im)
      TYPE(qs_ot_local_work_type), INTENT(INOUT)         :: w
      TYPE(qs_ot_type), INTENT(IN)                       :: qs_ot_env
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_hc, matrix_hc_im

      INTEGER                                            :: k, n, o

      CALL dbcsr_get_info(qs_ot_env%matrix_x, nfullrows_total=n, nfullcols_total=k)
      CALL qs_ot_allocate_local_work(w, n, k)
      IF (ALLOCATED(w%hc)) THEN
         IF (SIZE(w%hc, 1) /= n .OR. SIZE(w%hc, 2) /= k) DEALLOCATE (w%hc)
      END IF
      IF (.NOT. ALLOCATED(w%hc)) ALLOCATE (w%hc(n, k, 2))
      CALL qs_ot_load_local_pair(matrix_hc, matrix_hc_im, w%hc)
      CALL qs_ot_load_local_pair(qs_ot_env%matrix_x, qs_ot_env%matrix_x_im, w%x)
      CALL qs_ot_load_local_pair(qs_ot_env%matrix_sx, qs_ot_env%matrix_sx_im, w%sx)
      IF (qs_ot_env%settings%ot_algorithm == 'REF') THEN
         CALL qs_ot_load_local_pair(qs_ot_env%matrix_ref_inv_sqrt, qs_ot_env%matrix_ref_inv_sqrt_im, w%sinp)
         RETURN
      END IF
      CALL dbcsr_get_info(qs_ot_env%matrix_sc0, nfullcols_total=o)
      IF (ALLOCATED(w%sc0)) THEN
         IF (SIZE(w%sc0, 1) /= n .OR. SIZE(w%buf1_ortho, 2) /= k .OR. o /= SIZE(w%sc0, 2)) THEN
            DEALLOCATE (w%sc0, w%target, w%os, w%buf1_ortho, w%buf2_ortho, w%tmp_ortho)
         END IF
      END IF
      IF (.NOT. ALLOCATED(w%sc0)) THEN
         ALLOCATE (w%sc0(n, o, 2), w%target(n, o, 2), w%os(o, o, 2), &
                   w%buf1_ortho(o, k, 2), w%buf2_ortho(o, k, 2), w%tmp_ortho(o, k))
      END IF
      CALL qs_ot_load_local_pair(qs_ot_env%matrix_c0, qs_ot_env%matrix_c0_im, w%c0)
      CALL qs_ot_load_local_pair(qs_ot_env%matrix_r, qs_ot_env%matrix_r_im, w%r)
      CALL qs_ot_load_local_pair(qs_ot_env%matrix_sinp, qs_ot_env%matrix_sinp_im, w%sinp)
      CALL qs_ot_load_local_pair(qs_ot_env%matrix_sc0, qs_ot_env%matrix_sc0_im, w%sc0)
      CALL qs_ot_load_local_pair(qs_ot_env%matrix_os, qs_ot_env%matrix_os_im, w%os)
      CALL copy_dbcsr_to_local(qs_ot_env%matrix_sinp_b, w%sinp_b)
      CALL copy_dbcsr_to_local(qs_ot_env%matrix_cosp_b, w%cosp_b)
      IF (ASSOCIATED(qs_ot_env%preconditioner)) THEN
         CALL qs_ot_load_local_pair(qs_ot_env%matrix_psc0, qs_ot_env%matrix_psc0_im, w%target)
      ELSE
         w%target(:, :, :) = w%sc0(:, :, :)
      END IF
   END SUBROUTINE qs_ot_load_local_work

! **************************************************************************************************
!> \brief Pack a complex pair, allowing different sparse patterns in its two components.
!> \param matrix_re real DBCSR matrix
!> \param matrix_im imaginary DBCSR matrix
!> \param dense dense real/imaginary pair
! **************************************************************************************************
   SUBROUTINE qs_ot_load_local_pair(matrix_re, matrix_im, dense)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_re, matrix_im
      REAL(dp), INTENT(OUT)                              :: dense(:, :, :)

      CALL copy_dbcsr_to_local(matrix_re, dense(:, :, 1))
      CALL copy_dbcsr_to_local(matrix_im, dense(:, :, 2))
   END SUBROUTINE qs_ot_load_local_pair

! **************************************************************************************************
!> \brief Multiply dense real/imaginary pairs in the same four-product order as DBCSR.
!> \param op_a N or C
!> \param op_b N or C
!> \param a left operand
!> \param b right operand
!> \param c result
!> \param tmp real workspace, at least as large as one component of C
!> \param gemm_ctx caller-owned GEMM context
! **************************************************************************************************
   SUBROUTINE qs_ot_complex_multiply_pair_local(op_a, op_b, a, b, c, tmp, gemm_ctx)
      CHARACTER(LEN=1), INTENT(IN)                       :: op_a, op_b
      REAL(dp), CONTIGUOUS, INTENT(IN)                   :: a(:, :, :), b(:, :, :)
      REAL(dp), CONTIGUOUS, INTENT(INOUT)                :: c(:, :, :), tmp(:, :)
      TYPE(local_gemm_ctxt_type)                         :: gemm_ctx

      CALL qs_ot_complex_multiply_parts_local(op_a, op_b, a(:, :, 1), a(:, :, 2), &
                                              b(:, :, 1), b(:, :, 2), c(:, :, 1), c(:, :, 2), tmp, gemm_ctx)
   END SUBROUTINE qs_ot_complex_multiply_pair_local

! **************************************************************************************************
!> \brief Multiply separate real/imaginary arrays, including cached overlap matrices.
!> \param op_a N or C
!> \param op_b N or C
!> \param a_re real left operand
!> \param a_im imaginary left operand
!> \param b_re real right operand
!> \param b_im imaginary right operand
!> \param c_re real result
!> \param c_im imaginary result
!> \param tmp real workspace
!> \param gemm_ctx caller-owned GEMM context
! **************************************************************************************************
   SUBROUTINE qs_ot_complex_multiply_parts_local(op_a, op_b, a_re, a_im, b_re, b_im, &
                                                 c_re, c_im, tmp, gemm_ctx)
      CHARACTER(LEN=1), INTENT(IN)                       :: op_a, op_b
      REAL(dp), CONTIGUOUS, INTENT(IN)                   :: a_re(:, :), a_im(:, :), b_re(:, :), &
                                                            b_im(:, :)
      REAL(dp), CONTIGUOUS, INTENT(INOUT)                :: c_re(:, :), c_im(:, :), tmp(:, :)
      TYPE(local_gemm_ctxt_type)                         :: gemm_ctx

      CHARACTER(LEN=1)                                   :: trans_a, trans_b
      INTEGER                                            :: k, m, n
      REAL(dp)                                           :: sign_a, sign_b

      CPASSERT(op_a == 'N' .OR. op_a == 'C')
      CPASSERT(op_b == 'N' .OR. op_b == 'C')
      trans_a = MERGE('N', 'T', op_a == 'N')
      trans_b = MERGE('N', 'T', op_b == 'N')
      sign_a = MERGE(1.0_dp, -1.0_dp, op_a == 'N')
      sign_b = MERGE(1.0_dp, -1.0_dp, op_b == 'N')
      m = SIZE(c_re, 1)
      n = SIZE(c_re, 2)
      k = MERGE(SIZE(a_re, 2), SIZE(a_re, 1), op_a == 'N')
      CALL gemm_ctx%gemm(trans_a, trans_b, m, n, k, 1.0_dp, a_re, SIZE(a_re, 1), &
                         b_re, SIZE(b_re, 1), 0.0_dp, c_re, SIZE(c_re, 1))
      CALL gemm_ctx%gemm(trans_a, trans_b, m, n, k, 1.0_dp, a_im, SIZE(a_im, 1), &
                         b_im, SIZE(b_im, 1), 0.0_dp, tmp, SIZE(tmp, 1))
      c_re(:, :) = c_re(:, :) - sign_a*sign_b*tmp(1:m, 1:n)
      CALL gemm_ctx%gemm(trans_a, trans_b, m, n, k, sign_b, a_re, SIZE(a_re, 1), &
                         b_im, SIZE(b_im, 1), 0.0_dp, c_im, SIZE(c_im, 1))
      CALL gemm_ctx%gemm(trans_a, trans_b, m, n, k, sign_a, a_im, SIZE(a_im, 1), &
                         b_re, SIZE(b_re, 1), 0.0_dp, tmp, SIZE(tmp, 1))
      c_im(:, :) = c_im(:, :) + tmp(1:m, 1:n)
   END SUBROUTINE qs_ot_complex_multiply_parts_local

! **************************************************************************************************
!> \brief Contract one STRICT term in the local eigenvector basis.
!> \param w channel matrices and scratch
!> \param orbital dense X or C0 pair
!> \param kernel real sinc or cos divided differences
!> \param RESULT contraction in the eigenvector basis
! **************************************************************************************************
   SUBROUTINE qs_ot_strict_frechet_local(w, orbital, kernel, RESULT)
      TYPE(qs_ot_local_work_type)                        :: w
      REAL(dp), CONTIGUOUS, INTENT(IN)                   :: orbital(:, :, :), kernel(:, :)
      REAL(dp), CONTIGUOUS, INTENT(INOUT)                :: result(:, :, :)

      CALL qs_ot_complex_multiply_local('C', 'N', w%hc, orbital, w%buf2, w%buf1(:, :, 1), w%gemm_ctx)
      CALL qs_ot_complex_multiply_local('N', 'N', w%buf2, w%r, w%buf1, w%buf4(:, :, 1), w%gemm_ctx)
      CALL qs_ot_complex_multiply_local('C', 'N', w%r, w%buf1, w%buf2, w%buf4(:, :, 1), w%gemm_ctx)
      RESULT(:, :, 1) = w%buf2(:, :, 1)*kernel(:, :)
      RESULT(:, :, 2) = w%buf2(:, :, 2)*kernel(:, :)
   END SUBROUTINE qs_ot_strict_frechet_local

! **************************************************************************************************
!> \brief Evaluate and project one local STRICT gradient without allocating or calling DBCSR.
!> \param w channel matrices and scratch, with a prepared inverse tangent metric
! **************************************************************************************************
   SUBROUTINE qs_ot_strict_gradient_local(w)
      TYPE(qs_ot_local_work_type)                        :: w

      INTEGER                                            :: i, j
      REAL(dp)                                           :: im, re

      ! Direct term, followed by the X*sinc and C0*cos Frechet contributions.
      CALL qs_ot_complex_multiply_local('N', 'N', w%hc, w%sinp, w%out, w%tmp, w%gemm_ctx)
      CALL qs_ot_strict_frechet_local(w, w%x, w%sinp_b, w%buf3)
      CALL qs_ot_strict_frechet_local(w, w%c0, w%cosp_b, w%buf4)
      w%buf3(:, :, :) = w%buf3(:, :, :) + w%buf4(:, :, :)
      CALL qs_ot_complex_multiply_local('N', 'N', w%r, w%buf3, w%buf1, w%buf2(:, :, 1), w%gemm_ctx)
      CALL qs_ot_complex_multiply_local('N', 'C', w%buf1, w%r, w%buf3, w%buf2(:, :, 1), w%gemm_ctx)

      ! Merge each conjugate pair in place, without a transpose temporary.
      DO j = 1, w%k
         w%buf3(j, j, 1) = w%buf3(j, j, 1) + w%buf3(j, j, 1)
         w%buf3(j, j, 2) = 0.0_dp
         DO i = j + 1, w%k
            re = w%buf3(i, j, 1) + w%buf3(j, i, 1)
            im = w%buf3(i, j, 2) - w%buf3(j, i, 2)
            w%buf3(i, j, 1) = re
            w%buf3(j, i, 1) = re
            w%buf3(i, j, 2) = im
            w%buf3(j, i, 2) = -im
         END DO
      END DO
      CALL qs_ot_complex_multiply_local('N', 'N', w%sx, w%buf3, w%buf_nk, w%tmp, w%gemm_ctx)
      w%out(:, :, :) = w%out(:, :, :) + w%buf_nk(:, :, :)

      ! Tangent projection with the same preconditioned reference and inverse metric.
      CALL qs_ot_complex_multiply_local('C', 'N', w%target, w%out, w%buf1_ortho, w%tmp_ortho, w%gemm_ctx)
      CALL qs_ot_complex_multiply_local('N', 'N', w%os, w%buf1_ortho, w%buf2_ortho, w%tmp_ortho, w%gemm_ctx)
      CALL qs_ot_complex_multiply_local('N', 'N', w%sc0, w%buf2_ortho, w%buf_nk, w%tmp, w%gemm_ctx)
      w%out(:, :, :) = w%out(:, :, :) - w%buf_nk(:, :, :)
   END SUBROUTINE qs_ot_strict_gradient_local

END MODULE qs_ot_kp

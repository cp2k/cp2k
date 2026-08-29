!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief orbital transformations
!> \par History
!>      Added Taylor expansion based computation of the matrix functions (01.2004)
!>      added additional rotation variables for non-equivalent occupied orbs (08.2004)
!> \author Joost VandeVondele (06.2002)
! **************************************************************************************************
MODULE qs_ot
   USE arnoldi_api,                     ONLY: arnoldi_extremal
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_distribution_type, dbcsr_filter, &
        dbcsr_finalize, dbcsr_get_block_p, dbcsr_get_info, dbcsr_get_occupation, dbcsr_init_p, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_release, dbcsr_release_p, &
        dbcsr_reserve_blocks, dbcsr_scale, dbcsr_set, dbcsr_transposed, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert,&
                                              cp_dbcsr_cholesky_restore
   USE cp_dbcsr_contrib,                ONLY: dbcsr_add_on_diag,&
                                              dbcsr_frobenius_norm,&
                                              dbcsr_gershgorin_norm,&
                                              dbcsr_hadamard_product,&
                                              dbcsr_scale_by_vector
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_heevd,&
                                              cp_dbcsr_syevd
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: diag_complex
   USE message_passing,                 ONLY: mp_comm_type
   USE preconditioner,                  ONLY: apply_preconditioner
   USE preconditioner_types,            ONLY: preconditioner_type
   USE qs_ot_types,                     ONLY: qs_ot_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC  :: qs_ot_get_p
   PUBLIC  :: qs_ot_get_p_complex
   PUBLIC  :: qs_ot_get_orbitals
   PUBLIC  :: qs_ot_get_orbitals_complex
   PUBLIC  :: qs_ot_get_derivative
   PUBLIC  :: qs_ot_get_derivative_complex
   PUBLIC  :: qs_ot_get_orbitals_ref
   PUBLIC  :: qs_ot_get_orbitals_ref_complex
   PUBLIC  :: qs_ot_get_derivative_ref
   PUBLIC  :: qs_ot_get_derivative_ref_complex
   PUBLIC  :: qs_ot_apply_complex_frechet_dbcsr
   PUBLIC  :: qs_ot_complex_exp_frechet_kernel
   PUBLIC  :: qs_ot_fixed_n_energy_gradient
   PUBLIC  :: qs_ot_fixed_n_energy_hessian
   PUBLIC  :: qs_ot_fixed_n_projector_frechet
   PUBLIC  :: qs_ot_fixed_n_response_mu_shift
   PUBLIC  :: qs_ot_generate_rotation
   PUBLIC  :: qs_ot_generate_rotation_complex
   PUBLIC  :: qs_ot_rot_mat_derivative
   PUBLIC  :: qs_ot_rot_mat_derivative_complex
   PUBLIC  :: qs_ot_new_preconditioner
   PRIVATE :: qs_ot_p2m_diag
   PRIVATE :: qs_ot_p2m_diag_complex
   PRIVATE :: qs_ot_complex_multiply
   PRIVATE :: qs_ot_sinc
   PRIVATE :: qs_ot_ref_poly
   PRIVATE :: qs_ot_ref_chol
   PRIVATE :: qs_ot_ref_lwdn
   PRIVATE :: qs_ot_ref_decide
   PRIVATE :: qs_ot_ref_update
   PRIVATE :: qs_ot_refine
   PRIVATE :: qs_ot_on_the_fly_localize

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot'

CONTAINS

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

      denominator = MAX(0.0_dp, fixed_n_curvature_sum)
      IF (denominator <= EPSILON(denominator)) denominator = MAX(0.0_dp, local_curvature_sum)
      mu_shift = 0.0_dp
      IF (denominator > EPSILON(denominator)) mu_shift = weighted_energy_response/denominator

   END FUNCTION qs_ot_fixed_n_response_mu_shift

! **************************************************************************************************
!> \brief fixed-N Mermin gradient in auxiliary-energy coordinates
!> \param rayleigh_energy diagonal expectation values of the current Hamiltonian
!> \param energy_coordinate auxiliary band energies controlling the occupations
!> \param response_weight positive weighted occupation susceptibilities chi_i
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
      IF (fixed_n_weight_sum > EPSILON(fixed_n_weight_sum)) THEN
         fixed_n_mean = fixed_n_weighted_residual/fixed_n_weight_sum
      END IF
      gradient(:) = response_weight(:)* &
                    (fixed_n_mean - (rayleigh_energy(:) - energy_coordinate(:)))

   END SUBROUTINE qs_ot_fixed_n_energy_gradient

! **************************************************************************************************
!> \brief dense fixed-N occupation Hessian in auxiliary-energy coordinates
!>
!>        H = diag(chi) - chi chi^T / sum(chi) is symmetric positive semidefinite and has the
!>        constant-energy gauge as an exact null vector.
!> \param response_weight positive weighted occupation susceptibilities chi_i
!> \param fixed_n_weight_sum global sum_i chi_i for the fixed-N group
!> \param hessian projected fixed-N Hessian
! **************************************************************************************************
   PURE SUBROUTINE qs_ot_fixed_n_energy_hessian(response_weight, fixed_n_weight_sum, hessian)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: response_weight
      REAL(KIND=dp), INTENT(IN)                          :: fixed_n_weight_sum
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: hessian

      INTEGER                                            :: i, j, n

      n = SIZE(response_weight)
      hessian(:, :) = 0.0_dp
      DO i = 1, n
         hessian(i, i) = MAX(0.0_dp, response_weight(i))
      END DO
      IF (fixed_n_weight_sum > EPSILON(fixed_n_weight_sum)) THEN
         DO j = 1, n
            DO i = 1, n
               hessian(i, j) = hessian(i, j) - &
                               MAX(0.0_dp, response_weight(i))* &
                               MAX(0.0_dp, response_weight(j))/fixed_n_weight_sum
            END DO
         END DO
      END IF
      hessian(:, :) = 0.5_dp*(hessian + TRANSPOSE(hessian))

   END SUBROUTINE qs_ot_fixed_n_energy_hessian

! **************************************************************************************************
!> \brief fixed-N Frechet derivative of a smooth occupation projector
!>
!>        The spectral divided-difference kernel is invariant under rotations inside a degenerate
!>        eigenspace. Its diagonal includes the chemical-potential response of the complete fixed-N
!>        group, while off-diagonal terms describe the physical change of the spectral projector.
!> \param chc projected Hermitian Hamiltonian
!> \param dchc Hermitian Hamiltonian perturbation
!> \param occupation canonical occupations associated with the eigenvalues of chc
!> \param kpoint_weight irreducible-k-point weight
!> \param response_weight positive weighted occupation susceptibilities for this channel
!> \param fixed_n_weight_sum susceptibility summed over the complete fixed-N group
!> \param projector_derivative derivative of the weighted occupation projector
!> \param density_factor optional representation-dependent density prefactor
! **************************************************************************************************
   SUBROUTINE qs_ot_fixed_n_projector_frechet( &
      chc, dchc, occupation, kpoint_weight, response_weight, fixed_n_weight_sum, &
      projector_derivative, density_factor)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN)      :: chc, dchc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: occupation
      REAL(KIND=dp), INTENT(IN)                          :: kpoint_weight
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: response_weight
      REAL(KIND=dp), INTENT(IN)                          :: fixed_n_weight_sum
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: projector_derivative
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: density_factor

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: eigenvectors, kernel, work
      INTEGER                                            :: i, j, n
      REAL(KIND=dp)                                      :: coefficient, denominator, factor, &
                                                            gap_tolerance, local_weight_sum, &
                                                            mu_numerator, mu_shift, scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, weighted_occupation

      n = SIZE(chc, 1)
      CPASSERT(n > 0)
      CPASSERT(ALL(SHAPE(chc) == [n, n]))
      CPASSERT(ALL(SHAPE(dchc) == [n, n]))
      CPASSERT(SIZE(occupation) == n)
      CPASSERT(SIZE(response_weight) == n)
      CPASSERT(ALL(SHAPE(projector_derivative) == [n, n]))

      factor = 1.0_dp
      IF (PRESENT(density_factor)) factor = density_factor
      projector_derivative(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      IF (factor == 0.0_dp .OR. kpoint_weight <= 0.0_dp) RETURN

      ALLOCATE (eigenvectors(n, n), kernel(n, n), work(n, n), eigenvalues(n), &
                weighted_occupation(n))
      CALL diag_complex(chc, eigenvectors, eigenvalues)
      work(:, :) = MATMUL(CONJG(TRANSPOSE(eigenvectors)), MATMUL(dchc, eigenvectors))
      weighted_occupation(:) = factor*kpoint_weight*occupation(:)

      local_weight_sum = SUM(MAX(0.0_dp, response_weight(:)))
      mu_numerator = 0.0_dp
      DO i = 1, n
         mu_numerator = mu_numerator + MAX(0.0_dp, response_weight(i))* &
                        REAL(work(i, i), KIND=dp)
      END DO
      mu_shift = qs_ot_fixed_n_response_mu_shift(mu_numerator, local_weight_sum, &
                                                 fixed_n_weight_sum)

      scale = MAX(1.0_dp, MAXVAL(ABS(eigenvalues(:))))
      gap_tolerance = SQRT(EPSILON(1.0_dp))*scale
      kernel(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO j = 1, n
         DO i = 1, n
            IF (i == j) THEN
               coefficient = -factor*MAX(0.0_dp, response_weight(i))
               kernel(i, i) = CMPLX(coefficient*(REAL(work(i, i), KIND=dp) - mu_shift), &
                                    0.0_dp, KIND=dp)
            ELSE
               denominator = eigenvalues(i) - eigenvalues(j)
               IF (ABS(denominator) > gap_tolerance) THEN
                  coefficient = (weighted_occupation(i) - weighted_occupation(j))/denominator
               ELSE
                  coefficient = -factor*SQRT(MAX(0.0_dp, response_weight(i))* &
                                             MAX(0.0_dp, response_weight(j)))
               END IF
               kernel(i, j) = coefficient*work(i, j)
            END IF
         END DO
      END DO

      projector_derivative(:, :) = MATMUL(eigenvectors, &
                                          MATMUL(kernel, CONJG(TRANSPOSE(eigenvectors))))
      projector_derivative(:, :) = 0.5_dp*(projector_derivative + &
                                           CONJG(TRANSPOSE(projector_derivative)))

      DEALLOCATE (eigenvectors, kernel, work, eigenvalues, weighted_occupation)

   END SUBROUTINE qs_ot_fixed_n_projector_frechet

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
!> \brief gets ready to use the preconditioner/ or renew the preconditioner
!>        only keeps a pointer to the preconditioner.
!>        If you change the preconditioner, you have to call this routine
!>        you remain responsible of proper deallocate of your preconditioner
!>        (or you can reuse it on the next step of the computation)
!> \param qs_ot_env ...
!> \param preconditioner ...
! **************************************************************************************************
   SUBROUTINE qs_ot_new_preconditioner(qs_ot_env, preconditioner)
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(preconditioner_type), POINTER                 :: preconditioner

      INTEGER                                            :: ncoef

      qs_ot_env%preconditioner => preconditioner
      qs_ot_env%os_valid = .FALSE.
      IF (.NOT. ASSOCIATED(qs_ot_env%matrix_psc0)) THEN
         CALL dbcsr_init_p(qs_ot_env%matrix_psc0)
         CALL dbcsr_copy(qs_ot_env%matrix_psc0, qs_ot_env%matrix_sc0, 'matrix_psc0')
      END IF
      IF (qs_ot_env%has_complex_kpoint_state .AND. &
          .NOT. ASSOCIATED(qs_ot_env%matrix_psc0_im)) THEN
         CALL dbcsr_init_p(qs_ot_env%matrix_psc0_im)
         CALL dbcsr_copy(qs_ot_env%matrix_psc0_im, qs_ot_env%matrix_sc0_im, 'matrix_psc0_im')
      END IF

      IF (.NOT. qs_ot_env%use_dx) THEN
         qs_ot_env%use_dx = .TRUE.
         CALL dbcsr_init_p(qs_ot_env%matrix_dx)
         CALL dbcsr_copy(qs_ot_env%matrix_dx, qs_ot_env%matrix_gx, 'matrix_dx')
         IF (qs_ot_env%has_complex_kpoint_state) THEN
            CALL dbcsr_init_p(qs_ot_env%matrix_dx_im)
            CALL dbcsr_copy(qs_ot_env%matrix_dx_im, qs_ot_env%matrix_gx_im, 'matrix_dx_im')
         END IF
         IF (qs_ot_env%settings%do_rotation) THEN
            CALL dbcsr_init_p(qs_ot_env%rot_mat_dx)
            CALL dbcsr_copy(qs_ot_env%rot_mat_dx, qs_ot_env%rot_mat_gx, 'rot_mat_dx')
            IF (qs_ot_env%has_complex_kpoint_state) THEN
               CALL dbcsr_init_p(qs_ot_env%rot_mat_dx_im)
               CALL dbcsr_copy(qs_ot_env%rot_mat_dx_im, qs_ot_env%rot_mat_gx_im, 'rot_mat_dx_im')
            END IF
         END IF
         IF (qs_ot_env%settings%do_ener) THEN
            ncoef = SIZE(qs_ot_env%ener_gx)
            ALLOCATE (qs_ot_env%ener_dx(ncoef))
            qs_ot_env%ener_dx = 0.0_dp
         END IF
      END IF

   END SUBROUTINE qs_ot_new_preconditioner

! **************************************************************************************************
!> \brief multiply paired real/imaginary DBCSR matrices, with optional conjugate transposes
!> \param op_a N or C
!> \param op_b N or C
!> \param a_re real part of A
!> \param a_im imaginary part of A
!> \param b_re real part of B
!> \param b_im imaginary part of B
!> \param c_re real part of A*B
!> \param c_im imaginary part of A*B
!> \param tmp real workspace shaped like C
! **************************************************************************************************
   SUBROUTINE qs_ot_complex_multiply(op_a, op_b, a_re, a_im, b_re, b_im, c_re, c_im, tmp)
      CHARACTER(LEN=1), INTENT(IN)                       :: op_a, op_b
      TYPE(dbcsr_type)                                   :: a_re, a_im, b_re, b_im, c_re, c_im, tmp

      CHARACTER(LEN=1)                                   :: db_op_a, db_op_b
      REAL(KIND=dp)                                      :: sign_a, sign_b

      SELECT CASE (op_a)
      CASE ('N')
         db_op_a = 'N'
         sign_a = 1.0_dp
      CASE ('C')
         db_op_a = 'T'
         sign_a = -1.0_dp
      CASE DEFAULT
         CPABORT("Complex matrix product expects N or C for op_a")
      END SELECT
      SELECT CASE (op_b)
      CASE ('N')
         db_op_b = 'N'
         sign_b = 1.0_dp
      CASE ('C')
         db_op_b = 'T'
         sign_b = -1.0_dp
      CASE DEFAULT
         CPABORT("Complex matrix product expects N or C for op_b")
      END SELECT

      CALL dbcsr_multiply(db_op_a, db_op_b, 1.0_dp, a_re, b_re, 0.0_dp, c_re)
      CALL dbcsr_multiply(db_op_a, db_op_b, 1.0_dp, a_im, b_im, 0.0_dp, tmp)
      CALL dbcsr_add(c_re, tmp, alpha_scalar=1.0_dp, beta_scalar=-sign_a*sign_b)

      CALL dbcsr_multiply(db_op_a, db_op_b, sign_b, a_re, b_im, 0.0_dp, c_im)
      CALL dbcsr_multiply(db_op_a, db_op_b, sign_a, a_im, b_re, 0.0_dp, tmp)
      CALL dbcsr_add(c_im, tmp, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

   END SUBROUTINE qs_ot_complex_multiply

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param C_NEW ...
!> \param SC ...
!> \param G_OLD ...
!> \param D ...
! **************************************************************************************************
   SUBROUTINE qs_ot_on_the_fly_localize(qs_ot_env, C_NEW, SC, G_OLD, D)
      !
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type), POINTER                          :: C_NEW, SC, G_OLD, D

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_on_the_fly_localize'
      INTEGER, PARAMETER                                 :: taylor_order = 50
      REAL(KIND=dp), PARAMETER                           :: alpha = 0.1_dp, f2_eps = 0.01_dp

      INTEGER                                            :: col, col_size, handle, i, k, n, p, row, &
                                                            row_size
      REAL(dp), DIMENSION(:, :), POINTER                 :: block
      REAL(KIND=dp)                                      :: expfactor, f2, norm_fro, norm_gct, tmp
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type), POINTER                          :: C, Gp1, Gp2, GU, U
      TYPE(mp_comm_type)                                 :: group

      CALL timeset(routineN, handle)
      !
      !
      CALL dbcsr_get_info(C_NEW, nfullrows_total=n, nfullcols_total=k)
      !
      ! C = C*expm(-G)
      GU => qs_ot_env%buf1_k_k_nosym ! a buffer
      U => qs_ot_env%buf2_k_k_nosym ! a buffer
      Gp1 => qs_ot_env%buf3_k_k_nosym ! a buffer
      Gp2 => qs_ot_env%buf4_k_k_nosym ! a buffer
      C => qs_ot_env%buf1_n_k ! a buffer
      !
      ! compute the derivative of the norm
      !-------------------------------------------------------------------
      ! (x^2+eps)^1/2
      f2 = 0.0_dp
      CALL dbcsr_copy(C, C_NEW)
      CALL dbcsr_iterator_start(iter, C)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, row_size=row_size, col_size=col_size)
         DO p = 1, col_size ! p
         DO i = 1, row_size ! i
            tmp = SQRT(block(i, p)**2 + f2_eps)
            f2 = f2 + tmp
            block(i, p) = block(i, p)/tmp
         END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL dbcsr_get_info(C, group=group)
      CALL group%sum(f2)
      !
      !
      CALL dbcsr_multiply('T', 'N', 1.0_dp, C, C_NEW, 0.0_dp, GU)
      !
      ! antisymetrize
      CALL dbcsr_get_info(GU, distribution=dist)
      CALL dbcsr_transposed(U, GU, shallow_data_copy=.FALSE., &
                            use_distribution=dist, &
                            transpose_distribution=.FALSE.)
      CALL dbcsr_add(GU, U, alpha_scalar=-0.5_dp, beta_scalar=0.5_dp)
      !-------------------------------------------------------------------
      !
      norm_fro = dbcsr_frobenius_norm(GU)
      norm_gct = dbcsr_gershgorin_norm(GU)
      !write(*,*) 'qs_ot_localize: ||P-I||_f=',norm_fro,' ||P-I||_GCT=',norm_gct
      !
      !kscale = CEILING(LOG(MIN(norm_fro,norm_gct))/LOG(2.0_dp))
      !scale  = LOG(MIN(norm_fro,norm_gct))/LOG(2.0_dp)
      !write(*,*) 'qs_ot_localize: scale=',scale,' kscale=',kscale
      !
      ! rescale for steepest descent
      CALL dbcsr_scale(GU, -alpha)
      !
      ! compute unitary transform
      ! zeroth and first order
      expfactor = 1.0_dp
      CALL dbcsr_copy(U, GU)
      CALL dbcsr_scale(U, expfactor)
      CALL dbcsr_add_on_diag(U, 1.0_dp)
      ! other orders
      CALL dbcsr_copy(Gp1, GU)
      DO i = 2, taylor_order
         ! new power of G
         CALL dbcsr_multiply('N', 'N', 1.0_dp, GU, Gp1, 0.0_dp, Gp2)
         CALL dbcsr_copy(Gp1, Gp2)
         ! add to the taylor expansion so far
         expfactor = expfactor/REAL(i, KIND=dp)
         CALL dbcsr_add(U, Gp1, alpha_scalar=1.0_dp, beta_scalar=expfactor)
         norm_fro = dbcsr_frobenius_norm(Gp1)
         !write(*,*) 'Taylor expansion i=',i,' norm(X^i)/i!=',norm_fro*expfactor
         IF (norm_fro*expfactor < 1.0E-10_dp) EXIT
      END DO
      !
      ! rotate MOs
      CALL dbcsr_multiply('N', 'N', 1.0_dp, C_NEW, U, 0.0_dp, C)
      CALL dbcsr_copy(C_NEW, C)
      !
      ! rotate SC
      CALL dbcsr_multiply('N', 'N', 1.0_dp, SC, U, 0.0_dp, C)
      CALL dbcsr_copy(SC, C)
      !
      ! rotate D_i
      CALL dbcsr_multiply('N', 'N', 1.0_dp, D, U, 0.0_dp, C)
      CALL dbcsr_copy(D, C)
      !
      ! rotate G_i-1
      IF (ASSOCIATED(G_OLD)) THEN
         CALL dbcsr_multiply('N', 'N', 1.0_dp, G_OLD, U, 0.0_dp, C)
         CALL dbcsr_copy(G_OLD, C)
      END IF
      !
      CALL timestop(handle)
   END SUBROUTINE qs_ot_on_the_fly_localize

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param C_OLD ...
!> \param C_TMP ...
!> \param C_NEW ...
!> \param P ...
!> \param SC ...
!> \param update ...
! **************************************************************************************************
   SUBROUTINE qs_ot_ref_chol(qs_ot_env, C_OLD, C_TMP, C_NEW, P, SC, update)
      !
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type)                                   :: C_OLD, C_TMP, C_NEW, P, SC
      LOGICAL, INTENT(IN)                                :: update

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_ot_ref_chol'

      INTEGER                                            :: handle, k, n

      CALL timeset(routineN, handle)
      !
      CALL dbcsr_get_info(C_NEW, nfullrows_total=n, nfullcols_total=k)
      !
      ! P = U'*U
      CALL cp_dbcsr_cholesky_decompose(P, k, qs_ot_env%para_env, qs_ot_env%blacs_env)
      !
      ! C_NEW = C_OLD*inv(U)
      CALL cp_dbcsr_cholesky_restore(C_OLD, k, P, C_NEW, op="SOLVE", pos="RIGHT", &
                                     transa="N", para_env=qs_ot_env%para_env, blacs_env=qs_ot_env%blacs_env)
      !
      ! Update SC if needed
      IF (update) THEN
         CALL cp_dbcsr_cholesky_restore(SC, k, P, C_TMP, op="SOLVE", pos="RIGHT", &
                                        transa="N", para_env=qs_ot_env%para_env, blacs_env=qs_ot_env%blacs_env)
         CALL dbcsr_copy(SC, C_TMP)
      END IF
      !
      CALL timestop(handle)
   END SUBROUTINE qs_ot_ref_chol

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param C_OLD ...
!> \param C_TMP ...
!> \param C_NEW ...
!> \param P ...
!> \param SC ...
!> \param update ...
! **************************************************************************************************
   SUBROUTINE qs_ot_ref_lwdn(qs_ot_env, C_OLD, C_TMP, C_NEW, P, SC, update)
      !
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type)                                   :: C_OLD, C_TMP, C_NEW, P, SC
      LOGICAL, INTENT(IN)                                :: update

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_ot_ref_lwdn'

      INTEGER                                            :: handle, i, k, n
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: eig, fun
      TYPE(dbcsr_type), POINTER                          :: V, W

      CALL timeset(routineN, handle)
      !
      CALL dbcsr_get_info(C_NEW, nfullrows_total=n, nfullcols_total=k)
      !
      V => qs_ot_env%buf1_k_k_nosym ! a buffer
      W => qs_ot_env%buf2_k_k_nosym ! a buffer
      ALLOCATE (eig(k), fun(k))
      !
      CALL cp_dbcsr_syevd(P, V, eig, qs_ot_env%para_env, qs_ot_env%blacs_env)
      !
      ! compute the P^(-1/2)
      DO i = 1, k
         IF (eig(i) <= 0.0_dp) THEN
            CPABORT("P not positive definite")
         END IF
         IF (eig(i) < 1.0E-8_dp) THEN
            fun(i) = 0.0_dp
         ELSE
            fun(i) = 1.0_dp/SQRT(eig(i))
         END IF
      END DO
      CALL dbcsr_copy(W, V)
      CALL dbcsr_scale_by_vector(V, alpha=fun, side='right')
      CALL dbcsr_multiply('N', 'T', 1.0_dp, W, V, 0.0_dp, P)
      !
      ! Update C
      CALL dbcsr_multiply('N', 'N', 1.0_dp, C_OLD, P, 0.0_dp, C_NEW)
      !
      ! Update SC if needed
      IF (update) THEN
         CALL dbcsr_multiply('N', 'N', 1.0_dp, SC, P, 0.0_dp, C_TMP)
         CALL dbcsr_copy(SC, C_TMP)
      END IF
      !
      DEALLOCATE (eig, fun)
      !
      CALL timestop(handle)
   END SUBROUTINE qs_ot_ref_lwdn

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param C_OLD ...
!> \param C_TMP ...
!> \param C_NEW ...
!> \param P ...
!> \param SC ...
!> \param norm_in ...
!> \param update ...
! **************************************************************************************************
   SUBROUTINE qs_ot_ref_poly(qs_ot_env, C_OLD, C_TMP, C_NEW, P, SC, norm_in, update)
      !
      TYPE(qs_ot_type)                                   :: qs_ot_env
      TYPE(dbcsr_type), POINTER                          :: C_OLD, C_TMP, C_NEW, P
      TYPE(dbcsr_type)                                   :: SC
      REAL(dp), INTENT(IN)                               :: norm_in
      LOGICAL, INTENT(IN)                                :: update

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_ot_ref_poly'

      INTEGER                                            :: handle, irefine, k, n
      LOGICAL                                            :: quick_exit
      REAL(dp)                                           :: norm, norm_fro, norm_gct, occ_in, &
                                                            occ_out, rescale
      TYPE(dbcsr_type), POINTER                          :: BUF1, BUF2, BUF_NOSYM, FT, FY

      CALL timeset(routineN, handle)
      !
      CALL dbcsr_get_info(C_NEW, nfullrows_total=n, nfullcols_total=k)
      !
      BUF_NOSYM => qs_ot_env%buf1_k_k_nosym ! a buffer
      BUF1 => qs_ot_env%buf1_k_k_sym ! a buffer
      BUF2 => qs_ot_env%buf2_k_k_sym ! a buffer
      FY => qs_ot_env%buf3_k_k_sym ! a buffer
      FT => qs_ot_env%buf4_k_k_sym ! a buffer
      !
      ! initialize the norm (already computed in qs_ot_get_orbitals_ref)
      norm = norm_in
      !
      ! can we do a quick exit?
      quick_exit = .FALSE.
      IF (norm < qs_ot_env%settings%eps_irac_quick_exit) quick_exit = .TRUE.
      !
      ! lets refine
      rescale = 1.0_dp
      DO irefine = 1, qs_ot_env%settings%max_irac
         !
         ! rescaling
         IF (norm > 1.0_dp) THEN
            CALL dbcsr_scale(P, 1.0_dp/norm)
            rescale = rescale/SQRT(norm)
         END IF
         !
         ! get the refinement polynomial
         CALL qs_ot_refine(P, FY, BUF1, BUF2, qs_ot_env%settings%irac_degree, &
                           qs_ot_env%settings%eps_irac_filter_matrix)
         !
         ! collect the transformation
         IF (irefine == 1) THEN
            CALL dbcsr_copy(FT, FY, name='FT')
         ELSE
            CALL dbcsr_multiply('N', 'N', 1.0_dp, FT, FY, 0.0_dp, BUF1)
            IF (qs_ot_env%settings%eps_irac_filter_matrix > 0.0_dp) THEN
               occ_in = dbcsr_get_occupation(buf1)
               CALL dbcsr_filter(buf1, qs_ot_env%settings%eps_irac_filter_matrix)
               occ_out = dbcsr_get_occupation(buf1)
            END IF
            CALL dbcsr_copy(FT, BUF1, name='FT')
         END IF
         !
         ! quick exit if possible
         IF (quick_exit) THEN
            EXIT
         END IF
         !
         ! P = FY^T * P * FY
         CALL dbcsr_multiply('N', 'N', 1.0_dp, P, FY, 0.0_dp, BUF_NOSYM)
         IF (qs_ot_env%settings%eps_irac_filter_matrix > 0.0_dp) THEN
            occ_in = dbcsr_get_occupation(buf_nosym)
            CALL dbcsr_filter(buf_nosym, qs_ot_env%settings%eps_irac_filter_matrix)
            occ_out = dbcsr_get_occupation(buf_nosym)
         END IF
         CALL dbcsr_multiply('N', 'N', 1.0_dp, FY, BUF_NOSYM, 0.0_dp, P)
         IF (qs_ot_env%settings%eps_irac_filter_matrix > 0.0_dp) THEN
            occ_in = dbcsr_get_occupation(p)
            CALL dbcsr_filter(p, qs_ot_env%settings%eps_irac_filter_matrix)
            occ_out = dbcsr_get_occupation(p)
         END IF
         !
         ! check ||P-1||_gct
         CALL dbcsr_add_on_diag(P, -1.0_dp)
         norm_fro = dbcsr_frobenius_norm(P)
         norm_gct = dbcsr_gershgorin_norm(P)
         CALL dbcsr_add_on_diag(P, 1.0_dp)
         norm = MIN(norm_gct, norm_fro)
         !
         ! printing
         !
         ! blows up
         IF (norm > 1.0E10_dp) THEN
            CALL cp_abort(__LOCATION__, &
                          "Refinement blows up! "// &
                          "We need you to improve the code, please post your input on "// &
                          "the forum https://www.cp2k.org/")
         END IF
         !
         ! can we do a quick exit next step?
         IF (norm < qs_ot_env%settings%eps_irac_quick_exit) quick_exit = .TRUE.
         !
         ! are we done?
         IF (norm < qs_ot_env%settings%eps_irac) EXIT
         !
      END DO
      !
      ! C_NEW = C_NEW * FT * rescale
      CALL dbcsr_multiply('N', 'N', rescale, C_OLD, FT, 0.0_dp, C_NEW)
      IF (qs_ot_env%settings%eps_irac_filter_matrix > 0.0_dp) THEN
         occ_in = dbcsr_get_occupation(c_new)
         CALL dbcsr_filter(c_new, qs_ot_env%settings%eps_irac_filter_matrix)
         occ_out = dbcsr_get_occupation(c_new)
      END IF
      !
      ! update SC = SC * FY * rescale
      IF (update) THEN
         CALL dbcsr_multiply('N', 'N', rescale, SC, FT, 0.0_dp, C_TMP)
         IF (qs_ot_env%settings%eps_irac_filter_matrix > 0.0_dp) THEN
            occ_in = dbcsr_get_occupation(c_tmp)
            CALL dbcsr_filter(c_tmp, qs_ot_env%settings%eps_irac_filter_matrix)
            occ_out = dbcsr_get_occupation(c_tmp)
         END IF
         CALL dbcsr_copy(SC, C_TMP)
      END IF
      !
      CALL timestop(handle)
   END SUBROUTINE qs_ot_ref_poly

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env1 ...
!> \return ...
! **************************************************************************************************
   FUNCTION qs_ot_ref_update(qs_ot_env1) RESULT(update)
      !
      TYPE(qs_ot_type)                                   :: qs_ot_env1
      LOGICAL                                            :: update

      update = .FALSE.
      SELECT CASE (qs_ot_env1%settings%ot_method)
      CASE ("CG")
         SELECT CASE (qs_ot_env1%settings%line_search_method)
         CASE ("2PNT")
            IF (qs_ot_env1%line_search_count == 2) update = .TRUE.
         CASE DEFAULT
            CPABORT("NYI")
         END SELECT
      CASE ("DIIS")
         update = .TRUE.
      CASE ("BROY", "LBFG")
         ! These minimizers retain positions or secants in one fixed REF chart.
         update = .FALSE.
      CASE DEFAULT
         CPABORT("NYI")
      END SELECT
   END FUNCTION qs_ot_ref_update

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env1 ...
!> \param norm_in ...
!> \param ortho_irac ...
! **************************************************************************************************
   SUBROUTINE qs_ot_ref_decide(qs_ot_env1, norm_in, ortho_irac)
      !
      TYPE(qs_ot_type)                                   :: qs_ot_env1
      REAL(dp), INTENT(IN)                               :: norm_in
      CHARACTER(LEN=*), INTENT(INOUT)                    :: ortho_irac

      ortho_irac = qs_ot_env1%settings%ortho_irac
      IF (norm_in < qs_ot_env1%settings%eps_irac_switch) ortho_irac = "POLY"
   END SUBROUTINE qs_ot_ref_decide

! **************************************************************************************************
!> \brief ...
!> \param matrix_c ...
!> \param matrix_s ...
!> \param matrix_x ...
!> \param matrix_sx ...
!> \param matrix_gx_old ...
!> \param matrix_dx ...
!> \param qs_ot_env ...
!> \param qs_ot_env1 ...
! **************************************************************************************************
   SUBROUTINE qs_ot_get_orbitals_ref(matrix_c, matrix_s, matrix_x, matrix_sx, &
                                     matrix_gx_old, matrix_dx, qs_ot_env, qs_ot_env1)
      !
      TYPE(dbcsr_type), POINTER                          :: matrix_c, matrix_s, matrix_x, matrix_sx, &
                                                            matrix_gx_old, matrix_dx
      TYPE(qs_ot_type)                                   :: qs_ot_env, qs_ot_env1

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_orbitals_ref'

      CHARACTER(LEN=4)                                   :: ortho_irac
      INTEGER                                            :: handle, k, n
      LOGICAL                                            :: on_the_fly_loc, update
      REAL(dp)                                           :: norm, norm_fro, norm_gct, occ_in, occ_out
      TYPE(dbcsr_type), POINTER                          :: C_NEW, C_OLD, C_TMP, D, G_OLD, P, S, SC

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix_c, nfullrows_total=n, nfullcols_total=k)
      !
      C_NEW => matrix_c
      C_OLD => matrix_x ! need to be carefully updated for the gradient !
      SC => matrix_sx ! need to be carefully updated for the gradient !
      G_OLD => matrix_gx_old ! need to be carefully updated for localization !
      D => matrix_dx ! need to be carefully updated for localization !
      S => matrix_s

      P => qs_ot_env%p_k_k_sym ! a buffer
      C_TMP => qs_ot_env%buf1_n_k ! a buffer
      !
      ! do we need to update C_OLD and SC?
      update = qs_ot_ref_update(qs_ot_env1)
      !
      ! do we want to on the fly localize?
      ! for the moment this is set from the input,
      ! later we might want to localize every n-step or
      ! when the sparsity increases...
      on_the_fly_loc = qs_ot_env1%settings%on_the_fly_loc
      !
      ! compute SC = S*C
      IF (ASSOCIATED(S)) THEN
         CALL dbcsr_multiply('N', 'N', 1.0_dp, S, C_OLD, 0.0_dp, SC)
         IF (qs_ot_env1%settings%eps_irac_filter_matrix > 0.0_dp) THEN
            occ_in = dbcsr_get_occupation(sc)
            CALL dbcsr_filter(sc, qs_ot_env1%settings%eps_irac_filter_matrix)
            occ_out = dbcsr_get_occupation(sc)
         END IF
      ELSE
         CALL dbcsr_copy(SC, C_OLD)
      END IF
      !
      ! compute P = C'*SC
      CALL dbcsr_multiply('T', 'N', 1.0_dp, C_OLD, SC, 0.0_dp, P)
      IF (qs_ot_env1%settings%eps_irac_filter_matrix > 0.0_dp) THEN
         occ_in = dbcsr_get_occupation(p)
         CALL dbcsr_filter(p, qs_ot_env1%settings%eps_irac_filter_matrix)
         occ_out = dbcsr_get_occupation(p)
      END IF
      !
      ! check ||P-1||_f and ||P-1||_gct
      CALL dbcsr_add_on_diag(P, -1.0_dp)
      norm_fro = dbcsr_frobenius_norm(P)
      norm_gct = dbcsr_gershgorin_norm(P)
      CALL dbcsr_add_on_diag(P, 1.0_dp)
      norm = MIN(norm_gct, norm_fro)
      CALL qs_ot_ref_decide(qs_ot_env1, norm, ortho_irac)
      !
      ! select the orthogonality method
      SELECT CASE (ortho_irac)
      CASE ("CHOL")
         CALL qs_ot_ref_chol(qs_ot_env, C_OLD, C_TMP, C_NEW, P, SC, update)
      CASE ("LWDN")
         CALL qs_ot_ref_lwdn(qs_ot_env, C_OLD, C_TMP, C_NEW, P, SC, update)
      CASE ("POLY")
         CALL qs_ot_ref_poly(qs_ot_env, C_OLD, C_TMP, C_NEW, P, SC, norm, update)
      CASE DEFAULT
         CPABORT("Wrong argument")
      END SELECT
      !
      ! We update the C_i+1 and localization
      IF (update) THEN
         IF (on_the_fly_loc) THEN
            CALL qs_ot_on_the_fly_localize(qs_ot_env, C_NEW, SC, G_OLD, D)
         END IF
         CALL dbcsr_copy(C_OLD, C_NEW)
      END IF

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL qs_ot_generate_rotation(qs_ot_env)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, C_NEW, qs_ot_env%rot_mat_u, &
                             0.0_dp, C_TMP)
         CALL dbcsr_copy(C_NEW, C_TMP)
      END IF
      !
      CALL timestop(handle)
   END SUBROUTINE qs_ot_get_orbitals_ref

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
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s, c_re, 0.0_dp, sc_re)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s_im, c_im, 0.0_dp, tmp_nk)
      CALL dbcsr_add(sc_re, tmp_nk, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s, c_im, 0.0_dp, sc_im)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s_im, c_re, 0.0_dp, tmp_nk)
      CALL dbcsr_add(sc_im, tmp_nk, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

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
      CALL dbcsr_multiply('N', 'N', 1.0_dp, c_re, p_re, 0.0_dp, matrix_c)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, c_im, p_im, 0.0_dp, tmp_nk)
      CALL dbcsr_add(matrix_c, tmp_nk, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      CALL dbcsr_multiply('N', 'N', 1.0_dp, c_re, p_im, 0.0_dp, matrix_c_im)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, c_im, p_re, 0.0_dp, tmp_nk)
      CALL dbcsr_add(matrix_c_im, tmp_nk, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      IF (update) THEN
         CALL dbcsr_copy(c_re, matrix_c)
         CALL dbcsr_copy(c_im, matrix_c_im)

         CALL dbcsr_set(f_re, 0.0_dp)
         CALL dbcsr_add_on_diag(f_re, alpha=1.0_dp)
         CALL dbcsr_set(f_im, 0.0_dp)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s, c_re, 0.0_dp, sc_re)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s_im, c_im, 0.0_dp, tmp_nk)
         CALL dbcsr_add(sc_re, tmp_nk, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s, c_im, 0.0_dp, sc_im)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s_im, c_re, 0.0_dp, tmp_nk)
         CALL dbcsr_add(sc_im, tmp_nk, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL qs_ot_generate_rotation_complex(qs_ot_env)

         CALL dbcsr_copy(rotated_re, matrix_c, name="rotated_re")
         CALL dbcsr_copy(rotated_im, matrix_c_im, name="rotated_im")
         CALL dbcsr_copy(rotation_tmp, matrix_c, name="rotation_tmp")

         ! C_out = Q*U for complex Q and unitary U.
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_c, qs_ot_env%rot_mat_u, &
                             0.0_dp, rotated_re)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_c_im, qs_ot_env%rot_mat_u_im, &
                             0.0_dp, rotation_tmp)
         CALL dbcsr_add(rotated_re, rotation_tmp, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_c, qs_ot_env%rot_mat_u_im, &
                             0.0_dp, rotated_im)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_c_im, qs_ot_env%rot_mat_u, &
                             0.0_dp, rotation_tmp)
         CALL dbcsr_add(rotated_im, rotation_tmp, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

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
!> \brief  refinement polynomial of degree 2,3 and 4 (PRB 70, 193102 (2004))
!> \param P ...
!> \param FY ...
!> \param P2 ...
!> \param T ...
!> \param irac_degree ...
!> \param eps_irac_filter_matrix ...
! **************************************************************************************************
   SUBROUTINE qs_ot_refine(P, FY, P2, T, irac_degree, eps_irac_filter_matrix)
      TYPE(dbcsr_type), INTENT(inout)                    :: P, FY, P2, T
      INTEGER, INTENT(in)                                :: irac_degree
      REAL(dp), INTENT(in)                               :: eps_irac_filter_matrix

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_ot_refine'

      INTEGER                                            :: handle, k
      REAL(dp)                                           :: occ_in, occ_out, r

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(P, nfullcols_total=k)
      SELECT CASE (irac_degree)
      CASE (2)
         ! C_out = C_in * ( 15/8 * I - 10/8 * P + 3/8 * P^2)
         r = 3.0_dp/8.0_dp
         CALL dbcsr_multiply('N', 'N', r, P, P, 0.0_dp, FY)
         IF (eps_irac_filter_matrix > 0.0_dp) THEN
            occ_in = dbcsr_get_occupation(fy)
            CALL dbcsr_filter(fy, eps_irac_filter_matrix)
            occ_out = dbcsr_get_occupation(fy)
         END IF
         r = -10.0_dp/8.0_dp
         CALL dbcsr_add(FY, P, alpha_scalar=1.0_dp, beta_scalar=r)
         r = 15.0_dp/8.0_dp
         CALL dbcsr_add_on_diag(FY, alpha=r)
      CASE (3)
         ! C_out = C_in * ( 35/16 * I - 35/16 * P + 21/16 * P^2 - 5/16 P^3)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, P, P, 0.0_dp, P2)
         IF (eps_irac_filter_matrix > 0.0_dp) THEN
            occ_in = dbcsr_get_occupation(p2)
            CALL dbcsr_filter(p2, eps_irac_filter_matrix)
            occ_out = dbcsr_get_occupation(p2)
         END IF
         r = -5.0_dp/16.0_dp
         CALL dbcsr_multiply('N', 'N', r, P2, P, 0.0_dp, FY)
         IF (eps_irac_filter_matrix > 0.0_dp) THEN
            occ_in = dbcsr_get_occupation(fy)
            CALL dbcsr_filter(fy, eps_irac_filter_matrix)
            occ_out = dbcsr_get_occupation(fy)
         END IF
         r = 21.0_dp/16.0_dp
         CALL dbcsr_add(FY, P2, alpha_scalar=1.0_dp, beta_scalar=r)
         r = -35.0_dp/16.0_dp
         CALL dbcsr_add(FY, P, alpha_scalar=1.0_dp, beta_scalar=r)
         r = 35.0_dp/16.0_dp
         CALL dbcsr_add_on_diag(FY, alpha=r)
      CASE (4)
         ! C_out = C_in * ( 315/128 * I - 420/128 * P + 378/128 * P^2 - 180/128 P^3 + 35/128 P^4 )
         !       = C_in * ( 315/128 * I - 420/128 * P + 378/128 * P^2 + ( - 180/128 * P + 35/128 * P^2 ) * P^2 )
         CALL dbcsr_multiply('N', 'N', 1.0_dp, P, P, 0.0_dp, P2) ! P^2
         IF (eps_irac_filter_matrix > 0.0_dp) THEN
            occ_in = dbcsr_get_occupation(p2)
            CALL dbcsr_filter(p2, eps_irac_filter_matrix)
            occ_out = dbcsr_get_occupation(p2)
         END IF
         r = -180.0_dp/128.0_dp
         CALL dbcsr_add(T, P, alpha_scalar=0.0_dp, beta_scalar=r) ! T=-180/128*P
         r = 35.0_dp/128.0_dp
         CALL dbcsr_add(T, P2, alpha_scalar=1.0_dp, beta_scalar=r) ! T=T+35/128*P^2
         CALL dbcsr_multiply('N', 'N', 1.0_dp, T, P2, 0.0_dp, FY) ! Y=T*P^2
         IF (eps_irac_filter_matrix > 0.0_dp) THEN
            occ_in = dbcsr_get_occupation(fy)
            CALL dbcsr_filter(fy, eps_irac_filter_matrix)
            occ_out = dbcsr_get_occupation(fy)
         END IF
         r = 378.0_dp/128.0_dp
         CALL dbcsr_add(FY, P2, alpha_scalar=1.0_dp, beta_scalar=r) ! Y=Y+378/128*P^2
         r = -420.0_dp/128.0_dp
         CALL dbcsr_add(FY, P, alpha_scalar=1.0_dp, beta_scalar=r) ! Y=Y-420/128*P
         r = 315.0_dp/128.0_dp
         CALL dbcsr_add_on_diag(FY, alpha=r) ! Y=Y+315/128*I
      CASE DEFAULT
         CPABORT("This irac_order NYI")
      END SELECT
      CALL timestop(handle)
   END SUBROUTINE qs_ot_refine

! **************************************************************************************************
!> \brief ...
!> \param matrix_hc ...
!> \param matrix_x ...
!> \param matrix_sx ...
!> \param matrix_gx ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE qs_ot_get_derivative_ref(matrix_hc, matrix_x, matrix_sx, matrix_gx, &
                                       qs_ot_env)
      TYPE(dbcsr_type), POINTER                          :: matrix_hc, matrix_x, matrix_sx, matrix_gx
      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative_ref'

      INTEGER                                            :: handle, k, n
      REAL(dp)                                           :: occ_in, occ_out
      TYPE(dbcsr_type), POINTER                          :: C, CHC, G, HC, HC_work, SC

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix_x, nfullrows_total=n, nfullcols_total=k)
      !
      C => matrix_x ! NBsf*NOcc
      SC => matrix_sx ! NBsf*NOcc need to be up2date
      HC => matrix_hc ! NBsf*NOcc
      G => matrix_gx ! NBsf*NOcc
      CHC => qs_ot_env%buf1_k_k_sym ! buffer

      IF (qs_ot_env%settings%do_rotation) THEN
         ! The physical orbitals are C_current=Q*U. Pull dE/dC_current
         ! back to the unrotated REF basis before projecting it.
         CALL qs_ot_rot_mat_derivative(qs_ot_env)
         HC_work => qs_ot_env%buf1_n_k
         CALL dbcsr_multiply('N', 'T', 1.0_dp, HC, qs_ot_env%rot_mat_u, &
                             0.0_dp, HC_work)
      ELSE
         HC_work => HC
      END IF

      ! C'*(H*C)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, C, HC_work, 0.0_dp, CHC)
      IF (qs_ot_env%settings%eps_irac_filter_matrix > 0.0_dp) THEN
         occ_in = dbcsr_get_occupation(chc)
         CALL dbcsr_filter(chc, qs_ot_env%settings%eps_irac_filter_matrix)
         occ_out = dbcsr_get_occupation(chc)
      END IF
      ! (S*C)*(C'*H*C)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, SC, CHC, 0.0_dp, G)
      IF (qs_ot_env%settings%eps_irac_filter_matrix > 0.0_dp) THEN
         occ_in = dbcsr_get_occupation(g)
         CALL dbcsr_filter(g, qs_ot_env%settings%eps_irac_filter_matrix)
         occ_out = dbcsr_get_occupation(g)
      END IF
      ! G = 2*(1-S*C*C')*H*C
      CALL dbcsr_add(G, HC_work, alpha_scalar=-1.0_dp, beta_scalar=1.0_dp)
      !
      CALL timestop(handle)
   END SUBROUTINE qs_ot_get_derivative_ref

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
      REAL(dp)                                           :: occ_in, occ_out
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
      CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%matrix_x, f_re, 0.0_dp, q_re)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%matrix_x_im, f_im, 0.0_dp, g_re)
      CALL dbcsr_add(q_re, g_re, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%matrix_x, f_im, 0.0_dp, q_im)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%matrix_x_im, f_re, 0.0_dp, g_re)
      CALL dbcsr_add(q_im, g_re, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

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
         occ_in = dbcsr_get_occupation(b_re)
         CALL dbcsr_filter(b_re, qs_ot_env%settings%eps_irac_filter_matrix)
         occ_out = dbcsr_get_occupation(b_re)
         occ_in = dbcsr_get_occupation(b_im)
         CALL dbcsr_filter(b_im, qs_ot_env%settings%eps_irac_filter_matrix)
         occ_out = dbcsr_get_occupation(b_im)
      END IF

      ! S*Q = (S*X)*F.  G is temporary storage for this pair.
      CALL dbcsr_multiply('N', 'N', 1.0_dp, sc_re, f_re, 0.0_dp, g_re)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, sc_im, f_im, 0.0_dp, q_re)
      CALL dbcsr_add(g_re, q_re, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      CALL dbcsr_multiply('N', 'N', 1.0_dp, sc_re, f_im, 0.0_dp, g_im)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, sc_im, f_re, 0.0_dp, q_re)
      CALL dbcsr_add(g_im, q_re, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

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
         occ_in = dbcsr_get_occupation(g_re)
         CALL dbcsr_filter(g_re, qs_ot_env%settings%eps_irac_filter_matrix)
         occ_out = dbcsr_get_occupation(g_re)
         occ_in = dbcsr_get_occupation(g_im)
         CALL dbcsr_filter(g_im, qs_ot_env%settings%eps_irac_filter_matrix)
         occ_out = dbcsr_get_occupation(g_im)
      END IF

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL dbcsr_release(hc_rot_re)
         CALL dbcsr_release(hc_rot_im)
         CALL dbcsr_release(tmp_nk)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_ot_get_derivative_ref_complex

! **************************************************************************************************
!> \brief Transpose a square DBCSR matrix while retaining its target distribution.
!> \param matrix source matrix
!> \param transposed transposed matrix
!> \param identity_template square matrix template
! **************************************************************************************************
   SUBROUTINE qs_ot_square_transpose(matrix, transposed, identity_template)
      TYPE(dbcsr_type)                                   :: matrix, transposed, identity_template

      TYPE(dbcsr_type)                                   :: identity

      CALL dbcsr_copy(transposed, matrix, name='square_transposed')
      CALL dbcsr_copy(identity, identity_template, name='transpose_identity')
      CALL dbcsr_set(identity, 0.0_dp)
      CALL dbcsr_add_on_diag(identity, 1.0_dp)
      CALL dbcsr_multiply('T', 'N', 1.0_dp, matrix, identity, 0.0_dp, transposed)
      CALL dbcsr_release(identity)

   END SUBROUTINE qs_ot_square_transpose

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
      CALL dbcsr_multiply('N', 'N', 1.0_dp, anti_re, v_re, 0.0_dp, work_re)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, anti_im, v_im, 0.0_dp, tmp_kk)
      CALL dbcsr_add(work_re, tmp_kk, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, anti_re, v_im, 0.0_dp, work_im)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, anti_im, v_re, 0.0_dp, tmp_kk)
      CALL dbcsr_add(work_im, tmp_kk, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

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
      CALL dbcsr_multiply('N', 'N', 1.0_dp, v_re, inner_re, 0.0_dp, work_re)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, v_im, inner_im, 0.0_dp, tmp_kk)
      CALL dbcsr_add(work_re, tmp_kk, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, v_re, inner_im, 0.0_dp, work_im)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, v_im, inner_re, 0.0_dp, tmp_kk)
      CALL dbcsr_add(work_im, tmp_kk, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

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
      CALL dbcsr_multiply('N', 'N', 1.0_dp, sx_re, f_re, 0.0_dp, sq_re)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, sx_im, f_im, 0.0_dp, tmp_nk)
      CALL dbcsr_add(sq_re, tmp_nk, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, sx_re, f_im, 0.0_dp, sq_im)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, sx_im, f_re, 0.0_dp, tmp_nk)
      CALL dbcsr_add(sq_im, tmp_nk, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

      CALL dbcsr_copy(vertical_re, sx_re, name='vertical_re')
      CALL dbcsr_copy(vertical_im, sx_im, name='vertical_im')
      CALL dbcsr_multiply('N', 'N', 1.0_dp, sq_re, z_re, 0.0_dp, vertical_re)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, sq_im, z_im, 0.0_dp, tmp_nk)
      CALL dbcsr_add(vertical_re, tmp_nk, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, sq_re, z_im, 0.0_dp, vertical_im)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, sq_im, z_re, 0.0_dp, tmp_nk)
      CALL dbcsr_add(vertical_im, tmp_nk, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
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
!> \brief computes p=x*S*x and the matrix functionals related matrices
!> \param matrix_x ...
!> \param matrix_sx ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE qs_ot_get_p(matrix_x, matrix_sx, qs_ot_env)

      TYPE(dbcsr_type), POINTER                          :: matrix_x, matrix_sx
      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_ot_get_p'
      REAL(KIND=dp), PARAMETER                           :: rone = 1.0_dp, rzero = 0.0_dp

      INTEGER                                            :: handle, k, max_iter, n
      LOGICAL                                            :: converged
      REAL(KIND=dp)                                      :: max_ev, min_ev, threshold

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix_x, nfullrows_total=n, nfullcols_total=k)

      ! get the overlap
      CALL dbcsr_multiply('T', 'N', rone, matrix_x, matrix_sx, rzero, &
                          qs_ot_env%matrix_p)

      ! get an upper bound for the largest eigenvalue
      ! try using lancos first and fall back to gershgorin norm if it fails
      max_iter = 30; threshold = 1.0E-03_dp
      CALL arnoldi_extremal(qs_ot_env%matrix_p, max_ev, min_ev, converged, threshold, max_iter)
      qs_ot_env%largest_eval_upper_bound = MAX(max_ev, ABS(min_ev))

      IF (.NOT. converged) qs_ot_env%largest_eval_upper_bound = dbcsr_gershgorin_norm(qs_ot_env%matrix_p)
      CALL decide_strategy(qs_ot_env)
      IF (qs_ot_env%do_taylor) THEN
         CALL qs_ot_p2m_taylor(qs_ot_env)
      ELSE
         CALL qs_ot_p2m_diag(qs_ot_env)
      END IF

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL qs_ot_generate_rotation(qs_ot_env)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_ot_get_p

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
         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_dedu, &
                             qs_ot_env%rot_mat_evec_re, 0.0_dp, work_re)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_dedu_im, &
                             qs_ot_env%rot_mat_evec_im, 0.0_dp, tmp)
         CALL dbcsr_add(work_re, tmp, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_dedu, &
                             qs_ot_env%rot_mat_evec_im, 0.0_dp, work_im)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_dedu_im, &
                             qs_ot_env%rot_mat_evec_re, 0.0_dp, tmp)
         CALL dbcsr_add(work_im, tmp, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
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

         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_evec_re, &
                             outer_deriv_re, 0.0_dp, work_re)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_evec_im, &
                             outer_deriv_im, 0.0_dp, tmp)
         CALL dbcsr_add(work_re, tmp, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_evec_re, &
                             outer_deriv_im, 0.0_dp, work_im)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, qs_ot_env%rot_mat_evec_im, &
                             outer_deriv_re, 0.0_dp, tmp)
         CALL dbcsr_add(work_im, tmp, alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

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
! **************************************************************************************************
   SUBROUTINE qs_ot_get_p_complex(matrix_x, matrix_x_im, matrix_sx, matrix_sx_im, qs_ot_env)
      TYPE(dbcsr_type), POINTER                          :: matrix_x, matrix_x_im, matrix_sx, &
                                                            matrix_sx_im
      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_p_complex'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(qs_ot_env%has_complex_kpoint_state)
      CPASSERT(qs_ot_env%settings%ot_algorithm == 'TOD')

      CALL qs_ot_complex_multiply('C', 'N', matrix_x, matrix_x_im, matrix_sx, matrix_sx_im, &
                                  qs_ot_env%matrix_p, qs_ot_env%matrix_p_im, qs_ot_env%matrix_buf1)
      qs_ot_env%do_taylor = .FALSE.
      CALL qs_ot_p2m_diag_complex(qs_ot_env)

      CALL timestop(handle)

   END SUBROUTINE qs_ot_get_p_complex

! **************************************************************************************************
!> \brief computes the rotation matrix rot_mat_u that is associated to a given
!>        rot_mat_x using rot_mat_u=exp(rot_mat_x)
!> \param qs_ot_env a valid qs_ot_env
!> \par History
!>      08.2004 created [Joost VandeVondele]
!>      12.2024 Rewrite to use only real matrices [Ole Schuett]
! **************************************************************************************************
   SUBROUTINE qs_ot_generate_rotation(qs_ot_env)

      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_generate_rotation'

      INTEGER                                            :: handle, k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: exp_evals_im, exp_evals_re
      TYPE(dbcsr_type)                                   :: buf_1, buf_2

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(qs_ot_env%rot_mat_x, nfullrows_total=k)

      IF (k /= 0) THEN
         ! We want to compute: rot_mat_u = exp(i*rot_mat_x)

         ! Diagonalize: matrix = i*rot_mat_x.
         ! Note that matrix is imaginary and hermitian because rot_mat_x is real and anti-symmetric.
         CALL cp_dbcsr_heevd(matrix_im=qs_ot_env%rot_mat_x, &  ! matrix_re omitted because it's zero
                             eigenvectors_re=qs_ot_env%rot_mat_evec_re, &
                             eigenvectors_im=qs_ot_env%rot_mat_evec_im, &
                             eigenvalues=qs_ot_env%rot_mat_evals, &
                             para_env=qs_ot_env%para_env, &
                             blacs_env=qs_ot_env%blacs_env)

         ! Compute: exp_evals = EXP(-i*rot_mat_evals)
         ALLOCATE (exp_evals_re(k), exp_evals_im(k))
         exp_evals_re(:) = COS(-qs_ot_env%rot_mat_evals(:))
         exp_evals_im(:) = SIN(-qs_ot_env%rot_mat_evals(:))

         ! Compute: rot_mat_u = \sum_ij exp_evals_ij * |rot_mat_evec_i> <rot_mat_evec_j|
         ! Note that we need only two matrix multiplications because rot_mat_u is real.
         CALL dbcsr_copy(buf_1, qs_ot_env%rot_mat_evec_re, name="buf_1")
         CALL dbcsr_scale_by_vector(buf_1, alpha=exp_evals_re, side='right')
         CALL dbcsr_copy(buf_2, qs_ot_env%rot_mat_evec_im, name="buf_2")
         CALL dbcsr_scale_by_vector(buf_2, alpha=exp_evals_im, side='right')
         CALL dbcsr_add(buf_1, buf_2, alpha_scalar=+1.0_dp, beta_scalar=-1.0_dp)
         CALL dbcsr_multiply('N', 'T', 1.0_dp, buf_1, qs_ot_env%rot_mat_evec_re, 0.0_dp, qs_ot_env%rot_mat_u)

         CALL dbcsr_copy(buf_1, qs_ot_env%rot_mat_evec_im)
         CALL dbcsr_scale_by_vector(buf_1, alpha=exp_evals_re, side='right')
         CALL dbcsr_copy(buf_2, qs_ot_env%rot_mat_evec_re)
         CALL dbcsr_scale_by_vector(buf_2, alpha=exp_evals_im, side='right')
         CALL dbcsr_add(buf_1, buf_2, alpha_scalar=+1.0_dp, beta_scalar=+1.0_dp)
         CALL dbcsr_multiply('N', 'T', 1.0_dp, buf_1, qs_ot_env%rot_mat_evec_im, 1.0_dp, qs_ot_env%rot_mat_u)

         ! Clean up.
         CALL dbcsr_release(buf_1)
         CALL dbcsr_release(buf_2)
         DEALLOCATE (exp_evals_re, exp_evals_im)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_ot_generate_rotation

! **************************************************************************************************
!> \brief computes the derivative fields with respect to rot_mat_x
!> \param qs_ot_env valid qs_ot_env. In particular qs_ot_generate_rotation has to be called before
!>                        and the rot_mat_dedu matrix has to be up to date
!> \par History
!>      08.2004 created [ Joost VandeVondele ]
!>      12.2024 Rewrite to use only real matrices [Ole Schuett]
! **************************************************************************************************
   SUBROUTINE qs_ot_rot_mat_derivative(qs_ot_env)
      TYPE(qs_ot_type)                         :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_rot_mat_derivative'

      INTEGER                                  :: col, handle, i, iblock, j, k, max_blocks, nblocks, row
      INTEGER, ALLOCATABLE, DIMENSION(:)       :: cols, rows
      INTEGER, DIMENSION(:), POINTER           :: col_blk_offset, col_blk_size, row_blk_offset, row_blk_size
      REAL(KIND=dp)                            :: e1, e2
      TYPE(dbcsr_type)                         :: outer_deriv_re, outer_deriv_im, mat_buf, &
                                                  inner_deriv_re, inner_deriv_im
      TYPE(dbcsr_distribution_type)            :: dist
      TYPE(dbcsr_iterator_type)                :: iter
      REAL(dp), DIMENSION(:, :), POINTER       :: block_in_re, block_in_im, block_out_re, block_out_im
      LOGICAL                                  :: duplicate, found_in_im, found_in_re, found_out_im, &
                                                  found_out_re
      REAL(KIND=dp)                            :: im_part, re_part
      COMPLEX(dp)                              :: cval_in, cval_out
      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(qs_ot_env%rot_mat_u, nfullrows_total=k)
      IF (k /= 0) THEN
         CALL dbcsr_copy(qs_ot_env%matrix_buf1, qs_ot_env%rot_mat_dedu)
         ! now we get to the derivative wrt the antisymmetric matrix rot_mat_x
         CALL dbcsr_copy(mat_buf, qs_ot_env%rot_mat_dedu, "mat_buf")

         ! inner_deriv_ij = <rot_mat_evec_i| rot_mat_dedu |rot_mat_evec_j>
         CALL dbcsr_copy(inner_deriv_re, qs_ot_env%rot_mat_dedu, "inner_deriv_re") ! TODO just create
         CALL dbcsr_copy(inner_deriv_im, qs_ot_env%rot_mat_dedu, "inner_deriv_im") ! TODO just create

         CALL dbcsr_multiply('T', 'N', +1.0_dp, qs_ot_env%rot_mat_dedu, qs_ot_env%rot_mat_evec_im, 0.0_dp, mat_buf)
         CALL dbcsr_multiply('T', 'N', +1.0_dp, qs_ot_env%rot_mat_evec_im, mat_buf, 0.0_dp, inner_deriv_re)
         CALL dbcsr_multiply('T', 'N', +1.0_dp, qs_ot_env%rot_mat_evec_re, mat_buf, 0.0_dp, inner_deriv_im)

         CALL dbcsr_multiply('T', 'N', +1.0_dp, qs_ot_env%rot_mat_dedu, qs_ot_env%rot_mat_evec_re, 0.0_dp, mat_buf)
         CALL dbcsr_multiply('T', 'N', +1.0_dp, qs_ot_env%rot_mat_evec_re, mat_buf, 1.0_dp, inner_deriv_re)
         CALL dbcsr_multiply('T', 'N', -1.0_dp, qs_ot_env%rot_mat_evec_im, mat_buf, 1.0_dp, inner_deriv_im)

         ! Real and imaginary products can have different sparse block patterns.
         ! Form their union explicitly and treat a missing partner block as zero.
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
            duplicate = .FALSE.
            DO iblock = 1, nblocks
               duplicate = duplicate .OR. (rows(iblock) == row .AND. cols(iblock) == col)
            END DO
            IF (.NOT. duplicate) THEN
               nblocks = nblocks + 1
               rows(nblocks) = row
               cols(nblocks) = col
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)
         CALL dbcsr_iterator_start(iter, inner_deriv_im)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row, col)
            duplicate = .FALSE.
            DO iblock = 1, nblocks
               duplicate = duplicate .OR. (rows(iblock) == row .AND. cols(iblock) == col)
            END DO
            IF (.NOT. duplicate) THEN
               nblocks = nblocks + 1
               rows(nblocks) = row
               cols(nblocks) = col
            END IF
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

         CALL dbcsr_get_info(outer_deriv_re, row_blk_offset=row_blk_offset, col_blk_offset=col_blk_offset)
         CALL dbcsr_iterator_start(iter, outer_deriv_re)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row, col)
            CALL dbcsr_get_block_p(inner_deriv_re, row, col, block_in_re, found_in_re)
            CALL dbcsr_get_block_p(inner_deriv_im, row, col, block_in_im, found_in_im)
            CALL dbcsr_get_block_p(outer_deriv_re, row, col, block_out_re, found_out_re)
            CALL dbcsr_get_block_p(outer_deriv_im, row, col, block_out_im, found_out_im)
            CPASSERT(found_out_re .AND. found_out_im)

            DO i = 1, SIZE(block_out_re, 1)
            DO j = 1, SIZE(block_out_re, 2)
               e1 = qs_ot_env%rot_mat_evals(row_blk_offset(row) + i - 1)
               e2 = qs_ot_env%rot_mat_evals(col_blk_offset(col) + j - 1)
               re_part = 0.0_dp
               im_part = 0.0_dp
               IF (found_in_re) re_part = block_in_re(i, j)
               IF (found_in_im) im_part = block_in_im(i, j)
               cval_in = CMPLX(re_part, im_part, dp)
               cval_out = cval_in*cint(e1, e2)
               block_out_re(i, j) = REAL(cval_out)
               block_out_im(i, j) = AIMAG(cval_out)
            END DO
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)
         DEALLOCATE (rows, cols)
         CALL dbcsr_release(inner_deriv_re)
         CALL dbcsr_release(inner_deriv_im)

         ! Compute: matrix_buf1 = \sum_i outer_deriv_ij * |rot_mat_evec_i> <rot_mat_evec_j|
         CALL dbcsr_multiply('N', 'N', +1.0_dp, qs_ot_env%rot_mat_evec_re, outer_deriv_re, 0.0_dp, mat_buf)
         CALL dbcsr_multiply('N', 'N', -1.0_dp, qs_ot_env%rot_mat_evec_im, outer_deriv_im, 1.0_dp, mat_buf)
         CALL dbcsr_multiply('N', 'T', +1.0_dp, mat_buf, qs_ot_env%rot_mat_evec_re, 0.0_dp, qs_ot_env%matrix_buf1)

         CALL dbcsr_multiply('N', 'N', +1.0_dp, qs_ot_env%rot_mat_evec_re, outer_deriv_im, 0.0_dp, mat_buf)
         CALL dbcsr_multiply('N', 'N', +1.0_dp, qs_ot_env%rot_mat_evec_im, outer_deriv_re, 1.0_dp, mat_buf)
         CALL dbcsr_multiply('N', 'T', +1.0_dp, mat_buf, qs_ot_env%rot_mat_evec_im, 1.0_dp, qs_ot_env%matrix_buf1)

         ! Account for anti-symmetry of rot_mat_x without relying on
         ! STRICT-only matrix buffers. REF uses the same finite chart.
         CALL qs_ot_square_transpose(qs_ot_env%matrix_buf1, mat_buf, qs_ot_env%rot_mat_u)
         CALL dbcsr_copy(qs_ot_env%rot_mat_gx, mat_buf)
         CALL dbcsr_add(qs_ot_env%rot_mat_gx, qs_ot_env%matrix_buf1, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

         CALL dbcsr_release(mat_buf)
         CALL dbcsr_release(outer_deriv_re)
         CALL dbcsr_release(outer_deriv_im)
      END IF
      CALL timestop(handle)
   CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param e1 ...
!> \param e2 ...
!> \return ...
! **************************************************************************************************
      FUNCTION cint(e1, e2)
      REAL(KIND=dp)                                      :: e1, e2
      COMPLEX(KIND=dp)                                   :: cint

      COMPLEX(KIND=dp)                                   :: l1, l2, x
      INTEGER                                            :: I

         l1 = (0.0_dp, -1.0_dp)*e1
         l2 = (0.0_dp, -1.0_dp)*e2
         IF (ABS(l1 - l2) > 0.5_dp) THEN
            cint = (EXP(l1) - EXP(l2))/(l1 - l2)
         ELSE
            x = 1.0_dp
            cint = 0.0_dp
            DO I = 1, 16
               cint = cint + x
               x = x*(l1 - l2)/REAL(I + 1, KIND=dp)
            END DO
            cint = cint*EXP(l2)
         END IF
      END FUNCTION cint
   END SUBROUTINE qs_ot_rot_mat_derivative

! **************************************************************************************************
!> \brief decide strategy
!>        tries to decide if the taylor expansion of cos(sqrt(xsx)) converges rapidly enough
!>        to make a taylor expansion of the functions cos(sqrt(xsx)) and sin(sqrt(xsx))/sqrt(xsx)
!>        and their derivatives faster than their computation based on diagonalization since xsx can
!>        be very small, especially during dynamics, only a few terms might indeed be needed we find
!>        the necessary order N to have largest_eval_upper_bound**(N+1)/(2(N+1))! < eps_taylor
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE decide_strategy(qs_ot_env)
      TYPE(qs_ot_type)                                   :: qs_ot_env

      INTEGER                                            :: N
      REAL(KIND=dp)                                      :: num_error

      qs_ot_env%do_taylor = .FALSE.
      N = 0
      num_error = qs_ot_env%largest_eval_upper_bound/(2.0_dp)
      DO WHILE (num_error > qs_ot_env%settings%eps_taylor .AND. N <= qs_ot_env%settings%max_taylor)
         N = N + 1
         num_error = num_error*qs_ot_env%largest_eval_upper_bound/REAL((2*N + 1)*(2*N + 2), KIND=dp)
      END DO
      qs_ot_env%taylor_order = N
      IF (qs_ot_env%taylor_order <= qs_ot_env%settings%max_taylor) THEN
         qs_ot_env%do_taylor = .TRUE.
      END IF

   END SUBROUTINE decide_strategy

! **************************************************************************************************
!> \brief c=(c0*cos(p^0.5)+x*sin(p^0.5)*p^(-0.5)) x rot_mat_u
!>        this assumes that x is already ortho to S*C0, and that p is x*S*x
!>        rot_mat_u is an optional rotation matrix
!> \param matrix_c ...
!> \param matrix_x ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE qs_ot_get_orbitals(matrix_c, matrix_x, qs_ot_env)

      TYPE(dbcsr_type), POINTER                          :: matrix_c, matrix_x
      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_orbitals'
      REAL(KIND=dp), PARAMETER                           :: rone = 1.0_dp, rzero = 0.0_dp

      INTEGER                                            :: handle, k, n
      TYPE(dbcsr_type), POINTER                          :: matrix_kk

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix_x, nfullrows_total=n, nfullcols_total=k)

      ! rotate the multiplying matrices cosp and sinp instead of the result,
      ! this should be cheaper for large basis sets
      IF (qs_ot_env%settings%do_rotation) THEN
         matrix_kk => qs_ot_env%matrix_buf1
         CALL dbcsr_multiply('N', 'N', rone, qs_ot_env%matrix_cosp, &
                             qs_ot_env%rot_mat_u, rzero, matrix_kk)
      ELSE
         matrix_kk => qs_ot_env%matrix_cosp
      END IF

      CALL dbcsr_multiply('N', 'N', rone, qs_ot_env%matrix_c0, matrix_kk, &
                          rzero, matrix_c)

      IF (qs_ot_env%settings%do_rotation) THEN
         matrix_kk => qs_ot_env%matrix_buf1
         CALL dbcsr_multiply('N', 'N', rone, qs_ot_env%matrix_sinp, &
                             qs_ot_env%rot_mat_u, rzero, matrix_kk)
      ELSE
         matrix_kk => qs_ot_env%matrix_sinp
      END IF
      CALL dbcsr_multiply('N', 'N', rone, matrix_x, matrix_kk, &
                          rone, matrix_c)

      CALL timestop(handle)

   END SUBROUTINE qs_ot_get_orbitals

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

      CALL timestop(handle)

   END SUBROUTINE qs_ot_get_orbitals_complex

! **************************************************************************************************
!> \brief this routines computes dE/dx=dx, with dx ortho to sc0
!>        needs dE/dC=hc,C0,X,SX,p
!>        if preconditioned it will not be the derivative, but the lagrangian multiplier
!>        is changed so that P*dE/dx is the right derivative (i.e. in the allowed subspace)
!> \param matrix_hc ...
!> \param matrix_x ...
!> \param matrix_sx ...
!> \param matrix_gx ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE qs_ot_get_derivative(matrix_hc, matrix_x, matrix_sx, matrix_gx, &
                                   qs_ot_env)
      TYPE(dbcsr_type), POINTER                          :: matrix_hc, matrix_x, matrix_sx, matrix_gx
      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative'
      REAL(KIND=dp), PARAMETER                           :: rone = 1.0_dp, rzero = 0.0_dp

      INTEGER                                            :: handle, k, n, ortho_k
      TYPE(dbcsr_type), POINTER                          :: matrix_hc_local, matrix_target

      CALL timeset(routineN, handle)

      NULLIFY (matrix_hc_local)

      CALL dbcsr_get_info(matrix_x, nfullrows_total=n, nfullcols_total=k)

      ! could in principle be taken inside qs_ot_get_derivative_* for increased efficiency
      ! create a local rotated version of matrix_hc leaving matrix_hc untouched (needed
      ! for lagrangian multipliers)
      IF (qs_ot_env%settings%do_rotation) THEN
         CALL dbcsr_copy(matrix_gx, matrix_hc) ! use gx as temporary
         CALL dbcsr_init_p(matrix_hc_local)
         CALL dbcsr_copy(matrix_hc_local, matrix_hc, name='matrix_hc_local')
         CALL dbcsr_set(matrix_hc_local, 0.0_dp)
         CALL dbcsr_multiply('N', 'T', rone, matrix_gx, qs_ot_env%rot_mat_u, rzero, matrix_hc_local)
      ELSE
         matrix_hc_local => matrix_hc
      END IF

      IF (qs_ot_env%do_taylor) THEN
         CALL qs_ot_get_derivative_taylor(matrix_hc_local, matrix_x, matrix_sx, matrix_gx, qs_ot_env)
      ELSE
         CALL qs_ot_get_derivative_diag(matrix_hc_local, matrix_x, matrix_sx, matrix_gx, qs_ot_env)
      END IF

      ! and make it orthogonal
      CALL dbcsr_get_info(qs_ot_env%matrix_sc0, nfullcols_total=ortho_k)

      IF (ASSOCIATED(qs_ot_env%preconditioner)) THEN
         matrix_target => qs_ot_env%matrix_psc0
      ELSE
         matrix_target => qs_ot_env%matrix_sc0
      END IF
      ! first make the matrix os if not yet valid
      IF (.NOT. qs_ot_env%os_valid) THEN
         ! this assumes that the preconditioner is a single matrix
         ! that maps sc0 onto psc0

         IF (ASSOCIATED(qs_ot_env%preconditioner)) THEN
            CALL apply_preconditioner(qs_ot_env%preconditioner, qs_ot_env%matrix_sc0, &
                                      qs_ot_env%matrix_psc0)
         END IF
         CALL dbcsr_multiply('T', 'N', rone, &
                             qs_ot_env%matrix_sc0, matrix_target, &
                             rzero, qs_ot_env%matrix_os)
         CALL cp_dbcsr_cholesky_decompose(qs_ot_env%matrix_os, &
                                          para_env=qs_ot_env%para_env, blacs_env=qs_ot_env%blacs_env)
         CALL cp_dbcsr_cholesky_invert(qs_ot_env%matrix_os, &
                                       para_env=qs_ot_env%para_env, blacs_env=qs_ot_env%blacs_env, &
                                       uplo_to_full=.TRUE.)
         qs_ot_env%os_valid = .TRUE.
      END IF
      CALL dbcsr_multiply('T', 'N', rone, matrix_target, matrix_gx, &
                          rzero, qs_ot_env%matrix_buf1_ortho)
      CALL dbcsr_multiply('N', 'N', rone, qs_ot_env%matrix_os, &
                          qs_ot_env%matrix_buf1_ortho, rzero, qs_ot_env%matrix_buf2_ortho)
      CALL dbcsr_multiply('N', 'N', -rone, qs_ot_env%matrix_sc0, &
                          qs_ot_env%matrix_buf2_ortho, rone, matrix_gx)
      ! also treat the rot_mat gradient here
      IF (qs_ot_env%settings%do_rotation) THEN
         CALL qs_ot_rot_mat_derivative(qs_ot_env)
      END IF

      IF (qs_ot_env%settings%do_rotation) THEN
         CALL dbcsr_release_p(matrix_hc_local)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_ot_get_derivative

! **************************************************************************************************
!> \brief finite complex STRICT derivative, projected onto C0^H*S*X=0
!> \param matrix_hc real part of H(k)*C(k)
!> \param matrix_hc_im imaginary part of H(k)*C(k)
!> \param qs_ot_env OT channel state
! **************************************************************************************************
   SUBROUTINE qs_ot_get_derivative_complex(matrix_hc, matrix_hc_im, qs_ot_env)
      TYPE(dbcsr_type), POINTER                          :: matrix_hc, matrix_hc_im
      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative_complex'

      INTEGER                                            :: handle, i, k
      REAL(KIND=dp)                                      :: eval_scale, eval_threshold
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_type), POINTER                          :: target_im, target_re

      CALL timeset(routineN, handle)

      CPASSERT(qs_ot_env%has_complex_kpoint_state)
      CPASSERT(qs_ot_env%settings%ot_algorithm == 'TOD')
      CPASSERT(.NOT. qs_ot_env%settings%do_rotation)

      ! Direct X contribution, H*C sinc(sqrt(P)).
      CALL qs_ot_complex_multiply('N', 'N', matrix_hc, matrix_hc_im, &
                                  qs_ot_env%matrix_sinp, qs_ot_env%matrix_sinp_im, &
                                  qs_ot_env%matrix_gx, qs_ot_env%matrix_gx_im, &
                                  qs_ot_env%matrix_tmp_nk)

      ! Frechet contribution from X sinc(sqrt(P)).
      CALL qs_ot_complex_multiply('C', 'N', matrix_hc, matrix_hc_im, &
                                  qs_ot_env%matrix_x, qs_ot_env%matrix_x_im, &
                                  qs_ot_env%matrix_buf2, qs_ot_env%matrix_buf2_im, &
                                  qs_ot_env%matrix_buf1)
      CALL qs_ot_complex_multiply('N', 'N', qs_ot_env%matrix_buf2, qs_ot_env%matrix_buf2_im, &
                                  qs_ot_env%matrix_r, qs_ot_env%matrix_r_im, &
                                  qs_ot_env%matrix_buf1, qs_ot_env%matrix_buf1_im, &
                                  qs_ot_env%matrix_buf4)
      CALL qs_ot_complex_multiply('C', 'N', qs_ot_env%matrix_r, qs_ot_env%matrix_r_im, &
                                  qs_ot_env%matrix_buf1, qs_ot_env%matrix_buf1_im, &
                                  qs_ot_env%matrix_buf2, qs_ot_env%matrix_buf2_im, &
                                  qs_ot_env%matrix_buf4)
      CALL dbcsr_hadamard_product(qs_ot_env%matrix_buf2, qs_ot_env%matrix_sinp_b, &
                                  qs_ot_env%matrix_buf3)
      CALL dbcsr_hadamard_product(qs_ot_env%matrix_buf2_im, qs_ot_env%matrix_sinp_b, &
                                  qs_ot_env%matrix_buf3_im)

      ! Frechet contribution from C0 cos(sqrt(P)).
      CALL qs_ot_complex_multiply('C', 'N', matrix_hc, matrix_hc_im, &
                                  qs_ot_env%matrix_c0, qs_ot_env%matrix_c0_im, &
                                  qs_ot_env%matrix_buf2, qs_ot_env%matrix_buf2_im, &
                                  qs_ot_env%matrix_buf1)
      CALL qs_ot_complex_multiply('N', 'N', qs_ot_env%matrix_buf2, qs_ot_env%matrix_buf2_im, &
                                  qs_ot_env%matrix_r, qs_ot_env%matrix_r_im, &
                                  qs_ot_env%matrix_buf1, qs_ot_env%matrix_buf1_im, &
                                  qs_ot_env%matrix_buf4)
      CALL qs_ot_complex_multiply('C', 'N', qs_ot_env%matrix_r, qs_ot_env%matrix_r_im, &
                                  qs_ot_env%matrix_buf1, qs_ot_env%matrix_buf1_im, &
                                  qs_ot_env%matrix_buf2, qs_ot_env%matrix_buf2_im, &
                                  qs_ot_env%matrix_buf4)
      CALL dbcsr_hadamard_product(qs_ot_env%matrix_buf2, qs_ot_env%matrix_cosp_b, &
                                  qs_ot_env%matrix_buf4)
      CALL dbcsr_hadamard_product(qs_ot_env%matrix_buf2_im, qs_ot_env%matrix_cosp_b, &
                                  qs_ot_env%matrix_buf4_im)
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

      ! Preconditioner-aware projection onto the complex STRICT tangent space.
      IF (ASSOCIATED(qs_ot_env%preconditioner)) THEN
         target_re => qs_ot_env%matrix_psc0
         target_im => qs_ot_env%matrix_psc0_im
      ELSE
         target_re => qs_ot_env%matrix_sc0
         target_im => qs_ot_env%matrix_sc0_im
      END IF
      IF (.NOT. qs_ot_env%os_valid) THEN
         IF (ASSOCIATED(qs_ot_env%preconditioner)) THEN
            CALL apply_preconditioner(qs_ot_env%preconditioner, &
                                      qs_ot_env%matrix_sc0, qs_ot_env%matrix_sc0_im, &
                                      qs_ot_env%matrix_psc0, qs_ot_env%matrix_psc0_im)
         END IF
         CALL qs_ot_complex_multiply('C', 'N', qs_ot_env%matrix_sc0, qs_ot_env%matrix_sc0_im, &
                                     target_re, target_im, qs_ot_env%matrix_os, qs_ot_env%matrix_os_im, &
                                     qs_ot_env%matrix_buf1)
         CALL dbcsr_get_info(qs_ot_env%matrix_os, nfullrows_total=k)
         CALL cp_dbcsr_heevd(matrix_re=qs_ot_env%matrix_os, matrix_im=qs_ot_env%matrix_os_im, &
                             eigenvectors_re=qs_ot_env%matrix_r, eigenvectors_im=qs_ot_env%matrix_r_im, &
                             eigenvalues=qs_ot_env%evals, para_env=qs_ot_env%para_env, &
                             blacs_env=qs_ot_env%blacs_env)
         eval_scale = MAX(1.0_dp, MAXVAL(ABS(qs_ot_env%evals(1:k))))
         eval_threshold = 100.0_dp*EPSILON(1.0_dp)*eval_scale
         IF (MINVAL(qs_ot_env%evals(1:k)) <= eval_threshold) THEN
            CPABORT("Complex STRICT tangent metric is not positive definite")
         END IF
         DO i = 1, k
            qs_ot_env%dum(i) = 1.0_dp/qs_ot_env%evals(i)
         END DO
         CALL dbcsr_copy(qs_ot_env%matrix_buf1, qs_ot_env%matrix_r)
         CALL dbcsr_copy(qs_ot_env%matrix_buf1_im, qs_ot_env%matrix_r_im)
         CALL dbcsr_scale_by_vector(qs_ot_env%matrix_buf1, alpha=qs_ot_env%dum, side='right')
         CALL dbcsr_scale_by_vector(qs_ot_env%matrix_buf1_im, alpha=qs_ot_env%dum, side='right')
         CALL qs_ot_complex_multiply('N', 'C', qs_ot_env%matrix_buf1, qs_ot_env%matrix_buf1_im, &
                                     qs_ot_env%matrix_r, qs_ot_env%matrix_r_im, &
                                     qs_ot_env%matrix_os, qs_ot_env%matrix_os_im, &
                                     qs_ot_env%matrix_buf2)
         qs_ot_env%os_valid = .TRUE.
      END IF

      CALL qs_ot_complex_multiply('C', 'N', target_re, target_im, &
                                  qs_ot_env%matrix_gx, qs_ot_env%matrix_gx_im, &
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
      CALL dbcsr_add(qs_ot_env%matrix_gx, qs_ot_env%matrix_buf_nk, &
                     alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      CALL dbcsr_add(qs_ot_env%matrix_gx_im, qs_ot_env%matrix_buf_nk_im, &
                     alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)

      CALL timestop(handle)

   END SUBROUTINE qs_ot_get_derivative_complex

! **************************************************************************************************
!> \brief ...
!> \param matrix_hc ...
!> \param matrix_x ...
!> \param matrix_sx ...
!> \param matrix_gx ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE qs_ot_get_derivative_diag(matrix_hc, matrix_x, matrix_sx, matrix_gx, &
                                        qs_ot_env)

      TYPE(dbcsr_type), POINTER                          :: matrix_hc, matrix_x, matrix_sx, matrix_gx
      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative_diag'
      REAL(KIND=dp), PARAMETER                           :: rone = 1.0_dp, rzero = 0.0_dp

      INTEGER                                            :: handle, k, n
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix_x, nfullrows_total=n, nfullcols_total=k)

      ! go for the derivative now
      ! this de/dc*(dX/dx)*sinp
      CALL dbcsr_multiply('N', 'N', rone, matrix_hc, qs_ot_env%matrix_sinp, rzero, matrix_gx)
      ! overlap hc*x
      CALL dbcsr_multiply('T', 'N', rone, matrix_hc, matrix_x, rzero, qs_ot_env%matrix_buf2)
      ! get it in the basis of the eigenvectors
      CALL dbcsr_multiply('N', 'N', rone, qs_ot_env%matrix_buf2, qs_ot_env%matrix_r, &
                          rzero, qs_ot_env%matrix_buf1)
      CALL dbcsr_multiply('T', 'N', rone, qs_ot_env%matrix_r, qs_ot_env%matrix_buf1, &
                          rzero, qs_ot_env%matrix_buf2)

      ! get the schur product of O_uv*B_uv
      CALL dbcsr_hadamard_product(qs_ot_env%matrix_buf2, qs_ot_env%matrix_sinp_b, &
                                  qs_ot_env%matrix_buf3)

      ! overlap hc*c0
      CALL dbcsr_multiply('T', 'N', rone, matrix_hc, qs_ot_env%matrix_c0, rzero, &
                          qs_ot_env%matrix_buf2)
      ! get it in the basis of the eigenvectors
      CALL dbcsr_multiply('N', 'N', rone, qs_ot_env%matrix_buf2, qs_ot_env%matrix_r, &
                          rzero, qs_ot_env%matrix_buf1)
      CALL dbcsr_multiply('T', 'N', rone, qs_ot_env%matrix_r, qs_ot_env%matrix_buf1, &
                          rzero, qs_ot_env%matrix_buf2)

      CALL dbcsr_hadamard_product(qs_ot_env%matrix_buf2, qs_ot_env%matrix_cosp_b, &
                                  qs_ot_env%matrix_buf4)

      ! add the two bs and compute b+b^T
      CALL dbcsr_add(qs_ot_env%matrix_buf3, qs_ot_env%matrix_buf4, &
                     alpha_scalar=rone, beta_scalar=rone)

      ! get the b in the eigenvector basis
      CALL dbcsr_multiply('N', 'T', rone, qs_ot_env%matrix_buf3, qs_ot_env%matrix_r, &
                          rzero, qs_ot_env%matrix_buf1)
      CALL dbcsr_multiply('N', 'N', rone, qs_ot_env%matrix_r, qs_ot_env%matrix_buf1, &
                          rzero, qs_ot_env%matrix_buf3)
      CALL dbcsr_get_info(qs_ot_env%matrix_buf3, distribution=dist)
      CALL dbcsr_transposed(qs_ot_env%matrix_buf1, qs_ot_env%matrix_buf3, &
                            shallow_data_copy=.FALSE., use_distribution=dist, &
                            transpose_distribution=.FALSE.)
      CALL dbcsr_add(qs_ot_env%matrix_buf3, qs_ot_env%matrix_buf1, &
                     alpha_scalar=rone, beta_scalar=rone)

      ! and add to the derivative
      CALL dbcsr_multiply('N', 'N', rone, matrix_sx, qs_ot_env%matrix_buf3, &
                          rone, matrix_gx)
      CALL timestop(handle)

   END SUBROUTINE qs_ot_get_derivative_diag

! **************************************************************************************************
!> \brief compute the derivative of the taylor expansion below
!> \param matrix_hc ...
!> \param matrix_x ...
!> \param matrix_sx ...
!> \param matrix_gx ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE qs_ot_get_derivative_taylor(matrix_hc, matrix_x, matrix_sx, matrix_gx, &
                                          qs_ot_env)

      TYPE(dbcsr_type), POINTER                          :: matrix_hc, matrix_x, matrix_sx, matrix_gx
      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative_taylor'
      REAL(KIND=dp), PARAMETER                           :: rone = 1.0_dp, rzero = 0.0_dp

      INTEGER                                            :: handle, i, k, n
      REAL(KIND=dp)                                      :: cosfactor, sinfactor
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_type), POINTER                          :: matrix_left, matrix_right

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(matrix_x, nfullrows_total=n, nfullcols_total=k)

      ! go for the derivative now
      ! this de/dc*(dX/dx)*sinp i.e. zeroth order
      CALL dbcsr_multiply('N', 'N', rone, matrix_hc, qs_ot_env%matrix_sinp, rzero, matrix_gx)

      IF (qs_ot_env%taylor_order <= 0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! we store the matrix that will multiply sx in matrix_r
      CALL dbcsr_set(qs_ot_env%matrix_r, rzero)

      ! just better names for matrix_cosp_b and matrix_sinp_b (they are buffer space here)
      matrix_left => qs_ot_env%matrix_cosp_b
      matrix_right => qs_ot_env%matrix_sinp_b

      ! overlap hc*x and add its transpose to matrix_left
      CALL dbcsr_multiply('T', 'N', rone, matrix_hc, matrix_x, rzero, matrix_left)
      CALL dbcsr_get_info(matrix_left, distribution=dist)
      CALL dbcsr_transposed(qs_ot_env%matrix_buf1, matrix_left, &
                            shallow_data_copy=.FALSE., use_distribution=dist, &
                            transpose_distribution=.FALSE.)
      CALL dbcsr_add(matrix_left, qs_ot_env%matrix_buf1, &
                     alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      CALL dbcsr_copy(matrix_right, matrix_left)

      ! first order
      sinfactor = -1.0_dp/(2.0_dp*3.0_dp)
      CALL dbcsr_add(qs_ot_env%matrix_r, matrix_left, &
                     alpha_scalar=1.0_dp, beta_scalar=sinfactor)

      !      M
      !    OM+MO
      ! OOM+OMO+MOO
      !   ...
      DO i = 2, qs_ot_env%taylor_order
         sinfactor = sinfactor*(-1.0_dp)/REAL(2*i*(2*i + 1), KIND=dp)
         CALL dbcsr_multiply('N', 'N', rone, qs_ot_env%matrix_p, matrix_left, rzero, qs_ot_env%matrix_buf1)
         CALL dbcsr_multiply('N', 'N', rone, matrix_right, qs_ot_env%matrix_p, rzero, matrix_left)
         CALL dbcsr_copy(matrix_right, matrix_left)
         CALL dbcsr_add(matrix_left, qs_ot_env%matrix_buf1, &
                        1.0_dp, 1.0_dp)
         CALL dbcsr_add(qs_ot_env%matrix_r, matrix_left, &
                        alpha_scalar=1.0_dp, beta_scalar=sinfactor)
      END DO

      ! overlap hc*c0 and add its transpose to matrix_left
      CALL dbcsr_multiply('T', 'N', rone, matrix_hc, qs_ot_env%matrix_c0, rzero, matrix_left)
      CALL dbcsr_get_info(matrix_left, distribution=dist)
      CALL dbcsr_transposed(qs_ot_env%matrix_buf1, matrix_left, &
                            shallow_data_copy=.FALSE., use_distribution=dist, &
                            transpose_distribution=.FALSE.)
      CALL dbcsr_add(matrix_left, qs_ot_env%matrix_buf1, 1.0_dp, 1.0_dp)
      CALL dbcsr_copy(matrix_right, matrix_left)

      ! first order
      cosfactor = -1.0_dp/(1.0_dp*2.0_dp)
      CALL dbcsr_add(qs_ot_env%matrix_r, matrix_left, &
                     alpha_scalar=1.0_dp, beta_scalar=cosfactor)

      !      M
      !    OM+MO
      ! OOM+OMO+MOO
      !   ...
      DO i = 2, qs_ot_env%taylor_order
         cosfactor = cosfactor*(-1.0_dp)/REAL(2*i*(2*i - 1), KIND=dp)
         CALL dbcsr_multiply('N', 'N', rone, qs_ot_env%matrix_p, matrix_left, rzero, qs_ot_env%matrix_buf1)
         CALL dbcsr_multiply('N', 'N', rone, matrix_right, qs_ot_env%matrix_p, rzero, matrix_left)
         CALL dbcsr_copy(matrix_right, matrix_left)
         CALL dbcsr_add(matrix_left, qs_ot_env%matrix_buf1, 1.0_dp, 1.0_dp)
         CALL dbcsr_add(qs_ot_env%matrix_r, matrix_left, &
                        alpha_scalar=1.0_dp, beta_scalar=cosfactor)
      END DO

      ! and add to the derivative
      CALL dbcsr_multiply('N', 'N', rone, matrix_sx, qs_ot_env%matrix_r, rone, matrix_gx)

      CALL timestop(handle)

   END SUBROUTINE qs_ot_get_derivative_taylor

! *************************************************************************************************
!> \brief computes a taylor expansion.
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE qs_ot_p2m_taylor(qs_ot_env)
      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_ot_p2m_taylor'
      REAL(KIND=dp), PARAMETER                           :: rone = 1.0_dp, rzero = 0.0_dp

      INTEGER                                            :: handle, i, k
      REAL(KIND=dp)                                      :: cosfactor, sinfactor

      CALL timeset(routineN, handle)

      ! zeroth order
      CALL dbcsr_set(qs_ot_env%matrix_cosp, rzero)
      CALL dbcsr_set(qs_ot_env%matrix_sinp, rzero)
      CALL dbcsr_add_on_diag(qs_ot_env%matrix_cosp, rone)
      CALL dbcsr_add_on_diag(qs_ot_env%matrix_sinp, rone)

      IF (qs_ot_env%taylor_order <= 0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! first order
      cosfactor = -1.0_dp/(1.0_dp*2.0_dp)
      sinfactor = -1.0_dp/(2.0_dp*3.0_dp)
      CALL dbcsr_add(qs_ot_env%matrix_cosp, qs_ot_env%matrix_p, alpha_scalar=1.0_dp, beta_scalar=cosfactor)
      CALL dbcsr_add(qs_ot_env%matrix_sinp, qs_ot_env%matrix_p, alpha_scalar=1.0_dp, beta_scalar=sinfactor)
      IF (qs_ot_env%taylor_order <= 1) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! other orders
      CALL dbcsr_get_info(qs_ot_env%matrix_p, nfullrows_total=k)
      CALL dbcsr_copy(qs_ot_env%matrix_r, qs_ot_env%matrix_p)

      DO i = 2, qs_ot_env%taylor_order
         ! new power of p
         CALL dbcsr_multiply('N', 'N', rone, qs_ot_env%matrix_p, qs_ot_env%matrix_r, &
                             rzero, qs_ot_env%matrix_buf1)
         CALL dbcsr_copy(qs_ot_env%matrix_r, qs_ot_env%matrix_buf1)
         ! add to the taylor expansion so far
         cosfactor = cosfactor*(-1.0_dp)/REAL(2*i*(2*i - 1), KIND=dp)
         sinfactor = sinfactor*(-1.0_dp)/REAL(2*i*(2*i + 1), KIND=dp)
         CALL dbcsr_add(qs_ot_env%matrix_cosp, qs_ot_env%matrix_r, &
                        alpha_scalar=1.0_dp, beta_scalar=cosfactor)
         CALL dbcsr_add(qs_ot_env%matrix_sinp, qs_ot_env%matrix_r, &
                        alpha_scalar=1.0_dp, beta_scalar=sinfactor)
      END DO

      CALL timestop(handle)

   END SUBROUTINE qs_ot_p2m_taylor

! **************************************************************************************************
!> \brief given p, computes  - eigenstuff (matrix_r,evals)
!>        - cos(p^0.5),p^(-0.5)*sin(p^0.5)
!>        - the real b matrices, needed for the derivatives of these guys
!>        cosp_b_ij=(1/(2pii) * int(cos(z^1/2)/((z-eval(i))*(z-eval(j))))
!>        sinp_b_ij=(1/(2pii) * int(z^(-1/2)*sin(z^1/2)/((z-eval(i))*(z-eval(j))))
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE qs_ot_p2m_diag(qs_ot_env)

      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_ot_p2m_diag'
      REAL(KIND=dp), PARAMETER                           :: rone = 1.0_dp, rzero = 0.0_dp

      INTEGER                                            :: col, col_offset, col_size, handle, i, j, &
                                                            k, row, row_offset, row_size
      REAL(dp), DIMENSION(:, :), POINTER                 :: block
      REAL(KIND=dp)                                      :: a, b
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(qs_ot_env%matrix_p, nfullrows_total=k)
      CALL dbcsr_copy(qs_ot_env%matrix_buf1, qs_ot_env%matrix_p)
      CALL cp_dbcsr_syevd(qs_ot_env%matrix_buf1, qs_ot_env%matrix_r, qs_ot_env%evals, &
                          qs_ot_env%para_env, qs_ot_env%blacs_env)
      DO i = 1, k
         qs_ot_env%evals(i) = MAX(0.0_dp, qs_ot_env%evals(i))
      END DO

      !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(k,qs_ot_env)
      DO i = 1, k
         qs_ot_env%dum(i) = COS(SQRT(qs_ot_env%evals(i)))
      END DO
      CALL dbcsr_copy(qs_ot_env%matrix_buf1, qs_ot_env%matrix_r)
      CALL dbcsr_scale_by_vector(qs_ot_env%matrix_buf1, alpha=qs_ot_env%dum, side='right')
      CALL dbcsr_multiply('N', 'T', rone, qs_ot_env%matrix_r, qs_ot_env%matrix_buf1, &
                          rzero, qs_ot_env%matrix_cosp)

      !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(k,qs_ot_env)
      DO i = 1, k
         qs_ot_env%dum(i) = qs_ot_sinc(SQRT(qs_ot_env%evals(i)))
      END DO
      CALL dbcsr_copy(qs_ot_env%matrix_buf1, qs_ot_env%matrix_r)
      CALL dbcsr_scale_by_vector(qs_ot_env%matrix_buf1, alpha=qs_ot_env%dum, side='right')
      CALL dbcsr_multiply('N', 'T', rone, qs_ot_env%matrix_r, qs_ot_env%matrix_buf1, &
                          rzero, qs_ot_env%matrix_sinp)

      CALL dbcsr_copy(qs_ot_env%matrix_cosp_b, qs_ot_env%matrix_cosp)
      CALL dbcsr_iterator_start(iter, qs_ot_env%matrix_cosp_b)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         DO j = 1, col_size
         DO i = 1, row_size
            a = (SQRT(qs_ot_env%evals(row_offset + i - 1)) &
                 - SQRT(qs_ot_env%evals(col_offset + j - 1)))/2.0_dp
            b = (SQRT(qs_ot_env%evals(row_offset + i - 1)) &
                 + SQRT(qs_ot_env%evals(col_offset + j - 1)))/2.0_dp
            block(i, j) = -0.5_dp*qs_ot_sinc(a)*qs_ot_sinc(b)
         END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_copy(qs_ot_env%matrix_sinp_b, qs_ot_env%matrix_sinp)
      CALL dbcsr_iterator_start(iter, qs_ot_env%matrix_sinp_b)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         DO j = 1, col_size
         DO i = 1, row_size
            a = SQRT(qs_ot_env%evals(row_offset + i - 1))
            b = SQRT(qs_ot_env%evals(col_offset + j - 1))
            block(i, j) = qs_ot_sincf(a, b)
         END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE qs_ot_p2m_diag

! **************************************************************************************************
!> \brief diagonalize Hermitian P and build cos(sqrt(P)), sinc(sqrt(P)), and Frechet kernels
!> \param qs_ot_env complex STRICT channel state
! **************************************************************************************************
   SUBROUTINE qs_ot_p2m_diag_complex(qs_ot_env)
      TYPE(qs_ot_type)                                   :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_p2m_diag_complex'

      INTEGER                                            :: col, col_offset, col_size, handle, i, j, &
                                                            k, row, row_offset, row_size
      REAL(dp), DIMENSION(:, :), POINTER                 :: block
      REAL(KIND=dp)                                      :: a, b
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(qs_ot_env%matrix_p, nfullrows_total=k)
      CALL cp_dbcsr_heevd(matrix_re=qs_ot_env%matrix_p, matrix_im=qs_ot_env%matrix_p_im, &
                          eigenvectors_re=qs_ot_env%matrix_r, eigenvectors_im=qs_ot_env%matrix_r_im, &
                          eigenvalues=qs_ot_env%evals, para_env=qs_ot_env%para_env, &
                          blacs_env=qs_ot_env%blacs_env)
      DO i = 1, k
         qs_ot_env%evals(i) = MAX(0.0_dp, qs_ot_env%evals(i))
      END DO

      DO i = 1, k
         qs_ot_env%dum(i) = COS(SQRT(qs_ot_env%evals(i)))
      END DO
      CALL dbcsr_copy(qs_ot_env%matrix_buf1, qs_ot_env%matrix_r)
      CALL dbcsr_copy(qs_ot_env%matrix_buf1_im, qs_ot_env%matrix_r_im)
      CALL dbcsr_scale_by_vector(qs_ot_env%matrix_buf1, alpha=qs_ot_env%dum, side='right')
      CALL dbcsr_scale_by_vector(qs_ot_env%matrix_buf1_im, alpha=qs_ot_env%dum, side='right')
      CALL qs_ot_complex_multiply('N', 'C', qs_ot_env%matrix_buf1, qs_ot_env%matrix_buf1_im, &
                                  qs_ot_env%matrix_r, qs_ot_env%matrix_r_im, &
                                  qs_ot_env%matrix_cosp, qs_ot_env%matrix_cosp_im, &
                                  qs_ot_env%matrix_buf2)

      DO i = 1, k
         qs_ot_env%dum(i) = qs_ot_sinc(SQRT(qs_ot_env%evals(i)))
      END DO
      CALL dbcsr_copy(qs_ot_env%matrix_buf1, qs_ot_env%matrix_r)
      CALL dbcsr_copy(qs_ot_env%matrix_buf1_im, qs_ot_env%matrix_r_im)
      CALL dbcsr_scale_by_vector(qs_ot_env%matrix_buf1, alpha=qs_ot_env%dum, side='right')
      CALL dbcsr_scale_by_vector(qs_ot_env%matrix_buf1_im, alpha=qs_ot_env%dum, side='right')
      CALL qs_ot_complex_multiply('N', 'C', qs_ot_env%matrix_buf1, qs_ot_env%matrix_buf1_im, &
                                  qs_ot_env%matrix_r, qs_ot_env%matrix_r_im, &
                                  qs_ot_env%matrix_sinp, qs_ot_env%matrix_sinp_im, &
                                  qs_ot_env%matrix_buf2)

      CALL dbcsr_copy(qs_ot_env%matrix_cosp_b, qs_ot_env%matrix_cosp)
      CALL dbcsr_iterator_start(iter, qs_ot_env%matrix_cosp_b)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         DO j = 1, col_size
         DO i = 1, row_size
            a = (SQRT(qs_ot_env%evals(row_offset + i - 1)) &
                 - SQRT(qs_ot_env%evals(col_offset + j - 1)))/2.0_dp
            b = (SQRT(qs_ot_env%evals(row_offset + i - 1)) &
                 + SQRT(qs_ot_env%evals(col_offset + j - 1)))/2.0_dp
            block(i, j) = -0.5_dp*qs_ot_sinc(a)*qs_ot_sinc(b)
         END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_copy(qs_ot_env%matrix_sinp_b, qs_ot_env%matrix_sinp)
      CALL dbcsr_iterator_start(iter, qs_ot_env%matrix_sinp_b)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, block, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         DO j = 1, col_size
         DO i = 1, row_size
            a = SQRT(qs_ot_env%evals(row_offset + i - 1))
            b = SQRT(qs_ot_env%evals(col_offset + j - 1))
            block(i, j) = qs_ot_sincf(a, b)
         END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE qs_ot_p2m_diag_complex

! **************************************************************************************************
!> \brief computes sin(x)/x for all values of the argument
!> \param x ...
!> \return ...
! **************************************************************************************************
   FUNCTION qs_ot_sinc(x)

      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp)                                      :: qs_ot_sinc

      REAL(KIND=dp), PARAMETER :: q1 = 1.0_dp, q2 = -q1/(2.0_dp*3.0_dp), q3 = -q2/(4.0_dp*5.0_dp), &
         q4 = -q3/(6.0_dp*7.0_dp), q5 = -q4/(8.0_dp*9.0_dp), q6 = -q5/(10.0_dp*11.0_dp), &
         q7 = -q6/(12.0_dp*13.0_dp), q8 = -q7/(14.0_dp*15.0_dp), q9 = -q8/(16.0_dp*17.0_dp), &
         q10 = -q9/(18.0_dp*19.0_dp)

      REAL(KIND=dp)                                      :: y

      IF (ABS(x) > 0.5_dp) THEN
         qs_ot_sinc = SIN(x)/x
      ELSE
         y = x*x
         qs_ot_sinc = q1 + y*(q2 + y*(q3 + y*(q4 + y*(q5 + y*(q6 + y*(q7 + y*(q8 + y*(q9 + y*(q10)))))))))
      END IF
   END FUNCTION qs_ot_sinc

! **************************************************************************************************
!> \brief computes (1/(x^2-y^2))*(sinc(x)-sinc(y)) for all positive values of the arguments
!> \param xa ...
!> \param ya ...
!> \return ...
! **************************************************************************************************
   FUNCTION qs_ot_sincf(xa, ya)

      REAL(KIND=dp), INTENT(IN)                          :: xa, ya
      REAL(KIND=dp)                                      :: qs_ot_sincf

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: a, b, rs, sf, x, xs, y, ybx, ybxs

      ! this is currently a limit of the routine, could be removed rather easily
      IF (xa < 0) CPABORT("x is negative")
      IF (ya < 0) CPABORT("y is negative")

      IF (xa < ya) THEN
         x = ya
         y = xa
      ELSE
         x = xa
         y = ya
      END IF

      IF (x < 0.5_dp) THEN ! use series, keeping in mind that x,y,x+y,x-y can all be zero

         qs_ot_sincf = 0.0_dp
         IF (x > 0.0_dp) THEN
            ybx = y/x
         ELSE ! should be irrelevant  !?
            ybx = 0.0_dp
         END IF

         sf = -1.0_dp/((1.0_dp + ybx)*6.0_dp)
         rs = 1.0_dp
         ybxs = ybx
         xs = 1.0_dp

         DO i = 1, 10
            qs_ot_sincf = qs_ot_sincf + sf*rs*xs*(1.0_dp + ybxs)
            sf = -sf/(REAL((2*i + 2), dp)*REAL((2*i + 3), dp))
            rs = rs + ybxs
            ybxs = ybxs*ybx
            xs = xs*x*x
         END DO

      ELSE ! no series expansion
         IF (x - y > 0.1_dp) THEN ! safe to use the normal form
            qs_ot_sincf = (qs_ot_sinc(x) - qs_ot_sinc(y))/((x + y)*(x - y))
         ELSE
            a = (x + y)/2.0_dp
            b = (x - y)/2.0_dp ! might be close to zero
            ! y (=(a-b)) can not be close to zero since it is close to x>0.5
            qs_ot_sincf = (qs_ot_sinc(b)*COS(a) - qs_ot_sinc(a)*COS(b))/(2*x*y)
         END IF
      END IF

   END FUNCTION qs_ot_sincf

END MODULE qs_ot

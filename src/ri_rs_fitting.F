!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shared numerical operations for computing the RI-RS matrix Z_lP.
!> \par History
!>      09.2026 created Jan Wilhelm
! **************************************************************************************************
MODULE ri_rs_fitting
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_put_block,&
                                              dbcsr_type
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_solve
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_rs_fitting'
   REAL(KIND=dp), PARAMETER, PRIVATE :: jacobi_floor = 1.0E-16_dp

   PUBLIC :: evaluate_ri_rs_fit, &
             solve_ri_rs_fit, &
             store_Z_lP_columns

CONTAINS

! **************************************************************************************************
!> \brief Computes the regularized RI-RS fit from an AO collocation matrix and a replicated RHS.
!>
!>          D_ll' = [Σ_μ ϕ_μ(r_l) ϕ_μ(r_l')]²,
!>          d_l = 1/sqrt(max(D_ll, jacobi_floor)),
!>          D'_ll' = d_l D_ll' d_l' + λ δ_ll',
!>          Σ_l' D'_ll' Z'_l'P = d_l d_lP,
!>          Z_lP = d_l Z'_lP.
!>
!>        The RI-RS least-squares construction and balancing are described by Duchemin and
!>        Blase, arXiv:1902.00066, Sec. II.3. Here jacobi_floor bounds D_ll from below in
!>        both the LAPACK and ScaLAPACK paths. The regularization λ is added after Jacobi
!>        scaling; in the unscaled system its diagonal contribution is λ/d_l².
!> \param phi_local AO values on the fitting grid, replicated within the solve subgroup
!> \param d_lp replicated RHS d_lP on entry and fitted coefficients Z_lP on return
!> \param tikhonov regularization parameter λ in the Jacobi-scaled space
!> \param para_env_sub optional subgroup communicator for the distributed solve
!> \param blacs_env_sub optional BLACS context, supplied together with para_env_sub
! **************************************************************************************************
   SUBROUTINE solve_ri_rs_fit(phi_local, d_lp, tikhonov, para_env_sub, blacs_env_sub)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi_local
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: d_lp
      REAL(KIND=dp), INTENT(IN)                          :: tikhonov
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_sub
      TYPE(cp_blacs_env_type), OPTIONAL, POINTER         :: blacs_env_sub

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'solve_ri_rs_fit'

      INTEGER                                            :: handle, handle2, i, info, n_grid, n_rhs
      REAL(KIND=dp), ALLOCATABLE                         :: gram(:, :), jacobi_scale(:)

      CALL timeset(routineN, handle)

      n_grid = SIZE(phi_local, 1)
      n_rhs = SIZE(d_lp, 2)
      CPASSERT(SIZE(d_lp, 1) == n_grid)
      CPASSERT(PRESENT(para_env_sub) .EQV. PRESENT(blacs_env_sub))
      IF (PRESENT(para_env_sub)) THEN
         ! Only the row norms are replicated; the distributed solve builds its own Gram tiles.
         ALLOCATE (jacobi_scale(n_grid))
         !$OMP PARALLEL DO DEFAULT(NONE) &
         !$OMP SHARED(n_grid, phi_local, jacobi_scale) &
         !$OMP PRIVATE(i) SCHEDULE(STATIC)
         DO i = 1, n_grid
            jacobi_scale(i) = 1.0_dp/MAX(DOT_PRODUCT(phi_local(i, :), phi_local(i, :)), SQRT(jacobi_floor))
         END DO
         !$OMP END PARALLEL DO
      ELSE
         CALL factor_ri_rs_fit(phi_local, tikhonov, gram, jacobi_scale, info)
         IF (info /= 0) CPABORT('RI-RS Cholesky factorization failed')
      END IF

      ! Pre-scale the right-hand side: d_lP <- d_l d_lP.
      CALL scale_rows_by_diag(d_lp, jacobi_scale)
      IF (PRESENT(para_env_sub)) THEN
         CALL solve_D_lp_distributed(phi_local, jacobi_scale, d_lp, tikhonov, &
                                     para_env_sub, blacs_env_sub)
      ELSE
         ! Solve Σ_l' D'_ll' Z'_l'P = d_l d_lP, overwriting d_lp with Z'_lP.
         CALL timeset(routineN//'_dpotrs', handle2)
         CALL dpotrs('L', n_grid, n_rhs, gram, n_grid, d_lp, n_grid, info)
         CALL timestop(handle2)
         IF (info /= 0) CPABORT('RI-RS Cholesky solve failed')
      END IF
      ! Undo the conditioning on the solution: Z_lP = d_l Z'_lP.
      CALL scale_rows_by_diag(d_lp, jacobi_scale)

      CALL timestop(handle)

   END SUBROUTINE solve_ri_rs_fit

! **************************************************************************************************
!> \brief Build and factor the Jacobi-scaled squared AO Gram matrix, K = d D d + λ I.
!>        The lower triangle is used throughout; the factor is retained for subsequent solves.
!>        See solve_ri_rs_fit for the RI-RS balancing convention and reference.
!> \param phi AO values indexed (grid point, AO)
!> \param tikhonov regularization in the Jacobi-scaled space
!> \param factor lower Cholesky factor of K on successful return
!> \param jacobi_scale pair-feature scales d_l
!> \param info LAPACK factorization status; the caller chooses how to handle failure
!> \param overlap optional lower triangle of the scaled AO overlap V V^T, V_l = sqrt(d_l) phi_l
! **************************************************************************************************
   SUBROUTINE factor_ri_rs_fit(phi, tikhonov, factor, jacobi_scale, info, overlap)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi
      REAL(KIND=dp), INTENT(IN)                          :: tikhonov
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: factor(:, :), jacobi_scale(:)
      INTEGER, INTENT(OUT)                               :: info
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT), OPTIONAL  :: overlap(:, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'factor_ri_rs_fit'

      INTEGER                                            :: handle, handle2, i, j, n_grid

      CALL timeset(routineN, handle)
      n_grid = SIZE(phi, 1)
      ALLOCATE (factor(n_grid, n_grid), jacobi_scale(n_grid))
      ! First form S_ll' = Σ_μ phi_lμ phi_l'μ; D_ll' = S_ll'^2.
      CALL timeset(routineN//'_dsyrk', handle2)
      CALL dsyrk('L', 'N', n_grid, SIZE(phi, 2), 1.0_dp, phi, n_grid, 0.0_dp, factor, n_grid)
      CALL timestop(handle2)
      DO i = 1, n_grid
         jacobi_scale(i) = 1.0_dp/MAX(factor(i, i), SQRT(jacobi_floor))
      END DO
      IF (PRESENT(overlap)) THEN
         ALLOCATE (overlap(n_grid, n_grid))
         DO j = 1, n_grid
            DO i = j, n_grid
               overlap(i, j) = factor(i, j)*SQRT(jacobi_scale(i))*SQRT(jacobi_scale(j))
            END DO
         END DO
      END IF
      ! Add λ after Jacobi scaling; LAPACK needs only the lower triangle.
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(n_grid, factor, jacobi_scale, tikhonov) &
      !$OMP PRIVATE(j, i) SCHEDULE(DYNAMIC)
      DO j = 1, n_grid
         DO i = j, n_grid
            factor(i, j) = factor(i, j)**2*jacobi_scale(i)*jacobi_scale(j)
         END DO
         factor(j, j) = factor(j, j) + tikhonov
      END DO
      !$OMP END PARALLEL DO
      CALL timeset(routineN//'_dpotrf', handle2)
      CALL dpotrf('L', n_grid, factor, n_grid, info)
      CALL timestop(handle2)
      CALL timestop(handle)
   END SUBROUTINE factor_ri_rs_fit

! **************************************************************************************************
!> \brief Evaluate the normalized cluster residual and its exact derivative in AO space.
!>
!>        V_l = sqrt(d_l) phi_l, S = V V^T, K = S∘S + λI, and K Z' = b, with
!>        b_lP = V_l^T B_P V_l. B_P is the symmetric matrix defined by the upper AO triangle.
!>        The objective is ||B - A Z'||²/||B||², A_:l = svec(V_l V_l^T), without a penalty term.
!>        For λ>0 the derivative includes the response of the regularized coefficients.
!>        Each RI-column panel reuses V B_P for the RHS and gradient before releasing the panel.
!>        The balancing convention and RI-RS reference are given in solve_ri_rs_fit.
!> \param phi unscaled AO values indexed (grid point, AO)
!> \param int_3c three-center integrals indexed (AO, AO, RI); only the upper AO triangle is used
!> \param tikhonov regularization in the Jacobi-scaled space
!> \param normalized_error normalized squared residual
!> \param fit_successful false if the target is empty or the regularized fit cannot be solved
!> \param ao_gradient optional derivative with respect to the unscaled AO values
!> \param max_abs_error optional largest absolute unweighted integral error
!> \param workspace_budget_bytes optional total budget for evaluator-owned gradient workspaces;
!>        the caller deducts its own arrays first. Reserve all fixed workspaces before sizing
!>        the product cache. An absent or nonpositive budget streams one RI column.
! **************************************************************************************************
   SUBROUTINE evaluate_ri_rs_fit(phi, int_3c, tikhonov, normalized_error, fit_successful, &
                                 ao_gradient, max_abs_error, workspace_budget_bytes)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: int_3c
      REAL(KIND=dp), INTENT(IN)                          :: tikhonov
      REAL(KIND=dp), INTENT(OUT)                         :: normalized_error
      LOGICAL, INTENT(OUT)                               :: fit_successful
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: ao_gradient
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: max_abs_error
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: workspace_budget_bytes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'evaluate_ri_rs_fit'
      INTEGER, PARAMETER                                 :: residual_block_size = 64

      INTEGER                                            :: first, handle, handle2, i, icache, imu, &
                                                            info, inu, ip, j, last, n_ao, &
                                                            n_cached, n_grid, n_panel, n_ri, &
                                                            panel_first, panel_last
      REAL(KIND=dp) :: coefficient_norm2, error_scale, fixed_bytes, integral_norm2, na, ng, &
         projection, residual_norm2, roundoff_estimate, row_norm2, value, weight
      REAL(KIND=dp), ALLOCATABLE :: coefficient_work(:, :), factor(:, :), jacobi_scale(:), &
         overlap(:, :), products(:, :, :), residual(:, :), v(:, :), z_prime(:, :)

      CALL timeset(routineN, handle)
      n_grid = SIZE(phi, 1)
      n_ao = SIZE(phi, 2)
      n_ri = SIZE(int_3c, 3)
      CPASSERT(SIZE(int_3c, 1) == n_ao .AND. SIZE(int_3c, 2) == n_ao)
      CPASSERT(tikhonov >= 0.0_dp)
      IF (PRESENT(ao_gradient)) THEN
         CPASSERT(ALL(SHAPE(ao_gradient) == SHAPE(phi)))
         ao_gradient(:, :) = 0.0_dp
      END IF
      normalized_error = HUGE(normalized_error)
      IF (PRESENT(max_abs_error)) max_abs_error = HUGE(max_abs_error)
      fit_successful = .FALSE.
      IF (n_grid < 1 .OR. n_ri < 1) THEN
         CALL timestop(handle)
         RETURN
      END IF
      integral_norm2 = 0.0_dp
      DO ip = 1, n_ri
         DO inu = 1, n_ao
            integral_norm2 = integral_norm2 + int_3c(inu, inu, ip)**2 + &
                             2.0_dp*DOT_PRODUCT(int_3c(1:inu - 1, inu, ip), int_3c(1:inu - 1, inu, ip))
         END DO
      END DO
      IF (integral_norm2 <= TINY(1.0_dp)) THEN
         CALL timestop(handle)
         RETURN
      END IF

      n_panel = n_ri
      n_cached = 1
      IF (PRESENT(ao_gradient)) THEN
         n_panel = 1
         IF (PRESENT(workspace_budget_bytes)) THEN
            IF (workspace_budget_bytes > 0.0_dp) THEN
               ng = REAL(n_grid, dp)
               na = REAL(n_ao, dp)
               ! Reserve factor, overlap, V, two coefficient arrays, Jacobi scales and the
               ! largest residual panel before assigning all remaining budget to V B_P.
               fixed_bytes = 8.0_dp*(2.0_dp*ng**2 + ng*na + 2.0_dp*ng*REAL(n_ri, dp) + ng + &
                                     na*MIN(na, REAL(residual_block_size, dp)))
               n_panel = INT(MIN(REAL(n_ri, dp), MAX(0.0_dp, &
                                                     (workspace_budget_bytes - fixed_bytes)/(8.0_dp*ng*na))))
               IF (n_panel < 1) THEN
                  CPABORT("RI-RS minimum fit workspace exceeds memory budget")
               END IF
            END IF
         END IF
         n_cached = n_panel
         CALL factor_ri_rs_fit(phi, tikhonov, factor, jacobi_scale, info, overlap)
      ELSE
         CALL factor_ri_rs_fit(phi, tikhonov, factor, jacobi_scale, info)
      END IF
      IF (info /= 0) THEN
         CALL timestop(handle)
         RETURN
      END IF
      ALLOCATE (v(n_grid, n_ao), products(n_grid, n_ao, n_cached), &
                coefficient_work(n_grid, n_ri), z_prime(n_grid, n_ri))
      DO imu = 1, n_ao
         v(:, imu) = SQRT(jacobi_scale)*phi(:, imu)
      END DO
      ! K is independent of P. Solve a panel of RHSs and consume its cached products immediately.
      ! Keep Z' and the coefficient response for the final grid-space symmetric product.
      value = 0.0_dp
      coefficient_norm2 = 0.0_dp
      CALL timeset(routineN//'_ri_panels', handle2)
      DO panel_first = 1, n_ri, n_panel
         panel_last = MIN(panel_first + n_panel - 1, n_ri)
         DO ip = panel_first, panel_last
            icache = 1
            IF (PRESENT(ao_gradient)) icache = ip - panel_first + 1
            CALL dsymm('R', 'U', n_grid, n_ao, 1.0_dp, int_3c(:, :, ip), n_ao, v, n_grid, &
                       0.0_dp, products(:, :, icache), n_grid)
            DO i = 1, n_grid
               coefficient_work(i, ip) = DOT_PRODUCT(v(i, :), products(i, :, icache))
            END DO
         END DO
         z_prime(:, panel_first:panel_last) = coefficient_work(:, panel_first:panel_last)
         CALL dpotrs('L', n_grid, panel_last - panel_first + 1, factor, n_grid, &
                     z_prime(:, panel_first:panel_last), n_grid, info)
         IF (info /= 0) THEN
            CALL timestop(handle2)
            CALL timestop(handle)
            RETURN
         END IF
         DO ip = panel_first, panel_last
            value = value + DOT_PRODUCT(z_prime(:, ip), coefficient_work(:, ip))
            coefficient_norm2 = coefficient_norm2 + DOT_PRODUCT(z_prime(:, ip), z_prime(:, ip))
         END DO

         IF (PRESENT(ao_gradient)) THEN
            ! Reuse the RHS for Y = K^-1 Z' and then C = Z' + λY; λ=0 needs no response solve.
            coefficient_work(:, panel_first:panel_last) = z_prime(:, panel_first:panel_last)
            IF (tikhonov > 0.0_dp) THEN
               CALL dpotrs('L', n_grid, panel_last - panel_first + 1, factor, n_grid, &
                           coefficient_work(:, panel_first:panel_last), n_grid, info)
               IF (info /= 0) THEN
                  CALL timestop(handle2)
                  CALL timestop(handle)
                  RETURN
               END IF
               coefficient_work(:, panel_first:panel_last) = z_prime(:, panel_first:panel_last) + &
                                                             tikhonov*coefficient_work(:, panel_first:panel_last)
            END IF
            ! Accumulate -Σ_P C_lP (V B_P)_l while the panel products are still available.
            DO ip = panel_first, panel_last
               icache = ip - panel_first + 1
               DO imu = 1, n_ao
                  ao_gradient(:, imu) = ao_gradient(:, imu) - coefficient_work(:, ip)*products(:, imu, icache)
               END DO
            END DO
            IF (tikhonov > 0.0_dp) THEN
               ! Retain T = C - Z'/2 = Z'/2 + λY for M = Z'T^T + TZ'^T below.
               coefficient_work(:, panel_first:panel_last) = coefficient_work(:, panel_first:panel_last) - &
                                                             0.5_dp*z_prime(:, panel_first:panel_last)
            END IF
         END IF
      END DO
      CALL timestop(handle2)
      ! K Z' = b implies ||B-AZ'||² = ||B||² - tr(Z'^T b) - λ ||Z'||².
      weight = tikhonov*coefficient_norm2
      residual_norm2 = integral_norm2 - value - weight
      error_scale = integral_norm2 + ABS(value) + weight
      ! Large coefficients can also amplify the factorization's backward error. The scaled
      ! Gram has diagonal <= 1, hence ||K||_2 <= n_grid + λ; use a conservative roundoff estimate.
      roundoff_estimate = 8.0_dp*EPSILON(1.0_dp)*REAL(n_grid, dp)* &
                          (REAL(n_grid, dp) + tikhonov)*coefficient_norm2
      ! Near an exact fit, use bounded AO-column panels to avoid cancellation. This path also
      ! supplies the final maximum-error diagnostic, without storing a pair-by-RI residual.
      IF (PRESENT(max_abs_error) .OR. &
          residual_norm2 <= MAX(SQRT(EPSILON(1.0_dp))*error_scale, roundoff_estimate)) THEN
         CALL timeset(routineN//'_residual', handle2)
         ALLOCATE (residual(n_ao, MIN(residual_block_size, n_ao)))
         residual_norm2 = 0.0_dp
         IF (PRESENT(max_abs_error)) max_abs_error = 0.0_dp
         DO ip = 1, n_ri
            DO first = 1, n_ao, residual_block_size
               last = MIN(first + residual_block_size - 1, n_ao)
               DO inu = first, last
                  products(:, inu - first + 1, 1) = z_prime(:, ip)*v(:, inu)
               END DO
               CALL dgemm('T', 'N', n_ao, last - first + 1, n_grid, -1.0_dp, v, n_grid, products(:, :, 1), n_grid, &
                          0.0_dp, residual, n_ao)
               DO inu = first, last
                  DO imu = 1, inu
                     value = int_3c(imu, inu, ip) + residual(imu, inu - first + 1)
                     weight = REAL(2 - MERGE(1, 0, imu == inu), dp)
                     residual_norm2 = residual_norm2 + weight*value**2
                     IF (PRESENT(max_abs_error)) max_abs_error = MAX(max_abs_error, ABS(value))
                  END DO
               END DO
            END DO
         END DO
         CALL timestop(handle2)
      END IF
      DEALLOCATE (products)
      normalized_error = residual_norm2/integral_norm2

      IF (PRESENT(ao_gradient)) THEN
         CALL timeset(routineN//'_gradient', handle2)
         ! All panel solves are complete. Reuse the factor for M = Z'Z'^T + λ(Z'Y^T + YZ'^T).
         IF (tikhonov > 0.0_dp) THEN
            CALL dsyr2k('L', 'N', n_grid, n_ri, 1.0_dp, z_prime, n_grid, coefficient_work, n_grid, &
                        0.0_dp, factor, n_grid)
         ELSE
            CALL dsyrk('L', 'N', n_grid, n_ri, 1.0_dp, z_prime, n_grid, 0.0_dp, factor, n_grid)
         END IF
         DO j = 1, n_grid
            DO i = j, n_grid
               factor(i, j) = factor(i, j)*overlap(i, j)
            END DO
         END DO
         ! Add (M∘S)V to the accumulated negative term, then normalize the AO gradient.
         CALL dsymm('L', 'L', n_grid, n_ao, 1.0_dp, factor, n_grid, v, n_grid, &
                    1.0_dp, ao_gradient, n_grid)
         ao_gradient(:, :) = (4.0_dp/integral_norm2)*ao_gradient
         ! V_l = sqrt(d_l) phi_l. Differentiate the normalization only above the floor;
         ! on the clamped branch d_l is constant, including at the branch boundary.
         DO i = 1, n_grid
            row_norm2 = DOT_PRODUCT(phi(i, :), phi(i, :))
            IF (row_norm2 > SQRT(jacobi_floor)) THEN
               projection = DOT_PRODUCT(phi(i, :), ao_gradient(i, :))/row_norm2
               ao_gradient(i, :) = ao_gradient(i, :) - projection*phi(i, :)
            END IF
            ao_gradient(i, :) = SQRT(jacobi_scale(i))*ao_gradient(i, :)
         END DO
         CALL timestop(handle2)
      END IF
      fit_successful = .TRUE.
      CALL timestop(handle)
   END SUBROUTINE evaluate_ri_rs_fit

! **************************************************************************************************
!> \brief Multiplies every matrix row by the corresponding diagonal entry:
!>        A(l, :) <- d_l A(l, :).
!> \param matrix ...
!> \param diagonal ...
! **************************************************************************************************
   SUBROUTINE scale_rows_by_diag(matrix, diagonal)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: diagonal

      CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_rows_by_diag'

      INTEGER                                            :: handle, i, j

      CALL timeset(routineN, handle)
      CPASSERT(SIZE(diagonal) == SIZE(matrix, 1))

      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(matrix, diagonal) &
      !$OMP PRIVATE(j, i) &
      !$OMP SCHEDULE(STATIC)
      DO j = 1, SIZE(matrix, 2)
         DO i = 1, SIZE(matrix, 1)
            matrix(i, j) = matrix(i, j)*diagonal(i)
         END DO
      END DO
      !$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE scale_rows_by_diag

! **************************************************************************************************
!> \brief Stores dense Z_lP columns in the distributed block-sparse matrix.
!> \param mat_Z_lP ...
!> \param Z_local ...
!> \param local_grid_idx ...
!> \param n_local_grid ...
!> \param n_loc_ri ...
!> \param atom_P ...
!> \param r_blk_sizes ...
!> \param row_offset ...
!> \param eps_filter ...
! **************************************************************************************************
   SUBROUTINE store_Z_lP_columns(mat_Z_lP, Z_local, local_grid_idx, n_local_grid, n_loc_ri, &
                                 atom_P, r_blk_sizes, row_offset, eps_filter)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_Z_lP
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Z_local
      INTEGER, DIMENSION(:), INTENT(IN)                  :: local_grid_idx
      INTEGER, INTENT(IN)                                :: n_local_grid, n_loc_ri, atom_P
      INTEGER, DIMENSION(:), INTENT(IN)                  :: r_blk_sizes, row_offset
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'store_Z_lP_columns'

      INTEGER                                            :: current_chunk_size, g_pt, handle, i_blk, &
                                                            loc_ptr, r_end, r_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Z_blk

      CALL timeset(routineN, handle)

      ALLOCATE (Z_blk(MAXVAL(r_blk_sizes), n_loc_ri))
      ! local_grid_idx is ascending, so one walking pointer maps sphere rows to global grid blocks.
      loc_ptr = 1

      DO i_blk = 1, SIZE(r_blk_sizes)
         r_start = row_offset(i_blk) + 1
         r_end = row_offset(i_blk) + r_blk_sizes(i_blk)
         current_chunk_size = r_blk_sizes(i_blk)
         Z_blk = 0.0_dp

         DO WHILE (loc_ptr <= n_local_grid)
            g_pt = local_grid_idx(loc_ptr)
            IF (g_pt > r_end) EXIT
            Z_blk(g_pt - r_start + 1, 1:n_loc_ri) = Z_local(loc_ptr, 1:n_loc_ri)
            loc_ptr = loc_ptr + 1
         END DO

         ! DBCSR routes each retained block to its global owner when the caller finalizes mat_Z_lP.
         IF (MAXVAL(ABS(Z_blk(1:current_chunk_size, 1:n_loc_ri))) > eps_filter) THEN
            CALL dbcsr_put_block(mat_Z_lP, row=i_blk, col=atom_P, &
                                 block=Z_blk(1:current_chunk_size, 1:n_loc_ri))
         END IF
      END DO

      DEALLOCATE (Z_blk)

      CALL timestop(handle)

   END SUBROUTINE store_Z_lP_columns

! **************************************************************************************************
!> \brief Solves the Jacobi-scaled, regularized RI-RS system with ScaLAPACK.
!>
!>        Each rank builds its block-cyclic part of
!>
!>          D'_ll' = d_l [Σ_μ ϕ_μ(r_l) ϕ_μ(r_l')]² d_l' + λ δ_ll'.
!>
!>        The Cholesky solution Z'_lP is gathered into the replicated right-hand side.
!>        The caller then restores Z_lP = d_l Z'_lP.
!> \param phi_local replicated AO values on the fitting grid
!> \param d_vec Jacobi scaling vector d_l, with the same floor as in the dense solve
!> \param d_lp replicated scaled RHS d_l d_lP on entry, intermediate solution Z'_lP on return
!> \param tikhonov ...
!> \param para_env_sub ...
!> \param blacs_env_sub ...
! **************************************************************************************************
   SUBROUTINE solve_D_lp_distributed(phi_local, d_vec, d_lp, tikhonov, para_env_sub, blacs_env_sub)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: phi_local
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: d_vec
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: d_lp
      REAL(KIND=dp), INTENT(IN)                          :: tikhonov
      TYPE(mp_para_env_type), POINTER                    :: para_env_sub
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_D_lp_distributed'

      INTEGER                                            :: handle, i_loc, ig, info, j_loc, jg, &
                                                            n_ao, n_loc, n_rhs, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_b, fm_struct_D
      TYPE(cp_fm_type)                                   :: fm_b, fm_D

      CALL timeset(routineN, handle)
      n_loc = SIZE(phi_local, 1)
      n_ao = SIZE(phi_local, 2)
      n_rhs = SIZE(d_lp, 2)

      NULLIFY (fm_struct_D, fm_struct_b)
      CALL cp_fm_struct_create(fm_struct_D, para_env=para_env_sub, &
                               context=blacs_env_sub, &
                               nrow_global=n_loc, ncol_global=n_loc)
      CALL cp_fm_struct_create(fm_struct_b, para_env=para_env_sub, &
                               context=blacs_env_sub, &
                               nrow_global=n_loc, ncol_global=n_rhs)
      CALL cp_fm_create(fm_D, fm_struct_D)
      CALL cp_fm_create(fm_b, fm_struct_b)

      CALL cp_fm_get_info(fm_D, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, &
                          local_data=local_data)

      ! Build only the rank-owned Gram tiles, avoiding a replicated n_grid × n_grid matrix.
      IF (nrow_local > 0 .AND. ncol_local > 0) THEN
         BLOCK
            INTEGER, PARAMETER :: ntile = 1024
            INTEGER :: ib, ie, jb, je, mb, kb, ti, tj, handle_dgemm
            REAL(KIND=dp), ALLOCATABLE :: gram_t(:, :), phi_cols_t(:, :), phi_rows_t(:, :)
            ALLOCATE (phi_rows_t(ntile, n_ao), phi_cols_t(n_ao, ntile), gram_t(ntile, ntile))
            ! Pack AO values for the grid rows/columns selected by the block-cyclic FM layout.
            DO ib = 1, nrow_local, ntile
               ie = MIN(ib + ntile - 1, nrow_local)
               mb = ie - ib + 1
               !$OMP PARALLEL DO DEFAULT(NONE) &
               !$OMP SHARED(mb, n_ao, phi_rows_t, phi_local, row_indices, ib) &
               !$OMP PRIVATE(ti, j_loc) SCHEDULE(STATIC)
               DO j_loc = 1, n_ao
                  DO ti = 1, mb
                     phi_rows_t(ti, j_loc) = phi_local(row_indices(ib + ti - 1), j_loc)
                  END DO
               END DO
               !$OMP END PARALLEL DO
               DO jb = 1, ncol_local, ntile
                  je = MIN(jb + ntile - 1, ncol_local)
                  kb = je - jb + 1
                  !$OMP PARALLEL DO DEFAULT(NONE) &
                  !$OMP SHARED(kb, n_ao, phi_cols_t, phi_local, col_indices, jb) &
                  !$OMP PRIVATE(tj, i_loc) SCHEDULE(STATIC)
                  DO tj = 1, kb
                     DO i_loc = 1, n_ao
                        phi_cols_t(i_loc, tj) = phi_local(col_indices(jb + tj - 1), i_loc)
                     END DO
                  END DO
                  !$OMP END PARALLEL DO
                  CALL timeset(routineN//'_dgemm', handle_dgemm)
                  CALL dgemm('N', 'N', mb, kb, n_ao, &
                             1.0_dp, phi_rows_t, ntile, phi_cols_t, n_ao, &
                             0.0_dp, gram_t, ntile)
                  CALL timestop(handle_dgemm)
                  ! Square the AO Gram entries, apply d_l d_l', and add λ on the global diagonal.
                  !$OMP PARALLEL DO DEFAULT(NONE) &
                  !$OMP SHARED(mb, kb, gram_t, d_vec, row_indices, col_indices, ib, jb) &
                  !$OMP SHARED(local_data, tikhonov) &
                  !$OMP PRIVATE(ti, tj, ig, jg) SCHEDULE(STATIC)
                  DO tj = 1, kb
                     jg = col_indices(jb + tj - 1)
                     DO ti = 1, mb
                        ig = row_indices(ib + ti - 1)
                        local_data(ib + ti - 1, jb + tj - 1) = &
                           gram_t(ti, tj)*gram_t(ti, tj)*d_vec(ig)*d_vec(jg)
                        IF (ig == jg) THEN
                           local_data(ib + ti - 1, jb + tj - 1) = &
                              local_data(ib + ti - 1, jb + tj - 1) + tikhonov
                        END IF
                     END DO
                  END DO
                  !$OMP END PARALLEL DO
               END DO
            END DO
            DEALLOCATE (phi_rows_t, phi_cols_t, gram_t)
         END BLOCK
      END IF

      ! Distribute the pre-scaled RHS and solve for Z'_lP within the BLACS subgroup.
      CALL cp_fm_set_submatrix(fm_b, d_lp)

      CALL cp_fm_cholesky_decompose(fm_D, n=n_loc, info_out=info)
      IF (info /= 0) CPABORT("pdpotrf failed in solve_D_lp_distributed")

      CALL cp_fm_cholesky_solve(fm_D, fm_b, n=n_loc, info_out=info)
      IF (info /= 0) CPABORT("pdpotrs failed in solve_D_lp_distributed")

      ! Replicate Z'_lP within the subgroup; the caller applies the final Jacobi rescaling.
      CALL cp_fm_get_submatrix(fm_b, d_lp)

      CALL cp_fm_release(fm_D)
      CALL cp_fm_release(fm_b)
      CALL cp_fm_struct_release(fm_struct_D)
      CALL cp_fm_struct_release(fm_struct_b)

      CALL timestop(handle)

   END SUBROUTINE solve_D_lp_distributed

END MODULE ri_rs_fitting

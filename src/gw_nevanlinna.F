!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Nevanlinna analytic continuation of the GW correlation self-energy from the imaginary
!>        to the real frequency axis, following J. Fei, C.-N. Yeh, E. Gull, PRL 126, 056402 (2021).
!>        Minus the correlation self-energy maps the upper half plane into itself. After a Moebius
!>        map onto the unit disk the interpolation becomes a Schur problem whose solution is
!>        causal by construction. The Schur recursion loses digits quickly and is therefore done
!>        in double-double arithmetic (about 31 significant digits, module double_double), so that
!>        no quadruple-precision support of the compiler is needed; the evaluation is stable in
!>        double precision.
!> \par History
!>      11.2026 created
! **************************************************************************************************
MODULE gw_nevanlinna
   USE double_double,                   ONLY: &
        zdd, zdd_abs2, zdd_add, zdd_conjg, zdd_div, zdd_i, zdd_mul, zdd_one, zdd_sub, zdd_to_dp, &
        zdd_type, zdd_zero
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathlib,                         ONLY: diag_complex
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_nevanlinna'

   ! imaginary shifts delta*max|Sigma| tried in turn until the samples admit an interpolant
   REAL(KIND=dp), DIMENSION(9), PARAMETER :: nevanlinna_deltas = [0.0_dp, 1.0E-14_dp, 1.0E-12_dp, &
                                                                  1.0E-10_dp, 1.0E-8_dp, 1.0E-6_dp, &
                                                                  1.0E-4_dp, 1.0E-3_dp, 1.0E-2_dp]
   REAL(KIND=dp), PARAMETER :: tiny_schur = 1.0E-40_dp

   PUBLIC :: gw_nevanlinna_env_type, gw_qp_summary_type, nevanlinna_type, nevanlinna_create, &
             nevanlinna_evaluate, nevanlinna_evaluate_real_grid, nevanlinna_release
   PUBLIC :: nevanlinna_window_grid_type, nevanlinna_prepare_window_grid, nevanlinna_create_windows

! **************************************************************************************************
!> \brief Parameters of the quasiparticle solution with the Nevanlinna continuation; they are set
!>        in set_heuristic_parameters (gw_utils.F) and stay zero for the Padé continuation
! **************************************************************************************************
   TYPE gw_nevanlinna_env_type
      ! relative SVD cutoff of the least-squares minimax weights of Σ^c(iτ) -> Σ^c(iω); the
      ! truncated SVD keeps Σ^c(iω) causal, which the Nevanlinna interpolant requires
      REAL(KIND=dp)                                      :: rcond_weights = 0.0_dp
      ! fit windows: the lowest k >= window_min_points imaginary frequencies, with the highest one
      ! in [window_freq_min, window_freq_max]; the quasiparticle energy is the median over windows
      INTEGER                                            :: window_min_points = 0
      REAL(KIND=dp)                                      :: window_freq_min = 0.0_dp, &
                                                            window_freq_max = 0.0_dp
      ! a level is marked unstable if its window solutions spread more than window_spread
      REAL(KIND=dp)                                      :: window_spread = 0.0_dp
      ! Newton search from the DFT eigenvalue with steps of at most qp_newton_max_step, then the
      ! solution of largest weight Z within qp_maxz_radius of the Newton solution
      REAL(KIND=dp)                                      :: qp_newton_max_step = 0.0_dp, &
                                                            qp_maxz_radius = 0.0_dp
   END TYPE gw_nevanlinna_env_type

! **************************************************************************************************
!> \brief Schur parameters of one Nevanlinna interpolant
! **************************************************************************************************
   TYPE nevanlinna_type
      LOGICAL                                            :: solvable = .FALSE.
      LOGICAL                                            :: mirrored = .FALSE.
      INTEGER                                            :: n = 0
      COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE        :: xis, phis
   END TYPE nevanlinna_type

! **************************************************************************************************
!> \brief Frequency-dependent Schur factors shared by the bands and fit windows of a k-point.
! **************************************************************************************************
   TYPE nevanlinna_window_grid_type
      PRIVATE
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: omega
      TYPE(zdd_type), ALLOCATABLE, DIMENSION(:)          :: xis
      TYPE(zdd_type), ALLOCATABLE, DIMENSION(:, :)       :: factors
   END TYPE nevanlinna_window_grid_type

! **************************************************************************************************
!> \brief Nevanlinna diagnostics for one GW cycle. n_kpoints is the number of target k-points per spin.
! **************************************************************************************************
   TYPE gw_qp_summary_type
      INTEGER                                            :: n_kpoints = 1
      INTEGER(KIND=int_8)                                :: n_reliable = 0_int_8, n_unstable = 0_int_8, &
                                                            n_multi = 0_int_8, n_pade = 0_int_8
      REAL(KIND=dp)                                      :: start_time = 0.0_dp
   END TYPE gw_qp_summary_type

CONTAINS

! **************************************************************************************************
!> \brief Builds the Nevanlinna interpolant of y(x). If the samples are not causal within the
!>        accuracy of the data, a small imaginary shift -i*delta*max|y| is added, with delta taken
!>        from a fixed ladder; nev%solvable is .FALSE. if no delta of the ladder works.
!> \param nev the interpolant
!> \param x sample points, all in the upper or all in the lower half plane
!> \param y values of the correlation self-energy at x
!> \param pick_min optional smallest relative eigenvalue of the unshifted Pick matrix
! **************************************************************************************************
   SUBROUTINE nevanlinna_create(nev, x, y, pick_min)
      TYPE(nevanlinna_type), INTENT(OUT)                 :: nev
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: x, y
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: pick_min

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nevanlinna_create'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: xx, yy
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: pick, pick_vec
      INTEGER                                            :: handle, i, i_delta, k, n, n_valid
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: order
      REAL(KIND=dp)                                      :: scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: pick_eig
      TYPE(zdd_type), ALLOCATABLE, DIMENSION(:)          :: phis, xis

      CALL timeset(routineN, handle)

      IF (PRESENT(pick_min)) pick_min = 0.0_dp
      n = SIZE(x)
      nev%n = n
      ALLOCATE (nev%xis(n), nev%phis(n), xis(n), phis(n), xx(n), yy(n), order(n))
      nev%xis(:) = (0.0_dp, 0.0_dp)
      nev%phis(:) = (0.0_dp, 0.0_dp)

      ! occupied levels are sampled in the lower half plane: continue the mirror image
      nev%mirrored = ALL(AIMAG(x) < 0.0_dp)
      IF (nev%mirrored) THEN
         xx(:) = CONJG(x)
         yy(:) = CONJG(y)
      ELSE
         xx(:) = x
         yy(:) = y
      END IF

      IF (n == 0 .OR. .NOT. (ALL(ABS(xx) < HUGE(1.0_dp)) .AND. ALL(ABS(yy) < HUGE(1.0_dp)))) THEN
         DEALLOCATE (xis, phis, xx, yy, order)
         CALL timestop(handle)
         RETURN
      END IF

      CALL sort_by_imag(xx, order)

      scale = MAXVAL(ABS(yy))
      IF (scale == 0.0_dp) scale = 1.0_dp

      DO i = 1, n
         xis(i) = moebius(zdd(xx(order(i))))
      END DO

      DO i_delta = 1, SIZE(nevanlinna_deltas)
         DO i = 1, n
            phis(i) = shifted_sample(yy(order(i)), nevanlinna_deltas(i_delta)*scale)
         END DO

         IF (i_delta == 1 .AND. PRESENT(pick_min)) THEN
            ! smallest eigenvalue of the Pick matrix of the unshifted samples, relative to the
            ! largest one; negative values quantify how far the data are from causal
            ALLOCATE (pick(n, n), pick_vec(n, n), pick_eig(n))
            DO k = 1, n
               DO i = 1, n
                  pick(i, k) = zdd_to_dp(zdd_div(zdd_sub(zdd_one, zdd_mul(phis(i), zdd_conjg(phis(k)))), &
                                                 zdd_sub(zdd_one, zdd_mul(xis(i), zdd_conjg(xis(k))))))
               END DO
            END DO
            pick(:, :) = 0.5_dp*(pick + CONJG(TRANSPOSE(pick)))
            CALL diag_complex(pick, pick_vec, pick_eig)
            pick_min = MINVAL(pick_eig)/MAX(MAXVAL(ABS(pick_eig)), TINY(1.0_dp))
            DEALLOCATE (pick, pick_vec, pick_eig)
         END IF

         n_valid = schur_prefix(xis, phis)
         nev%solvable = n_valid == n
         IF (nev%solvable) EXIT
      END DO

      DO i = 1, n
         nev%xis(i) = zdd_to_dp(xis(i))
         IF (nev%solvable) nev%phis(i) = zdd_to_dp(phis(i))
      END DO

      DEALLOCATE (xis, phis, xx, yy, order)

      CALL timestop(handle)

   END SUBROUTINE nevanlinna_create

! **************************************************************************************************
!> \brief Prepare the Moebius nodes and triangular table of Blaschke factors for a prefix grid.
!> \param grid shared frequency data
!> \param omega imaginary frequencies in window order
! **************************************************************************************************
   SUBROUTINE nevanlinna_prepare_window_grid(grid, omega)
      TYPE(nevanlinna_window_grid_type), INTENT(OUT)     :: grid
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: omega

      INTEGER                                            :: j, k, n

      n = SIZE(omega)
      ALLOCATE (grid%omega(n))
      grid%omega(:) = ABS(omega)
      ! With a zero frequency the occupied-state mirror convention changes. An unordered
      ! grid can also give different sorted prefixes. Use the single-window path in either case.
      IF (.NOT. ALL(grid%omega > 0.0_dp .AND. grid%omega < HUGE(1.0_dp))) RETURN
      IF (ANY(grid%omega(2:n) < grid%omega(1:n - 1))) RETURN
      ALLOCATE (grid%xis(n))
      grid%xis(:) = moebius(zdd(CMPLX(0.0_dp, grid%omega, KIND=dp)))
      IF (ANY(zdd_abs2(grid%xis) >= 1.0_dp)) RETURN
      ALLOCATE (grid%factors(n, n))
      DO k = 1, n
         DO j = 1, k - 1
            grid%factors(j, k) = blaschke_factor(grid%xis(j), grid%xis(k))
         END DO
      END DO
   END SUBROUTINE nevanlinna_prepare_window_grid

! **************************************************************************************************
!> \brief Fit nested windows for a single band, reusing Schur parameters at equal absolute shifts.
!> \param fits one interpolant per window, with solvable false if all shifts failed
!> \param grid frequency data shared across bands
!> \param y self-energy samples in the same order as the grid
!> \param occupied whether the samples were taken on the negative imaginary axis
!> \param sizes positive window lengths in increasing order
! **************************************************************************************************
   SUBROUTINE nevanlinna_create_windows(fits, grid, y, occupied, sizes)
      TYPE(nevanlinna_type), DIMENSION(:), INTENT(OUT)   :: fits
      TYPE(nevanlinna_window_grid_type), INTENT(IN)      :: grid
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: y
      LOGICAL, INTENT(IN)                                :: occupied
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'nevanlinna_create_windows'
      COMPLEX(KIND=dp), PARAMETER                        :: ii = (0.0_dp, 1.0_dp)

      COMPLEX(KIND=dp), DIMENSION(SIZE(y))               :: yy
      INTEGER                                            :: handle, i, i_delta, n, n_valid
      REAL(KIND=dp)                                      :: previous_shift, shift, sign_freq
      REAL(KIND=dp), DIMENSION(SIZE(sizes))              :: window_scale
      REAL(KIND=dp), DIMENSION(SIZE(y))                  :: abs_y
      TYPE(zdd_type), DIMENSION(SIZE(y))                 :: phis

      CALL timeset(routineN, handle)
      abs_y(:) = ABS(y)
      IF (.NOT. ALLOCATED(grid%factors) .OR. .NOT. ALL(abs_y < HUGE(1.0_dp))) THEN
         ! A bad sample in a longer window must not invalidate its shorter prefixes.
         sign_freq = MERGE(-1.0_dp, 1.0_dp, occupied)
         DO i = 1, SIZE(sizes)
            n = sizes(i)
            CALL nevanlinna_create(fits(i), ii*(sign_freq*grid%omega(1:n)), y(1:n))
         END DO
         CALL timestop(handle)
         RETURN
      END IF

      yy(:) = y
      IF (occupied) yy(:) = CONJG(y)
      DO i = 1, SIZE(sizes)
         n = sizes(i)
         fits(i)%n = n
         fits(i)%mirrored = occupied
         window_scale(i) = MAXVAL(abs_y(1:n))
         IF (window_scale(i) == 0.0_dp) window_scale(i) = 1.0_dp
      END DO

      DO i_delta = 1, SIZE(nevanlinna_deltas)
         previous_shift = -1.0_dp
         n_valid = 0
         DO i = SIZE(sizes), 1, -1
            IF (fits(i)%solvable) CYCLE
            n = sizes(i)
            shift = nevanlinna_deltas(i_delta)*window_scale(i)
            ! A shorter window at this exact shift has the same parameters up to the
            ! first failed Schur step, even if the longer window needs another shift.
            IF (shift /= previous_shift) THEN
               phis(1:n) = shifted_sample(yy(1:n), shift)
               n_valid = schur_prefix(grid%xis(1:n), phis(1:n), grid%factors)
               previous_shift = shift
            END IF
            IF (n > n_valid) CYCLE
            fits(i)%solvable = .TRUE.
            ALLOCATE (fits(i)%xis(n), fits(i)%phis(n))
            fits(i)%xis(:) = zdd_to_dp(grid%xis(1:n))
            fits(i)%phis(:) = zdd_to_dp(phis(1:n))
         END DO
         IF (ALL(fits%solvable)) EXIT
      END DO
      CALL timestop(handle)
   END SUBROUTINE nevanlinna_create_windows

! **************************************************************************************************
!> \brief Evaluates the interpolant on the real axis
!> \param nev the interpolant
!> \param x_val real frequency
!> \param func_val continued correlation self-energy at x_val
! **************************************************************************************************
   SUBROUTINE nevanlinna_evaluate(nev, x_val, func_val)
      TYPE(nevanlinna_type), INTENT(IN)                  :: nev
      REAL(KIND=dp), INTENT(IN)                          :: x_val
      COMPLEX(KIND=dp), INTENT(OUT)                      :: func_val

      COMPLEX(KIND=dp), PARAMETER                        :: ii = (0.0_dp, 1.0_dp), &
                                                            one = (1.0_dp, 0.0_dp), &
                                                            zero = (0.0_dp, 0.0_dp)

      COMPLEX(KIND=dp)                                   :: aa, bb, bj, cc, cphi, dd, na, nb, nc, &
                                                            nd, phi, th, xi, xz
      INTEGER                                            :: j
      REAL(KIND=dp)                                      :: sc

      ! product of the Moebius transformations of the Schur recursion, applied to theta = 0
      xz = (CMPLX(x_val, 0.0_dp, KIND=dp) - ii)/(CMPLX(x_val, 0.0_dp, KIND=dp) + ii)
      aa = one
      bb = zero
      cc = zero
      dd = one
      DO j = 1, nev%n
         xi = nev%xis(j)
         phi = nev%phis(j)
         bj = (xz - xi)/(one - CONJG(xi)*xz)
         cphi = CONJG(phi)
         na = aa*bj + bb*cphi*bj
         nb = aa*phi + bb
         nc = cc*bj + dd*cphi*bj
         nd = cc*phi + dd
         ! rescale to avoid overflow; only the ratio bb/dd is used
         sc = MAX(ABS(na), ABS(nb), ABS(nc), ABS(nd))
         IF (sc == 0.0_dp) sc = 1.0_dp
         aa = na/sc
         bb = nb/sc
         cc = nc/sc
         dd = nd/sc
      END DO
      th = bb/dd
      ! inverse Moebius map, and Sigma = -(Nevanlinna function)
      func_val = -ii*(one + th)/(one - th)
      IF (nev%mirrored) func_val = CONJG(func_val)

   END SUBROUTINE nevanlinna_evaluate

! **************************************************************************************************
!> \brief Evaluates Re Sigma on a grid of real frequencies
!> \param nev the interpolant
!> \param x real frequencies
!> \param values real part of the continued correlation self-energy at x
! **************************************************************************************************
   SUBROUTINE nevanlinna_evaluate_real_grid(nev, x, values)
      TYPE(nevanlinna_type), INTENT(IN)                  :: nev
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: x
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: values

      COMPLEX(KIND=dp), PARAMETER                        :: ii = (0.0_dp, 1.0_dp)
      INTEGER, PARAMETER                                 :: block_size = 64, rescale_interval = 8
      REAL(KIND=dp), PARAMETER :: scalar_tol = SQRT(EPSILON(1.0_dp))

      COMPLEX(KIND=dp)                                   :: alpha, beta, c, p, phi, q, sigma
      COMPLEX(KIND=dp), DIMENSION(block_size)            :: u, v
      INTEGER                                            :: first, i, j, n
      REAL(KIND=dp)                                      :: scale
      REAL(KIND=dp), DIMENSION(block_size)               :: inv_scale, x_scaled

      IF (nev%n == 0) THEN
         values(:) = 0.0_dp
         RETURN
      END IF
      ! Near a node on the unit circle the two forms of the Blaschke factor can lose
      ! different digits. Use the scalar arithmetic for these ill-conditioned fits.
      IF (ANY(ABS(nev%xis(1:nev%n - 1)) > 1.0_dp - scalar_tol)) THEN
         DO i = 1, SIZE(x)
            CALL nevanlinna_evaluate(nev, x(i), sigma)
            values(i) = REAL(sigma, KIND=dp)
         END DO
         RETURN
      END IF
      DO first = 1, SIZE(x), block_size
         n = MIN(block_size, SIZE(x) - first + 1)
         ! Represent x by (x_scaled,inv_scale) to keep the factors bounded at large x.
         DO i = 1, n
            scale = MAX(1.0_dp, ABS(x(first + i - 1)))
            x_scaled(i) = x(first + i - 1)/scale
            inv_scale(i) = 1.0_dp/scale
         END DO
         ! Apply the product to (0,1) from right to left. The last factor gives (phi_n,1).
         u(1:n) = nev%phis(nev%n)
         v(1:n) = 1.0_dp
         DO j = nev%n - 1, 1, -1
            alpha = 1.0_dp - nev%xis(j)
            beta = -ii*(1.0_dp + nev%xis(j))
            phi = nev%phis(j)
            ! At real x the Blaschke factor is c/conjg(c), c = (1-xi)*x - i*(1+xi).
            ! A common real scale cancels from u/v. Bound c before multiplying the vector.
            !$OMP SIMD PRIVATE(c,scale,p,q)
            DO i = 1, n
               c = alpha*x_scaled(i) + beta*inv_scale(i)
               scale = MAX(ABS(REAL(c, KIND=dp)), ABS(AIMAG(c)), TINY(1.0_dp))
               c = c/scale
               p = c*u(i)
               q = CONJG(c)*v(i)
               u(i) = p + phi*q
               v(i) = q + CONJG(phi)*p
            END DO
            ! The common scale drops out of the final ratio; bound its drift in long fits.
            IF (MOD(j, rescale_interval) == 0) THEN
               !$OMP SIMD PRIVATE(scale)
               DO i = 1, n
                  scale = MAX(ABS(REAL(u(i), KIND=dp)), ABS(AIMAG(u(i))), &
                              ABS(REAL(v(i), KIND=dp)), ABS(AIMAG(v(i))), TINY(1.0_dp))
                  u(i) = u(i)/scale
                  v(i) = v(i)/scale
               END DO
            END IF
         END DO
         DO i = 1, n
            ! Near a pole use the scalar evaluator's order of operations for the sign test.
            IF (ABS(v(i) - u(i)) <= scalar_tol*ABS(v(i))) THEN
               CALL nevanlinna_evaluate(nev, x(first + i - 1), sigma)
               values(first + i - 1) = REAL(sigma, KIND=dp)
            ELSE
               ! Re[-i*(v+u)/(v-u)] is unchanged by the occupied-state conjugation.
               values(first + i - 1) = AIMAG((v(i) + u(i))/(v(i) - u(i)))
            END IF
         END DO
      END DO

   END SUBROUTINE nevanlinna_evaluate_real_grid

! **************************************************************************************************
!> \brief Releases the interpolant
!> \param nev the interpolant
! **************************************************************************************************
   SUBROUTINE nevanlinna_release(nev)
      TYPE(nevanlinna_type), INTENT(INOUT)               :: nev

      IF (ALLOCATED(nev%xis)) DEALLOCATE (nev%xis)
      IF (ALLOCATED(nev%phis)) DEALLOCATE (nev%phis)
      nev%solvable = .FALSE.
      nev%n = 0

   END SUBROUTINE nevanlinna_release

! **************************************************************************************************
!> \brief Moebius map of the upper half plane onto the unit disk
!> \param w point in the upper half plane
!> \return (w - i)/(w + i)
! **************************************************************************************************
   ELEMENTAL FUNCTION moebius(w) RESULT(h)
      TYPE(zdd_type), INTENT(IN)                         :: w
      TYPE(zdd_type)                                     :: h

      h = zdd_div(zdd_sub(w, zdd_i), zdd_add(w, zdd_i))

   END FUNCTION moebius

! **************************************************************************************************
!> \brief Map the shifted -Sigma sample to the unit disk.
!> \param y self-energy sample in the upper half plane
!> \param shift absolute imaginary shift
!> \return Moebius image of the shifted sample
! **************************************************************************************************
   ELEMENTAL FUNCTION shifted_sample(y, shift) RESULT(lam)
      COMPLEX(KIND=dp), INTENT(IN)                       :: y
      REAL(KIND=dp), INTENT(IN)                          :: shift
      TYPE(zdd_type)                                     :: lam

      lam = moebius(zdd(-(y - CMPLX(0.0_dp, shift, KIND=dp))))
   END FUNCTION shifted_sample

! **************************************************************************************************
!> \brief Blaschke factor at a pair of Moebius nodes.
!> \param xi Schur node
!> \param x evaluation node
!> \return Blaschke factor
! **************************************************************************************************
   PURE FUNCTION blaschke_factor(xi, x) RESULT(bfac)
      TYPE(zdd_type), INTENT(IN)                         :: xi, x
      TYPE(zdd_type)                                     :: bfac

      bfac = zdd_div(zdd_sub(x, xi), zdd_sub(zdd_one, zdd_mul(zdd_conjg(xi), x)))
   END FUNCTION blaschke_factor

! **************************************************************************************************
!> \brief Schur recursion through the first sample outside the unit disk.
!> \param xis Moebius images of the sample points, sorted by imaginary part
!> \param phis Moebius images of the function values on input, Schur parameters of the valid prefix on output
!> \param factors optional precomputed Blaschke factors
!> \return length of the valid prefix
! **************************************************************************************************
   FUNCTION schur_prefix(xis, phis, factors) RESULT(n_valid)
      TYPE(zdd_type), DIMENSION(:), INTENT(IN)           :: xis
      TYPE(zdd_type), DIMENSION(:), INTENT(INOUT)        :: phis
      TYPE(zdd_type), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: factors
      INTEGER                                            :: n_valid

      INTEGER                                            :: j, k
      TYPE(zdd_type)                                     :: bfac, den, phi

      n_valid = 0
      ! Finish one node before starting the next so a rejected prefix never updates the tail.
      DO k = 1, SIZE(xis)
         phi = phis(k)
         DO j = 1, k - 1
            IF (PRESENT(factors)) THEN
               bfac = factors(j, k)
            ELSE
               bfac = blaschke_factor(xis(j), xis(k))
            END IF
            den = zdd_sub(zdd_one, zdd_mul(zdd_conjg(phis(j)), phi))
            IF (zdd_abs2(den) < tiny_schur**2 .OR. zdd_abs2(bfac) < tiny_schur**2) THEN
               phi = zdd_zero
            ELSE
               phi = zdd_div(zdd_div(zdd_sub(phi, phis(j)), den), bfac)
            END IF
         END DO
         IF (zdd_abs2(phi) >= 1.0_dp) EXIT
         phis(k) = phi
         n_valid = k
      END DO

   END FUNCTION schur_prefix

! **************************************************************************************************
!> \brief Insertion sort of sample points by increasing imaginary part
!> \param z sample points
!> \param order permutation that sorts z
! **************************************************************************************************
   SUBROUTINE sort_by_imag(z, order)
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: z
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: order

      INTEGER                                            :: i, j, tmp

      DO i = 1, SIZE(z)
         order(i) = i
      END DO
      DO i = 2, SIZE(z)
         tmp = order(i)
         j = i - 1
         DO WHILE (j >= 1)
            IF (AIMAG(z(order(j))) <= AIMAG(z(tmp))) EXIT
            order(j + 1) = order(j)
            j = j - 1
         END DO
         order(j + 1) = tmp
      END DO

   END SUBROUTINE sort_by_imag

END MODULE gw_nevanlinna

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
   USE bibliography,                    ONLY: Fei2021,&
                                              cite_reference
   USE double_double,                   ONLY: &
        zdd, zdd_abs2, zdd_add, zdd_conjg, zdd_div, zdd_i, zdd_mul, zdd_one, zdd_sub, zdd_to_dp, &
        zdd_type, zdd_zero
   USE kinds,                           ONLY: dp
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

   PUBLIC :: gw_nevanlinna_env_type, nevanlinna_type, nevanlinna_create, nevanlinna_evaluate, &
             nevanlinna_release

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
      REAL(KIND=dp)                                      :: delta = 0.0_dp
      REAL(KIND=dp)                                      :: pick_min = 0.0_dp
      COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE        :: xis, phis
   END TYPE nevanlinna_type

CONTAINS

! **************************************************************************************************
!> \brief Builds the Nevanlinna interpolant of y(x). If the samples are not causal within the
!>        accuracy of the data, a small imaginary shift -i*delta*max|y| is added, with delta taken
!>        from a fixed ladder; nev%solvable is .FALSE. if no delta of the ladder works.
!> \param nev the interpolant
!> \param x sample points, all in the upper or all in the lower half plane
!> \param y values of the correlation self-energy at x
! **************************************************************************************************
   SUBROUTINE nevanlinna_create(nev, x, y)
      TYPE(nevanlinna_type), INTENT(OUT)                 :: nev
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: x, y

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nevanlinna_create'

      COMPLEX(KIND=dp)                                   :: s
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: xx, yy
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: pick, pick_vec
      INTEGER                                            :: handle, i, i_delta, k, n
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: order
      REAL(KIND=dp)                                      :: scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: pick_eig
      TYPE(zdd_type), ALLOCATABLE, DIMENSION(:)          :: lams, phis, xis

      CALL timeset(routineN, handle)

      CALL cite_reference(Fei2021)

      n = SIZE(x)
      nev%n = n
      ALLOCATE (nev%xis(n), nev%phis(n), xis(n), phis(n), lams(n), xx(n), yy(n), order(n))
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
         DEALLOCATE (xis, phis, lams, xx, yy, order)
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
            s = yy(order(i)) - CMPLX(0.0_dp, nevanlinna_deltas(i_delta)*scale, KIND=dp)
            ! -Sigma is the function that maps the upper half plane into itself
            lams(i) = moebius(zdd(-s))
         END DO

         IF (i_delta == 1) THEN
            ! smallest eigenvalue of the Pick matrix of the unshifted samples, relative to the
            ! largest one; negative values quantify how far the data are from causal
            ALLOCATE (pick(n, n), pick_vec(n, n), pick_eig(n))
            DO k = 1, n
               DO i = 1, n
                  pick(i, k) = zdd_to_dp(zdd_div(zdd_sub(zdd_one, zdd_mul(lams(i), zdd_conjg(lams(k)))), &
                                                 zdd_sub(zdd_one, zdd_mul(xis(i), zdd_conjg(xis(k))))))
               END DO
            END DO
            pick(:, :) = 0.5_dp*(pick + CONJG(TRANSPOSE(pick)))
            CALL diag_complex(pick, pick_vec, pick_eig)
            nev%pick_min = MINVAL(pick_eig)/MAX(MAXVAL(ABS(pick_eig)), TINY(1.0_dp))
            DEALLOCATE (pick, pick_vec, pick_eig)
         END IF

         IF (schur_parameters(xis, lams, phis)) THEN
            nev%solvable = .TRUE.
            nev%delta = nevanlinna_deltas(i_delta)
            EXIT
         END IF
      END DO

      DO i = 1, n
         nev%xis(i) = zdd_to_dp(xis(i))
         IF (nev%solvable) nev%phis(i) = zdd_to_dp(phis(i))
      END DO

      DEALLOCATE (xis, phis, lams, xx, yy, order)

      CALL timestop(handle)

   END SUBROUTINE nevanlinna_create

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
!> \brief Schur recursion; fails if an intermediate value leaves the unit disk, i.e. if no
!>        contractive interpolant of the samples exists
!> \param xis Moebius images of the sample points, sorted by imaginary part
!> \param lams Moebius images of the function values
!> \param phis Schur parameters
!> \return .TRUE. if the recursion stayed inside the unit disk
! **************************************************************************************************
   FUNCTION schur_parameters(xis, lams, phis) RESULT(ok)
      TYPE(zdd_type), DIMENSION(:), INTENT(IN)           :: xis, lams
      TYPE(zdd_type), DIMENSION(:), INTENT(OUT)          :: phis
      LOGICAL                                            :: ok

      INTEGER                                            :: j, k
      TYPE(zdd_type)                                     :: bfac, den, phi
      TYPE(zdd_type), DIMENSION(SIZE(xis))               :: vals

      vals(:) = lams
      phis(:) = zdd_zero
      ok = .TRUE.
      DO j = 1, SIZE(xis)
         phi = vals(j)
         IF (zdd_abs2(phi) >= 1.0_dp) THEN
            ok = .FALSE.
            RETURN
         END IF
         phis(j) = phi
         DO k = j + 1, SIZE(xis)
            bfac = zdd_div(zdd_sub(xis(k), xis(j)), zdd_sub(zdd_one, zdd_mul(zdd_conjg(xis(j)), xis(k))))
            den = zdd_sub(zdd_one, zdd_mul(zdd_conjg(phi), vals(k)))
            IF (zdd_abs2(den) < tiny_schur**2 .OR. zdd_abs2(bfac) < tiny_schur**2) THEN
               vals(k) = zdd_zero
               CYCLE
            END IF
            vals(k) = zdd_div(zdd_div(zdd_sub(vals(k), phi), den), bfac)
         END DO
      END DO

   END FUNCTION schur_parameters

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

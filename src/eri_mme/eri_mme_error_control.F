!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods aiming for error estimate and automatic cutoff calibration.
!>        integrals.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE eri_mme_error_control
   USE ao_util,                         ONLY: exp_radius
   USE eri_mme_gaussian,                ONLY: get_minimax_coeff_v_gspace,&
                                              hermite_gauss_norm
   USE eri_mme_lattice_summation,       ONLY: pgf_sum_2c_gspace_1d_deltal
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE message_passing,                 ONLY: mp_para_env_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eri_mme_error_control'

   PUBLIC :: calibrate_cutoff, cutoff_minimax_error, minimax_error, cutoff_error
CONTAINS

! **************************************************************************************************
!> \brief Find optimal cutoff minimizing errors due to minimax approximation and
!>        due to finite cutoff using bisection on the difference of the errors
!> \param hmat ...
!> \param h_inv ...
!> \param G_min ...
!> \param vol ...
!> \param zet_min   Minimum exponent
!> \param l_mm      Total ang. mom. quantum number
!> \param zet_max     Max. exponents to estimate cutoff error
!> \param l_max_zet       Max. total ang. mom. quantum numbers to estimate cutoff error
!> \param n_minimax Number of terms in minimax approximation
!> \param cutoff_l  Initial guess of lower bound for cutoff
!> \param cutoff_r  Initial guess of upper bound for cutoff
!> \param tol       Tolerance (cutoff precision)
!> \param delta     to modify initial guess interval
!> \param cutoff    Best cutoff
!> \param err_mm    Minimax error
!> \param err_c     Cutoff error
!> \param C_mm      Scaling constant to generalize AM-GM upper bound estimate to
!>                  minimax approx.
!> \param para_env ...
!> \param print_calib ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE calibrate_cutoff(hmat, h_inv, G_min, vol, zet_min, l_mm, zet_max, l_max_zet, &
                               n_minimax, cutoff_l, cutoff_r, tol, delta, &
                               cutoff, err_mm, err_c, C_mm, para_env, print_calib, unit_nr)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat, h_inv
      REAL(KIND=dp), INTENT(IN)                          :: G_min
      REAL(KIND=dp)                                      :: vol
      REAL(KIND=dp), INTENT(IN)                          :: zet_min
      INTEGER, INTENT(IN)                                :: l_mm
      REAL(KIND=dp), INTENT(IN)                          :: zet_max
      INTEGER, INTENT(IN)                                :: l_max_zet, n_minimax
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_l, cutoff_r, tol, delta
      REAL(KIND=dp), INTENT(OUT)                         :: cutoff, err_mm, err_c, C_mm
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      LOGICAL, INTENT(IN)                                :: print_calib
      INTEGER, INTENT(IN)                                :: unit_nr

      INTEGER                                            :: i, iter1, iter2, max_iter
      LOGICAL                                            :: do_print, valid_initial
      REAL(KIND=dp)                                      :: cutoff_mid, delta_c_mid, delta_mm_mid
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: minimax_aw
      REAL(KIND=dp), DIMENSION(2)                        :: cutoff_lr, delta_c, delta_mm

      do_print = unit_nr > 0 .AND. print_calib
      IF (do_print) THEN
         WRITE (unit_nr, '(/T2, A)') "ERI_MME| Basis set parameters for estimating minimax error"
         WRITE (unit_nr, '(T2, A, T67, ES12.2, 1X, I1)') "ERI_MME|   exp, l:", zet_min, l_mm
         WRITE (unit_nr, '(T2, A)') "ERI_MME| Basis set parameters for estimating cutoff error"
         WRITE (unit_nr, '(T2, A, T67, ES12.2, 1X, I1)') "ERI_MME|   exp, l:", zet_max, l_max_zet
      END IF

      max_iter = 100

      IF ((cutoff_r - cutoff_l)/(0.5_dp*(cutoff_r + cutoff_l)) .LE. tol) &
         CALL cp_abort(__LOCATION__, "difference of boundaries for cutoff "// &
                       "(MAX - MIN) must be greater than cutoff precision.")

      IF ((delta .GE. 1.0_dp) .OR. (delta .LE. 0.0_dp)) &
         CALL cp_abort(__LOCATION__, &
                       "relative delta to modify initial cutoff interval (DELTA) must be in (0, 1)")

      cutoff_lr(1) = cutoff_l
      cutoff_lr(2) = cutoff_r

      ALLOCATE (minimax_aw(2*n_minimax))

      IF (do_print) THEN
         WRITE (unit_nr, '(/T2, A)') "ERI_MME| Calibrating cutoff by bisecting error(minimax) - error(cutoff)"
         WRITE (unit_nr, '(T2, A, T72, ES9.2)') "ERI_MME| Rel. cutoff precision", tol
         WRITE (unit_nr, '(T2, A, T77, F4.1)') "ERI_MME| Rel. cutoff delta to modify initial interval", delta
      END IF

      ! 1) find valid initial values for bisection
      DO iter1 = 1, max_iter + 1
         IF (iter1 .GT. max_iter) &
            CALL cp_abort(__LOCATION__, &
                          "Maximum number of iterations in bisection to determine initial "// &
                          "cutoff interval has been exceeded.")

         cutoff_lr(1) = MAX(cutoff_lr(1), 0.5_dp*G_min**2)
         ! approx.) is hit

         DO i = 1, 2
            CALL cutoff_minimax_error(cutoff_lr(i), hmat, h_inv, vol, G_min, zet_min, l_mm, zet_max, l_max_zet, &
                                      n_minimax, minimax_aw, delta_mm(i), delta_c(i), C_mm, para_env)
         END DO

         valid_initial = .TRUE.
         IF ((delta_mm(1) - delta_c(1)) .GT. 0) THEN
            cutoff_lr(1) = cutoff_lr(1)*(1.0_dp - ABS(delta))
            valid_initial = .FALSE.
         END IF
         IF ((delta_mm(2) - delta_c(2)) .LT. 0) THEN
            cutoff_lr(2) = cutoff_lr(2)*(1.0_dp + ABS(delta))
            valid_initial = .FALSE.
         END IF

         IF (valid_initial) EXIT
      END DO

      ! 2) bisection to find cutoff s.t. err_minimax(cutoff) - err_cutoff(cutoff) = 0
      IF (do_print) WRITE (unit_nr, '(/T2, A)') &
         "ERI_MME| Step, cutoff (min, max, mid), err(minimax), err(cutoff), err diff"

      DO iter2 = 1, max_iter + 1
         IF (iter2 .GT. max_iter) &
            CALL cp_abort(__LOCATION__, &
                          "Maximum number of iterations in bisection to determine cutoff has been exceeded")

         cutoff_mid = 0.5_dp*(cutoff_lr(1) + cutoff_lr(2))
         CALL cutoff_minimax_error(cutoff_mid, hmat, h_inv, vol, G_min, zet_min, l_mm, zet_max, l_max_zet, &
                                   n_minimax, minimax_aw, delta_mm_mid, delta_c_mid, C_mm, para_env)
         IF (do_print) WRITE (unit_nr, '(T11, I2, F11.1, F11.1, F11.1, 3X, ES9.2, 3X, ES9.2, 3X, ES9.2)') &
            iter2, cutoff_lr(1), cutoff_lr(2), cutoff_mid, &
            delta_mm_mid, delta_c_mid, delta_mm_mid - delta_c_mid

         IF ((cutoff_lr(2) - cutoff_lr(1))/cutoff_mid .LT. tol) EXIT
         IF (delta_mm_mid - delta_c_mid .GT. 0) THEN
            cutoff_lr(2) = cutoff_mid
            delta_mm(2) = delta_mm_mid
            delta_c(2) = delta_c_mid
         ELSE
            cutoff_lr(1) = cutoff_mid
            delta_mm(1) = delta_mm_mid
            delta_c(1) = delta_c_mid
         END IF
      END DO
      err_mm = delta_mm_mid
      err_c = delta_c_mid
      cutoff = cutoff_mid

      IF (do_print) THEN
         WRITE (unit_nr, '(/T2, A)') "ERI_MME| Cutoff calibration number of steps:"
         WRITE (unit_nr, '(T2, A, T79, I2)') "ERI_MME|   Steps for initial interval", iter1 - 1
         WRITE (unit_nr, '(T2, A, T79, I2/)') "ERI_MME|   Bisection iteration steps", iter2 - 1
      END IF

   END SUBROUTINE calibrate_cutoff

! **************************************************************************************************
!> \brief Compute upper bounds for the errors of 2-center ERI's (P|P) due
!>        to minimax approximation and due to finite cutoff, where P is a
!>        normalized Hermite Gaussian.
!> \param cutoff ...
!> \param hmat ...
!> \param h_inv ...
!> \param vol ...
!> \param G_min ...
!> \param zet_min     Exponent of P to estimate minimax error
!> \param l_mm       total ang. mom. quantum number of P to estimate minimax error
!> \param zet_max   Max. exponents of P to estimate cutoff error
!> \param l_max_zet     Max. total ang. mom. quantum numbers of P to estimate cutoff error
!> \param n_minimax  Number of terms in minimax approximation
!> \param minimax_aw Minimax coefficients
!> \param err_mm     Minimax error
!> \param err_ctff   Cutoff error
!> \param C_mm       Scaling constant to generalize AM-GM upper bound estimate to
!>                   minimax approx.
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE cutoff_minimax_error(cutoff, hmat, h_inv, vol, G_min, zet_min, l_mm, zet_max, l_max_zet, &
                                   n_minimax, minimax_aw, err_mm, err_ctff, C_mm, para_env)
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat, h_inv
      REAL(KIND=dp), INTENT(IN)                          :: vol, G_min, zet_min
      INTEGER, INTENT(IN)                                :: l_mm
      REAL(KIND=dp), INTENT(IN)                          :: zet_max
      INTEGER, INTENT(IN)                                :: l_max_zet, n_minimax
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: minimax_aw
      REAL(KIND=dp), INTENT(OUT)                         :: err_mm, err_ctff, C_mm
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      REAL(KIND=dp)                                      :: delta_mm

      CALL minimax_error(cutoff, hmat, vol, G_min, zet_min, l_mm, &
                         n_minimax, minimax_aw, err_mm, delta_mm)
      CALL cutoff_error(cutoff, h_inv, G_min, zet_max, l_max_zet, &
                        n_minimax, minimax_aw, err_ctff, C_mm, para_env)

   END SUBROUTINE cutoff_minimax_error

! **************************************************************************************************
!> \brief   Minimax error, simple analytical formula
!>          Note minimax error may blow up for small exponents. This is also observed numerically,
!>          but in this case, error estimate is no upper bound.
!> \param cutoff ...
!> \param hmat ...
!> \param vol ...
!> \param G_min ...
!> \param zet_min    Exponent of P to estimate minimax error
!> \param l_mm       total ang. mom. quantum number of P to estimate minimax error
!> \param n_minimax  Number of terms in minimax approximation
!> \param minimax_aw Minimax coefficients
!> \param err_mm     Minimax error
!> \param delta_mm ...
!> \param potential ...
!> \param pot_par ...
! **************************************************************************************************
   SUBROUTINE minimax_error(cutoff, hmat, vol, G_min, zet_min, l_mm, &
                            n_minimax, minimax_aw, err_mm, delta_mm, potential, pot_par)
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat
      REAL(KIND=dp), INTENT(IN)                          :: vol, G_min, zet_min
      INTEGER, INTENT(IN)                                :: l_mm, n_minimax
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: minimax_aw
      REAL(KIND=dp), INTENT(OUT)                         :: err_mm, delta_mm
      INTEGER, INTENT(IN), OPTIONAL                      :: potential
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: pot_par

      INTEGER                                            :: i_xyz
      REAL(KIND=dp)                                      :: prod_mm_k

      CALL get_minimax_coeff_v_gspace(n_minimax, cutoff, G_min, minimax_aw(:), &
                                      potential=potential, pot_par=pot_par, err_minimax=delta_mm)

      prod_mm_k = 1.0_dp
      DO i_xyz = 1, 3
         prod_mm_k = prod_mm_k*(ABS(hmat(i_xyz, i_xyz))/twopi + &
                                MERGE(SQRT(2.0_dp/(zet_min*pi))*EXP(-1.0_dp), 0.0_dp, l_mm .GT. 0))
      END DO
      err_mm = 32*pi**4/vol*delta_mm*prod_mm_k

   END SUBROUTINE

! **************************************************************************************************
!> \brief Cutoff error, estimating G > G_c part of Ewald sum by using C/3 * 1/(Gx^2*Gy^2*Gz^2)^1/3 as an
!>        upper bound for 1/G^2 (AM-GM inequality) and its minimax approximation (factor C).
!>
!> Note: usually, minimax approx. falls off faster than 1/G**2, so C should be approximately 1.
!> The error is calculated for all l up to l_max and golden section search algorithm is
!> applied to find the exponent that maximizes cutoff error.
!> \param cutoff ...
!> \param h_inv ...
!> \param G_min ...
!> \param zet_max   Max. exponents of P to estimate cutoff error
!> \param l_max_zet     Max. total ang. mom. quantum numbers of P to estimate cutoff error
!> \param n_minimax  Number of terms in minimax approximation
!> \param minimax_aw Minimax coefficients
!> \param err_ctff   Cutoff error
!> \param C_mm       Scaling constant to generalize AM-GM upper bound estimate to
!>                   minimax approx.
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE cutoff_error(cutoff, h_inv, G_min, zet_max, l_max_zet, &
                           n_minimax, minimax_aw, err_ctff, C_mm, para_env)
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: h_inv
      REAL(KIND=dp), INTENT(IN)                          :: G_min, zet_max
      INTEGER, INTENT(IN)                                :: l_max_zet, n_minimax
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: minimax_aw
      REAL(KIND=dp), INTENT(OUT)                         :: err_ctff, C_mm
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      INTEGER                                            :: i_aw, iG, iter, max_iter, nG
      REAL(KIND=dp) :: C, dG, eps_zet, err0, err1, err_c, err_ctff_curr, err_ctff_prev, err_d, G, &
         G_1, G_c, gr, zet_a, zet_b, zet_c, zet_d, zet_div, zet_max_tmp

      ! parameters for finding exponent maximizing cutoff error

      eps_zet = 1.0E-05_dp ! tolerance for exponent
      zet_div = 2.0_dp ! sampling constant for finding initial values of exponents
      max_iter = 100 ! maximum number of iterations in golden section search
      G_c = SQRT(2.0*cutoff)

      zet_max_tmp = zet_max

      ! 2) Cutoff error, estimating G > G_c part of Ewald sum by using C/3 * 1/(Gx^2*Gy^2*Gz^2)^1/3 as an
      !                  upper bound for 1/G^2 (AM-GM inequality) and its minimax approximation (factor C).
      !                  Note: usually, minimax approx. falls off faster than 1/G**2, so C should be approximately 1.
      !                  The error is calculated for all l up to l_max and golden section search algorithm is
      !                  applied to find the exponent that maximizes cutoff error.
      G_1 = SQRT(1.0_dp/(3.0_dp*MINVAL(minimax_aw(1:n_minimax))))

      C_mm = 0.0_dp
      IF (G_1 .GT. G_c) THEN
         nG = 1000
         dG = (G_1 - G_c)/nG
         G = G_c
         DO iG = 1, nG
            G = MIN(G, G_c)
            C = 0.0_dp
            DO i_aw = 1, n_minimax
               C = C + 3.0_dp*minimax_aw(n_minimax + i_aw)*EXP(-3.0_dp*minimax_aw(i_aw)*G**2)*G**2
            END DO
            C_mm = MAX(C, C_mm)
            G = G + dG
         END DO
      ELSE
         DO i_aw = 1, n_minimax
            C_mm = C_mm + 3.0_dp*minimax_aw(n_minimax + i_aw)*EXP(-3.0_dp*minimax_aw(i_aw)*G_c**2)*G_c**2
         END DO
      END IF
      C = MAX(1.0_dp, C_mm)

      err_ctff_prev = 0.0_dp
      gr = 0.5_dp*(SQRT(5.0_dp) - 1.0_dp) ! golden ratio
      ! Find valid starting values for golden section search
      DO iter = 1, max_iter + 1
         IF (iter .GT. max_iter) &
            CALL cp_abort(__LOCATION__, "Maximum number of iterations for finding "// &
                          "exponent maximizing cutoff error has been exceeded.")

         CALL cutoff_error_fixed_exp(cutoff, h_inv, G_min, l_max_zet, zet_max_tmp, C, err_ctff_curr, para_env)
         IF (err_ctff_prev .GE. err_ctff_curr) THEN
            zet_a = zet_max_tmp
            zet_b = MIN(zet_max_tmp*zet_div**2, zet_max)
            EXIT
         ELSE
            err_ctff_prev = err_ctff_curr
         END IF
         zet_max_tmp = zet_max_tmp/zet_div
      END DO

      ! Golden section search
      zet_c = zet_b - gr*(zet_b - zet_a)
      zet_d = zet_a + gr*(zet_b - zet_a)
      DO iter = 1, max_iter + 1
         IF (ABS(zet_c - zet_d) .LT. eps_zet*(zet_a + zet_b)) THEN
            CALL cutoff_error_fixed_exp(cutoff, h_inv, G_min, l_max_zet, zet_a, C, err0, para_env)
            CALL cutoff_error_fixed_exp(cutoff, h_inv, G_min, l_max_zet, zet_b, C, err1, para_env)
            err_ctff_curr = MAX(err0, err1)
            EXIT
         END IF
         CALL cutoff_error_fixed_exp(cutoff, h_inv, G_min, l_max_zet, zet_c, C, err_c, para_env)
         CALL cutoff_error_fixed_exp(cutoff, h_inv, G_min, l_max_zet, zet_d, C, err_d, para_env)
         IF (err_c .GT. err_d) THEN
            zet_b = zet_d
            zet_d = zet_c
            zet_c = zet_b - gr*(zet_b - zet_a)
         ELSE
            zet_a = zet_c
            zet_c = zet_d
            zet_d = zet_a + gr*(zet_b - zet_a)
         END IF
      END DO
      err_ctff = err_ctff_curr

   END SUBROUTINE

! **************************************************************************************************
!> \brief Calculate cutoff error estimate by using C_mm/3 * 1/(Gx^2*Gy^2*Gz^2)^1/3
!>        as an upper bound for 1/G^2 (and its minimax approximation) for |G| > G_c.
!>        Error is referring to a basis function P with fixed exponent zet_max and
!>        max. angular momentum l_max_zet.
!> \param cutoff ...
!> \param h_inv ...
!> \param G_min ...
!> \param l_max_zet ...
!> \param zet_max ...
!> \param C_mm ...
!> \param err_c ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE cutoff_error_fixed_exp(cutoff, h_inv, G_min, l_max_zet, zet_max, C_mm, err_c, para_env)
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: h_inv
      REAL(KIND=dp), INTENT(IN)                          :: G_min
      INTEGER, INTENT(IN)                                :: l_max_zet
      REAL(KIND=dp), INTENT(IN)                          :: zet_max, C_mm
      REAL(KIND=dp), INTENT(OUT)                         :: err_c
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      INTEGER                                            :: ax, ay, az, G_l, G_u, Gl_first, Gl_last, &
                                                            Gu_first, Gu_last, i_xyz, l, my_p, &
                                                            n_Gl, n_Gl_left, n_Gl_p, n_Gu, &
                                                            n_Gu_left, n_Gu_p, n_p
      REAL(KIND=dp)                                      :: alpha_G, eps_G, err_c_l, G_c, G_rad, &
                                                            G_res, inv_lgth, prefactor, sum_G_diff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: S_G_l, S_G_u

      G_c = SQRT(2.0_dp*cutoff)
      eps_G = TINY(eps_G) ! sum up to machine precision
      G_res = 0.5_dp*G_min ! resolution for screening

      err_c = 0.0_dp
      alpha_G = 1.0_dp/(2.0_dp*zet_max)
      prefactor = 1.0_dp/zet_max

      ALLOCATE (S_G_l(0:2*l_max_zet, 3))
      ALLOCATE (S_G_u(0:2*l_max_zet, 3))

      G_rad = exp_radius(2*l_max_zet, alpha_G, eps_G, prefactor, epsabs=G_res)

      ! Parallelization of sum over G vectors
      my_p = para_env%mepos ! mpi rank
      n_p = para_env%num_pe ! total number of processes

      DO i_xyz = 1, 3
         inv_lgth = ABS(h_inv(i_xyz, i_xyz))

         G_l = FLOOR(G_c/(inv_lgth*twopi))
         G_u = FLOOR(G_rad/(inv_lgth*twopi))

         IF (G_u .LT. G_l) G_u = G_l

         ! Serial code:
         ! !Sum |G| <= G_c
         ! CALL pgf_sum_2c_gspace_1d_deltal(S_G_l(:, i_xyz), alpha_G, inv_lgth, -G_l, G_l, &
         !                               2.0_dp/3.0_dp, prefactor)
         ! !Sum |G| > G_c
         ! CALL pgf_sum_2c_gspace_1d_deltal(S_G_u(:, i_xyz), alpha_G, inv_lgth, G_l + 1, G_u, &
         !                               2.0_dp/3.0_dp, prefactor)

         ! Parallel code:
         n_Gu = MAX((G_u - G_l), 0)
         n_Gl = 2*G_l + 1
         n_Gu_p = n_Gu/n_p
         n_Gl_p = n_Gl/n_p
         n_Gu_left = MOD(n_Gu, n_p)
         n_Gl_left = MOD(n_Gl, n_p)

         IF (my_p .LT. n_Gu_left) THEN
            Gu_first = G_l + 1 + (n_Gu_p + 1)*my_p
            Gu_last = G_l + 1 + (n_Gu_p + 1)*(my_p + 1) - 1
         ELSE
            Gu_first = G_l + 1 + n_Gu_left + n_Gu_p*my_p
            Gu_last = G_l + 1 + n_Gu_left + n_Gu_p*(my_p + 1) - 1
         END IF

         IF (my_p .LT. n_Gl_left) THEN
            Gl_first = -G_l + (n_Gl_p + 1)*my_p
            Gl_last = -G_l + (n_Gl_p + 1)*(my_p + 1) - 1
         ELSE
            Gl_first = -G_l + n_Gl_left + n_Gl_p*my_p
            Gl_last = -G_l + n_Gl_left + n_Gl_p*(my_p + 1) - 1
         END IF

         ! Sum |G| <= G_c
         CALL pgf_sum_2c_gspace_1d_deltal(S_G_l(:, i_xyz), alpha_G, inv_lgth, Gl_first, Gl_last, &
                                          2.0_dp/3.0_dp, prefactor)

         ! Sum |G| > G_c
         CALL pgf_sum_2c_gspace_1d_deltal(S_G_u(:, i_xyz), alpha_G, inv_lgth, Gu_first, Gu_last, &
                                          2.0_dp/3.0_dp, prefactor)
      END DO

      CALL para_env%sum(S_G_l)
      CALL para_env%sum(S_G_u)

      S_G_u = S_G_u*2.0_dp ! to include negative values of G

      DO l = 0, l_max_zet
      DO ax = 0, l
      DO ay = 0, l - ax
         az = l - ax - ay

         ! Compute prod_k (S_G_l(l_k,k) + S_G_u(l_k,k)) - prod_k (S_G_l(l_k,k)) with k in {x, y, z}
         ! Note: term by term multiplication to avoid subtraction for numerical stability
         sum_G_diff = S_G_u(2*ax, 1)*S_G_u(2*ay, 2)*S_G_u(2*az, 3) + &
                      S_G_u(2*ax, 1)*S_G_u(2*ay, 2)*S_G_l(2*az, 3) + &
                      S_G_u(2*ax, 1)*S_G_l(2*ay, 2)*S_G_u(2*az, 3) + &
                      S_G_l(2*ax, 1)*S_G_u(2*ay, 2)*S_G_u(2*az, 3) + &
                      S_G_u(2*ax, 1)*S_G_l(2*ay, 2)*S_G_l(2*az, 3) + &
                      S_G_l(2*ax, 1)*S_G_u(2*ay, 2)*S_G_l(2*az, 3) + &
                      S_G_l(2*ax, 1)*S_G_l(2*ay, 2)*S_G_u(2*az, 3)

         err_c_l = 4.0_dp*pi**4*hermite_gauss_norm(zet_max, [ax, ay, az])**2* &
                   C_mm/3.0_dp*sum_G_diff

         err_c = MAX(err_c, err_c_l)
      END DO
      END DO
      END DO

      DEALLOCATE (S_G_u, S_G_l)

   END SUBROUTINE cutoff_error_fixed_exp

END MODULE eri_mme_error_control

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Ewald sums to represent integrals in direct and reciprocal lattice.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE eri_mme_lattice_summation
   #:include "eri_mme_unroll.fypp"

   USE ao_util, ONLY: exp_radius
   USE eri_mme_gaussian, ONLY: create_gaussian_overlap_dist_to_hermite, &
                               create_hermite_to_cartesian, &
                               eri_mme_coulomb, eri_mme_yukawa, &
                               eri_mme_longrange
   USE kinds, ONLY: dp, &
                    int_8
   USE mathconstants, ONLY: gaussi, &
                            pi, &
                            twopi
   USE orbital_pointers, ONLY: coset, &
                               indco, &
                               ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eri_mme_lattice_summation'

   PUBLIC :: &
      ellipsoid_bounds, &
      eri_mme_2c_get_bounds, &
      eri_mme_2c_get_rads, &
      eri_mme_3c_get_bounds, &
      eri_mme_3c_get_rads, &
      get_l, &
      pgf_sum_2c_gspace_1d, &
      pgf_sum_2c_gspace_1d_deltal, &
      pgf_sum_2c_gspace_3d, &
      pgf_sum_2c_rspace_1d, &
      pgf_sum_2c_rspace_3d, &
      pgf_sum_3c_1d, &
      pgf_sum_3c_3d

   INTEGER, PARAMETER, PRIVATE :: exp_w = 50, div_w = 10

CONTAINS

! **************************************************************************************************
!> \brief Get summation radii for 2c integrals
!> \param la_max ...
!> \param lb_max ...
!> \param zeta ...
!> \param zetb ...
!> \param a_mm ...
!> \param G_min ...
!> \param R_min ...
!> \param sum_precision ...
!> \param G_rad ...
!> \param R_rad ...
! **************************************************************************************************
   SUBROUTINE eri_mme_2c_get_rads(la_max, lb_max, zeta, zetb, a_mm, G_min, R_min, sum_precision, G_rad, R_rad)
      INTEGER, INTENT(IN)                                :: la_max, lb_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb, a_mm, G_min, R_min, &
                                                            sum_precision
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: G_rad, R_rad

      INTEGER                                            :: l_max
      REAL(KIND=dp)                                      :: alpha_G, alpha_R, G_res, R_res

      l_max = la_max + lb_max
      alpha_G = a_mm + 0.25_dp/zeta + 0.25_dp/zetb
      alpha_R = 0.25_dp/alpha_G

      G_res = 0.5_dp*G_min
      R_res = 0.5_dp*R_min

      IF (PRESENT(G_rad)) G_rad = exp_radius(l_max, alpha_G, sum_precision, 1.0_dp, epsabs=G_res)
      IF (PRESENT(R_rad)) R_rad = exp_radius(l_max, alpha_R, sum_precision, 1.0_dp, epsabs=R_res)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Get summation radii for 3c integrals
!> \param la_max ...
!> \param lb_max ...
!> \param lc_max ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param a_mm ...
!> \param G_min ...
!> \param R_min ...
!> \param sum_precision ...
!> \param G_rads_1 ...
!> \param R_rads_2 ...
!> \param R_rads_3 ...
! **************************************************************************************************
   SUBROUTINE eri_mme_3c_get_rads(la_max, lb_max, lc_max, zeta, zetb, zetc, a_mm, G_min, R_min, &
                                  sum_precision, G_rads_1, R_rads_2, R_rads_3)
      INTEGER, INTENT(IN)                                :: la_max, lb_max, lc_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb, zetc, a_mm
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: G_min, R_min
      REAL(KIND=dp), INTENT(IN)                          :: sum_precision
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT), OPTIONAL :: G_rads_1, R_rads_2
      REAL(KIND=dp), DIMENSION(2), INTENT(OUT), OPTIONAL :: R_rads_3

      REAL(KIND=dp)                                      :: alpha, alpha_R, beta, G_res, gamma, R_res

      ! exponents in G space
      alpha = 0.25_dp/zeta
      beta = 0.25_dp/zetb
      gamma = 0.25_dp/zetc + a_mm

      ! Summation radii and number of summands for all lattice summation methods
      ! sum method 1
      IF (PRESENT(G_rads_1)) THEN
         G_res = 0.5_dp*G_min
         G_rads_1(1) = exp_radius(la_max, alpha, sum_precision, 1.0_dp, G_res)
         G_rads_1(2) = exp_radius(lb_max, beta, sum_precision, 1.0_dp, G_res)
         G_rads_1(3) = exp_radius(lc_max, gamma, sum_precision, 1.0_dp, G_res)
      END IF

      ! sum method 2
      IF (PRESENT(R_rads_2)) THEN
         R_res = 0.5_dp*R_min
         R_rads_2(1) = exp_radius(lb_max + lc_max, 0.25_dp/(beta + gamma), sum_precision, 1.0_dp, epsabs=R_res)
         R_rads_2(2) = exp_radius(lc_max + la_max, 0.25_dp/(alpha + gamma), sum_precision, 1.0_dp, epsabs=R_res)
         R_rads_2(3) = exp_radius(lb_max + la_max, 0.25_dp/(alpha + beta), sum_precision, 1.0_dp, epsabs=R_res)
      END IF

      ! sum method 3

      IF (PRESENT(R_rads_3)) THEN
         R_res = 0.5_dp*R_min
         alpha_R = 1.0_dp/((zeta + zetb + zetc)/((zeta + zetb)*zetc) + 4.0_dp*a_mm)
         R_rads_3(1) = exp_radius(la_max + lb_max, zeta*zetb/(zeta + zetb), sum_precision, 1.0_dp, R_res)
         R_rads_3(2) = exp_radius(la_max + lb_max + lc_max, alpha_R, sum_precision, 1.0_dp, R_res)
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief Get summation bounds for 2c integrals
!> \param hmat ...
!> \param h_inv ...
!> \param vol ...
!> \param is_ortho ...
!> \param G_min ...
!> \param R_min ...
!> \param la_max ...
!> \param lb_max ...
!> \param zeta ...
!> \param zetb ...
!> \param a_mm ...
!> \param sum_precision ...
!> \param n_sum_1d ...
!> \param n_sum_3d ...
!> \param G_bounds ...
!> \param G_rad ...
!> \param R_bounds ...
!> \param R_rad ...
! **************************************************************************************************
   SUBROUTINE eri_mme_2c_get_bounds(hmat, h_inv, vol, is_ortho, G_min, R_min, la_max, lb_max, &
                                    zeta, zetb, a_mm, sum_precision, n_sum_1d, n_sum_3d, &
                                    G_bounds, G_rad, R_bounds, R_rad)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat, h_inv
      REAL(KIND=dp), INTENT(IN)                          :: vol
      LOGICAL, INTENT(IN)                                :: is_ortho
      REAL(KIND=dp), INTENT(IN)                          :: G_min, R_min
      INTEGER, INTENT(IN)                                :: la_max, lb_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb, a_mm, sum_precision
      INTEGER(KIND=int_8), DIMENSION(2, 3), INTENT(OUT)  :: n_sum_1d
      INTEGER(KIND=int_8), DIMENSION(2), INTENT(OUT)     :: n_sum_3d
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: G_bounds
      REAL(KIND=dp), INTENT(OUT)                         :: G_rad
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: R_bounds
      REAL(KIND=dp), INTENT(OUT)                         :: R_rad

      INTEGER                                            :: i_xyz
      REAL(KIND=dp)                                      :: ns_G, ns_R, vol_G

      CALL eri_mme_2c_get_rads(la_max, lb_max, zeta, zetb, a_mm, G_min, R_min, sum_precision, G_rad, R_rad)

      G_bounds = ellipsoid_bounds(G_rad, TRANSPOSE(hmat)/(2.0_dp*pi))
      R_bounds = ellipsoid_bounds(R_rad, h_inv)

      vol_G = twopi**3/vol

      IF (is_ortho) THEN
         DO i_xyz = 1, 3
            ns_G = 2.0_dp*G_bounds(i_xyz)
            ns_R = 2.0_dp*R_bounds(i_xyz)
            n_sum_1d(1, i_xyz) = nsum_2c_gspace_1d(ns_G, la_max, lb_max)
            n_sum_1d(2, i_xyz) = nsum_2c_rspace_1d(ns_R, la_max, lb_max)
         END DO
      ELSE
         ns_G = 4._dp/3._dp*pi*G_rad**3/vol_G
         ns_R = 4._dp/3._dp*pi*R_rad**3/vol
         n_sum_3d(1) = nsum_2c_gspace_3d(ns_G, la_max, lb_max)
         n_sum_3d(2) = nsum_2c_rspace_3d(ns_R, la_max, lb_max)
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief Get summation bounds for 3c integrals
!> \param hmat ...
!> \param h_inv ...
!> \param vol ...
!> \param is_ortho ...
!> \param G_min ...
!> \param R_min ...
!> \param la_max ...
!> \param lb_max ...
!> \param lc_max ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param a_mm ...
!> \param sum_precision ...
!> \param n_sum_1d ...
!> \param n_sum_3d ...
!> \param G_bounds_1 ...
!> \param G_rads_1 ...
!> \param R_bounds_2 ...
!> \param R_rads_2 ...
!> \param R_bounds_3 ...
!> \param R_rads_3 ...
! **************************************************************************************************
   SUBROUTINE eri_mme_3c_get_bounds(hmat, h_inv, vol, is_ortho, G_min, R_min, la_max, lb_max, lc_max, &
                                    zeta, zetb, zetc, a_mm, sum_precision, n_sum_1d, n_sum_3d, &
                                    G_bounds_1, G_rads_1, R_bounds_2, R_rads_2, R_bounds_3, R_rads_3)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat, h_inv
      REAL(KIND=dp), INTENT(IN)                          :: vol
      LOGICAL, INTENT(IN)                                :: is_ortho
      REAL(KIND=dp), INTENT(IN)                          :: G_min, R_min
      INTEGER, INTENT(IN)                                :: la_max, lb_max, lc_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb, zetc, a_mm, sum_precision
      INTEGER(KIND=int_8), DIMENSION(3, 3), INTENT(OUT)  :: n_sum_1d
      INTEGER(KIND=int_8), DIMENSION(3), INTENT(OUT)     :: n_sum_3d
      REAL(KIND=dp), DIMENSION(3, 3)                     :: G_bounds_1
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: G_rads_1
      REAL(KIND=dp), DIMENSION(3, 3)                     :: R_bounds_2
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: R_rads_2
      REAL(KIND=dp), DIMENSION(2, 3)                     :: R_bounds_3
      REAL(KIND=dp), DIMENSION(2), INTENT(OUT)           :: R_rads_3

      INTEGER                                            :: i, i_xyz, order_1, order_2
      REAL(KIND=dp)                                      :: ns1_G1, ns1_G2, ns2_G, ns2_R, ns3_R1, &
                                                            ns3_R2, vol_G
      CALL eri_mme_3c_get_rads(la_max, lb_max, lc_max, zeta, zetb, zetc, a_mm, G_min, R_min, sum_precision, &
                               G_rads_1=G_rads_1, R_rads_2=R_rads_2, R_rads_3=R_rads_3)

      vol_G = twopi**3/vol

      order_1 = MAXLOC(G_rads_1, DIM=1)
      order_2 = MINLOC(G_rads_1, DIM=1)

      DO i = 1, 3
         G_bounds_1(i, :) = ellipsoid_bounds(G_rads_1(i), TRANSPOSE(hmat)/(2.0_dp*pi))
      END DO

      DO i = 1, 3
         R_bounds_2(i, :) = ellipsoid_bounds(R_rads_2(i), h_inv)
      END DO

      DO i = 1, 2
         R_bounds_3(i, :) = ellipsoid_bounds(R_rads_3(i), h_inv)
      END DO

      IF (is_ortho) THEN
         DO i_xyz = 1, 3

            ns3_R1 = 2.0_dp*R_bounds_3(1, i_xyz)
            ns3_R2 = 2.0_dp*R_bounds_3(2, i_xyz)

            n_sum_1d(3, i_xyz) = nsum_3c_rspace_1d(ns3_R1, ns3_R2)
         END DO

      ELSE

         order_1 = MAXLOC(G_rads_1, DIM=1)
         order_2 = MINLOC(G_rads_1, DIM=1)

         SELECT CASE (order_1)
         CASE (1)
            ns1_G1 = 4._dp/3._dp*pi*G_rads_1(2)**3/vol_G
            ns1_G2 = 4._dp/3._dp*pi*G_rads_1(3)**3/vol_G
         CASE (2)
            ns1_G1 = 4._dp/3._dp*pi*G_rads_1(1)**3/vol_G
            ns1_G2 = 4._dp/3._dp*pi*G_rads_1(3)**3/vol_G
         CASE (3)
            ns1_G1 = 4._dp/3._dp*pi*G_rads_1(1)**3/vol_G
            ns1_G2 = 4._dp/3._dp*pi*G_rads_1(2)**3/vol_G
         END SELECT

         n_sum_3d(1) = nsum_3c_gspace_3d(ns1_G1, ns1_G2, la_max, lb_max, lc_max)

         ns2_G = 4._dp/3._dp*pi*G_rads_1(order_2)**3/vol_G
         ns2_R = 4._dp/3._dp*pi*R_rads_2(order_2)**3/vol

         ns3_R1 = 4._dp/3._dp*pi*R_rads_3(1)**3/vol
         ns3_R2 = 4._dp/3._dp*pi*R_rads_3(2)**3/vol

         SELECT CASE (order_2)
         CASE (1)
            n_sum_3d(2) = nsum_product_3c_gspace_3d(ns2_G, ns2_R, la_max, lb_max, lc_max)
         CASE (2)
            n_sum_3d(2) = nsum_product_3c_gspace_3d(ns2_G, ns2_R, lb_max, la_max, lc_max)
         CASE (3)
            n_sum_3d(2) = nsum_product_3c_gspace_3d(ns2_G, ns2_R, lc_max, lb_max, la_max)
         END SELECT

         n_sum_3d(3) = nsum_3c_rspace_3d(ns3_R1, ns3_R2, la_max, lb_max, lc_max)
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief Roughly estimated number of floating point operations
!> \param ns_G ...
!> \param l ...
!> \param m ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION nsum_2c_gspace_1d(ns_G, l, m)
      REAL(KIND=dp), INTENT(IN)                          :: ns_G
      INTEGER, INTENT(IN)                                :: l, m
      INTEGER(KIND=int_8)                                :: nsum_2c_gspace_1d

      nsum_2c_gspace_1d = NINT(ns_G*(2*exp_w + (l + m + 1)*5), KIND=int_8)
   END FUNCTION

! **************************************************************************************************
!> \brief Compute Ewald-like sum for 2-center ERIs in G space in 1 dimension
!>        S_G(l, alpha) = (-i)^l*inv_lgth*sum_G( C(l, alpha, G) exp(iGR) ), with
!>                        C(l, alpha, r) = r^l exp(-alpha*r^2),
!>        dG = inv_lgth*twopi and G = -G_bound*dG, (-G_bound + 1)*dG, ..., G_bound*dG
!>             for all l < = l_max.
!> \param S_G ...
!> \param R ...
!> \param alpha ...
!> \param inv_lgth ...
!> \param G_c ...
!> \note  S_G is real.
! **************************************************************************************************
   PURE SUBROUTINE pgf_sum_2c_gspace_1d(S_G, R, alpha, inv_lgth, G_c)
      REAL(KIND=dp), DIMENSION(0:), INTENT(OUT)          :: S_G
      REAL(KIND=dp), INTENT(IN)                          :: R, alpha, inv_lgth, G_c

      COMPLEX(KIND=dp)                                   :: exp_tot
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: S_G_c
      INTEGER                                            :: gg, l, l_max
      REAL(KIND=dp)                                      :: dG, G, G_pow_l

      dG = inv_lgth*twopi
      l_max = UBOUND(S_G, 1)

      ALLOCATE (S_G_c(0:l_max))
      S_G_c(:) = 0.0_dp
      DO gg = -FLOOR(G_c), FLOOR(G_c)
         G = gg*dG
         exp_tot = EXP(-alpha*G**2)*EXP(gaussi*G*R) ! cost: 2 exp_w flops
         G_pow_l = 1.0_dp
         DO l = 0, l_max
            S_G_c(l) = S_G_c(l) + G_pow_l*(-1.0_dp)**l*exp_tot ! cost: 4 flops
            G_pow_l = G_pow_l*G ! cost: 1 flop
         END DO
      END DO

      S_G(:) = REAL(S_G_c(0:l_max)*i_pow((/(l, l=0, l_max)/)))*inv_lgth
   END SUBROUTINE pgf_sum_2c_gspace_1d

! **************************************************************************************************
!> \brief Roughly estimated number of floating point operations
!> \param ns_G ...
!> \param l ...
!> \param m ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION nsum_2c_gspace_3d(ns_G, l, m)
      REAL(KIND=dp), INTENT(IN)                          :: ns_G
      INTEGER, INTENT(IN)                                :: l, m
      INTEGER(KIND=int_8)                                :: nsum_2c_gspace_3d

      nsum_2c_gspace_3d = NINT(ns_G*(2*exp_w + ncoset(l + m)*7), KIND=int_8)
   END FUNCTION

! **************************************************************************************************
!> \brief As pgf_sum_2c_gspace_1d but 3d sum required for non-orthorhombic cells
!> \param S_G ...
!> \param l_max ...
!> \param R ...
!> \param alpha ...
!> \param h_inv ...
!> \param G_c ...
!> \param G_rad ...
!> \param vol ...
!> \param coulomb ...
!> \note  MMME Method is not very efficient for non-orthorhombic cells
! **************************************************************************************************
   PURE SUBROUTINE pgf_sum_2c_gspace_3d(S_G, l_max, R, alpha, h_inv, G_c, G_rad, vol, potential, pot_par)
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: S_G
      INTEGER, INTENT(IN)                                :: l_max
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: R
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: h_inv
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: G_c
      REAL(KIND=dp), INTENT(IN)                          :: G_rad, vol
      INTEGER, INTENT(IN), OPTIONAL                      :: potential
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: pot_par

      COMPLEX(KIND=dp)                                   :: exp_tot
      INTEGER, DIMENSION(3)                              :: l_xyz
      INTEGER                                            :: gx, gy, gz, k, l, lco, lx, ly, lz
      COMPLEX(KIND=dp), DIMENSION(ncoset(l_max))         :: Ig
      REAL(KIND=dp)                                      :: G_rads_sq, G_sq
      REAL(KIND=dp), DIMENSION(3)                        :: G, G_x, G_y, G_z
      REAL(KIND=dp), DIMENSION(3, 0:l_max)               :: G_pow_l
      REAL(KIND=dp), DIMENSION(3, 3)                     :: ht

      ht = twopi*TRANSPOSE(h_inv)
      Ig(:) = 0.0_dp

      G_rads_sq = G_rad**2

      DO gx = -FLOOR(G_c(1)), FLOOR(G_c(1))
         G_x = ht(:, 1)*gx
         DO gy = -FLOOR(G_c(2)), FLOOR(G_c(2))
            G_y = ht(:, 2)*gy
            DO gz = -FLOOR(G_c(3)), FLOOR(G_c(3))
               G_z = ht(:, 3)*gz

               G = G_x + G_y + G_z
               G_sq = G(1)**2 + G(2)**2 + G(3)**2
               IF (G_sq .GT. G_rads_sq) CYCLE

               IF (PRESENT(potential)) THEN
                  IF (gx .EQ. 0 .AND. gy .EQ. 0 .AND. gz .EQ. 0) CYCLE
               END IF

               exp_tot = EXP(-alpha*G_sq)*EXP(gaussi*DOT_PRODUCT(G, R)) ! cost: 2 exp_w flops
               IF (PRESENT(potential)) THEN
                  SELECT CASE (potential)
                  CASE (eri_mme_coulomb)
                     exp_tot = exp_tot/G_sq
                  CASE (eri_mme_yukawa)
                     exp_tot = exp_tot/(G_sq + pot_par**2)
                     !exp_tot = exp_tot/G_sq
                  CASE (eri_mme_longrange)
                     exp_tot = exp_tot/G_sq*EXP(-G_sq/pot_par**2)
                  END SELECT
               END IF
               DO k = 1, 3
                  G_pow_l(k, 0) = 1.0_dp
                  DO l = 1, l_max
                     G_pow_l(k, l) = G_pow_l(k, l - 1)*G(k)
                  END DO
               END DO
               DO lco = 1, ncoset(l_max)
                  CALL get_l(lco, l, lx, ly, lz)
                  l_xyz = [lx, ly, lz]
                  Ig(coset(lx, ly, lz)) = Ig(coset(lx, ly, lz)) + & ! cost: 7 flops
                                          G_pow_l(1, lx)*G_pow_l(2, ly)*G_pow_l(3, lz)* &
                                          exp_tot*(-1.0_dp)**l*i_pow(l)
               END DO
            END DO
         END DO
      END DO
      S_G(:) = REAL(Ig(:), KIND=dp)/vol
   END SUBROUTINE pgf_sum_2c_gspace_3d

! **************************************************************************************************
!> \brief Roughly estimated number of floating point operations
!> \param ns_R ...
!> \param l ...
!> \param m ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION nsum_2c_rspace_1d(ns_R, l, m)
      REAL(KIND=dp), INTENT(IN)                          :: ns_R
      INTEGER, INTENT(IN)                                :: l, m
      INTEGER(KIND=int_8)                                :: nsum_2c_rspace_1d

      nsum_2c_rspace_1d = NINT(ns_R*(exp_w + (l + m + 1)*3), KIND=int_8)
   END FUNCTION

! **************************************************************************************************
!> \brief Compute Ewald-like sum for 2-center ERIs in R space in 1 dimension
!>        S_R(l, alpha) = SQRT(alpha/pi) sum_R'( H(l, alpha, R-R') ),
!>        with H(l, alpha, R) = (-d/dR)^l exp(-alpha*R^2),
!>        dR = lgth and R' = -R_min*dR, (-R_min + 1)*dR, ..., R_max*dR,
!>        for all l < = l_max.
!> \param S_R ...
!> \param R ...
!> \param alpha ...
!> \param lgth ...
!> \param R_c ...
!> \note  result is equivalent to pgf_sum_2c_gspace_1d with
!>              S_R(l, alpha) = S_G(l, 1/(4*alpha))
! **************************************************************************************************
   PURE SUBROUTINE pgf_sum_2c_rspace_1d(S_R, R, alpha, lgth, R_c)
      REAL(KIND=dp), DIMENSION(0:), INTENT(OUT)          :: S_R
      REAL(KIND=dp), INTENT(IN)                          :: R, alpha, lgth, R_c

      INTEGER                                            :: l, l_max, rr
      REAL(KIND=dp)                                      :: dR, exp_tot, R_pow_l, Rp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: h_to_c

      dR = lgth
      l_max = UBOUND(S_R, 1)

      ! 1) compute sum over C(l, alpha, R - R') instead of H(l, alpha, R - R')
      S_R(:) = 0.0_dp
      Rp = R - R_c*dR
      DO rr = CEILING(-R_c - R/dR), FLOOR(R_c - R/dR)
         Rp = R + rr*dR
         exp_tot = EXP(-alpha*Rp**2) ! cost: exp_w flops
         R_pow_l = 1.0_dp
         DO l = 0, l_max
            S_R(l) = S_R(l) + R_pow_l*exp_tot ! cost: 2 flops
            R_pow_l = R_pow_l*Rp ! cost: 1 flop
         END DO
      END DO

      ! 2) C --> H
      CALL create_hermite_to_cartesian(alpha, l_max, h_to_c)
      S_R = MATMUL(TRANSPOSE(h_to_c(0:l_max, 0:l_max)), S_R)*SQRT(alpha/pi)
   END SUBROUTINE pgf_sum_2c_rspace_1d

! **************************************************************************************************
!> \brief Roughly estimated number of floating point operations
!> \param ns_R ...
!> \param l ...
!> \param m ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION nsum_2c_rspace_3d(ns_R, l, m)
      REAL(KIND=dp), INTENT(IN)                          :: ns_R
      INTEGER, INTENT(IN)                                :: l, m
      INTEGER(KIND=int_8)                                :: nsum_2c_rspace_3d

      nsum_2c_rspace_3d = NINT(ns_R*(exp_w + ncoset(l + m)*(4 + ncoset(l + m)*4)), KIND=int_8)
   END FUNCTION

! **************************************************************************************************
!> \brief As pgf_sum_2c_rspace_1d but 3d sum required for non-orthorhombic cells
!> \param S_R ...
!> \param l_max ...
!> \param R ...
!> \param alpha ...
!> \param hmat ...
!> \param h_inv ...
!> \param R_c ...
!> \param R_rad ...
!> \note  MMME Method is not very efficient for non-orthorhombic cells
! **************************************************************************************************
   PURE SUBROUTINE pgf_sum_2c_rspace_3d(S_R, l_max, R, alpha, hmat, h_inv, R_c, R_rad)
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: S_R
      INTEGER, INTENT(IN)                                :: l_max
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: R
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat, h_inv
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: R_c
      REAL(KIND=dp), INTENT(IN)                          :: R_rad

      INTEGER                                            :: k, l, lco, llx, lly, llz, lx, ly, lz, &
                                                            sx, sy, sz
      REAL(KIND=dp)                                      :: exp_tot, R_rad_sq, R_sq
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: h_to_c
      REAL(KIND=dp), DIMENSION(3)                        :: R_l, R_r, Rp, Rx, Ry, Rz, s_shift
      REAL(KIND=dp), DIMENSION(3, 0:l_max)               :: R_pow_l
      REAL(KIND=dp), DIMENSION(ncoset(l_max))            :: S_R_C

      S_R(:) = 0.0_dp
      S_R_C(:) = 0.0_dp

      s_shift = MATMUL(h_inv, -R)
      R_l = -R_c + s_shift
      R_r = R_c + s_shift

      R_rad_sq = R_rad**2

      DO sx = CEILING(R_l(1)), FLOOR(R_r(1))
         Rx = hmat(:, 1)*sx
         DO sy = CEILING(R_l(2)), FLOOR(R_r(2))
            Ry = hmat(:, 2)*sy
            DO sz = CEILING(R_l(3)), FLOOR(R_r(3))
               Rz = hmat(:, 3)*sz
               Rp = Rx + Ry + Rz
               R_sq = (Rp(1) + R(1))**2 + (Rp(2) + R(2))**2 + (Rp(3) + R(3))**2
               IF (R_sq .GT. R_rad_sq) CYCLE
               exp_tot = EXP(-alpha*R_sq) ! cost: exp_w flops
               DO k = 1, 3
                  R_pow_l(k, 0) = 1.0_dp
                  DO l = 1, l_max
                     R_pow_l(k, l) = R_pow_l(k, l - 1)*(Rp(k) + R(k))
                  END DO
               END DO
               DO lco = 1, ncoset(l_max)
                  CALL get_l(lco, l, lx, ly, lz)
                  S_R_C(coset(lx, ly, lz)) = S_R_C(coset(lx, ly, lz)) + R_pow_l(1, lx)*R_pow_l(2, ly)*R_pow_l(3, lz)*exp_tot ! cost: 4 flops
               END DO
            END DO
         END DO
      END DO

      CALL create_hermite_to_cartesian(alpha, l_max, h_to_c)

      DO lco = 1, ncoset(l_max)
         CALL get_l(lco, l, lx, ly, lz)

         DO llx = 0, lx
         DO lly = 0, ly
         DO llz = 0, lz
            S_R(coset(lx, ly, lz)) = S_R(coset(lx, ly, lz)) + & ! cost: 4 flops
                                     h_to_c(llx, lx)*h_to_c(lly, ly)*h_to_c(llz, lz)* &
                                     S_R_C(coset(llx, lly, llz))
         END DO
         END DO
         END DO
      END DO
      S_R(:) = S_R(:)*(alpha/pi)**1.5_dp

   END SUBROUTINE pgf_sum_2c_rspace_3d

! **************************************************************************************************
!> \brief Compute 1d sum
!>        S_G(l, alpha) = inv_lgth*sum_G( C(l, alpha, delta_l, G) ) with
!>          C(l, alpha, delta_l, G) = prefactor*|G|^(l-delta_l) exp(-alpha*G^2)
!>          if G not equal 0
!>          C(l = 0, alpha, delta_l, 0) = 1, C(l>0, alpha, delta_l, 0) = 0
!>        dG = inv_lgth*twopi and G = -G_bound*dG, (-G_bound + 1)*dG, ..., G_bound*dG
!>        for all l < = l_max.
!> \param S_G ...
!> \param alpha ...
!> \param inv_lgth ...
!> \param G_min ...
!> \param G_c ...
!> \param delta_l ...
!> \param prefactor ...
!> \note  needed for cutoff error estimate
! **************************************************************************************************
   PURE SUBROUTINE pgf_sum_2c_gspace_1d_deltal(S_G, alpha, inv_lgth, G_min, G_c, delta_l, prefactor)
      REAL(KIND=dp), DIMENSION(0:), INTENT(OUT)          :: S_G
      REAL(KIND=dp), INTENT(IN)                          :: alpha, inv_lgth
      INTEGER, INTENT(IN)                                :: G_min, G_c
      REAL(KIND=dp), INTENT(IN)                          :: delta_l, prefactor

      INTEGER                                            :: k, k0, k1, l, l_max
      REAL(KIND=dp)                                      :: dG, exp_tot, G, prefac

      prefac = prefactor*inv_lgth
      dG = inv_lgth*twopi ! positive
      l_max = UBOUND(S_G, 1)

      S_G(:) = 0.0_dp
      IF (0 .LT. G_min) THEN
         k0 = G_min; k1 = 0
      ELSE IF (G_c .LT. 0) THEN
         k0 = 0; k1 = G_c
      ELSE ! Gmin <= 0 /\ 0 <= Gc
         S_G(0) = prefac
         k0 = 1; k1 = -1
      END IF
      ! positive range
      IF (0 .LT. k0) THEN
         DO k = k0, G_c
            G = k*dG; exp_tot = EXP(-alpha*G**2)*prefac
            DO l = 0, l_max
               S_G(l) = S_G(l) + G**(l - delta_l)*exp_tot
            END DO
         END DO
      END IF
      ! negative range
      IF (k1 .LT. 0) THEN
         DO k = G_min, k1
            G = k*dG; exp_tot = EXP(-alpha*G**2)*prefac
            DO l = 0, l_max
               S_G(l) = S_G(l) + ABS(G)**(l - delta_l)*exp_tot
            END DO
         END DO
      END IF
   END SUBROUTINE pgf_sum_2c_gspace_1d_deltal

! **************************************************************************************************
!> \brief Compute Ewald-like sum for 3-center integrals in 1 dimension
!>        S_G(l, m, n, alpha, beta, gamma) = i^(l+m+n)*(-1)^(l+m)*inv_lgth^2*
!>                                           sum_G sum_G'( exp(i G R1)
!>                                           C(l,alpha,G) C(m,beta,G'-G) C(n,gamma,G') exp(i G' R2) )
!>        for all l < = l_max, m <= m_max, n <= n_max.
!>        a_mm is the minimax exponent.
!>        alpha =  1/(4 zeta), beta = 1/(4 zetb), gamma = 1/(4 zetc) + a_mm
!>        R1 = RB-RA; R2 = RC-RB
!>        Note on method / order arguments:
!>        Three equivalent methods (Poisson summation) to compute this sum over
!>        Cartesian Gaussians C or Hermite Gaussians H and
!>        reciprocal lattice vectors G or direct lattice vectors R:
!>        - method 1: sum_G sum_G' C(G) C(G,G') C(G')
!>        - method 2: sum_G sum_R C(G) C(R)
!>        - method 3: sum_R sum_R' H(R, R')
!>        The order parameter selects the Gaussian functions over which the sum is performed
!>        method 1: order = 1, 2, 3
!>        method 2: order = 1, 2, 3
!>        method 3: order = 1
!>        If method and order are not present, the method / order that converges fastest is
!>        automatically chosen.
!> \param S_G ...
!> \param RA ...
!> \param RB ...
!> \param RC ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param a_mm ...
!> \param lgth ...
!> \param G_bounds_1 ...
!> \param R_bounds_2 ...
!> \param R_bounds_3 ...
!> \param method ...
!> \param method_out ...
!> \param order ...
! **************************************************************************************************
   SUBROUTINE pgf_sum_3c_1d(S_G, RA, RB, RC, zeta, zetb, zetc, a_mm, lgth, R_bounds_3)
      REAL(KIND=dp), DIMENSION(0:, 0:, 0:), INTENT(OUT)  :: S_G
      REAL(KIND=dp), INTENT(IN)                          :: RA, RB, RC, zeta, zetb, zetc, a_mm, lgth
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: R_bounds_3

      INTEGER                                            :: l_max, m_max, n_max
      INTEGER                                            :: nR1, nR2
      INTEGER                                            :: prop_exp

      l_max = UBOUND(S_G, 1)
      m_max = UBOUND(S_G, 2)
      n_max = UBOUND(S_G, 3)

      nR1 = 2*FLOOR(R_bounds_3(1)) + 1
      nR2 = 2*FLOOR(R_bounds_3(2)) + 1

      IF (nR1*nR2 > 1 + nR1*2) THEN
         prop_exp = 1
      ELSE
         prop_exp = 0
      END IF

      IF (MAXVAL([l_max, m_max, n_max]) .GT. ${lmax_unroll}$) THEN
         CALL pgf_sum_3c_rspace_1d_generic(S_G, RA, RB, RC, zeta, zetb, zetc, a_mm, lgth, R_bounds_3)
      ELSE
         #:for l_max in range(0,lmax_unroll+1)
            IF (l_max == ${l_max}$) THEN
               #:for m_max in range(0,lmax_unroll+1)
                  IF (m_max == ${m_max}$) THEN
                     #:for n_max in range(0, lmax_unroll+1)
                        IF (n_max == ${n_max}$) THEN
                           #:for prop_exp in range(0,2)
                              IF (prop_exp == ${prop_exp}$) THEN
                                 CALL pgf_sum_3c_rspace_1d_${l_max}$_${m_max}$_${n_max}$_exp_${prop_exp}$ (S_G, RA, RB, RC, &
                                                                                           zeta, zetb, zetc, a_mm, lgth, R_bounds_3)
                                 RETURN
                              END IF
                           #:endfor
                        END IF
                     #:endfor
                  END IF
               #:endfor
            END IF
         #:endfor
      END IF

   END SUBROUTINE pgf_sum_3c_1d

! **************************************************************************************************
!> \brief Roughly estimated number of floating point operations
!> \param ns_G1 ...
!> \param ns_G2 ...
!> \param l ...
!> \param m ...
!> \param n ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION nsum_3c_gspace_1d()
      INTEGER(KIND=int_8)                                :: nsum_3c_gspace_1d

      nsum_3c_gspace_1d = 15
   END FUNCTION

! **************************************************************************************************
!> \brief Roughly estimated number of floating point operations
!> \param ns_G ...
!> \param ns_R ...
!> \param l ...
!> \param m ...
!> \param n ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION nsum_product_3c_gspace_1d(ns_G, ns_R)
      REAL(KIND=dp), INTENT(IN) :: ns_G, ns_R
      INTEGER(KIND=int_8)                                :: nsum_product_3c_gspace_1d

      nsum_product_3c_gspace_1d = MIN(19, NINT(ns_G*(3 + ns_R*2)))
   END FUNCTION

! **************************************************************************************************
!> \brief Roughly estimated number of floating point operations
!> \param ns_R1 ...
!> \param ns_R2 ...
!> \param l ...
!> \param m ...
!> \param n ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION nsum_3c_rspace_1d(ns_R1, ns_R2)
      REAL(KIND=dp), INTENT(IN)                          :: ns_R1, ns_R2
      INTEGER(KIND=int_8)                                :: nsum_3c_rspace_1d

      nsum_3c_rspace_1d = NINT(MIN((4 + ns_R1*2), ns_R1*(ns_R2 + 1)), KIND=int_8)
   END FUNCTION

! **************************************************************************************************
!> \brief Helper routine: compute SQRT(alpha/pi) (-1)^n sum_(R, R') sum_{t=0}^{l+m} E(t,l,m) H(RC - P(R) - R', t + n, alpha)
!> with alpha = 1.0_dp/((a + b + c)/((a + b)*c) + 4.0_dp*a_mm),
!> P(R) = (a*(RA + R) + b*RB)/(a + b)
!> \param S_R ...
!> \param RA ...
!> \param RB ...
!> \param RC ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param a_mm ...
!> \param lgth ...
!> \param R_c ...
! **************************************************************************************************
   PURE SUBROUTINE pgf_sum_3c_rspace_1d_generic(S_R, RA, RB, RC, zeta, zetb, zetc, a_mm, lgth, R_c)
      REAL(KIND=dp), DIMENSION(0:, 0:, 0:), INTENT(OUT)  :: S_R
      REAL(KIND=dp), INTENT(IN)                          :: RA, RB, RC, zeta, zetb, zetc, a_mm, lgth
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: R_c

      INTEGER                                            :: ll, mm, l, k, l_max, m, m_max, n, n_max, rr1, rr2, t, l_tot_max
      REAL(KIND=dp)                                      :: alpha, dR, exp_tot, R, R1, R2, R_offset, &
                                                            R_pow_t, R_tmp, rr1_delta, rr2_delta, c1, c2, c3
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: S_R_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: h_to_c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: E

      dR = lgth
      alpha = 1.0_dp/((zeta + zetb + zetc)/((zeta + zetb)*zetc) + 4.0_dp*a_mm)
      l_max = UBOUND(S_R, 1)
      m_max = UBOUND(S_R, 2)
      n_max = UBOUND(S_R, 3)
      l_tot_max = l_max + m_max + n_max

      ALLOCATE (S_R_t(0:l_max + m_max + n_max))
      ALLOCATE (E(-1:l_max + m_max + 1, -1:l_max, -1:m_max))

      S_R(:, :, :) = 0.0_dp

      R_offset = RC - (zeta*RA + zetb*RB)/(zeta + zetb)

      ! inline CALL create_hermite_to_cartesian(alpha, l_tot_max, h_to_c)
      ALLOCATE (h_to_c(-1:l_tot_max + 1, 0:l_tot_max))
      h_to_c(:, :) = 0.0_dp
      h_to_c(0, 0) = 1.0_dp
      DO l = 0, l_tot_max - 1
         DO k = 0, l + 1
            h_to_c(k, l + 1) = -(k + 1)*h_to_c(k + 1, l) + 2.0_dp*alpha*h_to_c(k - 1, l)
         END DO
      END DO

      rr1_delta = (RA - RB)/dR
      DO rr1 = CEILING(-R_c(1) + rr1_delta), FLOOR(R_c(1) + rr1_delta)
         S_R_t(:) = 0.0_dp
         R1 = rr1*dR
         R_tmp = R_offset + R1*zeta/(zeta + zetb)
         rr2_delta = -R_tmp/dR
         DO rr2 = CEILING(-R_c(2) + rr2_delta), FLOOR(R_c(2) + rr2_delta)
            R2 = rr2*dR
            R = R_tmp + R2
            exp_tot = EXP(-alpha*R**2) ! cost: exp_w flops
            R_pow_t = 1.0_dp
            DO t = 0, l_max + m_max + n_max
               S_R_t(t) = S_R_t(t) + R_pow_t*exp_tot ! cost: 2 flops
               R_pow_t = R_pow_t*R ! cost: 1 flop
            END DO
         END DO

         ! C --> H
         S_R_t(:) = MATMUL(TRANSPOSE(h_to_c(0:l_max + m_max + n_max, 0:l_max + m_max + n_max)), S_R_t)*SQRT(alpha/pi)

         ! H --> HH
         !inline CALL create_gaussian_overlap_dist_to_hermite(l_max, m_max, zeta, zetb, RA-R1, RB, 2, E)

         E(:, :, :) = 0.0_dp
         E(0, 0, 0) = EXP(-zeta*zetb/(zeta + zetb)*(RA - R1 - RB)**2)

         c1 = 1.0_dp/(zeta + zetb)
         c2 = 2.0_dp*(zetb/(zeta + zetb))*(RB - (RA - R1))
         c3 = 2.0_dp*(zeta/(zeta + zetb))*(RA - R1 - RB)

         DO mm = 0, m_max - 1
            DO ll = 0, l_max - 1
               DO t = 0, ll + mm + 1
                  E(t, ll + 1, mm) = zeta*(c1*E(t - 1, ll, mm) + &
                                           c2*E(t, ll, mm) + &
                                           2*(t + 1)*E(t + 1, ll, mm) - &
                                           2*ll*E(t, ll - 1, mm))
                  E(t, ll, mm + 1) = zetb*(c1*E(t - 1, ll, mm) + &
                                           c3*E(t, ll, mm) + &
                                           2*(t + 1)*E(t + 1, ll, mm) - &
                                           2*mm*E(t, ll, mm - 1))
               END DO
            END DO
         END DO

         DO ll = 0, l_max - 1
            DO t = 0, ll + m_max + 1
               E(t, ll + 1, m_max) = zeta*(c1*E(t - 1, ll, m_max) + &
                                           c2*E(t, ll, m_max) + &
                                           2*(t + 1)*E(t + 1, ll, m_max) - &
                                           2*ll*E(t, ll - 1, m_max))
            END DO
         END DO

         DO mm = 0, m_max - 1
            DO t = 0, l_max + mm + 1
               E(t, l_max, mm + 1) = zetb*(c1*E(t - 1, l_max, mm) + &
                                           c3*E(t, l_max, mm) + &
                                           2*(t + 1)*E(t + 1, l_max, mm) - &
                                           2*mm*E(t, l_max, mm - 1))
            END DO
         END DO

         DO n = 0, n_max
            DO m = 0, m_max
               DO l = 0, l_max
                  DO t = 0, l + m
                     S_R(l, m, n) = S_R(l, m, n) + E(t, l, m)*(-1)**n*S_R_t(t + n) ! cost: 5 flops
                  END DO
               END DO
            END DO
         END DO
      END DO

      S_R = S_R*pi**(-0.5_dp)*((zeta + zetb)/(zeta*zetb))**(-0.5_dp)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Helper routine: compute SQRT(alpha/pi) (-1)^n sum_(R, R') sum_{t=0}^{l+m} E(t,l,m) H(RC - P(R) - R', t + n, alpha)
!> with alpha = 1.0_dp/((a + b + c)/((a + b)*c) + 4.0_dp*a_mm),
!> P(R) = (a*(RA + R) + b*RB)/(a + b)
!> \param S_R ...
!> \param RA ...
!> \param RB ...
!> \param RC ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param a_mm ...
!> \param lgth ...
!> \param R_c ...
! **************************************************************************************************
   #:for prop_exp in range(0,2)
      #:for l_max in range(0, lmax_unroll+1)
         #:for m_max in range(0, lmax_unroll+1)
            #:for n_max in range(0, lmax_unroll+1)
               #:set l_tot_max = l_max + m_max + n_max
               PURE SUBROUTINE pgf_sum_3c_rspace_1d_${l_max}$_${m_max}$_${n_max}$_exp_${prop_exp}$ ( &
                  S_R, RA, RB, RC, zeta, zetb, zetc, a_mm, lgth, R_c)
                  REAL(KIND=dp), DIMENSION(0:, 0:, 0:), INTENT(OUT)  :: S_R
                  REAL(KIND=dp), INTENT(IN)                          :: RA, RB, RC, zeta, zetb, zetc, a_mm, lgth
                  REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: R_c
                  INTEGER                                            :: rr1, rr2, rr2_l, rr2_r, rr1_l, rr1_r
                  REAL(KIND=dp)                                      :: alpha, alpha_E, dR, exp2_Rsq, R, R1, R_offset, &
                                                                        R_pow_t, R_tmp, rr1_delta, rr2_delta

                  #:if l_tot_max > 0
                     REAL(KIND=dp)                                      :: c1, c2, c3
                  #:endif
                  REAL(KIND=dp) :: ${", ".join(["S_R_t_{}".format(t) for t in range(0,l_tot_max+1)])}$
                  REAL(KIND=dp) :: ${", ".join(["S_R_t2_{}".format(t) for t in range(0,l_tot_max+1)])}$
  REAL(KIND=dp) :: ${", ".join([", ".join(["h_to_c_{}_{}".format(l1,l2) for l1 in range(0,l2+1)]) for l2 in range(0,l_tot_max+1)])}$
      REAL(KIND=dp) :: ${", ".join([", ".join([", ".join(["E_{}_{}_{}".format(t,l,m) for t in range(0,l+m+1)]) for l in range(0,l_max+1)]) for m in range(0,m_max+1)])}$

                  #:if prop_exp
                     REAL(KIND=dp) :: exp2_2RdR, exp_dRsq, exp_2dRsq !exp_E_dRsq, exp_E_2dRsq, exp_E_2RdR, exp_E_Rsq
                  #:endif

                  dR = lgth
                  alpha = 1.0_dp/((zeta + zetb + zetc)/((zeta + zetb)*zetc) + 4.0_dp*a_mm)

                  S_R(:, :, :) = 0.0_dp

                  R_offset = RC - (zeta*RA + zetb*RB)/(zeta + zetb)

                  h_to_c_0_0 = SQRT(alpha/pi)

                  #:for l in range(0, l_tot_max)
                     #:for k in range(0, l+2)
                        #:if k<l or k>0
                           h_to_c_${k}$_${l+1}$ = #{if k<l}#-${k+1}$*h_to_c_${k+1}$_${l}$#{endif}# #{if k > 0}#+2*alpha*h_to_c_${k-1}$_${l}$#{endif}#
                        #:else
                           h_to_c_${k}$_${l+1}$ = 0.0_dp
                        #:endif
                     #:endfor
                  #:endfor

                  #:if prop_exp
                     exp_dRsq = exp(-alpha*dR*dR)
                     exp_2dRsq = exp_dRsq*exp_dRsq
                  #:endif

                  rr1_delta = (RA - RB)/dR

                  rr1_l = CEILING(-R_c(1) + rr1_delta)
                  rr1_r = FLOOR(R_c(1) + rr1_delta)

                  R1 = rr1_l*dR

                  alpha_E = zeta*zetb/(zeta + zetb)

                  DO rr1 = rr1_l, rr1_r
                     #:for t in range(0, l_tot_max + 1)
                        S_R_t_${t}$ = 0.0_dp
                        S_R_t2_${t}$ = 0.0_dp
                     #:endfor
                     R_tmp = R_offset + R1*zeta/(zeta + zetb)
                     rr2_delta = -R_tmp/dR

                     rr2_l = CEILING(-R_c(2) + rr2_delta)
                     rr2_r = FLOOR(R_c(2) + rr2_delta)

                     R = R_tmp + (rr2_l)*dR

                     #:if prop_exp
                        exp2_2RdR = exp(-2*alpha*R*dR)
                        exp2_Rsq = exp(-alpha*R*R)
                     #:endif

                     DO rr2 = rr2_l, rr2_r
                        R_pow_t = 1.0_dp
                        #:if not prop_exp
                           exp2_Rsq = exp(-alpha*R*R)
                        #:endif
                        #:for t in range(0, l_tot_max + 1)
                           S_R_t_${t}$ = S_R_t_${t}$+R_pow_t*exp2_Rsq
                           #:if t<l_tot_max
                              R_pow_t = R_pow_t*R
                           #:endif
                        #:endfor

                        #:if prop_exp
                           exp2_Rsq = exp2_Rsq*exp_dRsq*exp2_2RdR
                           exp2_2RdR = exp2_2RdR*exp_2dRsq
                        #:endif
                        R = R + dR
                     END DO

                     ! C --> H
                     #:for l in range(0, l_tot_max+1)
                        #:for k in range(0, l+1)
                           S_R_t2_${l}$ = S_R_t2_${l}$+h_to_c_${k}$_${l}$*S_R_t_${k}$
                        #:endfor
                     #:endfor

                     ! H --> HH
                     E_0_0_0 = exp(-alpha_E*(RA - RB - R1)*(RA - RB - R1))

                     #:if l_tot_max > 0
                        c1 = 1.0_dp/(zeta + zetb)
                        c2 = 2.0_dp*(zetb/(zeta + zetb))*(RB - (RA - R1))
                        c3 = 2.0_dp*(zeta/(zeta + zetb))*(RA - R1 - RB)
                     #:endif

                     #:for m in range(0,m_max+1)
                        #:for l in range(0,l_max+1)
                           #:for t in range(0,l+m+2)
                              #:if l < l_max
                                 E_${t}$_${l+1}$_${m}$ = zeta*(#{if t>0}# c1*E_${t-1}$_${l}$_${m}$#{endif}# &
                                 #{if t<=l+m}# +c2*E_${t}$_${l}$_${m}$&#{endif}#
                                 #{if t<l+m}# +${2*(t+1)}$*E_${t+1}$_${l}$_${m}$ &#{endif}#
                                 #{if l>0 and t<=l-1+m}#-${2*l}$*E_${t}$_${l-1}$_${m}$#{endif}#)
                              #:endif
                              #:if m < m_max
                                 E_${t}$_${l}$_${m+1}$ = zetb*(#{if t>0}# c1*E_${t-1}$_${l}$_${m}$#{endif}# &
                                 #{if t<=l+m}#+c3*E_${t}$_${l}$_${m}$&#{endif}#
                                 #{if t<l+m}# +${2*(t+1)}$*E_${t+1}$_${l}$_${m}$ &#{endif}#
                                 #{if m>0 and t<=m-1+l}#-${2*m}$*E_${t}$_${l}$_${m-1}$#{endif}#)
                              #:endif
                           #:endfor
                        #:endfor
                     #:endfor

                     #:for n in range(0, n_max+1)
                        #:for m in range(0, m_max+1)
                           #:for l in range(0, l_max+1)
                              #:for t in range(0, l+m+1)
              S_R(${l}$, ${m}$, ${n}$) = S_R(${l}$, ${m}$, ${n}$) + E_${t}$_${l}$_${m}$*(${(-1)**n}$)*S_R_t2_${t+n}$ ! cost: 5 flops
                              #:endfor
                           #:endfor
                        #:endfor
                     #:endfor
                     R1 = R1 + dR
                  END DO

                  S_R = S_R*pi**(-0.5_dp)*((zeta + zetb)/(zeta*zetb))**(-0.5_dp)
               END SUBROUTINE
            #:endfor
         #:endfor
      #:endfor
   #:endfor

! **************************************************************************************************
!> \brief As pgf_sum_3c_1d but 3d sum required for non-orthorhombic cells
!> \param S_G ...
!> \param la_max ...
!> \param lb_max ...
!> \param lc_max ...
!> \param RA ...
!> \param RB ...
!> \param RC ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param a_mm ...
!> \param hmat ...
!> \param h_inv ...
!> \param vol ...
!> \param G_bounds_1 ...
!> \param R_bounds_2 ...
!> \param R_bounds_3 ...
!> \param G_rads_1 ...
!> \param R_rads_2 ...
!> \param R_rads_3 ...
!> \param method ...
!> \param method_out ...
!> \param order ...
! **************************************************************************************************
   SUBROUTINE pgf_sum_3c_3d(S_G, la_max, lb_max, lc_max, RA, RB, RC, &
                            zeta, zetb, zetc, a_mm, hmat, h_inv, vol, &
                            G_bounds_1, R_bounds_2, R_bounds_3, &
                            G_rads_1, R_rads_2, R_rads_3, &
                            method, method_out, order)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: S_G
      INTEGER, INTENT(IN)                                :: la_max, lb_max, lc_max
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: RA, RB, RC
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb, zetc, a_mm
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat, h_inv
      REAL(KIND=dp), INTENT(IN)                          :: vol
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: G_bounds_1, R_bounds_2
      REAL(KIND=dp), DIMENSION(2, 3), INTENT(IN)         :: R_bounds_3
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: G_rads_1, R_rads_2
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: R_rads_3
      INTEGER, INTENT(IN)                                :: method
      INTEGER, INTENT(OUT), OPTIONAL                     :: method_out
      INTEGER, INTENT(IN), OPTIONAL                      :: order

      INTEGER                                            :: l_max, m_max, n_max, sum_method, &
                                                            sum_order
      LOGICAL                                            :: assert_stm
      REAL(KIND=dp)                                      :: alpha, beta, G_rad, gamma, R_rad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: S_G_tmp
      REAL(KIND=dp), DIMENSION(3)                        :: G_bound, R1, R2, R_bound
      REAL(KIND=dp), DIMENSION(3, 3)                     :: ht

      IF (PRESENT(order)) THEN
         sum_order = order
      ELSE
         sum_order = 0
      END IF

      sum_method = method

      alpha = 0.25_dp/zeta
      beta = 0.25_dp/zetb
      gamma = 0.25_dp/zetc + a_mm

      l_max = la_max
      m_max = lb_max
      n_max = lc_max

      R1 = RB - RA
      R2 = RC - RB

      ht = twopi*TRANSPOSE(h_inv)

      SELECT CASE (sum_method)
      CASE (1) ! sum_G sum_G' C(G) C(G,G') C(G')

         IF (sum_order .EQ. 0) THEN
            sum_order = MAXLOC(G_bounds_1(:, 1), DIM=1)
            assert_stm = MINLOC(G_bounds_1(:, 1), DIM=1) == &
               MINLOC(G_bounds_1(:, 2), DIM=1) .AND. &
               MINLOC(G_bounds_1(:, 1), DIM=1) == &
               MINLOC(G_bounds_1(:, 3), DIM=1)
            CPASSERT(assert_stm)
         END IF

         CALL pgf_sum_3c_gspace_3d(S_G, l_max, m_max, n_max, R1, R2, alpha, beta, gamma, ht, vol, G_bounds_1, G_rads_1, sum_order)

      CASE (2) ! sum_G sum_R C(G) C(R)
         IF (sum_order .EQ. 0) THEN
            sum_order = MINLOC(G_bounds_1(:, 1), DIM=1)
            assert_stm = MINLOC(G_bounds_1(:, 1), DIM=1) == &
               MINLOC(G_bounds_1(:, 2), DIM=1) .AND. &
               MINLOC(G_bounds_1(:, 1), DIM=1) == &
               MINLOC(G_bounds_1(:, 3), DIM=1)
            CPASSERT(assert_stm)
         END IF
         R_rad = R_rads_2(sum_order)
         G_rad = G_rads_1(sum_order)
         R_bound = R_bounds_2(sum_order, :)
         G_bound = G_bounds_1(sum_order, :)
         SELECT CASE (sum_order)
         CASE (1)
            ALLOCATE (S_G_tmp(ncoset(l_max), ncoset(m_max), ncoset(n_max)))
            CALL pgf_sum_product_3c_gspace_3d(S_G_tmp, l_max, m_max, n_max, R1, R2, alpha, beta, gamma, hmat, h_inv, vol, &
                                              R_bound, G_bound, R_rad, G_rad)
            S_G = RESHAPE(S_G_tmp, SHAPE(S_G), order=[1, 2, 3])
         CASE (2)
            ALLOCATE (S_G_tmp(ncoset(m_max), ncoset(l_max), ncoset(n_max)))
            CALL pgf_sum_product_3c_gspace_3d(S_G_tmp, m_max, l_max, n_max, -R1, R1 + R2, beta, alpha, gamma, hmat, h_inv, vol, &
                                              R_bound, G_bound, R_rad, G_rad)
            S_G = RESHAPE(S_G_tmp, SHAPE(S_G), order=[2, 1, 3])
         CASE (3)
            ALLOCATE (S_G_tmp(ncoset(n_max), ncoset(m_max), ncoset(l_max)))
            CALL pgf_sum_product_3c_gspace_3d(S_G_tmp, n_max, m_max, l_max, -R2, -R1, gamma, beta, alpha, hmat, h_inv, vol, &
                                              R_bound, G_bound, R_rad, G_rad)
            S_G = RESHAPE(S_G_tmp, SHAPE(S_G), order=[3, 2, 1])
         END SELECT
      CASE (3) ! sum_R sum_R' H(R, R')
         CALL pgf_sum_3c_rspace_3d(S_G, l_max, m_max, n_max, RA, RB, RC, zeta, zetb, zetc, a_mm, hmat, h_inv, &
                                   R_bounds_3, R_rads_3)
         S_G = S_G*pi**(-1.5_dp)*((zeta + zetb)/(zeta*zetb))**(-1.5_dp)
      END SELECT

      IF (PRESENT(method_out)) method_out = sum_method

   END SUBROUTINE pgf_sum_3c_3d

! **************************************************************************************************
!> \brief Roughly estimated number of floating point operations
!> \param ns_G1 ...
!> \param ns_G2 ...
!> \param l ...
!> \param m ...
!> \param n ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION nsum_3c_gspace_3d(ns_G1, ns_G2, l, m, n)
      REAL(KIND=dp), INTENT(IN)                          :: ns_G1, ns_G2
      INTEGER, INTENT(IN)                                :: l, m, n
      INTEGER(KIND=int_8)                                :: nsum_3c_gspace_3d

      nsum_3c_gspace_3d = NINT(ns_G1*ns_G2*(5*exp_w + ncoset(l)*ncoset(m)*ncoset(n)*4), KIND=int_8)

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param S_G ...
!> \param l_max ...
!> \param m_max ...
!> \param n_max ...
!> \param R1 ...
!> \param R2 ...
!> \param alpha ...
!> \param beta ...
!> \param gamma ...
!> \param ht ...
!> \param vol ...
!> \param G_c ...
!> \param G_rad ...
!> \param sum_order ...
!> \param coulomb ...
! **************************************************************************************************
   PURE SUBROUTINE pgf_sum_3c_gspace_3d(S_G, l_max, m_max, n_max, R1, R2, alpha, beta, gamma, &
                                        ht, vol, G_c, G_rad, sum_order, coulomb)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: S_G
      INTEGER, INTENT(IN)                                :: l_max, m_max, n_max
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: R1, R2
      REAL(KIND=dp), INTENT(IN)                          :: alpha, beta, gamma
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: ht
      REAL(KIND=dp), INTENT(IN)                          :: vol
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: G_c
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: G_rad
      INTEGER, INTENT(IN)                                :: sum_order
      LOGICAL, INTENT(IN), OPTIONAL                      :: coulomb

      INTEGER, DIMENSION(3)                              :: G1c, G2c, G3c
      INTEGER                                            :: g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, &
                                                            g3z
      COMPLEX(KIND=dp), DIMENSION(ncoset(l_max), ncoset( &
                                  m_max), ncoset(n_max))                          :: S_G_c
      LOGICAL                                            :: use_coulomb
      REAL(KIND=dp)                                      :: G1_sq, G2_sq, G3_sq
      REAL(KIND=dp), DIMENSION(3)                        :: G1, G1_x, G1_y, G1_z, G2, G2_x, G2_y, &
                                                            G2_z, G3, G3_x, G3_y, G3_z, G_rads_sq

      S_G_c(:, :, :) = 0.0_dp

      G1c = FLOOR(G_c(1, :))
      G2c = FLOOR(G_c(2, :))
      G3c = FLOOR(G_c(3, :))

      ! we can not precompute exponentials for general cell
      ! Perform double G sum
      G_rads_sq = G_rad**2

      IF (PRESENT(coulomb)) THEN
         use_coulomb = coulomb
      ELSE
         use_coulomb = .FALSE.
      END IF

      SELECT CASE (sum_order)
      CASE (1)
         DO g2x = -G2c(1), G2c(1)
            G2_x = ht(:, 1)*g2x
            DO g2y = -G2c(2), G2c(2)
               G2_y = ht(:, 2)*g2y
               DO g2z = -G2c(3), G2c(3)
                  G2_z = ht(:, 3)*g2z
                  G2 = G2_x + G2_y + G2_z
                  G2_sq = G2(1)**2 + G2(2)**2 + G2(3)**2
                  IF (G2_sq .GT. G_rads_sq(2)) CYCLE
                  DO g3x = -G3c(1), G3c(1)
                     G3_x = ht(:, 1)*g3x
                     DO g3y = -G3c(2), G3c(2)
                        G3_y = ht(:, 2)*g3y
                        DO g3z = -G3c(3), G3c(3)
                           G3_z = ht(:, 3)*g3z
                           G3 = G3_x + G3_y + G3_z
                           G3_sq = G3(1)**2 + G3(2)**2 + G3(3)**2
                           IF (G3_sq .GT. G_rads_sq(3)) CYCLE
                           G1 = G3 - G2
                           G1_sq = G1(1)**2 + G1(2)**2 + G1(3)**2
                           IF (G1_sq .GT. G_rads_sq(1)) CYCLE
                           IF (use_coulomb) THEN
                              IF (g3x == 0 .AND. g3y == 0 .AND. g3z == 0) CYCLE
                           END IF
                           CALL pgf_product_3c_gspace_3d(S_G_c, G1, G1_sq, G2, G2_sq, G3, G3_sq, l_max, m_max, n_max, &
                                                         alpha, beta, gamma, R1, R2, use_coulomb)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      CASE (2)
         DO g1x = -G1c(1), G1c(1)
            G1_x = ht(:, 1)*g1x
            DO g1y = -G1c(2), G1c(2)
               G1_y = ht(:, 2)*g1y
               DO g1z = -G1c(3), G1c(3)
                  G1_z = ht(:, 3)*g1z
                  G1 = G1_x + G1_y + G1_z
                  G1_sq = G1(1)**2 + G1(2)**2 + G1(3)**2
                  IF (G1_sq .GT. G_rads_sq(1)) CYCLE
                  DO g3x = -G3c(1), G3c(1)
                     G3_x = ht(:, 1)*g3x
                     DO g3y = -G3c(2), G3c(2)
                        G3_y = ht(:, 2)*g3y
                        DO g3z = -G3c(3), G3c(3)
                           G3_z = ht(:, 3)*g3z
                           G3 = G3_x + G3_y + G3_z
                           G3_sq = G3(1)**2 + G3(2)**2 + G3(3)**2
                           IF (G3_sq .GT. G_rads_sq(3)) CYCLE
                           G2 = G3 - G1
                           G2_sq = G2(1)**2 + G2(2)**2 + G2(3)**2
                           IF (G2_sq .GT. G_rads_sq(2)) CYCLE
                           IF (use_coulomb) THEN
                              IF (g3x == 0 .AND. g3y == 0 .AND. g3z == 0) CYCLE
                           END IF
                           CALL pgf_product_3c_gspace_3d(S_G_c, G1, G1_sq, G2, G2_sq, G3, G3_sq, l_max, m_max, n_max, &
                                                         alpha, beta, gamma, R1, R2, use_coulomb)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      CASE (3)
         DO g1x = -G1c(1), G1c(1)
            G1_x = ht(:, 1)*g1x
            DO g1y = -G1c(2), G1c(2)
               G1_y = ht(:, 2)*g1y
               DO g1z = -G1c(3), G1c(3)
                  G1_z = ht(:, 3)*g1z
                  G1 = G1_x + G1_y + G1_z
                  G1_sq = G1(1)**2 + G1(2)**2 + G1(3)**2
                  IF (G1_sq .GT. G_rads_sq(1)) CYCLE
                  DO g2x = -G2c(1), G2c(1)
                     G2_x = ht(:, 1)*g2x
                     DO g2y = -G2c(2), G2c(2)
                        G2_y = ht(:, 2)*g2y
                        DO g2z = -G2c(3), G2c(3)
                           G2_z = ht(:, 3)*g2z
                           G2 = G2_x + G2_y + G2_z
                           G2_sq = G2(1)**2 + G2(2)**2 + G2(3)**2
                           IF (G2_sq .GT. G_rads_sq(2)) CYCLE
                           G3 = G1 + G2
                           G3_sq = G3(1)**2 + G3(2)**2 + G3(3)**2
                           IF (G3_sq .GT. G_rads_sq(3)) CYCLE
                           IF (use_coulomb) THEN
                              IF (g1x + g2x == 0 .AND. g1y + g2y == 0 .AND. g1z + g2z == 0) CYCLE
                           END IF
                           CALL pgf_product_3c_gspace_3d(S_G_c, G1, G1_sq, G2, G2_sq, G3, G3_sq, l_max, m_max, n_max, &
                                                         alpha, beta, gamma, R1, R2, use_coulomb)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END SELECT
      S_G = REAL(S_G_c, KIND=dp)/vol**2

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param S_G ...
!> \param G1 ...
!> \param G1_sq ...
!> \param G2 ...
!> \param G2_sq ...
!> \param G3 ...
!> \param G3_sq ...
!> \param l_max ...
!> \param m_max ...
!> \param n_max ...
!> \param alpha ...
!> \param beta ...
!> \param gamma ...
!> \param R1 ...
!> \param R2 ...
!> \param coulomb ...
! **************************************************************************************************
   PURE SUBROUTINE pgf_product_3c_gspace_3d(S_G, G1, G1_sq, G2, G2_sq, G3, G3_sq, l_max, m_max, n_max, &
                                            alpha, beta, gamma, R1, R2, coulomb)

      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: S_G
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: G1
      REAL(KIND=dp), INTENT(IN)                          :: G1_sq
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: G2
      REAL(KIND=dp), INTENT(IN)                          :: G2_sq
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: G3
      REAL(KIND=dp), INTENT(IN)                          :: G3_sq
      INTEGER, INTENT(IN)                                :: l_max, m_max, n_max
      REAL(KIND=dp), INTENT(IN)                          :: alpha, beta, gamma
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: R1, R2
      LOGICAL, INTENT(IN)                                :: coulomb

      COMPLEX(KIND=dp)                                   :: exp_tot
      INTEGER                                            :: k, l, lco, lx, ly, lz, m, mco, mx, my, &
                                                            mz, n, nco, nx, ny, nz
      COMPLEX(KIND=dp), DIMENSION(ncoset(n_max))         :: S_G_n
      COMPLEX(KIND=dp), DIMENSION(ncoset(m_max))         :: S_G_m
      COMPLEX(KIND=dp), DIMENSION(ncoset(l_max))         :: S_G_l
      REAL(KIND=dp), DIMENSION(3, 0:l_max)               :: G1_pow_l
      REAL(KIND=dp), DIMENSION(3, 0:m_max)               :: G2_pow_m
      REAL(KIND=dp), DIMENSION(3, 0:n_max)               :: G3_pow_n

      exp_tot = EXP(gaussi*DOT_PRODUCT(G1, R1))*EXP(-alpha*G1_sq)* & ! cost: 5 exp_w flops
                EXP(-beta*G2_sq)* &
                EXP(-gamma*G3_sq)*EXP(gaussi*DOT_PRODUCT(G3, R2))

      IF (coulomb) exp_tot = exp_tot/G3_sq

      DO k = 1, 3
         G1_pow_l(k, 0) = 1.0_dp
         DO l = 1, l_max
            G1_pow_l(k, l) = G1_pow_l(k, l - 1)*G1(k)
         END DO
         G2_pow_m(k, 0) = 1.0_dp
         DO m = 1, m_max
            G2_pow_m(k, m) = G2_pow_m(k, m - 1)*G2(k)
         END DO
         G3_pow_n(k, 0) = 1.0_dp
         DO n = 1, n_max
            G3_pow_n(k, n) = G3_pow_n(k, n - 1)*G3(k)
         END DO
      END DO

      DO lco = 1, ncoset(l_max)
         CALL get_l(lco, l, lx, ly, lz)
         S_G_l(lco) = G1_pow_l(1, lx)*G1_pow_l(2, ly)*G1_pow_l(3, lz)*(-1.0_dp)**l*i_pow(l)
      END DO

      DO mco = 1, ncoset(m_max)
         CALL get_l(mco, m, mx, my, mz)
         S_G_m(mco) = G2_pow_m(1, mx)*G2_pow_m(2, my)*G2_pow_m(3, mz)*(-1.0_dp)**m*i_pow(m)
      END DO

      DO nco = 1, ncoset(n_max)
         CALL get_l(nco, n, nx, ny, nz)
         S_G_n(nco) = G3_pow_n(1, nx)*G3_pow_n(2, ny)*G3_pow_n(3, nz)*i_pow(n)
      END DO

      DO nco = 1, ncoset(n_max)
         DO mco = 1, ncoset(m_max)
            DO lco = 1, ncoset(l_max)
               S_G(lco, mco, nco) = S_G(lco, mco, nco) + & ! cost: 4 flops
                                    S_G_l(lco)*S_G_m(mco)*S_G_n(nco)* &
                                    exp_tot
            END DO
         END DO
      END DO

   END SUBROUTINE pgf_product_3c_gspace_3d

! **************************************************************************************************
!> \brief Roughly estimated number of floating point operations
!> \param ns_G ...
!> \param ns_R ...
!> \param l ...
!> \param m ...
!> \param n ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION nsum_product_3c_gspace_3d(ns_G, ns_R, l, m, n)
      REAL(KIND=dp), INTENT(IN)                          :: ns_G, ns_R
      INTEGER, INTENT(IN)                                :: l, m, n
      INTEGER(KIND=int_8)                                :: nsum_product_3c_gspace_3d

      nsum_product_3c_gspace_3d = &
         NINT( &
         ns_G*( &
         (exp_w*2) + &
         ns_R*(exp_w*2 + ncoset(l + m)*7) + &
         3*nsum_gaussian_overlap(l, m, 1) + &
         ncoset(l)*ncoset(m)*(ncoset(l + m)*4 + ncoset(n)*8)), &
         KIND=int_8)
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param S_G ...
!> \param l_max ...
!> \param m_max ...
!> \param n_max ...
!> \param R1 ...
!> \param R2 ...
!> \param alpha ...
!> \param beta ...
!> \param gamma ...
!> \param hmat ...
!> \param h_inv ...
!> \param vol ...
!> \param R_c ...
!> \param G_c ...
!> \param R_rad ...
!> \param G_rad ...
! **************************************************************************************************
   PURE SUBROUTINE pgf_sum_product_3c_gspace_3d(S_G, l_max, m_max, n_max, R1, R2, alpha, beta, gamma, &
                                                hmat, h_inv, vol, R_c, G_c, R_rad, G_rad)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: S_G
      INTEGER, INTENT(IN)                                :: l_max, m_max, n_max
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: R1, R2
      REAL(KIND=dp), INTENT(IN)                          :: alpha, beta, gamma
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat, h_inv
      REAL(KIND=dp), INTENT(IN)                          :: vol
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: R_c, G_c
      REAL(KIND=dp), INTENT(IN)                          :: R_rad, G_rad

      COMPLEX(KIND=dp)                                   :: exp_tot
      INTEGER                                            :: gx, gy, gz, k, l, lco, lnco, lx, ly, lz, &
                                                            m, mco, n, nco
      COMPLEX(KIND=dp), &
         DIMENSION(ncoset(m_max), ncoset(n_max))         :: S_R
      COMPLEX(KIND=dp), DIMENSION(ncoset(l_max), ncoset( &
                                  m_max), ncoset(n_max))                          :: S_G_c
      REAL(KIND=dp)                                      :: G_rad_sq, G_sq, R_rad_sq
      REAL(KIND=dp), DIMENSION(3)                        :: G, G_x, G_y, G_z
      REAL(KIND=dp), DIMENSION(3, 0:l_max)               :: G_pow_l
      REAL(KIND=dp), DIMENSION(3, 3)                     :: ht
      REAL(KIND=dp), DIMENSION(ncoset(l_max))            :: S_G_c_l

      S_G_c(:, :, :) = 0.0_dp

      G_rad_sq = G_rad**2
      R_rad_sq = R_rad**2

      lnco = ncoset(l_max)

      ht = twopi*TRANSPOSE(h_inv)
      DO gx = -FLOOR(G_c(1)), FLOOR(G_c(1))
         G_x = ht(:, 1)*gx
         DO gy = -FLOOR(G_c(2)), FLOOR(G_c(2))
            G_y = ht(:, 2)*gy
            DO gz = -FLOOR(G_c(3)), FLOOR(G_c(3))
               G_z = ht(:, 3)*gz
               G = G_x + G_y + G_z
               G_sq = G(1)**2 + G(2)**2 + G(3)**2
               IF (G_sq .GT. G_rad_sq) CYCLE

               exp_tot = EXP(-alpha*G_sq)*EXP(gaussi*DOT_PRODUCT(G, R1)) ! cost: exp_w*2 flops

               DO k = 1, 3
                  G_pow_l(k, 0) = 1.0_dp
                  DO l = 1, l_max
                     G_pow_l(k, l) = G_pow_l(k, l - 1)*G(k)
                  END DO
               END DO

               CALL pgf_sum_product_3c_rspace_3d(S_R, m_max, n_max, G, R2, beta, gamma, hmat, h_inv, vol, R_c, R_rad_sq)

               DO lco = 1, ncoset(l_max)
                  CALL get_l(lco, l, lx, ly, lz)
                  S_G_c_l(lco) = G_pow_l(1, lx)*G_pow_l(2, ly)*G_pow_l(3, lz)*(-1.0_dp)**l
               END DO

               DO nco = 1, ncoset(n_max)
                  CALL get_l(nco, n)
                  DO mco = 1, ncoset(m_max)
                     CALL get_l(mco, m)
                     DO lco = 1, ncoset(l_max)
                        CALL get_l(lco, l)
                        S_G_c(lco, mco, nco) = & ! cost: 8 flops
                           S_G_c(lco, mco, nco) + &
                           S_G_c_l(lco)* &
                           exp_tot*i_pow(l + m + n)*(-1.0_dp)**m*S_R(mco, nco)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      S_G = REAL(S_G_c, KIND=dp)/vol**2

   END SUBROUTINE pgf_sum_product_3c_gspace_3d

! **************************************************************************************************
!> \brief ...
!> \param S_R ...
!> \param l_max ...
!> \param m_max ...
!> \param G ...
!> \param R ...
!> \param alpha ...
!> \param beta ...
!> \param hmat ...
!> \param h_inv ...
!> \param vol ...
!> \param R_c ...
!> \param R_rad_sq ...
! **************************************************************************************************
   PURE SUBROUTINE pgf_sum_product_3c_rspace_3d(S_R, l_max, m_max, G, R, alpha, beta, hmat, h_inv, vol, R_c, R_rad_sq)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(OUT)     :: S_R
      INTEGER, INTENT(IN)                                :: l_max, m_max
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: G, R
      REAL(KIND=dp), INTENT(IN)                          :: alpha, beta
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat, h_inv
      REAL(KIND=dp), INTENT(IN)                          :: vol
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: R_c
      REAL(KIND=dp), INTENT(IN)                          :: R_rad_sq

      COMPLEX(KIND=dp)                                   :: exp_tot
      INTEGER                                            :: k, l, lco, lx, ly, lz, m, mco, mx, my, &
                                                            mz, sx, sy, sz, t, tco, tx, ty, tz
      COMPLEX(KIND=dp), DIMENSION(ncoset(l_max + m_max))   :: S_R_t
      REAL(KIND=dp)                                      :: c1, c2, Rp_sq
      REAL(KIND=dp), &
         DIMENSION(-1:l_max + m_max + 1, -1:l_max, -1:m_max) :: E1, E2, E3
      REAL(KIND=dp), DIMENSION(3)                        :: R_l, R_r, Rp, Rx, Ry, Rz, s_shift
      REAL(KIND=dp), DIMENSION(3, 0:l_max + m_max)         :: R_pow_t

      c1 = 0.25_dp/(alpha + beta)
      c2 = alpha/(alpha + beta)

      S_R_t(:) = 0.0_dp
      S_R(:, :) = 0.0_dp

      s_shift = MATMUL(h_inv, R)
      R_l = -R_c + s_shift
      R_r = R_c + s_shift

      DO sx = CEILING(R_l(1)), FLOOR(R_r(1))
         Rx = hmat(:, 1)*sx
         DO sy = CEILING(R_l(2)), FLOOR(R_r(2))
            Ry = hmat(:, 2)*sy
            DO sz = CEILING(R_l(3)), FLOOR(R_r(3))
               Rz = hmat(:, 3)*sz
               Rp = Rx + Ry + Rz - R
               Rp_sq = Rp(1)**2 + Rp(2)**2 + Rp(3)**2
               IF (Rp_sq .GT. R_rad_sq) CYCLE

               exp_tot = EXP(-c1*Rp_sq)*EXP(-gaussi*c2*DOT_PRODUCT(Rp, G)) ! cost: exp_w*2 flops
               DO k = 1, 3
                  R_pow_t(k, 0) = 1.0_dp
                  DO t = 1, l_max + m_max
                     R_pow_t(k, t) = R_pow_t(k, t - 1)*Rp(k)
                  END DO
               END DO

               DO tco = 1, ncoset(l_max + m_max)
                  CALL get_l(tco, t, tx, ty, tz)
                  S_R_t(tco) = S_R_t(tco) + & ! cost: 7 flops
                               R_pow_t(1, tx)*R_pow_t(2, ty)*R_pow_t(3, tz)* &
                               (-1.0_dp)**t*i_pow(t)*exp_tot
               END DO

            END DO
         END DO
      END DO

      CALL create_gaussian_overlap_dist_to_hermite(l_max, m_max, alpha, beta, G(1), 0.0_dp, 1, E1)
      CALL create_gaussian_overlap_dist_to_hermite(l_max, m_max, alpha, beta, G(2), 0.0_dp, 1, E2)
      CALL create_gaussian_overlap_dist_to_hermite(l_max, m_max, alpha, beta, G(3), 0.0_dp, 1, E3)

      DO mco = 1, ncoset(m_max)
         CALL get_l(mco, m, mx, my, mz)
         DO lco = 1, ncoset(l_max)
            CALL get_l(lco, l, lx, ly, lz)
            DO tx = 0, lx + mx
            DO ty = 0, ly + my
            DO tz = 0, lz + mz
               tco = coset(tx, ty, tz)
               S_R(lco, mco) = S_R(lco, mco) + & ! cost: 4 flops
                               E1(tx, lx, mx)*E2(ty, ly, my)*E3(tz, lz, mz)*S_R_t(tco)

            END DO
            END DO
            END DO
         END DO
      END DO

      S_R(:, :) = S_R(:, :)*vol/(twopi)**3*(pi/(alpha + beta))**1.5_dp

   END SUBROUTINE pgf_sum_product_3c_rspace_3d

! **************************************************************************************************
!> \brief Roughly estimated number of floating point operations
!> \param ns_R1 ...
!> \param ns_R2 ...
!> \param l ...
!> \param m ...
!> \param n ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION nsum_3c_rspace_3d(ns_R1, ns_R2, l, m, n)
      REAL(KIND=dp), INTENT(IN)                          :: ns_R1, ns_R2
      INTEGER, INTENT(IN)                                :: l, m, n
      INTEGER(KIND=int_8)                                :: nsum_3c_rspace_3d

      nsum_3c_rspace_3d = &
         NINT( &
         ns_R1*( &
         ns_R2*(exp_w + ncoset(l + m + n)*4) + &
         3*nsum_gaussian_overlap(l, m, 2) + &
         ncoset(m)*ncoset(l)*( &
         ncoset(l + m)*2 + ncoset(n)*ncoset(l + m)*4)), &
         KIND=int_8)

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param S_R ...
!> \param l_max ...
!> \param m_max ...
!> \param n_max ...
!> \param RA ...
!> \param RB ...
!> \param RC ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param a_mm ...
!> \param hmat ...
!> \param h_inv ...
!> \param R_c ...
!> \param R_rad ...
! **************************************************************************************************
   SUBROUTINE pgf_sum_3c_rspace_3d(S_R, l_max, m_max, n_max, RA, RB, RC, zeta, zetb, zetc, a_mm, &
                                   hmat, h_inv, R_c, R_rad)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: S_R
      INTEGER, INTENT(IN)                                :: l_max, m_max, n_max
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: RA, RB, RC
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb, zetc, a_mm
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat, h_inv
      REAL(KIND=dp), DIMENSION(2, 3), INTENT(IN)         :: R_c
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: R_rad

      INTEGER                                            :: k, l, lco, lx, ly, lz, m, mco, mx, my, &
                                                            mz, n, nco, nx, ny, nz, s1x, s1y, s1z, &
                                                            s2x, s2y, s2z, t, tco, tnco, ttco, &
                                                            ttx, tty, ttz, tx, ty, tz
      REAL(KIND=dp)                                      :: alpha, exp_tot, R1_sq, R_sq
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: h_to_c
      REAL(KIND=dp), &
         DIMENSION(-1:l_max + m_max + 1, -1:l_max, -1:m_max) :: E1, E2, E3
      REAL(KIND=dp), DIMENSION(2)                        :: R_rad_sq
      REAL(KIND=dp), DIMENSION(3)                        :: R, R1, R1_l, R1_r, R1_tmp, R1x, R1y, &
                                                            R1z, R2_l, R2_r, R2x, R2y, R2z, &
                                                            R_offset, R_tmp, s1_shift, s2_shift
      REAL(KIND=dp), DIMENSION(3, 0:l_max + m_max + n_max)   :: R_pow_t
      REAL(KIND=dp), DIMENSION(ncoset(l_max + m_max), &
                               ncoset(l_max), ncoset(m_max))                   :: Eco
      REAL(KIND=dp), &
         DIMENSION(ncoset(l_max + m_max + n_max))            :: S_R_t
      REAL(KIND=dp), DIMENSION(ncoset(l_max + m_max + n_max) &
                               , ncoset(l_max + m_max + n_max))                    :: h_to_c_tco

      alpha = 1.0_dp/((zeta + zetb + zetc)/((zeta + zetb)*zetc) + 4.0_dp*a_mm)

      Eco(:, :, :) = 0.0_dp
      S_R(:, :, :) = 0.0_dp
      h_to_c_tco(:, :) = 0.0_dp

      R_offset = RC - (zeta*RA + zetb*RB)/(zeta + zetb)

      CALL create_hermite_to_cartesian(alpha, l_max + m_max + n_max, h_to_c)

      DO tco = 1, ncoset(l_max + m_max + n_max)
         CALL get_l(tco, t, tx, ty, tz)
         DO ttx = 0, tx
         DO tty = 0, ty
         DO ttz = 0, tz
            ttco = coset(ttx, tty, ttz)
            h_to_c_tco(ttco, tco) = h_to_c(ttx, tx)*h_to_c(tty, ty)*h_to_c(ttz, tz)
         END DO
         END DO
         END DO
      END DO

      s1_shift = MATMUL(h_inv, RA - RB)
      R1_l = -R_c(1, :) + s1_shift
      R1_r = R_c(1, :) + s1_shift

      R_rad_sq = R_rad**2

      DO s1x = CEILING(R1_l(1)), FLOOR(R1_r(1))
         R1x = hmat(:, 1)*s1x
         DO s1y = CEILING(R1_l(2)), FLOOR(R1_r(2))
            R1y = hmat(:, 2)*s1y
            DO s1z = CEILING(R1_l(3)), FLOOR(R1_r(3))
               R1z = hmat(:, 3)*s1z
               R1 = R1x + R1y + R1z
               S_R_t(:) = 0.0_dp
               R1_tmp = R1 - (RA - RB)
               R1_sq = R1_tmp(1)**2 + R1_tmp(2)**2 + R1_tmp(3)**2

               IF (R1_sq .GT. R_rad_sq(1)) CYCLE

               R_tmp = R_offset + R1*zeta/(zeta + zetb)
               s2_shift = -MATMUL(h_inv, R_tmp)
               R2_l = -R_c(2, :) + s2_shift
               R2_r = R_c(2, :) + s2_shift
               DO s2x = CEILING(R2_l(1)), FLOOR(R2_r(1))
                  R2x = hmat(:, 1)*s2x
                  DO s2y = CEILING(R2_l(2)), FLOOR(R2_r(2))
                     R2y = hmat(:, 2)*s2y
                     DO s2z = CEILING(R2_l(3)), FLOOR(R2_r(3))
                        R2z = hmat(:, 3)*s2z
                        R = R_tmp + R2x + R2y + R2z
                        R_sq = R(1)**2 + R(2)**2 + R(3)**2

                        IF (R_sq .GT. R_rad_sq(2)) CYCLE

                        exp_tot = EXP(-alpha*R_sq) ! cost: exp_w flops
                        DO k = 1, 3
                           R_pow_t(k, 0) = 1.0_dp
                           DO t = 1, l_max + m_max + n_max
                              R_pow_t(k, t) = R_pow_t(k, t - 1)*R(k)
                           END DO
                        END DO
                        DO tco = 1, ncoset(l_max + m_max + n_max)
                           CALL get_l(tco, t, tx, ty, tz)
                           S_R_t(tco) = S_R_t(tco) + R_pow_t(1, tx)*R_pow_t(2, ty)*R_pow_t(3, tz)*exp_tot ! cost: 4 flops
                        END DO
                     END DO
                  END DO
               END DO

               S_R_t(:) = MATMUL(TRANSPOSE(h_to_c_tco), S_R_t)*(alpha/pi)**1.5_dp

               CALL create_gaussian_overlap_dist_to_hermite(l_max, m_max, zeta, zetb, RA(1) - R1(1), RB(1), 2, E1)
               CALL create_gaussian_overlap_dist_to_hermite(l_max, m_max, zeta, zetb, RA(2) - R1(2), RB(2), 2, E2)
               CALL create_gaussian_overlap_dist_to_hermite(l_max, m_max, zeta, zetb, RA(3) - R1(3), RB(3), 2, E3)

               DO mco = 1, ncoset(m_max)
                  CALL get_l(mco, m, mx, my, mz)
                  DO lco = 1, ncoset(l_max)
                     CALL get_l(lco, l, lx, ly, lz)
                     DO tx = 0, lx + mx
                     DO ty = 0, ly + my
                     DO tz = 0, lz + mz
                        tco = coset(tx, ty, tz)
                        Eco(tco, lco, mco) = E1(tx, lx, mx)*E2(ty, ly, my)*E3(tz, lz, mz) ! cost: 2 flops
                     END DO
                     END DO
                     END DO
                  END DO
               END DO

               DO nco = 1, ncoset(n_max)
                  CALL get_l(nco, n, nx, ny, nz)
                  DO tco = 1, ncoset(l_max + m_max)
                     CALL get_l(tco, t, tx, ty, tz)
                     tnco = coset(tx + nx, ty + ny, tz + nz)
                     S_R(:, :, nco) = S_R(:, :, nco) + & ! cost: 4 flops
                                      Eco(tco, :, :)* &
                                      (-1)**n*S_R_t(tnco)

                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE pgf_sum_3c_rspace_3d

! **************************************************************************************************
!> \brief Compute bounding box for ellipsoid. This is needed in order to find summation bounds for
!>        sphere for sums over non-orthogonal lattice vectors.
!> \param s_rad sphere radius
!> \param s_to_e sphere to ellipsoid trafo
!> \return ...
! **************************************************************************************************
   PURE FUNCTION ellipsoid_bounds(s_rad, s_to_e)
      REAL(KIND=dp), INTENT(IN)                          :: s_rad
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: s_to_e
      REAL(KIND=dp), DIMENSION(3)                        :: ellipsoid_bounds

      INTEGER                                            :: i_xyz

      DO i_xyz = 1, 3
         ellipsoid_bounds(i_xyz) = SQRT(s_to_e(i_xyz, 1)**2 + s_to_e(i_xyz, 2)**2 + s_to_e(i_xyz, 3)**2)*s_rad
      END DO

   END FUNCTION ellipsoid_bounds

! **************************************************************************************************
!> \brief Roughly estimated number of floating point operations
!> \param l ...
!> \param m ...
!> \param H_or_C_product ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION nsum_gaussian_overlap(l, m, H_or_C_product)
      INTEGER, INTENT(IN)                                :: l, m, H_or_C_product
      INTEGER                                            :: nsum_gaussian_overlap

      INTEGER                                            :: loop

      nsum_gaussian_overlap = exp_w
      loop = (m + 1)*(l + 1)*(m + l + 2)
      IF (H_or_C_product .EQ. 1) THEN
         nsum_gaussian_overlap = nsum_gaussian_overlap + loop*16
      ELSE
         nsum_gaussian_overlap = nsum_gaussian_overlap + loop*32
      END IF
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param lco ...
!> \param l ...
!> \param lx ...
!> \param ly ...
!> \param lz ...
! **************************************************************************************************
   PURE ELEMENTAL SUBROUTINE get_l(lco, l, lx, ly, lz)
      INTEGER, INTENT(IN)                                :: lco
      INTEGER, INTENT(OUT)                               :: l
      INTEGER, INTENT(OUT), OPTIONAL                     :: lx, ly, lz

      l = SUM(indco(:, lco))
      IF (PRESENT(lx)) lx = indco(1, lco)
      IF (PRESENT(ly)) ly = indco(2, lco)
      IF (PRESENT(lz)) lz = indco(3, lco)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param i ...
!> \return ...
! **************************************************************************************************
   PURE ELEMENTAL FUNCTION i_pow(i)
      INTEGER, INTENT(IN)                                :: i
      COMPLEX(KIND=dp)                                   :: i_pow

      COMPLEX(KIND=dp), DIMENSION(0:3), PARAMETER :: &
         ip = (/(1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), (-1.0_dp, 0.0_dp), (0.0_dp, -1.0_dp)/)

      i_pow = ip(MOD(i, 4))

   END FUNCTION

END MODULE eri_mme_lattice_summation

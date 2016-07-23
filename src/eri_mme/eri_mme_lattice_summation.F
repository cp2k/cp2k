!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Ewald sums to represent integrals in direct and reciprocal lattice.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE eri_mme_lattice_summation

   USE ao_util,                         ONLY: exp_radius
   USE eri_mme_gaussian,                ONLY: create_gaussian_overlap_dist_to_hermite,&
                                              create_hermite_to_cartesian
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: gaussi,&
                                              pi,&
                                              twopi
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eri_mme_lattice_summation'

   ! array storing some powers of imaginary unit (in lack of a more elegant way)
   COMPLEX(KIND=dp), DIMENSION(0:23), PARAMETER :: &
      i_pow = [(1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), (-1.0_dp, 0.0_dp), (0.0_dp, -1.0_dp), &
               (1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), (-1.0_dp, 0.0_dp), (0.0_dp, -1.0_dp), &
               (1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), (-1.0_dp, 0.0_dp), (0.0_dp, -1.0_dp), &
               (1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), (-1.0_dp, 0.0_dp), (0.0_dp, -1.0_dp), &
               (1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), (-1.0_dp, 0.0_dp), (0.0_dp, -1.0_dp), &
               (1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), (-1.0_dp, 0.0_dp), (0.0_dp, -1.0_dp)]

   PUBLIC :: pgf_sum_2c_gspace_1d, &
             pgf_sum_2c_rspace_1d, &
             pgf_sum_2c_gspace_3d, &
             pgf_sum_2c_rspace_3d, &
             pgf_sum_2c_gspace_1d_deltal, &
             pgf_sum_3c_1d

   REAL(KIND=dp), PARAMETER, PUBLIC :: lattice_sum_prec = 1.0E-20_dp ! truncation threshold for summing PGFs

CONTAINS

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
!> \param G_min ...
!> \param G_c ...
!> \note  S_G is real.
! **************************************************************************************************
   PURE SUBROUTINE pgf_sum_2c_gspace_1d(S_G, R, alpha, inv_lgth, G_min, G_c)
      REAL(KIND=dp), DIMENSION(0:), INTENT(INOUT)        :: S_G
      REAL(KIND=dp), INTENT(IN)                          :: R, alpha, inv_lgth
      INTEGER, INTENT(IN)                                :: G_min, G_c

      CHARACTER(LEN=*), PARAMETER :: routineN = 'pgf_sum_2c_gspace_1d', &
         routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: exp_tot
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: S_G_c
      INTEGER                                            :: k, l, l_max
      REAL(KIND=dp)                                      :: dG, G, G_pow_l

      dG = inv_lgth*twopi
      l_max = UBOUND(S_G, 1)

      ALLOCATE (S_G_c(0:l_max))
      S_G_c(:) = 0.0_dp
      DO k = G_min, G_c
         G = k*dG
         exp_tot = EXP(-alpha*G**2)*EXP(gaussi*G*R)
         G_pow_l = 1.0_dp
         DO l = 0, l_max
            S_G_c(l) = S_G_c(l)+G_pow_l*(-1.0_dp)**l*exp_tot
            G_pow_l = G_pow_l*G
         ENDDO
      ENDDO

      S_G(:) = REAL(S_G_c(0:l_max)*i_pow(0:l_max))*inv_lgth
   END SUBROUTINE pgf_sum_2c_gspace_1d

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
!> \param R_min ...
!> \param R_max ...
!> \note  result is equivalent to pgf_sum_2c_gspace_1d with
!>              S_R(l, alpha) = S_G(l, 1/(4*alpha))
! **************************************************************************************************
   PURE SUBROUTINE pgf_sum_2c_rspace_1d(S_R, R, alpha, lgth, R_min, R_max)
      REAL(KIND=dp), DIMENSION(0:), INTENT(INOUT)        :: S_R
      REAL(KIND=dp), INTENT(IN)                          :: R, alpha, lgth
      INTEGER, INTENT(IN)                                :: R_min, R_max

      CHARACTER(LEN=*), PARAMETER :: routineN = 'pgf_sum_2c_rspace_1d', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: k, l, l_max
      REAL(KIND=dp)                                      :: dR, exp_tot, R_pow_l, Rp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: h_to_c

      dR = lgth
      l_max = UBOUND(S_R, 1)

      ! 1) compute sum over C(l, alpha, R - R') instead of H(l, alpha, R - R')
      S_R(:) = 0.0_dp
      Rp = R+R_min*dR
      DO k = R_min, R_max
         Rp = R+k*dR
         exp_tot = EXP(-alpha*Rp**2)
         R_pow_l = 1.0_dp
         DO l = 0, l_max
            S_R(l) = S_R(l)+R_pow_l*exp_tot
            R_pow_l = R_pow_l*Rp
         ENDDO
      ENDDO

      ! 2) C --> H
      CALL create_hermite_to_cartesian(alpha, l_max, h_to_c)
      S_R = MATMUL(TRANSPOSE(h_to_c), S_R)*SQRT(alpha/pi)
   END SUBROUTINE pgf_sum_2c_rspace_1d

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
      REAL(KIND=dp), DIMENSION(0:), INTENT(INOUT)        :: S_G
      REAL(KIND=dp), INTENT(IN)                          :: alpha, inv_lgth
      INTEGER, INTENT(IN)                                :: G_min, G_c
      REAL(KIND=dp), INTENT(IN)                          :: delta_l, prefactor

      CHARACTER(LEN=*), PARAMETER :: routineN = 'pgf_sum_2c_gspace_1d_deltal', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: k, l, l_max
      REAL(KIND=dp)                                      :: dG, exp_tot, G, prefac

      prefac = prefactor*inv_lgth
      dG = inv_lgth*twopi
      l_max = UBOUND(S_G, 1)

      S_G(:) = 0.0_dp
      DO k = G_min, G_c
         IF (k .NE. 0) THEN
            G = k*dG
            exp_tot = EXP(-alpha*G**2)*prefac
            DO l = 0, l_max
               S_G(l) = S_G(l)+ABS(G)**(l-delta_l)*exp_tot
            ENDDO
         ELSE
            S_G(0) = S_G(0)+prefac
         ENDIF
      ENDDO

   END SUBROUTINE pgf_sum_2c_gspace_1d_deltal

! **************************************************************************************************
!> \brief As pgf_sum_2c_gspace_1d but 3d sum required for non-orthorhombic cells
!> \param S_G ...
!> \param R ...
!> \param alpha ...
!> \param h_inv ...
!> \param G_min ...
!> \param G_c ...
!> \param vol ...
!> \note  MMME Method is not very efficient for non-orthorhombic cells
! **************************************************************************************************
   PURE SUBROUTINE pgf_sum_2c_gspace_3d(S_G, R, alpha, h_inv, G_min, G_c, vol)
      REAL(KIND=dp), DIMENSION(0:, 0:, 0:), &
         INTENT(INOUT)                                   :: S_G
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: R
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: h_inv
      INTEGER, DIMENSION(3), INTENT(IN)                  :: G_min, G_c
      REAL(KIND=dp), INTENT(IN)                          :: vol

      COMPLEX(KIND=dp)                                   :: exp_tot
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: Ig
      INTEGER                                            :: gx, gy, gz, k, l, l_max, lx, ly, lz
      INTEGER, DIMENSION(3)                              :: g_xyz, l_xyz
      REAL(KIND=dp)                                      :: G_sq, Ik
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: G_pow_l
      REAL(KIND=dp), DIMENSION(3)                        :: G
      REAL(KIND=dp), DIMENSION(3, 3)                     :: ht

      l_max = UBOUND(S_G, 1)
      ht = twopi*TRANSPOSE(h_inv)
      ALLOCATE (G_pow_l(3, 0:l_max))
      ALLOCATE (Ig(0:l_max, 0:l_max, 0:l_max))
      Ig(:, :, :) = 0.0_dp
      DO gx = G_min(1), G_c(1)
      DO gy = G_min(2), G_c(2)
      DO gz = G_min(3), G_c(3)
         g_xyz = [gx, gy, gz]
         G = MATMUL(ht, g_xyz)
         G_sq = G(1)**2+G(2)**2+G(3)**2
         exp_tot = EXP(-alpha*G_sq)*EXP(-gaussi*DOT_PRODUCT(G, -R))
         DO k = 1, 3
            G_pow_l(k, 0) = 1.0_dp
            DO l = 1, l_max
               G_pow_l(k, l) = G_pow_l(k, l-1)*G(k)
            ENDDO
         ENDDO
         IF (.NOT. (gx .EQ. 0 .AND. gy .EQ. 0 .AND. gz .EQ. 0)) THEN
            DO l = 0, l_max
            DO lx = 0, l
            DO ly = 0, l-lx
               lz = l-lx-ly
               l_xyz = [lx, ly, lz]
               Ik = 1.0_dp
               DO k = 1, 3
                  Ik = Ik*G_pow_l(k, l_xyz(k))
               ENDDO
               Ig(lx, ly, lz) = Ig(lx, ly, lz)+exp_tot/G_sq*Ik*(-gaussi)**l/vol
            ENDDO
            ENDDO
            ENDDO
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      S_G(:, :, :) = REAL(Ig(:, :, :), KIND=dp)
   END SUBROUTINE pgf_sum_2c_gspace_3d

! **************************************************************************************************
!> \brief As pgf_sum_2c_rspace_1d but 3d sum required for non-orthorhombic cells
!> \param S_R ...
!> \param R ...
!> \param alpha ...
!> \param hmat ...
!> \param R_min ...
!> \param R_max ...
!> \note  MMME Method is not very efficient for non-orthorhombic cells
! **************************************************************************************************
   PURE SUBROUTINE pgf_sum_2c_rspace_3d(S_R, R, alpha, hmat, R_min, R_max)
      REAL(KIND=dp), DIMENSION(0:, 0:, 0:), &
         INTENT(INOUT)                                   :: S_R
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: R
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: hmat
      INTEGER, DIMENSION(3), INTENT(IN)                  :: R_min, R_max

      INTEGER                                            :: k, l, l_max, ll, llx, lly, llz, lx, ly, &
                                                            lz, sx, sy, sz
      INTEGER, DIMENSION(3)                              :: l_xyz, s_xyz
      REAL(KIND=dp)                                      :: exp_tot, Ik, R_sq
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: h_to_c, R_pow_l
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Is, S_R_C
      REAL(KIND=dp), DIMENSION(3)                        :: Rp

      l_max = UBOUND(S_R, 1)
      ALLOCATE (S_R_C(0:l_max, 0:l_max, 0:l_max))
      ALLOCATE (R_pow_l(3, 0:l_max))
      ALLOCATE (Is(0:l_max, 0:l_max, 0:l_max))

      S_R(:, :, :) = 0.0_dp
      S_R_C(:, :, :) = 0.0_dp
      Is(:, :, :) = 0.0_dp

      DO sx = R_min(1), R_max(1)
      DO sy = R_min(2), R_max(2)
      DO sz = R_min(3), R_max(3)
         s_xyz = [sx, sy, sz]
         Rp = MATMUL(hmat, s_xyz)
         R_sq = (Rp(1)+R(1))**2+(Rp(2)+R(2))**2+(Rp(3)+R(3))**2
         exp_tot = EXP(-alpha*R_sq)
         DO k = 1, 3
            R_pow_l(k, 0) = 1.0_dp
            DO l = 1, l_max
               R_pow_l(k, l) = R_pow_l(k, l-1)*(Rp(k)+R(k))
            ENDDO
         ENDDO
         DO l = 0, l_max
         DO lx = 0, l
         DO ly = 0, l-lx
            lz = l-lx-ly
            l_xyz = [lx, ly, lz]
            Ik = 1.0_dp
            DO k = 1, 3
               Ik = Ik*R_pow_l(k, l_xyz(k))
            ENDDO
            Is(lx, ly, lz) = Is(lx, ly, lz)+exp_tot*Ik
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      S_R_C(:, :, :) = Is(:, :, :)

      CALL create_hermite_to_cartesian(alpha, l_max, h_to_c)
      DO l = 0, l_max
      DO lx = 0, l
      DO ly = 0, l-lx
         lz = l-lx-ly
         DO ll = 0, l_max
         DO llx = 0, ll
         DO lly = 0, ll-llx
            llz = ll-llx-lly
            S_R(lx, ly, lz) = S_R(lx, ly, lz)+ &
                              h_to_c(llx, lx)*h_to_c(lly, ly)*h_to_c(llz, lz)* &
                              S_R_C(llx, lly, llz)
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      S_R(:, :, :) = S_R(:, :, :)*(alpha/pi)**1.5_dp

   END SUBROUTINE pgf_sum_2c_rspace_3d

! **************************************************************************************************
!> \brief helper method to determine bounds for lattice sums for 3-center integrals in 1 dimension.
!> \param lgth ...
!> \param G_min ...
!> \param R_min ...
!> \param la_max ...
!> \param lb_max ...
!> \param lc_max ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param a_mm ...
!> \param n_sum ...
!> \param G_rads_1 ...
!> \param R_rads_2 ...
!> \param R_rads_3 ...
! **************************************************************************************************
   SUBROUTINE eri_mme_3c_get_bounds(lgth, G_min, R_min, la_max, lb_max, lc_max, &
                                    zeta, zetb, zetc, a_mm, n_sum, G_rads_1, R_rads_2, R_rads_3)
      REAL(KIND=dp)                                      :: lgth
      REAL(KIND=dp), INTENT(IN)                          :: G_min, R_min
      INTEGER, INTENT(IN)                                :: la_max, lb_max, lc_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb, zetc, a_mm
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: n_sum
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: G_rads_1, R_rads_2
      REAL(KIND=dp), DIMENSION(2), INTENT(OUT)           :: R_rads_3

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_3c_get_bounds', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: order_1, order_2
      REAL(KIND=dp)                                      :: alpha, alpha_R, beta, G_res, gamma, R_res

      ! resolution for Gaussian width
      G_res = 0.5_dp*G_min
      R_res = 0.5_dp*R_min

      ! exponents in G space
      alpha = 0.25_dp/zeta
      beta = 0.25_dp/zetb
      gamma = 0.25_dp/zetc+a_mm

      ! Summation radii and number of summands for all lattice summation methods
      ! sum method 1
      G_rads_1(1) = exp_radius(la_max, alpha, lattice_sum_prec, 1.0_dp, epsin=G_res)
      G_rads_1(2) = exp_radius(lb_max, beta, lattice_sum_prec, 1.0_dp, epsin=G_res)
      G_rads_1(3) = exp_radius(lc_max, gamma, lattice_sum_prec, 1.0_dp, epsin=G_res)

      order_1 = MAXLOC(G_rads_1, DIM=1)

      ! note: n_sum calculated here is an upper bound for the actual number of summands
      SELECT CASE (order_1)
      CASE (1)
         n_sum(1) = (FLOOR(2.0_dp*G_rads_1(2)*lgth/twopi)+1)*(FLOOR(2.0_dp*G_rads_1(3)*lgth/twopi)+1)
      CASE (2)
         n_sum(1) = (FLOOR(2.0_dp*G_rads_1(1)*lgth/twopi)+1)*(FLOOR(2.0_dp*G_rads_1(3)*lgth/twopi)+1)
      CASE (3)
         n_sum(1) = (FLOOR(2.0_dp*G_rads_1(1)*lgth/twopi)+1)*(FLOOR(2.0_dp*G_rads_1(2)*lgth/twopi)+1)
      END SELECT

      ! sum method 2
      order_2 = MINLOC(G_rads_1, DIM=1)
      R_rads_2(1) = exp_radius(lb_max+lc_max, 0.25_dp/(beta+gamma), lattice_sum_prec, 1.0_dp, epsin=R_res)
      R_rads_2(2) = exp_radius(lc_max+la_max, 0.25_dp/(alpha+gamma), lattice_sum_prec, 1.0_dp, epsin=R_res)
      R_rads_2(3) = exp_radius(lb_max+la_max, 0.25_dp/(alpha+beta), lattice_sum_prec, 1.0_dp, epsin=R_res)
      n_sum(2) = (FLOOR(2.0_dp*G_rads_1(order_2)*lgth/twopi)+1)*(FLOOR(2.0_dp*R_rads_2(order_2)/lgth)+1)

      ! sum method 3
      alpha_R = 1.0_dp/((zeta+zetb+zetc)/((zeta+zetb)*zetc)+4.0_dp*a_mm)
      R_rads_3(1) = exp_radius(la_max+lb_max, zeta*zetb/(zeta+zetb), lattice_sum_prec, 1.0_dp, epsin=R_res)
      R_rads_3(2) = exp_radius(la_max+lb_max+lc_max, alpha_R, lattice_sum_prec, 1.0_dp, epsin=R_res)
      n_sum(3) = (FLOOR(2.0_dp*R_rads_3(1)/lgth)+1)*(FLOOR(2.0_dp*R_rads_3(2)/lgth)+1)

   END SUBROUTINE eri_mme_3c_get_bounds

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
!> \param G_min ...
!> \param R_min ...
!> \param RA ...
!> \param RB ...
!> \param RC ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param a_mm ...
!> \param lgth ...
!> \param method ...
!> \param method_out ...
!> \param order ...
! **************************************************************************************************
   SUBROUTINE pgf_sum_3c_1d(S_G, G_min, R_min, RA, RB, RC, zeta, zetb, zetc, a_mm, lgth, method, method_out, order)
      REAL(KIND=dp), DIMENSION(0:, 0:, 0:), &
         INTENT(INOUT)                                   :: S_G
      REAL(KIND=dp), INTENT(IN)                          :: G_min, R_min, RA, RB, RC, zeta, zetb, &
                                                            zetc, a_mm, lgth
      INTEGER, INTENT(IN), OPTIONAL                      :: method, order

      CHARACTER(LEN=*), PARAMETER :: routineN = 'pgf_sum_3c_1d', routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: exp1, exp3
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: S_G_c
      INTEGER                                            :: G1c, G2c, G3c, i_pgf, k, k1, k2, &
                                                            k3, l_max, m_max, n_max, R_bound1, &
                                                            R_bound2, sum_method, sum_order, &
                                                            n_sum_predicted, handle
      INTEGER, DIMENSION(3)                              :: G_bounds, n_sum
      REAL(KIND=dp)                                      :: alpha, beta, dG, G, gamma, inv_lgth, R1, &
                                                            R2, R_rad_2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: exp2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: S_G_tmp
      REAL(KIND=dp), DIMENSION(2)                        :: R_rads_3
      REAL(KIND=dp), DIMENSION(3)                        :: G_rads_1, R_rads_2
      INTEGER, INTENT(OUT), OPTIONAL                     :: method_out
      INTEGER                                            :: n_sum_count
      LOGICAL, PARAMETER                                 :: debug = .FALSE.

      n_sum_count = -1
      IF (debug) n_sum_count = 0

      IF (PRESENT(order)) THEN
         sum_order = order
      ELSE
         sum_order = 0
      ENDIF

      IF (PRESENT(method)) THEN
         sum_method = method
      ELSE
         sum_method = 0
      ENDIF

      alpha = 0.25_dp/zeta
      beta = 0.25_dp/zetb
      gamma = 0.25_dp/zetc+a_mm

      inv_lgth = 1.0_dp/lgth
      dG = inv_lgth*twopi

      l_max = UBOUND(S_G, 1)
      m_max = UBOUND(S_G, 2)
      n_max = UBOUND(S_G, 3)

      R1 = RB-RA
      R2 = RC-RB

      CALL eri_mme_3c_get_bounds(lgth, G_min, R_min, l_max, m_max, n_max, zeta, zetb, &
                                 zetc, a_mm, n_sum, G_rads_1, R_rads_2, R_rads_3)

      IF (sum_method .EQ. 0) sum_method = MINLOC(n_sum, DIM=1)
      IF (debug) n_sum_predicted = n_sum(sum_method)

      ALLOCATE (S_G_c(0:l_max, 0:m_max, 0:n_max))
      S_G_c(:, :, :) = 0.0_dp

      IF (.NOT. sum_method .EQ. 3) THEN
         DO i_pgf = 1, 3
            G_bounds(i_pgf) = FLOOR(G_rads_1(i_pgf)/(inv_lgth*twopi))
         ENDDO
         G1c = G_bounds(1)
         G2c = G_bounds(2)
         G3c = G_bounds(3)
      ENDIF

      SELECT CASE (sum_method)
      CASE (1) ! sum_G sum_G' C(G) C(G,G') C(G)
         CALL timeset(routineN//"_GG", handle)

         ! rearrange the sum s.t. G and G' run over the most localized functions
         IF (sum_order .EQ. 0) THEN
            sum_order = MAXLOC(G_bounds, DIM=1)
         ENDIF

         ! Precompute exponentials
         SELECT CASE (sum_order)
         CASE (1)
            G1c = G_bounds(2)+G_bounds(3)
         CASE (2)
            G2c = G_bounds(1)+G_bounds(3)
         CASE (3)
            G3c = G_bounds(1)+G_bounds(2)
         END SELECT

         ALLOCATE (exp1(-G1c:G1c))
         ALLOCATE (exp2(-G2c:G2c))
         ALLOCATE (exp3(-G3c:G3c))

         DO k = -G1c, G1c
            G = k*dG
            exp1(k) = EXP(-alpha*G**2)*EXP(gaussi*G*R1)
         ENDDO

         DO k = -G2c, G2c
            G = k*dG
            exp2(k) = EXP(-beta*G**2)
         ENDDO

         DO k = -G3c, G3c
            G = k*dG
            exp3(k) = EXP(-gamma*G**2)*EXP(gaussi*G*R2)
         ENDDO

         ! Perform double G sum
         SELECT CASE (sum_order)
         CASE (1)
            DO k2 = -G2c, G2c
               DO k3 = -G3c, G3c
                  k1 = k3-k2
                  IF (debug) n_sum_count = n_sum_count+1
                  CALL pgf_product_gspace_3c_1d(S_G_c, k1, k2, k3, dG, l_max, m_max, n_max, G1c, G2c, G3c, exp1, exp2, exp3)
               ENDDO
            ENDDO
         CASE (2)
            DO k1 = -G1c, G1c
               DO k3 = -G3c, G3c
                  k2 = k3-k1
                  IF (debug) n_sum_count = n_sum_count+1
                  CALL pgf_product_gspace_3c_1d(S_G_c, k1, k2, k3, dG, l_max, m_max, n_max, G1c, G2c, G3c, exp1, exp2, exp3)
               ENDDO
            ENDDO
         CASE (3)
            DO k1 = -G1c, G1c
               DO k2 = -G2c, G2c
                  k3 = k1+k2
                  IF (debug) n_sum_count = n_sum_count+1
                  CALL pgf_product_gspace_3c_1d(S_G_c, k1, k2, k3, dG, l_max, m_max, n_max, G1c, G2c, G3c, exp1, exp2, exp3)
               ENDDO
            ENDDO
         END SELECT

         S_G = REAL(S_G_c, KIND=dp)*inv_lgth**2
         CALL timestop(handle)
      CASE (2) ! sum_G sum_R C(G) C(R)
         CALL timeset(routineN//"_GR", handle)
         IF (sum_order .EQ. 0) sum_order = MINLOC(G_bounds, DIM=1) ! 3 cases of rearranging the G sum runs over the most localized function
         R_rad_2 = R_rads_2(sum_order)
         SELECT CASE (sum_order)
         CASE (1)
            ALLOCATE (S_G_tmp(0:l_max, 0:m_max, 0:n_max))
            R_bound1 = CEILING((R2-R_rad_2)/lgth)
            R_bound2 = FLOOR((R2+R_rad_2)/lgth)
            CALL pgf_sum_product_gspace_3c_1d(S_G_tmp, R1, R2, alpha, beta, gamma, lgth, R_bound1, R_bound2, -G1c, G1c, &
                                              n_sum_count)
            S_G = RESHAPE(S_G_tmp, SHAPE(S_G), order=[1, 2, 3])
         CASE (2)
            ALLOCATE (S_G_tmp(0:m_max, 0:l_max, 0:n_max))
            R_bound1 = CEILING(((R1+R2)-R_rad_2)/lgth)
            R_bound2 = FLOOR(((R1+R2)+R_rad_2)/lgth)
            CALL pgf_sum_product_gspace_3c_1d(S_G_tmp, -R1, R1+R2, beta, alpha, gamma, lgth, R_bound1, R_bound2, -G2c, G2c, &
                                              n_sum_count)
            S_G = RESHAPE(S_G_tmp, SHAPE(S_G), order=[2, 1, 3])
         CASE (3)
            ALLOCATE (S_G_tmp(0:n_max, 0:m_max, 0:l_max))
            R_bound1 = CEILING((-R1-R_rad_2)/lgth)
            R_bound2 = FLOOR((-R1+R_rad_2)/lgth)
            CALL pgf_sum_product_gspace_3c_1d(S_G_tmp, -R2, -R1, gamma, beta, alpha, lgth, R_bound1, R_bound2, -G3c, G3c, &
                                              n_sum_count)
            S_G = RESHAPE(S_G_tmp, SHAPE(S_G), order=[3, 2, 1])
         END SELECT

         CALL timestop(handle)
      CASE (3) ! sum_R sum_R' H(R, R')
         CALL timeset(routineN//"_RR", handle)
         CALL pgf_sum_rspace_3c_1d(S_G, RA, RB, RC, zeta, zetb, zetc, a_mm, lgth, R_rads_3, n_sum_count)
         S_G = S_G*pi**(-0.5_dp)*((zeta+zetb)/(zeta*zetb))**(-0.5_dp)
         CALL timestop(handle)
      END SELECT

      IF (debug) THEN
         CPASSERT(n_sum_count .LE. n_sum_predicted)
      ENDIF
      IF (PRESENT(method_out)) method_out = sum_method

   CONTAINS

! **************************************************************************************************
!> \brief Helper routine: compute (-1)^(l+m)*(i)^(l+m+n)*G1^l*exp_G1*G2^m*exp_G2*G3^n*exp_G3
!> \param S_G ...
!> \param k1 ...
!> \param k2 ...
!> \param k3 ...
!> \param dG ...
!> \param l_max ...
!> \param m_max ...
!> \param n_max ...
!> \param G1c ...
!> \param G2c ...
!> \param G3c ...
!> \param exp1 ...
!> \param exp2 ...
!> \param exp3 ...
! **************************************************************************************************
      PURE SUBROUTINE pgf_product_gspace_3c_1d(S_G, k1, k2, k3, dG, l_max, m_max, n_max, G1c, G2c, G3c, exp1, exp2, exp3)
      COMPLEX(KIND=dp), DIMENSION(0:, 0:, 0:), &
         INTENT(INOUT)                                   :: S_G
      INTEGER, INTENT(IN)                                :: k1, k2, k3
      REAL(KIND=dp), INTENT(IN)                          :: dG
      INTEGER, INTENT(IN)                                :: l_max, m_max, n_max, G1c, G2c, G3c
      COMPLEX(KIND=dp), DIMENSION(-G1c:G1c), INTENT(IN)  :: exp1
      REAL(KIND=dp), DIMENSION(-G2c:G2c), INTENT(IN)     :: exp2
      COMPLEX(KIND=dp), DIMENSION(-G3c:G3c), INTENT(IN)  :: exp3

      INTEGER                                            :: l, m, n
      REAL(KIND=dp)                                      :: G1, G1_pow_l, G2, G2_pow_m, G3, G3_pow_n

         G1 = k1*dG
         G2 = k2*dG
         G3 = k3*dG

         G3_pow_n = 1.0_dp
         DO n = 0, n_max
            G2_pow_m = 1.0_dp
            DO m = 0, m_max
               G1_pow_l = 1.0_dp
               DO l = 0, l_max
                  S_G(l, m, n) = S_G(l, m, n)+G1_pow_l*exp1(k1)*G2_pow_m*exp2(k2)*G3_pow_n*exp3(k3)*i_pow(l+m+n)*(-1.0_dp)**(l+m)
                  G1_pow_l = G1_pow_l*G1
               ENDDO
               G2_pow_m = G2_pow_m*G2
            ENDDO
            G3_pow_n = G3_pow_n*G3
         ENDDO

      END SUBROUTINE pgf_product_gspace_3c_1d

! **************************************************************************************************
!> \brief Helper routine: compute (-1)^(l+m) i^(l+m+n)/lgth**2 sum(G) e^(i G R1) C(G, l, alpha) S_R(l, m)
!>        (S_R obtained with pgf_sum_product_rspace_3c_1d)
!> \param S_G ...
!> \param R1 ...
!> \param R2 ...
!> \param alpha ...
!> \param beta ...
!> \param gamma ...
!> \param lgth ...
!> \param R_min ...
!> \param R_c ...
!> \param G_min ...
!> \param G_c ...
!> \param n_sum_count ...
! **************************************************************************************************
      PURE SUBROUTINE pgf_sum_product_gspace_3c_1d(S_G, R1, R2, alpha, beta, gamma, lgth, R_min, R_c, G_min, G_c, n_sum_count)
      REAL(KIND=dp), DIMENSION(0:, 0:, 0:), &
         INTENT(INOUT)                                   :: S_G
      REAL(KIND=dp), INTENT(IN)                          :: R1, R2, alpha, beta, gamma, lgth
      INTEGER, INTENT(IN)                                :: R_min, R_c, G_min, G_c
      INTEGER, INTENT(INOUT), OPTIONAL                   :: n_sum_count

      COMPLEX(KIND=dp)                                   :: exp_tot
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: S_R
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: S_G_c
      INTEGER                                            :: k, l, l_max, m, m_max, n, n_max
      REAL(KIND=dp)                                      :: dG, G, G_pow_l

         dG = 1.0_dp/lgth*twopi
         l_max = UBOUND(S_G, 1)
         m_max = UBOUND(S_G, 2)
         n_max = UBOUND(S_G, 3)
         ALLOCATE (S_G_c(0:l_max, 0:m_max, 0:n_max))

         S_G_c(:, :, :) = 0.0_dp
         ALLOCATE (S_R(0:m_max, 0:n_max))
         DO k = G_min, G_c
            G = k*dG
            exp_tot = EXP(-alpha*G**2)*EXP(gaussi*G*R1)
            CALL pgf_sum_product_rspace_3c_1d(S_R, G, R2, beta, gamma, lgth, R_min, R_c, n_sum_count)
            DO n = 0, n_max
               DO m = 0, m_max
                  G_pow_l = 1.0_dp
                  DO l = 0, l_max
                     S_G_c(l, m, n) = S_G_c(l, m, n)+G_pow_l*exp_tot*i_pow(l+m+n)*(-1.0_dp)**(l+m)*S_R(m, n)
                     G_pow_l = G_pow_l*G
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         S_G = REAL(S_G_c, KIND=dp)/lgth**2

      END SUBROUTINE pgf_sum_product_gspace_3c_1d

! **************************************************************************************************
!> \brief Helper routine: compute sum(R') sum_{t=0}^{l+m} E(t,l,m) (-i)^t C(R'+R,t,alpha) * e^(i beta (R'+R) G)
!>        with alpha = 1/(4(a+b)), beta = a/(a+b).
!>        E(t,l,m) are the coefficients for expanding the product
!>        C(G'- G, l, a) C(G', m, b) into Hermite functions.
!> \param S_R ...
!> \param G ...
!> \param R ...
!> \param alpha ...
!> \param beta ...
!> \param lgth ...
!> \param R_min ...
!> \param R_c ...
!> \param n_sum_count ...
! **************************************************************************************************
      PURE SUBROUTINE pgf_sum_product_rspace_3c_1d(S_R, G, R, alpha, beta, lgth, R_min, R_c, n_sum_count)
      COMPLEX(KIND=dp), DIMENSION(0:, 0:), INTENT(INOUT) :: S_R
      REAL(KIND=dp), INTENT(IN)                          :: G, R, alpha, beta, lgth
      INTEGER, INTENT(IN)                                :: R_min, R_c
      INTEGER, INTENT(INOUT), OPTIONAL                   :: n_sum_count

      COMPLEX(KIND=dp)                                   :: exp_tot
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: S_R_t
      INTEGER                                            :: k, l, l_max, m, m_max, t
      REAL(KIND=dp)                                      :: c1, c2, dG, dR, R_pow_t, Rp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: E

         dR = lgth
         dG = twopi/dR
         l_max = UBOUND(S_R, 1)
         m_max = UBOUND(S_R, 2)
         c1 = 0.25_dp/(alpha+beta)
         c2 = alpha/(alpha+beta)

         ALLOCATE (S_R_t(0:l_max+m_max))

         S_R_t(:) = 0.0_dp
         S_R(:, :) = 0.0_dp

         DO k = R_min, R_c
            Rp = k*dR-R
            exp_tot = EXP(-c1*Rp**2)*EXP(-gaussi*c2*Rp*G)
            R_pow_t = 1.0_dp
            DO t = 0, l_max+m_max
               S_R_t(t) = S_R_t(t)+i_pow(t)*(-1)**t*R_pow_t*exp_tot
               R_pow_t = R_pow_t*Rp
            ENDDO
            IF (n_sum_count .GE. 0) n_sum_count = n_sum_count+1
         ENDDO

         ! H --> CC
         CALL create_gaussian_overlap_dist_to_hermite(l_max, m_max, alpha, beta, G, 0.0_dp, 1, E)
         DO m = 0, m_max
            DO l = 0, l_max
               DO t = 0, l+m
                  S_R(l, m) = S_R(l, m)+E(t, l, m)*S_R_t(t)
               ENDDO
            ENDDO
         ENDDO

         S_R(:, :) = S_R(:, :)*lgth/(twopi)*(pi/(alpha+beta))**0.5_dp

      END SUBROUTINE pgf_sum_product_rspace_3c_1d

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
!> \param R_rads ...
!> \param n_sum_count ...
! **************************************************************************************************
      PURE SUBROUTINE pgf_sum_rspace_3c_1d(S_R, RA, RB, RC, zeta, zetb, zetc, a_mm, lgth, R_rads, n_sum_count)
      REAL(KIND=dp), DIMENSION(0:, 0:, 0:), &
         INTENT(INOUT)                                   :: S_R
      REAL(KIND=dp), INTENT(IN)                          :: RA, RB, RC, zeta, zetb, zetc, a_mm, lgth
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: R_rads
      INTEGER, INTENT(INOUT), OPTIONAL                   :: n_sum_count

      INTEGER                                            :: k1, k2, l, l_max, m, m_max, n, n_max, &
                                                            R1_bound1, R1_bound2, R2_bound1, &
                                                            R2_bound2, t
      REAL(KIND=dp)                                      :: alpha, dR1, dR2, exp_tot, R, R1, R2, &
                                                            R_offset, R_pow_t, R_tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: S_R_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: h_to_c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: E

         dR1 = lgth
         dR2 = lgth
         alpha = 1.0_dp/((zeta+zetb+zetc)/((zeta+zetb)*zetc)+4.0_dp*a_mm)
         l_max = UBOUND(S_R, 1)
         m_max = UBOUND(S_R, 2)
         n_max = UBOUND(S_R, 3)

         ALLOCATE (S_R_t(0:l_max+m_max+n_max))

         S_R(:, :, :) = 0.0_dp

         R_offset = RC-(zeta*RA+zetb*RB)/(zeta+zetb)
         R1_bound1 = CEILING((RA-RB-R_rads(1))/dR1)
         R1_bound2 = FLOOR((RA-RB+R_rads(1))/dR1)

         CALL create_hermite_to_cartesian(alpha, l_max+m_max+n_max, h_to_c)

         DO k1 = R1_bound1, R1_bound2

            S_R_t(:) = 0.0_dp
            R1 = k1*dR1
            R_tmp = R_offset+R1*zeta/(zeta+zetb)
            R2_bound1 = CEILING((-R_tmp-R_rads(2))/dR2)
            R2_bound2 = FLOOR((-R_tmp+R_rads(2))/dR2)
            DO k2 = R2_bound1, R2_bound2
               R2 = k2*dR2
               R = R_tmp+R2
               exp_tot = EXP(-alpha*R**2)
               R_pow_t = 1.0_dp
               DO t = 0, l_max+m_max+n_max
                  S_R_t(t) = S_R_t(t)+R_pow_t*exp_tot
                  R_pow_t = R_pow_t*R
               ENDDO
               IF (n_sum_count .GE. 0) n_sum_count = n_sum_count+1
            ENDDO

            ! C --> H
            S_R_t(:) = MATMUL(TRANSPOSE(h_to_c), S_R_t)*SQRT(alpha/pi)

            ! H --> HH
            CALL create_gaussian_overlap_dist_to_hermite(l_max, m_max, zeta, zetb, RA-R1, RB, 2, E)

            DO n = 0, n_max
               DO m = 0, m_max
                  DO l = 0, l_max
                     DO t = 0, l+m
                        S_R(l, m, n) = S_R(l, m, n)+E(t, l, m)*(-1)**n*S_R_t(t+n)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            DEALLOCATE (E)
         ENDDO
      END SUBROUTINE pgf_sum_rspace_3c_1d

   END SUBROUTINE pgf_sum_3c_1d

END MODULE eri_mme_lattice_summation

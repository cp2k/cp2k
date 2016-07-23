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

   USE eri_mme_gaussian,                ONLY: create_hermite_to_cartesian
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: gaussi,&
                                              pi,&
                                              twopi
#include "./base/base_uses.f90"

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
             pgf_sum_2c_gspace_1d_deltal

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

END MODULE eri_mme_lattice_summation

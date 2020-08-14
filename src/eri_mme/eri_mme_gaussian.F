!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods related to properties of Hermite and Cartesian Gaussian functions.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE eri_mme_gaussian
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: gamma1
   USE minimax_exp,                     ONLY: get_exp_minimax_coeff
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eri_mme_gaussian'

   INTEGER, PARAMETER, PUBLIC :: eri_mme_coulomb = 1, &
                                 eri_mme_yukawa = 2, &
                                 eri_mme_longrange = 3

   PUBLIC :: &
      create_gaussian_overlap_dist_to_hermite, &
      create_hermite_to_cartesian, &
      get_minimax_coeff_v_gspace, &
      hermite_gauss_norm

CONTAINS

! **************************************************************************************************
!> \brief Create matrix to transform between cartesian and hermite gaussian
!>        basis functions.
!> \param zet    exponent
!> \param l_max ...
!> \param h_to_c transformation matrix with dimensions (0:l_max, 0:l_max)
!> \note  is idempotent, so transformation is the same
!>        in both directions.
! **************************************************************************************************
   PURE SUBROUTINE create_hermite_to_cartesian(zet, l_max, h_to_c)
      REAL(KIND=dp), INTENT(IN)                          :: zet
      INTEGER, INTENT(IN)                                :: l_max
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: h_to_c

      INTEGER                                            :: k, l

      ALLOCATE (h_to_c(-1:l_max + 1, 0:l_max))
      h_to_c(:, :) = 0.0_dp
      h_to_c(0, 0) = 1.0_dp
      DO l = 0, l_max - 1
         DO k = 0, l + 1
            h_to_c(k, l + 1) = -(k + 1)*h_to_c(k + 1, l) + 2.0_dp*zet*h_to_c(k - 1, l)
         ENDDO
      ENDDO

   END SUBROUTINE create_hermite_to_cartesian

! **************************************************************************************************
!> \brief Norm of 1d Hermite-Gauss functions
!> \param zet ...
!> \param l ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION hermite_gauss_norm(zet, l) RESULT(norm)
      REAL(KIND=dp), INTENT(IN)                          :: zet
      INTEGER, DIMENSION(3), INTENT(IN)                  :: l
      REAL(KIND=dp)                                      :: norm

      norm = 1.0_dp/SQRT((2.0_dp*zet)**(SUM(l) - 1.5_dp)*(gamma1(l(1))*gamma1(l(2))*gamma1(l(3))))

   END FUNCTION hermite_gauss_norm

! **************************************************************************************************
!> \brief Get minimax coefficient a_i and w_i for approximating
!>        1/G^2 by sum_i w_i exp(-a_i G^2)
!> \param n_minimax   Number of minimax terms
!> \param cutoff      Plane Wave cutoff
!> \param G_min       Minimum absolute value of G
!> \param minimax_aw  Minimax coefficients a_i, w_i
!> \param potential   potential to use. Accepts the following values:
!>                    1: coulomb potential V(r)=1/r
!>                    2: yukawa potential V(r)=e(-a*r)/r
!>                    3: long-range coulomb erf(a*r)/r
!> \param pot_par     potential parameter a for yukawa V(r)=e(-a*r)/r or long-range coulomb V(r)=erf(a*r)/r
!> \param err_minimax Maximum error MAX (|1/G^2-\sum_i w_i exp(-a_i G^2)|)
! **************************************************************************************************
   SUBROUTINE get_minimax_coeff_v_gspace(n_minimax, cutoff, G_min, minimax_aw, potential, pot_par, err_minimax)
      INTEGER, INTENT(IN)                                :: n_minimax
      REAL(KIND=dp), INTENT(IN)                          :: cutoff, G_min
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: minimax_aw
      INTEGER, INTENT(IN), OPTIONAL                      :: potential
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: pot_par
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: err_minimax

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_minimax_coeff_v_gspace'

      INTEGER                                            :: potential_prv
      REAL(KIND=dp)                                      :: dG, G_max, minimax_Rc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: a, w

      IF (PRESENT(potential)) THEN
         potential_prv = potential
      ELSE
         potential_prv = eri_mme_coulomb
      ENDIF

      IF (potential_prv > 3) THEN
         CPABORT("unknown potential")
      ENDIF

      IF ((potential_prv >= 2) .AND. .NOT. PRESENT(pot_par)) THEN
         CPABORT("potential parameter pot_par required for yukawa or long-range Coulomb")
      ENDIF

      dG = 1.0E-3 ! Resolution in G to determine error of minimax approximation

      ! Note: G_c = SQRT(2*cutoff) cutoff in 1 cartesian direction
      ! G_max = SQRT(3*G_c**2) maximum absolute value of G vector
      ! Minimax approx. needs to be valid in range [G_min, G_max]

      ! 1) compute minimax coefficients

      G_max = SQRT(3.0_dp*2.0_dp*cutoff)
      CPASSERT(G_max .GT. G_min)
      IF (potential_prv == eri_mme_coulomb .OR. potential_prv == eri_mme_longrange) THEN
         minimax_Rc = (G_max/G_min)**2
      ELSEIF (potential_prv == eri_mme_yukawa) THEN
         minimax_Rc = (G_max**2 + pot_par**2)/(G_min**2 + pot_par**2)
      ENDIF

      CALL get_exp_minimax_coeff(n_minimax, minimax_Rc, minimax_aw, err_minimax)

      ALLOCATE (a(n_minimax)); ALLOCATE (w(n_minimax))
      a(:) = minimax_aw(:n_minimax)
      w(:) = minimax_aw(n_minimax + 1:)
      SELECT CASE (potential_prv)
         ! Scale minimax coefficients to incorporate different Fourier transforms
      CASE (eri_mme_coulomb)
         ! FT = 1/G**2
         a(:) = a/G_min**2
         w(:) = w/G_min**2
      CASE (eri_mme_yukawa)
         ! FT = 1/(G**2 + pot_par**2)
         w(:) = w*EXP((-a*pot_par**2)/(G_min**2 + pot_par**2))/(G_min**2 + pot_par**2)
         a(:) = a/(G_min**2 + pot_par**2)
      CASE (eri_mme_longrange)
         ! FT = exp(-(G/pot_par)**2)/G**2
         ! approximating 1/G**2 as for Coulomb:
         a(:) = a/G_min**2
         w(:) = w/G_min**2
         ! incorporate exponential factor:
         a(:) = a + 1.0_dp/pot_par**2
      END SELECT
      minimax_aw = [a(:), w(:)]

      IF (PRESENT(err_minimax)) THEN
         IF (potential_prv == eri_mme_coulomb) THEN
            err_minimax = err_minimax/G_min**2
         ELSEIF (potential_prv == eri_mme_yukawa) THEN
            err_minimax = err_minimax/(G_min**2 + pot_par**2)
         ELSEIF (potential_prv == eri_mme_longrange) THEN
            err_minimax = err_minimax/G_min**2 ! approx. of Coulomb
            err_minimax = err_minimax*EXP(-G_min**2/pot_par**2) ! exponential factor
         ENDIF
      ENDIF

   END SUBROUTINE get_minimax_coeff_v_gspace

! **************************************************************************************************
!> \brief Expand 1d product of cartesian (or hermite) gaussians into single hermite gaussians:
!>        Find E_t^{lm} s.t.
!>        F(l, a, r-R1) * F(m, b, r-R2) = sum_{t=0}^{l+m} E_t^{lm} H(t, p, r-R_P)
!>        with p = a + b, R_P = (a*R1 + b*R2)/p. The function F can be either Cartesian
!>        Gaussian or Hermite Gaussian.
!> \param l ...
!> \param m ...
!> \param a ...
!> \param b ...
!> \param R1 ...
!> \param R2 ...
!> \param H_or_C_product 1: cartesian product, 2: hermite product
!> \param E ...
! **************************************************************************************************
   PURE SUBROUTINE create_gaussian_overlap_dist_to_hermite(l, m, a, b, R1, R2, H_or_C_product, E)
      INTEGER, INTENT(IN)                                :: l, m
      REAL(KIND=dp), INTENT(IN)                          :: a, b, R1, R2
      INTEGER, INTENT(IN)                                :: H_or_C_product
      REAL(KIND=dp), DIMENSION(-1:l+m+1, -1:l, -1:m), &
         INTENT(OUT)                                     :: E

      INTEGER                                            :: ll, mm, t
      REAL(KIND=dp)                                      :: c1, c2, c3

      E(:, :, :) = 0.0_dp
      E(0, 0, 0) = EXP(-a*b/(a + b)*(R1 - R2)**2) ! cost: exp_w flops

      c1 = 0.5_dp/(a + b)
      c2 = (b/(a + b))*(R2 - R1)
      c3 = (a/(a + b))*(R1 - R2)

      IF (H_or_C_product .EQ. 1) THEN ! Cartesian overlap dist
         DO mm = 0, m
            DO ll = 0, l
               DO t = 0, ll + mm + 1
                  IF (ll .LT. l) THEN
                     E(t, ll + 1, mm) = c1*E(t - 1, ll, mm) + & ! cost: 8 flops
                                        c2*E(t, ll, mm) + &
                                        (t + 1)*E(t + 1, ll, mm)
                  ENDIF
                  IF (mm .LT. m) THEN
                     E(t, ll, mm + 1) = c1*E(t - 1, ll, mm) + & ! cost: 8 flops
                                        c3*E(t, ll, mm) + &
                                        (t + 1)*E(t + 1, ll, mm)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ELSE ! Hermite overlap dist
         DO mm = 0, m
            DO ll = 0, l
               DO t = 0, ll + mm + 1
                  IF (ll .LT. l) THEN
                     E(t, ll + 1, mm) = a*(2*c1*E(t - 1, ll, mm) + & ! cost: 16 flops
                                           2*c2*E(t, ll, mm) + &
                                           2*(t + 1)*E(t + 1, ll, mm) - &
                                           2*ll*E(t, ll - 1, mm))
                  ENDIF
                  IF (mm .LT. m) THEN
                     E(t, ll, mm + 1) = b*(2*c1*E(t - 1, ll, mm) + & ! cost: 16 flops
                                           2*c3*E(t, ll, mm) + &
                                           2*(t + 1)*E(t + 1, ll, mm) - &
                                           2*mm*E(t, ll, mm - 1))

                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF

   END SUBROUTINE create_gaussian_overlap_dist_to_hermite
END MODULE eri_mme_gaussian

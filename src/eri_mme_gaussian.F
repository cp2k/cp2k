!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Low-level routines for Minimax-Ewald (MME) method for electron repulsion
!>        integrals.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE eri_mme_gaussian
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: gamma1,&
                                              maxfac
   USE minimax,                         ONLY: get_minimax_coeff
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eri_mme_gaussian'

   PUBLIC :: create_hermite_to_cartesian, hermite_gauss_norm, get_minimax_coeff_v_gspace

CONTAINS

! **************************************************************************************************
!> \brief Create matrix to transform between cartesian and hermite gaussian
!>        basis functions.
!> \param zet    exponent
!> \param h_to_c transformation matrix with dimensions (0:l_max, 0:l_max)
!> \note  is idempotent, so transformation is the same
!>        in both directions.
! **************************************************************************************************
   SUBROUTINE create_hermite_to_cartesian(zet, h_to_c)
      REAL(KIND=dp), INTENT(IN)                          :: zet
      REAL(KIND=dp), DIMENSION(0:, 0:), INTENT(INOUT)    :: h_to_c

      INTEGER                                            :: k, l, l_max

      l_max = UBOUND(h_to_c, 1)

      CPASSERT((UBOUND(h_to_c, 2) .EQ. l_max))

      h_to_c(:, :) = 0.0_dp
      h_to_c(0, 0) = 1.0_dp
      DO l = 0, l_max-1
         DO k = 0, l+1
            IF (k .EQ. 0) THEN
               h_to_c(l+1, k) = (k+1)*h_to_c(l, k+1)
            ELSE
               IF (k .LT. l) h_to_c(l+1, k) = (k+1)*h_to_c(l, k+1)-2.0_dp*zet*h_to_c(l, k-1)
               IF (k .GE. l) h_to_c(l+1, k) = -2.0_dp*zet*h_to_c(l, k-1)
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE create_hermite_to_cartesian

! **************************************************************************************************
!> \brief Norm of 1d Hermite-Gauss functions
!> \param zet ...
!> \param l ...
!> \retval norm ...
! **************************************************************************************************
   FUNCTION hermite_gauss_norm(zet, l) RESULT(norm)
      REAL(KIND=dp), INTENT(IN)                          :: zet
      INTEGER, DIMENSION(3), INTENT(IN)                  :: l
      REAL(KIND=dp)                                      :: norm

      CPASSERT(MAXVAL(l) .LE. maxfac)
      norm = 1.0_dp/SQRT((2.0_dp*zet)**(SUM(l)-1.5_dp)*(gamma1(l(1))*gamma1(l(2))*gamma1(l(3))))

   END FUNCTION hermite_gauss_norm

! **************************************************************************************************
!> \brief Get minimax coefficient a_i and w_i for approximating
!>        1/G^2 by sum_i w_i exp(-a_i G^2)
!> \param n_minimax   Number of minimax terms
!> \param cutoff      Plane Wave cutoff
!> \param G_min       Minimum absolute value of G
!> \param minimax_aw  Minimax coefficients a_i, w_i
!> \param err_minimax Maximum error MAX (|1/G^2-\sum_i w_i exp(-a_i G^2)|)
! **************************************************************************************************
   SUBROUTINE get_minimax_coeff_v_gspace(n_minimax, cutoff, G_min, minimax_aw, err_minimax)
      INTEGER, INTENT(IN)                                :: n_minimax
      REAL(KIND=dp), INTENT(INOUT)                       :: cutoff
      REAL(KIND=dp), INTENT(IN)                          :: G_min
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: minimax_aw
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: err_minimax

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_minimax_coeff_v_gspace', &
         routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: dG, G_max, minimax_Rc

      dG = 1.0E-3 ! Resolution in G to determine error of minimax approximation
      !min_cutoff = .FALSE.

      ! Note: G_c = SQRT(2*cutoff) cutoff in 1 cartesian direction
      ! G_max = SQRT(3*G_c**2) maximum absolute value of G vector
      ! Minimax approx. needs to be valid in range [G_min, G_max]

      ! 1) compute minimax coefficients
      G_max = SQRT(3.0_dp*2.0_dp*cutoff)
      minimax_Rc = (G_max/G_min)**2

      G_max = SQRT(3.0_dp*2.0_dp*cutoff)
      CPASSERT(G_max .GT. G_min)

      CALL get_minimax_coeff(n_minimax, minimax_Rc, minimax_aw, err_minimax)

      minimax_aw = minimax_aw/G_min**2
      IF (PRESENT(err_minimax)) err_minimax = err_minimax/G_min**2

   END SUBROUTINE get_minimax_coeff_v_gspace

END MODULE eri_mme_gaussian

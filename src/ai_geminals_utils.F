!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of Coulomb integrals over Correlated Cartesian Gaussian
!>        Geminals (electron repulsion integrals, ERIs).
!>        Utility routines
!> \par Literature
!>      S.L. Saito and Y. Suzuki, J. Chem. Phys. 114 1115 (2001)
!> \par History
!>      none
!> \author Juerg Hutter (07.2009)
! *****************************************************************************
MODULE ai_geminals_utils
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fac
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_geminals_utils'

  PRIVATE

  PUBLIC  :: matv, acoff

! *****************************************************************************

  CONTAINS

! *****************************************************************************
  FUNCTION matv(a,b) RESULT(c)
    REAL(dp), DIMENSION(2, 2)                :: a
    REAL(dp), DIMENSION(6)                   :: b, c

     c(1) = a(1,1)*b(1) + a(1,2)*b(4)
     c(2) = a(1,1)*b(2) + a(1,2)*b(5)
     c(3) = a(1,1)*b(3) + a(1,2)*b(6)
     c(4) = a(1,2)*b(1) + a(2,2)*b(4)
     c(5) = a(1,2)*b(2) + a(2,2)*b(5)
     c(6) = a(1,2)*b(3) + a(2,2)*b(6)

  END FUNCTION matv

! *****************************************************************************

  SUBROUTINE acoff(ac,nmax)
    REAL(dp), DIMENSION(0:, 0:, 0:)          :: ac
    INTEGER, INTENT(IN)                      :: nmax

    INTEGER                                  :: i, j, k, q

    ac = 0._dp
    DO i=0,nmax
      DO j=0,nmax
        DO k=0,i+j
          DO q=MAX(0,k-j),MIN(i,k)
            ac(i,j,k) = ac(i,j,k) + (-1)**(k-q)*choose(i,q)*choose(j,k-q)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE acoff

! *****************************************************************************

  FUNCTION choose(n,k) RESULT(c)
    INTEGER, INTENT(IN)                      :: n, k
    REAL(dp)                                 :: c

     IF(n==0) THEN
       c = 1._dp
     ELSE
       c = fac(n)/fac(k)/fac(n-k)
     END IF
  END FUNCTION choose

! *****************************************************************************

END MODULE ai_geminals_utils

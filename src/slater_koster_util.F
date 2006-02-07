!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
!!****** cp2k/slater_koster_util [1.0] *
!!
!!   NAME
!!     slater_koster_util
!!
!!   FUNCTION
!!     Some utility functions for the Slater-Koster Modules
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
!
    MODULE slater_koster_util
!------------------------------------------------------------------------------!
  USE kinds,                           ONLY: dp
!

      IMPLICIT NONE
!
      PRIVATE
      PUBLIC :: sph, dsph, out_prod, out_dprod, dpro
!
!!*****
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
      SUBROUTINE sph(l,r,gsl)
! Spherical Harmonics up to f functions
! There is a factor sqrt(4*PI)/sqrt(2*l+1) omitted
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), INTENT(IN)                :: r(3)
    REAL(KIND=dp), INTENT(OUT)               :: gsl(0:6)

    REAL(KIND=dp), PARAMETER :: sq15 = 3.872983346207417_dp , &
      sq3 = 1.732050807568877_dp , sq38 = 0.612372435695794_dp , &
      sq58 = 0.790569415042095_dp 

! sqrt(3.0_dp)
! sqrt(15.0_dp)
! sqrt(0.375_dp)
! sqrt(0.625_dp)

        SELECT CASE (l)
        CASE (0)
          gsl(0) = 1.0_dp
        CASE (1)
          gsl(0) = r(3)
          gsl(1) = r(1)
          gsl(2) = r(2)
        CASE (2)
          gsl(0) = 0.5_dp*(3.0_dp*r(3)*r(3)-1.0_dp)
          gsl(1) = sq3*r(1)*r(3)
          gsl(2) = sq3*r(2)*r(3)
          gsl(3) = 0.5_dp*sq3*(r(1)*r(1)-r(2)*r(2))
          gsl(4) = sq3*r(1)*r(2)
        CASE (3)
          gsl(0) = 0.5_dp*(5.0_dp*r(3)*r(3)-3.0_dp)*r(3)
          gsl(1) = sq38*(5.0_dp*r(3)*r(3)-1.0_dp)*r(1)
          gsl(2) = sq38*(5.0_dp*r(3)*r(3)-1.0_dp)*r(2)
          gsl(3) = 0.5_dp*sq15*(r(1)*r(1)-r(2)*r(2))*r(3)
          gsl(4) = sq15*r(1)*r(2)*r(3)
          gsl(5) = sq58*(r(1)*r(1)-3.0_dp*r(2)*r(2))*r(1)
          gsl(6) = sq58*(3.0_dp*r(1)*r(1)-r(2)*r(2))*r(2)
        END SELECT

      END SUBROUTINE sph
!------------------------------------------------------------------------------!
      SUBROUTINE dsph(l,r,dgsl)
! Derivatives of Spherical Harmonics up to f functions
! There is a factor sqrt(4*PI)/sqrt(2*l+1)/R omitted
    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), INTENT(IN)                :: r(3)
    REAL(KIND=dp), INTENT(OUT)               :: dgsl(0:6,1:3)

    REAL(KIND=dp), PARAMETER :: sq15 = 3.872983346207417_dp , &
      sq3 = 1.732050807568877_dp , sq38 = 0.612372435695794_dp , &
      sq58 = 0.790569415042095_dp 

    REAL(KIND=dp)                            :: dc(0:6)

! sqrt(3.0_dp)
! sqrt(15.0_dp)
! sqrt(0.375_dp)
! sqrt(0.625_dp)

        SELECT CASE (l)
        CASE (0)
          dgsl(0,:) = 0.0_dp
        CASE (1)
          dgsl(0:2,1:3) = 0.0_dp
          dgsl(1,1) = 1.0_dp
          dgsl(2,2) = 1.0_dp
          dgsl(0,3) = 1.0_dp
        CASE (2)
          dgsl(0,1) = 0.0_dp
          dgsl(1,1) = sq3*r(3)
          dgsl(2,1) = 0.0_dp
          dgsl(3,1) = sq3*r(1)
          dgsl(4,1) = sq3*r(2)
          dgsl(0,2) = 0.0_dp
          dgsl(1,2) = 0.0_dp
          dgsl(2,2) = sq3*r(3)
          dgsl(3,2) = -sq3*r(2)
          dgsl(4,2) = sq3*r(1)
          dgsl(0,3) = 3.0_dp*r(3)
          dgsl(1,3) = sq3*r(1)
          dgsl(2,3) = sq3*r(2)
          dgsl(3,3) = 0.0_dp
          dgsl(4,3) = 0.0_dp
        CASE (3)
          dgsl(0,1) = 0.0_dp
          dgsl(1,1) = sq38*(5.0_dp*r(3)*r(3)-1.0_dp)
          dgsl(2,1) = 0.0_dp
          dgsl(3,1) = sq15*r(1)*r(3)
          dgsl(4,1) = sq15*r(2)*r(3)
          dgsl(5,1) = 3.0_dp*sq58*(r(1)*r(1)-r(2)*r(2))
          dgsl(6,1) = 6.0_dp*sq58*r(1)*r(2)
          dgsl(0,2) = 0.0_dp
          dgsl(1,2) = 0.0_dp
          dgsl(2,2) = sq38*(5.0_dp*r(3)*r(3)-1.0_dp)
          dgsl(3,2) = -sq15*r(2)*r(3)
          dgsl(4,2) = sq15*r(1)*r(3)
          dgsl(5,2) = -6.0_dp*sq58*r(2)*r(1)
          dgsl(6,2) = 3.0_dp*sq58*(r(1)*r(1)-r(2)*r(2))
          dgsl(0,3) = 7.5_dp*r(3)*r(3)-1.5_dp
          dgsl(1,3) = 10.0_dp*sq38*r(3)*r(1)
          dgsl(2,3) = 10.0_dp*sq38*r(3)*r(2)
          dgsl(3,3) = 0.5_dp*sq15*(r(1)*r(1)-r(2)*r(2))
          dgsl(4,3) = sq15*r(1)*r(2)
          dgsl(5,3) = 0.0_dp
          dgsl(6,3) = 0.0_dp
        END SELECT
        dc(0:2*l) = r(1)*dgsl(0:2*l,1) + r(2)*dgsl(0:2*l,2) + r(3)*dgsl(0:2*l,3)
        dgsl(0:2*l,1) = dgsl(0:2*l,1) - r(1)*dc(0:2*l)
        dgsl(0:2*l,2) = dgsl(0:2*l,2) - r(2)*dc(0:2*l)
        dgsl(0:2*l,3) = dgsl(0:2*l,3) - r(3)*dc(0:2*l)

      END SUBROUTINE dsph
!------------------------------------------------------------------------------!
      SUBROUTINE out_prod(mat,v1,v2)
    REAL(KIND=dp), INTENT(OUT)               :: mat(:,:)
    REAL(KIND=dp), INTENT(IN)                :: v1(:), v2(:)

    INTEGER                                  :: i, j, n1, n2

        n1 = SIZE(v1)
        n2 = SIZE(v2)
        DO j = 1, n2
          DO i = 1, n1
            mat(i,j) = v1(i)*v2(j)
          END DO
        END DO
      END SUBROUTINE out_prod
!------------------------------------------------------------------------------!
      SUBROUTINE out_dprod(dmat,v1,v2,dv1,dv2)
    REAL(KIND=dp), INTENT(OUT)               :: dmat(:,:)
    REAL(KIND=dp), INTENT(IN)                :: v1(:), v2(:), dv1(:), dv2(:)

    INTEGER                                  :: i, j, n1, n2

        n1 = SIZE(v1)
        n2 = SIZE(v2)
        DO j = 1, n2
          DO i = 1, n1
            dmat(i,j) = dv1(i)*v2(j) + v1(i)*dv2(j)
          END DO
        END DO
      END SUBROUTINE out_dprod
!------------------------------------------------------------------------------!
      FUNCTION dpro(g1,g2,l1,l2)
    REAL(KIND=dp), INTENT(IN)                :: g1(0:6,0:6), g2(0:6,0:6)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(KIND=dp)                            :: dpro

    INTEGER                                  :: i, j

        dpro = 0.0_dp
        DO i = 0, 2*l2
          DO j = 0, 2*l1
            dpro = dpro + g1(j,i)*g2(j,i)
          END DO
        END DO

      END FUNCTION dpro
!------------------------------------------------------------------------------!
    END MODULE slater_koster_util
!------------------------------------------------------------------------------!

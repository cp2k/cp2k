!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
!!****** cp2k/slater_koster_matr [1.0] *
!!
!!   NAME
!!     slater_koster_matr
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

    MODULE slater_koster_matr
!------------------------------------------------------------------------------!
  USE kinds,                           ONLY: dp
  USE slater_koster_util,              ONLY: dsph,&
                                             out_dprod,&
                                             out_prod,&
                                             sph
!

      IMPLICIT NONE
!
      PRIVATE
      PUBLIC :: gmat, dgmat
!
      REAL(KIND=dp), PARAMETER :: zero = 0.0_dp
      REAL(KIND=dp), PARAMETER :: one = 1.0_dp

!!*****
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
      SUBROUTINE gmat(l1,l2,r,gmu,dll)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(KIND=dp), INTENT(IN)                :: r(3)
    REAL(KIND=dp), INTENT(OUT)               :: gmu(0:6,0:6,0:3), &
                                                dll(0:6,0:6,0:6)

    INTEGER                                  :: i, m1, m2, mo, mu
    REAL(KIND=dp)                            :: gsl(0:6,0:3)

        mu = MIN(l1,l2)
        mo = MAX(l1,l2)
        m1 = 2*l1
        m2 = 2*l2
        gsl(0:6,0:3) = zero
        DO i = 0, mo
          CALL sph(i,r,gsl(0:6,i))
        END DO
        dll(0:6,0:6,0) = zero
        DO i = 0, 2*mo
          dll(i,i,0) = 1.0_dp
        END DO

        IF (l1==0) THEN
!s-l
          gmu(0,0:m2,0) = gsl(0:m2,l2)
        ELSE IF (l2==0) THEN
!l-s
          gmu(0:m1,0,0) = gsl(0:m1,l1)
        ELSE
          SELECT CASE (l1)
          CASE (1)
            SELECT CASE (l2)
            CASE (1)
!p-p
              CALL out_prod(dll(0:2,0:2,2),gsl(0:2,1),gsl(0:2,1))
              gmu(0:2,0:2,0) = dll(0:2,0:2,2)
              gmu(0:2,0:2,1) = -dll(0:2,0:2,2) + dll(0:2,0:2,0)
            CASE (2)
!p-d
              CALL out_prod(dll(0:2,0:4,3),gsl(0:2,1),gsl(0:4,2))
              CALL set_mat(dll(0:2,0:4,1),r,1)
              gmu(0:2,0:4,0) = dll(0:2,0:4,3)
              gmu(0:2,0:4,1) = -2.0_dp/SQRT(3.0_dp)*dll(0:2,0:4,3) + &
                dll(0:2,0:4,1)
            CASE (3)
!p-f
              CALL out_prod(dll(0:2,0:6,4),gsl(0:2,1),gsl(0:6,3))
              CALL set_mat(dll(0:2,0:6,2),r,2)
              gmu(0:2,0:6,0) = dll(0:2,0:6,4)
              gmu(0:2,0:6,1) = -SQRT(1.5_dp)*dll(0:2,0:6,4) + &
                0.5_dp*dll(0:2,0:6,2)
            END SELECT
          CASE (2)
            SELECT CASE (l2)
            CASE (1)
!d-p
              CALL out_prod(dll(0:4,0:2,3),gsl(0:4,2),gsl(0:2,1))
              CALL set_mat(dll(0:4,0:2,1),r,1)
              gmu(0:4,0:2,0) = dll(0:4,0:2,3)
              gmu(0:4,0:2,1) = -2.0_dp/SQRT(3.0_dp)*dll(0:4,0:2,3) + &
                dll(0:4,0:2,1)
            CASE (2)
!d-d
              CALL out_prod(dll(0:4,0:4,4),gsl(0:4,2),gsl(0:4,2))
              CALL set_mat(dll(0:4,0:4,2),r,2)
              gmu(0:4,0:4,0) = dll(0:4,0:4,4)
              gmu(0:4,0:4,1) = -4.0_dp/3.0_dp*dll(0:4,0:4,4) + &
                dll(0:4,0:4,2) + dll(0:4,0:4,0)
              gmu(0:4,0:4,2) = 1.0_dp/3.0_dp*dll(0:4,0:4,4) - dll(0:4,0:4,2)
            CASE (3)
!d-f
              CALL out_prod(dll(0:4,0:6,5),gsl(0:4,2),gsl(0:6,3))
              CALL set_mat(dll(0:4,0:6,3),r,3)
              CALL set_mat(dll(0:4,0:6,1),r,1)
              gmu(0:4,0:6,0) = dll(0:4,0:6,5)
              gmu(0:4,0:6,1) = -SQRT(2.0_dp)*dll(0:4,0:6,5) + &
                0.5_dp*dll(0:4,0:6,3) + 0.5_dp*dll(0:4,0:6,1)
              gmu(0:4,0:6,2) = SQRT(0.2_dp)*dll(0:4,0:6,5) - &
                SQRT(0.4_dp)*dll(0:4,0:6,3)
            END SELECT
          CASE (3)
            SELECT CASE (l2)
            CASE (1)
!f-p
              CALL out_prod(dll(0:6,0:2,4),gsl(0:6,3),gsl(0:2,1))
              CALL set_mat(dll(0:6,0:2,2),r,2)
              gmu(0:6,0:2,0) = dll(0:6,0:2,4)
              gmu(0:6,0:2,1) = -SQRT(1.5_dp)*dll(0:6,0:2,4) + &
                0.5_dp*dll(0:6,0:2,2)
            CASE (2)
!f-d
              CALL out_prod(dll(0:6,0:4,5),gsl(0:6,3),gsl(0:4,2))
              CALL set_mat(dll(0:6,0:4,3),r,3)
              CALL set_mat(dll(0:6,0:4,1),r,1)
              gmu(0:6,0:4,0) = dll(0:6,0:4,5)
              gmu(0:6,0:4,1) = -SQRT(2.0_dp)*dll(0:6,0:4,5) + &
                0.5_dp*dll(0:6,0:4,3) + 0.5_dp*dll(0:6,0:4,1)
              gmu(0:6,0:4,2) = SQRT(0.2_dp)*dll(0:6,0:4,5) - &
                SQRT(0.4_dp)*dll(0:6,0:4,3)
            CASE (3)
!f-f
              CALL out_prod(dll(0:6,0:6,6),gsl(0:6,3),gsl(0:6,3))
              CALL set_mat(dll(0:6,0:6,4),r,4)
              CALL set_mat(dll(0:6,0:6,2),r,2)
              gmu(0:6,0:6,0) = dll(0:6,0:6,6)
              gmu(0:6,0:6,1) = -1.5_dp*dll(0:6,0:6,6) + &
                0.625_dp*dll(0:6,0:6,4) + 0.625_dp*dll(0:6,0:6,2) + &
                0.625_dp*dll(0:6,0:6,0)
              gmu(0:6,0:6,2) = 0.6_dp*dll(0:6,0:6,6) - dll(0:6,0:6,4)
              gmu(0:6,0:6,3) = -0.1_dp*dll(0:6,0:6,6) + &
                0.375_dp*dll(0:6,0:6,4) - 0.625_dp*dll(0:6,0:6,2) + &
                0.375_dp*dll(0:6,0:6,0)
            END SELECT
          END SELECT
        END IF

      END SUBROUTINE gmat
!------------------------------------------------------------------------------!
      SUBROUTINE dgmat(l1,l2,r,dgmu,ddll)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(KIND=dp), INTENT(IN)                :: r(3)
    REAL(KIND=dp), INTENT(OUT)               :: dgmu(0:6,0:6,0:3,1:3), &
                                                ddll(0:6,0:6,0:6,1:3)

    INTEGER                                  :: i, m1, m2, mo, mu
    REAL(KIND=dp)                            :: dgsl(0:6,0:3,1:3), &
                                                gsl(0:6,0:3)

        mu = MIN(l1,l2)
        mo = MAX(l1,l2)
        m1 = 2*l1
        m2 = 2*l2
        gsl(0:6,0:3) = zero
        dgsl(0:6,0:3,1:3) = zero
        DO i = 0, mo
          CALL sph(i,r,gsl(0:6,i))
          CALL dsph(i,r,dgsl(0:6,i,1:3))
        END DO

        IF (l1==0) THEN
!s-l
          dgmu(0,0:m2,0,1:3) = dgsl(0:m2,l2,1:3)
        ELSE IF (l2==0) THEN
!l-s
          dgmu(0:m1,0,0,1:3) = dgsl(0:m1,l1,1:3)
        ELSE
          SELECT CASE (l1)
          CASE (1)
            SELECT CASE (l2)
            CASE (1)
!p-p
              CALL out_dprod(ddll(0:2,0:2,2,1),gsl(0:2,1),gsl(0:2,1), &
                             dgsl(0:2,1,1),dgsl(0:2,1,1))
              CALL out_dprod(ddll(0:2,0:2,2,2),gsl(0:2,1),gsl(0:2,1), &
                             dgsl(0:2,1,2),dgsl(0:2,1,2))
              CALL out_dprod(ddll(0:2,0:2,2,3),gsl(0:2,1),gsl(0:2,1), &
                             dgsl(0:2,1,3),dgsl(0:2,1,3))
              dgmu(0:2,0:2,0,1:3) = ddll(0:2,0:2,2,1:3)
              dgmu(0:2,0:2,1,1:3) = -ddll(0:2,0:2,2,1:3)
            CASE (2)
!p-d
              CALL out_dprod(ddll(0:2,0:4,3,1),gsl(0:2,1),gsl(0:4,2), &
                            dgsl(0:2,1,1),dgsl(0:4,2,1))
              CALL out_dprod(ddll(0:2,0:4,3,2),gsl(0:2,1),gsl(0:4,2), &
                            dgsl(0:2,1,2),dgsl(0:4,2,2))
              CALL out_dprod(ddll(0:2,0:4,3,3),gsl(0:2,1),gsl(0:4,2), &
                            dgsl(0:2,1,3),dgsl(0:4,2,3))
              CALL set_dmat(ddll(0:2,0:4,1,1:3),r,1)
              dgmu(0:2,0:4,0,1:3) = ddll(0:2,0:4,3,1:3)
              dgmu(0:2,0:4,1,1:3) = -2.0_dp/SQRT(3.0_dp)*ddll(0:2,0:4,3,1:3) + &
                ddll(0:2,0:4,1,1:3)
            CASE (3)
!p-f
              CALL out_dprod(ddll(0:2,0:6,4,1),gsl(0:2,1),gsl(0:6,3), &
                             dgsl(0:2,1,1),dgsl(0:6,3,1))
              CALL out_dprod(ddll(0:2,0:6,4,2),gsl(0:2,1),gsl(0:6,3), &
                             dgsl(0:2,1,2),dgsl(0:6,3,2))
              CALL out_dprod(ddll(0:2,0:6,4,3),gsl(0:2,1),gsl(0:6,3), &
                             dgsl(0:2,1,3),dgsl(0:6,3,3))
              CALL set_dmat(ddll(0:2,0:6,2,1:3),r,2)
              dgmu(0:2,0:6,0,1:3) = ddll(0:2,0:6,4,1:3)
              dgmu(0:2,0:6,1,1:3) = -SQRT(1.5_dp)*ddll(0:2,0:6,4,1:3) + &
                0.5_dp*ddll(0:2,0:6,2,1:3)
            END SELECT
          CASE (2)
            SELECT CASE (l2)
            CASE (1)
!d-p
              CALL out_dprod(ddll(0:4,0:2,3,1),gsl(0:4,2),gsl(0:2,1), &
                             dgsl(0:4,2,1),dgsl(0:2,1,1))
              CALL out_dprod(ddll(0:4,0:2,3,2),gsl(0:4,2),gsl(0:2,1), &
                             dgsl(0:4,2,2),dgsl(0:2,1,2))
              CALL out_dprod(ddll(0:4,0:2,3,3),gsl(0:4,2),gsl(0:2,1), &
                             dgsl(0:4,2,3),dgsl(0:2,1,3))
              CALL set_dmat(ddll(0:4,0:2,1,1:3),r,1)
              dgmu(0:4,0:2,0,1:3) = ddll(0:4,0:2,3,1:3)
              dgmu(0:4,0:2,1,1:3) = -2.0_dp/SQRT(3.0_dp)*ddll(0:4,0:2,3,1:3) + &
                ddll(0:4,0:2,1,1:3)
            CASE (2)
!d-d
              CALL out_dprod(ddll(0:4,0:4,4,1),gsl(0:4,2),gsl(0:4,2), &
                             dgsl(0:4,2,1),dgsl(0:4,2,1))
              CALL out_dprod(ddll(0:4,0:4,4,2),gsl(0:4,2),gsl(0:4,2), &
                             dgsl(0:4,2,2),dgsl(0:4,2,2))
              CALL out_dprod(ddll(0:4,0:4,4,3),gsl(0:4,2),gsl(0:4,2), &
                             dgsl(0:4,2,3),dgsl(0:4,2,3))
              CALL set_dmat(ddll(0:4,0:4,2,1:3),r,2)
              dgmu(0:4,0:4,0,1:3) = ddll(0:4,0:4,4,1:3)
              dgmu(0:4,0:4,1,1:3) = -4.0_dp/3.0_dp*ddll(0:4,0:4,4,1:3) + &
                ddll(0:4,0:4,2,1:3)
              dgmu(0:4,0:4,2,1:3) = 1.0_dp/3.0_dp*ddll(0:4,0:4,4,1:3) -  &
                ddll(0:4,0:4,2,1:3)
            CASE (3)
!d-f
              CALL out_dprod(ddll(0:4,0:6,5,1),gsl(0:4,2),gsl(0:6,3), &
                             dgsl(0:4,2,1),dgsl(0:6,3,1))
              CALL out_dprod(ddll(0:4,0:6,5,2),gsl(0:4,2),gsl(0:6,3), &
                             dgsl(0:4,2,2),dgsl(0:6,3,2))
              CALL out_dprod(ddll(0:4,0:6,5,3),gsl(0:4,2),gsl(0:6,3), &
                             dgsl(0:4,2,3),dgsl(0:6,3,3))
              CALL set_dmat(ddll(0:4,0:6,3,1:3),r,3)
              CALL set_dmat(ddll(0:4,0:6,1,1:3),r,1)
              dgmu(0:4,0:6,0,1:3) = ddll(0:4,0:6,5,1:3)
              dgmu(0:4,0:6,1,1:3) = -SQRT(2.0_dp)*ddll(0:4,0:6,5,1:3) + &
                0.5_dp*ddll(0:4,0:6,3,1:3) + 0.5_dp*ddll(0:4,0:6,1,1:3)
              dgmu(0:4,0:6,2,1:3) = SQRT(0.2_dp)*ddll(0:4,0:6,5,1:3) - &
                SQRT(0.4_dp)*ddll(0:4,0:6,3,1:3)
            END SELECT
          CASE (3)
            SELECT CASE (l2)
            CASE (1)
!f-p
              CALL out_dprod(ddll(0:6,0:2,4,1),gsl(0:6,3),gsl(0:2,1), &
                             dgsl(0:6,3,1),dgsl(0:2,1,1))
              CALL out_dprod(ddll(0:6,0:2,4,2),gsl(0:6,3),gsl(0:2,1), &
                             dgsl(0:6,3,2),dgsl(0:2,1,2))
              CALL out_dprod(ddll(0:6,0:2,4,3),gsl(0:6,3),gsl(0:2,1), &
                             dgsl(0:6,3,3),dgsl(0:2,1,3))
              CALL set_dmat(ddll(0:6,0:2,2,1:3),r,2)
              dgmu(0:6,0:2,0,1:3) = ddll(0:6,0:2,4,1:3)
              dgmu(0:6,0:2,1,1:3) = -SQRT(1.5_dp)*ddll(0:6,0:2,4,1:3) + &
                0.5_dp*ddll(0:6,0:2,2,1:3)
            CASE (2)
!f-d
              CALL out_dprod(ddll(0:6,0:4,5,1),gsl(0:6,3),gsl(0:4,2), &
                            dgsl(0:6,3,1),dgsl(0:4,2,1))
              CALL out_dprod(ddll(0:6,0:4,5,2),gsl(0:6,3),gsl(0:4,2), &
                            dgsl(0:6,3,2),dgsl(0:4,2,2))
              CALL out_dprod(ddll(0:6,0:4,5,3),gsl(0:6,3),gsl(0:4,2), &
                            dgsl(0:6,3,3),dgsl(0:4,2,3))
              CALL set_dmat(ddll(0:6,0:4,3,1:3),r,3)
              CALL set_dmat(ddll(0:6,0:4,1,1:3),r,1)
              dgmu(0:6,0:4,0,1:3) = ddll(0:6,0:4,5,1:3)
              dgmu(0:6,0:4,1,1:3) = -SQRT(2.0_dp)*ddll(0:6,0:4,5,1:3) + &
                0.5_dp*ddll(0:6,0:4,3,1:3) + 0.5_dp*ddll(0:6,0:4,1,1:3)
              dgmu(0:6,0:4,2,1:3) = SQRT(0.2_dp)*ddll(0:6,0:4,5,1:3) - &
                SQRT(0.4_dp)*ddll(0:6,0:4,3,1:3)
            CASE (3)
!f-f
              CALL out_dprod(ddll(0:6,0:6,6,1),gsl(0:6,3),gsl(0:6,3), &
                             dgsl(0:6,3,1),dgsl(0:6,3,1))
              CALL out_dprod(ddll(0:6,0:6,6,2),gsl(0:6,3),gsl(0:6,3), &
                             dgsl(0:6,3,2),dgsl(0:6,3,2))
              CALL out_dprod(ddll(0:6,0:6,6,3),gsl(0:6,3),gsl(0:6,3), &
                             dgsl(0:6,3,3),dgsl(0:6,3,3))
              CALL set_dmat(ddll(0:6,0:6,4,1:3),r,4)
              CALL set_dmat(ddll(0:6,0:6,2,1:3),r,2)
              dgmu(0:6,0:6,0,1:3) = ddll(0:6,0:6,6,1:3)
              dgmu(0:6,0:6,1,1:3) = -1.5_dp*ddll(0:6,0:6,6,1:3) + &
                0.625_dp*ddll(0:6,0:6,4,1:3) + 0.625_dp*ddll(0:6,0:6,2,1:3)
              dgmu(0:6,0:6,2,1:3) = 0.6_dp*ddll(0:6,0:6,6,1:3) - ddll(0:6,0:6,4,1:3)
              dgmu(0:6,0:6,3,1:3) = -0.1_dp*ddll(0:6,0:6,6,1:3) + &
                0.375_dp*ddll(0:6,0:6,4,1:3) - 0.625_dp*ddll(0:6,0:6,2,1:3)
            END SELECT
          END SELECT
        END IF

      END SUBROUTINE dgmat
!------------------------------------------------------------------------------!
      SUBROUTINE set_mat(mat,r,num)
    REAL(KIND=dp), INTENT(OUT)               :: mat(0:,0:)
    REAL(KIND=dp), INTENT(IN)                :: r(3)
    INTEGER, INTENT(IN)                      :: num

    REAL(KIND=dp), PARAMETER                 :: s3 = 0.577350269189626_dp 

    INTEGER                                  :: l1, l2
    REAL(KIND=dp)                            :: x, xx, xy, xz, y, yy, yz, z, &
                                                zz

! sqrt(1.0_dp/3.0_dp)

        l1 = SIZE(mat,1)
        l2 = SIZE(mat,2)
        x = r(1)
        y = r(2)
        z = r(3)
        IF (l1==5 .AND. l2==3 .AND. num==1) THEN
! d1(d,p)
          mat(0:4,0:2) = zero
          mat(0,0) = 2.0_dp*s3*z
          mat(1,0) = x
          mat(2,0) = y
          mat(0,1) = -s3*x
          mat(1,1) = z
          mat(3,1) = x
          mat(4,1) = y
          mat(0,2) = -s3*y
          mat(2,2) = z
          mat(3,2) = -y
          mat(4,2) = x
        ELSE IF (l1==3 .AND. l2==5 .AND. num==1) THEN
! d1(p,d)
          mat(0:2,0:4) = zero
          mat(0,0) = 2.0_dp*s3*z
          mat(0,1) = x
          mat(0,2) = y
          mat(1,0) = -s3*x
          mat(1,1) = z
          mat(1,3) = x
          mat(1,4) = y
          mat(2,0) = -s3*y
          mat(2,2) = z
          mat(2,3) = -y
          mat(2,4) = x
        ELSE
          xx = x*x
          yy = y*y
          zz = z*z
          xy = x*y
          xz = x*z
          yz = y*z
          IF (l1==5 .AND. l2==5 .AND. num==2) THEN
! d2(d,d)
            mat(0,0) = zz - 2.0_dp/3.0_dp
            mat(0,1) = s3*xz
            mat(0,2) = s3*yz
            mat(0,3) = -s3*(xx-yy)
            mat(0,4) = -2.0_dp*s3*xy
            mat(1,0) = mat(0,1)
            mat(2,0) = mat(0,2)
            mat(3,0) = mat(0,3)
            mat(4,0) = mat(0,4)
            mat(1,1) = -yy
            mat(1,2) = xy
            mat(1,3) = xz
            mat(1,4) = yz
            mat(2,1) = mat(1,2)
            mat(3,1) = mat(1,3)
            mat(4,1) = mat(1,4)
            mat(2,2) = -xx
            mat(2,3) = -yz
            mat(2,4) = xz
            mat(3,2) = mat(2,3)
            mat(4,2) = mat(2,4)
            mat(3,3) = -zz
            mat(3,4) = zero
            mat(4,3) = zero
            mat(4,4) = -zz
          ELSE IF (l1==7 .AND. l2==3 .AND. num==2) THEN
! d2(f,p)
            mat(0,0) = SQRT(1.5_dp)*(3.0_dp*zz-1.0_dp)
            mat(0,1) = -SQRT(6.0_dp)*xz
            mat(0,2) = -SQRT(6.0_dp)*yz
            mat(1,0) = 4.0_dp*xz
            mat(1,1) = 0.5_dp*(5.0_dp*zz-2.0_dp*xx-1.0_dp)
            mat(1,2) = -xy
            mat(2,0) = 4.0_dp*yz
            mat(2,1) = -xy
            mat(2,2) = 0.5_dp*(5.0_dp*zz-2.0_dp*yy-1.0_dp)
            mat(3,0) = SQRT(2.5_dp)*(xx-yy)
            mat(3,1) = SQRT(10.0_dp)*xz
            mat(3,2) = -SQRT(10.0_dp)*yz
            mat(4,0) = SQRT(10.0_dp)*xy
            mat(4,1) = SQRT(10.0_dp)*yz
            mat(4,2) = SQRT(10.0_dp)*xz
            mat(5,0) = zero
            mat(5,1) = SQRT(3.75_dp)*(xx-yy)
            mat(5,2) = -SQRT(15.0_dp)*xy
            mat(6,0) = zero
            mat(6,1) = SQRT(15.0_dp)*xy
            mat(6,2) = SQRT(3.75_dp)*(xx-yy)
          ELSE IF (l1==3 .AND. l2==7 .AND. num==2) THEN
! d2(p,f)
            mat(0,0) = SQRT(1.5_dp)*(3.0_dp*zz-1.0_dp)
            mat(1,0) = -SQRT(6.0_dp)*xz
            mat(2,0) = -SQRT(6.0_dp)*yz
            mat(0,1) = 4.0_dp*xz
            mat(1,1) = 0.5_dp*(5.0_dp*zz-2.0_dp*xx-1.0_dp)
            mat(2,1) = -xy
            mat(0,2) = 4.0_dp*yz
            mat(1,2) = -xy
            mat(2,2) = 0.5_dp*(5.0_dp*zz-2.0_dp*yy-1.0_dp)
            mat(0,3) = SQRT(2.5_dp)*(xx-yy)
            mat(1,3) = SQRT(10.0_dp)*xz
            mat(2,3) = -SQRT(10.0_dp)*yz
            mat(0,4) = SQRT(10.0_dp)*xy
            mat(1,4) = SQRT(10.0_dp)*yz
            mat(2,4) = SQRT(10.0_dp)*xz
            mat(0,5) = zero
            mat(1,5) = SQRT(3.75_dp)*(xx-yy)
            mat(2,5) = -SQRT(15.0_dp)*xy
            mat(0,6) = zero
            mat(1,6) = SQRT(15.0_dp)*xy
            mat(2,6) = SQRT(3.75_dp)*(xx-yy)
          ELSE IF (l1==7 .AND. l2==5 .AND. num==1) THEN
! d1(f,d)
            mat(0,0) = SQRT(4.5_dp)*z
            mat(0,1) = -SQRT(1.5_dp)*x
            mat(0,2) = -SQRT(1.5_dp)*y
            mat(0,3) = zero
            mat(0,4) = zero
            mat(1,0) = SQRT(3.0_dp)*x
            mat(1,1) = 2.0_dp*z
            mat(1,2) = zero
            mat(1,3) = -0.5_dp*x
            mat(1,4) = -0.5_dp*y
            mat(2,0) = SQRT(3.0_dp)*y
            mat(2,1) = zero
            mat(2,2) = 2.0_dp*z
            mat(2,3) = 0.5_dp*y
            mat(2,4) = -0.5_dp*x
            mat(3,0) = zero
            mat(3,1) = SQRT(2.5_dp)*x
            mat(3,2) = -SQRT(2.5_dp)*y
            mat(3,3) = SQRT(2.5_dp)*z
            mat(3,4) = zero
            mat(4,0) = zero
            mat(4,1) = SQRT(2.5_dp)*y
            mat(4,2) = SQRT(2.5_dp)*x
            mat(4,3) = zero
            mat(4,4) = SQRT(2.5_dp)*z
            mat(5,0) = zero
            mat(5,1) = zero
            mat(5,2) = zero
            mat(5,3) = SQRT(3.75_dp)*x
            mat(5,4) = -SQRT(3.75_dp)*y
            mat(6,0) = zero
            mat(6,1) = zero
            mat(6,2) = zero
            mat(6,3) = SQRT(3.75_dp)*y
            mat(6,4) = SQRT(3.75_dp)*x
          ELSE IF (l1==5 .AND. l2==7 .AND. num==1) THEN
! d1(d,f)
            mat(0,0) = SQRT(4.5_dp)*z
            mat(1,0) = -SQRT(1.5_dp)*x
            mat(2,0) = -SQRT(1.5_dp)*y
            mat(3,0) = zero
            mat(4,0) = zero
            mat(0,1) = SQRT(3.0_dp)*x
            mat(1,1) = 2.0_dp*z
            mat(2,1) = zero
            mat(3,1) = -0.5_dp*x
            mat(4,1) = -0.5_dp*y
            mat(0,2) = SQRT(3.0_dp)*y
            mat(1,2) = zero
            mat(2,2) = 2.0_dp*z
            mat(3,2) = 0.5_dp*y
            mat(4,2) = -0.5_dp*x
            mat(0,3) = zero
            mat(1,3) = SQRT(2.5_dp)*x
            mat(2,3) = -SQRT(2.5_dp)*y
            mat(3,3) = SQRT(2.5_dp)*z
            mat(4,3) = zero
            mat(0,4) = zero
            mat(1,4) = SQRT(2.5_dp)*y
            mat(2,4) = SQRT(2.5_dp)*x
            mat(3,4) = zero
            mat(4,4) = SQRT(2.5_dp)*z
            mat(0,5) = zero
            mat(1,5) = zero
            mat(2,5) = zero
            mat(3,5) = SQRT(3.75_dp)*x
            mat(4,5) = -SQRT(3.75_dp)*y
            mat(0,6) = zero
            mat(1,6) = zero
            mat(2,6) = zero
            mat(3,6) = SQRT(3.75_dp)*y
            mat(4,6) = SQRT(3.75_dp)*x
          ELSE IF (l1==7 .AND. l2==5 .AND. num==3) THEN
! d3(f,d)
            mat(0,0) = SQRT(0.5_dp)*(4.0_dp*zz-3.0_dp)*z
            mat(0,1) = SQRT(1.5_dp)*x*zz
            mat(0,2) = SQRT(1.5_dp)*y*zz
            mat(0,3) = -SQRT(6.0_dp)*(xx-yy)*z
            mat(0,4) = -2.0_dp*SQRT(6.0_dp)*x*yz
            mat(1,0) = SQRT(0.75_dp)*(3.0_dp*zz-1.0_dp)*x
            mat(1,1) = 0.5_dp*(xx-5.0_dp*yy)*z
            mat(1,2) = 3.0_dp*x*yz
            mat(1,3) = (2.5_dp*zz-xx+yy)*x
            mat(1,4) = 0.5_dp*(5.0_dp*zz-4.0_dp*xx)*y
            mat(2,0) = SQRT(0.75_dp)*(3.0_dp*zz-1.0_dp)*y
            mat(2,1) = 3.0_dp*x*yz
            mat(2,2) = 0.5_dp*(yy-5.0_dp*xx)*z
            mat(2,3) = -(2.5_dp*zz+xx-yy)*y
            mat(2,4) = 0.5_dp*(5.0_dp*zz-4.0_dp*yy)*x
            mat(3,0) = zero
            mat(3,1) = SQRT(2.5_dp)*(zz-2.0_dp*yy)*x
            mat(3,2) = -SQRT(2.5_dp)*(zz-2.0_dp*xx)*y
            mat(3,3) = SQRT(2.5_dp)*(1.0_dp-2.0_dp*zz)*z
            mat(3,4) = zero
            mat(4,0) = zero
            mat(4,1) = SQRT(2.5_dp)*(1.0_dp-2.0_dp*yy)*y
            mat(4,2) = SQRT(2.5_dp)*(1.0_dp-2.0_dp*xx)*x
            mat(4,3) = zero
            mat(4,4) = SQRT(2.5_dp)*(1.0_dp-2.0_dp*zz)*z
            mat(5,0) = SQRT(1.25_dp)*(3.0_dp*yy-xx)*x
            mat(5,1) = SQRT(3.75_dp)*(xx-yy)*z
            mat(5,2) = -SQRT(15.0_dp)*x*yz
            mat(5,3) = -SQRT(3.75_dp)*x*zz
            mat(5,4) = SQRT(3.75_dp)*y*zz
            mat(6,0) = -SQRT(1.25_dp)*(3.0_dp*xx-yy)*y
            mat(6,1) = SQRT(15.0_dp)*x*yz
            mat(6,2) = SQRT(3.75_dp)*(xx-yy)*z
            mat(6,3) = -SQRT(3.75_dp)*y*zz
            mat(6,4) = -SQRT(3.75_dp)*x*zz
          ELSE IF (l1==5 .AND. l2==7 .AND. num==3) THEN
! d3(d,f)
            mat(0,0) = SQRT(0.5_dp)*(4.0_dp*zz-3.0_dp)*z
            mat(1,0) = SQRT(1.5_dp)*x*zz
            mat(2,0) = SQRT(1.5_dp)*y*zz
            mat(3,0) = -SQRT(6.0_dp)*(xx-yy)*z
            mat(4,0) = -2.0_dp*SQRT(6.0_dp)*x*yz
            mat(0,1) = SQRT(0.75_dp)*(3.0_dp*zz-1.0_dp)*x
            mat(1,1) = 0.5_dp*(xx-5.0_dp*yy)*z
            mat(2,1) = 3.0_dp*x*yz
            mat(3,1) = (2.5_dp*zz-xx+yy)*x
            mat(4,1) = 0.5_dp*(5.0_dp*zz-4.0_dp*xx)*y
            mat(0,2) = SQRT(0.75_dp)*(3.0_dp*zz-1.0_dp)*y
            mat(1,2) = 3.0_dp*x*yz
            mat(2,2) = 0.5_dp*(yy-5.0_dp*xx)*z
            mat(3,2) = -(2.5_dp*zz+xx-yy)*y
            mat(4,2) = 0.5_dp*(5.0_dp*zz-4.0_dp*yy)*x
            mat(0,3) = zero
            mat(1,3) = SQRT(2.5_dp)*(zz-2.0_dp*yy)*x
            mat(2,3) = -SQRT(2.5_dp)*(zz-2.0_dp*xx)*y
            mat(3,3) = SQRT(2.5_dp)*(1.0_dp-2.0_dp*zz)*z
            mat(4,3) = zero
            mat(0,4) = zero
            mat(1,4) = SQRT(2.5_dp)*(1.0_dp-2.0_dp*yy)*y
            mat(2,4) = SQRT(2.5_dp)*(1.0_dp-2.0_dp*xx)*x
            mat(3,4) = zero
            mat(4,4) = SQRT(2.5_dp)*(1.0_dp-2.0_dp*zz)*z
            mat(0,5) = SQRT(1.25_dp)*(3.0_dp*yy-xx)*x
            mat(1,5) = SQRT(3.75_dp)*(xx-yy)*z
            mat(2,5) = -SQRT(15.0_dp)*x*yz
            mat(3,5) = -SQRT(3.75_dp)*x*zz
            mat(4,5) = SQRT(3.75_dp)*y*zz
            mat(0,6) = -SQRT(1.25_dp)*(3.0_dp*xx-yy)*y
            mat(1,6) = SQRT(15.0_dp)*x*yz
            mat(2,6) = SQRT(3.75_dp)*(xx-yy)*z
            mat(3,6) = -SQRT(3.75_dp)*y*zz
            mat(4,6) = -SQRT(3.75_dp)*x*zz
          ELSE IF (l1==7 .AND. l2==7 .AND. num==2) THEN
! d2(f,f)
            mat(0,0) = 0.4_dp*(3.0_dp*zz-1.0_dp)
            mat(0,1) = SQRT(0.24_dp)*xz
            mat(0,2) = SQRT(0.24_dp)*yz
            mat(0,3) = -SQRT(0.6_dp)*(xx-yy)
            mat(0,4) = -SQRT(2.4_dp)*xy
            mat(0,5) = zero
            mat(0,6) = zero
            mat(1,0) = mat(0,1)
            mat(2,0) = mat(0,2)
            mat(3,0) = mat(0,3)
            mat(4,0) = mat(0,4)
            mat(5,0) = mat(0,5)
            mat(6,0) = mat(0,6)
            mat(1,1) = 0.3_dp*(1.0_dp-4.0_dp*yy+zz)
            mat(1,2) = 1.2_dp*xy
            mat(1,3) = SQRT(0.9_dp)*xz
            mat(1,4) = SQRT(0.9_dp)*yz
            mat(1,5) = -SQRT(0.15_dp)*(xx-yy)
            mat(1,6) = -SQRT(0.6_dp)*xy
            mat(2,1) = mat(1,2)
            mat(3,1) = mat(1,3)
            mat(4,1) = mat(1,4)
            mat(5,1) = mat(1,5)
            mat(6,1) = mat(1,6)
            mat(2,2) = 0.3_dp*(1.0_dp-4.0_dp*xx+zz)
            mat(2,3) = -SQRT(0.9_dp)*yz
            mat(2,4) = SQRT(0.9_dp)*xz
            mat(2,5) = SQRT(0.6_dp)*xy
            mat(2,6) = -SQRT(0.15_dp)*(xx-yy)
            mat(3,2) = mat(2,3)
            mat(4,2) = mat(2,4)
            mat(5,2) = mat(2,5)
            mat(6,2) = mat(2,6)
            mat(3,3) = zero
            mat(3,4) = zero
            mat(3,5) = SQRT(1.5_dp)*xz
            mat(3,6) = SQRT(1.5_dp)*yz
            mat(4,3) = mat(3,4)
            mat(5,3) = mat(3,5)
            mat(6,3) = mat(3,6)
            mat(4,4) = zero
            mat(4,5) = -SQRT(1.5_dp)*yz
            mat(4,6) = SQRT(1.5_dp)*xz
            mat(5,4) = mat(4,5)
            mat(6,4) = mat(4,6)
            mat(5,5) = -0.5_dp*(3.0_dp*zz-1.0_dp)
            mat(5,6) = zero
            mat(6,5) = zero
            mat(6,6) = -0.5_dp*(3.0_dp*zz-1.0_dp)
          ELSE IF (l1==7 .AND. l2==7 .AND. num==4) THEN
! d4(f,f)
            mat(0,0) = 0.6_dp*(5.0_dp*zz-4.0_dp)*zz
            mat(0,1) = SQRT(0.24_dp)*(5.0_dp*zz-2.0_dp)*xz
            mat(0,2) = SQRT(0.24_dp)*(5.0_dp*zz-2.0_dp)*yz
            mat(0,3) = -SQRT(0.6_dp)*(xx-yy)*zz
            mat(0,4) = -SQRT(2.4_dp)*xy*zz
            mat(0,5) = SQRT(3.6_dp)*(3.0_dp*yy-xx)*xz
            mat(0,6) = SQRT(3.6_dp)*(yy-3.0_dp*xx)*yz
            mat(1,0) = mat(0,1)
            mat(2,0) = mat(0,2)
            mat(3,0) = mat(0,3)
            mat(4,0) = mat(0,4)
            mat(5,0) = mat(0,5)
            mat(6,0) = mat(0,6)
            mat(1,1) = -0.4_dp*xx - 0.5_dp*(5.0_dp*yy-3.0_dp*xx)*zz
            mat(1,2) = 0.4_dp*(10.0_dp*zz-1.0_dp)*xy
            mat(1,3) = SQRT(0.1_dp)*(6.0_dp*zz-8.0_dp*yy-1.0_dp)*xz
            mat(1,4) = SQRT(0.1_dp)*(2.0_dp*zz-8.0_dp*yy+3.0_dp)*yz
            mat(1,5) = SQRT(3.75_dp)*(xx-yy-1.4_dp*xx*xx+1.2_dp*xx*yy+yy*yy &
              )
            mat(1,6) = SQRT(0.6_dp)*(4.0_dp*zz-4.0_dp*xx+1.0_dp)*xy
            mat(2,1) = mat(1,2)
            mat(3,1) = mat(1,3)
            mat(4,1) = mat(1,4)
            mat(5,1) = mat(1,5)
            mat(6,1) = mat(1,6)
            mat(2,2) = -0.4_dp*yy - 0.5_dp*(5.0_dp*xx-3.0_dp*yy)*zz
            mat(2,3) = -SQRT(0.1_dp)*(6.0_dp*zz-8.0_dp*xx-1.0_dp)*yz
            mat(2,4) = SQRT(0.1_dp)*(2.0_dp*zz-8.0_dp*xx+3.0_dp)*xz
            mat(2,5) = SQRT(0.6_dp)*(4.0_dp*yy-4.0_dp*zz-1.0_dp)*xy
            mat(2,6) = SQRT(3.75_dp)*(xx-yy-xx*xx-1.2_dp*xx*yy+1.4_dp*yy*yy &
              )
            mat(3,2) = mat(2,3)
            mat(4,2) = mat(2,4)
            mat(5,2) = mat(2,5)
            mat(6,2) = mat(2,6)
            mat(3,3) = (2.0_dp-3.0_dp*zz)*zz - 4.0_dp*xx*yy
            mat(3,4) = 2.0_dp*(xx-yy)*xy
            mat(3,5) = -SQRT(1.5_dp)*(2.0_dp*zz-1.0_dp)*xz
            mat(3,6) = -SQRT(1.5_dp)*(2.0_dp*zz-1.0_dp)*yz
            mat(4,3) = mat(3,4)
            mat(5,3) = mat(3,5)
            mat(6,3) = mat(3,6)
            mat(4,4) = -(2.0_dp*zz-1.0_dp)**2 + 4.0_dp*xx*yy
            mat(4,5) = SQRT(1.5_dp)*(2.0_dp*zz-1.0_dp)*yz
            mat(4,6) = -SQRT(1.5_dp)*(2.0_dp*zz-1.0_dp)*xz
            mat(5,4) = mat(4,5)
            mat(6,4) = mat(4,6)
            mat(5,5) = 1.5_dp*(zz-1.0_dp)*zz
            mat(5,6) = zero
            mat(6,5) = zero
            mat(6,6) = 1.5_dp*(zz-1.0_dp)*zz
          END IF
        END IF
      END SUBROUTINE set_mat
!------------------------------------------------------------------------------!
      SUBROUTINE set_dmat(dmat,r,num)
    REAL(KIND=dp), INTENT(OUT)               :: dmat(0:,0:,1:)
    REAL(KIND=dp), INTENT(IN)                :: r(3)
    INTEGER, INTENT(IN)                      :: num

    REAL(KIND=dp), PARAMETER                 :: s3 = 0.577350269189626_dp 

    INTEGER                                  :: l1, l2
    REAL(KIND=dp)                            :: dc(0:6,0:6), x, xx, xy, xz, &
                                                y, yy, yz, z, zz

! sqrt(1.0_dp/3.0_dp)

        l1 = SIZE(dmat,1)
        l2 = SIZE(dmat,2)
        x = r(1)
        y = r(2)
        z = r(3)
        IF (l1==5 .AND. l2==3 .AND. num==1) THEN
! d1(d,p)
          dmat(0:4,0:2,1:3) = zero
          dmat(1,0,1) = one
          dmat(0,1,1) = -s3
          dmat(3,1,1) = one
          dmat(4,2,1) = one
          dmat(2,0,2) = one
          dmat(0,2,2) = -s3
          dmat(3,2,2) = -one
          dmat(4,1,2) = one
          dmat(0,0,3) = 2.0_dp*s3
          dmat(1,1,3) = one
          dmat(2,2,3) = one
        ELSE IF (l1==3 .AND. l2==5 .AND. num==1) THEN
! d1(p,d)
          dmat(0:2,0:4,1:3) = zero
          dmat(0,1,1) = one
          dmat(1,0,1) = -s3
          dmat(1,3,1) = one
          dmat(2,4,1) = one
          dmat(0,2,2) = one
          dmat(2,0,2) = -s3
          dmat(2,3,2) = -one
          dmat(1,4,2) = one
          dmat(0,0,3) = 2.0_dp*s3
          dmat(1,1,3) = one
          dmat(2,2,3) = one
        ELSEIF (l1==5 .AND. l2==5 .AND. num==2) THEN
! d2(d,d)
          dmat(0:4,0:4,1:3) = zero
          dmat(0,1,1) = s3*z
          dmat(0,3,1) = -2.0_dp*s3*x
          dmat(0,4,1) = -2.0_dp*s3*y
          dmat(1,0,1) = dmat(0,1,1)
          dmat(3,0,1) = dmat(0,3,1)
          dmat(4,0,1) = dmat(0,4,1)
          dmat(1,2,1) = y
          dmat(1,3,1) = z
          dmat(2,1,1) = y
          dmat(3,1,1) = z
          dmat(2,2,1) = -2.0_dp*x
          dmat(2,4,1) = z
          dmat(4,2,1) = z
          dmat(0,2,2) = s3*z
          dmat(0,3,2) = 2.0_dp*s3*y
          dmat(0,4,2) = -2.0_dp*s3*x
          dmat(2,0,2) = dmat(0,2,2)
          dmat(3,0,2) = dmat(0,3,2)
          dmat(4,0,2) = dmat(0,4,2)
          dmat(1,1,2) = -2.0_dp*y
          dmat(1,2,2) = x
          dmat(1,4,2) = z
          dmat(2,1,2) = x
          dmat(4,1,2) = z
          dmat(2,3,2) = -z
          dmat(3,2,2) = -z
          dmat(0,0,3) = 2.0_dp*z
          dmat(0,1,3) = s3*x
          dmat(0,2,3) = s3*y
          dmat(1,0,3) = dmat(0,1,3)
          dmat(2,0,3) = dmat(0,2,3)
          dmat(1,3,3) = x
          dmat(1,4,3) = y
          dmat(3,1,3) = x
          dmat(4,1,3) = y
          dmat(2,3,3) = -y
          dmat(2,4,3) = x
          dmat(3,2,3) = -y
          dmat(4,2,3) = x
          dmat(3,3,3) = -2.0_dp*z
          dmat(4,4,3) = -2.0_dp*z
        ELSE IF (l1==7 .AND. l2==3 .AND. num==2) THEN
! d2(f,p)
          dmat(0,0,1) = zero
          dmat(0,1,1) = -SQRT(6.0_dp)*z
          dmat(0,2,1) = zero
          dmat(1,0,1) = 4.0_dp*z
          dmat(1,1,1) = -2.0_dp*x
          dmat(1,2,1) = -y
          dmat(2,0,1) = zero
          dmat(2,1,1) = -y
          dmat(2,2,1) = zero
          dmat(3,0,1) = 2.0_dp*SQRT(2.5_dp)*x
          dmat(3,1,1) = SQRT(10.0_dp)*z
          dmat(3,2,1) = zero
          dmat(4,0,1) = SQRT(10.0_dp)*y
          dmat(4,1,1) = zero
          dmat(4,2,1) = SQRT(10.0_dp)*z
          dmat(5,0,1) = zero
          dmat(5,1,1) = SQRT(15.0_dp)*x
          dmat(5,2,1) = -SQRT(15.0_dp)*y
          dmat(6,0,1) = zero
          dmat(6,1,1) = SQRT(15.0_dp)*y
          dmat(6,2,1) = SQRT(15.0_dp)*x

          dmat(0,0,2) = zero
          dmat(0,1,2) = zero
          dmat(0,2,2) = -SQRT(6.0_dp)*z
          dmat(1,0,2) = zero
          dmat(1,1,2) = zero
          dmat(1,2,2) = -x
          dmat(2,0,2) = 4.0_dp*z
          dmat(2,1,2) = -x
          dmat(2,2,2) = -2.0_dp*y
          dmat(3,0,2) = -2.0_dp*SQRT(2.5_dp)*y
          dmat(3,1,2) = zero
          dmat(3,2,2) = -SQRT(10.0_dp)*z
          dmat(4,0,2) = SQRT(10.0_dp)*x
          dmat(4,1,2) = SQRT(10.0_dp)*z
          dmat(4,2,2) = zero
          dmat(5,0,2) = zero
          dmat(5,1,2) = -SQRT(15.0_dp)*y
          dmat(5,2,2) = -SQRT(15.0_dp)*x
          dmat(6,0,2) = zero
          dmat(6,1,2) = SQRT(15.0_dp)*x
          dmat(6,2,2) = -SQRT(15.0_dp)*y

          dmat(0,0,3) = 6.0_dp*SQRT(1.5_dp)*z
          dmat(0,1,3) = -SQRT(6.0_dp)*x
          dmat(0,2,3) = -SQRT(6.0_dp)*y
          dmat(1,0,3) = 4.0_dp*x
          dmat(1,1,3) = 5.0_dp*z
          dmat(1,2,3) = zero
          dmat(2,0,3) = 4.0_dp*y
          dmat(2,1,3) = zero
          dmat(2,2,3) = 5.0_dp*z
          dmat(3,0,3) = zero
          dmat(3,1,3) = SQRT(10.0_dp)*x
          dmat(3,2,3) = -SQRT(10.0_dp)*y
          dmat(4,0,3) = zero
          dmat(4,1,3) = SQRT(10.0_dp)*y
          dmat(4,2,3) = SQRT(10.0_dp)*x
          dmat(5,0,3) = zero
          dmat(5,1,3) = zero
          dmat(5,2,3) = zero
          dmat(6,0,3) = zero
          dmat(6,1,3) = zero
          dmat(6,2,3) = zero
        ELSE IF (l1==3 .AND. l2==7 .AND. num==2) THEN
! d2(p,f)
          dmat(0,0,1) = zero
          dmat(1,0,1) = -SQRT(6.0_dp)*z
          dmat(2,0,1) = zero
          dmat(0,1,1) = 4.0_dp*z
          dmat(1,1,1) = -2.0_dp*x
          dmat(2,1,1) = -y
          dmat(0,2,1) = zero
          dmat(1,2,1) = -y
          dmat(2,2,1) = zero
          dmat(0,3,1) = 2.0_dp*SQRT(2.5_dp)*x
          dmat(1,3,1) = SQRT(10.0_dp)*z
          dmat(2,3,1) = zero
          dmat(0,4,1) = SQRT(10.0_dp)*y
          dmat(1,4,1) = zero
          dmat(2,4,1) = SQRT(10.0_dp)*z
          dmat(0,5,1) = zero
          dmat(1,5,1) = SQRT(15.0_dp)*x
          dmat(2,5,1) = -SQRT(15.0_dp)*y
          dmat(0,6,1) = zero
          dmat(1,6,1) = SQRT(15.0_dp)*y
          dmat(2,6,1) = SQRT(15.0_dp)*x

          dmat(0,0,2) = zero
          dmat(1,0,2) = zero
          dmat(2,0,2) = -SQRT(6.0_dp)*z
          dmat(0,1,2) = zero
          dmat(1,1,2) = zero
          dmat(2,1,2) = -x
          dmat(0,2,2) = 4.0_dp*z
          dmat(1,2,2) = -x
          dmat(2,2,2) = -2.0_dp*y
          dmat(0,3,2) = -2.0_dp*SQRT(2.5_dp)*y
          dmat(1,3,2) = zero
          dmat(2,3,2) = -SQRT(10.0_dp)*z
          dmat(0,4,2) = SQRT(10.0_dp)*x
          dmat(1,4,2) = SQRT(10.0_dp)*z
          dmat(2,4,2) = zero
          dmat(0,5,2) = zero
          dmat(1,5,2) = -SQRT(15.0_dp)*y
          dmat(2,5,2) = -SQRT(15.0_dp)*x
          dmat(0,6,2) = zero
          dmat(1,6,2) = SQRT(15.0_dp)*x
          dmat(2,6,2) = -SQRT(15.0_dp)*y

          dmat(0,0,3) = 6.0_dp*SQRT(1.5_dp)*z
          dmat(1,0,3) = -SQRT(6.0_dp)*x
          dmat(2,0,3) = -SQRT(6.0_dp)*y
          dmat(0,1,3) = 4.0_dp*x
          dmat(1,1,3) = 5.0_dp*z
          dmat(2,1,3) = zero
          dmat(0,2,3) = 4.0_dp*y
          dmat(1,2,3) = zero
          dmat(2,2,3) = 5.0_dp*z
          dmat(0,3,3) = zero
          dmat(1,3,3) = SQRT(10.0_dp)*x
          dmat(2,3,3) = -SQRT(10.0_dp)*y
          dmat(0,4,3) = zero
          dmat(1,4,3) = SQRT(10.0_dp)*y
          dmat(2,4,3) = SQRT(10.0_dp)*x
          dmat(0,5,3) = zero
          dmat(1,5,3) = zero
          dmat(2,5,3) = zero
          dmat(0,6,3) = zero
          dmat(1,6,3) = zero
          dmat(2,6,3) = zero
        ELSE IF (l1==7 .AND. l2==5 .AND. num==1) THEN
! d1(f,d)
          dmat(0:6,0:4,1:3) = zero
          dmat(0,1,1) = -SQRT(1.5_dp)
          dmat(1,0,1) = SQRT(3.0_dp)
          dmat(1,3,1) = -0.5_dp
          dmat(2,4,1) = -0.5_dp
          dmat(3,1,1) = SQRT(2.5_dp)
          dmat(4,2,1) = SQRT(2.5_dp)
          dmat(5,3,1) = SQRT(3.75_dp)
          dmat(6,4,1) = SQRT(3.75_dp)

          dmat(0,2,2) = -SQRT(1.5_dp)
          dmat(1,4,2) = -0.5_dp
          dmat(2,0,2) = SQRT(3.0_dp)
          dmat(2,3,2) = 0.5_dp
          dmat(3,2,2) = -SQRT(2.5_dp)
          dmat(4,1,2) = SQRT(2.5_dp)
          dmat(5,4,2) = -SQRT(3.75_dp)
          dmat(6,3,2) = SQRT(3.75_dp)

          dmat(0,0,3) = SQRT(4.5_dp)
          dmat(1,1,3) = 2.0_dp
          dmat(2,2,3) = 2.0_dp
          dmat(3,3,3) = SQRT(2.5_dp)
          dmat(4,4,3) = SQRT(2.5_dp)
        ELSE IF (l1==5 .AND. l2==7 .AND. num==1) THEN
! d1(d,f)
          dmat(0:4,0:6,1:3) = zero
          dmat(1,0,1) = -SQRT(1.5_dp)
          dmat(0,1,1) = SQRT(3.0_dp)
          dmat(3,1,1) = -0.5_dp
          dmat(4,2,1) = -0.5_dp
          dmat(1,3,1) = SQRT(2.5_dp)
          dmat(2,4,1) = SQRT(2.5_dp)
          dmat(3,5,1) = SQRT(3.75_dp)
          dmat(4,6,1) = SQRT(3.75_dp)

          dmat(2,0,2) = -SQRT(1.5_dp)
          dmat(4,1,2) = -0.5_dp
          dmat(0,2,2) = SQRT(3.0_dp)
          dmat(3,2,2) = 0.5_dp
          dmat(2,3,2) = -SQRT(2.5_dp)
          dmat(1,4,2) = SQRT(2.5_dp)
          dmat(4,5,2) = -SQRT(3.75_dp)
          dmat(3,6,2) = SQRT(3.75_dp)

          dmat(0,0,3) = SQRT(4.5_dp)
          dmat(1,1,3) = 2.0_dp
          dmat(2,2,3) = 2.0_dp
          dmat(3,3,3) = SQRT(2.5_dp)
          dmat(4,4,3) = SQRT(2.5_dp)
        ELSE IF (l1==7 .AND. l2==5 .AND. num==3) THEN
          xx = x*x
          yy = y*y
          zz = z*z
          xy = x*y
          xz = x*z
          yz = y*z
! d3(f,d)
          dmat(0,0,1) = zero
          dmat(0,1,1) = SQRT(1.5_dp)*zz
          dmat(0,2,1) = zero
          dmat(0,3,1) = -2.0_dp*SQRT(6.0_dp)*xz
          dmat(0,4,1) = -2.0_dp*SQRT(6.0_dp)*yz
          dmat(1,0,1) = SQRT(0.75_dp)*(3.0_dp*zz-1.0_dp)
          dmat(1,1,1) = xz
          dmat(1,2,1) = 3.0_dp*yz
          dmat(1,3,1) = 2.5_dp*zz-3.0_dp*xx+yy
          dmat(1,4,1) = -4.0_dp*xy
          dmat(2,0,1) = zero
          dmat(2,1,1) = 3.0_dp*yz
          dmat(2,2,1) = -5.0_dp*xz
          dmat(2,3,1) = -2.0_dp*xy
          dmat(2,4,1) = 0.5_dp*(5.0_dp*zz-4.0_dp*yy)
          dmat(3,0,1) = zero
          dmat(3,1,1) = SQRT(2.5_dp)*(zz-2.0_dp*yy)
          dmat(3,2,1) = SQRT(2.5_dp)*4.0_dp*xy
          dmat(3,3,1) = zero
          dmat(3,4,1) = zero
          dmat(4,0,1) = zero
          dmat(4,1,1) = zero
          dmat(4,2,1) = SQRT(2.5_dp)*(1.0_dp-6.0_dp*xx)
          dmat(4,3,1) = zero
          dmat(4,4,1) = zero
          dmat(5,0,1) = 3.0_dp*SQRT(1.25_dp)*(yy-xx)
          dmat(5,1,1) = SQRT(15.0_dp)*xz
          dmat(5,2,1) = -SQRT(15.0_dp)*yz
          dmat(5,3,1) = -SQRT(3.75_dp)*zz
          dmat(5,4,1) = zero
          dmat(6,0,1) = -SQRT(5.0_dp)*3.0_dp*xy
          dmat(6,1,1) = SQRT(15.0_dp)*yz
          dmat(6,2,1) = SQRT(15.0_dp)*xz
          dmat(6,3,1) = zero
          dmat(6,4,1) = -SQRT(3.75_dp)*zz

          dmat(0,0,2) = zero
          dmat(0,1,2) = zero
          dmat(0,2,2) = SQRT(1.5_dp)*zz
          dmat(0,3,2) = 2.0_dp*SQRT(6.0_dp)*yz
          dmat(0,4,2) = -2.0_dp*SQRT(6.0_dp)*xz
          dmat(1,0,2) = zero
          dmat(1,1,2) = -5.0_dp*yz
          dmat(1,2,2) = 3.0_dp*xz
          dmat(1,3,2) = 2.0_dp*xy
          dmat(1,4,2) = 0.5_dp*(5.0_dp*zz-4.0_dp*xx)
          dmat(2,0,2) = SQRT(0.75_dp)*(3.0_dp*zz-1.0_dp)
          dmat(2,1,2) = 3.0_dp*xz
          dmat(2,2,2) = yz
          dmat(2,3,2) = -(2.5_dp*zz+xx-3.0_dp*yy)
          dmat(2,4,2) = -4.0_dp*xy
          dmat(3,0,2) = zero
          dmat(3,1,2) = -SQRT(2.5_dp)*4.0_dp*xy
          dmat(3,2,2) = -SQRT(2.5_dp)*(zz-2.0_dp*xx)
          dmat(3,3,2) = zero
          dmat(3,4,2) = zero
          dmat(4,0,2) = zero
          dmat(4,1,2) = SQRT(2.5_dp)*(1.0_dp-6.0_dp*yy)
          dmat(4,2,2) = zero
          dmat(4,3,2) = zero
          dmat(4,4,2) = zero
          dmat(5,0,2) = 3.0_dp*SQRT(5.0_dp)*xy
          dmat(5,1,2) = -SQRT(15.0_dp)*yz
          dmat(5,2,2) = -SQRT(15.0_dp)*xz
          dmat(5,3,2) = zero
          dmat(5,4,2) = SQRT(3.75_dp)*zz
          dmat(6,0,2) = -SQRT(1.25_dp)*(3.0_dp*xx-3.0_dp*yy)
          dmat(6,1,2) = SQRT(15.0_dp)*xz
          dmat(6,2,2) = -SQRT(15.0_dp)*yz
          dmat(6,3,2) = -SQRT(3.75_dp)*zz
          dmat(6,4,2) = zero

          dmat(0,0,3) = SQRT(0.5_dp)*(12.0_dp*zz-3.0_dp)
          dmat(0,1,3) = 2.0_dp*SQRT(1.5_dp)*xz
          dmat(0,2,3) = 2.0_dp*SQRT(1.5_dp)*yz
          dmat(0,3,3) = -SQRT(6.0_dp)*(xx-yy)
          dmat(0,4,3) = -2.0_dp*SQRT(6.0_dp)*xy
          dmat(1,0,3) = 3.0_dp*SQRT(3.0_dp)*xz
          dmat(1,1,3) = 0.5_dp*(xx-5.0_dp*yy)
          dmat(1,2,3) = 3.0_dp*xy
          dmat(1,3,3) = 5.0_dp*xz
          dmat(1,4,3) = 5.0_dp*yz
          dmat(2,0,3) = 3.0_dp*SQRT(3.0_dp)*yz
          dmat(2,1,3) = 3.0_dp*xy
          dmat(2,2,3) = 0.5_dp*(yy-5.0_dp*xx)
          dmat(2,3,3) = -5.0_dp*yz
          dmat(2,4,3) = 5.0_dp*xz
          dmat(3,0,3) = zero
          dmat(3,1,3) = 2.0_dp*SQRT(2.5_dp)*xz
          dmat(3,2,3) = -2.0_dp*SQRT(2.5_dp)*yz
          dmat(3,3,3) = SQRT(2.5_dp)*(1.0_dp-6.0_dp*zz)
          dmat(3,4,3) = zero
          dmat(4,0,3) = zero
          dmat(4,1,3) = zero
          dmat(4,2,3) = zero
          dmat(4,3,3) = zero
          dmat(4,4,3) = SQRT(2.5_dp)*(1.0_dp-6.0_dp*zz)
          dmat(5,0,3) = zero
          dmat(5,1,3) = SQRT(3.75_dp)*(xx-yy)
          dmat(5,2,3) = -SQRT(15.0_dp)*xy
          dmat(5,3,3) = -SQRT(15.0_dp)*xz
          dmat(5,4,3) = SQRT(15.0_dp)*yz
          dmat(6,0,3) = zero
          dmat(6,1,3) = SQRT(15.0_dp)*xy
          dmat(6,2,3) = SQRT(3.75_dp)*(xx-yy)
          dmat(6,3,3) = -SQRT(15.0_dp)*yz
          dmat(6,4,3) = -SQRT(15.0_dp)*xz

        ELSE IF (l1==5 .AND. l2==7 .AND. num==3) THEN
          xx = x*x
          yy = y*y
          zz = z*z
          xy = x*y
          xz = x*z
          yz = y*z
! d3(d,f)
          dmat(0,0,1) = zero
          dmat(1,0,1) = SQRT(1.5_dp)*zz
          dmat(2,0,1) = zero
          dmat(3,0,1) = -2.0_dp*SQRT(6.0_dp)*xz
          dmat(4,0,1) = -2.0_dp*SQRT(6.0_dp)*yz
          dmat(0,1,1) = SQRT(0.75_dp)*(3.0_dp*zz-1.0_dp)
          dmat(1,1,1) = xz
          dmat(2,1,1) = 3.0_dp*yz
          dmat(3,1,1) = 2.5_dp*zz-3.0_dp*xx+yy
          dmat(4,1,1) = -4.0_dp*xy
          dmat(0,2,1) = zero
          dmat(1,2,1) = 3.0_dp*yz
          dmat(2,2,1) = -5.0_dp*xz
          dmat(3,2,1) = -2.0_dp*xy
          dmat(4,2,1) = 0.5_dp*(5.0_dp*zz-4.0_dp*yy)
          dmat(0,3,1) = zero
          dmat(1,3,1) = SQRT(2.5_dp)*(zz-2.0_dp*yy)
          dmat(2,3,1) = SQRT(2.5_dp)*4.0_dp*xy
          dmat(3,3,1) = zero
          dmat(4,3,1) = zero
          dmat(0,4,1) = zero
          dmat(1,4,1) = zero
          dmat(2,4,1) = SQRT(2.5_dp)*(1.0_dp-6.0_dp*xx)
          dmat(3,4,1) = zero
          dmat(4,4,1) = zero
          dmat(0,5,1) = 3.0_dp*SQRT(1.25_dp)*(yy-xx)
          dmat(1,5,1) = SQRT(15.0_dp)*xz
          dmat(2,5,1) = -SQRT(15.0_dp)*yz
          dmat(3,5,1) = -SQRT(3.75_dp)*zz
          dmat(4,5,1) = zero
          dmat(0,6,1) = -SQRT(5.0_dp)*3.0_dp*xy
          dmat(1,6,1) = SQRT(15.0_dp)*yz
          dmat(2,6,1) = SQRT(15.0_dp)*xz
          dmat(3,6,1) = zero
          dmat(4,6,1) = -SQRT(3.75_dp)*zz

          dmat(0,0,2) = zero
          dmat(1,0,2) = zero
          dmat(2,0,2) = SQRT(1.5_dp)*zz
          dmat(3,0,2) = 2.0_dp*SQRT(6.0_dp)*yz
          dmat(4,0,2) = -2.0_dp*SQRT(6.0_dp)*xz
          dmat(0,1,2) = zero
          dmat(1,1,2) = -5.0_dp*yz
          dmat(2,1,2) = 3.0_dp*xz
          dmat(3,1,2) = 2.0_dp*xy
          dmat(4,1,2) = 0.5_dp*(5.0_dp*zz-4.0_dp*xx)
          dmat(0,2,2) = SQRT(0.75_dp)*(3.0_dp*zz-1.0_dp)
          dmat(1,2,2) = 3.0_dp*xz
          dmat(2,2,2) = yz
          dmat(3,2,2) = -(2.5_dp*zz+xx-3.0_dp*yy)
          dmat(4,2,2) = -4.0_dp*xy
          dmat(0,3,2) = zero
          dmat(1,3,2) = -SQRT(2.5_dp)*4.0_dp*xy
          dmat(2,3,2) = -SQRT(2.5_dp)*(zz-2.0_dp*xx)
          dmat(3,3,2) = zero
          dmat(4,3,2) = zero
          dmat(0,4,2) = zero
          dmat(1,4,2) = SQRT(2.5_dp)*(1.0_dp-6.0_dp*yy)
          dmat(2,4,2) = zero
          dmat(3,4,2) = zero
          dmat(4,4,2) = zero
          dmat(0,5,2) = SQRT(5.0_dp)*3.0_dp*xy
          dmat(1,5,2) = -SQRT(15.0_dp)*yz
          dmat(2,5,2) = -SQRT(15.0_dp)*xz
          dmat(3,5,2) = zero
          dmat(4,5,2) = SQRT(3.75_dp)*zz
          dmat(0,6,2) = -SQRT(1.25_dp)*(3.0_dp*xx-3.0_dp*yy)
          dmat(1,6,2) = SQRT(15.0_dp)*xz
          dmat(2,6,2) = -SQRT(15.0_dp)*yz
          dmat(3,6,2) = -SQRT(3.75_dp)*zz
          dmat(4,6,2) = zero

          dmat(0,0,3) = SQRT(0.5_dp)*(12.0_dp*zz-3.0_dp)
          dmat(1,0,3) = 2.0_dp*SQRT(1.5_dp)*xz
          dmat(2,0,3) = 2.0_dp*SQRT(1.5_dp)*yz
          dmat(3,0,3) = -SQRT(6.0_dp)*(xx-yy)
          dmat(4,0,3) = -2.0_dp*SQRT(6.0_dp)*xy
          dmat(0,1,3) = 3.0_dp*SQRT(3.0_dp)*xz
          dmat(1,1,3) = 0.5_dp*(xx-5.0_dp*yy)
          dmat(2,1,3) = 3.0_dp*xy
          dmat(3,1,3) = 5.0_dp*xz
          dmat(4,1,3) = 5.0_dp*yz
          dmat(0,2,3) = 3.0_dp*SQRT(3.0_dp)*yz
          dmat(1,2,3) = 3.0_dp*xy
          dmat(2,2,3) = 0.5_dp*(yy-5.0_dp*xx)
          dmat(3,2,3) = -5.0_dp*yz
          dmat(4,2,3) = 5.0_dp*xz
          dmat(0,3,3) = zero
          dmat(1,3,3) = 2.0_dp*SQRT(2.5_dp)*xz
          dmat(2,3,3) = -2.0_dp*SQRT(2.5_dp)*yz
          dmat(3,3,3) = SQRT(2.5_dp)*(1.0_dp-6.0_dp*zz)
          dmat(4,3,3) = zero
          dmat(0,4,3) = zero
          dmat(1,4,3) = zero
          dmat(2,4,3) = zero
          dmat(3,4,3) = zero
          dmat(4,4,3) = SQRT(2.5_dp)*(1.0_dp-6.0_dp*zz)
          dmat(0,5,3) = zero
          dmat(1,5,3) = SQRT(3.75_dp)*(xx-yy)
          dmat(2,5,3) = -SQRT(15.0_dp)*xy
          dmat(3,5,3) = -SQRT(15.0_dp)*xz
          dmat(4,5,3) = SQRT(15.0_dp)*yz
          dmat(0,6,3) = zero
          dmat(1,6,3) = SQRT(15.0_dp)*xy
          dmat(2,6,3) = SQRT(3.75_dp)*(xx-yy)
          dmat(3,6,3) = -SQRT(15.0_dp)*yz
          dmat(4,6,3) = -SQRT(15.0_dp)*xz
        ELSE IF (l1==7 .AND. l2==7 .AND. num==2) THEN
! d2(f,f)
          dmat(0,0,1) = zero
          dmat(0,1,1) = SQRT(0.24_dp)*z
          dmat(0,2,1) = zero
          dmat(0,3,1) = -2.0_dp*SQRT(0.6_dp)*x
          dmat(0,4,1) = -SQRT(2.4_dp)*y
          dmat(0,5,1) = zero
          dmat(0,6,1) = zero
          dmat(1,0,1) = dmat(0,1,1)
          dmat(2,0,1) = dmat(0,2,1)
          dmat(3,0,1) = dmat(0,3,1)
          dmat(4,0,1) = dmat(0,4,1)
          dmat(5,0,1) = dmat(0,5,1)
          dmat(6,0,1) = dmat(0,6,1)
          dmat(1,1,1) = zero
          dmat(1,2,1) = 1.2_dp*y
          dmat(1,3,1) = SQRT(0.9_dp)*z
          dmat(1,4,1) = zero
          dmat(1,5,1) = -2.0_dp*SQRT(0.15_dp)*x
          dmat(1,6,1) = -SQRT(0.6_dp)*y
          dmat(2,1,1) = dmat(1,2,1)
          dmat(3,1,1) = dmat(1,3,1)
          dmat(4,1,1) = dmat(1,4,1)
          dmat(5,1,1) = dmat(1,5,1)
          dmat(6,1,1) = dmat(1,6,1)
          dmat(2,2,1) = -2.4_dp*x
          dmat(2,3,1) = zero
          dmat(2,4,1) = SQRT(0.9_dp)*z
          dmat(2,5,1) = SQRT(0.6_dp)*y
          dmat(2,6,1) = -2.0_dp*SQRT(0.15_dp)*x
          dmat(3,2,1) = dmat(2,3,1)
          dmat(4,2,1) = dmat(2,4,1)
          dmat(5,2,1) = dmat(2,5,1)
          dmat(6,2,1) = dmat(2,6,1)
          dmat(3,3,1) = zero
          dmat(3,4,1) = zero
          dmat(3,5,1) = SQRT(1.5_dp)*z
          dmat(3,6,1) = zero
          dmat(4,3,1) = dmat(3,4,1)
          dmat(5,3,1) = dmat(3,5,1)
          dmat(6,3,1) = dmat(3,6,1)
          dmat(4,4,1) = zero
          dmat(4,5,1) = zero
          dmat(4,6,1) = SQRT(1.5_dp)*z
          dmat(5,4,1) = dmat(4,5,1)
          dmat(6,4,1) = dmat(4,6,1)
          dmat(5,5,1) = zero
          dmat(5,6,1) = zero
          dmat(6,5,1) = zero
          dmat(6,6,1) = zero

          dmat(0,0,2) = zero
          dmat(0,1,2) = zero
          dmat(0,2,2) = SQRT(0.24_dp)*z
          dmat(0,3,2) = 2.0_dp*SQRT(0.6_dp)*y
          dmat(0,4,2) = -SQRT(2.4_dp)*x
          dmat(0,5,2) = zero
          dmat(0,6,2) = zero
          dmat(1,0,2) = dmat(0,1,2)
          dmat(2,0,2) = dmat(0,2,2)
          dmat(3,0,2) = dmat(0,3,2)
          dmat(4,0,2) = dmat(0,4,2)
          dmat(5,0,2) = dmat(0,5,2)
          dmat(6,0,2) = dmat(0,6,2)
          dmat(1,1,2) = -2.4_dp*y
          dmat(1,2,2) = 1.2_dp*x
          dmat(1,3,2) = zero
          dmat(1,4,2) = SQRT(0.9_dp)*z
          dmat(1,5,2) = 2.0_dp*SQRT(0.15_dp)*y
          dmat(1,6,2) = -SQRT(0.6_dp)*x
          dmat(2,1,2) = dmat(1,2,2)
          dmat(3,1,2) = dmat(1,3,2)
          dmat(4,1,2) = dmat(1,4,2)
          dmat(5,1,2) = dmat(1,5,2)
          dmat(6,1,2) = dmat(1,6,2)
          dmat(2,2,2) = zero
          dmat(2,3,2) = -SQRT(0.9_dp)*z
          dmat(2,4,2) = zero
          dmat(2,5,2) = SQRT(0.6_dp)*x
          dmat(2,6,2) = 2.0_dp*SQRT(0.15_dp)*y
          dmat(3,2,2) = dmat(2,3,2)
          dmat(4,2,2) = dmat(2,4,2)
          dmat(5,2,2) = dmat(2,5,2)
          dmat(6,2,2) = dmat(2,6,2)
          dmat(3,3,2) = zero
          dmat(3,4,2) = zero
          dmat(3,5,2) = zero
          dmat(3,6,2) = SQRT(1.5_dp)*z
          dmat(4,3,2) = dmat(3,4,2)
          dmat(5,3,2) = dmat(3,5,2)
          dmat(6,3,2) = dmat(3,6,2)
          dmat(4,4,2) = zero
          dmat(4,5,2) = -SQRT(1.5_dp)*z
          dmat(4,6,2) = zero
          dmat(5,4,2) = dmat(4,5,2)
          dmat(6,4,2) = dmat(4,6,2)
          dmat(5,5,2) = zero
          dmat(5,6,2) = zero
          dmat(6,5,2) = zero
          dmat(6,6,2) = zero

          dmat(0,0,3) = 2.4_dp*z
          dmat(0,1,3) = SQRT(0.24_dp)*x
          dmat(0,2,3) = SQRT(0.24_dp)*y
          dmat(0,3,3) = zero
          dmat(0,4,3) = zero
          dmat(0,5,3) = zero
          dmat(0,6,3) = zero
          dmat(1,0,3) = dmat(0,1,3)
          dmat(2,0,3) = dmat(0,2,3)
          dmat(3,0,3) = dmat(0,3,3)
          dmat(4,0,3) = dmat(0,4,3)
          dmat(5,0,3) = dmat(0,5,3)
          dmat(6,0,3) = dmat(0,6,3)
          dmat(1,1,3) = 0.6_dp*z
          dmat(1,2,3) = zero
          dmat(1,3,3) = SQRT(0.9_dp)*x
          dmat(1,4,3) = SQRT(0.9_dp)*y
          dmat(1,5,3) = zero
          dmat(1,6,3) = zero
          dmat(2,1,3) = dmat(1,2,3)
          dmat(3,1,3) = dmat(1,3,3)
          dmat(4,1,3) = dmat(1,4,3)
          dmat(5,1,3) = dmat(1,5,3)
          dmat(6,1,3) = dmat(1,6,3)
          dmat(2,2,3) = 0.6_dp*z
          dmat(2,3,3) = -SQRT(0.9_dp)*y
          dmat(2,4,3) = SQRT(0.9_dp)*x
          dmat(2,5,3) = zero
          dmat(2,6,3) = zero
          dmat(3,2,3) = dmat(2,3,3)
          dmat(4,2,3) = dmat(2,4,3)
          dmat(5,2,3) = dmat(2,5,3)
          dmat(6,2,3) = dmat(2,6,3)
          dmat(3,3,3) = zero
          dmat(3,4,3) = zero
          dmat(3,5,3) = SQRT(1.5_dp)*x
          dmat(3,6,3) = SQRT(1.5_dp)*y
          dmat(4,3,3) = dmat(3,4,3)
          dmat(5,3,3) = dmat(3,5,3)
          dmat(6,3,3) = dmat(3,6,3)
          dmat(4,4,3) = zero
          dmat(4,5,3) = -SQRT(1.5_dp)*y
          dmat(4,6,3) = SQRT(1.5_dp)*x
          dmat(5,4,3) = dmat(4,5,3)
          dmat(6,4,3) = dmat(4,6,3)
          dmat(5,5,3) = -3.0_dp*z
          dmat(5,6,3) = zero
          dmat(6,5,3) = zero
          dmat(6,6,3) = -3.0_dp*z
        ELSE IF (l1==7 .AND. l2==7 .AND. num==4) THEN
          xx = x*x
          yy = y*y
          zz = z*z
          xy = x*y
          xz = x*z
          yz = y*z
! d4(f,f)
          dmat(0,0,1) = zero
          dmat(0,1,1) = SQRT(0.24_dp)*(5.0_dp*zz-2.0_dp)*z
          dmat(0,2,1) = zero
          dmat(0,3,1) = -SQRT(0.6_dp)*2.0_dp*x*zz
          dmat(0,4,1) = -SQRT(2.4_dp)*y*zz
          dmat(0,5,1) = SQRT(3.6_dp)*(3.0_dp*yy-3.0_dp*xx)*z
          dmat(0,6,1) = -SQRT(3.6_dp)*6.0_dp*x*yz
          dmat(1,0,1) = dmat(0,1,1)
          dmat(2,0,1) = dmat(0,2,1)
          dmat(3,0,1) = dmat(0,3,1)
          dmat(4,0,1) = dmat(0,4,1)
          dmat(5,0,1) = dmat(0,5,1)
          dmat(6,0,1) = dmat(0,6,1)
          dmat(1,1,1) = -0.8_dp*x + 3.0_dp*x*zz
          dmat(1,2,1) = 0.4_dp*(10.0_dp*zz-1.0_dp)*y
          dmat(1,3,1) = SQRT(0.1_dp)*(6.0_dp*zz-8.0_dp*yy-1.0_dp)*z
          dmat(1,4,1) = zero
          dmat(1,5,1) = SQRT(3.75_dp)*(2.0_dp*x-5.6_dp*x*xx+2.4_dp*x*yy)
          dmat(1,6,1) = SQRT(0.6_dp)*(4.0_dp*zz-12.0_dp*xx+1.0_dp)*y
          dmat(2,1,1) = dmat(1,2,1)
          dmat(3,1,1) = dmat(1,3,1)
          dmat(4,1,1) = dmat(1,4,1)
          dmat(5,1,1) = dmat(1,5,1)
          dmat(6,1,1) = dmat(1,6,1)
          dmat(2,2,1) = -5.0_dp*x*zz
          dmat(2,3,1) = SQRT(0.1_dp)*16.0_dp*x*yz
          dmat(2,4,1) = SQRT(0.1_dp)*(2.0_dp*zz-24.0_dp*xx+3.0_dp)*z
          dmat(2,5,1) = SQRT(0.6_dp)*(4.0_dp*yy-4.0_dp*zz-1.0_dp)*y
          dmat(2,6,1) = SQRT(3.75_dp)*(2.0_dp*x-4.0_dp*x*xx-2.4_dp*x*yy)
          dmat(3,2,1) = dmat(2,3,1)
          dmat(4,2,1) = dmat(2,4,1)
          dmat(5,2,1) = dmat(2,5,1)
          dmat(6,2,1) = dmat(2,6,1)
          dmat(3,3,1) = -8.0_dp*x*yy
          dmat(3,4,1) = 2.0_dp*(3.0_dp*xx-yy)*y
          dmat(3,5,1) = -SQRT(1.5_dp)*(2.0_dp*zz-1.0_dp)*z
          dmat(3,6,1) = zero
          dmat(4,3,1) = dmat(3,4,1)
          dmat(5,3,1) = dmat(3,5,1)
          dmat(6,3,1) = dmat(3,6,1)
          dmat(4,4,1) = 8.0_dp*x*yy
          dmat(4,5,1) = zero
          dmat(4,6,1) = -SQRT(1.5_dp)*(2.0_dp*zz-1.0_dp)*z
          dmat(5,4,1) = dmat(4,5,1)
          dmat(6,4,1) = dmat(4,6,1)
          dmat(5,5,1) = zero
          dmat(5,6,1) = zero
          dmat(6,5,1) = zero
          dmat(6,6,1) = zero

          dmat(0,0,2) = zero
          dmat(0,1,2) = zero
          dmat(0,2,2) = SQRT(0.24_dp)*(5.0_dp*zz-2.0_dp)*z
          dmat(0,3,2) = SQRT(0.6_dp)*2.0_dp*y*zz
          dmat(0,4,2) = -SQRT(2.4_dp)*x*zz
          dmat(0,5,2) = SQRT(3.6_dp)*6.0_dp*y*xz
          dmat(0,6,2) = SQRT(3.6_dp)*(3.0_dp*yy-3.0_dp*xx)*z
          dmat(1,0,2) = dmat(0,1,2)
          dmat(2,0,2) = dmat(0,2,2)
          dmat(3,0,2) = dmat(0,3,2)
          dmat(4,0,2) = dmat(0,4,2)
          dmat(5,0,2) = dmat(0,5,2)
          dmat(6,0,2) = dmat(0,6,2)
          dmat(1,1,2) = -5.0_dp*y*zz
          dmat(1,2,2) = 0.4_dp*(10.0_dp*zz-1.0_dp)*x
          dmat(1,3,2) = SQRT(0.1_dp)*(-16.0_dp*y)*xz
          dmat(1,4,2) = SQRT(0.1_dp)*(2.0_dp*zz-24.0_dp*yy+3.0_dp)*z
          dmat(1,5,2) = SQRT(3.75_dp)*(-2.0_dp*y+2.4_dp*xx*y+4.0_dp*y*yy)
          dmat(1,6,2) = SQRT(0.6_dp)*(4.0_dp*zz-4.0_dp*xx+1.0_dp)*x
          dmat(2,1,2) = dmat(1,2,2)
          dmat(3,1,2) = dmat(1,3,2)
          dmat(4,1,2) = dmat(1,4,2)
          dmat(5,1,2) = dmat(1,5,2)
          dmat(6,1,2) = dmat(1,6,2)
          dmat(2,2,2) = -0.8_dp*y + 3.0_dp*y*zz
          dmat(2,3,2) = -SQRT(0.1_dp)*(6.0_dp*zz-8.0_dp*xx-1.0_dp)*z
          dmat(2,4,2) = zero
          dmat(2,5,2) = SQRT(0.6_dp)*(12.0_dp*yy-4.0_dp*zz-1.0_dp)*x
          dmat(2,6,2) = SQRT(3.75_dp)*(-2.0_dp*y-2.4_dp*xx*y+5.6_dp*y*yy)
          dmat(3,2,2) = dmat(2,3,2)
          dmat(4,2,2) = dmat(2,4,2)
          dmat(5,2,2) = dmat(2,5,2)
          dmat(6,2,2) = dmat(2,6,2)
          dmat(3,3,2) = -8.0_dp*xx*y
          dmat(3,4,2) = 2.0_dp*(xx-3.0_dp*yy)*x
          dmat(3,5,2) = zero
          dmat(3,6,2) = -SQRT(1.5_dp)*(2.0_dp*zz-1.0_dp)*z
          dmat(4,3,2) = dmat(3,4,2)
          dmat(5,3,2) = dmat(3,5,2)
          dmat(6,3,2) = dmat(3,6,2)
          dmat(4,4,2) = 8.0_dp*xx*y
          dmat(4,5,2) = SQRT(1.5_dp)*(2.0_dp*zz-1.0_dp)*z
          dmat(4,6,2) = zero
          dmat(5,4,2) = dmat(4,5,2)
          dmat(6,4,2) = dmat(4,6,2)
          dmat(5,5,2) = zero
          dmat(5,6,2) = zero
          dmat(6,5,2) = zero
          dmat(6,6,2) = zero

          dmat(0,0,3) = 0.6_dp*(20.0_dp*zz-8.0_dp)*z
          dmat(0,1,3) = SQRT(0.24_dp)*(15.0_dp*zz-2.0_dp)*x
          dmat(0,2,3) = SQRT(0.24_dp)*(15.0_dp*zz-2.0_dp)*y
          dmat(0,3,3) = -SQRT(0.6_dp)*(xx-yy)*2.0_dp*z
          dmat(0,4,3) = -SQRT(2.4_dp)*xy*2.0_dp*z
          dmat(0,5,3) = SQRT(3.6_dp)*(3.0_dp*yy-xx)*x
          dmat(0,6,3) = SQRT(3.6_dp)*(yy-3.0_dp*xx)*y
          dmat(1,0,3) = dmat(0,1,3)
          dmat(2,0,3) = dmat(0,2,3)
          dmat(3,0,3) = dmat(0,3,3)
          dmat(4,0,3) = dmat(0,4,3)
          dmat(5,0,3) = dmat(0,5,3)
          dmat(6,0,3) = dmat(0,6,3)
          dmat(1,1,3) = -(5.0_dp*yy-3.0_dp*xx)*z
          dmat(1,2,3) = 8.0_dp*z*xy
          dmat(1,3,3) = SQRT(0.1_dp)*(18.0_dp*zz-8.0_dp*yy-1.0_dp)*x
          dmat(1,4,3) = SQRT(0.1_dp)*(6.0_dp*zz-8.0_dp*yy+3.0_dp)*y
          dmat(1,5,3) = zero
          dmat(1,6,3) = SQRT(0.6_dp)*8.0_dp*z*xy
          dmat(2,1,3) = dmat(1,2,3)
          dmat(3,1,3) = dmat(1,3,3)
          dmat(4,1,3) = dmat(1,4,3)
          dmat(5,1,3) = dmat(1,5,3)
          dmat(6,1,3) = dmat(1,6,3)
          dmat(2,2,3) = -(5.0_dp*xx-3.0_dp*yy)*z
          dmat(2,3,3) = -SQRT(0.1_dp)*(18.0_dp*zz-8.0_dp*xx-1.0_dp)*y
          dmat(2,4,3) = SQRT(0.1_dp)*(6.0_dp*zz-8.0_dp*xx+3.0_dp)*x
          dmat(2,5,3) = SQRT(0.6_dp)*(-8.0_dp*z)*xy
          dmat(2,6,3) = zero
          dmat(3,2,3) = dmat(2,3,3)
          dmat(4,2,3) = dmat(2,4,3)
          dmat(5,2,3) = dmat(2,5,3)
          dmat(6,2,3) = dmat(2,6,3)
          dmat(3,3,3) = (4.0_dp-12.0_dp*zz)*z
          dmat(3,4,3) = zero
          dmat(3,5,3) = -SQRT(1.5_dp)*(6.0_dp*zz-1.0_dp)*x
          dmat(3,6,3) = -SQRT(1.5_dp)*(6.0_dp*zz-1.0_dp)*y
          dmat(4,3,3) = dmat(3,4,3)
          dmat(5,3,3) = dmat(3,5,3)
          dmat(6,3,3) = dmat(3,6,3)
          dmat(4,4,3) = -2.0_dp*(2.0_dp*zz-1.0_dp)*4.0_dp*z
          dmat(4,5,3) = SQRT(1.5_dp)*(6.0_dp*zz-1.0_dp)*y
          dmat(4,6,3) = -SQRT(1.5_dp)*(6.0_dp*zz-1.0_dp)*x
          dmat(5,4,3) = dmat(4,5,3)
          dmat(6,4,3) = dmat(4,6,3)
          dmat(5,5,3) = 1.5_dp*(4.0_dp*zz-2.0_dp)*z
          dmat(5,6,3) = zero
          dmat(6,5,3) = zero
          dmat(6,6,3) = 1.5_dp*(4.0_dp*zz-2.0_dp)*z
        END IF
        dc(0:l1-1,0:l2-1) = r(1)*dmat(:,:,1) + r(2)*dmat(:,:,2) + &
                            r(3)*dmat(:,:,3)
        dmat(:,:,1) = dmat(:,:,1) - r(1)*dc(0:l1-1,0:l2-1)
        dmat(:,:,2) = dmat(:,:,2) - r(2)*dc(0:l1-1,0:l2-1)
        dmat(:,:,3) = dmat(:,:,3) - r(3)*dc(0:l1-1,0:l2-1)
      END SUBROUTINE set_dmat
!------------------------------------------------------------------------------!
    END MODULE slater_koster_matr
!------------------------------------------------------------------------------!

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
  USE kinds,                           ONLY: dbl
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
      REAL (dbl), PARAMETER :: zero = 0._dbl
      REAL (dbl), PARAMETER :: one = 1._dbl

!!*****
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
      SUBROUTINE gmat(l1,l2,r,gmu,dll)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(dbl), INTENT(IN)                    :: r(3)
    REAL(dbl), INTENT(OUT)                   :: gmu(0:6,0:6,0:3), &
                                                dll(0:6,0:6,0:6)

    INTEGER                                  :: i, m1, m2, mo, mu
    REAL(dbl)                                :: gsl(0:6,0:3)

        mu = min(l1,l2)
        mo = max(l1,l2)
        m1 = 2*l1
        m2 = 2*l2
        gsl(0:6,0:3) = zero
        DO i = 0, mo
          CALL sph(i,r,gsl(0:6,i))
        END DO
        dll(0:6,0:6,0) = zero
        DO i = 0, 2*mo
          dll(i,i,0) = 1._dbl
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
              gmu(0:2,0:4,1) = -2._dbl/sqrt(3._dbl)*dll(0:2,0:4,3) + &
                dll(0:2,0:4,1)
            CASE (3)
!p-f
              CALL out_prod(dll(0:2,0:6,4),gsl(0:2,1),gsl(0:6,3))
              CALL set_mat(dll(0:2,0:6,2),r,2)
              gmu(0:2,0:6,0) = dll(0:2,0:6,4)
              gmu(0:2,0:6,1) = -sqrt(1.5_dbl)*dll(0:2,0:6,4) + &
                0.5_dbl*dll(0:2,0:6,2)
            END SELECT
          CASE (2)
            SELECT CASE (l2)
            CASE (1)
!d-p
              CALL out_prod(dll(0:4,0:2,3),gsl(0:4,2),gsl(0:2,1))
              CALL set_mat(dll(0:4,0:2,1),r,1)
              gmu(0:4,0:2,0) = dll(0:4,0:2,3)
              gmu(0:4,0:2,1) = -2._dbl/sqrt(3._dbl)*dll(0:4,0:2,3) + &
                dll(0:4,0:2,1)
            CASE (2)
!d-d
              CALL out_prod(dll(0:4,0:4,4),gsl(0:4,2),gsl(0:4,2))
              CALL set_mat(dll(0:4,0:4,2),r,2)
              gmu(0:4,0:4,0) = dll(0:4,0:4,4)
              gmu(0:4,0:4,1) = -4._dbl/3._dbl*dll(0:4,0:4,4) + &
                dll(0:4,0:4,2) + dll(0:4,0:4,0)
              gmu(0:4,0:4,2) = 1._dbl/3._dbl*dll(0:4,0:4,4) - dll(0:4,0:4,2)
            CASE (3)
!d-f
              CALL out_prod(dll(0:4,0:6,5),gsl(0:4,2),gsl(0:6,3))
              CALL set_mat(dll(0:4,0:6,3),r,3)
              CALL set_mat(dll(0:4,0:6,1),r,1)
              gmu(0:4,0:6,0) = dll(0:4,0:6,5)
              gmu(0:4,0:6,1) = -sqrt(2._dbl)*dll(0:4,0:6,5) + &
                0.5_dbl*dll(0:4,0:6,3) + 0.5_dbl*dll(0:4,0:6,1)
              gmu(0:4,0:6,2) = sqrt(0.2_dbl)*dll(0:4,0:6,5) - &
                sqrt(0.4_dbl)*dll(0:4,0:6,3)
            END SELECT
          CASE (3)
            SELECT CASE (l2)
            CASE (1)
!f-p
              CALL out_prod(dll(0:6,0:2,4),gsl(0:6,3),gsl(0:2,1))
              CALL set_mat(dll(0:6,0:2,2),r,2)
              gmu(0:6,0:2,0) = dll(0:6,0:2,4)
              gmu(0:6,0:2,1) = -sqrt(1.5_dbl)*dll(0:6,0:2,4) + &
                0.5_dbl*dll(0:6,0:2,2)
            CASE (2)
!f-d
              CALL out_prod(dll(0:6,0:4,5),gsl(0:6,3),gsl(0:4,2))
              CALL set_mat(dll(0:6,0:4,3),r,3)
              CALL set_mat(dll(0:6,0:4,1),r,1)
              gmu(0:6,0:4,0) = dll(0:6,0:4,5)
              gmu(0:6,0:4,1) = -sqrt(2._dbl)*dll(0:6,0:4,5) + &
                0.5_dbl*dll(0:6,0:4,3) + 0.5_dbl*dll(0:6,0:4,1)
              gmu(0:6,0:4,2) = sqrt(0.2_dbl)*dll(0:6,0:4,5) - &
                sqrt(0.4_dbl)*dll(0:6,0:4,3)
            CASE (3)
!f-f
              CALL out_prod(dll(0:6,0:6,6),gsl(0:6,3),gsl(0:6,3))
              CALL set_mat(dll(0:6,0:6,4),r,4)
              CALL set_mat(dll(0:6,0:6,2),r,2)
              gmu(0:6,0:6,0) = dll(0:6,0:6,6)
              gmu(0:6,0:6,1) = -1.5_dbl*dll(0:6,0:6,6) + &
                0.625_dbl*dll(0:6,0:6,4) + 0.625_dbl*dll(0:6,0:6,2) + &
                0.625_dbl*dll(0:6,0:6,0)
              gmu(0:6,0:6,2) = 0.6_dbl*dll(0:6,0:6,6) - dll(0:6,0:6,4)
              gmu(0:6,0:6,3) = -0.1_dbl*dll(0:6,0:6,6) + &
                0.375_dbl*dll(0:6,0:6,4) - 0.625_dbl*dll(0:6,0:6,2) + &
                0.375_dbl*dll(0:6,0:6,0)
            END SELECT
          END SELECT
        END IF

      END SUBROUTINE gmat
!------------------------------------------------------------------------------!
      SUBROUTINE dgmat(l1,l2,r,dgmu,ddll)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(dbl), INTENT(IN)                    :: r(3)
    REAL(dbl), INTENT(OUT)                   :: dgmu(0:6,0:6,0:3,1:3), &
                                                ddll(0:6,0:6,0:6,1:3)

    INTEGER                                  :: i, m1, m2, mo, mu
    REAL(dbl)                                :: dgsl(0:6,0:3,1:3), &
                                                gsl(0:6,0:3)

        mu = min(l1,l2)
        mo = max(l1,l2)
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
              dgmu(0:2,0:4,1,1:3) = -2._dbl/sqrt(3._dbl)*ddll(0:2,0:4,3,1:3) + &
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
              dgmu(0:2,0:6,1,1:3) = -sqrt(1.5_dbl)*ddll(0:2,0:6,4,1:3) + &
                0.5_dbl*ddll(0:2,0:6,2,1:3)
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
              dgmu(0:4,0:2,1,1:3) = -2._dbl/sqrt(3._dbl)*ddll(0:4,0:2,3,1:3) + &
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
              dgmu(0:4,0:4,1,1:3) = -4._dbl/3._dbl*ddll(0:4,0:4,4,1:3) + &
                ddll(0:4,0:4,2,1:3)
              dgmu(0:4,0:4,2,1:3) = 1._dbl/3._dbl*ddll(0:4,0:4,4,1:3) -  &
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
              dgmu(0:4,0:6,1,1:3) = -sqrt(2._dbl)*ddll(0:4,0:6,5,1:3) + &
                0.5_dbl*ddll(0:4,0:6,3,1:3) + 0.5_dbl*ddll(0:4,0:6,1,1:3)
              dgmu(0:4,0:6,2,1:3) = sqrt(0.2_dbl)*ddll(0:4,0:6,5,1:3) - &
                sqrt(0.4_dbl)*ddll(0:4,0:6,3,1:3)
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
              dgmu(0:6,0:2,1,1:3) = -sqrt(1.5_dbl)*ddll(0:6,0:2,4,1:3) + &
                0.5_dbl*ddll(0:6,0:2,2,1:3)
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
              dgmu(0:6,0:4,1,1:3) = -sqrt(2._dbl)*ddll(0:6,0:4,5,1:3) + &
                0.5_dbl*ddll(0:6,0:4,3,1:3) + 0.5_dbl*ddll(0:6,0:4,1,1:3)
              dgmu(0:6,0:4,2,1:3) = sqrt(0.2_dbl)*ddll(0:6,0:4,5,1:3) - &
                sqrt(0.4_dbl)*ddll(0:6,0:4,3,1:3)
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
              dgmu(0:6,0:6,1,1:3) = -1.5_dbl*ddll(0:6,0:6,6,1:3) + &
                0.625_dbl*ddll(0:6,0:6,4,1:3) + 0.625_dbl*ddll(0:6,0:6,2,1:3)
              dgmu(0:6,0:6,2,1:3) = 0.6_dbl*ddll(0:6,0:6,6,1:3) - ddll(0:6,0:6,4,1:3)
              dgmu(0:6,0:6,3,1:3) = -0.1_dbl*ddll(0:6,0:6,6,1:3) + &
                0.375_dbl*ddll(0:6,0:6,4,1:3) - 0.625_dbl*ddll(0:6,0:6,2,1:3)
            END SELECT
          END SELECT
        END IF

      END SUBROUTINE dgmat
!------------------------------------------------------------------------------!
      SUBROUTINE set_mat(mat,r,num)
    REAL(dbl), INTENT(OUT)                   :: mat(0:,0:)
    REAL(dbl), INTENT(IN)                    :: r(3)
    INTEGER, INTENT(IN)                      :: num

    REAL(dbl), PARAMETER                     :: s3 = 0.577350269189626_dbl 

    INTEGER                                  :: l1, l2
    REAL(dbl)                                :: x, xx, xy, xz, y, yy, yz, z, &
                                                zz

! sqrt(1._dbl/3._dbl)

        l1 = size(mat,1)
        l2 = size(mat,2)
        x = r(1)
        y = r(2)
        z = r(3)
        IF (l1==5 .AND. l2==3 .AND. num==1) THEN
! d1(d,p)
          mat(0:4,0:2) = zero
          mat(0,0) = 2._dbl*s3*z
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
          mat(0,0) = 2._dbl*s3*z
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
            mat(0,0) = zz - 2._dbl/3._dbl
            mat(0,1) = s3*xz
            mat(0,2) = s3*yz
            mat(0,3) = -s3*(xx-yy)
            mat(0,4) = -2._dbl*s3*xy
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
            mat(0,0) = sqrt(1.5_dbl)*(3._dbl*zz-1._dbl)
            mat(0,1) = -sqrt(6._dbl)*xz
            mat(0,2) = -sqrt(6._dbl)*yz
            mat(1,0) = 4._dbl*xz
            mat(1,1) = 0.5_dbl*(5._dbl*zz-2._dbl*xx-1._dbl)
            mat(1,2) = -xy
            mat(2,0) = 4._dbl*yz
            mat(2,1) = -xy
            mat(2,2) = 0.5_dbl*(5._dbl*zz-2._dbl*yy-1._dbl)
            mat(3,0) = sqrt(2.5_dbl)*(xx-yy)
            mat(3,1) = sqrt(10._dbl)*xz
            mat(3,2) = -sqrt(10._dbl)*yz
            mat(4,0) = sqrt(10._dbl)*xy
            mat(4,1) = sqrt(10._dbl)*yz
            mat(4,2) = sqrt(10._dbl)*xz
            mat(5,0) = zero
            mat(5,1) = sqrt(3.75_dbl)*(xx-yy)
            mat(5,2) = -sqrt(15._dbl)*xy
            mat(6,0) = zero
            mat(6,1) = sqrt(15._dbl)*xy
            mat(6,2) = sqrt(3.75_dbl)*(xx-yy)
          ELSE IF (l1==3 .AND. l2==7 .AND. num==2) THEN
! d2(p,f)
            mat(0,0) = sqrt(1.5_dbl)*(3._dbl*zz-1._dbl)
            mat(1,0) = -sqrt(6._dbl)*xz
            mat(2,0) = -sqrt(6._dbl)*yz
            mat(0,1) = 4._dbl*xz
            mat(1,1) = 0.5_dbl*(5._dbl*zz-2._dbl*xx-1._dbl)
            mat(2,1) = -xy
            mat(0,2) = 4._dbl*yz
            mat(1,2) = -xy
            mat(2,2) = 0.5_dbl*(5._dbl*zz-2._dbl*yy-1._dbl)
            mat(0,3) = sqrt(2.5_dbl)*(xx-yy)
            mat(1,3) = sqrt(10._dbl)*xz
            mat(2,3) = -sqrt(10._dbl)*yz
            mat(0,4) = sqrt(10._dbl)*xy
            mat(1,4) = sqrt(10._dbl)*yz
            mat(2,4) = sqrt(10._dbl)*xz
            mat(0,5) = zero
            mat(1,5) = sqrt(3.75_dbl)*(xx-yy)
            mat(2,5) = -sqrt(15._dbl)*xy
            mat(0,6) = zero
            mat(1,6) = sqrt(15._dbl)*xy
            mat(2,6) = sqrt(3.75_dbl)*(xx-yy)
          ELSE IF (l1==7 .AND. l2==5 .AND. num==1) THEN
! d1(f,d)
            mat(0,0) = sqrt(4.5_dbl)*z
            mat(0,1) = -sqrt(1.5_dbl)*x
            mat(0,2) = -sqrt(1.5_dbl)*y
            mat(0,3) = zero
            mat(0,4) = zero
            mat(1,0) = sqrt(3._dbl)*x
            mat(1,1) = 2._dbl*z
            mat(1,2) = zero
            mat(1,3) = -0.5_dbl*x
            mat(1,4) = -0.5_dbl*y
            mat(2,0) = sqrt(3._dbl)*y
            mat(2,1) = zero
            mat(2,2) = 2._dbl*z
            mat(2,3) = 0.5_dbl*y
            mat(2,4) = -0.5_dbl*x
            mat(3,0) = zero
            mat(3,1) = sqrt(2.5_dbl)*x
            mat(3,2) = -sqrt(2.5_dbl)*y
            mat(3,3) = sqrt(2.5_dbl)*z
            mat(3,4) = zero
            mat(4,0) = zero
            mat(4,1) = sqrt(2.5_dbl)*y
            mat(4,2) = sqrt(2.5_dbl)*x
            mat(4,3) = zero
            mat(4,4) = sqrt(2.5_dbl)*z
            mat(5,0) = zero
            mat(5,1) = zero
            mat(5,2) = zero
            mat(5,3) = sqrt(3.75_dbl)*x
            mat(5,4) = -sqrt(3.75_dbl)*y
            mat(6,0) = zero
            mat(6,1) = zero
            mat(6,2) = zero
            mat(6,3) = sqrt(3.75_dbl)*y
            mat(6,4) = sqrt(3.75_dbl)*x
          ELSE IF (l1==5 .AND. l2==7 .AND. num==1) THEN
! d1(d,f)
            mat(0,0) = sqrt(4.5_dbl)*z
            mat(1,0) = -sqrt(1.5_dbl)*x
            mat(2,0) = -sqrt(1.5_dbl)*y
            mat(3,0) = zero
            mat(4,0) = zero
            mat(0,1) = sqrt(3._dbl)*x
            mat(1,1) = 2._dbl*z
            mat(2,1) = zero
            mat(3,1) = -0.5_dbl*x
            mat(4,1) = -0.5_dbl*y
            mat(0,2) = sqrt(3._dbl)*y
            mat(1,2) = zero
            mat(2,2) = 2._dbl*z
            mat(3,2) = 0.5_dbl*y
            mat(4,2) = -0.5_dbl*x
            mat(0,3) = zero
            mat(1,3) = sqrt(2.5_dbl)*x
            mat(2,3) = -sqrt(2.5_dbl)*y
            mat(3,3) = sqrt(2.5_dbl)*z
            mat(4,3) = zero
            mat(0,4) = zero
            mat(1,4) = sqrt(2.5_dbl)*y
            mat(2,4) = sqrt(2.5_dbl)*x
            mat(3,4) = zero
            mat(4,4) = sqrt(2.5_dbl)*z
            mat(0,5) = zero
            mat(1,5) = zero
            mat(2,5) = zero
            mat(3,5) = sqrt(3.75_dbl)*x
            mat(4,5) = -sqrt(3.75_dbl)*y
            mat(0,6) = zero
            mat(1,6) = zero
            mat(2,6) = zero
            mat(3,6) = sqrt(3.75_dbl)*y
            mat(4,6) = sqrt(3.75_dbl)*x
          ELSE IF (l1==7 .AND. l2==5 .AND. num==3) THEN
! d3(f,d)
            mat(0,0) = sqrt(0.5_dbl)*(4._dbl*zz-3._dbl)*z
            mat(0,1) = sqrt(1.5_dbl)*x*zz
            mat(0,2) = sqrt(1.5_dbl)*y*zz
            mat(0,3) = -sqrt(6._dbl)*(xx-yy)*z
            mat(0,4) = -2._dbl*sqrt(6._dbl)*x*yz
            mat(1,0) = sqrt(0.75_dbl)*(3._dbl*zz-1._dbl)*x
            mat(1,1) = 0.5_dbl*(xx-5._dbl*yy)*z
            mat(1,2) = 3._dbl*x*yz
            mat(1,3) = (2.5_dbl*zz-xx+yy)*x
            mat(1,4) = 0.5_dbl*(5._dbl*zz-4._dbl*xx)*y
            mat(2,0) = sqrt(0.75_dbl)*(3._dbl*zz-1._dbl)*y
            mat(2,1) = 3._dbl*x*yz
            mat(2,2) = 0.5_dbl*(yy-5._dbl*xx)*z
            mat(2,3) = -(2.5_dbl*zz+xx-yy)*y
            mat(2,4) = 0.5_dbl*(5._dbl*zz-4._dbl*yy)*x
            mat(3,0) = zero
            mat(3,1) = sqrt(2.5_dbl)*(zz-2._dbl*yy)*x
            mat(3,2) = -sqrt(2.5_dbl)*(zz-2._dbl*xx)*y
            mat(3,3) = sqrt(2.5_dbl)*(1._dbl-2._dbl*zz)*z
            mat(3,4) = zero
            mat(4,0) = zero
            mat(4,1) = sqrt(2.5_dbl)*(1._dbl-2._dbl*yy)*y
            mat(4,2) = sqrt(2.5_dbl)*(1._dbl-2._dbl*xx)*x
            mat(4,3) = zero
            mat(4,4) = sqrt(2.5_dbl)*(1._dbl-2._dbl*zz)*z
            mat(5,0) = sqrt(1.25_dbl)*(3._dbl*yy-xx)*x
            mat(5,1) = sqrt(3.75_dbl)*(xx-yy)*z
            mat(5,2) = -sqrt(15._dbl)*x*yz
            mat(5,3) = -sqrt(3.75_dbl)*x*zz
            mat(5,4) = sqrt(3.75_dbl)*y*zz
            mat(6,0) = -sqrt(1.25_dbl)*(3._dbl*xx-yy)*y
            mat(6,1) = sqrt(15._dbl)*x*yz
            mat(6,2) = sqrt(3.75_dbl)*(xx-yy)*z
            mat(6,3) = -sqrt(3.75_dbl)*y*zz
            mat(6,4) = -sqrt(3.75_dbl)*x*zz
          ELSE IF (l1==5 .AND. l2==7 .AND. num==3) THEN
! d3(d,f)
            mat(0,0) = sqrt(0.5_dbl)*(4._dbl*zz-3._dbl)*z
            mat(1,0) = sqrt(1.5_dbl)*x*zz
            mat(2,0) = sqrt(1.5_dbl)*y*zz
            mat(3,0) = -sqrt(6._dbl)*(xx-yy)*z
            mat(4,0) = -2._dbl*sqrt(6._dbl)*x*yz
            mat(0,1) = sqrt(0.75_dbl)*(3._dbl*zz-1._dbl)*x
            mat(1,1) = 0.5_dbl*(xx-5._dbl*yy)*z
            mat(2,1) = 3._dbl*x*yz
            mat(3,1) = (2.5_dbl*zz-xx+yy)*x
            mat(4,1) = 0.5_dbl*(5._dbl*zz-4._dbl*xx)*y
            mat(0,2) = sqrt(0.75_dbl)*(3._dbl*zz-1._dbl)*y
            mat(1,2) = 3._dbl*x*yz
            mat(2,2) = 0.5_dbl*(yy-5._dbl*xx)*z
            mat(3,2) = -(2.5_dbl*zz+xx-yy)*y
            mat(4,2) = 0.5_dbl*(5._dbl*zz-4._dbl*yy)*x
            mat(0,3) = zero
            mat(1,3) = sqrt(2.5_dbl)*(zz-2._dbl*yy)*x
            mat(2,3) = -sqrt(2.5_dbl)*(zz-2._dbl*xx)*y
            mat(3,3) = sqrt(2.5_dbl)*(1._dbl-2._dbl*zz)*z
            mat(4,3) = zero
            mat(0,4) = zero
            mat(1,4) = sqrt(2.5_dbl)*(1._dbl-2._dbl*yy)*y
            mat(2,4) = sqrt(2.5_dbl)*(1._dbl-2._dbl*xx)*x
            mat(3,4) = zero
            mat(4,4) = sqrt(2.5_dbl)*(1._dbl-2._dbl*zz)*z
            mat(0,5) = sqrt(1.25_dbl)*(3._dbl*yy-xx)*x
            mat(1,5) = sqrt(3.75_dbl)*(xx-yy)*z
            mat(2,5) = -sqrt(15._dbl)*x*yz
            mat(3,5) = -sqrt(3.75_dbl)*x*zz
            mat(4,5) = sqrt(3.75_dbl)*y*zz
            mat(0,6) = -sqrt(1.25_dbl)*(3._dbl*xx-yy)*y
            mat(1,6) = sqrt(15._dbl)*x*yz
            mat(2,6) = sqrt(3.75_dbl)*(xx-yy)*z
            mat(3,6) = -sqrt(3.75_dbl)*y*zz
            mat(4,6) = -sqrt(3.75_dbl)*x*zz
          ELSE IF (l1==7 .AND. l2==7 .AND. num==2) THEN
! d2(f,f)
            mat(0,0) = 0.4_dbl*(3._dbl*zz-1._dbl)
            mat(0,1) = sqrt(0.24_dbl)*xz
            mat(0,2) = sqrt(0.24_dbl)*yz
            mat(0,3) = -sqrt(0.6_dbl)*(xx-yy)
            mat(0,4) = -sqrt(2.4_dbl)*xy
            mat(0,5) = zero
            mat(0,6) = zero
            mat(1,0) = mat(0,1)
            mat(2,0) = mat(0,2)
            mat(3,0) = mat(0,3)
            mat(4,0) = mat(0,4)
            mat(5,0) = mat(0,5)
            mat(6,0) = mat(0,6)
            mat(1,1) = 0.3_dbl*(1._dbl-4._dbl*yy+zz)
            mat(1,2) = 1.2_dbl*xy
            mat(1,3) = sqrt(0.9_dbl)*xz
            mat(1,4) = sqrt(0.9_dbl)*yz
            mat(1,5) = -sqrt(0.15_dbl)*(xx-yy)
            mat(1,6) = -sqrt(0.6_dbl)*xy
            mat(2,1) = mat(1,2)
            mat(3,1) = mat(1,3)
            mat(4,1) = mat(1,4)
            mat(5,1) = mat(1,5)
            mat(6,1) = mat(1,6)
            mat(2,2) = 0.3_dbl*(1._dbl-4._dbl*xx+zz)
            mat(2,3) = -sqrt(0.9_dbl)*yz
            mat(2,4) = sqrt(0.9_dbl)*xz
            mat(2,5) = sqrt(0.6_dbl)*xy
            mat(2,6) = -sqrt(0.15_dbl)*(xx-yy)
            mat(3,2) = mat(2,3)
            mat(4,2) = mat(2,4)
            mat(5,2) = mat(2,5)
            mat(6,2) = mat(2,6)
            mat(3,3) = zero
            mat(3,4) = zero
            mat(3,5) = sqrt(1.5_dbl)*xz
            mat(3,6) = sqrt(1.5_dbl)*yz
            mat(4,3) = mat(3,4)
            mat(5,3) = mat(3,5)
            mat(6,3) = mat(3,6)
            mat(4,4) = zero
            mat(4,5) = -sqrt(1.5_dbl)*yz
            mat(4,6) = sqrt(1.5_dbl)*xz
            mat(5,4) = mat(4,5)
            mat(6,4) = mat(4,6)
            mat(5,5) = -0.5_dbl*(3._dbl*zz-1._dbl)
            mat(5,6) = zero
            mat(6,5) = zero
            mat(6,6) = -0.5_dbl*(3._dbl*zz-1._dbl)
          ELSE IF (l1==7 .AND. l2==7 .AND. num==4) THEN
! d4(f,f)
            mat(0,0) = 0.6_dbl*(5._dbl*zz-4._dbl)*zz
            mat(0,1) = sqrt(0.24_dbl)*(5._dbl*zz-2._dbl)*xz
            mat(0,2) = sqrt(0.24_dbl)*(5._dbl*zz-2._dbl)*yz
            mat(0,3) = -sqrt(0.6_dbl)*(xx-yy)*zz
            mat(0,4) = -sqrt(2.4_dbl)*xy*zz
            mat(0,5) = sqrt(3.6_dbl)*(3._dbl*yy-xx)*xz
            mat(0,6) = sqrt(3.6_dbl)*(yy-3._dbl*xx)*yz
            mat(1,0) = mat(0,1)
            mat(2,0) = mat(0,2)
            mat(3,0) = mat(0,3)
            mat(4,0) = mat(0,4)
            mat(5,0) = mat(0,5)
            mat(6,0) = mat(0,6)
            mat(1,1) = -0.4_dbl*xx - 0.5_dbl*(5._dbl*yy-3._dbl*xx)*zz
            mat(1,2) = 0.4_dbl*(10._dbl*zz-1._dbl)*xy
            mat(1,3) = sqrt(0.1_dbl)*(6._dbl*zz-8._dbl*yy-1._dbl)*xz
            mat(1,4) = sqrt(0.1_dbl)*(2._dbl*zz-8._dbl*yy+3._dbl)*yz
            mat(1,5) = sqrt(3.75_dbl)*(xx-yy-1.4_dbl*xx*xx+1.2_dbl*xx*yy+yy*yy &
              )
            mat(1,6) = sqrt(0.6_dbl)*(4._dbl*zz-4._dbl*xx+1._dbl)*xy
            mat(2,1) = mat(1,2)
            mat(3,1) = mat(1,3)
            mat(4,1) = mat(1,4)
            mat(5,1) = mat(1,5)
            mat(6,1) = mat(1,6)
            mat(2,2) = -0.4_dbl*yy - 0.5_dbl*(5._dbl*xx-3._dbl*yy)*zz
            mat(2,3) = -sqrt(0.1_dbl)*(6._dbl*zz-8._dbl*xx-1._dbl)*yz
            mat(2,4) = sqrt(0.1_dbl)*(2._dbl*zz-8._dbl*xx+3._dbl)*xz
            mat(2,5) = sqrt(0.6_dbl)*(4._dbl*yy-4._dbl*zz-1._dbl)*xy
            mat(2,6) = sqrt(3.75_dbl)*(xx-yy-xx*xx-1.2_dbl*xx*yy+1.4_dbl*yy*yy &
              )
            mat(3,2) = mat(2,3)
            mat(4,2) = mat(2,4)
            mat(5,2) = mat(2,5)
            mat(6,2) = mat(2,6)
            mat(3,3) = (2._dbl-3._dbl*zz)*zz - 4._dbl*xx*yy
            mat(3,4) = 2._dbl*(xx-yy)*xy
            mat(3,5) = -sqrt(1.5_dbl)*(2._dbl*zz-1._dbl)*xz
            mat(3,6) = -sqrt(1.5_dbl)*(2._dbl*zz-1._dbl)*yz
            mat(4,3) = mat(3,4)
            mat(5,3) = mat(3,5)
            mat(6,3) = mat(3,6)
            mat(4,4) = -(2._dbl*zz-1._dbl)**2 + 4._dbl*xx*yy
            mat(4,5) = sqrt(1.5_dbl)*(2._dbl*zz-1._dbl)*yz
            mat(4,6) = -sqrt(1.5_dbl)*(2._dbl*zz-1._dbl)*xz
            mat(5,4) = mat(4,5)
            mat(6,4) = mat(4,6)
            mat(5,5) = 1.5_dbl*(zz-1._dbl)*zz
            mat(5,6) = zero
            mat(6,5) = zero
            mat(6,6) = 1.5_dbl*(zz-1._dbl)*zz
          END IF
        END IF
      END SUBROUTINE set_mat
!------------------------------------------------------------------------------!
      SUBROUTINE set_dmat(dmat,r,num)
    REAL(dbl), INTENT(OUT)                   :: dmat(0:,0:,1:)
    REAL(dbl), INTENT(IN)                    :: r(3)
    INTEGER, INTENT(IN)                      :: num

    REAL(dbl), PARAMETER                     :: s3 = 0.577350269189626_dbl 

    INTEGER                                  :: l1, l2
    REAL(dbl)                                :: dc(0:6,0:6), x, xx, xy, xz, &
                                                y, yy, yz, z, zz

! sqrt(1._dbl/3._dbl)

        l1 = size(dmat,1)
        l2 = size(dmat,2)
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
          dmat(0,0,3) = 2._dbl*s3
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
          dmat(0,0,3) = 2._dbl*s3
          dmat(1,1,3) = one
          dmat(2,2,3) = one
        ELSEIF (l1==5 .AND. l2==5 .AND. num==2) THEN
! d2(d,d)
          dmat(0:4,0:4,1:3) = zero
          dmat(0,1,1) = s3*z
          dmat(0,3,1) = -2._dbl*s3*x
          dmat(0,4,1) = -2._dbl*s3*y
          dmat(1,0,1) = dmat(0,1,1)
          dmat(3,0,1) = dmat(0,3,1)
          dmat(4,0,1) = dmat(0,4,1)
          dmat(1,2,1) = y
          dmat(1,3,1) = z
          dmat(2,1,1) = y
          dmat(3,1,1) = z
          dmat(2,2,1) = -2._dbl*x
          dmat(2,4,1) = z
          dmat(4,2,1) = z
          dmat(0,2,2) = s3*z
          dmat(0,3,2) = 2._dbl*s3*y
          dmat(0,4,2) = -2._dbl*s3*x
          dmat(2,0,2) = dmat(0,2,2)
          dmat(3,0,2) = dmat(0,3,2)
          dmat(4,0,2) = dmat(0,4,2)
          dmat(1,1,2) = -2._dbl*y
          dmat(1,2,2) = x
          dmat(1,4,2) = z
          dmat(2,1,2) = x
          dmat(4,1,2) = z
          dmat(2,3,2) = -z
          dmat(3,2,2) = -z
          dmat(0,0,3) = 2._dbl*z
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
          dmat(3,3,3) = -2._dbl*z
          dmat(4,4,3) = -2._dbl*z
        ELSE IF (l1==7 .AND. l2==3 .AND. num==2) THEN
! d2(f,p)
          dmat(0,0,1) = zero
          dmat(0,1,1) = -sqrt(6._dbl)*z
          dmat(0,2,1) = zero
          dmat(1,0,1) = 4._dbl*z
          dmat(1,1,1) = -2._dbl*x
          dmat(1,2,1) = -y
          dmat(2,0,1) = zero
          dmat(2,1,1) = -y
          dmat(2,2,1) = zero
          dmat(3,0,1) = 2._dbl*sqrt(2.5_dbl)*x
          dmat(3,1,1) = sqrt(10._dbl)*z
          dmat(3,2,1) = zero
          dmat(4,0,1) = sqrt(10._dbl)*y
          dmat(4,1,1) = zero
          dmat(4,2,1) = sqrt(10._dbl)*z
          dmat(5,0,1) = zero
          dmat(5,1,1) = sqrt(15._dbl)*x
          dmat(5,2,1) = -sqrt(15._dbl)*y
          dmat(6,0,1) = zero
          dmat(6,1,1) = sqrt(15._dbl)*y
          dmat(6,2,1) = sqrt(15._dbl)*x

          dmat(0,0,2) = zero
          dmat(0,1,2) = zero
          dmat(0,2,2) = -sqrt(6._dbl)*z
          dmat(1,0,2) = zero
          dmat(1,1,2) = zero
          dmat(1,2,2) = -x
          dmat(2,0,2) = 4._dbl*z
          dmat(2,1,2) = -x
          dmat(2,2,2) = -2._dbl*y
          dmat(3,0,2) = -2._dbl*sqrt(2.5_dbl)*y
          dmat(3,1,2) = zero
          dmat(3,2,2) = -sqrt(10._dbl)*z
          dmat(4,0,2) = sqrt(10._dbl)*x
          dmat(4,1,2) = sqrt(10._dbl)*z
          dmat(4,2,2) = zero
          dmat(5,0,2) = zero
          dmat(5,1,2) = -sqrt(15._dbl)*y
          dmat(5,2,2) = -sqrt(15._dbl)*x
          dmat(6,0,2) = zero
          dmat(6,1,2) = sqrt(15._dbl)*x
          dmat(6,2,2) = -sqrt(15._dbl)*y

          dmat(0,0,3) = 6._dbl*sqrt(1.5_dbl)*z
          dmat(0,1,3) = -sqrt(6._dbl)*x
          dmat(0,2,3) = -sqrt(6._dbl)*y
          dmat(1,0,3) = 4._dbl*x
          dmat(1,1,3) = 5._dbl*z
          dmat(1,2,3) = zero
          dmat(2,0,3) = 4._dbl*y
          dmat(2,1,3) = zero
          dmat(2,2,3) = 5._dbl*z
          dmat(3,0,3) = zero
          dmat(3,1,3) = sqrt(10._dbl)*x
          dmat(3,2,3) = -sqrt(10._dbl)*y
          dmat(4,0,3) = zero
          dmat(4,1,3) = sqrt(10._dbl)*y
          dmat(4,2,3) = sqrt(10._dbl)*x
          dmat(5,0,3) = zero
          dmat(5,1,3) = zero
          dmat(5,2,3) = zero
          dmat(6,0,3) = zero
          dmat(6,1,3) = zero
          dmat(6,2,3) = zero
        ELSE IF (l1==3 .AND. l2==7 .AND. num==2) THEN
! d2(p,f)
          dmat(0,0,1) = zero
          dmat(1,0,1) = -sqrt(6._dbl)*z
          dmat(2,0,1) = zero
          dmat(0,1,1) = 4._dbl*z
          dmat(1,1,1) = -2._dbl*x
          dmat(2,1,1) = -y
          dmat(0,2,1) = zero
          dmat(1,2,1) = -y
          dmat(2,2,1) = zero
          dmat(0,3,1) = 2._dbl*sqrt(2.5_dbl)*x
          dmat(1,3,1) = sqrt(10._dbl)*z
          dmat(2,3,1) = zero
          dmat(0,4,1) = sqrt(10._dbl)*y
          dmat(1,4,1) = zero
          dmat(2,4,1) = sqrt(10._dbl)*z
          dmat(0,5,1) = zero
          dmat(1,5,1) = sqrt(15._dbl)*x
          dmat(2,5,1) = -sqrt(15._dbl)*y
          dmat(0,6,1) = zero
          dmat(1,6,1) = sqrt(15._dbl)*y
          dmat(2,6,1) = sqrt(15._dbl)*x

          dmat(0,0,2) = zero
          dmat(1,0,2) = zero
          dmat(2,0,2) = -sqrt(6._dbl)*z
          dmat(0,1,2) = zero
          dmat(1,1,2) = zero
          dmat(2,1,2) = -x
          dmat(0,2,2) = 4._dbl*z
          dmat(1,2,2) = -x
          dmat(2,2,2) = -2._dbl*y
          dmat(0,3,2) = -2._dbl*sqrt(2.5_dbl)*y
          dmat(1,3,2) = zero
          dmat(2,3,2) = -sqrt(10._dbl)*z
          dmat(0,4,2) = sqrt(10._dbl)*x
          dmat(1,4,2) = sqrt(10._dbl)*z
          dmat(2,4,2) = zero
          dmat(0,5,2) = zero
          dmat(1,5,2) = -sqrt(15._dbl)*y
          dmat(2,5,2) = -sqrt(15._dbl)*x
          dmat(0,6,2) = zero
          dmat(1,6,2) = sqrt(15._dbl)*x
          dmat(2,6,2) = -sqrt(15._dbl)*y

          dmat(0,0,3) = 6._dbl*sqrt(1.5_dbl)*z
          dmat(1,0,3) = -sqrt(6._dbl)*x
          dmat(2,0,3) = -sqrt(6._dbl)*y
          dmat(0,1,3) = 4._dbl*x
          dmat(1,1,3) = 5._dbl*z
          dmat(2,1,3) = zero
          dmat(0,2,3) = 4._dbl*y
          dmat(1,2,3) = zero
          dmat(2,2,3) = 5._dbl*z
          dmat(0,3,3) = zero
          dmat(1,3,3) = sqrt(10._dbl)*x
          dmat(2,3,3) = -sqrt(10._dbl)*y
          dmat(0,4,3) = zero
          dmat(1,4,3) = sqrt(10._dbl)*y
          dmat(2,4,3) = sqrt(10._dbl)*x
          dmat(0,5,3) = zero
          dmat(1,5,3) = zero
          dmat(2,5,3) = zero
          dmat(0,6,3) = zero
          dmat(1,6,3) = zero
          dmat(2,6,3) = zero
        ELSE IF (l1==7 .AND. l2==5 .AND. num==1) THEN
! d1(f,d)
          dmat(0:6,0:4,1:3) = zero
          dmat(0,1,1) = -sqrt(1.5_dbl)
          dmat(1,0,1) = sqrt(3._dbl)
          dmat(1,3,1) = -0.5_dbl
          dmat(2,4,1) = -0.5_dbl
          dmat(3,1,1) = sqrt(2.5_dbl)
          dmat(4,2,1) = sqrt(2.5_dbl)
          dmat(5,3,1) = sqrt(3.75_dbl)
          dmat(6,4,1) = sqrt(3.75_dbl)

          dmat(0,2,2) = -sqrt(1.5_dbl)
          dmat(1,4,2) = -0.5_dbl
          dmat(2,0,2) = sqrt(3._dbl)
          dmat(2,3,2) = 0.5_dbl
          dmat(3,2,2) = -sqrt(2.5_dbl)
          dmat(4,1,2) = sqrt(2.5_dbl)
          dmat(5,4,2) = -sqrt(3.75_dbl)
          dmat(6,3,2) = sqrt(3.75_dbl)

          dmat(0,0,3) = sqrt(4.5_dbl)
          dmat(1,1,3) = 2._dbl
          dmat(2,2,3) = 2._dbl
          dmat(3,3,3) = sqrt(2.5_dbl)
          dmat(4,4,3) = sqrt(2.5_dbl)
        ELSE IF (l1==5 .AND. l2==7 .AND. num==1) THEN
! d1(d,f)
          dmat(0:4,0:6,1:3) = zero
          dmat(1,0,1) = -sqrt(1.5_dbl)
          dmat(0,1,1) = sqrt(3._dbl)
          dmat(3,1,1) = -0.5_dbl
          dmat(4,2,1) = -0.5_dbl
          dmat(1,3,1) = sqrt(2.5_dbl)
          dmat(2,4,1) = sqrt(2.5_dbl)
          dmat(3,5,1) = sqrt(3.75_dbl)
          dmat(4,6,1) = sqrt(3.75_dbl)

          dmat(2,0,2) = -sqrt(1.5_dbl)
          dmat(4,1,2) = -0.5_dbl
          dmat(0,2,2) = sqrt(3._dbl)
          dmat(3,2,2) = 0.5_dbl
          dmat(2,3,2) = -sqrt(2.5_dbl)
          dmat(1,4,2) = sqrt(2.5_dbl)
          dmat(4,5,2) = -sqrt(3.75_dbl)
          dmat(3,6,2) = sqrt(3.75_dbl)

          dmat(0,0,3) = sqrt(4.5_dbl)
          dmat(1,1,3) = 2._dbl
          dmat(2,2,3) = 2._dbl
          dmat(3,3,3) = sqrt(2.5_dbl)
          dmat(4,4,3) = sqrt(2.5_dbl)
        ELSE IF (l1==7 .AND. l2==5 .AND. num==3) THEN
          xx = x*x
          yy = y*y
          zz = z*z
          xy = x*y
          xz = x*z
          yz = y*z
! d3(f,d)
          dmat(0,0,1) = zero
          dmat(0,1,1) = sqrt(1.5_dbl)*zz
          dmat(0,2,1) = zero
          dmat(0,3,1) = -2._dbl*sqrt(6._dbl)*xz
          dmat(0,4,1) = -2._dbl*sqrt(6._dbl)*yz
          dmat(1,0,1) = sqrt(0.75_dbl)*(3._dbl*zz-1._dbl)
          dmat(1,1,1) = xz
          dmat(1,2,1) = 3._dbl*yz
          dmat(1,3,1) = 2.5_dbl*zz-3._dbl*xx+yy
          dmat(1,4,1) = -4._dbl*xy
          dmat(2,0,1) = zero
          dmat(2,1,1) = 3._dbl*yz
          dmat(2,2,1) = -5._dbl*xz
          dmat(2,3,1) = -2._dbl*xy
          dmat(2,4,1) = 0.5_dbl*(5._dbl*zz-4._dbl*yy)
          dmat(3,0,1) = zero
          dmat(3,1,1) = sqrt(2.5_dbl)*(zz-2._dbl*yy)
          dmat(3,2,1) = sqrt(2.5_dbl)*4._dbl*xy
          dmat(3,3,1) = zero
          dmat(3,4,1) = zero
          dmat(4,0,1) = zero
          dmat(4,1,1) = zero
          dmat(4,2,1) = sqrt(2.5_dbl)*(1._dbl-6._dbl*xx)
          dmat(4,3,1) = zero
          dmat(4,4,1) = zero
          dmat(5,0,1) = 3._dbl*sqrt(1.25_dbl)*(yy-xx)
          dmat(5,1,1) = sqrt(15._dbl)*xz
          dmat(5,2,1) = -sqrt(15._dbl)*yz
          dmat(5,3,1) = -sqrt(3.75_dbl)*zz
          dmat(5,4,1) = zero
          dmat(6,0,1) = -sqrt(5._dbl)*3._dbl*xy
          dmat(6,1,1) = sqrt(15._dbl)*yz
          dmat(6,2,1) = sqrt(15._dbl)*xz
          dmat(6,3,1) = zero
          dmat(6,4,1) = -sqrt(3.75_dbl)*zz

          dmat(0,0,2) = zero
          dmat(0,1,2) = zero
          dmat(0,2,2) = sqrt(1.5_dbl)*zz
          dmat(0,3,2) = 2._dbl*sqrt(6._dbl)*yz
          dmat(0,4,2) = -2._dbl*sqrt(6._dbl)*xz
          dmat(1,0,2) = zero
          dmat(1,1,2) = -5._dbl*yz
          dmat(1,2,2) = 3._dbl*xz
          dmat(1,3,2) = 2._dbl*xy
          dmat(1,4,2) = 0.5_dbl*(5._dbl*zz-4._dbl*xx)
          dmat(2,0,2) = sqrt(0.75_dbl)*(3._dbl*zz-1._dbl)
          dmat(2,1,2) = 3._dbl*xz
          dmat(2,2,2) = yz
          dmat(2,3,2) = -(2.5_dbl*zz+xx-3._dbl*yy)
          dmat(2,4,2) = -4._dbl*xy
          dmat(3,0,2) = zero
          dmat(3,1,2) = -sqrt(2.5_dbl)*4._dbl*xy
          dmat(3,2,2) = -sqrt(2.5_dbl)*(zz-2._dbl*xx)
          dmat(3,3,2) = zero
          dmat(3,4,2) = zero
          dmat(4,0,2) = zero
          dmat(4,1,2) = sqrt(2.5_dbl)*(1._dbl-6._dbl*yy)
          dmat(4,2,2) = zero
          dmat(4,3,2) = zero
          dmat(4,4,2) = zero
          dmat(5,0,2) = 3._dbl*sqrt(5._dbl)*xy
          dmat(5,1,2) = -sqrt(15._dbl)*yz
          dmat(5,2,2) = -sqrt(15._dbl)*xz
          dmat(5,3,2) = zero
          dmat(5,4,2) = sqrt(3.75_dbl)*zz
          dmat(6,0,2) = -sqrt(1.25_dbl)*(3._dbl*xx-3._dbl*yy)
          dmat(6,1,2) = sqrt(15._dbl)*xz
          dmat(6,2,2) = -sqrt(15._dbl)*yz
          dmat(6,3,2) = -sqrt(3.75_dbl)*zz
          dmat(6,4,2) = zero

          dmat(0,0,3) = sqrt(0.5_dbl)*(12._dbl*zz-3._dbl)
          dmat(0,1,3) = 2._dbl*sqrt(1.5_dbl)*xz
          dmat(0,2,3) = 2._dbl*sqrt(1.5_dbl)*yz
          dmat(0,3,3) = -sqrt(6._dbl)*(xx-yy)
          dmat(0,4,3) = -2._dbl*sqrt(6._dbl)*xy
          dmat(1,0,3) = 3._dbl*sqrt(3._dbl)*xz
          dmat(1,1,3) = 0.5_dbl*(xx-5._dbl*yy)
          dmat(1,2,3) = 3._dbl*xy
          dmat(1,3,3) = 5._dbl*xz
          dmat(1,4,3) = 5._dbl*yz
          dmat(2,0,3) = 3._dbl*sqrt(3._dbl)*yz
          dmat(2,1,3) = 3._dbl*xy
          dmat(2,2,3) = 0.5_dbl*(yy-5._dbl*xx)
          dmat(2,3,3) = -5._dbl*yz
          dmat(2,4,3) = 5._dbl*xz
          dmat(3,0,3) = zero
          dmat(3,1,3) = 2._dbl*sqrt(2.5_dbl)*xz
          dmat(3,2,3) = -2._dbl*sqrt(2.5_dbl)*yz
          dmat(3,3,3) = sqrt(2.5_dbl)*(1._dbl-6._dbl*zz)
          dmat(3,4,3) = zero
          dmat(4,0,3) = zero
          dmat(4,1,3) = zero
          dmat(4,2,3) = zero
          dmat(4,3,3) = zero
          dmat(4,4,3) = sqrt(2.5_dbl)*(1._dbl-6._dbl*zz)
          dmat(5,0,3) = zero
          dmat(5,1,3) = sqrt(3.75_dbl)*(xx-yy)
          dmat(5,2,3) = -sqrt(15._dbl)*xy
          dmat(5,3,3) = -sqrt(15._dbl)*xz
          dmat(5,4,3) = sqrt(15._dbl)*yz
          dmat(6,0,3) = zero
          dmat(6,1,3) = sqrt(15._dbl)*xy
          dmat(6,2,3) = sqrt(3.75_dbl)*(xx-yy)
          dmat(6,3,3) = -sqrt(15._dbl)*yz
          dmat(6,4,3) = -sqrt(15._dbl)*xz

        ELSE IF (l1==5 .AND. l2==7 .AND. num==3) THEN
          xx = x*x
          yy = y*y
          zz = z*z
          xy = x*y
          xz = x*z
          yz = y*z
! d3(d,f)
          dmat(0,0,1) = zero
          dmat(1,0,1) = sqrt(1.5_dbl)*zz
          dmat(2,0,1) = zero
          dmat(3,0,1) = -2._dbl*sqrt(6._dbl)*xz
          dmat(4,0,1) = -2._dbl*sqrt(6._dbl)*yz
          dmat(0,1,1) = sqrt(0.75_dbl)*(3._dbl*zz-1._dbl)
          dmat(1,1,1) = xz
          dmat(2,1,1) = 3._dbl*yz
          dmat(3,1,1) = 2.5_dbl*zz-3._dbl*xx+yy
          dmat(4,1,1) = -4._dbl*xy
          dmat(0,2,1) = zero
          dmat(1,2,1) = 3._dbl*yz
          dmat(2,2,1) = -5._dbl*xz
          dmat(3,2,1) = -2._dbl*xy
          dmat(4,2,1) = 0.5_dbl*(5._dbl*zz-4._dbl*yy)
          dmat(0,3,1) = zero
          dmat(1,3,1) = sqrt(2.5_dbl)*(zz-2._dbl*yy)
          dmat(2,3,1) = sqrt(2.5_dbl)*4._dbl*xy
          dmat(3,3,1) = zero
          dmat(4,3,1) = zero
          dmat(0,4,1) = zero
          dmat(1,4,1) = zero
          dmat(2,4,1) = sqrt(2.5_dbl)*(1._dbl-6._dbl*xx)
          dmat(3,4,1) = zero
          dmat(4,4,1) = zero
          dmat(0,5,1) = 3._dbl*sqrt(1.25_dbl)*(yy-xx)
          dmat(1,5,1) = sqrt(15._dbl)*xz
          dmat(2,5,1) = -sqrt(15._dbl)*yz
          dmat(3,5,1) = -sqrt(3.75_dbl)*zz
          dmat(4,5,1) = zero
          dmat(0,6,1) = -sqrt(5._dbl)*3._dbl*xy
          dmat(1,6,1) = sqrt(15._dbl)*yz
          dmat(2,6,1) = sqrt(15._dbl)*xz
          dmat(3,6,1) = zero
          dmat(4,6,1) = -sqrt(3.75_dbl)*zz

          dmat(0,0,2) = zero
          dmat(1,0,2) = zero
          dmat(2,0,2) = sqrt(1.5_dbl)*zz
          dmat(3,0,2) = 2._dbl*sqrt(6._dbl)*yz
          dmat(4,0,2) = -2._dbl*sqrt(6._dbl)*xz
          dmat(0,1,2) = zero
          dmat(1,1,2) = -5._dbl*yz
          dmat(2,1,2) = 3._dbl*xz
          dmat(3,1,2) = 2._dbl*xy
          dmat(4,1,2) = 0.5_dbl*(5._dbl*zz-4._dbl*xx)
          dmat(0,2,2) = sqrt(0.75_dbl)*(3._dbl*zz-1._dbl)
          dmat(1,2,2) = 3._dbl*xz
          dmat(2,2,2) = yz
          dmat(3,2,2) = -(2.5_dbl*zz+xx-3._dbl*yy)
          dmat(4,2,2) = -4._dbl*xy
          dmat(0,3,2) = zero
          dmat(1,3,2) = -sqrt(2.5_dbl)*4._dbl*xy
          dmat(2,3,2) = -sqrt(2.5_dbl)*(zz-2._dbl*xx)
          dmat(3,3,2) = zero
          dmat(4,3,2) = zero
          dmat(0,4,2) = zero
          dmat(1,4,2) = sqrt(2.5_dbl)*(1._dbl-6._dbl*yy)
          dmat(2,4,2) = zero
          dmat(3,4,2) = zero
          dmat(4,4,2) = zero
          dmat(0,5,2) = sqrt(5._dbl)*3._dbl*xy
          dmat(1,5,2) = -sqrt(15._dbl)*yz
          dmat(2,5,2) = -sqrt(15._dbl)*xz
          dmat(3,5,2) = zero
          dmat(4,5,2) = sqrt(3.75_dbl)*zz
          dmat(0,6,2) = -sqrt(1.25_dbl)*(3._dbl*xx-3._dbl*yy)
          dmat(1,6,2) = sqrt(15._dbl)*xz
          dmat(2,6,2) = -sqrt(15._dbl)*yz
          dmat(3,6,2) = -sqrt(3.75_dbl)*zz
          dmat(4,6,2) = zero

          dmat(0,0,3) = sqrt(0.5_dbl)*(12._dbl*zz-3._dbl)
          dmat(1,0,3) = 2._dbl*sqrt(1.5_dbl)*xz
          dmat(2,0,3) = 2._dbl*sqrt(1.5_dbl)*yz
          dmat(3,0,3) = -sqrt(6._dbl)*(xx-yy)
          dmat(4,0,3) = -2._dbl*sqrt(6._dbl)*xy
          dmat(0,1,3) = 3._dbl*sqrt(3._dbl)*xz
          dmat(1,1,3) = 0.5_dbl*(xx-5._dbl*yy)
          dmat(2,1,3) = 3._dbl*xy
          dmat(3,1,3) = 5._dbl*xz
          dmat(4,1,3) = 5._dbl*yz
          dmat(0,2,3) = 3._dbl*sqrt(3._dbl)*yz
          dmat(1,2,3) = 3._dbl*xy
          dmat(2,2,3) = 0.5_dbl*(yy-5._dbl*xx)
          dmat(3,2,3) = -5._dbl*yz
          dmat(4,2,3) = 5._dbl*xz
          dmat(0,3,3) = zero
          dmat(1,3,3) = 2._dbl*sqrt(2.5_dbl)*xz
          dmat(2,3,3) = -2._dbl*sqrt(2.5_dbl)*yz
          dmat(3,3,3) = sqrt(2.5_dbl)*(1._dbl-6._dbl*zz)
          dmat(4,3,3) = zero
          dmat(0,4,3) = zero
          dmat(1,4,3) = zero
          dmat(2,4,3) = zero
          dmat(3,4,3) = zero
          dmat(4,4,3) = sqrt(2.5_dbl)*(1._dbl-6._dbl*zz)
          dmat(0,5,3) = zero
          dmat(1,5,3) = sqrt(3.75_dbl)*(xx-yy)
          dmat(2,5,3) = -sqrt(15._dbl)*xy
          dmat(3,5,3) = -sqrt(15._dbl)*xz
          dmat(4,5,3) = sqrt(15._dbl)*yz
          dmat(0,6,3) = zero
          dmat(1,6,3) = sqrt(15._dbl)*xy
          dmat(2,6,3) = sqrt(3.75_dbl)*(xx-yy)
          dmat(3,6,3) = -sqrt(15._dbl)*yz
          dmat(4,6,3) = -sqrt(15._dbl)*xz
        ELSE IF (l1==7 .AND. l2==7 .AND. num==2) THEN
! d2(f,f)
          dmat(0,0,1) = zero
          dmat(0,1,1) = sqrt(0.24_dbl)*z
          dmat(0,2,1) = zero
          dmat(0,3,1) = -2._dbl*sqrt(0.6_dbl)*x
          dmat(0,4,1) = -sqrt(2.4_dbl)*y
          dmat(0,5,1) = zero
          dmat(0,6,1) = zero
          dmat(1,0,1) = dmat(0,1,1)
          dmat(2,0,1) = dmat(0,2,1)
          dmat(3,0,1) = dmat(0,3,1)
          dmat(4,0,1) = dmat(0,4,1)
          dmat(5,0,1) = dmat(0,5,1)
          dmat(6,0,1) = dmat(0,6,1)
          dmat(1,1,1) = zero
          dmat(1,2,1) = 1.2_dbl*y
          dmat(1,3,1) = sqrt(0.9_dbl)*z
          dmat(1,4,1) = zero
          dmat(1,5,1) = -2._dbl*sqrt(0.15_dbl)*x
          dmat(1,6,1) = -sqrt(0.6_dbl)*y
          dmat(2,1,1) = dmat(1,2,1)
          dmat(3,1,1) = dmat(1,3,1)
          dmat(4,1,1) = dmat(1,4,1)
          dmat(5,1,1) = dmat(1,5,1)
          dmat(6,1,1) = dmat(1,6,1)
          dmat(2,2,1) = -2.4_dbl*x
          dmat(2,3,1) = zero
          dmat(2,4,1) = sqrt(0.9_dbl)*z
          dmat(2,5,1) = sqrt(0.6_dbl)*y
          dmat(2,6,1) = -2._dbl*sqrt(0.15_dbl)*x
          dmat(3,2,1) = dmat(2,3,1)
          dmat(4,2,1) = dmat(2,4,1)
          dmat(5,2,1) = dmat(2,5,1)
          dmat(6,2,1) = dmat(2,6,1)
          dmat(3,3,1) = zero
          dmat(3,4,1) = zero
          dmat(3,5,1) = sqrt(1.5_dbl)*z
          dmat(3,6,1) = zero
          dmat(4,3,1) = dmat(3,4,1)
          dmat(5,3,1) = dmat(3,5,1)
          dmat(6,3,1) = dmat(3,6,1)
          dmat(4,4,1) = zero
          dmat(4,5,1) = zero
          dmat(4,6,1) = sqrt(1.5_dbl)*z
          dmat(5,4,1) = dmat(4,5,1)
          dmat(6,4,1) = dmat(4,6,1)
          dmat(5,5,1) = zero
          dmat(5,6,1) = zero
          dmat(6,5,1) = zero
          dmat(6,6,1) = zero

          dmat(0,0,2) = zero
          dmat(0,1,2) = zero
          dmat(0,2,2) = sqrt(0.24_dbl)*z
          dmat(0,3,2) = 2._dbl*sqrt(0.6_dbl)*y
          dmat(0,4,2) = -sqrt(2.4_dbl)*x
          dmat(0,5,2) = zero
          dmat(0,6,2) = zero
          dmat(1,0,2) = dmat(0,1,2)
          dmat(2,0,2) = dmat(0,2,2)
          dmat(3,0,2) = dmat(0,3,2)
          dmat(4,0,2) = dmat(0,4,2)
          dmat(5,0,2) = dmat(0,5,2)
          dmat(6,0,2) = dmat(0,6,2)
          dmat(1,1,2) = -2.4_dbl*y
          dmat(1,2,2) = 1.2_dbl*x
          dmat(1,3,2) = zero
          dmat(1,4,2) = sqrt(0.9_dbl)*z
          dmat(1,5,2) = 2._dbl*sqrt(0.15_dbl)*y
          dmat(1,6,2) = -sqrt(0.6_dbl)*x
          dmat(2,1,2) = dmat(1,2,2)
          dmat(3,1,2) = dmat(1,3,2)
          dmat(4,1,2) = dmat(1,4,2)
          dmat(5,1,2) = dmat(1,5,2)
          dmat(6,1,2) = dmat(1,6,2)
          dmat(2,2,2) = zero
          dmat(2,3,2) = -sqrt(0.9_dbl)*z
          dmat(2,4,2) = zero
          dmat(2,5,2) = sqrt(0.6_dbl)*x
          dmat(2,6,2) = 2._dbl*sqrt(0.15_dbl)*y
          dmat(3,2,2) = dmat(2,3,2)
          dmat(4,2,2) = dmat(2,4,2)
          dmat(5,2,2) = dmat(2,5,2)
          dmat(6,2,2) = dmat(2,6,2)
          dmat(3,3,2) = zero
          dmat(3,4,2) = zero
          dmat(3,5,2) = zero
          dmat(3,6,2) = sqrt(1.5_dbl)*z
          dmat(4,3,2) = dmat(3,4,2)
          dmat(5,3,2) = dmat(3,5,2)
          dmat(6,3,2) = dmat(3,6,2)
          dmat(4,4,2) = zero
          dmat(4,5,2) = -sqrt(1.5_dbl)*z
          dmat(4,6,2) = zero
          dmat(5,4,2) = dmat(4,5,2)
          dmat(6,4,2) = dmat(4,6,2)
          dmat(5,5,2) = zero
          dmat(5,6,2) = zero
          dmat(6,5,2) = zero
          dmat(6,6,2) = zero

          dmat(0,0,3) = 2.4_dbl*z
          dmat(0,1,3) = sqrt(0.24_dbl)*x
          dmat(0,2,3) = sqrt(0.24_dbl)*y
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
          dmat(1,1,3) = 0.6_dbl*z
          dmat(1,2,3) = zero
          dmat(1,3,3) = sqrt(0.9_dbl)*x
          dmat(1,4,3) = sqrt(0.9_dbl)*y
          dmat(1,5,3) = zero
          dmat(1,6,3) = zero
          dmat(2,1,3) = dmat(1,2,3)
          dmat(3,1,3) = dmat(1,3,3)
          dmat(4,1,3) = dmat(1,4,3)
          dmat(5,1,3) = dmat(1,5,3)
          dmat(6,1,3) = dmat(1,6,3)
          dmat(2,2,3) = 0.6_dbl*z
          dmat(2,3,3) = -sqrt(0.9_dbl)*y
          dmat(2,4,3) = sqrt(0.9_dbl)*x
          dmat(2,5,3) = zero
          dmat(2,6,3) = zero
          dmat(3,2,3) = dmat(2,3,3)
          dmat(4,2,3) = dmat(2,4,3)
          dmat(5,2,3) = dmat(2,5,3)
          dmat(6,2,3) = dmat(2,6,3)
          dmat(3,3,3) = zero
          dmat(3,4,3) = zero
          dmat(3,5,3) = sqrt(1.5_dbl)*x
          dmat(3,6,3) = sqrt(1.5_dbl)*y
          dmat(4,3,3) = dmat(3,4,3)
          dmat(5,3,3) = dmat(3,5,3)
          dmat(6,3,3) = dmat(3,6,3)
          dmat(4,4,3) = zero
          dmat(4,5,3) = -sqrt(1.5_dbl)*y
          dmat(4,6,3) = sqrt(1.5_dbl)*x
          dmat(5,4,3) = dmat(4,5,3)
          dmat(6,4,3) = dmat(4,6,3)
          dmat(5,5,3) = -3._dbl*z
          dmat(5,6,3) = zero
          dmat(6,5,3) = zero
          dmat(6,6,3) = -3._dbl*z
        ELSE IF (l1==7 .AND. l2==7 .AND. num==4) THEN
          xx = x*x
          yy = y*y
          zz = z*z
          xy = x*y
          xz = x*z
          yz = y*z
! d4(f,f)
          dmat(0,0,1) = zero
          dmat(0,1,1) = sqrt(0.24_dbl)*(5._dbl*zz-2._dbl)*z
          dmat(0,2,1) = zero
          dmat(0,3,1) = -sqrt(0.6_dbl)*2._dbl*x*zz
          dmat(0,4,1) = -sqrt(2.4_dbl)*y*zz
          dmat(0,5,1) = sqrt(3.6_dbl)*(3._dbl*yy-3._dbl*xx)*z
          dmat(0,6,1) = -sqrt(3.6_dbl)*6._dbl*x*yz
          dmat(1,0,1) = dmat(0,1,1)
          dmat(2,0,1) = dmat(0,2,1)
          dmat(3,0,1) = dmat(0,3,1)
          dmat(4,0,1) = dmat(0,4,1)
          dmat(5,0,1) = dmat(0,5,1)
          dmat(6,0,1) = dmat(0,6,1)
          dmat(1,1,1) = -0.8_dbl*x + 3._dbl*x*zz
          dmat(1,2,1) = 0.4_dbl*(10._dbl*zz-1._dbl)*y
          dmat(1,3,1) = sqrt(0.1_dbl)*(6._dbl*zz-8._dbl*yy-1._dbl)*z
          dmat(1,4,1) = zero
          dmat(1,5,1) = sqrt(3.75_dbl)*(2._dbl*x-5.6_dbl*x*xx+2.4_dbl*x*yy)
          dmat(1,6,1) = sqrt(0.6_dbl)*(4._dbl*zz-12._dbl*xx+1._dbl)*y
          dmat(2,1,1) = dmat(1,2,1)
          dmat(3,1,1) = dmat(1,3,1)
          dmat(4,1,1) = dmat(1,4,1)
          dmat(5,1,1) = dmat(1,5,1)
          dmat(6,1,1) = dmat(1,6,1)
          dmat(2,2,1) = -5._dbl*x*zz
          dmat(2,3,1) = sqrt(0.1_dbl)*16._dbl*x*yz
          dmat(2,4,1) = sqrt(0.1_dbl)*(2._dbl*zz-24._dbl*xx+3._dbl)*z
          dmat(2,5,1) = sqrt(0.6_dbl)*(4._dbl*yy-4._dbl*zz-1._dbl)*y
          dmat(2,6,1) = sqrt(3.75_dbl)*(2._dbl*x-4._dbl*x*xx-2.4_dbl*x*yy)
          dmat(3,2,1) = dmat(2,3,1)
          dmat(4,2,1) = dmat(2,4,1)
          dmat(5,2,1) = dmat(2,5,1)
          dmat(6,2,1) = dmat(2,6,1)
          dmat(3,3,1) = -8._dbl*x*yy
          dmat(3,4,1) = 2._dbl*(3._dbl*xx-yy)*y
          dmat(3,5,1) = -sqrt(1.5_dbl)*(2._dbl*zz-1._dbl)*z
          dmat(3,6,1) = zero
          dmat(4,3,1) = dmat(3,4,1)
          dmat(5,3,1) = dmat(3,5,1)
          dmat(6,3,1) = dmat(3,6,1)
          dmat(4,4,1) = 8._dbl*x*yy
          dmat(4,5,1) = zero
          dmat(4,6,1) = -sqrt(1.5_dbl)*(2._dbl*zz-1._dbl)*z
          dmat(5,4,1) = dmat(4,5,1)
          dmat(6,4,1) = dmat(4,6,1)
          dmat(5,5,1) = zero
          dmat(5,6,1) = zero
          dmat(6,5,1) = zero
          dmat(6,6,1) = zero

          dmat(0,0,2) = zero
          dmat(0,1,2) = zero
          dmat(0,2,2) = sqrt(0.24_dbl)*(5._dbl*zz-2._dbl)*z
          dmat(0,3,2) = sqrt(0.6_dbl)*2._dbl*y*zz
          dmat(0,4,2) = -sqrt(2.4_dbl)*x*zz
          dmat(0,5,2) = sqrt(3.6_dbl)*6._dbl*y*xz
          dmat(0,6,2) = sqrt(3.6_dbl)*(3._dbl*yy-3._dbl*xx)*z
          dmat(1,0,2) = dmat(0,1,2)
          dmat(2,0,2) = dmat(0,2,2)
          dmat(3,0,2) = dmat(0,3,2)
          dmat(4,0,2) = dmat(0,4,2)
          dmat(5,0,2) = dmat(0,5,2)
          dmat(6,0,2) = dmat(0,6,2)
          dmat(1,1,2) = -5._dbl*y*zz
          dmat(1,2,2) = 0.4_dbl*(10._dbl*zz-1._dbl)*x
          dmat(1,3,2) = sqrt(0.1_dbl)*(-16._dbl*y)*xz
          dmat(1,4,2) = sqrt(0.1_dbl)*(2._dbl*zz-24._dbl*yy+3._dbl)*z
          dmat(1,5,2) = sqrt(3.75_dbl)*(-2._dbl*y+2.4_dbl*xx*y+4._dbl*y*yy)
          dmat(1,6,2) = sqrt(0.6_dbl)*(4._dbl*zz-4._dbl*xx+1._dbl)*x
          dmat(2,1,2) = dmat(1,2,2)
          dmat(3,1,2) = dmat(1,3,2)
          dmat(4,1,2) = dmat(1,4,2)
          dmat(5,1,2) = dmat(1,5,2)
          dmat(6,1,2) = dmat(1,6,2)
          dmat(2,2,2) = -0.8_dbl*y + 3._dbl*y*zz
          dmat(2,3,2) = -sqrt(0.1_dbl)*(6._dbl*zz-8._dbl*xx-1._dbl)*z
          dmat(2,4,2) = zero
          dmat(2,5,2) = sqrt(0.6_dbl)*(12._dbl*yy-4._dbl*zz-1._dbl)*x
          dmat(2,6,2) = sqrt(3.75_dbl)*(-2._dbl*y-2.4_dbl*xx*y+5.6_dbl*y*yy)
          dmat(3,2,2) = dmat(2,3,2)
          dmat(4,2,2) = dmat(2,4,2)
          dmat(5,2,2) = dmat(2,5,2)
          dmat(6,2,2) = dmat(2,6,2)
          dmat(3,3,2) = -8._dbl*xx*y
          dmat(3,4,2) = 2._dbl*(xx-3._dbl*yy)*x
          dmat(3,5,2) = zero
          dmat(3,6,2) = -sqrt(1.5_dbl)*(2._dbl*zz-1._dbl)*z
          dmat(4,3,2) = dmat(3,4,2)
          dmat(5,3,2) = dmat(3,5,2)
          dmat(6,3,2) = dmat(3,6,2)
          dmat(4,4,2) = 8._dbl*xx*y
          dmat(4,5,2) = sqrt(1.5_dbl)*(2._dbl*zz-1._dbl)*z
          dmat(4,6,2) = zero
          dmat(5,4,2) = dmat(4,5,2)
          dmat(6,4,2) = dmat(4,6,2)
          dmat(5,5,2) = zero
          dmat(5,6,2) = zero
          dmat(6,5,2) = zero
          dmat(6,6,2) = zero

          dmat(0,0,3) = 0.6_dbl*(20._dbl*zz-8._dbl)*z
          dmat(0,1,3) = sqrt(0.24_dbl)*(15._dbl*zz-2._dbl)*x
          dmat(0,2,3) = sqrt(0.24_dbl)*(15._dbl*zz-2._dbl)*y
          dmat(0,3,3) = -sqrt(0.6_dbl)*(xx-yy)*2._dbl*z
          dmat(0,4,3) = -sqrt(2.4_dbl)*xy*2._dbl*z
          dmat(0,5,3) = sqrt(3.6_dbl)*(3._dbl*yy-xx)*x
          dmat(0,6,3) = sqrt(3.6_dbl)*(yy-3._dbl*xx)*y
          dmat(1,0,3) = dmat(0,1,3)
          dmat(2,0,3) = dmat(0,2,3)
          dmat(3,0,3) = dmat(0,3,3)
          dmat(4,0,3) = dmat(0,4,3)
          dmat(5,0,3) = dmat(0,5,3)
          dmat(6,0,3) = dmat(0,6,3)
          dmat(1,1,3) = -(5._dbl*yy-3._dbl*xx)*z
          dmat(1,2,3) = 8._dbl*z*xy
          dmat(1,3,3) = sqrt(0.1_dbl)*(18._dbl*zz-8._dbl*yy-1._dbl)*x
          dmat(1,4,3) = sqrt(0.1_dbl)*(6._dbl*zz-8._dbl*yy+3._dbl)*y
          dmat(1,5,3) = zero
          dmat(1,6,3) = sqrt(0.6_dbl)*8._dbl*z*xy
          dmat(2,1,3) = dmat(1,2,3)
          dmat(3,1,3) = dmat(1,3,3)
          dmat(4,1,3) = dmat(1,4,3)
          dmat(5,1,3) = dmat(1,5,3)
          dmat(6,1,3) = dmat(1,6,3)
          dmat(2,2,3) = -(5._dbl*xx-3._dbl*yy)*z
          dmat(2,3,3) = -sqrt(0.1_dbl)*(18._dbl*zz-8._dbl*xx-1._dbl)*y
          dmat(2,4,3) = sqrt(0.1_dbl)*(6._dbl*zz-8._dbl*xx+3._dbl)*x
          dmat(2,5,3) = sqrt(0.6_dbl)*(-8._dbl*z)*xy
          dmat(2,6,3) = zero
          dmat(3,2,3) = dmat(2,3,3)
          dmat(4,2,3) = dmat(2,4,3)
          dmat(5,2,3) = dmat(2,5,3)
          dmat(6,2,3) = dmat(2,6,3)
          dmat(3,3,3) = (4._dbl-12._dbl*zz)*z
          dmat(3,4,3) = zero
          dmat(3,5,3) = -sqrt(1.5_dbl)*(6._dbl*zz-1._dbl)*x
          dmat(3,6,3) = -sqrt(1.5_dbl)*(6._dbl*zz-1._dbl)*y
          dmat(4,3,3) = dmat(3,4,3)
          dmat(5,3,3) = dmat(3,5,3)
          dmat(6,3,3) = dmat(3,6,3)
          dmat(4,4,3) = -2._dbl*(2._dbl*zz-1._dbl)*4._dbl*z
          dmat(4,5,3) = sqrt(1.5_dbl)*(6._dbl*zz-1._dbl)*y
          dmat(4,6,3) = -sqrt(1.5_dbl)*(6._dbl*zz-1._dbl)*x
          dmat(5,4,3) = dmat(4,5,3)
          dmat(6,4,3) = dmat(4,6,3)
          dmat(5,5,3) = 1.5_dbl*(4._dbl*zz-2._dbl)*z
          dmat(5,6,3) = zero
          dmat(6,5,3) = zero
          dmat(6,6,3) = 1.5_dbl*(4._dbl*zz-2._dbl)*z
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

!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
!!****** cp2k/slater_koster [1.0] *
!!
!!   NAME
!!     slater_koster
!!
!!   FUNCTION
!!     Calculation of two-center s-f Slater-Koster integrals
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     A.K. McMahan, Phys. Rev. B, 58, p4293 (1998)
!!
!!   SOURCE
!******************************************************************************

    MODULE slater_koster
!------------------------------------------------------------------------------!
  USE kinds,                           ONLY: dbl
  USE slater_koster_matr,              ONLY: dgmat,&
                                             gmat
!


      IMPLICIT NONE
!
      PRIVATE
      PUBLIC :: m_integrals
!
!!*****
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!!****** slater_koster/m_integrals [1.0] *
!!
!!   NAME
!!     m_integrals
!!
!!   FUNCTION
!!     Calculation of two-center s-f Slater-Koster integrals
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     A.K. McMahan, Phys. Rev. B, 58, p4293 (1998)
!!
!!     Sign convention: t(lp,lq)=(-1)**(lp+lq)*(lp m,lq m)
!!     where (lp m,lq m) is the integral in the diatomic coordinate system
!!
!!   SOURCE
!******************************************************************************
      SUBROUTINE m_integrals(l1,l2,r,tll,hmm)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(dbl), INTENT(IN)                    :: r(3), tll(0:2*min(l1,l2))
    REAL(dbl), INTENT(OUT)                   :: hmm(0:2*l1,0:2*l2)

    INTEGER                                  :: i, m1, m2, mu
    REAL(dbl)                                :: dll(0:6,0:6,0:6), &
                                                gmu(0:6,0:6,0:3)

        IF (l1+l2==0) THEN
          hmm(0,0) = tll(0)
        ELSE IF (l1==0 .AND. l2==1) THEN
          hmm(0,0) = tll(0)*r(3)
          hmm(0,1) = tll(0)*r(1)
          hmm(0,2) = tll(0)*r(2)
        ELSE IF (l1==1 .AND. l2==0) THEN
          hmm(0,0) = tll(0)*r(3)
          hmm(1,0) = tll(0)*r(1)
          hmm(2,0) = tll(0)*r(2)
        ELSE IF (l1==1 .AND. l2==1) THEN
          hmm(0,0) = tll(0)*r(3)*r(3) + tll(1)*(1._dbl-r(3)*r(3))
          hmm(1,0) = (tll(0)-tll(1))*r(3)*r(1)
          hmm(2,0) = (tll(0)-tll(1))*r(3)*r(2)
          hmm(0,1) = hmm(1,0)
          hmm(0,2) = hmm(2,0)
          hmm(1,1) = tll(0)*r(1)*r(1) + tll(1)*(1._dbl-r(1)*r(1))
          hmm(2,1) = (tll(0)-tll(1))*r(1)*r(2)
          hmm(1,2) = hmm(2,1)
          hmm(2,2) = tll(0)*r(2)*r(2) + tll(1)*(1._dbl-r(2)*r(2))
        ELSE
          CALL gmat(l1,l2,r,gmu,dll)
          mu = min(l1,l2)
          m1 = 2*l1
          m2 = 2*l2
          hmm(0:m1,0:m2) = 0._dbl
          DO i = 0, mu
            hmm(0:m1,0:m2) = hmm(0:m1,0:m2) + gmu(0:m1,0:m2,i)*tll(i)
          END DO
        END IF
      END SUBROUTINE m_integrals
!!*****
!------------------------------------------------------------------------------!
      SUBROUTINE dm_integrals(l1,l2,r,tllr,dtll,pmm,grad)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(dbl), INTENT(IN)                    :: r(3), tllr(0:2*min(l1,l2)), &
                                                dtll(0:2*min(l1,l2)), &
                                                pmm(0:2*l1,0:2*l2)
    REAL(dbl), INTENT(OUT)                   :: grad(1:3)

    INTEGER                                  :: i, j, m1, m2, mu
    REAL(dbl)                                :: ddll(0:6,0:6,0:6,1:3), &
                                                dgmu(0:6,0:6,0:3,1:3), &
                                                hmm(0:6,0:6), ph

        mu = min(l1,l2)
        m1 = 2*l1
        m2 = 2*l2

        CALL m_integrals(l1,l2,r,dtll,hmm(0:m1,0:m2))
        ph = trace_AB(hmm(0:m1,0:m2),pmm(0:m1,0:m2))
        grad = ph * r

        IF (l1+l2==0) THEN
        ELSE
          CALL dgmat(l1,l2,r,dgmu,ddll)
          DO j = 1, 3
            hmm(0:m1,0:m2) = 0._dbl
            DO i = 0, mu
              hmm(0:m1,0:m2) = hmm(0:m1,0:m2) + dgmu(0:m1,0:m2,i,j)*tllr(i)
            END DO
            grad(j) = grad(j) + trace_AB(hmm(0:m1,0:m2),pmm(0:m1,0:m2))
          END DO
        END IF

      END SUBROUTINE dm_integrals
!------------------------------------------------------------------------------!
      FUNCTION trace_AB(A,B) RESULT(T)
    REAL(dbl), INTENT(IN)                    :: A(:,:), B(:,:)
    REAL(dbl)                                :: T

    INTEGER                                  :: i, m, n

        m = SIZE(A,1)
        n = SIZE(A,2)
        T = 0._dbl
        DO i = 1, m
          T = T + dot_product(A(:,i),B(:,i))
        END DO

      END FUNCTION trace_AB
!------------------------------------------------------------------------------!
    END MODULE slater_koster
!------------------------------------------------------------------------------!

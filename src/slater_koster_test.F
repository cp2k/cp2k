!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
!!****** cp2k/slater_koster_test [1.0] *
!!
!!   NAME
!!     slater_koster_test
!!
!!   FUNCTION
!!     Test for Slater-Koster matrices
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!   SOURCE
!!
!******************************************************************************
!
    MODULE slater_koster_test
!------------------------------------------------------------------------------!
  USE kinds,                           ONLY: dp
  USE slater_koster_matr,              ONLY: dgmat,&
                                             gmat
!

  USE slater_koster_util,              ONLY: dpro

      IMPLICIT NONE
!
      PRIVATE
      PUBLIC :: sk_test, dsk_test
!
!!***
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
      SUBROUTINE sk_test(l1,l2,emax)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(KIND=dp), INTENT(OUT)                   :: emax

    INTEGER                                  :: lm, lu
    REAL(KIND=dp) :: chi, dll1(0:6,0:6,0:6), dll2(0:6,0:6,0:6), &
      gmu1(0:6,0:6,0:3), gmu2(0:6,0:6,0:3), r1(3), r2(3), rr

        CALL random_number(r1)
        rr = dot_product(r1,r1)
        r1 = r1/sqrt(rr)
        CALL random_number(r2)
        rr = dot_product(r2,r2)
        r2 = r2/sqrt(rr)
        chi = dot_product(r1,r2)

        gmu1 = 0.0_dp
        gmu2 = 0.0_dp
        dll1 = 0.0_dp
        dll2 = 0.0_dp
        CALL gmat(l1,l2,r1,gmu1,dll1)
        CALL gmat(l1,l2,r2,gmu2,dll2)

        emax = 0.0_dp
        lm = max(l1,l2)
        lu = min(l1,l2)
        SELECT CASE (lu)
        CASE (0)
          CALL s_test(chi,gmu1,gmu2,l1,l2,emax)
        CASE (1)
          CALL p_test(chi,gmu1,gmu2,dll1,dll2,l1,l2,emax)
        CASE (2)
          CALL d_test(chi,gmu1,gmu2,dll1,dll2,l1,l2,emax)
        CASE (3)
          CALL f_test(chi,gmu1,gmu2,dll1,dll2,l1,l2,emax)
        END SELECT

      END SUBROUTINE sk_test
!------------------------------------------------------------------------------!
      SUBROUTINE dsk_test(l1,l2,emax)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(KIND=dp), INTENT(OUT)                   :: emax

    INTEGER                                  :: lm, lu
    REAL(KIND=dp) :: chi, dchi1(3), dchi2(3), ddll1(0:6,0:6,0:6,1:3), &
      ddll2(0:6,0:6,0:6,1:3), dgmu1(0:6,0:6,0:3,1:3), dgmu2(0:6,0:6,0:3,1:3), &
      dll1(0:6,0:6,0:6), dll2(0:6,0:6,0:6), gmu1(0:6,0:6,0:3), &
      gmu2(0:6,0:6,0:3), r1(3), r2(3), rr

        CALL random_number(r1)
        rr = dot_product(r1,r1)
        r1 = r1/sqrt(rr)
        CALL random_number(r2)
        rr = dot_product(r2,r2)
        r2 = r2/sqrt(rr)
        chi = dot_product(r1,r2)
        dchi1(1) = r2(1) - r1(1)*chi
        dchi1(2) = r2(2) - r1(2)*chi
        dchi1(3) = r2(3) - r1(3)*chi
        dchi2(1) = r1(1) - r2(1)*chi
        dchi2(2) = r1(2) - r2(2)*chi
        dchi2(3) = r1(3) - r2(3)*chi

        gmu1 = 0.0_dp
        gmu2 = 0.0_dp
        dll1 = 0.0_dp
        dll2 = 0.0_dp
        dgmu1 = 0.0_dp
        dgmu2 = 0.0_dp
        ddll1 = 0.0_dp
        ddll2 = 0.0_dp
        CALL gmat(l1,l2,r1,gmu1,dll1)
        CALL gmat(l1,l2,r2,gmu2,dll2)
        CALL dgmat(l1,l2,r1,dgmu1,ddll1)
        CALL dgmat(l1,l2,r2,dgmu2,ddll2)

        emax = 0.0_dp
        lm = max(l1,l2)
        lu = min(l1,l2)
        SELECT CASE (lu)
        CASE (0)
          CALL ds_test(chi,dchi1,dchi2,gmu1,gmu2,dgmu1,dgmu2,l1,l2,emax)
        CASE (1)
          CALL dp_test(chi,dchi1,dchi2,gmu1,gmu2,dgmu1,dgmu2, &
                         dll1,dll2,ddll1,ddll2,l1,l2,emax)
        CASE (2)
          CALL dd_test(chi,dchi1,dchi2,gmu1,gmu2,dgmu1,dgmu2, &
                         dll1,dll2,ddll1,ddll2,l1,l2,emax)
        CASE (3)
          CALL df_test(chi,dchi1,dchi2,gmu1,gmu2,dgmu1,dgmu2, &
                         dll1,dll2,ddll1,ddll2,l1,l2,emax)
        END SELECT

      END SUBROUTINE dsk_test
!------------------------------------------------------------------------------!
      SUBROUTINE s_test(chi,gmu1,gmu2,l1,l2,emax)
    REAL(KIND=dp), INTENT(IN)                    :: chi, gmu1(0:6,0:6,0:3), &
                                                gmu2(0:6,0:6,0:3)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(KIND=dp), INTENT(OUT)                   :: emax

    INTEGER                                  :: lm, lu
    REAL(KIND=dp)                                :: xt

        lm = max(l1,l2)
        lu = min(l1,l2)
        xt = dpro(gmu1(:,:,0),gmu2(:,:,0),l1,l2)
        SELECT CASE (lm)
        CASE (0)
          xt = abs(xt - 1.0_dp)
        CASE (1)
          xt = abs(xt - chi)
        CASE (2)
          xt = abs(xt - 0.5_dp*(3.0_dp*chi*chi-1.0_dp))
        CASE (3)
          xt = abs(xt - 0.5_dp*(5.0_dp*chi*chi-3.0_dp)*chi)
        END SELECT
        IF (xt>1.E-12) WRITE (*,*) ' Error in test for gmu ', l1, l2, xt
        emax = max(emax,xt)
      END SUBROUTINE s_test

      SUBROUTINE p_test(chi,gmu1,gmu2,dll1,dll2,l1,l2,emax)
    REAL(KIND=dp), INTENT(IN)                    :: chi, gmu1(0:6,0:6,0:3), &
                                                gmu2(0:6,0:6,0:3), &
                                                dll1(0:6,0:6,0:6), &
                                                dll2(0:6,0:6,0:6)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(KIND=dp), INTENT(OUT)                   :: emax

    INTEGER                                  :: lm, lu
    REAL(KIND=dp)                                :: dt1, dt2, dt3, xt1, xt2, xt3

        lm = max(l1,l2)
        lu = min(l1,l2)
        xt1 = dpro(gmu1(:,:,0),gmu2(:,:,0),l1,l2)
        xt2 = dpro(gmu1(:,:,0),gmu2(:,:,1),l1,l2)
        xt3 = dpro(gmu1(:,:,1),gmu2(:,:,1),l1,l2)
        SELECT CASE (lm)
        CASE (1)
          xt1 = abs(xt1 - chi*chi)
          xt2 = abs(xt2 - (1.0_dp-chi*chi))
          xt3 = abs(xt3 - (1.0_dp+chi*chi))
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 ', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 ', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 ', l1, l2, xt3
          emax = max(emax,xt1,xt2,xt3)
          dt1 = dpro(dll1(:,:,0),dll2(:,:,0),l1,l2)
          dt2 = dpro(dll1(:,:,0),dll2(:,:,2),l1,l2)
          dt3 = dpro(dll1(:,:,2),dll2(:,:,2),l1,l2)
          dt1 = abs(dt1-3.0_dp)
          dt2 = abs(dt2-1.0_dp)
          dt3 = abs(dt3-chi*chi)
          IF (dt1>1.E-12) WRITE (*,*) ' Error in test for dll #1 ', l1, l2, dt1
          IF (dt2>1.E-12) WRITE (*,*) ' Error in test for dll #2 ', l1, l2, dt2
          IF (dt3>1.E-12) WRITE (*,*) ' Error in test for dll #3 ', l1, l2, dt3
          emax = max(emax,dt1,dt2,dt3)
        CASE (2)
          xt1 = abs(xt1 - 0.5_dp*chi*(3.0_dp*chi*chi-1.0_dp))
          xt2 = abs(xt2 + sqrt(3.0_dp)*chi*(chi*chi-1.0_dp))
          xt3 = abs(xt3 - 2.0_dp*chi*chi*chi)
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 ', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 ', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 ', l1, l2, xt3
          emax = max(emax,xt1,xt2,xt3)
          dt1 = dpro(dll1(:,:,1),dll2(:,:,1),l1,l2)
          dt2 = dpro(dll1(:,:,1),dll2(:,:,3),l1,l2)
          dt3 = dpro(dll1(:,:,3),dll2(:,:,3),l1,l2)
          dt1 = abs(dt1-10.0_dp/3.0_dp*chi)
          dt2 = abs(dt2-2.0_dp*chi/sqrt(3.0_dp))
          dt3 = abs(dt3-0.5_dp*chi*(3.0_dp*chi*chi-1.0_dp))
          IF (dt1>1.E-12) WRITE (*,*) ' Error in test for dll #1 ', l1, l2, dt1
          IF (dt2>1.E-12) WRITE (*,*) ' Error in test for dll #2 ', l1, l2, dt2
          IF (dt3>1.E-12) WRITE (*,*) ' Error in test for dll #3 ', l1, l2, dt3
          emax = max(emax,dt1,dt2,dt3)
        CASE (3)
          xt1 = abs(xt1 - 0.5_dp*(5.0_dp*chi*chi-3.0_dp)*chi*chi)
          xt2 = abs(xt2 + sqrt(0.375_dp)*(5.0_dp*chi*chi-1.0_dp)*(chi*chi- &
              1.0_dp))
          xt3 = abs(xt3 - 0.25_dp*(15.0_dp*chi*chi*chi*chi-6.0_dp*chi*chi- &
              1.0_dp))
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 ', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 ', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 ', l1, l2, xt3
          emax = max(emax,xt1,xt2,xt3)
          dt1 = dpro(dll1(:,:,2),dll2(:,:,2),l1,l2)
          dt2 = dpro(dll1(:,:,2),dll2(:,:,4),l1,l2)
          dt3 = dpro(dll1(:,:,4),dll2(:,:,4),l1,l2)
          dt1 = abs(dt1-7.0_dp*(3.0_dp*chi*chi-1.0_dp))
          dt2 = abs(dt2-sqrt(1.5_dp)*(3.0_dp*chi*chi-1.0_dp))
          dt3 = abs(dt3-0.5_dp*chi*chi*(5.0_dp*chi*chi-3.0_dp))
          IF (dt1>1.E-12) WRITE (*,*) ' Error in test for dll #1 ', l1, l2, dt1
          IF (dt2>1.E-12) WRITE (*,*) ' Error in test for dll #2 ', l1, l2, dt2
          IF (dt3>1.E-12) WRITE (*,*) ' Error in test for dll #3 ', l1, l2, dt3
          emax = max(emax,dt1,dt2,dt3)
        END SELECT
      END SUBROUTINE p_test

      SUBROUTINE d_test(chi,gmu1,gmu2,dll1,dll2,l1,l2,emax)
    REAL(KIND=dp), INTENT(IN)                    :: chi, gmu1(0:6,0:6,0:3), &
                                                gmu2(0:6,0:6,0:3), &
                                                dll1(0:6,0:6,0:6), &
                                                dll2(0:6,0:6,0:6)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(KIND=dp), INTENT(OUT)                   :: emax

    INTEGER                                  :: lm, lu
    REAL(KIND=dp)                                :: c2, c4, dt1, dt2, dt3, dt4, &
                                                dt5, dt6, xt1, xt2, xt3, xt4, &
                                                xt5, xt6

        lm = max(l1,l2)
        lu = min(l1,l2)
        xt1 = dpro(gmu1(:,:,0),gmu2(:,:,0),l1,l2)
        xt2 = dpro(gmu1(:,:,0),gmu2(:,:,1),l1,l2)
        xt3 = dpro(gmu1(:,:,0),gmu2(:,:,2),l1,l2)
        xt4 = dpro(gmu1(:,:,1),gmu2(:,:,1),l1,l2)
        xt5 = dpro(gmu1(:,:,1),gmu2(:,:,2),l1,l2)
        xt6 = dpro(gmu1(:,:,2),gmu2(:,:,2),l1,l2)
        SELECT CASE (lm)
        CASE (2)
          c2 = chi*chi
          c4 = c2*c2
          xt1 = abs(xt1 - 0.25_dp*(3.0_dp*c2-1.0_dp)**2)
          xt2 = abs(xt2 + 3.0_dp*c2*(c2-1.0_dp))
          xt3 = abs(xt3 - 0.75_dp*(c2-1.0_dp)**2)
          xt4 = abs(xt4 - (4.0_dp*c4-3.0_dp*c2+1.0_dp))
          xt5 = abs(xt5 - (-c4+1.0_dp))
          xt6 = abs(xt6 - 0.25_dp*(c4+6.0_dp*c2+1.0_dp))
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 ', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 ', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 ', l1, l2, xt3
          IF (xt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 ', l1, l2, xt4
          IF (xt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 ', l1, l2, xt5
          IF (xt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 ', l1, l2, xt6
          emax = max(emax,xt1,xt2,xt3,xt4,xt5,xt6)
          dt1 = dpro(dll1(:,:,0),dll2(:,:,0),l1,l2)
          dt2 = dpro(dll1(:,:,0),dll2(:,:,2),l1,l2)
          dt3 = dpro(dll1(:,:,0),dll2(:,:,4),l1,l2)
          dt4 = dpro(dll1(:,:,2),dll2(:,:,2),l1,l2)
          dt5 = dpro(dll1(:,:,2),dll2(:,:,4),l1,l2)
          dt6 = dpro(dll1(:,:,4),dll2(:,:,4),l1,l2)
          dt1 = abs(dt1-5.0_dp)
          dt2 = abs(dt2+5.0_dp/3.0_dp)
          dt3 = abs(dt3-1.0_dp)
          dt4 = abs(dt4-1.0_dp/9.0_dp*(21.0_dp*c2-2.0_dp))
          dt5 = abs(dt5-1.0_dp/3.0_dp*(3.0_dp*c2-2.0_dp))
          dt6 = abs(dt6-1.0_dp/4.0_dp*(3.0_dp*c2-1.0_dp)**2)
          IF (dt1>1.E-12) WRITE (*,*) ' Error in test for dll #1 ', l1, l2, dt1
          IF (dt2>1.E-12) WRITE (*,*) ' Error in test for dll #2 ', l1, l2, dt2
          IF (dt3>1.E-12) WRITE (*,*) ' Error in test for dll #3 ', l1, l2, dt3
          IF (dt4>1.E-12) WRITE (*,*) ' Error in test for dll #4 ', l1, l2, dt4
          IF (dt5>1.E-12) WRITE (*,*) ' Error in test for dll #5 ', l1, l2, dt5
          IF (dt6>1.E-12) WRITE (*,*) ' Error in test for dll #6 ', l1, l2, dt6
          emax = max(emax,dt1,dt2,dt3,dt4,dt5,dt6)
        CASE (3)
          c2 = chi*chi
          c4 = c2*c2
          xt1 = abs(xt1 - 0.25_dp*(5.0_dp*c2-3.0_dp)*(3.0_dp*c2-1.0_dp)*chi)
          xt2 = abs(xt2 + 3.0_dp/sqrt(8.0_dp)*(5.0_dp*c2-1.0_dp)*(c2-1.0_dp)*chi)
          xt3 = abs(xt3 - sqrt(45.0_dp)/4.0_dp*chi*(c2-1.0_dp)**2)
          xt4 = abs(xt4 - 0.5_dp*(15.0_dp*c4-16.0_dp*c2+5.0_dp)*chi)
          xt5 = abs(xt5 + sqrt(0.625_dp)*(3.0_dp*c2+1.0_dp)*(c2-1.0_dp)*chi)
          xt6 = abs(xt6 - 0.25_dp*(3.0_dp*c4+10.0_dp*c2-5.0_dp)*chi)
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 ', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 ', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 ', l1, l2, xt3
          IF (xt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 ', l1, l2, xt4
          IF (xt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 ', l1, l2, xt5
          IF (xt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 ', l1, l2, xt6
          emax = max(emax,xt1,xt2,xt3,xt4,xt5,xt6)
          dt1 = dpro(dll1(:,:,1),dll2(:,:,1),l1,l2)
          dt2 = dpro(dll1(:,:,1),dll2(:,:,3),l1,l2)
          dt3 = dpro(dll1(:,:,1),dll2(:,:,5),l1,l2)
          dt4 = dpro(dll1(:,:,3),dll2(:,:,3),l1,l2)
          dt5 = dpro(dll1(:,:,3),dll2(:,:,5),l1,l2)
          dt6 = dpro(dll1(:,:,5),dll2(:,:,5),l1,l2)
          dt1 = abs(dt1-17.5_dp*chi)
          dt2 = abs(dt2+3.5_dp*chi)
          dt3 = abs(dt3-3.0_dp/sqrt(2.0_dp)*chi)
          dt4 = abs(dt4-0.5_dp*(24.0_dp*c2-13.0_dp)*chi)
          dt5 = abs(dt5-1.0_dp/sqrt(2.0_dp)*chi*(4.0_dp*c2-3.0_dp))
          dt6 = abs(dt6-0.25_dp*chi*(5.0_dp*c2-3.0_dp)*(3.0_dp*c2-1.0_dp))
          IF (dt1>1.E-12) WRITE (*,*) ' Error in test for dll #1 ', l1, l2, dt1
          IF (dt2>1.E-12) WRITE (*,*) ' Error in test for dll #2 ', l1, l2, dt2
          IF (dt3>1.E-12) WRITE (*,*) ' Error in test for dll #3 ', l1, l2, dt3
          IF (dt4>1.E-12) WRITE (*,*) ' Error in test for dll #4 ', l1, l2, dt4
          IF (dt5>1.E-12) WRITE (*,*) ' Error in test for dll #5 ', l1, l2, dt5
          IF (dt6>1.E-12) WRITE (*,*) ' Error in test for dll #6 ', l1, l2, dt6
          emax = max(emax,dt1,dt2,dt3,dt4,dt5,dt6)
        END SELECT
      END SUBROUTINE d_test

      SUBROUTINE f_test(chi,gmu1,gmu2,dll1,dll2,l1,l2,emax)
    REAL(KIND=dp), INTENT(IN)                    :: chi, gmu1(0:6,0:6,0:3), &
                                                gmu2(0:6,0:6,0:3), &
                                                dll1(0:6,0:6,0:6), &
                                                dll2(0:6,0:6,0:6)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(KIND=dp), INTENT(OUT)                   :: emax

    INTEGER                                  :: lm, lu
    REAL(KIND=dp)                                :: c2, c4, c6, dt0, dt1, dt2, &
                                                dt3, dt4, dt5, dt6, dt7, dt8, &
                                                dt9, xt0, xt1, xt2, xt3, xt4, &
                                                xt5, xt6, xt7, xt8, xt9

        lm = max(l1,l2)
        lu = min(l1,l2)
        xt1 = dpro(gmu1(:,:,0),gmu2(:,:,0),l1,l2)
        xt2 = dpro(gmu1(:,:,0),gmu2(:,:,1),l1,l2)
        xt3 = dpro(gmu1(:,:,0),gmu2(:,:,2),l1,l2)
        xt4 = dpro(gmu1(:,:,0),gmu2(:,:,3),l1,l2)
        xt5 = dpro(gmu1(:,:,1),gmu2(:,:,1),l1,l2)
        xt6 = dpro(gmu1(:,:,1),gmu2(:,:,2),l1,l2)
        xt7 = dpro(gmu1(:,:,1),gmu2(:,:,3),l1,l2)
        xt8 = dpro(gmu1(:,:,2),gmu2(:,:,2),l1,l2)
        xt9 = dpro(gmu1(:,:,2),gmu2(:,:,3),l1,l2)
        xt0 = dpro(gmu1(:,:,3),gmu2(:,:,3),l1,l2)
        c2 = chi*chi
        c4 = c2*c2
        c6 = c4*c2
        xt1 = abs(xt1 - 0.25_dp*(5.0_dp*c2-3.0_dp)**2*c2)
        xt2 = abs(xt2 + 0.375_dp*(5.0_dp*c2-1.0_dp)**2*(c2-1.0_dp))
        xt3 = abs(xt3 - 3.75_dp*c2*(c2-1.0_dp)**2)
        xt4 = abs(xt4 + 0.625_dp*(c2-1.0_dp)**3)
        xt5 = abs(xt5 - 0.0625_dp*(225.0_dp*c6-305.0_dp*c4+111.0_dp*c2+1.0_dp))
        xt6 = abs(xt6 + 0.625_dp*(9.0_dp*c4-2.0_dp*c2+1.0_dp)*(c2-1.0_dp))
        xt7 = abs(xt7 - 0.9375_dp*(c2+1.0_dp)*(c2-1.0_dp)**2)
        xt8 = abs(xt8 - 0.25_dp*(9.0_dp*c6+10.0_dp*c4-15.0_dp*c2+4.0_dp))
        xt9 = abs(xt9 + 0.375_dp*(c4+6.0_dp*c2+1.0_dp)*(c2-1.0_dp))
        xt0 = abs(xt0 - 0.0625_dp*(c4+14.0_dp*c2+1.0_dp)*(c2+1.0_dp))
        IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 ', l1, l2, xt1
        IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 ', l1, l2, xt2
        IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 ', l1, l2, xt3
        IF (xt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 ', l1, l2, xt4
        IF (xt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 ', l1, l2, xt5
        IF (xt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 ', l1, l2, xt6
        IF (xt7>1.E-12) WRITE (*,*) ' Error in test for gmu #7 ', l1, l2, xt7
        IF (xt8>1.E-12) WRITE (*,*) ' Error in test for gmu #8 ', l1, l2, xt8
        IF (xt9>1.E-12) WRITE (*,*) ' Error in test for gmu #9 ', l1, l2, xt9
        IF (xt0>1.E-12) WRITE (*,*) ' Error in test for gmu #0 ', l1, l2, xt0
        emax = max(emax,xt1,xt2,xt3,xt4,xt5,xt6,xt7,xt8,xt9,xt0)
        dt1 = dpro(dll1(:,:,0),dll2(:,:,0),l1,l2)
        dt2 = dpro(dll1(:,:,0),dll2(:,:,2),l1,l2)
        dt3 = dpro(dll1(:,:,0),dll2(:,:,4),l1,l2)
        dt4 = dpro(dll1(:,:,0),dll2(:,:,6),l1,l2)
        dt5 = dpro(dll1(:,:,2),dll2(:,:,2),l1,l2)
        dt6 = dpro(dll1(:,:,2),dll2(:,:,4),l1,l2)
        dt7 = dpro(dll1(:,:,2),dll2(:,:,6),l1,l2)
        dt8 = dpro(dll1(:,:,4),dll2(:,:,4),l1,l2)
        dt9 = dpro(dll1(:,:,4),dll2(:,:,6),l1,l2)
        dt0 = dpro(dll1(:,:,6),dll2(:,:,6),l1,l2)
        dt1 = abs(dt1-7.0_dp)
        dt2 = abs(dt2)
        dt3 = abs(dt3+1.4_dp)
        dt4 = abs(dt4-1.0_dp)
        dt5 = abs(dt5-1.68_dp*(3.0_dp*c2-1.0_dp))
        dt6 = abs(dt6-0.24_dp*(3.0_dp*c2-1.0_dp))
        dt7 = abs(dt7-0.4_dp*(3.0_dp*c2-1.0_dp))
        dt8 = abs(dt8-0.04_dp*(220.0_dp*c4-186.0_dp*c2+25.0_dp))
        dt9 = abs(dt9-0.6_dp*(5.0_dp*c2-4.0_dp)*c2)
        dt0 = abs(dt0-0.25_dp*(5.0_dp*c2-3.0_dp)**2*c2)
        IF (dt1>1.E-12) WRITE (*,*) ' Error in test for dll #1 ', l1, l2, dt1
        IF (dt2>1.E-12) WRITE (*,*) ' Error in test for dll #2 ', l1, l2, dt2
        IF (dt3>1.E-12) WRITE (*,*) ' Error in test for dll #3 ', l1, l2, dt3
        IF (dt4>1.E-12) WRITE (*,*) ' Error in test for dll #4 ', l1, l2, dt4
        IF (dt5>1.E-12) WRITE (*,*) ' Error in test for dll #5 ', l1, l2, dt5
        IF (dt6>1.E-12) WRITE (*,*) ' Error in test for dll #6 ', l1, l2, dt6
        IF (dt7>1.E-12) WRITE (*,*) ' Error in test for dll #7 ', l1, l2, dt7
        IF (dt8>1.E-12) WRITE (*,*) ' Error in test for dll #8 ', l1, l2, dt8
        IF (dt9>1.E-12) WRITE (*,*) ' Error in test for dll #9 ', l1, l2, dt9
        IF (dt0>1.E-12) WRITE (*,*) ' Error in test for dll #0 ', l1, l2, dt9
        IF (dt0>1.E-12) WRITE (*,*) ' Error in test for dll #0 ', l1, l2, dt0
        emax = max(emax,dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9,dt0)
      END SUBROUTINE f_test
!------------------------------------------------------------------------------!
      SUBROUTINE ds_test(chi,dchi1,dchi2,gmu1,gmu2,dgmu1,dgmu2,l1,l2,emax)
    REAL(KIND=dp), INTENT(IN) :: chi, dchi1(3), dchi2(3), gmu1(0:6,0:6,0:3), &
      gmu2(0:6,0:6,0:3), dgmu1(0:6,0:6,0:3,1:3), dgmu2(0:6,0:6,0:3,1:3)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(KIND=dp), INTENT(OUT)                   :: emax

    INTEGER                                  :: lm, lu
    REAL(KIND=dp)                                :: xt, yt, zt

        emax = 0.0_dp
        lm = max(l1,l2)
        lu = min(l1,l2)
        xt = dpro(dgmu1(:,:,0,1),gmu2(:,:,0),l1,l2)
        yt = dpro(dgmu1(:,:,0,2),gmu2(:,:,0),l1,l2)
        zt = dpro(dgmu1(:,:,0,3),gmu2(:,:,0),l1,l2)
        SELECT CASE (lm)
        CASE (0)
! ds-s
        CASE (1)
! ds-p
          xt = abs(xt - dchi1(1))
          yt = abs(yt - dchi1(2))
          zt = abs(zt - dchi1(3))
        CASE (2)
! ds-d
          xt = abs(xt - 3.0_dp*chi*dchi1(1))
          yt = abs(yt - 3.0_dp*chi*dchi1(2))
          zt = abs(zt - 3.0_dp*chi*dchi1(3))
        CASE (3)
! ds-f
          xt = abs(xt - 0.5_dp*(15.0_dp*chi*chi-3.0_dp)*dchi1(1))
          yt = abs(yt - 0.5_dp*(15.0_dp*chi*chi-3.0_dp)*dchi1(2))
          zt = abs(zt - 0.5_dp*(15.0_dp*chi*chi-3.0_dp)*dchi1(3))
        END SELECT
        IF (xt>1.E-12) WRITE (*,*) ' Error in test for dgmu x0 ', l1, l2, xt
        IF (yt>1.E-12) WRITE (*,*) ' Error in test for dgmu y0 ', l1, l2, yt
        IF (zt>1.E-12) WRITE (*,*) ' Error in test for dgmu z0 ', l1, l2, zt
        emax = max(xt,yt,zt)

        xt = dpro(gmu1(:,:,0),dgmu2(:,:,0,1),l1,l2)
        yt = dpro(gmu1(:,:,0),dgmu2(:,:,0,2),l1,l2)
        zt = dpro(gmu1(:,:,0),dgmu2(:,:,0,3),l1,l2)
        SELECT CASE (lm)
        CASE (0)
! s-ds
        CASE (1)
! s-dp
          xt = abs(xt - dchi2(1))
          yt = abs(yt - dchi2(2))
          zt = abs(zt - dchi2(3))
        CASE (2)
! s-dd
          xt = abs(xt - 3.0_dp*chi*dchi2(1))
          yt = abs(yt - 3.0_dp*chi*dchi2(2))
          zt = abs(zt - 3.0_dp*chi*dchi2(3))
        CASE (3)
! s-df
          xt = abs(xt - 0.5_dp*(15.0_dp*chi*chi-3.0_dp)*dchi2(1))
          yt = abs(yt - 0.5_dp*(15.0_dp*chi*chi-3.0_dp)*dchi2(2))
          zt = abs(zt - 0.5_dp*(15.0_dp*chi*chi-3.0_dp)*dchi2(3))
        END SELECT
        IF (xt>1.E-12) WRITE (*,*) ' Error in test for dgmu 0x', l1, l2, xt
        IF (yt>1.E-12) WRITE (*,*) ' Error in test for dgmu 0y', l1, l2, yt
        IF (zt>1.E-12) WRITE (*,*) ' Error in test for dgmu 0z', l1, l2, zt
        emax = max(emax,xt,yt,zt)
      END SUBROUTINE ds_test
      SUBROUTINE dp_test(chi,dchi1,dchi2,gmu1,gmu2,dgmu1,dgmu2, &
                         dll1,dll2,ddll1,ddll2,l1,l2,emax)
    REAL(KIND=dp), INTENT(IN) :: chi, dchi1(3), dchi2(3), gmu1(0:6,0:6,0:3), &
      gmu2(0:6,0:6,0:3), dgmu1(0:6,0:6,0:3,1:3), dgmu2(0:6,0:6,0:3,1:3), &
      dll1(0:6,0:6,0:6), dll2(0:6,0:6,0:6), ddll1(0:6,0:6,0:6,1:3), &
      ddll2(0:6,0:6,0:6,1:3)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(KIND=dp), INTENT(OUT)                   :: emax

    INTEGER                                  :: lm, lu
    REAL(KIND=dp)                                :: dx1, dx2, dx3, dy1, dy2, dy3, &
                                                dz1, dz2, dz3, xt1, xt2, xt3, &
                                                yt1, yt2, yt3, zt1, zt2, zt3

        emax = 0.0_dp
        lm = max(l1,l2)
        lu = min(l1,l2)
        xt1 = dpro(dgmu1(:,:,0,1),gmu2(:,:,0),l1,l2)
        xt2 = dpro(dgmu1(:,:,0,1),gmu2(:,:,1),l1,l2)
        xt3 = dpro(dgmu1(:,:,1,1),gmu2(:,:,1),l1,l2)
        yt1 = dpro(dgmu1(:,:,0,2),gmu2(:,:,0),l1,l2)
        yt2 = dpro(dgmu1(:,:,0,2),gmu2(:,:,1),l1,l2)
        yt3 = dpro(dgmu1(:,:,1,2),gmu2(:,:,1),l1,l2)
        zt1 = dpro(dgmu1(:,:,0,3),gmu2(:,:,0),l1,l2)
        zt2 = dpro(dgmu1(:,:,0,3),gmu2(:,:,1),l1,l2)
        zt3 = dpro(dgmu1(:,:,1,3),gmu2(:,:,1),l1,l2)
        SELECT CASE (lm)
        CASE (1)
! dp-p
          xt1 = abs(xt1 - 2.0_dp*chi*dchi1(1))
          xt2 = abs(xt2 + 2.0_dp*chi*dchi1(1))
          xt3 = abs(xt3 - 2.0_dp*chi*dchi1(1))
          yt1 = abs(yt1 - 2.0_dp*chi*dchi1(2))
          yt2 = abs(yt2 + 2.0_dp*chi*dchi1(2))
          yt3 = abs(yt3 - 2.0_dp*chi*dchi1(2))
          zt1 = abs(zt1 - 2.0_dp*chi*dchi1(3))
          zt2 = abs(zt2 + 2.0_dp*chi*dchi1(3))
          zt3 = abs(zt3 - 2.0_dp*chi*dchi1(3))
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 x0', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 x0', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 x0', l1, l2, xt3
          emax = max(emax,xt1,xt2,xt3)
          IF (yt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 y0', l1, l2, yt1
          IF (yt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 y0', l1, l2, yt2
          IF (yt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 y0', l1, l2, yt3
          emax = max(emax,yt1,yt2,yt3)
          IF (zt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 z0', l1, l2, zt1
          IF (zt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 z0', l1, l2, zt2
          IF (zt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 z0', l1, l2, zt3
          emax = max(emax,zt1,zt2,zt3)
          dx1 = dpro(ddll1(:,:,0,1),dll2(:,:,0),l1,l2)
          dx2 = dpro(ddll1(:,:,0,1),dll2(:,:,2),l1,l2)
          dx3 = dpro(ddll1(:,:,2,1),dll2(:,:,2),l1,l2)
          dy1 = dpro(ddll1(:,:,0,2),dll2(:,:,0),l1,l2)
          dy2 = dpro(ddll1(:,:,0,2),dll2(:,:,2),l1,l2)
          dy3 = dpro(ddll1(:,:,2,2),dll2(:,:,2),l1,l2)
          dz1 = dpro(ddll1(:,:,0,3),dll2(:,:,0),l1,l2)
          dz2 = dpro(ddll1(:,:,0,3),dll2(:,:,2),l1,l2)
          dz3 = dpro(ddll1(:,:,2,3),dll2(:,:,2),l1,l2)
          dx1 = abs(dx1)
          dx2 = abs(dx2)
          dx3 = abs(dx3-2.0_dp*chi*dchi1(1))
          dy1 = abs(dy1)
          dy2 = abs(dy2)
          dy3 = abs(dy3-2.0_dp*chi*dchi1(2))
          dz1 = abs(dz1)
          dz2 = abs(dz2)
          dz3 = abs(dz3-2.0_dp*chi*dchi1(3))
          IF (dx1>1.E-12) WRITE (*,*) ' Error in test for dll #1 x0', l1, l2, dx1
          IF (dx2>1.E-12) WRITE (*,*) ' Error in test for dll #2 x0', l1, l2, dx2
          IF (dx3>1.E-12) WRITE (*,*) ' Error in test for dll #3 x0', l1, l2, dx3
          emax = max(emax,dx1,dx2,dx3)
          IF (dy1>1.E-12) WRITE (*,*) ' Error in test for dll #1 y0', l1, l2, dy1
          IF (dy2>1.E-12) WRITE (*,*) ' Error in test for dll #2 y0', l1, l2, dy2
          IF (dy3>1.E-12) WRITE (*,*) ' Error in test for dll #3 y0', l1, l2, dy3
          emax = max(emax,dy1,dy2,dy3)
          IF (dz1>1.E-12) WRITE (*,*) ' Error in test for dll #1 z0', l1, l2, dz1
          IF (dz2>1.E-12) WRITE (*,*) ' Error in test for dll #2 z0', l1, l2, dz2
          IF (dz3>1.E-12) WRITE (*,*) ' Error in test for dll #3 z0', l1, l2, dz3
          emax = max(emax,dz1,dz2,dz3)
        CASE (2)
! dp-d
          xt1 = abs(xt1 - 0.5_dp*(9.0_dp*chi*chi-1.0_dp)*dchi1(1))
          xt2 = abs(xt2 + sqrt(3.0_dp)*(3.0_dp*chi*chi-1.0_dp)*dchi1(1))
          xt3 = abs(xt3 - 6.0_dp*chi*chi*dchi1(1))
          yt1 = abs(yt1 - 0.5_dp*(9.0_dp*chi*chi-1.0_dp)*dchi1(2))
          yt2 = abs(yt2 + sqrt(3.0_dp)*(3.0_dp*chi*chi-1.0_dp)*dchi1(2))
          yt3 = abs(yt3 - 6.0_dp*chi*chi*dchi1(2))
          zt1 = abs(zt1 - 0.5_dp*(9.0_dp*chi*chi-1.0_dp)*dchi1(3))
          zt2 = abs(zt2 + sqrt(3.0_dp)*(3.0_dp*chi*chi-1.0_dp)*dchi1(3))
          zt3 = abs(zt3 - 6.0_dp*chi*chi*dchi1(3))
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 x0', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 x0', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 x0', l1, l2, xt3
          emax = max(emax,xt1,xt2,xt3)
          IF (yt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 y0', l1, l2, yt1
          IF (yt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 y0', l1, l2, yt2
          IF (yt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 y0', l1, l2, yt3
          emax = max(emax,yt1,yt2,yt3)
          IF (zt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 z0', l1, l2, zt1
          IF (zt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 z0', l1, l2, zt2
          IF (zt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 z0', l1, l2, zt3
          emax = max(emax,zt1,zt2,zt3)
          dx1 = dpro(ddll1(:,:,1,1),dll2(:,:,1),l1,l2)
          dx2 = dpro(ddll1(:,:,1,1),dll2(:,:,3),l1,l2)
          dx3 = dpro(ddll1(:,:,3,1),dll2(:,:,3),l1,l2)
          dy1 = dpro(ddll1(:,:,1,2),dll2(:,:,1),l1,l2)
          dy2 = dpro(ddll1(:,:,1,2),dll2(:,:,3),l1,l2)
          dy3 = dpro(ddll1(:,:,3,2),dll2(:,:,3),l1,l2)
          dz1 = dpro(ddll1(:,:,1,3),dll2(:,:,1),l1,l2)
          dz2 = dpro(ddll1(:,:,1,3),dll2(:,:,3),l1,l2)
          dz3 = dpro(ddll1(:,:,3,3),dll2(:,:,3),l1,l2)
          dx1 = abs(dx1-10.0_dp/3.0_dp*dchi1(1))
          dx2 = abs(dx2-2.0_dp*dchi1(1)/sqrt(3.0_dp))
          dx3 = abs(dx3-0.5_dp*dchi1(1)*(9.0_dp*chi*chi-1.0_dp))
          dy1 = abs(dy1-10.0_dp/3.0_dp*dchi1(2))
          dy2 = abs(dy2-2.0_dp*dchi1(2)/sqrt(3.0_dp))
          dy3 = abs(dy3-0.5_dp*dchi1(2)*(9.0_dp*chi*chi-1.0_dp))
          dz1 = abs(dz1-10.0_dp/3.0_dp*dchi1(3))
          dz2 = abs(dz2-2.0_dp*dchi1(3)/sqrt(3.0_dp))
          dz3 = abs(dz3-0.5_dp*dchi1(3)*(9.0_dp*chi*chi-1.0_dp))
          IF (dx1>1.E-12) WRITE (*,*) ' Error in test for dll #1 x0', l1, l2, dx1
          IF (dx2>1.E-12) WRITE (*,*) ' Error in test for dll #2 x0', l1, l2, dx2
          IF (dx3>1.E-12) WRITE (*,*) ' Error in test for dll #3 x0', l1, l2, dx3
          emax = max(emax,dx1,dx2,dx3)
          IF (dy1>1.E-12) WRITE (*,*) ' Error in test for dll #1 y0', l1, l2, dy1
          IF (dy2>1.E-12) WRITE (*,*) ' Error in test for dll #2 y0', l1, l2, dy2
          IF (dy3>1.E-12) WRITE (*,*) ' Error in test for dll #3 y0', l1, l2, dy3
          emax = max(emax,dy1,dy2,dy3)
          IF (dz1>1.E-12) WRITE (*,*) ' Error in test for dll #1 z0', l1, l2, dz1
          IF (dz2>1.E-12) WRITE (*,*) ' Error in test for dll #2 z0', l1, l2, dz2
          IF (dz3>1.E-12) WRITE (*,*) ' Error in test for dll #3 z0', l1, l2, dz3
          emax = max(emax,dz1,dz2,dz3)
        CASE (3)
! dp-f
          xt1 = abs(xt1 - (10.0_dp*chi*chi-3.0_dp)*chi*dchi1(1))
          xt2 = abs(xt2 + sqrt(0.375_dp)*(20.0_dp*chi*chi-12.0_dp)*chi*dchi1(1))
          xt3 = abs(xt3 - (15.0_dp*chi*chi-3.0_dp)*chi*dchi1(1))
          yt1 = abs(yt1 - (10.0_dp*chi*chi-3.0_dp)*chi*dchi1(2))
          yt2 = abs(yt2 + sqrt(0.375_dp)*(20.0_dp*chi*chi-12.0_dp)*chi*dchi1(2))
          yt3 = abs(yt3 - (15.0_dp*chi*chi-3.0_dp)*chi*dchi1(2))
          zt1 = abs(zt1 - (10.0_dp*chi*chi-3.0_dp)*chi*dchi1(3))
          zt2 = abs(zt2 + sqrt(0.375_dp)*(20.0_dp*chi*chi-12.0_dp)*chi*dchi1(3))
          zt3 = abs(zt3 - (15.0_dp*chi*chi-3.0_dp)*chi*dchi1(3))
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 x0', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 x0', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 x0', l1, l2, xt3
          emax = max(emax,xt1,xt2,xt3)
          IF (yt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 y0', l1, l2, yt1
          IF (yt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 y0', l1, l2, yt2
          IF (yt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 y0', l1, l2, yt3
          emax = max(emax,yt1,yt2,yt3)
          IF (zt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 z0', l1, l2, zt1
          IF (zt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 z0', l1, l2, zt2
          IF (zt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 z0', l1, l2, zt3
          emax = max(emax,zt1,zt2,zt3)
          dx1 = dpro(ddll1(:,:,2,1),dll2(:,:,2),l1,l2)
          dx2 = dpro(ddll1(:,:,2,1),dll2(:,:,4),l1,l2)
          dx3 = dpro(ddll1(:,:,4,1),dll2(:,:,4),l1,l2)
          dy1 = dpro(ddll1(:,:,2,2),dll2(:,:,2),l1,l2)
          dy2 = dpro(ddll1(:,:,2,2),dll2(:,:,4),l1,l2)
          dy3 = dpro(ddll1(:,:,4,2),dll2(:,:,4),l1,l2)
          dz1 = dpro(ddll1(:,:,2,3),dll2(:,:,2),l1,l2)
          dz2 = dpro(ddll1(:,:,2,3),dll2(:,:,4),l1,l2)
          dz3 = dpro(ddll1(:,:,4,3),dll2(:,:,4),l1,l2)
          dx1 = abs(dx1-42.0_dp*chi*dchi1(1))
          dx2 = abs(dx2-sqrt(1.5_dp)*6.0_dp*chi*dchi1(1))
          dx3 = abs(dx3-chi*dchi1(1)*(10.0_dp*chi*chi-3.0_dp))
          dy1 = abs(dy1-42.0_dp*chi*dchi1(2))
          dy2 = abs(dy2-sqrt(1.5_dp)*6.0_dp*chi*dchi1(2))
          dy3 = abs(dy3-chi*dchi1(2)*(10.0_dp*chi*chi-3.0_dp))
          dz1 = abs(dz1-42.0_dp*chi*dchi1(3))
          dz2 = abs(dz2-sqrt(1.5_dp)*6.0_dp*chi*dchi1(3))
          dz3 = abs(dz3-chi*dchi1(3)*(10.0_dp*chi*chi-3.0_dp))
          IF (dx1>1.E-12) WRITE (*,*) ' Error in test for dll #1 x0', l1, l2, dx1
          IF (dx2>1.E-12) WRITE (*,*) ' Error in test for dll #2 x0', l1, l2, dx2
          IF (dx3>1.E-12) WRITE (*,*) ' Error in test for dll #3 x0', l1, l2, dx3
          emax = max(emax,dx1,dx2,dx3)
          IF (dy1>1.E-12) WRITE (*,*) ' Error in test for dll #1 y0', l1, l2, dy1
          IF (dy2>1.E-12) WRITE (*,*) ' Error in test for dll #2 y0', l1, l2, dy2
          IF (dy3>1.E-12) WRITE (*,*) ' Error in test for dll #3 y0', l1, l2, dy3
          emax = max(emax,dy1,dy2,dy3)
          IF (dz1>1.E-12) WRITE (*,*) ' Error in test for dll #1 z0', l1, l2, dz1
          IF (dz2>1.E-12) WRITE (*,*) ' Error in test for dll #2 z0', l1, l2, dz2
          IF (dz3>1.E-12) WRITE (*,*) ' Error in test for dll #3 z0', l1, l2, dz3
          emax = max(emax,dz1,dz2,dz3)
        END SELECT
        xt1 = dpro(gmu1(:,:,0),dgmu2(:,:,0,1),l1,l2)
        xt2 = dpro(gmu1(:,:,0),dgmu2(:,:,1,1),l1,l2)
        xt3 = dpro(gmu1(:,:,1),dgmu2(:,:,1,1),l1,l2)
        yt1 = dpro(gmu1(:,:,0),dgmu2(:,:,0,2),l1,l2)
        yt2 = dpro(gmu1(:,:,0),dgmu2(:,:,1,2),l1,l2)
        yt3 = dpro(gmu1(:,:,1),dgmu2(:,:,1,2),l1,l2)
        zt1 = dpro(gmu1(:,:,0),dgmu2(:,:,0,3),l1,l2)
        zt2 = dpro(gmu1(:,:,0),dgmu2(:,:,1,3),l1,l2)
        zt3 = dpro(gmu1(:,:,1),dgmu2(:,:,1,3),l1,l2)
        SELECT CASE (lm)
        CASE (1)
! p-dp
          xt1 = abs(xt1 - 2.0_dp*chi*dchi2(1))
          xt2 = abs(xt2 + 2.0_dp*chi*dchi2(1))
          xt3 = abs(xt3 - 2.0_dp*chi*dchi2(1))
          yt1 = abs(yt1 - 2.0_dp*chi*dchi2(2))
          yt2 = abs(yt2 + 2.0_dp*chi*dchi2(2))
          yt3 = abs(yt3 - 2.0_dp*chi*dchi2(2))
          zt1 = abs(zt1 - 2.0_dp*chi*dchi2(3))
          zt2 = abs(zt2 + 2.0_dp*chi*dchi2(3))
          zt3 = abs(zt3 - 2.0_dp*chi*dchi2(3))
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0x', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0x', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0x', l1, l2, xt3
          emax = max(emax,xt1,xt2,xt3)
          IF (yt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0y', l1, l2, yt1
          IF (yt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0y', l1, l2, yt2
          IF (yt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0y', l1, l2, yt3
          emax = max(emax,yt1,yt2,yt3)
          IF (zt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0z', l1, l2, zt1
          IF (zt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0z', l1, l2, zt2
          IF (zt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0z', l1, l2, zt3
          emax = max(emax,zt1,zt2,zt3)
          dx1 = dpro(dll1(:,:,0),ddll2(:,:,0,1),l1,l2)
          dx2 = dpro(dll1(:,:,0),ddll2(:,:,2,1),l1,l2)
          dx3 = dpro(dll1(:,:,2),ddll2(:,:,2,1),l1,l2)
          dy1 = dpro(dll1(:,:,0),ddll2(:,:,0,2),l1,l2)
          dy2 = dpro(dll1(:,:,0),ddll2(:,:,2,2),l1,l2)
          dy3 = dpro(dll1(:,:,2),ddll2(:,:,2,2),l1,l2)
          dz1 = dpro(dll1(:,:,0),ddll2(:,:,0,3),l1,l2)
          dz2 = dpro(dll1(:,:,0),ddll2(:,:,2,3),l1,l2)
          dz3 = dpro(dll1(:,:,2),ddll2(:,:,2,3),l1,l2)
          dx1 = abs(dx1)
          dx2 = abs(dx2)
          dx3 = abs(dx3-2.0_dp*chi*dchi2(1))
          dy1 = abs(dy1)
          dy2 = abs(dy2)
          dy3 = abs(dy3-2.0_dp*chi*dchi2(2))
          dz1 = abs(dz1)
          dz2 = abs(dz2)
          dz3 = abs(dz3-2.0_dp*chi*dchi2(3))
          IF (dx1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0x', l1, l2, dx1
          IF (dx2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0x', l1, l2, dx2
          IF (dx3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0x', l1, l2, dx3
          emax = max(emax,dx1,dx2,dx3)
          IF (dy1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0y', l1, l2, dy1
          IF (dy2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0y', l1, l2, dy2
          IF (dy3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0y', l1, l2, dy3
          emax = max(emax,dy1,dy2,dy3)
          IF (dz1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0z', l1, l2, dz1
          IF (dz2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0z', l1, l2, dz2
          IF (dz3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0z', l1, l2, dz3
          emax = max(emax,dz1,dz2,dz3)
        CASE (2)
! p-dd
          xt1 = abs(xt1 - 0.5_dp*(9.0_dp*chi*chi-1.0_dp)*dchi2(1))
          xt2 = abs(xt2 + sqrt(3.0_dp)*(3.0_dp*chi*chi-1.0_dp)*dchi2(1))
          xt3 = abs(xt3 - 6.0_dp*chi*chi*dchi2(1))
          yt1 = abs(yt1 - 0.5_dp*(9.0_dp*chi*chi-1.0_dp)*dchi2(2) )
          yt2 = abs(yt2 + sqrt(3.0_dp)*(3.0_dp*chi*chi-1.0_dp)*dchi2(2))
          yt3 = abs(yt3 - 6.0_dp*chi*chi*dchi2(2))
          zt1 = abs(zt1 - 0.5_dp*(9.0_dp*chi*chi-1.0_dp)*dchi2(3) )
          zt2 = abs(zt2 + sqrt(3.0_dp)*(3.0_dp*chi*chi-1.0_dp)*dchi2(3))
          zt3 = abs(zt3 - 6.0_dp*chi*chi*dchi2(3))
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0x', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0x', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0x', l1, l2, xt3
          emax = max(emax,xt1,xt2,xt3)
          IF (yt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0y', l1, l2, yt1
          IF (yt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0y', l1, l2, yt2
          IF (yt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0y', l1, l2, yt3
          emax = max(emax,yt1,yt2,yt3)
          IF (zt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0z', l1, l2, zt1
          IF (zt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0z', l1, l2, zt2
          IF (zt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0z', l1, l2, zt3
          emax = max(emax,zt1,zt2,zt3)
          dx1 = dpro(dll1(:,:,1),ddll2(:,:,1,1),l1,l2)
          dx2 = dpro(dll1(:,:,1),ddll2(:,:,3,1),l1,l2)
          dx3 = dpro(dll1(:,:,3),ddll2(:,:,3,1),l1,l2)
          dy1 = dpro(dll1(:,:,1),ddll2(:,:,1,2),l1,l2)
          dy2 = dpro(dll1(:,:,1),ddll2(:,:,3,2),l1,l2)
          dy3 = dpro(dll1(:,:,3),ddll2(:,:,3,2),l1,l2)
          dz1 = dpro(dll1(:,:,1),ddll2(:,:,1,3),l1,l2)
          dz2 = dpro(dll1(:,:,1),ddll2(:,:,3,3),l1,l2)
          dz3 = dpro(dll1(:,:,3),ddll2(:,:,3,3),l1,l2)
          dx1 = abs(dx1-10.0_dp/3.0_dp*dchi2(1))
          dx2 = abs(dx2-2.0_dp*dchi2(1)/sqrt(3.0_dp))
          dx3 = abs(dx3-0.5_dp*dchi2(1)*(9.0_dp*chi*chi-1.0_dp))
          dy1 = abs(dy1-10.0_dp/3.0_dp*dchi2(2))
          dy2 = abs(dy2-2.0_dp*dchi2(2)/sqrt(3.0_dp))
          dy3 = abs(dy3-0.5_dp*dchi2(2)*(9.0_dp*chi*chi-1.0_dp))
          dz1 = abs(dz1-10.0_dp/3.0_dp*dchi2(3))
          dz2 = abs(dz2-2.0_dp*dchi2(3)/sqrt(3.0_dp))
          dz3 = abs(dz3-0.5_dp*dchi2(3)*(9.0_dp*chi*chi-1.0_dp))
          IF (dx1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0x', l1, l2, dx1
          IF (dx2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0x', l1, l2, dx2
          IF (dx3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0x', l1, l2, dx3
          emax = max(emax,dx1,dx2,dx3)
          IF (dy1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0y', l1, l2, dy1
          IF (dy2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0y', l1, l2, dy2
          IF (dy3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0y', l1, l2, dy3
          emax = max(emax,dy1,dy2,dy3)
          IF (dz1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0z', l1, l2, dz1
          IF (dz2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0z', l1, l2, dz2
          IF (dz3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0z', l1, l2, dz3
          emax = max(emax,dz1,dz2,dz3)
        CASE (3)
! p-df
          xt1 = abs(xt1 - (10.0_dp*chi*chi-3.0_dp)*chi*dchi2(1))
          xt2 = abs(xt2 + sqrt(0.375_dp)*(20.0_dp*chi*chi-12.0_dp)*chi*dchi2(1))
          xt3 = abs(xt3 - (15.0_dp*chi*chi-3.0_dp)*chi*dchi2(1))
          yt1 = abs(yt1 - (10.0_dp*chi*chi-3.0_dp)*chi*dchi2(2))
          yt2 = abs(yt2 + sqrt(0.375_dp)*(20.0_dp*chi*chi-12.0_dp)*chi*dchi2(2))
          yt3 = abs(yt3 - (15.0_dp*chi*chi-3.0_dp)*chi*dchi2(2))
          zt1 = abs(zt1 - (10.0_dp*chi*chi-3.0_dp)*chi*dchi2(3))
          zt2 = abs(zt2 + sqrt(0.375_dp)*(20.0_dp*chi*chi-12.0_dp)*chi*dchi2(3))
          zt3 = abs(zt3 - (15.0_dp*chi*chi-3.0_dp)*chi*dchi2(3))
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0x', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0x', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0x', l1, l2, xt3
          emax = max(emax,xt1,xt2,xt3)
          IF (yt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0y', l1, l2, yt1
          IF (yt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0y', l1, l2, yt2
          IF (yt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0y', l1, l2, yt3
          emax = max(emax,yt1,yt2,yt3)
          IF (zt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0z', l1, l2, zt1
          IF (zt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0z', l1, l2, zt2
          IF (zt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0z', l1, l2, zt3
          emax = max(emax,zt1,zt2,zt3)
          dx1 = dpro(dll1(:,:,2),ddll2(:,:,2,1),l1,l2)
          dx2 = dpro(dll1(:,:,2),ddll2(:,:,4,1),l1,l2)
          dx3 = dpro(dll1(:,:,4),ddll2(:,:,4,1),l1,l2)
          dy1 = dpro(dll1(:,:,2),ddll2(:,:,2,2),l1,l2)
          dy2 = dpro(dll1(:,:,2),ddll2(:,:,4,2),l1,l2)
          dy3 = dpro(dll1(:,:,4),ddll2(:,:,4,2),l1,l2)
          dz1 = dpro(dll1(:,:,2),ddll2(:,:,2,3),l1,l2)
          dz2 = dpro(dll1(:,:,2),ddll2(:,:,4,3),l1,l2)
          dz3 = dpro(dll1(:,:,4),ddll2(:,:,4,3),l1,l2)
          dx1 = abs(dx1-42.0_dp*chi*dchi2(1))
          dx2 = abs(dx2-sqrt(1.5_dp)*6.0_dp*chi*dchi2(1))
          dx3 = abs(dx3-chi*dchi2(1)*(10.0_dp*chi*chi-3.0_dp))
          dy1 = abs(dy1-42.0_dp*chi*dchi2(2))
          dy2 = abs(dy2-sqrt(1.5_dp)*6.0_dp*chi*dchi2(2))
          dy3 = abs(dy3-chi*dchi2(2)*(10.0_dp*chi*chi-3.0_dp))
          dz1 = abs(dz1-42.0_dp*chi*dchi2(3))
          dz2 = abs(dz2-sqrt(1.5_dp)*6.0_dp*chi*dchi2(3))
          dz3 = abs(dz3-chi*dchi2(3)*(10.0_dp*chi*chi-3.0_dp))
          IF (dx1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0x', l1, l2, dx1
          IF (dx2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0x', l1, l2, dx2
          IF (dx3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0x', l1, l2, dx3
          emax = max(emax,dx1,dx2,dx3)
          IF (dy1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0y', l1, l2, dy1
          IF (dy2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0y', l1, l2, dy2
          IF (dy3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0y', l1, l2, dy3
          emax = max(emax,dy1,dy2,dy3)
          IF (dz1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0z', l1, l2, dz1
          IF (dz2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0z', l1, l2, dz2
          IF (dz3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0z', l1, l2, dz3
          emax = max(emax,dz1,dz2,dz3)
        END SELECT
      END SUBROUTINE dp_test
      SUBROUTINE dd_test(chi,dchi1,dchi2,gmu1,gmu2,dgmu1,dgmu2, &
                         dll1,dll2,ddll1,ddll2,l1,l2,emax)
    REAL(KIND=dp), INTENT(IN) :: chi, dchi1(3), dchi2(3), gmu1(0:6,0:6,0:3), &
      gmu2(0:6,0:6,0:3), dgmu1(0:6,0:6,0:3,1:3), dgmu2(0:6,0:6,0:3,1:3), &
      dll1(0:6,0:6,0:6), dll2(0:6,0:6,0:6), ddll1(0:6,0:6,0:6,1:3), &
      ddll2(0:6,0:6,0:6,1:3)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(KIND=dp), INTENT(OUT)                   :: emax

    INTEGER                                  :: lm, lu
    REAL(KIND=dp) :: c2, c4, dx1, dx2, dx3, dx4, dx5, dx6, dy1, dy2, dy3, dy4, &
      dy5, dy6, dz1, dz2, dz3, dz4, dz5, dz6, xt1, xt2, xt3, xt4, xt5, xt6, &
      yt1, yt2, yt3, yt4, yt5, yt6, zt1, zt2, zt3, zt4, zt5, zt6

        emax = 0.0_dp
        lm = max(l1,l2)
        lu = min(l1,l2)

        xt1 = dpro(dgmu1(:,:,0,1),gmu2(:,:,0),l1,l2)
        xt2 = dpro(dgmu1(:,:,0,1),gmu2(:,:,1),l1,l2)
        xt3 = dpro(dgmu1(:,:,0,1),gmu2(:,:,2),l1,l2)
        xt4 = dpro(dgmu1(:,:,1,1),gmu2(:,:,1),l1,l2)
        xt5 = dpro(dgmu1(:,:,1,1),gmu2(:,:,2),l1,l2)
        xt6 = dpro(dgmu1(:,:,2,1),gmu2(:,:,2),l1,l2)
        yt1 = dpro(dgmu1(:,:,0,2),gmu2(:,:,0),l1,l2)
        yt2 = dpro(dgmu1(:,:,0,2),gmu2(:,:,1),l1,l2)
        yt3 = dpro(dgmu1(:,:,0,2),gmu2(:,:,2),l1,l2)
        yt4 = dpro(dgmu1(:,:,1,2),gmu2(:,:,1),l1,l2)
        yt5 = dpro(dgmu1(:,:,1,2),gmu2(:,:,2),l1,l2)
        yt6 = dpro(dgmu1(:,:,2,2),gmu2(:,:,2),l1,l2)
        zt1 = dpro(dgmu1(:,:,0,3),gmu2(:,:,0),l1,l2)
        zt2 = dpro(dgmu1(:,:,0,3),gmu2(:,:,1),l1,l2)
        zt3 = dpro(dgmu1(:,:,0,3),gmu2(:,:,2),l1,l2)
        zt4 = dpro(dgmu1(:,:,1,3),gmu2(:,:,1),l1,l2)
        zt5 = dpro(dgmu1(:,:,1,3),gmu2(:,:,2),l1,l2)
        zt6 = dpro(dgmu1(:,:,2,3),gmu2(:,:,2),l1,l2)
        SELECT CASE (lm)
        CASE (2)
! dd-d
          xt1 = abs(xt1 - 3.0_dp*chi*(3.0_dp*chi*chi-1.0_dp)*dchi1(1))
          xt2 = abs(xt2 + 6.0_dp*chi*(2.0_dp*chi*chi-1.0_dp)*dchi1(1))
          xt3 = abs(xt3 - 3.0_dp*chi*(chi*chi-1.0_dp)*dchi1(1))
          xt4 = abs(xt4 - 2.0_dp*chi*(8.0_dp*chi*chi-3.0_dp)*dchi1(1))
          xt5 = abs(xt5 + 4.0_dp*chi*chi*chi*dchi1(1))
          xt6 = abs(xt6 - chi*(chi*chi+3.0_dp)*dchi1(1))
          yt1 = abs(yt1 - 3.0_dp*chi*(3.0_dp*chi*chi-1.0_dp)*dchi1(2))
          yt2 = abs(yt2 + 6.0_dp*chi*(2.0_dp*chi*chi-1.0_dp)*dchi1(2))
          yt3 = abs(yt3 - 3.0_dp*chi*(chi*chi-1.0_dp)*dchi1(2))
          yt4 = abs(yt4 - 2.0_dp*chi*(8.0_dp*chi*chi-3.0_dp)*dchi1(2))
          yt5 = abs(yt5 + 4.0_dp*chi*chi*chi*dchi1(2))
          yt6 = abs(yt6 - chi*(chi*chi+3.0_dp)*dchi1(2))
          zt1 = abs(zt1 - 3.0_dp*chi*(3.0_dp*chi*chi-1.0_dp)*dchi1(3))
          zt2 = abs(zt2 + 6.0_dp*chi*(2.0_dp*chi*chi-1.0_dp)*dchi1(3))
          zt3 = abs(zt3 - 3.0_dp*chi*(chi*chi-1.0_dp)*dchi1(3))
          zt4 = abs(zt4 - 2.0_dp*chi*(8.0_dp*chi*chi-3.0_dp)*dchi1(3))
          zt5 = abs(zt5 + 4.0_dp*chi*chi*chi*dchi1(3))
          zt6 = abs(zt6 - chi*(chi*chi+3.0_dp)*dchi1(3))
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 x0', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 x0', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 x0', l1, l2, xt3
          IF (xt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 x0', l1, l2, xt4
          IF (xt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 x0', l1, l2, xt5
          IF (xt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 x0', l1, l2, xt6
          emax = max(emax,xt1,xt2,xt3,xt4,xt5,xt6)
          IF (yt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 y0', l1, l2, yt1
          IF (yt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 y0', l1, l2, yt2
          IF (yt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 y0', l1, l2, yt3
          IF (yt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 y0', l1, l2, yt4
          IF (yt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 y0', l1, l2, yt5
          IF (yt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 y0', l1, l2, yt6
          emax = max(emax,yt1,yt2,yt3,yt4,yt5,yt6)
          IF (zt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 z0', l1, l2, zt1
          IF (zt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 z0', l1, l2, zt2
          IF (zt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 z0', l1, l2, zt3
          IF (zt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 z0', l1, l2, zt4
          IF (zt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 z0', l1, l2, zt5
          IF (zt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 z0', l1, l2, zt6
          emax = max(emax,zt1,zt2,zt3,zt4,zt5,zt6)
          dx1 = dpro(ddll1(:,:,0,1),dll2(:,:,0),l1,l2)
          dx2 = dpro(ddll1(:,:,0,1),dll2(:,:,2),l1,l2)
          dx3 = dpro(ddll1(:,:,0,1),dll2(:,:,4),l1,l2)
          dx4 = dpro(ddll1(:,:,2,1),dll2(:,:,2),l1,l2)
          dx5 = dpro(ddll1(:,:,2,1),dll2(:,:,4),l1,l2)
          dx6 = dpro(ddll1(:,:,4,1),dll2(:,:,4),l1,l2)
          dy1 = dpro(ddll1(:,:,0,2),dll2(:,:,0),l1,l2)
          dy2 = dpro(ddll1(:,:,0,2),dll2(:,:,2),l1,l2)
          dy3 = dpro(ddll1(:,:,0,2),dll2(:,:,4),l1,l2)
          dy4 = dpro(ddll1(:,:,2,2),dll2(:,:,2),l1,l2)
          dy5 = dpro(ddll1(:,:,2,2),dll2(:,:,4),l1,l2)
          dy6 = dpro(ddll1(:,:,4,2),dll2(:,:,4),l1,l2)
          dz1 = dpro(ddll1(:,:,0,3),dll2(:,:,0),l1,l2)
          dz2 = dpro(ddll1(:,:,0,3),dll2(:,:,2),l1,l2)
          dz3 = dpro(ddll1(:,:,0,3),dll2(:,:,4),l1,l2)
          dz4 = dpro(ddll1(:,:,2,3),dll2(:,:,2),l1,l2)
          dz5 = dpro(ddll1(:,:,2,3),dll2(:,:,4),l1,l2)
          dz6 = dpro(ddll1(:,:,4,3),dll2(:,:,4),l1,l2)
          dx1 = abs(dx1)
          dx2 = abs(dx2)
          dx3 = abs(dx3)
          dx4 = abs(dx4-14.0_dp/3.0_dp*chi*dchi1(1))
          dx5 = abs(dx5-2.0_dp*chi*dchi1(1))
          dx6 = abs(dx6-3.0_dp*(3.0_dp*chi*chi-1.0_dp)*chi*dchi1(1))
          dy1 = abs(dy1)
          dy2 = abs(dy2)
          dy3 = abs(dy3)
          dy4 = abs(dy4-14.0_dp/3.0_dp*chi*dchi1(2))
          dy5 = abs(dy5-2.0_dp*chi*dchi1(2))
          dy6 = abs(dy6-3.0_dp*(3.0_dp*chi*chi-1.0_dp)*chi*dchi1(2))
          dz1 = abs(dz1)
          dz2 = abs(dz2)
          dz3 = abs(dz3)
          dz4 = abs(dz4-14.0_dp/3.0_dp*chi*dchi1(3))
          dz5 = abs(dz5-2.0_dp*chi*dchi1(3))
          dz6 = abs(dz6-3.0_dp*(3.0_dp*chi*chi-1.0_dp)*chi*dchi1(3))
          IF (dx1>1.E-12) WRITE (*,*) ' Error in test for dll #1 x0', l1, l2, dx1
          IF (dx2>1.E-12) WRITE (*,*) ' Error in test for dll #2 x0', l1, l2, dx2
          IF (dx3>1.E-12) WRITE (*,*) ' Error in test for dll #3 x0', l1, l2, dx3
          IF (dx4>1.E-12) WRITE (*,*) ' Error in test for dll #4 x0', l1, l2, dx4
          IF (dx5>1.E-12) WRITE (*,*) ' Error in test for dll #5 x0', l1, l2, dx5
          IF (dx6>1.E-12) WRITE (*,*) ' Error in test for dll #6 x0', l1, l2, dx6
          emax = max(emax,dx1,dx2,dx3,dx4,dx5,dx6)
          IF (dy1>1.E-12) WRITE (*,*) ' Error in test for dll #1 y0', l1, l2, dy1
          IF (dy2>1.E-12) WRITE (*,*) ' Error in test for dll #2 y0', l1, l2, dy2
          IF (dy3>1.E-12) WRITE (*,*) ' Error in test for dll #3 y0', l1, l2, dy3
          IF (dy4>1.E-12) WRITE (*,*) ' Error in test for dll #4 y0', l1, l2, dy4
          IF (dy5>1.E-12) WRITE (*,*) ' Error in test for dll #5 y0', l1, l2, dy5
          IF (dy6>1.E-12) WRITE (*,*) ' Error in test for dll #6 y0', l1, l2, dy6
          emax = max(emax,dy1,dy2,dy3,dy4,dy5,dy6)
          IF (dz1>1.E-12) WRITE (*,*) ' Error in test for dll #1 z0', l1, l2, dz1
          IF (dz2>1.E-12) WRITE (*,*) ' Error in test for dll #2 z0', l1, l2, dz2
          IF (dz3>1.E-12) WRITE (*,*) ' Error in test for dll #3 z0', l1, l2, dz3
          IF (dz4>1.E-12) WRITE (*,*) ' Error in test for dll #4 z0', l1, l2, dz4
          IF (dz5>1.E-12) WRITE (*,*) ' Error in test for dll #5 z0', l1, l2, dz5
          IF (dz6>1.E-12) WRITE (*,*) ' Error in test for dll #6 z0', l1, l2, dz6
          emax = max(emax,dx1,dx2,dx3,dx4,dx5,dx6)
        CASE (3)
! dd-f
          c2 = chi*chi
          c4 = c2*c2
          xt1 = abs(xt1 - 0.25_dp*(75.0_dp*c4-42.0_dp*c2+3.0_dp)*dchi1(1))
          xt2 = abs(xt2 + 3.0_dp/sqrt(8.0_dp)*(25.0_dp*c4-18.0_dp*c2+1.0_dp)*dchi1(1))
          xt3 = abs(xt3 - sqrt(45.0_dp)/4.0_dp*(5.0_dp*c4-6.0_dp*c2+1.0_dp)*dchi1(1))
          xt4 = abs(xt4 - 0.5_dp*(75.0_dp*c4-48.0_dp*c2+5.0_dp)*dchi1(1))
          xt5 = abs(xt5 + sqrt(0.625_dp)*(15.0_dp*c4-6.0_dp*c2-1.0_dp)*dchi1(1))
          xt6 = abs(xt6 - 0.25_dp*(15.0_dp*c4+30.0_dp*c2-5.0_dp)*dchi1(1))
          yt1 = abs(yt1 - 0.25_dp*(75.0_dp*c4-42.0_dp*c2+3.0_dp)*dchi1(2))
          yt2 = abs(yt2 + 3.0_dp/sqrt(8.0_dp)*(25.0_dp*c4-18.0_dp*c2+1.0_dp)*dchi1(2))
          yt3 = abs(yt3 - sqrt(45.0_dp)/4.0_dp*(5.0_dp*c4-6.0_dp*c2+1.0_dp)*dchi1(2))
          yt4 = abs(yt4 - 0.5_dp*(75.0_dp*c4-48.0_dp*c2+5.0_dp)*dchi1(2))
          yt5 = abs(yt5 + sqrt(0.625_dp)*(15.0_dp*c4-6.0_dp*c2-1.0_dp)*dchi1(2))
          yt6 = abs(yt6 - 0.25_dp*(15.0_dp*c4+30.0_dp*c2-5.0_dp)*dchi1(2))
          zt1 = abs(zt1 - 0.25_dp*(75.0_dp*c4-42.0_dp*c2+3.0_dp)*dchi1(3))
          zt2 = abs(zt2 + 3.0_dp/sqrt(8.0_dp)*(25.0_dp*c4-18.0_dp*c2+1.0_dp)*dchi1(3))
          zt3 = abs(zt3 - sqrt(45.0_dp)/4.0_dp*(5.0_dp*c4-6.0_dp*c2+1.0_dp)*dchi1(3))
          zt4 = abs(zt4 - 0.5_dp*(75.0_dp*c4-48.0_dp*c2+5.0_dp)*dchi1(3))
          zt5 = abs(zt5 + sqrt(0.625_dp)*(15.0_dp*c4-6.0_dp*c2-1.0_dp)*dchi1(3))
          zt6 = abs(zt6 - 0.25_dp*(15.0_dp*c4+30.0_dp*c2-5.0_dp)*dchi1(3))
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 x0', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 x0', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 x0', l1, l2, xt3
          IF (xt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 x0', l1, l2, xt4
          IF (xt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 x0', l1, l2, xt5
          IF (xt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 x0', l1, l2, xt6
          emax = max(emax,xt1,xt2,xt3,xt4,xt5,xt6)
          IF (yt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 y0', l1, l2, yt1
          IF (yt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 y0', l1, l2, yt2
          IF (yt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 y0', l1, l2, yt3
          IF (yt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 y0', l1, l2, yt4
          IF (yt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 y0', l1, l2, yt5
          IF (yt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 y0', l1, l2, yt6
          emax = max(emax,yt1,yt2,yt3,yt4,yt5,yt6)
          IF (zt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 z0', l1, l2, zt1
          IF (zt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 z0', l1, l2, zt2
          IF (zt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 z0', l1, l2, zt3
          IF (zt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 z0', l1, l2, zt4
          IF (zt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 z0', l1, l2, zt5
          IF (zt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 z0', l1, l2, zt6
          emax = max(emax,zt1,zt2,zt3,zt4,zt5,zt6)
          dx1 = dpro(ddll1(:,:,1,1),dll2(:,:,1),l1,l2)
          dx2 = dpro(ddll1(:,:,1,1),dll2(:,:,3),l1,l2)
          dx3 = dpro(ddll1(:,:,1,1),dll2(:,:,5),l1,l2)
          dx4 = dpro(ddll1(:,:,3,1),dll2(:,:,3),l1,l2)
          dx5 = dpro(ddll1(:,:,3,1),dll2(:,:,5),l1,l2)
          dx6 = dpro(ddll1(:,:,5,1),dll2(:,:,5),l1,l2)
          dy1 = dpro(ddll1(:,:,1,2),dll2(:,:,1),l1,l2)
          dy2 = dpro(ddll1(:,:,1,2),dll2(:,:,3),l1,l2)
          dy3 = dpro(ddll1(:,:,1,2),dll2(:,:,5),l1,l2)
          dy4 = dpro(ddll1(:,:,3,2),dll2(:,:,3),l1,l2)
          dy5 = dpro(ddll1(:,:,3,2),dll2(:,:,5),l1,l2)
          dy6 = dpro(ddll1(:,:,5,2),dll2(:,:,5),l1,l2)
          dz1 = dpro(ddll1(:,:,1,3),dll2(:,:,1),l1,l2)
          dz2 = dpro(ddll1(:,:,1,3),dll2(:,:,3),l1,l2)
          dz3 = dpro(ddll1(:,:,1,3),dll2(:,:,5),l1,l2)
          dz4 = dpro(ddll1(:,:,3,3),dll2(:,:,3),l1,l2)
          dz5 = dpro(ddll1(:,:,3,3),dll2(:,:,5),l1,l2)
          dz6 = dpro(ddll1(:,:,5,3),dll2(:,:,5),l1,l2)
          dx1 = abs(dx1-17.5_dp*dchi1(1))
          dx2 = abs(dx2+3.5_dp*dchi1(1))
          dx3 = abs(dx3-3.0_dp/sqrt(2.0_dp)*dchi1(1))
          dx4 = abs(dx4-0.5_dp*(72.0_dp*c2-13.0_dp)*dchi1(1))
          dx5 = abs(dx5-1.0_dp/sqrt(2.0_dp)*dchi1(1)*(12.0_dp*c2-3.0_dp))
          dx6 = abs(dx6-0.25_dp*dchi1(1)*(75.0_dp*c4-42.0_dp*c2+3.0_dp))
          dy1 = abs(dy1-17.5_dp*dchi1(2))
          dy2 = abs(dy2+3.5_dp*dchi1(2))
          dy3 = abs(dy3-3.0_dp/sqrt(2.0_dp)*dchi1(2))
          dy4 = abs(dy4-0.5_dp*(72.0_dp*c2-13.0_dp)*dchi1(2))
          dy5 = abs(dy5-1.0_dp/sqrt(2.0_dp)*dchi1(2)*(12.0_dp*c2-3.0_dp))
          dy6 = abs(dy6-0.25_dp*dchi1(2)*(75.0_dp*c4-42.0_dp*c2+3.0_dp))
          dz1 = abs(dz1-17.5_dp*dchi1(3))
          dz2 = abs(dz2+3.5_dp*dchi1(3))
          dz3 = abs(dz3-3.0_dp/sqrt(2.0_dp)*dchi1(3))
          dz4 = abs(dz4-0.5_dp*(72.0_dp*c2-13.0_dp)*dchi1(3))
          dz5 = abs(dz5-1.0_dp/sqrt(2.0_dp)*dchi1(3)*(12.0_dp*c2-3.0_dp))
          dz6 = abs(dz6-0.25_dp*dchi1(3)*(75.0_dp*c4-42.0_dp*c2+3.0_dp))
          IF (dx1>1.E-12) WRITE (*,*) ' Error in test for dll #1 x0', l1, l2, dx1
          IF (dx2>1.E-12) WRITE (*,*) ' Error in test for dll #2 x0', l1, l2, dx2
          IF (dx3>1.E-12) WRITE (*,*) ' Error in test for dll #3 x0', l1, l2, dx3
          IF (dx4>1.E-12) WRITE (*,*) ' Error in test for dll #4 x0', l1, l2, dx4
          IF (dx5>1.E-12) WRITE (*,*) ' Error in test for dll #5 x0', l1, l2, dx5
          IF (dx6>1.E-12) WRITE (*,*) ' Error in test for dll #6 x0', l1, l2, dx6
          emax = max(emax,dx1,dx2,dx3,dx4,dx5,dx6)
          IF (dy1>1.E-12) WRITE (*,*) ' Error in test for dll #1 y0', l1, l2, dy1
          IF (dy2>1.E-12) WRITE (*,*) ' Error in test for dll #2 y0', l1, l2, dy2
          IF (dy3>1.E-12) WRITE (*,*) ' Error in test for dll #3 y0', l1, l2, dy3
          IF (dy4>1.E-12) WRITE (*,*) ' Error in test for dll #4 y0', l1, l2, dy4
          IF (dy5>1.E-12) WRITE (*,*) ' Error in test for dll #5 y0', l1, l2, dy5
          IF (dy6>1.E-12) WRITE (*,*) ' Error in test for dll #6 y0', l1, l2, dy6
          emax = max(emax,dy1,dy2,dy3,dy4,dy5,dy6)
          IF (dz1>1.E-12) WRITE (*,*) ' Error in test for dll #1 z0', l1, l2, dz1
          IF (dz2>1.E-12) WRITE (*,*) ' Error in test for dll #2 z0', l1, l2, dz2
          IF (dz3>1.E-12) WRITE (*,*) ' Error in test for dll #3 z0', l1, l2, dz3
          IF (dz4>1.E-12) WRITE (*,*) ' Error in test for dll #4 z0', l1, l2, dz4
          IF (dz5>1.E-12) WRITE (*,*) ' Error in test for dll #5 z0', l1, l2, dz5
          IF (dz6>1.E-12) WRITE (*,*) ' Error in test for dll #6 z0', l1, l2, dz6
          emax = max(emax,dx1,dx2,dx3,dx4,dx5,dx6)
        END SELECT
        xt1 = dpro(gmu1(:,:,0),dgmu2(:,:,0,1),l1,l2)
        xt2 = dpro(gmu1(:,:,0),dgmu2(:,:,1,1),l1,l2)
        xt3 = dpro(gmu1(:,:,0),dgmu2(:,:,2,1),l1,l2)
        xt4 = dpro(gmu1(:,:,1),dgmu2(:,:,1,1),l1,l2)
        xt5 = dpro(gmu1(:,:,1),dgmu2(:,:,2,1),l1,l2)
        xt6 = dpro(gmu1(:,:,2),dgmu2(:,:,2,1),l1,l2)
        yt1 = dpro(gmu1(:,:,0),dgmu2(:,:,0,2),l1,l2)
        yt2 = dpro(gmu1(:,:,0),dgmu2(:,:,1,2),l1,l2)
        yt3 = dpro(gmu1(:,:,0),dgmu2(:,:,2,2),l1,l2)
        yt4 = dpro(gmu1(:,:,1),dgmu2(:,:,1,2),l1,l2)
        yt5 = dpro(gmu1(:,:,1),dgmu2(:,:,2,2),l1,l2)
        yt6 = dpro(gmu1(:,:,2),dgmu2(:,:,2,2),l1,l2)
        zt1 = dpro(gmu1(:,:,0),dgmu2(:,:,0,3),l1,l2)
        zt2 = dpro(gmu1(:,:,0),dgmu2(:,:,1,3),l1,l2)
        zt3 = dpro(gmu1(:,:,0),dgmu2(:,:,2,3),l1,l2)
        zt4 = dpro(gmu1(:,:,1),dgmu2(:,:,1,3),l1,l2)
        zt5 = dpro(gmu1(:,:,1),dgmu2(:,:,2,3),l1,l2)
        zt6 = dpro(gmu1(:,:,2),dgmu2(:,:,2,3),l1,l2)
        SELECT CASE (lm)
        CASE (2)
! d-dd
          xt1 = abs(xt1 - 3.0_dp*chi*(3.0_dp*chi*chi-1.0_dp)*dchi2(1))
          xt2 = abs(xt2 + 6.0_dp*chi*(2.0_dp*chi*chi-1.0_dp)*dchi2(1))
          xt3 = abs(xt3 - 3.0_dp*chi*(chi*chi-1.0_dp)*dchi2(1))
          xt4 = abs(xt4 - 2.0_dp*chi*(8.0_dp*chi*chi-3.0_dp)*dchi2(1))
          xt5 = abs(xt5 + 4.0_dp*chi*chi*chi*dchi2(1))
          xt6 = abs(xt6 - chi*(chi*chi+3.0_dp)*dchi2(1))
          yt1 = abs(yt1 - 3.0_dp*chi*(3.0_dp*chi*chi-1.0_dp)*dchi2(2))
          yt2 = abs(yt2 + 6.0_dp*chi*(2.0_dp*chi*chi-1.0_dp)*dchi2(2))
          yt3 = abs(yt3 - 3.0_dp*chi*(chi*chi-1.0_dp)*dchi2(2))
          yt4 = abs(yt4 - 2.0_dp*chi*(8.0_dp*chi*chi-3.0_dp)*dchi2(2))
          yt5 = abs(yt5 + 4.0_dp*chi*chi*chi*dchi2(2))
          yt6 = abs(yt6 - chi*(chi*chi+3.0_dp)*dchi2(2))
          zt1 = abs(zt1 - 3.0_dp*chi*(3.0_dp*chi*chi-1.0_dp)*dchi2(3))
          zt2 = abs(zt2 + 6.0_dp*chi*(2.0_dp*chi*chi-1.0_dp)*dchi2(3))
          zt3 = abs(zt3 - 3.0_dp*chi*(chi*chi-1.0_dp)*dchi2(3))
          zt4 = abs(zt4 - 2.0_dp*chi*(8.0_dp*chi*chi-3.0_dp)*dchi2(3))
          zt5 = abs(zt5 + 4.0_dp*chi*chi*chi*dchi2(3))
          zt6 = abs(zt6 - chi*(chi*chi+3.0_dp)*dchi2(3))
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0x', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0x', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0x', l1, l2, xt3
          IF (xt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 0x', l1, l2, xt4
          IF (xt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 0x', l1, l2, xt5
          IF (xt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 0x', l1, l2, xt6
          emax = max(emax,xt1,xt2,xt3,xt4,xt5,xt6)
          IF (yt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0y', l1, l2, yt1
          IF (yt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0y', l1, l2, yt2
          IF (yt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0y', l1, l2, yt3
          IF (yt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 0y', l1, l2, yt4
          IF (yt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 0y', l1, l2, yt5
          IF (yt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 0y', l1, l2, yt6
          emax = max(emax,yt1,yt2,yt3,yt4,yt5,yt6)
          IF (zt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0z', l1, l2, zt1
          IF (zt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0z', l1, l2, zt2
          IF (zt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0z', l1, l2, zt3
          IF (zt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 0z', l1, l2, zt4
          IF (zt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 0z', l1, l2, zt5
          IF (zt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 0z', l1, l2, zt6
          emax = max(emax,zt1,zt2,zt3,zt4,zt5,zt6)
          dx1 = dpro(dll1(:,:,0),ddll2(:,:,0,1),l1,l2)
          dx2 = dpro(dll1(:,:,0),ddll2(:,:,2,1),l1,l2)
          dx3 = dpro(dll1(:,:,0),ddll2(:,:,4,1),l1,l2)
          dx4 = dpro(dll1(:,:,2),ddll2(:,:,2,1),l1,l2)
          dx5 = dpro(dll1(:,:,2),ddll2(:,:,4,1),l1,l2)
          dx6 = dpro(dll1(:,:,4),ddll2(:,:,4,1),l1,l2)
          dy1 = dpro(dll1(:,:,0),ddll2(:,:,0,2),l1,l2)
          dy2 = dpro(dll1(:,:,0),ddll2(:,:,2,2),l1,l2)
          dy3 = dpro(dll1(:,:,0),ddll2(:,:,4,2),l1,l2)
          dy4 = dpro(dll1(:,:,2),ddll2(:,:,2,2),l1,l2)
          dy5 = dpro(dll1(:,:,2),ddll2(:,:,4,2),l1,l2)
          dy6 = dpro(dll1(:,:,4),ddll2(:,:,4,2),l1,l2)
          dz1 = dpro(dll1(:,:,0),ddll2(:,:,0,3),l1,l2)
          dz2 = dpro(dll1(:,:,0),ddll2(:,:,2,3),l1,l2)
          dz3 = dpro(dll1(:,:,0),ddll2(:,:,4,3),l1,l2)
          dz4 = dpro(dll1(:,:,2),ddll2(:,:,2,3),l1,l2)
          dz5 = dpro(dll1(:,:,2),ddll2(:,:,4,3),l1,l2)
          dz6 = dpro(dll1(:,:,4),ddll2(:,:,4,3),l1,l2)
          dx1 = abs(dx1)
          dx2 = abs(dx2)
          dx3 = abs(dx3)
          dx4 = abs(dx4-14.0_dp/3.0_dp*chi*dchi2(1))
          dx5 = abs(dx5-2.0_dp*chi*dchi2(1))
          dx6 = abs(dx6-3.0_dp*(3.0_dp*chi*chi-1.0_dp)*chi*dchi2(1))
          dy1 = abs(dy1)
          dy2 = abs(dy2)
          dy3 = abs(dy3)
          dy4 = abs(dy4-14.0_dp/3.0_dp*chi*dchi2(2))
          dy5 = abs(dy5-2.0_dp*chi*dchi2(2))
          dy6 = abs(dy6-3.0_dp*(3.0_dp*chi*chi-1.0_dp)*chi*dchi2(2))
          dz1 = abs(dz1)
          dz2 = abs(dz2)
          dz3 = abs(dz3)
          dz4 = abs(dz4-14.0_dp/3.0_dp*chi*dchi2(3))
          dz5 = abs(dz5-2.0_dp*chi*dchi2(3))
          dz6 = abs(dz6-3.0_dp*(3.0_dp*chi*chi-1.0_dp)*chi*dchi2(3))
          IF (dx1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0x', l1, l2, dx1
          IF (dx2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0x', l1, l2, dx2
          IF (dx3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0x', l1, l2, dx3
          IF (dx4>1.E-12) WRITE (*,*) ' Error in test for dll #4 0x', l1, l2, dx4
          IF (dx5>1.E-12) WRITE (*,*) ' Error in test for dll #5 0x', l1, l2, dx5
          IF (dx6>1.E-12) WRITE (*,*) ' Error in test for dll #6 0x', l1, l2, dx6
          emax = max(emax,dx1,dx2,dx3,dx4,dx5,dx6)
          IF (dy1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0y', l1, l2, dy1
          IF (dy2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0y', l1, l2, dy2
          IF (dy3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0y', l1, l2, dy3
          IF (dy4>1.E-12) WRITE (*,*) ' Error in test for dll #4 0y', l1, l2, dy4
          IF (dy5>1.E-12) WRITE (*,*) ' Error in test for dll #5 0y', l1, l2, dy5
          IF (dy6>1.E-12) WRITE (*,*) ' Error in test for dll #6 0y', l1, l2, dy6
          emax = max(emax,dy1,dy2,dy3,dy4,dy5,dy6)
          IF (dz1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0z', l1, l2, dz1
          IF (dz2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0z', l1, l2, dz2
          IF (dz3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0z', l1, l2, dz3
          IF (dz4>1.E-12) WRITE (*,*) ' Error in test for dll #4 0z', l1, l2, dz4
          IF (dz5>1.E-12) WRITE (*,*) ' Error in test for dll #5 0z', l1, l2, dz5
          IF (dz6>1.E-12) WRITE (*,*) ' Error in test for dll #6 0z', l1, l2, dz6
          emax = max(emax,dx1,dx2,dx3,dx4,dx5,dx6)
        CASE (3)
! d-df
          c2 = chi*chi
          c4 = c2*c2
          xt1 = abs(xt1 - 0.25_dp*(75.0_dp*c4-42.0_dp*c2+3.0_dp)*dchi2(1))
          xt2 = abs(xt2 + 3.0_dp/sqrt(8.0_dp)*(25.0_dp*c4-18.0_dp*c2+1.0_dp)*dchi2(1))
          xt3 = abs(xt3 - sqrt(45.0_dp)/4.0_dp*(5.0_dp*c4-6.0_dp*c2+1.0_dp)*dchi2(1))
          xt4 = abs(xt4 - 0.5_dp*(75.0_dp*c4-48.0_dp*c2+5.0_dp)*dchi2(1))
          xt5 = abs(xt5 + sqrt(0.625_dp)*(15.0_dp*c4-6.0_dp*c2-1.0_dp)*dchi2(1))
          xt6 = abs(xt6 - 0.25_dp*(15.0_dp*c4+30.0_dp*c2-5.0_dp)*dchi2(1))
          yt1 = abs(yt1 - 0.25_dp*(75.0_dp*c4-42.0_dp*c2+3.0_dp)*dchi2(2))
          yt2 = abs(yt2 + 3.0_dp/sqrt(8.0_dp)*(25.0_dp*c4-18.0_dp*c2+1.0_dp)*dchi2(2))
          yt3 = abs(yt3 - sqrt(45.0_dp)/4.0_dp*(5.0_dp*c4-6.0_dp*c2+1.0_dp)*dchi2(2))
          yt4 = abs(yt4 - 0.5_dp*(75.0_dp*c4-48.0_dp*c2+5.0_dp)*dchi2(2))
          yt5 = abs(yt5 + sqrt(0.625_dp)*(15.0_dp*c4-6.0_dp*c2-1.0_dp)*dchi2(2))
          yt6 = abs(yt6 - 0.25_dp*(15.0_dp*c4+30.0_dp*c2-5.0_dp)*dchi2(2))
          zt1 = abs(zt1 - 0.25_dp*(75.0_dp*c4-42.0_dp*c2+3.0_dp)*dchi2(3))
          zt2 = abs(zt2 + 3.0_dp/sqrt(8.0_dp)*(25.0_dp*c4-18.0_dp*c2+1.0_dp)*dchi2(3))
          zt3 = abs(zt3 - sqrt(45.0_dp)/4.0_dp*(5.0_dp*c4-6.0_dp*c2+1.0_dp)*dchi2(3))
          zt4 = abs(zt4 - 0.5_dp*(75.0_dp*c4-48.0_dp*c2+5.0_dp)*dchi2(3))
          zt5 = abs(zt5 + sqrt(0.625_dp)*(15.0_dp*c4-6.0_dp*c2-1.0_dp)*dchi2(3))
          zt6 = abs(zt6 - 0.25_dp*(15.0_dp*c4+30.0_dp*c2-5.0_dp)*dchi2(3))
          IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0x', l1, l2, xt1
          IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0x', l1, l2, xt2
          IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0x', l1, l2, xt3
          IF (xt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 0x', l1, l2, xt4
          IF (xt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 0x', l1, l2, xt5
          IF (xt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 0x', l1, l2, xt6
          emax = max(emax,xt1,xt2,xt3,xt4,xt5,xt6)
          IF (yt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0y', l1, l2, yt1
          IF (yt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0y', l1, l2, yt2
          IF (yt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0y', l1, l2, yt3
          IF (yt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 0y', l1, l2, yt4
          IF (yt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 0y', l1, l2, yt5
          IF (yt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 0y', l1, l2, yt6
          emax = max(emax,yt1,yt2,yt3,yt4,yt5,yt6)
          IF (zt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0z', l1, l2, zt1
          IF (zt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0z', l1, l2, zt2
          IF (zt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0z', l1, l2, zt3
          IF (zt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 0z', l1, l2, zt4
          IF (zt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 0z', l1, l2, zt5
          IF (zt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 0z', l1, l2, zt6
          emax = max(emax,zt1,zt2,zt3,zt4,zt5,zt6)
          dx1 = dpro(dll1(:,:,1),ddll2(:,:,1,1),l1,l2)
          dx2 = dpro(dll1(:,:,1),ddll2(:,:,3,1),l1,l2)
          dx3 = dpro(dll1(:,:,1),ddll2(:,:,5,1),l1,l2)
          dx4 = dpro(dll1(:,:,3),ddll2(:,:,3,1),l1,l2)
          dx5 = dpro(dll1(:,:,3),ddll2(:,:,5,1),l1,l2)
          dx6 = dpro(dll1(:,:,5),ddll2(:,:,5,1),l1,l2)
          dy1 = dpro(dll1(:,:,1),ddll2(:,:,1,2),l1,l2)
          dy2 = dpro(dll1(:,:,1),ddll2(:,:,3,2),l1,l2)
          dy3 = dpro(dll1(:,:,1),ddll2(:,:,5,2),l1,l2)
          dy4 = dpro(dll1(:,:,3),ddll2(:,:,3,2),l1,l2)
          dy5 = dpro(dll1(:,:,3),ddll2(:,:,5,2),l1,l2)
          dy6 = dpro(dll1(:,:,5),ddll2(:,:,5,2),l1,l2)
          dz1 = dpro(dll1(:,:,1),ddll2(:,:,1,3),l1,l2)
          dz2 = dpro(dll1(:,:,1),ddll2(:,:,3,3),l1,l2)
          dz3 = dpro(dll1(:,:,1),ddll2(:,:,5,3),l1,l2)
          dz4 = dpro(dll1(:,:,3),ddll2(:,:,3,3),l1,l2)
          dz5 = dpro(dll1(:,:,3),ddll2(:,:,5,3),l1,l2)
          dz6 = dpro(dll1(:,:,5),ddll2(:,:,5,3),l1,l2)
          dx1 = abs(dx1-17.5_dp*dchi2(1))
          dx2 = abs(dx2+3.5_dp*dchi2(1))
          dx3 = abs(dx3-3.0_dp/sqrt(2.0_dp)*dchi2(1))
          dx4 = abs(dx4-0.5_dp*(72.0_dp*c2-13.0_dp)*dchi2(1))
          dx5 = abs(dx5-1.0_dp/sqrt(2.0_dp)*dchi2(1)*(12.0_dp*c2-3.0_dp))
          dx6 = abs(dx6-0.25_dp*dchi2(1)*(75.0_dp*c4-42.0_dp*c2+3.0_dp))
          dy1 = abs(dy1-17.5_dp*dchi2(2))
          dy2 = abs(dy2+3.5_dp*dchi2(2))
          dy3 = abs(dy3-3.0_dp/sqrt(2.0_dp)*dchi2(2))
          dy4 = abs(dy4-0.5_dp*(72.0_dp*c2-13.0_dp)*dchi2(2))
          dy5 = abs(dy5-1.0_dp/sqrt(2.0_dp)*dchi2(2)*(12.0_dp*c2-3.0_dp))
          dy6 = abs(dy6-0.25_dp*dchi2(2)*(75.0_dp*c4-42.0_dp*c2+3.0_dp))
          dz1 = abs(dz1-17.5_dp*dchi2(3))
          dz2 = abs(dz2+3.5_dp*dchi2(3))
          dz3 = abs(dz3-3.0_dp/sqrt(2.0_dp)*dchi2(3))
          dz4 = abs(dz4-0.5_dp*(72.0_dp*c2-13.0_dp)*dchi2(3))
          dz5 = abs(dz5-1.0_dp/sqrt(2.0_dp)*dchi2(3)*(12.0_dp*c2-3.0_dp))
          dz6 = abs(dz6-0.25_dp*dchi2(3)*(75.0_dp*c4-42.0_dp*c2+3.0_dp))
          IF (dx1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0x', l1, l2, dx1
          IF (dx2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0x', l1, l2, dx2
          IF (dx3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0x', l1, l2, dx3
          IF (dx4>1.E-12) WRITE (*,*) ' Error in test for dll #4 0x', l1, l2, dx4
          IF (dx5>1.E-12) WRITE (*,*) ' Error in test for dll #5 0x', l1, l2, dx5
          IF (dx6>1.E-12) WRITE (*,*) ' Error in test for dll #6 0x', l1, l2, dx6
          emax = max(emax,dx1,dx2,dx3,dx4,dx5,dx6)
          IF (dy1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0y', l1, l2, dy1
          IF (dy2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0y', l1, l2, dy2
          IF (dy3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0y', l1, l2, dy3
          IF (dy4>1.E-12) WRITE (*,*) ' Error in test for dll #4 0y', l1, l2, dy4
          IF (dy5>1.E-12) WRITE (*,*) ' Error in test for dll #5 0y', l1, l2, dy5
          IF (dy6>1.E-12) WRITE (*,*) ' Error in test for dll #6 0y', l1, l2, dy6
          emax = max(emax,dy1,dy2,dy3,dy4,dy5,dy6)
          IF (dz1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0z', l1, l2, dz1
          IF (dz2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0z', l1, l2, dz2
          IF (dz3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0z', l1, l2, dz3
          IF (dz4>1.E-12) WRITE (*,*) ' Error in test for dll #4 0z', l1, l2, dz4
          IF (dz5>1.E-12) WRITE (*,*) ' Error in test for dll #5 0z', l1, l2, dz5
          IF (dz6>1.E-12) WRITE (*,*) ' Error in test for dll #6 0z', l1, l2, dz6
          emax = max(emax,dx1,dx2,dx3,dx4,dx5,dx6)
        END SELECT
      END SUBROUTINE dd_test
      SUBROUTINE df_test(chi,dchi1,dchi2,gmu1,gmu2,dgmu1,dgmu2, &
                         dll1,dll2,ddll1,ddll2,l1,l2,emax)
    REAL(KIND=dp), INTENT(IN) :: chi, dchi1(3), dchi2(3), gmu1(0:6,0:6,0:3), &
      gmu2(0:6,0:6,0:3), dgmu1(0:6,0:6,0:3,1:3), dgmu2(0:6,0:6,0:3,1:3), &
      dll1(0:6,0:6,0:6), dll2(0:6,0:6,0:6), ddll1(0:6,0:6,0:6,1:3), &
      ddll2(0:6,0:6,0:6,1:3)
    INTEGER, INTENT(IN)                      :: l1, l2
    REAL(KIND=dp), INTENT(OUT)                   :: emax

    INTEGER                                  :: lm, lu
    REAL(KIND=dp) :: c2, c4, dx0, dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8, dx9, &
      dy0, dy1, dy2, dy3, dy4, dy5, dy6, dy7, dy8, dy9, dz0, dz1, dz2, dz3, &
      dz4, dz5, dz6, dz7, dz8, dz9, xt0, xt1, xt2, xt3, xt4, xt5, xt6, xt7, &
      xt8, xt9, yt0, yt1, yt2, yt3, yt4, yt5, yt6, yt7, yt8, yt9, zt0, zt1, &
      zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9

        emax = 0.0_dp
        lm = max(l1,l2)
        lu = min(l1,l2)
        c2 = chi*chi
        c4 = c2*c2
! df-f
        xt1 = dpro(dgmu1(:,:,0,1),gmu2(:,:,0),l1,l2)
        xt2 = dpro(dgmu1(:,:,0,1),gmu2(:,:,1),l1,l2)
        xt3 = dpro(dgmu1(:,:,0,1),gmu2(:,:,2),l1,l2)
        xt4 = dpro(dgmu1(:,:,0,1),gmu2(:,:,3),l1,l2)
        xt5 = dpro(dgmu1(:,:,1,1),gmu2(:,:,1),l1,l2)
        xt6 = dpro(dgmu1(:,:,1,1),gmu2(:,:,2),l1,l2)
        xt7 = dpro(dgmu1(:,:,1,1),gmu2(:,:,3),l1,l2)
        xt8 = dpro(dgmu1(:,:,2,1),gmu2(:,:,2),l1,l2)
        xt9 = dpro(dgmu1(:,:,2,1),gmu2(:,:,3),l1,l2)
        xt0 = dpro(dgmu1(:,:,3,1),gmu2(:,:,3),l1,l2)
        yt1 = dpro(dgmu1(:,:,0,2),gmu2(:,:,0),l1,l2)
        yt2 = dpro(dgmu1(:,:,0,2),gmu2(:,:,1),l1,l2)
        yt3 = dpro(dgmu1(:,:,0,2),gmu2(:,:,2),l1,l2)
        yt4 = dpro(dgmu1(:,:,0,2),gmu2(:,:,3),l1,l2)
        yt5 = dpro(dgmu1(:,:,1,2),gmu2(:,:,1),l1,l2)
        yt6 = dpro(dgmu1(:,:,1,2),gmu2(:,:,2),l1,l2)
        yt7 = dpro(dgmu1(:,:,1,2),gmu2(:,:,3),l1,l2)
        yt8 = dpro(dgmu1(:,:,2,2),gmu2(:,:,2),l1,l2)
        yt9 = dpro(dgmu1(:,:,2,2),gmu2(:,:,3),l1,l2)
        yt0 = dpro(dgmu1(:,:,3,2),gmu2(:,:,3),l1,l2)
        zt1 = dpro(dgmu1(:,:,0,3),gmu2(:,:,0),l1,l2)
        zt2 = dpro(dgmu1(:,:,0,3),gmu2(:,:,1),l1,l2)
        zt3 = dpro(dgmu1(:,:,0,3),gmu2(:,:,2),l1,l2)
        zt4 = dpro(dgmu1(:,:,0,3),gmu2(:,:,3),l1,l2)
        zt5 = dpro(dgmu1(:,:,1,3),gmu2(:,:,1),l1,l2)
        zt6 = dpro(dgmu1(:,:,1,3),gmu2(:,:,2),l1,l2)
        zt7 = dpro(dgmu1(:,:,1,3),gmu2(:,:,3),l1,l2)
        zt8 = dpro(dgmu1(:,:,2,3),gmu2(:,:,2),l1,l2)
        zt9 = dpro(dgmu1(:,:,2,3),gmu2(:,:,3),l1,l2)
        zt0 = dpro(dgmu1(:,:,3,3),gmu2(:,:,3),l1,l2)
        xt1 = abs(xt1 - 0.25_dp*chi*(150.0_dp*c4-120.0_dp*c2+18.0_dp)*dchi1(1))
        xt2 = abs(xt2 + 0.375_dp*chi*(150.0_dp*c4-140.0_dp*c2+22.0_dp)*dchi1(1))
        xt3 = abs(xt3 - 3.75_dp*chi*(6.0_dp*c4-8.0_dp*c2+2.0_dp)*dchi1(1))
        xt4 = abs(xt4 + 0.625_dp*chi*(6.0_dp*c4-12.0_dp*c2+6.0_dp)*dchi1(1))
        xt5 = abs(xt5 - 0.0625_dp*chi*(1350.0_dp*c4-1220.0_dp*c2+222.0_dp)*dchi1(1))
        xt6 = abs(xt6 + 0.625_dp*chi*(54.0_dp*c4-44.0_dp*c2+6.0_dp)*dchi1(1))
        xt7 = abs(xt7 - 0.9375_dp*chi*(6.0_dp*c4-4.0_dp*c2-2.0_dp)*dchi1(1))
        xt8 = abs(xt8 - 0.25_dp*chi*(54.0_dp*c4+40.0_dp*c2-30.0_dp)*dchi1(1))
        xt9 = abs(xt9 + 0.375_dp*chi*(6.0_dp*c4+20.0_dp*c2-10.0_dp)*dchi1(1))
        xt0 = abs(xt0 - 0.0625_dp*chi*(6.0_dp*c4+60.0_dp*c2+30.0_dp)*dchi1(1))
        yt1 = abs(yt1 - 0.25_dp*chi*(150.0_dp*c4-120.0_dp*c2+18.0_dp)*dchi1(2))
        yt2 = abs(yt2 + 0.375_dp*chi*(150.0_dp*c4-140.0_dp*c2+22.0_dp)*dchi1(2))
        yt3 = abs(yt3 - 3.75_dp*chi*(6.0_dp*c4-8.0_dp*c2+2.0_dp)*dchi1(2))
        yt4 = abs(yt4 + 0.625_dp*chi*(6.0_dp*c4-12.0_dp*c2+6.0_dp)*dchi1(2))
        yt5 = abs(yt5 - 0.0625_dp*chi*(1350.0_dp*c4-1220.0_dp*c2+222.0_dp)*dchi1(2))
        yt6 = abs(yt6 + 0.625_dp*chi*(54.0_dp*c4-44.0_dp*c2+6.0_dp)*dchi1(2))
        yt7 = abs(yt7 - 0.9375_dp*chi*(6.0_dp*c4-4.0_dp*c2-2.0_dp)*dchi1(2))
        yt8 = abs(yt8 - 0.25_dp*chi*(54.0_dp*c4+40.0_dp*c2-30.0_dp)*dchi1(2))
        yt9 = abs(yt9 + 0.375_dp*chi*(6.0_dp*c4+20.0_dp*c2-10.0_dp)*dchi1(2))
        yt0 = abs(yt0 - 0.0625_dp*chi*(6.0_dp*c4+60.0_dp*c2+30.0_dp)*dchi1(2))
        zt1 = abs(zt1 - 0.25_dp*chi*(150.0_dp*c4-120.0_dp*c2+18.0_dp)*dchi1(3))
        zt2 = abs(zt2 + 0.375_dp*chi*(150.0_dp*c4-140.0_dp*c2+22.0_dp)*dchi1(3))
        zt3 = abs(zt3 - 3.75_dp*chi*(6.0_dp*c4-8.0_dp*c2+2.0_dp)*dchi1(3))
        zt4 = abs(zt4 + 0.625_dp*chi*(6.0_dp*c4-12.0_dp*c2+6.0_dp)*dchi1(3))
        zt5 = abs(zt5 - 0.0625_dp*chi*(1350.0_dp*c4-1220.0_dp*c2+222.0_dp)*dchi1(3))
        zt6 = abs(zt6 + 0.625_dp*chi*(54.0_dp*c4-44.0_dp*c2+6.0_dp)*dchi1(3))
        zt7 = abs(zt7 - 0.9375_dp*chi*(6.0_dp*c4-4.0_dp*c2-2.0_dp)*dchi1(3))
        zt8 = abs(zt8 - 0.25_dp*chi*(54.0_dp*c4+40.0_dp*c2-30.0_dp)*dchi1(3))
        zt9 = abs(zt9 + 0.375_dp*chi*(6.0_dp*c4+20.0_dp*c2-10.0_dp)*dchi1(3))
        zt0 = abs(zt0 - 0.0625_dp*chi*(6.0_dp*c4+60.0_dp*c2+30.0_dp)*dchi1(3))
        IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 x0', l1, l2, xt1
        IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 x0', l1, l2, xt2
        IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 x0', l1, l2, xt3
        IF (xt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 x0', l1, l2, xt4
        IF (xt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 x0', l1, l2, xt5
        IF (xt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 x0', l1, l2, xt6
        IF (xt7>1.E-12) WRITE (*,*) ' Error in test for gmu #7 x0', l1, l2, xt7
        IF (xt8>1.E-12) WRITE (*,*) ' Error in test for gmu #8 x0', l1, l2, xt8
        IF (xt9>1.E-12) WRITE (*,*) ' Error in test for gmu #9 x0', l1, l2, xt9
        IF (xt0>1.E-12) WRITE (*,*) ' Error in test for gmu #0 x0', l1, l2, xt0
        emax = max(emax,xt1,xt2,xt3,xt4,xt5,xt6,xt7,xt8,xt9,xt0)
        IF (yt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 y0', l1, l2, yt1
        IF (yt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 y0', l1, l2, yt2
        IF (yt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 y0', l1, l2, yt3
        IF (yt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 y0', l1, l2, yt4
        IF (yt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 y0', l1, l2, yt5
        IF (yt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 y0', l1, l2, yt6
        IF (yt7>1.E-12) WRITE (*,*) ' Error in test for gmu #7 y0', l1, l2, yt7
        IF (yt8>1.E-12) WRITE (*,*) ' Error in test for gmu #8 y0', l1, l2, yt8
        IF (yt9>1.E-12) WRITE (*,*) ' Error in test for gmu #9 y0', l1, l2, yt9
        IF (yt0>1.E-12) WRITE (*,*) ' Error in test for gmu #0 y0', l1, l2, yt0
        emax = max(emax,yt1,yt2,yt3,yt4,yt5,yt6,yt7,yt8,yt9,yt0)
        IF (zt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 z0', l1, l2, zt1
        IF (zt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 z0', l1, l2, zt2
        IF (zt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 z0', l1, l2, zt3
        IF (zt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 z0', l1, l2, zt4
        IF (zt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 z0', l1, l2, zt5
        IF (zt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 z0', l1, l2, zt6
        IF (zt7>1.E-12) WRITE (*,*) ' Error in test for gmu #7 z0', l1, l2, zt7
        IF (zt8>1.E-12) WRITE (*,*) ' Error in test for gmu #8 z0', l1, l2, zt8
        IF (zt9>1.E-12) WRITE (*,*) ' Error in test for gmu #9 z0', l1, l2, zt9
        IF (zt0>1.E-12) WRITE (*,*) ' Error in test for gmu #0 z0', l1, l2, zt0
        emax = max(emax,zt1,zt2,zt3,zt4,zt5,zt6,zt7,zt8,zt9,zt0)
        dx1 = dpro(ddll1(:,:,0,1),dll2(:,:,0),l1,l2)
        dx2 = dpro(ddll1(:,:,0,1),dll2(:,:,2),l1,l2)
        dx3 = dpro(ddll1(:,:,0,1),dll2(:,:,4),l1,l2)
        dx4 = dpro(ddll1(:,:,0,1),dll2(:,:,6),l1,l2)
        dx5 = dpro(ddll1(:,:,2,1),dll2(:,:,2),l1,l2)
        dx6 = dpro(ddll1(:,:,2,1),dll2(:,:,4),l1,l2)
        dx7 = dpro(ddll1(:,:,2,1),dll2(:,:,6),l1,l2)
        dx8 = dpro(ddll1(:,:,4,1),dll2(:,:,4),l1,l2)
        dx9 = dpro(ddll1(:,:,4,1),dll2(:,:,6),l1,l2)
        dx0 = dpro(ddll1(:,:,6,1),dll2(:,:,6),l1,l2)
        dy1 = dpro(ddll1(:,:,0,2),dll2(:,:,0),l1,l2)
        dy2 = dpro(ddll1(:,:,0,2),dll2(:,:,2),l1,l2)
        dy3 = dpro(ddll1(:,:,0,2),dll2(:,:,4),l1,l2)
        dy4 = dpro(ddll1(:,:,0,2),dll2(:,:,6),l1,l2)
        dy5 = dpro(ddll1(:,:,2,2),dll2(:,:,2),l1,l2)
        dy6 = dpro(ddll1(:,:,2,2),dll2(:,:,4),l1,l2)
        dy7 = dpro(ddll1(:,:,2,2),dll2(:,:,6),l1,l2)
        dy8 = dpro(ddll1(:,:,4,2),dll2(:,:,4),l1,l2)
        dy9 = dpro(ddll1(:,:,4,2),dll2(:,:,6),l1,l2)
        dy0 = dpro(ddll1(:,:,6,2),dll2(:,:,6),l1,l2)
        dz1 = dpro(ddll1(:,:,0,3),dll2(:,:,0),l1,l2)
        dz2 = dpro(ddll1(:,:,0,3),dll2(:,:,2),l1,l2)
        dz3 = dpro(ddll1(:,:,0,3),dll2(:,:,4),l1,l2)
        dz4 = dpro(ddll1(:,:,0,3),dll2(:,:,6),l1,l2)
        dz5 = dpro(ddll1(:,:,2,3),dll2(:,:,2),l1,l2)
        dz6 = dpro(ddll1(:,:,2,3),dll2(:,:,4),l1,l2)
        dz7 = dpro(ddll1(:,:,2,3),dll2(:,:,6),l1,l2)
        dz8 = dpro(ddll1(:,:,4,3),dll2(:,:,4),l1,l2)
        dz9 = dpro(ddll1(:,:,4,3),dll2(:,:,6),l1,l2)
        dz0 = dpro(ddll1(:,:,6,3),dll2(:,:,6),l1,l2)
        dx1 = abs(dx1)
        dx2 = abs(dx2)
        dx3 = abs(dx3)
        dx4 = abs(dx4)
        dx5 = abs(dx5-10.08_dp*chi*dchi1(1))
        dx6 = abs(dx6-1.44_dp*chi*dchi1(1))
        dx7 = abs(dx7-2.4_dp*chi*dchi1(1))
        dx8 = abs(dx8-0.04_dp*chi*(880.0_dp*c2-372.0_dp)*dchi1(1))
        dx9 = abs(dx9-0.6_dp*chi*(20.0_dp*c2-8.0_dp)*dchi1(1))
        dx0 = abs(dx0-0.25_dp*chi*(150.0_dp*c4-120.0_dp*c2+18.0_dp)*dchi1(1))
        dy1 = abs(dy1)
        dy2 = abs(dy2)
        dy3 = abs(dy3)
        dy4 = abs(dy4)
        dy5 = abs(dy5-10.08_dp*chi*dchi1(2))
        dy6 = abs(dy6-1.44_dp*chi*dchi1(2))
        dy7 = abs(dy7-2.4_dp*chi*dchi1(2))
        dy8 = abs(dy8-0.04_dp*chi*(880.0_dp*c2-372.0_dp)*dchi1(2))
        dy9 = abs(dy9-0.6_dp*chi*(20.0_dp*c2-8.0_dp)*dchi1(2))
        dy0 = abs(dy0-0.25_dp*chi*(150.0_dp*c4-120.0_dp*c2+18.0_dp)*dchi1(2))
        dz1 = abs(dz1)
        dz2 = abs(dz2)
        dz3 = abs(dz3)
        dz4 = abs(dz4)
        dz5 = abs(dz5-10.08_dp*chi*dchi1(3))
        dz6 = abs(dz6-1.44_dp*chi*dchi1(3))
        dz7 = abs(dz7-2.4_dp*chi*dchi1(3))
        dz8 = abs(dz8-0.04_dp*chi*(880.0_dp*c2-372.0_dp)*dchi1(3))
        dz9 = abs(dz9-0.6_dp*chi*(20.0_dp*c2-8.0_dp)*dchi1(3))
        dz0 = abs(dz0-0.25_dp*chi*(150.0_dp*c4-120.0_dp*c2+18.0_dp)*dchi1(3))
        IF (dx1>1.E-12) WRITE (*,*) ' Error in test for dll #1 x0', l1, l2, dx1
        IF (dx2>1.E-12) WRITE (*,*) ' Error in test for dll #2 x0', l1, l2, dx2
        IF (dx3>1.E-12) WRITE (*,*) ' Error in test for dll #3 x0', l1, l2, dx3
        IF (dx4>1.E-12) WRITE (*,*) ' Error in test for dll #4 x0', l1, l2, dx4
        IF (dx5>1.E-12) WRITE (*,*) ' Error in test for dll #5 x0', l1, l2, dx5
        IF (dx6>1.E-12) WRITE (*,*) ' Error in test for dll #6 x0', l1, l2, dx6
        IF (dx7>1.E-12) WRITE (*,*) ' Error in test for dll #7 x0', l1, l2, dx7
        IF (dx8>1.E-12) WRITE (*,*) ' Error in test for dll #8 x0', l1, l2, dx8
        IF (dx9>1.E-12) WRITE (*,*) ' Error in test for dll #9 x0', l1, l2, dx9
        IF (dx0>1.E-12) WRITE (*,*) ' Error in test for dll #0 x0', l1, l2, dx0
        emax = max(emax,dx1,dx2,dx3,dx4,dx5,dx6,dx7,dx8,dx9,dx0)
        IF (dy1>1.E-12) WRITE (*,*) ' Error in test for dll #1 y0', l1, l2, dy1
        IF (dy2>1.E-12) WRITE (*,*) ' Error in test for dll #2 y0', l1, l2, dy2
        IF (dy3>1.E-12) WRITE (*,*) ' Error in test for dll #3 y0', l1, l2, dy3
        IF (dy4>1.E-12) WRITE (*,*) ' Error in test for dll #4 y0', l1, l2, dy4
        IF (dy5>1.E-12) WRITE (*,*) ' Error in test for dll #5 y0', l1, l2, dy5
        IF (dy6>1.E-12) WRITE (*,*) ' Error in test for dll #6 y0', l1, l2, dy6
        IF (dy7>1.E-12) WRITE (*,*) ' Error in test for dll #7 y0', l1, l2, dy7
        IF (dy8>1.E-12) WRITE (*,*) ' Error in test for dll #8 y0', l1, l2, dy8
        IF (dy9>1.E-12) WRITE (*,*) ' Error in test for dll #9 y0', l1, l2, dy9
        IF (dy0>1.E-12) WRITE (*,*) ' Error in test for dll #0 y0', l1, l2, dy0
        emax = max(emax,dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8,dy9,dy0)
        IF (dz1>1.E-12) WRITE (*,*) ' Error in test for dll #1 z0', l1, l2, dz1
        IF (dz2>1.E-12) WRITE (*,*) ' Error in test for dll #2 z0', l1, l2, dz2
        IF (dz3>1.E-12) WRITE (*,*) ' Error in test for dll #3 z0', l1, l2, dz3
        IF (dz4>1.E-12) WRITE (*,*) ' Error in test for dll #4 z0', l1, l2, dz4
        IF (dz5>1.E-12) WRITE (*,*) ' Error in test for dll #5 z0', l1, l2, dz5
        IF (dz6>1.E-12) WRITE (*,*) ' Error in test for dll #6 z0', l1, l2, dz6
        IF (dz7>1.E-12) WRITE (*,*) ' Error in test for dll #7 z0', l1, l2, dz7
        IF (dz8>1.E-12) WRITE (*,*) ' Error in test for dll #8 z0', l1, l2, dz8
        IF (dz9>1.E-12) WRITE (*,*) ' Error in test for dll #9 z0', l1, l2, dz9
        IF (dz0>1.E-12) WRITE (*,*) ' Error in test for dll #0 z0', l1, l2, dz0
        emax = max(emax,dz1,dz2,dz3,dz4,dz5,dz6,dz7,dz8,dz9,dz0)
! f-df
        xt1 = dpro(gmu1(:,:,0),dgmu2(:,:,0,1),l1,l2)
        xt2 = dpro(gmu1(:,:,0),dgmu2(:,:,1,1),l1,l2)
        xt3 = dpro(gmu1(:,:,0),dgmu2(:,:,2,1),l1,l2)
        xt4 = dpro(gmu1(:,:,0),dgmu2(:,:,3,1),l1,l2)
        xt5 = dpro(gmu1(:,:,1),dgmu2(:,:,1,1),l1,l2)
        xt6 = dpro(gmu1(:,:,1),dgmu2(:,:,2,1),l1,l2)
        xt7 = dpro(gmu1(:,:,1),dgmu2(:,:,3,1),l1,l2)
        xt8 = dpro(gmu1(:,:,2),dgmu2(:,:,2,1),l1,l2)
        xt9 = dpro(gmu1(:,:,2),dgmu2(:,:,3,1),l1,l2)
        xt0 = dpro(gmu1(:,:,3),dgmu2(:,:,3,1),l1,l2)
        yt1 = dpro(gmu1(:,:,0),dgmu2(:,:,0,2),l1,l2)
        yt2 = dpro(gmu1(:,:,0),dgmu2(:,:,1,2),l1,l2)
        yt3 = dpro(gmu1(:,:,0),dgmu2(:,:,2,2),l1,l2)
        yt4 = dpro(gmu1(:,:,0),dgmu2(:,:,3,2),l1,l2)
        yt5 = dpro(gmu1(:,:,1),dgmu2(:,:,1,2),l1,l2)
        yt6 = dpro(gmu1(:,:,1),dgmu2(:,:,2,2),l1,l2)
        yt7 = dpro(gmu1(:,:,1),dgmu2(:,:,3,2),l1,l2)
        yt8 = dpro(gmu1(:,:,2),dgmu2(:,:,2,2),l1,l2)
        yt9 = dpro(gmu1(:,:,2),dgmu2(:,:,3,2),l1,l2)
        yt0 = dpro(gmu1(:,:,3),dgmu2(:,:,3,2),l1,l2)
        zt1 = dpro(gmu1(:,:,0),dgmu2(:,:,0,3),l1,l2)
        zt2 = dpro(gmu1(:,:,0),dgmu2(:,:,1,3),l1,l2)
        zt3 = dpro(gmu1(:,:,0),dgmu2(:,:,2,3),l1,l2)
        zt4 = dpro(gmu1(:,:,0),dgmu2(:,:,3,3),l1,l2)
        zt5 = dpro(gmu1(:,:,1),dgmu2(:,:,1,3),l1,l2)
        zt6 = dpro(gmu1(:,:,1),dgmu2(:,:,2,3),l1,l2)
        zt7 = dpro(gmu1(:,:,1),dgmu2(:,:,3,3),l1,l2)
        zt8 = dpro(gmu1(:,:,2),dgmu2(:,:,2,3),l1,l2)
        zt9 = dpro(gmu1(:,:,2),dgmu2(:,:,3,3),l1,l2)
        zt0 = dpro(gmu1(:,:,3),dgmu2(:,:,3,3),l1,l2)
        xt1 = abs(xt1 - 0.25_dp*chi*(150.0_dp*c4-120.0_dp*c2+18.0_dp)*dchi2(1))
        xt2 = abs(xt2 + 0.375_dp*chi*(150.0_dp*c4-140.0_dp*c2+22.0_dp)*dchi2(1))
        xt3 = abs(xt3 - 3.75_dp*chi*(6.0_dp*c4-8.0_dp*c2+2.0_dp)*dchi2(1))
        xt4 = abs(xt4 + 0.625_dp*chi*(6.0_dp*c4-12.0_dp*c2+6.0_dp)*dchi2(1))
        xt5 = abs(xt5 - 0.0625_dp*chi*(1350.0_dp*c4-1220.0_dp*c2+222.0_dp)*dchi2(1))
        xt6 = abs(xt6 + 0.625_dp*chi*(54.0_dp*c4-44.0_dp*c2+6.0_dp)*dchi2(1))
        xt7 = abs(xt7 - 0.9375_dp*chi*(6.0_dp*c4-4.0_dp*c2-2.0_dp)*dchi2(1))
        xt8 = abs(xt8 - 0.25_dp*chi*(54.0_dp*c4+40.0_dp*c2-30.0_dp)*dchi2(1))
        xt9 = abs(xt9 + 0.375_dp*chi*(6.0_dp*c4+20.0_dp*c2-10.0_dp)*dchi2(1))
        xt0 = abs(xt0 - 0.0625_dp*chi*(6.0_dp*c4+60.0_dp*c2+30.0_dp)*dchi2(1))
        yt1 = abs(yt1 - 0.25_dp*chi*(150.0_dp*c4-120.0_dp*c2+18.0_dp)*dchi2(2))
        yt2 = abs(yt2 + 0.375_dp*chi*(150.0_dp*c4-140.0_dp*c2+22.0_dp)*dchi2(2))
        yt3 = abs(yt3 - 3.75_dp*chi*(6.0_dp*c4-8.0_dp*c2+2.0_dp)*dchi2(2))
        yt4 = abs(yt4 + 0.625_dp*chi*(6.0_dp*c4-12.0_dp*c2+6.0_dp)*dchi2(2))
        yt5 = abs(yt5 - 0.0625_dp*chi*(1350.0_dp*c4-1220.0_dp*c2+222.0_dp)*dchi2(2))
        yt6 = abs(yt6 + 0.625_dp*chi*(54.0_dp*c4-44.0_dp*c2+6.0_dp)*dchi2(2))
        yt7 = abs(yt7 - 0.9375_dp*chi*(6.0_dp*c4-4.0_dp*c2-2.0_dp)*dchi2(2))
        yt8 = abs(yt8 - 0.25_dp*chi*(54.0_dp*c4+40.0_dp*c2-30.0_dp)*dchi2(2))
        yt9 = abs(yt9 + 0.375_dp*chi*(6.0_dp*c4+20.0_dp*c2-10.0_dp)*dchi2(2))
        yt0 = abs(yt0 - 0.0625_dp*chi*(6.0_dp*c4+60.0_dp*c2+30.0_dp)*dchi2(2))
        zt1 = abs(zt1 - 0.25_dp*chi*(150.0_dp*c4-120.0_dp*c2+18.0_dp)*dchi2(3))
        zt2 = abs(zt2 + 0.375_dp*chi*(150.0_dp*c4-140.0_dp*c2+22.0_dp)*dchi2(3))
        zt3 = abs(zt3 - 3.75_dp*chi*(6.0_dp*c4-8.0_dp*c2+2.0_dp)*dchi2(3))
        zt4 = abs(zt4 + 0.625_dp*chi*(6.0_dp*c4-12.0_dp*c2+6.0_dp)*dchi2(3))
        zt5 = abs(zt5 - 0.0625_dp*chi*(1350.0_dp*c4-1220.0_dp*c2+222.0_dp)*dchi2(3))
        zt6 = abs(zt6 + 0.625_dp*chi*(54.0_dp*c4-44.0_dp*c2+6.0_dp)*dchi2(3))
        zt7 = abs(zt7 - 0.9375_dp*chi*(6.0_dp*c4-4.0_dp*c2-2.0_dp)*dchi2(3))
        zt8 = abs(zt8 - 0.25_dp*chi*(54.0_dp*c4+40.0_dp*c2-30.0_dp)*dchi2(3))
        zt9 = abs(zt9 + 0.375_dp*chi*(6.0_dp*c4+20.0_dp*c2-10.0_dp)*dchi2(3))
        zt0 = abs(zt0 - 0.0625_dp*chi*(6.0_dp*c4+60.0_dp*c2+30.0_dp)*dchi2(3))
        IF (xt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0x', l1, l2, xt1
        IF (xt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0x', l1, l2, xt2
        IF (xt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0x', l1, l2, xt3
        IF (xt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 0x', l1, l2, xt4
        IF (xt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 0x', l1, l2, xt5
        IF (xt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 0x', l1, l2, xt6
        IF (xt7>1.E-12) WRITE (*,*) ' Error in test for gmu #7 0x', l1, l2, xt7
        IF (xt8>1.E-12) WRITE (*,*) ' Error in test for gmu #8 0x', l1, l2, xt8
        IF (xt9>1.E-12) WRITE (*,*) ' Error in test for gmu #9 0x', l1, l2, xt9
        IF (xt0>1.E-12) WRITE (*,*) ' Error in test for gmu #0 0x', l1, l2, xt0
        emax = max(emax,xt1,xt2,xt3,xt4,xt5,xt6,xt7,xt8,xt9,xt0)
        IF (yt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0y', l1, l2, yt1
        IF (yt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0y', l1, l2, yt2
        IF (yt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0y', l1, l2, yt3
        IF (yt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 0y', l1, l2, yt4
        IF (yt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 0y', l1, l2, yt5
        IF (yt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 0y', l1, l2, yt6
        IF (yt7>1.E-12) WRITE (*,*) ' Error in test for gmu #7 0y', l1, l2, yt7
        IF (yt8>1.E-12) WRITE (*,*) ' Error in test for gmu #8 0y', l1, l2, yt8
        IF (yt9>1.E-12) WRITE (*,*) ' Error in test for gmu #9 0y', l1, l2, yt9
        IF (yt0>1.E-12) WRITE (*,*) ' Error in test for gmu #0 0y', l1, l2, yt0
        emax = max(emax,yt1,yt2,yt3,yt4,yt5,yt6,yt7,yt8,yt9,yt0)
        IF (zt1>1.E-12) WRITE (*,*) ' Error in test for gmu #1 0z', l1, l2, zt1
        IF (zt2>1.E-12) WRITE (*,*) ' Error in test for gmu #2 0z', l1, l2, zt2
        IF (zt3>1.E-12) WRITE (*,*) ' Error in test for gmu #3 0z', l1, l2, zt3
        IF (zt4>1.E-12) WRITE (*,*) ' Error in test for gmu #4 0z', l1, l2, zt4
        IF (zt5>1.E-12) WRITE (*,*) ' Error in test for gmu #5 0z', l1, l2, zt5
        IF (zt6>1.E-12) WRITE (*,*) ' Error in test for gmu #6 0z', l1, l2, zt6
        IF (zt7>1.E-12) WRITE (*,*) ' Error in test for gmu #7 0z', l1, l2, zt7
        IF (zt8>1.E-12) WRITE (*,*) ' Error in test for gmu #8 0z', l1, l2, zt8
        IF (zt9>1.E-12) WRITE (*,*) ' Error in test for gmu #9 0z', l1, l2, zt9
        IF (zt0>1.E-12) WRITE (*,*) ' Error in test for gmu #0 0z', l1, l2, zt0
        emax = max(emax,zt1,zt2,zt3,zt4,zt5,zt6,zt7,zt8,zt9,zt0)
        dx1 = dpro(dll1(:,:,0),ddll2(:,:,0,1),l1,l2)
        dx2 = dpro(dll1(:,:,0),ddll2(:,:,2,1),l1,l2)
        dx3 = dpro(dll1(:,:,0),ddll2(:,:,4,1),l1,l2)
        dx4 = dpro(dll1(:,:,0),ddll2(:,:,6,1),l1,l2)
        dx5 = dpro(dll1(:,:,2),ddll2(:,:,2,1),l1,l2)
        dx6 = dpro(dll1(:,:,2),ddll2(:,:,4,1),l1,l2)
        dx7 = dpro(dll1(:,:,2),ddll2(:,:,6,1),l1,l2)
        dx8 = dpro(dll1(:,:,4),ddll2(:,:,4,1),l1,l2)
        dx9 = dpro(dll1(:,:,4),ddll2(:,:,6,1),l1,l2)
        dx0 = dpro(dll1(:,:,6),ddll2(:,:,6,1),l1,l2)
        dy1 = dpro(dll1(:,:,0),ddll2(:,:,0,2),l1,l2)
        dy2 = dpro(dll1(:,:,0),ddll2(:,:,2,2),l1,l2)
        dy3 = dpro(dll1(:,:,0),ddll2(:,:,4,2),l1,l2)
        dy4 = dpro(dll1(:,:,0),ddll2(:,:,6,2),l1,l2)
        dy5 = dpro(dll1(:,:,2),ddll2(:,:,2,2),l1,l2)
        dy6 = dpro(dll1(:,:,2),ddll2(:,:,4,2),l1,l2)
        dy7 = dpro(dll1(:,:,2),ddll2(:,:,6,2),l1,l2)
        dy8 = dpro(dll1(:,:,4),ddll2(:,:,4,2),l1,l2)
        dy9 = dpro(dll1(:,:,4),ddll2(:,:,6,2),l1,l2)
        dy0 = dpro(dll1(:,:,6),ddll2(:,:,6,2),l1,l2)
        dz1 = dpro(dll1(:,:,0),ddll2(:,:,0,3),l1,l2)
        dz2 = dpro(dll1(:,:,0),ddll2(:,:,2,3),l1,l2)
        dz3 = dpro(dll1(:,:,0),ddll2(:,:,4,3),l1,l2)
        dz4 = dpro(dll1(:,:,0),ddll2(:,:,6,3),l1,l2)
        dz5 = dpro(dll1(:,:,2),ddll2(:,:,2,3),l1,l2)
        dz6 = dpro(dll1(:,:,2),ddll2(:,:,4,3),l1,l2)
        dz7 = dpro(dll1(:,:,2),ddll2(:,:,6,3),l1,l2)
        dz8 = dpro(dll1(:,:,4),ddll2(:,:,4,3),l1,l2)
        dz9 = dpro(dll1(:,:,4),ddll2(:,:,6,3),l1,l2)
        dz0 = dpro(dll1(:,:,6),ddll2(:,:,6,3),l1,l2)
        dx1 = abs(dx1)
        dx2 = abs(dx2)
        dx3 = abs(dx3)
        dx4 = abs(dx4)
        dx5 = abs(dx5-10.08_dp*chi*dchi2(1))
        dx6 = abs(dx6-1.44_dp*chi*dchi2(1))
        dx7 = abs(dx7-2.4_dp*chi*dchi2(1))
        dx8 = abs(dx8-0.04_dp*chi*(880.0_dp*c2-372.0_dp)*dchi2(1))
        dx9 = abs(dx9-0.6_dp*chi*(20.0_dp*c2-8.0_dp)*dchi2(1))
        dx0 = abs(dx0-0.25_dp*chi*(150.0_dp*c4-120.0_dp*c2+18.0_dp)*dchi2(1))
        dy1 = abs(dy1)
        dy2 = abs(dy2)
        dy3 = abs(dy3)
        dy4 = abs(dy4)
        dy5 = abs(dy5-10.08_dp*chi*dchi2(2))
        dy6 = abs(dy6-1.44_dp*chi*dchi2(2))
        dy7 = abs(dy7-2.4_dp*chi*dchi2(2))
        dy8 = abs(dy8-0.04_dp*chi*(880.0_dp*c2-372.0_dp)*dchi2(2))
        dy9 = abs(dy9-0.6_dp*chi*(20.0_dp*c2-8.0_dp)*dchi2(2))
        dy0 = abs(dy0-0.25_dp*chi*(150.0_dp*c4-120.0_dp*c2+18.0_dp)*dchi2(2))
        dz1 = abs(dz1)
        dz2 = abs(dz2)
        dz3 = abs(dz3)
        dz4 = abs(dz4)
        dz5 = abs(dz5-10.08_dp*chi*dchi2(3))
        dz6 = abs(dz6-1.44_dp*chi*dchi2(3))
        dz7 = abs(dz7-2.4_dp*chi*dchi2(3))
        dz8 = abs(dz8-0.04_dp*chi*(880.0_dp*c2-372.0_dp)*dchi2(3))
        dz9 = abs(dz9-0.6_dp*chi*(20.0_dp*c2-8.0_dp)*dchi2(3))
        dz0 = abs(dz0-0.25_dp*chi*(150.0_dp*c4-120.0_dp*c2+18.0_dp)*dchi2(3))
        IF (dx1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0x', l1, l2, dx1
        IF (dx2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0x', l1, l2, dx2
        IF (dx3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0x', l1, l2, dx3
        IF (dx4>1.E-12) WRITE (*,*) ' Error in test for dll #4 0x', l1, l2, dx4
        IF (dx5>1.E-12) WRITE (*,*) ' Error in test for dll #5 0x', l1, l2, dx5
        IF (dx6>1.E-12) WRITE (*,*) ' Error in test for dll #6 0x', l1, l2, dx6
        IF (dx7>1.E-12) WRITE (*,*) ' Error in test for dll #7 0x', l1, l2, dx7
        IF (dx8>1.E-12) WRITE (*,*) ' Error in test for dll #8 0x', l1, l2, dx8
        IF (dx9>1.E-12) WRITE (*,*) ' Error in test for dll #9 0x', l1, l2, dx9
        IF (dx0>1.E-12) WRITE (*,*) ' Error in test for dll #0 0x', l1, l2, dx0
        emax = max(emax,dx1,dx2,dx3,dx4,dx5,dx6,dx7,dx8,dx9,dx0)
        IF (dy1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0y', l1, l2, dy1
        IF (dy2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0y', l1, l2, dy2
        IF (dy3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0y', l1, l2, dy3
        IF (dy4>1.E-12) WRITE (*,*) ' Error in test for dll #4 0y', l1, l2, dy4
        IF (dy5>1.E-12) WRITE (*,*) ' Error in test for dll #5 0y', l1, l2, dy5
        IF (dy6>1.E-12) WRITE (*,*) ' Error in test for dll #6 0y', l1, l2, dy6
        IF (dy7>1.E-12) WRITE (*,*) ' Error in test for dll #7 0y', l1, l2, dy7
        IF (dy8>1.E-12) WRITE (*,*) ' Error in test for dll #8 0y', l1, l2, dy8
        IF (dy9>1.E-12) WRITE (*,*) ' Error in test for dll #9 0y', l1, l2, dy9
        IF (dy0>1.E-12) WRITE (*,*) ' Error in test for dll #0 0y', l1, l2, dy0
        emax = max(emax,dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8,dy9,dy0)
        IF (dz1>1.E-12) WRITE (*,*) ' Error in test for dll #1 0z', l1, l2, dz1
        IF (dz2>1.E-12) WRITE (*,*) ' Error in test for dll #2 0z', l1, l2, dz2
        IF (dz3>1.E-12) WRITE (*,*) ' Error in test for dll #3 0z', l1, l2, dz3
        IF (dz4>1.E-12) WRITE (*,*) ' Error in test for dll #4 0z', l1, l2, dz4
        IF (dz5>1.E-12) WRITE (*,*) ' Error in test for dll #5 0z', l1, l2, dz5
        IF (dz6>1.E-12) WRITE (*,*) ' Error in test for dll #6 0z', l1, l2, dz6
        IF (dz7>1.E-12) WRITE (*,*) ' Error in test for dll #7 0z', l1, l2, dz7
        IF (dz8>1.E-12) WRITE (*,*) ' Error in test for dll #8 0z', l1, l2, dz8
        IF (dz9>1.E-12) WRITE (*,*) ' Error in test for dll #9 0z', l1, l2, dz9
        IF (dz0>1.E-12) WRITE (*,*) ' Error in test for dll #0 0z', l1, l2, dz0
        emax = max(emax,dz1,dz2,dz3,dz4,dz5,dz6,dz7,dz8,dz9,dz0)
      END SUBROUTINE df_test
!------------------------------------------------------------------------------!
    END MODULE slater_koster_test
!------------------------------------------------------------------------------!

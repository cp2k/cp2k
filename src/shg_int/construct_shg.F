!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of the integrals over solid harmonic Gaussian(SHG) functions.
!>        Routines for (a|O(r12)|b) and overlap integrals (ab), (aba) and (abb).
!> \par Literature (partly)
!>      T.J. Giese and D. M. York, J. Chem. Phys, 128, 064104 (2008)
!>      T. Helgaker, P Joergensen, J. Olsen, Molecular Electronic-Structure
!>                                           Theory, Wiley
!> \par History
!>      created [04.2015]
!> \author Dorothea Golze
! **************************************************************************************************
MODULE construct_shg
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: dfac,&
                                              fac
   USE orbital_pointers,                ONLY: indso,&
                                              indso_inv,&
                                              nsoset
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'construct_shg'

! *** Public subroutines ***
   PUBLIC :: get_real_scaled_solid_harmonic, get_W_matrix, get_dW_matrix, &
             construct_int_shg_ab, construct_dev_shg_ab, construct_overlap_shg_aba, &
             dev_overlap_shg_aba, construct_overlap_shg_abb, dev_overlap_shg_abb

CONTAINS

! **************************************************************************************************
!> \brief computes the real scaled solid harmonics Rlm up to a given l
!> \param Rlm_c cosine part of real scaled soldi harmonics
!> \param Rlm_s sine part of real scaled soldi harmonics
!> \param l maximal l quantum up to where Rlm is calculated
!> \param r distance vector between a and b
!> \param r2 square of distance vector
! **************************************************************************************************
   SUBROUTINE get_real_scaled_solid_harmonic(Rlm_c, Rlm_s, l, r, r2)

      INTEGER, INTENT(IN)                                :: l
      REAL(KIND=dp), DIMENSION(0:l, -2*l:2*l), &
         INTENT(OUT)                                     :: Rlm_s, Rlm_c
      REAL(KIND=dp), DIMENSION(3)                        :: r
      REAL(KIND=dp)                                      :: r2

      CHARACTER(len=*), PARAMETER :: routineN = 'get_real_scaled_solid_harmonic', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, li, mi, prefac
      REAL(KIND=dp)                                      :: Rc, Rc_00, Rlm, Rmlm, Rplm, Rs, Rs_00, &
                                                            temp_c

      CALL timeset(routineN, handle)

      Rc_00 = 1.0_dp
      Rs_00 = 0.0_dp

      Rlm_c(0, 0) = Rc_00
      Rlm_s(0, 0) = Rs_00

      ! generate elements Rmm
      ! start
      IF (l > 0) THEN
         Rc = -0.5_dp*r(1)*Rc_00
         Rs = -0.5_dp*r(2)*Rc_00
         Rlm_c(1, 1) = Rc
         Rlm_s(1, 1) = Rs
         Rlm_c(1, -1) = -Rc
         Rlm_s(1, -1) = Rs
      ENDIF
      DO li = 2, l
         temp_c = (-r(1)*Rc+r(2)*Rs)/(REAL(2*(li-1)+2, dp))
         Rs = (-r(2)*Rc-r(1)*Rs)/(REAL(2*(li-1)+2, dp))
         Rc = temp_c
         Rlm_c(li, li) = Rc
         Rlm_s(li, li) = Rs
         IF (MODULO(li, 2) /= 0) THEN
            Rlm_c(li, -li) = -Rc
            Rlm_s(li, -li) = Rs
         ELSE
            Rlm_c(li, -li) = Rc
            Rlm_s(li, -li) = -Rs
         ENDIF
      ENDDO

      DO mi = 0, l-1
         Rmlm = Rlm_c(mi, mi)
         Rlm = r(3)*Rlm_c(mi, mi)
         Rlm_c(mi+1, mi) = Rlm
         IF (MODULO(mi, 2) /= 0) THEN
            Rlm_c(mi+1, -mi) = -Rlm
         ELSE
            Rlm_c(mi+1, -mi) = Rlm
         ENDIF
         DO li = mi+2, l
            prefac = (li+mi)*(li-mi)
            Rplm = (REAL(2*li-1, dp)*r(3)*Rlm-r2*Rmlm)/REAL(prefac, dp)
            Rmlm = Rlm
            Rlm = Rplm
            Rlm_c(li, mi) = Rlm
            IF (MODULO(mi, 2) /= 0) THEN
               Rlm_c(li, -mi) = -Rlm
            ELSE
               Rlm_c(li, -mi) = Rlm
            ENDIF
         ENDDO
      ENDDO
      DO mi = 1, l-1
         Rmlm = Rlm_s(mi, mi)
         Rlm = r(3)*Rlm_s(mi, mi)
         Rlm_s(mi+1, mi) = Rlm
         IF (MODULO(mi, 2) /= 0) THEN
            Rlm_s(mi+1, -mi) = Rlm
         ELSE
            Rlm_s(mi+1, -mi) = -Rlm
         ENDIF
         DO li = mi+2, l
            prefac = (li+mi)*(li-mi)
            Rplm = (REAL(2*li-1, dp)*r(3)*Rlm-r2*Rmlm)/REAL(prefac, dp)
            Rmlm = Rlm
            Rlm = Rplm
            Rlm_s(li, mi) = Rlm
            IF (MODULO(mi, 2) /= 0) THEN
               Rlm_s(li, -mi) = Rlm
            ELSE
               Rlm_s(li, -mi) = -Rlm
            ENDIF
         ENDDO
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE get_real_scaled_solid_harmonic

! **************************************************************************************************
!> \brief Calculate the prefactor A(l,m) = (-1)^m \sqrt[(2-delta(m,0))(l+m)!(l-m)!]
!> \param lmax maximal l quantum number
!> \param A matrix storing the prefactor for a given l and m
! **************************************************************************************************
   SUBROUTINE get_Alm(lmax, A)

      INTEGER, INTENT(IN)                                :: lmax
      REAL(KIND=dp), DIMENSION(0:lmax, 0:lmax), &
         INTENT(INOUT)                                   :: A

      CHARACTER(len=*), PARAMETER :: routineN = 'get_Alm', routineP = moduleN//':'//routineN

      INTEGER                                            :: l, m
      REAL(KIND=dp)                                      :: temp

      DO l = 0, lmax
         DO m = 0, l
            temp = SQRT(fac(l+m)*fac(l-m))
            IF (MODULO(m, 2) /= 0) temp = -temp
            IF (m /= 0) temp = temp*SQRT(2.0_dp)
            A(l, m) = temp
         ENDDO
      ENDDO

   END SUBROUTINE get_Alm

! **************************************************************************************************
!> \brief calculates the prefactors for the derivatives of the W matrix
!> \param lmax maximal l quantum number
!> \param dA_p = A(l,m)/A(l-1,m+1)
!> \param dA_m = A(l,m)/A(l-1,m-1)
!> \param dA = A(l,m)/A(l-1,m)
!> \note for m=0, W_l-1,-1 can't be read from Waux_mat, but we use
!>       W_l-1,-1 = -W_l-1,1 [cc(1), cs(2)] or W_l-1,-1 = W_l-1,1 [[sc(3), ss(4)], i.e.
!>       effectively we multiply dA_p by 2
! **************************************************************************************************
   SUBROUTINE get_dA_prefactors(lmax, dA_p, dA_m, dA)

      INTEGER, INTENT(IN)                                :: lmax
      REAL(KIND=dp), DIMENSION(0:lmax, 0:lmax), &
         INTENT(INOUT)                                   :: dA_p, dA_m, dA

      CHARACTER(len=*), PARAMETER :: routineN = 'get_dA_prefactors', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: l, m
      REAL(KIND=dp)                                      :: bm, bm_m, bm_p

      DO l = 0, lmax
         DO m = 0, l
            bm = 1.0_dp
            bm_m = 1.0_dp
            bm_p = 1.0_dp
            IF (m /= 0) bm = SQRT(2.0_dp)
            IF (m-1 /= 0) bm_m = SQRT(2.0_dp)
            IF (m+1 /= 0) bm_p = SQRT(2.0_dp)
            dA_p(l, m) = -bm/bm_p*SQRT(REAL((l-m)*(l-m-1), dp))
            dA_m(l, m) = -bm/bm_m*SQRT(REAL((l+m)*(l+m-1), dp))
            dA(l, m) = 2.0_dp*SQRT(REAL((l+m)*(l-m), dp))
            IF (m == 0) dA_p(l, m) = 2.0_dp*dA_p(l, m)
         ENDDO
      ENDDO
   END SUBROUTINE get_dA_prefactors

! **************************************************************************************************
!> \brief calculates the angular dependent-part of the SHG integrals,
!>        transformation matrix W, see literature above
!> \param lamax array of maximal l quantum number on a;
!>        lamax(lb) with lb= 0..lbmax
!> \param lbmax maximal l quantum number on b
!> \param lmax maximal l quantum number
!> \param Rc cosine part of real scaled solid harmonics
!> \param Rs sine part of real scaled solid harmonics
!> \param Waux_mat stores the angular-dependent part of the SHG integrals
! **************************************************************************************************
   SUBROUTINE get_W_matrix(lamax, lbmax, lmax, Rc, Rs, Waux_mat)

      INTEGER, DIMENSION(:), POINTER                     :: lamax
      INTEGER, INTENT(IN)                                :: lbmax, lmax
      REAL(KIND=dp), DIMENSION(0:lmax, -2*lmax:2*lmax), &
         INTENT(IN)                                      :: Rc, Rs
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: Waux_mat

      CHARACTER(len=*), PARAMETER :: routineN = 'get_W_matrix', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, j, k, la, labmin, laj, lb, lbj, &
                                                            ma, ma_m, ma_p, mb, mb_m, mb_p, nla, &
                                                            nlb
      REAL(KIND=dp) :: A_jk, A_lama, A_lbmb, Alm_fac, delta_k, prefac, Rca_m, Rca_p, Rcb_m, Rcb_p, &
         Rsa_m, Rsa_p, Rsb_m, Rsb_p, sign_fac, Wa(4), Wb(4), Wmat(4)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: A

      CALL timeset(routineN, handle)

      Wa(:) = 0.0_dp
      Wb(:) = 0.0_dp
      Wmat(:) = 0.0_dp

      ALLOCATE (A(0:lmax, 0:lmax))
      CALL get_Alm(lmax, A)

      DO lb = 0, lbmax
         nlb = nsoset(lb-1)
         DO la = 0, lamax(lb)
            nla = nsoset(la-1)
            labmin = MIN(la, lb)
            DO mb = 0, lb
               A_lbmb = A(lb, mb)
               IF (MODULO(lb, 2) /= 0) A_lbmb = -A_lbmb
               DO ma = 0, la
                  A_lama = A(la, ma)
                  Alm_fac = A_lama*A_lbmb
                  DO j = 0, labmin
                     laj = la-j
                     lbj = lb-j
                     prefac = Alm_fac*REAL(2**(la+lb-j), dp)*dfac(2*j-1)
                     delta_k = 0.5_dp
                     Wmat = 0.0_dp
                     DO k = 0, j
                        ma_m = ma-k
                        ma_p = ma+k
                        IF (laj < ABS(ma_m) .AND. laj < ABS(ma_p)) CYCLE
                        mb_m = mb-k
                        mb_p = mb+k
                        IF (lbj < ABS(mb_m) .AND. lbj < ABS(mb_p)) CYCLE
                        IF (k /= 0) delta_k = 1.0_dp
                        A_jk = fac(j+k)*fac(j-k)
                        IF (k /= 0) A_jk = 2.0_dp*A_jk
                        IF (MODULO(k, 2) /= 0) THEN
                           sign_fac = -1.0_dp
                        ELSE
                           sign_fac = 1.0_dp
                        ENDIF
                        Rca_m = Rc(laj, ma_m)
                        Rsa_m = Rs(laj, ma_m)
                        Rca_p = Rc(laj, ma_p)
                        Rsa_p = Rs(laj, ma_p)
                        Rcb_m = Rc(lbj, mb_m)
                        Rsb_m = Rs(lbj, mb_m)
                        Rcb_p = Rc(lbj, mb_p)
                        Rsb_p = Rs(lbj, mb_p)
                        Wa(1) = delta_k*(Rca_m+sign_fac*Rca_p)
                        Wb(1) = delta_k*(Rcb_m+sign_fac*Rcb_p)
                        Wa(2) = -Rsa_m+sign_fac*Rsa_p
                        Wb(2) = -Rsb_m+sign_fac*Rsb_p
                        Wmat(1) = Wmat(1)+prefac/A_jk*(Wa(1)*Wb(1)+Wa(2)*Wb(2))
                        IF (mb > 0) THEN
                           Wb(3) = delta_k*(Rsb_m+sign_fac*Rsb_p)
                           Wb(4) = Rcb_m-sign_fac*Rcb_p
                           Wmat(2) = Wmat(2)+prefac/A_jk*(Wa(1)*Wb(3)+Wa(2)*Wb(4))
                        ENDIF
                        IF (ma > 0) THEN
                           Wa(3) = delta_k*(Rsa_m+sign_fac*Rsa_p)
                           Wa(4) = Rca_m-sign_fac*Rca_p
                           Wmat(3) = Wmat(3)+prefac/A_jk*(Wa(3)*Wb(1)+Wa(4)*Wb(2))
                        ENDIF
                        IF (ma > 0 .AND. mb > 0) THEN
                           Wmat(4) = Wmat(4)+prefac/A_jk*(Wa(3)*Wb(3)+Wa(4)*Wb(4))
                        ENDIF
                     ENDDO
                     Waux_mat(j+1, nla+la+1+ma, nlb+lb+1+mb) = Wmat(1)
                     IF (mb > 0) Waux_mat(j+1, nla+la+1+ma, nlb+lb+1-mb) = Wmat(2)
                     IF (ma > 0) Waux_mat(j+1, nla+la+1-ma, nlb+lb+1+mb) = Wmat(3)
                     IF (ma > 0 .AND. mb > 0) Waux_mat(j+1, nla+la+1-ma, nlb+lb+1-mb) = Wmat(4)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE (A)

      CALL timestop(handle)

   END SUBROUTINE get_W_matrix

! **************************************************************************************************
!> \brief calculates derivatives of transformation matrix W,
!> \param lamax array of maximal l quantum number on a;
!>        lamax(lb) with lb= 0..lbmax
!> \param lbmax maximal l quantum number on b
!> \param Waux_mat stores the angular-dependent part of the SHG integrals
!> \param dWaux_mat stores the derivatives of the angular-dependent part of
!>        the SHG integrals
! **************************************************************************************************
   SUBROUTINE get_dW_matrix(lamax, lbmax, Waux_mat, dWaux_mat)

      INTEGER, DIMENSION(:), POINTER                     :: lamax
      INTEGER, INTENT(IN)                                :: lbmax
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: Waux_mat
      REAL(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(INOUT)                                   :: dWaux_mat

      CHARACTER(len=*), PARAMETER :: routineN = 'get_dW_matrix', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ima, imam, imb, imbm, ipa, ipam, &
                                                            ipb, ipbm, j, jmax, la, labm, labmin, &
                                                            lamb, lb, lmax, ma, mb, nla, nlam, &
                                                            nlb, nlbm
      REAL(KIND=dp)                                      :: dAa, dAa_m, dAa_p, dAb, dAb_m, dAb_p
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dA, dA_m, dA_p, Wam, Wamm, Wamp, Wbm, &
                                                            Wbmm, Wbmp

      CALL timeset(routineN, handle)

      jmax = MIN(MAXVAL(lamax), lbmax)
      ALLOCATE (Wam(0:jmax, 4), Wamm(0:jmax, 4), Wamp(0:jmax, 4))
      ALLOCATE (Wbm(0:jmax, 4), Wbmm(0:jmax, 4), Wbmp(0:jmax, 4))

      !*** get dA_p=A(l,m)/A(l-1,m+1)
      !*** get dA_m=A(l,m)/A(l-1,m-1)
      !*** get dA=2*A(l,m)/A(l-1,m)
      lmax = MAX(MAXVAL(lamax), lbmax)
      ALLOCATE (dA_p(0:lmax, 0:lmax), dA_m(0:lmax, 0:lmax), dA(0:lmax, 0:lmax))
      CALL get_dA_prefactors(lmax, dA_p, dA_m, dA)

      DO lb = 0, lbmax
         nlb = nsoset(lb-1)
         nlbm = 0
         IF (lb > 0) nlbm = nsoset(lb-2)
         DO la = 0, lamax(lb)
            nla = nsoset(la-1)
            nlam = 0
            IF (la > 0) nlam = nsoset(la-2)
            labmin = MIN(la, lb)
            lamb = MIN(la-1, lb)
            labm = MIN(la, lb-1)
            DO mb = 0, lb
               dAb = dA(lb, mb)
               dAb_p = dA_p(lb, mb)
               dAb_m = dA_m(lb, mb)
               ipb = nlb+lb+mb+1
               imb = nlb+lb-mb+1
               ipbm = nlbm+lb+mb
               imbm = nlbm+lb-mb
               DO ma = 0, la
                  dAa = dA(la, ma)
                  dAa_p = dA_p(la, ma)
                  dAa_m = dA_m(la, ma)
                  ipa = nla+la+ma+1
                  ima = nla+la-ma+1
                  ipam = nlam+la+ma
                  imam = nlam+la-ma
                  Wam(:, :) = 0.0_dp
                  Wamm(:, :) = 0.0_dp
                  Wamp(:, :) = 0.0_dp
                  !*** Wam: la-1, ma
                  IF (ma <= la-1) THEN
                     Wam(0:lamb, 1) = Waux_mat(1:lamb+1, ipam, ipb)
                     IF (mb > 0) Wam(0:lamb, 2) = Waux_mat(1:lamb+1, ipam, imb)
                     IF (ma > 0) Wam(0:lamb, 3) = Waux_mat(1:lamb+1, imam, ipb)
                     IF (ma > 0 .AND. mb > 0) Wam(0:lamb, 4) = Waux_mat(1:lamb+1, imam, imb)
                  ENDIF
                  !*** Wamm: la-1, ma-1
                  IF (ma-1 >= 0) THEN
                     Wamm(0:lamb, 1) = Waux_mat(1:lamb+1, ipam-1, ipb)
                     IF (mb > 0) Wamm(0:lamb, 2) = Waux_mat(1:lamb+1, ipam-1, imb)
                     IF (ma-1 > 0) Wamm(0:lamb, 3) = Waux_mat(1:lamb+1, imam+1, ipb) !order: e.g. -1 0 1, if < 0 |m|, -1 means  -m+1
                     IF (ma-1 > 0 .AND. mb > 0) Wamm(0:lamb, 4) = Waux_mat(1:lamb+1, imam+1, imb)
                  ENDIF
                  !*** Wamp: la-1, ma+1
                  IF (ma+1 <= la-1) THEN
                     Wamp(0:lamb, 1) = Waux_mat(1:lamb+1, ipam+1, ipb)
                     IF (mb > 0) Wamp(0:lamb, 2) = Waux_mat(1:lamb+1, ipam+1, imb)
                     IF (ma+1 > 0) Wamp(0:lamb, 3) = Waux_mat(1:lamb+1, imam-1, ipb)
                     IF (ma+1 > 0 .AND. mb > 0) Wamp(0:lamb, 4) = Waux_mat(1:lamb+1, imam-1, imb)
                  ENDIF
                  Wbm(:, :) = 0.0_dp
                  Wbmm(:, :) = 0.0_dp
                  Wbmp(:, :) = 0.0_dp
                  !*** Wbm: lb-1, mb
                  IF (mb <= lb-1) THEN
                     Wbm(0:labm, 1) = Waux_mat(1:labm+1, ipa, ipbm)
                     IF (mb > 0) Wbm(0:labm, 2) = Waux_mat(1:labm+1, ipa, imbm)
                     IF (ma > 0) Wbm(0:labm, 3) = Waux_mat(1:labm+1, ima, ipbm)
                     IF (ma > 0 .AND. mb > 0) Wbm(0:labm, 4) = Waux_mat(1:labm+1, ima, imbm)
                  ENDIF
                  !*** Wbmm: lb-1, mb-1
                  IF (mb-1 >= 0) THEN
                     Wbmm(0:labm, 1) = Waux_mat(1:labm+1, ipa, ipbm-1)
                     IF (mb-1 > 0) Wbmm(0:labm, 2) = Waux_mat(1:labm+1, ipa, imbm+1)
                     IF (ma > 0) Wbmm(0:labm, 3) = Waux_mat(1:labm+1, ima, ipbm-1)
                     IF (ma > 0 .AND. mb-1 > 0) Wbmm(0:labm, 4) = Waux_mat(1:labm+1, ima, imbm+1)
                  ENDIF
                  !*** Wbmp: lb-1, mb+1
                  IF (mb+1 <= lb-1) THEN
                     Wbmp(0:labm, 1) = Waux_mat(1:labm+1, ipa, ipbm+1)
                     IF (mb+1 > 0) Wbmp(0:labm, 2) = Waux_mat(1:labm+1, ipa, imbm-1)
                     IF (ma > 0) Wbmp(0:labm, 3) = Waux_mat(1:labm+1, ima, ipbm+1)
                     IF (ma > 0 .AND. mb+1 > 0) Wbmp(0:labm, 4) = Waux_mat(1:labm+1, ima, imbm-1)
                  ENDIF
                  DO j = 0, labmin
                     !*** x component
                     dWaux_mat(1, j+1, ipa, ipb) = dAa_p*Wamp(j, 1)-dAa_m*Wamm(j, 1) &
                                                   -dAb_p*Wbmp(j, 1)+dAb_m*Wbmm(j, 1)
                     IF (mb > 0) THEN
                        dWaux_mat(1, j+1, ipa, imb) = dAa_p*Wamp(j, 2)-dAa_m*Wamm(j, 2) &
                                                      -dAb_p*Wbmp(j, 2)+dAb_m*Wbmm(j, 2)
                     ENDIF
                     IF (ma > 0) THEN
                        dWaux_mat(1, j+1, ima, ipb) = dAa_p*Wamp(j, 3)-dAa_m*Wamm(j, 3) &
                                                      -dAb_p*Wbmp(j, 3)+dAb_m*Wbmm(j, 3)
                     ENDIF
                     IF (ma > 0 .AND. mb > 0) THEN
                        dWaux_mat(1, j+1, ima, imb) = dAa_p*Wamp(j, 4)-dAa_m*Wamm(j, 4) &
                                                      -dAb_p*Wbmp(j, 4)+dAb_m*Wbmm(j, 4)
                     ENDIF

                     !**** y component
                     dWaux_mat(2, j+1, ipa, ipb) = dAa_p*Wamp(j, 3)+dAa_m*Wamm(j, 3) &
                                                   -dAb_p*Wbmp(j, 2)-dAb_m*Wbmm(j, 2)
                     IF (mb > 0) THEN
                        dWaux_mat(2, j+1, ipa, imb) = dAa_p*Wamp(j, 4)+dAa_m*Wamm(j, 4) &
                                                      +dAb_p*Wbmp(j, 1)+dAb_m*Wbmm(j, 1)
                     ENDIF
                     IF (ma > 0) THEN
                        dWaux_mat(2, j+1, ima, ipb) = -dAa_p*Wamp(j, 1)-dAa_m*Wamm(j, 1) &
                                                      -dAb_p*Wbmp(j, 4)-dAb_m*Wbmm(j, 4)
                     ENDIF
                     IF (ma > 0 .AND. mb > 0) THEN
                        dWaux_mat(2, j+1, ima, imb) = -dAa_p*Wamp(j, 2)-dAa_m*Wamm(j, 2) &
                                                      +dAb_p*Wbmp(j, 3)+dAb_m*Wbmm(j, 3)
                     ENDIF
                     !**** z compnent
                     dWaux_mat(3, j+1, ipa, ipb) = dAa*Wam(j, 1)-dAb*Wbm(j, 1)
                     IF (mb > 0) THEN
                        dWaux_mat(3, j+1, ipa, imb) = dAa*Wam(j, 2)-dAb*Wbm(j, 2)
                     ENDIF
                     IF (ma > 0) THEN
                        dWaux_mat(3, j+1, ima, ipb) = dAa*Wam(j, 3)-dAb*Wbm(j, 3)
                     ENDIF
                     IF (ma > 0 .AND. mb > 0) THEN
                        dWaux_mat(3, j+1, ima, imb) = dAa*Wam(j, 4)-dAb*Wbm(j, 4)
                     ENDIF

                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE (Wam, Wamm, Wamp)
      DEALLOCATE (Wbm, Wbmm, Wbmp)
      DEALLOCATE (dA, dA_p, dA_m)

      CALL timestop(handle)

   END SUBROUTINE get_dW_matrix

! **************************************************************************************************
!> \brief calculates [ab] SHG overlap integrals using precomputed angular-
!>        dependent part
!> \param la set of l quantum number on a
!> \param first_sgfa indexing
!> \param nshella number of shells for a
!> \param lb set of l quantum number on b
!> \param first_sgfb indexing
!> \param nshellb number of shells for b
!> \param swork_cont contracted and normalized [s|s] integrals
!> \param Waux_mat precomputed angular-dependent part
!> \param sab contracted integral of spherical harmonic Gaussianslm
! **************************************************************************************************
   SUBROUTINE construct_int_shg_ab(la, first_sgfa, nshella, lb, first_sgfb, nshellb, &
                                   swork_cont, Waux_mat, sab)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: la, first_sgfa
      INTEGER, INTENT(IN)                                :: nshella
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lb, first_sgfb
      INTEGER, INTENT(IN)                                :: nshellb
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: swork_cont, Waux_mat
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: sab

      CHARACTER(len=*), PARAMETER :: routineN = 'construct_int_shg_ab', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: fnla, fnlb, fsgfa, fsgfb, handle, &
                                                            ishella, j, jshellb, labmin, lai, lbj, &
                                                            lnla, lnlb, lsgfa, lsgfb, mai, mbj
      REAL(KIND=dp)                                      :: prefac

      CALL timeset(routineN, handle)

      DO jshellb = 1, nshellb
         lbj = lb(jshellb)
         fnlb = nsoset(lbj-1)+1
         lnlb = nsoset(lbj)
         fsgfb = first_sgfb(jshellb)
         lsgfb = fsgfb+2*lbj
         DO ishella = 1, nshella
            lai = la(ishella)
            fnla = nsoset(lai-1)+1
            lnla = nsoset(lai)
            fsgfa = first_sgfa(ishella)
            lsgfa = fsgfa+2*lai
            labmin = MIN(lai, lbj)
            DO mbj = 0, 2*lbj
               DO mai = 0, 2*lai
                  DO j = 0, labmin
                     prefac = swork_cont(lai+lbj-j+1, ishella, jshellb)
                     sab(fsgfa+mai, fsgfb+mbj) = sab(fsgfa+mai, fsgfb+mbj) &
                                                 +prefac*Waux_mat(j+1, fnla+mai, fnlb+mbj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE construct_int_shg_ab

! **************************************************************************************************
!> \brief calculates derivatives of [ab] SHG overlap integrals using precomputed
!>        angular-dependent part
!> \param la set of l quantum number on a
!> \param first_sgfa indexing
!> \param nshella number of shells for a
!> \param lb set of l quantum number on b
!> \param first_sgfb indexing
!> \param nshellb number of shells for b
!> \param rab distance vector Ra-Rb
!> \param swork_cont contracted and normalized [s|s] integrals
!> \param Waux_mat precomputed angular-dependent part
!> \param dWaux_mat ...
!> \param dsab derivative of contracted integral of spherical harmonic Gaussians
! **************************************************************************************************
   SUBROUTINE construct_dev_shg_ab(la, first_sgfa, nshella, lb, first_sgfb, nshellb, rab, &
                                   swork_cont, Waux_mat, dWaux_mat, dsab)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: la, first_sgfa
      INTEGER, INTENT(IN)                                :: nshella
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lb, first_sgfb
      INTEGER, INTENT(IN)                                :: nshellb
      REAL(KIND=dp), INTENT(IN)                          :: rab(3)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: swork_cont, Waux_mat
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: dWaux_mat
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: dsab

      CHARACTER(len=*), PARAMETER :: routineN = 'construct_dev_shg_ab', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: fnla, fnlb, fsgfa, fsgfb, handle, i, &
                                                            ishella, j, jshellb, labmin, lai, lbj, &
                                                            lnla, lnlb, lsgfa, lsgfb
      REAL(KIND=dp)                                      :: dprefac, prefac, rabx2(3)

      CALL timeset(routineN, handle)

      rabx2(:) = 2.0_dp*rab
      DO jshellb = 1, nshellb
         lbj = lb(jshellb)
         fnlb = nsoset(lbj-1)+1
         lnlb = nsoset(lbj)
         fsgfb = first_sgfb(jshellb)
         lsgfb = fsgfb+2*lbj
         DO ishella = 1, nshella
            lai = la(ishella)
            fnla = nsoset(lai-1)+1
            lnla = nsoset(lai)
            fsgfa = first_sgfa(ishella)
            lsgfa = fsgfa+2*lai
            labmin = MIN(lai, lbj)
            DO j = 0, labmin
               prefac = swork_cont(lai+lbj-j+1, ishella, jshellb)
               dprefac = swork_cont(lai+lbj-j+2, ishella, jshellb) !j+1
               DO i = 1, 3
                  dsab(fsgfa:lsgfa, fsgfb:lsgfb, i) = dsab(fsgfa:lsgfa, fsgfb:lsgfb, i) &
                                                      +rabx2(i)*dprefac*Waux_mat(j+1, fnla:lnla, fnlb:lnlb) &
                                                      +prefac*dWaux_mat(i, j+1, fnla:lnla, fnlb:lnlb)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE construct_dev_shg_ab

! **************************************************************************************************
!> \brief calculates [aba] SHG overlap integrals using precomputed angular-
!>        dependent part
!> \param la set of l quantum number on a, orbital basis
!> \param first_sgfa indexing
!> \param nshella number of shells for a, orbital basis
!> \param lb set of l quantum number on b. orbital basis
!> \param first_sgfb indexing
!> \param nshellb number of shells for b, orbital basis
!> \param lca of l quantum number on a, aux basis
!> \param first_sgfca indexing
!> \param nshellca  number of shells for a, aux basis
!> \param cg_coeff Clebsch-Gordon coefficients
!> \param cg_none0_list list of none-zero Clebsch-Gordon coefficients
!> \param ncg_none0 number of non-zero Clebsch-Gordon coefficients
!> \param swork_cont contracted and normalized [s|ra^n|s] integrals
!> \param Waux_mat precomputed angular-dependent part
!> \param saba contracted overlap [aba] of spherical harmonic Gaussians
! **************************************************************************************************
   SUBROUTINE construct_overlap_shg_aba(la, first_sgfa, nshella, lb, first_sgfb, nshellb, &
                                        lca, first_sgfca, nshellca, cg_coeff, cg_none0_list, &
                                        ncg_none0, swork_cont, Waux_mat, saba)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: la, first_sgfa
      INTEGER, INTENT(IN)                                :: nshella
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lb, first_sgfb
      INTEGER, INTENT(IN)                                :: nshellb
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lca, first_sgfca
      INTEGER, INTENT(IN)                                :: nshellca
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: cg_coeff
      INTEGER, DIMENSION(:, :, :), INTENT(IN)            :: cg_none0_list
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: ncg_none0
      REAL(KIND=dp), DIMENSION(:, 0:, :, :, :), &
         INTENT(IN)                                      :: swork_cont
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Waux_mat
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: saba

      CHARACTER(len=*), PARAMETER :: routineN = 'construct_overlap_shg_aba', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, ia, il, ilist, ishella, isoa1, isoa2, isoaa, j, jb, jshellb, ka, kshella, &
         laa, labmin, lai, lak, lbj, maa, mai, mak, mbj, nla, nlb, sgfa, sgfb, sgfca
      REAL(KIND=dp)                                      :: prefac, stemp

      CALL timeset(routineN, handle)

      DO kshella = 1, nshellca
         lak = lca(kshella)
         sgfca = first_sgfca(kshella)
         ka = sgfca+lak
         DO jshellb = 1, nshellb
            lbj = lb(jshellb)
            nlb = nsoset(lbj-1)+lbj+1
            sgfb = first_sgfb(jshellb)
            jb = sgfb+lbj
            DO ishella = 1, nshella
               lai = la(ishella)
               sgfa = first_sgfa(ishella)
               ia = sgfa+lai
               DO mai = -lai, lai, 1
                  DO mak = -lak, lak, 1
                     isoa1 = indso_inv(lai, mai)
                     isoa2 = indso_inv(lak, mak)
                     DO mbj = -lbj, lbj, 1
                        DO ilist = 1, ncg_none0(isoa1, isoa2)
                           isoaa = cg_none0_list(isoa1, isoa2, ilist)
                           laa = indso(1, isoaa)
                           maa = indso(2, isoaa)
                           nla = nsoset(laa-1)+laa+1
                           labmin = MIN(laa, lbj)
                           il = INT((lai+lak-laa)/2)
                           stemp = 0.0_dp
                           DO j = 0, labmin
                              prefac = swork_cont(laa+lbj-j+1, il, ishella, jshellb, kshella)
                              stemp = stemp+prefac*Waux_mat(j+1, nla+maa, nlb+mbj)
                           ENDDO
                           saba(ia+mai, jb+mbj, ka+mak) = saba(ia+mai, jb+mbj, ka+mak)+cg_coeff(isoa1, isoa2, isoaa)*stemp
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE construct_overlap_shg_aba

! **************************************************************************************************
!> \brief calculates derivatives of [aba] SHG overlap integrals using
!>        precomputed angular-dependent part
!> \param la set of l quantum number on a, orbital basis
!> \param first_sgfa indexing
!> \param nshella number of shells for a, orbital basis
!> \param lb set of l quantum number on b. orbital basis
!> \param first_sgfb indexing
!> \param nshellb number of shells for b, orbital basis
!> \param lca of l quantum number on a, aux basis
!> \param first_sgfca indexing
!> \param nshellca  number of shells for a, aux basis
!> \param cg_coeff Clebsch-Gordon coefficients
!> \param cg_none0_list list of none-zero Clebsch-Gordon coefficients
!> \param ncg_none0 number of non-zero Clebsch-Gordon coefficients
!> \param rab distance vector Ra-Rb
!> \param swork_cont contracted and normalized [s|ra^n|s] integrals
!> \param Waux_mat precomputed angular-dependent part
!> \param dWaux_mat derivatives of precomputed angular-dependent part
!> \param dsaba derivative of contracted overlap [aba] of spherical harmonic
!>              Gaussians
! **************************************************************************************************
   SUBROUTINE dev_overlap_shg_aba(la, first_sgfa, nshella, lb, first_sgfb, nshellb, &
                                  lca, first_sgfca, nshellca, cg_coeff, cg_none0_list, &
                                  ncg_none0, rab, swork_cont, Waux_mat, dWaux_mat, dsaba)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: la, first_sgfa
      INTEGER, INTENT(IN)                                :: nshella
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lb, first_sgfb
      INTEGER, INTENT(IN)                                :: nshellb
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lca, first_sgfca
      INTEGER, INTENT(IN)                                :: nshellca
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: cg_coeff
      INTEGER, DIMENSION(:, :, :), INTENT(IN)            :: cg_none0_list
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: ncg_none0
      REAL(KIND=dp), INTENT(IN)                          :: rab(3)
      REAL(KIND=dp), DIMENSION(:, 0:, :, :, :), &
         INTENT(IN)                                      :: swork_cont
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Waux_mat
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: dWaux_mat
      REAL(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(INOUT)                                   :: dsaba

      CHARACTER(len=*), PARAMETER :: routineN = 'dev_overlap_shg_aba', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, i, ia, il, ilist, ishella, isoa1, isoa2, isoaa, j, jb, jshellb, ka, &
         kshella, laa, labmin, lai, lak, lbj, maa, mai, mak, mbj, nla, nlb, sgfa, sgfb, sgfca
      REAL(KIND=dp)                                      :: dprefac, dtemp(3), prefac, rabx2(3)

      CALL timeset(routineN, handle)

      rabx2(:) = 2.0_dp*rab

      DO kshella = 1, nshellca
         lak = lca(kshella)
         sgfca = first_sgfca(kshella)
         ka = sgfca+lak
         DO jshellb = 1, nshellb
            lbj = lb(jshellb)
            nlb = nsoset(lbj-1)+lbj+1
            sgfb = first_sgfb(jshellb)
            jb = sgfb+lbj
            DO ishella = 1, nshella
               lai = la(ishella)
               sgfa = first_sgfa(ishella)
               ia = sgfa+lai
               DO mai = -lai, lai, 1
                  DO mak = -lak, lak, 1
                     isoa1 = indso_inv(lai, mai)
                     isoa2 = indso_inv(lak, mak)
                     DO mbj = -lbj, lbj, 1
                        DO ilist = 1, ncg_none0(isoa1, isoa2)
                           isoaa = cg_none0_list(isoa1, isoa2, ilist)
                           laa = indso(1, isoaa)
                           maa = indso(2, isoaa)
                           nla = nsoset(laa-1)+laa+1
                           labmin = MIN(laa, lbj)
                           il = (lai+lak-laa)/2 ! lai+lak-laa always even
                           dtemp = 0.0_dp
                           DO j = 0, labmin
                              prefac = swork_cont(laa+lbj-j+1, il, ishella, jshellb, kshella)
                              dprefac = swork_cont(laa+lbj-j+2, il, ishella, jshellb, kshella)
                              DO i = 1, 3
                                 dtemp(i) = dtemp(i)+rabx2(i)*dprefac*Waux_mat(j+1, nla+maa, nlb+mbj) &
                                            +prefac*dWaux_mat(i, j+1, nla+maa, nlb+mbj)
                              ENDDO
                           ENDDO
                           DO i = 1, 3
                              dsaba(ia+mai, jb+mbj, ka+mak, i) = dsaba(ia+mai, jb+mbj, ka+mak, i) &
                                                                 +cg_coeff(isoa1, isoa2, isoaa)*dtemp(i)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE dev_overlap_shg_aba

! **************************************************************************************************
!> \brief calculates [abb] SHG overlap integrals using precomputed angular-
!>        dependent part
!> \param la set of l quantum number on a, orbital basis
!> \param first_sgfa indexing
!> \param nshella number of shells for a, orbital basis
!> \param lb set of l quantum number on b. orbital basis
!> \param first_sgfb indexing
!> \param nshellb number of shells for b, orbital basis
!> \param lcb l quantum number on b, aux basis
!> \param first_sgfcb indexing
!> \param nshellcb number of shells for b, aux basis
!> \param cg_coeff Clebsch-Gordon coefficients
!> \param cg_none0_list list of none-zero Clebsch-Gordon coefficients
!> \param ncg_none0 number of non-zero Clebsch-Gordon coefficients
!> \param swork_cont contracted and normalized [s|rb^n|s] integrals
!> \param Waux_mat precomputed angular-dependent part
!> \param sabb contracted overlap [abb] of spherical harmonic Gaussians
! **************************************************************************************************
   SUBROUTINE construct_overlap_shg_abb(la, first_sgfa, nshella, lb, first_sgfb, nshellb, &
                                        lcb, first_sgfcb, nshellcb, cg_coeff, cg_none0_list, &
                                        ncg_none0, swork_cont, Waux_mat, sabb)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: la, first_sgfa
      INTEGER, INTENT(IN)                                :: nshella
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lb, first_sgfb
      INTEGER, INTENT(IN)                                :: nshellb
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lcb, first_sgfcb
      INTEGER, INTENT(IN)                                :: nshellcb
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: cg_coeff
      INTEGER, DIMENSION(:, :, :), INTENT(IN)            :: cg_none0_list
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: ncg_none0
      REAL(KIND=dp), DIMENSION(:, 0:, :, :, :), &
         INTENT(IN)                                      :: swork_cont
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Waux_mat
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: sabb

      CHARACTER(len=*), PARAMETER :: routineN = 'construct_overlap_shg_abb', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, ia, il, ilist, ishella, isob1, isob2, isobb, j, jb, jshellb, kb, kshellb, &
         labmin, lai, lbb, lbj, lbk, mai, mbb, mbj, mbk, nla, nlb, sgfa, sgfb, sgfcb
      REAL(KIND=dp)                                      :: prefac, stemp, tsign

      CALL timeset(routineN, handle)

      DO kshellb = 1, nshellcb
         lbk = lcb(kshellb)
         sgfcb = first_sgfcb(kshellb)
         kb = sgfcb+lbk
         DO jshellb = 1, nshellb
            lbj = lb(jshellb)
            sgfb = first_sgfb(jshellb)
            jb = sgfb+lbj
            DO ishella = 1, nshella
               lai = la(ishella)
               nla = nsoset(lai-1)+lai+1
               sgfa = first_sgfa(ishella)
               ia = sgfa+lai
               DO mbj = -lbj, lbj, 1
                  DO mbk = -lbk, lbk, 1
                     isob1 = indso_inv(lbj, mbj)
                     isob2 = indso_inv(lbk, mbk)
                     DO mai = -lai, lai, 1
                        DO ilist = 1, ncg_none0(isob1, isob2)
                           isobb = cg_none0_list(isob1, isob2, ilist)
                           lbb = indso(1, isobb)
                           mbb = indso(2, isobb)
                           nlb = nsoset(lbb-1)+lbb+1
                           ! tsgin: because we take the transpose of auxmat (calculated for (la,lb))
                           tsign = 1.0_dp
                           IF (MODULO(lbb-lai, 2) /= 0) tsign = -1.0_dp
                           labmin = MIN(lai, lbb)
                           il = INT((lbj+lbk-lbb)/2)
                           stemp = 0.0_dp
                           DO j = 0, labmin
                              prefac = swork_cont(lai+lbb-j+1, il, ishella, jshellb, kshellb)
                              stemp = stemp+prefac*Waux_mat(j+1, nlb+mbb, nla+mai)
                           ENDDO
                           sabb(ia+mai, jb+mbj, kb+mbk) = sabb(ia+mai, jb+mbj, kb+mbk)+tsign*cg_coeff(isob1, isob2, isobb)*stemp
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE construct_overlap_shg_abb

! **************************************************************************************************
!> \brief calculates derivatives of [abb] SHG overlap integrals using
!>        precomputed angular-dependent part
!> \param la set of l quantum number on a, orbital basis
!> \param first_sgfa indexing
!> \param nshella number of shells for a, orbital basis
!> \param lb set of l quantum number on b. orbital basis
!> \param first_sgfb indexing
!> \param nshellb number of shells for b, orbital basis
!> \param lcb l quantum number on b, aux basis
!> \param first_sgfcb indexing
!> \param nshellcb number of shells for b, aux basis
!> \param cg_coeff Clebsch-Gordon coefficients
!> \param cg_none0_list list of none-zero Clebsch-Gordon coefficients
!> \param ncg_none0 number of non-zero Clebsch-Gordon coefficients
!> \param rab distance vector Ra-Rb
!> \param swork_cont contracted and normalized [s|rb^n|s] integrals
!> \param Waux_mat precomputed angular-dependent part
!> \param dWaux_mat derivatives of precomputed angular-dependent part
!> \param dsabb derivative of contracted overlap [abb] of spherical harmonic
!>        Gaussians
! **************************************************************************************************
   SUBROUTINE dev_overlap_shg_abb(la, first_sgfa, nshella, lb, first_sgfb, nshellb, &
                                  lcb, first_sgfcb, nshellcb, cg_coeff, cg_none0_list, &
                                  ncg_none0, rab, swork_cont, Waux_mat, dWaux_mat, dsabb)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: la, first_sgfa
      INTEGER, INTENT(IN)                                :: nshella
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lb, first_sgfb
      INTEGER, INTENT(IN)                                :: nshellb
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lcb, first_sgfcb
      INTEGER, INTENT(IN)                                :: nshellcb
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: cg_coeff
      INTEGER, DIMENSION(:, :, :), INTENT(IN)            :: cg_none0_list
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: ncg_none0
      REAL(KIND=dp), INTENT(IN)                          :: rab(3)
      REAL(KIND=dp), DIMENSION(:, 0:, :, :, :), &
         INTENT(IN)                                      :: swork_cont
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Waux_mat
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: dWaux_mat
      REAL(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(INOUT)                                   :: dsabb

      CHARACTER(len=*), PARAMETER :: routineN = 'dev_overlap_shg_abb', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, i, ia, il, ilist, ishella, isob1, isob2, isobb, j, jb, jshellb, kb, &
         kshellb, labmin, lai, lbb, lbj, lbk, mai, mbb, mbj, mbk, nla, nlb, sgfa, sgfb, sgfcb
      REAL(KIND=dp)                                      :: dprefac, dtemp(3), prefac, rabx2(3), &
                                                            tsign

      CALL timeset(routineN, handle)

      rabx2(:) = 2.0_dp*rab

      DO kshellb = 1, nshellcb
         lbk = lcb(kshellb)
         sgfcb = first_sgfcb(kshellb)
         kb = sgfcb+lbk
         DO jshellb = 1, nshellb
            lbj = lb(jshellb)
            sgfb = first_sgfb(jshellb)
            jb = sgfb+lbj
            DO ishella = 1, nshella
               lai = la(ishella)
               nla = nsoset(lai-1)+lai+1
               sgfa = first_sgfa(ishella)
               ia = sgfa+lai
               DO mbj = -lbj, lbj, 1
                  DO mbk = -lbk, lbk, 1
                     isob1 = indso_inv(lbj, mbj)
                     isob2 = indso_inv(lbk, mbk)
                     DO mai = -lai, lai, 1
                        DO ilist = 1, ncg_none0(isob1, isob2)
                           isobb = cg_none0_list(isob1, isob2, ilist)
                           lbb = indso(1, isobb)
                           mbb = indso(2, isobb)
                           nlb = nsoset(lbb-1)+lbb+1
                           ! tsgin: because we take the transpose of auxmat (calculated for (la,lb))
                           tsign = 1.0_dp
                           IF (MODULO(lbb-lai, 2) /= 0) tsign = -1.0_dp
                           labmin = MIN(lai, lbb)
                           il = (lbj+lbk-lbb)/2
                           dtemp = 0.0_dp
                           DO j = 0, labmin
                              prefac = swork_cont(lai+lbb-j+1, il, ishella, jshellb, kshellb)
                              dprefac = swork_cont(lai+lbb-j+2, il, ishella, jshellb, kshellb)
                              DO i = 1, 3
                                 dtemp(i) = dtemp(i)+rabx2(i)*dprefac*Waux_mat(j+1, nlb+mbb, nla+mai) &
                                            +prefac*dWaux_mat(i, j+1, nlb+mbb, nla+mai)
                              ENDDO
                           ENDDO
                           DO i = 1, 3
                              dsabb(ia+mai, jb+mbj, kb+mbk, i) = dsabb(ia+mai, jb+mbj, kb+mbk, i) &
                                                                 +tsign*cg_coeff(isob1, isob2, isobb)*dtemp(i)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE dev_overlap_shg_abb

END MODULE construct_shg


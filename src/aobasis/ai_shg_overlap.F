!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the overlap integrals over solid harmonic Gaussian 
!>        (SHG) functions.
!> \par Literature (partly)
!>      T.J. Giese and D. M. York, J. Chem. Phys, 128, 064104 (2008)
!>      T. Helgaker, P Joergensen, J. Olsen, Molecular Electronic-Structure
!>                                           Theory, Wiley 
!> \par History
!>      created [04.2015]
!> \author Dorothea Golze
! *****************************************************************************
MODULE ai_shg_overlap
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: dfac,&
                                             fac,&
                                             pi
  USE orbital_pointers,                ONLY: indso,&
                                             indso_inv,&
                                             nsoset_pm
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_shg_overlap'

! *** Public subroutines ***
  PUBLIC :: get_dW_matrix, get_W_matrix, get_real_scaled_solid_harmonic,&
            overlap_shg_ab, dev_overlap_shg_ab, overlap_shg_aba, dev_overlap_shg_aba,&
            overlap_shg_abb, dev_overlap_shg_abb, s_overlap_ab, s_overlap_abx

CONTAINS

! *****************************************************************************
!> \brief calculates the uncontracted, not normalized [s|s] overlap 
!> \param la_max maximal l quantum number on a
!> \param npgfa number of primitive Gaussian on a 
!> \param zeta set of exponents on a 
!> \param lb_max maximal l quantum number on b
!> \param npgfb number of primitive Gaussian on a
!> \param zetb set of exponents on a 
!> \param rab distance vector between a and b
!> \param s uncontracted overlap of s functions
!> \param calculate_forces ...
! *****************************************************************************
  SUBROUTINE s_overlap_ab(la_max,npgfa,zeta,lb_max,npgfb,zetb,rab,s,calculate_forces)

    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta
    INTEGER, INTENT(IN)                      :: lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: s
    LOGICAL, INTENT(IN)                      :: calculate_forces

    CHARACTER(len=*), PARAMETER :: routineN = 's_overlap_ab', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ids, ipgfa, jpgfb, &
                                                ndev
    REAL(KIND=dp)                            :: a, b, rab2, xhi, zet

    CALL timeset(routineN,handle)

   ! Distance of the centers a and b
    rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
    ndev = 0
    IF(calculate_forces) ndev = 1
    ! Loops over all pairs of primitive Gaussian-type functions
    DO ipgfa=1,npgfa
       DO jpgfb=1,npgfb

          ! Distance Screening   !maybe later 

          ! Calculate some prefactors
          a = zeta(ipgfa)
          b = zetb(jpgfb)
          zet =  a+b
          xhi =  a*b/zet

          ! [s|s] integral
          s(ipgfa,jpgfb,1) = (pi/zet)**(1.5_dp)*EXP(-xhi*rab2) 

          DO ids=2,la_max + lb_max + ndev + 1
            s(ipgfa,jpgfb,ids) = -xhi*s(ipgfa,jpgfb,ids-1)
          ENDDO

       END DO
    END DO

   CALL timestop(handle)

  END SUBROUTINE s_overlap_ab

! *****************************************************************************
!> \brief calculates [s|ra^n|s] integrals for [aba] and the [s|rb^n|s] 
!>        integrals for [abb]
!>        abx ... x stands either for a and b respectively
!> \param la_max maximal l quantum number on a, orbital basis
!> \param npgfa number of primitive Gaussian on a, orbital basis 
!> \param zeta set of exponents on a, orbital basis 
!> \param lb_max maximal l quantum number on b, orbital basis
!> \param npgfb number of primitive Gaussian on a, orbital basis
!> \param zetb set of exponents on b, orbital basis 
!> \param lx_max maximal l quantum number of aux basis on x (a or b) 
!> \param npgfx number of primitive Gaussian on x, aux basis
!> \param zetx set of exponents on x, aux basis
!> \param rab distance vector between a and b
!> \param s uncontracted [s|r^n|s] integrals
!> \param calculate_forces ...
!> \param calc_aba flag if [aba] or [abb] is calculated
! *****************************************************************************
  SUBROUTINE s_overlap_abx(la_max,npgfa,zeta,lb_max,npgfb,zetb,lx_max,npgfx,zetx,&
                           rab,s,calculate_forces,calc_aba)

    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta
    INTEGER, INTENT(IN)                      :: lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb
    INTEGER, INTENT(IN)                      :: lx_max, npgfx
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetx
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: s
    LOGICAL, INTENT(IN)                      :: calculate_forces, calc_aba

    CHARACTER(len=*), PARAMETER :: routineN = 's_overlap_abx', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ids, il, ipgfa, &
                                                jpgfb, kpgfx, laa_max, &
                                                lbb_max, lmax, n, ndev, nds, &
                                                nl, stat
    REAL(KIND=dp)                            :: a, b, exp_rab2, k, pfac, &
                                                prefac, rab2, sqrt_pi3, &
                                                sqrt_zet, sr_int, temp, x, &
                                                xhi, zet
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dsr_int, dtemp

    CALL timeset(routineN,handle)

    ! Distance of the centers a and b
    rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
    ndev = 0
    IF(calculate_forces) ndev = 1

    IF(calc_aba) THEN 
      laa_max =  la_max +  lx_max
      nl = INT(laa_max/2)  
      IF(la_max ==0 .OR. lx_max ==0)  nl = 0
      lmax = laa_max + lb_max
    ELSE
      lbb_max =  lb_max +  lx_max
      nl = INT(lbb_max/2)  
      IF(lb_max ==0 .OR. lx_max ==0)  nl = 0
      lmax = la_max + lbb_max
    ENDIF

    ALLOCATE(dtemp(nl+1), dsr_int(nl+1),stat=STAT)
    CPASSERT(stat==0)

    sqrt_pi3 = SQRT(pi**3)

    ! Loops over all pairs of primitive Gaussian-type functions
    DO ipgfa=1,npgfa
       DO jpgfb=1,npgfb
         DO kpgfx=1,npgfx

            ! Calculate some prefactors
            IF(calc_aba) THEN
              a = zeta(ipgfa) + zetx(kpgfx) 
              b = zetb(jpgfb)
              x = b
            ELSE
              a = zeta(ipgfa)  
              b = zetb(jpgfb) + zetx(kpgfx)
              x = a
            ENDIF

            zet =  a+b
            xhi =  a*b/zet
            exp_rab2 = EXP(-xhi*rab2)

            pfac = x**2 / zet
            sqrt_zet = SQRT(zet)

            DO il=0,nl
               nds =  lmax - 2*il + ndev + 1
               SELECT CASE(il)
               CASE(0)
                 ! [s|s] integral
                 s(ipgfa,jpgfb,kpgfx,il,1) = (pi/zet)**(1.5_dp)*exp_rab2
                 DO ids = 2, nds
                   n = ids-1
                   s(ipgfa,jpgfb,kpgfx,il,ids) = (-xhi)**n*s(ipgfa,jpgfb,kpgfx,il,1)
                 ENDDO
               CASE(1)
                 ![s|r^2|s] integral
                 sr_int = sqrt_pi3/sqrt_zet**5 * (3.0_dp + 2.0_dp*pfac*rab2)/2.0_dp 
                 s(ipgfa,jpgfb,kpgfx,il,1) = exp_rab2 *  sr_int
                 k = sqrt_pi3 * x**2 / sqrt_zet**7
                 DO ids = 2, nds
                   n = ids-1
                   s(ipgfa,jpgfb,kpgfx,il,ids) = (-xhi)**n*exp_rab2*sr_int &
                                                + n * (-xhi)**(n-1)* k * exp_rab2
                 ENDDO 
               CASE(2)
                 ![s|r^4|s] integral
                  prefac = sqrt_pi3/4.0_dp/sqrt_zet**7
                  temp = 15.0_dp + 20.0_dp * pfac *rab2 + 4.0_dp * (pfac*rab2)**2
                  sr_int = prefac * temp
                  s(ipgfa,jpgfb,kpgfx,il,1) = exp_rab2 * sr_int
                  !** derivatives
                  k = sqrt_pi3 * x**4 / sqrt_zet**11
                  dsr_int(1) = prefac * (20.0_dp * pfac + 8.0_dp * pfac**2 * rab2)
                  DO ids = 2, nds
                    n = ids-1
                    dtemp(1) = (-xhi)**n*exp_rab2*sr_int
                    dtemp(2) = n * (-xhi)**(n-1) * exp_rab2 * dsr_int(1)
                    dtemp(3) = (n**2 - n) * (-xhi)**(n-2) * k * exp_rab2  
                    s(ipgfa,jpgfb,kpgfx,il,ids) = dtemp(1) +  dtemp(2) + dtemp(3) 
                  ENDDO
               CASE(3)
                ![s|r^6|s] integral
                 prefac = sqrt_pi3/8.0_dp/sqrt_zet**9
                 temp = 105.0_dp + 210.0_dp * pfac * rab2 
                 temp = temp + 84.0_dp * (pfac*rab2)**2 + 8.0_dp * (pfac*rab2)**3
                 sr_int = prefac * temp
                 s(ipgfa,jpgfb,kpgfx,il,1) = exp_rab2 * sr_int
                 !** derivatives
                 k = sqrt_pi3 * x**6 / sqrt_zet**15
                 dsr_int(1) = prefac * (210.0_dp*pfac + 168.0_dp*pfac**2*rab2 &
                                        + 24.0_dp * pfac**3 * rab2**2)
                 dsr_int(2) = prefac * (168.0_dp*pfac**2 + 48.0_dp * pfac**3 * rab2)
                 DO ids = 2, nds
                   n = ids-1
                   dtemp(1) = (-xhi)**n*exp_rab2*sr_int
                   dtemp(2) = REAL(n,dp) * (-xhi)**(n-1) * exp_rab2 * dsr_int(1)
                   dtemp(3) = REAL(n**2 - n,dp)/2.0_dp * (-xhi)**(n-2) &
                              * exp_rab2 * dsr_int(2)  
                   dtemp(4) = REAL(n*(n-1)*(n-2),dp) *  (-xhi)**(n-3) * k * exp_rab2 
                   s(ipgfa,jpgfb,kpgfx,il,ids) =   dtemp(1) + dtemp(2)&
                                                 + dtemp(3) + dtemp(4) 
                 ENDDO
               CASE(4)
                ![s|r^8|s] integral
                 prefac = sqrt_pi3/16.0_dp/sqrt_zet**11
                 temp = 945.0_dp + 2520.0_dp * pfac * rab2 + 1512.0_dp * (pfac*rab2)**2
                 temp = temp + 288.0_dp * (pfac*rab2)**3 + 16.0_dp * (pfac*rab2)**4
                 sr_int =  prefac * temp
                 s(ipgfa,jpgfb,kpgfx,il,1) = exp_rab2 * sr_int
                 !** derivatives
                 k = sqrt_pi3 * x**8 / sqrt_zet**19
                 dsr_int(1) = 2520.0_dp * pfac + 3024.0_dp * pfac**2 * rab2
                 dsr_int(1) = dsr_int(1) + 864.0_dp * pfac**3 * rab2**2&
                              + 64.0_dp * pfac**4 * rab2**3
                 dsr_int(1) = prefac * dsr_int(1)
                 dsr_int(2) = 3024.0_dp * pfac**2 + 1728.0_dp * pfac**3 * rab2 
                 dsr_int(2) = dsr_int(2) + 192.0_dp * pfac**4 * rab2**2
                 dsr_int(2) = prefac * dsr_int(2)
                 dsr_int(3) = 1728.0_dp * pfac**3 + 384.0_dp * pfac**4 * rab2
                 dsr_int(3) = prefac * dsr_int(3)
                 DO ids = 2, nds
                   n = ids-1
                   dtemp(1) = (-xhi)**n*exp_rab2*sr_int
                   dtemp(2) = REAL(n,dp) * (-xhi)**(n-1) * exp_rab2 * dsr_int(1)
                   dtemp(3) = REAL(n**2 - n,dp)/2.0_dp * (-xhi)**(n-2)&
                              * exp_rab2 * dsr_int(2)  
                   dtemp(4) = REAL(n*(n-1)*(n-2),dp)/6.0_dp * (-xhi)**(n-3)&
                              * exp_rab2 * dsr_int(3)
                   dtemp(5) = REAL(n*(n-1)*(n-2)*(n-3),dp) * (-xhi)**(n-4) &
                              * k * exp_rab2
                   s(ipgfa,jpgfb,kpgfx,il,ids) =  dtemp(1) + dtemp(2) + dtemp(3)&
                                                + dtemp(4) + dtemp(5) 
                 ENDDO
               CASE(5)
                ![s|r^10|s] integral
                 prefac = sqrt_pi3/32.0_dp/sqrt_zet**13
                 temp = 10395.0_dp + 34650.0_dp * pfac * rab2 
                 temp = temp + 27720.0_dp * (pfac*rab2)**2 + 7920.0_dp * (pfac*rab2)**3
                 temp = temp + 880.0_dp * (pfac*rab2)**4 + 32.0_dp * (pfac*rab2)**5
                 sr_int = prefac * temp 
                 s(ipgfa,jpgfb,kpgfx,il,1) = exp_rab2 * sr_int
                 !** derivatives
                 k = sqrt_pi3 * x**10 / sqrt_zet**23
                 dsr_int(1) = 34650.0_dp * pfac + 55440.0_dp * pfac**2 * rab2
                 dsr_int(1) = dsr_int(1) + 23760.0_dp * pfac**3 * rab2**2
                 dsr_int(1) = dsr_int(1) + 3520.0_dp * pfac**4 * rab2**3 
                 dsr_int(1) = dsr_int(1) + 160.0_dp * pfac**5 * rab2**4 
                 dsr_int(1) = prefac * dsr_int(1) 
                 dsr_int(2) = 55440.0_dp * pfac**2 + 47520.0_dp * pfac**3 * rab2 
                 dsr_int(2) = dsr_int(2) + 10560.0_dp * pfac**4 * rab2**2 
                 dsr_int(2) = dsr_int(2) + 640.0_dp * pfac**5 * rab2**3 
                 dsr_int(2) = prefac * dsr_int(2) 
                 dsr_int(3) = 47520.0_dp * pfac**3 + 21120.0_dp * pfac**4 *rab2
                 dsr_int(3) = dsr_int(3) + 1920.0_dp * pfac**5 * rab2**2
                 dsr_int(3) = prefac * dsr_int(3) 
                 dsr_int(4) = 21120.0_dp * pfac**4 + 3840.0_dp * pfac**5 *rab2
                 dsr_int(4) = prefac * dsr_int(4)
                 DO ids = 2, nds
                   n = ids-1
                   dtemp(1) = (-xhi)**n*exp_rab2*sr_int
                   dtemp(2) = REAL(n,dp) * (-xhi)**(n-1) * exp_rab2 * dsr_int(1)
                   dtemp(3) = REAL(n**2 - n,dp)/2.0_dp * (-xhi)**(n-2)&
                              * exp_rab2 * dsr_int(2)  
                   dtemp(4) = REAL(n*(n-1)*(n-2),dp)/6.0_dp * (-xhi)**(n-3)&
                              * exp_rab2 * dsr_int(3)
                   dtemp(5) = REAL(n*(n-1)*(n-2)*(n-3),dp)/24.0_dp * (-xhi)**(n-4)&
                              * exp_rab2 * dsr_int(4)
                   dtemp(6) = REAL(n*(n-1)*(n-2)*(n-3)*(n-4),dp) * (-xhi)**(n-5)&
                             * k *  exp_rab2
                   s(ipgfa,jpgfb,kpgfx,il,ids) =   dtemp(1) + dtemp(2) + dtemp(3)&
                                                 + dtemp(4) + dtemp(5) + dtemp(6)
                 ENDDO
               CASE DEFAULT
                 CALL cp_unimplemented_error(fromWhere=routineP, &
                      message="SHG integrals not implemented when sum of l quantum"//& 
                      " number of orbital and ri basis larger than 11", &
                      error_level=cp_failure_level)
               END SELECT
               
                
            ENDDO

          END DO
       END DO
    END DO

    DEALLOCATE(dtemp, dsr_int, stat=STAT)
    CPASSERT(stat==0)

    CALL timestop(handle)

  END SUBROUTINE s_overlap_abx
 
! *****************************************************************************
!> \brief computes the real scaled solid harmonics Rlm up to a given l
!> \param Rlm_c cosine part of real scaled soldi harmonics 
!> \param Rlm_s sine part of real scaled soldi harmonics 
!> \param l maximal l quantum up to where Rlm is calculated
!> \param r distance vector between a and b
!> \param r2 square of distance vector
! *****************************************************************************
  SUBROUTINE get_real_scaled_solid_harmonic(Rlm_c,Rlm_s,l,r,r2)

    INTEGER, INTENT(IN)                      :: l
    REAL(KIND=dp), &
      DIMENSION(0:l, -2*l:2*l), INTENT(OUT)  :: Rlm_s, Rlm_c
    REAL(KIND=dp), DIMENSION(3)              :: r
    REAL(KIND=dp)                            :: r2

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'get_real_scaled_solid_harmonic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, li, mi, prefac
    REAL(KIND=dp)                            :: Rc, Rc_00, Rlm, Rmlm, Rplm, &
                                                Rs, Rs_00, temp_c

    CALL timeset(routineN,handle)

    Rc_00 = 1.0_dp
    Rs_00 = 0.0_dp
  
    Rlm_c(0,0) = Rc_00 
    Rlm_s(0,0) = Rs_00
 
 
    ! generate elements Rmm 
    ! start
    IF (l > 0) THEN
      Rc = -0.5_dp*r(1)*Rc_00 
      Rs = -0.5_dp*r(2)*Rc_00
      Rlm_c(1,1) = Rc
      Rlm_s(1,1) = Rs
      Rlm_c(1,-1) = -Rc
      Rlm_s(1,-1) =  Rs
    ENDIF
    DO li=2,l
       temp_c = (-r(1)*Rc + r(2)*Rs)/(REAL(2*(li-1) + 2,dp))
       Rs = (-r(2)*Rc - r(1)*Rs)/(REAL(2*(li-1) + 2,dp))
       Rc = temp_c
       Rlm_c(li,li) = Rc 
       Rlm_s(li,li) = Rs
       IF(MODULO(li,2) /= 0) THEN
         Rlm_c(li,-li) = -Rc
         Rlm_s(li,-li) =  Rs
       ELSE
         Rlm_c(li,-li) =  Rc
         Rlm_s(li,-li) = -Rs
       ENDIF 
    ENDDO

    DO mi =0,l-1
      Rmlm = Rlm_c(mi,mi)
      Rlm = r(3)*Rlm_c(mi,mi)
      Rlm_c(mi+1,mi) = Rlm
      IF(MODULO(mi,2) /=0) THEN
        Rlm_c(mi+1,-mi) = -Rlm
      ELSE
        Rlm_c(mi+1,-mi) =  Rlm
      ENDIF
      DO li = mi+2,l
         prefac = (li+mi)*(li-mi) 
         Rplm = (REAL(2*li-1,dp)*r(3)*Rlm - r2*Rmlm)/REAL(prefac,dp) 
         Rmlm =  Rlm
         Rlm  = Rplm
         Rlm_c(li,mi) = Rlm 
         IF(MODULO(mi,2) /= 0) THEN
           Rlm_c(li,-mi) = -Rlm
         ELSE
           Rlm_c(li,-mi) =  Rlm
         ENDIF 
      ENDDO
    ENDDO
    DO mi =1,l-1
       Rmlm = Rlm_s(mi,mi)
       Rlm = r(3)*Rlm_s(mi,mi)
       Rlm_s(mi+1,mi) = Rlm
       IF(MODULO(mi,2) /= 0) THEN
         Rlm_s(mi+1,-mi) =  Rlm
       ELSE
         Rlm_s(mi+1,-mi) = -Rlm
       ENDIF 
       DO li = mi+2,l
          prefac = (li+mi)*(li-mi) 
          Rplm = (REAL(2*li-1,dp)*r(3)*Rlm - r2*Rmlm)/REAL(prefac,dp) 
          Rmlm =  Rlm
          Rlm  = Rplm
          Rlm_s(li,mi) = Rlm 
          IF(MODULO(mi,2) /= 0) THEN
            Rlm_s(li,-mi) =  Rlm
          ELSE
            Rlm_s(li,-mi) = -Rlm
          ENDIF 
       ENDDO
    ENDDO

   CALL timestop(handle)

  END SUBROUTINE get_real_scaled_solid_harmonic

! *****************************************************************************
!> \brief calculates the angular dependent-part of the SHG integrals,
!>        transformation matrix W, see literature above
!> \param lamax array of maximal l quantum number on a;
!>        lamax(lb) with lb= 0..lbmax
!> \param lbmax maximal l quantum number on b
!> \param lmax maximal l quantum number
!> \param Rc cosine part of real scaled solid harmonics
!> \param Rs sine part of real scaled solid harmonics
!> \param Waux_mat stores the angular-dependent part of the SHG integrals
!>        last dimension is (1:4): cc(1), cs(2), sc(3), ss(4)
! *****************************************************************************
  SUBROUTINE get_W_matrix(lamax,lbmax,lmax,Rc,Rs,Waux_mat)

    INTEGER, DIMENSION(:), POINTER           :: lamax
    INTEGER, INTENT(IN)                      :: lbmax, lmax
    REAL(KIND=dp), &
      DIMENSION(0:lmax, -2*lmax:2*lmax), &
      INTENT(IN)                             :: Rc, Rs
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: Waux_mat

    CHARACTER(len=*), PARAMETER :: routineN = 'get_W_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, j, k, la, labmin, &
                                                laj, lb, lbj, ma, ma_m, ma_p, &
                                                mb, mb_m, mb_p, nla, nlb
    REAL(KIND=dp) :: A_jk, A_lama, A_lbmb, Alm_fac, delta_k, prefac, Rca_m, &
      Rca_p, Rcb_m, Rcb_p, Rsa_m, Rsa_p, Rsb_m, Rsb_p, sign_fac, Wa(4), &
      Wb(4), Wmat(4)

    CALL timeset(routineN,handle)

    Wa(:) = 0.0_dp
    Wb(:) = 0.0_dp
    Wmat(:) = 0.0_dp

    DO lb = 0,lbmax
       nlb = nsoset_pm(lb-1)
       DO la = 0,lamax(lb)
          nla = nsoset_pm(la-1)
          labmin = MIN(la,lb)
          DO mb = 0,lb
             A_lbmb = SQRT(fac(lb+mb)*fac(lb-mb))
             IF(MODULO(mb,2) /= 0 )  A_lbmb = - A_lbmb 
             IF(mb /= 0) A_lbmb = A_lbmb * SQRT(2.0_dp)
             IF(MODULO(lb,2) /= 0) A_lbmb = - A_lbmb
             DO ma = 0,la
                A_lama = SQRT(fac(la+ma)*fac(la-ma))
                IF(MODULO(ma,2) /= 0 )  A_lama = - A_lama 
                IF(ma /= 0) A_lama = A_lama * SQRT(2.0_dp)
                Alm_fac = A_lama * A_lbmb
                DO j=0,labmin
                   laj = la - j 
                   lbj = lb - j
                   prefac = Alm_fac * REAL(2**(la+lb-j),dp) * dfac(2*j-1) 
                   delta_k = 0.5_dp
                   Wmat = 0.0_dp
                   DO k=0,j
                      ma_m = ma - k 
                      ma_p = ma + k 
                      IF(laj< ABS(ma_m) .AND. laj < ABS(ma_p)) CYCLE
                      mb_m = mb - k 
                      mb_p = mb + k 
                      IF(lbj < ABS(mb_m) .AND. lbj < ABS(mb_p)) CYCLE
                      IF(k /= 0 ) delta_k = 1.0_dp
                      A_jk = fac(j+k)*fac(j-k)
                      IF(k /= 0) A_jk = 2.0_dp*A_jk 
                      IF(MODULO(k,2) /= 0) THEN
                        sign_fac = -1.0_dp   
                      ELSE
                        sign_fac = 1.0_dp   
                      ENDIF
                      Rca_m = Rc(laj,ma_m)
                      Rsa_m = Rs(laj,ma_m)
                      Rca_p = Rc(laj,ma_p)
                      Rsa_p = Rs(laj,ma_p)
                      Rcb_m = Rc(lbj,mb_m)
                      Rsb_m = Rs(lbj,mb_m)
                      Rcb_p = Rc(lbj,mb_p)
                      Rsb_p = Rs(lbj,mb_p)
                      Wa(1) = delta_k*(Rca_m + sign_fac * Rca_p) 
                      Wb(1) = delta_k*(Rcb_m + sign_fac * Rcb_p)
                      Wa(2) = -Rsa_m + sign_fac * Rsa_p
                      Wb(2) = -Rsb_m + sign_fac * Rsb_p
                      Wmat(1) = Wmat(1) + prefac/A_jk * (Wa(1) * Wb(1) + Wa(2) * Wb(2))
                      IF(mb > 0) THEN 
                        Wb(3) = delta_k*(Rsb_m + sign_fac * Rsb_p)
                        Wb(4) = Rcb_m - sign_fac * Rcb_p
                        Wmat(2) = Wmat(2) + prefac/A_jk * (Wa(1) * Wb(3) + Wa(2) * Wb(4))
                      ENDIF
                      IF(ma > 0) THEN 
                        Wa(3) = delta_k*(Rsa_m + sign_fac * Rsa_p)
                        Wa(4) =  Rca_m - sign_fac * Rca_p
                        Wmat(3) = Wmat(3) + prefac/A_jk * (Wa(3) * Wb(1) + Wa(4) * Wb(2))
                      ENDIF
                      IF(ma > 0 .AND. mb > 0) THEN
                       Wmat(4) = Wmat(4)  + prefac/A_jk * (Wa(3) * Wb(3) + Wa(4) * Wb(4))
                      ENDIF
                   ENDDO
                   Waux_mat(nla+ma+1,nlb+mb+1,j+1,1) = Wmat(1)
                   IF(mb > 0) Waux_mat(nla+ma+1,nlb+mb+1,j+1,2) = Wmat(2)
                   IF(ma > 0) Waux_mat(nla+ma+1,nlb+mb+1,j+1,3) = Wmat(3)
                   IF(ma > 0.AND. mb > 0) Waux_mat(nla+ma+1,nlb+mb+1,j+1,4) = Wmat(4)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE get_W_matrix

! *****************************************************************************
!> \brief calculates derivatives of transformation matrix W,
!> \param lamax array of maximal l quantum number on a;
!>        lamax(lb) with lb= 0..lbmax
!> \param lbmax maximal l quantum number on b
!> \param Waux_mat stores the angular-dependent part of the SHG integrals 
!> \param dWaux_mat stores the derivatives of the angular-dependent part of 
!>        the SHG integrals 
!>        last dimension is (1:4): cc(1), cs(2), sc(3), ss(4)
! *****************************************************************************
  SUBROUTINE get_dW_matrix(lamax,lbmax,Waux_mat,dWaux_mat)

    INTEGER, DIMENSION(:), POINTER           :: lamax
    INTEGER, INTENT(IN)                      :: lbmax
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: Waux_mat
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: dWaux_mat

    CHARACTER(len=*), PARAMETER :: routineN = 'get_dW_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ia, iam, ib, ibm, j, &
                                                jmax, la, labm, labmin, lamb, &
                                                lb, ma, mb, nla, nlam, nlb, &
                                                nlbm, stat
    REAL(KIND=dp)                            :: bma, bma_m, bma_p, bmb, &
                                                bmb_m, bmb_p, dAa, dAa_m, &
                                                dAa_p, dAb, dAb_m, dAb_p
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Wam, Wamm, Wamp, Wbm, Wbmm, &
                                                Wbmp

    CALL timeset(routineN,handle)

    jmax = MIN(MAXVAL(lamax),lbmax)
    ALLOCATE(Wam(0:jmax,4),Wamm(0:jmax,4),Wamp(0:jmax,4),STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(Wbm(0:jmax,4),Wbmm(0:jmax,4),Wbmp(0:jmax,4),STAT=stat)
    CPASSERT(stat==0)

    DO lb = 0,lbmax
       nlb = nsoset_pm(lb-1)
       nlbm = 0
       IF( lb > 0) nlbm = nsoset_pm(lb-2)
       DO la = 0,lamax(lb)
          nla = nsoset_pm(la-1)
          nlam = 0
          IF(la > 0) nlam = nsoset_pm(la-2)
          labmin = MIN(la,lb)
          lamb = MIN(la-1,lb)
          labm = MIN(la,lb-1)
          DO mb = 0,lb
             bmb   = 1.0_dp
             bmb_m = 1.0_dp
             bmb_p = 1.0_dp
             IF(mb /= 0) bmb = SQRT(2.0_dp)
             IF(mb-1 /= 0) bmb_m = SQRT(2.0_dp)
             IF(mb+1 /= 0) bmb_p = SQRT(2.0_dp)
             dAb_p = -bmb/bmb_p*SQRT(REAL((lb-mb)*(lb-mb-1),dp))
             dAb_m = -bmb/bmb_m*SQRT(REAL((lb+mb)*(lb+mb-1),dp))
             dAb   = 2.0_dp*SQRT(REAL((lb+mb)*(lb-mb),dp))
             !*** for m=0 W_l-1,-1 can't be read from Waux_mat, but we use
             !*** W_l-1,-1 = -W_l-1,1 [cc(1), cs(2)] or W_l-1,-1 = W_l-1,1 [[sc(3), ss(4)]
             IF(mb ==0) dAb_p = 2.0_dp * dAb_p 
             ib = nlb + mb + 1
             ibm = nlbm + mb + 1
             DO ma = 0,la
                bma   = 1.0_dp
                bma_m = 1.0_dp
                bma_p = 1.0_dp
                IF(ma /=0) bma = SQRT(2.0_dp)
                IF(ma-1 /= 0) bma_m = SQRT(2.0_dp)
                IF(ma+1 /= 0) bma_p = SQRT(2.0_dp)
                dAa_p = -bma/bma_p*SQRT(REAL((la-ma)*(la-ma-1),dp))
                dAa_m = -bma/bma_m*SQRT(REAL((la+ma)*(la+ma-1),dp))
                dAa = 2.0_dp*SQRT(REAL((la+ma)*(la-ma),dp))
                IF(ma ==0) dAa_p = 2.0_dp * dAa_p 
                ia = nla + ma + 1
                iam = nlam + ma + 1
                Wam(:,:)  = 0.0_dp
                Wamm(:,:) = 0.0_dp
                Wamp(:,:) = 0.0_dp
                !*** Wam: la-1, ma
                IF(ma <= la-1)  THEN
                   Wam(0:lamb,1)  = Waux_mat(iam,ib,1:lamb+1,1)
                   IF(mb > 0) Wam(0:lamb,2)  = Waux_mat(iam,ib,1:lamb+1,2)
                   IF(ma > 0) Wam(0:lamb,3)  = Waux_mat(iam,ib,1:lamb+1,3)
                   IF(ma > 0.AND.mb >0) Wam(0:lamb,4)  = Waux_mat(iam,ib,1:lamb+1,4)
                ENDIF
                !*** Wamm: la-1, ma-1 
                IF(ma-1 >= 0)  THEN
                   Wamm(0:lamb,1) = Waux_mat(iam-1,ib,1:lamb+1,1)
                   IF(mb > 0) Wamm(0:lamb,2) = Waux_mat(iam-1,ib,1:lamb+1,2)
                   IF(ma-1 > 0 ) Wamm(0:lamb,3) = Waux_mat(iam-1,ib,1:lamb+1,3)
                   IF(ma-1 > 0.AND.mb > 0) Wamm(0:lamb,4) = Waux_mat(iam-1,ib,1:lamb+1,4)
                ENDIF
                !*** Wamp: la-1, ma+1
                IF(ma+1 <= la-1) THEN
                   Wamp(0:lamb,1) = Waux_mat(iam+1,ib,1:lamb+1,1)
                   IF(mb > 0) Wamp(0:lamb,2) = Waux_mat(iam+1,ib,1:lamb+1,2)
                   IF(ma+1 > 0) Wamp(0:lamb,3) = Waux_mat(iam+1,ib,1:lamb+1,3)
                   IF(ma+1 > 0 .AND. mb >0) Wamp(0:lamb,4) = Waux_mat(iam+1,ib,1:lamb+1,4)
                ENDIF
                Wbm(:,:)  = 0.0_dp
                Wbmm(:,:) = 0.0_dp
                Wbmp(:,:) = 0.0_dp
                !*** Wbm: lb-1, mb
                IF(mb <= lb-1) THEN
                   Wbm(0:labm,1)  = Waux_mat(ia,ibm,1:labm+1,1)
                   IF(mb > 0) Wbm(0:labm,2)  = Waux_mat(ia,ibm,1:labm+1,2)
                   IF(ma > 0) Wbm(0:labm,3)  = Waux_mat(ia,ibm,1:labm+1,3)
                   IF(ma > 0 .AND. mb > 0) Wbm(0:labm,4)  = Waux_mat(ia,ibm,1:labm+1,4)
                ENDIF
                !*** Wbmm: lb-1, mb-1
                IF(mb-1 >= 0) THEN
                   Wbmm(0:labm,1) = Waux_mat(ia,ibm-1,1:labm+1,1)
                   IF(mb-1 > 0) Wbmm(0:labm,2) = Waux_mat(ia,ibm-1,1:labm+1,2)
                   IF(ma > 0) Wbmm(0:labm,3) = Waux_mat(ia,ibm-1,1:labm+1,3)
                   IF(ma > 0.AND.mb-1 > 0)  Wbmm(0:labm,4) = Waux_mat(ia,ibm-1,1:labm+1,4)
                ENDIF
                !*** Wbmp: lb-1, mb+1
                IF(mb+1 <= lb-1) THEN
                   Wbmp(0:labm,1) = Waux_mat(ia,ibm+1,1:labm+1,1)
                   IF(mb+1 > 0) Wbmp(0:labm,2) = Waux_mat(ia,ibm+1,1:labm+1,2)
                   IF(ma > 0) Wbmp(0:labm,3) = Waux_mat(ia,ibm+1,1:labm+1,3)
                   IF(ma > 0.AND. mb+1 > 0) Wbmp(0:labm,4) = Waux_mat(ia,ibm+1,1:labm+1,4)
                ENDIF
                DO j=0,labmin
                   !*** x component
                   dWaux_mat(1,ia,ib,j+1,1) =  dAa_p * Wamp(j,1) - dAa_m * Wamm(j,1) &
                                             - dAb_p * Wbmp(j,1) + dAb_m * Wbmm(j,1)
                   IF(mb > 0) THEN
                     dWaux_mat(1,ia,ib,j+1,2) =  dAa_p * Wamp(j,2) - dAa_m * Wamm(j,2) &
                                               - dAb_p * Wbmp(j,2) + dAb_m * Wbmm(j,2) 
                   ENDIF
                   IF(ma > 0) THEN
                     dWaux_mat(1,ia,ib,j+1,3) =  dAa_p * Wamp(j,3) - dAa_m * Wamm(j,3) &
                                               - dAb_p * Wbmp(j,3) + dAb_m * Wbmm(j,3) 
                   ENDIF
                   IF(ma > 0.AND. mb > 0) THEN
                     dWaux_mat(1,ia,ib,j+1,4) =  dAa_p * Wamp(j,4) - dAa_m * Wamm(j,4) &
                                               - dAb_p * Wbmp(j,4) + dAb_m * Wbmm(j,4)
                   ENDIF
                   !**** y component
                   dWaux_mat(2,ia,ib,j+1,1) =  dAa_p * Wamp(j,3) + dAa_m * Wamm(j,3) &
                                             - dAb_p * Wbmp(j,2) - dAb_m * Wbmm(j,2) 
                   IF(mb > 0) THEN
                     dWaux_mat(2,ia,ib,j+1,2) =  dAa_p * Wamp(j,4) + dAa_m * Wamm(j,4) &
                                               + dAb_p * Wbmp(j,1) + dAb_m * Wbmm(j,1) 
                   ENDIF
                   IF(ma > 0) THEN
                     dWaux_mat(2,ia,ib,j+1,3) = - dAa_p * Wamp(j,1) - dAa_m * Wamm(j,1) &
                                                - dAb_p * Wbmp(j,4) - dAb_m * Wbmm(j,4) 
                   ENDIF
                   IF(ma > 0.AND. mb > 0) THEN
                     dWaux_mat(2,ia,ib,j+1,4) = - dAa_p * Wamp(j,2) - dAa_m * Wamm(j,2) &
                                                + dAb_p * Wbmp(j,3) + dAb_m * Wbmm(j,3) 
                   ENDIF
                   !**** z compnent
                   dWaux_mat(3,ia,ib,j+1,1) = dAa*Wam(j,1) - dAb*Wbm(j,1)
                   IF(mb > 0) THEN
                      dWaux_mat(3,ia,ib,j+1,2) = dAa*Wam(j,2) - dAb*Wbm(j,2)
                   ENDIF
                   IF(ma > 0) THEN
                      dWaux_mat(3,ia,ib,j+1,3) = dAa*Wam(j,3) - dAb*Wbm(j,3)
                   ENDIF
                   IF(ma > 0.AND. mb > 0) THEN
                      dWaux_mat(3,ia,ib,j+1,4) = dAa*Wam(j,4) - dAb*Wbm(j,4)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(Wam,Wamm,Wamp,STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(Wbm,Wbmm,Wbmp,STAT=stat)
    CPASSERT(stat==0)

    CALL timestop(handle)

  END SUBROUTINE get_dW_matrix

! *****************************************************************************
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
!> \param sab contracted overlap of spherical harmonic Gaussians
! *****************************************************************************
  SUBROUTINE overlap_shg_ab(la,first_sgfa,nshella,lb,first_sgfb,nshellb,&
                            swork_cont,Waux_mat,sab)

    INTEGER, DIMENSION(:), INTENT(IN)        :: la, first_sgfa
    INTEGER, INTENT(IN)                      :: nshella
    INTEGER, DIMENSION(:), INTENT(IN)        :: lb, first_sgfb
    INTEGER, INTENT(IN)                      :: nshellb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: swork_cont
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: Waux_mat
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: sab

    CHARACTER(len=*), PARAMETER :: routineN = 'overlap_shg_ab', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ia, iaa, ishella, j, &
                                                jb, jbb, jshellb, labmin, &
                                                lai, lbj, mai, mbj, nla, nlb, &
                                                sgfa, sgfb
    REAL(KIND=dp)                            :: prefac

    CALL timeset(routineN,handle)

    DO ishella = 1, nshella
       lai = la(ishella)
       nla = nsoset_pm(lai-1)
       sgfa = first_sgfa(ishella)
       ia = sgfa + lai
       DO jshellb = 1, nshellb
          lbj = lb(jshellb)
          nlb = nsoset_pm(lbj-1)
          sgfb = first_sgfb(jshellb)
          jb = sgfb + lbj
          labmin = MIN(lai,lbj)
          DO j=0,labmin
             prefac = swork_cont(ishella,jshellb,lai+lbj-j+1)  
             ! mai = 0 and mbj = 0
             iaa = nla+1
             jbb = nlb+1
             sab(ia,jb) = sab(ia,jb) + prefac * Waux_mat(iaa,jbb,j+1,1)
             ! mai = 0 and mbj /= 0
             DO mbj=1,lbj 
                jbb = nlb+mbj+1
                sab(ia,jb+mbj) = sab(ia,jb+mbj) + prefac * Waux_mat(iaa,jbb,j+1,1)
                sab(ia,jb-mbj) = sab(ia,jb-mbj) + prefac * Waux_mat(iaa,jbb,j+1,2) 
             ENDDO
             ! mai /= 0 and mbj = 0
             jbb = nlb+1
             DO mai=1,lai
                iaa = nla+mai+1
                sab(ia+mai,jb) = sab(ia+mai,jb) + prefac * Waux_mat(iaa,jbb,j+1,1)
                sab(ia-mai,jb) = sab(ia-mai,jb) + prefac * Waux_mat(iaa,jbb,j+1,3)
             ENDDO
             ! mai /= 0 and mbj /= 0
             DO mai=1,lai
                iaa = nla+mai+1
                DO mbj=1,lbj
                  jbb = nlb+mbj+1
                   sab(ia+mai,jb+mbj) = sab(ia+mai,jb+mbj) + prefac * Waux_mat(iaa,jbb,j+1,1)
                   sab(ia+mai,jb-mbj) = sab(ia+mai,jb-mbj) + prefac * Waux_mat(iaa,jbb,j+1,2)
                   sab(ia-mai,jb+mbj) = sab(ia-mai,jb+mbj) + prefac * Waux_mat(iaa,jbb,j+1,3)
                   sab(ia-mai,jb-mbj) = sab(ia-mai,jb-mbj) + prefac * Waux_mat(iaa,jbb,j+1,4)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    CALL timestop(handle)

  END SUBROUTINE overlap_shg_ab

! *****************************************************************************
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
!> \param dWaux_mat derivatives of precomputed angular-dependent part
!> \param dsab derivative of contracted overlap of spherical harmonic Gaussians 
! *****************************************************************************
  SUBROUTINE dev_overlap_shg_ab(la,first_sgfa,nshella,lb,first_sgfb,nshellb,rab,&
                            swork_cont,Waux_mat,dWaux_mat,dsab)

    INTEGER, DIMENSION(:), INTENT(IN)        :: la, first_sgfa
    INTEGER, INTENT(IN)                      :: nshella
    INTEGER, DIMENSION(:), INTENT(IN)        :: lb, first_sgfb
    INTEGER, INTENT(IN)                      :: nshellb
    REAL(KIND=dp), INTENT(IN)                :: rab(3)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: swork_cont
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: Waux_mat
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), INTENT(IN)   :: dWaux_mat
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: dsab

    CHARACTER(len=*), PARAMETER :: routineN = 'dev_overlap_shg_ab', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ia, iaa, ishella, &
                                                j, jb, jbb, jshellb, labmin, &
                                                lai, lbj, mai, mbj, nla, nlb, &
                                                sgfa, sgfb
    REAL(KIND=dp)                            :: dprefac, prefac, rabx2(3)

    CALL timeset(routineN,handle)

    rabx2(:) = 2.0_dp * rab
    DO ishella = 1, nshella
       lai = la(ishella)
       nla = nsoset_pm(lai-1)
       sgfa = first_sgfa(ishella)
       ia = sgfa + lai
       DO jshellb = 1, nshellb
          lbj = lb(jshellb)
          nlb = nsoset_pm(lbj-1)
          sgfb = first_sgfb(jshellb)
          jb = sgfb + lbj
          labmin = MIN(lai,lbj)
          DO j=0,labmin
             prefac =  swork_cont(ishella,jshellb,lai+lbj-j+1)  
             dprefac = swork_cont(ishella,jshellb,lai+lbj-j+2) !j+1 
             ! mai = 0 and mbj = 0
             iaa = nla+1
             jbb = nlb+1
             DO i=1,3
              dsab(ia,jb,i) = dsab(ia,jb,i) + rabx2(i) * dprefac * Waux_mat(iaa,jbb,j+1,1)&
                                            + prefac * dWaux_mat(i,iaa,jbb,j+1,1)
             ENDDO
             ! mai = 0 and mbj /= 0
             DO mbj=1,lbj 
                jbb = nlb+mbj+1
                DO i=1,3
                  dsab(ia,jb+mbj,i) = dsab(ia,jb+mbj,i) + rabx2(i) * dprefac * Waux_mat(iaa,jbb,j+1,1)&
                                                        + prefac * dWaux_mat(i,iaa,jbb,j+1,1)
                  dsab(ia,jb-mbj,i) = dsab(ia,jb-mbj,i) + rabx2(i) * dprefac * Waux_mat(iaa,jbb,j+1,2)& 
                                                        + prefac * dWaux_mat(i,iaa,jbb,j+1,2)
                ENDDO
             ENDDO
             ! mai /= 0 and mbj = 0
             jbb = nlb+1
             DO mai=1,lai
                iaa = nla+mai+1
                DO i=1,3
                  dsab(ia+mai,jb,i) = dsab(ia+mai,jb,i) + rabx2(i) * dprefac * Waux_mat(iaa,jbb,j+1,1)&
                                                        + prefac * dWaux_mat(i,iaa,jbb,j+1,1)
                  dsab(ia-mai,jb,i) = dsab(ia-mai,jb,i) + rabx2(i) * dprefac * Waux_mat(iaa,jbb,j+1,3)&
                                                        + prefac * dWaux_mat(i,iaa,jbb,j+1,3)
                ENDDO
             ENDDO
             ! mai /= 0 and mbj /= 0
             DO mai=1,lai
                iaa = nla+mai+1
                DO mbj=1,lbj
                  jbb = nlb+mbj+1
                  DO i=1,3
                   dsab(ia+mai,jb+mbj,i) = dsab(ia+mai,jb+mbj,i) + rabx2(i) * dprefac * Waux_mat(iaa,jbb,j+1,1)&
                                                                 + prefac * dWaux_mat(i,iaa,jbb,j+1,1)
                   dsab(ia+mai,jb-mbj,i) = dsab(ia+mai,jb-mbj,i) + rabx2(i) * dprefac * Waux_mat(iaa,jbb,j+1,2)&
                                                                 + prefac * dWaux_mat(i,iaa,jbb,j+1,2)
                   dsab(ia-mai,jb+mbj,i) = dsab(ia-mai,jb+mbj,i) + rabx2(i) * dprefac * Waux_mat(iaa,jbb,j+1,3)&
                                                                 + prefac * dWaux_mat(i,iaa,jbb,j+1,3)
                   dsab(ia-mai,jb-mbj,i) = dsab(ia-mai,jb-mbj,i) + rabx2(i) * dprefac * Waux_mat(iaa,jbb,j+1,4)&
                                                                 + prefac * dWaux_mat(i,iaa,jbb,j+1,4)
                  ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    CALL timestop(handle)

  END SUBROUTINE dev_overlap_shg_ab

! *****************************************************************************
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
! *****************************************************************************
  SUBROUTINE overlap_shg_aba(la,first_sgfa,nshella,lb,first_sgfb,nshellb,&
                             lca,first_sgfca,nshellca,cg_coeff,cg_none0_list,&
                             ncg_none0,swork_cont,Waux_mat,saba)

    INTEGER, DIMENSION(:), INTENT(IN)        :: la, first_sgfa
    INTEGER, INTENT(IN)                      :: nshella
    INTEGER, DIMENSION(:), INTENT(IN)        :: lb, first_sgfb
    INTEGER, INTENT(IN)                      :: nshellb
    INTEGER, DIMENSION(:), INTENT(IN)        :: lca, first_sgfca
    INTEGER, INTENT(IN)                      :: nshellca
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: cg_coeff
    INTEGER, DIMENSION(:, :, :), INTENT(IN)  :: cg_none0_list
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: ncg_none0
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: swork_cont
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: Waux_mat
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: saba

    CHARACTER(len=*), PARAMETER :: routineN = 'overlap_shg_aba', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ia, il, ilist, ishella, isoa1, isoa2, isoaa, j, jb, &
      jshellb, ka, kshella, laa, labmin, lai, lak, lbj, maa, mai, mak, mbj, &
      nla, nlb, sgfa, sgfb, sgfca
    REAL(KIND=dp)                            :: prefac, stemp

    CALL timeset(routineN,handle)

    DO ishella = 1, nshella
       lai = la(ishella)
       sgfa = first_sgfa(ishella)
       ia = sgfa + lai
       DO jshellb = 1, nshellb
          lbj = lb(jshellb)
          nlb = nsoset_pm(lbj-1)
          sgfb = first_sgfb(jshellb)
          jb = sgfb + lbj
          DO kshella = 1, nshellca
             lak = lca(kshella)
             sgfca = first_sgfca(kshella)
             ka = sgfca + lak
             DO mai = -lai,lai,1
                DO mbj=-lbj,lbj,1
                   DO mak=-lak,lak,1
                      isoa1 = indso_inv(lai,mai)
                      isoa2 = indso_inv(lak,mak)
                      DO ilist = 1, ncg_none0(isoa1,isoa2)
                         isoaa = cg_none0_list(isoa1,isoa2,ilist) 
                         laa = indso(1,isoaa) 
                         maa = indso(2,isoaa) 
                         nla = nsoset_pm(laa-1)
                         labmin = MIN(laa,lbj)
                         il = INT((lai + lak - laa)/2) 
                         stemp = 0.0_dp
                         DO j=0,labmin
                            prefac = swork_cont(ishella,jshellb,kshella,il,laa+lbj-j+1) 
                            IF(maa >= 0.AND.mbj >= 0) stemp = stemp + prefac * Waux_mat(nla+ABS(maa)+1,nlb+ABS(mbj)+1,j+1,1)
                            IF(maa >= 0.AND.mbj  < 0) stemp = stemp + prefac * Waux_mat(nla+ABS(maa)+1,nlb+ABS(mbj)+1,j+1,2)
                            IF(maa <  0.AND.mbj >= 0) stemp = stemp + prefac * Waux_mat(nla+ABS(maa)+1,nlb+ABS(mbj)+1,j+1,3)
                            IF(maa <  0.AND.mbj <  0) stemp = stemp + prefac * Waux_mat(nla+ABS(maa)+1,nlb+ABS(mbj)+1,j+1,4)
                         ENDDO
                         saba(ia+mai,jb+mbj,ka+mak) = saba(ia+mai,jb+mbj,ka+mak) +  cg_coeff(isoa1,isoa2,isoaa)*stemp
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
 
    CALL timestop(handle)

  END SUBROUTINE overlap_shg_aba

! *****************************************************************************
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
! *****************************************************************************
  SUBROUTINE dev_overlap_shg_aba(la,first_sgfa,nshella,lb,first_sgfb,nshellb,&
                             lca,first_sgfca,nshellca,cg_coeff,cg_none0_list,&
                             ncg_none0,rab,swork_cont,Waux_mat,dWaux_mat,dsaba)

    INTEGER, DIMENSION(:), INTENT(IN)        :: la, first_sgfa
    INTEGER, INTENT(IN)                      :: nshella
    INTEGER, DIMENSION(:), INTENT(IN)        :: lb, first_sgfb
    INTEGER, INTENT(IN)                      :: nshellb
    INTEGER, DIMENSION(:), INTENT(IN)        :: lca, first_sgfca
    INTEGER, INTENT(IN)                      :: nshellca
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: cg_coeff
    INTEGER, DIMENSION(:, :, :), INTENT(IN)  :: cg_none0_list
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: ncg_none0
    REAL(KIND=dp), INTENT(IN)                :: rab(3)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: swork_cont
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: Waux_mat
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), INTENT(IN)   :: dWaux_mat
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(INOUT)                          :: dsaba

    CHARACTER(len=*), PARAMETER :: routineN = 'dev_overlap_shg_aba', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ia, il, ilist, ishella, isoa1, isoa2, isoaa, j, jb, &
      jshellb, ka, kshella, laa, labmin, lai, lak, lbj, maa, mai, mak, mbj, &
      nla, nlb, sgfa, sgfb, sgfca
    REAL(KIND=dp)                            :: dprefac, dtemp(3), prefac, &
                                                rabx2(3)

    CALL timeset(routineN,handle)

    rabx2(:) = 2.0_dp * rab

    DO ishella = 1, nshella
       lai = la(ishella)
       sgfa = first_sgfa(ishella)
       ia = sgfa + lai
       DO jshellb = 1, nshellb
          lbj = lb(jshellb)
          nlb = nsoset_pm(lbj-1)
          sgfb = first_sgfb(jshellb)
          jb = sgfb + lbj
          DO kshella = 1, nshellca
             lak = lca(kshella)
             sgfca = first_sgfca(kshella)
             ka = sgfca + lak
             DO mai = -lai,lai,1
                DO mbj=-lbj,lbj,1
                   DO mak=-lak,lak,1
                      isoa1 = indso_inv(lai,mai)
                      isoa2 = indso_inv(lak,mak)
                      DO ilist = 1, ncg_none0(isoa1,isoa2)
                         isoaa = cg_none0_list(isoa1,isoa2,ilist) 
                         laa = indso(1,isoaa) 
                         maa = indso(2,isoaa) 
                         nla = nsoset_pm(laa-1)
                         labmin = MIN(laa,lbj)
                         il = INT((lai + lak - laa)/2) 
                         dtemp = 0.0_dp
                         DO j=0,labmin
                            prefac = swork_cont(ishella,jshellb,kshella,il,laa+lbj-j+1) 
                            dprefac = swork_cont(ishella,jshellb,kshella,il,laa+lbj-j+2)
                            DO i=1,3 
                               IF(maa >= 0.AND.mbj >= 0) THEN
                                  dtemp(i) = dtemp(i) + rabx2(i) * dprefac * Waux_mat(nla+ABS(maa)+1,nlb+ABS(mbj)+1,j+1,1)&
                                                      + prefac * dWaux_mat(i,nla+ABS(maa)+1,nlb+ABS(mbj)+1,j+1,1)
                               ELSEIF(maa >= 0.AND.mbj  < 0) THEN
                                  dtemp(i) = dtemp(i) + rabx2(i) * dprefac * Waux_mat(nla+ABS(maa)+1,nlb+ABS(mbj)+1,j+1,2)&
                                                      + prefac * dWaux_mat(i,nla+ABS(maa)+1,nlb+ABS(mbj)+1,j+1,2)
                               ELSEIF(maa <  0.AND.mbj >= 0) THEN
                                  dtemp(i) = dtemp(i) + rabx2(i) * dprefac * Waux_mat(nla+ABS(maa)+1,nlb+ABS(mbj)+1,j+1,3)&
                                                      + prefac * dWaux_mat(i,nla+ABS(maa)+1,nlb+ABS(mbj)+1,j+1,3)
                               ELSEIF(maa <  0.AND.mbj <  0) THEN
                                  dtemp(i) = dtemp(i) + rabx2(i) * dprefac * Waux_mat(nla+ABS(maa)+1,nlb+ABS(mbj)+1,j+1,4)&
                                                      + prefac * dWaux_mat(i,nla+ABS(maa)+1,nlb+ABS(mbj)+1,j+1,4)
                               ENDIF
                            ENDDO
                         ENDDO
                         DO i=1,3
                            dsaba(ia+mai,jb+mbj,ka+mak,i) = dsaba(ia+mai,jb+mbj,ka+mak,i) &
                                                            + cg_coeff(isoa1,isoa2,isoaa)*dtemp(i)
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

! *****************************************************************************
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
! *****************************************************************************
  SUBROUTINE overlap_shg_abb(la,first_sgfa,nshella,lb,first_sgfb,nshellb,&
                             lcb,first_sgfcb,nshellcb,cg_coeff,cg_none0_list,&
                             ncg_none0,swork_cont,Waux_mat,sabb)

    INTEGER, DIMENSION(:), INTENT(IN)        :: la, first_sgfa
    INTEGER, INTENT(IN)                      :: nshella
    INTEGER, DIMENSION(:), INTENT(IN)        :: lb, first_sgfb
    INTEGER, INTENT(IN)                      :: nshellb
    INTEGER, DIMENSION(:), INTENT(IN)        :: lcb, first_sgfcb
    INTEGER, INTENT(IN)                      :: nshellcb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: cg_coeff
    INTEGER, DIMENSION(:, :, :), INTENT(IN)  :: cg_none0_list
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: ncg_none0
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: swork_cont
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: Waux_mat
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: sabb

    CHARACTER(len=*), PARAMETER :: routineN = 'overlap_shg_abb', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ia, il, ilist, ishella, isob1, isob2, isobb, j, jb, &
      jshellb, kb, kshellb, labmin, lai, lbb, lbj, lbk, mai, mbb, mbj, mbk, &
      nla, nlb, sgfa, sgfb, sgfcb
    REAL(KIND=dp)                            :: prefac, stemp, tsign

    CALL timeset(routineN,handle)

    DO ishella = 1, nshella
       lai = la(ishella)
       nla = nsoset_pm(lai-1) + 1
       sgfa = first_sgfa(ishella)
       ia = sgfa + lai
       DO jshellb = 1, nshellb
          lbj = lb(jshellb)
          nlb = nsoset_pm(lbj-1)
          sgfb = first_sgfb(jshellb)
          jb = sgfb + lbj
          DO kshellb = 1, nshellcb
             lbk = lcb(kshellb)
             sgfcb = first_sgfcb(kshellb)
             kb = sgfcb + lbk
             DO mai = -lai,lai,1
                DO mbj=-lbj,lbj,1
                   DO mbk=-lbk,lbk,1
                      isob1 = indso_inv(lbj,mbj)
                      isob2 = indso_inv(lbk,mbk)
                      DO ilist = 1, ncg_none0(isob1,isob2)
                         isobb = cg_none0_list(isob1,isob2,ilist) 
                         lbb = indso(1,isobb) 
                         mbb = indso(2,isobb) 
                         nlb = nsoset_pm(lbb-1) + 1
                         tsign = 1.0_dp
                         IF(MODULO(lbb-lai,2) /= 0 ) tsign = -1.0_dp  
                         labmin = MIN(lai,lbb)
                         il = INT((lbj + lbk - lbb)/2) 
                         stemp = 0.0_dp
                         DO j=0,labmin
                            prefac = swork_cont(ishella,jshellb,kshellb,il,lai+lbb-j+1) 
                            IF(mai >= 0 .AND. mbb >= 0) stemp = stemp + prefac * Waux_mat(nlb+ABS(mbb),nla+ABS(mai),j+1,1)
                            ! take the transpose, since auxmat calculated for (lbb,lai)
                            IF(mai >= 0.AND.mbb  < 0) stemp = stemp + prefac * Waux_mat(nlb+ABS(mbb),nla+ABS(mai),j+1,3)
                            IF(mai <  0.AND.mbb >= 0) stemp = stemp + prefac * Waux_mat(nlb+ABS(mbb),nla+ABS(mai),j+1,2)
                            IF(mai <  0.AND.mbb <  0) stemp = stemp + prefac * Waux_mat(nlb+ABS(mbb),nla+ABS(mai),j+1,4)
                         ENDDO
                         sabb(ia+mai,jb+mbj,kb+mbk) = sabb(ia+mai,jb+mbj,kb+mbk) +  tsign*cg_coeff(isob1,isob2,isobb)*stemp
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
 
    CALL timestop(handle)

  END SUBROUTINE overlap_shg_abb

! *****************************************************************************
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
! *****************************************************************************
  SUBROUTINE dev_overlap_shg_abb(la,first_sgfa,nshella,lb,first_sgfb,nshellb,&
                                 lcb,first_sgfcb,nshellcb,cg_coeff,cg_none0_list,&
                                 ncg_none0,rab,swork_cont,Waux_mat,dWaux_mat,dsabb)

    INTEGER, DIMENSION(:), INTENT(IN)        :: la, first_sgfa
    INTEGER, INTENT(IN)                      :: nshella
    INTEGER, DIMENSION(:), INTENT(IN)        :: lb, first_sgfb
    INTEGER, INTENT(IN)                      :: nshellb
    INTEGER, DIMENSION(:), INTENT(IN)        :: lcb, first_sgfcb
    INTEGER, INTENT(IN)                      :: nshellcb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: cg_coeff
    INTEGER, DIMENSION(:, :, :), INTENT(IN)  :: cg_none0_list
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: ncg_none0
    REAL(KIND=dp), INTENT(IN)                :: rab(3)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: swork_cont
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: Waux_mat
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), INTENT(IN)   :: dWaux_mat
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(INOUT)                          :: dsabb

    CHARACTER(len=*), PARAMETER :: routineN = 'dev_overlap_shg_abb', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ia, il, ilist, ishella, isob1, isob2, isobb, j, jb, &
      jshellb, kb, kshellb, labmin, lai, lbb, lbj, lbk, mai, mbb, mbj, mbk, &
      nla, nlb, sgfa, sgfb, sgfcb
    REAL(KIND=dp)                            :: dprefac, dtemp(3), prefac, &
                                                rabx2(3), tsign

    CALL timeset(routineN,handle)

    rabx2(:) = 2.0_dp * rab

    DO ishella = 1, nshella
       lai = la(ishella)
       nla = nsoset_pm(lai-1) + 1
       sgfa = first_sgfa(ishella)
       ia = sgfa + lai
       DO jshellb = 1, nshellb
          lbj = lb(jshellb)
          nlb = nsoset_pm(lbj-1)
          sgfb = first_sgfb(jshellb)
          jb = sgfb + lbj
          DO kshellb = 1, nshellcb
             lbk = lcb(kshellb)
             sgfcb = first_sgfcb(kshellb)
             kb = sgfcb + lbk
             DO mai = -lai,lai,1
                DO mbj=-lbj,lbj,1
                   DO mbk=-lbk,lbk,1
                      isob1 = indso_inv(lbj,mbj)
                      isob2 = indso_inv(lbk,mbk)
                      DO ilist = 1, ncg_none0(isob1,isob2)
                         isobb = cg_none0_list(isob1,isob2,ilist) 
                         lbb = indso(1,isobb) 
                         mbb = indso(2,isobb) 
                         nlb = nsoset_pm(lbb-1) + 1
                         tsign = 1.0_dp
                         IF(MODULO(lbb-lai,2) /= 0 ) tsign = -1.0_dp  
                         labmin = MIN(lai,lbb)
                         il = INT((lbj + lbk - lbb)/2) 
                         dtemp = 0.0_dp
                         DO j=0,labmin
                            prefac = swork_cont(ishella,jshellb,kshellb,il,lai+lbb-j+1) 
                            dprefac = swork_cont(ishella,jshellb,kshellb,il,lai+lbb-j+2) 
                            DO i=1,3 
                               IF(mai >= 0 .AND. mbb >= 0) THEN
                                 dtemp(i) = dtemp(i) + rabx2(i) * dprefac * Waux_mat(nlb+ABS(mbb),nla+ABS(mai),j+1,1)&
                                                     + prefac * dWaux_mat(i,nlb+ABS(mbb),nla+ABS(mai),j+1,1)
                               ENDIF
                               ! take the transpose, since auxmat calculated for (lbb,lai)
                               IF(mai >= 0.AND.mbb  < 0) THEN
                                 dtemp(i) = dtemp(i) + rabx2(i) * dprefac * Waux_mat(nlb+ABS(mbb),nla+ABS(mai),j+1,3)&
                                                     + prefac *  dWaux_mat(i,nlb+ABS(mbb),nla+ABS(mai),j+1,3)
                               ENDIF
                               IF(mai <  0.AND.mbb >= 0) THEN
                                 dtemp(i) = dtemp(i) + rabx2(i) * dprefac * Waux_mat(nlb+ABS(mbb),nla+ABS(mai),j+1,2)&
                                                     + prefac * dWaux_mat(i,nlb+ABS(mbb),nla+ABS(mai),j+1,2)
                               ENDIF
                               IF(mai <  0.AND.mbb <  0) THEN
                                 dtemp(i) = dtemp(i) + rabx2(i) * dprefac * Waux_mat(nlb+ABS(mbb),nla+ABS(mai),j+1,4)&
                                                     + prefac * dWaux_mat(i,nlb+ABS(mbb),nla+ABS(mai),j+1,4)
                               ENDIF 
                            ENDDO
                         ENDDO
                         DO i=1,3 
                            dsabb(ia+mai,jb+mbj,kb+mbk,i) = dsabb(ia+mai,jb+mbj,kb+mbk,i) &
                                                            +  tsign*cg_coeff(isob1,isob2,isobb)*dtemp(i)
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

END MODULE ai_shg_overlap


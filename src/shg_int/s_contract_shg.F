!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for calculating the s-integrals and their scalar derivatives with respect to rab2
!>        over solid harmonic Gaussian (SHG) functions + contraction routines for these integrals
!>        i)  (s|O(r12)|s) where O(r12) is the overlap, coulomb operator etc.
!>        ii) (aba) and (abb) s-overlaps
!> \par Literature (partly)
!>      T.J. Giese and D. M. York, J. Chem. Phys, 128, 064104 (2008)
!>      T. Helgaker, P Joergensen, J. Olsen, Molecular Electronic-Structure
!>                                           Theory, Wiley
!>      R. Ahlrichs, PCCP, 8, 3072 (2006)
!> \par History
!>      created [05.2016]
!> \author Dorothea Golze
! **************************************************************************************************
MODULE s_contract_shg
   USE gamma,                           ONLY: fgamma=>fgamma_0
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: dfac,&
                                              fac,&
                                              pi
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 's_contract_shg'

! *** Public subroutines ***
   PUBLIC :: s_overlap_ab, s_overlap_abb, s_coulomb_ab, s_verf_ab, s_verfc_ab, &
             s_vgauss_ab, s_gauss_ab, s_ra2m_ab

   PUBLIC :: contract_sint_ab_clow, contract_sint_ab_chigh, contract_s_overlap_aba, &
             contract_s_overlap_abb, contract_s_ra2m_ab

CONTAINS

! **************************************************************************************************
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
! **************************************************************************************************
   SUBROUTINE s_overlap_ab(la_max, npgfa, zeta, lb_max, npgfb, zetb, rab, s, calculate_forces)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: s
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 's_overlap_ab', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ids, ipgfa, jpgfb, ndev
      REAL(KIND=dp)                                      :: a, b, rab2, xhi, zet

      CALL timeset(routineN, handle)

      ! Distance of the centers a and b
      rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
      ndev = 0
      IF (calculate_forces) ndev = 1
      ! Loops over all pairs of primitive Gaussian-type functions
      DO ipgfa = 1, npgfa
         DO jpgfb = 1, npgfb

            ! Distance Screening   !maybe later

            ! Calculate some prefactors
            a = zeta(ipgfa)
            b = zetb(jpgfb)
            zet = a+b
            xhi = a*b/zet

            ! [s|s] integral
            s(ipgfa, jpgfb, 1) = (pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

            DO ids = 2, la_max+lb_max+ndev+1
               s(ipgfa, jpgfb, ids) = -xhi*s(ipgfa, jpgfb, ids-1)
            ENDDO

         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE s_overlap_ab

! **************************************************************************************************
!> \brief calculates [s|ra^n|s] integrals for [aba] and the [s|rb^n|s]
!>        integrals for [abb]
!> \param la_max maximal l quantum number on a, orbital basis
!> \param npgfa number of primitive Gaussian on a, orbital basis
!> \param zeta set of exponents on a, orbital basis
!> \param lb_max maximal l quantum number on b, orbital basis
!> \param npgfb number of primitive Gaussian on a, orbital basis
!> \param zetb set of exponents on b, orbital basis
!> \param lcb_max maximal l quantum number of aux basis on b
!> \param npgfcb number of primitive Gaussian on b.  aux basis
!> \param zetcb set of exponents on b, aux basis
!> \param rab distance vector between a and b
!> \param s uncontracted [s|r^n|s] integrals
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE s_overlap_abb(la_max, npgfa, zeta, lb_max, npgfb, zetb, lcb_max, npgfcb, zetcb, &
                            rab, s, calculate_forces)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb
      INTEGER, INTENT(IN)                                :: lcb_max, npgfcb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetcb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :, :, :, :), &
         INTENT(INOUT)                                   :: s
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 's_overlap_abb', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ids, il, ipgfa, j, jpgfb, &
                                                            kpgfb, lbb_max, lmax, n, ndev, nds, &
                                                            nfac, nl
      REAL(KIND=dp)                                      :: a, b, dfsr_int, exp_rab2, k, pfac, &
                                                            prefac, rab2, sqrt_pi3, sqrt_zet, &
                                                            sr_int, temp, xhi, zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dsr_int, dtemp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: coeff_srs

      CALL timeset(routineN, handle)

      ! Distance of the centers a and b
      rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
      ndev = 0
      IF (calculate_forces) ndev = 1

      lbb_max = lb_max+lcb_max
      nl = INT(lbb_max/2)
      IF (lb_max == 0 .OR. lcb_max == 0) nl = 0
      lmax = la_max+lbb_max

      ALLOCATE (dtemp(nl+1), dsr_int(nl+1))
      ALLOCATE (coeff_srs(nl+1, nl+1, nl+1))
      IF (nl > 5) CALL get_prefac_sabb(nl, coeff_srs)

      sqrt_pi3 = SQRT(pi**3)

      ! Loops over all pairs of primitive Gaussian-type functions
      DO kpgfb = 1, npgfcb
         DO jpgfb = 1, npgfb
            DO ipgfa = 1, npgfa

               !Calculate some prefactors
               a = zeta(ipgfa)
               b = zetb(jpgfb)+zetcb(kpgfb)

               zet = a+b
               xhi = a*b/zet
               exp_rab2 = EXP(-xhi*rab2)

               pfac = a**2/zet
               sqrt_zet = SQRT(zet)

               DO il = 0, nl
                  nds = lmax-2*il+ndev+1
                  SELECT CASE (il)
                  CASE (0)
                     ! [s|s] integral
                     s(ipgfa, jpgfb, kpgfb, il+1, 1) = (pi/zet)**(1.5_dp)*exp_rab2
                     DO ids = 2, nds
                        n = ids-1
                        s(ipgfa, jpgfb, kpgfb, il+1, ids) = (-xhi)**n*s(ipgfa, jpgfb, kpgfb, il+1, 1)
                     ENDDO
                  CASE (1)
                     ![s|r^2|s] integral
                     sr_int = sqrt_pi3/sqrt_zet**5*(3.0_dp+2.0_dp*pfac*rab2)/2.0_dp
                     s(ipgfa, jpgfb, kpgfb, il+1, 1) = exp_rab2*sr_int
                     k = sqrt_pi3*a**2/sqrt_zet**7
                     DO ids = 2, nds
                        n = ids-1
                        s(ipgfa, jpgfb, kpgfb, il+1, ids) = (-xhi)**n*exp_rab2*sr_int &
                                                            +n*(-xhi)**(n-1)*k*exp_rab2
                     ENDDO
                  CASE (2)
                     ![s|r^4|s] integral
                     prefac = sqrt_pi3/4.0_dp/sqrt_zet**7
                     temp = 15.0_dp+20.0_dp*pfac*rab2+4.0_dp*(pfac*rab2)**2
                     sr_int = prefac*temp
                     s(ipgfa, jpgfb, kpgfb, il+1, 1) = exp_rab2*sr_int
                     !** derivatives
                     k = sqrt_pi3*a**4/sqrt_zet**11
                     dsr_int(1) = prefac*(20.0_dp*pfac+8.0_dp*pfac**2*rab2)
                     DO ids = 2, nds
                        n = ids-1
                        dtemp(1) = (-xhi)**n*exp_rab2*sr_int
                        dtemp(2) = n*(-xhi)**(n-1)*exp_rab2*dsr_int(1)
                        dtemp(3) = (n**2-n)*(-xhi)**(n-2)*k*exp_rab2
                        s(ipgfa, jpgfb, kpgfb, il+1, ids) = dtemp(1)+dtemp(2)+dtemp(3)
                     ENDDO
                  CASE (3)
                     ![s|r^6|s] integral
                     prefac = sqrt_pi3/8.0_dp/sqrt_zet**9
                     temp = 105.0_dp+210.0_dp*pfac*rab2
                     temp = temp+84.0_dp*(pfac*rab2)**2+8.0_dp*(pfac*rab2)**3
                     sr_int = prefac*temp
                     s(ipgfa, jpgfb, kpgfb, il+1, 1) = exp_rab2*sr_int
                     !** derivatives
                     k = sqrt_pi3*a**6/sqrt_zet**15
                     dsr_int(1) = prefac*(210.0_dp*pfac+168.0_dp*pfac**2*rab2 &
                                          +24.0_dp*pfac**3*rab2**2)
                     dsr_int(2) = prefac*(168.0_dp*pfac**2+48.0_dp*pfac**3*rab2)
                     DO ids = 2, nds
                        n = ids-1
                        dtemp(1) = (-xhi)**n*exp_rab2*sr_int
                        dtemp(2) = REAL(n, dp)*(-xhi)**(n-1)*exp_rab2*dsr_int(1)
                        dtemp(3) = REAL(n**2-n, dp)/2.0_dp*(-xhi)**(n-2) &
                                   *exp_rab2*dsr_int(2)
                        dtemp(4) = REAL(n*(n-1)*(n-2), dp)*(-xhi)**(n-3)*k*exp_rab2
                        s(ipgfa, jpgfb, kpgfb, il+1, ids) = dtemp(1)+dtemp(2) &
                                                            +dtemp(3)+dtemp(4)
                     ENDDO
                  CASE (4)
                     ![s|r^8|s] integral
                     prefac = sqrt_pi3/16.0_dp/sqrt_zet**11
                     temp = 945.0_dp+2520.0_dp*pfac*rab2+1512.0_dp*(pfac*rab2)**2
                     temp = temp+288.0_dp*(pfac*rab2)**3+16.0_dp*(pfac*rab2)**4
                     sr_int = prefac*temp
                     s(ipgfa, jpgfb, kpgfb, il+1, 1) = exp_rab2*sr_int
                     !** derivatives
                     k = sqrt_pi3*a**8/sqrt_zet**19
                     dsr_int(1) = 2520.0_dp*pfac+3024.0_dp*pfac**2*rab2
                     dsr_int(1) = dsr_int(1)+864.0_dp*pfac**3*rab2**2 &
                                  +64.0_dp*pfac**4*rab2**3
                     dsr_int(1) = prefac*dsr_int(1)
                     dsr_int(2) = 3024.0_dp*pfac**2+1728.0_dp*pfac**3*rab2
                     dsr_int(2) = dsr_int(2)+192.0_dp*pfac**4*rab2**2
                     dsr_int(2) = prefac*dsr_int(2)
                     dsr_int(3) = 1728.0_dp*pfac**3+384.0_dp*pfac**4*rab2
                     dsr_int(3) = prefac*dsr_int(3)
                     DO ids = 2, nds
                        n = ids-1
                        dtemp(1) = (-xhi)**n*exp_rab2*sr_int
                        dtemp(2) = REAL(n, dp)*(-xhi)**(n-1)*exp_rab2*dsr_int(1)
                        dtemp(3) = REAL(n**2-n, dp)/2.0_dp*(-xhi)**(n-2) &
                                   *exp_rab2*dsr_int(2)
                        dtemp(4) = REAL(n*(n-1)*(n-2), dp)/6.0_dp*(-xhi)**(n-3) &
                                   *exp_rab2*dsr_int(3)
                        dtemp(5) = REAL(n*(n-1)*(n-2)*(n-3), dp)*(-xhi)**(n-4) &
                                   *k*exp_rab2
                        s(ipgfa, jpgfb, kpgfb, il+1, ids) = dtemp(1)+dtemp(2)+dtemp(3) &
                                                            +dtemp(4)+dtemp(5)
                     ENDDO
                  CASE (5)
                     ![s|r^10|s] integral
                     prefac = sqrt_pi3/32.0_dp/sqrt_zet**13
                     temp = 10395.0_dp+34650.0_dp*pfac*rab2
                     temp = temp+27720.0_dp*(pfac*rab2)**2+7920.0_dp*(pfac*rab2)**3
                     temp = temp+880.0_dp*(pfac*rab2)**4+32.0_dp*(pfac*rab2)**5
                     sr_int = prefac*temp
                     s(ipgfa, jpgfb, kpgfb, il+1, 1) = exp_rab2*sr_int
                     !** derivatives
                     k = sqrt_pi3*a**10/sqrt_zet**23
                     dsr_int(1) = 34650.0_dp*pfac+55440.0_dp*pfac**2*rab2
                     dsr_int(1) = dsr_int(1)+23760.0_dp*pfac**3*rab2**2
                     dsr_int(1) = dsr_int(1)+3520.0_dp*pfac**4*rab2**3
                     dsr_int(1) = dsr_int(1)+160.0_dp*pfac**5*rab2**4
                     dsr_int(1) = prefac*dsr_int(1)
                     dsr_int(2) = 55440.0_dp*pfac**2+47520.0_dp*pfac**3*rab2
                     dsr_int(2) = dsr_int(2)+10560.0_dp*pfac**4*rab2**2
                     dsr_int(2) = dsr_int(2)+640.0_dp*pfac**5*rab2**3
                     dsr_int(2) = prefac*dsr_int(2)
                     dsr_int(3) = 47520.0_dp*pfac**3+21120.0_dp*pfac**4*rab2
                     dsr_int(3) = dsr_int(3)+1920.0_dp*pfac**5*rab2**2
                     dsr_int(3) = prefac*dsr_int(3)
                     dsr_int(4) = 21120.0_dp*pfac**4+3840.0_dp*pfac**5*rab2
                     dsr_int(4) = prefac*dsr_int(4)
                     DO ids = 2, nds
                        n = ids-1
                        dtemp(1) = (-xhi)**n*exp_rab2*sr_int
                        dtemp(2) = REAL(n, dp)*(-xhi)**(n-1)*exp_rab2*dsr_int(1)
                        dtemp(3) = REAL(n**2-n, dp)/2.0_dp*(-xhi)**(n-2) &
                                   *exp_rab2*dsr_int(2)
                        dtemp(4) = REAL(n*(n-1)*(n-2), dp)/6.0_dp*(-xhi)**(n-3) &
                                   *exp_rab2*dsr_int(3)
                        dtemp(5) = REAL(n*(n-1)*(n-2)*(n-3), dp)/24.0_dp*(-xhi)**(n-4) &
                                   *exp_rab2*dsr_int(4)
                        dtemp(6) = REAL(n*(n-1)*(n-2)*(n-3)*(n-4), dp)*(-xhi)**(n-5) &
                                   *k*exp_rab2
                        s(ipgfa, jpgfb, kpgfb, il+1, ids) = dtemp(1)+dtemp(2)+dtemp(3) &
                                                            +dtemp(4)+dtemp(5)+dtemp(6)
                     ENDDO
                  CASE DEFAULT
                     !*** general formula; factor 1.5-2 slower than explicit expressions
                     prefac = exp_rab2/sqrt_zet**(2*il+3)
                     sr_int = 0.0_dp
                     DO i = 0, il
                        sr_int = sr_int+coeff_srs(i+1, 1, il+1)*(pfac)**i*rab2**i
                     ENDDO
                     s(ipgfa, jpgfb, kpgfb, il+1, 1) = prefac*sr_int
                     DO ids = 2, nds
                        n = ids-1
                        nfac = 1
                        dfsr_int = (-xhi)**n*sr_int
                        DO j = 1, il
                           temp = 0.0_dp
                           DO i = j, il
                              temp = temp+coeff_srs(i+1, j+1, il+1)*(pfac)**i*rab2**(i-j)
                           ENDDO
                           nfac = nfac*(n-j+1)
                           dfsr_int = dfsr_int+temp*REAL(nfac, dp)/fac(j)*(-xhi)**(n-j)
                        ENDDO
                        s(ipgfa, jpgfb, kpgfb, il+1, ids) = prefac*dfsr_int
                     ENDDO

                  END SELECT

               ENDDO

            END DO
         END DO
      END DO

      DEALLOCATE (dtemp, dsr_int)
      DEALLOCATE (coeff_srs)

      CALL timestop(handle)

   END SUBROUTINE s_overlap_abb

! **************************************************************************************************
!> \brief calculates the uncontracted, not normalized [0a|ra^(2m)|0b] two-center integral,
!>        where ra = r-Ra and Ra center of a
!> \param la_max maximal l quantum number on a
!> \param npgfa number of primitive Gaussian on a
!> \param zeta set of exponents on a
!> \param lb_max maximal l quantum number on b
!> \param npgfb number of primitive Gaussian on a
!> \param zetb set of exponents on a
!> \param m exponent of the ra operator
!> \param rab distance vector between a and b
!> \param s ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE s_ra2m_ab(la_max, npgfa, zeta, lb_max, npgfb, zetb, m, rab, s, calculate_forces)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb
      INTEGER, INTENT(IN)                                :: m
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(INOUT)                                   :: s
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 's_ra2m_ab', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ids, il, ipgfa, j, jpgfb, n, &
                                                            ndev, nds, nfac
      REAL(KIND=dp)                                      :: a, b, dfsr_int, exp_rab2, k, pfac, &
                                                            prefac, rab2, sqrt_pi3, sqrt_zet, &
                                                            sr_int, temp, xhi, zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dsr_int, dtemp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: coeff_srs

      CALL timeset(routineN, handle)

      ! Distance of the centers a and b
      rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
      ndev = 0
      IF (calculate_forces) ndev = 1

      ALLOCATE (dtemp(m+1), dsr_int(m+1))
      ALLOCATE (coeff_srs(m+1, m+1, m+1))
      CALL get_prefac_sabb(m, coeff_srs)
      !IF(m > 5) CALL get_prefac_sabb(m, coeff_srs)
      sqrt_pi3 = SQRT(pi**3)

      ! Loops over all pairs of primitive Gaussian-type functions
      DO ipgfa = 1, npgfa
         DO jpgfb = 1, npgfb

            ! Calculate some prefactors
            a = zeta(ipgfa)
            b = zetb(jpgfb)
            zet = a+b
            xhi = a*b/zet
            exp_rab2 = EXP(-xhi*rab2)

            sqrt_zet = SQRT(zet)
            pfac = b**2/zet

            nds = la_max+lb_max+ndev+1
            DO il = 0, m
               SELECT CASE (il)
               CASE (0)
                  ! [s|s] integral
                  s(ipgfa, jpgfb, m-il+1, 1) = (pi/zet)**(1.5_dp)*exp_rab2
                  DO ids = 2, nds
                     n = ids-1
                     s(ipgfa, jpgfb, m-il+1, ids) = (-xhi)**n*s(ipgfa, jpgfb, m-il+1, 1)
                  ENDDO
               CASE (1)
                  ![s|r^2|s] integral
                  sr_int = sqrt_pi3/sqrt_zet**5*(3.0_dp+2.0_dp*pfac*rab2)/2.0_dp
                  s(ipgfa, jpgfb, m-il+1, 1) = exp_rab2*sr_int
                  k = sqrt_pi3*b**2/sqrt_zet**7
                  DO ids = 2, nds
                     n = ids-1
                     s(ipgfa, jpgfb, m-il+1, ids) = (-xhi)**n*exp_rab2*sr_int &
                                                    +n*(-xhi)**(n-1)*k*exp_rab2
                  ENDDO
               CASE (2)
                  ![s|r^4|s] integral
                  prefac = sqrt_pi3/4.0_dp/sqrt_zet**7
                  temp = 15.0_dp+20.0_dp*pfac*rab2+4.0_dp*(pfac*rab2)**2
                  sr_int = prefac*temp
                  s(ipgfa, jpgfb, m-il+1, 1) = exp_rab2*sr_int
                  !** derivatives
                  k = sqrt_pi3*b**4/sqrt_zet**11
                  dsr_int(1) = prefac*(20.0_dp*pfac+8.0_dp*pfac**2*rab2)
                  DO ids = 2, nds
                     n = ids-1
                     dtemp(1) = (-xhi)**n*exp_rab2*sr_int
                     dtemp(2) = n*(-xhi)**(n-1)*exp_rab2*dsr_int(1)
                     dtemp(3) = (n**2-n)*(-xhi)**(n-2)*k*exp_rab2
                     s(ipgfa, jpgfb, m-il+1, ids) = dtemp(1)+dtemp(2)+dtemp(3)
                  ENDDO
               CASE (3)
                  ![s|r^6|s] integral
                  prefac = sqrt_pi3/8.0_dp/sqrt_zet**9
                  temp = 105.0_dp+210.0_dp*pfac*rab2
                  temp = temp+84.0_dp*(pfac*rab2)**2+8.0_dp*(pfac*rab2)**3
                  sr_int = prefac*temp
                  s(ipgfa, jpgfb, m-il+1, 1) = exp_rab2*sr_int
                  !** derivatives
                  k = sqrt_pi3*b**6/sqrt_zet**15
                  dsr_int(1) = prefac*(210.0_dp*pfac+168.0_dp*pfac**2*rab2 &
                                       +24.0_dp*pfac**3*rab2**2)
                  dsr_int(2) = prefac*(168.0_dp*pfac**2+48.0_dp*pfac**3*rab2)
                  DO ids = 2, nds
                     n = ids-1
                     dtemp(1) = (-xhi)**n*exp_rab2*sr_int
                     dtemp(2) = REAL(n, dp)*(-xhi)**(n-1)*exp_rab2*dsr_int(1)
                     dtemp(3) = REAL(n**2-n, dp)/2.0_dp*(-xhi)**(n-2) &
                                *exp_rab2*dsr_int(2)
                     dtemp(4) = REAL(n*(n-1)*(n-2), dp)*(-xhi)**(n-3)*k*exp_rab2
                     s(ipgfa, jpgfb, m-il+1, ids) = dtemp(1)+dtemp(2) &
                                                    +dtemp(3)+dtemp(4)
                  ENDDO
               CASE (4)
                  ![s|r^8|s] integral
                  prefac = sqrt_pi3/16.0_dp/sqrt_zet**11
                  temp = 945.0_dp+2520.0_dp*pfac*rab2+1512.0_dp*(pfac*rab2)**2
                  temp = temp+288.0_dp*(pfac*rab2)**3+16.0_dp*(pfac*rab2)**4
                  sr_int = prefac*temp
                  s(ipgfa, jpgfb, m-il+1, 1) = exp_rab2*sr_int
                  !** derivatives
                  k = sqrt_pi3*b**8/sqrt_zet**19
                  dsr_int(1) = 2520.0_dp*pfac+3024.0_dp*pfac**2*rab2
                  dsr_int(1) = dsr_int(1)+864.0_dp*pfac**3*rab2**2 &
                               +64.0_dp*pfac**4*rab2**3
                  dsr_int(1) = prefac*dsr_int(1)
                  dsr_int(2) = 3024.0_dp*pfac**2+1728.0_dp*pfac**3*rab2
                  dsr_int(2) = dsr_int(2)+192.0_dp*pfac**4*rab2**2
                  dsr_int(2) = prefac*dsr_int(2)
                  dsr_int(3) = 1728.0_dp*pfac**3+384.0_dp*pfac**4*rab2
                  dsr_int(3) = prefac*dsr_int(3)
                  DO ids = 2, nds
                     n = ids-1
                     dtemp(1) = (-xhi)**n*exp_rab2*sr_int
                     dtemp(2) = REAL(n, dp)*(-xhi)**(n-1)*exp_rab2*dsr_int(1)
                     dtemp(3) = REAL(n**2-n, dp)/2.0_dp*(-xhi)**(n-2) &
                                *exp_rab2*dsr_int(2)
                     dtemp(4) = REAL(n*(n-1)*(n-2), dp)/6.0_dp*(-xhi)**(n-3) &
                                *exp_rab2*dsr_int(3)
                     dtemp(5) = REAL(n*(n-1)*(n-2)*(n-3), dp)*(-xhi)**(n-4) &
                                *k*exp_rab2
                     s(ipgfa, jpgfb, m-il+1, ids) = dtemp(1)+dtemp(2)+dtemp(3) &
                                                    +dtemp(4)+dtemp(5)
                  ENDDO
               CASE (5)
                  ![s|r^10|s] integral
                  prefac = sqrt_pi3/32.0_dp/sqrt_zet**13
                  temp = 10395.0_dp+34650.0_dp*pfac*rab2
                  temp = temp+27720.0_dp*(pfac*rab2)**2+7920.0_dp*(pfac*rab2)**3
                  temp = temp+880.0_dp*(pfac*rab2)**4+32.0_dp*(pfac*rab2)**5
                  sr_int = prefac*temp
                  s(ipgfa, jpgfb, m-il+1, 1) = exp_rab2*sr_int
                  !** derivatives
                  k = sqrt_pi3*b**10/sqrt_zet**23
                  dsr_int(1) = 34650.0_dp*pfac+55440.0_dp*pfac**2*rab2
                  dsr_int(1) = dsr_int(1)+23760.0_dp*pfac**3*rab2**2
                  dsr_int(1) = dsr_int(1)+3520.0_dp*pfac**4*rab2**3
                  dsr_int(1) = dsr_int(1)+160.0_dp*pfac**5*rab2**4
                  dsr_int(1) = prefac*dsr_int(1)
                  dsr_int(2) = 55440.0_dp*pfac**2+47520.0_dp*pfac**3*rab2
                  dsr_int(2) = dsr_int(2)+10560.0_dp*pfac**4*rab2**2
                  dsr_int(2) = dsr_int(2)+640.0_dp*pfac**5*rab2**3
                  dsr_int(2) = prefac*dsr_int(2)
                  dsr_int(3) = 47520.0_dp*pfac**3+21120.0_dp*pfac**4*rab2
                  dsr_int(3) = dsr_int(3)+1920.0_dp*pfac**5*rab2**2
                  dsr_int(3) = prefac*dsr_int(3)
                  dsr_int(4) = 21120.0_dp*pfac**4+3840.0_dp*pfac**5*rab2
                  dsr_int(4) = prefac*dsr_int(4)
                  DO ids = 2, nds
                     n = ids-1
                     dtemp(1) = (-xhi)**n*exp_rab2*sr_int
                     dtemp(2) = REAL(n, dp)*(-xhi)**(n-1)*exp_rab2*dsr_int(1)
                     dtemp(3) = REAL(n**2-n, dp)/2.0_dp*(-xhi)**(n-2) &
                                *exp_rab2*dsr_int(2)
                     dtemp(4) = REAL(n*(n-1)*(n-2), dp)/6.0_dp*(-xhi)**(n-3) &
                                *exp_rab2*dsr_int(3)
                     dtemp(5) = REAL(n*(n-1)*(n-2)*(n-3), dp)/24.0_dp*(-xhi)**(n-4) &
                                *exp_rab2*dsr_int(4)
                     dtemp(6) = REAL(n*(n-1)*(n-2)*(n-3)*(n-4), dp)*(-xhi)**(n-5) &
                                *k*exp_rab2
                     s(ipgfa, jpgfb, m-il+1, ids) = dtemp(1)+dtemp(2)+dtemp(3) &
                                                    +dtemp(4)+dtemp(5)+dtemp(6)
                  ENDDO
               CASE DEFAULT
                  prefac = exp_rab2/sqrt_zet**(2*il+3)
                  sr_int = 0.0_dp
                  DO i = 0, il
                     sr_int = sr_int+coeff_srs(i+1, 1, il+1)*(pfac)**i*rab2**i
                  ENDDO
                  s(ipgfa, jpgfb, m-il+1, 1) = prefac*sr_int
                  DO ids = 2, nds
                     n = ids-1
                     nfac = 1
                     dfsr_int = (-xhi)**n*sr_int
                     DO j = 1, il
                        temp = 0.0_dp
                        DO i = j, il
                           temp = temp+coeff_srs(i+1, j+1, il+1)*(pfac)**i*rab2**(i-j)
                        ENDDO
                        nfac = nfac*(n-j+1)
                        dfsr_int = dfsr_int+temp*REAL(nfac, dp)/fac(j)*(-xhi)**(n-j)
                     ENDDO
                     s(ipgfa, jpgfb, m-il+1, ids) = prefac*dfsr_int
                  ENDDO
               END SELECT
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE (coeff_srs)
      DEALLOCATE (dtemp, dsr_int)
      CALL timestop(handle)

   END SUBROUTINE s_ra2m_ab

! **************************************************************************************************
!> \brief prefactor for the general formula to calculate the (0a|0b|0b) overlap integrals
!> \param nl ...
!> \param prefac ...
! **************************************************************************************************
   SUBROUTINE get_prefac_sabb(nl, prefac)
      INTEGER, INTENT(IN)                                :: nl
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: prefac

      CHARACTER(len=*), PARAMETER :: routineN = 'get_prefac_sabb', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: il, j, k
      REAL(KIND=dp)                                      :: sqrt_pi3, temp

      sqrt_pi3 = SQRT(pi**3)

      DO il = 0, nl
         temp = dfac(2*il+1)*sqrt_pi3*fac(il)/2.0_dp**il
         DO j = 0, il
            DO k = j, il
               prefac(k+1, j+1, il+1) = temp*2.0_dp**k/dfac(2*k+1)/fac(il-k)/fac(k-j)
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE get_prefac_sabb

! **************************************************************************************************
!> \brief calculates the uncontracted, not normalized [s|1/r12|s] two-center coulomb integral
!> \param la_max maximal l quantum number on a
!> \param npgfa number of primitive Gaussian on a
!> \param zeta set of exponents on a
!> \param lb_max maximal l quantum number on b
!> \param npgfb number of primitive Gaussian on a
!> \param zetb set of exponents on a
!> \param omega parameter not needed, but given for the sake of the abstract interface
!> \param rab distance vector between a and b
!> \param v uncontracted coulomb integral of s functions
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE s_coulomb_ab(la_max, npgfa, zeta, lb_max, npgfb, zetb, omega, rab, v, calculate_forces)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 's_coulomb_ab', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ids, ipgfa, jpgfb, n, ndev, nmax
      REAL(KIND=dp)                                      :: a, b, dummy, prefac, rab2, T, xhi, zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: f

      CALL timeset(routineN, handle)

      dummy = omega

      ! Distance of the centers a and b
      rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
      ndev = 0
      IF (calculate_forces) ndev = 1
      nmax = la_max+lb_max+ndev+1
      ALLOCATE (f(0:nmax))
      ! Loops over all pairs of primitive Gaussian-type functions
      DO ipgfa = 1, npgfa
         DO jpgfb = 1, npgfb

            ! Calculate some prefactors
            a = zeta(ipgfa)
            b = zetb(jpgfb)
            zet = a+b
            xhi = a*b/zet
            prefac = 2.0_dp*SQRT(pi**5)/(a*b)/SQRT(zet)
            T = xhi*rab2
            CALL fgamma(nmax-1, T, f)

            DO ids = 1, nmax
               n = ids-1
               v(ipgfa, jpgfb, ids) = prefac*(-xhi)**n*f(n)
            ENDDO

         END DO
      END DO
      DEALLOCATE (f)

      CALL timestop(handle)

   END SUBROUTINE s_coulomb_ab

! **************************************************************************************************
!> \brief calculates the uncontracted, not normalized [s|1/erf(omega*r12)/r12|s] two-center integral
!> \param la_max maximal l quantum number on a
!> \param npgfa number of primitive Gaussian on a
!> \param zeta set of exponents on a
!> \param lb_max maximal l quantum number on b
!> \param npgfb number of primitive Gaussian on a
!> \param zetb set of exponents on a
!> \param omega parameter in the operator
!> \param rab distance vector between a and b
!> \param v uncontracted erf(r)/r integral of s functions
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE s_verf_ab(la_max, npgfa, zeta, lb_max, npgfb, zetb, omega, rab, v, calculate_forces)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 's_verf_ab', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ids, ipgfa, jpgfb, n, ndev, nmax
      REAL(KIND=dp)                                      :: a, Arg, b, comega, prefac, rab2, T, xhi, &
                                                            zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: f

      CALL timeset(routineN, handle)

      ! Distance of the centers a and b
      rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
      ndev = 0
      IF (calculate_forces) ndev = 1
      nmax = la_max+lb_max+ndev+1
      ALLOCATE (f(0:nmax))
      ! Loops over all pairs of primitive Gaussian-type functions
      DO ipgfa = 1, npgfa
         DO jpgfb = 1, npgfb

            ! Calculate some prefactors
            a = zeta(ipgfa)
            b = zetb(jpgfb)
            zet = a+b
            xhi = a*b/zet
            comega = omega**2/(omega**2+xhi)
            prefac = 2.0_dp*SQRT(pi**5)*SQRT(comega)/(a*b)/SQRT(zet)
            T = xhi*rab2
            Arg = comega*T
            CALL fgamma(nmax-1, Arg, f)

            DO ids = 1, nmax
               n = ids-1
               v(ipgfa, jpgfb, ids) = prefac*(-xhi*comega)**n*f(n)
            ENDDO

         END DO
      END DO
      DEALLOCATE (f)

      CALL timestop(handle)

   END SUBROUTINE s_verf_ab

! **************************************************************************************************
!> \brief calculates the uncontracted, not normalized [s|1/erfc(omega*r12)/r12|s] two-center
!>        integral
!> \param la_max maximal l quantum number on a
!> \param npgfa number of primitive Gaussian on a
!> \param zeta set of exponents on a
!> \param lb_max maximal l quantum number on b
!> \param npgfb number of primitive Gaussian on a
!> \param zetb set of exponents on a
!> \param omega parameter in the operator
!> \param rab distance vector between a and b
!> \param v uncontracted erf(r)/r integral of s functions
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE s_verfc_ab(la_max, npgfa, zeta, lb_max, npgfb, zetb, omega, rab, v, calculate_forces)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 's_verfc_ab', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ids, ipgfa, jpgfb, n, ndev, nmax
      REAL(KIND=dp)                                      :: a, b, comega, comegaT, prefac, rab2, T, &
                                                            xhi, zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: fv, fverf

      CALL timeset(routineN, handle)

      ! Distance of the centers a and b
      rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
      ndev = 0
      IF (calculate_forces) ndev = 1
      nmax = la_max+lb_max+ndev+1
      ALLOCATE (fv(0:nmax), fverf(0:nmax))
      ! Loops over all pairs of primitive Gaussian-type functions
      DO ipgfa = 1, npgfa
         DO jpgfb = 1, npgfb

            ! Calculate some prefactors
            a = zeta(ipgfa)
            b = zetb(jpgfb)
            zet = a+b
            xhi = a*b/zet
            comega = omega**2/(omega**2+xhi)
            prefac = 2.0_dp*SQRT(pi**5)/(a*b)/SQRT(zet)
            T = xhi*rab2
            comegaT = comega*T
            CALL fgamma(nmax-1, T, fv)
            CALL fgamma(nmax-1, comegaT, fverf)

            DO ids = 1, nmax
               n = ids-1
               v(ipgfa, jpgfb, ids) = prefac*(-xhi)**n*(fv(n)-SQRT(comega)*comega**n*fverf(n))
            ENDDO

         END DO
      END DO
      DEALLOCATE (fv, fverf)

      CALL timestop(handle)

   END SUBROUTINE s_verfc_ab

! **************************************************************************************************
!> \brief calculates the uncontracted, not normalized [s|exp(-omega*r12**2)/r12|s] two-center
!>        integral
!> \param la_max maximal l quantum number on a
!> \param npgfa number of primitive Gaussian on a
!> \param zeta set of exponents on a
!> \param lb_max maximal l quantum number on b
!> \param npgfb number of primitive Gaussian on a
!> \param zetb set of exponents on a
!> \param omega parameter in the operator
!> \param rab distance vector between a and b
!> \param v uncontracted exp(-omega*r**2)/r integral of s functions
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE s_vgauss_ab(la_max, npgfa, zeta, lb_max, npgfb, zetb, omega, rab, v, calculate_forces)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 's_vgauss_ab', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ids, ipgfa, j, jpgfb, n, ndev, &
                                                            nmax
      REAL(KIND=dp)                                      :: a, b, eta, etaT, expT, oeta, prefac, &
                                                            rab2, T, xeta, xhi, zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: f

      CALL timeset(routineN, handle)

      ! Distance of the centers a and b
      rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
      ndev = 0
      IF (calculate_forces) ndev = 1
      nmax = la_max+lb_max+ndev+1
      ALLOCATE (f(0:nmax))
      ! Loops over all pairs of primitive Gaussian-type functions
      v = 0.0_dp
      DO ipgfa = 1, npgfa
         DO jpgfb = 1, npgfb

            ! Calculate some prefactors
            a = zeta(ipgfa)
            b = zetb(jpgfb)
            zet = a+b
            xhi = a*b/zet
            eta = xhi/(xhi+omega)
            oeta = omega*eta
            xeta = xhi*eta
            T = xhi*rab2
            expT = EXP(-omega/(omega+xhi)*T)
            prefac = 2.0_dp*SQRT(pi**5/zet**3)/(xhi+omega)*expT
            etaT = eta*T
            CALL fgamma(nmax-1, etaT, f)

            DO ids = 1, nmax
               n = ids-1
               DO j = 0, n
                  v(ipgfa, jpgfb, ids) = v(ipgfa, jpgfb, ids) &
                                         +prefac*fac(n)/fac(j)/fac(n-j)*(-oeta)**(n-j)*(-xeta)**j*f(j)
               ENDDO
            ENDDO

         END DO
      END DO
      DEALLOCATE (f)

      CALL timestop(handle)

   END SUBROUTINE s_vgauss_ab

! **************************************************************************************************
!> \brief calculates the uncontracted, not normalized [s|exp(-omega*r12**2)|s] two-center
!>        integral
!> \param la_max maximal l quantum number on a
!> \param npgfa number of primitive Gaussian on a
!> \param zeta set of exponents on a
!> \param lb_max maximal l quantum number on b
!> \param npgfb number of primitive Gaussian on a
!> \param zetb set of exponents on a
!> \param omega parameter in the operator
!> \param rab distance vector between a and b
!> \param v uncontracted exp(-omega*r**2) integral of s functions
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE s_gauss_ab(la_max, npgfa, zeta, lb_max, npgfb, zetb, omega, rab, v, calculate_forces)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 's_gauss_ab', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ids, ipgfa, jpgfb, n, ndev, nmax
      REAL(KIND=dp)                                      :: a, b, eta, expT, oeta, prefac, rab2, T, &
                                                            xhi, zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: f

      CALL timeset(routineN, handle)

      ! Distance of the centers a and b
      rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
      ndev = 0
      IF (calculate_forces) ndev = 1
      nmax = la_max+lb_max+ndev+1
      ALLOCATE (f(0:nmax))
      ! Loops over all pairs of primitive Gaussian-type functions
      DO ipgfa = 1, npgfa
         DO jpgfb = 1, npgfb

            ! Calculate some prefactors
            a = zeta(ipgfa)
            b = zetb(jpgfb)
            zet = a+b
            xhi = a*b/zet
            eta = xhi/(xhi+omega)
            oeta = omega*eta
            T = xhi*rab2
            expT = EXP(-omega/(omega+xhi)*T)
            prefac = pi**3/SQRT(zet**3)/SQRT((xhi+omega)**3)*expT

            DO ids = 1, nmax
               n = ids-1
               v(ipgfa, jpgfb, ids) = prefac*(-oeta)**n
            ENDDO

         END DO
      END DO
      DEALLOCATE (f)

      CALL timestop(handle)

   END SUBROUTINE s_gauss_ab

! **************************************************************************************************
!> \brief Contraction and normalization of the [s|O(r12)|s] integrals and their scalar derivatives;
!>        this routine is more efficient for uncontracted basis sets (clow), e.g. for ri basis sets
!> \param la set of l quantum numbers for a
!> \param npgfa number of primitive Gaussian on a
!> \param nshella number of shells for a
!> \param scona_shg SHG contraction/normalization matrix for a
!> \param lb set of l quantum numbers for b
!> \param npgfb number of primitive Gaussian on b
!> \param nshellb number of shells for b
!> \param sconb_shg SHG contraction/normalization matrix for b
!> \param swork matrix storing the uncontracted and unnormalized s-integrals and derivatives
!> \param swork_cont matrix storing the contracted and normalized s-integrals and derivatives
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE contract_sint_ab_clow(la, npgfa, nshella, scona_shg, lb, npgfb, nshellb, sconb_shg, &
                                    swork, swork_cont, calculate_forces)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: la
      INTEGER, INTENT(IN)                                :: npgfa, nshella
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: scona_shg
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lb
      INTEGER, INTENT(IN)                                :: npgfb, nshellb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sconb_shg
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: swork
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: swork_cont
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_sint_ab_clow', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ids, ids_start, ipgfa, ishella, &
                                                            jpgfb, jshellb, lai, lbj, ndev, nds

      CALL timeset(routineN, handle)

      ndev = 0
      IF (calculate_forces) ndev = 1

      swork_cont = 0.0_dp
      DO ishella = 1, nshella
         lai = la(ishella)
         DO jshellb = 1, nshellb
            lbj = lb(jshellb)
            nds = lai+lbj+1
            ids_start = nds-MIN(lai, lbj)
            DO ipgfa = 1, npgfa
               DO jpgfb = 1, npgfb
                  DO ids = ids_start, nds+ndev
                     swork_cont(ids, ishella, jshellb) = swork_cont(ids, ishella, jshellb) &
                                                         +scona_shg(ipgfa, ishella) &
                                                         *sconb_shg(jpgfb, jshellb) &
                                                         *swork(ipgfa, jpgfb, ids)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE contract_sint_ab_clow

! **************************************************************************************************
!> \brief Contraction and normalization of the [s|O(r12)|s] integrals and their scalar derivatives;
!>        this routine is more efficient for highly contracted basis sets (chigh)
!> \param npgfa number of primitive Gaussian on a
!> \param nshella number of shells for a
!> \param scona SHG contraction/normalization matrix for a
!> \param npgfb number of primitive Gaussian on b
!> \param nshellb number of shells for b
!> \param sconb SHG contraction/normalization matrix for b
!> \param nds maximal derivative of [s|O(r12)|s] with respect to rab2
!> \param swork matrix storing the uncontracted and unnormalized s-integrals and derivatives
!> \param swork_cont matrix storing the contracted and normalized s-integrals and derivatives
! **************************************************************************************************
   SUBROUTINE contract_sint_ab_chigh(npgfa, nshella, scona, &
                                     npgfb, nshellb, sconb, &
                                     nds, swork, swork_cont)

      INTEGER, INTENT(IN)                                :: npgfa, nshella
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: scona
      INTEGER, INTENT(IN)                                :: npgfb, nshellb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sconb
      INTEGER, INTENT(IN)                                :: nds
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: swork
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: swork_cont

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_sint_ab_chigh', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: work_cc, work_pc

      CALL timeset(routineN, handle)

      swork_cont = 0.0_dp
      ALLOCATE (work_pc(npgfb, nds, nshella))
      ALLOCATE (work_cc(nds, nshella, nshellb))

      CALL dgemm("T", "N", npgfb*nds, nshella, npgfa, 1._dp, swork, npgfa, scona, npgfa, &
                 0.0_dp, work_pc, npgfb*nds)
      CALL dgemm("T", "N", nds*nshella, nshellb, npgfb, 1._dp, work_pc, npgfb, sconb, npgfb, &
                 0.0_dp, swork_cont, nds*nshella)

      DEALLOCATE (work_pc, work_cc)
      CALL timestop(handle)

   END SUBROUTINE contract_sint_ab_chigh

! **************************************************************************************************
!> \brief Contraction, normalization and combinatiorial combination of the [0a|(r-Ra)^(2m)|0b]
!>         integrals and their scalar derivatives
!> \param npgfa number of primitive Gaussian on a
!> \param nshella number of shells for a
!> \param scon_ra2m contraction matrix on a containg the combinatorial factors
!> \param npgfb number of primitive Gaussian on b
!> \param nshellb number of shells for b
!> \param sconb SHG contraction/normalization matrix for b
!> \param swork matrix storing the uncontracted and unnormalized s-integrals and derivatives
!> \param swork_cont matrix storing the contracted and normalized s-integrals and derivatives
!> \param m exponent in operator (r-Ra)^(2m)
!> \param nds_max maximal derivative with respect to rab2
! **************************************************************************************************
   SUBROUTINE contract_s_ra2m_ab(npgfa, nshella, scon_ra2m, npgfb, nshellb, sconb, swork, swork_cont, &
                                 m, nds_max)

      INTEGER, INTENT(IN)                                :: npgfa, nshella
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: scon_ra2m
      INTEGER, INTENT(IN)                                :: npgfb, nshellb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sconb
      REAL(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(INOUT)                                   :: swork
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: swork_cont
      INTEGER, INTENT(IN)                                :: m, nds_max

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_s_ra2m_ab', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, my_m
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: work_pc

      CALL timeset(routineN, handle)
      my_m = m+1
      ALLOCATE (work_pc(npgfb, nds_max, nshella))
      work_pc = 0.0_dp
      swork_cont = 0.0_dp
      DO i = 1, my_m
         CALL dgemm("T", "N", npgfb*nds_max, nshella, npgfa, 1.0_dp, swork(1:npgfa, 1:npgfb, i, 1:nds_max), npgfa, &
                    scon_ra2m(1:npgfa, i, 1:nshella), npgfa, 1.0_dp, work_pc, npgfb*nds_max)
      ENDDO
      CALL dgemm("T", "N", nds_max*nshella, nshellb, npgfb, 1.0_dp, work_pc, npgfb, sconb, npgfb, 0.0_dp, &
                 swork_cont, nds_max*nshella)

      DEALLOCATE (work_pc)

      CALL timestop(handle)

   END SUBROUTINE contract_s_ra2m_ab

! **************************************************************************************************
!> \brief Contraction and normalization of the (abb) overlap
!> \param la set of l quantum numbers on a
!> \param npgfa number of primitive Gaussians functions on a; orbital basis
!> \param nshella number of shells for a; orbital basis
!> \param scona_shg SHG contraction/normalization matrix for a; orbital basis
!> \param lb set of l quantum numbers on b; orbital basis
!> \param npgfb number of primitive Gaussians functions on b; orbital basis
!> \param nshellb number of shells for b; orbital basis
!> \param lcb set of l quantum numbers on b; ri basis
!> \param npgfcb number of primitive Gaussians functions on b; ri basis
!> \param nshellcb number of shells for b; ri basis
!> \param orbb_index index for orbital basis at B for sconb_mix
!> \param rib_index index for ri basis at B for sconb_mix
!> \param sconb_mix  mixed contraction matrix for orbital and ri basis at B
!> \param nl_max related to the parameter m in (a|rb^(2m)|b)
!> \param nds_max derivative with respect to rab**2
!> \param swork matrix of storing the uncontracted and unnormalized s-integrals and derivatives
!> \param swork_cont matrix storing the contracted and normalized s-integrals and derivatives
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE contract_s_overlap_abb(la, npgfa, nshella, scona_shg, lb, npgfb, nshellb, &
                                     lcb, npgfcb, nshellcb, orbb_index, rib_index, sconb_mix, &
                                     nl_max, nds_max, swork, swork_cont, calculate_forces)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: la
      INTEGER, INTENT(IN)                                :: npgfa, nshella
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: scona_shg
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lb
      INTEGER, INTENT(IN)                                :: npgfb, nshellb
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lcb
      INTEGER, INTENT(IN)                                :: npgfcb, nshellcb
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: orbb_index, rib_index
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: sconb_mix
      INTEGER, INTENT(IN)                                :: nl_max, nds_max
      REAL(KIND=dp), DIMENSION(:, :, :, :, :), &
         INTENT(IN)                                      :: swork
      REAL(KIND=dp), DIMENSION(:, 0:, :, :, :), &
         INTENT(INOUT)                                   :: swork_cont
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'contract_s_overlap_abb', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: forb, fri, handle, ids, ids_start, iil, &
                                                            il, ishella, jpgfb, jshellb, kpgfb, &
                                                            kshellb, lai, lbb, lbj, lbk, ndev, &
                                                            nds, nl
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :)                        :: work_ppc

      CALL timeset(routineN, handle)

      ndev = 0
      IF (calculate_forces) ndev = 1

      ALLOCATE (work_ppc(npgfb, npgfcb, nl_max, nds_max, nshella))
      work_ppc = 0.0_dp
      CALL dgemm("T", "N", npgfb*npgfcb*nl_max*nds_max, nshella, npgfa, 1._dp, swork, npgfa, &
                 scona_shg, npgfa, 0.0_dp, work_ppc, npgfb*npgfcb*nl_max*nds_max)
      swork_cont = 0.0_dp
      DO kshellb = 1, nshellcb
         lbk = lcb(kshellb)
         DO jshellb = 1, nshellb
            lbj = lb(jshellb)
            nl = INT((lbj+lbk)/2)
            IF (lbj == 0 .OR. lbk == 0) nl = 0
            DO ishella = 1, nshella
               lai = la(ishella)
               DO il = 0, nl
                  lbb = lbj+lbk-2*il
                  nds = lai+lbb+1
                  ids_start = nds-MIN(lai, lbb)
                  DO jpgfb = 1, npgfb
                     forb = orbb_index(jpgfb, jshellb)
                     DO kpgfb = 1, npgfcb
                        fri = rib_index(kpgfb, kshellb)
                        DO ids = ids_start, nds+ndev
                           DO iil = 0, il
                              swork_cont(ids, il, ishella, jshellb, kshellb) = &
                                 swork_cont(ids, il, ishella, jshellb, kshellb) &
                                 +sconb_mix(iil+1, fri, forb, il+1)*work_ppc(jpgfb, kpgfb, il-iil+1, ids, ishella)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DEALLOCATE (work_ppc)
      CALL timestop(handle)

   END SUBROUTINE contract_s_overlap_abb

! **************************************************************************************************
!> \brief Contraction and normalization of the (aba) overlap
!> \param la set of l quantum numbers on a; orbital basis
!> \param npgfa number of primitive Gaussians functions on a; orbital basis
!> \param nshella number of shells for a; orbital basis
!> \param lb set of l quantum numbers on b; orbital basis
!> \param npgfb number of primitive Gaussians functions on b; orbital basis
!> \param nshellb number of shells for b; orbital basis
!> \param sconb_shg SHG contraction/normalization matrix for b; orbital basis
!> \param lca set of l quantum numbers on a; ri basis
!> \param npgfca number of primitive Gaussians functions on a; ri basis
!> \param nshellca number of shells for a; ri basis
!> \param orba_index index for orbital basis at A for scona_mix
!> \param ria_index index for ri basis at A for scona_mix
!> \param scona_mix mixed contraction matrix for orbital and ri basis at A
!> \param nl_max related to the parameter m in (a|ra^(2m)|b)
!> \param nds_max maximal derivative with respect to rab**2
!> \param swork matrix of storing the uncontracted and unnormalized s-integrals and derivatives
!> \param swork_cont matrix storing the contracted and normalized s-integrals and derivatives
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE contract_s_overlap_aba(la, npgfa, nshella, lb, npgfb, nshellb, sconb_shg, &
                                     lca, npgfca, nshellca, orba_index, ria_index, scona_mix, &
                                     nl_max, nds_max, swork, swork_cont, calculate_forces)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: la
      INTEGER, INTENT(IN)                                :: npgfa, nshella
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lb
      INTEGER, INTENT(IN)                                :: npgfb, nshellb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sconb_shg
      INTEGER, DIMENSION(:), INTENT(IN)                  :: lca
      INTEGER, INTENT(IN)                                :: npgfca, nshellca
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: orba_index, ria_index
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: scona_mix
      INTEGER, INTENT(IN)                                :: nl_max, nds_max
      REAL(KIND=dp), DIMENSION(:, :, :, :, :), &
         INTENT(IN)                                      :: swork
      REAL(KIND=dp), DIMENSION(:, 0:, :, :, :), &
         INTENT(INOUT)                                   :: swork_cont
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'contract_s_overlap_aba', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: forb, fri, handle, ids, ids_start, iil, &
                                                            il, ipgfa, ishella, jshellb, kpgfa, &
                                                            kshella, laa, lai, lak, lbj, ndev, &
                                                            nds, nl
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :)                        :: work_ppc

      CALL timeset(routineN, handle)

      ndev = 0
      IF (calculate_forces) ndev = 1

      ALLOCATE (work_ppc(npgfa, npgfca, nl_max, nds_max, nshellb))
      work_ppc = 0.0_dp
      CALL dgemm("T", "N", npgfa*npgfca*nl_max*nds_max, nshellb, npgfb, 1._dp, swork, npgfb, &
                 sconb_shg, npgfb, 0.0_dp, work_ppc, npgfa*npgfca*nl_max*nds_max)
      swork_cont = 0.0_dp
      DO kshella = 1, nshellca
         lak = lca(kshella)
         DO jshellb = 1, nshellb
            lbj = lb(jshellb)
            DO ishella = 1, nshella
               lai = la(ishella)
               nl = INT((lai+lak)/2)
               IF (lai == 0 .OR. lak == 0) nl = 0
               DO il = 0, nl
                  laa = lai+lak-2*il
                  nds = laa+lbj+1
                  ids_start = nds-MIN(laa, lbj)
                  DO kpgfa = 1, npgfca
                     fri = ria_index(kpgfa, kshella)
                     DO ipgfa = 1, npgfa
                        forb = orba_index(ipgfa, ishella)
                        DO ids = ids_start, nds+ndev
                           DO iil = 0, il
                              swork_cont(ids, il, ishella, jshellb, kshella) = &
                                 swork_cont(ids, il, ishella, jshellb, kshella) &
                                 +scona_mix(iil+1, fri, forb, il+1)*work_ppc(ipgfa, kpgfa, il-iil+1, ids, jshellb)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE (work_ppc)
      CALL timestop(handle)

   END SUBROUTINE contract_s_overlap_aba

END MODULE s_contract_shg

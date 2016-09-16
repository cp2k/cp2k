!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of integrals over Cartesian Gaussian-type functions for [a|(r-Ra)^(2m)|b]
!>        Ra is the postion of center a
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      none
!> \par Parameters
!>       - ax,ay,az    : Angular momentum index numbers of orbital a.
!>       - cx,cy,cz    : Angular momentum index numbers of orbital c.
!>       - coset       : Cartesian orbital set pointer.
!>       - dac         : Distance between the atomic centers a and c.
!>       - l{a,c}      : Angular momentum quantum number of shell a or c.
!>       - l{a,c}_max  : Maximum angular momentum quantum number of shell a or c.
!>       - l{a,c}_min  : Minimum angular momentum quantum number of shell a or c.
!>       - ncoset      : Number of orbitals in a Cartesian orbital set.
!>       - npgf{a,c}   : Degree of contraction of shell a or c.
!>       - rac         : Distance vector between the atomic centers a and c.
!>       - rac2        : Square of the distance between the atomic centers a and c.
!>       - zet{a,c}    : Exponents of the Gaussian-type functions a or c.
!>       - zetp        : Reciprocal of the sum of the exponents of orbital a and b.
!>       - zetw        : Reciprocal of the sum of the exponents of orbital a and c.
!> \author Dorothea Golze (08.2016)
! **************************************************************************************************
MODULE ai_operator_ra2m

   USE ai_os_rr,                        ONLY: os_rr_ovlp
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fac,&
                                              pi
   USE orbital_pointers,                ONLY: coset,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_operator_ra2m'

   PRIVATE

   ! *** Public subroutines ***

   PUBLIC :: operator_ra2m

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Calculation of the primitive two-center [a|(r-Ra)^(2m)|b] integrals over Cartesian
!>        Gaussian-type functions; operator is here (r-Ra)^(2m)
!> \param la_max ...
!> \param la_min ...
!> \param npgfa ...
!> \param zeta ...
!> \param lb_max ...
!> \param lb_min ...
!> \param npgfb ...
!> \param zetb ...
!> \param m exponent in (r-Ra)^(2m) operator
!> \param rab ...
!> \param sab ...
!> \param dsab ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE operator_ra2m(la_max, la_min, npgfa, zeta, &
                            lb_max, lb_min, npgfb, zetb, &
                            m, rab, sab, dsab, calculate_forces)
      INTEGER, INTENT(IN)                                :: la_max, la_min, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta
      INTEGER, INTENT(IN)                                :: lb_max, lb_min, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb
      INTEGER, INTENT(IN)                                :: m
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: sab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: dsab
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'operator_ra2m', routineP = moduleN//':'//routineN

      INTEGER                                            :: ax, ay, az, bx, by, bz, coa, cob, &
                                                            handle, i, ia, ib, ipgf, j, jpgf, k, &
                                                            la, lb, ldrr, lma, lmb, ma, mb, na, &
                                                            nb, ofa, ofb
      REAL(KIND=dp)                                      :: a, b, dumx, dumy, dumz, f0, prefac, &
                                                            rab2, tab, xhi, zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: rr
      REAL(KIND=dp), DIMENSION(3)                        :: rap, rbp

      CALL timeset(routineN, handle)

      sab = 0.0_dp
      IF (calculate_forces) dsab = 0.0_dp

      ! Distance of the centers a and b

      rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
      tab = SQRT(rab2)

      ! Maximum l for auxiliary integrals
      lma = la_max+2*m
      lmb = lb_max
      IF (calculate_forces) lma = lma+1
      ldrr = MAX(lma, lmb)+1

      ! Allocate space for auxiliary integrals
      ALLOCATE (rr(0:ldrr-1, 0:ldrr-1, 3))

      ! Number of integrals, check size of arrays
      ofa = ncoset(la_min-1)
      ofb = ncoset(lb_min-1)
      na = ncoset(la_max)-ofa
      nb = ncoset(lb_max)-ofb
      CPASSERT((SIZE(sab, 1) >= na*npgfa))
      CPASSERT((SIZE(sab, 2) >= nb*npgfb))

      ! Loops over all pairs of primitive Gaussian-type functions
      ma = 0
      DO ipgf = 1, npgfa
         mb = 0
         DO jpgf = 1, npgfb

            ! Calculate some prefactors
            a = zeta(ipgf)
            b = zetb(jpgf)
            zet = a+b
            xhi = a*b/zet
            rap = b*rab/zet
            rbp = -a*rab/zet

            ! [s|s] integral
            f0 = (pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

            ! Calculate the reccurence relation
            CALL os_rr_ovlp(rap, lma, rbp, lmb, zet, ldrr, rr)

            DO lb = lb_min, lb_max
            DO bx = 0, lb
            DO by = 0, lb-bx
               bz = lb-bx-by
               cob = coset(bx, by, bz)-ofb
               ib = mb+cob
               DO la = la_min, la_max
               DO ax = 0, la
               DO ay = 0, la-ax
                  az = la-ax-ay
                  coa = coset(ax, ay, az)-ofa
                  ia = ma+coa
                  DO i = 0, m
                  DO j = 0, m
                  DO k = 0, m
                     IF (i+j+k /= m) CYCLE
                     prefac = fac(m)/fac(i)/fac(j)/fac(k)
                     sab(ia, ib) = sab(ia, ib)+prefac*f0 &
                                   *rr(ax+2*i, bx, 1)*rr(ay+2*j, by, 2)*rr(az+2*k, bz, 3)
                     IF (calculate_forces) THEN
                        ! (da|b) = 2*a*(a+1|b) - N(a)*(a-1|b)
                        ! dx
                        dumx = 2.0_dp*a*rr(ax+2*i+1, bx, 1)
                        IF (ax+2*i > 0) dumx = dumx-REAL(ax+2*i, dp)*rr(ax+2*i-1, bx, 1)
                        dsab(ia, ib, 1) = dsab(ia, ib, 1)+prefac*f0*dumx*rr(ay+2*j, by, 2)*rr(az+2*k, bz, 3)
                        ! dy
                        dumy = 2.0_dp*a*rr(ay+2*j+1, by, 2)
                        IF (ay+2*j > 0) dumy = dumy-REAL(ay+2*j, dp)*rr(ay+2*j-1, by, 2)
                        dsab(ia, ib, 2) = dsab(ia, ib, 2)+prefac*f0*rr(ax+2*i, bx, 1)*dumy*rr(az+2*k, bz, 3)
                        ! dz
                        dumz = 2.0_dp*a*rr(az+2*k+1, bz, 3)
                        IF (az+2*k > 0) dumz = dumz-REAL(az+2*k, dp)*rr(az+2*k-1, bz, 3)
                        dsab(ia, ib, 3) = dsab(ia, ib, 3)+prefac*f0*rr(ax+2*i, bx, 1)*rr(ay+2*j, by, 2)*dumz
                     ENDIF
                  ENDDO
                  ENDDO
                  ENDDO
                  !
               ENDDO
               ENDDO
               ENDDO !la
            ENDDO
            ENDDO
            ENDDO !lb

            mb = mb+nb
         END DO
         ma = ma+na
      END DO

      DEALLOCATE (rr)

      CALL timestop(handle)

   END SUBROUTINE operator_ra2m

END MODULE ai_operator_ra2m

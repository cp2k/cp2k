!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Calculate the first derivative of an integral block.
!> \author  Matthias Krack
!> \date    05.10.2000
!> \version 1.0
!> \par Literature
!>          S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par Parameters
!>      - ax,ay,az  : Angular momentum index numbers of orbital a.
!>      - bx,by,bz  : Angular momentum index numbers of orbital b.
!>      - coset     : Cartesian orbital set pointer.
!>      - l{a,b}    : Angular momentum quantum number of shell a or b.
!>      - l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
!>      - l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
!>      - ncoset    : Number of orbitals in a Cartesian orbital set.
!>      - npgf{a,b} : Degree of contraction of shell a or b.
!>      - rab       : Distance vector between the atomic centers a and b.
!>      - rab2      : Square of the distance between the atomic centers a and b.
!>      - rac       : Distance vector between the atomic centers a and c.
!>      - rac2      : Square of the distance between the atomic centers a and c.
!>      - rbc       : Distance vector between the atomic centers b and c.
!>      - rbc2      : Square of the distance between the atomic centers b and c.
!>      - rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
!>      - zet{a,b}  : Exponents of the Gaussian-type functions a or b.
!>      - zetp      : Reciprocal of the sum of the exponents of orbital a and b.
! *****************************************************************************
MODULE ai_derivatives
  
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***
  PUBLIC :: adbdr, dabdr

CONTAINS

! *****************************************************************************
!> \brief   Calculate the first derivative of an integral block.
!>          This takes the derivative  with respect to
!>          the atomic position Ra, i.e. the center of the primitive on the left.
!>          To get the derivative of the left primitive with respect to r (orbital
!>           coordinate), take the opposite sign.
!>          To get the derivative with respect to the center of the primitive on
!>          the right Rb, take the opposite sign.
!>          To get the derivative of the right primitive with respect to r,
!>          do not change the sign.
!>          [da/dRi|O|b] = 2*zeta*[a+1i|O|b] - Ni(a)[a-1i|O|b]
!>
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param dab ...
!> \param ab ...
!> \param dabdx ...
!> \param dabdy ...
!> \param dabdz ...
!> \date    05.10.2000
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE dabdr(la_max,npgfa,zeta,rpgfa,la_min,lb_max,npgfb,rpgfb,lb_min,&
       dab,ab,dabdx,dabdy,dabdz)
    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: la_min, lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), INTENT(IN)                :: dab
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: ab
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: dabdx, dabdy, dabdz

    INTEGER :: ax, ay, az, bx, by, bz, coa, coamx, coamy, coamz, coapx, &
      coapy, coapz, cob, codb, i, ipgf, j, jpgf, la, lb, na, nb, nda, ndb
    REAL(KIND=dp)                            :: fa, fx, fy, fz

!   *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0
    nda = 0

    dabdx=0.0_dp
    dabdy=0.0_dp
    dabdz=0.0_dp

    DO ipgf=1,npgfa

       fa = 2.0_dp*zeta(ipgf)

       nb = 0
       ndb = 0

       DO jpgf=1,npgfb

          !       *** Screening ***

          IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
             DO j=nb+ncoset(lb_min-1)+1,nb+ncoset(lb_max)
                DO i=na+ncoset(la_min-1)+1,na+ncoset(la_max)
                   dabdx(i,j) = 0.0_dp
                   dabdy(i,j) = 0.0_dp
                   dabdz(i,j) = 0.0_dp
                END DO
             END DO
             nb = nb + ncoset(lb_max)
             ndb = ndb + ncoset(lb_max+1)
             CYCLE
          END IF

          !       *** [da/dRi|O|b] = 2*zeta*[a+1i|O|b] - Ni(a)[a-1i|O|b] ***

          DO la=0,la_max !MAX(0,la_min),la_max

             IF (la == 0) THEN

                coa = na + 1
                coapx = nda + 2
                coapy = nda + 3
                coapz = nda + 4

                DO lb= 0,lb_max !lb_min,lb_max
                   DO bx=0,lb
                      DO by=0,lb-bx
                         bz = lb - bx - by
                         cob = nb + coset(bx,by,bz)
                         codb = ndb + coset(bx,by,bz)
                         dabdx(coa,cob) = fa*ab(coapx,codb)
                         dabdy(coa,cob) = fa*ab(coapy,codb)
                         dabdz(coa,cob) = fa*ab(coapz,codb)
                      END DO
                   END DO
                END DO

             ELSE

                DO ax=0,la
                   DO ay=0,la-ax
                      az = la - ax - ay

                      coa = na + coset(ax,ay,az)
                      coamx = nda + coset(MAX(0,ax-1),ay,az)
                      coamy = nda + coset(ax,MAX(0,ay-1),az)
                      coamz = nda + coset(ax,ay,MAX(0,az-1))
                      coapx = nda + coset(ax+1,ay,az)
                      coapy = nda + coset(ax,ay+1,az)
                      coapz = nda + coset(ax,ay,az+1)

                      fx = REAL(ax,dp)
                      fy = REAL(ay,dp)
                      fz = REAL(az,dp)

                      DO lb=0,lb_max !lb_min,lb_max
                         DO bx=0,lb
                            DO by=0,lb-bx
                               bz = lb - bx - by
                               cob = nb + coset(bx,by,bz)
                               codb = ndb + coset(bx,by,bz)
                               dabdx(coa,cob) = fa*ab(coapx,codb) - fx*ab(coamx,codb)
                               dabdy(coa,cob) = fa*ab(coapy,codb) - fy*ab(coamy,codb)
                               dabdz(coa,cob) = fa*ab(coapz,codb) - fz*ab(coamz,codb)
                            END DO
                         END DO
                      END DO

                   END DO
                END DO

             END IF

          END DO

          nb = nb + ncoset(lb_max)
          ndb = ndb + ncoset(lb_max+1)

       END DO

       na = na + ncoset(la_max)
       nda = nda + ncoset(la_max+1)

    END DO

  END SUBROUTINE dabdr

! *****************************************************************************
!> \brief   Calculate the first derivative of an integral block.
!>          This takes the derivative  with respect to
!>          the atomic position Rb, i.e. the center of the primitive on the right.
!>          To get the derivative of the left primitive with respect to r
!>          (orbital coordinate), take the opposite sign.
!>          [a|O|db/dRi] = 2*zetb*[a|O|b+1i] - Ni(b)[a|O|b-1i]
!> \param la_max ...
!> \param npgfa ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param dab ...
!> \param ab ...
!> \param adbdx ...
!> \param adbdy ...
!> \param adbdz ...
!> \date    05.10.2000
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE adbdr(la_max,npgfa,rpgfa,la_min,lb_max,npgfb,zetb,rpgfb,lb_min,&
       dab,ab,adbdx,adbdy,adbdz)
    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa
    INTEGER, INTENT(IN)                      :: la_min, lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb, rpgfb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), INTENT(IN)                :: dab
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: ab
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: adbdx, adbdy, adbdz

    INTEGER :: ax, ay, az, bx, by, bz, coa, cob, cobmx, cobmy, cobmz, cobpx, &
      cobpy, cobpz, coda, i, ipgf, j, jpgf, la, lb, na, nb, nda, ndb
    REAL(KIND=dp)                            :: fb, fx, fy, fz

    na = 0
    nda = 0

    adbdx=0.0_dp
    adbdy=0.0_dp
    adbdz=0.0_dp
    !   *** Loop over all pairs of primitive Gaussian-type functions ***
    DO ipgf=1,npgfa

       nb = 0
       ndb = 0

       DO jpgf=1,npgfb

          fb = 2.0_dp*zetb(jpgf)

          !       *** Screening ***

          IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
             DO j=nb+ncoset(lb_min-1)+1,nb+ncoset(lb_max)
                DO i=na+ncoset(la_min-1)+1,na+ncoset(la_max)
                   adbdx(i,j) = 0.0_dp
                   adbdy(i,j) = 0.0_dp
                   adbdz(i,j) = 0.0_dp
                END DO
             END DO
             nb = nb + ncoset(lb_max)
             ndb = ndb + ncoset(lb_max+1)
             CYCLE
          END IF

          !       *** [a|O|db/dRi] = 2*zetb*[a|O|b+1i] - Ni(b)[a|O|b-1i] ***

          DO lb = 0,lb_max

             IF(lb==0) THEN

                cob = nb + 1
                cobpx = ndb + 2
                cobpy = ndb + 3
                cobpz = ndb + 4

                DO la=0,la_max !la_min,la_max
                   DO ax=0,la
                      DO ay=0,la-ax
                         az = la -ax -ay
                         coa = na + coset(ax,ay,az)
                         coda = nda + coset(ax,ay,az)
                         adbdx(coa,cob) = fb*ab(coda,cobpx)
                         adbdy(coa,cob) = fb*ab(coda,cobpy)
                         adbdz(coa,cob) = fb*ab(coda,cobpz)
                      END DO
                   END DO
                END DO
             ELSE

                DO bx=0,lb
                   DO by=0,lb-bx
                      bz = lb - bx - by

                      cob = nb + coset(bx,by,bz)
                      cobmx = ndb + coset(MAX(0,bx-1),by,bz)
                      cobmy = ndb + coset(bx,MAX(0,by-1),bz)
                      cobmz = ndb + coset(bx,by,MAX(0,bz-1))
                      cobpx = ndb + coset(bx+1,by,bz)
                      cobpy = ndb + coset(bx,by+1,bz)
                      cobpz = ndb + coset(bx,by,bz+1)

                      fx = REAL(bx,dp)
                      fy = REAL(by,dp)
                      fz = REAL(bz,dp)

                      DO la=0,la_max !la_min,la_max
                         DO ax=0,la
                            DO ay=0,la-ax
                               az = la -ax -ay
                               coa = na + coset(ax,ay,az)
                               coda = nda + coset(ax,ay,az)
                               adbdx(coa,cob) = fb*ab(coda,cobpx) - fx*ab(coda,cobmx)
                               adbdy(coa,cob) = fb*ab(coda,cobpy) - fy*ab(coda,cobmy)
                               adbdz(coa,cob) = fb*ab(coda,cobpz) - fz*ab(coda,cobmz)
                            END DO
                         END DO
                      END DO

                   END DO
                END DO

             END IF

          END DO

          nb = nb + ncoset(lb_max)
          ndb = ndb + ncoset(lb_max+1)

       END DO

       na = na + ncoset(la_max)
       nda = nda + ncoset(la_max+1)

    END DO

  END SUBROUTINE adbdr

END MODULE ai_derivatives


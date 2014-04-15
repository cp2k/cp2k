!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the kinetic energy integrals over Cartesian
!>      Gaussian-type functions.
!>
!>      [a|T|b] = [a|-nabla**2/2|b]
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      - Derivatives added (10.05.2002,MK)
!> \par Parameters
!>      - ax,ay,az  : Angular momentum index numbers of orbital a.
!>      - bx,by,bz  : Angular momentum index numbers of orbital b.
!>      - coset     : Cartesian orbital set pointer.
!>      - dab       : Distance between the atomic centers a and b.
!>      - l{a,b}    : Angular momentum quantum number of shell a or b.
!>      - l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
!>      - l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
!>      - ncoset    : Number of Cartesian orbitals up to l.
!>      - rab       : Distance vector between the atomic centers a and b.
!>      - rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
!>      - sab       : Shell set of overlap integrals.
!>      - tab       : Shell set of kinetic energy integrals.
!>      - zet{a,b}  : Exponents of the Gaussian functions.
!> \author Matthias Krack (31.07.2000)
! *****************************************************************************
MODULE ai_kinetic
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_kinetic'

! *** Public subroutines ***

  PUBLIC :: kinetic

CONTAINS

! *****************************************************************************
!> \brief    Calculation of the two-center kinetic energy integrals [a|T|b]
!>           over Cartesian Gaussian-type functions.
!> \author  Matthias Krack
!> \date    22.09.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE kinetic(la_max_set,la_min_set,npgfa,rpgfa,zeta,&
                     lb_max_set,lb_min_set,npgfb,rpgfb,zetb,&
                     rab,dab,sdab,tab,da_max,return_derivatives,&
                     t,ldt,pab,force_a)
    INTEGER, INTENT(IN)                      :: la_max_set, la_min_set, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa, zeta
    INTEGER, INTENT(IN)                      :: lb_max_set, lb_min_set, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb, zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: dab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: sdab
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: tab
    INTEGER, INTENT(IN)                      :: da_max
    LOGICAL, INTENT(IN)                      :: return_derivatives
    INTEGER, INTENT(IN)                      :: ldt
    REAL(KIND=dp), DIMENSION(ldt, ldt, *), &
      INTENT(INOUT)                          :: t
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: pab
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(OUT), OPTIONAL                  :: force_a

    CHARACTER(len=*), PARAMETER :: routineN = 'kinetic', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ay, az, bx, by, bz, cda, cdax, cday, cdaz, coa, coam2x, &
      coam2y, coam2z, coamx, coamy, coamz, coapx, coapy, coapz, cob, cobm2x, &
      cobm2y, cobm2z, cobmx, cobmy, cobmz, da, dax, day, daz, i, ipgf, j, jk, &
      jpgf, jstart, k, la, la_max, la_min, la_start, lb, lb_max, lb_min, &
      lb_start, na, nb, nb1, nda, nda1
    LOGICAL                                  :: calculate_force_a
    REAL(KIND=dp)                            :: f0, f1, f2, f3, f4, fax, fay, &
                                                faz, ftz, xi, zetp
    REAL(KIND=dp), DIMENSION(3)              :: rap, rbp

    IF (PRESENT(pab).AND.PRESENT(force_a)) THEN
      calculate_force_a = .TRUE.
      force_a(:) = 0.0_dp
    ELSE
      calculate_force_a = .FALSE.
    END IF

    la_max = la_max_set + da_max
    la_min = MAX(0,la_min_set-da_max)

    lb_max = lb_max_set
    lb_min = lb_min_set

!   *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0
    nda = 0

    DO ipgf=1,npgfa

      nda1 = nda + 1
      nb = 0

      DO jpgf=1,npgfb

        nb1 = nb + 1

!       *** Screening ***

        IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
          DO j=nb1,nb+ncoset(lb_max_set)
            DO i=na+1,na+ncoset(la_max_set)
              tab(i,j) = 0.0_dp
            END DO
          END DO
          IF (return_derivatives) THEN
            DO k=2,ncoset(da_max)
              jstart = (k - 1)*SIZE(tab,1)
              DO j=jstart+nb1,jstart+nb+ncoset(lb_max_set)
                DO i=na+1,na+ncoset(la_max_set)
                  tab(i,j) = 0.0_dp
                END DO
              END DO
            END DO
          END IF
          nb = nb + ncoset(lb_max_set)
          CYCLE
        END IF

!       *** Calculate some prefactors ***

        zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))

        f0 = zeta(ipgf)*zetp
        f1 = zetb(jpgf)*zetp
        f2 = 0.5_dp*zetp
        xi = zeta(ipgf)*f1

!       *** Calculate the basic two-center kinetic energy integral [s|T|s] ***

        t(1,1,1) = xi*(3.0_dp - 2.0_dp*xi*dab*dab)*sdab(nda1,nb1,1)

!       *** Recurrence steps: [s|T|s] -> [a|T|b] ***

        IF (la_max > 0) THEN

!         *** Vertical recurrence steps: [s|T|s] -> [a|T|s] ***

          rap(:) = f1*rab(:)

!         *** [p|T|s] = (Pi - Ai)*[s|T|s] + f1*[ds/dAi|s]  (i = x,y,z) ***

          t(2,1,1) = rap(1)*t(1,1,1) + f1*sdab(nda1,nb1,2) ! [px|T|s]
          t(3,1,1) = rap(2)*t(1,1,1) + f1*sdab(nda1,nb1,3) ! [py|T|s]
          t(4,1,1) = rap(3)*t(1,1,1) + f1*sdab(nda1,nb1,4) ! [pz|T|s]

          IF (la_max > 1) THEN

!           *** [d|T|s] ***

            f3 = f2*t(1,1,1)

            t( 5,1,1) = rap(1)*t(2,1,1) + f1*sdab(nda+2,nb1,2) + f3 ! [dx2|T|s]
            t( 6,1,1) = rap(1)*t(3,1,1) + f1*sdab(nda+3,nb1,2)      ! [dxy|T|s]
            t( 7,1,1) = rap(1)*t(4,1,1) + f1*sdab(nda+4,nb1,2)      ! [dxz|T|s]
            t( 8,1,1) = rap(2)*t(3,1,1) + f1*sdab(nda+3,nb1,3) + f3 ! [dy2|T|s]
            t( 9,1,1) = rap(2)*t(4,1,1) + f1*sdab(nda+4,nb1,3)      ! [dyz|T|s]
            t(10,1,1) = rap(3)*t(4,1,1) + f1*sdab(nda+4,nb1,4) + f3 ! [dz2|T|s]

            IF (la_max > 2) THEN

!             *** [f|T|s] ***

              f3 = 2.0_dp*f2

              t(11,1,1) = rap(1)*t( 5,1,1) + f1*sdab(nda+ 5,nb1,2) + f3*t(2,1,1) ! [fx3 |T|s]
              t(12,1,1) = rap(1)*t( 6,1,1) + f1*sdab(nda+ 6,nb1,2) + f2*t(3,1,1) ! [fx2y|T|s]
              t(13,1,1) = rap(1)*t( 7,1,1) + f1*sdab(nda+ 7,nb1,2) + f2*t(4,1,1) ! [fx2z|T|s]
              t(14,1,1) = rap(1)*t( 8,1,1) + f1*sdab(nda+ 8,nb1,2)               ! [fxy2|T|s]
              t(15,1,1) = rap(1)*t( 9,1,1) + f1*sdab(nda+ 9,nb1,2)               ! [fxyz|T|s]
              t(16,1,1) = rap(1)*t(10,1,1) + f1*sdab(nda+10,nb1,2)               ! [fxz2|T|s]
              t(17,1,1) = rap(2)*t( 8,1,1) + f1*sdab(nda+ 8,nb1,3) + f3*t(3,1,1) ! [fy3 |T|s]
              t(18,1,1) = rap(2)*t( 9,1,1) + f1*sdab(nda+ 9,nb1,3) + f2*t(4,1,1) ! [fy2z|T|s]
              t(19,1,1) = rap(2)*t(10,1,1) + f1*sdab(nda+10,nb1,3)               ! [fyz2|T|s]
              t(20,1,1) = rap(3)*t(10,1,1) + f1*sdab(nda+10,nb1,4) + f3*t(4,1,1) ! [fz3 |T|s]

              IF (la_max > 3) THEN

!               *** [g|T|s] ***

                f4 = 3.0_dp*f2

                t(21,1,1) = rap(1)*t(11,1,1) + f1*sdab(nda+11,nb1,2) + f4*t( 5,1,1) ! [gx4  |T|s]
                t(22,1,1) = rap(1)*t(12,1,1) + f1*sdab(nda+12,nb1,2) + f3*t( 6,1,1) ! [gx3y |T|s]
                t(23,1,1) = rap(1)*t(13,1,1) + f1*sdab(nda+13,nb1,2) + f3*t( 7,1,1) ! [gx3z |T|s]
                t(24,1,1) = rap(1)*t(14,1,1) + f1*sdab(nda+14,nb1,2) + f2*t( 8,1,1) ! [gx2y2|T|s]
                t(25,1,1) = rap(1)*t(15,1,1) + f1*sdab(nda+15,nb1,2) + f2*t( 9,1,1) ! [gx2yz|T|s]
                t(26,1,1) = rap(1)*t(16,1,1) + f1*sdab(nda+16,nb1,2) + f2*t(10,1,1) ! [gx2z2|T|s]
                t(27,1,1) = rap(1)*t(17,1,1) + f1*sdab(nda+17,nb1,2)                ! [gxy3 |T|s]
                t(28,1,1) = rap(1)*t(18,1,1) + f1*sdab(nda+18,nb1,2)                ! [gxy2z|T|s]
                t(29,1,1) = rap(1)*t(19,1,1) + f1*sdab(nda+19,nb1,2)                ! [gxyz2|T|s]
                t(30,1,1) = rap(1)*t(20,1,1) + f1*sdab(nda+20,nb1,2)                ! [gxz3 |T|s]
                t(31,1,1) = rap(2)*t(17,1,1) + f1*sdab(nda+17,nb1,3) + f4*t( 8,1,1) ! [gy4  |T|s]
                t(32,1,1) = rap(2)*t(18,1,1) + f1*sdab(nda+18,nb1,3) + f3*t( 9,1,1) ! [gy3z |T|s]
                t(33,1,1) = rap(2)*t(19,1,1) + f1*sdab(nda+19,nb1,3) + f2*t(10,1,1) ! [gy2z2|T|s]
                t(34,1,1) = rap(2)*t(20,1,1) + f1*sdab(nda+20,nb1,3)                ! [gyz3 |T|s]
                t(35,1,1) = rap(3)*t(20,1,1) + f1*sdab(nda+20,nb1,4) + f4*t(10,1,1) ! [gz4  |T|s]

!               *** [a|T|s] = (Pi - Ai)*[a-1i|T|s] + ***
!               ***           f1*[d(a-1i)/dAi|s] +   ***
!               ***           f2*Ni(a-1i)*[a-2i|T|s] ***

                DO la=5,la_max

!                 *** Increase the angular momentum component z of function a ***

                  coa = coset(0,0,la)
                  coamz = coset(0,0,la-1)
                  coam2z = coset(0,0,la-2)
                  t(coa,1,1) = rap(3)*t(coamz,1,1) +&
                               f1*sdab(nda+coamz,nb1,4) +&
                               f2*REAL(la-1,dp)*t(coam2z,1,1)

!                 *** Increase the angular momentum component y of function a ***

                  az = la - 1
                  coa = coset(0,1,az)
                  coamy = coset(0,0,az)
                  t(coa,1,1) = rap(2)*t(coamy,1,1) +&
                               f1*sdab(nda+coamy,nb1,3)
                  DO ay=2,la
                    az = la - ay
                    coa = coset(0,ay,az)
                    coamy = coset(0,ay-1,az)
                    coam2y = coset(0,ay-2,az)
                    t(coa,1,1) = rap(2)*t(coamy,1,1) +&
                                 f1*sdab(nda+coamy,nb1,3) +&
                                 f2*REAL(ay-1,dp)*t(coam2y,1,1)
                  END DO

!                 *** Increase the angular momentum component x of function a ***

                  DO ay=0,la-1
                    az = la - 1 - ay
                    coa = coset(1,ay,az)
                    coamx = coset(0,ay,az)
                    t(coa,1,1) = rap(1)*t(coamx,1,1) +&
                                 f1*sdab(nda+coamx,nb1,2)
                  END DO
                  DO ax=2,la
                    f3 = f2*REAL(ax-1,dp)
                    DO ay=0,la-ax
                      az = la - ax - ay
                      coa = coset(ax,ay,az)
                      coamx = coset(ax-1,ay,az)
                      coam2x = coset(ax-2,ay,az)
                      t(coa,1,1) = rap(1)*t(coamx,1,1) +&
                                   f1*sdab(nda+coamx,nb1,2) +&
                                   f3*t(coam2x,1,1)
                    END DO
                  END DO

                END DO

              END IF

            END IF

          END IF

!         *** Recurrence steps: [a|T|s] -> [a|T|b] ***

          IF (lb_max > 0) THEN

            DO j=2,ncoset(lb_max)
              DO i=1,ncoset(la_min)
                t(i,j,1) = 0.0_dp
              END DO
            END DO

!           *** Horizontal recurrence steps ***

            rbp(:) = rap(:) - rab(:)

!           *** [a|T|p] = [a+1i|T|s] - (Bi - Ai)*[a|T|s] - [da/dAi|s] ***

            IF (lb_max == 1) THEN
              la_start = la_min
            ELSE
              la_start = MAX(0,la_min-1)
            END IF

            DO la=la_start,la_max-1
              DO ax=0,la
                DO ay=0,la-ax
                  az = la - ax - ay
                  coa = coset(ax,ay,az)
                  coapx = coset(ax+1,ay,az)
                  coapy = coset(ax,ay+1,az)
                  coapz = coset(ax,ay,az+1)
                  t(coa,2,1) = t(coapx,1,1) - rab(1)*t(coa,1,1) -&
                               sdab(nda+coa,nb1,2)
                  t(coa,3,1) = t(coapy,1,1) - rab(2)*t(coa,1,1) -&
                               sdab(nda+coa,nb1,3)
                  t(coa,4,1) = t(coapz,1,1) - rab(3)*t(coa,1,1) -&
                               sdab(nda+coa,nb1,4)
                END DO
              END DO
            END DO

!           *** Vertical recurrence step ***

!           *** [a|T|p] = (Pi - Bi)*[a|T|s] + f2*Ni(a)*[a-1i|T|s] - ***
!           ***           f0*[da/dAi|s]                             ***

            DO ax=0,la_max
              fax = f2*REAL(ax,dp)
              DO ay=0,la_max-ax
                fay = f2*REAL(ay,dp)
                az = la_max - ax - ay
                faz = f2*REAL(az,dp)
                coa = coset(ax,ay,az)
                coamx = coset(ax-1,ay,az)
                coamy = coset(ax,ay-1,az)
                coamz = coset(ax,ay,az-1)
                t(coa,2,1) = rbp(1)*t(coa,1,1) + fax*t(coamx,1,1) -&
                             f0*sdab(nda+coa,nb1,2)
                t(coa,3,1) = rbp(2)*t(coa,1,1) + fay*t(coamy,1,1) -&
                             f0*sdab(nda+coa,nb1,3)
                t(coa,4,1) = rbp(3)*t(coa,1,1) + faz*t(coamz,1,1) -&
                             f0*sdab(nda+coa,nb1,4)
              END DO
            END DO

!           *** Recurrence steps: [a|T|p] -> [a|T|b] ***

            DO lb=2,lb_max

!             *** Horizontal recurrence steps ***

!             *** [a|T|b] = [a+1i|T|b-1i] - (Bi - Ai)*[a|T|b-1i] - ***
!             ***           [da/dAi|b-1i]                          ***

              IF (lb == lb_max) THEN
                la_start = la_min
              ELSE
                la_start = MAX(0,la_min-1)
              END IF

              DO la=la_start,la_max-1
                DO ax=0,la
                  DO ay=0,la-ax
                    az = la - ax - ay
                    coa = coset(ax,ay,az)
                    coapx = coset(ax+1,ay,az)
                    coapy = coset(ax,ay+1,az)
                    coapz = coset(ax,ay,az+1)

!                   *** Shift of angular momentum component z from a to b ***

                    cob = coset(0,0,lb)
                    cobmz = coset(0,0,lb-1)
                    t(coa,cob,1) = t(coapz,cobmz,1) -&
                                   rab(3)*t(coa,cobmz,1) -&
                                   sdab(nda+coa,nb+cobmz,4)

!                   *** Shift of angular momentum component y from a to b ***

                    bz = lb - 1
                    cob = coset(0,1,bz)
                    cobmy = coset(0,0,bz)
                    t(coa,cob,1) = t(coapy,cobmy,1) -&
                                   rab(2)*t(coa,cobmy,1) -&
                                   sdab(nda+coa,nb+cobmy,3)
                    DO by=2,lb
                      bz = lb - by
                      cob = coset(0,by,bz)
                      cobmy = coset(0,by-1,bz)
                      t(coa,cob,1) = t(coapy,cobmy,1) -&
                                     rab(2)*t(coa,cobmy,1) -&
                                     sdab(nda+coa,nb+cobmy,3)
                    END DO

!                   *** Shift of angular momentum component x from a to b ***

                    DO by=0,lb-1
                      bz = lb - 1 - by
                      cob = coset(1,by,bz)
                      cobmx = coset(0,by,bz)
                      t(coa,cob,1) = t(coapx,cobmx,1) -&
                                     rab(1)*t(coa,cobmx,1) -&
                                     sdab(nda+coa,nb+cobmx,2)
                    END DO
                    DO bx=2,lb
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        cob = coset(bx,by,bz)
                        cobmx = coset(bx-1,by,bz)
                        t(coa,cob,1) = t(coapx,cobmx,1) -&
                                       rab(1)*t(coa,cobmx,1) -&
                                       sdab(nda+coa,nb+cobmx,2)
                      END DO
                    END DO

                  END DO
                END DO
              END DO

!             *** Vertical recurrence step ***

!             *** [a|T|b] = (Pi - Bi)*[a|T|b-1i] +   ***
!             ***           f2*Ni(a)*[a-1i|T|b-1i] + ***
!             ***           f2*Ni(b-1i)*[a|T|b-2i] - ***
!             ***           f0*[da/dAi|b-1i]         ***

              DO ax=0,la_max
                fax = f2*REAL(ax,dp)
                DO ay=0,la_max-ax
                  fay = f2*REAL(ay,dp)
                  az = la_max - ax - ay
                  faz = f2*REAL(az,dp)
                  coa = coset(ax,ay,az)
                  coamx = coset(ax-1,ay,az)
                  coamy = coset(ax,ay-1,az)
                  coamz = coset(ax,ay,az-1)

!                 *** Increase the angular momentum component z of b ***

                  f3 = f2*REAL(lb-1,dp)
                  cob = coset(0,0,lb)
                  cobmz = coset(0,0,lb-1)
                  cobm2z = coset(0,0,lb-2)
                  t(coa,cob,1) = rbp(3)*t(coa,cobmz,1) +&
                                 faz*t(coamz,cobmz,1) +&
                                 f3*t(coa,cobm2z,1) -&
                                 f0*sdab(nda+coa,nb+cobmz,4)

!                 *** Increase the angular momentum component y of b ***

                  bz = lb - 1
                  cob = coset(0,1,bz)
                  cobmy = coset(0,0,bz)
                  t(coa,cob,1) = rbp(2)*t(coa,cobmy,1) +&
                                 fay*t(coamy,cobmy,1) -&
                                 f0*sdab(nda+coa,nb+cobmy,3)
                  DO by=2,lb
                    bz = lb - by
                    f3 = f2*REAL(by-1,dp)
                    cob = coset(0,by,bz)
                    cobmy = coset(0,by-1,bz)
                    cobm2y = coset(0,by-2,bz)
                    t(coa,cob,1) = rbp(2)*t(coa,cobmy,1) +&
                                   fay*t(coamy,cobmy,1) +&
                                   f3*t(coa,cobm2y,1) -&
                                   f0*sdab(nda+coa,nb+cobmy,3)
                  END DO

!                 *** Increase the angular momentum component x of b ***

                  DO by=0,lb-1
                    bz = lb - 1 - by
                    cob = coset(1,by,bz)
                    cobmx = coset(0,by,bz)
                    t(coa,cob,1) = rbp(1)*t(coa,cobmx,1) +&
                                   fax*t(coamx,cobmx,1) -&
                                   f0*sdab(nda+coa,nb+cobmx,2)
                  END DO
                  DO bx=2,lb
                    f3 = f2*REAL(bx-1,dp)
                    DO by=0,lb-bx
                      bz = lb - bx - by
                      cob = coset(bx,by,bz)
                      cobmx = coset(bx-1,by,bz)
                      cobm2x = coset(bx-2,by,bz)
                      t(coa,cob,1) = rbp(1)*t(coa,cobmx,1) +&
                                     fax*t(coamx,cobmx,1) +&
                                     f3*t(coa,cobm2x,1) -&
                                     f0*sdab(nda+coa,nb+cobmx,2)
                    END DO
                  END DO

                END DO
              END DO

            END DO

          END IF

        ELSE

          IF (lb_max > 0) THEN

!           *** Vertical recurrence steps: [s|T|s] -> [s|T|b] ***

            rbp(:) = -f0*rab(:)

!           *** [s|T|p] = (Pi - Bi)*[s|T|s] - f0*[ds/dAi|s] ***

            t(1,2,1) = rbp(1)*t(1,1,1) - f0*sdab(nda1,nb1,2) ! [s|T|px]
            t(1,3,1) = rbp(2)*t(1,1,1) - f0*sdab(nda1,nb1,3) ! [s|T|py]
            t(1,4,1) = rbp(3)*t(1,1,1) - f0*sdab(nda1,nb1,4) ! [s|T|pz]

            IF (lb_max > 1) THEN

!             *** [d|T|s] ***

              f3 = f2*t(1,1,1)

              t(1, 5,1) = rbp(1)*t(1,2,1) - f0*sdab(nda1,nb+2,2) + f3 ! [s|T|dx2]
              t(1, 6,1) = rbp(1)*t(1,3,1) - f0*sdab(nda1,nb+3,2)      ! [s|T|dxy]
              t(1, 7,1) = rbp(1)*t(1,4,1) - f0*sdab(nda1,nb+4,2)      ! [s|T|dxz]
              t(1, 8,1) = rbp(2)*t(1,3,1) - f0*sdab(nda1,nb+3,3) + f3 ! [s|T|dy2]
              t(1, 9,1) = rbp(2)*t(1,4,1) - f0*sdab(nda1,nb+4,3)      ! [s|T|dyz]
              t(1,10,1) = rbp(3)*t(1,4,1) - f0*sdab(nda1,nb+4,4) + f3 ! [s|T|dz2]

!             *** [s|T|b] = (Pi - Bi)*[s|T|b-1i] - ***
!             ***           f0*[ds/dAi|b-1i] +     ***
!             ***           f2*Ni(b-1i)*[s|T|b-2i] ***

              DO lb=2,lb_max

!               *** Increase the angular momentum component z of function b ***

                cob = coset(0,0,lb)
                cobmz = coset(0,0,lb-1)
                cobm2z = coset(0,0,lb-2)
                t(1,cob,1) = rbp(3)*t(1,cobmz,1) -&
                             f0*sdab(nda1,nb+cobmz,4) +&
                             f2*REAL(lb-1,dp)*t(1,cobm2z,1)

!               *** Increase the angular momentum component y of function b ***

                bz = lb - 1
                cob = coset(0,1,bz)
                cobmy = coset(0,0,bz)
                t(1,cob,1) = rbp(2)*t(1,cobmy,1) -&
                             f0*sdab(nda1,nb+cobmy,3)
                DO by=2,lb
                  bz = lb - by
                  cob = coset(0,by,bz)
                  cobmy = coset(0,by-1,bz)
                  cobm2y = coset(0,by-2,bz)
                  t(1,cob,1) = rbp(2)*t(1,cobmy,1) -&
                               f0*sdab(nda1,nb+cobmy,3) +&
                               f2*REAL(by-1,dp)*t(1,cobm2y,1)
                END DO

!               *** Increase the angular momentum component x of function b ***

                DO by=0,lb-1
                  bz = lb - 1 - by
                  cob = coset(1,by,bz)
                  cobmx = coset(0,by,bz)
                  t(1,cob,1) = rbp(1)*t(1,cobmx,1) -&
                               f0*sdab(nda1,nb+cobmx,2)
                END DO
                DO bx=2,lb
                  f3 = f2*REAL(bx-1,dp)
                  DO by=0,lb-bx
                    bz = lb - bx - by
                    cob = coset(bx,by,bz)
                    cobmx = coset(bx-1,by,bz)
                    cobm2x = coset(bx-2,by,bz)
                    t(1,cob,1) = rbp(1)*t(1,cobmx,1) -&
                                 f0*sdab(nda1,nb+cobmx,2) +&
                                 f3*t(1,cobm2x,1)
                  END DO
                END DO

              END DO

            END IF

          END IF

        END IF

!       *** Store the primitive kinetic energy integrals ***

        DO j=1,ncoset(lb_max_set)
          DO i=1,ncoset(la_max_set)
            tab(na+i,nb+j) = t(i,j,1)
          END DO
        END DO

!       *** Calculate the requested derivatives with respect  ***
!       *** to the nuclear coordinates of the atomic center a ***

        IF (return_derivatives) THEN
          la_start = 0
          lb_start = 0
        ELSE
          la_start = la_min_set
          lb_start = lb_min_set
        END IF

        DO da=0,da_max-1
          ftz = 2.0_dp*zeta(ipgf)
          DO dax=0,da
            DO day=0,da-dax
              daz = da - dax - day
              cda = coset(dax,day,daz)
              cdax = coset(dax+1,day,daz)
              cday = coset(dax,day+1,daz)
              cdaz = coset(dax,day,daz+1)

!             *** [da/dAi|T|b] = 2*zeta*[a+1i|T|b] - Ni(a)[a-1i|T|b] ***

              DO la=la_start,la_max-da-1
                DO ax=0,la
                  fax = REAL(ax,dp)
                  DO ay=0,la-ax
                    fay = REAL(ay,dp)
                    az = la - ax - ay
                    faz = REAL(az,dp)
                    coa = coset(ax,ay,az)
                    coamx = coset(ax-1,ay,az)
                    coamy = coset(ax,ay-1,az)
                    coamz = coset(ax,ay,az-1)
                    coapx = coset(ax+1,ay,az)
                    coapy = coset(ax,ay+1,az)
                    coapz = coset(ax,ay,az+1)
                    DO lb=lb_start,lb_max_set
                      DO bx=0,lb
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          cob = coset(bx,by,bz)
                          t(coa,cob,cdax) = ftz*t(coapx,cob,cda) -&
                                            fax*t(coamx,cob,cda)
                          t(coa,cob,cday) = ftz*t(coapy,cob,cda) -&
                                            fay*t(coamy,cob,cda)
                          t(coa,cob,cdaz) = ftz*t(coapz,cob,cda) -&
                                            faz*t(coamz,cob,cda)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO

            END DO
          END DO
        END DO

!       *** Return all the calculated derivatives of   ***
!       *** the kinetic energy integrals, if requested ***

        IF (return_derivatives) THEN
          DO k=2,ncoset(da_max)
            jstart = (k - 1)*SIZE(tab,1)
            DO j=1,ncoset(lb_max_set)
              jk = jstart + j
              DO i=1,ncoset(la_max_set)
                tab(na+i,nb+jk) = t(i,j,k)
              END DO
            END DO
          END DO
        END IF

!       *** Calculate the force contribution for the atomic center a ***

        IF (calculate_force_a) THEN
          DO k=1,3
            DO j=ncoset(lb_min_set-1)+1,ncoset(lb_max_set)
              DO i=ncoset(la_min_set-1)+1,ncoset(la_max_set)
                force_a(k) = force_a(k) + pab(na+i,nb+j)*t(i,j,k+1)
              END DO
            END DO
          END DO
        END IF

        nb = nb + ncoset(lb_max_set)

      END DO

      na = na + ncoset(la_max_set)
      nda = nda + ncoset(la_max)

    END DO

  END SUBROUTINE kinetic

END MODULE ai_kinetic

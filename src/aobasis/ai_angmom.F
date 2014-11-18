!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the angular momentum integrals over
!>      Cartesian Gaussian-type functions.
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      none
!> \author JGH (20.02.2005)
! *****************************************************************************
MODULE ai_angmom

! ax,ay,az  : Angular momentum index numbers of orbital a.
! bx,by,bz  : Angular momentum index numbers of orbital b.
! coset     : Cartesian orbital set pointer.
! dab       : Distance between the atomic centers a and b.
! l{a,b}    : Angular momentum quantum number of shell a or b.
! l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
! l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
! rac       : Distance vector between the atomic center a and C.
! rbc       : Distance vector between the atomic center b and C.
! rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
! angab     : Shell set of angular momentum integrals.
! zet{a,b}  : Exponents of the Gaussian-type functions a or b.
! zetp      : Reciprocal of the sum of the exponents of orbital a and b.

  USE ai_derivatives,                  ONLY: dabdr
  USE ai_moments,                      ONLY: moment
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: angmom

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param rac ...
!> \param rbc ...
!> \param angab ...
! *****************************************************************************
  SUBROUTINE angmom(la_max,npgfa,zeta,rpgfa,la_min,&
                    lb_max,npgfb,zetb,rpgfb, rac,rbc,angab)

    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: la_min, lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb, rpgfb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac, rbc
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(OUT)                            :: angab

    INTEGER                                  :: ax, ay, az, bx, by, bz, i, &
                                                ipgf, j, jpgf, la, la_start, &
                                                lb, na, nb
    REAL(KIND=dp)                            :: dab, f0, f1, f2, f3, fx, fy, &
                                                fz, rab2, zetp
    REAL(KIND=dp), DIMENSION(3)              :: ac1, ac2, ac3, bc1, bc2, bc3, &
                                                rab, rap, rbp
    REAL(KIND=dp), DIMENSION(ncoset(la_max),&
      ncoset(lb_max))                        :: s
    REAL(KIND=dp), DIMENSION(ncoset(la_max),&
      ncoset(lb_max), 3)                     :: as

!   *** Calculate the distance between the centers a and b ***

    rab = rbc - rac
    rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
    dab = SQRT(rab2)

!   *** Loop over all pairs of primitive Gaussian-type functions ***
    angab = 0.0_dp
    s = 0.0_dp
    as = 0.0_dp

    na = 0

    DO ipgf=1,npgfa

      nb = 0

      DO jpgf=1,npgfb

        s=0.0_dp
        as=0.0_dp

!       *** Screening ***

        IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
          DO j=nb+1,nb+ncoset(lb_max)
            DO i=na+1,na+ncoset(la_max)
              angab(i,j,1) = 0.0_dp
              angab(i,j,2) = 0.0_dp
              angab(i,j,3) = 0.0_dp
            END DO
          END DO
          nb = nb + ncoset(lb_max)
          CYCLE
        END IF

!       *** Calculate some prefactors ***

        zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))

        f0 = (pi*zetp)**1.5_dp
        f1 = zetb(jpgf)*zetp
        f2 = 0.5_dp*zetp
        !
        bc1(1) =  0._dp
        bc1(2) = -f1*rbc(3)
        bc1(3) =  f1*rbc(2)
        !
        bc2(1) =  f1*rbc(3)
        bc2(2) =  0._dp
        bc2(3) = -f1*rbc(1)
        !
        bc3(1) = -f1*rbc(2)
        bc3(2) =  f1*rbc(1)
        bc3(3) =  0._dp
        !
        ac1(1) =  0._dp
        ac1(2) =  zeta(ipgf)*zetp*rac(3)
        ac1(3) = -zeta(ipgf)*zetp*rac(2)
        !
        ac2(1) = -zeta(ipgf)*zetp*rac(3)
        ac2(2) =  0._dp
        ac2(3) =  zeta(ipgf)*zetp*rac(1)
        !
        ac3(1) =  zeta(ipgf)*zetp*rac(2)
        ac3(2) = -zeta(ipgf)*zetp*rac(1)
        ac3(3) =  0._dp
        !
!       *** Calculate the basic two-center overlap integral [s|s] ***
!       *** Calculate the basic two-center angular momentum integral [s|L|s] ***

        s(1,1) = f0*EXP(-zeta(ipgf)*f1*rab2)
        as(1,1,1) = 2._dp*zeta(ipgf)*f1*(rac(2)*rbc(3)-rac(3)*rbc(2))*s(1,1)
        as(1,1,2) = 2._dp*zeta(ipgf)*f1*(rac(3)*rbc(1)-rac(1)*rbc(3))*s(1,1)
        as(1,1,3) = 2._dp*zeta(ipgf)*f1*(rac(1)*rbc(2)-rac(2)*rbc(1))*s(1,1)

!       *** Recurrence steps: [s|L|s] -> [a|L|b] ***

        IF (la_max > 0) THEN

!         *** Vertical recurrence steps: [s|L|s] -> [a|L|s] ***

          rap(:) = f1*rab(:)

!         *** [p|s] = (Pi - Ai)*[s|s]  (i = x,y,z) ***
!         *** [p|Ln|s] = (Pi - Ai)*[s|Ln|s]+xb/(xa+xb)*(1i x BC)*[s|s]  (i = x,y,z) ***

          s(2,1) = rap(1)*s(1,1)
          s(3,1) = rap(2)*s(1,1)
          s(4,1) = rap(3)*s(1,1)
          as(2,1,1) = rap(1)*as(1,1,1) + bc1(1)*s(1,1)
          as(2,1,2) = rap(1)*as(1,1,2) + bc1(2)*s(1,1)
          as(2,1,3) = rap(1)*as(1,1,3) + bc1(3)*s(1,1)
          as(3,1,1) = rap(2)*as(1,1,1) + bc2(1)*s(1,1)
          as(3,1,2) = rap(2)*as(1,1,2) + bc2(2)*s(1,1)
          as(3,1,3) = rap(2)*as(1,1,3) + bc2(3)*s(1,1)
          as(4,1,1) = rap(3)*as(1,1,1) + bc3(1)*s(1,1)
          as(4,1,2) = rap(3)*as(1,1,2) + bc3(2)*s(1,1)
          as(4,1,3) = rap(3)*as(1,1,3) + bc3(3)*s(1,1)

!         *** [a|s] = (Pi - Ai)*[a-1i|s] + f2*Ni(a-1i)*[a-2i|s]          ***
!         *** [a|Ln|s] = (Pi - Ai)*[a-1i|Ln|s] + f2*Ni(a-1i)*[a-2i|Ln|s] ***
!         ***           + xb/(xa+xb)*{(1i x BC)}_n*[a-1i|s]                  ***

          DO la=2,la_max

!           *** Increase the angular momentum component z of function a ***

            s(coset(0,0,la),1) = rap(3)*s(coset(0,0,la-1),1) +&
                                 f2*REAL(la-1,dp)*s(coset(0,0,la-2),1)
            as(coset(0,0,la),1,1) = rap(3)*as(coset(0,0,la-1),1,1) +&
                                 f2*REAL(la-1,dp)*as(coset(0,0,la-2),1,1) +&
                                 bc3(1)*s(coset(0,0,la-1),1)
            as(coset(0,0,la),1,2) = rap(3)*as(coset(0,0,la-1),1,2) +&
                                 f2*REAL(la-1,dp)*as(coset(0,0,la-2),1,2) +&
                                 bc3(2)*s(coset(0,0,la-1),1)
            as(coset(0,0,la),1,3) = rap(3)*as(coset(0,0,la-1),1,3) +&
                                 f2*REAL(la-1,dp)*as(coset(0,0,la-2),1,3) +&
                                 bc3(3)*s(coset(0,0,la-1),1)

!           *** Increase the angular momentum component y of function a ***

            az = la - 1
            s(coset(0,1,az),1) = rap(2)*s(coset(0,0,az),1)
            as(coset(0,1,az),1,1) = rap(2)*as(coset(0,0,az),1,1) +&
                                    bc2(1)*s(coset(0,0,az),1)
            as(coset(0,1,az),1,2) = rap(2)*as(coset(0,0,az),1,2) +&
                                    bc2(2)*s(coset(0,0,az),1)
            as(coset(0,1,az),1,3) = rap(2)*as(coset(0,0,az),1,3) +&
                                    bc2(3)*s(coset(0,0,az),1)

            DO ay=2,la
              az = la - ay
              s(coset(0,ay,az),1) = rap(2)*s(coset(0,ay-1,az),1) +&
                                    f2*REAL(ay-1,dp)*s(coset(0,ay-2,az),1)
              as(coset(0,ay,az),1,1) = rap(2)*as(coset(0,ay-1,az),1,1) +&
                                    f2*REAL(ay-1,dp)*as(coset(0,ay-2,az),1,1) +&
                                    bc2(1)*s(coset(0,ay-1,az),1)
              as(coset(0,ay,az),1,2) = rap(2)*as(coset(0,ay-1,az),1,2) +&
                                    f2*REAL(ay-1,dp)*as(coset(0,ay-2,az),1,2) +&
                                    bc2(2)*s(coset(0,ay-1,az),1)
              as(coset(0,ay,az),1,3) = rap(2)*as(coset(0,ay-1,az),1,3) +&
                                    f2*REAL(ay-1,dp)*as(coset(0,ay-2,az),1,3) +&
                                    bc2(3)*s(coset(0,ay-1,az),1)
            END DO

!           *** Increase the angular momentum component x of function a ***

            DO ay=0,la-1
              az = la - 1 - ay
              s(coset(1,ay,az),1) = rap(1)*s(coset(0,ay,az),1)
              as(coset(1,ay,az),1,1) = rap(1)*as(coset(0,ay,az),1,1) +&
                                       bc1(1)*s(coset(0,ay,az),1)
              as(coset(1,ay,az),1,2) = rap(1)*as(coset(0,ay,az),1,2) +&
                                       bc1(2)*s(coset(0,ay,az),1)
              as(coset(1,ay,az),1,3) = rap(1)*as(coset(0,ay,az),1,3) +&
                                       bc1(3)*s(coset(0,ay,az),1)
            END DO

            DO ax=2,la
              f3 = f2*REAL(ax-1,dp)
              DO ay=0,la-ax
                az = la - ax - ay
                s(coset(ax,ay,az),1) = rap(1)*s(coset(ax-1,ay,az),1) +&
                                       f3*s(coset(ax-2,ay,az),1)
                as(coset(ax,ay,az),1,1) = rap(1)*as(coset(ax-1,ay,az),1,1) +&
                                       f3*as(coset(ax-2,ay,az),1,1) +&
                                       bc1(1)*s(coset(ax-1,ay,az),1)
                as(coset(ax,ay,az),1,2) = rap(1)*as(coset(ax-1,ay,az),1,2) +&
                                       f3*as(coset(ax-2,ay,az),1,2) +&
                                       bc1(2)*s(coset(ax-1,ay,az),1)
                as(coset(ax,ay,az),1,3) = rap(1)*as(coset(ax-1,ay,az),1,3) +&
                                       f3*as(coset(ax-2,ay,az),1,3) +&
                                       bc1(3)*s(coset(ax-1,ay,az),1)
              END DO
            END DO

          END DO

!         *** Recurrence steps: [a|L|s] -> [a|L|b] ***

          IF (lb_max > 0) THEN

            DO j=2,ncoset(lb_max)
              DO i=1,ncoset(la_max)
                s(i,j) = 0.0_dp
                as(i,j,1) = 0.0_dp
                as(i,j,2) = 0.0_dp
                as(i,j,3) = 0.0_dp
              END DO
            END DO

!           *** Horizontal recurrence steps ***

            rbp(:) = rap(:) - rab(:)

!           *** [a|L|p] = [a+1i|Lm|s] - (Bi - Ai)*[a|Lm|s] ***
!           ***         + [a+1k|s] + (Ak - Ck)*[a|s]  eps(i,m,k)

            IF (lb_max == 1) THEN
              la_start = la_min
            ELSE
              la_start = MAX(0,la_min-1)
            END IF

            DO la=la_start,la_max-1
              DO ax=0,la
                DO ay=0,la-ax
                  az = la - ax - ay
                  s(coset(ax,ay,az),2) = s(coset(ax+1,ay,az),1) - &
                                                 rab(1)*s(coset(ax,ay,az),1)
                  s(coset(ax,ay,az),3) = s(coset(ax,ay+1,az),1) - &
                                                 rab(2)*s(coset(ax,ay,az),1)
                  s(coset(ax,ay,az),4) = s(coset(ax,ay,az+1),1) - &
                                                 rab(3)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),2,1) = as(coset(ax+1,ay,az),1,1) - &
                                                 rab(1)*as(coset(ax,ay,az),1,1)
                  as(coset(ax,ay,az),3,1) = as(coset(ax,ay+1,az),1,1) - &
                                                 rab(2)*as(coset(ax,ay,az),1,1) &
                                - s(coset(ax,ay,az+1),1) - rac(3)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),4,1) = as(coset(ax,ay,az+1),1,1) - &
                                                 rab(3)*as(coset(ax,ay,az),1,1) &
                                + s(coset(ax,ay+1,az),1) + rac(2)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),2,2) = as(coset(ax+1,ay,az),1,2) - &
                                                 rab(1)*as(coset(ax,ay,az),1,2) &
                                + s(coset(ax,ay,az+1),1) + rac(3)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),3,2) = as(coset(ax,ay+1,az),1,2) - &
                                                 rab(2)*as(coset(ax,ay,az),1,2)
                  as(coset(ax,ay,az),4,2) = as(coset(ax,ay,az+1),1,2) - &
                                                 rab(3)*as(coset(ax,ay,az),1,2) &
                                - s(coset(ax+1,ay,az),1) - rac(1)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),2,3) = as(coset(ax+1,ay,az),1,3) - &
                                                 rab(1)*as(coset(ax,ay,az),1,3) &
                                - s(coset(ax,ay+1,az),1) - rac(2)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),3,3) = as(coset(ax,ay+1,az),1,3) - &
                                                 rab(2)*as(coset(ax,ay,az),1,3) &
                                + s(coset(ax+1,ay,az),1) + rac(1)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),4,3) = as(coset(ax,ay,az+1),1,3) - &
                                                 rab(3)*as(coset(ax,ay,az),1,3)
                END DO
              END DO
            END DO

!           *** Vertical recurrence step ***

!           *** [a|p] = (Pi - Bi)*[a|s] + f2*Ni(a)*[a-1i|s]       ***
!           *** [a|L|p] = (Pi - Bi)*[a|L|s] + f2*Ni(a)*[a-1i|L|s] ***
!           ***           + xa/(xa+xb)*(AC x 1i)*[a|s]            ***
!           ***           + 0.5*zetp*SUM_k Nk(a)*(1k x 1i)*[a-1k|s]   ***

            DO ax=0,la_max
              fx = f2*REAL(ax,dp)
              DO ay=0,la_max-ax
                fy = f2*REAL(ay,dp)
                az = la_max - ax - ay
                fz = f2*REAL(az,dp)
                IF (ax == 0) THEN
                  s(coset(ax,ay,az),2) = rbp(1)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),2,1) = rbp(1)*as(coset(ax,ay,az),1,1) +&
                                            ac1(1)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),2,2) = rbp(1)*as(coset(ax,ay,az),1,2) +&
                                            ac1(2)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),2,3) = rbp(1)*as(coset(ax,ay,az),1,3) +&
                                            ac1(3)*s(coset(ax,ay,az),1)
                ELSE
                  s(coset(ax,ay,az),2) = rbp(1)*s(coset(ax,ay,az),1) +&
                                         fx*s(coset(ax-1,ay,az),1)
                  as(coset(ax,ay,az),2,1) = rbp(1)*as(coset(ax,ay,az),1,1) +&
                                         fx*as(coset(ax-1,ay,az),1,1) +&
                                         ac1(1)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),2,2) = rbp(1)*as(coset(ax,ay,az),1,2) +&
                                         fx*as(coset(ax-1,ay,az),1,2) +&
                                         ac1(2)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),2,3) = rbp(1)*as(coset(ax,ay,az),1,3) +&
                                         fx*as(coset(ax-1,ay,az),1,3) +&
                                         ac1(3)*s(coset(ax,ay,az),1)
                END IF
                IF ( az > 0 ) as(coset(ax,ay,az),2,2) = as(coset(ax,ay,az),2,2) +&
                              f2*REAL(az,dp)*s(coset(ax,ay,az-1),1)
                IF ( ay > 0 ) as(coset(ax,ay,az),2,3) = as(coset(ax,ay,az),2,3) -&
                              f2*REAL(ay,dp)*s(coset(ax,ay-1,az),1)
                IF (ay == 0) THEN
                  s(coset(ax,ay,az),3) = rbp(2)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),3,1) = rbp(2)*as(coset(ax,ay,az),1,1) +&
                                            ac2(1)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),3,2) = rbp(2)*as(coset(ax,ay,az),1,2) +&
                                            ac2(2)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),3,3) = rbp(2)*as(coset(ax,ay,az),1,3) +&
                                            ac2(3)*s(coset(ax,ay,az),1)
                ELSE
                  s(coset(ax,ay,az),3) = rbp(2)*s(coset(ax,ay,az),1) +&
                                         fy*s(coset(ax,ay-1,az),1)
                  as(coset(ax,ay,az),3,1) = rbp(2)*as(coset(ax,ay,az),1,1) +&
                                         fy*as(coset(ax,ay-1,az),1,1) +&
                                         ac2(1)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),3,2) = rbp(2)*as(coset(ax,ay,az),1,2) +&
                                         fy*as(coset(ax,ay-1,az),1,2) +&
                                         ac2(2)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),3,3) = rbp(2)*as(coset(ax,ay,az),1,3) +&
                                         fy*as(coset(ax,ay-1,az),1,3) +&
                                         ac2(3)*s(coset(ax,ay,az),1)
                END IF
                IF ( az > 0 ) as(coset(ax,ay,az),3,1) = as(coset(ax,ay,az),3,1) -&
                              f2*REAL(az,dp)*s(coset(ax,ay,az-1),1)
                IF ( ax > 0 ) as(coset(ax,ay,az),3,3) = as(coset(ax,ay,az),3,3) +&
                              f2*REAL(ax,dp)*s(coset(ax-1,ay,az),1)
                IF (az == 0) THEN
                  s(coset(ax,ay,az),4) = rbp(3)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),4,1) = rbp(3)*as(coset(ax,ay,az),1,1) +&
                                            ac3(1)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),4,2) = rbp(3)*as(coset(ax,ay,az),1,2) +&
                                            ac3(2)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),4,3) = rbp(3)*as(coset(ax,ay,az),1,3) +&
                                            ac3(3)*s(coset(ax,ay,az),1)
                ELSE
                  s(coset(ax,ay,az),4) = rbp(3)*s(coset(ax,ay,az),1) +&
                                         fz*s(coset(ax,ay,az-1),1)
                  as(coset(ax,ay,az),4,1) = rbp(3)*as(coset(ax,ay,az),1,1) +&
                                         fz*as(coset(ax,ay,az-1),1,1) +&
                                         ac3(1)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),4,2) = rbp(3)*as(coset(ax,ay,az),1,2) +&
                                         fz*as(coset(ax,ay,az-1),1,2) +&
                                         ac3(2)*s(coset(ax,ay,az),1)
                  as(coset(ax,ay,az),4,3) = rbp(3)*as(coset(ax,ay,az),1,3) +&
                                         fz*as(coset(ax,ay,az-1),1,3) +&
                                         ac3(3)*s(coset(ax,ay,az),1)
                END IF
                IF ( ay > 0 ) as(coset(ax,ay,az),4,1) = as(coset(ax,ay,az),4,1) +&
                              f2*REAL(ay,dp)*s(coset(ax,ay-1,az),1)
                IF ( ax > 0 ) as(coset(ax,ay,az),4,2) = as(coset(ax,ay,az),4,2) -&
                              f2*REAL(ax,dp)*s(coset(ax-1,ay,az),1)
              END DO
            END DO

!           *** Recurrence steps: [a|L|p] -> [a|L|b] ***

            DO lb=2,lb_max

!             *** Horizontal recurrence steps ***

!             *** [a|Lm|b] = [a+1i|Lm|b-1i] - (Bi - Ai)*[a|Lm|b-1i] ***
!             ***         + [a+1k|b-1i] + (Ak - Ck)*[a|b-1i]  eps(i,m,k)

              IF (lb == lb_max) THEN
                la_start = la_min
              ELSE
                la_start = MAX(0,la_min-1)
              END IF

              DO la=la_start,la_max-1
                DO ax=0,la
                  DO ay=0,la-ax
                    az = la - ax - ay

!                   *** Shift of angular momentum component z from a to b ***

                    s(coset(ax,ay,az),coset(0,0,lb)) =&
                      s(coset(ax,ay,az+1),coset(0,0,lb-1)) -&
                      rab(3)*s(coset(ax,ay,az),coset(0,0,lb-1))
                    as(coset(ax,ay,az),coset(0,0,lb),1) =&
                      as(coset(ax,ay,az+1),coset(0,0,lb-1),1) -&
                      rab(3)*as(coset(ax,ay,az),coset(0,0,lb-1),1) &
                      +s(coset(ax,ay+1,az),coset(0,0,lb-1)) &
                      +rac(2)*s(coset(ax,ay,az),coset(0,0,lb-1))
                    as(coset(ax,ay,az),coset(0,0,lb),2) =&
                      as(coset(ax,ay,az+1),coset(0,0,lb-1),2) -&
                      rab(3)*as(coset(ax,ay,az),coset(0,0,lb-1),2) &
                      -s(coset(ax+1,ay,az),coset(0,0,lb-1)) &
                      -rac(1)*s(coset(ax,ay,az),coset(0,0,lb-1))
                    as(coset(ax,ay,az),coset(0,0,lb),3) =&
                      as(coset(ax,ay,az+1),coset(0,0,lb-1),3) -&
                      rab(3)*as(coset(ax,ay,az),coset(0,0,lb-1),3)

!                   *** Shift of angular momentum component y from a to b ***

                    DO by=1,lb
                      bz = lb - by
                      s(coset(ax,ay,az),coset(0,by,bz)) =&
                        s(coset(ax,ay+1,az),coset(0,by-1,bz)) -&
                        rab(2)*s(coset(ax,ay,az),coset(0,by-1,bz))
                      as(coset(ax,ay,az),coset(0,by,bz),1) =&
                        as(coset(ax,ay+1,az),coset(0,by-1,bz),1) -&
                        rab(2)*as(coset(ax,ay,az),coset(0,by-1,bz),1) &
                        -s(coset(ax,ay,az+1),coset(0,by-1,bz)) &
                        -rac(3)*s(coset(ax,ay,az),coset(0,by-1,bz))
                      as(coset(ax,ay,az),coset(0,by,bz),2) =&
                        as(coset(ax,ay+1,az),coset(0,by-1,bz),2) -&
                        rab(2)*as(coset(ax,ay,az),coset(0,by-1,bz),2)
                      as(coset(ax,ay,az),coset(0,by,bz),3) =&
                        as(coset(ax,ay+1,az),coset(0,by-1,bz),3) -&
                        rab(2)*as(coset(ax,ay,az),coset(0,by-1,bz),3) &
                        +s(coset(ax+1,ay,az),coset(0,by-1,bz)) &
                        +rac(1)*s(coset(ax,ay,az),coset(0,by-1,bz))
                    END DO

!                   *** Shift of angular momentum component x from a to b ***

                    DO bx=1,lb
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        s(coset(ax,ay,az),coset(bx,by,bz)) =&
                          s(coset(ax+1,ay,az),coset(bx-1,by,bz)) -&
                          rab(1)*s(coset(ax,ay,az),coset(bx-1,by,bz))
                        as(coset(ax,ay,az),coset(bx,by,bz),1) =&
                          as(coset(ax+1,ay,az),coset(bx-1,by,bz),1) -&
                          rab(1)*as(coset(ax,ay,az),coset(bx-1,by,bz),1)
                        as(coset(ax,ay,az),coset(bx,by,bz),2) =&
                          as(coset(ax+1,ay,az),coset(bx-1,by,bz),2) -&
                          rab(1)*as(coset(ax,ay,az),coset(bx-1,by,bz),2) &
                          +s(coset(ax,ay,az+1),coset(bx-1,by,bz)) &
                          +rac(3)*s(coset(ax,ay,az),coset(bx-1,by,bz))
                        as(coset(ax,ay,az),coset(bx,by,bz),3) =&
                          as(coset(ax+1,ay,az),coset(bx-1,by,bz),3) -&
                          rab(1)*as(coset(ax,ay,az),coset(bx-1,by,bz),3) &
                          -s(coset(ax,ay+1,az),coset(bx-1,by,bz)) &
                          -rac(2)*s(coset(ax,ay,az),coset(bx-1,by,bz))
                      END DO
                    END DO

                  END DO
                END DO
              END DO

!             *** Vertical recurrence step ***

!             *** [a|b] = (Pi - Bi)*[a|b-1i] + f2*Ni(a)*[a-1i|b-1i] +       ***
!             ***         f2*Ni(b-1i)*[a|b-2i]                              ***
!             *** [a|L|b] = (Pi - Bi)*[a|L|b-1i] + f2*Ni(a)*[a-1i|L|b-1i] + ***
!             ***         f2*Ni(b-1i)*[a|L|b-2i]                            ***
!             ***         + xa/(xa+xb)*(AC x 1i)*[a|b-1i]                   ***
!             ***         + 0.5*zetp*SUM_k Nk(a)*(1k x 1i)*[a-1k|b-1i]      ***

              DO ax=0,la_max
                fx = f2*REAL(ax,dp)
                DO ay=0,la_max-ax
                  fy = f2*REAL(ay,dp)
                  az = la_max - ax - ay
                  fz = f2*REAL(az,dp)

!                 *** Increase the angular momentum component z of function b ***

                  f3 = f2*REAL(lb-1,dp)

                  IF (az == 0) THEN
                    s(coset(ax,ay,az),coset(0,0,lb)) =&
                      rbp(3)*s(coset(ax,ay,az),coset(0,0,lb-1)) +&
                      f3*s(coset(ax,ay,az),coset(0,0,lb-2))
                    as(coset(ax,ay,az),coset(0,0,lb),1) =&
                      rbp(3)*as(coset(ax,ay,az),coset(0,0,lb-1),1) +&
                      f3*as(coset(ax,ay,az),coset(0,0,lb-2),1) +&
                      ac3(1)*s(coset(ax,ay,az),coset(0,0,lb-1))
                    as(coset(ax,ay,az),coset(0,0,lb),2) =&
                      rbp(3)*as(coset(ax,ay,az),coset(0,0,lb-1),2) +&
                      f3*as(coset(ax,ay,az),coset(0,0,lb-2),2) +&
                      ac3(2)*s(coset(ax,ay,az),coset(0,0,lb-1))
                    as(coset(ax,ay,az),coset(0,0,lb),3) =&
                      rbp(3)*as(coset(ax,ay,az),coset(0,0,lb-1),3) +&
                      f3*as(coset(ax,ay,az),coset(0,0,lb-2),3) +&
                      ac3(3)*s(coset(ax,ay,az),coset(0,0,lb-1))
                  ELSE
                    s(coset(ax,ay,az),coset(0,0,lb)) =&
                      rbp(3)*s(coset(ax,ay,az),coset(0,0,lb-1)) +&
                      fz*s(coset(ax,ay,az-1),coset(0,0,lb-1)) +&
                      f3*s(coset(ax,ay,az),coset(0,0,lb-2))
                    as(coset(ax,ay,az),coset(0,0,lb),1) =&
                      rbp(3)*as(coset(ax,ay,az),coset(0,0,lb-1),1) +&
                      fz*as(coset(ax,ay,az-1),coset(0,0,lb-1),1) +&
                      f3*as(coset(ax,ay,az),coset(0,0,lb-2),1) +&
                      ac3(1)*s(coset(ax,ay,az),coset(0,0,lb-1))
                    as(coset(ax,ay,az),coset(0,0,lb),2) =&
                      rbp(3)*as(coset(ax,ay,az),coset(0,0,lb-1),2) +&
                      fz*as(coset(ax,ay,az-1),coset(0,0,lb-1),2) +&
                      f3*as(coset(ax,ay,az),coset(0,0,lb-2),2) +&
                      ac3(2)*s(coset(ax,ay,az),coset(0,0,lb-1))
                    as(coset(ax,ay,az),coset(0,0,lb),3) =&
                      rbp(3)*as(coset(ax,ay,az),coset(0,0,lb-1),3) +&
                      fz*as(coset(ax,ay,az-1),coset(0,0,lb-1),3) +&
                      f3*as(coset(ax,ay,az),coset(0,0,lb-2),3) +&
                      ac3(3)*s(coset(ax,ay,az),coset(0,0,lb-1))
                  END IF
                  IF ( ay > 0 ) as(coset(ax,ay,az),coset(0,0,lb),1) =&
                                as(coset(ax,ay,az),coset(0,0,lb),1) +&
                                f2*REAL(ay,dp)*s(coset(ax,ay-1,az),coset(0,0,lb-1))
                  IF ( ax > 0 ) as(coset(ax,ay,az),coset(0,0,lb),2) =&
                                as(coset(ax,ay,az),coset(0,0,lb),2) -&
                                f2*REAL(ax,dp)*s(coset(ax-1,ay,az),coset(0,0,lb-1))

!                 *** Increase the angular momentum component y of function b ***

                  IF (ay == 0) THEN
                    bz = lb - 1
                    s(coset(ax,ay,az),coset(0,1,bz)) =&
                      rbp(2)*s(coset(ax,ay,az),coset(0,0,bz))
                    as(coset(ax,ay,az),coset(0,1,bz),1) =&
                      rbp(2)*as(coset(ax,ay,az),coset(0,0,bz),1) +&
                      ac2(1)*s(coset(ax,ay,az),coset(0,0,bz))
                    as(coset(ax,ay,az),coset(0,1,bz),2) =&
                      rbp(2)*as(coset(ax,ay,az),coset(0,0,bz),2) +&
                      ac2(2)*s(coset(ax,ay,az),coset(0,0,bz))
                    as(coset(ax,ay,az),coset(0,1,bz),3) =&
                      rbp(2)*as(coset(ax,ay,az),coset(0,0,bz),3) +&
                      ac2(3)*s(coset(ax,ay,az),coset(0,0,bz))
                    IF ( az > 0 ) as(coset(ax,ay,az),coset(0,1,bz),1) =&
                                  as(coset(ax,ay,az),coset(0,1,bz),1) -&
                                  f2*REAL(az,dp)*s(coset(ax,ay,az-1),coset(0,0,bz))
                    IF ( ax > 0 ) as(coset(ax,ay,az),coset(0,1,bz),3) =&
                                  as(coset(ax,ay,az),coset(0,1,bz),3) +&
                                  f2*REAL(ax,dp)*s(coset(ax-1,ay,az),coset(0,0,bz))
                    DO by=2,lb
                      bz = lb - by
                      f3 = f2*REAL(by-1,dp)
                      s(coset(ax,ay,az),coset(0,by,bz)) =&
                        rbp(2)*s(coset(ax,ay,az),coset(0,by-1,bz)) +&
                        f3*s(coset(ax,ay,az),coset(0,by-2,bz))
                      as(coset(ax,ay,az),coset(0,by,bz),1) =&
                        rbp(2)*as(coset(ax,ay,az),coset(0,by-1,bz),1) +&
                        f3*as(coset(ax,ay,az),coset(0,by-2,bz),1) +&
                        ac2(1)*s(coset(ax,ay,az),coset(0,by-1,bz))
                      as(coset(ax,ay,az),coset(0,by,bz),2) =&
                        rbp(2)*as(coset(ax,ay,az),coset(0,by-1,bz),2) +&
                        f3*as(coset(ax,ay,az),coset(0,by-2,bz),2) +&
                        ac2(2)*s(coset(ax,ay,az),coset(0,by-1,bz))
                      as(coset(ax,ay,az),coset(0,by,bz),3) =&
                        rbp(2)*as(coset(ax,ay,az),coset(0,by-1,bz),3) +&
                        f3*as(coset(ax,ay,az),coset(0,by-2,bz),3) +&
                        ac2(3)*s(coset(ax,ay,az),coset(0,by-1,bz))
                      IF ( az > 0 ) as(coset(ax,ay,az),coset(0,by,bz),1) =&
                                    as(coset(ax,ay,az),coset(0,by,bz),1) -&
                                    f2*REAL(az,dp)*s(coset(ax,ay,az-1),coset(0,by-1,bz))
                      IF ( ax > 0 ) as(coset(ax,ay,az),coset(0,by,bz),3) =&
                                    as(coset(ax,ay,az),coset(0,by,bz),3) +&
                                    f2*REAL(ax,dp)*s(coset(ax-1,ay,az),coset(0,by-1,bz))
                    END DO
                  ELSE
                    bz = lb - 1
                    s(coset(ax,ay,az),coset(0,1,bz)) =&
                      rbp(2)*s(coset(ax,ay,az),coset(0,0,bz)) +&
                      fy*s(coset(ax,ay-1,az),coset(0,0,bz))
                    as(coset(ax,ay,az),coset(0,1,bz),1) =&
                      rbp(2)*as(coset(ax,ay,az),coset(0,0,bz),1) +&
                      fy*as(coset(ax,ay-1,az),coset(0,0,bz),1) +&
                      ac2(1)*s(coset(ax,ay,az),coset(0,0,bz))
                    as(coset(ax,ay,az),coset(0,1,bz),2) =&
                      rbp(2)*as(coset(ax,ay,az),coset(0,0,bz),2) +&
                      fy*as(coset(ax,ay-1,az),coset(0,0,bz),2) +&
                      ac2(2)*s(coset(ax,ay,az),coset(0,0,bz))
                    as(coset(ax,ay,az),coset(0,1,bz),3) =&
                      rbp(2)*as(coset(ax,ay,az),coset(0,0,bz),3) +&
                      fy*as(coset(ax,ay-1,az),coset(0,0,bz),3) +&
                      ac2(3)*s(coset(ax,ay,az),coset(0,0,bz))
                    IF ( az > 0 ) as(coset(ax,ay,az),coset(0,1,bz),1) =&
                                  as(coset(ax,ay,az),coset(0,1,bz),1) -&
                                  f2*REAL(az,dp)*s(coset(ax,ay,az-1),coset(0,0,bz))
                    IF ( ax > 0 ) as(coset(ax,ay,az),coset(0,1,bz),3) =&
                                  as(coset(ax,ay,az),coset(0,1,bz),3) +&
                                  f2*REAL(ax,dp)*s(coset(ax-1,ay,az),coset(0,0,bz))
                    DO by=2,lb
                      bz = lb - by
                      f3 = f2*REAL(by-1,dp)
                      s(coset(ax,ay,az),coset(0,by,bz)) =&
                        rbp(2)*s(coset(ax,ay,az),coset(0,by-1,bz)) +&
                        fy*s(coset(ax,ay-1,az),coset(0,by-1,bz)) +&
                        f3*s(coset(ax,ay,az),coset(0,by-2,bz))
                      as(coset(ax,ay,az),coset(0,by,bz),1) =&
                        rbp(2)*as(coset(ax,ay,az),coset(0,by-1,bz),1) +&
                        fy*as(coset(ax,ay-1,az),coset(0,by-1,bz),1) +&
                        f3*as(coset(ax,ay,az),coset(0,by-2,bz),1) +&
                        ac2(1)*s(coset(ax,ay,az),coset(0,by-1,bz))
                      as(coset(ax,ay,az),coset(0,by,bz),2) =&
                        rbp(2)*as(coset(ax,ay,az),coset(0,by-1,bz),2) +&
                        fy*as(coset(ax,ay-1,az),coset(0,by-1,bz),2) +&
                        f3*as(coset(ax,ay,az),coset(0,by-2,bz),2) +&
                        ac2(2)*s(coset(ax,ay,az),coset(0,by-1,bz))
                      as(coset(ax,ay,az),coset(0,by,bz),3) =&
                        rbp(2)*as(coset(ax,ay,az),coset(0,by-1,bz),3) +&
                        fy*as(coset(ax,ay-1,az),coset(0,by-1,bz),3) +&
                        f3*as(coset(ax,ay,az),coset(0,by-2,bz),3) +&
                        ac2(3)*s(coset(ax,ay,az),coset(0,by-1,bz))
                      IF ( az > 0 ) as(coset(ax,ay,az),coset(0,by,bz),1) =&
                                    as(coset(ax,ay,az),coset(0,by,bz),1) -&
                                    f2*REAL(az,dp)*s(coset(ax,ay,az-1),coset(0,by-1,bz))
                      IF ( ax > 0 ) as(coset(ax,ay,az),coset(0,by,bz),3) =&
                                    as(coset(ax,ay,az),coset(0,by,bz),3) +&
                                    f2*REAL(ax,dp)*s(coset(ax-1,ay,az),coset(0,by-1,bz))
                    END DO
                  END IF

!                 *** Increase the angular momentum component x of function b ***

                  IF (ax == 0) THEN
                    DO by=0,lb-1
                      bz = lb - 1 - by
                      s(coset(ax,ay,az),coset(1,by,bz)) =&
                        rbp(1)*s(coset(ax,ay,az),coset(0,by,bz))
                      as(coset(ax,ay,az),coset(1,by,bz),1) =&
                        rbp(1)*as(coset(ax,ay,az),coset(0,by,bz),1) +&
                        ac1(1)*s(coset(ax,ay,az),coset(0,by,bz))
                      as(coset(ax,ay,az),coset(1,by,bz),2) =&
                        rbp(1)*as(coset(ax,ay,az),coset(0,by,bz),2) +&
                        ac1(2)*s(coset(ax,ay,az),coset(0,by,bz))
                      as(coset(ax,ay,az),coset(1,by,bz),3) =&
                        rbp(1)*as(coset(ax,ay,az),coset(0,by,bz),3) +&
                        ac1(3)*s(coset(ax,ay,az),coset(0,by,bz))
                      IF ( az > 0 ) as(coset(ax,ay,az),coset(1,by,bz),2) =&
                                    as(coset(ax,ay,az),coset(1,by,bz),2) +&
                                    f2*REAL(az,dp)*s(coset(ax,ay,az-1),coset(0,by,bz))
                      IF ( ay > 0 ) as(coset(ax,ay,az),coset(1,by,bz),3) =&
                                    as(coset(ax,ay,az),coset(1,by,bz),3) -&
                                    f2*REAL(ay,dp)*s(coset(ax,ay-1,az),coset(0,by,bz))
                    END DO
                    DO bx=2,lb
                      f3 = f2*REAL(bx-1,dp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        s(coset(ax,ay,az),coset(bx,by,bz)) =&
                          rbp(1)*s(coset(ax,ay,az),coset(bx-1,by,bz)) +&
                          f3*s(coset(ax,ay,az),coset(bx-2,by,bz))
                        as(coset(ax,ay,az),coset(bx,by,bz),1) =&
                          rbp(1)*as(coset(ax,ay,az),coset(bx-1,by,bz),1) +&
                          f3*as(coset(ax,ay,az),coset(bx-2,by,bz),1) +&
                          ac1(1)*s(coset(ax,ay,az),coset(bx-1,by,bz))
                        as(coset(ax,ay,az),coset(bx,by,bz),2) =&
                          rbp(1)*as(coset(ax,ay,az),coset(bx-1,by,bz),2) +&
                          f3*as(coset(ax,ay,az),coset(bx-2,by,bz),2) +&
                          ac1(2)*s(coset(ax,ay,az),coset(bx-1,by,bz))
                        as(coset(ax,ay,az),coset(bx,by,bz),3) =&
                          rbp(1)*as(coset(ax,ay,az),coset(bx-1,by,bz),3) +&
                          f3*as(coset(ax,ay,az),coset(bx-2,by,bz),3) +&
                          ac1(3)*s(coset(ax,ay,az),coset(bx-1,by,bz))
                        IF ( az > 0 ) as(coset(ax,ay,az),coset(bx,by,bz),2) =&
                                      as(coset(ax,ay,az),coset(bx,by,bz),2) +&
                                      f2*REAL(az,dp)*s(coset(ax,ay,az-1),coset(bx-1,by,bz))
                        IF ( ay > 0 ) as(coset(ax,ay,az),coset(bx,by,bz),3) =&
                                      as(coset(ax,ay,az),coset(bx,by,bz),3) -&
                                      f2*REAL(ay,dp)*s(coset(ax,ay-1,az),coset(bx-1,by,bz))
                      END DO
                    END DO
                  ELSE
                    DO by=0,lb-1
                      bz = lb - 1 - by
                      s(coset(ax,ay,az),coset(1,by,bz)) =&
                        rbp(1)*s(coset(ax,ay,az),coset(0,by,bz)) +&
                        fx*s(coset(ax-1,ay,az),coset(0,by,bz))
                      as(coset(ax,ay,az),coset(1,by,bz),1) =&
                        rbp(1)*as(coset(ax,ay,az),coset(0,by,bz),1) +&
                        fx*as(coset(ax-1,ay,az),coset(0,by,bz),1) +&
                        ac1(1)*s(coset(ax,ay,az),coset(0,by,bz))
                      as(coset(ax,ay,az),coset(1,by,bz),2) =&
                        rbp(1)*as(coset(ax,ay,az),coset(0,by,bz),2) +&
                        fx*as(coset(ax-1,ay,az),coset(0,by,bz),2) +&
                        ac1(2)*s(coset(ax,ay,az),coset(0,by,bz))
                      as(coset(ax,ay,az),coset(1,by,bz),3) =&
                        rbp(1)*as(coset(ax,ay,az),coset(0,by,bz),3) +&
                        fx*as(coset(ax-1,ay,az),coset(0,by,bz),3) +&
                        ac1(3)*s(coset(ax,ay,az),coset(0,by,bz))
                      IF ( az > 0 ) as(coset(ax,ay,az),coset(1,by,bz),2) =&
                                    as(coset(ax,ay,az),coset(1,by,bz),2) +&
                                    f2*REAL(az,dp)*s(coset(ax,ay,az-1),coset(0,by,bz))
                      IF ( ay > 0 ) as(coset(ax,ay,az),coset(1,by,bz),3) =&
                                    as(coset(ax,ay,az),coset(1,by,bz),3) -&
                                    f2*REAL(ay,dp)*s(coset(ax,ay-1,az),coset(0,by,bz))
                    END DO
                    DO bx=2,lb
                      f3 = f2*REAL(bx-1,dp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        s(coset(ax,ay,az),coset(bx,by,bz)) =&
                          rbp(1)*s(coset(ax,ay,az),coset(bx-1,by,bz)) +&
                          fx*s(coset(ax-1,ay,az),coset(bx-1,by,bz)) +&
                          f3*s(coset(ax,ay,az),coset(bx-2,by,bz))
                        as(coset(ax,ay,az),coset(bx,by,bz),1) =&
                          rbp(1)*as(coset(ax,ay,az),coset(bx-1,by,bz),1) +&
                          fx*as(coset(ax-1,ay,az),coset(bx-1,by,bz),1) +&
                          f3*as(coset(ax,ay,az),coset(bx-2,by,bz),1) +&
                          ac1(1)*s(coset(ax,ay,az),coset(bx-1,by,bz))
                        as(coset(ax,ay,az),coset(bx,by,bz),2) =&
                          rbp(1)*as(coset(ax,ay,az),coset(bx-1,by,bz),2) +&
                          fx*as(coset(ax-1,ay,az),coset(bx-1,by,bz),2) +&
                          f3*as(coset(ax,ay,az),coset(bx-2,by,bz),2) +&
                          ac1(2)*s(coset(ax,ay,az),coset(bx-1,by,bz))
                        as(coset(ax,ay,az),coset(bx,by,bz),3) =&
                          rbp(1)*as(coset(ax,ay,az),coset(bx-1,by,bz),3) +&
                          fx*as(coset(ax-1,ay,az),coset(bx-1,by,bz),3) +&
                          f3*as(coset(ax,ay,az),coset(bx-2,by,bz),3) +&
                          ac1(3)*s(coset(ax,ay,az),coset(bx-1,by,bz))
                      IF ( az > 0 ) as(coset(ax,ay,az),coset(bx,by,bz),2) =&
                                    as(coset(ax,ay,az),coset(bx,by,bz),2) +&
                                    f2*REAL(az,dp)*s(coset(ax,ay,az-1),coset(bx-1,by,bz))
                      IF ( ay > 0 ) as(coset(ax,ay,az),coset(bx,by,bz),3) =&
                                    as(coset(ax,ay,az),coset(bx,by,bz),3) -&
                                    f2*REAL(ay,dp)*s(coset(ax,ay-1,az),coset(bx-1,by,bz))
                      END DO
                    END DO
                  END IF

                END DO
              END DO

            END DO

          END IF

        ELSE

          IF (lb_max > 0) THEN

!           *** Vertical recurrence steps: [s|L|s] -> [s|L|b] ***

            rbp(:) = (f1 - 1.0_dp)*rab(:)

!           *** [s|p] = (Pi - Bi)*[s|s]                                  ***
!           *** [s|L|p] = (Pi - Bi)*[s|L|s] + xa/(xa+xb)*(AC x 1i)*[s|s] ***

            s(1,2) = rbp(1)*s(1,1)
            s(1,3) = rbp(2)*s(1,1)
            s(1,4) = rbp(3)*s(1,1)
            as(1,2,1) = rbp(1)*as(1,1,1) + ac1(1)*s(1,1)
            as(1,2,2) = rbp(1)*as(1,1,2) + ac1(2)*s(1,1)
            as(1,2,3) = rbp(1)*as(1,1,3) + ac1(3)*s(1,1)
            as(1,3,1) = rbp(2)*as(1,1,1) + ac2(1)*s(1,1)
            as(1,3,2) = rbp(2)*as(1,1,2) + ac2(2)*s(1,1)
            as(1,3,3) = rbp(2)*as(1,1,3) + ac2(3)*s(1,1)
            as(1,4,1) = rbp(3)*as(1,1,1) + ac3(1)*s(1,1)
            as(1,4,2) = rbp(3)*as(1,1,2) + ac3(2)*s(1,1)
            as(1,4,3) = rbp(3)*as(1,1,3) + ac3(3)*s(1,1)

!           *** [s|b] = (Pi - Bi)*[s|b-1i] + f2*Ni(b-1i)*[s|b-2i]       ***
!           *** [s|L|b] = (Pi - Bi)*[s|L|b-1i] + f2*Ni(b-1i)*[s|L|b-2i] ***
!           ***           + xa/(xa+xb)*(AC x 1i)*[s|s-1i]               ***

            DO lb=2,lb_max

!             *** Increase the angular momentum component z of function b ***

              s(1,coset(0,0,lb)) = rbp(3)*s(1,coset(0,0,lb-1)) +&
                                   f2*REAL(lb-1,dp)*s(1,coset(0,0,lb-2))
              as(1,coset(0,0,lb),1) = rbp(3)*as(1,coset(0,0,lb-1),1) +&
                                   f2*REAL(lb-1,dp)*as(1,coset(0,0,lb-2),1) +&
                                   ac3(1)*s(1,coset(0,0,lb-1))
              as(1,coset(0,0,lb),2) = rbp(3)*as(1,coset(0,0,lb-1),2) +&
                                   f2*REAL(lb-1,dp)*as(1,coset(0,0,lb-2),2) +&
                                   ac3(2)*s(1,coset(0,0,lb-1))
              as(1,coset(0,0,lb),3) = rbp(3)*as(1,coset(0,0,lb-1),3) +&
                                   f2*REAL(lb-1,dp)*as(1,coset(0,0,lb-2),3) +&
                                   ac3(3)*s(1,coset(0,0,lb-1))

!             *** Increase the angular momentum component y of function b ***

              bz = lb - 1
              s(1,coset(0,1,bz)) = rbp(2)*s(1,coset(0,0,bz))
              as(1,coset(0,1,bz),1) = rbp(2)*as(1,coset(0,0,bz),1) +&
                                      ac2(1)*s(1,coset(0,0,bz))
              as(1,coset(0,1,bz),2) = rbp(2)*as(1,coset(0,0,bz),2) +&
                                      ac2(2)*s(1,coset(0,0,bz))
              as(1,coset(0,1,bz),3) = rbp(2)*as(1,coset(0,0,bz),3) +&
                                      ac2(3)*s(1,coset(0,0,bz))

              DO by=2,lb
                bz = lb - by
                s(1,coset(0,by,bz)) = rbp(2)*s(1,coset(0,by-1,bz)) +&
                                      f2*REAL(by-1,dp)*s(1,coset(0,by-2,bz))
                as(1,coset(0,by,bz),1) = rbp(2)*as(1,coset(0,by-1,bz),1) +&
                                      f2*REAL(by-1,dp)*as(1,coset(0,by-2,bz),1) +&
                                      ac2(1)*s(1,coset(0,by-1,bz))
                as(1,coset(0,by,bz),2) = rbp(2)*as(1,coset(0,by-1,bz),2) +&
                                      f2*REAL(by-1,dp)*as(1,coset(0,by-2,bz),2) +&
                                      ac2(2)*s(1,coset(0,by-1,bz))
                as(1,coset(0,by,bz),3) = rbp(2)*as(1,coset(0,by-1,bz),3) +&
                                      f2*REAL(by-1,dp)*as(1,coset(0,by-2,bz),3) +&
                                      ac2(3)*s(1,coset(0,by-1,bz))
              END DO

!             *** Increase the angular momentum component x of function b ***

              DO by=0,lb-1
                bz = lb - 1 - by
                s(1,coset(1,by,bz)) = rbp(1)*s(1,coset(0,by,bz))
                as(1,coset(1,by,bz),1) = rbp(1)*as(1,coset(0,by,bz),1) +&
                                         ac1(1)*s(1,coset(0,by,bz))
                as(1,coset(1,by,bz),2) = rbp(1)*as(1,coset(0,by,bz),2) +&
                                         ac1(2)*s(1,coset(0,by,bz))
                as(1,coset(1,by,bz),3) = rbp(1)*as(1,coset(0,by,bz),3) +&
                                         ac1(3)*s(1,coset(0,by,bz))
              END DO

              DO bx=2,lb
                f3 = f2*REAL(bx-1,dp)
                DO by=0,lb-bx
                  bz = lb - bx - by
                  s(1,coset(bx,by,bz)) = rbp(1)*s(1,coset(bx-1,by,bz)) +&
                                         f3*s(1,coset(bx-2,by,bz))
                  as(1,coset(bx,by,bz),1) = rbp(1)*as(1,coset(bx-1,by,bz),1) +&
                                         f3*as(1,coset(bx-2,by,bz),1) +&
                                         ac1(1)*s(1,coset(bx-1,by,bz))
                  as(1,coset(bx,by,bz),2) = rbp(1)*as(1,coset(bx-1,by,bz),2) +&
                                         f3*as(1,coset(bx-2,by,bz),2) +&
                                         ac1(2)*s(1,coset(bx-1,by,bz))
                  as(1,coset(bx,by,bz),3) = rbp(1)*as(1,coset(bx-1,by,bz),3) +&
                                         f3*as(1,coset(bx-2,by,bz),3) +&
                                         ac1(3)*s(1,coset(bx-1,by,bz))
                END DO
              END DO

            END DO

          END IF

        END IF

        DO j=1,ncoset(lb_max)
          DO i=1,ncoset(la_max)
            angab(na+i,nb+j,1) = as(i,j,1)
            angab(na+i,nb+j,2) = as(i,j,2)
            angab(na+i,nb+j,3) = as(i,j,3)
          END DO
        END DO

        nb = nb + ncoset(lb_max)

      END DO

      na = na + ncoset(la_max)

    END DO

  END SUBROUTINE angmom

! *****************************************************************************
!> \brief ...
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param rac ...
!> \param rbc ...
!> \param angab ...
! *****************************************************************************
  SUBROUTINE angmom2(la_max,npgfa,zeta,rpgfa,la_min,&
                     lb_max,npgfb,zetb,rpgfb,lb_min,&
                     rac,rbc,angab)

    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: la_min, lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb, rpgfb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac, rbc
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: angab

    INTEGER                                  :: i, ia, ib, ja, jb, lda_min, &
                                                ldb_min, lds, lmax, na, nb
    REAL(KIND=dp)                            :: dab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: mab
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: dmab

    dab = SQRT( SUM ((rac - rbc)**2) )

    lda_min = MAX ( 0, la_min-1 )
    ldb_min = MAX ( 0, lb_min-1 )
    lmax = MAX(la_max+1,lb_max+1)
    lds = ncoset(lmax+1)*MAX(npgfa,npgfb)
    ALLOCATE(mab(lds,lds,3))
    ALLOCATE(dmab(lds,lds,3,3))
    mab = 0.0_dp
    CALL moment(lb_max+1,npgfb,zetb,rpgfb,ldb_min,&
                la_max,npgfa,zeta,rpgfa,&
                1,rbc,rac,mab)

    DO i=1,3
      CALL dabdr(lb_max,npgfb,zetb,rpgfb,lb_min,&
                 la_max,npgfa,rpgfa,la_min,&
                 dab,mab(:,:,i),dmab(:,:,1,i),dmab(:,:,2,i),dmab(:,:,3,i))
    END DO

    na=0
    DO ia=1,npgfa
      nb=0
      DO ib=1,npgfb
        DO ja=na+1,na+ncoset(la_max)
          DO jb=nb+1,nb+ncoset(lb_max)
            angab(ja,jb,1) = -dmab(jb,ja,3,2) + dmab(jb,ja,2,3)
            angab(ja,jb,2) = -dmab(jb,ja,1,3) + dmab(jb,ja,3,1)
            angab(ja,jb,3) = -dmab(jb,ja,2,1) + dmab(jb,ja,1,2)
          ENDDO
        ENDDO
        nb = nb + ncoset(lb_max)
      END DO
      na = na + ncoset(la_max)
    END DO

    DEALLOCATE(mab,dmab)

  END SUBROUTINE angmom2

END MODULE ai_angmom

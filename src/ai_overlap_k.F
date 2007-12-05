!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/ai_overlap_k [1.0] *
!!
!!   NAME
!!     ai_overlap_k
!!
!!   FUNCTION
!!     Calculation of the k-point overlap integrals over 
!!     Cartesian Gaussian-type functions.
!!
!!   AUTHOR
!!     Juerg Hutter (19.02.2007)
!!
!!   LITERATURE
!!     S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE ai_overlap_k

! *****************************************************************************

! ax,ay,az  : Angular momentum index numbers of orbital a.
! bx,by,bz  : Angular momentum index numbers of orbital b.
! coset     : Cartesian orbital set pointer.
! dab       : Distance between the atomic centers a and b.
! l{a,b}    : Angular momentum quantum number of shell a or b.
! l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
! l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
! rab       : Distance vector between the atomic centers a and b.
! rab2      : Square of the distance between the atomic centers a and b.
! rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
! sab       : Shell set of overlap integrals.
! zet{a,b}  : Exponents of the Gaussian-type functions a or b.
! zetp      : Reciprocal of the sum of the exponents of orbital a and b.
! kpoint    : k-point vector.

! *****************************************************************************

  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: gaussi,&
                                             pi
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: overlap_k

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE overlap_k(la_max,zeta,rpgfa,la_min,&
                       lb_max,zetb,rpgfb,lb_min,&
                       kpoint,rab,rab2,sabk)

!------------------------------------------------------------------------------

    INTEGER, INTENT(IN)                      :: la_max
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: la_min, lb_max
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb, rpgfb
    INTEGER, INTENT(IN)                      :: lb_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: kpoint, rab
    REAL(KIND=dp), INTENT(IN)                :: rab2
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: sabk

    INTEGER                                  :: ax, ay, az, bx, by, bz, i, &
                                                ipgf, j, jpgf, la, la_start, &
                                                lb, na, nb
    REAL(KIND=dp)                            :: dab, f0, f1, f2, f3, fx, fy, &
                                                fz, k2, zetp
    REAL(KIND=dp), DIMENSION(3)              :: rap, rbp
    REAL(KIND=dp), DIMENSION(ncoset(la_max),&
      ncoset(lb_max))                        :: s

!------------------------------------------------------------------------------

    dab = SQRT(rab2)
    k2 = SUM(kpoint(:)**2)

!   Loop over all pairs of primitive Gaussian-type functions

    na = 0

    DO ipgf=1,SIZE(zeta)

      nb = 0

      DO jpgf=1,SIZE(zetb)

!       Screening

        IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
          DO j=nb+1,nb+ncoset(lb_max)
            DO i=na+1,na+ncoset(la_max)
              sabk(i,j) = CMPLX(0.0_dp,0.0_dp,KIND=dp)
            END DO
          END DO
          nb = nb + ncoset(lb_max)
          CYCLE
        END IF

!       Calculate some prefactors

        zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))

        f0 = (pi*zetp)**1.5_dp
        f1 = zetb(jpgf)*zetp
        f2 = 0.5_dp*zetp

!       Calculate the basic overlap integral [s|k|s]

        s(1,1) = f0*EXP(-zeta(ipgf)*f1*rab2)*EXP(-0.25_dp*k2*zetp)

!       Recurrence steps: [s|k|s] -> [a|k|b]

        IF (la_max > 0) THEN

!         Vertical recurrence steps: [s|k|s] -> [a|k|s]

          rap(:) = f1*rab(:)

!         [p|k|s] = (Pi - Ai)*[s|k|s] + f2*I*ki*[s|k|s] (i = x,y,z)

          s(2,1) = rap(1)*s(1,1) + f2*kpoint(1)*s(1,1)*gaussi
          s(3,1) = rap(2)*s(1,1) + f2*kpoint(2)*s(1,1)*gaussi
          s(4,1) = rap(3)*s(1,1) + f2*kpoint(3)*s(1,1)*gaussi

!         [a|k|s] = (Pi - Ai)*[a-1i|k|s] + f2*I*ki*[a-1i|k|s] 
!                   + f2*Ni(a-1i)*[a-2i|k|s]

          DO la=2,la_max

!           Increase the angular momentum component z of function a

            s(coset(0,0,la),1) = rap(3)*s(coset(0,0,la-1),1) +&
                                 f2*kpoint(3)*s(coset(0,0,la-1),1)*gaussi +&
                                 f2*REAL(la-1,dp)*s(coset(0,0,la-2),1)

!           Increase the angular momentum component y of function a

            az = la - 1
            s(coset(0,1,az),1) = rap(2)*s(coset(0,0,az),1) +&
                                 f2*kpoint(2)*s(coset(0,0,az),1)*gaussi

            DO ay=2,la
              az = la - ay
              s(coset(0,ay,az),1) = rap(2)*s(coset(0,ay-1,az),1) +&
                                    f2*kpoint(2)*s(coset(0,ay-1,az),1)*gaussi +&
                                    f2*REAL(ay-1,dp)*s(coset(0,ay-2,az),1)
            END DO

!           Increase the angular momentum component x of function a

            DO ay=0,la-1
              az = la - 1 - ay
              s(coset(1,ay,az),1) = rap(1)*s(coset(0,ay,az),1) +&
                                    f2*kpoint(1)*s(coset(0,ay,az),1)*gaussi
            END DO

            DO ax=2,la
              f3 = f2*REAL(ax-1,dp)
              DO ay=0,la-ax
                az = la - ax - ay
                s(coset(ax,ay,az),1) = rap(1)*s(coset(ax-1,ay,az),1) +&
                                       f2*kpoint(1)*s(coset(ax-1,ay,az),1)*gaussi +&
                                       f3*s(coset(ax-2,ay,az),1)
              END DO
            END DO

          END DO

!         Recurrence steps: [a|k|s] -> [a|k|b]

          IF (lb_max > 0) THEN

            DO j=2,ncoset(lb_max)
              DO i=1,ncoset(la_max)
                s(i,j) = CMPLX(0.0_dp,0.0_dp,KIND=dp)
              END DO
            END DO

!           Horizontal recurrence steps

            rbp(:) = rap(:) - rab(:)

!           [a|k|p] = [a+1i|k|s] - (Bi - Ai)*[a|k|s]

            IF (lb_max == 1) THEN
              la_start = la_min
            ELSE
              la_start = MAX(0,la_min-1)
            END IF

            DO la=la_start,la_max-1
              DO ax=0,la
                DO ay=0,la-ax
                  az = la - ax - ay
                  s(coset(ax,ay,az),2) = s(coset(ax+1,ay,az),1) -&
                                         rab(1)*s(coset(ax,ay,az),1)
                  s(coset(ax,ay,az),3) = s(coset(ax,ay+1,az),1) -&
                                         rab(2)*s(coset(ax,ay,az),1)
                  s(coset(ax,ay,az),4) = s(coset(ax,ay,az+1),1) -&
                                         rab(3)*s(coset(ax,ay,az),1)
                END DO
              END DO
            END DO

!           Vertical recurrence step

!           [a|k|p] = (Pi - Bi)*[a|k|s] + f2*ki*I*[a|k|s] + f2*Ni(a)*[a-1i|k|s]

            DO ax=0,la_max
              fx = f2*REAL(ax,dp)
              DO ay=0,la_max-ax
                fy = f2*REAL(ay,dp)
                az = la_max - ax - ay
                fz = f2*REAL(az,dp)
                IF (ax == 0) THEN
                  s(coset(ax,ay,az),2) = rbp(1)*s(coset(ax,ay,az),1) +&
                                         f2*kpoint(1)*s(coset(ax,ay,az),1)*gaussi
                ELSE
                  s(coset(ax,ay,az),2) = rbp(1)*s(coset(ax,ay,az),1) +&
                                         f2*kpoint(1)*s(coset(ax,ay,az),1)*gaussi +&
                                         fx*s(coset(ax-1,ay,az),1)
                END IF
                IF (ay == 0) THEN
                  s(coset(ax,ay,az),3) = rbp(2)*s(coset(ax,ay,az),1) +&
                                         f2*kpoint(2)*s(coset(ax,ay,az),1)*gaussi
                ELSE
                  s(coset(ax,ay,az),3) = rbp(2)*s(coset(ax,ay,az),1) +&
                                         f2*kpoint(2)*s(coset(ax,ay,az),1)*gaussi +&
                                         fy*s(coset(ax,ay-1,az),1)
                END IF
                IF (az == 0) THEN
                  s(coset(ax,ay,az),4) = rbp(3)*s(coset(ax,ay,az),1) +&
                                         f2*kpoint(3)*s(coset(ax,ay,az),1)*gaussi
                ELSE
                  s(coset(ax,ay,az),4) = rbp(3)*s(coset(ax,ay,az),1) +&
                                         f2*kpoint(3)*s(coset(ax,ay,az),1)*gaussi +&
                                         fz*s(coset(ax,ay,az-1),1)
                END IF
              END DO
            END DO

!           Recurrence steps: [a|k|p] -> [a|k|b]

            DO lb=2,lb_max

!             Horizontal recurrence steps

!             [a|k|b] = [a+1i|k|b-1i] - (Bi - Ai)*[a|k|b-1i]

              IF (lb == lb_max) THEN
                la_start = la_min
              ELSE
                la_start = MAX(0,la_min-1)
              END IF

              DO la=la_start,la_max-1
                DO ax=0,la
                  DO ay=0,la-ax
                    az = la - ax - ay

!                   Shift of angular momentum component z from a to b

                    s(coset(ax,ay,az),coset(0,0,lb)) =&
                      s(coset(ax,ay,az+1),coset(0,0,lb-1)) -&
                      rab(3)*s(coset(ax,ay,az),coset(0,0,lb-1))

!                   Shift of angular momentum component y from a to b

                    DO by=1,lb
                      bz = lb - by
                      s(coset(ax,ay,az),coset(0,by,bz)) =&
                        s(coset(ax,ay+1,az),coset(0,by-1,bz)) -&
                        rab(2)*s(coset(ax,ay,az),coset(0,by-1,bz))
                    END DO

!                   Shift of angular momentum component x from a to b

                    DO bx=1,lb
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        s(coset(ax,ay,az),coset(bx,by,bz)) =&
                          s(coset(ax+1,ay,az),coset(bx-1,by,bz)) -&
                          rab(1)*s(coset(ax,ay,az),coset(bx-1,by,bz))
                      END DO
                    END DO

                  END DO
                END DO
              END DO

!             Vertical recurrence step

!             [a|k|b] = (Pi - Bi)*[a|k|b-1i] + f2*ki*I*[a|k|b-1i] 
!                       + f2*Ni(a)*[a-1i|k|b-1i] + f2*Ni(b-1i)*[a|k|b-2i]

              DO ax=0,la_max
                fx = f2*REAL(ax,dp)
                DO ay=0,la_max-ax
                  fy = f2*REAL(ay,dp)
                  az = la_max - ax - ay
                  fz = f2*REAL(az,dp)

!                 Increase the angular momentum component z of function b

                  f3 = f2*REAL(lb-1,dp)

                  IF (az == 0) THEN
                    s(coset(ax,ay,az),coset(0,0,lb)) =&
                      rbp(3)*s(coset(ax,ay,az),coset(0,0,lb-1)) +&
                      f2*kpoint(3)*s(coset(ax,ay,az),coset(0,0,lb-1))*gaussi +&
                      f3*s(coset(ax,ay,az),coset(0,0,lb-2))
                  ELSE
                    s(coset(ax,ay,az),coset(0,0,lb)) =&
                      rbp(3)*s(coset(ax,ay,az),coset(0,0,lb-1)) +&
                      f2*kpoint(3)*s(coset(ax,ay,az),coset(0,0,lb-1))*gaussi +&
                      fz*s(coset(ax,ay,az-1),coset(0,0,lb-1)) +&
                      f3*s(coset(ax,ay,az),coset(0,0,lb-2))
                  END IF

!                 Increase the angular momentum component y of function b

                  IF (ay == 0) THEN
                    bz = lb - 1
                    s(coset(ax,ay,az),coset(0,1,bz)) =&
                      rbp(2)*s(coset(ax,ay,az),coset(0,0,bz)) +&
                      f2*kpoint(2)*s(coset(ax,ay,az),coset(0,0,bz))*gaussi
                    DO by=2,lb
                      bz = lb - by
                      f3 = f2*REAL(by-1,dp)
                      s(coset(ax,ay,az),coset(0,by,bz)) =&
                        rbp(2)*s(coset(ax,ay,az),coset(0,by-1,bz)) +&
                        f2*kpoint(2)*s(coset(ax,ay,az),coset(0,by-1,bz))*gaussi +&
                        f3*s(coset(ax,ay,az),coset(0,by-2,bz))
                    END DO
                  ELSE
                    bz = lb - 1
                    s(coset(ax,ay,az),coset(0,1,bz)) =&
                      rbp(2)*s(coset(ax,ay,az),coset(0,0,bz)) +&
                      f2*kpoint(2)*s(coset(ax,ay,az),coset(0,0,bz))*gaussi +&
                      fy*s(coset(ax,ay-1,az),coset(0,0,bz))
                    DO by=2,lb
                      bz = lb - by
                      f3 = f2*REAL(by-1,dp)
                      s(coset(ax,ay,az),coset(0,by,bz)) =&
                        rbp(2)*s(coset(ax,ay,az),coset(0,by-1,bz)) +&
                        f2*kpoint(2)*s(coset(ax,ay,az),coset(0,by-1,bz))*gaussi +&
                        fy*s(coset(ax,ay-1,az),coset(0,by-1,bz)) +&
                        f3*s(coset(ax,ay,az),coset(0,by-2,bz))
                    END DO
                  END IF

!                 Increase the angular momentum component x of function b

                  IF (ax == 0) THEN
                    DO by=0,lb-1
                      bz = lb - 1 - by
                      s(coset(ax,ay,az),coset(1,by,bz)) =&
                        rbp(1)*s(coset(ax,ay,az),coset(0,by,bz)) +&
                        f2*kpoint(1)*s(coset(ax,ay,az),coset(0,by,bz))*gaussi
                    END DO
                    DO bx=2,lb
                      f3 = f2*REAL(bx-1,dp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        s(coset(ax,ay,az),coset(bx,by,bz)) =&
                          rbp(1)*s(coset(ax,ay,az),coset(bx-1,by,bz)) +&
                          f2*kpoint(1)*s(coset(ax,ay,az),coset(bx-1,by,bz))*gaussi +&
                          f3*s(coset(ax,ay,az),coset(bx-2,by,bz))
                      END DO
                    END DO
                  ELSE
                    DO by=0,lb-1
                      bz = lb - 1 - by
                      s(coset(ax,ay,az),coset(1,by,bz)) =&
                        rbp(1)*s(coset(ax,ay,az),coset(0,by,bz)) +&
                        f2*kpoint(1)*s(coset(ax,ay,az),coset(0,by,bz))*gaussi +&
                        fx*s(coset(ax-1,ay,az),coset(0,by,bz))
                    END DO
                    DO bx=2,lb
                      f3 = f2*REAL(bx-1,dp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        s(coset(ax,ay,az),coset(bx,by,bz)) =&
                          rbp(1)*s(coset(ax,ay,az),coset(bx-1,by,bz)) +&
                          f2*kpoint(1)*s(coset(ax,ay,az),coset(bx-1,by,bz))*gaussi +&
                          fx*s(coset(ax-1,ay,az),coset(bx-1,by,bz)) +&
                          f3*s(coset(ax,ay,az),coset(bx-2,by,bz))
                      END DO
                    END DO
                  END IF

                END DO
              END DO

            END DO

          END IF

        ELSE

          IF (lb_max > 0) THEN

!           Vertical recurrence steps: [s|k|s] -> [s|k|b]

            rbp(:) = (f1 - 1.0_dp)*rab(:)

!           [s|k|p] = (Pi - Bi)*[s|k|s] + f2*ki*I*[s|k|s]

            s(1,2) = rbp(1)*s(1,1) + f2*kpoint(1)*s(1,1)*gaussi
            s(1,3) = rbp(2)*s(1,1) + f2*kpoint(2)*s(1,1)*gaussi
            s(1,4) = rbp(3)*s(1,1) + f2*kpoint(3)*s(1,1)*gaussi

!           [s|k|b] = (Pi - Bi)*[s|k|b-1i] + f2*ki*I*[s|k|b-1i] + 
!                     f2*Ni(b-1i)*[s|k|b-2i]

            DO lb=2,lb_max

!             Increase the angular momentum component z of function b

              s(1,coset(0,0,lb)) = rbp(3)*s(1,coset(0,0,lb-1)) +&
                                   f2*kpoint(3)*s(1,coset(0,0,lb-1))*gaussi +&
                                   f2*REAL(lb-1,dp)*s(1,coset(0,0,lb-2))

!             Increase the angular momentum component y of function b

              bz = lb - 1
              s(1,coset(0,1,bz)) = rbp(2)*s(1,coset(0,0,bz)) +&
                                   f2*kpoint(2)*s(1,coset(0,0,bz))*gaussi

              DO by=2,lb
                bz = lb - by
                s(1,coset(0,by,bz)) = rbp(2)*s(1,coset(0,by-1,bz)) +&
                                      f2*kpoint(2)*s(1,coset(0,by-1,bz))*gaussi +&
                                      f2*REAL(by-1,dp)*s(1,coset(0,by-2,bz))
              END DO

!             Increase the angular momentum component x of function b

              DO by=0,lb-1
                bz = lb - 1 - by
                s(1,coset(1,by,bz)) = rbp(1)*s(1,coset(0,by,bz)) +&
                                      f2*kpoint(1)*s(1,coset(0,by,bz))*gaussi
              END DO

              DO bx=2,lb
                f3 = f2*REAL(bx-1,dp)
                DO by=0,lb-bx
                  bz = lb - bx - by
                  s(1,coset(bx,by,bz)) = rbp(1)*s(1,coset(bx-1,by,bz)) +&
                                         f2*kpoint(1)*s(1,coset(bx-1,by,bz))*gaussi +&
                                         f3*s(1,coset(bx-2,by,bz))
                END DO
              END DO

            END DO

          END IF

        END IF

        DO j=1,ncoset(lb_max)
          DO i=1,ncoset(la_max)
            sabk(na+i,nb+j) = s(i,j)
          END DO
        END DO

        nb = nb + ncoset(lb_max)

      END DO

      na = na + ncoset(la_max)

    END DO

  END SUBROUTINE overlap_k

! *****************************************************************************

END MODULE ai_overlap_k

! *****************************************************************************

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/ai_overlap3 [1.0] *
!!
!!   NAME
!!     ai_overlap
!!
!!   FUNCTION
!!     Calculation of three-center overlap integrals over Cartesian
!!     Gaussian-type functions.
!!
!!   AUTHOR
!!     Matthias Krack (26.06.2001)
!!
!!   LITERATURE
!!     S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!!
!!   MODIFICATION HISTORY
!!     - Derivatives added (02.05.2002,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE ai_overlap3

! *****************************************************************************

! ax,ay,az   : Angular momentum index numbers of orbital a.
! bx,by,bz   : Angular momentum index numbers of orbital b.
! coset      : Cartesian orbital set pointer.
! dab        : Distance between the atomic centers a and b.
! dac        : Distance between the atomic centers a and c.
! dbc        : Distance between the atomic centers b and c.
! l{a,b,c}   : Angular momentum quantum number of shell a, b or c.
! l{a,b}_max : Maximum angular momentum quantum number of shell a, b or c.
! ncoset     : Number of Cartesian orbitals up to l.
! rab        : Distance vector between the atomic centers a and b.
! rac        : Distance vector between the atomic centers a and c.
! rbc        : Distance vector between the atomic centers b and c.
! rpgf{a,b,c}: Radius of the primitive Gaussian-type function a or b.
! zet{a,b,c} : Exponents of the Gaussian-type functions a or b.
! zetg       : Reciprocal of the sum of the exponents of orbital a, b and c.
! zetp       : Reciprocal of the sum of the exponents of orbital a and b.

! *****************************************************************************

  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             indco,&
                                             ncoset
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: overlap3, overlap3_old

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE overlap3_old(la_max,npgfa,zeta,rpgfa,la_min,&
                      lb_max,npgfb,zetb,rpgfb,lb_min,&
                      lc_max,npgfc,zetc,rpgfc,lc_min,&
                      rab,dab,rac,dac,rbc,dbc,sabc)

!   Purpose: Calculation of three-center overlap integrals <a|c|b> over
!            Cartesian Gaussian functions. The Gaussian-type functions a and b
!            may be contracted and c is a primitive Gaussian-type function.

!   History: - Creation (26.06.2001, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: la_min, lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb, rpgfb
    INTEGER, INTENT(IN)                      :: lb_min, lc_max, npgfc
    REAL(KIND=dp), INTENT(IN)                :: zetc, rpgfc
    INTEGER, INTENT(IN)                      :: lc_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: dab
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac
    REAL(KIND=dp), INTENT(IN)                :: dac
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rbc
    REAL(KIND=dp), INTENT(IN)                :: dbc
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: sabc

    INTEGER :: ax, ay, az, bx, by, bz, coa, coax, coay, coaz, coc, cocx, &
      cocy, cocz, cx, cy, cz, i, ipgf, j, jpgf, k, la, la_start, lb, lc, na, &
      nb
    REAL(KIND=dp)                            :: f0, f1, f2, f3, fcx, fcy, &
                                                fcz, fx, fy, fz, rcp2, zetg, &
                                                zetp
    REAL(KIND=dp), DIMENSION(3)              :: rag, rbg, rcg, rcp
    REAL(KIND=dp), DIMENSION(ncoset(la_max),&
      ncoset(lb_max), ncoset(lc_max))        :: s

!   ---------------------------------------------------------------------------
!   *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0

    DO ipgf=1,npgfa

!     *** Screening ***

      IF (rpgfa(ipgf) + rpgfc < dac) THEN
        na = na + ncoset(la_max)
        CYCLE
      END IF

      nb = 0

      DO jpgf=1,npgfb

!       *** Screening ***

        IF ((rpgfb(jpgf) + rpgfc < dbc).OR.&
            (rpgfa(ipgf) + rpgfb(jpgf) < dab)) THEN
          nb = nb + ncoset(lb_max)
          CYCLE
        END IF

!       *** Calculate some prefactors ***

        zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))
        zetg = 1.0_dp/(zeta(ipgf) + zetb(jpgf) + zetc)

        f0 = (pi*zetg)**1.5_dp
        f1 = zetb(jpgf)*zetp
        f2 = 0.5_dp*zetg

        rcp(:) = f1*rab(:) - rac(:)

        rcp2 = rcp(1)*rcp(1) + rcp(2)*rcp(2) + rcp(3)*rcp(3)

!       *** Calculate the basic three-center overlap integral [s|s|s] ***

        s(1,1,1) = f0*EXP(-(zeta(ipgf)*f1*dab*dab + zetc*zetg*rcp2/zetp))

!       *** Recurrence steps: [s|s|s] -> [a|s|s] ***

        IF (la_max > 0) THEN

!         *** Vertical recurrence steps: [s|s|s] -> [a|s|s] ***

          rag(:) = zetg*(zetb(jpgf)*rab(:) + zetc*rac(:))

!         *** [p|s|s] = (Gi - Ai)*[s|s|s]  (i = x,y,z) ***

          s(2,1,1) = rag(1)*s(1,1,1)
          s(3,1,1) = rag(2)*s(1,1,1)
          s(4,1,1) = rag(3)*s(1,1,1)

!         *** [a|s|s] = (Gi - Ai)*[a-1i|s|s] + f2*Ni(a-1i)*[a-2i|s|s] ***

          DO la=2,la_max

!           *** Increase the angular momentum component z of function a ***

            s(coset(0,0,la),1,1) = rag(3)*s(coset(0,0,la-1),1,1) +&
                                   f2*REAL(la-1,dp)*s(coset(0,0,la-2),1,1)

!           *** Increase the angular momentum component y of function a ***

            az = la - 1
            s(coset(0,1,az),1,1) = rag(2)*s(coset(0,0,az),1,1)

            DO ay=2,la
              az = la - ay
              s(coset(0,ay,az),1,1) = rag(2)*s(coset(0,ay-1,az),1,1) +&
                                      f2*REAL(ay-1,dp)*s(coset(0,ay-2,az),1,1)
            END DO

!           *** Increase the angular momentum component x of function a ***

            DO ay=0,la-1
              az = la - 1 - ay
              s(coset(1,ay,az),1,1) = rag(1)*s(coset(0,ay,az),1,1)
            END DO

            DO ax=2,la
              f3 = f2*REAL(ax-1,dp)
              DO ay=0,la-ax
                az = la - ax - ay
                s(coset(ax,ay,az),1,1) = rag(1)*s(coset(ax-1,ay,az),1,1) +&
                                         f3*s(coset(ax-2,ay,az),1,1)
              END DO
            END DO

          END DO

!         *** Recurrence steps: [a|s|s] -> [a|s|b] ***

          IF (lb_max > 0) THEN

!           *** Horizontal recurrence steps ***

            rbg(:) = rag(:) - rab(:)

!           *** [a|s|p] = [a+1i|s|s] - (Bi - Ai)*[a|s|s] ***

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
                  coax = coset(ax+1,ay,az)
                  coay = coset(ax,ay+1,az)
                  coaz = coset(ax,ay,az+1)
                  s(coset(ax,ay,az),2,1) = s(coax,1,1) - rab(1)*s(coa,1,1)
                  s(coset(ax,ay,az),3,1) = s(coay,1,1) - rab(2)*s(coa,1,1)
                  s(coset(ax,ay,az),4,1) = s(coaz,1,1) - rab(3)*s(coa,1,1)
                END DO
              END DO
            END DO

!           *** Vertical recurrence step ***

!           *** [a|s|p] = (Gi - Bi)*[a|s|s] + f2*Ni(a)*[a-1i|s|s] ***

            DO ax=0,la_max
              fx = f2*REAL(ax,dp)
              DO ay=0,la_max-ax
                fy = f2*REAL(ay,dp)
                az = la_max - ax - ay
                fz = f2*REAL(az,dp)
                coa = coset(ax,ay,az)
                IF (ax == 0) THEN
                  s(coa,2,1) = rbg(1)*s(coa,1,1)
                ELSE
                  s(coa,2,1) = rbg(1)*s(coa,1,1) + fx*s(coset(ax-1,ay,az),1,1)
                END IF
                IF (ay == 0) THEN
                  s(coa,3,1) = rbg(2)*s(coa,1,1)
                ELSE
                  s(coa,3,1) = rbg(2)*s(coa,1,1) + fy*s(coset(ax,ay-1,az),1,1)
                END IF
                IF (az == 0) THEN
                  s(coa,4,1) = rbg(3)*s(coa,1,1)
                ELSE
                  s(coa,4,1) = rbg(3)*s(coa,1,1) + fz*s(coset(ax,ay,az-1),1,1)
                END IF
              END DO
            END DO

!           *** Recurrence steps: [a|s|p] -> [a|s|b] ***

            DO lb=2,lb_max

!             *** Horizontal recurrence steps ***

!             *** [a|s|b] = [a+1i|s|b-1i] - (Bi - Ai)*[a|s|b-1i] ***

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
                    coax = coset(ax+1,ay,az)
                    coay = coset(ax,ay+1,az)
                    coaz = coset(ax,ay,az+1)

!                   *** Shift of angular momentum component z from a to b ***

                    s(coa,coset(0,0,lb),1) =&
                      s(coaz,coset(0,0,lb-1),1) -&
                      rab(3)*s(coa,coset(0,0,lb-1),1)

!                   *** Shift of angular momentum component y from a to b ***

                    DO by=1,lb
                      bz = lb - by
                      s(coa,coset(0,by,bz),1) =&
                        s(coay,coset(0,by-1,bz),1) -&
                        rab(2)*s(coa,coset(0,by-1,bz),1)
                    END DO

!                   *** Shift of angular momentum component x from a to b ***

                    DO bx=1,lb
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        s(coa,coset(bx,by,bz),1) =&
                          s(coax,coset(bx-1,by,bz),1) -&
                          rab(1)*s(coa,coset(bx-1,by,bz),1)
                      END DO
                    END DO

                  END DO
                END DO
              END DO

!             *** Vertical recurrence step ***

!             *** [a|s|b] = (Gi - Bi)*[a|s|b-1i] +   ***
!             ***           f2*Ni(a)*[a-1i|s|b-1i] + ***
!             ***           f2*Ni(b-1i)*[a|s|b-2i]   ***

              DO ax=0,la_max
                fx = f2*REAL(ax,dp)
                DO ay=0,la_max-ax
                  fy = f2*REAL(ay,dp)
                  az = la_max - ax - ay
                  fz = f2*REAL(az,dp)

                  coa = coset(ax,ay,az)

                  f3 = f2*REAL(lb-1,dp)

!                 *** Shift of angular momentum component z from a to b ***

                  IF (az == 0) THEN
                    s(coa,coset(0,0,lb),1) =&
                      rbg(3)*s(coa,coset(0,0,lb-1),1) +&
                      f3*s(coa,coset(0,0,lb-2),1)
                  ELSE
                    coaz = coset(ax,ay,az-1)
                    s(coa,coset(0,0,lb),1) =&
                      rbg(3)*s(coa,coset(0,0,lb-1),1) +&
                      fz*s(coaz,coset(0,0,lb-1),1) +&
                      f3*s(coa,coset(0,0,lb-2),1)
                  END IF

!                 *** Shift of angular momentum component y from a to b ***

                  IF (ay == 0) THEN
                    bz = lb - 1
                    s(coa,coset(0,1,bz),1) =&
                      rbg(2)*s(coa,coset(0,0,bz),1)
                    DO by=2,lb
                      bz = lb - by
                      f3 = f2*REAL(by-1,dp)
                      s(coa,coset(0,by,bz),1) =&
                        rbg(2)*s(coa,coset(0,by-1,bz),1) +&
                        f3*s(coa,coset(0,by-2,bz),1)
                    END DO
                  ELSE
                    coay = coset(ax,ay-1,az)
                    bz = lb - 1
                    s(coa,coset(0,1,bz),1) =&
                      rbg(2)*s(coa,coset(0,0,bz),1) +&
                      fy*s(coay,coset(0,0,bz),1)
                    DO by=2,lb
                      bz = lb - by
                      f3 = f2*REAL(by-1,dp)
                      s(coa,coset(0,by,bz),1) =&
                        rbg(2)*s(coa,coset(0,by-1,bz),1) +&
                        fy*s(coay,coset(0,by-1,bz),1) +&
                        f3*s(coa,coset(0,by-2,bz),1)
                    END DO
                  END IF

!                 *** Shift of angular momentum component x from a to b ***

                  IF (ax == 0) THEN
                    DO by=0,lb-1
                      bz = lb - 1 - by
                      s(coa,coset(1,by,bz),1) =&
                        rbg(1)*s(coa,coset(0,by,bz),1)
                    END DO
                    DO bx=2,lb
                      f3 = f2*REAL(bx-1,dp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        s(coa,coset(bx,by,bz),1) =&
                          rbg(1)*s(coa,coset(bx-1,by,bz),1) +&
                          f3*s(coa,coset(bx-2,by,bz),1)
                      END DO
                    END DO
                  ELSE
                    coax = coset(ax-1,ay,az)
                    DO by=0,lb-1
                      bz = lb - 1 - by
                      s(coa,coset(1,by,bz),1) =&
                        rbg(1)*s(coa,coset(0,by,bz),1) +&
                        fx*s(coax,coset(0,by,bz),1)
                    END DO
                    DO bx=2,lb
                      f3 = f2*REAL(bx-1,dp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        s(coa,coset(bx,by,bz),1) =&
                          rbg(1)*s(coa,coset(bx-1,by,bz),1) +&
                          fx*s(coax,coset(bx-1,by,bz),1) +&
                          f3*s(coa,coset(bx-2,by,bz),1)
                      END DO
                    END DO
                  END IF

                END DO
              END DO

            END DO

          END IF

        ELSE

          IF (lb_max > 0) THEN

!           *** Vertical recurrence steps: [s|s|s] -> [s|s|b] ***

            rbg(:) = -zetg*(zeta(ipgf)*rab(:) - zetc*rbc(:))

!           *** [s|s|p] = (Gi - Bi)*[s|s|s] ***

            s(1,2,1) = rbg(1)*s(1,1,1)
            s(1,3,1) = rbg(2)*s(1,1,1)
            s(1,4,1) = rbg(3)*s(1,1,1)

!           *** [s|s|b] = (Gi - Bi)*[s|s|b-1i] + f2*Ni(b-1i)*[s|s|b-2i] ***

            DO lb=2,lb_max

!             *** Increase the angular momentum component z of function b ***

              s(1,coset(0,0,lb),1) = rbg(3)*s(1,coset(0,0,lb-1),1) +&
                                     f2*REAL(lb-1,dp)*s(1,coset(0,0,lb-2),1)

!             *** Increase the angular momentum component y of function b ***

              bz = lb - 1
              s(1,coset(0,1,bz),1) = rbg(2)*s(1,coset(0,0,bz),1)

              DO by=2,lb
                bz = lb - by
                s(1,coset(0,by,bz),1) =&
                  rbg(2)*s(1,coset(0,by-1,bz),1) +&
                  f2*REAL(by-1,dp)*s(1,coset(0,by-2,bz),1)
              END DO

!             *** Increase the angular momentum component x of function b ***

              DO by=0,lb-1
                bz = lb - 1 - by
                s(1,coset(1,by,bz),1) = rbg(1)*s(1,coset(0,by,bz),1)
              END DO

              DO bx=2,lb
                f3 = f2*REAL(bx-1,dp)
                DO by=0,lb-bx
                  bz = lb - bx - by
                  s(1,coset(bx,by,bz),1) = rbg(1)*s(1,coset(bx-1,by,bz),1) +&
                                           f3*s(1,coset(bx-2,by,bz),1)
                END DO
              END DO

            END DO

          END IF

        END IF

!       *** Recurrence steps: [a|s|b] -> [a|c|b] ***

        IF (lc_max > 0) THEN

!         *** Vertical recurrence steps: [s|s|s] -> [s|c|s] ***

          rcg(:) = -zetg*(zeta(ipgf)*rac(:) + zetb(jpgf)*rbc(:))

!         *** [s|p|s] = (Gi - Ci)*[s|s|s]  (i = x,y,z) ***

          s(1,1,2) = rcg(1)*s(1,1,1)
          s(1,1,3) = rcg(2)*s(1,1,1)
          s(1,1,4) = rcg(3)*s(1,1,1)

!         *** [s|c|s] = (Gi - Ci)*[s|c-1i|s] + f2*Ni(c-1i)*[s|c-2i|s] ***

          DO lc=2,lc_max

!           *** Increase the angular momentum component z of function c ***

            s(1,1,coset(0,0,lc)) = rcg(3)*s(1,1,coset(0,0,lc-1)) +&
                                   f2*REAL(lc-1,dp)*s(1,1,coset(0,0,lc-2))

!           *** Increase the angular momentum component y of function c ***

            cz = lc - 1
            s(1,1,coset(0,1,cz)) = rcg(2)*s(1,1,coset(0,0,cz))

            DO cy=2,lc
              cz = lc - cy
              s(1,1,coset(0,cy,cz)) = rcg(2)*s(1,1,coset(0,cy-1,cz)) +&
                                      f2*REAL(cy-1,dp)*s(1,1,coset(0,cy-2,cz))
            END DO

!           *** Increase the angular momentum component x of function c ***

            DO cy=0,lc-1
              cz = lc - 1 - cy
              s(1,1,coset(1,cy,cz)) = rcg(1)*s(1,1,coset(0,cy,cz))
            END DO

            DO cx=2,lc
              f3 = f2*REAL(cx-1,dp)
              DO cy=0,lc-cx
                cz = lc - cx - cy
                s(1,1,coset(cx,cy,cz)) = rcg(1)*s(1,1,coset(cx-1,cy,cz)) +&
                                         f3*s(1,1,coset(cx-2,cy,cz))
              END DO
            END DO

          END DO

!         *** Recurrence steps: [s|c|s] -> [a|c|b] ***

          DO lc=1,lc_max

            DO cx=0,lc
              DO cy=0,lc-cx
                cz = lc - cx - cy

                coc = coset(cx,cy,cz)
                cocx = coset(MAX(0,cx-1),cy,cz)
                cocy = coset(cx,MAX(0,cy-1),cz)
                cocz = coset(cx,cy,MAX(0,cz-1))

                fcx = f2*REAL(cx,dp)
                fcy = f2*REAL(cy,dp)
                fcz = f2*REAL(cz,dp)

!               *** Recurrence steps: [s|c|s] -> [a|c|s] ***

                IF (la_max > 0) THEN

!                 *** Vertical recurrence steps: [s|c|s] -> [a|c|s] ***

                  rag(:) = rcg(:) + rac(:)

!                 *** [p|c|s] = (Gi - Ai)*[s|c|s] + f2*Ni(c)*[s|c-1i|s] ***

                  s(2,1,coc) = rag(1)*s(1,1,coc) + fcx*s(1,1,cocx)
                  s(3,1,coc) = rag(2)*s(1,1,coc) + fcy*s(1,1,cocy)
                  s(4,1,coc) = rag(3)*s(1,1,coc) + fcz*s(1,1,cocz)

!                 *** [a|c|s] = (Gi - Ai)*[a-1i|c|s] +   ***
!                 ***           f2*Ni(a-1i)*[a-2i|c|s] + ***
!                 ***           f2*Ni(c)*[a-1i|c-1i|s]   ***

                  DO la=2,la_max

!                   *** Increase the angular momentum component z of a ***

                    s(coset(0,0,la),1,coc) =&
                      rag(3)*s(coset(0,0,la-1),1,coc) +&
                      f2*REAL(la-1,dp)*s(coset(0,0,la-2),1,coc) +&
                      fcz*s(coset(0,0,la-1),1,cocz)

!                   *** Increase the angular momentum component y of a ***

                    az = la - 1
                    s(coset(0,1,az),1,coc) =&
                      rag(2)*s(coset(0,0,az),1,coc) +&
                      fcy*s(coset(0,0,az),1,cocy)

                    DO ay=2,la
                      az = la - ay
                      s(coset(0,ay,az),1,coc) =&
                        rag(2)*s(coset(0,ay-1,az),1,coc) +&
                        f2*REAL(ay-1,dp)*s(coset(0,ay-2,az),1,coc) +&
                        fcy*s(coset(0,ay-1,az),1,cocy)
                    END DO

!                   *** Increase the angular momentum component x of a ***

                    DO ay=0,la-1
                      az = la - 1 - ay
                      s(coset(1,ay,az),1,coc) =&
                        rag(1)*s(coset(0,ay,az),1,coc) +&
                        fcx*s(coset(0,ay,az),1,cocx)
                    END DO

                    DO ax=2,la
                      f3 = f2*REAL(ax-1,dp)
                      DO ay=0,la-ax
                        az = la - ax - ay
                        s(coset(ax,ay,az),1,coc) =&
                          rag(1)*s(coset(ax-1,ay,az),1,coc) +&
                          f3*s(coset(ax-2,ay,az),1,coc) +&
                          fcx*s(coset(ax-1,ay,az),1,cocx)
                      END DO
                    END DO

                  END DO

!                 *** Recurrence steps: [a|c|s] -> [a|c|b] ***

                  IF (lb_max > 0) THEN

!                   *** Horizontal recurrence steps ***

                    rbg(:) = rag(:) - rab(:)

!                   *** [a|c|p] = [a+1i|c|s] - (Bi - Ai)*[a|c|s] ***

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
                          coax = coset(ax+1,ay,az)
                          coay = coset(ax,ay+1,az)
                          coaz = coset(ax,ay,az+1)
                          s(coa,2,coc) = s(coax,1,coc) - rab(1)*s(coa,1,coc)
                          s(coa,3,coc) = s(coay,1,coc) - rab(2)*s(coa,1,coc)
                          s(coa,4,coc) = s(coaz,1,coc) - rab(3)*s(coa,1,coc)
                        END DO
                      END DO
                    END DO

!                   *** Vertical recurrence step ***

!                   *** [a|c|p] = (Gi - Bi)*[a|c|s] +   ***
!                                 f2*Ni(a)*[a-1i|c|s] + ***
!                                 f2*Ni(c)*[a|c-1i|s]   ***

                    DO ax=0,la_max
                      fx = f2*REAL(ax,dp)
                      DO ay=0,la_max-ax
                        fy = f2*REAL(ay,dp)
                        az = la_max - ax - ay
                        fz = f2*REAL(az,dp)
                        coa = coset(ax,ay,az)
                        IF (ax == 0) THEN
                          s(coa,2,coc) = rbg(1)*s(coa,1,coc) +&
                                         fcx*s(coa,1,cocx)
                        ELSE
                          s(coa,2,coc) = rbg(1)*s(coa,1,coc) +&
                                         fx*s(coset(ax-1,ay,az),1,coc) +&
                                         fcx*s(coa,1,cocx)
                        END IF
                        IF (ay == 0) THEN
                          s(coa,3,coc) = rbg(2)*s(coa,1,coc) +&
                                         fcy*s(coa,1,cocy)
                        ELSE
                          s(coa,3,coc) = rbg(2)*s(coa,1,coc) +&
                                         fy*s(coset(ax,ay-1,az),1,coc) +&
                                         fcy*s(coa,1,cocy)
                        END IF
                        IF (az == 0) THEN
                          s(coa,4,coc) = rbg(3)*s(coa,1,coc) +&
                                         fcz*s(coa,1,cocz)
                        ELSE
                          s(coa,4,coc) = rbg(3)*s(coa,1,coc) +&
                                         fz*s(coset(ax,ay,az-1),1,coc) +&
                                         fcz*s(coa,1,cocz)
                        END IF
                      END DO
                    END DO

!                   *** Recurrence steps: [a|c|p] -> [a|c|b] ***

                    DO lb=2,lb_max

!                     *** Horizontal recurrence steps ***

!                     *** [a|c|b] = [a+1i|c|b-1i] - (Bi - Ai)*[a|c|b-1i] ***

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
                            coax = coset(ax+1,ay,az)
                            coay = coset(ax,ay+1,az)
                            coaz = coset(ax,ay,az+1)

!                           *** Shift of angular momentum ***
!                           *** component z from a to b   ***

                            s(coa,coset(0,0,lb),coc) =&
                              s(coaz,coset(0,0,lb-1),coc) -&
                              rab(3)*s(coa,coset(0,0,lb-1),coc)

!                           *** Shift of angular momentum ***
!                           *** component y from a to b   ***

                            DO by=1,lb
                              bz = lb - by
                              s(coa,coset(0,by,bz),coc) =&
                                s(coay,coset(0,by-1,bz),coc) -&
                                rab(2)*s(coa,coset(0,by-1,bz),coc)
                            END DO

!                           *** Shift of angular momentum ***
!                           *** component x from a to b   ***

                            DO bx=1,lb
                              DO by=0,lb-bx
                                bz = lb - bx - by
                                s(coa,coset(bx,by,bz),coc) =&
                                  s(coax,coset(bx-1,by,bz),coc) -&
                                  rab(1)*s(coa,coset(bx-1,by,bz),coc)
                              END DO
                            END DO

                          END DO
                        END DO
                      END DO

!                     *** Vertical recurrence step ***

!                     *** [a|c|b] = (Gi - Bi)*[a|c|b-1i] +   ***
!                     ***           f2*Ni(a)*[a-1i|c|b-1i] + ***
!                     ***           f2*Ni(b-1i)*[a|c|b-2i] + ***
!                     ***           f2*Ni(c)*[a|c-1i|b-1i]   ***

                      DO ax=0,la_max
                        fx = f2*REAL(ax,dp)
                        DO ay=0,la_max-ax
                          fy = f2*REAL(ay,dp)
                          az = la_max - ax - ay
                          fz = f2*REAL(az,dp)

                          coa = coset(ax,ay,az)
                          coax = coset(MAX(0,ax-1),ay,az)
                          coay = coset(ax,MAX(0,ay-1),az)
                          coaz = coset(ax,ay,MAX(0,az-1))

                          f3 = f2*REAL(lb-1,dp)

!                         *** Shift of angular momentum ***
!                         *** component z from a to b   ***

                          IF (az == 0) THEN
                            s(coa,coset(0,0,lb),coc) =&
                              rbg(3)*s(coa,coset(0,0,lb-1),coc) +&
                              f3*s(coa,coset(0,0,lb-2),coc) +&
                              fcz*s(coa,coset(0,0,lb-1),cocz)
                          ELSE
                            s(coa,coset(0,0,lb),coc) =&
                              rbg(3)*s(coa,coset(0,0,lb-1),coc) +&
                              fz*s(coaz,coset(0,0,lb-1),coc) +&
                              f3*s(coa,coset(0,0,lb-2),coc) +&
                              fcz*s(coa,coset(0,0,lb-1),cocz)
                          END IF

!                         *** Shift of angular momentum ***
!                         *** component y from a to b   ***

                          IF (ay == 0) THEN
                            bz = lb - 1
                            s(coa,coset(0,1,bz),coc) =&
                              rbg(2)*s(coa,coset(0,0,bz),coc) +&
                              fcy*s(coa,coset(0,0,bz),cocy)
                            DO by=2,lb
                              bz = lb - by
                              f3 = f2*REAL(by-1,dp)
                              s(coa,coset(0,by,bz),coc) =&
                                rbg(2)*s(coa,coset(0,by-1,bz),coc) +&
                                f3*s(coa,coset(0,by-2,bz),coc) +&
                                fcy*s(coa,coset(0,by-1,bz),cocy)
                            END DO
                          ELSE
                            bz = lb - 1
                            s(coa,coset(0,1,bz),coc) =&
                              rbg(2)*s(coa,coset(0,0,bz),coc) +&
                              fy*s(coay,coset(0,0,bz),coc) +&
                              fcy*s(coa,coset(0,0,bz),cocy)
                            DO by=2,lb
                              bz = lb - by
                              f3 = f2*REAL(by-1,dp)
                              s(coa,coset(0,by,bz),coc) =&
                                rbg(2)*s(coa,coset(0,by-1,bz),coc) +&
                                fy*s(coay,coset(0,by-1,bz),coc) +&
                                f3*s(coa,coset(0,by-2,bz),coc) +&
                                fcy*s(coa,coset(0,by-1,bz),cocy)
                            END DO
                          END IF

!                         *** Shift of angular momentum ***
!                         *** component x from a to b   ***

                          IF (ax == 0) THEN
                            DO by=0,lb-1
                              bz = lb - 1 - by
                              s(coa,coset(1,by,bz),coc) =&
                                rbg(1)*s(coa,coset(0,by,bz),coc) +&
                                fcx*s(coa,coset(0,by,bz),cocx)
                            END DO
                            DO bx=2,lb
                              f3 = f2*REAL(bx-1,dp)
                              DO by=0,lb-bx
                                bz = lb - bx - by
                                s(coa,coset(bx,by,bz),coc) =&
                                  rbg(1)*s(coa,coset(bx-1,by,bz),coc) +&
                                  f3*s(coa,coset(bx-2,by,bz),coc) +&
                                  fcx*s(coa,coset(bx-1,by,bz),cocx)
                              END DO
                            END DO
                          ELSE
                            DO by=0,lb-1
                              bz = lb - 1 - by
                              s(coa,coset(1,by,bz),coc) =&
                                rbg(1)*s(coa,coset(0,by,bz),coc) +&
                                fx*s(coax,coset(0,by,bz),coc) +&
                                fcx*s(coa,coset(0,by,bz),cocx)
                            END DO
                            DO bx=2,lb
                              f3 = f2*REAL(bx-1,dp)
                              DO by=0,lb-bx
                                bz = lb - bx - by
                                s(coa,coset(bx,by,bz),coc) =&
                                  rbg(1)*s(coa,coset(bx-1,by,bz),coc) +&
                                  fx*s(coax,coset(bx-1,by,bz),coc) +&
                                  f3*s(coa,coset(bx-2,by,bz),coc) +&
                                  fcx*s(coa,coset(bx-1,by,bz),cocx)
                              END DO
                            END DO
                          END IF

                        END DO
                      END DO

                    END DO

                  END IF

                ELSE

                  IF (lb_max > 0) THEN

!                   *** Vertical recurrence steps: [s|c|s] -> [s|c|b] ***

                    rbg(:) = rcg(:) + rbc(:)

!                   *** [s|c|p] = (Gi - Bi)*[s|c|s] + f2*Ni(c)*[s|c-1i|s] ***

                    s(1,2,coc) = rbg(1)*s(1,1,coc) + fcx*s(1,1,cocx)
                    s(1,3,coc) = rbg(2)*s(1,1,coc) + fcy*s(1,1,cocy)
                    s(1,4,coc) = rbg(3)*s(1,1,coc) + fcz*s(1,1,cocz)

!                   *** [s|c|b] = (Gi - Bi)*[s|c|b-1i] + ***
!                   ***           f2*Ni(b-1i)*[s|c|b-2i] ***
!                   ***           f2*Ni(c)*[s|c-1i|b-1i] ***

                    DO lb=2,lb_max

!                     *** Increase the angular momentum component z of b ***

                      s(1,coset(0,0,lb),coc) =&
                        rbg(3)*s(1,coset(0,0,lb-1),coc) +&
                        f2*REAL(lb-1,dp)*s(1,coset(0,0,lb-2),coc) +&
                        fcz*s(1,coset(0,0,lb-1),cocz)

!                     *** Increase the angular momentum component y of b ***

                      bz = lb - 1
                      s(1,coset(0,1,bz),coc) =&
                        rbg(2)*s(1,coset(0,0,bz),coc) +&
                        fcy*s(1,coset(0,0,bz),cocy)

                      DO by=2,lb
                        bz = lb - by
                        s(1,coset(0,by,bz),coc) =&
                          rbg(2)*s(1,coset(0,by-1,bz),coc) +&
                          f2*REAL(by-1,dp)*s(1,coset(0,by-2,bz),coc) +&
                          fcy*s(1,coset(0,by-1,bz),cocy)
                      END DO

!                     *** Increase the angular momentum component x of b ***

                      DO by=0,lb-1
                        bz = lb - 1 - by
                        s(1,coset(1,by,bz),coc) =&
                          rbg(1)*s(1,coset(0,by,bz),coc) +&
                          fcx*s(1,coset(0,by,bz),cocx)
                      END DO

                      DO bx=2,lb
                        f3 = f2*REAL(bx-1,dp)
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          s(1,coset(bx,by,bz),coc) =&
                            rbg(1)*s(1,coset(bx-1,by,bz),coc) +&
                            f3*s(1,coset(bx-2,by,bz),coc) +&
                            fcx*s(1,coset(bx-1,by,bz),cocx)
                        END DO
                      END DO

                    END DO

                  END IF

                END IF

              END DO
            END DO

          END DO

        END IF

!       *** Add the contribution of the current pair ***
!       *** of primitive Gaussian-type functions     ***

        DO k=ncoset(lc_min-1)+1,ncoset(lc_max)
          DO j=ncoset(lb_min-1)+1,ncoset(lb_max)
            DO i=ncoset(la_min-1)+1,ncoset(la_max)
              sabc(na+i,nb+j,k) = s(i,j,k)
            END DO
          END DO
        END DO

        nb = nb + ncoset(lb_max)

      END DO

      na = na + ncoset(la_max)

    END DO

  END SUBROUTINE overlap3_old

!   ***************************************************************************

  SUBROUTINE overlap3(la_max,npgfa,zeta,rpgfa,la_min,&
                      lb_max,npgfb,zetb,rpgfb,lb_min,&
                      lc_max,npgfc,zetc,rpgfc,lc_min,&
                      rab,dab,rac,dac,rbc,dbc,sabc)

!   Purpose: Calculation of three-center overlap integrals <a|b|c> over
!            Cartesian Gaussian functions. 

!   History: - Creation (08.2011, JGH)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: la_min, lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb, rpgfb
    INTEGER, INTENT(IN)                      :: lb_min, lc_max, npgfc
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetc, rpgfc
    INTEGER, INTENT(IN)                      :: lc_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: dab
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac
    REAL(KIND=dp), INTENT(IN)                :: dac
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rbc
    REAL(KIND=dp), INTENT(IN)                :: dbc
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: sabc

    INTEGER                                  :: i, ipgf, ix, j, jpgf, k, &
                                                kpgf, la, lam, lamm, lb, lbm, &
                                                lbmm, lc, lcm, lcmm, na, nb, &
                                                nc
    REAL(KIND=dp)                            :: aai, aax, aay, f0, f1, f2, &
                                                rcp2, zetg, zetp
    REAL(KIND=dp), DIMENSION(3)              :: rag, rbg, rcg, rcp
    REAL(KIND=dp), DIMENSION(ncoset(la_max),&
      ncoset(lb_max), ncoset(lc_max))        :: s

!   ---------------------------------------------------------------------------
!   *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0
    DO ipgf=1,npgfa

      nb = 0
      DO jpgf=1,npgfb

        ! *** Screening ***
        IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
          nb = nb + ncoset(lb_max)
          CYCLE
        END IF

        nc = 0
        DO kpgf=1,npgfc

          ! *** Screening ***
          IF ((rpgfb(jpgf) + rpgfc(kpgf) < dbc).OR.&
              (rpgfa(ipgf) + rpgfc(kpgf) < dac)) THEN
            nc = nc + ncoset(lc_max)
            CYCLE
          END IF

          ! *** Calculate some prefactors ***
          zetg = 1.0_dp/(zeta(ipgf) + zetb(jpgf) + zetc(kpgf))
          zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))
          f0 = (pi*zetg)**1.5_dp
          f1 = zetb(jpgf)*zetp
          f2 = 0.5_dp*zetg
          rcp(:) = f1*rab(:) - rac(:)
          rcp2 = rcp(1)*rcp(1) + rcp(2)*rcp(2) + rcp(3)*rcp(3)
          rag(:) =  zetg*(zetb(jpgf)*rab(:) + zetc(kpgf)*rac(:))
          rbg(:) = -zetg*(zeta(ipgf)*rab(:) - zetc(kpgf)*rbc(:))
          rcg(:) = -zetg*(zeta(ipgf)*rac(:) + zetb(jpgf)*rbc(:))

          ! *** Calculate the basic three-center overlap integral [s|s|s] ***
          s(1,1,1) = f0*EXP(-(zeta(ipgf)*f1*dab*dab + zetc(kpgf)*zetg*rcp2/zetp))

          ! *** Recurrence steps: [s|s|s] -> [s|s|c] ***
          IF (lc_max > 0) THEN
            ! *** [s|s|p] = (Gi - Ci)*[s|s|s]  (i = x,y,z) ***
            s(1,1,2) = rcg(1)*s(1,1,1)
            s(1,1,3) = rcg(2)*s(1,1,1)
            s(1,1,4) = rcg(3)*s(1,1,1)
            ! *** [s|s|c] = (Gi - Ci)*[s|s|c-1i] + f2*Ni(c-1i)*[s|s|c-2i] ***
            DO lc=5,ncoset(lc_max)
               CALL indices(lc,ix,lcm,lcmm,aai)
               !
               s(1,1,lc) = rcg(ix)*s(1,1,lcm) + f2*aai*s(1,1,lcmm)
               !
            END DO
          END IF
          ! *** Recurrence steps: [s|s|c] -> [s|b|c] ***
          IF (lb_max > 0) THEN
            DO lc=1,ncoset(lc_max)
              DO lb=2,ncoset(lb_max)
                 CALL indices(lb,ix,lbm,lbmm,aai)
                 CALL reduce_index(lc,ix,lcm,aax)
                 !
                 s(1,lb,lc) = rbg(ix)*s(1,lbm,lc) + f2*aai*s(1,lbmm,lc) + f2*aax*s(1,lbm,lcm)
                 !
              END DO
            END DO
          END IF
          ! *** Recurrence steps: [s|b|c] -> [a|b|c] ***
          IF (la_max > 0) THEN
            DO lc=1,ncoset(lc_max)
              DO lb=1,ncoset(lb_max)
                DO la=2,ncoset(la_max)
                   CALL indices(la,ix,lam,lamm,aai)
                   CALL reduce_index(lb,ix,lbm,aax)
                   CALL reduce_index(lc,ix,lcm,aay)
                   !
                   s(la,lb,lc) = rag(ix)*s(lam,lb,lc) + &
                        f2*aai*s(lamm,lb,lc) + f2*aax*s(lam,lbm,lc) + f2*aay*s(lam,lb,lcm)
                   !
                END DO
              END DO
            END DO
          END IF

          DO k=ncoset(lc_min-1)+1,ncoset(lc_max)
            DO j=ncoset(lb_min-1)+1,ncoset(lb_max)
              DO i=ncoset(la_min-1)+1,ncoset(la_max)
                sabc(na+i,nb+j,nc+k) = s(i,j,k)
              END DO
            END DO
          END DO

          nc = nc + ncoset(lc_max)
        END DO

        nb = nb + ncoset(lb_max)
      END DO

      na = na + ncoset(la_max)
    END DO

  END SUBROUTINE overlap3

  SUBROUTINE indices(lc,ix,lcm,lcmm,aai)
    INTEGER, INTENT(IN)                      :: lc
    INTEGER, INTENT(OUT)                     :: ix, lcm, lcmm
    REAL(KIND=dp), INTENT(OUT)               :: aai

    INTEGER, DIMENSION(1)                    :: ixx
    INTEGER, DIMENSION(3)                    :: ii

    ii(:) = indco(1:3,lc)
    ixx = MAXLOC(ii)
    ix = ixx(1)
    ii(ix) = ii(ix) - 1
    aai = REAL(ii(ix),dp)
    lcm = coset(ii(1),ii(2),ii(3))
    ii(ix) = ii(ix) - 1
    lcmm = coset(ii(1),ii(2),ii(3))

  END SUBROUTINE indices

  SUBROUTINE reduce_index(lx,ix,lxm,aax)
    INTEGER, INTENT(IN)                      :: lx, ix
    INTEGER, INTENT(OUT)                     :: lxm
    REAL(KIND=dp), INTENT(OUT)               :: aax

    INTEGER, DIMENSION(3)                    :: ii

    ii(:) = indco(1:3,lx)
    aax = REAL(ii(ix),dp)
    ii(ix) = ii(ix) - 1
    lxm = coset(ii(1),ii(2),ii(3))

  END SUBROUTINE reduce_index 

! *****************************************************************************

END MODULE ai_overlap3

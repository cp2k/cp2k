!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/ai_overlap_ppl [1.0] *
!!
!!   NAME
!!     ai_overlap_ppl
!!
!!   FUNCTION
!!     Calculation of three-center overlap integrals over Cartesian
!!     Gaussian-type functions for the local part of the Goedecker
!!     pseudopotential (GTH).
!!
!!     <a|V(local)|b> = <a|V(erf) + V(ppl)|b>
!!                    = <a|V(erf)|b> + <a|V(ppl)|b>
!!                    = <a|-Z(eff)*erf(SQRT(2)*alpha*r)/r +
!!                      (C1 + C2*(alpha*r)**2 + C3*(alpha*r)**4 +
!!                       C4*(alpha*r)**6)*exp(-(alpha*r)**2/2))|b>
!!
!!   AUTHOR
!!     Matthias Krack (04.10.2000)
!!
!!   LITERATURE
!!     S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!!     S. Goedecker, M. Teter and J. Hutter, Phys. Rev. B 54, 1703 (1996)
!!     C. Hartwigsen, S. Goedecker and J. Hutter, Phys. Rev. B 58, 3641 (1998)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE ai_overlap_ppl

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
! rab2       : Square of the distance between the atomic centers a and b.
! rac        : Distance vector between the atomic centers a and c.
! rac2       : Square of the distance between the atomic centers a and c.
! rbc        : Distance vector between the atomic centers b and c.
! rbc2       : Square of the distance between the atomic centers b and c.
! rpgf{a,b,c}: Radius of the primitive Gaussian-type function a or b.
! zet{a,b,c} : Exponents of the Gaussian-type functions a or b.
! zetg       : Reciprocal of the sum of the exponents of orbital a, b and c.
! zetp       : Reciprocal of the sum of the exponents of orbital a and b.

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE mathconstants,    ONLY: pi
  USE orbital_pointers, ONLY: coset,ncoset

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: overlap_ppl

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE overlap_ppl(la_max,zeta,rpgfa,la_min,&
                         lb_max,zetb,rpgfb,lb_min,&
                         cexp_ppl,zetc,rpgfc,&
                         rab,rab2,rac,rac2,rbc,rbc2,vab)

!   Purpose: Calculation of three-center overlap integrals <a|c|b> over
!            Cartesian Gaussian functions for the local part of the Goedecker
!            pseudopotential (GTH). c is a primitive Gaussian-type function
!            with a set of even angular momentum indices.

!            <a|V(ppl)|b> = <a| (C1 + C2*(alpha*r)**2 + C3*(alpha*r)**4 +
!                                C4*(alpha*r)**6)*exp(-(alpha*r)**2/2))|b>

!   History: - Creation (04.10.2000, Matthias Krack)

!   ***************************************************************************

    REAL(wp), INTENT(IN)                      :: rab2,rac2,rbc2,rpgfc,zetc
    INTEGER, INTENT(IN)                       :: la_max,la_min,lb_max,lb_min
    REAL(wp), DIMENSION(3), INTENT(IN)        :: rab,rac,rbc
    REAL(wp), DIMENSION(:), INTENT(IN)        :: cexp_ppl,rpgfa,rpgfb,zeta,zetb
    REAL(wp), DIMENSION(:,:), INTENT(INOUT)   :: vab

!   *** Local variables ***

    REAL(wp) :: dab,dac,dbc,f0,f1,f2,f3,fc,fx,fy,fz,rcg2,rcp2,zetg,zetp
    INTEGER  :: ax,ay,az,bx,by,bz,c,i,ipgf,j,jpgf,la,la_start,lb,lc,na,nb

    REAL(wp), DIMENSION(3) :: rag,rbg,rcg,rcp

    REAL(wp), DIMENSION(ncoset(la_max),ncoset(lb_max),SIZE(cexp_ppl)) :: sabc

!   ---------------------------------------------------------------------------

!   *** Calculate the distances of the centers a, b and c ***

    dab = SQRT(rab2)
    dac = SQRT(rac2)
    dbc = SQRT(rbc2)

!   *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0

    DO ipgf=1,SIZE(zeta)

!     *** Screening ***

      IF (rpgfa(ipgf) + rpgfc < dac) THEN
        na = na + ncoset(la_max)
        CYCLE
      END IF

      nb = 0

      DO jpgf=1,SIZE(zetb)

!       *** Screening ***

        IF ((rpgfb(jpgf) + rpgfc < dbc).OR.&
            (rpgfa(ipgf) + rpgfb(jpgf) < dab)) THEN
          nb = nb + ncoset(lb_max)
          CYCLE
        END IF

!       *** Calculate some prefactors ***

        zetp = 1.0_wp/(zeta(ipgf) + zetb(jpgf))
        zetg = 1.0_wp/(zeta(ipgf) + zetb(jpgf) + zetc)

        f0 = (pi*zetg)**1.5_wp
        f1 = zetb(jpgf)*zetp
        f2 = 0.5_wp*zetg

        rcp(:) = f1*rab(:) - rac(:)

        rcp2 = rcp(1)*rcp(1) + rcp(2)*rcp(2) + rcp(3)*rcp(3)

!       *** Calculate the basic three-center overlap integral [s|s|s] ***

        sabc(1,1,1) = f0*EXP(-(zeta(ipgf)*f1*rab2 + zetc*zetg*rcp2/zetp))

!       *** Recurrence steps: [s|s|s] -> [a|s|s] ***

        IF (la_max > 0) THEN

!         *** Vertical recurrence steps: [s|s|s] -> [a|s|s] ***

          rag(:) = zetg*(zetb(jpgf)*rab(:) + zetc*rac(:))

!         *** [p|s|s] = (Gi - Ai)*[s|s|s]  (i = x,y,z) ***

          sabc(2,1,1) = rag(1)*sabc(1,1,1)
          sabc(3,1,1) = rag(2)*sabc(1,1,1)
          sabc(4,1,1) = rag(3)*sabc(1,1,1)

!         *** [a|s|s] = (Gi - Ai)*[a-1i|s|s] + f2*Ni(a-1i)*[a-2i|s|s] ***

          DO la=2,la_max

!           *** Increase the angular momentum component z of function a ***

            sabc(coset(0,0,la),1,1) =&
              rag(3)*sabc(coset(0,0,la-1),1,1) +&
              f2*REAL(la-1,wp)*sabc(coset(0,0,la-2),1,1)

!           *** Increase the angular momentum component y of function a ***

            az = la - 1
            sabc(coset(0,1,az),1,1) = rag(2)*sabc(coset(0,0,az),1,1)

            DO ay=2,la
              az = la - ay
              sabc(coset(0,ay,az),1,1) =&
                rag(2)*sabc(coset(0,ay-1,az),1,1) +&
                f2*REAL(ay-1,wp)*sabc(coset(0,ay-2,az),1,1)
            END DO

!           *** Increase the angular momentum component x of function a ***

            DO ay=0,la-1
              az = la - 1 - ay
              sabc(coset(1,ay,az),1,1) = rag(1)*sabc(coset(0,ay,az),1,1)
            END DO

            DO ax=2,la
              f3 = f2*REAL(ax-1,wp)
              DO ay=0,la-ax
                az = la - ax - ay
                sabc(coset(ax,ay,az),1,1) =&
                  rag(1)*sabc(coset(ax-1,ay,az),1,1) +&
                  f3*sabc(coset(ax-2,ay,az),1,1)
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
                  sabc(coset(ax,ay,az),2,1) =&
                    sabc(coset(ax+1,ay,az),1,1) -&
                    rab(1)*sabc(coset(ax,ay,az),1,1)
                  sabc(coset(ax,ay,az),3,1) =&
                    sabc(coset(ax,ay+1,az),1,1) -&
                    rab(2)*sabc(coset(ax,ay,az),1,1)
                  sabc(coset(ax,ay,az),4,1) =&
                    sabc(coset(ax,ay,az+1),1,1) -&
                    rab(3)*sabc(coset(ax,ay,az),1,1)
                END DO
              END DO
            END DO

!           *** Vertical recurrence step ***

!           *** [a|s|p] = (Gi - Bi)*[a|s|s] + f2*Ni(a)*[a-1i|s|s] ***

            DO ax=0,la_max
              fx = f2*REAL(ax,wp)
              DO ay=0,la_max-ax
                fy = f2*REAL(ay,wp)
                az = la_max - ax - ay
                fz = f2*REAL(az,wp)
                IF (ax == 0) THEN
                  sabc(coset(ax,ay,az),2,1) =&
                    rbg(1)*sabc(coset(ax,ay,az),1,1)
                ELSE
                  sabc(coset(ax,ay,az),2,1) =&
                    rbg(1)*sabc(coset(ax,ay,az),1,1) +&
                    fx*sabc(coset(ax-1,ay,az),1,1)
                END IF
                IF (ay == 0) THEN
                  sabc(coset(ax,ay,az),3,1) =&
                    rbg(2)*sabc(coset(ax,ay,az),1,1)
                ELSE
                  sabc(coset(ax,ay,az),3,1) =&
                    rbg(2)*sabc(coset(ax,ay,az),1,1) +&
                    fy*sabc(coset(ax,ay-1,az),1,1)
                END IF
                IF (az == 0) THEN
                  sabc(coset(ax,ay,az),4,1) =&
                    rbg(3)*sabc(coset(ax,ay,az),1,1)
                ELSE
                  sabc(coset(ax,ay,az),4,1) =&
                    rbg(3)*sabc(coset(ax,ay,az),1,1) +&
                    fz*sabc(coset(ax,ay,az-1),1,1)
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

!                   *** Shift of angular momentum component z from a to b ***

                    sabc(coset(ax,ay,az),coset(0,0,lb),1) =&
                      sabc(coset(ax,ay,az+1),coset(0,0,lb-1),1) -&
                      rab(3)*sabc(coset(ax,ay,az),coset(0,0,lb-1),1)

!                   *** Shift of angular momentum component y from a to b ***

                    DO by=1,lb
                      bz = lb - by
                      sabc(coset(ax,ay,az),coset(0,by,bz),1) =&
                        sabc(coset(ax,ay+1,az),coset(0,by-1,bz),1) -&
                        rab(2)*sabc(coset(ax,ay,az),coset(0,by-1,bz),1)
                    END DO

!                   *** Shift of angular momentum component x from a to b ***

                    DO bx=1,lb
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        sabc(coset(ax,ay,az),coset(bx,by,bz),1) =&
                          sabc(coset(ax+1,ay,az),coset(bx-1,by,bz),1) -&
                          rab(1)*sabc(coset(ax,ay,az),coset(bx-1,by,bz),1)
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
                fx = f2*REAL(ax,wp)
                DO ay=0,la_max-ax
                  fy = f2*REAL(ay,wp)
                  az = la_max - ax - ay
                  fz = f2*REAL(az,wp)

!                 *** Shift of angular momentum component z from a to b ***

                  f3 = f2*REAL(lb-1,wp)

                  IF (az == 0) THEN
                    sabc(coset(ax,ay,az),coset(0,0,lb),1) =&
                      rbg(3)*sabc(coset(ax,ay,az),coset(0,0,lb-1),1) +&
                      f3*sabc(coset(ax,ay,az),coset(0,0,lb-2),1)
                  ELSE
                    sabc(coset(ax,ay,az),coset(0,0,lb),1) =&
                      rbg(3)*sabc(coset(ax,ay,az),coset(0,0,lb-1),1) +&
                      fz*sabc(coset(ax,ay,az-1),coset(0,0,lb-1),1) +&
                      f3*sabc(coset(ax,ay,az),coset(0,0,lb-2),1)
                  END IF

!                 *** Shift of angular momentum component y from a to b ***

                  IF (ay == 0) THEN
                    bz = lb - 1
                    sabc(coset(ax,ay,az),coset(0,1,bz),1) =&
                      rbg(2)*sabc(coset(ax,ay,az),coset(0,0,bz),1)
                    DO by=2,lb
                      bz = lb - by
                      f3 = f2*REAL(by-1,wp)
                      sabc(coset(ax,ay,az),coset(0,by,bz),1) =&
                        rbg(2)*sabc(coset(ax,ay,az),coset(0,by-1,bz),1) +&
                        f3*sabc(coset(ax,ay,az),coset(0,by-2,bz),1)
                    END DO
                  ELSE
                    bz = lb - 1
                    sabc(coset(ax,ay,az),coset(0,1,bz),1) =&
                      rbg(2)*sabc(coset(ax,ay,az),coset(0,0,bz),1) +&
                      fy*sabc(coset(ax,ay-1,az),coset(0,0,bz),1)
                    DO by=2,lb
                      bz = lb - by
                      f3 = f2*REAL(by-1,wp)
                      sabc(coset(ax,ay,az),coset(0,by,bz),1) =&
                        rbg(2)*sabc(coset(ax,ay,az),coset(0,by-1,bz),1) +&
                        fy*sabc(coset(ax,ay-1,az),coset(0,by-1,bz),1) +&
                        f3*sabc(coset(ax,ay,az),coset(0,by-2,bz),1)
                    END DO
                  END IF

!                 *** Shift of angular momentum component x from a to b ***

                  IF (ax == 0) THEN
                    DO by=0,lb-1
                      bz = lb - 1 - by
                      sabc(coset(ax,ay,az),coset(1,by,bz),1) =&
                        rbg(1)*sabc(coset(ax,ay,az),coset(0,by,bz),1)
                    END DO
                    DO bx=2,lb
                      f3 = f2*REAL(bx-1,wp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        sabc(coset(ax,ay,az),coset(bx,by,bz),1) =&
                          rbg(1)*sabc(coset(ax,ay,az),coset(bx-1,by,bz),1) +&
                          f3*sabc(coset(ax,ay,az),coset(bx-2,by,bz),1)
                      END DO
                    END DO
                  ELSE
                    DO by=0,lb-1
                      bz = lb - 1 - by
                      sabc(coset(ax,ay,az),coset(1,by,bz),1) =&
                        rbg(1)*sabc(coset(ax,ay,az),coset(0,by,bz),1) +&
                        fx*sabc(coset(ax-1,ay,az),coset(0,by,bz),1)
                    END DO
                    DO bx=2,lb
                      f3 = f2*REAL(bx-1,wp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        sabc(coset(ax,ay,az),coset(bx,by,bz),1) =&
                          rbg(1)*sabc(coset(ax,ay,az),coset(bx-1,by,bz),1) +&
                          fx*sabc(coset(ax-1,ay,az),coset(bx-1,by,bz),1) +&
                          f3*sabc(coset(ax,ay,az),coset(bx-2,by,bz),1)
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

            sabc(1,2,1) = rbg(1)*sabc(1,1,1)
            sabc(1,3,1) = rbg(2)*sabc(1,1,1)
            sabc(1,4,1) = rbg(3)*sabc(1,1,1)

!           *** [s|s|b] = (Gi - Bi)*[s|s|b-1i] + f2*Ni(b-1i)*[s|s|b-2i] ***

            DO lb=2,lb_max

!             *** Increase the angular momentum component z of function b ***

              sabc(1,coset(0,0,lb),1) =&
                rbg(3)*sabc(1,coset(0,0,lb-1),1) +&
                f2*REAL(lb-1,wp)*sabc(1,coset(0,0,lb-2),1)

!             *** Increase the angular momentum component y of function b ***

              bz = lb - 1
              sabc(1,coset(0,1,bz),1) = rbg(2)*sabc(1,coset(0,0,bz),1)

              DO by=2,lb
                bz = lb - by
                sabc(1,coset(0,by,bz),1) =&
                  rbg(2)*sabc(1,coset(0,by-1,bz),1) +&
                  f2*REAL(by-1,wp)*sabc(1,coset(0,by-2,bz),1)
              END DO

!             *** Increase the angular momentum component x of function b ***

              DO by=0,lb-1
                bz = lb - 1 - by
                sabc(1,coset(1,by,bz),1) = rbg(1)*sabc(1,coset(0,by,bz),1)
              END DO

              DO bx=2,lb
                f3 = f2*REAL(bx-1,wp)
                DO by=0,lb-bx
                  bz = lb - bx - by
                  sabc(1,coset(bx,by,bz),1) =&
                    rbg(1)*sabc(1,coset(bx-1,by,bz),1) +&
                    f3*sabc(1,coset(bx-2,by,bz),1)
                END DO
              END DO

            END DO

          END IF

        END IF

!       *** Calculate the contributions for lc > 0 ***

        DO c=2,SIZE(cexp_ppl)

          lc = 2*(c - 1)
          fc = f2*REAL(lc,wp)

!         *** Calculate the basic three-center overlap integrals [s|c|s] ***

          SELECT CASE (c)
          CASE (2)
            rcg(:) = -zetg*(zeta(ipgf)*rac(:) + zetb(jpgf)*rbc(:))
            rcg2 = rcg(1)*rcg(1) + rcg(2)*rcg(2) + rcg(3)*rcg(3)
            sabc(1,1,2) = (1.5_wp*zetg + rcg2)*sabc(1,1,1)
          CASE (3)
            sabc(1,1,3) = (2.5_wp*zetg + rcg2)*sabc(1,1,2) +&
                          rcg2*zetg*sabc(1,1,1)
          CASE (4)
            sabc(1,1,4) = (3.5_wp*zetg + rcg2)*sabc(1,1,3) +&
                          2.0_wp*rcg2*zetg*(sabc(1,1,2) + zetg*sabc(1,1,1))
          END SELECT

!         *** Recurrence steps: [s|c|s] -> [a|c|s] ***

          IF (la_max > 0) THEN

!           *** Vertical recurrence steps: [s|c|s] -> [a|c|s] ***

!           *** [p|c|s] = (Gi - Ai)*[s|c|s] +                          ***
!           ***           f2*Ni(c)*([p|c-2i|s] - (Ci - Ai)*[s|c-2i|s]) ***

            sabc(2,1,c) = rag(1)*sabc(1,1,c) +&
                          fc*(sabc(2,1,c-1) - rac(1)*sabc(1,1,c-1))
            sabc(3,1,c) = rag(2)*sabc(1,1,c) +&
                          fc*(sabc(3,1,c-1) - rac(2)*sabc(1,1,c-1))
            sabc(4,1,c) = rag(3)*sabc(1,1,c) +&
                          fc*(sabc(4,1,c-1) - rac(3)*sabc(1,1,c-1))

!           *** [a|c|s] = (Gi - Ai)*[a-1i|c|s] +                          ***
!           ***           f2*Ni(a-1i)*[a-2i|c|s] +                        ***
!           ***           f2*Ni(c)*([a|c-2i|s] - (Ci - Ai)*[a-1i|c-2i|s]) ***

            DO la=2,la_max

!             *** Increase the angular momentum component z of function a ***

              sabc(coset(0,0,la),1,c) =&
                rag(3)*sabc(coset(0,0,la-1),1,c) +&
                f2*REAL(la-1,wp)*sabc(coset(0,0,la-2),1,c) +&
                fc*(sabc(coset(0,0,la),1,c-1) -&
                    rac(3)*sabc(coset(0,0,la-1),1,c-1))

!             *** Increase the angular momentum component y of function a ***

              az = la - 1
              sabc(coset(0,1,az),1,c) =&
                rag(2)*sabc(coset(0,0,az),1,c) +&
                fc*(sabc(coset(0,1,az),1,c-1) -&
                    rac(2)*sabc(coset(0,0,az),1,c-1))

              DO ay=2,la
                az = la - ay
                sabc(coset(0,ay,az),1,c) =&
                  rag(2)*sabc(coset(0,ay-1,az),1,c) +&
                  f2*REAL(ay-1,wp)*sabc(coset(0,ay-2,az),1,c) +&
                  fc*(sabc(coset(0,ay,az),1,c-1) -&
                      rac(2)*sabc(coset(0,ay-1,az),1,c-1))
              END DO

!             *** Increase the angular momentum component x of function a ***

              DO ay=0,la-1
                az = la - 1 - ay
                sabc(coset(1,ay,az),1,c) =&
                  rag(1)*sabc(coset(0,ay,az),1,c) +&
                  fc*(sabc(coset(1,ay,az),1,c-1) -&
                      rac(1)*sabc(coset(0,ay,az),1,c-1))
              END DO

              DO ax=2,la
                f3 = f2*REAL(ax-1,wp)
                DO ay=0,la-ax
                  az = la - ax - ay
                  sabc(coset(ax,ay,az),1,c) =&
                    rag(1)*sabc(coset(ax-1,ay,az),1,c) +&
                    f3*sabc(coset(ax-2,ay,az),1,c) +&
                    fc*(sabc(coset(ax,ay,az),1,c-1) -&
                        rac(1)*sabc(coset(ax-1,ay,az),1,c-1))
                END DO
              END DO

            END DO

!           *** Recurrence steps: [a|c|s] -> [a|c|b] ***

            IF (lb_max > 0) THEN

!             *** Horizontal recurrence steps ***

!             *** [a|c|p] = [a+1i|c|s] - (Bi - Ai)*[a|c|s] ***

              IF (lb_max == 1) THEN
                la_start = la_min
              ELSE
                la_start = MAX(0,la_min-1)
              END IF

              DO la=la_start,la_max-1
                DO ax=0,la
                  DO ay=0,la-ax
                    az = la - ax - ay
                    sabc(coset(ax,ay,az),2,c) =&
                      sabc(coset(ax+1,ay,az),1,c) -&
                      rab(1)*sabc(coset(ax,ay,az),1,c)
                    sabc(coset(ax,ay,az),3,c) =&
                      sabc(coset(ax,ay+1,az),1,c) -&
                      rab(2)*sabc(coset(ax,ay,az),1,c)
                    sabc(coset(ax,ay,az),4,c) =&
                      sabc(coset(ax,ay,az+1),1,c) -&
                      rab(3)*sabc(coset(ax,ay,az),1,c)
                  END DO
                END DO
              END DO

!             *** Vertical recurrence step ***

!             *** [a|c|p] = (Gi - Bi)*[a|c|s] +                          ***
!             ***           f2*Ni(a)*[a-1i|c|s] +                        ***
!             ***           f2*Ni(c)*([a|c-2i|p] - (Ci - Bi)*[a|c-2i|s]) ***

              DO ax=0,la_max
                fx = f2*REAL(ax,wp)
                DO ay=0,la_max-ax
                  fy = f2*REAL(ay,wp)
                  az = la_max - ax - ay
                  fz = f2*REAL(az,wp)
                  IF (ax == 0) THEN
                    sabc(coset(ax,ay,az),2,c) =&
                      rbg(1)*sabc(coset(ax,ay,az),1,c) +&
                      fc*(sabc(coset(ax,ay,az),2,c-1) -&
                          rbc(1)*sabc(coset(ax,ay,az),1,c-1))
                  ELSE
                    sabc(coset(ax,ay,az),2,c) =&
                      rbg(1)*sabc(coset(ax,ay,az),1,c) +&
                      fx*sabc(coset(ax-1,ay,az),1,c) +&
                      fc*(sabc(coset(ax,ay,az),2,c-1) -&
                          rbc(1)*sabc(coset(ax,ay,az),1,c-1))
                  END IF
                  IF (ay == 0) THEN
                    sabc(coset(ax,ay,az),3,c) =&
                      rbg(2)*sabc(coset(ax,ay,az),1,c) +&
                      fc*(sabc(coset(ax,ay,az),3,c-1) -&
                          rbc(2)*sabc(coset(ax,ay,az),1,c-1))
                  ELSE
                    sabc(coset(ax,ay,az),3,c) =&
                      rbg(2)*sabc(coset(ax,ay,az),1,c) +&
                      fy*sabc(coset(ax,ay-1,az),1,c) +&
                      fc*(sabc(coset(ax,ay,az),3,c-1) -&
                          rbc(2)*sabc(coset(ax,ay,az),1,c-1))
                  END IF
                  IF (az == 0) THEN
                    sabc(coset(ax,ay,az),4,c) =&
                      rbg(3)*sabc(coset(ax,ay,az),1,c) +&
                      fc*(sabc(coset(ax,ay,az),4,c-1) -&
                          rbc(3)*sabc(coset(ax,ay,az),1,c-1))
                  ELSE
                    sabc(coset(ax,ay,az),4,c) =&
                      rbg(3)*sabc(coset(ax,ay,az),1,c) +&
                      fz*sabc(coset(ax,ay,az-1),1,c) +&
                      fc*(sabc(coset(ax,ay,az),4,c-1) -&
                          rbc(3)*sabc(coset(ax,ay,az),1,c-1))
                  END IF
                END DO
              END DO

!             *** Recurrence steps: [a|s|p] -> [a|s|b] ***

              DO lb=2,lb_max

!               *** Horizontal recurrence steps ***

!               *** [a|c|b] = [a+1i|c|b-1i] - (Bi - Ai)*[a|c|b-1i] ***

                IF (lb == lb_max) THEN
                  la_start = la_min
                ELSE
                  la_start = MAX(0,la_min-1)
                END IF

                DO la=la_start,la_max-1
                  DO ax=0,la
                    DO ay=0,la-ax
                      az = la - ax - ay

!                     *** Shift of angular momentum component x from a to b ***

                      sabc(coset(ax,ay,az),coset(0,0,lb),c) =&
                        sabc(coset(ax,ay,az+1),coset(0,0,lb-1),c) -&
                        rab(3)*sabc(coset(ax,ay,az),coset(0,0,lb-1),c)

!                     *** Shift of angular momentum component y from a to b ***

                      DO by=1,lb
                        bz = lb - by
                        sabc(coset(ax,ay,az),coset(0,by,bz),c) =&
                          sabc(coset(ax,ay+1,az),coset(0,by-1,bz),c) -&
                          rab(2)*sabc(coset(ax,ay,az),coset(0,by-1,bz),c)
                      END DO

!                     *** Shift of angular momentum component z from a to b ***

                      DO bx=1,lb
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          sabc(coset(ax,ay,az),coset(bx,by,bz),c) =&
                            sabc(coset(ax+1,ay,az),coset(bx-1,by,bz),c) -&
                            rab(1)*sabc(coset(ax,ay,az),coset(bx-1,by,bz),c)
                        END DO
                      END DO

                    END DO
                  END DO
                END DO

!               *** Vertical recurrence step ***

!               *** [a|c|b] = (Gi - Bi)*[a|c|b-1i] +             ***
!               ***           f2*Ni(a)*[a-1i|c|b-1i] +           ***
!               ***           f2*Ni(b-1i)*[a|c|b-2i] +           ***
!               ***           f2*Ni(c)*([a|c-2i|b] -             ***
!               ***                     (Ci - Bi)*[a|c-2i|b-1i]) ***

                DO ax=0,la_max
                  fx = f2*REAL(ax,wp)
                  DO ay=0,la_max-ax
                    fy = f2*REAL(ay,wp)
                    az = la_max - ax - ay
                    fz = f2*REAL(az,wp)

!                   *** Shift of angular momentum component z from a to b ***

                    f3 = f2*REAL(lb-1,wp)

                    IF (az == 0) THEN
                      sabc(coset(ax,ay,az),coset(0,0,lb),c) =&
                        rbg(3)*sabc(coset(ax,ay,az),coset(0,0,lb-1),c) +&
                        f3*sabc(coset(ax,ay,az),coset(0,0,lb-2),c) +&
                        fc*(sabc(coset(ax,ay,az),coset(0,0,lb),c-1) -&
                            rbc(3)*sabc(coset(ax,ay,az),coset(0,0,lb-1),c-1))
                    ELSE
                      sabc(coset(ax,ay,az),coset(0,0,lb),c) =&
                        rbg(3)*sabc(coset(ax,ay,az),coset(0,0,lb-1),c) +&
                        fz*sabc(coset(ax,ay,az-1),coset(0,0,lb-1),c) +&
                        f3*sabc(coset(ax,ay,az),coset(0,0,lb-2),c) +&
                        fc*(sabc(coset(ax,ay,az),coset(0,0,lb),c-1) -&
                            rbc(3)*sabc(coset(ax,ay,az),coset(0,0,lb-1),c-1))
                    END IF

!                   *** Shift of angular momentum component y from a to b ***

                    IF (ay == 0) THEN
                      bz = lb - 1
                      sabc(coset(ax,ay,az),coset(0,1,bz),c) =&
                        rbg(2)*sabc(coset(ax,ay,az),coset(0,0,bz),c) +&
                        fc*(sabc(coset(ax,ay,az),coset(0,1,bz),c-1) -&
                            rbc(2)*sabc(coset(ax,ay,az),coset(0,0,bz),c-1))
                      DO by=2,lb
                        bz = lb - by
                        f3 = f2*REAL(by-1,wp)
                        sabc(coset(ax,ay,az),coset(0,by,bz),c) =&
                          rbg(2)*sabc(coset(ax,ay,az),coset(0,by-1,bz),c) +&
                          f3*sabc(coset(ax,ay,az),coset(0,by-2,bz),c) +&
                          fc*(sabc(coset(ax,ay,az),coset(0,by,bz),c-1) -&
                              rbc(2)*sabc(coset(ax,ay,az),&
                                          coset(0,by-1,bz),c-1))
                      END DO
                    ELSE
                      bz = lb - 1
                      sabc(coset(ax,ay,az),coset(0,1,bz),c) =&
                        rbg(2)*sabc(coset(ax,ay,az),coset(0,0,bz),c) +&
                        fy*sabc(coset(ax,ay-1,az),coset(0,0,bz),c) +&
                        fc*(sabc(coset(ax,ay,az),coset(0,1,bz),c-1) -&
                            rbc(2)*sabc(coset(ax,ay,az),coset(0,0,bz),c-1))
                      DO by=2,lb
                        bz = lb - by
                        f3 = f2*REAL(by-1,wp)
                        sabc(coset(ax,ay,az),coset(0,by,bz),c) =&
                          rbg(2)*sabc(coset(ax,ay,az),coset(0,by-1,bz),c) +&
                          fy*sabc(coset(ax,ay-1,az),coset(0,by-1,bz),c) +&
                          f3*sabc(coset(ax,ay,az),coset(0,by-2,bz),c) +&
                          fc*(sabc(coset(ax,ay,az),coset(0,by,bz),c-1) -&
                              rbc(2)*sabc(coset(ax,ay,az),&
                                          coset(0,by-1,bz),c-1))
                      END DO
                    END IF

!                   *** Shift of angular momentum component x from a to b ***

                    IF (ax == 0) THEN
                      DO by=0,lb-1
                        bz = lb - 1 - by
                        sabc(coset(ax,ay,az),coset(1,by,bz),c) =&
                          rbg(1)*sabc(coset(ax,ay,az),coset(0,by,bz),c) +&
                          fc*(sabc(coset(ax,ay,az),coset(1,by,bz),c-1) -&
                              rbc(1)*sabc(coset(ax,ay,az),coset(0,by,bz),c-1))
                      END DO
                      DO bx=2,lb
                        f3 = f2*REAL(bx-1,wp)
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          sabc(coset(ax,ay,az),coset(bx,by,bz),c) =&
                            rbg(1)*sabc(coset(ax,ay,az),coset(bx-1,by,bz),c) +&
                            f3*sabc(coset(ax,ay,az),coset(bx-2,by,bz),c) +&
                            fc*(sabc(coset(ax,ay,az),coset(bx,by,bz),c-1) -&
                                rbc(1)*sabc(coset(ax,ay,az),&
                                            coset(bx-1,by,bz),c-1))
                        END DO
                      END DO
                    ELSE
                      DO by=0,lb-1
                        bz = lb - 1 - by
                        sabc(coset(ax,ay,az),coset(1,by,bz),c) =&
                          rbg(1)*sabc(coset(ax,ay,az),coset(0,by,bz),c) +&
                          fx*sabc(coset(ax-1,ay,az),coset(0,by,bz),c) +&
                          fc*(sabc(coset(ax,ay,az),coset(1,by,bz),c-1) -&
                                rbc(1)*sabc(coset(ax,ay,az),&
                                            coset(0,by,bz),c-1))
                      END DO
                      DO bx=2,lb
                        f3 = f2*REAL(bx-1,wp)
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          sabc(coset(ax,ay,az),coset(bx,by,bz),c) =&
                            rbg(1)*sabc(coset(ax,ay,az),coset(bx-1,by,bz),c) +&
                            fx*sabc(coset(ax-1,ay,az),coset(bx-1,by,bz),c) +&
                            f3*sabc(coset(ax,ay,az),coset(bx-2,by,bz),c) +&
                            fc*(sabc(coset(ax,ay,az),coset(bx,by,bz),c-1) -&
                                rbc(1)*sabc(coset(ax,ay,az),&
                                            coset(bx-1,by,bz),c-1))
                        END DO
                      END DO
                    END IF

                  END DO
                END DO

              END DO

            END IF

          ELSE

            IF (lb_max > 0) THEN

!             *** Vertical recurrence steps: [s|c|s] -> [s|c|b] ***

!             *** [s|c|p] = (Gi - Bi)*[s|c|s] +                          ***
!             ***           f2*Ni(c)*([s|c-2i|p] - (Ci - Bi)*[s|c-2i|s]) ***

              sabc(1,2,c) = rbg(1)*sabc(1,1,c) +&
                            fc*(sabc(1,2,c-1) - rbc(1)*sabc(1,1,c-1))
              sabc(1,3,c) = rbg(2)*sabc(1,1,c) +&
                            fc*(sabc(1,3,c-1) - rbc(2)*sabc(1,1,c-1))
              sabc(1,4,c) = rbg(3)*sabc(1,1,c) +&
                            fc*(sabc(1,4,c-1) - rbc(3)*sabc(1,1,c-1))

!             *** [s|c|b] = (Gi - Bi)*[s|c|b-1i] +                          ***
!             ***           f2*Ni(b-1i)*[s|c|b-2i] +                        ***
!             ***           f2*Ni(c)*([s|c-2i|b] - (Ci - Bi)*[s|c-2i|b-1i]) ***

              DO lb=2,lb_max

!               *** Increase the angular momentum component z of function b ***

                sabc(1,coset(0,0,lb),c) =&
                  rbg(3)*sabc(1,coset(0,0,lb-1),c) +&
                  f2*REAL(lb-1,wp)*sabc(1,coset(0,0,lb-2),c) +&
                  fc*(sabc(1,coset(0,0,lb),c-1) -&
                      rbc(3)*sabc(1,coset(0,0,lb-1),c-1))

!               *** Increase the angular momentum component y of function b ***

                bz = lb - 1
                sabc(1,coset(0,1,bz),c) =&
                  rbg(2)*sabc(1,coset(0,0,bz),c) +&
                  fc*(sabc(1,coset(0,1,bz),c-1) -&
                      rbc(2)*sabc(1,coset(0,0,bz),c-1))

                DO by=2,lb
                  bz = lb - by
                  sabc(1,coset(0,by,bz),c) =&
                    rbg(2)*sabc(1,coset(0,by-1,bz),c) +&
                    f2*REAL(by-1,wp)*sabc(1,coset(0,by-2,bz),c) +&
                    fc*(sabc(1,coset(0,by,bz),c-1) -&
                        rbc(2)*sabc(1,coset(0,by-1,bz),c-1))
                END DO

!               *** Increase the angular momentum component x of function b ***

                DO by=0,lb-1
                  bz = lb - 1 - by
                  sabc(1,coset(1,by,bz),c) =&
                    rbg(1)*sabc(1,coset(0,by,bz),c) +&
                    fc*(sabc(1,coset(1,by,bz),c-1) -&
                        rbc(1)*sabc(1,coset(0,by,bz),c-1))
                END DO

                DO bx=2,lb
                  f3 = f2*REAL(bx-1,wp)
                  DO by=0,lb-bx
                    bz = lb - bx - by
                    sabc(1,coset(bx,by,bz),c) =&
                      rbg(1)*sabc(1,coset(bx-1,by,bz),c) +&
                      f3*sabc(1,coset(bx-2,by,bz),c) +&
                      fc*(sabc(1,coset(bx,by,bz),c-1) -&
                          rbc(1)*sabc(1,coset(bx-1,by,bz),c-1))
                  END DO
                END DO

              END DO

            END IF

          END IF

        END DO

!       *** Add the contribution of the current pair ***
!       *** of primitive Gaussian-type functions     ***

        DO c=1,SIZE(cexp_ppl)
          DO j=ncoset(lb_min-1)+1,ncoset(lb_max)
            DO i=ncoset(la_min-1)+1,ncoset(la_max)
              vab(na+i,nb+j) = vab(na+i,nb+j) + cexp_ppl(c)*sabc(i,j,c)
            END DO
          END DO
        END DO

        nb = nb + ncoset(lb_max)

      END DO

      na = na + ncoset(la_max)

    END DO

  END SUBROUTINE overlap_ppl

! *****************************************************************************

END MODULE ai_overlap_ppl

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/ai_kinetic [1.0] *
!!
!!   NAME
!!     ai_kinetic
!!
!!   FUNCTION
!!     Calculation of the kinetic energy integrals over Cartesian
!!     Gaussian-type functions.
!!
!!     [a|T|b] = [a|-nabla**2/2|b]
!!
!!   AUTHOR
!!     Matthias Krack (31.07.2000)
!!
!!   LITERATURE
!!     S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE ai_kinetic

! *****************************************************************************

! ax,ay,az  : Angular momentum index numbers of orbital a.
! bx,by,bz  : Angular momentum index numbers of orbital b.
! coset     : Cartesian orbital set pointer.
! dab       : Distance between the atomic centers a and b.
! l{a,b}    : Angular momentum quantum number of shell a or b.
! l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
! l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
! ncoset    : Number of Cartesian orbitals up to l.
! rab       : Distance vector between the atomic centers a and b.
! rab2      : Square of the distance between the atomic centers a and b.
! rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
! sab       : Shell set of overlap integrals.
! tab       : Shell set of kinetic energy integrals.
! zet{a,b}  : Exponents of the Gaussian functions.

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE mathconstants,    ONLY: pi
  USE orbital_pointers, ONLY: coset,ncoset

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: kinetic

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE kinetic(la_max,zeta,rpgfa,la_min,&
                     lb_max,zetb,rpgfb,lb_min,&
                     rab,rab2,sab,tab)

!   Purpose: Calculation of the two-center kinetic energy integrals [a|T|b]
!            over Cartesian Gaussian-type functions.

!   History: - Creation (22.09.2000, Matthias Krack)

!   ***************************************************************************

    REAL(wp), INTENT(IN)                  :: rab2
    INTEGER, INTENT(IN)                   :: la_max,la_min,lb_max,lb_min
    REAL(wp), DIMENSION(3), INTENT(IN)    :: rab
    REAL(wp), DIMENSION(:), INTENT(IN)    :: rpgfa,rpgfb,zeta,zetb
    REAL(wp), DIMENSION(:,:), INTENT(IN)  :: sab
    REAL(wp), DIMENSION(:,:), INTENT(OUT) :: tab

!   *** Local variables ***

    REAL(wp) :: dab,f1,f2,f3,f4,f5,f6,f7,fx,fy,fz,xi,zetp
    INTEGER  :: ax,ay,az,bx,by,bz,i,ipgf,j,jpgf,la,la_start,lb,na,nb

    REAL(wp), DIMENSION(3) :: rap,rbp

    REAL(wp), DIMENSION(ncoset(la_max),ncoset(lb_max)) :: t

!   ---------------------------------------------------------------------------

!   *** Calculate the distance between the centers a and b ***

    dab = SQRT(rab2)

!   *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0

    DO ipgf=1,SIZE(zeta)

      nb = 0

      DO jpgf=1,SIZE(zetb)

!       *** Screening ***

        IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
          DO j=nb+ncoset(lb_min-1)+1,nb+ncoset(lb_max)
            DO i=na+ncoset(la_min-1)+1,na+ncoset(la_max)
              tab(i,j) = 0.0_wp
            END DO
          END DO
          nb = nb + ncoset(lb_max)
          CYCLE
        END IF

!       *** Calculate some prefactors ***

        zetp = 1.0_wp/(zeta(ipgf) + zetb(jpgf))

        f1 = zetb(jpgf)*zetp
        f2 = 0.5_wp*zetp
        xi = zeta(ipgf)*f1
        f3 = 2.0_wp*xi

!       *** Calculate the basic two-center kinetic energy integral [s|T|s] ***

        t(1,1) = xi*(3.0_wp - f3*rab2)*sab(na+1,nb+1)

!       *** Recurrence steps: [s|T|s] -> [a|T|b] ***

        IF (la_max > 0) THEN

!         *** Vertical recurrence steps: [s|T|s] -> [a|T|s] ***

          rap(:) = f1*rab(:)

!         *** [p|T|s] = (Pi - Ai)*[s|T|s] + f3*[p|s]  (i = x,y,z) ***

          t(2,1) = rap(1)*t(1,1) + f3*sab(na+2,nb+1)
          t(3,1) = rap(2)*t(1,1) + f3*sab(na+3,nb+1)
          t(4,1) = rap(3)*t(1,1) + f3*sab(na+4,nb+1)

!         *** [a|T|s] = (Pi - Ai)*[a-1i|T|s] +   ***
!         ***           f2*Ni(a-1i)*[a-2i|T|s] + ***
!         ***           f3*[a|s] -               ***
!         ***           f1*Ni(a-1i)*[a-2i|s])    ***

          DO la=2,la_max

!           *** Increase the angular momentum component z of function a ***

            t(coset(0,0,la),1) =&
              rap(3)*t(coset(0,0,la-1),1) +&
              f2*REAL(la-1,wp)*t(coset(0,0,la-2),1) +&
              f3*sab(na+coset(0,0,la),nb+1) -&
              f1*REAL(la-1,wp)*sab(na+coset(0,0,la-2),nb+1)

!           *** Increase the angular momentum component y of function a ***

            az = la - 1
            t(coset(0,1,az),1) =&
              rap(2)*t(coset(0,0,az),1) +&
              f3*sab(na+coset(0,1,az),nb+1)

            DO ay=2,la
              az = la - ay
              t(coset(0,ay,az),1) =&
                rap(2)*t(coset(0,ay-1,az),1) +&
                f2*REAL(ay-1,wp)*t(coset(0,ay-2,az),1) +&
                f3*sab(na+coset(0,ay,az),nb+1) -&
                f1*REAL(ay-1,wp)*sab(na+coset(0,ay-2,az),nb+1)
            END DO

!           *** Increase the angular momentum component x of function a ***

            DO ay=0,la-1
              az = la - 1 - ay
              t(coset(1,ay,az),1) =&
                rap(1)*t(coset(0,ay,az),1) +&
                f3*sab(na+coset(1,ay,az),nb+1)
            END DO

            DO ax=2,la
              f6 = f2*REAL(ax-1,wp)
              f7 = f1*REAL(ax-1,wp)
              DO ay=0,la-ax
                az = la - ax - ay
                t(coset(ax,ay,az),1) =&
                  rap(1)*t(coset(ax-1,ay,az),1) +&
                  f6*t(coset(ax-2,ay,az),1) +&
                  f3*sab(na+coset(ax,ay,az),nb+1) -&
                  f7*sab(na+coset(ax-2,ay,az),nb+1)
              END DO
            END DO

          END DO

!         *** Recurrence steps: [a|T|s] -> [a|T|b] ***

          IF (lb_max > 0) THEN

!           *** Horizontal recurrence steps ***

            f4 = 2.0_wp*zetb(jpgf)
            f5 = -zeta(ipgf)*zetp

            rbp(:) = rap(:) - rab(:)

!           *** [a|T|p] = [a+1i|T|s] - (Bi - Ai)*[a|T|s] +  f4*[a|p] ***

            IF (lb_max == 1) THEN
              la_start = la_min
            ELSE
              la_start = MAX(0,la_min-1)
            END IF

            DO la=la_start,la_max-1
              DO ax=0,la
                DO ay=0,la-ax
                  az = la - ax - ay
                  t(coset(ax,ay,az),2) = t(coset(ax+1,ay,az),1) -&
                                         rab(1)*t(coset(ax,ay,az),1) +&
                                         f4*sab(na+coset(ax,ay,az),nb+2)
                  t(coset(ax,ay,az),3) = t(coset(ax,ay+1,az),1) -&
                                         rab(2)*t(coset(ax,ay,az),1) +&
                                         f4*sab(na+coset(ax,ay,az),nb+3)
                  t(coset(ax,ay,az),4) = t(coset(ax,ay,az+1),1) -&
                                         rab(3)*t(coset(ax,ay,az),1) +&
                                         f4*sab(na+coset(ax,ay,az),nb+4)
                END DO
              END DO
            END DO

!           *** Vertical recurrence step ***

!           *** [a|T|p] = (Pi - Bi)*[a|T|s] + f2*Ni(a)*[a-1i|T|s] + ***
!           ***           f3*[a|p]                                  ***

            IF (lb_max == 1) THEN
              la_start = la_min
            ELSE
              la_start = MAX(0,la_min-1)
            END IF

            DO ax=0,la_max
              fx = f2*REAL(ax,wp)
              DO ay=0,la_max-ax
                fy = f2*REAL(ay,wp)
                az = la_max - ax - ay
                fz = f2*REAL(az,wp)
                IF (ax == 0) THEN
                  t(coset(ax,ay,az),2) = rbp(1)*t(coset(ax,ay,az),1) +&
                                         f3*sab(na+coset(ax,ay,az),nb+2)
                ELSE
                  t(coset(ax,ay,az),2) = rbp(1)*t(coset(ax,ay,az),1) +&
                                         fx*t(coset(ax-1,ay,az),1) +&
                                         f3*sab(na+coset(ax,ay,az),nb+2)
                END IF
                IF (ay == 0) THEN
                  t(coset(ax,ay,az),3) = rbp(2)*t(coset(ax,ay,az),1) +&
                                         f3*sab(na+coset(ax,ay,az),nb+3)
                ELSE
                  t(coset(ax,ay,az),3) = rbp(2)*t(coset(ax,ay,az),1) +&
                                         fy*t(coset(ax,ay-1,az),1) +&
                                         f3*sab(na+coset(ax,ay,az),nb+3)
                END IF
                IF (az == 0) THEN
                  t(coset(ax,ay,az),4) = rbp(3)*t(coset(ax,ay,az),1) +&
                                         f3*sab(na+coset(ax,ay,az),nb+4)
                ELSE
                  t(coset(ax,ay,az),4) = rbp(3)*t(coset(ax,ay,az),1) +&
                                         fz*t(coset(ax,ay,az-1),1) +&
                                         f3*sab(na+coset(ax,ay,az),nb+4)
                END IF
              END DO
            END DO

!           *** Recurrence steps: [a|T|p] -> [a|T|b] ***

            DO lb=2,lb_max

!             *** Horizontal recurrence steps ***

!             *** [a|T|b] = [a+1i|T|b-1i] - (Bi - Ai)*[a|T|b-1i] + ***
!             ***           f4*[a|b] - Ni(b-1i)*[a|b-2i]           ***

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

                    f6 = REAL(lb-1,wp)

                    t(coset(ax,ay,az),coset(0,0,lb)) =&
                      t(coset(ax,ay,az+1),coset(0,0,lb-1)) -&
                      rab(3)*t(coset(ax,ay,az),coset(0,0,lb-1)) +&
                      f4*sab(na+coset(ax,ay,az),nb+coset(0,0,lb)) -&
                      f6*sab(na+coset(ax,ay,az),nb+coset(0,0,lb-2))

!                   *** Shift of angular momentum component y from a to b ***

                    bz = lb - 1
                    t(coset(ax,ay,az),coset(0,1,bz)) =&
                      t(coset(ax,ay+1,az),coset(0,0,bz)) -&
                      rab(2)*t(coset(ax,ay,az),coset(0,0,bz)) +&
                      f4*sab(na+coset(ax,ay,az),nb+coset(0,1,bz))

                    DO by=2,lb
                      bz = lb - by
                      f6 = REAL(by-1,wp)
                      t(coset(ax,ay,az),coset(0,by,bz)) =&
                        t(coset(ax,ay+1,az),coset(0,by-1,bz)) -&
                        rab(2)*t(coset(ax,ay,az),coset(0,by-1,bz)) +&
                        f4*sab(na+coset(ax,ay,az),nb+coset(0,by,bz)) -&
                        f6*sab(na+coset(ax,ay,az),nb+coset(0,by-2,bz))
                    END DO

!                   *** Shift of angular momentum component x from a to b ***

                    DO by=0,lb-1
                      bz = lb - 1 - by
                        t(coset(ax,ay,az),coset(1,by,bz)) =&
                          t(coset(ax+1,ay,az),coset(0,by,bz)) -&
                          rab(1)*t(coset(ax,ay,az),coset(0,by,bz)) +&
                          f4*sab(na+coset(ax,ay,az),nb+coset(1,by,bz))
                    END DO

                    DO bx=2,lb
                      f6 = REAL(bx-1,wp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        t(coset(ax,ay,az),coset(bx,by,bz)) =&
                          t(coset(ax+1,ay,az),coset(bx-1,by,bz)) -&
                          rab(1)*t(coset(ax,ay,az),coset(bx-1,by,bz)) +&
                          f4*sab(na+coset(ax,ay,az),nb+coset(bx,by,bz)) -&
                          f6*sab(na+coset(ax,ay,az),nb+coset(bx-2,by,bz))
                      END DO
                    END DO

                  END DO
                END DO
              END DO

!             *** Vertical recurrence step ***

!             *** [a|T|b] = (Pi - Bi)*[a|T|b-1i] +   ***
!             ***           f2*Ni(a)*[a-1i|T|b-1i] + ***
!             ***           f2*Ni(b-1i)*[a|T|b-2i] + ***
!             ***           f3*[a|b] +               ***
!             ***           f5*Ni(b-1i)*[a|b-2i]     ***

              IF (lb == lb_max) THEN
                la_start = la_min
              ELSE
                la_start = MAX(0,la_min-1)
              END IF

              DO ax=0,la_max
                fx = f2*REAL(ax,wp)
                DO ay=0,la_max-ax
                  fy = f2*REAL(ay,wp)
                  az = la_max - ax - ay
                  fz = f2*REAL(az,wp)

!                 *** Shift of angular momentum component z from a to b ***

                  f6 = f2*REAL(lb-1,wp)
                  f7 = f5*REAL(lb-1,wp)

                  IF (az == 0) THEN
                    t(coset(ax,ay,az),coset(0,0,lb)) =&
                      rbp(3)*t(coset(ax,ay,az),coset(0,0,lb-1)) +&
                      f6*t(coset(ax,ay,az),coset(0,0,lb-2)) +&
                      f3*sab(na+coset(ax,ay,az),nb+coset(0,0,lb)) +&
                      f7*sab(na+coset(ax,ay,az),nb+coset(0,0,lb-2))
                  ELSE
                    t(coset(ax,ay,az),coset(0,0,lb)) =&
                      rbp(3)*t(coset(ax,ay,az),coset(0,0,lb-1)) +&
                      fz*t(coset(ax,ay,az-1),coset(0,0,lb-1)) +&
                      f6*t(coset(ax,ay,az),coset(0,0,lb-2)) +&
                      f3*sab(na+coset(ax,ay,az),nb+coset(0,0,lb)) +&
                      f7*sab(na+coset(ax,ay,az),nb+coset(0,0,lb-2))
                  END IF

!                 *** Shift of angular momentum component y from a to b ***

                  IF (ay == 0) THEN
                    bz = lb - 1
                    t(coset(ax,ay,az),coset(0,1,bz)) =&
                      rbp(2)*t(coset(ax,ay,az),coset(0,0,bz)) +&
                      f3*sab(na+coset(ax,ay,az),nb+coset(0,1,bz))
                    DO by=2,lb
                      bz = lb - by
                      f6 = f2*REAL(by-1,wp)
                      f7 = f5*REAL(by-1,wp)
                      t(coset(ax,ay,az),coset(0,by,bz)) =&
                        rbp(2)*t(coset(ax,ay,az),coset(0,by-1,bz)) +&
                        f6*t(coset(ax,ay,az),coset(0,by-2,bz)) +&
                        f3*sab(na+coset(ax,ay,az),nb+coset(0,by,bz)) +&
                        f7*sab(na+coset(ax,ay,az),nb+coset(0,by-2,bz))
                    END DO
                  ELSE
                    bz = lb - 1
                    t(coset(ax,ay,az),coset(0,1,bz)) =&
                      rbp(2)*t(coset(ax,ay,az),coset(0,0,bz)) +&
                      fy*t(coset(ax,ay-1,az),coset(0,0,bz)) +&
                      f3*sab(na+coset(ax,ay,az),nb+coset(0,1,bz))
                    DO by=2,lb
                      bz = lb - by
                      f6 = f2*REAL(by-1,wp)
                      f7 = f5*REAL(by-1,wp)
                      t(coset(ax,ay,az),coset(0,by,bz)) =&
                        rbp(2)*t(coset(ax,ay,az),coset(0,by-1,bz)) +&
                        fy*t(coset(ax,ay-1,az),coset(0,by-1,bz)) +&
                        f6*t(coset(ax,ay,az),coset(0,by-2,bz)) +&
                        f3*sab(na+coset(ax,ay,az),nb+coset(0,by,bz)) +&
                        f7*sab(na+coset(ax,ay,az),nb+coset(0,by-2,bz))
                    END DO
                  END IF

!                 *** Shift of angular momentum component x from a to b ***

                  IF (ax == 0) THEN
                    DO by=0,lb-1
                      bz = lb - 1 - by
                      t(coset(ax,ay,az),coset(1,by,bz)) =&
                        rbp(1)*t(coset(ax,ay,az),coset(0,by,bz)) +&
                        f3*sab(na+coset(ax,ay,az),nb+coset(1,by,bz))
                    END DO
                    DO bx=2,lb
                      f6 = f2*REAL(bx-1,wp)
                      f7 = f5*REAL(bx-1,wp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        t(coset(ax,ay,az),coset(bx,by,bz)) =&
                          rbp(1)*t(coset(ax,ay,az),coset(bx-1,by,bz)) +&
                          f6*t(coset(ax,ay,az),coset(bx-2,by,bz)) +&
                          f3*sab(na+coset(ax,ay,az),nb+coset(bx,by,bz)) +&
                          f7*sab(na+coset(ax,ay,az),nb+coset(bx-2,by,bz))
                      END DO
                    END DO
                  ELSE
                    DO by=0,lb-1
                      bz = lb - 1 - by
                      t(coset(ax,ay,az),coset(1,by,bz)) =&
                        rbp(1)*t(coset(ax,ay,az),coset(0,by,bz)) +&
                        fx*t(coset(ax-1,ay,az),coset(0,by,bz)) +&
                        f3*sab(na+coset(ax,ay,az),nb+coset(1,by,bz))
                    END DO
                    DO bx=2,lb
                      f6 = f2*REAL(bx-1,wp)
                      f7 = f5*REAL(bx-1,wp)
                      DO by=0,lb-bx
                        bz = lb - bx - by
                        t(coset(ax,ay,az),coset(bx,by,bz)) =&
                          rbp(1)*t(coset(ax,ay,az),coset(bx-1,by,bz)) +&
                          fx*t(coset(ax-1,ay,az),coset(bx-1,by,bz)) +&
                          f6*t(coset(ax,ay,az),coset(bx-2,by,bz)) +&
                          f3*sab(na+coset(ax,ay,az),nb+coset(bx,by,bz)) +&
                          f7*sab(na+coset(ax,ay,az),nb+coset(bx-2,by,bz))
                      END DO
                    END DO
                  END IF

                END DO
              END DO

            END DO

          END IF

        ELSE

          IF (lb_max > 0) THEN

!           *** Vertical recurrence steps: [s|T|s] -> [s|T|b] ***

            f5 = -zeta(ipgf)*zetp

            rbp(:) = f5*rab(:)

!           *** [s|T|p] = (Pi - Bi)*[s|T|s] + f3*[s|p] ***

            t(1,2) = rbp(1)*t(1,1) + f3*sab(na+1,nb+2)
            t(1,3) = rbp(2)*t(1,1) + f3*sab(na+1,nb+3)
            t(1,4) = rbp(3)*t(1,1) + f3*sab(na+1,nb+4)

!           *** [s|T|b] = (Pi - Bi)*[s|T|b-1i] +   ***
!           ***           f2*Ni(b-1i)*[s|T|b-2i] + ***
!           ***           f3*[s|b] +               ***
!           ***           f5*Ni(b-1i)*[s|b-2i])    ***

            DO lb=2,lb_max

!             *** Increase the angular momentum component z of function b ***

              t(1,coset(0,0,lb)) =&
                rbp(3)*t(1,coset(0,0,lb-1)) +&
                f2*REAL(lb-1,wp)*t(1,coset(0,0,lb-2)) +&
                f3*sab(na+1,nb+coset(0,0,lb)) +&
                f5*REAL(lb-1,wp)*sab(na+1,nb+coset(0,0,lb-2))

!             *** Increase the angular momentum component y of function b ***

              bz = lb - 1
              t(1,coset(0,1,bz)) =&
                rbp(2)*t(1,coset(0,0,bz)) +&
                f3*sab(na+1,nb+coset(0,1,bz))

              DO by=2,lb
                bz = lb - by
                t(1,coset(0,by,bz)) =&
                  rbp(2)*t(1,coset(0,by-1,bz)) +&
                  f2*REAL(by-1,wp)*t(1,coset(0,by-2,bz)) +&
                  f3*sab(na+1,nb+coset(0,by,bz)) +&
                  f5*REAL(by-1,wp)*sab(na+1,nb+coset(0,by-2,bz))
              END DO

!             *** Increase the angular momentum component x of function b ***

              DO by=0,lb-1
                bz = lb - 1 - by
                t(1,coset(1,by,bz)) =&
                  rbp(1)*t(1,coset(0,by,bz)) +&
                  f3*sab(na+1,nb+coset(1,by,bz))
              END DO

              DO bx=2,lb
                f6 = f2*REAL(bx-1,wp)
                f7 = f5*REAL(bx-1,wp)
                DO by=0,lb-bx
                  bz = lb - bx - by
                  t(1,coset(bx,by,bz)) =&
                    rbp(1)*t(1,coset(bx-1,by,bz)) +&
                    f6*t(1,coset(bx-2,by,bz)) +&
                    f3*sab(na+1,nb+coset(bx,by,bz)) +&
                    f7*sab(na+1,nb+coset(bx-2,by,bz))
                END DO
              END DO

            END DO

          END IF

        END IF

        DO j=ncoset(lb_min-1)+1,ncoset(lb_max)
          DO i=ncoset(la_min-1)+1,ncoset(la_max)
            tab(na+i,nb+j) = t(i,j)
          END DO
        END DO

        nb = nb + ncoset(lb_max)

      END DO

      na = na + ncoset(la_max)

    END DO

  END SUBROUTINE kinetic

! *****************************************************************************

END MODULE ai_kinetic

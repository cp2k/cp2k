!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/ai_nuclear [1.0] *
!!
!!   NAME
!!     ai_nuclear
!!
!!   FUNCTION
!!     Calculation of the two-center nuclear attraction integrals over
!!     Cartesian Gaussian-type functions.
!!
!!     V(nuc) = <a|-Z/R(c)|b>
!!
!!   AUTHOR
!!     Matthias Krack (26.09.2000)
!!
!!   LITERATURE
!!     S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE ai_nuclear

! *****************************************************************************

! ax,ay,az  : Angular momentum index numbers of orbital a.
! bx,by,bz  : Angular momentum index numbers of orbital b.
! coset     : Cartesian orbital set pointer.
! l{a,b}    : Angular momentum quantum number of shell a or b.
! l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
! l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
! ncoset    : Number of orbitals in a Cartesian orbital set.
! npgf{a,b} : Degree of contraction of shell a or b.
! rab       : Distance vector between the atomic centers a and b.
! rab2      : Square of the distance between the atomic centers a and b.
! rac       : Distance vector between the atomic centers a and c.
! rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
! zet{a,b}  : Exponents of the Gaussian-type functions a or b.
! zetp      : Reciprocal of the sum of the exponents of orbital a and b.

! *****************************************************************************

  USE kinds, ONLY: dp

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: nuclear

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE nuclear(la_max,zeta,rpgfa,la_min,&
                     lb_max,zetb,rpgfb,lb_min,&
                     zc,rab,rab2,rac,vab)

!   Purpose: Calculation of the primitive two-center nuclear attraction
!            integrals over Cartesian Gaussian-type functions.

!   History: - Creation (26.09.2000, Matthias Krack)

!   ***************************************************************************

    USE gamma,            ONLY: fgamma
    USE mathconstants,    ONLY: pi
    USE orbital_pointers, ONLY: coset,ncoset

    REAL(KIND = dp), INTENT(IN)                    :: rab2,zc
    INTEGER, INTENT(IN)                     :: la_max,la_min,lb_max,lb_min
    REAL(KIND = dp), DIMENSION(3), INTENT(IN)      :: rab,rac
    REAL(KIND = dp), DIMENSION(:), INTENT(IN)      :: rpgfa,rpgfb,zeta,zetb
    REAL(KIND = dp), DIMENSION(:,:), INTENT(INOUT) :: vab

!   *** Local variables ***

    REAL(KIND = dp) :: f0,f1,f2,f3,fx,fy,fz,t,zetp
    INTEGER  :: ax,ay,az,bx,by,bz,i,ipgf,j,jpgf,la,la_start,lb,n,na,nb,nmax

    REAL(KIND = dp), DIMENSION(3)               :: rap,rbp,rcp
    REAL(KIND = dp), DIMENSION(0:la_max+lb_max) :: f
    REAL(KIND = dp), DIMENSION(ncoset(la_max),&
                        ncoset(lb_max),&
                        la_max+lb_max+1) :: v

!   ---------------------------------------------------------------------------

    nmax = la_max + lb_max + 1

!   *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0

    DO ipgf=1,SIZE(zeta)

      nb = 0

      DO jpgf=1,SIZE(zetb)

!       *** Calculate some prefactors ***

        zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))

        f0 = 2.0_dp*pi*zetp
        f1 = zetb(jpgf)*zetp
        f2 = 0.5_dp*zetp

        f0 = f0*EXP(-zeta(ipgf)*f1*rab2)

        rap(:) = f1*rab(:)
        rcp(:) = rap(:) - rac(:)

!       *** Calculate the incomplete Gamma function ***

        t = (rcp(1)*rcp(1) + rcp(2)*rcp(2) + rcp(3)*rcp(3))/zetp

        CALL fgamma(nmax-1,t,f)

!       *** Calculate the basic nuclear attraction integrals [s|A(0)|s]{n} ***

        DO n=1,nmax
          v(1,1,n) = f0*f(n-1)
        END DO

!       *** Recurrence steps: [s|A(0)|s] -> [a|A(0)|b] ***

        IF (la_max > 0) THEN

!         *** Vertical recurrence steps: [s|A(0)|s] -> [a|A(0)|s] ***

!         *** [p|A(0)|s]{n} = (Pi - Ai)*[s|A(0)|s]{n} -              ***
!         ***                 (Pi - Ci)*[s|A(0)|s]{n+1}  (i = x,y,z) ***

          DO n=1,nmax-1
            v(2,1,n) = rap(1)*v(1,1,n) - rcp(1)*v(1,1,n+1)
            v(3,1,n) = rap(2)*v(1,1,n) - rcp(2)*v(1,1,n+1)
            v(4,1,n) = rap(3)*v(1,1,n) - rcp(3)*v(1,1,n+1)
          END DO

!         *** [a|A(0)|s]{n} = (Pi - Ai)*[a-1i|A(0)|s]{n} -    ***
!         ***                 (Pi - Ci)*[a-1i|A(0)|s]{n+1} +  ***
!         ***                 f2*Ni(a-1i)*([a-2i|A(0)|s]{n} - ***
!         ***                              [a-2i|A(0)|s]{n+1} ***

          DO la=2,la_max

            DO n=1,nmax-la

!             *** Increase the angular momentum component z of function a ***

              v(coset(0,0,la),1,n) =&
                rap(3)*v(coset(0,0,la-1),1,n) -&
                rcp(3)*v(coset(0,0,la-1),1,n+1) +&
                f2*REAL(la-1,dp)*(v(coset(0,0,la-2),1,n) -&
                                  v(coset(0,0,la-2),1,n+1))

!             *** Increase the angular momentum component y of function a ***

              az = la - 1
              v(coset(0,1,az),1,n) =&
                rap(2)*v(coset(0,0,az),1,n) -&
                rcp(2)*v(coset(0,0,az),1,n+1)

              DO ay=2,la
                az = la - ay
                v(coset(0,ay,az),1,n) =&
                  rap(2)*v(coset(0,ay-1,az),1,n) -&
                  rcp(2)*v(coset(0,ay-1,az),1,n+1) +&
                  f2*REAL(ay-1,dp)*(v(coset(0,ay-2,az),1,n) -&
                                    v(coset(0,ay-2,az),1,n+1))
              END DO

!             *** Increase the angular momentum component x of function a ***

              DO ay=0,la-1
                az = la - 1 - ay
                v(coset(1,ay,az),1,n) =&
                  rap(1)*v(coset(0,ay,az),1,n) -&
                  rcp(1)*v(coset(0,ay,az),1,n+1)
              END DO

              DO ax=2,la
                f3 = f2*REAL(ax-1,dp)
                DO ay=0,la-ax
                  az = la - ax - ay
                  v(coset(ax,ay,az),1,n) =&
                    rap(1)*v(coset(ax-1,ay,az),1,n) -&
                    rcp(1)*v(coset(ax-1,ay,az),1,n+1) +&
                    f3*(v(coset(ax-2,ay,az),1,n) -&
                        v(coset(ax-2,ay,az),1,n+1))
                END DO
              END DO

            END DO

          END DO

!         *** Recurrence steps: [a|A(0)|s] -> [a|A(0)|b] ***

          IF (lb_max > 0) THEN

!           *** Horizontal recurrence steps ***

            rbp(:) = rap(:) - rab(:)

!           *** [a||A(0)|p]{n} = [a+1i|A(0)|s]{n} -      ***
!           ***                  (Bi - Ai)*[a|A(0)|s]{n} ***

            IF (lb_max == 1) THEN
              la_start = la_min
            ELSE
              la_start = MAX(0,la_min-1)
            END IF

            DO la=la_start,la_max-1
              DO n=1,nmax-la-1
                DO ax=0,la
                  DO ay=0,la-ax
                    az = la - ax - ay
                    v(coset(ax,ay,az),2,n) =&
                      v(coset(ax+1,ay,az),1,n) -&
                      rab(1)*v(coset(ax,ay,az),1,n)
                    v(coset(ax,ay,az),3,n) =&
                      v(coset(ax,ay+1,az),1,n) -&
                      rab(2)*v(coset(ax,ay,az),1,n)
                    v(coset(ax,ay,az),4,n) =&
                      v(coset(ax,ay,az+1),1,n) -&
                      rab(3)*v(coset(ax,ay,az),1,n)
                  END DO
                END DO
              END DO
            END DO

!           *** Vertical recurrence step ***

!           *** [a|A(0)|p]{n} = (Pi - Bi)*[a|A(0)|s]{n} -     ***
!           ***                 (Pi - Ci)*[a|A(0)|s]{n+1} +   ***
!           ***                 f2*Ni(a)*([a-1i|A(0)|s]{n} -  ***
!           ***                           [a-1i|A(0)|s]{n+1}) ***

            DO n=1,nmax-la_max-1
              DO ax=0,la_max
                fx = f2*REAL(ax,dp)
                DO ay=0,la_max-ax
                  fy = f2*REAL(ay,dp)
                  az = la_max - ax - ay
                  fz = f2*REAL(az,dp)

                  IF (ax == 0) THEN
                    v(coset(ax,ay,az),2,n) =&
                      rbp(1)*v(coset(ax,ay,az),1,n) -&
                      rcp(1)*v(coset(ax,ay,az),1,n+1)
                  ELSE
                    v(coset(ax,ay,az),2,n) =&
                      rbp(1)*v(coset(ax,ay,az),1,n) -&
                      rcp(1)*v(coset(ax,ay,az),1,n+1) +&
                      fx*(v(coset(ax-1,ay,az),1,n) -&
                          v(coset(ax-1,ay,az),1,n+1))
                  END IF

                  IF (ay == 0) THEN
                    v(coset(ax,ay,az),3,n) =&
                      rbp(2)*v(coset(ax,ay,az),1,n) -&
                      rcp(2)*v(coset(ax,ay,az),1,n+1)
                  ELSE
                    v(coset(ax,ay,az),3,n) =&
                      rbp(2)*v(coset(ax,ay,az),1,n) -&
                      rcp(2)*v(coset(ax,ay,az),1,n+1) +&
                      fy*(v(coset(ax,ay-1,az),1,n) -&
                          v(coset(ax,ay-1,az),1,n+1))
                  END IF

                  IF (az == 0) THEN
                    v(coset(ax,ay,az),4,n) =&
                      rbp(3)*v(coset(ax,ay,az),1,n) -&
                      rcp(3)*v(coset(ax,ay,az),1,n+1)
                  ELSE
                    v(coset(ax,ay,az),4,n) =&
                      rbp(3)*v(coset(ax,ay,az),1,n) -&
                      rcp(3)*v(coset(ax,ay,az),1,n+1) +&
                      fz*(v(coset(ax,ay,az-1),1,n) -&
                          v(coset(ax,ay,az-1),1,n+1))
                  END IF

                END DO
              END DO
            END DO

!           *** Recurrence steps: [a|A(0)|p] -> [a|A(0)|b] ***

            DO lb=2,lb_max

!             *** Horizontal recurrence steps ***

!             *** [a||A(0)|b]{n} = [a+1i|A(0)|b-1i]{n} -      ***
!             ***                  (Bi - Ai)*[a|A(0)|b-1i]{n} ***

              IF (lb == lb_max) THEN
                la_start = la_min
              ELSE
                la_start = MAX(0,la_min-1)
              END IF

              DO la=la_start,la_max-1
                DO n=1,nmax-la-lb
                  DO ax=0,la
                    DO ay=0,la-ax
                      az = la - ax - ay

!                     *** Shift of angular momentum component z from a to b ***

                      v(coset(ax,ay,az),coset(0,0,lb),n) =&
                        v(coset(ax,ay,az+1),coset(0,0,lb-1),n) -&
                        rab(3)*v(coset(ax,ay,az),coset(0,0,lb-1),n)

!                     *** Shift of angular momentum component y from a to b ***

                      DO by=1,lb
                        bz = lb - by
                        v(coset(ax,ay,az),coset(0,by,bz),n) =&
                          v(coset(ax,ay+1,az),coset(0,by-1,bz),n) -&
                          rab(2)*v(coset(ax,ay,az),coset(0,by-1,bz),n)
                      END DO

!                     *** Shift of angular momentum component x from a to b ***

                      DO bx=1,lb
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          v(coset(ax,ay,az),coset(bx,by,bz),n) =&
                            v(coset(ax+1,ay,az),coset(bx-1,by,bz),n) -&
                            rab(1)*v(coset(ax,ay,az),coset(bx-1,by,bz),n)
                        END DO
                      END DO

                    END DO
                  END DO
                END DO
              END DO

!             *** Vertical recurrence step ***

!             *** [a|A(0)|b]{n} = (Pi - Bi)*[a|A(0)|b-1i]{n} -       ***
!             ***                 (Pi - Ci)*[a|A(0)|b-1i]{n+1} +     ***
!             ***                 f2*Ni(a)*([a-1i|A(0)|b-1i]{n} -    ***
!             ***                           [a-1i|A(0)|b-1i]{n+1}) + ***
!             ***                 f2*Ni(b-1i)*([a|A(0)|b-2i]{n} -    ***
!             ***                              [a|A(0)|b-2i]{n+1})   ***

              DO n=1,nmax-la_max-lb
                DO ax=0,la_max
                  fx = f2*REAL(ax,dp)
                  DO ay=0,la_max-ax
                    fy = f2*REAL(ay,dp)
                    az = la_max - ax - ay
                    fz = f2*REAL(az,dp)

!                   *** Shift of angular momentum component z from a to b ***

                    f3 = f2*REAL(lb-1,dp)

                    IF (az == 0) THEN
                      v(coset(ax,ay,az),coset(0,0,lb),n) =&
                        rbp(3)*v(coset(ax,ay,az),coset(0,0,lb-1),n) -&
                        rcp(3)*v(coset(ax,ay,az),coset(0,0,lb-1),n+1) +&
                        f3*(v(coset(ax,ay,az),coset(0,0,lb-2),n) -&
                            v(coset(ax,ay,az),coset(0,0,lb-2),n+1))
                    ELSE
                      v(coset(ax,ay,az),coset(0,0,lb),n) =&
                        rbp(3)*v(coset(ax,ay,az),coset(0,0,lb-1),n) -&
                        rcp(3)*v(coset(ax,ay,az),coset(0,0,lb-1),n+1) +&
                        fz*(v(coset(ax,ay,az-1),coset(0,0,lb-1),n) -&
                            v(coset(ax,ay,az-1),coset(0,0,lb-1),n+1)) +&
                        f3*(v(coset(ax,ay,az),coset(0,0,lb-2),n) -&
                            v(coset(ax,ay,az),coset(0,0,lb-2),n+1))
                    END IF

!                   *** Shift of angular momentum component y from a to b ***

                    IF (ay == 0) THEN
                      bz = lb - 1
                      v(coset(ax,ay,az),coset(0,1,bz),n) =&
                        rbp(2)*v(coset(ax,ay,az),coset(0,0,bz),n) -&
                        rcp(2)*v(coset(ax,ay,az),coset(0,0,bz),n+1)
                      DO by=2,lb
                        bz = lb - by
                        f3 = f2*REAL(by-1,dp)
                        v(coset(ax,ay,az),coset(0,by,bz),n) =&
                          rbp(2)*v(coset(ax,ay,az),coset(0,by-1,bz),n) -&
                          rcp(2)*v(coset(ax,ay,az),coset(0,by-1,bz),n+1) +&
                          f3*(v(coset(ax,ay,az),coset(0,by-2,bz),n) -&
                              v(coset(ax,ay,az),coset(0,by-2,bz),n+1))
                      END DO
                    ELSE
                      bz = lb - 1
                      v(coset(ax,ay,az),coset(0,1,bz),n) =&
                        rbp(2)*v(coset(ax,ay,az),coset(0,0,bz),n) -&
                        rcp(2)*v(coset(ax,ay,az),coset(0,0,bz),n+1) +&
                        fy*(v(coset(ax,ay-1,az),coset(0,0,bz),n) -&
                            v(coset(ax,ay-1,az),coset(0,0,bz),n+1))
                      DO by=2,lb
                        bz = lb - by
                        f3 = f2*REAL(by-1,dp)
                        v(coset(ax,ay,az),coset(0,by,bz),n) =&
                          rbp(2)*v(coset(ax,ay,az),coset(0,by-1,bz),n) -&
                          rcp(2)*v(coset(ax,ay,az),coset(0,by-1,bz),n+1) +&
                          fy*(v(coset(ax,ay-1,az),coset(0,by-1,bz),n) -&
                              v(coset(ax,ay-1,az),coset(0,by-1,bz),n+1)) +&
                          f3*(v(coset(ax,ay,az),coset(0,by-2,bz),n) -&
                              v(coset(ax,ay,az),coset(0,by-2,bz),n+1))
                      END DO
                    END IF

!                   *** Shift of angular momentum component x from a to b ***

                    IF (ax == 0) THEN
                      DO by=0,lb-1
                        bz = lb - 1 - by
                        v(coset(ax,ay,az),coset(1,by,bz),n) =&
                          rbp(1)*v(coset(ax,ay,az),coset(0,by,bz),n) -&
                          rcp(1)*v(coset(ax,ay,az),coset(0,by,bz),n+1)
                      END DO
                      DO bx=2,lb
                        f3 = f2*REAL(bx-1,dp)
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          v(coset(ax,ay,az),coset(bx,by,bz),n) =&
                            rbp(1)*v(coset(ax,ay,az),coset(bx-1,by,bz),n) -&
                            rcp(1)*v(coset(ax,ay,az),coset(bx-1,by,bz),n+1) +&
                            f3*(v(coset(ax,ay,az),coset(bx-2,by,bz),n) -&
                                v(coset(ax,ay,az),coset(bx-2,by,bz),n+1))
                        END DO
                      END DO
                    ELSE
                      DO by=0,lb-1
                        bz = lb - 1 - by
                        v(coset(ax,ay,az),coset(1,by,bz),n) =&
                          rbp(1)*v(coset(ax,ay,az),coset(0,by,bz),n) -&
                          rcp(1)*v(coset(ax,ay,az),coset(0,by,bz),n+1) +&
                          fx*(v(coset(ax-1,ay,az),coset(0,by,bz),n) -&
                              v(coset(ax-1,ay,az),coset(0,by,bz),n+1))
                      END DO
                      DO bx=2,lb
                        f3 = f2*REAL(bx-1,dp)
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          v(coset(ax,ay,az),coset(bx,by,bz),n) =&
                            rbp(1)*v(coset(ax,ay,az),coset(bx-1,by,bz),n) -&
                            rcp(1)*v(coset(ax,ay,az),coset(bx-1,by,bz),n+1) +&
                            fx*(v(coset(ax-1,ay,az),coset(bx-1,by,bz),n) -&
                                v(coset(ax-1,ay,az),coset(bx-1,by,bz),n+1)) +&
                            f3*(v(coset(ax,ay,az),coset(bx-2,by,bz),n) -&
                                v(coset(ax,ay,az),coset(bx-2,by,bz),n+1))
                        END DO
                      END DO
                    END IF

                  END DO
                END DO
              END DO

            END DO

          END IF

        ELSE

          IF (lb_max > 0) THEN

!           *** Vertical recurrence steps: [s|A(0)|s] -> [s|A(0)|b] ***

            rbp(:) = rap(:) - rab(:)

!           *** [s|A(0)|p]{n} = (Pi - Bi)*[s|A(0)|s]{n} - ***
!           ***                 (Pi - Ci)*[s|A(0)|s]{n+1} ***

            DO n=1,nmax-1
              v(1,2,n) = rbp(1)*v(1,1,n) - rcp(1)*v(1,1,n+1)
              v(1,3,n) = rbp(2)*v(1,1,n) - rcp(2)*v(1,1,n+1)
              v(1,4,n) = rbp(3)*v(1,1,n) - rcp(3)*v(1,1,n+1)
            END DO

!           *** [s|A(0)|b]{n} = (Pi - Bi)*[s|A(0)|b-1i]{n} -    ***
!           ***                 (Pi - Ci)*[s|A(0)|b-1i]{n+1} +  ***
!           ***                 f2*Ni(b-1i)*([s|A(0)|b-2i]{n} - ***
!           ***                              [s|A(0)|b-2i]{n+1} ***

            DO lb=2,lb_max

              DO n=1,nmax-lb

!               *** Increase the angular momentum component z of function b ***

                v(1,coset(0,0,lb),n) =&
                  rbp(3)*v(1,coset(0,0,lb-1),n) -&
                  rcp(3)*v(1,coset(0,0,lb-1),n+1) +&
                  f2*REAL(lb-1,dp)*(v(1,coset(0,0,lb-2),n) -&
                                    v(1,coset(0,0,lb-2),n+1))

!               *** Increase the angular momentum component y of function b ***

                bz = lb - 1
                v(1,coset(0,1,bz),n) =&
                  rbp(2)*v(1,coset(0,0,bz),n) -&
                  rcp(2)*v(1,coset(0,0,bz),n+1)

                DO by=2,lb
                  bz = lb - by
                  v(1,coset(0,by,bz),n) =&
                    rbp(2)*v(1,coset(0,by-1,bz),n) -&
                    rcp(2)*v(1,coset(0,by-1,bz),n+1) +&
                    f2*REAL(by-1,dp)*(v(1,coset(0,by-2,bz),n) -&
                                      v(1,coset(0,by-2,bz),n+1))
                END DO

!               *** Increase the angular momentum component x of function b ***

                DO by=0,lb-1
                  bz = lb - 1 - by
                  v(1,coset(1,by,bz),n) =&
                    rbp(1)*v(1,coset(0,by,bz),n) -&
                    rcp(1)*v(1,coset(0,by,bz),n+1)
                END DO

                DO bx=2,lb
                  f3 = f2*REAL(bx-1,dp)
                  DO by=0,lb-bx
                    bz = lb - bx - by
                    v(1,coset(bx,by,bz),n) =&
                      rbp(1)*v(1,coset(bx-1,by,bz),n) -&
                      rcp(1)*v(1,coset(bx-1,by,bz),n+1) +&
                      f3*(v(1,coset(bx-2,by,bz),n) -&
                          v(1,coset(bx-2,by,bz),n+1))
                  END DO
                END DO

              END DO

            END DO

          END IF

        END IF

!       *** Add the contribution of the current pair ***
!       *** of primitive Gaussian-type functions     ***

        DO j=ncoset(lb_min-1)+1,ncoset(lb_max)
          DO i=ncoset(la_min-1)+1,ncoset(la_max)
            vab(na+i,nb+j) = vab(na+i,nb+j) - zc*v(i,j,1)
          END DO
        END DO

        nb = nb + ncoset(lb_max)

      END DO

      na = na + ncoset(la_max)

    END DO

  END SUBROUTINE nuclear

! *****************************************************************************

END MODULE ai_nuclear

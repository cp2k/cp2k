!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of Coulomb integrals over Cartesian Gaussian-type functions
!>      (electron repulsion integrals, ERIs).
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      none
!> \par Parameters
!>       - ax,ay,az    : Angular momentum index numbers of orbital a.
!>       - bx,by,bz    : Angular momentum index numbers of orbital b.
!>       - cx,cy,cz    : Angular momentum index numbers of orbital c.
!>       - coset       : Cartesian orbital set pointer.
!>       - dab         : Distance between the atomic centers a and b.
!>       - dac         : Distance between the atomic centers a and c.
!>       - dbc         : Distance between the atomic centers b and c.
!>       - gccc        : Prefactor of the primitive Gaussian function c.
!>       - l{a,b,c}    : Angular momentum quantum number of shell a, b or c.
!>       - l{a,b,c}_max: Maximum angular momentum quantum number of shell a, b or c.
!>       - l{a,b,c}_min: Minimum angular momentum quantum number of shell a, b or c.
!>       - ncoset      : Number of orbitals in a Cartesian orbital set.
!>       - npgf{a,b}   : Degree of contraction of shell a or b.
!>       - rab         : Distance vector between the atomic centers a and b.
!>       - rab2        : Square of the distance between the atomic centers a and b.
!>       - rac         : Distance vector between the atomic centers a and c.
!>       - rac2        : Square of the distance between the atomic centers a and c.
!>       - rbc         : Distance vector between the atomic centers b and c.
!>       - rbc2        : Square of the distance between the atomic centers b and c.
!>       - rpgf{a,b,c} : Radius of the primitive Gaussian-type function a, b or c.
!>       - zet{a,b,c}  : Exponents of the Gaussian-type functions a, b or c.
!>       - zetp        : Reciprocal of the sum of the exponents of orbital a and b.
!>       - zetw        : Reciprocal of the sum of the exponents of orbital a, b and c.
!> \author Matthias Krack (22.08.2000)
! *****************************************************************************
MODULE ai_coulomb
  USE f77_blas
  USE gamma,                           ONLY: fgamma=>fgamma_0
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
#include "../common/cp_common_uses.f90"

   IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_coulomb'
   PRIVATE

 ! *** Public subroutines ***

  PUBLIC :: coulomb2, coulomb2_new

 CONTAINS

! *****************************************************************************
!> \brief   Calculation of the primitive two-center Coulomb integrals over
!>          Cartesian Gaussian-type functions.
!> \author  Matthias Krack
!> \date    05.12.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE coulomb2(la_max,npgfa,zeta,rpgfa,la_min,lc_max,npgfc,zetc,rpgfc,lc_min,&
       rac,rac2,vac,v,f,maxder,vac_plus)
    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: la_min, lc_max, npgfc
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetc, rpgfc
    INTEGER, INTENT(IN)                      :: lc_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac
    REAL(KIND=dp), INTENT(IN)                :: rac2
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: vac
    REAL(KIND=dp), DIMENSION(:, :, :)        :: v
    REAL(KIND=dp), DIMENSION(0:)             :: f
    INTEGER, INTENT(IN), OPTIONAL            :: maxder
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL :: vac_plus

    CHARACTER(len=*), PARAMETER :: routineN = 'coulomb2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ax, ay, az, coc, cocx, cocy, &
                                                cocz, cx, cy, cz, i, ipgf, j, &
                                                jpgf, la, lc, maxder_local, &
                                                n, na, nap, nc, nmax
    REAL(KIND=dp)                            :: dac, f0, f1, f2, f3, f4, f5, &
                                                f6, fcx, fcy, fcz, rho, t, &
                                                zetp, zetq, zetw
    REAL(KIND=dp), DIMENSION(3)              :: raw, rcw

    v = 0.0_dp

    maxder_local = 0
    IF(PRESENT(maxder)) THEN
       maxder_local = maxder
       vac_plus = 0.0_dp
    END IF

    nmax = la_max + lc_max + 1

    ! *** Calculate the distance of the centers a and c ***

    dac = SQRT(rac2)

    ! *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0
    nap = 0

    DO ipgf=1,npgfa

       nc = 0

       DO jpgf=1,npgfc

          ! *** Screening ***

          IF (rpgfa(ipgf) + rpgfc(jpgf) < dac) THEN
             DO j=nc+ncoset(lc_min-1)+1,nc+ncoset(lc_max)
                DO i=na+ncoset(la_min-1)+1,na+ncoset(la_max-maxder_local)
                   vac(i,j) = 0.0_dp
                END DO
             END DO
             nc = nc + ncoset(lc_max)
             CYCLE
          END IF

          ! *** Calculate some prefactors ***

          zetp = 1.0_dp/zeta(ipgf)
          zetq = 1.0_dp/zetc(jpgf)
          zetw = 1.0_dp/(zeta(ipgf) + zetc(jpgf))

          rho = zeta(ipgf)*zetc(jpgf)*zetw

          f0 = 2.0_dp*SQRT(pi**5*zetw)*zetp*zetq

          ! *** Calculate the incomplete Gamma function ***

          t = rho*rac2

          CALL fgamma(nmax-1,t,f)

          ! *** Calculate the basic two-center Coulomb integrals [s||s]{n} ***

          DO n=1,nmax
             v(1,1,n) = f0*f(n-1)
          END DO

          ! *** Vertical recurrence steps: [s||s] -> [s||c] ***

          IF (lc_max > 0) THEN

             f1 = 0.5_dp*zetq
             f2 = -rho*zetq

             rcw(:) = -zeta(ipgf)*zetw*rac(:)

             ! *** [s||p]{n} = (Wi - Ci)*[s||s]{n+1}  (i = x,y,z) ***

             DO n=1,nmax-1
                v(1,2,n) = rcw(1)*v(1,1,n+1)
                v(1,3,n) = rcw(2)*v(1,1,n+1)
                v(1,4,n) = rcw(3)*v(1,1,n+1)
             END DO

             ! ** [s||c]{n} = (Wi - Ci)*[s||c-1i]{n+1} +     ***
             ! **             f1*Ni(c-1i)*(   [s||c-2i]{n} + ***
             ! **                          f2*[s||c-2i]{n+1} ***

             DO lc=2,lc_max

                DO n=1,nmax-lc

                   ! **** Increase the angular momentum component z of c ***

                   v(1,coset(0,0,lc),n) =&
                        rcw(3)*v(1,coset(0,0,lc-1),n+1) +&
                        f1*REAL(lc-1,dp)*(   v(1,coset(0,0,lc-2),n) +&
                        f2*v(1,coset(0,0,lc-2),n+1))

                   ! *** Increase the angular momentum component y of c ***

                   cz = lc - 1
                   v(1,coset(0,1,cz),n) = rcw(2)*v(1,coset(0,0,cz),n+1)

                   DO cy=2,lc
                      cz = lc - cy
                      v(1,coset(0,cy,cz),n) =&
                           rcw(2)*v(1,coset(0,cy-1,cz),n+1) +&
                           f1*REAL(cy-1,dp)*(   v(1,coset(0,cy-2,cz),n) +&
                           f2*v(1,coset(0,cy-2,cz),n+1))
                   END DO

                   ! *** Increase the angular momentum component x of c ***

                   DO cy=0,lc-1
                      cz = lc - 1 - cy
                      v(1,coset(1,cy,cz),n) = rcw(1)*v(1,coset(0,cy,cz),n+1)
                   END DO

                   DO cx=2,lc
                      f6 = f1*REAL(cx-1,dp)
                      DO cy=0,lc-cx
                         cz = lc - cx - cy
                         v(1,coset(cx,cy,cz),n) =&
                              rcw(1)*v(1,coset(cx-1,cy,cz),n+1) +&
                              f6*(   v(1,coset(cx-2,cy,cz),n) +&
                              f2*v(1,coset(cx-2,cy,cz),n+1))
                      END DO
                   END DO

                END DO

             END DO

          END IF

          ! *** Vertical recurrence steps: [s||c] -> [a||c] ***

          IF (la_max > 0) THEN

             f3 = 0.5_dp*zetp
             f4 = -rho*zetp
             f5 = 0.5_dp*zetw

             raw(:) = zetc(jpgf)*zetw*rac(:)

             ! *** [p||s]{n} = (Wi - Ai)*[s||s]{n+1}  (i = x,y,z) ***

             DO n=1,nmax-1
                v(2,1,n) = raw(1)*v(1,1,n+1)
                v(3,1,n) = raw(2)*v(1,1,n+1)
                v(4,1,n) = raw(3)*v(1,1,n+1)
             END DO

             ! *** [a||s]{n} = (Wi - Ai)*[a-1i||s]{n+1} +      ***
             ! ***             f3*Ni(a-1i)*(   [a-2i||s]{n} +  ***
             ! ***                          f4*[a-2i||s]{n+1}) ***

             DO la=2,la_max

                DO n=1,nmax-la

                   ! *** Increase the angular momentum component z of a ***

                   v(coset(0,0,la),1,n) =&
                        raw(3)*v(coset(0,0,la-1),1,n+1) +&
                        f3*REAL(la-1,dp)*(   v(coset(0,0,la-2),1,n) +&
                        f4*v(coset(0,0,la-2),1,n+1))

                   ! *** Increase the angular momentum component y of a ***

                   az = la - 1
                   v(coset(0,1,az),1,n) = raw(2)*v(coset(0,0,az),1,n+1)

                   DO ay=2,la
                      az = la - ay
                      v(coset(0,ay,az),1,n) =&
                           raw(2)*v(coset(0,ay-1,az),1,n+1) +&
                           f3*REAL(ay-1,dp)*(   v(coset(0,ay-2,az),1,n) +&
                           f4*v(coset(0,ay-2,az),1,n+1))
                   END DO

                   ! *** Increase the angular momentum component x of a ***

                   DO ay=0,la-1
                      az = la - 1 - ay
                      v(coset(1,ay,az),1,n) = raw(1)*v(coset(0,ay,az),1,n+1)
                   END DO

                   DO ax=2,la
                      f6 = f3*REAL(ax-1,dp)
                      DO ay=0,la-ax
                         az = la - ax - ay
                         v(coset(ax,ay,az),1,n) =&
                              raw(1)*v(coset(ax-1,ay,az),1,n+1) +&
                              f6*(   v(coset(ax-2,ay,az),1,n) +&
                              f4*v(coset(ax-2,ay,az),1,n+1))
                      END DO
                   END DO

                END DO

             END DO

             DO lc=1,lc_max

                DO cx=0,lc
                   DO cy=0,lc-cx
                      cz = lc - cx - cy

                      coc = coset(cx,cy,cz)
                      cocx = coset(MAX(0,cx-1),cy,cz)
                      cocy = coset(cx,MAX(0,cy-1),cz)
                      cocz = coset(cx,cy,MAX(0,cz-1))

                      fcx = f5*REAL(cx,dp)
                      fcy = f5*REAL(cy,dp)
                      fcz = f5*REAL(cz,dp)

                      ! *** [p||c]{n} = (Wi - Ai)*[s||c]{n+1} + ***
                      ! ***             f5*Ni(c)*[s||c-1i]{n+1} ***

                      DO n=1,nmax-1-lc
                         v(2,coc,n) = raw(1)*v(1,coc,n+1) + fcx*v(1,cocx,n+1)
                         v(3,coc,n) = raw(2)*v(1,coc,n+1) + fcy*v(1,cocy,n+1)
                         v(4,coc,n) = raw(3)*v(1,coc,n+1) + fcz*v(1,cocz,n+1)
                      END DO

                      ! *** [a||c]{n} = (Wi - Ai)*[a-1i||c]{n+1} +        ***
                      ! ***             f3*Ni(a-1i)*(   [a-2i||c]{n} +    ***
                      ! ***                          f4*[a-2i||c]{n+1}) + ***
                      ! ***             f5*Ni(c)*[a-1i||c-1i]{n+1}        ***

                      DO la=2,la_max

                         DO n=1,nmax-la-lc

                            ! *** Increase the angular momentum component z of a ***

                            v(coset(0,0,la),coc,n) =&
                                 raw(3)*v(coset(0,0,la-1),coc,n+1) +&
                                 f3*REAL(la-1,dp)*(   v(coset(0,0,la-2),coc,n) +&
                                 f4*v(coset(0,0,la-2),coc,n+1)) +&
                                 fcz*v(coset(0,0,la-1),cocz,n+1)

                            ! *** Increase the angular momentum component y of a ***

                            az = la - 1
                            v(coset(0,1,az),coc,n) =&
                                 raw(2)*v(coset(0,0,az),coc,n+1) +&
                                 fcy*v(coset(0,0,az),cocy,n+1)

                            DO ay=2,la
                               az = la - ay
                               v(coset(0,ay,az),coc,n) =&
                                    raw(2)*v(coset(0,ay-1,az),coc,n+1) +&
                                    f3*REAL(ay-1,dp)*(   v(coset(0,ay-2,az),coc,n) +&
                                    f4*v(coset(0,ay-2,az),coc,n+1)) +&
                                    fcy*v(coset(0,ay-1,az),cocy,n+1)
                            END DO

                            ! *** Increase the angular momentum component x of a ***

                            DO ay=0,la-1
                               az = la - 1 - ay
                               v(coset(1,ay,az),coc,n) =&
                                    raw(1)*v(coset(0,ay,az),coc,n+1) +&
                                    fcx*v(coset(0,ay,az),cocx,n+1)
                            END DO

                            DO ax=2,la
                               f6 = f3*REAL(ax-1,dp)
                               DO ay=0,la-ax
                                  az = la - ax - ay
                                  v(coset(ax,ay,az),coc,n) =&
                                       raw(1)*v(coset(ax-1,ay,az),coc,n+1) +&
                                       f6*(   v(coset(ax-2,ay,az),coc,n) +&
                                       f4*v(coset(ax-2,ay,az),coc,n+1)) +&
                                       fcx*v(coset(ax-1,ay,az),cocx,n+1)
                               END DO
                            END DO

                         END DO

                      END DO

                   END DO
                END DO

             END DO

          END IF

          DO j=ncoset(lc_min-1)+1,ncoset(lc_max)
             DO i=ncoset(la_min-1)+1,ncoset(la_max-maxder_local)
                vac(na+i,nc+j) = v(i,j,1)
             END DO
          END DO

          IF(PRESENT(maxder)) THEN
             DO j=ncoset(lc_min-1)+1,ncoset(lc_max)
                DO i=1,ncoset(la_max)
                   vac_plus(nap+i,nc+j) = v(i,j,1)
                END DO
             END DO
          END IF

          nc = nc + ncoset(lc_max)

       END DO

       na = na + ncoset(la_max - maxder_local)
       nap = nap + ncoset(la_max)

    END DO

  END SUBROUTINE coulomb2
! *****************************************************************************
!> \brief   Calculation of the primitive two-center Coulomb integrals over
!>          Cartesian Gaussian-type functions. Same as coulomb2 treats derivative
!>          different (now vac_plus is symmetric)
!> \author  CJM
!> \date    6.02.2008
!> \version 1.0
! *****************************************************************************

  SUBROUTINE coulomb2_new(la_max,npgfa,zeta,la_min,lc_max,npgfc,zetc,lc_min,&
       rac,rac2,vac,v,f,maxder,vac_plus)
    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta
    INTEGER, INTENT(IN)                      :: la_min, lc_max, npgfc
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetc
    INTEGER, INTENT(IN)                      :: lc_min
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac
    REAL(KIND=dp), INTENT(IN)                :: rac2
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: vac
    REAL(KIND=dp), DIMENSION(:, :, :)        :: v
    REAL(KIND=dp), DIMENSION(0:)             :: f
    INTEGER, INTENT(IN), OPTIONAL            :: maxder
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL :: vac_plus

    CHARACTER(len=*), PARAMETER :: routineN = 'coulomb2_new', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ax, ay, az, coc, cocx, cocy, &
                                                cocz, cx, cy, cz, i, ipgf, j, &
                                                jpgf, la, lc, maxder_local, &
                                                n, na, nap, nc, ncp, nmax
    REAL(KIND=dp)                            :: dac, f0, f1, f2, f3, f4, f5, &
                                                f6, fcx, fcy, fcz, rho, t, &
                                                zetp, zetq, zetw
    REAL(KIND=dp), DIMENSION(3)              :: raw, rcw

    v = 0.0_dp

    maxder_local = 0
    IF(PRESENT(maxder)) THEN
       maxder_local = maxder
       vac_plus = 0.0_dp
    END IF

    nmax = la_max + lc_max + 1

    ! *** Calculate the distance of the centers a and c ***

    dac = SQRT(rac2)

    ! *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0
    nap = 0

    DO ipgf=1,npgfa

       nc = 0
       ncp = 0

       DO jpgf=1,npgfc

          ! *** Calculate some prefactors ***

          zetp = 1.0_dp/zeta(ipgf)
          zetq = 1.0_dp/zetc(jpgf)
          zetw = 1.0_dp/(zeta(ipgf) + zetc(jpgf))

          rho = zeta(ipgf)*zetc(jpgf)*zetw

          f0 = 2.0_dp*SQRT(pi**5*zetw)*zetp*zetq

          ! *** Calculate the incomplete Gamma function ***

          t = rho*rac2

          CALL fgamma(nmax-1,t,f)

          ! *** Calculate the basic two-center Coulomb integrals [s||s]{n} ***

          DO n=1,nmax
             v(1,1,n) = f0*f(n-1)
          END DO

          ! *** Vertical recurrence steps: [s||s] -> [s||c] ***

          IF (lc_max > 0) THEN

             f1 = 0.5_dp*zetq
             f2 = -rho*zetq

             rcw(:) = -zeta(ipgf)*zetw*rac(:)

             ! *** [s||p]{n} = (Wi - Ci)*[s||s]{n+1}  (i = x,y,z) ***

             DO n=1,nmax-1
                v(1,2,n) = rcw(1)*v(1,1,n+1)
                v(1,3,n) = rcw(2)*v(1,1,n+1)
                v(1,4,n) = rcw(3)*v(1,1,n+1)
             END DO

             ! ** [s||c]{n} = (Wi - Ci)*[s||c-1i]{n+1} +     ***
             ! **             f1*Ni(c-1i)*(   [s||c-2i]{n} + ***
             ! **                          f2*[s||c-2i]{n+1} ***

             DO lc=2,lc_max

                DO n=1,nmax-lc

                   ! **** Increase the angular momentum component z of c ***

                   v(1,coset(0,0,lc),n) =&
                        rcw(3)*v(1,coset(0,0,lc-1),n+1) +&
                        f1*REAL(lc-1,dp)*(   v(1,coset(0,0,lc-2),n) +&
                        f2*v(1,coset(0,0,lc-2),n+1))

                   ! *** Increase the angular momentum component y of c ***

                   cz = lc - 1
                   v(1,coset(0,1,cz),n) = rcw(2)*v(1,coset(0,0,cz),n+1)

                   DO cy=2,lc
                      cz = lc - cy
                      v(1,coset(0,cy,cz),n) =&
                           rcw(2)*v(1,coset(0,cy-1,cz),n+1) +&
                           f1*REAL(cy-1,dp)*(   v(1,coset(0,cy-2,cz),n) +&
                           f2*v(1,coset(0,cy-2,cz),n+1))
                   END DO

                   ! *** Increase the angular momentum component x of c ***

                   DO cy=0,lc-1
                      cz = lc - 1 - cy
                      v(1,coset(1,cy,cz),n) = rcw(1)*v(1,coset(0,cy,cz),n+1)
                   END DO

                   DO cx=2,lc
                      f6 = f1*REAL(cx-1,dp)
                      DO cy=0,lc-cx
                         cz = lc - cx - cy
                         v(1,coset(cx,cy,cz),n) =&
                              rcw(1)*v(1,coset(cx-1,cy,cz),n+1) +&
                              f6*(   v(1,coset(cx-2,cy,cz),n) +&
                              f2*v(1,coset(cx-2,cy,cz),n+1))
                      END DO
                   END DO

                END DO

             END DO

          END IF

          ! *** Vertical recurrence steps: [s||c] -> [a||c] ***

          IF (la_max > 0) THEN

             f3 = 0.5_dp*zetp
             f4 = -rho*zetp
             f5 = 0.5_dp*zetw

             raw(:) = zetc(jpgf)*zetw*rac(:)

             ! *** [p||s]{n} = (Wi - Ai)*[s||s]{n+1}  (i = x,y,z) ***

             DO n=1,nmax-1
                v(2,1,n) = raw(1)*v(1,1,n+1)
                v(3,1,n) = raw(2)*v(1,1,n+1)
                v(4,1,n) = raw(3)*v(1,1,n+1)
             END DO

             ! *** [a||s]{n} = (Wi - Ai)*[a-1i||s]{n+1} +      ***
             ! ***             f3*Ni(a-1i)*(   [a-2i||s]{n} +  ***
             ! ***                          f4*[a-2i||s]{n+1}) ***

             DO la=2,la_max

                DO n=1,nmax-la

                   ! *** Increase the angular momentum component z of a ***

                   v(coset(0,0,la),1,n) =&
                        raw(3)*v(coset(0,0,la-1),1,n+1) +&
                        f3*REAL(la-1,dp)*(   v(coset(0,0,la-2),1,n) +&
                        f4*v(coset(0,0,la-2),1,n+1))

                   ! *** Increase the angular momentum component y of a ***

                   az = la - 1
                   v(coset(0,1,az),1,n) = raw(2)*v(coset(0,0,az),1,n+1)

                   DO ay=2,la
                      az = la - ay
                      v(coset(0,ay,az),1,n) =&
                           raw(2)*v(coset(0,ay-1,az),1,n+1) +&
                           f3*REAL(ay-1,dp)*(   v(coset(0,ay-2,az),1,n) +&
                           f4*v(coset(0,ay-2,az),1,n+1))
                   END DO

                   ! *** Increase the angular momentum component x of a ***

                   DO ay=0,la-1
                      az = la - 1 - ay
                      v(coset(1,ay,az),1,n) = raw(1)*v(coset(0,ay,az),1,n+1)
                   END DO

                   DO ax=2,la
                      f6 = f3*REAL(ax-1,dp)
                      DO ay=0,la-ax
                         az = la - ax - ay
                         v(coset(ax,ay,az),1,n) =&
                              raw(1)*v(coset(ax-1,ay,az),1,n+1) +&
                              f6*(   v(coset(ax-2,ay,az),1,n) +&
                              f4*v(coset(ax-2,ay,az),1,n+1))
                      END DO
                   END DO

                END DO

             END DO

             DO lc=1,lc_max

                DO cx=0,lc
                   DO cy=0,lc-cx
                      cz = lc - cx - cy

                      coc = coset(cx,cy,cz)
                      cocx = coset(MAX(0,cx-1),cy,cz)
                      cocy = coset(cx,MAX(0,cy-1),cz)
                      cocz = coset(cx,cy,MAX(0,cz-1))

                      fcx = f5*REAL(cx,dp)
                      fcy = f5*REAL(cy,dp)
                      fcz = f5*REAL(cz,dp)

                      ! *** [p||c]{n} = (Wi - Ai)*[s||c]{n+1} + ***
                      ! ***             f5*Ni(c)*[s||c-1i]{n+1} ***

                      DO n=1,nmax-1-lc
                         v(2,coc,n) = raw(1)*v(1,coc,n+1) + fcx*v(1,cocx,n+1)
                         v(3,coc,n) = raw(2)*v(1,coc,n+1) + fcy*v(1,cocy,n+1)
                         v(4,coc,n) = raw(3)*v(1,coc,n+1) + fcz*v(1,cocz,n+1)
                      END DO

                      ! *** [a||c]{n} = (Wi - Ai)*[a-1i||c]{n+1} +        ***
                      ! ***             f3*Ni(a-1i)*(   [a-2i||c]{n} +    ***
                      ! ***                          f4*[a-2i||c]{n+1}) + ***
                      ! ***             f5*Ni(c)*[a-1i||c-1i]{n+1}        ***

                      DO la=2,la_max

                         DO n=1,nmax-la-lc

                            ! *** Increase the angular momentum component z of a ***

                            v(coset(0,0,la),coc,n) =&
                                 raw(3)*v(coset(0,0,la-1),coc,n+1) +&
                                 f3*REAL(la-1,dp)*(   v(coset(0,0,la-2),coc,n) +&
                                 f4*v(coset(0,0,la-2),coc,n+1)) +&
                                 fcz*v(coset(0,0,la-1),cocz,n+1)

                            ! *** Increase the angular momentum component y of a ***

                            az = la - 1
                            v(coset(0,1,az),coc,n) =&
                                 raw(2)*v(coset(0,0,az),coc,n+1) +&
                                 fcy*v(coset(0,0,az),cocy,n+1)

                            DO ay=2,la
                               az = la - ay
                               v(coset(0,ay,az),coc,n) =&
                                    raw(2)*v(coset(0,ay-1,az),coc,n+1) +&
                                    f3*REAL(ay-1,dp)*(   v(coset(0,ay-2,az),coc,n) +&
                                    f4*v(coset(0,ay-2,az),coc,n+1)) +&
                                    fcy*v(coset(0,ay-1,az),cocy,n+1)
                            END DO

                            ! *** Increase the angular momentum component x of a ***

                            DO ay=0,la-1
                               az = la - 1 - ay
                               v(coset(1,ay,az),coc,n) =&
                                    raw(1)*v(coset(0,ay,az),coc,n+1) +&
                                    fcx*v(coset(0,ay,az),cocx,n+1)
                            END DO

                            DO ax=2,la
                               f6 = f3*REAL(ax-1,dp)
                               DO ay=0,la-ax
                                  az = la - ax - ay
                                  v(coset(ax,ay,az),coc,n) =&
                                       raw(1)*v(coset(ax-1,ay,az),coc,n+1) +&
                                       f6*(   v(coset(ax-2,ay,az),coc,n) +&
                                       f4*v(coset(ax-2,ay,az),coc,n+1)) +&
                                       fcx*v(coset(ax-1,ay,az),cocx,n+1)
                               END DO
                            END DO

                         END DO

                      END DO

                   END DO
                END DO

             END DO

          END IF

          DO j=ncoset(lc_min-1)+1,ncoset(lc_max-maxder_local)
             DO i=ncoset(la_min-1)+1,ncoset(la_max-maxder_local)
                vac(na+i,nc+j) = v(i,j,1)
             END DO
          END DO

          IF(PRESENT(maxder)) THEN
             DO j=1,ncoset(lc_max)
                DO i=1,ncoset(la_max)
                   vac_plus(nap+i,ncp+j) = v(i,j,1)
                END DO
             END DO
          END IF

          nc = nc + ncoset ( lc_max - maxder_local )
          ncp = ncp + ncoset ( lc_max )

       END DO

       na = na + ncoset(la_max - maxder_local)
       nap = nap + ncoset(la_max)

    END DO

  END SUBROUTINE coulomb2_new


! *****************************************************************************
!> \brief   Calculation of the primitive three-center Coulomb integrals over
!>          Cartesian Gaussian-type functions (electron repulsion integrals,
!>          ERIs).
!> \author  Matthias Krack
!> \date    06.11.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE coulomb3(la_max,npgfa,zeta,rpgfa,la_min,lb_max,npgfb,zetb,rpgfb,lb_min,&
       lc_max,zetc,rpgfc,lc_min,gccc,rab,rab2,rac,rac2,rbc2,vabc,int_abc,&
       v,f,maxder,vabc_plus)
    INTEGER, INTENT(IN)                      :: la_max, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: la_min, lb_max, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb, rpgfb
    INTEGER, INTENT(IN)                      :: lb_min, lc_max
    REAL(KIND=dp), INTENT(IN)                :: zetc, rpgfc
    INTEGER, INTENT(IN)                      :: lc_min
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: gccc
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: rab2
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac
    REAL(KIND=dp), INTENT(IN)                :: rac2, rbc2
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: vabc
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(OUT)                            :: int_abc
    REAL(KIND=dp), DIMENSION(:, :, :, :)     :: v
    REAL(KIND=dp), DIMENSION(0:)             :: f
    INTEGER, INTENT(IN), OPTIONAL            :: maxder
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL :: vabc_plus

    CHARACTER(len=*), PARAMETER :: routineN = 'coulomb3', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ay, az, bx, by, bz, coc, cocx, cocy, cocz, cx, cy, cz, i, &
      ipgf, j, jpgf, k, kk, la, la_start, lb, lc, maxder_local, n, na, nap, &
      nb, nmax
    REAL(KIND=dp)                            :: dab, dac, dbc, f0, f1, f2, &
                                                f3, f4, f5, f6, f7, fcx, fcy, &
                                                fcz, fx, fy, fz, t, zetp, &
                                                zetq, zetw
    REAL(KIND=dp), DIMENSION(3)              :: rap, rbp, rcp, rcw, rpw

    v = 0.0_dp

    maxder_local = 0
    IF(PRESENT(maxder)) THEN
       maxder_local = maxder
    END IF

    nmax = la_max + lb_max + lc_max + 1

    ! *** Calculate the distances of the centers a, b and c ***

    dab = SQRT(rab2)
    dac = SQRT(rac2)
    dbc = SQRT(rbc2)

    ! *** Initialize integrals array
    int_abc = 0.0_dp

    ! *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0
    nap = 0

    DO ipgf=1,npgfa

       ! *** Screening ***
       IF (rpgfa(ipgf) + rpgfc < dac) THEN
          na = na + ncoset(la_max - maxder_local)
          nap = nap + ncoset(la_max)
          CYCLE
       END IF

       nb = 0

       DO jpgf=1,npgfb

          ! *** Screening ***
          IF (&
               (rpgfb(jpgf) + rpgfc < dbc).OR.&
               (rpgfa(ipgf) + rpgfb(jpgf) < dab)) THEN
             nb = nb + ncoset(lb_max)
             CYCLE
          END IF

          ! *** Calculate some prefactors ***

          zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))
          zetq = 1.0_dp/zetc
          zetw = 1.0_dp/(zeta(ipgf) + zetb(jpgf) + zetc)

          f0 = 2.0_dp*SQRT(pi**5*zetw)*zetp*zetq
          f1 = zetb(jpgf)*zetp
          f2 = 0.5_dp*zetp
          f4 = -zetc*zetw

          f0 = f0*EXP(-zeta(ipgf)*f1*rab2)

          rap(:) = f1*rab(:)
          rcp(:) = rap(:) - rac(:)
          rpw(:) = f4*rcp(:)

          ! *** Calculate the incomplete Gamma function ***

          t = -f4*(rcp(1)*rcp(1) + rcp(2)*rcp(2) + rcp(3)*rcp(3))/zetp

          CALL fgamma(nmax-1,t,f)

          ! *** Calculate the basic three-center Coulomb integrals [ss||s]{n} ***

          DO n=1,nmax
             v(1,1,1,n) = f0*f(n-1)
          END DO

          ! *** Recurrence steps: [ss||s] -> [as||s] ***

          IF (la_max > 0) THEN

             ! *** Vertical recurrence steps: [ss||s] -> [as||s] ***

             ! *** [ps||s]{n} = (Pi - Ai)*[ss||s]{n} +              ***
             ! ***              (Wi - Pi)*[ss||s]{n+1}  (i = x,y,z) ***

             DO n=1,nmax-1
                v(2,1,1,n) = rap(1)*v(1,1,1,n) + rpw(1)*v(1,1,1,n+1)
                v(3,1,1,n) = rap(2)*v(1,1,1,n) + rpw(2)*v(1,1,1,n+1)
                v(4,1,1,n) = rap(3)*v(1,1,1,n) + rpw(3)*v(1,1,1,n+1)
             END DO

             ! *** [as||s]{n} = (Pi - Ai)*[(a-1i)s||s]{n} +        ***
             ! ***              (Wi - Pi)*[(a-1i)s||s]{n+1} +      ***
             ! ***              f2*Ni(a-1i)*(   [(a-2i)s||s]{n} +  ***
             ! ***                           f4*[(a-2i)s||s]{n+1}) ***

             DO la=2,la_max

                DO n=1,nmax-la

                   ! *** Increase the angular momentum component z of a ***

                   v(coset(0,0,la),1,1,n) =&
                        rap(3)*v(coset(0,0,la-1),1,1,n) +&
                        rpw(3)*v(coset(0,0,la-1),1,1,n+1) +&
                        f2*REAL(la-1,dp)*(   v(coset(0,0,la-2),1,1,n) +&
                        f4*v(coset(0,0,la-2),1,1,n+1))

                   ! *** Increase the angular momentum component y of a ***

                   az = la - 1
                   v(coset(0,1,az),1,1,n) =&
                        rap(2)*v(coset(0,0,az),1,1,n) +&
                        rpw(2)*v(coset(0,0,az),1,1,n+1)

                   DO ay=2,la
                      az = la - ay
                      v(coset(0,ay,az),1,1,n) =&
                           rap(2)*v(coset(0,ay-1,az),1,1,n) +&
                           rpw(2)*v(coset(0,ay-1,az),1,1,n+1) +&
                           f2*REAL(ay-1,dp)*(   v(coset(0,ay-2,az),1,1,n) +&
                           f4*v(coset(0,ay-2,az),1,1,n+1))
                   END DO

                   ! *** Increase the angular momentum component x of a ***

                   DO ay=0,la-1
                      az = la - 1 - ay
                      v(coset(1,ay,az),1,1,n) =&
                           rap(1)*v(coset(0,ay,az),1,1,n) +&
                           rpw(1)*v(coset(0,ay,az),1,1,n+1)
                   END DO

                   DO ax=2,la
                      f3 = f2*REAL(ax-1,dp)
                      DO ay=0,la-ax
                         az = la - ax - ay
                         v(coset(ax,ay,az),1,1,n) =&
                              rap(1)*v(coset(ax-1,ay,az),1,1,n) +&
                              rpw(1)*v(coset(ax-1,ay,az),1,1,n+1) +&
                              f3*(   v(coset(ax-2,ay,az),1,1,n) +&
                              f4*v(coset(ax-2,ay,az),1,1,n+1))
                      END DO
                   END DO

                END DO

             END DO

             ! *** Recurrence steps: [as||s] -> [ab||s] ***

             IF (lb_max > 0) THEN

                ! *** Horizontal recurrence steps ***

                rbp(:) = rap(:) - rab(:)

                ! *** [ap||s]{n} = [(a+1i)s||s]{n} - (Bi - Ai)*[as||s]{n} ***

                la_start = MAX(0,la_min-1)

                DO la=la_start,la_max-1
                   DO n=1,nmax-la-1
                      DO ax=0,la
                         DO ay=0,la-ax
                            az = la - ax - ay
                            v(coset(ax,ay,az),2,1,n) =&
                                 v(coset(ax+1,ay,az),1,1,n) -&
                                 rab(1)*v(coset(ax,ay,az),1,1,n)
                            v(coset(ax,ay,az),3,1,n) =&
                                 v(coset(ax,ay+1,az),1,1,n) -&
                                 rab(2)*v(coset(ax,ay,az),1,1,n)
                            v(coset(ax,ay,az),4,1,n) =&
                                 v(coset(ax,ay,az+1),1,1,n) -&
                                 rab(3)*v(coset(ax,ay,az),1,1,n)
                         END DO
                      END DO
                   END DO
                END DO

                ! *** Vertical recurrence step ***

                ! *** [ap||s]{n} = (Pi - Bi)*[as||s]{n} +          ***
                ! ***              (Wi - Pi)*[as||s]{n+1} +        ***
                ! ***              f2*Ni(a)*(   [(a-1i)s||s]{n} +  ***
                ! ***                        f4*[(a-1i)s||s]{n+1}) ***

                DO n=1,nmax-la_max-1
                   DO ax=0,la_max
                      fx = f2*REAL(ax,dp)
                      DO ay=0,la_max-ax
                         fy = f2*REAL(ay,dp)
                         az = la_max - ax - ay
                         fz = f2*REAL(az,dp)

                         IF (ax == 0) THEN
                            v(coset(ax,ay,az),2,1,n) =&
                                 rbp(1)*v(coset(ax,ay,az),1,1,n) +&
                                 rpw(1)*v(coset(ax,ay,az),1,1,n+1)
                         ELSE
                            v(coset(ax,ay,az),2,1,n) =&
                                 rbp(1)*v(coset(ax,ay,az),1,1,n) +&
                                 rpw(1)*v(coset(ax,ay,az),1,1,n+1) +&
                                 fx*(   v(coset(ax-1,ay,az),1,1,n) +&
                                 f4*v(coset(ax-1,ay,az),1,1,n+1))
                         END IF

                         IF (ay == 0) THEN
                            v(coset(ax,ay,az),3,1,n) =&
                                 rbp(2)*v(coset(ax,ay,az),1,1,n) +&
                                 rpw(2)*v(coset(ax,ay,az),1,1,n+1)
                         ELSE
                            v(coset(ax,ay,az),3,1,n) =&
                                 rbp(2)*v(coset(ax,ay,az),1,1,n) +&
                                 rpw(2)*v(coset(ax,ay,az),1,1,n+1) +&
                                 fy*(   v(coset(ax,ay-1,az),1,1,n) +&
                                 f4*v(coset(ax,ay-1,az),1,1,n+1))
                         END IF

                         IF (az == 0) THEN
                            v(coset(ax,ay,az),4,1,n) =&
                                 rbp(3)*v(coset(ax,ay,az),1,1,n) +&
                                 rpw(3)*v(coset(ax,ay,az),1,1,n+1)
                         ELSE
                            v(coset(ax,ay,az),4,1,n) =&
                                 rbp(3)*v(coset(ax,ay,az),1,1,n) +&
                                 rpw(3)*v(coset(ax,ay,az),1,1,n+1) +&
                                 fz*(   v(coset(ax,ay,az-1),1,1,n) +&
                                 f4*v(coset(ax,ay,az-1),1,1,n+1))
                         END IF

                      END DO
                   END DO
                END DO

                ! *** Recurrence steps: [ap||s] -> [ab||s] ***

                DO lb=2,lb_max

                   ! *** Horizontal recurrence steps ***

                   ! *** [ab||s]{n} = [(a+1i)(b-1i)||s]{n} -    ***
                   ! ***              (Bi - Ai)*[a(b-1i)||s]{n} ***

                   la_start = MAX(0,la_min-1)

                   DO la=la_start,la_max-1
                      DO n=1,nmax-la-lb
                         DO ax=0,la
                            DO ay=0,la-ax
                               az = la - ax - ay

                               ! *** Shift of angular momentum component z from a to b ***

                               v(coset(ax,ay,az),coset(0,0,lb),1,n) =&
                                    v(coset(ax,ay,az+1),coset(0,0,lb-1),1,n) -&
                                    rab(3)*v(coset(ax,ay,az),coset(0,0,lb-1),1,n)

                               ! *** Shift of angular momentum component y from a to b ***

                               DO by=1,lb
                                  bz = lb - by
                                  v(coset(ax,ay,az),coset(0,by,bz),1,n) =&
                                       v(coset(ax,ay+1,az),coset(0,by-1,bz),1,n) -&
                                       rab(2)*v(coset(ax,ay,az),coset(0,by-1,bz),1,n)
                               END DO

                               ! *** Shift of angular momentum component x from a to b ***

                               DO bx=1,lb
                                  DO by=0,lb-bx
                                     bz = lb - bx - by
                                     v(coset(ax,ay,az),coset(bx,by,bz),1,n) =&
                                          v(coset(ax+1,ay,az),coset(bx-1,by,bz),1,n) -&
                                          rab(1)*v(coset(ax,ay,az),coset(bx-1,by,bz),1,n)
                                  END DO
                               END DO

                            END DO
                         END DO
                      END DO
                   END DO

                   ! *** Vertical recurrence step ***

                   ! *** [ab||s]{n} = (Pi - Bi)*[a(b-1i)||s]{n} +          ***
                   ! ***              (Wi - Pi)*[a(b-1i)||s]{n+1} +        ***
                   ! ***              f2*Ni(a)*(   [(a-1i)(b-1i)||s]{n} +  ***
                   ! ***                        f4*[(a-1i)(b-1i)||s]{n+1}) ***
                   ! ***              f2*Ni(b-1i)*(   [a(b-2i)||s]{n} +    ***
                   ! ***                           f4*[a(b-2i)||s]{n+1})   ***

                   DO n=1,nmax-la_max-lb
                      DO ax=0,la_max
                         fx = f2*REAL(ax,dp)
                         DO ay=0,la_max-ax
                            fy = f2*REAL(ay,dp)
                            az = la_max - ax - ay
                            fz = f2*REAL(az,dp)

                            ! *** Shift of angular momentum component z from a to b ***

                            f3 = f2*REAL(lb-1,dp)

                            IF (az == 0) THEN
                               v(coset(ax,ay,az),coset(0,0,lb),1,n) =&
                                    rbp(3)*v(coset(ax,ay,az),coset(0,0,lb-1),1,n) +&
                                    rpw(3)*v(coset(ax,ay,az),coset(0,0,lb-1),1,n+1) +&
                                    f3*(   v(coset(ax,ay,az),coset(0,0,lb-2),1,n) +&
                                    f4*v(coset(ax,ay,az),coset(0,0,lb-2),1,n+1))
                            ELSE
                               v(coset(ax,ay,az),coset(0,0,lb),1,n) =&
                                    rbp(3)*v(coset(ax,ay,az),coset(0,0,lb-1),1,n) +&
                                    rpw(3)*v(coset(ax,ay,az),coset(0,0,lb-1),1,n+1) +&
                                    fz*(   v(coset(ax,ay,az-1),coset(0,0,lb-1),1,n) +&
                                    f4*v(coset(ax,ay,az-1),coset(0,0,lb-1),1,n+1)) +&
                                    f3*(   v(coset(ax,ay,az),coset(0,0,lb-2),1,n) +&
                                    f4*v(coset(ax,ay,az),coset(0,0,lb-2),1,n+1))
                            END IF

                            ! *** Shift of angular momentum component y from a to b ***

                            IF (ay == 0) THEN
                               bz = lb - 1
                               v(coset(ax,ay,az),coset(0,1,bz),1,n) =&
                                    rbp(2)*v(coset(ax,ay,az),coset(0,0,bz),1,n) +&
                                    rpw(2)*v(coset(ax,ay,az),coset(0,0,bz),1,n+1)
                               DO by=2,lb
                                  bz = lb - by
                                  f3 = f2*REAL(by-1,dp)
                                  v(coset(ax,ay,az),coset(0,by,bz),1,n) =&
                                       rbp(2)*v(coset(ax,ay,az),coset(0,by-1,bz),1,n) +&
                                       rpw(2)*v(coset(ax,ay,az),coset(0,by-1,bz),1,n+1) +&
                                       f3*(   v(coset(ax,ay,az),coset(0,by-2,bz),1,n) +&
                                       f4*v(coset(ax,ay,az),coset(0,by-2,bz),1,n+1))
                               END DO
                            ELSE
                               bz = lb - 1
                               v(coset(ax,ay,az),coset(0,1,bz),1,n) =&
                                    rbp(2)*v(coset(ax,ay,az),coset(0,0,bz),1,n) +&
                                    rpw(2)*v(coset(ax,ay,az),coset(0,0,bz),1,n+1) +&
                                    fy*(   v(coset(ax,ay-1,az),coset(0,0,bz),1,n) +&
                                    f4*v(coset(ax,ay-1,az),coset(0,0,bz),1,n+1))
                               DO by=2,lb
                                  bz = lb - by
                                  f3 = f2*REAL(by-1,dp)
                                  v(coset(ax,ay,az),coset(0,by,bz),1,n) =&
                                       rbp(2)*v(coset(ax,ay,az),coset(0,by-1,bz),1,n) +&
                                       rpw(2)*v(coset(ax,ay,az),coset(0,by-1,bz),1,n+1) +&
                                       fy*(   v(coset(ax,ay-1,az),coset(0,by-1,bz),1,n) +&
                                       f4*v(coset(ax,ay-1,az),&
                                       coset(0,by-1,bz),1,n+1)) +&
                                       f3*(   v(coset(ax,ay,az),coset(0,by-2,bz),1,n) +&
                                       f4*v(coset(ax,ay,az),coset(0,by-2,bz),1,n+1))
                               END DO
                            END IF

                            ! *** Shift of angular momentum component x from a to b ***

                            IF (ax == 0) THEN
                               DO by=0,lb-1
                                  bz = lb - 1 - by
                                  v(coset(ax,ay,az),coset(1,by,bz),1,n) =&
                                       rbp(1)*v(coset(ax,ay,az),coset(0,by,bz),1,n) +&
                                       rpw(1)*v(coset(ax,ay,az),coset(0,by,bz),1,n+1)
                               END DO
                               DO bx=2,lb
                                  f3 = f2*REAL(bx-1,dp)
                                  DO by=0,lb-bx
                                     bz = lb - bx - by
                                     v(coset(ax,ay,az),coset(bx,by,bz),1,n) =&
                                          rbp(1)*v(coset(ax,ay,az),coset(bx-1,by,bz),1,n) +&
                                          rpw(1)*v(coset(ax,ay,az),&
                                          coset(bx-1,by,bz),1,n+1) +&
                                          f3*(   v(coset(ax,ay,az),coset(bx-2,by,bz),1,n) +&
                                          f4*v(coset(ax,ay,az),coset(bx-2,by,bz),1,n+1))
                                  END DO
                               END DO
                            ELSE
                               DO by=0,lb-1
                                  bz = lb - 1 - by
                                  v(coset(ax,ay,az),coset(1,by,bz),1,n) =&
                                       rbp(1)*v(coset(ax,ay,az),coset(0,by,bz),1,n) +&
                                       rpw(1)*v(coset(ax,ay,az),coset(0,by,bz),1,n+1) +&
                                       fx*(   v(coset(ax-1,ay,az),coset(0,by,bz),1,n) +&
                                       f4*v(coset(ax-1,ay,az),coset(0,by,bz),1,n+1))
                               END DO
                               DO bx=2,lb
                                  f3 = f2*REAL(bx-1,dp)
                                  DO by=0,lb-bx
                                     bz = lb - bx - by
                                     v(coset(ax,ay,az),coset(bx,by,bz),1,n) =&
                                          rbp(1)*v(coset(ax,ay,az),coset(bx-1,by,bz),1,n) +&
                                          rpw(1)*v(coset(ax,ay,az),&
                                          coset(bx-1,by,bz),1,n+1) +&
                                          fx*(   v(coset(ax-1,ay,az),&
                                          coset(bx-1,by,bz),1,n) +&
                                          f4*v(coset(ax-1,ay,az),&
                                          coset(bx-1,by,bz),1,n+1)) +&
                                          f3*(   v(coset(ax,ay,az),coset(bx-2,by,bz),1,n) +&
                                          f4*v(coset(ax,ay,az),coset(bx-2,by,bz),1,n+1))
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

                ! *** Vertical recurrence steps: [ss||s] -> [sb||s] ***

                rbp(:) = rap(:) - rab(:)

                ! *** [sp||s]{n} = (Pi - Bi)*[ss||s]{n} + ***
                ! ***              (Wi - Pi)*[ss||s]{n+1} ***

                DO n=1,nmax-1
                   v(1,2,1,n) = rbp(1)*v(1,1,1,n) + rpw(1)*v(1,1,1,n+1)
                   v(1,3,1,n) = rbp(2)*v(1,1,1,n) + rpw(2)*v(1,1,1,n+1)
                   v(1,4,1,n) = rbp(3)*v(1,1,1,n) + rpw(3)*v(1,1,1,n+1)
                END DO

                ! *** [sb||s]{n} = (Pi - Bi)*[s(b-1i)||s]{n} +        ***
                ! ***              (Wi - Pi)*[s(b-1i)||s]{n+1} +      ***
                ! ***              f2*Ni(b-1i)*(   [s(b-2i)||s]{n} +  ***
                ! ***                           f4*[s(b-2i)||s]{n+1}) ***

                DO lb=2,lb_max

                   DO n=1,nmax-lb

                      ! *** Increase the angular momentum component z of b ***

                      v(1,coset(0,0,lb),1,n) =&
                           rbp(3)*v(1,coset(0,0,lb-1),1,n) +&
                           rpw(3)*v(1,coset(0,0,lb-1),1,n+1) +&
                           f2*REAL(lb-1,dp)*(   v(1,coset(0,0,lb-2),1,n) +&
                           f4*v(1,coset(0,0,lb-2),1,n+1))

                      ! *** Increase the angular momentum component y of b ***

                      bz = lb - 1
                      v(1,coset(0,1,bz),1,n) =&
                           rbp(2)*v(1,coset(0,0,bz),1,n) +&
                           rpw(2)*v(1,coset(0,0,bz),1,n+1)

                      DO by=2,lb
                         bz = lb - by
                         v(1,coset(0,by,bz),1,n) =&
                              rbp(2)*v(1,coset(0,by-1,bz),1,n) +&
                              rpw(2)*v(1,coset(0,by-1,bz),1,n+1) +&
                              f2*REAL(by-1,dp)*(   v(1,coset(0,by-2,bz),1,n) +&
                              f4*v(1,coset(0,by-2,bz),1,n+1))
                      END DO

                      ! *** Increase the angular momentum component x of b ***

                      DO by=0,lb-1
                         bz = lb - 1 - by
                         v(1,coset(1,by,bz),1,n) =&
                              rbp(1)*v(1,coset(0,by,bz),1,n) +&
                              rpw(1)*v(1,coset(0,by,bz),1,n+1)
                      END DO

                      DO bx=2,lb
                         f3 = f2*REAL(bx-1,dp)
                         DO by=0,lb-bx
                            bz = lb - bx - by
                            v(1,coset(bx,by,bz),1,n) =&
                                 rbp(1)*v(1,coset(bx-1,by,bz),1,n) +&
                                 rpw(1)*v(1,coset(bx-1,by,bz),1,n+1) +&
                                 f3*(   v(1,coset(bx-2,by,bz),1,n) +&
                                 f4*v(1,coset(bx-2,by,bz),1,n+1))
                         END DO
                      END DO

                   END DO

                END DO

             END IF

          END IF

          ! *** Recurrence steps: [ab||s] -> [ab||c] ***

          IF (lc_max > 0) THEN

             ! *** Vertical recurrence steps: [ss||s] -> [ss||c] ***

             f5 = -zetw/zetp
             f6 = 0.5_dp*zetw
             f7 = 0.5_dp*zetq

             rcw(:) = rcp(:) + rpw(:)

             ! *** [ss||p]{n} = (Wi - Ci)*[ss||s]{n+1}  (i = x,y,z) ***

             DO n=1,nmax-1
                v(1,1,2,n) = rcw(1)*v(1,1,1,n+1)
                v(1,1,3,n) = rcw(2)*v(1,1,1,n+1)
                v(1,1,4,n) = rcw(3)*v(1,1,1,n+1)
             END DO

             ! *** [ss||c]{n} = (Wi - Ci)*[ss||c-1i]{n+1} + ***
             ! ***              f7*Ni(c-1i)*[ss||c-2i]{n} + ***
             ! ***              f5*[ss||c-2i]{n+1}          ***

             DO lc=2,lc_max

                DO n=1,nmax-lc

                   ! *** Increase the angular momentum component z of c ***

                   v(1,1,coset(0,0,lc),n) =&
                        rcw(3)*v(1,1,coset(0,0,lc-1),n+1) +&
                        f7*REAL(lc-1,dp)*(   v(1,1,coset(0,0,lc-2),n) +&
                        f5*v(1,1,coset(0,0,lc-2),n+1))

                   ! *** Increase the angular momentum component y of c ***

                   cz = lc - 1
                   v(1,1,coset(0,1,cz),n) = rcw(2)*v(1,1,coset(0,0,cz),n+1)

                   DO cy=2,lc
                      cz = lc - cy
                      v(1,1,coset(0,cy,cz),n) =&
                           rcw(2)*v(1,1,coset(0,cy-1,cz),n+1) +&
                           f7*REAL(cy-1,dp)*(   v(1,1,coset(0,cy-2,cz),n) +&
                           f5*v(1,1,coset(0,cy-2,cz),n+1))
                   END DO

                   ! *** Increase the angular momentum component x of c ***

                   DO cy=0,lc-1
                      cz = lc - 1 - cy
                      v(1,1,coset(1,cy,cz),n) = rcw(1)*v(1,1,coset(0,cy,cz),n+1)
                   END DO

                   DO cx=2,lc
                      DO cy=0,lc-cx
                         cz = lc - cx - cy
                         v(1,1,coset(cx,cy,cz),n) =&
                              rcw(1)*v(1,1,coset(cx-1,cy,cz),n+1) +&
                              f7*REAL(cx-1,dp)*(   v(1,1,coset(cx-2,cy,cz),n) +&
                              f5*v(1,1,coset(cx-2,cy,cz),n+1))
                      END DO
                   END DO

                END DO

             END DO

             ! *** Recurrence steps: [ss||c] -> [ab||c] ***

             DO lc=1,lc_max

                DO cx=0,lc
                   DO cy=0,lc-cx
                      cz = lc - cx - cy

                      coc = coset(cx,cy,cz)
                      cocx = coset(MAX(0,cx-1),cy,cz)
                      cocy = coset(cx,MAX(0,cy-1),cz)
                      cocz = coset(cx,cy,MAX(0,cz-1))

                      fcx = f6*REAL(cx,dp)
                      fcy = f6*REAL(cy,dp)
                      fcz = f6*REAL(cz,dp)

                      ! *** Recurrence steps: [ss||c] -> [as||c] ***

                      IF (la_max > 0) THEN

                         ! *** Vertical recurrence steps: [ss||c] -> [as||c] ***

                         ! *** [ps||c]{n} = (Pi - Ai)*[ss||c]{n} +                ***
                         ! ***              (Wi - Pi)*[ss||c]{n+1} +              ***
                         ! ***              f6*Ni(c)*[ss||c-1i]{n+1}  (i = x,y,z) ***

                         DO n=1,nmax-1-lc
                            v(2,1,coc,n) = rap(1)*v(1,1,coc,n) +&
                                 rpw(1)*v(1,1,coc,n+1) +&
                                 fcx*v(1,1,cocx,n+1)
                            v(3,1,coc,n) = rap(2)*v(1,1,coc,n) +&
                                 rpw(2)*v(1,1,coc,n+1) +&
                                 fcy*v(1,1,cocy,n+1)
                            v(4,1,coc,n) = rap(3)*v(1,1,coc,n) +&
                                 rpw(3)*v(1,1,coc,n+1) +&
                                 fcz*v(1,1,cocz,n+1)
                         END DO

                         ! *** [as||c]{n} = (Pi - Ai)*[(a-1i)s||c]{n} +          ***
                         ! ***              (Wi - Pi)*[(a-1i)s||c]{n+1} +        ***
                         ! ***              f2*Ni(a-1i)*(   [(a-2i)s||c]{n} +    ***
                         ! ***                           f4*[(a-2i)s||c]{n+1}) + ***
                         ! ***              f6*Ni(c)*[(a-1i)s||c-1i]{n+1}        ***

                         DO la=2,la_max

                            DO n=1,nmax-la-lc

                               ! *** Increase the angular momentum component z of a ***

                               v(coset(0,0,la),1,coc,n) =&
                                    rap(3)*v(coset(0,0,la-1),1,coc,n) +&
                                    rpw(3)*v(coset(0,0,la-1),1,coc,n+1) +&
                                    f2*REAL(la-1,dp)*(   v(coset(0,0,la-2),1,coc,n) +&
                                    f4*v(coset(0,0,la-2),1,coc,n+1)) +&
                                    fcz*v(coset(0,0,la-1),1,cocz,n+1)

                               ! *** Increase the angular momentum component y of a ***

                               az = la - 1
                               v(coset(0,1,az),1,coc,n) =&
                                    rap(2)*v(coset(0,0,az),1,coc,n) +&
                                    rpw(2)*v(coset(0,0,az),1,coc,n+1) +&
                                    fcy*v(coset(0,0,az),1,cocy,n+1)

                               DO ay=2,la
                                  f3 = f2*REAL(ay-1,dp)
                                  az = la - ay
                                  v(coset(0,ay,az),1,coc,n) =&
                                       rap(2)*v(coset(0,ay-1,az),1,coc,n) +&
                                       rpw(2)*v(coset(0,ay-1,az),1,coc,n+1) +&
                                       f3*(   v(coset(0,ay-2,az),1,coc,n) +&
                                       f4*v(coset(0,ay-2,az),1,coc,n+1)) +&
                                       fcy*v(coset(0,ay-1,az),1,cocy,n+1)
                               END DO

                               ! *** Increase the angular momentum component x of a ***

                               DO ay=0,la-1
                                  az = la - 1 - ay
                                  v(coset(1,ay,az),1,coc,n) =&
                                       rap(1)*v(coset(0,ay,az),1,coc,n) +&
                                       rpw(1)*v(coset(0,ay,az),1,coc,n+1) +&
                                       fcx*v(coset(0,ay,az),1,cocx,n+1)
                               END DO

                               DO ax=2,la
                                  f3 = f2*REAL(ax-1,dp)
                                  DO ay=0,la-ax
                                     az = la - ax - ay
                                     v(coset(ax,ay,az),1,coc,n) =&
                                          rap(1)*v(coset(ax-1,ay,az),1,coc,n) +&
                                          rpw(1)*v(coset(ax-1,ay,az),1,coc,n+1) +&
                                          f3*(   v(coset(ax-2,ay,az),1,coc,n) +&
                                          f4*v(coset(ax-2,ay,az),1,coc,n+1)) +&
                                          fcx*v(coset(ax-1,ay,az),1,cocx,n+1)
                                  END DO
                               END DO

                            END DO

                         END DO

                         ! *** Recurrence steps: [as||c] -> [ab||c] ***

                         IF (lb_max > 0) THEN

                            ! *** Horizontal recurrence steps ***

                            ! *** [ap||c]{n} = [(a+1i)s||c]{n} - (Bi - Ai)*[as||c]{n} ***

                            la_start = MAX(0,la_min-1)

                            DO la=la_start,la_max-1
                               DO n=1,nmax-la-1-lc
                                  DO ax=0,la
                                     DO ay=0,la-ax
                                        az = la - ax - ay
                                        v(coset(ax,ay,az),2,coc,n) =&
                                             v(coset(ax+1,ay,az),1,coc,n) -&
                                             rab(1)*v(coset(ax,ay,az),1,coc,n)
                                        v(coset(ax,ay,az),3,coc,n) =&
                                             v(coset(ax,ay+1,az),1,coc,n) -&
                                             rab(2)*v(coset(ax,ay,az),1,coc,n)
                                        v(coset(ax,ay,az),4,coc,n) =&
                                             v(coset(ax,ay,az+1),1,coc,n) -&
                                             rab(3)*v(coset(ax,ay,az),1,coc,n)
                                     END DO
                                  END DO
                               END DO
                            END DO

                            ! *** Vertical recurrence step ***

                            ! *** [ap||c]{n} = (Pi - Bi)*[as||c]{n} +            ***
                            ! ***              (Wi - Pi)*[as||c]{n+1} +          ***
                            ! ***              f2*Ni(a)*(   [(a-1i)s||c]{n} +    ***
                            ! ***                        f4*[(a-1i)s||c]{n+1}) + ***
                            ! ***              f6*Ni(c)*[(as||c-1i]{n+1})        ***

                            DO n=1,nmax-la_max-1-lc
                               DO ax=0,la_max
                                  fx = f2*REAL(ax,dp)
                                  DO ay=0,la_max-ax
                                     fy = f2*REAL(ay,dp)
                                     az = la_max - ax - ay
                                     fz = f2*REAL(az,dp)

                                     IF (ax == 0) THEN
                                        v(coset(ax,ay,az),2,coc,n) =&
                                             rbp(1)*v(coset(ax,ay,az),1,coc,n) +&
                                             rpw(1)*v(coset(ax,ay,az),1,coc,n+1) +&
                                             fcx*v(coset(ax,ay,az),1,cocx,n+1)
                                     ELSE
                                        v(coset(ax,ay,az),2,coc,n) =&
                                             rbp(1)*v(coset(ax,ay,az),1,coc,n) +&
                                             rpw(1)*v(coset(ax,ay,az),1,coc,n+1) +&
                                             fx*(   v(coset(ax-1,ay,az),1,coc,n) +&
                                             f4*v(coset(ax-1,ay,az),1,coc,n+1)) +&
                                             fcx*v(coset(ax,ay,az),1,cocx,n+1)
                                     END IF

                                     IF (ay == 0) THEN
                                        v(coset(ax,ay,az),3,coc,n) =&
                                             rbp(2)*v(coset(ax,ay,az),1,coc,n) +&
                                             rpw(2)*v(coset(ax,ay,az),1,coc,n+1) +&
                                             fcy*v(coset(ax,ay,az),1,cocy,n+1)
                                     ELSE
                                        v(coset(ax,ay,az),3,coc,n) =&
                                             rbp(2)*v(coset(ax,ay,az),1,coc,n) +&
                                             rpw(2)*v(coset(ax,ay,az),1,coc,n+1) +&
                                             fy*(   v(coset(ax,ay-1,az),1,coc,n) +&
                                             f4*v(coset(ax,ay-1,az),1,coc,n+1)) +&
                                             fcy*v(coset(ax,ay,az),1,cocy,n+1)
                                     END IF

                                     IF (az == 0) THEN
                                        v(coset(ax,ay,az),4,coc,n) =&
                                             rbp(3)*v(coset(ax,ay,az),1,coc,n) +&
                                             rpw(3)*v(coset(ax,ay,az),1,coc,n+1) +&
                                             fcz*v(coset(ax,ay,az),1,cocz,n+1)
                                     ELSE
                                        v(coset(ax,ay,az),4,coc,n) =&
                                             rbp(3)*v(coset(ax,ay,az),1,coc,n) +&
                                             rpw(3)*v(coset(ax,ay,az),1,coc,n+1) +&
                                             fz*(   v(coset(ax,ay,az-1),1,coc,n) +&
                                             f4*v(coset(ax,ay,az-1),1,coc,n+1)) +&
                                             fcz*v(coset(ax,ay,az),1,cocz,n+1)
                                     END IF

                                  END DO
                               END DO
                            END DO

                            ! *** Recurrence steps: [ap||c] -> [ab||c] ***

                            DO lb=2,lb_max

                               ! *** Horizontal recurrence steps ***

                               ! *** [ab||c]{n} = [(a+1i)(b-1i)||c]{n} -    ***
                               ! ***              (Bi - Ai)*[a(b-1i)||c]{n} ***

                               la_start = MAX(0,la_min-1)

                               DO la=la_start,la_max-1
                                  DO n=1,nmax-la-lb-lc
                                     DO ax=0,la
                                        DO ay=0,la-ax
                                           az = la - ax - ay

                                           ! *** Shift of angular momentum component z ***

                                           v(coset(ax,ay,az),coset(0,0,lb),coc,n) =&
                                                v(coset(ax,ay,az+1),&
                                                coset(0,0,lb-1),coc,n) -&
                                                rab(3)*v(coset(ax,ay,az),&
                                                coset(0,0,lb-1),coc,n)

                                           ! *** Shift of angular momentum component y ***

                                           DO by=1,lb
                                              bz = lb - by
                                              v(coset(ax,ay,az),coset(0,by,bz),coc,n) =&
                                                   v(coset(ax,ay+1,az),&
                                                   coset(0,by-1,bz),coc,n) -&
                                                   rab(2)*v(coset(ax,ay,az),&
                                                   coset(0,by-1,bz),coc,n)
                                           END DO

                                           ! *** Shift of angular momentum component x ***

                                           DO bx=1,lb
                                              DO by=0,lb-bx
                                                 bz = lb - bx - by
                                                 v(coset(ax,ay,az),coset(bx,by,bz),coc,n) =&
                                                      v(coset(ax+1,ay,az),&
                                                      coset(bx-1,by,bz),coc,n) -&
                                                      rab(1)*v(coset(ax,ay,az),&
                                                      coset(bx-1,by,bz),coc,n)
                                              END DO
                                           END DO

                                        END DO
                                     END DO
                                  END DO
                               END DO

                               ! *** Vertical recurrence step ***

                               ! *** [ab||c]{n} = (Pi - Bi)*[a(b-1i)||c]{n} +          ***
                               ! ***              (Wi - Pi)*[a(b-1i)||c]{n+1} +        ***
                               ! ***              f2*Ni(a)*(   [(a-1i)(b-1i)||c]{n} +  ***
                               ! ***                        f4*[(a-1i)(b-1i)||c]{n+1}) ***
                               ! ***              f2*Ni(b-1i)*(   [a(b-2i)||c]{n} +    ***
                               ! ***                           f4*[a(b-2i)||c]{n+1}) + ***
                               ! ***              f6*Ni(c)*[a(b-1i)||c-1i]{n+1})       ***

                               DO n=1,nmax-la_max-lb-lc
                                  DO ax=0,la_max
                                     fx = f2*REAL(ax,dp)
                                     DO ay=0,la_max-ax
                                        fy = f2*REAL(ay,dp)
                                        az = la_max - ax - ay
                                        fz = f2*REAL(az,dp)

                                        ! *** Shift of angular momentum component z from a to b ***

                                        f3 = f2*REAL(lb-1,dp)

                                        IF (az == 0) THEN
                                           v(coset(ax,ay,az),coset(0,0,lb),coc,n) =&
                                                rbp(3)*v(coset(ax,ay,az),&
                                                coset(0,0,lb-1),coc,n) +&
                                                rpw(3)*v(coset(ax,ay,az),&
                                                coset(0,0,lb-1),coc,n+1) +&
                                                f3*(   v(coset(ax,ay,az),&
                                                coset(0,0,lb-2),coc,n) +&
                                                f4*v(coset(ax,ay,az),&
                                                coset(0,0,lb-2),coc,n+1)) +&
                                                fcz*v(coset(ax,ay,az),&
                                                coset(0,0,lb-1),cocz,n+1)
                                        ELSE
                                           v(coset(ax,ay,az),coset(0,0,lb),coc,n) =&
                                                rbp(3)*v(coset(ax,ay,az),&
                                                coset(0,0,lb-1),coc,n) +&
                                                rpw(3)*v(coset(ax,ay,az),&
                                                coset(0,0,lb-1),coc,n+1) +&
                                                fz*(   v(coset(ax,ay,az-1),&
                                                coset(0,0,lb-1),coc,n) +&
                                                f4*v(coset(ax,ay,az-1),&
                                                coset(0,0,lb-1),coc,n+1)) +&
                                                f3*(   v(coset(ax,ay,az),&
                                                coset(0,0,lb-2),coc,n) +&
                                                f4*v(coset(ax,ay,az),&
                                                coset(0,0,lb-2),coc,n+1)) +&
                                                fcz*v(coset(ax,ay,az),&
                                                coset(0,0,lb-1),cocz,n+1)
                                        END IF

                                        ! *** Shift of angular momentum component y from a to b ***

                                        IF (ay == 0) THEN
                                           bz = lb - 1
                                           v(coset(ax,ay,az),coset(0,1,bz),coc,n) =&
                                                rbp(2)*v(coset(ax,ay,az),&
                                                coset(0,0,bz),coc,n) +&
                                                rpw(2)*v(coset(ax,ay,az),&
                                                coset(0,0,bz),coc,n+1) +&
                                                fcy*v(coset(ax,ay,az),&
                                                coset(0,0,bz),cocy,n+1)
                                           DO by=2,lb
                                              bz = lb - by
                                              f3 = f2*REAL(by-1,dp)
                                              v(coset(ax,ay,az),coset(0,by,bz),coc,n) =&
                                                   rbp(2)*v(coset(ax,ay,az),&
                                                   coset(0,by-1,bz),coc,n) +&
                                                   rpw(2)*v(coset(ax,ay,az),&
                                                   coset(0,by-1,bz),coc,n+1) +&
                                                   f3*(   v(coset(ax,ay,az),&
                                                   coset(0,by-2,bz),coc,n) +&
                                                   f4*v(coset(ax,ay,az),&
                                                   coset(0,by-2,bz),coc,n+1)) +&
                                                   fcy*v(coset(ax,ay,az),&
                                                   coset(0,by-1,bz),cocy,n+1)
                                           END DO
                                        ELSE
                                           bz = lb - 1
                                           v(coset(ax,ay,az),coset(0,1,bz),coc,n) =&
                                                rbp(2)*v(coset(ax,ay,az),&
                                                coset(0,0,bz),coc,n) +&
                                                rpw(2)*v(coset(ax,ay,az),&
                                                coset(0,0,bz),coc,n+1) +&
                                                fy*(   v(coset(ax,ay-1,az),&
                                                coset(0,0,bz),coc,n) +&
                                                f4*v(coset(ax,ay-1,az),&
                                                coset(0,0,bz),coc,n+1)) +&
                                                fcy*v(coset(ax,ay,az),&
                                                coset(0,0,bz),cocy,n+1)
                                           DO by=2,lb
                                              bz = lb - by
                                              f3 = f2*REAL(by-1,dp)
                                              v(coset(ax,ay,az),coset(0,by,bz),coc,n) =&
                                                   rbp(2)*v(coset(ax,ay,az),&
                                                   coset(0,by-1,bz),coc,n) +&
                                                   rpw(2)*v(coset(ax,ay,az),&
                                                   coset(0,by-1,bz),coc,n+1) +&
                                                   fy*(   v(coset(ax,ay-1,az),&
                                                   coset(0,by-1,bz),coc,n) +&
                                                   f4*v(coset(ax,ay-1,az),&
                                                   coset(0,by-1,bz),coc,n+1)) +&
                                                   f3*(   v(coset(ax,ay,az),&
                                                   coset(0,by-2,bz),coc,n) +&
                                                   f4*v(coset(ax,ay,az),&
                                                   coset(0,by-2,bz),coc,n+1)) +&
                                                   fcy*v(coset(ax,ay,az),&
                                                   coset(0,by-1,bz),cocy,n+1)
                                           END DO
                                        END IF

                                        ! *** Shift of angular momentum component x from a to b ***

                                        IF (ax == 0) THEN
                                           DO by=0,lb-1
                                              bz = lb - 1 - by
                                              v(coset(ax,ay,az),coset(1,by,bz),coc,n) =&
                                                   rbp(1)*v(coset(ax,ay,az),&
                                                   coset(0,by,bz),coc,n) +&
                                                   rpw(1)*v(coset(ax,ay,az),&
                                                   coset(0,by,bz),coc,n+1) +&
                                                   fcx*v(coset(ax,ay,az),&
                                                   coset(0,by,bz),cocx,n+1)
                                           END DO
                                           DO bx=2,lb
                                              f3 = f2*REAL(bx-1,dp)
                                              DO by=0,lb-bx
                                                 bz = lb - bx - by
                                                 v(coset(ax,ay,az),coset(bx,by,bz),coc,n) =&
                                                      rbp(1)*v(coset(ax,ay,az),&
                                                      coset(bx-1,by,bz),coc,n) +&
                                                      rpw(1)*v(coset(ax,ay,az),&
                                                      coset(bx-1,by,bz),coc,n+1) +&
                                                      f3*(   v(coset(ax,ay,az),&
                                                      coset(bx-2,by,bz),coc,n) +&
                                                      f4*v(coset(ax,ay,az),&
                                                      coset(bx-2,by,bz),coc,n+1)) +&
                                                      fcx*v(coset(ax,ay,az),&
                                                      coset(bx-1,by,bz),cocx,n+1)
                                              END DO
                                           END DO
                                        ELSE
                                           DO by=0,lb-1
                                              bz = lb - 1 - by
                                              v(coset(ax,ay,az),coset(1,by,bz),coc,n) =&
                                                   rbp(1)*v(coset(ax,ay,az),&
                                                   coset(0,by,bz),coc,n) +&
                                                   rpw(1)*v(coset(ax,ay,az),&
                                                   coset(0,by,bz),coc,n+1) +&
                                                   fx*(   v(coset(ax-1,ay,az),&
                                                   coset(0,by,bz),coc,n) +&
                                                   f4*v(coset(ax-1,ay,az),&
                                                   coset(0,by,bz),coc,n+1)) +&
                                                   fcx*v(coset(ax,ay,az),&
                                                   coset(0,by,bz),cocx,n+1)
                                           END DO
                                           DO bx=2,lb
                                              f3 = f2*REAL(bx-1,dp)
                                              DO by=0,lb-bx
                                                 bz = lb - bx - by
                                                 v(coset(ax,ay,az),coset(bx,by,bz),coc,n) =&
                                                      rbp(1)*v(coset(ax,ay,az),&
                                                      coset(bx-1,by,bz),coc,n) +&
                                                      rpw(1)*v(coset(ax,ay,az),&
                                                      coset(bx-1,by,bz),coc,n+1) +&
                                                      fx*(   v(coset(ax-1,ay,az),&
                                                      coset(bx-1,by,bz),coc,n) +&
                                                      f4*v(coset(ax-1,ay,az),&
                                                      coset(bx-1,by,bz),coc,n+1)) +&
                                                      f3*(   v(coset(ax,ay,az),&
                                                      coset(bx-2,by,bz),coc,n) +&
                                                      f4*v(coset(ax,ay,az),&
                                                      coset(bx-2,by,bz),coc,n+1)) +&
                                                      fcx*v(coset(ax,ay,az),&
                                                      coset(bx-1,by,bz),cocx,n+1)
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

                            ! *** Vertical recurrence steps: [ss||c] -> [sb||c] ***

                            ! *** [sp||c]{n} = (Pi - Bi)*[ss||c]{n} +    ***
                            ! ***              (Wi - Pi)*[ss||c]{n+1} +  ***
                            ! ***              f6*Ni(c)**[ss||c-1i]{n+1} ***

                            DO n=1,nmax-1-lc
                               v(1,2,coc,n) = rbp(1)*v(1,1,coc,n) +&
                                    rpw(1)*v(1,1,coc,n+1) +&
                                    fcx*v(1,1,cocx,n+1)
                               v(1,3,coc,n) = rbp(2)*v(1,1,coc,n) +&
                                    rpw(2)*v(1,1,coc,n+1) +&
                                    fcy*v(1,1,cocy,n+1)
                               v(1,4,coc,n) = rbp(3)*v(1,1,coc,n) +&
                                    rpw(3)*v(1,1,coc,n+1) +&
                                    fcz*v(1,1,cocz,n+1)
                            END DO

                            ! *** [sb||c]{n} = (Pi - Bi)*[s(b-1i)||c]{n} +          ***
                            ! ***              (Wi - Pi)*[s(b-1i)||c]{n+1} +        ***
                            ! ***              f2*Ni(b-1i)*(   [s(b-2i)||c]{n} +    ***
                            ! ***                           f4*[s(b-2i)||c]{n+1}) + ***
                            ! ***              f6*Ni(c)**[s(b-1i)||c-1i]{n+1}       ***

                            DO lb=2,lb_max

                               DO n=1,nmax-lb-lc

                                  ! *** Increase the angular momentum component z of b ***

                                  v(1,coset(0,0,lb),coc,n) =&
                                       rbp(3)*v(1,coset(0,0,lb-1),coc,n) +&
                                       rpw(3)*v(1,coset(0,0,lb-1),coc,n+1) +&
                                       f2*REAL(lb-1,dp)*(   v(1,coset(0,0,lb-2),coc,n) +&
                                       f4*v(1,coset(0,0,lb-2),coc,n+1)) +&
                                       fcz*v(1,coset(0,0,lb-1),cocz,n+1)

                                  ! *** Increase the angular momentum component y of b ***

                                  bz = lb - 1
                                  v(1,coset(0,1,bz),coc,n) =&
                                       rbp(2)*v(1,coset(0,0,bz),coc,n) +&
                                       rpw(2)*v(1,coset(0,0,bz),coc,n+1) +&
                                       fcy*v(1,coset(0,0,bz),cocy,n+1)

                                  DO by=2,lb
                                     f3 = f2*REAL(by-1,dp)
                                     bz = lb - by
                                     v(1,coset(0,by,bz),coc,n) =&
                                          rbp(2)*v(1,coset(0,by-1,bz),coc,n) +&
                                          rpw(2)*v(1,coset(0,by-1,bz),coc,n+1) +&
                                          f3*(   v(1,coset(0,by-2,bz),coc,n) +&
                                          f4*v(1,coset(0,by-2,bz),coc,n+1)) +&
                                          fcy*v(1,coset(0,by-1,bz),cocy,n+1)
                                  END DO

                                  ! *** Increase the angular momentum component x of b ***

                                  DO by=0,lb-1
                                     bz = lb - 1 - by
                                     v(1,coset(1,by,bz),coc,n) =&
                                          rbp(1)*v(1,coset(0,by,bz),coc,n) +&
                                          rpw(1)*v(1,coset(0,by,bz),coc,n+1) +&
                                          fcx*v(1,coset(0,by,bz),cocx,n+1)
                                  END DO

                                  DO bx=2,lb
                                     f3 = f2*REAL(bx-1,dp)
                                     DO by=0,lb-bx
                                        bz = lb - bx - by
                                        v(1,coset(bx,by,bz),coc,n) =&
                                             rbp(1)*v(1,coset(bx-1,by,bz),coc,n) +&
                                             rpw(1)*v(1,coset(bx-1,by,bz),coc,n+1) +&
                                             f3*(   v(1,coset(bx-2,by,bz),coc,n) +&
                                             f4*v(1,coset(bx-2,by,bz),coc,n+1)) +&
                                             fcx*v(1,coset(bx-1,by,bz),cocx,n+1)
                                     END DO
                                  END DO

                               END DO

                            END DO

                         END IF

                      END IF

                   END DO
                END DO

             END DO

          END IF

          ! *** Add the contribution of the current pair ***
          ! *** of primitive Gaussian-type functions     ***

          DO k=ncoset(lc_min-1)+1,ncoset(lc_max)
             kk = k - ncoset(lc_min-1)
             DO j=ncoset(lb_min-1)+1,ncoset(lb_max)
                DO i=ncoset(la_min-1)+1,ncoset(la_max-maxder_local)
                   vabc(na+i,nb+j) = vabc(na+i,nb+j) + gccc(kk)*v(i,j,k,1)
                   int_abc(na+i,nb+j,kk) = v(i,j,k,1)
                END DO
             END DO
          END DO

          IF(PRESENT(maxder)) THEN
             DO k=ncoset(lc_min-1)+1,ncoset(lc_max)
                kk = k - ncoset(lc_min-1)
                DO j=1,ncoset(lb_max)
                   DO i=1,ncoset(la_max)
                      vabc_plus(nap+i,nb+j) = vabc_plus(nap+i,nb+j) + gccc(kk)*v(i,j,k,1)
                   END DO
                END DO
             END DO
          END IF

          nb = nb + ncoset(lb_max)

       END DO

       na = na + ncoset(la_max - maxder_local)
       nap = nap + ncoset(la_max)

    END DO

  END SUBROUTINE coulomb3

END MODULE ai_coulomb

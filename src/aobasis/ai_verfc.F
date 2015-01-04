!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Build up the nuclear potential part of the core Hamiltonian matrix
!>      in the case of an allelectron calculation by computing the nuclear
!>      attraction integral matrix <a|-Z/r|b> and the integral matrix
!>      <a|Z*erf(r)/r|b>.  The integrals <a|Z*erf(r)/r|b> can be rewritten
!>      as the three-center Coulomb integrals <ab||c> with a  primitive
!>      s-type Gaussian function c multiplied by a factor N.
!>
!>                       erfc(r)
!>      <a|V|b> = -Z*<a|---------|b>
!>                          r
!>
!>                       1 - erf(r)
!>              = -Z*<a|------------|b>
!>                           r
!>
!>                        1           erf(r)
!>              = -Z*(<a|---|b> - <a|--------|b>)
!>                        r             r
!>
!>                        1
!>              = -Z*(<a|---|b> - N*<ab||c>)
!>                        r
!>
!>                       1
!>              = -Z*<a|---|b> + Z*N*<ab||c>
!>                       r
!>                   \_______/       \_____/
!>                       |              |
!>                    nuclear        coulomb
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par Parameters
!>      -  ax,ay,az  : Angular momentum index numbers of orbital a.
!>      -  bx,by,bz  : Angular momentum index numbers of orbital b.
!>      -  coset     : Cartesian orbital set pointer.
!>      -  dab       : Distance between the atomic centers a and b.
!>      -  dac       : Distance between the atomic centers a and c.
!>      -  dbc       : Distance between the atomic centers b and c.
!>      -  l{a,b}    : Angular momentum quantum number of shell a or b.
!>      -  l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
!>      -  l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
!>      -  ncoset    : Number of orbitals in a Cartesian orbital set.
!>      -  npgf{a,b} : Degree of contraction of shell a or b.
!>      -  rab       : Distance vector between the atomic centers a and b.
!>      -  rab2      : Square of the distance between the atomic centers a and b.
!>      -  rac       : Distance vector between the atomic centers a and c.
!>      -  rac2      : Square of the distance between the atomic centers a and c.
!>      -  rbc2      : Square of the distance between the atomic centers b and c.
!>      -  rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
!>      -  zet{a,b}  : Exponents of the Gaussian-type functions a or b.
!>      -  zetp      : Reciprocal of the sum of the exponents of orbital a and b.
!>      -  zetw        : Reciprocal of the sum of the exponents of orbital a, b and c.
!> Indices for pVp matrices
!>      -  ax,ay,az                : Angular momentum index numbers of orbital a.
!>      -  n{x1,y1,z1}{x2,y2,z2}   : indices for the pVp matrix elements
!>                                   <p_{x1,y1,z1}a|V|p_{x2,y2,z2}>
!>      -  co{a,b}{m,p}{x,y,z}     : coefficients of the pVp matrix elements: co == coset(lx,ly,lz)
!>                                   m == li - 1; p == li +1; li= x,y,z (which l quantum number is modified); l = a,b
!>                                   coset ( ... -1 ...) = 1  , should be zero, lmi = 1.0 for li >=1 and lmi = 0.0 for li = 0
!>                                  guarantees this as a prefactor
!>      -  f{a,b}{x,y,z},ft{a,b},z : f(t)li: prefactors as a result of derivations of l with respect to i
!>      -  pVpout                  : pVp matrix elements
!>      -  pVp                     : local pVp matrix elements
!>      -  lamin,lbmin             : minimal angular quantum number used in pVp matrices
!>      _  vnucp                   : potential of the nuclei
!>      -  vnuc                    : potential of the electrons
!>      _  pVpa, pVpb              : pVpl=coset(lx,ly,lz)
!>      -  na_pgf,nb_pgf           : indice for primitive gaussian functions
!> \par History
!>      10.2008 added pVp matrix elements (jens)
!> \author Matthias Krack (04.10.2000)
! *****************************************************************************


MODULE ai_verfc
  
  USE gamma,                           ONLY: fgamma=>fgamma_0
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_verfc'
  PRIVATE

! *** Public subroutines ***

  PUBLIC :: verfc

CONTAINS

! *****************************************************************************
!> \brief   Calculation of the primitive three-center nuclear potential
!>          integrals <a|Z*erfc(r)/r|b> over Cartesian Gaussian-type
!>          functions.
!> \param la_max1 ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min1 ...
!> \param lb_max1 ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min1 ...
!> \param zetc ...
!> \param rpgfc ...
!> \param zc ...
!> \param cerf ...
!> \param rab ...
!> \param rab2 ...
!> \param rac ...
!> \param rac2 ...
!> \param rbc2 ...
!> \param vabc ...
!> \param verf ...
!> \param vnuc ...
!> \param f ...
!> \param maxder ...
!> \param vabc_plus ...
!> \param vnabc ...
!> \param pVp_sum ...
!> \param pVp ...
!> \param dkh_erfc ...
!> \date    04.10.2000
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  SUBROUTINE verfc(la_max1,npgfa,zeta,rpgfa,la_min1,lb_max1,npgfb,zetb,rpgfb,lb_min1,&
       zetc,rpgfc,zc,cerf,rab,rab2,rac,rac2,rbc2,vabc,verf,vnuc,f,&
       maxder,vabc_plus,vnabc,pVp_sum,pVp,dkh_erfc)                                             !JT
    INTEGER, INTENT(IN)                      :: la_max1, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zeta, rpgfa
    INTEGER, INTENT(IN)                      :: la_min1, lb_max1, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zetb, rpgfb
    INTEGER, INTENT(IN)                      :: lb_min1
    REAL(KIND=dp), INTENT(IN)                :: zetc, rpgfc, zc, cerf
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(IN)                :: rab2
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac
    REAL(KIND=dp), INTENT(IN)                :: rac2, rbc2
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: vabc
    REAL(KIND=dp), DIMENSION(:, :, :)        :: verf, vnuc
    REAL(KIND=dp), DIMENSION(0:)             :: f
    INTEGER, INTENT(IN), OPTIONAL            :: maxder
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL :: vabc_plus, vnabc, pVp_sum, pVp
    LOGICAL, OPTIONAL                        :: dkh_erfc

    CHARACTER(len=*), PARAMETER :: routineN = 'verfc', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ay, az, bx, by, bz, coamx, coamy, coamz, coapx, coapy, &
      coapz, cobmx, cobmy, cobmz, cobpx, cobpy, cobpz, i, ipgf, j, jpgf, la, &
      la_max, la_min, la_start, lb, lb_max, lb_min, maxder_local, n, na, nap, &
      nb, nmax, nxx, nxy, nxz, nyx, nyy, nyz, nzx, nzy, nzz, pVpa, pVpb
    LOGICAL                                  :: do_dkh
    REAL(KIND=dp) :: dab, dac, dbc, f0, f1, f2, f3, f4, fax, fay, faz, fbx, &
      fby, fbz, ferf, fnuc, ftaz, ftbz, fx, fy, fz, rcp2, t, zetp, zetq, zetw
    REAL(KIND=dp), DIMENSION(3)              :: rap, rbp, rcp, rpw

!   ---------------------------------------------------------------------------

    IF(PRESENT(maxder)) THEN
      verf(1:ncoset(la_max1),1:ncoset(lb_max1),1) = 0._dp
      vnuc(1:ncoset(la_max1),1:ncoset(lb_max1),1) = 0._dp
    END IF

    do_dkh = .FALSE.

    IF(PRESENT(pVp)) THEN
      do_dkh=.TRUE.
      pVp = 0.0_dp
      ! pVp matrices requires angular momentum quantum numbers l+1 and l-1

      la_max = la_max1+1
      IF(PRESENT(maxder))THEN
        IF(maxder > 0) THEN
          la_max = la_max1
        END IF
      END IF
      lb_max = lb_max1+1

      la_min = MAX(0,la_min1-1)
      lb_min = MAX(0,lb_min1-1)
    ELSE
      do_dkh = .FALSE.
      la_max = la_max1
      lb_max = lb_max1
      la_min = la_min1
      lb_min = lb_min1
    END IF

    nmax = la_max + lb_max + 1

    maxder_local = 0
    IF(PRESENT(maxder)) maxder_local = maxder

!   *** Calculate the distances of the centers a, b and c ***

    dab = SQRT(rab2)
    dac = SQRT(rac2)
    dbc = SQRT(rbc2)

!   *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0
    nap = 0

    DO ipgf=1,npgfa

!     *** Screening ***

      IF(do_dkh)THEN
        IF(dkh_erfc)THEN
            IF (rpgfa(ipgf) + rpgfc < dac) THEN
            na = na + ncoset(la_max1-maxder_local)             !JT
            nap = nap + ncoset(la_max1)                        !JT
            CYCLE
          END IF
        END IF
      ELSE
        IF (rpgfa(ipgf) + rpgfc < dac) THEN
          na = na + ncoset(la_max1-maxder_local)             !JT
          nap = nap + ncoset(la_max1)                        !JT
          CYCLE
        END IF
      END IF

      nb = 0

      DO jpgf=1,npgfb

!       *** Screening ***
       IF(do_dkh)THEN
         IF(dkh_erfc)THEN
           IF ((rpgfb(jpgf) + rpgfc < dbc).OR.&
               (rpgfa(ipgf) + rpgfb(jpgf) < dab)) THEN
             nb = nb + ncoset(lb_max1)                        !JT
             CYCLE
           END IF
         END IF
       ELSE
         IF ((rpgfb(jpgf) + rpgfc < dbc).OR.&
             (rpgfa(ipgf) + rpgfb(jpgf) < dab)) THEN
           nb = nb + ncoset(lb_max1)                        !JT
           CYCLE
         END IF
       END IF
!       *** Calculate some prefactors ***

        zetp = 1.0_dp/(zeta(ipgf) + zetb(jpgf))
        zetq = 1.0_dp/zetc
        zetw = 1.0_dp/(zeta(ipgf) + zetb(jpgf) + zetc)

        f1 = zetb(jpgf)*zetp
        f2 = 0.5_dp*zetp
        f4 = -zetc*zetw

        f0 = EXP(-zeta(ipgf)*f1*rab2)

        fnuc = 2.0_dp*pi*zetp*f0
        ferf = 2.0_dp*SQRT(pi**5*zetw)*zetp*zetq*f0

        rap(:) = f1*rab(:)
        rcp(:) = rap(:) - rac(:)
        rpw(:) = f4*rcp(:)

!       *** Calculate the incomplete Gamma function values ***

        rcp2 = rcp(1)*rcp(1) + rcp(2)*rcp(2) + rcp(3)*rcp(3)

        t = rcp2/zetp

        CALL fgamma(nmax-1,t,f)

!       *** Calculate the basic nuclear attraction integrals [s|A(0)|s]{n} ***

        DO n=1,nmax
          vnuc(1,1,n) = fnuc*f(n-1)
        END DO

!       *** Calculate the incomplete Gamma function values ***

        t = -f4*rcp2/zetp

        CALL fgamma(nmax-1,t,f)

!       *** Calculate the basic three-center Coulomb integrals [ss||s]{n} ***

        DO n=1,nmax
          verf(1,1,n) = ferf*f(n-1)
        END DO

!       *** Recurrence steps: [s|A(0)|s] -> [a|A(0)|b] ***
!       ***                   [ss||s]    -> [ab||s]    ***

        IF (la_max > 0) THEN

!         *** Vertical recurrence steps: [s|A(0)|s] -> [a|A(0)|s] ***
!         ***                            [ss||s]    -> [as||s]    ***

!         *** [p|A(0)|s]{n} = (Pi - Ai)*[s|A(0)|s]{n} -           ***
!         ***                 (Pi - Ci)*[s|A(0)|s]{n+1}           ***
!         *** [ps||s]{n}    = (Pi - Ai)*[ss||s]{n} +              ***
!         ***                 (Wi - Pi)*[ss||s]{n+1}  (i = x,y,z) ***

          DO n=1,nmax-1
            vnuc(2,1,n) = rap(1)*vnuc(1,1,n) - rcp(1)*vnuc(1,1,n+1)
            verf(2,1,n) = rap(1)*verf(1,1,n) + rpw(1)*verf(1,1,n+1)
            vnuc(3,1,n) = rap(2)*vnuc(1,1,n) - rcp(2)*vnuc(1,1,n+1)
            verf(3,1,n) = rap(2)*verf(1,1,n) + rpw(2)*verf(1,1,n+1)
            vnuc(4,1,n) = rap(3)*vnuc(1,1,n) - rcp(3)*vnuc(1,1,n+1)
            verf(4,1,n) = rap(3)*verf(1,1,n) + rpw(3)*verf(1,1,n+1)
          END DO

!         *** [a|A(0)|s]{n} = (Pi - Ai)*[a-1i|A(0)|s]{n} -       ***
!         ***                 (Pi - Ci)*[a-1i|A(0)|s]{n+1} +     ***
!         ***                 f2*Ni(a-1i)*([a-2i|A(0)|s]{n} -    ***
!         ***                              [a-2i|A(0)|s]{n+1}    ***
!         *** [as||s]{n}    = (Pi - Ai)*[(a-1i)s||s]{n} +        ***
!         ***                 (Wi - Pi)*[(a-1i)s||s]{n+1} +      ***
!         ***                 f2*Ni(a-1i)*(   [(a-2i)s||s]{n} +  ***
!         ***                              f4*[(a-2i)s||s]{n+1}) ***

          DO la=2,la_max

            DO n=1,nmax-la

!             *** Increase the angular momentum component z of function a ***

              vnuc(coset(0,0,la),1,n) =&
                rap(3)*vnuc(coset(0,0,la-1),1,n) -&
                rcp(3)*vnuc(coset(0,0,la-1),1,n+1) +&
                f2*REAL(la-1,dp)*(vnuc(coset(0,0,la-2),1,n) -&
                                  vnuc(coset(0,0,la-2),1,n+1))
              verf(coset(0,0,la),1,n) =&
                rap(3)*verf(coset(0,0,la-1),1,n) +&
                rpw(3)*verf(coset(0,0,la-1),1,n+1) +&
                f2*REAL(la-1,dp)*(   verf(coset(0,0,la-2),1,n) +&
                                  f4*verf(coset(0,0,la-2),1,n+1))

!             *** Increase the angular momentum component y of function a ***

              az = la - 1
              vnuc(coset(0,1,az),1,n) =&
                rap(2)*vnuc(coset(0,0,az),1,n) -&
                rcp(2)*vnuc(coset(0,0,az),1,n+1)
              verf(coset(0,1,az),1,n) =&
                rap(2)*verf(coset(0,0,az),1,n) +&
                rpw(2)*verf(coset(0,0,az),1,n+1)

              DO ay=2,la
                az = la - ay
                vnuc(coset(0,ay,az),1,n) =&
                  rap(2)*vnuc(coset(0,ay-1,az),1,n) -&
                  rcp(2)*vnuc(coset(0,ay-1,az),1,n+1) +&
                  f2*REAL(ay-1,dp)*(vnuc(coset(0,ay-2,az),1,n) -&
                                    vnuc(coset(0,ay-2,az),1,n+1))
                verf(coset(0,ay,az),1,n) =&
                  rap(2)*verf(coset(0,ay-1,az),1,n) +&
                  rpw(2)*verf(coset(0,ay-1,az),1,n+1) +&
                  f2*REAL(ay-1,dp)*(   verf(coset(0,ay-2,az),1,n) +&
                                    f4*verf(coset(0,ay-2,az),1,n+1))
              END DO

!             *** Increase the angular momentum component x of function a ***

              DO ay=0,la-1
                az = la - 1 - ay
                vnuc(coset(1,ay,az),1,n) =&
                  rap(1)*vnuc(coset(0,ay,az),1,n) -&
                  rcp(1)*vnuc(coset(0,ay,az),1,n+1)
                verf(coset(1,ay,az),1,n) =&
                  rap(1)*verf(coset(0,ay,az),1,n) +&
                  rpw(1)*verf(coset(0,ay,az),1,n+1)
              END DO

              DO ax=2,la
                f3 = f2*REAL(ax-1,dp)
                DO ay=0,la-ax
                  az = la - ax - ay
                  vnuc(coset(ax,ay,az),1,n) =&
                    rap(1)*vnuc(coset(ax-1,ay,az),1,n) -&
                    rcp(1)*vnuc(coset(ax-1,ay,az),1,n+1) +&
                    f3*(vnuc(coset(ax-2,ay,az),1,n) -&
                        vnuc(coset(ax-2,ay,az),1,n+1))
                  verf(coset(ax,ay,az),1,n) =&
                    rap(1)*verf(coset(ax-1,ay,az),1,n) +&
                    rpw(1)*verf(coset(ax-1,ay,az),1,n+1) +&
                    f3*(   verf(coset(ax-2,ay,az),1,n) +&
                        f4*verf(coset(ax-2,ay,az),1,n+1))
                END DO
              END DO

            END DO

          END DO

!         *** Recurrence steps: [a|A(0)|s] -> [a|A(0)|b] ***
!         ***                   [as||s]    -> [ab||s]    ***

          IF (lb_max > 0) THEN

!           *** Horizontal recurrence steps ***

            rbp(:) = rap(:) - rab(:)

!           *** [a||A(0)|p]{n} = [a+1i|A(0)|s]{n} - (Bi - Ai)*[a|A(0)|s]{n} ***
!           *** [ap||s]{n}     = [(a+1i)s||s]{n}  - (Bi - Ai)*[as||s]{n}    ***

            la_start = MAX(0,la_min-1)

            DO la=la_start,la_max-1
              DO n=1,nmax-la-1
                DO ax=0,la
                  DO ay=0,la-ax
                    az = la - ax - ay
                    vnuc(coset(ax,ay,az),2,n) =&
                      vnuc(coset(ax+1,ay,az),1,n) -&
                      rab(1)*vnuc(coset(ax,ay,az),1,n)
                    verf(coset(ax,ay,az),2,n) =&
                      verf(coset(ax+1,ay,az),1,n) -&
                      rab(1)*verf(coset(ax,ay,az),1,n)
                    vnuc(coset(ax,ay,az),3,n) =&
                      vnuc(coset(ax,ay+1,az),1,n) -&
                      rab(2)*vnuc(coset(ax,ay,az),1,n)
                    verf(coset(ax,ay,az),3,n) =&
                      verf(coset(ax,ay+1,az),1,n) -&
                      rab(2)*verf(coset(ax,ay,az),1,n)
                    vnuc(coset(ax,ay,az),4,n) =&
                      vnuc(coset(ax,ay,az+1),1,n) -&
                      rab(3)*vnuc(coset(ax,ay,az),1,n)
                    verf(coset(ax,ay,az),4,n) =&
                      verf(coset(ax,ay,az+1),1,n) -&
                      rab(3)*verf(coset(ax,ay,az),1,n)
                  END DO
                END DO
              END DO
            END DO

!           *** Vertical recurrence step ***

!           *** [a|A(0)|p]{n} = (Pi - Bi)*[a|A(0)|s]{n} -        ***
!           ***                 (Pi - Ci)*[a|A(0)|s]{n+1} +      ***
!           ***                 f2*Ni(a)*([a-1i|A(0)|s]{n} -     ***
!           ***                           [a-1i|A(0)|s]{n+1})    ***
!           *** [ap||s]{n}    = (Pi - Bi)*[as||s]{n} +           ***
!           ***                 (Wi - Pi)*[as||s]{n+1} +         ***
!           ***                  f2*Ni(a)*(   [(a-1i)s||s]{n} +  ***
!           ***                            f4*[(a-1i)s||s]{n+1}) ***

            DO n=1,nmax-la_max-1
              DO ax=0,la_max
                fx = f2*REAL(ax,dp)
                DO ay=0,la_max-ax
                  fy = f2*REAL(ay,dp)
                  az = la_max - ax - ay
                  fz = f2*REAL(az,dp)

                  IF (ax == 0) THEN
                    vnuc(coset(ax,ay,az),2,n) =&
                      rbp(1)*vnuc(coset(ax,ay,az),1,n) -&
                      rcp(1)*vnuc(coset(ax,ay,az),1,n+1)
                    verf(coset(ax,ay,az),2,n) =&
                      rbp(1)*verf(coset(ax,ay,az),1,n) +&
                      rpw(1)*verf(coset(ax,ay,az),1,n+1)
                  ELSE
                    vnuc(coset(ax,ay,az),2,n) =&
                      rbp(1)*vnuc(coset(ax,ay,az),1,n) -&
                      rcp(1)*vnuc(coset(ax,ay,az),1,n+1) +&
                      fx*(vnuc(coset(ax-1,ay,az),1,n) -&
                          vnuc(coset(ax-1,ay,az),1,n+1))
                    verf(coset(ax,ay,az),2,n) =&
                      rbp(1)*verf(coset(ax,ay,az),1,n) +&
                      rpw(1)*verf(coset(ax,ay,az),1,n+1) +&
                      fx*(   verf(coset(ax-1,ay,az),1,n) +&
                          f4*verf(coset(ax-1,ay,az),1,n+1))
                  END IF

                  IF (ay == 0) THEN
                    vnuc(coset(ax,ay,az),3,n) =&
                      rbp(2)*vnuc(coset(ax,ay,az),1,n) -&
                      rcp(2)*vnuc(coset(ax,ay,az),1,n+1)
                    verf(coset(ax,ay,az),3,n) =&
                      rbp(2)*verf(coset(ax,ay,az),1,n) +&
                      rpw(2)*verf(coset(ax,ay,az),1,n+1)
                  ELSE
                    vnuc(coset(ax,ay,az),3,n) =&
                      rbp(2)*vnuc(coset(ax,ay,az),1,n) -&
                      rcp(2)*vnuc(coset(ax,ay,az),1,n+1) +&
                      fy*(vnuc(coset(ax,ay-1,az),1,n) -&
                          vnuc(coset(ax,ay-1,az),1,n+1))
                    verf(coset(ax,ay,az),3,n) =&
                      rbp(2)*verf(coset(ax,ay,az),1,n) +&
                      rpw(2)*verf(coset(ax,ay,az),1,n+1) +&
                      fy*(   verf(coset(ax,ay-1,az),1,n) +&
                          f4*verf(coset(ax,ay-1,az),1,n+1))
                  END IF

                  IF (az == 0) THEN
                    vnuc(coset(ax,ay,az),4,n) =&
                      rbp(3)*vnuc(coset(ax,ay,az),1,n) -&
                      rcp(3)*vnuc(coset(ax,ay,az),1,n+1)
                    verf(coset(ax,ay,az),4,n) =&
                      rbp(3)*verf(coset(ax,ay,az),1,n) +&
                      rpw(3)*verf(coset(ax,ay,az),1,n+1)
                  ELSE
                    vnuc(coset(ax,ay,az),4,n) =&
                      rbp(3)*vnuc(coset(ax,ay,az),1,n) -&
                      rcp(3)*vnuc(coset(ax,ay,az),1,n+1) +&
                      fz*(vnuc(coset(ax,ay,az-1),1,n) -&
                          vnuc(coset(ax,ay,az-1),1,n+1))
                    verf(coset(ax,ay,az),4,n) =&
                      rbp(3)*verf(coset(ax,ay,az),1,n) +&
                      rpw(3)*verf(coset(ax,ay,az),1,n+1) +&
                      fz*(   verf(coset(ax,ay,az-1),1,n) +&
                          f4*verf(coset(ax,ay,az-1),1,n+1))
                  END IF

                END DO
              END DO
            END DO

!           *** Recurrence steps: [a|A(0)|p] -> [a|A(0)|b] ***
!           ***                   [ap||s]    -> [ab||s]    ***

            DO lb=2,lb_max

!             *** Horizontal recurrence steps ***

!             *** [a||A(0)|b]{n} = [a+1i|A(0)|b-1i]{n} -      ***
!             ***                  (Bi - Ai)*[a|A(0)|b-1i]{n} ***
!             *** [ab||s]{n}     = [(a+1i)(b-1i)||s]{n} -     ***
!             ***                  (Bi - Ai)*[a(b-1i)||s]{n}  ***

              la_start = MAX(0,la_min-1)

              DO la=la_start,la_max-1
                DO n=1,nmax-la-lb
                  DO ax=0,la
                    DO ay=0,la-ax
                      az = la - ax - ay

!                     *** Shift of angular momentum component z from a to b ***

                      vnuc(coset(ax,ay,az),coset(0,0,lb),n) =&
                        vnuc(coset(ax,ay,az+1),coset(0,0,lb-1),n) -&
                        rab(3)*vnuc(coset(ax,ay,az),coset(0,0,lb-1),n)
                      verf(coset(ax,ay,az),coset(0,0,lb),n) =&
                        verf(coset(ax,ay,az+1),coset(0,0,lb-1),n) -&
                        rab(3)*verf(coset(ax,ay,az),coset(0,0,lb-1),n)

!                     *** Shift of angular momentum component y from a to b ***

                      DO by=1,lb
                        bz = lb - by
                        vnuc(coset(ax,ay,az),coset(0,by,bz),n) =&
                          vnuc(coset(ax,ay+1,az),coset(0,by-1,bz),n) -&
                          rab(2)*vnuc(coset(ax,ay,az),coset(0,by-1,bz),n)
                        verf(coset(ax,ay,az),coset(0,by,bz),n) =&
                          verf(coset(ax,ay+1,az),coset(0,by-1,bz),n) -&
                          rab(2)*verf(coset(ax,ay,az),coset(0,by-1,bz),n)
                      END DO

!                     *** Shift of angular momentum component x from a to b ***

                      DO bx=1,lb
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          vnuc(coset(ax,ay,az),coset(bx,by,bz),n) =&
                            vnuc(coset(ax+1,ay,az),coset(bx-1,by,bz),n) -&
                            rab(1)*vnuc(coset(ax,ay,az),coset(bx-1,by,bz),n)
                          verf(coset(ax,ay,az),coset(bx,by,bz),n) =&
                            verf(coset(ax+1,ay,az),coset(bx-1,by,bz),n) -&
                            rab(1)*verf(coset(ax,ay,az),coset(bx-1,by,bz),n)
                        END DO
                      END DO

                    END DO
                  END DO
                END DO
              END DO

!             *** Vertical recurrence step ***

!             *** [a|A(0)|b]{n} = (Pi - Bi)*[a|A(0)|b-1i]{n} -         ***
!             ***                 (Pi - Ci)*[a|A(0)|b-1i]{n+1} +       ***
!             ***                 f2*Ni(a)*([a-1i|A(0)|b-1i]{n} -      ***
!             ***                           [a-1i|A(0)|b-1i]{n+1}) +   ***
!             ***                 f2*Ni(b-1i)*([a|A(0)|b-2i]{n} -      ***
!             ***                              [a|A(0)|b-2i]{n+1})     ***
!             *** [ab||s]{n}    = (Pi - Bi)*[a(b-1i)||s]{n} +          ***
!             ***                 (Wi - Pi)*[a(b-1i)||s]{n+1} +        ***
!             ***                 f2*Ni(a)*(   [(a-1i)(b-1i)||s]{n} +  ***
!             ***                           f4*[(a-1i)(b-1i)||s]{n+1}) ***
!             ***                 f2*Ni(b-1i)*(   [a(b-2i)||s]{n} +    ***
!             ***                              f4*[a(b-2i)||s]{n+1})   ***

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
                      vnuc(coset(ax,ay,az),coset(0,0,lb),n) =&
                        rbp(3)*vnuc(coset(ax,ay,az),coset(0,0,lb-1),n) -&
                        rcp(3)*vnuc(coset(ax,ay,az),coset(0,0,lb-1),n+1) +&
                        f3*(vnuc(coset(ax,ay,az),coset(0,0,lb-2),n) -&
                            vnuc(coset(ax,ay,az),coset(0,0,lb-2),n+1))
                      verf(coset(ax,ay,az),coset(0,0,lb),n) =&
                        rbp(3)*verf(coset(ax,ay,az),coset(0,0,lb-1),n) +&
                        rpw(3)*verf(coset(ax,ay,az),coset(0,0,lb-1),n+1) +&
                        f3*(   verf(coset(ax,ay,az),coset(0,0,lb-2),n) +&
                            f4*verf(coset(ax,ay,az),coset(0,0,lb-2),n+1))
                    ELSE
                      vnuc(coset(ax,ay,az),coset(0,0,lb),n) =&
                        rbp(3)*vnuc(coset(ax,ay,az),coset(0,0,lb-1),n) -&
                        rcp(3)*vnuc(coset(ax,ay,az),coset(0,0,lb-1),n+1) +&
                        fz*(vnuc(coset(ax,ay,az-1),coset(0,0,lb-1),n) -&
                            vnuc(coset(ax,ay,az-1),coset(0,0,lb-1),n+1)) +&
                        f3*(vnuc(coset(ax,ay,az),coset(0,0,lb-2),n) -&
                            vnuc(coset(ax,ay,az),coset(0,0,lb-2),n+1))
                      verf(coset(ax,ay,az),coset(0,0,lb),n) =&
                        rbp(3)*verf(coset(ax,ay,az),coset(0,0,lb-1),n) +&
                        rpw(3)*verf(coset(ax,ay,az),coset(0,0,lb-1),n+1) +&
                        fz*(   verf(coset(ax,ay,az-1),coset(0,0,lb-1),n) +&
                            f4*verf(coset(ax,ay,az-1),coset(0,0,lb-1),n+1)) +&
                        f3*(   verf(coset(ax,ay,az),coset(0,0,lb-2),n) +&
                            f4*verf(coset(ax,ay,az),coset(0,0,lb-2),n+1))
                    END IF

!                   *** Shift of angular momentum component y from a to b ***

                    IF (ay == 0) THEN
                      bz = lb - 1
                      vnuc(coset(ax,ay,az),coset(0,1,bz),n) =&
                        rbp(2)*vnuc(coset(ax,ay,az),coset(0,0,bz),n) -&
                        rcp(2)*vnuc(coset(ax,ay,az),coset(0,0,bz),n+1)
                      verf(coset(ax,ay,az),coset(0,1,bz),n) =&
                        rbp(2)*verf(coset(ax,ay,az),coset(0,0,bz),n) +&
                        rpw(2)*verf(coset(ax,ay,az),coset(0,0,bz),n+1)
                      DO by=2,lb
                        bz = lb - by
                        f3 = f2*REAL(by-1,dp)
                        vnuc(coset(ax,ay,az),coset(0,by,bz),n) =&
                          rbp(2)*vnuc(coset(ax,ay,az),coset(0,by-1,bz),n) -&
                          rcp(2)*vnuc(coset(ax,ay,az),coset(0,by-1,bz),n+1) +&
                          f3*(vnuc(coset(ax,ay,az),coset(0,by-2,bz),n) -&
                              vnuc(coset(ax,ay,az),coset(0,by-2,bz),n+1))
                        verf(coset(ax,ay,az),coset(0,by,bz),n) =&
                          rbp(2)*verf(coset(ax,ay,az),coset(0,by-1,bz),n) +&
                          rpw(2)*verf(coset(ax,ay,az),coset(0,by-1,bz),n+1) +&
                          f3*(   verf(coset(ax,ay,az),coset(0,by-2,bz),n) +&
                              f4*verf(coset(ax,ay,az),coset(0,by-2,bz),n+1))
                      END DO
                    ELSE
                      bz = lb - 1
                      vnuc(coset(ax,ay,az),coset(0,1,bz),n) =&
                        rbp(2)*vnuc(coset(ax,ay,az),coset(0,0,bz),n) -&
                        rcp(2)*vnuc(coset(ax,ay,az),coset(0,0,bz),n+1) +&
                        fy*(vnuc(coset(ax,ay-1,az),coset(0,0,bz),n) -&
                            vnuc(coset(ax,ay-1,az),coset(0,0,bz),n+1))
                      verf(coset(ax,ay,az),coset(0,1,bz),n) =&
                        rbp(2)*verf(coset(ax,ay,az),coset(0,0,bz),n) +&
                        rpw(2)*verf(coset(ax,ay,az),coset(0,0,bz),n+1) +&
                        fy*(   verf(coset(ax,ay-1,az),coset(0,0,bz),n) +&
                            f4*verf(coset(ax,ay-1,az),coset(0,0,bz),n+1))
                      DO by=2,lb
                        bz = lb - by
                        f3 = f2*REAL(by-1,dp)
                        vnuc(coset(ax,ay,az),coset(0,by,bz),n) =&
                          rbp(2)*vnuc(coset(ax,ay,az),coset(0,by-1,bz),n) -&
                          rcp(2)*vnuc(coset(ax,ay,az),coset(0,by-1,bz),n+1) +&
                          fy*(vnuc(coset(ax,ay-1,az),coset(0,by-1,bz),n) -&
                              vnuc(coset(ax,ay-1,az),coset(0,by-1,bz),n+1)) +&
                          f3*(vnuc(coset(ax,ay,az),coset(0,by-2,bz),n) -&
                              vnuc(coset(ax,ay,az),coset(0,by-2,bz),n+1))
                        verf(coset(ax,ay,az),coset(0,by,bz),n) =&
                          rbp(2)*verf(coset(ax,ay,az),coset(0,by-1,bz),n) +&
                          rpw(2)*verf(coset(ax,ay,az),coset(0,by-1,bz),n+1) +&
                          fy*(   verf(coset(ax,ay-1,az),coset(0,by-1,bz),n) +&
                              f4*verf(coset(ax,ay-1,az),&
                                      coset(0,by-1,bz),n+1)) +&
                          f3*(   verf(coset(ax,ay,az),coset(0,by-2,bz),n) +&
                              f4*verf(coset(ax,ay,az),coset(0,by-2,bz),n+1))
                      END DO
                    END IF

!                   *** Shift of angular momentum component x from a to b ***

                    IF (ax == 0) THEN
                      DO by=0,lb-1
                        bz = lb - 1 - by
                        vnuc(coset(ax,ay,az),coset(1,by,bz),n) =&
                          rbp(1)*vnuc(coset(ax,ay,az),coset(0,by,bz),n) -&
                          rcp(1)*vnuc(coset(ax,ay,az),coset(0,by,bz),n+1)
                        verf(coset(ax,ay,az),coset(1,by,bz),n) =&
                          rbp(1)*verf(coset(ax,ay,az),coset(0,by,bz),n) +&
                          rpw(1)*verf(coset(ax,ay,az),coset(0,by,bz),n+1)
                      END DO
                      DO bx=2,lb
                        f3 = f2*REAL(bx-1,dp)
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          vnuc(coset(ax,ay,az),coset(bx,by,bz),n) =&
                            rbp(1)*vnuc(coset(ax,ay,az),coset(bx-1,by,bz),n) -&
                            rcp(1)*vnuc(coset(ax,ay,az),&
                                        coset(bx-1,by,bz),n+1) +&
                            f3*(vnuc(coset(ax,ay,az),coset(bx-2,by,bz),n) -&
                                vnuc(coset(ax,ay,az),coset(bx-2,by,bz),n+1))
                          verf(coset(ax,ay,az),coset(bx,by,bz),n) =&
                            rbp(1)*verf(coset(ax,ay,az),coset(bx-1,by,bz),n) +&
                            rpw(1)*verf(coset(ax,ay,az),&
                                        coset(bx-1,by,bz),n+1) +&
                            f3*(   verf(coset(ax,ay,az),coset(bx-2,by,bz),n) +&
                                f4*verf(coset(ax,ay,az),coset(bx-2,by,bz),n+1))
                        END DO
                      END DO
                    ELSE
                      DO by=0,lb-1
                        bz = lb - 1 - by
                        vnuc(coset(ax,ay,az),coset(1,by,bz),n) =&
                          rbp(1)*vnuc(coset(ax,ay,az),coset(0,by,bz),n) -&
                          rcp(1)*vnuc(coset(ax,ay,az),coset(0,by,bz),n+1) +&
                          fx*(vnuc(coset(ax-1,ay,az),coset(0,by,bz),n) -&
                              vnuc(coset(ax-1,ay,az),coset(0,by,bz),n+1))
                        verf(coset(ax,ay,az),coset(1,by,bz),n) =&
                          rbp(1)*verf(coset(ax,ay,az),coset(0,by,bz),n) +&
                          rpw(1)*verf(coset(ax,ay,az),coset(0,by,bz),n+1) +&
                          fx*(   verf(coset(ax-1,ay,az),coset(0,by,bz),n) +&
                              f4*verf(coset(ax-1,ay,az),coset(0,by,bz),n+1))
                      END DO
                      DO bx=2,lb
                        f3 = f2*REAL(bx-1,dp)
                        DO by=0,lb-bx
                          bz = lb - bx - by
                          vnuc(coset(ax,ay,az),coset(bx,by,bz),n) =&
                            rbp(1)*vnuc(coset(ax,ay,az),coset(bx-1,by,bz),n) -&
                            rcp(1)*vnuc(coset(ax,ay,az),&
                                        coset(bx-1,by,bz),n+1) +&
                            fx*(vnuc(coset(ax-1,ay,az),coset(bx-1,by,bz),n) -&
                                vnuc(coset(ax-1,ay,az),&
                                     coset(bx-1,by,bz),n+1)) +&
                            f3*(vnuc(coset(ax,ay,az),coset(bx-2,by,bz),n) -&
                                vnuc(coset(ax,ay,az),coset(bx-2,by,bz),n+1))
                          verf(coset(ax,ay,az),coset(bx,by,bz),n) =&
                            rbp(1)*verf(coset(ax,ay,az),coset(bx-1,by,bz),n) +&
                            rpw(1)*verf(coset(ax,ay,az),&
                                        coset(bx-1,by,bz),n+1) +&
                            fx*(   verf(coset(ax-1,ay,az),&
                                        coset(bx-1,by,bz),n) +&
                                f4*verf(coset(ax-1,ay,az),&
                                        coset(bx-1,by,bz),n+1)) +&
                            f3*(   verf(coset(ax,ay,az),coset(bx-2,by,bz),n) +&
                                f4*verf(coset(ax,ay,az),coset(bx-2,by,bz),n+1))
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
!           ***                            [ss||s]    -> [sb||s]    ***

            rbp(:) = rap(:) - rab(:)

!           *** [s|A(0)|p]{n} = (Pi - Bi)*[s|A(0)|s]{n} - ***
!           ***                 (Pi - Ci)*[s|A(0)|s]{n+1} ***
!           *** [sp||s]{n}    = (Pi - Bi)*[ss||s]{n} +    ***
!           ***                 (Wi - Pi)*[ss||s]{n+1}    ***

            DO n=1,nmax-1
              vnuc(1,2,n) = rbp(1)*vnuc(1,1,n) - rcp(1)*vnuc(1,1,n+1)
              verf(1,2,n) = rbp(1)*verf(1,1,n) + rpw(1)*verf(1,1,n+1)
              vnuc(1,3,n) = rbp(2)*vnuc(1,1,n) - rcp(2)*vnuc(1,1,n+1)
              verf(1,3,n) = rbp(2)*verf(1,1,n) + rpw(2)*verf(1,1,n+1)
              vnuc(1,4,n) = rbp(3)*vnuc(1,1,n) - rcp(3)*vnuc(1,1,n+1)
              verf(1,4,n) = rbp(3)*verf(1,1,n) + rpw(3)*verf(1,1,n+1)
            END DO

!           *** [s|A(0)|b]{n} = (Pi - Bi)*[s|A(0)|b-1i]{n} -       ***
!           ***                 (Pi - Ci)*[s|A(0)|b-1i]{n+1} +     ***
!           ***                 f2*Ni(b-1i)*([s|A(0)|b-2i]{n} -    ***
!           ***                              [s|A(0)|b-2i]{n+1}    ***
!           *** [sb||s]{n}    = (Pi - Bi)*[s(b-1i)||s]{n} +        ***
!           ***                 (Wi - Pi)*[s(b-1i)||s]{n+1} +      ***
!           ***                 f2*Ni(b-1i)*(   [s(b-2i)||s]{n} +  ***
!           ***                              f4*[s(b-2i)||s]{n+1}) ***

            DO lb=2,lb_max

              DO n=1,nmax-lb

!               *** Increase the angular momentum component z of function b ***

                vnuc(1,coset(0,0,lb),n) =&
                  rbp(3)*vnuc(1,coset(0,0,lb-1),n) -&
                  rcp(3)*vnuc(1,coset(0,0,lb-1),n+1) +&
                  f2*REAL(lb-1,dp)*(vnuc(1,coset(0,0,lb-2),n) -&
                                    vnuc(1,coset(0,0,lb-2),n+1))
                verf(1,coset(0,0,lb),n) =&
                  rbp(3)*verf(1,coset(0,0,lb-1),n) +&
                  rpw(3)*verf(1,coset(0,0,lb-1),n+1) +&
                  f2*REAL(lb-1,dp)*(   verf(1,coset(0,0,lb-2),n) +&
                                    f4*verf(1,coset(0,0,lb-2),n+1))

!               *** Increase the angular momentum component y of function b ***

                bz = lb - 1
                vnuc(1,coset(0,1,bz),n) =&
                  rbp(2)*vnuc(1,coset(0,0,bz),n) -&
                  rcp(2)*vnuc(1,coset(0,0,bz),n+1)
                verf(1,coset(0,1,bz),n) =&
                  rbp(2)*verf(1,coset(0,0,bz),n) +&
                  rpw(2)*verf(1,coset(0,0,bz),n+1)

                DO by=2,lb
                  bz = lb - by
                  vnuc(1,coset(0,by,bz),n) =&
                    rbp(2)*vnuc(1,coset(0,by-1,bz),n) -&
                    rcp(2)*vnuc(1,coset(0,by-1,bz),n+1) +&
                    f2*REAL(by-1,dp)*(vnuc(1,coset(0,by-2,bz),n) -&
                                      vnuc(1,coset(0,by-2,bz),n+1))
                  verf(1,coset(0,by,bz),n) =&
                    rbp(2)*verf(1,coset(0,by-1,bz),n) +&
                    rpw(2)*verf(1,coset(0,by-1,bz),n+1) +&
                    f2*REAL(by-1,dp)*(   verf(1,coset(0,by-2,bz),n) +&
                                      f4*verf(1,coset(0,by-2,bz),n+1))
                END DO

!               *** Increase the angular momentum component x of function b ***

                DO by=0,lb-1
                  bz = lb - 1 - by
                  vnuc(1,coset(1,by,bz),n) =&
                    rbp(1)*vnuc(1,coset(0,by,bz),n) -&
                    rcp(1)*vnuc(1,coset(0,by,bz),n+1)
                  verf(1,coset(1,by,bz),n) =&
                    rbp(1)*verf(1,coset(0,by,bz),n) +&
                    rpw(1)*verf(1,coset(0,by,bz),n+1)
                END DO

                DO bx=2,lb
                  f3 = f2*REAL(bx-1,dp)
                  DO by=0,lb-bx
                    bz = lb - bx - by
                    vnuc(1,coset(bx,by,bz),n) =&
                      rbp(1)*vnuc(1,coset(bx-1,by,bz),n) -&
                      rcp(1)*vnuc(1,coset(bx-1,by,bz),n+1) +&
                      f3*(vnuc(1,coset(bx-2,by,bz),n) -&
                          vnuc(1,coset(bx-2,by,bz),n+1))
                    verf(1,coset(bx,by,bz),n) =&
                      rbp(1)*verf(1,coset(bx-1,by,bz),n) +&
                      rpw(1)*verf(1,coset(bx-1,by,bz),n+1) +&
                      f3*(   verf(1,coset(bx-2,by,bz),n) +&
                          f4*verf(1,coset(bx-2,by,bz),n+1))
                  END DO
                END DO

              END DO

    END DO

  END IF

END IF

!       *** Add the contribution of the current pair ***
!       *** of primitive Gaussian-type functions     ***

!JT
        IF(do_dkh)THEN
          DO j=1,ncoset(lb_max1)
            DO i=ncoset(la_min1-1)+1,ncoset(la_max1-maxder_local)
              vnabc(na+i,nb+j) = vnabc(na+i,nb+j) + cerf*verf(i,j,1)
            END DO
          END DO
          DO j=1,ncoset(lb_max1)
            DO i=ncoset(la_min1-1)+1,ncoset(la_max1-maxder_local)
              vabc(na+i,nb+j) = vabc(na+i,nb+j) - zc*vnuc(i,j,1)
            END DO
          END DO
        ELSE
          DO j=1,ncoset(lb_max1)                                        !JT
            DO i=ncoset(la_min1-1)+1,ncoset(la_max1-maxder_local)       !JT
              vabc(na+i,nb+j) = vabc(na+i,nb+j) -&
                                zc*vnuc(i,j,1) + cerf*verf(i,j,1)
            END DO
          END DO
        END IF
!JT

        IF(PRESENT(maxder)) THEN
          DO j=1,ncoset(lb_max1)                                      !JT
            DO i=1,ncoset(la_max1)                                    !JT
              vabc_plus(nap+i,nb+j) = vabc_plus(nap+i,nb+j) -&
                              zc*vnuc(i,j,1) + cerf*verf(i,j,1)
            END DO
          END DO
        END IF
!JTs
        IF(do_dkh)THEN
          ! *********** Calculate pVp matrix **************
          ! [pa|V|pb]=4*zeta*zetb*[a+i|V|b+i]-2*zetb*Ni(a)*[a-i|V|b+i]-2*zeta*Ni(b)*[a+i|V|b-i]+Ni(a)*Ni(b)*[a-i|V|b-i]
          ! Every integral has been calculated before (except [-1|V|-1],[a|V|-1] and [-1|V|b],
          ! which do not contribute, because their prefactor Ni(a) or Ni(b) is zero)
          ! and is given by verf/vnuc(coset(ax,ay,az),coset(bx,by,bz),1)
          ! ***********************************************
          ftaz = 2.0_dp*zeta(ipgf)
          ftbz = 2.0_dp*zetb(jpgf)
          nxx=1
          nyy=2
          nzz=3
          nxy=4
          nxz=5
          nyx=6
          nyz=7
          nzx=8
          nzy=9
          DO la=0,la_max1
            DO ax=0,la
              fax = REAL(ax,dp)
              DO ay=0,la-ax
                fay = REAL(ay,dp)
                az = la - ax - ay
                faz = REAL(az,dp)
                pVpa  = coset(ax,ay,az)
                coamx = coset(ax-1,ay,az)
                coamy = coset(ax,ay-1,az)
                coamz = coset(ax,ay,az-1)
                coapx = coset(ax+1,ay,az)
                coapy = coset(ax,ay+1,az)
                coapz = coset(ax,ay,az+1)
                DO lb=0,lb_max1
                  DO bx=0,lb
                    fbx=REAL(bx,dp)
                    DO by=0,lb-bx
                      fby=REAL(by,dp)
                      bz = lb - bx - by
                      fbz=REAL(bz,dp)
                      pVpb = coset(bx,by,bz)
                      cobmx = coset(bx-1,by,bz)
                      cobmy = coset(bx,by-1,bz)
                      cobmz = coset(bx,by,bz-1)
                      cobpx = coset(bx+1,by,bz)
                      cobpy = coset(bx,by+1,bz)
                      cobpz = coset(bx,by,bz+1)
                      IF(dkh_erfc)THEN
                        pVp(pVpa,pVpb)=ftaz*ftbz*    (-zc*vnuc(coapx,cobpx,1)+cerf*verf(coapx,cobpx,1))-&
                                           ftaz* fbx*(-zc*vnuc(coapx,cobmx,1)+cerf*verf(coapx,cobmx,1))-&
                                           ftbz* fax*(-zc*vnuc(coamx,cobpx,1)+cerf*verf(coamx,cobpx,1))+&
                                           fax*  fbx*(-zc*vnuc(coamx,cobmx,1)+cerf*verf(coamx,cobmx,1))+&
                                          ftaz*ftbz* (-zc*vnuc(coapy,cobpy,1)+cerf*verf(coapy,cobpy,1))-&
                                           ftaz* fby*(-zc*vnuc(coapy,cobmy,1)+cerf*verf(coapy,cobmy,1))-&
                                           ftbz* fay*(-zc*vnuc(coamy,cobpy,1)+cerf*verf(coamy,cobpy,1))+&
                                           fay*  fby*(-zc*vnuc(coamy,cobmy,1)+cerf*verf(coamy,cobmy,1))+&
                                          ftaz*ftbz* (-zc*vnuc(coapz,cobpz,1)+cerf*verf(coapz,cobpz,1))-&
                                           ftaz* fbz*(-zc*vnuc(coapz,cobmz,1)+cerf*verf(coapz,cobmz,1))-&
                                           ftbz* faz*(-zc*vnuc(coamz,cobpz,1)+cerf*verf(coamz,cobpz,1))+&
                                           faz*  fbz*(-zc*vnuc(coamz,cobmz,1)+cerf*verf(coamz,cobmz,1))
                      ELSE
                        pVp(pVpa,pVpb)=ftaz*ftbz*(-zc*vnuc(coapx,cobpx,1))-&
                                           ftaz* fbx*(-zc*vnuc(coapx,cobmx,1))-&
                                           ftbz* fax*(-zc*vnuc(coamx,cobpx,1))+&
                                           fax*  fbx*(-zc*vnuc(coamx,cobmx,1))+&
                                          ftaz*ftbz*(-zc*vnuc(coapy,cobpy,1))-&
                                           ftaz* fby*(-zc*vnuc(coapy,cobmy,1))-&
                                           ftbz* fay*(-zc*vnuc(coamy,cobpy,1))+&
                                           fay*  fby*(-zc*vnuc(coamy,cobmy,1))+&
                                          ftaz*ftbz*(-zc*vnuc(coapz,cobpz,1))-&
                                           ftaz* fbz*(-zc*vnuc(coapz,cobmz,1))-&
                                           ftbz* faz*(-zc*vnuc(coamz,cobpz,1))+&
                                           faz*  fbz*(-zc*vnuc(coamz,cobmz,1))
                      END IF
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO

          DO j=1,ncoset(lb_max1)
            DO i=ncoset(la_min1-1)+1,ncoset(la_max1-maxder_local)
              pVp_sum(na+i,nb+j) = pVp_sum(na+i,nb+j) +pVp(i,j)
            END DO
          END DO
        END IF
!JTe
        nb = nb + ncoset(lb_max1)                                      !JT

      END DO

      na = na + ncoset(la_max1 - maxder_local)                         !JT
      nap = nap + ncoset(la_max1)                                      !JT

    END DO

  END SUBROUTINE verfc

END MODULE ai_verfc

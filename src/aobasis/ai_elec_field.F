!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
!>       - coset       : Cartesian orbital set pointer.
!>       - dab         : Distance between the atomic centers a and b.
!>       - dac         : Distance between the atomic centers a and c.
!>       - dbc         : Distance between the atomic centers b and c.
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
!> \author VW
! *****************************************************************************
MODULE ai_elec_field
  USE ai_os_rr,                        ONLY: os_rr_coul
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_elec_field'
  PRIVATE

 ! *** Public subroutines ***

  PUBLIC :: efg

CONTAINS

! *****************************************************************************
!> \brief   Calculation of the primitive electric field integrals over
!>          Cartesian Gaussian-type functions.
!> \param la_max ...
!> \param la_min ...
!> \param npgfa ...
!> \param rpgfa ...
!> \param zeta ...
!> \param lb_max ...
!> \param lb_min ...
!> \param npgfb ...
!> \param rpgfb ...
!> \param zetb ...
!> \param rac ...
!> \param rbc ...
!> \param rab ...
!> \param vab ...
!> \param ldrr1 ...
!> \param ldrr2 ...
!> \param rr ...
!> \date    02.03.2009
!> \author  VW
!> \version 1.0
! *****************************************************************************
  SUBROUTINE efg(la_max,la_min,npgfa,rpgfa,zeta,&
                 lb_max,lb_min,npgfb,rpgfb,zetb,&
                 rac,rbc,rab,vab,ldrr1,ldrr2,rr)
    INTEGER, INTENT(IN)                      :: la_max, la_min, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa, zeta
    INTEGER, INTENT(IN)                      :: lb_max, lb_min, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb, zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac, rbc, rab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: vab
    INTEGER, INTENT(IN)                      :: ldrr1, ldrr2
    REAL(KIND=dp), &
      DIMENSION(0:ldrr1-1, ldrr2, *), &
      INTENT(INOUT)                          :: rr

    CHARACTER(len=*), PARAMETER :: routineN = 'efg', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ay, az, bx, by, bz, coa, coam1x, coam1y, coam1z, coam2x, &
      coam2y, coam2z, coamxpy, coamxpz, coamxy, coamxz, coamypx, coamypz, &
      coamyz, coamzpx, coamzpy, coap1x, coap1y, coap1z, coap2x, coap2y, &
      coap2z, coapxy, coapxz, coapyz, cob, cobm1x, cobm1y, cobm1z, cobm2x, &
      cobm2y, cobm2z, cobmxpy, cobmxpz, cobmxy, cobmxz, cobmypx, cobmypz, &
      cobmyz, cobmzpx, cobmzpy, cobp1x, cobp1y, cobp1z, cobp2x, cobp2y, &
      cobp2z, cobpxy, cobpxz, cobpyz, i, ipgf, j, jpgf, la, lb, ma, mb, na, nb
    REAL(KIND=dp)                            :: dab, dum, dumxx, dumxy, &
                                                dumxz, dumyy, dumyz, dumzz, &
                                                f0, rab2, xhi, za, zb, zet
    REAL(KIND=dp), DIMENSION(3)              :: rap, rbp, rcp

! *** Calculate the distance of the centers a and c ***

    rab2 = rab(1)**2+rab(2)**2+rab(3)**2
    dab = SQRT(rab2)

    ! *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0

    DO ipgf=1,npgfa

       nb = 0

       DO jpgf=1,npgfb

          ! *** Screening ***

          IF(rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
             DO j=nb+1,nb+ncoset(lb_max)
                DO i=na+1,na+ncoset(la_max)
                   vab(i,j,1) = 0.0_dp
                   vab(i,j,2) = 0.0_dp
                   vab(i,j,3) = 0.0_dp
                   vab(i,j,4) = 0.0_dp
                   vab(i,j,5) = 0.0_dp
                   vab(i,j,6) = 0.0_dp
                ENDDO
             ENDDO
             nb = nb + ncoset(lb_max)
             CYCLE
          ENDIF

          ! *** Calculate some prefactors ***

          za = zeta(ipgf)
          zb = zetb(jpgf)
          zet  =  za+zb
          xhi  =  za*zb/zet
          rap  =  zb*rab/zet
          rbp  = -za*rab/zet
          rcp  = -(za*rac+zb*rbc)/zet
          f0 = 2.0_dp*SQRT(zet/pi)*(pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

          ! *** Calculate the reccurence relation

          CALL os_rr_coul(rap,la_max+2,rbp,lb_max+2,rcp,zet,ldrr1,ldrr2,rr)

          ! *** Calculate the primitive electric field gradient integrals ***

          DO lb = lb_min,lb_max
          DO bx = 0,lb
          DO by = 0,lb-bx
             bz = lb - bx - by
             cob    = coset(bx,by,bz)
             cobm1x = coset(MAX(bx-1,0),by,bz)
             cobm1y = coset(bx,MAX(by-1,0),bz)
             cobm1z = coset(bx,by,MAX(bz-1,0))
             cobm2x = coset(MAX(bx-2,0),by,bz)
             cobm2y = coset(bx,MAX(by-2,0),bz)
             cobm2z = coset(bx,by,MAX(bz-2,0))
             cobmxy = coset(MAX(bx-1,0),MAX(by-1,0),bz)
             cobmxz = coset(MAX(bx-1,0),by,MAX(bz-1,0))
             cobmyz = coset(bx,MAX(by-1,0),MAX(bz-1,0))
             cobp1x = coset(bx+1,by,bz)
             cobp1y = coset(bx,by+1,bz)
             cobp1z = coset(bx,by,bz+1)
             cobp2x = coset(bx+2,by,bz)
             cobp2y = coset(bx,by+2,bz)
             cobp2z = coset(bx,by,bz+2)
             cobpxy = coset(bx+1,by+1,bz)
             cobpxz = coset(bx+1,by,bz+1)
             cobpyz = coset(bx,by+1,bz+1)
             cobmxpy= coset(MAX(bx-1,0),by+1,bz)
             cobmxpz= coset(MAX(bx-1,0),by,bz+1)
             cobmypx= coset(bx+1,MAX(by-1,0),bz)
             cobmypz= coset(bx,MAX(by-1,0),bz+1)
             cobmzpx= coset(bx+1,by,MAX(bz-1,0))
             cobmzpy= coset(bx,by+1,MAX(bz-1,0))
             mb = nb + cob
             DO la = la_min,la_max
             DO ax = 0,la
             DO ay = 0,la-ax
                az = la - ax - ay
                coa    = coset(ax,ay,az)
                coap1x = coset(ax+1,ay,az)
                coap1y = coset(ax,ay+1,az)
                coap1z = coset(ax,ay,az+1)
                coap2x = coset(ax+2,ay,az)
                coap2y = coset(ax,ay+2,az)
                coap2z = coset(ax,ay,az+2)
                coapxy = coset(ax+1,ay+1,az)
                coapxz = coset(ax+1,ay,az+1)
                coapyz = coset(ax,ay+1,az+1)
                coam1x = coset(MAX(ax-1,0),ay,az)
                coam1y = coset(ax,MAX(ay-1,0),az)
                coam1z = coset(ax,ay,MAX(az-1,0))
                coam2x = coset(MAX(ax-2,0),ay,az)
                coam2y = coset(ax,MAX(ay-2,0),az)
                coam2z = coset(ax,ay,MAX(az-2,0))
                coamxy = coset(MAX(ax-1,0),MAX(ay-1,0),az)
                coamxz = coset(MAX(ax-1,0),ay,MAX(az-1,0))
                coamyz = coset(ax,MAX(ay-1,0),MAX(az-1,0))
                coamxpy= coset(MAX(ax-1,0),ay+1,az)
                coamxpz= coset(MAX(ax-1,0),ay,az+1)
                coamypx= coset(ax+1,MAX(ay-1,0),az)
                coamypz= coset(ax,MAX(ay-1,0),az+1)
                coamzpx= coset(ax+1,ay,MAX(az-1,0))
                coamzpy= coset(ax,ay+1,MAX(az-1,0))
                ma = na + coa
                !
                ! (a|xx|b)
                dum = 4.0_dp * ( za**2 * rr(0,coap2x,cob) + zb**2 * rr(0,coa,cobp2x) &
                     &         + 2.0_dp * za * zb * rr(0,coap1x,cobp1x) ) &
                     -2.0_dp * rr(0,coa,cob) * ( za * REAL(2*ax+1,dp) + zb * REAL(2*bx+1,dp) )
                IF(ax.GT.1) dum = dum + REAL(ax*(ax-1),dp) * rr(0,coam2x,cob)
                IF(bx.GT.1) dum = dum + REAL(bx*(bx-1),dp) * rr(0,coa,cobm2x)
                IF(ax.GT.0) dum = dum - 4.0_dp * zb * REAL(ax,dp) * rr(0,coam1x,cobp1x)
                IF(bx.GT.0) dum = dum - 4.0_dp * za * REAL(bx,dp) * rr(0,coap1x,cobm1x)
                IF(ax.GT.0.AND.bx.GT.0) dum = dum + 2.0_dp * REAL(ax*bx,dp) * rr(0,coam1x,cobm1x)
                dumxx = f0 * dum
                !
                ! (a|yy|b)
                dum = 4.0_dp * ( za**2 * rr(0,coap2y,cob) + zb**2 * rr(0,coa,cobp2y) &
                     &         + 2.0_dp * za * zb * rr(0,coap1y,cobp1y) ) &
                     -2.0_dp * rr(0,coa,cob) * ( za * REAL(2*ay+1,dp) + zb * REAL(2*by+1,dp) )
                IF(ay.GT.1) dum = dum + REAL(ay*(ay-1),dp) * rr(0,coam2y,cob)
                IF(by.GT.1) dum = dum + REAL(by*(by-1),dp) * rr(0,coa,cobm2y)
                IF(ay.GT.0) dum = dum - 4.0_dp * zb * REAL(ay,dp) * rr(0,coam1y,cobp1y)
                IF(by.GT.0) dum = dum - 4.0_dp * za * REAL(by,dp) * rr(0,coap1y,cobm1y)
                IF(ay.GT.0.AND.by.GT.0) dum = dum + 2.0_dp * REAL(ay*by,dp) * rr(0,coam1y,cobm1y)
                dumyy = f0 * dum
                !
                ! (a|zz|b)
                dum = 4.0_dp * ( za**2 * rr(0,coap2z,cob) + zb**2 * rr(0,coa,cobp2z) &
                     &         + 2.0_dp * za * zb * rr(0,coap1z,cobp1z) ) &
                     -2.0_dp * rr(0,coa,cob) * ( za * REAL(2*az+1,dp) + zb * REAL(2*bz+1,dp) )
                IF(az.GT.1) dum = dum + REAL(az*(az-1),dp) * rr(0,coam2z,cob)
                IF(bz.GT.1) dum = dum + REAL(bz*(bz-1),dp) * rr(0,coa,cobm2z)
                IF(az.GT.0) dum = dum - 4.0_dp * zb * REAL(az,dp) * rr(0,coam1z,cobp1z)
                IF(bz.GT.0) dum = dum - 4.0_dp * za * REAL(bz,dp) * rr(0,coap1z,cobm1z)
                IF(az.GT.0.AND.bz.GT.0) dum = dum + 2.0_dp * REAL(az*bz,dp) * rr(0,coam1z,cobm1z)
                dumzz = f0 * dum
                !
                ! (a|xy|b)
                dum = 4.0_dp * ( za**2 * rr(0,coapxy,cob) + zb**2 * rr(0,coa,cobpxy) &
                     &         + za * zb * ( rr(0,coap1x,cobp1y) + rr(0,coap1y,cobp1x) ) )
                IF(ax.GT.0) dum = dum - 2.0_dp * REAL(ax,dp) * &
                     &  ( za * rr(0,coamxpy,cob) + zb * rr(0,coam1x,cobp1y) )
                IF(ay.GT.0) dum = dum - 2.0_dp * REAL(ay,dp) * &
                     &  ( za * rr(0,coamypx,cob) + zb * rr(0,coam1y,cobp1x) )
                IF(ax.GT.0.AND.ay.GT.0) dum = dum + REAL(ax*ay,dp) * rr(0,coamxy,cob)
                IF(bx.GT.0) dum = dum - 2.0_dp * REAL(bx,dp) * &
                     &  ( zb * rr(0,coa,cobmxpy) + za * rr(0,coap1y,cobm1x) )
                IF(by.GT.0) dum = dum - 2.0_dp * REAL(by,dp) * &
                     &  ( zb * rr(0,coa,cobmypx) + za * rr(0,coap1x,cobm1y) )
                IF(bx.GT.0.AND.by.GT.0) dum = dum + REAL(bx*by,dp) * rr(0,coa,cobmxy)
                IF(ax.GT.0.AND.by.GT.0) dum = dum + REAL(ax*by,dp) * rr(0,coam1x,cobm1y)
                IF(ay.GT.0.AND.bx.GT.0) dum = dum + REAL(ay*bx,dp) * rr(0,coam1y,cobm1x)
                dumxy = f0 * dum
                !
                ! (a|xz|b)
                dum = 4.0_dp * ( za**2 * rr(0,coapxz,cob) + zb**2 * rr(0,coa,cobpxz) &
                     &         + za * zb * ( rr(0,coap1x,cobp1z) + rr(0,coap1z,cobp1x) ) )
                IF(ax.GT.0) dum = dum - 2.0_dp * REAL(ax,dp) * &
                     &  ( za * rr(0,coamxpz,cob) + zb * rr(0,coam1x,cobp1z) )
                IF(az.GT.0) dum = dum - 2.0_dp * REAL(az,dp) * &
                     &  ( za * rr(0,coamzpx,cob) + zb * rr(0,coam1z,cobp1x) )
                IF(ax.GT.0.AND.az.GT.0) dum = dum + REAL(ax*az,dp) * rr(0,coamxz,cob)
                IF(bx.GT.0) dum = dum - 2.0_dp * REAL(bx,dp) * &
                     &  ( zb * rr(0,coa,cobmxpz) + za * rr(0,coap1z,cobm1x) )
                IF(bz.GT.0) dum = dum - 2.0_dp * REAL(bz,dp) * &
                     &  ( zb * rr(0,coa,cobmzpx) + za * rr(0,coap1x,cobm1z) )
                IF(bx.GT.0.AND.bz.GT.0) dum = dum + REAL(bx*bz,dp) * rr(0,coa,cobmxz)
                IF(ax.GT.0.AND.bz.GT.0) dum = dum + REAL(ax*bz,dp) * rr(0,coam1x,cobm1z)
                IF(az.GT.0.AND.bx.GT.0) dum = dum + REAL(az*bx,dp) * rr(0,coam1z,cobm1x)
                dumxz = f0 * dum
                !
                ! (a|yz|b)
                dum = 4.0_dp * ( za**2 * rr(0,coapyz,cob) + zb**2 * rr(0,coa,cobpyz) &
                     &         + za * zb * ( rr(0,coap1y,cobp1z) + rr(0,coap1z,cobp1y) ) )
                IF(ay.GT.0) dum = dum - 2.0_dp * REAL(ay,dp) * &
                     &  ( za * rr(0,coamypz,cob) + zb * rr(0,coam1y,cobp1z) )
                IF(az.GT.0) dum = dum - 2.0_dp * REAL(az,dp) * &
                     &  ( za * rr(0,coamzpy,cob) + zb * rr(0,coam1z,cobp1y) )
                IF(ay.GT.0.AND.az.GT.0) dum = dum + REAL(ay*az,dp) * rr(0,coamyz,cob)
                IF(by.GT.0) dum = dum - 2.0_dp * REAL(by,dp) * &
                     &  ( zb * rr(0,coa,cobmypz) + za * rr(0,coap1z,cobm1y) )
                IF(bz.GT.0) dum = dum - 2.0_dp * REAL(bz,dp) * &
                     &  ( zb * rr(0,coa,cobmzpy) + za * rr(0,coap1y,cobm1z) )
                IF(by.GT.0.AND.bz.GT.0) dum = dum + REAL(by*bz,dp) * rr(0,coa,cobmyz)
                IF(ay.GT.0.AND.bz.GT.0) dum = dum + REAL(ay*bz,dp) * rr(0,coam1y,cobm1z)
                IF(az.GT.0.AND.by.GT.0) dum = dum + REAL(az*by,dp) * rr(0,coam1z,cobm1y)
                dumyz = f0 * dum
                !
                !
                vab(ma,mb,1) = ( 2.0_dp * dumxx - dumyy - dumzz ) / 3.0_dp !xx
                vab(ma,mb,2) = dumxy !xy
                vab(ma,mb,3) = dumxz !xz
                vab(ma,mb,4) = ( 2.0_dp * dumyy - dumzz - dumxx ) / 3.0_dp !yy
                vab(ma,mb,5) = dumyz !yz
                vab(ma,mb,6) = ( 2.0_dp * dumzz - dumxx - dumyy ) / 3.0_dp !zz
             ENDDO
             ENDDO
             ENDDO !la

          ENDDO
          ENDDO
          ENDDO !lb

          nb = nb + ncoset(lb_max)

       ENDDO

       na = na + ncoset(la_max)

    ENDDO

  END SUBROUTINE efg

END MODULE ai_elec_field

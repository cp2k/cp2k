!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xc_tpss_util [1.0] *
!!
!!   NAME
!!     xc_tpss_u
!!
!!   FUNCTION
!!     Functions needed by TPSS
!!
!!   AUTHOR
!!     JGH (04.01.2004)
!!
!!   MODIFICATION HISTORY
!!
!******************************************************************************

MODULE xc_tpss_util

! *****************************************************************************

  USE kinds,                           ONLY: dbl

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl),PARAMETER  :: small = 1.e-14_dbl
  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f83 = 8._dbl*f13, &
                          f1081 = 10._dbl/81._dbl


  PUBLIC :: tpssux, pz_u, pz_p, ccfun, cgefun, tpbec, tpbeca, tpbecs

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

      SUBROUTINE tpssux(rho,grho,tau,sx,v1x,v2x,vttx)
!       exchange part of TPSS functional 

    REAL(dbl), INTENT(IN)                    :: rho, grho, tau
    REAL(dbl), INTENT(OUT)                   :: sx, v1x, v2x, vttx

    REAL(dbl), PARAMETER :: b = 0.40e0_dbl, c = 1.59096e0_dbl, &
      e = 1.537e0_dbl, kappa = 0.804e0_dbl, mu = 0.21951e0_dbl, &
      se = 1.2397580409096_dbl

    REAL(dbl) :: alpha, dfxdg, dfxdr, dfxdt, dfxdx, dpdg, dpdr, dqbtda, &
      dqbtdp, dqbtdz, dxdp, dxdqbt, dxdz, dzdg, dzdr, dzdt, exunif, fx, p, &
      qbt, x, z

!sqrt(e)
! -----------------------------------------------------------------------------
!

        exunif = -3._dbl/(4._dbl*pi) * (3._dbl*pi*pi*rho)**f13
!
        IF ( ABS(tau) > small .AND. grho > small ) THEN
          p=grho/(4._dbl*(3._dbl*pi*pi)**f23*rho**f83)
          z=0.125_dbl*grho/rho/tau
          alpha = 5._dbl*f13*p * (1._dbl/z - 1._dbl)
          qbt = 0.45_dbl*(alpha-1._dbl)/SQRT(1._dbl+b*alpha*(alpha-1._dbl))+f23*p
          x = ( f1081 + c*z*z/(1._dbl+z*z)**2 ) * &
              p + 146._dbl/2025._dbl*qbt**2 &
              - 73._dbl/405._dbl * qbt * &
              SQRT ( 0.5_dbl*(0.6_dbl*z)**2+0.5_dbl*p*p) &
              + f1081**2*p*p/kappa+2._dbl*se*f1081*(0.6_dbl*z)**2+e*mu*p**3
          x = x/(1._dbl+se*p)**2
          fx = 1._dbl + kappa - kappa/(1._dbl + x/kappa)
        ELSE
          fx = 1._dbl
        ENDIF
        sx = rho*exunif*fX
!
        IF ( ABS(tau) > small .AND. grho > small ) THEN
          dfxdx  = 1._dbl/(1._dbl + x/kappa)**2
          dqbtda = 0.225_dbl*(2._dbl + b*alpha - b)/ &
                           (1._dbl + b*alpha*alpha - b*alpha)**1.5_dbl
          dqbtdp = dqbtda*5._dbl*f13*(1._dbl/z-1._dbl) + f23
          dqbtdz = dqbtda*5._dbl*f13*p*(-1._dbl/(z*z))
          dxdqbt = (146._dbl/2025._dbl * 2._dbl*qbt - &
                   73._dbl/405._dbl * SQRT ( 0.5_dbl*(0.6_dbl*z)**2+0.5_dbl*p*p ))/ &
                   (1._dbl+se*p)**2
          dxdz   = dxdqbt*dqbtdz + ( c*p*(2._dbl*z/(1._dbl+z*z)**2 - &
                   4._dbl*z**3/(1._dbl+z*z)**3 ) - 73._dbl/405._dbl * qbt*0.5_dbl * &
                   0.36_dbl*z / SQRT ( 0.5_dbl*(0.6_dbl*z)**2 + 0.5_dbl*p*p ) + &
                   4._dbl*se*f1081*0.36_dbl*z ) / (1._dbl+se*p)**2
          dxdp   = -2._dbl*se*x/(1._dbl+se*p) + dxdqbt*dqbtdp + &
                   ( ( f1081 + c*z*z/(1._dbl+z*z)**2 ) - 73._dbl/405._dbl*qbt * &
                   0.5_dbl*p/SQRT ( 0.5_dbl*(0.6_dbl*z)**2 + 0.5_dbl*p*p ) + &
                   2._dbl*f1081**2*p/kappa + 3._dbl*e*mu*p**2 )/(1._dbl+se*p)**2
          dpdr   = -f83*p/rho
          dpdg   = 2._dbl*p/SQRT(grho)
          dzdr   = -z/rho
          dzdg   = 2._dbl*z/SQRT(grho)
          dzdt   = -z/tau
          dfxdr  = dfxdx*(dxdp*dpdr + dxdz*dzdr)
          dfxdg  = dfxdx*(dxdp*dpdg + dxdz*dzdg)
          dfxdt  = dfxdx*dxdz*dzdt
!
          v1x  = exunif*(fx + f13*fx + rho*dfxdr)
          v2x  = rho*exunif*dfxdg/SQRT(grho)
          vttx = rho*exunif*dfxdt
        ELSE
          v1x  = exunif + f13*sx/rho
          v2x  = 0._dbl
          vttx = 0._dbl
        ENDIF

      END SUBROUTINE tpssux

! *****************************************************************************

      SUBROUTINE tpbec(rho,grho,sc,v1c,v2c)
!       PBE Correlation functional

    REAL(dbl), INTENT(IN)                    :: rho, grho
    REAL(dbl), INTENT(OUT)                   :: sc, v1c, v2c

    REAL(dbl), PARAMETER :: be = 0.06672455060314922_dbl, cx = -0.001667_dbl, &
      cxc0 = 0.002568_dbl, cc0 = -cx + cxc0, ga = 0.031090690869654895_dbl

    INTEGER                                  :: iflg
    REAL(dbl)                                :: a, aa, af, dadr, decdr, dhdr, &
                                                dhdt, dsda, dsdr, dsdt, dtdr, &
                                                ec, expe, h0, rs, s1, t, vc, &
                                                xkf, xks, xy, y

! -----------------------------------------------------------------------------

        rs    = (3._dbl/(4._dbl*pi*rho))**f13
        iflg  = 2
        IF(rs < 1.0_dbl) iflg=1
        CALL pz_u(rs,ec,vc,iflg)
        decdr = (vc-ec)/rho
        aa    = grho
        a     = SQRT(aa)
        xkf   = (9._dbl*pi/4._dbl)**f13/rs
        xks   = SQRT(4._dbl*xkf/pi)
        t     = a/(2._dbl*xks*rho)
        expe  = EXP(-ec/ga)
        af    = be/ga * (1._dbl/(expe-1._dbl))
        y     = af*t*t
        xy    = (1._dbl+y)/(1._dbl+y+y*y)
        s1    = 1._dbl+be/ga*t*t*xy
        h0    = ga * LOG(s1)
        dtdr  = -t*7._dbl/(6._dbl*rho)
        dadr  = af*af*expe/be*decdr
        dsda  = -be/ga * af * t**6 * (2._dbl+y) / (1._dbl+y+y*y)**2
        dsdt  = 2._dbl*be/ga * t * (1._dbl+2._dbl*y) / (1._dbl+y+y*y)**2
        dsdr  = dsda*dadr + dsdt*dtdr
        dhdt  = ga/s1*dsdt
        dhdr  = ga/s1*dsdr
        sc    = ec + h0
        v1c   = decdr + dhdr
        v2c   = dhdt*t/a

      END SUBROUTINE tpbec

! *****************************************************************************

      SUBROUTINE tpbeca(rhoa,grhoaa,sc,v1ca,v2ca)

    REAL(dbl), INTENT(IN)                    :: rhoa, grhoaa
    REAL(dbl), INTENT(OUT)                   :: sc, v1ca, v2ca

    REAL(dbl), PARAMETER :: be = 0.06672455060314922_dbl, cx = -0.001667_dbl, &
      cxc0 = 0.002568_dbl, cc0 = -cx + cxc0, ga = 0.031090690869654895_dbl

    INTEGER                                  :: iflg
    REAL(dbl) :: a, aa, af, dadra, decdra, dhdra, dhdt, dsda, dsdra, dsdt, &
      dtdra, ec, expe, grho, h0, phi, phi3, rho, rs, s1, t, vca, vcb, xkf, &
      xks, xy, y

! -----------------------------------------------------------------------------

        rho   = rhoa
        IF ( rho > small ) THEN
          phi   = 0.5_dbl*2._dbl**(2._dbl/3._dbl)
          phi3  = phi*phi*phi
          grho  = grhoaa
          rs    = (3._dbl/(4._dbl*pi*rho))**f13
          iflg=2
          IF(rs < 1.0_dbl) iflg=1
          CALL pz_p(rs,1._dbl,ec,vca,vcb,iflg)
!
          aa    = DMAX1(grho,small)
          a     = SQRT(aa)
          xkf   = (9._dbl*pi/4._dbl)**f13/rS
          xks   = SQRT(4._dbl*xkf/pi)
          t     = a/(2._dbl*xks*rho*phi)
          expe  = EXP(-ec/(phi3*ga))
          af    = be/ga * (1._dbl/(expe-1._dbl))
          y     = af*t*t
          xy    = (1._dbl+y)/(1._dbl+y+y*y)
          s1    = 1._dbl+be/ga*t*t*xy
          h0    = ga*phi3 * LOG(s1)
!
          decdra= (vca-ec)/rho
          dtdra = -t*(7._dbl/(6._dbl*rho))
          dadra = af*af*expe/(-be*phi3)*(-decdra)
          dsda  = -be/ga * af * t**6 * (2._dbl+y) / (1._dbl+y+y*y)**2
          dsdt  = 2._dbl*be/ga * t * (1._dbl+2._dbl*y) / (1._dbl+y+y*y)**2
          dsdra = dsda*dadra + dsdt*dtdra
          dhdt  = ga*phi3/s1*dsdt
          dhdra = ga*phi3/s1*dsdra
!
          sc    = ec + h0
          v1ca  = decdra + dhdra
          v2ca  = dhdt*t/a
        ELSE
          sc    = 0._dbl
          v1ca  = 0._dbl
          v2ca  = 0._dbl
        ENDIF

      END SUBROUTINE tpbeca

! *****************************************************************************

      SUBROUTINE tpbecs(rhoa,rhob,grhoaa,grhoab,grhobb,sc, &
                        v1ca,v2ca,v1cb,v2cb,v2cab)
!       PBE Correlation functional (spin polarized)

    REAL(dbl), INTENT(IN)                    :: rhoa, rhob, grhoaa, grhoab, &
                                                grhobb
    REAL(dbl), INTENT(OUT)                   :: sc, v1ca, v2ca, v1cb, v2cb, &
                                                v2cab

    REAL(dbl), PARAMETER                     :: be = 0.06672455060314922_dbl, &
                                                ga = 0.031090690869654895_dbl

    INTEGER                                  :: iflg
    REAL(dbl) :: a, aa, af, dadra, dadrb, decdra, decdrb, dedra, dedrb, &
      dhdra, dhdrb, dhdt, dphida, dphidb, dphide, dsda, dsdra, dsdrb, dsdt, &
      dtdphi, dtdra, dtdrb, ec, eta, expe, grho, h0, phi, phi3, rho, rs, s1, &
      t, vca, vcb, xkf, xks, xy, y

! -----------------------------------------------------------------------------

        rho=rhoa+rhob
        eta=(rhoa-rhob)/rho
        grho=grhoaa+2._dbl*grhoab+grhobb
        IF(ABS(eta) > 1._dbl) eta=SIGN(1.0_dbl,eta)
        rs = (0.75_dbl/(pi*rho))**f13
        iflg=2
        IF(rs < 1.0_dbl) IFLG=1
        CALL pz_p(rs,eta,ec,vca,vcb,iflg)
        decdra= (vca-ec)/rho
        decdrb= (vcb-ec)/rho
        phi   = 0.5_dbl*((1._dbl+eta)**f23+(1._dbl-eta)**f23)
        phi3  = phi*phi*phi
        aa    = DMAX1(grho,small)
        a     = SQRT(aa)
        xkf   = (9._dbl*pi/4._dbl)**f13/rs
        xks   = SQRT(4._dbl*xkf/pi)
        t     = a/(2._dbl*xks*rho*phi)
        expe  = EXP(-ec/(phi3*ga))
        af    = be/ga * (1._dbl/(expe-1._dbl))
        y     = af*t*t
        xy    = (1._dbl+y)/(1._dbl+y+y*y)
        s1    = 1._dbl+be/ga*t*t*xy
        h0    = ga*phi3 * LOG(s1)
        sc    = ec + h0
!
        IF(eta < 0.999999999999_dbl .AND. eta > -0.999999999999_dbl) THEN
          dtdphi= -t/phi
          dphide= 1._dbl/(3._dbl*(1._dbl+eta)**f13) - 1._dbl/(3._dbl*(1._dbl-eta)**f13)
          dedra = 2._dbl*rhob/(rho*rho)
          dedrb = -2._dbl*rhoa/(rho*rho)
          dphida= dphide*dedra
          dphidb= dphide*dedrb
          dtdra = -t*(dphida/phi+7._dbl/(6._dbl*rho))
          dtdrb = -t*(dphidb/phi+7._dbl/(6._dbl*rho))
          dadra = af*af*expe/(-be*phi3)*(3._dbl*ec/phi*dphida-decdra)
          dadrb = af*af*expe/(-be*phi3)*(3._dbl*ec/phi*dphidb-decdrb)
          dsda  = -be/ga * af * t**6 * (2._dbl+y) / (1._dbl+y+y*y)**2
          dsdt  = 2._dbl*be/ga * t * (1._dbl+2._dbl*y) / (1._dbl+y+y*y)**2
          dsdra = dsda*dadra + dsdt*dtdra
          dsdrb = dsda*dadrb + dsdt*dtdrb
          dhdt  = ga*phi3/s1*dsdt
          dhdra = 3._dbl*h0/phi*dphida + ga*phi3/s1*dsdra
          dhdrb = 3._dbl*h0/phi*dphidb + ga*phi3/s1*dsdrb
        ELSE
          dtdra = -t*(7._dbl/(6._dbl*rho))
          dtdrb = -t*(7._dbl/(6._dbl*rho))
          dadra = af*af*expe/(-be*phi3)*(-decdrA)
          dadrb = af*af*expe/(-be*phi3)*(-decdrb)
          dsda  = -be/ga * af * t**6 * (2._dbl+y) / (1._dbl+y+y*y)**2
          dsdt  = 2._dbl*be/ga * t * (1._dbl+2._dbl*y) / (1._dbl+y+y*y)**2
          dsdra = dsda*dadra + dsdt*dtdra
          dsdrb = dsda*dadrb + dsdt*dtdrb
          dhdt  = ga*phi3/s1*dsdt
          dhdra = ga*phi3/s1*dsdra
          dhdrb = ga*phi3/s1*dsdrb
        ENDIF
!
        v1ca  = decdra + dhdra
        v2ca  = dhdt*t/aa
        v1cb  = decdrb + dhdrb
        v2cb  = dhdt*t/aa
        v2cab = dhdt*t/aa

      END SUBROUTINE tpbecs

! *****************************************************************************

      SUBROUTINE pz_u(rs,epz,vpz,iflg)
!       J.P. Perdew and Alex Zunger PRB 23, 5048 (1981)

    REAL(dbl), INTENT(IN)                    :: rs
    REAL(dbl), INTENT(OUT)                   :: epz, vpz
    INTEGER, INTENT(IN)                      :: iflg

    REAL(dbl), PARAMETER :: a = 0.0311_dbl, b = -0.048_dbl, b1 = 1.0529_dbl, &
      b2 = 0.3334_dbl, c = 0.0020_dbl, d = -0.0116_dbl, gc = -0.1423_dbl

    REAL(dbl)                                :: dox, ox, rs1, rs2, xln

! -----------------------------------------------------------------------------

        epz=0.0_dbl
        vpz=0.0_dbl
        IF(iflg == 1) THEN
!         High density formula
          xln=LOG(rs)
          epz=a*xln+b+c*rs*xln+d*rs
          vpz=a*xln+(b-a/3._dbl)+2._dbl/3._dbl*c*rs*xln+(2._dbl*d-c)/3._dbl*rs
        ELSE IF(iflg == 2) THEN
!         Interpolation formula
          rs1=SQRT(rs)
          rs2=rs
          ox=1._dbl+b1*rs1+b2*rs2
          dox=1._dbl+7._dbl/6._dbl*b1*rs1+4._dbl/3._dbl*b2*rs2
          epz=gc/ox
          vpz=epz*dox/ox
        ENDIF

      END SUBROUTINE pz_u

! *****************************************************************************

      SUBROUTINE pz_p(rs,eta,epz,vpza,vpzb,iflg)
!       J.P. Perdew and Alex Zunger PRB 23, 5048 (1981)

    REAL(dbl), INTENT(IN)                    :: rs, eta
    REAL(dbl), INTENT(OUT)                   :: epz, vpza, vpzb
    INTEGER, INTENT(IN)                      :: iflg

    REAL(dbl), PARAMETER :: ap = 0.01555_dbl, au = 0.0311_dbl, &
      b1p = 1.3981_dbl, b1u = 1.0529_dbl, b2p = 0.2611_dbl, b2u = 0.3334_dbl, &
      bp = -0.0269_dbl, bu = -0.048_dbl, cp = 0.0007_dbl, cu = 0.0020_dbl, &
      dp = -0.0048_dbl, du = -0.0116_dbl, gcp = -0.0843_dbl, &
      gcu = -0.1423_dbl, small = 1.e-16_dbl

    REAL(dbl)                                :: dfeta, doxp, doxu, ep, eu, &
                                                f43, feta, oxp, oxu, rs1, &
                                                rs2, vp, vu, xln

! -----------------------------------------------------------------------------

      IF(iflg == 1) THEN
!       High density formula
        xln=LOG(rs)
        eu=au*xln+bu+cu*rs*xln+du*rs
        ep=ap*xln+bp+cp*rs*xln+dp*rs
        vu=au*xln+(bu-au/3._dbl)+2._dbl/3._dbl*cu*rs*xln+(2._dbl*du-cu)/3._dbl*rs
        vp=ap*xln+(bp-ap/3._dbl)+2._dbl/3._dbl*cp*rs*xln+(2._dbl*dp-cp)/3._dbl*rs
      ELSE IF(iflg == 2) THEN
!       Interpolation formula
        f43=4._dbl/3._dbl
        rs1=SQRT(rs)
        rs2=rs
        oxu=1._dbl+b1u*rs1+b2u*rs2
        doxu=1._dbl+7._dbl/6._dbl*b1u*rs1+f43*b2u*rs2
        eu=gcu/oxu
        vu=eu*doxu/oxu
        oxp=1._dbl+b1p*rs1+b2p*rs2
        doxp=1._dbl+7._dbl/6._dbl*b1p*rs1+f43*b2p*rs2
        ep=gcp/oxp
        vp=ep*doxp/oxp
      ENDIF
      IF (eta>(1._dbl-small)) THEN
         feta=0._dbl
         dfeta=0._dbl
      ELSE
         feta=((1._dbl+eta)**f43+(1._dbl-eta)**f43-2._dbl)/(2._dbl**f43-2._dbl)
         dfeta=f43*((1._dbl+eta)**f13-(1._dbl-eta)**f13)/(2._dbl**f43-2._dbl)
      END IF
      epz=eu+feta*(ep-eu)
      vpza=vu+feta*(vp-vu)+(ep-eu)*(1._dbl-eta)*dfeta
      vpzb=vu+feta*(vp-vu)+(ep-eu)*(-1._dbl-eta)*dfeta

      END SUBROUTINE pz_p

! *****************************************************************************

      SUBROUTINE cgefun(xsi,eta2,cxe,dcdx,dcde)

    REAL(dbl), INTENT(IN)                    :: xsi, eta2
    REAL(dbl), INTENT(OUT)                   :: cxe, dcdx, dcde

    REAL(dbl)                                :: cx0, dno, dx0, dxx, f43, f73, &
                                                x2, xx

        IF(ABS(xsi) < small .AND. ABS(eta2) < small*small) THEN
          cxe  = 0.53_dbl
          dcdx = 0._dbl
          dcde = 0._dbl
        ELSE IF(ABS(ABS(xsi)-1._dbl) < small .AND. &
                ABS(eta2) < small*small) THEN
          x2  = xsi*xsi
          cx0 = 0.53_dbl + 0.87_dbl*x2 + 0.50_dbl * x2*x2 + &
                           2.26_dbl * x2*x2*x2
          dx0 = 2._dbl*xsi * ( 0.87_dbl + x2 + 3._dbl**2.26_dbl * x2*x2 )
          cxe  = cx0
          dcdx = dx0
          dcde = 0._dbl
        ELSE
          f43 = -4._dbl/3._dbl
          f73 = -7._dbl/3._dbl
          x2  = xsi*xsi
          cx0 = 0.53_dbl + 0.87_dbl*x2 + 0.50_dbl * x2*x2 + 2.26_dbl * x2*x2*x2
          dx0 = 2._dbl*xsi * ( 0.87_dbl + x2 + 3._dbl**2.26_dbl * x2*x2 )
          IF(ABS(ABS(xsi)-1._dbl) < small) THEN
            xx  = 2._dbl**f43
            dxx = f43*2._dbl**f73
            dxx = dxx*SIGN(1._dbl,xsi)
          ELSE
            xx  = (1._dbl + xsi)**f43 + (1._dbl - xsi)**f43
            dxx = f43*(1._dbl + xsi)**f73 - f43*(1._dbl - xsi)**f73
          ENDIF
          dno = (1._dbl + 0.5_dbl*eta2*xx)
          cxe = cx0/dno**4
          dcdx = dx0/dno**4 - 4._dbl*cxe*0.5_dbl*eta2*dxx/dno
          dcde = -2._dbl*cxe*xx/dno
        ENDIF

      END SUBROUTINE cgefun

! *****************************************************************************

      SUBROUTINE ccfun(rhoa,rhob,grhoaa,grhobb,grhoab,cxe,dcdra,dcdrb,&
                       dcdgaa,dcdgbb,dcdgab)

    REAL(dbl), INTENT(IN)                    :: rhoa, rhob, grhoaa, grhobb, &
                                                grhoab
    REAL(dbl), INTENT(OUT)                   :: cxe, dcdra, dcdrb, dcdgaa, &
                                                dcdgbb, dcdgab

    REAL(dbl), PARAMETER                     :: f43 = 4._dbl/3._dbl

    REAL(dbl)                                :: cnn, dcde, dcdx, dedgaa, &
                                                dedgab, dedgbb, dedra, dedrb, &
                                                dedx, dxdra, dxdrb, dxsi0, &
                                                eta2, oo, rho, xsi

        rho    = rhoa + rhob
        xsi    = ( rhoa - rhob ) / rho
        IF(ABS(xsi) > 1._dbl) xsi=SIGN(1._dbl,xsi)

!-----------------------------------------------------------------------
! for some reason this (correct version) is not working
!
!      dxsi0  = (xsi-1._dbl)**2 * grhoaa + 2._dbl*(xsi*xsi-1._dbl)*grhoab &
!                + (xsi+1._dbl)**2 * grhobb
!
! this simplified version works (but is not correct)
!
       dxsi0  = (xsi-1._dbl)**2 * grhoaa + (xsi+1._dbl)**2 * grhobb
!-----------------------------------------------------------------------

        oo     = 0.25_dbl/(rho*rho)/(3._dbl*pi*pi*rho)**f23
        eta2   = dxsi0*oo
        CALL cgefun(xsi,eta2,cnn,dcdx,dcde)
        dxdra  =  2._dbl*rhob/(rho*rho)
        dxdrb  = -2._dbl*rhoa/(rho*rho)
        IF(dxsi0 < small*small) THEN
          dedx   = 0._dbl
          dedgaa = 0._dbl
          dedgbb = 0._dbl
          dedgab = 0._dbl
        ELSE

!-----------------------------------------------------------------------
!         dedx   = (2._dbl*(xsi-1._dbl) * grhoaa + 4._dbl*xsi*grhoab &
!                   + 2._dbl*(xsi+1._dbl) * grhobb)*oO
!         dedgaa = (xsi-1._dbl)**2 * oo
!         dedgbb = (xsi+1._dbl)**2 * oo
!         dedgab = 2._dbl*oo*(xsi*xsi-1._dbl)

          dedx   = (2._dbl*(xsi-1._dbl)*grhoaa+2._dbl*(xsi+1._dbl)*grhobb) * oo
          dedgaa = (xsi-1._dbl)**2 * oo
          dedgbb = (xsi+1._dbl)**2 * oo
          dedgab = 0._dbl
!-----------------------------------------------------------------------

        ENDIF

        dedra  = dedx*dxdra - 2._dbl*f43*eta2/rhO
        dedrb  = dedx*dxdrb - 2._dbl*f43*eta2/rho
!
        dcdra  = dcdx*dxdra + dcde*dedra
        dcdrb  = dcdx*dxdrb + dcde*dedrb
        dcdgaa = dcde*dedgaa
        dcdgbb = dcde*dedgbb
        dcdgab = dcde*dedgab
        cxe    = cnn

      END SUBROUTINE CCFUN

! *****************************************************************************

END MODULE xc_tpss_util

! *****************************************************************************

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xc_tpss_p [1.0] *
!!
!!   NAME
!!     xc_tpss_p
!!
!!   FUNCTION
!!     Calculate the TPSS Functional (spin polarized)
!!
!!   AUTHOR
!!     JGH (04.01.2004)
!!
!!   MODIFICATION HISTORY
!!
!******************************************************************************

MODULE xc_tpss_p

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE xc_tpss_util, ONLY: tpssux, pz_p, tpbec, tpbeca, tpbecs, ccfun

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl),PARAMETER  :: small = 1.e-14_dbl
  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f83 = 8._dbl*f13, &
                          f1081 = 10._dbl/81._dbl

  PUBLIC :: tpss_p

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

      SUBROUTINE tpss_p(rhoa,rhob,grhoaa,grhobb,grhoab,taua,taub, &
                        sx,v1xa,v2xa,v1xb,v2xb,vttxa,vttxb,sc,v1ca,v2ca, &
                        v1cb,v2cb,v2xab,v2cab,vttca,vttcb)

        REAL(dbl), INTENT(IN) :: rhoa,rhob,grhoaa,grhobb,grhoab,taua,taub
        REAL(dbl), INTENT(OUT):: sx,v1xa,v2xa,v1xb,v2xb,vttxa,vttxb,sc,v1ca,v2ca, &
                                 v1cb,v2cb,v2xab,v2cab,vttca,vttcb

        REAL(dbl) :: rho,grho,tau,sxx,v1x,v2x,vtt

        IF(rhoa > small .OR. rhob > small) THEN
!         Exchange part
!         alpha spin
          rho=2._dbl*rhoa
          grho=4._dbl*grhoaa
          tau=2._dbl*taua
          IF(rho > small) THEN
            CALL tpssux(rho,grho,tau,sxx,v1x,v2x,vtt)
            sx=0.5_dbl*sxx
            v1xa=v1x
            v2xa=2._dbl*v2x
            vttxa=vtt
          ELSE
            sx=0._dbl
            v1xa=0._dbl
            v2xa=0._dbl
            vttxa=0._dbl
          ENDIF
!         beta spin
          rho=2._dbl*rhob
          grho=4._dbl*grhobb
          tau=2._dbl*taub
          IF(RHO > SMALL) THEN
            CALL tpssux(rho,grho,tau,sxx,v1x,v2x,vtt)
            sx=sx + 0.5_dbl*sxx
            v1xb=v1x
            v2xb=2._dbl*v2x
            vttxb=vtt
          ELSE
            v1xb=0._dbl
            v2xb=0._dbl
            vttxb=0._dbl
          ENDIF
          v2xab=0._dbl
!         Correlation
          CALL tpsspc(rhoa,rhob,grhoaa,grhobb,grhoab,taua,taub,&
                      sc,v1ca,v2ca,v1cb,v2cb,v2cab,vttca,vttcb)
        ELSE
          sx=0._dbl
          sc=0._dbl
          v1xa=0._dbl
          v1xb=0._dbl
          v2xa=0._dbl
          v2xb=0._dbl
          v1ca=0._dbl
          v1cb=0._dbl
          v2ca=0._dbl
          v2cb=0._dbl
          v2xab=0._dbl
          v2cab=0._dbl
          vttxa=0._dbl
          vttxb=0._dbl
          vttca=0._dbl
          vttcb=0._dbl
        ENDIF

      END SUBROUTINE tpss_p

! *****************************************************************************

      SUBROUTINE tpsspc(rhoa,rhob,grhoaa,grhobb,grhoab,taua,taub, &
                        sc,v1ca,v2ca,v1cb,v2cb,v2cab,vtta,vttb)

        REAL(dbl), INTENT(IN) :: rhoa,rhob,grhoaa,grhobb,grhoab,taua,taub
        REAL(dbl), INTENT(OUT):: sc,v1ca,v2ca,v1cb,v2cb,v2cab,vtta,vttb

        REAL(dbl), PARAMETER :: d = 2.80_dbl

        REAL(dbl) :: rho,grho,tau,z,dzdr,rdzdt,dzdg,e,&
                     dedra,dedrb,dedga,dedgb,dedgab,dedz
        REAL(dbl) :: op,edaz,edbz,ro,to,eta,rs,ec,vca,vcb
        INTEGER   :: iflg

        rho = rhoa + rhob
        grho = grhoaa + 2._dbl*grhoab + grhobb
        tau = taua + taub
        IF( ABS(tau) > small .AND. grho > small) THEN
          z    = 0.125_dbl*grho/rho/tau
          dzdr = -z/rho
          rdzdt = -z*rho/tau
          dzdg = z/grho
          CALL revpkzbp(rhoa,rhob,grhoaa,grhobb,grhoab,z,e, &
                         dedra,dedrb,dedga,dedgb,dedgab,dedz)
          op   = 1._dbl + d * e * z**3
          sc   = rho * e * op
          edaz = dedra + dedz*dzdr
          v1ca = e*op + rho*edaz*op + rho*e*d*z*z*(edaz*z+3._dbl*e*dzdr)
          edbz = dedrb + dedz*dzdr
          v1cb = e*op + rho*edbz*op + rho*e*d*z*z*(edbz*z+3._dbl*e*dzdr)
          ro   = e*d*e*3._dbl*z*z
          to   = rho * ( 1._dbl + 2._dbl*d*e*z**3 )
          v2ca = 2._dbl*rho*ro*dzdg + 2._dbl*to*dedz*dzdg + to*dedga
          v2cb = 2._dbl*rho*ro*dzdg + 2._dbl*to*dedz*dzdg + to*dedgb
          v2cab= 2._dbl*rho*ro*dzdg + 2._dbl*to*dedz*dzdg + to*dedgab
          vtta = vtta + dedz*rdzdt*op + e*d*dedz*rdzdt*z**3 + ro*rdzdt
          vttb = vttb + dedz*rdzdt*op + e*d*dedz*rdzdt*z**3 + ro*rdzdt
      ELSE
        eta=(rhoa-rhob)/rho
        IF(ABS(eta) > 1._dbl) eta=SIGN(1._dbl,eta)
        rs = (0.75_dbl/(pi*rho))**f13
        iflg=2
        IF(rs < 1.0_dbl) iflg=1
        CALL pz_p(rs,eta,ec,vca,vcb,iflg)
        sc   = rho * ec
        v1ca = vca
        v2ca = 0._dbl
        v1cb = vcb
        v2cb = 0._dbl
        v2cab= 0._dbl
      ENDIF

      END SUBROUTINE tpsspc

! *****************************************************************************

      SUBROUTINE revpkzbp(rhoa,rhob,grhoaa,grhobb,grhoab,z,e, &
                         dedra,dedrb,dedga,dedgb,dedgab,dedz)

        REAL(dbl), INTENT(IN) :: rhoa,rhob,grhoaa,grhobb,grhoab,z
        REAL(dbl), INTENT(OUT):: e,dedra,dedrb,dedga,dedgb,dedgab,dedz

        REAL(dbl),PARAMETER :: C00=0.54_dbl
  
        REAL(dbl) :: sc,v1ca,v2ca,v1cb,v2cb,v2cab,sca,tv1ca,tv2ca,scb,tv1cb,tv2cb
        REAL(dbl) :: eca,t1a,t2a,t2aba,txaba,tyaba
        REAL(dbl) :: ecb,t1b,t2b,t2abb,txabb,tyabb
        REAL(dbl) :: cxe,dcdra,dcdrb,dcdgaa,dcdgbb,dcdgab
        REAL(dbl) :: rho,ef1,ef2,e2

        CALL tpbecs(rhoa,rhob,grhoaa,grhoab,grhobb,sc,v1ca,v2ca,v1cb,v2cb,v2cab)
        CALL tpbeca(rhoa,grhoaa,sca,tv1ca,tv2ca)
        CALL tpbeca(rhob,grhobb,scb,tv1cb,tv2cb)

        IF(sca > sc) THEN
          eca = sca
          t1a = tv1ca
          IF(grhoaa > small) THEN
             t2a = tv2ca/SQRT(grhoaa)
          ELSE
             t2a = 0._dbl
          ENDIF
          t2aba = 0._dbl
          txaba = 0._dbl
          tyaba = 0._dbl
        ELSE
          eca = sc
          t1a = v1ca
          t2a = v2ca
          t2aba = v2cab
          txaba = v1cb
          tyaba = v2cB
        ENDIF
        IF(scb > sc) THEN
          ecb = scb
          t1b = tv1cb
          IF(grhobb > small) THEN
             t2b = tv2cb/SQRT(grhobb)
          ELSE
             t2b = 0._dbl
          ENDIF
          t2abb = 0._dbl
          txabb = 0._dbl
          tyabb = 0._dbl
        ELSE
          ecb = sc
          t1b = v1cb
          t2b = v2cb
          t2abb = v2cab
          txabb = v1ca
          tyabb = v2cA
        ENDIF
!
        CALL ccfun(rhoa,rhob,grhoaa,grhobb,grhoab,cxe,dcdra,dcdrb, &
                   dcdgaa,dcdgbb,dcdgab)
!
        rho    = rhoa + rhob
        ef1    = 1._dbl + cxe * z * z
        ef2    = ( 1._dbl + cxe ) * z * z
        e2     = rhoa/rho * eca + rhob/rho * ecb
        e      = sc * ef1 - ef2 * e2
        dedra  = v1ca * ef1 - (-e2+eca + rhoa*t1a + rhob*txabb)/rho * ef2 &
                 + dcdra*z*z * ( sc - e2 )
        dedrb  = v1cb * ef1 - (-e2+ecb + rhob*t1b + rhoa*txaba)/rho * ef2 &
                 + dcdrb*z*z * ( sc - e2 )
        dedga  = v2ca * ef1 - ef2 * ( rhoa/rho*t2a + rhob/rho*tyabb ) &
                 + 2._dbl*dcdgaa*z*z * ( sc - e2 )
        dedgb  = v2cb * ef1 - ef2 * ( rhob/rho*t2b + rhoa/rho*tyaba ) &
                 + 2._dbl*dcdgbb*z*z * ( sc - e2 )
        dedgab = v2cab * ef1 - ef2 * (rhoa/rho*t2aba + rhob/rho*t2abb) &
                 + 2._dbl*dcdgab*z*z * ( sc - e2 )
        dedz   = 2._dbl * sc * cxe * z - 2._dbl * ( 1._dbl + cxe ) * z * e2

      END SUBROUTINE revpkzbp

! *****************************************************************************

END MODULE xc_tpss_p

! *****************************************************************************

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/xc_tpss_u [1.0] *
!!
!!   NAME
!!     xc_tpss_u
!!
!!   FUNCTION
!!     Calculate the TPSS Functional (spin unpolarized)
!!
!!   AUTHOR
!!     JGH (04.01.2004)
!!
!!   MODIFICATION HISTORY
!!
!******************************************************************************

MODULE xc_tpss_u

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE xc_tpss_util, ONLY: tpssux, pz_u, tpbec, tpbeca

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(dbl),PARAMETER  :: small = 1.e-14_dbl
  REAL(dbl), PARAMETER :: pi = 3.14159265358979323846264338_dbl
  REAL(dbl), PARAMETER :: f13 = 1._dbl/3._dbl, &
                          f23 = 2._dbl*f13, &
                          f83 = 8._dbl*f13, &
                          f1081 = 10._dbl/81._dbl

  PUBLIC :: tpss_u

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

      SUBROUTINE tpss_u(rho,grho,tau,sx,v1x,v2x,vttx,sc,v1c,v2c,vttc)

        REAL(dbl), INTENT(IN) :: rho,grho,tau
        REAL(dbl), INTENT(OUT):: sx,v1x,v2x,vttx,sc,v1c,v2c,vttc

        IF ( rho > small ) THEN
!         Exchange part
          CALL tpssux(rho,grho,tau,sx,v1x,v2x,vttx)

!         Correlation
          CALL tpssuc(rho,grho,tau,sc,v1c,v2c,vttc)
        ELSE
          sx   = 0._dbl
          v1x  = 0._dbl
          v2x  = 0._dbl
          vttx = 0._dbl
          sc   = 0._dbl
          v1c  = 0._dbl
          v2c  = 0._dbl
          vttc = 0._dbl
        END IF

      END SUBROUTINE tpss_u

! *****************************************************************************

      SUBROUTINE tpssuc(rho,grho,tau,sc,v1c,v2c,vttc)
!       correlation part of TPSS functional 

        REAL(dbl), INTENT(IN) :: rho,grho,tau
        REAL(dbl), INTENT(OUT):: sc,v1c,v2c,vttc

        REAL(dbl),PARAMETER   :: d = 2.8e0_dbl

        REAL(dbl)  :: z,dzdr,dzdg,dzdt,e,dedr,dedg,dedz,dedt,edrt,edgt,op,ro,rs
        INTEGER    :: iflg

! -----------------------------------------------------------------------------

        IF ( ABS(tau) > small .AND. grho > small ) THEN
          z    = 0.125_dbl*grho/rho/tau
          dzdr = -z/rho
          dzdg = 2._dbl*z/sqrt(grho)
          dzdt = -z/tau
          CALL revpkzb(rho,grho,z,e,dedr,dedg,dedz)
          dedt = dedz*dzdt
          edrt = dedr + dedz*dzdr
          edgt = dedg + dedz*dzdg
          op  = 1._dbl + d * e * z**3
          ro  = rho * e * d * z**2
          sc  = rho * e * op
          v1c = (e + rho*edrt)*op + ro * (edrt*z+3._dbl*e*dzdr)
          v2c = rho*edgt*op + ro * (edgt*z+3._dbl*e*dzdg)
          v2c = v2c/sqrt(grho)
          vttc = rho*dedt*op + ro * (dedt*z+3._dbl*e*dzdt)
        ELSE
          rs   = (3._dbl/(4._dbl*pi*rho))**(1._dbl/3._dbl)
          iflg = 2
          if(rs.lt.1.0_dbl) iflg=1
          CALL pz_u(rs,e,v1c,iflg)
          sc   = rho * e
          v2c  = 0._dbl
          vttc = 0._dbl
        ENDIF

      END SUBROUTINE tpssuc

! *****************************************************************************

      SUBROUTINE revpkzb(rho,grho,z,e,dedr,dedg,dedz)
!       Revised PKZB tau correlation functional

        REAL(dbl), INTENT(IN) :: rho,grho,z
        REAL(dbl), INTENT(OUT):: e,dedr,dedg,dedz

        REAL(dbl),PARAMETER   :: c00 = 0.53_dbl
        REAL(dbl),PARAMETER   :: c0p = 1.53_dbl

        REAL(dbl)  :: ecpbe,v1c,v2c,ecpbea,v1ca,v2ca,tecpbe,tv1c,tv2c,z2,t1

! -----------------------------------------------------------------------------

        CALL tpbec(rho,grho,ecpbe,v1c,v2c)
        CALL tpbeca(0.5_dbl*rho,0.25_dbl*grho,ecpbea,v1ca,v2ca)
        IF ( ecpbea > ecpbe ) THEN
          tecpbe = ecpbea
          tv1c   = 0.5_dbl*v1ca
          tv2c   = 0.5_dbl*v2ca
        ELSE
          tecpbe = ecpbe
          tv1c   = v1c
          tv2c   = v2c
        END IF
        z2 = z*z
        t1 = 1._dbl + c00*z2
        e    = ecpbe*t1 - c0p*z2*tecpbe
        dedr = v1c*t1 - c0p*z2*tv1c
        dedg = v2c*t1 - c0p*z2*tv2c
        dedz = 2._dbl*z * (ecpbe*c00 - c0p*tecpbe)

      END SUBROUTINE revpkzb

! *****************************************************************************

END MODULE xc_tpss_u

! *****************************************************************************

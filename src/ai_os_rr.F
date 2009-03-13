!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE ai_os_rr
  USE f77_blas
  USE gamma,                           ONLY: fgamma
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: coset
#include "cp_common_uses.h"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_os_rr'
  PRIVATE

 ! *** Public subroutines ***

  PUBLIC :: os_rr_coul

CONTAINS

! *****************************************************************************
!> \brief   Calculation of the Obara-Saika reccurence relation for 1/r_C
!> \author  VW
!> \date    02.03.2009
!> \version 1.0
! *****************************************************************************
  SUBROUTINE os_rr_coul(rap,la_max,rbp,lb_max,rcp,zet,rr,ldrr1,ldrr2)
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rap
    INTEGER, INTENT(IN)                      :: la_max
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rbp
    INTEGER, INTENT(IN)                      :: lb_max
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rcp
    REAL(dp), INTENT(IN)                     :: zet
    INTEGER, INTENT(IN)                      :: ldrr1, ldrr2
    REAL(dp), DIMENSION(0:ldrr1-1, ldrr2, *), &
      INTENT(INOUT)                          :: rr

    INTEGER :: ax, ay, az, bx, by, bz, coa, coa1x, coa1y, coa1z, coa2x, &
      coa2y, coa2z, cob, cob1x, cob1y, cob1z, cob2x, cob2y, cob2z, la, lb, m, &
      mmax
    REAL(dp)                                 :: g, rcp2, t

!
! initialize
!     

    mmax = la_max+lb_max
    g = 0.5_dp/zet
    !
    ! rr(0:mmax) should be initialized before
    !
    rcp2 = rcp(1)**2+rcp(2)**2+rcp(3)**2
    t = zet*rcp2
    CALL fgamma(mmax,t,rr(0:mmax,1,1))
    !
    ! recursion in la with lb=0
    !
    DO la=1,la_max
       DO ax=0,la
       DO ay=0,la-ax
          az=la-ax-ay
          coa   = coset(ax,ay,az)
          coa1x = coset(MAX(ax-1,0),ay,az)
          coa1y = coset(ax,MAX(ay-1,0),az)
          coa1z = coset(ax,ay,MAX(az-1,0))
          coa2x = coset(MAX(ax-2,0),ay,az)
          coa2y = coset(ax,MAX(ay-2,0),az)
          coa2z = coset(ax,ay,MAX(az-2,0))
          IF(az.GT.0) THEN
             DO m=0,mmax-la
                rr(m,coa,1)=rap(3)*rr(m,coa1z,1)-rcp(3)*rr(m+1,coa1z,1)
             ENDDO
             IF(az.GT.1) THEN                                                        
                DO m=0,mmax-la
                   rr(m,coa,1)=rr(m,coa,1)+g*REAL(az-1,dp)*(rr(m,coa2z,1)-rr(m+1,coa2z,1))
                ENDDO
             ENDIF
          ELSEIF(ay.GT.0) THEN
             DO m=0,mmax-la
                rr(m,coa,1)=rap(2)*rr(m,coa1y,1)-rcp(2)*rr(m+1,coa1y,1)
             ENDDO
             IF(ay.GT.1) THEN
                DO m=0,mmax-la
                   rr(m,coa,1)=rr(m,coa,1)+g*REAL(ay-1,dp)*(rr(m,coa2y,1)-rr(m+1,coa2y,1))
                ENDDO
             ENDIF
          ELSEIF(ax.GT.0) THEN
             DO m=0,mmax-la
                rr(m,coa,1)=rap(1)*rr(m,coa1x,1)-rcp(1)*rr(m+1,coa1x,1)
             ENDDO
             IF(ax.GT.1) THEN
                DO m=0,mmax-la
                   rr(m,coa,1)=rr(m,coa,1)+g*REAL(ax-1,dp)*(rr(m,coa2x,1)-rr(m+1,coa2x,1))
                ENDDO
             ENDIF
          ELSE
             STOP
          ENDIF
       ENDDO
       ENDDO
    ENDDO
    !
    ! recursion in lb with all possible la
    !
    DO la=0,la_max
       DO ax=0,la
       DO ay=0,la-ax
          az=la-ax-ay
          coa   = coset(ax,ay,az)
          coa1x = coset(MAX(ax-1,0),ay,az)
          coa1y = coset(ax,MAX(ay-1,0),az)
          coa1z = coset(ax,ay,MAX(az-1,0))
          coa2x = coset(MAX(ax-2,0),ay,az)
          coa2y = coset(ax,MAX(ay-2,0),az)
          coa2z = coset(ax,ay,MAX(az-2,0))
          DO lb=1,lb_max
             DO bx=0,lb
             DO by=0,lb-bx
                bz=lb-bx-by
                cob   = coset(bx,by,bz)
                cob1x = coset(MAX(bx-1,0),by,bz)
                cob1y = coset(bx,MAX(by-1,0),bz)
                cob1z = coset(bx,by,MAX(bz-1,0))
                cob2x = coset(MAX(bx-2,0),by,bz)
                cob2y = coset(bx,MAX(by-2,0),bz)
                cob2z = coset(bx,by,MAX(bz-2,0))
                IF(bz.GT.0) THEN
                   DO m=0,mmax-la-lb
                      rr(m,coa,cob)=rbp(3)*rr(m,coa,cob1z)-rcp(3)*rr(m+1,coa,cob1z)
                   ENDDO
                   IF(bz.GT.1) THEN
                      DO m=0,mmax-la-lb
                         rr(m,coa,cob)=rr(m,coa,cob)+g*REAL(bz-1,dp)*(rr(m,coa,cob2z)-rr(m+1,coa,cob2z))
                      ENDDO
                   ENDIF
                   IF(az.GT.0) THEN
                      DO m=0,mmax-la-lb
                         rr(m,coa,cob)=rr(m,coa,cob)+g*REAL(az,dp)*(rr(m,coa1z,cob1z)-rr(m+1,coa1z,cob1z))
                      ENDDO
                   ENDIF
                ELSEIF(by.GT.0) THEN
                   DO m=0,mmax-la-lb
                      rr(m,coa,cob)=rbp(2)*rr(m,coa,cob1y)-rcp(2)*rr(m+1,coa,cob1y)
                   ENDDO
                   IF(by.GT.1) THEN
                      DO m=0,mmax-la-lb
                         rr(m,coa,cob)=rr(m,coa,cob)+g*REAL(by-1,dp)*(rr(m,coa,cob2y)-rr(m+1,coa,cob2y))
                      ENDDO
                   ENDIF
                   IF(ay.GT.0) THEN
                      DO m=0,mmax-la-lb
                         rr(m,coa,cob)=rr(m,coa,cob)+g*REAL(ay,dp)*(rr(m,coa1y,cob1y)-rr(m+1,coa1y,cob1y))
                      ENDDO
                   ENDIF
                ELSEIF(bx.GT.0) THEN
                   DO m=0,mmax-la-lb
                      rr(m,coa,cob)=rbp(1)*rr(m,coa,cob1x)-rcp(1)*rr(m+1,coa,cob1x)
                   ENDDO
                   IF(bx.GT.1) THEN
                      DO m=0,mmax-la-lb
                         rr(m,coa,cob)=rr(m,coa,cob)+g*REAL(bx-1,dp)*(rr(m,coa,cob2x)-rr(m+1,coa,cob2x))
                      ENDDO
                   ENDIF
                   IF(ax.GT.0) THEN
                      DO m=0,mmax-la-lb
                         rr(m,coa,cob)=rr(m,coa,cob)+g*REAL(ax,dp)*(rr(m,coa1x,cob1x)-rr(m+1,coa1x,cob1x))
                      ENDDO
                   ENDIF
                ELSE
                   STOP
                ENDIF
             ENDDO
             ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    !
  END SUBROUTINE os_rr_coul

END MODULE ai_os_rr

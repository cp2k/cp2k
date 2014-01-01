!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the Fermi contact integrals over Cartesian
!>        Gaussian-type functions.
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!> \par Parameters
!>      - ax,ay,az  : Angular momentum index numbers of orbital a.
!>      - bx,by,bz  : Angular momentum index numbers of orbital b.
!>      - coset     : Cartesian orbital set pointer.
!>      - dab       : Distance between the atomic centers a and b.
!>      - l{a,b}    : Angular momentum quantum number of shell a or b.
!>      - l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
!>      - l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
!>      - rab       : Distance vector between the atomic centers a and b.
!>      - rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
!>      - sab       : Shell set of overlap integrals.
!>      - zet{a,b}  : Exponents of the Gaussian-type functions a or b.
!>      - zetp      : Reciprocal of the sum of the exponents of orbital a and b.
!> \author Matthias Krack (08.10.1999)
! *****************************************************************************
MODULE ai_fermi_contact
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_fermi_contact'

! *** Public subroutines ***
  PUBLIC :: fermi_contact

CONTAINS

! *****************************************************************************
!> \brief   Purpose: Calculation of the two-center Fermi contact integrals
!>          4/3*pi*[a|delta(r-c)|b] over Cartesian Gaussian-type functions.
!> \author  VW
!> \date    27.02.2009
!> \version 1.0
! *****************************************************************************
  SUBROUTINE fermi_contact(la_max,la_min,npgfa,rpgfa,zeta,&
                           lb_max,lb_min,npgfb,rpgfb,zetb,&
                           rac,rbc,dab,fcab,ldfc)
    INTEGER, INTENT(IN)                      :: la_max, la_min, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa, zeta
    INTEGER, INTENT(IN)                      :: lb_max, lb_min, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb, zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac, rbc
    REAL(KIND=dp)                            :: dab
    INTEGER, INTENT(IN)                      :: ldfc
    REAL(KIND=dp), DIMENSION(ldfc, *), &
      INTENT(INOUT)                          :: fcab

    CHARACTER(len=*), PARAMETER :: routineN = 'fermi_contact', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ax, ay, az, bx, by, bz, coa, &
                                                cob, i, ipgf, j, jpgf, la, &
                                                lb, ma, mb, na, nb
    REAL(KIND=dp)                            :: dac2, dbc2, f0, fax, fay, &
                                                faz, fbx, fby, fbz

! *** Calculate some prefactors ***

    dac2 = rac(1)**2 + rac(2)**2 + rac(3)**2
    dbc2 = rbc(1)**2 + rbc(2)**2 + rbc(3)**2

    ! *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0

    DO ipgf=1,npgfa

       nb = 0

       DO jpgf=1,npgfb

          ! *** Screening ***

          IF (rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
             DO j=nb+1,nb+ncoset(lb_max)
                DO i=na+1,na+ncoset(la_max)
                   fcab(i,j) = 0.0_dp
                ENDDO
             ENDDO
             nb = nb + ncoset(lb_max)
             CYCLE
          ENDIF

          ! *** Calculate some prefactors ***

          f0 = 4.0_dp/3.0_dp*pi*EXP(-zeta(ipgf)*dac2-zetb(jpgf)*dbc2)

          ! *** Calculate the primitive Fermi contact integrals ***

          DO lb = lb_min,lb_max
          DO bx = 0,lb
             fbx = 1.0_dp
             IF(bx.GT.0) fbx = (rbc(1))**bx
          DO by = 0,lb-bx
             bz = lb - bx - by
             cob = coset(bx,by,bz)
             mb = nb + cob
             fby = 1.0_dp
             IF(by.GT.0) fby = (rbc(2))**by
             fbz = 1.0_dp
             IF(bz.GT.0) fbz = (rbc(3))**bz
             DO la = la_min,la_max
             DO ax = 0,la
                fax = fbx
                IF(ax.GT.0) fax = fbx * (rac(1))**ax
             DO ay = 0,la-ax
                az = la - ax - ay
                coa = coset(ax,ay,az)
                ma = na + coa
                fay = fby
                IF(ay.GT.0) fay = fby * (rac(2))**ay
                faz = fbz
                IF(az.GT.0) faz = fbz * (rac(3))**az

                fcab(ma,mb) = f0 * fax * fay * faz

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

  END SUBROUTINE fermi_contact

END MODULE ai_fermi_contact

!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of spin orbit integrals over Cartesian Gaussian-type functions
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
MODULE ai_spin_orbit
  USE ai_os_rr,                        ONLY: os_rr_coul
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_spin_orbit'
  PRIVATE

 ! *** Public subroutines ***

  PUBLIC :: pso

CONTAINS

! *****************************************************************************
!> \brief   Calculation of the primitive paramagnetic spin orbit integrals over
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
  SUBROUTINE pso(la_max,la_min,npgfa,rpgfa,zeta,lb_max,lb_min,npgfb,rpgfb,zetb,&
                rac,rbc,rab,vab,ldrr1,ldrr2,rr)
    INTEGER, INTENT(IN)                      :: la_max, la_min, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa, zeta
    INTEGER, INTENT(IN)                      :: lb_max, lb_min, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb, zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rac, rbc, rab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: vab
    INTEGER, INTENT(IN)                      :: ldrr1, ldrr2
    REAL(dp), &
      DIMENSION(0:ldrr1-1, ldrr2, *), &
      INTENT(INOUT)                          :: rr

    CHARACTER(len=*), PARAMETER :: routineN = 'pso', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ay, az, bx, by, bz, coa, coam1x, coam1y, coam1z, coap1x, &
      coap1y, coap1z, cob, cobm1x, cobm1y, cobm1z, cobp1x, cobp1y, cobp1z, i, &
      ipgf, j, jpgf, la, lb, ma, mb, na, nb
    REAL(dp)                                 :: dab, dum1, dum2, f0, rab2, &
                                                xhi, zet, zetab
    REAL(dp), DIMENSION(3)                   :: rap, rbp, rcp

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
                ENDDO
             ENDDO
             nb = nb + ncoset(lb_max)
             CYCLE
          ENDIF

          ! *** Calculate some prefactors ***

          zetab=  zeta(ipgf)*zetb(jpgf)
          zet  =  zeta(ipgf)+zetb(jpgf)
          xhi  =  zetab/zet
          rap  =  zetb(jpgf)*rab/zet
          rbp  = -zeta(ipgf)*rab/zet
          rcp  = -(zeta(ipgf)*rac+zetb(jpgf)*rbc)/zet

          f0 = 2.0_dp*SQRT(zet/pi)*(pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

          ! *** Calculate the reccurence relation ***

          CALL os_rr_coul(rap,la_max+1,rbp,lb_max+1,rcp,zet,ldrr1,ldrr2,rr)

          ! *** Calculate the primitive Fermi contact integrals ***

          DO lb = lb_min,lb_max
          DO bx = 0,lb
          DO by = 0,lb-bx
             bz = lb - bx - by
             cob    = coset(bx,by,bz)
             cobm1x = coset(MAX(bx-1,0),by,bz)
             cobm1y = coset(bx,MAX(by-1,0),bz)
             cobm1z = coset(bx,by,MAX(bz-1,0))
             cobp1x = coset(bx+1,by,bz)
             cobp1y = coset(bx,by+1,bz)
             cobp1z = coset(bx,by,bz+1)
             mb = nb + cob
             DO la = la_min,la_max
             DO ax = 0,la
             DO ay = 0,la-ax
                az = la - ax - ay
                coa    = coset(ax,ay,az)
                coam1x = coset(MAX(ax-1,0),ay,az)
                coam1y = coset(ax,MAX(ay-1,0),az)
                coam1z = coset(ax,ay,MAX(az-1,0))
                coap1x = coset(ax+1,ay,az)
                coap1y = coset(ax,ay+1,az)
                coap1z = coset(ax,ay,az+1)
                ma = na + coa
                !
                !
                ! (a|pso_x|b) = (4*zeta*zetb*(a+y||b+z)
                !               -2*zeta*Nz(b)*(a+y||b-z)-2*zetb*Ny(a)*(a-y||b+z)
                !               +Ny(a)*Nz(b)*(a-y||b-z))
                !              -(4*zeta*zetb*(a+z||b+y)
                !               -2*zeta*Ny(b)*(a+z||b-y)-2*zetb*Nz(a)*(a-z||b+y)
                !               +Nz(a)*Ny(b)*(a-z||b-y))
                dum1 = 4.0_dp * zeta(ipgf) * zetb(jpgf) * rr(0,coap1y,cobp1z)
                IF(bz.GT.0) dum1 = dum1 - 2.0_dp * zeta(ipgf) * REAL(bz,dp) * rr(0,coap1y,cobm1z)
                IF(ay.GT.0) dum1 = dum1 - 2.0_dp * zetb(jpgf) * REAL(ay,dp) * rr(0,coam1y,cobp1z)
                IF(ay.GT.0.AND.bz.GT.0) dum1 = dum1 + REAL(ay,dp) * REAL(bz,dp) * rr(0,coam1y,cobm1z)
                !
                dum2 = 4.0_dp * zeta(ipgf) * zetb(jpgf) * rr(0,coap1z,cobp1y)
                IF(by.GT.0) dum2 = dum2 - 2.0_dp * zeta(ipgf) * REAL(by,dp) * rr(0,coap1z,cobm1y)
                IF(az.GT.0) dum2 = dum2 - 2.0_dp * zetb(jpgf) * REAL(az,dp) * rr(0,coam1z,cobp1y)
                IF(az.GT.0.AND.by.GT.0) dum2 = dum2 + REAL(az,dp) * REAL(by,dp) * rr(0,coam1z,cobm1y)
                vab(ma,mb,1) = f0 * (dum1 - dum2)
                !
                !
                ! (a|pso_y|b) = (4*zeta*zetb*(a+z||b+x)
                !               -2*zeta*Nx(b)*(a+z||b-x)-2*zetb*Nz(a)*(a-z||b+x)
                !               +Nz(a)*Nx(b)*(a-z||b-x))
                !              -(4*zeta*zetb*(a+x||b+z)
                !               -2*zeta*Nz(b)*(a+x||b-z)-2*zetb*Nx(a)*(a-x||b+z)
                !               +Nx(a)*Nz(b)*(a-x||b-z))
                dum1 = 4.0_dp * zeta(ipgf) * zetb(jpgf) * rr(0,coap1z,cobp1x)
                IF(bx.GT.0) dum1 = dum1 - 2.0_dp * zeta(ipgf) * REAL(bx,dp) * rr(0,coap1z,cobm1x)
                IF(az.GT.0) dum1 = dum1 - 2.0_dp * zetb(jpgf) * REAL(az,dp) * rr(0,coam1z,cobp1x)
                IF(az.GT.0.AND.bx.GT.0) dum1 = dum1 + REAL(az,dp) * REAL(bx,dp) * rr(0,coam1z,cobm1x)
                !
                dum2 = 4.0_dp * zeta(ipgf) * zetb(jpgf) * rr(0,coap1x,cobp1z)
                IF(bz.GT.0) dum2 = dum2 - 2.0_dp * zeta(ipgf) * REAL(bz,dp) * rr(0,coap1x,cobm1z)
                IF(ax.GT.0) dum2 = dum2 - 2.0_dp * zetb(jpgf) * REAL(ax,dp) * rr(0,coam1x,cobp1z)
                IF(ax.GT.0.AND.bz.GT.0) dum2 = dum2 + REAL(ax,dp) * REAL(bz,dp) * rr(0,coam1x,cobm1z)
                vab(ma,mb,2) = f0 * (dum1 - dum2)
                !
                !
                ! (a|pso_z|b) = (4*zeta*zetb*(a+x||b+y)
                !               -2*zeta*Ny(b)*(a+x||b-y)-2*zetb*Nx(a)*(a-x||b+y)
                !               +Nx(a)*Ny(b)*(a-x||b-y))
                !              -(4*zeta*zetb*(a+y||b+x)
                !               -2*zeta*Nx(b)*(a+y||b-x)-2*zetb*Ny(a)*(a-y||b+x)
                !               +Ny(a)*Nx(b)*(a-y||b-x))
                dum1 = 4.0_dp * zeta(ipgf) * zetb(jpgf) * rr(0,coap1x,cobp1y)
                IF(by.GT.0) dum1 = dum1 - 2.0_dp * zeta(ipgf) * REAL(by,dp) * rr(0,coap1x,cobm1y)
                IF(ax.GT.0) dum1 = dum1 - 2.0_dp * zetb(jpgf) * REAL(ax,dp) * rr(0,coam1x,cobp1y)
                IF(ax.GT.0.AND.by.GT.0) dum1 = dum1 + REAL(ax,dp) * REAL(by,dp) * rr(0,coam1x,cobm1y)
                !
                dum2 = 4.0_dp * zeta(ipgf) * zetb(jpgf) * rr(0,coap1y,cobp1x)
                IF(bx.GT.0) dum2 = dum2 - 2.0_dp * zeta(ipgf) * REAL(bx,dp) * rr(0,coap1y,cobm1x)
                IF(ay.GT.0) dum2 = dum2 - 2.0_dp * zetb(jpgf) * REAL(ay,dp) * rr(0,coam1y,cobp1x)
                IF(ay.GT.0.AND.bx.GT.0) dum2 = dum2 + REAL(ay,dp) * REAL(bx,dp) * rr(0,coam1y,cobm1x)
                vab(ma,mb,3) = f0 * (dum1 - dum2)
                !
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

  END SUBROUTINE pso

END MODULE ai_spin_orbit

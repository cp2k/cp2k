!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/functionals [1.1] *
!!
!!   NAME
!!     functionals
!!
!!   FUNCTION
!!     Calculate the exchange and/or correlation energy and potential for
!!     closed shell or open shell densities.
!!
!!     ***  WARNING !!!  ***
!!
!!     *** DEPRECATED!!! ***
!!
!!     use xc_functionals
!!
!!   AUTHOR
!!     Matthias Krack (03.03.2000)
!!
!!   MODIFICATION HISTORY
!!     - HCTH functional added (01.08.2001, MK)
!!
!!   SOURCE
!******************************************************************************

MODULE functionals

! *****************************************************************************

! Index:

! SUBROUTINE init_functionals(rho_threshold,drho_threshold,xalpha)
! SUBROUTINE becke88_0(rho,drho,ex,vx,vxg)
! SUBROUTINE becke88_1(rho,drho,ex,vx,vxg,weight)
! SUBROUTINE becke88_2(rho,drho,ex,vx,vxg,weight)
! SUBROUTINE becke88_3(rho,drho,ex,vx,vxg,weight)
! SUBROUTINE becke88_sp(rhos,drhos,exs,vxs,vxgs)
! SUBROUTINE becke88_sp_0(rhoa,rhob,drhoa,drhob,ex,vxa,vxb,vxga,vxgb)
! SUBROUTINE becke88_sp_1(rhoa,rhob,drhoa,drhob,ex,vxa,vxb,vxga,vxgb,weight)
! SUBROUTINE becke88_sp_2(rhoa,rhob,drhoa,drhob,ex,vxa,vxb,vxga,vxgb,weight)
! SUBROUTINE becke88_sp_3(rhoa,rhob,drhoa,drhob,ex,vxa,vxb,vxga,vxgb,weight)
! SUBROUTINE hcth_0(iparset,rho,drho,exc,vxc,vxcg)
! SUBROUTINE hcth_1(iparset,rho,drho,exc,vxc,vxcg,weight)
! SUBROUTINE hcth_2(iparset,rho,drho,exc,vxc,vxcg,weight)
! SUBROUTINE hcth_3(iparset,rho,drho,exc,vxc,vxcg,weight)
! SUBROUTINE lyp_0(rho,drho,ec,vc,vcg)
! SUBROUTINE lyp_1(rho,drho,ec,vc,vcg,weight)
! SUBROUTINE lyp_2(rho,drho,ec,vc,vcg,weight)
! SUBROUTINE lyp_3(rho,drho,ec,vc,vcg,weight)
! SUBROUTINE pade_0(rho,exc,vxc)
! SUBROUTINE pade_1(rho,exc,vxc,weight)
! SUBROUTINE pade_2(rho,exc,vxc,weight)
! SUBROUTINE pade_3(rho,exc,vxc,weight)
! SUBROUTINE pbe_c_0(rho,drho,ec,vc,vcg)
! SUBROUTINE pbe_c_1(rho,drho,ec,vc,vcg,weight)
! SUBROUTINE pbe_c_2(rho,drho,ec,vc,vcg,weight)
! SUBROUTINE pbe_c_3(rho,drho,ec,vc,vcg,weight)
! SUBROUTINE pbe_x_0(rho,drho,ex,vx,vxg)
! SUBROUTINE pbe_x_1(rho,drho,ex,vx,vxg,weight)
! SUBROUTINE pbe_x_2(rho,drho,ex,vx,vxg,weight)
! SUBROUTINE pbe_x_3(rho,drho,ex,vx,vxg,weight)
! SUBROUTINE perdew86_c_0(rho,drho,ec,vc,vcg)
! SUBROUTINE perdew86_c_1(rho,drho,ec,vc,vcg,weight)
! SUBROUTINE perdew86_c_2(rho,drho,ec,vc,vcg,weight)
! SUBROUTINE perdew86_c_3(rho,drho,ec,vc,vcg,weight)
! SUBROUTINE perdew86_x_0(rho,drho,ex,vx,vxg)
! SUBROUTINE perdew86_x_1(rho,drho,ex,vx,vxg,weight)
! SUBROUTINE perdew86_x_2(rho,drho,ex,vx,vxg,weight)
! SUBROUTINE perdew86_x_3(rho,drho,ex,vx,vxg,weight)
! SUBROUTINE slater_0(rho,ex,vx)
! SUBROUTINE slater_1(rho,ex,vx,weight)
! SUBROUTINE slater_2(rho,ex,vx,weight)
! SUBROUTINE slater_3(rho,ex,vx,weight)
! SUBROUTINE vwn_c_0(rho,ec,vc)
! SUBROUTINE vwn_c_1(rho,ec,vc,weight)
! SUBROUTINE vwn_c_2(rho,ec,vc,weight)
! SUBROUTINE vwn_c_3(rho,ec,vc,weight)
! SUBROUTINE vwn_x_0(rho,ex,vx)
! SUBROUTINE vwn_x_1(rho,ex,vx,weight)
! SUBROUTINE vwn_x_2(rho,ex,vx,weight)
! SUBROUTINE vwn_x_3(rho,ex,vx,weight)

! *****************************************************************************

  USE kinds, ONLY: dp

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(KIND = dp), PARAMETER :: f13 = 1.0_dp/3.0_dp,&
                         f16 = 0.5_dp*f13,&
                         f23 = 2.0_dp*f13,&
                         f43 = 4.0_dp*f13,&
                         f76 = 7.0_dp*f16,&
                         f83 = 8.0_dp*f13

! *** Type to hold XC functional infos ***

  TYPE xc_type
    CHARACTER(LEN=60) :: c_functional,x_functional
    LOGICAL           :: gradient_functional,smooth_gradient
  END TYPE xc_type

! *** Global variables ***

  REAL(KIND = dp) :: beta_pbe,cf,c_perdew86,cx_slater_e,cx_slater_v,cx_vwn_e,&
              cx_vwn_v,eps_drho,eps_rho,gamma_pbe,pi,r2kf,r2ks,rsfac,two13

! *** Public subroutines ***

  PUBLIC :: becke88,hcth,init_functionals,lyp,pade,pbe_c,pbe_x,perdew86_x,&
            perdew86_c,slater,vwn_c,vwn_x

! *** Public data types ***

  PUBLIC :: xc_type

! *****************************************************************************

  INTERFACE becke88
    MODULE PROCEDURE becke88_0,becke88_1,becke88_2,becke88_3,&
                     becke88_sp_0,becke88_sp_1,becke88_sp_2,becke88_sp_3
  END INTERFACE

  INTERFACE hcth
    MODULE PROCEDURE hcth_0,hcth_1,hcth_2,hcth_3
  END INTERFACE

  INTERFACE lyp
    MODULE PROCEDURE lyp_0,lyp_1,lyp_2,lyp_3
  END INTERFACE

  INTERFACE pade
    MODULE PROCEDURE pade_0,pade_1,pade_2,pade_3
  END INTERFACE

  INTERFACE pbe_c
    MODULE PROCEDURE pbe_c_0,pbe_c_1,pbe_c_2,pbe_c_3
  END INTERFACE

  INTERFACE pbe_x
    MODULE PROCEDURE pbe_x_0,pbe_x_1,pbe_x_2,pbe_x_3
  END INTERFACE

  INTERFACE perdew86_c
    MODULE PROCEDURE perdew86_c_0,perdew86_c_1,perdew86_c_2,perdew86_c_3
  END INTERFACE

  INTERFACE perdew86_x
    MODULE PROCEDURE perdew86_x_0,perdew86_x_1,perdew86_x_2,perdew86_x_3
  END INTERFACE

  INTERFACE slater
    MODULE PROCEDURE slater_0,slater_1,slater_2,slater_3
  END INTERFACE

  INTERFACE vwn_c
    MODULE PROCEDURE vwn_c_0,vwn_c_1,vwn_c_2,vwn_c_3
  END INTERFACE

  INTERFACE vwn_x
    MODULE PROCEDURE vwn_x_0,vwn_x_1,vwn_x_2,vwn_x_3
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE init_functionals(rho_threshold,drho_threshold,xalpha)

!   Purpose: Initialize the constants used for the calculation of the
!            exchange-correlation functionals. The default value for
!            xalpha is 2/3.

!   History: - Creation (16.03.2000, Matthias Krack)
!            - Optional parameter rho_threshold added (02.08.2001, MK)

!   ***************************************************************************

    REAL(KIND = dp), OPTIONAL, INTENT(IN) :: drho_threshold,rho_threshold,xalpha

!   ---------------------------------------------------------------------------

    pi = 2.0_dp*ACOS(0.0_dp)
    two13 = 2.0_dp**f13

!   *** Density values smaller than the threshold value are set to zero ***

    IF (PRESENT(rho_threshold)) THEN
      eps_rho = rho_threshold
    ELSE
      eps_rho = 1.0E-12_dp
    END IF

    IF (PRESENT(drho_threshold)) THEN
      eps_drho = drho_threshold
    ELSE
      eps_drho = 1.0E-10_dp
    END IF

!   *** Coefficient of the Thomas-Fermi kinetic-energy density ***

    cf = 0.3_dp*(3.0_dp*pi**2)**f23

!   *** Prefactor for the calculation of the effective electron radius ***

    rsfac = (f43*pi)**(-f13)

!   *** alpha = 2/3 (VWN exchange with spin polarization) ***

    cx_vwn_e = -0.75_dp*(3.0_dp/pi)**f13
    cx_vwn_v = f43*cx_vwn_e

!   *** Prefactors of the X_alpha exchange functional (Slater) ***

    IF (PRESENT(xalpha)) THEN
      cx_slater_e = -1.125_dp*xalpha*(3.0_dp/pi)**f13
      cx_slater_v = f43*cx_slater_e
    ELSE
      cx_slater_e = cx_vwn_e
      cx_slater_v = cx_vwn_v
    END IF

!   *** Prefactor in the Perdew correlation functional (1986) ***

    c_perdew86 = (9.0_dp*pi)**f16

!   *** The reciprocal of two times the local Fermi wave vector ***

    r2kf = 0.5_dp*(3.0_dp*pi**2)**(-f13)

!   *** Prefactors in the PBE correlation functional (1996) ***

    r2ks = 0.25_dp*(3.0_dp/pi)**(-f16)

    beta_pbe = 0.004235_dp*(16.0_dp/pi)*(3.0_dp*pi**2)**f13
    gamma_pbe = (1.0_dp - LOG(2.0_dp))/pi**2

  END SUBROUTINE init_functionals

! *****************************************************************************

  SUBROUTINE becke88_0(rho,drho,ex,vx,vxg)

!   Purpose: Calculate the gradient-corrected exchange energy and potential
!            of Becke (1988) for a closed shell density.

!   Literature: A. D. Becke, Phys. Rev. A 38, 3098 (1988)

!   History: - Creation (03.03.2000, Matthias Krack)

!   ***************************************************************************

    REAL(KIND = dp), INTENT(IN)  :: drho,rho
    REAL(KIND = dp), INTENT(OUT) :: ex,vx,vxg

!   *** Local parameters ***

    REAL(KIND = dp), PARAMETER :: beta = 0.0042_dp

!   *** Local variables ***

    REAL(KIND = dp) :: drhos,f,rho13,rho43,rhos13,rhos43,rhos,s,s2,y,z

!   ---------------------------------------------------------------------------

    IF (rho < eps_rho) THEN

      ex = 0.0_dp
      vx = 0.0_dp
      vxg = 0.0_dp

    ELSE

      rhos = 0.5_dp*rho
      drhos = 0.5_dp*drho

      rhos13 = rhos**f13
      rhos43 = rhos13*rhos

      rho13 = rho**f13
      rho43 = rho13*rho

!     *** LSDA exchange part (VWN) ***

      ex = cx_vwn_e*rho43
      vx = cx_vwn_v*rho13

!     *** GGA part ***

      s = drhos/rhos43
      s2 = s*s

!     *** Calculate the area sinus hyperbolicus: y = arsinh(s) ***

      y = LOG(s + SQRT(s2 + 1.0_dp))

      f = 1.0_dp/(1.0_dp + 6.0_dp*beta*s*y)
      z = 6.0_dp*beta*s*f*(y + s/SQRT(1.0_dp + s2))

      ex = ex - 2.0_dp*beta*drhos*s*f
      vx = vx - f43*beta*rhos13*s2*f*(z - 1.0_dp)
      vxg = -0.5_dp*beta*f*(2.0_dp - z)/rhos43

    END IF

  END SUBROUTINE becke88_0

! *****************************************************************************

  SUBROUTINE becke88_1(rho,drho,ex,vx,vxg,weight)

    REAL(KIND = dp), INTENT(INOUT)                      :: ex
    REAL(KIND = dp), DIMENSION(:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:), INTENT(INOUT)        :: vx,vxg
    REAL(KIND = dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0,vxg0
    INTEGER  :: i,n1

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)

    IF (PRESENT(weight)) THEN
      DO i=1,n1
        CALL becke88_0(rho(i),drho(i),ex0,vx0,vxg0)
        ex = ex + weight(i)*ex0
        vx(i) = vx(i) + weight(i)*vx0
        vxg(i) = vxg(i) + weight(i)*vxg0
      END DO
    ELSE
      DO i=1,n1
        CALL becke88_0(rho(i),drho(i),ex0,vx0,vxg0)
        ex = ex + ex0
        vx(i) = vx(i) + vx0
        vxg(i) = vxg(i) + vxg0
      END DO
    END IF

  END SUBROUTINE becke88_1

! *****************************************************************************

  SUBROUTINE becke88_2(rho,drho,ex,vx,vxg,weight)

    REAL(KIND = dp), INTENT(INOUT)                        :: ex
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:,:), INTENT(inOUT)        :: vx,vxg
    REAL(KIND = dp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0,vxg0
    INTEGER  :: i,j,n1,n2

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)

    IF (PRESENT(weight)) THEN
      DO i=1,n2
        DO j=1,n1
          CALL becke88_0(rho(j,i),drho(j,i),ex0,vx0,vxg0)
          ex = ex + weight(j,i)*ex0
          vx(j,i) = vx(j,i) + weight(j,i)*vx0
          vxg(j,i) = vxg(j,i) + weight(j,i)*vxg0
        END DO
      END DO
    ELSE
      DO i=1,n2
        DO j=1,n1
          CALL becke88_0(rho(j,i),drho(j,i),ex0,vx0,vxg0)
          ex = ex + ex0
          vx(j,i) = vx(j,i) + vx0
          vxg(j,i) = vxg(j,i) + vxg0
        END DO
      END DO
    END IF

  END SUBROUTINE becke88_2

! *****************************************************************************

  SUBROUTINE becke88_3(rho,drho,ex,vx,vxg,weight)

    REAL(KIND = dp), INTENT(INOUT)                          :: ex
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(INOUT)        :: vx,vxg
    REAL(KIND = dp), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0,vxg0
    INTEGER  :: i,j,k,n1,n2,n3

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)
    n3 = SIZE(rho,3)

    IF (PRESENT(weight)) THEN
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL becke88_0(rho(k,j,i),drho(k,j,i),ex0,vx0,vxg0)
            ex = ex + weight(k,j,i)*ex0
            vx(k,j,i) = vx(k,j,i) + weight(k,j,i)*vx0
            vxg(k,j,i) = vxg(k,j,i) + weight(k,j,i)*vxg0
          END DO
        END DO
      END DO
    ELSE
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL becke88_0(rho(k,j,i),drho(k,j,i),ex0,vx0,vxg0)
            ex = ex + ex0
            vx(k,j,i) = vx(k,j,i) + vx0
            vxg(k,j,i) = vxg(k,j,i) + vxg0
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE becke88_3

! *****************************************************************************

  SUBROUTINE becke88_sp(rhos,drhos,exs,vxs,vxgs)

!   Purpose: Calculate the gradient-corrected exchange energy and potential
!            of Becke (1988) for an open shell density (spin-polarized).

!   Literature: A. D. Becke, Phys. Rev. A 38, 3098 (1988)

!   History: - Creation (08.03.2000, Matthias Krack)

!   ***************************************************************************

    REAL(KIND = dp), INTENT(IN)  :: drhos,rhos
    REAL(KIND = dp), INTENT(OUT) :: exs,vxs,vxgs

!   *** Local parameters ***

    REAL(KIND = dp), PARAMETER :: beta = 0.0042_dp

!   *** Local variables ***

    REAL(KIND = dp) :: f,rhos13,rhos43,s,s2,y,z

!   ---------------------------------------------------------------------------

    rhos13 = rhos**f13
    rhos43 = rhos13*rhos

    s = drhos/rhos43
    s2 = s*s

!   *** Calculate the area sinus hyperbolicus: y = arsinh(s) ***

    y = LOG(s + SQRT(s2 + 1.0_dp))

    f = 1.0_dp/(1.0_dp + 6.0_dp*beta*s*y)
    z = 6.0_dp*beta*s*f*(y + s/SQRT(1.0_dp + s2))

    exs = -beta*drhos*s*f
    vxs = -f43*beta*rhos13*s2*f*(z - 1.0_dp)
    vxgs = -beta*f*(2.0_dp - z)/rhos43

  END SUBROUTINE becke88_sp

! *****************************************************************************

  SUBROUTINE becke88_sp_0(rhoa,rhob,drhoa,drhob,ex,vxa,vxb,vxga,vxgb)

    REAL(KIND = dp), INTENT(IN)  :: drhoa,drhob,rhoa,rhob
    REAL(KIND = dp), INTENT(OUT) :: ex,vxa,vxb,vxga,vxgb

!   *** Local variables ***

    REAL(KIND = dp) :: exa,exb

!   ---------------------------------------------------------------------------

    IF (rhoa < eps_rho) THEN
      exa = 0.0_dp
      vxa = 0.0_dp
      vxga = 0.0_dp
    ELSE
      CALL becke88_sp(rhoa,drhoa,exa,vxa,vxga)
    END IF

    IF (rhob < eps_rho) THEN
      exb = 0.0_dp
      vxb = 0.0_dp
      vxgb = 0.0_dp
    ELSE
      CALL becke88_sp(rhob,drhob,exb,vxb,vxgb)
    END IF

    ex = exa + exb

  END SUBROUTINE becke88_sp_0

! *****************************************************************************

  SUBROUTINE becke88_sp_1(rhoa,rhob,drhoa,drhob,ex,vxa,vxb,vxga,vxgb,weight)

    REAL(KIND = dp), INTENT(INOUT)                      :: ex
    REAL(KIND = dp), DIMENSION(:), INTENT(IN)           :: drhoa,drhob,rhoa,rhob
    REAL(KIND = dp), DIMENSION(:), INTENT(INOUT)        :: vxa,vxb,vxga,vxgb
    REAL(KIND = dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vxa0,vxb0,vxga0,vxgb0
    INTEGER  :: i,n1

!   ---------------------------------------------------------------------------

    n1 = MIN(SIZE(rhoa,1),SIZE(rhob,1),&
             SIZE(drhoa,1),SIZE(drhob,1))

    IF (PRESENT(weight)) THEN
      DO i=1,n1
        CALL becke88_sp_0(rhoa(i),rhob(i),drhoa(i),drhob(i),&
                          ex0,vxa0,vxb0,vxga0,vxgb0)
        ex = ex + weight(i)*ex0
        vxa(i) = vxa(i) + weight(i)*vxa0
        vxb(i) = vxb(i) + weight(i)*vxb0
        vxga(i) = vxga(i) + weight(i)*vxga0
        vxgb(i) = vxgb(i) + weight(i)*vxgb0
      END DO
    ELSE
      DO i=1,n1
        CALL becke88_sp_0(rhoa(i),rhob(i),drhoa(i),drhob(i),&
                          ex0,vxa0,vxb0,vxga0,vxgb0)
        ex = ex + ex0
        vxa(i) = vxa(i) + vxa0
        vxb(i) = vxb(i) + vxb0
        vxga(i) = vxga(i) + vxga0
        vxgb(i) = vxgb(i) + vxgb0
      END DO
    END IF

  END SUBROUTINE becke88_sp_1

! *****************************************************************************

  SUBROUTINE becke88_sp_2(rhoa,rhob,drhoa,drhob,ex,vxa,vxb,vxga,vxgb,weight)

    REAL(KIND = dp), INTENT(INOUT)                        :: ex
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)           :: drhoa,drhob,rhoa,rhob
    REAL(KIND = dp), DIMENSION(:,:), INTENT(INOUT)        :: vxa,vxb,vxga,vxgb
    REAL(KIND = dp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vxa0,vxb0,vxga0,vxgb0
    INTEGER  :: i,j,n1,n2

!   ---------------------------------------------------------------------------

    n1 = MIN(SIZE(rhoa,1),SIZE(rhob,1),&
             SIZE(drhoa,1),SIZE(drhob,1))
    n2 = MIN(SIZE(rhoa,2),SIZE(rhob,2),&
             SIZE(drhoa,2),SIZE(drhob,2))

    IF (PRESENT(weight)) THEN
      DO i=1,n2
        DO j=1,n1
          CALL becke88_sp_0(rhoa(j,i),rhob(j,i),drhoa(j,i),drhob(j,i),&
                            ex0,vxa0,vxb0,vxga0,vxgb0)
          ex = ex + weight(j,i)*ex0
          vxa(j,i) = vxa(j,i) + weight(j,i)*vxa0
          vxb(j,i) = vxb(j,i) + weight(j,i)*vxb0
          vxga(j,i) = vxga(j,i) + weight(j,i)*vxga0
          vxgb(j,i) = vxgb(j,i) + weight(j,i)*vxgb0
        END DO
      END DO
    ELSE
      DO i=1,n2
        DO j=1,n1
          CALL becke88_sp_0(rhoa(j,i),rhob(j,i),drhoa(j,i),drhob(j,i),&
                            ex0,vxa0,vxb0,vxga0,vxgb0)
          ex = ex + ex0
          vxa(j,i) = vxa(j,i) + vxa0
          vxb(j,i) = vxb(j,i) + vxb0
          vxga(j,i) = vxga(j,i) + vxga0
          vxgb(j,i) = vxgb(j,i) + vxgb0
        END DO
      END DO
    END IF

  END SUBROUTINE becke88_sp_2

! *****************************************************************************

  SUBROUTINE becke88_sp_3(rhoa,rhob,drhoa,drhob,ex,vxa,vxb,vxga,vxgb,weight)

    REAL(KIND = dp), INTENT(INOUT)                          :: ex
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(IN)           :: drhoa,drhob,rhoa,rhob
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(INOUT)        :: vxa,vxb,vxga,vxgb
    REAL(KIND = dp), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vxa0,vxb0,vxga0,vxgb0
    INTEGER  :: i,j,k,n1,n2,n3

!   ---------------------------------------------------------------------------

    n1 = MIN(SIZE(rhoa,1),SIZE(rhob,1),&
             SIZE(drhoa,1),SIZE(drhob,1))
    n2 = MIN(SIZE(rhoa,2),SIZE(rhob,2),&
             SIZE(drhoa,2),SIZE(drhob,2))
    n3 = MIN(SIZE(rhoa,3),SIZE(rhob,3),&
             SIZE(drhoa,3),SIZE(drhob,3))

    IF (PRESENT(weight)) THEN
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL becke88_sp_0(rhoa(k,j,i),rhob(k,j,i),drhoa(k,j,i),&
                              drhob(k,j,i),ex0,vxa0,vxb0,vxga0,vxgb0)
            ex = ex + weight(k,j,i)*ex0
            vxa(k,j,i) = vxa(k,j,i) + weight(k,j,i)*vxa0
            vxb(k,j,i) = vxb(k,j,i) + weight(k,j,i)*vxb0
            vxga(k,j,i) = vxga(k,j,i) + weight(k,j,i)*vxga0
            vxgb(k,j,i) = vxgb(k,j,i) + weight(k,j,i)*vxgb0
          END DO
        END DO
      END DO
    ELSE
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL becke88_sp_0(rhoa(k,j,i),rhob(k,j,i),drhoa(k,j,i),&
                              drhob(k,j,i),ex0,vxa0,vxb0,vxga0,vxgb0)
            ex = ex + ex0
            vxa(k,j,i) = vxa(k,j,i) + vxa0
            vxb(k,j,i) = vxb(k,j,i) + vxb0
            vxga(k,j,i) = vxga(k,j,i) + vxga0
            vxgb(k,j,i) = vxgb(k,j,i) + vxgb0
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE becke88_sp_3

! *****************************************************************************

  SUBROUTINE hcth_0(iparset,rho,drho,exc,vxc,vxcg)

!   Purpose: Calculate the gradient-corrected exchange energy and potential
!            of Hamprecht, Cohen, Tozer, and Handy (HCTH) for a closed shell
!            density.

!   Literature: - F. A. Hamprecht, A. J. Cohen, D. J. Tozer, and N. C. Handy,
!                 J. Chem. Phys. 109, 6264 (1998) -> HCTH/93
!               - A. D. Boese, N. L. Doltsinis, N. C. Handy, and M. Sprik,
!                 J. Chem. Phys. 112, 1670 (2000) -> HCTH/120 and HCTH/147
!               - A. D. Boese and N. C. Handy,
!                 J. Chem. Phys. 114, 5497 (2001) -> HCTH/407
!               - J. P. Perdew and Y. Wang,
!                 Phys. Rev. B 45, 13244 (1992) -> PW92

!   History: - Creation (27.07.2001, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN)   :: iparset
    REAL(KIND = dp), INTENT(IN)  :: drho,rho
    REAL(KIND = dp), INTENT(OUT) :: exc,vxc,vxcg

!   *** Local parameters ***

!   *** LSDA correlation parametrisation (PW92) ***

    REAL(KIND = dp), PARAMETER :: a0 = 0.031091_dp,&
                           a1 = 0.015545_dp,&
                           alpha0 = 0.21370_dp,&
                           alpha1 = 0.20548_dp

    REAL(KIND = dp), DIMENSION(4), PARAMETER :: beta0 = (/ 7.59570_dp,&
                                                    3.58760_dp,&
                                                    1.63820_dp,&
                                                    0.49294_dp/),&
                                         beta1 = (/14.11890_dp,&
                                                    6.19770_dp,&
                                                    3.36620_dp,&
                                                    0.62517_dp/)

!   *** GGA parametrisation (HCTH/iparset) ***

    REAL(KIND = dp), PARAMETER :: gamma_cab = 0.006_dp,&
                           gamma_css = 0.200_dp,&
                           gamma_xss = 0.004_dp

!   *** Local variables ***

    REAL(KIND = dp) :: dgcabddrho,dgcabdrho,dgcabds,dgcssddrho,dgcssdrho,dgcssds,&
                dgdrs,dgxssddrho,dgxssdrho,dgxssds,drsdrho,drhos,ecab,ecss,&
                exss,g,gcab,gcss,gs2,gxss,p,q,rho13,rho43,rhos13,rhos43,rhos,&
                rs,rs12,s,s2,u,vcab,vcss,vxss,x,y

    REAL(KIND = dp), DIMENSION(0:4) :: ccab,ccss,cxss

!   ---------------------------------------------------------------------------

    IF (rho < eps_rho) THEN

      exc = 0.0_dp
      vxc = 0.0_dp
      vxcg = 0.0_dp

    ELSE

!     *** Load the HCTH parameter set HCTH/iparset ***

      SELECT CASE (iparset)
      CASE (93)
        cxss(0) =  0.109320E+01_dp
        ccss(0) =  0.222601E+00_dp
        ccab(0) =  0.729974E+00_dp
        cxss(1) = -0.744056E+00_dp
        ccss(1) = -0.338622E-01_dp
        ccab(1) =  0.335287E+01_dp
        cxss(2) =  0.559920E+01_dp
        ccss(2) = -0.125170E-01_dp
        ccab(2) = -0.115430E+02_dp
        cxss(3) = -0.678549E+01_dp
        ccss(3) = -0.802496E+00_dp
        ccab(3) =  0.808564E+01_dp
        cxss(4) =  0.449357E+01_dp
        ccss(4) =  0.155396E+01_dp
        ccab(4) = -0.447857E+01_dp
      CASE (120)
        cxss(0) =  0.109163E+01_dp
        ccss(0) =  0.489508E+00_dp
        ccab(0) =  0.514730E+00_dp
        cxss(1) = -0.747215E+00_dp
        ccss(1) = -0.260699E+00_dp
        ccab(1) =  0.692982E+01_dp
        cxss(2) =  0.507833E+01_dp
        ccss(2) =  0.432917E+00_dp
        ccab(2) = -0.247073E+02_dp
        cxss(3) = -0.410746E+01_dp
        ccss(3) = -0.199247E+01_dp
        ccab(3) =  0.231098E+02_dp
        cxss(4) =  0.117173E+01_dp
        ccss(4) =  0.248531E+01_dp
        ccab(4) = -0.113234E+02_dp
      CASE (147)
        cxss(0) =  0.109025E+01_dp
        ccss(0) =  0.562576E+00_dp
        ccab(0) =  0.542352E+00_dp
        cxss(1) = -0.799194E+00_dp
        ccss(1) =  0.171436E-01_dp
        ccab(1) =  0.701464E+01_dp
        cxss(2) =  0.557212E+01_dp
        ccss(2) = -0.130636E+01_dp
        ccab(2) = -0.283822E+02_dp
        cxss(3) = -0.586760E+01_dp
        ccss(3) =  0.105747E+01_dp
        ccab(3) =  0.350329E+02_dp
        cxss(4) =  0.304544E+01_dp
        ccss(4) =  0.885429E+00_dp
        ccab(4) = -0.204284E+02_dp
      CASE (407)
        cxss(0) =  0.108184E+01_dp
        ccss(0) =  0.118777E+01_dp
        ccab(0) =  0.589076E+00_dp
        cxss(1) = -0.518339E+00_dp
        ccss(1) = -0.240292E+01_dp
        ccab(1) =  0.442374E+01_dp
        cxss(2) =  0.342562E+01_dp
        ccss(2) =  0.561741E+01_dp
        ccab(2) = -0.192218E+02_dp
        cxss(3) = -0.262901E+01_dp
        ccss(3) = -0.917923E+01_dp
        ccab(3) =  0.425721E+02_dp
        cxss(4) =  0.228855E+01_dp
        ccss(4) =  0.624798E+01_dp
        ccab(4) = -0.420052E+02_dp
      CASE DEFAULT
        STOP "in SUBROUTINE hcth_0: Invalid HCTH parameter set requested"
      END SELECT

!     *** rho_sigma = rho/2 = rho_alpha = rho_beta (same for |nabla rho|) ***

      rhos = 0.5_dp*rho
      drhos = 0.5_dp*drho

      rhos13 = rhos**f13
      rhos43 = rhos13*rhos

      rho13 = two13*rhos13
      rho43 = rho13*rho

!     *** LSDA exchange part (VWN) ***

      exss = cx_vwn_e*rho43
      vxss = cx_vwn_v*rho13

!     *** LSDA correlation part (PW92) ***

!     *** G(rho_sigma,0) => spin polarisation zeta = 1 ***

      rs = rsfac/rhos13
      rs12 = SQRT(rs)
      q = 2.0_dp*a1*(beta1(1) + (beta1(2) + (beta1(3) +&
                     beta1(4)*rs12)*rs12)*rs12)*rs12
      p = 1.0_dp + 1.0_dp/q
      x = -2.0_dp*a1*(1.0_dp + alpha1*rs)
      y = LOG(p)
      g = x*y
      dgdrs = -2.0_dp*a1*alpha1*y -&
               x*a1*(beta1(1)/rs12 + 2.0_dp*beta1(2) +&
                     3.0_dp*beta1(3)*rs12 + 4.0_dp*beta1(4)*rs)/(p*q*q)
      drsdrho = -f13*rs/rho
      ecss = rho*g
      vcss = g + rho*dgdrs*drsdrho

!     *** G(rho_alpha,rho_beta) => spin polarisation zeta = 0 ***

      rs = rsfac/rho13
      rs12 = SQRT(rs)
      q = 2.0_dp*a0*(beta0(1) + (beta0(2) + (beta0(3) +&
                     beta0(4)*rs12)*rs12)*rs12)*rs12
      p = 1.0_dp + 1.0_dp/q
      x = -2.0_dp*a0*(1.0_dp + alpha0*rs)
      y = LOG(p)
      g = x*y
      dgdrs = -2.0_dp*a0*alpha0*y -&
               x*a0*(beta0(1)/rs12 + 2.0_dp*beta0(2) +&
                     3.0_dp*beta0(3)*rs12 + 4.0_dp*beta0(4)*rs)/(p*q*q)
      drsdrho = -f13*rs/rho
      ecab = rho*g - ecss
      vcab = g + rho*dgdrs*drsdrho - vcss

!     *** GGA part (HCTH) ***

      s = drhos/rhos43
      s2 = s*s
      x = -f83/rho
      y = 2.0_dp/(drho*drho)

!     *** g_x(rho_sigma,rho_sigma) ***

      gs2 = gamma_xss*s2
      q = 1.0_dp/(1.0_dp + gs2)
      u = gs2*q
      gxss = cxss(0) + (cxss(1) + (cxss(2) + (cxss(3) + cxss(4)*u)*u)*u)*u
      dgxssds = q*(cxss(1) + (2.0_dp*cxss(2) + (3.0_dp*cxss(3) +&
                   4.0_dp*cxss(4)*u)*u)*u)*u
      dgxssdrho = x*dgxssds
      dgxssddrho = y*dgxssds

!     *** g_c(rho_sigma,rho_sigma) ***

      gs2 = gamma_css*s2
      q = 1.0_dp/(1.0_dp + gs2)
      u = gs2*q
      gcss = ccss(0) + (ccss(1) + (ccss(2) + (ccss(3) + ccss(4)*u)*u)*u)*u
      dgcssds = q*(ccss(1) + (2.0_dp*ccss(2) + (3.0_dp*ccss(3) +&
                   4.0_dp*ccss(4)*u)*u)*u)*u
      dgcssdrho = x*dgcssds
      dgcssddrho = y*dgcssds

!     *** g_c(rho_alpha,rho_beta) ***

      gs2 = gamma_cab*s2
      q = 1.0_dp/(1.0_dp + gs2)
      u = gs2*q
      gcab = ccab(0) + (ccab(1) + (ccab(2) + (ccab(3) + ccab(4)*u)*u)*u)*u
      dgcabds = q*(ccab(1) + (2.0_dp*ccab(2) + (3.0_dp*ccab(3) +&
                   4.0_dp*ccab(4)*u)*u)*u)*u
      dgcabdrho = x*dgcabds
      dgcabddrho = y*dgcabds

!     *** Finally collect all contributions ***

      exc = exss*gxss + ecss*gcss + ecab*gcab
      vxc = vxss*gxss + exss*dgxssdrho +&
            vcss*gcss + ecss*dgcssdrho +&
            vcab*gcab + ecab*dgcabdrho
      vxcg = exss*dgxssddrho + ecss*dgcssddrho + ecab*dgcabddrho

    END IF

  END SUBROUTINE hcth_0

! *****************************************************************************

  SUBROUTINE hcth_1(iparset,rho,drho,exc,vxc,vxcg,weight)

    INTEGER, INTENT(IN)                          :: iparset
    REAL(KIND = dp), INTENT(INOUT)                      :: exc
    REAL(KIND = dp), DIMENSION(:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:), INTENT(INOUT)        :: vxc,vxcg
    REAL(KIND = dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: exc0,vxc0,vxcg0
    INTEGER  :: i,n1

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)

    IF (PRESENT(weight)) THEN
      DO i=1,n1
        CALL hcth_0(iparset,rho(i),drho(i),exc0,vxc0,vxcg0)
        exc = exc + weight(i)*exc0
        vxc(i) = vxc(i) + weight(i)*vxc0
        vxcg(i) = vxcg(i) + weight(i)*vxcg0
      END DO
    ELSE
      DO i=1,n1
        CALL hcth_0(iparset,rho(i),drho(i),exc0,vxc0,vxcg0)
        exc = exc + exc0
        vxc(i) = vxc(i) + vxc0
        vxcg(i) = vxcg(i) + vxcg0
      END DO
    END IF

  END SUBROUTINE hcth_1

! *****************************************************************************

  SUBROUTINE hcth_2(iparset,rho,drho,exc,vxc,vxcg,weight)

    INTEGER, INTENT(IN)                            :: iparset
    REAL(KIND = dp), INTENT(INOUT)                        :: exc
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:,:), INTENT(INOUT)        :: vxc,vxcg
    REAL(KIND = dp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: exc0,vxc0,vxcg0
    INTEGER  :: i,j,n1,n2

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)

    IF (PRESENT(weight)) THEN
      DO i=1,n2
        DO j=1,n1
          CALL hcth_0(iparset,rho(j,i),drho(j,i),exc0,vxc0,vxcg0)
          exc = exc + weight(j,i)*exc0
          vxc(j,i) = vxc(j,i) + weight(j,i)*vxc0
          vxcg(j,i) = vxcg(j,i) + weight(j,i)*vxcg0
        END DO
      END DO
    ELSE
      DO i=1,n2
        DO j=1,n1
          CALL hcth_0(iparset,rho(j,i),drho(j,i),exc0,vxc0,vxcg0)
          exc = exc + exc0
          vxc(j,i) = vxc(j,i) + vxc0
          vxcg(j,i) = vxcg(j,i) + vxcg0
        END DO
      END DO
    END IF

  END SUBROUTINE hcth_2

! *****************************************************************************

  SUBROUTINE hcth_3(iparset,rho,drho,exc,vxc,vxcg,weight)

    INTEGER, INTENT(IN)                              :: iparset
    REAL(KIND = dp), INTENT(INOUT)                          :: exc
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(INOUT)        :: vxc,vxcg
    REAL(KIND = dp), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: exc0,vxc0,vxcg0
    INTEGER  :: i,j,k,n1,n2,n3

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)
    n3 = SIZE(rho,3)

    IF (PRESENT(weight)) THEN
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL hcth_0(iparset,rho(k,j,i),drho(k,j,i),exc0,vxc0,vxcg0)
            exc = exc + weight(k,j,i)*exc0
            vxc(k,j,i) = vxc(k,j,i) + weight(k,j,i)*vxc0
            vxcg(k,j,i) = vxcg(k,j,i) + weight(k,j,i)*vxcg0
          END DO
        END DO
      END DO
    ELSE
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL hcth_0(iparset,rho(k,j,i),drho(k,j,i),exc0,vxc0,vxcg0)
            exc = exc + exc0
            vxc(k,j,i) = vxc(k,j,i) + vxc0
            vxcg(k,j,i) = vxcg(k,j,i) + vxcg0
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE hcth_3

! *****************************************************************************

  SUBROUTINE lyp_0(rho,drho,ec,vc,vcg)

!   Purpose: Calculate the gradient-corrected correlation energy and potential
!            of Lee, Yang, and Parr (LYP, 1988) for a closed shell density.

!   Literature: - C. Lee, W. Yang, R. G. Parr,
!                 Phys. Rev. B 37, 785 (1988)
!               - B. Miehlich, A. Savin, H. Stoll, H. Preuss,
!                 Chem. Phys. Lett. 157, 200 (1989)

!   History: - Creation (11.03.2000, Matthias Krack)

!   ***************************************************************************

    REAL(KIND = dp), INTENT(IN)  :: drho,rho
    REAL(KIND = dp), INTENT(OUT) :: ec,vc,vcg

!   *** Local parameters ***

    REAL(KIND = dp), PARAMETER :: a = 0.04918_dp,&
                           b = 0.13200_dp,&
                           c = 0.25330_dp,&
                           d = 0.34900_dp

!   *** Local variables ***

    REAL(KIND = dp) :: crhom13,ddelta,delta,domega,doq,drho2,f,omega,q,rho13,rho53,&
                rho83,rhom13,rhom23

!   ---------------------------------------------------------------------------

    IF (rho < eps_rho) THEN

      ec = 0.0_dp
      vc = 0.0_dp
      vcg = 0.0_dp

    ELSE

      rho13 = rho**f13
      rho53 = rho13*rho13*rho
      rho83 = rho53*rho

      rhom13 = 1.0_dp/rho13
      rhom23 = rhom13*rhom13

      drho2 = drho*drho

      crhom13 = c*rhom13
      q = 1.0_dp + d*rhom13
      doq = d/q

      delta = crhom13 + doq*rhom13
      ddelta = f13*(doq*doq*rhom23 - delta)/rho

      f = 3.0_dp + 7.0_dp*delta

      omega = EXP(-crhom13)/(q*rho53)
      domega = f13*omega*((c + doq)*rhom13 - 5.0_dp)/rho

      ec = a*(b*omega*(f*drho2/72.0_dp - cf*rho83) - rho/q)
      vc = -a*((96.0_dp*d*rhom13 + 72.0_dp)/(q*q) -&
               b*((7.0_dp*omega*ddelta + f*domega)*drho2 -&
                  8.0_dp*cf*rho53*(24.0_dp*omega + 9.0_dp*rho*domega)))/72.0_dp
      vcg = a*b*omega*f/36.0_dp

    END IF

  END SUBROUTINE lyp_0

! *****************************************************************************

  SUBROUTINE lyp_1(rho,drho,ec,vc,vcg,weight)

    REAL(KIND = dp), INTENT(INOUT)                      :: ec
    REAL(KIND = dp), DIMENSION(:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:), INTENT(INOUT)        :: vc,vcg
    REAL(KIND = dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ec0,vc0,vcg0
    INTEGER  :: i,n1

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)

    IF (PRESENT(weight)) THEN
      DO i=1,n1
        CALL lyp_0(rho(i),drho(i),ec0,vc0,vcg0)
        ec = ec + weight(i)*ec0
        vc(i) = vc(i) + weight(i)*vc0
        vcg(i) = vcg(i) + weight(i)*vcg0
      END DO
    ELSE
      DO i=1,n1
        CALL lyp_0(rho(i),drho(i),ec0,vc0,vcg0)
        ec = ec + ec0
        vc(i) = vc(i) + vc0
        vcg(i) = vcg(i) + vcg0
      END DO
    END IF

  END SUBROUTINE lyp_1

! *****************************************************************************

  SUBROUTINE lyp_2(rho,drho,ec,vc,vcg,weight)

    REAL(KIND = dp), INTENT(INOUT)                        :: ec
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:,:), INTENT(INOUT)        :: vc,vcg
    REAL(KIND = dp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ec0,vc0,vcg0
    INTEGER  :: i,j,n1,n2

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)

    IF (PRESENT(weight)) THEN
      DO i=1,n2
        DO j=1,n1
          CALL lyp_0(rho(j,i),drho(j,i),ec0,vc0,vcg0)
          ec = ec + weight(j,i)*ec0
          vc(j,i) = vc(j,i) + weight(j,i)*vc0
          vcg(j,i) = vcg(j,i) + weight(j,i)*vcg0
        END DO
      END DO
    ELSE
      DO i=1,n2
        DO j=1,n1
          CALL lyp_0(rho(j,i),drho(j,i),ec0,vc0,vcg0)
          ec = ec + ec0
          vc(j,i) = vc(j,i) + vc0
          vcg(j,i) = vcg(j,i) + vcg0
        END DO
      END DO
    END IF

  END SUBROUTINE lyp_2

! *****************************************************************************

  SUBROUTINE lyp_3(rho,drho,ec,vc,vcg,weight)

    REAL(KIND = dp), INTENT(INOUT)                          :: ec
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(INOUT)        :: vc,vcg
    REAL(KIND = dp), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ec0,vc0,vcg0
    INTEGER  :: i,j,k,n1,n2,n3

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)
    n3 = SIZE(rho,3)

    IF (PRESENT(weight)) THEN
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL lyp_0(rho(k,j,i),drho(k,j,i),ec0,vc0,vcg0)
            ec = ec + weight(k,j,i)*ec0
            vc(k,j,i) = vc(k,j,i) + weight(k,j,i)*vc0
            vcg(k,j,i) = vcg(k,j,i) + weight(k,j,i)*vcg0
          END DO
        END DO
      END DO
    ELSE
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL lyp_0(rho(k,j,i),drho(k,j,i),ec0,vc0,vcg0)
            ec = ec + ec0
            vc(k,j,i) = vc(k,j,i) + vc0
            vcg(k,j,i) = vcg(k,j,i) + vcg0
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE lyp_3

! *****************************************************************************

  SUBROUTINE pade_0(rho,exc,vxc)

!   Purpose: Calculate the exchange-correlation energy and potential using the
!            Pade approximation of Goedecker, Teter and Hutter for LDA.

!   Literature: S. Goedecker, M. Teter and J. Hutter,
!               Phys. Rev. B 54, 1703 (1996)

!   History: - Creation (16.03.2000, Matthias Krack)

!   ***************************************************************************

    REAL(KIND = dp), INTENT(IN)  :: rho
    REAL(KIND = dp), INTENT(OUT) :: exc,vxc

!   *** Local parameters ***

    REAL(KIND = dp), PARAMETER :: a0 = 0.4581652932831429E+0_dp,&
                           a1 = 0.2217058676663745E+1_dp,&
                           a2 = 0.7405551735357053E+0_dp,&
                           a3 = 0.1968227878617998E-1_dp,&
                           b1 = 1.0000000000000000E+0_dp,&
                           b2 = 0.4504130959426697E+1_dp,&
                           b3 = 0.1110667363742916E+1_dp,&
                           b4 = 0.2359291751427506E-1_dp

!   *** Local variables ***

    REAL(KIND = dp) :: mydp,dq,epade,rs,p,q,rho13

!   ---------------------------------------------------------------------------

    IF (rho < eps_rho) THEN

      exc = 0.0_dp
      vxc = 0.0_dp

    ELSE

      rho13 = rho**f13

      rs = rsfac/rho13

      p = a0 + (a1 + (a2 + a3*rs)*rs)*rs
      q = (b1 + (b2 + (b3 + b4*rs)*rs)*rs)*rs
      epade = -p/q

      mydp = a1 + (2.0_dp*a2 + 3.0_dp*a3*rs)*rs
      dq = b1 + (2.0_dp*b2 + (3.0_dp*b3 + 4.0_dp*b4*rs)*rs)*rs

      exc = epade*rho
      vxc = epade + f13*rs*(mydp*q - p*dq)/q**2

    END IF

  END SUBROUTINE pade_0

! *****************************************************************************

  SUBROUTINE pade_1(rho,exc,vxc,weight)

    REAL(KIND = dp), INTENT(INOUT)                      :: exc
    REAL(KIND = dp), DIMENSION(:), INTENT(IN)           :: rho
    REAL(KIND = dp), DIMENSION(:), INTENT(INOUT)        :: vxc
    REAL(KIND = dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: exc0,vxc0
    INTEGER  :: i,n1

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)

    IF (PRESENT(weight)) THEN
      DO i=1,n1
        CALL pade_0(rho(i),exc0,vxc0)
        exc = exc + weight(i)*exc0
        vxc(i) = vxc(i) + weight(i)*vxc0
      END DO
    ELSE
      DO i=1,n1
        CALL pade_0(rho(i),exc0,vxc0)
        exc = exc + exc0
        vxc(i) = vxc(i) + vxc0
      END DO
    END IF

  END SUBROUTINE pade_1

! *****************************************************************************

  SUBROUTINE pade_2(rho,exc,vxc,weight)

    REAL(KIND = dp), INTENT(INOUT)                        :: exc
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)           :: rho
    REAL(KIND = dp), DIMENSION(:,:), INTENT(INOUT)        :: vxc
    REAL(KIND = dp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: exc0,vxc0
    INTEGER  :: i,j,n1,n2

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)

    IF (PRESENT(weight)) THEN
      DO i=1,n2
        DO j=1,n1
          CALL pade_0(rho(j,i),exc0,vxc0)
          exc = exc + weight(j,i)*exc0
          vxc(j,i) = vxc(j,i) + weight(j,i)*vxc0
        END DO
      END DO
    ELSE
      DO i=1,n2
        DO j=1,n1
          CALL pade_0(rho(j,i),exc0,vxc0)
          exc = exc + exc0
          vxc(j,i) = vxc(j,i) + vxc0
        END DO
      END DO
    ENDIF

  END SUBROUTINE pade_2

! *****************************************************************************

  SUBROUTINE pade_3(rho,exc,vxc,weight)

    REAL(KIND = dp), INTENT(INOUT)                          :: exc
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(IN)           :: rho
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(INOUT)        :: vxc
    REAL(KIND = dp), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: exc0,vxc0
    INTEGER  :: i,j,k,n1,n2,n3

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)
    n3 = SIZE(rho,3)

    IF (PRESENT(weight)) THEN
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL pade_0(rho(k,j,i),exc0,vxc0)
            exc = exc + weight(k,j,i)*exc0
            vxc(k,j,i) = vxc(k,j,i) + weight(k,j,i)*vxc0
          END DO
        END DO
      END DO
    ELSE
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL pade_0(rho(k,j,i),exc0,vxc0)
            exc = exc + exc0
            vxc(k,j,i) = vxc(k,j,i) + vxc0
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE pade_3

! *****************************************************************************

  SUBROUTINE pbe_c_0(rho,drho,ec,vc,vcg)

!   Purpose: Calculate the gradient-corrected correlation energy and potential
!            of Perdew, Burke and Ernzerhof (1996) for a closed shell density.

!   Literature: J. P. Perdew, K. Burke and M. Ernzerhof,
!               Phys. Rev. Lett., 77, 3865 (1996)

!   History: - Creation (22.03.2000, Matthias Krack)

!   ***************************************************************************

    REAL(KIND = dp), INTENT(IN)  :: drho,rho
    REAL(KIND = dp), INTENT(OUT) :: ec,vc,vcg

!   *** Local variables ***

    REAL(KIND = dp) :: a,at2,bog,dade,dedr,dhddr,dhdr,dpda,dpdr,dpdt,dtddr,dtdr,e,&
                expe,h,p,q,q2,rhom76,t,t2

!   ---------------------------------------------------------------------------

    CALL vwn_c_0(rho,ec,vc)

    IF ((rho < eps_rho).OR.(drho < eps_drho)) THEN

      vcg = 0.0_dp

    ELSE

      e = ec/rho

      rhom76 = rho**(-f76)

      t = r2ks*drho*rhom76
      t2 = t*t

      expe = EXP(-e/gamma_pbe)

      bog = beta_pbe/gamma_pbe

      a = bog/(expe - 1.0_dp)
      at2 = a*t2

      q = 1.0_dp/(1.0_dp + at2 + at2*at2)
      q2 = q*q

      p = 1.0_dp + bog*t2*(1.0_dp + at2)*q

      h = gamma_pbe*LOG(p)

      dpdt = 2.0_dp*bog*t*(1.0_dp + 2.0_dp*at2)*q2
      dtdr = -f76*t/rho
      dpda = -bog*at2*t2*t2*(2.0_dp + at2)*q2
      dade = a*a*expe/beta_pbe
      dedr = (vc - e)/rho
      dpdr = dpdt*dtdr + dpda*dade*dedr
      dhdr = gamma_pbe*dpdr/p
      dtddr = t/drho**2
      dhddr = gamma_pbe*dpdt*dtddr/p

      ec = ec + rho*h
      vc = vc + h + rho*dhdr
      vcg = rho*dhddr

    END IF

  END SUBROUTINE pbe_c_0

! *****************************************************************************

  SUBROUTINE pbe_c_1(rho,drho,ec,vc,vcg,weight)

    REAL(KIND = dp), INTENT(INOUT)                      :: ec
    REAL(KIND = dp), DIMENSION(:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:), INTENT(INOUT)        :: vc,vcg
    REAL(KIND = dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ec0,vc0,vcg0
    INTEGER  :: i,n1

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)

    IF (PRESENT(weight)) THEN
      DO i=1,n1
        CALL pbe_c_0(rho(i),drho(i),ec0,vc0,vcg0)
        ec = ec + weight(i)*ec0
        vc(i) = vc(i) + weight(i)*vc0
        vcg(i) = vcg(i) + weight(i)*vcg0
      END DO
    ELSE
      DO i=1,n1
        CALL pbe_c_0(rho(i),drho(i),ec0,vc0,vcg0)
        ec = ec + ec0
        vc(i) = vc(i) + vc0
        vcg(i) = vcg(i) + vcg0
      END DO
    END IF

  END SUBROUTINE pbe_c_1

! *****************************************************************************

  SUBROUTINE pbe_c_2(rho,drho,ec,vc,vcg,weight)

    REAL(KIND = dp), INTENT(INOUT)                        :: ec
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:,:), INTENT(INOUT)        :: vc,vcg
    REAL(KIND = dp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ec0,vc0,vcg0
    INTEGER  :: i,j,n1,n2

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)

    IF (PRESENT(weight)) THEN
      DO i=1,n2
        DO j=1,n1
          CALL pbe_c_0(rho(j,i),drho(j,i),ec0,vc0,vcg0)
          ec = ec + weight(j,i)*ec0
          vc(j,i) = vc(j,i) + weight(j,i)*vc0
          vcg(j,i) = vcg(j,i) + weight(j,i)*vcg0
        END DO
      END DO
    ELSE
      DO i=1,n2
        DO j=1,n1
          CALL pbe_c_0(rho(j,i),drho(j,i),ec0,vc0,vcg0)
          ec = ec + ec0
          vc(j,i) = vc(j,i) + vc0
          vcg(j,i) = vcg(j,i) + vcg0
        END DO
      END DO
    END IF

  END SUBROUTINE pbe_c_2

! *****************************************************************************

  SUBROUTINE pbe_c_3(rho,drho,ec,vc,vcg,weight)

    REAL(KIND = dp), INTENT(INOUT)                          :: ec
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(INOUT)        :: vc,vcg
    REAL(KIND = dp), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ec0,vc0,vcg0
    INTEGER  :: i,j,k,n1,n2,n3

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)
    n3 = SIZE(rho,3)

    IF (PRESENT(weight)) THEN
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL pbe_c_0(rho(k,j,i),drho(k,j,i),ec0,vc0,vcg0)
            ec = ec + weight(k,j,i)*ec0
            vc(k,j,i) = vc(k,j,i) + weight(k,j,i)*vc0
            vcg(k,j,i) = vcg(k,j,i) + weight(k,j,i)*vcg0
          END DO
        END DO
      END DO
    ELSE
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL pbe_c_0(rho(k,j,i),drho(k,j,i),ec0,vc0,vcg0)
            ec = ec + ec0
            vc(k,j,i) = vc(k,j,i) + vc0
            vcg(k,j,i) = vcg(k,j,i) + vcg0
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE pbe_c_3

! *****************************************************************************

  SUBROUTINE pbe_x_0(rho,drho,ex,vx,vxg)

!   Purpose: Calculate the gradient-corrected exchange energy and potential
!            of Perdew, Burke and Ernzerhof (1996) for a closed shell density.

!   Literature: J. P. Perdew, K. Burke and M. Ernzerhof,
!               Phys. Rev. Lett., 77, 3865 (1996)

!   History: - Creation (20.03.2000, Matthias Krack)

!   ***************************************************************************

    REAL(KIND = dp), INTENT(IN)  :: drho,rho
    REAL(KIND = dp), INTENT(OUT) :: ex,vx,vxg

!   *** Local parameters ***

    REAL(KIND = dp), PARAMETER :: kappa = 0.804_dp,&
                           mu = 0.2195149727645171_dp

!   *** Local variables ***

    REAL(KIND = dp) :: dfx,fx,p,rho13,rho43,s

!   ---------------------------------------------------------------------------

    IF ((rho < eps_rho).OR.(drho < eps_drho)) THEN

      ex = 0.0_dp
      vx = 0.0_dp
      vxg = 0.0_dp

    ELSE

      rho13 = rho**f13
      rho43 = rho13*rho

      s = r2kf*drho/rho43

      p = 1.0_dp/(1.0_dp + mu*s*s/kappa)

      fx = 1.0_dp + kappa*(1.0_dp - p)
      dfx = 2.0_dp*mu*s*p*p

      ex = cx_vwn_e*rho43*fx
      vx = cx_vwn_v*rho13*(fx - s*dfx)
      vxg = cx_vwn_e*r2kf*dfx/drho

    END IF

  END SUBROUTINE pbe_x_0

! *****************************************************************************

  SUBROUTINE pbe_x_1(rho,drho,ex,vx,vxg,weight)

    REAL(KIND = dp), INTENT(INOUT)                      :: ex
    REAL(KIND = dp), DIMENSION(:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:), INTENT(INOUT)        :: vx,vxg
    REAL(KIND = dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0,vxg0
    INTEGER  :: i,n1

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)

    IF (PRESENT(weight)) THEN
      DO i=1,n1
        CALL pbe_x_0(rho(i),drho(i),ex0,vx0,vxg0)
        ex = ex + weight(i)*ex0
        vx(i) = vx(i) + weight(i)*vx0
        vxg(i) = vxg(i) + weight(i)*vxg0
      END DO
    ELSE
      DO i=1,n1
        CALL pbe_x_0(rho(i),drho(i),ex0,vx0,vxg0)
        ex = ex + ex0
        vx(i) = vx(i) + vx0
        vxg(i) = vxg(i) + vxg0
      END DO
    END IF

  END SUBROUTINE pbe_x_1

! *****************************************************************************

  SUBROUTINE pbe_x_2(rho,drho,ex,vx,vxg,weight)

    REAL(KIND = dp), INTENT(INOUT)                        :: ex
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:,:), INTENT(INOUT)        :: vx,vxg
    REAL(KIND = dp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0,vxg0
    INTEGER  :: i,j,n1,n2

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)

    IF (PRESENT(weight)) THEN
      DO i=1,n2
        DO j=1,n1
          CALL pbe_x_0(rho(j,i),drho(j,i),ex0,vx0,vxg0)
          ex = ex + weight(j,i)*ex0
          vx(j,i) = vx(j,i) + weight(j,i)*vx0
          vxg(j,i) = vxg(j,i) + weight(j,i)*vxg0
        END DO
      END DO
    ELSE
      DO i=1,n2
        DO j=1,n1
          CALL pbe_x_0(rho(j,i),drho(j,i),ex0,vx0,vxg0)
          ex = ex + ex0
          vx(j,i) = vx(j,i) + vx0
          vxg(j,i) = vxg(j,i) + vxg0
        END DO
      END DO
    END IF

  END SUBROUTINE pbe_x_2

! *****************************************************************************

  SUBROUTINE pbe_x_3(rho,drho,ex,vx,vxg,weight)

    REAL(KIND = dp), INTENT(INOUT)                          :: ex
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(INOUT)        :: vx,vxg
    REAL(KIND = dp), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0,vxg0
    INTEGER  :: i,j,k,n1,n2,n3

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)
    n3 = SIZE(rho,3)

    IF (PRESENT(weight)) THEN
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL pbe_x_0(rho(k,j,i),drho(k,j,i),ex0,vx0,vxg0)
            ex = ex + weight(k,j,i)*ex0
            vx(k,j,i) = vx(k,j,i) + weight(k,j,i)*vx0
            vxg(k,j,i) = vxg(k,j,i) + weight(k,j,i)*vxg0
          END DO
        END DO
      END DO
    ELSE
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL pbe_x_0(rho(k,j,i),drho(k,j,i),ex0,vx0,vxg0)
            ex = ex + ex0
            vx(k,j,i) = vx(k,j,i) + vx0
            vxg(k,j,i) = vxg(k,j,i) + vxg0
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE pbe_x_3

! *****************************************************************************

  SUBROUTINE perdew86_c_0(rho,drho,ec,vc,vcg)

!   Purpose: Calculate the gradient-corrected correlation energy and potential
!            of Perdew (1986) for a closed shell density.

!   Literature: J. P. Perdew, Phys. Rev. B, 33, 8822 (1986)
!               Erratum: Phys. Rev. B, 34, 7406 (1986)

!   History: - Creation (16.03.2000, Matthias Krack)

!   ***************************************************************************

    REAL(KIND = dp), INTENT(IN)  :: drho,rho
    REAL(KIND = dp), INTENT(OUT) :: ec,vc,vcg

!   *** Local parameters ***

    REAL(KIND = dp), PARAMETER :: c1 = 0.001667_dp,&
                           c2 = 0.002568_dp,&
                           c3 = 0.023266_dp,&
                           c4 = 7.389E-6_dp,&
                           c5 = 8.723_dp,&
                           c6 = 0.472_dp,&
                           c7 = 1.0E4_dp*c4,&
                           f = 0.11_dp

!   *** Local variables ***

    REAL(KIND = dp) :: crho,dcrho,mydp,dq,drho2,ephi,rs,p,phi,q,rhom43,rhom76

!   ---------------------------------------------------------------------------

    IF (rho < eps_rho) THEN

      vcg = 0.0_dp

    ELSE

      rhom76 = rho**(-f76)

      rs = rsfac*rho**(-f13)

      p = c2 + (c3 + c4*rs)*rs
      q = 1.0_dp + (c5 + (c6 + c7*rs)*rs)*rs

      crho = c1 + p/q

      phi = c_perdew86*f*(c1 + c2)*drho*rhom76/crho

      IF (phi < 80.0_dp) THEN

        rhom43 = rho**(-f43)
        ephi = EXP(-phi)
        drho2 = drho*drho

        ec = ec + crho*EXP(-phi)*rhom43*drho2

        mydp = c3 + 2.0_dp*c4*rs
        dq = c5 + (2.0_dp*c6 + 3.0_dp*c7*rs)*rs

        dcrho = (mydp*q - p*dq)/(q*q)

        vc = vc - ephi*drho2*rhom76*rhom76*((f43 - f76*phi)*crho +&
                                            f13*rs*(phi + 1.0_dp)*dcrho)
        vcg = -ephi*rhom43*(phi - 2.0_dp)*crho

      ELSE

        vcg = 0.0_dp

      END IF

    END IF

  END SUBROUTINE perdew86_c_0

! *****************************************************************************

  SUBROUTINE perdew86_c_1(rho,drho,ec,vc,vcg,weight)

    REAL(KIND = dp), INTENT(INOUT)                      :: ec
    REAL(KIND = dp), DIMENSION(:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:), INTENT(INOUT)        :: vc,vcg
    REAL(KIND = dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ec0,vc0,vcg0
    INTEGER  :: i,n1

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)

    IF (PRESENT(weight)) THEN
      DO i=1,n1
        CALL perdew86_c_0(rho(i),drho(i),ec0,vc0,vcg0)
        ec = ec + weight(i)*ec0
        vc(i) = vc(i) + weight(i)*vc0
        vcg(i) = vcg(i) + weight(i)*vcg0
      END DO
    ELSE
      DO i=1,n1
        CALL perdew86_c_0(rho(i),drho(i),ec0,vc0,vcg0)
        ec = ec + ec0
        vc(i) = vc(i) + vc0
        vcg(i) = vcg(i) + vcg0
      END DO
    END IF

  END SUBROUTINE perdew86_c_1

! *****************************************************************************

  SUBROUTINE perdew86_c_2(rho,drho,ec,vc,vcg,weight)

    REAL(KIND = dp), INTENT(INOUT)                        :: ec
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:,:), INTENT(INOUT)        :: vc,vcg
    REAL(KIND = dp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ec0,vc0,vcg0
    INTEGER  :: i,j,n1,n2

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)

    IF (PRESENT(weight)) THEN
      DO i=1,n2
        DO j=1,n1
          CALL perdew86_c_0(rho(j,i),drho(j,i),ec0,vc0,vcg0)
          ec = ec + weight(j,i)*ec0
          vc(j,i) = vc(j,i) + weight(j,i)*vc0
          vcg(j,i) = vcg(j,i) + weight(j,i)*vcg0
        END DO
      END DO
    ELSE
      DO i=1,n2
        DO j=1,n1
          CALL perdew86_c_0(rho(j,i),drho(j,i),ec0,vc0,vcg0)
          ec = ec + ec0
          vc(j,i) = vc(j,i) + vc0
          vcg(j,i) = vcg(j,i) + vcg0
        END DO
      END DO
    END IF

  END SUBROUTINE perdew86_c_2

! *****************************************************************************

  SUBROUTINE perdew86_c_3(rho,drho,ec,vc,vcg,weight)

    REAL(KIND = dp), INTENT(INOUT)                          :: ec
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(INOUT)        :: vc,vcg
    REAL(KIND = dp), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ec0,vc0,vcg0
    INTEGER  :: i,j,k,n1,n2,n3

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)
    n3 = SIZE(rho,3)

    IF (PRESENT(weight)) THEN
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL perdew86_c_0(rho(k,j,i),drho(k,j,i),ec0,vc0,vcg0)
            ec = ec + weight(k,j,i)*ec0
            vc(k,j,i) = vc(k,j,i) + weight(k,j,i)*vc0
            vcg(k,j,i) = vcg(k,j,i) + weight(k,j,i)*vcg0
          END DO
        END DO
      END DO
    ELSE
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL perdew86_c_0(rho(k,j,i),drho(k,j,i),ec0,vc0,vcg0)
            ec = ec + ec0
            vc(k,j,i) = vc(k,j,i) + vc0
            vcg(k,j,i) = vcg(k,j,i) + vcg0
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE perdew86_c_3

! *****************************************************************************

  SUBROUTINE perdew86_x_0(rho,drho,ex,vx,vxg)

!   Purpose: Calculate the gradient-corrected exchange energy and potential
!            of Perdew and Wang (1986) for a closed shell density.

!   Literature: J. P. Perdew and Y. Wang, Phys. Rev. B, 33, 8800 (1986)

!   History: - Creation (17.03.2000, Matthias Krack)

!   ***************************************************************************

    REAL(KIND = dp), INTENT(IN)  :: drho,rho
    REAL(KIND = dp), INTENT(OUT) :: ex,vx,vxg

!   *** Local parameters ***

    REAL(KIND = dp), PARAMETER :: b = 14.0_dp,&
                           c = 0.2_dp,&
                           d = 0.0864_dp,&
                           m = 1.0_dp/15.0_dp,&
                           a = d/m

!   *** Local variables ***

    REAL(KIND = dp) :: dfx,mydp,fx,p,rho13,rho43,s,s2,t,tau

!   ---------------------------------------------------------------------------

    IF (rho < eps_rho) THEN

      ex = 0.0_dp
      vx = 0.0_dp
      vxg = 0.0_dp

    ELSE

      rho13 = rho**f13
      rho43 = rho13*rho

      s = r2kf*drho/rho43
      s2 = s*s

      p = 1.0_dp + (a + (b + c*s2)*s2)*s2
      mydp = 2.0_dp*s*(a + (2.0_dp*b + 3.0_dp*c*s2)*s2)

      fx = p**m
      dfx = m*p**(m - 1.0_dp)*mydp

!     *** Correction to fx for large s values ***

      IF (s > 5.0_dp) THEN
        tau = EXP(0.25_dp - 0.05_dp*s)
        t = tau*(1.0_dp - fx)
        fx = 1.0_dp - t
        dfx = 0.05_dp*t + tau*dfx
      END IF

      ex = cx_vwn_e*rho43*fx
      vx = cx_vwn_v*rho13*(fx - s*dfx)
      vxg = cx_vwn_e*r2kf*dfx/drho

    END IF

  END SUBROUTINE perdew86_x_0

! *****************************************************************************

  SUBROUTINE perdew86_x_1(rho,drho,ex,vx,vxg,weight)

    REAL(KIND = dp), INTENT(INOUT)                      :: ex
    REAL(KIND = dp), DIMENSION(:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:), INTENT(INOUT)        :: vx,vxg
    REAL(KIND = dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0,vxg0
    INTEGER  :: i,n1

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)

    IF (PRESENT(weight)) THEN
      DO i=1,n1
        CALL perdew86_x_0(rho(i),drho(i),ex0,vx0,vxg0)
        ex = ex + weight(i)*ex0
        vx(i) = vx(i) + weight(i)*vx0
        vxg(i) = vxg(i) + weight(i)*vxg0
      END DO
    ELSE
      DO i=1,n1
        CALL perdew86_x_0(rho(i),drho(i),ex0,vx0,vxg0)
        ex = ex + ex0
        vx(i) = vx(i) + vx0
        vxg(i) = vxg(i) + vxg0
      END DO
    END IF

  END SUBROUTINE perdew86_x_1

! *****************************************************************************

  SUBROUTINE perdew86_x_2(rho,drho,ex,vx,vxg,weight)

    REAL(KIND = dp), INTENT(INOUT)                        :: ex
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:,:), INTENT(INOUT)        :: vx,vxg
    REAL(KIND = dp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0,vxg0
    INTEGER  :: i,j,n1,n2

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)

    IF (PRESENT(weight)) THEN
      DO i=1,n2
        DO j=1,n1
          CALL perdew86_x_0(rho(j,i),drho(j,i),ex0,vx0,vxg0)
          ex = ex + weight(j,i)*ex0
          vx(j,i) = vx(j,i) + weight(j,i)*vx0
          vxg(j,i) = vxg(j,i) + weight(j,i)*vxg0
        END DO
      END DO
    ELSE
      DO i=1,n2
        DO j=1,n1
          CALL perdew86_x_0(rho(j,i),drho(j,i),ex0,vx0,vxg0)
          ex = ex + ex0
          vx(j,i) = vx(j,i) + vx0
          vxg(j,i) = vxg(j,i) + vxg0
        END DO
      END DO
    END IF

  END SUBROUTINE perdew86_x_2

! *****************************************************************************

  SUBROUTINE perdew86_x_3(rho,drho,ex,vx,vxg,weight)

    REAL(KIND = dp), INTENT(INOUT)                          :: ex
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(IN)           :: drho,rho
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(INOUT)        :: vx,vxg
    REAL(KIND = dp), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0,vxg0
    INTEGER  :: i,j,k,n1,n2,n3

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)
    n3 = SIZE(rho,3)

    IF (PRESENT(weight)) THEN
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL perdew86_x_0(rho(k,j,i),drho(k,j,i),ex0,vx0,vxg0)
            ex = ex + weight(k,j,i)*ex0
            vx(k,j,i) = vx(k,j,i) + weight(k,j,i)*vx0
            vxg(k,j,i) = vxg(k,j,i) + weight(k,j,i)*vxg0
          END DO
        END DO
      END DO
    ELSE
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL perdew86_x_0(rho(k,j,i),drho(k,j,i),ex0,vx0,vxg0)
            ex = ex + ex0
            vx(k,j,i) = vx(k,j,i) + vx0
            vxg(k,j,i) = vxg(k,j,i) + vxg0
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE perdew86_x_3

! *****************************************************************************

  SUBROUTINE slater_0(rho,ex,vx)

!   Purpose: Calculate the exchange energy and potential using the X_alpha
!            approximation of Slater for a closed shell density.

!   Literature: J. C. Slater, Phys. Rev. 81, 385 (1951)

!   History: - Creation (11.03.2000, Matthias Krack)

!   ***************************************************************************

    REAL(KIND = dp), INTENT(IN)  :: rho
    REAL(KIND = dp), INTENT(OUT) :: ex,vx

!   *** Local variables ***

    REAL(KIND = dp) :: rho13

!   ---------------------------------------------------------------------------

    IF (rho < eps_rho) THEN

      ex = 0.0_dp
      vx = 0.0_dp

    ELSE

      rho13 = rho**f13

      ex = cx_slater_e*rho13*rho
      vx = cx_slater_v*rho13

    END IF

  END SUBROUTINE slater_0

! *****************************************************************************

  SUBROUTINE slater_1(rho,ex,vx,weight)

    REAL(KIND = dp), INTENT(INOUT)                      :: ex
    REAL(KIND = dp), DIMENSION(:), INTENT(IN)           :: rho
    REAL(KIND = dp), DIMENSION(:), INTENT(INOUT)        :: vx
    REAL(KIND = dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0
    INTEGER  :: i,n1

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)

    IF (PRESENT(weight)) THEN
      DO i=1,n1
        CALL slater_0(rho(i),ex0,vx0)
        ex = ex + weight(i)*ex0
        vx(i) = vx(i) + weight(i)*vx0
      END DO
    ELSE
      DO i=1,n1
        CALL slater_0(rho(i),ex0,vx0)
        ex = ex + ex0
        vx(i) = vx(i) + vx0
      END DO
    END IF

  END SUBROUTINE slater_1

! *****************************************************************************

  SUBROUTINE slater_2(rho,ex,vx,weight)

    REAL(KIND = dp), INTENT(INOUT)                        :: ex
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)           :: rho
    REAL(KIND = dp), DIMENSION(:,:), INTENT(INOUT)        :: vx
    REAL(KIND = dp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0
    INTEGER  :: i,j,n1,n2

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)

    IF (PRESENT(weight)) THEN
      DO i=1,n2
        DO j=1,n1
          CALL slater_0(rho(j,i),ex0,vx0)
          ex = ex + weight(j,i)*ex0
          vx(j,i) = vx(j,i) + weight(j,i)*vx0
        END DO
      END DO
    ELSE
      DO i=1,n2
        DO j=1,n1
          CALL slater_0(rho(j,i),ex0,vx0)
          ex = ex + ex0
          vx(j,i) = vx(j,i) + vx0
        END DO
      END DO
    END IF

  END SUBROUTINE slater_2

! *****************************************************************************

  SUBROUTINE slater_3(rho,ex,vx,weight)

    REAL(KIND = dp), INTENT(INOUT)                          :: ex
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(IN)           :: rho
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(INOUT)        :: vx
    REAL(KIND = dp), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0
    INTEGER  :: i,j,k,n1,n2,n3

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)
    n3 = SIZE(rho,3)

    IF (PRESENT(weight)) THEN
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL slater_0(rho(k,j,i),ex0,vx0)
            ex = ex + weight(k,j,i)*ex0
            vx(k,j,i) = vx(k,j,i) + weight(k,j,i)*vx0
          END DO
        END DO
      END DO
    ELSE
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL slater_0(rho(k,j,i),ex0,vx0)
            ex = ex + ex0
            vx(k,j,i) = vx(k,j,i) + vx0
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE slater_3

! *****************************************************************************

  SUBROUTINE vwn_c_0(rho,ec,vc)

!   Purpose: Calculate the correlation energy and potential of Vosko, Wilk and
!            Nusair (VWN, 1980) for a closed shell density.

!   Literature: S. H. Vosko, L. Wilk and M. Nusair,
!               Can. J. Phys. 58, 1200 (1980)

!   History: - Creation (12.03.2000, Matthias Krack)

!   ***************************************************************************

    REAL(KIND = dp), INTENT(IN)  :: rho
    REAL(KIND = dp), INTENT(OUT) :: ec,vc

!   *** Local parameters (paramagnetic parameter set of VWN) ***

    REAL(KIND = dp), PARAMETER :: ap =  0.0621814_dp,&
                           bp =  3.72744_dp,&
                           cp = 12.9352_dp,&
                           xp = -0.10498_dp

!   *** Local variables ***

    REAL(KIND = dp) :: ecp,decp,qp,pp,x,x2,xmxp,xmxp2,xp2

!   ---------------------------------------------------------------------------

    IF (rho < eps_rho) THEN

      ec = 0.0_dp
      vc = 0.0_dp

    ELSE

!     *** Calculate the effective electron radius ***

      x2 = rsfac*rho**(-f13)

      x = SQRT(x2)

      xp2 = xp*xp
      xmxp = x - xp
      xmxp2 = xmxp*xmxp
      pp = x*x + bp*x + cp
      qp = SQRT(4.0_dp*cp - bp*bp)

      ecp = 0.5_dp*ap*(LOG(x2/xmxp2) -&
                       ((xp2 + cp)*LOG(pp/xmxp2) +&
                        2.0_dp*bp*(xp2 - cp)*ATAN(qp/(2.0_dp*x + bp))/qp)/&
                       (xp2 + bp*xp + cp))
      decp = 0.5_dp*f13*ap*((1.0_dp + bp/xmxp)*x2/pp - 1.0)

      ec = ecp*rho
      vc = ecp + decp

    END IF

  END SUBROUTINE vwn_c_0

! *****************************************************************************

  SUBROUTINE vwn_c_1(rho,ec,vc,weight)

    REAL(KIND = dp), INTENT(INOUT)                      :: ec
    REAL(KIND = dp), DIMENSION(:), INTENT(IN)           :: rho
    REAL(KIND = dp), DIMENSION(:), INTENT(INOUT)        :: vc
    REAL(KIND = dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ec0,vc0
    INTEGER  :: i,n1

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)

    IF (PRESENT(weight)) THEN
      DO i=1,n1
        CALL vwn_c_0(rho(i),ec0,vc0)
        ec = ec + weight(i)*ec0
        vc(i) = vc(i) + weight(i)*vc0
      END DO
    ELSE
      DO i=1,n1
        CALL vwn_c_0(rho(i),ec0,vc0)
        ec = ec + ec0
        vc(i) = vc(i) + vc0
      END DO
    END IF

  END SUBROUTINE vwn_c_1

! *****************************************************************************

  SUBROUTINE vwn_c_2(rho,ec,vc,weight)

    REAL(KIND = dp), INTENT(INOUT)                        :: ec
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)           :: rho
    REAL(KIND = dp), DIMENSION(:,:), INTENT(INOUT)        :: vc
    REAL(KIND = dp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ec0,vc0
    INTEGER  :: i,j,n1,n2

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)

    IF (PRESENT(weight)) THEN
      DO i=1,n2
        DO j=1,n1
          CALL vwn_c_0(rho(j,i),ec0,vc0)
          ec = ec + weight(j,i)*ec0
          vc(j,i) = vc(j,i) + weight(j,i)*vc0
        END DO
      END DO
    ELSE
      DO i=1,n2
        DO j=1,n1
          CALL vwn_c_0(rho(j,i),ec0,vc0)
          ec = ec + ec0
          vc(j,i) = vc(j,i) + vc0
        END DO
      END DO
    END IF

  END SUBROUTINE vwn_c_2

! *****************************************************************************

  SUBROUTINE vwn_c_3(rho,ec,vc,weight)

    REAL(KIND = dp), INTENT(INOUT)                          :: ec
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(IN)           :: rho
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(INOUT)        :: vc
    REAL(KIND = dp), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ec0,vc0
    INTEGER  :: i,j,k,n1,n2,n3

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)
    n3 = SIZE(rho,3)

    IF (PRESENT(weight)) THEN
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL vwn_c_0(rho(k,j,i),ec0,vc0)
            ec = ec + weight(k,j,i)*ec0
            vc(k,j,i) = vc(k,j,i) + weight(k,j,i)*vc0
          END DO
        END DO
      END DO
    ELSE
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL vwn_c_0(rho(k,j,i),ec0,vc0)
            ec = ec + ec0
            vc(k,j,i) = vc(k,j,i) + vc0
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE vwn_c_3

! *****************************************************************************

  SUBROUTINE vwn_x_0(rho,ex,vx)

!   Purpose: Calculate the exchange energy and potential of Vosko, Wilk and
!            Nusair (VWN, 1980) for a closed shell density.

!   Literature: S. H. Vosko, L. Wilk and M. Nusair,
!               Can. J. Phys. 58, 1200 (1980)

!   History: - Creation (22.03.2000, Matthias Krack)

!   ***************************************************************************

    REAL(KIND = dp), INTENT(IN)  :: rho
    REAL(KIND = dp), INTENT(OUT) :: ex,vx

!   *** Local variables ***

    REAL(KIND = dp) :: rho13

!   ---------------------------------------------------------------------------

    IF (rho < eps_rho) THEN

      ex = 0.0_dp
      vx = 0.0_dp

    ELSE

      rho13 = rho**f13

      ex = cx_vwn_e*rho13*rho
      vx = cx_vwn_v*rho13

    END IF

  END SUBROUTINE vwn_x_0

! *****************************************************************************

  SUBROUTINE vwn_x_1(rho,ex,vx,weight)

    REAL(KIND = dp), INTENT(INOUT)                      :: ex
    REAL(KIND = dp), DIMENSION(:), INTENT(IN)           :: rho
    REAL(KIND = dp), DIMENSION(:), INTENT(INOUT)        :: vx
    REAL(KIND = dp), DIMENSION(:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0
    INTEGER  :: i,n1

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)

    IF (PRESENT(weight)) THEN
      DO i=1,n1
        CALL vwn_x_0(rho(i),ex0,vx0)
        ex = ex + weight(i)*ex0
        vx(i) = vx(i) + weight(i)*vx0
      END DO
    ELSE
      DO i=1,n1
        CALL vwn_x_0(rho(i),ex0,vx0)
        ex = ex + ex0
        vx(i) = vx(i) + vx0
      END DO
    END IF

  END SUBROUTINE vwn_x_1

! *****************************************************************************

  SUBROUTINE vwn_x_2(rho,ex,vx,weight)

    REAL(KIND = dp), INTENT(INOUT)                        :: ex
    REAL(KIND = dp), DIMENSION(:,:), INTENT(IN)           :: rho
    REAL(KIND = dp), DIMENSION(:,:), INTENT(INOUT)        :: vx
    REAL(KIND = dp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0
    INTEGER  :: i,j,n1,n2

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)

    IF (PRESENT(weight)) THEN
      DO i=1,n2
        DO j=1,n1
          CALL vwn_x_0(rho(j,i),ex0,vx0)
          ex = ex + weight(j,i)*ex0
          vx(j,i) = vx(j,i) + weight(j,i)*vx0
        END DO
      END DO
    ELSE
      DO i=1,n2
        DO j=1,n1
          CALL vwn_x_0(rho(j,i),ex0,vx0)
          ex = ex + ex0
          vx(j,i) = vx(j,i) + vx0
        END DO
      END DO
    END IF

  END SUBROUTINE vwn_x_2

! *****************************************************************************

  SUBROUTINE vwn_x_3(rho,ex,vx,weight)

    REAL(KIND = dp), INTENT(INOUT)                          :: ex
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(IN)           :: rho
    REAL(KIND = dp), DIMENSION(:,:,:), INTENT(INOUT)        :: vx
    REAL(KIND = dp), DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: weight

!   *** Local variables ***

    REAL(KIND = dp) :: ex0,vx0
    INTEGER  :: i,j,k,n1,n2,n3

!   ---------------------------------------------------------------------------

    n1 = SIZE(rho,1)
    n2 = SIZE(rho,2)
    n3 = SIZE(rho,3)

    IF (PRESENT(weight)) THEN
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL vwn_x_0(rho(k,j,i),ex0,vx0)
            ex = ex + weight(k,j,i)*ex0
            vx(k,j,i) = vx(k,j,i) + weight(k,j,i)*vx0
          END DO
        END DO
      END DO
    ELSE
      DO i=1,n3
        DO j=1,n2
          DO k=1,n1
            CALL vwn_x_0(rho(k,j,i),ex0,vx0)
            ex = ex + ex0
            vx(k,j,i) = vx(k,j,i) + vx0
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE vwn_x_3

! *****************************************************************************

END MODULE functionals

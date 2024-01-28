!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief calculates fxc in the spirit of the b97 exchange/correlation functional
!> \author jgh
! **************************************************************************************************
MODULE xc_b97_fxc
   USE kinds,                           ONLY: dp
   USE pw_types,                        ONLY: pw_r3d_rs_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: b97_fxc_eval, b97_fcc_eval

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param rhos ...
!> \param norm_drhos ...
!> \param fxc ...
!> \param gx ...
!> \param cx ...
!> \param eps_rho ...
! **************************************************************************************************
   SUBROUTINE b97_fxc_eval(rhos, norm_drhos, fxc, gx, cx, eps_rho)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rhos, norm_drhos
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: fxc
      REAL(KIND=dp), INTENT(IN)                          :: gx
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cx
      REAL(KIND=dp), INTENT(IN)                          :: eps_rho

      CHARACTER(len=*), PARAMETER                        :: routineN = 'b97_fxc_eval'

      INTEGER                                            :: handle, i, io, j, k, norder
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(KIND=dp)                                      :: drho, gval, rho, sval, ux

      CALL timeset(routineN, handle)

      norder = SIZE(cx)
      bo(1:2, 1:3) = rhos%pw_grid%bounds_local(1:2, 1:3)
!$OMP PARALLEL DO PRIVATE(i,j,k,io,rho,drho,sval,gval,ux) DEFAULT(NONE)&
!$OMP SHARED(bo,rhos,norm_drhos,fxc,gx,cx,eps_rho,norder)
      DO k = bo(1, 3), bo(2, 3)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 1), bo(2, 1)

               rho = rhos%array(i, j, k)
               drho = norm_drhos%array(i, j, k)
               IF (rho > eps_rho) THEN
                  sval = gx*(drho/rho**1.33333333333333_dp)**2
                  ux = sval/(1._dp + sval)
                  gval = 0.0_dp
                  DO io = 0, norder - 1
                     gval = gval + cx(io + 1)*(ux**io)
                  END DO
                  fxc%array(i, j, k) = fxc%array(i, j, k)*gval
               END IF

            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE b97_fxc_eval

! **************************************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param rhob ...
!> \param norm_drhoa ...
!> \param norm_drhob ...
!> \param fcc ...
!> \param gcc ...
!> \param cco ...
!> \param eps_rho ...
! **************************************************************************************************
   SUBROUTINE b97_fcc_eval(rhoa, rhob, norm_drhoa, norm_drhob, fcc, gcc, cco, eps_rho)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rhoa, rhob, norm_drhoa, norm_drhob
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: fcc
      REAL(KIND=dp), INTENT(IN)                          :: gcc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cco
      REAL(KIND=dp), INTENT(IN)                          :: eps_rho

      CHARACTER(len=*), PARAMETER                        :: routineN = 'b97_fcc_eval'

      INTEGER                                            :: handle, i, io, j, k, norder
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(KIND=dp)                                      :: dra, drb, gval, ra, rb, sa, sb, sval, ux

      CALL timeset(routineN, handle)

      norder = SIZE(cco)
      bo(1:2, 1:3) = rhoa%pw_grid%bounds_local(1:2, 1:3)
!$OMP PARALLEL DO PRIVATE(i,j,k,ra,rb,dra,drb,sa,sb,sval,gval,ux,io) DEFAULT(NONE)&
!$OMP SHARED(bo,rhoa,rhob,norm_drhoa,norm_drhob,fcc,gcc,cco,norder,eps_rho)
      DO k = bo(1, 3), bo(2, 3)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 1), bo(2, 1)

               ra = rhoa%array(i, j, k)
               rb = rhob%array(i, j, k)
               dra = norm_drhoa%array(i, j, k)
               drb = norm_drhob%array(i, j, k)
               IF (ra > eps_rho .AND. rb > eps_rho) THEN
                  sa = (dra/ra**1.33333333333333_dp)**2
                  sb = (drb/rb**1.33333333333333_dp)**2
                  sval = 0.5_dp*gcc*(sa + sb)
                  ux = sval/(1._dp + sval)
                  gval = 0.0_dp
                  DO io = 0, norder - 1
                     gval = gval + cco(io + 1)*(ux**io)
                  END DO
                  fcc%array(i, j, k) = fcc%array(i, j, k)*gval
               END IF

            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE b97_fcc_eval

END MODULE xc_b97_fxc

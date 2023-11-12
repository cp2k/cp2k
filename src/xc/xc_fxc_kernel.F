!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Exchange and Correlation kernel functionals
!> \author JGH
! **************************************************************************************************
MODULE xc_fxc_kernel
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_scale,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
   USE xc_b97_fxc,                      ONLY: b97_fcc_eval,&
                                              b97_fxc_eval
   USE xc_input_constants,              ONLY: xc_deriv_pw
   USE xc_pade,                         ONLY: pade_fxc_eval,&
                                              pade_init
   USE xc_perdew_wang,                  ONLY: perdew_wang_fxc_calc
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_setall,&
                                              xc_rho_cflags_type
   USE xc_util,                         ONLY: xc_pw_gradient
   USE xc_xalpha,                       ONLY: xalpha_fxc_eval
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: calc_fxc_kernel

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_fxc_kernel'

CONTAINS

! **************************************************************************************************
!> \brief Exchange and Correlation kernel functional calculations
!> \param fxc_rspace ...
!> \param rho_r the value of the density in the real space
!> \param rho_g value of the density in the g space (needs to be associated
!>        only for gradient corrections)
!> \param tau_r value of the kinetic density tau on the grid (can be null,
!>        used only with meta functionals)
!> \param xc_kernel which functional to calculate, and how to do it
!> \param triplet ...
!> \param pw_pool the pool for the grids
!> \author JGH
! **************************************************************************************************
   SUBROUTINE calc_fxc_kernel(fxc_rspace, rho_r, rho_g, tau_r, xc_kernel, triplet, pw_pool)
      TYPE(pw_type), DIMENSION(:)                        :: fxc_rspace, rho_r, rho_g, tau_r
      TYPE(section_vals_type), POINTER                   :: xc_kernel
      LOGICAL, INTENT(IN)                                :: triplet
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      CHARACTER(len=*), PARAMETER                        :: routineN = 'calc_fxc_kernel'
      REAL(KIND=dp), PARAMETER                           :: eps_rho = 1.E-10_dp

      CHARACTER(len=20)                                  :: fxc_name
      INTEGER                                            :: handle, i, idir, j, k, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: lsd
      REAL(KIND=dp)                                      :: scalec, scalex
      REAL(KIND=dp), DIMENSION(3)                        :: ccaa, ccab, cxaa, g_ab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rvals
      TYPE(pw_type)                                      :: fxa, fxb, norm_drhoa, norm_drhob, rhoa, &
                                                            rhob, rhog, tmpg
      TYPE(pw_type), DIMENSION(3)                        :: drhoa
      TYPE(xc_rho_cflags_type)                           :: needs

      CPASSERT(ASSOCIATED(xc_kernel))
      CPASSERT(ASSOCIATED(pw_pool))

      CALL timeset(routineN, handle)

      nspins = SIZE(rho_r)
      lsd = (nspins == 2)
      IF (triplet) THEN
         CPASSERT(nspins == 1)
      END IF

      CALL section_vals_val_get(xc_kernel, "_SECTION_PARAMETERS_", c_val=fxc_name)
      CALL section_vals_val_get(xc_kernel, "SCALE_X", r_val=scalex)
      CALL section_vals_val_get(xc_kernel, "SCALE_C", r_val=scalec)

      CALL xc_rho_cflags_setall(needs, .FALSE.)
      CALL fxc_kernel_info(fxc_name, needs, lsd)

      CALL pw_pool%create_pw(rhoa, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool%create_pw(rhob, use_data=REALDATA3D, in_space=REALSPACE)
      IF (lsd) THEN
         CALL pw_copy(rho_r(1), rhoa)
         CALL pw_copy(rho_r(2), rhob)
      ELSE IF (triplet) THEN
         CALL pw_copy(rho_r(1), rhoa)
         CALL pw_copy(rho_r(1), rhob)
      ELSE
         CALL pw_copy(rho_r(1), rhoa)
         CALL pw_copy(rho_r(1), rhob)
         CALL pw_scale(rhoa, 0.5_dp)
         CALL pw_scale(rhob, 0.5_dp)
      END IF
      IF (needs%norm_drho) THEN
         ! deriv rho
         DO idir = 1, 3
            CALL pw_pool%create_pw(drhoa(idir), use_data=REALDATA3D, in_space=REALSPACE)
         END DO
         CALL pw_pool%create_pw(norm_drhoa, use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool%create_pw(norm_drhob, use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool%create_pw(rhog, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         CALL pw_pool%create_pw(tmpg, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         IF (lsd) THEN
            CALL pw_copy(rho_g(1), rhog)
         ELSE IF (triplet) THEN
            CALL pw_copy(rho_g(1), rhog)
         ELSE
            CALL pw_copy(rho_g(1), rhog)
            CALL pw_scale(rhog, 0.5_dp)
         END IF
         CALL xc_pw_gradient(rhoa, rhog, tmpg, drhoa(:), xc_deriv_pw)
         bo(1:2, 1:3) = rhoa%pw_grid%bounds_local(1:2, 1:3)
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(bo,norm_drhoa,drhoa)
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  norm_drhoa%cr3d(i, j, k) = SQRT(drhoa(1)%cr3d(i, j, k)**2 + &
                                                  drhoa(2)%cr3d(i, j, k)**2 + &
                                                  drhoa(3)%cr3d(i, j, k)**2)
               END DO
            END DO
         END DO
         IF (lsd) THEN
            CALL pw_copy(rho_g(2), rhog)
            CALL xc_pw_gradient(rhob, rhog, tmpg, drhoa(:), xc_deriv_pw)
            bo(1:2, 1:3) = rhob%pw_grid%bounds_local(1:2, 1:3)
!$OMP       PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(bo,norm_drhob,drhoa)
            DO k = bo(1, 3), bo(2, 3)
               DO j = bo(1, 2), bo(2, 2)
                  DO i = bo(1, 1), bo(2, 1)
                     norm_drhob%cr3d(i, j, k) = SQRT(drhoa(1)%cr3d(i, j, k)**2 + &
                                                     drhoa(2)%cr3d(i, j, k)**2 + &
                                                     drhoa(3)%cr3d(i, j, k)**2)
                  END DO
               END DO
            END DO
         ELSE
            norm_drhob%cr3d(:, :, :) = norm_drhoa%cr3d(:, :, :)
         END IF
         CALL pw_pool%give_back_pw(rhog)
         CALL pw_pool%give_back_pw(tmpg)
      END IF
      IF (needs%tau) THEN
         MARK_USED(tau_r)
         CPABORT("Meta functionals not available.")
      END IF

      SELECT CASE (TRIM(fxc_name))
      CASE ("PADEFXC")
         IF (scalec == scalex) THEN
            CALL pade_init(eps_rho)
            CALL pade_fxc_eval(rhoa, rhob, fxc_rspace(1), fxc_rspace(2), fxc_rspace(3))
            IF (scalex /= 1.0_dp) THEN
               CALL pw_scale(fxc_rspace(1), scalex)
               CALL pw_scale(fxc_rspace(2), scalex)
               CALL pw_scale(fxc_rspace(3), scalex)
            END IF
         ELSE
            CPABORT("PADE Fxc Kernel functional needs SCALE_X==SCALE_C")
         END IF
      CASE ("LDAFXC")
         CALL pw_zero(fxc_rspace(1))
         CALL pw_zero(fxc_rspace(2))
         CALL pw_zero(fxc_rspace(3))
         CALL xalpha_fxc_eval(rhoa, rhob, fxc_rspace(1), fxc_rspace(3), scalex, eps_rho)
         CALL perdew_wang_fxc_calc(rhoa, rhob, fxc_rspace(1), fxc_rspace(2), fxc_rspace(3), &
                                   scalec, eps_rho)
      CASE ("GGAFXC")
         ! get parameter
         CALL section_vals_val_get(xc_kernel, "GAMMA", r_vals=rvals)
         g_ab(1:3) = rvals(1:3)
         CALL section_vals_val_get(xc_kernel, "C_XAA", r_vals=rvals)
         cxaa(1:3) = rvals(1:3)
         CALL section_vals_val_get(xc_kernel, "C_CAA", r_vals=rvals)
         ccaa(1:3) = rvals(1:3)
         CALL section_vals_val_get(xc_kernel, "C_CAB", r_vals=rvals)
         ccab(1:3) = rvals(1:3)
         ! correlation
         CALL pw_zero(fxc_rspace(1))
         CALL pw_zero(fxc_rspace(2))
         CALL pw_zero(fxc_rspace(3))
         CALL perdew_wang_fxc_calc(rhoa, rhob, fxc_rspace(1), fxc_rspace(2), fxc_rspace(3), &
                                   scalec, eps_rho)
         CALL b97_fxc_eval(rhoa, norm_drhoa, fxc_rspace(1), g_ab(1), ccaa, eps_rho)
         CALL b97_fxc_eval(rhob, norm_drhob, fxc_rspace(3), g_ab(3), ccaa, eps_rho)
         CALL b97_fcc_eval(rhoa, rhob, norm_drhoa, norm_drhob, fxc_rspace(2), g_ab(2), ccab, eps_rho)
         ! exchange
         CALL pw_pool%create_pw(fxa, use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool%create_pw(fxb, use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_zero(fxa)
         CALL pw_zero(fxb)
         CALL xalpha_fxc_eval(rhoa, rhob, fxa, fxb, scalex, eps_rho)
         CALL b97_fxc_eval(rhoa, norm_drhoa, fxa, g_ab(1), cxaa, eps_rho)
         CALL b97_fxc_eval(rhob, norm_drhob, fxb, g_ab(1), cxaa, eps_rho)
         CALL pw_axpy(fxa, fxc_rspace(1))
         CALL pw_axpy(fxb, fxc_rspace(3))
         CALL pw_pool%give_back_pw(fxa)
         CALL pw_pool%give_back_pw(fxb)
      CASE ("NONE")
         CALL pw_zero(fxc_rspace(1))
         CALL pw_zero(fxc_rspace(2))
         CALL pw_zero(fxc_rspace(3))
      CASE default
         CPABORT("Fxc Kernel functional is defined incorrectly")
      END SELECT

      CALL pw_pool%give_back_pw(rhoa)
      CALL pw_pool%give_back_pw(rhob)
      IF (needs%norm_drho) THEN
         CALL pw_pool%give_back_pw(norm_drhoa)
         CALL pw_pool%give_back_pw(norm_drhob)
         DO idir = 1, 3
            CALL pw_pool%give_back_pw(drhoa(idir))
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE calc_fxc_kernel

! **************************************************************************************************
!> \brief ...
!> \param fxc_name ...
!> \param needs ...
!> \param lsd ...
! **************************************************************************************************
   SUBROUTINE fxc_kernel_info(fxc_name, needs, lsd)
      CHARACTER(len=20), INTENT(IN)                      :: fxc_name
      TYPE(xc_rho_cflags_type), INTENT(INOUT)            :: needs
      LOGICAL, INTENT(IN)                                :: lsd

      SELECT CASE (TRIM(fxc_name))
      CASE ("PADEFXC", "LDAFXC")
         IF (lsd) THEN
            needs%rho_spin = .TRUE.
         ELSE
            needs%rho = .TRUE.
         END IF
      CASE ("GGAFXC")
         IF (lsd) THEN
            needs%rho_spin = .TRUE.
            needs%norm_drho_spin = .TRUE.
            needs%norm_drho = .TRUE.
         ELSE
            needs%rho = .TRUE.
            needs%norm_drho = .TRUE.
         END IF
      CASE ("NONE")
      CASE default
         CPABORT("Fxc Kernel functional is defined incorrectly")
      END SELECT

   END SUBROUTINE fxc_kernel_info

END MODULE xc_fxc_kernel


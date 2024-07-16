!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_kernel_types
   USE admm_types,                      ONLY: admm_env_release,&
                                              admm_type
   USE cp_control_types,                ONLY: admm_control_release,&
                                              admm_control_type
   USE hfx_types,                       ONLY: hfx_release,&
                                              hfx_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: lri_density_type,&
                                              lri_environment_type
   USE qs_tddfpt2_stda_types,           ONLY: stda_env_type
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_release
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_release,&
                                              xc_rho_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kernel_types'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   ! number of first derivative components (3: d/dx, d/dy, d/dz)
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC :: full_kernel_env_type, kernel_env_type
   PUBLIC :: release_kernel_env

! **************************************************************************************************
!> \brief Collection of variables required to evaluate adiabatic TDDFPT kernel.
!> \par History
!>   * 12.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE full_kernel_env_type
      ! ground state electron density
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set => Null()
      ! response density
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho1_set => Null()
      !> first and second derivatives of exchange-correlation functional
      TYPE(xc_derivative_set_type)                       :: xc_deriv_set = xc_derivative_set_type()
      !> XC input section
      LOGICAL                                            :: do_exck = .FALSE.
      TYPE(section_vals_type), POINTER                   :: xc_section => Null()
      !> flags which indicate required components of the exchange-correlation functional
      !> (density, gradient, etc)
      TYPE(xc_rho_cflags_type)                           :: xc_rho1_cflags = xc_rho_cflags_type()
      !> the method used to compute position derivatives of densities and potentials
      INTEGER                                            :: deriv_method_id = -1
      !> the density smoothing method
      INTEGER                                            :: rho_smooth_id = -1
      !> scaling coefficients in the linear combination:
      !> K = alpha * K_{\alpha,\alpha} + beta * K_{\alpha,\beta}
      REAL(kind=dp)                                      :: alpha = -1.0_dp, beta = -1.0_dp
      !> flags for finite differences/analytic XC kernels
      LOGICAL                                            :: deriv2_analytic = .FALSE.
      LOGICAL                                            :: deriv3_analytic = .FALSE.
      ! Local resolution of the identity for Coulomb
      TYPE(lri_environment_type), POINTER                :: lri_env => Null()
      TYPE(lri_density_type), POINTER                    :: lri_density => Null()
      ! Short range HFX integral environment
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data => Null()
      TYPE(section_vals_type), POINTER                   :: hfxsr_section => Null()
      TYPE(admm_type), POINTER                           :: admm_env => Null()
      TYPE(admm_control_type), POINTER                   :: admm_control => NULL()
   END TYPE full_kernel_env_type
! **************************************************************************************************
!> \brief Type to hold environments for the different kernels
!> \par History
!>   * 04.2019 created [JHU]
! **************************************************************************************************
   TYPE kernel_env_type
      TYPE(full_kernel_env_type), POINTER                :: full_kernel => Null()
      TYPE(full_kernel_env_type), POINTER                :: admm_kernel => Null()
      TYPE(stda_env_type), POINTER                       :: stda_kernel => Null()
   END TYPE kernel_env_type

CONTAINS

! **************************************************************************************************
!> \brief Release kernel environment.
!> \param kernel_env  kernel environment (destroyed on exit)
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE release_kernel_env(kernel_env)
      TYPE(full_kernel_env_type), POINTER                :: kernel_env

      IF (ASSOCIATED(kernel_env)) THEN
         IF (ASSOCIATED(kernel_env%xc_rho1_set)) THEN
            CALL xc_rho_set_release(kernel_env%xc_rho1_set)
            DEALLOCATE (kernel_env%xc_rho1_set)
         END IF
         CALL xc_dset_release(kernel_env%xc_deriv_set)
         IF (ASSOCIATED(kernel_env%xc_rho_set)) THEN
            CALL xc_rho_set_release(kernel_env%xc_rho_set)
            DEALLOCATE (kernel_env%xc_rho_set)
         END IF
         IF (ASSOCIATED(kernel_env%x_data)) THEN
            CALL hfx_release(kernel_env%x_data)
         END IF
         IF (ASSOCIATED(kernel_env%admm_env)) THEN
            CALL admm_env_release(kernel_env%admm_env)
            DEALLOCATE (kernel_env%admm_env)
         END IF
         IF (ASSOCIATED(kernel_env%admm_control)) THEN
            CALL admm_control_release(kernel_env%admm_control)
         END IF
      END IF

   END SUBROUTINE release_kernel_env

END MODULE qs_kernel_types

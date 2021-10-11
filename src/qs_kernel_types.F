!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_kernel_types
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
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
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set
      ! response density
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho1_set
      !> first and second derivatives of exchange-correlation functional
      TYPE(xc_derivative_set_type), POINTER              :: xc_deriv_set
      !> XC input section
      TYPE(section_vals_type), POINTER                   :: xc_section
      !> flags which indicate required components of the exchange-correlation functional
      !> (density, gradient, etc)
      TYPE(xc_rho_cflags_type)                           :: xc_rho1_cflags
      !> the method used to compute position derivatives of densities and potentials
      INTEGER                                            :: deriv_method_id
      !> the density smoothing method
      INTEGER                                            :: rho_smooth_id
      !> scaling coefficients in the linear combination:
      !> K = alpha * K_{\alpha,\alpha} + beta * K_{\alpha,\beta}
      REAL(kind=dp)                                      :: alpha, beta
      !> flags for finite differences/analytic XC kernels
      LOGICAL                                            :: deriv2_analytic
      LOGICAL                                            :: deriv3_analytic
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'release_kernel_env', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(kernel_env)) THEN
         CALL xc_rho_set_release(kernel_env%xc_rho1_set)
         CALL xc_dset_release(kernel_env%xc_deriv_set)
         CALL xc_rho_set_release(kernel_env%xc_rho_set)
      END IF

   END SUBROUTINE release_kernel_env

END MODULE qs_kernel_types

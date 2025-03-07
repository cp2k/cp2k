!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Thermal regions type: to initialize and control the temperature of
!>        different regions
!> \par History
!>   - Added support for langevin regions (2014/01/08, LT)
!> \author MI
! **************************************************************************************************
MODULE thermal_region_types

   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: thermal_regions_type, &
             thermal_region_type, &
             allocate_thermal_regions, &
             release_thermal_regions

   TYPE thermal_regions_type
      INTEGER :: nregions = 0
      LOGICAL :: force_rescaling = .FALSE.
      REAL(KIND=dp) :: temp_reg0 = 0.0_dp
      LOGICAL, DIMENSION(:), POINTER                   :: do_langevin => NULL()
      TYPE(section_vals_type), POINTER                 :: section => NULL()
      TYPE(thermal_region_type), DIMENSION(:), POINTER :: thermal_region => NULL()
   END TYPE thermal_regions_type

   TYPE thermal_region_type
      INTEGER :: region_index = 0, npart = 0
      INTEGER, DIMENSION(:), POINTER :: part_index => NULL()
      REAL(KIND=dp) :: ekin = 0.0_dp, noisy_gamma_region = 0.0_dp, temperature = 0.0_dp, temp_expected = 0.0_dp, temp_tol = 0.0_dp
   END TYPE thermal_region_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermal_region_types'
CONTAINS

! **************************************************************************************************
!> \brief allocate thermal_regions
!> \param thermal_regions ...
!> \author
! **************************************************************************************************
   SUBROUTINE allocate_thermal_regions(thermal_regions)
      TYPE(thermal_regions_type), INTENT(OUT)            :: thermal_regions

      MARK_USED(thermal_regions)
   END SUBROUTINE allocate_thermal_regions

! **************************************************************************************************
!> \brief release thermal_regions
!> \param thermal_regions ...
!> \author
! **************************************************************************************************
   SUBROUTINE release_thermal_regions(thermal_regions)

      TYPE(thermal_regions_type), INTENT(INOUT)          :: thermal_regions

      INTEGER                                            :: ireg

      IF (ASSOCIATED(thermal_regions%thermal_region)) THEN
         DO ireg = 1, SIZE(thermal_regions%thermal_region)
            DEALLOCATE (thermal_regions%thermal_region(ireg)%part_index)
         END DO
         DEALLOCATE (thermal_regions%thermal_region)
      END IF
      IF (ASSOCIATED(thermal_regions%do_langevin)) THEN
         DEALLOCATE (thermal_regions%do_langevin)
      END IF

   END SUBROUTINE release_thermal_regions

END MODULE thermal_region_types

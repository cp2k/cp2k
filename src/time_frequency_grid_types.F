!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shared time/frequency quadrature data used by RPA and GW paths.
MODULE time_frequency_grid_types

   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'time_frequency_grid_types'

   TYPE time_frequency_grid_type
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: frequency
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: frequency_weights
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: imaginary_time
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: time_weights_at_zero_frequency
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: cosine_time_to_frequency_weights
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: cosine_frequency_to_time_weights
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: sine_time_to_frequency_weights
   END TYPE time_frequency_grid_type

   PUBLIC :: time_frequency_grid_type, time_frequency_grid_release

CONTAINS

! **************************************************************************************************
!> \brief Release all data owned by a time_frequency_grid_type object.
!> \param grid The grid object to release.
! **************************************************************************************************
   SUBROUTINE time_frequency_grid_release(grid)

      TYPE(time_frequency_grid_type), INTENT(INOUT)      :: grid

      IF (ALLOCATED(grid%frequency)) DEALLOCATE (grid%frequency)
      IF (ALLOCATED(grid%frequency_weights)) DEALLOCATE (grid%frequency_weights)
      IF (ALLOCATED(grid%imaginary_time)) DEALLOCATE (grid%imaginary_time)
      IF (ALLOCATED(grid%time_weights_at_zero_frequency)) DEALLOCATE (grid%time_weights_at_zero_frequency)
      IF (ALLOCATED(grid%cosine_time_to_frequency_weights)) DEALLOCATE (grid%cosine_time_to_frequency_weights)
      IF (ALLOCATED(grid%cosine_frequency_to_time_weights)) DEALLOCATE (grid%cosine_frequency_to_time_weights)
      IF (ALLOCATED(grid%sine_time_to_frequency_weights)) DEALLOCATE (grid%sine_time_to_frequency_weights)

   END SUBROUTINE time_frequency_grid_release

END MODULE time_frequency_grid_types

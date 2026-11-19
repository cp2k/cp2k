!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Canonical names and order of the SKALA TorchScript feature dictionary.
!>
!> All code that builds a SKALA input dictionary and the code that reads the
!> model's exported feature metadata must agree on these keys; keep them here
!> instead of repeating string literals.
! **************************************************************************************************
MODULE skala_feature_keys
   USE kinds,                           ONLY: default_string_length
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_feature_keys'

   INTEGER, PARAMETER :: skala_n_features = 9

   CHARACTER(len=*), PARAMETER :: &
      skala_key_grid_coords = "grid_coords", &
      skala_key_grid_weights = "grid_weights", &
      skala_key_atomic_grid_weights = "atomic_grid_weights", &
      skala_key_atomic_grid_sizes = "atomic_grid_sizes", &
      skala_key_atomic_grid_size_bound_shape = "atomic_grid_size_bound_shape", &
      skala_key_density = "density", &
      skala_key_grad = "grad", &
      skala_key_kin = "kin", &
      skala_key_coarse_0_atomic_coords = "coarse_0_atomic_coords"

   ! Canonical order of the feature dictionary, matching the model's exported
   ! metadata order ("features" extra file of the TorchScript archive).
   CHARACTER(len=default_string_length), PARAMETER :: skala_feature_names(skala_n_features) = [ &
                                                      CHARACTER(len=default_string_length) :: &
                                                      skala_key_density, skala_key_kin, skala_key_grad, skala_key_grid_coords, &
                                               skala_key_grid_weights, skala_key_atomic_grid_weights, skala_key_atomic_grid_sizes, &
                                                      skala_key_coarse_0_atomic_coords, skala_key_atomic_grid_size_bound_shape]

   PUBLIC :: skala_key_grid_coords, skala_key_grid_weights, skala_key_atomic_grid_weights
   PUBLIC :: skala_key_atomic_grid_sizes, skala_key_atomic_grid_size_bound_shape
   PUBLIC :: skala_key_density, skala_key_grad, skala_key_kin, skala_key_coarse_0_atomic_coords
   PUBLIC :: skala_validate_feature_list

CONTAINS

! **************************************************************************************************
!> \brief Abort if the model's exported feature list does not match the protocol.
!> \param features parsed "features" metadata of the TorchScript model
! **************************************************************************************************
   SUBROUTINE skala_validate_feature_list(features)
      CHARACTER(len=*), DIMENSION(:), INTENT(IN)         :: features

      INTEGER                                            :: i

      IF (SIZE(features) /= skala_n_features) THEN
         CALL cp_abort(__LOCATION__, moduleN//': SKALA TorchScript model exports an '// &
                       'unexpected number of features; CP2K expects the SKALA feature protocol.')
      END IF
      DO i = 1, skala_n_features
         IF (TRIM(features(i)) /= TRIM(skala_feature_names(i))) THEN
            CALL cp_abort(__LOCATION__, moduleN//': SKALA TorchScript model feature mismatch: '// &
                          'expected "'//TRIM(skala_feature_names(i))//'", got "'// &
                          TRIM(features(i))//'".')
         END IF
      END DO

   END SUBROUTINE skala_validate_feature_list

END MODULE skala_feature_keys

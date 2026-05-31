!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Small CP2K wrapper around the SKALA TorchScript functional protocol.
! **************************************************************************************************
MODULE skala_torch_api
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE string_utilities,                ONLY: uppercase
   USE torch_api,                       ONLY: &
        torch_dict_type, torch_model_forward_mol_tensor, torch_model_load, &
        torch_model_read_metadata, torch_model_release, torch_model_type, &
        torch_tensor_item_double, torch_tensor_release, torch_tensor_type, &
        torch_tensor_weighted_sum
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_torch_api'

   PUBLIC :: skala_torch_model_type, skala_torch_model_load, skala_torch_model_release
   PUBLIC :: skala_torch_model_get_exc, skala_torch_model_get_exc_density
   PUBLIC :: skala_torch_model_needs_feature, skala_torch_model_protocol_version

   TYPE skala_torch_model_type
      PRIVATE
      INTEGER                                            :: protocol_version = -1
      CHARACTER(len=default_string_length), ALLOCATABLE, &
         DIMENSION(:)                                    :: features
      TYPE(torch_model_type)                             :: torch_model
   END TYPE skala_torch_model_type

CONTAINS

! **************************************************************************************************
!> \brief Load a SKALA TorchScript model and its feature metadata.
!> \param model ...
!> \param filename ...
! **************************************************************************************************
   SUBROUTINE skala_torch_model_load(model, filename)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model
      CHARACTER(len=*), INTENT(IN)                       :: filename

      CHARACTER(:), ALLOCATABLE                          :: features_json, protocol_string
      INTEGER                                            :: ios

      CALL torch_model_load(model%torch_model, filename)
      protocol_string = torch_model_read_metadata(filename, "protocol_version")
      features_json = torch_model_read_metadata(filename, "features")
      READ (protocol_string, *, IOSTAT=ios) model%protocol_version
      IF (ios /= 0) CPABORT("Could not parse SKALA TorchScript protocol_version metadata")
      IF (model%protocol_version /= 2) THEN
         CPABORT("Unsupported SKALA TorchScript protocol version")
      END IF

      CALL parse_feature_list(features_json, model%features)

   END SUBROUTINE skala_torch_model_load

! **************************************************************************************************
!> \brief Release a loaded SKALA TorchScript model.
!> \param model ...
! **************************************************************************************************
   SUBROUTINE skala_torch_model_release(model)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model

      CALL torch_model_release(model%torch_model)
      IF (ALLOCATED(model%features)) DEALLOCATE (model%features)
      model%protocol_version = -1

   END SUBROUTINE skala_torch_model_release

! **************************************************************************************************
!> \brief Check whether a loaded SKALA model requests a feature.
!> \param model ...
!> \param feature ...
!> \return ...
! **************************************************************************************************
   FUNCTION skala_torch_model_needs_feature(model, feature) RESULT(needs_feature)
      TYPE(skala_torch_model_type), INTENT(IN)           :: model
      CHARACTER(len=*), INTENT(IN)                       :: feature
      LOGICAL                                            :: needs_feature

      CHARACTER(len=default_string_length)               :: feature_key, model_feature
      INTEGER                                            :: i

      feature_key = ADJUSTL(feature)
      CALL uppercase(feature_key)

      needs_feature = .FALSE.
      IF (.NOT. ALLOCATED(model%features)) RETURN

      DO i = 1, SIZE(model%features)
         model_feature = ADJUSTL(model%features(i))
         CALL uppercase(model_feature)
         IF (TRIM(model_feature) == TRIM(feature_key)) THEN
            needs_feature = .TRUE.
            RETURN
         END IF
      END DO

   END FUNCTION skala_torch_model_needs_feature

! **************************************************************************************************
!> \brief Return the loaded SKALA TorchScript protocol version.
!> \param model ...
!> \return ...
! **************************************************************************************************
   FUNCTION skala_torch_model_protocol_version(model) RESULT(protocol_version)
      TYPE(skala_torch_model_type), INTENT(IN)           :: model
      INTEGER                                            :: protocol_version

      protocol_version = model%protocol_version

   END FUNCTION skala_torch_model_protocol_version

! **************************************************************************************************
!> \brief Evaluate the SKALA exchange-correlation energy density.
!> \param model ...
!> \param inputs ...
!> \param exc_density ...
! **************************************************************************************************
   SUBROUTINE skala_torch_model_get_exc_density(model, inputs, exc_density)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model
      TYPE(torch_dict_type), INTENT(IN)                  :: inputs
      TYPE(torch_tensor_type), INTENT(INOUT)             :: exc_density

      CALL torch_model_forward_mol_tensor(model%torch_model, "get_exc_density", inputs, exc_density)

   END SUBROUTINE skala_torch_model_get_exc_density

! **************************************************************************************************
!> \brief Evaluate the weighted SKALA exchange-correlation energy.
!> \param model ...
!> \param inputs ...
!> \param grid_weights ...
!> \param exc_tensor ...
!> \param exc ...
! **************************************************************************************************
   SUBROUTINE skala_torch_model_get_exc(model, inputs, grid_weights, exc_tensor, exc)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model
      TYPE(torch_dict_type), INTENT(IN)                  :: inputs
      TYPE(torch_tensor_type), INTENT(IN)                :: grid_weights
      TYPE(torch_tensor_type), INTENT(INOUT)             :: exc_tensor
      REAL(KIND=dp), INTENT(OUT)                         :: exc

      TYPE(torch_tensor_type)                            :: exc_density

      CALL skala_torch_model_get_exc_density(model, inputs, exc_density)
      CALL torch_tensor_weighted_sum(exc_density, grid_weights, exc_tensor)
      exc = torch_tensor_item_double(exc_tensor)
      CALL torch_tensor_release(exc_density)

   END SUBROUTINE skala_torch_model_get_exc

! **************************************************************************************************
!> \brief Parse a TorchScript extra_files JSON list of feature names.
!> \param features_json ...
!> \param features ...
! **************************************************************************************************
   SUBROUTINE parse_feature_list(features_json, features)
      CHARACTER(len=*), INTENT(IN)                       :: features_json
      CHARACTER(len=default_string_length), &
         ALLOCATABLE, DIMENSION(:), INTENT(OUT)          :: features

      INTEGER                                            :: end_pos, feature_count, i, pos, quote1, &
                                                            quote2, start_pos

      feature_count = 0
      pos = 1
      DO
         quote1 = INDEX(features_json(pos:), '"')
         IF (quote1 == 0) EXIT
         start_pos = pos + quote1
         quote2 = INDEX(features_json(start_pos:), '"')
         IF (quote2 == 0) EXIT
         feature_count = feature_count + 1
         pos = start_pos + quote2
      END DO

      IF (feature_count == 0) CPABORT("SKALA TorchScript model does not list any features")
      ALLOCATE (features(feature_count))
      features = ""

      pos = 1
      DO i = 1, feature_count
         quote1 = INDEX(features_json(pos:), '"')
         start_pos = pos + quote1
         quote2 = INDEX(features_json(start_pos:), '"')
         end_pos = start_pos + quote2 - 2
         features(i) = features_json(start_pos:end_pos)
         pos = start_pos + quote2
      END DO

   END SUBROUTINE parse_feature_list

END MODULE skala_torch_api

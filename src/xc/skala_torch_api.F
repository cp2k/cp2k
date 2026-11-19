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
#if defined (__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: ieee_all, &
                              ieee_get_halting_mode, &
                              ieee_set_halting_mode
#endif
   USE kinds, ONLY: default_string_length, &
                    dp, &
                    int_8
   USE skala_feature_keys, ONLY: skala_validate_feature_list
   USE torch_api, ONLY: &
      torch_dict_type, torch_model_disable_parameter_gradients, torch_model_forward_mol_tensor, &
      torch_model_freeze_preserving_method, torch_model_load_with_metadata, torch_model_release, &
      torch_model_remap_device_constants, torch_model_promote_float32, torch_model_type, &
      torch_tensor_item_double, torch_tensor_release, torch_tensor_type, &
      torch_tensor_weighted_sum, torch_tensor_weighted_sum_segments
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_torch_api'

   PUBLIC :: skala_torch_model_type, skala_torch_model_load, skala_torch_model_release
   PUBLIC :: skala_torch_model_get_exc, skala_torch_model_get_exc_device
   PUBLIC :: skala_torch_model_get_exc_segments

   TYPE skala_torch_model_type
      PRIVATE
      TYPE(torch_model_type)                             :: torch_model
   END TYPE skala_torch_model_type

CONTAINS

! **************************************************************************************************
!> \brief Load a SKALA TorchScript model and its feature metadata.
!> \param model ...
!> \param filename ...
!> \param fp64 Promote float32 numerical state and casts; default true for differential consistency.
! **************************************************************************************************
   SUBROUTINE skala_torch_model_load(model, filename, fp64)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model
      CHARACTER(len=*), INTENT(IN)                       :: filename
      LOGICAL, INTENT(IN), OPTIONAL                      :: fp64

      CHARACTER(:), ALLOCATABLE                          :: features_json, protocol_string
      INTEGER                                            :: ios, protocol_version
      LOGICAL                                            :: use_fp64

      ! The model is loaded on the CPU (even for a CUDA target) so that all
      ! float32 -> float64 promotion and constant folding happen in host memory.
      CALL torch_model_load_with_metadata(model%torch_model, filename, &
                                          "protocol_version", protocol_string, &
                                          "features", features_json)
      CALL torch_model_disable_parameter_gradients(model%torch_model)
      READ (protocol_string, *, IOSTAT=ios) protocol_version
      IF (ios /= 0) CPABORT("Could not parse SKALA TorchScript protocol_version metadata")
      IF (protocol_version /= 2) THEN
         CPABORT("Unsupported SKALA TorchScript protocol version")
      END IF

      BLOCK
         CHARACTER(len=default_string_length), ALLOCATABLE :: features(:)
         CALL parse_feature_list(features_json, features)
         CALL skala_validate_feature_list(features)
         DEALLOCATE (features)
      END BLOCK
      use_fp64 = .TRUE.
      IF (PRESENT(fp64)) use_fp64 = fp64
      IF (use_fp64) CALL torch_model_promote_float32(model%torch_model)
      ! Move to the resolved device (a no-op on the CPU target) and remap the
      ! serialized device constants after promotion.
      CALL torch_model_remap_device_constants(model%torch_model)
      ! Preserve the exported SKALA entry point while folding constant model state.
      CALL torch_model_freeze_preserving_method(model%torch_model, "get_exc_density")

   END SUBROUTINE skala_torch_model_load

! **************************************************************************************************
!> \brief Release a loaded SKALA TorchScript model.
!> \param model ...
! **************************************************************************************************
   SUBROUTINE skala_torch_model_release(model)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model

      CALL torch_model_release(model%torch_model)

   END SUBROUTINE skala_torch_model_release

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

      CALL skala_torch_model_get_exc_impl(model, inputs, grid_weights, exc_tensor, .TRUE., exc)

   END SUBROUTINE skala_torch_model_get_exc

! **************************************************************************************************
!> \brief Evaluate the weighted SKALA energy on the device without reading it back.
!> \param model ...
!> \param inputs ...
!> \param grid_weights ...
!> \param exc_tensor ...
! **************************************************************************************************
   SUBROUTINE skala_torch_model_get_exc_device(model, inputs, grid_weights, exc_tensor)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model
      TYPE(torch_dict_type), INTENT(IN)                  :: inputs
      TYPE(torch_tensor_type), INTENT(IN)                :: grid_weights
      TYPE(torch_tensor_type), INTENT(INOUT)             :: exc_tensor

      REAL(KIND=dp)                                      :: exc

      CALL skala_torch_model_get_exc_impl(model, inputs, grid_weights, exc_tensor, .FALSE., exc)

   END SUBROUTINE skala_torch_model_get_exc_device

! **************************************************************************************************
!> \brief Shared implementation of the weighted SKALA energy; optionally read it back.
!> \param model ...
!> \param inputs ...
!> \param grid_weights ...
!> \param exc_tensor ...
!> \param read_exc whether to read the scalar energy back to the host
!> \param exc scalar energy (0 when read_exc is false)
! **************************************************************************************************
   SUBROUTINE skala_torch_model_get_exc_impl(model, inputs, grid_weights, exc_tensor, read_exc, exc)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model
      TYPE(torch_dict_type), INTENT(IN)                  :: inputs
      TYPE(torch_tensor_type), INTENT(IN)                :: grid_weights
      TYPE(torch_tensor_type), INTENT(INOUT)             :: exc_tensor
      LOGICAL, INTENT(IN)                                :: read_exc
      REAL(KIND=dp), INTENT(OUT)                         :: exc

      TYPE(torch_tensor_type)                            :: exc_density

#if defined (__HAS_IEEE_EXCEPTIONS)
      LOGICAL, DIMENSION(5)                              :: ieee_halt

      CALL ieee_get_halting_mode(IEEE_ALL, ieee_halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL torch_model_forward_mol_tensor(model%torch_model, "get_exc_density", inputs, exc_density)
      CALL torch_tensor_weighted_sum(exc_density, grid_weights, exc_tensor)
      CALL torch_tensor_release(exc_density)
      exc = 0.0_dp
      IF (read_exc) exc = torch_tensor_item_double(exc_tensor)
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, ieee_halt)
#endif

   END SUBROUTINE skala_torch_model_get_exc_impl

! **************************************************************************************************
!> \brief Weighted SKALA energy per contiguous flat segment, sharing one model forward pass.
!> \param model ...
!> \param inputs ...
!> \param grid_weights ...
!> \param seg_offsets zero-based flat offsets of the segments
!> \param seg_lengths segment lengths in flat elements
!> \param exc_tensor total weighted energy (differentiable, used for backward)
!> \param exc_segments per-segment weighted energies
! **************************************************************************************************
   SUBROUTINE skala_torch_model_get_exc_segments(model, inputs, grid_weights, seg_offsets, &
                                                 seg_lengths, exc_tensor, exc_segments)
      TYPE(skala_torch_model_type), INTENT(INOUT)        :: model
      TYPE(torch_dict_type), INTENT(IN)                  :: inputs
      TYPE(torch_tensor_type), INTENT(IN)                :: grid_weights
      INTEGER(KIND=int_8), DIMENSION(:), INTENT(IN)      :: seg_offsets, seg_lengths
      TYPE(torch_tensor_type), INTENT(INOUT)             :: exc_tensor
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: exc_segments

      TYPE(torch_tensor_type)                            :: exc_density
#if defined (__HAS_IEEE_EXCEPTIONS)
      LOGICAL, DIMENSION(5)                              :: ieee_halt
#endif

      CPASSERT(SIZE(seg_offsets) == SIZE(seg_lengths))
      CPASSERT(SIZE(exc_segments) == SIZE(seg_offsets))

#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, ieee_halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL torch_model_forward_mol_tensor(model%torch_model, "get_exc_density", inputs, exc_density)
      CALL torch_tensor_weighted_sum(exc_density, grid_weights, exc_tensor)
      ! All segment sums are reduced on device and copied back with a single
      ! host synchronization instead of one per segment.
      CALL torch_tensor_weighted_sum_segments(exc_density, grid_weights, seg_offsets, &
                                              seg_lengths, exc_segments)
      CALL torch_tensor_release(exc_density)
#if defined (__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, ieee_halt)
#endif

   END SUBROUTINE skala_torch_model_get_exc_segments

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

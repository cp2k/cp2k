!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Owned inverse USPP data, independent of the KS methods that consume it.
!>        Bases remain borrowed. Explicit SCF activation requires a prepared common frame.
! **************************************************************************************************
MODULE qs_upf_inverse_types
   USE basis_set_types,                 ONLY: gto_basis_set_p_type
   USE kinds,                           ONLY: dp
   USE kpoint_transitional,             ONLY: kpoint_transitional_release,&
                                              kpoint_transitional_type
   USE message_passing,                 ONLY: mp_comm_type
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
   USE upf_inverse,                     ONLY: upf_inverse_atomic_type,&
                                              upf_inverse_type
   USE upf_radial_fields,               ONLY: upf_radial_fields_type
   USE uspp_augmentation,               ONLY: uspp_augmentation_type

   IMPLICIT NONE
   PRIVATE

   TYPE qs_upf_inverse_type
      TYPE(upf_inverse_type) :: inverse
      LOGICAL :: atom_local = .FALSE., direct = .FALSE.
      TYPE(upf_inverse_atomic_type), ALLOCATABLE :: atomic_cache(:)
      TYPE(mp_comm_type) :: point_comm
      ! Bases are borrowed from the live Quickstep environment; radial caches and fields are owned.
      TYPE(gto_basis_set_p_type), ALLOCATABLE :: bases(:)
      TYPE(upf_radial_fields_type), ALLOCATABLE :: beta_cache(:), delta_cache(:)
      INTEGER, ALLOCATABLE :: ao_first(:), ao_last(:), projector_first(:), projector_last(:), center_atoms(:), center_kinds(:)
      INTEGER, ALLOCATABLE :: images(:, :)
      REAL(dp) :: cell(3, 3) = 0.0_dp, reduced_k(3) = 0.0_dp
      REAL(dp), ALLOCATABLE :: positions(:, :), points(:, :), weights(:)
      ! Supplied PP_Q and the reference metric's radial Q integral are separate data.
      COMPLEX(dp), ALLOCATABLE :: coefficients(:, :), q(:, :), radial_q(:, :), dion(:, :), beta(:, :), delta(:, :), &
                                  delta_gradient(:, :, :), gaussian(:, :), gaussian_gradient(:, :, :), &
                                  soft(:, :), soft_gradient(:, :, :)
   END TYPE qs_upf_inverse_type

   TYPE qs_upf_inverse_grid_type
      TYPE(qs_upf_inverse_type) :: ao
      INTEGER :: bounds(2, 3) = 0, bounds_local(2, 3) = 0, npts(3) = 0, mode = PW_MODE_LOCAL
      REAL(dp) :: dh(3, 3) = 0.0_dp, dvol = 0.0_dp
   END TYPE qs_upf_inverse_grid_type

   TYPE qs_upf_inverse_augmentation_type
      TYPE(uspp_augmentation_type), ALLOCATABLE :: kinds(:)
   END TYPE qs_upf_inverse_augmentation_type

   TYPE qs_upf_inverse_scf_type
      LOGICAL :: enabled = .FALSE.
      INTEGER :: density_evaluations = 0, operator_evaluations = 0
      TYPE(qs_upf_inverse_grid_type) :: frame
      TYPE(qs_upf_inverse_augmentation_type) :: augmentation
      TYPE(kpoint_transitional_type) :: metric
      COMPLEX(dp), ALLOCATABLE :: gaussian_kinetic(:, :)
   END TYPE qs_upf_inverse_scf_type

   PUBLIC :: qs_upf_inverse_type, qs_upf_inverse_grid_type, qs_upf_inverse_augmentation_type
   PUBLIC :: qs_upf_inverse_scf_type, qs_upf_inverse_scf_release

CONTAINS

! **************************************************************************************************
!> \brief Release the owned metric and fields before their borrowed KS basis owner is destroyed.
!> \param state Owned SCF state; pointer is nullified.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_release(state)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      IF (.NOT. ASSOCIATED(state)) RETURN
      CALL kpoint_transitional_release(state%metric)
      DEALLOCATE (state)
      NULLIFY (state)
   END SUBROUTINE qs_upf_inverse_scf_release

END MODULE qs_upf_inverse_types

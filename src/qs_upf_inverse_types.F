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
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_release
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE kinds,                           ONLY: dp
   USE kpoint_transform,                ONLY: kp_transform_plan_type
   USE kpoint_transitional,             ONLY: kpoint_transitional_release,&
                                              kpoint_transitional_type
   USE kpoint_types,                    ONLY: kpoint_clear_scf_smat,&
                                              kpoint_native_density_clear,&
                                              kpoint_native_operator_clear,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_comm_type
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE upf_inverse,                     ONLY: upf_inverse_atomic_type,&
                                              upf_inverse_type
   USE upf_projector_integrals,         ONLY: upf_projector_type
   USE upf_radial_fields,               ONLY: upf_radial_fields_type
   USE uspp_augmentation,               ONLY: uspp_augmentation_type

   IMPLICIT NONE
   PRIVATE

   TYPE qs_upf_inverse_atom_grid_ref_type
      TYPE(grid_atom_type), POINTER :: grid => NULL()
      TYPE(upf_projector_type) :: projectors
      REAL(dp), ALLOCATABLE :: directions(:, :), angular_weights(:), home_coupling(:, :)
   END TYPE qs_upf_inverse_atom_grid_ref_type
   PUBLIC :: qs_upf_inverse_atom_grid_ref_type

   TYPE qs_upf_inverse_type
      TYPE(qs_upf_inverse_atom_grid_ref_type), ALLOCATABLE :: defining_atom_grids(:)
      TYPE(upf_inverse_type) :: inverse
      LOGICAL :: atom_local = .FALSE., direct = .FALSE.
      TYPE(upf_inverse_atomic_type), ALLOCATABLE :: atomic_cache(:)
      INTEGER :: coarse_points = -1, patch_points = 0
      INTEGER, ALLOCATABLE :: patch_center(:)
      REAL(dp), ALLOCATABLE :: metric_weights(:)
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
      LOGICAL :: native_gaussian = .FALSE.
      TYPE(qs_upf_inverse_type) :: ao
      INTEGER, ALLOCATABLE :: patch_index(:, :, :), patch_offset(:, :)
      REAL(dp), ALLOCATABLE :: stencil(:, :, :, :), dstencil(:, :, :, :)
      INTEGER :: bounds(2, 3) = 0, bounds_local(2, 3) = 0, npts(3) = 0, mode = PW_MODE_LOCAL
      REAL(dp) :: dh(3, 3) = 0.0_dp, dvol = 0.0_dp
   END TYPE qs_upf_inverse_grid_type

   TYPE qs_upf_inverse_augmentation_type
      TYPE(uspp_augmentation_type), ALLOCATABLE :: kinds(:)
   END TYPE qs_upf_inverse_augmentation_type

   TYPE qs_upf_inverse_gaussian_type
      TYPE(kp_transform_plan_type) :: plan
      TYPE(mp_comm_type) :: comm
      TYPE(dbcsr_p_type), POINTER :: density(:) => NULL()
      REAL(dp), ALLOCATABLE :: reduced_k(:, :), weights(:)
      INTEGER :: nao = 0
   END TYPE qs_upf_inverse_gaussian_type

   TYPE qs_upf_inverse_kpoints_type
      ! The live initialized group environment is borrowed; all fields below are owned.
      TYPE(kpoint_type), POINTER :: kpoints => NULL()
      TYPE(qs_upf_inverse_grid_type), ALLOCATABLE :: frames(:)
      TYPE(qs_upf_inverse_augmentation_type) :: augmentation
      TYPE(qs_upf_inverse_gaussian_type) :: gaussian
      REAL(dp), ALLOCATABLE :: reduced_k(:, :), weights(:)
      INTEGER :: nao = 0, kp_range(2) = 0
   END TYPE qs_upf_inverse_kpoints_type

   TYPE qs_upf_inverse_scf_type
      LOGICAL :: native_gaussian = .FALSE.
      TYPE(dbcsr_p_type), POINTER :: direct_density(:) => NULL(), direct_potential(:) => NULL()
      LOGICAL :: enabled = .FALSE.
      INTEGER :: density_evaluations = 0, operator_evaluations = 0
      TYPE(qs_upf_inverse_grid_type) :: frame
      TYPE(qs_upf_inverse_augmentation_type) :: augmentation
      TYPE(kpoint_transitional_type) :: metric
      COMPLEX(dp), ALLOCATABLE :: gaussian_kinetic(:, :)
      TYPE(qs_upf_inverse_kpoints_type) :: bloch
      COMPLEX(dp), ALLOCATABLE :: metric_k(:, :, :), core_k(:, :, :), kinetic_k(:, :, :), ks_k(:, :, :, :)
      COMPLEX(dp), ALLOCATABLE :: gaussian_kinetic_k(:, :, :)
   END TYPE qs_upf_inverse_scf_type

   PUBLIC :: qs_upf_inverse_type, qs_upf_inverse_grid_type, qs_upf_inverse_augmentation_type
   PUBLIC :: qs_upf_inverse_kpoints_type, qs_upf_inverse_scf_type, qs_upf_inverse_scf_release
   PUBLIC :: qs_upf_inverse_gaussian_type, qs_upf_inverse_gaussian_release

CONTAINS

! **************************************************************************************************
!> \brief Release the owned metric and fields before their borrowed KS basis owner is destroyed.
!> \param state Owned SCF state; pointer is nullified.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_release(state)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      INTEGER                                            :: slot

      IF (.NOT. ASSOCIATED(state)) RETURN
      IF (state%enabled .AND. ASSOCIATED(state%bloch%kpoints)) THEN
         DO slot = 1, SIZE(state%bloch%kpoints%native_operators)
            CALL kpoint_native_operator_clear(state%bloch%kpoints%native_operators(slot))
         END DO
         CALL kpoint_native_density_clear(state%bloch%kpoints)
         CALL kpoint_clear_scf_smat(state%bloch%kpoints)
      END IF
      CALL kpoint_transitional_release(state%metric)
      CALL qs_upf_inverse_gaussian_release(state%bloch%gaussian)
      IF (ASSOCIATED(state%direct_density)) CALL dbcsr_deallocate_matrix_set(state%direct_density)
      IF (ASSOCIATED(state%direct_potential)) CALL dbcsr_deallocate_matrix_set(state%direct_potential)
      DEALLOCATE (state)
      NULLIFY (state)
   END SUBROUTINE qs_upf_inverse_scf_release

! **************************************************************************************************
!> \brief Release the owned native image matrices and phase map.
!> \param map ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_gaussian_release(map)
      TYPE(qs_upf_inverse_gaussian_type), INTENT(INOUT)  :: map

      INTEGER                                            :: image

      IF (ASSOCIATED(map%density)) THEN
         DO image = 1, SIZE(map%density)
            CALL dbcsr_release(map%density(image)%matrix)
            DEALLOCATE (map%density(image)%matrix)
         END DO
         DEALLOCATE (map%density)
      END IF
      map = qs_upf_inverse_gaussian_type()
   END SUBROUTINE qs_upf_inverse_gaussian_release

END MODULE qs_upf_inverse_types

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Build SKALA TorchScript feature dictionaries from CP2K GPW real-space grids.
! **************************************************************************************************
MODULE skala_gpw_features
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_array_utils,                  ONLY: cp_3d_r_cp_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE torch_api,                       ONLY: torch_dict_create,&
                                              torch_dict_insert,&
                                              torch_dict_release,&
                                              torch_dict_type,&
                                              torch_tensor_from_array,&
                                              torch_tensor_release,&
                                              torch_tensor_to_device_leaf,&
                                              torch_tensor_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                              xc_rho_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_gpw_features'

   PUBLIC :: skala_gpw_feature_type, skala_gpw_feature_build, skala_gpw_feature_release

   TYPE skala_gpw_feature_type
      TYPE(torch_dict_type)                             :: inputs
      TYPE(torch_tensor_type)                           :: atomic_grid_size_bound_shape_t
      TYPE(torch_tensor_type)                           :: atomic_grid_sizes_t
      TYPE(torch_tensor_type)                           :: atomic_grid_weights_t
      TYPE(torch_tensor_type)                           :: coarse_0_atomic_coords_t
      TYPE(torch_tensor_type)                           :: density_t
      TYPE(torch_tensor_type)                           :: grad_t
      TYPE(torch_tensor_type)                           :: grid_coords_t
      TYPE(torch_tensor_type)                           :: grid_weights_t
      TYPE(torch_tensor_type)                           :: kin_t
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)          :: feature_index
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)    :: atomic_grid_sizes
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :) :: atomic_grid_size_bound_shape
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)          :: atomic_grid_weights, grid_weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)       :: coarse_0_atomic_coords, density, &
                                                           grid_coords, kin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)    :: grad
      REAL(KIND=dp)                                      :: electron_count = 0.0_dp, &
                                                            grid_weight_sum = 0.0_dp, &
                                                            spin_moment = 0.0_dp
      LOGICAL                                            :: active = .FALSE.
   END TYPE skala_gpw_feature_type

CONTAINS

! **************************************************************************************************
!> \brief Build a flat SKALA molecular feature dictionary from a local GPW grid.
!> \param features ...
!> \param rho_set ...
!> \param rho_r ...
!> \param particle_set ...
!> \param cell ...
!> \param requires_grad ...
!> \param weights ...
! **************************************************************************************************
   SUBROUTINE skala_gpw_feature_build(features, rho_set, rho_r, particle_set, cell, &
                                      requires_grad, weights)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(pw_r3d_rs_type), DIMENSION(:), INTENT(IN)     :: rho_r
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN), OPTIONAL                      :: requires_grad
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: weights

      INTEGER, PARAMETER                                 :: nreal_per_point = 14

      INTEGER                                            :: i, iatom, ipt, ispin, j, k, local_row, &
                                                            max_grid_size, natom, nflat, &
                                                            nflat_local, nproc, nspins, owner, pe, &
                                                            pe_index, real_base
      INTEGER, ALLOCATABLE, DIMENSION(:) :: atom_offset, atom_position, feature_counts, &
         feature_displs, global_owner, local_owner, local_to_global, real_counts, real_displs
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: my_requires_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: global_feature, local_feature
      REAL(KIND=dp), DIMENSION(3)                        :: grid_point, owner_coord
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho, rhoa, rhob, tau_a, tau_b, tau_total
      TYPE(cp_3d_r_cp_type), DIMENSION(3)                :: drho, drhoa, drhob
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      my_requires_grad = .FALSE.
      IF (PRESENT(requires_grad)) my_requires_grad = requires_grad

      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(SIZE(rho_r) == 1 .OR. SIZE(rho_r) == 2)
      CPASSERT(ASSOCIATED(rho_r(1)%pw_grid))
      pw_grid => rho_r(1)%pw_grid

      natom = SIZE(particle_set)
      nspins = SIZE(rho_r)
      bo = pw_grid%bounds_local
      nflat_local = pw_grid%ngpts_local
      nproc = pw_grid%para%group%num_pe
      pe_index = pw_grid%para%group%mepos + 1

      CALL skala_gpw_feature_release(features)
      ALLOCATE (local_owner(nflat_local), local_feature(nreal_per_point*nflat_local), &
                feature_counts(nproc), feature_displs(nproc), real_counts(nproc), real_displs(nproc))
      ALLOCATE (features%feature_index(bo(1, 1):bo(2, 1), &
                                       bo(1, 2):bo(2, 2), &
                                       bo(1, 3):bo(2, 3)))
      features%feature_index = 0
      local_feature = 0.0_dp

      IF (nspins == 1) THEN
         CALL xc_rho_set_get(rho_set, rho=rho, drho=drho, tau=tau_total)
      ELSE
         CALL xc_rho_set_get(rho_set, rhoa=rhoa, rhob=rhob, drhoa=drhoa, drhob=drhob, &
                             tau_a=tau_a, tau_b=tau_b)
      END IF

      local_row = 0
      DO k = bo(1, 3), bo(2, 3)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 1), bo(2, 1)
               local_row = local_row + 1
               real_base = nreal_per_point*(local_row - 1)
               grid_point = grid_coordinate(pw_grid, [i, j, k])
               owner = nearest_atom(grid_point, particle_set, cell)
               local_owner(local_row) = owner
               features%feature_index(i, j, k) = local_row

               owner_coord = pbc(particle_set(owner)%r, cell, positive_range=.TRUE.)
               local_feature(real_base + 11:real_base + 13) = owner_coord + &
                                                              pbc(owner_coord, grid_point, cell)
               local_feature(real_base + 14) = pw_grid%dvol
               IF (PRESENT(weights)) THEN
                  IF (ASSOCIATED(weights)) THEN
                     local_feature(real_base + 14) = pw_grid%dvol*weights%array(i, j, k)
                  END IF
               END IF

               IF (nspins == 1) THEN
                  local_feature(real_base + 1) = 0.5_dp*rho(i, j, k)
                  local_feature(real_base + 2) = 0.5_dp*rho(i, j, k)
                  DO ispin = 1, 2
                     local_feature(real_base + 2 + 3*(ispin - 1) + 1) = 0.5_dp*drho(1)%array(i, j, k)
                     local_feature(real_base + 2 + 3*(ispin - 1) + 2) = 0.5_dp*drho(2)%array(i, j, k)
                     local_feature(real_base + 2 + 3*(ispin - 1) + 3) = 0.5_dp*drho(3)%array(i, j, k)
                     local_feature(real_base + 8 + ispin) = 0.5_dp*tau_total(i, j, k)
                  END DO
               ELSE
                  local_feature(real_base + 1) = rhoa(i, j, k)
                  local_feature(real_base + 2) = rhob(i, j, k)
                  local_feature(real_base + 3) = drhoa(1)%array(i, j, k)
                  local_feature(real_base + 4) = drhoa(2)%array(i, j, k)
                  local_feature(real_base + 5) = drhoa(3)%array(i, j, k)
                  local_feature(real_base + 6) = drhob(1)%array(i, j, k)
                  local_feature(real_base + 7) = drhob(2)%array(i, j, k)
                  local_feature(real_base + 8) = drhob(3)%array(i, j, k)
                  local_feature(real_base + 9) = tau_a(i, j, k)
                  local_feature(real_base + 10) = tau_b(i, j, k)
               END IF
            END DO
         END DO
      END DO

      ! SKALA groups all grid points by atom. For distributed GPW grids, each rank
      ! reconstructs this same global feature ordering and later maps only its local
      ! grid slice back to VXC.
      CALL pw_grid%para%group%allgather(nflat_local, feature_counts)
      feature_displs(1) = 0
      DO pe = 2, nproc
         feature_displs(pe) = feature_displs(pe - 1) + feature_counts(pe - 1)
      END DO
      real_counts = nreal_per_point*feature_counts
      real_displs = nreal_per_point*feature_displs
      nflat = SUM(feature_counts)
      ALLOCATE (global_owner(nflat), global_feature(nreal_per_point*nflat))
      CALL pw_grid%para%group%allgatherv(local_owner, global_owner, feature_counts, feature_displs)
      CALL pw_grid%para%group%allgatherv(local_feature, global_feature, real_counts, real_displs)

      ALLOCATE (features%atomic_grid_sizes(natom), atom_offset(natom + 1), atom_position(natom), &
                local_to_global(nflat_local))
      features%atomic_grid_sizes = 0_int_8
      DO ipt = 1, nflat
         features%atomic_grid_sizes(global_owner(ipt)) = features%atomic_grid_sizes(global_owner(ipt)) + 1_int_8
      END DO
      atom_offset(1) = 1
      DO iatom = 1, natom
         atom_offset(iatom + 1) = atom_offset(iatom) + INT(features%atomic_grid_sizes(iatom))
      END DO
      atom_position = atom_offset(:natom)
      max_grid_size = MAXVAL(INT(features%atomic_grid_sizes))

      ALLOCATE (features%density(nflat, 2), features%grad(nflat, 3, 2), features%kin(nflat, 2), &
                features%grid_coords(3, nflat), features%grid_weights(nflat), &
                features%atomic_grid_weights(nflat), features%coarse_0_atomic_coords(3, natom), &
                features%atomic_grid_size_bound_shape(0, max_grid_size))
      features%density = 0.0_dp
      features%grad = 0.0_dp
      features%kin = 0.0_dp
      features%grid_coords = 0.0_dp
      features%grid_weights = pw_grid%dvol
      features%atomic_grid_weights = pw_grid%dvol

      DO iatom = 1, natom
         features%coarse_0_atomic_coords(:, iatom) = pbc(particle_set(iatom)%r, cell, positive_range=.TRUE.)
      END DO

      DO ipt = 1, nflat
         owner = global_owner(ipt)
         local_row = atom_position(owner)
         atom_position(owner) = atom_position(owner) + 1
         real_base = nreal_per_point*(ipt - 1)
         features%density(local_row, :) = global_feature(real_base + 1:real_base + 2)
         features%grad(local_row, 1, 1) = global_feature(real_base + 3)
         features%grad(local_row, 2, 1) = global_feature(real_base + 4)
         features%grad(local_row, 3, 1) = global_feature(real_base + 5)
         features%grad(local_row, 1, 2) = global_feature(real_base + 6)
         features%grad(local_row, 2, 2) = global_feature(real_base + 7)
         features%grad(local_row, 3, 2) = global_feature(real_base + 8)
         features%kin(local_row, :) = global_feature(real_base + 9:real_base + 10)
         features%grid_coords(:, local_row) = global_feature(real_base + 11:real_base + 13)
         features%grid_weights(local_row) = global_feature(real_base + 14)
         features%atomic_grid_weights(local_row) = features%grid_weights(local_row)
         IF (ipt > feature_displs(pe_index) .AND. ipt <= feature_displs(pe_index) + nflat_local) THEN
            local_to_global(ipt - feature_displs(pe_index)) = local_row
         END IF
      END DO

      DO k = bo(1, 3), bo(2, 3)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 1), bo(2, 1)
               features%feature_index(i, j, k) = local_to_global(features%feature_index(i, j, k))
            END DO
         END DO
      END DO

      features%electron_count = SUM((features%density(:, 1) + features%density(:, 2))* &
                                    features%grid_weights)
      features%spin_moment = SUM((features%density(:, 1) - features%density(:, 2))* &
                                 features%grid_weights)
      features%grid_weight_sum = SUM(features%grid_weights)

      CALL torch_dict_create(features%inputs)
      CALL add_feature_tensors(features, my_requires_grad)
      features%active = .TRUE.

      DEALLOCATE (atom_offset, atom_position, feature_counts, feature_displs, global_feature, &
                  global_owner, local_feature, local_owner, local_to_global, real_counts, real_displs)

   END SUBROUTINE skala_gpw_feature_build

! **************************************************************************************************
!> \brief Release Torch objects and backing arrays owned by a feature bundle.
!> \param features ...
! **************************************************************************************************
   SUBROUTINE skala_gpw_feature_release(features)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features

      IF (features%active) THEN
         CALL torch_tensor_release(features%density_t)
         CALL torch_tensor_release(features%grad_t)
         CALL torch_tensor_release(features%kin_t)
         CALL torch_tensor_release(features%grid_coords_t)
         CALL torch_tensor_release(features%grid_weights_t)
         CALL torch_tensor_release(features%atomic_grid_weights_t)
         CALL torch_tensor_release(features%atomic_grid_sizes_t)
         CALL torch_tensor_release(features%coarse_0_atomic_coords_t)
         CALL torch_tensor_release(features%atomic_grid_size_bound_shape_t)
         CALL torch_dict_release(features%inputs)
         features%active = .FALSE.
      END IF

      IF (ALLOCATED(features%density)) DEALLOCATE (features%density)
      IF (ALLOCATED(features%grad)) DEALLOCATE (features%grad)
      IF (ALLOCATED(features%kin)) DEALLOCATE (features%kin)
      IF (ALLOCATED(features%feature_index)) DEALLOCATE (features%feature_index)
      IF (ALLOCATED(features%grid_coords)) DEALLOCATE (features%grid_coords)
      IF (ALLOCATED(features%grid_weights)) DEALLOCATE (features%grid_weights)
      IF (ALLOCATED(features%atomic_grid_weights)) DEALLOCATE (features%atomic_grid_weights)
      IF (ALLOCATED(features%atomic_grid_sizes)) DEALLOCATE (features%atomic_grid_sizes)
      IF (ALLOCATED(features%coarse_0_atomic_coords)) DEALLOCATE (features%coarse_0_atomic_coords)
      IF (ALLOCATED(features%atomic_grid_size_bound_shape)) &
         DEALLOCATE (features%atomic_grid_size_bound_shape)

   END SUBROUTINE skala_gpw_feature_release

! **************************************************************************************************
!> \brief Insert all SKALA feature tensors into the Torch dictionary.
!> \param features ...
!> \param requires_grad ...
! **************************************************************************************************
   SUBROUTINE add_feature_tensors(features, requires_grad)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      LOGICAL, INTENT(IN)                                :: requires_grad

      CALL torch_tensor_from_array(features%density_t, features%density, requires_grad=requires_grad)
      CALL torch_tensor_to_device_leaf(features%density_t, requires_grad)
      CALL torch_dict_insert(features%inputs, "density", features%density_t)
      CALL torch_tensor_from_array(features%grad_t, features%grad, requires_grad=requires_grad)
      CALL torch_tensor_to_device_leaf(features%grad_t, requires_grad)
      CALL torch_dict_insert(features%inputs, "grad", features%grad_t)
      CALL torch_tensor_from_array(features%kin_t, features%kin, requires_grad=requires_grad)
      CALL torch_tensor_to_device_leaf(features%kin_t, requires_grad)
      CALL torch_dict_insert(features%inputs, "kin", features%kin_t)
      CALL torch_tensor_from_array(features%grid_coords_t, features%grid_coords)
      CALL torch_tensor_to_device_leaf(features%grid_coords_t, .FALSE.)
      CALL torch_dict_insert(features%inputs, "grid_coords", features%grid_coords_t)
      CALL torch_tensor_from_array(features%grid_weights_t, features%grid_weights)
      CALL torch_tensor_to_device_leaf(features%grid_weights_t, .FALSE.)
      CALL torch_dict_insert(features%inputs, "grid_weights", features%grid_weights_t)
      CALL torch_tensor_from_array(features%atomic_grid_weights_t, features%atomic_grid_weights)
      CALL torch_tensor_to_device_leaf(features%atomic_grid_weights_t, .FALSE.)
      CALL torch_dict_insert(features%inputs, "atomic_grid_weights", features%atomic_grid_weights_t)
      CALL torch_tensor_from_array(features%atomic_grid_sizes_t, features%atomic_grid_sizes)
      CALL torch_tensor_to_device_leaf(features%atomic_grid_sizes_t, .FALSE.)
      CALL torch_dict_insert(features%inputs, "atomic_grid_sizes", features%atomic_grid_sizes_t)
      CALL torch_tensor_from_array(features%coarse_0_atomic_coords_t, features%coarse_0_atomic_coords)
      CALL torch_tensor_to_device_leaf(features%coarse_0_atomic_coords_t, .FALSE.)
      CALL torch_dict_insert(features%inputs, "coarse_0_atomic_coords", features%coarse_0_atomic_coords_t)
      CALL torch_tensor_from_array(features%atomic_grid_size_bound_shape_t, &
                                   features%atomic_grid_size_bound_shape)
      CALL torch_tensor_to_device_leaf(features%atomic_grid_size_bound_shape_t, .FALSE.)
      CALL torch_dict_insert(features%inputs, "atomic_grid_size_bound_shape", &
                             features%atomic_grid_size_bound_shape_t)

   END SUBROUTINE add_feature_tensors

! **************************************************************************************************
!> \brief Return the Cartesian coordinate of a regular GPW grid point.
!> \param pw_grid ...
!> \param index ...
!> \return ...
! **************************************************************************************************
   FUNCTION grid_coordinate(pw_grid, index) RESULT(coord)
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      INTEGER, DIMENSION(3), INTENT(IN)                  :: index
      REAL(KIND=dp), DIMENSION(3)                        :: coord

      INTEGER, DIMENSION(3)                              :: relative_index

      relative_index = index - pw_grid%bounds(1, :)
      coord = REAL(relative_index(1), KIND=dp)*pw_grid%dh(:, 1) + &
              REAL(relative_index(2), KIND=dp)*pw_grid%dh(:, 2) + &
              REAL(relative_index(3), KIND=dp)*pw_grid%dh(:, 3)

   END FUNCTION grid_coordinate

! **************************************************************************************************
!> \brief Assign a grid point to the nearest periodic atom.
!> \param grid_point ...
!> \param particle_set ...
!> \param cell ...
!> \return ...
! **************************************************************************************************
   FUNCTION nearest_atom(grid_point, particle_set, cell) RESULT(owner)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: grid_point
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER                                            :: owner

      INTEGER                                            :: iatom
      REAL(KIND=dp)                                      :: best_r2, r2
      REAL(KIND=dp), DIMENSION(3)                        :: rij

      owner = 1
      best_r2 = HUGE(1.0_dp)
      DO iatom = 1, SIZE(particle_set)
         rij = pbc(grid_point, particle_set(iatom)%r, cell)
         r2 = SUM(rij**2)
         IF (r2 < best_r2) THEN
            best_r2 = r2
            owner = iatom
         END IF
      END DO

   END FUNCTION nearest_atom

END MODULE skala_gpw_features

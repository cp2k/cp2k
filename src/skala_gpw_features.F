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
   USE message_passing,                 ONLY: mp_comm_type
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE torch_api,                       ONLY: &
        torch_dict_clone, torch_dict_create, torch_dict_insert, torch_dict_release, &
        torch_dict_type, torch_tensor_expand_dim, torch_tensor_from_array, torch_tensor_narrow, &
        torch_tensor_release, torch_tensor_reset_from_array, torch_tensor_to_device_leaf, &
        torch_tensor_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                              xc_rho_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_gpw_features'
   REAL(KIND=dp), PARAMETER, PRIVATE    :: layout_tol = 1.0E-12_dp
   INTEGER, PARAMETER, PRIVATE          :: ndynamic_per_point = 10, nrks_dynamic_per_point = 5, &
                                           nstatic_per_point = 5, ngrad_per_point = 10
   INTEGER, PARAMETER, PUBLIC           :: skala_gpw_atom_partition_hard = 1, &
                                           skala_gpw_atom_partition_smooth = 2
   REAL(KIND=dp), PARAMETER, PRIVATE    :: smooth_partition_eps = 1.0E-12_dp

   PUBLIC :: skala_gpw_atom_subchunk_count, skala_gpw_feature_build, &
             skala_gpw_feature_build_atom_subchunk, skala_gpw_feature_release, &
             skala_gpw_feature_type, skala_gpw_smooth_partition_derivatives

   TYPE skala_gpw_layout_cache_type
      INTEGER                                            :: chunk_atom_begin = 1, chunk_atom_end = 0, &
                                                            chunk_feature_begin = 1, &
                                                            chunk_feature_count = 0, chunk_natom = 0, &
                                                            natom = 0, nflat = 0, nflat_local = 0, &
                                                            npoint = 0, nproc = 0, &
                                                            atom_partition = skala_gpw_atom_partition_hard
      INTEGER, DIMENSION(2, 3)                           :: bo = 0, bounds = 0
      INTEGER, DIMENSION(3)                              :: npts = 0
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dynamic_counts, dynamic_displs, &
                                                            chunk_feature_counts, chunk_feature_displs, &
                                                            chunk_grad_counts, chunk_grad_displs, &
                                                            feature_counts, feature_displs, &
                                                            feature_source_points, global_to_feature, &
                                                            local_feature_counts, local_feature_offsets, &
                                                            local_feature_points, local_feature_rows, &
                                                            route_grad_return_recv_counts, &
                                                            route_grad_return_recv_displs, &
                                                            route_grad_return_send_counts, &
                                                            route_grad_return_send_displs, &
                                                            route_local_dest, chunk_return_positions, &
                                                            route_point_recv_counts, &
                                                            route_point_recv_displs, &
                                                            route_point_send_counts, &
                                                            route_point_send_displs, &
                                                            route_send_local_rows
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: feature_index
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: atomic_grid_sizes, chunk_atomic_grid_sizes, &
                                                            chunk_feature_indices
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: local_feature_indices
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :)  :: atomic_grid_size_bound_shape, &
                                                            chunk_atomic_grid_size_bound_shape
      TYPE(torch_dict_type)                              :: chunk_inputs
      TYPE(torch_dict_type)                              :: chunk_static_inputs
      TYPE(torch_dict_type)                              :: inputs
      TYPE(torch_dict_type)                              :: static_inputs
      TYPE(torch_tensor_type)                            :: atomic_grid_size_bound_shape_t
      TYPE(torch_tensor_type)                            :: atomic_grid_sizes_t
      TYPE(torch_tensor_type)                            :: atomic_grid_weights_t
      TYPE(torch_tensor_type)                            :: chunk_atomic_grid_size_bound_shape_t
      TYPE(torch_tensor_type)                            :: chunk_atomic_grid_sizes_t
      TYPE(torch_tensor_type)                            :: chunk_atomic_grid_weights_t
      TYPE(torch_tensor_type)                            :: chunk_coarse_0_atomic_coords_t
      TYPE(torch_tensor_type)                            :: chunk_density_t
      TYPE(torch_tensor_type)                            :: chunk_density_input_t
      TYPE(torch_tensor_type)                            :: chunk_feature_indices_t
      TYPE(torch_tensor_type)                            :: chunk_grad_t
      TYPE(torch_tensor_type)                            :: chunk_grad_input_t
      TYPE(torch_tensor_type)                            :: chunk_grid_coords_t
      TYPE(torch_tensor_type)                            :: chunk_grid_weights_t
      TYPE(torch_tensor_type)                            :: chunk_kin_t
      TYPE(torch_tensor_type)                            :: chunk_kin_input_t
      TYPE(torch_tensor_type)                            :: coarse_0_atomic_coords_t
      TYPE(torch_tensor_type)                            :: density_t
      TYPE(torch_tensor_type)                            :: grid_coords_t
      TYPE(torch_tensor_type)                            :: grid_weights_t
      TYPE(torch_tensor_type)                            :: grad_t
      TYPE(torch_tensor_type)                            :: kin_t
      TYPE(torch_tensor_type)                            :: local_feature_indices_t
      REAL(KIND=dp)                                      :: dvol = 0.0_dp, weight_sum = 0.0_dp, &
                                                            weight_sumsq = 0.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)                     :: cell_hmat = 0.0_dp, dh = 0.0_dp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: atomic_grid_weights, chunk_atomic_grid_weights, &
                                                            chunk_grid_weights, grid_weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_coords, chunk_coarse_0_atomic_coords, &
                                                            chunk_grid_coords, coarse_0_atomic_coords, &
                                                            grid_coords
      LOGICAL                                            :: active = .FALSE., has_weights = .FALSE., &
                                                            chunk_dynamic_input_views_active = .FALSE., &
                                                            chunk_dynamic_tensors_active = .FALSE., &
                                                            chunk_inputs_active = .FALSE., &
                                                            chunk_inputs_use_collapsed_rks = .FALSE., &
                                                            chunk_static_tensors_active = .FALSE., &
                                                            dynamic_tensors_active = .FALSE., &
                                                            inputs_active = .FALSE., &
                                                            static_tensors_active = .FALSE.
   END TYPE skala_gpw_layout_cache_type

   TYPE skala_gpw_feature_type
      INTEGER                                            :: chunk_feature_count = 0, nflat = 0, &
                                                            nflat_local = 0, &
                                                            atom_partition = skala_gpw_atom_partition_hard
      TYPE(torch_dict_type)                             :: inputs
      TYPE(torch_tensor_type)                           :: atomic_grid_size_bound_shape_t
      TYPE(torch_tensor_type)                           :: atomic_grid_sizes_t
      TYPE(torch_tensor_type)                           :: atomic_grid_weights_t
      TYPE(torch_tensor_type)                           :: coarse_0_atomic_coords_t
      TYPE(torch_tensor_type)                           :: density_input_t
      TYPE(torch_tensor_type)                           :: density_t
      TYPE(torch_tensor_type)                           :: grad_t
      TYPE(torch_tensor_type)                           :: grad_input_t
      TYPE(torch_tensor_type)                           :: grid_coords_t
      TYPE(torch_tensor_type)                           :: grid_weights_t
      TYPE(torch_tensor_type)                           :: kin_input_t
      TYPE(torch_tensor_type)                           :: kin_t
      TYPE(torch_tensor_type)                           :: local_feature_indices_t
      INTEGER, ALLOCATABLE, DIMENSION(:)                :: chunk_grad_counts, chunk_grad_displs, &
                                                           local_feature_counts, local_feature_offsets, &
                                                           local_feature_rows, &
                                                           chunk_return_positions, &
                                                           route_grad_return_recv_counts, &
                                                           route_grad_return_recv_displs, &
                                                           route_grad_return_send_counts, &
                                                           route_grad_return_send_displs, &
                                                           route_point_recv_counts, &
                                                           route_point_recv_displs, &
                                                           route_point_send_counts, &
                                                           route_point_send_displs, &
                                                           route_send_local_rows
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)          :: feature_index
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)    :: atomic_grid_sizes
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :) :: atomic_grid_size_bound_shape
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)          :: atomic_grid_weights, grid_weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)       :: chunk_density, chunk_kin, &
                                                           coarse_0_atomic_coords, density, &
                                                           grid_coords, kin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)    :: chunk_grad, grad
      REAL(KIND=dp)                                      :: electron_count = 0.0_dp, &
                                                            grid_weight_sum = 0.0_dp, &
                                                            spin_moment = 0.0_dp
      LOGICAL                                            :: active = .FALSE., owns_coordinate_tensor = .FALSE., &
                                                            owns_grid_coordinate_tensor = .FALSE., &
                                                            owns_weight_tensors = .FALSE., &
                                                            owns_dynamic_tensors = .TRUE., &
                                                            owns_inputs = .TRUE., &
                                                            owns_static_tensors = .TRUE., &
                                                            uses_atom_chunk_routing = .FALSE., &
                                                            uses_atom_chunks = .FALSE., &
                                                            uses_collapsed_rks_dynamic = .FALSE.
   END TYPE skala_gpw_feature_type

   TYPE(skala_gpw_layout_cache_type), SAVE               :: cached_layout

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
!> \param requires_coordinate_grad ...
!> \param requires_stress_grad ...
!> \param use_atom_chunks ...
!> \param route_atom_chunks ...
!> \param atom_partition ...
! **************************************************************************************************
   SUBROUTINE skala_gpw_feature_build(features, rho_set, rho_r, particle_set, cell, &
                                      requires_grad, weights, requires_coordinate_grad, &
                                      requires_stress_grad, use_atom_chunks, route_atom_chunks, &
                                      atom_partition)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(pw_r3d_rs_type), DIMENSION(:), INTENT(IN)     :: rho_r
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN), OPTIONAL                      :: requires_grad
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: weights
      LOGICAL, INTENT(IN), OPTIONAL                      :: requires_coordinate_grad, &
                                                            requires_stress_grad, use_atom_chunks, &
                                                            route_atom_chunks
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_partition

      INTEGER :: handle, i, ipt, ispin, j, k, local_row, my_atom_partition, &
         ndynamic_local_per_point, nflat, nflat_local, nspins, phase_handle, real_base, row
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL :: collapse_spin_dynamics, my_requires_coordinate_grad, my_requires_grad, &
         my_requires_stress_grad, my_route_atom_chunks, my_use_atom_chunks, &
         use_atom_chunk_protocol, use_atom_chunk_routing
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: global_dynamic, local_dynamic
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho, rhoa, rhob, tau_a, tau_b, tau_total
      TYPE(cp_3d_r_cp_type), DIMENSION(3)                :: drho, drhoa, drhob
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      CALL timeset("skala_gpw_feature_build", handle)

      my_requires_grad = .FALSE.
      IF (PRESENT(requires_grad)) my_requires_grad = requires_grad
      my_requires_coordinate_grad = .FALSE.
      IF (PRESENT(requires_coordinate_grad)) THEN
         my_requires_coordinate_grad = requires_coordinate_grad
      END IF
      my_requires_stress_grad = .FALSE.
      IF (PRESENT(requires_stress_grad)) my_requires_stress_grad = requires_stress_grad
      my_use_atom_chunks = .FALSE.
      IF (PRESENT(use_atom_chunks)) my_use_atom_chunks = use_atom_chunks
      my_route_atom_chunks = .FALSE.
      IF (PRESENT(route_atom_chunks)) my_route_atom_chunks = route_atom_chunks
      my_atom_partition = skala_gpw_atom_partition_hard
      IF (PRESENT(atom_partition)) my_atom_partition = atom_partition
      IF (my_atom_partition /= skala_gpw_atom_partition_hard .AND. &
          my_atom_partition /= skala_gpw_atom_partition_smooth) THEN
         CALL cp_abort(__LOCATION__, "Unknown native SKALA atom-partition mode.")
      END IF
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(SIZE(rho_r) == 1 .OR. SIZE(rho_r) == 2)
      CPASSERT(ASSOCIATED(rho_r(1)%pw_grid))
      pw_grid => rho_r(1)%pw_grid

      nspins = SIZE(rho_r)
      bo = pw_grid%bounds_local
      nflat_local = pw_grid%ngpts_local

      CALL timeset("skala_gpw_pre_release", phase_handle)
      CALL skala_gpw_feature_release(features)
      CALL timestop(phase_handle)

      CALL timeset("skala_gpw_layout_cache", phase_handle)
      CALL ensure_layout_cache(pw_grid, particle_set, cell, weights, my_atom_partition)
      CALL timestop(phase_handle)
      nflat = cached_layout%nflat
      use_atom_chunk_protocol = my_use_atom_chunks .AND. &
                                .NOT. (my_requires_coordinate_grad .OR. my_requires_stress_grad)
      use_atom_chunk_routing = use_atom_chunk_protocol .AND. my_route_atom_chunks
      collapse_spin_dynamics = nspins == 1 .AND. use_atom_chunk_routing
      ndynamic_local_per_point = ndynamic_per_point
      IF (collapse_spin_dynamics) ndynamic_local_per_point = nrks_dynamic_per_point
      ALLOCATE (local_dynamic(ndynamic_local_per_point*nflat_local))
      local_dynamic = 0.0_dp

      CALL timeset("skala_gpw_pack_local", phase_handle)
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
               real_base = ndynamic_local_per_point*(local_row - 1)

               IF (nspins == 1) THEN
                  IF (collapse_spin_dynamics) THEN
                     local_dynamic(real_base + 1) = 0.5_dp*rho(i, j, k)
                     local_dynamic(real_base + 2) = 0.5_dp*drho(1)%array(i, j, k)
                     local_dynamic(real_base + 3) = 0.5_dp*drho(2)%array(i, j, k)
                     local_dynamic(real_base + 4) = 0.5_dp*drho(3)%array(i, j, k)
                     local_dynamic(real_base + 5) = 0.5_dp*tau_total(i, j, k)
                  ELSE
                     local_dynamic(real_base + 1) = 0.5_dp*rho(i, j, k)
                     local_dynamic(real_base + 2) = 0.5_dp*rho(i, j, k)
                     DO ispin = 1, 2
                        local_dynamic(real_base + 2 + 3*(ispin - 1) + 1) = &
                           0.5_dp*drho(1)%array(i, j, k)
                        local_dynamic(real_base + 2 + 3*(ispin - 1) + 2) = &
                           0.5_dp*drho(2)%array(i, j, k)
                        local_dynamic(real_base + 2 + 3*(ispin - 1) + 3) = &
                           0.5_dp*drho(3)%array(i, j, k)
                        local_dynamic(real_base + 8 + ispin) = 0.5_dp*tau_total(i, j, k)
                     END DO
                  END IF
               ELSE
                  local_dynamic(real_base + 1) = rhoa(i, j, k)
                  local_dynamic(real_base + 2) = rhob(i, j, k)
                  local_dynamic(real_base + 3) = drhoa(1)%array(i, j, k)
                  local_dynamic(real_base + 4) = drhoa(2)%array(i, j, k)
                  local_dynamic(real_base + 5) = drhoa(3)%array(i, j, k)
                  local_dynamic(real_base + 6) = drhob(1)%array(i, j, k)
                  local_dynamic(real_base + 7) = drhob(2)%array(i, j, k)
                  local_dynamic(real_base + 8) = drhob(3)%array(i, j, k)
                  local_dynamic(real_base + 9) = tau_a(i, j, k)
                  local_dynamic(real_base + 10) = tau_b(i, j, k)
               END IF
            END DO
         END DO
      END DO
      CALL timestop(phase_handle)

      CALL timeset("skala_gpw_copy_layout", phase_handle)
      CALL copy_cached_layout(features, my_requires_coordinate_grad .OR. my_requires_stress_grad, &
                              my_requires_stress_grad .OR. &
                              (my_atom_partition == skala_gpw_atom_partition_smooth .AND. &
                               (my_requires_coordinate_grad .OR. my_requires_stress_grad)))
      CALL timestop(phase_handle)

      IF (use_atom_chunk_routing) THEN
         CALL timeset("skala_gpw_route_dyn", phase_handle)
         CALL route_atom_chunk_dynamics(features, local_dynamic, pw_grid%para%group, &
                                        collapse_spin_dynamics)
         features%uses_atom_chunk_routing = .TRUE.
         features%uses_atom_chunks = .TRUE.
         CALL timestop(phase_handle)
      ELSE
         ALLOCATE (global_dynamic(ndynamic_per_point*cached_layout%npoint))
         CALL timeset("skala_gpw_allgatherv", phase_handle)
         CALL pw_grid%para%group%allgatherv(local_dynamic, global_dynamic, &
                                            cached_layout%dynamic_counts, &
                                            cached_layout%dynamic_displs)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_reorder_dyn", phase_handle)
         ALLOCATE (features%density(nflat, 2), features%grad(nflat, 3, 2), &
                   features%kin(nflat, 2))
         features%density = 0.0_dp
         features%grad = 0.0_dp
         features%kin = 0.0_dp

         DO row = 1, nflat
            ipt = cached_layout%feature_source_points(row)
            real_base = ndynamic_per_point*(ipt - 1)
            features%density(row, :) = global_dynamic(real_base + 1:real_base + 2)
            features%grad(row, 1, 1) = global_dynamic(real_base + 3)
            features%grad(row, 2, 1) = global_dynamic(real_base + 4)
            features%grad(row, 3, 1) = global_dynamic(real_base + 5)
            features%grad(row, 1, 2) = global_dynamic(real_base + 6)
            features%grad(row, 2, 2) = global_dynamic(real_base + 7)
            features%grad(row, 3, 2) = global_dynamic(real_base + 8)
            features%kin(row, :) = global_dynamic(real_base + 9:real_base + 10)
         END DO
         CALL timestop(phase_handle)
      END IF

      CALL timeset("skala_gpw_feature_sums", phase_handle)
      IF (features%uses_atom_chunks) THEN
         features%electron_count = 0.0_dp
         features%spin_moment = 0.0_dp
         IF (features%chunk_feature_count > 0) THEN
            IF (features%uses_collapsed_rks_dynamic) THEN
               features%electron_count = SUM(2.0_dp*features%chunk_density(:, 1)* &
                                             cached_layout%chunk_grid_weights)
            ELSE
               features%electron_count = SUM((features%chunk_density(:, 1) + &
                                              features%chunk_density(:, 2))* &
                                             cached_layout%chunk_grid_weights)
               features%spin_moment = SUM((features%chunk_density(:, 1) - &
                                           features%chunk_density(:, 2))* &
                                          cached_layout%chunk_grid_weights)
            END IF
         END IF
         CALL pw_grid%para%group%sum(features%electron_count)
         CALL pw_grid%para%group%sum(features%spin_moment)
      ELSE
         features%electron_count = SUM((features%density(:, 1) + features%density(:, 2))* &
                                       features%grid_weights)
         features%spin_moment = SUM((features%density(:, 1) - features%density(:, 2))* &
                                    features%grid_weights)
      END IF
      features%grid_weight_sum = SUM(features%grid_weights)
      CALL timestop(phase_handle)

      CALL timeset("skala_gpw_tensor_update", phase_handle)
      IF (use_atom_chunk_protocol .AND. .NOT. features%uses_atom_chunks) THEN
         IF (features%chunk_feature_count > 0) CALL extract_atom_chunk_dynamics(features)
         features%uses_atom_chunks = .TRUE.
      END IF
      IF (.NOT. features%uses_atom_chunks .OR. features%chunk_feature_count > 0) THEN
         CALL add_feature_tensors(features, my_requires_grad, my_requires_coordinate_grad, &
                                  my_requires_stress_grad, &
                                  features%uses_atom_chunks, &
                                  requires_weight_grad= &
                                  (my_atom_partition == skala_gpw_atom_partition_smooth .AND. &
                                   (my_requires_coordinate_grad .OR. my_requires_stress_grad)))
      ELSE
         ! This rank participates in atom-chunk communication but owns no model input rows.
         features%owns_coordinate_tensor = .FALSE.
         features%owns_grid_coordinate_tensor = .FALSE.
         features%owns_weight_tensors = .FALSE.
         features%owns_dynamic_tensors = .FALSE.
         features%owns_inputs = .FALSE.
         features%owns_static_tensors = .FALSE.
      END IF
      CALL timestop(phase_handle)
      features%active = .TRUE.

      IF (ALLOCATED(global_dynamic)) DEALLOCATE (global_dynamic)
      DEALLOCATE (local_dynamic)
      CALL timestop(handle)

   END SUBROUTINE skala_gpw_feature_build

! **************************************************************************************************
!> \brief Ensure that static grid-to-atom layout data is cached for the current grid/geometry.
!> \param pw_grid ...
!> \param particle_set ...
!> \param cell ...
!> \param weights ...
!> \param atom_partition ...
! **************************************************************************************************
   SUBROUTINE ensure_layout_cache(pw_grid, particle_set, cell, weights, atom_partition)
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: weights
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_partition

      INTEGER                                            :: my_atom_partition, phase_handle
      LOGICAL                                            :: cache_matches

      my_atom_partition = skala_gpw_atom_partition_hard
      IF (PRESENT(atom_partition)) my_atom_partition = atom_partition
      IF (PRESENT(weights)) THEN
         CALL timeset("skala_gpw_layout_match", phase_handle)
         cache_matches = layout_cache_matches(pw_grid, particle_set, cell, weights, &
                                              my_atom_partition)
         CALL timestop(phase_handle)
         IF (cache_matches) RETURN
         CALL timeset("skala_gpw_layout_rebuild", phase_handle)
         CALL rebuild_layout_cache(pw_grid, particle_set, cell, weights, my_atom_partition)
         CALL timestop(phase_handle)
      ELSE
         CALL timeset("skala_gpw_layout_match", phase_handle)
         cache_matches = layout_cache_matches(pw_grid, particle_set, cell, &
                                              atom_partition=my_atom_partition)
         CALL timestop(phase_handle)
         IF (cache_matches) RETURN
         CALL timeset("skala_gpw_layout_rebuild", phase_handle)
         CALL rebuild_layout_cache(pw_grid, particle_set, cell, &
                                   atom_partition=my_atom_partition)
         CALL timestop(phase_handle)
      END IF

   END SUBROUTINE ensure_layout_cache

! **************************************************************************************************
!> \brief Check whether the current static layout cache can be reused.
!> \param pw_grid ...
!> \param particle_set ...
!> \param cell ...
!> \param weights ...
!> \param atom_partition ...
!> \return ...
! **************************************************************************************************
   FUNCTION layout_cache_matches(pw_grid, particle_set, cell, weights, atom_partition) RESULT(matches)
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: weights
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_partition
      LOGICAL                                            :: matches

      INTEGER                                            :: iatom, my_atom_partition
      LOGICAL                                            :: weights_match

      my_atom_partition = skala_gpw_atom_partition_hard
      IF (PRESENT(atom_partition)) my_atom_partition = atom_partition
      matches = .FALSE.
      IF (.NOT. cached_layout%active) RETURN
      IF (cached_layout%atom_partition /= my_atom_partition) RETURN
      IF (cached_layout%natom /= SIZE(particle_set)) RETURN
      IF (cached_layout%nflat_local /= pw_grid%ngpts_local) RETURN
      IF (cached_layout%nproc /= pw_grid%para%group%num_pe) RETURN
      IF (ANY(cached_layout%bo /= pw_grid%bounds_local)) RETURN
      IF (ANY(cached_layout%bounds /= pw_grid%bounds)) RETURN
      IF (ANY(cached_layout%npts /= pw_grid%npts)) RETURN
      IF (ABS(cached_layout%dvol - pw_grid%dvol) > layout_tol) RETURN
      IF (ANY(ABS(cached_layout%dh - pw_grid%dh) > layout_tol)) RETURN
      IF (ANY(ABS(cached_layout%cell_hmat - cell%hmat) > layout_tol)) RETURN
      IF (.NOT. ALLOCATED(cached_layout%atom_coords)) RETURN

      DO iatom = 1, SIZE(particle_set)
         IF (ANY(ABS(cached_layout%atom_coords(:, iatom) - particle_set(iatom)%r) > layout_tol)) RETURN
      END DO

      IF (PRESENT(weights)) THEN
         weights_match = layout_weights_match(pw_grid, weights)
      ELSE
         weights_match = layout_weights_match(pw_grid)
      END IF
      IF (.NOT. weights_match) RETURN

      matches = .TRUE.

   END FUNCTION layout_cache_matches

! **************************************************************************************************
!> \brief Check whether current optional integration weights match the cached static tensors.
!> \param pw_grid ...
!> \param weights ...
!> \return ...
! **************************************************************************************************
   FUNCTION layout_weights_match(pw_grid, weights) RESULT(matches)
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: weights
      LOGICAL                                            :: matches

      LOGICAL                                            :: has_weights
      REAL(KIND=dp)                                      :: weight_sum, weight_sumsq

      matches = .FALSE.
      MARK_USED(pw_grid)
      IF (PRESENT(weights)) THEN
         CALL weights_signature(weights, has_weights, weight_sum, weight_sumsq)
      ELSE
         CALL weights_signature(has_weights=has_weights, weight_sum=weight_sum, &
                                weight_sumsq=weight_sumsq)
      END IF

      IF (cached_layout%has_weights .NEQV. has_weights) RETURN
      IF (ABS(cached_layout%weight_sum - weight_sum) > layout_tol) RETURN
      IF (ABS(cached_layout%weight_sumsq - weight_sumsq) > layout_tol) RETURN

      matches = .TRUE.

   END FUNCTION layout_weights_match

! **************************************************************************************************
!> \brief Build the static SKALA layout cache.
!> \param pw_grid ...
!> \param particle_set ...
!> \param cell ...
!> \param weights ...
!> \param atom_partition ...
! **************************************************************************************************
   SUBROUTINE rebuild_layout_cache(pw_grid, particle_set, cell, weights, atom_partition)
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: weights
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_partition

      INTEGER :: feature_local, i, iatom, ipt, j, k, local_row, max_grid_size, max_local_features, &
         my_atom_partition, natom, nfeature_local, nflat, nflat_local, npoint, nproc, owner, pe, &
         pe_index, phase_handle, row, source_global, source_local, static_base
      INTEGER, ALLOCATABLE, DIMENSION(:) :: atom_offset, atom_position, chunk_atom_begin, &
         chunk_atom_end, cursor, feature_counts, feature_displs, global_owner, &
         global_source_points, local_feature_counts_tmp, local_owner, local_source_global, &
         local_source_points, point_counts, point_displs, static_counts, static_displs
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: has_weights
      REAL(KIND=dp)                                      :: base_weight, included_sum, &
                                                            partition_weight, weight_sum, &
                                                            weight_sumsq
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: distances, global_static, local_static, &
                                                            partition_weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_coords_pbc, atom_image_coords
      REAL(KIND=dp), DIMENSION(3)                        :: grid_point, owner_coord

      CALL release_layout_cache(cached_layout)

      my_atom_partition = skala_gpw_atom_partition_hard
      IF (PRESENT(atom_partition)) my_atom_partition = atom_partition
      natom = SIZE(particle_set)
      bo = pw_grid%bounds_local
      nflat_local = pw_grid%ngpts_local
      nproc = pw_grid%para%group%num_pe
      pe_index = pw_grid%para%group%mepos + 1

      IF (PRESENT(weights)) THEN
         CALL weights_signature(weights, has_weights, weight_sum, weight_sumsq)
      ELSE
         CALL weights_signature(has_weights=has_weights, weight_sum=weight_sum, &
                                weight_sumsq=weight_sumsq)
      END IF

      max_local_features = nflat_local
      IF (my_atom_partition == skala_gpw_atom_partition_smooth) THEN
         max_local_features = nflat_local*natom
      END IF
      ALLOCATE (local_owner(max_local_features), &
                local_source_points(max_local_features), &
                local_static(nstatic_per_point*max_local_features), &
                local_feature_counts_tmp(nflat_local), feature_counts(nproc), &
                feature_displs(nproc), point_counts(nproc), point_displs(nproc), &
                static_counts(nproc), static_displs(nproc), atom_coords_pbc(3, natom), &
                atom_image_coords(3, natom), distances(natom), partition_weights(natom))
      ALLOCATE (cached_layout%feature_index(bo(1, 1):bo(2, 1), &
                                            bo(1, 2):bo(2, 2), &
                                            bo(1, 3):bo(2, 3)))
      cached_layout%feature_index = 0
      local_static = 0.0_dp
      local_feature_counts_tmp = 0
      DO iatom = 1, natom
         atom_coords_pbc(:, iatom) = pbc(particle_set(iatom)%r, cell, positive_range=.TRUE.)
      END DO

      CALL timeset("skala_gpw_layout_local", phase_handle)
      local_row = 0
      nfeature_local = 0
      DO k = bo(1, 3), bo(2, 3)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 1), bo(2, 1)
               local_row = local_row + 1
               grid_point = grid_coordinate(pw_grid, [i, j, k])
               base_weight = pw_grid%dvol
               IF (PRESENT(weights)) THEN
                  IF (ASSOCIATED(weights)) base_weight = base_weight*weights%array(i, j, k)
               END IF
               cached_layout%feature_index(i, j, k) = local_row

               IF (my_atom_partition == skala_gpw_atom_partition_hard) THEN
                  owner = nearest_atom(grid_point, atom_coords_pbc, cell)
                  owner_coord = atom_coords_pbc(:, owner)
                  nfeature_local = nfeature_local + 1
                  local_feature_counts_tmp(local_row) = 1
                  local_owner(nfeature_local) = owner
                  local_source_points(nfeature_local) = local_row
                  static_base = nstatic_per_point*(nfeature_local - 1)
                  local_static(static_base + 1:static_base + 3) = grid_point
                  local_static(static_base + 4) = base_weight
                  local_static(static_base + 5) = base_weight
               ELSE
                  CALL smooth_atom_partition(grid_point, atom_coords_pbc, cell, &
                                             partition_weights, atom_image_coords, distances)
                  included_sum = SUM(partition_weights, MASK=partition_weights > smooth_partition_eps)
                  IF (included_sum <= 0.0_dp) THEN
                     owner = nearest_atom(grid_point, atom_coords_pbc, cell)
                     partition_weights = 0.0_dp
                     partition_weights(owner) = 1.0_dp
                     included_sum = 1.0_dp
                  END IF
                  DO iatom = 1, natom
                     IF (partition_weights(iatom) <= smooth_partition_eps) CYCLE
                     partition_weight = partition_weights(iatom)/included_sum
                     nfeature_local = nfeature_local + 1
                     local_feature_counts_tmp(local_row) = &
                        local_feature_counts_tmp(local_row) + 1
                     local_owner(nfeature_local) = iatom
                     local_source_points(nfeature_local) = local_row
                     static_base = nstatic_per_point*(nfeature_local - 1)
                     local_static(static_base + 1:static_base + 3) = grid_point
                     local_static(static_base + 4) = base_weight*partition_weight
                     local_static(static_base + 5) = base_weight
                  END DO
               END IF
            END DO
         END DO
      END DO
      CALL timestop(phase_handle)

      ! SKALA groups all grid points by atom. This ordering is static while the
      ! grid, cell, atom positions, and optional integration weights are unchanged.
      CALL timeset("skala_gpw_layout_gather", phase_handle)
      CALL pw_grid%para%group%allgather(nflat_local, point_counts)
      CALL counts_to_displs(point_counts, point_displs)
      npoint = SUM(point_counts)
      CALL pw_grid%para%group%allgather(nfeature_local, feature_counts)
      CALL counts_to_displs(feature_counts, feature_displs)
      DO pe = 1, nproc
         static_counts(pe) = nstatic_per_point*feature_counts(pe)
         static_displs(pe) = nstatic_per_point*feature_displs(pe)
      END DO
      nflat = SUM(feature_counts)
      ALLOCATE (global_owner(nflat), global_source_points(nflat), &
                global_static(nstatic_per_point*nflat), local_source_global(nfeature_local))
      DO feature_local = 1, nfeature_local
         local_source_global(feature_local) = point_displs(pe_index) + local_source_points(feature_local)
      END DO
      CALL pw_grid%para%group%allgatherv(local_owner(1:nfeature_local), global_owner, feature_counts, &
                                         feature_displs)
      CALL pw_grid%para%group%allgatherv(local_source_global, global_source_points, feature_counts, &
                                         feature_displs)
      CALL pw_grid%para%group%allgatherv(local_static(1:nstatic_per_point*nfeature_local), &
                                         global_static, static_counts, &
                                         static_displs)
      CALL timestop(phase_handle)

      ALLOCATE (cached_layout%chunk_feature_counts(nproc), &
                cached_layout%chunk_feature_displs(nproc), &
                cached_layout%chunk_grad_counts(nproc), cached_layout%chunk_grad_displs(nproc), &
                cached_layout%feature_counts(nproc), cached_layout%feature_displs(nproc), &
                cached_layout%dynamic_counts(nproc), cached_layout%dynamic_displs(nproc), &
                cached_layout%route_grad_return_recv_counts(nproc), &
                cached_layout%route_grad_return_recv_displs(nproc), &
                cached_layout%route_grad_return_send_counts(nproc), &
                cached_layout%route_grad_return_send_displs(nproc), &
                cached_layout%route_point_recv_counts(nproc), &
                cached_layout%route_point_recv_displs(nproc), &
                cached_layout%route_point_send_counts(nproc), &
                cached_layout%route_point_send_displs(nproc), &
                cached_layout%feature_source_points(nflat), &
                cached_layout%global_to_feature(npoint), cached_layout%atomic_grid_sizes(natom), &
                cached_layout%local_feature_counts(nflat_local), &
                cached_layout%local_feature_offsets(nflat_local + 1), &
                cached_layout%local_feature_rows(nfeature_local), &
                cached_layout%local_feature_points(nfeature_local), &
                cached_layout%local_feature_indices(nfeature_local), atom_offset(natom + 1), &
                atom_position(natom), chunk_atom_begin(nproc), chunk_atom_end(nproc), &
                cursor(nflat_local))
      cached_layout%feature_counts(:) = feature_counts
      cached_layout%feature_displs(:) = feature_displs
      cached_layout%dynamic_counts(:) = ndynamic_per_point*point_counts
      cached_layout%dynamic_displs(:) = ndynamic_per_point*point_displs
      cached_layout%atomic_grid_sizes = 0_int_8
      cached_layout%global_to_feature = 0
      cached_layout%local_feature_counts(:) = local_feature_counts_tmp
      cached_layout%local_feature_offsets(1) = 1
      DO local_row = 1, nflat_local
         cached_layout%local_feature_offsets(local_row + 1) = &
            cached_layout%local_feature_offsets(local_row) + &
            cached_layout%local_feature_counts(local_row)
      END DO
      cursor(:) = cached_layout%local_feature_offsets(1:nflat_local)

      CALL timeset("skala_gpw_layout_atom_sort", phase_handle)
      DO ipt = 1, nflat
         cached_layout%atomic_grid_sizes(global_owner(ipt)) = &
            cached_layout%atomic_grid_sizes(global_owner(ipt)) + 1_int_8
      END DO
      atom_offset(1) = 1
      DO iatom = 1, natom
         atom_offset(iatom + 1) = atom_offset(iatom) + INT(cached_layout%atomic_grid_sizes(iatom))
      END DO
      DO iatom = 1, natom
         atom_position(iatom) = atom_offset(iatom)
      END DO
      max_grid_size = MAXVAL(INT(cached_layout%atomic_grid_sizes))
      CALL build_atom_chunks(cached_layout%atomic_grid_sizes, atom_offset, nproc, &
                             chunk_atom_begin, chunk_atom_end, &
                             cached_layout%chunk_feature_counts, &
                             cached_layout%chunk_feature_displs)
      cached_layout%chunk_grad_counts(:) = ngrad_per_point*cached_layout%chunk_feature_counts
      cached_layout%chunk_grad_displs(:) = ngrad_per_point*cached_layout%chunk_feature_displs
      cached_layout%chunk_atom_begin = chunk_atom_begin(pe_index)
      cached_layout%chunk_atom_end = chunk_atom_end(pe_index)
      cached_layout%chunk_feature_begin = cached_layout%chunk_feature_displs(pe_index) + 1
      cached_layout%chunk_feature_count = cached_layout%chunk_feature_counts(pe_index)
      cached_layout%chunk_natom = cached_layout%chunk_atom_end - &
                                  cached_layout%chunk_atom_begin + 1

      ALLOCATE (cached_layout%grid_coords(3, nflat), cached_layout%grid_weights(nflat), &
                cached_layout%atomic_grid_weights(nflat), &
                cached_layout%coarse_0_atomic_coords(3, natom), &
                cached_layout%atomic_grid_size_bound_shape(0, max_grid_size), &
                cached_layout%atom_coords(3, natom))
      cached_layout%grid_coords = 0.0_dp
      cached_layout%grid_weights = 0.0_dp
      cached_layout%atomic_grid_weights = 0.0_dp
      cached_layout%atomic_grid_size_bound_shape = 0_int_8

      DO iatom = 1, natom
         cached_layout%atom_coords(:, iatom) = particle_set(iatom)%r
         cached_layout%coarse_0_atomic_coords(:, iatom) = atom_coords_pbc(:, iatom)
      END DO

      DO ipt = 1, nflat
         owner = global_owner(ipt)
         row = atom_position(owner)
         atom_position(owner) = atom_position(owner) + 1
         source_global = global_source_points(ipt)
         cached_layout%feature_source_points(row) = source_global
         IF (cached_layout%global_to_feature(source_global) == 0) THEN
            cached_layout%global_to_feature(source_global) = row
         END IF
         static_base = nstatic_per_point*(ipt - 1)
         cached_layout%grid_coords(:, row) = global_static(static_base + 1:static_base + 3)
         cached_layout%grid_weights(row) = global_static(static_base + 4)
         cached_layout%atomic_grid_weights(row) = global_static(static_base + 5)
         source_local = source_global - point_displs(pe_index)
         IF (source_local >= 1 .AND. source_local <= nflat_local) THEN
            feature_local = cursor(source_local)
            cursor(source_local) = cursor(source_local) + 1
            cached_layout%local_feature_rows(feature_local) = row
            cached_layout%local_feature_points(feature_local) = source_local
         END IF
      END DO

      CPASSERT(ALL(cached_layout%global_to_feature > 0))
      CPASSERT(ALL(cached_layout%local_feature_rows > 0))
      CPASSERT(ALL(cached_layout%local_feature_points > 0))
      DO k = bo(1, 3), bo(2, 3)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 1), bo(2, 1)
               local_row = cached_layout%feature_index(i, j, k)
               cached_layout%feature_index(i, j, k) = &
                  cached_layout%local_feature_rows(cached_layout%local_feature_offsets(local_row))
            END DO
         END DO
      END DO
      DO feature_local = 1, nfeature_local
         cached_layout%local_feature_indices(feature_local) = &
            INT(cached_layout%local_feature_rows(feature_local) - 1, KIND=int_8)
      END DO
      CALL timestop(phase_handle)
      CALL timeset("skala_gpw_layout_chunk_routes", phase_handle)
      CALL build_atom_chunk_routes(cached_layout, cached_layout%local_feature_rows, &
                                   pw_grid%para%group)
      CALL build_atom_chunk_layout(cached_layout)
      CALL timestop(phase_handle)

      cached_layout%natom = natom
      cached_layout%nflat = nflat
      cached_layout%nflat_local = nflat_local
      cached_layout%npoint = npoint
      cached_layout%nproc = nproc
      cached_layout%atom_partition = my_atom_partition
      cached_layout%bo = bo
      cached_layout%bounds = pw_grid%bounds
      cached_layout%npts = pw_grid%npts
      cached_layout%dvol = pw_grid%dvol
      cached_layout%dh = pw_grid%dh
      cached_layout%cell_hmat = cell%hmat
      cached_layout%weight_sum = weight_sum
      cached_layout%weight_sumsq = weight_sumsq
      cached_layout%has_weights = has_weights
      CALL timeset("skala_gpw_layout_tensors", phase_handle)
      CALL build_static_layout_tensors(cached_layout)
      CALL timestop(phase_handle)
      cached_layout%active = .TRUE.

      DEALLOCATE (atom_coords_pbc, atom_image_coords, atom_offset, atom_position, &
                  chunk_atom_begin, chunk_atom_end, cursor, feature_counts, feature_displs, &
                  global_owner, global_source_points, global_static, local_feature_counts_tmp, &
                  distances, local_owner, local_source_global, local_source_points, &
                  local_static, partition_weights, point_counts, point_displs, static_counts, &
                  static_displs)

   END SUBROUTINE rebuild_layout_cache

! **************************************************************************************************
!> \brief Build cached Torch tensors for static SKALA inputs.
!> \param cache ...
! **************************************************************************************************
   SUBROUTINE build_static_layout_tensors(cache)
      TYPE(skala_gpw_layout_cache_type), INTENT(INOUT)   :: cache

      CPASSERT(.NOT. cache%static_tensors_active)

      CALL torch_tensor_from_array(cache%grid_coords_t, cache%grid_coords)
      CALL torch_tensor_to_device_leaf(cache%grid_coords_t, .FALSE.)
      CALL torch_tensor_from_array(cache%grid_weights_t, cache%grid_weights)
      CALL torch_tensor_to_device_leaf(cache%grid_weights_t, .FALSE.)
      CALL torch_tensor_from_array(cache%atomic_grid_weights_t, cache%atomic_grid_weights)
      CALL torch_tensor_to_device_leaf(cache%atomic_grid_weights_t, .FALSE.)
      CALL torch_tensor_from_array(cache%atomic_grid_sizes_t, cache%atomic_grid_sizes)
      CALL torch_tensor_to_device_leaf(cache%atomic_grid_sizes_t, .FALSE.)
      CALL torch_tensor_from_array(cache%coarse_0_atomic_coords_t, cache%coarse_0_atomic_coords)
      CALL torch_tensor_to_device_leaf(cache%coarse_0_atomic_coords_t, .FALSE.)
      CALL torch_tensor_from_array(cache%atomic_grid_size_bound_shape_t, &
                                   cache%atomic_grid_size_bound_shape)
      CALL torch_tensor_to_device_leaf(cache%atomic_grid_size_bound_shape_t, .FALSE.)
      CALL torch_tensor_from_array(cache%local_feature_indices_t, cache%local_feature_indices)
      CALL torch_tensor_to_device_leaf(cache%local_feature_indices_t, .FALSE.)

      CALL torch_dict_create(cache%static_inputs)
      CALL torch_dict_insert(cache%static_inputs, "grid_coords", cache%grid_coords_t)
      CALL torch_dict_insert(cache%static_inputs, "grid_weights", cache%grid_weights_t)
      CALL torch_dict_insert(cache%static_inputs, "atomic_grid_weights", &
                             cache%atomic_grid_weights_t)
      CALL torch_dict_insert(cache%static_inputs, "atomic_grid_sizes", &
                             cache%atomic_grid_sizes_t)
      CALL torch_dict_insert(cache%static_inputs, "atomic_grid_size_bound_shape", &
                             cache%atomic_grid_size_bound_shape_t)
      cache%static_tensors_active = .TRUE.

      IF (cache%chunk_feature_count > 0) THEN
         CPASSERT(.NOT. cache%chunk_static_tensors_active)
         CALL torch_tensor_from_array(cache%chunk_grid_coords_t, cache%chunk_grid_coords)
         CALL torch_tensor_to_device_leaf(cache%chunk_grid_coords_t, .FALSE.)
         CALL torch_tensor_from_array(cache%chunk_grid_weights_t, cache%chunk_grid_weights)
         CALL torch_tensor_to_device_leaf(cache%chunk_grid_weights_t, .FALSE.)
         CALL torch_tensor_from_array(cache%chunk_atomic_grid_weights_t, &
                                      cache%chunk_atomic_grid_weights)
         CALL torch_tensor_to_device_leaf(cache%chunk_atomic_grid_weights_t, .FALSE.)
         CALL torch_tensor_from_array(cache%chunk_atomic_grid_sizes_t, &
                                      cache%chunk_atomic_grid_sizes)
         CALL torch_tensor_to_device_leaf(cache%chunk_atomic_grid_sizes_t, .FALSE.)
         CALL torch_tensor_from_array(cache%chunk_coarse_0_atomic_coords_t, &
                                      cache%chunk_coarse_0_atomic_coords)
         CALL torch_tensor_to_device_leaf(cache%chunk_coarse_0_atomic_coords_t, .FALSE.)
         CALL torch_tensor_from_array(cache%chunk_atomic_grid_size_bound_shape_t, &
                                      cache%chunk_atomic_grid_size_bound_shape)
         CALL torch_tensor_to_device_leaf(cache%chunk_atomic_grid_size_bound_shape_t, .FALSE.)
         CALL torch_tensor_from_array(cache%chunk_feature_indices_t, cache%chunk_feature_indices)
         CALL torch_tensor_to_device_leaf(cache%chunk_feature_indices_t, .FALSE.)

         CALL torch_dict_create(cache%chunk_static_inputs)
         CALL torch_dict_insert(cache%chunk_static_inputs, "grid_coords", &
                                cache%chunk_grid_coords_t)
         CALL torch_dict_insert(cache%chunk_static_inputs, "grid_weights", &
                                cache%chunk_grid_weights_t)
         CALL torch_dict_insert(cache%chunk_static_inputs, "atomic_grid_weights", &
                                cache%chunk_atomic_grid_weights_t)
         CALL torch_dict_insert(cache%chunk_static_inputs, "atomic_grid_sizes", &
                                cache%chunk_atomic_grid_sizes_t)
         CALL torch_dict_insert(cache%chunk_static_inputs, "atomic_grid_size_bound_shape", &
                                cache%chunk_atomic_grid_size_bound_shape_t)
         cache%chunk_static_tensors_active = .TRUE.
      END IF

   END SUBROUTINE build_static_layout_tensors

! **************************************************************************************************
!> \brief Copy static cached layout arrays into a feature bundle.
!> \param features ...
!> \param needs_coordinate_array ...
!> \param needs_grid_coordinate_array ...
! **************************************************************************************************
   SUBROUTINE copy_cached_layout(features, needs_coordinate_array, needs_grid_coordinate_array)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      LOGICAL, INTENT(IN)                                :: needs_coordinate_array, &
                                                            needs_grid_coordinate_array

      CPASSERT(cached_layout%active)

      ALLOCATE (features%feature_index(LBOUND(cached_layout%feature_index, 1): &
                                       UBOUND(cached_layout%feature_index, 1), &
                                       LBOUND(cached_layout%feature_index, 2): &
                                       UBOUND(cached_layout%feature_index, 2), &
                                       LBOUND(cached_layout%feature_index, 3): &
                                       UBOUND(cached_layout%feature_index, 3)))
      ALLOCATE (features%grid_weights(cached_layout%nflat))
      ALLOCATE (features%local_feature_counts(cached_layout%nflat_local), &
                features%local_feature_offsets(cached_layout%nflat_local + 1), &
                features%local_feature_rows(SIZE(cached_layout%local_feature_rows)))

      features%feature_index(:, :, :) = cached_layout%feature_index
      features%grid_weights(:) = cached_layout%grid_weights
      features%local_feature_counts(:) = cached_layout%local_feature_counts
      features%local_feature_offsets(:) = cached_layout%local_feature_offsets
      features%local_feature_rows(:) = cached_layout%local_feature_rows
      features%nflat = cached_layout%nflat
      features%nflat_local = cached_layout%nflat_local
      features%chunk_feature_count = cached_layout%chunk_feature_count
      features%atom_partition = cached_layout%atom_partition
      ALLOCATE (features%atomic_grid_sizes(cached_layout%natom))
      features%atomic_grid_sizes(:) = cached_layout%atomic_grid_sizes
      IF (needs_grid_coordinate_array) THEN
         ALLOCATE (features%grid_coords(3, cached_layout%nflat))
         ALLOCATE (features%atomic_grid_weights(cached_layout%nflat))
         features%grid_coords(:, :) = cached_layout%grid_coords
         features%atomic_grid_weights(:) = cached_layout%atomic_grid_weights
      END IF
      ALLOCATE (features%chunk_grad_counts(cached_layout%nproc), &
                features%chunk_grad_displs(cached_layout%nproc), &
                features%route_grad_return_recv_counts(cached_layout%nproc), &
                features%route_grad_return_recv_displs(cached_layout%nproc), &
                features%route_grad_return_send_counts(cached_layout%nproc), &
                features%route_grad_return_send_displs(cached_layout%nproc), &
                features%route_point_recv_counts(cached_layout%nproc), &
                features%route_point_recv_displs(cached_layout%nproc), &
                features%route_point_send_counts(cached_layout%nproc), &
                features%route_point_send_displs(cached_layout%nproc), &
                features%route_send_local_rows(SIZE(cached_layout%route_send_local_rows)))
      features%chunk_grad_counts(:) = cached_layout%chunk_grad_counts
      features%chunk_grad_displs(:) = cached_layout%chunk_grad_displs
      features%route_grad_return_recv_counts(:) = cached_layout%route_grad_return_recv_counts
      features%route_grad_return_recv_displs(:) = cached_layout%route_grad_return_recv_displs
      features%route_grad_return_send_counts(:) = cached_layout%route_grad_return_send_counts
      features%route_grad_return_send_displs(:) = cached_layout%route_grad_return_send_displs
      features%route_point_recv_counts(:) = cached_layout%route_point_recv_counts
      features%route_point_recv_displs(:) = cached_layout%route_point_recv_displs
      features%route_point_send_counts(:) = cached_layout%route_point_send_counts
      features%route_point_send_displs(:) = cached_layout%route_point_send_displs
      features%route_send_local_rows(:) = cached_layout%route_send_local_rows
      IF (needs_coordinate_array) THEN
         ALLOCATE (features%coarse_0_atomic_coords(3, cached_layout%natom))
         features%coarse_0_atomic_coords(:, :) = cached_layout%coarse_0_atomic_coords
      END IF

   END SUBROUTINE copy_cached_layout

! **************************************************************************************************
!> \brief Split the atom-ordered feature rows into contiguous atom chunks.
!> \param atomic_grid_sizes ...
!> \param atom_offset ...
!> \param nproc ...
!> \param chunk_atom_begin ...
!> \param chunk_atom_end ...
!> \param chunk_feature_counts ...
!> \param chunk_feature_displs ...
! **************************************************************************************************
   SUBROUTINE build_atom_chunks(atomic_grid_sizes, atom_offset, nproc, chunk_atom_begin, &
                                chunk_atom_end, chunk_feature_counts, chunk_feature_displs)
      INTEGER(KIND=int_8), DIMENSION(:), INTENT(IN)      :: atomic_grid_sizes
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_offset
      INTEGER, INTENT(IN)                                :: nproc
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: chunk_atom_begin, chunk_atom_end, &
                                                            chunk_feature_counts, &
                                                            chunk_feature_displs

      INTEGER :: best_limit, count, displ, end_atom, lower_limit, max_end_atom, midpoint, natom, &
         next_atom, next_count, pe, ranks_left, target_chunks, total_count, upper_limit

      natom = SIZE(atomic_grid_sizes)
      chunk_atom_begin = natom + 1
      chunk_atom_end = natom
      chunk_feature_counts = 0
      chunk_feature_displs = 0
      IF (natom == 0) RETURN

      target_chunks = MIN(nproc, natom)
      total_count = atom_offset(natom + 1) - 1
      lower_limit = MAXVAL(INT(atomic_grid_sizes))
      lower_limit = MAX(lower_limit, (total_count + target_chunks - 1)/target_chunks)
      upper_limit = total_count
      best_limit = upper_limit
      DO WHILE (lower_limit <= upper_limit)
         midpoint = (lower_limit + upper_limit)/2
         IF (atom_chunks_fit_limit(atomic_grid_sizes, midpoint, target_chunks)) THEN
            best_limit = midpoint
            upper_limit = midpoint - 1
         ELSE
            lower_limit = midpoint + 1
         END IF
      END DO

      displ = 0
      next_atom = 1
      DO pe = 1, nproc
         chunk_feature_displs(pe) = displ
         IF (pe > target_chunks .OR. next_atom > natom) CYCLE

         ranks_left = target_chunks - pe + 1
         chunk_atom_begin(pe) = next_atom
         max_end_atom = natom - ranks_left + 1
         end_atom = next_atom
         count = INT(atomic_grid_sizes(end_atom))
         DO WHILE (end_atom < max_end_atom)
            next_count = count + INT(atomic_grid_sizes(end_atom + 1))
            IF (next_count > best_limit) EXIT
            end_atom = end_atom + 1
            count = next_count
         END DO

         chunk_atom_end(pe) = end_atom
         chunk_feature_counts(pe) = atom_offset(end_atom + 1) - atom_offset(next_atom)
         displ = displ + chunk_feature_counts(pe)
         next_atom = end_atom + 1
      END DO

      CPASSERT(displ == atom_offset(natom + 1) - 1)

   END SUBROUTINE build_atom_chunks

! **************************************************************************************************
!> \brief Check if contiguous atom chunks can stay below a feature-count limit.
!> \param atomic_grid_sizes ...
!> \param limit ...
!> \param nchunks ...
!> \return ...
! **************************************************************************************************
   FUNCTION atom_chunks_fit_limit(atomic_grid_sizes, limit, nchunks) RESULT(fits)
      INTEGER(KIND=int_8), DIMENSION(:), INTENT(IN)      :: atomic_grid_sizes
      INTEGER, INTENT(IN)                                :: limit, nchunks
      LOGICAL                                            :: fits

      INTEGER                                            :: atom_count, chunk_count, iatom, &
                                                            used_chunks

      fits = .FALSE.
      IF (SIZE(atomic_grid_sizes) == 0) THEN
         fits = .TRUE.
         RETURN
      END IF

      used_chunks = 1
      chunk_count = 0
      DO iatom = 1, SIZE(atomic_grid_sizes)
         atom_count = INT(atomic_grid_sizes(iatom))
         IF (atom_count > limit) RETURN
         IF (chunk_count + atom_count > limit) THEN
            used_chunks = used_chunks + 1
            chunk_count = atom_count
         ELSE
            chunk_count = chunk_count + atom_count
         END IF
      END DO
      fits = used_chunks <= nchunks

   END FUNCTION atom_chunks_fit_limit

! **************************************************************************************************
!> \brief Return the MPI rank owning an atom-ordered feature row.
!> \param row ...
!> \param counts ...
!> \param displs ...
!> \return ...
! **************************************************************************************************
   FUNCTION feature_row_chunk_owner(row, counts, displs) RESULT(owner)
      INTEGER, INTENT(IN)                                :: row
      INTEGER, DIMENSION(:), INTENT(IN)                  :: counts, displs
      INTEGER                                            :: owner

      INTEGER                                            :: pe

      owner = 0
      DO pe = 1, SIZE(counts)
         IF (row > displs(pe) .AND. row <= displs(pe) + counts(pe)) THEN
            owner = pe
            RETURN
         END IF
      END DO

   END FUNCTION feature_row_chunk_owner

! **************************************************************************************************
!> \brief Build zero-based displacement arrays from per-rank counts.
!> \param counts ...
!> \param displs ...
! **************************************************************************************************
   SUBROUTINE counts_to_displs(counts, displs)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: counts
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: displs

      INTEGER                                            :: pe

      displs(1) = 0
      DO pe = 2, SIZE(counts)
         displs(pe) = displs(pe - 1) + counts(pe - 1)
      END DO

   END SUBROUTINE counts_to_displs

! **************************************************************************************************
!> \brief Precompute all-to-all routing between local grid rows and atom chunks.
!> \param cache ...
!> \param local_to_global ...
!> \param group ...
! **************************************************************************************************
   SUBROUTINE build_atom_chunk_routes(cache, local_to_global, group)
      TYPE(skala_gpw_layout_cache_type), INTENT(INOUT)   :: cache
      INTEGER, DIMENSION(:), INTENT(IN)                  :: local_to_global

      CLASS(mp_comm_type), INTENT(IN)                    :: group

      INTEGER                                            :: chunk_row, dest, local_feature, point_pos, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cursor, recv_meta, send_meta

      ALLOCATE (cache%route_local_dest(SIZE(local_to_global)), &
                cache%route_send_local_rows(SIZE(local_to_global)), &
                cache%chunk_return_positions(cache%chunk_feature_count), &
                cursor(SIZE(cache%route_point_send_counts)))
      cache%route_point_send_counts = 0
      cache%route_send_local_rows = 0
      cache%chunk_return_positions = 0
      DO local_feature = 1, SIZE(local_to_global)
         dest = feature_row_chunk_owner(local_to_global(local_feature), &
                                        cache%chunk_feature_counts, &
                                        cache%chunk_feature_displs)
         CPASSERT(dest > 0)
         cache%route_local_dest(local_feature) = dest
         cache%route_point_send_counts(dest) = cache%route_point_send_counts(dest) + 1
      END DO
      CALL counts_to_displs(cache%route_point_send_counts, cache%route_point_send_displs)
      cursor(:) = cache%route_point_send_displs + 1
      DO local_feature = 1, SIZE(local_to_global)
         dest = cache%route_local_dest(local_feature)
         point_pos = cursor(dest)
         cursor(dest) = cursor(dest) + 1
         cache%route_send_local_rows(point_pos) = cache%local_feature_points(local_feature)
      END DO
      CALL group%alltoall(cache%route_point_send_counts, cache%route_point_recv_counts, 1)
      CALL counts_to_displs(cache%route_point_recv_counts, cache%route_point_recv_displs)

      ALLOCATE (send_meta(SIZE(local_to_global)), recv_meta(cache%chunk_feature_count))
      cursor(:) = cache%route_point_send_displs + 1
      DO local_feature = 1, SIZE(local_to_global)
         dest = cache%route_local_dest(local_feature)
         point_pos = cursor(dest)
         cursor(dest) = cursor(dest) + 1
         send_meta(point_pos) = local_to_global(local_feature)
      END DO
      CALL group%alltoall(send_meta, cache%route_point_send_counts, &
                          cache%route_point_send_displs, recv_meta, &
                          cache%route_point_recv_counts, &
                          cache%route_point_recv_displs)
      DO point_pos = 1, cache%chunk_feature_count
         row = recv_meta(point_pos)
         chunk_row = row - cache%chunk_feature_begin + 1
         CPASSERT(chunk_row >= 1 .AND. chunk_row <= cache%chunk_feature_count)
         cache%chunk_return_positions(chunk_row) = point_pos
      END DO

      cache%route_grad_return_send_counts(:) = ngrad_per_point*cache%route_point_recv_counts
      cache%route_grad_return_send_displs(:) = ngrad_per_point*cache%route_point_recv_displs
      cache%route_grad_return_recv_counts(:) = ngrad_per_point*cache%route_point_send_counts
      cache%route_grad_return_recv_displs(:) = ngrad_per_point*cache%route_point_send_displs

      CPASSERT(SUM(cache%route_point_send_counts) == SIZE(local_to_global))
      CPASSERT(SUM(cache%route_point_recv_counts) == cache%chunk_feature_count)
      CPASSERT(ALL(cache%route_send_local_rows > 0))
      CPASSERT(ALL(cache%chunk_return_positions > 0))

      DEALLOCATE (cursor, recv_meta, send_meta)

   END SUBROUTINE build_atom_chunk_routes

! **************************************************************************************************
!> \brief Materialize the current rank's atom chunk static layout.
!> \param cache ...
! **************************************************************************************************
   SUBROUTINE build_atom_chunk_layout(cache)
      TYPE(skala_gpw_layout_cache_type), INTENT(INOUT)   :: cache

      INTEGER                                            :: irow, max_grid_size, row_begin, row_end

      IF (cache%chunk_feature_count <= 0 .OR. cache%chunk_natom <= 0) RETURN

      row_begin = cache%chunk_feature_begin
      row_end = row_begin + cache%chunk_feature_count - 1
      ALLOCATE (cache%chunk_grid_coords(3, cache%chunk_feature_count), &
                cache%chunk_grid_weights(cache%chunk_feature_count), &
                cache%chunk_atomic_grid_weights(cache%chunk_feature_count), &
                cache%chunk_atomic_grid_sizes(cache%chunk_natom), &
                cache%chunk_coarse_0_atomic_coords(3, cache%chunk_natom), &
                cache%chunk_feature_indices(cache%chunk_feature_count))
      cache%chunk_grid_coords(:, :) = cache%grid_coords(:, row_begin:row_end)
      cache%chunk_grid_weights(:) = cache%grid_weights(row_begin:row_end)
      cache%chunk_atomic_grid_weights(:) = cache%atomic_grid_weights(row_begin:row_end)
      cache%chunk_atomic_grid_sizes(:) = &
         cache%atomic_grid_sizes(cache%chunk_atom_begin:cache%chunk_atom_end)
      cache%chunk_coarse_0_atomic_coords(:, :) = &
         cache%coarse_0_atomic_coords(:, cache%chunk_atom_begin:cache%chunk_atom_end)

      max_grid_size = MAXVAL(INT(cache%chunk_atomic_grid_sizes))
      ALLOCATE (cache%chunk_atomic_grid_size_bound_shape(0, max_grid_size))
      cache%chunk_atomic_grid_size_bound_shape = 0_int_8
      DO irow = 1, cache%chunk_feature_count
         cache%chunk_feature_indices(irow) = INT(irow - 1, KIND=int_8)
      END DO

   END SUBROUTINE build_atom_chunk_layout

! **************************************************************************************************
!> \brief Send local dynamic feature rows to their atom-chunk owner ranks.
!> \param features ...
!> \param local_dynamic ...
!> \param group ...
!> \param collapse_spin_dynamics ...
! **************************************************************************************************
   SUBROUTINE route_atom_chunk_dynamics(features, local_dynamic, group, collapse_spin_dynamics)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: local_dynamic

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      LOGICAL, INTENT(IN)                                :: collapse_spin_dynamics

      INTEGER                                            :: chunk_row, dest, dyn_base, local_feature, local_row, &
                                                            ndynamic_route_per_point, nrecv, nsend, &
                                                            point_pos, src_base
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cursor, recv_counts, recv_displs, &
                                                            send_counts, send_displs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: recv_dynamic, send_dynamic

      nsend = SIZE(cached_layout%route_local_dest)
      nrecv = SUM(cached_layout%route_point_recv_counts)
      CPASSERT(nsend == SIZE(cached_layout%local_feature_rows))
      CPASSERT(nrecv == cached_layout%chunk_feature_count)
      ndynamic_route_per_point = ndynamic_per_point
      IF (collapse_spin_dynamics) ndynamic_route_per_point = nrks_dynamic_per_point

      ALLOCATE (send_dynamic(MAX(1, ndynamic_route_per_point*nsend)), &
                recv_dynamic(MAX(1, ndynamic_route_per_point*nrecv)), &
                cursor(cached_layout%nproc), send_counts(cached_layout%nproc), &
                send_displs(cached_layout%nproc), recv_counts(cached_layout%nproc), &
                recv_displs(cached_layout%nproc))
      send_counts(:) = ndynamic_route_per_point*cached_layout%route_point_send_counts
      send_displs(:) = ndynamic_route_per_point*cached_layout%route_point_send_displs
      recv_counts(:) = ndynamic_route_per_point*cached_layout%route_point_recv_counts
      recv_displs(:) = ndynamic_route_per_point*cached_layout%route_point_recv_displs
      cursor(:) = cached_layout%route_point_send_displs + 1
      DO local_feature = 1, nsend
         dest = cached_layout%route_local_dest(local_feature)
         point_pos = cursor(dest)
         cursor(dest) = cursor(dest) + 1
         dyn_base = ndynamic_route_per_point*(point_pos - 1)
         local_row = cached_layout%local_feature_points(local_feature)
         src_base = ndynamic_route_per_point*(local_row - 1)
         send_dynamic(dyn_base + 1:dyn_base + ndynamic_route_per_point) = &
            local_dynamic(src_base + 1:src_base + ndynamic_route_per_point)
      END DO

      CALL group%alltoall(send_dynamic, send_counts, send_displs, recv_dynamic, recv_counts, &
                          recv_displs)

      features%uses_collapsed_rks_dynamic = collapse_spin_dynamics
      IF (cached_layout%chunk_feature_count > 0) THEN
         IF (collapse_spin_dynamics) THEN
            ALLOCATE (features%chunk_density(cached_layout%chunk_feature_count, 1), &
                      features%chunk_grad(cached_layout%chunk_feature_count, 3, 1), &
                      features%chunk_kin(cached_layout%chunk_feature_count, 1), &
                      features%chunk_return_positions(cached_layout%chunk_feature_count))
         ELSE
            ALLOCATE (features%chunk_density(cached_layout%chunk_feature_count, 2), &
                      features%chunk_grad(cached_layout%chunk_feature_count, 3, 2), &
                      features%chunk_kin(cached_layout%chunk_feature_count, 2), &
                      features%chunk_return_positions(cached_layout%chunk_feature_count))
         END IF
         features%chunk_return_positions(:) = cached_layout%chunk_return_positions

         DO chunk_row = 1, cached_layout%chunk_feature_count
            point_pos = cached_layout%chunk_return_positions(chunk_row)
            CPASSERT(point_pos >= 1 .AND. point_pos <= cached_layout%chunk_feature_count)
            dyn_base = ndynamic_route_per_point*(point_pos - 1)
            IF (collapse_spin_dynamics) THEN
               features%chunk_density(chunk_row, 1) = recv_dynamic(dyn_base + 1)
               features%chunk_grad(chunk_row, 1, 1) = recv_dynamic(dyn_base + 2)
               features%chunk_grad(chunk_row, 2, 1) = recv_dynamic(dyn_base + 3)
               features%chunk_grad(chunk_row, 3, 1) = recv_dynamic(dyn_base + 4)
               features%chunk_kin(chunk_row, 1) = recv_dynamic(dyn_base + 5)
            ELSE
               features%chunk_density(chunk_row, :) = recv_dynamic(dyn_base + 1:dyn_base + 2)
               features%chunk_grad(chunk_row, 1, 1) = recv_dynamic(dyn_base + 3)
               features%chunk_grad(chunk_row, 2, 1) = recv_dynamic(dyn_base + 4)
               features%chunk_grad(chunk_row, 3, 1) = recv_dynamic(dyn_base + 5)
               features%chunk_grad(chunk_row, 1, 2) = recv_dynamic(dyn_base + 6)
               features%chunk_grad(chunk_row, 2, 2) = recv_dynamic(dyn_base + 7)
               features%chunk_grad(chunk_row, 3, 2) = recv_dynamic(dyn_base + 8)
               features%chunk_kin(chunk_row, :) = recv_dynamic(dyn_base + 9:dyn_base + 10)
            END IF
         END DO
         CPASSERT(ALL(features%chunk_return_positions > 0))
      END IF

      DEALLOCATE (cursor, recv_counts, recv_displs, recv_dynamic, send_counts, send_displs, &
                  send_dynamic)

   END SUBROUTINE route_atom_chunk_dynamics

! **************************************************************************************************
!> \brief Extract the current rank's atom chunk from the global dynamic feature arrays.
!> \param features ...
! **************************************************************************************************
   SUBROUTINE extract_atom_chunk_dynamics(features)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features

      INTEGER                                            :: row_begin, row_end

      CPASSERT(cached_layout%chunk_feature_count > 0)
      row_begin = cached_layout%chunk_feature_begin
      row_end = row_begin + cached_layout%chunk_feature_count - 1
      ALLOCATE (features%chunk_density(cached_layout%chunk_feature_count, 2), &
                features%chunk_grad(cached_layout%chunk_feature_count, 3, 2), &
                features%chunk_kin(cached_layout%chunk_feature_count, 2))
      features%chunk_density(:, :) = features%density(row_begin:row_end, :)
      features%chunk_grad(:, :, :) = features%grad(row_begin:row_end, :, :)
      features%chunk_kin(:, :) = features%kin(row_begin:row_end, :)

   END SUBROUTINE extract_atom_chunk_dynamics

! **************************************************************************************************
!> \brief Compute a local signature for optional integration weights.
!> \param weights ...
!> \param has_weights ...
!> \param weight_sum ...
!> \param weight_sumsq ...
! **************************************************************************************************
   SUBROUTINE weights_signature(weights, has_weights, weight_sum, weight_sumsq)
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: weights
      LOGICAL, INTENT(OUT)                               :: has_weights
      REAL(KIND=dp), INTENT(OUT)                         :: weight_sum, weight_sumsq

      has_weights = .FALSE.
      weight_sum = 0.0_dp
      weight_sumsq = 0.0_dp
      IF (PRESENT(weights)) THEN
         IF (ASSOCIATED(weights)) THEN
            has_weights = .TRUE.
            weight_sum = SUM(weights%array)
            weight_sumsq = SUM(weights%array*weights%array)
         END IF
      END IF

   END SUBROUTINE weights_signature

! **************************************************************************************************
!> \brief Release cached layout arrays.
!> \param cache ...
! **************************************************************************************************
   SUBROUTINE release_layout_cache(cache)
      TYPE(skala_gpw_layout_cache_type), INTENT(INOUT)   :: cache

      IF (cache%inputs_active) THEN
         CALL torch_dict_release(cache%inputs)
         cache%inputs_active = .FALSE.
      END IF

      IF (cache%chunk_inputs_active) THEN
         CALL torch_dict_release(cache%chunk_inputs)
         cache%chunk_inputs_active = .FALSE.
      END IF

      IF (cache%dynamic_tensors_active) THEN
         CALL torch_tensor_release(cache%density_t)
         CALL torch_tensor_release(cache%grad_t)
         CALL torch_tensor_release(cache%kin_t)
         cache%dynamic_tensors_active = .FALSE.
      END IF

      IF (cache%chunk_dynamic_tensors_active) THEN
         IF (cache%chunk_dynamic_input_views_active) THEN
            CALL torch_tensor_release(cache%chunk_density_input_t)
            CALL torch_tensor_release(cache%chunk_grad_input_t)
            CALL torch_tensor_release(cache%chunk_kin_input_t)
            cache%chunk_dynamic_input_views_active = .FALSE.
         END IF
         CALL torch_tensor_release(cache%chunk_density_t)
         CALL torch_tensor_release(cache%chunk_grad_t)
         CALL torch_tensor_release(cache%chunk_kin_t)
         cache%chunk_dynamic_tensors_active = .FALSE.
      END IF

      IF (cache%static_tensors_active) THEN
         CALL torch_tensor_release(cache%grid_coords_t)
         CALL torch_tensor_release(cache%grid_weights_t)
         CALL torch_tensor_release(cache%atomic_grid_weights_t)
         CALL torch_tensor_release(cache%atomic_grid_sizes_t)
         CALL torch_tensor_release(cache%coarse_0_atomic_coords_t)
         CALL torch_tensor_release(cache%atomic_grid_size_bound_shape_t)
         CALL torch_tensor_release(cache%local_feature_indices_t)
         CALL torch_dict_release(cache%static_inputs)
         cache%static_tensors_active = .FALSE.
      END IF

      IF (cache%chunk_static_tensors_active) THEN
         CALL torch_tensor_release(cache%chunk_grid_coords_t)
         CALL torch_tensor_release(cache%chunk_grid_weights_t)
         CALL torch_tensor_release(cache%chunk_atomic_grid_weights_t)
         CALL torch_tensor_release(cache%chunk_atomic_grid_sizes_t)
         CALL torch_tensor_release(cache%chunk_coarse_0_atomic_coords_t)
         CALL torch_tensor_release(cache%chunk_atomic_grid_size_bound_shape_t)
         CALL torch_tensor_release(cache%chunk_feature_indices_t)
         CALL torch_dict_release(cache%chunk_static_inputs)
         cache%chunk_static_tensors_active = .FALSE.
      END IF

      IF (ALLOCATED(cache%chunk_feature_counts)) DEALLOCATE (cache%chunk_feature_counts)
      IF (ALLOCATED(cache%chunk_feature_displs)) DEALLOCATE (cache%chunk_feature_displs)
      IF (ALLOCATED(cache%chunk_grad_counts)) DEALLOCATE (cache%chunk_grad_counts)
      IF (ALLOCATED(cache%chunk_grad_displs)) DEALLOCATE (cache%chunk_grad_displs)
      IF (ALLOCATED(cache%route_grad_return_recv_counts)) THEN
         DEALLOCATE (cache%route_grad_return_recv_counts)
      END IF
      IF (ALLOCATED(cache%route_grad_return_recv_displs)) THEN
         DEALLOCATE (cache%route_grad_return_recv_displs)
      END IF
      IF (ALLOCATED(cache%route_grad_return_send_counts)) THEN
         DEALLOCATE (cache%route_grad_return_send_counts)
      END IF
      IF (ALLOCATED(cache%route_grad_return_send_displs)) THEN
         DEALLOCATE (cache%route_grad_return_send_displs)
      END IF
      IF (ALLOCATED(cache%route_local_dest)) DEALLOCATE (cache%route_local_dest)
      IF (ALLOCATED(cache%chunk_return_positions)) DEALLOCATE (cache%chunk_return_positions)
      IF (ALLOCATED(cache%route_point_recv_counts)) DEALLOCATE (cache%route_point_recv_counts)
      IF (ALLOCATED(cache%route_point_recv_displs)) DEALLOCATE (cache%route_point_recv_displs)
      IF (ALLOCATED(cache%route_point_send_counts)) DEALLOCATE (cache%route_point_send_counts)
      IF (ALLOCATED(cache%route_point_send_displs)) DEALLOCATE (cache%route_point_send_displs)
      IF (ALLOCATED(cache%route_send_local_rows)) DEALLOCATE (cache%route_send_local_rows)
      IF (ALLOCATED(cache%dynamic_counts)) DEALLOCATE (cache%dynamic_counts)
      IF (ALLOCATED(cache%dynamic_displs)) DEALLOCATE (cache%dynamic_displs)
      IF (ALLOCATED(cache%feature_counts)) DEALLOCATE (cache%feature_counts)
      IF (ALLOCATED(cache%feature_displs)) DEALLOCATE (cache%feature_displs)
      IF (ALLOCATED(cache%feature_source_points)) DEALLOCATE (cache%feature_source_points)
      IF (ALLOCATED(cache%global_to_feature)) DEALLOCATE (cache%global_to_feature)
      IF (ALLOCATED(cache%feature_index)) DEALLOCATE (cache%feature_index)
      IF (ALLOCATED(cache%atomic_grid_sizes)) DEALLOCATE (cache%atomic_grid_sizes)
      IF (ALLOCATED(cache%chunk_atomic_grid_sizes)) DEALLOCATE (cache%chunk_atomic_grid_sizes)
      IF (ALLOCATED(cache%chunk_feature_indices)) DEALLOCATE (cache%chunk_feature_indices)
      IF (ALLOCATED(cache%local_feature_counts)) DEALLOCATE (cache%local_feature_counts)
      IF (ALLOCATED(cache%local_feature_indices)) DEALLOCATE (cache%local_feature_indices)
      IF (ALLOCATED(cache%local_feature_offsets)) DEALLOCATE (cache%local_feature_offsets)
      IF (ALLOCATED(cache%local_feature_points)) DEALLOCATE (cache%local_feature_points)
      IF (ALLOCATED(cache%local_feature_rows)) DEALLOCATE (cache%local_feature_rows)
      IF (ALLOCATED(cache%atomic_grid_size_bound_shape)) THEN
         DEALLOCATE (cache%atomic_grid_size_bound_shape)
      END IF
      IF (ALLOCATED(cache%chunk_atomic_grid_size_bound_shape)) THEN
         DEALLOCATE (cache%chunk_atomic_grid_size_bound_shape)
      END IF
      IF (ALLOCATED(cache%atomic_grid_weights)) DEALLOCATE (cache%atomic_grid_weights)
      IF (ALLOCATED(cache%chunk_atomic_grid_weights)) DEALLOCATE (cache%chunk_atomic_grid_weights)
      IF (ALLOCATED(cache%chunk_grid_weights)) DEALLOCATE (cache%chunk_grid_weights)
      IF (ALLOCATED(cache%grid_weights)) DEALLOCATE (cache%grid_weights)
      IF (ALLOCATED(cache%atom_coords)) DEALLOCATE (cache%atom_coords)
      IF (ALLOCATED(cache%chunk_coarse_0_atomic_coords)) THEN
         DEALLOCATE (cache%chunk_coarse_0_atomic_coords)
      END IF
      IF (ALLOCATED(cache%coarse_0_atomic_coords)) DEALLOCATE (cache%coarse_0_atomic_coords)
      IF (ALLOCATED(cache%chunk_grid_coords)) DEALLOCATE (cache%chunk_grid_coords)
      IF (ALLOCATED(cache%grid_coords)) DEALLOCATE (cache%grid_coords)

      cache%chunk_atom_begin = 1
      cache%chunk_atom_end = 0
      cache%chunk_feature_begin = 1
      cache%chunk_feature_count = 0
      cache%chunk_natom = 0
      cache%natom = 0
      cache%nflat = 0
      cache%nflat_local = 0
      cache%npoint = 0
      cache%nproc = 0
      cache%atom_partition = skala_gpw_atom_partition_hard
      cache%bo = 0
      cache%bounds = 0
      cache%npts = 0
      cache%dvol = 0.0_dp
      cache%weight_sum = 0.0_dp
      cache%weight_sumsq = 0.0_dp
      cache%cell_hmat = 0.0_dp
      cache%dh = 0.0_dp
      cache%active = .FALSE.
      cache%has_weights = .FALSE.
      cache%chunk_dynamic_tensors_active = .FALSE.
      cache%chunk_dynamic_input_views_active = .FALSE.
      cache%chunk_inputs_active = .FALSE.
      cache%chunk_inputs_use_collapsed_rks = .FALSE.
      cache%chunk_static_tensors_active = .FALSE.
      cache%dynamic_tensors_active = .FALSE.
      cache%inputs_active = .FALSE.
      cache%static_tensors_active = .FALSE.

   END SUBROUTINE release_layout_cache

! **************************************************************************************************
!> \brief Release Torch objects and backing arrays owned by a feature bundle.
!> \param features ...
! **************************************************************************************************
   SUBROUTINE skala_gpw_feature_release(features)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features

      IF (features%active) THEN
         IF (features%owns_dynamic_tensors) THEN
            IF (features%uses_collapsed_rks_dynamic) THEN
               CALL torch_tensor_release(features%density_input_t)
               CALL torch_tensor_release(features%grad_input_t)
               CALL torch_tensor_release(features%kin_input_t)
            END IF
            CALL torch_tensor_release(features%density_t)
            CALL torch_tensor_release(features%grad_t)
            CALL torch_tensor_release(features%kin_t)
         END IF
         IF (features%owns_static_tensors) THEN
            CALL torch_tensor_release(features%grid_coords_t)
            CALL torch_tensor_release(features%grid_weights_t)
            CALL torch_tensor_release(features%atomic_grid_weights_t)
            CALL torch_tensor_release(features%atomic_grid_sizes_t)
            CALL torch_tensor_release(features%atomic_grid_size_bound_shape_t)
         END IF
         IF (features%owns_grid_coordinate_tensor) THEN
            CALL torch_tensor_release(features%grid_coords_t)
         END IF
         IF (features%owns_weight_tensors) THEN
            CALL torch_tensor_release(features%grid_weights_t)
            CALL torch_tensor_release(features%atomic_grid_weights_t)
         END IF
         IF (features%owns_static_tensors .OR. features%owns_coordinate_tensor) THEN
            CALL torch_tensor_release(features%coarse_0_atomic_coords_t)
         END IF
         IF (features%owns_inputs) CALL torch_dict_release(features%inputs)
         features%active = .FALSE.
         features%owns_coordinate_tensor = .FALSE.
         features%owns_grid_coordinate_tensor = .FALSE.
         features%owns_weight_tensors = .FALSE.
         features%owns_dynamic_tensors = .TRUE.
         features%owns_inputs = .TRUE.
         features%owns_static_tensors = .TRUE.
         features%uses_atom_chunk_routing = .FALSE.
         features%uses_atom_chunks = .FALSE.
         features%uses_collapsed_rks_dynamic = .FALSE.
      END IF

      IF (ALLOCATED(features%chunk_density)) DEALLOCATE (features%chunk_density)
      IF (ALLOCATED(features%chunk_grad)) DEALLOCATE (features%chunk_grad)
      IF (ALLOCATED(features%chunk_kin)) DEALLOCATE (features%chunk_kin)
      IF (ALLOCATED(features%density)) DEALLOCATE (features%density)
      IF (ALLOCATED(features%grad)) DEALLOCATE (features%grad)
      IF (ALLOCATED(features%kin)) DEALLOCATE (features%kin)
      IF (ALLOCATED(features%chunk_grad_counts)) DEALLOCATE (features%chunk_grad_counts)
      IF (ALLOCATED(features%chunk_grad_displs)) DEALLOCATE (features%chunk_grad_displs)
      IF (ALLOCATED(features%chunk_return_positions)) DEALLOCATE (features%chunk_return_positions)
      IF (ALLOCATED(features%route_grad_return_recv_counts)) THEN
         DEALLOCATE (features%route_grad_return_recv_counts)
      END IF
      IF (ALLOCATED(features%route_grad_return_recv_displs)) THEN
         DEALLOCATE (features%route_grad_return_recv_displs)
      END IF
      IF (ALLOCATED(features%route_grad_return_send_counts)) THEN
         DEALLOCATE (features%route_grad_return_send_counts)
      END IF
      IF (ALLOCATED(features%route_grad_return_send_displs)) THEN
         DEALLOCATE (features%route_grad_return_send_displs)
      END IF
      IF (ALLOCATED(features%route_point_recv_counts)) THEN
         DEALLOCATE (features%route_point_recv_counts)
      END IF
      IF (ALLOCATED(features%route_point_recv_displs)) THEN
         DEALLOCATE (features%route_point_recv_displs)
      END IF
      IF (ALLOCATED(features%route_point_send_counts)) THEN
         DEALLOCATE (features%route_point_send_counts)
      END IF
      IF (ALLOCATED(features%route_point_send_displs)) THEN
         DEALLOCATE (features%route_point_send_displs)
      END IF
      IF (ALLOCATED(features%route_send_local_rows)) DEALLOCATE (features%route_send_local_rows)
      IF (ALLOCATED(features%feature_index)) DEALLOCATE (features%feature_index)
      IF (ALLOCATED(features%local_feature_counts)) DEALLOCATE (features%local_feature_counts)
      IF (ALLOCATED(features%local_feature_offsets)) DEALLOCATE (features%local_feature_offsets)
      IF (ALLOCATED(features%local_feature_rows)) DEALLOCATE (features%local_feature_rows)
      IF (ALLOCATED(features%grid_coords)) DEALLOCATE (features%grid_coords)
      IF (ALLOCATED(features%grid_weights)) DEALLOCATE (features%grid_weights)
      IF (ALLOCATED(features%atomic_grid_weights)) DEALLOCATE (features%atomic_grid_weights)
      IF (ALLOCATED(features%atomic_grid_sizes)) DEALLOCATE (features%atomic_grid_sizes)
      IF (ALLOCATED(features%coarse_0_atomic_coords)) DEALLOCATE (features%coarse_0_atomic_coords)
      IF (ALLOCATED(features%atomic_grid_size_bound_shape)) THEN
         DEALLOCATE (features%atomic_grid_size_bound_shape)
      END IF
      features%chunk_feature_count = 0
      features%nflat = 0
      features%nflat_local = 0
      features%atom_partition = skala_gpw_atom_partition_hard
      features%uses_atom_chunk_routing = .FALSE.
      features%uses_collapsed_rks_dynamic = .FALSE.

   END SUBROUTINE skala_gpw_feature_release

! **************************************************************************************************
!> \brief Return how many atom-contiguous subchunks the cached rank chunk needs.
!> \param max_rows ...
!> \return ...
! **************************************************************************************************
   FUNCTION skala_gpw_atom_subchunk_count(max_rows) RESULT(nsubchunks)
      INTEGER, INTENT(IN)                                :: max_rows
      INTEGER                                            :: nsubchunks

      INTEGER                                            :: atom_rows, iatom, rows

      nsubchunks = 0
      IF (.NOT. cached_layout%active) RETURN
      IF (cached_layout%chunk_natom <= 0) RETURN
      IF (max_rows <= 0) THEN
         nsubchunks = 1
         RETURN
      END IF

      rows = 0
      DO iatom = 1, cached_layout%chunk_natom
         atom_rows = INT(cached_layout%chunk_atomic_grid_sizes(iatom))
         IF (rows > 0 .AND. rows + atom_rows > max_rows) THEN
            nsubchunks = nsubchunks + 1
            rows = 0
         END IF
         rows = rows + atom_rows
      END DO
      IF (rows > 0) nsubchunks = nsubchunks + 1
      nsubchunks = MAX(1, nsubchunks)

   END FUNCTION skala_gpw_atom_subchunk_count

! **************************************************************************************************
!> \brief Build an atom-contiguous subchunk feature bundle from a rank-local atom chunk.
!> \param parent ...
!> \param features ...
!> \param subchunk_index ...
!> \param max_rows ...
!> \param requires_grad ...
! **************************************************************************************************
   SUBROUTINE skala_gpw_feature_build_atom_subchunk(parent, features, subchunk_index, &
                                                    max_rows, requires_grad)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: parent
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      INTEGER, INTENT(IN)                                :: subchunk_index, max_rows
      LOGICAL, INTENT(IN)                                :: requires_grad

      INTEGER                                            :: atom_begin, atom_count, atom_end, &
                                                            max_grid_size, row_begin, row_count, &
                                                            row_end

      CALL skala_gpw_feature_release(features)
      CPASSERT(parent%uses_atom_chunks)
      CALL atom_subchunk_bounds(subchunk_index, max_rows, atom_begin, atom_end, &
                                row_begin, row_end)
      atom_count = atom_end - atom_begin + 1
      row_count = row_end - row_begin + 1
      CPASSERT(atom_count > 0)
      CPASSERT(row_count > 0)
      MARK_USED(requires_grad)
      max_grid_size = MAXVAL(INT(cached_layout%chunk_atomic_grid_sizes(atom_begin:atom_end)))

      ALLOCATE (features%atomic_grid_size_bound_shape(0, max_grid_size))
      features%atomic_grid_size_bound_shape = 0_int_8

      features%chunk_feature_count = row_count
      features%nflat = parent%nflat
      features%nflat_local = parent%nflat_local
      features%grid_weight_sum = SUM(cached_layout%chunk_grid_weights(row_begin:row_end))
      features%uses_atom_chunks = .TRUE.
      features%uses_atom_chunk_routing = parent%uses_atom_chunk_routing
      CALL add_subchunk_feature_tensors(parent, features, atom_begin, atom_count, row_begin, &
                                        row_count)
      features%active = .TRUE.

   END SUBROUTINE skala_gpw_feature_build_atom_subchunk

! **************************************************************************************************
!> \brief Return atom and row bounds for an atom-contiguous rank-local subchunk.
!> \param subchunk_index ...
!> \param max_rows ...
!> \param atom_begin ...
!> \param atom_end ...
!> \param row_begin ...
!> \param row_end ...
! **************************************************************************************************
   SUBROUTINE atom_subchunk_bounds(subchunk_index, max_rows, atom_begin, atom_end, &
                                   row_begin, row_end)
      INTEGER, INTENT(IN)                                :: subchunk_index, max_rows
      INTEGER, INTENT(OUT)                               :: atom_begin, atom_end, row_begin, row_end

      INTEGER                                            :: atom_rows, current_subchunk, iatom, &
                                                            row_cursor, rows

      CPASSERT(subchunk_index > 0)
      CPASSERT(max_rows > 0)
      CPASSERT(cached_layout%chunk_natom > 0)

      atom_begin = 1
      atom_end = 0
      row_begin = 1
      row_end = 0
      current_subchunk = 1
      row_cursor = 1
      rows = 0
      DO iatom = 1, cached_layout%chunk_natom
         atom_rows = INT(cached_layout%chunk_atomic_grid_sizes(iatom))
         IF (rows > 0 .AND. rows + atom_rows > max_rows) THEN
            IF (current_subchunk == subchunk_index) THEN
               atom_end = iatom - 1
               row_end = row_cursor - 1
               RETURN
            END IF
            current_subchunk = current_subchunk + 1
            atom_begin = iatom
            row_begin = row_cursor
            rows = 0
         END IF
         rows = rows + atom_rows
         row_cursor = row_cursor + atom_rows
      END DO

      IF (current_subchunk == subchunk_index) THEN
         atom_end = cached_layout%chunk_natom
         row_end = row_cursor - 1
         RETURN
      END IF

      CPABORT("Requested native SKALA atom subchunk does not exist.")

   END SUBROUTINE atom_subchunk_bounds

! **************************************************************************************************
!> \brief Insert a subchunk into a Torch dictionary using static views of the cached chunk tensors.
!> \param parent ...
!> \param features ...
!> \param atom_begin ...
!> \param atom_count ...
!> \param row_begin ...
!> \param row_count ...
! **************************************************************************************************
   SUBROUTINE add_subchunk_feature_tensors(parent, features, atom_begin, atom_count, row_begin, &
                                           row_count)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: parent
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      INTEGER, INTENT(IN)                                :: atom_begin, atom_count, row_begin, &
                                                            row_count

      CPASSERT(cached_layout%chunk_static_tensors_active)
      CPASSERT(parent%active)
      CPASSERT(ALLOCATED(features%atomic_grid_size_bound_shape))

      features%owns_coordinate_tensor = .FALSE.
      features%owns_dynamic_tensors = .TRUE.
      features%owns_inputs = .TRUE.
      features%owns_static_tensors = .TRUE.
      features%uses_collapsed_rks_dynamic = parent%uses_collapsed_rks_dynamic

      CALL torch_tensor_narrow(cached_layout%chunk_grid_coords_t, 0, row_begin - 1, &
                               row_count, features%grid_coords_t)
      CALL torch_tensor_narrow(cached_layout%chunk_grid_weights_t, 0, row_begin - 1, &
                               row_count, features%grid_weights_t)
      CALL torch_tensor_narrow(cached_layout%chunk_atomic_grid_weights_t, 0, row_begin - 1, &
                               row_count, features%atomic_grid_weights_t)
      CALL torch_tensor_narrow(cached_layout%chunk_atomic_grid_sizes_t, 0, atom_begin - 1, &
                               atom_count, features%atomic_grid_sizes_t)
      CALL torch_tensor_narrow(cached_layout%chunk_coarse_0_atomic_coords_t, 0, &
                               atom_begin - 1, atom_count, features%coarse_0_atomic_coords_t)
      CALL torch_tensor_from_array(features%atomic_grid_size_bound_shape_t, &
                                   features%atomic_grid_size_bound_shape)
      CALL torch_tensor_to_device_leaf(features%atomic_grid_size_bound_shape_t, .FALSE.)
      CALL torch_tensor_narrow(parent%density_t, 1, row_begin - 1, row_count, &
                               features%density_t)
      CALL torch_tensor_narrow(parent%grad_t, 2, row_begin - 1, row_count, features%grad_t)
      CALL torch_tensor_narrow(parent%kin_t, 1, row_begin - 1, row_count, features%kin_t)
      IF (features%uses_collapsed_rks_dynamic) THEN
         CALL torch_tensor_expand_dim(features%density_t, 0, 2, features%density_input_t)
         CALL torch_tensor_expand_dim(features%grad_t, 0, 2, features%grad_input_t)
         CALL torch_tensor_expand_dim(features%kin_t, 0, 2, features%kin_input_t)
      END IF

      CALL torch_dict_create(features%inputs)
      CALL torch_dict_insert(features%inputs, "grid_coords", features%grid_coords_t)
      CALL torch_dict_insert(features%inputs, "grid_weights", features%grid_weights_t)
      CALL torch_dict_insert(features%inputs, "atomic_grid_weights", &
                             features%atomic_grid_weights_t)
      CALL torch_dict_insert(features%inputs, "atomic_grid_sizes", &
                             features%atomic_grid_sizes_t)
      CALL torch_dict_insert(features%inputs, "atomic_grid_size_bound_shape", &
                             features%atomic_grid_size_bound_shape_t)
      IF (features%uses_collapsed_rks_dynamic) THEN
         CALL torch_dict_insert(features%inputs, "density", features%density_input_t)
         CALL torch_dict_insert(features%inputs, "grad", features%grad_input_t)
         CALL torch_dict_insert(features%inputs, "kin", features%kin_input_t)
      ELSE
         CALL torch_dict_insert(features%inputs, "density", features%density_t)
         CALL torch_dict_insert(features%inputs, "grad", features%grad_t)
         CALL torch_dict_insert(features%inputs, "kin", features%kin_t)
      END IF
      CALL torch_dict_insert(features%inputs, "coarse_0_atomic_coords", &
                             features%coarse_0_atomic_coords_t)

   END SUBROUTINE add_subchunk_feature_tensors

! **************************************************************************************************
!> \brief Insert owned subchunk arrays into a Torch dictionary.
!> \param features ...
!> \param requires_grad ...
! **************************************************************************************************
   SUBROUTINE add_owned_feature_tensors(features, requires_grad)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      LOGICAL, INTENT(IN)                                :: requires_grad

      CPASSERT(ALLOCATED(features%chunk_density))
      CPASSERT(ALLOCATED(features%chunk_grad))
      CPASSERT(ALLOCATED(features%chunk_kin))
      CPASSERT(ALLOCATED(features%grid_coords))
      CPASSERT(ALLOCATED(features%grid_weights))
      CPASSERT(ALLOCATED(features%atomic_grid_weights))
      CPASSERT(ALLOCATED(features%atomic_grid_sizes))
      CPASSERT(ALLOCATED(features%atomic_grid_size_bound_shape))
      CPASSERT(ALLOCATED(features%coarse_0_atomic_coords))

      features%owns_coordinate_tensor = .FALSE.
      features%owns_dynamic_tensors = .TRUE.
      features%owns_inputs = .TRUE.
      features%owns_static_tensors = .TRUE.

      CALL torch_tensor_from_array(features%grid_coords_t, features%grid_coords)
      CALL torch_tensor_to_device_leaf(features%grid_coords_t, .FALSE.)
      CALL torch_tensor_from_array(features%grid_weights_t, features%grid_weights)
      CALL torch_tensor_to_device_leaf(features%grid_weights_t, .FALSE.)
      CALL torch_tensor_from_array(features%atomic_grid_weights_t, features%atomic_grid_weights)
      CALL torch_tensor_to_device_leaf(features%atomic_grid_weights_t, .FALSE.)
      CALL torch_tensor_from_array(features%atomic_grid_sizes_t, features%atomic_grid_sizes)
      CALL torch_tensor_to_device_leaf(features%atomic_grid_sizes_t, .FALSE.)
      CALL torch_tensor_from_array(features%coarse_0_atomic_coords_t, &
                                   features%coarse_0_atomic_coords)
      CALL torch_tensor_to_device_leaf(features%coarse_0_atomic_coords_t, .FALSE.)
      CALL torch_tensor_from_array(features%atomic_grid_size_bound_shape_t, &
                                   features%atomic_grid_size_bound_shape)
      CALL torch_tensor_to_device_leaf(features%atomic_grid_size_bound_shape_t, .FALSE.)
      CALL torch_tensor_from_array(features%density_t, features%chunk_density)
      CALL torch_tensor_to_device_leaf(features%density_t, requires_grad)
      CALL torch_tensor_from_array(features%grad_t, features%chunk_grad)
      CALL torch_tensor_to_device_leaf(features%grad_t, requires_grad)
      CALL torch_tensor_from_array(features%kin_t, features%chunk_kin)
      CALL torch_tensor_to_device_leaf(features%kin_t, requires_grad)

      CALL torch_dict_create(features%inputs)
      CALL torch_dict_insert(features%inputs, "grid_coords", features%grid_coords_t)
      CALL torch_dict_insert(features%inputs, "grid_weights", features%grid_weights_t)
      CALL torch_dict_insert(features%inputs, "atomic_grid_weights", &
                             features%atomic_grid_weights_t)
      CALL torch_dict_insert(features%inputs, "atomic_grid_sizes", &
                             features%atomic_grid_sizes_t)
      CALL torch_dict_insert(features%inputs, "atomic_grid_size_bound_shape", &
                             features%atomic_grid_size_bound_shape_t)
      CALL torch_dict_insert(features%inputs, "density", features%density_t)
      CALL torch_dict_insert(features%inputs, "grad", features%grad_t)
      CALL torch_dict_insert(features%inputs, "kin", features%kin_t)
      CALL torch_dict_insert(features%inputs, "coarse_0_atomic_coords", &
                             features%coarse_0_atomic_coords_t)

   END SUBROUTINE add_owned_feature_tensors

! **************************************************************************************************
!> \brief Insert all SKALA feature tensors into the Torch dictionary.
!> \param features ...
!> \param requires_grad ...
!> \param requires_coordinate_grad ...
!> \param requires_stress_grad ...
!> \param use_atom_chunks ...
!> \param requires_weight_grad ...
! **************************************************************************************************
   SUBROUTINE add_feature_tensors(features, requires_grad, requires_coordinate_grad, &
                                  requires_stress_grad, use_atom_chunks, requires_weight_grad)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      LOGICAL, INTENT(IN)                                :: requires_grad, requires_coordinate_grad, &
                                                            requires_stress_grad, use_atom_chunks
      LOGICAL, INTENT(IN), OPTIONAL                      :: requires_weight_grad

      LOGICAL                                            :: my_requires_weight_grad

      my_requires_weight_grad = .FALSE.
      IF (PRESENT(requires_weight_grad)) my_requires_weight_grad = requires_weight_grad

      CPASSERT(cached_layout%static_tensors_active)
      features%owns_static_tensors = .FALSE.
      features%owns_coordinate_tensor = .FALSE.
      features%owns_grid_coordinate_tensor = .FALSE.
      features%owns_weight_tensors = .FALSE.
      features%owns_dynamic_tensors = .FALSE.
      features%owns_inputs = .TRUE.
      IF (use_atom_chunks) THEN
         CPASSERT(.NOT. requires_coordinate_grad)
         CPASSERT(.NOT. requires_stress_grad)
         CPASSERT(.NOT. my_requires_weight_grad)
         CPASSERT(cached_layout%chunk_static_tensors_active)
         features%grid_coords_t = cached_layout%chunk_grid_coords_t
         features%grid_weights_t = cached_layout%chunk_grid_weights_t
         features%atomic_grid_weights_t = cached_layout%chunk_atomic_grid_weights_t
         features%atomic_grid_sizes_t = cached_layout%chunk_atomic_grid_sizes_t
         features%atomic_grid_size_bound_shape_t = &
            cached_layout%chunk_atomic_grid_size_bound_shape_t
         features%local_feature_indices_t = cached_layout%chunk_feature_indices_t

         IF (cached_layout%chunk_inputs_active .AND. &
             (cached_layout%chunk_inputs_use_collapsed_rks .NEQV. &
              features%uses_collapsed_rks_dynamic)) THEN
            CALL torch_dict_release(cached_layout%chunk_inputs)
            cached_layout%chunk_inputs_active = .FALSE.
         END IF
         IF (.NOT. features%uses_collapsed_rks_dynamic .AND. &
             cached_layout%chunk_dynamic_input_views_active) THEN
            CALL torch_tensor_release(cached_layout%chunk_density_input_t)
            CALL torch_tensor_release(cached_layout%chunk_grad_input_t)
            CALL torch_tensor_release(cached_layout%chunk_kin_input_t)
            cached_layout%chunk_dynamic_input_views_active = .FALSE.
         END IF

         CALL torch_tensor_reset_from_array(cached_layout%chunk_density_t, &
                                            features%chunk_density, requires_grad=requires_grad)
         features%density_t = cached_layout%chunk_density_t
         CALL torch_tensor_reset_from_array(cached_layout%chunk_grad_t, features%chunk_grad, &
                                            requires_grad=requires_grad)
         features%grad_t = cached_layout%chunk_grad_t
         CALL torch_tensor_reset_from_array(cached_layout%chunk_kin_t, features%chunk_kin, &
                                            requires_grad=requires_grad)
         features%kin_t = cached_layout%chunk_kin_t
         cached_layout%chunk_dynamic_tensors_active = .TRUE.

         IF (features%uses_collapsed_rks_dynamic .AND. &
             .NOT. cached_layout%chunk_dynamic_input_views_active) THEN
            CALL torch_tensor_expand_dim(cached_layout%chunk_density_t, 0, 2, &
                                         cached_layout%chunk_density_input_t)
            CALL torch_tensor_expand_dim(cached_layout%chunk_grad_t, 0, 2, &
                                         cached_layout%chunk_grad_input_t)
            CALL torch_tensor_expand_dim(cached_layout%chunk_kin_t, 0, 2, &
                                         cached_layout%chunk_kin_input_t)
            cached_layout%chunk_dynamic_input_views_active = .TRUE.
         END IF
         IF (features%uses_collapsed_rks_dynamic) THEN
            features%density_input_t = cached_layout%chunk_density_input_t
            features%grad_input_t = cached_layout%chunk_grad_input_t
            features%kin_input_t = cached_layout%chunk_kin_input_t
         END IF

         IF (.NOT. cached_layout%chunk_inputs_active) THEN
            CALL torch_dict_clone(cached_layout%chunk_static_inputs, cached_layout%chunk_inputs)
            IF (features%uses_collapsed_rks_dynamic) THEN
               CALL torch_dict_insert(cached_layout%chunk_inputs, "density", &
                                      features%density_input_t)
               CALL torch_dict_insert(cached_layout%chunk_inputs, "grad", &
                                      features%grad_input_t)
               CALL torch_dict_insert(cached_layout%chunk_inputs, "kin", &
                                      features%kin_input_t)
            ELSE
               CALL torch_dict_insert(cached_layout%chunk_inputs, "density", &
                                      cached_layout%chunk_density_t)
               CALL torch_dict_insert(cached_layout%chunk_inputs, "grad", &
                                      cached_layout%chunk_grad_t)
               CALL torch_dict_insert(cached_layout%chunk_inputs, "kin", &
                                      cached_layout%chunk_kin_t)
            END IF
            CALL torch_dict_insert(cached_layout%chunk_inputs, "coarse_0_atomic_coords", &
                                   cached_layout%chunk_coarse_0_atomic_coords_t)
            cached_layout%chunk_inputs_use_collapsed_rks = features%uses_collapsed_rks_dynamic
            cached_layout%chunk_inputs_active = .TRUE.
         END IF
         features%inputs = cached_layout%chunk_inputs
         features%owns_inputs = .FALSE.
         features%coarse_0_atomic_coords_t = cached_layout%chunk_coarse_0_atomic_coords_t
      ELSE
         IF (.NOT. requires_stress_grad .AND. .NOT. my_requires_weight_grad) THEN
            features%grid_coords_t = cached_layout%grid_coords_t
            features%grid_weights_t = cached_layout%grid_weights_t
            features%atomic_grid_weights_t = cached_layout%atomic_grid_weights_t
         END IF
         features%atomic_grid_sizes_t = cached_layout%atomic_grid_sizes_t
         features%atomic_grid_size_bound_shape_t = cached_layout%atomic_grid_size_bound_shape_t
         features%local_feature_indices_t = cached_layout%local_feature_indices_t

         CALL torch_tensor_reset_from_array(cached_layout%density_t, features%density, &
                                            requires_grad=requires_grad)
         features%density_t = cached_layout%density_t
         CALL torch_tensor_reset_from_array(cached_layout%grad_t, features%grad, &
                                            requires_grad=requires_grad)
         features%grad_t = cached_layout%grad_t
         CALL torch_tensor_reset_from_array(cached_layout%kin_t, features%kin, &
                                            requires_grad=requires_grad)
         features%kin_t = cached_layout%kin_t
         cached_layout%dynamic_tensors_active = .TRUE.

         IF (requires_coordinate_grad .OR. requires_stress_grad .OR. my_requires_weight_grad) THEN
            IF (requires_stress_grad .OR. my_requires_weight_grad) THEN
               CALL torch_dict_create(features%inputs)
               IF (requires_stress_grad) THEN
                  CALL torch_tensor_from_array(features%grid_coords_t, features%grid_coords)
                  CALL torch_tensor_to_device_leaf(features%grid_coords_t, .TRUE.)
                  CALL torch_dict_insert(features%inputs, "grid_coords", features%grid_coords_t)
                  features%owns_grid_coordinate_tensor = .TRUE.
               ELSE
                  features%grid_coords_t = cached_layout%grid_coords_t
                  CALL torch_dict_insert(features%inputs, "grid_coords", features%grid_coords_t)
               END IF
               CALL torch_tensor_from_array(features%grid_weights_t, features%grid_weights)
               CALL torch_tensor_to_device_leaf(features%grid_weights_t, .TRUE.)
               CALL torch_tensor_from_array(features%atomic_grid_weights_t, &
                                            features%atomic_grid_weights)
               CALL torch_tensor_to_device_leaf(features%atomic_grid_weights_t, .TRUE.)
               CALL torch_dict_insert(features%inputs, "grid_weights", features%grid_weights_t)
               CALL torch_dict_insert(features%inputs, "atomic_grid_weights", &
                                      features%atomic_grid_weights_t)
               CALL torch_dict_insert(features%inputs, "atomic_grid_sizes", &
                                      features%atomic_grid_sizes_t)
               CALL torch_dict_insert(features%inputs, "atomic_grid_size_bound_shape", &
                                      features%atomic_grid_size_bound_shape_t)
               features%owns_weight_tensors = .TRUE.
            ELSE
               CALL torch_dict_clone(cached_layout%static_inputs, features%inputs)
            END IF
            CALL torch_dict_insert(features%inputs, "density", features%density_t)
            CALL torch_dict_insert(features%inputs, "grad", features%grad_t)
            CALL torch_dict_insert(features%inputs, "kin", features%kin_t)
         ELSE
            IF (.NOT. cached_layout%inputs_active) THEN
               CALL torch_dict_clone(cached_layout%static_inputs, cached_layout%inputs)
               CALL torch_dict_insert(cached_layout%inputs, "density", cached_layout%density_t)
               CALL torch_dict_insert(cached_layout%inputs, "grad", cached_layout%grad_t)
               CALL torch_dict_insert(cached_layout%inputs, "kin", cached_layout%kin_t)
               CALL torch_dict_insert(cached_layout%inputs, "coarse_0_atomic_coords", &
                                      cached_layout%coarse_0_atomic_coords_t)
               cached_layout%inputs_active = .TRUE.
            END IF
            features%inputs = cached_layout%inputs
            features%owns_inputs = .FALSE.
            features%coarse_0_atomic_coords_t = cached_layout%coarse_0_atomic_coords_t
         END IF
      END IF

      IF (requires_coordinate_grad .OR. requires_stress_grad) THEN
         CPASSERT(.NOT. use_atom_chunks)
         CALL torch_tensor_from_array(features%coarse_0_atomic_coords_t, &
                                      features%coarse_0_atomic_coords)
         CALL torch_tensor_to_device_leaf(features%coarse_0_atomic_coords_t, .TRUE.)
         CALL torch_dict_insert(features%inputs, "coarse_0_atomic_coords", &
                                features%coarse_0_atomic_coords_t)
         features%owns_coordinate_tensor = .TRUE.
      END IF

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
!> \brief Build Becke-like smooth atom weights for one native-grid point.
!> \param grid_point ...
!> \param atom_coords ...
!> \param cell ...
!> \param weights ...
!> \param atom_image_coords ...
!> \param distances ...
! **************************************************************************************************
   SUBROUTINE smooth_atom_partition(grid_point, atom_coords, cell, weights, atom_image_coords, &
                                    distances)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: grid_point
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: atom_coords
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: weights
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: atom_image_coords
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: distances

      INTEGER                                            :: iatom, jatom, natom
      REAL(KIND=dp)                                      :: mu, rab, rsum, switch, total
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      REAL(KIND=dp), DIMENSION(3, SIZE(atom_coords, 2))  :: partition_atom_coords

      natom = SIZE(atom_coords, 2)
      CPASSERT(SIZE(weights) == natom)
      CPASSERT(SIZE(atom_image_coords, 1) == 3)
      CPASSERT(SIZE(atom_image_coords, 2) == natom)
      CPASSERT(SIZE(distances) == natom)

      DO iatom = 1, natom
         atom_image_coords(:, iatom) = &
            nearest_image_coordinate(atom_coords(:, iatom), grid_point, cell)
         partition_atom_coords(:, iatom) = &
            nearest_atom_image_coordinate(atom_coords(:, iatom), grid_point, cell)
         rij = grid_point - partition_atom_coords(:, iatom)
         distances(iatom) = SQRT(SUM(rij**2))
      END DO

      weights = 1.0_dp
      DO iatom = 1, natom - 1
         DO jatom = iatom + 1, natom
            rij = partition_atom_coords(:, iatom) - partition_atom_coords(:, jatom)
            rab = SQRT(SUM(rij**2))
            IF (rab <= layout_tol) CYCLE
            mu = (distances(iatom) - distances(jatom))/rab
            mu = MAX(-1.0_dp, MIN(1.0_dp, mu))
            switch = 0.5_dp*(1.0_dp - becke_shape(mu))
            weights(iatom) = weights(iatom)*switch
            weights(jatom) = weights(jatom)*(1.0_dp - switch)
         END DO
      END DO

      total = SUM(weights)
      IF (total > 0.0_dp) THEN
         weights = weights/total
      ELSE
         rsum = HUGE(1.0_dp)
         jatom = 1
         DO iatom = 1, natom
            IF (distances(iatom) < rsum) THEN
               rsum = distances(iatom)
               jatom = iatom
            END IF
         END DO
         weights = 0.0_dp
         weights(jatom) = 1.0_dp
      END IF

   END SUBROUTINE smooth_atom_partition

! **************************************************************************************************
!> \brief Build smooth atom weights and their atom/cell deformation derivatives.
!> \param grid_point ...
!> \param atom_coords ...
!> \param cell ...
!> \param weights ...
!> \param included ...
!> \param dweights_datom ...
!> \param dweights_dstrain ...
! **************************************************************************************************
   SUBROUTINE skala_gpw_smooth_partition_derivatives(grid_point, atom_coords, cell, &
                                                     weights, included, dweights_datom, &
                                                     dweights_dstrain)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: grid_point
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: atom_coords
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: weights
      LOGICAL, DIMENSION(:), INTENT(OUT)                 :: included
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: dweights_datom, dweights_dstrain

      INTEGER                                            :: iatom, idir, jatom, jdir, natom
      REAL(KIND=dp)                                      :: dist_diff, ds_dmu, included_sum, mu, &
                                                            mu_raw, one_minus_switch, rab, rsum, &
                                                            switch, total
      REAL(KIND=dp), DIMENSION(3)                        :: dmu_atom_i, dmu_atom_j, ds_atom_i, &
                                                            ds_atom_j, pair, unit_pair
      REAL(KIND=dp), DIMENSION(3, 3)                     :: dmu_strain, ds_strain, mean_strain
      REAL(KIND=dp), DIMENSION(3, SIZE(atom_coords, 2), &
         SIZE(atom_coords, 2))                           :: log_weight_atom
      REAL(KIND=dp), DIMENSION(3, SIZE(atom_coords, 2))  :: mean_atom, partition_atom_coords, rvecs, &
                                                            unit_rvecs
      REAL(KIND=dp), &
         DIMENSION(3, 3, SIZE(atom_coords, 2))           :: log_weight_strain
      REAL(KIND=dp), DIMENSION(SIZE(atom_coords, 2))     :: distances, normalized_weights, &
                                                            raw_weights

      natom = SIZE(atom_coords, 2)
      CPASSERT(SIZE(weights) == natom)
      CPASSERT(SIZE(included) == natom)
      CPASSERT(SIZE(dweights_datom, 1) == 3)
      CPASSERT(SIZE(dweights_datom, 2) == natom)
      CPASSERT(SIZE(dweights_datom, 3) == natom)
      CPASSERT(SIZE(dweights_dstrain, 1) == 3)
      CPASSERT(SIZE(dweights_dstrain, 2) == 3)
      CPASSERT(SIZE(dweights_dstrain, 3) == natom)

      weights = 0.0_dp
      included = .FALSE.
      dweights_datom = 0.0_dp
      dweights_dstrain = 0.0_dp
      raw_weights = 1.0_dp
      log_weight_atom = 0.0_dp
      log_weight_strain = 0.0_dp

      DO iatom = 1, natom
         partition_atom_coords(:, iatom) = &
            nearest_atom_image_coordinate(atom_coords(:, iatom), grid_point, cell)
         rvecs(:, iatom) = grid_point - partition_atom_coords(:, iatom)
         distances(iatom) = SQRT(SUM(rvecs(:, iatom)**2))
         IF (distances(iatom) > layout_tol) THEN
            unit_rvecs(:, iatom) = rvecs(:, iatom)/distances(iatom)
         ELSE
            unit_rvecs(:, iatom) = 0.0_dp
         END IF
      END DO

      DO iatom = 1, natom - 1
         DO jatom = iatom + 1, natom
            pair = partition_atom_coords(:, iatom) - partition_atom_coords(:, jatom)
            rab = SQRT(SUM(pair**2))
            IF (rab <= layout_tol) CYCLE
            unit_pair = pair/rab
            dist_diff = distances(iatom) - distances(jatom)
            mu_raw = dist_diff/rab
            mu = MAX(-1.0_dp, MIN(1.0_dp, mu_raw))
            switch = 0.5_dp*(1.0_dp - becke_shape(mu))
            one_minus_switch = 1.0_dp - switch

            IF (ABS(mu_raw) < 1.0_dp) THEN
               ds_dmu = -0.5_dp*becke_shape_derivative(mu)
            ELSE
               ds_dmu = 0.0_dp
            END IF
            IF (ABS(ds_dmu) > 0.0_dp .AND. switch > TINY(1.0_dp) .AND. &
                one_minus_switch > TINY(1.0_dp)) THEN
               dmu_atom_i = (-unit_rvecs(:, iatom)*rab - dist_diff*unit_pair)/rab**2
               dmu_atom_j = (unit_rvecs(:, jatom)*rab + dist_diff*unit_pair)/rab**2
               ds_atom_i = ds_dmu*dmu_atom_i
               ds_atom_j = ds_dmu*dmu_atom_j
               log_weight_atom(:, iatom, iatom) = &
                  log_weight_atom(:, iatom, iatom) + ds_atom_i/switch
               log_weight_atom(:, iatom, jatom) = &
                  log_weight_atom(:, iatom, jatom) - ds_atom_i/one_minus_switch
               log_weight_atom(:, jatom, iatom) = &
                  log_weight_atom(:, jatom, iatom) + ds_atom_j/switch
               log_weight_atom(:, jatom, jatom) = &
                  log_weight_atom(:, jatom, jatom) - ds_atom_j/one_minus_switch

               DO idir = 1, 3
                  DO jdir = 1, 3
                     dmu_strain(idir, jdir) = &
                        ((unit_rvecs(idir, iatom)*rvecs(jdir, iatom) - &
                          unit_rvecs(idir, jatom)*rvecs(jdir, jatom))*rab - &
                         dist_diff*unit_pair(idir)*pair(jdir))/rab**2
                  END DO
               END DO
               ds_strain = ds_dmu*dmu_strain
               log_weight_strain(:, :, iatom) = &
                  log_weight_strain(:, :, iatom) + ds_strain/switch
               log_weight_strain(:, :, jatom) = &
                  log_weight_strain(:, :, jatom) - ds_strain/one_minus_switch
            END IF

            raw_weights(iatom) = raw_weights(iatom)*switch
            raw_weights(jatom) = raw_weights(jatom)*one_minus_switch
         END DO
      END DO

      total = SUM(raw_weights)
      IF (total > 0.0_dp) THEN
         normalized_weights = raw_weights/total
         included = normalized_weights > smooth_partition_eps
      ELSE
         rsum = HUGE(1.0_dp)
         jatom = 1
         DO iatom = 1, natom
            IF (distances(iatom) < rsum) THEN
               rsum = distances(iatom)
               jatom = iatom
            END IF
         END DO
         included(jatom) = .TRUE.
         weights(jatom) = 1.0_dp
         RETURN
      END IF

      included_sum = SUM(raw_weights, MASK=included)
      IF (included_sum <= 0.0_dp) THEN
         rsum = HUGE(1.0_dp)
         jatom = 1
         DO iatom = 1, natom
            IF (distances(iatom) < rsum) THEN
               rsum = distances(iatom)
               jatom = iatom
            END IF
         END DO
         included = .FALSE.
         included(jatom) = .TRUE.
         weights = 0.0_dp
         weights(jatom) = 1.0_dp
         RETURN
      END IF

      DO iatom = 1, natom
         IF (included(iatom)) weights(iatom) = raw_weights(iatom)/included_sum
      END DO

      mean_atom = 0.0_dp
      mean_strain = 0.0_dp
      DO iatom = 1, natom
         IF (.NOT. included(iatom)) CYCLE
         mean_strain = mean_strain + weights(iatom)*log_weight_strain(:, :, iatom)
         DO jatom = 1, natom
            mean_atom(:, jatom) = mean_atom(:, jatom) + &
                                  weights(iatom)*log_weight_atom(:, jatom, iatom)
         END DO
      END DO

      DO iatom = 1, natom
         IF (.NOT. included(iatom)) CYCLE
         dweights_dstrain(:, :, iatom) = &
            weights(iatom)*(log_weight_strain(:, :, iatom) - mean_strain)
         DO jatom = 1, natom
            dweights_datom(:, jatom, iatom) = &
               weights(iatom)*(log_weight_atom(:, jatom, iatom) - mean_atom(:, jatom))
         END DO
      END DO

   END SUBROUTINE skala_gpw_smooth_partition_derivatives

! **************************************************************************************************
!> \brief Becke fuzzy-cell shape function.
!> \param mu ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION becke_shape(mu) RESULT(val)
      REAL(KIND=dp), INTENT(IN)                          :: mu
      REAL(KIND=dp)                                      :: val

      INTEGER                                            :: iter

      val = mu
      DO iter = 1, 3
         val = 0.5_dp*val*(3.0_dp - val*val)
      END DO

   END FUNCTION becke_shape

! **************************************************************************************************
!> \brief Derivative of the Becke fuzzy-cell shape function.
!> \param mu ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION becke_shape_derivative(mu) RESULT(val)
      REAL(KIND=dp), INTENT(IN)                          :: mu
      REAL(KIND=dp)                                      :: val

      INTEGER                                            :: iter
      REAL(KIND=dp)                                      :: x

      x = mu
      val = 1.0_dp
      DO iter = 1, 3
         val = val*1.5_dp*(1.0_dp - x*x)
         x = 0.5_dp*x*(3.0_dp - x*x)
      END DO

   END FUNCTION becke_shape_derivative

! **************************************************************************************************
!> \brief Return the atom image nearest to a regular-grid point.
!> \param atom_coord ...
!> \param grid_point ...
!> \param cell ...
!> \return ...
! **************************************************************************************************
   FUNCTION nearest_atom_image_coordinate(atom_coord, grid_point, cell) RESULT(coord)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: atom_coord, grid_point
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(3)                        :: coord

      REAL(KIND=dp)                                      :: dx, dy, dz

      IF (cell%orthorhombic) THEN
         dx = atom_coord(1) - grid_point(1)
         dy = atom_coord(2) - grid_point(2)
         dz = atom_coord(3) - grid_point(3)
         dx = dx - cell%hmat(1, 1)*cell%perd(1)*ANINT(cell%h_inv(1, 1)*dx)
         dy = dy - cell%hmat(2, 2)*cell%perd(2)*ANINT(cell%h_inv(2, 2)*dy)
         dz = dz - cell%hmat(3, 3)*cell%perd(3)*ANINT(cell%h_inv(3, 3)*dz)
         coord = grid_point + [dx, dy, dz]
      ELSE
         coord = grid_point + pbc(grid_point, atom_coord, cell)
      END IF

   END FUNCTION nearest_atom_image_coordinate

! **************************************************************************************************
!> \brief Return the grid-point image nearest to the owning atom coordinate.
!> \param owner_coord ...
!> \param grid_point ...
!> \param cell ...
!> \return ...
! **************************************************************************************************
   FUNCTION nearest_image_coordinate(owner_coord, grid_point, cell) RESULT(coord)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: owner_coord, grid_point
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(3)                        :: coord

      REAL(KIND=dp)                                      :: dx, dy, dz

      IF (cell%orthorhombic) THEN
         dx = grid_point(1) - owner_coord(1)
         dy = grid_point(2) - owner_coord(2)
         dz = grid_point(3) - owner_coord(3)
         dx = dx - cell%hmat(1, 1)*cell%perd(1)*ANINT(cell%h_inv(1, 1)*dx)
         dy = dy - cell%hmat(2, 2)*cell%perd(2)*ANINT(cell%h_inv(2, 2)*dy)
         dz = dz - cell%hmat(3, 3)*cell%perd(3)*ANINT(cell%h_inv(3, 3)*dz)
         coord = owner_coord + [dx, dy, dz]
      ELSE
         coord = owner_coord + pbc(owner_coord, grid_point, cell)
      END IF

   END FUNCTION nearest_image_coordinate

! **************************************************************************************************
!> \brief Assign a grid point to the nearest periodic atom.
!> \param grid_point ...
!> \param atom_coords ...
!> \param cell ...
!> \return ...
! **************************************************************************************************
   FUNCTION nearest_atom(grid_point, atom_coords, cell) RESULT(owner)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: grid_point
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: atom_coords
      TYPE(cell_type), POINTER                           :: cell
      INTEGER                                            :: owner

      INTEGER                                            :: iatom
      REAL(KIND=dp)                                      :: best_r2, dx, dy, dz, r2
      REAL(KIND=dp), DIMENSION(3)                        :: rij

      owner = 1
      best_r2 = HUGE(1.0_dp)
      IF (cell%orthorhombic) THEN
         DO iatom = 1, SIZE(atom_coords, 2)
            dx = grid_point(1) - atom_coords(1, iatom)
            dy = grid_point(2) - atom_coords(2, iatom)
            dz = grid_point(3) - atom_coords(3, iatom)
            dx = dx - cell%hmat(1, 1)*cell%perd(1)*ANINT(cell%h_inv(1, 1)*dx)
            dy = dy - cell%hmat(2, 2)*cell%perd(2)*ANINT(cell%h_inv(2, 2)*dy)
            dz = dz - cell%hmat(3, 3)*cell%perd(3)*ANINT(cell%h_inv(3, 3)*dz)
            r2 = dx*dx + dy*dy + dz*dz
            IF (r2 < best_r2) THEN
               best_r2 = r2
               owner = iatom
            END IF
         END DO
      ELSE
         DO iatom = 1, SIZE(atom_coords, 2)
            rij = pbc(grid_point, atom_coords(:, iatom), cell)
            r2 = SUM(rij**2)
            IF (r2 < best_r2) THEN
               best_r2 = r2
               owner = iatom
            END IF
         END DO
      END IF

   END FUNCTION nearest_atom

END MODULE skala_gpw_features

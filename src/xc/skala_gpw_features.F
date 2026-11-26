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
   USE skala_env_types,                 ONLY: &
        atom_chunk_layout_cache_type, skala_atom_ws_type, skala_execution_target_type, &
        skala_gapw_atom_cuda_tensor_cache_type, skala_gpw_feature_type, &
        skala_gpw_layout_cache_type, skala_gpw_scratch_type, skala_gpw_state_type, &
        skala_grad_scratch_type, skala_model_context_type, skala_state, &
        skala_weight_grad_cache_type
   USE skala_feature_keys,              ONLY: &
        skala_key_atomic_grid_size_bound_shape, skala_key_atomic_grid_sizes, &
        skala_key_atomic_grid_weights, skala_key_coarse_0_atomic_coords, skala_key_density, &
        skala_key_grad, skala_key_grid_coords, skala_key_grid_weights, skala_key_kin
   USE skala_periodic_partition,        ONLY: build_periodic_atom_image_layout,&
                                              periodic_atom_image_partition_from_layout,&
                                              skala_periodic_partition_type
   USE torch_api,                       ONLY: &
        staged_upload, torch_dict_clone, torch_dict_create, torch_dict_insert, torch_dict_release, &
        torch_dict_set, torch_pinned_finish, torch_pinned_reserve, torch_stream_synchronize, &
        torch_tensor_expand_dim, torch_tensor_from_array, torch_tensor_narrow, &
        torch_tensor_release, torch_tensor_set_requires_grad, torch_tensor_to_device_leaf, &
        torch_tensor_type
   USE xc_input_constants,              ONLY: skala_gpw_atom_partition_hard,&
                                              skala_gpw_atom_partition_smooth
   USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                              xc_rho_set_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   REAL(KIND=dp), PARAMETER, PRIVATE    :: layout_tol = 1.0E-12_dp
   INTEGER, PARAMETER, PRIVATE          :: ndynamic_per_point = 10, nrks_dynamic_per_point = 5, &
                                           nstatic_per_point = 5
   REAL(KIND=dp), PARAMETER, PRIVATE    :: smooth_partition_eps = 1.0E-12_dp
   REAL(KIND=dp), PARAMETER, PRIVATE    :: smooth_partition_taper_end = &
                                           10.0_dp*smooth_partition_eps

   PUBLIC :: skala_gpw_atom_subchunk_count, skala_gpw_atom_subchunk_layout, &
             skala_state, skala_gpw_state_type, skala_execution_target_type, &
             skala_gpw_layout_cache_type, skala_gpw_feature_type, skala_gpw_scratch_type, &
             skala_gapw_atom_cuda_tensor_cache_type, skala_atom_ws_type, &
             skala_model_context_type, atom_chunk_layout_cache_type, &
             skala_grad_scratch_type, skala_weight_grad_cache_type, &
             sort_chunks_by_padded_rows, grid_coordinate, allocate_empty_bound_shape, &
             skala_gpw_features_release_cache, &
             skala_gpw_feature_build, skala_gpw_feature_build_atom_subchunk_bounds, &
             skala_gpw_feature_release, &
             smooth_partition_derivative_action, &
             smooth_partition_weights, &
             build_periodic_atom_image_layout, &
             periodic_atom_image_partition_from_layout, skala_periodic_partition_type, &
             smooth_partition_atomic_weight_scale, smooth_partition_atomic_weight_scale_derivative

   INTERFACE ensure_real_1d
      MODULE PROCEDURE ensure_real_1d_alloc, ensure_real_1d_ptr
   END INTERFACE ensure_real_1d
   INTERFACE ensure_real_2d
      MODULE PROCEDURE ensure_real_2d_alloc, ensure_real_2d_ptr
   END INTERFACE ensure_real_2d
   INTERFACE ensure_int_1d
      MODULE PROCEDURE ensure_int_1d_alloc, ensure_int_1d_ptr
   END INTERFACE ensure_int_1d

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
!> \param requires_weight_grad ...
!> \param compute_diagnostics ...
!> \param atom_partition ...
!> \param cuda_device selected Torch device (-1 for CPU, >=0 for CUDA)
! **************************************************************************************************
   SUBROUTINE skala_gpw_feature_build(features, rho_set, rho_r, particle_set, cell, &
                                      requires_grad, weights, requires_coordinate_grad, &
                                      requires_stress_grad, use_atom_chunks, requires_weight_grad, &
                                      compute_diagnostics, atom_partition, cuda_device)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(pw_r3d_rs_type), DIMENSION(:), INTENT(IN)     :: rho_r
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN), OPTIONAL                      :: requires_grad
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: weights
      LOGICAL, INTENT(IN), OPTIONAL :: requires_coordinate_grad, requires_stress_grad, &
         use_atom_chunks, requires_weight_grad, compute_diagnostics
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_partition
      INTEGER, INTENT(IN)                                :: cuda_device

      INTEGER :: handle, i, ipt, ispin, j, k, local_row, my_atom_partition, &
         ndynamic_local_per_point, nflat, nflat_local, nspins, phase_handle, real_base, row
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL :: collapse_spin_dynamics, my_compute_diagnostics, my_requires_coordinate_grad, &
         my_requires_grad, my_requires_stress_grad, my_requires_weight_grad, my_use_atom_chunks
      REAL(KIND=dp)                                      :: grid_weight_local
      REAL(KIND=dp), DIMENSION(3)                        :: integrals
      REAL(KIND=dp), DIMENSION(:), POINTER               :: global_dynamic, local_dynamic
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
      my_requires_weight_grad = .FALSE.
      IF (PRESENT(requires_weight_grad)) my_requires_weight_grad = requires_weight_grad
      my_use_atom_chunks = .FALSE.
      IF (PRESENT(use_atom_chunks)) my_use_atom_chunks = use_atom_chunks
      my_compute_diagnostics = .FALSE.
      IF (PRESENT(compute_diagnostics)) my_compute_diagnostics = compute_diagnostics
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
      CALL skala_gpw_feature_release(features, keep_host_arrays=.TRUE.)
      CALL timestop(phase_handle)

      CALL timeset("skala_gpw_layout_cache", phase_handle)
      CALL ensure_layout_cache(pw_grid, particle_set, cell, weights, my_atom_partition, &
                               needs_full_static_tensors=.NOT. my_use_atom_chunks, &
                               needs_local_static=my_use_atom_chunks .AND. &
                               (my_requires_coordinate_grad .OR. my_requires_stress_grad .OR. &
                                my_requires_weight_grad), &
                               cuda_device=cuda_device)
      CALL timestop(phase_handle)
      nflat = skala_state%layout%nflat
      collapse_spin_dynamics = nspins == 1 .AND. my_use_atom_chunks
      ndynamic_local_per_point = ndynamic_per_point
      IF (collapse_spin_dynamics) ndynamic_local_per_point = nrks_dynamic_per_point
      CALL ensure_real_1d(skala_state%scratch%local_dynamic, ndynamic_local_per_point*nflat_local)
      skala_state%scratch%local_dynamic = 0.0_dp
      local_dynamic => skala_state%scratch%local_dynamic

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
      CALL copy_cached_layout(features, my_use_atom_chunks, &
                              my_requires_coordinate_grad .OR. my_requires_stress_grad .OR. &
                              my_requires_weight_grad, &
                              my_requires_stress_grad .OR. &
                              (my_atom_partition == skala_gpw_atom_partition_smooth .AND. &
                               (my_requires_coordinate_grad .OR. my_requires_stress_grad .OR. &
                                my_requires_weight_grad)))
      CALL timestop(phase_handle)

      IF (my_use_atom_chunks) THEN
         CALL timeset("skala_gpw_route_dyn", phase_handle)
         CALL route_atom_chunk_dynamics(features, local_dynamic, pw_grid%para%group, &
                                        collapse_spin_dynamics)
         features%uses_atom_chunks = .TRUE.
         CALL timestop(phase_handle)
      ELSE
         CALL ensure_real_1d(skala_state%scratch%global_dynamic, ndynamic_per_point*skala_state%layout%npoint)
         global_dynamic => skala_state%scratch%global_dynamic
         CALL timeset("skala_gpw_allgatherv", phase_handle)
         CALL pw_grid%para%group%allgatherv(local_dynamic, global_dynamic, &
                                            skala_state%layout%dynamic_counts, &
                                            skala_state%layout%dynamic_displs)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_reorder_dyn", phase_handle)
         CALL ensure_real_2d(features%density, nflat, 2)
         CALL ensure_real_3d(features%grad, nflat, 3, 2)
         CALL ensure_real_2d(features%kin, nflat, 2)

         DO row = 1, nflat
            ipt = skala_state%layout%feature_source_points(row)
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

      IF (my_compute_diagnostics) THEN
         CALL timeset("skala_gpw_feature_sums", phase_handle)
         IF (features%uses_atom_chunks) THEN
            features%electron_count = 0.0_dp
            features%kinetic_integral = 0.0_dp
            features%spin_moment = 0.0_dp
            IF (features%chunk_feature_count > 0) THEN
               IF (features%uses_collapsed_rks_dynamic) THEN
                  features%electron_count = SUM(2.0_dp*features%chunk_density(:, 1)* &
                                                skala_state%layout%chunk_grid_weights)
                  features%kinetic_integral = SUM(2.0_dp*features%chunk_kin(:, 1)* &
                                                  skala_state%layout%chunk_grid_weights)
               ELSE
                  features%electron_count = SUM((features%chunk_density(:, 1) + &
                                                 features%chunk_density(:, 2))* &
                                                skala_state%layout%chunk_grid_weights)
                  features%kinetic_integral = SUM((features%chunk_kin(:, 1) + &
                                                   features%chunk_kin(:, 2))* &
                                                  skala_state%layout%chunk_grid_weights)
                  features%spin_moment = SUM((features%chunk_density(:, 1) - &
                                              features%chunk_density(:, 2))* &
                                             skala_state%layout%chunk_grid_weights)
               END IF
            END IF
            integrals = [features%electron_count, features%kinetic_integral, features%spin_moment]
            CALL pw_grid%para%group%sum(integrals)
            features%electron_count = integrals(1)
            features%kinetic_integral = integrals(2)
            features%spin_moment = integrals(3)
         ELSE
            features%electron_count = SUM((features%density(:, 1) + features%density(:, 2))* &
                                          features%grid_weights)
            features%kinetic_integral = SUM((features%kin(:, 1) + features%kin(:, 2))* &
                                            features%grid_weights)
            features%spin_moment = SUM((features%density(:, 1) - features%density(:, 2))* &
                                       features%grid_weights)
         END IF
         IF (my_use_atom_chunks) THEN
            grid_weight_local = 0.0_dp
            IF (ALLOCATED(skala_state%layout%chunk_grid_weights)) THEN
               grid_weight_local = SUM(skala_state%layout%chunk_grid_weights)
            END IF
            CALL pw_grid%para%group%sum(grid_weight_local)
            features%grid_weight_sum = grid_weight_local
         ELSE
            features%grid_weight_sum = SUM(features%grid_weights)
         END IF
         CALL timestop(phase_handle)
      END IF

      CALL timeset("skala_gpw_tensor_update", phase_handle)
      IF (.NOT. features%uses_atom_chunks .OR. features%chunk_feature_count > 0) THEN
         CALL add_feature_tensors(features, my_requires_grad, my_requires_coordinate_grad, &
                                  my_requires_stress_grad, &
                                  features%uses_atom_chunks, &
                                  requires_weight_grad= &
                                  (my_atom_partition == skala_gpw_atom_partition_smooth .AND. &
                                   (my_requires_coordinate_grad .OR. my_requires_stress_grad .OR. &
                                    my_requires_weight_grad)))
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

      NULLIFY (local_dynamic, global_dynamic)
      CALL timestop(handle)

   END SUBROUTINE skala_gpw_feature_build

! **************************************************************************************************
!> \brief Ensure that static grid-to-atom layout data is cached for the current grid/geometry.
!> \param pw_grid ...
!> \param particle_set ...
!> \param cell ...
!> \param weights ...
!> \param atom_partition ...
!> \param needs_full_static_tensors ...
!> \param needs_local_static ...
!> \param cuda_device ...
! **************************************************************************************************
   SUBROUTINE ensure_layout_cache(pw_grid, particle_set, cell, weights, atom_partition, &
                                  needs_full_static_tensors, needs_local_static, cuda_device)
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: weights
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_partition
      LOGICAL, INTENT(IN)                                :: needs_full_static_tensors, &
                                                            needs_local_static
      INTEGER, INTENT(IN)                                :: cuda_device

      INTEGER                                            :: my_atom_partition, phase_handle
      LOGICAL                                            :: cache_matches

      my_atom_partition = skala_gpw_atom_partition_hard
      IF (PRESENT(atom_partition)) my_atom_partition = atom_partition

      IF (skala_state%layout_mode_known) THEN
         IF (skala_state%layout_mode_full .NEQV. needs_full_static_tensors) THEN
            IF (needs_full_static_tensors) THEN
               CALL cp_abort(__LOCATION__, &
                             "SKALA: the native layout was established in atom-chunk mode, "// &
                             "but a full-grid evaluation was requested. Mixing full-grid and "// &
                             "atom-chunk modes within one run is not supported.")
            ELSE
               CALL cp_abort(__LOCATION__, &
                             "SKALA: the native layout was established in full-grid mode, "// &
                             "but an atom-chunk evaluation was requested. Mixing full-grid and "// &
                             "atom-chunk modes within one run is not supported.")
            END IF
         END IF
      ELSE
         skala_state%layout_mode_known = .TRUE.
         skala_state%layout_mode_full = needs_full_static_tensors
      END IF

      IF (PRESENT(weights)) THEN
         CALL timeset("skala_gpw_layout_match", phase_handle)
         cache_matches = layout_cache_matches(pw_grid, particle_set, cell, weights, &
                                              my_atom_partition, cuda_device, needs_local_static)
         CALL timestop(phase_handle)
      ELSE
         CALL timeset("skala_gpw_layout_match", phase_handle)
         cache_matches = layout_cache_matches(pw_grid, particle_set, cell, &
                                              atom_partition=my_atom_partition, &
                                              cuda_device=cuda_device, &
                                              needs_local_static=needs_local_static)
         CALL timestop(phase_handle)
      END IF

      IF (cache_matches) THEN
         IF (needs_full_static_tensors .AND. .NOT. skala_state%layout%static_tensors_active) THEN
            CALL timeset("skala_gpw_layout_tensors", phase_handle)
            CALL build_full_static_layout_tensors(skala_state%layout)
            CALL timestop(phase_handle)
         END IF
         RETURN
      END IF

      CALL timeset("skala_gpw_layout_rebuild", phase_handle)
      IF (PRESENT(weights)) THEN
         CALL rebuild_layout_cache(pw_grid, particle_set, cell, weights, my_atom_partition, &
                                   needs_full_static_tensors, needs_local_static, cuda_device)
      ELSE
         CALL rebuild_layout_cache(pw_grid, particle_set, cell, &
                                   atom_partition=my_atom_partition, &
                                   needs_full_static_tensors=needs_full_static_tensors, &
                                   needs_local_static=needs_local_static, &
                                   cuda_device=cuda_device)
      END IF
      CALL timestop(phase_handle)

   END SUBROUTINE ensure_layout_cache

! **************************************************************************************************
!> \brief Check whether the geometry/periodicity/communicator part of the layout key matches.
!> \param cached_hmat ...
!> \param cached_perd ...
!> \param cached_comm ...
!> \param cell_hmat ...
!> \param cell_perd ...
!> \param comm ...
!> \return ...
! **************************************************************************************************
   FUNCTION skala_gpw_layout_metadata_matches(cached_hmat, cached_perd, cached_comm, &
                                              cell_hmat, cell_perd, comm) RESULT(matches)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: cached_hmat
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cached_perd
      INTEGER, INTENT(IN)                                :: cached_comm
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: cell_hmat
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell_perd
      INTEGER, INTENT(IN)                                :: comm
      LOGICAL                                            :: matches

      matches = .NOT. ANY(ABS(cached_hmat - cell_hmat) > layout_tol)
      matches = matches .AND. ALL(cached_perd == cell_perd)
      matches = matches .AND. cached_comm == comm

   END FUNCTION skala_gpw_layout_metadata_matches

! **************************************************************************************************
!> \brief Check whether the current static layout cache can be reused.
!>
!> The decision must be identical on every MPI rank: a reuse skips the layout
!> rebuild collectives, so a rank-local match while another rank mismatches
!> would deadlock. The local test (which includes the distributed integration
!> weights) is therefore reduced with a logical AND across the communicator.
!> \param pw_grid ...
!> \param particle_set ...
!> \param cell ...
!> \param weights ...
!> \param atom_partition ...
!> \param cuda_device ...
!> \param needs_local_static ...
!> \return ...
! **************************************************************************************************
   FUNCTION layout_cache_matches(pw_grid, particle_set, cell, weights, atom_partition, &
                                 cuda_device, needs_local_static) RESULT(matches)
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: weights
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_partition
      INTEGER, INTENT(IN)                                :: cuda_device
      LOGICAL, INTENT(IN), OPTIONAL                      :: needs_local_static
      LOGICAL                                            :: matches

      INTEGER                                            :: hit

      hit = MERGE(1, 0, layout_cache_matches_local(pw_grid, particle_set, cell, weights, &
                                                   atom_partition, cuda_device, needs_local_static))
      CALL pw_grid%para%group%min(hit)
      matches = (hit == 1)

   END FUNCTION layout_cache_matches

! **************************************************************************************************
!> \brief Rank-local part of layout_cache_matches (no collectives).
!> \param pw_grid ...
!> \param particle_set ...
!> \param cell ...
!> \param weights ...
!> \param atom_partition ...
!> \param cuda_device ...
!> \param needs_local_static ...
!> \return ...
! **************************************************************************************************
   FUNCTION layout_cache_matches_local(pw_grid, particle_set, cell, weights, atom_partition, &
                                       cuda_device, needs_local_static) RESULT(matches)
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: weights
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_partition
      INTEGER, INTENT(IN)                                :: cuda_device
      LOGICAL, INTENT(IN), OPTIONAL                      :: needs_local_static
      LOGICAL                                            :: matches

      INTEGER                                            :: iatom, my_atom_partition
      LOGICAL                                            :: weights_match

      my_atom_partition = skala_gpw_atom_partition_hard
      IF (PRESENT(atom_partition)) my_atom_partition = atom_partition
      matches = .FALSE.
      IF (.NOT. skala_state%layout%active) RETURN
      IF (PRESENT(needs_local_static)) THEN
         IF (needs_local_static .AND. .NOT. ALLOCATED(skala_state%layout%local_feature_grid_coords)) RETURN
      END IF
      IF (skala_state%layout%cuda_device /= cuda_device) RETURN
      IF (skala_state%layout%atom_partition /= my_atom_partition) RETURN
      IF (skala_state%layout%natom /= SIZE(particle_set)) RETURN
      IF (skala_state%layout%nproc /= pw_grid%para%group%num_pe) RETURN
      IF (ANY(skala_state%layout%bo /= pw_grid%bounds_local)) RETURN
      IF (ANY(skala_state%layout%bounds /= pw_grid%bounds)) RETURN
      IF (ANY(skala_state%layout%npts /= pw_grid%npts)) RETURN
      IF (ANY(ABS(skala_state%layout%dh - pw_grid%dh) > layout_tol)) RETURN
      IF (.NOT. skala_gpw_layout_metadata_matches(skala_state%layout%cell_hmat, &
                                                  skala_state%layout%cell_perd, &
                                                  skala_state%layout%comm_handle, cell%hmat, &
                                                  cell%perd, pw_grid%para%group%get_handle())) RETURN
      IF (.NOT. ALLOCATED(skala_state%layout%atom_coords)) RETURN

      DO iatom = 1, SIZE(particle_set)
         IF (ANY(ABS(skala_state%layout%atom_coords(:, iatom) - particle_set(iatom)%r) > layout_tol)) RETURN
      END DO

      IF (PRESENT(weights)) THEN
         weights_match = layout_weights_match(weights)
      ELSE
         weights_match = layout_weights_match()
      END IF
      IF (.NOT. weights_match) RETURN

      matches = .TRUE.

   END FUNCTION layout_cache_matches_local

! **************************************************************************************************
!> \brief Check whether current optional integration weights match the cached static tensors.
!> \param weights ...
!> \return ...
! **************************************************************************************************
   FUNCTION layout_weights_match(weights) RESULT(matches)
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: weights
      LOGICAL                                            :: matches

      LOGICAL                                            :: has_weights
      REAL(KIND=dp)                                      :: weight_sum, weight_sumsq

      matches = .FALSE.

      has_weights = .FALSE.
      IF (PRESENT(weights)) THEN
         IF (ASSOCIATED(weights)) has_weights = .TRUE.
      END IF
      IF (skala_state%layout%has_weights .NEQV. has_weights) RETURN

      IF (PRESENT(weights)) THEN
         CALL weights_signature(weights, has_weights, weight_sum, weight_sumsq)
      ELSE
         CALL weights_signature(has_weights=has_weights, weight_sum=weight_sum, &
                                weight_sumsq=weight_sumsq)
      END IF
      IF (ABS(skala_state%layout%weight_sum - weight_sum) > layout_tol) RETURN
      IF (ABS(skala_state%layout%weight_sumsq - weight_sumsq) > layout_tol) RETURN

      matches = .TRUE.

   END FUNCTION layout_weights_match

! **************************************************************************************************
!> \brief Build the static SKALA layout cache.
!> \param pw_grid ...
!> \param particle_set ...
!> \param cell ...
!> \param weights ...
!> \param atom_partition ...
!> \param needs_full_static_tensors ...
!> \param needs_local_static ...
!> \param cuda_device ...
! **************************************************************************************************
   SUBROUTINE rebuild_layout_cache(pw_grid, particle_set, cell, weights, atom_partition, &
                                   needs_full_static_tensors, needs_local_static, cuda_device)
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER            :: weights
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_partition
      LOGICAL, INTENT(IN)                                :: needs_full_static_tensors, &
                                                            needs_local_static
      INTEGER, INTENT(IN)                                :: cuda_device

      INTEGER :: chunk_row, dest, dest_base, dist_pos, feature_local, feature_slot, i, iatom, ipt, &
         j, jatom, k, local_feature, local_row, max_grid_size, max_local_features, &
         my_atom_partition, natom, nfeature_local, nflat, nflat_local, npoint, nproc, nx, ny, &
         owner, pe, pe_index, phase_handle, row, source_global, source_local, static_base
      INTEGER, ALLOCATABLE, DIMENSION(:) :: atom_chunk_owner, atom_offset, atom_position, &
         chunk_atom_begin, chunk_atom_end, chunk_feature_counts, chunk_feature_displs, cursor, &
         dist_counts, dist_cursor, dist_displs, dist_recv_counts, dist_recv_displs, &
         dist_send_counts, dist_send_displs, feature_counts, feature_displs, global_owner, &
         global_source_points, local_feature_counts_tmp, local_owner, local_source_global, &
         local_source_points, point_counts, point_displs, static_counts, static_displs
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: has_weights, nonperiodic
      REAL(KIND=dp)                                      :: atomic_base_weight, base_weight, &
                                                            included_sum, partition_weight, &
                                                            weight_sum, weight_sumsq
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dist_recv_buffer, dist_send_buffer, &
                                                            distances, global_static, &
                                                            local_static, partition_weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_coords_pbc, pair_distances, &
                                                            partition_atom_coords
      REAL(KIND=dp), DIMENSION(3)                        :: grid_point

      CALL release_layout_cache(skala_state%layout)

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
                local_static(nstatic_per_point*max_local_features), &
                local_feature_counts_tmp(nflat_local), feature_counts(nproc), &
                feature_displs(nproc), static_counts(nproc), static_displs(nproc), &
                atom_coords_pbc(3, natom))
      IF (needs_full_static_tensors) THEN
         ALLOCATE (local_source_points(max_local_features), point_counts(nproc), &
                   point_displs(nproc))
      END IF
      local_feature_counts_tmp = 0
      DO iatom = 1, natom
         atom_coords_pbc(:, iatom) = pbc(particle_set(iatom)%r, cell, positive_range=.TRUE.)
      END DO
      nonperiodic = ALL(cell%perd == 0)
      IF (my_atom_partition == skala_gpw_atom_partition_smooth .AND. nonperiodic) THEN
         ALLOCATE (pair_distances(natom, natom))
         pair_distances = 0.0_dp
         DO iatom = 1, natom - 1
            DO jatom = iatom + 1, natom
               pair_distances(iatom, jatom) = &
                  SQRT(SUM((atom_coords_pbc(:, iatom) - atom_coords_pbc(:, jatom))**2))
               pair_distances(jatom, iatom) = pair_distances(iatom, jatom)
            END DO
         END DO
      ELSE
         ALLOCATE (pair_distances(0, 0))
      END IF

      CALL timeset("skala_gpw_layout_local", phase_handle)
      nx = bo(2, 1) - bo(1, 1) + 1
      ny = bo(2, 2) - bo(1, 2) + 1
      IF (my_atom_partition == skala_gpw_atom_partition_hard) THEN
         nfeature_local = nflat_local
!$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(3) &
!$OMP SHARED(atom_coords_pbc, bo, skala_state, cell, local_feature_counts_tmp, local_owner, &
!$OMP        local_source_points, local_static, needs_full_static_tensors, nx, ny, pw_grid, &
!$OMP        weights) &
!$OMP PRIVATE(atomic_base_weight, base_weight, grid_point, i, j, k, local_row, owner, static_base)
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  local_row = i - bo(1, 1) + 1 + &
                              nx*(j - bo(1, 2) + ny*(k - bo(1, 3)))
                  grid_point = grid_coordinate(pw_grid, [i, j, k])
                  atomic_base_weight = pw_grid%dvol
                  base_weight = pw_grid%dvol
                  IF (PRESENT(weights)) THEN
                     IF (ASSOCIATED(weights)) base_weight = base_weight*weights%array(i, j, k)
                  END IF
                  owner = nearest_atom(grid_point, atom_coords_pbc, cell)
                  local_feature_counts_tmp(local_row) = 1
                  local_owner(local_row) = owner
                  IF (needs_full_static_tensors) local_source_points(local_row) = local_row
                  static_base = nstatic_per_point*(local_row - 1)
                  local_static(static_base + 1:static_base + 3) = grid_point
                  local_static(static_base + 4) = base_weight
                  local_static(static_base + 5) = atomic_base_weight
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
      ELSE
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(atom_coords_pbc, bo, skala_state, cell, local_feature_counts_tmp, local_owner, &
!$OMP        local_source_points, local_static, natom, needs_full_static_tensors, nonperiodic, &
!$OMP        nx, ny, pair_distances, pw_grid, weights) &
!$OMP PRIVATE(atomic_base_weight, base_weight, distances, feature_slot, grid_point, i, iatom, &
!$OMP         included_sum, j, k, local_feature, local_row, owner, partition_weight, &
!$OMP         partition_atom_coords, partition_weights, static_base)
         ! Keep point-local scratch outside the hot loop and private to each thread.
         ALLOCATE (distances(natom), partition_atom_coords(3, natom), partition_weights(natom))
!$OMP DO COLLAPSE(3)
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  local_row = i - bo(1, 1) + 1 + &
                              nx*(j - bo(1, 2) + ny*(k - bo(1, 3)))
                  grid_point = grid_coordinate(pw_grid, [i, j, k])
                  atomic_base_weight = pw_grid%dvol
                  base_weight = pw_grid%dvol
                  IF (PRESENT(weights)) THEN
                     IF (ASSOCIATED(weights)) base_weight = base_weight*weights%array(i, j, k)
                  END IF
                  IF (nonperiodic) THEN
                     CALL smooth_atom_partition(grid_point, atom_coords_pbc, cell, &
                                                partition_weights, partition_atom_coords, distances, &
                                                pair_distances)
                  ELSE
                     CALL smooth_atom_partition(grid_point, atom_coords_pbc, cell, &
                                                partition_weights, partition_atom_coords, distances)
                  END IF
                  included_sum = SUM(partition_weights, MASK=partition_weights > smooth_partition_eps)
                  IF (included_sum <= 0.0_dp) THEN
                     owner = nearest_atom(grid_point, atom_coords_pbc, cell)
                     partition_weights = 0.0_dp
                     partition_weights(owner) = 1.0_dp
                     included_sum = 1.0_dp
                  END IF
                  local_feature = 0
                  DO iatom = 1, natom
                     IF (partition_weights(iatom) <= smooth_partition_eps) CYCLE
                     partition_weight = partition_weights(iatom)/included_sum
                     local_feature = local_feature + 1
                     feature_slot = natom*(local_row - 1) + local_feature
                     local_owner(feature_slot) = iatom
                     IF (needs_full_static_tensors) local_source_points(feature_slot) = local_row
                     static_base = nstatic_per_point*(feature_slot - 1)
                     local_static(static_base + 1:static_base + 3) = grid_point
                     local_static(static_base + 4) = base_weight*partition_weight
                     local_static(static_base + 5) = atomic_base_weight* &
                                                     smooth_partition_atomic_weight_scale( &
                                                     partition_weight)
                  END DO
                  local_feature_counts_tmp(local_row) = local_feature
               END DO
            END DO
         END DO
!$OMP END DO
         DEALLOCATE (distances, partition_atom_coords, partition_weights)
!$OMP END PARALLEL

         ! Compact the conflict-free row slots while preserving the serial row/atom order.
         nfeature_local = 0
         DO local_row = 1, nflat_local
            CPASSERT(local_feature_counts_tmp(local_row) > 0)
            DO local_feature = 1, local_feature_counts_tmp(local_row)
               nfeature_local = nfeature_local + 1
               feature_slot = natom*(local_row - 1) + local_feature
               IF (feature_slot == nfeature_local) CYCLE
               local_owner(nfeature_local) = local_owner(feature_slot)
               IF (needs_full_static_tensors) local_source_points(nfeature_local) = &
                  local_source_points(feature_slot)
               local_static(nstatic_per_point*(nfeature_local - 1) + 1: &
                            nstatic_per_point*nfeature_local) = &
                  local_static(nstatic_per_point*(feature_slot - 1) + 1: &
                               nstatic_per_point*feature_slot)
            END DO
         END DO
      END IF
      CALL timestop(phase_handle)

      ! SKALA groups all grid points by atom. This ordering is static while the
      ! grid, cell, atom positions, and optional integration weights are unchanged.
      CALL timeset("skala_gpw_layout_gather", phase_handle)
      npoint = 0
      IF (needs_full_static_tensors) THEN
         point_counts = 0
         point_displs = 0
         CALL pw_grid%para%group%allgather(nflat_local, point_counts)
         CALL counts_to_displs(point_counts, point_displs)
         npoint = SUM(point_counts)
      END IF
      CALL pw_grid%para%group%allgather(nfeature_local, feature_counts)
      CALL counts_to_displs(feature_counts, feature_displs)
      DO pe = 1, nproc
         static_counts(pe) = nstatic_per_point*feature_counts(pe)
         static_displs(pe) = nstatic_per_point*feature_displs(pe)
      END DO
      nflat = SUM(feature_counts)
      IF (needs_full_static_tensors) THEN
         ALLOCATE (global_owner(nflat), global_source_points(nflat), &
                   global_static(nstatic_per_point*nflat), local_source_global(nfeature_local))
         DO feature_local = 1, nfeature_local
            local_source_global(feature_local) = point_displs(pe_index) + &
                                                 local_source_points(feature_local)
         END DO
         CALL pw_grid%para%group%allgatherv(local_owner(1:nfeature_local), global_owner, &
                                            feature_counts, feature_displs)
         CALL pw_grid%para%group%allgatherv(local_source_global, global_source_points, &
                                            feature_counts, feature_displs)
         CALL pw_grid%para%group%allgatherv(local_static(1:nstatic_per_point*nfeature_local), &
                                            global_static, static_counts, static_displs)
      ELSE
         ALLOCATE (global_owner(0), global_source_points(0), global_static(0))
      END IF
      CALL timestop(phase_handle)

      ALLOCATE (chunk_feature_counts(nproc), chunk_feature_displs(nproc), &
                skala_state%layout%route_point_recv_counts(nproc), &
                skala_state%layout%route_point_recv_displs(nproc), &
                skala_state%layout%route_point_send_counts(nproc), &
                skala_state%layout%route_point_send_displs(nproc), &
                skala_state%layout%atomic_grid_sizes(natom), &
                skala_state%layout%local_feature_offsets(nflat_local + 1), &
                skala_state%layout%local_feature_rows(nfeature_local), &
                atom_offset(natom + 1), chunk_atom_begin(nproc), chunk_atom_end(nproc))
      skala_state%layout%atomic_grid_sizes = 0_int_8
      skala_state%layout%local_feature_offsets(1) = 1
      DO local_row = 1, nflat_local
         skala_state%layout%local_feature_offsets(local_row + 1) = &
            skala_state%layout%local_feature_offsets(local_row) + &
            local_feature_counts_tmp(local_row)
      END DO

      IF (needs_local_static) THEN
         ALLOCATE (skala_state%layout%local_feature_grid_coords(3, nfeature_local), &
                   skala_state%layout%local_feature_grid_weights(nfeature_local), &
                   skala_state%layout%local_feature_atomic_grid_weights(nfeature_local))
         skala_state%layout%local_feature_grid_coords = 0.0_dp
         skala_state%layout%local_feature_grid_weights = 0.0_dp
         skala_state%layout%local_feature_atomic_grid_weights = 0.0_dp
      END IF

      CALL timeset("skala_gpw_layout_atom_sort", phase_handle)
      IF (needs_full_static_tensors) THEN
         ALLOCATE (skala_state%layout%dynamic_counts(nproc), skala_state%layout%dynamic_displs(nproc), &
                   skala_state%layout%feature_source_points(nflat), &
                   atom_position(natom), cursor(nflat_local))
         skala_state%layout%dynamic_counts(:) = ndynamic_per_point*point_counts
         skala_state%layout%dynamic_displs(:) = ndynamic_per_point*point_displs
         cursor(:) = skala_state%layout%local_feature_offsets(1:nflat_local)
         DO ipt = 1, nflat
            skala_state%layout%atomic_grid_sizes(global_owner(ipt)) = &
               skala_state%layout%atomic_grid_sizes(global_owner(ipt)) + 1_int_8
         END DO
         atom_offset(1) = 1
         DO iatom = 1, natom
            atom_offset(iatom + 1) = atom_offset(iatom) + INT(skala_state%layout%atomic_grid_sizes(iatom))
         END DO
         DO iatom = 1, natom
            atom_position(iatom) = atom_offset(iatom)
         END DO
      ELSE
         BLOCK
            INTEGER, ALLOCATABLE, DIMENSION(:) :: atom_before, atom_cursor, atom_inclusive, &
                                                  atom_total, local_atom_counts
            INTEGER, ALLOCATABLE, DIMENSION(:, :) :: inclusive_scan, local_scan

            ALLOCATE (local_atom_counts(natom), atom_total(natom), atom_inclusive(natom), &
                      atom_before(natom), atom_cursor(natom), local_scan(natom, 1), &
                      inclusive_scan(natom, 1))
            local_atom_counts = 0
            DO feature_local = 1, nfeature_local
               local_atom_counts(local_owner(feature_local)) = &
                  local_atom_counts(local_owner(feature_local)) + 1
            END DO
            atom_total(:) = local_atom_counts
            CALL pw_grid%para%group%sum(atom_total)
            local_scan(:, 1) = local_atom_counts
            CALL pw_grid%para%group%sum_partial(local_scan, inclusive_scan)
            atom_inclusive(:) = inclusive_scan(:, 1)
            atom_before(:) = atom_inclusive - local_atom_counts
            DO iatom = 1, natom
               skala_state%layout%atomic_grid_sizes(iatom) = INT(atom_total(iatom), KIND=int_8)
            END DO
            atom_offset(1) = 1
            DO iatom = 1, natom
               atom_offset(iatom + 1) = atom_offset(iatom) + &
                                        INT(skala_state%layout%atomic_grid_sizes(iatom))
            END DO
            atom_cursor = 0
            DO feature_local = 1, nfeature_local
               owner = local_owner(feature_local)
               row = atom_offset(owner) + atom_before(owner) + atom_cursor(owner)
               atom_cursor(owner) = atom_cursor(owner) + 1
               skala_state%layout%local_feature_rows(feature_local) = row
               IF (needs_local_static) THEN
                  static_base = nstatic_per_point*(feature_local - 1)
                  skala_state%layout%local_feature_grid_coords(:, feature_local) = &
                     local_static(static_base + 1:static_base + 3)
                  skala_state%layout%local_feature_grid_weights(feature_local) = &
                     local_static(static_base + 4)
                  skala_state%layout%local_feature_atomic_grid_weights(feature_local) = &
                     local_static(static_base + 5)
               END IF
            END DO
            DEALLOCATE (local_atom_counts, atom_total, atom_inclusive, atom_before, atom_cursor, &
                        local_scan, inclusive_scan)
         END BLOCK
      END IF
      max_grid_size = MAXVAL(INT(skala_state%layout%atomic_grid_sizes))
      IF (.NOT. needs_full_static_tensors) THEN
         CALL build_atom_chunks(skala_state%layout%atomic_grid_sizes, atom_offset, nproc, &
                                chunk_atom_begin, chunk_atom_end, &
                                chunk_feature_counts, chunk_feature_displs)
         skala_state%layout%chunk_atom_begin = chunk_atom_begin(pe_index)
         skala_state%layout%chunk_atom_end = chunk_atom_end(pe_index)
         skala_state%layout%chunk_feature_begin = chunk_feature_displs(pe_index) + 1
         skala_state%layout%chunk_feature_count = chunk_feature_counts(pe_index)
         skala_state%layout%chunk_natom = skala_state%layout%chunk_atom_end - &
                                          skala_state%layout%chunk_atom_begin + 1

         ALLOCATE (atom_chunk_owner(natom))
         atom_chunk_owner = 0
         DO pe = 1, nproc
            DO iatom = chunk_atom_begin(pe), chunk_atom_end(pe)
               IF (iatom >= 1 .AND. iatom <= natom) atom_chunk_owner(iatom) = pe
            END DO
         END DO
      END IF

      ALLOCATE (skala_state%layout%coarse_0_atomic_coords(3, natom), &
                skala_state%layout%atom_coords(3, natom))
      IF (needs_full_static_tensors) THEN
         CALL allocate_empty_bound_shape(skala_state%layout%atomic_grid_size_bound_shape, max_grid_size)
         ALLOCATE (skala_state%layout%grid_coords(3, nflat), skala_state%layout%grid_weights(nflat), &
                   skala_state%layout%atomic_grid_weights(nflat))
         skala_state%layout%grid_coords = 0.0_dp
         skala_state%layout%grid_weights = 0.0_dp
         skala_state%layout%atomic_grid_weights = 0.0_dp
      END IF

      DO iatom = 1, natom
         skala_state%layout%atom_coords(:, iatom) = particle_set(iatom)%r
         skala_state%layout%coarse_0_atomic_coords(:, iatom) = atom_coords_pbc(:, iatom)
      END DO

      IF (needs_full_static_tensors) THEN
         DO ipt = 1, nflat
            owner = global_owner(ipt)
            row = atom_position(owner)
            atom_position(owner) = atom_position(owner) + 1
            source_global = global_source_points(ipt)
            skala_state%layout%feature_source_points(row) = source_global
            static_base = nstatic_per_point*(ipt - 1)
            skala_state%layout%grid_coords(:, row) = global_static(static_base + 1:static_base + 3)
            skala_state%layout%grid_weights(row) = global_static(static_base + 4)
            skala_state%layout%atomic_grid_weights(row) = global_static(static_base + 5)
            source_local = source_global - point_displs(pe_index)
            IF (source_local >= 1 .AND. source_local <= nflat_local) THEN
               feature_local = cursor(source_local)
               cursor(source_local) = cursor(source_local) + 1
               skala_state%layout%local_feature_rows(feature_local) = row
               IF (needs_local_static) THEN
                  static_base = nstatic_per_point*(feature_local - 1)
                  skala_state%layout%local_feature_grid_coords(:, feature_local) = &
                     local_static(static_base + 1:static_base + 3)
                  skala_state%layout%local_feature_grid_weights(feature_local) = &
                     local_static(static_base + 4)
                  skala_state%layout%local_feature_atomic_grid_weights(feature_local) = &
                     local_static(static_base + 5)
               END IF
            END IF
         END DO
         CPASSERT(ALL(skala_state%layout%local_feature_rows > 0))
      ELSE
         CPASSERT(ALL(skala_state%layout%local_feature_rows > 0))
      END IF

      IF (.NOT. needs_full_static_tensors) THEN
         IF (skala_state%layout%chunk_feature_count > 0) THEN
            ALLOCATE (skala_state%layout%chunk_grid_coords(3, skala_state%layout%chunk_feature_count), &
                      skala_state%layout%chunk_grid_weights(skala_state%layout%chunk_feature_count), &
                      skala_state%layout%chunk_atomic_grid_weights(skala_state%layout%chunk_feature_count))
            skala_state%layout%chunk_grid_coords = 0.0_dp
            skala_state%layout%chunk_grid_weights = 0.0_dp
            skala_state%layout%chunk_atomic_grid_weights = 0.0_dp
         END IF
         ALLOCATE (dist_counts(nproc), dist_displs(nproc), &
                   dist_recv_counts(nproc), dist_recv_displs(nproc), dist_cursor(nproc), &
                   dist_send_counts(nproc), dist_send_displs(nproc))
         dist_counts = 0
         DO feature_local = 1, nfeature_local
            dest = atom_chunk_owner(local_owner(feature_local))
            CPASSERT(dest > 0)
            dist_counts(dest) = dist_counts(dest) + 1
         END DO
         CALL counts_to_displs(dist_counts, dist_displs)
         CALL pw_grid%para%group%alltoall(dist_counts, dist_recv_counts, 1)
         CALL counts_to_displs(dist_recv_counts, dist_recv_displs)
         ALLOCATE (dist_send_buffer(MAX(1, 6*nfeature_local)), &
                   dist_recv_buffer(MAX(1, 6*SUM(dist_recv_counts))))
         dist_send_buffer = 0.0_dp
         dist_recv_buffer = 0.0_dp
         dist_cursor(:) = dist_displs + 1
         DO feature_local = 1, nfeature_local
            dest = atom_chunk_owner(local_owner(feature_local))
            dist_pos = dist_cursor(dest)
            dist_cursor(dest) = dist_cursor(dest) + 1
            static_base = nstatic_per_point*(feature_local - 1)
            dest_base = 6*(dist_pos - 1)
            dist_send_buffer(dest_base + 1) = &
               REAL(skala_state%layout%local_feature_rows(feature_local), KIND=dp)
            dist_send_buffer(dest_base + 2:dest_base + 4) = &
               local_static(static_base + 1:static_base + 3)
            dist_send_buffer(dest_base + 5) = local_static(static_base + 4)
            dist_send_buffer(dest_base + 6) = local_static(static_base + 5)
         END DO
         dist_send_counts(:) = 6*dist_counts
         dist_send_displs(:) = 6*dist_displs
         dist_recv_counts(:) = 6*dist_recv_counts
         dist_recv_displs(:) = 6*dist_recv_displs
         CALL pw_grid%para%group%alltoall(dist_send_buffer, dist_send_counts, dist_send_displs, &
                                          dist_recv_buffer, dist_recv_counts, dist_recv_displs)
         IF (skala_state%layout%chunk_feature_count > 0) THEN
            DO i = 1, SUM(dist_recv_counts)/6
               dest_base = 6*(i - 1)
               row = NINT(dist_recv_buffer(dest_base + 1))
               chunk_row = row - skala_state%layout%chunk_feature_begin + 1
               CPASSERT(chunk_row >= 1)
               CPASSERT(chunk_row <= skala_state%layout%chunk_feature_count)
               skala_state%layout%chunk_grid_coords(:, chunk_row) = &
                  dist_recv_buffer(dest_base + 2:dest_base + 4)
               skala_state%layout%chunk_grid_weights(chunk_row) = dist_recv_buffer(dest_base + 5)
               skala_state%layout%chunk_atomic_grid_weights(chunk_row) = dist_recv_buffer(dest_base + 6)
            END DO
         END IF
         DEALLOCATE (dist_counts, dist_displs, dist_recv_counts, &
                     dist_recv_displs, dist_cursor, dist_send_counts, dist_send_displs, &
                     dist_send_buffer, dist_recv_buffer)
      END IF

      CALL timestop(phase_handle)
      IF (.NOT. needs_full_static_tensors) THEN
         CALL timeset("skala_gpw_layout_chunk_routes", phase_handle)
         CALL build_atom_chunk_routes(skala_state%layout, skala_state%layout%local_feature_rows, &
                                      local_owner, atom_chunk_owner, pw_grid%para%group)
         CALL build_atom_chunk_layout(skala_state%layout)
         CALL timestop(phase_handle)
      END IF

      skala_state%layout%natom = natom
      skala_state%layout%nflat = nflat
      skala_state%layout%nflat_local = nflat_local
      skala_state%layout%npoint = npoint
      skala_state%layout%nproc = nproc
      skala_state%layout%atom_partition = my_atom_partition
      skala_state%layout%cuda_device = cuda_device
      skala_state%layout%bo = bo
      skala_state%layout%bounds = pw_grid%bounds
      skala_state%layout%npts = pw_grid%npts
      skala_state%layout%dh = pw_grid%dh
      skala_state%layout%cell_hmat = cell%hmat
      skala_state%layout%cell_perd = cell%perd
      skala_state%layout%comm_handle = pw_grid%para%group%get_handle()
      skala_state%layout%weight_sum = weight_sum
      skala_state%layout%weight_sumsq = weight_sumsq
      skala_state%layout%has_weights = has_weights
      CALL timeset("skala_gpw_layout_tensors", phase_handle)
      IF (needs_full_static_tensors) THEN
         CALL build_full_static_layout_tensors(skala_state%layout)
      ELSE IF (skala_state%layout%chunk_feature_count > 0) THEN
         CALL build_chunk_static_layout_tensors(skala_state%layout)
      END IF
      CALL timestop(phase_handle)
      skala_state%layout%active = .TRUE.

      DEALLOCATE (atom_coords_pbc, atom_offset, &
                  chunk_atom_begin, chunk_atom_end, feature_counts, feature_displs, &
                  global_owner, global_source_points, global_static, local_feature_counts_tmp, &
                  local_owner, local_static, pair_distances, static_counts, static_displs)
      IF (ALLOCATED(local_source_points)) DEALLOCATE (local_source_points)
      IF (ALLOCATED(point_counts)) DEALLOCATE (point_counts)
      IF (ALLOCATED(point_displs)) DEALLOCATE (point_displs)
      IF (ALLOCATED(atom_position)) DEALLOCATE (atom_position)
      IF (ALLOCATED(cursor)) DEALLOCATE (cursor)
      IF (ALLOCATED(local_source_global)) DEALLOCATE (local_source_global)

   END SUBROUTINE rebuild_layout_cache

! **************************************************************************************************
!> \brief Build cached full-layout Torch tensors for static SKALA inputs.
!> \param cache ...
! **************************************************************************************************
   SUBROUTINE build_full_static_layout_tensors(cache)
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

      CALL torch_dict_create(cache%static_inputs)
      CALL torch_dict_insert(cache%static_inputs, skala_key_grid_coords, cache%grid_coords_t)
      CALL torch_dict_insert(cache%static_inputs, skala_key_grid_weights, cache%grid_weights_t)
      CALL torch_dict_insert(cache%static_inputs, skala_key_atomic_grid_weights, &
                             cache%atomic_grid_weights_t)
      CALL torch_dict_insert(cache%static_inputs, skala_key_atomic_grid_sizes, &
                             cache%atomic_grid_sizes_t)
      CALL torch_dict_insert(cache%static_inputs, skala_key_atomic_grid_size_bound_shape, &
                             cache%atomic_grid_size_bound_shape_t)
      cache%static_tensors_active = .TRUE.

   END SUBROUTINE build_full_static_layout_tensors

! **************************************************************************************************
!> \brief Build cached atom-chunk Torch tensors for static SKALA inputs.
!> \param cache ...
! **************************************************************************************************
   SUBROUTINE build_chunk_static_layout_tensors(cache)
      TYPE(skala_gpw_layout_cache_type), INTENT(INOUT)   :: cache

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

      CALL torch_dict_create(cache%chunk_static_inputs)
      CALL torch_dict_insert(cache%chunk_static_inputs, skala_key_grid_coords, &
                             cache%chunk_grid_coords_t)
      CALL torch_dict_insert(cache%chunk_static_inputs, skala_key_grid_weights, &
                             cache%chunk_grid_weights_t)
      CALL torch_dict_insert(cache%chunk_static_inputs, skala_key_atomic_grid_weights, &
                             cache%chunk_atomic_grid_weights_t)
      CALL torch_dict_insert(cache%chunk_static_inputs, skala_key_atomic_grid_sizes, &
                             cache%chunk_atomic_grid_sizes_t)
      CALL torch_dict_insert(cache%chunk_static_inputs, skala_key_atomic_grid_size_bound_shape, &
                             cache%chunk_atomic_grid_size_bound_shape_t)
      cache%chunk_static_tensors_active = .TRUE.

   END SUBROUTINE build_chunk_static_layout_tensors

! **************************************************************************************************
!> \brief Ensure the chunk static Torch tensors expose the requested differentiability.
!> \param grad_mode ...
! **************************************************************************************************
   SUBROUTINE ensure_chunk_static_tensor_grad(grad_mode)
      INTEGER, INTENT(IN)                                :: grad_mode

      CPASSERT(skala_state%layout%chunk_static_tensors_active)
      CALL torch_tensor_set_requires_grad(skala_state%layout%chunk_grid_coords_t, grad_mode >= 2)
      CALL torch_tensor_set_requires_grad(skala_state%layout%chunk_grid_weights_t, grad_mode >= 1)
      CALL torch_tensor_set_requires_grad(skala_state%layout%chunk_atomic_grid_weights_t, &
                                          grad_mode >= 2)
      CALL torch_tensor_set_requires_grad(skala_state%layout%chunk_coarse_0_atomic_coords_t, &
                                          grad_mode >= 2)

   END SUBROUTINE ensure_chunk_static_tensor_grad

! **************************************************************************************************
!> \brief Ensure the full-grid static Torch tensors expose the requested differentiability,
!>        reusing the cached leaves instead of re-uploading them every evaluation.
!> \param grad_mode ...
! **************************************************************************************************
   SUBROUTINE ensure_full_static_tensor_grad(grad_mode)
      INTEGER, INTENT(IN)                                :: grad_mode

      CPASSERT(skala_state%layout%static_tensors_active)
      CALL torch_tensor_set_requires_grad(skala_state%layout%grid_coords_t, grad_mode >= 2)
      CALL torch_tensor_set_requires_grad(skala_state%layout%grid_weights_t, grad_mode >= 1)
      CALL torch_tensor_set_requires_grad(skala_state%layout%atomic_grid_weights_t, grad_mode >= 2)
   END SUBROUTINE ensure_full_static_tensor_grad

! **************************************************************************************************
!> \brief Copy static cached layout arrays into a feature bundle.
!> \param features ...
!> \param use_atom_chunk_routing ...
!> \param needs_coordinate_array ...
!> \param needs_grid_coordinate_array ...
! **************************************************************************************************
   SUBROUTINE copy_cached_layout(features, use_atom_chunk_routing, needs_coordinate_array, &
                                 needs_grid_coordinate_array)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      LOGICAL, INTENT(IN)                                :: use_atom_chunk_routing, &
                                                            needs_coordinate_array, &
                                                            needs_grid_coordinate_array

      CPASSERT(skala_state%layout%active)

      IF (use_atom_chunk_routing) THEN
         CALL copy_routing_layout(features, needs_coordinate_array)
      ELSE
         CALL copy_full_layout(features, needs_grid_coordinate_array)
      END IF

      features%nflat = skala_state%layout%nflat
      features%nflat_local = skala_state%layout%nflat_local
      features%chunk_atom_begin = skala_state%layout%chunk_atom_begin
      features%chunk_feature_count = skala_state%layout%chunk_feature_count
      features%chunk_natom = skala_state%layout%chunk_natom
      features%atom_partition = skala_state%layout%atom_partition
      CALL ensure_i8_1d(features%atomic_grid_sizes, skala_state%layout%natom)
      features%atomic_grid_sizes(:) = skala_state%layout%atomic_grid_sizes
      IF (needs_coordinate_array) THEN
         CALL ensure_real_2d(features%coarse_0_atomic_coords, 3, skala_state%layout%natom)
         features%coarse_0_atomic_coords(:, :) = skala_state%layout%coarse_0_atomic_coords
      END IF

   END SUBROUTINE copy_cached_layout

! **************************************************************************************************
!> \brief Copy the full (non-routing) cached layout into a feature bundle.
!> \param features ...
!> \param needs_grid_coordinate_array ...
! **************************************************************************************************
   SUBROUTINE copy_full_layout(features, needs_grid_coordinate_array)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      LOGICAL, INTENT(IN)                                :: needs_grid_coordinate_array

      features%grid_weights => skala_state%layout%grid_weights
      features%local_feature_offsets => skala_state%layout%local_feature_offsets
      features%local_feature_rows => skala_state%layout%local_feature_rows
      IF (needs_grid_coordinate_array) THEN
         features%grid_coords => skala_state%layout%grid_coords
         features%atomic_grid_weights => skala_state%layout%atomic_grid_weights
      ELSE
         NULLIFY (features%grid_coords, features%atomic_grid_weights)
      END IF
      features%owns_layout_arrays = .FALSE.

   END SUBROUTINE copy_full_layout

! **************************************************************************************************
!> \brief Copy the per-rank routing layout into a feature bundle.
!> \param features ...
!> \param needs_coordinate_array ...
! **************************************************************************************************
   SUBROUTINE copy_routing_layout(features, needs_coordinate_array)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      LOGICAL, INTENT(IN)                                :: needs_coordinate_array

      INTEGER                                            :: feature_pos

      IF (.NOT. features%owns_layout_arrays) THEN
         NULLIFY (features%grid_coords, features%grid_weights, features%atomic_grid_weights, &
                  features%local_feature_offsets, features%local_feature_rows)
      END IF
      features%owns_layout_arrays = .TRUE.

      IF (needs_coordinate_array) THEN
         CALL ensure_int_1d(features%local_feature_offsets, skala_state%layout%nflat_local + 1)
         CALL ensure_int_1d(features%local_feature_rows, SIZE(skala_state%layout%local_feature_rows))
         CALL ensure_real_2d(features%grid_coords, 3, SIZE(skala_state%layout%local_feature_rows))
         CALL ensure_real_1d(features%grid_weights, SIZE(skala_state%layout%local_feature_rows))
         CALL ensure_real_1d(features%atomic_grid_weights, SIZE(skala_state%layout%local_feature_rows))
         features%local_feature_offsets(:) = skala_state%layout%local_feature_offsets
         DO feature_pos = 1, SIZE(skala_state%layout%local_feature_rows)
            features%grid_coords(:, feature_pos) = &
               skala_state%layout%local_feature_grid_coords(:, feature_pos)
            features%grid_weights(feature_pos) = &
               skala_state%layout%local_feature_grid_weights(feature_pos)
            features%atomic_grid_weights(feature_pos) = &
               skala_state%layout%local_feature_atomic_grid_weights(feature_pos)
            features%local_feature_rows(feature_pos) = feature_pos
         END DO
      END IF

      features%route_point_recv_counts => skala_state%layout%route_point_recv_counts
      features%route_point_recv_displs => skala_state%layout%route_point_recv_displs
      features%route_point_send_counts => skala_state%layout%route_point_send_counts
      features%route_point_send_displs => skala_state%layout%route_point_send_displs
      features%route_chunk_offsets => skala_state%layout%route_chunk_offsets
      features%route_chunk_rows => skala_state%layout%route_chunk_rows
      features%route_row_offsets => skala_state%layout%route_row_offsets
      features%route_row_positions => skala_state%layout%route_row_positions
      features%route_send_local_rows => skala_state%layout%route_send_local_rows
      features%route_feature_send_counts => skala_state%layout%route_feature_send_counts
      features%route_feature_send_displs => skala_state%layout%route_feature_send_displs
      features%route_feature_recv_counts => skala_state%layout%route_feature_recv_counts
      features%route_feature_recv_displs => skala_state%layout%route_feature_recv_displs
      features%route_send_local_features => skala_state%layout%route_send_local_features
      features%route_recv_chunk_rows => skala_state%layout%route_recv_chunk_rows

   END SUBROUTINE copy_routing_layout

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
!> \param local_atoms ...
!> \param atom_chunk_owner ...
!> \param group ...
! **************************************************************************************************
   SUBROUTINE build_atom_chunk_routes(cache, local_to_global, local_atoms, atom_chunk_owner, group)
      TYPE(skala_gpw_layout_cache_type), INTENT(INOUT)   :: cache
      INTEGER, DIMENSION(:), INTENT(IN)                  :: local_to_global, local_atoms, &
                                                            atom_chunk_owner

      CLASS(mp_comm_type), INTENT(IN)                    :: group

      INTEGER                                            :: chunk_row, dest, feature_begin, feature_end, &
                                                            feature_pos, local_feature, local_row, meta_pos, &
                                                            nflat_local, nproc, nrecv, nsend, pe, point_pos, &
                                                            route_index, row, row_position
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: chunk_cursor, chunk_route_counts, &
                                                            feature_cursor, feature_recv_counts, &
                                                            feature_recv_displs, feature_send_counts, &
                                                            feature_send_displs, meta_recv_counts, &
                                                            meta_recv_displs, meta_send_counts, &
                                                            meta_send_displs, recv_meta, route_cursor, &
                                                            route_feature_dest, route_last_row, &
                                                            route_last_position, send_meta

      nflat_local = SIZE(cache%local_feature_offsets) - 1
      nproc = SIZE(cache%route_point_send_counts)
      ALLOCATE (route_feature_dest(SIZE(local_to_global)), &
                feature_send_counts(nproc), feature_send_displs(nproc), &
                feature_recv_counts(nproc), feature_recv_displs(nproc), &
                route_last_row(nproc), route_last_position(nproc))
      cache%route_point_send_counts = 0
      feature_send_counts = 0
      route_feature_dest = 0
      route_last_row = 0
      ! Route each local grid row only once to each atom-chunk owner.
      DO local_row = 1, nflat_local
         feature_begin = cache%local_feature_offsets(local_row)
         feature_end = cache%local_feature_offsets(local_row + 1) - 1
         DO local_feature = feature_begin, feature_end
            dest = atom_chunk_owner(local_atoms(local_feature))
            CPASSERT(dest > 0)
            route_feature_dest(local_feature) = dest
            feature_send_counts(dest) = feature_send_counts(dest) + 1
            IF (route_last_row(dest) /= local_row) THEN
               cache%route_point_send_counts(dest) = cache%route_point_send_counts(dest) + 1
               route_last_row(dest) = local_row
            END IF
         END DO
      END DO
      CALL counts_to_displs(cache%route_point_send_counts, cache%route_point_send_displs)
      CALL counts_to_displs(feature_send_counts, feature_send_displs)
      CALL group%alltoall(cache%route_point_send_counts, cache%route_point_recv_counts, 1)
      CALL counts_to_displs(cache%route_point_recv_counts, cache%route_point_recv_displs)
      CALL group%alltoall(feature_send_counts, feature_recv_counts, 1)
      CALL counts_to_displs(feature_recv_counts, feature_recv_displs)

      nsend = SUM(cache%route_point_send_counts)
      nrecv = SUM(cache%route_point_recv_counts)
      ALLOCATE (cache%route_feature_send_counts(nproc), &
                cache%route_feature_send_displs(nproc), &
                cache%route_feature_recv_counts(nproc), &
                cache%route_feature_recv_displs(nproc), &
                cache%route_send_local_features(SIZE(local_to_global)), &
                cache%route_recv_chunk_rows(cache%chunk_feature_count))
      cache%route_feature_send_counts(:) = feature_send_counts
      cache%route_feature_send_displs(:) = feature_send_displs
      cache%route_feature_recv_counts(:) = feature_recv_counts
      cache%route_feature_recv_displs(:) = feature_recv_displs
      cache%route_send_local_features = 0
      cache%route_recv_chunk_rows = 0
      ALLOCATE (cache%route_send_local_rows(nsend), &
                cache%route_row_offsets(nflat_local + 1), &
                cache%route_row_positions(nsend), &
                route_cursor(nproc), feature_cursor(nproc), &
                send_meta(2*SIZE(local_to_global)), &
                recv_meta(2*cache%chunk_feature_count), &
                meta_send_counts(nproc), meta_send_displs(nproc), &
                meta_recv_counts(nproc), meta_recv_displs(nproc))
      cache%route_send_local_rows = 0
      cache%route_row_positions = 0
      send_meta = 0
      route_cursor(:) = cache%route_point_send_displs + 1
      feature_cursor(:) = feature_send_displs + 1
      route_last_row = 0
      route_last_position = 0
      cache%route_row_offsets(1) = 1
      row_position = 1
      DO local_row = 1, nflat_local
         feature_begin = cache%local_feature_offsets(local_row)
         feature_end = cache%local_feature_offsets(local_row + 1) - 1
         DO local_feature = feature_begin, feature_end
            dest = route_feature_dest(local_feature)
            IF (route_last_row(dest) /= local_row) THEN
               point_pos = route_cursor(dest)
               route_cursor(dest) = route_cursor(dest) + 1
               route_last_row(dest) = local_row
               route_last_position(dest) = point_pos
               cache%route_send_local_rows(point_pos) = local_row
               cache%route_row_positions(row_position) = point_pos
               row_position = row_position + 1
            ELSE
               point_pos = route_last_position(dest)
            END IF
            feature_pos = feature_cursor(dest)
            feature_cursor(dest) = feature_cursor(dest) + 1
            cache%route_send_local_features(feature_pos) = local_feature
            send_meta(2*feature_pos - 1) = local_to_global(local_feature)
            send_meta(2*feature_pos) = point_pos - cache%route_point_send_displs(dest)
         END DO
         cache%route_row_offsets(local_row + 1) = row_position
      END DO
      CPASSERT(row_position == nsend + 1)

      ! Tell each chunk owner which feature rows share a compressed route point.
      meta_send_counts(:) = 2*feature_send_counts
      meta_send_displs(:) = 2*feature_send_displs
      meta_recv_counts(:) = 2*feature_recv_counts
      meta_recv_displs(:) = 2*feature_recv_displs
      CALL group%alltoall(send_meta, meta_send_counts, meta_send_displs, recv_meta, &
                          meta_recv_counts, meta_recv_displs)

      ! Store the route-to-feature expansion in CSR form for gradient reduction.
      ALLOCATE (cache%chunk_return_positions(cache%chunk_feature_count), &
                cache%route_chunk_offsets(nrecv + 1), &
                cache%route_chunk_rows(cache%chunk_feature_count), &
                chunk_route_counts(nrecv), chunk_cursor(nrecv))
      cache%chunk_return_positions = 0
      cache%route_chunk_rows = 0
      chunk_route_counts = 0
      DO pe = 1, nproc
         DO feature_pos = feature_recv_displs(pe) + 1, &
            feature_recv_displs(pe) + feature_recv_counts(pe)
            meta_pos = 2*feature_pos
            row = recv_meta(meta_pos - 1)
            route_index = recv_meta(meta_pos)
            point_pos = cache%route_point_recv_displs(pe) + route_index
            CPASSERT(point_pos >= 1 .AND. point_pos <= nrecv)
            chunk_row = row - cache%chunk_feature_begin + 1
            CPASSERT(chunk_row >= 1 .AND. chunk_row <= cache%chunk_feature_count)
            cache%chunk_return_positions(chunk_row) = point_pos
            cache%route_recv_chunk_rows(feature_pos) = chunk_row
            chunk_route_counts(point_pos) = chunk_route_counts(point_pos) + 1
         END DO
      END DO
      cache%route_chunk_offsets(1) = 1
      DO point_pos = 1, nrecv
         cache%route_chunk_offsets(point_pos + 1) = &
            cache%route_chunk_offsets(point_pos) + chunk_route_counts(point_pos)
      END DO
      chunk_cursor(:) = cache%route_chunk_offsets(1:nrecv)
      DO pe = 1, nproc
         DO feature_pos = feature_recv_displs(pe) + 1, &
            feature_recv_displs(pe) + feature_recv_counts(pe)
            meta_pos = 2*feature_pos
            row = recv_meta(meta_pos - 1)
            chunk_row = row - cache%chunk_feature_begin + 1
            point_pos = cache%chunk_return_positions(chunk_row)
            cache%route_chunk_rows(chunk_cursor(point_pos)) = chunk_row
            chunk_cursor(point_pos) = chunk_cursor(point_pos) + 1
         END DO
      END DO

      CPASSERT(SUM(feature_send_counts) == SIZE(local_to_global))
      CPASSERT(SUM(feature_recv_counts) == cache%chunk_feature_count)
      CPASSERT(SUM(cache%route_point_send_counts) == nsend)
      CPASSERT(SUM(cache%route_point_recv_counts) == nrecv)
      CPASSERT(ALL(cache%chunk_return_positions > 0))
      CPASSERT(ALL(cache%route_chunk_rows > 0))
      CPASSERT(ALL(cache%route_row_positions > 0))
      CPASSERT(ALL(cache%route_send_local_rows > 0))

      DEALLOCATE (chunk_cursor, chunk_route_counts, feature_cursor, feature_recv_counts, &
                  feature_recv_displs, feature_send_counts, feature_send_displs, &
                  meta_recv_counts, meta_recv_displs, meta_send_counts, meta_send_displs, &
                  recv_meta, route_cursor, route_feature_dest, route_last_position, route_last_row, &
                  send_meta)

   END SUBROUTINE build_atom_chunk_routes

! **************************************************************************************************
!> \brief Materialize the current rank's atom chunk static layout.
!> \param cache ...
! **************************************************************************************************
   SUBROUTINE build_atom_chunk_layout(cache)
      TYPE(skala_gpw_layout_cache_type), INTENT(INOUT)   :: cache

      INTEGER                                            :: max_grid_size

      IF (cache%chunk_feature_count <= 0 .OR. cache%chunk_natom <= 0) RETURN

      ALLOCATE (cache%chunk_atomic_grid_sizes(cache%chunk_natom), &
                cache%chunk_coarse_0_atomic_coords(3, cache%chunk_natom))
      cache%chunk_atomic_grid_sizes(:) = &
         cache%atomic_grid_sizes(cache%chunk_atom_begin:cache%chunk_atom_end)
      cache%chunk_coarse_0_atomic_coords(:, :) = &
         cache%coarse_0_atomic_coords(:, cache%chunk_atom_begin:cache%chunk_atom_end)

      max_grid_size = MAXVAL(INT(cache%chunk_atomic_grid_sizes))
      CALL allocate_empty_bound_shape(cache%chunk_atomic_grid_size_bound_shape, max_grid_size)

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

      INTEGER                                            :: chunk_row, dyn_base, local_row, &
                                                            ndynamic_route_per_point, nrecv, nsend, &
                                                            point_pos, src_base
      INTEGER, DIMENSION(:), POINTER                     :: recv_counts, recv_displs, send_counts, &
                                                            send_displs
      REAL(KIND=dp), DIMENSION(:), POINTER               :: recv_dynamic, send_dynamic

      nsend = SIZE(skala_state%layout%route_send_local_rows)
      nrecv = SUM(skala_state%layout%route_point_recv_counts)
      CPASSERT(nsend == SUM(skala_state%layout%route_point_send_counts))
      ndynamic_route_per_point = ndynamic_per_point
      IF (collapse_spin_dynamics) ndynamic_route_per_point = nrks_dynamic_per_point

      CALL ensure_real_1d(skala_state%scratch%route_send_dynamic, &
                          ndynamic_route_per_point*nsend)
      CALL ensure_real_1d(skala_state%scratch%route_recv_dynamic, &
                          ndynamic_route_per_point*nrecv)
      CALL ensure_int_1d(skala_state%scratch%route_send_counts, skala_state%layout%nproc)
      CALL ensure_int_1d(skala_state%scratch%route_send_displs, skala_state%layout%nproc)
      CALL ensure_int_1d(skala_state%scratch%route_recv_counts, skala_state%layout%nproc)
      CALL ensure_int_1d(skala_state%scratch%route_recv_displs, skala_state%layout%nproc)
      send_dynamic => skala_state%scratch%route_send_dynamic
      recv_dynamic => skala_state%scratch%route_recv_dynamic
      send_counts => skala_state%scratch%route_send_counts
      send_displs => skala_state%scratch%route_send_displs
      recv_counts => skala_state%scratch%route_recv_counts
      recv_displs => skala_state%scratch%route_recv_displs
      send_counts(:) = ndynamic_route_per_point*skala_state%layout%route_point_send_counts
      send_displs(:) = ndynamic_route_per_point*skala_state%layout%route_point_send_displs
      recv_counts(:) = ndynamic_route_per_point*skala_state%layout%route_point_recv_counts
      recv_displs(:) = ndynamic_route_per_point*skala_state%layout%route_point_recv_displs
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(skala_state, local_dynamic, ndynamic_route_per_point, nsend, send_dynamic) &
!$OMP PRIVATE(dyn_base, local_row, point_pos, src_base)
      DO point_pos = 1, nsend
         dyn_base = ndynamic_route_per_point*(point_pos - 1)
         local_row = skala_state%layout%route_send_local_rows(point_pos)
         src_base = ndynamic_route_per_point*(local_row - 1)
         send_dynamic(dyn_base + 1:dyn_base + ndynamic_route_per_point) = &
            local_dynamic(src_base + 1:src_base + ndynamic_route_per_point)
      END DO
!$OMP END PARALLEL DO

      CALL group%alltoall(send_dynamic, send_counts, send_displs, recv_dynamic, recv_counts, &
                          recv_displs)

      features%uses_collapsed_rks_dynamic = collapse_spin_dynamics
      IF (skala_state%layout%chunk_feature_count > 0) THEN
         IF (collapse_spin_dynamics) THEN
            CALL ensure_real_2d(features%chunk_density, skala_state%layout%chunk_feature_count, 1)
            CALL ensure_real_3d(features%chunk_grad, skala_state%layout%chunk_feature_count, 3, 1)
            CALL ensure_real_2d(features%chunk_kin, skala_state%layout%chunk_feature_count, 1)
         ELSE
            CALL ensure_real_2d(features%chunk_density, skala_state%layout%chunk_feature_count, 2)
            CALL ensure_real_3d(features%chunk_grad, skala_state%layout%chunk_feature_count, 3, 2)
            CALL ensure_real_2d(features%chunk_kin, skala_state%layout%chunk_feature_count, 2)
         END IF

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(skala_state, collapse_spin_dynamics, features, ndynamic_route_per_point, nrecv, recv_dynamic) &
!$OMP PRIVATE(chunk_row, dyn_base, point_pos)
         DO chunk_row = 1, skala_state%layout%chunk_feature_count
            point_pos = skala_state%layout%chunk_return_positions(chunk_row)
            CPASSERT(point_pos >= 1 .AND. point_pos <= nrecv)
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
!$OMP END PARALLEL DO
      END IF

      NULLIFY (recv_counts, recv_displs, recv_dynamic, send_counts, send_displs, send_dynamic)

   END SUBROUTINE route_atom_chunk_dynamics

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
!> \brief Release the module-level layout cache.
! **************************************************************************************************
   SUBROUTINE skala_gpw_features_release_cache()
      CALL release_layout_cache(skala_state%layout)
      IF (ALLOCATED(skala_state%scratch%local_dynamic)) THEN
         DEALLOCATE (skala_state%scratch%local_dynamic)
      END IF
      IF (ALLOCATED(skala_state%scratch%global_dynamic)) THEN
         DEALLOCATE (skala_state%scratch%global_dynamic)
      END IF
      IF (ALLOCATED(skala_state%scratch%route_send_dynamic)) THEN
         DEALLOCATE (skala_state%scratch%route_send_dynamic)
      END IF
      IF (ALLOCATED(skala_state%scratch%route_recv_dynamic)) THEN
         DEALLOCATE (skala_state%scratch%route_recv_dynamic)
      END IF
      IF (ALLOCATED(skala_state%scratch%route_send_counts)) THEN
         DEALLOCATE (skala_state%scratch%route_send_counts)
      END IF
      IF (ALLOCATED(skala_state%scratch%route_send_displs)) THEN
         DEALLOCATE (skala_state%scratch%route_send_displs)
      END IF
      IF (ALLOCATED(skala_state%scratch%route_recv_counts)) THEN
         DEALLOCATE (skala_state%scratch%route_recv_counts)
      END IF
      IF (ALLOCATED(skala_state%scratch%route_recv_displs)) THEN
         DEALLOCATE (skala_state%scratch%route_recv_displs)
      END IF
      skala_state%layout_mode_known = .FALSE.
      skala_state%layout_mode_full = .FALSE.
   END SUBROUTINE skala_gpw_features_release_cache

! **************************************************************************************************
!> \brief Release cached layout arrays.
!> \param cache ...
! **************************************************************************************************
   SUBROUTINE release_layout_cache(cache)
      TYPE(skala_gpw_layout_cache_type), INTENT(INOUT)   :: cache

      IF (cache%grad_inputs_active) THEN
         CALL torch_dict_release(cache%grad_inputs)
         cache%grad_inputs_active = .FALSE.
      END IF

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
         CALL torch_dict_release(cache%chunk_static_inputs)
         cache%chunk_static_tensors_active = .FALSE.
      END IF

      IF (ALLOCATED(cache%chunk_return_positions)) DEALLOCATE (cache%chunk_return_positions)
      IF (ALLOCATED(cache%route_chunk_offsets)) DEALLOCATE (cache%route_chunk_offsets)
      IF (ALLOCATED(cache%route_chunk_rows)) DEALLOCATE (cache%route_chunk_rows)
      IF (ALLOCATED(cache%route_point_recv_counts)) DEALLOCATE (cache%route_point_recv_counts)
      IF (ALLOCATED(cache%route_point_recv_displs)) DEALLOCATE (cache%route_point_recv_displs)
      IF (ALLOCATED(cache%route_point_send_counts)) DEALLOCATE (cache%route_point_send_counts)
      IF (ALLOCATED(cache%route_point_send_displs)) DEALLOCATE (cache%route_point_send_displs)
      IF (ALLOCATED(cache%route_row_offsets)) DEALLOCATE (cache%route_row_offsets)
      IF (ALLOCATED(cache%route_row_positions)) DEALLOCATE (cache%route_row_positions)
      IF (ALLOCATED(cache%route_send_local_rows)) DEALLOCATE (cache%route_send_local_rows)
      IF (ALLOCATED(cache%route_feature_send_counts)) DEALLOCATE (cache%route_feature_send_counts)
      IF (ALLOCATED(cache%route_feature_send_displs)) DEALLOCATE (cache%route_feature_send_displs)
      IF (ALLOCATED(cache%route_feature_recv_counts)) DEALLOCATE (cache%route_feature_recv_counts)
      IF (ALLOCATED(cache%route_feature_recv_displs)) DEALLOCATE (cache%route_feature_recv_displs)
      IF (ALLOCATED(cache%route_send_local_features)) DEALLOCATE (cache%route_send_local_features)
      IF (ALLOCATED(cache%route_recv_chunk_rows)) DEALLOCATE (cache%route_recv_chunk_rows)
      IF (ALLOCATED(cache%dynamic_counts)) DEALLOCATE (cache%dynamic_counts)
      IF (ALLOCATED(cache%dynamic_displs)) DEALLOCATE (cache%dynamic_displs)
      IF (ALLOCATED(cache%feature_source_points)) DEALLOCATE (cache%feature_source_points)
      IF (ALLOCATED(cache%atomic_grid_sizes)) DEALLOCATE (cache%atomic_grid_sizes)
      IF (ALLOCATED(cache%chunk_atomic_grid_sizes)) DEALLOCATE (cache%chunk_atomic_grid_sizes)
      IF (ALLOCATED(cache%local_feature_offsets)) DEALLOCATE (cache%local_feature_offsets)
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
      IF (ALLOCATED(cache%local_feature_grid_coords)) DEALLOCATE (cache%local_feature_grid_coords)
      IF (ALLOCATED(cache%local_feature_grid_weights)) DEALLOCATE (cache%local_feature_grid_weights)
      IF (ALLOCATED(cache%local_feature_atomic_grid_weights)) THEN
         DEALLOCATE (cache%local_feature_atomic_grid_weights)
      END IF

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
!> \param keep_host_arrays keep the reusable host arrays for the next SCF step
! **************************************************************************************************
   SUBROUTINE skala_gpw_feature_release(features, keep_host_arrays)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_host_arrays

      LOGICAL                                            :: keep

      keep = .FALSE.
      IF (PRESENT(keep_host_arrays)) keep = keep_host_arrays

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
         features%inputs_active = .FALSE.
         features%active = .FALSE.
         features%owns_coordinate_tensor = .FALSE.
         features%owns_grid_coordinate_tensor = .FALSE.
         features%owns_weight_tensors = .FALSE.
         features%owns_dynamic_tensors = .TRUE.
         features%owns_inputs = .TRUE.
         features%owns_static_tensors = .TRUE.
         features%uses_atom_chunks = .FALSE.
         features%uses_collapsed_rks_dynamic = .FALSE.
      END IF

      IF (features%owns_layout_arrays) THEN
         IF (.NOT. keep) THEN
            IF (ASSOCIATED(features%local_feature_offsets)) DEALLOCATE (features%local_feature_offsets)
            IF (ASSOCIATED(features%local_feature_rows)) DEALLOCATE (features%local_feature_rows)
            IF (ASSOCIATED(features%grid_coords)) DEALLOCATE (features%grid_coords)
            IF (ASSOCIATED(features%grid_weights)) DEALLOCATE (features%grid_weights)
            IF (ASSOCIATED(features%atomic_grid_weights)) DEALLOCATE (features%atomic_grid_weights)
            NULLIFY (features%local_feature_offsets, features%local_feature_rows, &
                     features%grid_coords, features%grid_weights, features%atomic_grid_weights)
            features%owns_layout_arrays = .FALSE.
         END IF
      ELSE
         NULLIFY (features%local_feature_offsets, features%local_feature_rows, &
                  features%grid_coords, features%grid_weights, features%atomic_grid_weights)
      END IF

      IF (.NOT. keep) THEN
         IF (ALLOCATED(features%chunk_density)) DEALLOCATE (features%chunk_density)
         IF (ALLOCATED(features%chunk_grad)) DEALLOCATE (features%chunk_grad)
         IF (ALLOCATED(features%chunk_kin)) DEALLOCATE (features%chunk_kin)
         IF (ALLOCATED(features%density)) DEALLOCATE (features%density)
         IF (ALLOCATED(features%grad)) DEALLOCATE (features%grad)
         IF (ALLOCATED(features%kin)) DEALLOCATE (features%kin)
      END IF
      NULLIFY (features%route_point_recv_counts, features%route_point_recv_displs, &
               features%route_point_send_counts, features%route_point_send_displs, &
               features%route_chunk_offsets, features%route_chunk_rows, &
               features%route_row_offsets, features%route_row_positions, &
               features%route_send_local_rows, features%route_feature_send_counts, &
               features%route_feature_send_displs, features%route_feature_recv_counts, &
               features%route_feature_recv_displs, features%route_send_local_features, &
               features%route_recv_chunk_rows)
      IF (.NOT. keep) THEN
         IF (ALLOCATED(features%atomic_grid_sizes)) DEALLOCATE (features%atomic_grid_sizes)
         IF (ALLOCATED(features%coarse_0_atomic_coords)) DEALLOCATE (features%coarse_0_atomic_coords)
         IF (ALLOCATED(features%atomic_grid_size_bound_shape)) THEN
            DEALLOCATE (features%atomic_grid_size_bound_shape)
         END IF
      END IF
      features%chunk_feature_count = 0
      features%nflat = 0
      features%nflat_local = 0
      features%atom_partition = skala_gpw_atom_partition_hard
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
      IF (.NOT. skala_state%layout%active) RETURN
      IF (skala_state%layout%chunk_natom <= 0) RETURN
      IF (max_rows <= 0) THEN
         nsubchunks = 1
         RETURN
      END IF

      rows = 0
      DO iatom = 1, skala_state%layout%chunk_natom
         atom_rows = INT(skala_state%layout%chunk_atomic_grid_sizes(iatom))
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
!> \brief Sort chunk layouts in place, largest padded block first.
!> \param padded_rows key used for the descending order
!> \param atom_begin ...
!> \param atom_count ...
!> \param row_begin ...
!> \param row_count ...
! **************************************************************************************************
   SUBROUTINE sort_chunks_by_padded_rows(padded_rows, atom_begin, atom_count, row_begin, row_count)
      INTEGER(KIND=int_8), DIMENSION(:), INTENT(IN)      :: padded_rows
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: atom_begin, atom_count, row_begin, &
                                                            row_count

      INTEGER                                            :: atom_begin_tmp, atom_count_tmp, chunk, &
                                                            insert_at, row_begin_tmp, row_count_tmp
      INTEGER(KIND=int_8)                                :: key
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: keys

      CPASSERT(SIZE(atom_begin) == SIZE(padded_rows))
      CPASSERT(SIZE(atom_count) == SIZE(padded_rows))
      CPASSERT(SIZE(row_begin) == SIZE(padded_rows))
      CPASSERT(SIZE(row_count) == SIZE(padded_rows))

      ALLOCATE (keys(SIZE(padded_rows)))
      keys(:) = padded_rows
      DO chunk = 2, SIZE(keys)
         key = keys(chunk)
         atom_begin_tmp = atom_begin(chunk)
         atom_count_tmp = atom_count(chunk)
         row_begin_tmp = row_begin(chunk)
         row_count_tmp = row_count(chunk)
         insert_at = chunk
         DO WHILE (insert_at > 1 .AND. keys(insert_at - 1) < key)
            keys(insert_at) = keys(insert_at - 1)
            atom_begin(insert_at) = atom_begin(insert_at - 1)
            atom_count(insert_at) = atom_count(insert_at - 1)
            row_begin(insert_at) = row_begin(insert_at - 1)
            row_count(insert_at) = row_count(insert_at - 1)
            insert_at = insert_at - 1
         END DO
         keys(insert_at) = key
         atom_begin(insert_at) = atom_begin_tmp
         atom_count(insert_at) = atom_count_tmp
         row_begin(insert_at) = row_begin_tmp
         row_count(insert_at) = row_count_tmp
      END DO
      DEALLOCATE (keys)

   END SUBROUTINE sort_chunks_by_padded_rows

! **************************************************************************************************
!> \brief Build atom and row bounds for all atom-contiguous subchunks in one pass.
!> \param max_rows ...
!> \param atom_begin ...
!> \param atom_count ...
!> \param row_begin ...
!> \param row_count ...
! **************************************************************************************************
   SUBROUTINE skala_gpw_atom_subchunk_layout(max_rows, atom_begin, atom_count, row_begin, &
                                             row_count)
      INTEGER, INTENT(IN)                                :: max_rows
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: atom_begin, atom_count, row_begin, &
                                                            row_count

      INTEGER                                            :: atom_rows, iatom, nsubchunks, rows, &
                                                            subchunk
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: padded_rows

      nsubchunks = skala_gpw_atom_subchunk_count(max_rows)
      ALLOCATE (atom_begin(nsubchunks), atom_count(nsubchunks), row_begin(nsubchunks), &
                row_count(nsubchunks))
      IF (nsubchunks == 0) RETURN

      CPASSERT(skala_state%layout%active)
      CPASSERT(skala_state%layout%chunk_natom > 0)

      atom_begin = 0
      atom_count = 0
      row_begin = 0
      row_count = 0

      IF (max_rows <= 0) THEN
         atom_begin(1) = 1
         atom_count(1) = skala_state%layout%chunk_natom
         row_begin(1) = 1
         row_count(1) = skala_state%layout%chunk_feature_count
         RETURN
      END IF

      subchunk = 1
      atom_begin(subchunk) = 1
      row_begin(subchunk) = 1
      rows = 0
      DO iatom = 1, skala_state%layout%chunk_natom
         atom_rows = INT(skala_state%layout%chunk_atomic_grid_sizes(iatom))
         IF (rows > 0 .AND. rows + atom_rows > max_rows) THEN
            atom_count(subchunk) = iatom - atom_begin(subchunk)
            row_count(subchunk) = rows
            subchunk = subchunk + 1
            atom_begin(subchunk) = iatom
            row_begin(subchunk) = row_begin(subchunk - 1) + row_count(subchunk - 1)
            rows = 0
         END IF
         rows = rows + atom_rows
      END DO
      atom_count(subchunk) = skala_state%layout%chunk_natom - atom_begin(subchunk) + 1
      row_count(subchunk) = rows

      CPASSERT(subchunk == nsubchunks)

      ALLOCATE (padded_rows(nsubchunks))
      DO subchunk = 1, nsubchunks
         iatom = atom_begin(subchunk) + atom_count(subchunk) - 1
         padded_rows(subchunk) = INT(atom_count(subchunk), KIND=int_8)* &
                                 MAXVAL(skala_state%layout%chunk_atomic_grid_sizes( &
                                        atom_begin(subchunk):iatom))
      END DO

      ! Evaluate the largest padded block first so Torch can reuse its CUDA allocations.
      CALL sort_chunks_by_padded_rows(padded_rows, atom_begin, atom_count, row_begin, row_count)
      DEALLOCATE (padded_rows)

   END SUBROUTINE skala_gpw_atom_subchunk_layout

! **************************************************************************************************
!> \brief Release the view tensors/dict of a reusable subchunk bundle.
!>        Host arrays and the parent-owned tensors are left untouched.
!> \param features ...
! **************************************************************************************************
   SUBROUTINE release_subchunk_tensors(features)
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features

      IF (.NOT. features%active) RETURN
      CALL torch_tensor_release(features%grid_coords_t)
      CALL torch_tensor_release(features%grid_weights_t)
      CALL torch_tensor_release(features%atomic_grid_weights_t)
      CALL torch_tensor_release(features%atomic_grid_sizes_t)
      CALL torch_tensor_release(features%coarse_0_atomic_coords_t)
      CALL torch_tensor_release(features%atomic_grid_size_bound_shape_t)
      CALL torch_tensor_release(features%density_t)
      CALL torch_tensor_release(features%grad_t)
      CALL torch_tensor_release(features%kin_t)
      IF (features%uses_collapsed_rks_dynamic) THEN
         CALL torch_tensor_release(features%density_input_t)
         CALL torch_tensor_release(features%grad_input_t)
         CALL torch_tensor_release(features%kin_input_t)
      END IF
      IF (features%inputs_active) THEN
         CALL torch_dict_release(features%inputs)
         features%inputs_active = .FALSE.
      END IF
      features%active = .FALSE.
      features%owns_coordinate_tensor = .FALSE.
      features%owns_grid_coordinate_tensor = .FALSE.
      features%owns_weight_tensors = .FALSE.
      features%owns_dynamic_tensors = .FALSE.
      features%owns_inputs = .TRUE.
      features%owns_static_tensors = .FALSE.

   END SUBROUTINE release_subchunk_tensors

! **************************************************************************************************
!> \brief Build an atom-contiguous subchunk feature bundle from precomputed bounds.
!> \param parent ...
!> \param features ...
!> \param atom_begin ...
!> \param atom_count ...
!> \param row_begin ...
!> \param row_count ...
! **************************************************************************************************
   SUBROUTINE skala_gpw_feature_build_atom_subchunk_bounds(parent, features, atom_begin, &
                                                           atom_count, row_begin, row_count)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: parent
      TYPE(skala_gpw_feature_type), INTENT(INOUT)        :: features
      INTEGER, INTENT(IN)                                :: atom_begin, atom_count, row_begin, &
                                                            row_count

      INTEGER                                            :: atom_end, max_grid_size

      IF (features%active) CALL release_subchunk_tensors(features)
      CPASSERT(parent%uses_atom_chunks)
      atom_end = atom_begin + atom_count - 1
      CPASSERT(atom_count > 0)
      CPASSERT(row_count > 0)
      max_grid_size = MAXVAL(INT(skala_state%layout%chunk_atomic_grid_sizes(atom_begin:atom_end)))

      IF (ALLOCATED(features%atomic_grid_size_bound_shape)) THEN
         IF (SIZE(features%atomic_grid_size_bound_shape, 2) /= max_grid_size) THEN
            DEALLOCATE (features%atomic_grid_size_bound_shape)
         END IF
      END IF
      IF (.NOT. ALLOCATED(features%atomic_grid_size_bound_shape)) THEN
         ALLOCATE (features%atomic_grid_size_bound_shape(0, max_grid_size))
      END IF

      CALL add_subchunk_feature_tensors(parent, features, atom_begin, atom_count, row_begin, &
                                        row_count)
      features%active = .TRUE.

   END SUBROUTINE skala_gpw_feature_build_atom_subchunk_bounds

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

      CPASSERT(skala_state%layout%chunk_static_tensors_active)
      CPASSERT(parent%active)
      CPASSERT(ALLOCATED(features%atomic_grid_size_bound_shape))

      features%owns_coordinate_tensor = .FALSE.
      features%owns_dynamic_tensors = .TRUE.
      features%owns_inputs = .TRUE.
      features%owns_static_tensors = .TRUE.
      features%uses_collapsed_rks_dynamic = parent%uses_collapsed_rks_dynamic

      CALL torch_tensor_narrow(skala_state%layout%chunk_grid_coords_t, 0, row_begin - 1, &
                               row_count, features%grid_coords_t)
      CALL torch_tensor_narrow(skala_state%layout%chunk_grid_weights_t, 0, row_begin - 1, &
                               row_count, features%grid_weights_t)
      CALL torch_tensor_narrow(skala_state%layout%chunk_atomic_grid_weights_t, 0, row_begin - 1, &
                               row_count, features%atomic_grid_weights_t)
      CALL torch_tensor_narrow(skala_state%layout%chunk_atomic_grid_sizes_t, 0, atom_begin - 1, &
                               atom_count, features%atomic_grid_sizes_t)
      CALL torch_tensor_narrow(skala_state%layout%chunk_coarse_0_atomic_coords_t, 0, &
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

      IF (.NOT. features%inputs_active) THEN
         CALL torch_dict_create(features%inputs)
         features%inputs_active = .TRUE.
      END IF
      CALL torch_dict_set(features%inputs, skala_key_grid_coords, features%grid_coords_t)
      CALL torch_dict_set(features%inputs, skala_key_grid_weights, features%grid_weights_t)
      CALL torch_dict_set(features%inputs, skala_key_atomic_grid_weights, &
                          features%atomic_grid_weights_t)
      CALL torch_dict_set(features%inputs, skala_key_atomic_grid_sizes, &
                          features%atomic_grid_sizes_t)
      CALL torch_dict_set(features%inputs, skala_key_atomic_grid_size_bound_shape, &
                          features%atomic_grid_size_bound_shape_t)
      IF (features%uses_collapsed_rks_dynamic) THEN
         CALL torch_dict_set(features%inputs, skala_key_density, features%density_input_t)
         CALL torch_dict_set(features%inputs, skala_key_grad, features%grad_input_t)
         CALL torch_dict_set(features%inputs, skala_key_kin, features%kin_input_t)
      ELSE
         CALL torch_dict_set(features%inputs, skala_key_density, features%density_t)
         CALL torch_dict_set(features%inputs, skala_key_grad, features%grad_t)
         CALL torch_dict_set(features%inputs, skala_key_kin, features%kin_t)
      END IF
      CALL torch_dict_set(features%inputs, skala_key_coarse_0_atomic_coords, &
                          features%coarse_0_atomic_coords_t)

   END SUBROUTINE add_subchunk_feature_tensors

! **************************************************************************************************
!> \brief Stage-upload the three dynamic feature tensors and report reallocation.
!> \param density_t ...
!> \param density ...
!> \param grad_t ...
!> \param grad ...
!> \param kin_t ...
!> \param kin ...
!> \param requires_grad ...
!> \param cuda_device ...
!> \param dynamic_reallocated ...
! **************************************************************************************************
   SUBROUTINE staged_upload_dynamic(density_t, density, grad_t, grad, kin_t, kin, &
                                    requires_grad, cuda_device, dynamic_reallocated)
      TYPE(torch_tensor_type), INTENT(INOUT)             :: density_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: density
      TYPE(torch_tensor_type), INTENT(INOUT)             :: grad_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: grad
      TYPE(torch_tensor_type), INTENT(INOUT)             :: kin_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: kin
      LOGICAL, INTENT(IN)                                :: requires_grad
      INTEGER, INTENT(IN)                                :: cuda_device
      LOGICAL, INTENT(OUT)                               :: dynamic_reallocated

      INTEGER(KIND=int_8)                                :: upload_offset, upload_total
      LOGICAL                                            :: reallocated_density, reallocated_grad, &
                                                            reallocated_kin

      upload_total = INT(SIZE(density), KIND=int_8) + INT(SIZE(grad), KIND=int_8) + &
                     INT(SIZE(kin), KIND=int_8)
      IF (cuda_device >= 0) CALL torch_pinned_reserve(upload_total)
      upload_offset = 0_int_8
      upload_offset = staged_upload(density_t, density, upload_offset, requires_grad, &
                                    reallocated_density)
      upload_offset = staged_upload(grad_t, grad, upload_offset, requires_grad, &
                                    reallocated_grad)
      upload_offset = staged_upload(kin_t, kin, upload_offset, requires_grad, &
                                    reallocated_kin)
      CALL torch_stream_synchronize()
      CALL torch_pinned_finish()
      dynamic_reallocated = reallocated_density .OR. reallocated_grad .OR. reallocated_kin

   END SUBROUTINE staged_upload_dynamic

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

      INTEGER                                            :: grad_mode
      LOGICAL                                            :: dynamic_reallocated, &
                                                            my_requires_weight_grad

      my_requires_weight_grad = .FALSE.
      IF (PRESENT(requires_weight_grad)) my_requires_weight_grad = requires_weight_grad

      features%owns_static_tensors = .FALSE.
      features%owns_coordinate_tensor = .FALSE.
      features%owns_grid_coordinate_tensor = .FALSE.
      features%owns_weight_tensors = .FALSE.
      features%owns_dynamic_tensors = .FALSE.
      features%owns_inputs = .TRUE.
      IF (use_atom_chunks) THEN
         IF (requires_coordinate_grad .OR. requires_stress_grad) THEN
            grad_mode = 2
         ELSE IF (my_requires_weight_grad) THEN
            grad_mode = 1
         ELSE
            grad_mode = 0
         END IF
         CALL ensure_chunk_static_tensor_grad(grad_mode)
         CPASSERT(skala_state%layout%chunk_static_tensors_active)
         features%grid_coords_t = skala_state%layout%chunk_grid_coords_t
         features%grid_weights_t = skala_state%layout%chunk_grid_weights_t
         features%atomic_grid_weights_t = skala_state%layout%chunk_atomic_grid_weights_t
         features%atomic_grid_sizes_t = skala_state%layout%chunk_atomic_grid_sizes_t
         features%atomic_grid_size_bound_shape_t = &
            skala_state%layout%chunk_atomic_grid_size_bound_shape_t

         IF (skala_state%layout%chunk_inputs_active .AND. &
             (skala_state%layout%chunk_inputs_use_collapsed_rks .NEQV. &
              features%uses_collapsed_rks_dynamic)) THEN
            CALL torch_dict_release(skala_state%layout%chunk_inputs)
            skala_state%layout%chunk_inputs_active = .FALSE.
         END IF
         IF (.NOT. features%uses_collapsed_rks_dynamic .AND. &
             skala_state%layout%chunk_dynamic_input_views_active) THEN
            CALL torch_tensor_release(skala_state%layout%chunk_density_input_t)
            CALL torch_tensor_release(skala_state%layout%chunk_grad_input_t)
            CALL torch_tensor_release(skala_state%layout%chunk_kin_input_t)
            skala_state%layout%chunk_dynamic_input_views_active = .FALSE.
         END IF

         CALL staged_upload_dynamic(skala_state%layout%chunk_density_t, features%chunk_density, &
                                    skala_state%layout%chunk_grad_t, features%chunk_grad, &
                                    skala_state%layout%chunk_kin_t, features%chunk_kin, &
                                    requires_grad, skala_state%layout%cuda_device, &
                                    dynamic_reallocated)
         features%density_t = skala_state%layout%chunk_density_t
         features%grad_t = skala_state%layout%chunk_grad_t
         features%kin_t = skala_state%layout%chunk_kin_t
         IF (dynamic_reallocated) THEN
            IF (skala_state%layout%chunk_inputs_active) THEN
               CALL torch_dict_release(skala_state%layout%chunk_inputs)
               skala_state%layout%chunk_inputs_active = .FALSE.
            END IF
            IF (skala_state%layout%chunk_dynamic_input_views_active) THEN
               CALL torch_tensor_release(skala_state%layout%chunk_density_input_t)
               CALL torch_tensor_release(skala_state%layout%chunk_grad_input_t)
               CALL torch_tensor_release(skala_state%layout%chunk_kin_input_t)
               skala_state%layout%chunk_dynamic_input_views_active = .FALSE.
            END IF
         END IF
         skala_state%layout%chunk_dynamic_tensors_active = .TRUE.

         IF (features%uses_collapsed_rks_dynamic .AND. &
             .NOT. skala_state%layout%chunk_dynamic_input_views_active) THEN
            CALL torch_tensor_expand_dim(skala_state%layout%chunk_density_t, 0, 2, &
                                         skala_state%layout%chunk_density_input_t)
            CALL torch_tensor_expand_dim(skala_state%layout%chunk_grad_t, 0, 2, &
                                         skala_state%layout%chunk_grad_input_t)
            CALL torch_tensor_expand_dim(skala_state%layout%chunk_kin_t, 0, 2, &
                                         skala_state%layout%chunk_kin_input_t)
            skala_state%layout%chunk_dynamic_input_views_active = .TRUE.
         END IF
         IF (features%uses_collapsed_rks_dynamic) THEN
            features%density_input_t = skala_state%layout%chunk_density_input_t
            features%grad_input_t = skala_state%layout%chunk_grad_input_t
            features%kin_input_t = skala_state%layout%chunk_kin_input_t
         END IF

         IF (.NOT. skala_state%layout%chunk_inputs_active) THEN
            CALL torch_dict_clone(skala_state%layout%chunk_static_inputs, skala_state%layout%chunk_inputs)
            IF (features%uses_collapsed_rks_dynamic) THEN
               CALL torch_dict_insert(skala_state%layout%chunk_inputs, skala_key_density, &
                                      features%density_input_t)
               CALL torch_dict_insert(skala_state%layout%chunk_inputs, skala_key_grad, &
                                      features%grad_input_t)
               CALL torch_dict_insert(skala_state%layout%chunk_inputs, skala_key_kin, &
                                      features%kin_input_t)
            ELSE
               CALL torch_dict_insert(skala_state%layout%chunk_inputs, skala_key_density, &
                                      skala_state%layout%chunk_density_t)
               CALL torch_dict_insert(skala_state%layout%chunk_inputs, skala_key_grad, &
                                      skala_state%layout%chunk_grad_t)
               CALL torch_dict_insert(skala_state%layout%chunk_inputs, skala_key_kin, &
                                      skala_state%layout%chunk_kin_t)
            END IF
            CALL torch_dict_insert(skala_state%layout%chunk_inputs, skala_key_coarse_0_atomic_coords, &
                                   skala_state%layout%chunk_coarse_0_atomic_coords_t)
            skala_state%layout%chunk_inputs_use_collapsed_rks = features%uses_collapsed_rks_dynamic
            skala_state%layout%chunk_inputs_active = .TRUE.
         END IF
         features%inputs = skala_state%layout%chunk_inputs
         features%owns_inputs = .FALSE.
         features%coarse_0_atomic_coords_t = skala_state%layout%chunk_coarse_0_atomic_coords_t
      ELSE
         CPASSERT(skala_state%layout%static_tensors_active)
         IF (.NOT. requires_stress_grad .AND. .NOT. my_requires_weight_grad) THEN
            features%grid_coords_t = skala_state%layout%grid_coords_t
            features%grid_weights_t = skala_state%layout%grid_weights_t
            features%atomic_grid_weights_t = skala_state%layout%atomic_grid_weights_t
         END IF
         features%atomic_grid_sizes_t = skala_state%layout%atomic_grid_sizes_t
         features%atomic_grid_size_bound_shape_t = skala_state%layout%atomic_grid_size_bound_shape_t

         CALL staged_upload_dynamic(skala_state%layout%density_t, features%density, &
                                    skala_state%layout%grad_t, features%grad, &
                                    skala_state%layout%kin_t, features%kin, &
                                    requires_grad, skala_state%layout%cuda_device, &
                                    dynamic_reallocated)
         features%density_t = skala_state%layout%density_t
         features%grad_t = skala_state%layout%grad_t
         features%kin_t = skala_state%layout%kin_t
         IF (dynamic_reallocated .AND. skala_state%layout%inputs_active) THEN
            CALL torch_dict_release(skala_state%layout%inputs)
            skala_state%layout%inputs_active = .FALSE.
         END IF
         skala_state%layout%dynamic_tensors_active = .TRUE.

         IF (requires_coordinate_grad .OR. requires_stress_grad .OR. my_requires_weight_grad) THEN
            grad_mode = 0
            IF (requires_stress_grad) THEN
               grad_mode = 2
            ELSE IF (my_requires_weight_grad) THEN
               grad_mode = 1
            END IF
            CALL ensure_full_static_tensor_grad(grad_mode)
            IF (requires_stress_grad .OR. my_requires_weight_grad) THEN
               IF (.NOT. skala_state%layout%grad_inputs_active) THEN
                  CALL torch_dict_clone(skala_state%layout%static_inputs, skala_state%layout%grad_inputs)
                  skala_state%layout%grad_inputs_active = .TRUE.
               END IF
               features%grid_coords_t = skala_state%layout%grid_coords_t
               features%grid_weights_t = skala_state%layout%grid_weights_t
               features%atomic_grid_weights_t = skala_state%layout%atomic_grid_weights_t
               CALL torch_dict_set(skala_state%layout%grad_inputs, skala_key_grid_coords, &
                                   features%grid_coords_t)
               CALL torch_dict_set(skala_state%layout%grad_inputs, skala_key_grid_weights, &
                                   features%grid_weights_t)
               CALL torch_dict_set(skala_state%layout%grad_inputs, skala_key_atomic_grid_weights, &
                                   features%atomic_grid_weights_t)
               features%inputs = skala_state%layout%grad_inputs
               features%owns_inputs = .FALSE.
               features%owns_weight_tensors = .FALSE.
               features%owns_grid_coordinate_tensor = .FALSE.
            ELSE
               CALL torch_dict_clone(skala_state%layout%static_inputs, features%inputs)
            END IF
            CALL torch_dict_set(features%inputs, skala_key_density, features%density_t)
            CALL torch_dict_set(features%inputs, skala_key_grad, features%grad_t)
            CALL torch_dict_set(features%inputs, skala_key_kin, features%kin_t)
         ELSE
            CALL ensure_full_static_tensor_grad(0)
            IF (.NOT. skala_state%layout%inputs_active) THEN
               CALL torch_dict_clone(skala_state%layout%static_inputs, skala_state%layout%inputs)
               CALL torch_dict_insert(skala_state%layout%inputs, skala_key_density, &
                                      skala_state%layout%density_t)
               CALL torch_dict_insert(skala_state%layout%inputs, skala_key_grad, skala_state%layout%grad_t)
               CALL torch_dict_insert(skala_state%layout%inputs, skala_key_kin, skala_state%layout%kin_t)
               CALL torch_dict_insert(skala_state%layout%inputs, skala_key_coarse_0_atomic_coords, &
                                      skala_state%layout%coarse_0_atomic_coords_t)
               skala_state%layout%inputs_active = .TRUE.
            END IF
            features%inputs = skala_state%layout%inputs
            features%owns_inputs = .FALSE.
            features%coarse_0_atomic_coords_t = skala_state%layout%coarse_0_atomic_coords_t
         END IF
      END IF

      IF (.NOT. use_atom_chunks) THEN
         IF (requires_coordinate_grad .OR. requires_stress_grad) THEN
            CALL torch_tensor_from_array(features%coarse_0_atomic_coords_t, &
                                         features%coarse_0_atomic_coords)
            CALL torch_tensor_to_device_leaf(features%coarse_0_atomic_coords_t, .TRUE.)
            CALL torch_dict_set(features%inputs, skala_key_coarse_0_atomic_coords, &
                                features%coarse_0_atomic_coords_t)
            features%owns_coordinate_tensor = .TRUE.
         ELSE IF (my_requires_weight_grad) THEN
            CALL torch_dict_set(features%inputs, skala_key_coarse_0_atomic_coords, &
                                skala_state%layout%coarse_0_atomic_coords_t)
         END IF
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
!> \param partition_atom_coords ...
!> \param distances ...
!> \param pair_distances ...
! **************************************************************************************************
   SUBROUTINE smooth_atom_partition(grid_point, atom_coords, cell, weights, partition_atom_coords, &
                                    distances, pair_distances)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: grid_point
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: atom_coords
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: weights
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: partition_atom_coords
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: distances
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: pair_distances

      INTEGER                                            :: iatom, jatom, natom
      REAL(KIND=dp)                                      :: mu, rab, rsum, switch, total
      REAL(KIND=dp), DIMENSION(3)                        :: rij

      natom = SIZE(atom_coords, 2)
      CPASSERT(SIZE(weights) == natom)
      CPASSERT(SIZE(partition_atom_coords, 1) == 3)
      CPASSERT(SIZE(partition_atom_coords, 2) == natom)
      CPASSERT(SIZE(distances) == natom)

      IF (PRESENT(pair_distances)) THEN
         DO iatom = 1, natom
            rij = grid_point - atom_coords(:, iatom)
            distances(iatom) = SQRT(SUM(rij**2))
         END DO
      ELSE
         DO iatom = 1, natom
            partition_atom_coords(:, iatom) = &
               nearest_atom_image_coordinate(atom_coords(:, iatom), grid_point, cell)
            rij = grid_point - partition_atom_coords(:, iatom)
            distances(iatom) = SQRT(SUM(rij**2))
         END DO
      END IF

      weights = 1.0_dp
      DO iatom = 1, natom - 1
         DO jatom = iatom + 1, natom
            IF (PRESENT(pair_distances)) THEN
               rab = pair_distances(iatom, jatom)
            ELSE
               rij = partition_atom_coords(:, iatom) - partition_atom_coords(:, jatom)
               rab = SQRT(SUM(rij**2))
            END IF
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
!> \brief Build normalized smooth atom partition weights.
!> \param grid_point ...
!> \param atom_coords ...
!> \param cell ...
!> \param weights ...
!> \param included ...
!> \param partition_atom_coords_out ...
!> \param distances_out ...
! **************************************************************************************************
   SUBROUTINE smooth_partition_weights(grid_point, atom_coords, cell, weights, included, &
                                       partition_atom_coords_out, distances_out)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: grid_point
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: atom_coords
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: weights
      LOGICAL, DIMENSION(:), INTENT(OUT)                 :: included
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: partition_atom_coords_out
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT), OPTIONAL :: distances_out

      INTEGER                                            :: iatom, jatom, natom
      REAL(KIND=dp)                                      :: included_sum, rsum
      REAL(KIND=dp), DIMENSION(3, SIZE(atom_coords, 2))  :: partition_atom_coords
      REAL(KIND=dp), DIMENSION(SIZE(atom_coords, 2))     :: distances

      natom = SIZE(atom_coords, 2)
      CPASSERT(SIZE(weights) == natom)
      CPASSERT(SIZE(included) == natom)

      CALL smooth_atom_partition(grid_point, atom_coords, cell, weights, &
                                 partition_atom_coords, distances)
      IF (PRESENT(partition_atom_coords_out)) partition_atom_coords_out = partition_atom_coords
      IF (PRESENT(distances_out)) distances_out = distances
      included = weights > smooth_partition_eps
      included_sum = SUM(weights, MASK=included)

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
         IF (included(iatom)) weights(iatom) = weights(iatom)/included_sum
      END DO

   END SUBROUTINE smooth_partition_weights

! **************************************************************************************************
!> \brief Contract the smooth-partition logarithmic weight derivatives directly onto the
!>        force and/or virial, avoiding an O(natom^2) derivative tensor.
!>
!> With L(j,i) = d ln w_i / dR_j, m_j = sum_i w_i L(j,i) and A = sum_i G_i w_i, the force is
!> F_j = sum_i G_i w_i L(j,i) - A m_j. The pair loop accumulates sum_i G_i w_i L(j,i) and m_j
!> without storing L(j,i).
!> \param grid_point ...
!> \param atom_coords ...
!> \param cell ...
!> \param weights ...
!> \param included ...
!> \param weight_grad ...
!> \param atom_force ...
!> \param virial_xc ...
!> \param partition_atom_coords_in ...
!> \param distances_in ...
! **************************************************************************************************
   SUBROUTINE smooth_partition_derivative_action(grid_point, atom_coords, cell, &
                                                 weights, included, weight_grad, &
                                                 atom_force, virial_xc, partition_atom_coords_in, &
                                                 distances_in)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: grid_point
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: atom_coords
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: weights
      LOGICAL, DIMENSION(:), INTENT(IN)                  :: included
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: weight_grad
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL                                        :: atom_force
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: virial_xc
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: partition_atom_coords_in
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: distances_in

      INTEGER                                            :: iatom, idir, jatom, jdir, natom
      LOGICAL                                            :: do_force, do_virial
      REAL(KIND=dp)                                      :: dist_diff, ds_dmu, mu, mu_raw, &
                                                            one_minus_switch, rab, switch, &
                                                            total_weight_grad
      REAL(KIND=dp), DIMENSION(3)                        :: dmu_atom_i, dmu_atom_j, ds_atom_i, &
                                                            ds_atom_j, pair, unit_pair
      REAL(KIND=dp), DIMENSION(3, 3)                     :: dmu_strain, ds_strain, strain_action, &
                                                            strain_mean
      REAL(KIND=dp), DIMENSION(3, SIZE(atom_coords, 2))  :: atom_action, atom_mean, &
                                                            partition_atom_coords, rvecs, &
                                                            unit_rvecs
      REAL(KIND=dp), DIMENSION(SIZE(atom_coords, 2))     :: distances

      natom = SIZE(atom_coords, 2)
      do_force = PRESENT(atom_force)
      do_virial = PRESENT(virial_xc)
      CPASSERT(SIZE(weights) == natom)
      CPASSERT(SIZE(weight_grad) == natom)
      CPASSERT(SIZE(included) == natom)
      IF (do_force) THEN
         CPASSERT(SIZE(atom_force, 1) == 3)
         CPASSERT(SIZE(atom_force, 2) == natom)
      END IF

      DO iatom = 1, natom
         IF (PRESENT(partition_atom_coords_in)) THEN
            partition_atom_coords(:, iatom) = partition_atom_coords_in(:, iatom)
         ELSE
            partition_atom_coords(:, iatom) = &
               nearest_atom_image_coordinate(atom_coords(:, iatom), grid_point, cell)
         END IF
         rvecs(:, iatom) = grid_point - partition_atom_coords(:, iatom)
         IF (PRESENT(distances_in)) THEN
            distances(iatom) = distances_in(iatom)
         ELSE
            distances(iatom) = SQRT(SUM(rvecs(:, iatom)**2))
         END IF
         IF (distances(iatom) > layout_tol) THEN
            unit_rvecs(:, iatom) = rvecs(:, iatom)/distances(iatom)
         ELSE
            unit_rvecs(:, iatom) = 0.0_dp
         END IF
      END DO

      atom_mean = 0.0_dp
      atom_action = 0.0_dp
      strain_mean = 0.0_dp
      strain_action = 0.0_dp
      total_weight_grad = 0.0_dp
      DO iatom = 1, natom
         IF (included(iatom)) total_weight_grad = total_weight_grad + &
                                                  weight_grad(iatom)*weights(iatom)
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
            IF (ABS(ds_dmu) <= 0.0_dp .OR. switch <= TINY(1.0_dp) .OR. &
                one_minus_switch <= TINY(1.0_dp)) CYCLE
            dmu_atom_i = (-unit_rvecs(:, iatom)*rab - dist_diff*unit_pair)/rab**2
            dmu_atom_j = (unit_rvecs(:, jatom)*rab + dist_diff*unit_pair)/rab**2
            ds_atom_i = ds_dmu*dmu_atom_i
            ds_atom_j = ds_dmu*dmu_atom_j

            IF (included(iatom)) THEN
               atom_mean(:, iatom) = atom_mean(:, iatom) + weights(iatom)*ds_atom_i/switch
               atom_mean(:, jatom) = atom_mean(:, jatom) + weights(iatom)*ds_atom_j/switch
               atom_action(:, iatom) = atom_action(:, iatom) + &
                                       weight_grad(iatom)*weights(iatom)*ds_atom_i/switch
               atom_action(:, jatom) = atom_action(:, jatom) + &
                                       weight_grad(iatom)*weights(iatom)*ds_atom_j/switch
            END IF
            IF (included(jatom)) THEN
               atom_mean(:, iatom) = atom_mean(:, iatom) - &
                                     weights(jatom)*ds_atom_i/one_minus_switch
               atom_mean(:, jatom) = atom_mean(:, jatom) - &
                                     weights(jatom)*ds_atom_j/one_minus_switch
               atom_action(:, iatom) = atom_action(:, iatom) - &
                                       weight_grad(jatom)*weights(jatom)*ds_atom_i/one_minus_switch
               atom_action(:, jatom) = atom_action(:, jatom) - &
                                       weight_grad(jatom)*weights(jatom)*ds_atom_j/one_minus_switch
            END IF

            IF (do_virial) THEN
               DO idir = 1, 3
                  DO jdir = 1, 3
                     dmu_strain(idir, jdir) = &
                        ((unit_rvecs(idir, iatom)*rvecs(jdir, iatom) - &
                          unit_rvecs(idir, jatom)*rvecs(jdir, jatom))*rab - &
                         dist_diff*unit_pair(idir)*pair(jdir))/rab**2
                  END DO
               END DO
               ds_strain = ds_dmu*dmu_strain
               IF (included(iatom)) THEN
                  strain_mean = strain_mean + weights(iatom)*ds_strain/switch
                  strain_action = strain_action + &
                                  weight_grad(iatom)*weights(iatom)*ds_strain/switch
               END IF
               IF (included(jatom)) THEN
                  strain_mean = strain_mean - weights(jatom)*ds_strain/one_minus_switch
                  strain_action = strain_action - &
                                  weight_grad(jatom)*weights(jatom)*ds_strain/one_minus_switch
               END IF
            END IF
         END DO
      END DO

      IF (do_force) THEN
         DO jatom = 1, natom
            atom_force(:, jatom) = atom_force(:, jatom) + &
                                   atom_action(:, jatom) - total_weight_grad*atom_mean(:, jatom)
         END DO
      END IF

      IF (do_virial) THEN
         strain_action = strain_action - total_weight_grad*strain_mean
         DO idir = 1, 3
            DO jdir = 1, idir
               IF (idir == jdir) THEN
                  virial_xc(idir, jdir) = virial_xc(idir, jdir) + strain_action(idir, jdir)
               ELSE
                  virial_xc(jdir, idir) = virial_xc(jdir, idir) + strain_action(idir, jdir)
                  virial_xc(idir, jdir) = virial_xc(idir, jdir) + strain_action(idir, jdir)
               END IF
            END DO
         END DO
      END IF

   END SUBROUTINE smooth_partition_derivative_action

! **************************************************************************************************
!> \brief Smoothly suppress a sparse atom row's internal quadrature weight at the layout cutoff.
!> \param weight normalized smooth-partition weight
!> \return scale in the interval [0, 1]
! **************************************************************************************************
   PURE FUNCTION smooth_partition_atomic_weight_scale(weight) RESULT(scale)
      REAL(KIND=dp), INTENT(IN)                          :: weight
      REAL(KIND=dp)                                      :: scale

      REAL(KIND=dp)                                      :: x

      IF (weight <= smooth_partition_eps) THEN
         scale = 0.0_dp
      ELSE IF (weight >= smooth_partition_taper_end) THEN
         scale = 1.0_dp
      ELSE
         x = (weight - smooth_partition_eps)/ &
             (smooth_partition_taper_end - smooth_partition_eps)
         scale = x**3*(10.0_dp + x*(-15.0_dp + 6.0_dp*x))
      END IF

   END FUNCTION smooth_partition_atomic_weight_scale

! **************************************************************************************************
!> \brief Derivative of the sparse atom-row quadrature taper with respect to partition weight.
!> \param weight normalized smooth-partition weight
!> \return derivative of the internal quadrature-weight scale
! **************************************************************************************************
   PURE FUNCTION smooth_partition_atomic_weight_scale_derivative(weight) RESULT(derivative)
      REAL(KIND=dp), INTENT(IN)                          :: weight
      REAL(KIND=dp)                                      :: derivative

      REAL(KIND=dp)                                      :: x

      IF (weight <= smooth_partition_eps .OR. weight >= smooth_partition_taper_end) THEN
         derivative = 0.0_dp
      ELSE
         x = (weight - smooth_partition_eps)/ &
             (smooth_partition_taper_end - smooth_partition_eps)
         derivative = 30.0_dp*x**2*(1.0_dp - x)**2/ &
                      (smooth_partition_taper_end - smooth_partition_eps)
      END IF

   END FUNCTION smooth_partition_atomic_weight_scale_derivative

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

! **************************************************************************************************
!> \brief Allocate the zero-row placeholder shape for atomic_grid_size_bound_shape.
!> \param buf ...
!> \param max_grid_size ...
! **************************************************************************************************
   SUBROUTINE allocate_empty_bound_shape(buf, max_grid_size)
      INTEGER(KIND=int_8), ALLOCATABLE, &
         DIMENSION(:, :), INTENT(OUT)                    :: buf
      INTEGER, INTENT(IN)                                :: max_grid_size

      ALLOCATE (buf(0, max_grid_size))

   END SUBROUTINE allocate_empty_bound_shape

! **************************************************************************************************
!> \brief Ensure a reusable 1D real scratch buffer has the requested size.
!> \param buf ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE ensure_real_1d_alloc(buf, n)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: buf
      INTEGER, INTENT(IN)                                :: n

      IF (ALLOCATED(buf)) THEN
         IF (SIZE(buf) == MAX(1, n)) RETURN
         DEALLOCATE (buf)
      END IF
      ALLOCATE (buf(MAX(1, n)))
      buf = 0.0_dp

   END SUBROUTINE ensure_real_1d_alloc

! **************************************************************************************************
!> \brief Ensure a reusable 1D integer scratch buffer has the requested size.
!> \param buf ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE ensure_int_1d_alloc(buf, n)
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: buf
      INTEGER, INTENT(IN)                                :: n

      IF (ALLOCATED(buf)) THEN
         IF (SIZE(buf) == MAX(1, n)) RETURN
         DEALLOCATE (buf)
      END IF
      ALLOCATE (buf(MAX(1, n)))
      buf = 0

   END SUBROUTINE ensure_int_1d_alloc

! **************************************************************************************************
!> \brief Pointer variant of ensure_real_1d (owns the allocation).
!> \param buf ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE ensure_real_1d_ptr(buf, n)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: buf
      INTEGER, INTENT(IN)                                :: n

      IF (ASSOCIATED(buf)) THEN
         IF (SIZE(buf) == MAX(1, n)) RETURN
         DEALLOCATE (buf)
      END IF
      ALLOCATE (buf(MAX(1, n)))
      buf = 0.0_dp

   END SUBROUTINE ensure_real_1d_ptr

! **************************************************************************************************
!> \brief Pointer variant of ensure_real_2d (owns the allocation).
!> \param buf ...
!> \param n ...
!> \param m ...
! **************************************************************************************************
   SUBROUTINE ensure_real_2d_ptr(buf, n, m)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         POINTER                                         :: buf
      INTEGER, INTENT(IN)                                :: n, m

      IF (ASSOCIATED(buf)) THEN
         IF (SIZE(buf, 1) == n .AND. SIZE(buf, 2) == m) RETURN
         DEALLOCATE (buf)
      END IF
      ALLOCATE (buf(n, m))
      buf = 0.0_dp

   END SUBROUTINE ensure_real_2d_ptr

! **************************************************************************************************
!> \brief Pointer variant of ensure_int_1d (owns the allocation).
!> \param buf ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE ensure_int_1d_ptr(buf, n)
      INTEGER, DIMENSION(:), INTENT(INOUT), POINTER      :: buf
      INTEGER, INTENT(IN)                                :: n

      IF (ASSOCIATED(buf)) THEN
         IF (SIZE(buf) == MAX(1, n)) RETURN
         DEALLOCATE (buf)
      END IF
      ALLOCATE (buf(MAX(1, n)))
      buf = 0

   END SUBROUTINE ensure_int_1d_ptr

! **************************************************************************************************
!> \brief Ensure a reusable 1D int_8 scratch buffer has the requested size.
!> \param buf ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE ensure_i8_1d(buf, n)
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: buf
      INTEGER, INTENT(IN)                                :: n

      IF (ALLOCATED(buf)) THEN
         IF (SIZE(buf) == MAX(1, n)) RETURN
         DEALLOCATE (buf)
      END IF
      ALLOCATE (buf(MAX(1, n)))
      buf = 0_int_8

   END SUBROUTINE ensure_i8_1d

! **************************************************************************************************
!> \brief Ensure a reusable 2D real scratch buffer has the requested shape.
!> \param buf ...
!> \param n ...
!> \param m ...
! **************************************************************************************************
   SUBROUTINE ensure_real_2d_alloc(buf, n, m)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: buf
      INTEGER, INTENT(IN)                                :: n, m

      IF (ALLOCATED(buf)) THEN
         IF (SIZE(buf, 1) == n .AND. SIZE(buf, 2) == m) RETURN
         DEALLOCATE (buf)
      END IF
      ALLOCATE (buf(n, m))
      buf = 0.0_dp

   END SUBROUTINE ensure_real_2d_alloc

! **************************************************************************************************
!> \brief Ensure a reusable 3D real scratch buffer has the requested shape.
!> \param buf ...
!> \param n ...
!> \param m ...
!> \param k ...
! **************************************************************************************************
   SUBROUTINE ensure_real_3d(buf, n, m, k)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: buf
      INTEGER, INTENT(IN)                                :: n, m, k

      IF (ALLOCATED(buf)) THEN
         IF (SIZE(buf, 1) == n .AND. SIZE(buf, 2) == m .AND. SIZE(buf, 3) == k) RETURN
         DEALLOCATE (buf)
      END IF
      ALLOCATE (buf(n, m, k))
      buf = 0.0_dp

   END SUBROUTINE ensure_real_3d

END MODULE skala_gpw_features

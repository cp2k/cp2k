!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Experimental CP2K-native GPW real-space-grid path for SKALA TorchScript models.
! **************************************************************************************************
MODULE skala_gpw_functional
   USE cell_types,                      ONLY: cell_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE input_section_types,             ONLY: section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_comm_type
   USE offload_api,                     ONLY: offload_get_device_count,&
                                              offload_set_chosen_device
   USE particle_types,                  ONLY: particle_type
   USE pw_methods,                      ONLY: pw_scale,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE skala_gpw_features,              ONLY: skala_gpw_feature_build,&
                                              skala_gpw_feature_release,&
                                              skala_gpw_feature_type
   USE skala_torch_api,                 ONLY: skala_torch_model_get_exc,&
                                              skala_torch_model_load,&
                                              skala_torch_model_release,&
                                              skala_torch_model_type
   USE string_utilities,                ONLY: uppercase
   USE torch_api,                       ONLY: torch_cuda_is_available,&
                                              torch_tensor_backward_scalar,&
                                              torch_tensor_data_ptr,&
                                              torch_tensor_grad,&
                                              torch_tensor_release,&
                                              torch_tensor_type,&
                                              torch_use_cuda
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
   USE xc_util,                         ONLY: xc_pw_divergence,&
                                              xc_requires_tmp_g
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_gpw_functional'
   INTEGER, PARAMETER, PRIVATE          :: ngrad_per_point = 10

   PUBLIC :: ensure_native_skala_grid_scope, skala_gpw_eval, xc_section_uses_native_skala_grid

   TYPE(skala_torch_model_type), SAVE                  :: cached_model
   CHARACTER(len=default_path_length), SAVE            :: cached_model_path = ""
   LOGICAL, SAVE                                       :: cached_model_loaded = .FALSE.
   INTEGER, SAVE                                       :: cached_model_cuda_device = -3
   INTEGER, SAVE                                       :: logged_cuda_device = -3, &
                                                          logged_cuda_device_count = -1, &
                                                          logged_cuda_nproc = -1, &
                                                          logged_cuda_request = -3

CONTAINS

! **************************************************************************************************
!> \brief Return true if the GAUXC subsection requests the CP2K-native GPW grid path.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION xc_section_uses_native_skala_grid(xc_section) RESULT(uses_native_grid)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: uses_native_grid

      TYPE(section_vals_type), POINTER                   :: gauxc_section

      uses_native_grid = .FALSE.
      gauxc_section => get_gauxc_section(xc_section)
      IF (ASSOCIATED(gauxc_section)) THEN
         CALL section_vals_val_get(gauxc_section, "NATIVE_GRID", l_val=uses_native_grid)
      END IF

   END FUNCTION xc_section_uses_native_skala_grid

! **************************************************************************************************
!> \brief Enforce the currently implemented native SKALA GPW input scope.
!> \param xc_section ...
! **************************************************************************************************
   SUBROUTINE ensure_native_skala_grid_scope(xc_section)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section

      CHARACTER(len=default_path_length)                 :: model_key, model_name
      CHARACTER(len=default_string_length)               :: functional_key, functional_name
      INTEGER                                            :: ifun, nfun
      LOGICAL                                            :: native_grid
      TYPE(section_vals_type), POINTER                   :: functionals, gauxc_section, xc_fun

      NULLIFY (gauxc_section)
      IF (.NOT. ASSOCIATED(xc_section)) THEN
         CPABORT("Native SKALA GPW requires an XC section")
      END IF

      functionals => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      IF (.NOT. ASSOCIATED(functionals)) THEN
         CPABORT("Native SKALA GPW requires an XC_FUNCTIONAL section")
      END IF

      nfun = 0
      ifun = 0
      DO
         ifun = ifun + 1
         xc_fun => section_vals_get_subs_vals2(functionals, i_section=ifun)
         IF (.NOT. ASSOCIATED(xc_fun)) EXIT
         nfun = nfun + 1
         IF (xc_fun%section%name == "GAUXC") gauxc_section => xc_fun
      END DO

      IF (.NOT. ASSOCIATED(gauxc_section)) THEN
         CPABORT("Native SKALA GPW requires an XC_FUNCTIONAL%GAUXC section")
      END IF
      IF (nfun /= 1) THEN
         CPABORT("Native SKALA GPW requires GAUXC to be the only XC functional")
      END IF

      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID", l_val=native_grid)
      IF (.NOT. native_grid) RETURN

      CALL section_vals_val_get(gauxc_section, "FUNCTIONAL", c_val=functional_name)
      functional_key = ADJUSTL(functional_name)
      CALL uppercase(functional_key)
      IF (TRIM(functional_key) /= "PBE") THEN
         CPABORT("Native SKALA GPW currently requires GAUXC%FUNCTIONAL PBE")
      END IF

      CALL section_vals_val_get(gauxc_section, "MODEL", c_val=model_name)
      model_key = ADJUSTL(model_name)
      CALL uppercase(model_key)
      IF (TRIM(model_key) == "NONE" .OR. TRIM(model_key) == "") THEN
         CPABORT("Native SKALA GPW requires GAUXC%MODEL SKALA or a TorchScript model path")
      END IF

   END SUBROUTINE ensure_native_skala_grid_scope

! **************************************************************************************************
!> \brief Evaluate SKALA energy and first derivatives on a CP2K GPW grid.
!> \param vxc_rho ...
!> \param vxc_tau ...
!> \param exc ...
!> \param rho_r ...
!> \param rho_g ...
!> \param tau ...
!> \param xc_section ...
!> \param weights ...
!> \param pw_pool ...
!> \param particle_set ...
!> \param cell ...
!> \param compute_virial ...
!> \param virial_xc ...
!> \param just_energy ...
!> \param atom_force ...
! **************************************************************************************************
   SUBROUTINE skala_gpw_eval(vxc_rho, vxc_tau, exc, rho_r, rho_g, tau, xc_section, &
                             weights, pw_pool, particle_set, cell, compute_virial, virial_xc, &
                             just_energy, atom_force)
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: vxc_rho, vxc_tau
      REAL(KIND=dp), INTENT(OUT)                         :: exc
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: tau
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN)                                :: compute_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: virial_xc
      LOGICAL, INTENT(IN), OPTIONAL                      :: just_energy
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: atom_force

      CHARACTER(len=default_path_length)                 :: model_path
      INTEGER                                            :: native_grid_cuda_device, nspins, &
                                                            phase_handle, selected_cuda_device, &
                                                            xc_deriv_method_id, xc_rho_smooth_id
      LOGICAL :: lsd, my_just_energy, native_grid_atom_chunk_routing, native_grid_atom_chunks, &
         native_grid_diagnostics, native_grid_use_cuda, needs_atom_force
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: density_grad, kin_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: grad_grad
      TYPE(section_vals_type), POINTER                   :: gauxc_section
      TYPE(skala_gpw_feature_type)                       :: features
      TYPE(torch_tensor_type)                            :: atom_coord_grad_t, exc_tensor
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set

      virial_xc = 0.0_dp
      exc = 0.0_dp
      my_just_energy = .FALSE.
      IF (PRESENT(just_energy)) my_just_energy = just_energy
      needs_atom_force = PRESENT(atom_force)
      IF (needs_atom_force) atom_force = 0.0_dp

      IF (compute_virial) THEN
         CALL cp_abort(__LOCATION__, &
                       "Native SKALA GPW stress/virial is not implemented yet.")
      END IF
      IF (.NOT. ASSOCIATED(rho_g)) THEN
         CALL cp_abort(__LOCATION__, &
                       "Native SKALA GPW requires the reciprocal-space density to form density gradients.")
      END IF
      IF (.NOT. ASSOCIATED(tau)) THEN
         CALL cp_abort(__LOCATION__, &
                       "Native SKALA GPW requires the kinetic-energy density.")
      END IF

      nspins = SIZE(rho_r)
      lsd = (nspins /= 1)
      CALL get_skala_model_path(xc_section, model_path)
      gauxc_section => get_gauxc_section(xc_section)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_USE_CUDA", l_val=native_grid_use_cuda)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_CUDA_DEVICE", &
                                i_val=native_grid_cuda_device)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_ATOM_CHUNKS", &
                                l_val=native_grid_atom_chunks)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_ATOM_CHUNK_ROUTING", &
                                l_val=native_grid_atom_chunk_routing)
      native_grid_atom_chunk_routing = native_grid_atom_chunk_routing .OR. native_grid_atom_chunks
      native_grid_atom_chunks = native_grid_atom_chunks .OR. native_grid_atom_chunk_routing
      IF (native_grid_atom_chunks .AND. needs_atom_force) THEN
         CALL cp_abort(__LOCATION__, &
                       "Native SKALA GPW atom chunks are not implemented for atom forces yet.")
      END IF
      ! The portable SKALA export used by the regtests builds ragged-index tensors on CPU.
      CALL torch_use_cuda(native_grid_use_cuda)
      selected_cuda_device = configure_native_grid_cuda( &
                             native_grid_use_cuda, native_grid_cuda_device, rho_r(1)%pw_grid%para%group)
      CALL ensure_model_loaded(model_path, selected_cuda_device)

      IF (lsd) THEN
         needs%rho_spin = .TRUE.
         needs%drho_spin = .TRUE.
         needs%tau_spin = .TRUE.
      ELSE
         needs%rho = .TRUE.
         needs%drho = .TRUE.
         needs%tau = .TRUE.
      END IF

      CALL section_vals_val_get(xc_section, "XC_GRID%XC_DERIV", i_val=xc_deriv_method_id)
      CALL section_vals_val_get(xc_section, "XC_GRID%XC_SMOOTH_RHO", i_val=xc_rho_smooth_id)

      CALL xc_rho_set_create(rho_set, &
                             rho_r(1)%pw_grid%bounds_local, &
                             rho_cutoff=section_get_rval(xc_section, "density_cutoff"), &
                             drho_cutoff=section_get_rval(xc_section, "gradient_cutoff"), &
                             tau_cutoff=section_get_rval(xc_section, "tau_cutoff"))
      CALL xc_rho_set_update(rho_set, rho_r, rho_g, tau, needs, &
                             xc_deriv_method_id, xc_rho_smooth_id, pw_pool)

      CALL skala_gpw_feature_build(features, rho_set, rho_r, particle_set, cell, &
                                   requires_grad=(.NOT. my_just_energy), weights=weights, &
                                   requires_coordinate_grad=needs_atom_force, &
                                   use_atom_chunks=native_grid_atom_chunks, &
                                   route_atom_chunks=native_grid_atom_chunk_routing)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_DIAGNOSTICS", l_val=native_grid_diagnostics)
      IF (native_grid_diagnostics) THEN
         CALL print_native_grid_diagnostics(features, rho_r(1)%pw_grid%para%group%mepos == 0)
      END IF
      CALL skala_torch_model_get_exc(cached_model, features%inputs, &
                                     features%grid_weights_t, exc_tensor, exc)
      IF (features%uses_atom_chunks) CALL rho_r(1)%pw_grid%para%group%sum(exc)

      IF (.NOT. my_just_energy) THEN
         CALL timeset("skala_gpw_backward", phase_handle)
         CALL torch_tensor_backward_scalar(exc_tensor)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_fetch", phase_handle)
         IF (features%uses_atom_chunks) THEN
            CALL fetch_and_gather_atom_chunk_grads(features, rho_r(1)%pw_grid%para%group, &
                                                   density_grad, grad_grad, kin_grad)
         ELSE
            CALL fetch_local_feature_grads(features, density_grad, grad_grad, kin_grad)
         END IF
         IF (needs_atom_force) THEN
            CALL add_explicit_coordinate_force(atom_force, features, atom_coord_grad_t, &
                                               rho_r(1)%pw_grid%para%group%mepos == 0)
         END IF
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_vxc_unpack", phase_handle)
         CALL build_vxc_from_feature_grads(vxc_rho, vxc_tau, rho_r, pw_pool, &
                                           density_grad, grad_grad, kin_grad, &
                                           xc_deriv_method_id)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_release", phase_handle)
         DEALLOCATE (density_grad, grad_grad, kin_grad)
         IF (needs_atom_force) CALL torch_tensor_release(atom_coord_grad_t)
         CALL timestop(phase_handle)
      END IF

      CALL timeset("skala_gpw_cleanup", phase_handle)
      CALL torch_tensor_release(exc_tensor)
      CALL skala_gpw_feature_release(features)
      CALL xc_rho_set_release(rho_set, pw_pool=pw_pool)
      CALL torch_use_cuda(.TRUE.)
      CALL timestop(phase_handle)

   END SUBROUTINE skala_gpw_eval

! **************************************************************************************************
!> \brief Add the explicit SKALA derivative with respect to atom-center coordinates.
!> \param atom_force ...
!> \param features ...
!> \param atom_coord_grad_t ...
!> \param root_rank ...
! **************************************************************************************************
   SUBROUTINE add_explicit_coordinate_force(atom_force, features, atom_coord_grad_t, root_rank)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: atom_force
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(torch_tensor_type), INTENT(INOUT)             :: atom_coord_grad_t
      LOGICAL, INTENT(IN)                                :: root_rank

      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: atom_coord_grad

      NULLIFY (atom_coord_grad)
      CALL torch_tensor_grad(features%coarse_0_atomic_coords_t, atom_coord_grad_t)
      IF (root_rank) THEN
         CALL torch_tensor_data_ptr(atom_coord_grad_t, atom_coord_grad)
         CPASSERT(SIZE(atom_force, 1) == SIZE(atom_coord_grad, 1))
         CPASSERT(SIZE(atom_force, 2) == SIZE(atom_coord_grad, 2))
         atom_force(:, :) = atom_force(:, :) + atom_coord_grad(:, :)
      END IF

   END SUBROUTINE add_explicit_coordinate_force

! **************************************************************************************************
!> \brief Map full Torch feature gradients back to this rank's local grid order.
!> \param features ...
!> \param density_grad ...
!> \param grad_grad ...
!> \param kin_grad ...
! **************************************************************************************************
   SUBROUTINE fetch_local_feature_grads(features, density_grad, grad_grad, kin_grad)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: density_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: grad_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: kin_grad

      INTEGER                                            :: i, j, k, local_row, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: density_grad_all, kin_grad_all
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grad_grad_all
      TYPE(torch_tensor_type)                            :: density_grad_t, grad_grad_t, kin_grad_t

      NULLIFY (density_grad_all, grad_grad_all, kin_grad_all)
      CALL get_feature_grad_views(features, density_grad_t, grad_grad_t, kin_grad_t, &
                                  density_grad_all, grad_grad_all, kin_grad_all)
      CPASSERT(SIZE(density_grad_all, 1) == features%nflat)
      CPASSERT(SIZE(density_grad_all, 2) == 2)
      CPASSERT(SIZE(grad_grad_all, 1) == features%nflat)
      CPASSERT(SIZE(grad_grad_all, 2) == 3)
      CPASSERT(SIZE(grad_grad_all, 3) == 2)
      CPASSERT(SIZE(kin_grad_all, 1) == features%nflat)
      CPASSERT(SIZE(kin_grad_all, 2) == 2)

      ALLOCATE (density_grad(features%nflat_local, 2), &
                grad_grad(features%nflat_local, 3, 2), &
                kin_grad(features%nflat_local, 2))
      local_row = 0
      DO k = LBOUND(features%feature_index, 3), UBOUND(features%feature_index, 3)
         DO j = LBOUND(features%feature_index, 2), UBOUND(features%feature_index, 2)
            DO i = LBOUND(features%feature_index, 1), UBOUND(features%feature_index, 1)
               local_row = local_row + 1
               row = features%feature_index(i, j, k)
               CPASSERT(row >= 1 .AND. row <= features%nflat)
               density_grad(local_row, :) = density_grad_all(row, :)
               grad_grad(local_row, :, :) = grad_grad_all(row, :, :)
               kin_grad(local_row, :) = kin_grad_all(row, :)
            END DO
         END DO
      END DO
      CPASSERT(local_row == features%nflat_local)

      CALL torch_tensor_release(density_grad_t)
      CALL torch_tensor_release(grad_grad_t)
      CALL torch_tensor_release(kin_grad_t)

   END SUBROUTINE fetch_local_feature_grads

! **************************************************************************************************
!> \brief Pack atom-chunk Torch gradients into CP2K communication buffers.
!> \param features ...
!> \param TARGET ...
!> \param route_to_return_positions ...
! **************************************************************************************************
   SUBROUTINE pack_atom_chunk_grads(features, TARGET, route_to_return_positions)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: target
      LOGICAL, INTENT(IN)                                :: route_to_return_positions

      INTEGER                                            :: base, irow, point_pos
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: chunk_density_grad, chunk_kin_grad
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: chunk_grad_grad
      TYPE(torch_tensor_type)                            :: density_grad_t, grad_grad_t, kin_grad_t

      NULLIFY (chunk_density_grad, chunk_grad_grad, chunk_kin_grad)
      CALL get_feature_grad_views(features, density_grad_t, grad_grad_t, kin_grad_t, &
                                  chunk_density_grad, chunk_grad_grad, chunk_kin_grad)
      CPASSERT(SIZE(TARGET) == ngrad_per_point*features%chunk_feature_count)
      CPASSERT(SIZE(chunk_density_grad, 1) == features%chunk_feature_count)
      CPASSERT(SIZE(chunk_density_grad, 2) == 2)
      CPASSERT(SIZE(chunk_grad_grad, 1) == features%chunk_feature_count)
      CPASSERT(SIZE(chunk_grad_grad, 2) == 3)
      CPASSERT(SIZE(chunk_grad_grad, 3) == 2)
      CPASSERT(SIZE(chunk_kin_grad, 1) == features%chunk_feature_count)
      CPASSERT(SIZE(chunk_kin_grad, 2) == 2)

      TARGET = 0.0_dp
      DO irow = 1, features%chunk_feature_count
         IF (route_to_return_positions) THEN
            point_pos = features%chunk_return_positions(irow)
            CPASSERT(point_pos >= 1 .AND. point_pos <= features%chunk_feature_count)
         ELSE
            point_pos = irow
         END IF
         base = ngrad_per_point*(point_pos - 1)
         TARGET(base + 1:base + 2) = chunk_density_grad(irow, :)
         TARGET(base + 3) = chunk_grad_grad(irow, 1, 1)
         TARGET(base + 4) = chunk_grad_grad(irow, 2, 1)
         TARGET(base + 5) = chunk_grad_grad(irow, 3, 1)
         TARGET(base + 6) = chunk_grad_grad(irow, 1, 2)
         TARGET(base + 7) = chunk_grad_grad(irow, 2, 2)
         TARGET(base + 8) = chunk_grad_grad(irow, 3, 2)
         TARGET(base + 9:base + 10) = chunk_kin_grad(irow, :)
      END DO

      CALL torch_tensor_release(density_grad_t)
      CALL torch_tensor_release(grad_grad_t)
      CALL torch_tensor_release(kin_grad_t)

   END SUBROUTINE pack_atom_chunk_grads

! **************************************************************************************************
!> \brief Return CPU views of autograd outputs for the SKALA dynamic feature tensors.
!> \param features ...
!> \param density_grad_t ...
!> \param grad_grad_t ...
!> \param kin_grad_t ...
!> \param density_grad ...
!> \param grad_grad ...
!> \param kin_grad ...
! **************************************************************************************************
   SUBROUTINE get_feature_grad_views(features, density_grad_t, grad_grad_t, kin_grad_t, &
                                     density_grad, grad_grad, kin_grad)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(torch_tensor_type), INTENT(INOUT)             :: density_grad_t, grad_grad_t, kin_grad_t
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: density_grad
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grad_grad
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: kin_grad

      NULLIFY (density_grad, grad_grad, kin_grad)
      CALL torch_tensor_grad(features%density_t, density_grad_t)
      CALL torch_tensor_grad(features%grad_t, grad_grad_t)
      CALL torch_tensor_grad(features%kin_t, kin_grad_t)
      CALL torch_tensor_data_ptr(density_grad_t, density_grad)
      CALL torch_tensor_data_ptr(grad_grad_t, grad_grad)
      CALL torch_tensor_data_ptr(kin_grad_t, kin_grad)

   END SUBROUTINE get_feature_grad_views

! **************************************************************************************************
!> \brief Fetch atom-chunk gradients and route them back to their local grid owners.
!> \param features ...
!> \param group ...
!> \param density_grad ...
!> \param grad_grad ...
!> \param kin_grad ...
! **************************************************************************************************
   SUBROUTINE fetch_and_gather_atom_chunk_grads(features, group, density_grad, grad_grad, &
                                                kin_grad)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: density_grad, kin_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: grad_grad

      INTEGER                                            :: base, i, j, k, local_row, &
                                                            nflat_local, phase_handle, &
                                                            point_pos, row
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: chunk_grad_buffer, global_grad_buffer, &
                                                            recv_grad_buffer, send_grad_buffer

      CPASSERT(features%uses_atom_chunks)
      CPASSERT(features%chunk_feature_count > 0)

      nflat_local = features%nflat_local
      IF (features%uses_atom_chunk_routing) THEN
         CPASSERT(SUM(features%route_point_recv_counts) == features%chunk_feature_count)
         CPASSERT(SUM(features%route_point_send_counts) == nflat_local)

         ALLOCATE (send_grad_buffer(ngrad_per_point*features%chunk_feature_count), &
                   recv_grad_buffer(ngrad_per_point*nflat_local))

         CALL timeset("skala_gpw_grad_torch_pack", phase_handle)
         CALL pack_atom_chunk_grads(features, send_grad_buffer, .TRUE.)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_route_comm", phase_handle)
         CALL group%alltoall(send_grad_buffer, features%route_grad_return_send_counts, &
                             features%route_grad_return_send_displs, recv_grad_buffer, &
                             features%route_grad_return_recv_counts, &
                             features%route_grad_return_recv_displs)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_route_scatter", phase_handle)
         ALLOCATE (density_grad(nflat_local, 2), grad_grad(nflat_local, 3, 2), &
                   kin_grad(nflat_local, 2))
         density_grad = 0.0_dp
         grad_grad = 0.0_dp
         kin_grad = 0.0_dp
         DO point_pos = 1, nflat_local
            local_row = features%route_send_local_rows(point_pos)
            CPASSERT(local_row >= 1 .AND. local_row <= nflat_local)
            base = ngrad_per_point*(point_pos - 1)
            density_grad(local_row, :) = recv_grad_buffer(base + 1:base + 2)
            grad_grad(local_row, 1, 1) = recv_grad_buffer(base + 3)
            grad_grad(local_row, 2, 1) = recv_grad_buffer(base + 4)
            grad_grad(local_row, 3, 1) = recv_grad_buffer(base + 5)
            grad_grad(local_row, 1, 2) = recv_grad_buffer(base + 6)
            grad_grad(local_row, 2, 2) = recv_grad_buffer(base + 7)
            grad_grad(local_row, 3, 2) = recv_grad_buffer(base + 8)
            kin_grad(local_row, :) = recv_grad_buffer(base + 9:base + 10)
         END DO
         CALL timestop(phase_handle)

         DEALLOCATE (recv_grad_buffer, send_grad_buffer)
      ELSE
         ALLOCATE (chunk_grad_buffer(ngrad_per_point*features%chunk_feature_count), &
                   global_grad_buffer(ngrad_per_point*features%nflat))
         CALL timeset("skala_gpw_grad_torch_pack", phase_handle)
         CALL pack_atom_chunk_grads(features, chunk_grad_buffer, .FALSE.)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_allgatherv", phase_handle)
         CALL group%allgatherv(chunk_grad_buffer, global_grad_buffer, &
                               features%chunk_grad_counts, features%chunk_grad_displs)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_scatter", phase_handle)
         ALLOCATE (density_grad(nflat_local, 2), grad_grad(nflat_local, 3, 2), &
                   kin_grad(nflat_local, 2))
         local_row = 0
         DO k = LBOUND(features%feature_index, 3), UBOUND(features%feature_index, 3)
            DO j = LBOUND(features%feature_index, 2), UBOUND(features%feature_index, 2)
               DO i = LBOUND(features%feature_index, 1), UBOUND(features%feature_index, 1)
                  local_row = local_row + 1
                  row = features%feature_index(i, j, k)
                  CPASSERT(row >= 1 .AND. row <= features%nflat)
                  base = ngrad_per_point*(row - 1)
                  density_grad(local_row, :) = global_grad_buffer(base + 1:base + 2)
                  grad_grad(local_row, 1, 1) = global_grad_buffer(base + 3)
                  grad_grad(local_row, 2, 1) = global_grad_buffer(base + 4)
                  grad_grad(local_row, 3, 1) = global_grad_buffer(base + 5)
                  grad_grad(local_row, 1, 2) = global_grad_buffer(base + 6)
                  grad_grad(local_row, 2, 2) = global_grad_buffer(base + 7)
                  grad_grad(local_row, 3, 2) = global_grad_buffer(base + 8)
                  kin_grad(local_row, :) = global_grad_buffer(base + 9:base + 10)
               END DO
            END DO
         END DO
         CALL timestop(phase_handle)
         DEALLOCATE (chunk_grad_buffer, global_grad_buffer)

      END IF

   END SUBROUTINE fetch_and_gather_atom_chunk_grads

! **************************************************************************************************
!> \brief Fill CP2K VXC real-space arrays from Torch feature gradients.
!> \param vxc_rho ...
!> \param vxc_tau ...
!> \param rho_r ...
!> \param pw_pool ...
!> \param density_grad ...
!> \param grad_grad ...
!> \param kin_grad ...
!> \param xc_deriv_method_id ...
! **************************************************************************************************
   SUBROUTINE build_vxc_from_feature_grads(vxc_rho, vxc_tau, rho_r, pw_pool, &
                                           density_grad, grad_grad, kin_grad, &
                                           xc_deriv_method_id)
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: vxc_rho, vxc_tau, rho_r
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: density_grad
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: grad_grad
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: kin_grad
      INTEGER, INTENT(IN)                                :: xc_deriv_method_id

      INTEGER                                            :: i, ipt, ispin, j, k, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(KIND=dp)                                      :: dvol_inv
      TYPE(pw_c1d_gs_type)                               :: tmp_g, vxc_g
      TYPE(pw_r3d_rs_type), DIMENSION(3)                 :: grad_pw

      nspins = SIZE(rho_r)
      bo = rho_r(1)%pw_grid%bounds_local
      dvol_inv = 1.0_dp/rho_r(1)%pw_grid%dvol

      ALLOCATE (vxc_rho(nspins), vxc_tau(nspins))
      DO ispin = 1, nspins
         CALL pw_pool%create_pw(vxc_rho(ispin))
         CALL pw_pool%create_pw(vxc_tau(ispin))
         CALL pw_zero(vxc_rho(ispin))
         CALL pw_zero(vxc_tau(ispin))
      END DO

      IF (xc_requires_tmp_g(xc_deriv_method_id) .OR. rho_r(1)%pw_grid%spherical) THEN
         CALL pw_pool%create_pw(vxc_g)
         IF (.NOT. rho_r(1)%pw_grid%spherical) CALL pw_pool%create_pw(tmp_g)
      END IF

      DO ispin = 1, nspins
         DO i = 1, 3
            CALL pw_pool%create_pw(grad_pw(i))
            CALL pw_zero(grad_pw(i))
         END DO

         ipt = 0
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  ipt = ipt + 1
                  IF (nspins == 1) THEN
                     vxc_rho(1)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                 (density_grad(ipt, 1) + density_grad(ipt, 2))
                     vxc_tau(1)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                 (kin_grad(ipt, 1) + kin_grad(ipt, 2))
                     grad_pw(1)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                 (grad_grad(ipt, 1, 1) + grad_grad(ipt, 1, 2))
                     grad_pw(2)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                 (grad_grad(ipt, 2, 1) + grad_grad(ipt, 2, 2))
                     grad_pw(3)%array(i, j, k) = 0.5_dp*dvol_inv* &
                                                 (grad_grad(ipt, 3, 1) + grad_grad(ipt, 3, 2))
                  ELSE
                     vxc_rho(ispin)%array(i, j, k) = dvol_inv*density_grad(ipt, ispin)
                     vxc_tau(ispin)%array(i, j, k) = dvol_inv*kin_grad(ipt, ispin)
                     grad_pw(1)%array(i, j, k) = dvol_inv*grad_grad(ipt, 1, ispin)
                     grad_pw(2)%array(i, j, k) = dvol_inv*grad_grad(ipt, 2, ispin)
                     grad_pw(3)%array(i, j, k) = dvol_inv*grad_grad(ipt, 3, ispin)
                  END IF
               END DO
            END DO
         END DO

         DO i = 1, 3
            CALL pw_scale(grad_pw(i), -1.0_dp)
         END DO
         CALL xc_pw_divergence(xc_deriv_method_id, grad_pw, tmp_g, vxc_g, vxc_rho(ispin))

         DO i = 1, 3
            CALL pw_pool%give_back_pw(grad_pw(i))
         END DO
      END DO

      IF (ASSOCIATED(vxc_g%pw_grid)) CALL pw_pool%give_back_pw(vxc_g)
      IF (ASSOCIATED(tmp_g%pw_grid)) CALL pw_pool%give_back_pw(tmp_g)

   END SUBROUTINE build_vxc_from_feature_grads

! **************************************************************************************************
!> \brief Print optional diagnostics for the CP2K-native SKALA GPW feature block.
!> \param features ...
!> \param print_active ...
! **************************************************************************************************
   SUBROUTINE print_native_grid_diagnostics(features, print_active)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      LOGICAL, INTENT(IN)                                :: print_active

      INTEGER                                            :: iw

      IF (.NOT. print_active) RETURN

      iw = cp_logger_get_default_io_unit()
      WRITE (UNIT=iw, FMT="(/,T2,A,1X,ES19.11)") &
         "SKALA_GPW| Native grid feature electrons", features%electron_count
      WRITE (UNIT=iw, FMT="(T2,A,1X,ES19.11)") &
         "SKALA_GPW| Native grid feature spin moment", features%spin_moment
      WRITE (UNIT=iw, FMT="(T2,A,1X,ES19.11)") &
         "SKALA_GPW| Native grid feature weight sum", features%grid_weight_sum
      IF (features%uses_atom_chunks) THEN
         WRITE (UNIT=iw, FMT="(T2,A,1X,I0,1X,A,1X,I0)") &
            "SKALA_GPW| Native grid atom chunk rows", features%chunk_feature_count, &
            "of", features%nflat
      END IF

   END SUBROUTINE print_native_grid_diagnostics

! **************************************************************************************************
!> \brief Configure CUDA device selection for the native SKALA GPW Torch path.
!> \param use_cuda ...
!> \param requested_device ...
!> \param group ...
!> \return selected CUDA device, or -1 for CPU fallback/no visible CUDA device
! **************************************************************************************************
   FUNCTION configure_native_grid_cuda(use_cuda, requested_device, group) RESULT(selected_device)
      LOGICAL, INTENT(IN)                                :: use_cuda
      INTEGER, INTENT(IN)                                :: requested_device

      CLASS(mp_comm_type), INTENT(IN)                    :: group

      INTEGER                                            :: cuda_device_count, iw, pe, selected_device
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: selected_devices

      selected_device = -1

      IF (.NOT. use_cuda) RETURN

      IF (.NOT. torch_cuda_is_available()) THEN
         cuda_device_count = 0
      ELSE
         cuda_device_count = offload_get_device_count()
      END IF
      IF (cuda_device_count > 0) THEN
         IF (requested_device < 0) THEN
            selected_device = MOD(group%mepos, cuda_device_count)
         ELSE
            selected_device = requested_device
         END IF
      END IF
      IF (selected_device >= cuda_device_count) THEN
         CALL cp_abort(__LOCATION__, &
                       "GAUXC%NATIVE_GRID_CUDA_DEVICE selects a CUDA device outside the visible "// &
                       "CP2K offload device range.")
      END IF
      IF (selected_device >= 0) CALL offload_set_chosen_device(selected_device)

      ALLOCATE (selected_devices(group%num_pe))
      CALL group%allgather(selected_device, selected_devices)

      IF (group%mepos /= 0) RETURN
      IF (selected_device == logged_cuda_device .AND. &
          cuda_device_count == logged_cuda_device_count .AND. &
          group%num_pe == logged_cuda_nproc .AND. &
          requested_device == logged_cuda_request) RETURN

      iw = cp_logger_get_default_io_unit()
      IF (selected_device >= 0) THEN
         WRITE (UNIT=iw, FMT="(/,T2,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)") &
            "SKALA_GPW| Native grid Torch CUDA device", selected_device, &
            "of", cuda_device_count, "requested", requested_device
      ELSE
         WRITE (UNIT=iw, FMT="(/,T2,A)") &
            "SKALA_GPW| Native grid Torch CUDA requested, but no CP2K offload device is visible"
      END IF
      WRITE (UNIT=iw, FMT="(T2,A)", ADVANCE="NO") &
         "SKALA_GPW| Native grid Torch CUDA rank devices"
      DO pe = 1, group%num_pe
         WRITE (UNIT=iw, FMT="(1X,I0,A,I0)", ADVANCE="NO") pe - 1, ":", selected_devices(pe)
      END DO
      WRITE (UNIT=iw, FMT=*)

      logged_cuda_device = selected_device
      logged_cuda_device_count = cuda_device_count
      logged_cuda_nproc = group%num_pe
      logged_cuda_request = requested_device

   END FUNCTION configure_native_grid_cuda

! **************************************************************************************************
!> \brief Load and cache the TorchScript SKALA model.
!> \param model_path ...
!> \param cuda_device ...
! **************************************************************************************************
   SUBROUTINE ensure_model_loaded(model_path, cuda_device)
      CHARACTER(len=*), INTENT(IN)                       :: model_path
      INTEGER, INTENT(IN)                                :: cuda_device

      IF (cached_model_loaded) THEN
         IF (TRIM(cached_model_path) == TRIM(model_path) .AND. &
             cached_model_cuda_device == cuda_device) RETURN
         CALL skala_torch_model_release(cached_model)
         cached_model_loaded = .FALSE.
      END IF

      CALL skala_torch_model_load(cached_model, TRIM(model_path))
      cached_model_path = model_path
      cached_model_cuda_device = cuda_device
      cached_model_loaded = .TRUE.

   END SUBROUTINE ensure_model_loaded

! **************************************************************************************************
!> \brief Resolve the SKALA TorchScript model path from the GAUXC subsection.
!> \param xc_section ...
!> \param model_path ...
! **************************************************************************************************
   SUBROUTINE get_skala_model_path(xc_section, model_path)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      CHARACTER(len=default_path_length), INTENT(OUT)    :: model_path

      CHARACTER(len=default_path_length)                 :: model_key
      INTEGER                                            :: env_status
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      gauxc_section => get_gauxc_section(xc_section)
      IF (.NOT. ASSOCIATED(gauxc_section)) THEN
         CPABORT("Native SKALA GPW requires an XC_FUNCTIONAL%GAUXC section")
      END IF

      CALL section_vals_val_get(gauxc_section, "MODEL", c_val=model_path)
      model_key = ADJUSTL(model_path)
      CALL uppercase(model_key)
      IF (TRIM(model_key) == "NONE" .OR. TRIM(model_key) == "") THEN
         CPABORT("Native SKALA GPW requires GAUXC%MODEL SKALA or a TorchScript model path")
      ELSE IF (TRIM(model_key) == "SKALA") THEN
         CALL GET_ENVIRONMENT_VARIABLE("GAUXC_SKALA_MODEL", model_path, STATUS=env_status)
         IF (env_status /= 0 .OR. LEN_TRIM(model_path) == 0) THEN
            CPABORT("MODEL SKALA requires the GAUXC_SKALA_MODEL environment variable")
         END IF
      END IF

   END SUBROUTINE get_skala_model_path

! **************************************************************************************************
!> \brief Return the first GAUXC functional subsection, if present.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_gauxc_section(xc_section) RESULT(gauxc_section)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      INTEGER                                            :: ifun
      TYPE(section_vals_type), POINTER                   :: functionals, xc_fun

      NULLIFY (gauxc_section)
      IF (.NOT. ASSOCIATED(xc_section)) RETURN

      functionals => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      IF (.NOT. ASSOCIATED(functionals)) RETURN

      ifun = 0
      DO
         ifun = ifun + 1
         xc_fun => section_vals_get_subs_vals2(functionals, i_section=ifun)
         IF (.NOT. ASSOCIATED(xc_fun)) EXIT
         IF (xc_fun%section%name == "GAUXC") THEN
            gauxc_section => xc_fun
            EXIT
         END IF
      END DO

   END FUNCTION get_gauxc_section

END MODULE skala_gpw_functional

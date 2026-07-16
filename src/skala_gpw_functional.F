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
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_array_utils,                  ONLY: cp_3d_r_cp_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE input_section_types,             ONLY: section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_comm_type
   USE offload_api,                     ONLY: offload_set_chosen_device
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_scale,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE skala_gpw_features,              ONLY: skala_gpw_atom_partition_hard,&
                                              skala_gpw_atom_partition_smooth,&
                                              skala_gpw_atom_subchunk_count,&
                                              skala_gpw_feature_build,&
                                              skala_gpw_feature_build_atom_subchunk,&
                                              skala_gpw_feature_release,&
                                              skala_gpw_feature_type,&
                                              skala_gpw_smooth_partition_derivatives
   USE skala_torch_api,                 ONLY: skala_torch_model_get_exc,&
                                              skala_torch_model_get_exc_density,&
                                              skala_torch_model_load,&
                                              skala_torch_model_release,&
                                              skala_torch_model_type
   USE string_utilities,                ONLY: uppercase
   USE torch_api,                       ONLY: &
        torch_cuda_device_count, torch_cuda_is_available, torch_dict_create, torch_dict_insert, &
        torch_dict_release, torch_dict_type, torch_tensor_backward_scalar, torch_tensor_data_ptr, &
        torch_tensor_from_array, torch_tensor_grad, torch_tensor_release, &
        torch_tensor_to_device_leaf, torch_tensor_type, torch_use_cuda
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_get,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
   USE xc_util,                         ONLY: xc_pw_divergence,&
                                              xc_requires_tmp_g
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'skala_gpw_functional'
   INTEGER, PARAMETER, PRIVATE          :: atom_chunk_auto_max_rows = 400000, &
                                           atom_chunk_auto_min_rows = 100000, &
                                           atom_chunk_auto_row_quantum = 100000, &
                                           ncollapsed_grad_per_point = 5, ngrad_per_point = 10
   INTEGER, PARAMETER, PUBLIC           :: skala_gapw_density_partition_hard_minus_soft = 1, &
                                           skala_gapw_density_partition_hard_only = 2, &
                                           skala_gapw_density_partition_soft_only = 3, &
                                           skala_gapw_density_partition_none = 4

   PUBLIC :: ensure_native_skala_grid_scope, get_gauxc_section, skala_gapw_atom_vxc_of_r, &
             native_skala_gapw_density_partition, skala_gpw_eval, skala_gpw_exc_density, &
             xc_section_uses_native_skala_grid, xc_section_uses_gauxc_model

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
!> \brief Return true if the GAUXC subsection requests a model evaluation.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION xc_section_uses_gauxc_model(xc_section) RESULT(uses_gauxc_model)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: uses_gauxc_model

      CHARACTER(len=default_path_length)                 :: model_key, model_name, xc_key, xc_name
      TYPE(section_vals_type), POINTER                   :: gauxc_section

      uses_gauxc_model = .FALSE.
      gauxc_section => get_gauxc_section(xc_section)
      IF (ASSOCIATED(gauxc_section)) THEN
         CALL section_vals_val_get(gauxc_section, "MODEL", c_val=model_name)
         CALL section_vals_val_get(gauxc_section, "FUNCTIONAL", c_val=xc_name)
         model_key = ADJUSTL(model_name)
         xc_key = ADJUSTL(xc_name)
         CALL uppercase(model_key)
         CALL uppercase(xc_key)
         uses_gauxc_model = (TRIM(model_key) /= "" .AND. TRIM(model_key) /= "NONE" .AND. &
                             TRIM(model_key) /= TRIM(xc_key))
      END IF

   END FUNCTION xc_section_uses_gauxc_model

! **************************************************************************************************
!> \brief Return the hard/soft GAPW one-center density partition for native SKALA.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION native_skala_gapw_density_partition(xc_section) RESULT(partition)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      INTEGER                                            :: partition

      TYPE(section_vals_type), POINTER                   :: gauxc_section

      partition = skala_gapw_density_partition_hard_minus_soft
      gauxc_section => get_gauxc_section(xc_section)
      IF (ASSOCIATED(gauxc_section)) THEN
         CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_GAPW_DENSITY_PARTITION", &
                                   i_val=partition)
      END IF

      SELECT CASE (partition)
      CASE (skala_gapw_density_partition_hard_minus_soft, &
            skala_gapw_density_partition_hard_only, &
            skala_gapw_density_partition_soft_only, &
            skala_gapw_density_partition_none)
         CONTINUE
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "Unknown GAUXC%NATIVE_GRID_GAPW_DENSITY_PARTITION value.")
      END SELECT

   END FUNCTION native_skala_gapw_density_partition

! **************************************************************************************************
!> \brief Enforce the currently implemented native SKALA GPW input scope.
!> \param xc_section ...
! **************************************************************************************************
   SUBROUTINE ensure_native_skala_grid_scope(xc_section)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section

      CHARACTER(len=default_path_length)                 :: model_key, model_name
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
      INTEGER :: iw, native_grid_atom_chunk_max_rows, native_grid_atom_partition, &
         native_grid_atom_subchunks, native_grid_cuda_device, nspins, phase_handle, &
         selected_cuda_device, xc_deriv_method_id, xc_rho_smooth_id
      LOGICAL :: has_atom_chunk_work, have_atom_coord_grad, lsd, my_just_energy, &
         native_grid_atom_chunk_routing, native_grid_atom_chunks, native_grid_diagnostics, &
         native_grid_use_cuda, needs_atom_force, use_atom_subchunks
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: density_grad, kin_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: grad_grad
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_before
      TYPE(section_vals_type), POINTER                   :: gauxc_section
      TYPE(skala_gpw_feature_type)                       :: features
      TYPE(torch_tensor_type)                            :: atom_coord_grad_t, &
                                                            atomic_grid_weight_grad_t, exc_tensor, &
                                                            grid_coord_grad_t, grid_weight_grad_t
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set

      virial_xc = 0.0_dp
      exc = 0.0_dp
      my_just_energy = .FALSE.
      IF (PRESENT(just_energy)) my_just_energy = just_energy
      needs_atom_force = PRESENT(atom_force)
      IF (needs_atom_force) atom_force = 0.0_dp
      have_atom_coord_grad = .FALSE.

      IF (compute_virial .AND. my_just_energy) THEN
         CALL cp_abort(__LOCATION__, &
                       "Native SKALA GPW stress/virial requires feature gradients.")
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
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_ATOM_CHUNK_MAX_ROWS", &
                                i_val=native_grid_atom_chunk_max_rows)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_ATOM_PARTITION", &
                                i_val=native_grid_atom_partition)
      SELECT CASE (native_grid_atom_partition)
      CASE (1)
         native_grid_atom_partition = skala_gpw_atom_partition_hard
      CASE (2)
         native_grid_atom_partition = skala_gpw_atom_partition_smooth
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "Unknown GAUXC%NATIVE_GRID_ATOM_PARTITION value.")
      END SELECT
      native_grid_atom_chunk_routing = native_grid_atom_chunk_routing .OR. native_grid_atom_chunks
      native_grid_atom_chunks = native_grid_atom_chunks .OR. native_grid_atom_chunk_routing
      IF (native_grid_atom_chunk_max_rows < -1) THEN
         CALL cp_abort(__LOCATION__, &
                       "GAUXC%NATIVE_GRID_ATOM_CHUNK_MAX_ROWS must be -1, zero, or positive.")
      END IF
      IF (needs_atom_force .OR. compute_virial) THEN
         IF (native_grid_atom_partition == skala_gpw_atom_partition_hard) THEN
            native_grid_atom_partition = skala_gpw_atom_partition_smooth
         END IF
         native_grid_atom_chunk_routing = .FALSE.
         native_grid_atom_chunks = .FALSE.
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
                                   requires_coordinate_grad=(needs_atom_force .OR. compute_virial), &
                                   requires_stress_grad=compute_virial, &
                                   use_atom_chunks=native_grid_atom_chunks, &
                                   route_atom_chunks=native_grid_atom_chunk_routing, &
                                   atom_partition=native_grid_atom_partition)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_DIAGNOSTICS", l_val=native_grid_diagnostics)
      IF (native_grid_diagnostics) THEN
         CALL print_native_grid_diagnostics(features, rho_r(1)%pw_grid%para%group%mepos == 0)
      END IF

      IF (features%uses_atom_chunks .AND. native_grid_atom_chunk_max_rows == -1) THEN
         IF (native_grid_use_cuda) THEN
            native_grid_atom_chunk_max_rows = auto_atom_chunk_max_rows(features, &
                                                                       rho_r(1)%pw_grid%para%group)
         ELSE
            native_grid_atom_chunk_max_rows = 0
         END IF
      END IF
      IF (native_grid_diagnostics .AND. features%uses_atom_chunks .AND. &
          rho_r(1)%pw_grid%para%group%mepos == 0) THEN
         iw = cp_logger_get_default_io_unit()
         IF (iw > 0) THEN
            WRITE (UNIT=iw, FMT="(T2,A,1X,I0)") &
               "SKALA_GPW| Native grid atom chunk max rows", native_grid_atom_chunk_max_rows
         END IF
      END IF
      native_grid_atom_subchunks = 1
      IF (features%uses_atom_chunks .AND. native_grid_atom_chunk_max_rows > 0) THEN
         native_grid_atom_subchunks = skala_gpw_atom_subchunk_count(native_grid_atom_chunk_max_rows)
         CALL rho_r(1)%pw_grid%para%group%max(native_grid_atom_subchunks)
      END IF
      use_atom_subchunks = features%uses_atom_chunks .AND. native_grid_atom_subchunks > 1
      has_atom_chunk_work = .NOT. features%uses_atom_chunks .OR. features%chunk_feature_count > 0
      exc = 0.0_dp
      IF (use_atom_subchunks) THEN
         CALL evaluate_atom_subchunks(features, rho_r(1)%pw_grid%para%group, &
                                      native_grid_atom_chunk_max_rows, &
                                      compute_grads=(.NOT. my_just_energy), exc=exc, &
                                      density_grad=density_grad, grad_grad=grad_grad, &
                                      kin_grad=kin_grad, collapse_spin_grads=(nspins == 1))
      ELSE IF (has_atom_chunk_work) THEN
         CALL skala_torch_model_get_exc(cached_model, features%inputs, &
                                        features%grid_weights_t, exc_tensor, exc)
      END IF
      IF (features%uses_atom_chunks) CALL rho_r(1)%pw_grid%para%group%sum(exc)

      IF (.NOT. my_just_energy) THEN
         IF (.NOT. use_atom_subchunks) THEN
            IF (has_atom_chunk_work) THEN
               CALL timeset("skala_gpw_backward", phase_handle)
               CALL torch_tensor_backward_scalar(exc_tensor)
               CALL timestop(phase_handle)

               IF (compute_virial) THEN
                  IF (native_grid_diagnostics) virial_before = virial_xc
                  CALL build_weight_virial(virial_xc, features, exc, grid_weight_grad_t, &
                                           atomic_grid_weight_grad_t, &
                                           rho_r(1)%pw_grid%para%group%mepos == 0, &
                                           native_grid_diagnostics)
                  IF (native_grid_diagnostics) THEN
                     CALL print_virial_delta("weight-residual", virial_xc - virial_before, &
                                             rho_r(1)%pw_grid%para%group%mepos == 0)
                  END IF
               END IF
            END IF

            CALL timeset("skala_gpw_grad_fetch", phase_handle)
            IF (features%uses_atom_chunks) THEN
               CALL fetch_and_gather_atom_chunk_grads(features, rho_r(1)%pw_grid%para%group, &
                                                      density_grad, grad_grad, kin_grad)
            ELSE
               CALL fetch_local_feature_grads(features, density_grad, grad_grad, kin_grad)
            END IF
            CALL timestop(phase_handle)
         END IF
         IF (needs_atom_force) THEN
            CALL add_explicit_coordinate_force(atom_force, features, atom_coord_grad_t, &
                                               rho_r(1)%pw_grid%para%group%mepos == 0)
            IF (features%atom_partition == skala_gpw_atom_partition_smooth) THEN
               CALL add_smooth_partition_force(atom_force, features, particle_set, cell, rho_r, &
                                               grid_weight_grad_t, atomic_grid_weight_grad_t)
            END IF
            have_atom_coord_grad = .TRUE.
         END IF

         CALL timeset("skala_gpw_vxc_unpack", phase_handle)
         IF (compute_virial) THEN
            IF (native_grid_diagnostics) virial_before = virial_xc
            CALL build_virial_from_feature_grads(virial_xc, rho_set, rho_r, grad_grad)
            IF (native_grid_diagnostics) THEN
               CALL print_virial_delta("feature-gradient", virial_xc - virial_before, &
                                       rho_r(1)%pw_grid%para%group%mepos == 0)
               virial_before = virial_xc
            END IF
            IF (.NOT. have_atom_coord_grad) THEN
               CALL torch_tensor_grad(features%coarse_0_atomic_coords_t, atom_coord_grad_t)
               have_atom_coord_grad = .TRUE.
            END IF
            CALL build_static_coordinate_virial(virial_xc, features, atom_coord_grad_t, &
                                                grid_coord_grad_t, &
                                                rho_r(1)%pw_grid%para%group%mepos == 0, &
                                                native_grid_diagnostics)
            IF (native_grid_diagnostics) THEN
               CALL print_virial_delta("static-coordinates", virial_xc - virial_before, &
                                       rho_r(1)%pw_grid%para%group%mepos == 0)
               virial_before = virial_xc
            END IF
            IF (features%atom_partition == skala_gpw_atom_partition_smooth) THEN
               CALL build_smooth_partition_virial(virial_xc, features, particle_set, cell, rho_r, &
                                                  grid_weight_grad_t, atomic_grid_weight_grad_t)
               IF (native_grid_diagnostics) THEN
                  CALL print_virial_delta("smooth-partition", virial_xc - virial_before, &
                                          rho_r(1)%pw_grid%para%group%mepos == 0)
                  virial_before = virial_xc
               END IF
            END IF
         END IF
         CALL build_vxc_from_feature_grads(vxc_rho, vxc_tau, rho_r, pw_pool, &
                                           density_grad, grad_grad, kin_grad, &
                                           xc_deriv_method_id)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_release", phase_handle)
         DEALLOCATE (density_grad, grad_grad, kin_grad)
         IF (have_atom_coord_grad) CALL torch_tensor_release(atom_coord_grad_t)
         CALL timestop(phase_handle)
      END IF

      CALL timeset("skala_gpw_cleanup", phase_handle)
      IF (.NOT. use_atom_subchunks .AND. has_atom_chunk_work) CALL torch_tensor_release(exc_tensor)
      CALL skala_gpw_feature_release(features)
      CALL xc_rho_set_release(rho_set, pw_pool=pw_pool)
      CALL torch_use_cuda(.TRUE.)
      CALL timestop(phase_handle)

   END SUBROUTINE skala_gpw_eval

! **************************************************************************************************
!> \brief Evaluate the native SKALA XC energy density on the CP2K PW grid.
!> \param exc_r ...
!> \param rho_r ...
!> \param rho_g ...
!> \param tau ...
!> \param xc_section ...
!> \param weights ...
!> \param pw_pool ...
!> \param particle_set ...
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE skala_gpw_exc_density(exc_r, rho_r, rho_g, tau, xc_section, weights, pw_pool, &
                                    particle_set, cell)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: exc_r
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: tau
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell

      CHARACTER(len=default_path_length)                 :: model_path
      INTEGER :: feature_pos, i, j, k, local_row, native_grid_atom_partition, &
         native_grid_cuda_device, nspins, row, selected_cuda_device, xc_deriv_method_id, &
         xc_rho_smooth_id
      LOGICAL                                            :: lsd, native_grid_atom_chunk_routing, &
                                                            native_grid_atom_chunks, &
                                                            native_grid_use_cuda
      REAL(KIND=dp)                                      :: local_exc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: exc_density
      TYPE(section_vals_type), POINTER                   :: gauxc_section
      TYPE(skala_gpw_feature_type)                       :: features
      TYPE(torch_tensor_type)                            :: exc_density_t
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set

      CPASSERT(ASSOCIATED(rho_r))
      CPASSERT(ASSOCIATED(rho_g))
      CPASSERT(ASSOCIATED(tau))
      CALL pw_zero(exc_r)

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
      native_grid_atom_chunks = .FALSE.
      native_grid_atom_chunk_routing = .FALSE.
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_ATOM_PARTITION", &
                                i_val=native_grid_atom_partition)
      SELECT CASE (native_grid_atom_partition)
      CASE (1)
         native_grid_atom_partition = skala_gpw_atom_partition_hard
      CASE (2)
         native_grid_atom_partition = skala_gpw_atom_partition_smooth
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "Unknown GAUXC%NATIVE_GRID_ATOM_PARTITION value.")
      END SELECT

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
                                   requires_grad=.FALSE., weights=weights, &
                                   requires_coordinate_grad=.FALSE., &
                                   requires_stress_grad=.FALSE., &
                                   use_atom_chunks=.FALSE., route_atom_chunks=.FALSE., &
                                   atom_partition=native_grid_atom_partition)
      CALL skala_torch_model_get_exc_density(cached_model, features%inputs, exc_density_t)
      NULLIFY (exc_density)
      CALL torch_tensor_data_ptr(exc_density_t, exc_density)

      local_row = 0
      DO k = LBOUND(features%feature_index, 3), UBOUND(features%feature_index, 3)
         DO j = LBOUND(features%feature_index, 2), UBOUND(features%feature_index, 2)
            DO i = LBOUND(features%feature_index, 1), UBOUND(features%feature_index, 1)
               local_row = local_row + 1
               local_exc = 0.0_dp
               DO feature_pos = features%local_feature_offsets(local_row), &
                  features%local_feature_offsets(local_row + 1) - 1
                  row = features%local_feature_rows(feature_pos)
                  local_exc = local_exc + exc_density(row)*features%grid_weights(row)
               END DO
               exc_r%array(i, j, k) = local_exc/rho_r(1)%pw_grid%dvol
            END DO
         END DO
      END DO
      CPASSERT(local_row == features%nflat_local)

      CALL torch_tensor_release(exc_density_t)
      CALL skala_gpw_feature_release(features)
      CALL xc_rho_set_release(rho_set, pw_pool=pw_pool)
      CALL torch_use_cuda(.TRUE.)

   END SUBROUTINE skala_gpw_exc_density

! **************************************************************************************************
!> \brief Evaluate SKALA on a GAPW one-center atomic grid.
!> \param xc_section ...
!> \param grid_atom ...
!> \param group ...
!> \param atom_coord ...
!> \param rho ...
!> \param drho ...
!> \param tau ...
!> \param weights ...
!> \param lsd ...
!> \param nspins ...
!> \param na ...
!> \param nr ...
!> \param exc ...
!> \param vxc ...
!> \param vxg ...
!> \param vtau ...
!> \param energy_only ...
!> \param atom_force ...
!> \param atom_virial ...
! **************************************************************************************************
   SUBROUTINE skala_gapw_atom_vxc_of_r(xc_section, grid_atom, group, atom_coord, &
                                       rho, drho, tau, weights, lsd, nspins, na, nr, &
                                       exc, vxc, vxg, vtau, energy_only, atom_force, atom_virial)
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(grid_atom_type), POINTER                      :: grid_atom

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: atom_coord
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho, tau, vxc, vtau
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: drho, vxg
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: weights
      LOGICAL, INTENT(IN)                                :: lsd
      INTEGER, INTENT(IN)                                :: nspins, na, nr
      REAL(KIND=dp), INTENT(OUT)                         :: exc
      LOGICAL, INTENT(IN), OPTIONAL                      :: energy_only
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT), &
         OPTIONAL                                        :: atom_force
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT), &
         OPTIONAL                                        :: atom_virial

      CHARACTER(len=default_path_length)                 :: model_path
      INTEGER                                            :: ia, idir, ir, native_grid_cuda_device, &
                                                            jdir, nflat, row, selected_cuda_device
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: atomic_grid_sizes
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :)  :: atomic_grid_size_bound_shape
      LOGICAL                                            :: need_coord_grad, my_energy_only, native_grid_use_cuda
      REAL(KIND=dp)                                      :: tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: atomic_grid_weights, grid_weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: coarse_0_atomic_coords, density, &
                                                            grid_coords, kin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: grad
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: atom_coord_grad, density_grad, &
                                                            grid_coord_grad, kin_grad
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grad_grad
      TYPE(section_vals_type), POINTER                   :: gauxc_section
      TYPE(torch_dict_type)                              :: inputs
      TYPE(torch_tensor_type)                            :: atomic_grid_size_bound_shape_t, &
                                                            atomic_grid_sizes_t, &
                                                            atomic_grid_weights_t, &
                                                            atom_coord_grad_t, &
                                                            coarse_0_atomic_coords_t, density_t, &
                                                            density_grad_t, exc_tensor, grad_t, &
                                                            grad_grad_t, grid_coord_grad_t, &
                                                            grid_coords_t, grid_weights_t, kin_t, &
                                                            kin_grad_t

      CPASSERT(ASSOCIATED(xc_section))
      CPASSERT(ASSOCIATED(grid_atom))
      CPASSERT(ASSOCIATED(rho))
      CPASSERT(ASSOCIATED(drho))
      CPASSERT(ASSOCIATED(tau))

      my_energy_only = .FALSE.
      IF (PRESENT(energy_only)) my_energy_only = energy_only
      need_coord_grad = PRESENT(atom_force) .OR. PRESENT(atom_virial)
      exc = 0.0_dp
      IF (PRESENT(atom_force)) atom_force = 0.0_dp
      IF (PRESENT(atom_virial)) atom_virial = 0.0_dp
      IF (.NOT. my_energy_only) THEN
         vxc = 0.0_dp
         vxg = 0.0_dp
         vtau = 0.0_dp
      END IF

      CALL get_skala_model_path(xc_section, model_path)
      gauxc_section => get_gauxc_section(xc_section)
      CPASSERT(ASSOCIATED(gauxc_section))
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_USE_CUDA", l_val=native_grid_use_cuda)
      CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_CUDA_DEVICE", &
                                i_val=native_grid_cuda_device)
      CALL torch_use_cuda(native_grid_use_cuda)
      selected_cuda_device = configure_native_grid_cuda( &
                             native_grid_use_cuda, native_grid_cuda_device, group)
      CALL ensure_model_loaded(model_path, selected_cuda_device)

      nflat = na*nr
      ALLOCATE (density(nflat, 2), grad(nflat, 3, 2), kin(nflat, 2), &
                grid_coords(3, nflat), grid_weights(nflat), &
                atomic_grid_weights(nflat), atomic_grid_sizes(1), &
                coarse_0_atomic_coords(3, 1), atomic_grid_size_bound_shape(0, nflat))
      density = 0.0_dp
      grad = 0.0_dp
      kin = 0.0_dp
      grid_coords = 0.0_dp
      grid_weights = 0.0_dp
      atomic_grid_weights = 0.0_dp
      atomic_grid_sizes(1) = INT(nflat, KIND=int_8)
      atomic_grid_size_bound_shape = 0_int_8
      coarse_0_atomic_coords(:, 1) = atom_coord

      row = 0
      DO ir = 1, nr
         DO ia = 1, na
            row = row + 1
            grid_coords(1, row) = atom_coord(1) + grid_atom%rad(ir)* &
                                  grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia)
            grid_coords(2, row) = atom_coord(2) + grid_atom%rad(ir)* &
                                  grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia)
            grid_coords(3, row) = atom_coord(3) + grid_atom%rad(ir)*grid_atom%cos_pol(ia)
            grid_weights(row) = weights(ia, ir)
            atomic_grid_weights(row) = weights(ia, ir)
            IF (nspins == 1) THEN
               density(row, :) = 0.5_dp*rho(ia, ir, 1)
               DO idir = 1, 3
                  grad(row, idir, :) = 0.5_dp*drho(idir, ia, ir, 1)
               END DO
               kin(row, :) = 0.5_dp*tau(ia, ir, 1)
            ELSE
               density(row, :) = rho(ia, ir, 1:2)
               DO idir = 1, 3
                  grad(row, idir, :) = drho(idir, ia, ir, 1:2)
               END DO
               kin(row, :) = tau(ia, ir, 1:2)
            END IF
         END DO
      END DO

      CALL torch_tensor_from_array(grid_coords_t, grid_coords)
      CALL torch_tensor_to_device_leaf(grid_coords_t, need_coord_grad)
      CALL torch_tensor_from_array(grid_weights_t, grid_weights)
      CALL torch_tensor_to_device_leaf(grid_weights_t, .FALSE.)
      CALL torch_tensor_from_array(atomic_grid_weights_t, atomic_grid_weights)
      CALL torch_tensor_to_device_leaf(atomic_grid_weights_t, .FALSE.)
      CALL torch_tensor_from_array(atomic_grid_sizes_t, atomic_grid_sizes)
      CALL torch_tensor_to_device_leaf(atomic_grid_sizes_t, .FALSE.)
      CALL torch_tensor_from_array(atomic_grid_size_bound_shape_t, &
                                   atomic_grid_size_bound_shape)
      CALL torch_tensor_to_device_leaf(atomic_grid_size_bound_shape_t, .FALSE.)
      CALL torch_tensor_from_array(coarse_0_atomic_coords_t, coarse_0_atomic_coords)
      CALL torch_tensor_to_device_leaf(coarse_0_atomic_coords_t, need_coord_grad)
      CALL torch_tensor_from_array(density_t, density)
      CALL torch_tensor_to_device_leaf(density_t,.NOT. my_energy_only)
      CALL torch_tensor_from_array(grad_t, grad)
      CALL torch_tensor_to_device_leaf(grad_t,.NOT. my_energy_only)
      CALL torch_tensor_from_array(kin_t, kin)
      CALL torch_tensor_to_device_leaf(kin_t,.NOT. my_energy_only)

      CALL torch_dict_create(inputs)
      CALL torch_dict_insert(inputs, "grid_coords", grid_coords_t)
      CALL torch_dict_insert(inputs, "grid_weights", grid_weights_t)
      CALL torch_dict_insert(inputs, "atomic_grid_weights", atomic_grid_weights_t)
      CALL torch_dict_insert(inputs, "atomic_grid_sizes", atomic_grid_sizes_t)
      CALL torch_dict_insert(inputs, "atomic_grid_size_bound_shape", &
                             atomic_grid_size_bound_shape_t)
      CALL torch_dict_insert(inputs, "density", density_t)
      CALL torch_dict_insert(inputs, "grad", grad_t)
      CALL torch_dict_insert(inputs, "kin", kin_t)
      CALL torch_dict_insert(inputs, "coarse_0_atomic_coords", coarse_0_atomic_coords_t)

      CALL skala_torch_model_get_exc(cached_model, inputs, grid_weights_t, exc_tensor, exc)

      IF (.NOT. my_energy_only) THEN
         NULLIFY (atom_coord_grad, density_grad, grad_grad, grid_coord_grad, kin_grad)
         CALL torch_tensor_backward_scalar(exc_tensor)
         IF (need_coord_grad) THEN
            CALL torch_tensor_grad(grid_coords_t, grid_coord_grad_t)
            CALL torch_tensor_grad(coarse_0_atomic_coords_t, atom_coord_grad_t)
            CALL torch_tensor_data_ptr(grid_coord_grad_t, grid_coord_grad)
            CALL torch_tensor_data_ptr(atom_coord_grad_t, atom_coord_grad)
            IF (PRESENT(atom_force)) THEN
               atom_force(:) = atom_coord_grad(:, 1)
               DO row = 1, nflat
                  atom_force(:) = atom_force(:) + grid_coord_grad(:, row)
               END DO
            END IF
            IF (PRESENT(atom_virial)) THEN
               DO row = 1, nflat
                  DO idir = 1, 3
                     DO jdir = 1, 3
                        tmp = grid_coord_grad(idir, row)*coarse_0_atomic_coords(jdir, 1)
                        atom_virial(idir, jdir) = atom_virial(idir, jdir) + tmp
                     END DO
                  END DO
               END DO
               DO idir = 1, 3
                  DO jdir = 1, 3
                     tmp = atom_coord_grad(idir, 1)*coarse_0_atomic_coords(jdir, 1)
                     atom_virial(idir, jdir) = atom_virial(idir, jdir) + tmp
                  END DO
               END DO
            END IF
         END IF
         CALL torch_tensor_grad(density_t, density_grad_t)
         CALL torch_tensor_grad(grad_t, grad_grad_t)
         CALL torch_tensor_grad(kin_t, kin_grad_t)
         CALL torch_tensor_data_ptr(density_grad_t, density_grad)
         CALL torch_tensor_data_ptr(grad_grad_t, grad_grad)
         CALL torch_tensor_data_ptr(kin_grad_t, kin_grad)

         row = 0
         DO ir = 1, nr
            DO ia = 1, na
               row = row + 1
               IF (lsd) THEN
                  vxc(ia, ir, 1:2) = density_grad(row, 1:2)
                  DO idir = 1, 3
                     vxg(idir, ia, ir, 1:2) = grad_grad(row, idir, 1:2)
                  END DO
                  vtau(ia, ir, 1:2) = kin_grad(row, 1:2)
               ELSE
                  vxc(ia, ir, 1) = 0.5_dp*(density_grad(row, 1) + density_grad(row, 2))
                  DO idir = 1, 3
                     vxg(idir, ia, ir, 1) = &
                        0.5_dp*(grad_grad(row, idir, 1) + grad_grad(row, idir, 2))
                  END DO
                  vtau(ia, ir, 1) = 0.5_dp*(kin_grad(row, 1) + kin_grad(row, 2))
               END IF
            END DO
         END DO

         CALL torch_tensor_release(density_grad_t)
         CALL torch_tensor_release(grad_grad_t)
         CALL torch_tensor_release(kin_grad_t)
         IF (need_coord_grad) THEN
            CALL torch_tensor_release(grid_coord_grad_t)
            CALL torch_tensor_release(atom_coord_grad_t)
         END IF
      END IF

      CALL torch_tensor_release(exc_tensor)
      CALL torch_tensor_release(density_t)
      CALL torch_tensor_release(grad_t)
      CALL torch_tensor_release(kin_t)
      CALL torch_tensor_release(grid_coords_t)
      CALL torch_tensor_release(grid_weights_t)
      CALL torch_tensor_release(atomic_grid_weights_t)
      CALL torch_tensor_release(atomic_grid_sizes_t)
      CALL torch_tensor_release(atomic_grid_size_bound_shape_t)
      CALL torch_tensor_release(coarse_0_atomic_coords_t)
      CALL torch_dict_release(inputs)
      DEALLOCATE (atomic_grid_size_bound_shape, atomic_grid_sizes, atomic_grid_weights, &
                  coarse_0_atomic_coords, density, grad, grid_coords, grid_weights, kin)
      CALL torch_use_cuda(.TRUE.)

   END SUBROUTINE skala_gapw_atom_vxc_of_r

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
!> \brief Add the force from SMOOTH native-grid atom partition weights.
!> \param atom_force ...
!> \param features ...
!> \param particle_set ...
!> \param cell ...
!> \param rho_r ...
!> \param grid_weight_grad_t ...
!> \param atomic_grid_weight_grad_t ...
! **************************************************************************************************
   SUBROUTINE add_smooth_partition_force(atom_force, features, particle_set, cell, rho_r, &
                                         grid_weight_grad_t, atomic_grid_weight_grad_t)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: atom_force
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(torch_tensor_type), INTENT(INOUT)             :: grid_weight_grad_t, &
                                                            atomic_grid_weight_grad_t

      INTEGER                                            :: feature_begin, feature_end, feature_pos, &
                                                            i, iatom, j, jatom, k, local_row, &
                                                            natom, row
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: included
      REAL(KIND=dp)                                      :: base_weight, weight_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_coords_pbc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dweights_datom, dweights_dstrain
      REAL(KIND=dp), DIMENSION(3)                        :: grid_point
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atomic_grid_weight_grad, grid_weight_grad

      NULLIFY (atomic_grid_weight_grad, grid_weight_grad)
      CALL torch_tensor_grad(features%grid_weights_t, grid_weight_grad_t)
      CALL torch_tensor_grad(features%atomic_grid_weights_t, atomic_grid_weight_grad_t)
      CALL torch_tensor_data_ptr(grid_weight_grad_t, grid_weight_grad)
      CALL torch_tensor_data_ptr(atomic_grid_weight_grad_t, atomic_grid_weight_grad)

      natom = SIZE(particle_set)
      CPASSERT(SIZE(atom_force, 1) == 3)
      CPASSERT(SIZE(atom_force, 2) == natom)
      ALLOCATE (atom_coords_pbc(3, natom), included(natom), weights(natom), &
                dweights_datom(3, natom, natom), dweights_dstrain(3, 3, natom))
      DO iatom = 1, natom
         atom_coords_pbc(:, iatom) = pbc(particle_set(iatom)%r, cell, positive_range=.TRUE.)
      END DO

      bo = rho_r(1)%pw_grid%bounds_local
      local_row = 0
      DO k = bo(1, 3), bo(2, 3)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 1), bo(2, 1)
               local_row = local_row + 1
               grid_point = native_grid_coordinate(rho_r(1)%pw_grid, [i, j, k])
               CALL skala_gpw_smooth_partition_derivatives(grid_point, atom_coords_pbc, cell, &
                                                           weights, included, dweights_datom, &
                                                           dweights_dstrain)
               feature_begin = features%local_feature_offsets(local_row)
               feature_end = features%local_feature_offsets(local_row + 1) - 1
               CPASSERT(feature_end - feature_begin + 1 == COUNT(included))
               base_weight = 0.0_dp
               DO feature_pos = feature_begin, feature_end
                  row = features%local_feature_rows(feature_pos)
                  base_weight = base_weight + features%grid_weights(row)
               END DO
               feature_pos = feature_begin
               DO iatom = 1, natom
                  IF (.NOT. included(iatom)) CYCLE
                  row = features%local_feature_rows(feature_pos)
                  weight_grad = grid_weight_grad(row)
                  DO jatom = 1, natom
                     atom_force(:, jatom) = atom_force(:, jatom) + &
                                            weight_grad*base_weight* &
                                            dweights_datom(:, jatom, iatom)
                  END DO
                  feature_pos = feature_pos + 1
               END DO
               CPASSERT(feature_pos == feature_end + 1)
            END DO
         END DO
      END DO
      CPASSERT(local_row == features%nflat_local)

      DEALLOCATE (atom_coords_pbc, dweights_datom, dweights_dstrain, included, weights)
      CALL torch_tensor_release(grid_weight_grad_t)
      CALL torch_tensor_release(atomic_grid_weight_grad_t)

   END SUBROUTINE add_smooth_partition_force

! **************************************************************************************************
!> \brief Add the virial from SMOOTH native-grid atom partition weights.
!> \param virial_xc ...
!> \param features ...
!> \param particle_set ...
!> \param cell ...
!> \param rho_r ...
!> \param grid_weight_grad_t ...
!> \param atomic_grid_weight_grad_t ...
! **************************************************************************************************
   SUBROUTINE build_smooth_partition_virial(virial_xc, features, particle_set, cell, rho_r, &
                                            grid_weight_grad_t, atomic_grid_weight_grad_t)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: virial_xc
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(torch_tensor_type), INTENT(INOUT)             :: grid_weight_grad_t, &
                                                            atomic_grid_weight_grad_t

      INTEGER                                            :: feature_begin, feature_end, feature_pos, &
                                                            i, iatom, idir, j, jdir, k, local_row, &
                                                            natom, row
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: included
      REAL(KIND=dp)                                      :: base_weight, tmp, weight_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_coords_pbc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dweights_datom, dweights_dstrain
      REAL(KIND=dp), DIMENSION(3)                        :: grid_point
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atomic_grid_weight_grad, grid_weight_grad

      NULLIFY (atomic_grid_weight_grad, grid_weight_grad)
      CALL torch_tensor_grad(features%grid_weights_t, grid_weight_grad_t)
      CALL torch_tensor_grad(features%atomic_grid_weights_t, atomic_grid_weight_grad_t)
      CALL torch_tensor_data_ptr(grid_weight_grad_t, grid_weight_grad)
      CALL torch_tensor_data_ptr(atomic_grid_weight_grad_t, atomic_grid_weight_grad)

      natom = SIZE(particle_set)
      ALLOCATE (atom_coords_pbc(3, natom), included(natom), weights(natom), &
                dweights_datom(3, natom, natom), dweights_dstrain(3, 3, natom))
      DO iatom = 1, natom
         atom_coords_pbc(:, iatom) = pbc(particle_set(iatom)%r, cell, positive_range=.TRUE.)
      END DO

      bo = rho_r(1)%pw_grid%bounds_local
      local_row = 0
      DO k = bo(1, 3), bo(2, 3)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 1), bo(2, 1)
               local_row = local_row + 1
               grid_point = native_grid_coordinate(rho_r(1)%pw_grid, [i, j, k])
               CALL skala_gpw_smooth_partition_derivatives(grid_point, atom_coords_pbc, cell, &
                                                           weights, included, dweights_datom, &
                                                           dweights_dstrain)
               feature_begin = features%local_feature_offsets(local_row)
               feature_end = features%local_feature_offsets(local_row + 1) - 1
               CPASSERT(feature_end - feature_begin + 1 == COUNT(included))
               base_weight = 0.0_dp
               DO feature_pos = feature_begin, feature_end
                  row = features%local_feature_rows(feature_pos)
                  base_weight = base_weight + features%grid_weights(row)
               END DO
               feature_pos = feature_begin
               DO iatom = 1, natom
                  IF (.NOT. included(iatom)) CYCLE
                  row = features%local_feature_rows(feature_pos)
                  weight_grad = grid_weight_grad(row)
                  DO idir = 1, 3
                     DO jdir = 1, idir
                        tmp = weight_grad*base_weight*dweights_dstrain(idir, jdir, iatom)
                        virial_xc(jdir, idir) = virial_xc(jdir, idir) + tmp
                        IF (idir /= jdir) virial_xc(idir, jdir) = virial_xc(idir, jdir) + tmp
                     END DO
                  END DO
                  feature_pos = feature_pos + 1
               END DO
               CPASSERT(feature_pos == feature_end + 1)
            END DO
         END DO
      END DO
      CPASSERT(local_row == features%nflat_local)

      DEALLOCATE (atom_coords_pbc, dweights_datom, dweights_dstrain, included, weights)
      CALL torch_tensor_release(grid_weight_grad_t)
      CALL torch_tensor_release(atomic_grid_weight_grad_t)

   END SUBROUTINE build_smooth_partition_virial

! **************************************************************************************************
!> \brief Return the Cartesian coordinate of a regular GPW grid point.
!> \param pw_grid ...
!> \param index ...
!> \return ...
! **************************************************************************************************
   FUNCTION native_grid_coordinate(pw_grid, index) RESULT(coord)
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      INTEGER, DIMENSION(3), INTENT(IN)                  :: index
      REAL(KIND=dp), DIMENSION(3)                        :: coord

      INTEGER, DIMENSION(3)                              :: relative_index

      relative_index = index - pw_grid%bounds(1, :)
      coord = REAL(relative_index(1), KIND=dp)*pw_grid%dh(:, 1) + &
              REAL(relative_index(2), KIND=dp)*pw_grid%dh(:, 2) + &
              REAL(relative_index(3), KIND=dp)*pw_grid%dh(:, 3)

   END FUNCTION native_grid_coordinate

! **************************************************************************************************
!> \brief Evaluate a rank-local atom chunk as multiple atom-contiguous Torch subchunks.
!> \param features ...
!> \param group ...
!> \param max_rows ...
!> \param compute_grads ...
!> \param exc ...
!> \param density_grad ...
!> \param grad_grad ...
!> \param kin_grad ...
!> \param collapse_spin_grads ...
! **************************************************************************************************
   SUBROUTINE evaluate_atom_subchunks(features, group, max_rows, compute_grads, exc, &
                                      density_grad, grad_grad, kin_grad, collapse_spin_grads)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      INTEGER, INTENT(IN)                                :: max_rows
      LOGICAL, INTENT(IN)                                :: compute_grads, collapse_spin_grads
      REAL(KIND=dp), INTENT(OUT)                         :: exc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: density_grad, kin_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: grad_grad

      INTEGER                                            :: base, isubchunk, local_row, nflat_local, &
                                                            nroute_grad_per_point, nroute_points, &
                                                            nsubchunks, phase_handle, point_pos, &
                                                            subphase_handle
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: route_grad_return_recv_counts, &
                                                            route_grad_return_recv_displs, &
                                                            route_grad_return_send_counts, &
                                                            route_grad_return_send_displs
      REAL(KIND=dp)                                      :: subchunk_exc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: recv_grad_buffer, send_grad_buffer
      TYPE(skala_gpw_feature_type)                       :: subchunk
      TYPE(torch_tensor_type)                            :: subchunk_exc_tensor

      CPASSERT(features%uses_atom_chunks)
      CPASSERT(max_rows > 0)
      nflat_local = features%nflat_local
      nsubchunks = skala_gpw_atom_subchunk_count(max_rows)

      exc = 0.0_dp
      IF (compute_grads) THEN
         CPASSERT(features%uses_atom_chunk_routing)
         CPASSERT(SUM(features%route_point_recv_counts) == features%chunk_feature_count)
         nroute_points = SIZE(features%route_send_local_rows)
         CPASSERT(SUM(features%route_point_send_counts) == nroute_points)
         nroute_grad_per_point = ngrad_per_point
         IF (collapse_spin_grads) nroute_grad_per_point = ncollapsed_grad_per_point
         ALLOCATE (send_grad_buffer(MAX(1, nroute_grad_per_point*features%chunk_feature_count)), &
                   recv_grad_buffer(MAX(1, nroute_grad_per_point*nroute_points)), &
                   route_grad_return_send_counts(SIZE(features%route_point_recv_counts)), &
                   route_grad_return_send_displs(SIZE(features%route_point_recv_displs)), &
                   route_grad_return_recv_counts(SIZE(features%route_point_send_counts)), &
                   route_grad_return_recv_displs(SIZE(features%route_point_send_displs)))
         route_grad_return_send_counts(:) = &
            nroute_grad_per_point*features%route_point_recv_counts
         route_grad_return_send_displs(:) = &
            nroute_grad_per_point*features%route_point_recv_displs
         route_grad_return_recv_counts(:) = &
            nroute_grad_per_point*features%route_point_send_counts
         route_grad_return_recv_displs(:) = &
            nroute_grad_per_point*features%route_point_send_displs
      END IF

      CALL timeset("skala_gpw_atom_subchunks", phase_handle)
      DO isubchunk = 1, nsubchunks
         CALL timeset("skala_gpw_atom_subchunk_build", subphase_handle)
         CALL skala_gpw_feature_build_atom_subchunk(features, subchunk, isubchunk, &
                                                    max_rows, compute_grads)
         CALL timestop(subphase_handle)
         CALL timeset("skala_gpw_atom_subchunk_forward", subphase_handle)
         CALL skala_torch_model_get_exc(cached_model, subchunk%inputs, &
                                        subchunk%grid_weights_t, subchunk_exc_tensor, &
                                        subchunk_exc)
         CALL timestop(subphase_handle)
         exc = exc + subchunk_exc
         IF (compute_grads) THEN
            CALL timeset("skala_gpw_atom_subchunk_backward", subphase_handle)
            CALL torch_tensor_backward_scalar(subchunk_exc_tensor)
            CALL timestop(subphase_handle)
         END IF
         CALL timeset("skala_gpw_atom_subchunk_release", subphase_handle)
         CALL torch_tensor_release(subchunk_exc_tensor)
         CALL skala_gpw_feature_release(subchunk)
         CALL timestop(subphase_handle)
      END DO
      IF (compute_grads .AND. features%chunk_feature_count > 0) THEN
         CALL timeset("skala_gpw_atom_subchunk_grad_pack", subphase_handle)
         CALL pack_atom_chunk_grads(features, send_grad_buffer, .TRUE., collapse_spin_grads)
         CALL timestop(subphase_handle)
      END IF
      CALL timestop(phase_handle)

      IF (compute_grads) THEN
         CALL timeset("skala_gpw_grad_route_comm", phase_handle)
         CALL group%alltoall(send_grad_buffer, route_grad_return_send_counts, &
                             route_grad_return_send_displs, recv_grad_buffer, &
                             route_grad_return_recv_counts, route_grad_return_recv_displs)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_route_scatter", phase_handle)
         ALLOCATE (density_grad(nflat_local, 2), grad_grad(nflat_local, 3, 2), &
                   kin_grad(nflat_local, 2))
         density_grad = 0.0_dp
         grad_grad = 0.0_dp
         kin_grad = 0.0_dp
         DO point_pos = 1, nroute_points
            local_row = features%route_send_local_rows(point_pos)
            CPASSERT(local_row >= 1 .AND. local_row <= nflat_local)
            base = nroute_grad_per_point*(point_pos - 1)
            IF (collapse_spin_grads) THEN
               density_grad(local_row, :) = density_grad(local_row, :) + &
                                            recv_grad_buffer(base + 1)
               grad_grad(local_row, 1, :) = grad_grad(local_row, 1, :) + &
                                            recv_grad_buffer(base + 2)
               grad_grad(local_row, 2, :) = grad_grad(local_row, 2, :) + &
                                            recv_grad_buffer(base + 3)
               grad_grad(local_row, 3, :) = grad_grad(local_row, 3, :) + &
                                            recv_grad_buffer(base + 4)
               kin_grad(local_row, :) = kin_grad(local_row, :) + recv_grad_buffer(base + 5)
            ELSE
               density_grad(local_row, :) = density_grad(local_row, :) + &
                                            recv_grad_buffer(base + 1:base + 2)
               grad_grad(local_row, 1, 1) = grad_grad(local_row, 1, 1) + &
                                            recv_grad_buffer(base + 3)
               grad_grad(local_row, 2, 1) = grad_grad(local_row, 2, 1) + &
                                            recv_grad_buffer(base + 4)
               grad_grad(local_row, 3, 1) = grad_grad(local_row, 3, 1) + &
                                            recv_grad_buffer(base + 5)
               grad_grad(local_row, 1, 2) = grad_grad(local_row, 1, 2) + &
                                            recv_grad_buffer(base + 6)
               grad_grad(local_row, 2, 2) = grad_grad(local_row, 2, 2) + &
                                            recv_grad_buffer(base + 7)
               grad_grad(local_row, 3, 2) = grad_grad(local_row, 3, 2) + &
                                            recv_grad_buffer(base + 8)
               kin_grad(local_row, :) = kin_grad(local_row, :) + &
                                        recv_grad_buffer(base + 9:base + 10)
            END IF
         END DO
         CALL timestop(phase_handle)

         DEALLOCATE (recv_grad_buffer, route_grad_return_recv_counts, &
                     route_grad_return_recv_displs, route_grad_return_send_counts, &
                     route_grad_return_send_displs, send_grad_buffer)
      END IF

   END SUBROUTINE evaluate_atom_subchunks

! **************************************************************************************************
!> \brief Select an automatic CUDA atom-subchunk row cap.
!> \param features ...
!> \param group ...
!> \return ...
! **************************************************************************************************
   FUNCTION auto_atom_chunk_max_rows(features, group) RESULT(max_rows)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      INTEGER                                            :: max_rows

      INTEGER                                            :: local_rows_max, target_rows

      local_rows_max = features%chunk_feature_count
      CALL group%max(local_rows_max)
      IF (local_rows_max <= 0) THEN
         max_rows = 0
         RETURN
      END IF

      IF (group%num_pe > 1) THEN
         target_rows = CEILING(REAL(local_rows_max, KIND=dp)/2.0_dp)
         max_rows = atom_chunk_auto_row_quantum* &
                    ((target_rows + atom_chunk_auto_row_quantum - 1)/atom_chunk_auto_row_quantum)
      ELSE
         target_rows = NINT(REAL(local_rows_max, KIND=dp)/4.0_dp)
         max_rows = atom_chunk_auto_row_quantum* &
                    MAX(1, NINT(REAL(target_rows, KIND=dp)/ &
                                REAL(atom_chunk_auto_row_quantum, KIND=dp)))
      END IF
      max_rows = MAX(atom_chunk_auto_min_rows, MIN(atom_chunk_auto_max_rows, max_rows))

   END FUNCTION auto_atom_chunk_max_rows

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

      INTEGER                                            :: feature_pos, i, j, k, local_row, row
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
      density_grad = 0.0_dp
      grad_grad = 0.0_dp
      kin_grad = 0.0_dp
      local_row = 0
      DO k = LBOUND(features%feature_index, 3), UBOUND(features%feature_index, 3)
         DO j = LBOUND(features%feature_index, 2), UBOUND(features%feature_index, 2)
            DO i = LBOUND(features%feature_index, 1), UBOUND(features%feature_index, 1)
               local_row = local_row + 1
               DO feature_pos = features%local_feature_offsets(local_row), &
                  features%local_feature_offsets(local_row + 1) - 1
                  row = features%local_feature_rows(feature_pos)
                  CPASSERT(row >= 1 .AND. row <= features%nflat)
                  density_grad(local_row, :) = density_grad(local_row, :) + &
                                               density_grad_all(row, :)
                  grad_grad(local_row, :, :) = grad_grad(local_row, :, :) + &
                                               grad_grad_all(row, :, :)
                  kin_grad(local_row, :) = kin_grad(local_row, :) + kin_grad_all(row, :)
               END DO
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
!> \param collapse_spin_grads ...
! **************************************************************************************************
   SUBROUTINE pack_atom_chunk_grads(features, TARGET, route_to_return_positions, &
                                    collapse_spin_grads)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: target
      LOGICAL, INTENT(IN)                                :: route_to_return_positions
      LOGICAL, INTENT(IN), OPTIONAL                      :: collapse_spin_grads

      INTEGER                                            :: base, irow, ngrad_buffer_per_point, &
                                                            point_pos, target_points
      LOGICAL                                            :: my_collapse_spin_grads
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: chunk_density_grad, chunk_kin_grad
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: chunk_grad_grad
      TYPE(torch_tensor_type)                            :: density_grad_t, grad_grad_t, kin_grad_t

      my_collapse_spin_grads = .FALSE.
      IF (PRESENT(collapse_spin_grads)) my_collapse_spin_grads = collapse_spin_grads
      ngrad_buffer_per_point = ngrad_per_point
      IF (my_collapse_spin_grads) ngrad_buffer_per_point = ncollapsed_grad_per_point

      NULLIFY (chunk_density_grad, chunk_grad_grad, chunk_kin_grad)
      CALL get_feature_grad_views(features, density_grad_t, grad_grad_t, kin_grad_t, &
                                  chunk_density_grad, chunk_grad_grad, chunk_kin_grad)
      CPASSERT(MOD(SIZE(TARGET), ngrad_buffer_per_point) == 0)
      target_points = SIZE(TARGET)/ngrad_buffer_per_point
      CPASSERT(target_points >= features%chunk_feature_count)
      CPASSERT(SIZE(chunk_density_grad, 1) == features%chunk_feature_count)
      CPASSERT(SIZE(chunk_grad_grad, 1) == features%chunk_feature_count)
      CPASSERT(SIZE(chunk_grad_grad, 2) == 3)
      CPASSERT(SIZE(chunk_kin_grad, 1) == features%chunk_feature_count)
      IF (features%uses_collapsed_rks_dynamic) THEN
         CPASSERT(my_collapse_spin_grads)
         CPASSERT(SIZE(chunk_density_grad, 2) == 1)
         CPASSERT(SIZE(chunk_grad_grad, 3) == 1)
         CPASSERT(SIZE(chunk_kin_grad, 2) == 1)
      ELSE
         CPASSERT(SIZE(chunk_density_grad, 2) == 2)
         CPASSERT(SIZE(chunk_grad_grad, 3) == 2)
         CPASSERT(SIZE(chunk_kin_grad, 2) == 2)
      END IF

      DO irow = 1, features%chunk_feature_count
         IF (route_to_return_positions) THEN
            point_pos = features%chunk_return_positions(irow)
            CPASSERT(point_pos >= 1 .AND. point_pos <= target_points)
         ELSE
            point_pos = irow
         END IF
         base = ngrad_buffer_per_point*(point_pos - 1)
         IF (my_collapse_spin_grads) THEN
            IF (features%uses_collapsed_rks_dynamic) THEN
               TARGET(base + 1) = 0.5_dp*chunk_density_grad(irow, 1)
               TARGET(base + 2) = 0.5_dp*chunk_grad_grad(irow, 1, 1)
               TARGET(base + 3) = 0.5_dp*chunk_grad_grad(irow, 2, 1)
               TARGET(base + 4) = 0.5_dp*chunk_grad_grad(irow, 3, 1)
               TARGET(base + 5) = 0.5_dp*chunk_kin_grad(irow, 1)
            ELSE
               TARGET(base + 1) = 0.5_dp*(chunk_density_grad(irow, 1) + &
                                          chunk_density_grad(irow, 2))
               TARGET(base + 2) = 0.5_dp*(chunk_grad_grad(irow, 1, 1) + &
                                          chunk_grad_grad(irow, 1, 2))
               TARGET(base + 3) = 0.5_dp*(chunk_grad_grad(irow, 2, 1) + &
                                          chunk_grad_grad(irow, 2, 2))
               TARGET(base + 4) = 0.5_dp*(chunk_grad_grad(irow, 3, 1) + &
                                          chunk_grad_grad(irow, 3, 2))
               TARGET(base + 5) = 0.5_dp*(chunk_kin_grad(irow, 1) + chunk_kin_grad(irow, 2))
            END IF
         ELSE
            TARGET(base + 1:base + 2) = chunk_density_grad(irow, :)
            TARGET(base + 3) = chunk_grad_grad(irow, 1, 1)
            TARGET(base + 4) = chunk_grad_grad(irow, 2, 1)
            TARGET(base + 5) = chunk_grad_grad(irow, 3, 1)
            TARGET(base + 6) = chunk_grad_grad(irow, 1, 2)
            TARGET(base + 7) = chunk_grad_grad(irow, 2, 2)
            TARGET(base + 8) = chunk_grad_grad(irow, 3, 2)
            TARGET(base + 9:base + 10) = chunk_kin_grad(irow, :)
         END IF
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

      INTEGER                                            :: base, feature_pos, i, j, k, local_row, &
                                                            nflat_local, nroute_grad_per_point, &
                                                            nroute_points, phase_handle, point_pos, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: route_grad_return_recv_counts, &
                                                            route_grad_return_recv_displs, &
                                                            route_grad_return_send_counts, &
                                                            route_grad_return_send_displs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: chunk_grad_buffer, global_grad_buffer, &
                                                            recv_grad_buffer, send_grad_buffer

      CPASSERT(features%uses_atom_chunks)

      nflat_local = features%nflat_local
      IF (features%uses_atom_chunk_routing) THEN
         CPASSERT(SUM(features%route_point_recv_counts) == features%chunk_feature_count)
         nroute_points = SIZE(features%route_send_local_rows)
         CPASSERT(SUM(features%route_point_send_counts) == nroute_points)

         nroute_grad_per_point = ngrad_per_point
         IF (features%uses_collapsed_rks_dynamic) THEN
            nroute_grad_per_point = ncollapsed_grad_per_point
         END IF
         ALLOCATE (send_grad_buffer(MAX(1, nroute_grad_per_point*features%chunk_feature_count)), &
                   recv_grad_buffer(MAX(1, nroute_grad_per_point*nroute_points)), &
                   route_grad_return_send_counts(SIZE(features%route_point_recv_counts)), &
                   route_grad_return_send_displs(SIZE(features%route_point_recv_displs)), &
                   route_grad_return_recv_counts(SIZE(features%route_point_send_counts)), &
                   route_grad_return_recv_displs(SIZE(features%route_point_send_displs)))
         route_grad_return_send_counts(:) = &
            nroute_grad_per_point*features%route_point_recv_counts
         route_grad_return_send_displs(:) = &
            nroute_grad_per_point*features%route_point_recv_displs
         route_grad_return_recv_counts(:) = &
            nroute_grad_per_point*features%route_point_send_counts
         route_grad_return_recv_displs(:) = &
            nroute_grad_per_point*features%route_point_send_displs

         IF (features%chunk_feature_count > 0) THEN
            CALL timeset("skala_gpw_grad_torch_pack", phase_handle)
            CALL pack_atom_chunk_grads(features, send_grad_buffer, .TRUE., &
                                       features%uses_collapsed_rks_dynamic)
            CALL timestop(phase_handle)
         END IF

         CALL timeset("skala_gpw_grad_route_comm", phase_handle)
         CALL group%alltoall(send_grad_buffer, route_grad_return_send_counts, &
                             route_grad_return_send_displs, recv_grad_buffer, &
                             route_grad_return_recv_counts, route_grad_return_recv_displs)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_route_scatter", phase_handle)
         ALLOCATE (density_grad(nflat_local, 2), grad_grad(nflat_local, 3, 2), &
                   kin_grad(nflat_local, 2))
         density_grad = 0.0_dp
         grad_grad = 0.0_dp
         kin_grad = 0.0_dp
         DO point_pos = 1, nroute_points
            local_row = features%route_send_local_rows(point_pos)
            CPASSERT(local_row >= 1 .AND. local_row <= nflat_local)
            base = nroute_grad_per_point*(point_pos - 1)
            IF (features%uses_collapsed_rks_dynamic) THEN
               density_grad(local_row, :) = density_grad(local_row, :) + &
                                            recv_grad_buffer(base + 1)
               grad_grad(local_row, 1, :) = grad_grad(local_row, 1, :) + &
                                            recv_grad_buffer(base + 2)
               grad_grad(local_row, 2, :) = grad_grad(local_row, 2, :) + &
                                            recv_grad_buffer(base + 3)
               grad_grad(local_row, 3, :) = grad_grad(local_row, 3, :) + &
                                            recv_grad_buffer(base + 4)
               kin_grad(local_row, :) = kin_grad(local_row, :) + recv_grad_buffer(base + 5)
            ELSE
               density_grad(local_row, :) = density_grad(local_row, :) + &
                                            recv_grad_buffer(base + 1:base + 2)
               grad_grad(local_row, 1, 1) = grad_grad(local_row, 1, 1) + &
                                            recv_grad_buffer(base + 3)
               grad_grad(local_row, 2, 1) = grad_grad(local_row, 2, 1) + &
                                            recv_grad_buffer(base + 4)
               grad_grad(local_row, 3, 1) = grad_grad(local_row, 3, 1) + &
                                            recv_grad_buffer(base + 5)
               grad_grad(local_row, 1, 2) = grad_grad(local_row, 1, 2) + &
                                            recv_grad_buffer(base + 6)
               grad_grad(local_row, 2, 2) = grad_grad(local_row, 2, 2) + &
                                            recv_grad_buffer(base + 7)
               grad_grad(local_row, 3, 2) = grad_grad(local_row, 3, 2) + &
                                            recv_grad_buffer(base + 8)
               kin_grad(local_row, :) = kin_grad(local_row, :) + &
                                        recv_grad_buffer(base + 9:base + 10)
            END IF
         END DO
         CALL timestop(phase_handle)

         DEALLOCATE (recv_grad_buffer, route_grad_return_recv_counts, &
                     route_grad_return_recv_displs, route_grad_return_send_counts, &
                     route_grad_return_send_displs, send_grad_buffer)
      ELSE
         ALLOCATE (chunk_grad_buffer(MAX(1, ngrad_per_point*features%chunk_feature_count)), &
                   global_grad_buffer(ngrad_per_point*features%nflat))
         IF (features%chunk_feature_count > 0) THEN
            CALL timeset("skala_gpw_grad_torch_pack", phase_handle)
            CALL pack_atom_chunk_grads(features, chunk_grad_buffer, .FALSE.)
            CALL timestop(phase_handle)
         END IF

         CALL timeset("skala_gpw_grad_allgatherv", phase_handle)
         CALL group%allgatherv(chunk_grad_buffer, global_grad_buffer, &
                               features%chunk_grad_counts, features%chunk_grad_displs)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_scatter", phase_handle)
         ALLOCATE (density_grad(nflat_local, 2), grad_grad(nflat_local, 3, 2), &
                   kin_grad(nflat_local, 2))
         density_grad = 0.0_dp
         grad_grad = 0.0_dp
         kin_grad = 0.0_dp
         local_row = 0
         DO k = LBOUND(features%feature_index, 3), UBOUND(features%feature_index, 3)
            DO j = LBOUND(features%feature_index, 2), UBOUND(features%feature_index, 2)
               DO i = LBOUND(features%feature_index, 1), UBOUND(features%feature_index, 1)
                  local_row = local_row + 1
                  DO feature_pos = features%local_feature_offsets(local_row), &
                     features%local_feature_offsets(local_row + 1) - 1
                     row = features%local_feature_rows(feature_pos)
                     CPASSERT(row >= 1 .AND. row <= features%nflat)
                     base = ngrad_per_point*(row - 1)
                     density_grad(local_row, :) = density_grad(local_row, :) + &
                                                  global_grad_buffer(base + 1:base + 2)
                     grad_grad(local_row, 1, 1) = grad_grad(local_row, 1, 1) + &
                                                  global_grad_buffer(base + 3)
                     grad_grad(local_row, 2, 1) = grad_grad(local_row, 2, 1) + &
                                                  global_grad_buffer(base + 4)
                     grad_grad(local_row, 3, 1) = grad_grad(local_row, 3, 1) + &
                                                  global_grad_buffer(base + 5)
                     grad_grad(local_row, 1, 2) = grad_grad(local_row, 1, 2) + &
                                                  global_grad_buffer(base + 6)
                     grad_grad(local_row, 2, 2) = grad_grad(local_row, 2, 2) + &
                                                  global_grad_buffer(base + 7)
                     grad_grad(local_row, 3, 2) = grad_grad(local_row, 3, 2) + &
                                                  global_grad_buffer(base + 8)
                     kin_grad(local_row, :) = kin_grad(local_row, :) + &
                                              global_grad_buffer(base + 9:base + 10)
                  END DO
               END DO
            END DO
         END DO
         CALL timestop(phase_handle)
         DEALLOCATE (chunk_grad_buffer, global_grad_buffer)

      END IF

   END SUBROUTINE fetch_and_gather_atom_chunk_grads

! **************************************************************************************************
!> \brief Build the native SKALA XC virial from feature gradients.
!> \param virial_xc ...
!> \param rho_set ...
!> \param rho_r ...
!> \param grad_grad ...
! **************************************************************************************************
   SUBROUTINE build_virial_from_feature_grads(virial_xc, rho_set, rho_r, grad_grad)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: virial_xc
      TYPE(xc_rho_set_type), INTENT(IN)                  :: rho_set
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: grad_grad

      INTEGER                                            :: i, idir, ipt, ispin, j, jdir, k, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(KIND=dp)                                      :: grad_i, tmp
      TYPE(cp_3d_r_cp_type), DIMENSION(3)                :: drho, drhoa, drhob

      nspins = SIZE(rho_r)
      bo = rho_r(1)%pw_grid%bounds_local
      ipt = 0

      IF (nspins == 1) THEN
         CALL xc_rho_set_get(rho_set, drho=drho)
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  ipt = ipt + 1
                  DO idir = 1, 3
                     grad_i = 0.5_dp*(grad_grad(ipt, idir, 1) + grad_grad(ipt, idir, 2))
                     DO jdir = 1, idir
                        tmp = -grad_i*drho(jdir)%array(i, j, k)
                        virial_xc(jdir, idir) = virial_xc(jdir, idir) + tmp
                        virial_xc(idir, jdir) = virial_xc(jdir, idir)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      ELSE
         CALL xc_rho_set_get(rho_set, drhoa=drhoa, drhob=drhob)
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  ipt = ipt + 1
                  DO idir = 1, 3
                     DO jdir = 1, idir
                        tmp = 0.0_dp
                        DO ispin = 1, 2
                           IF (ispin == 1) THEN
                              tmp = tmp - grad_grad(ipt, idir, ispin)*drhoa(jdir)%array(i, j, k)
                           ELSE
                              tmp = tmp - grad_grad(ipt, idir, ispin)*drhob(jdir)%array(i, j, k)
                           END IF
                        END DO
                        virial_xc(jdir, idir) = virial_xc(jdir, idir) + tmp
                        virial_xc(idir, jdir) = virial_xc(jdir, idir)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END IF

   END SUBROUTINE build_virial_from_feature_grads

! **************************************************************************************************
!> \brief Print a native SKALA XC virial contribution for diagnostics.
!> \param label ...
!> \param delta ...
!> \param root_rank ...
! **************************************************************************************************
   SUBROUTINE print_virial_delta(label, delta, root_rank)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: delta
      LOGICAL, INTENT(IN)                                :: root_rank

      INTEGER                                            :: i, iw

      IF (.NOT. root_rank) RETURN
      iw = cp_logger_get_default_io_unit()
      IF (iw <= 0) RETURN
      WRITE (iw, "(T2,A,1X,A)") "SKALA_GPW| XC virial contribution", TRIM(label)
      DO i = 1, 3
         WRITE (iw, "(T2,A,1X,3ES20.10)") "SKALA_GPW|", delta(i, 1:3)
      END DO

   END SUBROUTINE print_virial_delta

! **************************************************************************************************
!> \brief Add explicit SKALA coordinate-feature contributions to the XC virial.
!> \param virial_xc ...
!> \param features ...
!> \param atom_coord_grad_t ...
!> \param grid_coord_grad_t ...
!> \param root_rank ...
!> \param print_components ...
! **************************************************************************************************
   SUBROUTINE build_static_coordinate_virial(virial_xc, features, atom_coord_grad_t, &
                                             grid_coord_grad_t, root_rank, print_components)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: virial_xc
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(torch_tensor_type), INTENT(INOUT)             :: atom_coord_grad_t, grid_coord_grad_t
      LOGICAL, INTENT(IN)                                :: root_rank
      LOGICAL, INTENT(IN), OPTIONAL                      :: print_components

      INTEGER                                            :: feature_pos, i, iatom, idir, iw, j, &
                                                            jdir, k, local_row, row
      LOGICAL                                            :: my_print_components
      REAL(KIND=dp)                                      :: tmp
      REAL(KIND=dp), DIMENSION(3, 3)                     :: atom_virial, grid_virial
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: atom_coord_grad, grid_coord_grad

      my_print_components = .FALSE.
      IF (PRESENT(print_components)) my_print_components = print_components

      NULLIFY (atom_coord_grad, grid_coord_grad)
      CALL torch_tensor_grad(features%grid_coords_t, grid_coord_grad_t)
      CALL torch_tensor_data_ptr(grid_coord_grad_t, grid_coord_grad)
      CALL torch_tensor_data_ptr(atom_coord_grad_t, atom_coord_grad)

      grid_virial = 0.0_dp
      atom_virial = 0.0_dp
      local_row = 0
      DO k = LBOUND(features%feature_index, 3), UBOUND(features%feature_index, 3)
         DO j = LBOUND(features%feature_index, 2), UBOUND(features%feature_index, 2)
            DO i = LBOUND(features%feature_index, 1), UBOUND(features%feature_index, 1)
               local_row = local_row + 1
               DO feature_pos = features%local_feature_offsets(local_row), &
                  features%local_feature_offsets(local_row + 1) - 1
                  row = features%local_feature_rows(feature_pos)
                  DO idir = 1, 3
                     DO jdir = 1, 3
                        tmp = grid_coord_grad(idir, row)*features%grid_coords(jdir, row)
                        grid_virial(idir, jdir) = grid_virial(idir, jdir) + tmp
                        virial_xc(idir, jdir) = virial_xc(idir, jdir) + tmp
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      CPASSERT(local_row == features%nflat_local)

      IF (root_rank) THEN
         DO iatom = 1, SIZE(features%coarse_0_atomic_coords, 2)
            DO idir = 1, 3
               DO jdir = 1, 3
                  tmp = atom_coord_grad(idir, iatom)*features%coarse_0_atomic_coords(jdir, iatom)
                  atom_virial(idir, jdir) = atom_virial(idir, jdir) + tmp
                  virial_xc(idir, jdir) = virial_xc(idir, jdir) + tmp
               END DO
            END DO
         END DO
      END IF

      IF (my_print_components .AND. root_rank) THEN
         iw = cp_logger_get_default_io_unit()
         IF (iw > 0) THEN
            CALL print_virial_delta("static-grid", grid_virial, .TRUE.)
            CALL print_virial_delta("static-atom", atom_virial, .TRUE.)
         END IF
      END IF

      CALL torch_tensor_release(grid_coord_grad_t)

   END SUBROUTINE build_static_coordinate_virial

! **************************************************************************************************
!> \brief Add residual SKALA weight-feature contributions to the XC virial.
!> \param virial_xc ...
!> \param features ...
!> \param exc ...
!> \param grid_weight_grad_t ...
!> \param atomic_grid_weight_grad_t ...
!> \param root_rank ...
!> \param print_components ...
! **************************************************************************************************
   SUBROUTINE build_weight_virial(virial_xc, features, exc, grid_weight_grad_t, &
                                  atomic_grid_weight_grad_t, root_rank, print_components)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: virial_xc
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      REAL(KIND=dp), INTENT(IN)                          :: exc
      TYPE(torch_tensor_type), INTENT(INOUT)             :: grid_weight_grad_t, &
                                                            atomic_grid_weight_grad_t
      LOGICAL, INTENT(IN)                                :: root_rank
      LOGICAL, INTENT(IN), OPTIONAL                      :: print_components

      INTEGER                                            :: feature_pos, i, idir, iw, j, k, &
                                                            local_row, row
      LOGICAL                                            :: my_print_components
      REAL(KIND=dp)                                      :: atomic_tmp, exc_tmp, grid_tmp, tmp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atomic_grid_weight_grad, grid_weight_grad

      my_print_components = .FALSE.
      IF (PRESENT(print_components)) my_print_components = print_components

      NULLIFY (atomic_grid_weight_grad, grid_weight_grad)
      CALL torch_tensor_grad(features%grid_weights_t, grid_weight_grad_t)
      CALL torch_tensor_grad(features%atomic_grid_weights_t, atomic_grid_weight_grad_t)
      CALL torch_tensor_data_ptr(grid_weight_grad_t, grid_weight_grad)
      CALL torch_tensor_data_ptr(atomic_grid_weight_grad_t, atomic_grid_weight_grad)

      grid_tmp = 0.0_dp
      atomic_tmp = 0.0_dp
      local_row = 0
      DO k = LBOUND(features%feature_index, 3), UBOUND(features%feature_index, 3)
         DO j = LBOUND(features%feature_index, 2), UBOUND(features%feature_index, 2)
            DO i = LBOUND(features%feature_index, 1), UBOUND(features%feature_index, 1)
               local_row = local_row + 1
               DO feature_pos = features%local_feature_offsets(local_row), &
                  features%local_feature_offsets(local_row + 1) - 1
                  row = features%local_feature_rows(feature_pos)
                  grid_tmp = grid_tmp + grid_weight_grad(row)*features%grid_weights(row)
                  atomic_tmp = atomic_tmp + &
                               atomic_grid_weight_grad(row)*features%atomic_grid_weights(row)
               END DO
            END DO
         END DO
      END DO
      CPASSERT(local_row == features%nflat_local)
      exc_tmp = 0.0_dp
      IF (root_rank) exc_tmp = -exc
      tmp = grid_tmp + atomic_tmp + exc_tmp

      IF (my_print_components .AND. root_rank) THEN
         iw = cp_logger_get_default_io_unit()
         IF (iw > 0) THEN
            WRITE (iw, "(T2,A,1X,ES20.10)") "SKALA_GPW| XC virial weight grid", grid_tmp
            WRITE (iw, "(T2,A,1X,ES20.10)") "SKALA_GPW| XC virial weight atomic", atomic_tmp
            WRITE (iw, "(T2,A,1X,ES20.10)") "SKALA_GPW| XC virial weight final", exc_tmp
            WRITE (iw, "(T2,A,1X,ES20.10)") "SKALA_GPW| XC virial weight residual", tmp
         END IF
      END IF

      DO idir = 1, 3
         virial_xc(idir, idir) = virial_xc(idir, idir) + tmp
      END DO

      CALL torch_tensor_release(grid_weight_grad_t)
      CALL torch_tensor_release(atomic_grid_weight_grad_t)

   END SUBROUTINE build_weight_virial

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

      INTEGER                                            :: atom_rows_max, atom_rows_min, &
                                                            chunk_rows_max, chunk_rows_min, iw
      REAL(KIND=dp)                                      :: chunk_imbalance

      IF (.NOT. print_active) RETURN

      iw = cp_logger_get_default_io_unit()
      IF (iw <= 0) RETURN
      WRITE (UNIT=iw, FMT="(/,T2,A,1X,ES19.11)") &
         "SKALA_GPW| Native grid feature electrons", features%electron_count
      WRITE (UNIT=iw, FMT="(T2,A,1X,ES19.11)") &
         "SKALA_GPW| Native grid feature spin moment", features%spin_moment
      WRITE (UNIT=iw, FMT="(T2,A,1X,ES19.11)") &
         "SKALA_GPW| Native grid feature weight sum", features%grid_weight_sum
      IF (ALLOCATED(features%atomic_grid_sizes)) THEN
         atom_rows_min = INT(MINVAL(features%atomic_grid_sizes))
         atom_rows_max = INT(MAXVAL(features%atomic_grid_sizes))
         WRITE (UNIT=iw, FMT="(T2,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)") &
            "SKALA_GPW| Native grid atom row range", atom_rows_min, "to", &
            atom_rows_max, "sum", INT(SUM(features%atomic_grid_sizes))
      END IF
      IF (features%uses_atom_chunks) THEN
         WRITE (UNIT=iw, FMT="(T2,A,1X,I0,1X,A,1X,I0)") &
            "SKALA_GPW| Native grid atom chunk rows", features%chunk_feature_count, &
            "of", features%nflat
         IF (ALLOCATED(features%chunk_grad_counts)) THEN
            chunk_rows_min = MINVAL(features%chunk_grad_counts)/ngrad_per_point
            chunk_rows_max = MAXVAL(features%chunk_grad_counts)/ngrad_per_point
            chunk_imbalance = REAL(chunk_rows_max, KIND=dp)/REAL(MAX(1, chunk_rows_min), KIND=dp)
            WRITE (UNIT=iw, FMT="(T2,A,1X,I0,1X,A,1X,I0,1X,A,1X,ES12.5)") &
               "SKALA_GPW| Native grid atom chunk row range", chunk_rows_min, &
               "to", chunk_rows_max, "imbalance", chunk_imbalance
         END IF
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
         cuda_device_count = torch_cuda_device_count()
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
                       "Torch CUDA device range.")
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
      IF (iw <= 0) RETURN
      IF (selected_device >= 0) THEN
         WRITE (UNIT=iw, FMT="(/,T2,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)") &
            "SKALA_GPW| Native grid Torch CUDA device", selected_device, &
            "of", cuda_device_count, "requested", requested_device
      ELSE
         WRITE (UNIT=iw, FMT="(/,T2,A)") &
            "SKALA_GPW| Native grid Torch CUDA requested, but no Torch CUDA device is visible"
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
      LOGICAL                                            :: native_grid_use_cuda
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
         CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_USE_CUDA", l_val=native_grid_use_cuda)
         IF (native_grid_use_cuda) THEN
            CALL GET_ENVIRONMENT_VARIABLE("GAUXC_SKALA_CUDA_MODEL", model_path, STATUS=env_status)
            IF (env_status == 0 .AND. LEN_TRIM(model_path) > 0) RETURN
         END IF
         CALL GET_ENVIRONMENT_VARIABLE("GAUXC_SKALA_MODEL", model_path, STATUS=env_status)
         IF (env_status /= 0 .OR. LEN_TRIM(model_path) == 0) THEN
            IF (native_grid_use_cuda) THEN
               CALL cp_abort(__LOCATION__, &
                             "MODEL SKALA CUDA path requires GAUXC_SKALA_CUDA_MODEL or GAUXC_SKALA_MODEL")
            ELSE
               CALL cp_abort(__LOCATION__, &
                             "MODEL SKALA requires the GAUXC_SKALA_MODEL environment variable")
            END IF
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

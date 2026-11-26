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
   USE OMP_LIB,                         ONLY: omp_get_max_threads
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_array_utils,                  ONLY: cp_3d_r_cp_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE input_constants,                 ONLY: do_gpu_auto
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_comm_type
   USE particle_types,                  ONLY: particle_type
   USE pw_methods,                      ONLY: pw_scale,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE skala_feature_keys,              ONLY: &
        skala_key_atomic_grid_size_bound_shape, skala_key_atomic_grid_sizes, &
        skala_key_atomic_grid_weights, skala_key_coarse_0_atomic_coords, skala_key_density, &
        skala_key_grad, skala_key_grid_coords, skala_key_grid_weights, skala_key_kin
   USE skala_gpw_features,              ONLY: &
        allocate_empty_bound_shape, grid_coordinate, skala_execution_target_type, &
        skala_gpw_atom_subchunk_count, skala_gpw_atom_subchunk_layout, skala_gpw_feature_build, &
        skala_gpw_feature_build_atom_subchunk_bounds, skala_gpw_feature_release, &
        skala_gpw_feature_type, skala_gpw_features_release_cache, skala_model_context_type, &
        skala_state, skala_weight_grad_cache_type, &
        smooth_partition_atomic_weight_scale_derivative, smooth_partition_derivative_action, &
        smooth_partition_weights, sort_chunks_by_padded_rows
   USE skala_model_selection,           ONLY: get_gauxc_section,&
                                              get_skala_native_section,&
                                              resolve_skala_execution_target,&
                                              resolve_skala_model_path,&
                                              select_cuda_target,&
                                              xc_section_uses_gauxc_model
   USE skala_torch_api,                 ONLY: skala_torch_model_get_exc,&
                                              skala_torch_model_get_exc_device,&
                                              skala_torch_model_get_exc_segments,&
                                              skala_torch_model_load,&
                                              skala_torch_model_release
   USE string_utilities,                ONLY: uppercase
   USE torch_api,                       ONLY: &
        staged_upload, torch_cuda_device_count, torch_dict_create, torch_dict_release, &
        torch_dict_set, torch_dict_type, torch_file_stamp, torch_get_num_threads, &
        torch_grad_pool_clear, torch_model_metadata_clear, torch_pinned_finish, &
        torch_pinned_release, torch_pinned_reserve, torch_set_default_threads, &
        torch_stream_synchronize, torch_tensor_accumulate, torch_tensor_backward_scalar, &
        torch_tensor_clone, torch_tensor_data_ptr, torch_tensor_from_array, torch_tensor_grad, &
        torch_tensor_grad_batch, torch_tensor_grad_batch3, torch_tensor_grad_batch3_deferred, &
        torch_tensor_grad_batch_deferred, torch_tensor_hessian_vector, torch_tensor_item_double, &
        torch_tensor_narrow, torch_tensor_release, torch_tensor_set_requires_grad, &
        torch_tensor_to_device_leaf, torch_tensor_type, torch_use_cuda
   USE xc_input_constants,              ONLY: skala_gapw_atom_composite_grid,&
                                              skala_gapw_common_grid,&
                                              skala_gapw_cp2k_default,&
                                              skala_gapw_direct_valence,&
                                              skala_gapw_paw_one_center,&
                                              skala_gapw_paw_one_center_split,&
                                              skala_gpw_atom_partition_hard,&
                                              skala_gpw_atom_partition_smooth
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_get,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
   USE xc_util,                         ONLY: xc_pw_divergence,&
                                              xc_requires_tmp_g
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   INTEGER, PARAMETER, PRIVATE          :: atom_chunk_auto_max_rows = 400000, &
                                           atom_chunk_auto_min_rows = 100000, &
                                           atom_chunk_auto_row_quantum = 100000, &
                                           ncollapsed_grad_per_point = 5, &
                                           ngrad_per_point = 10
   REAL(KIND=dp), PARAMETER, PRIVATE    :: weight_cache_cell_tol = 1.0E-12_dp
   INTEGER, PARAMETER, PUBLIC           :: skala_gapw_density_partition_hard_minus_soft = 1, &
                                           skala_gapw_density_partition_hard_only = 2, &
                                           skala_gapw_density_partition_soft_only = 3, &
                                           skala_gapw_density_partition_none = 4
   PUBLIC :: build_vxc_from_feature_grads, ensure_native_skala_grid_scope, &
             skala_gapw_atom_composite_energy, skala_gapw_atom_vxc_of_r_pair, &
             native_skala_gapw_atom_composite_requested, &
             native_skala_uses_atom_composite_grid, &
             native_skala_gapw_composite_direct_ao, native_skala_gapw_composite_reference, &
             native_skala_gapw_density_partition, skala_gpw_eval, skala_gpw_weight_derivative, &
             skala_execution_use_cuda, skala_gpw_release_caches, skala_gapw_representation, &
             skala_gpw_fxc, skala_gpw_fxc_prep

CONTAINS

! **************************************************************************************************
!> \brief Align the Torch CPU intra-op thread count with CP2K's per-process OpenMP budget.
! **************************************************************************************************
   SUBROUTINE ensure_skala_torch_threads()
      INTEGER                                            :: iw, omp_threads, torch_threads

      IF (skala_state%logged_torch_threads) RETURN
      skala_state%logged_torch_threads = .TRUE.
      omp_threads = omp_get_max_threads()
      CALL torch_set_default_threads(omp_threads)
      torch_threads = torch_get_num_threads()
      iw = cp_logger_get_default_io_unit()
      IF (iw > 0) THEN
         WRITE (UNIT=iw, FMT="(T2,A,I0,A,I0,A)") &
            "SKALA| Torch intra-op threads: ", torch_threads, &
            " (CP2K OpenMP threads per rank: ", omp_threads, ")"
      END IF

   END SUBROUTINE ensure_skala_torch_threads
! **************************************************************************************************
!> \brief Read the shared SKALA_NATIVE chunk/partition controls in one place.
!> \param xc_section ...
!> \param atom_chunks ...
!> \param max_rows ...
!> \param atom_partition optional mapped partition id
! **************************************************************************************************
   SUBROUTINE read_skala_native_controls(xc_section, atom_chunks, max_rows, atom_partition)
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(OUT)                               :: atom_chunks
      INTEGER, INTENT(OUT)                               :: max_rows
      INTEGER, INTENT(OUT), OPTIONAL                     :: atom_partition

      INTEGER                                            :: partition
      TYPE(section_vals_type), POINTER                   :: native_section

      atom_chunks = .TRUE.
      max_rows = -1
      partition = skala_gpw_atom_partition_smooth
      native_section => get_skala_native_section(xc_section)
      IF (ASSOCIATED(native_section)) THEN
         CALL section_vals_val_get(native_section, "ATOM_CHUNKS", l_val=atom_chunks)
         CALL section_vals_val_get(native_section, "ATOM_CHUNK_MAX_ROWS", i_val=max_rows)
         CALL section_vals_val_get(native_section, "ATOM_PARTITION", i_val=partition)
      END IF
      IF (max_rows < -1) THEN
         CALL cp_abort(__LOCATION__, &
                       "SKALA_NATIVE%ATOM_CHUNK_MAX_ROWS must be -1, zero, or positive.")
      END IF
      IF (PRESENT(atom_partition)) atom_partition = partition

   END SUBROUTINE read_skala_native_controls
! **************************************************************************************************
!> \brief Return the pseudopotential GAPW representation selected for an active model.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION skala_gapw_representation(xc_section) RESULT(representation)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      INTEGER                                            :: representation

      TYPE(section_vals_type), POINTER                   :: gauxc_section, native_section

      representation = skala_gapw_cp2k_default
      native_section => get_skala_native_section(xc_section)
      IF (ASSOCIATED(native_section)) THEN
         CALL section_vals_val_get(native_section, "PSEUDOPOTENTIAL_GAPW_REPRESENTATION", &
                                   i_val=representation)
      ELSE
         IF (.NOT. xc_section_uses_gauxc_model(xc_section)) RETURN
         gauxc_section => get_gauxc_section(xc_section)
         IF (.NOT. ASSOCIATED(gauxc_section)) RETURN
         CALL section_vals_val_get(gauxc_section, "PSEUDOPOTENTIAL_GAPW_REPRESENTATION", &
                                   i_val=representation)
      END IF

      SELECT CASE (representation)
      CASE (skala_gapw_direct_valence, skala_gapw_paw_one_center, &
            skala_gapw_cp2k_default, skala_gapw_paw_one_center_split)
         CONTINUE
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, "Unknown pseudopotential GAPW representation.")
      END SELECT

   END FUNCTION skala_gapw_representation

! **************************************************************************************************
!> \brief Return true if native SKALA should use the full GAPW ORB density on one common grid.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION native_skala_gapw_composite_reference(xc_section) RESULT(use_composite_reference)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: use_composite_reference

      TYPE(section_vals_type), POINTER                   :: native_section

      use_composite_reference = .FALSE.
      native_section => get_skala_native_section(xc_section)
      IF (ASSOCIATED(native_section)) THEN
         CALL section_vals_val_get(native_section, "GAPW_COMPOSITE_REFERENCE", &
                                   l_val=use_composite_reference)
      END IF

   END FUNCTION native_skala_gapw_composite_reference

! **************************************************************************************************
!> \brief Return true when the explicit atom-centered composite reference is requested.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION native_skala_gapw_atom_composite_requested(xc_section) RESULT(use_atom_composite)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: use_atom_composite

      TYPE(section_vals_type), POINTER                   :: native_section

      use_atom_composite = .FALSE.
      native_section => get_skala_native_section(xc_section)
      IF (ASSOCIATED(native_section)) THEN
         CALL section_vals_val_get(native_section, &
                                   "GAPW_ATOM_COMPOSITE_REFERENCE", &
                                   l_val=use_atom_composite)
      END IF

   END FUNCTION native_skala_gapw_atom_composite_requested

! **************************************************************************************************
!> \brief Return true when native Skala uses atom-centered grids.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION native_skala_uses_atom_composite_grid(xc_section) RESULT(use_atom_composite)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: use_atom_composite

      INTEGER                                            :: composite_grid
      TYPE(section_vals_type), POINTER                   :: native_section

      use_atom_composite = .FALSE.
      composite_grid = skala_gapw_atom_composite_grid
      native_section => get_skala_native_section(xc_section)
      IF (ASSOCIATED(native_section)) THEN
         CALL section_vals_val_get(native_section, "LAYOUT", &
                                   i_val=composite_grid)
      END IF
      SELECT CASE (composite_grid)
      CASE (skala_gapw_common_grid)
         CONTINUE
      CASE (skala_gapw_atom_composite_grid)
         use_atom_composite = .TRUE.
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, "Unknown native-grid layout.")
      END SELECT

   END FUNCTION native_skala_uses_atom_composite_grid

! **************************************************************************************************
!> \brief Return true if the GAPW composite reference uses direct full-ORB collocation.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION native_skala_gapw_composite_direct_ao(xc_section) RESULT(use_direct_ao)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: use_direct_ao

      TYPE(section_vals_type), POINTER                   :: native_section

      use_direct_ao = .FALSE.
      native_section => get_skala_native_section(xc_section)
      IF (ASSOCIATED(native_section)) THEN
         CALL section_vals_val_get(native_section, "GAPW_COMPOSITE_DIRECT_AO", &
                                   l_val=use_direct_ao)
      END IF
      use_direct_ao = native_skala_gapw_composite_reference(xc_section) .AND. use_direct_ao

   END FUNCTION native_skala_gapw_composite_direct_ao

! **************************************************************************************************
!> \brief Return the hard/soft GAPW one-center density partition for native SKALA.
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION native_skala_gapw_density_partition(xc_section) RESULT(partition)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      INTEGER                                            :: partition

      TYPE(section_vals_type), POINTER                   :: native_section

      partition = skala_gapw_density_partition_hard_minus_soft
      native_section => get_skala_native_section(xc_section)
      IF (ASSOCIATED(native_section)) THEN
         CALL section_vals_val_get(native_section, "GAPW_DENSITY_PARTITION", &
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
                       "Unknown SKALA_NATIVE%GAPW_DENSITY_PARTITION value.")
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
      TYPE(section_vals_type), POINTER                   :: functionals, native_section, xc_fun

      NULLIFY (native_section)
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
         IF (xc_fun%section%name == "SKALA_NATIVE") native_section => xc_fun
      END DO

      IF (.NOT. ASSOCIATED(native_section)) THEN
         CPABORT("Native SKALA GPW requires an XC_FUNCTIONAL%SKALA_NATIVE section")
      END IF
      IF (nfun /= 1) THEN
         CPABORT("SKALA_NATIVE must be the only XC functional")
      END IF

      CALL section_vals_val_get(native_section, "MODEL", c_val=model_name)
      model_key = ADJUSTL(model_name)
      CALL uppercase(model_key)
      IF (TRIM(model_key) == "NONE" .OR. TRIM(model_key) == "") THEN
         CPABORT("SKALA_NATIVE requires MODEL SKALA or a TorchScript model path")
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

      INTEGER :: i, ipt, ispin, iw, j, k, native_grid_atom_chunk_max_rows, &
         native_grid_atom_partition, native_grid_atom_subchunks, nspins, phase_handle, &
         selected_cuda_device, static_grad_atom_slot, static_grad_atomic_weight_slot, &
         static_grad_count, static_grad_grid_slot, static_grad_weight_slot, xc_deriv_method_id, &
         xc_rho_smooth_id
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL :: has_atom_chunk_work, lsd, my_just_energy, native_grid_atom_chunks, &
         native_grid_diagnostics, needs_atom_force, use_atom_subchunks, use_chunk_static_grad, &
         use_full_static_grad
      REAL(KIND=dp)                                      :: density_contraction, tau_contraction
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET   :: atomic_grid_weight_grad_g, &
                                                            grid_weight_grad_g
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         TARGET                                          :: atom_coord_grad_g, grid_coord_grad_g
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_before
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atomic_grid_weight_grad_p, &
                                                            grid_weight_grad_p
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: atom_coord_grad_p, density_grad, &
                                                            grid_coord_grad_p, kin_grad
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grad_grad
      TYPE(section_vals_type), POINTER                   :: native_section
      TYPE(skala_gpw_feature_type), POINTER              :: features
      TYPE(torch_tensor_type)                            :: exc_tensor
      TYPE(torch_tensor_type), DIMENSION(4)              :: static_grad_batch_t, &
                                                            static_grad_tensor_batch
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set

      virial_xc = 0.0_dp
      exc = 0.0_dp
      my_just_energy = .FALSE.
      IF (PRESENT(just_energy)) my_just_energy = just_energy
      needs_atom_force = PRESENT(atom_force)
      IF (needs_atom_force) atom_force = 0.0_dp

      IF (compute_virial .AND. my_just_energy) THEN
         CALL cp_abort(__LOCATION__, &
                       "Native SKALA GPW stress/virial requires feature gradients.")
      END IF

      IF (.NOT. my_just_energy) THEN
         IF (compute_virial) THEN
            skala_state%stress_mode_entered = .TRUE.
         ELSE IF (skala_state%stress_mode_entered) THEN
            CALL cp_abort(__LOCATION__, &
                          "SKALA: a stress/virial evaluation was performed earlier in this "// &
                          "run, but a non-stress evaluation was requested. Mixing stress and "// &
                          "non-stress modes within one run is not supported.")
         END IF
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
      CALL ensure_model_context(xc_section, rho_r(1)%pw_grid%para%group, &
                                cuda_device=selected_cuda_device)
      CALL read_skala_native_controls(xc_section, native_grid_atom_chunks, &
                                      native_grid_atom_chunk_max_rows, native_grid_atom_partition)
      native_section => get_skala_native_section(xc_section)
      IF (needs_atom_force .OR. compute_virial) THEN
         IF (native_grid_atom_partition == skala_gpw_atom_partition_hard) THEN
            CALL cp_abort(__LOCATION__, &
                          "SKALA_NATIVE%ATOM_PARTITION HARD supports energy/VXC only. "// &
                          "Set ATOM_PARTITION SMOOTH to compute forces or the virial.")
         END IF
      END IF

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

      CALL section_vals_val_get(native_section, "DIAGNOSTICS", l_val=native_grid_diagnostics)
      features => skala_state%features
      CALL skala_gpw_feature_build(features, rho_set, rho_r, particle_set, cell, &
                                   requires_grad=(.NOT. my_just_energy), weights=weights, &
                                   requires_coordinate_grad=(needs_atom_force .OR. compute_virial), &
                                   requires_stress_grad=compute_virial, &
                                   use_atom_chunks=native_grid_atom_chunks, &
                                   compute_diagnostics=native_grid_diagnostics, &
                                   atom_partition=native_grid_atom_partition, &
                                   cuda_device=selected_cuda_device)
      CALL ensure_grad_scratch(features%nflat_local)
      density_grad => skala_state%grad_scratch%density_grad
      grad_grad => skala_state%grad_scratch%grad_grad
      kin_grad => skala_state%grad_scratch%kin_grad
      use_chunk_static_grad = features%uses_atom_chunks .AND. &
                              (needs_atom_force .OR. compute_virial)
      IF (native_grid_diagnostics) THEN
         CALL print_native_grid_diagnostics(features, rho_r(1)%pw_grid%para%group%mepos == 0, &
                                            rho_r(1)%pw_grid%para%group)
      END IF

      IF (features%uses_atom_chunks .AND. native_grid_atom_chunk_max_rows == -1) THEN
         native_grid_atom_chunk_max_rows = auto_atom_chunk_max_rows(features, &
                                                                    rho_r(1)%pw_grid%para%group)
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
         CALL skala_torch_model_get_exc(skala_state%model, features%inputs, &
                                        features%grid_weights_t, exc_tensor, exc)
      END IF
      IF (features%uses_atom_chunks) CALL rho_r(1)%pw_grid%para%group%sum(exc)

      IF (.NOT. my_just_energy) THEN
         IF (.NOT. use_atom_subchunks) THEN
            IF (has_atom_chunk_work) THEN
               CALL timeset("skala_gpw_backward", phase_handle)
               CALL torch_tensor_backward_scalar(exc_tensor)
               CALL timestop(phase_handle)
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
         IF (use_chunk_static_grad) THEN
            CALL timeset("skala_gpw_coord_grad_fetch", phase_handle)
            CALL fetch_and_route_atom_chunk_grads( &
               features, rho_r(1)%pw_grid%para%group, weight_only=.FALSE., &
               local_grid_coord_grad=grid_coord_grad_g, &
               local_grid_weight_grad=grid_weight_grad_g, &
               local_atomic_grid_weight_grad=atomic_grid_weight_grad_g, &
               atom_coord_grad=atom_coord_grad_g)
            grid_coord_grad_p => grid_coord_grad_g
            grid_weight_grad_p => grid_weight_grad_g
            atomic_grid_weight_grad_p => atomic_grid_weight_grad_g
            atom_coord_grad_p => atom_coord_grad_g
            CALL timestop(phase_handle)
         END IF
         use_full_static_grad = .NOT. features%uses_atom_chunks .AND. &
                                (needs_atom_force .OR. compute_virial)
         IF (use_full_static_grad) THEN
            CALL timeset("skala_gpw_static_grad_fetch", phase_handle)
            NULLIFY (grid_coord_grad_p, grid_weight_grad_p, atomic_grid_weight_grad_p, &
                     atom_coord_grad_p)
            static_grad_count = 0
            static_grad_atom_slot = 0
            static_grad_grid_slot = 0
            static_grad_weight_slot = 0
            static_grad_atomic_weight_slot = 0
            IF (needs_atom_force .OR. compute_virial) THEN
               static_grad_count = static_grad_count + 1
               static_grad_atom_slot = static_grad_count
               static_grad_tensor_batch(static_grad_count) = features%coarse_0_atomic_coords_t
            END IF
            IF (compute_virial) THEN
               static_grad_count = static_grad_count + 1
               static_grad_grid_slot = static_grad_count
               static_grad_tensor_batch(static_grad_count) = features%grid_coords_t
            END IF
            IF (compute_virial .OR. &
                (needs_atom_force .AND. &
                 features%atom_partition == skala_gpw_atom_partition_smooth)) THEN
               static_grad_count = static_grad_count + 1
               static_grad_weight_slot = static_grad_count
               static_grad_tensor_batch(static_grad_count) = features%grid_weights_t
               static_grad_count = static_grad_count + 1
               static_grad_atomic_weight_slot = static_grad_count
               static_grad_tensor_batch(static_grad_count) = features%atomic_grid_weights_t
            END IF
            IF (static_grad_count > 0) THEN
               CALL torch_tensor_grad_batch(static_grad_tensor_batch(1:static_grad_count), &
                                            static_grad_batch_t(1:static_grad_count))
               IF (static_grad_atom_slot > 0) THEN
                  CALL torch_tensor_data_ptr(static_grad_batch_t(static_grad_atom_slot), &
                                             atom_coord_grad_p)
               END IF
               IF (static_grad_grid_slot > 0) THEN
                  CALL torch_tensor_data_ptr(static_grad_batch_t(static_grad_grid_slot), &
                                             grid_coord_grad_p)
               END IF
               IF (static_grad_weight_slot > 0) THEN
                  CALL torch_tensor_data_ptr(static_grad_batch_t(static_grad_weight_slot), &
                                             grid_weight_grad_p)
               END IF
               IF (static_grad_atomic_weight_slot > 0) THEN
                  CALL torch_tensor_data_ptr(static_grad_batch_t(static_grad_atomic_weight_slot), &
                                             atomic_grid_weight_grad_p)
               END IF
            END IF
            CALL timestop(phase_handle)
         END IF
         IF (features%atom_partition == skala_gpw_atom_partition_smooth .AND. &
             (needs_atom_force .OR. compute_virial) .AND. ASSOCIATED(grid_weight_grad_p)) THEN
            CALL store_weight_grad_cache(grid_weight_grad_p, features, particle_set, cell, &
                                         rho_r, tau)
         END IF
         IF (compute_virial) THEN
            IF (native_grid_diagnostics) virial_before = virial_xc
            CALL build_weight_virial(virial_xc, features, exc, &
                                     rho_r(1)%pw_grid%para%group%mepos == 0, &
                                     native_grid_diagnostics, &
                                     grid_weight_grad_in=grid_weight_grad_p, &
                                     atomic_grid_weight_grad_in=atomic_grid_weight_grad_p)
            IF (native_grid_diagnostics) THEN
               CALL print_virial_delta("weight-residual", virial_xc - virial_before, &
                                       rho_r(1)%pw_grid%para%group%mepos == 0)
            END IF
         END IF
         IF (needs_atom_force) THEN
            CALL add_explicit_coordinate_force(atom_force, &
                                               rho_r(1)%pw_grid%para%group%mepos == 0, &
                                               atom_coord_grad_in=atom_coord_grad_p)
         END IF

         IF ((needs_atom_force .OR. compute_virial) .AND. &
             features%atom_partition == skala_gpw_atom_partition_smooth) THEN
            IF (compute_virial .AND. native_grid_diagnostics) virial_before = virial_xc
            IF (needs_atom_force .AND. compute_virial) THEN
               CALL add_smooth_partition_force_virial( &
                  atom_force=atom_force, virial_xc=virial_xc, features=features, &
                  particle_set=particle_set, cell=cell, rho_r=rho_r, &
                  grid_weight_grad_in=grid_weight_grad_p, &
                  atomic_grid_weight_grad_in=atomic_grid_weight_grad_p)
            ELSE IF (needs_atom_force) THEN
               CALL add_smooth_partition_force_virial( &
                  atom_force=atom_force, features=features, particle_set=particle_set, &
                  cell=cell, rho_r=rho_r, &
                  grid_weight_grad_in=grid_weight_grad_p, &
                  atomic_grid_weight_grad_in=atomic_grid_weight_grad_p)
            ELSE
               CALL add_smooth_partition_force_virial( &
                  virial_xc=virial_xc, features=features, particle_set=particle_set, &
                  cell=cell, rho_r=rho_r, &
                  grid_weight_grad_in=grid_weight_grad_p, &
                  atomic_grid_weight_grad_in=atomic_grid_weight_grad_p)
            END IF
            IF (compute_virial .AND. native_grid_diagnostics) THEN
               CALL print_virial_delta("smooth-partition", virial_xc - virial_before, &
                                       rho_r(1)%pw_grid%para%group%mepos == 0)
               virial_before = virial_xc
            END IF
         END IF

         CALL timeset("skala_gpw_vxc_unpack", phase_handle)
         IF (native_grid_diagnostics) THEN
            bo = rho_r(1)%pw_grid%bounds_local
            density_contraction = 0.0_dp
            tau_contraction = 0.0_dp
            ipt = 0
            DO k = bo(1, 3), bo(2, 3)
               DO j = bo(1, 2), bo(2, 2)
                  DO i = bo(1, 1), bo(2, 1)
                     ipt = ipt + 1
                     IF (nspins == 1) THEN
                        density_contraction = density_contraction + rho_r(1)%array(i, j, k)* &
                                              0.5_dp*(density_grad(ipt, 1) + density_grad(ipt, 2))
                        tau_contraction = tau_contraction + tau(1)%array(i, j, k)* &
                                          0.5_dp*(kin_grad(ipt, 1) + kin_grad(ipt, 2))
                     ELSE
                        DO ispin = 1, nspins
                           density_contraction = density_contraction + &
                                                 rho_r(ispin)%array(i, j, k)*density_grad(ipt, ispin)
                           tau_contraction = tau_contraction + &
                                             tau(ispin)%array(i, j, k)*kin_grad(ipt, ispin)
                        END DO
                     END IF
                  END DO
               END DO
            END DO
            CALL rho_r(1)%pw_grid%para%group%sum(density_contraction)
            CALL rho_r(1)%pw_grid%para%group%sum(tau_contraction)
            IF (rho_r(1)%pw_grid%para%group%mepos == 0) THEN
               iw = cp_logger_get_default_io_unit()
               IF (iw > 0) THEN
                  WRITE (iw, "(T2,A,1X,ES20.10)") &
                     "SKALA_GPW| XC density-gradient contraction", density_contraction
                  WRITE (iw, "(T2,A,1X,ES20.10)") &
                     "SKALA_GPW| XC kinetic-gradient contraction", tau_contraction
               END IF
            END IF
         END IF
         IF (compute_virial) THEN
            IF (native_grid_diagnostics) virial_before = virial_xc
            CALL build_virial_from_feature_grads(virial_xc, rho_set, rho_r, grad_grad)
            IF (native_grid_diagnostics) THEN
               CALL print_virial_delta("feature-gradient", virial_xc - virial_before, &
                                       rho_r(1)%pw_grid%para%group%mepos == 0)
               virial_before = virial_xc
            END IF
            CALL build_static_coordinate_virial(virial_xc, features, &
                                                rho_r(1)%pw_grid%para%group%mepos == 0, &
                                                native_grid_diagnostics, &
                                                atom_coord_grad_in=atom_coord_grad_p, &
                                                grid_coord_grad_in=grid_coord_grad_p)
            IF (native_grid_diagnostics) THEN
               CALL print_virial_delta("static-coordinates", virial_xc - virial_before, &
                                       rho_r(1)%pw_grid%para%group%mepos == 0)
               virial_before = virial_xc
            END IF
         END IF
         CALL build_vxc_from_feature_grads(vxc_rho, vxc_tau, rho_r, pw_pool, &
                                           density_grad, grad_grad, kin_grad, &
                                           xc_deriv_method_id)
         CALL timestop(phase_handle)

         CALL timeset("skala_gpw_grad_release", phase_handle)
         IF (use_full_static_grad) THEN
            DO i = 1, static_grad_count
               CALL torch_tensor_release(static_grad_batch_t(i))
            END DO
         END IF
         IF (ALLOCATED(atom_coord_grad_g)) DEALLOCATE (atom_coord_grad_g)
         IF (ALLOCATED(grid_coord_grad_g)) DEALLOCATE (grid_coord_grad_g)
         IF (ALLOCATED(grid_weight_grad_g)) DEALLOCATE (grid_weight_grad_g)
         IF (ALLOCATED(atomic_grid_weight_grad_g)) DEALLOCATE (atomic_grid_weight_grad_g)
         CALL timestop(phase_handle)
      END IF

      CALL timeset("skala_gpw_cleanup", phase_handle)
      IF (.NOT. use_atom_subchunks .AND. has_atom_chunk_work) CALL torch_tensor_release(exc_tensor)
      CALL skala_gpw_feature_release(features, keep_host_arrays=.TRUE.)
      CALL xc_rho_set_release(rho_set, pw_pool=pw_pool)
      CALL timestop(phase_handle)

   END SUBROUTINE skala_gpw_eval

! **************************************************************************************************
!> \brief Evaluate the native SKALA XC energy density on the real-space grid.
!>        On return weight_deriv_r holds the unweighted energy density
!>        sum_atom w_atom * eps_atom, i.e. the quantity xc_density expects before
!>        the caller multiplies by dvol. Do not interpret it as dE/d(multiplier):
!>        the optional weights pointer only enters the feature grid weights and is
!>        expected to be unassociated for the accurate-XCINT force path.
!> \param weight_deriv_r ...
!> \param rho_r ...
!> \param rho_g ...
!> \param tau ...
!> \param xc_section ...
!> \param weights optional external integration weights (may be unassociated)
!> \param pw_pool ...
!> \param particle_set ...
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE skala_gpw_weight_derivative(weight_deriv_r, rho_r, rho_g, tau, xc_section, &
                                          weights, pw_pool, particle_set, cell)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: weight_deriv_r
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: tau
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell

      INTEGER :: chunk_max_rows, feature_begin, feature_end, feature_pos, i, iatom, isub, j, k, &
         local_row, native_grid_atom_partition, natom, nspins, nsubchunks, nvals, row, &
         selected_cuda_device, xc_deriv_method_id, xc_rho_smooth_id
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: subchunk_atom_begin, &
                                                            subchunk_atom_count, &
                                                            subchunk_row_begin, subchunk_row_count
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: cache_hit, have_exc_tensor, &
                                                            have_weight_grad_t, lsd, &
                                                            native_grid_atom_chunks
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: included
      REAL(KIND=dp)                                      :: local_derivative
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: local_weight_grad, partition_weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_coords_pbc
      REAL(KIND=dp), DIMENSION(3)                        :: grid_point
      REAL(KIND=dp), DIMENSION(:), POINTER               :: grid_weight_grad
      TYPE(skala_gpw_feature_type)                       :: subchunk
      TYPE(skala_gpw_feature_type), POINTER              :: features
      TYPE(torch_tensor_type)                            :: exc_tensor, grid_weight_grad_t
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set

      CPASSERT(ASSOCIATED(rho_r))
      CPASSERT(ASSOCIATED(rho_g))
      CPASSERT(ASSOCIATED(tau))
      CALL pw_zero(weight_deriv_r)

      nspins = SIZE(rho_r)
      lsd = (nspins /= 1)
      CALL ensure_model_context(xc_section, rho_r(1)%pw_grid%para%group, &
                                cuda_device=selected_cuda_device)
      CALL read_skala_native_controls(xc_section, native_grid_atom_chunks, chunk_max_rows, &
                                      native_grid_atom_partition)
      IF (.NOT. native_grid_atom_chunks) chunk_max_rows = 0
      IF (native_grid_atom_partition == skala_gpw_atom_partition_hard) THEN
         CALL cp_abort(__LOCATION__, &
                       "SKALA_NATIVE%ATOM_PARTITION HARD supports energy/VXC only. "// &
                       "Set ATOM_PARTITION SMOOTH to compute forces or the virial.")
      END IF

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

      features => skala_state%features

      have_exc_tensor = .FALSE.
      have_weight_grad_t = .FALSE.
      NULLIFY (grid_weight_grad)

      ! Try the cache before building/uploading the feature tensors: the key
      ! (positions, cell, density/tau signatures, spin/atom mode and partition)
      ! uniquely fixes the layout, so a hit can reuse the cached gradient and skip
      ! the feature build, upload and model forward/backward.
      cache_hit = .FALSE.
      IF (skala_state%layout%active .AND. skala_state%weight_grad%valid) THEN
         cache_hit = weight_grad_cache_matches(rho_r(1)%pw_grid%para%group, particle_set, cell, &
                                               rho_r, tau, nspins, native_grid_atom_chunks, &
                                               native_grid_atom_partition, &
                                               skala_state%weight_grad%nvals)
      END IF
      IF (cache_hit) THEN
         nvals = skala_state%weight_grad%nvals
      ELSE
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
                                      requires_weight_grad=.TRUE., &
                                      use_atom_chunks=native_grid_atom_chunks, &
                                      atom_partition=native_grid_atom_partition, &
                                      cuda_device=selected_cuda_device)
         IF (features%uses_atom_chunks) THEN
            nvals = SIZE(features%route_send_local_features)
         ELSE
            nvals = features%nflat
         END IF
      END IF

      IF (cache_hit) THEN
         IF (native_grid_atom_chunks) THEN
            ALLOCATE (local_weight_grad(nvals))
            local_weight_grad(:) = skala_state%weight_grad%gradient
         ELSE
            grid_weight_grad => skala_state%weight_grad%gradient
         END IF
      ELSE IF (features%uses_atom_chunks) THEN
         IF (chunk_max_rows == -1) THEN
            chunk_max_rows = auto_atom_chunk_max_rows(features, rho_r(1)%pw_grid%para%group)
         END IF
         nsubchunks = 1
         IF (chunk_max_rows > 0) THEN
            IF (features%chunk_feature_count > 0) THEN
               nsubchunks = skala_gpw_atom_subchunk_count(chunk_max_rows)
            END IF
            CALL rho_r(1)%pw_grid%para%group%max(nsubchunks)
         END IF
         IF (nsubchunks > 1) THEN
            CALL skala_gpw_atom_subchunk_layout(chunk_max_rows, subchunk_atom_begin, &
                                                subchunk_atom_count, subchunk_row_begin, &
                                                subchunk_row_count)
            DO isub = 1, SIZE(subchunk_atom_begin)
               CALL skala_gpw_feature_build_atom_subchunk_bounds( &
                  features, subchunk, subchunk_atom_begin(isub), subchunk_atom_count(isub), &
                  subchunk_row_begin(isub), subchunk_row_count(isub))
               CALL skala_torch_model_get_exc_device(skala_state%model, subchunk%inputs, &
                                                     subchunk%grid_weights_t, exc_tensor)
               CALL torch_tensor_backward_scalar(exc_tensor)
               CALL torch_tensor_release(exc_tensor)
               CALL skala_gpw_feature_release(subchunk)
            END DO
            DEALLOCATE (subchunk_atom_begin, subchunk_atom_count, subchunk_row_begin, &
                        subchunk_row_count)
         ELSE IF (features%chunk_feature_count > 0) THEN
            CALL skala_torch_model_get_exc_device(skala_state%model, features%inputs, &
                                                  features%grid_weights_t, exc_tensor)
            CALL torch_tensor_backward_scalar(exc_tensor)
            have_exc_tensor = .TRUE.
         END IF
         CALL fetch_and_route_atom_chunk_grads(features, rho_r(1)%pw_grid%para%group, &
                                               weight_only=.TRUE., &
                                               local_grid_weight_grad=local_weight_grad)
      ELSE
         CALL skala_torch_model_get_exc_device(skala_state%model, features%inputs, &
                                               features%grid_weights_t, exc_tensor)
         CALL torch_tensor_backward_scalar(exc_tensor)
         CALL torch_tensor_grad(features%grid_weights_t, grid_weight_grad_t)
         CALL torch_tensor_data_ptr(grid_weight_grad_t, grid_weight_grad)
         have_exc_tensor = .TRUE.
         have_weight_grad_t = .TRUE.
      END IF

      natom = SIZE(particle_set)
      ALLOCATE (atom_coords_pbc(3, natom), included(natom), partition_weights(natom))
      DO iatom = 1, natom
         atom_coords_pbc(:, iatom) = pbc(particle_set(iatom)%r, cell, positive_range=.TRUE.)
      END DO

      bo = rho_r(1)%pw_grid%bounds_local
      local_row = 0
      DO k = bo(1, 3), bo(2, 3)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 1), bo(2, 1)
               local_row = local_row + 1
               feature_begin = skala_state%layout%local_feature_offsets(local_row)
               feature_end = skala_state%layout%local_feature_offsets(local_row + 1) - 1
               local_derivative = 0.0_dp
               grid_point = grid_coordinate(rho_r(1)%pw_grid, [i, j, k])
               CALL smooth_partition_weights(grid_point, atom_coords_pbc, cell, &
                                             partition_weights, included)
               CPASSERT(feature_end - feature_begin + 1 == COUNT(included))
               feature_pos = feature_begin
               DO iatom = 1, natom
                  IF (.NOT. included(iatom)) CYCLE
                  IF (native_grid_atom_chunks) THEN
                     ! Atom-chunk gradients are indexed in identity feature order.
                     local_derivative = local_derivative + &
                                        partition_weights(iatom)*local_weight_grad(feature_pos)
                  ELSE
                     row = skala_state%layout%local_feature_rows(feature_pos)
                     local_derivative = local_derivative + &
                                        partition_weights(iatom)*grid_weight_grad(row)
                  END IF
                  feature_pos = feature_pos + 1
               END DO
               CPASSERT(feature_pos == feature_end + 1)
               weight_deriv_r%array(i, j, k) = local_derivative
            END DO
         END DO
      END DO
      CPASSERT(local_row == skala_state%layout%nflat_local)

      DEALLOCATE (atom_coords_pbc, included, partition_weights)
      IF (have_weight_grad_t) CALL torch_tensor_release(grid_weight_grad_t)
      IF (have_exc_tensor) CALL torch_tensor_release(exc_tensor)
      IF (ALLOCATED(local_weight_grad)) DEALLOCATE (local_weight_grad)
      IF (.NOT. cache_hit) THEN
         CALL skala_gpw_feature_release(features, keep_host_arrays=.TRUE.)
         CALL xc_rho_set_release(rho_set, pw_pool=pw_pool)
      END IF

   END SUBROUTINE skala_gpw_weight_derivative

! **************************************************************************************************
!> \brief Evaluate the hard and soft one-center densities of one atom in a single model pass.
!>        Both blocks share the atom grid but keep separate energies, potentials and grads.
!> \param xc_section ...
!> \param grid_atom ...
!> \param group ...
!> \param atom_coord ...
!> \param rho_h ...
!> \param drho_h ...
!> \param tau_h ...
!> \param rho_s ...
!> \param drho_s ...
!> \param tau_s ...
!> \param weight_h ...
!> \param weight_s ...
!> \param lsd ...
!> \param nspins ...
!> \param na ...
!> \param nr ...
!> \param evaluate_hard ...
!> \param evaluate_soft ...
!> \param exc_h ...
!> \param vxc_h ...
!> \param vxg_h ...
!> \param vtau_h ...
!> \param exc_s ...
!> \param vxc_s ...
!> \param vxg_s ...
!> \param vtau_s ...
!> \param energy_only ...
!> \param atom_force_h ...
!> \param atom_virial_h ...
!> \param atom_force_s ...
!> \param atom_virial_s ...
! **************************************************************************************************
   SUBROUTINE skala_gapw_atom_vxc_of_r_pair(xc_section, grid_atom, group, atom_coord, &
                                            rho_h, drho_h, tau_h, rho_s, drho_s, tau_s, &
                                            weight_h, weight_s, lsd, nspins, na, nr, &
                                            evaluate_hard, evaluate_soft, &
                                            exc_h, vxc_h, vxg_h, vtau_h, &
                                            exc_s, vxc_s, vxg_s, vtau_s, &
                                            energy_only, atom_force_h, atom_virial_h, &
                                            atom_force_s, atom_virial_s)
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(grid_atom_type), POINTER                      :: grid_atom

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: atom_coord
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho_h, tau_h, vxc_h, vtau_h, &
                                                            rho_s, tau_s, vxc_s, vtau_s
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: drho_h, vxg_h, drho_s, vxg_s
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: weight_h, weight_s
      LOGICAL, INTENT(IN)                                :: lsd
      INTEGER, INTENT(IN)                                :: nspins, na, nr
      LOGICAL, INTENT(IN)                                :: evaluate_hard, evaluate_soft
      REAL(KIND=dp), INTENT(OUT)                         :: exc_h, exc_s
      LOGICAL, INTENT(IN), OPTIONAL                      :: energy_only
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT), &
         OPTIONAL                                        :: atom_force_h, atom_force_s
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT), &
         OPTIONAL                                        :: atom_virial_h, atom_virial_s

      INTEGER                                            :: ia, idir, iblock, ir, jdir, nblocks, &
                                                            nflat, nflat_tot, &
                                                            row, row_in_block, selected_cuda_device
      INTEGER(KIND=int_8)                                :: upload_offset, upload_total
      LOGICAL :: is_hard, need_coord_grad, my_energy_only, native_grid_use_cuda
      REAL(KIND=dp)                                      :: tmp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: atom_coord_grad, density_grad, &
                                                            grid_coord_grad, kin_grad
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho_b, tau_b, vtau_b, vxc_b, grad_grad
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: drho_b, vxg_b
      TYPE(torch_tensor_type), DIMENSION(2)              :: grad_batch_t
      TYPE(torch_tensor_type)                            :: density_grad_t, exc_tensor, grad_grad_t, &
                                                            kin_grad_t
      INTEGER                                            :: iblock_h, iblock_s

      CPASSERT(ASSOCIATED(xc_section))
      CPASSERT(ASSOCIATED(grid_atom))
      CPASSERT(ASSOCIATED(rho_h))
      CPASSERT(ASSOCIATED(drho_h))
      CPASSERT(ASSOCIATED(tau_h))
      CPASSERT(ASSOCIATED(rho_s))
      CPASSERT(ASSOCIATED(drho_s))
      CPASSERT(ASSOCIATED(tau_s))

      my_energy_only = .FALSE.
      IF (PRESENT(energy_only)) my_energy_only = energy_only
      need_coord_grad = PRESENT(atom_force_h) .OR. PRESENT(atom_virial_h) .OR. &
                        PRESENT(atom_force_s) .OR. PRESENT(atom_virial_s)
      exc_h = 0.0_dp
      exc_s = 0.0_dp
      IF (PRESENT(atom_force_h)) atom_force_h = 0.0_dp
      IF (PRESENT(atom_force_s)) atom_force_s = 0.0_dp
      IF (PRESENT(atom_virial_h)) atom_virial_h = 0.0_dp
      IF (PRESENT(atom_virial_s)) atom_virial_s = 0.0_dp
      IF (.NOT. my_energy_only) THEN
         vxc_h = 0.0_dp
         vxg_h = 0.0_dp
         vtau_h = 0.0_dp
         vxc_s = 0.0_dp
         vxg_s = 0.0_dp
         vtau_s = 0.0_dp
      END IF

      iblock_h = 0
      iblock_s = 0
      nblocks = 0
      IF (evaluate_hard) THEN
         nblocks = nblocks + 1
         iblock_h = nblocks
      END IF
      IF (evaluate_soft) THEN
         nblocks = nblocks + 1
         iblock_s = nblocks
      END IF
      IF (nblocks == 0) RETURN

      CALL ensure_model_context(xc_section, group, native_grid_use_cuda, selected_cuda_device)

      nflat = na*nr
      nflat_tot = nblocks*nflat
      IF (skala_state%atom_ws%nflat_tot /= nflat_tot .OR. skala_state%atom_ws%na /= na .OR. skala_state%atom_ws%nr /= nr .OR. &
          skala_state%atom_ws%nspins /= nspins .OR. skala_state%atom_ws%nblocks /= nblocks .OR. &
          skala_state%atom_ws%device /= selected_cuda_device) THEN
         CALL release_atom_ws()
         ALLOCATE (skala_state%atom_ws%density(nflat_tot, 2), skala_state%atom_ws%grad(nflat_tot, 3, 2), &
                   skala_state%atom_ws%kin(nflat_tot, 2), skala_state%atom_ws%grid_coords(3, nflat_tot), &
                   skala_state%atom_ws%grid_weights(nflat_tot), skala_state%atom_ws%atomic_grid_weights(nflat_tot), &
                   skala_state%atom_ws%atomic_grid_sizes(nblocks), &
                   skala_state%atom_ws%coarse_0_atomic_coords(3, nblocks), &
                   skala_state%atom_ws%atomic_grid_size_bound_shape(0, nflat), skala_state%atom_ws%weight_b(na, nr), &
                   skala_state%atom_ws%seg_offsets(nblocks), skala_state%atom_ws%seg_lengths(nblocks), &
                   skala_state%atom_ws%exc_segments(nblocks))
         skala_state%atom_ws%na = na
         skala_state%atom_ws%nr = nr
         skala_state%atom_ws%nspins = nspins
         skala_state%atom_ws%nblocks = nblocks
         skala_state%atom_ws%nflat_tot = nflat_tot
         skala_state%atom_ws%device = selected_cuda_device
         CALL torch_tensor_from_array(skala_state%atom_ws%atomic_grid_sizes_t, skala_state%atom_ws%atomic_grid_sizes)
         CALL torch_tensor_to_device_leaf(skala_state%atom_ws%atomic_grid_sizes_t, .FALSE.)
         CALL torch_tensor_from_array(skala_state%atom_ws%atomic_grid_size_bound_shape_t, &
                                      skala_state%atom_ws%atomic_grid_size_bound_shape)
         CALL torch_tensor_to_device_leaf(skala_state%atom_ws%atomic_grid_size_bound_shape_t, .FALSE.)
         skala_state%atom_ws%int_tensors_active = .TRUE.
      END IF
      skala_state%atom_ws%density = 0.0_dp
      skala_state%atom_ws%grad = 0.0_dp
      skala_state%atom_ws%kin = 0.0_dp
      skala_state%atom_ws%grid_coords = 0.0_dp
      skala_state%atom_ws%grid_weights = 0.0_dp
      skala_state%atom_ws%atomic_grid_weights = 0.0_dp

      DO iblock = 1, nblocks
         is_hard = (evaluate_hard .AND. iblock == iblock_h)
         IF (is_hard) THEN
            rho_b => rho_h
            drho_b => drho_h
            tau_b => tau_h
            skala_state%atom_ws%weight_b(:, :) = weight_h
         ELSE
            rho_b => rho_s
            drho_b => drho_s
            tau_b => tau_s
            skala_state%atom_ws%weight_b(:, :) = weight_s
         END IF
         skala_state%atom_ws%coarse_0_atomic_coords(:, iblock) = atom_coord
         skala_state%atom_ws%atomic_grid_sizes(iblock) = INT(nflat, KIND=int_8)
         DO ir = 1, nr
            DO ia = 1, na
               row_in_block = (ir - 1)*na + ia
               row = (iblock - 1)*nflat + row_in_block
               skala_state%atom_ws%grid_coords(1, row) = atom_coord(1) + grid_atom%rad(ir)* &
                                                         grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia)
               skala_state%atom_ws%grid_coords(2, row) = atom_coord(2) + grid_atom%rad(ir)* &
                                                         grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia)
               skala_state%atom_ws%grid_coords(3, row) = atom_coord(3) + grid_atom%rad(ir)*grid_atom%cos_pol(ia)
               skala_state%atom_ws%grid_weights(row) = skala_state%atom_ws%weight_b(ia, ir)
               skala_state%atom_ws%atomic_grid_weights(row) = skala_state%atom_ws%weight_b(ia, ir)
               IF (nspins == 1) THEN
                  skala_state%atom_ws%density(row, :) = 0.5_dp*rho_b(ia, ir, 1)
                  DO idir = 1, 3
                     skala_state%atom_ws%grad(row, idir, :) = 0.5_dp*drho_b(idir, ia, ir, 1)
                  END DO
                  skala_state%atom_ws%kin(row, :) = 0.5_dp*tau_b(ia, ir, 1)
               ELSE
                  skala_state%atom_ws%density(row, :) = rho_b(ia, ir, 1:2)
                  DO idir = 1, 3
                     skala_state%atom_ws%grad(row, idir, :) = drho_b(idir, ia, ir, 1:2)
                  END DO
                  skala_state%atom_ws%kin(row, :) = tau_b(ia, ir, 1:2)
               END IF
            END DO
         END DO
      END DO

      IF (native_grid_use_cuda) THEN
         upload_total = INT(SIZE(skala_state%atom_ws%grid_coords), KIND=int_8) + &
                        INT(SIZE(skala_state%atom_ws%grid_weights), KIND=int_8) + &
                        INT(SIZE(skala_state%atom_ws%atomic_grid_weights), KIND=int_8) + &
                        INT(SIZE(skala_state%atom_ws%coarse_0_atomic_coords), KIND=int_8) + &
                        INT(SIZE(skala_state%atom_ws%density), KIND=int_8) + INT(SIZE(skala_state%atom_ws%grad), KIND=int_8) + &
                        INT(SIZE(skala_state%atom_ws%kin), KIND=int_8)
         CALL torch_pinned_reserve(upload_total)
         upload_offset = 0_int_8
         upload_offset = staged_upload(skala_state%atom_ws%grid_coords_t, skala_state%atom_ws%grid_coords, upload_offset, &
                                       need_coord_grad)
         upload_offset = staged_upload(skala_state%atom_ws%grid_weights_t, skala_state%atom_ws%grid_weights, upload_offset, .FALSE.)
         upload_offset = staged_upload(skala_state%atom_ws%atomic_grid_weights_t, skala_state%atom_ws%atomic_grid_weights, &
                                       upload_offset, .FALSE.)
         upload_offset = staged_upload(skala_state%atom_ws%coarse_0_atomic_coords_t, skala_state%atom_ws%coarse_0_atomic_coords, &
                                       upload_offset, need_coord_grad)
         upload_offset = staged_upload(skala_state%atom_ws%density_t, skala_state%atom_ws%density, upload_offset, &
                                       .NOT. my_energy_only)
         upload_offset = staged_upload(skala_state%atom_ws%grad_t, skala_state%atom_ws%grad, upload_offset,.NOT. my_energy_only)
         upload_offset = staged_upload(skala_state%atom_ws%kin_t, skala_state%atom_ws%kin, upload_offset,.NOT. my_energy_only)
         CALL torch_stream_synchronize()
         CALL torch_pinned_finish()
         skala_state%atom_ws%tensors_active = .TRUE.
      ELSE
         CALL torch_tensor_from_array(skala_state%atom_ws%grid_coords_t, skala_state%atom_ws%grid_coords)
         CALL torch_tensor_to_device_leaf(skala_state%atom_ws%grid_coords_t, need_coord_grad)
         CALL torch_tensor_from_array(skala_state%atom_ws%grid_weights_t, skala_state%atom_ws%grid_weights)
         CALL torch_tensor_to_device_leaf(skala_state%atom_ws%grid_weights_t, .FALSE.)
         CALL torch_tensor_from_array(skala_state%atom_ws%atomic_grid_weights_t, skala_state%atom_ws%atomic_grid_weights)
         CALL torch_tensor_to_device_leaf(skala_state%atom_ws%atomic_grid_weights_t, .FALSE.)
         CALL torch_tensor_from_array(skala_state%atom_ws%coarse_0_atomic_coords_t, skala_state%atom_ws%coarse_0_atomic_coords)
         CALL torch_tensor_to_device_leaf(skala_state%atom_ws%coarse_0_atomic_coords_t, need_coord_grad)
         CALL torch_tensor_from_array(skala_state%atom_ws%density_t, skala_state%atom_ws%density)
         CALL torch_tensor_to_device_leaf(skala_state%atom_ws%density_t,.NOT. my_energy_only)
         CALL torch_tensor_from_array(skala_state%atom_ws%grad_t, skala_state%atom_ws%grad)
         CALL torch_tensor_to_device_leaf(skala_state%atom_ws%grad_t,.NOT. my_energy_only)
         CALL torch_tensor_from_array(skala_state%atom_ws%kin_t, skala_state%atom_ws%kin)
         CALL torch_tensor_to_device_leaf(skala_state%atom_ws%kin_t,.NOT. my_energy_only)
      END IF

      IF (.NOT. skala_state%atom_ws%dict_active) THEN
         CALL torch_dict_create(skala_state%atom_ws%inputs)
         skala_state%atom_ws%dict_active = .TRUE.
      END IF
      CALL torch_dict_set(skala_state%atom_ws%inputs, skala_key_grid_coords, skala_state%atom_ws%grid_coords_t)
      CALL torch_dict_set(skala_state%atom_ws%inputs, skala_key_grid_weights, skala_state%atom_ws%grid_weights_t)
      CALL torch_dict_set(skala_state%atom_ws%inputs, skala_key_atomic_grid_weights, &
                          skala_state%atom_ws%atomic_grid_weights_t)
      CALL torch_dict_set(skala_state%atom_ws%inputs, skala_key_atomic_grid_sizes, skala_state%atom_ws%atomic_grid_sizes_t)
      CALL torch_dict_set(skala_state%atom_ws%inputs, skala_key_atomic_grid_size_bound_shape, &
                          skala_state%atom_ws%atomic_grid_size_bound_shape_t)
      CALL torch_dict_set(skala_state%atom_ws%inputs, skala_key_density, skala_state%atom_ws%density_t)
      CALL torch_dict_set(skala_state%atom_ws%inputs, skala_key_grad, skala_state%atom_ws%grad_t)
      CALL torch_dict_set(skala_state%atom_ws%inputs, skala_key_kin, skala_state%atom_ws%kin_t)
      CALL torch_dict_set(skala_state%atom_ws%inputs, skala_key_coarse_0_atomic_coords, &
                          skala_state%atom_ws%coarse_0_atomic_coords_t)

      DO iblock = 1, nblocks
         skala_state%atom_ws%seg_offsets(iblock) = INT((iblock - 1)*nflat, KIND=int_8)
         skala_state%atom_ws%seg_lengths(iblock) = INT(nflat, KIND=int_8)
      END DO
      CALL skala_torch_model_get_exc_segments(skala_state%model, skala_state%atom_ws%inputs, &
                                              skala_state%atom_ws%grid_weights_t, &
                                              skala_state%atom_ws%seg_offsets, skala_state%atom_ws%seg_lengths, &
                                              exc_tensor, skala_state%atom_ws%exc_segments)
      IF (evaluate_hard) exc_h = skala_state%atom_ws%exc_segments(iblock_h)
      IF (evaluate_soft) exc_s = skala_state%atom_ws%exc_segments(iblock_s)

      IF (.NOT. my_energy_only) THEN
         NULLIFY (atom_coord_grad, density_grad, grad_grad, grid_coord_grad, kin_grad)
         CALL torch_tensor_backward_scalar(exc_tensor)
         IF (need_coord_grad) THEN
            CALL torch_tensor_grad_batch_deferred( &
               [skala_state%atom_ws%grid_coords_t, skala_state%atom_ws%coarse_0_atomic_coords_t], grad_batch_t)
            CALL torch_tensor_data_ptr(grad_batch_t(1), grid_coord_grad)
            CALL torch_tensor_data_ptr(grad_batch_t(2), atom_coord_grad)
         END IF
         CALL torch_tensor_grad_batch3_deferred(skala_state%atom_ws%density_t, &
                                                skala_state%atom_ws%grad_t, skala_state%atom_ws%kin_t, &
                                                density_grad_t, grad_grad_t, kin_grad_t)
         CALL torch_tensor_data_ptr(density_grad_t, density_grad)
         CALL torch_tensor_data_ptr(grad_grad_t, grad_grad)
         CALL torch_tensor_data_ptr(kin_grad_t, kin_grad)
         CALL torch_stream_synchronize()

         DO iblock = 1, nblocks
            is_hard = (evaluate_hard .AND. iblock == iblock_h)
            IF (is_hard) THEN
               vxc_b => vxc_h
               vxg_b => vxg_h
               vtau_b => vtau_h
            ELSE
               vxc_b => vxc_s
               vxg_b => vxg_s
               vtau_b => vtau_s
            END IF
            IF (need_coord_grad) THEN
               IF (is_hard .AND. PRESENT(atom_force_h)) THEN
                  atom_force_h(:) = atom_coord_grad(:, iblock)
                  DO row_in_block = 1, nflat
                     atom_force_h(:) = atom_force_h(:) + &
                                       grid_coord_grad(:, (iblock - 1)*nflat + row_in_block)
                  END DO
               ELSE IF (.NOT. is_hard .AND. PRESENT(atom_force_s)) THEN
                  atom_force_s(:) = atom_coord_grad(:, iblock)
                  DO row_in_block = 1, nflat
                     atom_force_s(:) = atom_force_s(:) + &
                                       grid_coord_grad(:, (iblock - 1)*nflat + row_in_block)
                  END DO
               END IF
               IF (is_hard .AND. PRESENT(atom_virial_h)) THEN
                  DO row_in_block = 1, nflat
                     row = (iblock - 1)*nflat + row_in_block
                     DO idir = 1, 3
                        DO jdir = 1, 3
                           tmp = grid_coord_grad(idir, row)*skala_state%atom_ws%coarse_0_atomic_coords(jdir, iblock)
                           atom_virial_h(idir, jdir) = atom_virial_h(idir, jdir) + tmp
                        END DO
                     END DO
                  END DO
                  DO idir = 1, 3
                     DO jdir = 1, 3
                        tmp = atom_coord_grad(idir, iblock)*skala_state%atom_ws%coarse_0_atomic_coords(jdir, iblock)
                        atom_virial_h(idir, jdir) = atom_virial_h(idir, jdir) + tmp
                     END DO
                  END DO
               ELSE IF (.NOT. is_hard .AND. PRESENT(atom_virial_s)) THEN
                  DO row_in_block = 1, nflat
                     row = (iblock - 1)*nflat + row_in_block
                     DO idir = 1, 3
                        DO jdir = 1, 3
                           tmp = grid_coord_grad(idir, row)*skala_state%atom_ws%coarse_0_atomic_coords(jdir, iblock)
                           atom_virial_s(idir, jdir) = atom_virial_s(idir, jdir) + tmp
                        END DO
                     END DO
                  END DO
                  DO idir = 1, 3
                     DO jdir = 1, 3
                        tmp = atom_coord_grad(idir, iblock)*skala_state%atom_ws%coarse_0_atomic_coords(jdir, iblock)
                        atom_virial_s(idir, jdir) = atom_virial_s(idir, jdir) + tmp
                     END DO
                  END DO
               END IF
            END IF

            DO ir = 1, nr
               DO ia = 1, na
                  row_in_block = (ir - 1)*na + ia
                  row = (iblock - 1)*nflat + row_in_block
                  IF (lsd) THEN
                     vxc_b(ia, ir, 1:2) = density_grad(row, 1:2)
                     DO idir = 1, 3
                        vxg_b(idir, ia, ir, 1:2) = grad_grad(row, idir, 1:2)
                     END DO
                     vtau_b(ia, ir, 1:2) = kin_grad(row, 1:2)
                  ELSE
                     vxc_b(ia, ir, 1) = 0.5_dp*(density_grad(row, 1) + density_grad(row, 2))
                     DO idir = 1, 3
                        vxg_b(idir, ia, ir, 1) = &
                           0.5_dp*(grad_grad(row, idir, 1) + grad_grad(row, idir, 2))
                     END DO
                     vtau_b(ia, ir, 1) = 0.5_dp*(kin_grad(row, 1) + kin_grad(row, 2))
                  END IF
               END DO
            END DO
         END DO

         CALL torch_tensor_release(density_grad_t)
         CALL torch_tensor_release(grad_grad_t)
         CALL torch_tensor_release(kin_grad_t)
         IF (need_coord_grad) THEN
            CALL torch_tensor_release(grad_batch_t(1))
            CALL torch_tensor_release(grad_batch_t(2))
         END IF
      END IF

      CALL torch_tensor_release(exc_tensor)
      IF (.NOT. skala_state%atom_ws%tensors_active) THEN
         CALL torch_tensor_release(skala_state%atom_ws%density_t)
         CALL torch_tensor_release(skala_state%atom_ws%grad_t)
         CALL torch_tensor_release(skala_state%atom_ws%kin_t)
         CALL torch_tensor_release(skala_state%atom_ws%grid_coords_t)
         CALL torch_tensor_release(skala_state%atom_ws%grid_weights_t)
         CALL torch_tensor_release(skala_state%atom_ws%atomic_grid_weights_t)
         CALL torch_tensor_release(skala_state%atom_ws%coarse_0_atomic_coords_t)
      END IF

   END SUBROUTINE skala_gapw_atom_vxc_of_r_pair

! **************************************************************************************************
!> \brief Evaluate a rank-local set of complete atom blocks and sum their SKALA energies.
!> \param xc_section ...
!> \param group ...
!> \param density ...
!> \param grad ...
!> \param kin ...
!> \param grid_coords ...
!> \param grid_weights ...
!> \param atomic_grid_weights ...
!> \param atomic_grid_sizes ...
!> \param atomic_coords ...
!> \param exc ...
!> \param density_grad_out ...
!> \param grad_grad_out ...
!> \param kin_grad_out ...
!> \param grid_coord_grad_out ...
!> \param grid_weight_grad_out ...
!> \param atomic_grid_weight_grad_out ...
!> \param atom_coord_grad_out ...
!> \param static_geometry_valid ...
! **************************************************************************************************
   SUBROUTINE skala_gapw_atom_composite_energy(xc_section, group, density, grad, kin, &
                                               grid_coords, grid_weights, atomic_grid_weights, &
                                               atomic_grid_sizes, atomic_coords, exc, &
                                               density_grad_out, grad_grad_out, kin_grad_out, &
                                               grid_coord_grad_out, grid_weight_grad_out, &
                                               atomic_grid_weight_grad_out, atom_coord_grad_out, &
                                               static_geometry_valid)
      TYPE(section_vals_type), POINTER                   :: xc_section

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(IN) :: density, grid_coords, &
                                                                 atomic_coords, kin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(IN) :: grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: grid_weights, atomic_grid_weights
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: atomic_grid_sizes
      REAL(KIND=dp), INTENT(OUT)                         :: exc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(OUT), OPTIONAL :: density_grad_out, &
                                                                            kin_grad_out, &
                                                                            grid_coord_grad_out, &
                                                                            atom_coord_grad_out
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(OUT), OPTIONAL :: grad_grad_out
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), INTENT(OUT), OPTIONAL :: grid_weight_grad_out, &
                                                                         atomic_grid_weight_grad_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: static_geometry_valid
      CHARACTER(len=default_path_length)                 :: model_path

      INTEGER :: atom_begin, atom_count, chunk_max_rows, chunk_row_count, chunk_row_start, &
                 ichunk, local_natom, local_nrow, max_grid_size, &
                 nchunks, phase_handle, selected_cuda_device
      INTEGER(KIND=int_8)                                :: upload_offset, upload_total
      INTEGER, DIMENSION(:), POINTER                     :: chunk_atom_begin, chunk_atom_count, &
                                                            chunk_row_begin, chunk_row_counts
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :)  :: atomic_grid_size_bound_shape
      LOGICAL :: active_rank, native_grid_atom_chunks, native_grid_use_cuda, &
                 need_coordinate_derivatives, need_derivatives, upload_static
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atomic_grid_weight_grad, grid_weight_grad
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: atom_coord_grad, density_grad, &
                                                            grid_coord_grad, kin_grad
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grad_grad
      TYPE(torch_dict_type)                              :: inputs
      TYPE(torch_tensor_type), DIMENSION(4)              :: grad_batch_t
      TYPE(torch_tensor_type)                            :: exc_total_t
      TYPE(torch_tensor_type)                            :: atomic_grid_size_bound_shape_t, &
                                                            atomic_grid_sizes_t, &
                                                            atomic_grid_sizes_chunk_t, &
                                                            atomic_grid_weights_chunk_t, &
                                                            atomic_coords_chunk_t, &
                                                            density_chunk_t, density_grad_t, &
                                                            exc_tensor, grad_chunk_t, grad_grad_t, &
                                                            grid_coords_chunk_t, grid_weights_chunk_t, &
                                                            kin_chunk_t, kin_grad_t
      TYPE(torch_tensor_type), TARGET                    :: atomic_coords_local_t, &
                                                            atomic_grid_weights_local_t, &
                                                            density_local_t, grad_local_t, &
                                                            grid_coords_local_t, grid_weights_local_t, &
                                                            kin_local_t
      TYPE(torch_tensor_type), POINTER                   :: atomic_coords_t, atomic_grid_weights_t, &
                                                            density_t, grad_t, grid_coords_t, &
                                                            grid_weights_t, kin_t

      CPASSERT(ASSOCIATED(xc_section))
      CPASSERT(SIZE(density, 1) == SIZE(grid_weights))
      CPASSERT(SIZE(density, 1) == SIZE(atomic_grid_weights))
      CPASSERT(SIZE(density, 1) == SIZE(grid_coords, 2))
      CPASSERT(SIZE(density, 1) == SIZE(grad, 1))
      CPASSERT(SIZE(density, 1) == SIZE(kin, 1))
      CPASSERT(SIZE(density, 2) == 2)
      CPASSERT(SIZE(grad, 2) == 3)
      CPASSERT(SIZE(grad, 3) == 2)
      CPASSERT(SIZE(kin, 2) == 2)
      CPASSERT(SIZE(atomic_grid_sizes) == SIZE(atomic_coords, 2))
      CPASSERT(SUM(atomic_grid_sizes) == INT(SIZE(density, 1), KIND=int_8))
      need_derivatives = PRESENT(density_grad_out)
      CPASSERT(PRESENT(grad_grad_out) .EQV. need_derivatives)
      CPASSERT(PRESENT(kin_grad_out) .EQV. need_derivatives)
      need_coordinate_derivatives = PRESENT(grid_coord_grad_out)
      CPASSERT(PRESENT(grid_weight_grad_out) .EQV. need_coordinate_derivatives)
      CPASSERT(PRESENT(atomic_grid_weight_grad_out) .EQV. need_coordinate_derivatives)
      CPASSERT(PRESENT(atom_coord_grad_out) .EQV. need_coordinate_derivatives)
      CPASSERT((.NOT. need_coordinate_derivatives) .OR. need_derivatives)

      local_nrow = SIZE(density, 1)
      local_natom = SIZE(atomic_grid_sizes)
      active_rank = local_natom > 0
      CPASSERT(active_rank .EQV. (local_nrow > 0))

      CALL get_skala_model_path(xc_section, model_path)
      native_grid_use_cuda = skala_execution_use_cuda(xc_section)
      CALL read_skala_native_controls(xc_section, native_grid_atom_chunks, chunk_max_rows)
      IF (.NOT. native_grid_atom_chunks) chunk_max_rows = 0
      IF (chunk_max_rows == -1) THEN
         chunk_max_rows = auto_atom_composite_chunk_max_rows(atomic_grid_sizes, group)
      END IF
      selected_cuda_device = configure_native_grid_cuda( &
                             native_grid_use_cuda, group)
      exc = 0.0_dp
      IF (active_rank) THEN
         CALL timeset("skala_atom_tensor_setup", phase_handle)
         CALL ensure_model_loaded(model_path, selected_cuda_device, xc_section)
         CALL ensure_atom_chunk_layout(atomic_grid_sizes, chunk_max_rows)
         chunk_atom_begin => skala_state%chunk_layout%atom_begin
         chunk_atom_count => skala_state%chunk_layout%atom_count
         chunk_row_begin => skala_state%chunk_layout%row_begin
         chunk_row_counts => skala_state%chunk_layout%row_counts
         nchunks = SIZE(chunk_atom_begin)
         CPASSERT(nchunks > 0)

         IF (native_grid_use_cuda) THEN
            grid_coords_t => skala_state%atom_cuda_tensors%grid_coords_t
            grid_weights_t => skala_state%atom_cuda_tensors%grid_weights_t
            atomic_grid_weights_t => skala_state%atom_cuda_tensors%atomic_grid_weights_t
            atomic_coords_t => skala_state%atom_cuda_tensors%atomic_coords_t
            density_t => skala_state%atom_cuda_tensors%density_t
            grad_t => skala_state%atom_cuda_tensors%grad_t
            kin_t => skala_state%atom_cuda_tensors%kin_t
            upload_static = .TRUE.
            IF (PRESENT(static_geometry_valid)) THEN
               IF (static_geometry_valid .AND. skala_state%atom_static_valid .AND. &
                   skala_state%atom_static_device == selected_cuda_device .AND. &
                   skala_state%atom_static_nrow == local_nrow .AND. &
                   skala_state%atom_static_natom == local_natom) upload_static = .FALSE.
            END IF
            IF (upload_static) THEN
               upload_total = INT(SIZE(grid_coords), KIND=int_8) + &
                              INT(SIZE(grid_weights), KIND=int_8) + &
                              INT(SIZE(atomic_grid_weights), KIND=int_8) + &
                              INT(SIZE(atomic_coords), KIND=int_8) + &
                              INT(SIZE(density), KIND=int_8) + &
                              INT(SIZE(grad), KIND=int_8) + &
                              INT(SIZE(kin), KIND=int_8)
            ELSE
               upload_total = INT(SIZE(density), KIND=int_8) + &
                              INT(SIZE(grad), KIND=int_8) + &
                              INT(SIZE(kin), KIND=int_8)
            END IF
            CALL torch_pinned_reserve(upload_total)
            upload_offset = 0_int_8
            IF (upload_static) THEN
               upload_offset = staged_upload(grid_coords_t, grid_coords, upload_offset, &
                                             need_coordinate_derivatives)
               upload_offset = staged_upload(grid_weights_t, grid_weights, upload_offset, &
                                             need_coordinate_derivatives)
               upload_offset = staged_upload(atomic_grid_weights_t, atomic_grid_weights, &
                                             upload_offset, need_coordinate_derivatives)
               upload_offset = staged_upload(atomic_coords_t, atomic_coords, upload_offset, &
                                             need_coordinate_derivatives)
            ELSE
               CALL torch_tensor_set_requires_grad(grid_coords_t, need_coordinate_derivatives)
               CALL torch_tensor_set_requires_grad(grid_weights_t, need_coordinate_derivatives)
               CALL torch_tensor_set_requires_grad(atomic_grid_weights_t, need_coordinate_derivatives)
               CALL torch_tensor_set_requires_grad(atomic_coords_t, need_coordinate_derivatives)
            END IF
            upload_offset = staged_upload(density_t, density, upload_offset, need_derivatives)
            upload_offset = staged_upload(grad_t, grad, upload_offset, need_derivatives)
            upload_offset = staged_upload(kin_t, kin, upload_offset, need_derivatives)
            CALL torch_stream_synchronize()
            CALL torch_pinned_finish()
            skala_state%atom_cuda_tensors_active = .TRUE.
            IF (upload_static) THEN
               skala_state%atom_static_valid = .TRUE.
               skala_state%atom_static_device = selected_cuda_device
               skala_state%atom_static_nrow = local_nrow
               skala_state%atom_static_natom = local_natom
            END IF
         ELSE
            grid_coords_t => grid_coords_local_t
            grid_weights_t => grid_weights_local_t
            atomic_grid_weights_t => atomic_grid_weights_local_t
            atomic_coords_t => atomic_coords_local_t
            density_t => density_local_t
            grad_t => grad_local_t
            kin_t => kin_local_t
            CALL torch_tensor_from_array(grid_coords_t, grid_coords)
            CALL torch_tensor_to_device_leaf(grid_coords_t, need_coordinate_derivatives)
            CALL torch_tensor_from_array(grid_weights_t, grid_weights)
            CALL torch_tensor_to_device_leaf(grid_weights_t, need_coordinate_derivatives)
            CALL torch_tensor_from_array(atomic_grid_weights_t, atomic_grid_weights)
            CALL torch_tensor_to_device_leaf(atomic_grid_weights_t, need_coordinate_derivatives)
            CALL torch_tensor_from_array(atomic_coords_t, atomic_coords)
            CALL torch_tensor_to_device_leaf(atomic_coords_t, need_coordinate_derivatives)
            CALL torch_tensor_from_array(density_t, density)
            CALL torch_tensor_to_device_leaf(density_t, need_derivatives)
            CALL torch_tensor_from_array(grad_t, grad)
            CALL torch_tensor_to_device_leaf(grad_t, need_derivatives)
            CALL torch_tensor_from_array(kin_t, kin)
            CALL torch_tensor_to_device_leaf(kin_t, need_derivatives)
         END IF
         CALL torch_tensor_from_array(atomic_grid_sizes_t, atomic_grid_sizes)
         CALL torch_tensor_to_device_leaf(atomic_grid_sizes_t, .FALSE.)
         CALL timestop(phase_handle)

         CALL torch_dict_create(inputs)
         DO ichunk = 1, nchunks
            atom_begin = chunk_atom_begin(ichunk)
            atom_count = chunk_atom_count(ichunk)
            chunk_row_start = chunk_row_begin(ichunk)
            chunk_row_count = chunk_row_counts(ichunk)
            max_grid_size = INT(MAXVAL(atomic_grid_sizes( &
                                       atom_begin:atom_begin + atom_count - 1)))
            CALL allocate_empty_bound_shape(atomic_grid_size_bound_shape, max_grid_size)

            CALL torch_tensor_narrow(grid_coords_t, 0, chunk_row_start - 1, &
                                     chunk_row_count, grid_coords_chunk_t)
            CALL torch_tensor_narrow(grid_weights_t, 0, chunk_row_start - 1, &
                                     chunk_row_count, grid_weights_chunk_t)
            CALL torch_tensor_narrow(atomic_grid_weights_t, 0, chunk_row_start - 1, &
                                     chunk_row_count, atomic_grid_weights_chunk_t)
            CALL torch_tensor_narrow(atomic_grid_sizes_t, 0, atom_begin - 1, &
                                     atom_count, atomic_grid_sizes_chunk_t)
            CALL torch_tensor_narrow(atomic_coords_t, 0, atom_begin - 1, &
                                     atom_count, atomic_coords_chunk_t)
            CALL torch_tensor_narrow(density_t, 1, chunk_row_start - 1, &
                                     chunk_row_count, density_chunk_t)
            CALL torch_tensor_narrow(grad_t, 2, chunk_row_start - 1, &
                                     chunk_row_count, grad_chunk_t)
            CALL torch_tensor_narrow(kin_t, 1, chunk_row_start - 1, &
                                     chunk_row_count, kin_chunk_t)
            CALL torch_tensor_from_array(atomic_grid_size_bound_shape_t, &
                                         atomic_grid_size_bound_shape)
            CALL torch_tensor_to_device_leaf(atomic_grid_size_bound_shape_t, .FALSE.)

            CALL torch_dict_set(inputs, skala_key_grid_coords, grid_coords_chunk_t)
            CALL torch_dict_set(inputs, skala_key_grid_weights, grid_weights_chunk_t)
            CALL torch_dict_set(inputs, skala_key_atomic_grid_weights, atomic_grid_weights_chunk_t)
            CALL torch_dict_set(inputs, skala_key_atomic_grid_sizes, atomic_grid_sizes_chunk_t)
            CALL torch_dict_set(inputs, skala_key_atomic_grid_size_bound_shape, &
                                atomic_grid_size_bound_shape_t)
            CALL torch_dict_set(inputs, skala_key_density, density_chunk_t)
            CALL torch_dict_set(inputs, skala_key_grad, grad_chunk_t)
            CALL torch_dict_set(inputs, skala_key_kin, kin_chunk_t)
            CALL torch_dict_set(inputs, skala_key_coarse_0_atomic_coords, atomic_coords_chunk_t)

            CALL timeset("skala_atom_torch_forward", phase_handle)
            CALL skala_torch_model_get_exc_device(skala_state%model, inputs, grid_weights_chunk_t, &
                                                  exc_tensor)
            CALL torch_tensor_accumulate(exc_total_t, exc_tensor)
            CALL timestop(phase_handle)
            IF (need_derivatives) THEN
               CALL timeset("skala_atom_torch_backward", phase_handle)
               CALL torch_tensor_backward_scalar(exc_tensor)
               CALL timestop(phase_handle)
            END IF

            CALL torch_tensor_release(exc_tensor)
            CALL torch_tensor_release(grid_coords_chunk_t)
            CALL torch_tensor_release(grid_weights_chunk_t)
            CALL torch_tensor_release(atomic_grid_weights_chunk_t)
            CALL torch_tensor_release(atomic_grid_sizes_chunk_t)
            CALL torch_tensor_release(atomic_coords_chunk_t)
            CALL torch_tensor_release(density_chunk_t)
            CALL torch_tensor_release(grad_chunk_t)
            CALL torch_tensor_release(kin_chunk_t)
            CALL torch_tensor_release(atomic_grid_size_bound_shape_t)
            DEALLOCATE (atomic_grid_size_bound_shape)
         END DO
         CALL torch_dict_release(inputs)
         exc = exc + torch_tensor_item_double(exc_total_t)
         CALL torch_tensor_release(exc_total_t)
      END IF
      CALL group%sum(exc)

      IF (need_derivatives) THEN
         ALLOCATE (density_grad_out(local_nrow, 2), grad_grad_out(local_nrow, 3, 2), &
                   kin_grad_out(local_nrow, 2))
         density_grad_out = 0.0_dp
         grad_grad_out = 0.0_dp
         kin_grad_out = 0.0_dp
         IF (active_rank) THEN
            NULLIFY (density_grad, grad_grad, kin_grad)
            CALL timeset("skala_atom_grad_fetch", phase_handle)
            CALL torch_tensor_grad_batch3_deferred(density_t, grad_t, kin_t, density_grad_t, &
                                                   grad_grad_t, kin_grad_t)
            CALL torch_tensor_data_ptr(density_grad_t, density_grad)
            CALL torch_tensor_data_ptr(grad_grad_t, grad_grad)
            CALL torch_tensor_data_ptr(kin_grad_t, kin_grad)
            CALL timestop(phase_handle)
         END IF

         IF (need_coordinate_derivatives) THEN
            ALLOCATE (grid_coord_grad_out(3, local_nrow), grid_weight_grad_out(local_nrow), &
                      atomic_grid_weight_grad_out(local_nrow), atom_coord_grad_out(3, local_natom))
            grid_coord_grad_out = 0.0_dp
            grid_weight_grad_out = 0.0_dp
            atomic_grid_weight_grad_out = 0.0_dp
            atom_coord_grad_out = 0.0_dp
            IF (active_rank) THEN
               NULLIFY (atomic_grid_weight_grad, atom_coord_grad, grid_coord_grad, grid_weight_grad)
               CALL timeset("skala_atom_coord_grad_fetch", phase_handle)
               CALL torch_tensor_grad_batch_deferred( &
                  [grid_coords_t, grid_weights_t, atomic_grid_weights_t, atomic_coords_t], &
                  grad_batch_t)
               CALL torch_tensor_data_ptr(grad_batch_t(1), grid_coord_grad)
               CALL torch_tensor_data_ptr(grad_batch_t(2), grid_weight_grad)
               CALL torch_tensor_data_ptr(grad_batch_t(3), atomic_grid_weight_grad)
               CALL torch_tensor_data_ptr(grad_batch_t(4), atom_coord_grad)
               CALL timestop(phase_handle)
            END IF
         END IF

         IF (active_rank) THEN
            CALL torch_stream_synchronize()
            density_grad_out(:, :) = density_grad
            grad_grad_out(:, :, :) = grad_grad
            kin_grad_out(:, :) = kin_grad
            IF (need_coordinate_derivatives) THEN
               grid_coord_grad_out(:, :) = grid_coord_grad
               grid_weight_grad_out(:) = grid_weight_grad
               atomic_grid_weight_grad_out(:) = atomic_grid_weight_grad
               atom_coord_grad_out(:, :) = atom_coord_grad
               CALL torch_tensor_release(grad_batch_t(1))
               CALL torch_tensor_release(grad_batch_t(2))
               CALL torch_tensor_release(grad_batch_t(3))
               CALL torch_tensor_release(grad_batch_t(4))
            END IF
            CALL torch_tensor_release(density_grad_t)
            CALL torch_tensor_release(grad_grad_t)
            CALL torch_tensor_release(kin_grad_t)
         END IF
      END IF

      IF (active_rank) THEN
         CALL timeset("skala_atom_tensor_release", phase_handle)
         IF (.NOT. native_grid_use_cuda) THEN
            CALL torch_tensor_release(density_t)
            CALL torch_tensor_release(grad_t)
            CALL torch_tensor_release(kin_t)
            CALL torch_tensor_release(grid_coords_t)
            CALL torch_tensor_release(grid_weights_t)
            CALL torch_tensor_release(atomic_grid_weights_t)
            CALL torch_tensor_release(atomic_coords_t)
         END IF
         CALL torch_tensor_release(atomic_grid_sizes_t)
         NULLIFY (chunk_atom_begin, chunk_atom_count, chunk_row_begin, chunk_row_counts)
         CALL timestop(phase_handle)
      END IF

   END SUBROUTINE skala_gapw_atom_composite_energy

! **************************************************************************************************
!> \brief Select an automatic row cap for atom-composite Torch calls.
!> \param atomic_grid_sizes ...
!> \param group ...
!> \return ...
! **************************************************************************************************
   FUNCTION auto_atom_composite_chunk_max_rows(atomic_grid_sizes, group) RESULT(max_rows)
      INTEGER(KIND=int_8), DIMENSION(:), INTENT(IN)      :: atomic_grid_sizes

      CLASS(mp_comm_type), INTENT(IN)                    :: group

      INTEGER                                            :: local_max_atom_rows, local_natom, &
                                                            local_padded_rows, max_rows

      local_padded_rows = 0
      local_natom = SIZE(atomic_grid_sizes)
      IF (local_natom > 0) THEN
         local_max_atom_rows = INT(MAXVAL(atomic_grid_sizes))
         IF (local_natom > atom_chunk_auto_max_rows .OR. &
             local_max_atom_rows > atom_chunk_auto_max_rows/local_natom) THEN
            local_padded_rows = atom_chunk_auto_max_rows + 1
         ELSE
            local_padded_rows = local_natom*local_max_atom_rows
         END IF
      END IF
      CALL group%max(local_padded_rows)
      IF (local_padded_rows <= atom_chunk_auto_max_rows) THEN
         max_rows = 0
         RETURN
      END IF
      max_rows = atom_chunk_auto_max_rows

   END FUNCTION auto_atom_composite_chunk_max_rows

! **************************************************************************************************
!> \brief Rebuild the atom-composite chunk layout only when the grid sizes change.
!> \param atomic_grid_sizes ...
!> \param max_rows ...
! **************************************************************************************************
   SUBROUTINE ensure_atom_chunk_layout(atomic_grid_sizes, max_rows)
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: atomic_grid_sizes
      INTEGER, INTENT(IN)                                :: max_rows

      LOGICAL                                            :: valid

      valid = skala_state%chunk_layout%max_rows == max_rows
      IF (valid) THEN
         valid = ALLOCATED(skala_state%chunk_layout%grid_sizes)
         IF (valid) valid = SIZE(skala_state%chunk_layout%grid_sizes) == SIZE(atomic_grid_sizes)
         IF (valid) valid = ALL(skala_state%chunk_layout%grid_sizes == atomic_grid_sizes)
      END IF
      IF (valid) RETURN

      IF (ALLOCATED(skala_state%chunk_layout%atom_begin)) THEN
         DEALLOCATE (skala_state%chunk_layout%atom_begin, skala_state%chunk_layout%atom_count, &
                     skala_state%chunk_layout%row_begin, skala_state%chunk_layout%row_counts)
      END IF
      CALL atom_composite_chunk_layout(atomic_grid_sizes, max_rows, &
                                       skala_state%chunk_layout%atom_begin, &
                                       skala_state%chunk_layout%atom_count, &
                                       skala_state%chunk_layout%row_begin, &
                                       skala_state%chunk_layout%row_counts)
      IF (ALLOCATED(skala_state%chunk_layout%grid_sizes)) THEN
         IF (SIZE(skala_state%chunk_layout%grid_sizes) /= SIZE(atomic_grid_sizes)) THEN
            DEALLOCATE (skala_state%chunk_layout%grid_sizes)
         END IF
      END IF
      IF (.NOT. ALLOCATED(skala_state%chunk_layout%grid_sizes)) THEN
         ALLOCATE (skala_state%chunk_layout%grid_sizes(SIZE(atomic_grid_sizes)))
      END IF
      skala_state%chunk_layout%grid_sizes(:) = atomic_grid_sizes
      skala_state%chunk_layout%max_rows = max_rows

   END SUBROUTINE ensure_atom_chunk_layout

! **************************************************************************************************
!> \brief Split contiguous atom-composite rows without dividing an atomic block.
!> \param atomic_grid_sizes ...
!> \param max_rows ...
!> \param atom_begin ...
!> \param atom_count ...
!> \param row_begin ...
!> \param row_count ...
! **************************************************************************************************
   SUBROUTINE atom_composite_chunk_layout(atomic_grid_sizes, max_rows, atom_begin, atom_count, &
                                          row_begin, row_count)
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: atomic_grid_sizes
      INTEGER, INTENT(IN)                                :: max_rows
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: atom_begin, atom_count, row_begin, &
                                                            row_count

      INTEGER                                            :: atom_rows, chunk_atoms, &
                                                            chunk_max_atom_rows, iatom, natom, &
                                                            nchunks, rows, subchunk
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: padded_rows

      natom = SIZE(atomic_grid_sizes)
      CPASSERT(natom > 0)
      CPASSERT(ALL(atomic_grid_sizes > 0_int_8))
      IF (max_rows <= 0) THEN
         ALLOCATE (atom_begin(1), atom_count(1), row_begin(1), row_count(1))
         atom_begin(1) = 1
         atom_count(1) = natom
         row_begin(1) = 1
         row_count(1) = INT(SUM(atomic_grid_sizes))
         RETURN
      END IF

      nchunks = 1
      chunk_atoms = 0
      chunk_max_atom_rows = 0
      DO iatom = 1, natom
         atom_rows = INT(atomic_grid_sizes(iatom))
         IF (chunk_atoms > 0 .AND. &
             MAX(chunk_max_atom_rows, atom_rows) > max_rows/(chunk_atoms + 1)) THEN
            nchunks = nchunks + 1
            chunk_atoms = 0
            chunk_max_atom_rows = 0
         END IF
         chunk_atoms = chunk_atoms + 1
         chunk_max_atom_rows = MAX(chunk_max_atom_rows, atom_rows)
      END DO
      ALLOCATE (atom_begin(nchunks), atom_count(nchunks), row_begin(nchunks), row_count(nchunks))

      subchunk = 1
      atom_begin(subchunk) = 1
      row_begin(subchunk) = 1
      chunk_atoms = 0
      chunk_max_atom_rows = 0
      rows = 0
      DO iatom = 1, natom
         atom_rows = INT(atomic_grid_sizes(iatom))
         IF (chunk_atoms > 0 .AND. &
             MAX(chunk_max_atom_rows, atom_rows) > max_rows/(chunk_atoms + 1)) THEN
            atom_count(subchunk) = iatom - atom_begin(subchunk)
            row_count(subchunk) = rows
            subchunk = subchunk + 1
            atom_begin(subchunk) = iatom
            row_begin(subchunk) = row_begin(subchunk - 1) + row_count(subchunk - 1)
            chunk_atoms = 0
            chunk_max_atom_rows = 0
            rows = 0
         END IF
         chunk_atoms = chunk_atoms + 1
         chunk_max_atom_rows = MAX(chunk_max_atom_rows, atom_rows)
         rows = rows + atom_rows
      END DO
      atom_count(subchunk) = natom - atom_begin(subchunk) + 1
      row_count(subchunk) = rows
      CPASSERT(subchunk == nchunks)

      ALLOCATE (padded_rows(nchunks))
      DO subchunk = 1, nchunks
         iatom = atom_begin(subchunk) + atom_count(subchunk) - 1
         padded_rows(subchunk) = INT(atom_count(subchunk), KIND=int_8)* &
                                 INT(MAXVAL(atomic_grid_sizes(atom_begin(subchunk):iatom)), &
                                     KIND=int_8)
      END DO

      CALL sort_chunks_by_padded_rows(padded_rows, atom_begin, atom_count, row_begin, row_count)
      DEALLOCATE (padded_rows)

   END SUBROUTINE atom_composite_chunk_layout

! **************************************************************************************************
!> \brief Add the explicit SKALA derivative with respect to atom-center coordinates.
!> \param atom_force ...
!> \param root_rank ...
!> \param atom_coord_grad_in ...
! **************************************************************************************************
   SUBROUTINE add_explicit_coordinate_force(atom_force, root_rank, atom_coord_grad_in)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: atom_force
      LOGICAL, INTENT(IN)                                :: root_rank
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: atom_coord_grad_in

      IF (root_rank) THEN
         CPASSERT(SIZE(atom_force, 1) == SIZE(atom_coord_grad_in, 1))
         CPASSERT(SIZE(atom_force, 2) == SIZE(atom_coord_grad_in, 2))
         atom_force(:, :) = atom_force(:, :) + atom_coord_grad_in(:, :)
      END IF

   END SUBROUTINE add_explicit_coordinate_force

! **************************************************************************************************
!> \brief Add the force and/or virial from SMOOTH native-grid atom partition weights.
!>        The logarithmic weight derivatives are contracted directly, so no
!>        O(natom^2) derivative tensor is formed.
!> \param atom_force ...
!> \param virial_xc ...
!> \param features ...
!> \param particle_set ...
!> \param cell ...
!> \param rho_r ...
!> \param grid_weight_grad_in ...
!> \param atomic_grid_weight_grad_in ...
! **************************************************************************************************
   SUBROUTINE add_smooth_partition_force_virial(atom_force, virial_xc, features, particle_set, cell, &
                                                rho_r, grid_weight_grad_in, atomic_grid_weight_grad_in)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL                                        :: atom_force
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: virial_xc
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), TARGET    :: grid_weight_grad_in, &
                                                            atomic_grid_weight_grad_in

      INTEGER                                            :: feature_begin, feature_end, feature_pos, &
                                                            i, iatom, j, k, local_row, natom, row
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: do_force, do_virial
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: included
      REAL(KIND=dp)                                      :: grid_base_weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: distances, weight_grad, weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_coords_pbc, partition_atom_coords
      REAL(KIND=dp), DIMENSION(3)                        :: grid_point
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atomic_grid_weight_grad, grid_weight_grad

      do_force = PRESENT(atom_force)
      do_virial = PRESENT(virial_xc)
      CPASSERT(do_force .OR. do_virial)

      grid_weight_grad => grid_weight_grad_in
      atomic_grid_weight_grad => atomic_grid_weight_grad_in

      natom = SIZE(particle_set)
      IF (do_force) THEN
         CPASSERT(SIZE(atom_force, 1) == 3)
         CPASSERT(SIZE(atom_force, 2) == natom)
      END IF
      ALLOCATE (atom_coords_pbc(3, natom), included(natom), weights(natom), weight_grad(natom), &
                distances(natom), partition_atom_coords(3, natom))
      DO iatom = 1, natom
         atom_coords_pbc(:, iatom) = pbc(particle_set(iatom)%r, cell, positive_range=.TRUE.)
      END DO

      bo = rho_r(1)%pw_grid%bounds_local
      local_row = 0
      DO k = bo(1, 3), bo(2, 3)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 1), bo(2, 1)
               local_row = local_row + 1
               grid_point = grid_coordinate(rho_r(1)%pw_grid, [i, j, k])
               CALL smooth_partition_weights(grid_point, atom_coords_pbc, cell, weights, included, &
                                             partition_atom_coords_out=partition_atom_coords, &
                                             distances_out=distances)
               feature_begin = features%local_feature_offsets(local_row)
               feature_end = features%local_feature_offsets(local_row + 1) - 1
               CPASSERT(feature_end - feature_begin + 1 == COUNT(included))
               grid_base_weight = 0.0_dp
               DO feature_pos = feature_begin, feature_end
                  row = features%local_feature_rows(feature_pos)
                  grid_base_weight = grid_base_weight + features%grid_weights(row)
               END DO
               weight_grad = 0.0_dp
               feature_pos = feature_begin
               DO iatom = 1, natom
                  IF (.NOT. included(iatom)) CYCLE
                  row = features%local_feature_rows(feature_pos)
                  weight_grad(iatom) = grid_base_weight*grid_weight_grad(row) + &
                                       rho_r(1)%pw_grid%dvol*atomic_grid_weight_grad(row)* &
                                       smooth_partition_atomic_weight_scale_derivative( &
                                       weights(iatom))
                  feature_pos = feature_pos + 1
               END DO
               CPASSERT(feature_pos == feature_end + 1)

               IF (do_force .AND. do_virial) THEN
                  CALL smooth_partition_derivative_action( &
                     grid_point, atom_coords_pbc, cell, weights, included, weight_grad, &
                     atom_force=atom_force, virial_xc=virial_xc, &
                     partition_atom_coords_in=partition_atom_coords, distances_in=distances)
               ELSE IF (do_force) THEN
                  CALL smooth_partition_derivative_action( &
                     grid_point, atom_coords_pbc, cell, weights, included, weight_grad, &
                     atom_force=atom_force, partition_atom_coords_in=partition_atom_coords, &
                     distances_in=distances)
               ELSE
                  CALL smooth_partition_derivative_action( &
                     grid_point, atom_coords_pbc, cell, weights, included, weight_grad, &
                     virial_xc=virial_xc, partition_atom_coords_in=partition_atom_coords, &
                     distances_in=distances)
               END IF
            END DO
         END DO
      END DO
      CPASSERT(local_row == features%nflat_local)

      DEALLOCATE (atom_coords_pbc, included, weights, weight_grad, distances, partition_atom_coords)

   END SUBROUTINE add_smooth_partition_force_virial

! **************************************************************************************************
!> \brief Return the Cartesian coordinate of a regular GPW grid point.
!> \brief Prepare scaled alltoall counts/displacements and exchange buffers.
!> \param send_src_counts ...
!> \param send_src_displs ...
!> \param recv_src_counts ...
!> \param recv_src_displs ...
!> \param factor ...
!> \param send_base_size unscaled send element count
!> \param recv_base_size unscaled receive element count
!> \param send_counts ...
!> \param send_displs ...
!> \param recv_counts ...
!> \param recv_displs ...
!> \param send_buffer ...
!> \param recv_buffer ...
!> \param zero_buffers ...
!> \return ...
! **************************************************************************************************
   SUBROUTINE prepare_route_exchange(send_src_counts, send_src_displs, recv_src_counts, &
                                     recv_src_displs, factor, send_base_size, recv_base_size, &
                                     send_counts, send_displs, recv_counts, recv_displs, &
                                     send_buffer, recv_buffer, zero_buffers)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: send_src_counts, send_src_displs, &
                                                            recv_src_counts, recv_src_displs
      INTEGER, INTENT(IN)                                :: factor, send_base_size, recv_base_size
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: send_counts, send_displs, recv_counts, &
                                                            recv_displs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: send_buffer, recv_buffer
      LOGICAL, INTENT(IN), OPTIONAL                      :: zero_buffers

      LOGICAL                                            :: zero

      zero = .FALSE.
      IF (PRESENT(zero_buffers)) zero = zero_buffers

      ALLOCATE (send_counts(SIZE(send_src_counts)), send_displs(SIZE(send_src_displs)), &
                recv_counts(SIZE(recv_src_counts)), recv_displs(SIZE(recv_src_displs)))
      send_counts(:) = factor*send_src_counts
      send_displs(:) = factor*send_src_displs
      recv_counts(:) = factor*recv_src_counts
      recv_displs(:) = factor*recv_src_displs

      ALLOCATE (send_buffer(MAX(1, factor*send_base_size)), &
                recv_buffer(MAX(1, factor*recv_base_size)))
      IF (zero) THEN
         send_buffer = 0.0_dp
         recv_buffer = 0.0_dp
      END IF

   END SUBROUTINE prepare_route_exchange

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
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: density_grad, kin_grad
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: grad_grad

      INTEGER                                            :: isubchunk, nroute_grad_per_point, &
                                                            nroute_recv_points, nroute_send_points, &
                                                            nsubchunks, phase_handle, subphase_handle
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: route_grad_return_recv_counts, &
                                                            route_grad_return_recv_displs, &
                                                            route_grad_return_send_counts, &
                                                            route_grad_return_send_displs, &
                                                            subchunk_atom_begin, &
                                                            subchunk_atom_count, &
                                                            subchunk_row_begin, &
                                                            subchunk_row_count
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: recv_grad_buffer, send_grad_buffer
      TYPE(skala_gpw_feature_type)                       :: subchunk
      TYPE(torch_tensor_type)                            :: subchunk_exc_tensor, subchunk_exc_total_t

      CPASSERT(features%uses_atom_chunks)
      CPASSERT(max_rows > 0)
      CALL skala_gpw_atom_subchunk_layout(max_rows, subchunk_atom_begin, subchunk_atom_count, &
                                          subchunk_row_begin, subchunk_row_count)
      nsubchunks = SIZE(subchunk_atom_begin)

      exc = 0.0_dp
      IF (compute_grads) THEN
         nroute_recv_points = SUM(features%route_point_recv_counts)
         nroute_send_points = SIZE(features%route_send_local_rows)
         CPASSERT(SUM(features%route_point_send_counts) == nroute_send_points)
         CPASSERT(SIZE(features%route_chunk_offsets) == nroute_recv_points + 1)
         nroute_grad_per_point = ngrad_per_point
         IF (collapse_spin_grads) nroute_grad_per_point = ncollapsed_grad_per_point
         CALL prepare_route_exchange( &
            features%route_point_recv_counts, features%route_point_recv_displs, &
            features%route_point_send_counts, features%route_point_send_displs, &
            nroute_grad_per_point, nroute_recv_points, nroute_send_points, &
            route_grad_return_send_counts, route_grad_return_send_displs, &
            route_grad_return_recv_counts, route_grad_return_recv_displs, &
            send_grad_buffer, recv_grad_buffer)
      END IF

      CALL timeset("skala_gpw_atom_subchunks", phase_handle)
      DO isubchunk = 1, nsubchunks
         CALL timeset("skala_gpw_atom_subchunk_build", subphase_handle)
         CALL skala_gpw_feature_build_atom_subchunk_bounds(features, subchunk, &
                                                           subchunk_atom_begin(isubchunk), &
                                                           subchunk_atom_count(isubchunk), &
                                                           subchunk_row_begin(isubchunk), &
                                                           subchunk_row_count(isubchunk))
         CALL timestop(subphase_handle)
         CALL timeset("skala_gpw_atom_subchunk_forward", subphase_handle)
         CALL skala_torch_model_get_exc_device(skala_state%model, subchunk%inputs, &
                                               subchunk%grid_weights_t, subchunk_exc_tensor)
         CALL torch_tensor_accumulate(subchunk_exc_total_t, subchunk_exc_tensor)
         CALL timestop(subphase_handle)
         IF (compute_grads) THEN
            CALL timeset("skala_gpw_atom_subchunk_backward", subphase_handle)
            CALL torch_tensor_backward_scalar(subchunk_exc_tensor)
            CALL timestop(subphase_handle)
         END IF
         CALL timeset("skala_gpw_atom_subchunk_release", subphase_handle)
         CALL torch_tensor_release(subchunk_exc_tensor)
         CALL timestop(subphase_handle)
      END DO
      IF (subchunk%active) CALL skala_gpw_feature_release(subchunk)
      IF (nsubchunks > 0) THEN
         exc = exc + torch_tensor_item_double(subchunk_exc_total_t)
         CALL torch_tensor_release(subchunk_exc_total_t)
      END IF
      IF (compute_grads .AND. features%chunk_feature_count > 0) THEN
         CALL timeset("skala_gpw_atom_subchunk_grad_pack", subphase_handle)
         CALL pack_atom_chunk_grads(features, send_grad_buffer, collapse_spin_grads)
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
         CALL scatter_routed_atom_chunk_grads(features, recv_grad_buffer, collapse_spin_grads, &
                                              density_grad, grad_grad, kin_grad)
         CALL timestop(phase_handle)

         DEALLOCATE (recv_grad_buffer, route_grad_return_recv_counts, &
                     route_grad_return_recv_displs, route_grad_return_send_counts, &
                     route_grad_return_send_displs, send_grad_buffer)
      END IF
      DEALLOCATE (subchunk_atom_begin, subchunk_atom_count, subchunk_row_begin, subchunk_row_count)

   END SUBROUTINE evaluate_atom_subchunks

! **************************************************************************************************
!> \brief Select an automatic atom-subchunk row cap.
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
!> \brief Compute the density/tau change-detection signature in a single pass per array,
!>        without materializing full-grid temporaries.
!> \param rho_r ...
!> \param tau ...
!> \param rho_sum ...
!> \param rho_sumsq ...
!> \param tau_sum ...
!> \param tau_sumsq ...
! **************************************************************************************************
   SUBROUTINE weight_grad_signature(rho_r, tau, rho_sum, rho_sumsq, tau_sum, tau_sumsq)
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r, tau
      REAL(KIND=dp), INTENT(OUT)                         :: rho_sum, rho_sumsq, tau_sum, tau_sumsq

      INTEGER                                            :: i, ispin, j, k
      REAL(KIND=dp)                                      :: x, y

      rho_sum = 0.0_dp
      rho_sumsq = 0.0_dp
      tau_sum = 0.0_dp
      tau_sumsq = 0.0_dp
      DO ispin = 1, SIZE(rho_r)
         DO k = 1, SIZE(rho_r(ispin)%array, 3)
            DO j = 1, SIZE(rho_r(ispin)%array, 2)
               DO i = 1, SIZE(rho_r(ispin)%array, 1)
                  x = rho_r(ispin)%array(i, j, k)
                  rho_sum = rho_sum + x
                  rho_sumsq = rho_sumsq + x*x
               END DO
            END DO
         END DO
         DO k = 1, SIZE(tau(ispin)%array, 3)
            DO j = 1, SIZE(tau(ispin)%array, 2)
               DO i = 1, SIZE(tau(ispin)%array, 1)
                  y = tau(ispin)%array(i, j, k)
                  tau_sum = tau_sum + y
                  tau_sumsq = tau_sumsq + y*y
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE weight_grad_signature

! **************************************************************************************************
!> \brief Store the smooth-partition weight gradient for later reuse.
!> \param gradient ...
!> \param features ...
!> \param particle_set ...
!> \param cell ...
!> \param rho_r ...
!> \param tau ...
! **************************************************************************************************
   SUBROUTINE store_weight_grad_cache(gradient, features, particle_set, cell, rho_r, tau)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: gradient
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r, tau

      INTEGER                                            :: iatom

      IF (ASSOCIATED(skala_state%weight_grad%gradient)) THEN
         IF (SIZE(skala_state%weight_grad%gradient) /= SIZE(gradient)) THEN
            DEALLOCATE (skala_state%weight_grad%gradient)
         END IF
      END IF
      IF (.NOT. ASSOCIATED(skala_state%weight_grad%gradient)) THEN
         ALLOCATE (skala_state%weight_grad%gradient(SIZE(gradient)))
      END IF
      skala_state%weight_grad%gradient = gradient

      IF (ALLOCATED(skala_state%weight_grad%positions)) THEN
         IF (SIZE(skala_state%weight_grad%positions, 2) /= SIZE(particle_set)) THEN
            DEALLOCATE (skala_state%weight_grad%positions)
         END IF
      END IF
      IF (.NOT. ALLOCATED(skala_state%weight_grad%positions)) THEN
         ALLOCATE (skala_state%weight_grad%positions(3, SIZE(particle_set)))
      END IF
      DO iatom = 1, SIZE(particle_set)
         skala_state%weight_grad%positions(:, iatom) = particle_set(iatom)%r
      END DO

      CALL weight_grad_signature(rho_r, tau, skala_state%weight_grad%rho_sum, &
                                 skala_state%weight_grad%rho_sumsq, &
                                 skala_state%weight_grad%tau_sum, skala_state%weight_grad%tau_sumsq)
      skala_state%weight_grad%nvals = SIZE(gradient)
      skala_state%weight_grad%nspins = SIZE(rho_r)
      skala_state%weight_grad%atom_chunks = features%uses_atom_chunks
      skala_state%weight_grad%atom_partition = features%atom_partition
      skala_state%weight_grad%cell_hmat = cell%hmat
      skala_state%weight_grad%cell_perd = cell%perd
      skala_state%weight_grad%valid = .TRUE.

   END SUBROUTINE store_weight_grad_cache

! **************************************************************************************************
!> \brief Whether the cached weight gradient matches the requested evaluation.
!>
!> The reuse decision must be identical on every MPI rank: a hit skips the
!> model backward and the atom-chunk gradient exchange, so a rank-local hit
!> while another rank misses would deadlock. The local test is therefore
!> reduced with a logical AND across the communicator.
!> \param group ...
!> \param particle_set ...
!> \param cell ...
!> \param rho_r ...
!> \param tau ...
!> \param nspins ...
!> \param atom_chunks ...
!> \param atom_partition ...
!> \param nvals ...
!> \return ...
! **************************************************************************************************
   FUNCTION weight_grad_cache_matches(group, particle_set, cell, rho_r, tau, nspins, &
                                      atom_chunks, atom_partition, nvals) RESULT(matches)
      CLASS(mp_comm_type), INTENT(IN)                    :: group
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r, tau
      INTEGER, INTENT(IN)                                :: nspins, atom_partition, nvals
      LOGICAL, INTENT(IN)                                :: atom_chunks
      LOGICAL                                            :: matches
      INTEGER                                            :: hit

      hit = MERGE(1, 0, weight_grad_cache_matches_local(particle_set, cell, rho_r, tau, &
                                                        nspins, atom_chunks, atom_partition, nvals))
      CALL group%min(hit)
      matches = (hit == 1)

   END FUNCTION weight_grad_cache_matches

! **************************************************************************************************
!> \brief Rank-local part of weight_grad_cache_matches (no collectives).
!> \param particle_set ...
!> \param cell ...
!> \param rho_r ...
!> \param tau ...
!> \param nspins ...
!> \param atom_chunks ...
!> \param atom_partition ...
!> \param nvals ...
!> \return ...
! **************************************************************************************************
   FUNCTION weight_grad_cache_matches_local(particle_set, cell, rho_r, tau, nspins, &
                                            atom_chunks, atom_partition, nvals) RESULT(matches)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r, tau
      INTEGER, INTENT(IN)                                :: nspins
      LOGICAL, INTENT(IN)                                :: atom_chunks
      INTEGER, INTENT(IN)                                :: atom_partition, nvals
      LOGICAL                                            :: matches

      INTEGER                                            :: iatom
      REAL(KIND=dp)                                      :: rho_sum, rho_sumsq, tau_sum, tau_sumsq

      matches = .FALSE.
      IF (nvals <= 0) RETURN
      IF (.NOT. skala_state%weight_grad%valid) RETURN
      IF (skala_state%weight_grad%nspins /= nspins) RETURN
      IF (skala_state%weight_grad%atom_chunks .NEQV. atom_chunks) RETURN
      IF (skala_state%weight_grad%atom_partition /= atom_partition) RETURN
      IF (skala_state%weight_grad%nvals /= nvals) RETURN
      IF (ANY(ABS(skala_state%weight_grad%cell_hmat - cell%hmat) > weight_cache_cell_tol)) RETURN
      IF (ANY(skala_state%weight_grad%cell_perd /= cell%perd)) RETURN
      IF (.NOT. ALLOCATED(skala_state%weight_grad%positions)) RETURN
      IF (SIZE(skala_state%weight_grad%positions, 2) /= SIZE(particle_set)) RETURN
      DO iatom = 1, SIZE(particle_set)
         IF (ANY(skala_state%weight_grad%positions(:, iatom) /= particle_set(iatom)%r)) RETURN
      END DO

      CALL weight_grad_signature(rho_r, tau, rho_sum, rho_sumsq, tau_sum, tau_sumsq)
      IF (rho_sum /= skala_state%weight_grad%rho_sum) RETURN
      IF (rho_sumsq /= skala_state%weight_grad%rho_sumsq) RETURN
      IF (tau_sum /= skala_state%weight_grad%tau_sum) RETURN
      IF (tau_sumsq /= skala_state%weight_grad%tau_sumsq) RETURN
      matches = .TRUE.

   END FUNCTION weight_grad_cache_matches_local

! **************************************************************************************************
!> \brief Ensure the reusable host gradient scratch has the local-grid shape.
!> \param nflat_local ...
! **************************************************************************************************
   SUBROUTINE ensure_grad_scratch(nflat_local)
      INTEGER, INTENT(IN)                                :: nflat_local

      IF (ASSOCIATED(skala_state%grad_scratch%density_grad)) THEN
         IF (SIZE(skala_state%grad_scratch%density_grad, 1) == nflat_local .AND. &
             SIZE(skala_state%grad_scratch%grad_grad, 1) == nflat_local .AND. &
             SIZE(skala_state%grad_scratch%kin_grad, 1) == nflat_local) RETURN
         DEALLOCATE (skala_state%grad_scratch%density_grad, skala_state%grad_scratch%grad_grad, &
                     skala_state%grad_scratch%kin_grad)
      END IF
      ALLOCATE (skala_state%grad_scratch%density_grad(nflat_local, 2), &
                skala_state%grad_scratch%grad_grad(nflat_local, 3, 2), &
                skala_state%grad_scratch%kin_grad(nflat_local, 2))
      skala_state%grad_scratch%density_grad = 0.0_dp
      skala_state%grad_scratch%grad_grad = 0.0_dp
      skala_state%grad_scratch%kin_grad = 0.0_dp

   END SUBROUTINE ensure_grad_scratch

! **************************************************************************************************
!> \brief Map full Torch feature gradients back to this rank's local grid order.
!> \param features ...
!> \param density_grad ...
!> \param grad_grad ...
!> \param kin_grad ...
! **************************************************************************************************
   SUBROUTINE fetch_local_feature_grads(features, density_grad, grad_grad, kin_grad)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: density_grad
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: grad_grad
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: kin_grad

      INTEGER                                            :: feature_pos, local_row, nrows_all, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: density_grad_all, kin_grad_all
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grad_grad_all
      TYPE(torch_tensor_type)                            :: density_grad_t, grad_grad_t, kin_grad_t

      NULLIFY (density_grad_all, grad_grad_all, kin_grad_all)
      nrows_all = features%nflat
      CALL get_feature_grad_views(features, density_grad_t, grad_grad_t, kin_grad_t, &
                                  density_grad_all, grad_grad_all, kin_grad_all)
      CPASSERT(SIZE(density_grad_all, 1) == nrows_all)
      CPASSERT(SIZE(density_grad_all, 2) == 2)
      CPASSERT(SIZE(grad_grad_all, 1) == nrows_all)
      CPASSERT(SIZE(grad_grad_all, 2) == 3)
      CPASSERT(SIZE(grad_grad_all, 3) == 2)
      CPASSERT(SIZE(kin_grad_all, 1) == nrows_all)
      CPASSERT(SIZE(kin_grad_all, 2) == 2)

      density_grad = 0.0_dp
      grad_grad = 0.0_dp
      kin_grad = 0.0_dp
      DO local_row = 1, features%nflat_local
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

      CALL torch_tensor_release(density_grad_t)
      CALL torch_tensor_release(grad_grad_t)
      CALL torch_tensor_release(kin_grad_t)

   END SUBROUTINE fetch_local_feature_grads

! **************************************************************************************************
!> \brief Cache the unmodified reference fields for the nonlocal Skala response.
!> \param rho_set Reference density, gradient and kinetic-energy density.
!> \param rho_r Real-space reference density.
!> \param tau_r Reference kinetic-energy density.
!> \param xc_section XC settings.
!> \param pw_pool Grid and Fourier transforms.
! **************************************************************************************************
   SUBROUTINE skala_gpw_fxc_prep(rho_set, rho_r, tau_r, xc_section, pw_pool)
      TYPE(xc_rho_set_type), INTENT(OUT)                 :: rho_set
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r, tau_r
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      INTEGER                                            :: ispin
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: tau_copy
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(xc_rho_cflags_type)                           :: needs

      IF (.NOT. ASSOCIATED(tau_r)) THEN
         CALL cp_abort(__LOCATION__, "Native Skala response requires kinetic-energy density.")
      END IF
      IF (SIZE(rho_r) == 1) THEN
         needs%rho = .TRUE.
         needs%drho = .TRUE.
         needs%tau = .TRUE.
      ELSE
         needs%rho_spin = .TRUE.
         needs%drho_spin = .TRUE.
         needs%tau_spin = .TRUE.
      END IF
      CALL xc_rho_set_create(rho_set, rho_r(1)%pw_grid%bounds_local, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))
      NULLIFY (rho_g)
      CALL xc_rho_set_update(rho_set, rho_r, rho_g, tau_r, needs, &
                             section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                             section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), pw_pool)
      DO ispin = 1, SIZE(tau_r)
         NULLIFY (tau_copy)
         CALL pw_pool%create_cr3d(tau_copy)
         tau_copy = tau_r(ispin)%array
         IF (SIZE(tau_r) == 1) THEN
            rho_set%tau => tau_copy
            rho_set%owns%tau = .TRUE.
         ELSE IF (ispin == 1) THEN
            rho_set%tau_a => tau_copy
         ELSE
            rho_set%tau_b => tau_copy
            rho_set%owns%tau_spin = .TRUE.
         END IF
      END DO
   END SUBROUTINE skala_gpw_fxc_prep

! **************************************************************************************************
!> \brief Apply the nonlocal Skala Hessian to a density and kinetic-density response.
!> \param rho0_set Reference fields.
!> \param rho1_r Density response.
!> \param tau1_r Kinetic-energy-density response.
!> \param xc_section XC settings.
!> \param weights Integration weights.
!> \param pw_pool Grid and Fourier transforms.
!> \param particle_set Fixed nuclei.
!> \param cell Fixed simulation cell.
!> \param fxc_rho XC potential response.
!> \param fxc_tau Kinetic-energy-density potential response.
! **************************************************************************************************
   SUBROUTINE skala_gpw_fxc(rho0_set, rho1_r, tau1_r, xc_section, weights, pw_pool, &
                            particle_set, cell, fxc_rho, fxc_tau)
      TYPE(xc_rho_set_type), INTENT(INOUT)               :: rho0_set
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho1_r, tau1_r
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rho, fxc_tau

      INTEGER :: i, isub, max_rows, nroute_grad_per_point, nroute_recv_points, nroute_send_points, &
         nsubchunks, partition, row0, selected_device, spin_dim
      INTEGER, ALLOCATABLE, DIMENSION(:) :: route_grad_return_recv_counts, &
         route_grad_return_recv_displs, route_grad_return_send_counts, &
         route_grad_return_send_displs, sub_atom_begin, sub_atom_count, sub_row_begin, &
         sub_row_count
      LOGICAL                                            :: collapse_spin_grads, fp64, use_cuda
      REAL(KIND=dp), ALLOCATABLE :: density_grad(:, :), grad_grad(:, :, :), kin_grad(:, :), &
         recv_grad_buffer(:), resp_density(:, :), resp_grad(:, :, :), resp_kin(:, :), &
         send_grad_buffer(:)
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sub_density, sub_kin
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: sub_grad
      TYPE(section_vals_type), POINTER                   :: native_section
      TYPE(skala_gpw_feature_type)                       :: dir_sub, direction, ref_sub, reference
      TYPE(torch_tensor_type)                            :: directions(3), exc_sub_t, sub_resp(3)
      TYPE(xc_rho_set_type)                              :: rho1_set

      fp64 = .TRUE.
      partition = skala_gpw_atom_partition_smooth
      native_section => get_skala_native_section(xc_section)
      IF (ASSOCIATED(native_section)) THEN
         CALL section_vals_val_get(native_section, "MODEL_FP64", l_val=fp64)
         CALL section_vals_val_get(native_section, "ATOM_PARTITION", i_val=partition)
      END IF
      IF (.NOT. fp64) THEN
         CALL cp_abort(__LOCATION__, "Native Skala response requires MODEL_FP64.")
      END IF
      CALL ensure_model_context(xc_section, rho1_r(1)%pw_grid%para%group, &
                                use_cuda, selected_device, higher_derivatives=.TRUE.)
      CALL torch_use_cuda(use_cuda)
      ! Response uses the atom-chunk layout. Release the cached layout only when
      ! it currently holds the full-grid variant (to bound memory and avoid mixing
      ! modes); otherwise reuse it so repeated response vectors do not rebuild the
      ! atom-chunk layout and static tensors on every call.
      IF (skala_state%layout_mode_known .AND. skala_state%layout_mode_full) THEN
         CALL skala_gpw_features_release_cache()
      END IF

      CALL skala_gpw_fxc_prep(rho1_set, rho1_r, tau1_r, xc_section, pw_pool)
      CALL skala_gpw_feature_build(direction, rho1_set, rho1_r, particle_set, cell, &
                                   requires_grad=.FALSE., weights=weights, atom_partition=partition, &
                                   use_atom_chunks=.TRUE., cuda_device=selected_device)
      IF (direction%chunk_feature_count > 0) THEN
         CALL torch_tensor_clone(direction%density_t, directions(1))
         CALL torch_tensor_clone(direction%grad_t, directions(2))
         CALL torch_tensor_clone(direction%kin_t, directions(3))
      END IF
      CALL skala_gpw_feature_build(reference, rho0_set, rho1_r, particle_set, cell, &
                                   requires_grad=.TRUE., weights=weights, atom_partition=partition, &
                                   use_atom_chunks=.TRUE., cuda_device=selected_device)
      collapse_spin_grads = reference%uses_collapsed_rks_dynamic
      spin_dim = 2
      IF (collapse_spin_grads) spin_dim = 1
      ! The response is the Hessian-vector of a sum of per-atom contributions, so it
      ! is block-diagonal per chunk row. Split this rank's chunk into row subchunks to
      ! bound peak memory and accumulate each subchunk's feature gradients.
      max_rows = -1
      IF (ASSOCIATED(native_section)) THEN
         CALL section_vals_val_get(native_section, "ATOM_CHUNK_MAX_ROWS", i_val=max_rows)
      END IF
      IF (max_rows < -1) THEN
         CALL cp_abort(__LOCATION__, &
                       "SKALA_NATIVE%ATOM_CHUNK_MAX_ROWS must be -1, zero, or positive.")
      END IF
      nsubchunks = 0
      IF (reference%chunk_feature_count > 0) THEN
         IF (max_rows == -1) THEN
            max_rows = auto_atom_chunk_max_rows(reference, rho1_r(1)%pw_grid%para%group)
         END IF
         CALL skala_gpw_atom_subchunk_layout(max_rows, sub_atom_begin, sub_atom_count, &
                                             sub_row_begin, sub_row_count)
         nsubchunks = SIZE(sub_atom_begin)
         ! Route the direction subchunks through the independent clones.
         direction%density_t = directions(1)
         direction%grad_t = directions(2)
         direction%kin_t = directions(3)
         direction%owns_dynamic_tensors = .FALSE.
         ALLOCATE (resp_density(reference%chunk_feature_count, spin_dim), &
                   resp_grad(reference%chunk_feature_count, 3, spin_dim), &
                   resp_kin(reference%chunk_feature_count, spin_dim))
         resp_density = 0.0_dp
         resp_grad = 0.0_dp
         resp_kin = 0.0_dp
         DO isub = 1, nsubchunks
            CALL skala_gpw_feature_build_atom_subchunk_bounds(reference, ref_sub, &
                                                              sub_atom_begin(isub), sub_atom_count(isub), &
                                                              sub_row_begin(isub), sub_row_count(isub))
            CALL skala_gpw_feature_build_atom_subchunk_bounds(direction, dir_sub, &
                                                              sub_atom_begin(isub), sub_atom_count(isub), &
                                                              sub_row_begin(isub), sub_row_count(isub))
            CALL skala_torch_model_get_exc_device(skala_state%model, ref_sub%inputs, &
                                                  ref_sub%grid_weights_t, exc_sub_t)
            CALL torch_tensor_hessian_vector(exc_sub_t, &
                                             [ref_sub%density_t, ref_sub%grad_t, ref_sub%kin_t], &
                                             [dir_sub%density_t, dir_sub%grad_t, dir_sub%kin_t], &
                                             sub_resp, to_host=.TRUE.)
            NULLIFY (sub_density, sub_grad, sub_kin)
            CALL torch_tensor_data_ptr(sub_resp(1), sub_density)
            CALL torch_tensor_data_ptr(sub_resp(2), sub_grad)
            CALL torch_tensor_data_ptr(sub_resp(3), sub_kin)
            row0 = sub_row_begin(isub) - 1
            resp_density(row0 + 1:row0 + sub_row_count(isub), :) = sub_density
            resp_grad(row0 + 1:row0 + sub_row_count(isub), :, :) = sub_grad
            resp_kin(row0 + 1:row0 + sub_row_count(isub), :) = sub_kin
            DO i = 1, SIZE(sub_resp)
               CALL torch_tensor_release(sub_resp(i))
            END DO
            CALL torch_tensor_release(exc_sub_t)
            CALL skala_gpw_feature_release(ref_sub)
            CALL skala_gpw_feature_release(dir_sub)
         END DO
      END IF

      ! Route the chunk feature gradients back to the grid-point owners.
      nroute_recv_points = SUM(reference%route_point_recv_counts)
      nroute_send_points = SIZE(reference%route_send_local_rows)
      nroute_grad_per_point = ngrad_per_point
      IF (collapse_spin_grads) nroute_grad_per_point = ncollapsed_grad_per_point
      CALL prepare_route_exchange( &
         reference%route_point_recv_counts, reference%route_point_recv_displs, &
         reference%route_point_send_counts, reference%route_point_send_displs, &
         nroute_grad_per_point, nroute_recv_points, nroute_send_points, &
         route_grad_return_send_counts, route_grad_return_send_displs, &
         route_grad_return_recv_counts, route_grad_return_recv_displs, &
         send_grad_buffer, recv_grad_buffer)
      IF (reference%chunk_feature_count > 0) THEN
         CALL pack_atom_chunk_values(reference, send_grad_buffer, collapse_spin_grads, &
                                     resp_density, resp_grad, resp_kin)
      END IF
      CALL rho1_r(1)%pw_grid%para%group%alltoall(send_grad_buffer, &
                                                 route_grad_return_send_counts, route_grad_return_send_displs, recv_grad_buffer, &
                                                 route_grad_return_recv_counts, route_grad_return_recv_displs)
      ALLOCATE (density_grad(reference%nflat_local, 2), &
                grad_grad(reference%nflat_local, 3, 2), &
                kin_grad(reference%nflat_local, 2))
      CALL scatter_routed_atom_chunk_grads(reference, recv_grad_buffer, collapse_spin_grads, &
                                           density_grad, grad_grad, kin_grad)
      CALL build_vxc_from_feature_grads(fxc_rho, fxc_tau, rho1_r, pw_pool, &
                                        density_grad, grad_grad, kin_grad, &
                                        section_get_ival(xc_section, "XC_GRID%XC_DERIV"))

      DEALLOCATE (density_grad, grad_grad, kin_grad, send_grad_buffer, recv_grad_buffer, &
                  route_grad_return_send_counts, route_grad_return_send_displs, &
                  route_grad_return_recv_counts, route_grad_return_recv_displs)
      IF (reference%chunk_feature_count > 0) THEN
         DEALLOCATE (resp_density, resp_grad, resp_kin, &
                     sub_atom_begin, sub_atom_count, sub_row_begin, sub_row_count)
      END IF
      IF (direction%chunk_feature_count > 0) THEN
         DO i = 1, SIZE(directions)
            CALL torch_tensor_release(directions(i))
         END DO
      END IF
      CALL skala_gpw_feature_release(direction)
      CALL skala_gpw_feature_release(reference)
      CALL xc_rho_set_release(rho1_set, pw_pool=pw_pool)
      CALL torch_use_cuda(use_cuda)
   END SUBROUTINE skala_gpw_fxc

! **************************************************************************************************
!> \brief Pack atom-chunk Torch gradients into CP2K communication buffers.
!> \param features ...
!> \param TARGET ...
!> \param collapse_spin_grads ...
! **************************************************************************************************
   SUBROUTINE pack_atom_chunk_grads(features, TARGET, collapse_spin_grads)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: target
      LOGICAL, INTENT(IN), OPTIONAL                      :: collapse_spin_grads

      LOGICAL                                            :: my_collapse_spin_grads
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: chunk_density_grad, chunk_kin_grad
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: chunk_grad_grad
      TYPE(torch_tensor_type)                            :: density_grad_t, grad_grad_t, kin_grad_t

      my_collapse_spin_grads = .FALSE.
      IF (PRESENT(collapse_spin_grads)) my_collapse_spin_grads = collapse_spin_grads

      NULLIFY (chunk_density_grad, chunk_grad_grad, chunk_kin_grad)
      CALL get_feature_grad_views(features, density_grad_t, grad_grad_t, kin_grad_t, &
                                  chunk_density_grad, chunk_grad_grad, chunk_kin_grad)
      CPASSERT(SIZE(chunk_density_grad, 1) == features%chunk_feature_count)
      CPASSERT(SIZE(chunk_grad_grad, 1) == features%chunk_feature_count)
      CPASSERT(SIZE(chunk_grad_grad, 2) == 3)
      CPASSERT(SIZE(chunk_kin_grad, 1) == features%chunk_feature_count)
      CPASSERT(.NOT. my_collapse_spin_grads .OR. features%uses_collapsed_rks_dynamic)
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

      CALL pack_atom_chunk_values(features, TARGET, my_collapse_spin_grads, &
                                  chunk_density_grad, chunk_grad_grad, chunk_kin_grad)

      CALL torch_tensor_release(density_grad_t)
      CALL torch_tensor_release(grad_grad_t)
      CALL torch_tensor_release(kin_grad_t)

   END SUBROUTINE pack_atom_chunk_grads

! **************************************************************************************************
!> \brief Pack chunk-level Hessian feature gradients into the point-routed buffer.
!> \param features ...
!> \param TARGET ...
!> \param collapse_spin_grads ...
!> \param density_values ...
!> \param grad_values ...
!> \param kin_values ...
! **************************************************************************************************
   SUBROUTINE pack_atom_chunk_values(features, TARGET, collapse_spin_grads, &
                                     density_values, grad_values, kin_values)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: target
      LOGICAL, INTENT(IN)                                :: collapse_spin_grads
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: density_values
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: grad_values
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: kin_values

      INTEGER                                            :: base, feature_pos, irow, &
                                                            ngrad_buffer_per_point, point_pos, &
                                                            target_points

      ngrad_buffer_per_point = ngrad_per_point
      IF (collapse_spin_grads) ngrad_buffer_per_point = ncollapsed_grad_per_point
      CPASSERT(MOD(SIZE(TARGET), ngrad_buffer_per_point) == 0)
      target_points = SIZE(TARGET)/ngrad_buffer_per_point
      CPASSERT(target_points == SIZE(features%route_chunk_offsets) - 1)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(collapse_spin_grads, density_values, features, grad_values, kin_values, &
!$OMP        ngrad_buffer_per_point, TARGET, target_points) &
!$OMP PRIVATE(base, feature_pos, irow, point_pos)
      DO point_pos = 1, target_points
         base = ngrad_buffer_per_point*(point_pos - 1)
         TARGET(base + 1:base + ngrad_buffer_per_point) = 0.0_dp
         DO feature_pos = features%route_chunk_offsets(point_pos), &
            features%route_chunk_offsets(point_pos + 1) - 1
            irow = features%route_chunk_rows(feature_pos)
            CPASSERT(irow >= 1 .AND. irow <= features%chunk_feature_count)
            IF (collapse_spin_grads) THEN
               TARGET(base + 1) = TARGET(base + 1) + 0.5_dp*density_values(irow, 1)
               TARGET(base + 2) = TARGET(base + 2) + 0.5_dp*grad_values(irow, 1, 1)
               TARGET(base + 3) = TARGET(base + 3) + 0.5_dp*grad_values(irow, 2, 1)
               TARGET(base + 4) = TARGET(base + 4) + 0.5_dp*grad_values(irow, 3, 1)
               TARGET(base + 5) = TARGET(base + 5) + 0.5_dp*kin_values(irow, 1)
            ELSE
               TARGET(base + 1:base + 2) = TARGET(base + 1:base + 2) + density_values(irow, :)
               TARGET(base + 3) = TARGET(base + 3) + grad_values(irow, 1, 1)
               TARGET(base + 4) = TARGET(base + 4) + grad_values(irow, 2, 1)
               TARGET(base + 5) = TARGET(base + 5) + grad_values(irow, 3, 1)
               TARGET(base + 6) = TARGET(base + 6) + grad_values(irow, 1, 2)
               TARGET(base + 7) = TARGET(base + 7) + grad_values(irow, 2, 2)
               TARGET(base + 8) = TARGET(base + 8) + grad_values(irow, 3, 2)
               TARGET(base + 9:base + 10) = TARGET(base + 9:base + 10) + kin_values(irow, :)
            END IF
         END DO
      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE pack_atom_chunk_values

! **************************************************************************************************
!> \brief Scatter routed atom-chunk gradients into local grid-row order.
!> \param features ...
!> \param recv_grad_buffer ...
!> \param collapse_spin_grads ...
!> \param density_grad ...
!> \param grad_grad ...
!> \param kin_grad ...
! **************************************************************************************************
   SUBROUTINE scatter_routed_atom_chunk_grads(features, recv_grad_buffer, collapse_spin_grads, &
                                              density_grad, grad_grad, kin_grad)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: recv_grad_buffer
      LOGICAL, INTENT(IN)                                :: collapse_spin_grads
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: density_grad
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: grad_grad
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: kin_grad

      INTEGER                                            :: base, local_row, nflat_local, &
                                                            nroute_grad_per_point, nroute_points, &
                                                            point_pos, row_route_pos

      nflat_local = features%nflat_local
      nroute_points = SIZE(features%route_send_local_rows)
      nroute_grad_per_point = ngrad_per_point
      IF (collapse_spin_grads) nroute_grad_per_point = ncollapsed_grad_per_point
      CPASSERT(SIZE(recv_grad_buffer) >= nroute_grad_per_point*nroute_points)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(collapse_spin_grads, density_grad, features, grad_grad, kin_grad, nflat_local, &
!$OMP        nroute_grad_per_point, nroute_points, recv_grad_buffer) &
!$OMP PRIVATE(base, local_row, point_pos, row_route_pos)
      DO local_row = 1, nflat_local
         density_grad(local_row, :) = 0.0_dp
         grad_grad(local_row, :, :) = 0.0_dp
         kin_grad(local_row, :) = 0.0_dp
         DO row_route_pos = features%route_row_offsets(local_row), &
            features%route_row_offsets(local_row + 1) - 1
            point_pos = features%route_row_positions(row_route_pos)
            CPASSERT(point_pos >= 1 .AND. point_pos <= nroute_points)
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
      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE scatter_routed_atom_chunk_grads

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
      CALL torch_tensor_grad_batch3(features%density_t, features%grad_t, features%kin_t, &
                                    density_grad_t, grad_grad_t, kin_grad_t)
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
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: density_grad, kin_grad
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: grad_grad

      INTEGER                                            :: nroute_grad_per_point, &
                                                            nroute_recv_points, nroute_send_points, &
                                                            phase_handle
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: route_grad_return_recv_counts, &
                                                            route_grad_return_recv_displs, &
                                                            route_grad_return_send_counts, &
                                                            route_grad_return_send_displs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: recv_grad_buffer, send_grad_buffer

      CPASSERT(features%uses_atom_chunks)

      nroute_recv_points = SUM(features%route_point_recv_counts)
      nroute_send_points = SIZE(features%route_send_local_rows)
      CPASSERT(SUM(features%route_point_send_counts) == nroute_send_points)
      CPASSERT(SIZE(features%route_chunk_offsets) == nroute_recv_points + 1)

      nroute_grad_per_point = ngrad_per_point
      IF (features%uses_collapsed_rks_dynamic) THEN
         nroute_grad_per_point = ncollapsed_grad_per_point
      END IF
      CALL prepare_route_exchange( &
         features%route_point_recv_counts, features%route_point_recv_displs, &
         features%route_point_send_counts, features%route_point_send_displs, &
         nroute_grad_per_point, nroute_recv_points, nroute_send_points, &
         route_grad_return_send_counts, route_grad_return_send_displs, &
         route_grad_return_recv_counts, route_grad_return_recv_displs, &
         send_grad_buffer, recv_grad_buffer)

      IF (features%chunk_feature_count > 0) THEN
         CALL timeset("skala_gpw_grad_torch_pack", phase_handle)
         CALL pack_atom_chunk_grads(features, send_grad_buffer, &
                                    features%uses_collapsed_rks_dynamic)
         CALL timestop(phase_handle)
      END IF

      CALL timeset("skala_gpw_grad_route_comm", phase_handle)
      CALL group%alltoall(send_grad_buffer, route_grad_return_send_counts, &
                          route_grad_return_send_displs, recv_grad_buffer, &
                          route_grad_return_recv_counts, route_grad_return_recv_displs)
      CALL timestop(phase_handle)

      CALL timeset("skala_gpw_grad_route_scatter", phase_handle)
      CALL scatter_routed_atom_chunk_grads(features, recv_grad_buffer, &
                                           features%uses_collapsed_rks_dynamic, &
                                           density_grad, grad_grad, kin_grad)
      CALL timestop(phase_handle)

      DEALLOCATE (recv_grad_buffer, route_grad_return_recv_counts, &
                  route_grad_return_recv_displs, route_grad_return_send_counts, &
                  route_grad_return_send_displs, send_grad_buffer)

   END SUBROUTINE fetch_and_gather_atom_chunk_grads

! **************************************************************************************************
!> \brief Fetch and route atom-chunk feature gradients.
!>
!> In full mode it routes grid-coordinate, grid-weight, atomic-grid-weight and
!> atomic-coordinate adjoints. In weight-only mode it fetches and routes only
!> dE/d(grid_weights), which is all the accurate-XCINT weight derivative needs.
!> \param features ...
!> \param group ...
!> \param weight_only ...
!> \param local_grid_coord_grad ...
!> \param local_grid_weight_grad ...
!> \param local_atomic_grid_weight_grad ...
!> \param atom_coord_grad ...
! **************************************************************************************************
   SUBROUTINE fetch_and_route_atom_chunk_grads(features, group, weight_only, &
                                               local_grid_coord_grad, local_grid_weight_grad, &
                                               local_atomic_grid_weight_grad, atom_coord_grad)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      LOGICAL, INTENT(IN)                                :: weight_only
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT), OPTIONAL                           :: local_grid_coord_grad, atom_coord_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT), OPTIONAL                           :: local_grid_weight_grad, &
                                                            local_atomic_grid_weight_grad

      INTEGER                                            :: base, comps, feature_pos, irow, &
                                                            local_feature, natom, nfeature_local, &
                                                            phase_handle
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: recv_counts, recv_displs, send_counts, &
                                                            send_displs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: recv_buffer, send_buffer
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: chunk_atom_coord_grad, chunk_coord_grad
      REAL(KIND=dp), DIMENSION(:), POINTER               :: chunk_atomic_grid_weight_grad, &
                                                            chunk_grid_weight_grad
      TYPE(torch_tensor_type), DIMENSION(4)              :: grad_batch_t
      TYPE(torch_tensor_type)                            :: grid_weight_grad_t

      CPASSERT(features%uses_atom_chunks)
      CPASSERT(PRESENT(local_grid_weight_grad))

      comps = 1
      IF (.NOT. weight_only) comps = 5
      nfeature_local = SIZE(features%route_send_local_features)

      ALLOCATE (local_grid_weight_grad(nfeature_local))
      local_grid_weight_grad = 0.0_dp
      IF (.NOT. weight_only) THEN
         CPASSERT(PRESENT(local_grid_coord_grad))
         CPASSERT(PRESENT(local_atomic_grid_weight_grad))
         CPASSERT(PRESENT(atom_coord_grad))
         natom = SIZE(features%atomic_grid_sizes)
         ALLOCATE (local_grid_coord_grad(3, nfeature_local), &
                   local_atomic_grid_weight_grad(nfeature_local), atom_coord_grad(3, natom))
         local_grid_coord_grad = 0.0_dp
         local_atomic_grid_weight_grad = 0.0_dp
         atom_coord_grad = 0.0_dp
      END IF

      CALL prepare_route_exchange( &
         features%route_feature_recv_counts, features%route_feature_recv_displs, &
         features%route_feature_send_counts, features%route_feature_send_displs, &
         comps, features%chunk_feature_count, nfeature_local, &
         send_counts, send_displs, recv_counts, recv_displs, &
         send_buffer, recv_buffer, zero_buffers=.TRUE.)

      IF (features%chunk_feature_count > 0) THEN
         CALL timeset("skala_gpw_feature_grad_torch", phase_handle)
         IF (weight_only) THEN
            NULLIFY (chunk_grid_weight_grad)
            CALL torch_tensor_grad(features%grid_weights_t, grid_weight_grad_t)
            CALL torch_tensor_data_ptr(grid_weight_grad_t, chunk_grid_weight_grad)
            DO feature_pos = 1, features%chunk_feature_count
               irow = features%route_recv_chunk_rows(feature_pos)
               send_buffer(feature_pos) = chunk_grid_weight_grad(irow)
            END DO
            CALL torch_tensor_release(grid_weight_grad_t)
         ELSE
            NULLIFY (chunk_coord_grad, chunk_grid_weight_grad, chunk_atomic_grid_weight_grad, &
                     chunk_atom_coord_grad)
            CALL torch_tensor_grad_batch( &
               [features%grid_coords_t, features%grid_weights_t, &
                features%atomic_grid_weights_t, features%coarse_0_atomic_coords_t], &
               grad_batch_t)
            CALL torch_tensor_data_ptr(grad_batch_t(1), chunk_coord_grad)
            CALL torch_tensor_data_ptr(grad_batch_t(2), chunk_grid_weight_grad)
            CALL torch_tensor_data_ptr(grad_batch_t(3), chunk_atomic_grid_weight_grad)
            CALL torch_tensor_data_ptr(grad_batch_t(4), chunk_atom_coord_grad)
            DO feature_pos = 1, features%chunk_feature_count
               irow = features%route_recv_chunk_rows(feature_pos)
               base = comps*(feature_pos - 1)
               send_buffer(base + 1:base + 3) = chunk_coord_grad(:, irow)
               send_buffer(base + 4) = chunk_grid_weight_grad(irow)
               send_buffer(base + 5) = chunk_atomic_grid_weight_grad(irow)
            END DO
            DO irow = 1, features%chunk_natom
               atom_coord_grad(:, features%chunk_atom_begin + irow - 1) = &
                  chunk_atom_coord_grad(:, irow)
            END DO
            CALL torch_tensor_release(grad_batch_t(1))
            CALL torch_tensor_release(grad_batch_t(2))
            CALL torch_tensor_release(grad_batch_t(3))
            CALL torch_tensor_release(grad_batch_t(4))
         END IF
         CALL timestop(phase_handle)
      END IF

      CALL group%alltoall(send_buffer, send_counts, send_displs, recv_buffer, recv_counts, &
                          recv_displs)
      IF (.NOT. weight_only) CALL group%sum(atom_coord_grad)

      DO feature_pos = 1, nfeature_local
         local_feature = features%route_send_local_features(feature_pos)
         IF (weight_only) THEN
            local_grid_weight_grad(local_feature) = recv_buffer(feature_pos)
         ELSE
            base = comps*(feature_pos - 1)
            local_grid_coord_grad(:, local_feature) = recv_buffer(base + 1:base + 3)
            local_grid_weight_grad(local_feature) = recv_buffer(base + 4)
            local_atomic_grid_weight_grad(local_feature) = recv_buffer(base + 5)
         END IF
      END DO

      DEALLOCATE (send_buffer, recv_buffer, send_counts, send_displs, recv_counts, recv_displs)

   END SUBROUTINE fetch_and_route_atom_chunk_grads

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
                     DO jdir = 1, 3
                        tmp = -grad_i*drho(jdir)%array(i, j, k)
                        virial_xc(jdir, idir) = virial_xc(jdir, idir) + tmp
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
                     DO jdir = 1, 3
                        tmp = 0.0_dp
                        DO ispin = 1, 2
                           IF (ispin == 1) THEN
                              tmp = tmp - grad_grad(ipt, idir, ispin)*drhoa(jdir)%array(i, j, k)
                           ELSE
                              tmp = tmp - grad_grad(ipt, idir, ispin)*drhob(jdir)%array(i, j, k)
                           END IF
                        END DO
                        virial_xc(jdir, idir) = virial_xc(jdir, idir) + tmp
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
!> \param root_rank ...
!> \param print_components ...
!> \param atom_coord_grad_in ...
!> \param grid_coord_grad_in ...
! **************************************************************************************************
   SUBROUTINE build_static_coordinate_virial(virial_xc, features, root_rank, &
                                             print_components, atom_coord_grad_in, grid_coord_grad_in)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: virial_xc
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      LOGICAL, INTENT(IN)                                :: root_rank
      LOGICAL, INTENT(IN), OPTIONAL                      :: print_components
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), TARGET :: atom_coord_grad_in, grid_coord_grad_in

      INTEGER                                            :: feature_pos, iatom, idir, iw, jdir, &
                                                            local_row, row
      LOGICAL                                            :: my_print_components
      REAL(KIND=dp)                                      :: tmp
      REAL(KIND=dp), DIMENSION(3)                        :: atom_grad_sum, grid_grad_sum
      REAL(KIND=dp), DIMENSION(3, 3)                     :: atom_virial, grid_virial
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: atom_coord_grad, grid_coord_grad

      my_print_components = .FALSE.
      IF (PRESENT(print_components)) my_print_components = print_components

      atom_coord_grad => atom_coord_grad_in
      grid_coord_grad => grid_coord_grad_in

      grid_virial = 0.0_dp
      atom_virial = 0.0_dp
      grid_grad_sum = 0.0_dp
      atom_grad_sum = 0.0_dp
      DO local_row = 1, features%nflat_local
         DO feature_pos = features%local_feature_offsets(local_row), &
            features%local_feature_offsets(local_row + 1) - 1
            row = features%local_feature_rows(feature_pos)
            grid_grad_sum(:) = grid_grad_sum(:) + grid_coord_grad(:, row)
            DO idir = 1, 3
               DO jdir = 1, 3
                  tmp = grid_coord_grad(idir, row)*features%grid_coords(jdir, row)
                  grid_virial(idir, jdir) = grid_virial(idir, jdir) + tmp
                  virial_xc(idir, jdir) = virial_xc(idir, jdir) + tmp
               END DO
            END DO
         END DO
      END DO

      IF (root_rank) THEN
         DO iatom = 1, SIZE(features%coarse_0_atomic_coords, 2)
            atom_grad_sum(:) = atom_grad_sum(:) + atom_coord_grad(:, iatom)
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
            WRITE (iw, "(T2,A,1X,3ES20.10)") "SKALA_GPW| XC coordinate gradient grid sum", &
               grid_grad_sum
            WRITE (iw, "(T2,A,1X,3ES20.10)") "SKALA_GPW| XC coordinate gradient atom sum", &
               atom_grad_sum
            WRITE (iw, "(T2,A,1X,3ES20.10)") "SKALA_GPW| XC coordinate gradient total sum", &
               grid_grad_sum + atom_grad_sum
         END IF
      END IF

   END SUBROUTINE build_static_coordinate_virial

! **************************************************************************************************
!> \brief Add residual SKALA weight-feature contributions to the XC virial.
!> \param virial_xc ...
!> \param features ...
!> \param exc ...
!> \param root_rank ...
!> \param print_components ...
!> \param grid_weight_grad_in ...
!> \param atomic_grid_weight_grad_in ...
! **************************************************************************************************
   SUBROUTINE build_weight_virial(virial_xc, features, exc, root_rank, print_components, &
                                  grid_weight_grad_in, atomic_grid_weight_grad_in)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: virial_xc
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      REAL(KIND=dp), INTENT(IN)                          :: exc
      LOGICAL, INTENT(IN)                                :: root_rank
      LOGICAL, INTENT(IN), OPTIONAL                      :: print_components
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), TARGET    :: grid_weight_grad_in, &
                                                            atomic_grid_weight_grad_in

      INTEGER                                            :: feature_pos, idir, iw, local_row, row
      LOGICAL                                            :: my_print_components
      REAL(KIND=dp)                                      :: atomic_tmp, exc_tmp, grid_tmp, tmp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atomic_grid_weight_grad, grid_weight_grad

      my_print_components = .FALSE.
      IF (PRESENT(print_components)) my_print_components = print_components

      grid_weight_grad => grid_weight_grad_in
      atomic_grid_weight_grad => atomic_grid_weight_grad_in

      grid_tmp = 0.0_dp
      atomic_tmp = 0.0_dp
      DO local_row = 1, features%nflat_local
         DO feature_pos = features%local_feature_offsets(local_row), &
            features%local_feature_offsets(local_row + 1) - 1
            row = features%local_feature_rows(feature_pos)
            grid_tmp = grid_tmp + grid_weight_grad(row)*features%grid_weights(row)
            atomic_tmp = atomic_tmp + &
                         atomic_grid_weight_grad(row)*features%atomic_grid_weights(row)
         END DO
      END DO
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
      LOGICAL                                            :: valid_spin_shape
      REAL(KIND=dp)                                      :: dvol_inv
      TYPE(pw_c1d_gs_type)                               :: tmp_g, vxc_g
      TYPE(pw_r3d_rs_type), DIMENSION(3)                 :: grad_pw

      nspins = SIZE(rho_r)
      bo = rho_r(1)%pw_grid%bounds_local
      dvol_inv = 1.0_dp/rho_r(1)%pw_grid%dvol
      CPASSERT(SIZE(kin_grad, 2) == SIZE(density_grad, 2))
      CPASSERT(SIZE(grad_grad, 3) == SIZE(density_grad, 2))
      CPASSERT(SIZE(grad_grad, 2) == 3)
      valid_spin_shape = SIZE(density_grad, 2) == nspins .OR. &
                         (nspins == 1 .AND. SIZE(density_grad, 2) == 2)
      CPASSERT(valid_spin_shape)

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
                     IF (SIZE(density_grad, 2) == 1) THEN
                        vxc_rho(1)%array(i, j, k) = 0.5_dp*dvol_inv*density_grad(ipt, 1)
                        vxc_tau(1)%array(i, j, k) = 0.5_dp*dvol_inv*kin_grad(ipt, 1)
                        grad_pw(1)%array(i, j, k) = 0.5_dp*dvol_inv*grad_grad(ipt, 1, 1)
                        grad_pw(2)%array(i, j, k) = 0.5_dp*dvol_inv*grad_grad(ipt, 2, 1)
                        grad_pw(3)%array(i, j, k) = 0.5_dp*dvol_inv*grad_grad(ipt, 3, 1)
                     ELSE
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
                     END IF
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
!> \param group ...
! **************************************************************************************************
   SUBROUTINE print_native_grid_diagnostics(features, print_active, group)
      TYPE(skala_gpw_feature_type), INTENT(IN)           :: features
      LOGICAL, INTENT(IN)                                :: print_active

      CLASS(mp_comm_type), INTENT(IN)                    :: group

      INTEGER                                            :: atom_rows_max, atom_rows_min, &
                                                            chunk_rows_max, chunk_rows_min, iw
      REAL(KIND=dp)                                      :: chunk_imbalance

      chunk_rows_min = 0
      chunk_rows_max = 0
      IF (features%uses_atom_chunks) THEN
         chunk_rows_min = features%chunk_feature_count
         chunk_rows_max = features%chunk_feature_count
         CALL group%min(chunk_rows_min)
         CALL group%max(chunk_rows_max)
      END IF

      IF (.NOT. print_active) RETURN

      iw = cp_logger_get_default_io_unit()
      IF (iw <= 0) RETURN
      WRITE (UNIT=iw, FMT="(/,T2,A,1X,ES19.11)") &
         "SKALA_GPW| Native grid feature electrons", features%electron_count
      WRITE (UNIT=iw, FMT="(T2,A,1X,ES19.11)") &
         "SKALA_GPW| Native grid feature spin moment", features%spin_moment
      WRITE (UNIT=iw, FMT="(T2,A,1X,ES19.11)") &
         "SKALA_GPW| Native grid feature kinetic integral", features%kinetic_integral
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
         chunk_imbalance = REAL(chunk_rows_max, KIND=dp)/REAL(MAX(1, chunk_rows_min), KIND=dp)
         WRITE (UNIT=iw, FMT="(T2,A,1X,I0,1X,A,1X,I0,1X,A,1X,ES12.5)") &
            "SKALA_GPW| Native grid atom chunk row range", chunk_rows_min, &
            "to", chunk_rows_max, "imbalance", chunk_imbalance
      END IF

   END SUBROUTINE print_native_grid_diagnostics

! **************************************************************************************************
!> \brief Resolve the GPU/CUDA execution mode for the SKALA_NATIVE path.
!> \param xc_section ...
!> \return true when Torch should run on CUDA
! **************************************************************************************************
   FUNCTION skala_execution_use_cuda(xc_section) RESULT(use_cuda)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL                                            :: use_cuda

      INTEGER                                            :: gpu_mode
      LOGICAL                                            :: backend_has_gpu
      TYPE(section_vals_type), POINTER                   :: gauxc_section, native_section
      TYPE(skala_execution_target_type)                  :: target

      gpu_mode = do_gpu_auto
      native_section => get_skala_native_section(xc_section)
      IF (ASSOCIATED(native_section)) THEN
         CALL section_vals_val_get(native_section, "USE_GPU", i_val=gpu_mode)
      ELSE
         gauxc_section => get_gauxc_section(xc_section)
         IF (ASSOCIATED(gauxc_section)) THEN
            CALL section_vals_val_get(gauxc_section, "USE_GPU", i_val=gpu_mode)
         END IF
      END IF

      backend_has_gpu = .FALSE.
#if defined(__LIBTORCH_CUDA)
      backend_has_gpu = .TRUE.
#endif
      CALL ensure_skala_torch_threads()
      CALL resolve_skala_execution_target(gpu_mode, backend_has_gpu, "SKALA", TARGET)
      use_cuda = target%use_cuda

   END FUNCTION skala_execution_use_cuda

! **************************************************************************************************
!> \brief Configure CUDA device selection for the native SKALA GPW Torch path.
!> \param use_cuda ...
!> \param group ...
!> \return selected CUDA device, or -1 when CUDA is disabled
! **************************************************************************************************
   FUNCTION configure_native_grid_cuda(use_cuda, group) RESULT(selected_device)
      LOGICAL, INTENT(IN)                                :: use_cuda

      CLASS(mp_comm_type), INTENT(IN)                    :: group

      INTEGER                                            :: cuda_device_count, iw, selected_device
      TYPE(skala_execution_target_type)                  :: target

      selected_device = -1

      IF (.NOT. use_cuda) RETURN

      CALL select_cuda_target("SKALA", TARGET)
      selected_device = target%device
      cuda_device_count = torch_cuda_device_count()

      IF (group%mepos /= 0) RETURN
      IF (selected_device == skala_state%logged_cuda_device .AND. &
          cuda_device_count == skala_state%logged_cuda_device_count .AND. &
          group%num_pe == skala_state%logged_cuda_nproc) RETURN

      iw = cp_logger_get_default_io_unit()
      IF (iw <= 0) RETURN
      WRITE (UNIT=iw, FMT="(/,T2,A,1X,I0,1X,A,1X,I0)") &
         "SKALA_GPW| Native grid Torch CUDA device", selected_device, &
         "of", cuda_device_count

      skala_state%logged_cuda_device = selected_device
      skala_state%logged_cuda_device_count = cuda_device_count
      skala_state%logged_cuda_nproc = group%num_pe

   END FUNCTION configure_native_grid_cuda

! **************************************************************************************************
!> \brief Resolve the SKALA model path, execution target and device once per XC section.
!>        Repeated per-atom section scans, environment lookups and model file
!>        stamps are avoided; the context is invalidated by the cache release.
!> \param xc_section XC settings, including the native model precision.
!> \param group ...
!> \param use_cuda ...
!> \param cuda_device ...
!> \param higher_derivatives ...
! **************************************************************************************************
   SUBROUTINE ensure_model_context(xc_section, group, use_cuda, cuda_device, higher_derivatives)
      TYPE(section_vals_type), POINTER                   :: xc_section

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      LOGICAL, INTENT(OUT), OPTIONAL                     :: use_cuda
      INTEGER, INTENT(OUT)                               :: cuda_device
      LOGICAL, INTENT(IN), OPTIONAL                      :: higher_derivatives

      CHARACTER(len=default_path_length)                 :: model_path
      LOGICAL                                            :: prepare_derivatives

      prepare_derivatives = .FALSE.
      IF (PRESENT(higher_derivatives)) prepare_derivatives = higher_derivatives

      IF (skala_state%model_context%valid) THEN
         IF (ASSOCIATED(skala_state%model_context%xc_section, xc_section)) THEN
            IF (prepare_derivatives .AND. .NOT. skala_state%model_higher_derivatives) THEN
               CALL get_skala_model_path(xc_section, model_path)
               CALL ensure_model_loaded(model_path, skala_state%model_context%device, xc_section, &
                                        higher_derivatives=.TRUE.)
            END IF
            IF (PRESENT(use_cuda)) use_cuda = skala_state%model_context%use_cuda
            cuda_device = skala_state%model_context%device
            RETURN
         END IF
      END IF

      CALL get_skala_model_path(xc_section, model_path)
      skala_state%model_context%use_cuda = skala_execution_use_cuda(xc_section)
      skala_state%model_context%device = configure_native_grid_cuda(skala_state%model_context%use_cuda, group)
      CALL ensure_model_loaded(model_path, skala_state%model_context%device, xc_section, &
                               higher_derivatives=prepare_derivatives)
      skala_state%model_context%xc_section => xc_section
      skala_state%model_context%valid = .TRUE.
      IF (PRESENT(use_cuda)) use_cuda = skala_state%model_context%use_cuda
      cuda_device = skala_state%model_context%device

   END SUBROUTINE ensure_model_context

! **************************************************************************************************
!> \brief Load and cache the TorchScript SKALA model.
!> \param model_path ...
!> \param cuda_device ...
!> \param xc_section XC settings, including the native model precision.
!> \param higher_derivatives ...
! **************************************************************************************************
   SUBROUTINE ensure_model_loaded(model_path, cuda_device, xc_section, higher_derivatives)
      CHARACTER(len=*), INTENT(IN)                       :: model_path
      INTEGER, INTENT(IN)                                :: cuda_device
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN), OPTIONAL                      :: higher_derivatives

      CHARACTER(len=default_path_length)                 :: model_stamp
      INTEGER                                            :: iw
      LOGICAL                                            :: fp64, prepare_derivatives
      TYPE(section_vals_type), POINTER                   :: native_section

      fp64 = .TRUE.
      prepare_derivatives = .FALSE.
      IF (PRESENT(higher_derivatives)) prepare_derivatives = higher_derivatives
      native_section => get_skala_native_section(xc_section)
      IF (ASSOCIATED(native_section)) THEN
         CALL section_vals_val_get(native_section, "MODEL_FP64", l_val=fp64)
      END IF
      fp64 = fp64 .OR. prepare_derivatives
      CALL torch_file_stamp(TRIM(model_path), model_stamp)

      IF (skala_state%model_loaded) THEN
         IF (TRIM(skala_state%model_path) == TRIM(model_path) .AND. &
             skala_state%model_cuda_device == cuda_device .AND. &
             (skala_state%model_fp64 .EQV. fp64) .AND. &
             (.NOT. prepare_derivatives .OR. skala_state%model_higher_derivatives) .AND. &
             TRIM(skala_state%model_stamp) == TRIM(model_stamp)) RETURN
         CALL skala_torch_model_release(skala_state%model)
         skala_state%model_loaded = .FALSE.
      END IF

      CALL skala_torch_model_load(skala_state%model, TRIM(model_path), fp64=fp64, &
                                  higher_derivatives=prepare_derivatives)
      skala_state%model_path = model_path
      skala_state%model_stamp = model_stamp
      skala_state%model_cuda_device = cuda_device
      skala_state%model_fp64 = fp64
      skala_state%model_higher_derivatives = prepare_derivatives
      skala_state%model_loaded = .TRUE.

      IF (cuda_device >= 0) THEN
         iw = cp_logger_get_default_io_unit()
         IF (iw > 0) THEN
            WRITE (UNIT=iw, FMT="(T2,A,1X,I0)") &
               "SKALA_GPW| LibTorch CUDA model active on device", cuda_device
         END IF
      END IF

   END SUBROUTINE ensure_model_loaded

! **************************************************************************************************
!> \brief Resolve the SKALA TorchScript model path from the SKALA_NATIVE subsection.
!> \param xc_section ...
!> \param model_path ...
! **************************************************************************************************
   SUBROUTINE get_skala_model_path(xc_section, model_path)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      CHARACTER(len=default_path_length), INTENT(OUT)    :: model_path

      CHARACTER(len=default_path_length)                 :: model_name
      LOGICAL                                            :: use_cuda
      TYPE(section_vals_type), POINTER                   :: gauxc_section, native_section

      use_cuda = skala_execution_use_cuda(xc_section)
      native_section => get_skala_native_section(xc_section)
      IF (ASSOCIATED(native_section)) THEN
         CALL section_vals_val_get(native_section, "MODEL", c_val=model_name)
      ELSE
         gauxc_section => get_gauxc_section(xc_section)
         IF (.NOT. ASSOCIATED(gauxc_section)) THEN
            CPABORT("Native SKALA GPW requires an XC_FUNCTIONAL%SKALA_NATIVE section")
         END IF
         CALL section_vals_val_get(gauxc_section, "MODEL", c_val=model_name)
      END IF
      CALL resolve_skala_model_path(model_name, use_cuda, model_path)

   END SUBROUTINE get_skala_model_path

! **************************************************************************************************
!> \brief Release the reusable per-atom workspace (host arrays, tensors, dict).
! **************************************************************************************************
   SUBROUTINE release_atom_ws()
      IF (skala_state%atom_ws%tensors_active) THEN
         CALL torch_tensor_release(skala_state%atom_ws%grid_coords_t)
         CALL torch_tensor_release(skala_state%atom_ws%grid_weights_t)
         CALL torch_tensor_release(skala_state%atom_ws%atomic_grid_weights_t)
         CALL torch_tensor_release(skala_state%atom_ws%coarse_0_atomic_coords_t)
         CALL torch_tensor_release(skala_state%atom_ws%density_t)
         CALL torch_tensor_release(skala_state%atom_ws%grad_t)
         CALL torch_tensor_release(skala_state%atom_ws%kin_t)
         skala_state%atom_ws%tensors_active = .FALSE.
      END IF
      IF (skala_state%atom_ws%int_tensors_active) THEN
         CALL torch_tensor_release(skala_state%atom_ws%atomic_grid_sizes_t)
         CALL torch_tensor_release(skala_state%atom_ws%atomic_grid_size_bound_shape_t)
         skala_state%atom_ws%int_tensors_active = .FALSE.
      END IF
      IF (skala_state%atom_ws%dict_active) THEN
         CALL torch_dict_release(skala_state%atom_ws%inputs)
         skala_state%atom_ws%dict_active = .FALSE.
      END IF
      IF (ALLOCATED(skala_state%atom_ws%density)) THEN
         DEALLOCATE (skala_state%atom_ws%density, skala_state%atom_ws%grad, &
                     skala_state%atom_ws%kin, skala_state%atom_ws%grid_coords, &
                     skala_state%atom_ws%grid_weights, skala_state%atom_ws%atomic_grid_weights, &
                     skala_state%atom_ws%atomic_grid_sizes, &
                     skala_state%atom_ws%coarse_0_atomic_coords, &
                     skala_state%atom_ws%atomic_grid_size_bound_shape, &
                     skala_state%atom_ws%weight_b, skala_state%atom_ws%seg_offsets, &
                     skala_state%atom_ws%seg_lengths, skala_state%atom_ws%exc_segments)
      END IF
      skala_state%atom_ws%na = -1
      skala_state%atom_ws%nr = -1
      skala_state%atom_ws%nspins = -1
      skala_state%atom_ws%nblocks = -1
      skala_state%atom_ws%nflat_tot = -1
      skala_state%atom_ws%device = -2
   END SUBROUTINE release_atom_ws

! **************************************************************************************************
!> \brief Release all module-level SKALA/Torch caches.
! **************************************************************************************************
   SUBROUTINE skala_gpw_release_caches()
      IF (.NOT. ASSOCIATED(skala_state)) RETURN
      IF (skala_state%model_loaded) THEN
         CALL skala_torch_model_release(skala_state%model)
         skala_state%model_loaded = .FALSE.
      END IF
      CALL torch_grad_pool_clear()
      CALL torch_pinned_release()
      CALL torch_model_metadata_clear()
      CALL skala_gpw_feature_release(skala_state%features)
      IF (skala_state%atom_cuda_tensors_active) THEN
         CALL torch_tensor_release(skala_state%atom_cuda_tensors%grid_coords_t)
         CALL torch_tensor_release(skala_state%atom_cuda_tensors%grid_weights_t)
         CALL torch_tensor_release(skala_state%atom_cuda_tensors%atomic_grid_weights_t)
         CALL torch_tensor_release(skala_state%atom_cuda_tensors%atomic_coords_t)
         CALL torch_tensor_release(skala_state%atom_cuda_tensors%density_t)
         CALL torch_tensor_release(skala_state%atom_cuda_tensors%grad_t)
         CALL torch_tensor_release(skala_state%atom_cuda_tensors%kin_t)
         skala_state%atom_cuda_tensors_active = .FALSE.
      END IF
      skala_state%atom_static_valid = .FALSE.
      skala_state%atom_static_device = -1
      skala_state%atom_static_nrow = -1
      skala_state%atom_static_natom = -1
      CALL release_atom_ws()
      skala_state%model_context%valid = .FALSE.
      NULLIFY (skala_state%model_context%xc_section)
      skala_state%stress_mode_entered = .FALSE.
      IF (ASSOCIATED(skala_state%grad_scratch%density_grad)) THEN
         DEALLOCATE (skala_state%grad_scratch%density_grad, skala_state%grad_scratch%grad_grad, skala_state%grad_scratch%kin_grad)
      END IF
      skala_state%weight_grad%valid = .FALSE.
      IF (ASSOCIATED(skala_state%weight_grad%gradient)) DEALLOCATE (skala_state%weight_grad%gradient)
      IF (ALLOCATED(skala_state%weight_grad%positions)) DEALLOCATE (skala_state%weight_grad%positions)
      skala_state%weight_grad = skala_weight_grad_cache_type()
      skala_state%model_context = skala_model_context_type()
      IF (ALLOCATED(skala_state%chunk_layout%atom_begin)) THEN
         DEALLOCATE (skala_state%chunk_layout%atom_begin)
      END IF
      IF (ALLOCATED(skala_state%chunk_layout%atom_count)) THEN
         DEALLOCATE (skala_state%chunk_layout%atom_count)
      END IF
      IF (ALLOCATED(skala_state%chunk_layout%row_begin)) THEN
         DEALLOCATE (skala_state%chunk_layout%row_begin)
      END IF
      IF (ALLOCATED(skala_state%chunk_layout%row_counts)) THEN
         DEALLOCATE (skala_state%chunk_layout%row_counts)
      END IF
      IF (ALLOCATED(skala_state%chunk_layout%grid_sizes)) THEN
         DEALLOCATE (skala_state%chunk_layout%grid_sizes)
      END IF
      skala_state%chunk_layout%max_rows = -1
      skala_state%logged_torch_threads = .FALSE.
      skala_state%logged_cuda_device = -3
      skala_state%logged_cuda_device_count = -1
      skala_state%logged_cuda_nproc = -1
      skala_state%model_path = ""
      skala_state%model_stamp = ""
      skala_state%model_fp64 = .FALSE.
      skala_state%model_cuda_device = -3
      skala_state%execution_resolved = .FALSE.
      skala_state%execution_target = skala_execution_target_type()
      CALL skala_gpw_features_release_cache()
   END SUBROUTINE skala_gpw_release_caches

END MODULE skala_gpw_functional

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Consolidated SKALA runtime state (layout, features, scratch, model context) and the
!>        cached XC section pointers, with create/get/release helpers.
! **************************************************************************************************
MODULE skala_env_types
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              int_8
   USE qs_native_grid_cache,            ONLY: native_grid_cache_type,&
                                              release_native_grid_cache
   USE skala_torch_api,                 ONLY: skala_torch_model_type
   USE torch_api,                       ONLY: torch_dict_type,&
                                              torch_tensor_type
   USE xc_input_constants,              ONLY: skala_gpw_atom_partition_hard
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: skala_env_create, skala_env_finalize, skala_env_get, skala_env_release, skala_state
   PUBLIC :: atom_chunk_layout_cache_type, skala_atom_ws_type, skala_execution_target_type
   PUBLIC :: skala_gapw_atom_cuda_tensor_cache_type, skala_gpw_feature_type
   PUBLIC :: skala_gpw_layout_cache_type, skala_gpw_scratch_type, skala_gpw_state_type
   PUBLIC :: skala_grad_scratch_type, skala_model_context_type, skala_weight_grad_cache_type
   PUBLIC :: native_atom_cache_type, native_need_entry_type, native_geometry_key_type

   TYPE :: skala_execution_target_type
      LOGICAL                                            :: use_cuda = .FALSE.
      INTEGER                                            :: device = -1
   END TYPE skala_execution_target_type

   TYPE :: skala_gpw_layout_cache_type
      INTEGER                                            :: chunk_atom_begin = 1, chunk_atom_end = 0, &
                                                            chunk_feature_begin = 1, &
                                                            chunk_feature_count = 0, chunk_natom = 0, &
                                                            natom = 0, nflat = 0, nflat_local = 0, &
                                                            npoint = 0, nproc = 0, &
                                                            cuda_device = -3, &
                                                            atom_partition = skala_gpw_atom_partition_hard
      INTEGER, DIMENSION(2, 3)                           :: bo = 0, bounds = 0
      INTEGER, DIMENSION(3)                              :: cell_perd = 0, npts = 0
      INTEGER                                            :: comm_handle = -1
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dynamic_counts, dynamic_displs, &
                                                            feature_source_points, &
                                                            local_feature_offsets, &
                                                            local_feature_rows, &
                                                            chunk_return_positions, route_chunk_offsets, &
                                                            route_chunk_rows, route_row_offsets, &
                                                            route_row_positions, route_send_local_rows, &
                                                            route_point_recv_counts, &
                                                            route_point_recv_displs, &
                                                            route_point_send_counts, &
                                                            route_point_send_displs, &
                                                            route_feature_send_counts, &
                                                            route_feature_send_displs, &
                                                            route_feature_recv_counts, &
                                                            route_feature_recv_displs, &
                                                            route_send_local_features, &
                                                            route_recv_chunk_rows
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: atomic_grid_sizes, chunk_atomic_grid_sizes
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :)  :: atomic_grid_size_bound_shape, &
                                                            chunk_atomic_grid_size_bound_shape
      TYPE(torch_dict_type)                              :: chunk_inputs
      TYPE(torch_dict_type)                              :: chunk_static_inputs
      TYPE(torch_dict_type)                              :: grad_inputs
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
      REAL(KIND=dp)                                      :: weight_sum = 0.0_dp, &
                                                            weight_sumsq = 0.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)                     :: cell_hmat = 0.0_dp, dh = 0.0_dp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: atomic_grid_weights, chunk_atomic_grid_weights, &
                                                            chunk_grid_weights, grid_weights, &
                                                            local_feature_atomic_grid_weights, &
                                                            local_feature_grid_weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: atom_coords, chunk_coarse_0_atomic_coords, &
                                                            chunk_grid_coords, coarse_0_atomic_coords, &
                                                            grid_coords, local_feature_grid_coords
      LOGICAL                                            :: active = .FALSE., has_weights = .FALSE., &
                                                            chunk_dynamic_input_views_active = .FALSE., &
                                                            chunk_dynamic_tensors_active = .FALSE., &
                                                            chunk_inputs_active = .FALSE., &
                                                            chunk_inputs_use_collapsed_rks = .FALSE., &
                                                            chunk_static_tensors_active = .FALSE., &
                                                            dynamic_tensors_active = .FALSE., &
                                                            grad_inputs_active = .FALSE., &
                                                            inputs_active = .FALSE., &
                                                            static_tensors_active = .FALSE.
   END TYPE skala_gpw_layout_cache_type

   TYPE :: skala_gpw_feature_type
      INTEGER                                            :: chunk_atom_begin = 1, &
                                                            chunk_feature_count = 0, chunk_natom = 0, &
                                                            nflat = 0, nflat_local = 0, &
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
      INTEGER, POINTER, DIMENSION(:)                    :: local_feature_offsets => NULL(), &
                                                           local_feature_rows => NULL()
      INTEGER, POINTER, DIMENSION(:)                    :: route_point_recv_counts => NULL(), &
                                                           route_point_recv_displs => NULL(), &
                                                           route_point_send_counts => NULL(), &
                                                           route_point_send_displs => NULL(), &
                                                           route_chunk_offsets => NULL(), &
                                                           route_chunk_rows => NULL(), &
                                                           route_row_offsets => NULL(), &
                                                           route_row_positions => NULL(), &
                                                           route_send_local_rows => NULL(), &
                                                           route_feature_send_counts => NULL(), &
                                                           route_feature_send_displs => NULL(), &
                                                           route_feature_recv_counts => NULL(), &
                                                           route_feature_recv_displs => NULL(), &
                                                           route_send_local_features => NULL(), &
                                                           route_recv_chunk_rows => NULL()
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)    :: atomic_grid_sizes
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :) :: atomic_grid_size_bound_shape
      REAL(KIND=dp), POINTER, DIMENSION(:)              :: atomic_grid_weights => NULL(), &
                                                           grid_weights => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:, :)           :: grid_coords => NULL()
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)       :: chunk_density, chunk_kin, &
                                                           coarse_0_atomic_coords, density, kin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)    :: chunk_grad, grad
      REAL(KIND=dp)                                      :: electron_count = 0.0_dp, &
                                                            grid_weight_sum = 0.0_dp, &
                                                            kinetic_integral = 0.0_dp, &
                                                            spin_moment = 0.0_dp
      LOGICAL                                            :: active = .FALSE., owns_coordinate_tensor = .FALSE., &
                                                            owns_grid_coordinate_tensor = .FALSE., &
                                                            owns_layout_arrays = .FALSE., &
                                                            owns_weight_tensors = .FALSE., &
                                                            owns_dynamic_tensors = .TRUE., &
                                                            owns_inputs = .TRUE., &
                                                            inputs_active = .FALSE., &
                                                            owns_static_tensors = .TRUE., &
                                                            uses_atom_chunks = .FALSE., &
                                                            uses_collapsed_rks_dynamic = .FALSE.
   END TYPE skala_gpw_feature_type

   TYPE :: skala_gpw_scratch_type
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: local_dynamic, global_dynamic, &
                                                         route_send_dynamic, route_recv_dynamic
      INTEGER, ALLOCATABLE, DIMENSION(:)              :: route_send_counts, route_send_displs, &
                                                         route_recv_counts, route_recv_displs
   END TYPE skala_gpw_scratch_type

   TYPE :: skala_gapw_atom_cuda_tensor_cache_type
      TYPE(torch_tensor_type)                            :: atomic_coords_t, atomic_grid_weights_t, &
                                                            density_t, grad_t, grid_coords_t, &
                                                            grid_weights_t, kin_t
   END TYPE skala_gapw_atom_cuda_tensor_cache_type

   TYPE :: skala_atom_ws_type
      INTEGER                                            :: na = -1
      INTEGER                                            :: nr = -1
      INTEGER                                            :: nspins = -1
      INTEGER                                            :: nblocks = -1
      INTEGER                                            :: nflat_tot = -1
      INTEGER                                            :: device = -2
      LOGICAL                                            :: tensors_active = .FALSE.
      LOGICAL                                            :: dict_active = .FALSE.
      LOGICAL                                            :: int_tensors_active = .FALSE.
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: atomic_grid_sizes, seg_lengths, &
                                                            seg_offsets
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :)  :: atomic_grid_size_bound_shape
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: atomic_grid_weights, exc_segments, &
                                                            grid_weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: coarse_0_atomic_coords, density, &
                                                            grid_coords, kin, weight_b
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: grad
      TYPE(torch_tensor_type)                            :: atomic_grid_size_bound_shape_t, &
                                                            atomic_grid_sizes_t, &
                                                            atomic_grid_weights_t, &
                                                            coarse_0_atomic_coords_t, density_t, &
                                                            grad_t, grid_coords_t, grid_weights_t, &
                                                            kin_t
      TYPE(torch_dict_type)                              :: inputs
   END TYPE skala_atom_ws_type

   TYPE :: skala_model_context_type
      LOGICAL                                            :: valid = .FALSE.
      TYPE(section_vals_type), POINTER                   :: xc_section => NULL()
      LOGICAL                                            :: use_cuda = .FALSE.
      INTEGER                                            :: device = -1
   END TYPE skala_model_context_type

   TYPE :: atom_chunk_layout_cache_type
      INTEGER                                            :: max_rows = -1
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_begin, atom_count, &
                                                            row_begin, row_counts
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: grid_sizes
   END TYPE atom_chunk_layout_cache_type

   TYPE :: skala_grad_scratch_type
      REAL(KIND=dp), POINTER, DIMENSION(:, :)    :: density_grad => NULL(), kin_grad => NULL()
      REAL(KIND=dp), POINTER, DIMENSION(:, :, :) :: grad_grad => NULL()
   END TYPE skala_grad_scratch_type

   TYPE :: skala_weight_grad_cache_type
      REAL(KIND=dp), POINTER, DIMENSION(:)       :: gradient => NULL()
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: positions
      REAL(KIND=dp)                              :: rho_sum = 0.0_dp, rho_sumsq = 0.0_dp, &
                                                    tau_sum = 0.0_dp, tau_sumsq = 0.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)             :: cell_hmat = 0.0_dp
      INTEGER, DIMENSION(3)                      :: cell_perd = 0
      INTEGER                                    :: nvals = -1, nspins = -1, atom_partition = -1
      LOGICAL                                    :: atom_chunks = .FALSE., valid = .FALSE.
   END TYPE skala_weight_grad_cache_type

   TYPE :: native_need_entry_type
      INTEGER(int_8) :: lin = 0_int_8
      INTEGER :: pos = 0
   END TYPE native_need_entry_type

   TYPE :: native_geometry_key_type
      REAL(dp) :: hmat(3, 3) = 0.0_dp
      INTEGER :: perd(3) = 0
      REAL(dp), ALLOCATABLE :: atoms(:, :)
      LOGICAL :: ready = .FALSE.
   END TYPE native_geometry_key_type

   TYPE :: native_atom_cache_type
      INTEGER :: native_need_nrows = -1, native_need_npts(3) = 0, native_need_nfg = 0, &
                 native_need_nneed = 0, native_need_nbucket = 0, native_need_nrecv = 0
      INTEGER :: native_need_group_handle = -1, native_need_num_pe = -1, &
                 native_need_mepos = -1, native_need_mode = -1
      INTEGER :: native_need_bounds_local(2, 3) = 0
      REAL(dp) :: native_need_dh_inv(3, 3) = 0.0_dp
      LOGICAL :: native_need_wrap = .FALSE., native_need_ready = .FALSE., &
                 native_need_replicated = .FALSE.
      INTEGER(int_8), ALLOCATABLE :: native_need_base(:, :), native_need_lin(:)
      INTEGER, ALLOCATABLE :: native_need_bucket_start(:), native_need_local(:), &
                              native_need_scount(:), native_need_sdispl(:), &
                              native_need_rcount(:), native_need_rdispl(:), &
                              native_need_sidx(:), native_need_ridx(:)
      REAL(dp), ALLOCATABLE :: native_need_values(:, :), native_need_reply(:)
      TYPE(native_need_entry_type), ALLOCATABLE :: native_need_bucket_entry(:)
      TYPE(native_geometry_key_type) :: native_need_geometry, native_comp_geometry
      INTEGER :: native_comp_nrows = -1
      LOGICAL :: native_comp_ready = .FALSE.
      INTEGER, ALLOCATABLE :: native_comp_local_atoms(:)
      REAL(dp), ALLOCATABLE :: native_comp_coords(:, :), native_comp_base_w(:), &
                               native_comp_grid_w(:)
      INTEGER :: adjoint_scratch_nthread = 0, adjoint_scratch_nfg = 0, adjoint_scratch_nneed = 0
      REAL(dp), ALLOCATABLE :: adjoint_thread_scratch(:, :, :)
   END TYPE native_atom_cache_type

   TYPE :: skala_gpw_state_type
      TYPE(skala_gpw_layout_cache_type)          :: layout
      LOGICAL                                    :: layout_mode_known = .FALSE.
      LOGICAL                                    :: layout_mode_full = .FALSE.
      TYPE(skala_gpw_scratch_type)               :: scratch
      TYPE(skala_gpw_feature_type)               :: features
      TYPE(skala_torch_model_type)               :: model
      CHARACTER(len=default_path_length)         :: model_path = ""
      CHARACTER(len=default_path_length)         :: model_stamp = ""
      LOGICAL                                    :: model_loaded = .FALSE.
      LOGICAL                                    :: model_fp64 = .FALSE.
      LOGICAL                                    :: model_higher_derivatives = .FALSE.
      INTEGER                                    :: model_cuda_device = -3
      LOGICAL                                    :: logged_torch_threads = .FALSE.
      INTEGER                                    :: logged_cuda_device = -3
      INTEGER                                    :: logged_cuda_device_count = -1
      INTEGER                                    :: logged_cuda_nproc = -1
      TYPE(skala_gapw_atom_cuda_tensor_cache_type) :: atom_cuda_tensors
      LOGICAL                                    :: atom_cuda_tensors_active = .FALSE.
      LOGICAL                                    :: atom_static_valid = .FALSE.
      INTEGER                                    :: atom_static_device = -1
      INTEGER                                    :: atom_static_nrow = -1
      INTEGER                                    :: atom_static_natom = -1
      TYPE(skala_atom_ws_type)                   :: atom_ws
      TYPE(skala_model_context_type)             :: model_context
      LOGICAL                                    :: stress_mode_entered = .FALSE.
      TYPE(atom_chunk_layout_cache_type)         :: chunk_layout
      TYPE(skala_grad_scratch_type)              :: grad_scratch
      TYPE(skala_weight_grad_cache_type)         :: weight_grad
      LOGICAL                                    :: execution_resolved = .FALSE.
      TYPE(skala_execution_target_type)          :: execution_target
      TYPE(native_grid_cache_type)                       :: native_grid_cache
      TYPE(section_vals_type), POINTER                   :: native_section => NULL()
      TYPE(section_vals_type), POINTER                   :: gauxc_section => NULL()
      TYPE(section_vals_type), POINTER                   :: resolved_xc_section => NULL()
      LOGICAL                                            :: sections_resolved = .FALSE.
      TYPE(native_atom_cache_type)                       :: native_atom
   END TYPE skala_gpw_state_type
   TYPE(skala_gpw_state_type), POINTER :: skala_state => NULL()

CONTAINS

! **************************************************************************************************
!> \brief Resolve and cache the XC_FUNCTIONAL subsections used by SKALA.
!> \param env ...
!> \param xc_section the XC input section
! **************************************************************************************************
   SUBROUTINE skala_env_create(env, xc_section)
      TYPE(skala_gpw_state_type), INTENT(INOUT), POINTER :: env
      TYPE(section_vals_type), INTENT(IN), OPTIONAL, &
         POINTER                                         :: xc_section

      IF (.NOT. ASSOCIATED(env)) ALLOCATE (env)
      skala_state => env
      NULLIFY (env%native_section, env%gauxc_section, env%resolved_xc_section)
      env%sections_resolved = .FALSE.
      IF (PRESENT(xc_section)) THEN
         IF (.NOT. ASSOCIATED(xc_section)) RETURN
         env%resolved_xc_section => xc_section
         env%native_section => find_named_section(xc_section, "SKALA_NATIVE")
         env%gauxc_section => find_named_section(xc_section, "GAUXC")
         env%sections_resolved = .TRUE.
      END IF
   END SUBROUTINE skala_env_create

! **************************************************************************************************
!> \brief Return the consolidated SKALA state instance.
!> \return pointer to the SKALA state
! **************************************************************************************************
   FUNCTION skala_env_get() RESULT(env)
      TYPE(skala_gpw_state_type), POINTER                :: env

      env => skala_state
   END FUNCTION skala_env_get

! **************************************************************************************************
!> \brief Release the SKALA section cache. Heavy Torch/layout caches are freed by
!>        skala_gpw_release_caches before this is called.
! **************************************************************************************************
   SUBROUTINE skala_env_release()
      IF (.NOT. ASSOCIATED(skala_state)) RETURN
      CALL release_native_grid_cache(skala_state%native_grid_cache)
      NULLIFY (skala_state%native_section, skala_state%gauxc_section, &
               skala_state%resolved_xc_section)
      skala_state%sections_resolved = .FALSE.
   END SUBROUTINE skala_env_release

! **************************************************************************************************
!> \brief Release the SKALA caches and deallocate the qs_env-owned state instance.
!> \param env the SKALA state owned by the qs_env
! **************************************************************************************************
   SUBROUTINE skala_env_finalize(env)
      TYPE(skala_gpw_state_type), INTENT(INOUT), POINTER :: env

      IF (.NOT. ASSOCIATED(env)) RETURN
      IF (ASSOCIATED(skala_state, env)) THEN
         CALL skala_env_release()
         NULLIFY (skala_state)
      END IF
      DEALLOCATE (env)
   END SUBROUTINE skala_env_finalize

! **************************************************************************************************
!> \brief Return the named XC_FUNCTIONAL subsection present in the input.
!> \param xc_section the XC input section
!> \param name subsection name
!> \return the subsection pointer or NULL
! **************************************************************************************************
   FUNCTION find_named_section(xc_section, name) RESULT(res)
      TYPE(section_vals_type), POINTER                   :: xc_section
      CHARACTER(len=*), INTENT(IN)                       :: name
      TYPE(section_vals_type), POINTER                   :: res

      INTEGER                                            :: i
      TYPE(section_vals_type), POINTER                   :: functionals, sub

      NULLIFY (res)
      IF (.NOT. ASSOCIATED(xc_section)) RETURN
      functionals => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      IF (.NOT. ASSOCIATED(functionals)) RETURN
      DO i = 1, functionals%section%n_subsections
         sub => functionals%subs_vals(i, 1)%section_vals
         IF (SIZE(sub%values, 2) == 0) CYCLE
         IF (TRIM(sub%section%name) == TRIM(name)) res => sub
      END DO
   END FUNCTION find_named_section

END MODULE skala_env_types

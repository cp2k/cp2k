!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines that build the integrals of the Vxc potential calculated
!>      for the atomic density in the basis set of spherical primitives
! **************************************************************************************************
MODULE qs_vxc_atom
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_array_utils,                  ONLY: cp_3d_r_cp_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE external_potential_types,        ONLY: gth_potential_type,&
                                              sgp_potential_type
   USE input_constants,                 ONLY: cdft_alpha_constraint,&
                                              cdft_beta_constraint,&
                                              cdft_charge_constraint,&
                                              cdft_magnetization_constraint,&
                                              xc_none
   USE input_section_types,             ONLY: section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                              get_number_of_lebedev_grid,&
                                              init_lebedev_grids,&
                                              lebedev_grid
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE orbital_pointers,                ONLY: indco,&
                                              indso,&
                                              nco,&
                                              ncoset,&
                                              nsoset
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL,&
                                              pw_grid_type
   USE pw_grids,                        ONLY: pw_grid_map_global_point
   USE pw_methods,                      ONLY: pw_axpy
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_cdft_grid,                    ONLY: cdft_point_context_create,&
                                              cdft_point_context_release,&
                                              cdft_point_context_type,&
                                              cdft_point_weights
   USE qs_cdft_types,                   ONLY: cdft_control_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                              create_grid_atom,&
                                              grid_atom_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              has_nlcc,&
                                              qs_kind_type
   USE qs_linres_types,                 ONLY: nablavks_atom_type
   USE qs_native_grid_cache,            ONLY: fetch_native_grid_stencil,&
                                              native_grid_interp_npts,&
                                              native_grid_interp_offset_max,&
                                              native_grid_interp_offset_min,&
                                              native_grid_interpolation_stencil_type,&
                                              prepare_native_grid_cache,&
                                              store_native_grid_stencil
   USE qs_rho_atom_methods,             ONLY: replicate_rho_atom_radial
   USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                              rho_atom_coeff,&
                                              rho_atom_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc_atom_utils,               ONLY: &
        calc_rho_angular, calc_rho_nlcc, calc_tau_atom, calc_weight_function, &
        create_tau_basis_cache, dgaVtaudgb, evaluate_nlcc_primitive_fields, gaVxcgb_GC, &
        gaVxcgb_noGC, gapw_atom_grid_support_radius, release_tau_basis_cache, tau_basis_cache_type
   USE skala_gpw_features,              ONLY: build_periodic_atom_image_layout,&
                                              periodic_atom_image_partition_from_layout,&
                                              skala_periodic_partition_type,&
                                              skala_state
   USE skala_gpw_functional,            ONLY: build_vxc_from_feature_grads,&
                                              native_skala_gapw_density_partition,&
                                              skala_gapw_atom_composite_energy,&
                                              skala_gapw_atom_vxc_of_r_pair,&
                                              skala_gapw_density_partition_hard_minus_soft,&
                                              skala_gapw_density_partition_hard_only,&
                                              skala_gapw_density_partition_none,&
                                              skala_gapw_density_partition_soft_only
   USE skala_model_selection,           ONLY: get_gauxc_section,&
                                              get_skala_native_section,&
                                              xc_section_uses_gauxc_model,&
                                              xc_section_uses_native_skala_grid
   USE spherical_harmonics,             ONLY: y_lm
   USE util,                            ONLY: get_limit
   USE virial_types,                    ONLY: virial_type
   USE xc_atom,                         ONLY: fill_rho_set,&
                                              vxc_of_r_epr,&
                                              vxc_of_r_new,&
                                              xc_rho_set_atom_update
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_create,&
                                              xc_dset_release,&
                                              xc_dset_zero_all
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_input_constants,              ONLY: skala_gapw_cp2k_default,&
                                              skala_gapw_direct_valence,&
                                              skala_gapw_paw_one_center,&
                                              skala_gapw_paw_one_center_split
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_get,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
   USE skala_env_types, ONLY: native_geometry_key_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_vxc_atom'

   PUBLIC :: calculate_vxc_atom, &
             calculate_vxc_atom_epr, &
             gapw_cdft_one_center, &
             add_gapw_atom_grid_interpolation_adjoint, &
             atom_grid_image_bounds, &
             interpolate_gapw_atom_grid_fields, &
             release_native_skala_caches

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Return true if a cached geometry snapshot matches the current cell and atoms.
!> \param key cached geometry snapshot
!> \param cell current cell
!> \param atom_coords all atom coordinates
!> \return ...
! **************************************************************************************************
   FUNCTION native_geometry_matches(key, cell, atom_coords) RESULT(matches)
      TYPE(native_geometry_key_type), INTENT(IN)         :: key
      TYPE(cell_type), POINTER                           :: cell
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: atom_coords
      LOGICAL                                            :: matches

      matches = key%ready .AND. ALLOCATED(key%atoms)
      IF (matches) matches = ALL(SHAPE(key%atoms) == SHAPE(atom_coords))
      IF (matches) matches = ALL(key%atoms == atom_coords)
      matches = matches .AND. ALL(key%hmat == cell%hmat) .AND. ALL(key%perd == cell%perd)
   END FUNCTION native_geometry_matches

! **************************************************************************************************
!> \brief Store the current cell and atom geometry into a cache snapshot.
!> \param key geometry snapshot to update
!> \param cell current cell
!> \param atom_coords all atom coordinates
! **************************************************************************************************
   SUBROUTINE store_native_geometry_key(key, cell, atom_coords)
      TYPE(native_geometry_key_type), INTENT(INOUT)      :: key
      TYPE(cell_type), POINTER                           :: cell
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: atom_coords

      IF (ALLOCATED(key%atoms)) THEN
         IF (.NOT. ALL(SHAPE(key%atoms) == SHAPE(atom_coords))) DEALLOCATE (key%atoms)
      END IF
      IF (.NOT. ALLOCATED(key%atoms)) THEN
         ALLOCATE (key%atoms(SIZE(atom_coords, 1), SIZE(atom_coords, 2)))
      END IF
      key%atoms(:, :) = atom_coords
      key%hmat = cell%hmat
      key%perd = cell%perd
      key%ready = .TRUE.
   END SUBROUTINE store_native_geometry_key

! **************************************************************************************************
!> \brief Ensure a grow-only per-thread adjoint scratch buffer of the given size.
!> \param nthread number of OpenMP threads
!> \param nfg number of adjoint channels
!> \param nneed number of local needed points
! **************************************************************************************************
   SUBROUTINE ensure_adjoint_thread_scratch(nthread, nfg, nneed)
      INTEGER, INTENT(IN)                                :: nthread, nfg, nneed

     IF (skala_state%native_atom%adjoint_scratch_nthread >= nthread .AND. skala_state%native_atom%adjoint_scratch_nfg >= nfg .AND. &
          skala_state%native_atom%adjoint_scratch_nneed >= nneed) RETURN
      IF (ALLOCATED(skala_state%native_atom%adjoint_thread_scratch)) DEALLOCATE (skala_state%native_atom%adjoint_thread_scratch)
      skala_state%native_atom%adjoint_scratch_nthread = MAX(nthread, skala_state%native_atom%adjoint_scratch_nthread)
      skala_state%native_atom%adjoint_scratch_nfg = MAX(nfg, skala_state%native_atom%adjoint_scratch_nfg)
      skala_state%native_atom%adjoint_scratch_nneed = MAX(nneed, skala_state%native_atom%adjoint_scratch_nneed)
      ALLOCATE (skala_state%native_atom%adjoint_thread_scratch( &
                skala_state%native_atom%adjoint_scratch_nfg, &
                skala_state%native_atom%adjoint_scratch_nneed, &
                skala_state%native_atom%adjoint_scratch_nthread))
   END SUBROUTINE ensure_adjoint_thread_scratch

! **************************************************************************************************
!> \brief Release the process-lifetime native SKALA caches (needed-point schedule,
!>        atom-composite coordinates and adjoint scratch).
! **************************************************************************************************
   SUBROUTINE release_native_skala_caches()
      IF (.NOT. ASSOCIATED(skala_state)) RETURN
      IF (ALLOCATED(skala_state%native_atom%native_need_base)) DEALLOCATE (skala_state%native_atom%native_need_base)
      CALL native_need_index_release()
      IF (ALLOCATED(skala_state%native_atom%native_need_local)) DEALLOCATE (skala_state%native_atom%native_need_local)
      IF (ALLOCATED(skala_state%native_atom%native_need_scount)) DEALLOCATE (skala_state%native_atom%native_need_scount)
      IF (ALLOCATED(skala_state%native_atom%native_need_sdispl)) DEALLOCATE (skala_state%native_atom%native_need_sdispl)
      IF (ALLOCATED(skala_state%native_atom%native_need_rcount)) DEALLOCATE (skala_state%native_atom%native_need_rcount)
      IF (ALLOCATED(skala_state%native_atom%native_need_rdispl)) DEALLOCATE (skala_state%native_atom%native_need_rdispl)
      IF (ALLOCATED(skala_state%native_atom%native_need_sidx)) DEALLOCATE (skala_state%native_atom%native_need_sidx)
      IF (ALLOCATED(skala_state%native_atom%native_need_ridx)) DEALLOCATE (skala_state%native_atom%native_need_ridx)
      IF (ALLOCATED(skala_state%native_atom%native_need_values)) DEALLOCATE (skala_state%native_atom%native_need_values)
      IF (ALLOCATED(skala_state%native_atom%native_need_reply)) DEALLOCATE (skala_state%native_atom%native_need_reply)
      IF (ALLOCATED(skala_state%native_atom%native_need_geometry%atoms)) THEN
         DEALLOCATE (skala_state%native_atom%native_need_geometry%atoms)
      END IF
      skala_state%native_atom%native_need_geometry%ready = .FALSE.
      skala_state%native_atom%native_need_ready = .FALSE.
      IF (ALLOCATED(skala_state%native_atom%native_comp_coords)) DEALLOCATE (skala_state%native_atom%native_comp_coords)
      IF (ALLOCATED(skala_state%native_atom%native_comp_base_w)) DEALLOCATE (skala_state%native_atom%native_comp_base_w)
      IF (ALLOCATED(skala_state%native_atom%native_comp_grid_w)) DEALLOCATE (skala_state%native_atom%native_comp_grid_w)
      IF (ALLOCATED(skala_state%native_atom%native_comp_local_atoms)) DEALLOCATE (skala_state%native_atom%native_comp_local_atoms)
      IF (ALLOCATED(skala_state%native_atom%native_comp_geometry%atoms)) THEN
         DEALLOCATE (skala_state%native_atom%native_comp_geometry%atoms)
      END IF
      skala_state%native_atom%native_comp_geometry%ready = .FALSE.
      skala_state%native_atom%native_comp_ready = .FALSE.
      skala_state%native_atom%native_comp_nrows = -1
      IF (ALLOCATED(skala_state%native_atom%adjoint_thread_scratch)) DEALLOCATE (skala_state%native_atom%adjoint_thread_scratch)
      skala_state%native_atom%adjoint_scratch_nthread = 0
      skala_state%native_atom%adjoint_scratch_nfg = 0
      skala_state%native_atom%adjoint_scratch_nneed = 0
   END SUBROUTINE release_native_skala_caches

! **************************************************************************************************
!> \brief Initialize atom-centered quadrature for native Skala layouts.
!> \param kind_set quantum kinds
!> \param dft_control DFT controls supplying the radial quadrature
! **************************************************************************************************
   SUBROUTINE ensure_native_skala_atom_grids(kind_set, dft_control)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control

      INTEGER                                            :: ikind, ll, na, nr, quadrature
      TYPE(grid_atom_type), POINTER                      :: grid_atom

      quadrature = dft_control%qs_control%gapw_control%quadrature
      CALL init_lebedev_grids()
      DO ikind = 1, SIZE(kind_set)
         NULLIFY (grid_atom)
         CALL get_qs_kind(kind_set(ikind), grid_atom=grid_atom, ngrid_ang=na, ngrid_rad=nr)
         IF (ASSOCIATED(grid_atom)) THEN
            IF (ASSOCIATED(grid_atom%weight) .AND. grid_atom%nr == nr) CYCLE
         ELSE
            CALL allocate_grid_atom(kind_set(ikind)%grid_atom)
            grid_atom => kind_set(ikind)%grid_atom
         END IF
         ll = get_number_of_lebedev_grid(n=na)
         na = lebedev_grid(ll)%n
         grid_atom%ng_sphere = na
         grid_atom%nr = nr
         CALL create_grid_atom(grid_atom, nr, na, 0, ll, quadrature)
      END DO
      CALL deallocate_lebedev_grids()

   END SUBROUTINE ensure_native_skala_atom_grids

! **************************************************************************************************
!> \brief Decide whether a kind contributes hard-minus-soft primitive fields.
!> \param paw_atom whether CP2K constructed a one-center representation for the kind
!> \param gapw_representation requested Skala pseudopotential GAPW representation
!> \param has_pseudopotential whether the kind uses a GTH or semi-global pseudopotential
!> \param zeff valence charge of the potential
!> \param zatom atomic number
!> \return true when hard-minus-soft fields contribute for this kind
! **************************************************************************************************
   PURE FUNCTION native_skala_uses_one_center_kind( &
      paw_atom, gapw_representation, has_pseudopotential, zeff, zatom) RESULT(use_one_center)
      LOGICAL, INTENT(IN)                                :: paw_atom
      INTEGER, INTENT(IN)                                :: gapw_representation
      LOGICAL, INTENT(IN)                                :: has_pseudopotential
      REAL(dp), INTENT(IN)                               :: zeff
      INTEGER, INTENT(IN)                                :: zatom
      LOGICAL                                            :: use_one_center

      use_one_center = paw_atom
      IF (.NOT. use_one_center) RETURN

      SELECT CASE (gapw_representation)
      CASE (skala_gapw_direct_valence)
         use_one_center = .NOT. has_pseudopotential
      CASE (skala_gapw_paw_one_center, skala_gapw_paw_one_center_split)
         CONTINUE
      CASE (skala_gapw_cp2k_default)
         IF (has_pseudopotential .AND. &
             ABS(zeff - REAL(zatom, dp)) <= 1.0E-10_dp) use_one_center = .FALSE.
      END SELECT

   END FUNCTION native_skala_uses_one_center_kind

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param energy_only ...
!> \param exc1 the on-body ex energy contribution
!> \param adiabatic_rescale_factor ...
!> \param kind_set_external provides a non-default kind_set to use
!> \param rho_atom_set_external provides a non-default atomic density set to use
!> \param xc_section_external provides an external non-default XC
!> \param calculate_forces ...
!> \param composite_vxc_rho ...
!> \param composite_vxc_tau ...
!> \param composite_reference_active ...
!> \param direct_valence_atom_grid evaluate the smooth valence fields on atom-centered grids
!> \param atom_composite_grid evaluate GAPW primitive fields on atom-centered composite grids
! **************************************************************************************************
   SUBROUTINE calculate_vxc_atom(qs_env, energy_only, exc1, &
                                 adiabatic_rescale_factor, kind_set_external, &
                                 rho_atom_set_external, xc_section_external, calculate_forces, &
                                 composite_vxc_rho, composite_vxc_tau, composite_reference_active, &
                                 direct_valence_atom_grid, atom_composite_grid)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: energy_only
      REAL(dp), INTENT(INOUT)                            :: exc1
      REAL(dp), INTENT(IN), OPTIONAL                     :: adiabatic_rescale_factor
      TYPE(qs_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: kind_set_external
      TYPE(rho_atom_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: rho_atom_set_external
      TYPE(section_vals_type), OPTIONAL, POINTER         :: xc_section_external
      LOGICAL, INTENT(IN), OPTIONAL                      :: calculate_forces
      TYPE(pw_r3d_rs_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: composite_vxc_rho, composite_vxc_tau
      LOGICAL, INTENT(OUT), OPTIONAL                     :: composite_reference_active
      LOGICAL, INTENT(IN), OPTIONAL                      :: direct_valence_atom_grid, &
                                                            atom_composite_grid

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_vxc_atom'
      INTEGER, PARAMETER :: lsd_grad_slots(10) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], &
         lsd_value_slots(6) = [5, 6, 7, 8, 9, 10], rks_grad_slots(5) = [1, 2, 3, 4, 5], &
         rks_value_slots(3) = [3, 4, 5]

      INTEGER :: adjoint_handle, adjoint_nchannels, adjoint_nn, adjoint_off, &
         atom_composite_components, bo(2), chunk_end_row, chunk_start_row, &
         composite_image_periodicity(3), composite_local_atom, composite_local_natom, &
         composite_nflat, composite_partition_target_image, composite_row, force_handle, &
         gapw_density_partition, gapw_representation, gather_handle, handle, ia, iat, iatom, &
         ichunk, icomponent, idir, ikind, image_i1, image_i2, image_i3, image_lower(3), &
         image_shell(3), image_shift(3), image_upper(3), inode_x, inode_y, inode_z, interp_handle, &
         ir, ispin, itask, iw, jdir
      INTEGER :: local_partition_target_image, myfun, na, natom, nchunks_per_atom, ngl, nr, &
         nrows_chunk, nspins, ntasks, nthread, num_pe, on, owner_rank, pos, row, source_atom, &
         target_atom, thread_id, xc_deriv_method_id, xc_rho_smooth_id, zatom
      INTEGER(int_8)                                     :: lin
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: composite_atomic_grid_sizes, &
                                                            composite_local_grid_sizes
      INTEGER, ALLOCATABLE, DIMENSION(:) :: chunk_atom, chunk_begin, chunk_end, &
         composite_atom_end, composite_atom_kind, composite_atom_kind_index, composite_atom_start, &
         composite_grid_atom, composite_local_atoms
      INTEGER, DIMENSION(2)                              :: owner_bounds
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL :: accint, atom_composite_active, atom_composite_diagnostic, &
         atom_composite_reference, comp_geometry_valid, composite_partition_layout_valid, &
         direct_valence_atom_composite, donlcc, evaluate_hard, evaluate_soft, gradient_f, lsd, &
         my_calculate_forces, native_grid_diagnostics, nlcc, one_center_kind, paw_atom, &
         paw_pseudopotentials, requested_atom_composite_grid, rho_g_valid, skala_atom_grid, &
         source_matrix_local, tau_f, tau_r_valid, use_atom_composite_density, &
         use_atom_composite_gradient, use_atom_composite_tau, use_virial
      REAL(dp) :: agr, alpha, atom_composite_exc, atom_composite_nelec, coefficient, &
         composite_cross_cutoff_max, composite_cross_density_max, composite_cross_grad_max, &
         composite_cross_kin_max, composite_density_max, composite_density_min, &
         composite_grad_max, composite_kin_max, composite_kin_min, composite_tau_integral, &
         cross_cutoff, density_cut, exc_h, exc_s, feature_vxc_analytic, feature_vxc_fd, &
         feature_vxc_minus, feature_vxc_plus, feature_vxc_step, gradient_cut, &
         local_partition_weight, my_adiabatic_rescale_factor, nlcc_density, nlcc_spin_factor, &
         one_center_density_field_contraction
      REAL(dp) :: one_center_density_matrix_contraction, one_center_field_contraction, &
         one_center_gradient_field_contraction, one_center_gradient_matrix_contraction, &
         one_center_matrix_contraction, one_center_rho_grad_field_contraction, &
         one_center_rho_grad_matrix_contraction, one_center_tau_field_contraction, &
         one_center_tau_matrix_contraction, one_center_tensor_contraction, partition_weight, &
         smooth_grid_contraction, smooth_input_contraction, target_partition_adjoint, tau_cut, zeff
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: adjoint_reduce_buf, &
         composite_atomic_grid_weight_grad, composite_base_grid_weights, &
         composite_grid_weight_grad, composite_grid_weights, fw
      REAL(dp), ALLOCATABLE, DIMENSION(:, :) :: adjoint_compact, adjoint_local_flat, &
         composite_atom_coord_grad, composite_atom_coords, composite_cross_density, &
         composite_cross_force, composite_cross_force_local, composite_cross_kin, &
         composite_density, composite_density_grad, composite_explicit_force, &
         composite_grid_coord_force, composite_grid_coord_force_local, composite_grid_coord_grad, &
         composite_grid_coords, composite_kin, composite_kin_grad, composite_local_atom_coords, &
         composite_model_atom_force, composite_moving_smooth_force, &
         composite_moving_smooth_force_local, composite_nlcc_center_force
      REAL(dp), ALLOCATABLE, DIMENSION(:, :) :: composite_nlcc_center_force_local, &
         composite_nlcc_target_force, composite_nlcc_target_force_local, &
         composite_partition_force, composite_partition_force_local, &
         composite_smooth_density_cache, composite_smooth_kin_cache, local_partition_datom
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), TARGET     :: smooth_density_adjoint_storage, &
                                                            smooth_kin_adjoint_storage
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :) :: composite_cross_grad, composite_grad, &
         composite_grad_grad, composite_int_h, composite_int_s, composite_smooth_gradient_cache, &
         vtau_h_local, vtau_s_local, vxc_h_local, vxc_s_local
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), TARGET  :: smooth_grad_adjoint_storage
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :, :)       :: vxg_h_local, vxg_s_local
      REAL(dp), DIMENSION(1, 1, 1)                       :: tau_d
      REAL(dp), DIMENSION(1, 1, 1, 1)                    :: rho_d
      REAL(dp), DIMENSION(10)                            :: composite_value_batch
      REAL(dp), DIMENSION(2) :: composite_smooth_density_adjoint_value, &
         composite_smooth_density_value, composite_smooth_kin_adjoint_value, &
         composite_smooth_kin_value, cross_density, cross_density_adjoint, cross_kin, &
         cross_kin_adjoint
      REAL(dp), DIMENSION(3) :: composite_point, cross_displacement, cross_spatial_derivative, &
         fractional, image_translation, nlcc_gradient, nlcc_spatial_derivative, &
         skala_atom_force_h, skala_atom_force_s, spatial_derivative
      REAL(dp), DIMENSION(3, 10)                         :: composite_grad_batch
      REAL(dp), DIMENSION(3, 2) :: composite_smooth_gradient_adjoint_value, &
         composite_smooth_gradient_value, cross_density_spatial, cross_grad, cross_grad_adjoint, &
         cross_kin_spatial
      REAL(dp), DIMENSION(3, 3) :: composite_cross_image_virial, &
         composite_cross_image_virial_local, composite_explicit_virial, composite_feature_virial, &
         composite_interpolation_virial, composite_partition_strain_virial, &
         local_partition_dstrain, nlcc_hessian, skala_atom_virial, skala_atom_virial_h, &
         skala_atom_virial_s
      REAL(dp), DIMENSION(3, 3, 2)                       :: cross_grad_spatial
      REAL(dp), DIMENSION(4)                             :: feature_component_analytic, &
                                                            feature_component_fd
      REAL(dp), DIMENSION(:, :), POINTER                 :: rho_nlcc, thread_adj, weight_h, weight_s
      REAL(dp), DIMENSION(:, :, :), POINTER :: rho_h, rho_s, smooth_rho, smooth_rhoa, smooth_rhob, &
         smooth_tau, smooth_tau_a, smooth_tau_b, tau_h, tau_s, vtau_h, vtau_s, vxc_h, vxc_s
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: drho_h, drho_s, vxg_h, vxg_s
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_3d_r_cp_type), DIMENSION(3)                :: smooth_drho, smooth_drhoa, smooth_drhob
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(native_grid_interpolation_stencil_type)       :: interpolation_stencil
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: smooth_rho_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: smooth_rho_r, smooth_tau_r, &
                                                            smooth_vxc_rho, smooth_vxc_tau
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: my_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: cpc_h, cpc_s, dr_h, dr_s, int_hh, &
                                                            int_ss, r_h, r_s
      TYPE(rho_atom_coeff), DIMENSION(:, :), POINTER     :: r_h_d, r_s_d
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: my_rho_atom_set
      TYPE(rho_atom_type), POINTER                       :: rho_atom
      TYPE(section_vals_type), POINTER                   :: input, model_section, my_xc_section, &
                                                            native_section, xc_fun_section
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential
      TYPE(skala_periodic_partition_type)                :: composite_partition_layout
      TYPE(tau_basis_cache_type)                         :: tau_basis_cache
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set_h, rho_set_s, smooth_rho_set

! -------------------------------------------------------------------------

      CALL timeset(routineN, handle)
      composite_partition_layout_valid = .FALSE.

      NULLIFY (atom_list)
      NULLIFY (auxbas_pw_pool)
      NULLIFY (my_kind_set)
      NULLIFY (atomic_kind_set)
      NULLIFY (cell)
      NULLIFY (grid_atom)
      NULLIFY (gth_potential)
      NULLIFY (force)
      NULLIFY (harmonics)
      NULLIFY (input)
      NULLIFY (para_env)
      NULLIFY (particle_set)
      NULLIFY (pw_env)
      NULLIFY (rho_atom)
      NULLIFY (rho_struct)
      NULLIFY (my_rho_atom_set)
      NULLIFY (rho_nlcc)
      NULLIFY (smooth_rho, smooth_rhoa, smooth_rhob, smooth_tau, smooth_tau_a, smooth_tau_b)
      NULLIFY (smooth_rho_g, smooth_rho_r, smooth_tau_r)
      NULLIFY (smooth_vxc_rho, smooth_vxc_tau)
      DO idir = 1, 3
         NULLIFY (smooth_drho(idir)%array, smooth_drhoa(idir)%array, smooth_drhob(idir)%array)
      END DO
      NULLIFY (sgp_potential)
      NULLIFY (virial)
      my_calculate_forces = .FALSE.
      IF (PRESENT(calculate_forces)) my_calculate_forces = calculate_forces
      IF (PRESENT(composite_reference_active)) composite_reference_active = .FALSE.
      direct_valence_atom_composite = .FALSE.
      IF (PRESENT(direct_valence_atom_grid)) THEN
         direct_valence_atom_composite = direct_valence_atom_grid
      END IF
      requested_atom_composite_grid = .FALSE.
      IF (PRESENT(atom_composite_grid)) requested_atom_composite_grid = atom_composite_grid

      IF (PRESENT(adiabatic_rescale_factor)) THEN
         my_adiabatic_rescale_factor = adiabatic_rescale_factor
      ELSE
         my_adiabatic_rescale_factor = 1.0_dp
      END IF

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      cell=cell, &
                      para_env=para_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=my_kind_set, &
                      input=input, &
                      particle_set=particle_set, &
                      pw_env=pw_env, &
                      virial=virial, &
                      rho_atom_set=my_rho_atom_set, &
                      force=force)

      IF (dft_control%qs_control%gapw_xc) THEN
         CALL get_qs_env(qs_env=qs_env, rho_xc=rho_struct)
      ELSE
         CALL get_qs_env(qs_env=qs_env, rho=rho_struct)
      END IF

      IF (PRESENT(kind_set_external)) my_kind_set => kind_set_external
      IF (PRESENT(rho_atom_set_external)) my_rho_atom_set => rho_atom_set_external

      nlcc = has_nlcc(my_kind_set)
      accint = dft_control%qs_control%gapw_control%accurate_xcint

      my_xc_section => section_vals_get_subs_vals(input, "DFT%XC")

      IF (PRESENT(xc_section_external)) my_xc_section => xc_section_external

      xc_fun_section => section_vals_get_subs_vals(my_xc_section, "XC_FUNCTIONAL")
      CALL section_vals_val_get(xc_fun_section, "_SECTION_PARAMETERS_", &
                                i_val=myfun)
      skala_atom_grid = xc_section_uses_gauxc_model(my_xc_section)
      gapw_representation = skala_gapw_cp2k_default
      atom_composite_diagnostic = .FALSE.
      atom_composite_reference = .FALSE.
      paw_pseudopotentials = .FALSE.
      native_grid_diagnostics = .FALSE.
      atom_composite_components = 1
      feature_vxc_step = 3.0E-3_dp
      IF (skala_atom_grid) THEN
         native_section => get_skala_native_section(my_xc_section)
         model_section => native_section
         IF (.NOT. ASSOCIATED(model_section)) model_section => get_gauxc_section(my_xc_section)
         CPASSERT(ASSOCIATED(model_section))
         CALL section_vals_val_get(model_section, "PSEUDOPOTENTIAL_GAPW_REPRESENTATION", &
                                   i_val=gapw_representation)
         IF (ASSOCIATED(native_section)) THEN
            CALL section_vals_val_get(native_section, "DIAGNOSTICS", &
                                      l_val=native_grid_diagnostics)
         END IF
      END IF
      IF (skala_atom_grid) THEN
         DO ikind = 1, SIZE(my_kind_set)
            NULLIFY (gth_potential, sgp_potential)
            CALL get_qs_kind(my_kind_set(ikind), paw_atom=paw_atom, &
                             gth_potential=gth_potential, sgp_potential=sgp_potential)
            paw_pseudopotentials = paw_pseudopotentials .OR. &
                                   (paw_atom .AND. (ASSOCIATED(gth_potential) .OR. &
                                                    ASSOCIATED(sgp_potential)))
         END DO
      END IF
      IF (skala_atom_grid .AND. xc_section_uses_native_skala_grid(my_xc_section)) THEN
         native_section => get_skala_native_section(my_xc_section)
         CPASSERT(ASSOCIATED(native_section))
         CALL section_vals_val_get(native_section, &
                                   "GAPW_ATOM_COMPOSITE_DIAGNOSTIC", &
                                   l_val=atom_composite_diagnostic)
         CALL section_vals_val_get(native_section, &
                                   "GAPW_ATOM_COMPOSITE_REFERENCE", &
                                   l_val=atom_composite_reference)
         CALL section_vals_val_get(native_section, &
                                   "GAPW_ATOM_COMPOSITE_COMPONENTS", &
                                   i_val=atom_composite_components)
         CALL section_vals_val_get(native_section, &
                                   "GAPW_ATOM_COMPOSITE_FD_STEP", &
                                   r_val=feature_vxc_step)
      END IF
      atom_composite_reference = atom_composite_reference .OR. &
                                 (gapw_representation == skala_gapw_paw_one_center .AND. &
                                  paw_pseudopotentials)
      atom_composite_reference = atom_composite_reference .OR. requested_atom_composite_grid
      atom_composite_reference = atom_composite_reference .OR. direct_valence_atom_composite
      atom_composite_active = atom_composite_diagnostic .OR. atom_composite_reference
      use_atom_composite_density = atom_composite_components <= 2
      use_atom_composite_gradient = atom_composite_components <= 2
      use_atom_composite_tau = atom_composite_components == 1 .OR. &
                               atom_composite_components == 3
      IF (atom_composite_active) THEN
         CALL ensure_native_skala_atom_grids(my_kind_set, dft_control)
      END IF
      IF (direct_valence_atom_composite) THEN
         use_atom_composite_density = .FALSE.
         use_atom_composite_gradient = .FALSE.
         use_atom_composite_tau = .FALSE.
      END IF
      ! CP2K's auxiliary PW fields are represented on one index-periodic cell even
      ! when the physical Poisson problem is isolated or partially periodic.
      composite_image_periodicity = 1
      IF (PRESENT(composite_reference_active)) composite_reference_active = atom_composite_reference
      gapw_density_partition = skala_gapw_density_partition_hard_minus_soft
      IF (skala_atom_grid) THEN
         gapw_density_partition = native_skala_gapw_density_partition(my_xc_section)
      END IF
      use_virial = ASSOCIATED(virial)
      IF (use_virial) use_virial = my_calculate_forces .AND. &
                                   virial%pv_calculate .AND. (.NOT. virial%pv_numer)

      IF (myfun == xc_none) THEN
         exc1 = 0.0_dp
         my_rho_atom_set(:)%exc_h = 0.0_dp
         my_rho_atom_set(:)%exc_s = 0.0_dp
      ELSE
         CALL section_vals_val_get(my_xc_section, "DENSITY_CUTOFF", &
                                   r_val=density_cut)
         CALL section_vals_val_get(my_xc_section, "GRADIENT_CUTOFF", &
                                   r_val=gradient_cut)
         CALL section_vals_val_get(my_xc_section, "TAU_CUTOFF", &
                                   r_val=tau_cut)

         lsd = dft_control%lsd
         nspins = dft_control%nspins
         needs = xc_functionals_get_needs(xc_fun_section, &
                                          lsd=lsd, &
                                          calc_potential=.TRUE.)

         gradient_f = (needs%drho .OR. needs%drho_spin) .OR. skala_atom_grid
         tau_f = (needs%tau .OR. needs%tau_spin) .OR. skala_atom_grid

         IF (atom_composite_active) THEN
            IF (lsd) THEN
               needs%rho_spin = .TRUE.
               needs%drho_spin = .TRUE.
               needs%tau_spin = .TRUE.
            ELSE
               needs%rho = .TRUE.
               needs%drho = .TRUE.
               needs%tau = .TRUE.
            END IF

            ALLOCATE (composite_atomic_grid_sizes(SIZE(particle_set)), &
                      composite_atom_kind(SIZE(particle_set)), &
                      composite_atom_kind_index(SIZE(particle_set)), &
                      composite_atom_start(SIZE(particle_set)), &
                      composite_atom_end(SIZE(particle_set)), &
                      composite_atom_coords(3, SIZE(particle_set)))
            composite_atomic_grid_sizes = 0_int_8
            composite_atom_kind = 0
            composite_atom_kind_index = 0
            composite_atom_start = 0
            composite_atom_end = 0
            DO iatom = 1, SIZE(particle_set)
               composite_atom_coords(:, iatom) = particle_set(iatom)%r
            END DO
            DO ikind = 1, SIZE(atomic_kind_set)
               CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
               NULLIFY (gth_potential, sgp_potential)
               CALL get_qs_kind(my_kind_set(ikind), paw_atom=paw_atom, &
                                gth_potential=gth_potential, grid_atom=grid_atom, &
                                sgp_potential=sgp_potential, zatom=zatom, zeff=zeff)
               DO iat = 1, natom
                  iatom = atom_list(iat)
                  composite_atomic_grid_sizes(iatom) = INT(grid_atom%nr*grid_atom%ng_sphere, KIND=int_8)
                  composite_atom_kind(iatom) = ikind
                  composite_atom_kind_index(iatom) = iat
               END DO
            END DO
            IF (ANY(composite_atomic_grid_sizes <= 0_int_8)) THEN
               CALL cp_abort(__LOCATION__, &
                             "The atom-composite diagnostic requires a GAPW one-center grid for every atom.")
            END IF

            composite_local_natom = 0
            DO ikind = 1, SIZE(atomic_kind_set)
               CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
               NULLIFY (gth_potential, sgp_potential)
               CALL get_qs_kind(my_kind_set(ikind), paw_atom=paw_atom, &
                                gth_potential=gth_potential, sgp_potential=sgp_potential, &
                                zatom=zatom, zeff=zeff)
               bo = get_limit(natom, para_env%num_pe, para_env%mepos)
               composite_local_natom = composite_local_natom + MAX(0, bo(2) - bo(1) + 1)
            END DO
            ALLOCATE (composite_local_atoms(composite_local_natom), &
                      composite_local_grid_sizes(composite_local_natom), &
                      composite_local_atom_coords(3, composite_local_natom))
            composite_local_atom = 0
            composite_nflat = 0
            DO ikind = 1, SIZE(atomic_kind_set)
               CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
               NULLIFY (gth_potential, sgp_potential)
               CALL get_qs_kind(my_kind_set(ikind), paw_atom=paw_atom, &
                                gth_potential=gth_potential, sgp_potential=sgp_potential, &
                                zatom=zatom, zeff=zeff)
               bo = get_limit(natom, para_env%num_pe, para_env%mepos)
               DO iat = bo(1), bo(2)
                  iatom = atom_list(iat)
                  composite_local_atom = composite_local_atom + 1
                  composite_local_atoms(composite_local_atom) = iatom
                  composite_local_grid_sizes(composite_local_atom) = &
                     composite_atomic_grid_sizes(iatom)
                  composite_local_atom_coords(:, composite_local_atom) = &
                     composite_atom_coords(:, iatom)
                  composite_atom_start(iatom) = composite_nflat + 1
                  composite_nflat = composite_nflat + INT(composite_atomic_grid_sizes(iatom))
                  composite_atom_end(iatom) = composite_nflat
               END DO
            END DO
            CPASSERT(composite_local_atom == composite_local_natom)
            ALLOCATE (composite_density(composite_nflat, 2), &
                      composite_grad(composite_nflat, 3, 2), &
                      composite_kin(composite_nflat, 2), &
                      composite_grid_atom(composite_nflat), &
                      composite_grid_coords(3, composite_nflat), &
                      composite_grid_weights(composite_nflat), &
                      composite_base_grid_weights(composite_nflat))
            IF (atom_composite_diagnostic) THEN
               ALLOCATE (composite_smooth_density_cache(composite_nflat, 2), &
                         composite_smooth_gradient_cache(composite_nflat, 3, 2), &
                         composite_smooth_kin_cache(composite_nflat, 2))
            END IF
            composite_density = 0.0_dp
            composite_grad = 0.0_dp
            composite_kin = 0.0_dp
            DO composite_local_atom = 1, composite_local_natom
               iatom = composite_local_atoms(composite_local_atom)
               composite_grid_atom(composite_atom_start(iatom):composite_atom_end(iatom)) = iatom
            END DO
            IF (atom_composite_diagnostic) THEN
               composite_smooth_density_cache = 0.0_dp
               composite_smooth_gradient_cache = 0.0_dp
               composite_smooth_kin_cache = 0.0_dp
            END IF
            composite_grid_coords = 0.0_dp
            composite_grid_weights = 0.0_dp
            composite_base_grid_weights = 0.0_dp

            CALL qs_rho_get(rho_struct, rho_r=smooth_rho_r, rho_g=smooth_rho_g, &
                            tau_r=smooth_tau_r, rho_g_valid=rho_g_valid, &
                            tau_r_valid=tau_r_valid)
            CPASSERT(rho_g_valid)
            CPASSERT(tau_r_valid)
            CPASSERT(ASSOCIATED(smooth_rho_r))
            CPASSERT(ASSOCIATED(smooth_rho_g))
            CPASSERT(ASSOCIATED(smooth_tau_r))
            CPASSERT(smooth_rho_r(1)%pw_grid%para%group%num_pe == para_env%num_pe)
            CPASSERT(smooth_rho_r(1)%pw_grid%para%group%mepos == para_env%mepos)
            CALL prepare_native_grid_cache(skala_state%native_grid_cache, smooth_rho_r(1)%pw_grid, &
                                           cell, composite_nflat, .TRUE.)
            CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
            CALL section_vals_val_get(my_xc_section, "XC_GRID%XC_DERIV", &
                                      i_val=xc_deriv_method_id)
            CALL section_vals_val_get(my_xc_section, "XC_GRID%XC_SMOOTH_RHO", &
                                      i_val=xc_rho_smooth_id)
            CALL xc_rho_set_create(smooth_rho_set, smooth_rho_r(1)%pw_grid%bounds_local, &
                                   rho_cutoff=section_get_rval(my_xc_section, "density_cutoff"), &
                                   drho_cutoff=section_get_rval(my_xc_section, "gradient_cutoff"), &
                                   tau_cutoff=section_get_rval(my_xc_section, "tau_cutoff"))
            CALL xc_rho_set_update(smooth_rho_set, smooth_rho_r, smooth_rho_g, smooth_tau_r, needs, &
                                   xc_deriv_method_id, xc_rho_smooth_id, auxbas_pw_pool)
            CALL timeset("calculate_vxc_atom_gather", gather_handle)
  comp_geometry_valid = skala_state%native_atom%native_comp_ready .AND. skala_state%native_atom%native_comp_nrows == composite_nflat
            IF (comp_geometry_valid) THEN
               comp_geometry_valid = native_geometry_matches(skala_state%native_atom%native_comp_geometry, cell, &
                                                             composite_atom_coords)
               ! The cached coordinates cover this rank's local atoms only, so a
               ! change in the local partition must invalidate them too.
               comp_geometry_valid = comp_geometry_valid .AND. ALLOCATED(skala_state%native_atom%native_comp_local_atoms)
               IF (comp_geometry_valid) THEN
                  comp_geometry_valid = SIZE(skala_state%native_atom%native_comp_local_atoms) == composite_local_natom
               END IF
               IF (comp_geometry_valid) THEN
                  comp_geometry_valid = ALL(skala_state%native_atom%native_comp_local_atoms == composite_local_atoms)
               END IF
            ELSE
               comp_geometry_valid = .FALSE.
            END IF
            IF (comp_geometry_valid) THEN
               composite_grid_coords(:, :) = skala_state%native_atom%native_comp_coords
               composite_base_grid_weights(:) = skala_state%native_atom%native_comp_base_w
               composite_grid_weights(:) = skala_state%native_atom%native_comp_grid_w
            ELSE
               IF (ALLOCATED(skala_state%native_atom%native_comp_coords)) THEN
                  DEALLOCATE (skala_state%native_atom%native_comp_coords, &
                              skala_state%native_atom%native_comp_base_w, &
                              skala_state%native_atom%native_comp_grid_w)
               END IF
        IF (ALLOCATED(skala_state%native_atom%native_comp_local_atoms)) DEALLOCATE (skala_state%native_atom%native_comp_local_atoms)
               ALLOCATE (skala_state%native_atom%native_comp_coords(3, composite_nflat), &
                         skala_state%native_atom%native_comp_base_w(composite_nflat), &
                         skala_state%native_atom%native_comp_grid_w(composite_nflat), &
                         skala_state%native_atom%native_comp_local_atoms(SIZE(composite_local_atoms)))
               ! The periodic image layout only depends on the geometry, not on
               ! the target atom, so build it once and share it across threads.
               CALL build_periodic_atom_image_layout( &
                  composite_atom_coords, cell, 1, composite_image_periodicity, &
                  composite_partition_layout, composite_partition_target_image)
               composite_partition_layout_valid = .TRUE.
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(NONE) &
!$OMP PRIVATE(composite_local_atom, iatom, grid_atom, nr, na, &
!$OMP         composite_partition_target_image, ir, ia, composite_row, composite_point, &
!$OMP         partition_weight) &
!$OMP SHARED(cell, composite_atom_coords, composite_atom_end, composite_atom_kind, &
!$OMP        composite_atom_start, composite_base_grid_weights, composite_grid_coords, &
!$OMP        composite_grid_weights, composite_image_periodicity, composite_local_atoms, &
!$OMP        composite_local_natom, composite_partition_layout, my_kind_set, particle_set)
               DO composite_local_atom = 1, composite_local_natom
                  iatom = composite_local_atoms(composite_local_atom)
                  composite_partition_target_image = iatom
                  CALL get_qs_kind(my_kind_set(composite_atom_kind(iatom)), grid_atom=grid_atom)
                  nr = grid_atom%nr
                  na = grid_atom%ng_sphere
                  DO ir = 1, nr
                     DO ia = 1, na
                        composite_row = composite_atom_start(iatom) + (ir - 1)*na + ia - 1
                        composite_point(1) = particle_set(iatom)%r(1) + grid_atom%rad(ir)* &
                                             grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia)
                        composite_point(2) = particle_set(iatom)%r(2) + grid_atom%rad(ir)* &
                                             grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia)
                        composite_point(3) = particle_set(iatom)%r(3) + &
                                             grid_atom%rad(ir)*grid_atom%cos_pol(ia)
                        CALL periodic_atom_image_partition_from_layout( &
                           composite_point, composite_partition_layout, &
                           composite_partition_target_image, partition_weight)
                        composite_grid_coords(:, composite_row) = composite_point
                        composite_base_grid_weights(composite_row) = grid_atom%weight(ia, ir)
                        composite_grid_weights(composite_row) = &
                           composite_base_grid_weights(composite_row)*partition_weight
                     END DO
                  END DO
               END DO
!$OMP END PARALLEL DO
               skala_state%native_atom%native_comp_coords(:, :) = composite_grid_coords
               skala_state%native_atom%native_comp_base_w(:) = composite_base_grid_weights
               skala_state%native_atom%native_comp_grid_w(:) = composite_grid_weights
               CALL store_native_geometry_key(skala_state%native_atom%native_comp_geometry, cell, composite_atom_coords)
               skala_state%native_atom%native_comp_local_atoms(:) = composite_local_atoms
               skala_state%native_atom%native_comp_nrows = composite_nflat
               skala_state%native_atom%native_comp_ready = .TRUE.
            END IF

            CALL native_grid_build_needed(smooth_rho_r(1)%pw_grid, cell, smooth_rho_r(1)%pw_grid%para%group, &
                                          .TRUE., composite_nflat, &
                                          composite_grid_coords, composite_atom_coords)
            IF (lsd) THEN
               CALL xc_rho_set_get(smooth_rho_set, rhoa=smooth_rhoa, rhob=smooth_rhob, &
                                   drhoa=smooth_drhoa, drhob=smooth_drhob, &
                                   tau_a=smooth_tau_a, tau_b=smooth_tau_b)
               CALL native_grid_alloc_fetch(10)
               CALL native_grid_extract_points(smooth_rho_r(1)%pw_grid, 1, 10, smooth_rhoa)
               CALL native_grid_extract_points(smooth_rho_r(1)%pw_grid, 2, 10, smooth_rhob)
               CALL native_grid_extract_points(smooth_rho_r(1)%pw_grid, 3, 10, smooth_tau_a)
               CALL native_grid_extract_points(smooth_rho_r(1)%pw_grid, 4, 10, smooth_tau_b)
               DO idir = 1, 3
                  CALL native_grid_extract_points(smooth_rho_r(1)%pw_grid, 4 + idir, 10, &
                                                  smooth_drhoa(idir)%array)
                  CALL native_grid_extract_points(smooth_rho_r(1)%pw_grid, 7 + idir, 10, &
                                                  smooth_drhob(idir)%array)
               END DO
            ELSE
               CALL xc_rho_set_get(smooth_rho_set, rho=smooth_rho, drho=smooth_drho, &
                                   tau=smooth_tau)
               CALL native_grid_alloc_fetch(5)
               CALL native_grid_extract_points(smooth_rho_r(1)%pw_grid, 1, 5, smooth_rho)
               CALL native_grid_extract_points(smooth_rho_r(1)%pw_grid, 2, 5, smooth_tau)
               DO idir = 1, 3
                  CALL native_grid_extract_points(smooth_rho_r(1)%pw_grid, 2 + idir, 5, &
                                                  smooth_drho(idir)%array)
               END DO
            END IF
            CALL native_grid_exchange_fetch(smooth_rho_r(1)%pw_grid%para%group)
            CALL timestop(gather_handle)
         END IF

         ! Initialize energy contribution from the one center XC terms to zero
         exc1 = 0.0_dp

         ! Nullify some pointers for work-arrays
         NULLIFY (rho_h, drho_h, rho_s, drho_s, weight_h, weight_s)
         NULLIFY (vxc_h, vxc_s, vxg_h, vxg_s)
         NULLIFY (tau_h, tau_s)
         NULLIFY (vtau_h, vtau_s)

         ! Here starts the loop over all the atoms

         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
            NULLIFY (gth_potential, sgp_potential)
            CALL get_qs_kind(my_kind_set(ikind), paw_atom=paw_atom, &
                             gth_potential=gth_potential, harmonics=harmonics, &
                             grid_atom=grid_atom, sgp_potential=sgp_potential, &
                             zatom=zatom, zeff=zeff)
            one_center_kind = .NOT. direct_valence_atom_composite
            IF (one_center_kind) THEN
               CALL get_qs_kind(my_kind_set(ikind), basis_set=basis_1c, basis_type="GAPW_1C")
               one_center_kind = paw_atom
               IF (skala_atom_grid) THEN
                  one_center_kind = native_skala_uses_one_center_kind( &
                                    paw_atom, gapw_representation, &
                                    ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential), &
                                    zeff, zatom)
               END IF
            END IF
            IF (.NOT. one_center_kind .AND. .NOT. atom_composite_active) CYCLE

            nr = grid_atom%nr
            na = grid_atom%ng_sphere

            IF (one_center_kind) THEN
               ! Prepare the structures needed to calculate and store the one-center XC derivatives.

               ! Array dimension: here anly one dimensional arrays are used,
               ! i.e. only the first column of deriv_data is read.
               ! The other to dimensions  are set to size equal 1
               bounds(1:2, 1:3) = 1
               bounds(2, 1) = na
               bounds(2, 2) = nr

               ! set integration weights
               weight_h => grid_atom%weight
               weight_s => grid_atom%weight
               IF (accint) THEN
                  on = dft_control%qs_control%gapw_control%oweights
                  alpha = dft_control%qs_control%gapw_control%aw(ikind)
                  IF (ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                     IF (grid_atom%gapw_weight_alpha /= alpha) DEALLOCATE (grid_atom%gapw_weight_s)
                  END IF
                  IF (.NOT. ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                     ALLOCATE (grid_atom%gapw_weight_s(na, nr))
                     ALLOCATE (fw(nr))
                     CALL calc_weight_function(fw, grid_atom%rad2, on, alpha)
                     DO ir = 1, nr
                        agr = 1.0_dp - fw(ir)
                        grid_atom%gapw_weight_s(:, ir) = agr*grid_atom%weight(:, ir)
                     END DO
                     DEALLOCATE (fw)
                     grid_atom%gapw_weight_alpha = alpha
                  END IF
                  weight_s => grid_atom%gapw_weight_s
               END IF

               IF (.NOT. skala_atom_grid) THEN
                  CALL xc_dset_create(deriv_set, local_bounds=bounds)
                  CALL xc_rho_set_create(rho_set_h, bounds, rho_cutoff=density_cut, &
                                         drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
                  CALL xc_rho_set_create(rho_set_s, bounds, rho_cutoff=density_cut, &
                                         drho_cutoff=gradient_cut, tau_cutoff=tau_cut)

                  CALL xc_rho_set_atom_update(rho_set_h, needs, nspins, bounds)
                  CALL xc_rho_set_atom_update(rho_set_s, needs, nspins, bounds)
               END IF

               CALL reallocate(rho_h, 1, na, 1, nr, 1, nspins)
               CALL reallocate(rho_s, 1, na, 1, nr, 1, nspins)
               CALL reallocate(vxc_h, 1, na, 1, nr, 1, nspins)
               CALL reallocate(vxc_s, 1, na, 1, nr, 1, nspins)
               !
               IF (gradient_f) THEN
                  CALL reallocate(drho_h, 1, 4, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(drho_s, 1, 4, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(vxg_h, 1, 3, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(vxg_s, 1, 3, 1, na, 1, nr, 1, nspins)
               END IF

               IF (tau_f) THEN
                  CALL create_tau_basis_cache(tau_basis_cache, grid_atom, basis_1c, harmonics)
                  CALL reallocate(tau_h, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(tau_s, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(vtau_h, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(vtau_s, 1, na, 1, nr, 1, nspins)
               END IF

               ! NLCC for separate hard and soft one-center densities.
               donlcc = .FALSE.
               IF (nlcc) THEN
                  NULLIFY (rho_nlcc)
                  rho_nlcc => my_kind_set(ikind)%nlcc_pot
                  IF (ASSOCIATED(rho_nlcc)) donlcc = .TRUE.
               END IF
            END IF

            ! Distribute the atoms of this kind

            num_pe = para_env%num_pe
            bo = get_limit(natom, para_env%num_pe, para_env%mepos)

            DO iat = bo(1), bo(2)
               iatom = atom_list(iat)

               IF (one_center_kind) THEN
                  my_rho_atom_set(iatom)%exc_h = 0.0_dp
                  my_rho_atom_set(iatom)%exc_s = 0.0_dp

                  rho_atom => my_rho_atom_set(iatom)
                  rho_h = 0.0_dp
                  rho_s = 0.0_dp
                  IF (gradient_f) THEN
                     NULLIFY (r_h, r_s, dr_h, dr_s, r_h_d, r_s_d)
                     CALL get_rho_atom(rho_atom=rho_atom, rho_rad_h=r_h, &
                                       rho_rad_s=r_s, drho_rad_h=dr_h, &
                                       drho_rad_s=dr_s, rho_rad_h_d=r_h_d, &
                                       rho_rad_s_d=r_s_d)
                     drho_h = 0.0_dp
                     drho_s = 0.0_dp
                  ELSE
                     NULLIFY (r_h, r_s)
                     CALL get_rho_atom(rho_atom=rho_atom, rho_rad_h=r_h, rho_rad_s=r_s)
                     rho_d = 0.0_dp
                  END IF
                  IF (tau_f) THEN
                     CALL calc_tau_atom(tau_h, tau_s, rho_atom, tau_basis_cache, nspins)
                  ELSE
                     tau_d = 0.0_dp
                  END IF

                  DO ir = 1, nr
                     CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_f, &
                                           ir, r_h, r_s, rho_h, rho_s, dr_h, dr_s, &
                                           r_h_d, r_s_d, drho_h, drho_s)
                     IF (donlcc) THEN
                        CALL calc_rho_nlcc(grid_atom, nspins, gradient_f, &
                                           ir, rho_nlcc(:, 1), rho_h, rho_s, &
                                           rho_nlcc(:, 2), drho_h, drho_s)
                     END IF
                  END DO
               END IF

               IF (atom_composite_active) THEN
                  CALL timeset("calculate_vxc_atom_interp", interp_handle)
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP PRIVATE(composite_row, composite_point, composite_smooth_density_value, &
!$OMP         composite_smooth_gradient_value, composite_smooth_kin_value, &
!$OMP         idir, ispin, &
!$OMP         gth_potential, nlcc_density, nlcc_gradient, nlcc_hessian, nlcc_spin_factor, &
!$OMP         interpolation_stencil, sgp_potential, source_atom) &
!$OMP SHARED(atom_composite_reference, cell, composite_atom_kind, &
!$OMP        composite_atom_start, &
!$OMP        composite_density, composite_grad, &
!$OMP        composite_grid_coords, composite_kin, &
!$OMP        composite_smooth_density_cache, composite_smooth_gradient_cache, &
!$OMP        composite_smooth_kin_cache, drho_h, drho_s, grid_atom, iatom, &
!$OMP        lsd, my_kind_set, na, nlcc, nr, one_center_kind, &
!$OMP        particle_set, qs_env, rho_h, rho_s, smooth_rho_r, tau_h, tau_s, &
!$OMP        atom_composite_diagnostic, use_atom_composite_density, &
!$OMP        use_atom_composite_gradient, use_atom_composite_tau)
                  DO ir = 1, nr
                     DO ia = 1, na
                        composite_row = composite_atom_start(iatom) + (ir - 1)*na + ia - 1
                        composite_point = composite_grid_coords(:, composite_row)
                        CALL get_native_grid_stencil(smooth_rho_r(1)%pw_grid, cell, &
                                                     composite_row, composite_point, &
                                                     interpolation_stencil, .TRUE.)
                        IF (lsd) THEN
                           CALL interpolate_native_grid_fields_compact_lsd( &
                              interpolation_stencil, composite_smooth_density_value, &
                              composite_smooth_gradient_value, composite_smooth_kin_value)
                           IF (atom_composite_diagnostic) THEN
                              composite_smooth_density_cache(composite_row, :) = &
                                 composite_smooth_density_value
                              composite_smooth_gradient_cache(composite_row, :, :) = &
                                 composite_smooth_gradient_value
                              composite_smooth_kin_cache(composite_row, :) = &
                                 composite_smooth_kin_value
                           END IF
                           DO ispin = 1, 2
                              composite_density(composite_row, ispin) = &
                                 composite_smooth_density_value(ispin)
                              composite_grad(composite_row, :, ispin) = &
                                 composite_smooth_gradient_value(:, ispin)
                              composite_kin(composite_row, ispin) = &
                                 composite_smooth_kin_value(ispin)
                              IF (one_center_kind .AND. use_atom_composite_density) THEN
                                 composite_density(composite_row, ispin) = &
                                    composite_density(composite_row, ispin) + &
                                    rho_h(ia, ir, ispin) - rho_s(ia, ir, ispin)
                              END IF
                              IF (one_center_kind .AND. use_atom_composite_gradient) THEN
                                 DO idir = 1, 3
                                    composite_grad(composite_row, idir, ispin) = &
                                       composite_grad(composite_row, idir, ispin) + &
                                       drho_h(idir, ia, ir, ispin) - drho_s(idir, ia, ir, ispin)
                                 END DO
                              END IF
                              IF (one_center_kind .AND. use_atom_composite_tau) THEN
                                 composite_kin(composite_row, ispin) = &
                                    composite_kin(composite_row, ispin) + &
                                    tau_h(ia, ir, ispin) - tau_s(ia, ir, ispin)
                              END IF
                           END DO
                        ELSE
                           CALL interpolate_native_grid_fields_compact( &
                              1, 3, 4, 5, 2, interpolation_stencil, &
                              composite_smooth_density_value(1), &
                              composite_smooth_gradient_value(:, 1), &
                              composite_smooth_kin_value(1))
                           IF (atom_composite_diagnostic) THEN
                              composite_smooth_density_cache(composite_row, 1) = &
                                 composite_smooth_density_value(1)
                              composite_smooth_gradient_cache(composite_row, :, 1) = &
                                 composite_smooth_gradient_value(:, 1)
                              composite_smooth_kin_cache(composite_row, 1) = &
                                 composite_smooth_kin_value(1)
                           END IF
                           composite_density(composite_row, :) = &
                              0.5_dp*composite_smooth_density_value(1)
                           DO idir = 1, 3
                              composite_grad(composite_row, idir, :) = &
                                 0.5_dp*composite_smooth_gradient_value(idir, 1)
                           END DO
                           composite_kin(composite_row, :) = &
                              0.5_dp*composite_smooth_kin_value(1)
                           IF (one_center_kind .AND. use_atom_composite_density) THEN
                              composite_density(composite_row, :) = &
                                 composite_density(composite_row, :) + &
                                 0.5_dp*(rho_h(ia, ir, 1) - rho_s(ia, ir, 1))
                           END IF
                           IF (one_center_kind .AND. use_atom_composite_gradient) THEN
                              DO idir = 1, 3
                                 composite_grad(composite_row, idir, :) = &
                                    composite_grad(composite_row, idir, :) + &
                                    0.5_dp*(drho_h(idir, ia, ir, 1) - &
                                            drho_s(idir, ia, ir, 1))
                              END DO
                           END IF
                           IF (one_center_kind .AND. use_atom_composite_tau) THEN
                              composite_kin(composite_row, :) = composite_kin(composite_row, :) + &
                                                                0.5_dp*(tau_h(ia, ir, 1) - &
                                                                        tau_s(ia, ir, 1))
                           END IF
                        END IF
                        IF (atom_composite_reference .AND. nlcc) THEN
                           nlcc_spin_factor = 0.5_dp
                           DO source_atom = 1, SIZE(particle_set)
                              NULLIFY (gth_potential, sgp_potential)
                              CALL get_qs_kind(my_kind_set(composite_atom_kind(source_atom)), &
                                               gth_potential=gth_potential, &
                                               sgp_potential=sgp_potential)
                              CALL evaluate_nlcc_primitive_fields( &
                                 composite_point, particle_set(source_atom)%r, &
                                 gth_potential, sgp_potential, nlcc_density, &
                                 nlcc_gradient, nlcc_hessian)
                              composite_density(composite_row, :) = &
                                 composite_density(composite_row, :) + &
                                 nlcc_spin_factor*nlcc_density
                              DO idir = 1, 3
                                 composite_grad(composite_row, idir, :) = &
                                    composite_grad(composite_row, idir, :) + &
                                    nlcc_spin_factor*nlcc_gradient(idir)
                              END DO
                           END DO
                        END IF
                     END DO
                  END DO
!$OMP END PARALLEL DO
                  CALL timestop(interp_handle)
                  CPASSERT(nr*na == composite_atom_end(iatom) - composite_atom_start(iatom) + 1)
               END IF

               IF (.NOT. one_center_kind) CYCLE

               IF (.NOT. skala_atom_grid) THEN
                  DO ir = 1, nr
                     IF (tau_f) THEN
                        CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, drho_h, tau_h, na, ir)
                        CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, drho_s, tau_s, na, ir)
                     ELSE IF (gradient_f) THEN
                        CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, drho_h, tau_d, na, ir)
                        CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, drho_s, tau_d, na, ir)
                     ELSE
                        CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, rho_d, tau_d, na, ir)
                        CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, rho_d, tau_d, na, ir)
                     END IF
                  END DO
               END IF

               evaluate_hard = .TRUE.
               evaluate_soft = .TRUE.
               skala_atom_force_h = 0.0_dp
               skala_atom_force_s = 0.0_dp
               skala_atom_virial_h = 0.0_dp
               skala_atom_virial_s = 0.0_dp
               IF (skala_atom_grid) THEN
                  SELECT CASE (gapw_density_partition)
                  CASE (skala_gapw_density_partition_hard_minus_soft)
                     CONTINUE
                  CASE (skala_gapw_density_partition_hard_only)
                     evaluate_soft = .FALSE.
                  CASE (skala_gapw_density_partition_soft_only)
                     evaluate_hard = .FALSE.
                  CASE (skala_gapw_density_partition_none)
                     evaluate_hard = .FALSE.
                     evaluate_soft = .FALSE.
                  CASE DEFAULT
                     CALL cp_abort(__LOCATION__, &
                                   "Unknown SKALA_NATIVE%GAPW_DENSITY_PARTITION value.")
                  END SELECT
               END IF
               IF (atom_composite_reference) THEN
                  evaluate_hard = .FALSE.
                  evaluate_soft = .FALSE.
               END IF

               IF (skala_atom_grid) THEN
                  IF (.NOT. evaluate_hard) THEN
                     exc_h = 0.0_dp
                     IF (.NOT. energy_only) THEN
                        vxc_h = 0.0_dp
                        IF (ASSOCIATED(vxg_h)) vxg_h = 0.0_dp
                        IF (ASSOCIATED(vtau_h)) vtau_h = 0.0_dp
                     END IF
                  END IF
                  IF (.NOT. evaluate_soft) THEN
                     exc_s = 0.0_dp
                     IF (.NOT. energy_only) THEN
                        vxc_s = 0.0_dp
                        IF (ASSOCIATED(vxg_s)) vxg_s = 0.0_dp
                        IF (ASSOCIATED(vtau_s)) vtau_s = 0.0_dp
                     END IF
                  END IF
                  IF (evaluate_hard .OR. evaluate_soft) THEN
                     CALL skala_gapw_atom_vxc_of_r_pair( &
                        my_xc_section, grid_atom, para_env, particle_set(iatom)%r, &
                        rho_h, drho_h, tau_h, rho_s, drho_s, tau_s, &
                        weight_h, weight_s, lsd, nspins, na, nr, &
                        evaluate_hard, evaluate_soft, &
                        exc_h, vxc_h, vxg_h, vtau_h, &
                        exc_s, vxc_s, vxg_s, vtau_s, &
                        energy_only=energy_only, &
                        atom_force_h=skala_atom_force_h, atom_virial_h=skala_atom_virial_h, &
                        atom_force_s=skala_atom_force_s, atom_virial_s=skala_atom_virial_s)
                  END IF
                  rho_atom%exc_h = rho_atom%exc_h + exc_h
                  rho_atom%exc_s = rho_atom%exc_s + exc_s
               ELSE
                  CALL xc_dset_zero_all(deriv_set)
                  IF (.NOT. evaluate_hard) THEN
                     exc_h = 0.0_dp
                     IF (.NOT. energy_only) THEN
                        vxc_h = 0.0_dp
                        IF (ASSOCIATED(vxg_h)) vxg_h = 0.0_dp
                        IF (ASSOCIATED(vtau_h)) vtau_h = 0.0_dp
                     END IF
                  ELSE
                     CALL vxc_of_r_new(xc_fun_section, rho_set_h, deriv_set, 1, needs, weight_h, &
                                       lsd, na, nr, exc_h, vxc_h, vxg_h, vtau_h, energy_only=energy_only, &
                                       adiabatic_rescale_factor=my_adiabatic_rescale_factor)
                  END IF
                  rho_atom%exc_h = rho_atom%exc_h + exc_h

                  CALL xc_dset_zero_all(deriv_set)
                  IF (.NOT. evaluate_soft) THEN
                     exc_s = 0.0_dp
                     IF (.NOT. energy_only) THEN
                        vxc_s = 0.0_dp
                        IF (ASSOCIATED(vxg_s)) vxg_s = 0.0_dp
                        IF (ASSOCIATED(vtau_s)) vtau_s = 0.0_dp
                     END IF
                  ELSE
                     CALL vxc_of_r_new(xc_fun_section, rho_set_s, deriv_set, 1, needs, weight_s, &
                                       lsd, na, nr, exc_s, vxc_s, vxg_s, vtau_s, energy_only=energy_only, &
                                       adiabatic_rescale_factor=my_adiabatic_rescale_factor)
                  END IF
                  rho_atom%exc_s = rho_atom%exc_s + exc_s
               END IF

               ! Add contributions to the exc energy

               exc1 = exc1 + rho_atom%exc_h - rho_atom%exc_s
               IF (skala_atom_grid .AND. my_calculate_forces .AND. ASSOCIATED(force)) THEN
                  force(ikind)%rho_elec(:, iat) = force(ikind)%rho_elec(:, iat) + &
                                                  skala_atom_force_h - skala_atom_force_s
               END IF
               IF (skala_atom_grid .AND. use_virial) THEN
                  skala_atom_virial = skala_atom_virial_h - skala_atom_virial_s
                  DO idir = 1, 3
                     DO jdir = 1, 3
                        virial%pv_gapw(idir, jdir) = virial%pv_gapw(idir, jdir) + &
                                                     skala_atom_virial(idir, jdir)
                        virial%pv_virial(idir, jdir) = virial%pv_virial(idir, jdir) + &
                                                       skala_atom_virial(idir, jdir)
                     END DO
                  END DO
               END IF

               ! Integration to get the matrix elements relative to the vxc_atom
               ! here the products with the primitives is done: gaVxcgb
               ! internal transformation to get the integral in cartesian Gaussians

               IF (.NOT. energy_only) THEN
                  NULLIFY (int_hh, int_ss)
                  CALL get_rho_atom(rho_atom=rho_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
                  IF (gradient_f) THEN
                     CALL gaVxcgb_GC(vxc_h, vxc_s, vxg_h, vxg_s, int_hh, int_ss, &
                                     grid_atom, basis_1c, harmonics, nspins)
                  ELSE
                     CALL gaVxcgb_noGC(vxc_h, vxc_s, int_hh, int_ss, &
                                       grid_atom, basis_1c, harmonics, nspins)
                  END IF
                  IF (tau_f) THEN
                     CALL dgaVtaudgb(vtau_h, vtau_s, int_hh, int_ss, tau_basis_cache, nspins)
                  END IF
               END IF ! energy_only
               NULLIFY (r_h, r_s, dr_h, dr_s)
            END DO ! iat

            IF (one_center_kind) THEN
               IF (tau_f) CALL release_tau_basis_cache(tau_basis_cache)

               IF (.NOT. skala_atom_grid) THEN
                  CALL xc_dset_release(deriv_set)
                  CALL xc_rho_set_release(rho_set_h)
                  CALL xc_rho_set_release(rho_set_s)
               END IF
            END IF
         END DO ! ikind

         IF (atom_composite_active) THEN
            ALLOCATE (composite_cross_density(composite_nflat, 2), &
                      composite_cross_grad(composite_nflat, 3, 2), &
                      composite_cross_kin(composite_nflat, 2))
            composite_cross_density = 0.0_dp
            composite_cross_grad = 0.0_dp
            composite_cross_kin = 0.0_dp
            composite_cross_cutoff_max = 0.0_dp
            IF (.NOT. direct_valence_atom_composite) THEN
               DO ikind = 1, SIZE(atomic_kind_set)
                  CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
                  NULLIFY (gth_potential, sgp_potential)
                  CALL get_qs_kind(my_kind_set(ikind), paw_atom=paw_atom, &
                                   gth_potential=gth_potential, harmonics=harmonics, &
                                   grid_atom=grid_atom, sgp_potential=sgp_potential, &
                                   zatom=zatom, zeff=zeff)
                  CALL get_qs_kind(my_kind_set(ikind), basis_set=basis_1c, basis_type="GAPW_1C")
                  IF (.NOT. native_skala_uses_one_center_kind( &
                      paw_atom, gapw_representation, &
                      ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential), &
                      zeff, zatom)) CYCLE

                  CALL replicate_rho_atom_radial( &
                     para_env, my_rho_atom_set, my_kind_set(ikind), atom_list, natom, nspins)
                  nr = grid_atom%nr
                  na = grid_atom%ng_sphere
                  CALL create_tau_basis_cache(tau_basis_cache, grid_atom, basis_1c, harmonics)
                  CALL reallocate(rho_h, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(rho_s, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(drho_h, 1, 4, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(drho_s, 1, 4, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(tau_h, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(tau_s, 1, na, 1, nr, 1, nspins)

                  ! The one-center density matrices are already globally reduced. Distribute the
                  ! overlap work by target atom so that every rank constructs only its model rows.
                  DO iat = 1, natom
                     source_atom = atom_list(iat)
                     rho_atom => my_rho_atom_set(source_atom)
                     NULLIFY (r_h, r_s, dr_h, dr_s, r_h_d, r_s_d)
                     CALL get_rho_atom(rho_atom=rho_atom, rho_rad_h=r_h, rho_rad_s=r_s, &
                                       drho_rad_h=dr_h, drho_rad_s=dr_s, &
                                       rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d)
                     rho_h = 0.0_dp
                     rho_s = 0.0_dp
                     drho_h = 0.0_dp
                     drho_s = 0.0_dp
                     CALL calc_tau_atom(tau_h, tau_s, rho_atom, tau_basis_cache, nspins)
                     DO ir = 1, nr
                        CALL calc_rho_angular(grid_atom, harmonics, nspins, .TRUE., &
                                              ir, r_h, r_s, rho_h, rho_s, dr_h, dr_s, &
                                              r_h_d, r_s_d, drho_h, drho_s)
                     END DO

                     cross_cutoff = gapw_atom_grid_support_radius( &
                                    grid_atom, rho_h, rho_s, drho_h, drho_s, tau_h, tau_s)
                     IF (cross_cutoff <= 0.0_dp) CYCLE
                     composite_cross_cutoff_max = MAX(composite_cross_cutoff_max, cross_cutoff)
                     image_shell = 0
                     DO idir = 1, 3
                        IF (cell%perd(idir) == 1) THEN
                           image_shell(idir) = CEILING( &
                                               cross_cutoff*SQRT(SUM(cell%h_inv(idir, :)**2))) + 1
                        END IF
                     END DO

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) &
!$OMP PRIVATE(image_lower, image_upper, composite_row, cross_density, cross_density_spatial, &
!$OMP         cross_displacement, cross_grad, cross_grad_spatial, cross_kin, &
!$OMP         cross_kin_spatial, fractional, idir, image_i1, image_i2, image_i3, jdir, &
!$OMP         image_shift, image_translation, target_atom) &
!$OMP SHARED(cell, composite_cross_density, composite_cross_grad, composite_cross_kin, &
!$OMP        composite_grid_atom, composite_grid_coords, composite_nflat, cross_cutoff, &
!$OMP        drho_h, drho_s, grid_atom, harmonics, image_shell, lsd, nspins, &
!$OMP        particle_set, rho_h, rho_s, source_atom, tau_h, tau_s, &
!$OMP        use_atom_composite_density, use_atom_composite_gradient, use_atom_composite_tau)
                     DO composite_row = 1, composite_nflat
                        target_atom = composite_grid_atom(composite_row)
                        fractional = 0.0_dp
                        DO idir = 1, 3
                           DO jdir = 1, 3
                              fractional(idir) = fractional(idir) + cell%h_inv(idir, jdir)* &
                                                 (composite_grid_coords(jdir, composite_row) - &
                                                  particle_set(source_atom)%r(jdir))
                           END DO
                        END DO
                        CALL atom_grid_image_bounds(cell, fractional, cross_cutoff, image_shell, &
                                                    image_lower, image_upper)
                        DO image_i3 = image_lower(3), image_upper(3)
                           DO image_i2 = image_lower(2), image_upper(2)
                              DO image_i1 = image_lower(1), image_upper(1)
                                 image_shift = [image_i1, image_i2, image_i3]
                                 IF (target_atom == source_atom .AND. &
                                     ALL(image_shift == 0)) CYCLE
                                 image_translation = MATMUL( &
                                                     cell%hmat, REAL(image_shift, dp))
                                 cross_displacement = composite_grid_coords(:, composite_row) - &
                                                      particle_set(source_atom)%r - &
                                                      image_translation
                                 IF (SQRT(SUM(cross_displacement**2)) > cross_cutoff) CYCLE
                                 CALL interpolate_gapw_atom_grid_fields( &
                                    grid_atom, harmonics, cross_displacement, cross_cutoff, nspins, &
                                    rho_h, rho_s, drho_h, drho_s, tau_h, tau_s, &
                                    cross_density, cross_grad, cross_kin, cross_density_spatial, &
                                    cross_grad_spatial, cross_kin_spatial, calculate_spatial=.FALSE.)
                                 IF (lsd) THEN
                                    IF (use_atom_composite_density) THEN
                                       composite_cross_density(composite_row, 1:2) = &
                                          composite_cross_density(composite_row, 1:2) + &
                                          cross_density(1:2)
                                    END IF
                                    IF (use_atom_composite_gradient) THEN
                                       composite_cross_grad(composite_row, :, 1:2) = &
                                          composite_cross_grad(composite_row, :, 1:2) + &
                                          cross_grad(:, 1:2)
                                    END IF
                                    IF (use_atom_composite_tau) THEN
                                       composite_cross_kin(composite_row, 1:2) = &
                                          composite_cross_kin(composite_row, 1:2) + cross_kin(1:2)
                                    END IF
                                 ELSE
                                    IF (use_atom_composite_density) THEN
                                       composite_cross_density(composite_row, :) = &
                                          composite_cross_density(composite_row, :) + &
                                          0.5_dp*cross_density(1)
                                    END IF
                                    IF (use_atom_composite_gradient) THEN
                                       DO idir = 1, 3
                                          composite_cross_grad(composite_row, idir, :) = &
                                             composite_cross_grad(composite_row, idir, :) + &
                                             0.5_dp*cross_grad(idir, 1)
                                       END DO
                                    END IF
                                    IF (use_atom_composite_tau) THEN
                                       composite_cross_kin(composite_row, :) = &
                                          composite_cross_kin(composite_row, :) + &
                                          0.5_dp*cross_kin(1)
                                    END IF
                                 END IF
                              END DO
                           END DO
                        END DO
                     END DO
!$OMP END PARALLEL DO
                  END DO

                  CALL release_tau_basis_cache(tau_basis_cache)
               END DO
            END IF

            IF (native_grid_diagnostics) THEN
               composite_cross_density_max = MAXVAL(ABS(composite_cross_density))
               composite_cross_grad_max = MAXVAL(ABS(composite_cross_grad))
               composite_cross_kin_max = MAXVAL(ABS(composite_cross_kin))
               CALL para_env%max(composite_cross_cutoff_max)
               CALL para_env%max(composite_cross_density_max)
               CALL para_env%max(composite_cross_grad_max)
               CALL para_env%max(composite_cross_kin_max)
               iw = cp_logger_get_default_io_unit()
               IF (iw > 0) THEN
                  WRITE (UNIT=iw, FMT="(T2,A,4(1X,ES20.12))") &
                     "SKALA_GPW| Atom-composite cross support/maxima", &
                     composite_cross_cutoff_max, composite_cross_density_max, &
                     composite_cross_grad_max, composite_cross_kin_max
               END IF
            END IF
            composite_density(:, :) = composite_density(:, :) + composite_cross_density(:, :)
            composite_grad(:, :, :) = composite_grad(:, :, :) + composite_cross_grad(:, :, :)
            composite_kin(:, :) = composite_kin(:, :) + composite_cross_kin(:, :)
            DEALLOCATE (composite_cross_density, composite_cross_grad, composite_cross_kin)

            CPASSERT(ALL(composite_grid_weights >= 0.0_dp))
            IF (atom_composite_diagnostic .OR. native_grid_diagnostics) THEN
               atom_composite_nelec = SUM(composite_grid_weights* &
                                          (composite_density(:, 1) + composite_density(:, 2)))
               CALL para_env%sum(atom_composite_nelec)
            END IF
            IF (atom_composite_reference .AND. my_calculate_forces) THEN
               CALL skala_gapw_atom_composite_energy( &
                  my_xc_section, para_env, composite_density, composite_grad, composite_kin, &
                  composite_grid_coords, composite_grid_weights, composite_base_grid_weights, &
                  composite_local_grid_sizes, composite_local_atom_coords, atom_composite_exc, &
                  composite_density_grad, composite_grad_grad, composite_kin_grad, &
                  composite_grid_coord_grad, composite_grid_weight_grad, &
                  composite_atomic_grid_weight_grad, composite_atom_coord_grad, &
                  static_geometry_valid=comp_geometry_valid)
               ALLOCATE (composite_cross_force(3, SIZE(particle_set)), &
                         composite_explicit_force(3, SIZE(particle_set)), &
                         composite_grid_coord_force(3, SIZE(particle_set)), &
                         composite_model_atom_force(3, SIZE(particle_set)), &
                         composite_moving_smooth_force(3, SIZE(particle_set)), &
                         composite_nlcc_center_force(3, SIZE(particle_set)), &
                         composite_nlcc_target_force(3, SIZE(particle_set)), &
                         composite_partition_force(3, SIZE(particle_set)))
               composite_cross_force = 0.0_dp
               composite_cross_image_virial = 0.0_dp
               composite_model_atom_force = 0.0_dp
               composite_grid_coord_force = 0.0_dp
               composite_moving_smooth_force = 0.0_dp
               composite_nlcc_center_force = 0.0_dp
               composite_nlcc_target_force = 0.0_dp
               composite_partition_force = 0.0_dp
               composite_explicit_virial = 0.0_dp
               composite_feature_virial = 0.0_dp
               composite_interpolation_virial = 0.0_dp
               composite_partition_strain_virial = 0.0_dp
               DO composite_local_atom = 1, composite_local_natom
                  composite_model_atom_force(:, composite_local_atoms(composite_local_atom)) = &
                     composite_atom_coord_grad(:, composite_local_atom)
               END DO
               nchunks_per_atom = 4
               ALLOCATE (chunk_atom(composite_local_natom*nchunks_per_atom), &
                         chunk_begin(composite_local_natom*nchunks_per_atom), &
                         chunk_end(composite_local_natom*nchunks_per_atom))
               itask = 0
               DO composite_local_atom = 1, composite_local_natom
                  iatom = composite_local_atoms(composite_local_atom)
                  nrows_chunk = composite_atom_end(iatom) - composite_atom_start(iatom) + 1
                  DO ichunk = 1, nchunks_per_atom
                     chunk_start_row = composite_atom_start(iatom) + &
                                       (ichunk - 1)*nrows_chunk/nchunks_per_atom
                     chunk_end_row = composite_atom_start(iatom) + &
                                     ichunk*nrows_chunk/nchunks_per_atom - 1
                     IF (chunk_end_row < chunk_start_row) CYCLE
                     itask = itask + 1
                     chunk_atom(itask) = iatom
                     chunk_begin(itask) = chunk_start_row
                     chunk_end(itask) = chunk_end_row
                  END DO
               END DO
               ntasks = itask
               CALL timeset("calculate_vxc_atom_force", force_handle)
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP PRIVATE(composite_local_atom, composite_row, composite_smooth_gradient_value, &
!$OMP         composite_value_batch, iatom, idir, ispin, jdir, interpolation_stencil) &
!$OMP SHARED(cell, composite_atom_end, composite_atom_start, composite_grad_grad, &
!$OMP        composite_grid_coords, composite_local_atoms, composite_local_natom, &
!$OMP        lsd, qs_env, smooth_rho_r) &
!$OMP REDUCTION(+:composite_feature_virial)
               DO composite_local_atom = 1, composite_local_natom
                  iatom = composite_local_atoms(composite_local_atom)
                  DO composite_row = composite_atom_start(iatom), composite_atom_end(iatom)
                     CALL get_native_grid_stencil(smooth_rho_r(1)%pw_grid, cell, &
                                                  composite_row, composite_grid_coords(:, composite_row), &
                                                  interpolation_stencil, .FALSE.)
                     IF (lsd) THEN
                        CALL interpolate_native_grid_values_compact(lsd_value_slots, 6, &
                                                                    interpolation_stencil, composite_value_batch)
                        DO jdir = 1, 3
                           composite_smooth_gradient_value(jdir, 1) = composite_value_batch(jdir)
                           composite_smooth_gradient_value(jdir, 2) = composite_value_batch(3 + jdir)
                        END DO
                     ELSE
                        CALL interpolate_native_grid_values_compact(rks_value_slots, 3, &
                                                                    interpolation_stencil, composite_value_batch)
                        DO jdir = 1, 3
                           composite_smooth_gradient_value(jdir, :) = 0.5_dp*composite_value_batch(jdir)
                        END DO
                     END IF
                     DO ispin = 1, 2
                        DO idir = 1, 3
                           DO jdir = 1, 3
                              composite_feature_virial(jdir, idir) = &
                                 composite_feature_virial(jdir, idir) - &
                                 composite_grad_grad(composite_row, idir, ispin)* &
                                 composite_smooth_gradient_value(jdir, ispin)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
!$OMP END PARALLEL DO
               ! The periodic image layout only depends on the geometry, not on
               ! the target atom, so build it once and share it across threads.
               IF (.NOT. composite_partition_layout_valid) THEN
                  CALL build_periodic_atom_image_layout(composite_atom_coords, cell, 1, &
                                                        composite_image_periodicity, &
                                                        composite_partition_layout, local_partition_target_image)
                  composite_partition_layout_valid = .TRUE.
               END IF
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(composite_row, itask, &
!$OMP         gth_potential, iatom, idir, ispin, jdir, &
!$OMP         nlcc_density, nlcc_gradient, nlcc_hessian, nlcc_spatial_derivative, &
!$OMP         nlcc_spin_factor, composite_nlcc_center_force_local, &
!$OMP         composite_partition_force_local, composite_grid_coord_force_local, &
!$OMP         composite_moving_smooth_force_local, composite_nlcc_target_force_local, &
!$OMP         local_partition_datom, local_partition_dstrain, &
!$OMP         local_partition_target_image, local_partition_weight, sgp_potential, source_atom, &
!$OMP         spatial_derivative, composite_grad_batch, target_atom, target_partition_adjoint, &
!$OMP         interpolation_stencil) &
!$OMP REDUCTION(+:composite_interpolation_virial, composite_partition_strain_virial) &
!$OMP SHARED(cell, chunk_atom, chunk_begin, chunk_end, composite_atom_coords, composite_atom_end, &
!$OMP        composite_atom_kind, composite_atom_start, composite_atomic_grid_weight_grad, &
!$OMP        composite_base_grid_weights, composite_density_grad, &
!$OMP        composite_grad_grad, composite_grid_coord_force, composite_grid_coord_grad, &
!$OMP        composite_grid_coords, composite_grid_weight_grad, composite_image_periodicity, &
!$OMP        composite_kin_grad, composite_local_atoms, composite_local_natom, &
!$OMP        composite_moving_smooth_force, &
!$OMP        composite_nlcc_center_force, composite_nlcc_target_force, composite_partition_force, &
!$OMP        composite_partition_layout, lsd, my_kind_set, ntasks, &
!$OMP        nlcc, particle_set, qs_env, smooth_rho_r)
               ALLOCATE (composite_nlcc_center_force_local(3, SIZE(particle_set)), &
                         composite_partition_force_local(3, SIZE(particle_set)), &
                         composite_grid_coord_force_local(3, SIZE(particle_set)), &
                         composite_moving_smooth_force_local(3, SIZE(particle_set)), &
                         composite_nlcc_target_force_local(3, SIZE(particle_set)), &
                         local_partition_datom(3, SIZE(particle_set)))
               composite_nlcc_center_force_local = 0.0_dp
               composite_partition_force_local = 0.0_dp
               composite_grid_coord_force_local = 0.0_dp
               composite_moving_smooth_force_local = 0.0_dp
               composite_nlcc_target_force_local = 0.0_dp
!$OMP DO SCHEDULE(DYNAMIC)
               DO itask = 1, ntasks
                  iatom = chunk_atom(itask)
                  local_partition_target_image = iatom
                  DO composite_row = chunk_begin(itask), chunk_end(itask)
                     composite_grid_coord_force_local(:, iatom) = &
                        composite_grid_coord_force_local(:, iatom) + &
                        composite_grid_coord_grad(:, composite_row)
                     CALL get_native_grid_stencil(smooth_rho_r(1)%pw_grid, cell, &
                                                  composite_row, composite_grid_coords(:, composite_row), &
                                                  interpolation_stencil, .FALSE.)
                     IF (lsd) THEN
                        CALL interpolate_native_grid_gradients_compact(lsd_grad_slots, 10, &
                                                                       smooth_rho_r(1)%pw_grid, &
                                                                       interpolation_stencil, composite_grad_batch)
                        spatial_derivative = &
                           composite_density_grad(composite_row, 1)*composite_grad_batch(:, 1) + &
                           composite_density_grad(composite_row, 2)*composite_grad_batch(:, 2) + &
                           composite_kin_grad(composite_row, 1)*composite_grad_batch(:, 3) + &
                           composite_kin_grad(composite_row, 2)*composite_grad_batch(:, 4)
                        DO idir = 1, 3
                           spatial_derivative = spatial_derivative + &
                                                composite_grad_grad(composite_row, idir, 1)* &
                                                composite_grad_batch(:, 4 + idir) + &
                                                composite_grad_grad(composite_row, idir, 2)* &
                                                composite_grad_batch(:, 7 + idir)
                        END DO
                     ELSE
                        CALL interpolate_native_grid_gradients_compact(rks_grad_slots, 5, &
                                                                       smooth_rho_r(1)%pw_grid, &
                                                                       interpolation_stencil, composite_grad_batch)
                        spatial_derivative = 0.5_dp*SUM(composite_density_grad(composite_row, :))* &
                                             composite_grad_batch(:, 1) + &
                                             0.5_dp*SUM(composite_kin_grad(composite_row, :))* &
                                             composite_grad_batch(:, 2)
                        DO idir = 1, 3
                           spatial_derivative = spatial_derivative + &
                                                0.5_dp*SUM(composite_grad_grad(composite_row, idir, :))* &
                                                composite_grad_batch(:, 2 + idir)
                        END DO
                     END IF
                     DO idir = 1, 3
                        DO jdir = 1, 3
                           composite_interpolation_virial(idir, jdir) = &
                              composite_interpolation_virial(idir, jdir) + &
                              spatial_derivative(idir)*( &
                              composite_grid_coords(jdir, composite_row) - &
                              particle_set(iatom)%r(jdir))
                        END DO
                     END DO
                     composite_moving_smooth_force_local(:, iatom) = &
                        composite_moving_smooth_force_local(:, iatom) + spatial_derivative
                     IF (nlcc) THEN
                        nlcc_spin_factor = 0.5_dp
                        DO source_atom = 1, SIZE(particle_set)
                           NULLIFY (gth_potential, sgp_potential)
                           CALL get_qs_kind(my_kind_set(composite_atom_kind(source_atom)), &
                                            gth_potential=gth_potential, &
                                            sgp_potential=sgp_potential)
                           CALL evaluate_nlcc_primitive_fields( &
                              composite_grid_coords(:, composite_row), &
                              particle_set(source_atom)%r, gth_potential, sgp_potential, &
                              nlcc_density, nlcc_gradient, nlcc_hessian)
                           nlcc_spatial_derivative = 0.0_dp
                           DO ispin = 1, 2
                              nlcc_spatial_derivative = nlcc_spatial_derivative + &
                                                        nlcc_spin_factor*composite_density_grad(composite_row, ispin)* &
                                                        nlcc_gradient
                              DO idir = 1, 3
                                 DO jdir = 1, 3
                                    nlcc_spatial_derivative(jdir) = &
                                       nlcc_spatial_derivative(jdir) + nlcc_spin_factor* &
                                       composite_grad_grad(composite_row, idir, ispin)* &
                                       nlcc_hessian(idir, jdir)
                                 END DO
                              END DO
                           END DO
                           composite_nlcc_target_force_local(:, iatom) = &
                              composite_nlcc_target_force_local(:, iatom) + nlcc_spatial_derivative
                           composite_nlcc_center_force_local(:, source_atom) = &
                              composite_nlcc_center_force_local(:, source_atom) - nlcc_spatial_derivative
                        END DO
                     END IF
                     CALL periodic_atom_image_partition_from_layout( &
                        composite_grid_coords(:, composite_row), composite_partition_layout, &
                        local_partition_target_image, local_partition_weight, &
                        local_partition_datom, local_partition_dstrain)
                     target_partition_adjoint = composite_base_grid_weights(composite_row)* &
                                                composite_grid_weight_grad(composite_row)
                     DO target_atom = 1, SIZE(particle_set)
                        composite_partition_force_local(:, target_atom) = &
                           composite_partition_force_local(:, target_atom) + &
                           target_partition_adjoint*local_partition_datom(:, target_atom)
                     END DO
                     composite_partition_force_local(:, iatom) = &
                        composite_partition_force_local(:, iatom) - target_partition_adjoint* &
                        SUM(local_partition_datom, DIM=2)
                     composite_partition_strain_virial = &
                        composite_partition_strain_virial - target_partition_adjoint* &
                        local_partition_dstrain
                  END DO
               END DO
!$OMP END DO
!$OMP CRITICAL(skala_atom_composite_force_reduction)
               composite_nlcc_center_force(:, :) = composite_nlcc_center_force(:, :) + &
                                                   composite_nlcc_center_force_local
               composite_partition_force(:, :) = composite_partition_force(:, :) + &
                                                 composite_partition_force_local
               composite_grid_coord_force(:, :) = composite_grid_coord_force(:, :) + &
                                                  composite_grid_coord_force_local
               composite_moving_smooth_force(:, :) = composite_moving_smooth_force(:, :) + &
                                                     composite_moving_smooth_force_local
               composite_nlcc_target_force(:, :) = composite_nlcc_target_force(:, :) + &
                                                   composite_nlcc_target_force_local
!$OMP END CRITICAL(skala_atom_composite_force_reduction)
               DEALLOCATE (composite_nlcc_center_force_local, composite_partition_force_local, &
                           composite_grid_coord_force_local, composite_moving_smooth_force_local, &
                           composite_nlcc_target_force_local, local_partition_datom)
!$OMP END PARALLEL
               DEALLOCATE (chunk_atom, chunk_begin, chunk_end)
               CALL timestop(force_handle)
               ! The target atom grids can overlap augmentation regions of other atoms
               ! and periodic images. Differentiate the same discrete interpolation used
               ! in the forward composite fields, including the explicit image strain.
               IF (.NOT. direct_valence_atom_composite) THEN
                  DO ikind = 1, SIZE(atomic_kind_set)
                     CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
                     NULLIFY (gth_potential, sgp_potential)
                     CALL get_qs_kind(my_kind_set(ikind), paw_atom=paw_atom, &
                                      gth_potential=gth_potential, harmonics=harmonics, &
                                      grid_atom=grid_atom, sgp_potential=sgp_potential, &
                                      zatom=zatom, zeff=zeff)
                     CALL get_qs_kind(my_kind_set(ikind), basis_set=basis_1c, basis_type="GAPW_1C")
                     IF (.NOT. native_skala_uses_one_center_kind( &
                         paw_atom, gapw_representation, &
                         ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential), &
                         zeff, zatom)) CYCLE

                     nr = grid_atom%nr
                     na = grid_atom%ng_sphere
                     CALL create_tau_basis_cache(tau_basis_cache, grid_atom, basis_1c, harmonics)
                     CALL reallocate(rho_h, 1, na, 1, nr, 1, nspins)
                     CALL reallocate(rho_s, 1, na, 1, nr, 1, nspins)
                     CALL reallocate(drho_h, 1, 4, 1, na, 1, nr, 1, nspins)
                     CALL reallocate(drho_s, 1, 4, 1, na, 1, nr, 1, nspins)
                     CALL reallocate(tau_h, 1, na, 1, nr, 1, nspins)
                     CALL reallocate(tau_s, 1, na, 1, nr, 1, nspins)

                     DO iat = 1, natom
                        source_atom = atom_list(iat)
                        rho_atom => my_rho_atom_set(source_atom)
                        NULLIFY (r_h, r_s, dr_h, dr_s, r_h_d, r_s_d)
                        CALL get_rho_atom(rho_atom=rho_atom, rho_rad_h=r_h, rho_rad_s=r_s, &
                                          drho_rad_h=dr_h, drho_rad_s=dr_s, &
                                          rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d)
                        rho_h = 0.0_dp
                        rho_s = 0.0_dp
                        drho_h = 0.0_dp
                        drho_s = 0.0_dp
                        CALL calc_tau_atom(tau_h, tau_s, rho_atom, tau_basis_cache, nspins)
                        DO ir = 1, nr
                           CALL calc_rho_angular(grid_atom, harmonics, nspins, .TRUE., &
                                                 ir, r_h, r_s, rho_h, rho_s, dr_h, dr_s, &
                                                 r_h_d, r_s_d, drho_h, drho_s)
                        END DO

                        cross_cutoff = gapw_atom_grid_support_radius( &
                                       grid_atom, rho_h, rho_s, drho_h, drho_s, tau_h, tau_s)
                        IF (cross_cutoff <= 0.0_dp) CYCLE
                        image_shell = 0
                        DO idir = 1, 3
                           IF (cell%perd(idir) == 1) THEN
                              image_shell(idir) = CEILING( &
                                                  cross_cutoff*SQRT(SUM(cell%h_inv(idir, :)**2))) + 1
                           END IF
                        END DO

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(image_lower, image_upper, composite_cross_force_local, composite_cross_image_virial_local, &
!$OMP         composite_row, cross_density, cross_density_adjoint, cross_density_spatial, &
!$OMP         cross_displacement, cross_grad, cross_grad_adjoint, cross_grad_spatial, &
!$OMP         cross_kin, cross_kin_adjoint, cross_kin_spatial, cross_spatial_derivative, &
!$OMP         fractional, idir, image_i1, image_i2, image_i3, image_shift, image_translation, &
!$OMP         ispin, jdir, target_atom) &
!$OMP SHARED(cell, composite_cross_force, composite_cross_image_virial, &
!$OMP        composite_density_grad, composite_grad_grad, composite_grid_atom, &
!$OMP        composite_grid_coords, composite_kin_grad, composite_nflat, cross_cutoff, &
!$OMP        drho_h, drho_s, grid_atom, harmonics, image_shell, lsd, nspins, &
!$OMP        particle_set, rho_h, rho_s, source_atom, tau_h, tau_s, &
!$OMP        use_atom_composite_density, use_atom_composite_gradient, use_atom_composite_tau)
                        ALLOCATE (composite_cross_force_local(3, SIZE(particle_set)))
                        composite_cross_force_local = 0.0_dp
                        composite_cross_image_virial_local = 0.0_dp
!$OMP DO SCHEDULE(STATIC)
                        DO composite_row = 1, composite_nflat
                           target_atom = composite_grid_atom(composite_row)
                           IF (lsd) THEN
                              cross_density_adjoint = 0.0_dp
                              cross_grad_adjoint = 0.0_dp
                              cross_kin_adjoint = 0.0_dp
                              IF (use_atom_composite_density) THEN
                                 cross_density_adjoint(1:2) = &
                                    composite_density_grad(composite_row, 1:2)
                              END IF
                              IF (use_atom_composite_gradient) THEN
                                 cross_grad_adjoint(:, 1:2) = &
                                    composite_grad_grad(composite_row, :, 1:2)
                              END IF
                              IF (use_atom_composite_tau) THEN
                                 cross_kin_adjoint(1:2) = &
                                    composite_kin_grad(composite_row, 1:2)
                              END IF
                           ELSE
                              cross_density_adjoint = 0.0_dp
                              cross_grad_adjoint = 0.0_dp
                              cross_kin_adjoint = 0.0_dp
                              IF (use_atom_composite_density) THEN
                                 cross_density_adjoint(1) = 0.5_dp* &
                                                            SUM(composite_density_grad(composite_row, :))
                              END IF
                              IF (use_atom_composite_gradient) THEN
                                 DO idir = 1, 3
                                    cross_grad_adjoint(idir, 1) = 0.5_dp* &
                                                                  SUM(composite_grad_grad(composite_row, idir, :))
                                 END DO
                              END IF
                              IF (use_atom_composite_tau) THEN
                                 cross_kin_adjoint(1) = 0.5_dp* &
                                                        SUM(composite_kin_grad(composite_row, :))
                              END IF
                           END IF

                           fractional = 0.0_dp
                           DO idir = 1, 3
                              DO jdir = 1, 3
                                 fractional(idir) = fractional(idir) + cell%h_inv(idir, jdir)* &
                                                    (composite_grid_coords(jdir, composite_row) - &
                                                     particle_set(source_atom)%r(jdir))
                              END DO
                           END DO
                           CALL atom_grid_image_bounds(cell, fractional, cross_cutoff, image_shell, &
                                                       image_lower, image_upper)
                           DO image_i3 = image_lower(3), image_upper(3)
                              DO image_i2 = image_lower(2), image_upper(2)
                                 DO image_i1 = image_lower(1), image_upper(1)
                                    image_shift = [image_i1, image_i2, image_i3]
                                    IF (target_atom == source_atom .AND. &
                                        ALL(image_shift == 0)) CYCLE
                                    image_translation = MATMUL( &
                                                        cell%hmat, REAL(image_shift, dp))
                                    cross_displacement = &
                                       composite_grid_coords(:, composite_row) - &
                                       particle_set(source_atom)%r - image_translation
                                    IF (SQRT(SUM(cross_displacement**2)) > cross_cutoff) CYCLE
                                    CALL interpolate_gapw_atom_grid_fields( &
                                       grid_atom, harmonics, cross_displacement, cross_cutoff, &
                                       nspins, rho_h, rho_s, drho_h, drho_s, tau_h, tau_s, &
                                       cross_density, cross_grad, cross_kin, cross_density_spatial, &
                                       cross_grad_spatial, cross_kin_spatial)
                                    cross_spatial_derivative = 0.0_dp
                                    DO ispin = 1, nspins
                                       DO idir = 1, 3
                                          cross_spatial_derivative(idir) = &
                                             cross_spatial_derivative(idir) + &
                                             cross_density_adjoint(ispin)* &
                                             cross_density_spatial(idir, ispin) + &
                                             cross_kin_adjoint(ispin)* &
                                             cross_kin_spatial(idir, ispin)
                                          DO jdir = 1, 3
                                             cross_spatial_derivative(idir) = &
                                                cross_spatial_derivative(idir) + &
                                                cross_grad_adjoint(jdir, ispin)* &
                                                cross_grad_spatial(jdir, idir, ispin)
                                          END DO
                                       END DO
                                    END DO
                                    composite_cross_force_local(:, target_atom) = &
                                       composite_cross_force_local(:, target_atom) + &
                                       cross_spatial_derivative
                                    composite_cross_force_local(:, source_atom) = &
                                       composite_cross_force_local(:, source_atom) - &
                                       cross_spatial_derivative
                                    DO idir = 1, 3
                                       DO jdir = 1, 3
                                          composite_cross_image_virial_local(idir, jdir) = &
                                             composite_cross_image_virial_local(idir, jdir) + &
                                             cross_spatial_derivative(idir)*image_translation(jdir)
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
!$OMP END DO
!$OMP CRITICAL(skala_atom_composite_cross_reduction)
                        composite_cross_force(:, :) = &
                           composite_cross_force(:, :) + composite_cross_force_local(:, :)
                        composite_cross_image_virial = composite_cross_image_virial + &
                                                       composite_cross_image_virial_local
!$OMP END CRITICAL(skala_atom_composite_cross_reduction)
                        DEALLOCATE (composite_cross_force_local)
!$OMP END PARALLEL
                     END DO
                     CALL release_tau_basis_cache(tau_basis_cache)
                  END DO
               END IF

               composite_explicit_force(:, :) = composite_model_atom_force(:, :) + &
                                                composite_grid_coord_force(:, :) + &
                                                composite_moving_smooth_force(:, :) + &
                                                composite_cross_force(:, :) + &
                                                composite_nlcc_center_force(:, :) + &
                                                composite_nlcc_target_force(:, :) + &
                                                composite_partition_force(:, :)
               ! CP2K stores +dE/dR in the electronic force components, while the
               ! virial is -dE/dstrain. Model coordinates, atom-grid centers,
               ! NLCC centers, and partition centers move affinely with their atoms.
               ! The smooth-field interpolation force is excluded here: its affine
               ! response is already in the PW stress and its non-affine local-grid
               ! correction is composite_interpolation_virial.
               DO iatom = 1, SIZE(particle_set)
                  DO idir = 1, 3
                     DO jdir = 1, 3
                        composite_explicit_virial(idir, jdir) = &
                           composite_explicit_virial(idir, jdir) - ( &
                           composite_model_atom_force(idir, iatom) + &
                           composite_grid_coord_force(idir, iatom) + &
                           composite_cross_force(idir, iatom) + &
                           composite_nlcc_center_force(idir, iatom) + &
                           composite_nlcc_target_force(idir, iatom) + &
                           composite_partition_force(idir, iatom))*particle_set(iatom)%r(jdir)
                     END DO
                  END DO
               END DO
               composite_explicit_virial = composite_explicit_virial + &
                                           composite_partition_strain_virial + &
                                           composite_cross_image_virial
               IF (ASSOCIATED(force)) THEN
                  DO iatom = 1, SIZE(particle_set)
                     ikind = composite_atom_kind(iatom)
                     iat = composite_atom_kind_index(iatom)
                     CPASSERT(ikind > 0 .AND. iat > 0)
                     force(ikind)%rho_elec(:, iat) = force(ikind)%rho_elec(:, iat) + &
                                                     composite_explicit_force(:, iatom)
                  END DO
               END IF
               IF (use_virial) THEN
                  ! The local radial vectors of the atom-centered quadrature do not
                  ! deform with the periodic cell. Correct the standard PW response,
                  ! which follows fixed fractional coordinates, by the corresponding
                  ! non-affine interpolation derivative.
                  virial%pv_xc = virial%pv_xc + composite_feature_virial - &
                                 composite_interpolation_virial
                  virial%pv_gapw = virial%pv_gapw + composite_explicit_virial
                  virial%pv_virial = virial%pv_virial + composite_explicit_virial
               END IF
               IF (native_grid_diagnostics) THEN
                  CALL para_env%sum(composite_cross_force)
                  CALL para_env%sum(composite_cross_image_virial)
                  CALL para_env%sum(composite_model_atom_force)
                  CALL para_env%sum(composite_grid_coord_force)
                  CALL para_env%sum(composite_moving_smooth_force)
                  CALL para_env%sum(composite_nlcc_center_force)
                  CALL para_env%sum(composite_nlcc_target_force)
                  CALL para_env%sum(composite_partition_force)
                  CALL para_env%sum(composite_explicit_force)
                  CALL para_env%sum(composite_explicit_virial)
                  CALL para_env%sum(composite_feature_virial)
                  CALL para_env%sum(composite_interpolation_virial)
                  iw = cp_logger_get_default_io_unit()
                  IF (iw > 0) THEN
                     DO iatom = 1, SIZE(particle_set)
                        WRITE (UNIT=iw, FMT="(T2,A,1X,I0,3(1X,ES20.12))") &
                           "SKALA_GPW| Atom-composite model-atom force", iatom, &
                           composite_model_atom_force(:, iatom)
                        WRITE (UNIT=iw, FMT="(T2,A,1X,I0,3(1X,ES20.12))") &
                           "SKALA_GPW| Atom-composite grid-coordinate force", iatom, &
                           composite_grid_coord_force(:, iatom)
                        WRITE (UNIT=iw, FMT="(T2,A,1X,I0,3(1X,ES20.12))") &
                           "SKALA_GPW| Atom-composite moving-smooth force", iatom, &
                           composite_moving_smooth_force(:, iatom)
                        WRITE (UNIT=iw, FMT="(T2,A,1X,I0,3(1X,ES20.12))") &
                           "SKALA_GPW| Atom-composite cross-region force", iatom, &
                           composite_cross_force(:, iatom)
                        WRITE (UNIT=iw, FMT="(T2,A,1X,I0,3(1X,ES20.12))") &
                           "SKALA_GPW| Atom-composite NLCC-center force", iatom, &
                           composite_nlcc_center_force(:, iatom)
                        WRITE (UNIT=iw, FMT="(T2,A,1X,I0,3(1X,ES20.12))") &
                           "SKALA_GPW| Atom-composite NLCC-target force", iatom, &
                           composite_nlcc_target_force(:, iatom)
                        WRITE (UNIT=iw, FMT="(T2,A,1X,I0,3(1X,ES20.12))") &
                           "SKALA_GPW| Atom-composite partition force", iatom, &
                           composite_partition_force(:, iatom)
                        WRITE (UNIT=iw, FMT="(T2,A,1X,I0,3(1X,ES20.12))") &
                           "SKALA_GPW| Atom-composite explicit force", iatom, &
                           composite_explicit_force(:, iatom)
                     END DO
                     WRITE (UNIT=iw, FMT="(T2,A)") &
                        "SKALA_GPW| Atom-composite explicit virial"
                     DO idir = 1, 3
                        WRITE (UNIT=iw, FMT="(T2,A,1X,3ES20.10)") &
                           "SKALA_GPW|", composite_explicit_virial(idir, :)
                     END DO
                     WRITE (UNIT=iw, FMT="(T2,A)") &
                        "SKALA_GPW| Atom-composite cross-image virial"
                     DO idir = 1, 3
                        WRITE (UNIT=iw, FMT="(T2,A,1X,3ES20.10)") &
                           "SKALA_GPW|", composite_cross_image_virial(idir, :)
                     END DO
                     WRITE (UNIT=iw, FMT="(T2,A)") &
                        "SKALA_GPW| Atom-composite feature virial"
                     DO idir = 1, 3
                        WRITE (UNIT=iw, FMT="(T2,A,1X,3ES20.10)") &
                           "SKALA_GPW|", composite_feature_virial(idir, :)
                     END DO
                     WRITE (UNIT=iw, FMT="(T2,A)") &
                        "SKALA_GPW| Atom-composite interpolation virial"
                     DO idir = 1, 3
                        WRITE (UNIT=iw, FMT="(T2,A,1X,3ES20.10)") &
                           "SKALA_GPW|", composite_interpolation_virial(idir, :)
                     END DO
                  END IF
               END IF
               DEALLOCATE (composite_atom_coord_grad, composite_atomic_grid_weight_grad, &
                           composite_cross_force, &
                           composite_explicit_force, composite_grid_coord_force, &
                           composite_grid_coord_grad, composite_grid_weight_grad, &
                           composite_model_atom_force, composite_moving_smooth_force, &
                           composite_nlcc_center_force, &
                           composite_nlcc_target_force, &
                           composite_partition_force)
            ELSE
               CALL skala_gapw_atom_composite_energy( &
                  my_xc_section, para_env, composite_density, composite_grad, composite_kin, &
                  composite_grid_coords, composite_grid_weights, composite_base_grid_weights, &
                  composite_local_grid_sizes, composite_local_atom_coords, atom_composite_exc, &
                  composite_density_grad, composite_grad_grad, composite_kin_grad, &
                  static_geometry_valid=comp_geometry_valid)
            END IF
            adjoint_nchannels = MERGE(2, 1, lsd)
            ALLOCATE (adjoint_compact(skala_state%native_atom%native_need_nfg, MAX(skala_state%native_atom%native_need_nneed, 1)))
            adjoint_compact = 0.0_dp
            nthread = 1
!$          nthread = omp_get_max_threads()
            CALL ensure_adjoint_thread_scratch(nthread, SIZE(adjoint_compact, 1), &
                                               SIZE(adjoint_compact, 2))
            CALL timeset("calculate_vxc_atom_adjoint", adjoint_handle)
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(composite_row, composite_smooth_density_adjoint_value, inode_x, inode_y, inode_z, &
!$OMP         composite_smooth_gradient_adjoint_value, composite_smooth_kin_adjoint_value, &
!$OMP         interpolation_stencil, thread_adj, thread_id, lin, pos, coefficient, idir) &
!$OMP SHARED(adjoint_compact, cell, composite_density_grad, skala_state, &
!$OMP        composite_grad_grad, composite_grid_coords, composite_kin_grad, composite_nflat, &
!$OMP        lsd, qs_env, &
!$OMP        smooth_rho_r)
            thread_id = 1
!$          thread_id = omp_get_thread_num() + 1
            thread_adj => skala_state%native_atom%adjoint_thread_scratch(1:SIZE(adjoint_compact, 1), &
                                                                         1:SIZE(adjoint_compact, 2), thread_id)
            thread_adj = 0.0_dp
!$OMP DO SCHEDULE(DYNAMIC)
            DO composite_row = 1, composite_nflat
               CALL get_native_grid_stencil(smooth_rho_r(1)%pw_grid, cell, &
                                            composite_row, composite_grid_coords(:, composite_row), &
                                            interpolation_stencil, .FALSE.)
               IF (.NOT. interpolation_stencil%active) CYCLE
               IF (lsd) THEN
                  composite_smooth_density_adjoint_value = composite_density_grad(composite_row, :)
                  composite_smooth_gradient_adjoint_value = composite_grad_grad(composite_row, :, :)
                  composite_smooth_kin_adjoint_value = composite_kin_grad(composite_row, :)
               ELSE
                  composite_smooth_density_adjoint_value = 0.0_dp
                  composite_smooth_gradient_adjoint_value = 0.0_dp
                  composite_smooth_kin_adjoint_value = 0.0_dp
                  composite_smooth_density_adjoint_value(1) = &
                     SUM(composite_density_grad(composite_row, :))
                  composite_smooth_gradient_adjoint_value(:, 1) = &
                     SUM(composite_grad_grad(composite_row, :, :), DIM=2)
                  composite_smooth_kin_adjoint_value(1) = SUM(composite_kin_grad(composite_row, :))
               END IF
               DO inode_z = 1, native_grid_interp_npts
                  IF (.NOT. interpolation_stencil%valid(inode_z, 3)) CYCLE
                  DO inode_y = 1, native_grid_interp_npts
                     IF (.NOT. interpolation_stencil%valid(inode_y, 2)) CYCLE
                     DO inode_x = 1, native_grid_interp_npts
                        IF (.NOT. interpolation_stencil%valid(inode_x, 1)) CYCLE
                        lin = native_need_linear_index(interpolation_stencil, inode_x, inode_y, &
                                                       inode_z, skala_state%native_atom%native_need_npts)
                        pos = native_need_lookup(lin)
                        CPASSERT(pos > 0)
                        coefficient = interpolation_stencil%weight(inode_x, 1)* &
                                      interpolation_stencil%weight(inode_y, 2)* &
                                      interpolation_stencil%weight(inode_z, 3)
                        IF (lsd) THEN
                           thread_adj(1, pos) = thread_adj(1, pos) + &
                                                coefficient*composite_smooth_density_adjoint_value(1)
                           thread_adj(2, pos) = thread_adj(2, pos) + &
                                                coefficient*composite_smooth_density_adjoint_value(2)
                           thread_adj(3, pos) = thread_adj(3, pos) + &
                                                coefficient*composite_smooth_kin_adjoint_value(1)
                           thread_adj(4, pos) = thread_adj(4, pos) + &
                                                coefficient*composite_smooth_kin_adjoint_value(2)
                           DO idir = 1, 3
                              thread_adj(4 + idir, pos) = thread_adj(4 + idir, pos) + &
                                                          coefficient*composite_smooth_gradient_adjoint_value(idir, 1)
                              thread_adj(7 + idir, pos) = thread_adj(7 + idir, pos) + &
                                                          coefficient*composite_smooth_gradient_adjoint_value(idir, 2)
                           END DO
                        ELSE
                           thread_adj(1, pos) = thread_adj(1, pos) + &
                                                coefficient*composite_smooth_density_adjoint_value(1)
                           thread_adj(2, pos) = thread_adj(2, pos) + &
                                                coefficient*composite_smooth_kin_adjoint_value(1)
                           DO idir = 1, 3
                              thread_adj(2 + idir, pos) = thread_adj(2 + idir, pos) + &
                                                          coefficient*composite_smooth_gradient_adjoint_value(idir, 1)
                           END DO
                        END IF
                     END DO
                  END DO
               END DO
            END DO
!$OMP END DO
!$OMP CRITICAL(skala_atom_composite_adjoint_compact)
            adjoint_compact(:, :) = adjoint_compact(:, :) + thread_adj
!$OMP END CRITICAL(skala_atom_composite_adjoint_compact)
!$OMP END PARALLEL
            CALL timestop(adjoint_handle)
            ALLOCATE (adjoint_local_flat(skala_state%native_atom%native_need_nfg, smooth_rho_r(1)%pw_grid%ngpts_local))
            CALL native_grid_exchange_scatter(smooth_rho_r(1)%pw_grid%para%group, adjoint_compact, adjoint_local_flat)
            ngl = smooth_rho_r(1)%pw_grid%ngpts_local
            ALLOCATE (smooth_density_adjoint_storage(ngl, adjoint_nchannels), &
                      smooth_grad_adjoint_storage(ngl, 3, adjoint_nchannels), &
                      smooth_kin_adjoint_storage(ngl, adjoint_nchannels))
            DO ispin = 1, adjoint_nchannels
               IF (lsd) THEN
                  smooth_density_adjoint_storage(:, ispin) = adjoint_local_flat(ispin, :)
                  smooth_kin_adjoint_storage(:, ispin) = adjoint_local_flat(2 + ispin, :)
                  DO idir = 1, 3
                     smooth_grad_adjoint_storage(:, idir, ispin) = &
                        adjoint_local_flat(4 + 3*(ispin - 1) + idir, :)
                  END DO
               ELSE
                  smooth_density_adjoint_storage(:, 1) = adjoint_local_flat(1, :)
                  smooth_kin_adjoint_storage(:, 1) = adjoint_local_flat(2, :)
                  DO idir = 1, 3
                     smooth_grad_adjoint_storage(:, idir, 1) = adjoint_local_flat(2 + idir, :)
                  END DO
               END IF
            END DO
            DEALLOCATE (adjoint_compact, adjoint_local_flat)
            IF (atom_composite_diagnostic) THEN
               smooth_input_contraction = 0.0_dp
!$OMP PARALLEL DO SCHEDULE(STATIC) REDUCTION(+:smooth_input_contraction) DEFAULT(NONE) &
!$OMP PRIVATE(composite_smooth_density_value, composite_smooth_gradient_value, &
!$OMP         composite_smooth_kin_value, idir, ispin) &
!$OMP SHARED(composite_density_grad, composite_grad_grad, composite_kin_grad, composite_nflat, &
!$OMP        composite_smooth_density_cache, composite_smooth_gradient_cache, &
!$OMP        composite_smooth_kin_cache, lsd)
               DO composite_row = 1, composite_nflat
                  composite_smooth_density_value = composite_smooth_density_cache(composite_row, :)
                  composite_smooth_gradient_value = &
                     composite_smooth_gradient_cache(composite_row, :, :)
                  composite_smooth_kin_value = composite_smooth_kin_cache(composite_row, :)
                  DO ispin = 1, 2
                     IF (lsd) THEN
                        smooth_input_contraction = smooth_input_contraction + &
                                                   composite_density_grad(composite_row, ispin)* &
                                                   composite_smooth_density_value(ispin) + &
                                                   composite_kin_grad(composite_row, ispin)* &
                                                   composite_smooth_kin_value(ispin)
                        DO idir = 1, 3
                           smooth_input_contraction = smooth_input_contraction + &
                                                      composite_grad_grad(composite_row, idir, ispin)* &
                                                      composite_smooth_gradient_value(idir, ispin)
                        END DO
                     ELSE
                        smooth_input_contraction = smooth_input_contraction + 0.5_dp*( &
                                                   composite_density_grad(composite_row, ispin)* &
                                                   composite_smooth_density_value(1) + &
                                                   composite_kin_grad(composite_row, ispin)* &
                                                   composite_smooth_kin_value(1))
                        DO idir = 1, 3
                           smooth_input_contraction = smooth_input_contraction + 0.5_dp* &
                                                      composite_grad_grad(composite_row, idir, ispin)* &
                                                      composite_smooth_gradient_value(idir, 1)
                        END DO
                     END IF
                  END DO
               END DO
!$OMP END PARALLEL DO
            END IF
            IF (atom_composite_diagnostic) CALL para_env%sum(smooth_input_contraction)
            CALL build_vxc_from_feature_grads( &
               smooth_vxc_rho, smooth_vxc_tau, smooth_rho_r, auxbas_pw_pool, &
               smooth_density_adjoint_storage, smooth_grad_adjoint_storage, &
               smooth_kin_adjoint_storage, xc_deriv_method_id)
            DEALLOCATE (smooth_density_adjoint_storage, smooth_grad_adjoint_storage, &
                        smooth_kin_adjoint_storage)
            IF (atom_composite_diagnostic) THEN
               smooth_grid_contraction = 0.0_dp
            END IF
            DO ispin = 1, nspins
               IF (atom_composite_diagnostic) THEN
                  smooth_grid_contraction = smooth_grid_contraction + smooth_rho_r(1)%pw_grid%dvol* &
                                            (SUM(smooth_vxc_rho(ispin)%array* &
                                                 smooth_rho_r(ispin)%array) + &
                                             SUM(smooth_vxc_tau(ispin)%array* &
                                                 smooth_tau_r(ispin)%array))
               END IF
               IF (atom_composite_reference) THEN
                  CPASSERT(PRESENT(composite_vxc_rho))
                  CPASSERT(PRESENT(composite_vxc_tau))
                  CPASSERT(ASSOCIATED(composite_vxc_rho))
                  CPASSERT(ASSOCIATED(composite_vxc_tau))
                  CPASSERT(SIZE(composite_vxc_rho) == nspins)
                  CPASSERT(SIZE(composite_vxc_tau) == nspins)
                  CALL pw_axpy(smooth_vxc_rho(ispin), composite_vxc_rho(ispin), 1.0_dp)
                  CALL pw_axpy(smooth_vxc_tau(ispin), composite_vxc_tau(ispin), 1.0_dp)
               END IF
               CALL auxbas_pw_pool%give_back_pw(smooth_vxc_rho(ispin))
               CALL auxbas_pw_pool%give_back_pw(smooth_vxc_tau(ispin))
            END DO
            IF (atom_composite_diagnostic) CALL para_env%sum(smooth_grid_contraction)
            DEALLOCATE (smooth_vxc_rho, smooth_vxc_tau)

            one_center_field_contraction = 0.0_dp
            one_center_matrix_contraction = 0.0_dp
            one_center_density_field_contraction = 0.0_dp
            one_center_density_matrix_contraction = 0.0_dp
            one_center_gradient_field_contraction = 0.0_dp
            one_center_gradient_matrix_contraction = 0.0_dp
            one_center_rho_grad_field_contraction = 0.0_dp
            one_center_rho_grad_matrix_contraction = 0.0_dp
            one_center_tau_field_contraction = 0.0_dp
            one_center_tau_matrix_contraction = 0.0_dp
            IF (.NOT. direct_valence_atom_composite) THEN
               DO ikind = 1, SIZE(atomic_kind_set)
                  CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
                  NULLIFY (gth_potential, sgp_potential)
                  CALL get_qs_kind(my_kind_set(ikind), paw_atom=paw_atom, &
                                   gth_potential=gth_potential, harmonics=harmonics, &
                                   grid_atom=grid_atom, sgp_potential=sgp_potential, &
                                   zatom=zatom, zeff=zeff)
                  CALL get_qs_kind(my_kind_set(ikind), basis_set=basis_1c, basis_type="GAPW_1C")
                  IF (.NOT. native_skala_uses_one_center_kind( &
                      paw_atom, gapw_representation, &
                      ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential), &
                      zeff, zatom)) CYCLE

                  nr = grid_atom%nr
                  na = grid_atom%ng_sphere
                  CALL reallocate(rho_h, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(rho_s, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(vxc_h, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(vxc_s, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(drho_h, 1, 4, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(drho_s, 1, 4, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(vxg_h, 1, 3, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(vxg_s, 1, 3, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(tau_h, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(tau_s, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(vtau_h, 1, na, 1, nr, 1, nspins)
                  CALL reallocate(vtau_s, 1, na, 1, nr, 1, nspins)
                  CALL create_tau_basis_cache(tau_basis_cache, grid_atom, basis_1c, harmonics)
                  ALLOCATE (adjoint_reduce_buf(10*na*nr*nspins))

                  bo = get_limit(natom, para_env%num_pe, para_env%mepos)
                  DO iat = 1, natom
                     iatom = atom_list(iat)
                     source_matrix_local = iat >= bo(1) .AND. iat <= bo(2)
                     rho_atom => my_rho_atom_set(iatom)
                     NULLIFY (cpc_h, cpc_s, r_h, r_s, dr_h, dr_s, r_h_d, r_s_d, int_hh, int_ss)
                     CALL get_rho_atom(rho_atom=rho_atom, cpc_h=cpc_h, cpc_s=cpc_s, &
                                       rho_rad_h=r_h, rho_rad_s=r_s, drho_rad_h=dr_h, &
                                       drho_rad_s=dr_s, rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d, &
                                       ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
                     rho_h = 0.0_dp
                     rho_s = 0.0_dp
                     drho_h = 0.0_dp
                     drho_s = 0.0_dp
                     CALL calc_tau_atom(tau_h, tau_s, rho_atom, tau_basis_cache, nspins)
                     DO ir = 1, nr
                        CALL calc_rho_angular(grid_atom, harmonics, nspins, .TRUE., &
                                              ir, r_h, r_s, rho_h, rho_s, dr_h, dr_s, &
                                              r_h_d, r_s_d, drho_h, drho_s)
                     END DO

                     vxc_h = 0.0_dp
                     vxc_s = 0.0_dp
                     vxg_h = 0.0_dp
                     vxg_s = 0.0_dp
                     vtau_h = 0.0_dp
                     vtau_s = 0.0_dp
                     IF (source_matrix_local) THEN
                        composite_row = composite_atom_start(iatom) - 1
                        DO ir = 1, nr
                           DO ia = 1, na
                              composite_row = composite_row + 1
                              IF (lsd) THEN
                                 IF (use_atom_composite_density) THEN
                                    vxc_h(ia, ir, 1:2) = composite_density_grad(composite_row, 1:2)
                                    vxc_s(ia, ir, 1:2) = composite_density_grad(composite_row, 1:2)
                                 END IF
                                 IF (use_atom_composite_gradient) THEN
                                    DO idir = 1, 3
                                       vxg_h(idir, ia, ir, 1:2) = &
                                          composite_grad_grad(composite_row, idir, 1:2)
                                       vxg_s(idir, ia, ir, 1:2) = &
                                          composite_grad_grad(composite_row, idir, 1:2)
                                    END DO
                                 END IF
                                 IF (use_atom_composite_tau) THEN
                                    vtau_h(ia, ir, 1:2) = composite_kin_grad(composite_row, 1:2)
                                    vtau_s(ia, ir, 1:2) = composite_kin_grad(composite_row, 1:2)
                                 END IF
                              ELSE
                                 IF (use_atom_composite_density) THEN
                                    vxc_h(ia, ir, 1) = 0.5_dp* &
                                                       SUM(composite_density_grad(composite_row, :))
                                    vxc_s(ia, ir, 1) = vxc_h(ia, ir, 1)
                                 END IF
                                 IF (use_atom_composite_gradient) THEN
                                    DO idir = 1, 3
                                       vxg_h(idir, ia, ir, 1) = 0.5_dp* &
                                                                SUM(composite_grad_grad( &
                                                                    composite_row, idir, :))
                                       vxg_s(idir, ia, ir, 1) = vxg_h(idir, ia, ir, 1)
                                    END DO
                                 END IF
                                 IF (use_atom_composite_tau) THEN
                                    vtau_h(ia, ir, 1) = 0.5_dp* &
                                                        SUM(composite_kin_grad(composite_row, :))
                                    vtau_s(ia, ir, 1) = vtau_h(ia, ir, 1)
                                 END IF
                              END IF
                           END DO
                        END DO
                        CPASSERT(composite_row == composite_atom_end(iatom))
                     END IF

                     cross_cutoff = gapw_atom_grid_support_radius( &
                                    grid_atom, rho_h, rho_s, drho_h, drho_s, tau_h, tau_s)
                     IF (cross_cutoff > 0.0_dp) THEN
                        image_shell = 0
                        DO idir = 1, 3
                           IF (cell%perd(idir) == 1) THEN
                              image_shell(idir) = CEILING( &
                                                  cross_cutoff*SQRT(SUM(cell%h_inv(idir, :)**2))) + 1
                           END IF
                        END DO
                        ! Each thread accumulates its own one-center potentials; combine them
                        ! before the existing MPI reduction over target-atom owners.
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(image_lower, image_upper, composite_row, cross_density_adjoint, cross_displacement, &
!$OMP         cross_grad_adjoint, cross_kin_adjoint, fractional, idir, image_i1, &
!$OMP         image_i2, image_i3, image_shift, image_translation, jdir, target_atom, &
!$OMP         vxc_h_local, vxc_s_local, vxg_h_local, vxg_s_local, vtau_h_local, vtau_s_local) &
!$OMP SHARED(cell, composite_density_grad, composite_grad_grad, composite_grid_atom, &
!$OMP        composite_grid_coords, composite_kin_grad, composite_nflat, cross_cutoff, &
!$OMP        grid_atom, harmonics, iatom, image_shell, lsd, nspins, particle_set, &
!$OMP        use_atom_composite_density, use_atom_composite_gradient, use_atom_composite_tau, &
!$OMP        vxc_h, vxc_s, vxg_h, vxg_s, vtau_h, vtau_s)
                        ALLOCATE (vxc_h_local, MOLD=vxc_h)
                        ALLOCATE (vxc_s_local, MOLD=vxc_s)
                        ALLOCATE (vxg_h_local, MOLD=vxg_h)
                        ALLOCATE (vxg_s_local, MOLD=vxg_s)
                        ALLOCATE (vtau_h_local, MOLD=vtau_h)
                        ALLOCATE (vtau_s_local, MOLD=vtau_s)
                        vxc_h_local = 0.0_dp
                        vxc_s_local = 0.0_dp
                        vxg_h_local = 0.0_dp
                        vxg_s_local = 0.0_dp
                        vtau_h_local = 0.0_dp
                        vtau_s_local = 0.0_dp
!$OMP DO SCHEDULE(STATIC)
                        DO composite_row = 1, composite_nflat
                           target_atom = composite_grid_atom(composite_row)
                           cross_density_adjoint = 0.0_dp
                           cross_grad_adjoint = 0.0_dp
                           cross_kin_adjoint = 0.0_dp
                           IF (lsd) THEN
                              IF (use_atom_composite_density) THEN
                                 cross_density_adjoint(1:2) = &
                                    composite_density_grad(composite_row, 1:2)
                              END IF
                              IF (use_atom_composite_gradient) THEN
                                 cross_grad_adjoint(:, 1:2) = &
                                    composite_grad_grad(composite_row, :, 1:2)
                              END IF
                              IF (use_atom_composite_tau) THEN
                                 cross_kin_adjoint(1:2) = &
                                    composite_kin_grad(composite_row, 1:2)
                              END IF
                           ELSE
                              IF (use_atom_composite_density) THEN
                                 cross_density_adjoint(1) = 0.5_dp* &
                                                            SUM(composite_density_grad(composite_row, :))
                              END IF
                              IF (use_atom_composite_gradient) THEN
                                 DO idir = 1, 3
                                    cross_grad_adjoint(idir, 1) = 0.5_dp* &
                                                                  SUM(composite_grad_grad(composite_row, idir, :))
                                 END DO
                              END IF
                              IF (use_atom_composite_tau) THEN
                                 cross_kin_adjoint(1) = 0.5_dp* &
                                                        SUM(composite_kin_grad(composite_row, :))
                              END IF
                           END IF
                           fractional = 0.0_dp
                           DO idir = 1, 3
                              DO jdir = 1, 3
                                 fractional(idir) = fractional(idir) + cell%h_inv(idir, jdir)* &
                                                    (composite_grid_coords(jdir, composite_row) - &
                                                     particle_set(iatom)%r(jdir))
                              END DO
                           END DO
                           CALL atom_grid_image_bounds(cell, fractional, cross_cutoff, image_shell, &
                                                       image_lower, image_upper)
                           DO image_i3 = image_lower(3), image_upper(3)
                              DO image_i2 = image_lower(2), image_upper(2)
                                 DO image_i1 = image_lower(1), image_upper(1)
                                    image_shift = [image_i1, image_i2, image_i3]
                                    IF (target_atom == iatom .AND. ALL(image_shift == 0)) CYCLE
                                    image_translation = MATMUL( &
                                                        cell%hmat, REAL(image_shift, dp))
                                    cross_displacement = &
                                       composite_grid_coords(:, composite_row) - &
                                       particle_set(iatom)%r - image_translation
                                    ! Reject inactive images before allocating interpolation scratch.
                                    IF (SQRT(SUM(cross_displacement**2)) > cross_cutoff) CYCLE
                                    CALL add_gapw_atom_grid_interpolation_adjoint( &
                                       grid_atom, harmonics, cross_displacement, cross_cutoff, &
                                       nspins, cross_density_adjoint, cross_grad_adjoint, &
                                       cross_kin_adjoint, vxc_h_local, vxc_s_local, vxg_h_local, vxg_s_local, &
                                       vtau_h_local, vtau_s_local)
                                 END DO
                              END DO
                           END DO
                        END DO
!$OMP END DO
!$OMP CRITICAL(skala_atom_composite_adjoint_reduction)
                        vxc_h = vxc_h + vxc_h_local
                        vxc_s = vxc_s + vxc_s_local
                        vxg_h = vxg_h + vxg_h_local
                        vxg_s = vxg_s + vxg_s_local
                        vtau_h = vtau_h + vtau_h_local
                        vtau_s = vtau_s + vtau_s_local
!$OMP END CRITICAL(skala_atom_composite_adjoint_reduction)
                        DEALLOCATE (vxc_h_local, vxc_s_local, vxg_h_local, vxg_s_local, vtau_h_local, vtau_s_local)
!$OMP END PARALLEL
                     END IF

                     ! Model rows are distributed by target atom, while CP2K stores each
                     ! one-center matrix on the rank owning its source atom. Reduce the exact
                     ! interpolation transpose to that owner in one collective (only the
                     ! owner consumes the result), packing the six arrays into one buffer.
                     owner_rank = para_env%mepos
                     DO row = 0, para_env%num_pe - 1
                        owner_bounds = get_limit(natom, para_env%num_pe, row)
                        IF (iat >= owner_bounds(1) .AND. iat <= owner_bounds(2)) THEN
                           owner_rank = row
                           EXIT
                        END IF
                     END DO
                     adjoint_nn = na*nr*nspins
                     adjoint_off = 0
                     adjoint_reduce_buf(adjoint_off + 1:adjoint_off + adjoint_nn) = &
                        RESHAPE(vxc_h, [adjoint_nn])
                     adjoint_off = adjoint_off + adjoint_nn
                     adjoint_reduce_buf(adjoint_off + 1:adjoint_off + adjoint_nn) = &
                        RESHAPE(vxc_s, [adjoint_nn])
                     adjoint_off = adjoint_off + adjoint_nn
                     adjoint_reduce_buf(adjoint_off + 1:adjoint_off + 3*adjoint_nn) = &
                        RESHAPE(vxg_h, [3*adjoint_nn])
                     adjoint_off = adjoint_off + 3*adjoint_nn
                     adjoint_reduce_buf(adjoint_off + 1:adjoint_off + 3*adjoint_nn) = &
                        RESHAPE(vxg_s, [3*adjoint_nn])
                     adjoint_off = adjoint_off + 3*adjoint_nn
                     adjoint_reduce_buf(adjoint_off + 1:adjoint_off + adjoint_nn) = &
                        RESHAPE(vtau_h, [adjoint_nn])
                     adjoint_off = adjoint_off + adjoint_nn
                     adjoint_reduce_buf(adjoint_off + 1:adjoint_off + adjoint_nn) = &
                        RESHAPE(vtau_s, [adjoint_nn])
                     CALL para_env%sum(adjoint_reduce_buf, owner_rank)
                     IF (.NOT. source_matrix_local) CYCLE
                     adjoint_off = 0
                     vxc_h = RESHAPE(adjoint_reduce_buf(adjoint_off + 1:adjoint_off + adjoint_nn), &
                                     SHAPE(vxc_h))
                     adjoint_off = adjoint_off + adjoint_nn
                     vxc_s = RESHAPE(adjoint_reduce_buf(adjoint_off + 1:adjoint_off + adjoint_nn), &
                                     SHAPE(vxc_s))
                     adjoint_off = adjoint_off + adjoint_nn
                     vxg_h = RESHAPE(adjoint_reduce_buf(adjoint_off + 1:adjoint_off + 3*adjoint_nn), &
                                     SHAPE(vxg_h))
                     adjoint_off = adjoint_off + 3*adjoint_nn
                     vxg_s = RESHAPE(adjoint_reduce_buf(adjoint_off + 1:adjoint_off + 3*adjoint_nn), &
                                     SHAPE(vxg_s))
                     adjoint_off = adjoint_off + 3*adjoint_nn
                     vtau_h = RESHAPE(adjoint_reduce_buf(adjoint_off + 1:adjoint_off + adjoint_nn), &
                                      SHAPE(vtau_h))
                     adjoint_off = adjoint_off + adjoint_nn
                     vtau_s = RESHAPE(adjoint_reduce_buf(adjoint_off + 1:adjoint_off + adjoint_nn), &
                                      SHAPE(vtau_s))

                     IF (atom_composite_diagnostic) THEN
                        one_center_rho_grad_field_contraction = &
                           one_center_rho_grad_field_contraction + SUM(vxc_h*(rho_h - rho_s)) + &
                           SUM(vxg_h*(drho_h(1:3, :, :, :) - drho_s(1:3, :, :, :)))
                        one_center_density_field_contraction = one_center_density_field_contraction + &
                                                               SUM(vxc_h*(rho_h - rho_s))
                        one_center_gradient_field_contraction = one_center_gradient_field_contraction + &
                                                                SUM(vxg_h*(drho_h(1:3, :, :, :) - &
                                                                           drho_s(1:3, :, :, :)))
                        one_center_tau_field_contraction = one_center_tau_field_contraction + &
                                                           SUM(vtau_h*(tau_h - tau_s))
                     END IF

                     ALLOCATE (composite_int_h(SIZE(int_hh(1)%r_coef, 1), &
                                               SIZE(int_hh(1)%r_coef, 2), nspins), &
                               composite_int_s(SIZE(int_ss(1)%r_coef, 1), &
                                               SIZE(int_ss(1)%r_coef, 2), nspins))
                     DO ispin = 1, nspins
                        composite_int_h(:, :, ispin) = int_hh(ispin)%r_coef
                        composite_int_s(:, :, ispin) = int_ss(ispin)%r_coef
                        int_hh(ispin)%r_coef = 0.0_dp
                        int_ss(ispin)%r_coef = 0.0_dp
                     END DO
                     CALL gaVxcgb_noGC(vxc_h, vxc_s, int_hh, int_ss, &
                                       grid_atom, basis_1c, harmonics, nspins)
                     DO ispin = 1, nspins
                        IF (atom_composite_diagnostic) THEN
                           one_center_density_matrix_contraction = &
                              one_center_density_matrix_contraction + &
                              contract_one_center_matrix(cpc_h(ispin)%r_coef, &
                                                         int_hh(ispin)%r_coef, &
                                                         tau_basis_cache%n2oindex) - &
                              contract_one_center_matrix(cpc_s(ispin)%r_coef, &
                                                         int_ss(ispin)%r_coef, &
                                                         tau_basis_cache%n2oindex)
                        END IF
                        int_hh(ispin)%r_coef = 0.0_dp
                        int_ss(ispin)%r_coef = 0.0_dp
                     END DO
                     CALL gaVxcgb_GC(vxc_h, vxc_s, vxg_h, vxg_s, int_hh, int_ss, &
                                     grid_atom, basis_1c, harmonics, nspins)
                     DO ispin = 1, nspins
                        IF (atom_composite_diagnostic) THEN
                           one_center_rho_grad_matrix_contraction = &
                              one_center_rho_grad_matrix_contraction + &
                              contract_one_center_matrix(cpc_h(ispin)%r_coef, &
                                                         int_hh(ispin)%r_coef, &
                                                         tau_basis_cache%n2oindex) - &
                              contract_one_center_matrix(cpc_s(ispin)%r_coef, &
                                                         int_ss(ispin)%r_coef, &
                                                         tau_basis_cache%n2oindex)
                        END IF
                     END DO
                     CALL dgaVtaudgb(vtau_h, vtau_s, int_hh, int_ss, tau_basis_cache, nspins)
                     DO ispin = 1, nspins
                        IF (atom_composite_diagnostic) THEN
                           one_center_matrix_contraction = one_center_matrix_contraction + &
                                                           contract_one_center_matrix(cpc_h(ispin)%r_coef, &
                                                                                      int_hh(ispin)%r_coef, &
                                                                                      tau_basis_cache%n2oindex) - &
                                                           contract_one_center_matrix(cpc_s(ispin)%r_coef, &
                                                                                      int_ss(ispin)%r_coef, &
                                                                                      tau_basis_cache%n2oindex)
                        END IF
                        IF (.NOT. atom_composite_reference) THEN
                           int_hh(ispin)%r_coef = composite_int_h(:, :, ispin)
                           int_ss(ispin)%r_coef = composite_int_s(:, :, ispin)
                        END IF
                     END DO
                     DEALLOCATE (composite_int_h, composite_int_s)
                  END DO
                  DEALLOCATE (adjoint_reduce_buf)
                  CALL release_tau_basis_cache(tau_basis_cache)
               END DO
            END IF
            IF (atom_composite_diagnostic) THEN
               CALL para_env%sum(one_center_density_field_contraction)
               CALL para_env%sum(one_center_density_matrix_contraction)
               CALL para_env%sum(one_center_gradient_field_contraction)
               CALL para_env%sum(one_center_matrix_contraction)
               CALL para_env%sum(one_center_rho_grad_field_contraction)
               CALL para_env%sum(one_center_rho_grad_matrix_contraction)
               CALL para_env%sum(one_center_tau_field_contraction)
               one_center_field_contraction = one_center_rho_grad_field_contraction + &
                                              one_center_tau_field_contraction
               one_center_gradient_matrix_contraction = one_center_rho_grad_matrix_contraction - &
                                                        one_center_density_matrix_contraction
               one_center_tau_matrix_contraction = one_center_matrix_contraction - &
                                                   one_center_rho_grad_matrix_contraction
               feature_component_analytic(1) = SUM(composite_density_grad*composite_density)
               feature_component_analytic(2) = SUM(composite_grad_grad*composite_grad)
               feature_component_analytic(3) = SUM(composite_kin_grad*composite_kin)
               feature_component_analytic(4) = SUM(feature_component_analytic(1:3))
               CALL para_env%sum(feature_component_analytic)
               one_center_tensor_contraction = feature_component_analytic(4) - &
                                               smooth_input_contraction
               DO icomponent = 1, 4
                  SELECT CASE (icomponent)
                  CASE (1)
                     composite_density = (1.0_dp + feature_vxc_step)*composite_density
                  CASE (2)
                     composite_grad = (1.0_dp + feature_vxc_step)*composite_grad
                  CASE (3)
                     composite_kin = (1.0_dp + feature_vxc_step)*composite_kin
                  CASE (4)
                     composite_density = (1.0_dp + feature_vxc_step)*composite_density
                     composite_grad = (1.0_dp + feature_vxc_step)*composite_grad
                     composite_kin = (1.0_dp + feature_vxc_step)*composite_kin
                  END SELECT
                  CALL skala_gapw_atom_composite_energy( &
                     my_xc_section, para_env, composite_density, composite_grad, composite_kin, &
                     composite_grid_coords, composite_grid_weights, composite_base_grid_weights, &
                     composite_local_grid_sizes, composite_local_atom_coords, feature_vxc_plus, &
                     static_geometry_valid=comp_geometry_valid)
                  SELECT CASE (icomponent)
                  CASE (1)
                     composite_density = ((1.0_dp - feature_vxc_step)/ &
                                          (1.0_dp + feature_vxc_step))*composite_density
                  CASE (2)
                     composite_grad = ((1.0_dp - feature_vxc_step)/ &
                                       (1.0_dp + feature_vxc_step))*composite_grad
                  CASE (3)
                     composite_kin = ((1.0_dp - feature_vxc_step)/ &
                                      (1.0_dp + feature_vxc_step))*composite_kin
                  CASE (4)
                     composite_density = ((1.0_dp - feature_vxc_step)/ &
                                          (1.0_dp + feature_vxc_step))*composite_density
                     composite_grad = ((1.0_dp - feature_vxc_step)/ &
                                       (1.0_dp + feature_vxc_step))*composite_grad
                     composite_kin = ((1.0_dp - feature_vxc_step)/ &
                                      (1.0_dp + feature_vxc_step))*composite_kin
                  END SELECT
                  CALL skala_gapw_atom_composite_energy( &
                     my_xc_section, para_env, composite_density, composite_grad, composite_kin, &
                     composite_grid_coords, composite_grid_weights, composite_base_grid_weights, &
                     composite_local_grid_sizes, composite_local_atom_coords, feature_vxc_minus, &
                     static_geometry_valid=comp_geometry_valid)
                  SELECT CASE (icomponent)
                  CASE (1)
                     composite_density = composite_density/(1.0_dp - feature_vxc_step)
                  CASE (2)
                     composite_grad = composite_grad/(1.0_dp - feature_vxc_step)
                  CASE (3)
                     composite_kin = composite_kin/(1.0_dp - feature_vxc_step)
                  CASE (4)
                     composite_density = composite_density/(1.0_dp - feature_vxc_step)
                     composite_grad = composite_grad/(1.0_dp - feature_vxc_step)
                     composite_kin = composite_kin/(1.0_dp - feature_vxc_step)
                  END SELECT
                  feature_component_fd(icomponent) = &
                     (feature_vxc_plus - feature_vxc_minus)/(2.0_dp*feature_vxc_step)
               END DO
               feature_vxc_analytic = feature_component_analytic(4)
               feature_vxc_fd = feature_component_fd(4)
               iw = cp_logger_get_default_io_unit()
               IF (iw > 0) THEN
                  WRITE (UNIT=iw, FMT="(T2,A,1X,I0)") &
                     "SKALA_GPW| Atom-composite reference components", atom_composite_components
                  WRITE (UNIT=iw, FMT="(T2,A,1X,ES24.16)") &
                     "SKALA_GPW| Atom-composite reference electrons", atom_composite_nelec
                  WRITE (UNIT=iw, FMT="(T2,A,1X,ES24.16)") &
                     "SKALA_GPW| Atom-composite reference XC energy", atom_composite_exc
                  WRITE (UNIT=iw, FMT="(T2,A,1X,ES24.16)") &
                     "SKALA_GPW| Atom-composite feature VXC contraction", feature_vxc_analytic
                  WRITE (UNIT=iw, FMT="(T2,A,1X,ES24.16)") &
                     "SKALA_GPW| Atom-composite feature finite difference", feature_vxc_fd
                  WRITE (UNIT=iw, FMT="(T2,A,1X,ES24.16)") &
                     "SKALA_GPW| Atom-composite feature VXC error", &
                     feature_vxc_analytic - feature_vxc_fd
                  WRITE (UNIT=iw, FMT="(T2,A,3(1X,ES24.16))") &
                     "SKALA_GPW| Atom-composite smooth PW adjoint", smooth_input_contraction, &
                     smooth_grid_contraction, smooth_grid_contraction - smooth_input_contraction
                  WRITE (UNIT=iw, FMT="(T2,A,4(1X,ES24.16))") &
                     "SKALA_GPW| Atom-composite one-center adjoint", one_center_tensor_contraction, &
                     one_center_field_contraction, one_center_matrix_contraction, &
                     one_center_matrix_contraction - one_center_tensor_contraction
                  WRITE (UNIT=iw, FMT="(T2,A,4(1X,ES24.16))") &
                     "SKALA_GPW| Atom-composite one-center channels", &
                     one_center_rho_grad_field_contraction, one_center_rho_grad_matrix_contraction, &
                     one_center_tau_field_contraction, one_center_tau_matrix_contraction
                  WRITE (UNIT=iw, FMT="(T2,A,4(1X,ES24.16))") &
                     "SKALA_GPW| Atom-composite one-center rho-gradient", &
                     one_center_density_field_contraction, one_center_density_matrix_contraction, &
                     one_center_gradient_field_contraction, one_center_gradient_matrix_contraction
                  DO icomponent = 1, 3
                     WRITE (UNIT=iw, FMT="(T2,A,1X,I0,3(1X,ES24.16))") &
                        "SKALA_GPW| Atom-composite component VXC", icomponent, &
                        feature_component_analytic(icomponent), feature_component_fd(icomponent), &
                        feature_component_analytic(icomponent) - feature_component_fd(icomponent)
                  END DO
               END IF
            END IF
            IF (atom_composite_reference) THEN
               exc1 = atom_composite_exc
               IF (native_grid_diagnostics) THEN
                  IF (composite_nflat > 0) THEN
                     composite_density_min = MINVAL(composite_density)
                     composite_density_max = MAXVAL(composite_density)
                     composite_kin_min = MINVAL(composite_kin)
                     composite_kin_max = MAXVAL(composite_kin)
                     composite_grad_max = MAXVAL(ABS(composite_grad))
                  ELSE
                     composite_density_min = HUGE(1.0_dp)
                     composite_density_max = -HUGE(1.0_dp)
                     composite_kin_min = HUGE(1.0_dp)
                     composite_kin_max = -HUGE(1.0_dp)
                     composite_grad_max = 0.0_dp
                  END IF
                  composite_tau_integral = &
                     SUM(composite_grid_weights*SUM(composite_kin, DIM=2))
                  CALL para_env%min(composite_density_min)
                  CALL para_env%max(composite_density_max)
                  CALL para_env%min(composite_kin_min)
                  CALL para_env%max(composite_kin_max)
                  CALL para_env%max(composite_grad_max)
                  CALL para_env%sum(composite_tau_integral)
               END IF
               iw = cp_logger_get_default_io_unit()
               IF (iw > 0) THEN
                  WRITE (UNIT=iw, FMT="(T2,A,1X,ES24.16)") &
                     "SKALA_GPW| Active atom-composite XC energy", atom_composite_exc
                  IF (native_grid_diagnostics) THEN
                     WRITE (UNIT=iw, FMT="(T2,A,1X,ES24.16)") &
                        "SKALA_GPW| Active atom-composite electrons", atom_composite_nelec
                     WRITE (UNIT=iw, FMT="(T2,A,2(1X,ES24.16))") &
                        "SKALA_GPW| Active atom-composite density range", &
                        composite_density_min, composite_density_max
                     WRITE (UNIT=iw, FMT="(T2,A,2(1X,ES24.16))") &
                        "SKALA_GPW| Active atom-composite tau range", &
                        composite_kin_min, composite_kin_max
                     WRITE (UNIT=iw, FMT="(T2,A,1X,ES24.16)") &
                        "SKALA_GPW| Active atom-composite tau integral", &
                        composite_tau_integral
                     WRITE (UNIT=iw, FMT="(T2,A,1X,ES24.16)") &
                        "SKALA_GPW| Active atom-composite max gradient", &
                        composite_grad_max
                  END IF
               END IF
            END IF
            CALL xc_rho_set_release(smooth_rho_set, pw_pool=auxbas_pw_pool)
            DEALLOCATE (composite_atomic_grid_sizes, composite_atom_kind, &
                        composite_atom_kind_index, composite_atom_start, composite_atom_end, &
                        composite_atom_coords, composite_local_atoms, composite_local_grid_sizes, &
                        composite_local_atom_coords, composite_grid_atom, &
                        composite_density, composite_grad, composite_kin, &
                        composite_density_grad, composite_grad_grad, composite_kin_grad, &
                        composite_grid_coords, composite_grid_weights, &
                        composite_base_grid_weights)
            IF (ALLOCATED(composite_smooth_density_cache)) THEN
               DEALLOCATE (composite_smooth_density_cache, composite_smooth_gradient_cache, &
                           composite_smooth_kin_cache)
            END IF
         END IF

         IF (.NOT. atom_composite_reference) CALL para_env%sum(exc1)

         IF (ASSOCIATED(rho_h)) DEALLOCATE (rho_h)
         IF (ASSOCIATED(rho_s)) DEALLOCATE (rho_s)
         IF (ASSOCIATED(vxc_h)) DEALLOCATE (vxc_h)
         IF (ASSOCIATED(vxc_s)) DEALLOCATE (vxc_s)

         IF (gradient_f) THEN
            IF (ASSOCIATED(drho_h)) DEALLOCATE (drho_h)
            IF (ASSOCIATED(drho_s)) DEALLOCATE (drho_s)
            IF (ASSOCIATED(vxg_h)) DEALLOCATE (vxg_h)
            IF (ASSOCIATED(vxg_s)) DEALLOCATE (vxg_s)
         END IF

         IF (tau_f) THEN
            IF (ASSOCIATED(tau_h)) DEALLOCATE (tau_h)
            IF (ASSOCIATED(tau_s)) DEALLOCATE (tau_s)
            IF (ASSOCIATED(vtau_h)) DEALLOCATE (vtau_h)
            IF (ASSOCIATED(vtau_s)) DEALLOCATE (vtau_s)
         END IF

      END IF !xc_none

      CALL timestop(handle)

   END SUBROUTINE calculate_vxc_atom

! **************************************************************************************************
!> \brief Add the GAPW one-center correction to CDFT values and operators.
!> \param qs_env Quickstep environment
!> \param energy_only skip construction of the CDFT one-center operator
!> \param calculate_forces evaluate explicit derivatives of the partition weights
!> \param values constraint values from the hard-minus-soft one-center densities
!> \param electronic_charge optional one-center corrections to atomic populations
!> \param operator_group optional group for which to build the unscaled weight operator
!> \param rho_atom_operator_set optional destination for the one-center operator integrals
! **************************************************************************************************
   SUBROUTINE gapw_cdft_one_center(qs_env, energy_only, calculate_forces, values, &
                                   electronic_charge, operator_group, rho_atom_operator_set)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: energy_only, calculate_forces
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: values
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: electronic_charge
      INTEGER, INTENT(IN), OPTIONAL                      :: operator_group
      TYPE(rho_atom_type), DIMENSION(:), POINTER, &
         OPTIONAL                                        :: rho_atom_operator_set

      INTEGER                                            :: atom, channel, ia, iat, igroup, ikind, &
                                                            ir, natom, natom_kind, nspins
      INTEGER, DIMENSION(2)                              :: atom_bounds
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: lsd, paw_atom
      REAL(KIND=dp)                                      :: delta_density, point_factor, spin_factor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: atomic_weights, group_weights
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: explicit_derivative, &
                                                            group_point_derivative
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: group_atom_derivative
      REAL(KIND=dp), DIMENSION(3)                        :: point
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho_h, rho_s
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: drho_h, drho_s
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: vlocal_h, vlocal_s
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(cdft_point_context_type)                      :: context
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: dr_h, dr_s, int_hh, int_ss, r_h, r_s
      TYPE(rho_atom_coeff), DIMENSION(:, :), POINTER     :: r_h_d, r_s_d
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: operator_atom_set, rho_atom_set
      TYPE(rho_atom_type), POINTER                       :: operator_atom, rho_atom

      NULLIFY (atom_list, atomic_kind_set, basis_1c, cdft_control, dft_control, force, &
               grid_atom, harmonics, int_hh, int_ss, para_env, particle_set, r_h, r_s, &
               dr_h, dr_s, r_h_d, r_s_d, rho_h, rho_s, drho_h, drho_s, &
               operator_atom, operator_atom_set, rho_atom, rho_atom_set, qs_kind_set, &
               vlocal_h, vlocal_s)
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, dft_control=dft_control, &
                      force=force, natom=natom, para_env=para_env, particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, rho_atom_set=rho_atom_set)
      CPASSERT(ASSOCIATED(atomic_kind_set))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(para_env))
      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(qs_kind_set))
      CPASSERT(ASSOCIATED(rho_atom_set))
      operator_atom_set => rho_atom_set
      IF (PRESENT(rho_atom_operator_set)) operator_atom_set => rho_atom_operator_set
      CPASSERT(ASSOCIATED(operator_atom_set))
      cdft_control => dft_control%qs_control%cdft_control
      CPASSERT(ASSOCIATED(cdft_control))
      nspins = dft_control%nspins
      lsd = dft_control%lsd
      CPASSERT(SIZE(values) == SIZE(cdft_control%group))
      IF (PRESENT(operator_group)) THEN
         CPASSERT(operator_group >= 1 .AND. operator_group <= SIZE(cdft_control%group))
      END IF
      IF (PRESENT(electronic_charge)) THEN
         CPASSERT(SIZE(electronic_charge, 1) == natom)
         CPASSERT(SIZE(electronic_charge, 2) == nspins)
      END IF

      CALL cdft_point_context_create(qs_env, context, calculate_forces)
      ALLOCATE (group_weights(context%ngroup), group_point_derivative(3, context%ngroup), &
                group_atom_derivative(3, natom, context%ngroup), &
                explicit_derivative(3, natom))
      IF (PRESENT(electronic_charge)) ALLOCATE (atomic_weights(natom))
      values = 0.0_dp
      explicit_derivative = 0.0_dp
      IF (PRESENT(electronic_charge)) electronic_charge = 0.0_dp

      DO ikind = 1, SIZE(atomic_kind_set)
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom_kind)
         CALL get_qs_kind(qs_kind_set(ikind), paw_atom=paw_atom, grid_atom=grid_atom, &
                          harmonics=harmonics)
         IF (.NOT. paw_atom) CYCLE
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_1c, basis_type="GAPW_1C")
         CPASSERT(ASSOCIATED(grid_atom))
         CPASSERT(ASSOCIATED(harmonics))
         CPASSERT(ASSOCIATED(basis_1c))
         ALLOCATE (vlocal_h(grid_atom%ng_sphere, grid_atom%nr, nspins), &
                   vlocal_s(grid_atom%ng_sphere, grid_atom%nr, nspins))

         atom_bounds = get_limit(natom_kind, para_env%num_pe, para_env%mepos)
         DO iat = atom_bounds(1), atom_bounds(2)
            atom = atom_list(iat)
            rho_atom => rho_atom_set(atom)
            NULLIFY (r_h, r_s)
            CALL get_rho_atom(rho_atom=rho_atom, rho_rad_h=r_h, rho_rad_s=r_s)
            CALL reallocate(rho_h, 1, grid_atom%ng_sphere, 1, grid_atom%nr, 1, nspins)
            CALL reallocate(rho_s, 1, grid_atom%ng_sphere, 1, grid_atom%nr, 1, nspins)
            rho_h = 0.0_dp
            rho_s = 0.0_dp
            DO ir = 1, grid_atom%nr
               CALL calc_rho_angular(grid_atom, harmonics, nspins, .FALSE., ir, r_h, r_s, &
                                     rho_h, rho_s, dr_h, dr_s, r_h_d, r_s_d, drho_h, drho_s)
            END DO
            vlocal_h = 0.0_dp
            vlocal_s = 0.0_dp

            DO ir = 1, grid_atom%nr
               DO ia = 1, grid_atom%ng_sphere
                  point(1) = particle_set(atom)%r(1) + grid_atom%rad(ir)* &
                             grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia)
                  point(2) = particle_set(atom)%r(2) + grid_atom%rad(ir)* &
                             grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia)
                  point(3) = particle_set(atom)%r(3) + grid_atom%rad(ir)*grid_atom%cos_pol(ia)
                  IF (PRESENT(electronic_charge)) THEN
                     CALL cdft_point_weights(context, point, group_weights, &
                                             group_point_derivative, group_atom_derivative, &
                                             atomic_weights)
                  ELSE
                     CALL cdft_point_weights(context, point, group_weights, &
                                             group_point_derivative, group_atom_derivative)
                  END IF
                  DO channel = 1, nspins
                     delta_density = rho_h(ia, ir, channel) - rho_s(ia, ir, channel)
                     IF (PRESENT(electronic_charge)) THEN
                        electronic_charge(:, channel) = electronic_charge(:, channel) + &
                                                        grid_atom%weight(ia, ir)*atomic_weights* &
                                                        delta_density
                     END IF
                     point_factor = 0.0_dp
                     DO igroup = 1, context%ngroup
                        spin_factor = cdft_spin_factor( &
                                      cdft_control%group(igroup)%constraint_type, channel, lsd)
                        IF (PRESENT(operator_group)) THEN
                           IF (igroup == operator_group) THEN
                              point_factor = point_factor + group_weights(igroup)
                           END IF
                        ELSE
                           point_factor = point_factor + cdft_control%strength(igroup)* &
                                          group_weights(igroup)*spin_factor
                        END IF
                        values(igroup) = values(igroup) + grid_atom%weight(ia, ir)* &
                                         group_weights(igroup)*delta_density*spin_factor
                        IF (calculate_forces) THEN
                           explicit_derivative(:, :) = &
                              explicit_derivative + grid_atom%weight(ia, ir)* &
                              cdft_control%strength(igroup)*delta_density* &
                              spin_factor*group_atom_derivative(:, :, igroup)
                           explicit_derivative(:, atom) = explicit_derivative(:, atom) + &
                                                          grid_atom%weight(ia, ir)* &
                                                          cdft_control%strength(igroup)*delta_density* &
                                                          spin_factor*group_point_derivative(:, igroup)
                        END IF
                     END DO
                     vlocal_h(ia, ir, channel) = grid_atom%weight(ia, ir)*point_factor
                     vlocal_s(ia, ir, channel) = vlocal_h(ia, ir, channel)
                  END DO
               END DO
            END DO

            IF (.NOT. energy_only) THEN
               NULLIFY (int_hh, int_ss)
               operator_atom => operator_atom_set(atom)
               CALL get_rho_atom(rho_atom=operator_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
               CALL gaVxcgb_noGC(vlocal_h, vlocal_s, int_hh, int_ss, grid_atom, &
                                 basis_1c, harmonics, nspins)
            END IF
         END DO
         DEALLOCATE (vlocal_h, vlocal_s)
      END DO

      CALL para_env%sum(values)
      CALL para_env%sum(explicit_derivative)
      IF (PRESENT(electronic_charge)) CALL para_env%sum(electronic_charge)
      IF (calculate_forces .AND. ASSOCIATED(force) .AND. para_env%is_source()) THEN
         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom_kind)
            DO iat = 1, natom_kind
               atom = atom_list(iat)
               force(ikind)%rho_elec(:, iat) = force(ikind)%rho_elec(:, iat) + &
                                               explicit_derivative(:, atom)
            END DO
         END DO
      END IF

      IF (ASSOCIATED(rho_h)) DEALLOCATE (rho_h)
      IF (ASSOCIATED(rho_s)) DEALLOCATE (rho_s)
      IF (ALLOCATED(atomic_weights)) DEALLOCATE (atomic_weights)
      DEALLOCATE (explicit_derivative, group_atom_derivative, group_point_derivative, group_weights)
      CALL cdft_point_context_release(context)

   CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param constraint_type ...
!> \param channel ...
!> \param lsd ...
!> \return ...
! **************************************************************************************************
      FUNCTION cdft_spin_factor(constraint_type, channel, lsd) RESULT(factor)
      INTEGER, INTENT(IN)                                :: constraint_type, channel
      LOGICAL, INTENT(IN)                                :: lsd
      REAL(KIND=dp)                                      :: factor

         SELECT CASE (constraint_type)
         CASE (cdft_charge_constraint)
            factor = 1.0_dp
         CASE (cdft_magnetization_constraint)
            CPASSERT(lsd)
            factor = MERGE(1.0_dp, -1.0_dp, channel == 1)
         CASE (cdft_alpha_constraint)
            CPASSERT(lsd)
            factor = MERGE(1.0_dp, 0.0_dp, channel == 1)
         CASE (cdft_beta_constraint)
            CPASSERT(lsd)
            factor = MERGE(1.0_dp, 0.0_dp, channel == 2)
         CASE DEFAULT
            CPABORT("Unknown CDFT constraint type.")
         END SELECT
      END FUNCTION cdft_spin_factor

   END SUBROUTINE gapw_cdft_one_center

! **************************************************************************************************
!> \brief Contract a compact one-center density matrix with an integral in the padded old basis.
!> \param density_matrix ...
!> \param integral_matrix ...
!> \param new_to_old ...
!> \return ...
! **************************************************************************************************
   FUNCTION contract_one_center_matrix(density_matrix, integral_matrix, new_to_old) RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: density_matrix, integral_matrix
      INTEGER, DIMENSION(:), INTENT(IN)                  :: new_to_old
      REAL(KIND=dp)                                      :: value

      INTEGER                                            :: ibas, jbas, nbas

      nbas = SIZE(density_matrix, 1)
      CPASSERT(SIZE(density_matrix, 2) == nbas)
      CPASSERT(SIZE(new_to_old) >= nbas)
      CPASSERT(MINVAL(new_to_old(1:nbas)) >= 1)
      CPASSERT(MAXVAL(new_to_old(1:nbas)) <= SIZE(integral_matrix, 1))
      CPASSERT(MAXVAL(new_to_old(1:nbas)) <= SIZE(integral_matrix, 2))

      value = 0.0_dp
      DO jbas = 1, nbas
         DO ibas = 1, nbas
            value = value + density_matrix(ibas, jbas)* &
                    integral_matrix(new_to_old(ibas), new_to_old(jbas))
         END DO
      END DO

   END FUNCTION contract_one_center_matrix

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE native_need_index_release()
      IF (ALLOCATED(skala_state%native_atom%native_need_lin)) DEALLOCATE (skala_state%native_atom%native_need_lin)
      IF (ALLOCATED(skala_state%native_atom%native_need_bucket_entry)) DEALLOCATE (skala_state%native_atom%native_need_bucket_entry)
      IF (ALLOCATED(skala_state%native_atom%native_need_bucket_start)) DEALLOCATE (skala_state%native_atom%native_need_bucket_start)
      skala_state%native_atom%native_need_nbucket = 0
   END SUBROUTINE native_need_index_release

! **************************************************************************************************
!> \brief Smallest power of two greater than or equal to n.
!> \param n ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION native_need_next_pow2(n) RESULT(p)
      INTEGER(int_8), INTENT(IN)                         :: n
      INTEGER(int_8)                                     :: p

      p = 1_int_8
      DO WHILE (p < n .AND. p < HUGE(p)/2_int_8)
         p = 2_int_8*p
      END DO
   END FUNCTION native_need_next_pow2

! **************************************************************************************************
!> \brief Bucket of a grid linear index (bit-mixed, power-of-two table).
!> \param lin global grid linear index
!> \param nbucket power-of-two bucket count
!> \return 1-based bucket
! **************************************************************************************************
   PURE FUNCTION native_need_bucket(lin, nbucket) RESULT(b)
      INTEGER(int_8), INTENT(IN)                         :: lin
      INTEGER, INTENT(IN)                                :: nbucket
      INTEGER                                            :: b

      INTEGER(int_8)                                     :: x

      x = IEOR(lin, ISHFT(lin, -21))
      b = INT(IAND(x, INT(nbucket - 1, int_8))) + 1
   END FUNCTION native_need_bucket

! **************************************************************************************************
!> \brief Map a grid linear index to its compact needed-point position (0 if absent).
!> \param lin global grid linear index
!> \return compact position
! **************************************************************************************************
   FUNCTION native_need_lookup(lin) RESULT(pos)
      INTEGER(int_8), INTENT(IN)                         :: lin
      INTEGER                                            :: pos

      INTEGER                                            :: b, k

      pos = 0
      IF (skala_state%native_atom%native_need_nbucket == 0) RETURN
      b = native_need_bucket(lin, skala_state%native_atom%native_need_nbucket)
      DO k = skala_state%native_atom%native_need_bucket_start(b), skala_state%native_atom%native_need_bucket_start(b + 1) - 1
         IF (skala_state%native_atom%native_need_bucket_entry(k)%lin == lin) THEN
            pos = skala_state%native_atom%native_need_bucket_entry(k)%pos
            RETURN
         END IF
      END DO
   END FUNCTION native_need_lookup

! **************************************************************************************************
!> \brief Linear index of a stencil node within the auxiliary PW grid.
!> \param stencil interpolation stencil
!> \param ix x node
!> \param iy y node
!> \param iz z node
!> \param npts grid point counts
!> \return one-based grid linear index
! **************************************************************************************************
   PURE FUNCTION native_need_linear_index(stencil, ix, iy, iz, npts) RESULT(lin)
      TYPE(native_grid_interpolation_stencil_type), &
         INTENT(IN)                                      :: stencil
      INTEGER, INTENT(IN)                                :: ix, iy, iz, npts(3)
      INTEGER(int_8)                                     :: lin

      lin = 1_int_8 + INT(stencil%relative_index(ix, 1), int_8) + INT(npts(1), int_8)*( &
            INT(stencil%relative_index(iy, 2), int_8) + INT(npts(2), int_8)* &
            INT(stencil%relative_index(iz, 3), int_8))
   END FUNCTION native_need_linear_index

! **************************************************************************************************
!> \brief Fetch a cached interpolation stencil, creating (and optionally storing) it on a miss.
!> \param pw_grid auxiliary real-space grid
!> \param cell current cell
!> \param row composite row index
!> \param point grid point
!> \param stencil returned stencil
!> \param store whether to store a newly built stencil in the cache
! **************************************************************************************************
   SUBROUTINE get_native_grid_stencil(pw_grid, cell, row, point, stencil, store)
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: row
      REAL(dp), INTENT(IN)                               :: point(3)
      TYPE(native_grid_interpolation_stencil_type), &
         INTENT(OUT)                                     :: stencil
      LOGICAL, INTENT(IN)                                :: store

      IF (fetch_native_grid_stencil(skala_state%native_grid_cache, row, point, stencil)) RETURN
      CALL create_native_grid_interpolation_stencil(stencil, pw_grid, cell, point, .TRUE.)
      IF (store) CALL store_native_grid_stencil(skala_state%native_grid_cache, row, point, stencil)
   END SUBROUTINE get_native_grid_stencil

! **************************************************************************************************
!> \brief Build the compact set of PW grid points touched by the local atom-composite stencils
!>        and the owner/local-index exchange schedule.
!> \param grid auxiliary real-space grid
!> \param cell current cell
!> \param group communicator over which the PW grid is distributed (usually para_env)
!> \param wrap wrap the auxiliary cell in nonperiodic directions
!> \param nrows number of local atom-composite rows
!> \param points local row coordinates
!> \param atom_coords all atom coordinates (identical on every rank)
! **************************************************************************************************
   SUBROUTINE native_grid_build_needed(grid, cell, group, wrap, nrows, points, atom_coords)
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell

      CLASS(mp_comm_type), INTENT(IN)                    :: group
      LOGICAL, INTENT(IN)                                :: wrap
      INTEGER, INTENT(IN)                                :: nrows
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: points, atom_coords

      INTEGER                                            :: b, capacity, i, ip, j, k, lindex, &
                                                            npe, npts(3), nneed, nbwords, owner, pos, row
      INTEGER(int_8)                                     :: bit, lin, lin0, mask, nbucket8, ngrid8, &
                                                            nx8, ny8, rel8(3), word
      INTEGER(int_8), ALLOCATABLE, DIMENSION(:)          :: cand, grown, present
      INTEGER(int_8), ALLOCATABLE, DIMENSION(:, :)       :: work_base
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bucket_count, cand_local, cand_owner, &
                                                            fill
      LOGICAL                                            :: valid
      REAL(dp), DIMENSION(3)                             :: relative
      TYPE(native_grid_interpolation_stencil_type)       :: stencil

      npts = grid%npts
      BLOCK
         INTEGER :: lim(2)
         lim = get_limit(npts(1), grid%para%group%num_pe_cart(1), grid%para%group%mepos_cart(1))
         lim = lim + grid%bounds(1, 1) - 1
         CPASSERT(lim(1) == grid%bounds_local(1, 1) .AND. lim(2) == grid%bounds_local(2, 1))
         lim = get_limit(npts(2), grid%para%group%num_pe_cart(2), grid%para%group%mepos_cart(2))
         lim = lim + grid%bounds(1, 2) - 1
         CPASSERT(lim(1) == grid%bounds_local(1, 2) .AND. lim(2) == grid%bounds_local(2, 2))
         CPASSERT(grid%bounds_local(1, 3) == grid%bounds(1, 3))
         CPASSERT(grid%bounds_local(2, 3) == grid%bounds(2, 3))
      END BLOCK
      valid = skala_state%native_atom%native_need_ready .AND. skala_state%native_atom%native_need_nrows == nrows
      IF (valid) THEN
  valid = ALL(skala_state%native_atom%native_need_npts == npts) .AND. ALL(skala_state%native_atom%native_need_dh_inv == grid%dh_inv)
         valid = valid .AND. (skala_state%native_atom%native_need_wrap .EQV. wrap)
         ! The exchange schedule is tied to the communicator and local grid
         ! partition, which are not captured by the geometry keys above.
         valid = valid .AND. skala_state%native_atom%native_need_num_pe == group%num_pe
         valid = valid .AND. skala_state%native_atom%native_need_mepos == group%mepos
         valid = valid .AND. skala_state%native_atom%native_need_group_handle == group%get_handle()
         valid = valid .AND. skala_state%native_atom%native_need_mode == grid%para%mode
         valid = valid .AND. ALL(skala_state%native_atom%native_need_bounds_local == grid%bounds_local)
      END IF
      IF (valid) valid = native_geometry_matches(skala_state%native_atom%native_need_geometry, cell, atom_coords)
      ! Only build the per-row fractional base when the cheap keys still match
      ! (to compare against the cached copy) or when a rebuild is required.
      IF (valid) THEN
         ALLOCATE (work_base(3, MAX(nrows, 1)))
         work_base = 0_int_8
         DO row = 1, nrows
            relative = MATMUL(grid%dh_inv, points(:, row))
            work_base(:, row) = INT(FLOOR(relative), int_8)
         END DO
         IF (ALLOCATED(skala_state%native_atom%native_need_base)) THEN
            valid = ALL(SHAPE(skala_state%native_atom%native_need_base) == SHAPE(work_base))
            IF (valid) valid = ALL(skala_state%native_atom%native_need_base == work_base)
         ELSE
            valid = .FALSE.
         END IF
         IF (valid) THEN
            DEALLOCATE (work_base)
            RETURN
         END IF
      END IF
      IF (.NOT. ALLOCATED(work_base)) THEN
         ALLOCATE (work_base(3, MAX(nrows, 1)))
         work_base = 0_int_8
         DO row = 1, nrows
            relative = MATMUL(grid%dh_inv, points(:, row))
            work_base(:, row) = INT(FLOOR(relative), int_8)
         END DO
      END IF

      IF (ALLOCATED(skala_state%native_atom%native_need_local)) DEALLOCATE (skala_state%native_atom%native_need_local)
      IF (ALLOCATED(skala_state%native_atom%native_need_scount)) DEALLOCATE (skala_state%native_atom%native_need_scount)
      IF (ALLOCATED(skala_state%native_atom%native_need_sdispl)) DEALLOCATE (skala_state%native_atom%native_need_sdispl)
      IF (ALLOCATED(skala_state%native_atom%native_need_rcount)) DEALLOCATE (skala_state%native_atom%native_need_rcount)
      IF (ALLOCATED(skala_state%native_atom%native_need_rdispl)) DEALLOCATE (skala_state%native_atom%native_need_rdispl)
      IF (ALLOCATED(skala_state%native_atom%native_need_sidx)) DEALLOCATE (skala_state%native_atom%native_need_sidx)
      IF (ALLOCATED(skala_state%native_atom%native_need_ridx)) DEALLOCATE (skala_state%native_atom%native_need_ridx)
      IF (ALLOCATED(skala_state%native_atom%native_need_base)) DEALLOCATE (skala_state%native_atom%native_need_base)
      IF (ALLOCATED(skala_state%native_atom%native_need_values)) DEALLOCATE (skala_state%native_atom%native_need_values)
      IF (ALLOCATED(skala_state%native_atom%native_need_reply)) DEALLOCATE (skala_state%native_atom%native_need_reply)
      CALL native_need_index_release()
      skala_state%native_atom%native_need_ready = .FALSE.
      skala_state%native_atom%native_need_npts = npts

      ngrid8 = INT(npts(1), int_8)*INT(npts(2), int_8)*INT(npts(3), int_8)
      nbwords = INT((ngrid8 + 63_int_8)/64_int_8)
      ALLOCATE (PRESENT(nbwords))
      present = 0_int_8
      nneed = 0
      capacity = 1024
      ALLOCATE (cand(capacity))
      cand = 0_int_8
      nx8 = INT(npts(1), int_8)
      ny8 = INT(npts(2), int_8)
      DO row = 1, nrows
         CALL create_native_grid_interpolation_stencil(stencil, grid, cell, points(:, row), wrap, &
                                                       indices_only=.TRUE.)
         IF (.NOT. stencil%active) CYCLE
         DO k = 1, native_grid_interp_npts
            IF (.NOT. stencil%valid(k, 3)) CYCLE
            DO j = 1, native_grid_interp_npts
               IF (.NOT. stencil%valid(j, 2)) CYCLE
               DO i = 1, native_grid_interp_npts
                  IF (.NOT. stencil%valid(i, 1)) CYCLE
                  lin = native_need_linear_index(stencil, i, j, k, npts)
                  word = (lin - 1_int_8)/64_int_8 + 1_int_8
                  bit = MOD(lin - 1_int_8, 64_int_8)
                  mask = ISHFT(1_int_8, INT(bit))
                  IF (IAND(PRESENT(word), mask) /= 0_int_8) CYCLE
                  PRESENT(word) = IOR(PRESENT(word), mask)
                  nneed = nneed + 1
                  IF (nneed > capacity) THEN
                     IF (capacity > HUGE(1)/2) THEN
                        CALL cp_abort(__LOCATION__, &
                                      "Native-grid candidate list would exceed the integer range.")
                     END IF
                     ALLOCATE (grown(2*capacity))
                     grown(1:capacity) = cand(1:capacity)
                     DEALLOCATE (cand)
                     CALL MOVE_ALLOC(grown, cand)
                     capacity = 2*capacity
                  END IF
                  cand(nneed) = lin
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE (present)
      skala_state%native_atom%native_need_nneed = nneed

      npe = group%num_pe
      skala_state%native_atom%native_need_replicated = (grid%para%mode == PW_MODE_LOCAL)
      ALLOCATE (skala_state%native_atom%native_need_lin(MAX(nneed, 1)))
      ALLOCATE (skala_state%native_atom%native_need_local(MAX(nneed, 1)))
      ALLOCATE (skala_state%native_atom%native_need_scount(npe), skala_state%native_atom%native_need_sdispl(npe))
      skala_state%native_atom%native_need_scount = 0
      skala_state%native_atom%native_need_lin = 0_int_8

      IF (skala_state%native_atom%native_need_replicated) THEN
         CPASSERT(group%num_pe == 1)
         DO ip = 1, nneed
            skala_state%native_atom%native_need_local(ip) = INT(cand(ip))
            skala_state%native_atom%native_need_lin(ip) = cand(ip)
         END DO
         skala_state%native_atom%native_need_scount(1) = nneed
         skala_state%native_atom%native_need_sdispl(1) = 0
         ALLOCATE (skala_state%native_atom%native_need_rcount(npe), skala_state%native_atom%native_need_rdispl(npe))
         skala_state%native_atom%native_need_rcount(:) = skala_state%native_atom%native_need_scount
         skala_state%native_atom%native_need_rdispl(:) = skala_state%native_atom%native_need_sdispl
         ALLOCATE (skala_state%native_atom%native_need_sidx(MAX(nneed, 1)), skala_state%native_atom%native_need_ridx(MAX(nneed, 1)))
         skala_state%native_atom%native_need_sidx(:) = skala_state%native_atom%native_need_local
         skala_state%native_atom%native_need_ridx(:) = skala_state%native_atom%native_need_local
      ELSE
         ALLOCATE (cand_owner(MAX(nneed, 1)), cand_local(MAX(nneed, 1)))
         DO ip = 1, nneed
            lin0 = cand(ip) - 1_int_8
            rel8(1) = MOD(lin0, nx8)
            lin0 = lin0/nx8
            rel8(2) = MOD(lin0, ny8)
            rel8(3) = lin0/ny8
            i = grid%bounds(1, 1) + INT(rel8(1))
            j = grid%bounds(1, 2) + INT(rel8(2))
            k = grid%bounds(1, 3) + INT(rel8(3))
            CALL pw_grid_map_global_point(grid, i, j, k, owner, lindex)
            cand_owner(ip) = owner
            cand_local(ip) = lindex
            skala_state%native_atom%native_need_scount(owner + 1) = skala_state%native_atom%native_need_scount(owner + 1) + 1
         END DO
         skala_state%native_atom%native_need_sdispl(1) = 0
         DO ip = 2, npe
            skala_state%native_atom%native_need_sdispl(ip) = &
               skala_state%native_atom%native_need_sdispl(ip - 1) + &
               skala_state%native_atom%native_need_scount(ip - 1)
         END DO
         ALLOCATE (fill(npe))
         fill = 0
         DO ip = 1, nneed
            owner = cand_owner(ip)
            pos = skala_state%native_atom%native_need_sdispl(owner + 1) + fill(owner + 1) + 1
            fill(owner + 1) = fill(owner + 1) + 1
            skala_state%native_atom%native_need_local(pos) = cand_local(ip)
            skala_state%native_atom%native_need_lin(pos) = cand(ip)
         END DO
         DEALLOCATE (fill, cand_owner, cand_local)
         ALLOCATE (skala_state%native_atom%native_need_sidx(MAX(nneed, 1)))
         IF (nneed > 0) skala_state%native_atom%native_need_sidx(1:nneed) = skala_state%native_atom%native_need_local(1:nneed)
         ALLOCATE (skala_state%native_atom%native_need_rcount(npe), skala_state%native_atom%native_need_rdispl(npe))
         CALL group%alltoall(skala_state%native_atom%native_need_scount, skala_state%native_atom%native_need_rcount, 1)
         skala_state%native_atom%native_need_rdispl(1) = 0
         DO ip = 2, npe
            skala_state%native_atom%native_need_rdispl(ip) = &
               skala_state%native_atom%native_need_rdispl(ip - 1) + &
               skala_state%native_atom%native_need_rcount(ip - 1)
         END DO
         ALLOCATE (skala_state%native_atom%native_need_ridx(MAX(SUM(skala_state%native_atom%native_need_rcount), 1)))
         CALL group%alltoall(skala_state%native_atom%native_need_sidx, &
                             skala_state%native_atom%native_need_scount, &
                             skala_state%native_atom%native_need_sdispl, &
                             skala_state%native_atom%native_need_ridx, &
                             skala_state%native_atom%native_need_rcount, &
                             skala_state%native_atom%native_need_rdispl)
      END IF
      skala_state%native_atom%native_need_nrecv = SUM(skala_state%native_atom%native_need_rcount)

      nbucket8 = native_need_next_pow2(MAX(2_int_8*INT(nneed, int_8), 64_int_8))
      IF (nbucket8 > INT(HUGE(1), int_8)) THEN
         CALL cp_abort(__LOCATION__, "Native-grid hash table would exceed the integer range.")
      END IF
      skala_state%native_atom%native_need_nbucket = INT(nbucket8)
      ALLOCATE (skala_state%native_atom%native_need_bucket_start(skala_state%native_atom%native_need_nbucket + 1))
      ALLOCATE (skala_state%native_atom%native_need_bucket_entry(MAX(nneed, 1)))
      ALLOCATE (bucket_count(skala_state%native_atom%native_need_nbucket))
      bucket_count = 0
      DO pos = 1, nneed
         b = native_need_bucket(skala_state%native_atom%native_need_lin(pos), skala_state%native_atom%native_need_nbucket)
         bucket_count(b) = bucket_count(b) + 1
      END DO
      skala_state%native_atom%native_need_bucket_start(1) = 1
      DO b = 1, skala_state%native_atom%native_need_nbucket
     skala_state%native_atom%native_need_bucket_start(b + 1) = skala_state%native_atom%native_need_bucket_start(b) + bucket_count(b)
      END DO
      bucket_count(:) = skala_state%native_atom%native_need_bucket_start(1:skala_state%native_atom%native_need_nbucket)
      DO pos = 1, nneed
         b = native_need_bucket(skala_state%native_atom%native_need_lin(pos), skala_state%native_atom%native_need_nbucket)
         skala_state%native_atom%native_need_bucket_entry(bucket_count(b))%lin = skala_state%native_atom%native_need_lin(pos)
         skala_state%native_atom%native_need_bucket_entry(bucket_count(b))%pos = pos
         bucket_count(b) = bucket_count(b) + 1
      END DO
      DEALLOCATE (bucket_count, cand, skala_state%native_atom%native_need_lin)

      ALLOCATE (skala_state%native_atom%native_need_base(SIZE(work_base, 1), SIZE(work_base, 2)))
      skala_state%native_atom%native_need_base(:, :) = work_base
      CALL store_native_geometry_key(skala_state%native_atom%native_need_geometry, cell, atom_coords)
      DEALLOCATE (work_base)
      skala_state%native_atom%native_need_dh_inv = grid%dh_inv
      skala_state%native_atom%native_need_wrap = wrap
      skala_state%native_atom%native_need_nrows = nrows
      skala_state%native_atom%native_need_num_pe = group%num_pe
      skala_state%native_atom%native_need_mepos = group%mepos
      skala_state%native_atom%native_need_group_handle = group%get_handle()
      skala_state%native_atom%native_need_mode = grid%para%mode
      skala_state%native_atom%native_need_bounds_local = grid%bounds_local
      skala_state%native_atom%native_need_ready = .TRUE.

   END SUBROUTINE native_grid_build_needed

! **************************************************************************************************
!> \brief Read one local smooth field only at the local indices requested by other ranks.
!> \param grid auxiliary real-space grid
!> \param slot destination field slot
!> \param nfg number of gathered fields (stride of the reply buffer)
!> \param values local field values indexed by the grid bounds
! **************************************************************************************************
   SUBROUTINE native_grid_extract_points(grid, slot, nfg, values)
      TYPE(pw_grid_type), INTENT(IN)                     :: grid
      INTEGER, INTENT(IN)                                :: slot, nfg
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: values

      INTEGER                                            :: bo(2, 3), i, idx, ix, iy, iz, j, k, &
                                                            lb(3), li, nlx, nly, rem

      bo = grid%bounds_local
      lb = LBOUND(values)
      nlx = bo(2, 1) - bo(1, 1) + 1
      nly = bo(2, 2) - bo(1, 2) + 1
      DO idx = 1, skala_state%native_atom%native_need_nrecv
         li = skala_state%native_atom%native_need_ridx(idx)
         rem = li - 1
         ix = MOD(rem, nlx)
         rem = rem/nlx
         iy = MOD(rem, nly)
         iz = rem/nly
         i = bo(1, 1) + ix
         j = bo(1, 2) + iy
         k = bo(1, 3) + iz
         skala_state%native_atom%native_need_reply((idx - 1)*nfg + slot) = &
            values(lb(1) + i - bo(1, 1), lb(2) + j - bo(1, 2), lb(3) + k - bo(1, 3))
      END DO

   END SUBROUTINE native_grid_extract_points

! **************************************************************************************************
!> \brief Allocate the compact gather buffers for the current field count.
!> \param nfg number of gathered fields (5 RKS, 10 LSD)
! **************************************************************************************************
   SUBROUTINE native_grid_alloc_fetch(nfg)
      INTEGER, INTENT(IN)                                :: nfg

      INTEGER                                            :: ncol, nreply

      skala_state%native_atom%native_need_nfg = nfg
      ncol = MAX(skala_state%native_atom%native_need_nneed, 1)
      IF (.NOT. ALLOCATED(skala_state%native_atom%native_need_values)) THEN
         ALLOCATE (skala_state%native_atom%native_need_values(nfg, ncol))
         skala_state%native_atom%native_need_values = 0.0_dp
      ELSE IF (SIZE(skala_state%native_atom%native_need_values, 1) /= nfg .OR. &
               SIZE(skala_state%native_atom%native_need_values, 2) /= ncol) THEN
         DEALLOCATE (skala_state%native_atom%native_need_values)
         ALLOCATE (skala_state%native_atom%native_need_values(nfg, ncol))
         skala_state%native_atom%native_need_values = 0.0_dp
      END IF
      nreply = MAX(skala_state%native_atom%native_need_nrecv*nfg, 1)
      IF (.NOT. ALLOCATED(skala_state%native_atom%native_need_reply)) THEN
         ALLOCATE (skala_state%native_atom%native_need_reply(nreply))
         skala_state%native_atom%native_need_reply = 0.0_dp
      ELSE IF (SIZE(skala_state%native_atom%native_need_reply) /= nreply) THEN
         DEALLOCATE (skala_state%native_atom%native_need_reply)
         ALLOCATE (skala_state%native_atom%native_need_reply(nreply))
         skala_state%native_atom%native_need_reply = 0.0_dp
      END IF

   END SUBROUTINE native_grid_alloc_fetch

! **************************************************************************************************
!> \brief Sparse all-to-all gather of the requested smooth-field values into the compact buffer.
!> \param group communicator over which the PW grid is distributed
! **************************************************************************************************
   SUBROUTINE native_grid_exchange_fetch(group)
      CLASS(mp_comm_type), INTENT(IN)                    :: group

      INTEGER                                            :: f, idx, k, nfg, npe, pos
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rcnt, rdisp, scnt, sdisp
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: recv

      nfg = skala_state%native_atom%native_need_nfg
      IF (skala_state%native_atom%native_need_replicated) THEN
         DO pos = 1, skala_state%native_atom%native_need_nneed
            DO f = 1, nfg
               skala_state%native_atom%native_need_values(f, pos) = skala_state%native_atom%native_need_reply((pos - 1)*nfg + f)
            END DO
         END DO
         RETURN
      END IF
      npe = group%num_pe
      ALLOCATE (scnt(npe), sdisp(npe), rcnt(npe), rdisp(npe))
      DO k = 1, npe
         scnt(k) = skala_state%native_atom%native_need_rcount(k)*nfg
         sdisp(k) = skala_state%native_atom%native_need_rdispl(k)*nfg
         rcnt(k) = skala_state%native_atom%native_need_scount(k)*nfg
         rdisp(k) = skala_state%native_atom%native_need_sdispl(k)*nfg
      END DO
      ALLOCATE (recv(MAX(SUM(rcnt), 1)))
      recv = 0.0_dp
      CALL group%alltoall(skala_state%native_atom%native_need_reply, scnt, sdisp, recv, rcnt, rdisp)
      DO k = 1, npe
         DO idx = 1, skala_state%native_atom%native_need_scount(k)
            pos = skala_state%native_atom%native_need_sdispl(k) + idx
            DO f = 1, nfg
        skala_state%native_atom%native_need_values(f, pos) = recv((skala_state%native_atom%native_need_sdispl(k) + idx - 1)*nfg + f)
            END DO
         END DO
      END DO
      DEALLOCATE (recv, scnt, sdisp, rcnt, rdisp)

   END SUBROUTINE native_grid_exchange_fetch

! **************************************************************************************************
!> \brief Sparse all-to-all scatter-add of the compact adjoint into local grid arrays.
!> \param group communicator over which the PW grid is distributed
!> \param adj compact adjoint in the same slot/owner order as the fetch
!> \param local_adj compact local adjoint, flattened in grid bounds order
! **************************************************************************************************
   SUBROUTINE native_grid_exchange_scatter(group, adj, local_adj)
      CLASS(mp_comm_type), INTENT(IN)                    :: group
      REAL(dp), DIMENSION(:, :), INTENT(IN), CONTIGUOUS  :: adj
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: local_adj

      INTEGER                                            :: f, idx, k, li, nfg, npe, nneed, pos
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rcnt, rdisp, scnt, sdisp
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: recv, send_flat

      nfg = skala_state%native_atom%native_need_nfg
      nneed = skala_state%native_atom%native_need_nneed
      local_adj = 0.0_dp
      IF (skala_state%native_atom%native_need_replicated) THEN
         DO pos = 1, nneed
            local_adj(:, skala_state%native_atom%native_need_local(pos)) = &
               local_adj(:, skala_state%native_atom%native_need_local(pos)) + adj(:, pos)
         END DO
         RETURN
      END IF
      npe = group%num_pe
      ALLOCATE (scnt(npe), sdisp(npe), rcnt(npe), rdisp(npe))
      DO k = 1, npe
         scnt(k) = skala_state%native_atom%native_need_scount(k)*nfg
         sdisp(k) = skala_state%native_atom%native_need_sdispl(k)*nfg
         rcnt(k) = skala_state%native_atom%native_need_rcount(k)*nfg
         rdisp(k) = skala_state%native_atom%native_need_rdispl(k)*nfg
      END DO
      ALLOCATE (send_flat(MAX(SUM(scnt), 1)))
      ALLOCATE (recv(MAX(SUM(rcnt), 1)))
      send_flat = 0.0_dp
      recv = 0.0_dp
      DO k = 1, npe
         DO idx = 1, skala_state%native_atom%native_need_scount(k)
            pos = skala_state%native_atom%native_need_sdispl(k) + idx
            DO f = 1, nfg
               send_flat((skala_state%native_atom%native_need_sdispl(k) + idx - 1)*nfg + f) = adj(f, pos)
            END DO
         END DO
      END DO
      CALL group%alltoall(send_flat, scnt, sdisp, recv, rcnt, rdisp)
      DO k = 1, npe
         DO idx = 1, skala_state%native_atom%native_need_rcount(k)
            li = skala_state%native_atom%native_need_ridx(skala_state%native_atom%native_need_rdispl(k) + idx)
            DO f = 1, nfg
               local_adj(f, li) = local_adj(f, li) + &
                                  recv((skala_state%native_atom%native_need_rdispl(k) + idx - 1)*nfg + f)
            END DO
         END DO
      END DO
      DEALLOCATE (recv, send_flat, scnt, sdisp, rcnt, rdisp)

   END SUBROUTINE native_grid_exchange_scatter

! **************************************************************************************************
!> \brief Interpolate the compact density/gradient/kinetic fields in one stencil traversal.
!> \param slot_rho density field slot
!> \param slot_gx x-gradient field slot
!> \param slot_gy y-gradient field slot
!> \param slot_gz z-gradient field slot
!> \param slot_kin kinetic-density field slot
!> \param stencil interpolation stencil
!> \param density_value ...
!> \param grad_value ...
!> \param kin_value ...
! **************************************************************************************************
   SUBROUTINE interpolate_native_grid_fields_compact(slot_rho, slot_gx, slot_gy, slot_gz, &
                                                     slot_kin, stencil, density_value, grad_value, &
                                                     kin_value)
      INTEGER, INTENT(IN)                                :: slot_rho, slot_gx, slot_gy, slot_gz, &
                                                            slot_kin
      TYPE(native_grid_interpolation_stencil_type), &
         INTENT(IN)                                      :: stencil
      REAL(dp), INTENT(OUT)                              :: density_value
      REAL(dp), DIMENSION(3), INTENT(OUT)                :: grad_value
      REAL(dp), INTENT(OUT)                              :: kin_value

      INTEGER                                            :: inode_x, inode_y, inode_z, pos
      INTEGER(int_8)                                     :: lin
      REAL(dp)                                           :: coefficient

      density_value = 0.0_dp
      grad_value = 0.0_dp
      kin_value = 0.0_dp
      IF (.NOT. stencil%active) RETURN
      DO inode_z = 1, native_grid_interp_npts
         IF (.NOT. stencil%valid(inode_z, 3)) CYCLE
         DO inode_y = 1, native_grid_interp_npts
            IF (.NOT. stencil%valid(inode_y, 2)) CYCLE
            DO inode_x = 1, native_grid_interp_npts
               IF (.NOT. stencil%valid(inode_x, 1)) CYCLE
               lin = native_need_linear_index(stencil, inode_x, inode_y, inode_z, skala_state%native_atom%native_need_npts)
               pos = native_need_lookup(lin)
               CPASSERT(pos > 0)
               coefficient = stencil%weight(inode_x, 1)* &
                             stencil%weight(inode_y, 2)* &
                             stencil%weight(inode_z, 3)
               density_value = density_value + coefficient*skala_state%native_atom%native_need_values(slot_rho, pos)
               grad_value(1) = grad_value(1) + coefficient*skala_state%native_atom%native_need_values(slot_gx, pos)
               grad_value(2) = grad_value(2) + coefficient*skala_state%native_atom%native_need_values(slot_gy, pos)
               grad_value(3) = grad_value(3) + coefficient*skala_state%native_atom%native_need_values(slot_gz, pos)
               kin_value = kin_value + coefficient*skala_state%native_atom%native_need_values(slot_kin, pos)
            END DO
         END DO
      END DO

   END SUBROUTINE interpolate_native_grid_fields_compact

! **************************************************************************************************
!> \brief Interpolate the compact density, gradient and kinetic-density fields for both spins of
!>        an LSD atom-composite point in a single stencil traversal.
!> \param stencil interpolation stencil
!> \param density_value per-spin density
!> \param grad_value per-spin density gradient
!> \param kin_value per-spin kinetic-density value
! **************************************************************************************************
   SUBROUTINE interpolate_native_grid_fields_compact_lsd(stencil, density_value, grad_value, kin_value)
      TYPE(native_grid_interpolation_stencil_type), &
         INTENT(IN)                                      :: stencil
      REAL(dp), DIMENSION(2), INTENT(OUT)                :: density_value
      REAL(dp), DIMENSION(3, 2), INTENT(OUT)             :: grad_value
      REAL(dp), DIMENSION(2), INTENT(OUT)                :: kin_value

      INTEGER                                            :: inode_x, inode_y, inode_z, pos
      INTEGER(int_8)                                     :: lin
      REAL(dp)                                           :: coefficient

      density_value = 0.0_dp
      grad_value = 0.0_dp
      kin_value = 0.0_dp
      IF (.NOT. stencil%active) RETURN
      DO inode_z = 1, native_grid_interp_npts
         IF (.NOT. stencil%valid(inode_z, 3)) CYCLE
         DO inode_y = 1, native_grid_interp_npts
            IF (.NOT. stencil%valid(inode_y, 2)) CYCLE
            DO inode_x = 1, native_grid_interp_npts
               IF (.NOT. stencil%valid(inode_x, 1)) CYCLE
               lin = native_need_linear_index(stencil, inode_x, inode_y, inode_z, skala_state%native_atom%native_need_npts)
               pos = native_need_lookup(lin)
               CPASSERT(pos > 0)
               coefficient = stencil%weight(inode_x, 1)* &
                             stencil%weight(inode_y, 2)* &
                             stencil%weight(inode_z, 3)
               density_value(1) = density_value(1) + &
                                  coefficient*skala_state%native_atom%native_need_values(1, pos)
               density_value(2) = density_value(2) + &
                                  coefficient*skala_state%native_atom%native_need_values(2, pos)
               kin_value(1) = kin_value(1) + &
                              coefficient*skala_state%native_atom%native_need_values(3, pos)
               kin_value(2) = kin_value(2) + &
                              coefficient*skala_state%native_atom%native_need_values(4, pos)
               grad_value(1, 1) = grad_value(1, 1) + &
                                  coefficient*skala_state%native_atom%native_need_values(5, pos)
               grad_value(2, 1) = grad_value(2, 1) + &
                                  coefficient*skala_state%native_atom%native_need_values(6, pos)
               grad_value(3, 1) = grad_value(3, 1) + &
                                  coefficient*skala_state%native_atom%native_need_values(7, pos)
               grad_value(1, 2) = grad_value(1, 2) + &
                                  coefficient*skala_state%native_atom%native_need_values(8, pos)
               grad_value(2, 2) = grad_value(2, 2) + &
                                  coefficient*skala_state%native_atom%native_need_values(9, pos)
               grad_value(3, 2) = grad_value(3, 2) + &
                                  coefficient*skala_state%native_atom%native_need_values(10, pos)
            END DO
         END DO
      END DO
   END SUBROUTINE interpolate_native_grid_fields_compact_lsd

! **************************************************************************************************
!> \brief Interpolate several compact field values from one precomputed stencil, traversing the
!>        stencil once for all requested slots.
!> \param slots field slots
!> \param nslots number of requested slots
!> \param stencil precomputed interpolation stencil
!> \param values interpolated values, one per slot
! **************************************************************************************************
   SUBROUTINE interpolate_native_grid_values_compact(slots, nslots, stencil, values)
      INTEGER, INTENT(IN)                                :: nslots
      INTEGER, DIMENSION(nslots), INTENT(IN)             :: slots
      TYPE(native_grid_interpolation_stencil_type), &
         INTENT(IN)                                      :: stencil
      REAL(dp), DIMENSION(nslots), INTENT(OUT)           :: values

      INTEGER                                            :: inode_x, inode_y, inode_z, islot, pos
      INTEGER(int_8)                                     :: lin
      REAL(dp)                                           :: coefficient

      values = 0.0_dp
      IF (.NOT. stencil%active) RETURN
      DO inode_z = 1, native_grid_interp_npts
         IF (.NOT. stencil%valid(inode_z, 3)) CYCLE
         DO inode_y = 1, native_grid_interp_npts
            IF (.NOT. stencil%valid(inode_y, 2)) CYCLE
            DO inode_x = 1, native_grid_interp_npts
               IF (.NOT. stencil%valid(inode_x, 1)) CYCLE
               lin = native_need_linear_index(stencil, inode_x, inode_y, inode_z, skala_state%native_atom%native_need_npts)
               pos = native_need_lookup(lin)
               CPASSERT(pos > 0)
               coefficient = stencil%weight(inode_x, 1)* &
                             stencil%weight(inode_y, 2)* &
                             stencil%weight(inode_z, 3)
               DO islot = 1, nslots
                  values(islot) = values(islot) + &
                                  coefficient*skala_state%native_atom%native_need_values(slots(islot), pos)
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE interpolate_native_grid_values_compact

! **************************************************************************************************
!> \brief Interpolate the Cartesian gradients of several compact fields from one precomputed
!>        stencil, traversing the stencil once for all requested slots.
!> \param slots field slots
!> \param nslots number of requested slots
!> \param pw_grid auxiliary grid (for the incremental cell matrix)
!> \param stencil precomputed interpolation stencil
!> \param gradients interpolated gradients, gradients(:, islot) per slot
! **************************************************************************************************
   SUBROUTINE interpolate_native_grid_gradients_compact(slots, nslots, pw_grid, stencil, gradients)
      INTEGER, INTENT(IN)                                :: nslots
      INTEGER, DIMENSION(nslots), INTENT(IN)             :: slots
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(native_grid_interpolation_stencil_type), &
         INTENT(IN)                                      :: stencil
      REAL(dp), DIMENSION(3, nslots), INTENT(OUT)        :: gradients

      INTEGER                                            :: idir, inode_x, inode_y, inode_z, islot, &
                                                            pos
      INTEGER(int_8)                                     :: lin
      REAL(dp)                                           :: value
      REAL(dp), DIMENSION(3)                             :: coefficient
      REAL(dp), DIMENSION(3, nslots)                     :: gradient_relative

      CPASSERT(ASSOCIATED(pw_grid))
      gradients = 0.0_dp
      IF (.NOT. stencil%active) RETURN
      gradient_relative = 0.0_dp
      DO inode_z = 1, native_grid_interp_npts
         IF (.NOT. stencil%valid(inode_z, 3)) CYCLE
         DO inode_y = 1, native_grid_interp_npts
            IF (.NOT. stencil%valid(inode_y, 2)) CYCLE
            DO inode_x = 1, native_grid_interp_npts
               IF (.NOT. stencil%valid(inode_x, 1)) CYCLE
               lin = native_need_linear_index(stencil, inode_x, inode_y, inode_z, skala_state%native_atom%native_need_npts)
               pos = native_need_lookup(lin)
               CPASSERT(pos > 0)
               coefficient(1) = stencil%derivative_weight(inode_x, 1)* &
                                stencil%weight(inode_y, 2)*stencil%weight(inode_z, 3)
               coefficient(2) = stencil%weight(inode_x, 1)* &
                                stencil%derivative_weight(inode_y, 2)*stencil%weight(inode_z, 3)
               coefficient(3) = stencil%weight(inode_x, 1)* &
                                stencil%weight(inode_y, 2)*stencil%derivative_weight(inode_z, 3)
               DO islot = 1, nslots
                  value = skala_state%native_atom%native_need_values(slots(islot), pos)
                  DO idir = 1, 3
                     gradient_relative(idir, islot) = gradient_relative(idir, islot) + &
                                                      coefficient(idir)*value
                  END DO
               END DO
            END DO
         END DO
      END DO
      gradients = MATMUL(TRANSPOSE(pw_grid%dh_inv), gradient_relative)

   END SUBROUTINE interpolate_native_grid_gradients_compact

! **************************************************************************************************
!> \brief Build the tensor-product interpolation stencil for one Cartesian point.
!> \param stencil ...
!> \param pw_grid ...
!> \param cell ...
!> \param point ...
!> \param wrap_auxiliary_cell wrap all auxiliary-grid directions
!> \param indices_only skip interpolation weights when only grid indices are needed
! **************************************************************************************************
   SUBROUTINE create_native_grid_interpolation_stencil(stencil, pw_grid, cell, point, &
                                                       wrap_auxiliary_cell, indices_only)
      TYPE(native_grid_interpolation_stencil_type), &
         INTENT(OUT)                                     :: stencil
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: point
      LOGICAL, INTENT(IN), OPTIONAL                      :: wrap_auxiliary_cell, indices_only

      INTEGER                                            :: idir, inode, relative_index
      INTEGER, DIMENSION(3)                              :: base
      LOGICAL                                            :: build_weights, wrap_grid
      REAL(KIND=dp), DIMENSION(3)                        :: fraction, relative

      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(ASSOCIATED(cell))
      stencil%active = .FALSE.
      stencil%valid = .FALSE.
      stencil%weight = 0.0_dp
      stencil%derivative_weight = 0.0_dp
      build_weights = .TRUE.
      IF (PRESENT(indices_only)) build_weights = .NOT. indices_only
      wrap_grid = ANY(cell%perd == 0)
      IF (PRESENT(wrap_auxiliary_cell)) wrap_grid = wrap_auxiliary_cell
      relative = MATMUL(pw_grid%dh_inv, point)
      DO idir = 1, 3
         IF (cell%perd(idir) == 1 .OR. wrap_grid) THEN
            relative(idir) = MODULO(relative(idir), REAL(pw_grid%npts(idir), KIND=dp))
         ELSE IF (relative(idir) <= -REAL(native_grid_interp_offset_max, dp) .OR. &
                  relative(idir) >= REAL(pw_grid%npts(idir) - &
                                         native_grid_interp_offset_min, KIND=dp)) THEN
            RETURN
         END IF
         base(idir) = FLOOR(relative(idir))
         FRACTION(idir) = relative(idir) - REAL(base(idir), KIND=dp)
         IF (build_weights) THEN
            CALL native_grid_lagrange_weights(FRACTION(idir), stencil%weight(:, idir), &
                                              stencil%derivative_weight(:, idir))
         END IF
         DO inode = 1, native_grid_interp_npts
            relative_index = base(idir) + native_grid_interp_offset_min + inode - 1
            IF (cell%perd(idir) == 1 .OR. wrap_grid) THEN
               relative_index = MODULO(relative_index, pw_grid%npts(idir))
               stencil%valid(inode, idir) = .TRUE.
            ELSE IF (relative_index >= 0 .AND. relative_index < pw_grid%npts(idir)) THEN
               stencil%valid(inode, idir) = .TRUE.
            END IF
            stencil%relative_index(inode, idir) = relative_index
         END DO
      END DO
      stencil%active = .TRUE.

   END SUBROUTINE create_native_grid_interpolation_stencil

! **************************************************************************************************
!> \brief Build the value and optional derivative weights for the native-grid interpolation.
!> \param fraction fractional coordinate between two grid points
!> \param weights interpolation weights
!> \param derivative_weights optional derivatives with respect to fraction
! **************************************************************************************************
   SUBROUTINE native_grid_lagrange_weights(fraction, weights, derivative_weights)
      REAL(KIND=dp), INTENT(IN)                          :: fraction
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: weights
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT), OPTIONAL :: derivative_weights

      INTEGER                                            :: inode, jnode, knode, npoints
      REAL(KIND=dp)                                      :: denominator, derivative, numerator
      REAL(KIND=dp), DIMENSION(SIZE(weights))            :: nodes

      npoints = SIZE(weights)
      IF (PRESENT(derivative_weights)) THEN
         CPASSERT(SIZE(derivative_weights) == npoints)
      END IF
      DO inode = 1, npoints
         nodes(inode) = REAL(native_grid_interp_offset_min + inode - 1, dp)
      END DO
      DO inode = 1, npoints
         denominator = 1.0_dp
         numerator = 1.0_dp
         DO jnode = 1, npoints
            IF (jnode == inode) CYCLE
            denominator = denominator*(nodes(inode) - nodes(jnode))
            numerator = numerator*(fraction - nodes(jnode))
         END DO
         weights(inode) = numerator/denominator
         IF (PRESENT(derivative_weights)) THEN
            derivative = 0.0_dp
            DO knode = 1, npoints
               IF (knode == inode) CYCLE
               numerator = 1.0_dp
               DO jnode = 1, npoints
                  IF (jnode == inode .OR. jnode == knode) CYCLE
                  numerator = numerator*(fraction - nodes(jnode))
               END DO
               derivative = derivative + numerator/denominator
            END DO
            derivative_weights(inode) = derivative
         END IF
      END DO

   END SUBROUTINE native_grid_lagrange_weights

! **************************************************************************************************
!> \brief Conservatively bound periodic images that can intersect a Cartesian support sphere.
!> \param cell current cell and periodicity
!> \param fractional target-minus-source displacement in lattice coordinates
!> \param cutoff current, density-dependent Cartesian support radius
!> \param image_shell original enclosing image shell
!> \param lower inclusive lower image indices
!> \param upper inclusive upper image indices (may be smaller than lower)
! **************************************************************************************************
   SUBROUTINE atom_grid_image_bounds(cell, fractional, cutoff, image_shell, lower, upper)
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(dp), INTENT(IN)                               :: fractional(3), cutoff
      INTEGER, INTENT(IN)                                :: image_shell(3)
      INTEGER, INTENT(OUT)                               :: lower(3), upper(3)

      INTEGER                                            :: center, idir
      REAL(dp)                                           :: condition_bound, padding, radius

      ! |(H^-1 d)_i| <= cutoff ||(H^-1)_i|| also holds for triclinic cells.
      ! Inflate the bounds for cancellation in fractional coordinates and the inverse cell.
      condition_bound = MAX(1.0_dp, SUM(ABS(cell%hmat))*SUM(ABS(cell%h_inv)))
      lower = 0
      upper = 0
      DO idir = 1, 3
         IF (cell%perd(idir) == 0) CYCLE
         center = NINT(fractional(idir))
         radius = cutoff*SQRT(SUM(cell%h_inv(idir, :)**2))
         padding = 128.0_dp*EPSILON(1.0_dp)*condition_bound* &
                   MAX(1.0_dp, MAXVAL(ABS(fractional)), radius)
         lower(idir) = MAX(center - image_shell(idir), CEILING(fractional(idir) - radius - padding))
         upper(idir) = MIN(center + image_shell(idir), FLOOR(fractional(idir) + radius + padding))
      END DO
   END SUBROUTINE atom_grid_image_bounds

! **************************************************************************************************
!> \brief Express one radial node's first and second derivatives as linear combinations of the
!>        local node values used by the C2-continuous quintic Hermite interpolant.
!> \param grid_atom radial source grid
!> \param descending whether radial nodes are stored in descending order
!> \param node logical index of the node whose derivatives are required
!> \param logical_start first logical index represented by the coefficient arrays
!> \param slope coefficients of the first derivative
!> \param curvature coefficients of the second derivative
! **************************************************************************************************
   SUBROUTINE radial_node_derivative_coefficients( &
      grid_atom, descending, node, logical_start, slope, curvature)
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      LOGICAL, INTENT(IN)                                :: descending
      INTEGER, INTENT(IN)                                :: node, logical_start
      REAL(dp), DIMENSION(4), INTENT(OUT)                :: slope, curvature

      INTEGER                                            :: center, hi, lo, n
      REAL(dp)                                           :: h_hi, h_lo, x_center, x_hi, x_lo

      n = grid_atom%nr
      slope = 0.0_dp
      curvature = 0.0_dp
      IF (node == 1) THEN
         lo = 1
         hi = 2
      ELSE IF (node == n) THEN
         lo = n - 1
         hi = n
      ELSE
         lo = node - 1
         hi = node + 1
      END IF
      x_lo = grid_atom%rad(MERGE(n + 1 - lo, lo, descending))
      x_hi = grid_atom%rad(MERGE(n + 1 - hi, hi, descending))
      slope(lo - logical_start + 1) = -1.0_dp/(x_hi - x_lo)
      slope(hi - logical_start + 1) = 1.0_dp/(x_hi - x_lo)

      IF (n < 3) RETURN
      center = MIN(MAX(node, 2), n - 1)
      lo = center - 1
      hi = center + 1
      x_lo = grid_atom%rad(MERGE(n + 1 - lo, lo, descending))
      x_center = grid_atom%rad(MERGE(n + 1 - center, center, descending))
      x_hi = grid_atom%rad(MERGE(n + 1 - hi, hi, descending))
      h_lo = x_center - x_lo
      h_hi = x_hi - x_center
      curvature(lo - logical_start + 1) = 2.0_dp/(h_lo*(h_lo + h_hi))
      curvature(center - logical_start + 1) = &
         -2.0_dp*(1.0_dp/h_lo + 1.0_dp/h_hi)/(h_lo + h_hi)
      curvature(hi - logical_start + 1) = 2.0_dp/(h_hi*(h_lo + h_hi))

   END SUBROUTINE radial_node_derivative_coefficients

! **************************************************************************************************
!> \brief Build mutually consistent value and Cartesian-derivative weights for interpolation from
!>        a GAPW radial/Lebedev atom grid to one arbitrary local displacement.
!> \param grid_atom source atom grid
!> \param harmonics source spherical harmonics
!> \param displacement point minus source-image coordinate
!> \param cutoff compact support radius
!> \param radial_indices active radial nodes
!> \param radial_weights interpolation weights on the active radial nodes
!> \param radial_derivative_weights optional radial derivatives of the interpolation weights
!> \param nradial number of active radial nodes
!> \param angular_weights angular interpolation weights
!> \param angular_derivative_weights optional Cartesian derivatives of the angular weights
!> \param active whether the point lies inside the compact support
!> \param direction_out ...
! **************************************************************************************************
   SUBROUTINE atom_grid_interpolation_weights( &
      grid_atom, harmonics, displacement, cutoff, radial_indices, radial_weights, &
      radial_derivative_weights, nradial, angular_weights, angular_derivative_weights, active, &
      direction_out)
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: displacement
      REAL(dp), INTENT(IN)                               :: cutoff
      INTEGER, DIMENSION(4), INTENT(OUT)                 :: radial_indices
      REAL(dp), DIMENSION(4), INTENT(OUT)                :: radial_weights
      REAL(dp), DIMENSION(4), INTENT(OUT), OPTIONAL      :: radial_derivative_weights
      INTEGER, INTENT(OUT)                               :: nradial
      REAL(dp), DIMENSION(:), INTENT(OUT)                :: angular_weights
      REAL(dp), DIMENSION(:, :), INTENT(OUT), OPTIONAL   :: angular_derivative_weights
      LOGICAL, INTENT(OUT)                               :: active
      REAL(dp), DIMENSION(3), INTENT(OUT), OPTIONAL      :: direction_out

      INTEGER                                            :: ia, ic, inode, iso, l, left, left_pos, &
                                                            logical_end, logical_start, lx, ly, &
                                                            lz, n, right_pos, shell_index
      LOGICAL                                            :: descending
      REAL(dp)                                           :: dh00, dh01, dh10, dh11, dh20, dh21, h00, &
                                                            h01, h10, h11, h20, h21, h_interval, &
                                                            monomial, radius, solid_derivative, t, &
                                                            x1, x2
      REAL(dp), DIMENSION(3)                             :: direction
      REAL(dp), DIMENSION(harmonics%max_s_harm)          :: angular_values
      REAL(dp), DIMENSION(4)                             :: curvature_left, curvature_right, &
                                                            slope_left, slope_right
      REAL(dp), DIMENSION(3, harmonics%max_s_harm)       :: angular_value_derivatives

      CPASSERT(ASSOCIATED(grid_atom))
      CPASSERT(ASSOCIATED(harmonics))
      CPASSERT(SIZE(angular_weights) == grid_atom%ng_sphere)
      IF (PRESENT(angular_derivative_weights)) THEN
         CPASSERT(SIZE(angular_derivative_weights, 1) == 3)
         CPASSERT(SIZE(angular_derivative_weights, 2) == grid_atom%ng_sphere)
      END IF

      radial_indices = 0
      radial_weights = 0.0_dp
      IF (PRESENT(radial_derivative_weights)) radial_derivative_weights = 0.0_dp
      angular_weights = 0.0_dp
      IF (PRESENT(angular_derivative_weights)) angular_derivative_weights = 0.0_dp
      IF (PRESENT(direction_out)) direction_out = 0.0_dp
      nradial = 0
      active = .FALSE.
      n = grid_atom%nr
      descending = grid_atom%rad(1) > grid_atom%rad(n)
      radius = SQRT(SUM(displacement**2))
      IF (radius > cutoff .OR. &
          radius > grid_atom%rad(MERGE(1, n, descending))) RETURN
      IF (radius <= 1.0E-12_dp) RETURN
      direction = displacement/radius
      IF (PRESENT(direction_out)) direction_out = direction

      CPASSERT(n >= 2)
      left = n - 1
      IF (radius <= grid_atom%rad(MERGE(n, 1, descending))) THEN
         left = 1
      ELSE
         DO inode = 1, n - 1
            IF (radius <= grid_atom%rad(MERGE(n - inode, inode + 1, descending))) THEN
               left = inode
               EXIT
            END IF
         END DO
      END IF

      logical_start = MAX(1, left - 1)
      logical_end = MIN(n, left + 2)
      nradial = logical_end - logical_start + 1
      DO inode = 1, nradial
         radial_indices(inode) = MERGE(n + 2 - logical_start - inode, &
                                       logical_start + inode - 1, descending)
      END DO
      left_pos = left - logical_start + 1
      right_pos = left_pos + 1
      x1 = grid_atom%rad(radial_indices(left_pos))
      x2 = grid_atom%rad(radial_indices(right_pos))
      h_interval = x2 - x1
      t = (radius - x1)/h_interval

      h00 = 1.0_dp - 10.0_dp*t**3 + 15.0_dp*t**4 - 6.0_dp*t**5
      h10 = t - 6.0_dp*t**3 + 8.0_dp*t**4 - 3.0_dp*t**5
      h20 = 0.5_dp*(t**2 - 3.0_dp*t**3 + 3.0_dp*t**4 - t**5)
      h01 = 10.0_dp*t**3 - 15.0_dp*t**4 + 6.0_dp*t**5
      h11 = -4.0_dp*t**3 + 7.0_dp*t**4 - 3.0_dp*t**5
      h21 = 0.5_dp*(t**3 - 2.0_dp*t**4 + t**5)

      CALL radial_node_derivative_coefficients( &
         grid_atom, descending, left, logical_start, slope_left, curvature_left)
      CALL radial_node_derivative_coefficients( &
         grid_atom, descending, left + 1, logical_start, slope_right, curvature_right)
      radial_weights(left_pos) = radial_weights(left_pos) + h00
      radial_weights(right_pos) = radial_weights(right_pos) + h01
      radial_weights(1:nradial) = radial_weights(1:nradial) + &
                                  h_interval*(h10*slope_left(1:nradial) + &
                                              h11*slope_right(1:nradial)) + &
                                  h_interval**2*(h20*curvature_left(1:nradial) + &
                                                 h21*curvature_right(1:nradial))
      IF (PRESENT(radial_derivative_weights)) THEN
         dh00 = -30.0_dp*t**2 + 60.0_dp*t**3 - 30.0_dp*t**4
         dh10 = 1.0_dp - 18.0_dp*t**2 + 32.0_dp*t**3 - 15.0_dp*t**4
         dh20 = 0.5_dp*(2.0_dp*t - 9.0_dp*t**2 + 12.0_dp*t**3 - 5.0_dp*t**4)
         dh01 = 30.0_dp*t**2 - 60.0_dp*t**3 + 30.0_dp*t**4
         dh11 = -12.0_dp*t**2 + 28.0_dp*t**3 - 15.0_dp*t**4
         dh21 = 0.5_dp*(3.0_dp*t**2 - 8.0_dp*t**3 + 5.0_dp*t**4)
         radial_derivative_weights(left_pos) = &
            radial_derivative_weights(left_pos) + dh00/h_interval
         radial_derivative_weights(right_pos) = &
            radial_derivative_weights(right_pos) + dh01/h_interval
         radial_derivative_weights(1:nradial) = radial_derivative_weights(1:nradial) + &
                                                dh10*slope_left(1:nradial) + &
                                                dh11*slope_right(1:nradial) + &
                                                h_interval*(dh20*curvature_left(1:nradial) + &
                                                            dh21*curvature_right(1:nradial))
      END IF

      angular_values = 0.0_dp
      IF (PRESENT(angular_derivative_weights)) angular_value_derivatives = 0.0_dp
      DO iso = 1, harmonics%max_s_harm
         l = indso(1, iso)
         CALL y_lm(direction, angular_values(iso), l, indso(2, iso))
         IF (l == 0 .OR. .NOT. PRESENT(angular_derivative_weights)) CYCLE
         shell_index = iso - nsoset(l - 1)
         DO ic = 1, nco(l)
            lx = indco(1, ic + ncoset(l - 1))
            ly = indco(2, ic + ncoset(l - 1))
            lz = indco(3, ic + ncoset(l - 1))
            IF (lx > 0) THEN
               monomial = REAL(lx, dp)*direction(1)**(lx - 1)* &
                          direction(2)**ly*direction(3)**lz
               angular_value_derivatives(1, iso) = angular_value_derivatives(1, iso) + &
                                                   orbtramat(l)%slm(shell_index, ic)*monomial
            END IF
            IF (ly > 0) THEN
               monomial = direction(1)**lx*REAL(ly, dp)*direction(2)**(ly - 1)* &
                          direction(3)**lz
               angular_value_derivatives(2, iso) = angular_value_derivatives(2, iso) + &
                                                   orbtramat(l)%slm(shell_index, ic)*monomial
            END IF
            IF (lz > 0) THEN
               monomial = direction(1)**lx*direction(2)**ly* &
                          REAL(lz, dp)*direction(3)**(lz - 1)
               angular_value_derivatives(3, iso) = angular_value_derivatives(3, iso) + &
                                                   orbtramat(l)%slm(shell_index, ic)*monomial
            END IF
         END DO
         DO inode = 1, 3
            solid_derivative = angular_value_derivatives(inode, iso)
            angular_value_derivatives(inode, iso) = (solid_derivative - &
                                                     REAL(l, dp)*angular_values(iso)* &
                                                     direction(inode))/radius
         END DO
      END DO

      DO ia = 1, grid_atom%ng_sphere
         angular_weights(ia) = grid_atom%wa(ia)* &
                               DOT_PRODUCT(harmonics%slm(ia, 1:harmonics%max_s_harm), &
                                           angular_values)
         IF (PRESENT(angular_derivative_weights)) THEN
            DO inode = 1, 3
               angular_derivative_weights(inode, ia) = grid_atom%wa(ia)* &
                                                       DOT_PRODUCT(harmonics%slm(ia, 1:harmonics%max_s_harm), &
                                                                   angular_value_derivatives(inode, :))
            END DO
         END IF
      END DO
      active = .TRUE.

   END SUBROUTINE atom_grid_interpolation_weights

! **************************************************************************************************
!> \brief Interpolate hard-minus-soft primitive fields and their spatial derivatives from one
!>        GAPW atom grid.
!> \param grid_atom radial and angular source grid
!> \param harmonics spherical-harmonic representation of the source grid
!> \param displacement target point relative to the source-atom image
!> \param cutoff compact support radius of the source fields
!> \param nspins number of spin channels
!> \param rho_h hard one-center density values
!> \param rho_s soft one-center density values
!> \param drho_h hard one-center density-gradient values
!> \param drho_s soft one-center density-gradient values
!> \param tau_h hard one-center kinetic-energy-density values
!> \param tau_s soft one-center kinetic-energy-density values
!> \param density interpolated hard-minus-soft density
!> \param gradient interpolated hard-minus-soft density gradient
!> \param kin interpolated hard-minus-soft kinetic-energy density
!> \param density_spatial Cartesian derivatives of density
!> \param gradient_spatial Cartesian derivatives of the density gradient
!> \param kin_spatial Cartesian derivatives of the kinetic-energy density
!> \param calculate_spatial evaluate spatial derivatives (default true), otherwise return zeros
! **************************************************************************************************
   SUBROUTINE interpolate_gapw_atom_grid_fields( &
      grid_atom, harmonics, displacement, cutoff, nspins, rho_h, rho_s, drho_h, drho_s, &
      tau_h, tau_s, density, gradient, kin, density_spatial, gradient_spatial, kin_spatial, calculate_spatial)
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: displacement
      REAL(dp), INTENT(IN)                               :: cutoff
      INTEGER, INTENT(IN)                                :: nspins
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: rho_h, rho_s
      REAL(dp), DIMENSION(:, :, :, :), INTENT(IN)        :: drho_h, drho_s
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: tau_h, tau_s
      REAL(dp), DIMENSION(2), INTENT(OUT)                :: density
      REAL(dp), DIMENSION(3, 2), INTENT(OUT)             :: gradient
      REAL(dp), DIMENSION(2), INTENT(OUT)                :: kin
      REAL(dp), DIMENSION(3, 2), INTENT(OUT)             :: density_spatial
      REAL(dp), DIMENSION(3, 3, 2), INTENT(OUT)          :: gradient_spatial
      REAL(dp), DIMENSION(3, 2), INTENT(OUT)             :: kin_spatial
      LOGICAL, INTENT(IN), OPTIONAL                      :: calculate_spatial

      INTEGER                                            :: ia, idir, inode, ir, ispin, jdir, nradial
      INTEGER, DIMENSION(4)                              :: radial_indices
      LOGICAL                                            :: active, need_spatial
      REAL(dp)                                           :: derivative_weight, weight
      REAL(dp), DIMENSION(3)                             :: direction
      REAL(dp), DIMENSION(grid_atom%ng_sphere)           :: angular_weights
      REAL(dp), DIMENSION(4)                             :: radial_derivative_weights, radial_weights
      REAL(dp), DIMENSION(3, grid_atom%ng_sphere)        :: angular_derivative_weights

      density = 0.0_dp
      gradient = 0.0_dp
      kin = 0.0_dp
      density_spatial = 0.0_dp
      gradient_spatial = 0.0_dp
      kin_spatial = 0.0_dp
      need_spatial = .TRUE.
      IF (PRESENT(calculate_spatial)) need_spatial = calculate_spatial
      IF (need_spatial) THEN
         CALL atom_grid_interpolation_weights( &
            grid_atom, harmonics, displacement, cutoff, radial_indices, radial_weights, &
            radial_derivative_weights, nradial, angular_weights, angular_derivative_weights, active, &
            direction_out=direction)
      ELSE
         CALL atom_grid_interpolation_weights( &
            grid_atom, harmonics, displacement, cutoff, radial_indices, radial_weights, &
            nradial=nradial, angular_weights=angular_weights, active=active, direction_out=direction)
      END IF
      IF (active) THEN
         DO inode = 1, nradial
            ir = radial_indices(inode)
            DO ia = 1, grid_atom%ng_sphere
               weight = radial_weights(inode)*angular_weights(ia)
               DO ispin = 1, nspins
                  density(ispin) = density(ispin) + &
                                   weight*(rho_h(ia, ir, ispin) - rho_s(ia, ir, ispin))
                  kin(ispin) = kin(ispin) + &
                               weight*(tau_h(ia, ir, ispin) - tau_s(ia, ir, ispin))
                  DO idir = 1, 3
                     gradient(idir, ispin) = gradient(idir, ispin) + &
                                             weight*(drho_h(idir, ia, ir, ispin) - &
                                                     drho_s(idir, ia, ir, ispin))
                     IF (.NOT. need_spatial) CYCLE
                     derivative_weight = radial_derivative_weights(inode)* &
                                         direction(idir)*angular_weights(ia) + &
                                         radial_weights(inode)*angular_derivative_weights(idir, ia)
                     density_spatial(idir, ispin) = density_spatial(idir, ispin) + &
                                                    derivative_weight*(rho_h(ia, ir, ispin) - rho_s(ia, ir, ispin))
                     kin_spatial(idir, ispin) = kin_spatial(idir, ispin) + &
                                                derivative_weight*(tau_h(ia, ir, ispin) - tau_s(ia, ir, ispin))
                     DO jdir = 1, 3
                        gradient_spatial(jdir, idir, ispin) = &
                           gradient_spatial(jdir, idir, ispin) + derivative_weight*( &
                           drho_h(jdir, ia, ir, ispin) - drho_s(jdir, ia, ir, ispin))
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END IF
   END SUBROUTINE interpolate_gapw_atom_grid_fields

! **************************************************************************************************
!> \brief Apply the exact transpose of interpolate_gapw_atom_grid_fields to one-center potentials.
!> \param grid_atom radial and angular source grid
!> \param harmonics spherical-harmonic representation of the source grid
!> \param displacement target point relative to the source-atom image
!> \param cutoff compact support radius of the source fields
!> \param nspins number of spin channels
!> \param density_adjoint model derivative with respect to density
!> \param gradient_adjoint model derivative with respect to the density gradient
!> \param kin_adjoint model derivative with respect to kinetic-energy density
!> \param vxc_h accumulated hard one-center density potential
!> \param vxc_s accumulated soft one-center density potential
!> \param vxg_h accumulated hard one-center density-gradient potential
!> \param vxg_s accumulated soft one-center density-gradient potential
!> \param vtau_h accumulated hard one-center kinetic-energy-density potential
!> \param vtau_s accumulated soft one-center kinetic-energy-density potential
! **************************************************************************************************
   SUBROUTINE add_gapw_atom_grid_interpolation_adjoint( &
      grid_atom, harmonics, displacement, cutoff, nspins, density_adjoint, gradient_adjoint, &
      kin_adjoint, vxc_h, vxc_s, vxg_h, vxg_s, vtau_h, vtau_s)
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: displacement
      REAL(dp), INTENT(IN)                               :: cutoff
      INTEGER, INTENT(IN)                                :: nspins
      REAL(dp), DIMENSION(2), INTENT(IN)                 :: density_adjoint
      REAL(dp), DIMENSION(3, 2), INTENT(IN)              :: gradient_adjoint
      REAL(dp), DIMENSION(2), INTENT(IN)                 :: kin_adjoint
      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: vxc_h, vxc_s
      REAL(dp), DIMENSION(:, :, :, :), INTENT(INOUT)     :: vxg_h, vxg_s
      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: vtau_h, vtau_s

      INTEGER                                            :: ia, idir, inode, ir, ispin, nradial
      INTEGER, DIMENSION(4)                              :: radial_indices
      LOGICAL                                            :: active
      REAL(dp)                                           :: value
      REAL(dp), DIMENSION(4)                             :: radial_weights
      REAL(dp), DIMENSION(grid_atom%ng_sphere)           :: angular_weights

      CALL atom_grid_interpolation_weights( &
         grid_atom, harmonics, displacement, cutoff, radial_indices, radial_weights, &
         nradial=nradial, angular_weights=angular_weights, active=active)
      IF (active) THEN
         DO inode = 1, nradial
            ir = radial_indices(inode)
            DO ia = 1, grid_atom%ng_sphere
               value = radial_weights(inode)*angular_weights(ia)
               DO ispin = 1, nspins
                  ! CP2K applies the hard-minus-soft sign when the two one-center
                  ! matrices are assembled, so both stored potentials carry the
                  ! same transpose-interpolation coefficient.
                  vxc_h(ia, ir, ispin) = vxc_h(ia, ir, ispin) + value*density_adjoint(ispin)
                  vxc_s(ia, ir, ispin) = vxc_s(ia, ir, ispin) + value*density_adjoint(ispin)
                  vtau_h(ia, ir, ispin) = vtau_h(ia, ir, ispin) + value*kin_adjoint(ispin)
                  vtau_s(ia, ir, ispin) = vtau_s(ia, ir, ispin) + value*kin_adjoint(ispin)
                  DO idir = 1, 3
                     vxg_h(idir, ia, ir, ispin) = vxg_h(idir, ia, ir, ispin) + &
                                                  value*gradient_adjoint(idir, ispin)
                     vxg_s(idir, ia, ir, ispin) = vxg_s(idir, ia, ir, ispin) + &
                                                  value*gradient_adjoint(idir, ispin)
                  END DO
               END DO
            END DO
         END DO
      END IF
   END SUBROUTINE add_gapw_atom_grid_interpolation_adjoint

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param exc1 the on-body ex energy contribution
!> \param gradient_atom_set ...
! **************************************************************************************************
   SUBROUTINE calculate_vxc_atom_epr(qs_env, exc1, gradient_atom_set)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), INTENT(INOUT)                            :: exc1
      TYPE(nablavks_atom_type), DIMENSION(:), POINTER    :: gradient_atom_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_vxc_atom_epr'

      INTEGER                                            :: bo(2), handle, ia, iat, iatom, idir, &
                                                            ikind, ir, ispin, myfun, na, natom, &
                                                            nr, nspins, num_pe, on
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: accint, donlcc, gradient_f, lsd, nlcc, &
                                                            paw_atom, tau_f
      REAL(dp)                                           :: agr, alpha, density_cut, exc_h, exc_s, &
                                                            gradient_cut, tau_cut
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: fw
      REAL(dp), DIMENSION(1, 1, 1)                       :: tau_d
      REAL(dp), DIMENSION(1, 1, 1, 1)                    :: rho_d
      REAL(dp), DIMENSION(:, :), POINTER                 :: rho_nlcc, weight_h, weight_s
      REAL(dp), DIMENSION(:, :, :), POINTER              :: rho_h, rho_s, tau_h, tau_s, vtau_h, &
                                                            vtau_s, vxc_h, vxc_s
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: drho_h, drho_s, vxg_h, vxg_s
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: my_kind_set
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: dr_h, dr_s, int_hh, int_ss, r_h, r_s
      TYPE(rho_atom_coeff), DIMENSION(:, :), POINTER     :: r_h_d, r_s_d
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: my_rho_atom_set
      TYPE(rho_atom_type), POINTER                       :: rho_atom
      TYPE(section_vals_type), POINTER                   :: input, my_xc_section, xc_fun_section
      TYPE(tau_basis_cache_type)                         :: tau_basis_cache
      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set_h, rho_set_s

! -------------------------------------------------------------------------

      CALL timeset(routineN, handle)

      NULLIFY (atom_list)
      NULLIFY (my_kind_set)
      NULLIFY (atomic_kind_set)
      NULLIFY (grid_atom)
      NULLIFY (harmonics)
      NULLIFY (input)
      NULLIFY (para_env)
      NULLIFY (rho_atom)
      NULLIFY (my_rho_atom_set)
      NULLIFY (rho_nlcc)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=my_kind_set, &
                      input=input, &
                      rho_atom_set=my_rho_atom_set)

      nlcc = has_nlcc(my_kind_set)
      accint = dft_control%qs_control%gapw_control%accurate_xcint

      my_xc_section => section_vals_get_subs_vals(input, &
                                                  "PROPERTIES%LINRES%EPR%PRINT%G_TENSOR%XC")
      xc_fun_section => section_vals_get_subs_vals(my_xc_section, "XC_FUNCTIONAL")
      CALL section_vals_val_get(xc_fun_section, "_SECTION_PARAMETERS_", &
                                i_val=myfun)

      IF (myfun == xc_none) THEN
         exc1 = 0.0_dp
         my_rho_atom_set(:)%exc_h = 0.0_dp
         my_rho_atom_set(:)%exc_s = 0.0_dp
      ELSE
         CALL section_vals_val_get(my_xc_section, "DENSITY_CUTOFF", &
                                   r_val=density_cut)
         CALL section_vals_val_get(my_xc_section, "GRADIENT_CUTOFF", &
                                   r_val=gradient_cut)
         CALL section_vals_val_get(my_xc_section, "TAU_CUTOFF", &
                                   r_val=tau_cut)

         lsd = dft_control%lsd
         nspins = dft_control%nspins
         needs = xc_functionals_get_needs(xc_fun_section, &
                                          lsd=lsd, &
                                          calc_potential=.TRUE.)

         ! whatever the xc, if epr_xc, drho_spin is needed
         needs%drho_spin = .TRUE.

         gradient_f = (needs%drho .OR. needs%drho_spin)
         tau_f = (needs%tau .OR. needs%tau_spin)

         ! Initialize energy contribution from the one center XC terms to zero
         exc1 = 0.0_dp

         ! Nullify some pointers for work-arrays
         NULLIFY (rho_h, drho_h, rho_s, drho_s, weight_h, weight_s)
         NULLIFY (vxc_h, vxc_s, vxg_h, vxg_s)
         NULLIFY (tau_h, tau_s)
         NULLIFY (vtau_h, vtau_s)

         ! Here starts the loop over all the atoms

         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
            CALL get_qs_kind(my_kind_set(ikind), paw_atom=paw_atom, &
                             harmonics=harmonics, grid_atom=grid_atom)
            CALL get_qs_kind(my_kind_set(ikind), basis_set=basis_1c, basis_type="GAPW_1C")

            IF (.NOT. paw_atom) CYCLE

            nr = grid_atom%nr
            na = grid_atom%ng_sphere

            ! Prepare the structures needed to calculate and store the xc derivatives

            ! Array dimension: here anly one dimensional arrays are used,
            ! i.e. only the first column of deriv_data is read.
            ! The other to dimensions  are set to size equal 1
            bounds(1:2, 1:3) = 1
            bounds(2, 1) = na
            bounds(2, 2) = nr

            ! set integration weights
            weight_h => grid_atom%weight
            weight_s => grid_atom%weight
            IF (accint) THEN
               on = dft_control%qs_control%gapw_control%oweights
               alpha = dft_control%qs_control%gapw_control%aw(ikind)
               IF (ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                  IF (grid_atom%gapw_weight_alpha /= alpha) DEALLOCATE (grid_atom%gapw_weight_s)
               END IF
               IF (.NOT. ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                  ALLOCATE (grid_atom%gapw_weight_s(na, nr))
                  ALLOCATE (fw(nr))
                  CALL calc_weight_function(fw, grid_atom%rad2, on, alpha)
                  DO ir = 1, nr
                     agr = 1.0_dp - fw(ir)
                     grid_atom%gapw_weight_s(:, ir) = agr*grid_atom%weight(:, ir)
                  END DO
                  DEALLOCATE (fw)
                  grid_atom%gapw_weight_alpha = alpha
               END IF
               weight_s => grid_atom%gapw_weight_s
            END IF

            ! create a place where to put the derivatives
            CALL xc_dset_create(deriv_set, local_bounds=bounds)
            ! create the place where to store the argument for the functionals
            CALL xc_rho_set_create(rho_set_h, bounds, rho_cutoff=density_cut, &
                                   drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
            CALL xc_rho_set_create(rho_set_s, bounds, rho_cutoff=density_cut, &
                                   drho_cutoff=gradient_cut, tau_cutoff=tau_cut)

            ! allocate the required 3d arrays where to store rho and drho
            CALL xc_rho_set_atom_update(rho_set_h, needs, nspins, bounds)
            CALL xc_rho_set_atom_update(rho_set_s, needs, nspins, bounds)

            CALL reallocate(rho_h, 1, na, 1, nr, 1, nspins)
            CALL reallocate(rho_s, 1, na, 1, nr, 1, nspins)
            CALL reallocate(vxc_h, 1, na, 1, nr, 1, nspins)
            CALL reallocate(vxc_s, 1, na, 1, nr, 1, nspins)
            !
            IF (gradient_f) THEN
               CALL reallocate(drho_h, 1, 4, 1, na, 1, nr, 1, nspins)
               CALL reallocate(drho_s, 1, 4, 1, na, 1, nr, 1, nspins)
               CALL reallocate(vxg_h, 1, 3, 1, na, 1, nr, 1, nspins)
               CALL reallocate(vxg_s, 1, 3, 1, na, 1, nr, 1, nspins)
            END IF

            IF (tau_f) THEN
               CALL create_tau_basis_cache(tau_basis_cache, grid_atom, basis_1c, harmonics)
               CALL reallocate(tau_h, 1, na, 1, nr, 1, nspins)
               CALL reallocate(tau_s, 1, na, 1, nr, 1, nspins)
               CALL reallocate(vtau_h, 1, na, 1, nr, 1, nspins)
               CALL reallocate(vtau_s, 1, na, 1, nr, 1, nspins)
            END IF

            ! NLCC: prepare rho and drho of the core charge for this KIND
            donlcc = .FALSE.
            IF (nlcc) THEN
               NULLIFY (rho_nlcc)
               rho_nlcc => my_kind_set(ikind)%nlcc_pot
               IF (ASSOCIATED(rho_nlcc)) donlcc = .TRUE.
            END IF

            ! Distribute the atoms of this kind

            num_pe = para_env%num_pe
            bo = get_limit(natom, para_env%num_pe, para_env%mepos)

            DO iat = bo(1), bo(2)
               iatom = atom_list(iat)

               my_rho_atom_set(iatom)%exc_h = 0.0_dp
               my_rho_atom_set(iatom)%exc_s = 0.0_dp

               rho_atom => my_rho_atom_set(iatom)
               rho_h = 0.0_dp
               rho_s = 0.0_dp
               IF (gradient_f) THEN
                  NULLIFY (r_h, r_s, dr_h, dr_s, r_h_d, r_s_d)
                  CALL get_rho_atom(rho_atom=rho_atom, rho_rad_h=r_h, &
                                    rho_rad_s=r_s, drho_rad_h=dr_h, &
                                    drho_rad_s=dr_s, rho_rad_h_d=r_h_d, &
                                    rho_rad_s_d=r_s_d)
                  drho_h = 0.0_dp
                  drho_s = 0.0_dp
               ELSE
                  NULLIFY (r_h, r_s)
                  CALL get_rho_atom(rho_atom=rho_atom, rho_rad_h=r_h, rho_rad_s=r_s)
                  rho_d = 0.0_dp
               END IF
               IF (tau_f) THEN
                  !compute tau on the grid all at once
                  CALL calc_tau_atom(tau_h, tau_s, rho_atom, tau_basis_cache, nspins)
               ELSE
                  tau_d = 0.0_dp
               END IF

               DO ir = 1, nr
                  CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_f, &
                                        ir, r_h, r_s, rho_h, rho_s, dr_h, dr_s, &
                                        r_h_d, r_s_d, drho_h, drho_s)
                  IF (donlcc) THEN
                     CALL calc_rho_nlcc(grid_atom, nspins, gradient_f, &
                                        ir, rho_nlcc(:, 1), rho_h, rho_s, rho_nlcc(:, 2), drho_h, drho_s)
                  END IF
               END DO
               DO ir = 1, nr
                  IF (tau_f) THEN
                     CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, drho_h, tau_h, na, ir)
                     CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, drho_s, tau_s, na, ir)
                  ELSE IF (gradient_f) THEN
                     CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, drho_h, tau_d, na, ir)
                     CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, drho_s, tau_d, na, ir)
                  ELSE
                     CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, rho_d, tau_d, na, ir)
                     CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, rho_d, tau_d, na, ir)
                  END IF
               END DO

               !-------------------!
               ! hard atom density !
               !-------------------!
               CALL xc_dset_zero_all(deriv_set)
               CALL vxc_of_r_epr(xc_fun_section, rho_set_h, deriv_set, needs, weight_h, &
                                 lsd, na, nr, exc_h, vxc_h, vxg_h, vtau_h)
               rho_atom%exc_h = rho_atom%exc_h + exc_h

               !-------------------!
               ! soft atom density !
               !-------------------!
               CALL xc_dset_zero_all(deriv_set)
               CALL vxc_of_r_epr(xc_fun_section, rho_set_s, deriv_set, needs, weight_s, &
                                 lsd, na, nr, exc_s, vxc_s, vxg_s, vtau_s)
               rho_atom%exc_s = rho_atom%exc_s + exc_s

               DO ispin = 1, nspins
                  DO idir = 1, 3
                     DO ir = 1, nr
                        DO ia = 1, na
                           gradient_atom_set(iatom)%nablavks_vec_rad_h(idir, ispin)%r_coef(ir, ia) = &
                              gradient_atom_set(iatom)%nablavks_vec_rad_h(idir, ispin)%r_coef(ir, ia) &
                              + vxg_h(idir, ia, ir, ispin)
                           gradient_atom_set(iatom)%nablavks_vec_rad_s(idir, ispin)%r_coef(ir, ia) = &
                              gradient_atom_set(iatom)%nablavks_vec_rad_s(idir, ispin)%r_coef(ir, ia) &
                              + vxg_s(idir, ia, ir, ispin)
                        END DO ! ia
                     END DO ! ir
                  END DO ! idir
               END DO ! ispin

               ! Add contributions to the exc energy

               exc1 = exc1 + rho_atom%exc_h - rho_atom%exc_s

               ! Integration to get the matrix elements relative to the vxc_atom
               ! here the products with the primitives is done: gaVxcgb
               ! internal transformation to get the integral in cartesian Gaussians

               NULLIFY (int_hh, int_ss)
               CALL get_rho_atom(rho_atom=rho_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
               IF (gradient_f) THEN
                  CALL gaVxcgb_GC(vxc_h, vxc_s, vxg_h, vxg_s, int_hh, int_ss, &
                                  grid_atom, basis_1c, harmonics, nspins)
               ELSE
                  CALL gaVxcgb_noGC(vxc_h, vxc_s, int_hh, int_ss, &
                                    grid_atom, basis_1c, harmonics, nspins)
               END IF
               IF (tau_f) THEN
                  CALL dgaVtaudgb(vtau_h, vtau_s, int_hh, int_ss, &
                                  tau_basis_cache, nspins)
               END IF
               NULLIFY (r_h, r_s, dr_h, dr_s)
            END DO ! iat

            IF (tau_f) CALL release_tau_basis_cache(tau_basis_cache)

            ! Release the xc structure used to store the xc derivatives
            CALL xc_dset_release(deriv_set)
            CALL xc_rho_set_release(rho_set_h)
            CALL xc_rho_set_release(rho_set_s)
         END DO ! ikind

         CALL para_env%sum(exc1)

         IF (ASSOCIATED(rho_h)) DEALLOCATE (rho_h)
         IF (ASSOCIATED(rho_s)) DEALLOCATE (rho_s)
         IF (ASSOCIATED(vxc_h)) DEALLOCATE (vxc_h)
         IF (ASSOCIATED(vxc_s)) DEALLOCATE (vxc_s)

         IF (gradient_f) THEN
            IF (ASSOCIATED(drho_h)) DEALLOCATE (drho_h)
            IF (ASSOCIATED(drho_s)) DEALLOCATE (drho_s)
            IF (ASSOCIATED(vxg_h)) DEALLOCATE (vxg_h)
            IF (ASSOCIATED(vxg_s)) DEALLOCATE (vxg_s)
         END IF

         IF (tau_f) THEN
            IF (ASSOCIATED(tau_h)) DEALLOCATE (tau_h)
            IF (ASSOCIATED(tau_s)) DEALLOCATE (tau_s)
            IF (ASSOCIATED(vtau_h)) DEALLOCATE (vtau_h)
            IF (ASSOCIATED(vtau_s)) DEALLOCATE (vtau_s)
         END IF

      END IF !xc_none

      CALL timestop(handle)

   END SUBROUTINE calculate_vxc_atom_epr

END MODULE qs_vxc_atom

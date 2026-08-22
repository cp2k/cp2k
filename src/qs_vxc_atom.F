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
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_array_utils,                  ONLY: cp_3d_r_cp_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type,&
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
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indco,&
                                              indso,&
                                              nco,&
                                              ncoset,&
                                              nsoset
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_types,                  ONLY: particle_type
   USE paw_basis_types,                 ONLY: get_paw_basis_info
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
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
   USE qs_harmonics_atom,               ONLY: get_none0_cg_list,&
                                              harmonics_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              has_nlcc,&
                                              qs_kind_type
   USE qs_linres_types,                 ONLY: nablavks_atom_type
   USE qs_rho_atom_methods,             ONLY: replicate_rho_atom_radial
   USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                              rho_atom_coeff,&
                                              rho_atom_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE skala_gpw_features,              ONLY: build_periodic_atom_image_layout,&
                                              periodic_atom_image_partition,&
                                              periodic_atom_image_partition_from_layout,&
                                              skala_gpw_smooth_partition_derivatives,&
                                              smooth_atom_partition,&
                                              smooth_partition_atomic_weight_scale,&
                                              smooth_partition_atomic_weight_scale_derivative
   USE skala_gpw_functional,            ONLY: &
        build_vxc_from_feature_grads, get_gauxc_section, native_skala_gapw_density_partition, &
        skala_gapw_atom_composite_energy, skala_gapw_atom_vxc_of_r, &
        skala_gapw_density_partition_hard_minus_soft, skala_gapw_density_partition_hard_only, &
        skala_gapw_density_partition_none, skala_gapw_density_partition_soft_only, &
        xc_section_uses_gauxc_model, xc_section_uses_native_skala_evaluator
   USE spherical_harmonics,             ONLY: y_lm
   USE util,                            ONLY: get_limit
   USE virial_types,                    ONLY: virial_type
   USE xc_atom,                         ONLY: fill_rho_set,&
                                              vxc_of_r_epr,&
                                              vxc_of_r_new,&
                                              xc_2nd_deriv_of_r,&
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
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_vxc_atom'

   ! The wider stencil suppresses mixed-kind eggbox forces while preserving an exact transpose.
   INTEGER, PARAMETER, PRIVATE :: native_grid_interp_offset_min = -5, &
                                  native_grid_interp_offset_max = 6, &
                                  native_grid_interp_npts = &
                                  native_grid_interp_offset_max - &
                                  native_grid_interp_offset_min + 1
   ! A wrapped stencil can touch both end tiles and one adjacent interior tile when the final
   ! tile is shorter than the stencil. Three tiles per direction are therefore sufficient.
   INTEGER, PARAMETER, PRIVATE :: native_grid_adjoint_tile_edge = 64, &
                                  native_grid_adjoint_max_tiles_per_direction = 3, &
                                  native_grid_adjoint_max_bins_per_row = &
                                  native_grid_adjoint_max_tiles_per_direction**3

   TYPE native_grid_interpolation_stencil_type
      INTEGER, DIMENSION(native_grid_interp_npts, 3)     :: relative_index = 0
      REAL(KIND=dp), DIMENSION(native_grid_interp_npts, 3) :: weight = 0.0_dp
      LOGICAL, DIMENSION(native_grid_interp_npts, 3)     :: valid = .FALSE.
      LOGICAL                                            :: active = .FALSE.
   END TYPE native_grid_interpolation_stencil_type

   TYPE tau_basis_cache_type
      INTEGER                                            :: maxso = 0, na = 0, nr = 0, nsatbas = 0, &
                                                            nset = 0
      INTEGER, DIMENSION(:), POINTER                     :: lmax => NULL(), lmin => NULL(), &
                                                            n2oindex => NULL(), npgf => NULL(), &
                                                            o2nindex => NULL()
      REAL(dp), DIMENSION(:, :), POINTER                 :: zet => NULL()
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: grad
   END TYPE tau_basis_cache_type

   PUBLIC :: calculate_vxc_atom, &
             calculate_vxc_atom_epr, &
             calculate_xc_2nd_deriv_atom, &
             calc_rho_angular, &
             calculate_gfxc_atom, &
             gapw_cdft_one_center, &
             gfxc_atom_diff, &
             gaVxcgb_noGC

CONTAINS

! **************************************************************************************************
!> \brief Evaluate an NLCC density and its first two Cartesian derivatives.
!> \param point evaluation point
!> \param center pseudopotential center
!> \param gth_potential optional GTH potential
!> \param sgp_potential optional separable Gaussian potential
!> \param rho core density
!> \param gradient Cartesian density gradient
!> \param hessian Cartesian density Hessian
! **************************************************************************************************
   SUBROUTINE evaluate_nlcc_primitive_fields(point, center, gth_potential, sgp_potential, &
                                             rho, gradient, hessian)
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: point, center
      TYPE(gth_potential_type), INTENT(IN), POINTER      :: gth_potential
      TYPE(sgp_potential_type), INTENT(IN), POINTER      :: sgp_potential
      REAL(dp), INTENT(OUT)                              :: rho
      REAL(dp), DIMENSION(3), INTENT(OUT)                :: gradient
      REAL(dp), DIMENSION(3, 3), INTENT(OUT)             :: hessian

      INTEGER                                            :: ic, idir, iexp, jdir, n_nlcc, nexp_nlcc, &
                                                            power
      INTEGER, DIMENSION(:), POINTER                     :: nct_nlcc
      LOGICAL                                            :: has_sgp_nlcc, nlcc_present
      REAL(dp)                                           :: alpha, beta, d2poly, dpoly, exponential, &
                                                            poly, r2, rho_x, rho_xx, scaled_r2
      REAL(dp), DIMENSION(3)                             :: displacement
      REAL(dp), DIMENSION(:), POINTER                    :: a_nlcc, alpha_nlcc, c_nlcc
      REAL(dp), DIMENSION(:, :), POINTER                 :: cval_nlcc

      NULLIFY (a_nlcc, alpha_nlcc, c_nlcc, cval_nlcc, nct_nlcc)
      rho = 0.0_dp
      rho_x = 0.0_dp
      rho_xx = 0.0_dp
      displacement = point - center
      r2 = DOT_PRODUCT(displacement, displacement)

      IF (ASSOCIATED(gth_potential)) THEN
         CALL get_potential(gth_potential, nlcc_present=nlcc_present, &
                            nexp_nlcc=nexp_nlcc, alpha_nlcc=alpha_nlcc, &
                            nct_nlcc=nct_nlcc, cval_nlcc=cval_nlcc)
         IF (nlcc_present) THEN
            DO iexp = 1, nexp_nlcc
               alpha = alpha_nlcc(iexp)
               beta = 0.5_dp/(alpha*alpha)
               scaled_r2 = r2/(alpha*alpha)
               exponential = EXP(-0.5_dp*scaled_r2)
               DO ic = 1, nct_nlcc(iexp)
                  power = ic - 1
                  poly = cval_nlcc(ic, iexp)*scaled_r2**power
                  dpoly = 0.0_dp
                  IF (power > 0) THEN
                     dpoly = cval_nlcc(ic, iexp)*REAL(power, dp)* &
                             scaled_r2**(power - 1)/(alpha*alpha)
                  END IF
                  d2poly = 0.0_dp
                  IF (power > 1) THEN
                     d2poly = cval_nlcc(ic, iexp)*REAL(power*(power - 1), dp)* &
                              scaled_r2**(power - 2)/(alpha**4)
                  END IF
                  rho = rho + exponential*poly
                  rho_x = rho_x + exponential*(dpoly - beta*poly)
                  rho_xx = rho_xx + exponential*(d2poly - 2.0_dp*beta*dpoly + beta*beta*poly)
               END DO
            END DO
         END IF
      ELSE IF (ASSOCIATED(sgp_potential)) THEN
         CALL get_potential(sgp_potential, has_nlcc=has_sgp_nlcc, n_nlcc=n_nlcc, &
                            a_nlcc=a_nlcc, c_nlcc=c_nlcc)
         IF (has_sgp_nlcc) THEN
            DO iexp = 1, n_nlcc
               exponential = EXP(-a_nlcc(iexp)*r2)
               rho = rho + c_nlcc(iexp)*exponential
               rho_x = rho_x - a_nlcc(iexp)*c_nlcc(iexp)*exponential
               rho_xx = rho_xx + a_nlcc(iexp)**2*c_nlcc(iexp)*exponential
            END DO
         END IF
      END IF

      gradient = 2.0_dp*rho_x*displacement
      DO idir = 1, 3
         DO jdir = 1, 3
            hessian(idir, jdir) = 4.0_dp*rho_xx*displacement(idir)*displacement(jdir)
         END DO
         hessian(idir, idir) = hessian(idir, idir) + 2.0_dp*rho_x
      END DO
   END SUBROUTINE evaluate_nlcc_primitive_fields

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

      INTEGER :: adjoint_bin, adjoint_entry, adjoint_nbins, adjoint_nchannels, &
         adjoint_tile_count(3), adjoint_tile_lower(3), adjoint_tile_upper(3), &
         atom_composite_components, base_shift(3), bo(2), composite_descriptor_target_image, &
         composite_image_periodicity(3), composite_local_atom, composite_local_natom, &
         composite_nflat, composite_partition_target_image, composite_pw_nflat, composite_row, &
         gapw_density_partition, gapw_representation, handle, ia, iat, iatom, icomponent, idir, &
         ikind, image_i1, image_i2, image_i3, image_shell(3), image_shift(3), ir, ispin, iw, jdir, &
         myfun, na, natom, nr, nspins
      INTEGER :: num_pe, source_atom, target_atom, xc_deriv_method_id, xc_rho_smooth_id, zatom
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: composite_atomic_grid_sizes, &
                                                            composite_local_grid_sizes
      INTEGER, ALLOCATABLE, DIMENSION(:) :: adjoint_bin_offsets, adjoint_bin_rows, &
         composite_atom_end, composite_atom_kind, composite_atom_kind_index, composite_atom_start, &
         composite_grid_atom, composite_local_atoms
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL :: accint, atom_composite_active, atom_composite_diagnostic, &
         atom_composite_reference, direct_valence_atom_composite, donlcc, evaluate_hard, &
         evaluate_soft, gradient_f, image_partition_atom_composite, lsd, my_calculate_forces, &
         native_grid_diagnostics, nlcc, one_center_kind, paw_atom, paw_pseudopotentials, &
         requested_atom_composite_grid, rho_g_valid, skala_atom_grid, source_matrix_local, tau_f, &
         tau_r_valid, use_atom_composite_density, use_atom_composite_gradient, &
         use_atom_composite_tau, use_virial
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: composite_partition_included
      REAL(dp) :: agr, alpha, atom_composite_exc, atom_composite_nelec, &
         composite_cross_cutoff_max, composite_cross_density_max, composite_cross_grad_max, &
         composite_cross_kin_max, composite_density_max, composite_density_min, &
         composite_grad_max, composite_kin_max, composite_kin_min, composite_tau_integral, &
         cross_cutoff, density_cut, descriptor_window_adjoint, descriptor_window_weight, exc_h, &
         exc_s, feature_vxc_analytic, feature_vxc_fd, feature_vxc_minus, feature_vxc_plus, &
         feature_vxc_step, gradient_cut, local_partition_weight, my_adiabatic_rescale_factor, &
         nlcc_density, nlcc_spin_factor
      REAL(dp) :: one_center_density_field_contraction, one_center_density_matrix_contraction, &
         one_center_field_contraction, one_center_gradient_field_contraction, &
         one_center_gradient_matrix_contraction, one_center_matrix_contraction, &
         one_center_rho_grad_field_contraction, one_center_rho_grad_matrix_contraction, &
         one_center_tau_field_contraction, one_center_tau_matrix_contraction, &
         one_center_tensor_contraction, partition_adjoint, partition_scale, partition_weight, &
         smooth_grid_contraction, smooth_input_contraction, target_partition_adjoint, tau_cut, zeff
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: composite_atomic_grid_weight_grad, &
         composite_atomic_grid_weights, composite_base_grid_weights, composite_distances, &
         composite_grid_weight_grad, composite_grid_weights, composite_partition_weights
      REAL(dp), ALLOCATABLE, DIMENSION(:, :) :: composite_atom_coord_grad, composite_atom_coords, &
         composite_cross_density, composite_cross_force, composite_cross_force_local, &
         composite_cross_kin, composite_density, composite_density_grad, &
         composite_descriptor_image_coords, composite_explicit_force, composite_grid_coord_force, &
         composite_grid_coord_grad, composite_grid_coords, composite_kin, composite_kin_grad, &
         composite_local_atom_coords, composite_model_atom_force, composite_moving_smooth_force, &
         composite_nlcc_center_force, composite_nlcc_center_force_local, &
         composite_nlcc_target_force
      REAL(dp), ALLOCATABLE, DIMENSION(:, :) :: composite_partition_atom_coords, &
         composite_partition_force, composite_partition_force_local, &
         composite_partition_image_coords, composite_smooth_density_cache, &
         composite_smooth_kin_cache, local_partition_datom
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), TARGET     :: smooth_density_adjoint_storage, &
                                                            smooth_kin_adjoint_storage
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :) :: composite_cross_grad, composite_grad, &
         composite_grad_grad, composite_int_h, composite_int_s, composite_partition_datom, &
         composite_partition_dstrain, composite_smooth_gradient_cache
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), TARGET  :: smooth_grad_adjoint_storage
      REAL(dp), DIMENSION(1, 1, 1)                       :: tau_d
      REAL(dp), DIMENSION(1, 1, 1, 1)                    :: rho_d
      REAL(dp), DIMENSION(2) :: composite_smooth_density_adjoint_value, &
         composite_smooth_density_value, composite_smooth_kin_adjoint_value, &
         composite_smooth_kin_value, cross_density, cross_density_adjoint, cross_kin, &
         cross_kin_adjoint
      REAL(dp), DIMENSION(3) :: composite_point, cross_displacement, cross_spatial_derivative, &
         fractional, image_translation, nlcc_gradient, nlcc_spatial_derivative, &
         skala_atom_force_h, skala_atom_force_s, spatial_derivative
      REAL(dp), DIMENSION(3, 1)                          :: local_descriptor_datom
      REAL(dp), DIMENSION(3, 2) :: composite_smooth_gradient_adjoint_value, &
         composite_smooth_gradient_value, cross_density_spatial, cross_grad, cross_grad_adjoint, &
         cross_kin_spatial
      REAL(dp), DIMENSION(3, 3) :: composite_cross_image_virial, &
         composite_cross_image_virial_local, composite_explicit_virial, composite_feature_virial, &
         composite_interpolation_virial, composite_partition_strain_virial, &
         local_descriptor_dstrain, local_partition_dstrain, nlcc_hessian, skala_atom_virial, &
         skala_atom_virial_h, skala_atom_virial_s
      REAL(dp), DIMENSION(3, 3, 2)                       :: cross_grad_spatial
      REAL(dp), DIMENSION(4)                             :: feature_component_analytic, &
                                                            feature_component_fd
      REAL(dp), DIMENSION(:, :), POINTER                 :: rho_nlcc, smooth_density_adjoint, &
                                                            smooth_kin_adjoint, weight_h, weight_s
      REAL(dp), DIMENSION(:, :, :), POINTER :: composite_smooth_rho, composite_smooth_rhoa, &
         composite_smooth_rhob, composite_smooth_tau, composite_smooth_tau_a, &
         composite_smooth_tau_b, rho_h, rho_s, smooth_grad_adjoint, smooth_rho, smooth_rhoa, &
         smooth_rhob, smooth_tau, smooth_tau_a, smooth_tau_b, tau_h, tau_s, vtau_h, vtau_s, vxc_h, &
         vxc_s
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: drho_h, drho_s, vxg_h, vxg_s
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_3d_r_cp_type), DIMENSION(3) :: composite_smooth_drho, composite_smooth_drhoa, &
         composite_smooth_drhob, smooth_drho, smooth_drhoa, smooth_drhob
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
      TYPE(section_vals_type), POINTER                   :: gauxc_section, input, my_xc_section, &
                                                            xc_fun_section
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential
      TYPE(tau_basis_cache_type)                         :: tau_basis_cache
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set_h, rho_set_s, smooth_rho_set

! -------------------------------------------------------------------------

      CALL timeset(routineN, handle)

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
      NULLIFY (composite_smooth_rho, composite_smooth_rhoa, composite_smooth_rhob, &
               composite_smooth_tau, composite_smooth_tau_a, composite_smooth_tau_b)
      NULLIFY (smooth_rho_g, smooth_rho_r, smooth_tau_r)
      NULLIFY (smooth_vxc_rho, smooth_vxc_tau)
      DO idir = 1, 3
         NULLIFY (smooth_drho(idir)%array, smooth_drhoa(idir)%array, smooth_drhob(idir)%array)
         NULLIFY (composite_smooth_drho(idir)%array, &
                  composite_smooth_drhoa(idir)%array, &
                  composite_smooth_drhob(idir)%array)
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
         gauxc_section => get_gauxc_section(my_xc_section)
         CPASSERT(ASSOCIATED(gauxc_section))
         CALL section_vals_val_get(gauxc_section, "PSEUDOPOTENTIAL_GAPW_REPRESENTATION", &
                                   i_val=gapw_representation)
         CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_DIAGNOSTICS", &
                                   l_val=native_grid_diagnostics)
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
      IF (skala_atom_grid .AND. xc_section_uses_native_skala_evaluator(my_xc_section)) THEN
         CALL section_vals_val_get(gauxc_section, &
                                   "NATIVE_GRID_GAPW_ATOM_COMPOSITE_DIAGNOSTIC", &
                                   l_val=atom_composite_diagnostic)
         CALL section_vals_val_get(gauxc_section, &
                                   "NATIVE_GRID_GAPW_ATOM_COMPOSITE_REFERENCE", &
                                   l_val=atom_composite_reference)
         CALL section_vals_val_get(gauxc_section, &
                                   "NATIVE_GRID_GAPW_ATOM_COMPOSITE_COMPONENTS", &
                                   i_val=atom_composite_components)
         CALL section_vals_val_get(gauxc_section, &
                                   "NATIVE_GRID_GAPW_ATOM_COMPOSITE_FD_STEP", &
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
      image_partition_atom_composite = atom_composite_active
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
                      composite_atom_coords(3, SIZE(particle_set)), &
                      composite_partition_weights(SIZE(particle_set)), &
                      composite_partition_atom_coords(3, SIZE(particle_set)), &
                      composite_distances(SIZE(particle_set)))
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
                      composite_smooth_density_cache(composite_nflat, 2), &
                      composite_smooth_gradient_cache(composite_nflat, 3, 2), &
                      composite_smooth_kin_cache(composite_nflat, 2), &
                      composite_grid_coords(3, composite_nflat), &
                      composite_grid_weights(composite_nflat), &
                      composite_base_grid_weights(composite_nflat), &
                      composite_atomic_grid_weights(composite_nflat))
            composite_density = 0.0_dp
            composite_grad = 0.0_dp
            composite_kin = 0.0_dp
            DO composite_local_atom = 1, composite_local_natom
               iatom = composite_local_atoms(composite_local_atom)
               composite_grid_atom(composite_atom_start(iatom):composite_atom_end(iatom)) = iatom
            END DO
            composite_smooth_density_cache = 0.0_dp
            composite_smooth_gradient_cache = 0.0_dp
            composite_smooth_kin_cache = 0.0_dp
            composite_grid_coords = 0.0_dp
            composite_grid_weights = 0.0_dp
            composite_base_grid_weights = 0.0_dp
            composite_atomic_grid_weights = 0.0_dp

            CALL qs_rho_get(rho_struct, rho_r=smooth_rho_r, rho_g=smooth_rho_g, &
                            tau_r=smooth_tau_r, rho_g_valid=rho_g_valid, &
                            tau_r_valid=tau_r_valid)
            CPASSERT(rho_g_valid)
            CPASSERT(tau_r_valid)
            CPASSERT(ASSOCIATED(smooth_rho_r))
            CPASSERT(ASSOCIATED(smooth_rho_g))
            CPASSERT(ASSOCIATED(smooth_tau_r))
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
            IF (lsd) THEN
               CALL xc_rho_set_get(smooth_rho_set, rhoa=smooth_rhoa, rhob=smooth_rhob, &
                                   drhoa=smooth_drhoa, drhob=smooth_drhob, &
                                   tau_a=smooth_tau_a, tau_b=smooth_tau_b)
               CALL gather_native_grid_field(smooth_rhoa, smooth_rho_r(1)%pw_grid, para_env, &
                                             composite_smooth_rhoa)
               CALL gather_native_grid_field(smooth_rhob, smooth_rho_r(1)%pw_grid, para_env, &
                                             composite_smooth_rhob)
               CALL gather_native_grid_field(smooth_tau_a, smooth_rho_r(1)%pw_grid, para_env, &
                                             composite_smooth_tau_a)
               CALL gather_native_grid_field(smooth_tau_b, smooth_rho_r(1)%pw_grid, para_env, &
                                             composite_smooth_tau_b)
               DO idir = 1, 3
                  CALL gather_native_grid_field(smooth_drhoa(idir)%array, &
                                                smooth_rho_r(1)%pw_grid, para_env, &
                                                composite_smooth_drhoa(idir)%array)
                  CALL gather_native_grid_field(smooth_drhob(idir)%array, &
                                                smooth_rho_r(1)%pw_grid, para_env, &
                                                composite_smooth_drhob(idir)%array)
               END DO
            ELSE
               CALL xc_rho_set_get(smooth_rho_set, rho=smooth_rho, drho=smooth_drho, &
                                   tau=smooth_tau)
               CALL gather_native_grid_field(smooth_rho, smooth_rho_r(1)%pw_grid, para_env, &
                                             composite_smooth_rho)
               CALL gather_native_grid_field(smooth_tau, smooth_rho_r(1)%pw_grid, para_env, &
                                             composite_smooth_tau)
               DO idir = 1, 3
                  CALL gather_native_grid_field(smooth_drho(idir)%array, &
                                                smooth_rho_r(1)%pw_grid, para_env, &
                                                composite_smooth_drho(idir)%array)
               END DO
            END IF
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
               IF (accint) THEN
                  weight_h => grid_atom%weight
                  alpha = dft_control%qs_control%gapw_control%aw(ikind)
                  IF (ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                     IF (grid_atom%gapw_weight_alpha /= alpha) DEALLOCATE (grid_atom%gapw_weight_s)
                  END IF
                  IF (.NOT. ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                     ALLOCATE (grid_atom%gapw_weight_s(na, nr))
                     DO ir = 1, nr
                        agr = 1.0_dp - EXP(-alpha*grid_atom%rad2(ir))
                        grid_atom%gapw_weight_s(:, ir) = grid_atom%weight(:, ir)*agr
                     END DO
                     grid_atom%gapw_weight_alpha = alpha
                  END IF
                  weight_s => grid_atom%gapw_weight_s
               ELSE
                  weight_h => grid_atom%weight
                  weight_s => grid_atom%weight
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
                  IF (image_partition_atom_composite) THEN
                     CALL build_periodic_atom_image_layout( &
                        composite_atom_coords, cell, iatom, composite_image_periodicity, &
                        composite_partition_image_coords, composite_partition_target_image)
                     CALL build_periodic_atom_image_layout( &
                        composite_atom_coords(:, iatom:iatom), cell, 1, &
                        composite_image_periodicity, composite_descriptor_image_coords, &
                        composite_descriptor_target_image)
                  END IF
!$OMP PARALLEL DO COLLAPSE(2) IF (image_partition_atom_composite) SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP PRIVATE(composite_row, composite_point, composite_smooth_density_value, &
!$OMP         composite_smooth_gradient_value, composite_smooth_kin_value, &
!$OMP         descriptor_window_weight, idir, ispin, &
!$OMP         gth_potential, nlcc_density, nlcc_gradient, nlcc_hessian, nlcc_spin_factor, &
!$OMP         interpolation_stencil, partition_scale, partition_weight, sgp_potential, source_atom) &
!$OMP SHARED(atom_composite_reference, cell, composite_atom_coords, composite_atom_kind, &
!$OMP        composite_atom_start, composite_atomic_grid_weights, composite_base_grid_weights, &
!$OMP        composite_density, composite_descriptor_image_coords, &
!$OMP        composite_descriptor_target_image, composite_distances, composite_grad, &
!$OMP        composite_grid_coords, composite_grid_weights, composite_kin, &
!$OMP        composite_partition_atom_coords, composite_partition_image_coords, &
!$OMP        composite_partition_target_image, composite_partition_weights, &
!$OMP        composite_smooth_density_cache, composite_smooth_drho, composite_smooth_drhoa, &
!$OMP        composite_smooth_drhob, composite_smooth_gradient_cache, &
!$OMP        composite_smooth_kin_cache, composite_smooth_rho, composite_smooth_rhoa, &
!$OMP        composite_smooth_rhob, composite_smooth_tau, composite_smooth_tau_a, &
!$OMP        composite_smooth_tau_b, drho_h, drho_s, grid_atom, iatom, &
!$OMP        image_partition_atom_composite, lsd, my_kind_set, na, nlcc, nr, one_center_kind, &
!$OMP        particle_set, rho_h, rho_s, smooth_rho_r, tau_h, tau_s, &
!$OMP        use_atom_composite_density, use_atom_composite_gradient, use_atom_composite_tau)
                  DO ir = 1, nr
                     DO ia = 1, na
                        composite_row = composite_atom_start(iatom) + (ir - 1)*na + ia - 1
                        composite_point(1) = particle_set(iatom)%r(1) + grid_atom%rad(ir)* &
                                             grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia)
                        composite_point(2) = particle_set(iatom)%r(2) + grid_atom%rad(ir)* &
                                             grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia)
                        composite_point(3) = particle_set(iatom)%r(3) + &
                                             grid_atom%rad(ir)*grid_atom%cos_pol(ia)
                        composite_grid_coords(:, composite_row) = composite_point
                        IF (image_partition_atom_composite) THEN
                           CALL periodic_atom_image_partition_from_layout( &
                              composite_point, composite_partition_image_coords, &
                              composite_partition_target_image, partition_weight)
                           ! The self-image partition defines a smooth atom-centered periodic
                           ! descriptor domain without truncating it at neighboring atoms.
                           CALL periodic_atom_image_partition_from_layout( &
                              composite_point, composite_descriptor_image_coords, &
                              composite_descriptor_target_image, descriptor_window_weight)
                           partition_scale = smooth_partition_atomic_weight_scale( &
                                             descriptor_window_weight)
                        ELSE
                           CALL smooth_atom_partition( &
                              composite_point, composite_atom_coords, cell, &
                              composite_partition_weights, composite_partition_atom_coords, &
                              composite_distances)
                           partition_weight = composite_partition_weights(iatom)
                           partition_scale = 1.0_dp
                        END IF
                        composite_base_grid_weights(composite_row) = grid_atom%weight(ia, ir)
                        composite_atomic_grid_weights(composite_row) = &
                           composite_base_grid_weights(composite_row)*partition_scale
                        composite_grid_weights(composite_row) = &
                           composite_base_grid_weights(composite_row)* &
                           partition_weight
                        CALL create_native_grid_interpolation_stencil( &
                           interpolation_stencil, smooth_rho_r(1)%pw_grid, cell, composite_point, &
                           image_partition_atom_composite)
                        IF (lsd) THEN
                           CALL interpolate_native_grid_fields( &
                              composite_smooth_rhoa, composite_smooth_drhoa(1)%array, &
                              composite_smooth_drhoa(2)%array, composite_smooth_drhoa(3)%array, &
                              composite_smooth_tau_a, interpolation_stencil, &
                              composite_smooth_density_value(1), &
                              composite_smooth_gradient_value(:, 1), &
                              composite_smooth_kin_value(1))
                           CALL interpolate_native_grid_fields( &
                              composite_smooth_rhob, composite_smooth_drhob(1)%array, &
                              composite_smooth_drhob(2)%array, composite_smooth_drhob(3)%array, &
                              composite_smooth_tau_b, interpolation_stencil, &
                              composite_smooth_density_value(2), &
                              composite_smooth_gradient_value(:, 2), &
                              composite_smooth_kin_value(2))
                           composite_smooth_density_cache(composite_row, :) = &
                              composite_smooth_density_value
                           composite_smooth_gradient_cache(composite_row, :, :) = &
                              composite_smooth_gradient_value
                           composite_smooth_kin_cache(composite_row, :) = &
                              composite_smooth_kin_value
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
                           CALL interpolate_native_grid_fields( &
                              composite_smooth_rho, composite_smooth_drho(1)%array, &
                              composite_smooth_drho(2)%array, composite_smooth_drho(3)%array, &
                              composite_smooth_tau, interpolation_stencil, &
                              composite_smooth_density_value(1), &
                              composite_smooth_gradient_value(:, 1), &
                              composite_smooth_kin_value(1))
                           composite_smooth_density_cache(composite_row, 1) = &
                              composite_smooth_density_value(1)
                           composite_smooth_gradient_cache(composite_row, :, 1) = &
                              composite_smooth_gradient_value(:, 1)
                           composite_smooth_kin_cache(composite_row, 1) = &
                              composite_smooth_kin_value(1)
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
                           nlcc_spin_factor = MERGE(1.0_dp, 0.5_dp, lsd)
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
                  IF (image_partition_atom_composite) THEN
                     DEALLOCATE (composite_descriptor_image_coords, composite_partition_image_coords)
                  END IF
                  CPASSERT(nr*na == composite_atom_end(iatom) - composite_atom_start(iatom) + 1)
               END IF

               IF (.NOT. one_center_kind) CYCLE

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
                                   "Unknown GAUXC%NATIVE_GRID_GAPW_DENSITY_PARTITION value.")
                  END SELECT
               END IF
               IF (atom_composite_reference) THEN
                  evaluate_hard = .FALSE.
                  evaluate_soft = .FALSE.
               END IF

               !-------------------!
               ! hard atom density !
               !-------------------!
               CALL xc_dset_zero_all(deriv_set)
               IF (.NOT. evaluate_hard) THEN
                  exc_h = 0.0_dp
                  IF (.NOT. energy_only) THEN
                     vxc_h = 0.0_dp
                     IF (ASSOCIATED(vxg_h)) vxg_h = 0.0_dp
                     IF (ASSOCIATED(vtau_h)) vtau_h = 0.0_dp
                  END IF
               ELSE IF (skala_atom_grid) THEN
                  CALL skala_gapw_atom_vxc_of_r( &
                     my_xc_section, grid_atom, para_env, particle_set(iatom)%r, &
                     rho_h, drho_h, tau_h, weight_h, lsd, nspins, na, nr, &
                     exc_h, vxc_h, vxg_h, vtau_h, energy_only=energy_only, &
                     atom_force=skala_atom_force_h, atom_virial=skala_atom_virial_h)
               ELSE
                  CALL vxc_of_r_new(xc_fun_section, rho_set_h, deriv_set, 1, needs, weight_h, &
                                    lsd, na, nr, exc_h, vxc_h, vxg_h, vtau_h, energy_only=energy_only, &
                                    adiabatic_rescale_factor=my_adiabatic_rescale_factor)
               END IF
               rho_atom%exc_h = rho_atom%exc_h + exc_h

               !-------------------!
               ! soft atom density !
               !-------------------!
               CALL xc_dset_zero_all(deriv_set)
               IF (.NOT. evaluate_soft) THEN
                  exc_s = 0.0_dp
                  IF (.NOT. energy_only) THEN
                     vxc_s = 0.0_dp
                     IF (ASSOCIATED(vxg_s)) vxg_s = 0.0_dp
                     IF (ASSOCIATED(vtau_s)) vtau_s = 0.0_dp
                  END IF
               ELSE IF (skala_atom_grid) THEN
                  CALL skala_gapw_atom_vxc_of_r( &
                     my_xc_section, grid_atom, para_env, particle_set(iatom)%r, &
                     rho_s, drho_s, tau_s, weight_s, lsd, nspins, na, nr, &
                     exc_s, vxc_s, vxg_s, vtau_s, energy_only=energy_only, &
                     atom_force=skala_atom_force_s, atom_virial=skala_atom_virial_s)
               ELSE
                  CALL vxc_of_r_new(xc_fun_section, rho_set_s, deriv_set, 1, needs, weight_s, &
                                    lsd, na, nr, exc_s, vxc_s, vxg_s, vtau_s, energy_only=energy_only, &
                                    adiabatic_rescale_factor=my_adiabatic_rescale_factor)
               END IF
               rho_atom%exc_s = rho_atom%exc_s + exc_s

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

               CALL xc_dset_release(deriv_set)
               CALL xc_rho_set_release(rho_set_h)
               CALL xc_rho_set_release(rho_set_s)
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
!$OMP PRIVATE(base_shift, composite_row, cross_density, cross_density_spatial, &
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
                        DO idir = 1, 3
                           base_shift(idir) = cell%perd(idir)*NINT(fractional(idir))
                        END DO
                        DO image_i3 = base_shift(3) - image_shell(3), &
                           base_shift(3) + image_shell(3)
                           DO image_i2 = base_shift(2) - image_shell(2), &
                              base_shift(2) + image_shell(2)
                              DO image_i1 = base_shift(1) - image_shell(1), &
                                 base_shift(1) + image_shell(1)
                                 image_shift = [image_i1, image_i2, image_i3]
                                 IF (target_atom == source_atom .AND. &
                                     ALL(image_shift == 0)) CYCLE
                                 image_translation = MATMUL( &
                                                     cell%hmat, REAL(image_shift, dp))
                                 cross_displacement = composite_grid_coords(:, composite_row) - &
                                                      particle_set(source_atom)%r - &
                                                      image_translation
                                 CALL interpolate_gapw_atom_grid_fields( &
                                    grid_atom, harmonics, cross_displacement, cross_cutoff, nspins, &
                                    rho_h, rho_s, drho_h, drho_s, tau_h, tau_s, &
                                    cross_density, cross_grad, cross_kin, cross_density_spatial, &
                                    cross_grad_spatial, cross_kin_spatial)
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
            atom_composite_nelec = SUM(composite_grid_weights* &
                                       (composite_density(:, 1) + composite_density(:, 2)))
            CALL para_env%sum(atom_composite_nelec)
            IF (atom_composite_reference .AND. my_calculate_forces) THEN
               CALL skala_gapw_atom_composite_energy( &
                  my_xc_section, para_env, composite_density, composite_grad, composite_kin, &
                  composite_grid_coords, composite_grid_weights, composite_atomic_grid_weights, &
                  composite_local_grid_sizes, composite_local_atom_coords, atom_composite_exc, &
                  composite_density_grad, composite_grad_grad, composite_kin_grad, &
                  composite_grid_coord_grad, composite_grid_weight_grad, &
                  composite_atomic_grid_weight_grad, composite_atom_coord_grad)
               ALLOCATE (composite_cross_force(3, SIZE(particle_set)), &
                         composite_explicit_force(3, SIZE(particle_set)), &
                         composite_grid_coord_force(3, SIZE(particle_set)), &
                         composite_model_atom_force(3, SIZE(particle_set)), &
                         composite_moving_smooth_force(3, SIZE(particle_set)), &
                         composite_nlcc_center_force(3, SIZE(particle_set)), &
                         composite_nlcc_target_force(3, SIZE(particle_set)), &
                         composite_partition_force(3, SIZE(particle_set)), &
                         composite_partition_included(SIZE(particle_set)), &
                         composite_partition_datom(3, SIZE(particle_set), SIZE(particle_set)), &
                         composite_partition_dstrain(3, 3, SIZE(particle_set)))
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
!$OMP PARALLEL DO IF (image_partition_atom_composite) SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP PRIVATE(composite_local_atom, composite_row, composite_smooth_gradient_value, &
!$OMP         iatom, idir, ispin, jdir) &
!$OMP SHARED(cell, composite_atom_end, composite_atom_start, composite_grad_grad, &
!$OMP        composite_grid_coords, composite_local_atoms, composite_local_natom, &
!$OMP        composite_smooth_drho, composite_smooth_drhoa, composite_smooth_drhob, &
!$OMP        image_partition_atom_composite, lsd, smooth_rho_r) &
!$OMP REDUCTION(+:composite_feature_virial)
               DO composite_local_atom = 1, composite_local_natom
                  iatom = composite_local_atoms(composite_local_atom)
                  DO composite_row = composite_atom_start(iatom), composite_atom_end(iatom)
                     IF (lsd) THEN
                        DO jdir = 1, 3
                           composite_smooth_gradient_value(jdir, 1) = &
                              interpolate_native_grid( &
                              composite_smooth_drhoa(jdir)%array, smooth_rho_r(1)%pw_grid, &
                              cell, composite_grid_coords(:, composite_row), &
                              image_partition_atom_composite)
                           composite_smooth_gradient_value(jdir, 2) = &
                              interpolate_native_grid( &
                              composite_smooth_drhob(jdir)%array, smooth_rho_r(1)%pw_grid, &
                              cell, composite_grid_coords(:, composite_row), &
                              image_partition_atom_composite)
                        END DO
                     ELSE
                        DO jdir = 1, 3
                           composite_smooth_gradient_value(jdir, :) = 0.5_dp* &
                                                                      interpolate_native_grid( &
                                                                      composite_smooth_drho(jdir)%array, smooth_rho_r(1)%pw_grid, &
                                                                      cell, composite_grid_coords(:, composite_row), &
                                                                      image_partition_atom_composite)
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
!$OMP PARALLEL IF (image_partition_atom_composite) DEFAULT(NONE) &
!$OMP PRIVATE(composite_local_atom, composite_row, descriptor_window_adjoint, &
!$OMP         descriptor_window_weight, gth_potential, iatom, idir, ispin, jdir, &
!$OMP         nlcc_density, nlcc_gradient, nlcc_hessian, nlcc_spatial_derivative, &
!$OMP         nlcc_spin_factor, partition_adjoint, composite_nlcc_center_force_local, &
!$OMP         composite_partition_force_local, local_descriptor_datom, &
!$OMP         local_descriptor_dstrain, local_partition_datom, local_partition_dstrain, &
!$OMP         local_partition_weight, sgp_potential, source_atom, &
!$OMP         spatial_derivative, target_atom, target_partition_adjoint) &
!$OMP REDUCTION(+:composite_interpolation_virial, composite_partition_strain_virial) &
!$OMP SHARED(cell, composite_atom_coord_grad, composite_atom_coords, composite_atom_end, &
!$OMP        composite_atom_kind, composite_atom_start, composite_atomic_grid_weight_grad, &
!$OMP        composite_atomic_grid_weights, composite_base_grid_weights, composite_density_grad, &
!$OMP        composite_grad_grad, composite_grid_coord_force, composite_grid_coord_grad, &
!$OMP        composite_grid_coords, composite_grid_weight_grad, composite_image_periodicity, &
!$OMP        composite_kin_grad, composite_local_atoms, composite_local_natom, &
!$OMP        composite_model_atom_force, composite_moving_smooth_force, &
!$OMP        composite_nlcc_center_force, composite_nlcc_target_force, composite_partition_datom, &
!$OMP        composite_partition_dstrain, composite_partition_force, composite_partition_included, &
!$OMP        composite_partition_weights, composite_smooth_drho, composite_smooth_drhoa, &
!$OMP        composite_smooth_drhob, composite_smooth_rho, composite_smooth_rhoa, &
!$OMP        composite_smooth_rhob, composite_smooth_tau, composite_smooth_tau_a, &
!$OMP        composite_smooth_tau_b, image_partition_atom_composite, lsd, my_kind_set, nlcc, &
!$OMP        particle_set, smooth_rho_r)
               ALLOCATE (composite_nlcc_center_force_local(3, SIZE(particle_set)), &
                         composite_partition_force_local(3, SIZE(particle_set)), &
                         local_partition_datom(3, SIZE(particle_set)))
               composite_nlcc_center_force_local = 0.0_dp
               composite_partition_force_local = 0.0_dp
!$OMP DO SCHEDULE(DYNAMIC)
               DO composite_local_atom = 1, composite_local_natom
                  iatom = composite_local_atoms(composite_local_atom)
                  composite_model_atom_force(:, iatom) = &
                     composite_atom_coord_grad(:, composite_local_atom)
                  DO composite_row = composite_atom_start(iatom), composite_atom_end(iatom)
                     composite_grid_coord_force(:, iatom) = &
                        composite_grid_coord_force(:, iatom) + &
                        composite_grid_coord_grad(:, composite_row)
                     IF (lsd) THEN
                        spatial_derivative = &
                           composite_density_grad(composite_row, 1)* &
                           interpolate_native_grid_gradient( &
                           composite_smooth_rhoa, smooth_rho_r(1)%pw_grid, cell, &
                           composite_grid_coords(:, composite_row), &
                           image_partition_atom_composite) + &
                           composite_density_grad(composite_row, 2)* &
                           interpolate_native_grid_gradient( &
                           composite_smooth_rhob, smooth_rho_r(1)%pw_grid, cell, &
                           composite_grid_coords(:, composite_row), &
                           image_partition_atom_composite) + &
                           composite_kin_grad(composite_row, 1)* &
                           interpolate_native_grid_gradient( &
                           composite_smooth_tau_a, smooth_rho_r(1)%pw_grid, cell, &
                           composite_grid_coords(:, composite_row), &
                           image_partition_atom_composite) + &
                           composite_kin_grad(composite_row, 2)* &
                           interpolate_native_grid_gradient( &
                           composite_smooth_tau_b, smooth_rho_r(1)%pw_grid, cell, &
                           composite_grid_coords(:, composite_row), &
                           image_partition_atom_composite)
                        DO idir = 1, 3
                           spatial_derivative = spatial_derivative + &
                                                composite_grad_grad(composite_row, idir, 1)* &
                                                interpolate_native_grid_gradient( &
                                                composite_smooth_drhoa(idir)%array, &
                                                smooth_rho_r(1)%pw_grid, cell, &
                                                composite_grid_coords(:, composite_row), &
                                                image_partition_atom_composite) + &
                                                composite_grad_grad(composite_row, idir, 2)* &
                                                interpolate_native_grid_gradient( &
                                                composite_smooth_drhob(idir)%array, &
                                                smooth_rho_r(1)%pw_grid, cell, &
                                                composite_grid_coords(:, composite_row), &
                                                image_partition_atom_composite)
                        END DO
                     ELSE
                        spatial_derivative = 0.5_dp*SUM( &
                                             composite_density_grad(composite_row, :))* &
                                             interpolate_native_grid_gradient( &
                                             composite_smooth_rho, smooth_rho_r(1)%pw_grid, cell, &
                                             composite_grid_coords(:, composite_row), &
                                             image_partition_atom_composite) + &
                                             0.5_dp*SUM(composite_kin_grad(composite_row, :))* &
                                             interpolate_native_grid_gradient( &
                                             composite_smooth_tau, smooth_rho_r(1)%pw_grid, cell, &
                                             composite_grid_coords(:, composite_row), &
                                             image_partition_atom_composite)
                        DO idir = 1, 3
                           spatial_derivative = spatial_derivative + 0.5_dp*SUM( &
                                                composite_grad_grad(composite_row, idir, :))* &
                                                interpolate_native_grid_gradient( &
                                                composite_smooth_drho(idir)%array, &
                                                smooth_rho_r(1)%pw_grid, &
                                                cell, composite_grid_coords(:, composite_row), &
                                                image_partition_atom_composite)
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
                     composite_moving_smooth_force(:, iatom) = &
                        composite_moving_smooth_force(:, iatom) + spatial_derivative
                     IF (nlcc) THEN
                        nlcc_spin_factor = MERGE(1.0_dp, 0.5_dp, lsd)
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
                           composite_nlcc_target_force(:, iatom) = &
                              composite_nlcc_target_force(:, iatom) + nlcc_spatial_derivative
                           composite_nlcc_center_force_local(:, source_atom) = &
                              composite_nlcc_center_force_local(:, source_atom) - nlcc_spatial_derivative
                        END DO
                     END IF
                     IF (image_partition_atom_composite) THEN
                        CALL periodic_atom_image_partition( &
                           composite_grid_coords(:, composite_row), composite_atom_coords, cell, &
                           iatom, local_partition_weight, local_partition_datom, &
                           local_partition_dstrain, composite_image_periodicity)
                        CALL periodic_atom_image_partition( &
                           composite_grid_coords(:, composite_row), &
                           composite_atom_coords(:, iatom:iatom), cell, 1, &
                           descriptor_window_weight, local_descriptor_datom, &
                           local_descriptor_dstrain, composite_image_periodicity)
                        ! Its atom derivative cancels against the moving target grid; periodic
                        ! image strain remains an explicit contribution to the virial.
                        target_partition_adjoint = composite_base_grid_weights(composite_row)* &
                                                   composite_grid_weight_grad(composite_row)
                        descriptor_window_adjoint = composite_base_grid_weights(composite_row)* &
                                                    composite_atomic_grid_weight_grad(composite_row)* &
                                                    smooth_partition_atomic_weight_scale_derivative( &
                                                    descriptor_window_weight)
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
                           local_partition_dstrain - descriptor_window_adjoint* &
                           local_descriptor_dstrain
                     ELSE
                        CALL skala_gpw_smooth_partition_derivatives( &
                           composite_grid_coords(:, composite_row), composite_atom_coords, cell, &
                           composite_partition_weights, composite_partition_included, &
                           composite_partition_datom, composite_partition_dstrain)
                        partition_adjoint = composite_grid_weight_grad(composite_row)* &
                                            composite_atomic_grid_weights(composite_row)
                        DO target_atom = 1, SIZE(particle_set)
                           composite_partition_force_local(:, target_atom) = &
                              composite_partition_force_local(:, target_atom) + &
                              partition_adjoint* &
                              composite_partition_datom(:, target_atom, iatom)
                        END DO
                        composite_partition_force_local(:, iatom) = &
                           composite_partition_force_local(:, iatom) - partition_adjoint* &
                           SUM(composite_partition_datom(:, :, iatom), DIM=2)
                     END IF
                  END DO
               END DO
!$OMP END DO
!$OMP CRITICAL(skala_atom_composite_force_reduction)
               composite_nlcc_center_force(:, :) = composite_nlcc_center_force(:, :) + &
                                                   composite_nlcc_center_force_local
               composite_partition_force(:, :) = composite_partition_force(:, :) + &
                                                 composite_partition_force_local
!$OMP END CRITICAL(skala_atom_composite_force_reduction)
               DEALLOCATE (composite_nlcc_center_force_local, composite_partition_force_local, &
                           local_partition_datom)
!$OMP END PARALLEL
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
!$OMP PRIVATE(base_shift, composite_cross_force_local, composite_cross_image_virial_local, &
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
                           DO idir = 1, 3
                              base_shift(idir) = cell%perd(idir)*NINT(fractional(idir))
                           END DO
                           DO image_i3 = base_shift(3) - image_shell(3), &
                              base_shift(3) + image_shell(3)
                              DO image_i2 = base_shift(2) - image_shell(2), &
                                 base_shift(2) + image_shell(2)
                                 DO image_i1 = base_shift(1) - image_shell(1), &
                                    base_shift(1) + image_shell(1)
                                    image_shift = [image_i1, image_i2, image_i3]
                                    IF (target_atom == source_atom .AND. &
                                        ALL(image_shift == 0)) CYCLE
                                    image_translation = MATMUL( &
                                                        cell%hmat, REAL(image_shift, dp))
                                    cross_displacement = &
                                       composite_grid_coords(:, composite_row) - &
                                       particle_set(source_atom)%r - image_translation
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
                  virial%pv_xc = composite_feature_virial - composite_interpolation_virial
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
                           composite_partition_datom, composite_partition_dstrain, &
                           composite_partition_force, composite_partition_included)
            ELSE
               CALL skala_gapw_atom_composite_energy( &
                  my_xc_section, para_env, composite_density, composite_grad, composite_kin, &
                  composite_grid_coords, composite_grid_weights, composite_atomic_grid_weights, &
                  composite_local_grid_sizes, composite_local_atom_coords, atom_composite_exc, &
                  composite_density_grad, composite_grad_grad, composite_kin_grad)
            END IF
            composite_pw_nflat = PRODUCT(smooth_rho_r(1)%pw_grid%npts)
            adjoint_nchannels = MERGE(2, 1, lsd)
            ALLOCATE (smooth_density_adjoint_storage(composite_pw_nflat, adjoint_nchannels), &
                      smooth_grad_adjoint_storage(composite_pw_nflat, 3, adjoint_nchannels), &
                      smooth_kin_adjoint_storage(composite_pw_nflat, adjoint_nchannels))
            smooth_density_adjoint_storage = 0.0_dp
            smooth_grad_adjoint_storage = 0.0_dp
            smooth_kin_adjoint_storage = 0.0_dp
            CALL build_native_grid_adjoint_bins( &
               smooth_rho_r(1)%pw_grid, cell, composite_grid_coords, &
               image_partition_atom_composite, adjoint_tile_count, adjoint_bin_offsets, &
               adjoint_bin_rows)
            adjoint_nbins = SIZE(adjoint_bin_offsets) - 1
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(NONE) &
!$OMP PRIVATE(adjoint_entry, composite_row, composite_smooth_density_adjoint_value, &
!$OMP         composite_smooth_gradient_adjoint_value, composite_smooth_kin_adjoint_value, &
!$OMP         adjoint_tile_lower, adjoint_tile_upper, interpolation_stencil) &
!$OMP SHARED(adjoint_bin_offsets, adjoint_bin_rows, adjoint_nbins, adjoint_nchannels, &
!$OMP        adjoint_tile_count, cell, composite_density_grad, composite_grad_grad, &
!$OMP        composite_grid_coords, composite_kin_grad, image_partition_atom_composite, lsd, &
!$OMP        smooth_density_adjoint_storage, smooth_grad_adjoint_storage, &
!$OMP        smooth_kin_adjoint_storage, smooth_rho_r)
            DO adjoint_bin = 1, adjoint_nbins
               CALL native_grid_adjoint_tile_bounds( &
                  smooth_rho_r(1)%pw_grid, adjoint_bin, adjoint_tile_count, &
                  adjoint_tile_lower, adjoint_tile_upper)
               DO adjoint_entry = adjoint_bin_offsets(adjoint_bin), &
                  adjoint_bin_offsets(adjoint_bin + 1) - 1
                  composite_row = adjoint_bin_rows(adjoint_entry)
                  CALL create_native_grid_interpolation_stencil( &
                     interpolation_stencil, smooth_rho_r(1)%pw_grid, cell, &
                     composite_grid_coords(:, composite_row), image_partition_atom_composite)
                  IF (lsd) THEN
                     composite_smooth_density_adjoint_value = &
                        composite_density_grad(composite_row, :)
                     composite_smooth_gradient_adjoint_value = &
                        composite_grad_grad(composite_row, :, :)
                     composite_smooth_kin_adjoint_value = composite_kin_grad(composite_row, :)
                  ELSE
                     composite_smooth_density_adjoint_value = 0.0_dp
                     composite_smooth_gradient_adjoint_value = 0.0_dp
                     composite_smooth_kin_adjoint_value = 0.0_dp
                     composite_smooth_density_adjoint_value(1) = &
                        SUM(composite_density_grad(composite_row, :))
                     composite_smooth_gradient_adjoint_value(:, 1) = &
                        SUM(composite_grad_grad(composite_row, :, :), DIM=2)
                     composite_smooth_kin_adjoint_value(1) = &
                        SUM(composite_kin_grad(composite_row, :))
                  END IF
                  CALL add_native_grid_fields_adjoint_tile( &
                     smooth_density_adjoint_storage, smooth_grad_adjoint_storage, &
                     smooth_kin_adjoint_storage, smooth_rho_r(1)%pw_grid, interpolation_stencil, &
                     composite_smooth_density_adjoint_value, &
                     composite_smooth_gradient_adjoint_value, &
                     composite_smooth_kin_adjoint_value, adjoint_nchannels, &
                     adjoint_tile_lower, adjoint_tile_upper)
               END DO
            END DO
!$OMP END PARALLEL DO
            DEALLOCATE (adjoint_bin_offsets, adjoint_bin_rows)
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
            smooth_density_adjoint => smooth_density_adjoint_storage
            smooth_grad_adjoint => smooth_grad_adjoint_storage
            smooth_kin_adjoint => smooth_kin_adjoint_storage
            ! The interpolation transpose still uses CP2K's global FFT-grid layout. Reduce only
            ! this PW adjoint; atom-grid feature rows and their model derivatives stay rank-local.
            CALL para_env%sum(smooth_density_adjoint)
            CALL para_env%sum(smooth_grad_adjoint)
            CALL para_env%sum(smooth_kin_adjoint)
            CALL para_env%sum(smooth_input_contraction)
            CALL build_vxc_from_feature_grads( &
               smooth_vxc_rho, smooth_vxc_tau, smooth_rho_r, auxbas_pw_pool, &
               smooth_density_adjoint, smooth_grad_adjoint, smooth_kin_adjoint, &
               xc_deriv_method_id, global_grid_layout=.TRUE.)
            NULLIFY (smooth_density_adjoint, smooth_grad_adjoint, smooth_kin_adjoint)
            DEALLOCATE (smooth_density_adjoint_storage, smooth_grad_adjoint_storage, &
                        smooth_kin_adjoint_storage)
            smooth_grid_contraction = 0.0_dp
            DO ispin = 1, nspins
               smooth_grid_contraction = smooth_grid_contraction + smooth_rho_r(1)%pw_grid%dvol* &
                                         (SUM(smooth_vxc_rho(ispin)%array* &
                                              smooth_rho_r(ispin)%array) + &
                                          SUM(smooth_vxc_tau(ispin)%array* &
                                              smooth_tau_r(ispin)%array))
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
            CALL para_env%sum(smooth_grid_contraction)
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
                        DO composite_local_atom = 1, composite_local_natom
                           target_atom = composite_local_atoms(composite_local_atom)
                           DO composite_row = composite_atom_start(target_atom), &
                              composite_atom_end(target_atom)
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
                              DO idir = 1, 3
                                 base_shift(idir) = cell%perd(idir)*NINT(fractional(idir))
                              END DO
                              DO image_i3 = base_shift(3) - image_shell(3), &
                                 base_shift(3) + image_shell(3)
                                 DO image_i2 = base_shift(2) - image_shell(2), &
                                    base_shift(2) + image_shell(2)
                                    DO image_i1 = base_shift(1) - image_shell(1), &
                                       base_shift(1) + image_shell(1)
                                       image_shift = [image_i1, image_i2, image_i3]
                                       IF (target_atom == iatom .AND. ALL(image_shift == 0)) CYCLE
                                       image_translation = MATMUL( &
                                                           cell%hmat, REAL(image_shift, dp))
                                       cross_displacement = &
                                          composite_grid_coords(:, composite_row) - &
                                          particle_set(iatom)%r - image_translation
                                       CALL add_gapw_atom_grid_interpolation_adjoint( &
                                          grid_atom, harmonics, cross_displacement, cross_cutoff, &
                                          nspins, cross_density_adjoint, cross_grad_adjoint, &
                                          cross_kin_adjoint, vxc_h, vxc_s, vxg_h, vxg_s, &
                                          vtau_h, vtau_s)
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END IF

                     ! Model rows are distributed by target atom, while CP2K stores each
                     ! one-center matrix on the rank owning its source atom. Sum the exact
                     ! interpolation transpose before forming that matrix.
                     CALL para_env%sum(vxc_h)
                     CALL para_env%sum(vxc_s)
                     CALL para_env%sum(vxg_h)
                     CALL para_env%sum(vxg_s)
                     CALL para_env%sum(vtau_h)
                     CALL para_env%sum(vtau_s)
                     IF (.NOT. source_matrix_local) CYCLE

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
                        one_center_density_matrix_contraction = &
                           one_center_density_matrix_contraction + &
                           contract_one_center_matrix(cpc_h(ispin)%r_coef, &
                                                      int_hh(ispin)%r_coef, &
                                                      tau_basis_cache%n2oindex) - &
                           contract_one_center_matrix(cpc_s(ispin)%r_coef, &
                                                      int_ss(ispin)%r_coef, &
                                                      tau_basis_cache%n2oindex)
                        int_hh(ispin)%r_coef = 0.0_dp
                        int_ss(ispin)%r_coef = 0.0_dp
                     END DO
                     CALL gaVxcgb_GC(vxc_h, vxc_s, vxg_h, vxg_s, int_hh, int_ss, &
                                     grid_atom, basis_1c, harmonics, nspins)
                     DO ispin = 1, nspins
                        one_center_rho_grad_matrix_contraction = &
                           one_center_rho_grad_matrix_contraction + &
                           contract_one_center_matrix(cpc_h(ispin)%r_coef, &
                                                      int_hh(ispin)%r_coef, &
                                                      tau_basis_cache%n2oindex) - &
                           contract_one_center_matrix(cpc_s(ispin)%r_coef, &
                                                      int_ss(ispin)%r_coef, &
                                                      tau_basis_cache%n2oindex)
                     END DO
                     CALL dgaVtaudgb(vtau_h, vtau_s, int_hh, int_ss, tau_basis_cache, nspins)
                     DO ispin = 1, nspins
                        one_center_matrix_contraction = one_center_matrix_contraction + &
                                                        contract_one_center_matrix(cpc_h(ispin)%r_coef, &
                                                                                   int_hh(ispin)%r_coef, &
                                                                                   tau_basis_cache%n2oindex) - &
                                                        contract_one_center_matrix(cpc_s(ispin)%r_coef, &
                                                                                   int_ss(ispin)%r_coef, &
                                                                                   tau_basis_cache%n2oindex)
                        IF (.NOT. atom_composite_reference) THEN
                           int_hh(ispin)%r_coef = composite_int_h(:, :, ispin)
                           int_ss(ispin)%r_coef = composite_int_s(:, :, ispin)
                        END IF
                     END DO
                     DEALLOCATE (composite_int_h, composite_int_s)
                  END DO
                  CALL release_tau_basis_cache(tau_basis_cache)
               END DO
            END IF
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
            IF (atom_composite_diagnostic) THEN
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
                     composite_grid_coords, composite_grid_weights, composite_atomic_grid_weights, &
                     composite_local_grid_sizes, composite_local_atom_coords, feature_vxc_plus)
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
                     composite_grid_coords, composite_grid_weights, composite_atomic_grid_weights, &
                     composite_local_grid_sizes, composite_local_atom_coords, feature_vxc_minus)
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
                        composite_partition_weights, composite_partition_atom_coords, &
                        composite_distances, composite_density, composite_grad, composite_kin, &
                        composite_smooth_density_cache, composite_smooth_gradient_cache, &
                        composite_smooth_kin_cache, &
                        composite_density_grad, composite_grad_grad, composite_kin_grad, &
                        composite_grid_coords, composite_grid_weights, &
                        composite_base_grid_weights, composite_atomic_grid_weights)
            IF (lsd) THEN
               DEALLOCATE (composite_smooth_rhoa, composite_smooth_rhob, &
                           composite_smooth_tau_a, composite_smooth_tau_b)
               DO idir = 1, 3
                  DEALLOCATE (composite_smooth_drhoa(idir)%array, &
                              composite_smooth_drhob(idir)%array)
               END DO
            ELSE
               DEALLOCATE (composite_smooth_rho, composite_smooth_tau)
               DO idir = 1, 3
                  DEALLOCATE (composite_smooth_drho(idir)%array)
               END DO
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
!> \brief Replicate a distributed real-space field for atom-grid interpolation.
!> \param local_values ...
!> \param pw_grid ...
!> \param group ...
!> \param global_values ...
! **************************************************************************************************
   SUBROUTINE gather_native_grid_field(local_values, pw_grid, group, global_values)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         POINTER                                         :: local_values
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(mp_para_env_type), INTENT(IN)                 :: group
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT), &
         POINTER                                         :: global_values

      INTEGER, DIMENSION(2, 3)                           :: bo

      CPASSERT(ASSOCIATED(local_values))
      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(.NOT. ASSOCIATED(global_values))
      bo = pw_grid%bounds_local
      ALLOCATE (global_values(pw_grid%bounds(1, 1):pw_grid%bounds(2, 1), &
                              pw_grid%bounds(1, 2):pw_grid%bounds(2, 2), &
                              pw_grid%bounds(1, 3):pw_grid%bounds(2, 3)))
      global_values = 0.0_dp
      global_values(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)) = &
         local_values
      CALL group%sum(global_values)

   END SUBROUTINE gather_native_grid_field

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
            CALL native_grid_lagrange_weights(FRACTION(idir), stencil%weight(:, idir))
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
!> \brief Interpolate density, gradient, and kinetic-density fields in one stencil traversal.
!> \param density ...
!> \param grad_x ...
!> \param grad_y ...
!> \param grad_z ...
!> \param kin ...
!> \param stencil ...
!> \param density_value ...
!> \param grad_value ...
!> \param kin_value ...
! **************************************************************************************************
   SUBROUTINE interpolate_native_grid_fields(density, grad_x, grad_y, grad_z, kin, stencil, &
                                             density_value, grad_value, kin_value)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: density, grad_x, grad_y, grad_z, kin
      TYPE(native_grid_interpolation_stencil_type), &
         INTENT(IN)                                      :: stencil
      REAL(KIND=dp), INTENT(OUT)                         :: density_value
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: grad_value
      REAL(KIND=dp), INTENT(OUT)                         :: kin_value

      INTEGER                                            :: inode_x, inode_y, inode_z
      INTEGER, DIMENSION(3)                              :: index, lower_bound
      REAL(KIND=dp)                                      :: coefficient

      CPASSERT(ALL(SHAPE(grad_x) == SHAPE(density)))
      CPASSERT(ALL(SHAPE(grad_y) == SHAPE(density)))
      CPASSERT(ALL(SHAPE(grad_z) == SHAPE(density)))
      CPASSERT(ALL(SHAPE(kin) == SHAPE(density)))
      density_value = 0.0_dp
      grad_value = 0.0_dp
      kin_value = 0.0_dp
      IF (.NOT. stencil%active) RETURN
      lower_bound = LBOUND(density)
      DO inode_z = 1, native_grid_interp_npts
         IF (.NOT. stencil%valid(inode_z, 3)) CYCLE
         INDEX(3) = lower_bound(3) + stencil%relative_index(inode_z, 3)
         DO inode_y = 1, native_grid_interp_npts
            IF (.NOT. stencil%valid(inode_y, 2)) CYCLE
            INDEX(2) = lower_bound(2) + stencil%relative_index(inode_y, 2)
            DO inode_x = 1, native_grid_interp_npts
               IF (.NOT. stencil%valid(inode_x, 1)) CYCLE
               INDEX(1) = lower_bound(1) + stencil%relative_index(inode_x, 1)
               coefficient = stencil%weight(inode_x, 1)* &
                             stencil%weight(inode_y, 2)* &
                             stencil%weight(inode_z, 3)
               density_value = density_value + coefficient*density(INDEX(1), INDEX(2), INDEX(3))
               grad_value(1) = grad_value(1) + coefficient*grad_x(INDEX(1), INDEX(2), INDEX(3))
               grad_value(2) = grad_value(2) + coefficient*grad_y(INDEX(1), INDEX(2), INDEX(3))
               grad_value(3) = grad_value(3) + coefficient*grad_z(INDEX(1), INDEX(2), INDEX(3))
               kin_value = kin_value + coefficient*kin(INDEX(1), INDEX(2), INDEX(3))
            END DO
         END DO
      END DO

   END SUBROUTINE interpolate_native_grid_fields

! **************************************************************************************************
!> \brief Group atom-grid rows by the disjoint PW tiles touched by their interpolation stencils.
!> \param pw_grid ...
!> \param cell ...
!> \param points atom-grid coordinates
!> \param wrap_auxiliary_cell wrap all auxiliary-grid directions
!> \param tile_count number of tiles along each PW-grid direction
!> \param bin_offsets CSR offsets into bin_rows
!> \param bin_rows atom-grid rows touching each tile
! **************************************************************************************************
   SUBROUTINE build_native_grid_adjoint_bins(pw_grid, cell, points, wrap_auxiliary_cell, &
                                             tile_count, bin_offsets, bin_rows)
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: points
      LOGICAL, INTENT(IN)                                :: wrap_auxiliary_cell
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: tile_count
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: bin_offsets, bin_rows

      INTEGER                                            :: candidate, ibin, idir, inode, irow, ix, &
                                                            iy, iz, nbin, ntouched, position
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bin_counts, bin_cursor, row_bin_count
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: row_bins
      INTEGER, DIMENSION(3)                              :: touched_count
      INTEGER, DIMENSION(&
         native_grid_adjoint_max_tiles_per_direction, 3) :: touched_tiles
      TYPE(native_grid_interpolation_stencil_type)       :: stencil

      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(ASSOCIATED(cell))
      CPASSERT(SIZE(points, 1) == 3)
      tile_count = (pw_grid%npts + native_grid_adjoint_tile_edge - 1)/ &
                   native_grid_adjoint_tile_edge
      nbin = PRODUCT(tile_count)
      ALLOCATE (row_bin_count(SIZE(points, 2)), &
                row_bins(native_grid_adjoint_max_bins_per_row, SIZE(points, 2)))
      row_bin_count = 0
      row_bins = 0
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP PRIVATE(candidate, idir, inode, ix, iy, iz, ntouched, stencil, &
!$OMP         touched_count, touched_tiles) &
!$OMP SHARED(cell, points, pw_grid, row_bin_count, row_bins, tile_count, wrap_auxiliary_cell)
      DO irow = 1, SIZE(points, 2)
         CALL create_native_grid_interpolation_stencil( &
            stencil, pw_grid, cell, points(:, irow), wrap_auxiliary_cell, indices_only=.TRUE.)
         IF (.NOT. stencil%active) CYCLE
         touched_count = 0
         touched_tiles = 0
         DO idir = 1, 3
            DO inode = 1, native_grid_interp_npts
               IF (.NOT. stencil%valid(inode, idir)) CYCLE
               candidate = stencil%relative_index(inode, idir)/native_grid_adjoint_tile_edge + 1
               IF (ALL(touched_tiles(1:touched_count(idir), idir) /= candidate)) THEN
                  touched_count(idir) = touched_count(idir) + 1
                  CPASSERT(touched_count(idir) <= SIZE(touched_tiles, 1))
                  touched_tiles(touched_count(idir), idir) = candidate
               END IF
            END DO
         END DO
         ntouched = 0
         DO iz = 1, touched_count(3)
            DO iy = 1, touched_count(2)
               DO ix = 1, touched_count(1)
                  ntouched = ntouched + 1
                  CPASSERT(ntouched <= native_grid_adjoint_max_bins_per_row)
                  row_bins(ntouched, irow) = 1 + touched_tiles(ix, 1) - 1 + tile_count(1)*( &
                                             touched_tiles(iy, 2) - 1 + tile_count(2)*( &
                                             touched_tiles(iz, 3) - 1))
               END DO
            END DO
         END DO
         row_bin_count(irow) = ntouched
      END DO
!$OMP END PARALLEL DO

      ALLOCATE (bin_counts(nbin), bin_cursor(nbin), bin_offsets(nbin + 1))
      bin_counts = 0
      DO irow = 1, SIZE(points, 2)
         DO ibin = 1, row_bin_count(irow)
            bin_counts(row_bins(ibin, irow)) = bin_counts(row_bins(ibin, irow)) + 1
         END DO
      END DO
      bin_offsets(1) = 1
      DO ibin = 1, nbin
         bin_offsets(ibin + 1) = bin_offsets(ibin) + bin_counts(ibin)
      END DO
      ALLOCATE (bin_rows(bin_offsets(nbin + 1) - 1))
      bin_cursor(:) = bin_offsets(1:nbin)
      DO irow = 1, SIZE(points, 2)
         DO ibin = 1, row_bin_count(irow)
            candidate = row_bins(ibin, irow)
            position = bin_cursor(candidate)
            bin_rows(position) = irow
            bin_cursor(candidate) = position + 1
         END DO
      END DO
      DEALLOCATE (bin_counts, bin_cursor, row_bin_count, row_bins)

   END SUBROUTINE build_native_grid_adjoint_bins

! **************************************************************************************************
!> \brief Return the inclusive PW-grid bounds owned by one linear tile index.
!> \param pw_grid ...
!> \param tile_index linear tile index
!> \param tile_count number of tiles along each PW-grid direction
!> \param tile_lower zero-based lower grid index
!> \param tile_upper zero-based upper grid index
! **************************************************************************************************
   SUBROUTINE native_grid_adjoint_tile_bounds(pw_grid, tile_index, tile_count, &
                                              tile_lower, tile_upper)
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      INTEGER, INTENT(IN)                                :: tile_index
      INTEGER, DIMENSION(3), INTENT(IN)                  :: tile_count
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: tile_lower, tile_upper

      INTEGER                                            :: linear_tile
      INTEGER, DIMENSION(3)                              :: expected_tile_count, tile_coord

      CPASSERT(ASSOCIATED(pw_grid))
      expected_tile_count = (pw_grid%npts + native_grid_adjoint_tile_edge - 1)/ &
                            native_grid_adjoint_tile_edge
      CPASSERT(ALL(tile_count == expected_tile_count))
      CPASSERT(tile_index >= 1 .AND. tile_index <= PRODUCT(tile_count))
      linear_tile = tile_index - 1
      tile_coord(1) = MOD(linear_tile, tile_count(1))
      linear_tile = linear_tile/tile_count(1)
      tile_coord(2) = MOD(linear_tile, tile_count(2))
      tile_coord(3) = linear_tile/tile_count(2)
      tile_lower = tile_coord*native_grid_adjoint_tile_edge
      tile_upper = MIN(tile_lower + native_grid_adjoint_tile_edge - 1, pw_grid%npts - 1)

   END SUBROUTINE native_grid_adjoint_tile_bounds

! **************************************************************************************************
!> \brief Apply the primitive-field interpolation transpose inside one disjoint PW tile.
!> \param density ...
!> \param grad ...
!> \param kin ...
!> \param pw_grid ...
!> \param stencil ...
!> \param density_value ...
!> \param grad_value ...
!> \param kin_value ...
!> \param nchannels number of spin channels accumulated by the adjoint
!> \param tile_lower zero-based lower grid index owned by this call
!> \param tile_upper zero-based upper grid index owned by this call
! **************************************************************************************************
   SUBROUTINE add_native_grid_fields_adjoint_tile( &
      density, grad, kin, pw_grid, stencil, density_value, grad_value, kin_value, nchannels, &
      tile_lower, tile_upper)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: density
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: grad
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: kin
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(native_grid_interpolation_stencil_type), &
         INTENT(IN)                                      :: stencil
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: density_value
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grad_value
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: kin_value
      INTEGER, INTENT(IN)                                :: nchannels
      INTEGER, DIMENSION(3), INTENT(IN)                  :: tile_lower, tile_upper

      INTEGER                                            :: idir, inode_x, inode_y, inode_z, ipt, &
                                                            ispin
      INTEGER, DIMENSION(3)                              :: relative_index
      REAL(KIND=dp)                                      :: coefficient

      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(SIZE(density, 1) == PRODUCT(pw_grid%npts))
      CPASSERT(SIZE(density, 2) == nchannels)
      CPASSERT(SIZE(grad, 1) == PRODUCT(pw_grid%npts))
      CPASSERT(SIZE(grad, 2) == 3)
      CPASSERT(SIZE(grad, 3) == nchannels)
      CPASSERT(SIZE(kin, 1) == PRODUCT(pw_grid%npts))
      CPASSERT(SIZE(kin, 2) == nchannels)
      CPASSERT(SIZE(density_value) == 2)
      CPASSERT(ALL(SHAPE(grad_value) == [3, 2]))
      CPASSERT(SIZE(kin_value) == 2)
      CPASSERT(nchannels >= 1 .AND. nchannels <= 2)
      CPASSERT(ALL(tile_lower >= 0))
      CPASSERT(ALL(tile_upper >= tile_lower))
      CPASSERT(ALL(tile_upper < pw_grid%npts))
      IF (.NOT. stencil%active) RETURN
      DO inode_z = 1, native_grid_interp_npts
         IF (.NOT. stencil%valid(inode_z, 3)) CYCLE
         relative_index(3) = stencil%relative_index(inode_z, 3)
         IF (relative_index(3) < tile_lower(3) .OR. relative_index(3) > tile_upper(3)) CYCLE
         DO inode_y = 1, native_grid_interp_npts
            IF (.NOT. stencil%valid(inode_y, 2)) CYCLE
            relative_index(2) = stencil%relative_index(inode_y, 2)
            IF (relative_index(2) < tile_lower(2) .OR. relative_index(2) > tile_upper(2)) CYCLE
            DO inode_x = 1, native_grid_interp_npts
               IF (.NOT. stencil%valid(inode_x, 1)) CYCLE
               relative_index(1) = stencil%relative_index(inode_x, 1)
               IF (relative_index(1) < tile_lower(1) .OR. relative_index(1) > tile_upper(1)) CYCLE
               ipt = 1 + relative_index(1) + pw_grid%npts(1)*( &
                     relative_index(2) + pw_grid%npts(2)*relative_index(3))
               coefficient = stencil%weight(inode_x, 1)* &
                             stencil%weight(inode_y, 2)* &
                             stencil%weight(inode_z, 3)
               DO ispin = 1, nchannels
                  density(ipt, ispin) = density(ipt, ispin) + coefficient*density_value(ispin)
                  DO idir = 1, 3
                     grad(ipt, idir, ispin) = grad(ipt, idir, ispin) + &
                                              coefficient*grad_value(idir, ispin)
                  END DO
                  kin(ipt, ispin) = kin(ipt, ispin) + coefficient*kin_value(ispin)
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE add_native_grid_fields_adjoint_tile

! **************************************************************************************************
!> \brief Interpolate a replicated native-grid field at a Cartesian point.
!> \param values ...
!> \param pw_grid ...
!> \param cell ...
!> \param point ...
!> \param wrap_auxiliary_cell wrap all auxiliary-grid directions
!> \return ...
! **************************************************************************************************
   FUNCTION interpolate_native_grid(values, pw_grid, cell, point, wrap_auxiliary_cell) RESULT(value)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: values
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: point
      LOGICAL, INTENT(IN), OPTIONAL                      :: wrap_auxiliary_cell
      REAL(KIND=dp)                                      :: value

      INTEGER                                            :: corner_x, corner_y, corner_z, idir
      INTEGER, DIMENSION(3)                              :: base, index, relative_index
      LOGICAL                                            :: wrap_grid
      REAL(KIND=dp)                                      :: coefficient
      REAL(KIND=dp), DIMENSION(3)                        :: fraction, relative
      REAL(KIND=dp), &
         DIMENSION(native_grid_interp_npts, 3)           :: weights

      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(ASSOCIATED(cell))
      DO idir = 1, 3
         CPASSERT(SIZE(values, idir) == pw_grid%npts(idir))
      END DO

      wrap_grid = ANY(cell%perd == 0)
      IF (PRESENT(wrap_auxiliary_cell)) wrap_grid = wrap_auxiliary_cell
      relative = MATMUL(pw_grid%dh_inv, point)
      DO idir = 1, 3
         IF (cell%perd(idir) == 1 .OR. wrap_grid) THEN
            relative(idir) = MODULO(relative(idir), REAL(pw_grid%npts(idir), KIND=dp))
         ELSE IF (relative(idir) <= -REAL(native_grid_interp_offset_max, dp) .OR. &
                  relative(idir) >= REAL(pw_grid%npts(idir) - &
                                         native_grid_interp_offset_min, KIND=dp)) THEN
            value = 0.0_dp
            RETURN
         END IF
         base(idir) = FLOOR(relative(idir))
         FRACTION(idir) = relative(idir) - REAL(base(idir), KIND=dp)
         CALL native_grid_lagrange_weights(FRACTION(idir), weights(:, idir))
      END DO

      value = 0.0_dp
      DO corner_z = native_grid_interp_offset_min, native_grid_interp_offset_max
         DO corner_y = native_grid_interp_offset_min, native_grid_interp_offset_max
            DO corner_x = native_grid_interp_offset_min, native_grid_interp_offset_max
               relative_index = base + [corner_x, corner_y, corner_z]
               coefficient = weights(corner_x - native_grid_interp_offset_min + 1, 1)* &
                             weights(corner_y - native_grid_interp_offset_min + 1, 2)* &
                             weights(corner_z - native_grid_interp_offset_min + 1, 3)
               DO idir = 1, 3
                  IF (cell%perd(idir) == 1 .OR. wrap_grid) THEN
                     relative_index(idir) = MODULO(relative_index(idir), pw_grid%npts(idir))
                  ELSE IF (relative_index(idir) < 0 .OR. &
                           relative_index(idir) >= pw_grid%npts(idir)) THEN
                     coefficient = 0.0_dp
                  END IF
               END DO
               IF (coefficient == 0.0_dp) CYCLE
               index = LBOUND(values) + relative_index
               value = value + coefficient*values(INDEX(1), INDEX(2), INDEX(3))
            END DO
         END DO
      END DO

   END FUNCTION interpolate_native_grid

! **************************************************************************************************
!> \brief Return the Cartesian gradient of native-grid interpolation at one point.
!> \param values ...
!> \param pw_grid ...
!> \param cell ...
!> \param point ...
!> \param wrap_auxiliary_cell wrap all auxiliary-grid directions
!> \return ...
! **************************************************************************************************
   FUNCTION interpolate_native_grid_gradient(values, pw_grid, cell, point, &
                                             wrap_auxiliary_cell) RESULT(gradient)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: values
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: point
      LOGICAL, INTENT(IN), OPTIONAL                      :: wrap_auxiliary_cell
      REAL(KIND=dp), DIMENSION(3)                        :: gradient

      INTEGER                                            :: corner_x, corner_y, corner_z, idir, jdir
      INTEGER, DIMENSION(3)                              :: base, corner_index, index, relative_index
      LOGICAL                                            :: wrap_grid
      REAL(KIND=dp)                                      :: coefficient
      REAL(KIND=dp), DIMENSION(3)                        :: fraction, gradient_relative, relative
      REAL(KIND=dp), &
         DIMENSION(native_grid_interp_npts, 3)           :: derivative_weights, weights

      CPASSERT(ASSOCIATED(pw_grid))
      CPASSERT(ASSOCIATED(cell))
      DO idir = 1, 3
         CPASSERT(SIZE(values, idir) == pw_grid%npts(idir))
      END DO

      wrap_grid = ANY(cell%perd == 0)
      IF (PRESENT(wrap_auxiliary_cell)) wrap_grid = wrap_auxiliary_cell
      relative = MATMUL(pw_grid%dh_inv, point)
      DO idir = 1, 3
         IF (cell%perd(idir) == 1 .OR. wrap_grid) THEN
            relative(idir) = MODULO(relative(idir), REAL(pw_grid%npts(idir), KIND=dp))
         ELSE IF (relative(idir) <= -REAL(native_grid_interp_offset_max, dp) .OR. &
                  relative(idir) >= REAL(pw_grid%npts(idir) - &
                                         native_grid_interp_offset_min, KIND=dp)) THEN
            gradient = 0.0_dp
            RETURN
         END IF
         base(idir) = FLOOR(relative(idir))
         FRACTION(idir) = relative(idir) - REAL(base(idir), KIND=dp)
         CALL native_grid_lagrange_weights( &
            FRACTION(idir), weights(:, idir), derivative_weights(:, idir))
      END DO

      gradient_relative = 0.0_dp
      DO corner_z = native_grid_interp_offset_min, native_grid_interp_offset_max
         DO corner_y = native_grid_interp_offset_min, native_grid_interp_offset_max
            DO corner_x = native_grid_interp_offset_min, native_grid_interp_offset_max
               relative_index = base + [corner_x, corner_y, corner_z]
               corner_index = [corner_x, corner_y, corner_z] - &
                              native_grid_interp_offset_min + 1
               DO idir = 1, 3
                  IF (cell%perd(idir) == 1 .OR. wrap_grid) THEN
                     relative_index(idir) = MODULO(relative_index(idir), pw_grid%npts(idir))
                  ELSE IF (relative_index(idir) < 0 .OR. &
                           relative_index(idir) >= pw_grid%npts(idir)) THEN
                     EXIT
                  END IF
               END DO
               IF (idir <= 3) CYCLE
               index = LBOUND(values) + relative_index
               DO idir = 1, 3
                  coefficient = 1.0_dp
                  DO jdir = 1, 3
                     IF (jdir == idir) THEN
                        coefficient = coefficient*derivative_weights( &
                                      corner_index(jdir), jdir)
                     ELSE
                        coefficient = coefficient*weights( &
                                      corner_index(jdir), jdir)
                     END IF
                  END DO
                  gradient_relative(idir) = gradient_relative(idir) + &
                                            coefficient*values(INDEX(1), INDEX(2), INDEX(3))
               END DO
            END DO
         END DO
      END DO
      gradient = MATMUL(TRANSPOSE(pw_grid%dh_inv), gradient_relative)

   END FUNCTION interpolate_native_grid_gradient

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
!> \brief Compact support radius of one hard-minus-soft atom-grid field.
!> \param grid_atom radial and angular source grid
!> \param rho_h hard one-center density values
!> \param rho_s soft one-center density values
!> \param drho_h hard one-center density-gradient values
!> \param drho_s soft one-center density-gradient values
!> \param tau_h hard one-center kinetic-energy-density values
!> \param tau_s soft one-center kinetic-energy-density values
!> \return outermost radius required by any primitive hard-minus-soft field
! **************************************************************************************************
   FUNCTION gapw_atom_grid_support_radius( &
      grid_atom, rho_h, rho_s, drho_h, drho_s, tau_h, tau_s) RESULT(cutoff)
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: rho_h, rho_s
      REAL(dp), DIMENSION(:, :, :, :), INTENT(IN)        :: drho_h, drho_s
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: tau_h, tau_s
      REAL(dp)                                           :: cutoff

      INTEGER                                            :: ir, support_ir
      LOGICAL                                            :: descending
      REAL(dp)                                           :: amplitude

      cutoff = 0.0_dp
      descending = grid_atom%rad(1) > grid_atom%rad(grid_atom%nr)
      IF (descending) THEN
         DO ir = 1, grid_atom%nr
            amplitude = MAX( &
                        MAXVAL(ABS(rho_h(:, ir, :) - rho_s(:, ir, :))), &
                        MAXVAL(ABS(drho_h(1:3, :, ir, :) - drho_s(1:3, :, ir, :))), &
                        MAXVAL(ABS(tau_h(:, ir, :) - tau_s(:, ir, :))))
            IF (amplitude > 1.0E-12_dp) THEN
               support_ir = MAX(1, ir - 2)
               cutoff = grid_atom%rad(support_ir)
               RETURN
            END IF
         END DO
      ELSE
         DO ir = grid_atom%nr, 1, -1
            amplitude = MAX( &
                        MAXVAL(ABS(rho_h(:, ir, :) - rho_s(:, ir, :))), &
                        MAXVAL(ABS(drho_h(1:3, :, ir, :) - drho_s(1:3, :, ir, :))), &
                        MAXVAL(ABS(tau_h(:, ir, :) - tau_s(:, ir, :))))
            IF (amplitude > 1.0E-12_dp) THEN
               support_ir = MIN(grid_atom%nr, ir + 2)
               cutoff = grid_atom%rad(support_ir)
               RETURN
            END IF
         END DO
      END IF

   END FUNCTION gapw_atom_grid_support_radius

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
!> \param radial_derivative_weights radial derivatives of the interpolation weights
!> \param nradial number of active radial nodes
!> \param angular_weights angular interpolation weights
!> \param angular_derivative_weights Cartesian derivatives of the angular weights
!> \param active whether the point lies inside the compact support
! **************************************************************************************************
   SUBROUTINE atom_grid_interpolation_weights( &
      grid_atom, harmonics, displacement, cutoff, radial_indices, radial_weights, &
      radial_derivative_weights, nradial, angular_weights, angular_derivative_weights, active)
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: displacement
      REAL(dp), INTENT(IN)                               :: cutoff
      INTEGER, DIMENSION(4), INTENT(OUT)                 :: radial_indices
      REAL(dp), DIMENSION(4), INTENT(OUT)                :: radial_weights, radial_derivative_weights
      INTEGER, INTENT(OUT)                               :: nradial
      REAL(dp), DIMENSION(:), INTENT(OUT)                :: angular_weights
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: angular_derivative_weights
      LOGICAL, INTENT(OUT)                               :: active

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
      CPASSERT(SIZE(angular_derivative_weights, 1) == 3)
      CPASSERT(SIZE(angular_derivative_weights, 2) == grid_atom%ng_sphere)

      radial_indices = 0
      radial_weights = 0.0_dp
      radial_derivative_weights = 0.0_dp
      angular_weights = 0.0_dp
      angular_derivative_weights = 0.0_dp
      nradial = 0
      active = .FALSE.
      n = grid_atom%nr
      descending = grid_atom%rad(1) > grid_atom%rad(n)
      radius = SQRT(SUM(displacement**2))
      IF (radius > cutoff .OR. &
          radius > grid_atom%rad(MERGE(1, n, descending))) RETURN
      IF (radius <= 1.0E-12_dp) RETURN
      direction = displacement/radius

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
      dh00 = -30.0_dp*t**2 + 60.0_dp*t**3 - 30.0_dp*t**4
      dh10 = 1.0_dp - 18.0_dp*t**2 + 32.0_dp*t**3 - 15.0_dp*t**4
      dh20 = 0.5_dp*(2.0_dp*t - 9.0_dp*t**2 + 12.0_dp*t**3 - 5.0_dp*t**4)
      dh01 = 30.0_dp*t**2 - 60.0_dp*t**3 + 30.0_dp*t**4
      dh11 = -12.0_dp*t**2 + 28.0_dp*t**3 - 15.0_dp*t**4
      dh21 = 0.5_dp*(3.0_dp*t**2 - 8.0_dp*t**3 + 5.0_dp*t**4)

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
      radial_derivative_weights(left_pos) = &
         radial_derivative_weights(left_pos) + dh00/h_interval
      radial_derivative_weights(right_pos) = &
         radial_derivative_weights(right_pos) + dh01/h_interval
      radial_derivative_weights(1:nradial) = radial_derivative_weights(1:nradial) + &
                                             dh10*slope_left(1:nradial) + &
                                             dh11*slope_right(1:nradial) + &
                                             h_interval*(dh20*curvature_left(1:nradial) + &
                                                         dh21*curvature_right(1:nradial))

      angular_values = 0.0_dp
      angular_value_derivatives = 0.0_dp
      DO iso = 1, harmonics%max_s_harm
         l = indso(1, iso)
         CALL y_lm(direction, angular_values(iso), l, indso(2, iso))
         IF (l == 0) CYCLE
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
         DO inode = 1, 3
            angular_derivative_weights(inode, ia) = grid_atom%wa(ia)* &
                                                    DOT_PRODUCT(harmonics%slm(ia, 1:harmonics%max_s_harm), &
                                                                angular_value_derivatives(inode, :))
         END DO
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
! **************************************************************************************************
   SUBROUTINE interpolate_gapw_atom_grid_fields( &
      grid_atom, harmonics, displacement, cutoff, nspins, rho_h, rho_s, drho_h, drho_s, &
      tau_h, tau_s, density, gradient, kin, density_spatial, gradient_spatial, kin_spatial)
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

      INTEGER                                            :: ia, idir, inode, ir, ispin, jdir, nradial
      INTEGER, DIMENSION(4)                              :: radial_indices
      LOGICAL                                            :: active
      REAL(dp)                                           :: derivative_weight, weight
      REAL(dp), DIMENSION(grid_atom%ng_sphere)           :: angular_weights
      REAL(dp), DIMENSION(4)                             :: radial_derivative_weights, radial_weights
      REAL(dp), DIMENSION(3, grid_atom%ng_sphere)        :: angular_derivative_weights

      density = 0.0_dp
      gradient = 0.0_dp
      kin = 0.0_dp
      density_spatial = 0.0_dp
      gradient_spatial = 0.0_dp
      kin_spatial = 0.0_dp
      CALL atom_grid_interpolation_weights( &
         grid_atom, harmonics, displacement, cutoff, radial_indices, radial_weights, &
         radial_derivative_weights, nradial, angular_weights, angular_derivative_weights, active)
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
                     derivative_weight = radial_derivative_weights(inode)* &
                                         displacement(idir)/SQRT(SUM(displacement**2))* &
                                         angular_weights(ia) + radial_weights(inode)* &
                                         angular_derivative_weights(idir, ia)
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
      REAL(dp), DIMENSION(grid_atom%ng_sphere)           :: angular_weights
      REAL(dp), DIMENSION(4)                             :: radial_derivative_weights, radial_weights
      REAL(dp), DIMENSION(3, grid_atom%ng_sphere)        :: angular_derivative_weights

      CALL atom_grid_interpolation_weights( &
         grid_atom, harmonics, displacement, cutoff, radial_indices, radial_weights, &
         radial_derivative_weights, nradial, angular_weights, angular_derivative_weights, active)
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
                                                            nr, nspins, num_pe
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: accint, donlcc, gradient_f, lsd, nlcc, &
                                                            paw_atom, tau_f
      REAL(dp)                                           :: agr, alpha, density_cut, exc_h, exc_s, &
                                                            gradient_cut, tau_cut
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
            IF (accint) THEN
               weight_h => grid_atom%weight
               alpha = dft_control%qs_control%gapw_control%aw(ikind)
               IF (ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                  IF (grid_atom%gapw_weight_alpha /= alpha) DEALLOCATE (grid_atom%gapw_weight_s)
               END IF
               IF (.NOT. ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                  ALLOCATE (grid_atom%gapw_weight_s(na, nr))
                  DO ir = 1, nr
                     agr = 1.0_dp - EXP(-alpha*grid_atom%rad2(ir))
                     grid_atom%gapw_weight_s(:, ir) = grid_atom%weight(:, ir)*agr
                  END DO
                  grid_atom%gapw_weight_alpha = alpha
               END IF
               weight_s => grid_atom%gapw_weight_s
            ELSE
               weight_h => grid_atom%weight
               weight_s => grid_atom%weight
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

! **************************************************************************************************
!> \brief ...
!> \param rho_atom_set ...
!> \param rho1_atom_set ...
!> \param qs_env ...
!> \param xc_section ...
!> \param para_env ...
!> \param do_tddfpt2 New implementation of TDDFT.
!> \param do_triplet ...
!> \param do_sf ...
!> \param kind_set_external ...
! **************************************************************************************************
   SUBROUTINE calculate_xc_2nd_deriv_atom(rho_atom_set, rho1_atom_set, qs_env, xc_section, para_env, &
                                          do_tddfpt2, do_triplet, do_sf, kind_set_external)

      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set, rho1_atom_set
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_tddfpt2, do_triplet, do_sf
      TYPE(qs_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: kind_set_external

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_xc_2nd_deriv_atom'

      INTEGER                                            :: atom, handle, iatom, ikind, ir, na, &
                                                            natom, nr, nspins
      INTEGER, DIMENSION(2)                              :: local_loop_limit
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: accint, donlcc, gradient_functional, &
                                                            lsd, my_do_sf, nlcc, paw_atom, &
                                                            scale_rho, tau_f
      REAL(KIND=dp)                                      :: agr, alpha, density_cut, gradient_cut, &
                                                            rtot, tau_cut
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: vtau_h, vtau_s, vxc_h, vxc_s
      REAL(KIND=dp), DIMENSION(1, 1, 1)                  :: rtau
      REAL(KIND=dp), DIMENSION(1, 1, 1, 1)               :: rrho
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rho_nlcc, weight_h, weight_s
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho1_h, rho1_s, rho_h, rho_s, tau1_h, &
                                                            tau1_s, tau_h, tau_s
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: drho1_h, drho1_s, drho_h, drho_s, vxg_h, &
                                                            vxg_s
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: my_kind_set, qs_kind_set
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: dr1_h, dr1_s, dr_h, dr_s, int_hh, &
                                                            int_ss, r1_h, r1_s, r_h, r_s
      TYPE(rho_atom_coeff), DIMENSION(:, :), POINTER     :: r1_h_d, r1_s_d, r_h_d, r_s_d
      TYPE(rho_atom_type), POINTER                       :: rho1_atom, rho_atom
      TYPE(section_vals_type), POINTER                   :: input, xc_fun_section
      TYPE(tau_basis_cache_type)                         :: tau_basis_cache
      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho1_set_h, rho1_set_s, rho_set_h, &
                                                            rho_set_s

! -------------------------------------------------------------------------

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set)
      NULLIFY (rho_h, rho_s, drho_h, drho_s, weight_h, weight_s)
      NULLIFY (rho1_h, rho1_s, drho1_h, drho1_s)
      NULLIFY (vxc_h, vxc_s, vxg_h, vxg_s)
      NULLIFY (tau_h, tau_s, tau1_h, tau1_s, vtau_h, vtau_s)
      NULLIFY (rho_nlcc)

      CALL get_qs_env(qs_env=qs_env, &
                      input=input, &
                      dft_control=dft_control, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set)

      IF (PRESENT(kind_set_external)) THEN
         my_kind_set => kind_set_external
      ELSE
         my_kind_set => qs_kind_set
      END IF
      nlcc = has_nlcc(my_kind_set)

      accint = dft_control%qs_control%gapw_control%accurate_xcint

      CALL section_vals_val_get(input, "DFT%LSD", l_val=lsd)
      CALL section_vals_val_get(xc_section, "DENSITY_CUTOFF", &
                                r_val=density_cut)
      CALL section_vals_val_get(xc_section, "GRADIENT_CUTOFF", &
                                r_val=gradient_cut)
      CALL section_vals_val_get(xc_section, "TAU_CUTOFF", &
                                r_val=tau_cut)

      my_do_sf = .FALSE.
      IF (PRESENT(do_sf)) my_do_sf = do_sf

      xc_fun_section => section_vals_get_subs_vals(xc_section, &
                                                   "XC_FUNCTIONAL")
      IF (lsd) THEN
         nspins = 2
      ELSE
         nspins = 1
      END IF

      scale_rho = .FALSE.
      IF (PRESENT(do_tddfpt2) .AND. PRESENT(do_triplet)) THEN
         IF (nspins == 1 .AND. do_triplet) THEN
            lsd = .TRUE.
            scale_rho = .TRUE.
         END IF
      ELSE IF (PRESENT(do_triplet)) THEN
         IF (nspins == 1 .AND. do_triplet) lsd = .TRUE.
      END IF

      needs = xc_functionals_get_needs(xc_fun_section, lsd=lsd, &
                                       calc_potential=.TRUE.)
      gradient_functional = needs%drho .OR. needs%drho_spin
      tau_f = (needs%tau .OR. needs%tau_spin)
      IF (.NOT. tau_f) rtau = 0.0_dp

      !  Here starts the loop over all the atoms
      DO ikind = 1, SIZE(atomic_kind_set)

         NULLIFY (atom_list, harmonics, grid_atom)
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
         CALL get_qs_kind(my_kind_set(ikind), paw_atom=paw_atom, &
                          harmonics=harmonics, grid_atom=grid_atom)
         CALL get_qs_kind(my_kind_set(ikind), basis_set=basis_1c, basis_type="GAPW_1C")
         IF (.NOT. paw_atom) CYCLE

         nr = grid_atom%nr
         na = grid_atom%ng_sphere

         ! set integration weights
         IF (accint) THEN
            weight_h => grid_atom%weight
            alpha = dft_control%qs_control%gapw_control%aw(ikind)
            IF (ASSOCIATED(grid_atom%gapw_weight_s)) THEN
               IF (grid_atom%gapw_weight_alpha /= alpha) DEALLOCATE (grid_atom%gapw_weight_s)
            END IF
            IF (.NOT. ASSOCIATED(grid_atom%gapw_weight_s)) THEN
               ALLOCATE (grid_atom%gapw_weight_s(na, nr))
               DO ir = 1, nr
                  agr = 1.0_dp - EXP(-alpha*grid_atom%rad2(ir))
                  grid_atom%gapw_weight_s(:, ir) = grid_atom%weight(:, ir)*agr
               END DO
               grid_atom%gapw_weight_alpha = alpha
            END IF
            weight_s => grid_atom%gapw_weight_s
         ELSE
            weight_h => grid_atom%weight
            weight_s => grid_atom%weight
         END IF

         ! Array dimension: here anly one dimensional arrays are used,
         ! i.e. only the first column of deriv_data is read.
         ! The other to dimensions  are set to size equal 1.
         bounds(1:2, 1:3) = 1
         bounds(2, 1) = na
         bounds(2, 2) = nr

         CALL xc_dset_create(deriv_set, local_bounds=bounds)
         CALL xc_rho_set_create(rho_set_h, bounds, rho_cutoff=density_cut, &
                                drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
         CALL xc_rho_set_create(rho_set_s, bounds, rho_cutoff=density_cut, &
                                drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
         CALL xc_rho_set_create(rho1_set_h, bounds, rho_cutoff=density_cut, &
                                drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
         CALL xc_rho_set_create(rho1_set_s, bounds, rho_cutoff=density_cut, &
                                drho_cutoff=gradient_cut, tau_cutoff=tau_cut)

         ! allocate the required 3d arrays where to store rho and drho
         IF (nspins == 1 .AND. .NOT. lsd) THEN
            CALL xc_rho_set_atom_update(rho_set_h, needs, 1, bounds)
            CALL xc_rho_set_atom_update(rho1_set_h, needs, 1, bounds)
            CALL xc_rho_set_atom_update(rho_set_s, needs, 1, bounds)
            CALL xc_rho_set_atom_update(rho1_set_s, needs, 1, bounds)
         ELSE
            CALL xc_rho_set_atom_update(rho_set_h, needs, 2, bounds)
            CALL xc_rho_set_atom_update(rho1_set_h, needs, 2, bounds)
            CALL xc_rho_set_atom_update(rho_set_s, needs, 2, bounds)
            CALL xc_rho_set_atom_update(rho1_set_s, needs, 2, bounds)
         END IF

         ALLOCATE (rho_h(1:na, 1:nr, 1:nspins), rho1_h(1:na, 1:nr, 1:nspins), &
                   rho_s(1:na, 1:nr, 1:nspins), rho1_s(1:na, 1:nr, 1:nspins))

         ALLOCATE (vxc_h(1:na, 1:nr, 1:nspins), vxc_s(1:na, 1:nr, 1:nspins))
         vxc_h = 0.0_dp
         vxc_s = 0.0_dp

         IF (tau_f) THEN
            CALL create_tau_basis_cache(tau_basis_cache, grid_atom, basis_1c, harmonics)
            ALLOCATE (tau_h(1:na, 1:nr, 1:nspins), tau1_h(1:na, 1:nr, 1:nspins), &
                      tau_s(1:na, 1:nr, 1:nspins), tau1_s(1:na, 1:nr, 1:nspins))
            ALLOCATE (vtau_h(1:na, 1:nr, 1:nspins), vtau_s(1:na, 1:nr, 1:nspins))
         END IF

         IF (gradient_functional) THEN
            ALLOCATE (drho_h(1:4, 1:na, 1:nr, 1:nspins), drho1_h(1:4, 1:na, 1:nr, 1:nspins), &
                      drho_s(1:4, 1:na, 1:nr, 1:nspins), drho1_s(1:4, 1:na, 1:nr, 1:nspins))
            ALLOCATE (vxg_h(1:3, 1:na, 1:nr, 1:nspins), vxg_s(1:3, 1:na, 1:nr, 1:nspins))
         ELSE
            ALLOCATE (drho_h(1, 1, 1, 1), drho1_h(1, 1, 1, 1), &
                      drho_s(1, 1, 1, 1), drho1_s(1, 1, 1, 1))
            ALLOCATE (vxg_h(1, 1, 1, 1), vxg_s(1, 1, 1, 1))
            rrho = 0.0_dp
         END IF
         vxg_h = 0.0_dp
         vxg_s = 0.0_dp

         ! NLCC: prepare rho and drho of the core charge for this KIND
         donlcc = .FALSE.
         IF (nlcc) THEN
            NULLIFY (rho_nlcc)
            rho_nlcc => my_kind_set(ikind)%nlcc_pot
            IF (ASSOCIATED(rho_nlcc)) donlcc = .TRUE.
         END IF

         ! parallelization
         local_loop_limit = get_limit(natom, para_env%num_pe, para_env%mepos)

         DO iatom = local_loop_limit(1), local_loop_limit(2) !1,natom
            atom = atom_list(iatom)

            rho_atom_set(atom)%exc_h = 0.0_dp
            rho_atom_set(atom)%exc_s = 0.0_dp
            rho1_atom_set(atom)%exc_h = 0.0_dp
            rho1_atom_set(atom)%exc_s = 0.0_dp

            rho_atom => rho_atom_set(atom)
            rho1_atom => rho1_atom_set(atom)
            NULLIFY (r_h, r_s, dr_h, dr_s, r_h_d, r_s_d)
            NULLIFY (r1_h, r1_s, dr1_h, dr1_s, r1_h_d, r1_s_d)
            rho_h = 0.0_dp
            rho_s = 0.0_dp
            rho1_h = 0.0_dp
            rho1_s = 0.0_dp
            IF (gradient_functional) THEN
               CALL get_rho_atom(rho_atom=rho_atom, &
                                 rho_rad_h=r_h, rho_rad_s=r_s, &
                                 drho_rad_h=dr_h, drho_rad_s=dr_s, &
                                 rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d)
               CALL get_rho_atom(rho_atom=rho1_atom, &
                                 rho_rad_h=r1_h, rho_rad_s=r1_s, &
                                 drho_rad_h=dr1_h, drho_rad_s=dr1_s, &
                                 rho_rad_h_d=r1_h_d, rho_rad_s_d=r1_s_d)
               drho_h = 0.0_dp; drho_s = 0.0_dp
               drho1_h = 0.0_dp; drho1_s = 0.0_dp
            ELSE
               CALL get_rho_atom(rho_atom=rho_atom, &
                                 rho_rad_h=r_h, rho_rad_s=r_s)
               CALL get_rho_atom(rho_atom=rho1_atom, &
                                 rho_rad_h=r1_h, rho_rad_s=r1_s)
            END IF

            rtot = 0.0_dp

            DO ir = 1, nr
               CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_functional, &
                                     ir, r_h, r_s, rho_h, rho_s, dr_h, dr_s, r_h_d, r_s_d, &
                                     drho_h, drho_s)
               IF (donlcc) THEN
                  CALL calc_rho_nlcc(grid_atom, nspins, gradient_functional, &
                                     ir, rho_nlcc(:, 1), rho_h, rho_s, rho_nlcc(:, 2), drho_h, drho_s)
               END IF
               CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_functional, &
                                     ir, r1_h, r1_s, rho1_h, rho1_s, dr1_h, dr1_s, r1_h_d, r1_s_d, &
                                     drho1_h, drho1_s)
            END DO
            IF (tau_f) THEN
               CALL calc_tau_atom(tau_h, tau_s, rho_atom, tau_basis_cache, nspins)
               CALL calc_tau_atom(tau1_h, tau1_s, rho1_atom, tau_basis_cache, nspins)
            END IF
            IF (scale_rho) THEN
               rho_h = 2.0_dp*rho_h
               rho_s = 2.0_dp*rho_s
               IF (gradient_functional) THEN
                  drho_h = 2.0_dp*drho_h
                  drho_s = 2.0_dp*drho_s
               END IF
               IF (tau_f) THEN
                  tau_h = 2.0_dp*tau_h
                  tau_s = 2.0_dp*tau_s
               END IF
            END IF

            DO ir = 1, nr
               IF (tau_f) THEN
                  CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, drho_h, tau_h, na, ir)
                  CALL fill_rho_set(rho1_set_h, lsd, nspins, needs, rho1_h, drho1_h, tau1_h, na, ir)
                  CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, drho_s, tau_s, na, ir)
                  CALL fill_rho_set(rho1_set_s, lsd, nspins, needs, rho1_s, drho1_s, tau1_s, na, ir)
               ELSE IF (gradient_functional) THEN
                  CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, drho_h, rtau, na, ir)
                  CALL fill_rho_set(rho1_set_h, lsd, nspins, needs, rho1_h, drho1_h, rtau, na, ir)
                  CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, drho_s, rtau, na, ir)
                  CALL fill_rho_set(rho1_set_s, lsd, nspins, needs, rho1_s, drho1_s, rtau, na, ir)
               ELSE
                  CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, rrho, rtau, na, ir)
                  CALL fill_rho_set(rho1_set_h, lsd, nspins, needs, rho1_h, rrho, rtau, na, ir)
                  CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, rrho, rtau, na, ir)
                  CALL fill_rho_set(rho1_set_s, lsd, nspins, needs, rho1_s, rrho, rtau, na, ir)
               END IF
            END DO

            CALL xc_2nd_deriv_of_r(xc_section=xc_section, &
                                   rho_set=rho_set_h, rho1_set=rho1_set_h, &
                                   deriv_set=deriv_set, &
                                   w=weight_h, vxc=vxc_h, vxg=vxg_h, vtau=vtau_h, do_triplet=do_triplet, &
                                   do_sf=my_do_sf)
            CALL xc_2nd_deriv_of_r(xc_section=xc_section, &
                                   rho_set=rho_set_s, rho1_set=rho1_set_s, &
                                   deriv_set=deriv_set, &
                                   w=weight_s, vxc=vxc_s, vxg=vxg_s, vtau=vtau_s, do_triplet=do_triplet, &
                                   do_sf=my_do_sf)

            CALL get_rho_atom(rho_atom=rho1_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
            IF (gradient_functional) THEN
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

         END DO

         ! some cleanup
         DEALLOCATE (rho_h, rho_s, rho1_h, rho1_s, vxc_h, vxc_s)
         DEALLOCATE (drho_h, drho_s, vxg_h, vxg_s)
         DEALLOCATE (drho1_h, drho1_s)
         IF (tau_f) THEN
            DEALLOCATE (tau_h, tau_s, tau1_h, tau1_s)
            DEALLOCATE (vtau_h, vtau_s)
            CALL release_tau_basis_cache(tau_basis_cache)
         END IF

         CALL xc_dset_release(deriv_set)
         CALL xc_rho_set_release(rho_set_h)
         CALL xc_rho_set_release(rho1_set_h)
         CALL xc_rho_set_release(rho_set_s)
         CALL xc_rho_set_release(rho1_set_s)
      END DO

      CALL timestop(handle)

   END SUBROUTINE calculate_xc_2nd_deriv_atom

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho0_atom_set ...
!> \param rho1_atom_set ...
!> \param rho2_atom_set ...
!> \param kind_set ...
!> \param xc_section ...
!> \param is_triplet ...
!> \param accuracy ...
! **************************************************************************************************
   SUBROUTINE calculate_gfxc_atom(qs_env, rho0_atom_set, rho1_atom_set, rho2_atom_set, &
                                  kind_set, xc_section, is_triplet, accuracy)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho0_atom_set, rho1_atom_set, &
                                                            rho2_atom_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: kind_set
      TYPE(section_vals_type), OPTIONAL, POINTER         :: xc_section
      LOGICAL, INTENT(IN)                                :: is_triplet
      INTEGER, INTENT(IN)                                :: accuracy

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_gfxc_atom'
      REAL(KIND=dp), PARAMETER                           :: epsrho = 5.e-4_dp

      INTEGER                                            :: bo(2), handle, iat, iatom, ikind, ir, &
                                                            istep, mspins, myfun, na, natom, nf, &
                                                            nr, ns, nspins, nstep, num_pe
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: accint, donlcc, gradient_f, lsd, nlcc, &
                                                            paw_atom, tau_f
      REAL(dp)                                           :: agr, alpha, beta, density_cut, exc_h, &
                                                            exc_s, gradient_cut, oeps1, oeps2, &
                                                            tau_cut
      REAL(dp), DIMENSION(1, 1, 1)                       :: tau_d
      REAL(dp), DIMENSION(1, 1, 1, 1)                    :: rho_d
      REAL(dp), DIMENSION(:, :), POINTER                 :: rho_nlcc, weight_h, weight_s
      REAL(dp), DIMENSION(:, :, :), POINTER              :: rho0_h, rho0_s, rho1_h, rho1_s, rho_h, &
                                                            rho_s, tau0_h, tau0_s, tau1_h, tau1_s, &
                                                            tau_h, tau_s, vtau_h, vtau_s, vxc_h, &
                                                            vxc_s
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: drho0_h, drho0_s, drho1_h, drho1_s, &
                                                            drho_h, drho_s, vxg_h, vxg_s
      REAL(KIND=dp), DIMENSION(-4:4)                     :: ak, bl
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: dr_h, dr_s, fint_hh, fint_ss, int_hh, &
                                                            int_ss, r_h, r_s
      TYPE(rho_atom_coeff), DIMENSION(:, :), POINTER     :: r_h_d, r_s_d
      TYPE(rho_atom_type), POINTER                       :: rho0_atom, rho1_atom, rho2_atom
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(tau_basis_cache_type)                         :: tau_basis_cache
      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set_h, rho_set_s

      CALL timeset(routineN, handle)

      NULLIFY (vtau_h, vtau_s)

      ak = 0.0_dp
      bl = 0.0_dp
      SELECT CASE (accuracy)
      CASE (:4)
         nstep = 2
         ak(-2:2) = [1.0_dp, -8.0_dp, 0.0_dp, 8.0_dp, -1.0_dp]/12.0_dp
         bl(-2:2) = [-1.0_dp, 16.0_dp, -30.0_dp, 16.0_dp, -1.0_dp]/12.0_dp
      CASE (5:7)
         nstep = 3
         ak(-3:3) = [-1.0_dp, 9.0_dp, -45.0_dp, 0.0_dp, 45.0_dp, -9.0_dp, 1.0_dp]/60.0_dp
         bl(-3:3) = [2.0_dp, -27.0_dp, 270.0_dp, -490.0_dp, 270.0_dp, -27.0_dp, 2.0_dp]/180.0_dp
      CASE (8:)
         nstep = 4
         ak(-4:4) = [1.0_dp, -32.0_dp/3.0_dp, 56.0_dp, -224.0_dp, 0.0_dp, &
                     224.0_dp, -56.0_dp, 32.0_dp/3.0_dp, -1.0_dp]/280.0_dp
         bl(-4:4) = [-1.0_dp, 128.0_dp/9.0_dp, -112.0_dp, 896.0_dp, -14350.0_dp/9.0_dp, &
                     896.0_dp, -112.0_dp, 128.0_dp/9.0_dp, -1.0_dp]/560.0_dp
      END SELECT
      oeps1 = 1.0_dp/epsrho
      oeps2 = 1.0_dp/(epsrho**2)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      atomic_kind_set=atomic_kind_set)

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      CALL section_vals_val_get(xc_fun_section, "_SECTION_PARAMETERS_", i_val=myfun)

      accint = dft_control%qs_control%gapw_control%accurate_xcint

      IF (myfun == xc_none) THEN
         ! no action needed?
      ELSE
         CALL section_vals_val_get(xc_section, "DENSITY_CUTOFF", r_val=density_cut)
         CALL section_vals_val_get(xc_section, "GRADIENT_CUTOFF", r_val=gradient_cut)
         CALL section_vals_val_get(xc_section, "TAU_CUTOFF", r_val=tau_cut)

         nlcc = has_nlcc(kind_set)
         lsd = dft_control%lsd
         nspins = dft_control%nspins
         mspins = nspins
         IF (is_triplet) THEN
            CPASSERT(nspins == 1)
            lsd = .TRUE.
            mspins = 2
         END IF
         needs = xc_functionals_get_needs(xc_fun_section, lsd=lsd, calc_potential=.TRUE.)
         gradient_f = (needs%drho .OR. needs%drho_spin)
         tau_f = (needs%tau .OR. needs%tau_spin)

         ! Here starts the loop over all the atoms
         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
            CALL get_qs_kind(kind_set(ikind), paw_atom=paw_atom, &
                             harmonics=harmonics, grid_atom=grid_atom)
            CALL get_qs_kind(kind_set(ikind), basis_set=basis_1c, basis_type="GAPW_1C")

            IF (.NOT. paw_atom) CYCLE

            nr = grid_atom%nr
            na = grid_atom%ng_sphere

            ! set integration weights
            IF (accint) THEN
               weight_h => grid_atom%weight
               alpha = dft_control%qs_control%gapw_control%aw(ikind)
               IF (ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                  IF (grid_atom%gapw_weight_alpha /= alpha) DEALLOCATE (grid_atom%gapw_weight_s)
               END IF
               IF (.NOT. ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                  ALLOCATE (grid_atom%gapw_weight_s(na, nr))
                  DO ir = 1, nr
                     agr = 1.0_dp - EXP(-alpha*grid_atom%rad2(ir))
                     grid_atom%gapw_weight_s(:, ir) = grid_atom%weight(:, ir)*agr
                  END DO
                  grid_atom%gapw_weight_alpha = alpha
               END IF
               weight_s => grid_atom%gapw_weight_s
            ELSE
               weight_h => grid_atom%weight
               weight_s => grid_atom%weight
            END IF

            ! Prepare the structures needed to calculate and store the xc derivatives

            ! Array dimension: here anly one dimensional arrays are used,
            ! i.e. only the first column of deriv_data is read.
            ! The other to dimensions  are set to size equal 1
            bounds(1:2, 1:3) = 1
            bounds(2, 1) = na
            bounds(2, 2) = nr

            ! create a place where to put the derivatives
            CALL xc_dset_create(deriv_set, local_bounds=bounds)
            ! create the place where to store the argument for the functionals
            CALL xc_rho_set_create(rho_set_h, bounds, rho_cutoff=density_cut, &
                                   drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
            CALL xc_rho_set_create(rho_set_s, bounds, rho_cutoff=density_cut, &
                                   drho_cutoff=gradient_cut, tau_cutoff=tau_cut)

            ! allocate the required 3d arrays where to store rho and drho
            CALL xc_rho_set_atom_update(rho_set_h, needs, mspins, bounds)
            CALL xc_rho_set_atom_update(rho_set_s, needs, mspins, bounds)

            ALLOCATE (rho_h(na, nr, mspins), rho_s(na, nr, mspins), &
                      rho0_h(na, nr, nspins), rho0_s(na, nr, nspins), &
                      rho1_h(na, nr, nspins), rho1_s(na, nr, nspins))
            ALLOCATE (vxc_h(na, nr, mspins), vxc_s(na, nr, mspins))
            IF (gradient_f) THEN
               ALLOCATE (drho_h(4, na, nr, mspins), drho_s(4, na, nr, mspins), &
                         drho0_h(4, na, nr, nspins), drho0_s(4, na, nr, nspins), &
                         drho1_h(4, na, nr, nspins), drho1_s(4, na, nr, nspins))
               ALLOCATE (vxg_h(3, na, nr, mspins), vxg_s(3, na, nr, mspins))
            END IF
            IF (tau_f) THEN
               CALL create_tau_basis_cache(tau_basis_cache, grid_atom, basis_1c, harmonics)
               ALLOCATE (tau_h(na, nr, mspins), tau_s(na, nr, mspins), &
                         tau0_h(na, nr, nspins), tau0_s(na, nr, nspins), &
                         tau1_h(na, nr, nspins), tau1_s(na, nr, nspins))
               ALLOCATE (vtau_h(na, nr, mspins), vtau_s(na, nr, mspins))
            END IF
            !
            ! NLCC: prepare rho and drho of the core charge for this KIND
            donlcc = .FALSE.
            IF (nlcc) THEN
               NULLIFY (rho_nlcc)
               rho_nlcc => kind_set(ikind)%nlcc_pot
               IF (ASSOCIATED(rho_nlcc)) donlcc = .TRUE.
            END IF

            ! Distribute the atoms of this kind
            num_pe = para_env%num_pe
            bo = get_limit(natom, num_pe, para_env%mepos)

            DO iat = bo(1), bo(2)
               iatom = atom_list(iat)
               !
               NULLIFY (int_hh, int_ss)
               rho0_atom => rho0_atom_set(iatom)
               CALL get_rho_atom(rho_atom=rho0_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
               ALLOCATE (fint_ss(nspins), fint_hh(nspins))
               DO ns = 1, nspins
                  nf = SIZE(int_ss(ns)%r_coef, 1)
                  ALLOCATE (fint_ss(ns)%r_coef(nf, nf))
                  nf = SIZE(int_hh(ns)%r_coef, 1)
                  ALLOCATE (fint_hh(ns)%r_coef(nf, nf))
               END DO

               ! RHO0
               rho0_h = 0.0_dp
               rho0_s = 0.0_dp
               rho0_atom => rho0_atom_set(iatom)
               IF (gradient_f) THEN
                  NULLIFY (r_h, r_s, dr_h, dr_s, r_h_d, r_s_d)
                  CALL get_rho_atom(rho_atom=rho0_atom, rho_rad_h=r_h, rho_rad_s=r_s, drho_rad_h=dr_h, &
                                    drho_rad_s=dr_s, rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d)
                  drho0_h = 0.0_dp
                  drho0_s = 0.0_dp
               ELSE
                  NULLIFY (r_h, r_s)
                  CALL get_rho_atom(rho_atom=rho0_atom, rho_rad_h=r_h, rho_rad_s=r_s)
                  rho_d = 0.0_dp
               END IF
               DO ir = 1, nr
                  CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_f, &
                                        ir, r_h, r_s, rho0_h, rho0_s, dr_h, dr_s, &
                                        r_h_d, r_s_d, drho0_h, drho0_s)
                  IF (donlcc) THEN
                     CALL calc_rho_nlcc(grid_atom, nspins, gradient_f, &
                                        ir, rho_nlcc(:, 1), rho0_h, rho0_s, rho_nlcc(:, 2), drho0_h, drho0_s)
                  END IF
               END DO
               IF (tau_f) THEN
                  !compute tau on the grid all at once
                  CALL calc_tau_atom(tau0_h, tau0_s, rho0_atom, tau_basis_cache, nspins)
               ELSE
                  tau_d = 0.0_dp
               END IF
               ! RHO1
               rho1_h = 0.0_dp
               rho1_s = 0.0_dp
               rho1_atom => rho1_atom_set(iatom)
               IF (gradient_f) THEN
                  NULLIFY (r_h, r_s, dr_h, dr_s, r_h_d, r_s_d)
                  CALL get_rho_atom(rho_atom=rho1_atom, rho_rad_h=r_h, rho_rad_s=r_s, drho_rad_h=dr_h, &
                                    drho_rad_s=dr_s, rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d)
                  drho1_h = 0.0_dp
                  drho1_s = 0.0_dp
               ELSE
                  NULLIFY (r_h, r_s)
                  CALL get_rho_atom(rho_atom=rho1_atom, rho_rad_h=r_h, rho_rad_s=r_s)
               END IF
               DO ir = 1, nr
                  CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_f, &
                                        ir, r_h, r_s, rho1_h, rho1_s, dr_h, dr_s, &
                                        r_h_d, r_s_d, drho1_h, drho1_s)
               END DO
               IF (tau_f) THEN
                  !compute tau on the grid all at once
                  CALL calc_tau_atom(tau1_h, tau1_s, rho1_atom, tau_basis_cache, nspins)
               END IF
               ! RHO2
               rho2_atom => rho2_atom_set(iatom)

               DO istep = -nstep, nstep

                  beta = REAL(istep, KIND=dp)*epsrho

                  IF (is_triplet) THEN
                     rho_h(:, :, 1) = rho0_h(:, :, 1) + beta*rho1_h(:, :, 1)
                     rho_h(:, :, 2) = rho0_h(:, :, 1)
                     rho_h = 0.5_dp*rho_h
                     rho_s(:, :, 1) = rho0_s(:, :, 1) + beta*rho1_s(:, :, 1)
                     rho_s(:, :, 2) = rho0_s(:, :, 1)
                     rho_s = 0.5_dp*rho_s
                     IF (gradient_f) THEN
                        drho_h(:, :, :, 1) = drho0_h(:, :, :, 1) + beta*drho1_h(:, :, :, 1)
                        drho_h(:, :, :, 2) = drho0_h(:, :, :, 1)
                        drho_h = 0.5_dp*drho_h
                        drho_s(:, :, :, 1) = drho0_s(:, :, :, 1) + beta*drho1_s(:, :, :, 1)
                        drho_s(:, :, :, 2) = drho0_s(:, :, :, 1)
                        drho_s = 0.5_dp*drho_s
                     END IF
                     IF (tau_f) THEN
                        tau_h(:, :, 1) = tau0_h(:, :, 1) + beta*tau1_h(:, :, 1)
                        tau_h(:, :, 2) = tau0_h(:, :, 1)
                        tau_h = 0.5_dp*tau0_h
                        tau_s(:, :, 1) = tau0_s(:, :, 1) + beta*tau1_s(:, :, 1)
                        tau_s(:, :, 2) = tau0_s(:, :, 1)
                        tau_s = 0.5_dp*tau0_s
                     END IF
                  ELSE
                     rho_h = rho0_h + beta*rho1_h
                     rho_s = rho0_s + beta*rho1_s
                     IF (gradient_f) THEN
                        drho_h = drho0_h + beta*drho1_h
                        drho_s = drho0_s + beta*drho1_s
                     END IF
                     IF (tau_f) THEN
                        tau_h = tau0_h + beta*tau1_h
                        tau_s = tau0_s + beta*tau1_s
                     END IF
                  END IF
                  !
                  IF (gradient_f) THEN
                     drho_h(4, :, :, :) = SQRT( &
                                          drho_h(1, :, :, :)*drho_h(1, :, :, :) + &
                                          drho_h(2, :, :, :)*drho_h(2, :, :, :) + &
                                          drho_h(3, :, :, :)*drho_h(3, :, :, :))

                     drho_s(4, :, :, :) = SQRT( &
                                          drho_s(1, :, :, :)*drho_s(1, :, :, :) + &
                                          drho_s(2, :, :, :)*drho_s(2, :, :, :) + &
                                          drho_s(3, :, :, :)*drho_s(3, :, :, :))
                  END IF

                  DO ir = 1, nr
                     IF (tau_f) THEN
                        CALL fill_rho_set(rho_set_h, lsd, mspins, needs, rho_h, drho_h, tau_h, na, ir)
                        CALL fill_rho_set(rho_set_s, lsd, mspins, needs, rho_s, drho_s, tau_s, na, ir)
                     ELSE IF (gradient_f) THEN
                        CALL fill_rho_set(rho_set_h, lsd, mspins, needs, rho_h, drho_h, tau_d, na, ir)
                        CALL fill_rho_set(rho_set_s, lsd, mspins, needs, rho_s, drho_s, tau_d, na, ir)
                     ELSE
                        CALL fill_rho_set(rho_set_h, lsd, mspins, needs, rho_h, rho_d, tau_d, na, ir)
                        CALL fill_rho_set(rho_set_s, lsd, mspins, needs, rho_s, rho_d, tau_d, na, ir)
                     END IF
                  END DO

                  ! hard atom density !
                  CALL xc_dset_zero_all(deriv_set)
                  CALL vxc_of_r_new(xc_fun_section, rho_set_h, deriv_set, 1, needs, weight_h, &
                                    lsd, na, nr, exc_h, vxc_h, vxg_h, vtau_h)
                  IF (is_triplet) THEN
                     vxc_h(:, :, 1) = vxc_h(:, :, 1) - vxc_h(:, :, 2)
                     IF (gradient_f) THEN
                        vxg_h(:, :, :, 1) = vxg_h(:, :, :, 1) - vxg_h(:, :, :, 2)
                     END IF
                     IF (tau_f) THEN
                        vtau_h(:, :, 1) = vtau_h(:, :, 1) - vtau_h(:, :, 2)
                     END IF
                  END IF
                  ! soft atom density !
                  CALL xc_dset_zero_all(deriv_set)
                  CALL vxc_of_r_new(xc_fun_section, rho_set_s, deriv_set, 1, needs, weight_s, &
                                    lsd, na, nr, exc_s, vxc_s, vxg_s, vtau_s)
                  IF (is_triplet) THEN
                     vxc_s(:, :, 1) = vxc_s(:, :, 1) - vxc_s(:, :, 2)
                     IF (gradient_f) THEN
                        vxg_s(:, :, :, 1) = vxg_s(:, :, :, 1) - vxg_s(:, :, :, 2)
                     END IF
                     IF (tau_f) THEN
                        vtau_s(:, :, 1) = vtau_s(:, :, 1) - vtau_s(:, :, 2)
                     END IF
                  END IF
                  ! potentials
                  DO ns = 1, nspins
                     fint_hh(ns)%r_coef(:, :) = 0.0_dp
                     fint_ss(ns)%r_coef(:, :) = 0.0_dp
                  END DO
                  IF (gradient_f) THEN
                     CALL gaVxcgb_GC(vxc_h, vxc_s, vxg_h, vxg_s, fint_hh, fint_ss, &
                                     grid_atom, basis_1c, harmonics, nspins)
                  ELSE
                     CALL gaVxcgb_noGC(vxc_h, vxc_s, fint_hh, fint_ss, &
                                       grid_atom, basis_1c, harmonics, nspins)
                  END IF
                  IF (tau_f) THEN
                     CALL dgaVtaudgb(vtau_h, vtau_s, fint_hh, fint_ss, &
                                     tau_basis_cache, nspins)
                  END IF
                  ! first derivative fxc
                  NULLIFY (int_hh, int_ss)
                  CALL get_rho_atom(rho_atom=rho1_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
                  DO ns = 1, nspins
                     int_ss(ns)%r_coef(:, :) = int_ss(ns)%r_coef(:, :) + oeps1*ak(istep)*fint_ss(ns)%r_coef(:, :)
                     int_hh(ns)%r_coef(:, :) = int_hh(ns)%r_coef(:, :) + oeps1*ak(istep)*fint_hh(ns)%r_coef(:, :)
                  END DO
                  ! second derivative gxc
                  NULLIFY (int_hh, int_ss)
                  CALL get_rho_atom(rho_atom=rho2_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
                  DO ns = 1, nspins
                     int_ss(ns)%r_coef(:, :) = int_ss(ns)%r_coef(:, :) + oeps2*bl(istep)*fint_ss(ns)%r_coef(:, :)
                     int_hh(ns)%r_coef(:, :) = int_hh(ns)%r_coef(:, :) + oeps2*bl(istep)*fint_hh(ns)%r_coef(:, :)
                  END DO
               END DO
               !
               DO ns = 1, nspins
                  DEALLOCATE (fint_ss(ns)%r_coef)
                  DEALLOCATE (fint_hh(ns)%r_coef)
               END DO
               DEALLOCATE (fint_ss, fint_hh)

            END DO ! iat

            ! Release the xc structure used to store the xc derivatives
            CALL xc_dset_release(deriv_set)
            CALL xc_rho_set_release(rho_set_h)
            CALL xc_rho_set_release(rho_set_s)

            DEALLOCATE (rho_h, rho_s, rho0_h, rho0_s, rho1_h, rho1_s)
            DEALLOCATE (vxc_h, vxc_s)
            IF (gradient_f) THEN
               DEALLOCATE (drho_h, drho_s, drho0_h, drho0_s, drho1_h, drho1_s)
               DEALLOCATE (vxg_h, vxg_s)
            END IF
            IF (tau_f) THEN
               DEALLOCATE (tau_h, tau_s, tau0_h, tau0_s, tau1_h, tau1_s)
               DEALLOCATE (vtau_h, vtau_s)
               CALL release_tau_basis_cache(tau_basis_cache)
            END IF
         END DO ! ikind

      END IF !xc_none

      CALL timestop(handle)

   END SUBROUTINE calculate_gfxc_atom

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho0_atom_set ...
!> \param rho1_atom_set ...
!> \param rho2_atom_set ...
!> \param kind_set ...
!> \param xc_section ...
!> \param is_triplet ...
!> \param accuracy ...
!> \param epsrho ...
! **************************************************************************************************
   SUBROUTINE gfxc_atom_diff(qs_env, rho0_atom_set, rho1_atom_set, rho2_atom_set, &
                             kind_set, xc_section, is_triplet, accuracy, epsrho)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho0_atom_set, rho1_atom_set, &
                                                            rho2_atom_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: kind_set
      TYPE(section_vals_type), OPTIONAL, POINTER         :: xc_section
      LOGICAL, INTENT(IN)                                :: is_triplet
      INTEGER, INTENT(IN)                                :: accuracy
      REAL(KIND=dp), INTENT(IN)                          :: epsrho

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'gfxc_atom_diff'

      INTEGER                                            :: bo(2), handle, iat, iatom, ikind, ir, &
                                                            istep, mspins, myfun, na, natom, nf, &
                                                            nr, ns, nspins, nstep, num_pe
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: accint, donlcc, gradient_f, lsd, nlcc, &
                                                            paw_atom, tau_f
      REAL(dp)                                           :: agr, alpha, beta, density_cut, &
                                                            gradient_cut, oeps1, tau_cut
      REAL(dp), CONTIGUOUS, DIMENSION(:, :, :), POINTER  :: vtau_h, vtau_s, vxc_h, vxc_s
      REAL(dp), DIMENSION(1, 1, 1)                       :: tau_d
      REAL(dp), DIMENSION(1, 1, 1, 1)                    :: rho_d
      REAL(dp), DIMENSION(:, :), POINTER                 :: rho_nlcc, weight_h, weight_s
      REAL(dp), DIMENSION(:, :, :), POINTER              :: rho0_h, rho0_s, rho1_h, rho1_s, rho_h, &
                                                            rho_s, tau0_h, tau0_s, tau1_h, tau1_s, &
                                                            tau_h, tau_s
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: drho0_h, drho0_s, drho1_h, drho1_s, &
                                                            drho_h, drho_s, vxg_h, vxg_s
      REAL(KIND=dp), DIMENSION(-4:4)                     :: ak
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: dr_h, dr_s, fint_hh, fint_ss, int_hh, &
                                                            int_ss, r_h, r_s
      TYPE(rho_atom_coeff), DIMENSION(:, :), POINTER     :: r_h_d, r_s_d
      TYPE(rho_atom_type), POINTER                       :: rho0_atom, rho1_atom, rho2_atom
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(tau_basis_cache_type)                         :: tau_basis_cache
      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho1_set_h, rho1_set_s, rho_set_h, &
                                                            rho_set_s

      CALL timeset(routineN, handle)

      NULLIFY (vtau_h, vtau_s)

      ak = 0.0_dp
      SELECT CASE (accuracy)
      CASE (:4)
         nstep = 2
         ak(-2:2) = [1.0_dp, -8.0_dp, 0.0_dp, 8.0_dp, -1.0_dp]/12.0_dp
      CASE (5:7)
         nstep = 3
         ak(-3:3) = [-1.0_dp, 9.0_dp, -45.0_dp, 0.0_dp, 45.0_dp, -9.0_dp, 1.0_dp]/60.0_dp
      CASE (8:)
         nstep = 4
         ak(-4:4) = [1.0_dp, -32.0_dp/3.0_dp, 56.0_dp, -224.0_dp, 0.0_dp, &
                     224.0_dp, -56.0_dp, 32.0_dp/3.0_dp, -1.0_dp]/280.0_dp
      END SELECT
      oeps1 = 1.0_dp/epsrho

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      atomic_kind_set=atomic_kind_set)

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      CALL section_vals_val_get(xc_fun_section, "_SECTION_PARAMETERS_", i_val=myfun)

      accint = dft_control%qs_control%gapw_control%accurate_xcint

      IF (myfun == xc_none) THEN
         ! no action needed?
      ELSE
         ! calculate fxc
         CALL calculate_xc_2nd_deriv_atom(rho0_atom_set, rho1_atom_set, qs_env, xc_section, para_env, &
                                          do_triplet=is_triplet, kind_set_external=kind_set)

         CALL section_vals_val_get(xc_section, "DENSITY_CUTOFF", r_val=density_cut)
         CALL section_vals_val_get(xc_section, "GRADIENT_CUTOFF", r_val=gradient_cut)
         CALL section_vals_val_get(xc_section, "TAU_CUTOFF", r_val=tau_cut)

         nlcc = has_nlcc(kind_set)
         lsd = dft_control%lsd
         nspins = dft_control%nspins
         mspins = nspins
         IF (is_triplet) THEN
            CPASSERT(nspins == 1)
            lsd = .TRUE.
            mspins = 2
         END IF
         needs = xc_functionals_get_needs(xc_fun_section, lsd=lsd, calc_potential=.TRUE.)
         gradient_f = (needs%drho .OR. needs%drho_spin)
         tau_f = (needs%tau .OR. needs%tau_spin)

         ! Here starts the loop over all the atoms
         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
            CALL get_qs_kind(kind_set(ikind), paw_atom=paw_atom, &
                             harmonics=harmonics, grid_atom=grid_atom)
            CALL get_qs_kind(kind_set(ikind), basis_set=basis_1c, basis_type="GAPW_1C")

            IF (.NOT. paw_atom) CYCLE

            nr = grid_atom%nr
            na = grid_atom%ng_sphere

            ! set integration weights
            IF (accint) THEN
               weight_h => grid_atom%weight
               alpha = dft_control%qs_control%gapw_control%aw(ikind)
               IF (ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                  IF (grid_atom%gapw_weight_alpha /= alpha) DEALLOCATE (grid_atom%gapw_weight_s)
               END IF
               IF (.NOT. ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                  ALLOCATE (grid_atom%gapw_weight_s(na, nr))
                  DO ir = 1, nr
                     agr = 1.0_dp - EXP(-alpha*grid_atom%rad2(ir))
                     grid_atom%gapw_weight_s(:, ir) = grid_atom%weight(:, ir)*agr
                  END DO
                  grid_atom%gapw_weight_alpha = alpha
               END IF
               weight_s => grid_atom%gapw_weight_s
            ELSE
               weight_h => grid_atom%weight
               weight_s => grid_atom%weight
            END IF

            ! Prepare the structures needed to calculate and store the xc derivatives

            ! Array dimension: here anly one dimensional arrays are used,
            ! i.e. only the first column of deriv_data is read.
            ! The other to dimensions  are set to size equal 1
            bounds(1:2, 1:3) = 1
            bounds(2, 1) = na
            bounds(2, 2) = nr

            ! create a place where to put the derivatives
            CALL xc_dset_create(deriv_set, local_bounds=bounds)
            ! create the place where to store the argument for the functionals
            CALL xc_rho_set_create(rho_set_h, bounds, rho_cutoff=density_cut, &
                                   drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
            CALL xc_rho_set_create(rho_set_s, bounds, rho_cutoff=density_cut, &
                                   drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
            CALL xc_rho_set_create(rho1_set_h, bounds, rho_cutoff=density_cut, &
                                   drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
            CALL xc_rho_set_create(rho1_set_s, bounds, rho_cutoff=density_cut, &
                                   drho_cutoff=gradient_cut, tau_cutoff=tau_cut)

            ! allocate the required 3d arrays where to store rho and drho
            CALL xc_rho_set_atom_update(rho_set_h, needs, mspins, bounds)
            CALL xc_rho_set_atom_update(rho_set_s, needs, mspins, bounds)
            CALL xc_rho_set_atom_update(rho1_set_h, needs, mspins, bounds)
            CALL xc_rho_set_atom_update(rho1_set_s, needs, mspins, bounds)

            ALLOCATE (rho_h(na, nr, nspins), rho_s(na, nr, nspins), &
                      rho0_h(na, nr, nspins), rho0_s(na, nr, nspins), &
                      rho1_h(na, nr, nspins), rho1_s(na, nr, nspins))
            ALLOCATE (vxc_h(na, nr, nspins), vxc_s(na, nr, nspins))
            IF (gradient_f) THEN
               ALLOCATE (drho_h(4, na, nr, nspins), drho_s(4, na, nr, nspins), &
                         drho0_h(4, na, nr, nspins), drho0_s(4, na, nr, nspins), &
                         drho1_h(4, na, nr, nspins), drho1_s(4, na, nr, nspins))
               ALLOCATE (vxg_h(3, na, nr, nspins), vxg_s(3, na, nr, nspins))
            END IF
            IF (tau_f) THEN
               CALL create_tau_basis_cache(tau_basis_cache, grid_atom, basis_1c, harmonics)
               ALLOCATE (tau_h(na, nr, nspins), tau_s(na, nr, nspins), &
                         tau0_h(na, nr, nspins), tau0_s(na, nr, nspins), &
                         tau1_h(na, nr, nspins), tau1_s(na, nr, nspins))
               ALLOCATE (vtau_h(na, nr, nspins), vtau_s(na, nr, nspins))
            END IF
            !
            ! NLCC: prepare rho and drho of the core charge for this KIND
            donlcc = .FALSE.
            IF (nlcc) THEN
               NULLIFY (rho_nlcc)
               rho_nlcc => kind_set(ikind)%nlcc_pot
               IF (ASSOCIATED(rho_nlcc)) donlcc = .TRUE.
            END IF

            ! Distribute the atoms of this kind
            num_pe = para_env%num_pe
            bo = get_limit(natom, num_pe, para_env%mepos)

            DO iat = bo(1), bo(2)
               iatom = atom_list(iat)
               !
               NULLIFY (int_hh, int_ss)
               rho0_atom => rho0_atom_set(iatom)
               CALL get_rho_atom(rho_atom=rho0_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
               ALLOCATE (fint_ss(nspins), fint_hh(nspins))
               DO ns = 1, nspins
                  nf = SIZE(int_ss(ns)%r_coef, 1)
                  ALLOCATE (fint_ss(ns)%r_coef(nf, nf))
                  nf = SIZE(int_hh(ns)%r_coef, 1)
                  ALLOCATE (fint_hh(ns)%r_coef(nf, nf))
               END DO

               ! RHO0
               rho0_h = 0.0_dp
               rho0_s = 0.0_dp
               rho0_atom => rho0_atom_set(iatom)
               IF (gradient_f) THEN
                  NULLIFY (r_h, r_s, dr_h, dr_s, r_h_d, r_s_d)
                  CALL get_rho_atom(rho_atom=rho0_atom, rho_rad_h=r_h, rho_rad_s=r_s, drho_rad_h=dr_h, &
                                    drho_rad_s=dr_s, rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d)
                  drho0_h = 0.0_dp
                  drho0_s = 0.0_dp
               ELSE
                  NULLIFY (r_h, r_s)
                  CALL get_rho_atom(rho_atom=rho0_atom, rho_rad_h=r_h, rho_rad_s=r_s)
                  rho_d = 0.0_dp
               END IF
               DO ir = 1, nr
                  CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_f, &
                                        ir, r_h, r_s, rho0_h, rho0_s, dr_h, dr_s, &
                                        r_h_d, r_s_d, drho0_h, drho0_s)
                  IF (donlcc) THEN
                     CALL calc_rho_nlcc(grid_atom, nspins, gradient_f, &
                                        ir, rho_nlcc(:, 1), rho0_h, rho0_s, rho_nlcc(:, 2), drho0_h, drho0_s)
                  END IF
               END DO
               IF (tau_f) THEN
                  !compute tau on the grid all at once
                  CALL calc_tau_atom(tau0_h, tau0_s, rho0_atom, tau_basis_cache, nspins)
               ELSE
                  tau_d = 0.0_dp
               END IF
               ! RHO1
               rho1_h = 0.0_dp
               rho1_s = 0.0_dp
               rho1_atom => rho1_atom_set(iatom)
               IF (gradient_f) THEN
                  NULLIFY (r_h, r_s, dr_h, dr_s, r_h_d, r_s_d)
                  CALL get_rho_atom(rho_atom=rho1_atom, rho_rad_h=r_h, rho_rad_s=r_s, drho_rad_h=dr_h, &
                                    drho_rad_s=dr_s, rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d)
                  drho1_h = 0.0_dp
                  drho1_s = 0.0_dp
               ELSE
                  NULLIFY (r_h, r_s)
                  CALL get_rho_atom(rho_atom=rho1_atom, rho_rad_h=r_h, rho_rad_s=r_s)
               END IF
               DO ir = 1, nr
                  CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_f, &
                                        ir, r_h, r_s, rho1_h, rho1_s, dr_h, dr_s, &
                                        r_h_d, r_s_d, drho1_h, drho1_s)
               END DO
               IF (tau_f) THEN
                  !compute tau on the grid all at once
                  CALL calc_tau_atom(tau1_h, tau1_s, rho1_atom, tau_basis_cache, nspins)
               END IF

               DO ir = 1, nr
                  IF (tau_f) THEN
                     CALL fill_rho_set(rho1_set_h, lsd, nspins, needs, rho1_h, drho1_h, tau1_h, na, ir)
                     CALL fill_rho_set(rho1_set_s, lsd, nspins, needs, rho1_s, drho1_s, tau1_s, na, ir)
                  ELSE IF (gradient_f) THEN
                     CALL fill_rho_set(rho1_set_h, lsd, nspins, needs, rho1_h, drho1_h, tau_d, na, ir)
                     CALL fill_rho_set(rho1_set_s, lsd, nspins, needs, rho1_s, drho1_s, tau_d, na, ir)
                  ELSE
                     CALL fill_rho_set(rho1_set_h, lsd, nspins, needs, rho1_h, rho_d, tau_d, na, ir)
                     CALL fill_rho_set(rho1_set_s, lsd, nspins, needs, rho1_s, rho_d, tau_d, na, ir)
                  END IF
               END DO

               ! RHO2
               rho2_atom => rho2_atom_set(iatom)

               DO istep = -nstep, nstep

                  beta = REAL(istep, KIND=dp)*epsrho

                  rho_h = rho0_h + beta*rho1_h
                  rho_s = rho0_s + beta*rho1_s
                  IF (gradient_f) THEN
                     drho_h = drho0_h + beta*drho1_h
                     drho_s = drho0_s + beta*drho1_s
                  END IF
                  IF (tau_f) THEN
                     tau_h = tau0_h + beta*tau1_h
                     tau_s = tau0_s + beta*tau1_s
                  END IF
                  !
                  IF (gradient_f) THEN
                     drho_h(4, :, :, :) = SQRT( &
                                          drho_h(1, :, :, :)*drho_h(1, :, :, :) + &
                                          drho_h(2, :, :, :)*drho_h(2, :, :, :) + &
                                          drho_h(3, :, :, :)*drho_h(3, :, :, :))

                     drho_s(4, :, :, :) = SQRT( &
                                          drho_s(1, :, :, :)*drho_s(1, :, :, :) + &
                                          drho_s(2, :, :, :)*drho_s(2, :, :, :) + &
                                          drho_s(3, :, :, :)*drho_s(3, :, :, :))
                  END IF

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

                  ! hard atom density !
                  CALL xc_dset_zero_all(deriv_set)
                  CALL xc_2nd_deriv_of_r(xc_section=xc_section, &
                                         rho_set=rho_set_h, rho1_set=rho1_set_h, &
                                         deriv_set=deriv_set, &
                                         w=weight_h, vxc=vxc_h, vxg=vxg_h, vtau=vtau_h, &
                                         do_triplet=is_triplet)
                  ! soft atom density !
                  CALL xc_dset_zero_all(deriv_set)
                  CALL xc_2nd_deriv_of_r(xc_section=xc_section, &
                                         rho_set=rho_set_s, rho1_set=rho1_set_s, &
                                         deriv_set=deriv_set, &
                                         w=weight_s, vxc=vxc_s, vxg=vxg_s, vtau=vtau_s, &
                                         do_triplet=is_triplet)
                  ! potentials
                  DO ns = 1, nspins
                     fint_hh(ns)%r_coef(:, :) = 0.0_dp
                     fint_ss(ns)%r_coef(:, :) = 0.0_dp
                  END DO
                  IF (gradient_f) THEN
                     CALL gaVxcgb_GC(vxc_h, vxc_s, vxg_h, vxg_s, fint_hh, fint_ss, &
                                     grid_atom, basis_1c, harmonics, nspins)
                  ELSE
                     CALL gaVxcgb_noGC(vxc_h, vxc_s, fint_hh, fint_ss, &
                                       grid_atom, basis_1c, harmonics, nspins)
                  END IF
                  IF (tau_f) THEN
                     CALL dgaVtaudgb(vtau_h, vtau_s, fint_hh, fint_ss, &
                                     tau_basis_cache, nspins)
                  END IF
                  ! second derivative gxc
                  NULLIFY (int_hh, int_ss)
                  CALL get_rho_atom(rho_atom=rho2_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
                  DO ns = 1, nspins
                     int_ss(ns)%r_coef(:, :) = int_ss(ns)%r_coef(:, :) + oeps1*ak(istep)*fint_ss(ns)%r_coef(:, :)
                     int_hh(ns)%r_coef(:, :) = int_hh(ns)%r_coef(:, :) + oeps1*ak(istep)*fint_hh(ns)%r_coef(:, :)
                  END DO
               END DO
               !
               DO ns = 1, nspins
                  DEALLOCATE (fint_ss(ns)%r_coef)
                  DEALLOCATE (fint_hh(ns)%r_coef)
               END DO
               DEALLOCATE (fint_ss, fint_hh)

            END DO ! iat

            ! Release the xc structure used to store the xc derivatives
            CALL xc_dset_release(deriv_set)
            CALL xc_rho_set_release(rho_set_h)
            CALL xc_rho_set_release(rho_set_s)
            CALL xc_rho_set_release(rho1_set_h)
            CALL xc_rho_set_release(rho1_set_s)

            DEALLOCATE (rho_h, rho_s, rho0_h, rho0_s, rho1_h, rho1_s)
            DEALLOCATE (vxc_h, vxc_s)
            IF (gradient_f) THEN
               DEALLOCATE (drho_h, drho_s, drho0_h, drho0_s, drho1_h, drho1_s)
               DEALLOCATE (vxg_h, vxg_s)
            END IF
            IF (tau_f) THEN
               DEALLOCATE (tau_h, tau_s, tau0_h, tau0_s, tau1_h, tau1_s)
               DEALLOCATE (vtau_h, vtau_s)
               CALL release_tau_basis_cache(tau_basis_cache)
            END IF
         END DO ! ikind

      END IF !xc_none

      CALL timestop(handle)

   END SUBROUTINE gfxc_atom_diff

! **************************************************************************************************
!> \brief ...
!> \param grid_atom ...
!> \param harmonics ...
!> \param nspins ...
!> \param grad_func ...
!> \param ir ...
!> \param r_h ...
!> \param r_s ...
!> \param rho_h ...
!> \param rho_s ...
!> \param dr_h ...
!> \param dr_s ...
!> \param r_h_d ...
!> \param r_s_d ...
!> \param drho_h ...
!> \param drho_s ...
! **************************************************************************************************
   SUBROUTINE calc_rho_angular(grid_atom, harmonics, nspins, grad_func, &
                               ir, r_h, r_s, rho_h, rho_s, &
                               dr_h, dr_s, r_h_d, r_s_d, drho_h, drho_s)

      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      INTEGER, INTENT(IN)                                :: nspins
      LOGICAL, INTENT(IN)                                :: grad_func
      INTEGER, INTENT(IN)                                :: ir
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: r_h, r_s
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho_h, rho_s
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: dr_h, dr_s
      TYPE(rho_atom_coeff), DIMENSION(:, :), POINTER     :: r_h_d, r_s_d
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: drho_h, drho_s

      INTEGER                                            :: ia, iso, ispin, na
      REAL(KIND=dp)                                      :: rad, urad

      CPASSERT(ASSOCIATED(r_h))
      CPASSERT(ASSOCIATED(r_s))
      CPASSERT(ASSOCIATED(rho_h))
      CPASSERT(ASSOCIATED(rho_s))
      IF (grad_func) THEN
         CPASSERT(ASSOCIATED(dr_h))
         CPASSERT(ASSOCIATED(dr_s))
         CPASSERT(ASSOCIATED(r_h_d))
         CPASSERT(ASSOCIATED(r_s_d))
         CPASSERT(ASSOCIATED(drho_h))
         CPASSERT(ASSOCIATED(drho_s))
      END IF

      na = grid_atom%ng_sphere
      rad = grid_atom%rad(ir)
      urad = grid_atom%oorad2l(ir, 1)
      DO ispin = 1, nspins
         DO iso = 1, harmonics%max_iso_not0
            DO ia = 1, na
               rho_h(ia, ir, ispin) = rho_h(ia, ir, ispin) + &
                                      r_h(ispin)%r_coef(ir, iso)*harmonics%slm(ia, iso)
               rho_s(ia, ir, ispin) = rho_s(ia, ir, ispin) + &
                                      r_s(ispin)%r_coef(ir, iso)*harmonics%slm(ia, iso)
            END DO ! ia
         END DO ! iso
      END DO ! ispin

      IF (grad_func) THEN
         DO ispin = 1, nspins
            DO iso = 1, harmonics%max_iso_not0
               DO ia = 1, na

                  ! components of the gradient of rho1 hard
                  drho_h(1, ia, ir, ispin) = drho_h(1, ia, ir, ispin) + &
                                             dr_h(ispin)%r_coef(ir, iso)* &
                                             harmonics%a(1, ia)*harmonics%slm(ia, iso) + &
                                             r_h_d(1, ispin)%r_coef(ir, iso)* &
                                             harmonics%slm(ia, iso)

                  drho_h(2, ia, ir, ispin) = drho_h(2, ia, ir, ispin) + &
                                             dr_h(ispin)%r_coef(ir, iso)* &
                                             harmonics%a(2, ia)*harmonics%slm(ia, iso) + &
                                             r_h_d(2, ispin)%r_coef(ir, iso)* &
                                             harmonics%slm(ia, iso)

                  drho_h(3, ia, ir, ispin) = drho_h(3, ia, ir, ispin) + &
                                             dr_h(ispin)%r_coef(ir, iso)* &
                                             harmonics%a(3, ia)*harmonics%slm(ia, iso) + &
                                             r_h_d(3, ispin)%r_coef(ir, iso)* &
                                             harmonics%slm(ia, iso)

                  ! components of the gradient of rho1 soft
                  drho_s(1, ia, ir, ispin) = drho_s(1, ia, ir, ispin) + &
                                             dr_s(ispin)%r_coef(ir, iso)* &
                                             harmonics%a(1, ia)*harmonics%slm(ia, iso) + &
                                             r_s_d(1, ispin)%r_coef(ir, iso)* &
                                             harmonics%slm(ia, iso)

                  drho_s(2, ia, ir, ispin) = drho_s(2, ia, ir, ispin) + &
                                             dr_s(ispin)%r_coef(ir, iso)* &
                                             harmonics%a(2, ia)*harmonics%slm(ia, iso) + &
                                             r_s_d(2, ispin)%r_coef(ir, iso)* &
                                             harmonics%slm(ia, iso)

                  drho_s(3, ia, ir, ispin) = drho_s(3, ia, ir, ispin) + &
                                             dr_s(ispin)%r_coef(ir, iso)* &
                                             harmonics%a(3, ia)*harmonics%slm(ia, iso) + &
                                             r_s_d(3, ispin)%r_coef(ir, iso)* &
                                             harmonics%slm(ia, iso)

               END DO ! ia
            END DO ! iso
            DO ia = 1, na
               drho_h(4, ia, ir, ispin) = SQRT( &
                                          drho_h(1, ia, ir, ispin)*drho_h(1, ia, ir, ispin) + &
                                          drho_h(2, ia, ir, ispin)*drho_h(2, ia, ir, ispin) + &
                                          drho_h(3, ia, ir, ispin)*drho_h(3, ia, ir, ispin))

               drho_s(4, ia, ir, ispin) = SQRT( &
                                          drho_s(1, ia, ir, ispin)*drho_s(1, ia, ir, ispin) + &
                                          drho_s(2, ia, ir, ispin)*drho_s(2, ia, ir, ispin) + &
                                          drho_s(3, ia, ir, ispin)*drho_s(3, ia, ir, ispin))
            END DO ! ia
         END DO ! ispin
      END IF

   END SUBROUTINE calc_rho_angular

! **************************************************************************************************
!> \brief Precompute radial and angular factors for GAPW meta-GGA tau contractions
!> \param tau_cache precomputed compact one-center gradient basis
!> \param grid_atom atom-centered integration grid
!> \param basis_1c GAPW one-center basis
!> \param harmonics spherical harmonics on the atom-centered grid
! **************************************************************************************************
   SUBROUTINE create_tau_basis_cache(tau_cache, grid_atom, basis_1c, harmonics)

      TYPE(tau_basis_cache_type), INTENT(INOUT)          :: tau_cache
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      INTEGER                                            :: dir, ia, igrid, ip, ipgf, ir, iset, iso, &
                                                            l, starti
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: a1, a2, gexp, r1, r2
      REAL(dp), DIMENSION(:, :), POINTER                 :: slm
      REAL(dp), DIMENSION(:, :, :), POINTER              :: dslm_dxyz

      NULLIFY (slm, dslm_dxyz)

      CALL release_tau_basis_cache(tau_cache)

      CALL get_gto_basis_set(gto_basis_set=basis_1c, lmax=tau_cache%lmax, &
                             lmin=tau_cache%lmin, maxso=tau_cache%maxso, &
                             npgf=tau_cache%npgf, nset=tau_cache%nset, &
                             zet=tau_cache%zet)
      CALL get_paw_basis_info(basis_1c, o2nindex=tau_cache%o2nindex, &
                              n2oindex=tau_cache%n2oindex, &
                              nsatbas=tau_cache%nsatbas)

      tau_cache%nr = grid_atom%nr
      tau_cache%na = grid_atom%ng_sphere
      slm => harmonics%slm
      dslm_dxyz => harmonics%dslm_dxyz

      ALLOCATE (tau_cache%grad(tau_cache%na*tau_cache%nr, tau_cache%nsatbas, 3))
      ALLOCATE (a1(tau_cache%na), a2(tau_cache%na), gexp(tau_cache%nr), &
                r1(tau_cache%nr), r2(tau_cache%nr))
      tau_cache%grad = 0.0_dp

      DO iset = 1, tau_cache%nset
         DO ipgf = 1, tau_cache%npgf(iset)
            starti = (iset - 1)*tau_cache%maxso + &
                     (ipgf - 1)*nsoset(tau_cache%lmax(iset))
            gexp(1:tau_cache%nr) = EXP(-tau_cache%zet(ipgf, iset)* &
                                       grid_atom%rad2(1:tau_cache%nr))
            DO iso = nsoset(tau_cache%lmin(iset) - 1) + 1, nsoset(tau_cache%lmax(iset))
               ip = tau_cache%o2nindex(starti + iso)
               IF (ip == 0) CYCLE
               l = indso(1, iso)

               r1(1:tau_cache%nr) = grid_atom%rad(1:tau_cache%nr)**(l - 1)*gexp(1:tau_cache%nr)
               r2(1:tau_cache%nr) = -2.0_dp*tau_cache%zet(ipgf, iset)* &
                                    grid_atom%rad2(1:tau_cache%nr)*r1(1:tau_cache%nr)

               DO dir = 1, 3
                  a1(1:tau_cache%na) = dslm_dxyz(dir, 1:tau_cache%na, iso)
                  a2(1:tau_cache%na) = harmonics%a(dir, 1:tau_cache%na)*slm(1:tau_cache%na, iso)
                  DO ir = 1, tau_cache%nr
                     DO ia = 1, tau_cache%na
                        igrid = ia + (ir - 1)*tau_cache%na
                        tau_cache%grad(igrid, ip, dir) = r1(ir)*a1(ia) + r2(ir)*a2(ia)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      DEALLOCATE (a1, a2, gexp, r1, r2)

   END SUBROUTINE create_tau_basis_cache

! **************************************************************************************************
!> \brief Release precomputed GAPW meta-GGA tau factors
!> \param tau_cache precomputed compact one-center gradient basis
! **************************************************************************************************
   SUBROUTINE release_tau_basis_cache(tau_cache)

      TYPE(tau_basis_cache_type), INTENT(INOUT)          :: tau_cache

      IF (ALLOCATED(tau_cache%grad)) DEALLOCATE (tau_cache%grad)
      IF (ASSOCIATED(tau_cache%n2oindex)) DEALLOCATE (tau_cache%n2oindex)
      IF (ASSOCIATED(tau_cache%o2nindex)) DEALLOCATE (tau_cache%o2nindex)
      NULLIFY (tau_cache%lmax, tau_cache%lmin, tau_cache%n2oindex, tau_cache%npgf, &
               tau_cache%zet, tau_cache%o2nindex)
      tau_cache%maxso = 0
      tau_cache%na = 0
      tau_cache%nr = 0
      tau_cache%nsatbas = 0
      tau_cache%nset = 0

   END SUBROUTINE release_tau_basis_cache

! **************************************************************************************************
!> \brief Computes tau hard and soft on the atomic grids for meta-GGA calculations
!> \param tau_h the hard part of tau
!> \param tau_s the soft part of tau
!> \param rho_atom atom-centered density matrices
!> \param tau_cache precomputed compact one-center gradient basis
!> \param nspins number of spin channels
!> \note This is a rewrite to correct a meta-GGA GAPW bug. This is more brute force than the original,
!>       which was done along in qs_rho_atom_methods.F, but makes sure that no corner is cut in
!>       terms of accuracy (A. Bussy)
! **************************************************************************************************
   SUBROUTINE calc_tau_atom(tau_h, tau_s, rho_atom, tau_cache, nspins)

      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: tau_h, tau_s
      TYPE(rho_atom_type), POINTER                       :: rho_atom
      TYPE(tau_basis_cache_type), INTENT(IN)             :: tau_cache
      INTEGER, INTENT(IN)                                :: nspins

      CHARACTER(len=*), PARAMETER                        :: routineN = 'calc_tau_atom'

      INTEGER                                            :: dir, handle, ia, ibas, igrid, ir, ispin, &
                                                            na, nbas, ngrid, nr
      REAL(dp)                                           :: tau_sum
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: work

      CALL timeset(routineN, handle)

      CPASSERT(ALLOCATED(tau_cache%grad))

      !zeroing tau, assuming it is already allocated
      tau_h = 0.0_dp
      tau_s = 0.0_dp

      nr = tau_cache%nr
      na = tau_cache%na
      nbas = tau_cache%nsatbas
      ngrid = na*nr
      ALLOCATE (work(ngrid, nbas))

      DO ispin = 1, nspins
         DO dir = 1, 3
            CALL dgemm('N', 'T', ngrid, nbas, nbas, 0.5_dp, tau_cache%grad(:, :, dir), &
                       ngrid, rho_atom%cpc_h(ispin)%r_coef, nbas, 0.0_dp, work, ngrid)
!$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(2) SCHEDULE(STATIC) &
!$OMP SHARED(dir, ispin, na, nbas, nr, tau_cache, tau_h, work) &
!$OMP PRIVATE(ia, ibas, igrid, ir, tau_sum)
            DO ir = 1, nr
               DO ia = 1, na
                  igrid = ia + (ir - 1)*na
                  tau_sum = 0.0_dp
                  DO ibas = 1, nbas
                     tau_sum = tau_sum + tau_cache%grad(igrid, ibas, dir)*work(igrid, ibas)
                  END DO
                  tau_h(ia, ir, ispin) = tau_h(ia, ir, ispin) + tau_sum
               END DO
            END DO
!$OMP END PARALLEL DO

            CALL dgemm('N', 'T', ngrid, nbas, nbas, 0.5_dp, tau_cache%grad(:, :, dir), &
                       ngrid, rho_atom%cpc_s(ispin)%r_coef, nbas, 0.0_dp, work, ngrid)
!$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(2) SCHEDULE(STATIC) &
!$OMP SHARED(dir, ispin, na, nbas, nr, tau_cache, tau_s, work) &
!$OMP PRIVATE(ia, ibas, igrid, ir, tau_sum)
            DO ir = 1, nr
               DO ia = 1, na
                  igrid = ia + (ir - 1)*na
                  tau_sum = 0.0_dp
                  DO ibas = 1, nbas
                     tau_sum = tau_sum + tau_cache%grad(igrid, ibas, dir)*work(igrid, ibas)
                  END DO
                  tau_s(ia, ir, ispin) = tau_s(ia, ir, ispin) + tau_sum
               END DO
            END DO
!$OMP END PARALLEL DO
         END DO
      END DO

      DEALLOCATE (work)

      CALL timestop(handle)

   END SUBROUTINE calc_tau_atom

! **************************************************************************************************
!> \brief ...
!> \param grid_atom ...
!> \param nspins ...
!> \param grad_func ...
!> \param ir ...
!> \param rho_nlcc ...
!> \param rho_h ...
!> \param rho_s ...
!> \param drho_nlcc ...
!> \param drho_h ...
!> \param drho_s ...
! **************************************************************************************************
   SUBROUTINE calc_rho_nlcc(grid_atom, nspins, grad_func, &
                            ir, rho_nlcc, rho_h, rho_s, drho_nlcc, drho_h, drho_s)

      TYPE(grid_atom_type), POINTER                      :: grid_atom
      INTEGER, INTENT(IN)                                :: nspins
      LOGICAL, INTENT(IN)                                :: grad_func
      INTEGER, INTENT(IN)                                :: ir
      REAL(KIND=dp), DIMENSION(:)                        :: rho_nlcc
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho_h, rho_s
      REAL(KIND=dp), DIMENSION(:)                        :: drho_nlcc
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: drho_h, drho_s

      INTEGER                                            :: ia, ispin, na
      REAL(KIND=dp)                                      :: drho, dx, dy, dz, rad, rho, urad, xsp

      CPASSERT(ASSOCIATED(rho_h))
      CPASSERT(ASSOCIATED(rho_s))
      IF (grad_func) THEN
         CPASSERT(ASSOCIATED(drho_h))
         CPASSERT(ASSOCIATED(drho_s))
      END IF

      na = grid_atom%ng_sphere
      rad = grid_atom%rad(ir)
      urad = grid_atom%oorad2l(ir, 1)

      xsp = REAL(nspins, KIND=dp)
      rho = rho_nlcc(ir)/xsp
      DO ispin = 1, nspins
         rho_h(1:na, ir, ispin) = rho_h(1:na, ir, ispin) + rho
         rho_s(1:na, ir, ispin) = rho_s(1:na, ir, ispin) + rho
      END DO ! ispin

      IF (grad_func) THEN
         drho = drho_nlcc(ir)/xsp
         DO ispin = 1, nspins
            DO ia = 1, na
               IF (grid_atom%azi(ia) == 0.0_dp) THEN
                  dx = 0.0_dp
                  dy = 0.0_dp
               ELSE
                  dx = grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia)
                  dy = grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia)
               END IF
               dz = grid_atom%cos_pol(ia)
               ! components of the gradient of rho1 hard
               drho_h(1, ia, ir, ispin) = drho_h(1, ia, ir, ispin) + drho*dx
               drho_h(2, ia, ir, ispin) = drho_h(2, ia, ir, ispin) + drho*dy
               drho_h(3, ia, ir, ispin) = drho_h(3, ia, ir, ispin) + drho*dz
               ! components of the gradient of rho1 soft
               drho_s(1, ia, ir, ispin) = drho_s(1, ia, ir, ispin) + drho*dx
               drho_s(2, ia, ir, ispin) = drho_s(2, ia, ir, ispin) + drho*dy
               drho_s(3, ia, ir, ispin) = drho_s(3, ia, ir, ispin) + drho*dz
               ! norm of gradient
               drho_h(4, ia, ir, ispin) = SQRT( &
                                          drho_h(1, ia, ir, ispin)*drho_h(1, ia, ir, ispin) + &
                                          drho_h(2, ia, ir, ispin)*drho_h(2, ia, ir, ispin) + &
                                          drho_h(3, ia, ir, ispin)*drho_h(3, ia, ir, ispin))

               drho_s(4, ia, ir, ispin) = SQRT( &
                                          drho_s(1, ia, ir, ispin)*drho_s(1, ia, ir, ispin) + &
                                          drho_s(2, ia, ir, ispin)*drho_s(2, ia, ir, ispin) + &
                                          drho_s(3, ia, ir, ispin)*drho_s(3, ia, ir, ispin))
            END DO ! ia
         END DO ! ispin
      END IF

   END SUBROUTINE calc_rho_nlcc

! **************************************************************************************************
!> \brief ...
!> \param vxc_h ...
!> \param vxc_s ...
!> \param int_hh ...
!> \param int_ss ...
!> \param grid_atom ...
!> \param basis_1c ...
!> \param harmonics ...
!> \param nspins ...
! **************************************************************************************************
   SUBROUTINE gaVxcgb_noGC(vxc_h, vxc_s, int_hh, int_ss, grid_atom, basis_1c, harmonics, nspins)

      REAL(dp), DIMENSION(:, :, :), POINTER              :: vxc_h, vxc_s
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: int_hh, int_ss
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      INTEGER, INTENT(IN)                                :: nspins

      CHARACTER(len=*), PARAMETER                        :: routineN = 'gaVxcgb_noGC'

      INTEGER :: handle, ia, ic, icg, ipgf1, ipgf2, ir, iset1, iset2, iso, iso1, iso2, ispin, l, &
         ld, lmax12, lmax_expansion, lmin12, m1, m2, max_iso_not0, max_iso_not0_local, max_s_harm, &
         maxl, maxso, n1, n2, na, ngau1, ngau2, nngau1, nr, nset, size1
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cg_n_list
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cg_list
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: g1, g2
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: gg, gVg_h, gVg_s, matso_h, matso_s, vx
      REAL(dp), DIMENSION(:, :), POINTER                 :: zet
      REAL(dp), DIMENSION(:, :, :), POINTER              :: my_CG

      CALL timeset(routineN, handle)

      NULLIFY (lmin, lmax, npgf, zet, my_CG)

      CALL get_gto_basis_set(gto_basis_set=basis_1c, lmax=lmax, lmin=lmin, &
                             maxso=maxso, maxl=maxl, npgf=npgf, &
                             nset=nset, zet=zet)

      nr = grid_atom%nr
      na = grid_atom%ng_sphere
      my_CG => harmonics%my_CG
      max_iso_not0 = harmonics%max_iso_not0
      lmax_expansion = indso(1, max_iso_not0)
      max_s_harm = harmonics%max_s_harm

      ALLOCATE (g1(nr), g2(nr), gg(nr, 0:2*maxl))
      ALLOCATE (gVg_h(na, 0:2*maxl), gVg_s(na, 0:2*maxl))
      ALLOCATE (matso_h(nsoset(maxl), nsoset(maxl)), &
                matso_s(nsoset(maxl), nsoset(maxl)))
      ALLOCATE (vx(na, nr))
      ALLOCATE (cg_list(2, nsoset(maxl)**2, max_s_harm), cg_n_list(max_s_harm))

      g1 = 0.0_dp
      g2 = 0.0_dp
      m1 = 0
      DO iset1 = 1, nset
         n1 = nsoset(lmax(iset1))
         m2 = 0
         DO iset2 = 1, nset
            CALL get_none0_cg_list(my_CG, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                   max_s_harm, lmax_expansion, cg_list, cg_n_list, max_iso_not0_local)
            CPASSERT(max_iso_not0_local <= max_iso_not0)

            n2 = nsoset(lmax(iset2))
            DO ipgf1 = 1, npgf(iset1)
               ngau1 = n1*(ipgf1 - 1) + m1
               size1 = nsoset(lmax(iset1)) - nsoset(lmin(iset1) - 1)
               nngau1 = nsoset(lmin(iset1) - 1) + ngau1

               g1(1:nr) = EXP(-zet(ipgf1, iset1)*grid_atom%rad2(1:nr))
               DO ipgf2 = 1, npgf(iset2)
                  ngau2 = n2*(ipgf2 - 1) + m2

                  g2(1:nr) = EXP(-zet(ipgf2, iset2)*grid_atom%rad2(1:nr))
                  lmin12 = lmin(iset1) + lmin(iset2)
                  lmax12 = lmax(iset1) + lmax(iset2)

                  ! reduce expansion local densities
                  IF (lmin12 <= lmax_expansion) THEN

                     gg = 0.0_dp
                     IF (lmin12 == 0) THEN
                        gg(1:nr, lmin12) = g1(1:nr)*g2(1:nr)
                     ELSE
                        gg(1:nr, lmin12) = grid_atom%rad2l(1:nr, lmin12)*g1(1:nr)*g2(1:nr)
                     END IF

                     ! limit the expansion of the local densities to a max L
                     IF (lmax12 > lmax_expansion) lmax12 = lmax_expansion

                     DO l = lmin12 + 1, lmax12
                        gg(1:nr, l) = grid_atom%rad(1:nr)*gg(:, l - 1)
                     END DO

                     DO ispin = 1, nspins
                        ld = lmax12 + 1
                        DO ir = 1, nr
                           vx(1:na, ir) = vxc_h(1:na, ir, ispin)
                        END DO
                        CALL dgemm('N', 'N', na, ld, nr, 1.0_dp, vx(1:na, 1:nr), na, &
                                   gg(1:nr, 0:lmax12), nr, 0.0_dp, gVg_h(1:na, 0:lmax12), na)
                        DO ir = 1, nr
                           vx(1:na, ir) = vxc_s(1:na, ir, ispin)
                        END DO
                        CALL dgemm('N', 'N', na, ld, nr, 1.0_dp, vx(1:na, 1:nr), na, &
                                   gg(1:nr, 0:lmax12), nr, 0.0_dp, gVg_s(1:na, 0:lmax12), na)

                        matso_h = 0.0_dp
                        matso_s = 0.0_dp
                        DO iso = 1, max_iso_not0_local
                           DO icg = 1, cg_n_list(iso)
                              iso1 = cg_list(1, icg, iso)
                              iso2 = cg_list(2, icg, iso)
                              l = indso(1, iso1) + indso(1, iso2)

                              CPASSERT(l <= lmax_expansion)
                              DO ia = 1, na
                                 matso_h(iso1, iso2) = matso_h(iso1, iso2) + &
                                                       gVg_h(ia, l)* &
                                                       my_CG(iso1, iso2, iso)* &
                                                       harmonics%slm(ia, iso)
                                 matso_s(iso1, iso2) = matso_s(iso1, iso2) + &
                                                       gVg_s(ia, l)* &
                                                       my_CG(iso1, iso2, iso)* &
                                                       harmonics%slm(ia, iso)
                              END DO
                           END DO
                        END DO

                        ! Write in the global matrix
                        DO ic = nsoset(lmin(iset2) - 1) + 1, nsoset(lmax(iset2))
                           iso1 = nsoset(lmin(iset1) - 1) + 1
                           iso2 = ngau2 + ic
                           CALL daxpy(size1, 1.0_dp, matso_h(iso1, ic), 1, &
                                      int_hh(ispin)%r_coef(nngau1 + 1, iso2), 1)
                           CALL daxpy(size1, 1.0_dp, matso_s(iso1, ic), 1, &
                                      int_ss(ispin)%r_coef(nngau1 + 1, iso2), 1)
                        END DO

                     END DO ! ispin

                  END IF ! lmax_expansion

               END DO ! ipfg2
            END DO ! ipfg1
            m2 = m2 + maxso
         END DO ! iset2
         m1 = m1 + maxso
      END DO ! iset1

      DEALLOCATE (g1, g2, gg, matso_h, matso_s, gVg_s, gVg_h, vx)

      DEALLOCATE (cg_list, cg_n_list)

      CALL timestop(handle)

   END SUBROUTINE gaVxcgb_noGC

! **************************************************************************************************
!> \brief ...
!> \param vxc_h ...
!> \param vxc_s ...
!> \param vxg_h ...
!> \param vxg_s ...
!> \param int_hh ...
!> \param int_ss ...
!> \param grid_atom ...
!> \param basis_1c ...
!> \param harmonics ...
!> \param nspins ...
! **************************************************************************************************
   SUBROUTINE gaVxcgb_GC(vxc_h, vxc_s, vxg_h, vxg_s, int_hh, int_ss, &
                         grid_atom, basis_1c, harmonics, nspins)

      REAL(dp), DIMENSION(:, :, :), POINTER              :: vxc_h, vxc_s
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: vxg_h, vxg_s
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: int_hh, int_ss
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      INTEGER, INTENT(IN)                                :: nspins

      CHARACTER(len=*), PARAMETER                        :: routineN = 'gaVxcgb_GC'

      INTEGER :: dmax_iso_not0_local, handle, ia, ic, icg, ipgf1, ipgf2, ir, iset1, iset2, iso, &
         iso1, iso2, ispin, l, lmax12, lmax_expansion, lmin12, m1, m2, max_iso_not0, &
         max_iso_not0_local, max_s_harm, maxl, maxso, n1, n2, na, ngau1, ngau2, nngau1, nr, nset, &
         size1
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cg_n_list, dcg_n_list
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cg_list, dcg_list
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf
      REAL(dp)                                           :: urad
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: g1, g2
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: dgg, gg, gVXCg_h, gVXCg_s, matso_h, &
                                                            matso_s
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: gVXGg_h, gVXGg_s
      REAL(dp), DIMENSION(:, :), POINTER                 :: zet
      REAL(dp), DIMENSION(:, :, :), POINTER              :: my_CG
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: my_CG_dxyz

      CALL timeset(routineN, handle)

      NULLIFY (lmin, lmax, npgf, zet, my_CG, my_CG_dxyz)

      CALL get_gto_basis_set(gto_basis_set=basis_1c, lmax=lmax, lmin=lmin, &
                             maxso=maxso, maxl=maxl, npgf=npgf, &
                             nset=nset, zet=zet)

      nr = grid_atom%nr
      na = grid_atom%ng_sphere
      my_CG => harmonics%my_CG
      my_CG_dxyz => harmonics%my_CG_dxyz
      max_iso_not0 = harmonics%max_iso_not0
      lmax_expansion = indso(1, max_iso_not0)
      max_s_harm = harmonics%max_s_harm

      ALLOCATE (g1(nr), g2(nr), gg(nr, 0:2*maxl), dgg(nr, 0:2*maxl))
      ALLOCATE (gVXCg_h(na, 0:2*maxl), gVXCg_s(na, 0:2*maxl))
      ALLOCATE (gVXGg_h(3, na, 0:2*maxl), gVXGg_s(3, na, 0:2*maxl))
      ALLOCATE (cg_list(2, nsoset(maxl)**2, max_s_harm), cg_n_list(max_s_harm), &
                dcg_list(2, nsoset(maxl)**2, max_s_harm), dcg_n_list(max_s_harm))

      ALLOCATE (matso_h(nsoset(maxl), nsoset(maxl)), &
                matso_s(nsoset(maxl), nsoset(maxl)))

      DO ispin = 1, nspins

         g1 = 0.0_dp
         g2 = 0.0_dp
         m1 = 0
         DO iset1 = 1, nset
            n1 = nsoset(lmax(iset1))
            m2 = 0
            DO iset2 = 1, nset
               CALL get_none0_cg_list(my_CG, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                      max_s_harm, lmax_expansion, cg_list, cg_n_list, max_iso_not0_local)
               CPASSERT(max_iso_not0_local <= max_iso_not0)
               CALL get_none0_cg_list(my_CG_dxyz, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                      max_s_harm, lmax_expansion, dcg_list, dcg_n_list, dmax_iso_not0_local)

               n2 = nsoset(lmax(iset2))
               DO ipgf1 = 1, npgf(iset1)
                  ngau1 = n1*(ipgf1 - 1) + m1
                  size1 = nsoset(lmax(iset1)) - nsoset(lmin(iset1) - 1)
                  nngau1 = nsoset(lmin(iset1) - 1) + ngau1

                  g1(1:nr) = EXP(-zet(ipgf1, iset1)*grid_atom%rad2(1:nr))
                  DO ipgf2 = 1, npgf(iset2)
                     ngau2 = n2*(ipgf2 - 1) + m2

                     g2(1:nr) = EXP(-zet(ipgf2, iset2)*grid_atom%rad2(1:nr))
                     lmin12 = lmin(iset1) + lmin(iset2)
                     lmax12 = lmax(iset1) + lmax(iset2)

                     !test reduce expansion local densities
                     IF (lmin12 <= lmax_expansion) THEN

                        gg = 0.0_dp
                        dgg = 0.0_dp

                        IF (lmin12 == 0) THEN
                           gg(1:nr, lmin12) = g1(1:nr)*g2(1:nr)
                        ELSE
                           gg(1:nr, lmin12) = grid_atom%rad2l(1:nr, lmin12)*g1(1:nr)*g2(1:nr)
                        END IF

                        !test reduce expansion local densities
                        IF (lmax12 > lmax_expansion) lmax12 = lmax_expansion

                        DO l = lmin12 + 1, lmax12
                           gg(1:nr, l) = grid_atom%rad(1:nr)*gg(:, l - 1)
                           dgg(1:nr, l - 1) = dgg(1:nr, l - 1) - 2.0_dp*(zet(ipgf1, iset1) + &
                                                                         zet(ipgf2, iset2))*gg(1:nr, l)
                        END DO
                        dgg(1:nr, lmax12) = dgg(1:nr, lmax12) - 2.0_dp*(zet(ipgf1, iset1) + &
                                                                        zet(ipgf2, iset2))*grid_atom%rad(1:nr)* &
                                            gg(1:nr, lmax12)

                        gVXCg_h = 0.0_dp
                        gVXCg_s = 0.0_dp
                        gVXGg_h = 0.0_dp
                        gVXGg_s = 0.0_dp

                        ! Cross Term
                        DO l = lmin12, lmax12
                           DO ia = 1, na
                              DO ir = 1, nr
                                 gVXCg_h(ia, l) = gVXCg_h(ia, l) + &
                                                  gg(ir, l)*vxc_h(ia, ir, ispin) + &
                                                  dgg(ir, l)* &
                                                  (vxg_h(1, ia, ir, ispin)*harmonics%a(1, ia) + &
                                                   vxg_h(2, ia, ir, ispin)*harmonics%a(2, ia) + &
                                                   vxg_h(3, ia, ir, ispin)*harmonics%a(3, ia))

                                 gVXCg_s(ia, l) = gVXCg_s(ia, l) + &
                                                  gg(ir, l)*vxc_s(ia, ir, ispin) + &
                                                  dgg(ir, l)* &
                                                  (vxg_s(1, ia, ir, ispin)*harmonics%a(1, ia) + &
                                                   vxg_s(2, ia, ir, ispin)*harmonics%a(2, ia) + &
                                                   vxg_s(3, ia, ir, ispin)*harmonics%a(3, ia))

                                 urad = grid_atom%oorad2l(ir, 1)

                                 gVXGg_h(1, ia, l) = gVXGg_h(1, ia, l) + &
                                                     vxg_h(1, ia, ir, ispin)* &
                                                     gg(ir, l)*urad

                                 gVXGg_h(2, ia, l) = gVXGg_h(2, ia, l) + &
                                                     vxg_h(2, ia, ir, ispin)* &
                                                     gg(ir, l)*urad

                                 gVXGg_h(3, ia, l) = gVXGg_h(3, ia, l) + &
                                                     vxg_h(3, ia, ir, ispin)* &
                                                     gg(ir, l)*urad

                                 gVXGg_s(1, ia, l) = gVXGg_s(1, ia, l) + &
                                                     vxg_s(1, ia, ir, ispin)* &
                                                     gg(ir, l)*urad

                                 gVXGg_s(2, ia, l) = gVXGg_s(2, ia, l) + &
                                                     vxg_s(2, ia, ir, ispin)* &
                                                     gg(ir, l)*urad

                                 gVXGg_s(3, ia, l) = gVXGg_s(3, ia, l) + &
                                                     vxg_s(3, ia, ir, ispin)* &
                                                     gg(ir, l)*urad

                              END DO ! ir
                           END DO ! ia
                        END DO ! l

                        matso_h = 0.0_dp
                        matso_s = 0.0_dp
                        DO iso = 1, max_iso_not0_local
                           DO icg = 1, cg_n_list(iso)
                              iso1 = cg_list(1, icg, iso)
                              iso2 = cg_list(2, icg, iso)

                              l = indso(1, iso1) + indso(1, iso2)

                              !test reduce expansion local densities
                              CPASSERT(l <= lmax_expansion)
                              DO ia = 1, na
                                 matso_h(iso1, iso2) = matso_h(iso1, iso2) + &
                                                       gVXCg_h(ia, l)* &
                                                       harmonics%slm(ia, iso)* &
                                                       my_CG(iso1, iso2, iso)
                                 matso_s(iso1, iso2) = matso_s(iso1, iso2) + &
                                                       gVXCg_s(ia, l)* &
                                                       harmonics%slm(ia, iso)* &
                                                       my_CG(iso1, iso2, iso)
                              END DO ! ia

                              !test reduce expansion local densities

                           END DO

                        END DO ! iso

                        DO iso = 1, dmax_iso_not0_local
                           DO icg = 1, dcg_n_list(iso)
                              iso1 = dcg_list(1, icg, iso)
                              iso2 = dcg_list(2, icg, iso)

                              l = indso(1, iso1) + indso(1, iso2)
                              !test reduce expansion local densities
                              CPASSERT(l <= lmax_expansion)
                              DO ia = 1, na
                                 matso_h(iso1, iso2) = matso_h(iso1, iso2) + &
                                                       (gVXGg_h(1, ia, l)*my_CG_dxyz(1, iso1, iso2, iso) + &
                                                        gVXGg_h(2, ia, l)*my_CG_dxyz(2, iso1, iso2, iso) + &
                                                        gVXGg_h(3, ia, l)*my_CG_dxyz(3, iso1, iso2, iso))* &
                                                       harmonics%slm(ia, iso)

                                 matso_s(iso1, iso2) = matso_s(iso1, iso2) + &
                                                       (gVXGg_s(1, ia, l)*my_CG_dxyz(1, iso1, iso2, iso) + &
                                                        gVXGg_s(2, ia, l)*my_CG_dxyz(2, iso1, iso2, iso) + &
                                                        gVXGg_s(3, ia, l)*my_CG_dxyz(3, iso1, iso2, iso))* &
                                                       harmonics%slm(ia, iso)

                              END DO ! ia

                              !test reduce expansion local densities

                           END DO ! icg
                        END DO ! iso
                        !test reduce expansion local densities
                     END IF ! lmax_expansion

                     !  Write in the global matrix
                     DO ic = nsoset(lmin(iset2) - 1) + 1, nsoset(lmax(iset2))
                        iso1 = nsoset(lmin(iset1) - 1) + 1
                        iso2 = ngau2 + ic
                        CALL daxpy(size1, 1.0_dp, matso_h(iso1, ic), 1, &
                                   int_hh(ispin)%r_coef(nngau1 + 1, iso2), 1)
                        CALL daxpy(size1, 1.0_dp, matso_s(iso1, ic), 1, &
                                   int_ss(ispin)%r_coef(nngau1 + 1, iso2), 1)
                     END DO

                  END DO ! ipfg2
               END DO ! ipfg1
               m2 = m2 + maxso
            END DO ! iset2
            m1 = m1 + maxso
         END DO ! iset1
      END DO ! ispin

      DEALLOCATE (g1, g2, gg, dgg, matso_h, matso_s, gVXCg_h, gVXCg_s, gVXGg_h, gVXGg_s)
      DEALLOCATE (cg_list, cg_n_list, dcg_list, dcg_n_list)

      CALL timestop(handle)

   END SUBROUTINE gaVxcgb_GC

! **************************************************************************************************
!> \brief Integrates 0.5 * grad_ga .dot. (V_tau * grad_gb) on the atomic grid for meta-GGA
!> \param vtau_h the hard tau potential
!> \param vtau_s the soft tau potential
!> \param int_hh hard one-center matrix contribution
!> \param int_ss soft one-center matrix contribution
!> \param tau_cache precomputed compact one-center gradient basis
!> \param nspins number of spin channels
!> \note This is a rewrite to correct meta-GGA GAPW bug. This is more brute force than the original
!>       but makes sure that no corner is cut in terms of accuracy (A. Bussy)
! **************************************************************************************************
   SUBROUTINE dgaVtaudgb(vtau_h, vtau_s, int_hh, int_ss, tau_cache, nspins)

      REAL(dp), DIMENSION(:, :, :), POINTER              :: vtau_h, vtau_s
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: int_hh, int_ss
      TYPE(tau_basis_cache_type), INTENT(IN)             :: tau_cache
      INTEGER, INTENT(IN)                                :: nspins

      CHARACTER(len=*), PARAMETER                        :: routineN = 'dgaVtaudgb'

      INTEGER                                            :: dir, handle, ia, ibas, igrid, iold, ir, &
                                                            ispin, jbas, jold, max_old_basis, na, &
                                                            nbas, ngrid, nr
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: int_h, int_s, weighted_grad

      CALL timeset(routineN, handle)

      CPASSERT(ALLOCATED(tau_cache%grad))
      CPASSERT(ASSOCIATED(tau_cache%n2oindex))

      nr = tau_cache%nr
      na = tau_cache%na
      nbas = tau_cache%nsatbas
      ngrid = na*nr
      max_old_basis = MAXVAL(tau_cache%n2oindex)
      ALLOCATE (int_h(nbas, nbas), int_s(nbas, nbas), weighted_grad(ngrid, nbas))

      DO ispin = 1, nspins
         CPASSERT(SIZE(int_hh(ispin)%r_coef, 1) >= max_old_basis)
         CPASSERT(SIZE(int_hh(ispin)%r_coef, 2) >= max_old_basis)
         CPASSERT(SIZE(int_ss(ispin)%r_coef, 1) >= max_old_basis)
         CPASSERT(SIZE(int_ss(ispin)%r_coef, 2) >= max_old_basis)
         int_h = 0.0_dp
         int_s = 0.0_dp
         DO dir = 1, 3
!$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(3) SCHEDULE(STATIC) &
!$OMP SHARED(dir, ispin, na, nbas, nr, tau_cache, vtau_h, weighted_grad) &
!$OMP PRIVATE(ia, ibas, igrid, ir)
            DO ibas = 1, nbas
               DO ir = 1, nr
                  DO ia = 1, na
                     igrid = ia + (ir - 1)*na
                     weighted_grad(igrid, ibas) = vtau_h(ia, ir, ispin)* &
                                                  tau_cache%grad(igrid, ibas, dir)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            CALL dgemm('T', 'N', nbas, nbas, ngrid, 0.5_dp, tau_cache%grad(:, :, dir), &
                       ngrid, weighted_grad, ngrid, 1.0_dp, int_h, nbas)

!$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(3) SCHEDULE(STATIC) &
!$OMP SHARED(dir, ispin, na, nbas, nr, tau_cache, vtau_s, weighted_grad) &
!$OMP PRIVATE(ia, ibas, igrid, ir)
            DO ibas = 1, nbas
               DO ir = 1, nr
                  DO ia = 1, na
                     igrid = ia + (ir - 1)*na
                     weighted_grad(igrid, ibas) = vtau_s(ia, ir, ispin)* &
                                                  tau_cache%grad(igrid, ibas, dir)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            CALL dgemm('T', 'N', nbas, nbas, ngrid, 0.5_dp, tau_cache%grad(:, :, dir), &
                       ngrid, weighted_grad, ngrid, 1.0_dp, int_s, nbas)
         END DO

!$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(2) SCHEDULE(STATIC) &
!$OMP SHARED(int_h, int_hh, int_s, int_ss, ispin, nbas, tau_cache) &
!$OMP PRIVATE(ibas, iold, jbas, jold)
         DO jbas = 1, nbas
            DO ibas = 1, nbas
               jold = tau_cache%n2oindex(jbas)
               iold = tau_cache%n2oindex(ibas)
               int_hh(ispin)%r_coef(iold, jold) = int_hh(ispin)%r_coef(iold, jold) + &
                                                  int_h(ibas, jbas)
               int_ss(ispin)%r_coef(iold, jold) = int_ss(ispin)%r_coef(iold, jold) + &
                                                  int_s(ibas, jbas)
            END DO
         END DO
!$OMP END PARALLEL DO
      END DO

      DEALLOCATE (int_h, int_s, weighted_grad)

      CALL timestop(handle)

   END SUBROUTINE dgaVtaudgb

END MODULE qs_vxc_atom

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Integration grid, AO images and radial tables for the periodic RI-RS integrals d_lP(R).
!> \author Jan Wilhelm
! **************************************************************************************************
MODULE ri_rs_kp_types
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_rs_kp_types'

   PUBLIC :: ao_image_type, block_class_cache_type, l_block_partition_type, &
             radial_interpolation_type, uniform_settings_type

   TYPE uniform_settings_type
      REAL(dp)                          :: accuracy = 0.0_dp, &
                                           spacing = 0.0_dp, padding = 0.0_dp, &
                                           l_block_edge = 0.0_dp
      INTEGER                           :: counts(3) = 0, subdivisions(3) = 0
      REAL(dp)                          :: lower(3) = 0.0_dp, &
                                           origin(3) = 0.0_dp, &
                                           step(3, 3) = 0.0_dp, &
                                           weight = 0.0_dp
   END TYPE uniform_settings_type

   TYPE l_block_partition_type
      INTEGER                           :: nblocks = 0
      INTEGER, ALLOCATABLE              :: first(:), indices(:)
      REAL(dp), ALLOCATABLE             :: lower(:, :), upper(:, :)
   END TYPE l_block_partition_type

   ! Values of atom-centred functions at r - R_A - hT on the points of one class of cell-aligned
   ! quadrature blocks: the RI potentials V_P (with point bounds) or the AO values phi_mu. A block
   ! shifted by the lattice vector hT' gives the same values at U = T - T', so all lattice copies
   ! of the class share them. Slots are indexed by (U, A).
   TYPE block_class_cache_type
      INTEGER                           :: block_class = -1, n_slots = 0, n_used = 0
      INTEGER, ALLOCATABLE              :: keys(:, :), slot(:, :, :, :)
      LOGICAL, ALLOCATABLE              :: evaluated(:, :)
      REAL(dp), ALLOCATABLE             :: point_bound(:, :), values(:, :, :)
   END TYPE block_class_cache_type

   TYPE ao_image_type
      TYPE(gto_basis_set_type), POINTER :: basis => NULL()
      ! Atom A and cell R of the image, centred at R_A + hR.
      INTEGER                           :: atom = 0, cell(3) = 0
      REAL(dp)                          :: center(3) = 0.0_dp
      REAL(dp)                          :: bound = 0.0_dp, decay = 0.0_dp
      REAL(dp), ALLOCATABLE             :: primary(:, :)
      INTEGER, ALLOCATABLE              :: primary_l(:)
      REAL(dp), ALLOCATABLE             :: envelope_exponents(:), &
                                           radial_envelope(:, :)
   END TYPE ao_image_type

   ! Per-atom interpolation tables for the erfc kernel.
   TYPE radial_interpolation_type
      REAL(dp), ALLOCATABLE             :: short_abscissae(:)
      REAL(dp), ALLOCATABLE             :: coefficients(:, :), &
                                           series(:, :, :), &
                                           short_series(:, :, :), &
                                           short_moments(:, :, :, :)
      ! Radial value and derivative, indexed (2, primitive/angular channel,
      ! squared-radius index).
      REAL(dp), ALLOCATABLE             :: mesh_moments(:, :, :)
      INTEGER, ALLOCATABLE              :: mesh_first(:, :), &
                                           spherical_mesh_first(:, :)
      ! An optional exact factorization of primitive Cartesian rows through
      ! spherical RI functions.
      REAL(dp), ALLOCATABLE             :: spherical_selection(:, :), &
                                           spherical_angular(:, :, :, :), &
                                           spherical_primitive(:, :, :)
      LOGICAL                           :: spherical_factorized = .FALSE.
      REAL(dp)                          :: mesh_spacing = 0.0_dp, &
                                           mesh_offset = 0.0_dp, &
                                           mesh_exponent = 0.0_dp
      REAL(dp)                          :: attenuation = 0.0_dp
      ! Bound on all |V_P(r)| of the atom on radial intervals, indexed from r = 0.
      REAL(dp), ALLOCATABLE             :: short_bounds(:)
   END TYPE radial_interpolation_type

END MODULE ri_rs_kp_types

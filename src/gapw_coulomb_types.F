!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Sparse Gaussian-product maps and geometry-dependent augmentation Coulomb tensors.
! **************************************************************************************************
MODULE gapw_coulomb_types
   USE kinds,                           ONLY: dp

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gapw_coulomb_term_type, gapw_coulomb_kind_type, gapw_coulomb_pair_type, gapw_coulomb_cache_type

   TYPE gapw_coulomb_term_type
      REAL(dp) :: alpha = 0.0_dp
      INTEGER :: power = 0, l = 0, offset = 0
      ! Upper triangle of the padded spherical primitive product matrix.
      INTEGER, ALLOCATABLE :: row(:), col(:), m(:)
      REAL(dp), ALLOCATABLE :: factor(:)
   END TYPE gapw_coulomb_term_type

   TYPE gapw_coulomb_kind_type
      INTEGER :: nprimitive = 0, ncomponent = 0
      TYPE(gapw_coulomb_term_type), ALLOCATABLE :: term(:)
   END TYPE gapw_coulomb_kind_type

   TYPE gapw_coulomb_pair_type
      ! Last index: potential, three displacement derivatives, nine strain derivatives.
      REAL(dp), ALLOCATABLE :: integral(:, :, :, :)
   END TYPE gapw_coulomb_pair_type

   TYPE gapw_coulomb_cache_type
      REAL(dp) :: eta = 0.0_dp, hmat(3, 3) = 0.0_dp
      INTEGER :: periodicity(3) = 1
      REAL(dp), ALLOCATABLE :: position(:, :)
      TYPE(gapw_coulomb_kind_type), ALLOCATABLE :: kind(:)
      TYPE(gapw_coulomb_pair_type), ALLOCATABLE :: pair(:, :)
   END TYPE gapw_coulomb_cache_type
END MODULE gapw_coulomb_types

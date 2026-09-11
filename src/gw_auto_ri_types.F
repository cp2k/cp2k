!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Input and persistent data for automatic RI basis optimization.
!> \par History
!>      09.2026 created [Jan Wilhelm]
! **************************************************************************************************
MODULE gw_auto_ri_types
   USE base_hooks,                      ONLY: timeset,&
                                              timestop
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_auto_ri_types'

   PUBLIC :: auto_ri_release, auto_ri_type

   TYPE auto_ri_type
      ! Input parameters.
      REAL(KIND=dp)                                      :: delta = -1.0_dp
      REAL(KIND=dp)                                      :: neighbor_radius = -1.0_dp
      LOGICAL                                            :: enabled = .FALSE.
      REAL(KIND=dp)                                      :: ri_ao_ratio = -1.0_dp

      ! Heuristic parameters.
      REAL(KIND=dp)                                      :: metric_null_threshold = 1.0E-10_dp
      REAL(KIND=dp)                                      :: metric_rank_threshold = 1.0E-6_dp
      REAL(KIND=dp)                                      :: min_score_ratio = 1.0E-12_dp
      REAL(KIND=dp)                                      :: pair_absolute_threshold = 1.0E-12_dp
      REAL(KIND=dp)                                      :: psd_tolerance = 1.0E-10_dp

      ! Optimized basis and metric matrices.
      TYPE(cp_fm_type)                                   :: coulomb_matrix = cp_fm_type()
      TYPE(cp_fm_type)                                   :: inverse_ri_metric = cp_fm_type()
      INTEGER                                            :: fit_group_count = 0
      LOGICAL                                            :: ready = .FALSE.
      INTEGER, ALLOCATABLE                               :: optimized_size(:)
      REAL(KIND=dp), ALLOCATABLE                         :: transform(:, :, :)

      ! AA and AB coefficient matrices, packed consecutively in column-major order.
      INTEGER, ALLOCATABLE                               :: fit_group_atom_a(:), fit_group_atom_b(:)
      INTEGER, ALLOCATABLE                               :: fit_group_first_a(:), fit_group_first_b(:)
      INTEGER, ALLOCATABLE                               :: fit_group_n_to_a(:)
      INTEGER, ALLOCATABLE                               :: fit_group_nrow(:), fit_group_ncolumn(:)
      INTEGER, ALLOCATABLE                               :: fit_group_offset(:)
      REAL(KIND=dp), ALLOCATABLE                         :: fit_group_coefficient(:)
   END TYPE auto_ri_type

CONTAINS

! **************************************************************************************************
!> \brief Releases all data owned by an optimized AUTO_RI basis.
!> \param auto_ri AUTO_RI input and optimized-basis data
! **************************************************************************************************
   SUBROUTINE auto_ri_release(auto_ri)
      TYPE(auto_ri_type), INTENT(INOUT)                  :: auto_ri

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'auto_ri_release'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ALLOCATED(auto_ri%optimized_size)) DEALLOCATE (auto_ri%optimized_size)
      IF (ALLOCATED(auto_ri%transform)) DEALLOCATE (auto_ri%transform)
      IF (ALLOCATED(auto_ri%fit_group_atom_a)) DEALLOCATE (auto_ri%fit_group_atom_a)
      IF (ALLOCATED(auto_ri%fit_group_atom_b)) DEALLOCATE (auto_ri%fit_group_atom_b)
      IF (ALLOCATED(auto_ri%fit_group_first_a)) DEALLOCATE (auto_ri%fit_group_first_a)
      IF (ALLOCATED(auto_ri%fit_group_first_b)) DEALLOCATE (auto_ri%fit_group_first_b)
      IF (ALLOCATED(auto_ri%fit_group_n_to_a)) DEALLOCATE (auto_ri%fit_group_n_to_a)
      IF (ALLOCATED(auto_ri%fit_group_nrow)) DEALLOCATE (auto_ri%fit_group_nrow)
      IF (ALLOCATED(auto_ri%fit_group_ncolumn)) DEALLOCATE (auto_ri%fit_group_ncolumn)
      IF (ALLOCATED(auto_ri%fit_group_offset)) DEALLOCATE (auto_ri%fit_group_offset)
      IF (ALLOCATED(auto_ri%fit_group_coefficient)) DEALLOCATE (auto_ri%fit_group_coefficient)
      CALL cp_fm_release(auto_ri%coulomb_matrix)
      CALL cp_fm_release(auto_ri%inverse_ri_metric)
      auto_ri%fit_group_count = 0
      auto_ri%ready = .FALSE.

      CALL timestop(handle)

   END SUBROUTINE auto_ri_release

END MODULE gw_auto_ri_types

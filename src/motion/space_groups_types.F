!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Space Group Symmetry Type Module  (version 1.0, Ferbruary 12, 2021)
!> \par History
!>      Pierre-André Cazade [pcazade] 02.2021 - University of Limerick
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
MODULE space_groups_types

   USE cell_types,                      ONLY: cell_release,&
                                              cell_type
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'space_groups_types'

   TYPE spgr_type
      LOGICAL                                     :: keep_space_group = .FALSE.
      LOGICAL                                     :: symlib = .FALSE.
      LOGICAL                                     :: print_atoms = .FALSE.
      INTEGER                                     :: iunit = -1
      INTEGER                                     :: istriz = -1
      REAL(KIND=dp)                               :: eps_symmetry = 1.0e-4_dp
      INTEGER                                     :: nparticle = 0
      INTEGER                                     :: nparticle_sym = 0
      INTEGER                                     :: n_atom = 0
      INTEGER                                     :: n_core = 0
      INTEGER                                     :: n_shell = 0
      INTEGER                                     :: n_atom_sym = 0
      INTEGER                                     :: n_core_sym = 0
      INTEGER                                     :: n_shell_sym = 0
      INTEGER, DIMENSION(:), ALLOCATABLE          :: atype
      INTEGER, DIMENSION(:, :), ALLOCATABLE       :: eqatom
      LOGICAL, DIMENSION(:), ALLOCATABLE          :: lop, lat
      REAL(KIND=dp), DIMENSION(3)                 :: pol = 0.0_dp
      !SPGLIB
      INTEGER                                     :: space_group_number = 0
      CHARACTER(len=11)                           :: international_symbol = ""
      CHARACTER(len=6)                            :: pointgroup_symbol = ""
      CHARACTER(len=7)                            :: schoenflies = ""
      INTEGER                                     :: n_operations = 0
      INTEGER                                     :: n_reduced_operations = 0
      INTEGER                                     :: n_operations_subset = 0
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE    :: rotations
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE    :: rotations_subset
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: translations
      TYPE(cell_type), POINTER                    :: cell_ref => NULL()
   END TYPE spgr_type

   PUBLIC :: spgr_type, release_spgr_type

CONTAINS

! **************************************************************************************************
!> \brief Release the SPGR type
!> \param spgr  The SPGR type
!> \par History
!>      01.2020 created [pcazade]
!> \author Pierre-André Cazade (first version)
! **************************************************************************************************
   SUBROUTINE release_spgr_type(spgr)

      TYPE(spgr_type), POINTER                           :: spgr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'release_spgr_type', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(spgr)) THEN

         IF (ALLOCATED(spgr%rotations)) THEN
            DEALLOCATE (spgr%rotations)
         END IF
         IF (ALLOCATED(spgr%rotations_subset)) THEN
            DEALLOCATE (spgr%rotations_subset)
         END IF
         IF (ALLOCATED(spgr%translations)) THEN
            DEALLOCATE (spgr%translations)
         END IF
         IF (ALLOCATED(spgr%atype)) THEN
            DEALLOCATE (spgr%atype)
         END IF
         IF (ALLOCATED(spgr%eqatom)) THEN
            DEALLOCATE (spgr%eqatom)
         END IF
         IF (ALLOCATED(spgr%lop)) THEN
            DEALLOCATE (spgr%lop)
         END IF
         IF (ALLOCATED(spgr%lat)) THEN
            DEALLOCATE (spgr%lat)
         END IF

         CALL cell_release(spgr%cell_ref)

         DEALLOCATE (spgr)
      END IF

   END SUBROUTINE release_spgr_type

END MODULE space_groups_types

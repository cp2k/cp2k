!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
MODULE integration_grid_types

   USE kinds,                           ONLY: dp
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'integration_grid_types'

   TYPE grid_batch_val_1d_type
      INTEGER                                                  :: np1
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE                 :: val1d
   END TYPE grid_batch_val_1d_type

   TYPE grid_batch_val_2d_type
      INTEGER                                                  :: np1, np2
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE              :: val2d
   END TYPE grid_batch_val_2d_type

   TYPE gnlist_type
      INTEGER, DIMENSION(:), ALLOCATABLE                       :: atom_list
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE              :: atom_pos
   END TYPE gnlist_type

   TYPE grid_batch_info_type
      INTEGER                                                  :: np
      INTEGER                                                  :: ref_atom
      INTEGER                                                  :: ibatch
      TYPE(gnlist_type)                                        :: gnlist
      REAL(KIND=dp), DIMENSION(3)                              :: rcenter
      REAL(KIND=dp)                                            :: radius
      REAL(dp), DIMENSION(:, :), ALLOCATABLE                   :: rco
      REAL(dp), DIMENSION(:), ALLOCATABLE                      :: weight
      REAL(dp), DIMENSION(:), ALLOCATABLE                      :: wref
      REAL(dp), DIMENSION(:), ALLOCATABLE                      :: wsum
   END TYPE grid_batch_info_type

   TYPE integration_grid_type
      INTEGER                                                  :: nbatch
      TYPE(grid_batch_info_type), DIMENSION(:), ALLOCATABLE    :: grid_batch
   END TYPE integration_grid_type

   TYPE integration_grid_value_type
      INTEGER                                                  :: nbatch
      TYPE(grid_batch_val_1d_type), DIMENSION(:), ALLOCATABLE  :: grid_val_1d
      TYPE(grid_batch_val_2d_type), DIMENSION(:), ALLOCATABLE  :: grid_val_2d
   END TYPE integration_grid_value_type

   PUBLIC :: integration_grid_type, allocate_intgrid, deallocate_intgrid
   PUBLIC :: integration_grid_value_type, allocate_intgrid_val, deallocate_intgrid_val

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief   Initialize integration_grid_type
!> \param int_grid ...
!> \date    02.2018
!> \param
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE allocate_intgrid(int_grid)

      TYPE(integration_grid_type), POINTER               :: int_grid

      IF (ASSOCIATED(int_grid)) CALL deallocate_intgrid(int_grid)
      ALLOCATE (int_grid)
      int_grid%nbatch = 0

   END SUBROUTINE allocate_intgrid

! **************************************************************************************************
!> \brief   Deallocate integration_grid_type
!> \param int_grid ...
!> \date    02.2018
!> \param
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE deallocate_intgrid(int_grid)
      TYPE(integration_grid_type), POINTER               :: int_grid

      INTEGER                                            :: i

      IF (ASSOCIATED(int_grid)) THEN
         IF (ALLOCATED(int_grid%grid_batch)) THEN
            DO i = 1, int_grid%nbatch
               IF (ALLOCATED(int_grid%grid_batch(i)%rco)) DEALLOCATE (int_grid%grid_batch(i)%rco)
               IF (ALLOCATED(int_grid%grid_batch(i)%weight)) DEALLOCATE (int_grid%grid_batch(i)%weight)
               IF (ALLOCATED(int_grid%grid_batch(i)%wref)) DEALLOCATE (int_grid%grid_batch(i)%wref)
               IF (ALLOCATED(int_grid%grid_batch(i)%wsum)) DEALLOCATE (int_grid%grid_batch(i)%wsum)
               !
               IF (ALLOCATED(int_grid%grid_batch(i)%gnlist%atom_list)) DEALLOCATE (int_grid%grid_batch(i)%gnlist%atom_list)
               IF (ALLOCATED(int_grid%grid_batch(i)%gnlist%atom_pos)) DEALLOCATE (int_grid%grid_batch(i)%gnlist%atom_pos)
            END DO
            DEALLOCATE (int_grid%grid_batch)
         END IF
         DEALLOCATE (int_grid)
      ELSE
         CALL cp_abort(__LOCATION__, &
                       "The pointer int_grid is not associated and "// &
                       "cannot be deallocated")
      END IF
   END SUBROUTINE deallocate_intgrid

! **************************************************************************************************
!> \brief   Initialize integration_grid_value_type
!> \param int_grid ...
!> \date    02.2018
!> \param
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE allocate_intgrid_val(int_grid)

      TYPE(integration_grid_value_type), POINTER         :: int_grid

      IF (ASSOCIATED(int_grid)) CALL deallocate_intgrid_val(int_grid)
      ALLOCATE (int_grid)
      int_grid%nbatch = 0

   END SUBROUTINE allocate_intgrid_val

! **************************************************************************************************
!> \brief   Deallocate integration_grid_value_type
!> \param int_grid ...
!> \date    02.2018
!> \param
!> \author  JGH
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE deallocate_intgrid_val(int_grid)
      TYPE(integration_grid_value_type), POINTER         :: int_grid

      INTEGER                                            :: i

      IF (ASSOCIATED(int_grid)) THEN
         IF (ALLOCATED(int_grid%grid_val_1d)) THEN
            DO i = 1, int_grid%nbatch
               IF (ALLOCATED(int_grid%grid_val_1d(i)%val1d)) DEALLOCATE (int_grid%grid_val_1d(i)%val1d)
            END DO
            DEALLOCATE (int_grid%grid_val_1d)
         END IF
         IF (ALLOCATED(int_grid%grid_val_2d)) THEN
            DO i = 1, int_grid%nbatch
               IF (ALLOCATED(int_grid%grid_val_2d(i)%val2d)) DEALLOCATE (int_grid%grid_val_2d(i)%val2d)
            END DO
            DEALLOCATE (int_grid%grid_val_2d)
         END IF
         DEALLOCATE (int_grid)
      ELSE
         CALL cp_abort(__LOCATION__, &
                       "The pointer int_grid is not associated and "// &
                       "cannot be deallocated")
      END IF
   END SUBROUTINE deallocate_intgrid_val

END MODULE integration_grid_types

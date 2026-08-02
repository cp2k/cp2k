!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM cell_types_TEST
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE kinds,                           ONLY: dp

   IMPLICIT NONE

   TYPE(cell_type), POINTER                           :: cell

   ALLOCATE (cell)
   cell%hmat = 0.0_dp
   cell%h_inv = 0.0_dp
   cell%hmat(1, 1) = 2.0_dp
   cell%hmat(2, 2) = 3.0_dp
   cell%hmat(3, 3) = 4.0_dp
   cell%h_inv(1, 1) = 0.5_dp
   cell%h_inv(2, 2) = 1.0_dp/3.0_dp
   cell%h_inv(3, 3) = 0.25_dp
   cell%perd = 1

   cell%orthorhombic = .TRUE.
   CALL check_half_cell_boundaries(cell)

   cell%orthorhombic = .FALSE.
   CALL check_half_cell_boundaries(cell)

   DEALLOCATE (cell)

CONTAINS

! **************************************************************************************************
!> \brief Check that roundoff around half-cell boundaries does not change the periodic image.
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE check_half_cell_boundaries(cell)

      TYPE(cell_type), POINTER                           :: cell

      INTEGER, DIMENSION(3), PARAMETER                   :: image_offset = [1, -1, 2]
      REAL(KIND=dp), PARAMETER                           :: delta = 8.0_dp*EPSILON(1.0_dp)

      REAL(KIND=dp), DIMENSION(3)                        :: r, scaled, wrapped, wrapped_scaled

      scaled = [0.5_dp - delta, -0.5_dp + delta, 1.5_dp - delta]
      r = scaled*[cell%hmat(1, 1), cell%hmat(2, 2), cell%hmat(3, 3)]
      wrapped = pbc(r, cell)
      wrapped_scaled = wrapped/[cell%hmat(1, 1), cell%hmat(2, 2), cell%hmat(3, 3)]
      IF (ANY(wrapped_scaled >= 0.0_dp)) THEN
         ERROR STOP "Roundoff changed the image at a half-cell boundary"
      END IF
      IF (ANY(ABS(ABS(wrapped_scaled) - 0.5_dp) > 2.0_dp*delta)) THEN
         ERROR STOP "Unexpected wrapped coordinate at a half-cell boundary"
      END IF

      scaled = [0.5_dp - 128.0_dp*EPSILON(1.0_dp), &
                -0.5_dp - 128.0_dp*EPSILON(1.0_dp), &
                1.5_dp + 128.0_dp*EPSILON(1.0_dp)]
      r = scaled*[cell%hmat(1, 1), cell%hmat(2, 2), cell%hmat(3, 3)]
      wrapped = pbc(r, cell)
      wrapped_scaled = wrapped/[cell%hmat(1, 1), cell%hmat(2, 2), cell%hmat(3, 3)]
      IF (ANY(SIGN(1.0_dp, wrapped_scaled) /= [1.0_dp, 1.0_dp, -1.0_dp])) THEN
         ERROR STOP "Coordinate outside the boundary tolerance changed image"
      END IF

      scaled = [0.5_dp - delta, -0.5_dp + delta, 1.5_dp - delta]
      r = scaled*[cell%hmat(1, 1), cell%hmat(2, 2), cell%hmat(3, 3)]
      wrapped = pbc(r, cell, image_offset)
      wrapped_scaled = wrapped/[cell%hmat(1, 1), cell%hmat(2, 2), cell%hmat(3, 3)]
      IF (ANY(ABS(wrapped_scaled - ([-0.5_dp, -0.5_dp, -0.5_dp] + &
                                    REAL(image_offset, KIND=dp))) > 2.0_dp*delta)) THEN
         ERROR STOP "Periodic image offset is inconsistent at a half-cell boundary"
      END IF

   END SUBROUTINE check_half_cell_boundaries

END PROGRAM cell_types_TEST

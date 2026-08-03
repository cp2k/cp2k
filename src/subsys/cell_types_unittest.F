!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM cell_types_TEST
   USE cell_types,                      ONLY: cell_type,&
                                              pbc_stable,&
                                              real_to_scaled,&
                                              scaled_to_real
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
   cell%hmat(1, 2) = 0.5_dp
   cell%h_inv(1, 2) = -1.0_dp/12.0_dp
   CALL check_half_cell_boundaries(cell)

   DEALLOCATE (cell)

CONTAINS

! **************************************************************************************************
!> \brief Check that roundoff around half-cell boundaries does not change the periodic image.
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE check_half_cell_boundaries(cell)

      TYPE(cell_type), POINTER                           :: cell

      INTEGER, DIMENSION(3), PARAMETER :: lattice_translation = [-1, 1, -2]
      REAL(KIND=dp), PARAMETER                           :: delta = 8.0_dp*EPSILON(1.0_dp)

      REAL(KIND=dp), DIMENSION(3)                        :: r, scaled, wrapped, wrapped_scaled, &
                                                            wrapped_translated, &
                                                            wrapped_translated_scaled

      scaled = [0.5_dp - delta, -0.5_dp + delta, 1.5_dp - delta]
      CALL scaled_to_real(r, scaled, cell)
      wrapped = pbc_stable(r, cell)
      CALL real_to_scaled(wrapped_scaled, wrapped, cell)
      IF (ANY(SIGN(1.0_dp, wrapped_scaled) /= [-1.0_dp, -1.0_dp, -1.0_dp])) THEN
         ERROR STOP "Roundoff changed the image at a half-cell boundary"
      END IF
      IF (ANY(ABS(ABS(wrapped_scaled) - 0.5_dp) > 2.0_dp*delta)) THEN
         ERROR STOP "Unexpected wrapped coordinate at a half-cell boundary"
      END IF

      CALL scaled_to_real(r, scaled + REAL(lattice_translation, KIND=dp), cell)
      wrapped_translated = pbc_stable(r, cell)
      CALL real_to_scaled(wrapped_translated_scaled, wrapped_translated, cell)
      IF (ANY(ABS(wrapped_translated_scaled - wrapped_scaled) > 4.0_dp*delta)) THEN
         ERROR STOP "Periodic wrapping changed under a lattice translation"
      END IF

      scaled = [0.5_dp - 128.0_dp*EPSILON(1.0_dp), &
                -0.5_dp - 128.0_dp*EPSILON(1.0_dp), &
                1.5_dp + 128.0_dp*EPSILON(1.0_dp)]
      CALL scaled_to_real(r, scaled, cell)
      wrapped = pbc_stable(r, cell)
      CALL real_to_scaled(wrapped_scaled, wrapped, cell)
      IF (ANY(SIGN(1.0_dp, wrapped_scaled) /= [1.0_dp, 1.0_dp, -1.0_dp])) THEN
         ERROR STOP "Coordinate outside the boundary tolerance changed image"
      END IF

   END SUBROUTINE check_half_cell_boundaries

END PROGRAM cell_types_TEST

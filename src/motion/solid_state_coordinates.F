!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Rotation-free atom and cell coordinates for solid-state transition searches
! **************************************************************************************************
MODULE solid_state_coordinates
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: det_3x3,&
                                              inv_3x3
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER, PUBLIC :: ss_cell_dof = 6
   INTEGER, PARAMETER :: row(ss_cell_dof) = [1, 2, 3, 1, 1, 2], col(ss_cell_dof) = [1, 2, 3, 2, 3, 3]

   TYPE, PUBLIC :: solid_state_frame_type
      REAL(KIND=dp) :: rotation(3, 3) = 0.0_dp, reference(3, 3) = 0.0_dp, width = 0.0_dp
   END TYPE solid_state_frame_type

   PUBLIC :: solid_state_frame, solid_state_pack, solid_state_unpack, solid_state_force

CONTAINS

! **************************************************************************************************
!> \brief Fix the orientation and length scale from the reference cell (columns are cell vectors).
!> \param frame reference coordinate frame
!> \param h reference cell matrix
!> \param natom number of atoms
!> \param weight dimensionless scale for cell motion
! **************************************************************************************************
   SUBROUTINE solid_state_frame(frame, h, natom, weight)
      TYPE(solid_state_frame_type), INTENT(OUT)          :: frame
      REAL(KIND=dp), INTENT(IN)                          :: h(3, 3)
      INTEGER, INTENT(IN)                                :: natom
      REAL(KIND=dp), INTENT(IN)                          :: weight

      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: volume

      volume = det_3x3(h)
      IF (volume <= 0.0_dp) CALL cp_abort(__LOCATION__, &
                                          "Solid-state transition searches require right-handed cells.")
      CPASSERT(natom > 0 .AND. weight > 0.0_dp)
      frame%reference = 0.0_dp
      DO j = 1, 3
         frame%rotation(:, j) = h(:, j)
         DO i = 1, j - 1
            frame%reference(i, j) = DOT_PRODUCT(frame%rotation(:, i), frame%rotation(:, j))
            frame%rotation(:, j) = frame%rotation(:, j) - frame%reference(i, j)*frame%rotation(:, i)
         END DO
         frame%reference(j, j) = NORM2(frame%rotation(:, j))
         IF (frame%reference(j, j) <= 1.0E-10_dp) THEN
            CPABORT("Degenerate solid-state reference cell.")
         END IF
         frame%rotation(:, j) = frame%rotation(:, j)/frame%reference(j, j)
      END DO
      frame%width = weight*SQRT(REAL(natom, dp))*(volume/natom)**(1.0_dp/3.0_dp)
   END SUBROUTINE solid_state_frame

! **************************************************************************************************
!> \brief Convert a cell and Cartesian positions to atom and cell coordinates in length units.
!> \param frame reference coordinate frame
!> \param h cell matrix
!> \param r Cartesian positions, flattened by atom
!> \param x joint coordinates: reference Cartesian positions and scaled cell deformation
! **************************************************************************************************
   SUBROUTINE solid_state_pack(frame, h, r, x)
      TYPE(solid_state_frame_type), INTENT(IN)           :: frame
      REAL(KIND=dp), INTENT(IN)                          :: h(3, 3), r(:)
      REAL(KIND=dp), INTENT(OUT)                         :: x(:)

      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: deformation(3, 3), inverse(3, 3)

      n = SIZE(r)
      CPASSERT(MOD(n, 3) == 0 .AND. SIZE(x) == n + ss_cell_dof)
      deformation = MATMUL(TRANSPOSE(frame%rotation), MATMUL(h, inv_3x3(frame%reference)))
      IF (MAX(ABS(deformation(2, 1)), ABS(deformation(3, 1)), ABS(deformation(3, 2))) > &
          1.0E-8_dp*MAX(1.0_dp, MAXVAL(ABS(deformation)))) THEN
         CALL cp_abort(__LOCATION__, &
                       "Solid-state replicas must have a common cell orientation and lattice correspondence.")
      END IF
      DO i = 1, 3
         IF (deformation(i, i) <= 0.0_dp) THEN
            CPABORT("Solid-state cell vectors must point consistently.")
         END IF
      END DO
      inverse = MATMUL(frame%reference, inv_3x3(h))
      DO i = 1, n, 3
         x(i:i + 2) = MATMUL(inverse, r(i:i + 2))
      END DO
      DO i = 1, ss_cell_dof
         x(n + i) = frame%width*deformation(row(i), col(i))
      END DO
      DO i = 1, 3
         x(n + i) = frame%width*LOG(deformation(i, i))
      END DO
   END SUBROUTINE solid_state_pack

! **************************************************************************************************
!> \brief Recover the cell and Cartesian positions without changing the reference orientation.
!> \param frame reference coordinate frame
!> \param x joint coordinates in length units
!> \param h cell matrix
!> \param r Cartesian positions, flattened by atom
! **************************************************************************************************
   SUBROUTINE solid_state_unpack(frame, x, h, r)
      TYPE(solid_state_frame_type), INTENT(IN)           :: frame
      REAL(KIND=dp), INTENT(IN)                          :: x(:)
      REAL(KIND=dp), INTENT(OUT)                         :: h(3, 3), r(:)

      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: transform(3, 3)

      n = SIZE(r)
      CPASSERT(MOD(n, 3) == 0 .AND. SIZE(x) == n + ss_cell_dof)
      transform = MATMUL(frame%rotation, cell_deformation(frame, x(n + 1:)))
      h = MATMUL(transform, frame%reference)
      DO i = 1, n, 3
         r(i:i + 2) = MATMUL(transform, x(i:i + 2))
      END DO
   END SUBROUTINE solid_state_unpack

! **************************************************************************************************
!> \brief Pull back minus the derivative of E+pV to the joint coordinates.
!> \param frame reference coordinate frame
!> \param x joint coordinates in length units
!> \param f Cartesian forces, flattened by atom
!> \param virial configurational virial, including the affine motion at fixed fractional coordinates
!> \param pv external hydrostatic pressure times current volume
!> \param force force in the joint coordinates
! **************************************************************************************************
   SUBROUTINE solid_state_force(frame, x, f, virial, pv, force)
      TYPE(solid_state_frame_type), INTENT(IN)           :: frame
      REAL(KIND=dp), INTENT(IN)                          :: x(:), f(:), virial(3, 3), pv
      REAL(KIND=dp), INTENT(OUT)                         :: force(:)

      INTEGER                                            :: i, n
      REAL(KIND=dp)                                      :: deformation(3, 3), transform(3, 3), &
                                                            w(3, 3)

      n = SIZE(f)
      CPASSERT(MOD(n, 3) == 0 .AND. SIZE(x) == n + ss_cell_dof)
      CPASSERT(SIZE(force) == SIZE(x))
      deformation = cell_deformation(frame, x(n + 1:))
      transform = MATMUL(frame%rotation, deformation)
      DO i = 1, n, 3
         force(i:i + 2) = MATMUL(TRANSPOSE(transform), f(i:i + 2))
      END DO
      ! The atomic coordinates are proportional to fractional coordinates. The
      ! virial already includes their affine motion. Do not add f*r a second time.
      w = MATMUL(TRANSPOSE(frame%rotation), MATMUL(virial, frame%rotation))
      DO i = 1, 3
         w(i, i) = w(i, i) - pv
      END DO
      w = MATMUL(w, TRANSPOSE(inv_3x3(deformation)))/frame%width
      DO i = 1, ss_cell_dof
         force(n + i) = w(row(i), col(i))
      END DO
      DO i = 1, 3
         force(n + i) = force(n + i)*deformation(i, i)
      END DO
   END SUBROUTINE solid_state_force

! **************************************************************************************************
!> \brief Construct the upper-triangular deformation from its six coordinates.
!> \param frame reference coordinate frame
!> \param q scaled cell deformation coordinates
!> \return the deformation matrix
! **************************************************************************************************
   FUNCTION cell_deformation(frame, q) RESULT(deformation)
      TYPE(solid_state_frame_type), INTENT(IN)           :: frame
      REAL(KIND=dp), INTENT(IN)                          :: q(ss_cell_dof)
      REAL(KIND=dp)                                      :: deformation(3, 3)

      INTEGER                                            :: i

      IF (ANY(ABS(q) > 30.0_dp*frame%width)) THEN
         CALL cp_abort(__LOCATION__, &
                       "Solid-state cell step too large. Reduce the optimization step size.")
      END IF
      deformation = 0.0_dp
      DO i = 1, ss_cell_dof
         deformation(row(i), col(i)) = q(i)/frame%width
      END DO
      DO i = 1, 3
         deformation(i, i) = EXP(deformation(i, i))
      END DO
   END FUNCTION cell_deformation
END MODULE solid_state_coordinates

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM rirs_grid_utils_unittest
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE rirs_grid_utils,                 ONLY: initialize_rirs_grid

   IMPLICIT NONE

   INTEGER                                             :: n_candidates
   REAL(KIND=dp), ALLOCATABLE                          :: points(:, :)
   TYPE(particle_type), DIMENSION(2), TARGET            :: particles
   TYPE(particle_type), DIMENSION(:), POINTER           :: particle_set

   particle_set => particles
   particle_set(1)%r = [0.0_dp, 0.0_dp, 0.0_dp]
   particle_set(2)%r = [2.0_dp, 0.0_dp, 0.0_dp]
   ALLOCATE (points(3, 4))
   points(:, :) = RESHAPE([ &
                          0.0_dp, 0.0_dp, 0.0_dp, &
                          0.5_dp, 0.0_dp, 0.0_dp, &
                          1.0_dp, 0.0_dp, 0.0_dp, &
                          1.5_dp, 0.0_dp, 0.0_dp], [3, 4])

   CALL initialize_rirs_grid(points, 2, 1, particle_set, n_candidates)
   IF (n_candidates /= 3) ERROR STOP "incorrect RI-RS Voronoi candidate count"
   IF (MAXVAL(ABS(points(:, 1) - [0.0_dp, 0.0_dp, 0.0_dp])) > 1.0E-14_dp) THEN
      ERROR STOP "RI-RS maximin selection did not start at the origin"
   END IF
   IF (MAXVAL(ABS(points(:, 2) - [1.0_dp, 0.0_dp, 0.0_dp])) > 1.0E-14_dp) THEN
      ERROR STOP "RI-RS maximin selection returned the wrong second point"
   END IF

   DEALLOCATE (points)
   ALLOCATE (points(3, 5))
   points(:, :) = RESHAPE([ &
                          0.0_dp, 0.0_dp, 0.0_dp, &
                          0.5_dp, 0.0_dp, 0.0_dp, &
                          1.0_dp, 0.0_dp, 0.0_dp, &
                          1.5_dp, 0.0_dp, 0.0_dp, &
                          2.0_dp, 0.0_dp, 0.0_dp], [3, 5])
   CALL initialize_rirs_grid(points, 4, 1, particle_set, n_candidates)
   IF (n_candidates /= 3) ERROR STOP "incorrect RI-RS Voronoi fallback candidate count"
   IF (SIZE(points, 2) /= 4) ERROR STOP "RI-RS Voronoi fallback returned the wrong grid size"
   IF (MAXVAL(ABS(points(:, 2) - [2.0_dp, 0.0_dp, 0.0_dp])) > 1.0E-14_dp) THEN
      ERROR STOP "RI-RS Voronoi fallback did not select from the full source grid"
   END IF

END PROGRAM rirs_grid_utils_unittest

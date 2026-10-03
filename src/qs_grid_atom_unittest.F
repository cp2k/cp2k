!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM qs_grid_atom_unittest
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: do_gapw_gcs,&
                                              do_gapw_gct,&
                                              do_gapw_log
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: pi
   USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                              create_grid_atom,&
                                              deallocate_grid_atom,&
                                              grid_atom_type

   IMPLICIT NONE

   INTEGER :: n, quadrature, size_index, type_index
   INTEGER, PARAMETER :: sizes(3) = [400, 6400, 12800]
   INTEGER, PARAMETER :: types(3) = [do_gapw_gcs, do_gapw_gct, do_gapw_log]
   REAL(dp) :: error, maximum_error
   TYPE(grid_atom_type), POINTER :: grid => NULL()

   CALL cp2k_init()
   maximum_error = 0.0_dp
   DO type_index = 1, SIZE(types)
      quadrature = types(type_index)
      DO size_index = 1, SIZE(sizes)
         n = sizes(size_index)
         CALL allocate_grid_atom(grid)
         CALL create_grid_atom(grid, n, 1, 2, 0, quadrature)
         IF (.NOT. ALL(ieee_is_finite(grid%rad))) ERROR STOP "Nonfinite radial points"
         IF (.NOT. ALL(ieee_is_finite(grid%wr))) ERROR STOP "Nonfinite radial weights"
         IF (ANY(grid%rad <= 0.0_dp)) ERROR STOP "Nonpositive radial points"
         IF (ANY(grid%wr <= 0.0_dp)) ERROR STOP "Nonpositive radial weights"
         IF (ANY(grid%rad(:n - 1) <= grid%rad(2:))) ERROR STOP "Unordered radial points"

         ! The weights include r**2; these three integrals are analytic on [0,infinity).
         error = ABS(SUM(grid%wr*EXP(-grid%rad)) - 2.0_dp)
         error = MAX(error, ABS(SUM(grid%wr*EXP(-grid%rad)/grid%rad) - 1.0_dp))
         error = MAX(error, ABS(SUM(grid%wr*EXP(-grid%rad2)) - SQRT(pi)/4.0_dp))
         WRITE (*, '(A,3I8,ES18.8)') "Radial quadrature: ", type_index, n, quadrature, error
         IF (error > 1.e-10_dp) ERROR STOP "Radial moment error"
         maximum_error = MAX(maximum_error, error)

         ! Independent 80-digit endpoint references: no clipping of rounded x to +/-1.
         IF (n == 6400 .AND. quadrature == do_gapw_gct) THEN
            IF (ABS(grid%rad(n)/4.83458313724687745e-18_dp - 1.0_dp) > 1.e-13_dp) THEN
               ERROR STOP "Transformed grid inner endpoint"
            END IF
            IF (ABS(grid%rad(1)/2.06843066219244770e17_dp - 1.0_dp) > 1.e-13_dp) THEN
               ERROR STOP "Transformed grid outer endpoint"
            END IF
         ELSE IF (n == 6400 .AND. quadrature == do_gapw_log) THEN
            IF (ABS(grid%rad(n)/6.97482911687147684e-18_dp - 1.0_dp) > 1.e-13_dp) THEN
               ERROR STOP "Logarithmic grid inner endpoint"
            END IF
            IF (ABS(grid%rad(1)/57.5213142095328827_dp - 1.0_dp) > 1.e-13_dp) THEN
               ERROR STOP "Logarithmic grid outer endpoint"
            END IF
         END IF
         CALL deallocate_grid_atom(grid)
      END DO
   END DO
   WRITE (*, '(A,ES18.8)') "Maximum radial moment error: ", maximum_error
   CALL cp2k_finalize()
END PROGRAM qs_grid_atom_unittest

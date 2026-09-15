! SPDX-License-Identifier: GPL-2.0-or-later
PROGRAM native_grid_cache_fixture
   USE cell_types, ONLY: cell_type
   USE kinds, ONLY: dp, int_8
   USE omp_lib, ONLY: omp_get_wtime, omp_set_num_threads
   USE pw_grid_types, ONLY: pw_grid_type
   USE qs_native_grid_cache
   USE stencil_kernel, ONLY: create_native_grid_interpolation_stencil
   IMPLICIT NONE
   TYPE(native_grid_cache_type) :: cache, other
   TYPE(native_grid_interpolation_stencil_type) :: expected, actual
   TYPE(pw_grid_type), POINTER :: grid
   TYPE(cell_type), POINTER :: cell
   INTEGER, PARAMETER :: nrows = 257
   REAL(dp) :: points(3, nrows), t0, elapsed, checksum
   INTEGER :: i, d, m, w, scenario, pass, rep, misses
   LOGICAL :: hit, wrap

   ALLOCATE (grid, cell)
   CALL omp_set_num_threads(4)
   DO m = 0, 7
      DO w = 0, 1
         wrap = w == 1
         grid%dh_inv = 0.0_dp
         cell%hmat = 0.0_dp
         DO d = 1, 3
            grid%dh_inv(d, d) = 2.0_dp
            cell%hmat(d, d) = 8.0_dp
            cell%perd(d) = MERGE(1, 0, BTEST(m, d - 1))
         END DO
         grid%npts = 16
         DO i = 1, nrows
            points(:, i) = [SIN(REAL(i, dp)), COS(0.7_dp*i), SIN(1.3_dp*i)]*10.0_dp
         END DO
         DO scenario = 0, 8
            SELECT CASE (scenario)
            CASE (1)
               points(1, :) = points(1, :) + 1.0E-6_dp
            CASE (2)
               points = points(:, nrows:1:-1)
            CASE (3)
               grid%dh_inv(1, 2) = 0.21_dp
            CASE (4)
               grid%npts(2) = 17
            CASE (5)
               cell%hmat(1, 2) = 0.03_dp
            CASE (6)
               cell%perd(1) = 1 - cell%perd(1)
            CASE (7)
               wrap = .NOT. wrap
            CASE (8)
               CALL release_native_grid_cache(cache)
            END SELECT
            CALL prepare_native_grid_cache(cache, grid, cell, nrows, wrap)
            DO pass = 1, 2
               misses = 0
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(cache, grid, cell, points, wrap) &
!$OMP PRIVATE(i, expected, actual, hit) REDUCTION(+:misses)
               DO i = 1, nrows
                  CALL create_native_grid_interpolation_stencil(expected, grid, cell, points(:, i), wrap)
                  hit = fetch_native_grid_stencil(cache, i, points(:, i), actual)
                  IF (.NOT. hit) THEN
                     misses = misses + 1
                     actual = expected
                     CALL store_native_grid_stencil(cache, i, points(:, i), actual)
                  END IF
                  CALL compare(expected, actual)
               END DO
!$OMP END PARALLEL DO
               IF (pass == 2 .AND. misses /= 0) ERROR STOP 'Unchanged geometry missed'
               IF (pass == 1 .AND. scenario > 0) THEN
                  ! Reversing an odd-sized row set leaves the middle point unchanged.
                  IF (misses /= nrows - MERGE(1, 0, scenario == 2)) ERROR STOP 'Changed geometry reused'
               END IF
            END DO
            ! Repeated readers of the same rows must never write the cache.
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(cache, grid, cell, points, wrap) PRIVATE(i, expected, actual, hit)
            DO i = 1, 3*nrows
               CALL create_native_grid_interpolation_stencil(expected, grid, cell, points(:, MOD(i, nrows) + 1), wrap)
               hit = fetch_native_grid_stencil(cache, MOD(i, nrows) + 1, points(:, MOD(i, nrows) + 1), actual)
               IF (.NOT. hit) ERROR STOP 'Concurrent read missed'
               CALL compare(expected, actual)
            END DO
!$OMP END PARALLEL DO
         END DO
         CALL release_native_grid_cache(cache)
      END DO
   END DO
   CALL prepare_native_grid_cache(cache, grid, cell, nrows, wrap, max_bytes=1024_int_8)
   CALL create_native_grid_interpolation_stencil(expected, grid, cell, points(:, 1), wrap)
   CALL store_native_grid_stencil(cache, 1, points(:, 1), expected)
   IF (.NOT. fetch_native_grid_stencil(cache, 1, points(:, 1), actual)) ERROR STOP 'Small cache did not store'
   IF (fetch_native_grid_stencil(cache, nrows, points(:, 1), actual)) ERROR STOP 'Memory cap ignored'
   IF (fetch_native_grid_stencil(cache, 0, points(:, 1), actual)) ERROR STOP 'Invalid row hit'
   IF (fetch_native_grid_stencil(other, 1, points(:, 1), actual)) ERROR STOP 'Environment alias'
   CALL prepare_native_grid_cache(cache, grid, cell, 0, wrap)
   IF (fetch_native_grid_stencil(cache, 1, points(:, 1), actual)) ERROR STOP 'Empty rank hit'
   CALL prepare_native_grid_cache(cache, grid, cell, nrows, wrap, max_bytes=0_int_8)
   CALL store_native_grid_stencil(cache, 1, points(:, 1), expected)
   IF (fetch_native_grid_stencil(cache, 1, points(:, 1), actual)) ERROR STOP 'Disabled cache hit'
   CALL release_native_grid_cache(cache)
   CALL release_native_grid_cache(cache)
   WRITE (*, '(A)') 'PASS: exact stencils, motion, reordered grids, cell/grid/periodicity/wrap changes, OpenMP, memory cap, empty rank'
   CALL prepare_native_grid_cache(cache, grid, cell, nrows, wrap)
   DO i = 1, nrows
      CALL create_native_grid_interpolation_stencil(expected, grid, cell, points(:, i), wrap)
      CALL store_native_grid_stencil(cache, i, points(:, i), expected)
   END DO
   DO pass = 1, 2
      checksum = 0.0_dp
      t0 = omp_get_wtime()
      DO rep = 1, 1000
         DO i = 1, nrows
            IF (pass == 1) THEN
               CALL create_native_grid_interpolation_stencil(actual, grid, cell, points(:, i), wrap)
            ELSE
               hit = fetch_native_grid_stencil(cache, i, points(:, i), actual)
               IF (.NOT. hit) ERROR STOP 'Benchmark cache miss'
            END IF
            checksum = checksum + SUM(actual%weight)
         END DO
      END DO
      elapsed = omp_get_wtime() - t0
      WRITE (*, '(A,I0,A,F12.6,A,ES20.12)') 'Stencil benchmark variant=', pass, ' seconds=', elapsed, ' checksum=', checksum
   END DO
CONTAINS
   SUBROUTINE compare(a, b)
      TYPE(native_grid_interpolation_stencil_type), INTENT(IN) :: a, b
      IF (a%active .NEQV. b%active) ERROR STOP 'Active differs'
      IF (ANY(a%valid .NEQV. b%valid)) ERROR STOP 'Valid differs'
      IF (ANY(a%weight /= b%weight)) ERROR STOP 'Weights differ'
      IF (ANY(a%relative_index /= b%relative_index)) ERROR STOP 'Indices differ'
   END SUBROUTINE compare
END PROGRAM native_grid_cache_fixture

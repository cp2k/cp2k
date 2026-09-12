! SPDX-License-Identifier: GPL-2.0-or-later
! Synthetic fields on nonuniform radial grids and real CP2K Lebedev grids.
PROGRAM gapw_atom_adjoint_fixture
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: init_lebedev_grids,&
                                              lebedev_grid
   USE omp_lib,                         ONLY: omp_get_wtime,&
                                              omp_set_num_threads
   USE orbital_pointers,                ONLY: indso,&
                                              init_orbital_pointers
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
   USE original_kernel,                 ONLY: old_fields => interpolate_gapw_atom_grid_fields,&
                                              old_project => project,&
                                              old_weights => atom_grid_interpolation_weights
   USE parallel_kernel,                 ONLY: parallel_project => project
   USE particle_types,                  ONLY: particle_type
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE screened_kernel,                 ONLY: screened_project => project
   USE spherical_harmonics,             ONLY: y_lm
   USE values_only_kernel,              ONLY: new_fields => interpolate_gapw_atom_grid_fields,&
                                              new_weights => atom_grid_interpolation_weights,&
                                              value_project => project

   IMPLICIT NONE
   TYPE(grid_atom_type), POINTER :: grid
   TYPE(harmonics_atom_type), POINTER :: harmonics
   TYPE(cell_type), POINTER :: cell
   TYPE(particle_type), POINTER :: particles(:)
   REAL(dp), POINTER :: vh(:, :, :), vs(:, :, :), th(:, :, :), ts(:, :, :)
   REAL(dp), POINTER :: gh(:, :, :, :), gs(:, :, :, :)
   REAL(dp), ALLOCATABLE :: coords(:, :), density(:, :), gradient(:, :, :), kin(:, :)
   INTEGER, ALLOCATABLE :: atoms(:), local_atoms(:), first(:), last(:)
   REAL(dp), ALLOCATABLE :: reference(:), actual(:)
   LOGICAL :: flags(3)
   INTEGER :: ns, i, n, order, periodic, mask, threads, variant, rep, np
   INTEGER, PARAMETER :: counts(4) = [2, 3, 7, 150], thread_counts(5) = [1, 2, 4, 8, 16]
   REAL(dp) :: t0, elapsed, best, error, max_error
   CHARACTER(32) :: mode, arg
   CHARACTER(12), PARAMETER :: names(4) = [CHARACTER(12) :: 'original', 'values_only', 'screened', 'parallel']

   CALL get_command_argument(1, mode)
   CALL get_command_argument(2, arg)
   READ (arg, *) np
   CALL init_orbital_pointers(6)
   i = -1
   CALL init_spherical_harmonics(6, i)
   CALL init_lebedev_grids()
   ALLOCATE (grid, harmonics, cell, particles(3))
   max_error = 0.0_dp
   DO n = 1, SIZE(counts)
      CALL make_grid(counts(n), 5, 4)
      DO order = 1, 2
         CALL check_weights()
         grid%rad = grid%rad(grid%nr:1:-1)
      END DO
      CALL free_grid()
   END DO
   WRITE (*, '(A)') 'PASS: value weights and full/optional derivatives agree exactly with baseline'

   CALL make_grid(7, 5, 4)
   DO ns = 1, 2
      CALL allocate_potentials()
      CALL check_forward_fields(.FALSE.)
      CALL check_transpose()
      DO periodic = 0, 2
         CALL make_cell(periodic)
         CALL make_points(37)
         DO mask = 0, 7
            flags = [BTEST(mask, 0), BTEST(mask, 1), BTEST(mask, 2)]
            CALL reset_potentials()
            CALL project_variant(1)
            reference = pack_potentials()
            DO threads = 1, SIZE(thread_counts)
               CALL omp_set_num_threads(thread_counts(threads))
               DO variant = 2, 4
                  CALL reset_potentials()
                  CALL project_variant(variant)
                  actual = pack_potentials()
                  error = MAXVAL(ABS(actual - reference))/MAX(1.0_dp, MAXVAL(ABS(reference)))
                  max_error = MAX(max_error, error)
                  IF (error > 2.0E-13_dp) ERROR STOP 'Parallel adjoint differs from original'
               END DO
            END DO
         END DO
         CALL check_target_partitions()
         CALL free_points()
      END DO
      CALL make_points(0)
      flags = .TRUE.
      CALL reset_potentials()
      reference = pack_potentials()
      CALL project_variant(4)
      actual = pack_potentials()
      IF (ANY(actual /= reference)) ERROR STOP 'Empty-rank potentials changed'
      CALL free_points()
      CALL free_potentials()
   END DO
   CALL free_grid()
   WRITE (*, '(A,ES12.4)') 'PASS: transpose, spin/masks, periodic/triclinic images, empty rank; max relative error=', max_error
   IF (TRIM(mode) == 'check') STOP

   CALL make_grid(150, 13, 6)
   CALL make_cell(2)
   CALL make_points(np)
   ns = 2
   flags = .TRUE.
   CALL allocate_potentials()
   CALL check_forward_fields(.TRUE.)
   WRITE (*, '(A,I0,A)') 'BENCHMARK: 150 radial x 770 Lebedev, lmax=6, two spins, ', np, ' rows; best of three seconds'
   WRITE (*, '(A,F8.3,A)') 'Private reduction storage per thread: ', REAL(20*150*770*8, dp)/1024**2, ' MiB'
   DO variant = 1, 4
      DO threads = 1, SIZE(thread_counts)
         IF (variant /= 4 .AND. threads /= 1) CYCLE
         CALL omp_set_num_threads(thread_counts(threads))
         best = HUGE(1.0_dp)
         DO rep = 1, 3
            CALL reset_potentials()
            t0 = omp_get_wtime()
            CALL project_variant(variant)
            elapsed = omp_get_wtime() - t0
            best = MIN(best, elapsed)
         END DO
         actual = pack_potentials()
         IF (variant == 1) reference = actual
         error = MAXVAL(ABS(actual - reference))/MAX(1.0_dp, MAXVAL(ABS(reference)))
         IF (error > 2.0E-12_dp) ERROR STOP 'Benchmark equivalence failed'
         WRITE (*, '(A,1X,I2,1X,F12.6,1X,ES12.4)') names(variant), thread_counts(threads), best, error
      END DO
   END DO
CONTAINS
   SUBROUTINE check_forward_fields(benchmark)
      LOGICAL, INTENT(IN)                                :: benchmark

      INTEGER                                            :: axis, ia, ir, p, pair, pass, repeats, s, &
                                                            variant
      REAL(dp) :: begin, checksum(2), d(2), displacement(3), dm1(2), dn(2), dp1(2), &
         drho(3, grid%ng_sphere, grid%nr, ns), ds(3, 2), dsn(3, 2), &
         dsoft(3, grid%ng_sphere, grid%nr, ns), elapsed(2), fd_error, g(3, 2), gm1(3, 2), &
         gn(3, 2), gp1(3, 2), gs2(3, 3, 2), gsn(3, 3, 2), h, normalization, &
         rho(grid%ng_sphere, grid%nr, ns), soft(grid%ng_sphere, grid%nr, ns), t(2), &
         tau(grid%ng_sphere, grid%nr, ns), tm1(2), tn(2), tp1(2), ts2(3, 2), tsn(3, 2), &
         tsoft(grid%ng_sphere, grid%nr, ns), x(3)

      DO s = 1, ns
         DO ir = 1, grid%nr
            DO ia = 1, grid%ng_sphere
               rho(ia, ir, s) = SIN(REAL(ia + ir + s, dp))
               soft(ia, ir, s) = COS(REAL(ia - ir + s, dp))
               drho(:, ia, ir, s) = rho(ia, ir, s)*[0.3_dp, -0.2_dp, 0.7_dp]
               dsoft(:, ia, ir, s) = soft(ia, ir, s)*[-0.7_dp, 0.4_dp, 0.1_dp]
            END DO
         END DO
      END DO
      tau = 0.4_dp*rho
      tsoft = 0.7_dp*soft
      DO p = 0, 105
         x = REAL(p, dp)*0.05_dp*[2.0_dp, 3.0_dp, 6.0_dp]/7.0_dp
         CALL old_fields(grid, harmonics, x, 5.0_dp, ns, rho, soft, drho, dsoft, tau, tsoft, d, g, t, ds, gs2, ts2)
         CALL new_fields(grid, harmonics, x, 5.0_dp, ns, rho, soft, drho, dsoft, tau, tsoft, dn, gn, tn, dsn, gsn, tsn)
         IF (ANY(d /= dn) .OR. ANY(g /= gn) .OR. ANY(t /= tn)) ERROR STOP 'Forward values differ'
         IF (ANY(ds /= dsn) .OR. ANY(gs2 /= gsn) .OR. ANY(ts2 /= tsn)) ERROR STOP 'Forward derivatives differ'
         CALL new_fields(grid, harmonics, x, 5.0_dp, ns, rho, soft, drho, dsoft, tau, tsoft, &
                         dn, gn, tn, dsn, gsn, tsn, calculate_spatial=.FALSE.)
         IF (ANY(d /= dn) .OR. ANY(g /= gn) .OR. ANY(t /= tn)) ERROR STOP 'Value-only forward differs'
         IF (ANY(dsn /= 0) .OR. ANY(gsn /= 0) .OR. ANY(tsn /= 0)) ERROR STOP 'Unused derivatives not zero'
      END DO
      WRITE (*, '(A,I0)') 'PASS: forward rho/gradient/tau and retained force derivatives are bitwise identical; spins=', ns
      x = [0.9_dp, -0.7_dp, 1.1_dp]
      CALL new_fields(grid, harmonics, x, 5.0_dp, ns, rho, soft, drho, dsoft, tau, tsoft, d, g, t, ds, gs2, ts2)
      normalization = MAX(1.0_dp, MAXVAL(ABS(ds)), MAXVAL(ABS(gs2)), MAXVAL(ABS(ts2)))
      fd_error = 0.0_dp
      DO axis = 1, 3
         DO p = 1, 2
            h = 1.0E-5_dp/REAL(p, dp)
            displacement = x
            displacement(axis) = x(axis) + h
            CALL new_fields(grid, harmonics, displacement, 5.0_dp, ns, rho, soft, drho, dsoft, tau, tsoft, &
                            dp1, gp1, tp1, dsn, gsn, tsn, calculate_spatial=.FALSE.)
            displacement(axis) = x(axis) - h
            CALL new_fields(grid, harmonics, displacement, 5.0_dp, ns, rho, soft, drho, dsoft, tau, tsoft, &
                            dm1, gm1, tm1, dsn, gsn, tsn, calculate_spatial=.FALSE.)
            fd_error = MAX(fd_error, MAXVAL(ABS((dp1 - dm1)/(2*h) - ds(axis, :)))/normalization, &
                           MAXVAL(ABS((gp1 - gm1)/(2*h) - gs2(:, axis, :)))/normalization, &
                           MAXVAL(ABS((tp1 - tm1)/(2*h) - ts2(axis, :)))/normalization)
         END DO
      END DO
      IF (fd_error > 2.0E-7_dp) ERROR STOP 'Value-only forward spatial finite difference failed'
      WRITE (*, '(A,ES12.4)') 'PASS: value-only forward spatial finite differences; max relative error=', fd_error
      IF (.NOT. benchmark) RETURN
      repeats = 10000
      DO pair = 1, 5
         checksum = 0.0_dp
         DO pass = 1, 2
            variant = MERGE(pass, 3 - pass, MOD(pair, 2) == 1)
            begin = omp_get_wtime()
            DO p = 1, repeats
               x = [0.9_dp + REAL(MOD(p, 17), dp)*0.013_dp, -0.7_dp, 1.1_dp]
               IF (variant == 1) THEN
                  CALL old_fields(grid, harmonics, x, 5.0_dp, ns, rho, soft, drho, dsoft, tau, tsoft, d, g, t, ds, gs2, ts2)
               ELSE
                  CALL new_fields(grid, harmonics, x, 5.0_dp, ns, rho, soft, drho, dsoft, tau, tsoft, &
                                  d, g, t, ds, gs2, ts2, calculate_spatial=.FALSE.)
               END IF
               checksum(variant) = checksum(variant) + SUM(d) + SUM(g) + SUM(t)
            END DO
            elapsed(variant) = omp_get_wtime() - begin
         END DO
         IF (checksum(1) /= checksum(2)) ERROR STOP 'Forward benchmark checksum differs'
         WRITE (*, '(A,3F12.6)') 'FORWARD seconds full/values-only and speedup: ', elapsed, elapsed(1)/elapsed(2)
      END DO
   END SUBROUTINE

   SUBROUTINE make_grid(nr, lg, lmax)
      INTEGER, INTENT(IN)                                :: nr, lg, lmax

      INTEGER                                            :: ia, ir, iso

      grid%nr = nr
      grid%ng_sphere = lebedev_grid(lg)%n
      harmonics%max_s_harm = (lmax + 1)**2
      ALLOCATE (grid%rad(nr), grid%wa(grid%ng_sphere))
      ALLOCATE (harmonics%slm(grid%ng_sphere, harmonics%max_s_harm))
      DO ir = 1, nr
         grid%rad(ir) = 0.05_dp + 4.95_dp*(REAL(ir - 1, dp)/REAL(nr - 1, dp))**2
      END DO
      grid%wa = 4.0_dp*ACOS(-1.0_dp)*lebedev_grid(lg)%w
      DO iso = 1, harmonics%max_s_harm
         DO ia = 1, grid%ng_sphere
            CALL y_lm(lebedev_grid(lg)%r(:, ia), harmonics%slm(ia, iso), indso(1, iso), indso(2, iso))
         END DO
      END DO
   END SUBROUTINE

   SUBROUTINE free_grid()
      DEALLOCATE (grid%rad, grid%wa, harmonics%slm)
   END SUBROUTINE

   SUBROUTINE check_weights()
      INTEGER                                            :: ii(4), jj(4), ni, nj, p
      LOGICAL                                            :: active_i, active_j
      REAL(dp) :: ai(grid%ng_sphere), aj(grid%ng_sphere), di(4), dj(4), gi(3, grid%ng_sphere), &
         gj(3, grid%ng_sphere), r, wi(4), wj(4), x(3)

      DO p = 0, 100 + grid%nr
         r = REAL(p, dp)*0.052_dp
         IF (p > 100) r = grid%rad(p - 100)
         x = r*[2.0_dp, 3.0_dp, 6.0_dp]/7.0_dp
         CALL old_weights(grid, harmonics, x, 5.0_dp, ii, wi, di, ni, ai, gi, active_i)
         CALL new_weights(grid, harmonics, x, 5.0_dp, jj, wj, dj, nj, aj, gj, active_j)
         IF (active_i .NEQV. active_j) ERROR STOP 'Support changed'
         IF (ni /= nj .OR. ANY(ii /= jj)) ERROR STOP 'Radial indices changed'
         IF (ANY(wi /= wj) .OR. ANY(ai /= aj)) ERROR STOP 'Value weights changed'
         IF (ANY(di /= dj) .OR. ANY(gi /= gj)) ERROR STOP 'Spatial derivatives changed'
         CALL new_weights(grid, harmonics, x, 5.0_dp, jj, wj, nradial=nj, angular_weights=aj, active=active_j)
         IF (ANY(wi /= wj) .OR. ANY(ai /= aj)) ERROR STOP 'Value-only weights changed'
         CALL new_weights(grid, harmonics, x, 5.0_dp, jj, wj, dj, nj, aj, active=active_j)
         IF (ANY(di /= dj)) ERROR STOP 'Radial-only derivatives changed'
         CALL new_weights(grid, harmonics, x, 5.0_dp, jj, wj, nradial=nj, angular_weights=aj, &
                          angular_derivative_weights=gj, active=active_j)
         IF (ANY(gi /= gj)) ERROR STOP 'Angular-only derivatives changed'
      END DO
   END SUBROUTINE

   SUBROUTINE make_cell(periodic)
      INTEGER, INTENT(IN)                                :: periodic

      cell%hmat = 0.0_dp
      cell%h_inv = 0.0_dp
      DO i = 1, 3
         cell%hmat(i, i) = 8.0_dp
         cell%h_inv(i, i) = 0.125_dp
      END DO
      cell%perd = MERGE(1, 0, periodic > 0)
      IF (periodic == 2) THEN
         cell%hmat(1, 2) = 1.6_dp
         cell%h_inv(1, 2) = -0.025_dp
      END IF
      particles(1)%r = [0.1_dp, 0.2_dp, 0.3_dp]
   END SUBROUTINE

   SUBROUTINE make_points(np)
      INTEGER, INTENT(IN)                                :: np

      INTEGER                                            :: a, p
      REAL(dp)                                           :: q

      ALLOCATE (coords(3, np), density(np, 2), gradient(np, 3, 2), kin(np, 2))
      ALLOCATE (atoms(np), local_atoms(3), first(3), last(3))
      local_atoms = [1, 2, 3]
      DO a = 1, 3
         first(a) = (a - 1)*np/3 + 1
         last(a) = a*np/3
         atoms(first(a):last(a)) = a
      END DO
      DO p = 1, np
         q = REAL((p + 1)/2, dp)
         coords(:, p) = particles(1)%r + 3.1_dp*[SIN(q), COS(0.7_dp*q), SIN(1.3_dp*q)]
         IF (MOD(p, 11) == 0) coords(:, p) = coords(:, p) + cell%hmat(:, 1)
         density(p, :) = [SIN(q), COS(q)]
         kin(p, :) = [COS(0.2_dp*q), SIN(0.3_dp*q)]
         DO a = 1, 3
            gradient(p, a, :) = [SIN(q + a), COS(q - a)]
         END DO
      END DO
   END SUBROUTINE

   SUBROUTINE free_points()
      DEALLOCATE (coords, density, gradient, kin, atoms, local_atoms, first, last)
   END SUBROUTINE

   SUBROUTINE check_target_partitions()
      ! Mimic target-atom MPI ownership and its subsequent sum without changing
      ! the actual caller loop: each partition owns one contiguous atom block.
      INTEGER                                            :: a, hi, lo, rank_first(3), rank_last(3)
      REAL(dp), ALLOCATABLE                              :: contribution(:), initial(:), total(:)

      CALL reset_potentials()
      initial = pack_potentials()
      total = initial
      DO a = 1, 3
         lo = first(a)
         hi = last(a)
         rank_first = 1
         rank_last = 0
         rank_last(a) = hi - lo + 1
         vh = 0; vs = 0; gh = 0; gs = 0; th = 0; ts = 0
         CALL parallel_project(cell, particles, grid, harmonics, ns, flags, &
                               coords(:, lo:hi), atoms(lo:hi), local_atoms(a:a), rank_first, rank_last, &
                               density(lo:hi, :), gradient(lo:hi, :, :), kin(lo:hi, :), 5.0_dp, vh, vs, gh, gs, th, ts)
         contribution = pack_potentials()
         total = total + contribution
      END DO
      IF (MAXVAL(ABS(total - reference)) > 2.0E-13_dp*MAX(1.0_dp, MAXVAL(ABS(reference)))) &
         ERROR STOP 'Target ownership partition sum differs from original'
   END SUBROUTINE

   SUBROUTINE allocate_potentials()
      INTEGER                                            :: na, nr

      na = grid%ng_sphere
      nr = grid%nr
      ALLOCATE (vh(na, nr, ns), vs(na, nr, ns), th(na, nr, ns), ts(na, nr, ns), gh(3, na, nr, ns), gs(3, na, nr, ns))
   END SUBROUTINE

   SUBROUTINE free_potentials()
      DEALLOCATE (vh, vs, th, ts, gh, gs)
   END SUBROUTINE

   SUBROUTINE reset_potentials()
      vh = 0.031_dp
      vs = -0.072_dp
      th = 0.11_dp
      ts = -0.03_dp
      gh = 0.07_dp
      gs = -0.012_dp
   END SUBROUTINE

   FUNCTION pack_potentials() RESULT(p)
      REAL(dp), ALLOCATABLE                              :: p(:)

      p = [RESHAPE(vh, [SIZE(vh)]), RESHAPE(vs, [SIZE(vs)]), RESHAPE(th, [SIZE(th)]), &
           RESHAPE(ts, [SIZE(ts)]), RESHAPE(gh, [SIZE(gh)]), RESHAPE(gs, [SIZE(gs)])]
   END FUNCTION

   SUBROUTINE project_variant(which)
      INTEGER, INTENT(IN)                                :: which

      SELECT CASE (which)
      CASE (1)
         CALL old_project(cell, particles, grid, harmonics, ns, flags, coords, atoms, local_atoms, first, last, &
                          density, gradient, kin, 5.0_dp, vh, vs, gh, gs, th, ts)
      CASE (2)
         CALL value_project(cell, particles, grid, harmonics, ns, flags, coords, atoms, local_atoms, first, last, &
                            density, gradient, kin, 5.0_dp, vh, vs, gh, gs, th, ts)
      CASE (3)
         CALL screened_project(cell, particles, grid, harmonics, ns, flags, coords, atoms, local_atoms, first, last, &
                               density, gradient, kin, 5.0_dp, vh, vs, gh, gs, th, ts)
      CASE (4)
         CALL parallel_project(cell, particles, grid, harmonics, ns, flags, coords, atoms, local_atoms, first, last, &
                               density, gradient, kin, 5.0_dp, vh, vs, gh, gs, th, ts)
      END SELECT
   END SUBROUTINE

   SUBROUTINE check_transpose()
      USE parallel_kernel, ONLY: add_gapw_atom_grid_interpolation_adjoint
      INTEGER                                            :: ia, ir, s
      REAL(dp) :: d(2), da(2), drho(3, grid%ng_sphere, grid%nr, ns), ds(3, 2), &
         dsoft(3, grid%ng_sphere, grid%nr, ns), g(3, 2), ga(3, 2), gs2(3, 3, 2), lhs, &
         rho(grid%ng_sphere, grid%nr, ns), rhs, soft(grid%ng_sphere, grid%nr, ns), t(2), ta(2), &
         tau(grid%ng_sphere, grid%nr, ns), tsoft(grid%ng_sphere, grid%nr, ns), tss(3, 2), x(3)

      DO s = 1, ns
         DO ir = 1, grid%nr
            DO ia = 1, grid%ng_sphere
               rho(ia, ir, s) = SIN(REAL(ia + ir + s, dp))
               soft(ia, ir, s) = COS(REAL(ia - ir + s, dp))
               drho(:, ia, ir, s) = rho(ia, ir, s)*[0.3_dp, -0.2_dp, 0.7_dp]
               dsoft(:, ia, ir, s) = soft(ia, ir, s)*[-0.7_dp, 0.4_dp, 0.1_dp]
            END DO
         END DO
      END DO
      tau = 0.4_dp*rho
      tsoft = 0.7_dp*soft
      da = [0.3_dp, -0.2_dp]
      ga = RESHAPE([0.1_dp, 0.7_dp, -0.3_dp, 0.8_dp, -0.2_dp, 0.3_dp], [3, 2])
      ta = [0.9_dp, -0.6_dp]
      x = [0.9_dp, -0.7_dp, 1.1_dp]
      CALL old_fields(grid, harmonics, x, 5.0_dp, ns, rho, soft, drho, dsoft, tau, tsoft, d, g, t, ds, gs2, tss)
      vh = 0; vs = 0; gh = 0; gs = 0; th = 0; ts = 0
      CALL add_gapw_atom_grid_interpolation_adjoint(grid, harmonics, x, 5.0_dp, ns, da, ga, ta, vh, vs, gh, gs, th, ts)
      lhs = SUM(d*da) + SUM(g*ga) + SUM(t*ta)
      rhs = SUM(vh*rho) - SUM(vs*soft) + SUM(gh*drho) - SUM(gs*dsoft) + SUM(th*tau) - SUM(ts*tsoft)
      IF (ABS(lhs - rhs) > 2.0E-13_dp*MAX(1.0_dp, ABS(lhs))) ERROR STOP 'Adjoint dot-product identity failed'
   END SUBROUTINE
END PROGRAM

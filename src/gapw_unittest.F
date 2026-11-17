!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM gapw_unittest
   USE basis_set_types,                 ONLY: allocate_gto_basis_set,&
                                              create_primitive_basis_set,&
                                              deallocate_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE gapw_1c_basis_set,               ONLY: create_1c_basis
   USE hartree_local_methods,           ONLY: gapw_background_potential
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: init_lebedev_grids,&
                                              lebedev_grid
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE mathlib,                         ONLY: inv_3x3
   USE orbital_pointers,                ONLY: indso,&
                                              init_orbital_pointers
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE qs_vxc_atom,                     ONLY: add_gapw_atom_grid_interpolation_adjoint,&
                                              atom_grid_image_bounds,&
                                              interpolate_gapw_atom_grid_fields
   USE spherical_harmonics,             ONLY: y_lm

   IMPLICIT NONE

   INTEGER, PARAMETER :: nr = 1601, ng = 32
   REAL(dp), PARAMETER :: alpha = 1.4_dp, beta = 0.55_dp, box = 12.0_dp, dr = 0.00625_dp
   REAL(dp), PARAMETER :: charges(2) = [-1.0_dp, 0.3_dp]
   REAL(dp) :: rad2(nr), weights(nr), delta_rho(nr), potential(nr), shifted(nr)
   REAL(dp) :: q, g2, periodic_difference, local_difference, correction, error, h, eplus, eminus
   INTEGER :: i, ix, iy, iz, io_unit

   CALL cp2k_init()
   io_unit = cp_logger_get_default_io_unit()
   CALL check_1c_basis()
   CALL check_atom_grid_adjoint()
   CALL check_atom_grid_images()
   DO i = 1, nr
      rad2(i) = ((i - 1)*dr)**2
   END DO
   weights = 2.0_dp*dr/3.0_dp
   weights(2:nr - 1:2) = 4.0_dp*dr/3.0_dp
   weights(1) = dr/3.0_dp
   weights(nr) = dr/3.0_dp
   weights = weights*rad2

   DO i = 1, SIZE(charges)
      q = charges(i)
      ! Electron-minus-nuclear l=0 density. Both normalized Gaussians carry q,
      ! so their difference has zero monopole and a nonzero second moment.
      delta_rho = q*SQRT(fourpi)*((alpha/pi)**1.5_dp*EXP(-alpha*rad2) - &
                                  (beta/pi)**1.5_dp*EXP(-beta*rad2))
      potential = gapw_background_potential(-q, box**3, rad2)
      correction = SUM(potential*delta_rho*weights)

      ! Independent reference: exact reciprocal-space Coulomb energies of the
      ! two periodic Gaussians with G=0 omitted, minus their free-space energies.
      periodic_difference = 0.0_dp
      DO iz = -ng, ng
         DO iy = -ng, ng
            DO ix = -ng, ng
               IF (ix == 0 .AND. iy == 0 .AND. iz == 0) CYCLE
               g2 = (2.0_dp*pi/box)**2*REAL(ix*ix + iy*iy + iz*iz, dp)
               periodic_difference = periodic_difference + &
                                     2.0_dp*pi*q*q/box**3/g2*(EXP(-g2/(2.0_dp*alpha)) - EXP(-g2/(2.0_dp*beta)))
            END DO
         END DO
      END DO
      local_difference = q*q*(SQRT(alpha) - SQRT(beta))/SQRT(2.0_dp*pi)
      error = ABS(periodic_difference - local_difference - correction)
      IF (io_unit > 0) WRITE (io_unit, '(A,F8.3,ES16.7)') "GAPW periodic Gaussian background:", q, error
      IF (error > 1.e-12_dp) ERROR STOP "GAPW background disagrees with periodic Gaussian electrostatics"

      ! The explicit volume derivative contributes E_background to each
      ! diagonal component of -dE/dstrain (local radial functions held fixed).
      h = 1.e-5_dp
      shifted = gapw_background_potential(-q, box**3*(1.0_dp + h), rad2)
      eplus = SUM(shifted*delta_rho*weights)
      shifted = gapw_background_potential(-q, box**3*(1.0_dp - h), rad2)
      eminus = SUM(shifted*delta_rho*weights)
      IF (ABS(-(eplus - eminus)/(2.0_dp*h) - correction) > 1.e-11_dp) THEN
         ERROR STOP "GAPW background volume derivative failed"
      END IF
   END DO
   potential = gapw_background_potential(0.0_dp, box**3, rad2)
   IF (ANY(potential /= 0.0_dp)) ERROR STOP "Neutral cell has a background potential"
   CALL cp2k_finalize()
CONTAINS

! **************************************************************************************************
!> \brief Check that one-centre basis extension retains diffuse exponents without duplicates.
! **************************************************************************************************
   SUBROUTINE check_1c_basis()
      INTEGER, PARAMETER                                 :: expected(4) = [3, 4, 5, 7]
      REAL(KIND=dp), PARAMETER :: spacing(4) = [2.5_dp, 2.0_dp, 1.75_dp, 1.5_dp]

      INTEGER                                            :: level
      REAL(KIND=dp)                                      :: diffuse, ratio
      TYPE(gto_basis_set_type), POINTER                  :: extended, orbitals, seed

      NULLIFY (extended, orbitals, seed)
      CALL init_orbital_pointers(5)
      CALL init_spherical_harmonics(5, 0)
      DO level = 1, 4
         ratio = SPACING(level)
         CALL allocate_gto_basis_set(seed)
         seed%name = 'two-s-primitives'
         seed%norm_type = 2
         seed%nset = 1
         ALLOCATE (seed%lmin(1), seed%lmax(1), seed%npgf(1), seed%zet(2, 1))
         seed%lmin = 0
         seed%lmax = 0
         seed%npgf = 2
         ! The second extension candidate duplicates the second soft exponent.
         ! More diffuse, well-separated candidates must still be admitted.
         seed%zet(:, 1) = [1.0_dp, 1.0_dp/SQRT(ratio)]
         CALL create_primitive_basis_set(seed, orbitals)
         orbitals%set_radius = 0.0_dp
         CALL create_1c_basis(orbitals, orbitals, extended, level)
         diffuse = 1.0_dp/(ratio*SQRT(ratio))
         IF (io_unit > 0) WRITE (io_unit, '(A,I2,A,I3)') &
            'GAPW extension level ', level, ': added s primitives ', extended%npgf(2)
         IF (extended%npgf(2) /= expected(level)) ERROR STOP 'GAPW extension lost admissible diffuse exponents'
         IF (MINVAL(ABS(extended%zet(:extended%npgf(2), 2) - diffuse)) > 1.e-12_dp) THEN
            ERROR STOP 'GAPW extension is missing the first admissible diffuse exponent'
         END IF
         IF (MINVAL(ABS(extended%zet(:extended%npgf(2), 2) - seed%zet(2, 1))) < 1.e-12_dp) THEN
            ERROR STOP 'GAPW extension duplicated an existing soft exponent'
         END IF
         CALL deallocate_gto_basis_set(extended)
         CALL deallocate_gto_basis_set(orbitals)
         CALL deallocate_gto_basis_set(seed)
      END DO
   END SUBROUTINE check_1c_basis

! **************************************************************************************************
!> \brief Check the hard-minus-soft density, gradient and tau interpolation transpose.
! **************************************************************************************************
   SUBROUTINE check_atom_grid_adjoint()
      INTEGER, PARAMETER                                 :: nrad = 7
      REAL(dp), PARAMETER :: radii(5) = [0.0_dp, 0.06_dp, 1.5_dp, 4.99_dp, 5.01_dp]

      INTEGER                                            :: a, ir, iso, na, ns, order, sample, spin
      REAL(dp)                                           :: d(2), da(2), ds(3, 2), g(3, 2), &
                                                            ga(3, 2), gs(3, 3, 2), lhs(3), &
                                                            max_error, rhs(3), t(2), ta(2), &
                                                            ts(3, 2), x(3)
      REAL(dp), ALLOCATABLE :: grad_h(:, :, :, :), grad_s(:, :, :, :), rho_h(:, :, :), &
         rho_s(:, :, :), tau_h(:, :, :), tau_s(:, :, :), v_h(:, :, :), v_s(:, :, :), &
         vg_h(:, :, :, :), vg_s(:, :, :, :), vt_h(:, :, :), vt_s(:, :, :)
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      CALL init_lebedev_grids()
      ALLOCATE (grid, harmonics)
      na = lebedev_grid(5)%n
      grid%nr = nrad
      grid%ng_sphere = na
      harmonics%max_s_harm = 25
      ALLOCATE (grid%rad(nrad), grid%wa(na), harmonics%slm(na, 25))
      grid%rad = [(0.05_dp + 4.95_dp*(REAL(ir - 1, dp)/REAL(nrad - 1, dp))**2, ir=1, nrad)]
      grid%wa = fourpi*lebedev_grid(5)%w
      DO iso = 1, 25
         DO a = 1, na
            CALL y_lm(lebedev_grid(5)%r(:, a), harmonics%slm(a, iso), indso(1, iso), indso(2, iso))
         END DO
      END DO
      ALLOCATE (rho_h(na, nrad, 2), rho_s(na, nrad, 2), tau_h(na, nrad, 2), tau_s(na, nrad, 2), &
                v_h(na, nrad, 2), v_s(na, nrad, 2), vt_h(na, nrad, 2), vt_s(na, nrad, 2), &
                grad_h(3, na, nrad, 2), grad_s(3, na, nrad, 2), vg_h(3, na, nrad, 2), vg_s(3, na, nrad, 2))
      DO spin = 1, 2
         DO ir = 1, nrad
            DO a = 1, na
               rho_h(a, ir, spin) = SIN(REAL(a + ir + spin, dp))
               rho_s(a, ir, spin) = COS(REAL(a - ir + spin, dp))
               grad_h(:, a, ir, spin) = rho_h(a, ir, spin)*[0.3_dp, -0.2_dp, 0.7_dp]
               grad_s(:, a, ir, spin) = rho_s(a, ir, spin)*[-0.7_dp, 0.4_dp, 0.1_dp]
            END DO
         END DO
      END DO
      tau_h(:, :, :) = 0.4_dp*rho_h
      tau_s(:, :, :) = 0.7_dp*rho_s
      da = [0.3_dp, -0.2_dp]
      ga = RESHAPE([0.1_dp, 0.7_dp, -0.3_dp, 0.8_dp, -0.2_dp, 0.3_dp], [3, 2])
      ta = [0.9_dp, -0.6_dp]
      max_error = 0.0_dp
      DO order = 1, 2
         DO ns = 1, 2
            DO sample = 1, SIZE(radii)
               x = radii(sample)*[2.0_dp, 3.0_dp, 6.0_dp]/7.0_dp
               CALL interpolate_gapw_atom_grid_fields(grid, harmonics, x, 5.0_dp, ns, &
                                                      rho_h, rho_s, grad_h, grad_s, tau_h, tau_s, &
                                                      d, g, t, ds, gs, ts, calculate_spatial=.FALSE.)
               v_h = 0.0_dp
               v_s = 0.0_dp
               vg_h = 0.0_dp
               vg_s = 0.0_dp
               vt_h = 0.0_dp
               vt_s = 0.0_dp
               CALL add_gapw_atom_grid_interpolation_adjoint(grid, harmonics, x, 5.0_dp, ns, &
                                                             da, ga, ta, v_h, v_s, vg_h, vg_s, vt_h, vt_s)
               lhs = [SUM(d*da), SUM(g*ga), SUM(t*ta)]
               rhs = [SUM(v_h*rho_h) - SUM(v_s*rho_s), SUM(vg_h*grad_h) - SUM(vg_s*grad_s), &
                      SUM(vt_h*tau_h) - SUM(vt_s*tau_s)]
               IF (sample == 3 .AND. MINVAL(ABS(lhs)) < 1.e-8_dp) ERROR STOP 'Atom-grid fixture is inactive'
               IF (sample == 1 .OR. sample == 5) THEN
                  IF (ANY(lhs /= 0.0_dp) .OR. ANY(rhs /= 0.0_dp)) ERROR STOP 'Inactive atom-grid point contributed'
               END IF
               max_error = MAX(max_error, MAXVAL(ABS(lhs - rhs)/MAX(1.0_dp, ABS(lhs), ABS(rhs))))
            END DO
         END DO
         grid%rad = grid%rad(nrad:1:-1)
      END DO
      IF (io_unit > 0) WRITE (io_unit, '(A,ES16.7)') 'GAPW atom-grid adjoint error: ', max_error
      IF (max_error > 2.e-13_dp) ERROR STOP 'GAPW interpolation transpose identity failed'
      DEALLOCATE (grid%rad, grid%wa, harmonics%slm)
      DEALLOCATE (grid, harmonics)
   END SUBROUTINE check_atom_grid_adjoint

! **************************************************************************************************
!> \brief Ensure periodic-image screening retains every image within the Cartesian cutoff.
! **************************************************************************************************
   SUBROUTINE check_atom_grid_images()
      INTEGER, PARAMETER                                 :: masks(3) = [0, 3, 7]

      INTEGER                                            :: center(3), hi(3), i1, i2, i3, lo(3), &
                                                            mask, p, sample, shell(3), shift(3)
      REAL(dp)                                           :: cutoff, fractional(3), point(3)
      TYPE(cell_type), POINTER                           :: cell

      ALLOCATE (cell)
      cell%hmat = RESHAPE([8.0_dp, 0.0_dp, 0.0_dp, 3.2_dp, 8.0_dp, 0.0_dp, 0.0_dp, 1.6_dp, 8.0_dp], [3, 3])
      cell%h_inv = inv_3x3(cell%hmat)
      DO mask = 1, SIZE(masks)
         cell%perd = [(MERGE(1, 0, BTEST(masks(mask), p)), p=0, 2)]
         DO sample = 1, 3
            fractional = [0.5_dp, -0.5_dp, 0.2_dp] + &
                         REAL(sample - 1, dp)*cell%perd*[100.0_dp, -51.0_dp, 77.0_dp]
            point = MATMUL(cell%hmat, fractional)
            center = cell%perd*NINT(fractional)
            cutoff = NORM2(point - MATMUL(cell%hmat, REAL(center, dp)))
            shell = cell%perd*(CEILING(cutoff*SQRT(SUM(cell%h_inv**2, DIM=2))) + 1)
            CALL atom_grid_image_bounds(cell, fractional, cutoff, shell, lo, hi)
            IF (ANY(lo < center - shell) .OR. ANY(hi > center + shell)) ERROR STOP 'Image bounds expanded'
            DO i3 = center(3) - shell(3), center(3) + shell(3)
               DO i2 = center(2) - shell(2), center(2) + shell(2)
                  DO i1 = center(1) - shell(1), center(1) + shell(1)
                     shift = [i1, i2, i3]
                     IF (NORM2(point - MATMUL(cell%hmat, REAL(shift, dp))) > cutoff) CYCLE
                     IF (ANY(shift < lo) .OR. ANY(shift > hi)) ERROR STOP 'Image within support omitted'
                  END DO
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE (cell)
   END SUBROUTINE check_atom_grid_images
END PROGRAM gapw_unittest

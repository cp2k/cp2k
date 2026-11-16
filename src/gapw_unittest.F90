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
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE gapw_1c_basis_set,               ONLY: create_1c_basis
   USE hartree_local_methods,           ONLY: gapw_background_potential
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE orbital_pointers,                ONLY: init_orbital_pointers
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics

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
END PROGRAM gapw_unittest

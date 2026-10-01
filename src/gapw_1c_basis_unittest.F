!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM gapw_1c_basis_unittest
   USE basis_set_types,                 ONLY: allocate_gto_basis_set,&
                                              create_primitive_basis_set,&
                                              deallocate_gto_basis_set,&
                                              gto_basis_set_type
   USE gapw_1c_basis_set,               ONLY: create_1c_basis
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE orbital_pointers,                ONLY: init_orbital_pointers
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics

   IMPLICIT NONE
   INTEGER :: level
   INTEGER, PARAMETER :: expected(4) = [3, 4, 5, 7]
   REAL(KIND=dp) :: ratio, diffuse
   REAL(KIND=dp), PARAMETER :: spacing(4) = [2.5_dp, 2.0_dp, 1.75_dp, 1.5_dp]
   TYPE(gto_basis_set_type), POINTER :: seed => NULL(), orbitals => NULL(), extended => NULL()

   CALL cp2k_init()
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
      PRINT '(A,I2,A,I3)', 'GAPW extension level ', level, ': added s primitives ', extended%npgf(2)
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
   CALL cp2k_finalize()
END PROGRAM gapw_1c_basis_unittest

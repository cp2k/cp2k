!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM skala_periodic_partition_unittest
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathlib,                         ONLY: inv_3x3
   USE skala_periodic_partition,        ONLY: build_periodic_atom_image_layout,&
                                              periodic_atom_image_partition,&
                                              periodic_atom_image_partition_from_layout,&
                                              skala_periodic_partition_type

   IMPLICIT NONE
   TYPE(cell_type), POINTER :: cell, perturbed
   TYPE(skala_periodic_partition_type) :: layout
   REAL(dp) :: atoms(3, 2), trial(3, 2), hmat(3, 3), point(3), shifted(3)
   REAL(dp) :: value, plus, minus, direct, total, datom(3, 2), dstrain(3, 3), mu
   REAL(dp) :: unity_error, derivative_error, translation_error
   REAL(dp), PARAMETER :: h = 2.e-6_dp
   INTEGER :: mask, periodicity(3), a, d, j, k, l, TARGET, sample

   CALL cp2k_init()
   ALLOCATE (cell, perturbed)
   hmat = RESHAPE([4.0_dp, 0.1_dp, 0.2_dp, 0.7_dp, 5.0_dp, 0.3_dp, -0.2_dp, 0.4_dp, 6.0_dp], [3, 3])
   atoms(:, 1) = [0.2_dp, 0.3_dp, -0.1_dp]
   atoms(:, 2) = [1.8_dp, 2.2_dp, 1.3_dp]
   unity_error = 0.0_dp
   derivative_error = 0.0_dp
   translation_error = 0.0_dp
   DO mask = 0, 7
      periodicity = [MERGE(1, 0, BTEST(mask, 0)), MERGE(1, 0, BTEST(mask, 1)), MERGE(1, 0, BTEST(mask, 2))]
      cell%hmat = hmat
      cell%h_inv = inv_3x3(hmat)
      cell%perd = periodicity
      perturbed%perd = periodicity
      CALL build_periodic_atom_image_layout(atoms, cell, 1, periodicity, layout, TARGET)
      DO sample = 1, 3
         point = [1.1_dp, 1.3_dp, 1.7_dp]
         IF (sample == 2) point = MATMUL(hmat, [0.5_dp, 0.2_dp, 0.1_dp])
         IF (sample == 3) point = MATMUL(hmat, [1.5_dp, 0.2_dp, 0.1_dp]) + atoms(:, 1)
         CALL periodic_atom_image_partition_from_layout(point, layout, TARGET, value, datom, dstrain)
         CALL periodic_atom_image_partition(point, atoms, cell, 1, direct, image_periodicity=periodicity)
         IF (ABS(direct - value) > 2.e-14_dp) ERROR STOP 'Cached partition changed the value'
         IF (mask == 0) THEN
            mu = (NORM2(point - atoms(:, 1)) - NORM2(point - atoms(:, 2)))/NORM2(atoms(:, 1) - atoms(:, 2))
            DO j = 1, 3
               mu = 1.5_dp*mu - 0.5_dp*mu**3
            END DO
            IF (ABS(value - 0.5_dp*(1.0_dp - mu)) > 2.e-14_dp) ERROR STOP 'Nonperiodic Becke weight changed'
            IF (MAXVAL(ABS(dstrain)) > 2.e-14_dp) ERROR STOP 'Molecular partition has an image-strain response'
         END IF
         total = 0.0_dp
         DO a = 1, 2
            DO l = -3*periodicity(3), 3*periodicity(3)
               DO k = -3*periodicity(2), 3*periodicity(2)
                  DO j = -3*periodicity(1), 3*periodicity(1)
                     shifted = point - MATMUL(hmat, REAL([j, k, l], dp))
                     CALL periodic_atom_image_partition(shifted, atoms, cell, a, direct, image_periodicity=periodicity)
                     total = total + direct
                  END DO
               END DO
            END DO
         END DO
         unity_error = MAX(unity_error, ABS(total - 1.0_dp))
         trial = atoms
         trial(:, 2) = atoms(:, 2) + MATMUL(hmat, REAL([1, -1, 2]*periodicity, dp))
         CALL periodic_atom_image_partition(point, trial, cell, 1, direct, image_periodicity=periodicity)
         translation_error = MAX(translation_error, ABS(direct - value))
         DO a = 1, 2
            DO d = 1, 3
               trial = atoms
               trial(d, a) = atoms(d, a) + h
               CALL periodic_atom_image_partition(point, trial, cell, 1, plus, image_periodicity=periodicity)
               trial(d, a) = atoms(d, a) - h
               CALL periodic_atom_image_partition(point, trial, cell, 1, minus, image_periodicity=periodicity)
               derivative_error = MAX(derivative_error, ABS((plus - minus)/(2*h) - datom(d, a)))
            END DO
         END DO
         DO d = 1, 3
            DO j = 1, 3
               perturbed%hmat = hmat
               perturbed%hmat(d, :) = hmat(d, :) + h*hmat(j, :)
               perturbed%h_inv = inv_3x3(perturbed%hmat)
               CALL periodic_atom_image_partition(point, atoms, perturbed, 1, plus, image_periodicity=periodicity)
               perturbed%hmat(d, :) = hmat(d, :) - h*hmat(j, :)
               perturbed%h_inv = inv_3x3(perturbed%hmat)
               CALL periodic_atom_image_partition(point, atoms, perturbed, 1, minus, image_periodicity=periodicity)
               derivative_error = MAX(derivative_error, ABS((plus - minus)/(2*h) - dstrain(d, j)))
            END DO
         END DO
      END DO
   END DO
   WRITE (*, '(A,ES16.7)') 'PERIODIC PARTITION UNITY ERROR ', unity_error
   WRITE (*, '(A,ES16.7)') 'PERIODIC PARTITION TRANSLATION ERROR ', translation_error
   WRITE (*, '(A,ES16.7)') 'PERIODIC PARTITION ATOM/STRAIN ERROR ', derivative_error
   IF (unity_error > 2.e-13_dp) ERROR STOP 'Periodic atom-image weights do not sum to one'
   IF (translation_error > 2.e-13_dp) ERROR STOP 'Moving an atom by a lattice vector changed its image family'
   IF (derivative_error > 3.e-8_dp) ERROR STOP 'Periodic partition derivatives fail finite differences'
   DEALLOCATE (cell, perturbed)
   CALL cp2k_finalize()
END PROGRAM skala_periodic_partition_unittest

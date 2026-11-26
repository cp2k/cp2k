!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM kpsym_unittest
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cryssym,                         ONLY: crys_sym_gen,&
                                              csym_type,&
                                              kpoint_gen,&
                                              release_csym_type
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE physcon,                         ONLY: bohr
   USE qs_kpoint_density_state,         ONLY: kpoint_density_metric_transport

   IMPLICIT NONE
   REAL(dp), PARAMETER :: ice_ix(3, 36) = RESHAPE([ &
                                                  -0.006830_dp, 0.340009_dp, 0.190161_dp, &
                                                  0.113981_dp, 0.164412_dp, 0.283330_dp, &
                                                  0.309270_dp, 0.370346_dp, 0.106284_dp, &
                                                  0.159991_dp, 0.493170_dp, 0.440161_dp, &
                                                  0.335588_dp, 0.613981_dp, 0.533330_dp, &
                                                  0.129654_dp, 0.809270_dp, 0.356284_dp, &
                                                  0.840009_dp, 0.506830_dp, -0.059839_dp, &
                                                  0.664412_dp, 0.386019_dp, 0.033330_dp, &
                                                  0.870346_dp, 0.190730_dp, 0.856284_dp, &
                                                  0.506830_dp, 0.840009_dp, 1.059839_dp, &
                                                  0.386019_dp, 0.664412_dp, 0.966670_dp, &
                                                  0.190730_dp, 0.870346_dp, 0.143716_dp, &
                                                  0.493170_dp, 0.159991_dp, 0.559839_dp, &
                                                  0.613981_dp, 0.335588_dp, 0.466670_dp, &
                                                  0.809270_dp, 0.129654_dp, 0.643716_dp, &
                                                  0.370346_dp, 0.309270_dp, -0.106284_dp, &
                                                  0.659991_dp, 1.006830_dp, 0.309839_dp, &
                                                  0.835588_dp, 0.886019_dp, 0.216670_dp, &
                                                  0.629654_dp, 0.690730_dp, 0.393716_dp, &
                                                  1.006830_dp, 0.659991_dp, 0.690161_dp, &
                                                  0.886019_dp, 0.835588_dp, 0.783330_dp, &
                                                  0.690730_dp, 0.629654_dp, 0.606284_dp, &
                                                  0.340009_dp, -0.006830_dp, 0.809839_dp, &
                                                  0.164412_dp, 0.113981_dp, 0.716670_dp, &
                                                  0.114088_dp, 0.310030_dp, 0.267366_dp, &
                                                  0.189970_dp, 0.614088_dp, 0.517366_dp, &
                                                  0.810030_dp, 0.385912_dp, 0.017365_dp, &
                                                  0.385912_dp, 0.810030_dp, 0.982634_dp, &
                                                  0.614088_dp, 0.189970_dp, 0.482634_dp, &
                                                  0.689970_dp, 0.885912_dp, 0.232634_dp, &
                                                  0.885912_dp, 0.689970_dp, 0.767366_dp, &
                                                  0.310030_dp, 0.114088_dp, 0.732634_dp, &
                                                  0.404256_dp, 0.404256_dp, 0.000000_dp, &
                                                  0.095744_dp, 0.904256_dp, 0.250000_dp, &
                                                  0.904256_dp, 0.095744_dp, 0.750000_dp, &
                                                  0.595744_dp, 0.595744_dp, 0.500000_dp], [3, 36])
   REAL(dp) :: cell(3, 3), coords(3, 36), shifted(3, 36)
   INTEGER :: i, io_unit
   INTEGER :: types(36)

   CALL cp2k_init()
   CALL check_metric_transport()
   io_unit = cp_logger_get_default_io_unit()
   cell = 0.0_dp
   cell(1, 1) = 6.784903_dp*bohr
   cell(2, 2) = 6.784903_dp*bohr
   cell(3, 3) = 6.807916_dp*bohr
   types(:24) = 1
   types(25:) = 2
   coords = ice_ix

   ! Rounded coordinates put seven operations within a few ulps of the boundary.
   CALL check_mesh(coords, cell, 1.0e-6_dp, 6, .FALSE., .TRUE.)
   CALL check_mesh(coords, cell, 1.0e-7_dp, 14, .FALSE., .TRUE.)
   shifted = coords + SPREAD([0.137_dp, -0.231_dp, 0.319_dp], 2, 36)
   CALL check_mesh(shifted, cell, 1.0e-6_dp, 6, .FALSE., .TRUE.)
   CALL check_mesh(coords, 1.3_dp*cell, 1.0e-6_dp, 6, .FALSE., .TRUE.)
   DO i = 1, 24
      shifted(:, i) = coords(:, 25 - i)
   END DO
   DO i = 25, 36
      shifted(:, i) = coords(:, 61 - i)
   END DO
   CALL check_mesh(shifted, cell, 1.0e-6_dp, 6, .FALSE., .TRUE.)

   ! Genuine tolerance-induced nonclosure must still trigger the orbit guard.
   coords = ice_ix
   coords(3, 27) = coords(3, 27) + 1.0e-6_dp
   coords(3, 33) = coords(3, 33) + 0.75e-6_dp
   coords(3, 36) = coords(3, 36) - 0.75e-6_dp
   CALL check_mesh(coords, cell, 1.0e-6_dp, 14, .TRUE., .FALSE.)
   IF (io_unit > 0) WRITE (io_unit, *) "K290 tolerance-boundary and orbit-guard tests passed"
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Restart transport must preserve a common spatial action and nonorthogonal AO gauges.
! **************************************************************************************************
   SUBROUTINE check_metric_transport()
      COMPLEX(dp) :: basis(2, 2), current(2, 2), expected(2, 2), inverse(2, 2), previous(2, 2), &
         rotation(2, 2), scaled_current(2, 2), scaled_previous(2, 2)
      COMPLEX(dp), ALLOCATABLE                           :: transport(:, :)
      INTEGER                                            :: info

      previous = CMPLX(RESHAPE([2.0_dp, 0.4_dp, 0.4_dp, 2.0_dp], [2, 2]), 0.0_dp, dp)
      current = CMPLX(RESHAPE([1.4_dp, 0.6_dp, 0.6_dp, 1.4_dp], [2, 2]), 0.0_dp, dp)
      rotation = CMPLX(RESHAPE([0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], [2, 2]), 0.0_dp, dp)
      expected(1, 1) = 0.5_dp*(SQRT(1.2_dp) + SQRT(2.0_dp))
      expected(1, 2) = 0.5_dp*(SQRT(1.2_dp) - SQRT(2.0_dp))
      expected(2, 1) = expected(1, 2)
      expected(2, 2) = expected(1, 1)
      CALL kpoint_density_metric_transport(previous, current, transport, info)
      IF (info /= 0) ERROR STOP 'Positive restart metric pair rejected'
      IF (MAXVAL(ABS(transport - expected)) > 1.e-12_dp) ERROR STOP 'Incorrect relative metric square root'
      IF (MAXVAL(ABS(MATMUL(CONJG(TRANSPOSE(transport)), MATMUL(current, transport)) - previous)) > 1.e-12_dp) THEN
         ERROR STOP 'Restart transport changed the charge metric'
      END IF
      IF (MAXVAL(ABS(MATMUL(transport, rotation) - MATMUL(rotation, transport))) > 1.e-12_dp) THEN
         ERROR STOP 'Restart transport broke a common spatial symmetry'
      END IF
      ! A complex nonunitary AO change makes both the relative root and spatial action non-Hermitian.
      basis = CMPLX(0.0_dp, 0.0_dp, dp)
      basis(1, 1) = 1.2_dp
      basis(2, 2) = 0.7_dp
      basis(1, 2) = CMPLX(0.0_dp, 0.2_dp, dp)
      inverse = CMPLX(0.0_dp, 0.0_dp, dp)
      inverse(1, 1) = 1.0_dp/1.2_dp
      inverse(2, 2) = 1.0_dp/0.7_dp
      inverse(1, 2) = -basis(1, 2)/(1.2_dp*0.7_dp)
      scaled_previous = MATMUL(CONJG(TRANSPOSE(basis)), MATMUL(previous, basis))
      scaled_current = MATMUL(CONJG(TRANSPOSE(basis)), MATMUL(current, basis))
      CALL kpoint_density_metric_transport(scaled_previous, scaled_current, transport, info)
      IF (info /= 0) ERROR STOP 'Nonorthogonal AO metric pair rejected'
      IF (MAXVAL(ABS(transport - MATMUL(inverse, MATMUL(expected, basis)))) > 1.e-12_dp) THEN
         ERROR STOP 'Restart transport depends on the AO coordinate frame'
      END IF
      CALL kpoint_density_metric_transport(previous, previous, transport, info)
      IF (info /= 0) ERROR STOP 'Unchanged restart metric rejected'
      transport(1, 1) = transport(1, 1) - 1.0_dp
      transport(2, 2) = transport(2, 2) - 1.0_dp
      IF (MAXVAL(ABS(transport)) > 1.e-12_dp) ERROR STOP 'Unchanged metric moved the restart state'
      previous(1, 1) = -1.0_dp
      CALL kpoint_density_metric_transport(previous, current, transport, info)
      IF (info == 0) ERROR STOP 'Nonpositive restart metric accepted'
   END SUBROUTINE check_metric_transport

! **************************************************************************************************
!> \brief Check mesh reduction and atomic operation consistency without an SCF calculation.
!> \param x fractional atomic coordinates
!> \param h cell matrix in bohr
!> \param delta K290 symmetry tolerance
!> \param expected expected irreducible k-point count
!> \param fallback expected inversion fallback
!> \param closed whether to check closure of the retained atomic operations
! **************************************************************************************************
   SUBROUTINE check_mesh(x, h, delta, expected, fallback, closed)
      REAL(dp), INTENT(IN)                               :: x(:, :), h(3, 3), delta
      INTEGER, INTENT(IN)                                :: expected
      LOGICAL, INTENT(IN)                                :: fallback, closed

      INTEGER                                            :: a, i, j, k
      LOGICAL                                            :: found
      REAL(dp)                                           :: product(3, 3)
      TYPE(csym_type)                                    :: csym

      CALL crys_sym_gen(csym, x, types, h, delta=delta, use_spglib=.FALSE.)
      CALL kpoint_gen(csym, [3, 3, 3], symm=.TRUE., shift=[0.0_dp, 0.0_dp, 0.0_dp])
      IF (csym%nkpoint /= expected) ERROR STOP "Incorrect irreducible mesh size"
      IF (csym%inversion_only .NEQV. fallback) ERROR STOP "Incorrect orbit-guard decision"
      IF (ABS(SUM(csym%wkpoint) - 27.0_dp) > 1.0e-12_dp) ERROR STOP "Incorrect weight sum"
      IF (ANY(csym%kplink <= 0)) ERROR STOP "Incomplete mesh mapping"
      DO k = 1, SIZE(csym%rt, 3)
         DO a = 1, SIZE(types)
            IF (COUNT(csym%f0(:, k) == a) /= 1) ERROR STOP "Nonbijective atom mapping"
         END DO
         IF (ANY(types(csym%f0(:, k)) /= types)) ERROR STOP "Atom types not preserved"
      END DO
      IF (closed) THEN
         DO i = 1, csym%nrtot
            DO j = 1, csym%nrtot
               product = MATMUL(csym%rt(:, :, i), csym%rt(:, :, j))
               found = .FALSE.
               DO k = 1, csym%nrtot
                  IF (MAXVAL(ABS(product - csym%rt(:, :, k))) > 1.0e-12_dp) CYCLE
                  IF (ANY(csym%f0(csym%f0(:, j), i) /= csym%f0(:, k))) CYCLE
                  found = .TRUE.
               END DO
               IF (.NOT. found) ERROR STOP "Operation group not closed"
            END DO
         END DO
      END IF
      CALL release_csym_type(csym)
   END SUBROUTINE check_mesh
END PROGRAM kpsym_unittest

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Accepted and trial complete complex densities in the native Bloch variational space.
!>        No inverse transformation to a finite real-space image set is used for mixing.
! **************************************************************************************************
MODULE qs_kpoint_density_state
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_set_submatrix
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: uspp_method_inverse,&
                                              uspp_method_inverse_1c
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE qs_kpoint_operators,             ONLY: kpoint_slot_s
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: kpoint_density_state_active
   PUBLIC :: kpoint_density_state_initialize, kpoint_density_state_capture
   PUBLIC :: kpoint_density_state_accept, kpoint_density_state_delta, kpoint_density_state_restore
   PUBLIC :: kpoint_density_state_write, kpoint_density_state_read
   PUBLIC :: kpoint_density_state_mulliken
   PUBLIC :: kpoint_density_metric_transport
CONTAINS

! **************************************************************************************************
!> \brief Test explicit complete Bloch density storage; null Gamma pointers are safe.
!> \param kpoints ...
!> \return ...
! **************************************************************************************************
   LOGICAL FUNCTION kpoint_density_state_active(kpoints) RESULT(active)
      TYPE(kpoint_type), POINTER                         :: kpoints

      active = .FALSE.
      IF (ASSOCIATED(kpoints)) active = ALLOCATED(kpoints%native_density)
   END FUNCTION kpoint_density_state_active

! **************************************************************************************************
!> \brief Explicitly attach group-local accepted/trial density storage to a native Bloch state.
!> \param kpoints ...
!> \param representation Identity of the orbital representation stored in restart files.
! **************************************************************************************************
   SUBROUTINE kpoint_density_state_initialize(kpoints, representation)
      TYPE(kpoint_type), INTENT(INOUT)                   :: kpoints
      INTEGER, INTENT(IN)                                :: representation

      INTEGER                                            :: ik, nspin, spin

      CPASSERT(.NOT. ALLOCATED(kpoints%native_density))
      nspin = SIZE(kpoints%kp_env(1)%kpoint_env%pmat, 2)
      ALLOCATE (kpoints%native_density(SIZE(kpoints%kp_env), nspin), &
                kpoints%native_trial(SIZE(kpoints%kp_env), nspin))
      DO spin = 1, nspin
         DO ik = 1, SIZE(kpoints%kp_env)
            CALL cp_cfm_create(kpoints%native_density(ik, spin), &
                               kpoints%kp_env(ik)%kpoint_env%pmat(1, spin)%matrix_struct)
            CALL cp_cfm_create(kpoints%native_trial(ik, spin), &
                               kpoints%kp_env(ik)%kpoint_env%pmat(1, spin)%matrix_struct)
            kpoints%native_density(ik, spin)%local_data = CMPLX(0.0_dp, 0.0_dp, dp)
            kpoints%native_trial(ik, spin)%local_data = CMPLX(0.0_dp, 0.0_dp, dp)
         END DO
      END DO
      kpoints%native_trial_valid = .FALSE.
      ! A representation-specific generalized solve seeds P(k) in the physical metric.
      kpoints%native_density_seeded = .FALSE.
      kpoints%native_density_representation = representation
   END SUBROUTINE kpoint_density_state_initialize

! **************************************************************************************************
!> \brief Write complete complex densities and their physical metrics to a k-point restart.
!>        Each owning group contributes once. Storage is independent of the MPI grouping and
!>        does not project the density onto a finite set of real-space images.
!> \param kpoints ...
!> \param unit Restart unit on the world source, nonpositive on other ranks.
! **************************************************************************************************
   SUBROUTINE kpoint_density_state_write(kpoints, unit)
      TYPE(kpoint_type), INTENT(IN)                      :: kpoints
      INTEGER, INTENT(IN)                                :: unit

      COMPLEX(dp), ALLOCATABLE                           :: matrix(:, :)
      INTEGER                                            :: ik, local, n, nspin, spin
      LOGICAL                                            :: owns

      CPASSERT(kpoints%native_density_seeded)
      CPASSERT(ALLOCATED(kpoints%native_density))
      n = kpoints%native_density(1, 1)%matrix_struct%nrow_global
      nspin = SIZE(kpoints%native_density, 2)
      ALLOCATE (matrix(n, n))
      IF (unit > 0) THEN
         WRITE (unit) kpoints%nkp, nspin, n, kpoints%native_density_representation
         WRITE (unit) kpoints%xkp(1:3, :), kpoints%wkp
      END IF
      DO ik = 1, kpoints%nkp
         local = ik - kpoints%kp_range(1) + 1
         owns = ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2)
         matrix = CMPLX(0.0_dp, 0.0_dp, dp)
         IF (owns) THEN
            CALL cp_cfm_get_submatrix(kpoints%native_operators(kpoint_slot_s)%matrices(local, 1), matrix)
            IF (.NOT. kpoints%para_env_kp%is_source()) matrix = CMPLX(0.0_dp, 0.0_dp, dp)
         END IF
         CALL kpoints%para_env%sum(matrix)
         IF (unit > 0) WRITE (unit) matrix
         DO spin = 1, nspin
            matrix = CMPLX(0.0_dp, 0.0_dp, dp)
            IF (owns) THEN
               CALL cp_cfm_get_submatrix(kpoints%native_density(local, spin), matrix)
               IF (.NOT. kpoints%para_env_kp%is_source()) matrix = CMPLX(0.0_dp, 0.0_dp, dp)
            END IF
            CALL kpoints%para_env%sum(matrix)
            IF (unit > 0) WRITE (unit) matrix
         END DO
      END DO
      DEALLOCATE (matrix)
   END SUBROUTINE kpoint_density_state_write

! **************************************************************************************************
!> \brief Restore complete Bloch densities, preserving occupations in the current physical metric.
!>        A = (M_new^-1 M_old)^1/2 gives P_new = A P_old A^H. This congruence preserves
!>        positivity, occupation eigenvalues, charge and common metric symmetries.
!>        The two inverse methods share reconstructed Gaussian orbital variables and can
!>        use this transport. Direct soft-orbital densities require a separate transformation.
!> \param kpoints ...
!> \param unit Restart unit on the world source, nonpositive on other ranks.
! **************************************************************************************************
   SUBROUTINE kpoint_density_state_read(kpoints, unit)
      TYPE(kpoint_type), INTENT(INOUT)                   :: kpoints
      INTEGER, INTENT(IN)                                :: unit

      COMPLEX(dp), ALLOCATABLE                           :: density(:, :), metric(:, :), &
                                                            previous_metric(:, :), transform(:, :)
      INTEGER                                            :: ik, info, ios, local, n, nao, nk, nspin, &
                                                            representation, spin
      LOGICAL                                            :: compatible_representation, owns
      REAL(dp)                                           :: charge, tolerance
      REAL(dp), ALLOCATABLE                              :: coordinates(:, :), weights(:)

      CPASSERT(ALLOCATED(kpoints%native_density))
      n = kpoints%native_density(1, 1)%matrix_struct%nrow_global
      IF (unit > 0) THEN
         READ (unit, IOSTAT=ios) nk, nspin, nao, representation
         IF (ios /= 0) CPABORT("Incomplete native k-point restart header")
      END IF
      CALL kpoints%para_env%bcast(nk)
      CALL kpoints%para_env%bcast(nspin)
      CALL kpoints%para_env%bcast(nao)
      CALL kpoints%para_env%bcast(representation)
      compatible_representation = representation == kpoints%native_density_representation
      IF (ANY(representation == [uspp_method_inverse, uspp_method_inverse_1c])) THEN
         compatible_representation = &
            ANY(kpoints%native_density_representation == [uspp_method_inverse, uspp_method_inverse_1c])
      END IF
      IF (nk /= kpoints%nkp .OR. nspin /= SIZE(kpoints%native_density, 2) .OR. nao /= n .OR. &
          .NOT. compatible_representation) THEN
         CPABORT("Native k-point restart mesh, spin, basis or representation mismatch")
      END IF
      ALLOCATE (coordinates(3, nk), weights(nk))
      IF (unit > 0) THEN
         READ (unit, IOSTAT=ios) coordinates, weights
         IF (ios /= 0) CPABORT("Incomplete native k-point restart mesh")
      END IF
      CALL kpoints%para_env%bcast(coordinates)
      CALL kpoints%para_env%bcast(weights)
      IF (MAXVAL(ABS(coordinates - kpoints%xkp(1:3, :))) > 1.0E-12_dp .OR. &
          MAXVAL(ABS(weights - kpoints%wkp)) > 1.0E-12_dp) THEN
         CPABORT("Native k-point restart has different coordinates or weights")
      END IF
      DEALLOCATE (coordinates, weights)
      ALLOCATE (density(n, n), metric(n, n), previous_metric(n, n))
      DO ik = 1, nk
         local = ik - kpoints%kp_range(1) + 1
         owns = ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2)
         IF (unit > 0) THEN
            READ (unit, IOSTAT=ios) previous_metric
            IF (ios /= 0) CPABORT("Incomplete native k-point restart metric")
         END IF
         CALL kpoints%para_env%bcast(previous_metric)
         IF (owns) THEN
            CALL cp_cfm_get_submatrix(kpoints%native_operators(kpoint_slot_s)%matrices(local, 1), metric)
            tolerance = 1.0E-10_dp*MAX(1.0_dp, MAXVAL(ABS(previous_metric)))
            IF (MAXVAL(ABS(previous_metric - CONJG(TRANSPOSE(previous_metric)))) > tolerance) THEN
               CPABORT("Non-Hermitian native k-point restart metric")
            END IF
            CALL kpoint_density_metric_transport(previous_metric, metric, transform, info)
            IF (info /= 0) CPABORT("Invalid native k-point restart metric transport")
         END IF
         DO spin = 1, nspin
            IF (unit > 0) THEN
               READ (unit, IOSTAT=ios) density
               IF (ios /= 0) CPABORT("Incomplete native k-point restart density")
            END IF
            CALL kpoints%para_env%bcast(density)
            IF (owns) THEN
               IF (.NOT. ALL(ieee_is_finite(REAL(density, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(density)))) THEN
                  CPABORT("Nonfinite native k-point restart density")
               END IF
               tolerance = 1.0E-10_dp*MAX(1.0_dp, MAXVAL(ABS(density)))
               IF (MAXVAL(ABS(density - CONJG(TRANSPOSE(density)))) > tolerance) THEN
                  CPABORT("Non-Hermitian native k-point restart density")
               END IF
               charge = REAL(SUM(density*CONJG(previous_metric)), dp)
               density(:, :) = MATMUL(transform, MATMUL(density, CONJG(TRANSPOSE(transform))))
               IF (ABS(SUM(density*CONJG(metric)) - charge) > 1.0E-9_dp*MAX(1.0_dp, ABS(charge))) THEN
                  CPABORT("Native k-point restart charge normalization failed")
               END IF
               CALL cp_cfm_set_submatrix(kpoints%native_density(local, spin), density)
            END IF
         END DO
         IF (ALLOCATED(transform)) DEALLOCATE (transform)
      END DO
      DEALLOCATE (density, metric, previous_metric)
      kpoints%native_density_seeded = .TRUE.
      kpoints%native_trial_valid = .FALSE.
      CALL kpoint_density_state_restore(kpoints)
   END SUBROUTINE kpoint_density_state_read

! **************************************************************************************************
!> \brief Positive relative square root A, with A^H M_new A = M_old.
!>        As a function of M_new^-1 M_old it commutes with every common metric symmetry,
!>        including nonunitary AO actions. A Cholesky-to-Cholesky map need not do so.
!> \param previous_metric Hermitian positive previous charge metric M_old
!> \param metric Hermitian positive current charge metric M_new
!> \param transform Principal square root of M_new^-1 M_old
!> \param info Zero on success, nonzero for invalid input or failed factorization
! **************************************************************************************************
   SUBROUTINE kpoint_density_metric_transport(previous_metric, metric, transform, info)
      COMPLEX(dp), INTENT(IN)                            :: previous_metric(:, :), metric(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: transform(:, :)
      INTEGER, INTENT(OUT)                               :: info

      COMPLEX(dp), ALLOCATABLE                           :: relative(:, :), upper(:, :), work(:)
      INTEGER                                            :: n
      REAL(dp), ALLOCATABLE                              :: real_work(:), values(:)

      n = SIZE(metric, 1)
      CPASSERT(n > 0 .AND. ALL(SHAPE(metric) == [n, n]))
      CPASSERT(ALL(SHAPE(previous_metric) == [n, n]))
      info = -1
      IF (.NOT. ALL(ieee_is_finite(REAL(metric, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(metric)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(previous_metric, dp))) .OR. &
          .NOT. ALL(ieee_is_finite(AIMAG(previous_metric)))) RETURN
      IF (MAXVAL(ABS(metric - CONJG(TRANSPOSE(metric)))) > 1.e-10_dp*MAX(1.0_dp, MAXVAL(ABS(metric)))) RETURN
      IF (MAXVAL(ABS(previous_metric - CONJG(TRANSPOSE(previous_metric)))) > &
          1.e-10_dp*MAX(1.0_dp, MAXVAL(ABS(previous_metric)))) RETURN
      ALLOCATE (upper, SOURCE=metric)
      CALL zpotrf('U', n, upper, n, info)
      IF (info /= 0) RETURN
      ! M_new = R^H R. B = R^-H M_old R^-1 is Hermitian positive.
      ALLOCATE (relative, SOURCE=previous_metric)
      CALL ztrsm('L', 'U', 'C', 'N', n, n, CMPLX(1.0_dp, 0.0_dp, dp), upper, n, relative, n)
      CALL ztrsm('R', 'U', 'N', 'N', n, n, CMPLX(1.0_dp, 0.0_dp, dp), upper, n, relative, n)
      relative(:, :) = 0.5_dp*(relative + CONJG(TRANSPOSE(relative)))
      IF (.NOT. ALL(ieee_is_finite(REAL(relative, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(relative)))) THEN
         info = -1
         RETURN
      END IF
      ALLOCATE (values(n), work(2*n), real_work(3*n))
      CALL zheev('V', 'U', n, relative, n, values, work, SIZE(work), real_work, info)
      IF (info /= 0) RETURN
      IF (.NOT. ALL(ieee_is_finite(values)) .OR. MINVAL(values) <= 0.0_dp) THEN
         info = n + 1
         RETURN
      END IF
      ALLOCATE (transform, SOURCE=MATMUL(relative*SPREAD(SQRT(values), 1, n), CONJG(TRANSPOSE(relative))))
      ! A = R^-1 sqrt(B) R. No inverse of either dense metric is formed.
      CALL ztrmm('R', 'U', 'N', 'N', n, n, CMPLX(1.0_dp, 0.0_dp, dp), upper, n, transform, n)
      CALL ztrsm('L', 'U', 'N', 'N', n, n, CMPLX(1.0_dp, 0.0_dp, dp), upper, n, transform, n)
      IF (.NOT. ALL(ieee_is_finite(REAL(transform, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(transform)))) info = -1
   END SUBROUTINE kpoint_density_metric_transport

! **************************************************************************************************
!> \brief Capture the newly diagonalized density while retaining the density that built H.
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE kpoint_density_state_capture(kpoints)
      TYPE(kpoint_type), INTENT(INOUT)                   :: kpoints

      INTEGER                                            :: ik, spin
      TYPE(kpoint_env_type), POINTER                     :: kp

      CPASSERT(ALLOCATED(kpoints%native_density))
      DO spin = 1, SIZE(kpoints%native_trial, 2)
         DO ik = 1, SIZE(kpoints%native_trial, 1)
            kp => kpoints%kp_env(ik)%kpoint_env
            kpoints%native_trial(ik, spin)%local_data(:, :) = CMPLX(kp%pmat(1, spin)%local_data, 0.0_dp, dp)
            IF (.NOT. kpoints%use_real_wfn) THEN
               kpoints%native_trial(ik, spin)%local_data(:, :) = &
                  kpoints%native_trial(ik, spin)%local_data + CMPLX(0.0_dp, kp%pmat(2, spin)%local_data, dp)
            END IF
         END DO
      END DO
      kpoints%native_trial_valid = .TRUE.
   END SUBROUTINE kpoint_density_state_capture

! **************************************************************************************************
!> \brief Maximum complex AO density change, reduced over the world once without k weights.
!> \param kpoints ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE kpoint_density_state_delta(kpoints, delta)
      TYPE(kpoint_type), INTENT(IN)                      :: kpoints
      REAL(dp), INTENT(OUT)                              :: delta

      INTEGER                                            :: ik, spin

      CPASSERT(kpoints%native_trial_valid)
      delta = 0.0_dp
      DO spin = 1, SIZE(kpoints%native_trial, 2)
         DO ik = 1, SIZE(kpoints%native_trial, 1)
            delta = MAX(delta, MAXVAL(ABS(kpoints%native_trial(ik, spin)%local_data - &
                                          kpoints%native_density(ik, spin)%local_data)))
         END DO
      END DO
      CALL kpoints%para_env%max(delta)
   END SUBROUTINE kpoint_density_state_delta

! **************************************************************************************************
!> \brief Accept a convex direct density step. The raw trial remains available for finalization.
!> \param kpoints ...
!> \param alpha ...
! **************************************************************************************************
   SUBROUTINE kpoint_density_state_accept(kpoints, alpha)
      TYPE(kpoint_type), INTENT(INOUT)                   :: kpoints
      REAL(dp), INTENT(IN)                               :: alpha

      INTEGER                                            :: ik, spin

      CPASSERT(kpoints%native_trial_valid)
      CPASSERT(alpha >= 0.0_dp .AND. alpha <= 1.0_dp)
      DO spin = 1, SIZE(kpoints%native_trial, 2)
         DO ik = 1, SIZE(kpoints%native_trial, 1)
            IF (alpha == 1.0_dp) THEN
               kpoints%native_density(ik, spin)%local_data(:, :) = kpoints%native_trial(ik, spin)%local_data
            ELSE
               kpoints%native_density(ik, spin)%local_data(:, :) = &
                  (1.0_dp - alpha)*kpoints%native_density(ik, spin)%local_data + &
                  alpha*kpoints%native_trial(ik, spin)%local_data
            END IF
         END DO
      END DO
      CALL kpoint_density_state_restore(kpoints)
   END SUBROUTINE kpoint_density_state_accept

! **************************************************************************************************
!> \brief Restore the accepted P(k) used by the next native commutator/DIIS evaluation.
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE kpoint_density_state_restore(kpoints)
      TYPE(kpoint_type), INTENT(INOUT)                   :: kpoints

      INTEGER                                            :: ik, spin
      TYPE(kpoint_env_type), POINTER                     :: kp

      CPASSERT(ALLOCATED(kpoints%native_density))
      DO spin = 1, SIZE(kpoints%native_density, 2)
         DO ik = 1, SIZE(kpoints%native_density, 1)
            kp => kpoints%kp_env(ik)%kpoint_env
            kp%pmat(1, spin)%local_data(:, :) = REAL(kpoints%native_density(ik, spin)%local_data, dp)
            IF (.NOT. kpoints%use_real_wfn) THEN
               kp%pmat(2, spin)%local_data(:, :) = AIMAG(kpoints%native_density(ik, spin)%local_data)
            END IF
         END DO
      END DO
   END SUBROUTINE kpoint_density_state_restore

! **************************************************************************************************
!> \brief Mulliken populations from the accepted complex density and its registered charge metric.
!>        Each group source contributes its actual k-point weights once. The caller sums world ranks.
!> \param kpoints ...
!> \param populations Partial orbital populations, overwritten.
! **************************************************************************************************
   SUBROUTINE kpoint_density_state_mulliken(kpoints, populations)
      TYPE(kpoint_type), INTENT(IN)                      :: kpoints
      REAL(dp), INTENT(OUT)                              :: populations(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: density(:, :), metric(:, :)
      INTEGER                                            :: ik, local, n, spin

      CPASSERT(ALLOCATED(kpoints%native_density))
      CPASSERT(ALLOCATED(kpoints%native_operators(kpoint_slot_s)%matrices))
      n = SIZE(populations, 1)
      CPASSERT(SIZE(populations, 2) == SIZE(kpoints%native_density, 2))
      ALLOCATE (density(n, n), metric(n, n))
      populations = 0.0_dp
      DO local = 1, SIZE(kpoints%native_density, 1)
         ik = kpoints%kp_range(1) + local - 1
         CALL cp_cfm_get_submatrix(kpoints%native_operators(kpoint_slot_s)%matrices(local, 1), metric)
         DO spin = 1, SIZE(populations, 2)
            CALL cp_cfm_get_submatrix(kpoints%native_density(local, spin), density)
            IF (kpoints%para_env_kp%is_source()) THEN
               populations(:, spin) = populations(:, spin) + &
                                      kpoints%wkp(ik)*SUM(REAL(density*CONJG(metric), dp), DIM=2)
            END IF
         END DO
      END DO
   END SUBROUTINE kpoint_density_state_mulliken
END MODULE qs_kpoint_density_state

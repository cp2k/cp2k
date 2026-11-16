!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Metric-covariant spectral flattening with existing distributed matrix iterations.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE spectral_localizer_flatten
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_create,&
                                              dbcsr_get_info,&
                                              dbcsr_multiply,&
                                              dbcsr_release,&
                                              dbcsr_type
   USE cp_dbcsr_contrib,                ONLY: dbcsr_add_on_diag,&
                                              dbcsr_frobenius_norm
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE iterate_matrix,                  ONLY: invert_Hotelling,&
                                              matrix_sign_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE spectral_localizer,              ONLY: localizer_matrix_evaluate
   USE spectral_localizer_dbcsr,        ONLY: localizer_dbcsr_check,&
                                              localizer_dbcsr_dense,&
                                              localizer_dbcsr_from_real,&
                                              localizer_dbcsr_gap_pencil,&
                                              localizer_dbcsr_realify,&
                                              localizer_dbcsr_release,&
                                              localizer_matrix_type

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: localizer_flatten_dbcsr, localizer_flatten_checked

CONTAINS

! **************************************************************************************************
!> \brief Establish a positive AO metric and electronic gap before iterating the matrix sign.
!> \param h Covariant physical Hamiltonian, unchanged
!> \param s Covariant physical metric, unchanged
!> \param energy Reference energy
!> \param scale Positive flattened-band energy
!> \param tolerance Matrix-sign residual and Hermiticity tolerance
!> \param metric_tolerance Minimum overlap eigenvalue ratio
!> \param gap_tolerance Absolute minimum original electronic gap
!> \param para_env Communicator
!> \param flattened Newly initialized shifted Hamiltonian only on success
!> \param gap Original electronic gap lower and upper bounds, not the flattened localizer gap
!> \param metric_ratio Overlap eigenvalue ratio
!> \param residuals Four matrix-sign residuals
!> \param status 0 success, 1 unresolved input gap, -1 input, -2 metric, -4 solver,
!>        -6 matrix-sign iteration or structure failure
! **************************************************************************************************
   SUBROUTINE localizer_flatten_checked(h, s, energy, scale, tolerance, metric_tolerance, gap_tolerance, &
                                        para_env, flattened, gap, metric_ratio, residuals, status)
      TYPE(localizer_matrix_type), INTENT(INOUT)         :: h, s
      REAL(dp), INTENT(IN)                               :: energy, scale, tolerance, &
                                                            metric_tolerance, gap_tolerance
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(localizer_matrix_type), INTENT(OUT)           :: flattened
      REAL(dp), INTENT(OUT)                              :: gap(2), metric_ratio, residuals(4)
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp), ALLOCATABLE                           :: dense(:, :), dense_metric(:, :)
      INTEGER                                            :: invariant, n, ns
      LOGICAL                                            :: valid
      TYPE(localizer_matrix_type)                        :: metric, pencil

      status = -1
      gap = 0.0_dp
      metric_ratio = 0.0_dp
      residuals = HUGE(1.0_dp)
      IF (.NOT. ALL(ieee_is_finite([energy, scale, tolerance, metric_tolerance, gap_tolerance]))) RETURN
      IF (scale <= 0.0_dp .OR. tolerance <= 0.0_dp .OR. tolerance >= 1.0_dp) RETURN
      IF (metric_tolerance <= 0.0_dp .OR. metric_tolerance >= 1.0_dp .OR. gap_tolerance < 0.0_dp) RETURN
      CALL dbcsr_get_info(h%real_part, nfullrows_total=n)
      CALL dbcsr_get_info(s%real_part, nfullrows_total=ns)
      IF (n < 1 .OR. n /= ns) RETURN
      CALL localizer_dbcsr_check(h, para_env, tolerance, valid)
      IF (.NOT. valid) RETURN
      CALL localizer_dbcsr_check(s, para_env, tolerance, valid)
      IF (.NOT. valid) RETURN
      CALL localizer_dbcsr_gap_pencil(h, s, energy, pencil, metric)
      CALL localizer_dbcsr_dense(pencil, para_env, dense)
      CALL localizer_dbcsr_dense(metric, para_env, dense_metric)
      IF (para_env%is_source()) THEN
         CALL localizer_matrix_evaluate(dense, dense_metric(:n, :n), .FALSE., tolerance, metric_tolerance, &
                                        gap_tolerance, invariant, gap(1), metric_ratio, status)
         IF (status > 1) status = -4
         gap(2) = gap(1)
      END IF
      CALL localizer_dbcsr_release(pencil)
      CALL localizer_dbcsr_release(metric)
      CALL para_env%bcast(status, para_env%source)
      CALL para_env%bcast(gap, para_env%source)
      CALL para_env%bcast(metric_ratio, para_env%source)
      IF (status /= 0) RETURN
      CALL para_env%bcast(invariant, para_env%source)
      IF (invariant /= 0) THEN
         status = -4
         RETURN
      END IF
      CALL localizer_flatten_dbcsr(h, s, energy, scale, tolerance, para_env, flattened, residuals, status)
      IF (status /= 0) status = -6
   END SUBROUTINE localizer_flatten_checked

! **************************************************************************************************
!> \brief Form scale*S*sign(S^-1*(H-energy*S)) without dense AO matrices or S^-1/2.
!> \param h Covariant Hermitian Hamiltonian; unchanged
!> \param s Positive definite AO metric; must be validated by the caller
!> \param energy Reference energy, separated from the generalized spectrum at input accuracy
!> \param scale Positive energy scale of the flattened bands
!> \param tolerance Residual tolerance, not a certificate of the original Hamiltonian gap
!> \param para_env DBCSR communicator
!> \param flattened Newly initialized covariant operator only on success
!> \param residuals Inverse, sign-square, complex-structure and Hermiticity residuals
!> \param status 0 success, -1 invalid input, 1 unconverged sign/inverse, 2 broken structure
!> \note Caller must establish a positive metric and a resolved electronic gap independently.
!>       Matrix-sign convergence alone cannot certify the occupied subspace near a gap closing.
! **************************************************************************************************
   SUBROUTINE localizer_flatten_dbcsr(h, s, energy, scale, tolerance, para_env, flattened, residuals, status)
      TYPE(localizer_matrix_type), INTENT(INOUT)         :: h, s
      REAL(dp), INTENT(IN)                               :: energy, scale, tolerance
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(localizer_matrix_type), INTENT(OUT)           :: flattened
      REAL(dp), INTENT(OUT)                              :: residuals(4)
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: n, ns
      LOGICAL                                            :: valid
      REAL(dp)                                           :: filter, norm
      TYPE(dbcsr_type)                                   :: a, b, inverse, sign_matrix, trial, work

      status = -1
      residuals = HUGE(1.0_dp)
      IF (.NOT. ALL(ieee_is_finite([energy, scale, tolerance]))) RETURN
      IF (scale <= 0.0_dp .OR. tolerance <= 0.0_dp .OR. tolerance >= 1.0_dp) RETURN
      CALL dbcsr_get_info(h%real_part, nfullrows_total=n)
      CALL dbcsr_get_info(s%real_part, nfullrows_total=ns)
      IF (n < 1 .OR. n /= ns) RETURN
      CALL localizer_dbcsr_check(h, para_env, tolerance, valid)
      IF (.NOT. valid) RETURN
      CALL localizer_dbcsr_check(s, para_env, tolerance, valid)
      IF (.NOT. valid) RETURN
      CALL localizer_dbcsr_realify(h, a)
      CALL localizer_dbcsr_realify(s, b)
      CALL dbcsr_add(a, b, 1.0_dp, -energy)
      CALL dbcsr_create(inverse, template=b)
      CALL dbcsr_create(sign_matrix, template=b)
      CALL dbcsr_create(trial, template=b)
      CALL dbcsr_create(work, template=b)
      filter = MIN(1.e-16_dp, 0.01_dp*tolerance**2)
      DO
         norm = dbcsr_frobenius_norm(a)
         status = 1
         IF (.NOT. ieee_is_finite(norm) .OR. norm == 0.0_dp) EXIT
         CALL invert_Hotelling(inverse, b, 0.1_dp*tolerance, filter_eps=filter, &
                               accelerator_order=0, silent=.TRUE.)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, b, inverse, 0.0_dp, work)
         CALL dbcsr_add_on_diag(work, -1.0_dp)
         residuals(1) = dbcsr_frobenius_norm(work)
         IF (.NOT. ieee_is_finite(residuals(1)) .OR. residuals(1) > tolerance) EXIT
         CALL dbcsr_multiply('N', 'N', 1.0_dp, inverse, a, 0.0_dp, trial)
         CALL matrix_sign_Newton_Schulz(sign_matrix, trial, filter, sign_order=2, iounit=-1)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, sign_matrix, sign_matrix, 0.0_dp, work)
         CALL dbcsr_add_on_diag(work, -1.0_dp)
         residuals(2) = dbcsr_frobenius_norm(work)
         IF (.NOT. ieee_is_finite(residuals(2)) .OR. residuals(2) > tolerance) EXIT
         CALL dbcsr_multiply('N', 'N', scale, b, sign_matrix, 0.0_dp, trial)
         CALL localizer_dbcsr_from_real(trial, h, flattened)
         CALL dbcsr_release(work)
         CALL localizer_dbcsr_realify(flattened, work)
         norm = MAX(dbcsr_frobenius_norm(trial), TINY(1.0_dp))
         CALL dbcsr_add(work, trial, 1.0_dp, -1.0_dp)
         residuals(3) = dbcsr_frobenius_norm(work)/norm
         CALL localizer_dbcsr_check(flattened, para_env, tolerance, valid, residuals(4))
         status = 2
         IF (.NOT. valid .OR. .NOT. ALL(ieee_is_finite(residuals)) .OR. MAXVAL(residuals) > tolerance) THEN
            CALL localizer_dbcsr_release(flattened)
            EXIT
         END IF
         status = 0
         EXIT
      END DO
      CALL dbcsr_release(a)
      CALL dbcsr_release(b)
      CALL dbcsr_release(inverse)
      CALL dbcsr_release(sign_matrix)
      CALL dbcsr_release(trial)
      CALL dbcsr_release(work)
   END SUBROUTINE localizer_flatten_dbcsr

END MODULE spectral_localizer_flatten

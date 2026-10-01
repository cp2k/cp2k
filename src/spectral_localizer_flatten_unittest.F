!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Complex metric-covariant matrix sign against exact spectral constructions.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
PROGRAM spectral_localizer_flatten_unittest
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_finalize_lib, dbcsr_init_lib, dbcsr_put_block, dbcsr_type_no_symmetry
   USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                              cp_logger_create,&
                                              cp_logger_release,&
                                              cp_logger_type,&
                                              cp_rm_default_logger
   USE fft_tools,                       ONLY: finalize_fft,&
                                              init_fft
   USE ieee_arithmetic,                 ONLY: ieee_quiet_nan,&
                                              ieee_value
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: default_output_unit
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE spectral_localizer_bloch,        ONLY: localizer_flatten_bloch
   USE spectral_localizer_dbcsr,        ONLY: localizer_dbcsr_dense,&
                                              localizer_dbcsr_release,&
                                              localizer_matrix_type
   USE spectral_localizer_flatten,      ONLY: localizer_flatten_checked,&
                                              localizer_flatten_dbcsr
#include "./base/base_uses.f90"

   IMPLICIT NONE
   INTEGER, PARAMETER :: n = 12, nb = 6, bs = 2
   TYPE(mp_para_env_type), POINTER :: para_env
   TYPE(cp_logger_type), POINTER :: logger
   TYPE(dbcsr_distribution_type) :: distribution
   TYPE(localizer_matrix_type) :: h, s, flat
   INTEGER, POINTER :: process_grid(:, :), row_dist(:), col_dist(:), sizes(:)
   COMPLEX(dp) :: transform(n, n), u(n, n), rotation(n, n), scaled(n, n), hdense(n, n), sdense(n, n), exact(n, n)
   COMPLEX(dp), ALLOCATABLE :: result(:, :)
   REAL(dp) :: energies(n), residuals(4), difference, maximum_error, energy, scale, gap(2), metric_ratio, exact_gap
   INTEGER :: i, j, mode, nr, nc, io, status
   CHARACTER(LEN=32) :: argument

   ALLOCATE (para_env)
   CALL mp_world_init(para_env)
   io = -1
   IF (para_env%is_source()) io = default_output_unit
   CALL cp_logger_create(logger, para_env=para_env, default_global_unit_nr=io, close_global_unit_on_dealloc=.FALSE.)
   CALL cp_add_default_logger(logger)
   CALL dbcsr_init_lib(para_env%get_handle(), io)
   CALL init_fft('FFTSG', alltoall=.FALSE., pool_limit=4, wisdom_file='', plan_style=1)
   nr = para_env%num_pe
   nc = 1
   DO i = 1, command_argument_count()
      CALL get_command_argument(i, argument)
      IF (TRIM(argument) == '--column-grid') THEN
         nr = 1
         nc = para_env%num_pe
      END IF
   END DO
   ALLOCATE (process_grid(0:nr - 1, 0:nc - 1), row_dist(nb), col_dist(nb), sizes(nb))
   process_grid = RESHAPE([(i, i=0, para_env%num_pe - 1)], [nr, nc])
   row_dist = [(MODULO(i - 1, nr), i=1, nb)]
   col_dist = [(MODULO(i - 1, nc), i=1, nb)]
   sizes = bs
   CALL dbcsr_distribution_new(distribution, group=para_env%get_handle(), pgrid=process_grid, &
                               row_dist=row_dist, col_dist=col_dist)
   energy = 0.002_dp
   scale = 1.7_dp
   maximum_error = 0.0_dp
   DO mode = 0, 5
      u = 0.0_dp
      transform = 0.0_dp
      DO i = 1, n
         energies(i) = 0.01_dp*(i - n/2) - 0.005_dp
         transform(i, i) = 1.0_dp
         DO j = 1, n
            u(i, j) = EXP(CMPLX(0.0_dp, twopi*REAL((i - 1)*(j - 1), dp)/n, dp))/SQRT(REAL(n, dp))
         END DO
         IF (mode == 0) CYCLE
         transform(i, i) = 0.5_dp + 0.1_dp*i
         IF (i < n) transform(i, i + 1) = CMPLX(0.13_dp, -0.08_dp, dp)
      END DO
      IF (mode == 2) transform(n, n) = 0.1_dp
      IF (mode == 4) energies = energies + 0.5_dp
      IF (mode == 5) energies = energies - 0.5_dp
      IF (mode == 3) THEN
         u = 0.0_dp
         rotation = 0.0_dp
         transform(:n/2, n/2 + 1:) = 0.0_dp
         transform(n/2 + 1:, :n/2) = 0.0_dp
         transform(n/2 + 1:, n/2 + 1:) = CONJG(transform(:n/2, :n/2))
         DO i = 1, n/2
            energies(i) = 0.01_dp*(i - n/4) - 0.005_dp
            energies(n/2 + i) = energies(i)
            rotation(i, i) = COS(0.37_dp)
            rotation(n/2 + i, n/2 + i) = COS(0.37_dp)
            rotation(i, n/2 + i) = SIN(0.37_dp)
            rotation(n/2 + i, i) = -SIN(0.37_dp)
            DO j = 1, n/2
               u(i, j) = EXP(CMPLX(0.0_dp, twopi*REAL((i - 1)*(j - 1), dp)/(n/2), dp))/SQRT(REAL(n/2, dp))
            END DO
         END DO
         u(n/2 + 1:, n/2 + 1:) = CONJG(u(:n/2, :n/2))
         u = MATMUL(rotation, u)
      END IF
      scaled = u*SPREAD(energies, 1, n)
      hdense = MATMUL(CONJG(TRANSPOSE(transform)), MATMUL(MATMUL(scaled, CONJG(TRANSPOSE(u))), transform))
      sdense = MATMUL(CONJG(TRANSPOSE(transform)), transform)
      scaled = u*SPREAD(SIGN(scale, energies - energy), 1, n)
      exact = MATMUL(CONJG(TRANSPOSE(transform)), MATMUL(MATMUL(scaled, CONJG(TRANSPOSE(u))), transform))
      CALL make_operator(hdense, h)
      CALL make_operator(sdense, s)
      CALL localizer_flatten_dbcsr(h, s, energy, scale, 1.e-9_dp, para_env, flat, residuals, status)
      IF (status /= 0) THEN
         IF (io > 0) WRITE (io, *) 'Flatten failure: ', mode, status, residuals
         CPABORT('Distributed spectral flattening failed')
      END IF
      CALL localizer_dbcsr_dense(flat, para_env, RESULT)
      difference = MAXVAL(ABS(RESULT - exact))
      maximum_error = MAX(maximum_error, difference)
      IF (difference > 1.e-9_dp) THEN
         CPABORT('Covariant sign differs from exact eigensystem')
      END IF
      IF (mode == 3) THEN
         IF (MAXVAL(ABS(RESULT(:n/2, :n/2) - CONJG(RESULT(n/2 + 1:, n/2 + 1:)))) > 1.e-9_dp) THEN
            CPABORT('Flattening broke physical TR diagonal blocks')
         END IF
         IF (MAXVAL(ABS(RESULT(:n/2, n/2 + 1:) + CONJG(RESULT(n/2 + 1:, :n/2)))) > 1.e-9_dp) THEN
            CPABORT('Flattening broke physical TR spin mixing')
         END IF
      END IF
      DEALLOCATE (RESULT)
      CALL localizer_dbcsr_dense(h, para_env, RESULT)
      IF (MAXVAL(ABS(RESULT - hdense)) > 1.e-14_dp) THEN
         CPABORT('Flattening modified input Hamiltonian')
      END IF
      DEALLOCATE (RESULT)
      CALL localizer_dbcsr_dense(s, para_env, RESULT)
      IF (MAXVAL(ABS(RESULT - sdense)) > 1.e-14_dp) THEN
         CPABORT('Flattening modified input metric')
      END IF
      DEALLOCATE (RESULT)
      CALL localizer_dbcsr_release(flat)
      CALL localizer_flatten_checked(h, s, energy, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                     para_env, flat, gap, metric_ratio, residuals, status)
      IF (status /= 0) CPABORT('Checked flattening failed')
      exact_gap = MINVAL(ABS(energies - energy))
      IF (gap(1) > exact_gap + 1.e-10_dp .OR. gap(2) < exact_gap - 1.e-10_dp) THEN
         CPABORT('Input electronic gap differs from reference')
      END IF
      IF (metric_ratio <= 1.e-10_dp) CPABORT('Invalid positive metric bound')
      CALL localizer_dbcsr_dense(flat, para_env, RESULT)
      difference = MAXVAL(ABS(RESULT - exact))
      maximum_error = MAX(maximum_error, difference)
      IF (difference > 1.e-9_dp) THEN
         CPABORT('Checked flattening changed the exact result')
      END IF
      DEALLOCATE (RESULT)
      CALL localizer_dbcsr_release(flat)
      CALL localizer_dbcsr_release(h)
      CALL localizer_dbcsr_release(s)
   END DO

   hdense = 0.0_dp
   sdense = 0.0_dp
   DO i = 1, n
      sdense(i, i) = 1.0_dp
      hdense(i, i) = 1.0_dp
   END DO
   hdense(1, 1) = 0.0_dp
   CALL make_operator(hdense, h)
   CALL make_operator(sdense, s)
   CALL localizer_flatten_dbcsr(h, s, 0.0_dp, scale, 1.e-9_dp, para_env, flat, residuals, status)
   IF (status /= 1) CPABORT('Singular spectral sign was accepted')
   CALL localizer_flatten_checked(h, s, 0.0_dp, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                  para_env, flat, gap, metric_ratio, residuals, status)
   IF (status /= 1) CPABORT('Checked flattening accepted a zero spectral mode')
   CALL localizer_flatten_dbcsr(h, s, energy, -scale, 1.e-9_dp, para_env, flat, residuals, status)
   IF (status /= -1) CPABORT('Negative flattening scale was accepted')
   CALL localizer_flatten_dbcsr(h, s, ieee_value(energy, ieee_quiet_nan), scale, 1.e-9_dp, &
                                para_env, flat, residuals, status)
   IF (status /= -1) CPABORT('Nonfinite energy was accepted')
   CALL localizer_dbcsr_release(h)
   hdense = 0.0_dp
   CALL make_operator(hdense, h)
   CALL localizer_flatten_dbcsr(h, s, 0.0_dp, scale, 1.e-9_dp, para_env, flat, residuals, status)
   IF (status /= 1) CPABORT('Zero Hamiltonian sign was accepted')
   CALL localizer_dbcsr_release(h)
   hdense(1, 2) = 1.0_dp
   CALL make_operator(hdense, h)
   CALL localizer_flatten_dbcsr(h, s, energy, scale, 1.e-9_dp, para_env, flat, residuals, status)
   IF (status /= -1) CPABORT('Non-Hermitian Hamiltonian was accepted')
   CALL localizer_dbcsr_release(h)
   CALL localizer_dbcsr_release(s)
   hdense = sdense
   hdense(1, 1) = 1.e-11_dp
   CALL make_operator(hdense, h)
   CALL make_operator(sdense, s)
   CALL localizer_flatten_checked(h, s, 0.0_dp, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                  para_env, flat, gap, metric_ratio, residuals, status)
   IF (status /= 1) CPABORT('Unresolved near-zero spectral mode was accepted')
   CALL localizer_dbcsr_release(s)
   sdense(1, 1) = -1.0_dp
   CALL make_operator(sdense, s)
   CALL localizer_flatten_checked(h, s, 0.0_dp, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                  para_env, flat, gap, metric_ratio, residuals, status)
   IF (status /= -2) CPABORT('Indefinite metric was accepted')
   CALL localizer_dbcsr_release(s)
   sdense(1, 1) = 1.e-12_dp
   CALL make_operator(sdense, s)
   CALL localizer_flatten_checked(h, s, 0.0_dp, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                  para_env, flat, gap, metric_ratio, residuals, status)
   IF (status /= -2) CPABORT('Ill-conditioned metric was accepted')
   CALL localizer_dbcsr_release(h)
   CALL localizer_dbcsr_release(s)
   CALL check_bloch()
   IF (io > 0) WRITE (io, '(A,ES14.5)') 'Spectral flattening exact-matrix maximum error: ', maximum_error
   CALL dbcsr_distribution_release(distribution)
   DEALLOCATE (process_grid, row_dist, col_dist, sizes)
   CALL dbcsr_finalize_lib()
   CALL finalize_fft(para_env, '')
   CALL cp_rm_default_logger()
   CALL cp_logger_release(logger)
   DEALLOCATE (para_env)
   CALL mp_world_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Complex nonorthogonal full-band Fourier references in spin/cell/atom order.
! **************************************************************************************************
   SUBROUTINE check_bloch()
      INTEGER                                            :: band, grid(3), k, kr(3), m, ncell, &
                                                            nscalar, nsmall, nspin, p, pnext, &
                                                            r(3), rnext, row, t, tnext
      REAL(dp)                                           :: bloch_error, bloch_gap, norm_error, &
                                                            phase, ref_gap

      bloch_error = 0.0_dp
      DO m = 1, 7
         sizes = bs
         grid = [3, 1, 1]
         nspin = 2
         IF (m == 2) grid = [1, 3, 1]
         IF (m == 3) grid = [1, 1, 3]
         IF (m == 4) THEN
            grid = [2, 3, 1]
            nspin = 1
         END IF
         IF (m == 5) THEN
            nspin = 1
            sizes = [1, 3, 1, 3, 1, 3]
         END IF
         IF (m == 6) grid = 1
         ncell = PRODUCT(grid)
         nscalar = n/(nspin*ncell)
         nsmall = nscalar*nspin
         u = 0.0_dp
         transform = 0.0_dp
         DO row = 1, n
            t = MODULO((row - 1)/nscalar, ncell)
            p = ((row - 1)/(ncell*nscalar))*nscalar + MODULO(row - 1, nscalar)
            r = [MODULO(t, grid(1)), MODULO(t/grid(1), grid(2)), t/(grid(1)*grid(2))]
            transform(row, row) = 1.0_dp + 0.03_dp*p
            pnext = MODULO(p + 1, nsmall)
            rnext = (pnext/nscalar)*ncell*nscalar + t*nscalar + MODULO(pnext, nscalar) + 1
            transform(row, rnext) = transform(row, rnext) + CMPLX(0.07_dp, -0.04_dp, dp)
            tnext = MODULO(t + 1, ncell)
            rnext = (p/nscalar)*ncell*nscalar + tnext*nscalar + MODULO(p, nscalar) + 1
            ! A 1D cyclic displacement is a lattice translation for these first-axis-fast grids.
            IF (m /= 4) transform(row, rnext) = transform(row, rnext) + CMPLX(0.05_dp, 0.02_dp, dp)
            DO k = 0, ncell - 1
               kr = [MODULO(k, grid(1)), MODULO(k/grid(1), grid(2)), k/(grid(1)*grid(2))]
               DO band = 0, nsmall - 1
                  phase = twopi*(SUM(REAL(kr*r, dp)/REAL(grid, dp)) + REAL(p*band, dp)/REAL(nsmall, dp))
                  u(row, k*nsmall + band + 1) = EXP(CMPLX(0.0_dp, phase, dp))/SQRT(REAL(n, dp))
                  energies(k*nsmall + band + 1) = 0.2_dp + 0.03_dp*band + 0.01_dp*COS(REAL(k, dp))
                  IF (m == 7) energies(k*nsmall + band + 1) = 0.25_dp + 0.01_dp*COS(REAL(k, dp))
                  IF (band < nsmall/2) energies(k*nsmall + band + 1) = -energies(k*nsmall + band + 1)
               END DO
            END DO
         END DO
         scaled = u*SPREAD(energies, 1, n)
         hdense = MATMUL(CONJG(TRANSPOSE(transform)), MATMUL(MATMUL(scaled, CONJG(TRANSPOSE(u))), transform))
         sdense = MATMUL(CONJG(TRANSPOSE(transform)), transform)
         scaled = u*SPREAD(SIGN(scale, energies - energy), 1, n)
         exact = MATMUL(CONJG(TRANSPOSE(transform)), MATMUL(MATMUL(scaled, CONJG(TRANSPOSE(u))), transform))
         CALL make_operator(hdense, h)
         CALL make_operator(sdense, s)
         CALL localizer_flatten_bloch(h, s, grid, nspin, energy, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                      para_env, flat, bloch_gap, metric_ratio, residuals, status)
         IF (status /= 0) CPABORT('Bloch flattening failed')
         CALL localizer_dbcsr_dense(flat, para_env, RESULT)
         norm_error = MAXVAL(ABS(RESULT - exact))
         bloch_error = MAX(bloch_error, norm_error)
         maximum_error = MAX(maximum_error, norm_error)
         ref_gap = MINVAL(ABS(energies - energy))
         IF (norm_error > 1.e-11_dp .OR. ABS(bloch_gap - ref_gap) > 1.e-12_dp) THEN
            CPABORT('Bloch flattening differs from exact full-torus reference')
         END IF
         CALL localizer_dbcsr_release(flat)
         CALL localizer_flatten_checked(h, s, energy, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                        para_env, flat, gap, metric_ratio, residuals, status)
         IF (status /= 0) CPABORT('Full-torus reference flattening failed')
         CALL localizer_dbcsr_dense(flat, para_env, RESULT)
         IF (MAXVAL(ABS(RESULT - exact)) > 1.e-10_dp) THEN
            CPABORT('Full-torus/Bloch sign mismatch')
         END IF
         CALL localizer_dbcsr_release(flat)
         CALL localizer_dbcsr_release(h)
         IF (m /= 6) THEN
            hdense(1, 1) = hdense(1, 1) + 0.01_dp
            CALL make_operator(hdense, h)
            CALL localizer_flatten_bloch(h, s, grid, nspin, energy, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                         para_env, flat, bloch_gap, metric_ratio, residuals, status)
            IF (status /= -7) CPABORT('Broken translations accepted')
            CALL localizer_dbcsr_release(h)
         END IF
         CALL make_operator(energy*sdense, h)
         CALL localizer_flatten_bloch(h, s, grid, nspin, energy, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                      para_env, flat, bloch_gap, metric_ratio, residuals, status)
         IF (status /= 1) CPABORT('Closed Bloch gap accepted')
         CALL localizer_dbcsr_release(s)
         CALL make_operator(-sdense, s)
         CALL localizer_flatten_bloch(h, s, grid, nspin, energy, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                      para_env, flat, bloch_gap, metric_ratio, residuals, status)
         IF (status /= -2) CPABORT('Bloch flattening accepted an indefinite metric')
         CALL localizer_dbcsr_release(h)
         CALL localizer_dbcsr_release(s)
      END DO
      sizes = bs
      ! Each Bloch metric is well conditioned, but their union has an unacceptable ratio.
      hdense = 0.0_dp
      sdense = 0.0_dp
      DO row = 1, n
         hdense(row, row) = 1.0_dp
         sdense(row, row) = 0.5_dp*(1.0_dp + 1.e-12_dp)
         rnext = 1 + MODULO(row - 1 + n/2, n)
         sdense(row, rnext) = 0.5_dp*(1.0_dp - 1.e-12_dp)
      END DO
      CALL make_operator(hdense, h)
      CALL make_operator(sdense, s)
      CALL localizer_flatten_bloch(h, s, [2, 1, 1], 1, 0.0_dp, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                   para_env, flat, bloch_gap, metric_ratio, residuals, status)
      IF (status /= -2) CPABORT('Global Bloch metric bound failed')
      CALL localizer_flatten_bloch(h, s, [0, 1, 1], 1, 0.0_dp, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                   para_env, flat, bloch_gap, metric_ratio, residuals, status)
      IF (status /= -1) CPABORT('Invalid Bloch grid accepted')
      CALL localizer_dbcsr_release(h)
      CALL localizer_dbcsr_release(s)
      IF (io > 0) WRITE (io, '(A,ES14.5)') 'Seven Bloch flattening pencils; maximum error: ', bloch_error
   END SUBROUTINE check_bloch

! **************************************************************************************************
!> \brief Insert only locally owned dense-reference blocks.
!> \param dense Reference matrix
!> \param matrix New complex DBCSR operator
! **************************************************************************************************
   SUBROUTINE make_operator(dense, matrix)
      COMPLEX(dp), INTENT(IN)                            :: dense(:, :)
      TYPE(localizer_matrix_type), INTENT(OUT)           :: matrix

      INTEGER                                            :: col, first, last, left, right, row

      CALL dbcsr_create(matrix%real_part, 'FLATTEN REAL', distribution, dbcsr_type_no_symmetry, sizes, sizes)
      CALL dbcsr_create(matrix%imaginary_part, template=matrix%real_part)
      DO row = 1, nb
         first = SUM(sizes(:row - 1)) + 1
         last = SUM(sizes(:row))
         DO col = 1, nb
            IF (process_grid(row_dist(row), col_dist(col)) /= para_env%mepos) CYCLE
            left = SUM(sizes(:col - 1)) + 1
            right = SUM(sizes(:col))
            CALL dbcsr_put_block(matrix%real_part, row, col, REAL(dense(first:last, left:right), dp))
            CALL dbcsr_put_block(matrix%imaginary_part, row, col, AIMAG(dense(first:last, left:right)))
         END DO
      END DO
      CALL dbcsr_finalize(matrix%real_part)
      CALL dbcsr_finalize(matrix%imaginary_part)
   END SUBROUTINE make_operator

END PROGRAM spectral_localizer_flatten_unittest

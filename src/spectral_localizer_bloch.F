!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Exact spectral flattening of a translation-invariant AO torus in primitive Bloch blocks.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE spectral_localizer_bloch
   USE cp_dbcsr_api, ONLY: dbcsr_add, dbcsr_create, dbcsr_distribution_get, &
                           dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_info, &
                           dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
                           dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_put_block, dbcsr_type
   USE cp_dbcsr_contrib, ONLY: dbcsr_frobenius_norm
   USE ieee_arithmetic, ONLY: ieee_is_finite
#if defined(__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: IEEE_ALL, ieee_get_halting_mode, ieee_set_halting_mode
#endif
   USE kinds, ONLY: dp
   USE kpoint_lattice_fft, ONLY: k_grid_to_cell_fft
   USE message_passing, ONLY: mp_para_env_type
   USE spectral_localizer_dbcsr, ONLY: localizer_dbcsr_check, localizer_dbcsr_release, localizer_matrix_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: localizer_flatten_bloch

CONTAINS

! **************************************************************************************************
!> \brief Flatten all Bloch bands, then return to the original distributed real-cell/spin AO basis.
!> \param h Original complex Hamiltonian in spin/cell/atom order, unchanged
!> \param s Original complex metric in the same order, unchanged
!> \param grid Full Gamma-centered torus dimensions
!> \param nspin Number of physical spin blocks (one or two)
!> \param energy Absolute query energy
!> \param scale Positive flattened energy
!> \param tolerance Structure, covariance and eigensystem tolerance
!> \param metric_tolerance Minimum global overlap eigenvalue ratio
!> \param gap_tolerance Minimum original electronic distance from the query
!> \param para_env Communicator; independent eigensystems are distributed across its ranks
!> \param flattened Newly allocated covariant shifted Hamiltonian only on success
!> \param gap Original electronic gap, not a flattened localizer gap
!> \param metric_ratio Global minimum/maximum overlap eigenvalue ratio across all k-points
!> \param residuals Covariance, metric orthonormality, eigenproblem and final Hermiticity residuals
!> \param status 0 success, 1 unresolved gap, -1 invalid layout/data, -2 metric,
!>        -4 eigensolver, -6 excessive residual, -7 broken translation covariance
!> \note No occupied-state projection, symmetry reduction or modification of coordinate operators.
! **************************************************************************************************
   SUBROUTINE localizer_flatten_bloch(h, s, grid, nspin, energy, scale, tolerance, metric_tolerance, &
                                      gap_tolerance, para_env, flattened, gap, metric_ratio, residuals, status)
      TYPE(localizer_matrix_type), INTENT(INOUT)         :: h, s
      INTEGER, INTENT(IN)                                :: grid(3), nspin
      REAL(dp), INTENT(IN)                               :: energy, scale, tolerance, &
                                                            metric_tolerance, gap_tolerance
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(localizer_matrix_type), INTENT(OUT)           :: flattened
      REAL(dp), INTENT(OUT)                              :: gap, metric_ratio, residuals(4)
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp), ALLOCATABLE                           :: fk(:, :, :), hc(:, :, :), hk(:, :, :), &
                                                            sc(:, :, :), sk(:, :, :)
      INTEGER                                            :: ik, nb, ncell, ns, nscalar, nsmall, &
                                                            substatus
      INTEGER, ALLOCATABLE                               :: cell_index(:), cells(:, :), starts(:)
      INTEGER, POINTER                                   :: sizes(:), ssize(:)
      LOGICAL                                            :: valid
      REAL(dp)                                           :: all_bounds(2), bounds(2), difference, &
                                                            errors(2), gap_threshold, local_gap, &
                                                            threshold
      REAL(dp), ALLOCATABLE                              :: weights(:), xkp(:, :)
      TYPE(localizer_matrix_type)                        :: reconstructed

      status = -1
      gap = 0.0_dp
      metric_ratio = 0.0_dp
      residuals = HUGE(1.0_dp)
      IF (ANY(grid < 1) .OR. (nspin /= 1 .AND. nspin /= 2)) RETURN
      IF (.NOT. ALL(ieee_is_finite([energy, scale, tolerance, metric_tolerance, gap_tolerance]))) RETURN
      IF (scale <= 0.0_dp .OR. tolerance <= 0.0_dp .OR. tolerance >= 1.0_dp) RETURN
      IF (metric_tolerance <= 0.0_dp .OR. metric_tolerance >= 1.0_dp .OR. gap_tolerance < 0.0_dp) RETURN
      CALL dbcsr_get_info(h%real_part, nblkrows_total=nb, row_blk_size=sizes)
      CALL dbcsr_get_info(s%real_part, nblkrows_total=ns, row_blk_size=ssize)
      IF (nb < 1 .OR. nb /= ns) RETURN
      IF (ANY(sizes /= ssize)) RETURN
      IF (ANY(sizes <= 0)) RETURN
      ! Bound the products before interpreting the repeated atom-block layout.
      ncell = 1
      DO ik = 1, 3
         IF (grid(ik) > nb/ncell) RETURN
         ncell = ncell*grid(ik)
      END DO
      IF (nspin > nb/ncell) RETURN
      IF (MOD(nb, nspin*ncell) /= 0) RETURN
      ns = nb/(nspin*ncell)
      DO ik = 1, nb
         IF (sizes(ik) /= sizes(1 + MODULO(ik - 1, ns))) RETURN
      END DO
      CALL localizer_dbcsr_check(h, para_env, tolerance, valid)
      IF (.NOT. valid) RETURN
      CALL localizer_dbcsr_check(s, para_env, tolerance, valid)
      IF (.NOT. valid) RETURN
      nscalar = SUM(sizes(:ns))
      nsmall = nscalar*nspin
      ALLOCATE (starts(nb), cell_index(nb), cells(3, ncell), xkp(3, ncell), weights(ncell))
      DO ik = 1, nb
         starts(ik) = 1 + SUM(sizes(:MODULO(ik - 1, ns))) + ((ik - 1)/(ns*ncell))*nscalar
         cell_index(ik) = 1 + MODULO((ik - 1)/ns, ncell)
      END DO
      DO ik = 1, ncell
         cells(:, ik) = [MODULO(ik - 1, grid(1)), MODULO((ik - 1)/grid(1), grid(2)), &
                         (ik - 1)/(grid(1)*grid(2))]
         xkp(:, ik) = REAL(cells(:, ik), dp)/REAL(grid, dp)
      END DO
      ALLOCATE (hc(nsmall, nsmall, ncell), sc(nsmall, nsmall, ncell), &
                hk(nsmall, nsmall, ncell), sk(nsmall, nsmall, ncell), fk(nsmall, nsmall, ncell))
      CALL torus_blocks(h, starts, cell_index, para_env, hc)
      CALL torus_blocks(s, starts, cell_index, para_env, sc)
      residuals = 0.0_dp
      ! Check the complete matrix, including missing blocks, before using one cell as representative.
      CALL expand_blocks(hc, h, grid, starts, cell_index, para_env, reconstructed)
      CALL matrix_difference(reconstructed, h, difference)
      CALL localizer_dbcsr_release(reconstructed)
      residuals(1) = difference
      CALL expand_blocks(sc, s, grid, starts, cell_index, para_env, reconstructed)
      CALL matrix_difference(reconstructed, s, difference)
      CALL localizer_dbcsr_release(reconstructed)
      residuals(1) = MAX(residuals(1), difference)
      status = -7
      IF (.NOT. ieee_is_finite(residuals(1)) .OR. residuals(1) > tolerance) RETURN
      ! B(R) -> B(k) has positive phase and no normalization; the inverse has negative phase / Nk.
      weights = 1.0_dp
      CALL k_grid_to_cell_fft(hc, -xkp, grid, cells, hk, weights=weights)
      CALL k_grid_to_cell_fft(sc, -xkp, grid, cells, sk, weights=weights)
      fk = 0.0_dp
      all_bounds = [HUGE(1.0_dp), 0.0_dp]
      gap = HUGE(1.0_dp)
      gap_threshold = gap_tolerance
      status = 0
      DO ik = 1, ncell
         IF (MODULO(ik - 1, para_env%num_pe) /= para_env%mepos) CYCLE
         CALL flatten_block(hk(:, :, ik), sk(:, :, ik), energy, scale, tolerance, metric_tolerance, fk(:, :, ik), &
                            local_gap, bounds, threshold, errors, substatus)
         IF (substatus /= 0) THEN
            status = MIN(status, substatus)
            CYCLE
         END IF
         all_bounds = [MIN(all_bounds(1), bounds(1)), MAX(all_bounds(2), bounds(2))]
         gap = MIN(gap, local_gap)
         gap_threshold = MAX(gap_threshold, threshold)
         residuals(2:3) = MAX(residuals(2:3), errors)
      END DO
      CALL para_env%min(status)
      IF (status /= 0) RETURN
      CALL para_env%min(all_bounds(1))
      CALL para_env%max(all_bounds(2))
      CALL para_env%min(gap)
      CALL para_env%max(gap_threshold)
      CALL para_env%max(residuals)
      metric_ratio = all_bounds(1)/all_bounds(2)
      status = -2
      IF (metric_ratio <= metric_tolerance) RETURN
      status = 1
      IF (gap <= gap_threshold) RETURN
      status = -6
      IF (MAXVAL(residuals) > tolerance .OR. .NOT. ALL(ieee_is_finite(residuals))) RETURN
      CALL para_env%sum(fk)
      CALL k_grid_to_cell_fft(fk, xkp, grid, cells, hc)
      CALL expand_blocks(hc, h, grid, starts, cell_index, para_env, flattened)
      CALL localizer_dbcsr_check(flattened, para_env, tolerance, valid, residuals(4))
      IF (.NOT. valid) THEN
         CALL localizer_dbcsr_release(flattened)
         RETURN
      END IF
      status = 0
   END SUBROUTINE localizer_flatten_bloch

! **************************************************************************************************
!> \brief Extract the distributed first bra cell for every physical spin; no full AO gathering.
!> \param matrix Original matrix
!> \param starts Primitive spin/AO offsets for each atomic block
!> \param cell_index Cell of each block
!> \param para_env Communicator
!> \param blocks One primitive spin matrix for every cell displacement
! **************************************************************************************************
   SUBROUTINE torus_blocks(matrix, starts, cell_index, para_env, blocks)
      TYPE(localizer_matrix_type), INTENT(INOUT) :: matrix
      INTEGER, INTENT(IN) :: starts(:), cell_index(:)
      TYPE(mp_para_env_type), POINTER :: para_env
      COMPLEX(dp), INTENT(OUT) :: blocks(:, :, :)

      CALL extract(matrix%real_part, CMPLX(1.0_dp, 0.0_dp, dp), .TRUE.)
      CALL extract(matrix%imaginary_part, CMPLX(0.0_dp, 1.0_dp, dp), .FALSE.)
      CALL para_env%sum(blocks)

   CONTAINS

! **************************************************************************************************
!> \brief Accumulate one real component into the complex block row.
!> \param component Real DBCSR component
!> \param factor Real or imaginary unit
!> \param initialize Reset output before accumulation
! **************************************************************************************************
      SUBROUTINE extract(component, factor, initialize)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: component
      COMPLEX(dp), INTENT(IN)                            :: factor
      LOGICAL, INTENT(IN)                                :: initialize

      INTEGER                                            :: c, col, i, j, row
      LOGICAL                                            :: transposed
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator

         IF (initialize) blocks = 0.0_dp
         CALL dbcsr_iterator_start(iterator, component)
         DO WHILE (dbcsr_iterator_blocks_left(iterator))
            CALL dbcsr_iterator_next_block(iterator, row, col, block, transposed=transposed)
            CPASSERT(.NOT. transposed)
            IF (cell_index(row) /= 1) CYCLE
            i = starts(row)
            j = starts(col)
            c = cell_index(col)
            blocks(i:i + SIZE(block, 1) - 1, j:j + SIZE(block, 2) - 1, c) = &
               blocks(i:i + SIZE(block, 1) - 1, j:j + SIZE(block, 2) - 1, c) + factor*block
         END DO
         CALL dbcsr_iterator_stop(iterator)
      END SUBROUTINE extract
   END SUBROUTINE torus_blocks

! **************************************************************************************************
!> \brief Expand one block row into the original DBCSR distribution without changing spin ordering.
!> \param blocks Primitive spin matrices at each displacement
!> \param template Original full-storage AO layout
!> \param grid Torus dimensions
!> \param starts Primitive AO starts
!> \param cell_index Cell of each atomic block
!> \param para_env Communicator
!> \param matrix Newly allocated full torus matrix
! **************************************************************************************************
   SUBROUTINE expand_blocks(blocks, template, grid, starts, cell_index, para_env, matrix)
      COMPLEX(dp), INTENT(IN) :: blocks(:, :, :)
      TYPE(localizer_matrix_type), INTENT(IN) :: template
      INTEGER, INTENT(IN) :: grid(3), starts(:), cell_index(:)
      TYPE(mp_para_env_type), POINTER :: para_env
      TYPE(localizer_matrix_type), INTENT(OUT) :: matrix

      INTEGER :: row, col, i, j, c, d(3), t, u
      INTEGER, POINTER :: sizes(:), rows(:), cols(:), pgrid(:, :)
      TYPE(dbcsr_distribution_type) :: distribution

      CALL dbcsr_get_info(template%real_part, distribution=distribution, row_blk_size=sizes)
      CALL dbcsr_distribution_get(distribution, row_dist=rows, col_dist=cols, pgrid=pgrid)
      CALL dbcsr_create(matrix%real_part, template=template%real_part)
      CALL dbcsr_create(matrix%imaginary_part, template=template%imaginary_part)
      DO row = 1, SIZE(starts)
         t = cell_index(row) - 1
         i = starts(row)
         DO col = 1, SIZE(starts)
            IF (pgrid(rows(row), cols(col)) /= para_env%mepos) CYCLE
            u = cell_index(col) - 1
            j = starts(col)
            d = MODULO([MODULO(u, grid(1)) - MODULO(t, grid(1)), &
                        MODULO(u/grid(1), grid(2)) - MODULO(t/grid(1), grid(2)), &
                        u/(grid(1)*grid(2)) - t/(grid(1)*grid(2))], grid)
            c = 1 + d(1) + grid(1)*(d(2) + grid(2)*d(3))
            CALL dbcsr_put_block(matrix%real_part, row, col, &
                                 REAL(blocks(i:i + sizes(row) - 1, j:j + sizes(col) - 1, c), dp))
            CALL dbcsr_put_block(matrix%imaginary_part, row, col, &
                                 AIMAG(blocks(i:i + sizes(row) - 1, j:j + sizes(col) - 1, c)))
         END DO
      END DO
      CALL dbcsr_finalize(matrix%real_part)
      CALL dbcsr_finalize(matrix%imaginary_part)
   END SUBROUTINE expand_blocks

! **************************************************************************************************
!> \brief Relative complex Frobenius difference; leaves the trial as the difference.
!> \param trial Reconstructed matrix
!> \param reference Original matrix
!> \param error Relative difference
! **************************************************************************************************
   SUBROUTINE matrix_difference(trial, reference, error)
      TYPE(localizer_matrix_type), INTENT(INOUT)         :: trial, reference
      REAL(dp), INTENT(OUT)                              :: error

      REAL(dp)                                           :: norm

      norm = HYPOT(dbcsr_frobenius_norm(reference%real_part), dbcsr_frobenius_norm(reference%imaginary_part))
      CALL dbcsr_add(trial%real_part, reference%real_part, 1.0_dp, -1.0_dp)
      CALL dbcsr_add(trial%imaginary_part, reference%imaginary_part, 1.0_dp, -1.0_dp)
      error = HYPOT(dbcsr_frobenius_norm(trial%real_part), dbcsr_frobenius_norm(trial%imaginary_part))/ &
              MAX(norm, TINY(1.0_dp))
   END SUBROUTINE matrix_difference

! **************************************************************************************************
!> \brief Full generalized Bloch eigensystem, with metric/gap and numerical residual diagnostics.
!> \param h Hermitian primitive Bloch Hamiltonian
!> \param s Hermitian positive metric
!> \param energy Absolute query energy
!> \param scale Flattened band energy
!> \param tolerance Structure tolerance
!> \param metric_tolerance Minimum local overlap eigenvalue ratio
!> \param flat Covariant shifted sign operator
!> \param gap Distance of the original spectrum from energy
!> \param bounds Smallest/largest original metric eigenvalues
!> \param threshold Roundoff gap threshold
!> \param errors Metric orthonormality and relative eigenproblem residual
!> \param status 0 success, -1 invalid data, -2 metric, -4 eigensolver
! **************************************************************************************************
   SUBROUTINE flatten_block(h, s, energy, scale, tolerance, metric_tolerance, flat, gap, bounds, threshold, errors, status)
      COMPLEX(dp), INTENT(IN) :: h(:, :), s(:, :)
      REAL(dp), INTENT(IN) :: energy, scale, tolerance, metric_tolerance
      COMPLEX(dp), INTENT(OUT) :: flat(:, :)
      REAL(dp), INTENT(OUT) :: gap, bounds(2), threshold, errors(2)
      INTEGER, INTENT(OUT) :: status

      COMPLEX(dp), ALLOCATABLE :: a(:, :), b(:, :), work(:), sc(:, :), residual(:, :)
      REAL(dp), ALLOCATABLE :: eigenvalues(:), rwork(:)
      REAL(dp) :: hnorm, snorm, herror, serror
      INTEGER :: n, info, i
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL :: halt(SIZE(IEEE_ALL))
#endif

      status = -1
      n = SIZE(h, 1)
      IF (.NOT. ALL(ieee_is_finite(REAL(h, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(h)))) RETURN
      IF (.NOT. ALL(ieee_is_finite(REAL(s, dp))) .OR. .NOT. ALL(ieee_is_finite(AIMAG(s)))) RETURN
      hnorm = MAX(MAXVAL(ABS(h)), TINY(1.0_dp))
      snorm = MAX(MAXVAL(ABS(s)), TINY(1.0_dp))
      herror = MAXVAL(ABS(h - CONJG(TRANSPOSE(h))))
      serror = MAXVAL(ABS(s - CONJG(TRANSPOSE(s))))
      IF (herror > tolerance*hnorm .OR. serror > tolerance*snorm) RETURN
      ALLOCATE (a(n, n), b(n, n), work(MAX(1, 2*n)), eigenvalues(n), rwork(MAX(1, 3*n)))
      a = s
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zheev('N', 'U', n, a, n, eigenvalues, work, SIZE(work), rwork, info)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      status = -4
      IF (info /= 0 .OR. .NOT. ALL(ieee_is_finite(eigenvalues))) RETURN
      bounds = [MINVAL(eigenvalues), MAXVAL(eigenvalues)]
      status = -2
      IF (bounds(2) <= 0.0_dp .OR. bounds(1) <= metric_tolerance*bounds(2)) RETURN
      a = h - energy*s
      b = s
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zhegv(1, 'V', 'U', n, a, n, b, n, eigenvalues, work, SIZE(work), rwork, info)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      status = -4
      IF (info /= 0 .OR. .NOT. ALL(ieee_is_finite(eigenvalues))) RETURN
      gap = MINVAL(ABS(eigenvalues))
      threshold = REAL(n, dp)*(herror + (ABS(energy) + MAXVAL(ABS(eigenvalues)))*serror)/bounds(1) + &
                  100.0_dp*EPSILON(1.0_dp)*REAL(n, dp)* &
                  MAX(MAXVAL(ABS(eigenvalues)), MAXVAL(ABS(h - energy*s))/bounds(1))
      sc = MATMUL(s, a)
      residual = MATMUL(CONJG(TRANSPOSE(a)), sc)
      DO i = 1, n
         residual(i, i) = residual(i, i) - 1.0_dp
      END DO
      errors(1) = SQRT(SUM(ABS(residual)**2))
      residual = MATMUL(h - energy*s, a) - sc*SPREAD(eigenvalues, 1, n)
      errors(2) = SQRT(SUM(ABS(residual)**2))/MAX(TINY(1.0_dp), &
                                                  SQRT(SUM(ABS(h - energy*s)**2))*SQRT(SUM(ABS(a)**2)))
      flat = MATMUL(sc*SPREAD(SIGN(scale, eigenvalues), 1, n), CONJG(TRANSPOSE(sc)))
      status = 0
   END SUBROUTINE flatten_block
END MODULE spectral_localizer_bloch

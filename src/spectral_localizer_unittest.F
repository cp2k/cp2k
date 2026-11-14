!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE spectral_localizer_test_cases
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: check_localizer_2d, check_localizer_3d, check_localizer_flattening
CONTAINS

! **************************************************************************************************
!> \brief Representative open/periodic Chern and AII models with analytic gap controls.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
   SUBROUTINE check_localizer_2d()
      USE kinds, ONLY: dp
      USE spectral_localizer, ONLY: localizer_evaluate, localizer_periodic_evaluate
      IMPLICIT NONE
      COMPLEX(dp) :: h0(2, 2), s0(2, 2), r0(2, 2, 2), zero(2, 2, 2)
      COMPLEX(dp), ALLOCATABLE :: h(:, :), s(:, :), r(:, :, :), hs(:, :), ss(:, :), &
                                  rs(:, :, :), cosine(:, :, :), sine(:, :, :)
      INTEGER :: i, j, n, index, status
      REAL(dp) :: gap, ratio, reference, two_pi

      ! Spin-degenerate atomic pencil: generalized gap sqrt((2-E)**2+kappa**2*|r-point|**2).
      s0 = 0.0_dp
      s0(1, 1) = 2.0_dp
      s0(2, 2) = 2.0_dp
      h0 = 2*s0
      r0(:, :, 1) = 3*s0
      r0(:, :, 2) = 4*s0
      CALL localizer_evaluate(h0, s0, r0, [0.0_dp, 0.0_dp], 1.0_dp, 0.2_dp, .TRUE., &
                              1.e-10_dp, 1.e-10_dp, 1.e-9_dp, index, gap, ratio, status)
      IF (status /= 0 .OR. index /= 0 .OR. ABS(gap - SQRT(2.0_dp)) > 1.e-12_dp) THEN
         ERROR STOP 'Analytic generalized localizer failed'
      END IF
      CALL localizer_evaluate(h0, s0, r0, [3.0_dp, 4.0_dp], 2.0_dp, 0.2_dp, .TRUE., &
                              1.e-10_dp, 1.e-10_dp, 1.e-9_dp, index, gap, ratio, status)
      IF (status /= 1) ERROR STOP 'Zero gap was assigned an index'
      s0(1, 1) = -2.0_dp
      CALL localizer_evaluate(h0, s0, r0, [0.0_dp, 0.0_dp], 0.0_dp, 0.2_dp, .FALSE., &
                              1.e-10_dp, 1.e-10_dp, 1.e-9_dp, index, gap, ratio, status)
      IF (status /= -2) ERROR STOP 'Indefinite metric accepted'

      two_pi = 8.0_dp*ATAN(1.0_dp)
      s0(1, 1) = 1.0_dp
      s0(2, 2) = 1.0_dp
      h0 = 3*s0
      r0(:, :, 1) = s0
      r0(:, :, 2) = s0
      zero = 0.0_dp
      ! At a pi/2 origin, the periodic atomic gap is sqrt((3-1)**2+1).
      CALL localizer_periodic_evaluate(h0, s0, r0, zero, [two_pi/4, 0.0_dp], 0.0_dp, 1.0_dp, .TRUE., &
                                       1.e-10_dp, 1.e-10_dp, 1.e-9_dp, index, gap, ratio, status)
      IF (status /= 0 .OR. index /= 0 .OR. ABS(gap - SQRT(5.0_dp)) > 1.e-12_dp) ERROR STOP 'Periodic atomic gap failed'

      CALL qwz_model(5, -1.0_dp, .FALSE., h, s, r)
      CALL localizer_evaluate(h, s, r, [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .FALSE., &
                              1.e-10_dp, 1.e-10_dp, 1.e-9_dp, index, reference, ratio, status)
      IF (status /= 0 .OR. index /= -1 .OR. ABS(reference - 0.448472312212568_dp) > 1.e-11_dp) THEN
         ERROR STOP 'Nontrivial open QWZ model failed'
      END IF
      CALL localizer_evaluate(h, s, r(:, :, [2, 1]), [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .FALSE., &
                              1.e-10_dp, 1.e-10_dp, 1.e-9_dp, index, gap, ratio, status)
      IF (status /= 0 .OR. index /= 1 .OR. ABS(gap - reference) > 1.e-11_dp) ERROR STOP 'Oriented-plane reversal failed'

      ! Time-reversed Chern copies with TR-preserving, spin-nonconserving mixing.
      n = SIZE(h, 1)
      ALLOCATE (hs(2*n, 2*n), ss(2*n, 2*n), rs(2*n, 2*n, 2))
      hs = 0.0_dp
      ss = 0.0_dp
      rs = 0.0_dp
      hs(:n, :n) = h
      hs(n + 1:, n + 1:) = CONJG(h)
      ss(:n, :n) = s
      ss(n + 1:, n + 1:) = s
      rs(:n, :n, :) = r
      rs(n + 1:, n + 1:, :) = r
      DO i = 1, n, 2
         hs(i, n + i + 1) = 0.1_dp
         hs(i + 1, n + i) = -0.1_dp
      END DO
      hs(n + 1:, :n) = CONJG(TRANSPOSE(hs(:n, n + 1:)))
      CALL localizer_evaluate(hs, ss, rs, [0.0_dp, 0.0_dp], 0.0_dp, 0.25_dp, .TRUE., &
                              1.e-10_dp, 1.e-10_dp, 1.e-9_dp, index, gap, ratio, status)
      IF (status /= 0 .OR. index /= 1) ERROR STOP 'Spin-mixed nontrivial AII model failed'

      CALL qwz_model(6, -1.0_dp, .TRUE., h, s, r)
      n = SIZE(h, 1)
      ALLOCATE (cosine(n, n, 2), sine(n, n, 2))
      cosine = 0.0_dp
      sine = 0.0_dp
      DO j = 1, 2
         DO i = 1, n
            cosine(i, i, j) = COS(two_pi*REAL(r(i, i, j), dp)/6)
            sine(i, i, j) = SIN(two_pi*REAL(r(i, i, j), dp)/6)
         END DO
      END DO
      CALL localizer_periodic_evaluate(h, s, cosine, sine, [0.0_dp, 0.0_dp], 0.0_dp, 1.0_dp, .FALSE., &
                                       1.e-10_dp, 1.e-10_dp, 1.e-9_dp, index, gap, ratio, status)
      IF (status /= 0 .OR. index /= -1) ERROR STOP 'Periodic QWZ Chern phase failed'
   CONTAINS

! **************************************************************************************************
!> \brief Qi-Wu-Zhang square in site-major, two-orbital order.
!> \param length Sites along each direction
!> \param mass Onsite mass
!> \param torus Whether to close the boundary hoppings
!> \param h Hamiltonian
!> \param s Identity metric
!> \param r Centered Cartesian positions
! **************************************************************************************************
      SUBROUTINE qwz_model(length, mass, torus, h, s, r)
      INTEGER, INTENT(IN)                                :: length
      REAL(dp), INTENT(IN)                               :: mass
      LOGICAL, INTENT(IN)                                :: torus
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: h(:, :), s(:, :), r(:, :, :)

      COMPLEX(dp)                                        :: tx(2, 2), ty(2, 2)
      INTEGER                                            :: i, j, n, orbital, x, y

         n = 2*length**2
         ALLOCATE (h(n, n), s(n, n), r(n, n, 2))
         h = 0.0_dp
         s = 0.0_dp
         r = 0.0_dp
         tx = 0.0_dp
         tx(1, 1) = 0.5_dp
         tx(2, 2) = -0.5_dp
         tx(1, 2) = CMPLX(0.0_dp, -0.5_dp, dp)
         tx(2, 1) = tx(1, 2)
         ty = 0.0_dp
         ty(1, 1) = 0.5_dp
         ty(2, 2) = -0.5_dp
         ty(1, 2) = -0.5_dp
         ty(2, 1) = 0.5_dp
         DO y = 1, length
            DO x = 1, length
               i = 2*((y - 1)*length + x) - 1
               h(i, i) = mass
               h(i + 1, i + 1) = -mass
               DO orbital = 0, 1
                  s(i + orbital, i + orbital) = 1.0_dp
                  r(i + orbital, i + orbital, 1) = REAL(x, dp) - REAL(length + 1, dp)/2
                  r(i + orbital, i + orbital, 2) = REAL(y, dp) - REAL(length + 1, dp)/2
               END DO
               IF (x < length .OR. torus) THEN
                  j = 2*((y - 1)*length + MODULO(x, length) + 1) - 1
                  h(i:i + 1, j:j + 1) = tx
                  h(j:j + 1, i:i + 1) = CONJG(TRANSPOSE(tx))
               END IF
               IF (y < length .OR. torus) THEN
                  j = 2*(MODULO(y, length)*length + x) - 1
                  h(i:i + 1, j:j + 1) = ty
                  h(j:j + 1, i:i + 1) = CONJG(TRANSPOSE(ty))
               END IF
            END DO
         END DO
      END SUBROUTINE qwz_model
   END SUBROUTINE check_localizer_2d

! **************************************************************************************************
!> \brief Strong/trivial Wilson-Dirac models, independent singular gaps and QR fallback.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
   SUBROUTINE check_localizer_3d()
#if defined(__HAS_IEEE_EXCEPTIONS)
      USE ieee_exceptions, ONLY: IEEE_ALL, ieee_get_halting_mode, ieee_set_halting_mode
#endif
      USE kinds, ONLY: dp
      USE spectral_localizer, ONLY: localizer_3d_evaluate, localizer_3d_matrix_evaluate
      IMPLICIT NONE
      COMPLEX(dp), PARAMETER :: ii = CMPLX(0.0_dp, 1.0_dp, dp)
      INTEGER, PARAMETER :: side = 3, sites = side**3, n = 4*sites
      REAL(dp), PARAMETER :: masses(2) = [4.0_dp, 2.0_dp]
      INTEGER, PARAMETER :: indices(2) = [0, 1]
      COMPLEX(dp) :: h(n, n), s(n, n), r(n, n, 3)
      COMPLEX(dp) :: h0(2, 2), s0(2, 2), r0(2, 2, 3)
      COMPLEX(dp) :: pauli(2, 2, 3), a(2*n, 2*n), dummy(1, 1), work(10*n)
      REAL(dp) :: singular(2*n), rwork(10*n), gap, reference, ratio, mass
      INTEGER :: i, j, axis, model, index, status, io
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL :: halt(SIZE(IEEE_ALL)), current_halt(SIZE(IEEE_ALL))
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
#endif

      io = cp_logger_get_default_io_unit()
      CALL determinant_controls()
      pauli = 0.0_dp
      pauli(1, 2, 1) = 1.0_dp
      pauli(2, 1, 1) = 1.0_dp
      pauli(1, 2, 2) = -ii
      pauli(2, 1, 2) = ii
      pauli(1, 1, 3) = 1.0_dp
      pauli(2, 2, 3) = -1.0_dp
      DO model = 1, SIZE(masses)
         mass = masses(model)
         CALL dirac_model(mass)
         CALL evaluate(h, s, r, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp)
         IF (status /= 0 .OR. index /= indices(model)) ERROR STOP '3D strong/trivial parity mismatch'
         reference = gap
      END DO
      ! Independently check the singular gap once, on the nontrivial mass-2 model.
      a = 0.0_dp
      a(:n, :n) = h
      a(n + 1:, n + 1:) = h
      DO axis = 1, 3
         DO j = 1, 2
            DO i = 1, 2
               a((i - 1)*n + 1:i*n, (j - 1)*n + 1:j*n) = &
                  a((i - 1)*n + 1:i*n, (j - 1)*n + 1:j*n) - 0.5_dp*ii*pauli(i, j, axis)*r(:, :, axis)
            END DO
         END DO
      END DO
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
      CALL zgesvd('N', 'N', 2*n, 2*n, a, 2*n, singular, dummy, 1, dummy, 1, work, SIZE(work), rwork, status)
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
      IF (status /= 0 .OR. ABS(reference - MINVAL(singular)) > 2.e-12_dp) ERROR STOP '3D singular gap mismatch'
      IF (io > 0) WRITE (io, '(A,F4.1,I3,2ES18.9)') '3D Dirac mass, index, gap, SVD error: ', &
         mass, index, reference, ABS(reference - MINVAL(singular))

      ! A spin-degenerate atomic pair has the analytic chiral gap |H|.
      s0 = 0.0_dp
      s0(1, 1) = 1.0_dp
      s0(2, 2) = 1.0_dp
      r0 = 0.0_dp
      h0 = 0.1_dp*s0
      CALL evaluate(h0, s0, r0, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp)
      IF (status /= 0 .OR. index /= 0 .OR. ABS(gap - 0.1_dp) > 1.e-14_dp) ERROR STOP '3D atomic gap failed'
      CALL localizer_3d_evaluate(h0, s0, r0, [0.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, 0.5_dp, &
                                 1.e-10_dp, 1.e-10_dp, 0.2_dp, index, gap, ratio, status)
      IF (status /= 1) ERROR STOP '3D unresolved gap assigned an index'
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, current_halt)
      IF (ANY(halt .NEQV. current_halt)) ERROR STOP '3D localizer changed IEEE exception modes'
#endif
      IF (io > 0) WRITE (io, *) '3D AII determinant and generalized-gap tests passed.'
   CONTAINS
! **************************************************************************************************
!> \brief Chiral-block determinant controls independent of the Wilson-Dirac model.
!> \note Wilkinson's matrix has exponential LU pivot growth despite a large singular gap.
! **************************************************************************************************
      SUBROUTINE determinant_controls()
      INTEGER, PARAMETER                                 :: m = 64, quarter = m/4

      COMPLEX(dp)                                        :: block(m, m), metric(m/2, m/2), &
                                                            unused(1, 1), workspace(5*m)
      INTEGER                                            :: col, col_sign, info, orientation, &
                                                            paired_col, paired_row, RESULT, row, &
                                                            row_sign
      REAL(dp)                                           :: matrix(m, m), measured_gap, &
                                                            metric_condition, real_work(5*m), &
                                                            svd_gap, values(m)

         metric = 0.0_dp
         DO row = 1, m/2
            metric(row, row) = 1.0_dp
         END DO
         matrix = 0.0_dp
         DO row = 1, m
            matrix(row, row) = 1.0_dp
            matrix(row, :row - 1) = -1.0_dp
            matrix(row, m) = 1.0_dp + 0.25_dp*MODULO(row, 3)
         END DO
         DO orientation = 0, 1
            ! Positive last-column entries give a positive determinant. One row reversal changes its sign.
            IF (orientation == 1) matrix(1, :) = -matrix(1, :)
            ! Exact dyadic entries in Q W Q^dagger retain the pivot-growth fixture after conversion.
            DO col = 1, m
               paired_col = (3 - (col - 1)/quarter)*quarter + MODULO(col - 1, quarter) + 1
               col_sign = MERGE(1, -1, col <= quarter .OR. col > 3*quarter)
               DO row = 1, m
                  paired_row = (3 - (row - 1)/quarter)*quarter + MODULO(row - 1, quarter) + 1
                  row_sign = MERGE(1, -1, row <= quarter .OR. row > 3*quarter)
                  block(row, col) = 0.5_dp*(matrix(row, col) + row_sign*col_sign*matrix(paired_row, paired_col) + &
                                            ii*(col_sign*matrix(row, paired_col) - row_sign*matrix(paired_row, col)))
               END DO
            END DO
            CALL localizer_3d_matrix_evaluate(block, metric, 1.e-10_dp, 1.e-10_dp, 1.e-10_dp, &
                                              RESULT, measured_gap, metric_condition, info)
            IF (info /= 0 .OR. RESULT /= orientation) ERROR STOP '3D high-growth determinant sign'
#if defined(__HAS_IEEE_EXCEPTIONS)
            CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif
            CALL zgesvd('N', 'N', m, m, block, m, values, unused, 1, unused, 1, workspace, &
                        SIZE(workspace), real_work, info)
#if defined(__HAS_IEEE_EXCEPTIONS)
            CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif
            svd_gap = MINVAL(values)
            IF (info /= 0 .OR. ABS(measured_gap - svd_gap) > 1.e-12_dp) ERROR STOP '3D high-growth singular gap'
            IF (io > 0) WRITE (io, '(A,I2,2ES18.9)') '3D high-growth determinant, gap, SVD error: ', &
               RESULT, measured_gap, ABS(measured_gap - svd_gap)
         END DO
      END SUBROUTINE determinant_controls

! **************************************************************************************************
!> \brief Fixed tolerances and localizer scale for the small Wilson-Dirac control.
!> \param hh Physical Hamiltonian
!> \param ss Physical overlap
!> \param rr Three Cartesian moments
!> \param point Query position
!> \param energy Query energy
! **************************************************************************************************
      SUBROUTINE evaluate(hh, ss, rr, point, energy)
      COMPLEX(dp), INTENT(IN)                            :: hh(:, :), ss(:, :), rr(:, :, :)
      REAL(dp), INTENT(IN)                               :: point(3), energy

         CALL localizer_3d_evaluate(hh, ss, rr, point, energy, 0.5_dp, 1.e-10_dp, 1.e-10_dp, 1.e-10_dp, &
                                    index, gap, ratio, status)
      END SUBROUTINE evaluate
! **************************************************************************************************
!> \brief Open 3x3x3 Wilson-Dirac cube in fixed physical up/down ordering.
!> \param m Onsite mass; periodic Bloch mass is m+cos(kx)+cos(ky)+cos(kz)
! **************************************************************************************************
      SUBROUTINE dirac_model(m)
      REAL(dp), INTENT(IN)                               :: m

      INTEGER                                            :: col, direction, neighbor, orbital, &
                                                            other, pos(3), row, site, spin, x, y, z

         h = 0.0_dp
         s = 0.0_dp
         r = 0.0_dp
         DO z = 0, side - 1
            DO y = 0, side - 1
               DO x = 0, side - 1
                  pos = [x, y, z]
                  site = 1 + x + side*y + side**2*z
                  DO spin = 1, 2
                     DO orbital = 1, 2
                        row = site + sites*(orbital - 1) + 2*sites*(spin - 1)
                        s(row, row) = 1.0_dp
                        h(row, row) = m*REAL(3 - 2*orbital, dp)
                        r(row, row, :) = REAL(pos - 1, dp)
                        DO direction = 1, 3
                           IF (pos(direction) == side - 1) CYCLE
                           neighbor = site + side**(direction - 1)
                           col = neighbor + sites*(orbital - 1) + 2*sites*(spin - 1)
                           h(row, col) = 0.5_dp*REAL(3 - 2*orbital, dp)
                           h(col, row) = CONJG(h(row, col))
                           DO other = 1, 2
                              col = neighbor + sites*(2 - orbital) + 2*sites*(other - 1)
                              h(row, col) = -0.5_dp*ii*pauli(spin, other, direction)
                              h(col, row) = CONJG(h(row, col))
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END SUBROUTINE dirac_model
   END SUBROUTINE check_localizer_3d

! **************************************************************************************************
!> \brief Distributed and Bloch flattening against shared exact spectral constructions.
!> \param para_env Parallel environment initialized by the unit-test driver
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
   SUBROUTINE check_localizer_flattening(para_env)
      USE cp_dbcsr_api, ONLY: dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, &
                              dbcsr_distribution_type, dbcsr_finalize, dbcsr_finalize_lib, &
                              dbcsr_init_lib, dbcsr_put_block, dbcsr_type_no_symmetry
      USE fft_api, ONLY: FFT_LIBRARY_BACKEND_DEFAULT, FFT_LIBRARY_FFTW_MODE_ESTIMATE
      USE fft_tools, ONLY: finalize_fft, init_fft
      USE kinds, ONLY: dp
      USE mathconstants, ONLY: twopi
      USE message_passing, ONLY: mp_para_env_type
      USE offload_api, ONLY: offload_get_device_count, &
                             offload_set_chosen_device
      USE spectral_localizer_bloch, ONLY: localizer_flatten_bloch
      USE spectral_localizer_dbcsr, ONLY: localizer_dbcsr_dense, localizer_dbcsr_release, localizer_matrix_type
      USE spectral_localizer_flatten, ONLY: localizer_flatten_checked
#include "./base/base_uses.f90"

      IMPLICIT NONE
      INTEGER, PARAMETER :: n = 12, nb = 6, bs = 2
      INTEGER, PARAMETER :: grids(3, 3) = RESHAPE([1, 3, 1, 2, 3, 1, 3, 1, 1], [3, 3]), spins(3) = [2, 1, 1]
      REAL(dp), PARAMETER :: energy = 0.002_dp, scale = 1.7_dp
      TYPE(mp_para_env_type), POINTER :: para_env
      TYPE(dbcsr_distribution_type) :: distribution
      TYPE(localizer_matrix_type) :: h, s, flat
      INTEGER, POINTER :: process_grid(:, :), row_dist(:), col_dist(:), sizes(:)
      COMPLEX(dp) :: transform(n, n), u(n, n), hdense(n, n), sdense(n, n), exact(n, n)
      COMPLEX(dp), ALLOCATABLE :: result(:, :)
      REAL(dp) :: energies(n), residuals(4), gap(2), metric_ratio, ref_gap, bloch_gap, phase, maximum_error
      INTEGER :: band, grid(3), i, io, k, kr(3), m, nc, ncell, nr, nscalar, nsmall, nspin, &
                 p, pnext, r(3), rnext, row, status, t, tnext
      CHARACTER(LEN=32) :: argument

      io = cp_logger_get_default_io_unit()
      ! Select the active offload device before initializing the FFT library.
      IF (offload_get_device_count() > 0) THEN
         CALL offload_set_chosen_device(MOD(para_env%mepos, offload_get_device_count()))
      END IF
      CALL dbcsr_init_lib(para_env%get_handle(), io)
      CALL init_fft(FFT_LIBRARY_BACKEND_DEFAULT, alltoall=.FALSE., pool_limit=4, wisdom_file='', &
                    plan_style=FFT_LIBRARY_FFTW_MODE_ESTIMATE, status=status)
      IF (status /= 0) CPABORT("FFT initalization step failed")
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
      CALL dbcsr_distribution_new(distribution, group=para_env%get_handle(), pgrid=process_grid, &
                                  row_dist=row_dist, col_dist=col_dist)
      maximum_error = 0.0_dp
      ! Spin/cell ordering, a 2D grid and nonuniform atom blocks share one reference generator.
      DO m = 1, SIZE(spins)
         grid = grids(:, m)
         nspin = spins(m)
         sizes = bs
         IF (m == 3) sizes = [1, 3, 1, 3, 1, 3]
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
            IF (m /= 2) transform(row, rnext) = transform(row, rnext) + CMPLX(0.05_dp, 0.02_dp, dp)
            DO k = 0, ncell - 1
               kr = [MODULO(k, grid(1)), MODULO(k/grid(1), grid(2)), k/(grid(1)*grid(2))]
               DO band = 0, nsmall - 1
                  phase = twopi*(SUM(REAL(kr*r, dp)/REAL(grid, dp)) + REAL(p*band, dp)/REAL(nsmall, dp))
                  u(row, k*nsmall + band + 1) = EXP(CMPLX(0.0_dp, phase, dp))/SQRT(REAL(n, dp))
                  energies(k*nsmall + band + 1) = 0.2_dp + 0.03_dp*band + 0.01_dp*COS(REAL(k, dp))
                  IF (band < nsmall/2) energies(k*nsmall + band + 1) = -energies(k*nsmall + band + 1)
               END DO
            END DO
         END DO
         ! B=U^dagger T gives H=B^dagger diag(e)B, S=B^dagger B and an exact covariant sign.
         transform = MATMUL(CONJG(TRANSPOSE(u)), transform)
         hdense = MATMUL(CONJG(TRANSPOSE(transform)), transform*SPREAD(energies, 2, n))
         sdense = MATMUL(CONJG(TRANSPOSE(transform)), transform)
         exact = MATMUL(CONJG(TRANSPOSE(transform)), transform*SPREAD(SIGN(scale, energies - energy), 2, n))
         ref_gap = MINVAL(ABS(energies - energy))
         CALL make_operator(hdense, h)
         CALL make_operator(sdense, s)
         IF (m == 1) THEN
            CALL localizer_flatten_checked(h, s, energy, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                           para_env, flat, gap, metric_ratio, residuals, status)
            IF (status /= 0) CPABORT('Checked flattening failed')
            IF (gap(1) > ref_gap + 1.e-10_dp .OR. gap(2) < ref_gap - 1.e-10_dp .OR. metric_ratio <= 1.e-10_dp) THEN
               CPABORT('Checked flattening metric/gap bound failed')
            END IF
            CALL localizer_dbcsr_dense(flat, para_env, RESULT)
            maximum_error = MAX(maximum_error, MAXVAL(ABS(RESULT - exact)))
            IF (MAXVAL(ABS(RESULT - exact)) > 1.e-9_dp) THEN
               CPABORT('Covariant sign differs from exact eigensystem')
            END IF
            DEALLOCATE (RESULT)
            CALL localizer_dbcsr_release(flat)
         END IF
         CALL localizer_flatten_bloch(h, s, grid, nspin, energy, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                      para_env, flat, bloch_gap, metric_ratio, residuals, status)
         IF (status /= 0) CPABORT('Bloch flattening failed')
         CALL localizer_dbcsr_dense(flat, para_env, RESULT)
         maximum_error = MAX(maximum_error, MAXVAL(ABS(RESULT - exact)))
         IF (MAXVAL(ABS(RESULT - exact)) > 1.e-11_dp .OR. ABS(bloch_gap - ref_gap) > 1.e-12_dp) THEN
            CPABORT('Bloch flattening differs from exact full-torus reference')
         END IF
         DEALLOCATE (RESULT)
         CALL localizer_dbcsr_release(flat)
         CALL localizer_dbcsr_release(h)
         IF (m == 1) THEN
            hdense(1, 1) = hdense(1, 1) + 0.01_dp
            CALL make_operator(hdense, h)
            CALL localizer_flatten_bloch(h, s, grid, nspin, energy, scale, 1.e-9_dp, 1.e-10_dp, 1.e-8_dp, &
                                         para_env, flat, bloch_gap, metric_ratio, residuals, status)
            IF (status /= -7) CPABORT('Broken translations accepted')
            CALL localizer_dbcsr_release(h)
         END IF
         CALL localizer_dbcsr_release(s)
      END DO
      IF (io > 0) WRITE (io, '(A,ES14.5)') 'Spectral flattening exact-matrix maximum error: ', maximum_error
      CALL dbcsr_distribution_release(distribution)
      DEALLOCATE (process_grid, row_dist, col_dist, sizes)
      CALL dbcsr_finalize_lib()
      CALL finalize_fft('')
   CONTAINS

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
   END SUBROUTINE check_localizer_flattening

END MODULE spectral_localizer_test_cases

PROGRAM spectral_localizer_unittest
   USE cp_log_handling, ONLY: cp_add_default_logger, &
                              cp_logger_create, &
                              cp_logger_release, &
                              cp_logger_type, &
                              cp_rm_default_logger
   USE machine, ONLY: default_output_unit
   USE message_passing, ONLY: mp_para_env_type, mp_world_finalize, mp_world_init
   USE spectral_localizer_test_cases, ONLY: check_localizer_2d, check_localizer_3d, check_localizer_flattening
   IMPLICIT NONE
   TYPE(mp_para_env_type), POINTER :: para_env
   TYPE(cp_logger_type), POINTER :: logger
   INTEGER :: io

   ALLOCATE (para_env)
   CALL mp_world_init(para_env)
   io = -1
   IF (para_env%is_source()) io = default_output_unit
   CALL cp_logger_create(logger, para_env=para_env, default_global_unit_nr=io, close_global_unit_on_dealloc=.FALSE.)
   CALL cp_add_default_logger(logger)
   CALL check_localizer_2d()
   CALL check_localizer_3d()
   CALL check_localizer_flattening(para_env)
   CALL cp_rm_default_logger()
   CALL cp_logger_release(logger)
   DEALLOCATE (para_env)
   CALL mp_world_finalize()
END PROGRAM spectral_localizer_unittest

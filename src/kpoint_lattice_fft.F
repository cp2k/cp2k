!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Batched lattice Fourier transforms for complete regular k-point grids.
!>
!>        Real-space image cells are folded modulo the reciprocal grid.  Shifted Monkhorst-Pack
!>        meshes are handled by a phase twist before the FFT.  Non-regular or incomplete k-point
!>        lists transparently use the direct phase sum.
! **************************************************************************************************
MODULE kpoint_lattice_fft
   USE fft_tools,                       ONLY: BWFFT,&
                                              FFT_RADIX_NEXT,&
                                              fft3d,&
                                              fft_alloc,&
                                              fft_dealloc,&
                                              fft_radix_operations,&
                                              fft_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: gaussi,&
                                              twopi,&
                                              z_zero
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_lattice_fft'

   PUBLIC :: cell_to_k_grid_fft, &
             regular_kpoint_grid

CONTAINS

! **************************************************************************************************
!> \brief Test and map a complete, uniformly shifted reciprocal grid.
!> \param xkp reciprocal coordinates, in units of reciprocal lattice vectors
!> \param nkp_grid regular-grid dimensions
!> \param grid_index FFT-grid index for every k point
!> \param k_offset common reciprocal-coordinate offset
!> \return true for a complete regular grid
! **************************************************************************************************
   LOGICAL FUNCTION regular_kpoint_grid(xkp, nkp_grid, grid_index, k_offset) RESULT(regular)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      INTEGER, DIMENSION(3), INTENT(IN)                  :: nkp_grid
      INTEGER, DIMENSION(:, :), INTENT(OUT), OPTIONAL    :: grid_index
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT), OPTIONAL :: k_offset

      REAL(KIND=dp), PARAMETER                           :: map_tolerance = 2.0E-10_dp

      INTEGER                                            :: d, ik, m, nkp
      INTEGER, DIMENSION(3)                              :: index
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :)           :: occupied
      REAL(KIND=dp)                                      :: offset(3), scaled

      regular = .FALSE.
      IF (SIZE(xkp, 1) < 3 .OR. ANY(nkp_grid <= 0)) RETURN
      nkp = SIZE(xkp, 2)
      IF (PRODUCT(nkp_grid) /= nkp) RETURN
      IF (PRESENT(grid_index)) THEN
         IF (SIZE(grid_index, 1) < 3 .OR. SIZE(grid_index, 2) < nkp) RETURN
      END IF

      DO d = 1, 3
         scaled = MODULO(REAL(nkp_grid(d), KIND=dp)*xkp(d, 1), 1.0_dp)
         IF (ABS(scaled) < map_tolerance .OR. ABS(scaled - 1.0_dp) < map_tolerance) scaled = 0.0_dp
         offset(d) = scaled/REAL(nkp_grid(d), KIND=dp)
      END DO

      ALLOCATE (occupied(nkp_grid(1), nkp_grid(2), nkp_grid(3)), source=.FALSE.)
      DO ik = 1, nkp
         DO d = 1, 3
            scaled = REAL(nkp_grid(d), KIND=dp)*(xkp(d, ik) - offset(d))
            m = NINT(scaled)
            IF (ABS(scaled - REAL(m, KIND=dp)) > map_tolerance) THEN
               DEALLOCATE (occupied)
               RETURN
            END IF
            INDEX(d) = MODULO(m, nkp_grid(d)) + 1
         END DO
         IF (occupied(INDEX(1), INDEX(2), INDEX(3))) THEN
            DEALLOCATE (occupied)
            RETURN
         END IF
         occupied(INDEX(1), INDEX(2), INDEX(3)) = .TRUE.
         IF (PRESENT(grid_index)) grid_index(1:3, ik) = index
      END DO
      regular = ALL(occupied)
      DEALLOCATE (occupied)

      IF (regular .AND. PRESENT(k_offset)) k_offset = offset

   END FUNCTION regular_kpoint_grid

! **************************************************************************************************
!> \brief Transform a batch of real matrices from image cells to every supplied k point.
!> \param values_rs real-space matrices; the last dimension enumerates image cells
!> \param index_to_cell integer lattice vector for every image cell
!> \param xkp reciprocal coordinates of all requested k points
!> \param nkp_grid dimensions of the candidate regular reciprocal grid
!> \param values_k complex matrices at every k point
!> \param used_fft reports whether the regular-grid FFT path was used
!> \param deriv_direction optional Cartesian k derivative (1, 2, or 3)
!> \param hmat direct-lattice cell matrix, required for a derivative
!> \param selected_kpoints optional indices of k points to retain in values_k
! **************************************************************************************************
   SUBROUTINE cell_to_k_grid_fft(values_rs, index_to_cell, xkp, nkp_grid, values_k, &
                                 used_fft, deriv_direction, hmat, selected_kpoints)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: values_rs
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      INTEGER, DIMENSION(3), INTENT(IN)                  :: nkp_grid
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: values_k
      LOGICAL, INTENT(OUT), OPTIONAL                     :: used_fft
      INTEGER, INTENT(IN), OPTIONAL                      :: deriv_direction
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN), &
         OPTIONAL                                        :: hmat
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: selected_kpoints

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: cell_factor
      COMPLEX(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: fft_in, fft_out
      INTEGER                                            :: attempt, d, handle, i, icell, ik, &
                                                            ik_out, j, nout, radix_length, stat
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: cell_index, grid_index
      INTEGER, DIMENSION(3)                              :: n, nfft
      LOGICAL                                            :: compatible, regular
      REAL(KIND=dp)                                      :: arg
      REAL(KIND=dp), DIMENSION(3)                        :: cell_vector, k_offset

      CALL timeset("cell_to_k_grid_fft", handle)

      IF (SIZE(index_to_cell, 1) < 3 .OR. &
          SIZE(index_to_cell, 2) /= SIZE(values_rs, 3)) THEN
         CPABORT("Inconsistent real-space image-cell mapping in lattice FFT")
      END IF
      nout = SIZE(xkp, 2)
      IF (PRESENT(selected_kpoints)) THEN
         nout = SIZE(selected_kpoints)
         IF (ANY(selected_kpoints < 1) .OR. ANY(selected_kpoints > SIZE(xkp, 2))) THEN
            CPABORT("Selected lattice-FFT k-point index is out of range")
         END IF
      END IF
      IF (SIZE(values_k, 1) /= SIZE(values_rs, 1) .OR. &
          SIZE(values_k, 2) /= SIZE(values_rs, 2) .OR. &
          SIZE(values_k, 3) /= nout) THEN
         CPABORT("Inconsistent input and output matrix batches in lattice FFT")
      END IF
      IF (PRESENT(deriv_direction)) THEN
         IF (.NOT. PRESENT(hmat)) CPABORT("Lattice-FFT derivative requested without a cell matrix")
         IF (deriv_direction < 1 .OR. deriv_direction > 3) THEN
            CPABORT("Lattice-FFT derivative direction must be 1, 2, or 3")
         END IF
      END IF

      ALLOCATE (grid_index(3, SIZE(xkp, 2)))
      regular = regular_kpoint_grid(xkp, nkp_grid, grid_index, k_offset)
      n = nkp_grid
      nfft = n
      IF (regular .AND. fft_type /= 3) THEN
         DO d = 1, 3
            ! FFTSG rejects some short or non-radix lengths.  Embed the N-point
            ! transform in a supported integer multiple L: real-space samples
            ! occupy every L/N-th slot, while the first N output frequencies
            ! are exactly the desired N-point DFT.  This keeps common 4x4x1
            ! and 2D/1D Monkhorst-Pack grids on the FFT path.
            compatible = .FALSE.
            DO attempt = 0, 15
               IF (nfft(d) >= 3) THEN
                  CALL fft_radix_operations(nfft(d), radix_length, FFT_RADIX_NEXT)
                  IF (radix_length == nfft(d)) THEN
                     compatible = .TRUE.
                     EXIT
                  END IF
               END IF
               nfft(d) = nfft(d) + n(d)
            END DO
            IF (.NOT. compatible) regular = .FALSE.
         END DO
      END IF
      IF (.NOT. regular) THEN
         CALL direct_cell_to_k(values_rs, index_to_cell, xkp, values_k, deriv_direction, hmat, &
                               selected_kpoints)
         IF (PRESENT(used_fft)) used_fft = .FALSE.
         DEALLOCATE (grid_index)
         CALL timestop(handle)
         RETURN
      END IF

      NULLIFY (fft_in, fft_out)
      CALL fft_alloc(fft_in, nfft)
      CALL fft_alloc(fft_out, nfft)
      ALLOCATE (cell_index(3, SIZE(values_rs, 3)), cell_factor(SIZE(values_rs, 3)))
      DO icell = 1, SIZE(values_rs, 3)
         DO d = 1, 3
            cell_index(d, icell) = MODULO(index_to_cell(d, icell), n(d))*(nfft(d)/n(d)) + 1
         END DO
         arg = SUM(k_offset*REAL(index_to_cell(1:3, icell), KIND=dp))
         cell_factor(icell) = EXP(gaussi*twopi*arg)
         IF (PRESENT(deriv_direction)) THEN
            cell_vector = MATMUL(hmat, REAL(index_to_cell(1:3, icell), KIND=dp))
            cell_factor(icell) = cell_factor(icell)*gaussi*cell_vector(deriv_direction)
         END IF
      END DO

      stat = 0
      DO j = 1, SIZE(values_rs, 2)
         DO i = 1, SIZE(values_rs, 1)
            fft_in = z_zero
            DO icell = 1, SIZE(values_rs, 3)
               fft_in(cell_index(1, icell), cell_index(2, icell), cell_index(3, icell)) = &
                  fft_in(cell_index(1, icell), cell_index(2, icell), cell_index(3, icell)) + &
                  cell_factor(icell)*values_rs(i, j, icell)
            END DO
            CALL fft3d(BWFFT, nfft, fft_in, fft_out, status=stat)
            IF (stat /= 0) EXIT
            DO ik_out = 1, nout
               ik = ik_out
               IF (PRESENT(selected_kpoints)) ik = selected_kpoints(ik_out)
               values_k(i, j, ik_out) = &
                  fft_out(grid_index(1, ik), grid_index(2, ik), grid_index(3, ik))
            END DO
         END DO
         IF (stat /= 0) EXIT
      END DO

      IF (stat /= 0) THEN
         CALL direct_cell_to_k(values_rs, index_to_cell, xkp, values_k, deriv_direction, hmat, &
                               selected_kpoints)
         regular = .FALSE.
      END IF
      IF (PRESENT(used_fft)) used_fft = regular

      CALL fft_dealloc(fft_in)
      CALL fft_dealloc(fft_out)
      DEALLOCATE (cell_factor, cell_index, grid_index)
      CALL timestop(handle)

   END SUBROUTINE cell_to_k_grid_fft

! **************************************************************************************************
!> \brief Direct phase-sum fallback for arbitrary reciprocal-point lists.
!> \param values_rs ...
!> \param index_to_cell ...
!> \param xkp ...
!> \param values_k ...
!> \param deriv_direction ...
!> \param hmat ...
!> \param selected_kpoints ...
! **************************************************************************************************
   SUBROUTINE direct_cell_to_k(values_rs, index_to_cell, xkp, values_k, &
                               deriv_direction, hmat, selected_kpoints)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: values_rs
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: values_k
      INTEGER, INTENT(IN), OPTIONAL                      :: deriv_direction
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN), &
         OPTIONAL                                        :: hmat
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: selected_kpoints

      COMPLEX(KIND=dp)                                   :: factor
      INTEGER                                            :: icell, ik, ik_out
      REAL(KIND=dp)                                      :: arg
      REAL(KIND=dp), DIMENSION(3)                        :: cell_vector

      values_k = z_zero
      DO ik_out = 1, SIZE(values_k, 3)
         ik = ik_out
         IF (PRESENT(selected_kpoints)) ik = selected_kpoints(ik_out)
         DO icell = 1, SIZE(values_rs, 3)
            arg = SUM(xkp(1:3, ik)*REAL(index_to_cell(1:3, icell), KIND=dp))
            factor = EXP(gaussi*twopi*arg)
            IF (PRESENT(deriv_direction)) THEN
               cell_vector = MATMUL(hmat, REAL(index_to_cell(1:3, icell), KIND=dp))
               factor = factor*gaussi*cell_vector(deriv_direction)
            END IF
            values_k(:, :, ik_out) = values_k(:, :, ik_out) + factor*values_rs(:, :, icell)
         END DO
      END DO

   END SUBROUTINE direct_cell_to_k

END MODULE kpoint_lattice_fft

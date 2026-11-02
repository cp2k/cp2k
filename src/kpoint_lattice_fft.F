!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Sparse direct and regular-grid lattice Fourier transforms.
!>
!>        Image cells are folded modulo the reciprocal grid; shifted meshes use a phase twist.
!>        Forward transforms can embed a regular-grid subset. The inverse transform can also
!>        embed an explicitly allowed subset with its supplied weights and zero missing points.
!>        Generic low-level wrappers retain a direct phase-sum fallback for non-RSKP
!>        callers. The selected RSKP FFT backend requires successful FFT execution.
! **************************************************************************************************
MODULE kpoint_lattice_fft
   USE fft_lib,                         ONLY: fft_supports_arbitrary_lengths
   USE fft_tools,                       ONLY: BWFFT,&
                                              FFT_RADIX_NEXT,&
                                              FWFFT,&
                                              fft_1d_many,&
                                              fft_alloc,&
                                              fft_dealloc,&
                                              fft_radix_operations
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathconstants,                   ONLY: gaussi,&
                                              twopi,&
                                              z_zero
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_lattice_fft'

   PUBLIC :: cell_to_k_grid_fft, cell_to_k_grid_prepare, cell_to_k_grid_release, cell_to_k_grid_work_type
   PUBLIC :: k_grid_to_cell_fft, k_grid_to_cell_prepare, k_grid_to_cell_execute, k_grid_to_cell_release
   PUBLIC :: k_grid_to_cell_work_type, lattice_fft_shape, regular_kpoint_grid

   INTEGER(KIND=int_8), PARAMETER :: fft_batch_bytes = 32_int_8*1024_int_8**2

   ! Reusable forward lattice geometry and bounded numerical FFT workspace.
   TYPE cell_to_k_grid_work_type
      LOGICAL :: regular = .FALSE.
      INTEGER :: nfft(3) = 0
      INTEGER, ALLOCATABLE, PRIVATE :: grid_index(:, :), cell_index(:, :)
      COMPLEX(KIND=dp), ALLOCATABLE, PRIVATE :: cell_factor(:)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER, PRIVATE :: fft_in(:) => NULL(), fft_out(:) => NULL()
   END TYPE cell_to_k_grid_work_type

   ! Call-scoped mathematical maps and buffers. Backend plans remain owned by fft_tools.
   TYPE k_grid_to_cell_work_type
      PRIVATE
      LOGICAL :: regular = .FALSE.
      INTEGER :: nfft(3) = 0
      INTEGER, ALLOCATABLE :: cells(:, :), grid_index(:, :), cell_index(:, :)
      REAL(KIND=dp), ALLOCATABLE :: xkp(:, :), weight(:)
      COMPLEX(KIND=dp), ALLOCATABLE :: cell_phase(:)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER :: fft_in(:) => NULL(), fft_out(:) => NULL()
   END TYPE k_grid_to_cell_work_type

CONTAINS

! **************************************************************************************************
!> \brief Test and map a uniformly shifted reciprocal grid or an explicitly allowed subset.
!> \param xkp reciprocal coordinates, in units of reciprocal lattice vectors
!> \param nkp_grid regular-grid dimensions
!> \param grid_index FFT-grid index for every k point
!> \param k_offset common reciprocal-coordinate offset
!> \param allow_incomplete accept a subset of a regular grid (missing points are zero)
!> \return true for a regular grid, complete unless allow_incomplete is set
! **************************************************************************************************
   LOGICAL FUNCTION regular_kpoint_grid(xkp, nkp_grid, grid_index, k_offset, allow_incomplete) RESULT(regular)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      INTEGER, DIMENSION(3), INTENT(IN)                  :: nkp_grid
      INTEGER, DIMENSION(:, :), INTENT(OUT), OPTIONAL    :: grid_index
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT), OPTIONAL :: k_offset
      LOGICAL, INTENT(IN), OPTIONAL                      :: allow_incomplete

      REAL(KIND=dp), PARAMETER                           :: map_tolerance = 2.0E-10_dp

      INTEGER                                            :: d, ik, m, nkp
      INTEGER, DIMENSION(3)                              :: index
      LOGICAL                                            :: require_complete
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :)           :: occupied
      REAL(KIND=dp)                                      :: offset(3), scaled

      regular = .FALSE.
      IF (SIZE(xkp, 1) < 3 .OR. ANY(nkp_grid <= 0)) RETURN
      nkp = SIZE(xkp, 2)
      require_complete = .TRUE.
      IF (PRESENT(allow_incomplete)) require_complete = .NOT. allow_incomplete
      IF (nkp < 1 .OR. nkp > PRODUCT(nkp_grid)) RETURN
      IF (require_complete .AND. PRODUCT(nkp_grid) /= nkp) RETURN
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
            IF (ABS(scaled - REAL(m, KIND=dp)) > map_tolerance) RETURN
            INDEX(d) = MODULO(m, nkp_grid(d)) + 1
         END DO
         IF (occupied(INDEX(1), INDEX(2), INDEX(3))) RETURN
         occupied(INDEX(1), INDEX(2), INDEX(3)) = .TRUE.
         IF (PRESENT(grid_index)) grid_index(1:3, ik) = index
      END DO
      regular = .TRUE.

      IF (regular .AND. PRESENT(k_offset)) k_offset = offset

   END FUNCTION regular_kpoint_grid

! **************************************************************************************************
!> \brief Find backend-supported FFT dimensions that preserve the original lattice grid.
!>        Sample the padded transform at integer multiples of the original cell indices.
!> \param n original grid dimensions
!> \param nfft supported integer multiples, if available
!> \param allow_arbitrary use exact dimensions when supported by the active FFT backend
!> \return whether all dimensions have a supported multiple
! **************************************************************************************************
   LOGICAL FUNCTION lattice_fft_shape(n, nfft, allow_arbitrary) RESULT(compatible)

      INTEGER, DIMENSION(3), INTENT(IN)                  :: n
      INTEGER, DIMENSION(3), INTENT(OUT)                 :: nfft
      LOGICAL, INTENT(IN), OPTIONAL                      :: allow_arbitrary

      INTEGER                                            :: attempt, d, radix_length

      compatible = .FALSE.
      nfft = n
      IF (ANY(n <= 0)) RETURN
      IF (PRESENT(allow_arbitrary)) THEN
         IF (allow_arbitrary .AND. fft_supports_arbitrary_lengths()) THEN
            ! For capable backends, the PW grid radix list is a cost heuristic, not a capability limit.
            compatible = .TRUE.
            RETURN
         END IF
      END IF
      DO d = 1, 3
         DO attempt = 0, 15
            IF (nfft(d) == 1) EXIT
            IF (nfft(d) >= 2) THEN
               CALL fft_radix_operations(nfft(d), radix_length, FFT_RADIX_NEXT)
               IF (radix_length == nfft(d)) EXIT
            END IF
            nfft(d) = nfft(d) + n(d)
         END DO
         IF (attempt > 15) RETURN
      END DO
      compatible = .TRUE.

   END FUNCTION lattice_fft_shape

! **************************************************************************************************
!> \brief Transform all lattice lines in three batched 1D FFTs.
!>        Each axis moves from first to last through the existing transpose interface.
!>        Input order is (x,y,z,batch), result order (batch,x,y,z), returned in buffer.
!> \param sign fft_tools transform direction
!> \param n lattice dimensions
!> \param nbatch number of independent matrix elements
!> \param buffer input, then transformed data; ownership swaps with scratch
!> \param scratch equally sized FFT workspace
!> \param stat nonzero on backend execution failure
! **************************************************************************************************
   SUBROUTINE lattice_fft_execute(sign, n, nbatch, buffer, scratch, stat)
      INTEGER, INTENT(IN)                                :: sign, n(3), nbatch
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: buffer(:), scratch(:)
      INTEGER, INTENT(OUT)                               :: stat

      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: swap(:)
      INTEGER                                            :: d, nlines

      stat = 0
      DO d = 1, 3
         IF (n(d) == 1) CYCLE
         nlines = nbatch*PRODUCT(n)/n(d)
         CALL fft_1d_many(sign, n(d), nlines, .FALSE., .TRUE., n(d), nlines, &
                          buffer, scratch, 1.0_dp, stat)
         IF (stat /= 0) RETURN
         swap => buffer
         buffer => scratch
         scratch => swap
      END DO
   END SUBROUTINE lattice_fft_execute

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
!> \param require_fft abort if FFT cannot be used; never execute a direct phase sum
!> \param work optional prepared geometry, unchanged throughout repeated calls
!> \param cell_first first prepared image-map column for this batch
!> \param max_workspace_bytes optional FFT workspace bound
! **************************************************************************************************
   SUBROUTINE cell_to_k_grid_fft(values_rs, index_to_cell, xkp, nkp_grid, values_k, &
                              used_fft, deriv_direction, hmat, selected_kpoints, require_fft, work, cell_first, max_workspace_bytes)

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
      LOGICAL, INTENT(IN), OPTIONAL                      :: require_fft
      TYPE(cell_to_k_grid_work_type), INTENT(INOUT), &
         OPTIONAL, TARGET                                :: work
      INTEGER, INTENT(IN), OPTIONAL                      :: cell_first
      INTEGER(KIND=int_8), INTENT(IN), OPTIONAL          :: max_workspace_bytes

      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: input(:, :, :, :), output(:, :, :, :)
      INTEGER                                            :: first, handle, i, icell, ik, ik_out, &
                                                            image, j, nbatch, nout, start, stat, &
                                                            width
      INTEGER(KIND=int_8)                                :: budget, needed, ngrid
      LOGICAL                                            :: success
      TYPE(cell_to_k_grid_work_type), POINTER            :: plan
      TYPE(cell_to_k_grid_work_type), TARGET             :: local_work

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
         IF (.NOT. PRESENT(hmat)) THEN
            CALL cp_abort(__LOCATION__, "Lattice-FFT derivative requested without a cell matrix")
         END IF
         IF (deriv_direction < 1 .OR. deriv_direction > 3) THEN
            CPABORT("Lattice-FFT derivative direction must be 1, 2, or 3")
         END IF
      END IF

      IF (PRESENT(work)) THEN
         plan => work
      ELSE
         plan => local_work
         CALL cell_to_k_grid_prepare(plan, index_to_cell, xkp, nkp_grid, deriv_direction, hmat)
      END IF
      IF (.NOT. plan%regular .AND. PRESENT(require_fft)) THEN
         IF (require_fft) CPABORT("Forced lattice FFT needs a supported complete regular grid")
      END IF
      success = .FALSE.
      IF (plan%regular) THEN
         start = 1
         IF (PRESENT(cell_first)) start = cell_first
         CPASSERT(start >= 1 .AND. start + SIZE(values_rs, 3) - 1 <= SIZE(plan%cell_factor))
         ngrid = PRODUCT(INT(plan%nfft, int_8))
         budget = fft_batch_bytes
         IF (PRESENT(max_workspace_bytes)) budget = MIN(budget, max_workspace_bytes)
         width = MIN(SIZE(values_rs, 1), INT(MAX(1_int_8, budget/(48_int_8*ngrid))))
         IF (width == 0) THEN
            success = .TRUE.
         ELSE
            needed = ngrid*width
            CPASSERT(needed <= HUGE(0))
            IF (ASSOCIATED(plan%fft_in)) THEN
               ! Reuse capacity while also allowing a smaller caller budget.
               IF (SIZE(plan%fft_in, KIND=int_8) < needed .OR. &
                   32_int_8*SIZE(plan%fft_in, KIND=int_8) + 16_int_8*needed > budget) THEN
                  CALL fft_dealloc(plan%fft_in)
                  CALL fft_dealloc(plan%fft_out)
               END IF
            END IF
            IF (.NOT. ASSOCIATED(plan%fft_in)) THEN
               CALL fft_alloc(plan%fft_in, [INT(needed)])
               CALL fft_alloc(plan%fft_out, [INT(needed)])
            END IF
            stat = 0
            DO j = 1, SIZE(values_rs, 2)
               DO first = 1, SIZE(values_rs, 1), width
                  nbatch = MIN(width, SIZE(values_rs, 1) - first + 1)
                  input(1:plan%nfft(1), 1:plan%nfft(2), 1:plan%nfft(3), 1:nbatch) => plan%fft_in(1:INT(ngrid)*nbatch)
                  input = z_zero
                  DO i = 1, nbatch
                     DO icell = 1, SIZE(values_rs, 3)
                        image = start + icell - 1
                        input(plan%cell_index(1, image), plan%cell_index(2, image), plan%cell_index(3, image), i) = &
                           input(plan%cell_index(1, image), plan%cell_index(2, image), plan%cell_index(3, image), i) + &
                           plan%cell_factor(image)*values_rs(first + i - 1, j, icell)
                     END DO
                  END DO
                  CALL lattice_fft_execute(BWFFT, plan%nfft, nbatch, plan%fft_in, plan%fft_out, stat)
                  IF (stat /= 0) EXIT
                  output(1:nbatch, 1:plan%nfft(1), 1:plan%nfft(2), 1:plan%nfft(3)) => plan%fft_in(1:INT(ngrid)*nbatch)
                  DO ik_out = 1, SIZE(values_k, 3)
                     ik = ik_out
                     IF (PRESENT(selected_kpoints)) ik = selected_kpoints(ik_out)
                     values_k(first:first + nbatch - 1, j, ik_out) = &
                        output(:, plan%grid_index(1, ik), plan%grid_index(2, ik), plan%grid_index(3, ik))
                  END DO
               END DO
               IF (stat /= 0) EXIT
            END DO
            success = stat == 0
         END IF
      END IF
      IF (.NOT. success) THEN
         IF (PRESENT(require_fft)) THEN
            IF (require_fft) CPABORT("Forced lattice FFT execution failed")
         END IF
         CALL direct_cell_to_k(values_rs, index_to_cell, xkp, values_k, deriv_direction, hmat, selected_kpoints)
      END IF
      IF (PRESENT(used_fft)) used_fft = success
      IF (.NOT. PRESENT(work)) CALL cell_to_k_grid_release(local_work)
      CALL timestop(handle)
   END SUBROUTINE cell_to_k_grid_fft

! **************************************************************************************************
!> \brief Prepare fixed forward lattice maps; the caller keeps geometry and k points unchanged.
!> \param work owned image/grid maps and reusable FFT buffers
!> \param index_to_cell image translations, including repeated cells from different atom pairs
!> \param xkp full reciprocal grid
!> \param nkp_grid reciprocal-grid dimensions
!> \param deriv_direction optional Cartesian derivative direction
!> \param hmat direct-lattice cell matrix for a derivative
! **************************************************************************************************
   SUBROUTINE cell_to_k_grid_prepare(work, index_to_cell, xkp, nkp_grid, deriv_direction, hmat)
      TYPE(cell_to_k_grid_work_type), INTENT(INOUT)      :: work
      INTEGER, INTENT(IN)                                :: index_to_cell(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:, :)
      INTEGER, INTENT(IN)                                :: nkp_grid(3)
      INTEGER, INTENT(IN), OPTIONAL                      :: deriv_direction
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: hmat(3, 3)

      INTEGER                                            :: d, icell
      REAL(KIND=dp)                                      :: arg, cell_vector(3), k_offset(3)

      CALL cell_to_k_grid_release(work)
      ALLOCATE (work%grid_index(3, SIZE(xkp, 2)))
      work%regular = regular_kpoint_grid(xkp, nkp_grid, work%grid_index, k_offset, &
                                         allow_incomplete=.TRUE.)
      IF (work%regular) work%regular = lattice_fft_shape(nkp_grid, work%nfft, allow_arbitrary=.TRUE.)
      IF (.NOT. work%regular) RETURN
      IF (PRESENT(deriv_direction)) THEN
         CPASSERT(PRESENT(hmat))
         CPASSERT(deriv_direction >= 1 .AND. deriv_direction <= 3)
      END IF
      ALLOCATE (work%cell_index(3, SIZE(index_to_cell, 2)), work%cell_factor(SIZE(index_to_cell, 2)))
      DO icell = 1, SIZE(index_to_cell, 2)
         DO d = 1, 3
            work%cell_index(d, icell) = MODULO(index_to_cell(d, icell), nkp_grid(d))* &
                                        (work%nfft(d)/nkp_grid(d)) + 1
         END DO
         arg = SUM(k_offset*REAL(index_to_cell(1:3, icell), dp))
         work%cell_factor(icell) = EXP(gaussi*twopi*arg)
         IF (PRESENT(deriv_direction)) THEN
            cell_vector = MATMUL(hmat, REAL(index_to_cell(1:3, icell), dp))
            work%cell_factor(icell) = work%cell_factor(icell)*gaussi*cell_vector(deriv_direction)
         END IF
      END DO
   END SUBROUTINE cell_to_k_grid_prepare

! **************************************************************************************************
!> \brief Release forward lattice maps and numerical buffers.
!> \param work caller-owned prepared transform
! **************************************************************************************************
   SUBROUTINE cell_to_k_grid_release(work)
      TYPE(cell_to_k_grid_work_type), INTENT(INOUT)      :: work

      IF (ASSOCIATED(work%fft_in)) CALL fft_dealloc(work%fft_in)
      IF (ASSOCIATED(work%fft_out)) CALL fft_dealloc(work%fft_out)
      IF (ALLOCATED(work%grid_index)) DEALLOCATE (work%grid_index)
      IF (ALLOCATED(work%cell_index)) DEALLOCATE (work%cell_index, work%cell_factor)
      work%regular = .FALSE.
   END SUBROUTINE cell_to_k_grid_release

! **************************************************************************************************
!> \brief Transform reciprocal-grid data back to image-cell values.
!> \param values_k complex matrices on the supplied reciprocal points
!> \param xkp reciprocal coordinates, in units of reciprocal lattice vectors
!> \param nkp_grid dimensions of the reciprocal grid
!> \param index_to_cell integer lattice vector for every requested image cell
!> \param values_rs reconstructed complex matrices in the requested image cells
!> \param used_fft reports whether the regular-grid FFT path was used
!> \param weights optional k weights; otherwise use the existing 1/number-of-points normalization
!> \param allow_incomplete embed a unique regular-grid subset, zero-filling missing points
! **************************************************************************************************
   SUBROUTINE k_grid_to_cell_fft(values_k, xkp, nkp_grid, index_to_cell, values_rs, used_fft, &
                                 weights, allow_incomplete)

      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: values_k
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      INTEGER, DIMENSION(3), INTENT(IN)                  :: nkp_grid
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: values_rs
      LOGICAL, INTENT(OUT), OPTIONAL                     :: used_fft
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: weights
      LOGICAL, INTENT(IN), OPTIONAL                      :: allow_incomplete

      TYPE(k_grid_to_cell_work_type)                     :: work

      CALL k_grid_to_cell_prepare(work, xkp, nkp_grid, index_to_cell, weights, allow_incomplete)
      CALL k_grid_to_cell_execute(work, values_k, values_rs, used_fft)
      CALL k_grid_to_cell_release(work)

   END SUBROUTINE k_grid_to_cell_fft

! **************************************************************************************************
!> \brief Prepare one inverse lattice transform for repeated matrix batches.
!> \param work call-scoped maps, weights and buffers; no backend plan ownership
!> \param xkp reciprocal coordinates
!> \param nkp_grid candidate reciprocal-grid dimensions
!> \param index_to_cell requested lattice translations
!> \param weights final k weights, otherwise 1/number-of-points
!> \param allow_incomplete allow missing grid points, with zero values
! **************************************************************************************************
   SUBROUTINE k_grid_to_cell_prepare(work, xkp, nkp_grid, index_to_cell, weights, allow_incomplete)

      TYPE(k_grid_to_cell_work_type), INTENT(INOUT)      :: work
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      INTEGER, DIMENSION(3), INTENT(IN)                  :: nkp_grid
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: weights
      LOGICAL, INTENT(IN), OPTIONAL                      :: allow_incomplete

      INTEGER                                            :: d, icell, ncell, nkp
      REAL(KIND=dp), DIMENSION(3)                        :: k_offset

      CALL k_grid_to_cell_release(work)
      nkp = SIZE(xkp, 2)
      ncell = SIZE(index_to_cell, 2)
      IF (SIZE(xkp, 1) < 3 .OR. SIZE(index_to_cell, 1) < 3 .OR. nkp == 0) THEN
         CALL cp_abort(__LOCATION__, "Inconsistent reciprocal-grid mapping in inverse lattice FFT")
      END IF
      ALLOCATE (work%xkp(3, nkp), work%cells(3, ncell), work%weight(nkp), work%grid_index(3, nkp))
      work%xkp(:, :) = xkp(1:3, :)
      work%cells(:, :) = index_to_cell(1:3, :)
      work%weight(:) = 1.0_dp/REAL(nkp, dp)
      IF (PRESENT(weights)) THEN
         CPASSERT(SIZE(weights) == nkp)
         work%weight(:) = weights
      END IF
      work%regular = regular_kpoint_grid(xkp, nkp_grid, work%grid_index, k_offset, allow_incomplete)
      IF (work%regular) work%regular = lattice_fft_shape(nkp_grid, work%nfft, allow_arbitrary=.TRUE.)
      IF (.NOT. work%regular) RETURN

      ALLOCATE (work%cell_index(3, ncell), work%cell_phase(ncell))
      DO icell = 1, ncell
         DO d = 1, 3
            work%cell_index(d, icell) = MODULO(index_to_cell(d, icell), nkp_grid(d))* &
                                        (work%nfft(d)/nkp_grid(d)) + 1
         END DO
         work%cell_phase(icell) = EXP(-gaussi*twopi*SUM(k_offset*REAL(index_to_cell(1:3, icell), dp)))
      END DO

   END SUBROUTINE k_grid_to_cell_prepare

! **************************************************************************************************
!> \brief Execute a matrix batch using the prepared inverse lattice transform.
!>        Call outside OpenMP worker regions: fft_tools owns its shared pool and planner.
!> \param work prepared transform; scratch is reused, never shared between concurrent callers
!> \param values_k matrices at the prepared reciprocal points
!> \param values_rs matrices at the prepared lattice translations
!> \param used_fft whether this batch used FFT; false denotes a direct fallback or no result
!> \param require_fft abort on an unsupported mapping or FFT failure, without a direct fallback
!> \param allow_direct_fallback execute the generic direct sum when FFT is unavailable; defaults to true
! **************************************************************************************************
   SUBROUTINE k_grid_to_cell_execute(work, values_k, values_rs, used_fft, require_fft, allow_direct_fallback)

      TYPE(k_grid_to_cell_work_type), INTENT(INOUT)      :: work
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: values_k
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: values_rs
      LOGICAL, INTENT(OUT), OPTIONAL                     :: used_fft
      LOGICAL, INTENT(IN), OPTIONAL                      :: require_fft, allow_direct_fallback

      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: input(:, :, :, :), output(:, :, :, :)
      INTEGER                                            :: first, handle, i, icell, ik, j, nbatch, &
                                                            stat, width
      INTEGER(KIND=int_8)                                :: ngrid

      CALL timeset("k_grid_to_cell_fft", handle)
      CPASSERT(ALLOCATED(work%xkp))
      IF (SIZE(values_k, 3) /= SIZE(work%xkp, 2) .OR. &
          SIZE(values_rs, 1) /= SIZE(values_k, 1) .OR. &
          SIZE(values_rs, 2) /= SIZE(values_k, 2) .OR. &
          SIZE(values_rs, 3) /= SIZE(work%cells, 2)) THEN
         CALL cp_abort(__LOCATION__, "Inconsistent matrix batches in inverse lattice FFT")
      END IF

      stat = 0
      IF (PRESENT(require_fft)) THEN
         IF (require_fft .AND. .NOT. work%regular) THEN
            CPABORT("Forced inverse lattice FFT needs a supported regular grid")
         END IF
      END IF
      IF (work%regular) THEN
         ngrid = PRODUCT(INT(work%nfft, int_8))
         width = MIN(SIZE(values_k, 1), INT(MAX(1_int_8, fft_batch_bytes/(48_int_8*ngrid))))
         IF (width == 0) THEN
            IF (PRESENT(used_fft)) used_fft = .TRUE.
            CALL timestop(handle)
            RETURN
         END IF
         CPASSERT(ngrid*width <= HUGE(0))
         IF (ASSOCIATED(work%fft_in)) THEN
            IF (SIZE(work%fft_in, KIND=int_8) < ngrid*width) THEN
               CALL fft_dealloc(work%fft_in)
               CALL fft_dealloc(work%fft_out)
            END IF
         END IF
         IF (.NOT. ASSOCIATED(work%fft_in)) THEN
            CALL fft_alloc(work%fft_in, [INT(ngrid)*width])
            CALL fft_alloc(work%fft_out, [INT(ngrid)*width])
         END IF
         DO j = 1, SIZE(values_k, 2)
            DO first = 1, SIZE(values_k, 1), width
               nbatch = MIN(width, SIZE(values_k, 1) - first + 1)
               input(1:work%nfft(1), 1:work%nfft(2), 1:work%nfft(3), 1:nbatch) => &
                  work%fft_in(1:INT(ngrid)*nbatch)
               input = z_zero
               DO i = 1, nbatch
                  DO ik = 1, SIZE(values_k, 3)
                     input(work%grid_index(1, ik), work%grid_index(2, ik), work%grid_index(3, ik), i) = &
                        work%weight(ik)*values_k(first + i - 1, j, ik)
                  END DO
               END DO
               CALL lattice_fft_execute(FWFFT, work%nfft, nbatch, work%fft_in, work%fft_out, stat)
               IF (stat /= 0) EXIT
               output(1:nbatch, 1:work%nfft(1), 1:work%nfft(2), 1:work%nfft(3)) => &
                  work%fft_in(1:INT(ngrid)*nbatch)
               DO icell = 1, SIZE(values_rs, 3)
                  ! The 1D passes are unnormalized; supplied k weights are already final.
                  values_rs(first:first + nbatch - 1, j, icell) = work%cell_phase(icell)* &
                                          output(:, work%cell_index(1, icell), work%cell_index(2, icell), work%cell_index(3, icell))
               END DO
            END DO
            IF (stat /= 0) EXIT
         END DO
      END IF
      IF (stat /= 0) work%regular = .FALSE.
      IF (.NOT. work%regular) THEN
         IF (PRESENT(require_fft)) THEN
            IF (require_fft) CPABORT("Forced inverse lattice FFT execution failed")
         END IF
         IF (PRESENT(allow_direct_fallback)) THEN
            IF (.NOT. allow_direct_fallback) THEN
               values_rs = z_zero
               IF (PRESENT(used_fft)) used_fft = .FALSE.
               CALL timestop(handle)
               RETURN
            END IF
         END IF
         CALL direct_k_to_cell(values_k, work%xkp, work%cells, values_rs, work%weight)
      END IF
      IF (PRESENT(used_fft)) used_fft = work%regular
      CALL timestop(handle)

   END SUBROUTINE k_grid_to_cell_execute

! **************************************************************************************************
!> \brief Release call-scoped inverse lattice maps and FFT buffers.
!> \param work ...
! **************************************************************************************************
   SUBROUTINE k_grid_to_cell_release(work)

      TYPE(k_grid_to_cell_work_type), INTENT(INOUT)      :: work

      IF (ASSOCIATED(work%fft_in)) CALL fft_dealloc(work%fft_in)
      IF (ASSOCIATED(work%fft_out)) CALL fft_dealloc(work%fft_out)
      IF (ALLOCATED(work%xkp)) DEALLOCATE (work%xkp, work%cells, work%weight, work%grid_index)
      IF (ALLOCATED(work%cell_index)) DEALLOCATE (work%cell_index, work%cell_phase)
      work%regular = .FALSE.

   END SUBROUTINE k_grid_to_cell_release

! **************************************************************************************************
!> \brief Direct inverse phase sum for arbitrary reciprocal-point lists.
!> \param values_k ...
!> \param xkp ...
!> \param index_to_cell ...
!> \param values_rs ...
!> \param weights optional k weights, otherwise 1/number-of-points
! **************************************************************************************************
   SUBROUTINE direct_k_to_cell(values_k, xkp, index_to_cell, values_rs, weights)

      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN)   :: values_k
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: values_rs
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: weights

      COMPLEX(KIND=dp)                                   :: factor
      INTEGER                                            :: icell, ik
      REAL(KIND=dp)                                      :: arg, normalization

      values_rs = z_zero
      normalization = 1.0_dp/REAL(SIZE(xkp, 2), KIND=dp)
      DO icell = 1, SIZE(index_to_cell, 2)
         DO ik = 1, SIZE(xkp, 2)
            IF (PRESENT(weights)) normalization = weights(ik)
            arg = SUM(xkp(1:3, ik)*REAL(index_to_cell(1:3, icell), KIND=dp))
            factor = normalization*EXP(-gaussi*twopi*arg)
            values_rs(:, :, icell) = values_rs(:, :, icell) + factor*values_k(:, :, ik)
         END DO
      END DO

   END SUBROUTINE direct_k_to_cell

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
